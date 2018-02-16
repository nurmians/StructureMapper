# This file is part of the StructureMapper algorithm.
# Please cite the authors if you find this software useful
#
# https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361
 
# MIT License
#
# Copyright 2018 Anssi Nurminen and Vesa P. Hytönen
# Faculty of Medicine and Life Sciences and BioMediTech, University of Tampere, Arvo Ylpön katu 34, 33520 Tampere, Finland 
# Fimlab Laboratories, Biokatu 4, 33520 Tampere, Finland
# Tampere, Finland
# 16-02-2018

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys
import os
import ConfigParser
from subprocess import Popen, PIPE
import time
import re

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append( "%s/../utils" % curdir)
from fileops import InsertFolderSeparator, GetFolder, DriveLetter, SystemSpecificFolderSeparators

class Timeout( Exception):
    pass

def RunWithIt( aCmd, aTimeout=15):
    proc = Popen( aCmd, bufsize=0, stdout=PIPE, stderr=PIPE, shell=True)
    poll_seconds = .250
    deadline = time.time()+aTimeout
    while time.time() < deadline and proc.poll() == None:
        time.sleep( poll_seconds)

    if proc.poll() == None:
        #Process did not finish in time
        if float(sys.version[:3]) >= 2.6:
            proc.terminate()
            proc.wait()
        #raise Timeout()
    else:
        proc.wait()

    stdout, stderr = proc.communicate()
    #print "STDOUT: %s" % stdout
    #print "STDERR:", stderr
    return stdout, stderr, proc.returncode


def FindFolderThatContainsFile( aWildCardFile, aFilePattern, aStartPath=None, aLookFromDriveC=False):
    return GetFolder( FindFile( aWildCardFile, aFilePattern, aStartPath, aLookFromDriveC))



#aFile need to be in the form of a regexp pattern
def FindFile( aWildCardFile, aFilePattern, aStartPath=None, aLookFromDriveC=False):


    if os.name == 'nt':
        if not aStartPath: aStartPath = "C:"

        print "Looking for file '%s' in '%s'." % (aWildCardFile, aStartPath)

        file_pattern = re.compile( " (%s)" % aFilePattern)
        dir_pattern = re.compile("\s([A-Z]:\\\\.*)$") #Look for drive letters
        try:
            #find_cmd = "dir /S %s%s > I:\\test.txt" % ( InsertFolderSeparator( aStartPath), aWildCardFile)
            find_cmd = "dir /S %s%s" % ( InsertFolderSeparator( aStartPath), aWildCardFile)
            #print "Looking with cmd: '%s'" % find_cmd
            stdout, stderr, returncode = RunWithIt( find_cmd)
        except Timeout:
            pass

        directory = ""
        file = ""
        for line in stdout.split("\n"):

            line = line.strip()
            #print "line: ", line

            m = dir_pattern.search( line)
            if m:
                directory = m.group( 0).strip() + "\\"
                continue

            if line.find("<DIR>") >= 0: continue

            m = file_pattern.search( line)
            if m and len( directory):
                file = m.group( 0).strip()
                filepath = directory + file
                if os.path.isfile( filepath):
                    print "File found: '%s'" % filepath
                    return filepath

        #Not found, look from drive C as well
        if aStartPath.find("C:") < 0 and aLookFromDriveC:
            return FindFile( aWildCardFile, aFilePattern, aStartPath="C:", aLookFromDriveC=False)

    else:

        aWildCardFile = aWildCardFile.replace( ".exe", "")
        aFilePattern = aFilePattern.replace( "\\.exe", "")

        if not aStartPath: aStartPath = "/"

        print "Looking for file '%s' in '%s'." % (aWildCardFile, aStartPath)

        stdout, stderr, returncode = RunWithIt("find / -name '%s' 2> /dev/null" % aWildCardFile)
        for line in stdout.split("\n"):

            filepath = line.strip()
            if os.path.isfile( filepath):
                    return filepath

    print "WARNING: Could not find file '%s' in '%s'" % (aWildCardFile, aStartPath)
    return ""



def SetConfig( aTuple, aCustomDict, aConfig):
    section, key, val = aTuple
    if key in aCustomDict: val = str( aCustomDict[ key])
    aConfig.set( section, key, val)

def WriteConfig( aConfigFile, aCustomSettings={}, aAutoConfig=False):

    config = ConfigParser.RawConfigParser( allow_no_value = True)
    config.optionxform = str

    try:

        config.add_section('Directories')

        config.set('Directories', '# Paths can be relative ("../") or absolute ("C:\\")')
        config.set('Directories', '# Use empty values to disable features')

        default_pdb = './pdb' if not aAutoConfig else FindFolderThatContainsFile( "pdb????.ent", "pdb....\\.ent", DriveLetter(), aLookFromDriveC=True)
        if aAutoConfig: print "PDB files folder set to '%s'." % default_pdb
        if aAutoConfig and not len( default_pdb): default_pdb = './pdb'
        config.set('Directories', '# File storage directory for downloading PDB database structure files')
        SetConfig(('Directories', 'PDB_DIR', default_pdb), aCustomSettings, config)

        default_val = InsertFolderSeparator( default_pdb) + "ASA"
        if aAutoConfig and len( default_val): print "PDB ASA files folder set to '%s'." % default_val
        config.set('Directories', '# File storage directory for precalculated surface area calculations (ASA per residue)')
        SetConfig(('Directories', 'PRECALC_ASA_DIR', default_val), aCustomSettings, config)

        default_val = InsertFolderSeparator( default_pdb) + "ASA_SC"
        if aAutoConfig: print "PDB ASA_SC files folder set to '%s'." % default_val
        config.set('Directories', '# File storage directory for precalculated surface area calculations (ASA per residue sidechain)')
        SetConfig(('Directories', 'PRECALC_ASA_SC_DIR', default_val), aCustomSettings, config)

        default_val = './blast' if not aAutoConfig else FindFolderThatContainsFile( "pdbaa.00.phr", "pdbaa\\.00\\.phr", DriveLetter(), aLookFromDriveC=True)
        config.set('Directories', '# Local Blast database directory, needs to contain database "pdbaa".')
        SetConfig(('Directories', 'BLAST_DB', default_val), aCustomSettings, config)

        default_val = './dssp' if not aAutoConfig else FindFolderThatContainsFile( "????.dssp", "....\\.dssp", DriveLetter(), aLookFromDriveC=True)
        if aAutoConfig: print "DSSP folder for precalculated dssp files set to '%s'." % default_val
        config.set('Directories', '# Directory where to store DSSP secondary structure predictions')
        SetConfig(('Directories', 'DEFAULT_DSSP_DIR', default_val), aCustomSettings, config)

        default_val = SystemSpecificFolderSeparators('../results')
        config.set('Directories', '# Directory for the results')
        if aAutoConfig: print "Results will be stored in '%s'." % default_val
        SetConfig(('Directories', 'DEFAULT_RESULT_DIR', default_val), aCustomSettings, config)

        default_val = '' if not aAutoConfig else FindFolderThatContainsFile( "blastp.exe", "blast\\.exe", DriveLetter(), aLookFromDriveC=True)
        if aAutoConfig and len( default_val): print "Blastp executable found in '%s'.  [OK]" % default_val
        elif aAutoConfig: print "WARNING: Could not locate Blastp executable. This setting need to be set before using the algorithm (if the executable is not in system path)."
        config.set('Directories', '# BLASTP (protein Blast) executable for finding holomolgous structures.')
        SetConfig(('Directories', 'DEFAULT_BLASTP', default_val), aCustomSettings, config)

        default_val = '' if not aAutoConfig else FindFile( "iupred.exe", "iupred\\.exe", DriveLetter(), aLookFromDriveC=True)
        if aAutoConfig and len( default_val): print "IUPRED executable found in '%s'.  [OK]" % default_val
        config.set('Directories', '# IUPRED executable for predicting disorder.')
        SetConfig(('Directories', 'DEFAULT_IUPRED', default_val), aCustomSettings, config)

        default_val = '' if not aAutoConfig else FindFile( "dssp*.exe", "dssp.*\\.exe", DriveLetter(), aLookFromDriveC=True)
        if aAutoConfig and len( default_val): print "DSSP executable found in '%s'.  [OK]" % default_val
        config.set('Directories', '# DSSP executable for calculating secondary structure.')
        SetConfig(('Directories', 'DEFAULT_DSSP', default_val), aCustomSettings, config)

        default_val = '' if (not aAutoConfig or os.name != 'nt') else FindFile( "PymolWin.exe", "PymolWin\\.exe", DriveLetter(), aLookFromDriveC=True)
        config.set('Directories', '# PyMol executable for viewing results. Separate multiple alternatives with ",".')
        SetConfig(('Directories', 'PYMOL_PATH', default_val), aCustomSettings, config)

        config.set('Directories', '# Path for storing temporary files')
        SetConfig(('Directories', 'WORK_DIR', SystemSpecificFolderSeparators('./')), aCustomSettings, config)

        #ONLINE
        config.add_section('Online')
        config.set('Online', '# FTP URL for DSSP predictions')
        SetConfig(('Online', 'DEFAULT_DSSP_FTP', "ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/"), aCustomSettings, config)

        #FILTERING
        config.add_section('Filtering')

        config.set('Filtering', '# Minimum bitscore for usig blast homology structure (0 for no filter)')
        SetConfig(('Filtering', 'BITSCORE_THRESHOLD', "50.0"), aCustomSettings, config)

        config.set('Filtering', '# Max number of BLAST homologs to consider in evalutation (0 for all available)')
        SetConfig(('Filtering', 'MAX_BLAST_HOMOLOGS', "0"), aCustomSettings, config)

        config.set('Filtering', '# Minimum alignment score between sequence and structure (0-100) for using structure')
        SetConfig(('Filtering', 'ALIGNMENT_THRESHOLD', "70"), aCustomSettings, config)

        #config.set('Filtering', '# How many amino acids at the POI are used for finding the POI in the structure (-1 for full seq alignment)')
        #config.set('Filtering', 'ALIGNMENT_WINDOW', "10")

        config.set('Filtering', '# IUPRED sequence window size')
        SetConfig(('Filtering', 'IUPRED_SEQWINDOW', '8'), aCustomSettings, config)

        config.set('Filtering', '# Default number of amino acids to consider at the N-terminal side of the POI (0 for both neighboring aas)')
        SetConfig(('Filtering', 'DEFAULT_POI_POSITIONS', '0'), aCustomSettings, config)

        config.set('Filtering', '# Set to 1 to consider only sidechains in calculations as the default mode')
        SetConfig(('Filtering', 'DEFAULT_ONLY_SIDECHAINS', '0'), aCustomSettings, config)

        print "Writing config file to '%s'..." % aConfigFile,
        with open( aConfigFile, 'wb') as configfile:
            config.write( configfile)
        print "Done."

    except Exception as ex:
       sys.stderr.write("ERROR: Could not write configuration file '%s'.\n" % aConfigFile)
       message = "An exception of type {0} occured. Arguments: {1!r}\n".format(type(ex).__name__, ex.args)
       sys.stderr.write( message)
       raise
       return -4

    return 0


def ReadConfig( aConfigFile):

    import ConfigParser

    config = ConfigParser.RawConfigParser( allow_no_value = True)
    config.optionxform = str

    configurations = {}

    #CREATE default config file if it does not exist
    if not os.path.isfile( aConfigFile) or int( os.stat( aConfigFile).st_size) <= 0:

        sys.stdout.write("INFO: Creating config file with default values...\n")
        sys.stdout.write("INFO: Please specify default directories in file '%s'.\n" % aConfigFile)

        sys.exit( WriteConfig( aConfigFile))
        #Exit( WriteConfig( aConfigFile)) #Write default config file

    else:

        all_settings = []
        #                   var_type,mandatory,  section,      option,             default
        all_settings.append( ('str',   True,  'Directories', 'PRECALC_ASA_DIR',        ""))
        all_settings.append( ('str',   True,  'Directories', 'PRECALC_ASA_SC_DIR',     ""))
        all_settings.append( ('str',   True,  'Directories', 'PDB_DIR',                ""))
        all_settings.append( ('str',   True,  'Directories', 'BLAST_DB',               ""))
        all_settings.append( ('str',   False, 'Directories', 'DEFAULT_BLASTP',         ""))
        all_settings.append( ('str',   True,  'Directories', 'DEFAULT_DSSP_DIR',       ""))
        all_settings.append( ('str',   True,  'Directories', 'DEFAULT_IUPRED',         ""))
        all_settings.append( ('str',   True,  'Directories', 'DEFAULT_DSSP',           ""))
        all_settings.append( ('str',   False, 'Directories', 'WORK_DIR',               "./"))
        #all_settings.append( ('str',   True,  'Directories', 'DEFAULT_HUM_PROTEOME',   ""))
        all_settings.append( ('str',   True,  'Directories', 'DEFAULT_RESULT_DIR',     ""))
        all_settings.append( ('str',   False, 'Directories', 'PYMOL_PATH',             ""))
        all_settings.append( ('str',   True,  'Online',      'DEFAULT_DSSP_FTP',       ""))
        all_settings.append( ('float', True,  'Filtering',   'BITSCORE_THRESHOLD',     60.0))
        all_settings.append( ('int', False,   'Filtering',   'MAX_BLAST_HOMOLOGS',     30)) #default 5
        all_settings.append( ('int',   True,  'Filtering',   'ALIGNMENT_THRESHOLD',    70))
        #all_settings.append( ('int',   False, 'Filtering',   'ALIGNMENT_WINDOW',       10))
        all_settings.append( ('int',   False, 'Filtering',   'IUPRED_SEQWINDOW',       8))
        all_settings.append( ('int',   False, 'Filtering',   'DEFAULT_POI_POSITIONS',  0)) #0 == between
        all_settings.append( ('int',   False, 'Filtering',   'DEFAULT_ONLY_SIDECHAINS',0))

        missing = False

        config.read( aConfigFile)

        for var_type, mandatory, section, option, default in all_settings:

            try:

                value = config.get( section, option).replace("\"", "")#Remove quotes, if used
                #Typecast
                if var_type == "int": value = int( value)
                elif var_type == "float": value = float( value)
                elif var_type == "list": value = value.split(",")

                configurations[ option] = value

            except ValueError:
                sys.stderr.write( "ERROR: Configuration file: %s.\n" % str( noe))
                sys.stderr.write( "Argument '%s' should be of type %s.\n" % (noe.args[ 0], var_type))
                missing = True
            except ConfigParser.NoOptionError as noe:

                if mandatory:
                    sys.stderr.write( "ERROR: Configuration file: %s.\n" % str( noe))
                    sys.stderr.write( "Add missing option '%s' or remove file '%s' to reset to default values.\n" % (noe.args[ 0], aConfigFile))
                    missing = True
                else:
                    configurations[ option] = default

            except ConfigParser.NoSectionError as nse:
                sys.stderr.write( "ERROR: Configuration file: %s.\n" % str( nse))
                sys.stderr.write( "Add missing Section '%s' or remove file '%s' to reset to default values.\n" % (nse.args[ 0], aConfigFile))
                missing = True
            except Exception as ex:
                sys.stderr.write( "ERROR: Could not read configuration file '%s'.\n" % aConfigFile)
                message = "An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args)
                sys.stderr.write( message)
                sys.exit( -5)


        if missing: sys.exit( -1)

    return configurations


if __name__ == "__main__":

    #Testing
    #FindFile( "I:\\", "iupred.exe", "iupred\.exe")
    #FindFile( "I:\\", "pdb*.ent", "pdb.*\.ent")
    WriteConfig( "test.ini", aAutoConfig=True)

