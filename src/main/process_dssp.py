#!/usr/bin/python

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
import urllib #dssp fetching
from subprocess import call

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append( "%s/../utils" % curdir)

from struct_methods import *
from fileops import *
import online_pdb

DSSP_EXECUTABLE = ""

def SetDSSP_Executable( aPathToExe):
    global DSSP_EXECUTABLE
    DSSP_EXECUTABLE = aPathToExe

def GetDSSP_Executable():
    global DSSP_EXECUTABLE
    return DSSP_EXECUTABLE

#Return 0 on success
#Run DSSP algorithm (executable) on local machine
def RunDSSP( aPdbFile, aOutputfile, aOverWriteExisting=True, aDSSP_Executable="", aLogFile=None):

    #global DSSP_EXECUTABLE
    #global CONFIG_VALUES

    #CHECK EXE
    if len( aDSSP_Executable) == 0:
        aDSSP_Executable = GetDSSP_Executable()

    if len( aDSSP_Executable) == 0:
        if not hasattr( RunDSSP, "exe_not_set_warned"):
            LogIt( "DSSP WARNING: DSSP binary executable not set in config.ini.", aLogFile, -2) #Show once
            RunDSSP.exe_not_set_warned = 1
        return -1 #Dssp exe not set

    dssp_exe    = SystemSpecificFolderSeparators( aDSSP_Executable) #CONFIG_VALUES['DEFAULT_DSSP'])
    aPdbFile    = SystemSpecificFolderSeparators( aPdbFile)
    aOutputfile = SystemSpecificFolderSeparators( aOutputfile)


    if len( aDSSP_Executable) != 0 and not os.path.isfile( dssp_exe):
        if not hasattr( RunDSSP, "exe_not_found_warned"):
            LogIt( "DSSP ERROR: DSSP binary executable set but not found at '%s'." % dssp_exe, aLogFile, -2) #Show once
            RunDSSP.exe_not_found_warned = 1
        return -50 #Dssp exe not set

    if not aOverWriteExisting and os.path.isfile( aOutputfile) and int( os.stat( aOutputfile).st_size) > 0:
        return 0 #File exists

    if not os.path.isfile( aPdbFile):
        LogIt( "DSSP ERROR: File '%s' does not exists.\n" % aPdbFile, aLogFile, -2)
        return -1

    cmd = dssp_exe + " -i %s > %s" % (aPdbFile, aOutputfile)
    #print "CMD:", cmd
    retval = call( cmd, shell=True)
    #print "DSSP_RETVAL:", retval
    return retval

def FetchDSSP( aPDB_ID, aDSSP_folder, aDSSP_filename, aPdbDir=None, aLogFile=None):

    if not hasattr( FetchDSSP, "fetched"): FetchDSSP.fetched = 0
    if not hasattr( FetchDSSP, "errors"): FetchDSSP.errors = 0
    if not hasattr( FetchDSSP, "consecutive_errors"): FetchDSSP.consecutive_errors = 0
    if not hasattr( FetchDSSP, "unfetchables"): FetchDSSP.unfetchables = {}
    if not hasattr( FetchDSSP, "timeouts"): FetchDSSP.timeouts = 0

    DSSPbaseUrl = FetchDSSP.DEFAULT_DSSP_FTP

    dssp_filename = aDSSP_filename
    dssp_path = InsertFolderSeparator( aDSSP_folder) + dssp_filename

    runOffline = aPdbDir and len( aPdbDir)

    if runOffline and not os.path.exists( aPdbDir):
        raise RuntimeError( "DSSP ERROR: Not a valid PDB folder path: '%s'." % aPdbDir)
        #sys.stderr.write( )
        #Exit( -5)

    #If file does not exist, create (run DSSP locally) or fetch (download) it
    if not os.path.isfile( dssp_path) or int( os.stat( dssp_path).st_size) <= 0:


        if not runOffline and (not DSSPbaseUrl or not len( DSSPbaseUrl)):
            #sys.stderr.write( "DSSP ERROR: DSSP URL not specified. Please specify download URL in config.ini with key 'DEFAULT_DSSP_FTP'.\n")
            raise RuntimeError( "DSSP ERROR: DSSP URL not specified. Please specify download URL in config.ini with key 'DEFAULT_DSSP_FTP'.")

        fetch_retries = 0 #0 == No retries
        fetch_success = False

        while not fetch_success and fetch_retries >= 0 and aPDB_ID.upper() not in FetchDSSP.unfetchables:

            #FETCH ONLINE
            try:
                fetch_retval = -1
                #If DSSP Executable set, run DSSP
                if runOffline:

                    #USE OFFLINE IF POSSIBLE
                    url = InsertFolderSeparator( aPdbDir)+PDBFileFromAccCode( aPDB_ID) #offline url (path)
                    fetch_retval = RunDSSP( url, dssp_path, aOverWriteExisting=False, aLogFile=aLogFile )

                    if fetch_retval == 0:
                        fetch_success = True
                        LogIt( "DSSP INFO: Calculated DSSP for structure %s to file '%s'" % (aPDB_ID.upper(), dssp_filename), aLogFile, 1)
                    elif fetch_retval == -50:
                        #DSSP executable set but not found
                        raise RuntimeError( "DSSP executable set at '%s' but not found." % GetDSSP_Executable())
                        #Exit( -50)
                    else:
                        LogIt( "DSSP INFO: Failed to calculate DSSP for file '%s' (Only backbone available?).\n" % dssp_filename, aLogFile, 2)
                        FetchDSSP.unfetchables[ aPDB_ID.upper()] = True
                else:
                    #ONLINE
                    url = InsertFolderSeparator( DSSPbaseUrl) + dssp_filename

                    import socket
                    socket.setdefaulttimeout( 7)

                    urllib.urlretrieve( url, dssp_path)
                    LogIt( "DSSP INFO: Fetched: %s" % dssp_filename, aLogFile, 1)
                    fetch_success = True
                    FetchDSSP.fetched += 1
                    FetchDSSP.consecutive_errors = 0

                    if FetchDSSP.fetched == 20:
                        print "DSSP INFO: You can download the full DSSP database from '%s' using FTP to speed up this step of the processing." % DSSPbaseUrl
                    break

            except IOError as ex:
                sys.stderr.write( "DSSP ERROR: %s: Error fetching DSSP file: '%s' from '%s'.\n" % (aPDB_ID.upper(), dssp_filename, url))
                sys.stderr.write( "            An exception of type {0} occured. Arguments: {1!r}\n".format(type(ex).__name__, ex.args))
                FetchDSSP.unfetchables[ aPDB_ID.upper()] = True
                FetchDSSP.timeouts += 1
            except Exception as ex:
                sys.stderr.write( "DSSP ERROR: %s: Error fetching DSSP file: '%s' from '%s'.\n" % (aPDB_ID.upper(), dssp_filename, url))
                sys.stderr.write( "            An exception of type {0} occured. Arguments: {1!r}\n".format(type(ex).__name__, ex.args))

                if fetch_retries > 0:

                    sys.stderr.write( "DSSP: Retrying DSSP Fetch of '%s'...\n" % aPDB_ID.upper())
                    fetch_retries -= 1
                    time.sleep( 1)
                else:

                    #No more retries
                    FetchDSSP.errors += 1
                    FetchDSSP.consecutive_errors += 1
                    FetchDSSP.unfetchables[ aPDB_ID.upper()] = True
                    LogIt( "DSSP WARNING: DSSP for %s could not be calculated or retrieved from '%s'.\n" % (aPDB_ID.upper(), url), aLogFile, -2)
                    LogIt( "              This could be due to incomplete atom information (Only CAs?).", aLogFile, -2)
                    if not runOffline and FetchDSSP.consecutive_errors >= 20: #Only online in use
                        raise RuntimeError( "DSSP ERROR: Too many consecutive errors fetching DSSP files. Quitting...")
                        #Exit( -7)
                    break

        if not fetch_success:
            return (False, "Unable to fetch DSSP file") #Return

        return (True, "OK")


def GetDSSPpath( pdb_id, DSSP_folder):

    dssp_filename = pdb_id.lower() + ".dssp"
    dssp_path = InsertFolderSeparator( DSSP_folder) + dssp_filename
    return dssp_path

#Returns file contents as list of lines
def ReadDSSP( pdb_id, DSSP_folder, aPdbFolder, fetch_if_missing=True, log_file=None):

    dssp_filename = pdb_id.lower() + ".dssp"
    dssp_path = InsertFolderSeparator( DSSP_folder) + dssp_filename

    if not os.path.isfile( dssp_path):

        if fetch_if_missing:
            retval, msg = FetchDSSP( pdb_id, DSSP_folder, dssp_filename, aPdbDir=aPdbFolder, aLogFile=log_file )

            #Fetch failed?
            if retval == False:
                LogIt( "DSSP WARNING: %s: Could not retrieve file '%s'. (%s)\n" % (pdb_id.upper(), dssp_filename, msg), log_file, -2)
                return []
        else:
            return []

    #Read DSSP File
    content = []

    try:
        with open( dssp_path) as f:
            content = f.readlines()
    except IOError as ex:
        LogIt( "DSSP ERROR: %s: Unable to read DSSP file: '%s' Error: %s\n" % ( pdb_id.upper(), dssp_path, ex.strerror), log_file, -2)
    except Exception as ex:
        LogIt( "DSSP ERROR: %s: Could not read file '%s'. (%s)\n" % (pdb_id.upper(), dssp_filename, str( ex)), log_file, -2)

    if len( content) == 0:
    	LogIt( "DSSP ERROR: %s: File '%s' has no content.\n" % ( pdb_id.upper(), dssp_path), log_file, -2)

    return content

#Returns seqwindow of DSSP output around selected residue
#Arg: aDSSP_Content list of lines in a dssp file
def ProcessDSSP( pdb_id, aDSSP_Content, get_chain, get_resnum, seqwindow=10, log_file=None, aFilePath=""):


    get_resnum = get_resnum.strip() #Can include spaces and an insertion code

    line_num = 0

    #no odd-size windows allowed
    halfwindow = (seqwindow+1) / 2

    #returns a "?" in a list for each residue whos DSSP assesment was not found
    #retval = (["?"] * (halfwindow*2+1))
    #retval[ halfwindow] = '*'

    get_residue_index = -1
    assesments = []
    lines_to_go = halfwindow

    #print "Looking for: '%s'" % get_resnum

    #Read DSSP File
    try:

        header_needed = True
        for line in aDSSP_Content:

            if lines_to_go <= 0: break

            line_num += 1

            if header_needed:
                line = line.strip()
                if len( line) and line[ 0] == "#":
                    header_needed = False
                continue

            # #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
            #    1   70 A R              0   0  165      0, 0.0     9,-0.2     0, 0.0    55,-0.1   0.000 360.0 360.0 360.0 138.8  300.3   -7.1   -2.6
            #   86        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0
            # 1625  736 B P     >  -     0   0   38      0, 0.0     4,-2.7     0, 0.0     3,-0.4  -0.424  29.7-120.3 -57.6 135.9   32.2  -49.2   10.1
            # 1626  736CB S  H  > S+     0   0   94      2,-0.2     4,-2.7     1,-0.2     5,-0.2   0.853 110.6  60.2 -46.4 -38.0   30.6  -47.3    7.2

            chain = line[ 11:12].strip()
            #print "chain: '%s', resnum '%s'" % (line[ 11:12], line[5:11])
            if get_chain != chain: #Wrong chain
                if get_residue_index == -1: continue #Residue has not been found yet, look further
                else: break                          #Residue has been found and chain has changed, stop

            if line[ 13] == "!": continue

            #window_lines.append( line)
            try:
                res_int = int( line[5:10])
            except:
                sys.stderr.write( "DSSP ERROR: Invalid integer '%s' on line %i in file %s.\n\n" % (line[5:10], line_num, aFilePath))
                raise

            dssp_assesment = line[14:17].strip()
            if not len( dssp_assesment): dssp_assesment = "_" #Using underscore instead of a space
            assesments.append( dssp_assesment)

            #Keep assesments at the size of a halfwindow before targeted residue has been found
            if get_residue_index == -1 and len( assesments) > halfwindow: del assesments[ 0]

            resnum = line[5:11].strip() #includes insertion codes and empty spaces

            if get_residue_index == -1:
                if resnum == get_resnum:
                    get_residue_index = len( assesments)-1 #This was the residue we were looking for
            else:
                lines_to_go -= 1

    except Exception as ex:
        sys.stderr.write( "ERROR: %s: Error reading DSSP file: '%s'.\n" % (pdb_id.upper(), aFilePath))
        template =        "       An exception of type {0} occured. Arguments: {1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message + "\n")
        retval = ["Unknown"]
        raise

    if get_residue_index == -1:
        if header_needed == False: #File had a header but no correct residue, files that have only CAs (DSSP not applicable) have no header
            sys.stderr.write( "DSSP ERROR: Residue '%s' in chain '%s' not found in DSSP file %s.\n\n" % (get_resnum, get_chain, aFilePath))
        retval = ["Unknown"]


    #print "DSSP: %s" % dssp_assesment
    #Set

    #Fill to seqwindow size with "?" if necessary
    retval = (["?"]* (max( 0, halfwindow - ( get_residue_index+1)))) + assesments + (["?"]* (max( 0, halfwindow - ( len( assesments)-( get_residue_index+1)))))
    retval = retval[ 0:halfwindow] + ["*"] + retval[ halfwindow:]

    #if resnum < get_resnum and resnum + halfwindow > get_resnum:
    #    retval[ halfwindow - (get_resnum - resnum) - 1] = dssp_assesment
    #elif resnum == get_resnum:
    #    retval[ halfwindow-1] = dssp_assesment
    #elif resnum == get_resnum+1:
    #    retval[ halfwindow+1] = dssp_assesment
    #elif resnum > get_resnum and resnum - halfwindow <= get_resnum:
    #    retval[ halfwindow + (resnum - get_resnum)] = dssp_assesment

    #aa = line[12:14]
    #retval[ get_resnums.index( resnum)] = dssp

    return retval



def main():

    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #print 'Argument List:', str(sys.argv)

    if len( sys.argv) < 2:
        sys.stderr.write( "ERROR: Please specify an input pdb file.\n")
        sys.exit( -1)

    inputFile = sys.argv[ 1]
    if not os.path.isfile( inputFile):
        sys.stderr.write( "ERROR: Input file '%s' not found.\n")
        sys.exit( -1)

    outputFile = None
    if len( sys.argv) > 2:
        outputFile = sys.argv[ 2]

    if not outputFile or len( outputFile) == 0: outputFile = inputFile + ".dssp"


    SetDSSP_Executable( "../../DSSP/dssp-2.0.4-win32.exe")

    fetch_retval = RunDSSP( inputFile, outputFile, aOverWriteExisting=True, aLogFile=None )

    if fetch_retval == 0:
        fetch_success = True
        LogIt( "DSSP INFO: Calculated DSSP for structure '%s' to file '%s'" % (inputFile, outputFile), None, 1)
    elif fetch_retval == -50:
        #DSSP executable set but not found
        print "ERROR: DSSP Executable set to '%s' but not found." % GetDSSP_Executable()
        sys.exit( -50)
    else:
        LogIt( "DSSP INFO: Failed to calculate DSSP for file '%s' (Only backbone available?).\n" % outputFile, None, 2)


if __name__ == "__main__":
    main()