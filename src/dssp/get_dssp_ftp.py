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


import ftplib
import os
import socket
import sys
import urllib

from subprocess import call

sys.path.append(  "%s/../utils" % os.path.realpath(__file__))
from fileops import *

HOST = 'ftp.cmbi.ru.nl'
DIRN = '/pub/molbio/data/dssp'


def RunDSSP( aExe, aPdbFile, aOutputfile, aOverWriteExisting=True):

    retval = call( aExe + " -i %s > %s" % (aPdbFile, aOutputfile), shell=True)
    #print "DSSPRETVAL:", return_code
    return retval



def FetchFTP_File( files, output_dir, host=HOST, dir=DIRN):

    output_dir = InsertFolderSeparator( output_dir)
    retvals = []

    if type( files) == str:
        arr = []
        arr.append( files)
        files = arr

    if len( files) < 1:
        print "Nothing to fetch."
        return retvals


    try:
        ftp = ftplib.FTP(HOST)
    except (socket.error, socket.gaierror), e:
        print 'ERROR: cannot reach "%s"' % HOST
        return retvals
    print '*** Connected to host "%s"' % HOST

    try:
        ftp.login()
    except ftplib.error_perm:
        print 'ERROR: cannot login anonymously'
        ftp.quit()
        return retvals
    print '*** Logged in as "anonymous"'

    try:
        ftp.cwd(DIRN)
    except ftplib.error_perm:
        print 'ERROR: cannot CD to "%s"' % DIRN
        ftp.quit()
        return retvals
    print '*** Changed to "%s" folder' % DIRN

    errors = 0
    max_errors = 1
    not_found = 0
    fetched = 0
    skipped = 0

    #files = []

    #try:
    #    files = ftp.nlst()
    #except ftplib.error_perm, resp:
    #    if str(resp) == "550 No files found":
    #        print "No files in this directory"
    #    else:
    #        raise
    #
    #for f in files:
    #    print f
    #
    #sys.exit(0)

    for f in files:

        f = ChangeToFolder( f.lower(), "")

        if f.find( "pdb") == 0:
            f = f[ 3:]

        if f.find( ".ent") >= 0:
            f = f[ :-4]

        if f.find( ".pdb") >= 0:
            f = f[ :-4]

        if len( f) == 4:
            f = f + ".dssp"

        output_file = ChangeToFolder( f, output_dir)
        if os.path.isfile( output_file) and int( os.stat( output_file).st_size) > 0:
            #print "File '%s' already exists. Skipping." % output_file
            skipped += 1
            continue

        try:
            print "Fetching FTP file: %s... " % f,
            with open(output_file, 'wb') as of:
                ftp.retrbinary('RETR %s' % f, of.write)
        except ftplib.error_perm:
            print 'ERROR: Cannot read "%s"' % ( f)
            #os.unlink( output_file)
            not_found += 1
            #errors += 1
            if errors >= max_errors:
                break
        except:
            print 'ERROR: Cannot write "%s"' % ( output_file)
            #os.unlink( output_file)
            errors += 1
            if errors >= max_errors:
                break

        else:
            retvals.append( output_file)
            fetched += 1
            print 'Done.  %i/%i (%.2f%%)' % (not_found+skipped+fetched,len( files), float(not_found+skipped+fetched) / len( files) *100.0)


    if errors >= max_errors:
        print "Too many errors encountered."

    print ""
    print ""
    print "All Done."
    print ""
    print "Retrieved %i files." % len( retvals)
    print "Skipped %i files." % skipped
    print "%i were not found." % not_found

    ftp.quit()

    return retvals



def GetDSSP( pdb_id, DSSP_folder, get_chain, get_resnum, seqwindow = 10 ):

    try:
        get_resnum = int( get_resnum)
    except:
        sys.stderr.write( "ERROR: %s: resnum '%s' is not an integer.\n" % (pdb_id.upper(), get_resnum))

    #no odd-size windows allowed
    halfwindow = (seqwindow+1) / 2

    #returns a "?" in a list for each residue whos DSSP assesment was not found
    retval = (["?"] * (halfwindow*2+1))
    retval[ halfwindow] = '*'


    dssp_filename = pdb_id.lower() + ".dssp"
    dssp_path = InsertFolderSeparator( DSSP_folder) + dssp_filename


    #File does not exist, fetch it
    if not os.path.isfile( dssp_path) or int( os.stat( dssp_path).st_size) <= 0:

        try:
            urllib.urlretrieve("ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/" + dssp_filename, dssp_path)
            print "Fetched: %s" % dssp_filename
        except Exception as ex:
            sys.stderr.write( "ERROR: %s: Error fetching DSSP file: '%s'.\n" % (pdb_id.upper(), dssp_filename))
            template =        "       An exception of type {0} occured. Arguments: {1!r}"
            message = template.format(type(ex).__name__, ex.args)
            sys.stderr.write( message + "\n")
            return ["Unknown"] #Return


    #Read DSSP File
    try:
        with open( dssp_path, "r" ) as f:

            header = True
            for line in f:

                if header:
                    line = line.strip()
                    if len( line) and line[ 0] == "#":
                        header = False
                else:

                    #  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
                    #    1   70 A R              0   0  165      0, 0.0     9,-0.2     0, 0.0    55,-0.1   0.000 360.0 360.0 360.0 138.8  300.3   -7.1   -2.6

                    chain = line[ 10:12].strip()
                    if get_chain != chain: continue #Wrong chain

                    resnum = int( line[5:10])
                    dssp_assesment = line[14:17].strip()

                    if not len( dssp_assesment): dssp_assesment = "_" #Using underscore instead of a space
                    #print "DSSP: %s" % dssp_assesment


                    if resnum < get_resnum and resnum + halfwindow > get_resnum:
                        retval[ halfwindow - (get_resnum - resnum) - 1] = dssp_assesment
                    elif resnum == get_resnum:
                        retval[ halfwindow-1] = dssp_assesment
                    elif resnum == get_resnum+1:
                        retval[ halfwindow+1] = dssp_assesment
                    elif resnum > get_resnum and resnum - halfwindow <= get_resnum:
                        retval[ halfwindow + (resnum - get_resnum)] = dssp_assesment

                    #aa = line[12:14]
                    #retval[ get_resnums.index( resnum)] = dssp

    except IOError as e:
        sys.stderr.write( "ERROR: %s: Unable to read DSSP file: '%s' Error: %s\n" % ( pdb_id.upper(), dssp_path, e.strerror))
        retval = ["Unknown"]
    except Exception as ex:
        sys.stderr.write( "ERROR: %s: Error reading DSSP file: '%s'.\n" % (pdb_id.upper(), dssp_path))
        template =        "       An exception of type {0} occured. Arguments: {1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message + "\n")
        retval = ["Unknown"]

    return retval



def main():

    #Example usage:
    RunDSSP( "../../dssp/dssp-2.0.4-win32.exe", "../../PDB/pdb1a0j.ent", "./test.dssp", False )


if __name__ == "__main__":
  main()
