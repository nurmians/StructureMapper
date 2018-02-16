#!/usr/bin/python
# coding=UTF-8
# -*- coding: UTF-8 -*-

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
import shutil
import time
from multiprocessing import Process

from Bio.PDB.PDBList import PDBList

from fileops import InsertFolderSeparator


def FetchPDBFile( PDB_Code, output_folder="", force_fetch=False, obsoletes=[], verbose=False, aAsync=True):

    PDB_Code = PDB_Code.upper()
    output_folder = InsertFolderSeparator( output_folder)

    if PDB_Code in obsoletes:
        print "PDB Code '%s' is obsolete." % PDB_Code
        return ""

    if not os.path.isdir( output_folder):
        sys.stderr.write( "'%s' does not specify a valid directory.\n" % output_folder)
        return ""

    filename = output_folder + "pdb" + PDB_Code.lower() +".ent"

    if not os.path.isfile( filename):
        force_fetch = True
    elif int( os.stat( filename).st_size) <= 0:
        force_fetch = True
        try:
            os.remove( filename)
        except:
            pass


    if not force_fetch:
        try:
            with open( filename, "r") as f:
                if verbose: print "File '%s' already exists." % filename
                return filename
        except IOError: pass

    retval = -100

    if aAsync == True:
        retval = RetrievePDBFile( PDB_Code, output_folder)
    else:
        #Returns file name on success
        retval = DoRetrievePDBFile( PDB_Code, output_folder) #Synchronous
        if len( retval) > 0: retval = 0
        else: retval = -4

    if retval != 0:
        sys.stderr.write( "ERROR: Could not fetch file with PDB code '%s'. (%s)\n" % (PDB_Code, str( retval)))

    return retval


def FetchObsoletes( write_to_file=True, folder=""):

    folder = InsertFolderSeparator( folder)
    obfile = folder + "obsoletes.txt"

    obs = []

    read_ok = True

    if not os.path.isfile( obfile):
        #Obsoletes file does not exist
        is_ok, obs = DoFetchObsoletes( obfile)
        if not is_ok:
            sys.stderr.write("No information on obsolete structures could be obtained. '%s'.\n" % obfile)
            sys.exit(-1) #EXIT
    else:

        fileCreation = os.path.getmtime( obfile) #modification time
        #fileCreation = os.path.getctime( obfile) #creating time, does not change when truncated
        now = time.time()
        oneday_ago = now - 1*60*60*24*1 #24h


        if fileCreation < oneday_ago:
            #Obsoletes file is more than 24h old (refresh)
            sys.stderr.write("Obsoletes file '%s' needs to be refreshed.\n" % obfile)
            is_ok, obs = DoFetchObsoletes( obfile)
            if not is_ok:
                sys.stderr.write("Existing obsoletes file could not be refreshed. '%s'.\n" % obfile)

    #Use existing or freshly fetched file
    try:
        with open( obfile, 'r') as f:
            obs = f.read().splitlines()
    except IOError:
        sys.stderr.write("Could not read obsoletes file: '%s'.\n" % obfil)
        sys.exit(-1) #EXIT


    return obs


def DoFetchObsoletes( filename ):

    obs = []
    fetch_ok = True

    try:
        pdblist = PDBList();
        sys.stdout.write("INFO: Fetching obsolete structure information online...\n")
        obs = pdblist.get_all_obsolete();
    except:
        fetch_ok = False
        sys.stderr.write("[FAILED]\nUnable to fetch obsolete structures information online.\n")

    if fetch_ok:
        sys.stdout.write("[OK].\n")

        try:
            with open( filename, 'w') as f:
                for ob in obs:
                    f.write("%s\n" % ob)
            fetch_ok = True

        except IOError:
            sys.stderr.write("ERROR: Could not write obsoletes into file: '%s'.\n" % filename)
            fetch_ok = False

    return fetch_ok, obs


#Async
def RetrievePDBFile( aPDB_Code, aFolder ):

    filename = aFolder + "pdb" + aPDB_Code.lower() +".ent"
    backup_made = False

    if os.path.isfile( filename):
        try:
            shutil.move( filename, filename + ".bak")
            backup_made = True
        except:
            sys.stderr.write( "ERROR: Could not create a backup of existing file '%s'\n" % filename)
            return -3

    timeout_retries = 5
    first = True

    while timeout_retries > 0 or first:

        first = False

        p = Process( target=DoRetrievePDBFile, args=(aPDB_Code, aFolder)) #multiprocessing
        p.start()

        # Wait for 20 seconds or until process finishes
        p.join(20)

        # If thread is still active
        if p.is_alive():
            sys.stderr.write( "ERROR: Retrieve timeout.\n")

            # Terminate
            p.terminate()
            p.join()
            #Timed out
            timeout_retries -= 1
            if timeout_retries > 0: sys.stderr.write( "INFO: Retrying download in %i seconds...\n" % 30)
            else: sys.stderr.write( "INFO: No timeout retries left.\n")
            continue

        break


    #File exists [OK]
    if os.path.isfile( filename) and int( os.stat( filename).st_size) > 0:
        return 0

    if backup_made:
        try:
            shutil.move( filename + ".bak", filename)
        except:
            sys.stderr.write( "ERROR: Could not replace original file '%s' with created backup file '%s'\n" % (filename, (filename + ".bak")))

    #File does not exist
    return -1


USE_ALT_PDB_SERVER = False

def DoRetrievePDBFile( aPDB_Code, aFolder ):

    global USE_ALT_PDB_SERVER

    done = False
    errors_before_quit = 20
    seconds_between_retries = 30
    fetchedfile = ""
    alt_server = "http://www.rcsb.org/pdb/files/"

    while done == False:

        pdblist = None
        if USE_ALT_PDB_SERVER: pdblist = PDBList( server=alt_server)
        else: pdblist = PDBList()
        #pdblist = PDBList( server='ftp://ftp.wwpdb.org')
        #server = 'ftp://ftp.rcsb.org'
        #server = "ftp.ebi.ac.uk/pub/databases/pdb/"

        try:
            #http://biopython.org/DIST/docs/api/Bio.PDB.PDBList%27-pysrc.html
            #fetchedfile = pdblist.retrieve_pdb_file( pdb_code=aPDB_Code, pdir=aFolder, file_format="pdb", obsolete=False)
            fetchedfile = pdblist.retrieve_pdb_file( pdb_code=aPDB_Code, pdir=aFolder, file_format="pdb", obsolete=False)
            done = True
            if fetchedfile and len( fetchedfile) and (fetchedfile.find( ".ent") > 0 or fetchedfile.find( ".pdb") > 0):
                #print "Structure fetched, PDB code: " + aPDB_Code
                print "INFO: Structure " + aPDB_Code + " fetched. [OK]"
                #io = PDBIO()
                #io.set_structure( s)
                #io.save( filename)
            else:
                print "WARNING: Fetch failed [FAIL]"

        except IOError as ex:
            sys.stderr.write( "WARNING: Could not download structure {0}. An exception of type {1} occured.\n       Arguments: {2!r}\n".format( aPDB_Code, type( ex).__name__, ex.args))
            sys.stderr.write( "INFO: Retrying connection in %i seconds...\n" % seconds_between_retries)


            for a in ex.args:
                #Downloading too many structures too fast?
                if str( a).lower().find("too many") >= 0:
                    seconds_between_retries += 10
                    break
                if str( a).lower().find("No such file") >= 0:
                    #No need to retry
                    return fetchedfile
                if str( a).lower().find("did not properly respond") >= 0:
                    #No need to retry
                    sys.stderr.write( "INFO: Switching download thread to alternative server '%s'.\n" % alt_server)
                    USE_ALT_PDB_SERVER = True

            time.sleep( seconds_between_retries)
            done = False
            errors_before_quit -= 1
            if errors_before_quit <= 0:
                sys.stderr.write( "ERROR: Failed too many times. Quitting...\n")
                break


    return fetchedfile

