# This file is part of the StructureMapper algorithm.
# Please cite the authors if you find this software useful
#
# https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361
 
# MIT License
#
# Copyright 2018 user Nurminen and Vesa P. Hytönen
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
from datetime import datetime

from sis_globals import *


def JobHashToDir( aJobHash):
    global RESULT_FOLDER
    return InsertFolderSeparator( RESULT_FOLDER + aJobHash)

def LogError( aStr, aPrint=2, aLock=None):
    global ERROR_LOG
    try:
        aStr = datetime.now().strftime('[%Y-%m-%d %H:%M] ') + aStr
        LogIt( aStr, aFilename=ERROR_LOG, print_it=aPrint, lock=aLock)
    except Exception as ex:
        sys.stderr.write( "ERROR: Exception when writing error log.\n")
        sys.stderr.write( "       An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))

def LogServer( aStr, aPrint=1, aLock=None):
    global JOB_LOG
    try:
        aStr = datetime.now().strftime('[%Y-%m-%d %H:%M] ') + aStr
        LogIt( aStr, aFilename=JOB_LOG, print_it=aPrint, lock=aLock)
    except Exception as ex:
        sys.stderr.write( "ERROR: Exception when writing server log.\n")
        sys.stderr.write( "       An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))

#print_it: 1=stdout, 2=stderr, 0=Do not print to screen, -1=Print if no file specified (stdout), -2=Print if no file specified (stderr),
def LogIt( aStr, aFilename=None, print_it=1, lock=None):

    if not len( aStr): return
    elif aStr[ -1] != "\n": aStr += "\n"

    #Print always
    if( print_it == 1): sys.stdout.write( aStr)
    elif( print_it == 2): sys.stderr.write( aStr)

    #If it's a filehandle
    if aFilename and type( aFilename) == file:
        if lock: lock.acquire()
        aFilename.write( aStr)
        if lock: lock.release()
        return #RETURN

    if aFilename == None or len( aFilename) == 0:
        #Print if there is nowhere else to print
        if( print_it == -1): sys.stdout.write( aStr)
        elif( print_it == -2): sys.stderr.write( aStr)
        return #RETURN

    if lock: lock.acquire()

    #Write to file
    try:
        with open( aFilename, 'a+') as log_file:
            log_file.write( aStr)
    except IOError as e:
        sys.stderr.write( "WARNING: Could not write log file: '%s'." % ( aFilename))
    except:
        sys.stderr.write( "WARNING: Could not write log file (general): '%s'." % ( aFilename))

    if lock: lock.release()



def InsertFolderSeparator( folder):

    if len( folder) == 0: return folder

    folder_sep = "/"
    if folder.find( "\\") >= 0: folder_sep = "\\"

    #Make sure folder ends with a folder separator character
    if folder[ -1] != folder_sep:
        folder += folder_sep

    return SystemSpecificFolderSeparators( folder)

def SystemSpecificFolderSeparators( aPath):

    if os.name == 'nt': #windows
        aPath = aPath.replace( "/", "\\")
        return aPath.replace( "\\\\", "\\")

    aPath = aPath.replace( "\\", "/")
    return aPath.replace( "//", "/")

def DjangoConfigPaths():

    retvals = {}

    if os.name == 'nt': #windows
        #DEV environment
        retvals['PRECALC_ASA_DIR'] = "I:\\user\\PDB\\ASA"
        retvals['PRECALC_ASA_SC_DIR'] = "I:\\user\\PDB\\ASA_SC"
        retvals['PDB_DIR'] = "I:\\user\\PDB"
        retvals['BLAST_DB'] = "I:\\user\\Blast_DB"
        retvals['DEFAULT_DSSP_DIR'] = "I:\\user\\dssp\\files"
        retvals['DEFAULT_DSSP'] = "I:\\user\\Proteases\\dssp\\dssp-2.0.4-win32.exe"
        retvals['DEFAULT_IUPRED'] = "I:\\user\\Proteases\\disorder\\iupred\\iupred.exe"
    else:
        #Online
        retvals['PRECALC_ASA_DIR'] = "/home/user/sis/PDB/ASA"
        retvals['PRECALC_ASA_SC_DIR'] = "/home/user/sis/PDB/ASA_SC"
        retvals['PDB_DIR'] = "/home/user/sis/PDB"
        retvals['BLAST_DIR'] = "/home/user/sis/BLAST/DB"
        retvals['BLAST_DB'] = "/home/user/sis/BLAST/DB"
        retvals['DEFAULT_BLASTP'] = "/home/user/sis/BLAST/ncbi-blast-2.5.0+/bin"
        retvals['DEFAULT_DSSP_DIR'] = "/home/user/sis/DSSP/DB"
        retvals['DEFAULT_DSSP'] = "/home/user/sis/DSSP/dssp-2.0.4-linux-i386"
        retvals['DEFAULT_IUPRED'] = "/home/user/sis/iupred/iupred"

    return retvals

def TextToHtml():
    pass


#Returs an array of record IDs
def GetRecordACC( aIdentifier, aWithPos=True ):

    """"Given a SeqRecord identifier string, return the accession number as a string.

    e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "2765613"
    e.g. "gi|2765613_Pos232|emb|Z78488.1|PTZ78488" -> "2765613_Pos232"

    If there are

    """
    aIdentifier = aIdentifier.strip()
    if aIdentifier[ 0] == ">": aIdentifier = aIdentifier[ 1:]
    cols = aIdentifier.split("|")

    acc = ""

    if len( cols) < 2: acc = aIdentifier[ 0:10]
    elif cols[ 1][0:4] == "POI:": acc = cols[ 0]
    else: acc = cols[ 1]

    if not aWithPos: acc = StripPosTag( acc)

    if len( acc) == 0:
        sys.stderr.write( "ERROR: Could not find identifier for seq '%s'.\n") % aIdentifier
        sys.exit( -32)

    return acc