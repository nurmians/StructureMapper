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


import time
import signal, os, sys
import re
import signal

from multiprocessing import Lock, Value, Array

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import is_aa

from fileops import InsertFolderSeparator

#Buffered SyncPrint allows for more controlled printing when
#running parallel processing. When Set to buffered, output is not printed directly to console,
#but to a buffer which can be flushed (printed) at desired points or intervals.
g_SyncPrintBuffer_max = 20000 #buffer size in number of chars
g_SyncPrintBuffered = Value( 'i', 0 ) #YES == 1, NO == 0
g_SyncPrintBuffer = Array( 'c', " "*g_SyncPrintBuffer_max )
g_SyncPrintBuffer.value = ""
g_SyncPrintBuffer_overflow = Value( 'i', 0) #Overflow indicator

#Pass the return value to a thread as an initializer argument
def GetSyncPrintBufferForThread():
    global g_SyncPrintBuffered, g_SyncPrintBuffer, g_SyncPrintBuffer_overflow
    return (g_SyncPrintBuffered, g_SyncPrintBuffer, g_SyncPrintBuffer_overflow)

#In the Thread, call this method with the return value of GetSyncPrintBufferForThread()
def SetSyncPrintBufferInThread( aBuffer):
    global g_SyncPrintBuffered, g_SyncPrintBuffer, g_SyncPrintBuffer_overflow
    g_SyncPrintBuffered = aBuffer[ 0]
    g_SyncPrintBuffer = aBuffer[ 1]
    g_SyncPrintBuffer_overflow = aBuffer[ 2]


def SetSyncPrintBuffered( aSetBuffered=True):
    global g_SyncPrintBuffered, g_SyncPrintBuffer

    SyncPrintFlushBuffer()
    g_SyncPrintBuffered.value = (1 if aSetBuffered else 0)


def IsSyncPrintBuffered():
    global g_SyncPrintBuffered
    return (g_SyncPrintBuffered.value == 1)

def SyncPrintWriteBuffer( aStr, aLock=False, aNewLine=True):
    global g_SyncPrintBuffer, g_SyncPrintBuffer_max, g_SyncPrintBuffer_overflow

    #SyncErrPrint( "Buffering: '%s'" % aStr)
    msg = aStr + ("\n" if aNewLine else "")
    len_msg = len( msg)

    if aLock: aLock.acquire()

    if len( g_SyncPrintBuffer.value) + len_msg <= g_SyncPrintBuffer_max:


        g_SyncPrintBuffer.value += msg
        #print "Buffering: '%s'. " % msg.strip(), aLock
        if aLock: aLock.release()

    elif not g_SyncPrintBuffer_overflow.value:
        g_SyncPrintBuffer_overflow.value = 1
        if aLock: aLock.release()
        SyncErrPrint( "WARNING: SyncPrintBuffer overflow...", aLock)
    else:
        if aLock: aLock.release()

def SyncPrintFlushBuffer( aLock=False, aTitle=None ):

    global g_SyncPrintBuffer, g_SyncPrintBuffer_overflow

    if aLock: aLock.acquire()
    if len( g_SyncPrintBuffer.value) == 0:
        if aLock: aLock.release()
        return #nothing to print

    #sys.stdout.write( "Clearing buffer:\n")
    if aTitle: sys.stdout.write( aTitle)
    sys.stdout.write( g_SyncPrintBuffer.value)
    g_SyncPrintBuffer.value = "" #Clear buffer
    g_SyncPrintBuffer_overflow.value = 0
    if aLock: aLock.release()

def SyncPrint( aStr, aLock=False, aNewLine=True, aBypassBuffer=False):

    if aLock and not aBypassBuffer and IsSyncPrintBuffered():
        SyncPrintWriteBuffer( aStr, aLock, aNewLine)
        return

    nl = "\n" if aNewLine else ""
    if aLock: aLock.acquire()
    sys.stdout.write( aStr + nl)
    if aLock: aLock.release()

def SyncErrPrint( aStr, aLock=False, aNewLine=True):
    nl = "\n" if aNewLine else ""
    if aLock: aLock.acquire()
    sys.stderr.write( aStr + nl)
    if aLock: aLock.release()

def CountRecords( aRecordFile):
    count = 0
    try:
        with open( aRecordFile, "r") as rf:
            for line in rf:
                if line[ 0] == ">": count += 1
    except Exception as ex:
        template = "ERROR: An exception of type {0} occured. Arguments: {1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message + "\n")
        return -1
    return count

#AA -> a
def RealChain( aChain):
    if aChain and len( aChain) == 2:
        return aChain[ 0].lower()
    return aChain

#a --> AA
def FilenameChain( aChain):
    if aChain and len( aChain) and aChain[ 0].upper() != aChain[ 0]:
        return (aChain[ 0]*2).upper()
    return aChain

#Generator
def FastaIter( aFastaFile, aJoinSeqLines=True, aRecordSelector=None):
    """
    given a fasta file. yield tuples of header, sequence
    """
    from itertools import groupby
    rec_index = -1
    fh = open( aFastaFile, "r")
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[ 1] for x in groupby(fh, lambda line: line[ 0] == ">"))
    for grouped in faiter:

        if aRecordSelector != None:
            rec_index += 1
            if aRecordSelector( rec_index) == False:
                faiter.next()
                continue
        # drop the ">"
        header = grouped.next()
        if header[ 0] == ";": continue
        header = header[1:].strip()
        #header = header[1:].strip()

        if aJoinSeqLines: seq = "".join( s.strip() for s in faiter.next())
        else: seq = "\n".join( s.strip() for s in faiter.next())

        #print "Yielding: %s" % header #debug
        yield header, seq
        #if aRecordSelector == None or aRecordSelector( rec_index) == True:
            #yield header, seq

    fh.close()

#Generator
def FastaHeaderIter( fasta_name):

    with open(fasta_name, "r") as fhf:

        for line in fhf:
            #stripped = line.strip()
            #if len( line) > 0 and line.strip()[ 0] == ">":
            if len( line) > 0 and line[ 0] == ">": #Needs to be at the beginning of the line
                yield line.rstrip()


def ParseDescriptionFromFastaHeader( aHeader):

    cols = aHeader.split("|")
    if len( cols) == 1: return aHeader

    acc = GetRecordACC( aHeader, aWithPos=True, aRemoveInvalidChars=False)
    acc_pos = aHeader.find( acc)

    cols = aHeader[ acc_pos+len( acc):].split("|")
    if len( cols) <= 1: return ""

    m = re.search( "(.+?) ..=", cols[ 1])
    if m: return m.group( 1)
    
    return cols[ 1] 



def ParseDescriptionFromFastaHeaderOLD( aHeader):

    cols = aHeader.split("|")

    desc_col = -1

    if hasattr( ParseDescriptionFromFastaHeader, "desc_col"):
        desc_col = ParseDescriptionFromFastaHeader.desc_col

    if desc_col < 0:
        c = -1
        for col in cols:
            c += 1
            m = re.search( "(.+?) ..=", col)
            if m:
                desc_col = c
                break

    if desc_col < 0:
        c = -1
        max_spaces = -1
        max_space_col = -1

        for col in cols:
            c += 1
            spaces = col.count( " ")
            if spaces > max_spaces:
                max_spaces = spaces
                max_space_col = c

        if max_space_col <= 0:
            desc_col = max_space_col

    if desc_col < 0: return aHeader #Could not parse description

    ParseDescriptionFromFastaHeader.desc_col = desc_col
    desc = cols[ desc_col]

    m = re.search( "(.+?) ..=", desc)
    if m: desc = m.group( 1)

    #Remove e.g. IDH1_HUMAN
    m2 = re.search( "(\S+_\S+\s)", desc)
    if m2: desc = desc.replace( m2.group( 1), "")

    return desc

def StripPosTag( aAcc):

    if not aAcc or not len( aAcc): return ""

    pospos = aAcc.find( "_Pos")
    if pospos >= 0: return aAcc[:pospos]
    return aAcc

#returns POI N-terminal aa position as str
def GetRecordPOI( aDescription):

    retval = -1

    matchObj = re.search( "POI:(\d+)", aDescription)
    if matchObj: retval = int( matchObj.group( 1))

    return retval

#Returs an array of record IDs
def GetRecordACC( aIdentifier, aWithPos=True, aRemoveInvalidChars=True ):

    """"Given a SeqRecord identifier string, return the accession number as a string.

    e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "2765613"
    e.g. "gi|2765613_Pos232|emb|Z78488.1|PTZ78488" -> "2765613_Pos232"

    This method needs to handle fasta headers correctly when they are first inputted
    and after they have been processed by fasta_formatter

    """
    aIdentifier = aIdentifier.strip()
    if aIdentifier[ 0] == ">": aIdentifier = aIdentifier[ 1:]
    cols = map( str.strip, aIdentifier.split("|"))

    acc = ""
    #print "COLS:", cols

    if len( cols) < 2: acc = aIdentifier[ 0:10]
    elif len( cols) == 2 and len( cols[ 0]) > 2: acc = cols[ 0]
    elif cols[ 1].startswith( "POI:"): acc = cols[ 0]
    else: acc = cols[ 1]

    if not aWithPos: acc = StripPosTag( acc)

    if len( acc) == 0:
        sys.stderr.write( "ERROR: Could not find identifier for seq '%s'.\n" % aIdentifier)
        sys.exit( -32)

    if aRemoveInvalidChars: acc = ''.join(c if c not in ";,: \t" else "_" for c in acc)
    return acc

def GetRecordIDsFromFile( aFile, aWithPos=True ):

    rec_ids = []
    seq_headers = FastaHeaderIter( aFile)
    for header in seq_headers:
        rec_ids.append( GetRecordACC( header, aWithPos))
    return rec_ids


#Returs an array of record IDs
def GetRecordIDs( records, aWithPos=True ):

    retarr = []
    ids = {}

    for rec in records:

        acc = GetRecordACC( rec.description, aWithPos)
        if aWithPos and acc in ids:
            sys.stderr.write( "ERROR: ID '%s' is not unique within the dataset. Sequence '%s'. Exitting.\n" % (acc, rec.description))
            sys.exit( -45)

        ids[ acc] = True
        retarr.append( acc)

        #cols = rec.description.split("|")
        #found = False
        #for col in cols:
        #    if col.isalpha(): continue
        #    else:
        #
        #        newid = AddUnique( col, ids)
        #        retarr.append( newid)
        #        found = True
        #        break
        #
        #if not found:
        #      retarr.append( AddUnique( rec.description[0:10].translate( None, "|"), ids))

    return retarr


#FASTA_HEADER_GENE_PATTERN = re.compile( "GN\s*\=\s*([A-Z0-9-]+)") #old: "GN\s*\=\s*(\S+)(?:\||\s|$)"
FASTA_HEADER_GENE_PATTERN = re.compile( "GN\s*\=\s*([^\s\|$]+)") #old: "GN\s*\=\s*(\S+)(?:\||\s|$)"

#Creates a dictionary of key:entry_id value:gene
#For sequences that have "GN=*" specified in their description
def GetGenesForRecords( records ):

    global FASTA_HEADER_GENE_PATTERN

    retval = {}

    for rec in records:

        #cols = rec.description.split("|")

        matchObj = FASTA_HEADER_GENE_PATTERN.search( rec.description, re.I)

        if matchObj:
            retval[ GetRecordACC( rec.id)] = matchObj.group( 1)
            #print "DBG: Gene in '%s' => '%s'" % (rec.description, matchObj.group( 1))
        else:
            pass
            #print "DBG: No gene found in '%s'" % rec.description

    return retval


def GetGeneFromFastaHeader( aHeader):

    global FASTA_HEADER_GENE_PATTERN

    matchObj = FASTA_HEADER_GENE_PATTERN.search( aHeader, re.I)
    if matchObj: return matchObj.group( 1)
    return ""



#Argument is array of sequence IDs
def GetSequenceInfoFromFastaHeaders( aRecordFile ):

    #records = list( SeqIO.parse( aRecordFile, "fasta"))
    seq_dict = dict()

    print "Reading record information from file '%s'...." % aRecordFile

    re_seqid = re.compile( "\|(.+)\|")
    re_gene = re.compile( "GN=\s*(\S+)")
    re_desc = re.compile( "\|\S+\s+(.+?)\s+\S+=")



    with open( aRecordFile, "r") as rf:

        linenum = 0

        for line in rf:
            linenum += 1
            if line[ 0] != ">": continue

            m = re_seqid.search( line)
            if not m:
                sys.stderr.write( "WARNING: No Sequence id found for header on line %i: '%s'.\n" % (linenum, line) )
                continue

            seq_id = m.group( 1)
            gene = ""
            m = re_gene.search( line)
            if m: gene = m.group( 1)

            desc = ""
            m = re_desc.search( line)
            if m: desc = m.group( 1)

            if seq_id in seq_dict:
                sys.stderr.write( "WARNING: Seq_ID '%s' found multiple times. Line %i ignored.\n" % (seq_id, linenum))
                continue

            if seq_id:
                seq_dict[ seq_id] = { 'GENE' : gene, 'DESC' : desc }


    print "Found %i sequence headers." % len( seq_dict)

    return seq_dict

def AddUnique( aId, aUsed ):

    test = aId
    i = 2

    while test in aUsed:
        test = aId + "_" + str( i)
        i += 1

    aUsed[ test] = True
    return test


AA_THREE =  {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'HSE': 'H', 'HSD': 'H', 'HSP': 'H',
             'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER': '*', 'ALA': 'A',
             'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', '---': '-'} # 'XAA': 'X'



def AAcode_3_to_1( seq):

    global AA_THREE

    '''Turn a three letter protein into a one letter protein.

    The 3 letter code can be upper, lower, or any mix of cases
    The seq input length should be a factor of 3 or else results
    in an error

    >>>AAcode_3_to_1('METHISARGARGMET')
    >>>'MHRRM'

    '''
    if len(seq) % 3 != 0:
        sys.stderr.write( "ERROR: Sequence is not a factor of 3 in length.\n")
        return None

    upper_seq = seq.upper()
    single_seq = ''

    for i in range( len( upper_seq) / 3):
        res = upper_seq[3*i:3*i+3]
        if res in AA_THREE:
            single_seq += AA_THREE[ res]
        else:
            single_seq += "-"

    return single_seq

def PDBFileFromAccCode( code):
    return "pdb" + code.lower() + ".ent"

def GetPDBFileTitle( pdb_filename ):

    title = ""
    found = False
    try:
        with open( pdb_filename, 'r') as f:
            for line in f:
                if line.find( "TITLE") == 0:
                    if found: title += " " #add space
                    title += line[ 10:].strip()
                    found = True #Title can be on multiple rows
                elif found: break
    except:
        pass

    return title

#Get sequence length from PDB file
#returns (length,sequence)
def GetSeqLengthFromPDBFile( structure_id, filename, chain_id, include_gaps = False, quiet = False ):

    sequence = ""

    ppar = PDBParser( PERMISSIVE=1, QUIET=True)

    if not quiet: print "Reading file: " + filename

    try:
        with open( filename): pass
    except IOError:
        print "FILE ERROR: File '%s' not found." % filename
        return -1, sequence


    structure = ppar.get_structure(structure_id, filename )

    #Error getting structure
    if len( structure) == 0: return  -1, sequence

    model = structure[ 0]

    all_chains = []

    for chain in model:
        all_chains.append( chain.id)

    if chain_id not in all_chains and len( chain_id) > 1 and not chain_id.isdigit():
        #ZZ -> z
        chain_id = chain_id[0:1].lower()

    chain_in_model = chain_id in model

    if chain_in_model and len( model[ chain_id]):
        chain = model[ chain_id]

        seq_length = 0
        first = 999999
        last = 0
        previous_num = -1

        for residue in chain:
            if is_aa( residue):

                res_num = -1
                try:
                    res_num = int( residue.id[ 1])
                except ValueError:
                    sys.stderr.write("WARNING: Invalid residue number '%s' in file '%s'\n" % (residue.id[ 1], filename))
                    continue

                if include_gaps and previous_num > 0 and res_num > previous_num+1:
                    sequence += ("-"*(3*((res_num-previous_num)-1))) #insert gaps if there is a gap in numbering

                previous_num = res_num
                sequence += residue.get_resname()

                #print "resnum %s %i" % (residue.resname, res_num)
                first = min( first, res_num)
                last = max( last, res_num)

        seq_length = last - first + 1

        #print structure_id + " sequence: " + sequence
        return seq_length, AAcode_3_to_1( sequence)

    #Error finding chain in structure
    return  -1, sequence

#errors 0=ignore, 1=report, 2=raise exception
def TypeCast( value, type_str, errors=0):

    try:
        if type_str == "s":
            return str( value)
        elif type_str == "i":
            return int( value)
        elif type_str == "f":
            return float( value)

    except ValueError:
        if errors == 1: sys.stderr.write( "value: %s is not of type %s\n" % (value, type_str))
        if errors == 2: raise


    return value


def RealSeqID( sid):

    sid = sid.strip()
    return re.split(':\d+$', sid)[ 0]


def printDictionary( dict):

    print ""

    for k in sorted( dict.iterkeys()):
        print k,
        print ":",
        print dict[ k]

    print ""


def ClearLog( filename):

    try: os.remove( filename)
    except: pass

#print_it: 1=stdout, 2=stderr, 0=Do not print to screen, -1=Print only if no aFilename == None specified (stdout), -2=Print only if no aFilename == None specified (stderr),
def LogIt( aStr, aFilename=None, print_it=1, lock=None):

    if not len( aStr): return
    elif aStr[ -1] != "\n": aStr += "\n"

    #Print always
    if( print_it == 1): SyncPrint( aStr, lock, False) #sys.stdout.write( aStr)
    elif( print_it == 2): SyncErrPrint( aStr, lock, False) #sys.stderr.write( aStr)

    #If it's a filehandle
    if aFilename and type( aFilename) == file:
        if lock: lock.acquire()
        aFilename.write( aStr)
        if lock: lock.release()
        return #RETURN

    if aFilename == None:# or len( aFilename) == 0:
        #Print if there is nowhere else to print
        if( print_it == -1): SyncPrint( aStr, lock, False)
        elif( print_it == -2): SyncErrPrint( aStr, lock, False)
        return #RETURN

    #Not logged
    if len( aFilename) == 0: return

    if lock: lock.acquire()

    #Write to file
    try:
        with open( aFilename, 'a') as log_file:
            log_file.write( aStr)
    except IOError as e:
        SyncErrPrint( "WARNING: Could not write log file: '%s'." % ( aFilename), lock)
    except:
        SyncErrPrint( "WARNING: Could not write log file (general): '%s'." % ( aFilename), lock)

    if lock: lock.release()


class DelayedKeyboardInterrupt( object):
    def __enter__( self):
        self.signal_received = False
        self.old_handler = signal.getsignal( signal.SIGINT)
        signal.signal( signal.SIGINT, self.handler)

    def handler( self, sig, frame):
        self.signal_received = (sig, frame)
        sys.stderr.write( 'KeyboardInterrupt received.\nPlease wait for critical process to finish...\n')

    def __exit__( self, type, value, traceback):
        signal.signal( signal.SIGINT, self.old_handler)
        if self.signal_received:
            self.old_handler( *self.signal_received)


def WrapSeq( aSeq, aLen=60, aEndWithNewLine=True): #chunkify

    #much slow
    #from textwrap import wrap
    #return "\n".join( wrap( aSeq.replace("\n", ""), aLen)) + ("\n" if aEndWithNewLine else "")

    aSeq = aSeq.replace('\n', '')
    return "\n".join( (aSeq[0+i:aLen+i] for i in range(0, len( aSeq), aLen))) + ("\n" if aEndWithNewLine else "")

