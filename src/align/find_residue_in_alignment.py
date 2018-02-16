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
import re
import time
import random
import signal
import glob
import textwrap
import getopt
import string
from subprocess import call

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(  "%s/../utils" % curdir)
from struct_methods import LogIt

g_not_found_retval = ([-1, "?", "----------"], [-1, "?", "----------"], -1)

counted_aas = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

#for only poisite windowed sequences
#seq_res_num starts from 1
#returns tuple of arrays: ( [pdb_res_pos,pdbres,pdbseq], [seqresnum,seqres,seqseq], align_score)
#pdb_res_pos is Nth residue in the pdb file (and chain) and the first N-terminal residue of the poisite
def FindResidueInPoiSiteAlignment( alignment_file, cutseq_length, score_cutoff=-1, log_file="", aPositions=0):

    global g_not_found_retval, counted_aas

    PDB = 0
    SEQ = 1
    SCORE = 2

    #Check alignment score if threshold specified
    #;ALIGNMENT_SCORE:100
    alignment_score = -1

    #Read alignment Score
    if score_cutoff >= 0:

        found = False
        re_score = re.compile('ALIGNMENT_SCORE:\s*(\d+)', re.IGNORECASE)

        try:
            with open( alignment_file, 'r') as f:

                while not found:

                    line = f.readline()
                    if len( line) == 0:
                        LogIt( "ERROR: End of file reached in file '%s'. (FindResidueInPoiSiteAlignment)" % alignment_file, log_file, 2 )
                        break #End of file reached

                    line = line.strip()

                    if len( line) == 0:
                        continue
                    elif line[ 0] == ";" or line[ 0] == "#":
                        r = re_score.search( line)
                        if r:
                            alignment_score = int( r.group(1))
                            found = True
                    else:
                        break
        except ValueError:
            LogIt( "WARNING: Invalid alignment score found in file '%s'." % alignment_file, log_file, 2  )
        except:
            pass

        if not found:
            LogIt( "WARNING: File '%s' did not have an alignment score.\n" % alignment_file, log_file, 2  )
            return g_not_found_retval #RETURN
        elif alignment_score < score_cutoff:
            return (g_not_found_retval[ 0], g_not_found_retval[ 1], alignment_score) #RETURN


    try:
        align_records = list( SeqIO.parse( alignment_file, "fasta"))
    except:
        LogIt(  "Could not find or open file '%s'.\n" % alignment_file, log_file, 2  )
        return g_not_found_retval


    if len( align_records) < 2:
        LogIt(  "ERROR: Alignment file '%s' contains less than 2 alignments.\n" % alignment_file, log_file, 2  )
        return g_not_found_retval
    if len( align_records) > 2:
        LogIt(  "WARNING: Alignment file '%s' contains more than 2 alignments.\n" % alignment_file, log_file, 2  )
        LogIt(  "         Only first two are used.\n", log_file, 2  )


    #PDB file is first
    re_cutnum = re.compile('POI:\s*(\d+)-(\d+)', re.IGNORECASE)
    re_cutseq = re.compile('\s*POISEQ:\s*([A-Z-]+)/([A-Z-]+)', re.IGNORECASE)


    seq_res_num =  -1
    res = re_cutnum.search( align_records[ SEQ].description)
    if not res:
        LogIt(  "ERROR: Alignment file '%s' has no POI information in fasta header.\n" % alignment_file, log_file, 2  )
        return g_not_found_retval #RETURN

    try:
        cut_1 = int( res.group( 1))
        cut_2 = int( res.group( 2))
        seq_res_num = cut_1
    except ValueError:
        LogIt(  "ERROR: Could not convert '%s' and '%s' to integer in file '%s'.\n" % (res.group( 1), res.group( 2), alignment_file), log_file, 2  )
        return g_not_found_retval #RETURN

    res = re_cutseq.search( align_records[ SEQ].description)
    if not res:
        LogIt(  "WARNING: Alignment file '%s' has no POISEQ information.\n" % alignment_file, log_file, 2 )
        return g_not_found_retval #RETURN

    before, after = GetSurroundings( cutseq_length, cut_1, align_records[ SEQ].seq)
    seq_res_aa = res.group( 1)[ -1]
    seq_res_seq = "%s*%s" % (res.group( 1), res.group( 2))


    pdb_res_pos = 1
    pdb_res_aa = "?"
    pdb_res_seq = "----------"

    pdb_seq = align_records[ PDB].seq
    seq_seq = align_records[ SEQ].seq

    found = False

    cutseq = ""
    n_of_residues = 0

    #Check sequence
    for s in range( len( seq_seq)):

        aa = seq_seq[ s]

        if aa in counted_aas:
            n_of_residues += 1
            cutseq += aa

    if n_of_residues < cutseq_length/2:
        LogIt(  "ALIGN WARNING: short cutsite sequence: '%s'.\n" % cutseq, log_file, 2  )

    if n_of_residues > cutseq_length:
        LogIt(  "ALIGN ERROR: cut site sequence too long (%i).\n" % n_of_residues, log_file, 2  )
        return g_not_found_retval #RETURN

    seq_to_go = ((cutseq_length+(n_of_residues%2)) / 2)
    #print "seq_to_go: %i " % seq_to_go

    #if positions if 0 there needs to be a C-terminal residue after POI
    range_limit = 0 if aPositions > 0 else 1
    range_end = len( pdb_seq)

    for s in range( len( seq_seq)-range_limit):

        aa = seq_seq[ s]

        if aa in counted_aas:
            seq_to_go -= 1

        if seq_to_go <= 0:

            if aPositions == 0: #Both sides of POI
                if pdb_seq[ s] == "-" or pdb_seq[ s+1] == "-":
                    LogIt( "ALIGN INFO: PDB structure has a gap (unstructured region) at the POI-site, in file '%s'.\n" % os.path.basename( alignment_file), log_file, 1 )
                    return g_not_found_retval
            else:
                if pdb_seq[ s] == "-":
                    LogIt( "ALIGN INFO: PDB structure has a gap (unstructured region) at the POI-site, in file '%s'.\n" % os.path.basename( alignment_file), log_file, 1 )
                    return g_not_found_retval


            n_side = str( pdb_seq[ max( 0, s-4):min(range_end, s+1)].upper())
            while len( n_side) < 5: n_side = "-" + n_side

            c_side = str( pdb_seq[ min(range_end, s+1):min(range_end, s+6)].upper())
            while len( c_side) < 5: c_side += "-"

            pdb_res_seq = n_side + "*" + c_side

            #pdb_res_seq = pdb_seq[ max(0,s-4):s].upper() + pdb_seq[ s].upper() + "*" + pdb_seq[ s+1].upper() + pdb_seq[s+2:s+6].upper()

            pdb_res_aa = pdb_seq[ s].upper()
            found = True
            break

        if pdb_seq[ s] in counted_aas:
            #print pdb_seq[ s],
            pdb_res_pos += 1

    if not found:
        return g_not_found_retval


    return ([pdb_res_pos, pdb_res_aa, str(pdb_res_seq)], [seq_res_num, seq_res_aa, seq_res_seq], alignment_score)


def GetSurroundings( aWindow, aPOI_pos, aSeq):

    before = ""
    after = ""

    aWindow = aWindow - (aWindow % 2) #Make it even

    cur_pos = aPOI_pos
    while len( before) < aWindow/2 and cur_pos <= 0:
        aa = aSeq[ cur_pos]
        if aa != "-": before += aa
        cur_pos -= 1

    before = before[::-1] #reverse
    before = ("-"*((aWindow/2)-len( before)))+before #pad with gaps if necessary

    cur_pos = aPOI_pos+1
    while len( after) < aWindow/2 and cur_pos < len( aSeq):
        aa = aSeq[ cur_pos]
        if aa != "-": after += aa
        cur_pos += 1

    after = after+("-"*((aWindow/2)-len( after))) #pad with gaps if necessary

    return (before, after)


#For full sequences
#seq_res_num starts from 1
#so does returned PDB res num
def FindResidueInAlignment( alignment_file, seq_res_num):

    global g_not_found_retval, counted_aas

    PDB = 0
    SEQ = 1

    try:
        align_records = list( SeqIO.parse( alignment_file, "fasta"))
    except:
        sys.stderr.write( "ERROR: Could not find or open file '%s'.\n" % alignment_file )
        return g_not_found_retval


    if len( align_records) < 2:
        sys.stderr.write( "ERROR: Alignment file '%s' contains less than 2 alignments.\n" % alignment_file )
        return g_not_found_retval
    if len( align_records) > 2:
        sys.stderr.write( "WARNING: Alignment file '%s' contains more than 2 alignments.\n" % alignment_file )
        sys.stderr.write( "         Only first two are used.\n")


    #PDB file is first

    pdb_res_num = 1
    pdb_res_aa = "?"
    pdb_res_seq = "----------"

    seq_res_aa = "?"
    seq_res_seq = "?"


    seq_to_go = seq_res_num

    pdb_seq = align_records[ PDB].seq
    seq_seq = align_records[ SEQ].seq



    found = False

    #Loop until to penultimate residue
    for s in range( len( seq_seq)-1):

        seq_res_aa = seq_seq[ s]

        if seq_res_aa in counted_aas:
            seq_to_go -= 1

        if seq_to_go <= 0:
            pdb_res_aa = pdb_seq[ s]
            #Residues between both sides of cutsite have to be matched in the alignment
            if pdb_seq[ s] == "-" or pdb_seq[ s+1] == "-":
                pdb_res_num = -1
            pdb_res_seq = pdb_seq[max(0,s-4):s].upper() + pdb_seq[ s].upper() + "*" + pdb_seq[ s+1].upper() + pdb_seq[s+2:s+6].upper()
            seq_res_seq = seq_seq[max(0,s-4):s].upper() + seq_seq[ s].upper() + "*" + seq_seq[ s+1].upper() + seq_seq[s+2:s+6].upper()
            found = True
            break
        elif pdb_seq[ s] in counted_aas:
            pdb_res_num += 1

    if not found:
        return g_not_found_retval

    return ([pdb_res_num, pdb_res_aa, str(pdb_res_seq)], [seq_res_num, seq_res_aa, str(seq_res_seq)])



def PrintAlignHelp():

    print """
    This script find the corresponding residue of a PDB file based on an sequence alignment
    Residue numbering starts from 1

    Usage: script.py alignment_file residue_number

    -h  --help             Print this message.
    -v  --version          1.0
    """

def main():

    align_file = ""
    aa_num = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hv', ['help', 'version'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)


    #Input & Output files
    for arg in args:
        if len( align_file) == 0:
            align_file = arg
        elif len( aa_num) == 0:
            aa_num = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)


    #Flags
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            PrintAlignHelp()
            sys.exit( 0)
        elif opt in ('-v', '--version'):
            PrintAlignHelp()
            sys.exit( 0)

    try:
        aa_num = int( aa_num)
    except ValueError:
        sys.stderr.write( "Second parameter needs to be a number (first aa is 1).\n")
        sys.exit( -2)

    res, orig = FindResidueInAlignment( align_file, aa_num)

    if res[ 0] < 0:
        print "No matching residue was found."
    else:
        print "Seq %i aa '%s', in \"%s\"." % tuple( orig)
        print "Mathcing PDB file residue is number %i, '%s' in seq \"%s\"." % tuple( res)


if __name__ == "__main__":
  main()

