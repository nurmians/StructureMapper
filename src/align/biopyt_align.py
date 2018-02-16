#!/usr/bin/python
# coding=UTF-8
# -*- coding: UTF-8 -*-

import sys
import os
import re
from subprocess import call


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

curdir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(  "%s/../utils" % curdir)
from struct_methods import *
from fileops import *
from obj_cache import ObjCache
import online_pdb

sys.path.append(  "%s/../asa" % curdir)
import pdbatoms

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

g_pairwise_plock = False

#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/
#g_matrix = matlist.blosum62
g_matrix = matlist.pam250
g_gap_open = -10
g_gap_extend = -2.0

#Cache recently used PDB sequences so that the same PDB does not get read over and over
SEQ_CACHE = ObjCache( 20)
LEN_CACHE = ObjCache( 20)

def SetSeqCacheSize( aSize):
    global SEQ_CACHE, LEN_CACHE
    SEQ_CACHE.SetSize( aSize)
    LEN_CACHE.SetSize( aSize)

def SetPairwisePrintLock( aLock):
	global g_pairwise_plock
	g_pairwise_plock = aLock

#Uses only sequence around cutsite
#Return 1 = OK, 0 = skipped, -1 = Error
def AlignCutSeqWindowWithPairwise( thread_id, cutsite, cutseqlen, record_identifier, fasta_record, PDB_code, PDB_filename, PDB_chain, output_folder, skip_existing=True, log_file=None):

    #Use only sequence around cutsite
    cutsite_start = max(0, (cutsite-(cutseqlen/2)))
    cutsite_end = cutsite_start+cutseqlen

    #Create new record with just the part of the sequence around the cutsite
    new_record = SeqRecord( fasta_record.seq[ cutsite_start:cutsite_end], id=fasta_record.id, name=fasta_record.name, description=fasta_record.description)

    #print "RECORD ID: '%s'" % new_record.id
    #print "RECORD DESC: '%s'" % new_record.description

    return AlignCutSeqWithPairwise( thread_id, record_identifier, new_record, PDB_code, PDB_filename, PDB_chain, output_folder, skip_existing, log_file=log_file)


#Return 1 = OK, 0 = skipped, -1 = Error
#Use full sequence for alignment
def AlignCutSeqWithPairwise( thread_id, record_identifier, fasta_record, PDB_code, PDB_filename, PDB_chain, output_folder, skip_existing=False, log_file=None):

    global g_pairwise_plock, g_matrix,g_gap_open, g_gap_extend, SEQ_CACHE, LEN_CACHE

    PDB_code = PDB_code.upper()
    PDB_chain = PDB_chain.strip()

    #In blast results double upper chain means lower case: FF -> f
    if len( PDB_chain) == 2:
        if PDB_chain[ 0] == PDB_chain[ 1]: PDB_chain = PDB_chain[ 0].lower()
        else: LogIt( "WARNING: Dubious chain identifier '%s' for record '%s' and structure %s." % (PDB_chain, record_identifier, PDB_code), log_file, 2, lock=g_pairwise_plock )

    #Replace chars that are unsuitable for filenames with underscores
    record_identifier = record_identifier[:100] #max length 100
    record_filename = ValidFilename( record_identifier)
    output_folder = InsertFolderSeparator( output_folder)
    record_path = output_folder + record_filename

    #Write input file for Clustal
    output_file = record_path + ".align"
    #output_final = record_path + ".align"
    #clustal_temp_scorefile  = output_folder + "clustal_score_%i.tmp" % thread_id
    error_log_file  = output_folder + "pairwise_error_log_%i.txt" % thread_id

    #print "Reading sequence from PDB file: '%s' " % PDB_filename
    #OLD
    #pdb_length, pdb_sequence = GetSeqLengthFromPDBFile( PDB_code, PDB_filename, PDB_chain, include_gaps=True, quiet=True)

    #NEW
    cache_key = PDB_code + PDB_chain
    found_in_cache = LEN_CACHE.IsInCache( cache_key)
    pdb_sequence = ""
    pdb_length = -1

    #Announce job str
    align_log_str = "[%s]%i Aligning: %s|%s <-> %s" % (("+" if found_in_cache else "-"), thread_id, PDB_code, PDB_chain, record_identifier)

    if skip_existing and os.path.isfile( output_file):
        #If Verbose?
        #SyncPrint( align_log_str + " (already exists, skipped)\n  File: '%s' " % clustal_output_filename, lock=g_pairwise_plock)
        return 0
    else:
        #LogIt( "{:<80}".format(align_log_str) + " [%i]" % thread_id, log_file, 1, lock=g_pairwise_plock )
        #LogIt( align_log_str.ljust( 80) + " [%i]" % thread_id, log_file, 1, lock=g_pairwise_plock )
        LogIt( align_log_str, log_file, 1, lock=g_pairwise_plock )

    if found_in_cache:
        #Use Cached
        pdb_length = LEN_CACHE.Retrieve( cache_key)
        pdb_sequence = SEQ_CACHE.Retrieve( cache_key)
        #print "INFO: Seq for %s found in thread alignment cache." % (PDB_code + "|" + PDB_chain)
    else:

        #Make sure pdb file exists
        if not os.path.isfile( PDB_filename):
            pdb_folder = GetFolder( PDB_filename)
            fetch_err = online_pdb.FetchPDBFile( PDB_code, output_folder=pdb_folder, obsoletes=online_pdb.FetchObsoletes( write_to_file=True, folder=pdb_folder), aAsync=False, verbose=False)
            if fetch_err == 0: LogIt( "ALIGN INFO: PDB file '%s' fetched for alignment. [OK]" % PDB_filename, log_file, 1, lock=g_pairwise_plock )
            else: LogIt( "ALIGN ERROR: PDB file '%s' could not be fetched for alignment." % PDB_filename, log_file, 2, lock=g_pairwise_plock )

        pdb = None
        try:
            with open( PDB_filename, "r" ) as pdb_fh:
                pdb = pdbatoms.ReadPDB( pdb_fh, PDB_chain, firstmodelonly=True, aQuiet=True )
        except Exception as ex:
            LogIt( "ERROR: Could not read file '%s'. (%s)" % (PDB_filename, str( ex)), log_file, 2, lock=g_pairwise_plock )
            return -3 #Fatal error

        if pdb["stats"]["n_ptm_aas"] > 0:
            LogIt( "INFO: Converting %i PTMs in %s." % (pdb["stats"]["n_ptm_aas"], PDB_filename), log_file, 1, lock=g_pairwise_plock )
            pdbatoms.ConvertPTMAAstoRegularAtoms( pdb, aRemovePTMs=False, aQuiet=True) #Convert without removal of extra atoms is slighly faster
        pdb_sequence = pdbatoms.FastaSequence( pdb["residues"], aAtomTypes=["A"], aInsertGaps=True, aQuiet=True) #Only amino acids are written as seq
        pdb_length = len( pdb["residues"])

        LEN_CACHE.Add( cache_key, pdb_length)
        SEQ_CACHE.Add( cache_key, pdb_sequence)


    if pdb_length < 1:
        LogIt( "ERROR: No PDB sequence found for %s <-> %s File: '%s'" % (record_identifier, PDB_code + "|" + PDB_chain, PDB_filename), log_file, 2, lock=g_pairwise_plock)
        return -98 #RETURN
    elif pdb_length < 10:
        LogIt( "WARNING: Short PDB alignment sequence (len:%i) for %s <-> %s File: '%s'" % (pdb_length, record_identifier, PDB_code + "|" + PDB_chain, PDB_filename), log_file, 2, lock=g_pairwise_plock)


    #ALIGN!
    alns = pairwise2.align.localds( str( pdb_sequence).replace( "-", ""), str( fasta_record.seq), g_matrix, g_gap_open, g_gap_extend, penalize_end_gaps=False, one_alignment_only=True)
    #alns = pairwise2.align.localms( pdb_sequence, fasta_record.seq, 1, -3, -7, -2, one_alignment_only=1, score_only=1)

    if len( alns) == 0:
        sys.stderr.write( "WARNING: Alignment failed for file '%s' and seq '%s'.\n" % (PDB_filename, fasta_record.id))

    try:
        top_aln = alns[ 0]
    except IndexError:
        sys.stderr.write( "ERROR: Alignmnet failed for file '%s' and seq '%s'.\n" % (PDB_filename, fasta_record.id))
        sys.stderr.write( "ERROR: '%s'\n<=>\n'%s'\n" % ( str( pdb_sequence).replace( "-", ""), str(fasta_record.seq)))
        raise


    aln_pdb, aln_seq, score, begin, end = top_aln
    #print aln_pdb+'\n'+aln_seq

    #Own scoring method
    score = ScorePairwise2Alignment( aln_pdb, aln_seq)

    #If verbose?
    #print "PDB:", aln_pdb
    #print "SEQ:", aln_seq

    #print alns

    out_records = []
    out_records.append( SeqRecord( Seq( aln_pdb), id="PDB_FILE|" + PDB_code + "|chain:" + PDB_chain, description=";" + GetPDBFileTitle( PDB_filename).capitalize()))
    out_records.append( SeqRecord( Seq( aln_seq), id=fasta_record.id, name=fasta_record.name, description=fasta_record.description))


    SeqIO.write(out_records, output_file, "fasta")

    #Write alignment score to file
    try:
        PrependToFile( output_file, ";ALIGNMENT_SCORE:%i\n" % int( round( score)))
    except IOError as ex:
        LogIt( "ERROR: Unable to write score to alignment file: '%s' Error: %s\n" % (clustal_output_filename, ex.strerror), log_file, 2, lock=g_pairwise_plock)
        LogIt( "       An exception of type {0} occured. Arguments: {1!r}".format( type( ex).__name__, ex.args), 2, lock=g_pairwise_plock )
        return -12


    #if score_it:
    #    score_file = output_folder + record_filename + ".alignscore"
    #    ScoreAlignment( record_filename, PDB_code + "|" + PDB_chain, pdb_length, clustal_output_filename, score_file)

    #If verbose?
    #print "SCORE: ", score

    return 1 #OK


def ScorePairwise2Alignment( aPdbSeq, aCutSeq):

	cutseqlen = 0
	gaps = 0
	identities = 0
	start = -1
	end = -1

	for i in range( len( aCutSeq)):
		aa = aCutSeq[ i]
		if aa != "-":
			if start < 0: start = i
			end = i
			cutseqlen += 1
			if aPdbSeq[ i] == aa : identities += 1

	#print "SEQM:", aCutSeq[ start:end]
	if start < 0 or end < 0 or cutseqlen < 0: return 0
	middle_gaps = aCutSeq[ start:end].count( "-") + aPdbSeq[ start:end].count( "-")

	percentage = (identities / float( cutseqlen)) * 100.0 #identities
	percentage -= middle_gaps *  (100.0/cutseqlen)        #penalties
	return int( max( 0, round( percentage)))
