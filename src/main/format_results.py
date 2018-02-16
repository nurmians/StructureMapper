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

#print "%s/../utils" % os.path.dirname(__file__)
#sys.exit( 0)
curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(  "%s/../utils" % curdir)

from struct_methods import *
from fileops import ChangeToFolder



#Dictionary (Acc without pos) of lists (blast result lines) of dictionaries (blast row vals)
#Blast results are not dependent on POI pos identifiers, the index uses acc numbers without "_Pos"
#keys = "rank", query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start",
#        "s_end", "evalue", "bitscore","pdb_id","chain","struct_id","rank"
def CreateBlastResultIndex( aBlastFile, pdb_folder, aObsoletes=[]):


    if not aObsoletes or len( aObsoletes) == 0:
        aObsoletes = FetchObsoletes( True, pdb_folder)

    #Return all blast result rows for each sequence
    blast_index = {}
    #seq_results = []
    #first = True
    error_count = 0
    #n_blast_results = 0


    # Fields: query id,         subject id,             % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    #sp|P31946-2|1433B_HUMAN    gi|67464627|pdb|2BQ0|A  99.58       238 1   0   1   238 3   240 2e-176   489
    blast_col_keys = [ "query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
    blast_col_types= [ "s",        "s",          "f",        "i",                "i",          "i",         "i",       "i",     "i",       "i",     "s",      "f"]

    try:
        blast_record_reader = FastaIter( aBlastFile, False, None) #Generator works on blast result file even though not a fasta file
    except:
        sys.stderr.write("ERROR: Could not open blast result file '%s'.\n" % blast_file)
        return {} #RETURN

    cur_rec = -1

    for blast_header, blast_entries in blast_record_reader:

        if error_count > 20:
            sys.stderr.write("20+ errors encountered. Quitting.\n")
            return {} #RETURN

        cur_rec += 1
        #seq_acc = GetRecordACC( blast_header, aWithPos=True)
        acc = GetRecordACC( blast_header, aWithPos=False)

        if acc in blast_index: continue #Already in the index

        blast_index[ acc] = [] # new list
        unique_structs_for_seq = {}
        rank = 0

        for blast_entry in blast_entries.splitlines():

            rank += 1 #rank advances even if duplicate or invalid
            blast_cols = blast_entry.split( "\t")

            if len( blast_cols) != 12:
                print "ERROR: Bad Blast result for entry '%s'. (%i/12 tab separated cols)" % (blast_header, len( blast_cols))
                error_count += 1
                continue
            else:

                blast_dict = {}

                struct_cols = blast_cols[ 1].split( '|')
                blast_dict["pdb_id"]  = struct_cols[ 3].strip()
                blast_dict["chain"]  = struct_cols[ 4]
                blast_dict["rank"] = rank
                blast_dict["struct_id"]  = "%s|%s" % ( blast_dict["pdb_id"], blast_dict["chain"])

                #Remove obsolete entries
                if blast_dict["pdb_id"] in aObsoletes: continue
                #Add PDBID + CHAIN pairs only once per sequence
                if blast_dict["struct_id"] in unique_structs_for_seq: continue


                #blast_dict["entry_id"] = "%s_%i" % (seq_acc, rank)
                blast_dict["pdb_filepath"] = PDBFileFromAccCode( blast_dict["pdb_id"], pdb_folder)
                #blast_dict["align_file"] = align_file_format % blast_dict["entry_id"]

                for c in range( len(blast_cols)):
                    blast_dict[ blast_col_keys[ c]] = TypeCast( blast_cols[ c], blast_col_types[ c], 1)

                #seq_results.append( blast_dict)
                unique_structs_for_seq[ blast_dict["struct_id"]] = True

                blast_index[ acc].append( blast_dict)


    return blast_index

# Fields: query id,         subject id,             % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#sp|P31946-2|1433B_HUMAN    gi|67464627|pdb|2BQ0|A  99.58       238 1   0   1   238 3   240 2e-176   489
BLAST_COL_KEYS = [ "query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
BLAST_COL_TYPES= [ "s",        "s",          "f",        "i",                "i",          "i",         "i",       "i",     "i",       "i",     "s",      "f"]


#Input Blast results separated by newlines
#Return list (blast result lines) of dictionaries (blast row vals)
def ParseBlastEntries( aBlastEntries, aPdbFolder="", aObsoletes=[]):

    global BLAST_COL_KEYS, BLAST_COL_TYPES

    rank = 0
    unique_structs_for_seq = {}
    entry_list = [] #retvals

    for blast_entry in aBlastEntries.splitlines():

        rank += 1 #rank advances even if duplicate or invalid
        blast_cols = blast_entry.split( "\t")

        if len( blast_cols) != 12:
            print "ERROR: Bad Blast result for entry '%s'. (%i/12 tab separated cols)" % (blast_header, len( blast_cols))
            error_count += 1
            continue
        else:

            blast_dict = {}
            chain = ""

            #Two possible formats in BLAST results:
            if blast_cols[ 1].find( "|") >= 0:
                #gn|P23919|DTYMK    gi|359545626|pdb|2XX3|A 100.00  212 0   0   1   212 21  232 7e-156   436
                struct_cols = blast_cols[ 1].split( '|')
                blast_dict["pdb_id"]  = struct_cols[ 3].strip()
                chain = struct_cols[ 4]
            else:
                #gn|P30613_Pos2|PKLR     4IP7_A  99.815  541     1       0       34      574     3       543     0.0     1096
                struct_cols = blast_cols[ 1].split( '_')
                blast_dict["pdb_id"] = struct_cols[ 0].strip() #3IKM
                chain = struct_cols[ 1]

            ##Blast results can have "GG" when it should be "g"
            if len( chain) == 2 and chain[ 0] == chain[ 1]: chain = chain[ 0].lower()
            blast_dict["chain"]  = chain

            blast_dict["rank"] = rank
            blast_dict["struct_id"]  = "%s|%s" % ( blast_dict["pdb_id"], blast_dict["chain"])

            #Remove obsolete entries
            if blast_dict["pdb_id"] in aObsoletes: continue
            #Add PDBID + CHAIN pairs only once per sequence
            if blast_dict["struct_id"] in unique_structs_for_seq: continue


            #blast_dict["entry_id"] = "%s_%i" % (seq_acc, rank)
            if len( aPdbFolder): blast_dict["pdb_filepath"] = ChangeToFolder( PDBFileFromAccCode( blast_dict["pdb_id"]), aPdbFolder)
            #blast_dict["align_file"] = align_file_format % blast_dict["entry_id"]

            for c in range( len(blast_cols)):
                blast_dict[ BLAST_COL_KEYS[ c]] = TypeCast( blast_cols[ c], BLAST_COL_TYPES[ c], 1)

            #seq_results.append( blast_dict)
            unique_structs_for_seq[ blast_dict["struct_id"]] = True

            entry_list.append( blast_dict)

    return entry_list
