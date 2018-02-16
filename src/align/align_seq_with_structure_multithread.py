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
import traceback
import re
import time
import random
import signal
import glob
import textwrap
import getopt
import string
import gc
import itertools
from subprocess import call
from datetime import datetime, timedelta
#from dateutil.relativedelta import relativedelta

from multiprocessing import Process, Pool, Array, Lock, Value, cpu_count, Manager

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList

#from needle_align import *
from biopyt_align import * #AlignCutSeqWindowWithPairwise

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(  "%s/../utils" % curdir)
from struct_methods import *
from fileops import *


#valid_filename_chars = "-_.()[]{}= %s%s" % (string.ascii_letters, string.digits)


#Shared Global variables
align_printlock = False
n_blast_entries_processed = Value('i', 0)
n_alignments_done = Value('i', 0)
n_alignments_skipped = Value('i', 0)
n_alignments_error = Value('i', 0)
fatal_errors = Value('i', 0)
obsoletes_dict = {}


#Init func when running only Alignments with Pool
def PoolInitializer( *args):
    global n_alignments_done, n_alignments_skipped, n_alignments_error, align_printlock, obsoletes_dict, n_blast_entries_processed
    n_alignments_done, n_alignments_skipped, n_alignments_error, align_printlock, obsoletes_dict, n_blast_entries_processed = args
    #SetSyncPrintBuffered( False)
    signal.signal( signal.SIGINT, signal.SIG_IGN)

#FUNCTIONS
def ReportProcessStart( file):
    sys.stdout.write("Processing file '%s'.\n" % (file))


def ReportDone( thread_tuple):

    global align_printlock, fatal_errors

    SyncPrint( "Thread %i: Done. (%i)\n" % thread_tuple, aLock=align_printlock)

    if thread_tuple[ 1] == -99: #User interrupted process
        fatal_errors.value = -99
        SyncPrint( "Shutting down processes...\n", aLock=align_printlock)
    elif thread_tuple[ 1] < 0:
        fatal_errors.value = 1
        SyncPrint( "Shutting down processes...\n", aLock=align_printlock)

    #global printlock, number
    #SyncPrint( "File processed.", printlock )
    #SyncPrint( "Files done: %i %s" % ( str(file)), printlock)



def GetFastaFilesInFolder( folder):

    folder = InsertFolderSeparator( folder)
    glob_path1 = folder + "*.fasta"

    file_list = []

    for globfile in glob.glob( glob_path1):
        file_list.append( globfile)


    return file_list


def ScoreAlignment( seq_id, PDB_ID, PDB_LEN, filename, scorefilename):

    global align_printlock
    #Read ClustalW output file
    #outfile = open( "clustal_aligned.fasta", 'r')

    #alignment = SeqIO.parse( filename, "fastq")
    alignment_records = list( SeqIO.parse( filename, "fasta"))

    if len( alignment_records) != 2:
        SyncPrint( "ERROR: Bad alignment file '%s': contains %i records" % (filename, len( alignment_records)), align_printlock )
        return #RETURN

    seq_pdb = alignment_records[ 0].seq
    seq_proteome = alignment_records[ 1].seq

    #print "Seq PDB: % s " % seq_pdb
    #print "Seq Proteome: % s" % seq_proteome

    align_first = -1
    align_last = -1
    unaligned = 0
    identity = 0
    gaps = 0
    gap_opens = 0
    maybe_gaps = 0


    scorefile = open( scorefilename, "w")



    #If PDB sequence is longer, set unaligned to excess amount
    if len( seq_pdb) > len( seq_proteome):
        unaligned = len( seq_pdb) - len( seq_proteome)

    for i in range( 0, min( len( seq_pdb), len( seq_proteome))):

        #Gap in Uniprot seq
        if seq_pdb[ i] != '-' and seq_proteome[ i] == '-':
            unaligned += 1
        #Same residue in both
        elif seq_pdb[ i] != '-' and seq_proteome[ i] == seq_pdb[ i]:
            identity += 1
        #Gap in PDB sequence
        elif align_first >= 0 and seq_pdb[ i] == '-' and seq_proteome[ i] != '-':
            maybe_gaps += 1

        #Aligned residues
        if seq_pdb[ i] != '-' and seq_proteome[ i] != '-':
            align_last = i
            gaps += maybe_gaps

            if maybe_gaps > 0:
                gap_opens += 1

            maybe_gaps = 0

            if align_first < 0:
                align_first = i


    alignment_length = 0

    if align_first >= 0:
       alignment_length = align_last - align_first + 1

    #hyphen_count = seq_pdb.count('-')

    #Write score
    # SEQ_ID SEQ_LEN, PDB_ID, PDB_LEN, Alignment length,
    scorefile.write( "%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % ( seq_id, len( seq_proteome), PDB_ID, PDB_LEN, alignment_length, align_first, align_last, identity, gaps, unaligned, gap_opens))

    #Debug
    #alignment_length_file.close()
    #sys.exit("debug stop")

    scorefile.write("\n")
    scorefile.close()

def CreateAlignFilename( aProjectPrefix, aRecordId, aPDB_code, aPDB_chain):

    filename = "%s_[%s=%s(%s)].align" % (aProjectPrefix[:100], aRecordId.upper(), aPDB_code.upper(), aPDB_chain.upper())
    return ValidFilename( filename)



def RecordSelector( aCurIndex):
    #global align_printlock
    if (aCurIndex / RecordSelector.chunksize) % RecordSelector.tc != RecordSelector.tid: return False
    #SyncErrPrint( "YES for %i! tid: %i tc: %i" % (aCurIndex, RecordSelector.tid, RecordSelector.tc), align_printlock)
    return True


def InitRecordSelector( aThreadCount, aThreadId, aChunkSize=10):
    RecordSelector.tc = aThreadCount
    RecordSelector.tid = aThreadId
    RecordSelector.chunksize = aChunkSize

def FillAlignmentQueue( blast_record_reader, blast_entries, seq_index, fasta_records, aQueueSize = 10 ):

    global align_printlock

    appended = 0

    for blast_header, blast_record in blast_record_reader:

        record_id = GetRecordACC( blast_header)

        try:
            fasta_records.append( seq_index[ record_id])
            blast_entries.append( (blast_header, blast_record))
        except Exception as ex:

            SyncErrPrint( "THREAD ERROR: No matching sequence entry found for blast entry '%s':\n" % (record_id), align_printlock)
            SyncErrPrint( "SeqIndex size: %i" % len( seq_index), align_printlock)

            message = "An exception of type {0} occured. Arguments: {1!r}".format( type(ex).__name__, ex.args)
            SyncErrPrint( message, align_printlock)

            return -2 #RETURN ERROR

        appended += 1

        if len( blast_entries) >= aQueueSize: break

    return appended

def ProcessBlastResultsAsync( aBlastResultFile, aSequenceFile, aOutputFolder, aProcessID, aPDB_Folder, aObsoletes, aScoreIt=False, aSkipExisting=True, aLock=False):

    args_str =( "%s\t%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i" % (aBlastResultFile, aSequenceFile, aOutputFolder, 1, aProcessID, aBlastResultFile + ".log", aPDB_Folder, aScoreIt, aSkipExisting))

    #p = Process( target=ProcessBlastMultiThreadInit, args=(args_str, aObsoletes, aLock, g_SyncPrintBuffer, g_SyncPrintBuffered, g_SyncPrintBuffer_overflow))
    p = Process( target=ProcessBlastMultiThreadInit, args=(args_str, aObsoletes, aLock, GetSyncPrintBufferForThread()))
    p.start()
    return p


#For parallel processing when blasting
#threaded
#def ProcessBlastMultiThreadInit( arg_str, aObsoletes, aLock, aSyncPrintBuffer, aSyncPrintBuffered, aSyncPrintOverFlow):
def ProcessBlastMultiThreadInit( arg_str, aObsoletes, aLock, aSyncPrintBuffer):
    global align_printlock, obsoletes_dict
    align_printlock = aLock
    obsoletes_dict = aObsoletes
    thread_id = arg_str.split( "\t")[ 4]
    blast_file = arg_str.split( "\t")[ 0]

    #SetSyncPrintBufferInThread( aSyncPrintBuffered, aSyncPrintBuffer, aSyncPrintOverFlow)
    SetSyncPrintBufferInThread( aSyncPrintBuffer)
    SetAlignPrintLock( aLock)

    #DEBUG
    #print "\n\n\nIsSyncPrintBuffered in thread %i: %i" % (int(thread_id), IsSyncPrintBuffered()), "LOCK:", aLock
    #print "\n\n\n"

    SyncPrint( "Align Progress: Aligning blast results '%s'. [%s]" % (blast_file, thread_id), aLock)

    sig = signal.signal(signal.SIGINT, signal.SIG_IGN)
    ProcessBlastMultiThread( arg_str)
    signal.signal(signal.SIGINT, sig)

    SyncPrint( "Align Progress: File '%s' Done. [%s]" % (blast_file, thread_id), aLock)

    SetAlignPrintLock( False)


#Threaded
#def ProcessBlastMultiThread( arg_str, aLock ):
def ProcessBlastMultiThread( arg_str ):
    global align_printlock
    retval = -1
    thread_id = -1

    try:
        thread_id = int( arg_str.split( "\t")[ 4])
        retval = DoProcessBlastMultiThread( arg_str)

    except Exception as ex:

        err_str = "ERROR: Thread %i:\n" % thread_id
        message = "       An exception of type {0} occured. Arguments: {1!r}".format( type(ex).__name__, ex.args)
        sys.stderr.write( "\n" + err_str + message + "\n")
        SyncPrint( "\n" + err_str + message + "\n", align_printlock)
        traceback.print_exc(file=sys.stderr)

    return (thread_id, retval)

#Threaded
def DoProcessBlastMultiThread( arg_str ):

    global n_alignments_done, n_alignments_skipped, n_alignments_error, align_printlock, obsoletes_dict, n_blast_entries_processed#, seq_index

    #read arguments
    cols = str(arg_str).split( "\t")

    blast_file = cols[ 0]
    sequence_file = cols[ 1]
    output_folder = cols[ 2]
    thread_count = int( cols[ 3])
    thread_id = int( cols[ 4])
    log_file = cols[ 5]
    PDB_folder = cols[ 6]
    score_it = True if int( cols[ 7]) else False
    skip_existing = True if int( cols[ 8]) else False
    #sys.stderr.write( "Align Progress: Debug Thread B '%i'.\n" % thread_id)
    #raise ValueError('A very specific bad thing happened')

    #SyncPrint( "Obsoletes %i: %i (%i)" % (thread_id, len(obsoletes_dict), obsoletes_dict['125D'] ))
    #SyncPrint( "Score it: %i " % (score_it), align_printlock)
    #SyncPrint( "skip_existing: %i " % (skip_existing), align_printlock)

    path, filename = os.path.split( sequence_file)
    name, suffix = SplitFileName( filename)

    seq_index = SeqIO.index( sequence_file, "fasta", key_function=GetRecordACC)
    n_records = len( seq_index)
    n_total_blast_results = 0
    n_errors = 0


    output_folder = InsertFolderSeparator( output_folder)
    log_filename = output_folder + "align_%i_log.txt" % thread_id
    align_scores_filename = output_folder + "align_%i_scores.txt" % thread_id
    re_cutnum = re.compile('POI:\s*(\d+)-(\d+)', re.IGNORECASE)

    no_struct_found = 0
    align_ok = 0
    align_skipped = 0
    align_error = 0

    #sys.stderr.write( "Align Progress: Debug Thread '%i'." % thread_id)

    #SyncPrint( "Align Progress: Thread %i starting..." % thread_id, align_printlock )

    with open( log_filename, "w") as log_file:

        LogIt( "INFO: Processing BLAST file '%s'." % blast_file, log_file, -1, lock=align_printlock )

        blast_record_reader = None

        if thread_count <= 1:
            blast_record_reader = FastaIter( blast_file, False, None)
        else:
            #Process only every Nth record in the blast file in this thread
            InitRecordSelector( thread_count, thread_id, 10)
            #.blast file is not a fasta file but FastaIter method still parses the records nicely, one by one
            blast_record_reader = FastaIter( blast_file, False, RecordSelector) #generator

        q_blast_entries = []
        q_fasta_records = []
        queue_retval = -1

        while True:

            #Fill queue, read 10 sequences at a time
            queue_retval = FillAlignmentQueue( blast_record_reader, q_blast_entries, seq_index, q_fasta_records, 10 )
            if queue_retval < 1: break #Exit loop condition

            assert( len( q_blast_entries) == len( q_fasta_records))

            #Process queue
            while len( q_blast_entries):

                #FIFO queue
                blast_header, blast_record = q_blast_entries.pop( 0)
                fasta_rec = q_fasta_records.pop( 0)

                record_id = GetRecordACC( blast_header)
                record_desc = fasta_rec.description

                cutsite =  -1
                res = re_cutnum.search( record_desc)
                #print records[ r].description

                if res:
                    try:
                        cut_1 = int( res.group( 1))
                        cut_2 = int( res.group( 2))
                        cutsite = cut_1
                        #print record_desc #+ " homologs: %i" % len( PDB_homologs[ r])
                    except ValueError:
                        pass

                else:
                    LogIt( "INFO: No POI-site specified for: %s" % (record_desc), log_file, 1, lock=align_printlock)

                n_total_blast_results += 1
                blast_rank = 0 #1st rank is 1, not 0
                n_ok_blast_results = 0
                used_structures = set()

                #For each structure in blast results
                for blast_result in blast_record.split( "\n"):

                    blast_result = blast_result.strip()
                    if len( blast_result) == 0: continue #Skip empty rows

                    blast_row = {}
                    blast_rank += 1
                    #print "Blast result: '''%s'''" % blast_result
                    cols = blast_result.split( "\t")

                    if len( cols) != 12:
                        LogIt( "ALIGN ERROR: Bad Blast result '%s', rank: %i" % (blast_result, blast_rank), log_file, 2, lock=align_printlock )
                        n_errors += 1

                        if n_errors >= 5:
                            SyncErrPrint( "ALIGN ERROR: %i errors encountered in blast results. Exitting..\n" % n_errors, align_printlock)
                            return -2 #RETURN

                        continue

                    blast_row[ 'entry_id'] = record_id #cols[ 0].split( '|')[ 1]

                    #Two possible formats in BLAST results:
                    if cols[ 1].find( "|") >= 0:
                        #gn|P23919|DTYMK    gi|359545626|pdb|2XX3|A 100.00  212 0   0   1   212 21  232 7e-156   436
                        blast_row[ 'pdb_code'] = cols[ 1].split( '|')[ 3].upper() #3IKM
                        blast_row[ 'chain_id'] = cols[ 1].split( '|')[ 4]         #A
                    else:
                        #gn|P30613_Pos2|PKLR     4IP7_A  99.815  541     1       0       34      574     3       543     0.0     1096
                        blast_row[ 'pdb_code'] = cols[ 1].split( '_')[ 0].upper() #3IKM
                        blast_row[ 'chain_id'] = cols[ 1].split( '_')[ 1]         #A

                    blast_row[ 'struct_id']= ( blast_row[ 'pdb_code'] + "|" + blast_row[ 'chain_id'])  #3IKM|A
                    blast_row[ 'pdb_file'] = InsertFolderSeparator( PDB_folder) + "pdb" + blast_row[ 'pdb_code'].lower() + ".ent"

                    if blast_row[ 'pdb_code'] in obsoletes_dict:
                        #Skip obsolete structures
                        LogIt( "INFO: Entry for '%s' contains obsolete structure %s. Skipping." % (blast_row[ 'entry_id'], blast_row[ 'pdb_code']), log_file, 1, lock=align_printlock)
                        continue

                    n_ok_blast_results += 1


                    #Blast results can contain the same structure+chain multiple times
                    if blast_row['struct_id'] not in used_structures:

                        identifier = "%s_%s_%i" % (name, record_id, blast_rank)

                        if cutsite < 0:
                            LogIt( "No POI-site found for seq '%s'. Exitting.\n" % record_id, log_file, 2, lock=align_printlock)
                            return -2 #RETURN
                            #AlignSeqWithPDBFile( identifier, records[ r], blast_entry[PDB_CODE], blast_entry[FILE], blast_entry[CHAIN], output_folder, False, score_it, skip_existing )
                        else:

                            #Needle alignment
                            #align_retval = AlignCutSeqWindowWithNeedle( thread_id=thread_id, cutsite=cutsite, cutseqlen=10, record_identifier=identifier, \
                            #                                            fasta_record=fasta_rec, PDB_code=blast_row[ 'pdb_code'], PDB_filename=blast_row[ 'pdb_file'], \
                            #                                            PDB_chain=blast_row[ 'chain_id'], output_folder=output_folder, skip_existing=skip_existing, log_file=log_file )

                            #Biopython Pairwise2
                            align_retval = AlignCutSeqWindowWithPairwise( thread_id=thread_id, cutsite=cutsite, cutseqlen=10, record_identifier=identifier, \
                                                                          fasta_record=fasta_rec, PDB_code=blast_row['pdb_code'], PDB_filename=blast_row['pdb_file'], \
                                                                          PDB_chain=blast_row['chain_id'], output_folder=output_folder, skip_existing=skip_existing, log_file=log_file )

                            #Align with ClustalW2
                            #align_retval = AlignCutSeqWithPDBFile( thread_id, cutsite, 10, identifier, fasta_rec, blast_row[ 'pdb_code'], blast_row[ 'pdb_file'], \
                            #                                       blast_row[ 'chain_id'], output_folder, False, score_it, skip_existing, log_file=log_file )

                            #debug
                            #print "NEEDLE: %i" % align_retval
                            #sys.exit( -69)

                            if align_retval == 1: align_ok += 1
                            elif align_retval == 0: align_skipped += 1
                            elif align_retval == -98: align_error += 1 #No sequence in PDB file found
                            else:
                                align_error += 1
                                return align_retval

                        #Mark as processed
                        used_structures.add( blast_row[ 'struct_id'])

                #Report progress
                align_printlock.acquire()
                n_blast_entries_processed.value += 1
                n_alignments_done.value += align_ok
                n_alignments_skipped.value += align_skipped
                n_alignments_error.value += align_error
                align_printlock.release()

                align_error = 0
                align_skipped = 0
                align_ok = 0


        #After all processed (or interrupted by error)
        if queue_retval < 0:

                SyncErrPrint( "THREAD ERROR: No matching sequence entry found for blast entry %i, '%s':\n" % (cur_blast_entry_index, record_id))
                SyncErrPrint( "SeqIndex size: %i" % len( seq_index))

                #Sanity check failed
                error_str = "Thread %s error:\n" % thread_id
                template =  "An exception of type {0} occured. Arguments: {1!r}"
                message = template.format( type(ex).__name__, ex.args)
                SyncErrPrint( error_str + message, align_printlock)

                return -2 #RETURN


        if n_ok_blast_results == 0:
            no_struct_found += 1

    return thread_id

def SetAlignPrintLock( aLock):
    global align_printlock
    align_printlock = aLock
    #SetNeedlePrintLock( aLock)
    SetPairwisePrintLock( aLock)

#Aligns sequences (or only cut site sequences) with sequences exported from PDB files.
def ProcessBlastResults( sequence_filename, blast_filename, output_folder, PDB_folder, obsoletes=None, skip_existing=False, score_it=False, threads=0 ):

    global n_alignments_done, n_alignments_skipped, n_alignments_error, align_printlock, obsoletes_dict, fatal_errors #, seq_index

    if not os.path.isdir( PDB_folder):
        sys.stderr.write("ERROR: '%s' in not a valid folder.\n" % PDB_folder)
        return -1 #RETURN

    if not obsoletes:
        obsoletes = FetchObsoletes( True)

    #Create Shared dictionary
    manager = Manager()
    obsoletes_dict = manager.dict()
    for o in obsoletes: obsoletes_dict[ o] = True

    #Open Blast file
    try:
        blast_file = open( blast_filename, 'r')
    except:
        sys.stderr.write("ERROR: Error opening file: '%s'.\n" % blast_filename)
        return -1 #RETURN


    #Read sequences
    #try:
    #    print "Building index of seq records from file '%s'... " % sequence_filename,
    #    seq_index = SeqIO.index( sequence_filename, "fasta", key_function=GetRecordACC)
    #    print "Done."
    #except Exception as ex:
    #    sys.stderr.write("Error getting sequences from file: '%s'.\n" % sequence_filename)
    #    template =       "An exception of type {0} occured. Arguments: {1!r}"
    #    message = template.format( type(ex).__name__, ex.args)
    #    sys.stderr.write( message + "\n")
    #    return -1 #RETURN

    n_records = CountRecords( sequence_filename) #len( seq_index)
    print "Align Progress: Found %i blast entries to align in file '%s'." % (n_records, blast_filename)

    #Figure out how many cpus to use
    n_cpus = cpu_count()
    used_cpus = n_cpus

    if threads < 0:
        used_cpus = n_cpus + threads
    elif threads > 0:
        used_cpus = threads


    #enforce cpu limits
    used_cpus = min( n_cpus, used_cpus) #cap at detected system cores
    if used_cpus > ( n_records / 10): used_cpus = (n_records / 10) #No need for additional cores beyond 10 alignments each
    used_cpus = max( 1, used_cpus) #minimum of 1 cores

    print "Align Progress: Aligning %i records with %i cores." % (n_records, used_cpus)

    output_folder = InsertFolderSeparator( output_folder)
    log_filename = output_folder + "align_log.txt"

    SetAlignPrintLock( Lock())
    SetSyncPrintBuffered( False)
    #manager = Manager()
    #manager.Value('i', 0)

    print "INFO: Skipping existing files: %s" % ("Yes" if skip_existing else "No")
    print "INFO: Creating alignment score files: %s" % ("Yes" if score_it else "No")
    print "Align Progress: Starting %i thread%s and creating sequence record indices..." % (used_cpus, "s" if used_cpus > 1 else "")

    #print "Reading blast file: '%s'" % blast_filename
    args = []
    for c in range( used_cpus):
        args.append( "%s\t%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i" % (blast_filename, sequence_filename, output_folder, used_cpus, c, log_filename, PDB_folder, score_it, skip_existing))

    #Create worker pool
    pool = Pool( processes=used_cpus, initializer=PoolInitializer, initargs=(n_alignments_done, n_alignments_skipped, n_alignments_error, align_printlock, obsoletes_dict, n_blast_entries_processed ))

    for arg in args:
        pool.apply_async( ProcessBlastMultiThread, (arg, ), callback=ReportDone)

    pool.close() # no more tasks

    total_processed = 0

    start_time = time.time()
    prev_progress = 0
    processing_speed_per_minute = []
    report_interval = 20
    progress_counter = max( 5, report_interval - 5) #first report after 5 secs

    try:

        #Trick to keep Ctrl+C interrupt working during processing
        while total_processed < n_records:
            time.sleep( 1)
            progress_counter += 1

            if fatal_errors.value > 0:
                raise Exception('Alignment Thread error.')
            elif fatal_errors.value == -99:
                raise KeyboardInterrupt('Ctrl+C detected.')

            if progress_counter >= report_interval:
                progress_counter = 0
                total_processed = n_blast_entries_processed.value
                #print "ISBUFFERED:", IsSyncPrintBuffered()
                #SyncPrint( "Prev: %i, cur: %i" % (prev_progress, total_processed), False)#align_printlock)

                n_processed = total_processed - prev_progress
                if n_processed == 0: continue
                speed_msg = ""
                if prev_progress != 0 and n_alignments_done.value > 100:

                    per_minute = round(n_processed * (60.0 / report_interval))

                    #Calc average
                    if len( processing_speed_per_minute) > 30: processing_speed_per_minute.pop( 0)
                    processing_speed_per_minute.append( per_minute)
                    avg_speed = int( sum( processing_speed_per_minute) / len( processing_speed_per_minute))

                    n_remaining = n_records - total_processed
                    if avg_speed == 0: avg_speed = 1
                    n_minutes_remaining = n_remaining / avg_speed
                    completion = datetime.now() + timedelta(minutes = n_minutes_remaining)
                    if n_minutes_remaining <= 60:
                        #              ALIGN PROGRESS:
                        speed_msg = "\n                Speed: %i/min. Estimated time of completion: %s (%i minutes remaining)" % (avg_speed, completion.strftime('%H:%M'), n_minutes_remaining)
                    elif n_minutes_remaining < 60*24: #in 24h
                        h = n_minutes_remaining / 60
                        m = n_minutes_remaining - (h*60)
                        speed_msg = "\n                Speed: %i/min. Estimated time of completion: %s (in %ih %imin)" % (avg_speed, completion.strftime('%H:%M'), h, m)
                    else:
                        h = n_minutes_remaining / 60
                        speed_msg = "\n                Speed: %i/min. Estimated time of completion: %s (in >%i hours)" % (avg_speed, completion.strftime('%a %b %d %Y %H:%M'), h)


                SyncPrint( "\n\nALIGN PROGRESS: Entries processed: %i/%i (%2.1f%%), Alignments: %i, Skipped: %i, Errors: %i%s\n\n" % ( total_processed, n_records, (float(total_processed) / n_records * 100.0), n_alignments_done.value, n_alignments_skipped.value, n_alignments_error.value, speed_msg), align_printlock )
                prev_progress = total_processed
            elif progress_counter % 5 == 0: #Update every 5 secs
                total_processed = n_blast_entries_processed.value


        print "\nAlign Progress: Work has finished."
    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        sys.stderr.write( "\n\nProcess interrupted by user.\n\n")
        align_printlock = False
        return -99
    except Exception as ex:
        pool.terminate()
        pool.join()
        sys.stderr.write( "\n\nERROR: Fatal error occurred during alignment.\n\n")
        template = "An exception of type {0} occured. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message +"\n")

        align_printlock = False
        return -1
    else:
        pool.join()  # wrap up current tasks

    print "Align Progress: Alignment took %i seconds" % (int(time.time() - start_time))

    SetAlignPrintLock( False)
    gc.collect()

    #n_found = n_total_blast_results - no_struct_found
    #LogIt( "Homologous PDB structures found for %i/%i sequences (%.1f%%)." % ( n_found, n_total_blast_results, (float(n_found)/ n_total_blast_results * 100.0)), log_file, 1)
    #LogIt( "\n\nINFO: Alignment results (%i total):\nALIGNED=%i\nSKIPPED=%i\nERROR  =%i\n" % (align_ok+align_skipped+align_error, align_ok, align_skipped, align_error), log_file, 1 )

    #Collect Thread logs into one log
    ConcatenateFiles( GetDirFiles( output_folder + "align_*_log.txt"), log_filename, aRemoveAfterWritten=True)


    return 0 #DONE




def PrintHelp():

    print """
    This script runs alignments for a pair of blast results and fasta sequences
    Usage: script.py fasta_sequence_file blast_results_file output_folder
     or: fasta_sequence_folder blast_results_folder output_folder

    -p  --pdb_folder       Folder where structure files are stored
    -h  --help             Print this message.
    -v  --version          1.0
    """

def main():

    seq_file = ""
    blast_file = ""
    output_folder = ""
    pdb_folder = "I:/anSsi/PDB/"

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:hv', ['pdb_folder', 'help', 'version'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)



    #Input & Output files
    for arg in args:
        if len( seq_file) == 0:
            seq_file = arg
        elif len( blast_file) == 0:
            blast_file = arg
        elif len( output_folder) == 0:
            output_folder = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)


    #Flags
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            PrintHelp()
            sys.exit( 0)
        elif opt in ('-v', '--version'):
            PrintHelp()
            sys.exit( 0)
        elif opt in ('-p', '--pdb_folder'):
            pdb_folder = arg
            sys.exit( 0)


    if len( sys.argv) < 3:
        sys.stderr.write("Too few arguments.\n")
        sys.exit( -2) #Exit

    seq_file_list = []
    blast_file_list = []

    if os.path.isdir( seq_file) and os.path.isdir( blast_file):
        #Both dirs
        seq_file_list = GetFastaFilesInFolder( seq_file)
        for f in seq_file_list:
            blast_file_list.append( f + ".blast")

        if len(seq_file_list) < len( blast_file_list):
            sys.stderr.write("Warning: Blast files not found for all sequence (fasta) files.\n")

    elif os.path.isfile( seq_file) and os.path.isfile( blast_file):
        #Both files
        seq_file_list.append( seq_file)
        blast_file_list.append( blast_file)
    elif os.path.isfile( seq_file) and os.path.isdir( blast_file):
        #File and dir
        seq_file_list.append( seq_file)
        blast_file_list.append( InsertFolderSeparator( blast_file) + ".blast")
    else:
        sys.stderr.write("Bad arguments given. Specify either files or folders.\n")


    if len( seq_file_list) > 1:
        print "Found %i files to create alignments for."

    obsoletes = FetchObsoletes( True)

    for s in range( len( seq_file_list)):

        path, filename = os.path.split( seq_file_list[ s])
        name, suffix = SplitFileName( filename)

        output_folder = InsertFolderSeparator( output_folder)

        print "Processing files: '%s' and '%s'." % (seq_file_list[ s], blast_file_list[ s])
        ProcessBlastResults( seq_file_list[ s], blast_file_list[ s], output_folder, pdb_folder, obsoletes )


    print ""
    print "All done."



if __name__ == "__main__":
  main()

