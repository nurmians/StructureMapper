#!/usr/bin/python
# coding=UTF-8
# -*- coding: UTF-8 -*-

#input file of folder of fasta files
#Uses BLAST to search for homologous PDB structures


#    This file is part of asa_uta.py.
#
#    asa_uta is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    asa_uta is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with asa_uta.  If not, see <http://www.gnu.org/licenses/>.

## Authors: Anssi Nurminen, Vesa P. Hytönen
## Institute of Biomedical Technology
## University of Tampere, Tampere, Finland
## and
## BioMediTech, Tampere, Finland
## 14-11-2013

# Please cite the authors if you find this script useful."

import sys
import os
import re
import time
import random
import signal
import glob
import textwrap
import getopt
import logging
import shutil
from datetime import datetime, timedelta
from multiprocessing import Lock, Process

try: from cStringIO import StringIO
except: from StringIO import StringIO

from Bio import SeqIO
from Bio.Application import ApplicationError
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.PDBList import PDBList

#from blast_repository import *
from blast_sqlite3_repository import *

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(  "%s/../utils" % curdir)
from fileops import *
from struct_methods import *#FetchObsoletes, GetRecordIDs, GetRecordACC, SyncPrint, SyncErrPrint
import online_pdb

sys.path.append(  "%s/../align" % curdir)
from align_seq_with_structure_multithread  import ProcessBlastResultsAsync

align_processes = []
blast_print_lock = Lock()
original_sigint = signal.getsignal( signal.SIGINT)

def sigIntHandler( aSignal, aFrame):
    global align_processes, blast_print_lock, original_sigint

    signal.signal( signal.SIGINT, original_sigint) #Restore

    SetSyncPrintBuffered( False)
    SyncPrint( "Process interrupted by user.", blast_print_lock)

    TerminateThreads()

    sys.exit( -99)
    #raise KeyboardInterrupt


def TerminateThreads():
    global align_processes, blast_print_lock

    SetSyncPrintBuffered( False)

    for p in align_processes:
        SyncPrint( "Terminating remaining alignment process %i... " % int( p[ 2]), blast_print_lock, False )
        try:
            p[ 0].terminate()
            SyncPrint( "[OK]", blast_print_lock)
        except:
            SyncPrint( "[Fail]", blast_print_lock)
            pass


def printDictionary( dict):

    print ""

    for k in sorted( dict.iterkeys()):
        print k,
        print ":",
        print dict[ k]

    print ""


def WriteArrayToFile( filehandle, linearray ):

    for line in linearray:
        filehandle.write( line + "\n")


def ReportSpeed( aCurrent, aTotal, aPrintLock=False):

    #Using function attributes
    if not hasattr( ReportSpeed, "prev_time"):
        ReportSpeed.prev_time = time.time()
        ReportSpeed.prev_count = aCurrent
        ReportSpeed.processing_speed_per_minute = []
        return "Speed: Calculating processing speed..."

    now = time.time()
    seconds_passed = int( now - ReportSpeed.prev_time)
    if seconds_passed < 10: return #report interval

    n_processed = aCurrent - ReportSpeed.prev_count
    entries_per_minute = round( n_processed * (60.0 / seconds_passed))

    #Calc average
    if len( ReportSpeed.processing_speed_per_minute) > 19: ReportSpeed.processing_speed_per_minute.pop( 0)
    ReportSpeed.processing_speed_per_minute.append( entries_per_minute)
    avg_speed = int( sum( ReportSpeed.processing_speed_per_minute) / len( ReportSpeed.processing_speed_per_minute))
    if avg_speed == 0: avg_speed = 1

    n_minutes_remaining = (aTotal - aCurrent) / avg_speed

    completion = datetime.now() + timedelta(minutes = n_minutes_remaining)
    if n_minutes_remaining <= 1:
        #            Blast Progress:
        speed_msg = "         Speed: %i/min. Estimated time of completion: %s (almost done)" % (avg_speed, completion.strftime('%H:%M'))
    elif n_minutes_remaining <= 60:
        speed_msg = "         Speed: %i/min. Estimated time of completion: %s (%i minutes remaining)" % (avg_speed, completion.strftime('%H:%M'), n_minutes_remaining)
    elif n_minutes_remaining < 60*24: #in 24h
        h = n_minutes_remaining / 60 #h is float without conversion?
        m = n_minutes_remaining - (h*60)
        speed_msg = "         Speed: %i/min. Estimated time of completion: %s (in %ih %imin)" % (avg_speed, completion.strftime('%H:%M'), h, m)
    else:
        h = n_minutes_remaining / 60
        speed_msg = "         Speed: %i/min. Estimated time of completion: %s (in >%i hours)" % (avg_speed, completion.strftime('%a %b %d %Y %H:%M'), h)

    ReportSpeed.prev_time = now
    ReportSpeed.prev_count = aCurrent

    if aPrintLock: aPrintLock.acquire()

    print "\n--------------------------------------------------------------------------"
    print "Blast Progress: Processed records %i/%i (%.1f%%)\n%s" % (aCurrent, aTotal, (float( aCurrent)/aTotal)*100.0, speed_msg)
    print   "--------------------------------------------------------------------------\n"

    if aPrintLock: aPrintLock.release()

    return speed_msg


def FindStructures( blast_exe_folder, input_file, output_file = "", PDB_save_folder="", BLAST_DB_folder="", parallel_align=True, max_threads=0, max_results=5, e_val_threshold=0.001, fetch_files=True, reset_progress=False, obsoletes=[], verbose=False, sequence_window=30):

    global blast_print_lock

    retval = -1
    SetSyncPrintBuffered( True)
    use_repository = len( BLAST_DB_folder) > 0

    try:
        retval = DoFindStructures( blast_exe_folder, input_file, output_file, PDB_save_folder, BLAST_DB_folder, parallel_align, max_threads, max_results, e_val_threshold, fetch_files, reset_progress, obsoletes, verbose, sequence_window=sequence_window)
    except KeyboardInterrupt:
        if use_repository: CloseRepository()
        SyncErrPrint( "BLAST ERROR: Process interrupted by user.", blast_print_lock)
        TerminateThreads()
    except Exception as ex:
        SyncErrPrint( "BLAST ERROR: Exception caught running FindStructures in '%s'." % (__file__), blast_print_lock)
        SyncErrPrint( "             An exception of type {0} occured. Arguments: {1!r}".format( type( ex).__name__, ex.args), blast_print_lock)
        TerminateThreads()
        #DEBUG
        #raise

        #exc_type, exc_obj, exc_tb = sys.exc_info()
        #fname = os.path.split( exc_tb.tb_frame.f_code.co_filename)[1]
        #print( exc_type, fname, exc_tb.tb_lineno)

    SetSyncPrintBuffered( False)

    return retval

#blastp has a -remote flag for doing queries in online DBs
def DoFindStructures( blast_exe_folder, input_file, output_file = "", PDB_save_folder="", BLAST_DB_folder="", parallel_align=True, max_threads=0, max_results=5, e_val_threshold=0.001, fetch_files=True, reset_progress=False, obsoletes=[], verbose=False, sequence_window=30):

    global align_processes, blast_print_lock

    use_repository = len( BLAST_DB_folder) > 0

    PDB_save_folder = InsertFolderSeparator( PDB_save_folder)
    BLAST_DB_folder = InsertFolderSeparator( BLAST_DB_folder)
    BLAST_EXE_FOLDER = InsertFolderSeparator( blast_exe_folder)
    repository_file = BLAST_DB_folder + "blast_repository.sqlite"

    if use_repository:
        print "INFO: Using repository of previous results: '%s'" % repository_file
        SetRepositoryFile( repository_file)

    sequence_window += sequence_window % 2 #Make window even
    if sequence_window > 0:
        if verbose: sys.stderr.write( "INFO: Using sequence window of %i amino acids around each POI for BLASTing.\n" % sequence_window )
        if sequence_window < 10: sys.stderr.write( "WARNING: Blasting with a small sequence window (%i) might yield unexpected results.\n" % sequence_window )
    elif verbose: sys.stderr.write( "INFO: Using full sequences for BLASTing.\n" % sequence_window )

    #if not specified, limit to 8 threads
    if max_threads == 0: max_threads = 8

    #0 or less means all possible results
    if max_results <= 0: max_results = 999

    if verbose and max_results > 0: sys.stderr.write( "INFO: Using maximum of %i homologous structures for each POI.\n" % max_results )
    elif verbose: sys.stderr.write( "INFO: Using all available homologous structures for each POI.\n" )

    #Reserve 1 thread for blasting and the rest for alignments
    if parallel_align and max_threads == 1:
        sys.stderr.write( "INFO: Cannot use only one thread for parallel alignments.\n" )
        sys.stderr.write( "      No parallel alignments are done.\n" )
        parallel_align = False
    elif max_threads > 1:
        max_threads -= 1

    if parallel_align: print "INFO: using %i threads for parallel alignments." % max_threads

    #Terminate parallel alignment threads in own handler
    if parallel_align:
        signal.signal(signal.SIGINT, sigIntHandler)

    #Format default output filename if necessary
    if len( output_file) == 0: output_file = input_file + ".blast"

    #Check folder
    if not os.path.isdir( PDB_save_folder):
        sys.stderr.write( "ERROR: '%s' does not specify a valid directory.\n" % PDB_save_folder)
        return -1

    print "BLAST progress: Indexing file '%s'..." % (input_file)
    record_index = SeqIO.index( input_file, "fasta", key_function=GetRecordACC)
    n_records = len( record_index)

    print "BLAST progress: Read %i sequence records from file '%s'" % (n_records, input_file)

    if n_records <= 0:
        sys.stderr.write( "ERROR: No records to process.\n")
        return -1

    if len( obsoletes) == 0:
        obsoletes = online_pdb.FetchObsoletes( True, PDB_save_folder);


    seq_record_ids = GetRecordIDsFromFile( input_file, aWithPos=True)
    missing_blast_record_ids = seq_record_ids
    n_progress = 0

    #Sanity check
    if n_records != len( missing_blast_record_ids):
        sys.stderr.write( "ERROR: Missmatching number of sequence IDs. Quitting.\n")
        return -4

    #Get current progress
    if not reset_progress:
        #Get sequense ids from blast output file
        #Check if work is already done
        if os.path.isfile( output_file):

            blast_record_ids = GetRecordIDsFromFile( output_file, aWithPos=True)

            #blast_record_ids = []
            #with open( output_file, "r") as opf:
            #    for line in opf:
            #        if line[ 0] == ">":
            #            blast_record_ids.append( GetRecordACC( line[1:])) #skip ">"
            missing_blast_record_ids = list( set( seq_record_ids) - set( blast_record_ids))

        n_progress = n_records-len( missing_blast_record_ids)

        if n_progress > 0:
            print "BLAST progress: %i/%i (%.1f%%) records already BLASTed." % (n_progress, n_records, (float( n_progress)/n_records*100.0))
            print "BLAST progress: Remove file '%s' to re-process." % output_file
    elif os.path.isfile( output_file):
        #Reset
        #Remove existing output file (if any)
        try:
            print "BLAST progress: Resetting any previous BLAST progress... ",
            os.remove( output_file)
            print "Done."
        except OSError:
            sys.stderr.write("\nERROR: Could not remove BLAST progress file: '%s'. Remove it to reset BLAST progress." % output_file)
            sys.exit(-2)


    n_missing = len( missing_blast_record_ids)

    if n_missing > 0:
        #Open output file for writing
        print "BLAST progress: Writing output to file '%s'." % output_file
    else:
        print "BLAST progress: All Done."
        RemoveParallelAlignFiles( GetFolder(output_file))
        return 0


    fetched_pdbs = {}

    cur_result_chunk_file = None
    cur_result_chunk_filename = None
    chunk_size = 20
    file_written = False
    #cur_chunk = 0
    #align_processes = []
    align_file_queue = []
    align_max_queue = max( 10, max_threads*2 + 1)
    align_errors = 0
    #blast_paralign_%i.blast

    #Create align folder
    align_folder = InsertFolderSeparator( GetFolder( output_file) + "align" )
    if not os.path.exists( align_folder):  os.makedirs( align_folder)

    #Remove previous parallel alignment files (if any)
    RemoveParallelAlignFiles( GetFolder(output_file))


    with open( output_file, "a+") as results_file:

        previous_sequence = ""
        previous_id_stub = ""

        lines = []
        n_lines = 0
        line_ok = []
        n_blasted = 0

        #re_pos = re.compile( "_Pos\d+")

        for i in range( max(0, n_progress - 10), n_records):

            #print "I:%i" % i
            seq_id = seq_record_ids[ i]
            #cur_record = records[ i]
            cur_record = record_index[ seq_id]

            poipos = GetRecordPOI( cur_record.description)

            cur_seq = ""

            if sequence_window > 0:
                #Blast with partial (windowed) sequence
                start_pos = poipos - (sequence_window / 2)
                end_pos = start_pos + sequence_window
                if start_pos < 0:
                    start_pos = 0
                    end_pos = sequence_window
                if end_pos > len( cur_record.seq):
                    end_pos = len( cur_record.seq)
                    start_pos = max( 0, end_pos - sequence_window)
                cur_seq = str( cur_record.seq[ start_pos:end_pos])
            else:
                #Blast with full sequence
                cur_seq = str( cur_record.seq)


            same_id_as_previous = False
            fetched_from_repository = False

            if sequence_window <= 0:
                #For full sequence blasting, same id means same results
                same_id_as_previous = len( previous_id_stub) and previous_id_stub == StripPosTag( seq_id) #re.sub("_Pos\d+", "", seq_id) #StripPosTag


            #Sanity check
            if i < n_progress and seq_id in missing_blast_record_ids:
                sys.stderr.write( "ERROR: Number of processed entries (%i) does not match with current blast file.\n" % n_progress)
                sys.stderr.write( "ERROR: Entry '%s' should already be processed according to its index in the file. \n" % seq_id)
                sys.stderr.write( "ERROR: Consider restarting BLAST processing by removing the file '%s'\n" % output_file)
                sys.stderr.write( "       or using the '--reset' flag.\n")
                sys.exit( -3)


            if seq_id not in missing_blast_record_ids:
                msg = "Entry %i: %s already blasted." % (i, seq_id)
                print msg.ljust( 55), "[OK]"
                continue #Already fetched

            #DEBUG
            #PrintEntry( cur_seq)
            #sys.exit( 0)

            #If the sequence is exactly the same as for the previous fasta entry, use those results
            if cur_seq and ( same_id_as_previous or (len( cur_seq) == len( previous_sequence) and str(cur_seq) == previous_sequence)):
                #Use previous results
                if verbose: SyncPrint( "Blast search already processed for %s (#%i), using existing results." % (cur_record.id, i), blast_print_lock)
            elif not cur_seq or len(cur_seq) == 0:
                SyncErrPrint( "ERROR: No sequence for %s (#%i)." % (cur_record.id, i), blast_print_lock)
                sys.exit( -2)
            else:

                #print "\nProcessing record %i/%i (%.1f%%): %s" % (n_progress+1, n_records, ((n_progress+1.0)/n_records)*100.0, cur_record.id)
                previous_sequence = str(cur_seq)
                previous_id_stub = StripPosTag( seq_id) #re.sub("_Pos\d+", "", seq_id)

                #blastp -query "A:\3IKM_a.fasta" -db pdbaa -outfmt 6
                lines = None
                if use_repository: lines = FetchFromRepository( cur_seq, -1 if max_results > 100 else max_results, 30)

                SyncPrint( "\nBLAST:\n[%s] Processing record '%s'" % ("+" if lines != None else "-", cur_record.id), blast_print_lock, aBypassBuffer=True)

                if lines != None:
                    #Found in repository
                    #lines = FetchFromRepository( seq_id, repository)[ 1:] #Skip header row
                    #SyncPrint(  "[+] Result for '%s' found in repository." % (seq_id), blast_print_lock, aBypassBuffer=True )

                    if verbose and sequence_window <= 100: SyncPrint( "    Seq: '%s'" % cur_seq )
                    fetched_from_repository = True
                else:

                    #WRITE CURRENT SEQUENCE INTO temp FILE for blasting
                    with open( 'current.fasta', 'w') as tempf:
                        tempf.write( ">%s\n" % cur_record.description)
                        tempf.write( "%s\n" % cur_seq)

                    #SeqIO.write( cur_record, 'current.fasta', "fasta")

                    BlastX_Command = NcbiblastxCommandline( cmd=blast_exe, query="-", out="-", outfmt=6, db=("I:/anSsi/Blast_DB/"+'pdbaa'), evalue=0.001, max_target_seqs=5)
                    blast_output = StringIO( BlastX_Command( stdin=">Test_seq\nTDKPVYTPDQSVKIRVYSLSDDLKPAKRETVLT\n")[ 0])
                    #blast_output = StringIO( BlastX_Command( stdin=">Test_seq\nAAAAAAAAA\n")[ 0])

                    #for record in NCBIXML.parse( blast_handle):
                    bstd = blast_output.getvalue()


                    #BLAST IT
                    #SyncPrint( "Blasting...", blast_print_lock, aBypassBuffer=True )
                    retry_blast = 3

                    while retry_blast:

                        errors = False
                        if verbose and sequence_window <= 100: SyncPrint( "Blasting seq: '%s'" % cur_seq )
                        #blastp_cline = NcbiblastxCommandline( cmd=BLAST_EXE_FOLDER+'blastp', out='out.blast', outfmt=6, query='current.fasta', db=(BLAST_DB_folder+'pdbaa'), evalue=e_val_threshold, max_target_seqs=max_results)
                        blastp_cline = NcbiblastxCommandline( cmd=BLAST_EXE_FOLDER+'blastp', out='out.blast', outfmt=6, query='current.fasta', db=(BLAST_DB_folder+'pdbaa'), evalue=e_val_threshold, max_target_seqs=max_results)

                        try:
                            blastp_cline()
                        except ApplicationError as ex:
                            #Sometimes blast cannot find "pdbaa.00.phr" even though it is there
                            SyncErrPrint( "ERROR: Exception caught running NcbiblastxCommandline in '%s'.\n" % (__file__), blast_print_lock)
                            SyncErrPrint( "       An exception of type {0} occured. Arguments: {1!r}".format(type(ex).__name__, ex.args), blast_print_lock )
                            errors = True
                            retry_blast -= 1
                            time.sleep( 3)
                        except Exception as ex:
                            errors = True
                            SyncErrPrint( "ERROR: Exception caught running NcbiblastxCommandline in '%s'.\n" % (__file__), blast_print_lock)
                            SyncErrPrint( "       An exception of type {0} occured. Arguments: {1!r}".format(type(ex).__name__, ex.args), blast_print_lock )
                            sys.exit( -66)

                        if errors:
                            if retry_blast > 0:
                                SyncErrPrint( "Retrying blast...", blast_print_lock )
                            elif retry_blast == 0:
                                SyncErrPrint( "Retries exhausted. Exitting.", blast_print_lock )
                                sys.exit(-67)
                        else:
                            retry_blast = 0 #exit loop

                    n_blasted += 1

                    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

                    #read blast output
                    f = open('out.blast')
                    lines = map( str.strip, f.readlines())
                    f.close()

                n_lines = len( lines)
                line_ok = [True]*n_lines

                if n_lines < 1:
                    SyncPrint( "NO RESULTS", blast_print_lock, aBypassBuffer=True )
                else:

                    #print lines
                    #for line in lines: print "'%s'" % line
                    for line in lines: print line


                    #PDB
                    #p = PDBParser( PERMISSIVE=1)

                    for l in range( n_lines):

                        cols = lines[ l].split( "\t" );

                        if len( cols) != 12:
                            SyncErrPrint( "BLAST ERROR: line error: '%s'" % lines[ l], blast_print_lock )
                            line_ok[ l] = False
                            continue

                        #Two possible formats in BLAST results
                        if cols[ 1].find( "|") >= 0:
                            #gn|P23919|DTYMK    gi|359545626|pdb|2XX3|A 100.00  212 0   0   1   212 21  232 7e-156   436
                            id_col = cols[ 1].split( "|")
                            pdb_id = id_col[ 3]
                        else:
                            #gn|P30613_Pos2|PKLR     4IP7_A  99.815  541     1       0       34      574     3       543     0.0     1096
                            id_col = cols[ 1].split( "_")
                            pdb_id = id_col[ 0]


                        if pdb_id in obsoletes:
                            SyncPrint( "%s is obsolete." % pdb_id, blast_print_lock, aBypassBuffer=True )
                            #Get rid of obsoletes
                            line_ok[ l] = False
                            continue
                        else:
                            if len( pdb_id) != 4:
                                SyncErrPrint( "PDB ID ERROR:" + pdb_id, blast_print_lock )
                            elif fetch_files and pdb_id not in fetched_pdbs:
                                #Get file from online RSCB database
                                online_pdb.FetchPDBFile( pdb_id, PDB_save_folder, False, obsoletes, verbose )
                                fetched_pdbs[ pdb_id] = True



            #For parallel align processing
            if parallel_align:

                #print "Clearing buffer..."
                #SyncPrintFlushBuffer( blast_print_lock, aTitle="Alignments:\n")

                try:
                    new_chunk_filename = GetFolder( output_file) + "blast_paralign_%i.blast" % (i - (i % chunk_size))
                    #New chunk of BLAST results
                    if new_chunk_filename != cur_result_chunk_filename:

                        #Close current file
                        if cur_result_chunk_file: cur_result_chunk_file.close()

                        #Remove processed entries
                        CleanParallelAlignFiles( align_processes)

                        queue_full = len( align_file_queue) >= align_max_queue

                        #Insert current file to queue, if queue not full
                        if cur_result_chunk_file:

                            if not queue_full and file_written: #Max size for file queue
                                SyncErrPrint( "INFO: Queueing file for alignment: '%s'" % cur_result_chunk_filename, blast_print_lock )
                                align_file_queue.append( cur_result_chunk_filename)
                            else:
                                #Queue full, or file empty?
                                #File not needed
                                try: os.remove( cur_result_chunk_filename)
                                except: SyncErrPrint( "WARNING: Could not remove previous blast result file '%s'." % cur_result_chunk_filename, blast_print_lock )


                        #Start new alignment processes for first file in queue
                        while len( align_file_queue) > 0 and len( align_processes) < max_threads:

                            next_align_file = align_file_queue.pop( 0)

                            if not os.path.isfile( next_align_file):
                                SyncErrPrint( "\n\nWARNING: File '%s' does not exist." % next_align_file, blast_print_lock )
                                continue
                            if int( os.stat( next_align_file).st_size) <= 0:
                                SyncErrPrint( "\n\nWARNING: File '%s' is empty." % next_align_file, blast_print_lock )
                                continue

                            process_id = 0
                            if len( align_processes): process_id = min( set( range(len( align_processes)+1))-set( zip(*align_processes)[ 2])) #Lowest unused integer
                            if verbose: SyncPrint( "Align Progress: Starting new alignment process for entries %i-%i... (id:%i)" % (( i-(i % chunk_size)), i-1, process_id), blast_print_lock, aBypassBuffer=True )
                            process = ProcessBlastResultsAsync( next_align_file, input_file, align_folder, process_id, PDB_save_folder, obsoletes, aScoreIt=False, aSkipExisting=True, aLock=blast_print_lock)
                            align_processes.append( (process, next_align_file, process_id))

                        #open file for next chunk
                        cur_result_chunk_filename = new_chunk_filename

                        #Stop writing new files if queue already full
                        if not queue_full: cur_result_chunk_file = open( new_chunk_filename, "a+")
                        else: cur_result_chunk_file = None
                        file_written = False #Making sure not to add empty files to queue

                except Exception as ex:
                    SyncErrPrint( "WARNING: Error with parallel alignment processing:", blast_print_lock )
                    SyncErrPrint( "         An exception of type {0} occured. Arguments: {1!r}\n".format( type(ex).__name__, ex.args), blast_print_lock )
                    cur_result_chunk_file = None
                    align_errors += 1
                    if align_errors > 3:
                        SyncErrPrint( "ERROR: Too many parallel alignment errors. Use '-n' to not use parallel aligning during BALST queries. Exitting.", blast_print_lock )
                        sys.exit(-68)



            #Disable keyboard interrupt during writing of result file(s)
            sig = signal.signal(signal.SIGINT, signal.SIG_IGN)

            blast_record = []
            blast_record_header = ">%s\n" % cur_record.id
            blast_record.append( blast_record_header)

            results_file.write( blast_record_header.rstrip() + (", [%i]\n" % i) )
            if cur_result_chunk_file: cur_result_chunk_file.write( blast_record_header.rstrip() + (", [%i]\n" % i) )

            structures_found = 0

            for l in range( n_lines):
                if line_ok[ l]: #Get rid of obsoletes
                    blast_record.append( lines[ l])
                    results_file.write( "%s\n" % lines[ l])
                    if cur_result_chunk_file: cur_result_chunk_file.write( "%s\n" % lines[ l])
                    structures_found += 1

            results_file.write( "\n" )
            if cur_result_chunk_file: cur_result_chunk_file.write( "\n" )
            results_file.flush()
            if cur_result_chunk_file: cur_result_chunk_file.flush()
            file_written = True

            if not same_id_as_previous and not fetched_from_repository and use_repository: AddToRepository( cur_seq, blast_record[1:], max_results ) #AddToRepository( seq_id, repository, repository_file, blast_record)

            #Log
            logging.info(('STRUCT_COUNT:%02i: ' % structures_found) + cur_record.id )

            #Re-enable keyboard interrupt
            signal.signal(signal.SIGINT, sig)

            n_progress += 1
            SyncPrintFlushBuffer( blast_print_lock, aTitle="\nAlignments:\n")
            if not fetched_from_repository: ReportSpeed( n_progress, n_records, blast_print_lock)


    CloseRepository()

    if cur_result_chunk_file: cur_result_chunk_file.close()
    CleanParallelAlignFiles( align_processes)

    if len( align_processes):
        SyncPrint( "Waiting for 180 secs for remaining %i alignment process%s to finnish..." % (len( align_processes), ("" if len( align_processes) == 1 else "es")), blast_print_lock, aBypassBuffer=True )
        timer = 0
        while timer < 180:
            time.sleep( 2)
            SyncPrintFlushBuffer( blast_print_lock)
            timer += 2
            alive = 0
            for ap in align_processes:
                if ap[ 0].is_alive(): alive += 1
            if alive == 0:
                print "INFO: All threads have finished."
                break


        CleanParallelAlignFiles( align_processes)
        for p in align_processes:
            SyncPrint( "Terminating process %i... " % p[ 2], blast_print_lock, False, aBypassBuffer=True )
            p[ 0].terminate()
            SyncPrint( "Done.", blast_print_lock, aBypassBuffer=True )

    print "BLASTED %i entries.\nFile: '%s' written.\n" % ( n_blasted, output_file )

    RemoveParallelAlignFiles( GetFolder(output_file))
    signal.signal( signal.SIGINT, original_sigint)


    return 0 #All OK

def RemoveParallelAlignFiles( aBlastFolder):
    global blast_print_lock

    #Remove previous parallel alignment files (if any)
    retval, files = GetFolderFiles( aBlastFolder, "blast_paralign_*.blast")

    if retval == 0 and len( files):
        #align_file_queue.append( ChangeToFolder( files.pop( 0), align_folder))
        for cf in files:
            try:
                os.remove( ChangeToFolder( cf, aBlastFolder))
            except:
                SyncErrPrint( "WARNING: could not remove previous blast result file '%s'." % cf, blast_print_lock)

def CleanParallelAlignFiles( align_processes):

    global blast_print_lock

    for ap in reversed( align_processes):
        if not ap[ 0].is_alive():
            try:
                thread_id = ap[ 2]
                file_to_be_removed = ap[ 1]
                exit_code = ap[ 0].exitcode
                align_processes.remove( ap)
                SyncPrint( "Align Progress: File '%s' processed in thread %i.  %s" % (file_to_be_removed, thread_id, ("[OK]" if exit_code == 0 else ("[FAIL-%i]" % exit_code))), blast_print_lock)

                if os.path.isfile( file_to_be_removed): os.remove( file_to_be_removed)
            except (OSError, WindowsError) as oex:
                SyncErrPrint( "\nWARNING: Could not delete parallel alignment file '%s'." % file_to_be_removed, blast_print_lock)
                SyncErrPrint( "         An exception of type {0} occured. Arguments: {1!r}".format( type( oex).__name__, oex.args), blast_print_lock)
            except Exception as ex:
                SyncErrPrint( "\nWARNING: Could not clean parallel align files.")
                SyncErrPrint( "         An exception of type {0} occured. Arguments: {1!r}".format( type( ex).__name__, ex.args), blast_print_lock)

def GetFolder( filepath):
    path, filename = os.path.split( filepath)
    return InsertFolderSeparator( path)

def ChangeToFolder( filepath, folder):

    path, filename = os.path.split( filepath)
    return InsertFolderSeparator( folder) + filename

def PrintHelp():

    print """
          Usage: blast_structure_finder.py [flags] input_file output_file


          -e  --evalue FLOAT     E-value threshold for BLAST queries. default = 0.001
          -m  --max INT          Max number of results per file.
          -d  --db_folder PATH   Path and name to local BLAST database.
          -s  --save_folder PATH Specify where fetched PDB files are saved.
          -r  --reset            Reset progress and reprocess files from the beginning
          -n  --nofetch          Do not fetch structures.
          -f  --fetch            Fetch result structures from PDB database [default]
          -h  --help             Print this message.
          -v  --version          1.0
          """

def main():

    input_file = ""
    output_file = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'e:m:d:s:rnfhv', ['evalue=', 'max=', 'db_folder=', 'save_folder=', 'reset', 'nofetch', 'fetch', 'help', 'version'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        PrintHelp()
        #sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)

    if len( sys.argv) < 2:
        sys.stderr.write("Too few arguments.\n")
        PrintHelp()
        sys.exit( -2) #Exit

    #Input & Output files
    for arg in args:
        if len( input_file) == 0:
            input_file = arg
        elif len( output_file) == 0:
            output_file = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)


    file_list = []
    outfile_list = []

    if os.path.isfile( input_file):

        #Input File
        print "MAIN: Processing file '%s' using BLAST to find homologous structures." % input_file
        file_list.append( input_file)

        if len( output_file) and os.path.isdir( output_file):
            outfile_list.append( ChangeToFolder( input_file, output_file) + ".blast")
        else:
            outfile_list.append( input_file + ".blast")

    elif os.path.isdir( input_file):

        input_file = InsertFolderSeparator( input_file)

        #Input Folder
        file_list = GetFastaFilesInFolder( input_file)
        print "MAIN: Processing all fasta files in folder: '%s'" % input_file
        print "MAIN: Found %i fasta files." % len( file_list)

        if len( output_file):
            if not os.path.isdir( output_file):
                sys.stderr.write("'%s' does not specify a valid folder.\n" % output_file)
                sys.exit( 2) #Exit
            else:
                output_file = InsertFolderSeparator( output_file)
                for f in file_list:
                    outfile_list.append( ChangeToFolder( f, output_file) + ".blast")
        else:
            for f in file_list:
                outfile_list.append( f + ".blast")
    else:
        sys.stderr.write("'%s' does not specify a valid file or folder.\n" % input_file)
        sys.exit( -2) #Exit

    #e:m:d:s:
    e_val = 0.001
    max_results = 5
    db="I:/anSsi/Blast_DB/"
    save_path="I:/anSsi/PDB/"
    reset=False
    fetch=True

    #Flags
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            PrintHelp()
            sys.exit( 0)
        elif opt in ('-v', '--version'):
            PrintHelp()
            sys.exit( 0)
        elif opt in ('-r', '--reset'): reset = True
        elif opt in ('-f', '--fetch'): fetch=True
        elif opt in ('-n', '--nofetch'): fetch=False
        elif opt in ('-d', '--db_folder'):
            db = arg
            if arg and arg.upper() == "NONE": db = ""
        elif opt in ('-s', '--save_folder'): save_path = arg
        elif opt in ('-e', '--evalue'):
            try:
                e_val = float( arg)
            except ValueError:
                sys.stderr.write( "Bad evalue given: '%s'. Using default: 0.001\n" % arg)
        elif opt in ('-m', '--max'):
            try:
                max_results = int( arg)
            except ValueError:
                sys.stderr.write( u"Bad max_results given: '%s'. Using default: 4\n" % arg)
        else:
            sys.stderr.write("Unknown option: '%s'.\n" % opt)
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)

    if len( file_list) < 1:
        print "No files to parse."
        sys.exit( 0)



    #Find structures
    for i in range( len( file_list)):

        #def FindStructures( blast_exe_folder, input_file, output_file = "", PDB_save_folder="", BLAST_DB_folder="", parallel_align=True, max_threads=0, max_results=5, e_val_threshold=0.001, fetch_files=True, reset_progress=False, obsoletes=[], verbose=False, sequence_window=30):
        FindStructures( blast_exe_folder="C:/blast-2.4.0+/bin", input_file=file_list[ i], output_file=outfile_list[ i], PDB_save_folder=save_path, BLAST_DB_folder=db, parallel_align=False,
                        max_results=max_results, e_val_threshold=e_val, fetch_files=fetch, reset_progress=reset, verbose=True, sequence_window=20)



    print ""
    print "All done."



if __name__ == "__main__":


    #try: from cStringIO import StringIO
    #except: from StringIO import StringIO

    from StringIO import StringIO

    out_file = "test.out"
    #handle = StringIO( "EIGEKGVNLSGGQK")
    #data = handle.getvalue()

    blast_exe = "C:/Program Files/NCBI/blast-2.2.28+/bin/blastp.exe"
    #blast_exe = "C:/Program Files/NCBI/blast-2.6.0+/bin/blastp.exe"
    #blast_exe = "C:/blast-2.2.28+.abk/bin/blastp.exe"
    #blast_exe = "C:/blast-2.4.0+/bin/blastp.exe"

    if not os.path.isfile( blast_exe):
        blast_exe = "C:/Program Files/NCBI/blast-2.6.0+/bin/blastp.exe"
    #   blast_exe = "C:/blast-2.4.0+/bin/blastp.exe"

    if not os.path.isfile( blast_exe):
        sys.stderr.write( "ERROR: Blast executable not found.\n")
        sys.exit( -1)

    try: os.remove( out_file)
    except: pass

    blast_output = StringIO()
    blast_input = StringIO("MGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\r\nMGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\r\n")
    inputseq = """>gn|A0A096P6L9_Pos74|C5 OS=rat GN=C5|POI:74-75|POISEQ:GYVNLS/PENKFQ
                  MGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGYTEAFDATISLKSYPDKKVTYSSGYVNLSPENKFQNSALLTLPPKQFPRDENPVSHVYLEVVSMHFSKSKKIPITYDNGFLFIHTDKPVYTPDQSVKIRVYSLSDDLKPAKRETVLTFVDPEGTEVDIVEENDYTGIISFPDFKIPSNPKYGVWTIKAKYKKDFTTTGTAYFEVKEYVLPRFSVSIEPESNFIGYKNFKNFEITVKARYFYNKMVPDAEVYIFFGLREDIKEDEKQMMHKAMQAATLMDGAAQTCWLAETEVREINYNMFEDRNNNYSYIACTVTESSGGFSEEAEIPGIKYVLSPYTLNLVATPLFLKPGIPFSIKVQVKDSLEQLVGGVPVTLMAQTVNVNQETSDLEPKRSITHSADGVASFVVNLPSEVTSLKFEVKTDAPELPEENQASKEYEAVTYSSLSQSYIYIGWTENYKPMLVGEYLNIIVTPKSPYIDKITHYNYLILSKGKIVQYGTKEKLLYSSYQNINIPVTQDMVPSARLLVYYIVTGEQTAELVADAVWINIEEKCGNQLQVHLSPDKDVYSPGQTVSLDMVTEADSWVALSAVDSAVYGVRGKAKRAMQRVFQAFDDKSDLGCGAGGGRDNVDVFHLAGLTFLTNANADDSQYHDDSCKEILRPKRDLQLLHQKVEEQAAKYKHRVPKKCCYDGARENKYETCEQRVARVTIGPHCIRAFNECCTIADKIRKESHHKGMLLGRIQIKALLPVMKAEIRSYFPESWLWEVHRVPKRNQLQVALPDSLTTWEIQGIGISDNGICVADTLKAKVFKDVFLEMNIPYSVVRGEQIQLKGTVYNYRTSGTMFCVKMSAVEGICTPGSSAASPQTSRSSRCVRQRIEGSSSHLVTFSLLPLEIGLHSINFSLETSFGKEILVKTLRVVPEGIKRESYAGVTLDPRGVYGIVNRRKEFPYRIPLDLVPKTNVKRILSVKGLLIGEFLSTVLSKEGIDILTHLPKGSAEAELMSIVPVFYVFHYLEAGNHWNIFHPDTLARKQSLQKKIKEGLVSVMSYRNADYSYSMWKGASSSAWLTAFALRVLGQVNKYVKQDQYSICNSLLWLIEKCQLENGSFKENSQYLPIKLQGTLPAEAQENTLYLTAFSVIGIRKAIGICPTEKIYTALAKADSFLLERTLPSKSTFTLAIVAYALSLGDRTHPKFRSIVSALKREALVKGDPPIYRFWRDTLQRPDSSAPNSGTAGMVETTAYALLTSLNLKETSYVNPIIKWLSEEQRYGGGFYSTQDTINAIEGLTEYSLLVKQLHLDMDINVSYKHKGDFYQYKVTEKNFLGRPVEVPLNDDLIVTTGYSSGLATVYVKTVVHKTSVAEEFCSFYLKIDTQEVEASSYLSYSDSGHKRIIACASYKPSKEESASGSSHAVMDILLPTGIGANQEDLRALVEGVDQLLTDYQIKDSHVILQLNSIPSRDFLCVRFRIFELFQVGFLNPATFTVYEYHRPDKQCTMIYSTSDTNLQRVCEGAACKCVEADCGQLQAELDLAISADTRKETACKPEIAYAYKVRITSATEENIFVKYTATLLDIYKTGEAAAEKDSEITFIKKISCTNANLVKGKQYLIMGKEALQIKHNFSFKYIYPLDSSTWIEYWPTDTTCPSCQAFVANLDEFAEDIFLNGCE"""


    ################
    #    TEST 1    #
    ################

    import subprocess
    from Bio.Blast.Applications import NcbiblastpCommandline

    #query = 'NNAGFLD\nSNLIIVLNDN'  #your string from some external source
    query = ">Test_seq\nMIVSDIEANALLESVTKFHCGVIYDYSTAEYVSYRPSDFGAYLDALEAEVARGGLIVFHNGHKYDVPALT\n"
    #query = "MIVSDIEANALLESVTKFHCGVIYDYSTAEYVSYRPSDFGAYLDALEAEVARGGLIVFHNGHKYDVPALT\nMIVSDIEANALLESVTKFHCGVIYDYSTAEYVSYRPSDFGAYLDALEAEVARGGLIVFHNGHKYDVPALT\n\n\r\n\r\n"
    #blastp_cline = NcbiblastpCommandline( db="pdbaa", outfmt=6) #format Blast command
    blastp_cmdline = NcbiblastxCommandline( cmd=blast_exe, query='test_input.fasta', db="pdbaa", outfmt=6) #format Blast command
    print "CMD:", str(blastp_cmdline)

    #s = StringIO( query)
    #sys.stdin = s
    #process = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) #setup the process
    process = subprocess.Popen( str(blastp_cmdline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=False, shell=True) #setup the process
    #process.stdin.write( 'MIVSDIEANALLESVTKFHCG\n\rMIVSDIEANALLESVTKFHCG\n\r')
    #process.stdout.read()

    out, err = process.communicate( input='MIVSDIEANALLESVTKFHCGVIYD') #run the process
    #process.wait()
    #print process.stdout.read()


    print "OUT (%i):" % len( out), out
    print "ERR (%i):" % len( err), err
    print "RETCODE:", process.returncode
    sys.exit( 0)


    ################
    #    TEST 2    #
    ################

    #query = "MIVSDIEANALLESVTKFHCGVIYDYSTAEYVSYRPSDFGAYLDALEAEVARGGLIVFHNGHKYDVPALT"
    query = ">Test_seq\nMIVSDIEANALLESVTKFHCGVIYDYSTAEYVSYRPSDFGAYLDALEAEVARGGLIVFHNGHKYDVPALT\n"

    #blastp_cline = NcbiblastxCommandline( cmd=blast_exe, query="test_input.fasta", out="-", outfmt=6, db=("I:/anSsi/Blast_DB/"+'pdbaa'), evalue=0.001, max_target_seqs=5)
    #blastp_cline = NcbiblastxCommandline( cmd=blast_exe, query="-", out="-", outfmt=6, db=("I:/anSsi/Blast_DB/"+'pdbaa'), evalue=0.001, max_target_seqs=5)
    blastp_cline = NcbiblastxCommandline( cmd=blast_exe, outfmt=6, query="-", out="-", db="I:/anSsi/Blast_DB/pdbaa", evalue=0.001, max_target_seqs=5)
    print "CMD:", blastp_cline


    print >>blast_output, blastp_cline( stdin=query)#[ 0]
    #print >>blast_output, blastp_cline( stdin=inputseq)[ 0]
    bstd = blast_output.getvalue()
    print "OUTPUT:"
    print bstd
    sys.exit( 0)


    ################
    #    TEST 3    #
    ################

    #BlastX_Command = NcbiblastxCommandline( cmd=blast_exe, query="-", out="-", outfmt=6, db=("I:/anSsi/Blast_DB/"+'pdbaa'), evalue=0.001, max_target_seqs=5)
    BlastX_Command = NcbiblastxCommandline( cmd=blast_exe, query="-", out="-", outfmt=6, db="I:/anSsi/Blast_DB/pdbaa", evalue=0.001, max_target_seqs=99)
    print "CMD:", BlastX_Command

    #print >>blast_output, BlastX_Command( stdin=">Test_seq\nTDKPVYTPDQSVKIRVYSLSDDLKPAKRETVLT\n")[ 0]
    #print >>blast_output, BlastX_Command( stdin=inputseq)[ 0]
    s = StringIO("MGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\nMGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\n")
    sys.stdin = s.getvalue()
    #r = raw_input('What you say?\n')


    #print >>blast_output, BlastX_Command( stdin=blast_input)[ 0]
    print >>blast_output, BlastX_Command(stdin="MGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\nMGLWGLLCLLIFLDKTWGQEQTYVISAPKIFRVGSSENVVIQAHGY\n")[ 0]
    sys.stdin = sys.__stdin__

    #print >>blast_output, BlastX_Command( stdin=">Test_seq\nHTGFLTEYVATRWYRAP\n")[ 0]

    #for record in NCBIXML.parse( blast_handle):
    bstd = blast_output.getvalue()
    #print blast_output.getvalue()
    print "OUTPUT:"
    print bstd
    #print blast_handle.
    #for r in blast_handle.getvalue():
    #    print r
    #for record in blast_output.getvalue().split("\n"):
    #    print record

    #if os.path.isfile( out_file):
    #    with open( out_file, "r") as tof:
    #        print tof.read()
    #else:
    #    print "No results"

    #child = Popen( str( blastp_cline),
    #               stdin=data, #PIPE,
    #               stdout=PIPE,
    #               stderr=PIPE,
    #               universal_newlines=True,
    #               shell=False)
#
#
    ##child.stdin.write("EIGEKGVNLSGGQK")
    #child.stdin.close()





    #main()

