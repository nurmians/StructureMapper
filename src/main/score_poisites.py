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
import getopt
import string
import math
import logging
import gc
import multiprocessing
import atexit

try:
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastxCommandline
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBList import PDBList
except ImportError
    stderr.write("Biopython installation not found. Try running: 'pip install biopython'.")  


from format_results import *
from collect_results import *
from write_results import *
from main_record import MainRecord

curdir = os.path.dirname(os.path.abspath(__file__))

sys.path.append( "%s/../utils" % curdir)
from struct_methods import *
from fasta_formatter import ConvertFastaLocatorFormat, ExtractCustomValuesFromFastaHeader
from fileops import *
import progress_bar as pb
import online_pdb

sys.path.append( "%s/../blast" % curdir)
from blast_structure_finder_w_align import FindStructures

sys.path.append( "%s/../disorder" % curdir)
from disorder_predictor import PredictDisorder, ReadDisorderResults

sys.path.append( "%s/../align" % curdir)
from find_residue_in_alignment import FindResidueInAlignment, FindResidueInPoiSiteAlignment
from align_seq_with_structure_multithread  import ProcessBlastResults
#from align_seq_with_structure  import ProcessBlastResults

sys.path.append( "%s/../asa" % curdir)
from pdbatoms import *
#from asa_uta import CalcASAForResidues #,runASAFiles

from config_writer import ReadConfig, WriteConfig

#valid_filename_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)

#Default folders for files and databases, change these.
DEFAULT_CONFIG_FILE = "./config.ini"
CONFIG_VALUES = {}
CUSTOM_COLS_FILE = "custom_cols.txt"
UNEXPECTED_EXIT = True

def SummarizeResults( aSeqFile, aBlastFile, aBitScoreThreshold, aOutputFolder, aGroupedResultsFile, aOutputFile, aDisorderFile, aObsoletes=[], aCustomCols=[] ):

    global CONFIG_VALUES

    summary_headers = ["ENTRY", "ENTRY_ID", "GENE", "POI", "POISEQ", "STRUCTURES", "ACCESSIBLE", "BLAST_ENTRIES", "BEST_BITSCORE", "BITSCORE_PASS", "BEST_ALIGNMENT", "ALIGN_PASS", "CONSIDERED_BLAST_ENTRIES","POI_FOUND_IN"] #, "DISORDER_SCORE"]

    #Add custom cols
    for custom_col in reversed( aCustomCols):
        if custom_col in summary_headers:
            LogIt( "SUMMARY ERROR: Cannot use custom column name '%s' because it is already in use." % custom_col, log_file, 2 )
            aCustomCols.remove( custom_col)
    summary_headers.extend( aCustomCols)


    cur_file = ""
    n_total = 0
    try:

        cur_file = aSeqFile

        with open( aSeqFile, "r") as sf:
            for line in sf:
                if line[ 0] == ">": n_total += 1

        seq_headers = FastaHeaderIter( aSeqFile)
        cur_file = aBlastFile
        blast_records = FastaIter( aBlastFile, aJoinSeqLines=False)
    except Exception as ex:
        sys.stderr.write( "SUMMARY ERROR: Error opening file: '%s'.\n" % (cur_file))
        sys.stderr.write( "               An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
        Exit( -1)

    first_line = True
    n_orig_records = 0
    n_struct_records = 0

    ProgressBar( aProgress=0, aTotal=n_total, aIntervalSeconds=1.0)

    #Open output
    with open( aOutputFile, "w") as opf:

        for seq_header in seq_headers:

            n_orig_records += 1
            summary_values = {}

            #EXTRACT DATA FROM SEQ FILE
            seq_acc = GetRecordACC( seq_header)

            summary_values["ENTRY"] = GetRecordACC( seq_header, aWithPos=False)
            summary_values["ENTRY_ID"] = GetRecordACC( seq_header, aWithPos=True)
            summary_values["GENE"] = GetGeneFromFastaHeader( seq_header)

            header_vals = ExtractCustomValuesFromFastaHeader( seq_header, aAllValues=True)

            #Collect values from fasta header
            for hv in header_vals.keys():
                if hv == "POI":      summary_values["POI"] = header_vals["POI"].split("-")[ 0] if "POI" in header_vals else ""
                elif hv == "POISEQ": summary_values["POISEQ"] = header_vals["POISEQ"] if "POISEQ" in header_vals else ""
                else: summary_values[ hv] = header_vals[ hv]

            #EXTRACT DATA FROM BLAST FILE
            blast_header, blast_entries = blast_records.next()
            blast_acc = GetRecordACC( blast_header)

            assert( seq_acc == blast_acc)
            #print "SEQ:", seq_acc , " BLAST:", blast_acc

            b_entries = ParseBlastEntries( blast_entries, aObsoletes=aObsoletes)
            summary_values['N_BLAST_ENTRIES'] = len( b_entries)
            summary_values['CONSIDERED_BLAST_ENTRIES'] = ",".join( ["%s|%s" % (x['pdb_id'],x['chain']) for x in b_entries])
            #POI_STRUCTURES

            #Get bitscore
            bitscores = []
            for b_entry in b_entries: bitscores.append( b_entry['bitscore'])
            summary_values["BEST_BITSCORE"] = max( bitscores) if len( bitscores) else -1

            if aBitScoreThreshold < 0.0001:
                summary_values["BITSCORE_PASS"] = "n/a"
            else:
                bitscore_pass = summary_values["BEST_BITSCORE"] >= aBitScoreThreshold
                summary_values["BITSCORE_PASS"] = "TRUE" if bitscore_pass else "FALSE"


            #EXTRACT DATA FROM RESULT PDB FILES
            poidir = InsertFolderSeparator( aOutputFolder) + "/results/poi/"

            search_acc = GetRecordACC( seq_header, aWithPos=True)
            search_str = SystemSpecificFolderSeparators( poidir + search_acc + "_*.poisite.pdb")

            pdbfiles = []
            rank = 0
            while True:
                rank += 1
                if rank > len( b_entries): break

                filename = search_str.replace( "*", "%i" % rank)
                if os.path.isfile( filename): pdbfiles.append( filename)

            #print "Looking for:", search_str,
            #pdbfiles = GetDirFiles( search_str)
            #print " Found:", len( pdbfiles)

            summary_values["STRUCTURES"] = len( pdbfiles)
            summary_values["ALIGN_PASS"] = "TRUE" if len( pdbfiles) > 0 else "FALSE"
            if len( pdbfiles) > 0: n_struct_records += 1

            align_scores = []
            accessibility = []
            poi_structures = []

            for file in pdbfiles:

                file_data = GetResultFileRemarkInfo( file)
                align_scores.append( int( file_data["ALIGN_SCORE"]))
                accessibility.append( file_data["ISOL_ASA"] )
                poi_structures.append( "%s|%s" % (file_data["PDBCODE"], file_data["POI_CHAIN"]) )

            summary_values["BEST_ALIGNMENT"] = max( align_scores) if len( align_scores) else "0"
            summary_values["ACCESSIBLE"] = "/".join( list( set( accessibility)))
            summary_values["POI_FOUND_IN"] = ",".join( poi_structures)

            #Write header
            if first_line:

                opf.write( "\t".join( summary_headers) + "\n")
                first_line = False

            #WRITE SUMMARY LINE
            for col in summary_headers:

                val = summary_values[ col] if col in summary_values else ""
                opf.write( "%s\t" % str( val))

            opf.write( "\n")

            ProgressBar( aProgress=n_orig_records, aTotal=n_total)



    ProgressBar( aProgress=n_orig_records, aTotal=n_total, aFinalize=True)

    LogIt( "SUMMARY INFO: Analysis was ran for %i sequences." % n_orig_records, aOutputFile, 1)
    LogIt( "SUMMARY INFO: %i (%.1f%%) of these sequences had structural information available.\n" % (n_struct_records, float(n_struct_records)/n_orig_records*100.0 ), aOutputFile, 1)
    print "SUMMARY INFO: Summary of results in '%s'" % aOutputFile

    #All Done
    return 0

#Remove blast entries that have less than "cutoff" bitscore values
def BitscoreFilter( main_records, cutoff=60.0, log_file=None):

    n_rejected = 0
    n_accepted = 0
    n_unknown  = 0

    if cutoff < 0.0001:
        LogIt( "BITSCORE INFO: No bitscore filter in use.", log_file)
        return

    LogIt( "BITSCORE INFO: Enforcing BLAST bitscore threshold %.1f..." % cutoff, log_file)
    info_str = ""

    #Get rid of structures with less than cutoff bitscore
    #Get rid of obsolete structures
    for m_rec in main_records:

        for b_entry in m_rec.BlastEntries()[:]: #Use a copy of the list to iterate through

            entry_id = m_rec.EntryIdentifier( b_entry["rank"])

            if type( b_entry["bitscore"]) != float:
                info_str = "BITSCORE ERROR: %s: '%s' bitscore unreadable." % (entry_id, b_entry["bitscore"])
                LogIt( info_str, log_file)
                n_unknown += 1
                m_rec.BlastEntries().remove( b_entry) #Remove
                continue

            if b_entry["bitscore"] < cutoff:
                info_str = "BITSCORE INFO: %s: '%.1f' is below accepted bitscore threshold." % (entry_id, b_entry["bitscore"])
                LogIt( info_str, log_file, -1)
                n_rejected += 1
                m_rec.BlastEntries().remove( b_entry) #Remove
                continue

            n_accepted += 1



    #Report stats
    n_total = n_accepted+n_rejected+n_unknown
    n_failed = n_rejected+n_unknown
    info_str  = "\n\n"
    #info_str += "%i structures were found obsolete.\n" % (n_obsolete)
    info_str += "BITSCORE INFO: %i structures did not have a readable bitscore.%s\n\n" % (n_unknown, ("     [OK]" if n_unknown == 0 else ""))
    info_str += "BITSCORE INFO: Bitscore filter results PASS: %i/%i (%.1f%%)\n" % (n_accepted, n_total, float( n_accepted) / n_total * 100.0)
    info_str += "BITSCORE INFO: Bitscore filter results FAIL: %i/%i (%.1f%%)\n\n" % (n_failed, n_total, float( n_failed) / n_total * 100.0)

    LogIt( info_str, log_file)

#Remove blast entries that have less than "cutoff" alignment scores
#Also fill blast entries with information about the alignment
#This filter must be run even with 0 cutoff
def AlignmentFilter( main_records, align_folder, cutoff=70, log_file=None, positions=0):

    n_site_found = 0
    n_site_total = 0

    if cutoff % 10 != 0:
        cutoff = cutoff - (cutoff % 10)
        cutoff = min( 100, cutoff)
        LogIt( "FILTER INFO: Setting ALIGNMENT cutoff to even %i. (%i/100 or above scores will PASS)" % (cutoff, cutoff), log_file)

    info_str = "FILTER INFO: Enforcing ALIGNMENT score threshold of %i/100..." % cutoff
    LogIt( info_str, log_file)

    align_folder = InsertFolderSeparator( align_folder)

    n_rejected = 0
    n_accepted = 0
    n_unknown  = 0
    n_alignment_files = 0

    n_cutsite_not_found = 0
    percentiles = [0]*11

    unknowns = []
    errorneus_pdbs = set()

    curnum = 0
    n_records = len( main_records)

    n_records_pass = 0
    n_records_fail = 0
    record_passed = False

    #Process alignment files
    for m_rec in main_records:

        curnum += 1
        record_passed = False

        if curnum % 2000 == 0:
            print "FILTER INFO: Processing record %i/%i (%.1f%%)" % ( curnum, n_records, float( curnum) / n_records* 100.0)

        for blast_entry in reversed( m_rec.BlastEntries()): #loop reversed to be able to remove

            n_site_total += 1
            entry_id = m_rec.EntryIdentifier( blast_entry["rank"])

            #DEBUG
            #if curnum > 70000:
            #    print "INFO: Processing record %s" % ( blast_entry["entry_id"])

            #print "ALiGNFILE : %s" % align_file
            align_file = ValidFilename( m_rec.AlignmentFile( blast_entry["rank"]))
            #print "ALiGNFILE2: %s" % align_file

            if not os.path.isfile( align_file):
                info_str = "FILTER WARNING: %s: Alignment file '%s' not found." % (entry_id, align_file)
                LogIt( info_str, log_file, 2)
                n_unknown += 1
                unknowns.append( align_file )
                errorneus_pdbs.add( blast_entry["pdb_filepath"])

                if len( errorneus_pdbs) > 5:
                    LogIt( "FILTER ERROR: Errors detected with alignments with more than 5 different pdb files. Exitting.", log_file, 2)
                    Exit( 0) #EXIT
                m_rec.BlastEntries().remove( blast_entry) #Remove
                continue #continue

            #Find POI site in structure
            #res, orig = FindResidueInAlignment( align_file, cut_sites[ i])
            #returns tuple of arrays: ( [pdbresnum,pdbres,pdbseq], [seqresnum,seqres,seqseq], align_score)
            res, orig, alignment_score = FindResidueInPoiSiteAlignment( align_file, 10, cutoff, log_file, aPositions=positions)
            n_alignment_files += 1


            if alignment_score > 100: #Maximum score is 100
                info_str = "FILTER WARNING: Alignment score over 100 in file %s. (%i)" % (align_file, alignment_score)
                LogIt( info_str, log_file, 2)
                alignment_score = 100

            if alignment_score < 0: #No alignment score was found

                n_rejected += 1
                n_cutsite_not_found += 1

                LogIt( "FILTER INFO: No aligned residues at the POI-site were found in structure %s for alignment '%s'." % (blast_entry["struct_id"], os.path.basename( align_file)), log_file, 1)

                m_rec.BlastEntries().remove( blast_entry) #Remove

            elif alignment_score < cutoff: #Below threshold

                n_rejected += 1

                percentiles[ (alignment_score / 10)] += 1

                info_str = "FAIL: %s: Alignment with pdb %s scored %i." % (entry_id, blast_entry["struct_id"] , alignment_score)
                LogIt( info_str, log_file, -1)

                m_rec.BlastEntries().remove( blast_entry) #Remove


            else:
                n_site_found += 1
                n_accepted += 1
                record_passed = True

                percentiles[ (alignment_score / 10)] += 1

                info_str = ""
                info_str += "PASS:%s: score %i " % (entry_id, alignment_score)
                info_str += "Seq %i aa '%s', in \"%s\".\n" % (tuple( orig))
                info_str += "%s Matching PDB '%s' " % (((6+len(entry_id))*" "), blast_entry["struct_id"])
                info_str += "residue is number %i, '%s' in seq \"%s\"." % ( tuple( res))

                LogIt( info_str, log_file, -1)

                blast_entry["alignment_score"] = alignment_score
                blast_entry["pdb_poipos"] = res[ 0]
                blast_entry["pdb_poires"] = res[ 1]
                blast_entry["pdb_poiseq"] = res[ 2]
                blast_entry["seq_poipos"] = orig[ 0]
                blast_entry["seq_poires"] = orig[ 1]
                blast_entry["seq_poiseq"] = orig[ 2]

        if record_passed: n_records_pass += 1
        else: n_records_fail += 1

    #Report stats
    n_total = n_accepted+n_rejected+n_unknown
    if n_total == 0:
        LogIt( "ALIGN ERROR: No alignments were located in structures. Exitting.", log_file, 2)
        return -1

    n_failed = n_rejected+n_unknown
    info_str  = "\n"
    info_str += "ALIGN INFO: Alignment files processed: %i/%i \n\n" % (n_alignment_files,n_total)
    LogIt( info_str, log_file)

    if n_unknown > 0:

        LogIt( "ALIGN WARNING: %i alignment%s did not have a readable score.\n" % (n_unknown, "s" if n_unknown > 1 else ""), log_file)
        LogIt( "\n".join( unknowns[:20]), log_file, 1)
        if n_unknown > 20: LogIt( "Only first 20 of %i logged." % n_unknown, log_file, 1)

    assert( n_records_pass+n_records_fail == len( main_records))

    info_str  = ""
    info_str += "Alignment score filter records PASS: %i/%i (%.1f%%)\n" % (n_records_pass, len( main_records), float( n_records_pass) / len( main_records) * 100.0)
    info_str += "Alignment score filter records FAIL: %i/%i (%.1f%%)\n\n" % (n_records_fail, len( main_records), float( n_records_fail) / len( main_records) * 100.0)

    info_str += "Alignment score filter blast entries PASS: %i/%i (%.1f%%)\n" % (n_accepted, n_total, float( n_accepted) / n_total * 100.0)
    info_str += "Alignment score filter blast entries FAIL: %i/%i (%.1f%%)\n\n" % (n_failed, n_total, float( n_failed) / n_total * 100.0)

    LogIt( info_str, log_file)

    info_str = "Distribution:\n"

    for p in range( len(percentiles)-1, -1, -1):
        if p*10 < cutoff and p*10+10 >= cutoff: info_str += "-----\n" #Print threshold
        info_str += "%3i: %4i (%.1f%%)\n" % (p*10,percentiles[p],(percentiles[p]/float((n_total-n_unknown)-n_cutsite_not_found))*100.0)


    info_str += "\n"
    info_str += "POI not found:  %i\n" % n_cutsite_not_found
    info_str += "File not found: %i\n" % n_unknown
    info_str += "Total files:    %i\n" % n_total
    LogIt( info_str, log_file)

    #print percentiles

    print ""
    print "Total of %i/%i POI sites found in structures by sequence alignment." % (n_site_found, n_site_total)

    return 0


#Processes POIs (points of interest) for sequences in seq_file
#blast needs to be done and results have to be in blast_folder (*.blast)
#alignment needs to be done and results need to be in align_folder
def ProcessPoiSites( seq_file, blast_folder, align_folder, output_folder, pdb_folder, precalc_asa_folder=None, dssp_folder="", obsoletes=[], positions=0, reset=False, \
                     bit_score_cutoff=60.0, align_score_cutoff=70, threads=-1, aCustomColsFile="", sidechain_only=False, aWriteUsedBiomolecules=False, aDSSPExecutable="" ):


    print "INFO: Processing POIs in file: '%s'..." % seq_file

    #Read sequence records
    #cut_records = list( SeqIO.parse( seq_file, "fasta"))
    #cut_record_ids = GetRecordIDs( cut_records )
    align_folder = InsertFolderSeparator( align_folder)
    blast_folder = InsertFolderSeparator( blast_folder)


    seq_record_index = SeqIO.index( seq_file, "fasta", key_function=GetRecordACC)
    n_seq_records = len( seq_record_index)


    if not n_seq_records:
        sys.stderr.write( "ERROR: No sequence record entries found in file '%s'.\n" % seq_file)
        return -1


    #Read blast results
    blast_file = ChangeToFolder( seq_file+".blast", blast_folder)

    if not os.path.isfile( blast_file):
        sys.stderr.write("ERROR: Blast result file '%s' not found.\n" % (blast_file))
        return -1 #RETURN
    else:
        print "INFO: Processing Blast file: '%s'" % blast_file


    #Construct alignment file path
    seq_path, seq_filename = os.path.split( seq_file)
    seq_name, seq_suffix = SplitFileName( seq_filename)


    try:
        cur_file = seq_file
        seq_headers = FastaHeaderIter( seq_file)
        cur_file = blast_file
        blast_records = FastaIter( blast_file, aJoinSeqLines=False)
    except Exception as ex:
        sys.stderr.write( "PP ERROR: Error opening file: '%s'.\n" % (cur_file))
        sys.stderr.write( "          An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
        Exit( -1)


    #Get all acc codes in same order as they are in the file
    main_records = []
    MainRecord.SetAlignPathAndPrefix( "%s%s_" % (align_folder, seq_name))

    #Create main_records
    for sh in seq_headers:

        acc = GetRecordACC( sh, aWithPos=True)
        poi = GetRecordPOI( sh)

        blast_header, blast_entries = blast_records.next()
        blast_acc = GetRecordACC( blast_header)

        assert( acc == blast_acc) #Sanity check

        mr = MainRecord( acc, poi, seq_record_index, ParseBlastEntries( blast_entries, pdb_folder, obsoletes))
        main_records.append( mr)

    #Fasta records need to contain 2 fields:
    #|POI:2053-2054|POISEQ:EKSA-ATWD
    #Removes entries with invalid fasta header formatting from main_records
    #n_cutsites = CheckRecords( main_records)

    #Sanity check
    if n_seq_records != len( main_records):
        sys.stderr.write("ERROR: The number of POIs does not match the number of records (%i vs. %i).\n" % (len( cut_sites), len( main_records)))
        return -3 #RETURN

    #BLAST result structure:
    #         0                 1                       2           3                 4           5          6         7       8         9       10      11
    # Fields: query id,         subject id,             % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    #sp|P31946-2|1433B_HUMAN    gi|67464627|pdb|2BQ0|A  99.58       238 1   0   1   238 3   240 2e-176   489


    if precalc_asa_folder != None: InsertFolderSeparator( precalc_asa_folder)
    output_folder = InsertFolderSeparator( output_folder)
    log_filename = output_folder + "poi_filter_log.txt"
    #ClearLog( log_filename) #Start new log_file
    print "Generating log file: '%s'" % log_filename

    retval = -1

    #Clears and opens a log file
    with open( log_filename, 'w') as poisite_log:

        #Filter by bitscore
        if bit_score_cutoff > 0.00001:
            BitscoreFilter( main_records, bit_score_cutoff, poisite_log) #Removes filtered entries from blast index
        else:
            LogIt( "\nINFO: Bitscore filtering disabled.", poisite_log)

        #Filter by alignment score
        retval = AlignmentFilter( main_records, align_folder, align_score_cutoff, poisite_log, positions)

        if len( main_records) == 0 or retval < 0:
            LogIt( "\nINFO: No cutsites were found in PDB structure data.", poisite_log)
            return -7

        #Write result files
        try:
            retval = WriteResultFiles( main_records, seq_file, positions, sidechain_only, output_folder, dssp_folder, threads, reset, poisite_log,
                                       precalc_asa_folder, aCustomColsFile=aCustomColsFile, aWriteUsedBiomolecules=aWriteUsedBiomolecules, aDSSPExecutable=aDSSPExecutable)
        except Exception as ex:
            sys.stderr.write( "MAIN: Error occurred during writing of result files. (%s)\n" % str( ex))
            raise
            retval = -5 #Error during WriteResultFiles

    return retval

#NOT USED
def ProcessASAQueue( asa_queue, sidechain_only, threads=-1 ):

    if len( asa_queue) == 0:
        print ""
        print "ASA already calculated for all files..."
        return

    print ""
    print "Calculating ASA for %i files..." % len( asa_queue)
    print ""

    #List indices
    PDB = 0 #input
    ASA = 1 #output
    CHAIN = 2
    #REMARKS = 3
    #PRECALC_FILE = 4

    input_file_list = []
    chain_list = []
    outfiles_list = []
    all_outfiles = []
    all_infiles = []

    uniques = {}
    outputs = {}

    for qi in reversed( range( len( asa_queue))):

        identifier = asa_queue[ qi][ PDB] + "|" + asa_queue[ qi][ CHAIN]
        output = asa_queue[ qi][ ASA]

        if output in outputs:
            sys.stderr.write( "ASA QUEUE ERROR: Output file '%s' multiple times in queue. Removed: '%s'\n" % (output, identifier.replace("\t", ",")))
            del asa_queue[ qi]
            continue

        outputs[ output] = True
        file_list.append( qitem[ PDB])
        outfiles_list.append( qitem[ ASA])
        chain_list.append( qitem[ CHAIN])


    #Full residues or sidechains only
    calc_mode = "R" if not sidechain_only else "S"

    #Run with default modification mode "convert and remove PTMs"
    processed_files = runASAFiles( file_list, chain_list, outfiles_list, threads=threads, skip_existing_files=False, firstmodelonly=True, calc_mode=calc_mode)


    print ""
    print "ASA calculations Done."
    print ""


#For reliability column
#Returns value 0.0 - 1.0, 1.0 being a good match
def CheckWordMatch( seq_desc, pdb_desc):

    seq_desc = seq_desc.lower()
    pdb_desc = pdb_desc.lower()

    #Remove special chars form comparison
    seq_desc = re.sub('[-/]', ' ', seq_desc)
    seq_desc = re.sub('[\(\)]', '', seq_desc)
    seq_desc = seq_desc.split(" ")

    pdb_desc = re.sub('[-/]', ' ', pdb_desc)
    pdb_desc = re.sub('[\(\)]', '', pdb_desc)
    pdb_desc = pdb_desc.split(" ")

    retval = 0.0
    matches = ""

    exclude = ["protein", "proteins", "subunit", "fragment", "isoform", "molecule", "chain", "domain", "family", "mutant", "factor", \
               "and", "of", "the", "a", "with", "in", "wild", "type", \
               "human", "mouse" ,"rabbit", "bovine", "cow", "yeast", "drosophila", "melanogaster", "xenopus", "laevis", "rattus", "norvegicus", "gallus" ]

    #Check for acronyms
    acronyms_seq = []
    letters_seq = []
    for word in seq_desc:
        if not len( word): continue
        letters_seq.append( word[ 0])
        m = re.match("([a-z]{3,})(\d+)", word)
        if m: acronyms_seq.append( m.group(1))

    acronyms_pdb = []
    letters_pdb = []
    for word in pdb_desc:
        if not len( word): continue
        letters_pdb.append( word[ 0])
        m = re.match("([a-z]{3,})(\d+)", word)
        if m: acronyms_pdb.append( m.group(1))


    for acronym in acronyms_seq:
        for l in range( len(letters_pdb)-len( acronym)):
            if acronym == letters_pdb[l:len( acronym)]:
                retval = 0.8
                if len( matches): matches += "," + acronym
                else: matches = acronym
                break

    for acronym in acronyms_pdb:
        for l in range( len(letters_seq)-len( acronym)):
            if acronym == letters_seq[l:len( acronym)]:
                retval = 0.8
                if len( matches): matches += "," + acronym
                else: matches = acronym
                break


    #Look for similar words in both desc and title
    for word in seq_desc:
        if not len( word): continue
        if word in exclude: continue

        #if len( word) <= 2 and word.islower(): continue
        #if word.isdigit() and float( word) < 5.0: continue

        if word in pdb_desc:

            if word.isdigit():
                if float( word) < 5.0:
                    retval += 0.1
                elif float( word) < 10.1:
                    retval += 0.4
                elif float( word) < 20.1:
                    retval += 0.6
                else:
                    retval += 1.0
            elif len( word) < 2:
                retval += 0.5
            else:
                retval = 1.0

            if len( matches): matches += "," + word
            else: matches = word

            if retval > 0.99: break

    #enforce limits
    retval = min( 1.0, retval)
    retval = max( 0.0, retval)

    return (retval, matches)


#Create link
#Coldict is for result_file_entry only
#Pymol link column works on Windows if WinPymol.exe is in system path (see/edit file pymol_launcher.bat if you want to launch without path)
#This method fills PYMOL and PYMOL_BIOMOL columns in the grouped_results_file.txt and creates clickable hyperlink for viewing result structures using Excel and Pymol
#For each result file row, a batch (.bat) file is created. When this batch file is launched, it reads a row from pymol_script_lib file that contains instructions
#for PyMol to load and show the resulting structures with the POI highlighted in red. The scripts are stored in a single file to avoid creating thousands of files that
#are slow to delete with larger runs. When the hyperlink is clicked, the one row of instructions is written to a temporary file that is read into Pymol
#when it is launched. These extra steps are necessary, because the length of the instructions can exceed what can be inputted directly from the commandline. Also, the maximum number of
#structures that are loaded is restricted to 20 structures, due to pymol load time and batch file variables not being able to handle longer strings.
def SetPymolLinkCol( g_entry, g_entry_map, grouped_result_file_row, grouped_row_map):

    global CONFIG_VALUES

    pymol_exe = ""
    grouped_result_file_row[ grouped_row_map['PYMOL']] = ""
    grouped_result_file_row[ grouped_row_map['PYMOL_BIOMOL']] = ""

    #Tested only on Windows
    if os.name != 'nt': return

    pymol_exe = CONFIG_VALUES["PYMOL_PATH"]

    if len( pymol_exe) == 0 or len( g_entry[ 1]) == 0: return

    #first_seq_id = g_entry[ grouped_row_map['SEQ_ID']][ 0]
    poifiles = []
    biomolfiles = []
    #print "\nGENTRY:"
    #print g_entry
    #print "\GENTRY COLDICT:"
    #print g_entry_map
    #print "\nRESFILE:"
    #print grouped_result_file_row
    #print "\nGROUPED COLDICT:"
    #print grouped_row_map
    #print "\nCHAINS:"
    #print g_entry[ g_entry_map['POI_CHAIN']]
    poi_chains = g_entry[ g_entry_map['POI_CHAIN']][0:20] #View only first 20

    for seq_id in g_entry[ g_entry_map['ENTRY']][0:20]:
        poifile = SystemSpecificFolderSeparators( InsertFolderSeparator(CONFIG_VALUES['DEFAULT_RESULT_DIR']) + "results/poi/" + seq_id + ".poisite.pdb")
        biomolfile = SystemSpecificFolderSeparators( InsertFolderSeparator(CONFIG_VALUES['DEFAULT_RESULT_DIR']) + "results/biomol/" + seq_id + ".biomol.pdb")

        if not os.path.isfile( poifile):
            sys.stderr.write( "PYMOL WARNING: POI-pdb file '%s' not found for entry '%s'.\n" % (poifile, seq_id))
            return

        poifiles.append( poifile)

        if os.path.isfile( biomolfile):
            biomolfiles.append( biomolfile)


    #poifile_for_batchfile = os.path.abspath( poifile)
    #biofile_for_batchfile = os.path.abspath( biomolfile)

    pymol_exe = SystemSpecificFolderSeparators( pymol_exe)
    view_dir  = SystemSpecificFolderSeparators( InsertFolderSeparator(CONFIG_VALUES['DEFAULT_RESULT_DIR']) + "results/view/")
    biomol_view_dir  = SystemSpecificFolderSeparators( InsertFolderSeparator(CONFIG_VALUES['DEFAULT_RESULT_DIR']) + "results/view_biomol/")

    scriptlib = view_dir + "script_library.txt"
    biomol_scriptlib = biomol_view_dir + "script_library.txt"
    script_file = view_dir + "launch_pymol.txt"
    biomol_script_file = biomol_view_dir + "launch_pymol.txt"
    #print grouped_result_file_row
    #print "DICT:", grouped_row_map
    #print "INDEX", grouped_row_map['SEQ_ID']
    #print "SEQID:", grouped_result_file_row[ grouped_row_map['SEQ_ID']]
    #sys.exit( 0)
    entry_id = grouped_result_file_row[ grouped_row_map['SEQ_ID']]
    batchFile = SystemSpecificFolderSeparators( view_dir + "view_" + entry_id + ".bat")
    bio_batchFile = SystemSpecificFolderSeparators( biomol_view_dir + "view_biomol_" + entry_id + ".bat")

    if not hasattr( SetPymolLinkCol, "batch_written"):
        #First time
        if not os.path.exists( view_dir): os.makedirs( view_dir)
        if not os.path.exists( biomol_view_dir): os.makedirs( biomol_view_dir)

        try: os.unlink( scriptlib)
        except: pass

        try: os.unlink( biomol_scriptlib)
        except: pass

        SetPymolLinkCol.batch_written = True

        import shutil
        #Copy launcher files
        try:
            utils_dir = "%s/../utils/" % os.path.dirname(os.path.abspath(__file__))
            launcher_file = utils_dir + "pymol_launcher.bat"
            if os.path.isfile( launcher_file):
                #if True or not os.path.isfile( view_dir + "pymol_launcher.bat"):        shutil.copyfile( utils_dir + "pymol_launcher.bat", view_dir + "pymol_launcher.bat")
                #if True or not os.path.isfile( biomol_view_dir + "pymol_launcher.bat"): shutil.copyfile( utils_dir + "pymol_launcher.bat", biomol_view_dir + "pymol_launcher.bat")
                SetPymolLinkCol.launcher = launcher_file
            else:
                sys.stderr.write("PYMOL ERROR: No Pymol launcher file found.\n")
                return
        except Exception as ex:
                sys.stderr.write("PYMOL ERROR: Failed to set up PYMOL launcher file.\n")
                sys.stderr.write( str( ex) + "\n")
                return


    elif not hasattr( SetPymolLinkCol, "launcher"):
        #Launcher bat file has not been found
        return


    pdbcodes = grouped_result_file_row[ grouped_row_map['PDB_FILES']].split(",")

    try:
        #Write batch file for viewing structure with relevant chains
        WritePymolScriptLib( scriptlib, entry_id, poifiles, poi_chains, aPDB_Codes=pdbcodes)
        with open( batchFile, "w+") as bf: bf.write( "%s %s %s %s\n" % ( os.path.abspath( SetPymolLinkCol.launcher), entry_id, os.path.abspath( scriptlib), os.path.abspath(script_file)) )
        grouped_result_file_row[ grouped_row_map['PYMOL']] = '=HYPERLINK("file://%s","view")' % ( os.path.abspath( batchFile).replace("\\","/"))

        #Write batch file for biomol
        if len( biomolfiles):
            WritePymolScriptLib( biomol_scriptlib, entry_id, biomolfiles, poi_chains, aPDB_Codes=pdbcodes)
            with open( bio_batchFile, "w+") as bmbf: bmbf.write( "%s %s %s %s\n" % ( os.path.abspath( SetPymolLinkCol.launcher), entry_id, os.path.abspath( biomol_scriptlib), os.path.abspath(biomol_script_file)) )
            grouped_result_file_row[ grouped_row_map['PYMOL_BIOMOL']] = '=HYPERLINK("file://%s","view")' % ( os.path.abspath( bio_batchFile).replace("\\","/"))

    except KeyboardInterrupt:
        print "Process interrupted by user. Exitting..."
        Exit( 0) #Exit
    except Exception:
        sys.stderr.write( "PYMOL VIEW ERROR: Unexpected error: %s\n" % sys.exc_info()[0])
        raise

    #DEBUG
    #sys.exit( 0)


#Write Pymol startup scripts into a single file to avoid creating thousands of extra files
#Each excel row in grouped results file needs to have own launcher .bat file, because Excel cannot pass arguments to files without VB
def WritePymolScriptLib( aLibFile, aEntryId, aPoiFiles, aPoiChains, aPDB_Codes=[]):

    #launch_cmd = '"%s" -d "' % aPymolExe
    launch_cmd = "%s " % aEntryId

    first_obj = None
    first_chain = None

    #print "POI_CHAINS:", aPoiChains

    for i in range( len( aPoiFiles)):

        poifile = os.path.abspath( aPoiFiles[ i])
        pdb_code = aPDB_Codes[ i]

        orig_sel = ""
        obj_name = aPoiFiles[ i].replace( ".poisite.pdb", "")
        obj_chain = aPoiChains[ i]

        if len( pdb_code) == 4:
            #orig_sel = "select %i_%s, poi;" % (i, pdb_code)
            obj_name = "%s" % (pdb_code.upper())

        launch_cmd += 'load %s, %s;select poi_%s, %s and b > 94.0;show spheres, poi_%s and name ca;%s' % (poifile, obj_name, obj_name, obj_name, obj_name, orig_sel)

        #Align structures
        if first_obj != None:
            launch_cmd += 'align %s and chain %s, %s and chain %s;' % (obj_name, obj_chain, first_obj, first_chain)
        else:
            first_obj = obj_name
            first_chain = obj_chain

    #Finally:
    launch_cmd += 'orient %s;hide lines;show cartoon;spectrum b;' % first_obj

    #Write one row to script lib file
    with open( aLibFile, "a+") as libf:
        libf.write( launch_cmd + "\n" )


def SetReliability( result_file_entry, col_dict):

    r_entry = result_file_entry
    r_entry[ col_dict['RELIABILITY']] = 100.0
    return

    #similar words in title and desc
    #Align score
    #bitscore
    #Organism
    #resolution

    #"RELIABILITY", "RE_BITSCORE", "RE_IDENTITY", "RE_DESC", "RE_RESOLUTION", "RE_ALIGNMENT"

    r_entry = result_file_entry

    #Matching words
    r_entry['RE_DESC'], r_entry['DEBUG_DESC'] = CheckWordMatch( r_entry["DESC"], r_entry["TITLE"] + " " + r_entry["COMPOUND"])


    align_pts = float(r_entry['ALIGN_SCORE']) / 100.0
    align_pts = min( 1.0, align_pts)
    r_entry['RE_ALIGNMENT'] = align_pts

    r_entry['RE_IDENTITY'] = 0.0
    if r_entry['ORGANISM_TAXID'].find("9606") > 0 or float(r_entry['IDENTITY']) > 79.99999:
        r_entry['RE_IDENTITY'] = 1.0

    bitscore_pts = float(r_entry['BITSCORE']) / 500.0
    bitscore_pts = min( 1.0, bitscore_pts)
    r_entry['RE_BITSCORE'] = bitscore_pts

    resolution_pts = 0.0
    resolution = float(r_entry['RESOLUTION'])
    if resolution < 0.0001:
        #NMR and solution scattering files do not have resolution reported
        resolution_pts = 0.5
    else:
        #points [x4] = y0
        #       [x2] = y1
        # 0 = a4 + c
        # 1 = a2 + c => c = -a4a => 1 = a2 - a4 => a = -0,5 => c = -2
        resolution_pts = -0.5*resolution + 2

    resolution_pts = min( 1.0, resolution_pts)
    resolution_pts = max( 0.0, resolution_pts)
    r_entry['RE_RESOLUTION'] = resolution_pts

    mult = [ 0.2, 0.2, 0.2, 0.1, 0.3 ]
    retval = r_entry['RE_DESC']       * mult[0] + \
             r_entry['RE_ALIGNMENT']  * mult[1] + \
             r_entry['RE_BITSCORE']   * mult[2] + \
             r_entry['RE_IDENTITY']   * mult[3] + \
             r_entry['RE_RESOLUTION'] * mult[4]


    #If there is no DSSP assesment available, there's something wrong with the structure
    if r_entry['DSSP'].find("Unknown") > 0: retval = retval / 2.0

    r_entry['RELIABILITY'] = retval


def FormatResults( input_file, output_file, disorder_file=None, aCustomColumns=[]):

    from collections import OrderedDict

    content = []

    try:
        with open( input_file, 'r') as f:
            content = f.readlines()
    except IOError as e:
        sys.stderr.write( "FORMAT ERROR: File '%s'\n" % input_file )
        sys.stderr.write( "FORMAT ERROR: File error({0}): {1}\n".format(e.errno, e.strerror))
        return -1
    except:
        sys.stderr.write( "FORMAT ERROR: Unexpected error: %s\n" % sys.exc_info()[0])
        return -2

    line = 0
    col_map = {} #Map from collected results to grouped results
    grouped_entries = []
    cur_id = ""
    prev_id = "?"
    cur_entry = []
    dataseries =  GetFilename( output_file).replace( "_grouped_results.txt", "")

    #Group entries together based on their seq_id
    for entry in content:

        line += 1
        cols = entry.split("\t")

        #First line contains column headers create a dict for mapping
        if line == 1:
            i = 0
            for col in cols:
                #Remove last trailing newline
                col_map[ col.strip()] = i
                i += 1
            #print col_map
            #Exit(0)
            continue

        cur_protein = cols[ 0]
        cur_id = re.sub(r'_\d+$', '', cols[ 1]) #Remove blast rank information from SEQ_ID to group entries "_1"

        if cur_id != prev_id or cur_protein != prev_protein:
            #New entry found
            #print "new entry: %s" % cur_id
            #Append previous to results (if it exists)
            if len( cur_entry):
                grouped_entries.append( cur_entry[:])

            cur_entry = []
            for col in cols:
                cur_entry.append( [col.rstrip()])

        else:
            #Add to previous entry
            #print "existing entry: %s" % cur_id
            for c in range( len(cols)):
                cur_entry[ c].append( cols[ c].rstrip())

        prev_protein = cur_protein
        prev_id = cur_id

    #Append last entry
    if len( cur_entry):
        grouped_entries.append( cur_entry[:])


    #Results are now formatted to a list grouped by SEQ_ID

    #Colum types: any, avg, min, max, 1st, list,special, count, delta, check_identity
    #FORMATTED (GROUPED RESULT FILE COLUMNS)
    #How to map result file columns in grouped file
    col_types = [ ("DATASERIES","check_identity",""), \
                  ("SEQ_ID","check_identity",""),\

                  #("RELIABILITY","avg","RELIABILITY"),\

                  ("GENE","any",""),\
                  ("DESC","set",""),\

                  ("N_STRUCTURES","count","SEQ_ID"),\

                  ("PYMOL","special",""),\
                  ("PYMOL_BIOMOL","special",""),\
                  ("ASA_AVG","special",""),\
                  ("TEMPF_SCORE","special",""),\
                  ("IUPRED","special",""),\
                  #("DSSP_SCORE","special",""),\
                  ("DSSP","list",""),\

                  ("INTERFACE","vote",""),\

                  #negative value means missing value
                  ("ISOL_ASA_DECISION","vote","ISOL_ASA"),\
                  ("ISOL_ASA_AVG","avg_no_negatives","ISOL_ASA_AVG"),\
                  ("BIOL_ASA_DECISION","vote","BIOL_ASA"),\
                  ("BIOL_ASA_AVG","avg_no_negatives","BIOL_ASA_AVG"), \
                  ("ASYM_ASA_DECISION","vote","ASYM_ASA"),\
                  ("ASYM_ASA_AVG","avg_no_negatives","ASYM_ASA_AVG"),\

                  ("BIOL_ASA_POI_MIN","list",""),\
                  ("BIOL_ASA_POI_MAX","list",""),\
                  ("BIOL_ASA_POI_DELTA","delta_no_negatives","BIOL_ASA_POI_MIN,BIOL_ASA_POI_MAX"),\

                  ("ISOL_ASA_POI_MIN","list",""),\
                  ("ISOL_ASA_POI_MAX","list",""),\
                  ("ISOL_ASA_POI_DELTA","delta_no_negatives","ISOL_ASA_POI_MIN,ISOL_ASA_POI_MAX"),\

                  ("ASYM_ASA_POI_MIN","list",""),\
                  ("ASYM_ASA_POI_MAX","list",""),\
                  ("ASYM_ASA_POI_DELTA","delta_no_negatives","ASYM_ASA_POI_MIN,ASYM_ASA_POI_MAX"),\

                  ("TEMPF_POI_AVG","avg_nozeroes","TEMPF_AVG"),\
                  ("TEMPF_POI_VALS","list","TEMPF_AVG"),\
                  ("TEMPF_POI_MIN","list","TEMPF_POI_MIN"),\
                  ("TEMPF_POI_MAX","list","TEMPF_POI_MAX"),\
                  ("TEMPF_FILE_MIN","list",""),\
                  ("TEMPF_FILE_MAX","list",""),\

                  ("PDB_FILES","list","PDBCODE"),\
                  ("PDB_CHAINS","list","POI_CHAIN"),\
                  ("PDB_SITES","list","POIRESIDUE_1"),\
                  ("PDB_PTMS","list","PDBPTM"),\

                  ("BLAST_RANKS","list","RANK"),\

                  ("ENTRY","matching","ENTRY"),\
                  ("SEQ_POIPOS","any","POIPOS"),\
                  ("SEQ_POISEQ","any","POISEQ"),\
                  ("PDBSEQ","list",""),\
                  ("ALIGNMENTS","list","ALIGN_SCORE"),\
                  ("ALIGN_AVG","avg","ALIGN_SCORE"),\
                  ("BITSCORE_AVG","avg","BITSCORE"),\
                  ("IDENTITY_AVG","avg","IDENTITY"),\
                  ("METHODS","list","METHOD"),\
                  #("RESOLUTION_AVG","avg_nozeroes","RESOLUTION"),\
                  ("R-FREE_AVG","avg_nozeroes","R-FREE"),\
                  ("RESOLUTIONS","list","RESOLUTION"),\
                  ("ORGANISM_SCIENTIFIC","list",""),\
                  ("ORGANISM_TAXID","list",""), \

                  ("INTERFACES","list","INTERFACE"),\
                  ("INTERFACE_ASA_DELTA","list","INTERFACE_ASA_DELTA"),\
                  ("INTERFACE_ASA_DELTA_AVG","avg","INTERFACE_ASA_DELTA"),\

                  #("CONSIDERED_HOMOLOGS","special",""),\

                  ("BIOL_ASA_DECISIONS","list","BIOL_ASA"),\
                  ("BIOL_ASA_VALS","list","BIOL_ASA_AVG"),\

                  ("ISOL_ASA_DECISIONS","list","ISOL_ASA"),\
                  ("ISOL_ASA_VALS","list","ISOL_ASA_AVG"),\

                  ("ASYM_ASA_DECISIONS","list","ASYM_ASA"),\
                  ("ASYM_ASA_VALS","list","ASYM_ASA_AVG")\

                ]

    #Custom columns
    colnames = zip(*col_types)[ 0] #first elem of tuples in list
    aCustomColumns = list( set( aCustomColumns)) #Make sure that values are unique
    for cfhv in reversed( aCustomColumns):
        if cfhv in colnames:
            sys.stderr.write( "FORMAT ERROR: Cannot use custom column name '%s' because it is already in use.\n" % cfhv)
            aCustomColumns.remove( cfhv)
        else:
            print "INFO: Adding custom column '%s' to grouped results." % cfhv
            col_types.append( (cfhv,"check_identity","")) #Custom cols should have same value between blast entries

    #Create dictinary for mapping cols
    col_dict = {}
    for i in range( len( col_types)):
        col_dict[ col_types[ i][ 0]] = i

    print "INFO: Number of columns: %i" % len( col_types)

    disorder_result = {}
    if disorder_file:

        if os.path.isfile( disorder_file):
            try:
                disorder_result = ReadDisorderResults( disorder_file)
                print "INFO: Read %i disorder results." % len( disorder_result)
            except:
                sys.stderr.write( "ERROR: Could not read disorder results in file ('%s')\n" % disorder_file)
        else:
            sys.stderr.write( "ERROR: Disorder file specified but not found. ('%s')\n" % disorder_file)
            sys.exit( -1)
    else:
        print "INFO: Disorder calculation file not specified."


    progress = 0
    total = len( grouped_entries)
    print "INFO: Processing grouped entries..."
    pb.Init( progress, total)

    formatted_entries = []
    for g_entry in grouped_entries:

        #g_entry contains rows from individual results file
        #new_entry is those rows combined to one row

        new_entry = ProcessGroupedEntry( g_entry, col_types, col_map)


        #Set disorder prediction
        if len( disorder_result):

            #print "SEQ_ID:", new_entry[ col_dict['SEQ_ID']]
            entry_id = new_entry[ col_dict['SEQ_ID']] #SEQ_ID includes Pos (if it exists) but not blast rank.
            #regex = re.search( "^(.+)_\d+\s*$", g_entry[ col_dict['ENTRY']])
            #
            #if regex:
            #    entry = regex.group( 1)
            #    #print "ENTRY:", entry
            if entry_id in disorder_result:
                new_entry[ col_dict['IUPRED']] = str( disorder_result[ entry_id])
            else:
                new_entry[ col_dict['IUPRED']] = "-1.0"
                pb.PrintAbove( "WARNING: No disorder prediction found for '%s'.\n" % ( entry_id))
                #sys.stderr.write( "WARNING: No disorder prediction found for '%s'.\n" % ( entry_id))

        new_entry[ col_dict['DATASERIES']] = dataseries

        #Special
        SetPymolLinkCol( g_entry, col_map, new_entry, col_dict)
        SetAsaAvg( new_entry, col_dict)
        SetTempfScore( new_entry, col_dict)


        #SetHomologs( new_entry, col_dict)

        #SetReliability( new_entry, col_dict)
        #SetScore( new_entry, col_dict, disorder_result)
        #SetASADecision( new_entry, col_dict)

        #Append
        formatted_entries.append( new_entry)
        progress += 1
        pb.Update( progress, total, aShowSpeed=True, aEstimateCompletion=True)


    pb.Finalize( progress, total)


    n_entries = len( formatted_entries)
    print "Writing %i entries to file: '%s'..." % (n_entries,output_file)

    try:
        with open( output_file, 'w') as f:
            #Write column headers
            headers = "\t".join( [i[0] for i in col_types])
            f.write( headers + "\n")
            #Write data
            for entry in formatted_entries:
                #print entry[ 0] + "\n"

                f.write( "\t".join( [str(i) for i in entry]) + "\n")

    except IOError as e:
        sys.stderr.write( "ERROR: Formatted file error({0}): {1}\n".format(e.errno, e.strerror))
        return -1
    except:
        sys.stderr.write( "ERROR: Unexpected error when writing formatted results file: %s\n" % sys.exc_info()[0])
        raise
        #return -2

    print "INFO: Wrote %i formatted entries." % n_entries
    #All done formatting
    return n_entries

#H = α-helix
#B = residue in isolated β-bridge
#E = extended strand, participates in β ladder
#G = 3-helix (310 helix)
#I = 5 helix (π-helix)
#T = hydrogen bonded turn
#S = bend
#_ = loop or irregular
#? = missing from structure
def ScoreDSSP( dssp_output ):

    if dssp_output.find("Unknown") >= 0:
        return 0.0

    n_H = dssp_output.count('H')
    n_B = dssp_output.count('B')
    n_E = dssp_output.count('E')
    n_G = dssp_output.count('G')
    n_I = dssp_output.count('I')
    n_T = dssp_output.count('T')
    n_S = dssp_output.count('S')
    n_irreg = dssp_output.count('_')
    n_missing = dssp_output.count('?')

    positions = len( dssp_output) - 1  #'*'
    if positions == 0: return 100.0

    scoreblock = 100.0 / positions

    score = scoreblock * (n_S + n_irreg + n_missing + n_T)
    return score



def SetAsaAvg( result_file_entry, col_dict):

    result_file_entry[ col_dict['ASA_AVG']] = ""

    asa_avg_list = [result_file_entry[ col_dict['BIOL_ASA_AVG']], \
                    result_file_entry[ col_dict['ISOL_ASA_AVG']], \
                    result_file_entry[ col_dict['ASYM_ASA_AVG']]]

    for preferred_val in asa_avg_list:
        val = -1.0
        try: val = float( preferred_val)
        except: pass
        if val >= 0.0:
            result_file_entry[ col_dict['ASA_AVG']] = "%.2f" % val
            break


def SetTempfScore( grouped_entry, col_dict):

    #TEMPF
    t_mins = map( float, grouped_entry[ col_dict['TEMPF_FILE_MIN']].split(",")) #[float(i) for i in grouped_entry[col_dict['TEMPF_FILE_MIN']].split(",")]
    t_maxs = map( float, grouped_entry[ col_dict['TEMPF_FILE_MAX']].split(",")) #[float(i) for i in grouped_entry[col_dict['TEMPF_FILE_MAX']].split(",")]
    t_vals = map( float, grouped_entry[ col_dict['TEMPF_POI_VALS']].split(",")) #[float(i) for i in grouped_entry[col_dict['TEMPF_POI_VALS']]] #Avg Tempf for POI residues
    t_100_vals = []

    #Normalize tempf within file
    for v in range( len(t_mins)):
        min_ = t_mins[ v]
        max_ = t_maxs[ v]
        val_ = t_vals[ v]
        if min_ < 0.0001 : continue
        if max_ < 0.0001 : continue
        if (max_ - min_) < 0.0001 : continue
        #Normalize between 0.0-100.0
        normalized = 0.0 + (((val_-min_)*(100.0-0.0))/(max_-min_))
        t_100_vals.append( normalized)

    tempf_score = -1.0

    if len( t_100_vals):
        tempf_score = (sum( t_100_vals) / float(len( t_100_vals))) #AVG


    grouped_entry[ col_dict['TEMPF_SCORE']] = "%.2f" % tempf_score


def SetScore( g_entry, col_dict, disorder_dict = {}):

    #ASA val = Surface accessible area
    #normalized Tempf
    #Cleavage score
    #DSSP

    #Weights for each score
    mult = { "ASA":0.3, "TEMPF":0.1, "SEQ":0.3, "DSSP":0.3 }

    #All individual scores are in the range of 0.0 - 100.0

    #ASA
    asa_avg = float(g_entry[ col_dict['ASA_AVG']])
    #"WHAT?"
    #asa_score = ((1.0/12.0)*asa_avg-1.5)
    asa_score = (asa_avg*2.5)

    #Set between 0.0 and 100.0
    asa_score = min( 100.0, asa_score)
    asa_score = max( 0.0, asa_score)

    #TEMPF
    t_mins = [float(i) for i in g_entry[col_dict['TEMPF_FILE_MIN']].split(",")]
    t_maxs = [float(i) for i in g_entry[col_dict['TEMPF_FILE_MAX']].split(",")]
    t_vals = [float(i) for i in g_entry[col_dict['TEMPF_POI_VALS']].split(",")]
    t_100_vals = []

    #Normalize tempf within file
    for v in range( len(t_mins)):
        min_ = t_mins[ v]
        max_ = t_maxs[ v]
        val_ = t_vals[ v]
        if min_ < 0.0001 : continue
        if max_ < 0.0001 : continue
        if (max_ - min_) < 0.0001 : continue
        #Normalize between 0.0-100.0
        normalized =  0.0 + (((val_-min_)*(100.0-0.0))/(max_-min_))
        t_100_vals.append( normalized)

    tempf_score = -1.0

    if len( t_100_vals):
        tempf_score = (sum( t_100_vals) / float(len( t_100_vals))) #AVG

    #Normalize
    #tempf_score = tempf_score

    #Set between 0.0 and 100.0
    #tempf_score = min( 100.0, tempf_score)
    #tempf_score = max( 0.0, tempf_score)

    #SEQ
    #seq_score = float(g_entry[ col_dict['CLV']])

    #DSSP
    dssp_list = [ScoreDSSP( i) for i in g_entry[col_dict['DSSP']].split(",")]
    dssp_score = 0.0

    if len( dssp_list):
        dssp_score = sum( dssp_list) / len( dssp_list)


    score = asa_score    * mult["ASA"] + \
            tempf_score  * mult["TEMPF"] + \
            1.0    * mult["SEQ"] + \
            dssp_score   * mult["DSSP"]


    g_entry[ col_dict['TEMPF_SCORE']] = "%.2f" % tempf_score
    g_entry[ col_dict['DSSP_SCORE']] = "%.2f" % dssp_score
    g_entry[ col_dict['SCORE']] = "%.2f" % (score / 100.0)


def ProcessGroupedEntry( g_entry, col_types, col_map ):

    result = [] #[""]*len( col_types)

    #print "Process1:", g_entry[ col_map[ "SEQ_ID"]]

    for i in range( len( col_types)):
        col_name = col_types[ i][ 0]
        col_operation = col_types[ i][ 1] #operation to perform when grouping

        if col_operation == "special":
            result.append("")
            continue

        #Original col names for values to be grouped
        col_args = col_types[ i][ 2].split(",")
        if len( col_args) == 1 and len( col_args[0]) == 0:
            col_args[ 0] = col_name #grouped col name same as original

        #Make a list of all the values that need to be processed
        value_list = []
        value_errors = []
        for arg in col_args:
            #result_file_column_name = arg
            if arg not in col_map.keys():
                value_errors.append( "'%s' is not a result file column." % arg)
                #DEBUG
                #if arg == "SEQ_ID": print "VALUE:", g_entry[ col_map[ arg]]
            else:
                value_list.extend( g_entry[ col_map[ arg]])

        if len( value_errors):
            raise ValueError( "\n".join( value_errors))


        #DEBUG
        #if col_operation == "check_identity": print "Check:", g_entry[ 1]

        try:
            #Group Operation
            result.append( GroupOperation( value_list, col_operation))
            #DEBUG
            #if col_operation == "check_identity": print "Group_OP:", GroupOperation( value_list, col_operation)
        except ValueError as ex:
            sys.stderr.write( "FORMAT ERROR: Value error raised when processing column: %s\n" % col_name)
            message =         "Arguments: {0!r}\n".format( ex.args)
            sys.stderr.write( message)
            Exit( -4)
            Exit( -1)
            print value_list
            raise

    #"\t".join( result) + "\n"
    return result

#Deprecated
#ScoreResult( 38.98, 50.0, 0.8, 800.0, 100)
def ScoreResult( asa, tempf, clv, bitscore, align):

    #       asa  tempf clv   bs   align
    mult = [0.3, 0.3,  0.25, 0.1, 0.05]

    #http://www.wessa.net/rwasp_fitdistrnorm.wasp
    #mean =  38.9841116751269
    #dev = 24.9294262111256
    #tempf_score = min( 100.0, StandartDev( tempf, mean, dev) * 100.0)
    if tempf < 0.0001 and tempf > 0.0001:
        tempf_score = 50.0 #undefined
    else:
        tempf_score = min( 100.0, (-1.0*(1.0/900.0)*math.pow((tempf-50),2)+1)*100.0)
        tempf_score = max( 0.0, tempf_score)

    #print "TempF: %.3f" % tempf_score

    #mean = 20.7578540772532
    #dev = 13.1619983003219
    #asa_score = min( 100.0, StandartDev( asa, mean, dev) * 100.0)

    asa_score = min( 100.0, ((1.0/12.0)*asa-1.5)* 100.0 )
    asa_score = max( 0.0, asa_score)

    #print "ASA: %.3f" % asa_score

    #Quadratic x=0.5,y=0.0  and x=1.0,y=1.0
    clv_score = min( 100.0, (3.0*(clv*clv) - clv + 0.0) * 100.0)
    clv_score = max( clv_score, 0.0)

    #print "CLV: %.3f" % clv_score

    bitscore_score = min( 100.0, bitscore / 5)
    bitscore_score = max( bitscore_score, 0.0)

    #print "BITSCORE: %.3f" % bitscore_score

    #Linear x=70,y=0.0  and x=100,y=1.0
    align_score = min( 100.0, ((1.0/30.0)*align-(2.0+(1.0/3.0)))*100.0 )
    align_score = max( align_score, 0.0)

    #print "ALIGN: %.3f" % align_score

    return (asa_score*mult[ 0]+tempf_score*mult[ 1]+clv_score*mult[ 2]+bitscore_score*mult[ 3]+align_score*mult[ 4] )


#Operation types: any, avg, avg_nozeroes, min, max, list,special, count, delta, sum
def GroupOperation( value_list, optype):

    if len( value_list) == 0:
        sys.stderr("WARNING: empty value list given in GroupOperation")
        return ""

    if optype == "special":
        return ""
    elif optype == "any" or optype == "1st":
        return str(value_list[ 0])
    elif optype == "avg":
        if len( value_list) == 0: return ""
        return "%.4f" % (sum([float(i) for i in value_list]) / float(len( value_list)))
    elif optype == "avg_nozeroes":
        val = 0.0
        n = 0
        for v in value_list:
            fv = float(v)
            if fv < 0.00001 and fv > -0.00001: continue #Skip 0.0 values
            val += fv
            n += 1
        return "0.0" if n == 0 else str(val / float( n))
    elif optype == "avg_no_negatives":
        val = 0.0
        n = 0
        for v in value_list:
            fv = float(v)
            if fv < 0.0: continue #Skip negative values
            val += fv
            n += 1
        return "" if n == 0 else str(val / float( n))
    elif optype == "min":
        return str(min( [float(i) for i in value_list]))
    elif optype == "max":
        return str(max( [float(i) for i in value_list]))
    elif optype == "list":
        return ",".join(value_list)
    elif optype == "vote":
        options = list( set( value_list))
        counts = []
        for o in options:
            counts.append( value_list.count( o))
        max_count = int( max( counts))
        retval = ""
        for i in range( len( options)):
            if counts[ i] == max_count:
                retval += "%s%s" % (("/" if len( retval) else ""), options[ i])
        return retval
    elif optype == "set": #No duplicates
        return ",".join( list(set(value_list)))
    elif optype == "count":
        return str(len( value_list))
    elif optype == "delta":
        if not value_list or type( value_list) != list or len( value_list) == 0: return "0"
        vmin = min( [float(i) for i in value_list])
        vmax = max( [float(i) for i in value_list])
        return str((vmax - vmin))
    elif optype == "delta_no_negatives":
        if not value_list or type( value_list) != list or len( value_list) == 0: return "0"
        filtered_value_list = filter( lambda a: float( a) >= 0.0, value_list )
        if len( filtered_value_list) == 0: return "0"
        vmin = min( map(float, filtered_value_list))
        vmax = max( map(float, filtered_value_list))
        return str((vmax - vmin))
    elif optype == "sum":
        return str(sum([float(i) for i in value_list]))
    elif optype == "matching":
        #Return string of the matching part of chars in value_list strings
        min_len = min( map( len, value_list))
        cut_off_pos = 0
        for c in range( min_len):
            matching = True
            for v in range( len( value_list)-1):
                if value_list[ v][ c] != value_list[ v+1][ c]:
                    matching = False
                    break
            if matching:
                cut_off_pos = c+1
            else:
                break
        while cut_off_pos > 0:
            last_char = value_list[ 0][ cut_off_pos-1]
            if last_char == " " or last_char == "_":
                cut_off_pos -= 1
            else:
                break
        return str( value_list[ 0][0:cut_off_pos])
    elif optype == "check_identity":
        for i in range( len( value_list)-1):
            if value_list[ i] != value_list[ i+1]:
                raise ValueError( "Identity check failed: '%s' vs. '%s'" % (value_list[ i], value_list[ i+1]))

        return str(value_list[ 0])

    raise ValueError( "Unknown optype encountered: %s" % optype)


def StandartDev( x, mean, dev):
#def sd( x, mean, dev):
#StandartDev( 38.9841116751269, 38.9841116751269, 24.9294262111256)

    power = -1.0* (((x-mean)*(x-mean)) / (2*dev*dev))
    return (1 / (dev*math.sqrt(2*math.pi)))*math.pow( math.e, power)


def PrintHelp():

    script_file = ""# os.path.basename(__file__)
    cores = 0

    try:
        script_file = os.path.basename(__file__)
        cores = multiprocessing.cpu_count()
    except: pass

    if not len( script_file): script_file = "score_poisites.py"

    print u"""
    This script finds and calculates metrics of structural data for marked points of interest (POI) in input fasta sequences.
    Default functionality is to include the amino acids on both sides of the marked POI.
    Usage: %s [-h|-a][-s][-v]|[-r][-n][-w INT][-p INT][-b INT][-t INT][-c FILE] fasta_sequence_file   [output_folder]
    or:    %s                                                                   fasta_sequence_folder [output_folder]

    POIs can be indicated with an asterisk ("*") in the sequence. Each sequence can contain multiple POIs. The fasta sequences
    should  have unique idenfiers separated with a "|" in the fasta sequence header. If the fasta header sequence description
    contains a GN tag ("GN=XXXX"), XXXX will be included in the GENE column in the results. Custom result file columns can be
    created by inserting values in the format of "|MYSCORE:1.23|" to the sequence fasta headers after the sequence descriptions.

    -s  --sidechain        Include only sidechains in the calculations. By default, backbone atoms are included in calculations.
                           If the marked POIs have glycines (no sidechain) they will appear always as buried residues with no
                           accessible surface are (ASA) when this option is enabled.
    -p  --positions INT    Consider INT amino acids before (N-terminal side) the marked POI in calculations. Default value (p=0)
                           is special, and considers the two amino acids on either side of the POI.
    -b  --blastres INT     Maximum number of blast results to use for each POI (0 for unlimited). [default: 30]
                           Using more resulting structures extends runtime due to more required alignments to locate POIs in
                           the structure. For large proteins, where only partial structures are available, limiting the number
                           of results with full sequence BLAST searches can result in not all of the available structural
                           information being used. Some structures can have more than 60 homologs deposited.
    -w  --blastwin INT     Number of amino acids around the POIs to use for BLASTing [default: 0, full sequences]
                           A minimum window of 20 amino acids is recommended for windowed BLAST searches.
                           If you have many POIs per sequence (>2), BLASTing with full sequences is significantly faster.
                           If you have large proteins that have many partial structures, you might miss some of them when
                           BLASTing with full sequences and non-zero --blastres param.
    -n  --no_cache         Do not cache BLAST results into sqlite database
    -i  --biomolecules     Write used biological assemblies to file in result folder
    -r  --reset            Reset progress and remove any existing files
    -l  --no_paralign      Do not create alignments parallelly with blast queries
    -c  --config           Path to config.ini file if in other than current dir
    -a  --autoconfig       Try to automatically configure config.ini paths
    -v  --verbose          More progress output
    -h  --help             Print this message.
    -t  --threads          How many threads to use for ASA calculations [default:0 (all)]
                           Negative values can be used to limit thread count from max cores.
                           This system has %i available cores.
    -f  --format           Only reformat results
    -d  --debug            Print additional information in case of errors

    In addition to these options, several values can be configured through the "config.ini" file. If a "config.ini" file is
    not found, a file with default values will be created.


    Authors                Anssi Nurminen, Vesa P. Hytönen
                           University of Tampere
    Version                1.0c; Feb 2018
    """ % (script_file," "*len(script_file), cores)



def main():

    #See beginning of file
    global DEFAULT_CONFIG_FILE, CONFIG_VALUES, CUSTOM_COLS_FILE

    atexit.register( ExitHandler)
    config_file = DEFAULT_CONFIG_FILE

    reset = False
    threads = 0 #all available
    only_format = False
    verbose = False
    parallel_align = True
    blast_caching = True
    sidechain_only = False
    seq_file = ""
    output_folder = ""
    positions = 0
    autoconfig_set = False
    blast_results = -1
    blast_window = 0 #zero = blast with full seq
    debug = False
    write_bioassemblies = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'p:t:c:b:w:srfhvdnladbi', ['positions', 'threads', 'config', 'blastres', 'blastwin', 'sidechain', 'reset', 'format', 'help', 'version', 'verbose', 'no_cache', 'no_paralign', 'autoconfig', 'debug', 'biomolecules', 'biomols'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        sys.stderr.write( "See -h for usage.\n")
        Exit( 2)


    #Input & Output files
    for arg in args:
        if len( seq_file) == 0:
            seq_file = arg
        elif len( output_folder) == 0:
            output_folder = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            Exit( 2)

    #Use current dir if it was not specified
    if len( output_folder) == 0: output_folder = "."

    #Flags
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            PrintHelp()
            Exit( 0)
        elif opt in ('-r', '--reset'):
            if not reset: print "RESET"
            reset = True
        elif opt in ('-s', '--sidechain'):
            sidechain_only = True
        elif opt in ('-b', '--blastres'):
            try:
                blast_results = int( arg)
                if blast_results < 0:
                    sys.stderr.write( "ERROR: Invalid blastres (-b) argument %s: '%i'.\n" % ( opt,blast_results))
                    Exit( -2)
            except ValueError:
                sys.stderr.write( "ERROR: Invalid blastres (-b) argument %s: '%s'. Positive integer or zero expected.\n" % (opt,arg))
                Exit( -2)
        elif opt in ('-w', '--blastwin'):
            try:
                blast_window = int( arg)
                if blast_window < 0:
                    sys.stderr.write( "ERROR: Invalid blastwin (-w) argument %s: '%i'.\n" % ( opt,blast_window))
                    Exit( -2)
                blast_window += blast_window % 2 #Make even if not zero
            except ValueError:
                sys.stderr.write( "ERROR: Invalid blastwin (-w) argument %s: '%s'. Positive integer or zero expected.\n" % (opt,arg))
                Exit( -2)
        elif opt in ('-v', '--verbose'):
            if not verbose: print "VERBOSE"
            verbose = True
        elif opt in ('-l', '--no_paralign'):
            parallel_align = False
        elif opt in ('-n', '--no_cache'):
            blast_caching = False
        elif opt in ('-i', '--biomolecules', '--biomols'):
            write_bioassemblies = True
        elif opt in ('-f', '--format'):
            only_format = True
        elif opt in ('-d', '--debug'):
            if not debug: print "DEBUG"
            debug = True
        elif opt in ('-t', '--threads'):
            try:
                threads = int(arg)
            except ValueError:
                sys.stderr.write( "WARNING: Invalid thread argument '%s'. Using default thread count.\n" % (arg))
            if threads == 0:
                available_cores = multiprocessing.cpu_count()
                if available_cores > 64:
                    sys.stdout.write( "INFO: This system is detected to have %i cores.\n" % available_cores)
                    sys.stdout.write( "INFO: To use more than 64 cores, specify exact number with \"-t CORES\".\n")
                    sys.stdout.write( "INFO: Using 64 cores.\n" )
                    threads = 64
                else:
                    sys.stdout.write( "INFO: Using all available %i cores.\n" % available_cores)
            elif threads < 0:
                if multiprocessing.cpu_count() + threads < 1:
                    sys.stdout.write( "WARNING: Only %i cores available, cannot limit by %i. Using only a single core.\n" % (multiprocessing.cpu_count(), (threads*-1)))
                    threads = 1
                else:
                    limit_cores_by = -1*threads
                    threads = multiprocessing.cpu_count()+threads
                    if limit_cores_by == 1:
                        sys.stdout.write( "INFO: Using all but %i core => %i cores.\n" % (limit_cores_by,threads))
                    else:
                        sys.stdout.write( "INFO: Using all but %i cores => %i cores.\n" % (limit_cores_by,threads))
            else:
                sys.stdout.write( "INFO: Using max %i threads.\n" % (threads))

        elif opt in ('-c', '--config'):
            config_file = arg
        elif opt in ('-a','--autoconfig'):
            autoconfig_set = True
        elif opt in ('-p', '--positions'):
            try:
                positions = int( arg)
                if positions < 0:
                    sys.stderr.write( "ERROR: Invalid position argument %s: '%i'.\n" % ( opt,positions))
                    Exit( -2)
            except ValueError:
                sys.stderr.write( "ERROR: Invalid position argument %s: '%s'.\n" % (opt,arg))
                Exit( -2)
        else:
            sys.stderr.write("WARNING: Unknown option ignored: '%s'.\n" % opt)

    #Autoconfigure config.ini
    if autoconfig_set:
        sys.stdout.write("INFO: Creating configuration file '%s'.\n" % config_file)
        sys.stdout.write("INFO: Make sure you have the necessary components installed for the auto-configuration to be able to find them.\n")
        sys.stdout.write("INFO: Search has a 15 second timeout for each file, which may not be enough for a thorough search.\n\n")
        WriteConfig( config_file, aAutoConfig=True )
        Exit( 0)

    if len( config_file) > 1 and not os.path.isfile( config_file):
        sys.stderr.write("INFO: No configuration file was found. Creating a default file: '%s'.\n" % config_file)
        sys.stderr.write("INFO: Try option --autoconfig after required components such as BLAST have been installed.\n" % config_file)
        WriteConfig( config_file, aAutoConfig=False )
        Exit( 0)

    if len( sys.argv) < 1:
        sys.stderr.write("ERROR: Too few arguments.\n")
        Exit( -2) #Exit

    if len( seq_file) == 0:
        sys.stderr.write("ERROR: No input sequence file specified.\n\n")
        PrintHelp()
        #sys.stderr.write("See -h for usage.\n")
        Exit( 0)

    if not os.path.isfile( seq_file):
        sys.stderr.write("ERROR: '%s' does not specify a valid file\n" % seq_file)
        Exit( -2)


    config = ReadConfig( config_file) #Exits if not ok
    CONFIG_VALUES = config

    precalc_asa_folder = config['PRECALC_ASA_DIR'] if not sidechain_only else config['PRECALC_ASA_SC_DIR']
    pdb_folder   = config['PDB_DIR']
    blast_db_dir = config['BLAST_DB']
    dssp_folder  = config['DEFAULT_DSSP_DIR']
    #human_proteome_file = config['DEFAULT_HUM_PROTEOME']
    FetchDSSP.DEFAULT_DSSP_FTP = config['DEFAULT_DSSP_FTP']

    if precalc_asa_folder and len( precalc_asa_folder) and not os.path.exists( precalc_asa_folder): os.makedirs( precalc_asa_folder)

    if verbose: sys.stdout.write( "INFO: Configuration file '%s' read successfully.\n" % config_file)

    #Positions
    if positions != config['DEFAULT_POI_POSITIONS']:
        sys.stdout.write( "INFO: Config file POI positions overridden on commandline (%i => %i).   [OK]\n" % (config['DEFAULT_POI_POSITIONS'], positions))
        config['DEFAULT_POI_POSITIONS'] = positions

    #Sidechains
    if sidechain_only: config['DEFAULT_ONLY_SIDECHAINS'] = True
    else: sidechain_only = config['DEFAULT_ONLY_SIDECHAINS']

    if verbose and sidechain_only: sys.stdout.write( "INFO: Considering only sidechains in structure calculations.\n")


    #Use default folder for output if it was not specified
    if len( output_folder) == 0:
        #path,filename,name,suffix = SplitFilePath( seq_file )
        #output_folder = InsertFolderSeparator( seq_path + seq_name)
        output_folder = InsertFolderSeparator( config['DEFAULT_RESULT_DIR'])

    if output_folder == "./" or output_folder == ".\\":
        output_folder = "./" + seq_name + "/"


    output_folder = InsertFolderSeparator( output_folder)
    config['DEFAULT_RESULT_DIR'] = output_folder

    seq_path,seq_filename,seq_name,seq_suffix = SplitFilePath( seq_file )

    align_folder = InsertFolderSeparator( output_folder + "align" )
    result_folder = InsertFolderSeparator( output_folder + "results" )
    disorder_file = output_folder + "disorder_predictions.txt"
    custom_cols_file = output_folder + CUSTOM_COLS_FILE
    if reset:
        try: os.remove( custom_cols_file)
        except: pass

    #Create directories if needed
    if not os.path.exists( output_folder): os.makedirs( output_folder)
    if not os.path.exists( align_folder):  os.makedirs( align_folder)
    if not os.path.exists( result_folder): os.makedirs( result_folder)

    WriteProgress( "MAIN: Process started.", aNewFile=True)

    #Make sure input sequences are formatted correctly
    formatted_input_sequences_file = output_folder + "input_sequences.fasta"
    if reset or not os.path.isfile( formatted_input_sequences_file):

        print ""
        WriteProgress( "MAIN: Formatting input sequences...")
        print ""
        retval = ConvertFastaLocatorFormat( seq_file, formatted_input_sequences_file, aStar=False, aHeader=True, aMarker="between")
        gc.collect()
        if retval != True:
            sys.stdout.write( "ERROR: Fatal error while formatting input sequences.\n")
            Exit( -42) #Exit

    seq_file = formatted_input_sequences_file

    blast_file = ChangeToFolder( seq_file + ".blast", output_folder)
    bit_score_threshold = config["BITSCORE_THRESHOLD"]
    if blast_window > 0: bit_score_threshold = 0.0


    if not only_format:

        if reset:
            print "INFO: Emptying folder '%s'." % result_folder
            ClearFolder( result_folder, aClearSubfolders=True)
            print "INFO: Emptying folder '%s'." % align_folder
            ClearFolder( align_folder, aClearSubfolders=True)
            #sys.exit( 0)
        print "OBSOLETES: %s" % pdb_folder
        obsoletes = online_pdb.FetchObsoletes( True, pdb_folder)


        #1. Find homologs using BLAST
        #2. Align structure sequences with fasta sequences
        #3. Identify poisites in the structures
        #4. Get temp factors and calc ASA for files

        print ""
        WriteProgress( "MAIN: Finding structure homologs using BLAST...")
        print ""

        logging.basicConfig( filename=output_folder+'blast.log',level=logging.DEBUG)
        #logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)

        if blast_results < 0:
            blast_results = config['MAX_BLAST_HOMOLOGS']

        if blast_results <= 0: sys.stdout.write( "INFO: Using all available blast results.\n")
        else: sys.stdout.write( "INFO: Using max %i blast results.\n" % (blast_results))

        if blast_window <= 0: sys.stdout.write( "INFO: Using full sequences for BLASTing.\n")
        else: sys.stdout.write( "INFO: Using %i amino acids around POIs for BLASTing.\n" % (blast_window))


        #BLAST EM! (and fetch PDB structures if necessary)
        retval = -1
        try:
            retval = FindStructures( config['DEFAULT_BLASTP'], seq_file, blast_file, pdb_folder, blast_db_dir, parallel_align=parallel_align,\
                                     max_threads=threads, max_results=blast_results, e_val_threshold=0.001, fetch_files=True, reset_progress=reset,\
                                     obsoletes=obsoletes, verbose=verbose, sequence_window=blast_window, use_repository=blast_caching, debug=debug)
        except KeyboardInterrupt:
            #print "Process interrupted by user. Exitting..."
            Exit( -99) #Exit
        except Exception as ex:
            sys.stderr.write( "FATAL ERROR: An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
            if debug: raise

        if retval < 0:
            sys.stderr.write("ERROR: Fatal error occurred when finding homologous PDB structures using BLAST.\n")
            Exit( retval) #Exit

        gc.collect()


        print ""
        WriteProgress( "MAIN: Creating alignment files...")
        print ""

        logging.basicConfig( filename=output_folder+'score_run.log',level=logging.DEBUG)

        #ALIGN 'EM!
        #Creates alignment files that include an alignment score in range 0-100
        #Argument: score_it=True is obsolete, used only with ClustalW
        try:
            retval = ProcessBlastResults( seq_file, blast_file, align_folder, pdb_folder, obsoletes=obsoletes, skip_existing=(not reset), score_it=False, threads=threads )
        except Exception as ex:
            sys.stderr.write( "FATAL ERROR: An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
            if debug: raise
            retval = -89

        if retval == -99:
            #Process interrupted by user
            Exit( retval) #Exit
        if retval < 0:
            sys.stderr.write("ERROR: Fatal error occurred when aligning POI-sites to PDB structures.\n")
            Exit( retval) #Exit

        gc.collect()

        print ""
        WriteProgress( "MAIN: Finding POI sites in structures...")
        print ""

        #Find POI sites in structures
        try:

            retval = ProcessPoiSites( seq_file, GetFolder( blast_file), align_folder, result_folder, pdb_folder, precalc_asa_folder, dssp_folder, \
                                      obsoletes, positions=positions, reset=reset, bit_score_cutoff=bit_score_threshold, align_score_cutoff=config["ALIGNMENT_THRESHOLD"],\
                                      threads=threads, aCustomColsFile=custom_cols_file, sidechain_only=sidechain_only, aWriteUsedBiomolecules=write_bioassemblies,
                                      aDSSPExecutable=CONFIG_VALUES['DEFAULT_DSSP']) #Writes custom cols file
        except Exception as ex:
            sys.stderr.write( "FATAL ERROR: An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
            if debug: raise
            retval = -89

        if retval == -7:
            sys.stderr.write("ERROR: No results.\n")
            Exit( retval) #Exit
        elif retval < 0:
            sys.stderr.write("ERROR: Fatal error occurred when processing POI-sites.\n")
            Exit( retval) #Exit

        gc.collect()

        if len( config["DEFAULT_IUPRED"]):
            if os.path.isfile( config["DEFAULT_IUPRED"]):
                print ""
                WriteProgress( "MAIN: Running disorder predictions")
                print ""
                retval = PredictDisorder( config["DEFAULT_IUPRED"], seq_file, disorder_file, aSeqWindow=config['IUPRED_SEQWINDOW'], aForceReset=reset, aWorkDir=config["WORK_DIR"])

                if retval < 0:
                    sys.stderr.write("ERROR: Fatal error occurred when running disorder predictions. Quitting.\n")
                    Exit( retval) #Exit
            else:
                sys.stderr.write( "MAIN: IUPred executable not found at '%s'.\n" % config["DEFAULT_IUPRED"])
        elif verbose:
            sys.stdout.write( "MAIN: IUPred executable not defined in file 'config.ini'. No disorder predictions are included in the results.\n")

    print ""
    WriteProgress( "MAIN: Collecting results...")
    print ""
    results_file = output_folder + seq_name + "_results.txt"
    collect_log_file = output_folder + seq_name + "_collect_log.txt"
    custom_cols = ReadListFromFile( custom_cols_file, aRemoveDuplicates=True)
    if len( custom_cols): WriteListToFile( custom_cols, custom_cols_file, aAppend=False) #Get rid of possible duplicates

    with open( collect_log_file, "w+") as collect_log:
        #Goes through result pdb files and collects the data to a single file
        CollectResults( result_folder, results_file, seq_file, collect_log, aCustomColumns=custom_cols)

    print ""
    WriteProgress( "MAIN: Formatting results (grouping entries)...")
    print ""
    grouped_results_file = output_folder + seq_name + "_grouped_results.txt"
    retval = FormatResults( input_file=results_file, output_file=grouped_results_file, disorder_file = disorder_file, aCustomColumns=custom_cols )

    if retval < 0:
        sys.stderr.write("ERROR: Fatal error occurred when formatting results.\n")
        Exit( retval) #Exit


    print ""
    WriteProgress( "MAIN: Writing results summary...")
    print ""
    summary_file = output_folder + seq_name + "_summary.txt"
    retval = -1
    retval = SummarizeResults( seq_file, blast_file, bit_score_threshold, output_folder, grouped_results_file, aOutputFile=summary_file, aDisorderFile = disorder_file, aObsoletes=online_pdb.FetchObsoletes( True, pdb_folder), aCustomCols=custom_cols )

    if retval < 0:
        sys.stderr.write("ERROR: Fatal error occurred when summarizing results.\n")
        Exit( retval) #Exit


    print ""
    #print "MAIN: All done."
    Exit( 0)

def WriteProgress( aMsg, aPrintIt=True, aNewFile=False):

    global CONFIG_VALUES

    if len( aMsg) == 0: return
    elif aMsg[ -1] != "\n": aMsg = aMsg + "\n"

    if aPrintIt: sys.stdout.write( aMsg)

    if 'DEFAULT_RESULT_DIR' in CONFIG_VALUES and len( CONFIG_VALUES['DEFAULT_RESULT_DIR']):
        progress_file = "progress.txt"
        progress_file = InsertFolderSeparator( CONFIG_VALUES['DEFAULT_RESULT_DIR']) + progress_file

        filemode = "w" if aNewFile else "a"

        try:
            with open( progress_file, filemode) as pf:
                pf.write( aMsg)
        except Exception:
            sys.stderr.write( "ERROR: Writing progress file '%s' failed.\n" % progress_file)
    else:
        sys.stderr.write( "ERROR: No result folder set for writing progress.\n" % progress_file)

def ExitHandler():
    global UNEXPECTED_EXIT
    if UNEXPECTED_EXIT:
        try: WriteProgress( "MAIN: Unexpected exit.")
        except: pass


def Exit( aExitCode):

    global UNEXPECTED_EXIT

    try:
        if aExitCode == 0:
            WriteProgress( "MAIN: DONE. [%s]" % str(aExitCode), aPrintIt=False)
        elif aExitCode == -99:
            WriteProgress( "MAIN: Process Interrupted by user. [%s]" % str(aExitCode))
        else:
            WriteProgress( "MAIN: FAILED. [%s]" % str(aExitCode))
    except:
        pass

    UNEXPECTED_EXIT = False
    sys.exit( aExitCode)


if __name__ == "__main__":

  main()

