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
from subprocess import call

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

sys.path.append(  "%s/../utils" % os.path.dirname(os.path.abspath(__file__)))
from fasta_formatter import GetMarkerLocation
from fileops import *
from struct_methods import *
import progress_bar as pb

#Returns dictionary with seq_id as key
def ReadDisorderResults( aDisorderFile):

    disorder_results = {}
    #re_header = re.compile( "^\s*>.+?\|([^\|]+)\|")
    re_result = re.compile( "^\s*disorder.+?:\s*([\d\.]+)")
    #re_cutpos = re.compile( "(?:POI):\s*(\d+)")

    cur_value = ""
    diso_records = None

    try:
        #Not a fasta file but ">" at the header lines works
        diso_records = FastaIter( aDisorderFile )
    except IOError as e:
        sys.stderr.write( "DISORDER ERROR: Could not open disorder results file '%s'.\n" % ( aDisorderFile))
        sys.stderr.write( "DISORDER ERROR: %s\n" % ( e.strerror))
        sys.exit(-2)
        #return disorder_results

    for d_header, d_result in diso_records:

        poipos = int( GetRecordPOI( d_header))

        if poipos < 0:
            sys.stderr.write( "DISORDER ERROR: Entry '%s' has no POI specified in file '%s'.\n" % (d_header, aDisorderFile))
            sys.exit(-2)

        m = re_result.search( d_result) #Is it a result line

        if m != None:
            cur_value = float( m.group( 1))

            #Insert _Pos\d+ if needed
            entry_id = GetRecordACC( d_header, aWithPos=True) #cur_seq if (cur_seq.find("_Pos") >= 0) else (cur_seq + "_Pos" + cur_pos)

            if entry_id in disorder_results:
                sys.stderr.write( "DISORDER WARNING: Duplicate sequence id: '%s' on line %i. Duplicate result ignored.\n" % (entry_id, linenum))
            else:
                disorder_results[ entry_id] = cur_value


    print "DISORDER INFO: Read disorder prediction results for %i entries." % len( disorder_results)

    return disorder_results



def PredictDisorder( aIUPred, aInputFile, aOutputFile, aSeqWindow = 8, aForceReset=False, aWorkDir=None):

    print "IUPRED: Initializing IUPRED..."

    #predictor_exe = "..\..\disorder\iupred\iupred.exe"
    predictor_exe = aIUPred
    abspath = os.path.abspath( predictor_exe)
    if os.name == 'nt':
        os.environ["IUPred_PATH"] = abspath[:abspath.rfind( "\\")]
    else:
        os.environ["IUPred_PATH"] = abspath[:abspath.rfind( "/")]
    #print "Path:" + abspath[:abspath.rfind( "\\")]
    #sys.exit(0)

    if aForceReset:
        try: os.remove( aOutputFile)
        except: pass

    if aSeqWindow <= 0:
        sys.stderr.write("IUPRED ERROR: invalid sequence window specified: %i.\n" % aSeqWindow)

    seq_records = list( SeqIO.parse( aInputFile, "fasta"))

    work_dir = os.path.dirname(os.path.abspath(__file__))
    if aWorkDir and len( aWorkDir): work_dir = aWorkDir

    temp_inputfile = work_dir + "/temp_in.fasta.tmp"
    temp_outputfile = work_dir + "/temp_out.fasta.tmp"

    temp_inputfile  = SystemSpecificFolderSeparators( temp_inputfile)
    temp_outputfile = SystemSpecificFolderSeparators( temp_outputfile)

    re_clv = re.compile( "(POI):\s*[\d\.]+")
    n_cleavage_sites = sum( [ (1 if re_clv.search( x.description) != None else 0) for x in seq_records ])

    if len( seq_records) == 0:
        print "IUPRED ERROR: File '%s' did not contain any fasta sequences." % aInputFile
        return -1
    elif n_cleavage_sites == 0:
        print "IUPRED ERROR: No sequences with a POI found in file '%s'." % aInputFile
        return -1

    print "IUPRED: File '%s' has %i/%i (%.1f%%) sequences with a known cleavage site." % (aInputFile, n_cleavage_sites, len(seq_records), float(n_cleavage_sites)/len(seq_records)*100.0)

    #Check if work already done
    if os.path.isfile( aOutputFile):

        n_seqs = 0

        seq_headers = FastaHeaderIter( aOutputFile)
        for seq_header in seq_headers: n_seqs += 1

        if n_seqs == n_cleavage_sites:
            print "IUPRED: Outputfile already exists with a correct number of sequences. (%i/%i matching)" % (n_seqs,n_cleavage_sites)
            print "IUPRED: Skipping disorder prediction. Use reset, or delete file '%s' to redo calculations." % aOutputFile
            return 0
        else:
            print "IUPRED: Overwriting exiting file. (%i/%i matching)" % (n_seqs,n_cleavage_sites)


    print "IUPRED: Starting disorder prediction for %i sequences..." % n_cleavage_sites #len( seq_records)
    print "IUPRED: If this process is interrupted it will be restarted from the beginning."
    print "IUPRED: Writing output to file '%s'." % aOutputFile
    progress = 0
    missing = 0
    total = n_cleavage_sites #len( seq_records)

    predictions = {}
    re_resultrow = re.compile("^\s*(\d+)\s(.)\s+([\d\.]+)$", re.IGNORECASE)

    pb.InitProgressBar( progress, total, aUpdateIntervalSeconds=2.0)

    with open( aOutputFile, "w+") as outf:

        for rec in seq_records:

            markers = GetMarkerLocation( rec)
            if len( markers) == 0:
                if rec.description.find( "No_information") < 0:
                    pb.PrintAboveProgressBar( "IUPRED WARNING: No sequence marker found for record '%s'.\n" % (rec.description))
                else:
                    missing += 1
                continue
            if len( markers) > 1:
                pb.PrintAboveProgressBar( "IUPRED WARNING: Multiple sequence markers found for record '%s'.\n" % (rec.description))
                continue

            cleaned_seq = str( rec.seq)
            cleaned_seq = cleaned_seq.replace( "-", "")
            cleaned_seq = cleaned_seq.replace( "*", "")

            location = markers[ 0]

            results = {}

            if cleaned_seq in predictions.keys():
                results = predictions[ cleaned_seq]
            else:
                rec.seq = Seq( cleaned_seq, IUPAC.protein)

                SeqIO.write( rec, temp_inputfile, "fasta")
                call_str = predictor_exe + " " + temp_inputfile + " long > " + temp_outputfile
                call_str = SystemSpecificFolderSeparators( call_str)
                pb.PrintAboveProgressBar( "Predicting: " + GetRecordACC( rec.id, aWithPos=False))

                #PREDICT
                return_code = -1

                try:
                    return_code = call( call_str, shell=True)
                except Exception as ex:
                    sys.stderr.write("IUPRED ERROR: Could not run iupred '%s'.\n" % call_str)
                    sys.stderr.write( type( ex).__name__ + ": '" + str( ex)+"'\n")
                    raise

                if return_code != 0:
                    pb.PrintAboveProgressBar( "IUPRED ERROR: '%s' returned with code: %i.\n" % (predictor_exe, return_code))
                    return return_code


                with open( temp_outputfile, "r") as analysis_file:
                    content = analysis_file.readlines()

                #Read results
                for line in content:
                    m_result = re_resultrow.match( line)
                    if m_result:
                        results[ int( m_result.group( 1))] = float( m_result.group( 3))

                predictions[ cleaned_seq] = results

            disorder_avg = 0.0
            range_start = max( 1, location - aSeqWindow/2) #Numbering starts from 1
            range_end = min( len(rec.seq), location + aSeqWindow/2 + aSeqWindow % 2)


            for i in range( range_start, range_end):
                try:
                    disorder_avg += results[ i]
                except KeyError:
                    pb.PrintAboveProgressBar( "IUPRED ERROR: Error processing sequence '%s'.\n" % (rec.id))
                    pb.PrintAboveProgressBar( "IUPRED ERROR: Results did not contain residue '%i'.\n" % (i))
                    pb.PrintAboveProgressBar( "IUPRED ERROR: Residue range: %i - %i.\n" % (range_start, range_end))
                    pb.FinalizeProgressBar( progress, total)
                    sys.exit( -1) #EXIT

            disorder_avg /= float( range_end - range_start)

            outf.write( ">%s\ndisorder_@%i: %.4f\n\n" % ( rec.description, location, disorder_avg))
            progress += 1

            #Report progress
            pb.ProgressBar( progress, total)
            #if progress % 20 == 0:
            #    print "DISORDER INFO: Progress %i/%i (%.1f%%)" % (progress, total, float(progress)/total*100.0)

    pb.FinalizeProgressBar( progress, total)
    print "IUPRED INFO: Disorder predictions finished."

    try: os.remove( temp_inputfile)
    except: pass

    try: os.remove( temp_outputfile)
    except: pass

    if missing > 0:
        print "IUPRED INFO: %i sequences were missing cutsite information." % missing
    else:
        print "IUPRED INFO: All cutsites were found. [OK]"

    return 0

