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



import re
import sys
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

#sys.path.append( "../utils" )


from struct_methods import *


AAS = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","U","X"]
#U == selenocysteine
RE_CLEAN_PATTERN = re_cleanpattern = re.compile(r'[\s-]+')

FASTA_HEADER_VALUES_PAT = re.compile( "\s*([A-Z]+):(.*)\s*")
RESERVED_FASTA_HEADER_COLS = [ "ENTRY", "ENTRY_ID", "GENE", "POI", "POISEQ", "STRUCTURES", "ACCESSIBLE", "FLEXIBLE", "BLAST_ENTRIES", "BEST_BITSCORE", "BITSCORE_PASS", "BEST_ALIGNMENT", "ALIGN_PASS", "DISORDER_SCORE"
                               "DATASERIES","SEQ_ID","RELIABILITY","GENE","DESC","N_STRUCTURES","PYMOL","IUPRED","DSSP","ASA_VALS","ASA_AVG","ASA_DECISION","ASA_POI_DELTA","ASA_POI_MIN","ASA_POI_MAX","TEMPF_POI_VALS",
                               "TEMPF_POI_MIN","TEMPF_POI_MAX","TEMPF_FILE_MIN","TEMPF_FILE_MAX","PDB_FILES","PDB_CHAINS","PDB_SITES","PDB_PTMS","BLAST_RANKS","ENTRY","SEQ_POIPOS","SEQ_POISEQ","PDBSEQ","ALIGNMENTS","ALIGN_AVG",
                               "BITSCORE_AVG","IDENTITY_AVG","METHODS","R-FREE_AVG","RESOLUTIONS","ORGANISM_SCIENTIFIC","ORGANISM_TAXID" ]

def ExtractCustomValuesFromFastaHeader( aHeader, aAllValues=False ):

    global FASTA_HEADER_VALUES_PAT

    retvals = {}
    cols = aHeader.split("|")

    #pattern: |XXX:0000|
    pattern = FASTA_HEADER_VALUES_PAT

    for c in cols:
        m = pattern.match( c)
        if m:
            if not aAllValues and (m.group( 1) == "POI" or m.group( 1) == "POISEQ"): continue
            retvals[ m.group( 1)] = m.group( 2)

    return retvals


def ReservedHeaderCustomCols( aCustomCols):

    global RESERVED_FASTA_HEADER_COLS

    retvals = []

    for custom_col in aCustomCols:
        if custom_col in RESERVED_FASTA_HEADER_COLS:
            retvals.append( custom_col)

    return retvals


def GetStrMarkerPositions( aSeq, aMarker="between"):

    global AAS, RE_CLEAN_PATTERN

    markers = []
    seq_pos = 0 #seq starts from 1
    #re_cleanpattern = re.compile(r'[\s-]+')
    cleaned_seq = re.sub( RE_CLEAN_PATTERN, '', str( aSeq))

    #Find markers
    for si in range( len( cleaned_seq)):
        aa = cleaned_seq[ si]
        if aa in AAS:
            seq_pos += 1
        elif aa != '*':
            sys.stderr.write( "WARNING: Sequence contains unknown letter '%s'.\n" % ( aa))
            seq_pos += 1
        else:
            #Marker found
            pos_mod = 0
            if aMarker == "before": pos_mod = -1

            markers.append( seq_pos + pos_mod) #is '*'

    return markers



def GetMarkerLocation( aSeqRecord, aMarker="between"):

    global AAS

    re_cutnum = re.compile('POI:\s*(\d+)-(\d+)', re.IGNORECASE)
    re_cleanpattern = re.compile(r'[\s-]+')

    res = re_cutnum.search( aSeqRecord.description)
    header_contains_location = (res != None)
    marked_location_in_header = -1

    if header_contains_location:
        cut_1 = int( res.group( 1))
        cut_2 = int( res.group( 2))
        marked_location_in_header = cut_1
        return [ int(marked_location_in_header)]


    markers = []
    seq_pos = 0 #seq starts from 1
    cleaned_seq = re.sub( re_cleanpattern, '', str(aSeqRecord.seq))

    #Find markers
    for si in range( len( cleaned_seq)):
        aa = cleaned_seq[ si]
        if aa in AAS:
            seq_pos += 1
        elif aa != '*':
            sys.stderr.write( "WARNING:'" + rec.description + "' contains unknown letter '%s'.\n" % ( aa))
            seq_pos += 1
        else:
            pos_mod = 0
            if aMarker == "before": pos_mod = -1

            markers.append( seq_pos + pos_mod) #is '*'

    return markers

def ConvertFastaLocatorFormat( aInputFile, aOutputFile, aStar=False, aHeader=True, aMarker="between"):

    global AAS

    records = list( SeqIO.parse( aInputFile, "fasta"))

    sys.stderr.write( "INFO: Read %i sequence records from file '%s'.\n" % (len( records), aInputFile))

    re_cutnum = re.compile('POI:\s*(\d+)-(\d+)', re.IGNORECASE)
    re_cutseq = re.compile('\s*POISEQ:\s*([A-Z]+)/([A-Z]+)', re.IGNORECASE)
    re_pattern = re.compile(r'[\s-]+') #to remove dashes and spaces
    re_id = re.compile(r'^([^\|]+\|)([^\|]+)\|')

    identifiers = {}

    #Find duplicate identifiers
    for rec in records:
        acc = GetRecordACC( rec.id)
        if acc in identifiers: identifiers[ acc] += 1
        else: identifiers[ acc] = 1

    outseqs = []
    written_identifiers = {}
    duplicate_warnings = 0

    custom_header_cols = set()

    for rec in records:

        custom_header_dict = ExtractCustomValuesFromFastaHeader( rec.description)
        for key in custom_header_dict.keys(): custom_header_cols.add( key)


        marked_location_in_header = -1

        cleaned_seq = re.sub( re_pattern, '', str(rec.seq))
        seq_res_num = -1
        res = re_cutnum.search( rec.description)

        header_contains_location = (res != None)
        marked_location_in_header = -1

        if header_contains_location:
            cut_1 = int( res.group( 1))
            cut_2 = int( res.group( 2))
            marked_location_in_header = cut_1


        markers = []
        seq_pos = 0 #seq starts from 1
        #Find markers
        for si in range( len( cleaned_seq)):
            aa = cleaned_seq[ si]
            if aa in AAS:
                seq_pos += 1
            elif aa != '*':
                sys.stderr.write( "FORMAT WARNING:'" + rec.description + "' contains unknown letter '%s'.\n" % ( aa))
                seq_pos += 1
            else:
                markers.append( seq_pos) #is '*'

        #[m.start() for m in re.finditer('\*', cleaned_seq)]
        #star_pos = cleaned_seq.find( '*')


        #No markers and no header location information
        if len( markers) == 0 and not header_contains_location:
            sys.stderr.write( "FORMAT ERROR:'" + rec.description + "' does not contain a valid locator. Sequence ignored.\n")
            continue

        #If no markers, use only header location information
        if len( markers) == 0: markers = [ marked_location_in_header]

        marker_index = 0
        #used_ids = {}

        acc = GetRecordACC( rec.id)

        for location in markers:

            marker_index += 1

            if aMarker == "before":
                location += -1

            #Sanity check
            if location < 0 or location > len( cleaned_seq):
                sys.stderr.write( "FORMAT ERROR:'" + rec.description + "' has an invalid locator.\n")
                continue

            new_description = rec.description.strip()
            new_acc = acc

            #Append index to identifier to keep sequences' indentifiers unique
            #Use _Pos identifier for all sequences
            #if True: #(len( markers) > 1 or identifiers[ acc] > 1) and new_acc.find("_Pos") < 0:
                #new_description = re.sub( re_id, r'\1\2_Pos%i|' % location, new_description)
            new_acc = new_acc + ("_Pos%i" % location)
            orig_acc = GetRecordACC( rec.id, aRemoveInvalidChars=False)
            #new_description = new_description.replace( orig_acc, new_acc)
            new_description = new_acc
            new_description = new_description.replace( ",", "_") #No commas allowed in description to make PDB remark parsing possible

            #DEBUG
            #print "ACC:", acc
            #print "NEWACC:", new_acc
            #print "NEW_DESC:", new_description
            #sys.exit( 0)


            #Has this sequences with same position already been written?
            if new_acc in written_identifiers:
                if duplicate_warnings == 10:
                    sys.stderr.write( "FORMAT WARNING: Only first 10 duplicates warned.\n")
                else:
                    sys.stderr.write( "FORMAT WARNING: Record '%s' position %i is a duplicate.\n" % (rec.description[:30], location))
                    sys.stderr.write( "                Duplicate skipped.\n")
                duplicate_warnings += 1
                continue #Skip it

            written_identifiers[ new_acc] = True
            cleaner_seq = cleaned_seq.replace( "*", "")

            #Insert header locator
            if aHeader and not header_contains_location:
                seq_before = cleaner_seq[max(0,location-6):location]
                seq_before = ("-"*(6-len( seq_before)))+seq_before #pad with "-"
                seq_after = cleaner_seq[location:location+6]
                seq_after = seq_after+("-"*(6-len( seq_after)))  #pad with "-"

                new_description = new_description + ("|POI:%i-%i|POISEQ:%s/%s" % (location, location+1, seq_before, seq_after))

            #Remove header locator
            if not aHeader and header_contains_location:
                #print "old: '%s'" % rec.description
                new_description = new_description[:new_description.find("|POI:")]
                #print "new: '%s'" % rec.description

            first_col = new_description.split("|")[ 0]
            if not first_col.islower() or not first_col.isalpha():
                new_description = "xx|" + new_description


            #Strip
            new_description = "|".join( map( str.strip, new_description.split( "|")))

            #Insert star
            if aStar:
                cleaned_seq = cleaner_seq[:location] + "*" + cleaner_seq[location:]

            #Remove star
            if not aStar:
                cleaned_seq = cleaner_seq



            #print "DESC: '%s'" % new_description
            #print "ID: '%s'" % new_description
            #print "NAME: '%s'" % new_description

            new_record = rec[:]
            new_record.id = new_description
            new_record.name = ""
            new_record.description = ""
            new_record.seq = Seq( cleaned_seq, IUPAC.protein)
            outseqs.append( new_record)
            #rec.seq = Seq( re.sub("(.{64})", "\\1\n", cleaned_seq, 0, re.DOTALL), IUPAC.protein)

    #print records[ 0].description
    #print records[ 0].seq

    output_handle = open( aOutputFile, "w")
    SeqIO.write( outseqs, output_handle, "fasta")
    output_handle.close()

    sys.stderr.write( "FORMAT INFO: Wrote %i sequence records to file '%s'.\n" % (len(outseqs), aOutputFile))

    if duplicate_warnings:
        sys.stderr.write( "FORMAT INFO: Found %i duplicate entries.\n" % (duplicate_warnings))

    reserved_cols = ReservedHeaderCustomCols( list( custom_header_cols))
    if len( reserved_cols):
        sys.stderr.write( "FORMAT WARNING: Custom column names '%s' are conflicting with reserved column names.\n" % ( ",".join(reserved_cols) ))

    return True
