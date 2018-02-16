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

curdir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(  "%s/../asa" % curdir)
from pdbatoms import ReadPDB

sys.path.append(  "%s/../utils" % curdir)
from struct_methods import *
from fileops import *
from progress_bar import *

#Return (AVG, MIN, MAX)
def GetBColValues( aFilename, aResindex, aChain, aPositions):

    if len( aChain) == 2 and aChain[ 0] == aChain[ 1]: aChain = aChain[ 0].lower()
    pdb_info = ReadPDB( aFilename, [ aChain], firstmodelonly=True, aQuiet=True, aConvertPTMs=True)
    n_residues = len( pdb_info["residues"])

    if n_residues == 0:
        sys.stderr.write("COLLECT ERROR: No residues in file %s for chain %s.\n" % (aFilename, aChain))
        return []

    if aPositions == 0: #special case for neighboring aas
        first_res_ind = aResindex
        last_res_ind = aResindex+1
    else:
        first_res_ind = aResindex - (aPositions-1)
        last_res_ind = aResindex

    first_res_ind = max( 0, first_res_ind)
    last_res_ind = min( n_residues-1, last_res_ind)

    bcol_vals = []

    for resid in pdb_info["residues"]:
        if resid.resindex >=  first_res_ind and resid.resindex <= last_res_ind: # and resid.Chain() == aChain:
            bcol_vals.append( resid.BColAvg())

    if len( bcol_vals) == 0:
        sys.stderr.write("COLLECT WARNING: No bcol values found in file %s for residue index %i:%s.\n" % (aFilename, aResindex, aChain))

    return bcol_vals

    #if not len( bcol_vals): return (0, 0.0,0.0,0.0) #Residues not found
    #return (len( bcol_vals), sum( bcol_vals) / len( bcol_vals), min( bcol_vals), max( bcol_vals))



def ParseFastaHeader( header, outputdict, gene_dict, log_file="" ):

    cols = header.split("|")

    for col in cols:
        if col == "sp" or col == "tr" or col == "gn":
            outputdict["DATABASE"] = col
        elif col.find("POI:") == 0:
            outputdict["SEQPOISITE"] = col[ len( "POI:"):]
        elif col.find("POISEQ:") == 0:
            outputdict["SEQPOISEQ"] = col[ len( "POISEQ:"):]
        elif col.find("CLV:") == 0:
            outputdict["CLV"] = col[ len( "CLV:"):]
        elif col.find("NETSURF:") == 0:
            outputdict["NETSURF"] = col[ len( "NETSURF:"):]

    entry_acc = GetRecordACC( header)
    outputdict["SEQ_ID"] = entry_acc
    gene = ""
    found_in_gene_dict = False

    if entry_acc in gene_dict:
        gene =  gene_dict[ entry_acc]
        found_in_gene_dict = True

    outputdict["GENE"] = gene

    #DEBUG
    #if header.find( "O75390") >= 0:
    #    print "\n\n\n\n\n\n"
    #    print "header:", header

    outputdict["DESC"] = ParseDescriptionFromFastaHeader( header)

    return found_in_gene_dict



#Parse REMARK 69 info
def GetResultFileRemarkInfo( aResult_file):

    #REMARK 69  RUNINFO: DATASERIES:input_sequences
    #REMARK 69  SEQINFO: ENTRY:A4D111_Pos77_3
    #REMARK 69  SEQINFO: DESC: gn|A4D111_Pos77|OS=Homo sapiens|POI_77-78|POISEQ_PEEVSS/IVLTKM
    #REMARK 69  PDBINFO: TITLE:Fragment-based screening of hsp70 sheds light on the functional role of atp-binding site residues
    #REMARK 69  PDBINFO: HEADER:Chaperone, COMPOUND:Heat shock cognate 71 kda protein, PUBLISHED: 22-SEP-15
    #REMARK 69  PDBINFO: KEYWORDS: Heat shock protein; hsp70; hsp72; hsc70; atpase; bag1; chaperone; 2 fragment
    #REMARK 69  PDBINFO: PDBCODE:5AQL, CHAIN:A, N_RESIDUES:997, FILE:..\..\PDB\pdb5aql.ent
    #REMARK 69  PDBINFO: RESOLUTION:1.69, R-VALUE:0.17, R-FREE:0.09, TEMPF_MIN:13.65, TEMPF_MAX:106.47, METHOD:X-RAY DIFFRACTION
    #REMARK 69  PDBINFO: ORGANISM_SCIENTIFIC:HOMO SAPIENS, ORGANISM_COMMON:HUMAN, ORGANISM_TAXID:9606
    #REMARK 69  PDBSITE: PDBPOI:SER A 121 ;
    #REMARK 69  PDBSITE: PDBSEQ: EEVSS*MVLTK, PDBPOS:121
    #REMARK 69  PDBSITE: PDBMOD:sidechains, PDBAAS:1
    #REMARK 69  ORIGSEQ: POISEQ:PEEVSS*IVLTKM, POIPOS:77
    #REMARK 69  BLASTIN: RANK:3, ALIGN_SCORE:90, BITSCORE:40.80, IDENTITY:95.00
    #REMARK 69  DSSP_IN: DSSP:HHHHH*HHHHH
    #REMARK 69  ASA_ALL: ASA_STR: Buried, ASA_AVG: 9.11, ASA_VALUES: 9.11
    #REMARK 69  ASA_ISO: ASA_ISO_STR: Buried, ASA_ISO_AVG: 9.11, ASA_ISO_VALUES: 9.11
    #REMARK 69  INTERFA: PROXIMITY_CHAINS: A, DUPLICATE_CHAINS: C, INTERFACE: No
    #REMARK 69  TEMPFAC: TEMPF_AVG: 21.19, TEMPF_VALUES: 21.19

    retvals = {}


    with open( aResult_file, "r") as file:

        for line in file:

            if line.startswith( "ATOM "): break #No more remarks
            if line.find( "REMARK 69  ") != 0: continue

            #Remove remark and line identifier (i.e. "PDBSITE:")
            line = line[ len( "REMARK 69  ") + len( "RUNINFO: "):]
            cols = line.split(",")

            for col in cols:
                keyval = col.split(":", 1)

                if len( keyval) < 2:
                    sys.stderr.write("PDB RESULT FILE ERROR: Bad format for key '%s' in file '%s'.\n" % ( keyval, aResult_file))
                    assert( len( keyval) == 2)

                key = keyval[ 0].strip()
                value = keyval[ 1].strip()

                if value.find(";") >= 0:
                    value = value.split( ";") #array
                    value = filter( None, value) #Remove empty values

                retvals[ key] = value
                #print "KEY:'%s' VAL:'%s'" % (str(keyval[ 0].strip()), str(":".join( keyval[ 1:]).strip()))
                try:
                    if key == "PDBPOI" or key == "PDBCLV":
                        index = 0
                        for poi_res in value:
                            index += 1
                            #retvals["POIRESNAME_%i" % index]  = value[ 0:3].strip()
                            #retvals["POICHAIN_%i"    % index]  = value[ 4]
                            #retvals["POIRESNUM_%i"   % index]  = value[ 5:10].strip()
                            retvals["POIRESIDUE_%i"   % index]  = poi_res
                except:
                    sys.stderr.write("PDB RESULT FILE ERROR: Bad format for key '%s' in file '%s'.\n" % ( key, result_file))
                    sys.exit(-33)

    return retvals


def GetBColVal( pdb_file, chain, resnum):

    try:
        with open( pdb_file, "r") as file:

            for line in file:

                if line.find( "ATOM") == 0 or line.find( "HETATM") == 0:

                    if line[21] == chain and line[22:27].strip() == resnum:
                        return float(line[60:66])

        #If the residue is not found in the file, it is most likely due to it being in an unstructured region
        #Unstructured regions are likely to be flexible regions in the surface of the protein and thefore
        #good targets for proteases. Assigning them a good score of 100.0 (tempf, ASA) will make them stand out in the results.
        #Other causes could be that the Original PDB file contains HETATM entries or unrecognised residues. Results scoring
        #100 for both ASA and temp should be checked manually. There should not be too many.
        sys.stderr.write( "INFO: Residue %s, chain:'%s', not found in file '%s'.\n" % ( resnum, chain, pdb_file))
        sys.stderr.write( "      Assigning default b-col value of 100.0 (most likely it is in an unstructured region).\n")
        return 100.0

    except: pass

    sys.stderr.write( "ERROR: Failed to read residue %i (chain %s) in file '%s'. Temp factor (b-col) could not be read.\n" % ( resnum, chain, pdb_file ))
    return -1.0


def CreateGeneDictionaryForEntries( seq_file, log_file):

    gene_dict = {}

    #Read sequence file for getting GENE names
    if seq_file and len( seq_file) and os.path.isfile( seq_file):

        LogIt( "INFO: Reading 'GN=' params from fasta headers in file '%s' for gene names" % seq_file, log_file, 1 )

        #seq_records = list( SeqIO.parse( seq_file, "fasta"))
        gene_dict = {}
        n_headers = 0

        try:
            record_headers = FastaHeaderIter( seq_file)

            for header in record_headers:
                n_headers += 1
                gene_dict[ GetRecordACC( header)] = GetGeneFromFastaHeader( header)

            if n_headers == 0:
                LogIt( "ERROR: No valid sequence records found in file: %s.\nStopping." % seq_file, log_file, 2)
                sys.exit( -3)

            if len( gene_dict):
                LogIt( "INFO: Found %i gene names for %i entries. (%.1f%%)" % ( len( gene_dict), n_headers, float(len( gene_dict))/n_headers*100.0), log_file, 1 )
            elif len( gene_dict) == 0:
                LogIt( "ERROR: No gene names were found in file %s. (>X|X|X|X GN=GENE_NAME)\nStopping." % seq_file, log_file, 2)
                sys.exit( -4)

        except IOError as e:
            LogIt( "ERROR: File error({0}): {1}\n".format(e.errno, e.strerror), log_file, 2)
            sys.exit( -1)
    elif not seq_file or len( seq_file) == 0:
        LogIt( "No sequence file for reading GENE names was specified.", log_file, 1 )
    else:
        LogIt( "Sequence file '%s' for reading GENE names was not found." % seq_file, log_file, 1 )

    return gene_dict


def ListIt( aVar):
    if type( aVar) != list: return [aVar]
    return aVar


#Goes through result pdb files and collects the data to a single (tab-delimited) file
def CollectResults( result_folder, output_file, seq_file="", log_file=None, aCustomColumns=[]):

    #Read gene names from seq_file (human_proteome)
    #File needs to be a fasta file that has "GN=" definitions in seq descriptors (">")
    gene_dict = CreateGeneDictionaryForEntries( seq_file, log_file)
    result_folder = InsertFolderSeparator( result_folder) + "poi/"

    #Find result files
    result_files = GetDirFiles( result_folder + "*.poisite.pdb")

    if len( result_files) == 0:
        LogIt( "COLLECT ERROR: No result files found in folder '%s'." % result_folder, log_file, 2 )
        sys.exit( -32)
    else:
        LogIt( "COLLECT INFO: Found %i result files to collect data from." % len( result_files), log_file, 1 )

    #RESULT FILE COLUMNS
    col_order = ["DATASERIES","ENTRY","SEQ_ID","DESC","GENE","TEMPF_AVG","TEMPF_VALUES", "TEMPF_POI_MIN", "TEMPF_POI_MAX", \
                 "ISOL_ASA", "ISOL_ASA_AVG","ISOL_ASA_POI_MIN","ISOL_ASA_POI_MAX", \
                 "BIOL_ASA", "BIOL_ASA_AVG","BIOL_ASA_POI_MIN","BIOL_ASA_POI_MAX", \
                 "ASYM_ASA", "ASYM_ASA_AVG","ASYM_ASA_POI_MIN","ASYM_ASA_POI_MAX", \
                 "POI_CHAIN", "PROXIMITY_CHAINS", "DUPLICATE_CHAINS", "INTERFACE", "INTERFACE_ASA_DELTA", "DSSP", \
                 "PDBCODE", "PDBAAS", "PDBPTM", "TEMPF_FILE_MIN", "TEMPF_FILE_MAX", "PDBMOD", "ALL_CHAINS", "TITLE", "COMPOUND", "KEYWORDS", "DATABASE", \
                 "POIPOS","POISEQ","PDBSEQ","ALIGN_SCORE", \
                 "RANK","BITSCORE","N_RESIDUES","IDENTITY", "RESOLUTION", "R-VALUE", "R-FREE", "METHOD", "ORGANISM_COMMON", \
                 "ORGANISM_SCIENTIFIC", "ORGANISM_TAXID", "PUBLISHED"]

    #Add custom columns from fasta header
    aCustomColumns = list( set( aCustomColumns)) #Make sure values are unique
    for custom_col in reversed( aCustomColumns):
        if custom_col in col_order:
            LogIt( "COLUMN ERROR: Cannot use custom column name '%s' because it is already in use." % custom_col, log_file, 2 )
            aCustomColumns.remove( custom_col)

    if len( aCustomColumns):
        print "INFO: Adding custom column%s '%s' to results." % ("s" if len( aCustomColumns) > 1 else "", ", ".join( aCustomColumns))
        col_order.extend( aCustomColumns)


    gene_found = 0
    n_genes = 0

    #Write result file
    print "COLLECT PROGRESS: Writing collected results file '%s'..." % output_file

    header_written = False

    ProgressBar( aProgress=0, aTotal=len( result_files), aIntervalSeconds=1.0)

    with open( output_file, "w") as f:


        for poi_file in result_files:

            file_data = GetResultFileRemarkInfo( poi_file)
            assert( "PDBPOS" in file_data and "PDBAAS" in file_data)

            try:
                resindex = int( file_data["PDBPOS"]) - 1 #position to index
                positions = file_data["PDBAAS"] #How many positions upstream to use (0 for neighboring)
                if positions == "between": positions = 0
                else: positions = int( positions)

                tempf_values = map( float, file_data["TEMPF_VALUES"])
                tempf_values = filter( lambda a: a >= 0.0, tempf_values)
                if not len( tempf_values): tempf_values.append( -1.0) # If no values remain after filter, insert negative value
                file_data["TEMPF_POI_MIN"] = "%.2f" % min( tempf_values)
                file_data["TEMPF_POI_MAX"] = "%.2f" % max( tempf_values)

                asa_values = map( float, file_data["ASYM_ASA_VALUES"])
                asa_values = filter( lambda a: a >= 0.0, asa_values)
                if not len( asa_values): asa_values.append( -1.0)
                file_data["ASYM_ASA_POI_MIN"] = "%.2f" % min( asa_values)
                file_data["ASYM_ASA_POI_MAX"] = "%.2f" % max( asa_values)

                asa_isol_values = map( float, file_data["ISOL_ASA_VALUES"])
                asa_isol_values = filter( lambda a: a >= 0.0, asa_isol_values)
                if not len( asa_isol_values): asa_isol_values.append( -1.0)
                file_data["ISOL_ASA_POI_MIN"] = "%.2f" % min( asa_isol_values)
                file_data["ISOL_ASA_POI_MAX"] = "%.2f" % max( asa_isol_values)

                asa_biol_values = map( float, file_data["BIOL_ASA_VALUES"])
                asa_biol_values = filter( lambda a: a >= 0.0, asa_biol_values)
                if not len( asa_biol_values): asa_biol_values.append( -1.0)
                file_data["BIOL_ASA_POI_MIN"] = "%.2f" % min( asa_biol_values)
                file_data["BIOL_ASA_POI_MAX"] = "%.2f" % max( asa_biol_values)

                positions = file_data["PDBAAS"]
                if positions.lower() == "between": positions = 2
                else: positions = int( positions)
                for p in range( 1, positions+1):
                    key = "POIRESIDUE_%i" % p
                    if key not in file_data: file_data[ key] = ""


                gene_found += ParseFastaHeader( file_data["DESC"], file_data, gene_dict, log_file)
                n_genes += 1

            except ValueError as ex:
                LogIt( "COLLECT ERROR: Corrupt remark header in file '%s'. (%s)\n" % (poi_file, str( ex)), log_file, 2)
                raise #re-raise

            #Add empty values if value not found in tempf file
            for co in col_order:
                if co not in file_data.keys(): file_data[ co] = ""

            #Check first file for all keys
            #Add data not in col_order
            if not header_written:
                unordered = list( set( file_data) - set( col_order))

                for uo in unordered:
                    if uo == "FILE" and os.name != "nt": continue #Do not add file col on online version of server
                    col_order.append( uo)
                #Write column headers
                f.write( "\t".join( col_order))
                f.write( "\n")
                header_written = True

            line = []

            for co in col_order:
                value = file_data[ co]
                if type( value) == list: line.append( ",".join([str( x) for x in value]))
                else: line.append( str( value))

            #Write result
            f.write( "\t".join( line))
            f.write( "\n")

            ProgressBar( aProgress=n_genes)

            #if n_genes % 1000 == 0:
            #    print "Progress: %i/%i (%.1f%%)" % (n_genes, len( result_files), float(n_genes)/len( result_files)*100.0)


        ProgressBar( aProgress=n_genes, aTotal=len( result_files), aFinalize=True)

        if n_genes == 0:
            sys.stderr.write( "COLLECT ERROR: No result genes found.\n")
            sys.exit( -33)

        LogIt( "\n\nGene information found reliably for %i/%i entries. (%.1f%%)" % (gene_found, n_genes, float(gene_found)/n_genes*100.0), log_file, 1)
        return col_order

