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
import gc
from datetime import datetime, timedelta

from main_record import MainRecord
from process_dssp import *

curdir = os.path.dirname(os.path.abspath(__file__))

sys.path.append( "%s/../utils" % curdir)
from struct_methods import *
from fileops import *
import online_pdb
from obj_cache import ObjCache
from fasta_formatter import ExtractCustomValuesFromFastaHeader

sys.path.append( "%s/../asa" % curdir)
from pdbatoms import *
from asa_uta import CalcASAForResidues #,runASAFiles


#Writes POISITE PDB files for all cleavage sites
#Information about the cleavage site analysis is inserted in the beginning of generated files as "REMARK 69" lines
#The cleavage site in the PDB file is located here
#DSSP prediction is fetched or calculated here
def WriteResultFiles( main_records, seq_file, positions, sidechain_only, output_folder, dssp_folder="", threads=-1, reset=False,  \
                      poisite_log=None, precalc_asa_folder=None, aCustomColsFile="", aWriteUsedBiomolecules=True, aDSSPExecutable=""):

    n_main_records = 0 if not main_records else len( main_records)

    if n_main_records == 0:
        sys.stderr.write( "WRITE ERROR: No records to process.\n")
        return -54 #Error

    SetDSSP_Executable( aDSSPExecutable)

    #Set folders
    seq_path, seq_filename, seq_name, seq_suffix = SplitFilePath( seq_file)
    output_folder = InsertFolderSeparator( output_folder)

    poi_folder = output_folder + "poi/"
    if not os.path.exists( poi_folder):
        try:
            os.makedirs( poi_folder)
        except Exception as ex:
            LogIt( "WRITE ERROR: Could not create folder '%s' for result files. (%s)\n" % (poi_folder, str( ex)), poisite_log, 2)
            sys.exit( -34)


    #Create data strcuture for processing the entries
    pdb_dict = {}
    n_entries = 0

    #Create a dict of all pdb files to be processed
    #Process everything sorted by required pdb file to minimize disk I/O
    for m_rec in main_records:

        blast_entries = m_rec.BlastEntries()
        #Count entries that need to be processed
        n_entries += len( blast_entries)

        for be in blast_entries:
            pdb_id = be["pdb_id"]
            if pdb_id not in pdb_dict: pdb_dict[ pdb_id] = []
            elem = (m_rec, be)
            pdb_dict[ pdb_id].append( elem)


    entry_progress = 0
    n_errors = 0


    #DSSP - secondary structure prediction
    dssp_log_filename = output_folder + "dssp_log.txt"
    dssp_success = 0
    dssp_fail = 0
    dssp_skipped = 0
    dssp_unknown = 0
    dssp_unknown_str = "Unknown"
    dssp_gapped_structures = 0
    missing_dssp_files = {}
    gapped_dssp_postitions = {}
    processing_started = False
    n_biomol = 0
    n_no_biomol = 0

    pdb_cache = ObjCache( 1)
    dssp_cache = ObjCache( 1)
    biomol_cache = ObjCache( 5) #Biomols can be very large in memory

    custom_cols = []
    n_processed_structures = -1
    n_structures_to_process = len( pdb_dict)
    LogIt( "WRITE INFO: Processing %i pdb files for %i POI entries." % (n_structures_to_process, n_entries), poisite_log, 1)

    #Process everything by required pdb file to minimize disk I/O
    for pdb_id, entry_list in pdb_dict.iteritems():
    #for pdb_id in pdb_dict.keys():

    	#entry_list = pdb_dict[ pdb_id]
        biomol_cache.Reset()
        pdb_cache.Reset()
        dssp_cache.Reset()

        #if n_processed_structures % 100 == 0: gc.collect() #Make sure we have enough memory

        #Loop for all entries that use the same pdb file
        for m_rec, b_entry in entry_list:

            assert( b_entry["pdb_id"] == pdb_id)

            if n_errors > 100:
                sys.stderr.write( "ERROR: Encountered too many errors, exitting.\n")
                return -1 #EXIT

            entry_progress += 1

            #Compose filenames
            entry_id = m_rec.EntryIdentifier( b_entry["rank"])
            pdb_filename = ValidFilename( entry_id)

            poisite_pdb_file = SystemSpecificFolderSeparators( poi_folder + pdb_filename + ".poisite.pdb")

            if not reset:
                #Already processed
                if os.path.isfile( poisite_pdb_file) and int( os.stat( poisite_pdb_file).st_size) > 0:

                    #Result file found with > 0 filesize, no need to process further
                    LogIt(  "%s: Files already processed. (%i/%i)" % ( entry_id, entry_progress, n_entries), poisite_log, -1 )
                    dssp_skipped += 1
                    continue

            pdb = None
            pdb_filepath = b_entry["pdb_filepath"]
            #pdb_id = b_entry["pdb_id"]
            pdb_chain = RealChain( b_entry["chain"])
            pdb_key = pdb_id
            pdb_is_cached = pdb_cache.IsInCache( pdb_key)
            processing_started = True

            if pdb_is_cached:
                #Use cached file
                pdb = pdb_cache.Retrieve( pdb_key)
                Residue.ResetBCols( pdb["residues"]) #Remove "coloring" based on previous POI
            else:

                #Read original PDB File
                try:
                    #Make sure pdb file exists
                    if not os.path.isfile( pdb_filepath) or int( os.stat( pdb_filepath).st_size) <= 0:
                        online_pdb.FetchPDBFile( pdb_id, output_folder=GetFolder( pdb_filepath), force_fetch=True, aAsync=False)

                    with open( pdb_filepath, "r" ) as pdb_fh:
                        #Read all chains, and first model only
                        #PTMs have been converted during alignment so they need to be converted here to make
                        #position (index) within the chain of residue match
                        pdb = ReadPDB( pdb_fh, chains=[], firstmodelonly=True, aConvertPTMs=True, aWarnUnknown=False)

                    pdb_cache.Add( pdb_key, pdb)

                except IOError as e:
                    LogIt(  "ERROR: %s: Unable to read pdb file: '%s' Error: %s\n" % (entry_id, b_entry["pdb_filepath"], e.strerror), poisite_log, 2 )
                    n_errors += 1
                    continue
                except Exception as ex:
                    LogIt( "ERROR: %s: Error reading PDB file: '%s'.\n" % (entry_id, b_entry["pdb_filepath"]), poisite_log, 2)
                    LogIt( "       An exception of type {0} occured. Arguments: {1!r}".format( type( ex).__name__, ex.args), poisite_log, 2 )
                    n_errors += 1
                    continue

                if len( pdb["residues"]) == 0 or pdb["stats"]["amino_acids"] == 0:
                    LogIt( "WRITE ERROR: PDB file '%s' contained no amino acid residues." % (os.path.basename( pdb_filepath)), poisite_log, 2 )
                    n_errors += 1
                    continue

            stats = pdb["stats"]
            #Report processing
            sys.stdout.write( "[%s]%s %s|%s : %s\n" % (("+" if pdb_is_cached else "-"), ("b " if stats[ "n_biomolecules"] > 0 else "  "), pdb_id, pdb_chain, entry_id))

            header_vals = ExtractCustomValuesFromFastaHeader( m_rec.Description())

            for hv in header_vals:
                if hv not in custom_cols:
                    custom_cols.append( hv)
                    WriteListToFile( custom_cols, aCustomColsFile, aAppend=True)

            b_entry = dict( b_entry.items() + header_vals.items()) #Combine 2 dictionaries
            b_entry_identifier = entry_id

            remarks = []
            poi_remark_index = -1


            #Replace all commas with ';' to keep remarks readable
            for key in stats:
                value = stats[ key]
                if type( value) == str: stats[ key] = value.replace(',','_')

            try:
                                #    keep length
                                #|                 |
                remarks.append( "REMARK 69  RUNINFO: DATASERIES:%s" % Remark( seq_name)) #0 remove chars that could complicate parsing later
                remarks.append( "REMARK 69  SEQINFO: ENTRY:%s"      % Remark( b_entry_identifier))     #1
                remarks.append( "REMARK 69  SEQINFO: DESC:%s"       % FormatFastaDescription( m_rec.Description()))
                remarks.append( "REMARK 69  ORIGSEQ: POISEQ:%12s, POIPOS:%s" % ( Remark( b_entry["seq_poiseq"]), Remark( b_entry["seq_poipos"])))
                remarks.append( "REMARK 69  PDBINFO: TITLE:%s"      % Remark( stats["title"].capitalize()))
                remarks.append( "REMARK 69  PDBINFO: HEADER:%s, COMPOUND:%s, PUBLISHED: %s" % ( Remark( stats["header"].capitalize()), Remark( stats["compound"].capitalize()), Remark( stats["published"])))
                remarks.append( "REMARK 69  PDBINFO: KEYWORDS: %s"  % Remark( stats["keywords"].capitalize()))
                remarks.append( "REMARK 69  PDBINFO: PDBCODE:%s, ALL_CHAINS:%s, N_RESIDUES:%i, FILE:%s" % (pdb_id, Remark( stats["chains"]), stats["n_residues"], b_entry["pdb_filepath"]))
                remarks.append( "REMARK 69  PDBINFO: RESOLUTION:%.2f, R-VALUE:%.2f, R-FREE:%.2f, TEMPF_FILE_MIN:%.2f, TEMPF_FILE_MAX:%.2f, METHOD:%s" % (stats["resolution"],stats["R-Value"],stats["R-Free"],stats["TempF-min"],stats["TempF-max"],Remark( stats["method"]) ))
                remarks.append( "REMARK 69  PDBINFO: ORGANISM_SCIENTIFIC:%s, ORGANISM_COMMON:%s, ORGANISM_TAXID:%i" % (stats["organism_scientific"], stats["organism_common"], stats["organism_taxid"]))
                remarks.append( "REMARK 69  BLASTIN: RANK:%i, ALIGN_SCORE:%i, BITSCORE:%.2f, IDENTITY:%.2f" % (b_entry["rank"], b_entry["alignment_score"], b_entry["bitscore"], b_entry["identity"] ))
                remarks.append( "REMARK 69  PDBSITE: PDBSEQ:%12s, PDBPOS:%s" % ( b_entry["pdb_poiseq"], b_entry["pdb_poipos"] )) #7
                remarks.append( "REMARK 69  PDBSITE: PDBMOD:%s, PDBAAS:%s, PDBPTM:" % ( "sidechains" if sidechain_only else "residue", "between" if positions == 0 else str( positions))) #7
                ptm_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  PDBSITE: POI_CHAIN: %s, PDBPOI:" % b_entry["chain"]) #filled later
                poi_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  DSSP_IN: DSSP:") #filled later
                dssp_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  ASAISOL: ") #filled later
                asa_iso_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  ASABIOL: ") #filled later
                asa_bio_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  ASAASYM: ") #filled later
                asa_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  INTERFA: ") #filled later
                interface_remark_index = len( remarks)-1
                remarks.append( "REMARK 69  POITEMP: ") #filled later
                tempf_remark_index = len( remarks)-1

                #Custom remarks
                if len( header_vals):
                    customs =   "REMARK 69  CUSTOMS: "
                    first = True
                    for hv in header_vals.keys():
                        if not first: customs += ", "
                        customs += "%s:%s" % ( hv, header_vals[ hv])
                        first = False

                    #remarks.append( "REMARK 69  CUSTOMV: NETSURF:%s, CLV:%s" % (b_entry["NETSURF"], b_entry["CLV"]))
                    remarks.append( customs)

            except Exception as ex:
                raise
                sys.stderr.write( "ERROR: %s: Unable to Format PDB file header.\n" % ( b_entry_identifier))
                sys.stderr.write( "       An exception of type '{0}' occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
                printDictionary( b_entry)

                n_errors += 1
                return -2

            #Find POI-site residues
            pos = 0
            POI_residues = []

            #Find POI residues
            for residue in pdb["residues"]:

                #If pdb data structure is from cache b-col values have been modified -> restore them
                if residue.Chain() != pdb_chain: continue
                elif residue.AtomType() != "A": continue
                elif residue.IsAlternateLocation(): continue

                #POIPOS determines the Nth amino acid, skip everything else.
                #When sequence is extracted from the PDB, everything but known aa codes
                #are marked as '-'. In the alignment file POIPOS is calculated by non-hyphen chars.
                resname = residue.ResName()

                #if len( resname) != 3 or AAcode_3_to_1( resname) == '-': continue

                pos += 1
                #if pos < 15: print ("RES %i:" % pos), residue.ResName(), residue.ResNum()

                if pos == b_entry["pdb_poipos"]:
                    #N-Terminal residue at POI found
                    POI_residues = FindResidues( pdb["residues"], residue.ResNum(), residue.Chain(), aAdjacentPositions=(1 if positions == 0 else -1*(positions-1)))

                    pdb_poi_res = ""
                    pdb_ptm_types = ""
                    for poi_res in POI_residues:
                        poi_resnum = poi_res.ResNum()
                        if poi_resnum.isdigit(): poi_resnum += " " #Fill place for insertion code
                        pdb_poi_res += "%-3s %s%5s%s" % ( poi_res.ResName(), residue.Chain(), poi_resnum, ("" if len( pdb_poi_res) == 0 else ";")) #separator
                        pdb_ptm_types += "%s%s" % ( GetShortModificationType( stats, residue), ("" if len( pdb_ptm_types) == 0 else ";")) #separator
                    remarks[ poi_remark_index] += pdb_poi_res + ";" #Fill PDBSITE: PDBPOI:
                    remarks[ ptm_remark_index] += pdb_ptm_types + ";"

                    break #POI residues found

            if len( POI_residues) == 0:
                sys.stderr.write( "ERROR: POI position %i/%i in chain %s was not found for entry '%s' in file %s\n" % (b_entry["pdb_poipos"], pos, b_entry["chain"], entry_id, pdb_filepath))
                n_errors += 1
                return -54
                #continue

            #Get duplicate (dimer) chains in file
            duplicate_chains = SimpleFindDuplicateChains( pdb, aChains=[], aChainsToReport=[pdb_chain], aQuiet=False )
            duplicate_chains = sorted( duplicate_chains[ pdb_chain]) #Make it a list

            #Get relevant chains by spatial proximity of POI
            relevant_chains = set( pdb_chain)
            in_range_residues = set()
            for poi_res in POI_residues:
                res_list = ResiduesInRange( poi_res, pdb["residues"], aDistance=10.0 )
                for r in res_list: in_range_residues.add( r)
            in_range_residues = list( in_range_residues)

            for irr in in_range_residues: relevant_chains.add( irr.Chain())
            relevant_chains = list( relevant_chains)

            #######################
            #Get Tempf
            tempf_values = []
            for poi_res in POI_residues: tempf_values.append( poi_res.BColAvg())

            remarks[ tempf_remark_index] += "TEMPF_AVG: %.2f, " % (-1.0 if len( tempf_values) == 0 else sum( tempf_values) / float( len( tempf_values)))
            remarks[ tempf_remark_index] += "TEMPF_VALUES: %s;, " % (";".join( ['{:.2f}'.format( x) for x in tempf_values]))
            remarks[ tempf_remark_index] += "TEMPF_POI_MIN: %.2f, " % (-1.0 if len( tempf_values) == 0 else min( tempf_values))
            remarks[ tempf_remark_index] += "TEMPF_POI_MAX: %.2f"   % (-1.0 if len( tempf_values) == 0 else max( tempf_values))

            #DEBUG
            #print "POI_RESIDUES:"
            #for r in POI_residues:
            #    print str( r)
            #sys.exit( 0)

            #######################
            #Calculate ASA
            asa_asym_values = CalcASAForResidues( POI_residues, pdb, relevant_chains, aCalcMode=("S" if sidechain_only else "R"), aQuiet=True)
            #asa_asym_values_f = filter( lambda a: a >= 0.0, asa_asym_values ) #Remove negative values that mean there were no atoms that were included (e.g. sidechain mode for Gly)

            asym_asa_decision, asym_asa_avg = ASADecision( FastaSequence( POI_residues), asa_asym_values, sidechain_only )
            #asym_asa_avg = -1.0 if len( asa_asym_values_f) == 0 else sum( asa_asym_values_f) / float( len( asa_asym_values_f))
            remarks[ asa_remark_index] += "ASYM_ASA: %s, " % asym_asa_decision
            remarks[ asa_remark_index] += "ASYM_ASA_AVG: %.2f, " % asym_asa_avg
            remarks[ asa_remark_index] += "ASYM_ASA_VALUES: %s;" % (";".join( ['{:.2f}'.format( x) for x in asa_asym_values])) #Add ; to make it an array even if one value

            #######################
            #POI chain in isolation
            asa_iso_values = CalcASAForResidues( POI_residues, pdb, POI_residues[ 0].Chain(), aCalcMode=("S" if sidechain_only else "R"), aQuiet=True) #isolated chain
            #asa_iso_values_f = filter( lambda a: a >= 0.0, asa_iso_values )

            asa_iso_decision, isol_asa_avg = ASADecision( FastaSequence( POI_residues), asa_iso_values, sidechain_only )
            #isol_asa_avg = -1.0 if len( asa_iso_values_f) == 0 else sum( asa_iso_values_f) / float( len( asa_iso_values_f))
            remarks[ asa_iso_remark_index] += "ISOL_ASA: %s, " % asa_iso_decision
            remarks[ asa_iso_remark_index] += "ISOL_ASA_AVG: %.2f, " % isol_asa_avg
            remarks[ asa_iso_remark_index] += "ISOL_ASA_VALUES: %s;" % (";".join( ['{:.2f}'.format( x) for x in asa_iso_values])) #Add ; to make it an array even if one value

            ######################
            #Calc ASA from POI chain in biological assembly (if specified in PDB file)
            asa_bio_decision = "None"
            asa_bio_values = []
            #asa_bio_values_f = [] #filtered
            biol_asa_avg = -1.0
            bio_pdb = None

            if stats[ "n_biomolecules"] > 0:

                #Construct biological assemblies
                bio_found = False
                assemblies = ReadBiologicalAssemblies( pdb["remarks"])
                for assembly in assemblies:
                    if pdb_chain not in assembly["chains"]: continue

                    biomol_id = pdb_id + str( assembly["num"])

                    if biomol_cache.IsInCache( biomol_id):
                        bio_pdb = biomol_cache.Retrieve( biomol_id)
                    else:
                        bio_pdb = NewPdb()
                        #sys.stderr.write( "Creating biological assembly %i... " % assembly["num"])
                        CreateBiologialAssembly( assembly, aPdbSource=pdb, aPdbTarget=bio_pdb, aQuiet=True, aCenterResidue=POI_residues[ -1], aRestrictAboveAtoms=100000, aDistance=150.0)
                        #sys.stderr.write( "Done.\n")

                        #Add only biomols that were not restricted around a POI to cache
                        #if True or not IsRestricted( bio_pdb["remarks"]): biomol_cache.Add( biomol_id, bio_pdb)
                        if not IsRestricted( bio_pdb["remarks"]): biomol_cache.Add( biomol_id, bio_pdb)

                    residue = POI_residues[ 0]
                    biol_POI_residues = FindResidues( bio_pdb["residues"], residue.ResNum(), residue.Chain(), aAdjacentPositions=(1 if positions == 0 else -1*(positions-1)))

                    asa_bio_values = CalcASAForResidues( biol_POI_residues, bio_pdb, aChains=[], aCalcMode=("S" if sidechain_only else "R"), aQuiet=True)
                    #asa_bio_values_f = filter( lambda a: a >= 0.0, asa_bio_values )
                    asa_bio_decision, biol_asa_avg = ASADecision( FastaSequence( biol_POI_residues), asa_bio_values, sidechain_only )
                    if aWriteUsedBiomolecules:
                        SetBColForResultFileResidues( bio_pdb["residues"], biol_POI_residues, duplicate_chains)
                    bio_found = True
                    n_biomol += 1
                    break #Create and use only first biomol chain is a part of (if in multiple)

                if not bio_found: LogIt( "INFO: Chain '%s' not found in biological assemblies.\n" % pdb_chain, poisite_log, 1)
            else: n_no_biomol += 1

            #biol_asa_avg = -1.0 if len( asa_bio_values_f) == 0 else sum( asa_bio_values_f) / float( len( asa_bio_values_f))
            remarks[ asa_bio_remark_index] += "BIOL_ASA: %s, " % asa_bio_decision
            remarks[ asa_bio_remark_index] += "BIOL_ASA_AVG: %.2f, " % biol_asa_avg
            remarks[ asa_bio_remark_index] += "BIOL_ASA_VALUES: %s;" % (";".join( ['{:.2f}'.format( x) for x in asa_bio_values])) #Add ; to make it an array even if one value

            #####################
            #Set interface
            #No, Homomer, Heteromer
            interface = "No"
            max_pairwise_delta = 0.0
            interface_asa_delta = 0.0

            if isol_asa_avg >= 0.0:
                if   biol_asa_avg >= 0.0: max_pairwise_delta = MaxPairwiseDelta( asa_iso_values, asa_bio_values) #interface_asa_delta = (isol_asa_avg - biol_asa_avg)
                elif asym_asa_avg >= 0.0: max_pairwise_delta = MaxPairwiseDelta( asa_iso_values, asa_asym_values) # (isol_asa_avg - asym_asa_avg)

            if max_pairwise_delta >= 0.1 or max_pairwise_delta <= -0.1: #0.1 percent of ASA difference
                if len( set( relevant_chains).intersection( set(duplicate_chains))) > 0: interface = "Homomer" #Dimer?
                else: interface = "Heteromer"

            #compare_asa_decision = (asa_bio_decision if asa_bio_decision != "None" else asym_asa_decision) #Use biological assembly for comparison if available
            #if (asa_iso_decision == "Surface" and (compare_asa_decision == "Buried" or compare_asa_decision == "Intermediate")) or \
            #   (asa_iso_decision == "Intermediate" and compare_asa_decision == "Buried"):
            #   if len( set( relevant_chains).intersection( set(duplicate_chains))) > 0: interface = "Homomer" #Dimer?
            #   else: interface = "Heteromer"

            proximity_chains = relevant_chains[:] #Make a copy
            proximity_chains.remove( pdb_chain)
            remarks[ interface_remark_index] += "PROXIMITY_CHAINS: %s, " % (";".join( proximity_chains)  if len( proximity_chains)  else "None;")
            remarks[ interface_remark_index] += "DUPLICATE_CHAINS: %s, " % (";".join( duplicate_chains) if len( duplicate_chains) else "None;")
            remarks[ interface_remark_index] += "INTERFACE: %s," % interface
            remarks[ interface_remark_index] += "INTERFACE_ASA_DELTA: %.2f" % interface_asa_delta

            ################
            #DSSP
            fetch_dssp = dssp_folder != None and len( dssp_folder) > 0 #Fetch files if folder specified
            #sys.stderr.write( "Processing DSSP... ")

            dssp_content = []
            only_ca_atoms = stats["single_atom_residue_ratio"] > 0.9 #Ratio of single atom residues in PDB file

            if fetch_dssp and not only_ca_atoms and pdb_id not in missing_dssp_files: #do not try refetching missing files

                if dssp_cache.IsInCache( pdb_id):
                    dssp_content = dssp_cache.Retrieve( pdb_id)
                else:
                    try:
                        dssp_content = ReadDSSP( pdb_id, dssp_folder, aPdbFolder=GetFolder( pdb_filepath), fetch_if_missing=fetch_dssp, log_file=None )
                        dssp_cache.Add( pdb_id, dssp_content)
                    except RuntimeError as ex:
                        LogIt( "ERROR: Error processing DSSP for %s: '%s'" %  (pdb_id, str( ex)), dssp_log_filename, 2)
                        if str( ex).find("Too many consecutive errors") >= 0: return -7 #Exit

            if not dssp_content or len( dssp_content) == 0 or only_ca_atoms:
                b_entry["dssp"] = dssp_unknown_str
                if fetch_dssp: missing_dssp_files[ pdb_id] = True
            else:

                dssp = ProcessDSSP( pdb_id, dssp_content, POI_residues[ 0].Chain(), POI_residues[ 0].ResNum(), max(10, positions), log_file=dssp_log_filename, aFilePath=GetDSSPpath(pdb_id, dssp_folder) )

                b_entry["dssp"] = ''.join( dssp)
                if b_entry["dssp"].find('*') >= 0:
                    dssp_success += 1
                    if b_entry["dssp"].count('?') > 0:
                        #Log structures with gaps
                        dssp_location_str = "%s:%s:%s" % ( pdb_id, POI_residues[ 0].Chain(), POI_residues[ 0].ResNum())
                        if dssp_location_str not in gapped_dssp_postitions:
                            #Log only once for each position
                            dssp_gapped_structures += 1
                            gapped_dssp_postitions[ dssp_location_str] = True
                            LogIt( "DSSP assesment for %s has gaps: %s" % (dssp_location_str.ljust( 10), b_entry["dssp"] ), dssp_log_filename, 0)
                else:
                    dssp_fail += 1
                    if b_entry["dssp"].find( 'not found') >= 0:
                        b_entry["dssp"] = dssp_unknown_str
                        #missing_dssp_files[ b_entry["pdb_id"]] = True

            if b_entry["dssp"] == dssp_unknown_str: dssp_unknown += 1
            remarks[ dssp_remark_index] += b_entry["dssp"]


            #Set chain "colors"
            #MARK POI (and duplicates) in Bcol
            SetBColForResultFileResidues( pdb["residues"], POI_residues, duplicate_chains)

            #################
            #Finalize remarks and write POI file
            orig_remarks = pdb["remarks"]
            pdb["remarks"] = remarks
            write_retval = WritePdbToFile( b_entry_identifier, pdb, poisite_pdb_file, aChains=list(set(relevant_chains)|set(duplicate_chains)))
            pdb["remarks"]= orig_remarks

            if write_retval == False: n_errors += 1

            if aWriteUsedBiomolecules and bio_pdb and len( bio_pdb["atoms"]):
                bio_pdb["remarks"] = remarks
                biomol_folder = output_folder + "biomol/"
                biomol_pdb_file = SystemSpecificFolderSeparators( biomol_folder + pdb_filename + ".biomol.pdb")
                if not os.path.exists( biomol_folder): os.makedirs( biomol_folder)
                write_retval = WritePdbToFile( b_entry_identifier, bio_pdb, biomol_pdb_file, aChains=[])


        #main records processed
        n_processed_structures += 1

        if processing_started: ReportWriteSpeed( entry_progress, n_entries)

    ##############
    #PROCESSING DONE

    #Log DSSP results


    total_processed = n_biomol+n_no_biomol
    LogIt( "\nBiological assemblies found for %i/%i (%.1f%%) POIs." % (n_biomol,total_processed,Percentage(n_biomol, total_processed)), poisite_log, 1)

    if dssp_success+dssp_fail > 0:
        n_dssp_total = dssp_success+dssp_fail+dssp_unknown
        dssp_percentage = float(dssp_success) / (n_dssp_total) * 100.0
        LogIt( "\nDSSP assesment found for %i/%i (%.1f%%) structures" % (dssp_success, n_dssp_total, dssp_percentage), dssp_log_filename, -1)
        LogIt( "\nDSSP assesment found for %i/%i (%.1f%%) structures" % (dssp_success, n_dssp_total, dssp_percentage), poisite_log, 1)
    LogIt( "\n%i files were skipped, because they had been already processed." % (dssp_skipped), dssp_log_filename, -1)
    LogIt( "\n%i files were skipped, because they had been already processed." % (dssp_skipped), poisite_log, 1)

    if len( missing_dssp_files):
        with open( dssp_log_filename, 'w') as dssp_file:
            LogIt( "Missing DSSP files (%i):\n" % len( missing_dssp_files), dssp_file, 1)
            for mf in missing_dssp_files.keys():
                dssp_file.write( mf + "\n")
    else:
        LogIt( "No missing DSSP files.", dssp_log_filename)

    return 0 #DONE

def Percentage( aVal, aTotal):
    if aTotal == 0.0: return 0.0
    return (aVal / float(aTotal)) * 100.0

def SetBColForResultFileResidues( aResidueList, aPoiResidues, aDuplicateChains=[]):

    poi_chain = aPoiResidues[ 0].Chain()
    for residue in aResidueList:
        cur_chain = residue.Chain()
        if cur_chain == poi_chain:
            if residue in aPoiResidues: residue.SetBCol( 100.0)
            else: residue.SetBCol( 0.0) #POI chain residues
        elif cur_chain in aDuplicateChains:
            if residue.IsMatchingResidue( aPoiResidues): residue.SetBCol( 95.0) #Duplicate chain matching POI-residue
            else: residue.SetBCol( 20.0 + (aDuplicateChains.index( cur_chain) % 4) * 10) #Duplicate chain residue 20 - 50
        else:
            residue.SetBCol( 70.0 + (ord( cur_chain) % 3) * 5) #Proxmity chain residues 70-80

def MaxPairwiseDelta( aList1, aList2):

    if len( aList1) != len( aList2):
        sys.stderr.write("ERROR: Lists should be of the same length: %i vs. %i\n" % (len(aList1), len(aList2)))
        assert(False)
    deltas = [0.0 if (x<0.0 or y<0.0) else (x - y) for x, y in zip(aList1, aList2)]
    return max( deltas)



def WritePdbToFile( aEntryId, aPdb, aFilename, aChains=[]):

    try:
        with open( aFilename, "w" ) as pdb_fh:
            pdb_fh.write( "\n".join( aPdb["remarks"]) + "\n")
            WriteToFile( pdb_fh, aPdb["atoms"], types=["A","D"], chains=aChains)
    except IOError as e:
        sys.stderr.write( "ERROR: %s: Unable to write pdb file: '%s' Error: %s\n" % (aEntryId, aFilename, e.strerror))
        if e.strerror.find("No space left on device") >= 0: sys.exit( -34)
        return False
    except Exception as ex:
        sys.stderr.write( "ERROR: %s: Unable to write PDB file: '%s'.\n" % (aEntryId, aFilename))
        template =        "       An exception of type {0} occured. Arguments: {1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message + "\n")
        return False

    return True



def ReportWriteSpeed( aCurrent, aTotal, aPrintLock=False):

    #Using function attributes
    if not hasattr( ReportWriteSpeed, "prev_time"):
        ReportWriteSpeed.prev_time = time.time()
        ReportWriteSpeed.prev_count = aCurrent
        ReportWriteSpeed.processing_speed_per_minute = []
        return "WRITE INFO: Calculating processing speed..."

    now = time.time()
    seconds_passed = int( now - ReportWriteSpeed.prev_time)
    if seconds_passed < 15: return #report interval

    n_processed = aCurrent - ReportWriteSpeed.prev_count
    entries_per_minute = round( n_processed * (60.0 / seconds_passed))

    #Calc average
    if len( ReportWriteSpeed.processing_speed_per_minute) > 19: ReportWriteSpeed.processing_speed_per_minute.pop( 0)
    ReportWriteSpeed.processing_speed_per_minute.append( entries_per_minute)
    avg_speed = int( sum( ReportWriteSpeed.processing_speed_per_minute) / len( ReportWriteSpeed.processing_speed_per_minute))
    if avg_speed == 0: avg_speed = 1

    n_minutes_remaining = (aTotal - aCurrent) / avg_speed

    completion = datetime.now() + timedelta(minutes = n_minutes_remaining)
    if n_minutes_remaining <= 1:
        #            Blast Progress:
        speed_msg = "   Speed: %i/min. Estimated time of completion: %s (almost done)" % (avg_speed, completion.strftime('%H:%M'))
    elif n_minutes_remaining <= 60:
        speed_msg = "   Speed: %i/min. Estimated time of completion: %s (%i minutes remaining)" % (avg_speed, completion.strftime('%H:%M'), n_minutes_remaining)
    elif n_minutes_remaining < 60*24: #in 24h
        h = n_minutes_remaining / 60
        m = n_minutes_remaining - (h*60)
        speed_msg = "   Speed: %i/min. Estimated time of completion: %s (in %ih %imin)" % (avg_speed, completion.strftime('%H:%M'), h, m)
    else:
        h = n_minutes_remaining / 60
        speed_msg = "   Speed: %i/min. Estimated time of completion: %s (in >%i hours)" % (avg_speed, completion.strftime('%a %b %d %Y %H:%M'), h)

    ReportWriteSpeed.prev_time = now
    ReportWriteSpeed.prev_count = aCurrent

    if aPrintLock: aPrintLock.acquire()

    print "\n----------------------------------------------------------------------------"
    print   "Progress: Processed records %i/%i (%.1f%%)\n%s" % (aCurrent, aTotal, (float( aCurrent)/aTotal)*100.0, speed_msg)
    print   "----------------------------------------------------------------------------\n"

    if aPrintLock: aPrintLock.release()

    return speed_msg



def Remark( aRemark):
    if type( aRemark) == list: return ";".join( [re.sub('[,;]', '_', str( s)) for s in aRemark]) + ";"
    return re.sub('[,;]', '_', str( aRemark))

def FormatFastaDescription( aFastaDescription):

    aFastaDesction = re.sub( r"\|POI:\d+-\d+", "", aFastaDescription)
    aFastaDesction = re.sub( r"\|POISEQ:[A-Z-]+/[A-Z-]+", "", aFastaDescription)
    return Remark( aFastaDesction)


ASA_BB_THRESHOLDS = { 'G':14.0, }
ASA_SS_THRESHOLDS = { 'G':-1.0, 'A':14.0, }

def ASADecision( aSequence, aAsaValues, aSideChainOnly):

    global ASA_BB_THRESHOLDS, ASA_SS_THRESHOLDS

    if len( aAsaValues) == 0: return "Unknown"

    if len( aSequence) != len( aAsaValues):
        sys.stderr.write("ERROR: Sequence '%s' length %i does not match %i ASA values ('%s').\n" % (aSequence, len( aSequence),len( aAsaValues), aAsaValues))
        assert( len( aSequence) == len( aAsaValues))

    asa_thresholds_array = ASA_SS_THRESHOLDS if aSideChainOnly else ASA_BB_THRESHOLDS
    asa_threshold_values = []

    asa_vals = []

    for i in range( len( aSequence)):
        if aAsaValues[ i] < 0.0: continue
        asa_vals.append( aAsaValues[ i])
        aa = aSequence[ i]
        threshold = 15.0
        if aa in asa_thresholds_array: threshold = asa_thresholds_array[ aa]
        if threshold < 0.0: continue
        asa_threshold_values.append( threshold)

    if len( asa_threshold_values) == 0: return ("Buried", -1.0) #All Gs in sidechain mode

    threshold_avg = sum( asa_threshold_values) / float( len( asa_threshold_values))
    asa_avg = sum( asa_vals) / float( len( asa_vals))

    delta = asa_avg - threshold_avg

    decision = "Surface"
    if delta < 0.0: decision = "Buried"
    elif delta < 10.0: decision = "Intermediate"

    return (decision, asa_avg)


