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
import operator
import re
import math

pdbatoms_version = "1.0b"

#These residue names in the PDB file will be recognized as AMINO ACIDS
pdb_aas = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU", "GLY",
           "HIS","HSD","HSP","HSE","ILE","LEU","LYS", "MET",
           "PHE","PRO","SER","THR","TRP","TYR","VAL"] #, "GLX","ASX"]

res_map = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
           "GLY":"G","HIS":"H","HSD":"H","HSP":"H","HSE":"H","ILE":"I","LEU":"L","LYS":"K",
           "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V" }

#These residue names in the PDB file will be recognized as DNA/RNA
pdb_dna = {"DA":"a","DC":"c","DT":"t","DG":"g","A":"a","C":"c","T":"t","G":"g","U":"u"}

#These residue names in the PDB file will be recognized as IONS (single atom ions)
pdb_ion = ["MG","ZN","FE","NA","CL","K","CA","MN","AL","AU","AG","CU"] #list is incomplete

#These residue names in the PDB file will be recognized as OTHER
pdb_other = [ "SO4","PO4","ATP","ADP","NAD"] #list is incomplete

#These residue names in the PDB file will be recognized as WATER
pdb_water = ["HOH","HHO","OHH","SOL"]

#Van Der Waals radii for different atoms
VDW_radii = { "H" : 1.20, "C" : 1.70, "N" : 1.55, "O" : 1.52, "F" : 1.47,
              "NA": 2.27, "MG": 1.73, "SI": 2.10, "P" : 1.80, "S" : 1.80,
              "CL": 1.75, "K" : 2.75, "CU": 1.40, "ZN": 1.39, "SE": 1.90,
              "BR": 1.85, "I" : 1.98, "MN": 2.05, "AL": 1.84, "AU": 1.66,
              "FE": 1.94, "CA": 2.31, "AG": 2.03 }

#Post-translationally modified amino acids
#http://www.ebi.ac.uk/pdbe-srv/pdbechem/
PTM_aas = {
              #Phosphorylated
              "PTR":"TYR", # PHOSPHOTYROSINE
              "NEP":"HIS", # N1-PHOSPHONOHISTIDINE
              "HIP":"HIS", # Nd1-PHOSPHONOHISTIDINE
              "SEP":"SER", # PHOSPHOSERINE; PHOSPHORYLATED SER
              "TPO":"THR", # PHOSPHOTHREONINE;PHOSPHORYLATED;THREONINE;PHOSPHORYLATION
              "PHD":"ASP", # ASPARTYL PHOSPHATE
              "PTM":"TYR", # ALPHA-METHYL-O-PHOSPHOTYROSINE
              "FTY":"TYR", # DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
              "SBG":"SER", # O-[(S)-HYDROXY(METHYL)PHOSPHORYL]-L-SERINE
              "SEN":"SER", # O-[N,N-dimethylphosphoramidate]-L-serine
              "SGB":"SER", # O-[(S)-METHYL(1-METHYLETHOXY)PHOSPHORYL]-L-SERINE
              "SGX":"SER", # O-[(S)-AMINO(METHOXY)PHOSPHORYL]-L-SERINE
              "SGR":"SER", # O-[(R)-AMINO(METHOXY)PHOSPHORYL]-L-SERINE
              "4TP":"THR", # 4-HYDROXY-L-THREONINE-5-MONOPHOSPHATE
              "28X":"THR", # O-[(R)-{[(3R)-3,4-dihydroxybutyl]oxy}(hydroxy)phosphoryl]-L-threonine
              #Acetylated
              "AYA":"ALA", #N-ACETYLALANINE; N-ACETYL ALANINE
              "TH5":"THR", #O-ACETYL-L-THREONINE
              "ALY":"LYS", #N(6)-ACETYLLYSINE
              "SAC":"SER", #N-ACETYL-SERINE
              "OAS":"SER", #O-ACETYLSERINE
              "SCY":"CYS", #S-ACETYL-CYSTEINE
              #Mono-methylation
              "MLZ":"LYS", #N-METHYL-LYSINE
              "MVA":"VAL", #N-METHYLVALINE
              "MLE":"LEU", #N-METHYLLEUCINE
              "MAA":"ALA", #N-METHYLALANINE
              "MEA":"PHE", #N-METHYLPHENYLALANINE
              "0A9":"PHE", #METHYL L-PHENYLALANINATE
              "MME":"MET", #N-METHYL METHIONINE
              "4PH":"PHE", #4-METHYL-L-PHENYLALANINE
              "SMC":"CYS", #S-METHYLCYSTEINE
              "CMT":"CYS", #O-METHYLCYSTEINE
              "HIC":"HIS", #4-METHYL-HISTIDINE
              "AGM":"ARG", #5-METHYL-ARGININE
              "IML":"ILE", #N-METHYL-ISOLEUCINE
              "MLL":"LEU", #METHYL L-LEUCINATE
              "MVA":"VAL", #METHYLVALINE
              "VME":"VAL", #O- METHYLVALINE
              "MEN":"ASN", #N-METHYL ASPARAGINE
              "IML":"ILE", #N-METHYL-ISOLEUCINE
              "OLT":"THR", #O-METHYL-L-THREONINE
              #SelenoMethionine (Not a PTM)
              "MSE":"MET"
              }


PTM_types = {"Phosphorylation" : ["PTR","NEP","HIP","SEP","TPO","PHD","PTM","FTY","SBG","SEN","SGB","SGX","SGR","4TP","28X"],
             "Acetylation":      ["AYA","TH5","ALY","SAC","OAS","SCY"],
             "Methylation":      ["MLZ","MVA","MLE","MAA","MEA","0A9","MME","4PH","SMC","CMT","HIC","AGM","IML","MLL","MVA","VME","MEN","IML","OLT"],
             "Selenomethionine": ["MSE"] }

atoms_of_residue = { "ALA": ["N","CA","C","O","CB"],
                     "ARG": ["N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"],
                     "ASN": ["N","CA","C","O","CB","CG","OD1","ND2"],
                     "ASP": ["N","CA","C","O","CB","CG","OD1","OD2"],
                     "ASX": ["N","CA","C","O","CB","CG","XD1","XD2"],
                     "CYS": ["N","CA","C","O","CB","SG"],
                     "GLN": ["N","CA","C","O","CB","CG","CD","OE1","NE2"],
                     "GLU": ["N","CA","C","O","CB","CG","CD","OE1","OE2"],
                     "GLX": ["N","CA","C","O","CB","CG","CD","XE1","XE2"],
                     "GLY": ["N","CA","C","O"],
                     "HIS": ["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"],
                     "HSD": ["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"],
                     "HSP": ["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"],
                     "HSE": ["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"],
                     "ILE": ["N","CA","C","O","CB","CG1","CG2","CD1"],
                     "LEU": ["N","CA","C","O","CB","CG","CD1","CD2"],
                     "LYS": ["N","CA","C","O","CB","CG","CD","CE","NZ"],
                     "MET": ["N","CA","C","O","CB","CG","SD","CE"],
                     "PHE": ["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"],
                     "PRO": ["N","CA","C","O","CB","CG","CD"],
                     "SER": ["N","CA","C","O","CB","OG"],
                     "THR": ["N","CA","C","O","CB","OG1","CG2"],
                     "TRP": ["N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"],
                     "TYR": ["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"],
                     "VAL": ["N","CA","C","O","CB","CG1","CG2"] }

#pdb_ion array needs to con tain all existing ion keys in VDW_radii

#CA (Calsium) would be confused with C-alphas which are also marked as CA
#HG (gamma-hydrogens) would be confused with Mercury (Hg)
#HE (epsilon-hydrogens) would be confused with Helium (He)
#PB (Lead) could be problematic as well

def PTM_Type( aResname):
    global PTM_types
    resname = aResname.upper()
    for ptm_type, names in PTM_types.iteritems():
        if resname in names: return ptm_type
    return ""

def ConvertibleToAA( aResname ):

    global PTM_aas, res_map
    resname = aResname.upper()
    if resname in PTM_aas: return PTM_aas[ resname]
    if resname in res_map: return resname
    return "?"

#Convert Post-translationally modified aas from HETATM entries to ATOM entries with regular aa residue names
def ConvertPTMAAstoRegularAtoms( aPdbInfo, aRemovePTMs=False, aQuiet=False):

    n_residues_converted = 0
    residue_list = aPdbInfo["residues"]
    atom_list = aPdbInfo["atoms"]
    stats = aPdbInfo["stats"]

    for residue in residue_list:
        orig_resname = residue.ResName()
        n_atoms_converted, n_atoms_removed = residue.ConvertToAAIfPossible( aRemovePTMs)
        #if n_atoms_converted and not aQuiet: print "Converted PTM: %s %s conv:%i, remove:%i" % (residue.ResNum(), residue.ResName(), n_atoms_converted, n_atoms_removed)
        if n_atoms_converted and not aQuiet: print "PDB INFO: Converted %s -> %s %s %s" % (orig_resname, residue.ResName(), residue.Chain(), residue.ResNum())
        if n_atoms_converted > 0:
            n_residues_converted += 1
            stats["modified_residues"] += 1
            stats["hetatms"] -= n_atoms_converted
            stats["hetatms"] -= n_atoms_removed
            stats["total_atom_records"] -= n_atoms_removed
            stats["atoms"] += n_atoms_converted
            stats["amino_acids"] += 1
            stats["other_residues"] -= 1
            UpdateStatsAtomRatios( stats, residue)

    if aRemovePTMs:
      for atom in reversed( atom_list):
          if atom.isptmaminoacid:
              removable = not atom.IsBBorSideChain()
              #print atom, ("- %i" % removable)
              if removable: atom_list.remove( atom)


    return n_residues_converted

def IsAA( lettercode3 ):

    global pdb_aas

    if lettercode3.upper() in pdb_aas:
        return True

    return False

def IsPTMAA( lettercode3 ):

    global PTM_aas
    return lettercode3.upper() in PTM_aas


def MaxVDWradius():
    return max( VDW_radii.iteritems(), key=operator.itemgetter(1))[1]

def chunkify( aStr, aLen):
    if aLen < 1: return aStr
    return ("\n".join( [aStr[i:i+aLen] for i in range(0, len(aStr), aLen)]) + "\n")

def rangify( aNumbers, aSeparator=",", aSortEm=True):

    if len( aNumbers) < 1: return ""
    elif len( aNumbers) == 1: return str( numbers[ 0])

    if aSortEm: aNumbers = sorted( aNumbers)

    ranges = sum((list(t) for t in zip(aNumbers, aNumbers[1:]) if t[0]+1 != t[1]), [])
    iranges = iter(aNumbers[0:1] + ranges + aNumbers[-1:])

    retvals = []

    for n in iranges:
        start = n
        end = next(iranges)
        if start != end: retvals.append(  "%i-%i" % (start, end))
        else: retvals.append( "%i" % (start))

    return aSeparator.join( retvals)

def derangify( aRangeStr):

    retvals = []
    ranges = aRangeStr.split(",")

    for cur_range in ranges:
        parts = cur_range.split( "-")
        try:
          parts = map(int, parts)
        except ValueError as ex:
          sys.stderr.write( "ERROR: %s.\n" % str( ex))
          continue
        start = min( parts)
        end = max( parts)

        for i in range( start, end+1):
            retvals.append( i)

    return sorted( retvals)

#aResNum should be str and the include insertion code if it exists
#aModel = -1 (first matching residue in file)
#aAdjacentPositions = other residues next to the queried residue to be returned, negative values = towards the N-terminus
#             if a non-zero number, residues are returned in an array from N- to C-terminus
#Returns an array of residues
def FindResidues( aResidueList, aResNum, aChain, aModel=-1, aAdjacentPositions=None, aIncludeAlternateLocations=False):

    if type( aResNum) == int:
        sys.stderr.write( "RESIDUE WARNING: Argument 'aResNum' should be of type str and may include insertion codes.\n")
        aResNum = str( aResNum)

    aResNum = aResNum.strip()

    for r in range( len( aResidueList)):
        res = aResidueList[ r]

        if not aIncludeAlternateLocations and res.IsAlternateLocation(): continue #Not necessary, but does no harm to check here

        #print "Resnum: '%s' c:'%s'" % (res.ResNum(), res.Chain())
        if aResNum == res.ResNum() and aChain == res.Chain():

            if aModel >= 0 and aModel != res.Model(): continue

            retarr = [ res]
            if aAdjacentPositions == None or aAdjacentPositions == 0: return retarr #RETURN, only queried res

            step = 1 if aAdjacentPositions > 0 else -1
            cur_pos = r
            filled_positions = 0

            while True:

                cur_pos += step
                pos_res = aResidueList[ cur_pos]

                if not pos_res.IsAminoAcid(): continue

                if cur_pos >= len( aResidueList) or cur_pos < 0:
                    #End of residue list reached
                    break

                if aChain == pos_res.Chain() and (aModel < 0 or aModel == pos_res.Model()):
                    if not aIncludeAlternateLocations and pos_res.IsAlternateLocation(): continue
                    if step > 0: retarr.append( pos_res)
                    else: retarr.insert( 0, pos_res)
                    if pos_res.IsPrimaryLocation(): filled_positions += 1
                else:
                    #Chain or model ended before all position were filled
                    break

                if filled_positions >= aAdjacentPositions:
                    #All residues found
                    break

            return retarr

            #for i in range( r+step, (r+step)+aPositions, step):
            #    if i < len( aResidueList) and i >= 0:
            #        pos_res = aResidueList[ i]
            #        if aChain == pos_res.Chain() and (aModel < 0 or aModel == pos_res.Model()):
            #            if step > 0: retarr.append( aResidueList[ i])
            #            else: retarr.insert( 0, aResidueList[ i])
            #        else: break
            #    else: break

            #return retarr

    sys.stderr.write( "RESIDUE ERROR: Could not find residue '%s' in chain '%s'. (model:'%s')\n" % (aResNum, aChain, aModel))
    return []


#aDistance between Alpha-Carbons in Ångströms
#Returns a list of residues that are in aDistance range
def ResiduesInRange( aResidue, aResidueList, aDistance=10.0):

    if not aResidue or aResidue == None: return []

    in_range_residues = []

    ca_coord = aResidue.CoordsCA()
    if ca_coord[ 0] == 0.0 and ca_coord[ 1] == 0.0 and ca_coord[ 2] == 0.0:
        sys.stderr.write( "RESIDUE ERROR: Could not find CA coordinates for %s %s. (model:%s)\n" % (aResidue.ResNum(), aResidue.Chain(), aResidue.Model()))
        return []

    for res_other in aResidueList:

        distance = Residue.CoordDistance( ca_coord, res_other.CoordsCA())
        #print "Distance:", ca_coord, " To ", res_other.CoordsCA(), "  == ", distance
        if distance <= aDistance:

            if aResidue != res_other: in_range_residues.append( res_other)

    return in_range_residues


#Checks which chains among aChains could be duplicates based on their sequences
#aChains = chains that will be compared together (empty for all)
#aChainsToReport = returned dict will contain the duplicate chains in the set of aChains for these chains (empty for all aChains)
#aPercentageThreshold = Required aa identity percentage for chains to be declared duplicates
#Returns dict with chain as key (aChainsToReport), duplicate chain list as value
def SimpleFindDuplicateChains( aPdbInfo, aChains = [], aChainsToReport = [], aPercentageThreshold=80.0, aQuiet=False, aPrintPercentages=False ):

    if len( aChains) == 0:
      aChains = aPdbInfo["stats"]["chains"] #All chains

    if len( aChainsToReport) == 0:
      aChainsToReport = aChains

    duplicates = {}
    not_duplicates = {}

    models = FastaSequence( aPdbInfo["residues"], aPerModelAndChain=True, aInsertGaps=True, aInitialGaps=True, aAtomTypes=["A","D"], aQuiet=False)
    if len( models) == 0:
        sys.stderr.write("ERROR: No models found for duplicate processing.\n")
        return {}

    #Only look at first model
    #seq_tuple = models[ sorted( models.keys())[ 0]]
    chain_list, sequences = models[ sorted( models.keys())[ 0]]

    warned = []

    for c in aChainsToReport:

        #print "Chain: ", c
        new_duplicates = set()
        skippable = set()

        for d in duplicates.keys():
            if c in duplicates[ d]: new_duplicates.add( d)

        for d in not_duplicates.keys():
            if c in not_duplicates[ d]: skippable.add( d)

        for other in aChains:
            if other == c: continue #Do not compare to self
            if other in new_duplicates or other in skippable:
                #print "Skipped %s: dup" % other, other in new_duplicates, " tested:", other in skippable
                continue #This pair has already been processed

            seq_c = ""
            seq_other = ""

            try:
                index_c = chain_list.index( c)
                seq_c = sequences[ index_c]
            except ValueError:
                if not aQuiet and c not in warned:
                    sys.stderr.write("WARNING: Chain '%s' not found in duplicate processing.\n" % c)
                    warned.append( c)
                continue

            try:
                index_other = chain_list.index( other)
                seq_other = sequences[ index_other]
            except ValueError:
                if not aQuiet and other not in warned:
                    sys.stderr.write("WARNING: Chain '%s' not found in duplicate processing.\n" % other)
                    warned.append( other)
                continue

            min_range = min( len( seq_c), len( seq_other))
            matched = 0
            unmatched = 0

            for i in range( min_range):
                aa_c = seq_c[ i]
                aa_other = seq_other[ i]
                if aa_c != '-' or aa_other != '-':
                    if aa_c == aa_other: matched += 1
                    else: unmatched += 1

            total = matched + unmatched
            if total == 0: continue

            percentage = (matched / float(total)) * 100.0
            if aPrintPercentages: sys.stdout.write( "%s vs. %s: %.2f%%\n" % (c, other, percentage))

            if percentage >= aPercentageThreshold-0.00001: new_duplicates.add( other)
            else: skippable.add( other)

        not_duplicates[ c] = list( skippable)
        duplicates[ c] = list( new_duplicates)

    return duplicates

#Returns either string sequence or dict of models[ modelnum] = ([chains],[sequences]) if aPerModelAndChain=True
def FastaSequence( aResidueList, aAtomTypes=["A","D","I"], aInsertGaps=False, aInitialGaps=False, aSplitToRowsOfLength=0, aPerModelAndChain=False, aQuiet=False ):

    global res_map, pdb_dna

    seq = ""
    cur_chain = None
    cur_model = None
    expected_residue_seqnum = 1
    prev_icode = ""
    prev_reseqnum = 0
    first = True
    prev_residue_is_hetam = False
    warnings = 0
    models = {}
    #chains = []
    #chain_seqs = []
    first_in_chain = True

    for residue in aResidueList:

        if first:

            if not aInitialGaps: expected_residue_seqnum = residue.ResSeqNum()
            #elif residue.ResSeqNum() > 0: seq += (residue.ResSeqNum()-1) * "-"
            prev_icode = residue.InsertionCode()

            models[ residue.Model()] = ([], [])
            #models[ residue.Model()][ 0].append( residue.Chain())

        if not residue.IsPrimaryLocation(): continue #Exclude alternative locations of residues from sequence

        #When model or chain has changed
        if (cur_chain and cur_chain != residue.Chain()) or (cur_model and cur_model != residue.Model()):

            prev_residue_is_hetam = False
            first_in_chain = True
            prev_icode = residue.InsertionCode()
            if aPerModelAndChain:
                #if len( seq):
                models[ cur_model][ 1].append( seq)
                models[ cur_model][ 0].append( cur_chain)
                seq = ""
            else:
                seq += "/" #chain or model separator

            if cur_model and cur_model != residue.Model() and residue.Model() not in models:
                #New model
                models[ residue.Model()] = ([], [])


            if aInitialGaps: expected_residue_seqnum = 1
            else: expected_residue_seqnum = residue.ResSeqNum() #Set to whatever the new number is
        elif not first:
            first_in_chain = False

        first = False

        #No gaps between HETATM or ATOM and HETATM entries
        het_atom = residue.IsHETATOM()
        if het_atom or prev_residue_is_hetam: expected_residue_seqnum = residue.ResSeqNum()
        prev_residue_is_hetam = het_atom

        #No gaps at icode changes
        icode = residue.InsertionCode()
        if icode != prev_icode: expected_residue_seqnum = residue.ResSeqNum()
        prev_icode = icode

        #Insert gaps if residue seq num is not consecutive
        if aInsertGaps and expected_residue_seqnum > 0:
            if expected_residue_seqnum+10000 < residue.ResSeqNum():
                seq += "-"*10000 #Curb to 10000 lines
            else:
                while expected_residue_seqnum < residue.ResSeqNum():
                    seq += "-"
                    expected_residue_seqnum += 1

        #Append to sequence if atom type is suitable
        if residue.AtomType() in aAtomTypes: seq += residue.OneLetterResName()

        #Sanity check
        if residue.ResSeqNum() < expected_residue_seqnum and not first_in_chain:
            warnings += 1
            if not aQuiet:
                if warnings <= 10: sys.stderr.write("WARNING: Residue numbering is not consistent for residue: %s %s %s\n" % (residue.ResNum(), residue.ResName(), residue.Chain() ))
                if warnings == 10: sys.stderr.write("         Only first 10 warned.\n")

        expected_residue_seqnum = residue.ResSeqNum() + 1

        cur_chain = residue.Chain()
        cur_model = residue.Model()

    if aPerModelAndChain:
        #if len( seq):
        models[ cur_model][ 0].append( cur_chain)
        models[ cur_model][ 1].append( seq)
        return models

    #All done
    return chunkify( seq, aSplitToRowsOfLength)


def PrintFastaPerModelAndChain( aFastaSeqPerModelAndChain, aOnlyFirstModel=True ):

  #print aFastaSeqPerModelAndChain
  #return

  model_nums = sorted( map( int, aFastaSeqPerModelAndChain.keys()))

  for model_num in model_nums:

      model = aFastaSeqPerModelAndChain[ str( model_num)]
      chains = model[ 0]
      seqs = model[ 1]
      if model_num >= 0: print "Model %s:" % str( model_num)

      for i in range( 0, len( max(chains, seqs))):
          print "Chain %s:" % chains[ i]
          print "Seq:\n%s" % seqs[ i]

      if aOnlyFirstModel: break


def GetModificationType( aPDBStats, aResidue ):

    short_desc = aResidue.ShortStrDescription()
    cols = short_desc.split(";")
    index = False
    for col in cols:
        if col.startswith( "I:"):
            index = col
            break

    assert( index != False)

    #Compare indices
    for key in aPDBStats["ptm_aas"]:
       if key.find( index) >= 0: return aPDBStats["ptm_aas"][ key]

    #if short_desc in aPDBStats["ptm_aas"]: return aPDBStats["ptm_aas"][ short_desc]
    return "None"

def GetShortModificationType( aPDBStats, aResidue ):
    ptm_type = GetModificationType( aPDBStats, aResidue)

    if ptm_type == "Phosphorylation": return "P"
    elif ptm_type == "Acetylation": return "A"
    elif ptm_type == "Methylation": return "M"
    elif ptm_type == "Selenomethionine": return "S"

    return "N" #None


#Writes a list of atoms to file
def WriteToFile( openfilehandle, atom_list, types = [], chains=[] ):
#def WriteAtomsToFile( openfilehandle, atom_list, types = [], chains=[] ):

    handle_opened = False
    if type( openfilehandle) is str:
      openfilehandle = open( openfilehandle, "w")
      handle_opened = True

    all_types = ["A","D","I","U","W","O"]

    if len( types) == 0:
        types = all_types

    for t in types:
        if t not in all_types: sys.stderr.write("WARNING: Unknown atom type: '%s' ignored.\n" % t)


    for atom in atom_list:

        if ("A" in types and atom.isaminoacid) or \
           ("W" in types and atom.iswater) or \
           ("D" in types and atom.isdna) or \
           ("U" in types and atom.isunknown) or \
           ("I" in types and atom.ision) or \
           ("O" in types and atom.isother):

            if len( chains) and atom.chain not in chains: continue
            openfilehandle.write( atom.FileDescription())


    openfilehandle.write( "END\n")

    if handle_opened: openfilehandle.close()


#Returns dict containing keys "atoms", "residues", "stats"
#under key "atoms" there is a list of the atoms of the pdb in the order they were read
#under key "residues" there is a list of the residues of the pdb in the order they were read
#under key "stats" there is a dict of different statistics of the PDB file
#aConvertPTMs==False : No modifications to Post translational modifcations
#aConvertPTMs=="C"   : Convert PTMs to regular Amino acids (remove PTM atoms)
#aConvertPTMs==True  : Convert PTMs to regular Amino acids (remove PTM atoms)
#aConvertPTMs=="A"   : Convert PTMs from HETATM to ATOM and keep all PTM atoms
def ReadPDB( infile, chains=[], models=[], firstmodelonly=False, aQuiet=False, aConvertPTMs=False, aWarnUnknown=True ):

    #Convert chain format if necessary
    #Chain "EE" means really chain e   (in BLAST results)
    for i in range( len( chains)):
        chain = chains[ i]
        if len( chain) == 2 and chain.upper() and chain.isalpha() and chain[ 0] == chain[ 1]:
            chains[ i] = chain[ 0].lower()
            #print " -> CHAIN__: %s" % chains[ i]

    #Make sure models are int
    for m in reversed( range( len( models))):
        if type( models[ m]) != int:
          try:
              models[ m] = int( models[ m])
          except Exception:
              sys.stderr.write( "ERROR: Could not convert model '%s' to integer." % models[ m])
              models.remove[ m]

    if type( infile) == file:
        return DoReadPDB( infile, chains, models, firstmodelonly, aQuiet=aQuiet, aConvertPTMs=aConvertPTMs, aWarnUnknown=aWarnUnknown)
    elif type( infile) == str:
        try:
            with open( infile, "r") as filehandle:
               return DoReadPDB( filehandle, chains, models, firstmodelonly, aQuiet=aQuiet, aConvertPTMs=aConvertPTMs, aWarnUnknown=aWarnUnknown)
        except IOError as ex:
            sys.stderr.write("PDB ERROR: %s\n" % (str( ex)))
    elif type( infile) == list:
        retvals = []
        for f in infile:
            retvals.append( ReadPDB( f, chains, models, firstmodelonly, aQuiet=aQuiet, aConvertPTMs=aConvertPTMs, aWarnUnknown=aWarnUnknown))
        return retvals
    else:
        raise TypeError("Invalid type for argument infile.")

#infile is an open filehandle
def DoReadPDB( infile, chains = [], models = [], firstmodelonly=False, aQuiet=False, aConvertPTMs=False, aWarnUnknown=True):

    global PTM_types

    atoms = []
    residues = []
    remarks = []

    stats = { "amino_acids":0,
              "dna":0,
              "ions":0,
              "ion_names":[],
              "waters":0,
              "n_chains":0,
              "n_models":0,
              "models":[],
              "ignored_models":[],
              "hydrogens":0,
              "chains":[],
              "n_ignored_chains":0,
              "ignored_chains":[],

              "n_anisou" : 0, #Anisotropic Temperature Factor -lines
              "non-atom_lines":0,
              "total_atom_records":0,  #atoms + hetatm
              "ignored_atom_records":0,
              "accepted_atom_records":0,
              "accepted_aa_atom_records":0,
              "hetatms":0,
              "atoms":0,
              "ignored_atoms":0,
              "ignored_hetatms":0,

              "n_residues":0,
              "other_residues":0,
              "unknown_residues":0,
              #"ignored_residues":0,
              "unknown_radius":0,
              "largest_radius":1.2,
              "atom_types": {},
              "atoms_per_chain":{},
              "resolution": 0.0,
              "method":"unknown",
              "R-Free":0.0,
              "R-Value":0.0,
              "TempF-avg":0.0,
              "TempF-min":999.0,
              "TempF-max":0.0,
              "title":"",
              "compound":"",
              "header":"",
              "published":"",
              "pdbcode":"",
              "keywords":"",
              "organism_scientific": "",
              "organism_common": "",
              "organism_taxid": 0,
              "erroneous_atom_lines": 0,
              "n_atoms_with_alt_locations": 0,
              "n_residues_with_altloc" : 0,
              "n_aas_with_altloc" : 0,
              "n_alternate_location_residues" : 0, #Number of non-primary locations
              "alternate_locations" : {}, #Contains residues woth alternate locations, see globvar res_id_templ for format
              "n_ptm_aas":0, #Post-translationally modified amino acids (including each alternate location for residues)
              "ptm_aas": {}, #Contains Post-translationally modified residues, see globvar res_id_templ for format
              "modified_residues":0,
              "n_biomolecules":0,
              "incomplete_residue_ratio":0.0, # Percentage of incomplete amino acid entries
              "single_atom_residue_ratio": 0.0,
              "backbone_only_ratio": 0.0
            }

    try:
      stats["filename"] = infile.name
    except:
      stats["filename"] = "unknown"

    prev_chain = "?"
    cur_model = -1
    prev_model = -1
    first_in_model = -1
    prev_atom_resname = "???"
    prev_atom_insert = "?"
    prev_resnum = -999999

    #Some pdb files separate residues by an "insertion code"
    #This means that the residue number does not change only the code: GLU 100 -> TYR 100A
    prev_insert_code = "?"
    prev_altloc_indicator = "?"
    insertions_reported = False
    alt_locations_reported = False

    residue_reported = -999999
    warned_atom_types = {}

    #running_atom_id = 0
    accepted_atoms = 0

    line_num = 0

    cur_resid = Residue()
    newAtom = None

    chains = filter( None, chains) #Remove empty strings
    chains_set = len( chains) > 0

    ##with open( filename, 'r') as f:

    readinginfo = True
    atominfos = ["MODEL ","ATOM  ","HETATM"]

    #Read line by line
    for line in infile:

        line_num += 1

        #Skip empty lines
        if len( line) < 2: continue

        litem = line[0:6] #line-item

        #Parse non-atom information until first atom (or model, hetmatm) line if found
        if readinginfo and litem not in atominfos:
            ParseRemark( line, stats)
            if line.startswith("REMARK"): remarks.append( line)
            continue

        readinginfo = False

        try:

            #Models
            if litem == "MODEL ":
                stats["n_models"] += 1
                try:
                    model_str = line[6:].strip()
                    prev_model = cur_model
                    cur_model = int( model_str)
                    first_in_model = cur_model

                    if cur_resid.ModelNum() == -1: cur_resid.SetModelnum( cur_model) #When first model found in file

                    #if first_model == -1: first_model = cur_model
                    if prev_model == -1 and firstmodelonly:
                        models = [cur_model] #Set argument to contain only the first model
                        stats["models"].append( cur_model)
                    elif cur_model >= 0 and cur_model not in stats["models"] and len( models) and cur_model in models:
                        stats["models"].append( cur_model)

                except ValueError:
                    if not aQuiet: sys.stderr.write( "WARNING: Bad model number found: '%s' on line %i.\n" % (model_str, line_num))
            #elif firstmodelonly and first_model >= 0 and cur_model != first_model:
            #    continue #Skip all lines but the first model in file if firstmodelonly argument enabled
            elif litem == "ATOM  " or litem == "HETATM":

                #print "Model %i in models: %i" % (cur_model, cur_model in models), models

                #Models
                if cur_model >= 0 and len( models) and cur_model not in models:
                    if cur_model not in stats["ignored_models"]: stats["ignored_models"].append( cur_model)
                    #Skip model
                    continue

                #running_atom_id += 1
                #newAtom = Atom( line, running_atom_id, cur_model)
                newAtom = Atom( line, -1, cur_model)

                stats["total_atom_records"] += 1

                #if newAtom.locind != " ":
                  #stats["n_atoms_with_alt_locations"] += 1
                  #if stats["n_atoms_with_alt_locations"] == 1:
                    #Warn once
                    #sys.stderr.write( "WARNING: file '%s' contains alternate locations for atoms.\n" % infile.name)

                if first_in_model >= 0:
                    newAtom.first_in_model = first_in_model
                    first_in_model = -1

                #Skip unwanted chains
                if chains_set and newAtom.chain not in chains:
                    if newAtom.chain not in stats["ignored_chains"]:
                        stats["ignored_chains"].append( newAtom.chain)
                        stats["n_ignored_chains"] += 1
                    continue #skip chain

                #skip hydrogens
                if newAtom.ishydrogen:
                    stats["hydrogens"] += 1
                    stats["ignored_atom_records"] += 1
                    if newAtom.isatom: stats["ignored_atoms"] += 1
                    elif newAtom.ishetatm: stats["ignored_hetatms"] += 1
                    continue

                #skip water
                if newAtom.iswater:
                    stats["waters"] += 1
                    stats["ignored_atom_records"] += 1
                    if newAtom.isatom: stats["ignored_atoms"] += 1
                    elif newAtom.ishetatm: stats["ignored_hetatms"] += 1
                    continue

                if newAtom.isother or newAtom.isptmaminoacid:
                    stats[ "other_residues"] += 1

                #unknown types
                if newAtom.isunknown:
                    stats[ "unknown_residues"] += 1

                    if residue_reported != newAtom.resnum:
                        if not aQuiet and aWarnUnknown: sys.stderr.write( "PDB INFO: Unknown residue name on line %04d: %s\n" % (line_num, newAtom.resname))
                        residue_reported = newAtom.resnum

                if newAtom.hasradius:

                    stats["accepted_atom_records"] += 1

                    if newAtom.isaminoacid:
                        stats["accepted_aa_atom_records"] += 1

                    if newAtom.isatom: stats["atoms"] += 1
                    elif newAtom.ishetatm: stats["hetatms"] += 1

                    if newAtom.VDW_type not in stats["atom_types"]: stats["atom_types"][ newAtom.VDW_type] = 0
                    stats["atom_types"][ newAtom.VDW_type] += 1

                    #set chain stats for new atom
                    if newAtom.chain not in stats["chains"]:
                        stats["chains"].append( newAtom.chain)
                        stats["atoms_per_chain"][newAtom.chain] = 0
                        stats["n_chains"] += 1

                    if newAtom.radius > stats["largest_radius"]: stats["largest_radius"] = newAtom.radius
                    stats["atoms_per_chain"][newAtom.chain] += 1

                else:
                    #Ignore and warn about atoms that have no defined radius
                    stats["unknown_radius"] += 1
                    stats["ignored_atom_records"] += 1
                    if newAtom.isatom: stats["ignored_atoms"] += 1
                    elif newAtom.ishetatm: stats["ignored_hetatms"] += 1

                    if newAtom.atomname not in warned_atom_types.keys():
                        if not aQuiet:
                            sys.stderr.write( "WARNING: No VDW radius could be determined for atom on line %04d: %s %i %s\n" % (line_num, newAtom.atomname, newAtom.resnum, newAtom.resname))
                            sys.stderr.write( "         Atom is ignored. This warning is displayed only once per atom type.\n")
                        warned_atom_types[ newAtom.atomname] = 1
                    continue

                #DEBUG
                #if cur_resid.NumberOfAtoms() > 0:
                #    print newAtom.resname, newAtom.resnum, ("'%s'" % newAtom.AltLocationIndicator())
                #    if int( cur_resid.ResNum()) > 184: sys.exit( 0)

                #Residue number or chain has changed
                if newAtom.resnum != prev_resnum or newAtom.chain != prev_chain or newAtom.Insertion() != prev_insert_code:

                    if cur_resid.NumberOfAtoms() > 0:
                        AddResidue( cur_resid, residues, atoms, stats)

                        #New Residue
                        cur_resid = Residue()
                        cur_resid.SetModelnum( cur_model)

                elif newAtom.resname != prev_atom_resname and newAtom.resnum == prev_resnum:

                    #Residue name has changed but the numbering did not
                    if newAtom.Insertion() != prev_insert_code and not insertions_reported:
                        if not aQuiet: sys.stderr.write( "INFO: File '%s' contains inserted residues, 1st on line: %04d\n" % (infile.name, line_num))
                        insertions_reported = True
                    if newAtom.AltLocationIndicator() != prev_altloc_indicator and not alt_locations_reported:
                        if not aQuiet: sys.stderr.write( "INFO: File '%s' contains alternative locations for residues, 1st on line: %04d\n" % (infile.name, line_num))
                        alt_locations_reported = True

                    if newAtom.Insertion() == prev_insert_code and newAtom.AltLocationIndicator() == prev_altloc_indicator and not aQuiet:
                        sys.stderr.write( "WARNING: residue numbering does not change between residues (%s->%s) in file '%s' on line %04d: '%s'\n" % (prev_atom_resname, newAtom.resname, infile.name, line_num, line[0:30]))


                prev_chain = newAtom.chain
                prev_resnum = newAtom.resnum
                prev_insert_code = newAtom.Insertion()
                prev_altloc_indicator = newAtom.AltLocationIndicator()
                prev_atom_resname = newAtom.resname
                prev_atom_insert = newAtom.insert

                #add atom to residue
                cur_resid.AddAtom( newAtom)
                accepted_atoms += 1




                #Cumulative avg
                stats[ "TempF-avg"] = (stats[ "TempF-avg"]*(accepted_atoms-1) + float(newAtom.b)) / accepted_atoms
                stats[ "TempF-min"] = min( float(newAtom.b), stats[ "TempF-min"])
                stats[ "TempF-max"] = max( float(newAtom.b), stats[ "TempF-max"])

            elif litem == "ANISOU":
                stats["n_anisou"] += 1
                stats["non-atom_lines"] += 1
                #Add ANISOU information to previous atom here
            else:
                stats["non-atom_lines"] += 1
                #print "LINE:'%s'" % line.rstrip()
                if litem[0:3] == "TER" or litem == "ENDMDL":
                    if cur_resid and len( cur_resid.atoms) > 0: cur_resid.SetTerminal( True) #cur_resid.atoms[ -1].terminal = True
                    #if len( atoms): atoms[ -1].terminal = True

        except ValueError as err:
            if not aQuiet: sys.stderr.write( "PDB ERROR: Bad ATOM or HETATM entry on line %04d: %s\n%s" % (line_num, err.args, line))
            stats[ "erroneous_atom_lines"] += 1


    ######################
    #All lines processed
    stats[ "total_lines"] = line_num

    #Append last residue in file
    if cur_resid.NumberOfAtoms() > 0: AddResidue( cur_resid, residues, atoms, stats)

    ret_dict = { "atoms":atoms, "residues":residues, "stats":stats, "remarks":remarks }

    #PTM conversions
    if aConvertPTMs and aConvertPTMs != False and stats["n_ptm_aas"] > 0:
        if aConvertPTMs == "C" or aConvertPTMs == True: ConvertPTMAAstoRegularAtoms( ret_dict, aRemovePTMs=True, aQuiet=aQuiet)
        elif aConvertPTMs == "A": ConvertPTMAAstoRegularAtoms( ret_dict, aRemovePTMs=False, aQuiet=aQuiet)

    return ret_dict

#Model;Chain;ResName;ResNum;Alternate location;Convertible to AA (Original)
#res_id_templ = "M:%s;C:%s;N:%s;R:%s;L:%s;O:%s;"

def AddResidue( aResidue, aResiduesList, aAtomList, aStats):

    #global res_id_templ

    residue_list = []
    if aResidue.HasAltLocation():
        residue_list = SplitResidueToAlternateLocations( aResidue, aAtomList)
    else:
        residue_list.append( aResidue)


    for cur_resid in residue_list:

        #Append residue to returned array
        cur_resid.SetResIndex( len( aResiduesList))
        aResiduesList.append( cur_resid)

        #APPEND atoms to AtomList and set unique ids
        for atom in cur_resid.atoms:
            atom.atom_id = len( aAtomList)+1
            aAtomList.append( atom)

        primary = cur_resid.IsPrimaryLocation()
        res_name = cur_resid.ResName()

        #Update stats
        if primary:
            aStats["n_residues"] += 1

        if cur_resid.IsAminoAcid() and primary:
            aStats["amino_acids"] += 1
            UpdateStatsAtomRatios( aStats, cur_resid)

        if cur_resid.HasAltLocation():
            if primary:
                aStats["n_residues_with_altloc"] += 1
                aStats["n_atoms_with_alt_locations"] += cur_resid.NumberOfAtoms()
                if cur_resid.IsAminoAcid(): aStats["n_aas_with_altloc"] += 1
            else:
                aStats["n_alternate_location_residues"] += 1
            res_id_str = cur_resid.ShortStrDescription() #res_id_templ % (cur_resid.Model(), cur_resid.Chain(), res_name, cur_resid.ResNum(), cur_resid.AltLoc(), ConvertibleToAA( res_name))
            aStats["alternate_locations"][ res_id_str] = cur_resid.AltLoc()

        if cur_resid.IsPTMAminoAcid():
            aStats["n_ptm_aas"] += 1
            res_id_str = cur_resid.ShortStrDescription() #res_id_templ % (cur_resid.Model(), cur_resid.Chain(), res_name, cur_resid.ResNum(), cur_resid.AltLoc(), ConvertibleToAA( res_name))
            aStats["ptm_aas"][ res_id_str] = PTM_Type( res_name)
            #aStats["ptm_aas"].append( cur_resid.Chain(), cur_resid.ResNum())
        elif cur_resid.IsIon():
            aStats["ions"] += 1
            if len( cur_resid.ResName()) > 0 and res_name not in aStats["ion_names"]: aStats["ion_names"].append( res_name)
        elif cur_resid.IsDNA(): aStats["dna"] += 1




def ParseRemark( line, stats):

    try:
        if line[0:6] == "TITLE ":
            if len( stats["title"]) > 0:
                stats["title"] += " " #add space, title can be on multiple rows
            stats["title"] += line[10:].strip()
        elif line[0:6] == "KEYWDS":
            if len( stats["keywords"]): stats["keywords"] += " "
            stats["keywords"] += line[6:].strip()
        elif line[0:10] == "COMPND   2":
            stats["compound"] = re.sub( r';', '', line[10:]).strip()
            stats["compound"] = re.sub( r'MOLECULE:', '', stats["compound"]).strip()
        elif line[0:6] == "HEADER":
            m = re.match("(.*)\s+(\d{2}-.{3}-\d{2})\s+(.*)", line[6:])
            if m:
                stats["header"] = m.group(1).strip().capitalize()
                stats["published"] = m.group(2).strip().upper()
                stats["pdbcode"] = m.group(3).strip().upper()
        elif line[0:6] == "EXPDTA":
            stats["method"] = line[6:].strip()
        elif line.find("RESOLUTION.") > 0:
            m = re.match(".*(\d+\.\d+).*", line)
            if m: stats["resolution"] = float( m.group(1))
        elif line.find("FREE R VALUE    ") > 0:
            m = re.match(".*(\d+\.\d+).*", line)
            if m: stats["R-Free"] = float( m.group(1))
        elif line.find("3   R VALUE") > 0 and line.find("WORKING SET") > 0:
            m = re.match(".*(\d+\.\d+).*", line)
            if m: stats["R-Value"] = float( m.group(1))
        elif line.find("ORGANISM_SCIENTIFIC:") > 0:
            m = re.match(".*ORGANISM_SCIENTIFIC:\s*(.*)", line)
            if m:
                stats["organism_scientific"] = m.group(1)
                stats["organism_scientific"] = re.sub( r';', '', stats["organism_scientific"]).strip()
        elif line.find("ORGANISM_COMMON:") > 0:
            m = re.match(".*ORGANISM_COMMON:\s*(.*)", line)
            if m:
                stats["organism_common"] = m.group(1)
                stats["organism_common"] = re.sub( r';', '', stats["organism_common"]).strip()
        elif line.find("ORGANISM_TAXID:") > 0:
            m = re.match(".*ORGANISM_TAXID:\s*(\d+).*", line)
            if m: stats["organism_taxid"] = int(m.group(1))
        elif line.find( "350 BIOMOLECULE:") > 0:
            stats["n_biomolecules"] += 1


    except ValueError:
        raise


class Residue( object):

    def __init__( self):
        self.atoms = []
        #resindex is an unique running identifier for each residue within the pdb file
        #Files with multiple models can contain residues with identical residue sequence numbers ("A1809") and atom serial numbers ("ATOM    2405")
        #that should not be confused with this number. The resindex is not read from or written to the pdb file.
        self.resindex = -999999
        self.modelnum = -1
        self.firstloc = True

    def SetResIndex( self, index):
        self.resindex = index
        for a in self.atoms:
            a.resindex = index

    def SetModelnum( self, modelnum):
        self.modelnum = modelnum

    def AddAtom( self, atom):

        if atom.resindex >= 0: sys.stderr.write( "Warning: Atom %i already assigned to a residue.\n" % atom.index )

        self.atoms.append( atom)
        atom.resindex = self.resindex

    def NumberOfAtoms( self):
        return len( self.atoms)

    def OtherAtoms( self, atom):
        retval = []
        for a in self.atoms:
            if a != atom: retval.append( a)
        return retval

    def IsHETATOM( self):
      return False if len( self.atoms) == 0 else self.atoms[ 0].ishetatm

    def IsAminoAcid( self):
        return False if len( self.atoms) == 0 else self.atoms[ 0].isaminoacid

    def IsPTMAminoAcid( self):
        return False if len( self.atoms) == 0 else self.atoms[ 0].isptmaminoacid

    #Model;Chain;ResName;ResNum;Alternate location;Convertible to AA
    short_str_desc = "M:%s;C:%s;N:%s;R:%s;L:%s;A:%s;I:%i"

    def ShortStrDescription( self):
        return Residue.short_str_desc % (self.Model(), self.Chain(), self.ResName(), self.ResNum(), self.AltLoc(), ConvertibleToAA( self.ResName()), self.resindex)

    def HasMissingAtoms( self):
        global atoms_of_residue
        resname = self.ResName()
        #print "RESNAME:", self.ResName(), "ATOMS:%i/%i" % (len( self.atoms), len( atoms_of_residue[ resname]))
        if resname not in atoms_of_residue: return None
        return True if len( self.atoms) < len( atoms_of_residue[ resname]) else False

    def HasSideChain( self):
        return True if len( self.atoms) > 4 else False

    def IsIon( self):
        return False if len( self.atoms) == 0 else self.atoms[ 0].ision

    def IsDNA( self):
        return False if len( self.atoms) == 0 else self.atoms[ 0].isdna

    def IsWater( self):
        return False if len( self.atoms) == 0 else self.atoms[ 0].iswater

    def AltLoc( self):
        return "" if len( self.atoms) == 0 else self.atoms[ 0].AltLocationIndicator()

    def HasAltLocation( self):
        for a in self.atoms:
            if a.HasAltLocation(): return True
        return False

    def IsPrimaryLocation( self):
        return self.firstloc

    def IsAlternateLocation( self):
        return not self.firstloc

    def AtomType( self):
        if   self.IsAminoAcid(): return "A"
        elif self.IsIon():       return "I"
        elif self.IsWater():     return "W"
        elif self.IsDNA():       return "D"
        elif self.atoms[ 0].isother: return "O"
        return "U" #unknown

        return False if len( self.atoms) == 0 else self.atoms[ 0].iswater

    def ResName( self):
        return "" if len( self.atoms) == 0 else self.atoms[ 0].resname

    def OneLetterResName( self):
        global res_map
        if len( self.atoms) == 0: return ""
        if self.IsAminoAcid(): return res_map[ self.ResName()]
        elif self.IsDNA(): return pdb_dna[ self.ResName()]
        return "."

    #WARNING: return value is str and can include insertion codes
    def ResNum( self):
        return "" if len( self.atoms) == 0 else ("%i%s" % (self.atoms[ 0].resnum, self.atoms[ 0].Insertion()))

    #Excludes insertion code if it exists (integer). Note: Not the same thing as resindex
    def ResSeqNum( self):
        return 0 if len( self.atoms) == 0 else self.atoms[ 0].resnum

    def InsertionCode( self):
        return "" if len( self.atoms) == 0 else ("%s" % self.atoms[ 0].Insertion())

    def Chain( self):
        return " " if len( self.atoms) == 0 else self.atoms[ 0].chain

    def Model( self):
        return str( self.modelnum)

    def ModelNum( self):
        return self.modelnum

    def CoordsCA( self):
        for a in self.atoms:
            if a.atomname == "CA": return a.coord
        return [0.0, 0.0, 0.0]

    @staticmethod
    def CoordDistance( coord1, coord2):
        dx = coord2[ 0]-coord1[ 0]
        dy = coord2[ 1]-coord1[ 1]
        dz = coord2[ 2]-coord1[ 2]
        return math.sqrt( dx*dx+dy*dy+dz*dz )

    @staticmethod
    def ResetBCols( aResidueList):
        for r in aResidueList: r.RestoreOriginalAtomBColValues()

    def BColAvg( self):
        retval = 0.0

        if len( self.atoms) == 0: return retval #Avoid division by zero
        for a in self.atoms: retval += a.b
        return (retval / float( len( self.atoms)))

    def SetBCol( self, val):
        for a in self.atoms:
            a.b = val

    def ResetBCol( self):
        for a in self.atoms:
            a.b = a.original_bcol

    def FirstAtomSerialNumber( self):
        return " " if len( self.atoms) == 0 else self.atoms[ 0].index

    def __str__( self):
        retval = ""
        for a in self.atoms:
            retval += a.Line() + "\n"
        return retval

    def __eq__( self, other): # == (Same chain, model, everything)
        if other: return self.resindex == other.resindex
        return False

    def __ne__( self, other): # !=
        return not self.__eq__( other)

    #File order
    def __lt__( self, other): # <
        if other: return self.resindex < other.resindex
        return False

    def __gt__( self, other): # >
        if other: return self.resindex > other.resindex
        return False

    def __le__( self, other): # <=
        if other: return self.resindex <= other.resindex
        return False

    def __ge__( self, other): # >=
        if other: return self.resindex >= other.resindex
        return False

    def __contains__(self, atom):
        return atom in self.atoms

    def __hash__( self):
        return hash( self.resindex)

    def IsMatchingResidue( self, other): #Matching residue in another chain and or model

        if other and type( other) == Residue and self.ResNum() == other.ResNum() and self.ResName() == other.ResName(): return True
        elif other and type( other) == list:
            for o in other:
                if self.IsMatchingResidue( o): return True
        return False

    def PrintCA( self):
        for a in self.atoms:
            if a.atomname == "CA":
              print a.Line()
              break
        return ""

    def Print( self):
        if len( self.atoms) == 0: print "Empty residue"
        for a in self.atoms: print a.Line()


    def SetTerminal( self, aIsTerminal):
        if len( self.atoms) > 0: self.atoms[ -1].terminal = True

    def RestoreOriginalAtomBColValues( self):
        for a in self.atoms: a.RestoreOriginalBcolValue()

    def ConvertToAAIfPossible( self, aRemovePTM=False):

        n_atoms_converted = 0
        n_atoms_removed = 0
        for a in reversed( self.atoms):

            converted = a.ConvertPTMtoATOMentry()
            if converted: n_atoms_converted += 1

            #if converted: print "Converted:", a

            if converted and aRemovePTM:
                resname = a.resname

                if a.atomname == "SE" and resname == "MET": a.atomname = "SD" #Selenomethionine: Selenine -> Sulfur

                #print "Resname", resname, ":", resname in atoms_of_residue.keys()
                elif resname in atoms_of_residue.keys() and a.atomname not in atoms_of_residue[ resname]:
                    self.atoms.remove( a)
                    n_atoms_removed += 1
                    n_atoms_converted -= 1 #was converted, now removed
                    #print "Removed:", a

        return n_atoms_converted, n_atoms_removed

    @staticmethod
    def TranslateCoord( aMatrix, aCoord):

        ret_vec = [0.0,0.0,0.0]

        for p in range( 3):
            ret_vec[ p] =  aMatrix[ p][ 0]*aCoord[ 0] + aMatrix[ p][ 1]*aCoord[ 1] + aMatrix[ p][ 2]*aCoord[ 2]
            ret_vec[ p] += aMatrix[ p][ 3]

        return ret_vec


    def Translate( self, aMatrix):

        for atom in self.atoms: atom.coord = Residue.TranslateCoord( aMatrix, atom.coord)


class Atom:

    backbone_atom_names = ['C','N','O','CA']

    def __init__( self, line, atom_id, modelnum=-1):

        global pdb_aas, pdb_dna, pdb_ion, pdb_water, pdb_other, VDW_radii

        l = len( line)
        if l < 54:
            raise ValueError("PDB line is too short.")

        self.atom_id = atom_id
        self.terminal = False
        self.first_in_model = -1
        self.modelnum = modelnum

        #mandatory
        self.atom = line[0:6]              #1-6
        self.index = int(line[6:11])       #7-11
        self.atomname = line[12:16].strip()#13-16
        self.locind = line[16]             #17
        self.resname = line[17:20].strip() #18-20
        self.chain = line[21]              #22
        self.resnum = int(line[22:26])     #23-26
        self.insert = line[26]             #27
        self.coord = [float(line[30:38]),  #31-38
                      float(line[38:46]),  #39-46
                      float(line[46:54])]  #47-54
        #optional
        self.occ = 100.0
        self.b = 0.0
        self.segid = ""
        self.elsym = ""
        self.charge = ""

        if l >= 60: self.occ = float(line[54:60])    #55-60
        if l >= 66: self.b = float(line[60:66])      #61-66
        if l >= 76: self.segid = line[72:76].strip() #73-76
        if l >= 78: self.elsym = line[76:78].strip() #77-78
        if l >= 80: self.charge = line[78:80].strip()#79-80

        #additional
        self.original_bcol = self.b

        self.isatom = line[0:6] == "ATOM  "
        self.ishetatm = line[0:6] == "HETATM"

        self.ishydrogen = self.atomname[ 0] == "H"

        self.isaminoacid = IsAA( self.resname )
        self.isptmaminoacid = IsPTMAA( self.resname ) #Post-translationally modified AA + Selenomethionine
        self.iswater = self.resname in pdb_water
        self.ision = self.resname in pdb_ion
        self.isother = self.resname in pdb_other

        #if self.ision: print "'%s'" % self.resname

        self.isdna = self.resname in pdb_dna

        self.isunknown = not (self.ishydrogen or self.iswater or self.isaminoacid or self.ision or self.isdna or self.isother or self.isptmaminoacid)

        #self.box = [-1,-1,-1]

        self.asa = 0.0 #Accessible surface area
        self.sa = 0.0  #Surface area (area that does not collide withing the residue)
        self.pasa = 0.0  #Percentage of the surface area that is solvent accessible

        self.resindex = -999999

        #Set VDW radius for atom
        self.VDW_type = "C"
        self.hasradius = False
        self.radius = 0.0

        self.SetVDW()
        #self.SetIndex()


    def __eq__(self, other):
        return self.atom_id == other.atom_id

    def __ne__(self, other):
        return self.atom_id != other.atom_id

    def __str__( self):
        return self.Line()

    def ShortLineFormat( self):
        #"12: CB SER C 161 "
        return "%04i: %3s%s %s %s %s" % (self.resindex, self.atomname, self.locind, self.resname, self.chain, self.resnum)

    def SameModel( self, other):
      return (self.modelnum == other.modelnum)

    def SameResidue( self, other):
        return (self.resindex == other.resindex)

    def SetIndex( self, aIndex=None): #Max space for index is 5 chars
        if aIndex: self.index = aIndex
        if self.index > 99999: self.index %= 100000

    def Insertion( self):
        return (self.insert if self.insert.isalpha() else "")

    def HasAltLocation( self):
        return self.locind != " "

    def AltLocationIndicator( self):
        return self.locind

    def IsBackBone( self):
        return self.isaminoacid and (self.atomname in Atom.backbone_atom_names)

    def IsBBorSideChain( self ):
        global atoms_of_residue
        return (self.resname in atoms_of_residue and self.atomname in atoms_of_residue[ self.resname])

    def FileDescription( self):

        model_str = ""
        ter_str = ""
        model_end_str = ""

        if self.first_in_model >= 0:
            model_str =     "MODEL%9i                                                                  \n" % self.first_in_model

        if self.terminal:
            ter_str =       "TER   %5i      %-3s %s%4i                                                      \n" % (self.index+1,self.resname[:3],self.chain[:1],self.resnum)

        if self.modelnum >= 0 and self.terminal:
            model_end_str = "ENDMDL                                                                          \n"

        return (model_str + self.Line() + "\n" + ter_str + model_end_str)


    def IsOfIncludedType( self, included_types ):

        acc = included_types

        if self.isaminoacid and "A" in acc: return True
        if self.isdna       and "D" in acc: return True
        if self.isunknown   and "U" in acc: return True
        if self.ision       and "I" in acc: return True
        if self.isother     and "O" in acc: return True

        return False

    #Post-translational modifications
    def ConvertPTMtoATOMentry( self):

        global PTM_aas

        if self.isptmaminoacid and not self.isaminoacid:
            self.resname = PTM_aas[ self.resname]
            self.ishetatm = False
            self.isaminoacid = True
            return True
        return False

    def RestoreOriginalBcolValue( self):
        self.b = self.original_bcol

    def SetVDW( self):

        global pdb_aas, pdb_dna, pdb_ion, pdb_water, pdb_other, VDW_radii

        if self.ision:
            #ION
            if self.atomname in VDW_radii and self.resname in VDW_radii:
                self.radius = VDW_radii[ self.atomname]
                self.hasradius = True
                self.VDW_type = self.atomname
        else:
            #Not an ion
            if len( self.elsym) and self.elsym in VDW_radii and self.elsym not in pdb_ion:
                self.radius = VDW_radii[ self.elsym]
                self.hasradius = True
                self.VDW_type = self.elsym
            elif self.atomname in VDW_radii and self.atomname not in pdb_ion:
                self.radius = VDW_radii[ self.atomname]
                self.hasradius = True
                self.VDW_type = self.atomname
            elif len( self.atomname) and self.atomname[ 0] in VDW_radii and self.atomname[ 0] not in pdb_ion:
                self.radius = VDW_radii[ self.atomname[ 0]]
                self.hasradius = True
                self.VDW_type = self.atomname[ 0]

            else:
                no_digit_name = ''.join([i for i in self.atomname if i.isalpha()])

                if len( no_digit_name) and no_digit_name in VDW_radii and no_digit_name not in pdb_ion:
                    self.radius = VDW_radii[ no_digit_name]
                    self.hasradius = True
                    self.VDW_type = no_digit_name
                elif len( no_digit_name) and no_digit_name[ 0] in VDW_radii and no_digit_name[ 0] not in pdb_ion:
                    self.radius = VDW_radii[ no_digit_name[ 0]]
                    self.hasradius = True
                    self.VDW_type = no_digit_name[ 0]


    def Line( self):
        line =  "%-6s%5i %-4s%s%-3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f       %-4s%-2s%-2s" % (
                "ATOM" if not self.ishetatm else "HETATM", #1-6
                self.index,      #7-11
                self.atomname[:4] if len( self.atomname) >= 4 else " " + self.atomname,  #13-16
                self.locind,     #17
                self.resname[:3],#18-20
                self.chain[:1],  #22
                self.resnum,     #23-26
                self.insert[:1], #27
                self.coord[ 0],
                self.coord[ 1],
                self.coord[ 2],  #31-38,39-46,47-54
                self.occ,        #55-60
                self.b,          #61-66
                self.segid,      #73-76
                self.elsym,      #77-78
                self.charge)     #79-80

        return line


#Returns list of biomolecules specified with REMARK 350
#list (each biomolecule) contains dict with keys "num", list:"chains", list:"matrices"
def ReadBiologicalAssemblies( aRemarks):

    matrices = []
    cur_matrix = -1
    cur_chains = []
    cur_molec = -1
    r350 = "REMARK 350 "
    biomolecules = []


    for remark in aRemarks:

        if not remark.startswith( r350): continue

        remark_data = remark[ len( r350):]
        cols = remark_data.split()
        if len( cols) == 0: continue

        if cols[ 0].startswith( "BIOMOLECULE"):
            try:
                pos = remark.rfind( ":")
                biomol_str = remark[ pos+1:].strip()
                cur_molec = int( biomol_str.split(",")[ 0]) #Sometimes more than one molecule "2XDV"
                biomolecules.append( { "num":cur_molec, "chains":[], "matrices":[] })
                cur_matrix = -1
            except Exception as ex:
                sys.stderr.write( "PDB ERROR: Could not interpret Biomolecule information from line:\n")
                sys.stderr.write( "          '%s'\n" % remark)
                sys.stderr.write( "           %s\n" % str( ex))

        elif remark_data.startswith("APPLY THE FOLLOWING TO CHAINS"):

            try:
                pos = remark.rfind( ":")
                cur_chains = map( str.strip, remark[ pos+1:].split(","))
                cur_chains = filter( None, cur_chains)
                if len( biomolecules): biomolecules[ -1]["chains"] = cur_chains
                #print "CHAINS: %s" % ",".join( cur_chains)
            except Exception as ex:
                sys.stderr.write( "PDB ERROR: Could not interpret chains information from line:\n")
                sys.stderr.write( "          '%s'" % remark)
                sys.stderr.write( "           %s" % str( ex))
        elif remark_data.find("       AND CHAINS:") >= 0:
            try:
                pos = remark.rfind( ":")
                cur_chains = map( str.strip, remark[ pos+1:].split(","))
                cur_chains = filter( None, cur_chains)
                if len( biomolecules): biomolecules[ -1]["chains"].extend(  cur_chains)
                #print "CHAINS: %s" % ",".join( cur_chains)
            except Exception as ex:
                sys.stderr.write( "PDB ERROR: Could not interpret chains information from line:\n")
                sys.stderr.write( "          '%s'" % remark)
                sys.stderr.write( "           %s" % str( ex))
        elif cols[ 0].startswith("BIOMT"):
            #REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
            matrix_num = cols[ 1]
            if matrix_num != cur_matrix:
                # new transformation matrix
                #print "Curmolec:", cur_molec
                if len( biomolecules): biomolecules[ -1]["matrices"].append([])
                #print biomolecules

            #matrices = biomolecules[ -1]["matrices"]
            #last_matrix = matrices[ -1]
            if len( biomolecules): biomolecules[ -1]["matrices"][-1].append( map( float, cols[ 2:]))
            cur_matrix = matrix_num

    return biomolecules

def NewPdb():
    new_pdb =  {"residues":[],"atoms":[],"remarks":[],"stats":{}}
    new_pdb["stats"]["atom_types"]  = {}
    new_pdb["stats"]["amino_acids"] = 0
    new_pdb["stats"]["largest_radius"]  = 1.2
    new_pdb["stats"]["filename"] = "new"
    new_pdb["stats"]["pdbcode"] = "Unknown"
    new_pdb["stats"]["incomplete_residue_ratio"] = 0.0
    new_pdb["stats"]["single_atom_residue_ratio"] = 0.0
    new_pdb["stats"]["backbone_only_ratio"] = 0.0
    return new_pdb


def UpdateStatsAtomRatios( aStats, aResid):

    incomplete_res_val = 1.0 if aResid.HasMissingAtoms() else 0.0
    cur_n_atoms = aResid.NumberOfAtoms()
    singleatom_res_val = 1.0 if cur_n_atoms == 1 else 0.0
    backbone_res_val = 1.0 if cur_n_atoms <= 4 else 0.0
    aStats["incomplete_residue_ratio"] = (aStats["incomplete_residue_ratio"] * (aStats["amino_acids"]-1) + incomplete_res_val) / aStats["amino_acids"]
    aStats["single_atom_residue_ratio"] = (aStats["single_atom_residue_ratio"] * (aStats["amino_acids"]-1) + singleatom_res_val) / aStats["amino_acids"]
    aStats["backbone_only_ratio"] = (aStats["backbone_only_ratio"] * (aStats["amino_acids"]-1) + backbone_res_val) / aStats["amino_acids"]


#aChains, write only biomolecules that contain one of these chains
def WriteBiologicalAssembliesToFiles( aPdb, aOutputFolder, aChains=[], aBiomolecules=[], aQuiet=False, aCenter=[]):

    #def ReadPDB( infile, chains=[], models=[], firstmodelonly=False, aQuiet=False, aConvertPTMs=False, aWarnUnknown=True ):
    pdb = aPdb#ReadPDB( aInputFile, aChains, firstmodelonly=True, aConvertPTMs=False, aWarnUnknown=False, aQuiet=aQuiet)
    assemblies_written = 0


    assemblies = ReadBiologicalAssemblies( pdb["remarks"])

    if len( assemblies) == 0:
        if not aQuiet: sys.stderr.write( "INFO: No biological assemblies found in file '%s'.\n" % pdb["stats"]["filename"])
        return 0

    if not aQuiet: sys.stderr.write( "INFO: %i biological assemblies found in file '%s'.\n" % (len( assemblies), pdb["stats"]["filename"]))

    for assembly in assemblies:

        if aBiomolecules and len( aBiomolecules) > 0 and assembly["num"] not in aBiomolecules: continue
        if aChains and len( aChains):
            common_chains = list( set(aChains).intersection( assembly["chains"]))
            if len( common_chains) == 0: continue

        bio_pdb = NewPdb()
        CreateBiologialAssembly( assembly, aPdbSource=pdb, aPdbTarget=bio_pdb, aQuiet=aQuiet, aCenterResidue=aCenter)

        bio_pdb["remarks"] = pdb["remarks"]

        input_file = pdb["stats"]["filename"]
        if input_file.find( ".pdb") < 0 and input_file.find( ".ent") < 0 and "pdbcode" in pdb["stats"]: input_file = pdb["stats"]["pdbcode"]


        outputfile = pdb["stats"]["filename"].replace( ".pdb", ".biomol%i.pdb" % assembly["num"] )
        outputfile = pdb["stats"]["filename"].replace( ".ent", ".biomol%i.ent" % assembly["num"] )

        if outputfile.find( ".pdb") < 0 and outputfile.find( ".ent") < 0: output_file += ".biomol%i.pdb" % assembly["num"]

        if not aOutputFolder or len( aOutputFolder) == 0: aOutputFolder = "./"
        elif len( aOutputFolder) and aOutputFolder[ -1] != "/" and aOutputFolder[ -1] != "\\": aOutputFolder += "/"

        #in_path, in_filename = os.path.split( input_file)
        out_path, out_filename = os.path.split( outputfile)
        outputfile = aOutputFolder + out_filename

        try:
            remark_source = pdb["stats"]["pdbcode"] if "pdbcode" in pdb["stats"] else pdb["stats"]["filename"]
            with open( outputfile, "w" ) as pdb_fh:
                if not aQuiet: print "PDB INFO: Writing assembly '%s' to file '%s'." % (assembly["num"], outputfile)
                #pdb_fh.write( "".join( pdb["remarks"]))
                pdb_fh.write( "REMARK 300 BIOLOGICAL ASSEMBLY '%s' OF %s\n" % ( assembly["num"], pdb["stats"]["pdbcode"]))
                WriteToFile( pdb_fh, bio_pdb["atoms"], types=["A","D"], chains=aChains)

        except IOError as e:
            sys.stderr.write( "ERROR: Unable to write pdb file: '%s' Error: %s\n" % ( outputfile, e.strerror))
            return False

        except Exception as ex:
            sys.stderr.write( "ERROR: Unable to write PDB file: '%s'.\n" % ( outputfile))
            template =        "       An exception of type {0} occured. Arguments: {1!r}"
            message = template.format(type(ex).__name__, ex.args)
            sys.stderr.write( message + "\n")

        assemblies_written += 1


    if not aQuiet and len( aChains): sys.stderr.write( "INFO: No biological assemblies found for chains '%s' in file '%s'." % (",".join( aChains), pdb["stats"]["filename"]))
    return assemblies_written

def IsRestricted( aRemarks):
    for remark in aRemarks:
        if remark.startswith( u"REMARK 300 Construct restricted"): return True
    return False

#aCenter will be translated to biological assembly coordinates
#With aRestrictAboveAtoms large biological assemblies can be limited to contain only residues within aDistance of ACenter
#Assemblies with less than aRestrictAboveAtoms will not be limited
#If the assembly has been limited, the file will contain a "REMARK 300" stating this
def CreateBiologialAssembly( aAssembly, aPdbSource, aPdbTarget=None, aQuiet=False, aCenterResidue=None, aRestrictAboveAtoms=-1, aDistance=200.0):

    if aPdbTarget == None: aPdbTarget = aPdbSource

    #aPdbTarget["stats"]["atom_types"] = aPdbSource["stats"]["atom_types"].copy()
    #aPdbTarget["stats"]["largest_radius"] = aPdbSource["stats"]["largest_radius"]

    #print "KEYS:", aPdbSource.keys()
    if not aQuiet:
        pdb_code = "Unknown"
        if "pdbcode" in aPdbSource["stats"]: pdb_code = aPdbSource["stats"]["pdbcode"]
        print "PDB INFO: Creating biological assembly '%s' for %s." % (aAssembly["num"], pdb_code)

    translation_matrices = aAssembly["matrices"]
    chains = sorted( list( set( aAssembly["chains"]))) #Make sure no duplicates

    #Note: It is unclear what 1QKU biomol definitions mean
    #"Chains from A to Z", like in 1QKU
    #if len( chains) == 2 and chains[ 0] == "A" and chains[ 1] == "Z" and "Z" not in aPdbSource["stats"]["chains"]:
    #    chains = aPdbSource["stats"]["chains"]
    #else:

    #Check that all chains for biol assembly actually exist. I'm looking at you, 1QKU >:|
    for c in chains:
        if c not in aPdbSource["stats"]["chains"]: sys.stderr.write("PDB WARNING: Biological assembly chain '%s' not found in file.\n" % c)
    chains = list( set( chains).intersection( set(aPdbSource["stats"]["chains"])))


    if len( translation_matrices) == 0: sys.stderr.write("PDB ERROR: No translation matrices specified for assembly.\n")
    if len( chains) == 0: sys.stderr.write("PDB ERROR: No chains specified for assembly.\n")

    #Prevent re-copying copied residues
    last_res_id = -1
    if aPdbTarget == aPdbSource: last_res_id = aPdbSource["residues"][ -1].resindex

    restrict_around = None
    number_of_atoms_to_generate = 0

    for chain in chains: number_of_atoms_to_generate += aPdbSource["stats"]["atoms_per_chain"][ chain]
    number_of_atoms_to_generate *= len( translation_matrices)
    if not aQuiet: print "INFO: Generating %i atoms for assembly." % number_of_atoms_to_generate

    if aCenterResidue and aRestrictAboveAtoms >= 0 and number_of_atoms_to_generate > aRestrictAboveAtoms:
        center = aCenterResidue.CoordsCA()
        #Use first matrix for translation of coords
        restrict_around = Residue.TranslateCoord( translation_matrices[ 0], center)
        aPdbTarget["remarks"].append( u"REMARK 300 Construct restricted within %.1fÅ of [%.2f, %.2f, %.2f]\n" % tuple([aDistance]+center))
        if not aQuiet: print u"Restricting construct within %.1fÅ of [%.2f, %.2f, %.2f]" % tuple([aDistance]+center)

    n_residues = 0

    for chain in chains:
        #print "Processing CHAIN:", chain
        for matrix in translation_matrices:
            #print "Processing MTX:", matrix
            for residue in aPdbSource["residues"]:

                #Prevent re-copying copied residues
                if last_res_id >= 0 and residue.resindex > last_res_id: break

                if residue.Chain() == chain:

                    if restrict_around:
                        #Calc CA distance
                        distance = Residue.CoordDistance( restrict_around, Residue.TranslateCoord( matrix, residue.CoordsCA()))
                        if distance > aDistance: continue #Skip residue

                    rescopy = CreateCopyOfResidue( aPdbSource, residue, aPdbTarget)
                    rescopy.Translate( matrix)
                    n_residues += 1

    if not aQuiet: print "INFO: Generated %i residues" % n_residues

#TODO order in aAtomList
def SplitResidueToAlternateLocations( aResidue, aAtomList ):


    alt_locations = []
    for atom in aResidue.atoms:
        if atom.locind != " " and atom.locind not in alt_locations: alt_locations.append( atom.locind)

    #alt_locations.append( " ") #Some residues contain common atoms with each other
    if len( alt_locations) <= 1: return [aResidue]
    alt_locations.append( alt_locations.pop( 0)) #primary location to last

    cur_residue = aResidue
    retvals = [aResidue]

    #atoms with loc_index " " stay with aResidue
    #Copies created for others

    for cur_location in alt_locations:

        #print "Processing Location:", cur_location
        primary_location = (cur_location == alt_locations[ -1])

        if not primary_location:
            cur_residue = Residue() #New residue
            #Residue index set later
            cur_residue.SetModelnum( aResidue.atoms[ 0].modelnum)
            cur_residue.firstloc = False
            retvals.append( cur_residue)
        #else:
        #    print "PRIMARY LOCATION"

        for ai in reversed( range( len( aResidue.atoms))):

            cur_atom = aResidue.atoms[ ai]

            if primary_location:
                #Last location in loop
                if cur_atom.locind == " ": cur_atom.locind = cur_location
            else:

                #Shared atoms between alternate locations
                if cur_atom.locind == " ":
                    #Create copy of atom
                    #print "Copying atom '%s' -> '%s': '%s'" % (cur_atom.locind, cur_location, cur_atom.Line())
                    new_atom = CopyAtom( cur_atom) #New atom with new index
                    new_atom.locind = cur_location
                    cur_residue.AddAtom( new_atom)
                    #aAtomList.append( new_atom)
                #Non-shared Atoms of alternate locations
                elif cur_atom.locind == cur_location:
                    #print "Transferring atom %s -> %s: '%s'" % (cur_atom.locind, cur_location, cur_atom.Line())
                    #Transfer atom from original residue
                    cur_atom.resindex = -1 #Reset
                    cur_residue.AddAtom( cur_atom)
                    #Remove from original
                    #aResidue.atoms.remove( cur_atom)
                    del aResidue.atoms[ ai]

        if not primary_location: cur_residue.atoms = list( reversed( cur_residue.atoms))


    #Atom indices set when residue added
    #print "aResidue:"
    #aResidue.Print()
    #cur_atom_index = len( aAtomList)+1
    #for r in retvals:
    #    for a in r.atoms:
    #        a.atom_id = cur_atom_index
    #        cur_atom_index += 1

    return retvals


#resindex not copied
def CopyAtom( aSourceAtom, aNewAtomId=-1 ):

    new_atom = Atom( aSourceAtom.Line(), aNewAtomId, aSourceAtom.modelnum)
    new_atom.index = aSourceAtom.index #Keep index

    #Copy values not directly on pdb line
    new_atom.first_in_model = aSourceAtom.first_in_model
    new_atom.terminal = aSourceAtom.terminal
    new_atom.original_bcol = aSourceAtom.original_bcol
    new_atom.asa = aSourceAtom.asa
    new_atom.sa = aSourceAtom.sa
    new_atom.pasa = aSourceAtom.pasa

    return new_atom


def CreateCopyOfResidue( aPdbSource, aResidue, aPdbTarget=None ):

    if not aResidue or aResidue.NumberOfAtoms() == 0  or len( aPdbSource["residues"]) == 0: return None

    #assert( aPdbSource and aPdbSource != None)
    if aPdbTarget == None: aPdbTarget = aPdbSource

    stats = aPdbTarget["stats"]
    n_target_pdb_atoms = len( aPdbTarget["atoms"])
    new_atoms = []
    pdb_running_atom_id = 1 if n_target_pdb_atoms == 0 else (aPdbTarget["atoms"][ -1].atom_id + 1)
    pdb_running_atom_index = 1 if n_target_pdb_atoms == 0 else (aPdbTarget["atoms"][ -1].index + 1)

    new_residue = Residue()
    new_residue.SetResIndex( len( aPdbTarget["residues"]))
    new_residue.SetModelnum( aResidue.atoms[ 0].modelnum)
    new_residue.firstloc = aResidue.firstloc

    for atom in aResidue.atoms:

        #Create atom with new atom_id and index number
        new_atom = Atom( atom.Line(), pdb_running_atom_id, atom.modelnum)
        new_atom.index = pdb_running_atom_index #Might clash with HOH residues after all aas in file

        #Copy values not directly on pdb line
        new_atom.first_in_model = atom.first_in_model
        new_atom.terminal = atom.terminal
        new_atom.original_bcol = atom.original_bcol
        new_atom.asa = atom.asa
        new_atom.sa = atom.sa
        new_atom.pasa = atom.pasa


        aPdbTarget["atoms"].append( new_atom)

        new_residue.AddAtom( new_atom)

        if new_atom.radius > stats["largest_radius"]: stats["largest_radius"] = new_atom.radius
        if new_atom.VDW_type not in stats["atom_types"]: stats["atom_types"][ new_atom.VDW_type] = 0
        stats["atom_types"][ new_atom.VDW_type] += 1


        pdb_running_atom_id += 1
        pdb_running_atom_index += 1
        if pdb_running_atom_index > 99999: pdb_running_atom_index = 1 #PDB file format only allows 5 chars for index

    aPdbTarget["residues"].append( new_residue)

    if new_residue.IsAminoAcid():

        stats["amino_acids"] += 1
        UpdateStatsAtomRatios( stats, new_residue)

    return new_residue


def GetVersion():
  global pdbatoms_version
  return pdbatoms_version


def help():

    print """
    This script reads PDB files and outputs the selected chains and models
    Usage: %s [-c CHAINS][-m MODELS][-p MODE][-f][-s][-h][-b][-i] input_file [output_file/folder]

    -c  --chains  CHAINS     Chain selection [default:ALL] (separate with commas)
                             Double capital letter chains are interpreted as lower
                             case chains: "AA" -> "a".
    -m  --models  MODELS     Model selection [default:ALL] (separate with commas)
                             If the file does not have any MODEL entries this
                             flag is ignored.
    -a  --atomtypes [AWDUIO] Included atom types as a string of letters i.e. "ADI"
                              A : amino acids, W : water, D : DNA/RNA
                              U : unknown, I : ion, O : other
    -f  --firstmodelonly     Select only the first model in file [default:False]
                             (--models flag will be ignored)
    -i  --info               Print file information (statistics)
    -k  --infokeys           Print only selected file information (separate with commas)
                             Use -i to see full list of keys
    -p  --ptms [N|C|A]       Post-translational modification conversions (+selenomethionine)
                             Modes: N = No conversions [default]
                                    C = Convert to regular amino acids
                                    A = Atomize, only change PTM atoms from HETATM to ATOM
    -b  --biomolecules       Write biological assemblies in input_file to output_folder.
                             If no output_folder in specified output is written to
                             current folder. Use the chains argument to write only molecules
                             for specific chains. When combined with the --info flag available
                             biomolecules are only listed.
    -q  --quiet              Less output
    -s  --sequence           Extract sequence (Default atom type selection is "ADI")
                             Other than atoms of type A and D are outputted as "."
    -g  --gapped_sequence    Extract sequence with gaps at missing residues (numbering)
    -h  --help               Print this message.

    Version: %s
    """ % (os.path.basename(__file__), GetVersion())



#Standalone usage
def main():

    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c:m:a:k:p:fihsgqb', ['chains', 'models', 'firstmodel','info', 'help','sequence','gapped_sequence', 'biomolecules'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)

    inputfile = ""
    outputfile = ""
    firstmodelonly = False
    print_stats = False
    statkeys = []
    chains = []
    models = []
    atomtypes = []
    sequence_selection = 0
    quiet = False
    conversions = "N"
    write_biomols = False

    #Flags
    for opt, arg in opts:
        #print "OPT", opt
        if opt in ('-f', '--firstmodel'):
            firstmodelonly = True
        elif opt in ('-i', '--info'):
            print_stats = True
        elif opt in ('-s', '--sequence'):
            sequence_selection = 1
        elif opt in ('-g', '--gapped_sequence'):
            sequence_selection = 2
        elif opt in ('-q', '--quiet'):
            quiet = True
        elif opt in ('-b', '--biomolecules'):
            write_biomols = True
        elif opt in ('-h', '--help'):
            help()
            sys.exit( 0)
        elif opt in ('-c', '--chains'):
            chains = arg.replace('"','').split(",")
        elif opt in ('-p', '--ptms'):
            conversions = arg
        elif opt in ('-m', '--models'):
            models = arg.replace('"','').split(",")
        elif opt in ('-k', '--infokeys'):
            statkeys = arg.replace('"','').split(",")
        elif opt in ('-a', '--atomtypes'):
            atomtypes = arg.replace('"','').split(",")
            if len( atomtypes) == 1: atomtypes = list( arg.replace('"',''))
        else:
            sys.stderr.write("WARNING: Unknown option ignored: '%s'.\n" % opt)

    #Input & Output files
    for arg in args:
        if len( inputfile) == 0:
            inputfile = arg
        elif len( outputfile) == 0:
            outputfile = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)

    if not len( inputfile):
        sys.stderr.write("No input file specified.)" )
        help()
        sys.exit( 0)

    if not os.path.isfile( inputfile):
        sys.stderr.write("ERROR: '%s' does not specify a valid file.\n" % inputfile)
        sys.exit( -2)

    #{ "atoms":atoms, "residues":residues, "stats":stats }
    info = ReadPDB( inputfile, chains, models, firstmodelonly, aQuiet=quiet, aConvertPTMs=conversions)

    #Warnings:
    for c in chains:
        if c not in info["stats"]["chains"]:
            sys.stderr.write("WARNING: Chain '%s' was not found.\n" % c)

    for m in models:
        if m not in info["stats"]["models"]:
            if not firstmodelonly: sys.stderr.write("WARNING: Model '%s' was not found.\n" % m)

    if write_biomols:

        #if len( outputfile) == 0:
            #sys.stderr.write( "ERROR: No output folder specified.\n" )
            #outputfile = os.path.split( inputfile)[ 0]

        #-b -i
        if print_stats:
            assemblies = ReadBiologicalAssemblies( info["remarks"])
            if len( assemblies) == 0:
                sys.stdout.write( "INFO: No biological assemblies found in file '%s'.\n" % info["stats"]["filename"])
            else:
                sys.stdout.write( "INFO: %i biological assembl%s found in file '%s'.\n" % (len( assemblies), ("ies" if len(assemblies) != 1 else "y"), info["stats"]["filename"]))
                sys.stdout.write( "INFO: Available assemblies: %s.\n" % ",".join([str(x["num"]) for x in assemblies]))
                if len( chains) > 0:
                    available_assemblies = []
                    for c in chains:
                        for a in assemblies:
                            common_chains = list( set( chains).intersection( a["chains"]))
                            if len( common_chains) > 0:
                                available_assemblies.append( str(a["num"]))
                    available_assemblies = list( set( available_assemblies))
                    if not len( available_assemblies): available_assemblies.append( "None")
                    sys.stdout.write( "INFO: Available assembl%s for chains %s: %s.\n" % ( ("ies" if len(assemblies) != 1 else "y"), ",".join( chains), ",".join(available_assemblies)))
            sys.exit( 0)


        if len( outputfile) == 0:
            outputfile = "./"

        if not os.path.isdir( outputfile):
            sys.stderr.write( "ERROR: '%s' does not specify a valid folder.\n" % outputfile)

        written = WriteBiologicalAssembliesToFiles( info, outputfile, aChains=chains, aBiomolecules=[], aQuiet=quiet)
        if not quiet: print "Wrote %i biological assembly file%s." % (written, ("s" if written != 1 else ""))
        sys.exit( 0)

    outfile_handle = None
    if len( outputfile):
      try:
        outfile_handle = open( outputfile, "w+")
      except Exception as ex:
          sys.stderr.write( "ERROR: Error writing output file file: '%s'.\n" % (cur_file))
          sys.stderr.write( "       An exception of type {0} occured. Arguments: {1!r}\n".format( type( ex).__name__, ex.args))
          sys.exit( -1)
    else:
        outfile_handle = sys.stdout #Print to stdout

    fasta_seq = ""
    if sequence_selection > 0:
        #Print sequence
        if len(atomtypes) == 0: atomtypes = ["A","D","I"] #default
        row_len = 60
        if print_stats or len( statkeys): row_len = 0 #do not split
        if sequence_selection == 1: fasta_seq = FastaSequence( info["residues"], atomtypes, aInsertGaps=False, aInitialGaps=False, aSplitToRowsOfLength=row_len, aPerModelAndChain=print_stats)
        if sequence_selection == 2: fasta_seq = FastaSequence( info["residues"], atomtypes, aInsertGaps=True,  aInitialGaps=False, aSplitToRowsOfLength=row_len, aPerModelAndChain=print_stats)

        #print fasta_seq

        #Only seq and no stats
        if not print_stats and len( statkeys) == 0:
            models_str = "" #Format models string if file contains models
            if len( info["stats"]["models"]): models_str = "|MODELS:%s" % (rangify( info["stats"]["models"]) if len( info["stats"]["models"]) > 1 else str(info["stats"]["models"][ 0]))

            #print fasta header
            print ">%s|CHAINS:%s%s|ATOMTYPES:%s" % (os.path.basename(inputfile), ",".join(info["stats"]["chains"]), models_str, ",".join( atomtypes))

            if len( info["residues"]) == 0: print ";EMPTY SEQUENCE"
            else: print fasta_seq
            sys.exit( 0) #all done


    #Print file info (stats) and now the included lines from file
    if print_stats or len( statkeys) > 0:

        if type( fasta_seq) == str:
            info["stats"]["sequence"] = zip( info["stats"]["chains"], fasta_seq.split("/")) #Insert sequence to stats
        else:
            info["stats"]["sequence"] = fasta_seq

        print "File '%s' information:" % inputfile

        if len( statkeys):
            for k in statkeys:
                if k not in  info["stats"]: sys.stderr.write( "WARNING: Key '%s' is not a valid file info key.\n" % k)
                else: print "%s : %s" % ( k, info["stats"][ k])
        else: #print all keys
            for k in sorted( info["stats"].iterkeys()):
                print k,
                print ":",
                print info["stats"][ k]

        sys.exit( 0) #all done

    if sequence_selection != 0: sys.exit( 0) #all done

    WriteToFile( outfile_handle, info["atoms"], atomtypes)
    sys.exit( 0) #all done



if __name__ == "__main__":
    main()

    #Example usage and testing:

    #pdb = ReadPDB( "pdb3ikm.ent")
    #pdb = ReadPDB( "pdb2vn9.ent")
    #pdb = ReadPDB( "pdb5afu.ent")
    #pdb = ReadPDB( "pdb3uei.ent")
    #pdb = ReadPDB( "pdb5hy8.ent")
    #pdb = ReadPDB( "pdb3dtp.ent")
    #pdb = ReadPDB( "pdb4gb3.ent")

    #print "missing ratio:", pdb["stats"]["incomplete_residue_ratio"]

    #WriteBiologicalAssembliesToFiles( pdb, aOutputFolder="../asa/", aChains=[], aBiomolecules=[], aCenter=[64.705,  68.537,  61.114])
    #WriteBiologicalAssembliesToFiles( pdb, aOutputFolder="../asa/", aChains=[], aBiomolecules=[], aCenter=[-15.100,  35.483, 237.376])

    #assemblies = ReadBiologicalAssemblies( pdb["remarks"])
    #for ass in assemblies:
    #    new_pdb = {"residues":[],"atoms":[]}
    #    CreateBiologialAssembly( ass, pdb, new_pdb)
    #    WriteToFile( "I:/biomol_%i.pdb" % ass["num"], new_pdb["atoms"], chains=[], types=[])


    #def FastaSequence( aResidueList, aAtomTypes=["A","D","I"], aInsertGaps=False, aInitialGaps=False, aSplitToRowsOfLength=0, aPerModelAndChain=False, aQuiet=False ):
    #print FastaSequence( pdb["residues"], aInsertGaps=True, aInitialGaps=True, aPerModelAndChain=True, aAtomTypes=["A","D"] )

    #a_res = FindResidues( pdb["residues"], aResNum="117", aChain="A")
    #for m_res in pdb["residues"]:
    #    print m_res.IsMatchingResidue( a_res), m_res.PrintCA(),

    #res = FindResidues( pdb["residues"], aResNum="458", aChain="A", aModel=-1, aPositions=6000)
    #count = 0
    #for r in res:
    #  count += 1
    #  print count, ":",
    #  r.PrintCA()

    #reslist = ResiduesInRange( FindResidue( pdb["residues"], aResNum="458", aChain="A", aModel=-1), pdb["residues"], aDistance=10.0 )
    #for r in reslist: r.PrintCA()
    #print SimpleFindDuplicateChains( pdb, aChainsToReport=["A"])