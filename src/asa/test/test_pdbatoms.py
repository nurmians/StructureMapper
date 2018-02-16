#!/usr/bin/python
# coding=UTF-8
# -*- coding: UTF-8 -*-

import sys
import os

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append( "%s/../" % curdir)

import pdbatoms


def main():

	test_pass = 0
	test_fail = 0

	result = TestFastaSeq()

	#PrintResult( result)



def TestFastaSeq():


	test_pass = 0
	test_fail = 0

	#def ReadPDB( infile, chains=[], models=[], firstmodelonly=False, aQuiet=False, aConvertPTMs=False, aWarnUnknown=True ):
	pdb = {}
	try:
		pdb = pdbatoms.ReadPDB( "pdb3vuw.ent", chains=["A","E"], aQuiet=True )
	except Exception as ex:
		print "Test failed to run:", str( ex)
		#PrintResult( None) #Test incomplete
		return




	print "TEST: Fasta sequence extraction 1",
	seq = pdbatoms.FastaSequence( filter( lambda x: x.Chain() == "A", pdb["residues"]), aQuiet=True )
	real_seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLR"
	result = (seq == real_seq)
	PrintResult( result)
	if not result:
		print "Missmatch of sequences: '%s'" % seq
		print "                      : '%s'" % real_seq

	print "TEST: Fasta sequence extraction 2",
	seq = pdbatoms.FastaSequence( filter( lambda x: x.Chain() == "E", pdb["residues"]), aQuiet=True )
	real_seq = "KQRCRAPACDHFGNAKCNGYCNECFQFKQMY.."
	result = (seq == real_seq)
	PrintResult( result)
	if not result:
		print "Missmatch of sequences: '%s'" % seq
		print "                      : '%s'" % real_seq

	print "TEST: Fasta sequence extraction 3",
	seq = pdbatoms.FastaSequence( filter( lambda x: x.Chain() == "E", pdb["residues"]), aAtomTypes=["A"], aQuiet=True )
	real_seq = "KQRCRAPACDHFGNAKCNGYCNECFQFKQMY"
	result = (seq == real_seq)
	PrintResult( result)
	if not result:
		print "Missmatch of sequences: '%s'" % seq
		print "                      : '%s'" % real_seq


	try:
		pdb = pdbatoms.ReadPDB( "pdb2uv2.ent", chains=["A"], aQuiet=True )
	except Exception as ex:
		print "Test failed to run:", str( ex)
		#PrintResult( None) #Test incomplete
		return

	print "TEST: Fasta sequence extraction 4",
	seq = pdbatoms.FastaSequence( filter( lambda x: x.Chain() == "A", pdb["residues"]), aAtomTypes=["A"], aInitialGaps=True, aInsertGaps=True, aQuiet=True )
	real_seq = "--------------------YEHVTRDLNPEDFWEIIGELGD--FGKVYKAQNKETSVLAAAKVIDTKSEEELED" \
			   "YMVEIDILASCDHPNIVKLLDAFYYENNLWILIEFCAGGAVDAVMLELERPLTESQIQVVCKQTLDALNYLHDNKIIHRD" \
			   "LKAGNILFTLDGDIKLADFGVSAKNTRTIQRRDSFIGTPYWMAPEVVMCETSKDRPYDYKADVWSLGITLIEMAEIEPPH" \
			   "HELNPMRVLLKIAKSEPPTLAQPSRWSSNFKDFLKKCLEKNVDARWTTSQLLQHPFVTVDSNKPIRELIAEAK"

	result = (seq == real_seq)
	PrintResult( result)
	if not result:
		print "Missmatch of sequences: '%s'" % seq
		print "                      : '%s'" % real_seq


def PrintResult( aResult):
	if aResult == None: print "[SKIP]"
	elif aResult: print "[PASS]"
	else: print "[FAIL]"


if __name__ == "__main__":
    main()

