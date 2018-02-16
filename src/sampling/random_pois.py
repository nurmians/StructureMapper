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
import time
import random
from random import SystemRandom

from datetime import datetime

from Bio import SeqIO

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append( "%s/../utils" % curdir)
from struct_methods import *
from fasta_formatter import ConvertFastaLocatorFormat, ExtractCustomValuesFromFastaHeader
from fileops import *
import progress_bar as pb


#INPUT / OUTPUT
selection_file = "Human_proteome_Uniprot.fasta"
selections_to_make = 5000
#output_file = "Human_proteome_%i_random_Y.fasta" % selections_to_make
#suitable_aas = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] #All
#suitable_aas = ["A","C","D","E","F","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] #No G
#suitable_aas = ["K"] #K
#suitable_aas = ["K","R"]
suitable_aas = ["S"]
#suitable_aas = ["S","T","Y"] #Phosphorylation
output_file = "Human_proteome_%i_random_%s.fasta" % (selections_to_make,"".join(suitable_aas))

print "Reading file '%s'..." % selection_file

seq_record_index = SeqIO.index( selection_file, "fasta", key_function=GetRecordACC)
n_available_selections = len( seq_record_index)

if n_available_selections == 0:
	print "No sequences available"
	sys.ext( -1)

print "Found %i protein sequences to select from." % n_available_selections


seq_record_indices = {}
i = -1
for k in seq_record_index.keys():
	i += 1
	seq_record_indices[ i] = k

dt = datetime.now()
#dt.microsecond

#print "Generating %i random numbers..." % selections_to_make,
#random.SystemRandom( dt.microsecond)
cryptogen = SystemRandom( dt.microsecond)


#seq_indices = [cryptogen.randrange( n_available_selections) for i in range( selections_to_make)]
#seq_indices = [cryptogen.randrange( 1) for i in range(10)]

print "Selecting %i sequences at random" % selections_to_make

n_pois = 0
n_notfound = 0
n_pois_written = 0
#GetFolder( )
selected_pois = {} #No duplicates
selected_aas = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0}

output_file = ChangeToFolder( output_file, curdir)
print "Writing output to file: '%s'" % output_file
#sys.exit( 0)

pb.InitProgressBar( n_pois_written, selections_to_make, 1.0)

with open( output_file, "w") as of:

	while n_pois_written < selections_to_make:

	#for header, seq in seqs:

		rand_index = cryptogen.randrange( n_available_selections)
		seq_key = seq_record_indices[ rand_index]
		cur_seq = str( seq_record_index[ seq_key].seq).replace( "*", "")

		if len( cur_seq) < 25: continue #Skip short sequences
		seq_len = len( cur_seq)
		header = seq_record_index[ seq_key].description

		poi_positions = []
		if seq_key in selected_pois: poi_positions = selected_pois[ seq_key]
		tries = 0

		while True:
			if tries >= 1500:
				pb.PrintAboveProgressBar( "No suitable poi found for seq: %s" % header )
				n_notfound += 1
				break
			new_pos = cryptogen.randrange( seq_len) #Index of aa
			aa = cur_seq[ new_pos]

			if new_pos not in poi_positions and aa in suitable_aas: #["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]: #All
				poi_positions.append( new_pos)
				n_pois += 1
				pb.PrintAboveProgressBar( "%i: %s: inserting POI %s@%i" % (n_pois, header[:20], cur_seq[ new_pos], new_pos))
				poi_seq = cur_seq[:new_pos+1] + "*" + cur_seq[new_pos+1:] #Insert * _AFTER_ aa
				of.write( ">" + header + "\n")
				of.write( WrapSeq( poi_seq))
				selected_pois[ seq_key] = poi_positions[:]
				n_pois_written += 1
				selected_aas[ aa] += 1
				pb.ProgressBar( n_pois_written, selections_to_make)
				break
			else:
				tries += 1

		if n_notfound > 1000000:
			print "ERROR: Could not find enough suitable amino acids"
			sys.exit( -2)

pb.FinalizeProgressBar( n_pois_written, selections_to_make)

print "Generated %i POIs." % n_pois_written

print "Results:"
for key, val in selected_aas.iteritems():
	percentage = (float( val) / n_pois_written) * 100.0
	sys.stdout.write( "%s\t%s\t%.1f%%\t%.4f\n" % (key, val, percentage, percentage/100))

print "All done."

