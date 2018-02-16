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
sys.path.append( "%s/../utils" % curdir)

from struct_methods import *


class MainRecord:

    align_path_and_prefix = ""

    @staticmethod
    def SetAlignPathAndPrefix( aAlignPath ):
        MainRecord.align_path_and_prefix = aAlignPath

    def __init__(self, acc, poi, record_index, blast_entries, exit_on_error=True):
        self.acc = acc #Should include "_Pos"
        self.record_index = record_index
        #self.blast_index  = blast_index
        self.nopos_acc = StripPosTag( acc)
        self.poi = poi

        #if self.nopos_acc not in blast_index:
        #    sys.stderr.write( "ERROR: Main record '%s' has no entry '%s' in the Blast Index.\n" %  (self.acc, self.nopos_acc))
        #    if exit_on_error: Exit( -87)
        if self.acc not in record_index:
            sys.stderr.write( "ERROR: Main record '%s' has no entry in the Sequence record Index.\n" %  (self.acc))
            if exit_on_error: Exit( -88)
        if self.poi < 0: #0 is pos before first amino
            sys.stderr.write( "ERROR: Main record '%s' has no POI defined.\n" %  (self.acc))
            if exit_on_error: Exit( -89)

        #Need to make a copy for filtering later
        #Same blast results can be filtered differently based on seq alignment
        self.blast_entries = blast_entries# blast_index[ self.nopos_acc][:]

    def __eq__( self, another):
        if type(another) == str: return self.acc == another
        return (self.acc == another.acc)
    def __hash__( self):
        return hash( self.acc)

    def Acc( self):
        return self.acc

    def BlastAcc( self):
        return self.nopos_acc

    def SequenceStr( self):
        return str( self.record_index[ self.acc].seq)

    def SequenceRecord( self):
        return self.record_index[ self.acc]

    def Description( self):
        #print "DESC:", self.record_index[ self.acc].description
        return self.record_index[ self.acc].description

    def BlastEntries( self):
        return self.blast_entries
        #return self.blast_index[ self.nopos_acc]

    def Pdbs( self):
        retvals = []
        for b in self.blast_entries:
            try: retvals.append( b["pdb_id"])
            except: pass
        return retvals

    def EntryIdentifier( self, aBlastRank):
         return self.acc + ("_%s" % str( aBlastRank))

    def AlignmentFile( self, aBlastRank):
        if len( MainRecord.align_path_and_prefix) == 0: raise ValueError( "Call SetAlignPathAndPrefix() once before calling this method")
        return MainRecord.align_path_and_prefix + self.EntryIdentifier( aBlastRank) + ".align"
