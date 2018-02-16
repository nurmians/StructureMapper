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
import operator

class ObjCache:

    def __init__( self, aSize = 10):

        self.SetSize( 10)


    def SetSize( self, aSize):

        self.obj_cache_size = aSize
        self.Reset()


    def Reset( self):

        self.obj_cache = [ None] * self.obj_cache_size
        self.obj_cache_keys = {}
        self.obj_cache_index = 0
        self.key_list = []



    def Add( self, aKey, aObj):

        if aKey not in self.obj_cache_keys:
            self.obj_cache[ self.obj_cache_index] = aObj
            self.obj_cache_keys[ aKey] = self.obj_cache_index #store key index
            self.obj_cache_index += 1
            if self.obj_cache_index >= self.obj_cache_size: self.obj_cache_index = 0

            self.key_list.append( aKey)

            if len( self.key_list) > self.obj_cache_size:
                #remove oldest
                del self.obj_cache_keys[ self.key_list[ 0]]
                del self.key_list[ 0]
            return True #added
        else:
            self.SetRenewed( aKey)


        return False

    def SetRenewed( self, aKey):

        if aKey in self.obj_cache_keys and self.key_list[ -1] != aKey:
            self.key_list.remove( aKey)
            self.key_list.append( aKey)


    def Retrieve( self, aKey):

        if aKey in self.obj_cache_keys:
            self.SetRenewed( aKey)
            return self.obj_cache[ self.obj_cache_keys[ aKey]]

        return None

    def PrintOut( self):

        #Sanity check
        assert( len( self.obj_cache_keys) == len( self.key_list))

        if len( self.obj_cache_keys) == 0:
            print "ObjCache is empty."
            return

        keys = sorted( self.obj_cache_keys.items(), key=operator.itemgetter( 1))
        last_inserted_index = self.obj_cache_index - 1
        if last_inserted_index < 0: last_inserted_index = self.obj_cache_size - 1

        for key, val in keys:

            last_inserted = "" if val != last_inserted_index else " [newest]"
            print "%i: key:'%s' val:'%s'%s" % (val, key, self.obj_cache[ val], last_inserted)


