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
import time

#USAGE:
#import progress_bar as pb
#pb.Init( progress, total)
#pb.PrintAbove("helloes")
#pb.Update( progress, total, aShowSpeed=True, aEstimateCompletion=True)
#pb.Finalize( progress, total)

def PrintAboveProgressBar( aMsg):

    aMsg = aMsg.rstrip()
    if not hasattr( ProgressBar, "progress_len") or not ProgressBar.isatty:
        sys.stderr.write( aMsg + "\n")
        return

    str_filler = int(max(0, ProgressBar.progress_len-len(aMsg)))*" "
    msg = "\r" + aMsg + str_filler + "\n"
    ProgressBar( aForcePrint=True, aWithMsg=msg)

def InitProgressBar( aProgress=0, aTotal=1, aBarLength=30, aUpdateIntervalSeconds=0.5):

    #Using function attributes
    ProgressBar.isatty = sys.stdout.isatty() #Ouputting to a console?
    if ProgressBar.isatty == False: return

    ProgressBar.prev_time = time.time()
    ProgressBar.progress_str = "\rProgress: %%%id/%%i %%s (%%.1f%%%%)%%s%%s " % len( str( aTotal))
    ProgressBar.progress_len = len( ProgressBar.progress_str) + aBarLength + 4
    ProgressBar.cur_prog = max( aProgress, 0)
    ProgressBar.total = max( 1, aTotal)
    ProgressBar.start_progress = ProgressBar.cur_prog
    ProgressBar.run_time = 0.0
    ProgressBar.update_int = 1.0 if aUpdateIntervalSeconds == None else float( aUpdateIntervalSeconds)
    ProgressBar.finalized = False

def FinalizeProgressBar( aProgress, aTotal):
    ProgressBar( aProgress, aTotal, aFinalize=True )

def ProgressBar( aProgress=-1, aTotal=-1, aIntervalSeconds=None, aBarLength=30, aFinalize=False, aForcePrint=False, aShowSpeed=True, aEstimateCompletion=True, aWithMsg=""):

    if not hasattr( ProgressBar, "prev_time"):
        InitProgressBar( aProgress, aTotal, aBarLength, aIntervalSeconds)

    if not ProgressBar.isatty: return #Print only to console

    if aProgress < 0: aProgress = ProgressBar.cur_prog
    else: ProgressBar.cur_prog = aProgress

    if aTotal < 0: aTotal = ProgressBar.total
    else: ProgressBar.total = aTotal

    now = time.time()
    speed_str = ""
    completion_str = ""

    if not aForcePrint and not aFinalize:
        #Check for passed time
        seconds_passed = now - ProgressBar.prev_time
        interval = aIntervalSeconds if aIntervalSeconds != None else ProgressBar.update_int
        if seconds_passed < interval: return #report interval

    if (aShowSpeed or aEstimateCompletion) and ProgressBar.run_time > 0.5 and not aFinalize:
        ProgressBar.run_time += now - ProgressBar.prev_time
        speed_fl = (aProgress - ProgressBar.start_progress) / ProgressBar.run_time
        if aShowSpeed:
            speed_str = " [%.1f/s]" % speed_fl
        if aEstimateCompletion and speed_fl > 0.00001:
            completion_fl = ((aTotal - aProgress) / speed_fl)
            full_m, s = divmod( int( round( completion_fl)), 60)
            h, m = divmod(full_m, 60)
            if h > 99:
                completion_str = " Complete in >4 days."
            else:
                h_str = "" if h == 0 else ("%ih " % h)
                m_str = "" if full_m == 0 else ("%im " % m)
                s_str = "" if full_m >= 10 else ("%02is" % s)
                completion_str = " Complete in %s%s%s" % (h_str,m_str,s_str)

    ProgressBar.prev_time = now

    ratio = min( 1.0, aProgress / float( aTotal))
    p_size = int( round( aBarLength * ratio))
    arrow = ((p_size-1)*"=")+">" if (p_size != 0 and aProgress != aTotal) else (p_size)*"="

    bar = "[%s%s]" % (arrow, ((aBarLength-p_size)*"_"))

    pb_str = ProgressBar.progress_str % ( aProgress, aTotal,  bar, (float( aProgress)/aTotal*100.0), completion_str, speed_str)
    pb_len = len( pb_str)
    if pb_len < ProgressBar.progress_len: pb_str += (" "*(ProgressBar.progress_len-pb_len))
    if pb_len != ProgressBar.progress_len: ProgressBar.progress_len = pb_len

    sys.stderr.write( aWithMsg + pb_str)
    if aFinalize and (not hasattr( ProgressBar, "finalized") or ProgressBar.finalized != True):
        sys.stderr.write( "\n")
        sys.stdout.flush()
        #sys.stdout.write( "\n")
        ProgressBar.finalized = True

#Ease of use
def Init( aProgress=0, aTotal=1, aBarLength=30, aUpdateIntervalSeconds=0.5): InitProgressBar( aProgress=aProgress, aTotal=aTotal, aBarLength=aBarLength, aUpdateIntervalSeconds=aUpdateIntervalSeconds)
def Update( aProgress=-1, aTotal=-1, aShowSpeed=True, aEstimateCompletion=True, aWithMsg=""): ProgressBar( aProgress=aProgress, aTotal=aTotal, aShowSpeed=aShowSpeed, aEstimateCompletion=aEstimateCompletion, aWithMsg=aWithMsg)
def Finalize( aProgress, aTotal): FinalizeProgressBar( aProgress, aTotal)
def PrintAbove( aMsg): PrintAboveProgressBar( aMsg)
