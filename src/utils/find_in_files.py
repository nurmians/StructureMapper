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
import re
import glob


def FindInFiles( aPaths, aKeywords):


    n_printed = 0

    for path in aPaths:

        filelist = []

        #print "path:", path
        #print "glob:,", glob.glob( path)
        for globfile in glob.glob( path):
            #print "file:", globfile
            filelist.append( globfile)

        #print "list:", filelist
        #print ""
        #print ""


        for read_file in filelist:

            print "Searching from file '%s'..." % read_file
            with open( read_file, "r") as rf:

                linenum = 0
                printed = {}

                for line in rf:

                    linenum += 1
                    words = line.rstrip().split("|")
                    for w in words:
                        w = w.strip()
                        if  len( w) == 4 and w in aKeywords:
                            acc_id = line[:10]
                            if acc_id not in printed:
                            #print "FILE: '%s', LINE %04d: %s" % (read_file, linenum, line.rstrip())
                                print "LINE %7d: %s" % ( linenum, line.rstrip())
                                printed[ acc_id] = True
                                n_printed += 1
                            break

                    #for keyword in aKeywords.keys():
                    #    if keyword.find( line) >= 0:
                    #        print "FILE: '%s', LINE %04d: %s" % (read_file, linenum, line)

    print "Keywords found on %i lines" % n_printed




def main():

    #Testing and example usage:
    files = [   "../set1/input_sequences.fasta.blast", \
                "../set2/input_sequences.fasta.blast", \
                "../set3/input_sequences.fasta.blast", \
                "../set4/input_sequences.fasta.blast" ]
    keywords = [ "1A31","1A81","1AD5","1AOT","1APM","1AYA","1BF5","1BG1","1BKX","1BM2","1BMB","1C8L","1CDK","1CM8","1CMK","1CSY","1CTP","1E9H","1EEN","1EUD","1F34","1F8A","1FBV","1FMK","1FMO","1FOT","1FU5","1FYR","1G1F","1G6G","1GNG","1GPA","1GXC","1GZ2","1H1P","1H1W","1I3Z","1I7W","1IB1","1IRS","1J3H","1J4X","1JBP","1JDY","1JST","1JSU","1JU5","1JYQ","1K3A","1K4S","1KA6","1KDX","1KHX","1KSW","1L3R","1LCJ","1LCK","1LF8","1LKK","1LKL","1MW4","1NEX","1NH3","1O6K","1O6L","1O9U","1OGU","1OHE","1OL5","1OVA","1P22","1PY1","1Q24","1Q61","1Q8W","1QCF","1QMZ","1QPE","1QWO","1R0Z","1R1P","1RDQ","1RQQ","1RR8","1RRJ","1SCU","1SEU","1SHA","1SHC","1SHD","1SMH","1STC","1SVH","1SYK","1SZM","1TCE","1TRN","1TZE","1U54","1U5Q","1U7F","1U7V","1UA2","1UEF","1UHG","1UU3","1UU9","1V54","1W98","1X27","1XH7","1XH9","1XJD","1YDT","1YGR","1YHS","1YJM","1YRK","1YRP","1YVJ","1YVL","1YWT","1Z5M","1Z8D","1ZFP","1ZRZ","2A19","2AFF","2AZM","2B05","2B2T","2BBU","2BFX","2BFY","2BIK","2BR9","2BTP","2BVA","2C1A","2C1J","2C30","2C63","2CBL","2CCI","2CDZ","2CI9","2CIA","2CJM","2CJZ","2CPK","2DVJ","2DXP","2ERK","2ERZ","2F57","2FCI","2FO0","2FP4","2G35","2G57","2G9X","2GCD","2GHQ","2GNF","2GNG","2GNJ","2GQG","2GU8","2H7D","2H8H","2HHJ","2HMH","2I0E","2IVT","2IVV","2IW6","2IW8","2J0L","2J90","2JDO","2JDS","2JDT","2JED","2JFL","2JFM","2JG8","2JGZ","2JQI","2K7L","2KB3","2KFU","2KMD","2L1C","2LAX","2LAY","2LCT","2LDR","2LNW","2LQW","2LXT","2M3B","2MAE","2MFQ","2MRK","2MX4","2N04","2NMB","2NRU","2NRY","2O8Y","2OBJ","2OIB","2OW3","2PIE","2PSG","2PT3","2PTK","2PVF","2Q0N","2QCS","2QG5","2QKW","2QLL","2QO7","2QON","2QUR","2QVS","2R5L","2RLT","2RSY","2RT5","2RVM","2UV2","2UVY","2V7A","2V7D","2V7O","2VAG","2VGO","2VIF","2VO0","2VRX","2VX3","2W1C","2W1I","2W3O","2WTV","2X4Z","2XCH","2XCK","2XIK","2XIX","2XIY","2Y69","2Y94","2YOA","2Z5Z","2ZM1","2ZM3","2ZOQ","3A62","3A77","3A7F","3A8W","3AGL","3AGM","3AJ4","3AL3","3AMA","3ANQ","3BEG","3BHB","3BHT","3BLH","3C4W","3CD3","3CI5","3CKX","3CLY","3COJ","3COM","3CQU","3CY3","3D42","3D5W","3D9P","3DCV","3DDP","3DDQ","3DK6","3DK7","3DND","3DOG","3E5A","3E87","3ERH","3EXH","3EYG","3F2A","3F3Z","3F69","3F7Z","3F88","3FBV","3FHI","3FXX","3FXZ","3GB2","3GQI","3H9F","3H9O","3HA6","3HGK","3HL2","3HRC","3I3W","3IAF","3IDC","3IFQ","3IOP","3IW4","3JBL","3JRW","3K05","3K0H","3K15","3K2L","3KK8","3KK9","3KUL","3KVW","3KXZ","3L4J","3L6F","3L9M","3LJ0","3LPB","3LW1","3LXN","3MA3","3MAZ","3MI9","3MK0","3ML4","3MVJ","3MXC","3N84","3NAY","3NNX","3NUN","3NUS","3NX8","3O2Q","3O7L","3O8I","3OCB","3OJY","3OLL","3OLR","3OP0","3ORX","3OVE","3P9Y","3PFQ","3PFV","3PPZ","3PVB","3PX7","3PXE","3PY3","3Q52","3QAL","3QAM","3QC4","3QD2","3QHR","3QIC","3R5Q","3RTX","3RVG","3RWP","3SAY","3SC1","3SDJ","3SQD","3SRV","3SZM","3TL0","3TL8","3TMP","3TNH","3TNP","3TNQ","3TUY","3TXO","3UAL","3UAT","3UBW","3UEC","3UIG","3UIM","3UOT","3UZD","3WA4","3WDZ","3WE4","3X2U","3ZC6","3ZDI","3ZEW","3ZH8","3ZNI","3ZO2","3ZRK","3ZUV","4A07","4A49","4A4B","4A4C","4A7C","4AE6","4AE9","4AFJ","4AQC","4AZE","4AZF","4B8M","4BCF","4BCK","4BCM","4BCN","4BCO","4BCQ","4BH6","4BJU","4BN1","4BZN","4C0O","4C0S","4C0T","4C2V","4C2W","4C36","4C38","4C4E","4CEG","4CFE","4CFH","4CFU","4CKI","4CKT","4CRS","4CSE","4CT1","4CXA","4D0W","4DAT","4DC2","4DFX","4DFY","4DG3","4DIN","4DIT","4E4L","4E4M","4E6D","4E7W","4EC8","4ELJ","4EOI","4EOJ","4EOK","4EOM","4EON","4EOO","4EOP","4EOQ","4EOS","4EQC","4EUU","4EWH","4EWQ","4EY0","4FIE","4FIF","4FJ3","4FXW","4GFH","4GFU","4GL9","4GPZ","4GV1","4HGE","4HKC","4I3Z","4I41","4IAA","4IAC","4IAN","4II5","4IMI","4IMY","4IUG","4IW0","4IZA","4J1R","4J6S","4JDJ","4JIZ","4JLU","4JMG","4JMH","4JQI","4JS8","4JXT","4K11","4KIK","4KJD","4KS7","4KUJ","4KXF","4L1J","4L1U","4L46","4L7E","4LPA","4M69","4MQ1","4MYG","4NCT","4NM3","4NM5","4NND","4NST","4NTS","4O0V","4O21","4O46","4OR5","4OTD","4PSI","4PSW","4Q9Z","4QBS","4QFG","4QML","4QOZ","4QPM","4QSY","4RA4","4RER","4RGW","4RHG","4RMZ","4ROJ","4U1P","4U8Z","4UN0","4V0G","4V11","4W8E","4WB5","4WB6","4WB7","4WB8","4WIH","4WJN","4WNO","4WRQ","4WZP","4X3F","4X6R","4X8N","4XBR","4XBU","4XHK","4XWH","4XX9","4XZ1","4Y5U","4Y73","4YK9","4YM4","4YMJ","4YOO","4YOS","4Z16","4Z84","4ZHX","4ZJV","4ZLO","4ZPZ","5AMN","5AP1","5AP3","5B72","5BRK","5BRM","5CAW","5CDW","5D39","5DF6","5DH3","5DLT","5DMV","5DMZ","5DN3","5DT3","5DYJ","5E4H","5E50","5ECG","5EFQ","5EG3","5EW9","5EYK","5EZV","5F9E","5FB0","5FC2","5FUR","5FWK","5GLS","5HES","5HLP","5HVK","5IIS","5IL0","5JEJ","5K1E","5K3Y","5K5N","5K9P","5L2W","5LI1","5LI5","5LI9","5LIH"]

    key_dict = dict((k,True) for k in keywords)
    print "Searching for %i keywords." % len( key_dict)

    FindInFiles( files, key_dict)

if __name__ == "__main__":
  main()