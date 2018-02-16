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

import sqlite3
import hashlib #import md5
from datetime import datetime, date

REP_CONN = None
REP_FILEPATH = "./blast_repository.sqlite"
REP_MAINTENANCE_DONE = False

def SetRepositoryFile( aFilepath):
    global REP_FILEPATH, REP_CONN
    REP_FILEPATH = aFilepath

    if REP_CONN: REP_CONN.close()
    REP_CONN = None

def CloseRepository( aDoMaintenance=True):
    global REP_CONN, REP_MAINTENANCE_DONE

    if REP_CONN != None:
        if not REP_MAINTENANCE_DONE and aDoMaintenance: DeleteOld()
        REP_CONN.close()
        REP_CONN = None

def ConnectToRepository():

    global REP_FILEPATH, REP_CONN

    if REP_CONN != None: return REP_CONN

    if os.path.isfile( REP_FILEPATH) and int( os.stat( REP_FILEPATH).st_size) >0:
        REP_CONN = sqlite3.connect( REP_FILEPATH, detect_types=sqlite3.PARSE_DECLTYPES)
        #conn.text_factory = bytes
        return REP_CONN

    return CreateRepository( REP_FILEPATH)


def CreateRepository( aDB_Filename):

    global REP_CONN

    if os.path.isfile( aDB_Filename) and int( os.stat( aDB_Filename).st_size) >0:
        sys.stdout.write( "INFO: Repository file '%s' already exists.\n")
        return

    folder = os.path.split( aDB_Filename)[ 0]
    if not os.path.exists( folder): os.makedirs( folder)

    REP_CONN = sqlite3.connect( aDB_Filename, detect_types=sqlite3.PARSE_DECLTYPES)
    #conn.text_factory = bytes

    REP_CONN.execute('''CREATE TABLE BLAST_RESULTS
                        (HASH CHAR(37) PRIMARY KEY NOT NULL,
                         COUNT INT NOT NULL,
                         RESULTS TEXT NOT NULL,
                         INSERTED DATE NOT NULL);''')

    return REP_CONN


def SeqHash( aSeq):
    #For queries of length 37 or below, use sequene directly
    if len( aSeq) <= 37: return aSeq
    #For longer seqs use md5 + length as hash
    return hashlib.md5( aSeq.strip()).hexdigest() + ("%05i" % min( 99999, len( aSeq)))

def PrintEntry( aSeq):

    global REP_CONN, REP_FILEPATH

    hashcode = SeqHash( aSeq)

    db = ConnectToRepository()

    cur = db.cursor()
    cur.execute( "SELECT * FROM BLAST_RESULTS WHERE HASH='%s'" % hashcode)
    row = cur.fetchone()

    if row:
        print "Hash:      ", row[ 0]
        print "MaxResults:", row[ 1]
        print "Inserted:  ", row[ 3]
        print "Results:   ",  "None" if (not row[ 2] or not len( row[ 2])) else row[ 2].replace("\n", "\n            ").rstrip()
        #print row
    else:
        print "No entry found"

def FetchFromRepository( aSeq, aCount=-1, aMaxDaysAgo=60 ):

    global REP_CONN, REP_FILEPATH

    hashcode = SeqHash( aSeq)

    db = ConnectToRepository()

    cur = db.cursor()
    cur.execute( "SELECT * FROM BLAST_RESULTS WHERE HASH='%s'" % hashcode)
    row = cur.fetchone()

    if row:

        days_ago = (date.today() - row[ 3]).days #Five days ago == 5
        #print "DAYS AGO:", days_ago
        if days_ago > aMaxDaysAgo: return None

        #Check if result are sufficient
        max_results = row[ 1]
        results = str(row[ 2]).split("\n")
        n_results = (0 if results[ 0] == "" else len( results))

        if n_results == 0: return [] #return empty array if no results (always sufficient)

        if max_results > 0 and max_results == n_results: #If saved output has been limited, and we might be missing results from full set
            if aCount <= 0: return None #We might be missing some
            elif aCount > max_results: return None #We might be missing some

        #If repository has more results that was asked for, need to limit the number of results
        if aCount > 0 and n_results > aCount: return LimitResults( results, aCount)

        #return as is
        return results #return BLAST results


    return None

def DeleteEmpty():

    global REP_MAINTENANCE_DONE

    db = ConnectToRepository()
    cur = db.cursor()

    print "Deleting empty result entries...",

    cur.execute( "DELETE FROM BLAST_RESULTS WHERE LENGTH( RESULTS) == 0")

    #for row in db.execute( "SELECT * FROM BLAST_RESULTS WHERE LENGTH( RESULTS) == 0 LIMIT 50"):
    #    print "hash: ", row[ 0], " max:", row[ 1], " date:", row[ 3]
    #    print "   ", row[ 2].replace("\n", "\n    ")#.split("\n")

    db.commit()
    print "Done."

    #REP_MAINTENANCE_DONE = True

def DeleteOld( aDayLimit=60):

    global REP_MAINTENANCE_DONE

    db = ConnectToRepository()
    cur = db.cursor()

    print "INFO: Removing old entries (>%i days) from BLAST repository..." % aDayLimit,
    cur.execute( "DELETE FROM BLAST_RESULTS WHERE INSERTED <= date('now','-%i day')" % aDayLimit)
    db.commit()
    print "Done."

    REP_MAINTENANCE_DONE = True

def DeleteNew( aDayLimit=1):

    global REP_MAINTENANCE_DONE

    db = ConnectToRepository()
    cur = db.cursor()

    print "INFO: Removing new entries (deposited in the last %i days) from BLAST repository..." % aDayLimit,
    cur.execute( "DELETE FROM BLAST_RESULTS WHERE INSERTED >= date('now','-%i day')" % aDayLimit)
    db.commit()
    print "Done."

    REP_MAINTENANCE_DONE = True

def SelectOld( aDayLimit=60):

    db = ConnectToRepository()

    print "Entries older than %i days:" % aDayLimit
    for row in db.execute( "SELECT * FROM BLAST_RESULTS WHERE INSERTED <= date('now','-%i day')" % aDayLimit):
        print "hash: ", row[ 0], " max:", row[ 1], " date:", row[ 3]
        print "   ", row[ 2].replace("\n", "\n    ")#.split("\n")
        #print ""

def SelectToday():

    db = ConnectToRepository()

    found = False

    print "Entries for today:"
    for row in db.execute( "SELECT * FROM BLAST_RESULTS WHERE INSERTED >= date('now','-1 day')"):
        #print "Type:", type( row)
        print "hash: ", row[ 0], " max:", row[ 1], " date:", row[ 3]
        print "   ", row[ 2].replace("\n", "\n    ")#.split("\n")
        #print ""
        #print "\n".join( row[ 2]) #.split("\n")
        found = True

    if not found: print "No Entries deposited today."

def LimitResults( aResults, aLimit):

    prev_chain = ""
    prev_pdbcode = ""

    results = []
    if len( aResults) == 0: return results

    #Need to check how many unique results
    for blast_entry in aResults:

        cur_chain = "?"
        cur_pdbcode = "????"

        blast_cols = blast_entry.split( "\t")

        #Two possible formats in BLAST results:
        if blast_cols[ 1].find( "|") >= 0:
            #gn|P23919|DTYMK    gi|359545626|pdb|2XX3|A 100.00  212 0   0   1   212 21  232 7e-156   436
            struct_cols = blast_cols[ 1].split( '|')
            cur_pdbcode = struct_cols[ 3].strip()
            cur_chain = struct_cols[ 4]
        else:
            #gn|P30613_Pos2|PKLR     4IP7_A  99.815  541     1       0       34      574     3       543     0.0     1096
            struct_cols = blast_cols[ 1].split( '_')
            cur_pdbcode = struct_cols[ 0].strip() #3IKM
            cur_chain = struct_cols[ 1]

        if cur_pdbcode != prev_pdbcode or cur_chain != prev_chain:
            results.append( blast_entry)
            if len( results) >= aLimit: break
            prev_chain = cur_chain
            prev_pdbcode = cur_pdbcode

    #return "\n".join( results)
    return results

#aMaxCount -1 for all possible
def AddToRepository( aSeq, aBlastResults, aMaxCount=-1, aForceCommit=False ):

    hashcode = SeqHash( aSeq)
    today = date.today()
    aMaxCount = max( 0, aMaxCount)

    db = ConnectToRepository()
    cur = db.cursor()

    cache_result_lines = [("blast_cache" + l[ l.index("\t"):]) for l in aBlastResults] #change first col to "blast_cache"

    #cur.execute( """INSERT INTO BLAST_RESULTS (HASH, COUNT, RESULTS, INSERTED) VALUES( ?,?,?,? )
    #                ON DUPLICATE KEY UPDATE COUNT=VALUES(COUNT), RESULTS=VALUES(RESULTS), INSERTED=VALUES(INSERTED);""", (hashcode, aMaxCount, aBlastResults, today))
    cur.execute( "INSERT OR REPLACE INTO BLAST_RESULTS (HASH, COUNT, RESULTS, INSERTED) VALUES( ?,?,?,? )", (hashcode, aMaxCount, "\n".join( cache_result_lines), today))
    if aForceCommit: db.commit()

    #print "Inserted: '%s'" % aSeq

def debug():

    #EXAMPLE USAGE AND TESTING:

    global REP_FILEPATH
    REP_FILEPATH = "../../Blast_DB/blast_repository.sqlite"

    #AddToRepository( "ABC", ["gn|P23919|DTYMK\tgi|359545626|pdb|2XX3|A\t100.00\t212\t0\t0\t1\t212\t21\t232\t7e-156\t436","gn|P23919|DTYMK\tgi|359545626|pdb|2XX3|B100.00\t212\t0\t0\t1\t212\t21\t232\t7e-156\t436"], -1)
    #AddToRepository( "CBA", [], 3)

    #DeleteEmpty()

    print FetchFromRepository( "TYSSGYVNLSPENKFQNSAL", -1, 4 )
    print FetchFromRepository( "HGILARRPSYRKILKDLSSE", -1, 6 )

    #KQFPRDENPVSHVYLEVVSMHFSKSKKIPITYDNGFLFIHTDKPVYTPDQSVKIRVYSLSDDLKP
    #print "'%s'" % FetchFromRepository( "MVALSLKISIGNVVKTMQFEPSTMVYDACRIIRERIPEAPAGPPSDFGLFLSDDDPKKGIWLEAGKALDYYMLRNGDTMEYRKKQRPLKIRMLDGTVKTIMVDDSKTVTDMLMTICARIGITNHDEYSLVRELMEEKKEEITGTLRKDKTLLRDEKKMEKLKQKLHTDDELNWLDHGRTLREQGVEEHETLLLRRKFFYSDQNVDSRDPVQLNLLYVQARDDILNGSHPVSFDKACEFAGFQCQIQFGPHNEQKHKAGFLDLKDFLPKEYVKQKGERKIFQAHKNCGQMSEIEAKVRYVKLARSLKTYGVSFFLVKEKMKGKNKLVPRLLGITKECVMRVDEKTKEVIQEWNLTNIKRWAASPKSFTLDFGDYQDGYYSVQTTEGEQIAQLIAGYIDIILKKKKSKDHFGLEGDEESTMLEDSVSPKKSTVLQQQYNRVGKVEHGSVALPAIMRSGASGPENFQVGSMPPAQQQITSGQMHRGHMPPLTSAQQALTGTINSSMQAVQAAQATLDDFDTLPPLGQDAASKAWRKNKMDESKHEIHSQVDAITAGTASVVNLTAGDPAETDYTAVGCAVTTISSNLTEMSRGVKLLAALLEDEGGSGRPLLQAAKGLAGAVSELLRSAQPASAEPRQNLLQAAGNVGQASGELLQQIGESDTDPHFQDALMQLAKAVASAAAALVLKAKSVAQRTEDSGLQTQVIAAATQCALSTSQLVACTKVVAPTISSPVCQEQLVEAGRLVAKAVEGCVSASQAATEDGQLLRGVGAAATAVTQALNELLQHVKAHATGAGPAGRYDQATDTILTVTENIFSSMGDAGEMVRQARILAQATSDLVNAIKADAEGESDLENSRKLLSAAKILADATAKMVEAAKGAAAHPDSEEQQQRLREAAEGLRMATNAAAQNAIKKKLVQRLEHAAKQAAASATQTIAAAQHAASTPKASAGPQPLLVQSCKAVAEQIPLLVQGVRGSQAQPDSPSAQLALIAASQSFLQPGGKMVAAAKASVPTIQDQASAMQLSQCAKNLGTALAELRTAAQKAQEACGPLEMDSALSVVQNLEKDLQEVKAAARDGKLKPLPGETMEKCTQDLGNSTKAVSSAIAQLLGEVAQGNENYAGIAARDVAGGLRSLAQAARGVAALTSDPAVQAIVLDTASDVLDKASSLIEEAKKAAGHPGDPESQQRLAQVAKAVTQALNRCVSCLPGQRDVDNALRAVGDASKRLLSDSLPPSTGTFQEAQSRLNEAAAGLNQAATELVQASRGTPQDLARASGRFGQDFSTFLEAGVEMAGQAPSQEDRAQVVSNLKGISMSSSKLLLAAKALSTDPAAPNLKSQLAAAARAVTDSINQLITMCTQQAPGQKECDNALRELETVRELLENPVQPINDMSYFGCLDSVMENSKVLGEAMTGISQNAKNGNLPEFGDAISTASKALCGFTEAAAQAAYLVGVSDPNSQAGQQGLVEPTQFARANQAIQMACQSLGEPGCTQAQVLSAATIVAKHTSALCNSCRLASARTTNPTAKRQFVQSAKEVANSTANLVKTIKALDGAFTEENRAQCRAATAPLLEAVDNLSAFASNPEFSSIPAQISPEGRAAMEPIVISAKTMLESAGGLIQTARALAVNPRDPPSWSVLAGHSRTVSDSIKKLITSMRDKAPGQLECETAIAALNSCLRDLDQASLAAVSQQLAPREGISQEALHTQMLTAVQEISHLIEPLANAARAEASQLGHKVSQMAQYFEPLTLAAVGAASKTLSHPQQMALLDQTKTLAESALQLLYTAKEAGGNPKQAAHTQEALEEAVQMMTEAVEDLTTTLNEAASAAGVVGGMVDSITQAINQLDEGPMGEPEGSFVDYQTTMVRTAKAIAVTVQEMVTKSNTSPEELGPLANQLTSDYGRLASEAKPAAVAAENEEIGSHIKHRVQELGHGCAALVTKAGALQCSPSDAYTKKELIECARRVSEKVSHVLAALQAGNRGTQACITAASAVSGIIADLDTTIMFATAGTLNREGTETFADHREGILKTAKVLVEDTKVLVQNAAGSQEKLAQAAQSSVATITRLADVVKLGAASLGAEDPETQVVLINAVKDVAKALGDLISATKAAAGKVGDDPAVWQLKNSAKVMVTNVTSLLKTVKAVEDEATKGTRALEATTEHIRQELAVFCSPEPPAKTSTPEDFIRMTKGITMATAKAVAAGNSCRQEDVIATANLSRRAIADMLRACKEAAYHPEVAPDVRLRALHYGRECANGYLELLDHVLLTLQKPSPELKQQLTGHSKRVAGSVTELIQAAEAMKGTEWVDPEDPTVIAENELLGAAAAIEAAAKKLEQLKPRAKPKEADESLNFEEQILEAAKSIAAATSALVKAASAAQRELVAQGKVGAIPANALDDGQWSQGLISAARMVAAATNNLCEAANAAVQGHASQEKLISSAKQVAASTAQLLVACKVKADQDSEAMKRLQAAGNAVKRASDNLVKAAQKAAAFEEQENETVVVKEKMVGGIAQIIAAQEEMLRKERELEEARKKLAQIRQQQYKFLPSELRDEH", -1)
    #DeleteOld()

    #SelectOld( 60)
    #SelectToday()
    #DeleteNew()
    #SelectToday()
    CloseRepository( aDoMaintenance=False)


if __name__ == "__main__":

  debug()
  sys.stdout.write("\n")
  sys.stderr.write("\nERROR: module '%s' is not standalone.\n" % __file__)

