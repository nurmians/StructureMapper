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
import getopt
import math
import time
import ntpath
import re
import signal

from multiprocessing import Process, Pool, Array, Lock, Value, cpu_count

import spherepoints
import pdbatoms
#import fileops

curdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append( "%s/../utils" % curdir)
from struct_methods import SyncPrint, SyncErrPrint
from fileops import *
from delayed_keyboard_interrupt import DelayedKeyboardInterrupt

#Global variables
globvar = { 'distribution_type' : "i",
            'distribution_parameter' : 4,
            'chains' : "",
            'verbose' : False,
            'filemode' : "single",
            'skip_existing' : False,
            'threads' : -1,
            'ASA_level' : "residue",
            'input_file' : "",
            'output_file' : "",
            'include': ["A"],
            #'firstmodelonly':False,
            'probe' : 1.5, #Ångstrom radius
            'mods' : ["C"] #Post-translational modifications: [C]onvert, [A]tomize, [N]one
          }

software_version = "1.0"

modnames = { "C" : "CONVERT_AND_REMOVE_PTMS", "A" : "CONVERT_PTMS", "N" : "None" }

mod_str = "REMARK  68 MODIFICATIONS:"
appl_mod_str = "REMARK  68 APPLIED_MODIFICATIONS:"

#Global variable
printlock = Lock()
number = Value('i', 0)
processed_files = []

def GetFileModifications( aFilePath):

    global modnames, mod_str, appl_mod_str

    mod_list = ["N"]
    appl_mods = 0

    #mod_str = "REMARK  68 MODIFICATIONS:"
    #appl_mod_str = "REMARK  68 APPLIED_MODIFICATIONS:"

    with open( aFilePath, "r") as f:

        for line in f:
            if line.find( mod_str) == 0:

                retval = ModNamesStrToModChars( line[ len( mod_str):])

            elif line.find( appl_mod_str) == 0:
                try:
                    appl_mods = int( line[ len( appl_mod_str):].strip())
                except:
                    sys.stderr.write( "ERROR: could not convert number of applied modifications to integer. ('%s')" % line[ len( appl_mod_str):].strip())

            elif line.find( "ATOM") == 0: break

    return mod_list, appl_mods

def ModNamesStrToModChars( aModNamesStr):

    global modnames

    mod_names_list = aModNamesStr.rstrip().split(",")
    retvals = []

    for mn in mod_names_list:
        mn = mn.strip()
        for key in modnames.keys():
            if modnames[ key] == mn: retvals.append( key)

    if len( retvals) == 0: retvals.append( "N")
    return retvals


def ModsToNames( aMods):

    global modnames
    retval = ""
    for m in aMods:
        if m in modnames:
            if len( retval): retval += ", "
            retval += modnames[ m]
    if len( retval) == 0: retval = modnames["N"]
    return retval

#API
#calc_mode: "R" = residues, "S" = sidechain, "A" = atoms
def ASA( infile, outfile, probe=1.5, sample_multiplier=4, chains = "", threads=-1, report=False, calc_mode="R", verbose=False, include=["A"], firstmodelonly=True, modifications=[] ):

    global globvar, printlock

    assert( calc_mode in ['R','S','A'])

    globvar['distribution_type'] = "i"
    globvar['distribution_parameter'] = sample_multiplier
    globvar['chains'] = chains
    globvar['verbose'] = verbose
    globvar['filemode'] = "single"
    globvar['skip_existing'] = False
    globvar['threads'] = threads
    globvar['ASA_level'] = "residue"
    if calc_mode == "A": globvar['ASA_level'] = "atoms"
    elif calc_mode == "S": globvar['ASA_level'] = "sidechain"
    globvar['input_file'] = infile
    globvar['output_file'] = outfile
    globvar['include'] = include
    globvar['probe'] = probe
    globvar['mods'] = modifications
    #globvar['firstmodelonly'] = firstmodelonly


    return runASA( infile, chains, outfile, threads, report, printlock, firstmodelonly, calc_mode=CalcModeFromStr( globvar['ASA_level']), modifications=globvar['mods']  )


#FUNCTIONS
#def SyncPrint( aStr, aLock=False, aNewLine=True):
#    nl = "\n" if aNewLine else ""
#    if aLock: aLock.acquire()
#    sys.stdout.write( aStr + nl)
#    if aLock: aLock.release()
#
#def SyncErrPrint( aStr, aLock=False, aNewLine=True):
#    nl = "\n" if aNewLine else ""
#    if aLock: aLock.acquire()
#    sys.stderr.write( aStr + nl)
#    if aLock: aLock.release()

def PoolInitializer( *args):
    global number, printlock, globvar
    number, printlock, globvar = args
    signal.signal( signal.SIGINT, signal.SIG_IGN)


def ReportProcessStart( file):
    sys.stdout.write("Processing file %s\n" % (file + ".asa"))


def ReportDone( file):

    global processed_files, printlock, number

    printlock.acquire()

    if len( file) >= 6 and file[0:6] == "Error:":
        sys.stderr.write( file + "\n")
    else:
        processed_files.append( file)

    printlock.release()

    #global printlock, number
    #SyncPrint( "File processed.", printlock )
    #SyncPrint( "Files done: %i %s" % ( str(file)), printlock)



def runASAFiles( files, chains = "", output_folder = "", threads = -1, skip_existing_files=False, firstmodelonly=False, calc_mode="R", modifications=[]):

    global number, printlock, processed_files, globvar

    if skip_existing_files: RemoveExistingASAfilesFromPDBFileList( files, output_folder )

    n_threads = cpu_count()

    if threads < 0:
        n_threads = max( 1, n_threads + threads)
    elif threads > 0:
        #Allow more than core count of threads
        n_threads = threads #min( n_threads, threads)


    #n_threads = cpu_count() if threads < 1 else threads
    n_files = len( files)
    if n_files < 1: return [] #avoid divide by zero
    n_processes_per_file = max( 1, n_threads / n_files)


    #AssertionError: "daemonic processes are not allowed to have children"
    #Small series would have to be processed differently
    if n_processes_per_file > 1:
        if n_files > 1:
            print "To improve performance time, process these %i files without the -e or -d options." % n_files
        else:
            print "To improve performance time, process this file without the -e or -d options."
    n_processes_per_file = 1


    result_files = []

    #Create output files in the same folder as input files if no folder is specified
    if type( output_folder) == str:

        if len( output_folder) == 0:
            output_folder = GetFolder( files[ 0])

        output_folder = InsertFolderSeparator( output_folder)

        result_files = [output_folder + InsertASAToFilePath( GetFilename( file)) for file in files]
    else:
        result_files = output_folder


    if type( chains) == str:
        chains = [chains] * n_files
    if len( chains) != n_files:
        sys.stderr.write( "ASA ERROR: Chains count did not match files count.")
        return [] #RETURN
    if len( result_files) != n_files:
        sys.stderr.write( "ASA ERROR: result files count did not match files count.")
        return [] #RETURN


    args = []
    for f in range( len(files)):
        args.append( "%s\t%s\t%s\t%i\t%i\t%s\t%s" % (files[f], chains[ f], result_files[ f], n_processes_per_file, int( firstmodelonly), calc_mode, "".join(modifications)))

    #Create worker pool
    pool = Pool( processes=n_threads, initializer=PoolInitializer, initargs=(number, printlock, globvar))

    #try:
    for arg in args:
        pool.apply_async( runASASingleArg, (arg,), callback=ReportDone)

    pool.close() # no more tasks

    job_count = len( args)
    progress_counter = 0
    processed_files = []

    try:

        #Trick to keep Ctrl+C interrupt working during processing
        while number.value < job_count:
            time.sleep( 1)
            progress_counter += 1

            if progress_counter == 10:
                progress_counter = 0
                SyncPrint( "ASA PROGRESS: Files processed: %i/%i  (%2.1f%%)" % (number.value, job_count, (float(number.value) / job_count * 100.0)), printlock )
                #SyncPrint( "NUM: %i" % number.value, printlock)

        SyncPrint( "ASA PROGRESS: Files processed: %i/%i  (%2.1f%%)" % (number.value, job_count, (float(number.value) / job_count * 100.0)), printlock )
        print "\nASA INFO: Work has finished."

    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        sys.stderr.write( "\n\nProcess interrupted by user.\n\n")
        return processed_files
    else:
        #pool.close() # no more tasks
        pool.join()  # wrap up current tasks

    return processed_files


def runASAFolder( folder, chains = "", output_folder = "", threads = -1, skip_existing_files=False, calc_mode="R", modifications=[]):


    ok, file_list = ListOfFolderFiles( folder)
    n_total_files = len( file_list)

    if skip_existing_files: RemoveExistingASAfilesFromPDBFileList( file_list, output_folder )

    if ok != 0:
        sys.stderr.write( "Failed to find files to process. Exiting.\n")
        sys.exit( 5)
    if len( file_list) < 1:
        if n_total_files > 0:
            sys.stderr.write( "All files in the directory already calculated.\n")
        else:
            sys.stderr.write( "No PDB files found in folder '%s'. Exiting.\n" % folder)

        sys.exit( 0)

    if len( output_folder) and not os.path.isdir( output_folder):
        sys.stderr.write( "Output directory '%s' is not a valid directory. Exiting." % output_folder)
        sys.exit( 6)


    sys.stdout.write("Found %i files. Starting processing...\n" % len( file_list))


    start = time.time()

    result_files = runASAFiles( file_list, chains, output_folder, threads, skip_existing_files, calc_mode=calc_mode, modifications=modifications)

    end = time.time()


    if len( result_files) == 0:
        sys.stdout.write( "No files were processed.\n")
    else:
        sys.stdout.write( "\n\nCalculations done for %i files in %.2f seconds\n\n" % (len( result_files), float( end - start)))




def runASASeries( input_filename, chains = "", output_folder = "", threads = -1, skip_existing_files=False, calc_mode="R", modifications=[] ):

    if not os.path.isfile( input_filename):
        sys.stderr.write("Warning: File '%s' does not exist.\n" % input_filename)

    output_folder = InsertFolderSeparator( output_folder)

    #Find series files
    basename, file_list = ListOfSeriesFiles( input_filename)

    if len( basename) == 0: basename = "VMD_heatmap"
    heatmap_file = output_folder + basename + ".hm.txt"

    n_total_files = len( file_list)

    if n_total_files > 0:
        sys.stdout.write("%i files found in series.\n" % len( file_list))

    skipped_asa_files = []

    if skip_existing_files:
        skipped_asa_files = RemoveExistingASAfilesFromPDBFileList( file_list, output_folder)

    #print files

    if len( file_list) < 1:

        if n_total_files > 0:
            sys.stderr.write("All files in the series already calculated.\n")

            #Create heatmap from skipped files only
            if len( skipped_asa_files) > 0: CreateHeatmap( skipped_asa_files, heatmap_file, basename + " ASA" )

            return
        else:
            sys.stderr.write("Error finding file series.\n")
            return

    if len( output_folder) and not os.path.isdir( output_folder):
        sys.stderr.write( "Output directory '%s' is not a valid directory. Exiting." % output_folder)
        return



    sys.stdout.write( "Starting processing...\n")

    #Start timer
    start = time.time()

    #Run calculations
    result_files = runASAFiles( file_list, chains, output_folder, threads, calc_mode=calc_mode, modifications=modifications)

    end = time.time()


    if len( result_files) == 0:
        sys.stdout.write( "No files were processed.\n")
    else:
        sys.stdout.write( "\n\nCalculations done for %i files in %.2f seconds\n\n" % (len( result_files), float( end - start)))

    #Include skipped files (previously calculated) in heatmap calculation
    result_files.extend( skipped_asa_files)

    #Create heatmap
    CreateHeatmap( result_files, heatmap_file, basename + " ASA" )




def CreateHeatmap( file_list, output_filename="", title="ASA" ):

    global globvar

    #Gather and create results file from .asa. files

    #Compatible output with VMD heatmapper plugin
    #Open output
    if len( output_filename) == 0: output_filename = "VMD_heatmap.hm.txt"

    try:
        outfile = open( output_filename, 'w')
        outfile_ok = True
    except IOError as e:
        sys.stderr.write( "Output file: %s Error: %s\n" % (output_filename, e.strerror))
        sys.stderr.write( "Cannot write VMD heatmap file.\n" % (output_filename, e.strerror))
        return #RETURN


    results = {}
    ordered_keys = []
    n_files = 0

    file_list = SortFilesByNumber( file_list )

    #Read results
    for result_file in file_list:

        #Try opening and reading RESULT file
        try:
            with open( result_file, 'r') as f:
                #Read PDB file
                if globvar["verbose"]: print "Reading file: %s" % result_file

                #{atoms, residues, stats}
                pdb = pdbatoms.ReadPDB( f )

                if len( pdb["atoms"]) < 1:
                    sys.stderr.write("Warning: file %s did not contain any atom records.\n" % result_file )
                    continue

                #Process by residue
                if globvar["ASA_level"] == "residue" or globvar["ASA_level"] == "sidechain":

                    for residue in pdb["residues"]:

                        #ResNum() may contain an insertion code (char) with the actual residue number
                        resnum = residue.ResNum()
                        if resnum not in results.keys(): results[ resnum] = []
                        results[ resnum].append( residue.atoms[ 0].b)

                        #Keep order
                        ordered_keys.append( resnum)

                #Process by atom
                else:

                    for atom in pdb["atoms"]:

                        str_index = "%i" % atom.index
                        if str_index not in results.keys(): results[ str_index] = []
                        results[ str_index].append( atom.b)

                        #Keep order
                        ordered_keys.append( str_index)

                n_files += 1

        except IOError as e:
            sys.stderr.write( "ASA file: '%s' Error: %s\n" % (result_file, e.strerror))


    if n_files < 1:
        sys.stderr.write( "No files found for generating a heatmap.\n")
        return #RETURN

    sys.stdout.write( "Results for heatmap gathered from %i files.\n" % n_files)


    #Write file
    if outfile_ok:

        outfile.write( '-min "0.00"\n' )
        outfile.write( '-max "100.0"\n' )
        outfile.write( '-xlabel "Series"\n' )
        outfile.write( '-xorigin "0"\n' )
        outfile.write( '-xstep "1"\n' )
        ylabel = "Resids"
        if globvar["ASA_level"] == "atom": ylabel = "Atoms"
        elif globvar["ASA_level"] == "sidechain": ylabel = "Resids (sidechains)"
        outfile.write( '-ylabel "%s"\n' % ylabel)
        outfile.write( '-title "%s"\n' % title  )
        outfile.write( '-numbering "%s"\n' % ('atom' if globvar["ASA_level"] == "atom" else 'resid')    )

        for key in ordered_keys:
            outfile.write( "%s: " % key)
            outfile.write( ';'.join(map(str, results[ key])))
            outfile.write( "\n")


    sys.stdout.write( "VMD heatmap written to '%s'.\n" % output_filename)

#Work around for passing arguments to multiprocessing pool workers
def runASASingleArg( arg_str):

    global printlock, number, globvar

    #SyncPrint( "Num1: %i " % number.value, printlock )
    #number.value += 1
    #print "Num2: %2 " % number.value

    args = arg_str.split("\t")
    #print "ARGS:" + arg_str
    try:
        retval = runASA( args[ 0], args[ 1], args[ 2], int(args[ 3]), False, printlock, bool(args[ 4]), args[ 5], list( args[ 6]))
    except Exception as ex:

        printlock.acquire()
        template = "An exception of type {0} occured. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        sys.stderr.write( message)
        number.value += 1
        printlock.release()
        return ("Error: failed to process file: '%s'" % args[ 0])

    printlock.acquire()
    number.value += 1
    printlock.release()

    return retval



def IsOfIncludedType( atom ):

    global globvar
    acc = globvar["include"]

    if atom.isaminoacid and "A" in acc: return True
    if atom.isdna       and "D" in acc: return True
    if atom.isunknown   and "U" in acc: return True
    if atom.ision       and "I" in acc: return True
    if atom.isother     and "O" in acc: return True

    return False

#chains format: chains:models
#either can be omitted
#model 0 means first model only
#models can contain ranges
#separate multiple by commas
def ExtractChainsAndModels( aStr):

    #Process chains and models to list format
    chain_arg = aStr.split(":")
    firstmodelonly = 0

    chain_list = chain_arg[ 0].split(",")
    if len( chain_list) == 1 and len( chain_list[ 0]) > 1: chain_list = list( chain_arg[ 0]) #No commas used to separate chains
    elif len( chain_list) == 1 and chain_list[ 0].strip() == "*": chain_list = [] #All

    model_list = []
    if len( chain_arg) > 1:
        if chain_arg[ 1] == "*": model_list = [] #All
        else:
            model_list = pdbatoms.derangify( chain_arg[ 1])
            if len( model_list) == 1 and model_list[ 0] == 0: #Model 0 means first model only
                firstmodelonly = 1
                model_list = []

    #filter = Remove empty strings
    return (filter( None, chain_list), filter( None, model_list), firstmodelonly)

#RUN ASA  #chains format: chains:models
def runASA( input_filename, chains = "", output_filename = "", threads = -1, report_progress=True, printlock = False, firstmodelonly=False, calc_mode="R", modifications=[]):

    chain_list, model_list, arg_firstmodelonly = ExtractChainsAndModels( chains)
    if arg_firstmodelonly: firstmodelonly = arg_firstmodelonly

    if not os.path.isfile( input_filename):
        sys.stderr.write( "ASA ERROR: Input file '%s' does not exist.\n" % input_filename)
        return ""

    pdb = []
    locked = False
    n_ptms = 0 #Post-translational modifications

    #Try opening and reading INPUT file
    try:

        with open( input_filename, 'r') as infile:
            #Read PDB file

            #{atoms, residues, stats}

            #Acquire lock to print all messages from processing the pdb file without interruptions
            if printlock:
                printlock.acquire()
                locked = True

            sys.stdout.write( "ASA INFO: Processing file: %s\n" % input_filename )

            pdb = pdbatoms.ReadPDB( infile, chain_list, model_list, firstmodelonly=firstmodelonly )

            if globvar["verbose"]:
                sys.stdout.write( "Read %i atom records.\n" % ( len( pdb[ "atoms"])))
                printDictionary( pdb[ "stats"])

            if locked:
                printlock.release()
                locked = False

            if "C" in modifications or "A" in modifications:
                remove_ptms = "C" in modifications
                n_ptms = pdbatoms.ConvertPTMAAstoRegularAtoms( pdb, remove_ptms)
                if n_ptms > 0:
                    sys.stderr.write( "ASA INFO: Converted %i post-translationally modified residues.\n" % n_ptms)


    except IOError as e:

        if locked: printlock.release()

        SyncErrPrint( "File: '%s' Error: %s" % (input_filename, e.strerror), printlock)
        SyncErrPrint( "ASA could not be calculated.", printlock)
        return ""


    #check that OUTPUT file is writable
    #If no file has been specified, use stdout
    outfile = None

    if len( output_filename):
        try:
            outfile = open( output_filename, 'w+')
            #SyncErrPrint( "Writing ASA output to file '%s'" % output_filename, printlock)
        except IOError as e:
            SyncErrPrint( "ASA ERROR: Could not write file '%s', error: %s" % (output_filename, e.strerror), printlock)
            return ""
    else:
        outfile = sys.stdout


    #Find atoms for calculations
    atom_list = []
    for a in pdb["atoms"]:
        if IsOfIncludedType( a):
            atom_list.append( a)


    n_atoms = len( atom_list)

    if n_atoms < 1:
        SyncErrPrint(  "ASA WARNING: No ATOM or HETATM records to process found in file: %s." % (input_filename), printlock)
        #outfile.write( "REMARK   TOTAL ASA: 0.0\n")
        WriteRemarks( outfile, 0.0, aFirstModelOnly=firstmodelonly )
        outfile.write( "END\n")
        outfile.close()
        return "" #RETURN


    probe_r = float( globvar["probe"])
    probe_d = 2 * probe_r


    #Fill dict with sample points on a spherical surface (either icosaedric or golder ratio)
    atom_spheres = CreateAtomSpheres( probe_r, pdb["stats"][ "atom_types"])
    n_spherepoints = len( atom_spheres[ atom_spheres.keys()[ 0]])


    #Sort atoms into "boxes" in 3D based on probe and max atom size
    max_atom_radius = pdb["stats"]["largest_radius"]
    n_boxes, limits, side_len, boxed_atoms = SortToBoxes( atom_list, max_atom_radius * 2 + probe_r * 2)


    #Use either the number of cores or specified amount of threads
    n_threads = cpu_count() if threads < 1 else threads

    #Make sure there are enough calculations to do to make
    #starting new processes worthwhile
    reduced = max( min( n_threads, n_atoms*n_spherepoints / 100000), 1)
    if reduced < n_threads and report_progress:
        print "Reducing thread count from %i to %i based on the number of required calculations." % (n_threads, reduced)
        n_threads = reduced

    #Limit thread count between 1-500 threads (processes)
    n_threads = min( 500, max( n_threads, 1))


    processes = []

    atoms_per_thread = n_atoms / n_threads + 1

    #Shared memory arrays for returning results
    #'d' for double
    arr_SA    = Array( 'd', [0.0]*n_atoms, lock=False)     #Atom's contribution to residue's Surface Area
    arr_ASA   = Array( 'd', [0.0]*n_atoms, lock=False)     #Atom's Accessible Surface Area
    arr_total = Array( 'd', [0.0]*n_threads, lock=False)

    total_ASA = 0.0
    start = time.time()


    if n_threads == 1:

        #Calculate only in this process

        if report_progress and globvar["verbose"]: SyncPrint( "Using a single process.", printlock)
        if report_progress: SyncPrint( "\nCalculating ASA...", printlock)

        CalcASA( 0, arr_SA, arr_ASA, arr_total, probe_r, atom_spheres, boxed_atoms, atom_list, pdb["residues"], 0, n_atoms, report_progress)

    else:

        #Calculate using multiple Processes
        if report_progress and globvar["verbose"]: SyncPrint(  "Using %i processes." % n_threads, printlock)

        if report_progress: SyncPrint( "\nCalculating ASA...", printlock)

        #Create processes
        for t in range( 1, n_threads):

            r_start = t*atoms_per_thread
            r_end = min( n_atoms, t*atoms_per_thread+atoms_per_thread)

            new_thread = Process( target=CalcASA, args=( t, arr_SA, arr_ASA, arr_total, probe_r, atom_spheres, boxed_atoms,
                                                         atom_list, pdb["residues"], r_start, r_end, report_progress))
            processes.append( new_thread)

        for p in processes:
            p.start()

        #use this process as well
        CalcASA( 0, arr_SA, arr_ASA, arr_total, probe_r, atom_spheres, boxed_atoms, atom_list, pdb["residues"], 0, atoms_per_thread, report_progress)

        for p in processes:
            p.join()
            #print "Joined!"


    end = time.time()

    #Calculate total ASA from each threads' individual results
    #i = 0
    #for result in arr_total:
    #    print "Results %i: %f" % (i, result)
    #    i+=1
    #    total_ASA += result

    total_ASA = sum( arr_total)

    if report_progress:

        print u"Done.\n\nTotal ASA: %8.3f Å²\n\n" % total_ASA

        processes = "" if n_threads == 1 else ("[%ip]" % n_threads)
        print u"ASA calculated for %i atoms (%i sample points each) in %.2fs, with probe radius %.2fÅ. %s\n" % (
                                                                                              n_atoms,
                                                                                              n_spherepoints,
                                                                                              float(end - start),
                                                                                              probe_r, processes),
    #Reset b-cols for ALL atoms
    for a in pdb["atoms"]:
        a.b = 0.0

    #Transfer results from shared memory to atom data structure
    for i in xrange( 0, n_atoms):
        atom = atom_list[ i]
        atom.asa = arr_ASA[ i]
        atom.sa = arr_SA[ i]
        atom.b = atom.asa

    #Process results by residue
    if calc_mode == "S":   ProcessBySidechain( pdb[ "residues"])
    elif calc_mode == "R": ProcessByResidue(   pdb[ "residues"])

    #Write results
    #Write only AMINO ACIDS and DNA
    with DelayedKeyboardInterrupt():

        WriteRemarks( outfile, total_ASA, n_ptms, aFirstModelOnly=firstmodelonly )
        #outfile.write( "REMARK  68 TOTAL ASA: %.3f\n" % total_ASA)
        pdbatoms.WriteToFile( outfile, pdb["atoms"], globvar["include"])
        outfile.close()

    return output_filename


#aCalcMode: R = Per residue, S = Per residue sidechain, A = per atom
#WARNING: Uses the aResidues atom data for the calculations and restores b.col values to original values after calculations are done
def CalcASAForResidues( aResidues, aPdb, aChains=[], aModels=[], aProbe_r=1.5, aIncluded_atom_types=["A","D"], aCalcMode="R", aQuiet=False ):

    if aCalcMode == "S" and aPdb["stats"]["backbone_only_ratio"] > 0.9:
        if not aQuiet: sys.stderr.write( "ASA ERROR: Structure does not contain sidechain information for ASA calculations.\n")
        return [-1.0] * len(aResidues)
    elif aCalcMode == "S" and aPdb["stats"]["incomplete_residue_ratio"] > 0.5:
        if not aQuiet: sys.stderr.write( "ASA WARNING: Structure has incomplete amino acid information for ASA sidechain calculations.\n")

    #Find atoms for calculations
    full_atom_list = []
    for a in aPdb["atoms"]:
        if len( aChains) and a.chain    not in aChains: continue
        if len( aModels) and a.modelnum not in aModels: continue
        if len( aIncluded_atom_types) and not a.IsOfIncludedType( aIncluded_atom_types): continue
        full_atom_list.append( a)

    #print "Full atom list", len( full_atom_list)

    #Fill dict with sample points on a spherical surface (either icosaedric or golden ratio)
    atom_spheres = CreateAtomSpheres( aProbe_r, aPdb["stats"][ "atom_types"])
    n_spherepoints = len( atom_spheres[ atom_spheres.keys()[ 0]])

    max_atom_radius = aPdb["stats"]["largest_radius"]
    n_boxes, limits, side_len, boxed_atoms = SortToBoxes( full_atom_list, max_atom_radius * 2 + aProbe_r * 2)

    atom_list = []
    chains = set()
    for res in aResidues:
        atom_list += res.atoms
        chains.add( res.Chain())
    n_atoms = len( atom_list)

    if len( aChains):
        for c in chains:
            if c not in aChains:
                if not aQuiet: sys.stderr.write( "ASA WARNING: Chain of queried residue ('%s') not included in calculations.\n" % c)
    #print "atom list", n_atoms

    arr_SA    = Array( 'd', [0.0]*n_atoms, lock=False)     #Atom's contribution to residue's Surface Area
    arr_ASA   = Array( 'd', [0.0]*n_atoms, lock=False)     #Atom's Accessible Surface Area
    arr_total = Array( 'd', [0.0]*1, lock=False)

    CalcASA( 0, arr_SA, arr_ASA, arr_total, aProbe_r, atom_spheres, boxed_atoms, atom_list, aPdb["residues"], 0, n_atoms, report_progress=False )

    #for res in aResidues: res.ResetBCol()

    #Transfer results from shared memory to atom data structure
    for i in xrange( 0, n_atoms):
        atom = atom_list[ i]
        atom.asa = arr_ASA[ i]
        atom.sa = arr_SA[ i]
        atom.b = atom.asa

    ret_array = []


    #Process results by residue
    if aCalcMode == "S":
        ProcessBySidechain( aResidues)
        for res in aResidues:
            ret_array.append( res.BColAvg())
            res.RestoreOriginalAtomBColValues()
    elif aCalcMode == "R":
        ProcessByResidue( aResidues)
        for res in aResidues:
            ret_array.append( res.BColAvg())
            res.RestoreOriginalAtomBColValues()
    else:
        #Results by Atom
        for res in aResidues:
            atom_arr = []
            for atom in res.atoms:
                atom_arr.append( atom.b)
                #print atom.atomname, atom.resname, atom.b
                a.RestoreOriginalBcolValue()
            ret_array.append( atom_arr)

    return ret_array



def CalcASA( t_index, arr_SA, arr_ASA, arr_total, probe_r, atom_spheres, boxed_atoms, atoms, residues, r_start, r_end, report_progress=True):

    probe_d = 2*probe_r
    total_ASA = 0.0
    #total_SA = 0.0

    n_spherepoints = len( atom_spheres[ atom_spheres.keys()[ 0]])
    #print "Spherepoints for thread %i\n" % n_spherepoints

    n_boxes = [ len(boxed_atoms), len(boxed_atoms[ 0]), len(boxed_atoms[ 0][ 0])]

    reported = None
    atoms_processed = 0

    #First thread (0) is a reporting thread
    if t_index == 0 and report_progress:
        reported = time.time()
        print "Progress: 0.1%",

    #For each atom in the specified range
    #for atom in atoms[ r_start:r_end]:
    for atom_index in xrange(r_start, r_end):

        atom = atoms[ atom_index]

        #if only_sidechains and atom.IsBackBone():
        #    #Skip backbone atoms
        #    #arr_ASA[ atom_index] = 0.0 #Array already initialized to zeros, no need to set
        #    #arr_SA[ atom_index] = 0.0  #Array already initialized to zeros, no need to set
        #    continue

        atom_accessible_pts = 0
        atom_collision_distance = probe_d + atom.radius

        #atom_pts are probe center points to be tested
        #atom_pts = []
        #if use_numpy: atom_pts = (atom_spheres[ atom.VDW_type] * numpy.array( atom.coord)).tolist()
        #else: atom_pts = TransposeSphere( atom_spheres[ atom.VDW_type], atom.coord )

        atom_pts = TransposeSphere( atom_spheres[ atom.VDW_type], atom.coord )

        #Check for collisisons within the residue
        #This slows down calculation time a little bit, but has the benefit of separating the
        #the surface area that collides within the residue from other collisions
        #sample points that collide withing the residue are removed from the array atom_pts
        #print "Atom: ", atom.ShortLineFormat() #DEBUG
        (atom_surface_pts, checks_made) = CheckForCollisions( atom_pts, residues[ atom.resindex].OtherAtoms( atom), probe_r, deleteColliding=True)

        proxmity_atoms = []
        #skipped_atoms = 0
        #atoms_in_adj_boxes = 0

        #print "ATOM:", atom.Line()
        adj_boxes = AdjacentBoxes( atom.box, n_boxes, True)

        #Check the same box where the current atom is in first
        for adj_box in adj_boxes:
            for box_atom in boxed_atoms[ adj_box[ 0]][ adj_box[ 1]][ adj_box[ 2]]:

                #Atoms of the same residue already checked
                if atom.SameResidue( box_atom): continue

                #Atoms from different models do not collide
                if not atom.SameModel( box_atom): continue

                #atoms_in_adj_boxes += 1

                atom_distance = Distance( atom.coord, box_atom.coord)
                if atom_distance > atom_collision_distance + box_atom.radius:
                    #skipped_atoms += 1
                    continue

                proxmity_atoms.append( box_atom)


        #print "atoms in adj: %i, skipped: %i" % ( atoms_in_adj_boxes, skipped_atoms )
        (atom_accessible_pts, checks_made) = CheckForCollisions( atom_pts, proxmity_atoms, probe_r)


        #Calculations for atom done, set b col
        #Sphere surface area: 4*pi*r^2
        #ASA surface is atom's VDW radius + probe_r spherical surface
        sphere_area = 4*math.pi*(atom.radius+probe_r)*(atom.radius+probe_r)
        ASA = sphere_area * (float(atom_accessible_pts) / n_spherepoints) #Atom's Accessible Surface Area
        SA =  sphere_area * (float(atom_surface_pts) / n_spherepoints)    #Atom's contribution to residue's Surface Area

        #Percentage of how much of the area that does not collide with the residue itself is solvent accessible
        #PASA = 0.00 if atom_surface_pts == 0 else ((float(atom_accessible_pts) / atom_surface_pts) * 100.0)

        #Report results through shared memory arrays
        arr_ASA[ atom_index] = ASA
        arr_SA[ atom_index] = SA

        total_ASA += ASA

        #Report progress at regular time intervals
        if t_index == 0 and report_progress:
            atoms_processed += 1
            now = time.time()
            if float( now - reported) > 4.0:
                reported = now
                print "\rProgress: %.1f%%" % (float( atoms_processed) / ( r_end - r_start) * 100.0),


    #Done!
    if t_index == 0 and report_progress: print "\rProgress: %.1f%%\n" % 100
    arr_total[ t_index] = total_ASA


def WriteRemarks( outfile, total_asa=0.0, aAppliedModifications=0, aFirstModelOnly=0):

    global globvar, modnames, mod_str, appl_mod_str, software_version

    chain_list, model_list, arg_firstmodelonly = ExtractChainsAndModels( globvar['chains'])
    if arg_firstmodelonly: aFirstModelOnly = arg_firstmodelonly

    chain_str = "All"
    model_str = "First" if aFirstModelOnly else "All"

    if len( chain_list) > 0: chain_str = ";".join( chain_list)
    if len( model_list) > 0: model_str = pdbatoms.rangify( model_list, aSeparator=";") #Do not use commas

    outfile.write( "REMARK  68 INCLUDED_CHAINS: %s, INCLUDED_MODELS: %s, INCLUDED_ATOMS:%s\n" % (chain_str, model_str, "".join(globvar['include'])))
    outfile.write( "REMARK  68 ASA_TYPE: %s\n" % globvar['ASA_level'])
    outfile.write( "REMARK  68 ASA_VERSION: %s, PDBATOMS_VERSION: %s\n" % (software_version, pdbatoms.GetVersion()))
    outfile.write( "%s %s\n" % (mod_str, ModsToNames( globvar['mods'])))

    if aAppliedModifications > 0:
        outfile.write( "%s %i\n" % (appl_mod_str, aAppliedModifications))

    sample_type = "Unknown"
    sample_points = -1

    if globvar['distribution_type'] == "i":
        sample_type = "Icosahedric"
        sample_points = spherepoints.NumberOfIcosahedricPoints( globvar['distribution_parameter'])
    else:
        sample_type = "Vogel"
        sample_points = globvar['distribution_parameter']

    outfile.write( "REMARK  68 SAMPLE_POINTS: %i, SAMPLE_TYPE: %s\n" % (sample_points, sample_type))

    outfile.write( "REMARK  68 PROBE_RADIUS: %.3f\n" % globvar['probe'])
    outfile.write( "REMARK  68 TOTAL_ASA: %.3f\n" % total_asa)
    outfile.write( "REMARK  68 ASA_UTA: Anssi Nurminen; Vesa P. Hytonen; University of Tampere; 2018\n")



def ProcessBySidechain( residues):

    return ProcessByResidue( residues, sidechains_only=True)


def ProcessByResidue( residues, sidechains_only=False):

    #For each residue
    for res in residues:

        #Calculate only for included types
        if not IsOfIncludedType( res.atoms[ 0]):
            res.ResetBCol()
            continue

        #average over residue
        #accessible = 0.0
        #total = 0.0
        #over_20 = False
        #maximum = 0.0

        residue_asa = 0.0 #accessible surface area
        residue_sa = 0.0  #surface area
        n_included_atoms = 0

        for atom in res.atoms:

            if sidechains_only and atom.IsBackBone(): continue #Skip backbone atoms if sidechains_only True
            residue_asa += atom.asa
            residue_sa += atom.sa
            n_included_atoms += 1

        percentage = 0.0

        if residue_sa > 0.00001:
            percentage = (residue_asa / residue_sa) * 100.0
        elif n_included_atoms == 0:
            percentage = -1.0

        #Write b-col for residue
        for atom in res.atoms:
            atom.b = percentage




def CheckForCollisions( test_points, box_atoms, probe_r, deleteColliding=False):

    accessible_points = 0
    calculations_made = 0
    #p = 0

    #For each spherepoint for atom...
    for p in reversed( range( len( test_points))):
    #while p < len( test_points):
    #for point in test_points:

        collision = False
        test_point = test_points[ p]

        for box_a in range( len( box_atoms)):

            #measure distance
            box_atom = box_atoms[ box_a]
            probe_collision_distance = box_atom.radius + probe_r
            point_distance = Distance( test_point, box_atom.coord)
            calculations_made += 1

            if point_distance < probe_collision_distance:
                collision = True
                box_atoms.insert( 0, box_atoms.pop( box_a)) #The next point is more likely to collide with same atom, move it to top
                break

        if not collision: accessible_points += 1
        elif deleteColliding: del test_points[ p]

    #print "Accessible points: %i" % accessible_points
    return (accessible_points, calculations_made)


def CreateAtomSpheres( probe_r, atom_types):

    global globvar

    #Fill array of spherical points (either icosaedric or golder ratio)
    sphere = [[0.0,0.0,0.0]]
    atom_spheres = {}

    if globvar["distribution_type"] == "i":
        sphere = spherepoints.IcosahedronPointsOnSphere( globvar["distribution_parameter"])
    else:
        sphere = spherepoints.GoldenSpiralPointsOnSphere( globvar["distribution_parameter"])

    #for atom_type in pdbatoms.VDW_radii.iterkeys():
    for atom_type in atom_types.iterkeys():
        #Include probe radius
        atom_spheres[ atom_type] = SetSphereRadius( sphere, float( pdbatoms.VDW_radii[ atom_type]) + probe_r )

    return atom_spheres


def Distance( coord1, coord2):

    dx = coord2[ 0]-coord1[ 0]
    dy = coord2[ 1]-coord1[ 1]
    dz = coord2[ 2]-coord1[ 2]

    return math.sqrt( dx*dx+dy*dy+dz*dz )


def BoxCoords( atom_coords, side_length, limits ):

    return [int( math.floor( (atom_coords[ 0] - limits[ 0][ 0]) / side_length)),
            int( math.floor( (atom_coords[ 1] - limits[ 1][ 0]) / side_length)),
            int( math.floor( (atom_coords[ 2] - limits[ 2][ 0]) / side_length))]


def AdjacentBoxes( box, box_limits, include_self = True):

    if( box[ 0] >= box_limits[ 0]): sys.stderr.write( "Warning: Box x-coordinate out of bounds\n")
    if( box[ 1] >= box_limits[ 1]): sys.stderr.write( "Warning: Box y-coordinate out of bounds\n")
    if( box[ 2] >= box_limits[ 2]): sys.stderr.write( "Warning: Box z-coordinate out of bounds\n")

    retval = []
    if include_self:
        retval.append( [ box[ 0],box[ 1],box[ 2]] )

    for x in xrange( max( 0, box[ 0]-1), min( box[ 0]+1, box_limits[ 0]-1)+1):
        for y in xrange( max( 0, box[ 1]-1), min( box[ 1]+1, box_limits[ 1]-1)+1):
            for z in xrange( max( 0, box[ 2]-1), min( box[ 2]+1, box_limits[ 2]-1)+1):
                #Do not add box itself
                if x == box[ 0] and y == box[ 1] and z == box[ 2]: continue
                #add adjacent box
                retval.append( [x,y,z])

    return retval


def SortToBoxes( atoms, side_length):

    if side_length < 0.001:
        sys.stderr.write( "Error: Cannot group atoms with side length %6.2f" % side_length )
        return [[[]]]

    if len( atoms) == 0: return []

    limits_x = [ atoms[ 0].coord[ 0], atoms[ 0].coord[ 0]]
    limits_y = [ atoms[ 0].coord[ 1], atoms[ 0].coord[ 1]]
    limits_z = [ atoms[ 0].coord[ 2], atoms[ 0].coord[ 2]]


    for atom in atoms:
        if   limits_x[ 0] > atom.coord[ 0]: limits_x[ 0] = atom.coord[ 0]
        elif limits_x[ 1] < atom.coord[ 0]: limits_x[ 1] = atom.coord[ 0]

        if   limits_y[ 0] > atom.coord[ 1]: limits_y[ 0] = atom.coord[ 1]
        elif limits_y[ 1] < atom.coord[ 1]: limits_y[ 1] = atom.coord[ 1]

        if   limits_z[ 0] > atom.coord[ 2]: limits_z[ 0] = atom.coord[ 2]
        elif limits_z[ 1] < atom.coord[ 2]: limits_z[ 1] = atom.coord[ 2]

    x_length = limits_x[ 1] - limits_x[ 0] + 0.0001
    y_length = limits_y[ 1] - limits_y[ 0] + 0.0001
    z_length = limits_z[ 1] - limits_z[ 0] + 0.0001

    x_boxes = int( math.ceil( x_length / side_length))
    y_boxes = int( math.ceil( y_length / side_length))
    z_boxes = int( math.ceil( z_length / side_length))

    limits = [limits_x, limits_y, limits_z]

    #3D list for holding atoms
    #boxed = x_boxes * [ y_boxes * [ z_boxes * []]]
    boxed = [[[[] for k in xrange(z_boxes)] for j in xrange(y_boxes)] for i in xrange(x_boxes)]
    #print "boxed = %i * [ %i * [ %i * ['a']]]" % (x_boxes, y_boxes, z_boxes)

    for atom in atoms:

        box_coords = BoxCoords( atom.coord, side_length, limits )
        #print "atom in box: %i,%i,%i" % (box_coords[ 0], box_coords[ 1], box_coords[ 2])

        atom.box = [ box_coords[ 0], box_coords[ 1], box_coords[ 2]]

        boxed[ box_coords[ 0]][ box_coords[ 1]][ box_coords[ 2]].append( atom)


    return ([x_boxes, y_boxes, z_boxes], [limits_x, limits_y, limits_z], side_length, boxed)




def SetSphereRadius( sphere, r ):

    retval = []
    for sp in sphere: retval.append( [sp[ 0]*r, sp[ 1]*r, sp[ 2]*r])
    return retval


def TransposeSphere( sphere, coords):

    retval = []
    for sp in sphere: retval.append( [sp[ 0] + coords[ 0], sp[ 1] + coords[ 1], sp[ 2] + coords[ 2]])
    return retval


def usage():

    print u"\n"
    print u"This script calculates accessible surface area (ASA)"
    print u"for amino acids in a PDB structure file using the Lee-Richards method."
    print u"Results of calculation are reported in the b-col of the outputted"
    print u"PDB file. For information on the file format, see: www.rcsb.org/"
    print u""
    print u"Lee B, Richards FM. J Mol Biol. 1971 Feb 14;55(3):379-400."
    print u""
    print u"Authors: Anssi nurminen, Vesa P. Hytönen. 2018, University of Tampere"
    print u""
    print u"Script version: 1.0c"
    print u"Released under the GNU General Public License. Use -l for details."
    print u"Please cite the authors if you find this script useful."
    print u""
    print u""
    print u"Usage: %s [-h|-l] [-v] [-a|-r|-s] [-e|-d] [-k]" % os.path.basename(__file__)
    print u"       [-p RADIUS] [-c CHAINS] [-i MIN|-v POINTS]     inputfile [outputfile]"
    print u""
    print u""
    print u"Arguments:"
    print u" inputfile             PDB-format structure file. (See: www.rcsb.org)"
    print u" outputfile [optional] If an outputfile is not specified the results will be written to"
    print u"                       a file named: 'inputfile.asa.pdb'."
    print u""
    print u"Flags:"
    print u" -a --atom             Calculate ASA for each amino acid atom."
    print u""
    print u" -r --residue          Calculate ASA percentage for each residue. [default]"
    print u"                       The percentage represents the portion of the surface area"
    print u"                       of the residue that can be in contact with the solvent."
    print u""
    print u" -s --sidechain        Calculate ASA for amino acid sidechain atoms only and reports"
    print u"                       the percentage just like with the --residue (-r) option."
    print u""
    print u" -n --include  TYPES   String to specify which atoms to include in ASA calculations."
    print u"                       A=Amino acids, D=DNA/RNA, I=Ions, U=Unknown, O=Other"
    print u"                       Example usage: \"-n ADU\", Default is only \"A\""
    print u""
    print u" -e --series           Input file is the first of a series of files to process."
    print u"                       Other files of the series are found based on a number that needs to"
    print u"                       be present in the filename. Example: \"file001.pdb\"."
    print u"                       In addition, an output file (.hm) is generated that combines results"
    print u"                       in a text file that is compatible with the VMD heatmapper plugin."
    print u"                       The results are stored as .asa.pdb files in the outputfolder."
    print u"                       If no output folder is specified results are stored in the input folder."
    print u""
    print u" -d --directory        Calculate ASA for all .pdb and .ent files in the given directory (inputfile)"
    print u"                       The results are stored as .asa.pdb (or .asa.ent) files in the outputfolder."
    print u"                       If no output folder is specified results are stored in the input folder."
    print u""
    print u" -k --skip             Skip files that already have a .asa.pdb file in the specified output"
    print u"                       folder. Default mode is to recalculate ASA for these files also."
    print u""
    print u" -p --probe RADIUS     Radius of the probe sphere (in Ångströms) that is used for calculating"
    print u"                       surface accessibility. Default: 1.5Å"
    print u""
    print u" -c --chains           Calculate ASA for specified chains and models only. The chains and models that "
    print u"                       are not included will be ignored completely in calculations. Multiple chains and "
    print u"                       models can be specified simultaneously when separated with commas."
    print u"                       Format: \"-c chains:models\". Use \"-c :0\" to process only the first model in file."
    print u"                       Example: \"-c A,C,D:1,2,5-10\". Default: All chains, all models"
    print u""
    print u" -m --mods OPTION      Convert post-translationally modified amino acids (PTMs) to regular amino"
    print u"                       acids. These include PTR, NEP, HIP, SEP, TPO, AYA, MLZ, ... (see pdbatoms.py)"
    print u"                       Options: C = convert and remove [default], converts PTM amino HETATM to ATOM"
    print u"                                    and removes non-aa atoms from them."
    print u"                                A = only convert from HETATM to ATOM but do not remove any atoms."
    print u"                                N = none, make no changes to original PDB."
    print u" -t --threads  INT     Specify the number of threads to use for calculations."
    print u"                       Default value is the same number as detected system cores."
    print u"                       With negative values, the number of cores is decreased from maximum cores."
    try:    print u"                       Default value: %i (detected number of system cores)." % cpu_count()
    except: print u"                       Error: Could not detect number of cores."
    print u""
    print u" -i --icosahedron MIN  Use a geometric method (icosahedron) for generating an even sphere"
    print u"                       sample point distribution."
    print u"                       MIN speficies the minimum number of sample points per atom."
    print u"                       Possible numbers of sample points per atom: 12, 42, 162, 642, 2562,"
    print u"                       10242, 40962, 163842, 655362, ...  [default: 162]"
    print u" -g --golden POINTS    Use the Vogel method (golden spiral) to generate a sample sphere"
    print u"                       point distribution."
    print u""
    print u"                       POINTS specifies the exact number of points that will be"
    print u"                       distributed on a spherical surface."
    print u""
    print u" -v --verbose          Print more detailed output during execution."
    print u" -l --license          Print out license information."
    print u" -h --help             Print out usage instructions."
    print u""
    print u"Results of the calculation will be reported in square Ångströms in atom-mode (-a) and"
    print u"as percentages in residue-mode (-r), and written to the b-factor columns (temperature factor,"
    print u"cols 61-66) of the outputted PDB file. Using less spherical sample points will speed up"
    print u"calculation time, but lower accuracy."


def license():
    print """
    This file is part of the StructureMapper algorithm.
    Please cite the authors if you find this software useful
    
    https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361
    
    MIT License
    
    Copyright 2018 Anssi Nurminen and Vesa P. Hytönen
    Faculty of Medicine and Life Sciences and BioMediTech, University of Tampere, Arvo Ylpön katu 34, 33520 Tampere, Finland 
    Fimlab Laboratories, Biokatu 4, 33520 Tampere, Finland
    Tampere, Finland
    16-02-2018

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation the
    rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all 
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    """


def printDictionary( dict):

    print ""

    for k in sorted( dict.iterkeys()):
        print k,
        print ":",
        print dict[ k]

    print ""

def CalcModeFromStr( aStr):
    if aStr == "sidechain": return "S"
    elif aStr == "atom": return "A"
    return "R"

def CountMatchingArgs( aList1, aList2):
    count = 0
    if len( aList1) == 0 or len( aList2) == 0: return 0
    l2 = zip( *aList2)[ 0]
    for e in aList1:
        if e in l2: count += 1
    return count

def main():

    global globvar

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:g:t:c:p:n:m:dksearhvl', ['icosahedron=', 'vogel=', 'threads=','chains=', 'probe=', 'include=', 'mods=', 'sidechain', 'series', 'atom', 'residue', 'help', 'verbose'])
    except getopt.GetoptError as err:
        sys.stderr.write( err.msg)
        sys.stderr.write( "\n")
        sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)

    #Input & Output files
    for arg in args:
        if len( globvar[ "input_file"]) == 0:
            globvar[ "input_file"] = arg
        elif len( globvar[ "output_file"]) == 0:
            globvar[ "output_file"] = arg
        else:
            sys.stderr.write("Too many arguments.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)

    exclusives = [ '-e', '--series', '-d', '--directory']
    if CountMatchingArgs( exclusives, opts) > 1:
        sys.stderr.write("ERROR: Please specify only one option:'-e' or '-d'.\n")
        sys.exit( -4)

    exclusives = [ '-a', '--atom', '-r', '--residue', '-s', '--sidechain']
    if CountMatchingArgs( exclusives, opts) > 1:
        sys.stderr.write("ERROR: Please specify only one option '-a', '-r' or '-s'.\n")
        sys.exit( -4)

    #Flags
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit( 2)
        elif opt in ('-v', '--verbose'):
            print "Verbose mode on"
            globvar["verbose"] = True
        elif opt in ('-e', '--series'):
            globvar["filemode"] = "series"
        elif opt in ('-d', '--directory'):
            globvar["filemode"] = "directory"
        elif opt in ('-k', '--skip'):
            globvar["skip_existing"] = True,
        elif opt in ('-l', '--license'):
            license()
            sys.exit( 0)
        elif opt in ('-i', '--icosahedron'):
            globvar["distribution_type"] = "i"
            ok = True
            try:
                min_samples = int( arg)
            except ValueError:
                ok = False

            if not ok or min_samples < 1 or min_samples > 10000000:
                sys.stderr.write( "Bad minimum amount of sample points: '%s'.\n" % arg)
                min_samples = spherepoints.NumberOfIcosahedricPoints( 4)

            #get appropriate multiplier
            multiplier = 1
            while min_samples > spherepoints.NumberOfIcosahedricPoints( multiplier):
                multiplier += 1

            globvar["distribution_parameter"] = multiplier

            if min_samples != spherepoints.NumberOfIcosahedricPoints( multiplier):
                #If used number of points is not exactly the requested number, report it
                sys.stdout.write( "Using %i evenly distributed sample points per atom.\n" % spherepoints.NumberOfIcosahedricPoints( multiplier))


        elif opt in ('-g', '--golden'):
            globvar["distribution_type"] = "g"
            try:
                globvar["distribution_parameter"] = int( arg)
            except ValueError:
                globvar["distribution_parameter"] = 920
                sys.stderr.write( "Bad number of sample points: '%s'.\nUsing a 920-point Vogel sample point distribution.\n" % arg)

            if globvar["distribution_parameter"] < 6:
                sys.stderr.write( "Minimum possible number of point for Vogel distribution is 6.\n")
                globvar["distribution_parameter"] = 6

        elif opt in ('-c', '--chains'):
            globvar["chains"] = arg

        elif opt in ('-p', '--probe'):
            try:
                globvar["probe"] = float( arg)
            except ValueError:
                globvar["probe"] = 1.5
                sys.stderr.write( u"Bad probe radius: '%s'. Using default probe radius: 1.5Å.\n" % arg)

        elif opt in ('-n', '--include'):
                globvar["include"] = []
                accepted_vals = ['A','D','I','U','O']
                for a in arg:
                    au = a.upper()
                    if au in accepted_vals and au not in globvar["include"]:
                        globvar["include"].append( au)
                    else:
                        sys.stderr.write( u"Bad include type ignored: '%s'.\n" % a)

                if len( globvar["include"]) == 0:
                    sys.stderr.write("No valid residue types specified. Exiting.\n")
                    sys.exit( 0)
        elif opt in ('-m', '--mods'):
                globvar["mods"] = []
                accepted_vals = ['C','A','N']
                for a in arg:
                    au = a.upper()
                    if au in accepted_vals and au not in globvar["mods"]:
                        globvar["mods"].append( au)
                    else:
                        sys.stderr.write( u"Bad mod type ignored: '%s'.\n" % a)

                if len( globvar["mods"]) == 0:
                    sys.stderr.write("No valid mod types specified. Exiting.\n")
                    sys.exit( 0)
        elif opt in ('-t', '--threads'):

            #Using actually processes not threads (but still callilng them threads).
            #0 for default; the number of cores on system
            try:
                globvar["threads"] = int( arg)

                if globvar["threads"] <= 0:

                    requested_threads = cpu_count() + globvar["threads"]
                    if requested_threads < 1:
                        sys.stderr.write( "WARNING: Bad thread count specified: '%s'. Using only a single thread.\n" % arg)
                        globvar["threads"] = 1
                    else:
                        globvar["threads"] = requested_threads
                        sys.stderr.write( "INFO: Using %i cores.\n" % globvar["threads"])

            except ValueError:
                sys.stderr.write( "WARNING: Bad thread count specified: '%s'. Using only a single thread.\n" % arg)
                globvar["threads"] = 1

        elif opt in ('-a', '--atom'):
            globvar["ASA_level"] = "atom"
        elif opt in ('-r', '--residue'):
            globvar["ASA_level"] = "residue"
        elif opt in ('-s', '--sidechain'):
            globvar["ASA_level"] = "sidechain"
            sys.stderr.write("INFO: Mode: sidechains only\n")
        else:
            sys.stderr.write("Bad arguments specified.\n")
            sys.stderr.write("See -h for usage.\n")
            sys.exit( 2)

    if globvar["verbose"]:
        printDictionary( globvar)


    if len( globvar[ "input_file"]) == 0:
        sys.stderr.write("No input file specified.\n")
        sys.stderr.write("See -h for usage.\n")
        sys.exit( 2)


    if globvar["filemode"] == "single":
        if len( globvar[ "output_file"]) == 0:
            globvar[ "output_file"] = InsertASAToFilePath( globvar[ "input_file"])
            sys.stdout.write("Writing output to file: '%s'\n" % globvar[ "output_file"])
        elif os.path.isdir( globvar[ "output_file"]):
            #Use default naming convention
            globvar[ "output_file"] = InsertASAToFilePath( ChangeToFolder( globvar[ "input_file"], globvar[ "output_file"]))


    try:

        #RUN
        if globvar["filemode"] == "series":
            runASASeries( globvar[ "input_file"], globvar[ "chains"], globvar[ "output_file"], globvar[ "threads"], globvar[ "skip_existing"], calc_mode=CalcModeFromStr( globvar["ASA_level"]), modifications=globvar["mods"])
        elif globvar["filemode"] == "directory":
            runASAFolder( globvar[ "input_file"], globvar[ "chains"], globvar[ "output_file"], globvar[ "threads"], globvar[ "skip_existing"], calc_mode=CalcModeFromStr( globvar["ASA_level"]), modifications=globvar["mods"])
        else:
            runASA( globvar[ "input_file"], globvar[ "chains"], globvar[ "output_file"], globvar[ "threads"], calc_mode=CalcModeFromStr( globvar["ASA_level"]), modifications=globvar["mods"] )

    except KeyboardInterrupt:

        sys.stderr.write( "\n\nProcess interrupted by user.\n\n")




if __name__ == "__main__":
  main()

 #Example usage:
 #pdb = pdbatoms.ReadPDB( "pdb3ikm.ent")
 #residues = []
 #residues.append( pdbatoms.FindResidue( pdb["residues"], aResNum="458", aChain="D", aModel=-1))
 #residues.append( pdbatoms.FindResidue( pdb["residues"], aResNum="264", aChain="F", aModel=-1))
 #residues.append( pdbatoms.FindResidue( pdb["residues"], aResNum="644", aChain="D", aModel=-1))

 #for r in residues: r.PrintCA()

 ##def CalcASAForResidues( aResidues, aPdb, aProbe_r=1.5, aChains=[], aModels=[], aIncluded_atom_types=["A","D"] ):
 #print CalcASAForResidues( residues, pdb, aChains=[], aCalcMode="R")

 #for r in residues: r.PrintCA()

