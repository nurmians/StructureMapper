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

import os
import sys
import glob
import re
import string

g_valid_filename_chars = "-_.()[]{}= %s%s" % (string.ascii_letters, string.digits)

def ClearFolder( aFolder, aClearSubfolders=False, aVerbose=False):

    subfolders = []
    if aVerbose: print "Clearing folder: '%s'" % aFolder

    for the_file in os.listdir( aFolder):
        file_path = os.path.join( aFolder, the_file)
        try:
            if os.path.isfile(file_path):
                if aVerbose: print "filepath:", file_path
                os.unlink(file_path)
            elif os.path.isdir( file_path):
                if aVerbose: print "subfolder:", file_path
                subfolders.append( the_file)
        except Exception, e:
            sys.stderr.write( "ERROR:" + str( e) + "\n")

    if aClearSubfolders:
        for subfolder in subfolders: ClearFolder( InsertFolderSeparator( InsertFolderSeparator( aFolder) + subfolder))


def PrependToFile( filename, line, newfilename=None):

    if newfilename == None:
        with open( filename,'r+') as f:
            content = f.read()
            f.seek( 0,0)
            f.write( line.rstrip('\r\n') + '\n' + content)
    else:
        with open( filename,'r') as source_file:
            with open( newfilename,'w') as dest_file:
                dest_file.write( line.rstrip('\r\n') + '\n')
                dest_file.write( source_file.read())

def InsertASAToFilePath( filepath ):

    path, filename = os.path.split( filepath)

    if filename.lower().find( ".asa.") >= 0: return filepath #RETURN

    name, suffix = SplitFileName( filename)
    return InsertFolderSeparator( path) + name + ".asa." + suffix


def IsASAFile( filename):
    if filename.lower().find( ".asa.") >= 0: return True
    return False


def ChangeToFolder( filepath, folder):

    path, filename = os.path.split( filepath)
    return InsertFolderSeparator( folder) + filename

def ConcatenateFiles( aFileList, aToFile="", aAppendIfExists=True, aRemoveAfterWritten=False ):

    if not len( aFileList): return

    outfile = None

    try:
        if not len( aToFile):
            outfile = sys.stdout
        elif aAppendIfExists:
            outfile = open( aToFile, 'a')
        else:
            outfile = open( aToFile, 'w')
    except:
        sys.stderr.write( "FILE ERROR: Could open file '%s' for writing." % aToFile)
        return

    filenames = aFileList[:]

    for fname in reversed( filenames):
        try:
            with open( fname) as infile:
                for line in infile: outfile.write( line)
        except:
            sys.stderr.write( "FILE ERROR: Could not read file '%s'. File skipped." % fname)
            filenames.remove( fname)

    outfile.close()

    if aRemoveAfterWritten:
        for fname in filenames:
            try: os.remove( fname)
            except: pass




def WriteListToFile( aList, aFile, aAppend=False):

    if aFile == None or len( aFile) == 0: return

    if len( aList) == 0 and aAppend: return #nothing to do

    try:
        if len( aList) == 0 and not aAppend:
                open(aFile, 'w').close() #clear file
        if not aAppend:
            with open( aFile, "w+") as f:   #reset
                f.write( "\n".join( aList) + "\n") #newline after last elem is important
        else:
            with open( aFile, "a+") as f:   #append
                f.write( "\n".join( aList) + "\n") #newline after last elem is important
    except:
        sys.stderr.write( "FILE ERROR: Could not write list to file '%s'.")

def ReadListFromFile( aFile, aRemoveDuplicates=False):
    lines = []

    if not os.path.isfile( aFile): return lines

    try:
        with open( aFile, "r") as f:
            lines = f.readlines()
    except:
        sys.stderr.write( "FILE ERROR: Could not read list from file '%s'." % aFile)

    if len( lines):
        lines = map( str.strip, lines)
        if lines[ -1] == "": del lines[ -1] #Last newline can cause an empty last element

        if aRemoveDuplicates: lines = list( set( lines))

    return lines



def SortFilesByNumber( file_list ):

    regex = re.compile( r'(\d+)(?!.*\d)')

    numvals = []
    errors = False

    for file in sorted( file_list):

        find = regex.findall( file)

        if len( find) < 1:
            numvals.append ( -1)
            errors = True
            continue

        try:
            num = int( find[ 0])
            numvals.append( num)
        except ValueError:
            numvals.append ( -1)
            errors = True
            continue

    if errors:
        sys.stderr.write( "Warning: Files could not be sorted without errors.\n")

    sorted_files = [x for (y,x) in sorted(zip(numvals,file_list))]
    return sorted_files


def RemoveExistingASAfilesFromPDBFileList( file_list, folder=""):

    skipped_asa_files = []

    i = 0
    while i < len( file_list):

        file = file_list[ i]
        if len( folder): file = ChangeToFolder( file, folder)

        asa_file = InsertASAToFilePath( file)
        #print "ASA_file:" + asa_file

        if not os.path.isfile( asa_file):
            #file does not exist, calculate
            i += 1
            continue

        file_size = int( os.stat( asa_file).st_size)

        if file_size == 0:
            #file empty, recalculate
            i += 1
            continue

        #file exists, remove from list
        sys.stdout.write( "Info: File '%s' already exists. Skipping.\n" % asa_file )
        skipped_asa_files.append( asa_file)
        del file_list[ i]

    return skipped_asa_files


def SplitFilePath( filepath ):

    path, filename = os.path.split( filepath)
    name, suffix = SplitFileName( filename)

    return (InsertFolderSeparator( path), filename, name, suffix)


def GetFileNameWithoutPathAndSuffix( aFilePath):
    a,b,c,d = SplitFilePath( aFilePath)
    return c

def SplitFileName( filename ):

    suffix = ''
    name = ''
    fileparts = filename.split( '.')

    l = len( fileparts)

    if l == 0:
        pass
    elif l == 1:
        name = fileparts[ 0]
    elif l == 2:
        name = fileparts[ 0]
        suffix = fileparts[ 1]
    else:
        name = fileparts[ 0]
        suffix = '.'.join( fileparts[1:])

        #suffix = fileparts[ -1]
        #name = ''.join( fileparts[:-1])

    return (name, suffix)


def InsertFolderSeparator( folder):

    if len( folder) == 0: return folder

    folder_sep = "/"
    if folder.find( "\\") >= 0: folder_sep = "\\"

    #Make sure folder ends with a folder separator character
    if folder[ -1] != folder_sep:
        folder += folder_sep

    return SystemSpecificFolderSeparators( folder)

def DriveLetter( aFilePath=None):
    if os.name == 'nt':
        if not aFilePath or len( aFilePath) == 0: aFilePath = __file__
        return os.path.abspath( aFilePath)[:3]
    return "/"


def TryInt( a):
    try: return int( a)
    except: return a

def AlphanumKey( a):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ TryInt( c) for c in re.split('([0-9]+)', a) ]

def NaturalSort( a):
    """ Sort the given list in the way that humans expect.
    """
    a.sort( key=AlphanumKey)

#path can include wild cards *
def GetDirFiles( path):

    filelist = []

    if path.find( "*") < 0 and os.path.isdir( path):
        path = InsertFolderSeparator( path) + "*"

    for globfile in glob.glob( path):
        filelist.append( globfile)

    NaturalSort( filelist)
    return filelist


def GetFolderFiles( aFolder, aWildCard="*"):

    #import glob

    if len( aFolder) < 1:
        sys.stderr.write( "No folder specified.\n")
        return (-1, [])

    #aFolder = LinuxFolderSeparators( aFolder)
    is_dir = os.path.isdir( aFolder)

    if not is_dir:
        sys.stderr.write( "Argument '%s' is not an existing folder.\n" % aFolder )
        return (-2, [])

    #Make sure folder ends with a folder separator character
    aFolder = InsertFolderSeparator( aFolder)

    file_path = aFolder + aWildCard

    file_list = []

    #print  "Globbing: '%s'" % file_path
    #Get PDB filenames
    for globfile in glob.glob( file_path):
        file_list.append( globfile)

    return (0, file_list)

def SystemSpecificFolderSeparators( aPath):

    if os.name == 'nt': return aPath.replace( "/", "\\")
    return aPath.replace( "\\", "/")



def ListOfFolderFiles( folder, aPattern=""):

    #import glob

    if len( folder) < 1:
        sys.stderr.write( "No folder specified.\n")
        return (-1, [])

    is_dir = os.path.isdir( folder)

    if not is_dir:
        sys.stderr.write( "Argument %s is not an existing folder.\n" % folder )
        return (-2, [])

    #Make sure folder ends with a folder separator character
    folder = InsertFolderSeparator( folder)


    file_list = []
    file_path = folder



    if len( aPattern): #Custom pattern

        file_path = folder + aPattern

        file_list = []

        #Get PDB filenames
        for globfile in glob.glob( file_path):
            if IsASAFile( globfile): continue #Do not list asa files
            file_list.append( globfile)

    else:

        file_path = folder + "*.pdb"

        #Get PDB filenames
        for globfile in glob.glob( file_path):
            if IsASAFile( globfile): continue #Do not list asa files
            file_list.append( globfile)

        #PDB files downloaded from the PDB database use suffix .ent
        file_path = folder + "*.ent"

        #Get ent filenames
        for globfile in glob.glob( file_path):
            if IsASAFile( globfile): continue #Do not list asa files
            file_list.append( globfile)


    return (0, file_list)

def GetFilename( filepath):
    path, filename = os.path.split( filepath)
    return filename

def GetFolder( filepath):
    path, filename = os.path.split( filepath)
    return InsertFolderSeparator( path)

def ListOfSeriesFiles( file):

    #file = file.replace( ".asa.", "")
    is_asa_file = IsASAFile( file)

    #check whether "/" or "\" is used
    folder_sep = "/"
    if file.find( "\\") >= 0: folder_sep = "\\"

    path, filename = os.path.split( file)
    name, suffix = SplitFileName( filename)


    if len( name) == 0:
        return []

    #Find the last number in filename
    regex = re.compile( r'(.*?)(\d+)(?!.*\d)')


    find = regex.findall( name)

    if len( find) < 1 or len( find[ 0]) < 2:
        sys.stderr.write( "Filename '%s.%s' does not contain a series number.\n" % (name, suffix) )
        return ("", []) #RETURN

    first_num_of_series = -1
    basename = find[ 0][ 0]

    try:
        first_num_of_series = int( find[ 0][ 1])
    except ValueError:
        stderr.write( "Filename '%s.%s' does not contain a valid series number.\n" % (name, suffix) )
        return ("", []) #RETURN


    numvals = []
    retvals = []
    regex2 = re.compile( r'(\d+)(?!.*\d)')

    skipped = 0



    file_path = ""
    if len( path): file_path += path + folder_sep
    file_path += basename + "*"
    if len( suffix): file_path += "." + suffix

    #print "GLOB:" + file_path
    #print "suffix:" + suffix
    #print glob.glob( file_path)


    for globfile in glob.glob( file_path):

        find = regex2.findall( globfile)

        if len( find) < 1:
            skipped += 1
            continue
        try:
            num = int( find[ 0])
        except ValueError:
            skipped += 1
            continue

        if num < first_num_of_series:
            skipped += 1
            continue
        elif num in numvals:
            sys.stderr.write( "Warning: file number %i in series found multiple times. Skipping file: %s\n" % (num, globfile) )
            skipped += 1
            continue

        #Do not list asa files if first file of the series is not an asa file
        #and vice versa
        if not is_asa_file and IsASAFile( globfile): continue
        elif is_asa_file and not IsASAFile( globfile): continue

        retvals.append( globfile)
        numvals.append( num)


    sorted_files = [x for (y,x) in sorted(zip(numvals,retvals))]
    return (basename, sorted_files)



def GetFastaFilesInFolder( folder):

    folder = InsertFolderSeparator( folder)
    glob_path1 = folder + "*.fasta"

    file_list = []

    for globfile in glob.glob( glob_path1):
        file_list.append( globfile)


    return file_list

def LinuxFolderSeparators( aPath):
    return aPath.replace( "\\", "/")

#Global
#g_valid_filename_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)

def ValidFilename( aFilepath, aReplacement='_'):

    global g_valid_filename_chars

    path, filename, name, suffix = SplitFilePath( aFilepath)

    retval = ''.join(c if c in g_valid_filename_chars else aReplacement for c in filename)
    return (path + retval)

