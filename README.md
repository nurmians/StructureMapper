# StructureMapper

The StructureMapper algorithm processes amino acid sequences of proteins with marked points of interest and finds and analyzes associated structural data. An online version of the algorithm is available at [structuremapper.uta.fi](http://structuremapper.uta.fi).

The offline version of the algorithm is designed to process points of interest in a high-throughput manner, effectively utilizing parallellization.

If you find the StructureMapper algorithm useful please cite the authors and the associated article:
"StructureMapper: a high-throughput algorithm for analyzing protein sequence locations in structural data", Bioinformatics 2018, open access.

https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361

## Installing the algorithm

Only Python 2.7 and [BLAST](https://blast.ncbi.nlm.nih.gov) installations are required for running the algorithm locally. For faster processing, it is recommended that additional executables are also installed locally. The algorithm has been run and tested on Windows and Linux operating systems.

## Python

Install Python first. The algorithm is written using Python 2.7.  
**Linux:** ```sudo apt-get install python2.7```  
**Windows:** https://www.python.org/downloads/windows/  

The biopython module might require also the package libpython2.7-dev  
**Linux:** ```sudo apt-get install libpython2.7-dev```  

Install required python modules using pip:  

```
pip install biopython  
pip install numpy  
```

## BLAST

Setup a local BLAST DB (<50MB) and the blast executable.

Instructions can be found [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

Or simply download the suitable executable files for your system:

Blast binaries (all operating systems):
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

The path of the BLAST executable folder ('bin') can be added to the system path or provided with the 'DEFAULT_BLASTP' param in the config.ini file (see below). Only BLAST database 'pdbaa' is required for running StructureMapper. The pdbaa database contains sequences found in the [PDB](http://www.rcsb.org).

### Downloading the BLAST pdbaa database:

First, create and enter the folder use wish to use for the database files.

**Linux:**
Create a folder and download and extract the file:  
```wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz```  
```tar -zxvf pdbaa.tar.gz```  

**Windows:**
You can use the update_blastdb.pl perl script that comes with the BLAST binaries package and, 7-zip (for example) to extract the packages (requires [Perl](https://www.perl.org/get.html#win32) and [7-zip](http://www.7-zip.org) to be installed).

```
perl update_blastdb.pl pdbaa
"C:\Program Files\7-Zip\7z.exe" e -aoa pdbaa.tar.gz  
"C:\Program Files\7-Zip\7z.exe" x -aoa pdbaa.tar
```  

A config.ini file is created on first run of the StructureMapper algorithm and it needs to be modified to contain the paths to the BLAST DB and executable. Example config.ini lines:   
```DEFAULT_BLASTP = C:/blast/blast-2.7.1+/bin```   
```BLAST_DB = C:/blast/DB/```

You can also try running StructureMapper with the --autoconfig flag. The structures required for the analysis are downloaded if and when they are needed. For large datasets, the processing can be speeded up by downloading the PDB structure files in advance and by using precalculated ASA files.

## DSSP

The algorithm uses DSSP for secondary structure evaluations. 
Available at: ftp://ftp.cmbi.ru.nl/pub/software/dssp/  
**Linux:** ```wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64```  
make sure to set ```chmod a+x dssp-2.0.4-linux-amd64``` to be able to run it.  

The algorithm can also download the dssp files if the DSSP algorithm is not installed locally.

## Disorder predictions

If you are interested in disordered regions you can download Iupred from http://iupred.enzim.hu and specify the executable in the algorithm's config file (config.ini, created on first run).

## Configuration

The algorithm is run by executing "python score_poisites.py" in the src/main directory (run with "python2.7 score_poisites.py" if python 3 is set as the default).
Example usage (run.bat):  
```python src/main/score_poisites.py -b 3 -t 1 src/main/test/test_input.fasta ./results```

A config.ini file is created on the first run, and it needs to be modified to contain all the required 
paths, such as BLAST (executable and DB).

## Input

StructureMapper takes as its input FASTA formatted amino acid sequence files (.fasta) with asterisks ('\*') as point of interest (POI) markers. The input file can contain multiple sequences and each sequence can contain multiple POIs. Each sequence should have a unique identifier in its header.   
Example fasta sequence header:  
>\>sp|**O75838**|CIB2_HUMAN Calcium and integrin-binding family member 2 OS=Homo sapiens GN=CIB2 PE=1 SV=1
