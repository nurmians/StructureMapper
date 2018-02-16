# StructureMapper

The StructureMapper algorithm processes amino acid sequences of proteins with marked points of interest and finds and analyzes associated structural data. An online version of the algorithm is available at http://structuremapper.uta.fi

The offline version of the algorithm is designed to process points of interest in a high-throughput manner, effectively utilizing parallellization.

If you find the StructureMapper algorithm useful please cite the authors and the associated article:
"StructureMapper: a high-throughput algorithm for analyzing protein sequence locations in structural data", Bioinformatics 2018, open access.

https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361

# Installing the algorithm

Setting up the algorithm for running it locally on either windows or a linux machine (Mac not tested).

# Python

Install Python first. The algorithm is written using Python 2.7.
(sudo apt-get install python2.7)

the biopython module might require also the package libpython2.7-dev
(sudo apt-get install libpython2.7-dev)

Install required python modules using pip:

pip install biopython

pip install numpy

# BLAST

Setup a local BLAST DB (<50MB) and the blast executable.

Instructions can be found here:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Or simply download the suitable executable files for your system:

Binaries:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Linux:
Create a folder and download and extract the file
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz
tar -zxvf pdbaa.tar.gz

Windows:
You can use the update_blastdb.pl perl script that comes with the BLAST binaries packaage and for example 7-zip to extract the packages

perl update_blastdb.pl pdbaa
"C:\Program Files\7-Zip\7z.exe" e -aoa pdbaa.tar.gz
"C:\Program Files\7-Zip\7z.exe" x -aoa pdbaa.tar

The homologous structures are downloaded if and when they are needed. For large datasets, the processing can be speeded up by downloading the PDB structure files and precalculating ASA and ASA_SC files.

# DSSP

The algorithm uses DSSP for secondary structure evaluations. 
Available at: ftp://ftp.cmbi.ru.nl/pub/software/dssp/
(wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64)
make sure to set chmod a+x dssp-2.0.4-linux-amd64 to be able to run it

The algorithm can also dowload the dssp files if no executable is provided

# Disorder predictions

If you are interested in disordered regions you can download Iupred from http://iupred.enzim.hu and specify the executable in the algorithm's config file (config.ini, created on first run).

# Configuration

The algorithm is run by running "score_poisites.py" in the src/main directory.
Example usage: "python score_poisites.py-t -1 -s -p 1 -b 4 "C:\myinputfile.fasta" ..\results\myinput_results"

A config.ini file is created on the first run, and it needs to be modified to contain all the required 
paths BLAST (executable and DB).
