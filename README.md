Bayesian-based LCA taxonomic classification method
--------------------------------------------------

## System Requirements
* Python (version 3 or above)
* Linux
* libyaml

## Python packages
* python3-dev
* setuptools
* Biopython
* pyyaml

## Additional executable files (required)
* BLAST binary (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/)
* MUSCLE (http://www.drive5.com/muscle/downloads.htm)

## Database (required)
* 16S Microbial database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz)

## Install
Checkout the source code: https://github.com/qunfengdong/BLCA
```python
$ git clone https://github.com/qunfengdong/BLCA.git
$ cd BLCA
$ python3 setup.py install
```

## Getting started
Go to the working directory, save the input FASTA file and config.py file in the same folder.

## Config file
Create config file "config.py" in the current folder to setup the parameters
```python
[terminal]$ cat config.py
##-----------------------------------------------------------
## Input file in FASTA format
FILENAME = 'inputfile.fasta'

## Output file
OUTFILE = 'inputfile.annotation.txt'

## BLAST parameters
BLAST_BINARY = 'path/to/dir/ncbi-blast-2.3.0+/bin'
## Path to the actual indexed BLAST database with 16SMicrobial.n* files
BLAST_DATABASE = 'path/to/file/16SMicrobial'
## MUSCLE executable file
MUSCLE_BINARY = '/path/to/bin/muscle'

## BLAST output filtering options
## Select hits with at least a score 100
BLAST_CUTOFF_SCORE = 100
## Select hits within less then 10% score of the top hit
BLAST_CUTOFF_PERCENT = 10
## Select hits with at least 95% coverage
BLAST_COVERAGE = 95
## Select hits with at least 95% percentage identity
BLAST_PERCENTAGE_IDENTITY = 95

## Mutiple sequence alignment parameters
## Select 10 base pairs of hit sequences
## upstream and downstream of the alignment
HIT_SEQUENCE_BPS = 10

## Alignment scoring values
ALIGNMENT_GAP = -2.5
ALIGNMENT_MISMATCH = -2
ALIGNMENT_MATCH = 1

## Bootstrap
BOOTSTRAP = 100
##-----------------------------------------------------------
```

## Run BLCA on command line
```python
[terminal]$ python3
>>> import blca
>>> import config
>>> blca.execute()
```

## Run BLCA by creating a file called analysis.py
```python
import blca
import config
blca.execute()
```

Execute the file 'analysis.py' at the terminal
```python
[terminal]$ python3 analysis.py
```

## Output




