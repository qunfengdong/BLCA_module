Bayesian-based LCA taxonomic classification method
--------------------------------------------------

## Requirements
* Python3

## Install
Checkout the source code: https://github.com/qunfengdong/BLCA
```python
python3 setup.py install
```
## Config file
Create config file "config.py" in the current folder to setup the parameters
```python
## Input file
FILENAME = '16SMicrobial_150.V4.subsample3.fasta'

## Output file
OUTFILE = '16SMicrobial_150.V4.subsample3.fasta.annotation'

## BLAST parameters
BLAST_BINARY = 'data/ncbi-blast-2.3.0+/bin'
BLAST_DATABASE = 'data/16SMicrobial/16SMicrobial'
BLAST_CUTOFF_SCORE = 100
BLAST_CUTOFF_PERCENT = 10
BLAST_COVERAGE = 95
BLAST_PERCENTAGE_IDENTITY = 95

## MSA parameters
HIT_SEQUENCE_BPS = 10
MUSCLE_BINARY = '/usr/bin/muscle'

## Scores
ALIGNMENT_GAP = -2.5
ALIGNMENT_MISMATCH = -2
ALIGNMENT_MATCH = 1

## Bootstrap
BOOTSTRAP = 100
```

## Run BLCA
Run BLCA in another file
```
import blca
import config

blca.runblca()
```

