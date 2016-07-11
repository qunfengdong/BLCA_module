##-----------------------------------------------------------
## Input file in FASTA format
FILENAME = 'example_input_file.fasta'

## Output file
OUTFILE = 'annotation_output.txt'

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

## Cutoff 
CUTOFF = 0

##-----------------------------------------------------------

