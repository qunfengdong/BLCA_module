# 2016_10_12
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz
tar zxvf 16SMicrobial.tar.gz

## Extract the fasta sequences
../ncbi-blast-2.3.0+/bin/blastdbcmd -entry all -db 16SMicrobial -dbtype nucl > 16SMicrobialseq.fasta
grep "^>" 16SMicrobialseq.fasta > 16SMicrobialseq.fasta.names
python3 extract_gi.py > 16SMicrobialseq.fasta.names.gi

## Get taxonomy info
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
unzip gi_taxid_nucl.zip

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxvf taxdump.tar.gz

## Parse the taxonomy files
python3 taxonomy.py
mv subset* ../.
