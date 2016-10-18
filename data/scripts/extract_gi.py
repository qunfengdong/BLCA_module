import re

with open("16SMicrobialseq.fasta.names", "r") as fh:
  for line in fh:
    line = line.rstrip()
    gi = re.findall(r'>gi\|([0-9]*)\|ref\|.*', line)
    print(gi[0])
