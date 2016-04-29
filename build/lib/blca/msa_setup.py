import os
import sys
from .helpers import *
from Bio import SearchIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline as muscle

def get_hit_seq(fastafile, filename):
	yamlfile = yaml_load_file('storage.yaml')
	blout = SearchIO.parse(filename, 'blast-text')
	for query in blout:
		seqid = query.id.split("\n")[0]
		#print(seqid)
		fh = open("multi_" + seqid + ".fasta", 'a')
		yamlfile[seqid]['hits'] = {}
		for hit in query.hits:
			gi = re.match(r"gi\|(.*)\|ref", hit.id).group(1)
			yamlfile[seqid]['hits'][gi] = {}
			#print(yamlfile[seqid]['hits'])
			for hsp in hit.hsps:
				#print(hsp.hit)
				#print(hsp.hit_strand)
				#print(hsp.hit_start)
				#print(hsp.hit_end)
				hitstart = hsp.hit_start + 1 - 10
				hitstart = 1 if hitstart < 0 else hitstart
				hitend = hsp.hit_end + 1 + 10
				hitstrand = "plus" if (hsp.hit_strand == 1) else "minus"
				#print(hsp.hit_end)
				out = os.popen("data/ncbi-blast-2.3.0+/bin/blastdbcmd -db data/16SMicrobial/16SMicrobial -dbtype nucl -entry " + str(gi) + " -range " + str(hitstart) + "-" + str(hitend) + " -strand " + str(hitstrand)).read()
				fh.write(out)
				#print(hsp.hit.seq)
				#print(hsp.query.seq)
				#print("-----")
		fh.close()
		#break
	yaml_dump_file(yamlfile)

def get_hit_seq_megan(fastafile, filename):
	yamlfile = yaml_load_file('storage.yaml')
	blout = SearchIO.parse(filename, 'blast-text')
	for query in blout:
		seqid = query.id.split("\n")[0]
		seqid = seqid
		#print(seqid)
		fh = open("multi_" + seqid + ".fasta", 'a')
		yamlfile[seqid]['hits'] = {}
		topscore = 0
		for hit in query.hits:
			gi = re.match(r"gi\|(.*)\|ref", hit.id).group(1)
			yamlfile[seqid]['hits'][gi] = {}
			#print(yamlfile[seqid]['hits'])
			for hsp in hit.hsps:
				if topscore == 0:
					topscore = hsp.bitscore
				if hsp.bitscore < (topscore * 0.95):
					#print(seqid, " not included: ", gi, " score: ", str(hsp.bitscore), " topscore: ", topscore)
					continue
				if hsp.bitscore < 100:
					continue
				#print(hsp.hit)
				#print(hsp.hit_strand)
				#print(hsp.hit_start)
				#print(hsp.hit_end)
				hitstart = hsp.hit_start + 1 - 10
				hitstart = 1 if hitstart <= 0 else hitstart
				hitend = hsp.hit_end + 1 + 10
				hitstrand = "plus" if (hsp.hit_strand == 1) else "minus"
				#print(hitstart, hitend)
				out = os.popen("data/ncbi-blast-2.3.0+/bin/blastdbcmd -db data/16SMicrobial/16SMicrobial -dbtype nucl -entry " + str(gi) + " -range " + str(hitstart) + "-" + str(hitend) + " -strand " + str(hitstrand)).read()
				fh.write(out)
				#print(hsp.hit.seq)
				#print(hsp.query.seq)
				#print("-----")
		fh.close()
		#break
	yaml_dump_file(yamlfile)

def setup_msa(fastafile):
	seqlist = {}
	#seqs = []
	print("INFO: Creating MSA input files")
	for seq in SeqIO.parse(fastafile, "fasta"):
		#seqs.append(seq.id)
		seqlist[seq.id] = {}
		fh = open("multi_" + seq.id + ".fasta", 'w')
		fh.write(">" + seq.id + "\n")
		fh.write(str(seq.seq) + "\n")
		fh.close()
	yaml_dump_file(seqlist)
	#get_hit_seq(fastafile, fastafile + ".blout")
	get_hit_seq_megan(fastafile, fastafile + ".blastn")
	print("DONE: MSA input files generated.")
	#return seqs

def run_muscle():
	yamlfile = yaml_load_file('storage.yaml')
	tot = len(yamlfile.keys())
	count = 1
	for seqid in yamlfile:
		#print(seqid)
		sys.stdout.write("Files: %d of %d   \r" % (count, tot))
		count += 1
		filename= "multi_" + seqid + ".fasta"
		outfile = filename + ".maln"
		muscle_cline = muscle(input=filename, out=outfile)()
		#print(outfile)
		sys.stdout.flush()
