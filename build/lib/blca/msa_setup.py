import os
import sys
from .helpers import *
from Bio import SearchIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline as muscle
from .helpers import my_module

def get_hit_seq_megan(filename):
	yamlfile = yaml_load_file()
	blout = SearchIO.parse(filename, 'blast-text')
	for query in blout:
		seqid = query.id.split("\n")[0]
		seqid = seqid
		#print(seqid)
		fh = open("multi_" + seqid + ".fasta", 'a')
		yamlfile[seqid]['hits'] = {}
		bcp = 1 - (my_module.BLAST_CUTOFF_PERCENT / 100)
		topscore = 0
		for hit in query.hits:
			gi = re.match(r"gi\|(.*)\|ref", hit.id).group(1)
			yamlfile[seqid]['hits'][gi] = {}
			#print(yamlfile[seqid]['hits'])
			for hsp in hit.hsps:
				if topscore == 0:
					topscore = hsp.bitscore
				if hsp.bitscore < (topscore * bcp):
					#print(seqid, " not included: ", gi, " score: ", str(hsp.bitscore), " topscore: ", topscore)
					continue
				if hsp.bitscore < my_module.BLAST_CUTOFF_SCORE:
					continue
				if int(100 * hsp.ident_num/hsp.aln_span) < my_module.BLAST_CUTOFF_PERCENT:
					continue
				if int(100 * (float(hsp.query_end - hsp.query_start + 1)/query.seq_len)) < my_module.BLAST_COVERAGE:
					continue
				#print(hsp.hit, "\t", hsp.hit_strand, "\t", hsp.hit_start, "\t", hsp.hit_end, "\t", hsp.bitscore)
				#print(hsp.query_end, "\t", hsp.query_start, "\t", query.seq_len , "\t", int(100 * (float(hsp.query_end - hsp.query_start + 1)/query.seq_len)))
				hitstart = hsp.hit_start + 1 - 10
				hitstart = 1 if hitstart <= 0 else hitstart
				hitend = hsp.hit_end + 1 + 10
				hitstrand = "plus" if (hsp.hit_strand == 1) else "minus"
				#print(hitstart, hitend)
				out = os.popen(my_module.BLAST_BINARY + "/blastdbcmd -db " + my_module.BLAST_DATABASE + " -dbtype nucl -entry " + str(gi) + " -range " + str(hitstart) + "-" + str(hitend) + " -strand " + str(hitstrand)).read()
				#print("out", out)
				fh.write(out)
				#print(hsp.hit.seq)
				#print(hsp.query.seq)
				#print("-----")
		fh.close()
		#break
	yaml_dump_file(yamlfile)

def setup_msa():
	seqlist = {}
	#seqs = []
	print("INFO: Creating MSA input files")
	for seq in SeqIO.parse(my_module.FILENAME, "fasta"):
		#seqs.append(seq.id)
		seqlist[seq.id] = {}
		fh = open("multi_" + seq.id + ".fasta", 'w')
		fh.write(">" + seq.id + "\n")
		fh.write(str(seq.seq) + "\n")
		fh.close()
	yaml_dump_file(seqlist)
	get_hit_seq_megan(my_module.FILENAME + ".blastn")
	print("DONE: MSA input files generated.")
	#return seqs

def run_muscle():
	yamlfile = yaml_load_file()
	tot = len(yamlfile.keys())
	count = 1
	print("INFO: Running MUSCLE for MSA")
	for seqid in yamlfile:
		#print(seqid)
		sys.stdout.write("Files: %d of %d   \r" % (count, tot))
		count += 1
		filename= "multi_" + seqid + ".fasta"
		outfile = filename + ".maln"
		muscle_cline = muscle(cmd=my_module.MUSCLE_BINARY, input=filename, out=outfile)()
		#print(outfile)
		sys.stdout.flush()
	print("DONE: MSA complete.")
