import random
import os
import sys
import itertools
from .helpers import *
from Bio import SeqIO
import importlib
from importlib.machinery import SourceFileLoader
try:
	my_module = importlib.import_module('config')
except:
	this_dir, this_filename = os.path.split(__file__)
	my_module = SourceFileLoader("settings", this_dir + "/settings.py").load_module()

def setup_query(seq):
	for i in range(len(seq)):
		if(seq[i] != "-"):
			break
		seq = seq[:i] + "#" + seq[i+1:]
	for i in reversed(range(len(seq))):
		if(seq[i] != "-"):
			break
		seq = seq[:i] + "#" + seq[i+1:]
	return seq


def setup_query_pos(seq):
	for i in range(len(seq)):
		if(seq[i] != "-"):
			front = i
			break
		seq = seq[:i] + "#" + seq[i+1:]
	for i in reversed(range(len(seq))):
		if(seq[i] != "-"):
			back = i
			break
		seq = seq[:i] + "#" + seq[i+1:]
	return seq, front, back

def compute_score(seq1, seq2):
	score = 0
	for i in range(0,len(seq1)):
		#print(seq1[i], seq2[i])
		if ((seq1[i] == "#") or (seq2[i] == "#")):
			score += 0
		if ((seq1[i] == "-") or (seq2[i] == "-")):
			score += my_module.ALIGNMENT_GAP
		elif (seq1[i] != seq2[i]):
			score += my_module.ALIGNMENT_MISMATCH
		elif ((seq1[i] == "-") and (seq2[i] == "-")):
			score += 0
		elif (seq1[i] == seq2[i]):
			score += my_module.ALIGNMENT_MATCH
	return score

def calculate_prob(seqsdic, query, queryseq):
	scoredic = {}
	sumscore = 0
	for seqid in seqsdic:
		if seqid == query:
			continue
		score = compute_score(queryseq, seqsdic[seqid])
		scoredic[seqid] = {}
		scoredic[seqid]['score'] = score
		sumscore += score
	for seqid in seqsdic:
		if seqid == query:
			continue
		prob = float(scoredic[seqid]['score']) / sumscore
		scoredic[seqid]['prob'] = prob
	return scoredic

def compute_pairwise_file(alignfile, query):
	data = yaml_load_file()
	seqsdic = {}
	for seq in SeqIO.parse(alignfile, "fasta"):
		seqid = get_gi(seq.id) if (seq.id != query) else seq.id
		seqsdic[seqid] = seq.seq
	#print(">>", query)
	#print(">>", seqsdic)
	queryseq = setup_query(seqsdic[query])
	data[query]['hits'] = calculate_prob(seqsdic, query, queryseq)
	yaml_dump_file(data)

def transpose_file(infile, outfile):
	indata = open(infile, "r").read()
	#lines = '\n'.join(''.join(i) for i in zip(*indata.split()))
	lines = '\n'.join(''.join(i) for i in itertools.zip_longest(*indata.split(), fillvalue='#'))
	fh = open(outfile, "w")
	fh.write(lines)
	fh.close()

def randomize_file(infile, outfile, queryregion):
	lines = open(infile).readlines()
	fh = open(outfile, "w")
	for i in range(queryregion):
		r = random.randrange(queryregion)
		#print(r)
		fh.write(lines[r].rstrip("\r|\n") + "\n")
	fh.write(''.join(lines[queryregion:]))
	fh.close()

def get_seqsinfo(infile, query):
	#print(infile)
	a = open(infile, 'r').readlines()
	seqsdic = {}
	for line in a:
		line = line.rstrip("\r|\n")
		cols = line.split("XXXXX")
		#print("(" + cols[0] + ") (" + cols[1] + ") (" + query + ")")
		geneid = cols[1].rstrip("#")
		seqsdic[geneid] = cols[0]
	return seqsdic

def prob_highest(scoredic):
	tophit = []
	toppro = 0
	for hitid in scoredic:
		if toppro == 0:
			tophit.append(hitid)
			toppro = scoredic[hitid]['prob']
		elif scoredic[hitid]['prob'] == toppro:
			tophit.append(hitid)
		elif scoredic[hitid]['prob'] > toppro:
			tophit = [hitid]
			toppro = scoredic[hitid]['prob']
	return tophit

def update_bootstrap(confidence, qid):
	data = yaml_load_file()
	data[qid]['bootstrap'] = confidence
	yaml_dump_file(data)

def bootstrap_muscle_alignment_new(alignfile, query):
	print(alignfile)

def bootstrap_muscle_alignment(alignfile, query):
	seqsdic = {}
	for seq in SeqIO.parse(alignfile, "fasta"):
		seqid = get_gi(seq.id) if (seq.id != query) else seq.id
		seqsdic[seqid] = seq.seq

	queryseq, querystart, queryend = setup_query_pos(seqsdic[query])
	queryregion = queryend - querystart + 1

	bootstrap_confidence = {}

	for i in range(my_module.BOOTSTRAP):
		fh = open("tmp_"+str(i), "w")
		for seqid in seqsdic:
			fh.write(str(seqsdic[seqid][querystart:queryend+1]) + "XXXXX" + seqid.ljust(100) + "\n")
		fh.close()

		transpose_file("tmp_"+ str(i), "tmp_"+ str(i) + "_t")
		randomize_file("tmp_"+ str(i) + "_t", "tmp_"+ str(i) + "_t_s", queryregion)
		transpose_file("tmp_"+ str(i) + "_t_s", "tmp_"+ str(i) + "_t_s_t")

		seqsinfo = get_seqsinfo("tmp_"+ str(i) + "_t_s_t", query)
		#print(">>>", query)
		#print(">>>", seqsinfo)
		queryseq = setup_query(seqsinfo[query])
		s = calculate_prob(seqsinfo, query, queryseq)
		hitid = prob_highest(s)
		#print(hitid)
		for hid in hitid:
			if hid in bootstrap_confidence:
				bootstrap_confidence[hid] += (float(1)/len(hitid)) * (100/my_module.BOOTSTRAP)
			else:
				bootstrap_confidence[hid] = (float(1)/len(hitid)) * (100/my_module.BOOTSTRAP)
		os.remove("tmp_"+str(i))
		os.remove("tmp_"+str(i) + "_t")
		os.remove("tmp_"+str(i) + "_t_s")
		os.remove("tmp_"+str(i) + "_t_s_t")
		#break
	#print(bootstrap_confidence)
	update_bootstrap(bootstrap_confidence, query)

def compute():
	yamlfile = yaml_load_file()
	tot = len(yamlfile.keys())
	count = 1
	print("INFO: Running BOOTSTRAP")
	for seqid in yamlfile:
		#print(seqid)
		sys.stdout.write("Files: %d of %d   \r" % (count, tot))
		count += 1
		#compute_pairwise_file("multi_" + seqid + ".fasta.maln", seqid )
		bootstrap_muscle_alignment("multi_" + seqid + ".fasta.maln", seqid )
		#os.remove("multi_" + seqid + ".fasta.maln")
		sys.stdout.flush()
	print("DONE: BOOTSTRAP complete.")
