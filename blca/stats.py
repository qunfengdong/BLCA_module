import random
import timeit
import sys
import itertools
from .helpers import *
from Bio import SeqIO
from .helpers import my_module

def setup_query(seq):
	"""
	Add # to the ends of the query sequence
	:param seq:
	:return: seq
	"""
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
	"""
	Get the position of the start and end position of the alignment
	:param seq:
	:return:
	"""
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
	"""
	Get the score
	:param seq1:
	:param seq2:
	:return:
	"""
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
	"""
	Calculate probability
	:param seqsdic:
	:param query:
	:param queryseq:
	:return: dic
	"""
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
    if sumscore == 0:
      prob = 0
    else:
      prob = float(scoredic[seqid]['score']) / sumscore
		scoredic[seqid]['prob'] = prob
	#print(scoredic)
	return scoredic

def compute_pairwise_file(yamlfile, alignfile, query):
	"""
	Computes the score
	:param yamlfile:
	:param alignfile:
	:param query:
	:return:
	"""
	seqsdic = {}
	for seq in SeqIO.parse(alignfile, "fasta"):
		seqid = get_gi(seq.id) if (seq.id != query) else seq.id
		seqsdic[seqid] = seq.seq
	#print(">>", query)
	#print(">>", seqsdic)
	queryseq = setup_query(seqsdic[query])
	yamlfile[query]['hits'] = calculate_prob(seqsdic, query, queryseq)


def prob_highest(scoredic):
	"""
	Get the highest probability hit
	:param scoredic:
	:return: tophit
	"""
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

def bootstrap_muscle_alignment_new(yamlfile, alignfile, query):
	#print(alignfile)
	seqsdic = {}
	for seq in SeqIO.parse(alignfile, "fasta"):
		seqid = get_gi(seq.id) if (seq.id != query) else seq.id
		seqsdic[seqid] = seq.seq

	queryseq, querystart, queryend = setup_query_pos(seqsdic[query])
	queryregion = queryend - querystart + 1

	bootstrap_confidence = {}

	for i in range(my_module.BOOTSTRAP):
		bootdic = {}
		for seqid in seqsdic:
			bootdic[seqid] = seqsdic[seqid]
		seqsinfo = randomizefile(bootdic, query)
		queryseq = setup_query(seqsinfo[query])
		s = calculate_prob(seqsinfo, query, queryseq)
		hitid = prob_highest(s)
		#print(hitid)
		for hid in hitid:
			if hid in bootstrap_confidence:
				bootstrap_confidence[hid] += (float(1)/len(hitid)) * (100/my_module.BOOTSTRAP)
			else:
				bootstrap_confidence[hid] = (float(1)/len(hitid)) * (100/my_module.BOOTSTRAP)
	#print(bootstrap_confidence)
	yamlfile[query]['bootstrap'] = bootstrap_confidence

def randomizefile(bootdic, query):
	"""
	Get random sequence dic
	:param bootdic:
	:param query:
	:return:
	"""
	count = len(bootdic[query])
	randomindex = []
	randombootdic = {}
	for i in range(count):
		r = random.randrange(count)
		randomindex.append(r)

	for id in bootdic.keys():
		seq = []
		for i in randomindex:
			seq.append(bootdic[id][i])
		randombootdic[id] = ''.join(seq)
	return randombootdic

def compute():
	yamlfile = yaml_load_file()
	tot = len(yamlfile.keys())
	count = 1
	print("INFO: Running BOOTSTRAP")
	start_time = timeit.default_timer()
	for seqid in yamlfile:
		#print(seqid)
		sys.stdout.write("Files: %d of %d   \r" % (count, tot))
		count += 1
		compute_pairwise_file(yamlfile, "multi_" + seqid + ".fasta.maln", seqid )
		bootstrap_muscle_alignment_new(yamlfile, "multi_" + seqid + ".fasta.maln", seqid )
		os.remove("multi_" + seqid + ".fasta.maln")
		os.remove("multi_" + seqid + ".fasta")
		sys.stdout.flush()
	yaml_dump_file(yamlfile)
	elapsed = timeit.default_timer() - start_time
	print("DONE: BOOTSTRAP complete. (time taken: %d secs)" % int(elapsed))
