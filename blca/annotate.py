import subprocess
from .helpers import *

def read_yaml(filename):
	read_data = yaml_load_file(filename)
	#print("Loading: storage")
	read_taxid = yaml_load('data/subset/subset_gi_taxid.yaml')
	#print("Loading: gi_taxid")
	read_names = yaml_load('data/subset/subset_names.yaml')
	#print("Loading: names")
	read_nodes = yaml_load('data/subset/subset_nodes.yaml')
	#print("Loading: nodes")
	#print("=========================")
	return read_data, read_taxid, read_names, read_nodes

def get_path(tid, score, read_taxid, read_names, read_nodes):
	lin = {}
	#print(tid)
	if tid == '1161941':
		return lin
	lin['organism'] = read_names[tid]
	lin['organism_score'] = score
	for i in range(10):
		sid = read_nodes[tid]
		#print(sid, read_names[tid])
		# #print(sid['parent'])
		if(sid['rank'] != 'no rank'):
			rank = sid['rank']
			scorekey = rank + '_score'
			lin[rank] = read_names[tid]
			lin[scorekey] = score
		tid = sid['parent']
		if tid == 1:
			break
		#print(lin)
	return lin

def get_path_from_gi(gi, score, read_taxid, read_names, read_nodes):
	#print(read_taxid)
	#print(read_names)
	gi_taxid = read_taxid[gi]
	path = get_path(gi_taxid, score, read_taxid, read_names, read_nodes)
	return path

def get_original_taxonomy(id):
	ok = subprocess.check_output(["grep", id, "SILVA_123_SSURef_100_V34.taxonomy"], shell=True)
	ok.rstrip("\r|\n")
	print(ok.split("\t")[1])

def read_16s_annotation(read_names, read_nodes):
	gilist = {}
	order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	for gi in read_nodes.keys():
		lin = {}
		tid = gi
		for i in range(10):
			sid = read_nodes[tid]
			#print(sid, read_names[tid])
			# #print(sid['parent'])
			if(sid['rank'] != 'no rank'):
				rank = sid['rank']
				scorekey = rank + '_score'
				lin[rank] = read_names[tid]
			tid = sid['parent']
			if tid == 1:
				break
		s = []
		for taxa in reversed(order):
			if taxa in lin:
				s.append(lin[taxa])
			else:
				s.append("XXXXX")
		gilist[gi] = ";".join(s)
		#print(lin)
	#print(">>", gilist)
	return gilist

def read_original_annotation():
	#fh = open("SILVA_123_SSURef_100_V34.taxonomy", "r")
	fh = open("SILVA_123_SSURef_151_V4.taxonomy", "r")
	ann = {}
	for line in fh:
		line = line.rstrip("\r|\n")
		cols = line.split("\t")
		ann[cols[0]] = cols[1]
	return ann

def read_megan_output(filename):
	'''
	[kreva@kadesktop:tools]$ ./blast2lca \
	-top 20 \
	-i /home/kreva/Desktop/Workspace/2016_04_04_MEGAN_upgrade/analysis_001/16SMicrobial_150.V13.subsample2.fasta.blastn \
	-o /home/kreva/Desktop/Workspace/2016_04_04_MEGAN_upgrade/analysis_001/16SMicrobial_150.V13.subsample2.fasta.blastn.megan6.taxonomy.txt \
	-ko /home/kreva/Desktop/Workspace/2016_04_04_MEGAN_upgrade/analysis_001/16SMicrobial_150.V13.subsample2.fasta.blastn.kegg.txt
	'''
	fh = open(filename + ".blastn.megan6.taxonomy.txt", "r")
	#fh = open('16SMicrobial_150.V13.subsample1.fasta.dft.blastn.0.1.taxonomy', 'r')
	ann = {}
	for line in fh:
		line = line.rstrip("\r|\n")
		cols = line.split("\t")
		ann[cols[0]] = cols[1]
	return ann

def read_rdp_output(filename):
	fh = open(filename + ".rdp", "r")
	ann = {}
	for line in fh:
		line = line.rstrip("\r|\n")
		cols = line.split("\t\t")
		ann[cols[0]] = cols[1]
	return ann

def annotate(filename):
	'''
	latest
	:param filename:
	:return:
	'''
	diclin = {}
	read_data, read_taxid, read_names, read_nodes = read_yaml()
	for seq_id in read_data:
		diclin[seq_id] = {}
		for gi in read_data[seq_id]['bootstrap'].keys():
			#print(seq_id, gi, read_data[seq_id]['bootstrap'][gi])
			diclin[seq_id][gi] = get_path_from_gi(gi, read_data[seq_id]['bootstrap'][gi], read_taxid, read_names, read_nodes)

	order = ['organism', 'subspecies', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	print("\t" + "\t".join(reversed(order)))
	#ann = read_16s_annotation(read_names, read_nodes)
	#megan = read_megan_output(filename)
	#rdp = read_rdp_output(filename)

	for seq_id in diclin:
		print(seq_id, end='')
		con = {}
		for gi in diclin[seq_id]:
			for taxa in order:
				if taxa not in con:
					con[taxa] = {}
				if taxa not in diclin[seq_id][gi]:
					continue
				t = diclin[seq_id][gi][taxa]
				s = diclin[seq_id][gi][taxa + '_score']
				#print(t,s)
				if t in con[taxa]:
					#con[taxa][t]['score'] += s
					con[taxa][t] += float(s)
				else:
					#con[taxa][t] = {}
					con[taxa][t] = float(s)
					#con[taxa][t]['score'] = s

		'''
		for taxa in order:
			print(taxa)
			for name in con[taxa]:
				print("\t", name, "\t", round(con[taxa][name]['score'], 1))
		print("---")
		'''
		for taxa in reversed(order):
			print("\t", end='')

			if con[taxa]:
				name = max(con[taxa], key=lambda i: con[taxa][i])
				if con[taxa][name] >= 0:
					print(name + " (" + str(round(con[taxa][name])) + ")", end='')
			else:
				print('', end='')
		print("")


def annotate_megan_blca_rdp(filename):
	diclin = {}
	read_data, read_taxid, read_names, read_nodes = read_yaml(filename)
	for seq_id in read_data:
		diclin[seq_id] = {}
		for gi in read_data[seq_id]['bootstrap'].keys():
			#print(seq_id, gi, read_data[seq_id]['bootstrap'][gi])
			diclin[seq_id][gi] = get_path_from_gi(gi, read_data[seq_id]['bootstrap'][gi], read_taxid, read_names, read_nodes)

	order = ['org', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	#print("\t\t\t" + "\t".join(order))
	ann = read_16s_annotation(read_names, read_nodes)
	megan = read_megan_output(filename)
	rdp = read_rdp_output(filename)

	for seq_id in diclin:
		print(seq_id, end='')
		con = {}
		for gi in diclin[seq_id]:
			for taxa in order:
				if taxa not in con:
					con[taxa] = {}
				if taxa not in diclin[seq_id][gi]:
					continue
				t = diclin[seq_id][gi][taxa]
				s = diclin[seq_id][gi][taxa + '_score']
				#print(t,s)
				if t in con[taxa]:
					#con[taxa][t]['score'] += s
					con[taxa][t] += float(s)
				else:
					#con[taxa][t] = {}
					con[taxa][t] = float(s)
					#con[taxa][t]['score'] = s

		'''
		for taxa in order:
			print(taxa)
			for name in con[taxa]:
				print("\t", name, "\t", round(con[taxa][name]['score'], 1))
		print("---")
		'''
		for taxa in order:
			print("\t", end='')

			if con[taxa]:
				name = max(con[taxa], key=lambda i: con[taxa][i])
				if con[taxa][name] >= 80:
					print(name + " (" + str(round(con[taxa][name])) + ")", end='')
			else:
				print('', end='')

		print("\t", end='')
		giid = re.match(r"gi_([0-9]*)_|ref", seq_id).group(1)
		giidtaxid = read_taxid[giid]
		print(ann[giidtaxid], end='')

		print("\t", end='')
		print(megan[seq_id], end='')

		print("\t", end='')
		print(rdp[seq_id], end='')

		print("")


def annotate_megan_blca():
	diclin = {}
	read_data, read_taxid, read_names, read_nodes = read_yaml()
	for seq_id in read_data:
		diclin[seq_id] = {}
		for gi in read_data[seq_id]['bootstrap'].keys():
			#print(seq_id, gi, read_data[seq_id]['bootstrap'][gi])
			diclin[seq_id][gi] = get_path_from_gi(gi, read_data[seq_id]['bootstrap'][gi], read_taxid, read_names, read_nodes)

	order = ['org', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	print("\t\t\t" + "\t".join(order))
	ann = read_original_annotation()
	megan = read_megan_output()

	for seq_id in diclin:
		print(seq_id, end='')
		print("\t", end='')
		print(ann[seq_id], end='')
		print("\t", end='')
		print(megan[seq_id], end='')
		con = {}
		for gi in diclin[seq_id]:
			for taxa in order:
				if taxa not in con:
					con[taxa] = {}
				if taxa not in diclin[seq_id][gi]:
					continue
				t = diclin[seq_id][gi][taxa]
				s = diclin[seq_id][gi][taxa + '_score']
				#print(t,s)
				if t in con[taxa]:
					#con[taxa][t]['score'] += s
					con[taxa][t] += float(s)
				else:
					#con[taxa][t] = {}
					con[taxa][t] = float(s)
					#con[taxa][t]['score'] = s

		'''
		for taxa in order:
			print(taxa)
			for name in con[taxa]:
				print("\t", name, "\t", round(con[taxa][name]['score'], 1))
		print("---")
		'''
		for taxa in order:
			print("\t", end='')
			'''
			for name in con[taxa]:
				if con[taxa][name] >= 80:
					#print(name + " (" + str(round(con[taxa][name]['score'])) + ")", end='')
					print(name + " (" + str(round(con[taxa][name])) + ")", end='')
			'''

			if con[taxa]:
				name = max(con[taxa], key=lambda i: con[taxa][i])
				print(name + " (" + str(round(con[taxa][name])) + ")", end='')
			else:
				print('', end='')

		print("")


def annotate_blca():
	diclin = {}
	read_data, read_taxid, read_names, read_nodes = read_yaml()
	for seq_id in read_data:
		diclin[seq_id] = {}
		for gi in read_data[seq_id]['bootstrap'].keys():
			#print(seq_id, gi, read_data[seq_id]['bootstrap'][gi])
			diclin[seq_id][gi] = get_path_from_gi(gi, read_data[seq_id]['bootstrap'][gi], read_taxid, read_names, read_nodes)

	order = ['org', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	print("\t\t\t" + "\t".join(order))
	ann = read_original_annotation()
	megan = read_megan_output()

	for seq_id in diclin:
		print(seq_id, end='')
		con = {}
		for gi in diclin[seq_id]:
			for taxa in order:
				if taxa not in con:
					con[taxa] = {}
				if taxa not in diclin[seq_id][gi]:
					continue
				t = diclin[seq_id][gi][taxa]
				s = diclin[seq_id][gi][taxa + '_score']
				#print(t,s)
				if t in con[taxa]:
					#con[taxa][t]['score'] += s
					con[taxa][t] += float(s)
				else:
					#con[taxa][t] = {}
					con[taxa][t] = float(s)
					#con[taxa][t]['score'] = s

		for taxa in order:
			print("\t", end='')
			if con[taxa]:
				name = max(con[taxa], key=lambda i: con[taxa][i])
				print(name + " (" + str(round(con[taxa][name])) + ")", end='')
			else:
				print('', end='')

		print("")




