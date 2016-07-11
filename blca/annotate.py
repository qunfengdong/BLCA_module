import subprocess
import os
from .helpers import *
from .helpers import my_module

def read_yaml():
	this_dir, this_filename = os.path.split(__file__)
	datadir = os.path.abspath(os.path.join(this_dir, os.pardir))
	read_data = yaml_load_file()
	print("INFO: Loading taxonomy information.")
	#print("Loading: storage")
	read_taxid = yaml_load(datadir + '/blca/subset_gi_taxid.yaml')
	#print("Loading: gi_taxid")
	read_names = yaml_load(datadir + '/blca/subset_names.yaml')
	#print("Loading: names")
	read_nodes = yaml_load(datadir + '/blca/subset_nodes.yaml')
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


def annotate():
	'''
	Get the annotation
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

	outfile = open(my_module.OUTFILE, "w")
	#print("\t" + "\t".join(reversed(order)))
	outfile.write("sequence_id\t" + "\t".join(reversed(order)) + "\n")
	#ann = read_16s_annotation(read_names, read_nodes)
	#megan = read_megan_output(filename)
	#rdp = read_rdp_output(filename)

	for seq_id in diclin:
		#print(seq_id, end='')
		outfile.write(seq_id)
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
			#print("\t", end='')
			outfile.write("\t")

			if (taxa in con) and con[taxa]:
				#print(con[taxa])
				name = max(con[taxa], key=lambda i: con[taxa][i])
				if con[taxa][name] >= my_module.CUTOFF:
					#print(name + " (" + str(round(con[taxa][name])) + ")", end='')
					outfile.write(name + " (" + str(round(con[taxa][name])) + ")")
				else:
					outfile.write('NA')
			else:
				#print('', end='')
				outfile.write('NA')
		outfile.write("\n")
	outfile.close()
	print("DONE: Species annotation.")
	print("RESULTS: Annotations are available in the file " + my_module.OUTFILE)






