import re
from .annotate import *

def getindexvalue(cols, pos):
	try:
		return cols[pos]
	except IndexError:
		return ''

def parse_original(line):
	cols = line.split(";")
	#print(cols)
	taxa = {}
	taxa['superkingdom'] = getindexvalue(cols, 0)
	taxa['phylum'] = getindexvalue(cols, 1)
	taxa['class'] = getindexvalue(cols, 2)
	taxa['order'] = getindexvalue(cols, 3)
	taxa['family'] = getindexvalue(cols, 4)
	taxa['genus'] = getindexvalue(cols, 5)
	taxa['species'] = getindexvalue(cols, 6)
	#print(taxa)
	return taxa

def parse_rdp(cols):
	#print(cols)
	taxa = {}
	taxa['superkingdom'] = getindexvalue(cols, 11)
	taxa['phylum'] = getindexvalue(cols, 14)
	taxa['class'] = getindexvalue(cols, 17)
	taxa['order'] = getindexvalue(cols, 20)
	taxa['family'] = getindexvalue(cols, 23)
	taxa['genus'] = getindexvalue(cols, 26)
	#print(taxa)
	return taxa

def parse_megan5(line):
	cols = line.split("XXXXX")
	taxa = {}
	for t in cols:
		items = re.search('\[(.*)\] (.*); ([0-9]*)', t)
		name = items.group(1)
		taxa[name.lower()] = items.group(2)
		taxa[name.lower() + "_score"] = items.group(3)
	#print(taxa)
	return taxa

def parse_megan6(line):
	## Modify the output file
	'''
	:%s/; ;/\t/g
	:%s/;p/XXXp/g
	'''
	shand = {'p': 'phylum', 'o': 'order', 'c': 'class', 'f': 'family', 'g': 'genus', 's': 'species', 'd': 'superkingdom'}
	cols = line.split("XXX")
	taxa = {}
	for t in cols:
		items = re.search('(.*)__(.*); ([0-9]*)', t)
		name = items.group(1)
		taxa[shand[name]] = items.group(2)
		taxa[shand[name] + "_score"] = items.group(3)
	#print(taxa)
	return taxa


def parse_blca_name(col):
	rcol = col.replace(' ', '')
	if rcol == '':
		return '', ''
	a = re.search('(.*) \(([0-9]*)\)', col)
	return a.group(1), a.group(2)

def parse_blca(cols):
	taxa = {}
	taxa['superkingdom'], taxa['superkingdom_score'] = parse_blca_name(cols[8])
	taxa['phylum'], taxa['phylum_score'] = parse_blca_name(cols[7])
	taxa['class'], taxa['class_score'] = parse_blca_name(cols[6])
	taxa['order'], taxa['order_score'] = parse_blca_name(cols[5])
	taxa['family'], taxa['family_score'] = parse_blca_name(cols[4])
	taxa['genus'], taxa['genus_score'] = parse_blca_name(cols[3])
	taxa['species'], taxa['species_score'] = parse_blca_name(cols[2])
	taxa['org'], taxa['org_score'] = parse_blca_name(cols[1])
	#print(taxa)
	return(taxa)

def validate_print(taxa, orig, blca, megan, rdp, count_blca, count_megan, count_rdp):
	print("(", orig[taxa], ")") if taxa in orig else print('')
	print("(", blca[taxa], ")") if taxa in blca else print('')
	print("(", megan[taxa], ")") if taxa in megan else print('')
	print("(", rdp[taxa], ")") if taxa in rdp else print('')
	print(count_blca, count_megan, count_rdp)
	print("-----")

def validate_megan(filename):
	order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	for taxa in order:
		fh = open(filename)
		print(taxa)
		count_megan = 0
		count_blca = 0
		count_rdp = 0
		#order = ['genus']
		for line in fh:
			line = line.rstrip("\r|\n")
			#print(line)
			cols = line.split("\t")
			id = cols[0]
			#print("\n\n", id)
			orig = parse_original(cols[9])
			megan = parse_megan6(cols[10])
			blca = parse_blca(cols)

			#validate(taxa, orig, megan, blca)

			if taxa in megan:
				if megan[taxa] in orig[taxa]:
					count_megan += 1

			if taxa in blca:
				if blca[taxa] in orig[taxa]:
					count_blca += 1

		fh.close()
		print(count_megan, count_blca)
		#break

def validate(filename):
	order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	for taxa in order:
		fh = open(filename)
		print(taxa)
		count_megan = 0
		count_blca = 0
		count_rdp = 0
		#order = ['genus']
		for line in fh:
			line = line.rstrip("\r|\n")
			#print(line)
			cols = line.split("\t")
			id = cols[0]
			#print("\n\n", id)
			blca = parse_blca(cols)
			orig = parse_original(cols[9])
			megan = parse_megan6(cols[10])
			rdp = parse_rdp(cols)

			if taxa in megan:
				if megan[taxa] == orig[taxa]:
					count_megan += 1

			if taxa in blca:
				if blca[taxa] == orig[taxa]:
					count_blca += 1

			if taxa in rdp:
				if rdp[taxa] == orig[taxa]:
					count_rdp += 1

			#if taxa == 'species':
			#	validate_print(taxa, orig, blca, megan, rdp, count_blca, count_megan, count_rdp)

		fh.close()
		print(count_blca, count_megan, count_rdp)
		#break

