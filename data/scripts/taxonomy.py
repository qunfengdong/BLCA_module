import pickle
import glob

def getlistoffiles():
	a = glob.glob("gi_taxid*")
	a.remove("gi_taxid_nucl.dmp")
	return a

def yaml_dump_file(filename, data):
	stream = open(filename, 'w')
	yaml.dump(data, stream)
	stream.close()

def yaml_load_file(filename):
	stream = open(filename)
	data = yaml.safe_load(stream)
	stream.close()
	return data

def get_gi_numbers(filename):
	gitaxid = []
	with open(filename, 'r') as fh:
		for gi in fh:
			gi = gi.rstrip("\r|\n")
			gitaxid.append(gi)
	return gitaxid

def get_taxids(gilist):
	taxids = {}
	print(len(gilist))
	for ff in getlistoffiles():
		print(ff, len(gilist))
		with open(ff, "r") as fh:
			for line in fh:
				line = line.rstrip("\r|\n")
				gi, taxid = line.split("\t")
				if gi in gilist:
					taxids[gi] = taxid
					gilist.remove(gi)
	print(len(gilist))
	return taxids

def lineage(arrids):
  lindic = {}
  #while(arrids != ['1161941', '1']):
  while(len(arrids) > 2):
  #for i in range(7):
    print(arrids)
    fh = open("nodes.dmp", "r")
    for line in fh:
      line = line.rstrip("\r|\n")
      cols = line.split("\t|\t")
      #print(cols[0], cols[1])
      if(cols[0] in arrids):
        lindic[cols[0]] = {}
        lindic[cols[0]]['parent'] = cols[1]
        lindic[cols[0]]['rank'] = cols[2]
        arrids.remove(cols[0])
        if(cols[1] not in arrids):
          arrids.append(cols[1])
    fh.close()
    #break
  return lindic

def taxnames(arrids):
  print(arrids)
  names = {}
  while(len(arrids) > 0):
  #for i in range(3):
    fh = open("names.dmp", "r")
    for line in fh:
      line = line.rstrip("\r|\n")
      cols = line.split("\t|\t")
      #print("("+cols[0]+")")
      if((cols[0] in arrids) and ('scientific' in cols[3])):
        names[cols[0]] = cols[1]
        arrids.remove(cols[0])
    fh.close()
    print(arrids)
  #print(names)
  return names

def parse_gi_taxid():
	gilist = get_gi_numbers('16SMicrobialseq.fasta.names.gi')
	taxids = get_taxids(gilist)
	pickle.out(taxids, open('subset_gi_taxid.yaml', "wb"))

def parse_nodes():
  gis = yaml_load_file('subset_gi_taxid.yaml')
  list_gis = list(gis.values())
  lin = lineage(list_gis)
  yaml_dump_file('subset_nodes.yaml', lin)

def parse_names():
  lis = yaml_load_file('subset_nodes.yaml')
  list_gis = list(lis.keys())
  lin = taxnames(list_gis)
  yaml_dump_file('subset_names.yaml', lin)

if __name__ == "__main__":
  parse_gi_taxid()
  parse_nodes()
  parse_names()




