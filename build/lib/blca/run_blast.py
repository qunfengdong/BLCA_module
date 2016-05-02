import os
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
import importlib
from importlib.machinery import SourceFileLoader
try:
	my_module = importlib.import_module('config')
except:
	this_dir, this_filename = os.path.split(__file__)
	my_module = SourceFileLoader("settings", this_dir + "/settings.py").load_module()

def blast_seq():
	outfile = my_module.FILENAME + ".blastn"
	try:
		print("INFO: Running BLAST")
		blastn(cmd=my_module.BLAST_BINARY + "/blastn", query=my_module.FILENAME, db=my_module.BLAST_DATABASE, evalue=1e-20, out=outfile, show_gis="true", dust = 'no', soft_masking = "false", num_descriptions=500, num_alignments=500)()
		print("INFO: Verifying BLAST")
		if(verify_blast(outfile)):
			raise "BLAST output count is not equal to input sequence count"
		print("INFO: BLAST complete")
	except BaseException as e:
		print("Unable to BLAST")
		print(e)

def verify_blast(outputfile):
	count_out = os.popen("grep -c 'Query=' " + outputfile).read()
	count_in = os.popen("grep -c '>' " + my_module.FILENAME).read()
	if(int(count_in.rstrip("\r|\n")) !=  int(count_out.rstrip("\r|\n"))):
		return True
	return False

