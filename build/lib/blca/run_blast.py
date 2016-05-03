import os
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from .helpers import my_module

def blast_seq():
	outfile = my_module.FILENAME + ".blastn"
	try:
		print("INFO: Running BLAST")
		blastn(cmd=my_module.BLAST_BINARY + "/blastn", query=my_module.FILENAME, db=my_module.BLAST_DATABASE, evalue=1e-20, out=outfile, show_gis="true", dust = 'no', soft_masking = "false", num_descriptions=5, num_alignments=5)()
		print("INFO: Verifying BLAST")
		if(verify_blast(outfile)):
			raise "BLAST output count is not equal to input sequence count"
		print("DONE: BLAST complete")
	except BaseException as e:
		raise "Unable to BLAST"

def verify_blast(outputfile):
	count_out = os.popen("grep -c 'Query=' " + outputfile).read()
	count_in = os.popen("grep -c '>' " + my_module.FILENAME).read()
	if(int(count_in.rstrip("\r|\n")) !=  int(count_out.rstrip("\r|\n"))):
		return True
	return False

