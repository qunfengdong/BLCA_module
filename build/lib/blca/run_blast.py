import os
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
#from config import *

def blast_seq():
	outfile = FILENAME + ".blastn"
	try:
		print("INFO: Running BLAST")
		blastn(cmd=BLAST_BINARY + "/blastn", query=FILENAME, db=BLAST_DATABASE, evalue=1e-20, out=outfile, show_gis="true", dust = 'no', soft_masking = "false", num_descriptions=500, num_alignments=500)()
		print("INFO: Verifying BLAST")
		if(verify_blast(outfile)):
			raise "BLAST output count is not equal to input sequence count"
		print("INFO: BLAST complete")
	except BaseException as e:
		print("Unable to BLAST")
		print(e)

def verify_blast(outputfile):
	count_out = os.popen("grep -c 'Query=' " + outputfile).read()
	count_in = os.popen("grep -c '>' " + FILENAME).read()
	if(int(count_in.rstrip("\r|\n")) !=  int(count_out.rstrip("\r|\n"))):
		return True
	return False

