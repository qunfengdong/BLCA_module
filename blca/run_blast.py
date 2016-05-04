import os
import timeit
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from .helpers import my_module

def blast_seq():
	outfile = my_module.FILENAME + ".blastn"
	print("INFO: Running BLAST")
	start_time = timeit.default_timer()
	blastn(cmd=my_module.BLAST_BINARY + "/blastn", query=my_module.FILENAME, db=my_module.BLAST_DATABASE, evalue=1e-20, out=outfile, show_gis="true", dust = 'no', soft_masking = "false", num_descriptions=500, num_alignments=500)()
	print("INFO: Verifying BLAST")
	if(verify_blast(outfile)):
		raise Exception(""" (%s) BLAST output count is not equal to input sequence count""" % __file__ )
	elapsed = timeit.default_timer() - start_time
	print("DONE: BLAST complete (time taken: %d secs)" % int(elapsed))

def verify_blast(outputfile):
	count_out = os.popen("grep -c 'Query=' " + outputfile).read()
	count_in = os.popen("grep -c '>' " + my_module.FILENAME).read()
	if(int(count_in.rstrip("\r|\n")) !=  int(count_out.rstrip("\r|\n"))):
		return True
	return False

