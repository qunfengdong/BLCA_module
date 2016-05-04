import os
import glob
from Bio import SeqIO
import importlib

my_module = ''

def config_file():
	path = os.getcwd()
	if not os.path.exists('config.py'):
		raise Exception("""Config file "config.py" is not present.

Please visit github (https://github.com/qunfengdong/BLCA)
Copy paste the contents of config file to config.py in the current folder (%s)""" % (path))

def config_fastafile():
	if my_module.FILENAME == '':
		raise Exception("""No input FASTA file provided in "config.py"

Please verify the contents of config.py
""")

def blast_binary():
	if my_module.BLAST_BINARY == '':
		raise Exception("""The directory for BLAST binary files are not provided in "config.py"

Please verify the contents of config.py
""")

def blast_database():
	if my_module.BLAST_DATABASE == '':
		raise Exception("""The database for BLAST is not provided in "config.py"

Please verify the contents of config.py
""")

def muscle_binary():
	if my_module.MUSCLE_BINARY == '':
		raise Exception("""The binary for MUSCLE is not provided in "config.py"

Please verify the contents of config.py
""")

def config_fastafile_exists():
	path = os.getcwd()
	if not os.path.exists(my_module.FILENAME):
		raise Exception("""FASTA file (%s) does not exist at (%s)""" % (my_module.FILENAME, path))

def config_fastafile_format():
	fasta_sequences = SeqIO.parse(open(my_module.FILENAME),'fasta')
	count = 0
	for fasta in fasta_sequences:
		count += 1
	if count == 0:
		raise Exception("""FASTA file (%s) is not in FASTA format""" % my_module.FILENAME)

def config_fastafile_empty():
	if os.stat(my_module.FILENAME).st_size == 0:
		raise Exception("""FASTA file (%s) is am empty file.""" % my_module.FILENAME)

def blast_binary_exists():
	if not os.path.exists(my_module.BLAST_BINARY):
		raise Exception("""BLAST binary files are not present at (%s)""" % my_module.BLAST_BINARY)

def blast_database_exists():
	f = glob.glob(my_module.BLAST_DATABASE + ".n*")
	#print(my_module.BLAST_DATABASE)
	#print(f)
	if len(f) == 0:
		raise Exception("""BLAST database files are not present at (%s)""" % my_module.BLAST_DATABASE)

def muscle_binary_exists():
	if not os.path.exists(my_module.MUSCLE_BINARY):
		raise Exception("""MUSCLE binary file is not present at (%s)""" % my_module.MUSCLE_BINARY)

def fasta_format():
	try:
		for seq in SeqIO.parse(my_module.FILENAME, "fasta"):
			pass
	except Exception as e:
		raise Exception(e)

def others_integers():
	if not isinstance(my_module.BLAST_CUTOFF_SCORE, int):
		raise Exception("""BLAST_CUTOFF_SCORE must be integer in config.py""")
	if not isinstance(my_module.BLAST_CUTOFF_PERCENT, int):
		raise Exception("""BLAST_CUTOFF_SCORE must be integer in config.py""")
	if (my_module.BLAST_CUTOFF_PERCENT > 100) or (my_module.BLAST_CUTOFF_PERCENT < 0):
		raise Exception("""BLAST_CUTOFF_PERCENT must be between 0-100 in config.py""")
	if not isinstance(my_module.BLAST_COVERAGE, int):
		raise Exception("""BLAST_COVERAGE must be integer in config.py""")
	if (my_module.BLAST_COVERAGE > 100) or (my_module.BLAST_COVERAGE < 0):
		raise Exception("""BLAST_COVERAGE must be between 0-100 in config.py""")
	if not isinstance(my_module.BLAST_PERCENTAGE_IDENTITY, int):
		raise Exception("""BLAST_PERCENTAGE_IDENTITY must be integer in config.py""")
	if (my_module.BLAST_PERCENTAGE_IDENTITY > 100) or (my_module.BLAST_PERCENTAGE_IDENTITY < 0):
		raise Exception("""BLAST_PERCENTAGE_IDENTITY must be between 0-100 in config.py""")
	if not isinstance(my_module.HIT_SEQUENCE_BPS, int):
		raise Exception("""HIT_SEQUENCE_BPS must be integer in config.py""")
	if not isinstance(my_module.ALIGNMENT_GAP, float) and not isinstance(my_module.ALIGNMENT_GAP, int):
		raise Exception("""ALIGNMENT_GAP must be integer in config.py""")
	if not isinstance(my_module.ALIGNMENT_MISMATCH, float) and not isinstance(my_module.ALIGNMENT_MISMATCH, int):
		raise Exception("""ALIGNMENT_MISMATCH must be integer in config.py""")
	if not isinstance(my_module.ALIGNMENT_MATCH, float) and not isinstance(my_module.ALIGNMENT_MATCH, int):
		raise Exception("""ALIGNMENT_MATCH must be integer in config.py""")
	if not isinstance(my_module.BOOTSTRAP, int):
		raise Exception("""BOOTSTRAP must be integer in config.py""")

def verify():
	global my_module
	config_file()
	my_module = importlib.import_module('config')
	config_fastafile()
	config_fastafile_exists()
	config_fastafile_empty()
	config_fastafile_format()
	#fasta_format()
	blast_binary()
	blast_binary_exists()
	blast_database()
	blast_database_exists()
	muscle_binary()
	muscle_binary_exists()
	others_integers()

