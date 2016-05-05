import blca

def info():
	return ('Bayesian-based LCA taxonomic classification method')

def execute():
	try:
		blca.verify()
		#blca.blast_seq()
		#blca.setup_msa()
		#blca.run_muscle()
		#blca.compute()
		#blca.annotate()
	except Exception as e:
		print("ERROR: ", e)