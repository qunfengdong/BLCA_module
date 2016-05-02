import blca

def info():
	return ('Bayesian-based LCA taxonomic classification method')

def execute():
	blca.blast_seq()
	blca.setup_msa()
	blca.run_muscle()
	blca.compute()
	blca.annotate()