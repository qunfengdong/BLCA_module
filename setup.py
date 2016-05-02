from setuptools import setup

def readme():
	with open("README.md") as f:
		return f.read()

setup(name='blca',
			version='0.6',
			description='Bayesian-based LCA taxonomic classification method',
			url='https://github.com/qunfengdong/BLCA',
			author='Kashi Reva',
			author_email='kashi.mail@gmail.com',
			license='MIT',
			packages=['blca'],
			data_files=[('blca', ['data/subset_gi_taxid.yaml', 'data/subset_names.yaml', 'data/subset_nodes.yaml'])],
			#include_package_data=True,
			package_data={'blca': ['data/*.yaml']},
			scripts=['bin/runblca.py'],
			install_requires=[
				'biopython', 'pyyaml',
			],
			test_suite='nose.collector',
			tests_require=['nose'],
			zip_safe=False
			)

