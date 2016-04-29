from setuptools import setup

def readme():
	with open("README.rst") as f:
		return f.read()

setup(name='blca',
			version='0.6',
			description='Bayesian-based LCA taxonomic classification method',
			url='https://github.com/qunfengdong/BLCA',
			author='Kashi Reva',
			author_email='kashi.mail@gmail.com',
			license='MIT',
			packages=['blca'],
			include_package_data=True,
			package_data={},
			install_requires=[
				'biopython', 'pyyaml', 'yaml',
			],
			test_suite='nose.collector',
			tests_require=['nose'],
			zip_safe=False
			)

