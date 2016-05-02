import re
import yaml
import importlib
from importlib.machinery import SourceFileLoader
try:
	my_module = importlib.import_module('config')
except:
	this_dir, this_filename = os.path.split(__file__)
	my_module = SourceFileLoader("settings", this_dir + "/settings.py").load_module()

def yaml_dump_file(data):
	stream = open(my_module.FILENAME + '.yaml', 'w')
	yaml.dump(data, stream)
	stream.close()

def yaml_load_file():
	stream = open(my_module.FILENAME + '.yaml')
	data = yaml.safe_load(stream)
	stream.close()
	return data

def yaml_load(filename):
	stream = open(filename)
	data = yaml.safe_load(stream)
	stream.close()
	return data

def get_gi(hit):
	return re.match(r"gi\|(.*)\|ref", hit).group(1)

def modify_id(id):
	return id.replace("|", "_")

