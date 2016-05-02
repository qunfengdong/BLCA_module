import re
import yaml

def yaml_dump_file(filename, data):
	stream = open(filename + '.yaml', 'w')
	yaml.dump(data, stream)
	stream.close()

def yaml_load_file(filename):
	stream = open(filename + '.yaml')
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

