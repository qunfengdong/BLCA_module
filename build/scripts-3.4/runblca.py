#!/usr/bin/python3
import os
import importlib
from importlib.machinery import SourceFileLoader
try:
	my_module = importlib.import_module('config')
except:
	this_dir, this_filename = os.path.split(__file__)
	my_module = SourceFileLoader("settings", this_dir + "/settings.py").load_module()

print('kashi')
