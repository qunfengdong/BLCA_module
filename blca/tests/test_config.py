from unittest import TestCase

import blca
'''
def foo():
	raise Exception("doh")
'''
class TestInfo(TestCase):
	def test_no_config_file(self):
		#self.assertRaises(Exception, foo)
		try:
			blca.verify()
		except Exception as e:
			self.assertTrue("ERROR" in str(e))
