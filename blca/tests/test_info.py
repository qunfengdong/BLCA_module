from unittest import TestCase

import blca

class TestInfo(TestCase):
	def test_is_string(self):
		s = blca.info()
		self.assertTrue(isinstance(s, str))
