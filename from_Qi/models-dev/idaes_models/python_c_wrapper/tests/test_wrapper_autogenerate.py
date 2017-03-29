# from context import setup_path; setup_path()
from idaes_models.python_c_wrapper.autogenerate import wrapper

import unittest
from unittest import TestCase
from idaes_models.core.util.misc import category

"""
Before running this test, set environment variables:
export CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7"
export MODELS_HOME=/home/jovyan/models/models-fork/models
export ADOLC_LIBS=/usr/local/lib64
"""
SKIP_ALL_TESTS = False

class TestAutogenerate(TestCase):

	def setUp(self):
		global SKIP_ALL_TESTS
		try:
			import os
			if 'MODELS_HOME' in os.environ:
				self._models_home = os.environ['MODELS_HOME']
			else:
				self._models_home = "/home/jovyan/models"
			if 'CPPFLAGS' in os.environ:
				self._cpp_flags = os.environ['CPPFLAGS']
			else:
				self._cpp_flags = "-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7"
			if 'ADOLC_LIBS' in os.environ:
				self._adolc_libs = os.environ['ADOLC_LIBS']
			else:
				self._adolc_libs = "/usr/local/src/pyadolc/PACKAGES/ADOL-C/ADOL-C/.libs"
			SKIP_ALL_TESTS = False
		except Exception:
			print ("Error in getting environment variables. Set MODELS_HOME,CPPFLAGS and ADOLC_LIBS variables.")
			print ("Skipping tests.")
			SKIP_ALL_TESTS = True

	@category('frequent', 'docker-only')
	def test_wrapper_compilation_failure(self):
		# Bad paths :
		global SKIP_ALL_TESTS
		if SKIP_ALL_TESTS:
			self.skipTest('Skipping.')
		CPPFLAGS = ADOLC_LIBS = MODELS_HOME = "/not/a/real/path"
		r_val = wrapper.compile_wrapper_library(CPPFLAGS,MODELS_HOME,ADOLC_LIBS)
		assert (r_val == False)

	@category('frequent', 'docker-only')
	def test_wrapper_compilation_success(self):
		global SKIP_ALL_TESTS
		if SKIP_ALL_TESTS:
			self.skipTest('Skipping.')

		CPPFLAGS , ADOLC_LIBS , MODELS_HOME = self._cpp_flags , self._adolc_libs , self._models_home

		import os.path
		import subprocess
		r_val = wrapper.compile_wrapper_library(CPPFLAGS,MODELS_HOME,ADOLC_LIBS)
		assert (r_val == True)
		# Commented temporarily.
		#assert(os.path.isfile("{0}/{1}".format(MODELS_HOME,"src/wrapper/.libs/libwrapper.so")))

	@category('frequent', 'docker-only')
	def test_get_method_list(self):
		global SKIP_ALL_TESTS
		if SKIP_ALL_TESTS:
			self.skipTest('Skipping.')
		import os
		import sys
		module_dir = os.path.abspath(self._models_home + "/idaes_models/python_c_wrapper/sampledir")
		sys.path.insert(0, module_dir)
		funcgen,addfunc = wrapper.get_method_list(module_dir)

		# Only 4 methods are defined that fit the Python restrictions out of possible 6 (2 don't match criteria):
		assert((len(funcgen) == 4) and (len(funcgen) == len(addfunc)))

	@category('frequent', 'docker-only')
	def test_generate_c_code (self):
		global SKIP_ALL_TESTS
		if SKIP_ALL_TESTS:
			self.skipTest('Skipping.')
		# dest, funcgen_statements,addfunc_statements
		dest = "{0}/prop_lib/".format(self._models_home + "/idaes_models/python_c_wrapper")
		import os
		module_dir = os.path.abspath(self._models_home + "/idaes_models/python_c_wrapper/sampledir")
		funcgen,addfunc = wrapper.get_method_list(module_dir)

		assert(len(funcgen) > 0 and len(addfunc) > 0)
		dest_dir = os.path.abspath(self._models_home + "/idaes_models/python_c_wrapper/prop_lib")
		return_val = wrapper.generate_c_code(dest_dir,funcgen,addfunc)
		print (return_val)

if __name__ == '__main__':
    unittest.main()

