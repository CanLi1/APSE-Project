#!/usr/bin/python
from __future__ import print_function
import argparse
import idaes_models.python_c_wrapper.autogenerate.wrapper
import os

def main():
#if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Python/C Wrapper",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--path', help='Directory of Python modules.',required=True)
	parser.add_argument('--dest', help='Directory to create C library.',required=True)
	
	parser.add_argument('--CPPFLAGS', help='-I[ADOL-C include, e.g: -I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL] \
-I[Python include, e.g: -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7] \
defaults to $CPPFLAGS environment variable.',default="")
	parser.add_argument('--MODELS_HOME', help='Home directory of models repo. \
	defaults to $MODELS_HOME environment variable.',default="")
	parser.add_argument('--ADOLC_LIBS', help='ADOL-C lib location, e.g: /usr/local/lib64., \
	defaults to $ADOLC_LIBS environment variable.',default="")
	
	args = parser.parse_args()

	# Check on args:
	if not args.CPPFLAGS or not args.MODELS_HOME or not args.ADOLC_LIBS:
		print ("Error: Make sure you either set the proper environment variabels for CPPFLAGS,MODELS_HOME,ADOLC_LIBS or pass them.")
		exit()

	# Compile wrapper library:
	compile_wrapper_lib_status = wrapper.compile_wrapper_library(args.CPPFLAGS,args.MODELS_HOME,args.ADOLC_LIBS)
	if not compile_wrapper_lib_status:
		print ("Wrapper library failed to build.")
		exit()
	else:
		# Compile properties library:
		funcgen_statements , addfunc_statements = wrapper.get_method_list(args.path)
		if (len(funcgen_statements) == 0 or len(addfunc_statements) == 0):
			print ("Error during parsing the Python directory {0}.".format(args.path))
			exit()	
		file_status = wrapper.generate_c_code(args.dest,funcgen_statements,addfunc_statements)
		if not file_status:
			print ("Error during generating C code. Do you have access to the directory {1}?".format(args.dest))
			exit()
		compile_properties_lib_status = wrapper.wrap(args.dest,args.CPPFLAGS,args.MODELS_HOME,args.ADOLC_LIBS)
		if (compile_properties_lib_status):
			print ("Properties library has been generated @ {0}/{1}.{2}".format(args.dest,wrapper.LIB_NAME,wrapper.SO_FILE))
			wrapper.generate_export_statements(args.dest,args.path,args.MODELS_HOME,args.ADOLC_LIBS,args.CPPFLAGS)
			print ("Set environment variables using 'source {0}/{1}'".format(args.dest,wrapper.HELPER_NAME))
		else:
			print ("Failed to generate property library.")
			exit()

