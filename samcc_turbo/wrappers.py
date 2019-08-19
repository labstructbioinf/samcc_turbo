# wrappers for external software

import subprocess

def run_dssp(pdbpath):
	"""Run DSSP software"""
	#FIXME expand

	subprocess.call('dssp -i ' + pdbpath + ' -o ' + pdbpath.split('.')[0] + '.dssp', shell=True)

	return pdbpath.split('.')[0] + '.dssp'

def run_socket(pdbpath, dssppath, socket_binpath):
	"""Run Socket using pdb and calculated dssp file"""
	#FIXME expand

	subprocess.call(socket_binpath + ' -c 7.4 -l -f ' + pdbpath + ' -s ' + dssppath + ' > ' + pdbpath.split('.')[0] + '.socket', shell=True)

	return pdbpath.split('.')[0] + '.socket'
