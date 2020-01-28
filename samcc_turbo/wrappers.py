# wrappers for external software

import subprocess

def run_dssp(pdbpath, dssp_binpath, outfile_path=None):
	"""Run DSSP software"""
	#FIXME expand

	if outfile_path:
		print(outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp')
		subprocess.call(dssp_binpath + ' -i ' + pdbpath + ' -o ' + outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp', shell=True)
		return outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp'
	else:
		print(pdbpath.split('/')[-1].split('.')[0] + '.dssp')
		subprocess.call(dssp_binpath + ' -i ' + pdbpath + ' -o ' + pdbpath.split('/')[-1].split('.')[0] + '.dssp', shell=True)
		return pdbpath.split('/')[-1].split('.')[0] + '.dssp'

def run_socket(pdbpath, dssppath, socket_binpath):
	"""Run Socket using pdb and calculated dssp file"""
	#FIXME expand

	subprocess.call(socket_binpath + ' -c 7.4 -l -f ' + pdbpath + ' -s ' + dssppath + ' > ' + pdbpath.split('.')[0] + '.socket', shell=True)

	return pdbpath.split('.')[0] + '.socket'
