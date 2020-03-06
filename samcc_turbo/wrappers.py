# wrappers for external software

import subprocess, tempfile

def run_dssp(pdbpath, dssp_binpath, outfile_path=None):
	"""Run DSSP software"""
	#FIXME expand

	if outfile_path:
		#print(outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp')
		subprocess.call(dssp_binpath + ' -i ' + pdbpath + ' -o ' + outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp', shell=True)
		return outfile_path + pdbpath.split('/')[-1].split('.')[0] + '.dssp'
	else:
		#print(pdbpath.split('/')[-1].split('.')[0] + '.dssp')
		subprocess.call(dssp_binpath + ' -i ' + pdbpath + ' -o ' + pdbpath.split('/')[-1].split('.')[0] + '.dssp', shell=True)
		return pdbpath.split('/')[-1].split('.')[0] + '.dssp'

def run_socket(pdbpath, dssppath, socket_binpath, cut_off=7.4):
	"""Run Socket using pdb and calculated dssp file"""
	#FIXME expand
	
	f = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
	f.close()
	
	outfile_name = f.name
	
	cmd = "%s -c %s -l -f %s -s %s > %s" % (socket_binpath, cut_off, pdbpath, dssppath, outfile_name)

	# socket_binpath + ' -c 7.4 -l -f ' + pdbpath + ' -s ' + dssppath + ' > ' + pdbpath.split('.')[0] + '.socket'
	subprocess.call(cmd, shell=True)
	
	return outfile_name

	#return pdbpath.split('.')[0] + '.socket'
