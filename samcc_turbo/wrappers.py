### wrappers for external software ###

import subprocess
import tempfile

def run_dssp(pdbpath, dssp_binpath):
	"""Run DSSP software.
	Only 2.2.x versions are compatibile with Socket.

	Arguments:
	pdbpath      - path to PDB file
	dssp_binpath - path to dssp bin file
	
	Returns:
	Filepath to created .dssp file
	"""

	f = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
	f.close()

	outfile_name = f.name

	cmd = f"{dssp_binpath} -i {pdbpath} -o {outfile_name}"

	try:
		subprocess.run(cmd, shell=True, check=True)
		return outfile_name
	except subprocess.CalledProcessError:
		print('SamCC-Turbo requires dssp 2.2.1 to run. Specify path in bin_paths')

def run_socket(pdbpath, dssppath, socket_binpath, cut_off=7.4):
	"""Run Socket using pdb and calculated dssp file

	Arguments:
	pdbpath        - path to PDB file
	dssppath       - path to .dssp file
	socket_binpath - path to Socket bin file
	cut_off        - Socket cutoff parameter (max distance where knobes-into-holes are detected)

	Returns:
	Filepath to created Socket output file

	"""

	f = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
	f.close()

	outfile_name = f.name

	cmd = "%s -c %s -l -f %s -s %s > %s" % (socket_binpath, cut_off, pdbpath, dssppath, outfile_name)

	try:
		subprocess.run(cmd, shell=True, check=True)
		return outfile_name
	except subprocess.CalledProcessError:
		print('SamCC-Turbo requires Socket 3.03 to run. Specify path in bin_paths')
