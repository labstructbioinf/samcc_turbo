### probably not needed imports
# import sys, itertools, glob
#
# from math import factorial
#
# from Bio.PDB import StructureBuilder
# from Bio.PDB import PDBIO
# from Bio.SeqUtils import seq1
#
# from Bio.PDB import Vector

# import warnings
#import importlib ### DEV
#importlib.reload(ld) ### DEV
#importlib.reload(pdl) ### DEV
# from io import StringIO
### end old imports




### DEV
from . import wrappers, socket_parser, socketClass, helper_functions
import importlib
importlib.reload(wrappers)
importlib.reload(socket_parser)
importlib.reload(socketClass)
importlib.reload(helper_functions)

### END-DEV

import pickle
import sys
from .wrappers import run_dssp
from .wrappers import run_socket
from .socket_parser import parse_socket_output
from .socketClass import socket_class


def run_samcc_turbo(pdbpath, mode='auto-detect', deffile=None, plot=True,
					save_df=True, bin_paths={'dssp':'dssp', 'socket':'socket'}):
	"""Main function for running samcc-turbo.

	Arguments:
	pdbpath   -- path to pdb file
	mode      -- mode of running SamCC (default 'auto-detect')
	           - auto-detect: automatic detection of layers in the bundle (requires
	                          installed dssp and Socket)
	           - deffile: use layer definition from file
	deffile   -- path to file defining layers used in deffile mode (default False)
	plot      -- plot results and also save plot to file (default True)
	save_df   -- save result DataFrame to pickle (default True)
	bin_paths -- dictionary of paths to binaries of dssp and socket
	             (default {'dssp': 'dssp', 'socket':'socket'})
	"""

	###DEV-TXT: minimalne użycie: mierzaczka(pdb) odpala dssp, socket, parsuje socket, mierzy, zwraca df, wykres i sesje pymol
	###DEV-TXT: wersja 0.1A: użycie pliku definiującego helisy (stare samcc)
	###DEV-TXT: wersja 0.2: kontrola parametrów: czy wykres, df, pymol, ścieżki

	# get pdbid
	pdbid = pdbpath.split('/')[-1].split('.')[0]

	# run dssp
	dssppath = run_dssp(pdbpath, bin_paths['dssp'])

	if mode == 'auto-detect':
		# run socket ('auto-detect' mode)
		socketpath     = run_socket(pdbpath, dssppath, bin_paths['socket'])
		socket_data    = parse_socket_output(socketpath)
		s 			   = socket_class(socket_data, pdbpath)

	elif mode == 'deffile':
		if deffile == None:
			print('Please provide file with bundle definition for measurement.')
			sys.exit(-1)
		s = socket_class(deffile, pdbpath)

	else:
		#FIXME add behaviour for unknown mode
		print('Unknown mode')

	bundles = s.get_bundles(mode=mode, res_num_layer_detection=5)

	for bid, bundle in enumerate(bundles):
		bundle.calc_bundleaxis()
		bundle.calc_periodicity()
		bundle.calc_radius()
		bundle.calc_crick()
		bundle.calc_crickdev(3.5, 7, optimal_ph1=19.5)
		bundle.calc_axialshift()
		if plot: # make plot and save it to file
			bundle.plot(pdbid + '.png', elements=['Periodicity', 'Radius', 'CrickDev', 'Shift'])

		if save_df: # dump pickle with dataframe of measured values
			pickle.dump(bundle.gendf(), open(pdbpath.split('.')[0] + '_coil_' + str(bid) + '.p', 'wb'))

		bundle.pymol_plot_layer(filename=pdbpath ,savepath='/'.join(pdbpath.split('/')[:-1]), suffix='coil_' + str(bid),
								pymol_version=2.0, color_selection=True, helix_order=bundle.helix_order)
