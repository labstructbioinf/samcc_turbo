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


### Exceptions

class NoResidue(Exception):
    pass

class TooShort(Exception):
    pass


### DEV
from . import wrappers, socket_parser, socketClass
import importlib
importlib.reload(wrappers)
importlib.reload(socket_parser)
importlib.reload(socketClass)
### END-DEV

import pickle
from .wrappers import run_dssp
from .wrappers import run_socket
from .socket_parser import parse_socket_output
from .socketClass import socket_class


def run_samcc_turbo(pdbpath):
	"""Main function for running samcc-turbo."""

	###DEV-TXT: minimalne użycie: mierzaczka(pdb) odpala dssp, socket, parsuje socket, mierzy, zwraca df, wykres i sesje pymol
	###DEV-TXT: wersja 0.2: kontrola parametrów: czy wykres, df, pymol, ścieżki, zmienne coilowe (np. calc_crickdev)

	# get pdbid
	pdbid = pdbpath.split('/')[-1].split('.')[0]

	# run dssp
	dssppath = run_dssp(pdbpath)

	# run socket
	socket_binpath = '/home/kszczepaniak/Apps/socket3.03/socket'
	socketpath     = run_socket(pdbpath, dssppath, socket_binpath)
	socket_data    = parse_socket_output(socketpath)

	s = socket_class(socket_data, pdbpath)
	bundles = s.get_bundles(mode='ks', res_num_layer_detection=5)

	for b in enumerate(bundles):
		b[1].calc_bundleaxis()
		b[1].calc_periodicity()
		b[1].calc_radius()
		b[1].calc_crick()
		b[1].calc_crickdev(3.5, 7, optimal_ph1=19.5)
		b[1].calc_axialshift()
		b[1].plot(pdbid + '.png', elements=['Periodicity', 'Radius', 'CrickDev', 'Shift'])

		#print(b[1].gendf())

		pickle.dump(b[1].gendf(), open(pdbpath.split('.')[0] + '_coil_' + str(b[0]) + '.p', 'wb'))

		b[1].pymol_plot_layer(filename=pdbpath ,savepath='/'.join(pdbpath.split('/')[:-1]), suffix='coil_' + str(b[0]),
								pymol_version=2.0, color_selection=True, helix_order=b[1].helix_order)
