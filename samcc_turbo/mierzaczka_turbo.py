### Main SamCC-Turbo file ###
DEBUG=False

import pickle
import sys
from .wrappers import run_dssp
from .wrappers import run_socket
from .socket_parser import parse_socket_output
from .socketClass import socket_class
from .bundleClass import bundleClass

def run_samcc_turbo(pdbpath, mode='auto-detect', defdata=None,
					plot=True, save_df=True, save_pse=True,
					bin_paths={'dssp':'dssp', 'socket':'socket'},
					layer_detect_n=5, max_dist='auto', search_set_n=9):
	"""Main function for running samcc-turbo.

	Arguments (general):
	pdbpath   -- path to pdb file (.pdb)
	mode      -- mode of running SamCC (default 'auto-detect')
	           - auto-detect: automatic detection of layers in the bundle (requires
	                          installed dssp and Socket)
			   - defdata: use layer definition from list
	plot      -- plot results and also save plot to file (default True)
	save_df   -- save result DataFrame to pickle (default True)
	save_pse  -- save pymol session with drawn layers (default True)
	bin_paths -- dictionary of paths to binaries of dssp and socket
	             (default {'dssp': 'dssp', 'socket':'socket'})

	Arguments (auto-detect mode specific):
	layer_detect_n -- number of residues from each helix that will be combined
					  into layers that will be considered in layer detection;
					  (default 5) [this is "r" parameter in the paper]
	max_dist       -- distance threshold [A] that will be used to indicate unusual
					  layer setting; if total distance between residues in layer
					  if greater that this number the layer will be discarded, if
					  no layer will pass this filter then layer detection will
					  be repeated with bigger layer_detect_n; if set to auto
					  then threshold will be assigned basing on oligomerization
					  (default auto)
	search_set_n   -- number of starting layer settings that will be selected
					  basing on minimal distance between residues; they will be compared
					  in terms of angle between all layers in structure (default 9)
	Arguments (defdata mode specific):
	defdata   -- list of parameters defining layer setting (default None)
	"""

	# get pdbid
	pdbid = pdbpath.split('/')[-1].split('.')[0]

	if mode == 'auto-detect':
		# run dssp and socket ('auto-detect' mode)
		dssppath       = run_dssp(pdbpath, bin_paths['dssp'])
		socketpath     = run_socket(pdbpath, dssppath, bin_paths['socket'])
		socket_data    = parse_socket_output(socketpath)
		s 			   = socket_class(socket_data, pdbpath)

		bundles = s.get_bundles(res_num_layer_detection=layer_detect_n,
								distance_threshold=max_dist,
								search_layer_setting_num=search_set_n)

		for bid, bundle in enumerate(bundles):
			bundle.calc_bundleaxis()
			bundle.get_helicesaxis()
			bundle.calc_periodicity()
			bundle.calc_radius()
			bundle.calc_crick()
			bundle.calc_crickdev(3.5, 7, optimal_ph1=19.5)
			bundle.calc_axialshift()
			bundle.assign_positions()
			if plot: # make plot and save it to file
				bundle.plot(pdbid + '.png', elements=['Periodicity', 'Radius', 'CrickDev', 'Shift'])

			if save_df: # dump pickle with dataframe of measured values
				pickle.dump(bundle.gendf(), open(pdbpath.split('/')[-1].split('.')[0] + '_coil_' + str(bid) + '.p', 'wb'))

			if save_pse:
				bundle.pymol_plot_layer(filename=pdbpath ,savepath='/'.join(pdbpath.split('/')[:-1]), suffix='coil_' + str(bid),
										pymol_version=2.0, color_selection=True, helix_order=bundle.helix_order, helices_axis=bundle.helices_axis)

	elif mode == 'defdata':

		bundle = bundleClass()
		bundle.from_defdata(pdbpath, *defdata)

		bundle.calc_bundleaxis()
		bundle.get_helicesaxis()
		bundle.calc_periodicity()
		bundle.calc_radius()
		bundle.calc_crick()
		bundle.calc_crickdev(3.5, 7, optimal_ph1=19.5)
		bundle.calc_axialshift()
		bundle.assign_positions()

		if plot: # make plot and save it to file
			bundle.plot(pdbid + '.png', elements=['Periodicity', 'Radius', 'CrickDev', 'Shift'])

		if save_df: # dump pickle with dataframe of measured values
			pickle.dump(bundle.gendf(), open(pdbpath.split('/')[-1].split('.')[0] + '_coil.p', 'wb'))

		if save_pse:
			bundle.pymol_plot_layer(filename=pdbpath ,savepath='/'.join(pdbpath.split('/')[:-1]), suffix='coil',
									pymol_version=2.0, color_selection=True, helix_order=bundle.helix_order, helices_axis=bundle.helices_axis)

	else:
		print('Unknown mode. Available modes: auto-detect and defdata')
