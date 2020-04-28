### socketClass definition ###

import numpy as np
from .helper_functions import detect_helices_orientation
from .searchLayerClass import helixAxisBundleClass
from .residueClass import residueClass
from .chainClass import chainClass
from .bundleClass import bundleClass
from .layerClass import layerClass
from Bio.PDB import PDBParser
from .layer_detection import select_minimal_distance_layer_set
from .layer_detection import select_minimal_angle_layer_set2
from .layer_detection import select_minimal_angle_layer_set
from .layer_detection import find_bundle_boundry_layer_from_all
from .layer_detection import create_pymol_selection_from_socket_results
from .layer_detection import select_minimal_total_distance_layer_set
from .layer_detection import select_minimal_dist_to_plane_layer_set
from .layer_detection import get_layers_set_with_best_ranks
from .layer_detection import find_closest_CA_to_point

### Exceptions

class NoResidue(Exception):
	pass

class TooShort(Exception):
	pass

class socket_class():
	"""Class responsible for handling Socket-derived data"""

	def __init__(self, socket_data, pdb_filename):

		self.socket_data = socket_data
		self.pdb_filename = pdb_filename

	def get_bundles(self, res_num_layer_detection=5, distance_threshold='auto', search_layer_setting_num=9):
		"""Given a SOCKET raw data calculates a proper bundle(s) in which all
		helices are of the same length and residues are arrranged in layers.

		Arguments (same description as in main SamCC-Turbo function):
		res_num_layer_detection   -- number of residues from each helix that will be combined
						  		     into layers that will be considered in layer detection;
						             (default 5) [this is "r" parameter in the paper]
		distance_threshold        -- distance threshold [A] that will be used to indicate unusual
							  		 layer setting; if total distance between residues in layer
									 if greater that this number the layer will be discarded, if
									 no layer will pass this filter then layer detection will
									 be repeated with bigger layer_detect_n; if set to auto
									 then threshold will be assigned basing on oligomerization
									 (default auto)
		search_layer_setting_num  -- number of starting layer settings that will be selected
									 basing on minimal distance between residues; they will be compared
									 in terms of angle between all layers in structure (default 9)

		Returns:
			bundles: list of bundleClass instances (one structure can contain multiple bundles)
		"""

		def parse_socket():
			"""Parse raw SOCKET data and generate a list of helices

			Arguments:
			socket_data -- raw SOCKET data

			Returns:
			bundles: list of bundles (i.e list of helices: from, to, chain, anti-parallel)
			"""

			bundles = []

			for cc_id in self.socket_data:
				cc = self.socket_data[cc_id]
				# detect helices orientation
				o = detect_helices_orientation(cc['indices'], cc['relations'])

				bundleDesc = []
				for helix_id in cc['indices']:
					helix = cc['indices'][helix_id]
					bundleDesc.append([helix['start'], helix['end'], helix['chain'], o[helix_id]==-1])

				# if all helices in a bundle are AP then we should
				# switch them all to P (parallel)

				if all([i[3] for i in bundleDesc]):
					new_bundleDesc=[]
					for i in bundleDesc:
						new_bundleDesc.append(i[:-1] + [False])
					bundleDesc = new_bundleDesc

				bundles.append(bundleDesc)

			return bundles

		def parse_pdb(input_helices, pdb_filename):
			"""
			Parses PDB file according to SOCKET definitions

			Arguments:
				input_helices: list of helices (from, to, chain, ap)

			Returns:
				chains: list of chainClass instances
						(note that these chains may have different lengths)
			"""

			parser=PDBParser(QUIET=True)
			structure=parser.get_structure('temp', pdb_filename)

			chains = []

			for helix in input_helices:
				start, end, chain, ap = helix

				# solve the problem with het atoms
				preres = {}
				try:
					for i in structure[0][chain]:
						j = i.get_full_id()
						preres[j[3][1]] = i
				except KeyError:
					raise NoResidue()

				try:
					residues = [residueClass(preres[r], chain) for r in range(start, end+1)]
				except KeyError:
					raise NoResidue

				if ap == True:
					residues = residues[::-1]

				#print(len(residues), ap, chain)
				chain = chainClass(residues, ap, chain)
				chain.new_starts=[]
				chains.append(chain)

			assert len(chains)>1
			return chains

		bundles = []

		for input_helices in parse_socket():

			try:
				chains = parse_pdb(input_helices, self.pdb_filename)
			except NoResidue:
				print('No residue error')
				continue

			# assign distance_threshold according to oligomerization state
			default_distance_threshold = {2:30, 3:40, 4:50, 5:60, 6:70, 7:80, 8:90, 9:100}
			if distance_threshold == 'auto':
				distance_threshold_set = default_distance_threshold[len(input_helices)]
			else:
				distance_threshold_set = distance_threshold

			# calculate axis of the bundle
			for c in chains:
				c.calc_axis(smooth=False)

			# initialize axis bundle from chain axis points
			helices_axis_all = helixAxisBundleClass(chains)

			# flag axis points that are too far away from other points
			helices_axis_all.verify_points_distances(max_distance=20)

			# initialize boundry_layers list and search scope
			boundry_layers                   = []
			res_num_layer_detection_asserted = res_num_layer_detection

			# condition: any layer (from all segments) has distance belowe threshold
			while(not boundry_layers):

				# determine helix order
				helices_axis_all.detect_helix_order(res_num_layer_detection_asserted)

				# subset axis points from start, middle and end of bundle (according to distance flag)
				helices_pts_first  = helices_axis_all[:res_num_layer_detection_asserted]
				helices_pts_last   = helices_axis_all[-res_num_layer_detection_asserted:]
				helices_pts_middle = helices_axis_all.get_middle_points(res_num_layer_detection_asserted)

				for helices_pts in [helices_pts_first, helices_pts_middle, helices_pts_last]:
					boundry_layers += helices_pts.get_all_bundle_boundry_layers()

				boundry_layers = find_bundle_boundry_layer_from_all(boundry_layers, distance_threshold=distance_threshold_set, search_layer_setting_num=search_layer_setting_num)

				# increase scope of searched residues
				# INCREMENT - to activate multiloop uncomment next line and comment out four following lines
				# res_num_layer_detection_asserted += 1
				if res_num_layer_detection_asserted == 15:
					break
				else:
					res_num_layer_detection_asserted = 15
				# INCREMENT END

			assert len(boundry_layers)>0

			# merge all boundry layers into one list - but only selected ones
			percentile_threshold = np.percentile([ bl.total_distance for bl in boundry_layers ], 50)
			boundry_layers       = [ bl for bl in boundry_layers if bl.total_distance < percentile_threshold ]

			layers_sets = helices_axis_all.find_all_layers_from_layer(boundry_layers)

			# dimer case
			if len(input_helices) == 2:
				best_layer_set = select_minimal_distance_layer_set(layers_sets)
			# all other cases
			else:
				best_layer_set = select_minimal_angle_layer_set(layers_sets, best_layer_nb=1)

			# convert best layer set to bundleClass()
			# -1/+1 to chain indexes to accomodate fact that chain edges are excluded from layer search
			# this is because it is not possible to calculate helix axis points for edge residues
			for c, b in zip(chains, best_layer_set.iterchain()):
				cut_start = b[0].point_id -1 # id of first residue in best layer chain
				cut_stop  = b[-1].point_id +1 # id of last residue in best layer chain
				c.res = c.res[cut_start:cut_stop+1]
				c.axis = c.axis[cut_start:cut_stop+1]

			bundle = bundleClass()

			# check whether all chains are of equal length
			temp = [len(i.res) for i in chains]
			assert all([i==temp[0] for i in temp])

			# final definition of helices
			bundle.chains = chains

			# generate layers
			bundle.layers = []
			for layer_residues in zip(*[c.res for c in bundle.chains]):
				bundle.layers.append(layerClass(layer_residues))

			# generate pymol-readable selection
			# this will cut range of shown helices to the residues that are within resonable distance from all other helices
			# it does not affect range of drawn layers in pymol
			input_helices_cut = []
			for ha, ih in zip(helices_axis_all[:].iterchain(), input_helices):
				if ih[3]:
					input_helices_cut.append([ha[-1].res.get_id()[1], ha[0].res.get_id()[1], ih[2], ih[3]])
				else:
					input_helices_cut.append([ha[0].res.get_id()[1], ha[-1].res.get_id()[1], ih[2], ih[3]])

			input_helices = input_helices_cut

			bundle.pymol_selection = create_pymol_selection_from_socket_results(input_helices)
			bundle.ppo = best_layer_set.ppo

			# save order of helices (helix indexes as in socket data)
			bundle.helix_order = best_layer_set.neighbour_interactions

			bundles.append(bundle)

		return bundles
