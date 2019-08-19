# DEV
from . import helper_functions, chainClass, layer_detection, bundleClass
import importlib
importlib.reload(helper_functions)
importlib.reload(chainClass)
importlib.reload(layer_detection)
importlib.reload(bundleClass)
# END-DEV

from .helper_functions import detect_helices_orientation
from .residueClass import residueClass
from .chainClass import chainClass
from .bundleClass import bundleClass
from .layerClass import layerClass
from Bio.PDB import PDBParser
from .layer_detection import detect_helix_order, find_bundle_boundry_layer_universal
from .layer_detection import find_closest_CA_to_point
from .layer_detection import find_all_layers_from_layer
from .layer_detection import select_minimal_distance_layer_set
from .layer_detection import select_minimal_angle_layer_set
from .layer_detection import create_pymol_selection_from_socket_results

class socket_class():
	""" """
	#FIXME add docs

	def __init__(self, socket_data, pdb_filename):
		self.socket_data = socket_data
		self.pdb_filename = pdb_filename


	def get_bundles(self, mode, res_num_layer_detection=5):
		"""
		Given a SOCKET raw data calculates a proper bundle(s) in which all helices are of the same length
		and residues are arrranged in layers.

		Arguments:
			mode: "ks" for Krzysiek's approach and "sdh" for Staszek's approach
			res_num_layer_detection: used only in the "ks" mode. states how many residues from each
									 end will be considered in search for the first layer (int)

		Returns:
			bundles: list of bundleClass instances (one structure can contain multiple bundles)

		"""


		def parse_socket():
			"""
			Parse raw SOCKET data and generate a list of helices

			Arguments:
				socket_data: raw SOCKET data

			Returns:
				bundles: list of bundles (i.e list of helices: from, to, chain, anti-parallel)
			"""

			bundles = []
			for cc in self.socket_data[0]:
				o = detect_helices_orientation(cc['indices'], self.socket_data[1])
				bundleDesc = []
				for helix in cc['indices']:
					bundleDesc.append((int(helix[1]), int(helix[2]), helix[3], o[int(helix[0])]==-1))
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

			parser=PDBParser()
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

		def cut2(ch1, ch2):

			#FIXME is this function used? if not: delete

			ch1_O = ch1.res[1].O
			ch2_O = ch2.res[1].O

			ch1_to_all_ch2_dist = [(ch1_O - r.O).norm() for r in ch2.res[1:-1]]
			ch2_to_all_ch1_dist = [(ch2_O - r.O).norm() for r in ch1.res[1:-1]]

			if min(ch1_to_all_ch2_dist) < min(ch2_to_all_ch1_dist):
				ch2.new_starts.append(str(np.argmin(ch1_to_all_ch2_dist))+"_"+ch1.chain)
			else:
				ch1.new_starts.append(str(np.argmin(ch2_to_all_ch1_dist))+"_"+ch2.chain)

			return ch1, ch2

		# --------------------------------------------------------------------------------

		bundles = []

		for input_helices in parse_socket():
			# print(input_helices)

			chains = parse_pdb(input_helices, self.pdb_filename)

			if mode == 'sdh':
				"""
				Layer detection implemented by Staszek
				"""

				for c in chains:
					c.calc_axis()

				# compare helices all-vs-all
				for pos1 in range(1, len(chains)):
					chain1 = chains[pos1]
					for pos2 in range(0,pos1):
						chain2 = chains[pos2]
						chain1, chain2 = cut2(chain1, chain2)

				# cut the helices
				for c in chains:
					if len(c.new_starts)>0:
						c.res = c.res[max([int(i.split("_")[0]) for i in c.new_starts]):]

				# cut to the shortest helix
				minlen = min([len(c.res) for c in chains])

				for c in chains:
					c.res = c.res[:minlen]
					if len(c.res)<4: raise TooShort

				assert all([len(c.res)==len(chains[0].res) for c in chains])

			elif mode == 'ks':
				#FIXME clean from devel, add short description how detection works

				for c in chains:
					c.calc_axis(smooth=False)

				helices_CA = [c.get_CA() for c in chains]

				# equivalent of "find_best_fit_line_to_helices_CAs"
				# (calculation of helices_axis, helices_axis_all, helices_pts_first, helices_pts_last

				helices_axis_all       = [[list(i.get_array()) for i in c.axis[1:-1]] for c in chains]
				helices_axis_all_store = [[list(i.get_array()) for i in c.axis[1:-1]] for c in chains]

				### DEVEL
				# print('HLX AXIS ALL')
				# hlx_temp = ([(c[0], c[-1]) for c in helices_axis_all])
				# for c in hlx_temp:
				# 	print('CHAIN', c)
				from scipy.spatial import distance
				helices_axis_all_cut     = [ [] for _x in range(len(helices_axis_all)) ]
				helices_axis_all_cut_pos = [ [] for _x in range(len(helices_axis_all)) ]

				for h1 in enumerate(helices_axis_all):
					# print('----------')
					# print('HLX', h1[0])
					for res1 in h1[1]:
						# print('RES')
						res_stat = []
						for h2 in enumerate(helices_axis_all):
							if h1[0] != h2[0]:
								for res2 in h2[1]:
									dst = distance.euclidean(res1, res2)
									if dst < 20:
										#print('FOUND', dst)
										res_stat.append(True)
										break
								else:
									res_stat.append(False)
						# print(res_stat)
						if all(res_stat):
							helices_axis_all_cut[h1[0]].append(res1)
							helices_axis_all_cut_pos[h1[0]].append(True)
						else:
							helices_axis_all_cut_pos[h1[0]].append(False)
						# print('---')

				helices_axis_all = helices_axis_all_cut
				# print('CUT')
				# for c in helices_axis_all:
				# 	print(len(c))
				### END-DEVEL

				helices_axis = [(c[0], c[-1]) for c in helices_axis_all]

				shortest_helix_len = min([len(c.axis)-2 for c in chains])

				if shortest_helix_len < res_num_layer_detection:
					raise TooShort
					#res_num_layer_detection_asserted = shortest_helix_len
					#print('Warning: one of the helices was too short to start layer search with res_num_layer_detection= ' + str(res_num_layer_detection) + '; searching with res_num_layer_detection=' + str(res_num_layer_detection_asserted))
				else:
					res_num_layer_detection_asserted = res_num_layer_detection

				helices_pts_first = [c[:res_num_layer_detection_asserted] for c in helices_axis_all]
				helices_pts_last = [c[-res_num_layer_detection_asserted:] for c in helices_axis_all]
				# middle-out helices_pts
				# middle_point = int(len(helices_axis_all[0]) / 2) # OLD VERSION
				# helices_pts_middle = [c[middle_point-int(res_num_layer_detection_asserted/2):middle_point+int(res_num_layer_detection_asserted/2)] for c in helices_axis_all] # OLD VERSION
				middle_points = [ int(len(helix) / 2) for helix in helices_axis_all ]
				helices_pts_middle = [c[1][middle_points[c[0]]-int(res_num_layer_detection_asserted/2):middle_points[c[0]]+int(res_num_layer_detection_asserted/2)] for c in enumerate(helices_axis_all)]

				# detect helix order - and propagate it to all boundry layer searches
				neighbour_interactions, ppo = detect_helix_order(helices_pts_middle, helices_pts_all=helices_axis_all)
				#neighbour_interactions = ld.detect_helix_order(helices_pts_first, helices_pts_last, helices_pts_middle)
				# print(neighbour_interactions)

				# find first layer starting from different segments of coil
				first_layer, first_distance       = find_bundle_boundry_layer_universal(helices_pts_first, neighbour_interactions)
				CA_first_layer, first_layer_ids   = find_closest_CA_to_point(first_layer, helices_CA)

				last_layer, last_distance         = find_bundle_boundry_layer_universal(helices_pts_last, neighbour_interactions)
				CA_last_layer, last_layer_ids     = find_closest_CA_to_point(last_layer, helices_CA)

				middle_layer, middle_distance     = find_bundle_boundry_layer_universal(helices_pts_middle, neighbour_interactions)
				CA_middle_layer, middle_layer_ids = find_closest_CA_to_point(middle_layer, helices_CA)

				### DEVEL
				### if none of layer settings return reasonable distance: search again with bigger part of the bundle

				if not any([first_distance, last_distance, middle_distance]):
					# print(first_distance, last_distance, middle_distance)
					res_num_layer_detection_asserted = 15
					if shortest_helix_len < res_num_layer_detection_asserted:
						res_num_layer_detection_asserted = shortest_helix_len

					# print(shortest_helix_len)

					helices_pts_first = [c[:res_num_layer_detection_asserted] for c in helices_axis_all]
					helices_pts_last = [c[-res_num_layer_detection_asserted:] for c in helices_axis_all]
					# middle-out helices_pts
					# middle_point = int(len(helices_axis_all[0]) / 2) # OLD VERSION
					# helices_pts_middle = [c[middle_point-int(res_num_layer_detection_asserted/2):middle_point+int(res_num_layer_detection_asserted/2)] for c in helices_axis_all] # OLD VERSION
					middle_points = [ int(len(helix) / 2) for helix in helices_axis_all ]
					helices_pts_middle = [c[1][middle_points[c[0]]-int(res_num_layer_detection_asserted/2):middle_points[c[0]]+int(res_num_layer_detection_asserted/2)] for c in enumerate(helices_axis_all)]

					# detect helix order - and propagate it to all boundry layer searches
					#neighbour_interactions = ld.detect_helix_order(helices_pts_first, helices_pts_last, helices_pts_middle)
					neighbour_interactions, ppo = detect_helix_order(helices_pts_middle, helices_pts_all=helices_axis_all)
					# print(neighbour_interactions)

					# find first layer starting from different segments of coil
					first_layer, first_distance       = find_bundle_boundry_layer_universal(helices_pts_first, neighbour_interactions)
					CA_first_layer, first_layer_ids   = find_closest_CA_to_point(first_layer, helices_CA)

					last_layer, last_distance         = find_bundle_boundry_layer_universal(helices_pts_last, neighbour_interactions)
					CA_last_layer, last_layer_ids     = find_closest_CA_to_point(last_layer, helices_CA)

					middle_layer, middle_distance     = find_bundle_boundry_layer_universal(helices_pts_middle, neighbour_interactions)
					CA_middle_layer, middle_layer_ids = find_closest_CA_to_point(middle_layer, helices_CA)

					# print('CORRECTED:', first_distance, last_distance, middle_distance)

				### end search again
				### END-DEVEL

				CA_layers                         = [CA_first_layer, CA_last_layer, CA_middle_layer]

				# find all layers starting from all found points
				CA_layers_from_first   = find_all_layers_from_layer(first_layer_ids, helices_CA)
				CA_layers_from_last    = find_all_layers_from_layer(last_layer_ids, helices_CA)
				CA_layers_from_middle  = find_all_layers_from_layer(middle_layer_ids, helices_CA)

				# check average angle between layers and select layer setting with smaller avg angle for samcc input
				CA_layers_all                  = [CA_layers_from_first, CA_layers_from_last, CA_layers_from_middle]
				layer_ids_all                  = [(first_layer_ids,'both'), (last_layer_ids, 'both'), (middle_layer_ids, 'both')]

				# case of dimers (only 2 helices)
				if len(CA_layers_from_first[0]) == 2:
					best_layer_set, best_layer_ids = select_minimal_distance_layer_set(CA_layers_all, layer_ids_all)
				# all other cases
				else:
					best_layer_set, best_layer_ids = select_minimal_angle_layer_set(CA_layers_all, layer_ids_all)

				# print ('BEST LAYER', best_layer_ids)

				bundle_length = len(best_layer_set)
				len_upstream = min(best_layer_ids[0])

				for c, b in zip(chains, best_layer_ids[0]):

					cut = b - len_upstream
					#print ("!", len(c.res), len(c.res[cut:cut+bundle_length]), len(c.axis))

					c.res = c.res[cut:cut+bundle_length]
					c.axis = c.axis[cut:cut+bundle_length]

			else:
				print ('unknown mode %s' % mode)

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

			### DEVEL
			input_helices_cut = []
			# print(input_helices)
			# for c in helices_axis_all_store:
			# 	print(len(c))
			for c in enumerate(helices_axis_all_cut_pos):
				# print(c[1])
				if input_helices[c[0]][3]:
					res_nums = [ r for r in range(input_helices[c[0]][0]+1, input_helices[c[0]][1])][::-1]
				else:
					res_nums = [ r for r in range(input_helices[c[0]][0], input_helices[c[0]][1]-1)]
				# print(res_nums)
				res_num_cut = [ r[1] for r in enumerate(res_nums) if c[1][r[0]] ]
				# print(res_num_cut, '\n')
				input_helices_cut.append((min(res_num_cut), max(res_num_cut), input_helices[c[0]][2], input_helices[c[0]][3]))

			# print(input_helices_cut)
			input_helices = input_helices_cut
			# for c in chains:
			# 	for r in c.res:
			# 		print(r)
			### END-DEVEL

			bundle.pymol_selection = create_pymol_selection_from_socket_results(input_helices)
			bundle.ppo = ppo

			# save order of helices (helix indexes as in socket data) [ks mode only now]
			if mode == 'ks':
				bundle.helix_order = neighbour_interactions

			bundles.append(bundle)

		return bundles
