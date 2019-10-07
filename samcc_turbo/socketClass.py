# DEV
from . import helper_functions, chainClass, layer_detection, bundleClass
import importlib
importlib.reload(helper_functions)
importlib.reload(chainClass)
importlib.reload(layer_detection)
importlib.reload(bundleClass)

from .mierzaczka_turbo import DEBUG
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

### Exceptions

class NoResidue(Exception):
	pass

class TooShort(Exception):
	pass

class socket_class():
	""" """
	#FIXME add docs

	def __init__(self, socket_data, pdb_filename):

		#FIXME self.socket_data is not always containing socket data - should fix name to reflect that

		self.socket_data = socket_data
		self.pdb_filename = pdb_filename


	def get_bundles(self, mode, res_num_layer_detection=5, distance_threshold=50, search_layer_setting_num=3):
		#FIXME update docstring once final version
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

			#print(self.socket_data[1])

			for cc in self.socket_data[0]:
				o = detect_helices_orientation(cc['indices'], self.socket_data[1])
				bundleDesc = []
				for helix in cc['indices']:
					bundleDesc.append([int(helix[1]), int(helix[2]), helix[3], o[int(helix[0])]==-1])
					
				# if all helices in a bundle are AP then we should
				# switch them all to P (parallel)
				
				if all([i[3] for i in bundleDesc]):
					new_bundleDesc=[]
					for i in bundleDesc:
						new_bundleDesc.append(i[:-1] + [False])
					bundleDesc = new_bundleDesc

				bundles.append(bundleDesc)

			return bundles

		def from_deffile():

			#FIXME probably not useful - check and if so then delete
			#FIXME elaborate on deffile format

			"""
			Reads bundle based on classic SamCC definition file.

			Arguments:

				deffile: name of the definition file.

			Example input:

				deffile:

				/path/to/pdb/file/file.pdb
				A,11,29,a
				A,60,42,c
				A,106,88,f
				A,117,135,e


			"""

			#FIXME fix naming, self.socket_data is misleading

			deffile = self.socket_data

			f=open(deffile, 'r')

			filename=f.readline().strip()
			line=f.readline().strip()

			chains=[]
			chains_names = []
			chains_ap = []
			chains_heptad = []

			while line:
				chain, fr, to, heptad = line.split(',')
				chains_heptad.append(heptad)
				fr,to=int(fr), int(to)
				if fr<to:
					#chains.append(range(fr,to+1))
					chains.append((fr,to))
					chains_ap.append(False)
				else:
					#chains.append(range(fr,to-1,-1))
					chains.append((to,fr))
					chains_ap.append(True)

				while chain in chains_names:
					chain += "'"
				chains_names.append(chain)

				line=f.readline().strip()

			f.close()

			#print(filename, chains, chains_names, chains_ap, chains_heptad)

			socket_like_data = [[ (list(indexes)[0], list(indexes)[-1], chain, ap) for indexes, chain, ap in zip(chains, chains_names, chains_ap) ]]
			#socket_like_data = [[ (list(indexes)[0], list(indexes)[-1], chain, ap) for indexes, chain, ap in zip(chains, chains_names, chains_ap) ]]

			return socket_like_data

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


		if mode == 'deffile':
			"""Do not determine order in automatic way but use bundle definition file."""
			#FIXME use directly Bundle class, do not wrap with Socket class
			#FIXME once it works delete this option

			for input_helices in from_deffile():
				# print('IH', input_helices)

				chains = parse_pdb(input_helices, self.pdb_filename)

			for c in chains:
				c.calc_axis(smooth=False)


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

			bundle.pymol_selection = create_pymol_selection_from_socket_results(input_helices)

			#FIXME without helix_order pymol session will not have layers
			# need some function that will exrapolate helix order from deffile/structure
			#bundle.helix_order = [(0, 1), (1, 2), (2, 3), (0, 3)]

			bundles.append(bundle)

		elif mode == 'auto-detect':

			for input_helices in parse_socket():
				if DEBUG:
					print('IH', input_helices)

				chains = parse_pdb(input_helices, self.pdb_filename)

				#FIXME clean from devel, add short description how detection works

				# calculate axis of the bundle
				for c in chains:
					c.calc_axis(smooth=False)

				# get list of CA coords for all chains(helices)
				helices_CA = [c.get_CA() for c in chains]
				if DEBUG:
					print('Helices CA')
					for h in helices_CA:
						print(len(h))
				# helices_CA is list of helices, each is list of CA coords (Biopython objects) [dev-doc]

				# equivalent of "find_best_fit_line_to_helices_CAs"
				# (calculation of helices_axis, helices_axis_all, helices_pts_first, helices_pts_last

				helices_axis_all       = [[list(i.get_array()) for i in c.axis[1:-1]] for c in chains]
				helices_axis_all_store = [[list(i.get_array()) for i in c.axis[1:-1]] for c in chains]

				if DEBUG:
					print('Helices axis all')
					for h in helices_axis_all:
						print(len(h))
				# helices_axis_all is list of helices, each is list of axis coords (3-element list of floats) [dev-doc]

				### DEVEL
				# print('HLX AXIS ALL')
				# hlx_temp = ([(c[0], c[-1]) for c in helices_axis_all])
				# for c in hlx_temp:
				# 	print('CHAIN', c)
				from scipy.spatial import distance
				helices_axis_all_cut     = [ [] for _x in range(len(helices_axis_all)) ]
				helices_axis_all_cut_pos = [ [] for _x in range(len(helices_axis_all)) ]

				# helices_axis_all_cut is created as list of helices, each helix empty list [dev-doc]
				# helices_axis_all_cut_pos is created as list of helices, each helix empty list [dev-doc]
				# helices_axis_all_cut_pos is used for pymol-readable selection generation [dev-docs]

				if DEBUG:
					print('Helices axis all cut')
					print(helices_axis_all_cut)
					for h in helices_axis_all_cut:
						print(len(h))

				#FIXME here the search is done by comparing distances of points on helices axis - change var names and docs to reflect that
				# code verifying if all residues have any residue from other helix within reasonable distance (20A) [dev-doc]
				# this code compares all vs all residues from all helices
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
										# print('FOUND', dst)
										# found residue within cutoff distance on this helix, break and search next helix
										res_stat.append(True)
										break
								else:
									# there is no residue within cutoff distance on this helix, search next one
									res_stat.append(False)
						# print(res_stat)
						if all(res_stat):
							# there is residue within cutoff distance on every other helix
							# add this residue coords to its helix (use empty-filled helices_axis_all_cut)
							# add info whether this residue was accepted to its helix (use empty-filled helices_axis_all_cut_pos)
							helices_axis_all_cut[h1[0]].append(res1)
							helices_axis_all_cut_pos[h1[0]].append(True)
						else:
							# just add info whether this residue was accepted to its helix (use empty-filled helices_axis_all_cut_pos)
							helices_axis_all_cut_pos[h1[0]].append(False)
						# print('---')

				if DEBUG:
					print('Helices axis all cut - after loop')
					for h in helices_axis_all_cut:
						print(len(h))
					print('Helices axis all cut pos - after loop')
					for h in helices_axis_all_cut_pos:
						print(len(h))

				# update helices_axis_all so it contains only residues accepted by the distance filter [dev-doc]
				# could this create helices of different lengths? [dev-issue] - yes it can
				helices_axis_all = helices_axis_all_cut
				# print('CUT')
				# for c in helices_axis_all:
				# 	print(len(c))
				### END-DEVEL

				#FIXME this probably was used to draw approximated helix axis in pymol, now not used, use or delete (it is still arg in drawing function)
				helices_axis = [(c[0], c[-1]) for c in helices_axis_all]

				# FIXME possible issue (next block too):
				# here shortest_helix_len is determined from chains not from helix_axis_all_cut so this value can be larger than actual measured helix [dev-doc]
				shortest_helix_len = min([len(c.axis)-2 for c in chains])
				if DEBUG:
					print('shortest len', shortest_helix_len)

				# this block checks if the length of shortest helix is bigger than number of residues on helix that will be searched [dev-doc]
				if shortest_helix_len < res_num_layer_detection:
					raise TooShort
					#res_num_layer_detection_asserted = shortest_helix_len
					#print('Warning: one of the helices was too short to start layer search with res_num_layer_detection= ' + str(res_num_layer_detection) + '; searching with res_num_layer_detection=' + str(res_num_layer_detection_asserted))
				else:
					res_num_layer_detection_asserted = res_num_layer_detection

				# from each helix n positions on its axis are extracted; n = res_num_layer_detection parameter [dev-doc]
				# since helix_axis_all is basis for extraction these helices are already cut using minimal distance between residues [dev-doc]
				# points are extracted from 3 segments: beginning, middle and end of a bundle [dev-doc]
				helices_pts_first = [c[:res_num_layer_detection_asserted] for c in helices_axis_all]
				helices_pts_last = [c[-res_num_layer_detection_asserted:] for c in helices_axis_all]
				# middle-out helices_pts
				# for each helix take n/n-1 points from middle of a given helix; n when n is even and n-1 when n is uneven
				middle_points = [ int(len(helix) / 2) for helix in helices_axis_all ]
				helices_pts_middle = [ c[1][middle_points[c[0]]-int(res_num_layer_detection_asserted/2):middle_points[c[0]]+int(res_num_layer_detection_asserted/2)] for c in enumerate(helices_axis_all) ]

				if DEBUG:
					print('helices pts first')
					for h in helices_pts_first:
						print(len(h))
					print('middle points', middle_points)
					print('helices pts middle')
					for h in helices_pts_middle:
						print(len(h))


				# detect helix order - and propagate it to all boundry layer searches
				if DEBUG:
					print('helices pts middle type', type(helices_pts_middle))
				#FIXME ppo is used only to draw middle points of helices in pymol, probably should be excluded in production version
				#FIXME examine detect_helix_order function and add docs
				neighbour_interactions, ppo = detect_helix_order(helices_pts_middle, helices_pts_all=helices_axis_all)
				#FIXME clean old code if no longer needed
				#neighbour_interactions = ld.detect_helix_order(helices_pts_first, helices_pts_last, helices_pts_middle)
				if DEBUG:
					print('neighbour interactions', neighbour_interactions)

				# find first layer starting from different segments of coil, these segments are defined from helix_axis_all after minimal distance filter (see above) [dev-doc]
				#FIXME set distance flag and check its value after each segment: if False then skip to layer detection rerun with higher res_num_layer_detection
				#FIXME it step may be computational costly so it will speed up script
				first_layer, first_distance       = find_bundle_boundry_layer_universal(helices_pts_first, neighbour_interactions,
													distance_threshold=distance_threshold, search_layer_setting_num=search_layer_setting_num)

				#FIXME find CA on given helix that is closest to the point on axis - maybe this should find CA related to this point not the closest one
				#FIXME i.e. CA that was used to compute this point on axis (if there was just one) - check how exactly axis is calculated
				if DEBUG:
					print('helices CA - to find closest CA')
					for h in helices_CA:
						print(len(h))

				#FIXME: n-distance: this should operate on first_layer as a list now
				#FIXME: single min-dist version:
				#CA_first_layer, first_layer_ids   = find_closest_CA_to_point(first_layer, helices_CA)
				CA_first_layer_ids   = [ find_closest_CA_to_point(layer, helices_CA) for layer in first_layer ]

				last_layer, last_distance         = find_bundle_boundry_layer_universal(helices_pts_last, neighbour_interactions,
													distance_threshold=distance_threshold, search_layer_setting_num=search_layer_setting_num)
				#CA_last_layer, last_layer_ids     = find_closest_CA_to_point(last_layer, helices_CA)
				CA_last_layer_ids   = [ find_closest_CA_to_point(layer, helices_CA) for layer in last_layer ]

				middle_layer, middle_distance     = find_bundle_boundry_layer_universal(helices_pts_middle, neighbour_interactions,
													distance_threshold=distance_threshold, search_layer_setting_num=search_layer_setting_num)
				#CA_middle_layer, middle_layer_ids = find_closest_CA_to_point(middle_layer, helices_CA)
				CA_middle_layer_ids = [ find_closest_CA_to_point(layer, helices_CA) for layer in middle_layer ]

				### DEVEL
				### if none of layer settings return reasonable distance: search again with bigger part of the bundle
				#FIXME should probably be iterative search with various res_num_layer_detection values

				if not any([first_distance, last_distance, middle_distance]):
					# print(first_distance, last_distance, middle_distance)
					res_num_layer_detection_asserted = 15
					if shortest_helix_len < res_num_layer_detection_asserted:
						res_num_layer_detection_asserted = shortest_helix_len

					# print(shortest_helix_len)

					helices_pts_first = [c[:res_num_layer_detection_asserted] for c in helices_axis_all]
					helices_pts_last = [c[-res_num_layer_detection_asserted:] for c in helices_axis_all]
					# middle-out helices_pts
					middle_points = [ int(len(helix) / 2) for helix in helices_axis_all ]
					helices_pts_middle = [c[1][middle_points[c[0]]-int(res_num_layer_detection_asserted/2):middle_points[c[0]]+int(res_num_layer_detection_asserted/2)] for c in enumerate(helices_axis_all)]

					# detect helix order - and propagate it to all boundry layer searches
					#neighbour_interactions = ld.detect_helix_order(helices_pts_first, helices_pts_last, helices_pts_middle)
					neighbour_interactions, ppo = detect_helix_order(helices_pts_middle, helices_pts_all=helices_axis_all)
					# print(neighbour_interactions)

					# find first layer starting from different segments of coil

					#FIXME: n-distance: this should operate on first_layer as a list now
					#FIXME: single min-dist version:
					#CA_first_layer, first_layer_ids   = find_closest_CA_to_point(first_layer, helices_CA)
					first_layer, first_distance       = find_bundle_boundry_layer_universal(helices_pts_first, neighbour_interactions,
														distance_threshold=distance_threshold, search_layer_setting_num=search_layer_setting_num)
					CA_first_layer_ids   = [ find_closest_CA_to_point(layer, helices_CA) for layer in first_layer ]

					last_layer, last_distance         = find_bundle_boundry_layer_universal(helices_pts_last, neighbour_interactions)
					#CA_last_layer, last_layer_ids     = find_closest_CA_to_point(last_layer, helices_CA)
					CA_last_layer_ids   = [ find_closest_CA_to_point(layer, helices_CA) for layer in last_layer ]

					middle_layer, middle_distance     = find_bundle_boundry_layer_universal(helices_pts_middle, neighbour_interactions)
					#CA_middle_layer, middle_layer_ids = find_closest_CA_to_point(middle_layer, helices_CA)
					CA_middle_layer_ids = [ find_closest_CA_to_point(layer, helices_CA) for layer in middle_layer ]

					# print('CORRECTED:', first_distance, last_distance, middle_distance)

				### end search again
				### END-DEVEL

				#FIXME does not seem to be used - examine and delete if not needed
				# CA_layers                         = [CA_first_layer, CA_last_layer, CA_middle_layer]

				# find all layers starting from all found points
				#FIXME: singular version:
				#CA_layers_from_first   = find_all_layers_from_layer(first_layer_ids, helices_CA)
				#FIXME n-distance:
				CA_layers_from_first   = [ find_all_layers_from_layer(layer[1], helices_CA) for layer in CA_first_layer_ids ]
				CA_layers_from_last    = [ find_all_layers_from_layer(layer[1], helices_CA) for layer in CA_last_layer_ids ]
				CA_layers_from_middle  = [ find_all_layers_from_layer(layer[1], helices_CA) for layer in CA_middle_layer_ids ]

				# check average angle between layers and select layer setting with smaller avg angle for samcc input
				#FIXME n-distance: now selects from multipe starting layers for 3 point in bundle
				CA_layers_all                  = [CA_layers_from_first, CA_layers_from_last, CA_layers_from_middle]
				#FIXME singular version:
				#layer_ids_all                  = [(first_layer_ids,'both'), (last_layer_ids, 'both'), (middle_layer_ids, 'both')]
				layer_ids_all                  = [ [ (layer[1],'both') for layer in CA_first_layer_ids ],
												   [ (layer[1],'both') for layer in CA_last_layer_ids ],
												   [ (layer[1],'both') for layer in CA_middle_layer_ids ] ]

				#FIXME implement new method for dimers
				#FIXME cast helices onto plane and measure angle between points on first helix and corresponding points on second helix
				#FIXME n-distance: now selects from multiple starting layers
				# case of dimers (only 2 helices)
				if len(CA_layers_from_first[0]) == 2:
					best_layer_set, best_layer_ids = select_minimal_distance_layer_set(CA_layers_all, layer_ids_all)
				# all other cases
				else:
					best_layer_set, best_layer_ids = select_minimal_angle_layer_set(CA_layers_all, layer_ids_all)

				# at this point best layer setting is selected and list of ids of first residues at each of its helices is stored in best_layer_ids [doc-docs]
				# best layer_ids will be then used to cut fragments of chains to pass them to the bundleClass instance [dev-docs]
				if DEBUG:
					print ('BEST LAYER', best_layer_ids)

				bundle_length = len(best_layer_set)
				len_upstream = min(best_layer_ids[0])

				pymol_temp = []

				for c, b in zip(chains, best_layer_ids[0]):

					cut = b - len_upstream
					#print ("!", len(c.res), len(c.res[cut:cut+bundle_length]), len(c.axis))

					if DEBUG:
						print('cut selected layers; cut:cut+bundle_len', cut, cut+bundle_length)
					c.res = c.res[cut:cut+bundle_length]
					c.axis = c.axis[cut:cut+bundle_length]

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
				#FIXME this only stores the area detected by Socket as CC - it does not affect range of drawn layers in pymol
				#FIXME there was idea to cut it in case of unreasonable Socket results (not paired helices)
				#FIXME decide what Socket detected range should be drawn and implement solution

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

				if DEBUG:
					print('pymol session data', input_helices_cut)
				input_helices = input_helices_cut
				# for c in chains:
				# 	for r in c.res:
				# 		print(r)
				### END-DEVEL

				#FIXME pymol session in final version should indicate Socket detected region and show layer bars in actually measured layers (from bundleClass chains)
				#FIXME ppo probably not needed in production
				bundle.pymol_selection = create_pymol_selection_from_socket_results(input_helices)
				bundle.ppo = ppo

				# save order of helices (helix indexes as in socket data) [ks mode only now]
				if mode == 'auto-detect':
					bundle.helix_order = neighbour_interactions

				bundles.append(bundle)

		return bundles
