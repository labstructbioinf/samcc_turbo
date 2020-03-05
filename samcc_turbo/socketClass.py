# DEV
from . import helper_functions, chainClass, layer_detection, bundleClass, searchLayerClass
from operator import attrgetter
import numpy as np
import importlib
import heapq
importlib.reload(helper_functions)
importlib.reload(chainClass)
importlib.reload(layer_detection)
importlib.reload(bundleClass)
importlib.reload(searchLayerClass)

from .mierzaczka_turbo import DEBUG
# END-DEV

from .helper_functions import detect_helices_orientation
from .searchLayerClass import helixAxisBundleClass
from .residueClass import residueClass
from .chainClass import chainClass
from .bundleClass import bundleClass
from .layerClass import layerClass
from Bio.PDB import PDBParser
from .layer_detection import select_minimal_distance_layer_set
from .layer_detection import select_minimal_angle_layer_set2, select_minimal_angle_layer_set
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
	""" """
	#FIXME add docs

	def __init__(self, socket_data, pdb_filename):

		#FIXME self.socket_data is not always containing socket data - should fix name to reflect that

		self.socket_data = socket_data
		self.pdb_filename = pdb_filename

	def get_bundles(self, mode, res_num_layer_detection=5, distance_threshold='auto', search_layer_setting_num=3):
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
				#if DEBUG:
				# print('IH', input_helices)
				# print('LEN IH', len(input_helices))

				try:
					chains = parse_pdb(input_helices, self.pdb_filename)
					# for c in chains:
					# 	print(c)
				except NoResidue:
					print('No residue error')
					continue

				#FIXME clean from devel, add short description how detection works

				# assign distance_threshold according to oligomerization state
				#default_distance_threshold = {2:30, 3:40, 4:50, 5:60, 6:70, 7:80, 8:90, 9:100}
				default_distance_threshold = {2:30, 3:40, 4:50, 5:60, 6:70, 7:80, 8:90, 9:100}
				if distance_threshold == 'auto':
					distance_threshold_set = default_distance_threshold[len(input_helices)]
				else:
					distance_threshold_set = distance_threshold

				# calculate axis of the bundle
				for c in chains:
					c.calc_axis(smooth=False)

				# get list of CA coords for all chains(helices)
				# helices_CA = [c.get_CA() for c in chains]
				# if DEBUG:
				# 	print('Helices CA')
				# 	for h in helices_CA:
				# 		print(len(h))
				# helices_CA is list of helices, each is list of CA coords (Biopython objects) [dev-doc]

				# equivalent of "find_best_fit_line_to_helices_CAs"
				# (calculation of helices_axis, helices_axis_all, helices_pts_first, helices_pts_last

				### DEVEL

				#imports - move up once finish devel


				# initialize axis bundle from chain axis points
				helices_axis_all = helixAxisBundleClass(chains)


				# helices_axis_all.show_points()
				# print('====='*15)

				# flag axis points that are too far away from other points
				helices_axis_all.verify_points_distances(max_distance=20)
				# helices_axis_all.show_points()
				# print('====='*15)

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

					# if DEBUG:
					# helices_pts_first.show_points()
					# helices_pts_middle.show_points()
					# helices_pts_last.show_points()

					# helices_pts_first.show_points() #FIXME devel-code

					# find bundle boundry layer (min distance layer)
					# set of possible best layers in form of list of searchLayer
					# if list is empty then increment search scope



					for helices_pts in [helices_pts_first, helices_pts_middle, helices_pts_last]:
					# for helices_pts in [helices_pts_first, helices_pts_middle]:
						boundry_layers += helices_pts.get_all_bundle_boundry_layers()
					boundry_layers = find_bundle_boundry_layer_from_all(boundry_layers, distance_threshold=distance_threshold_set, search_layer_setting_num=search_layer_setting_num)

					# old non-merge version, delete to production
					# for helices_pts in [helices_pts_first, helices_pts_middle, helices_pts_last]:
					# # for helices_pts in [helices_pts_first, helices_pts_middle]:
					# 	boundry_layers += helices_pts.find_bundle_boundry_layer(distance_threshold=distance_threshold_set, search_layer_setting_num=search_layer_setting_num)
					# end-del


					# increase scope of searched residues
					# INCREMENT - to activate multiloop uncomment next line and comment out four following lines
					# res_num_layer_detection_asserted += 1
					if res_num_layer_detection_asserted == 15:
						break
					else:
						res_num_layer_detection_asserted = 15
					# INCREMENT END

				# convert boundry_layers layers to closest-CA (this might work better with ap bundles)
				# boundry_layers = find_closest_CA_to_point(boundry_layers, helices_axis_all)

				# find all layers from designated layers
				# return list of helixAxisBundleClass objects truncated by defined first layers

				# print('===BEFORE===')
				# for bl in boundry_layers:
				# 	print(bl)
				# print('====='*15)

				# merge all boundry layers into one list - but only selected ones
				percentile_threshold = np.percentile([ bl.total_distance for bl in boundry_layers ], 50)
				boundry_layers       = [ bl for bl in boundry_layers if bl.total_distance < percentile_threshold ]
				#boundry_layers = heapq.nsmallest(res_num_layer_detection_asserted, boundry_layers, key=attrgetter('total_distance'))

				layers_sets = helices_axis_all.find_all_layers_from_layer(boundry_layers)
				# print('BOUNDRY LAYERS')
				# for bl in boundry_layers:
				# 	print(bl)
				# print('====='*15)
				#layers_sets = [layers_sets[0]]

				# first make filter for tot distance of layers based on given boundry layer
				# best_layer_set_distances  = select_minimal_total_distance_layer_set(layers_sets)
				# best_layer_set_dist2plane = select_minimal_dist_to_plane_layer_set(layers_sets)
				# for layer_set in best_layer_set_distances:
				# 	layer_set.show_points()
				# layers_sets = best_layer_set_distances
				# layers_sets = best_layer_set_dist2plane

				# select best layer set (for dimer and trimer+)
				# function should be keep outside of class

				# dimer case
				if len(input_helices) == 2:
					best_layer_set = select_minimal_distance_layer_set(layers_sets)
				# all other cases
				else:
					best_layer_set = select_minimal_angle_layer_set(layers_sets, best_layer_nb=1)
					# best_layer_set = select_minimal_angle_layer_set2(layers_sets, best_layer_nb=1)
					# check 3 params and get layer with min(rank, rank, rank) [devel-idea]
					# layers_sets = select_minimal_angle_layer_set(layers_sets, best_layer_nb='rank')
					# best_layer_set = select_minimal_total_distance_layer_set(layers_sets, best_layer_nb=1)
					# layers_sets = select_minimal_dist_to_plane_layer_set(layers_sets, best_layer_nb='rank')
					# best_layer_set = get_layers_set_with_best_ranks(layers_sets)

				#best_layer_set = select_minimal_total_distance_layer_set(layers_sets)[0]

				# for ls in layers_sets:
				# 	print(ls.average_dist_angle)
				# 	print(ls.average_dist_to_plane)

				# best_layer_set.show_points()


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

				# show selected set [dev]
				# for ha, ih in zip(best_layer_set[:].iterchain(), input_helices):
				# 	if ih[3]:
				# 		print([ha[-1].res.get_id()[1], ha[0].res.get_id()[1], ih[2], ih[3]])
				# 	else:
				# 		print([ha[0].res.get_id()[1], ha[-1].res.get_id()[1], ih[2], ih[3]])


				# generate pymol-readable selection
				# yellow layers between helices: all residues indicated by best_layer_set [dev-docs]
				# red helices coloring: all residues within reasonable distance from other helices (axisPoint.distance_flag=True) [dev-docs]
				# FIXME: blue helices coloring: all residues indictaed by Socket as CC but excluded because of distance_flag=False [dev-docs]
				# green: all other residues in a structure [dev-docs]

				# this will cut range of shown helices to the residues that are within resonable distance from all other helices
				# it does not affect range of drawn layers in pymol
				input_helices_cut = []
				for ha, ih in zip(helices_axis_all[:].iterchain(), input_helices):
					if ih[3]:
						input_helices_cut.append([ha[-1].res.get_id()[1], ha[0].res.get_id()[1], ih[2], ih[3]])
					else:
						input_helices_cut.append([ha[0].res.get_id()[1], ha[-1].res.get_id()[1], ih[2], ih[3]])
				if DEBUG:
					print('pymol session data', input_helices_cut)

				# helices_axis_all.show_points()
				# print('pymol session data', input_helices_cut)
				input_helices = input_helices_cut

				#FIXME pymol session in final version should indicate Socket detected region and show layer bars in actually measured layers (from bundleClass chains)
				#FIXME ppo probably not needed in production
				bundle.pymol_selection = create_pymol_selection_from_socket_results(input_helices)
				bundle.ppo = best_layer_set.ppo

				# save order of helices (helix indexes as in socket data) [ks mode only now]
				if mode == 'auto-detect':
					bundle.helix_order = best_layer_set.neighbour_interactions

				bundles.append(bundle)

		return bundles
