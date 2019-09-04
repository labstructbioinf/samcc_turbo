import itertools
import heapq
import numpy as np
import scipy.spatial.distance as distance
from Bio import PDB
from .bundle import get_local_axis

def create_pymol_selection_from_socket_results(indices):

	#FIXME format docstring

	''' create pymol-readable selection from socket data (converted to bundleDesc) '''
	''' input: list of helices, each helix represented as a tuple in format (first_residue_number, last_residue_number, chain_id) '''
	''' output: list of pymol-readable selections, each as a string '''

	selection = []

	for idx in indices:
		selection.append('(chain {} and resi {}-{})'.format(idx[2], idx[0], idx[1]))

	return selection

def detect_helices_orientation(indices, orientation_data):

	#FIXME check if not redundant with function from helper_functions
	#FIXME format docstring in function and subfunctions

	''' detect orientation of all helices from socket data '''
	''' input: list of helices, each as dict of their features, '''
	''' list of all helix orientations, each element one relation in form [helix1_id, helix2_id, orientation] '''
	''' output: list of helix orientation in given bundle, each element cooresponds to one helix (p/ap) '''

	def convert_orientation_data_to_dict(orientation_data):

		''' convert data about orientation of helices from list ot dictionary '''
		''' if there is missing data about relations - it is extrapolated from other relations '''
		''' input: list of all helix orientations, each element one relation in form [helix1_id, helix2_id, orientation]'''
		''' output: dictionary of orientation with tuples of (helix1_id, helix2_id) format as keys and orientations (str) as values '''

		def find_missing_relations(orientation_dict):

			# find missing relations
			hids              = set([ pair[0] for pair in orientation_dict.keys() ] + [ pair[1] for pair in orientation_dict.keys() ])
			all_relations     = set([ hid for hid in itertools.permutations(hids, 2) ])
			missing_relations = all_relations - set(orientation_dict.keys())
			for relation in missing_relations:
				for hid in hids:
					if (hid != relation[0] and hid != relation[1] and (hid, relation[0]) not in missing_relations and (hid, relation[1]) not in missing_relations):
						orientation_hid_to_rel0 = orientation_dict[(hid, relation[0])]
						orientation_hid_to_rel1 = orientation_dict[(hid, relation[1])]
						if (orientation_hid_to_rel0 == orientation_hid_to_rel1):
							orientation_dict[relation[0], relation[1]] = 'parallel'
						else:
							orientation_dict[relation[0], relation[1]] = 'antiparallel'
						break

			return missing_relations

		orientation_dict = {}
		for relation in orientation_data:
			orientation_dict[(relation[0], relation[1])] = relation[2]
			orientation_dict[(relation[1], relation[0])] = relation[2]

		missing_relations = find_missing_relations(orientation_dict)

		### some relations require additional iterations of algorithm
		### some relations are impossible to extrapolate but those are not necessarry for samcc (hence break possibility)_
		while(missing_relations):
			print('Searching for missing relations...')
			if missing_relations == find_missing_relations(orientation_dict):
				break
			else:
				missing_relations = find_missing_relations(orientation_dict)

		return orientation_dict

	def convert_indices_to_string(indices):
		''' convert int indices in coil description to string '''
		''' assures data compatibility '''

		return [ (str(coil[0]), str(coil[1]), str(coil[2]), coil[3]) for coil in indices ]

	# convert int to string to unify keys in further dict
	indices = convert_indices_to_string(indices)
	orientation      = [ '' for _s in indices]
	orientation_dict = convert_orientation_data_to_dict(orientation_data)
	helices_ids      = [ hid[0] for hid in indices ]

	# initialize first helix as parallel
	main_helix       = helices_ids[0]
	hid              = helices_ids.index(main_helix)
	orientation[hid] = 'p'
	check_later      = []

	# detect all other helices orientation in relation to the first one
	for other_helix in helices_ids[1:]:
		pair_orientation = orientation_dict[(main_helix, other_helix)]
		hid_other        = helices_ids.index(other_helix)
		if pair_orientation == 'parallel':
			orientation[hid_other] = 'p'
		elif pair_orientation == 'antiparallel':
			orientation[hid_other] = 'ap'

	return orientation

def get_dssp_data(model, pdbpath, selection, orientation, verbose=False, dssp='mkdssp'):

	#FIXME check if still needed - if yes: format docstring

	''' get dssp data for given selection of pdb file '''
	''' helices that are in antiparallel direction are reversed (residue order) '''
	''' returns list of helices, each helix is list of its residues in format (chain, ('', resnum, '')) '''
	''' additionally returns dict of format helix_id:pymol_selection_and_orientation '''

	dssp = PDB.DSSP(model, pdbpath, dssp=dssp)

	helices           = [] # each helix is list of dssp keys
	helices2selection = {} # each record is helix_id:selection orientation

	for helix in enumerate(selection):
		chain_id      = helix[1].strip('()').split(' ')[1]
		range_select  = helix[1].strip('()').split(' ')[-1]
		res_range     = range(int(range_select.split('-')[0]), int(range_select.split('-')[1])+1)
		current_helix = []

		for res in list(dssp.keys()):
			if (res[0] == chain_id and res[1][1] in res_range):
				current_helix.append(res)

		# if helix is in antiparallel orientation: reverse residue order
		if orientation[helix[0]] == 'ap':
			current_helix.reverse()
		helices.append(current_helix)
		helices2selection[helix[0]] = helix[1] + ' ' + orientation[helix[0]]

	if verbose:
		print('Detected "helices":')
		for helix in enumerate(helices):
			print('')
			print('HELIX ' + str(helix[0]) + ' length: ' + str(len(helix[1])))
			for res in helix[1]:
				print(dssp[res], res)

	return helices, helices2selection

def extract_CA_coords(helices, model, verbose=False):

	#FIXME format docstring

	''' get CA coordinates for rach residue in each helix '''
	''' input: helix list, each helix is list of its residues in format (chain, ('', resnum, '')) '''
	''' output: list of helices, each helix is list of CA coords (Biopython Vector class) of its residues '''

	print('Extracting CA coordinates...')
	helices_CA = []
	for helix in enumerate(helices):

		#FIXME
		###!!! TEMPORARY SOLUTION: ignore KeyError structures (i.e. those with HETATM)
		##
		try:
			helix_CA = [ model[res[0]][res[1]]['CA'].get_vector() for res in helix[1] ]
		except KeyError:
			return KeyError
		##
		###!!!

		if verbose:
			print('CA COORDS: HELIX ' + str(helix[0]))
			for h in helix_CA:
				print(h)
		helices_CA.append(helix_CA)
	print('Extracted.')

	return helices_CA

def convert_Vector_data_to_np_array(helices_CA_Vector):

	#FIXME format docstring

	''' converts list of helices format from Biopython Vector to lists of 3 points '''
	''' input: list of helices, each helix is list of Vectors, each is coords for one residue '''
	''' output: list of helices, each helix is list of 3-element lists describing residue coords '''

	return [ np.array([ (v[0], v[1], v[2]) for v in helix_CA ]) for helix_CA in helices_CA_Vector ]

def find_best_fit_line_to_helices_CAs(helices_CA, mode, res_num_layer_detection=5):

	#FIXME verify if still used - if yes: format docstring

	''' find line that will fit best to a set of helix CA atoms coords '''
	''' algorithm used is the algorithm used in SOCKET '''
	''' pts_per_res option is used only in SVD, pts_res_num is used in both, those params are explained in SVD/SOCKET fundtions '''

	''' input: list of helices, each helix is list of Vectors, each is coords for one residue '''
	''' input: mode of helix axis fit (currently only socket),  res_num_layer_detection is described in socket algorithm'''

	''' output: list of helix axis (HAP) points(start/stop), list of pts_res_num first HAP, list of pts_res_num last HAP '''
	''' output: described in detail in compute_helix_axis_and_layer_points_SOCKET() '''

	def compute_helix_axis_and_layer_points_SOCKET(helices_CA_array, res_num_ld):

		''' computes helix axis as a set of points, one per layer (layer = set of residues, one from each helix) '''
		''' input: list of helices, each helix is list of 3-element lists describing residue CA coords  '''
		''' input: res_num_l(ayer)d(detection) states how many residues from each end will be considered in search for first layer (int) '''
		''' output: (1)list of 2-element tuples (each tuple represents one helix), both elements of tuple are 3-element list of helix axis point coords '''
		''' first element represents first residue (helix axis) of the helix, second the last - connecting them with line in pymol will give rough estimate of helix axis '''
		''' output: (2)list of helices, each helix represented by list of 3-element coords of helix axis points, this is full, non truncated list of axis points '''
		''' output: (3/4)list of helices, each helix represented by list of 3-element coords of helix axis points (1 point/residue) '''
		''' both output 3 and 4 are truncated according to pts_res_num parameter '''
		''' if pts_res_num is set to n the output 3 will give only first n points on helix axis and output 4 will give last n points '''
		''' if any of helices is too short to get n layers then n is set to length of the shortest helix '''

		helices_axis      = []
		helices_axis_all  = []
		helices_pts_first = []
		helices_pts_last  = []

		for helix_CA in helices_CA:
			helix_vector     = get_local_axis(helix_CA)
			helix_list       = [ [v[0], v[1], v[2]] for v in helix_vector ]

			# get first and last residue coords for line drawing
			helices_axis.append((helix_list[0], helix_list[-1]))
			# store entire helix axis
			helices_axis_all.append(helix_list)

			# get points for layer search
			helices_pts_first.append(helix_list[:res_num_ld])
			helices_pts_last.append(helix_list[-res_num_ld:])

		return helices_axis, helices_axis_all, helices_pts_first, helices_pts_last

	helices_CA_array   = convert_Vector_data_to_np_array(helices_CA)
	shortest_helix_len = min([ len(h) for h in helices_CA_array ])
	if shortest_helix_len < res_num_layer_detection:
		res_num_layer_detection_asserted = shortest_helix_len - 2
		print('Warning: one of the helices was to short for start layer search with res_num_layer_detection= ' + str(res_num_layer_detection) + '; searching with res_num_layer_detection=' + str(res_num_layer_detection_asserted))
	else:
		res_num_layer_detection_asserted = res_num_layer_detection

	if mode == 'socket':
		return compute_helix_axis_and_layer_points_SOCKET(helices_CA_array, res_num_layer_detection_asserted)

def detect_helix_order(*helices_pts_set, helices_pts_all):

	#FIXME add/format docstrings, clean from dev code

	def filter_out_diagonals_angles(first_points_on_helix):

		def calc_plane_bis(x, y, z):
			a = np.column_stack((x, y, z))
			return np.linalg.lstsq(a, np.ones_like(x), rcond=None)[0]

		def project_points(x, y, z, a, b, c):
			"""
			Projects the points with coordinates x, y, z onto the plane
			defined by a*x + b*y + c*z = 1
			"""
			vector_norm = a*a + b*b + c*c
			normal_vector = np.array([a, b, c]) / np.sqrt(vector_norm)
			point_in_plane = np.array([a, b, c]) / vector_norm

			points = np.column_stack((x, y, z))
			points_from_point_in_plane = points - point_in_plane
			proj_onto_normal_vector = np.dot(points_from_point_in_plane,
											 normal_vector)
			proj_onto_plane = (points_from_point_in_plane -
							   proj_onto_normal_vector[:, None]*normal_vector)

			# print('POINT-IN-PLANE')
			# print(point_in_plane)
			# print('PROJ-ONTO-PLANE')
			# print(proj_onto_plane)

			return point_in_plane + proj_onto_plane

		def unit_vector(vector):
			""" Returns the unit vector of the vector.  """
			return vector / np.linalg.norm(vector)

		def angle_between(v1, v2):
			""" Returns the angle in radians between vectors 'v1' and 'v2'::

					>>> angle_between((1, 0, 0), (0, 1, 0))
					1.5707963267948966
					>>> angle_between((1, 0, 0), (1, 0, 0))
					0.0
					>>> angle_between((1, 0, 0), (-1, 0, 0))
					3.141592653589793
			"""
			v1_u = unit_vector(v1)
			v2_u = unit_vector(v2)
			return np.rad2deg(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

		def determine_helix_order(points_projected):

			helix_order          = [0] # set first helix (from socket data) as first helix
			current_helix_vector = points_projected[0,:] # get point representing first helix

			while(len(helix_order) < len(points_projected)):
				angles = {}
				for other_vector in enumerate(points_projected):
					if other_vector[0] not in helix_order:
						angles[other_vector[0]] = angle_between(current_helix_vector, other_vector[1])
				# print('ANGLES', angles)
				min_angle = min(angles, key=angles.get)
				helix_order.append(min_angle)
				current_helix_vector = points_projected[min_angle,:]

			return helix_order

		def make_list_of_neighbour_interactions(helix_order):

			neighbour_interactions = []

			for i in range(len(helix_order) - 1):
				neighbour_interactions.append((helix_order[i], helix_order[i+1]))
			neighbour_interactions.append((helix_order[0], helix_order[-1]))

			# correct order of helices (smaller id must be always first - this is format used later)
			neighbour_interactions = [ (i[0], i[1]) if i[0] < i[1] else (i[1], i[0]) for i in neighbour_interactions ]

			return neighbour_interactions

		# process points data
		x = [ point[0] for point in first_points_on_helix ]
		y = [ point[1] for point in first_points_on_helix ]
		z = [ point[2] for point in first_points_on_helix ]

		points_projected = project_points(x, y, z, *calc_plane_bis(x, y, z))

		ppo = project_points(x, y, z, *calc_plane_bis(x, y, z))

		# set center of coordinate system at centroid point
		points_projected = points_projected[:-1,:] - points_projected[-1,:]

		# print('POINTS PROJECTED CORRECTED TO ORIGIN')
		# print(points_projected)
        #
		# print('CDIST')
		points_distances    = distance.cdist(points_projected, points_projected)
		max_points_distance = points_distances[points_distances>0].max() - points_distances[points_distances>0].min()

		helix_order            = determine_helix_order(points_projected)
		# print('HELIX ORDER')
		# print(helix_order)
		neighbour_interactions = make_list_of_neighbour_interactions(helix_order)

		### test plotting - remove once working
		# import matplotlib.pyplot as plt
		# from mpl_toolkits.mplot3d import Axes3D
        #
		# # plot raw data
		# plt.figure()
		# ax = plt.subplot(111, projection='3d')
		# ax.scatter(points_projected[:,0], points_projected[:,1], points_projected[:,2], color='b')
		# plt.show()
		###

		return neighbour_interactions, max_points_distance, ppo

	def calculate_and_add_centroid(first_points_on_helix):

		average_x = np.mean([point[0] for point in first_points_on_helix])
		average_y = np.mean([point[1] for point in first_points_on_helix])
		average_z = np.mean([point[2] for point in first_points_on_helix])

		centroid = [average_x, average_y, average_z]
		first_points_on_helix.append(centroid)

		return first_points_on_helix

	def find_closest_first_points_on_helix_to_centroid(centroid, helices_pts_all):
		''' find a set of points so each point will be on different helix and distance from this point to centroid will be smallest for entire helix '''

		first_points_on_helix = []
		for helix in enumerate(helices_pts_all):
			dist = distance.cdist([centroid], helix[1])
			# print(dist[0])
			# print('CDIST-CENTROID', min(dist[0]), dist[0].argmin())
			first_points_on_helix.append(helices_pts_all[helix[0]][dist[0].argmin()])

		first_points_on_helix.append(centroid)
		return first_points_on_helix

	largest_points_distance     = []
	neighbour_interactions_list = []
	ppo = [] # points used to determine helix order, FIXME to better var name

	for helices_pts in helices_pts_set:
		pts_distances = {} # dict of distances between points on helices axes, key=(h_id1,h_id2)

		for h1 in enumerate(helices_pts):
			for h2 in enumerate(helices_pts):
				if (h1[0] != h2[0] and not (h2[0],h1[0]) in pts_distances.keys()):
					dist = distance.cdist(h1[1], h2[1]) # Compute distance between each pair of the two collections of inputs.
					pts_distances[(h1[0],h2[0])] = dist
		dist_matrices_number = len(pts_distances.keys()) # number of all possible distances (including diagonals)

		# mapping of distance id to record in pts_distances
		distance_to_pos = { pos:list(pts_distances.keys())[pos] for pos in range(dist_matrices_number)}
		# select only relations between neighbouring helices - new angle-based method
		# print('LEN H')
		# for h in helices_pts:
		# 	print(len(h))

		first_points_on_helix = [ helix[int(len(helix)/2)] for helix in helices_pts ]
		#print('FIRST POINT ON HELIX', first_points_on_helix)
		#first_points_on_helix = [ helix[0] for helix in helices_pts ]
		first_points_on_helix = calculate_and_add_centroid(first_points_on_helix)

		# find first_points_on_helix from centroid
		first_points_on_helix = find_closest_first_points_on_helix_to_centroid(first_points_on_helix[-1], helices_pts_all)

		neighbour_interactions, max_points_distance, ppo = filter_out_diagonals_angles(first_points_on_helix)
		ppo = first_points_on_helix
		# print(neighbour_interactions)
		# print(max_points_distance)
		largest_points_distance.append(max_points_distance)
		neighbour_interactions_list.append(neighbour_interactions)

	best_order = largest_points_distance.index(min(largest_points_distance))

	return neighbour_interactions_list[best_order], ppo

def find_bundle_boundry_layer_universal(helices_pts, neighbour_interactions, distance_threshold=50):

	#FIXME format docstrings, clean from devel code

	''' find bundle layer with minimal distance between helices (i.e. points on their axes) '''
	''' distance is calculated only between edges, diagonals are not included '''
	''' input: list of helices, each helix is list of n 3-element lists of coords, n=res_num_layer_detection '''
	''' input: distance_threshold(int): maximal sum of distances in layer that is accepted as properly detected layer '''
	''' output: layer with minimal distance between helices and list of helix interactions cosidered in calculations '''
	''' output format(first_layer): list of helices, each helix is one 3-element list with coords of point on its axis'''
	''' output format(neighbour_interactions):  list of 2-element tuples, first element is id of first helix, second element is id of second helix '''

	def get_distances(layer, distance_to_pos, pts_distances):

		''' calculate sum of distances between helices in layer '''
		''' input: tuple of ids of points on helices axis (tuple length = number of helices), '''
		''' dict mapping id of distance to tuple indicating helix interaction, '''
		''' dict of distances between point on helices axis '''
		''' output: sum of distances between helices in layer(float) '''

		# special case for dimers, only one distance per layer so immediate return
		if len(distance_to_pos) == 1:
			distance_id = distance_to_pos[0]
			first_pos = layer[distance_id[0]]
			second_pos = layer[distance_id[1]]
			return pts_distances[distance_id][first_pos][second_pos]

		# all other cases (3+ helices in bundle)
		distance_sum = 0
		for dist in enumerate(layer):
			# each layer is tuple of residue indexes, e.g. (0,0,1,1) means that first residues from helix 1 and 2 and second from 3 and 4 will be considered
			# dist is tuple of numerators and residue indexes, e.g. (0,1), (1,0), (2,1), (3,1)
			# number of helices is equal to number of considered interactions between them (diagonal exclusion) - therefore enumeration allows to consider all interactions
			distance_id = distance_to_pos[dist[0]]
			# distance_id is identifier of interaction, e.g. (0,1) is interaction between helix 0 and 1
			first_pos = layer[distance_id[0]]
			# this will give the id of residue in first of helices considered in interaction
			second_pos = layer[distance_id[1]]
			# same for residue from second interaction
			distance_sum += pts_distances[distance_id][first_pos][second_pos]
			# access calculated distance data between selected residues

		return distance_sum

	pts_distances = {} # dict of distances between points on helices axes, key=(h_id1,h_id2)

	for h1 in enumerate(helices_pts):
		for h2 in enumerate(helices_pts):
			if (h1[0] != h2[0] and not (h2[0],h1[0]) in pts_distances.keys()):
				dist = distance.cdist(h1[1], h2[1]) # Compute distance between each pair of the two collections of inputs.
				pts_distances[(h1[0],h2[0])] = dist
	dist_matrices_number = len(pts_distances.keys()) # number of all possible distances (including diagonals)

	# mapping of distance id to record in pts_distances
	distance_to_pos = { pos:list(pts_distances.keys())[pos] for pos in range(dist_matrices_number)}

	# select only relations between neighbouring helices
	distance_to_pos = { dist[0]:dist[1] for dist in distance_to_pos.items() if dist[1] in neighbour_interactions }
	distance_to_pos = { dist[0]:dist[1] for dist in enumerate(distance_to_pos.values()) }

	pts_number  = len(list(pts_distances.values())[0]) # number of layers considered (= res_num_layer_detection)
	points_list = [ p for p in range(pts_number) ] # list of layer ids (length = pts_number = res_num_layer_detection)
	# compare all distances and select the smallest one
	if len(distance_to_pos) == 1:
		# special case for dimers
		distances = { layer:get_distances(layer, distance_to_pos, pts_distances) for layer in itertools.product(points_list, repeat=2) }
	else:
		distances = { layer:get_distances(layer, distance_to_pos, pts_distances) for layer in itertools.product(points_list, repeat=len(distance_to_pos)) }
	min_distance = min(distances, key=distances.get)

	# extract coords of the layer that has minimal distance between helices
	first_layer = []
	for helix_point in enumerate(min_distance):
		first_layer.append([ pt for pt in helices_pts[helix_point[0]][helix_point[1]] ])

	# check if smallest distance between layers to large for reasonable layer
	if min(distances.values()) > distance_threshold:
		distance_check = False
	else:
		distance_check = True

	return first_layer, distance_check

def find_closest_CA_to_point(layer_pts, helices_CA):

	#FIXME format docstring
	#FIXME this should relate point on axis to residue that this point was calculated from

	''' find CA closest to the point on helix axis '''
	''' input: (layer_pts) list of helices, each represented as 3-element list of coords of point on their axis '''
	''' input: (helices_CA) list of helices, each represented as list of Biopython Vectors, each Vector represents one CA of helix '''
	''' output: (CA_coords) list of helices, each represented by 3-element list of coords of CA closest to desired point on helix axis '''
	''' output: (CA_ids) list of identifiers (int) of CAs that are closest to desired points, each identifier represents one helix '''

	helices_CA_array  = convert_Vector_data_to_np_array(helices_CA)
	CA_coords         = []
	CA_ids            = []

	for point in enumerate(layer_pts):
		dist       = distance.cdist(np.array([point[1]]), helices_CA_array[point[0]])
		closest_id = np.argmin(dist)

		closest_coords = [ coord for coord in helices_CA_array[point[0]][closest_id] ]
		CA_coords.append(closest_coords)
		CA_ids.append(closest_id)

	return CA_coords, CA_ids

def find_all_layers_from_layer(layer_ids, helices_CA):

	#FIXME format docstring and comments

	''' get all layers starting from selected layer (set of CAs) '''
	''' input: (layer_ids) first layer defined as list of identifiers (int) of CAs that comprise it, each identifier represents one helix '''
	''' input: (helices_CA) list of helices, each represented as list of Biopython Vectors, each Vector represents one CA of helix '''
	''' output: list of layers, each layer represented as n-element list of 3-element lists of coords (float) of CAs, n=number of helices '''

	helices_CA_array  = convert_Vector_data_to_np_array(helices_CA)
	layers_all        = []

	# check how long is shortest helix - this is stop downstream
	helices_length    = [ len(helix) for helix in helices_CA_array ]
	max_layers        = min(np.array(helices_length) - np.array(layer_ids))
	# check which id layer_id is smallest: layer_id-1 is start of layers "on top"
	closest_id_to_1   = min(layer_ids)
	max_layers        = max_layers + closest_id_to_1
	# adjust layer_ids so each element points to id of helix residue forming first layer upstream
	layer_ids         = list(map(lambda x: x - (closest_id_to_1), layer_ids))

	# take every layer starting from start to stop
	for layer in range(max_layers):
		layer_coords = []
		for helix in range(len(helices_CA_array)):
			layer_coords.append([ coord for coord in helices_CA_array[helix][layer+layer_ids[helix]] ])
		layers_all.append(layer_coords)

	return layers_all

def select_minimal_angle_layer_set_OLD(CA_layers_list, layer_ids_list):

	#FIXME archive/delete

	''' === this version of the function is deprecated === '''

	''' find layer set with minimal angle between layers '''
	''' input: (CA_layers_list) list of lists of layers, each layer represented as n-element list of 3-element lists of coords (float) of CAs '''
	''' n=number of helices, topmost list comprise of lists of layers generated by scanning with different methods (from top or bottom of the bundle) '''
	''' input: (layer_ids_list) list of lists of identifiers (int) of CAs used to define first layer, each identifier represents one helix '''
	''' topmost list elements correspond to the elements of topmost list from CA_layers_list '''
	''' output: (best_layer_set) list of layers with minimal angle between them in format identical to elements of topmost list of CA_layers_list '''
	''' output: (best_layer_ids) element of layer_ids_list corresponding to selected best_layer_set element'''

	def check_layers_shape(CA_layers):

		''' compute average angle between layers in bundle '''
		''' for each layer find all possible planes and compute angles between planes in vertical groups '''
		''' result is averaged angle from list of average angles from groups of planes '''
		''' input: list of layers, each layer is n-element list of 3-element lists (coords), n=oligomerization state '''
		''' output: average angle between layers (float) '''

		def calculate_plane_equation_from_points(layer):

			''' calculate plane equation in form ax+bx+cx+d=0 '''
			''' input: list of three points defining layer, each point represented as 3-element list of coords '''
			''' output: plane equation in form of dicionary where keys correspond to coefficients of equation '''

			a = (layer[1][1] - layer[0][1]) * (layer[2][2] - layer[0][2]) - (layer[2][1] - layer[0][1]) * (layer[1][2] - layer[0][2])
			b = (layer[1][2] - layer[0][2]) * (layer[2][0] - layer[0][0]) - (layer[2][2] - layer[0][2]) * (layer[1][0] - layer[0][0])
			c = (layer[1][0] - layer[0][0]) * (layer[2][1] - layer[0][1]) - (layer[2][0] - layer[0][0]) * (layer[1][1] - layer[0][1])
			d = -(a * layer[0][0] + b * layer[0][1] + c * layer[0][2])
			plane_equation = {'a':a, 'b':b, 'c':c, 'd':d}

			return plane_equation

		def calculate_angle_between_planes(plane1, plane2, silent_warning=True):

			''' calculate angle between two planes defined with equation ax+bx+cx+d '''
			''' cos = |a1a2 + b1b2 + c1c2| / (sqrt(a1^2 + b1^2 + c1^2) * sqrt(a2^2 + b2^2 + c2^2)) '''
			''' input: (plane1/2) plane equation in form of dicionary where keys correspond to coefficients of equation '''
			''' output: angle between planes in degrees (float) '''

			numerator  = abs(sum([ plane1[coef]*plane2[coef] for coef in 'abc' ]))
			len_plane1 = np.sqrt(sum([ np.power(plane1[coef], 2) for coef in 'abc' ]))
			len_plane2 = np.sqrt(sum([ np.power(plane2[coef], 2) for coef in 'abc' ]))

			cos_phi    = numerator / (len_plane1 * len_plane2)
			if (1+1e-8 > cos_phi > 1.0):
				cos_phi = 1.0
				if not silent_warning:
					print('Warning: calculated cosine slightly beyound function codomain(max 1e-8 beyound) - corrected to 1.0.')
			phi_deg    = np.rad2deg(np.arccos(cos_phi))

			return phi_deg

		angle_layers_frame = []

		for layer_frame in itertools.combinations([ l for l in range(len(CA_layers[0]))], 3):

			#print('FRAME', layer_frame)
			plane_layers = []
			angle_layers = []

			for layer in CA_layers:
				layer_from_frame = [layer[layer_frame[0]], layer[layer_frame[1]], layer[layer_frame[2]]]
				plane_layers.append(calculate_plane_equation_from_points(layer_from_frame))

			for layer in range(len(plane_layers)-1):
				angle_layers.append(calculate_angle_between_planes(plane_layers[layer], plane_layers[layer+1]))

			angle_layers_frame.append(np.mean(angle_layers))

		return np.mean(angle_layers_frame)

	layer_set_angles     = (list(map(check_layers_shape, CA_layers_list)))
	best_layer_set_angle = min(layer_set_angles)

	best_layer_set = CA_layers_list[layer_set_angles.index(best_layer_set_angle)]
	best_layer_ids = layer_ids_list[layer_set_angles.index(best_layer_set_angle)]

	return best_layer_set, best_layer_ids

def select_minimal_angle_layer_set(CA_layers_list, layer_ids_list):

	#FIXME format docstring, clean from devel code

	''' find layer set with minimal angle between layers '''
	''' input: (CA_layers_list) list of lists of layers, each layer represented as n-element list of 3-element lists of coords (float) of CAs '''
	''' n=number of helices, topmost list comprise of lists of layers generated by scanning with different methods (from top or bottom of the bundle) '''
	''' input: (layer_ids_list) list of lists of identifiers (int) of CAs used to define first layer, each identifier represents one helix '''
	''' topmost list elements correspond to the elements of topmost list from CA_layers_list '''
	''' output: (best_layer_set) list of layers with minimal angle between them in format identical to elements of topmost list of CA_layers_list '''
	''' output: (best_layer_ids) element of layer_ids_list corresponding to selected best_layer_set element'''

	def check_layers_shape(CA_layers):

		''' compute average angle between layers in bundle '''
		''' for each layer find all possible planes and compute angles between planes in vertical groups '''
		''' result is averaged angle from list of average angles from groups of planes '''
		''' input: list of layers, each layer is n-element list of 3-element lists (coords), n=oligomerization state '''
		''' output: average angle between layers (float) '''

		def calculate_plane_equation_from_points(x, y, z):
			a = np.column_stack((x, y, z))
			return np.linalg.lstsq(a, np.ones_like(x), rcond=None)[0]

		def calculate_angle_between_planes(plane1, plane2, silent_warning=True):

			''' calculate angle between two planes defined with equation ax+bx+cx+d '''
			''' cos = |a1a2 + b1b2 + c1c2| / (sqrt(a1^2 + b1^2 + c1^2) * sqrt(a2^2 + b2^2 + c2^2)) '''
			''' input: (plane1/2) plane equation in form of dicionary where keys correspond to coefficients of equation '''
			''' output: angle between planes in degrees (float) '''

			numerator  = abs(sum([ plane1[coef]*plane2[coef] for coef in 'abc' ]))
			len_plane1 = np.sqrt(sum([ np.power(plane1[coef], 2) for coef in 'abc' ]))
			len_plane2 = np.sqrt(sum([ np.power(plane2[coef], 2) for coef in 'abc' ]))

			cos_phi    = numerator / (len_plane1 * len_plane2)
			if (1+1e-8 > cos_phi > 1.0):
				cos_phi = 1.0
				if not silent_warning:
					print('Warning: calculated cosine slightly beyound function codomain(max 1e-8 beyound) - corrected to 1.0.')
			phi_deg    = np.rad2deg(np.arccos(cos_phi))

			return phi_deg

		layer_equations = []
		layer_angles    = []

		for layer in CA_layers:
			# process points data
			x = [ point[0] for point in layer ]
			y = [ point[1] for point in layer ]
			z = [ point[2] for point in layer ]

			plane_equation = calculate_plane_equation_from_points(x,y,z)
			layer_equations.append({'a':plane_equation[0], 'b':plane_equation[1], 'c':plane_equation[2]})

		for layer in range(len(layer_equations)-1):
			# calculate angle between two planes
			layer_angles.append(calculate_angle_between_planes(layer_equations[layer], layer_equations[layer+1]))


		return np.mean(layer_angles)

	layer_set_angles     = (list(map(check_layers_shape, CA_layers_list)))
	# print(layer_ids_list)
	# print(layer_set_angles)
	best_layer_set_angle = min(layer_set_angles)

	best_layer_set = CA_layers_list[layer_set_angles.index(best_layer_set_angle)]
	best_layer_ids = layer_ids_list[layer_set_angles.index(best_layer_set_angle)]

	return best_layer_set, best_layer_ids

def select_minimal_distance_layer_set(CA_layers_list, layer_ids_list):

	#FIXME rewrite to use new idea for dimers: cast on plane and measure angles between lines connecting helices
	#FIXME format docstring

	''' find set of layers with minimal distance between 2 helices (dimers only) '''
	''' input/output same as described in select_minimal_angle_layer_set function '''
	''' only difference is that this function search for minimal distance between helices, not angle between layers '''

	def check_layers_distances(CA_layers):

		layer_distances = []

		for layer in CA_layers:
			point1 = layer[0]
			point2 = layer[1]

			layer_distances.append(distance.euclidean(point1, point2))

		return np.mean(layer_distances)

	layer_set_distances     = (list(map(check_layers_distances, CA_layers_list)))
	best_layer_set_distance = min(layer_set_distances)

	best_layer_set = CA_layers_list[layer_set_distances.index(best_layer_set_distance)]
	best_layer_ids = layer_ids_list[layer_set_distances.index(best_layer_set_distance)]

	return best_layer_set, best_layer_ids

def map_layers_to_resi(best_layer_ids, selection, orientation, bundle_length):

	#FIXME verify if still needed - if yes: format docstring

	''' convert detected layer definition to SamCC-readable selection '''
	''' input: (best_layer_ids) list of identifiers (int) of CAs used to define first layer, each identifier represents one helix '''
	''' input: (selection) list of pymol-readable selections, each as a string '''
	''' input: (orientation) list of helix orientation in given bundle, each element cooresponds to one helix (p/ap) '''
	''' input: (bundle_length) length of the bundle in current layer definition (int) '''
	''' output: samcc-ready selection, list of helices, each helix a list of format [ chain_id(string), start_residue(int), stop_residue(int) ] '''

	samcc_selection = []

	# find start,end residue for each helix in selection
	for helix in enumerate(selection):
		res_range  = helix[1].strip('()').split('resi')[1].split('-')
		res_start  = int(res_range[0])
		res_end    = int(res_range[1])
		chain_id   = helix[1].strip('()').split(' ')[1]

		len_upstream = min(best_layer_ids[0])
		if orientation[helix[0]] == 'p':
			samcc_start = res_start + best_layer_ids[0][helix[0]] - len_upstream
			samcc_end   = samcc_start + bundle_length - 1
			samcc_data  = [chain_id, samcc_start, samcc_end]
		elif orientation[helix[0]] == 'ap':
			samcc_start = res_end - best_layer_ids[0][helix[0]] + len_upstream
			samcc_end   = samcc_start - bundle_length + 1
			samcc_data  = [chain_id, samcc_start, samcc_end]
		samcc_selection.append(samcc_data)

	return samcc_selection

def convert_to_bundleClass_format(samcc_selection):

	#FIXME format docstring

	''' converts samcc selection from this library to original format used by bundleClass '''
	''' input: samcc-ready selection, list of helices, each helix a list of format [ chain_id(string), start_residue(int), stop_residue(int) ] '''
	''' output: BundleClass-ready selection, list of helices, each helix a tuple of format ( start_residue(int), stop_residue(int), chain_id(string), antiparallel(bool) ) '''

	input_helices = []
	for helix in samcc_selection:
		if helix[1] < helix[2]:
			input_helices.append((helix[1], helix[2], helix[0], False))
		else:
			input_helices.append((helix[2], helix[1], helix[0], True))

	return input_helices
