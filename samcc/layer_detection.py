### Main layer detection engine ###

import itertools
import heapq
import numpy as np
import scipy.spatial.distance as distance
import scipy.optimize
import functools
from operator import attrgetter
from Bio import PDB
from .bundle import get_local_axis

def create_pymol_selection_from_socket_results(indices):
	"""create pymol-readable selection from socket data (converted to bundleDesc)
	input: list of helices, each helix represented as a tuple in format (first_residue_number, last_residue_number, chain_id)
	output: list of pymol-readable selections, each as a string
	"""

	selection = []

	for idx in indices:
		selection.append('(chain {} and resi {}-{})'.format(idx[2], idx[0], idx[1]))

	return selection

def extract_CA_coords(helices, model, verbose=False):
	"""get CA coordinates for rach residue in each helix
	input: helix list, each helix is list of its residues in format (chain, ('', resnum, ''))
	output: list of helices, each helix is list of CA coords (Biopython Vector class) of its residues
	"""

	print('Extracting CA coordinates...')
	helices_CA = []
	for helix in enumerate(helices):

		#ignore KeyError structures (i.e. those with HETATM)
		try:
			helix_CA = [ model[res[0]][res[1]]['CA'].get_vector() for res in helix[1] ]
		except KeyError:
			return KeyError

		if verbose:
			print('CA COORDS: HELIX ' + str(helix[0]))
			for h in helix_CA:
				print(h)
		helices_CA.append(helix_CA)

	return helices_CA

def convert_Vector_data_to_np_array(helices_CA_Vector):
	"""converts list of helices format from Biopython Vector to lists of 3 points
	input: list of helices, each helix is list of Vectors, each is coords for one residue
	output: list of helices, each helix is list of 3-element lists describing residue coords
	"""

	return [ np.array([ (v[0], v[1], v[2]) for v in helix_CA ]) for helix_CA in helices_CA_Vector ]

def select_minimal_angle_layer_set2(layers_sets, best_layer_nb=1):

	def check_layers_shape(layers_set):

		def calculate_plane_equation_from_points(points):

			def plane(x, y, params):
				a = params[0]
				b = params[1]
				c = params[2]
				z = a*x + b*y + c
				return z

			def error(params, points):
				result = 0
				for (x,y,z) in points:
					plane_z = plane(x, y, params)
					diff = abs(plane_z - z)
					result += diff**2
				return result

			def cross(a, b):
				return [a[1]*b[2] - a[2]*b[1],
						a[2]*b[0] - a[0]*b[2],
						a[0]*b[1] - a[1]*b[0]]

			fun = functools.partial(error, points=points)
			params0 = [0, 0, 0]
			res = scipy.optimize.minimize(fun, params0)

			return res

		def calculate_angle_between_planes(plane1, plane2, silent_warning=True):

			"""calculate angle between two planes defined with equation ax+bx+cx+d
			cos = |a1a2 + b1b2 + c1c2| / (sqrt(a1^2 + b1^2 + c1^2) * sqrt(a2^2 + b2^2 + c2^2))
			input: (plane1/2) plane equation in form of dicionary where keys correspond to coefficients of equation
			output: angle between planes in degrees (float)
			"""

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

		for layer in layers_set.iterlayer():

			points = [ tuple(point.coords) for point in layer ]

			plane_equation = calculate_plane_equation_from_points(points)
			layer_equations.append({'a':plane_equation.x[0], 'b':plane_equation.x[1], 'c':plane_equation.x[2]})

		for layer in range(len(layer_equations)-1):
			# calculate angle between two planes
			layer_angles.append(calculate_angle_between_planes(layer_equations[layer], layer_equations[layer+1]))

		layers_set.average_dist_angle = np.median(layer_angles) #FIXME np.mean(layer_angles)
		print('V2', layers_set.average_dist_angle)

		return layers_set

	# this will return same layer sets list but with layers with set avg angle attribute
	layer_set_angles = (list(map(check_layers_shape, layers_sets)))

	if best_layer_nb == 1:
		best_layer_set_angle = min(layer_set_angles, key=attrgetter('average_dist_angle'))
	elif best_layer_nb == 'rank':
		best_layer_set_angle = sorted(layer_set_angles, key=attrgetter('average_dist_angle'))
	else:
		best_layer_set_angle = heapq.nsmallest(best_layer_nb, layer_set_angles, key=attrgetter('average_dist_angle'))

	return best_layer_set_angle

def select_minimal_angle_layer_set(layers_sets, best_layer_nb=1):
	"""find layer set with minimal angle between layers
	input: (CA_layers_list) list of lists of layers, each layer represented as n-element list of 3-element lists of coords (float) of CAs
	n=number of helices, topmost list comprise of lists of layers generated by scanning with different methods (from top or bottom of the bundle)
	input: (layer_ids_list) list of lists of identifiers (int) of CAs used to define first layer, each identifier represents one helix
	topmost list elements correspond to the elements of topmost list from CA_layers_list
	output: (best_layer_set) list of layers with minimal angle between them in format identical to elements of topmost list of CA_layers_list
	output: (best_layer_ids) element of layer_ids_list corresponding to selected best_layer_set element
	"""

	def check_layers_shape(layers_set):
		"""compute average angle between layers in bundle
		for each layer find all possible planes and compute angles between planes in vertical groups
		result is averaged angle from list of average angles from groups of planes
		input: list of layers, each layer is n-element list of 3-element lists (coords), n=oligomerization state
		output: average angle between layers (float)
		"""

		def calculate_plane_equation_from_points(x, y, z):
			a = np.column_stack((x, y, z))
			return np.linalg.lstsq(a, np.ones_like(x), rcond=None)[0]

		def calculate_angle_between_planes(plane1, plane2, silent_warning=True):

			"""calculate angle between two planes defined with equation ax+bx+cx+d
			cos = |a1a2 + b1b2 + c1c2| / (sqrt(a1^2 + b1^2 + c1^2) * sqrt(a2^2 + b2^2 + c2^2))
			input: (plane1/2) plane equation in form of dicionary where keys correspond to coefficients of equation
			output: angle between planes in degrees (float)
			"""

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

		for layer in layers_set.iterlayer():

			x = [ point.coords[0] for point in layer ] # all x coords from layer
			y = [ point.coords[1] for point in layer ]
			z = [ point.coords[2] for point in layer ]

			plane_equation = calculate_plane_equation_from_points(x,y,z)
			layer_equations.append({'a':plane_equation[0], 'b':plane_equation[1], 'c':plane_equation[2]})

		for layer in range(len(layer_equations)-1):
			# calculate angle between two planes
			layer_angles.append(calculate_angle_between_planes(layer_equations[layer], layer_equations[layer+1]))

		layers_set.average_dist_angle = np.median(layer_angles)

		return layers_set

	# this will return same layer sets list but with layers with set avg angle attribute
	layer_set_angles = (list(map(check_layers_shape, layers_sets)))

	if best_layer_nb == 1:
		best_layer_set_angle = min(layer_set_angles, key=attrgetter('average_dist_angle'))
	elif best_layer_nb == 'rank':
		best_layer_set_angle = sorted(layer_set_angles, key=attrgetter('average_dist_angle'))
	else:
		best_layer_set_angle = heapq.nsmallest(best_layer_nb, layer_set_angles, key=attrgetter('average_dist_angle'))

	return best_layer_set_angle

def find_bundle_boundry_layer_from_all(boundry_layers, distance_threshold, search_layer_setting_num):

	min_distance_set = heapq.nsmallest(search_layer_setting_num, boundry_layers)

	return [ layer for layer in min_distance_set if layer.total_distance <= distance_threshold ]

def select_minimal_distance_layer_set(layers_sets):
	"""find set of layers with minimal distance between 2 helices (dimers only)
	input/output same as described in select_minimal_angle_layer_set function
	only difference is that this function search for minimal distance between helices, not angle between layers
	"""

	def check_layers_distances(layers_set):

		layer_distances = []

		for layer in layers_set.iterlayer():
			point1 = layer[0].coords
			point2 = layer[1].coords

			layer_distances.append(distance.euclidean(point1, point2))

		layers_set.average_dist_angle = np.mean(layer_distances)

		return layers_set

	layer_set_distances     = (list(map(check_layers_distances, layers_sets)))
	best_layer_set_distance = min(layer_set_distances, key=attrgetter('average_dist_angle'))

	return best_layer_set_distance

def convert_to_bundleClass_format(samcc_selection):
	"""converts samcc selection from this library to original format used by bundleClass
	input: samcc-ready selection, list of helices, each helix a list of format [ chain_id(string), start_residue(int), stop_residue(int) ]
	output: BundleClass-ready selection, list of helices, each helix a tuple of format ( start_residue(int), stop_residue(int), chain_id(string), antiparallel(bool) )
	"""

	input_helices = []
	for helix in samcc_selection:
		if helix[1] < helix[2]:
			input_helices.append((helix[1], helix[2], helix[0], False))
		else:
			input_helices.append((helix[2], helix[1], helix[0], True))

	return input_helices

def select_minimal_dist_to_plane_layer_set(layers_sets, best_layer_nb=1):

	def check_layers_dist_to_plane(layers_set):

		def calculate_plane_equation_from_points(x, y, z):
			a = np.column_stack((x, y, z))
			return np.linalg.lstsq(a, np.ones_like(x), rcond=None)[0]

		def calculate_point_to_plane_distance(plane_equation, point):

			# plane coefficients
			A = plane_equation[0]
			B = plane_equation[1]
			C = plane_equation[2]

			# point coords
			x = point.coords[0]
			y = point.coords[1]
			z = point.coords[2]

			numerator   = abs(A*x + B*y + C*z)
			denominator = np.sqrt(A**2 + B**2 + C**2)

			return numerator / denominator

		layer_distances = []

		for layer in layers_set.iterlayer():

			x = [ point.coords[0] for point in layer ] # all x coords from layer
			y = [ point.coords[1] for point in layer ]
			z = [ point.coords[2] for point in layer ]

			plane_equation  = calculate_plane_equation_from_points(x,y,z)
			points_in_layer = []
			for point in layer:
				points_in_layer.append(calculate_point_to_plane_distance(plane_equation, point))

			layer_distances.append(np.mean(points_in_layer))

		layers_set.average_dist_to_plane = np.mean(layer_distances)

		return layers_set

	layer_set_dist_to_plane = (list(map(check_layers_dist_to_plane, layers_sets)))

	if best_layer_nb == 1:
		best_layer_set_dist_to_plane = min(layer_set_dist_to_plane, key=attrgetter('average_dist_to_plane'))
	elif best_layer_nb == 'rank':
		best_layer_set_dist_to_plane = sorted(layer_set_dist_to_plane, key=attrgetter('average_dist_to_plane'))
	else:
		best_layer_set_dist_to_plane = heapq.nsmallest(best_layer_nb, layer_set_dist_to_plane, key=attrgetter('average_dist_to_plane'))

	return best_layer_set_dist_to_plane

def select_minimal_total_distance_layer_set(layers_sets, best_layer_nb=1):

	def check_layers_total_distances(layers_set):

		def calculate_total_distance(layer, neighbour_interactions):
			"""Calculate total distance between neighbouring points in layer
			neighbourhood is determined according to neighbour_interactions list
			"""

			total_distance = 0

			for axis_point1 in layer:
				for axis_point2 in layer:
					if (axis_point1.helix_id, axis_point2.helix_id) in neighbour_interactions:
						total_distance += distance.euclidean(axis_point1.CA_coords, axis_point2.CA_coords)

			return total_distance

		layers_all_distances = []
		for layer in layers_set.iterlayer():
			layers_all_distances.append(calculate_total_distance(layer, layers_set.neighbour_interactions))

		layers_set.average_layers_dist = np.mean(layers_all_distances)

		return layers_set

	layer_set_total_distances = list(map(check_layers_total_distances, layers_sets))

	if best_layer_nb == 1:
		best_layer_set_total_distances = min(layer_set_total_distances, key=attrgetter('average_layers_dist'))
	elif best_layer_nb == 'rank':
		best_layer_set_total_distances = sorted(layer_set_total_distances, key=attrgetter('average_layers_dist'))
	else:
		best_layer_set_total_distances = heapq.nsmallest(best_layer_nb, layer_set_total_distances, key=attrgetter('average_layers_dist'))

	return best_layer_set_total_distances

def get_layers_set_with_best_ranks(layers_sets):

	best_layer_set_angle = sorted(layers_sets, key=attrgetter('average_dist_angle'))
	best_layer_set_dist_to_plane = sorted(layers_sets, key=attrgetter('average_dist_to_plane'))
	best_layer_set_total_distances = sorted(layers_sets, key=attrgetter('average_layers_dist'))

	for pos, layer_set in enumerate(best_layer_set_angle):
		layer_set.ranks += pos
	for pos, layer_set in enumerate(best_layer_set_total_distances):
		layer_set.ranks += pos

	return min(layers_sets, key=attrgetter('ranks'))

def find_closest_CA_to_point(boundry_layers, helices_axis_all):

	def get_closest_CA_to_axis_point(point, helix_pts):

		helix_point = helix_pts[0]
		dst_best    = distance.euclidean(helix_pts[0].CA_coords, point.coords)
		for helix_pt in helix_pts:
			dst = distance.euclidean(helix_pt.CA_coords, point.coords)
			if dst < dst_best:
				helix_point = helix_pt
				dst_best    = dst
		return helix_point

	boundry_layers_CA = []

	for layer, helix_axis in zip(boundry_layers, helices_axis_all):
		layer_adjusted = []
		for point in layer:
			closest_CA_helix_point = get_closest_CA_to_axis_point(point, [helix_axis[point.point_id+i] for i in range(-2,3) if point.point_id+i>0])
			layer_adjusted.append(closest_CA_helix_point)
		layer.axis_points = layer_adjusted
		boundry_layers_CA.append(layer)

	return boundry_layers_CA
