### classes for storing bundle as a collection of helix axis and points on helix axis ###

import itertools
import heapq
import numpy as np
from scipy.spatial import distance
from .chainClass import chainClass

class NotEqualAxisLen(Exception):
	pass

class helixAxisBundleClass():
	### class for storing all helices axis

	def __init__(self, chains, mode='from_chain', neighbour_interactions=None, ppo=None, average_dist_angle=None,
				 average_layers_dist=None, average_dist_to_plane=None, ranks=0):
		if mode == 'from_chain':
			self.helix_axis_all = [ helixAxisClass(helix_axis, helix_id=helix_id) for helix_id, helix_axis in enumerate(chains) ]
			self.neighbour_interactions = neighbour_interactions
			self.ppo = ppo
			self.average_dist_angle = average_dist_angle
			self.average_layers_dist = average_layers_dist
			self.average_dist_to_plane = average_dist_to_plane
			self.ranks = ranks
		elif mode == 'from_slice':
			self.helix_axis_all = chains
			self.neighbour_interactions = neighbour_interactions
			self.ppo = ppo
			self.average_dist_angle = average_dist_angle
			self.average_layers_dist = average_layers_dist
			self.average_dist_to_plane = average_dist_to_plane
			self.ranks = ranks

	def __getitem__(self, index):
		"""return new helixAxisBundleClass object but truncated
		when slicing only points with distance_flag = True are considered
		"""

		if isinstance(index, slice):
			helix_axis_sliced = []
			for helix_axis in self.helix_axis_all:
				# select only point with distance_flag = True and slice according to index
				helix_axis_sliced.append(helixAxisClass([ p for p in helix_axis.axis_points if p.distance_flag ][index]))

			return helixAxisBundleClass(helix_axis_sliced, mode='from_slice', neighbour_interactions=self.neighbour_interactions,
										ppo=self.ppo, average_dist_angle=self.average_dist_angle, average_layers_dist=self.average_layers_dist,
										average_dist_to_plane=self.average_dist_to_plane, ranks=self.ranks)
		else:
			return self.helix_axis_all[index]

	def iterchain(self):
		return self.helix_axis_all

	def iterlayer(self):
		"""Iterate over layers of the bundle, each item is one layer (list of points)
		Works only if all axis of equal length, else raise error
		"""

		# check if all helix axis of equal len, else raise error
		if not all(i.length() == self.helix_axis_all[0].length() for i in self.helix_axis_all):
			raise NotEqualAxisLen('Not all helix axis in the bundle are of equal length')

		# construct list of layers
		layers_all = []
		for l in range(self.helix_axis_all[0].length()):
			layer = []
			for ha in range(len(self.helix_axis_all)):
				layer.append(self.helix_axis_all[ha][l])
			layers_all.append(layer)

		return layers_all

	def get_middle_points(self, slice_size):
		"""for each helix take n/n-1 points from middle of a given helix; n when n is even and n-1 when n is uneven"""

		helix_axis_sliced = []
		for helix_axis in self.helix_axis_all:
			# select only point with distance_flag = True and slice according to index
			points_flag  = [ p for p in helix_axis.axis_points if p.distance_flag ]
			middle_point = int(len(points_flag) / 2)
			helix_axis_sliced.append(helixAxisClass(points_flag[middle_point - int(slice_size/2):middle_point + int(slice_size/2)]))

		return helixAxisBundleClass(helix_axis_sliced, mode='from_slice', neighbour_interactions=self.neighbour_interactions,
									ppo=self.ppo, average_dist_angle=self.average_dist_angle, average_layers_dist=self.average_layers_dist,
									average_dist_to_plane=self.average_dist_to_plane, ranks=self.ranks)

	def convert_to_coords_list(self):
		"""Return list of lists[helix axis] of lists[point coords]"""

		helix_axis_convert = []
		for helix_axis in self.helix_axis_all:
			helix_axis_convert.append([ [ c for c in p.coords ] for p in helix_axis.axis_points if (p.distance_flag and p.coords != []) ])

		return helix_axis_convert

	def detect_helix_order(self, res_num_layer_detection_asserted):

		def filter_out_diagonals_angles(first_points_on_helix):

			def calc_plane_bis(x, y, z):
				a = np.column_stack((x, y, z))
				return np.linalg.lstsq(a, np.ones_like(x), rcond=None)[0]

			def project_points(x, y, z, a, b, c):
				"""Projects the points with coordinates x, y, z onto the plane
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

			points_distances    = distance.cdist(points_projected, points_projected)
			max_points_distance = points_distances[points_distances>0].max() - points_distances[points_distances>0].min()

			helix_order            = determine_helix_order(points_projected)
			neighbour_interactions = make_list_of_neighbour_interactions(helix_order)

			return neighbour_interactions, max_points_distance, ppo

		def calculate_and_add_centroid(first_points_on_helix):

			average_x = np.mean([point[0] for point in first_points_on_helix])
			average_y = np.mean([point[1] for point in first_points_on_helix])
			average_z = np.mean([point[2] for point in first_points_on_helix])

			centroid = [average_x, average_y, average_z]
			first_points_on_helix.append(centroid)

			return first_points_on_helix

		def find_closest_first_points_on_helix_to_centroid(centroid, helices_pts_all):
			"""find a set of points so each point will be on different helix and distance
			from this point to centroid will be smallest for entire helix
			"""

			first_points_on_helix = []
			for helix in enumerate(helices_pts_all):
				dist = distance.cdist([centroid], helix[1])
				first_points_on_helix.append(helices_pts_all[helix[0]][dist[0].argmin()])

			first_points_on_helix.append(centroid)
			return first_points_on_helix

		middle_pts = self.get_middle_points(res_num_layer_detection_asserted)

		helices_pts_set = (middle_pts.convert_to_coords_list(),)
		helices_pts_all = self[:].convert_to_coords_list()

		largest_points_distance     = []
		neighbour_interactions_list = []
		ppo = [] # points used to determine helix order

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

			first_points_on_helix = [ helix[int(len(helix)/2)] for helix in helices_pts ]
			first_points_on_helix = calculate_and_add_centroid(first_points_on_helix)

			# find first_points_on_helix from centroid
			first_points_on_helix = find_closest_first_points_on_helix_to_centroid(first_points_on_helix[-1], helices_pts_all)

			neighbour_interactions, max_points_distance, ppo = filter_out_diagonals_angles(first_points_on_helix)
			ppo = first_points_on_helix
			largest_points_distance.append(max_points_distance)
			neighbour_interactions_list.append(neighbour_interactions)

		best_order = largest_points_distance.index(min(largest_points_distance))

		# add to self:
		self.neighbour_interactions = neighbour_interactions_list[best_order]
		self.ppo = ppo

	def get_all_bundle_boundry_layers(self):
		"""Only get list of searchLayers - do not evaluate best ones at this point"""

		return [ searchLayer(l, self.neighbour_interactions) for l in itertools.product(*self.helix_axis_all, repeat=1) ]

	def find_bundle_boundry_layer(self, distance_threshold, search_layer_setting_num):

		# use searchLayer class to create layer settings and calculate distance between their points
		layers_all = [ searchLayer(l, self.neighbour_interactions) for l in itertools.product(*self.helix_axis_all, repeat=1) ]

		# get set of layers with n smallest distances between their points
		min_distance_set = heapq.nsmallest(search_layer_setting_num, layers_all)

		# check if smallest distance between layers to large for reasonable layer
		# return only layers that meet reasonable distance criteria
		return [ layer for layer in min_distance_set if layer.total_distance <= distance_threshold ]

	def find_all_layers_from_layer(self, boundry_layers):

		layer_sets = []

		for layer in boundry_layers:
			helix_data = []
			for helix_axis, layer_point in zip(self.helix_axis_all, layer):
				# search helix axis for boundry point and at the same time calculate how many points upstream and downstream from it
				pts_before = 0
				pts_after  = 0
				pt_found   = False
				for axis_point in helix_axis:

					if axis_point.distance_flag == False:
						continue
					if axis_point.point_id == layer_point.point_id:
						pt_found = True
					# point before investigated one
					elif pt_found == False:
						pts_before += 1
					elif pt_found == True:
						pts_after += 1
				helix_data.append([pts_before, pts_after])
			max_upstream   = min([h[0] for h in helix_data])
			max_downstream = min([h[1] for h in helix_data])

			# extract axis points and construct helixAxisBundleClass:
			helix_axis_sliced = []
			for helix_axis, layer_point in zip(self.helix_axis_all, layer):
				idx_start = layer_point.point_id - max_upstream
				idx_stop  = layer_point.point_id + max_downstream + 1
				helix_axis_sliced.append(helix_axis[idx_start:idx_stop])
			layer_sets.append(helixAxisBundleClass(helix_axis_sliced, mode='from_slice', neighbour_interactions=self.neighbour_interactions,
												   ppo=self.ppo, average_dist_angle=self.average_dist_angle, average_layers_dist=self.average_layers_dist,
												   average_dist_to_plane=self.average_dist_to_plane, ranks=self.ranks))

		return layer_sets

	def verify_points_distances(self, max_distance=20):
		"""flag axis points that does not meet distance criteria
		code verifying if all helix axis points have any axis points
		from other helix within reasonable distance (20A) [dev-doc]
		this code compares all vs all axis points from all helices
		"""

		for h1_idx, helix_axis1 in enumerate(self.helix_axis_all):
			for point1 in helix_axis1.axis_points:
				# if edge (first or last on helix) flag as False
				if point1.coords == []:
					point1.distance_flag = False
					continue
				point_stat = []
				for h2_idx, helix_axis2 in enumerate(self.helix_axis_all):
					if h1_idx != h2_idx:
						for point2 in helix_axis2.axis_points:
							if (any(point1.CA_coords) and any(point2.CA_coords)):
								dst = distance.euclidean(point1.CA_coords, point2.CA_coords)
								if dst < max_distance:
									# found residue within cutoff distance on this helix, break and search next helix
									point_stat.append(True)
									break
						else:
							# there is no residue within cutoff distance on this helix, search next one
							point_stat.append(False)

				if all(point_stat):
					# there is residue within cutoff distance on every other helix - set distance_flag to True
					point1.distance_flag = True
				else:
					# set flag to False: at least one helix is too far away from this point
					point1.distance_flag = False

	def show_points(self):
		print('Bundle of', len(self.helix_axis_all), 'axis')
		print('Helix order:', self.neighbour_interactions)
		print('Average distance/angle of layers:', self.average_dist_angle)
		print('Average total distance of layers:', self.average_layers_dist)
		print('Average distance to plane of points in layers:', self.average_dist_to_plane)
		print('Layer set aggregated rank:', self.ranks)
		for hid, helix in enumerate(self.helix_axis_all):
			print('Axis', hid)
			for r in helix.axis_points:
				print('Point', r.point_id  ,'from:', r.res, ', dist_flag:', r.distance_flag, ', chain:', r.chain_name, '(' + str(r.helix_id) + ')')
			print('---'*5)

class helixAxisClass():
	"""class for storing helix axis - composed of helix points"""

	def __init__(self, helix_axis, helix_id=None):
		if isinstance(helix_axis, chainClass):
			self.axis_points = [ axisPoint(pid,data[0],data[1],helix_id) for pid, data in enumerate(zip(helix_axis.axis, helix_axis.res)) ]
		elif isinstance(helix_axis, list):
			self.axis_points = helix_axis

	def __repr__(self):
		return 'Chain of len ' + str(len(self.axis_points))

	def __getitem__(self, index):
		if isinstance(index, slice):
			return helixAxisClass(self.axis_points[index])
		else:
			return self.axis_points[index]

	def length(self, dst_flag=False):
		if dst_flag:
			return len([ p for p in self.axis_points if p.distance_flag ])
		else:
			return len(self.axis_points)

class axisPoint():
	"""class for storing helix points"""

	def __init__(self, pid, axis_point, residue, helix_id):
		self.point_id = pid
		if axis_point:
			self.coords = axis_point.get_array()
		else:
			self.coords = []
		self.CA_coords = residue.Ca.get_array()
		self.distance_flag = None
		self.res = residue.res
		self.chain_name = residue.chain_name
		self.helix_id   = helix_id

	def __repr__(self):
		if self.coords != []:
			return 'Axis point {0} at {1:.2f},{2:.2f},{3:.2f}'.format(self.point_id, *[float(c) for c in self.coords])
		else:
			return 'Axis point {0} at [no coords]'.format(self.point_id)

class searchLayer():
	"""class for storing layer data in search for min distance layers"""

	def __init__(self, layer_points, neighbour_interactions):
		self.axis_points    = [ l for l in layer_points ]
		self.total_distance = self.calculate_total_distance(neighbour_interactions)

	def __repr__(self):
		return 'Layer# Tot dist: {0:.2f}  Point ids: {1} Res ids: {2}'.format(self.total_distance, str([ p.point_id for p in self.axis_points ]), str([ p.res.get_id()[1] for p in self.axis_points ]))

	def __getitem__(self, index):
		return self.axis_points[index]

	def __lt__(self, other):
		return self.total_distance < other.total_distance

	def calculate_total_distance(self, neighbour_interactions):
		"""Calculate total distance between neighbouring points in layer
		neighbourhood is determined according to neighbour_interactions list
		"""

		total_distance = 0

		for axis_point1 in self.axis_points:
			for axis_point2 in self.axis_points:
				# if (axis_point1.helix_id != axis_point2.helix_id):
				if (axis_point1.helix_id, axis_point2.helix_id) in neighbour_interactions:
					total_distance += distance.euclidean(axis_point1.coords, axis_point2.coords)

		return total_distance
