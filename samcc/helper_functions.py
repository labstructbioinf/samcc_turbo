### Helper functions ###

import networkx as nx
import itertools
import numpy as np
from math import factorial

def gen_expected_crick_angles(P, rep_len, start_ph1, ap=False):

	step = 360 / P
	if ap:
		sign=-1
	else:
		sign=1
	return [adj(start_ph1+(sign * i * float(step))) for i in range(rep_len)]

def adj2(ang):

	ang = adj(ang)
	if ang < -180:
		ang += 360
	return ang

def adj(ang):

	ang += 180.0
	if ang>360:
		c = int(ang / 360.0)
		return ang - (360*c) - 180.0
	else:
		return ang - 180.0

def adjustangle(angle):

	if abs(angle)>180 and angle > 0:
		angle = angle - 360.0
	elif abs(angle)>180 and angle < 0:
		angle += 360.0
	return angle

def diffangle(targetA, sourceA):

	a = targetA - sourceA
	a = (a + 180) % 360 - 180
	return a

def crick_to_pos(start_Cr_ang, exp_helix_crick):

	diff = [abs(adj(start_Cr_ang-i)) for i in exp_helix_crick]
	mindiff = min(diff)

	start_Cr_ang_pos = diff.index(mindiff)

	try:
		name = chr(97+start_Cr_ang_pos)
	except ValueError:
		name = '?'

	return start_Cr_ang_pos, name, mindiff

def calc_crick_ang_dev(twister, exp_helix_crick, firstpos, lastpos, force_start_Cr_ang_pos=None):

	if firstpos == lastpos:
		assert lastpos == None

	# define Crick angle of the starting position
	if force_start_Cr_ang_pos==None:
		if lastpos==None:
			firstpos = 0
		if type(twister) is list:
			start_Cr_ang = twister[firstpos]
		else:
			start_Cr_ang = twister.iloc[firstpos]['Cr_ang']
		start_Cr_ang_pos, name, _ = crick_to_pos(start_Cr_ang, exp_helix_crick)
	else:
		start_Cr_ang_pos = force_start_Cr_ang_pos

	cr = itertools.cycle(exp_helix_crick)

	# skip to the first repeat pos
	for n in range(start_Cr_ang_pos):
		next(cr)

	if lastpos==None:
		laspos = len(twister)
	if type(twister) is list:
		data = twister[firstpos:lastpos]
	else:
		data = twister['Cr_ang'].iloc[firstpos:lastpos]

	Cr_ang_dev = [adj2(c - next(cr)) for c in data]

	m,b = np.polyfit(range(len(Cr_ang_dev)), Cr_ang_dev, 1)

	return Cr_ang_dev, m, b

def window(seq, n=2):
	"""Returns a sliding window (of width n) over data from the iterable
	s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
	"""

	it = iter(seq)
	result = tuple(itertools.islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
	The Savitzky-Golay filter removes high frequency noise from data.
	It has the advantage of preserving the original shape and
	features of the signal better than other types of filtering
	approaches, such as moving averages techniques.
	Parameters
	----------
	y : array_like, shape (N,)
		the values of the time history of the signal.
	window_size : int
		the length of the window. Must be an odd integer number.
	order : int
		the order of the polynomial used in the filtering.
		Must be less then `window_size` - 1.
	deriv: int
		the order of the derivative to compute (default = 0 means only smoothing)
	Returns
	-------
	ys : ndarray, shape (N)
		the smoothed signal (or it's n-th derivative).
	"""

	try:
		window_size = np.abs(np.int(window_size))
		order = np.abs(np.int(order))
	except ValueError as msg:
		raise ValueError("window_size and order have to be of type int")
	if window_size % 2 != 1 or window_size < 1:
		raise TypeError("window_size size must be a positive odd number")
	if window_size < order + 2:
		raise TypeError("window_size is too small for the polynomials order")
	order_range = range(order+1)
	half_window = (window_size -1) // 2
	# precompute coefficients
	b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
	m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
	# pad the signal at the extremes with
	# values taken from the signal itself
	firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
	lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
	y = np.concatenate((firstvals, y, lastvals))
	return np.convolve( m[::-1], y, mode='valid')

# SOCKET helper functions

def detect_helices_orientation(indices, orientation_data):

	"""
	Arguments

	indices: dict of bundle helices from Socket
	orientation_data: relative orientation of the helices from indices

	Returns

	a dict helix_id -> orientation

	this function uses a graph to convert relative orientations between the helices
	to parallel/anti-parallel labels for each helix

	"""

	G=nx.Graph()
	for o in orientation_data:
		G.add_edge(o[0], o[1], o=o[2])

	res = {}

	for g in nx.connected_component_subgraphs(G):

		first=True
		for e in list(nx.dfs_edges(g)):
			if first:
				g.node[e[0]]['o'] = 1
				first = False

			if g[e[0]][e[1]]['o']=='antiparallel':
				g.node[e[1]]['o'] = g.node[e[0]]['o']*-1
			else:
				g.node[e[1]]['o'] = g.node[e[0]]['o']

		for n in g.nodes():
				res[n] = g.node[n]['o']

	return res
