import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import pandas as pd
import warnings
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from .layerClass import layerClass
from .residueClass import residueClass
from .chainClass import chainClass
from .pymol_draw_layers import save_layers_to_pymol
#plt.switch_backend('SVG')

class bundleClass():
	"""Class describing single bundle of the structure"""

	def __init__(self):
		self.layers = []
		self.pymol_selection = []
		self.helix_order = []
		self.helices_axis = []
		self.ppo = []

	def calc_bundleaxis(self):
		"""Calculate axis of entire bundle"""

		axis=[]

		for l in self.layers[1:-1]:
			res_in_layer = [r.O for r in l.res]
			C = np.sum(res_in_layer) ** (1.0/len(self.chains))
			for r in l.res:
				r.C = C
			axis.append(C)

		self.axis = [None]+axis+[None]

	def get_helicesaxis(self):
		"""Get axis of individual helices"""

		for c in self.chains:
			self.helices_axis.append(c.axis)

	def calc_pitch_angle(self):
		"""Calculate pitch angle"""

		for c in self.chains:
			c.calc_pitch_angle()

	def calc_crick(self):
		"""Calculates Crick angles for all chains in the bundle."""

		for c in self.chains:
			c.calc_crick()

	def calc_crickdev(self, P, REP, optimal_ph1=19.5, smooth=False):
		"""Calculates Crick's angle deviations for all chains in a bundle
		P - periodicity (e.g. 3.5 for 7/2)
		REP - repeat length (eg. 7 for 7/2)
		"""
		for c in self.chains:
			c.calc_crickdev(P, REP, optimal_ph1=optimal_ph1, smooth=smooth)

	def calc_periodicity(self, smooth=False, helix_p=3.63):
		"""Calculates periodicity for all chains in a bundle."""
		for c in self.chains:
			c.calc_periodicity(smooth=smooth, helix_p=helix_p)

	def calc_radius(self):
		"""Calculates radius for all chains in a bundle."""
		for c in self.chains:
			c.calc_radius()

	def calc_axialshift(self):
		"""Calculates axial shift for all chains in a bundle."""
		for c in self.chains:
			c.calc_axialshift()

	def assign_positions(self):
		"""Assign heptad positions using TWISTER algorithm."""
		for c in self.chains:
			c.assign_positions()

	def write_helical_axes(self, filename):
		"""Writes helical axes in PDB format."""

		sb = StructureBuilder()
		sb.init_structure('')
		sb.init_model('')

		for cpos, chain in enumerate(self.chains):
			sb.init_chain(str(cpos))
			sb.init_seg('')

			for pos, i in enumerate(chain.res[1:-1]):
				sb.init_residue('ALA',' ',pos,' ')
				sb.init_atom('CA', i.O._ar,  0, 0, " ", ' CA ', 1)

		io=PDBIO()
		io.set_structure(sb.structure)
		io.save(filename)

	def from_defdata(self, filename, chains, chains_names, chains_ap, chains_heptad):
		"""Reads bundle based on provided parameters.

			Arguments:

				filename: name of the PDB file
		"""

		parser=PDBParser()
		structure=parser.get_structure('temp', filename)

		assert len(set([len(i) for i in chains]))==1, "Provided chains are not of equal length! Aborting!"

		# read layers
		for layer in zip(*chains):
			l=layerClass([ residueClass(structure[0][c][r], c) for r, c in zip(layer, [x[0] for x in chains_names]) ])
			self.layers.append(l)

		# read chains
		self.chains = [list([]) for _ in range(len(chains_names))]

		for l in self.layers:
			for pos, r in enumerate(l.res):
				self.chains[pos].append(r)

		for pos, c in enumerate(self.chains):
			self.chains[pos] = chainClass(self.chains[pos], chains_ap[pos], chains_names[pos], chains_heptad[pos])
			self.chains[pos].calc_axis()

	def gendf(self, cols=['radius', 'crick', 'P', 'shift', 'crdev', 'positions_twister', 'p', 'positions_samcc']):
		"""Generates pandas df containg all bundle parameters."""

		index, data = [], []
		warned = []

		for chainpos, chain in enumerate(self.chains):
			for rpos in range(len(chain.radius)):

				index.append((rpos, chainpos))

				row_dict = {'res_name':chain.res[rpos].res.get_resname(),
							'res_number':chain.res[rpos].res.full_id[3][1],
							'chain_name':chain.res[rpos].chain_name,
							'chain_ap':chain.ap}

				for c in cols:
					try:
						row_dict[c] = getattr(chain, c)[rpos]
					except AttributeError:
						if not c in warned:
							warnings.warn("%s data missing" % c)
							warned.append(c)

				data.append(row_dict)

		index = pd.MultiIndex.from_tuples(index, names=['layer', 'chain'])

		return pd.DataFrame(data, index=index).sort_index()

	def plot(self, filename="", elements=[]):

		"""Plot data measured by the layer"""

		assert len(elements)>1
		assert type(elements)==list

		sns.set_style("ticks")

		TICKFONTSIZE = 15
		AXESFONTSIZE = 15

		fig, axes = plt.subplots(nrows=len(elements), sharey=False, sharex=False, figsize=(15, len(elements)*3))
		a=0

		plt.subplots_adjust(hspace=0.5)

		colgen = itertools.cycle('gbmr')

		for e in elements:
			for chain, color in zip(self.chains, colgen):

				if e=='Periodicity':
					# Plot bundle periodicity
					axes[a].plot(chain.P[1:-1], c=color, linewidth=1.0)
					#for i in [17/5., 7/2., 25/7., 18/5., 11/3., 15/4., 19/5.]:
					#	axes[0].axhline(ls='--', color='0.5', y=i, linewidth=1.0)
					#axes[a].axhline(ls='--', color='0.5', y=np.mean(chain.P[1:-1]), linewidth=2.0)
					axes[a].set_ylabel(r'Periodicity', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(3.3, 3.9)
					axes[a].yaxis.set_ticks(np.arange(3.3, 3.9, 0.1))

				elif e=='Radius':
					# Plot radius
					axes[a].plot(chain.radius[1:-1], c=color, linewidth=1.0)
					axes[a].set_ylabel(r'Radius', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(6, 9)
					axes[a].yaxis.set_ticks(np.arange(5, 9, 0.5))

				elif e=='CrickDev':
					# Plot axial rotation
					axes[a].plot(chain.crdev[1:-1], c=color, linewidth=1.0)
					axes[a].set_ylabel(r'Helix axial rotation', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(-26, 26)
					axes[a].yaxis.set_ticks(np.arange(-26, 26, 10))

				elif e=='Shift':
					# Plot axial shift
					axes[a].plot(chain.shift[1:-1], c=color, linewidth=1.0)
					axes[a].set_ylabel(r'Axial shift', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(-5, 5)
					axes[a].yaxis.set_ticks(np.arange(-5, 5, 1))

				elif e=='w0':
					axes[a].plot([i.w0 for i in chain.res[1:-1]], c=color, linewidth=1.0)
					axes[a].set_ylabel(r'w0', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(-5, 5)
					axes[a].yaxis.set_ticks(np.arange(-5, 5, 1))

				elif e=='w1':
					axes[a].plot([i.w1 for i in chain.res[1:-1]], c=color, linewidth=1.0)
					axes[a].set_ylabel(r'w1', fontsize=AXESFONTSIZE)
					axes[a].set_ylim(90, 110)
					axes[a].yaxis.set_ticks(np.arange(90, 120, 5))

				else:
					pass

				if a == (len(elements)-1):
					axes[a].set_xlabel(r'Position in the sequence', fontsize=AXESFONTSIZE)
			a+=1

		for ax in axes:
			ax.xaxis.set_ticks(np.arange(0,len(self.layers)-2,1))

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(TICKFONTSIZE)

		if filename=="":
			plt.show()
		else:
			plt.savefig(filename, dpi=150)

	def pymol_plot_layer(self, filename, savepath, suffix, pymol_version, color_selection=False, helix_order=False, helices_axis=None):
		"""Save layers to pymol session.

		Arguments:
			filename: path to pdb file
			savepath: path to where pse file will be saved
			suffix: suffix appended to the name of saved file (pdb_id)
			pymol_version: version of pymol used to generate session (differences in syntax between 1.x and 2.x)
			color_selection: should bundles be highlighted in pymol session

		Returns:
			None, session is saved to file
		"""

		# plot only interaction between helices that were used to determine layers
		if helix_order:
			h_order = self.helix_order
		else:
			h_order = None

		# plot layers by axis points
		layer_points = [ layer.get_layer_axis() for layer in self.layers if (len(layer.get_layer_axis()) == len(layer.res)) ]

		if color_selection:
			save_layers_to_pymol(filename, layer_points, savepath, suffix, pymol_version, color_selection=self.pymol_selection,
								 helix_order=h_order, ppo=self.ppo, bundle_axis=self.axis, helices_axis=self.helices_axis)
		else:
			save_layers_to_pymol(filename, layer_points, savepath, suffix, pymol_version, color_selection=False,
								 helix_order=h_order, ppo=self.ppo, bundle_axis=self.axis, helices_axis=self.helices_axis)
