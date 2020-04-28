import numpy as np
import itertools
from Bio.PDB import calc_dihedral
from Bio.PDB import vector_to_axis
from collections import deque
from .bundle import get_local_axis
from .helper_functions import diffangle
from .helper_functions import gen_expected_crick_angles
from .helper_functions import crick_to_pos
from .helper_functions import calc_crick_ang_dev
from .helper_functions import savitzky_golay

class chainClass():
	"""Class describing single chain (helix) of the structure"""

	def __init__(self, residues, ap, chain, heptad='x'):
		self.res = residues
		for r in self.res:
			r.parentChain = self
		self.ap = ap
		self.chain = chain
		self.heptad = heptad

	def __str__(self):
		return "%s [%s-%s]" % (self.chain, self.res[0].res.id[1], self.res[-1].res.id[1])

	def get_seq(self):
		"""Get chain sequence"""

		s = "".join([seq1(r.res.resname) for r in self.res])
		if self.ap: s=s[::-1]
		return s

	def get_CA(self):
		"""Returns a list of helices. Each helix is a list of CA atoms.
		Each CA atom is represented as 3-dimensional Biopython Vector object.
		"""

		return [r.Ca for r in self.res]

	def calc_axis(self, smooth=False):
		"""Calculate axis"""

		temp_axis = get_local_axis([i.Ca for i in self.res])

		if smooth:
			self.axis = [temp_axis[0]] + [np.sum(w[0]+Vector(w[1]._ar*2)+w[2])**(1.0/4) for w in window(temp_axis, n=3)] + [temp_axis[-1]]
		else:
			self.axis = temp_axis

		# None values at edges to assure equal length of parameters tables
		self.axis = [None] + self.axis + [None]

		for pos, ax in enumerate(self.axis):
			self.res[pos].O = ax

	def calc_radius(self):
		"""Calculates distance between helical axis points and bundle axis points"""

		for pos in range(1, len(self.res)-1):
			self.res[pos].radius = (self.res[pos].O - self.res[pos].C).norm()

		self.res[0].radius = None
		self.res[-1].radius = None

		radius = [r.radius for r in self.res[1:-1]]
		self.radius = [None]+radius+[None]

	def calc_periodicity(self, smooth=False, helix_p=3.63):
		"""Calculates bundle periodicity from the perspective of a helix

		Arguments:
			helix_p: None - use local helix periodicity for caluclations
				     value - use a given fixed value of helix periodicity (default=3.63; alpha-helix periodicity)
		"""

		for pos, r in enumerate(self.res[1:-1]):
			# real position in the array
			realpos=pos+1

			### Calculate helical phase yield per residue (w1)
			# help on dihedral angles:
			# https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
			temp_w1 = []
			for other in [1, -1]:

				if (realpos+other>0) and (realpos+other<len(self.res)-1):
					temp_w1.append(np.degrees(calc_dihedral(r.Ca,\
															r.O,\
															self.res[realpos+other].O,\
															self.res[realpos+other].Ca)))

			self.res[realpos].w1 = np.mean(temp_w1)
			self.res[realpos].p = 360.0/self.res[realpos].w1

			### Calculate coiled coil phase yield per residue (w0)
			# https://github.com/numpy/numpy/issues/8529
			# Superhelical frequency (negative means left-handed superhelix)
			temp_w0 = []

			for other in [1, -1]:
				if (realpos+other>0) and (realpos+other<len(self.res)-1):
					with np.errstate(invalid="ignore"):
						temp_w0.append(np.degrees(calc_dihedral(r.O,\
																r.C,\
																self.res[realpos+other].C,\
																self.res[realpos+other].O)))

			self.res[realpos].w0 = np.mean(temp_w0)

			### Calculate bundle periodicity
			if helix_p==None: # use local helix periodicity (p) to calculate bundle periodicity (P)
				self.res[realpos].P = self.res[realpos].p / (1-(self.res[realpos].p*self.res[realpos].w0)/360)
			else: # use fixed p value
				self.res[realpos].P = helix_p / (1-(helix_p * self.res[realpos].w0)/360)

		self.res[0].P, self.res[-1].P = None, None
		self.p = [None]+[i.p for i in self.res[1:-1]]+[None]

		if smooth:
			self.P = [None]+list(savitzky_golay([r.P for r in self.res[1:-1]],5,1))+[None]
			# Update P in the residue objects
			for pos, i in enumerate(self.P):
				self.res[pos].P=i
		else:
			self.P = [None]+[r.P for r in self.res[1:-1]]+[None]

	def calc_crick(self):
		"""Calculates Crick angles for the chain."""

		def ClosestPointOnLine(a, b, p):
			ap = p-a
			ab = b-a
			result = a + np.dot(ap,ab)/np.dot(ab,ab) * ab
			return result

		def angle(v1, v2):
			return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

		for pos, r in enumerate(self.res[1:-1]):
			realpos=pos+1

			isFirstLayer = realpos==1
			isLastLayer = realpos==len(self.res)-2

			# parallel
			if not self.ap:
				if isFirstLayer:
					O_next = self.res[realpos+1].O.get_array()
					C_next = self.res[realpos+1].C.get_array()
				else:
					O_next = self.res[realpos-1].O.get_array()
					C_next = self.res[realpos-1].C.get_array()

			# anti-parallel
			else:
				if isLastLayer:
					O_next = self.res[realpos-1].O.get_array()
					C_next = self.res[realpos-1].C.get_array()
				else:
					O_next = self.res[realpos+1].O.get_array()
					C_next = self.res[realpos+1].C.get_array()

			O = r.O.get_array()
			Ca = r.Ca.get_array()
			C = r.C.get_array()

			# Project Ca onto the helix axis (O_new)
			O_new = ClosestPointOnLine(O, O_next, Ca)

			# Define a plane perpendicular to the helix axis
			# and find intersection between this plane and the bundle axis (C_new)
			n = O-O_new
			V0 = O_new
			w = C - O_new
			u = C_next - C
			N = -np.dot(n, w)
			D = np.dot(n, u)
			sI = N / D
			C_new = C+ sI * u

			# Define sign of the Crick angle
			mixed = np.dot(np.cross(O_next - O_new, Ca - O_new), C_new - O_new)

			if mixed<0:
				sign=1
			else:
				sign=-1

			if not self.ap:
				if not isFirstLayer: sign = sign * -1
			else:
				if not isLastLayer: sign = sign * -1

			r.crick = np.degrees(angle(Ca-O_new, C_new-O_new)) * sign

		self.res[0].crick=None
		self.res[-1].crick=None

		crick = [r.crick for r in self.res[1:-1]]
		self.crick = [None]+crick+[None]


	def old_calc_crick(self, fix_Crick=True):

		"""
		Currently not used!
		"""

		pass

		"""

		f_problem = []

		for pos, r in enumerate(self.res[1:-1]):

			realpos=pos+1

			OCa = r.Ca - r.O
			OC = r.C - r.O

			r.OCa = OCa
			r.OC = OC

			#Calculate the mixed product of vectors (O-Ca,O-C,O-Ox).
			#Just consider the vector between the current helical axis point
			#(O in your notation) and the next one (Ox).
			#If the mixed product is negative, reverse the sign of the Crick angle
			#(or vice versa, depending on definitions).

			if not self.ap:
				if realpos==1:
					OY = (self.res[realpos+1].O - r.O)
				else:
					OY = (self.res[realpos-1].O - r.O)
			else:
				if realpos==len(self.res)-2:
					OY = (self.res[realpos-1].O - r.O)
				else:
					OY = (self.res[realpos+1].O - r.O)

			mixed = (OY ** OCa) * OC

			if mixed<0:
				sign=1
			else:
				sign=-1

			if not self.ap:
				if realpos!=1: sign = sign * -1
			else:
				if realpos!=len(self.res)-2: sign = sign * -1

			r.crick = np.degrees(OCa.angle(OC)) * sign

			# defining the sign of the Crick angle close to the +/-180 deg. may be problematic
			# such cases are stored for further check-up
			if abs(mixed)<5 and abs(r.crick)>150:
				f_problem.append(realpos)

		self.res[0].crick=None
		self.res[-1].crick=None

		print(f_problem)

		if not fix_Crick: f_problem = []

		#iterate over potentially problematic Crick angles
		for fp in f_problem:

			expected = []

			# calculate w1 values in the neighbouring residue pairs
			if fp>2:
				expected.append(abs(diffangle(self.res[fp-1].crick, self.res[fp-2].crick)))
			if fp<len(self.res)-3:
				expected.append(abs(diffangle(self.res[fp+1].crick, self.res[fp+2].crick)))

			if len(expected)==0: continue

			expected = np.mean(expected)

			current = []
			alternative = []

			if fp>1:
				current.append(abs(diffangle(self.res[fp].crick, self.res[fp-1].crick)))
				alternative.append(abs(diffangle(-self.res[fp].crick, self.res[fp-1].crick)))
			if fp<len(self.res)-2:
				current.append(abs(diffangle(self.res[fp].crick, self.res[fp+1].crick)))
				alternative.append(abs(diffangle(-self.res[fp].crick, self.res[fp+1].crick)))

			if max(np.absolute(current - expected)) > max(np.absolute(alternative - expected)):
				print(fp, 'fixed!')
				self.res[fp].crick = -self.res[fp].crick

		crick = [r.crick for r in self.res[1:-1]]
		self.crick = [None]+crick+[None]

		"""

	def calc_pitch_angle(self):
		"""Calculate pitch angle"""

		for pos in range(len(self.res)-2):
			realpos=pos+1

			if realpos==1:
				other=1
			else:
				other=-1

			O = self.res[realpos].O - self.res[realpos+other].O
			C = self.res[realpos].C - self.res[realpos+other].C

			self.res[realpos].A = np.degrees(C.angle(O))
		self.res[0].A=None
		self.res[-1].A=None
		self.A = [None]+[r.A for r in self.res[1:-1]]+[None]

	def assign_positions(self):
		"""Assigns heptad positions using TWISTER algorithm."""

		hpos=['?', '?']
		pos = 2

		if self.ap:
			crick = self.crick[::-1]
		else:
			crick = self.crick

		while pos < len(self.res)-2:

			if crick[pos-1] < 0 and crick[pos] > 0 and abs(crick[pos-1]) > abs(crick[pos]):
				hpos.extend(['a', 'b', 'c'])
				pos+=3
			elif crick[pos] < 0 and crick[pos+1] > 0 and abs(crick[pos]) < abs(crick[pos+1]):
				hpos.extend(['d', 'e', 'f', 'g'])
				pos+=4
			else:
				hpos.append('?')
				pos+=1


		hpos = hpos[:len(self.res)-2]

		# add missing positions at N end
		fpos=0
		try:
			while hpos[fpos]=='?': fpos+=1
		except IndexError:
			pass
		else:
			hep = itertools.cycle('gfedcba')
			while not next(hep) == hpos[fpos]: pass
			fadd = [next(hep) for i in range(fpos)][::-1]
			hpos = fadd + hpos[fpos:]

		# add missing positions at C end
		epos=len(hpos)-1
		try:
			while hpos[epos]=='?': epos-=1
		except IndexError:
			hpos.extend(['?', '?'])
		else:
			hep = itertools.cycle('abcdefg')
			while not next(hep) == hpos[epos]: pass
			eadd = [next(hep) for i in range(len(self.res)-epos-1)]
			hpos = hpos[:epos+1] + eadd

		assert len(hpos) == len(self.res)

		# store assigned heptad positions in the object
		if self.ap:
			hpos = hpos[::-1]
		self.positions_twister = hpos

		# Old approach:
		#angles = gen_expected_crick_angles(P, REP, optimal_ph1)
		#self.positions = ['?'] + [crick_to_pos(c, angles)[1] for c in self.crick[1:-1]] + ['?']

	def calc_crickdev(self, P, REP, optimal_ph1=19.5, smooth=False):

		angles = gen_expected_crick_angles(P, REP, optimal_ph1)
		posnames = [chr(97+i) for i in range(len(angles))]

		g = deque(zip(angles, posnames))

		best_fit=1000000

		for i in range(len(g)):
			unpack_g = np.asarray(list(g))

			heptad_cycle = itertools.cycle(unpack_g[:,0].astype(float))

			measured_circk = [r.crick for r in self.res[1:-1]]
			exp_crick = [next(heptad_cycle) for _ in range(len(measured_circk))]

			if self.ap:
				exp_crick = exp_crick[::-1]

			crdev = calc_crick_ang_dev(measured_circk, exp_crick, 0, len(self.res)-2, force_start_Cr_ang_pos=0)[0]

			fit = np.mean([abs(i) for i in crdev])

			if fit < best_fit:
				best_crdev = crdev
				posname_cycle = itertools.cycle(unpack_g[:,1].astype(str))
				best_posnames = [next(posname_cycle) for _ in range(len(best_crdev)+1)]

				first_name = ord(best_posnames[0])-1
				if first_name<97:
					assert first_name==96
					first_name = posnames[-1]
				else:
					first_name = chr(first_name)

				best_posnames = [first_name] + best_posnames

				if self.ap:
					best_posnames = best_posnames[::-1]
				best_fit = fit

			g.rotate(1)

		if smooth:
			self.crdev = [None]+list(savitzky_golay(best_crdev,5,1))+[None]
		else:
			self.crdev = [None]+best_crdev+[None]

		self.positions_samcc = best_posnames

		for pos, i in enumerate(self.crdev):
			self.res[pos].crdev = i



	def calc_crickdev_old(self, P, REP, optimal_ph1=19.5, smooth=False):
		"""
		Calculates Crick Angle deviation.

		Arguments:
			P (int): bundle periodicity
			REP (int): repeat length
		"""

		angles = gen_expected_crick_angles(P, REP, optimal_ph1)

		# No starting heptad for index=0 (no data!) was defined. Define it based on the Crick angle
		if self.heptad == 'x':

			ref_crick = self.res[1].crick
			ref_crick_pos = 1

			if abs(ref_crick)>=165:
				#print(ref_crick, self.res[2].crick)
				ref_crick = self.res[2].crick
				ref_crick_pos = 2

			start_Cr_ang_pos, name, _ = crick_to_pos(ref_crick, angles)
			if not self.ap:
				hpos = start_Cr_ang_pos-ref_crick_pos
			else:
				hpos = start_Cr_ang_pos+ref_crick_pos

			if hpos<0: hpos=REP+hpos  # was -1 -> +hpos
			if hpos>=REP: hpos=hpos-REP   # was == -> >=  ; 0 -> hpos-REP
		else:
			hpos = ord(self.heptad)-97

		# Define heptad positions for all residues
		heptads = ""
		h = hpos
		for _ in range(len(self.res)):
			heptads += chr(h+97)
			if self.ap:
				h-=1
				if h<0: h = REP-1
			else:
				h+=1
				if h==REP: h = 0

		self.heptads = list(heptads)

		hep2angle = dict([(chr(l+97), a) for l,a in enumerate(angles)])

		if self.ap:
			hpos-=1
		else:
			hpos+=1

		if hpos<0:
			hpos = REP-1
		elif hpos==REP:
			hpos = 0

		heptad = chr(hpos+97)

		exp_ahelix_crick = gen_expected_crick_angles(P, REP, hep2angle[heptad], ap=self.ap)

		crdev = calc_crick_ang_dev([r.crick for r in self.res[1:-1]], exp_ahelix_crick, 0, len(self.res)-2,\
										   force_start_Cr_ang_pos=0)[0]


		if smooth:
			self.crdev = [None]+list(savitzky_golay(crdev,5,1))+[None]
		else:
			self.crdev = [None]+crdev+[None]

		for pos, i in enumerate(self.crdev):
			self.res[pos].crdev = i

	def calc_axialshift(self):
		"""C = prevnewC - newC # os CC (nowo wyliczone wartosci z helis, ktorych osie sa przesuniete)
        O = positions[pos]['alphaaxis'][chainnr] - newC # nowa os helisy
        Omove = vector_to_axis(C, O) - O
		"""

		for pos, r in enumerate(self.res[1:-1]):
			realpos=pos+1


			if realpos==1:
				pre = self.res[realpos+1].C
				C = pre - r.C
			else:
				pre = self.res[realpos-1].C
				C = pre - r.C

			O = r.O - r.C

			preOmove = vector_to_axis(C, O)
			Omove = preOmove - O

			temp_newO = r.O + Omove
			temp_distance = (r.O - pre).norm() - (temp_newO - pre).norm()

			r.shift = Omove.norm() * np.sign(temp_distance)
			if realpos==1:
				r.shift*=-1

		self.res[0].shift = None
		self.res[-1].shift = None

		self.shift = [i.shift for i in self.res]
