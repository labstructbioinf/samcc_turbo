import sys, math
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import StructureBuilder
from Bio.PDB import PDBIO
from Bio.PDB import vector_to_axis
from Bio.PDB import Vector


DEBUG = False

def get_local_axis(calphas):

	#FIXME add docs, consider more meaningful variable names or add comments

	unit_local_axis = {}
	origin = {}
	if len(calphas)>=4:
		for j in range(len(calphas)-3):
			vec12 = calphas[j+1] - calphas[j]
			vec23 = calphas[j+2] - calphas[j+1]
			vec34 = calphas[j+3] - calphas[j+2]
			dv13 = vec12 - vec23
			dv24 = vec23 - vec34
			cross_product = (dv13 ** dv24).normalized()
			unit_local_axis[j] = cross_product

			if DEBUG:
				print(j,unit_local_axis[j])

			dmag = math.sqrt(dv13._ar[0]*dv13._ar[0] + dv13._ar[1]*dv13._ar[1] + dv13._ar[2]*dv13._ar[2])
			emag = math.sqrt(dv24._ar[0]*dv24._ar[0] + dv24._ar[1]*dv24._ar[1] + dv24._ar[2]*dv24._ar[2])
			dot_product = dv13 * dv24

			costheta = dot_product/(dmag*emag)
			costheta1 = 1.0 - costheta
			radmag = math.sqrt(dmag*emag) / (costheta1*2.0)
			dv13 = dv13 ** (1/dmag)
			rad = dv13**radmag
			origin[j] = calphas[j+1] - rad

			if DEBUG:
				print('\torigin',j, origin[j])

			dv24 = dv24 ** (1/emag)
			rad = dv24**radmag
			origin[j+1] = calphas[j+2] - rad

			if DEBUG:
				print('\torigin',j+1, origin[j+1])
				print()

	else:
		print("Helix is less than 4 residues and has no axis\n\n")
		sys.exit(-1)

	return list(origin.values())

def smooth(ls, av=1):

	#FIXME is this function still used? delete if not, otherwise add docs

	finalls = []
	for l in range(len(ls)):
		div = 0
		suma = None
		try:
			for x in range(-av, av+1):

				if l+x<0 or l+x>len(ls)-1:
					raise IndexError
				multi = av+1-abs(x)



				if suma == None:
					suma = Vector(ls[l+x]._ar*(multi))
				else:
					#if x==av:
					suma = suma + Vector(ls[l+x]._ar* ( multi ) )

				div += multi
					#else:
					#suma = suma + ls[l+x]


		except IndexError:
			div = div * 1.0
			suma = ls[l]
		else:
			div = div * 1.0
			suma = suma / div



		finalls.append(suma)
	return finalls

###FIXME are these classes used? delete if not

class mydict(dict):
	def __getattr__(self, name):
		return self[name]
	def __setattr__(self, key, value):
		self[key] = value

class bundle(mydict):

	def __init__(self, id):
		self.id = id
		self.chains = mydict()

	def neighbourslist(self):
		m = max(self.chains.keys())
		l = [(i, i+1) for i in range(m)]
		l.append((m, 0))
		return l

	def writeaxes(self):

		sb = StructureBuilder.StructureBuilder()
		sb.init_structure('')
		sb.init_model('')

		for chain in self.chains.keys():
			sb.init_chain(str(chain))
			sb.init_seg('')

			for pos, i in enumerate(self.chains[chain]['axis']):

				sb.init_residue('ALA',' ',pos,' ')
				sb.init_atom('CA', i,  0, 0, " ", ' CA ', 1)

				#print i


		io=PDBIO()
		io.set_structure(sb.structure)
		io.save("%s.axes.pdb" % self.id)

class chain(mydict):
	pass

	"""
	def __init__(self):
		self.data = dict()

	def __setitem__(self, key, data):
		self.data[key] = data

	def __getitem__(self, key):
		return self.data[key]

	"""

###

def readbundle(rangefile):

	#FIXME is this function still used?

	"""
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser()
	"""


	parser = PDBParser(PERMISSIVE=1)

	ID = rangefile.split('.')[0]

	structure = parser.get_structure(ID, ID+'.pdb')

	b = bundle(ID)


	lastlen = None

	for pos, n in enumerate(open(rangefile)):
		f, t, chainid = n.strip().split(',')
		f, t = int(f), int(t)

		print("chain %s length: %s" % (pos, f-t))

		assert f!=t
		assert lastlen==None or f-t==lastlen
		lastlen == f-t

		c = chain()

		if f>t:
			r = range(f,t-1,-1)
			c.ap = True
		else:
			r = range(f,t+1)
			c.ap = False

		#print r


		c['id'] = chainid.upper()
		c['str_range'] = "%s-%s" % (f,t)
		c['range'] = r
		c['ca'] = [structure[0][c['id']][res]['CA'].get_vector() for res in c['range']]

		#print len(c.ca)



		#print [structure[0][c['id']][res] for res in c['range']]

		#if not c.ap:
		#	c['axis'] = get_local_axis(c['ca'])[:-1]
		#else:
		c['axis'] = get_local_axis(c['ca'])


		#print len(c.axis)

		#c['ca'] = c['ca'][:-2]

		#print len(c['ca'])
		#print len(c['axis'])

		b.chains[pos] = c

	return b



if __name__ == "__main__":

	#FIXME is this still useful?

	# http://paulbourke.net/geometry/pointlineplane/
	from L3D import LineLineIntersect3D, Point
	from sympy import Point3D, Line3D, Segment3D

	bundle = readbundle(sys.argv[1])

	#raise IOError

	bundle.writeaxes()

	for chain1, chain2 in bundle.neighbourslist():


		dist= []

		for i in range(len(bundle.chains[chain1].axis)-1):


			p1a = bundle.chains[chain1].axis[i].get_array()
			p1b = bundle.chains[chain1].axis[i+1].get_array()

			p2a = bundle.chains[chain2].axis[i].get_array()
			p2b = bundle.chains[chain2].axis[i+1].get_array()


			a = LineLineIntersect3D(Point(*p1a), Point(*p1b), Point(*p2a), Point(*p2b))

			if a.on_segment1 == 1 and a.on_segment2 == 1:
				dist.append(a.inters_dist)
			else:

				p1a = Point3D(*p1a, evaluate = True)
				p1b = Point3D(*p1b, evaluate = True)
				p2a = Point3D(*p2a, evaluate = True)
				p2b = Point3D(*p2b, evaluate = True)

				s1 = Segment3D(p1a, p1b, evaluate=True)
				s2 = Segment3D(p2a, p2b, evaluate=True)




				temp = []

				P = s1.projection(p2a)
				if s1.contains(P): temp.append(p2a.distance(P))

				P = s1.projection(p2b)
				if s1.contains(P): temp.append(p2b.distance(P))

				P = s2.projection(p1a)
				if s2.contains(P): temp.append(p1a.distance(P))

				P = s2.projection(p1b)
				if s2.contains(P): temp.append(p1b.distance(P))



				temp.extend([(p1a.distance(p2a)),\
						     (p1a.distance(p2b)),\
						     (p1b.distance(p2a)),\
						     (p1b.distance(p2b))])


				dist.append(float(min(temp)))

			"""
			a = LineLineIntersect3D(makepoint(bundle.chains[chain1].axis[0]), makepoint(bundle.chains[chain1].axis[1]), makepoint(bundle.chains[chain2].axis[0]), makepoint(bundle.chains[chain2].axis[1]))

			print a.inters_dist
			print




			print "Member 1 apparent intersection point = ", a.Pmem1
			print "Member 2 apparent intersection point = ", a.Pmem2
			print "Unit vector length of segment connecting apparent intersection points = ", a.uv
			print "Absolute distance between apparent intersection points = ", a.inters_dist
			print "Member 1 Left End 'X' distance to apparent intersection point = ", a.left_dist
			print "Member 1 Right End 'X' distance to apparent intersection point = ", a.right_dist
			print "Apparent intersection point position with respect to Member 1 = ", a.position
			if a.on_segment1 == 1:
				print "Point Pa lies on Member 1"
			else:
				print "Point Pa does not lie on Member 1"
			if a.on_segment2 == 1:
				print "Point Pb lies on Member 2"
			else:
				print "Point Pb does not lie on Member 2"

			"""


			#sys.exit(-1)
		#print

		#dist = [(a1-a2).norm() for a1, a2 in zip(bundle.chains[chain1]['axis'], bundle.chains[chain2]['axis'])]
		#for d in dist:
		#	print d

		#print "%s/%s/%s" % (round(np.mean(dist),1), round(np.median(dist),1), round(min(dist),1))

		print("%s [%s] - %s [%s]: %s" % (bundle.chains[chain1].id, bundle.chains[chain1].str_range, bundle.chains[chain2].id, bundle.chains[chain2].str_range, round(min(dist),1)))
		#print "%s-%s: %s" % (bundle.chains[chain1].id, bundle.chains[chain2].id, round(min(dist),1))


		#print
