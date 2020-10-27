import sys
import math
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import StructureBuilder
from Bio.PDB import PDBIO
from Bio.PDB import vector_to_axis
from Bio.PDB import Vector

def get_local_axis(calphas):
	"""Get local axis"""

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

			dmag = math.sqrt(dv13._ar[0]*dv13._ar[0] + dv13._ar[1]*dv13._ar[1] + dv13._ar[2]*dv13._ar[2])
			emag = math.sqrt(dv24._ar[0]*dv24._ar[0] + dv24._ar[1]*dv24._ar[1] + dv24._ar[2]*dv24._ar[2])
			dot_product = dv13 * dv24

			costheta = dot_product/(dmag*emag)
			costheta1 = 1.0 - costheta
			radmag = math.sqrt(dmag*emag) / (costheta1*2.0)
			dv13 = dv13 ** (1/dmag)
			rad = dv13**radmag
			origin[j] = calphas[j+1] - rad

			dv24 = dv24 ** (1/emag)
			rad = dv24**radmag
			origin[j+1] = calphas[j+2] - rad

	else:
		print("Helix is less than 4 residues and has no axis\n\n")
		sys.exit(-1)

	return list(origin.values())
