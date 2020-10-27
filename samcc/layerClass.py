class layerClass():
	"""Class describing single layer of the bundle"""

	def __init__(self, residues):
		self.res = residues

	def __str__(self):
		return " ".join(["%s:%s" % (r.chain_name, r.res.id[1]) for r in self.res])

	def get_layer_CA(self):
		return [ [res.Ca[0], res.Ca[1], res.Ca[2] ] for res in self.res ]

	def get_layer_axis(self):
		return [ [res.O[0], res.O[1], res.O[2] ] for res in self.res if res.O != None ]

	def get_layer_axis_temp(self):
		return [ type(res.O) for res in self.res ]
