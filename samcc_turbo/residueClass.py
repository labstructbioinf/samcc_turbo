class residueClass():
	
	#FIXME add docs

	def __init__(self, res, chain):
		self.res = res
		self.Ca = self.res['CA'].get_vector()
		self.chain_name = chain
		self.C = None
		self.O = None
