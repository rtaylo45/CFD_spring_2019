"""
Author: Zack Taylor
Date: 2/7/19

The Laplace package used by the physics driver
"""

class Laplace(object):

	"""
	@brief Initializes the Laplace object

	@param Mesh		The Mesh object
	"""
	def __init__(self, mesh):
		self.mesh = mesh
		self.CoeffA = None
		self.CoeffB = None
		self.CoeffC = None
		self.CoeffD = None
		self.CoeffE = None

		self.__setACoefficients()
		self.__setNodeSource()

	"""
	@brief Sets the coefficients for the A matrix
	"""
	def __setACoefficients(self):
		self.CoeffA = 1./(self.mesh.dy**2.)
		self.CoeffB = 1./(self.mesh.dx**2.)
		self.CoeffC = -(2./(self.mesh.dx**2.) + 2./(self.mesh.dy**2.))
		self.CoeffD = 1./(self.mesh.dx**2.)
		self.CoeffE = 1./(self.mesh.dy**2.)

	"""
	@Brief Sets the coeffients for the b matrix 
	"""
	def __setNodeSource(self):
		for node in self.mesh.nodes:
			if not node.solved:
				node.source = 0.0
				if node.east.solved:
					node.source += -node.east.solution/self.mesh.dx**2.
				if node.west.solved:
					node.source += -node.west.solution/self.mesh.dx**2.
				if node.north.solved:
					node.source += -node.north.solution/self.mesh.dy**2.
				if node.south.solved:
					node.source += -node.south.solution/self.mesh.dy**2.
