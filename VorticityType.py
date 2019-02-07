"""
Author: Zack Taylor
Date: 2/7/19

The Vorticity package used by the physics driver
"""

class Vorticity(object):

	"""
	@brief Initializes the Vorticity object

	@param Mesh		The Mesh object
	"""
	def __init__(self, mesh, dt=0.1, Re=10.0):
		self.mesh = mesh
		self.dt = dt
		self.Re = Re
		self.CoeffA = None
		self.CoeffB = None
		self.CoeffC = None
		self.CoeffD = None
		self.CoeffE = None

	"""
	@brief Sets the coefficients for the A matrix
	"""
	def __setACoefficients(self, node):
		v = node.xVelocity
		u = node.yVelocity
		dy = self.mesh.dy
		dx = self.mesh.dx
		Re = self.Re
		dt = self.dt

		self.CoeffA = v/(2.*dy) - 1./(Re*dy**2)
		self.CoeffB = u/(2.*dx) - 1./(Re*dx**2)
		self.CoeffC = 1./dt + 2./(Re*dx**2) + 2./(Re*dy**2)
		self.CoeffD = -u/(2.*dx) - 1./(Re*dx**2)
		self.CoeffE = -v/(2*dy) - 1./(Re*dy**2)

	"""
	@Brief Sets the coeffients for the b matrix 
	"""
	def __setNodeSource(self):
		for node in self.mesh.nodes:
			if not node.solved:
				node.source = node.solutionVorticity/self.dt

	"""
	@Brief Loops through the model to update the velocity
	"""
	def __updateVelocities(self):
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(K)
			u = (node.north.Vorticity-node.south.Vorticity)/(2*self.mesh.dy)
			v = (node.east.Vorticity-node.west.Vorticity)/(2*self.mesh.dx)
			node.yVelocity = u
			node.xVelocity = v

	"""
    @Brief Sets up the A matrix for a 5 point grid
	"""
	def getAMatrix(self):
		self.__updateVelocities()
		numOfColumns = self.mesh.numOfxNodes-2
		numOfRows = self.mesh.numOfyNodes-2
		numOfUnknowns = numOfColumns*numOfRows
		A = np.zeros((numOfUnknowns,numOfUnknowns))

		# diagonal
		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeBySolIndex(i)
			self.__setACoefficients(node)
			# main diagonal
			A[i,i] = self.CoeffC

			# the off off diagonals
			if i+numOfRows < numOfUnknowns:
				A[i,i+numOfRows] = self.CoeffE
			if i+1 > numOfRows:
				A[i,i-numOfRows] = self.CoeffA

			# the off diagonals
			if (i+1)%numOfRows:
				A[i+1,i] = self.CoeffD
				A[i,i+1] = self.CoeffB

		return A

	"""
	@Brief Builds and returns the b vector 
	"""
	def getbVector(self):
		self.__setNodeSource()
		B = np.zeros((self.mesh.maxSolIndex,1))
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			B[k] = node.source

		return B 
