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
	@brief runs the presolve to apply BC
	"""
	def runPresolve(self):
		self.__applyBC()

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
	@Brief Sets the boundary conditions for each face
	"""
	def __applyBC(self):
	
		# the velocities at each boundary	
		vWest = self.mesh.westBC
		vEast = self.mesh.eastBC
		uNorth = self.mesh.northBC
		uSouth = self.mesh.southBC

		for j in xrange(self.mesh.numOfyNodes):
			for i in xrange(self.mesh.numOfxNodes):

				node = self.mesh.getNodeByLoc(i,j)
				# First row. South BC
				if j == 0:
					temp = ((7*node.Vorticity - 8*node.north.Vorticity - 
							node.north.north.Vorticity)/(2.*self.mesh.dy**2) - 
							3.*uSouth/self.mesh.dy)
					node.Vorticity = temp
					node.StreamFunct = 0.0
					node.yVelocity = ySouth
					node.solved = True 

				# Last Row. North BC
				elif j == self.mesh.numOfyNodes-1:
					temp = ((7*node.Vorticity - 8*node.south.Vorticity - 
							node.south.south.Vorticity)/(2.*self.mesh.dy**2) - 
							3.*uNorth/self.mesh.dy)
					node.Vorticity = temp
					node.StreamFunct = 0.0
					node.yVelocity = yNorth
					node.solved = True

				# First column. West BC
				elif i == 0:
					temp = ((7*node.Vorticity - 8*node.east.Vorticity - 
							node.east.east.Vorticity)/(2.*self.mesh.dx**2) + 
							3.*vWest/self.mesh.dx)
					node.Vorticity = temp
					node.StreamFunct = 0.0
					node.xVelocity = vWest
					node.solved = True

				# Last column. East BC
				elif i == self.mesh.numOfxNodes-1:
					temp = ((7*node.Vorticity - 8*node.west.Vorticity - 
							node.west.west.Vorticity)/(2.*self.mesh.dx**2) + 
							3.*vEast/self.mesh.dx)
					node.Vorticity = temp
					node.StreamFunct = 0.0
					node.xVelocity = vEast
					node.solved = True
            
				else:
					pass

