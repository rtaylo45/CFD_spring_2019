"""
Author: Zack Taylor
Date: 2/7/19

The Vorticity package used by the physics driver
"""
import numpy as np

class Vorticity(object):

	"""
	@brief Initializes the Vorticity object

	@param Mesh		The Mesh object
	"""
	def __init__(self, mesh, dt=1.0, Re=100.0):
		self.mesh = mesh
		self.dt = dt
		self.Re = Re
		self.CoeffA = None
		self.CoeffB = None
		self.CoeffC = None
		self.CoeffD = None
		self.CoeffE = None
        # BC on north side of mesh
		self.northBC = None
        # BC type on the north side of mesh
		self.northBCType = None
        # BC on south side of mesh
		self.southBC = None
        # BC type on the south side of mesh
		self.southBCType = None
        # BC on east side of mesh
		self.eastBC = None
        # BC type on the east side of mesh
		self.eastBCType = None
        # BC on west side of mesh
		self.westBC = None
        # BC type on the west side of mesh
		self.westBCType = None

	"""
	@Brief Setter for boundary condition

	@param side     Face where BC is applied
	@param BC       The boundary condition value
	@param BCType   The type of boundary condition being applied
	"""
	def setBC(self, side, BC, BCType=0):
		if side == "north":
			self.northBC = BC
			self.northBCType = BCType
		elif side == "south":
			self.southBC = BC
			self.southBCType = BCType
		elif side == "east":
			self.eastBC = BC
			self.eastBCType = BCType
		elif side == "west":
			self.westBC = BC
			self.westBCType = BCType
		else:
			print "Invalid BC"

	"""
	@brief runs the presolve to apply BC
	"""
	def runPresolve(self):
		self.__applyBC()
		self.__setNodeSource()

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
				A[i,i+numOfRows] = self.CoeffB
			if i+1 > numOfRows:
				A[i,i-numOfRows] = self.CoeffD

			# the off diagonals
			if (i+1)%numOfRows:
				A[i,i+1] = self.CoeffA
			
			if (i)%(numOfRows):
				A[i,i-1] = self.CoeffE

		return A

	"""
	@Brief Builds and returns the b vector 
	"""
	def getbVector(self):
		self.__applyBC()
		self.__setNodeSource()
		B = np.zeros((self.mesh.maxSolIndex,1))
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			B[k] = node.VorticitySource
		return B 

	"""
	@brief Sets the coefficients for the A matrix
	"""
	def __setACoefficients(self, node):
		v = node.yVelocity
		u = node.xVelocity
		dy = self.mesh.dy
		dx = self.mesh.dx
		Re = self.Re
		dt = self.dt

		self.CoeffA = v/(2.*dy) - 1./(Re*dy**2)
		self.CoeffB = u/(2.*dx) - 1./(Re*dx**2)
		self.CoeffC = 1./dt + 2./(Re*dx**2) + 2./(Re*dy**2)
		self.CoeffD = -u/(2.*dx) - 1./(Re*dx**2)
		self.CoeffE = -v/(2.*dy) - 1./(Re*dy**2)

	"""
	@Brief Sets the coeffients for the b matrix 
	"""
	def __setNodeSource(self):
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			self.__setACoefficients(node)
			node.VorticitySource = node.VorticitySolution/self.dt

			if node.south.solved:
				node.VorticitySource -= node.south.VorticitySolution*self.CoeffE
			if node.north.solved:
				node.VorticitySource -= node.north.VorticitySolution*self.CoeffA
			if node.east.solved:
				node.VorticitySource -= node.east.VorticitySolution*self.CoeffB
			if node.west.solved:
				node.VorticitySource -= node.west.VorticitySolution*self.CoeffD
	"""
	@Brief Loops through the model to update the velocity
	"""
	def __updateVelocities(self):
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			node.yVelocity = ((node.north.LaplaceSolution - 
			node.south.LaplaceSolution)/(2.*self.mesh.dy))

			node.xVelocity = (-(node.east.LaplaceSolution -
			node.west.LaplaceSolution)/(2.*self.mesh.dx))

	"""
	@Brief Sets the boundary conditions for each face
	"""
	def __applyBC(self):
		#vRight = self.eastBC
		#vLeft = self.westBC
		#uTop = self.northBC
		#uBottom = self.southBC

		vRight = 0.0
		vLeft = 0.0
		uTop = 1.0
		uBottom = 0.0

		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):

				node = self.mesh.getNodeByLoc(i,j)
				# Last Row. North BC
				if j == self.mesh.numOfyNodes-1:
					vort = ((7.*node.LaplaceSolution - 8.*node.south.LaplaceSolution
					+ node.south.south.LaplaceSolution)/(2.*self.mesh.dy**2) - 
					(3./self.mesh.dy)*uTop)
					node.VorticitySolution = vort
					node.solved = True

				# First row. South BC
				elif j == 0:
					vort = ((7.*node.LaplaceSolution - 8.*node.north.LaplaceSolution
					+ node.north.north.LaplaceSolution)/(2.*self.mesh.dy**2) + 
					(3./self.mesh.dy)*uBottom)
					node.VorticitySolution = vort
					node.solved = True 


				# First column. West BC
				elif i == 0:
					vort = ((7.*node.LaplaceSolution - 8.*node.east.LaplaceSolution
					+ node.east.east.LaplaceSolution)/(2.*self.mesh.dx**2) - 
					(3./self.mesh.dx)*vLeft)
					node.VorticitySolution = vort
					node.solved = True

				# Last column. East BC
				elif i == self.mesh.numOfxNodes-1:
					vort = ((7.*node.LaplaceSolution - 8.*node.west.LaplaceSolution
					+ node.west.west.LaplaceSolution)/(2.*self.mesh.dx**2) + 
					(3./self.mesh.dx)*vRight)
					node.VorticitySolution = vort
					node.solved = True
				else:
					pass































