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
	def __init__(self, mesh, dt=0.01, Re=100.0):
		self.mesh = mesh
		self.dt = dt
		self.Re = Re
		self.CoeffA = None
		self.CoeffB = None
		self.CoeffC = None
		self.CoeffD = None
		self.CoeffE = None
		self.northBC = None
		self.southBC = None
		self.eastBC = None
		self.westBC = None

	"""
	@Brief Setter for boundary condition

	@param side     Face where BC is applied
	@param BC       The boundary condition value
	"""
	def setBC(self, side, BC, BCType=0):
		if side == "north":
			self.northBC = BC
		elif side == "south":
			self.southBC = BC
		elif side == "east":
			self.eastBC = BC
		elif side == "west":
			self.westBC = BC
		else:
			print "Invalid BC"

	"""
    @Brief Sets up the A matrix for a 5 point grid
	"""
	def getAMatrix(self):
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
		B = np.zeros((self.mesh.maxSolIndex,1))
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			self.__setACoefficients(node)
			B[k] = B[k] + node.VorticitySolution/self.dt
			
			if node.west.solved:
				B[k] = B[k] - node.west.VorticitySolution*self.CoeffD
			if node.east.solved:
				B[k] = B[k] - node.east.VorticitySolution*self.CoeffB
			if node.south.solved:
				B[k] = B[k] - node.south.VorticitySolution*self.CoeffE
			if node.north.solved:
				B[k] = B[k] - node.north.VorticitySolution*self.CoeffA
		return B 

	"""
	@brief Sets the coefficients for the A matrix
	"""
	def __setACoefficients(self, node):
		v = node.vVelocity
		u = node.uVelocity
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
	@Brief Loops through the model to update the velocity
	"""
	def updateVelocities(self):
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				
				if not node.solved:
					node.uVelocity = ((node.north.LaplaceSolution - 
					node.south.LaplaceSolution)/(2.*self.mesh.dy))

					node.vVelocity = -((node.west.LaplaceSolution -
					node.east.LaplaceSolution)/(2.*self.mesh.dx))

	"""
	@Brief Sets the boundary conditions for each face
	"""
	def upDateBC(self):
		vEast = self.eastBC
		vWest = self.westBC
		vNorth = self.northBC
		vSouth = self.southBC

		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
			
				# west BC	
				if i == 0:
					vort = ((7.*node.LaplaceSolution - 
					8.*node.east.LaplaceSolution + 
					node.east.east.LaplaceSolution)/(2.*self.mesh.dx**2.))	
					node.VorticitySolution = vort
					node.solved = True

				# east BC
				if i == (self.mesh.numOfxNodes-1):
					vort = ((7.*node.LaplaceSolution - 
					8.*node.west.LaplaceSolution + 
					node.west.west.LaplaceSolution)/(2.*self.mesh.dx**2.))	
					node.VorticitySolution = vort
					node.solved = True
						
				# south BC	
				if j == 0:
					vort = ((7.*node.LaplaceSolution - 
					8.*node.north.LaplaceSolution + 
					node.north.north.LaplaceSolution)/(2.*self.mesh.dy**2.))	
					node.VorticitySolution = vort
					node.solved = True

				# north BC
				if j == (self.mesh.numOfyNodes-1):
					vort = ((7.*node.LaplaceSolution - 
					8.*node.south.LaplaceSolution + 
					node.south.south.LaplaceSolution)/(2.*(self.mesh.dy**2.))
					- 3.*vNorth/(self.mesh.dy))	
					node.VorticitySolution = vort
					node.solved = True
















