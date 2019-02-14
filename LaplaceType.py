"""
Author: Zack Taylor
Date: 2/7/19

The Laplace package used by the physics driver
"""
import numpy as np

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
		self.__setACoefficients()
		self.__setNodeSource()

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
		b = np.zeros((self.mesh.maxSolIndex,1))
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			b[k] = node.LaplaceSource
		return b
			
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
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			node.LaplaceSource = -node.VorticitySolution

			if node.south.solved:
				node.LaplaceSource += node.south.LaplaceSolution*self.CoeffD
			if node.north.solved:
				node.LaplaceSource += node.north.LaplaceSolution*self.CoeffA
			if node.east.solved:
				node.LaplaceSource += node.east.LaplaceSolution*self.CoeffB
			if node.west.solved:
				node.LaplaceSource += node.west.LaplaceSolution*self.CoeffD
		
	"""
	@Brief Sets the boundary conditions for each face
	"""
	def __applyBC(self):
		for j in xrange(self.mesh.numOfyNodes):
			for i in xrange(self.mesh.numOfxNodes):

				node = self.mesh.getNodeByLoc(i,j)
				# First row. South BC
				if j == 0:
					node.LaplaceSolution = self.southBC 
					node.solved = True 

				# Last Row. North BC
				elif j == self.mesh.numOfyNodes-1:
					node.LaplaceSolution = self.northBC
					node.solved = True

				# First column. West BC
				elif i == 0:
					node.LaplaceSolution = self.westBC
					node.solved = True

				# Last column. East BC
				elif i == self.mesh.numOfxNodes-1:
					node.LaplaceSolution = self.eastBC
					node.solved = True
            
				else:
					pass

