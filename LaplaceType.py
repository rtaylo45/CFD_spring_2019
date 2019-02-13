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
		pass
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
		pass	
	"""
	@Brief Sets the boundary conditions for each face
	"""
	def __applyBC(self):
		for j in xrange(self.mesh.numOfyNodes):
			for i in xrange(self.mesh.numOfxNodes):

				node = self.mesh.getNodeByLoc(i,j)
				# First row. South BC
				if j == 0:
					node.solution = self.mesh.southBC 
					node.solved = True 

				# Last Row. North BC
				elif j == self.mesh.numOfyNodes-1:
					node.solution = self.mesh.northBC
					node.solved = True

				# First column. West BC
				elif i == 0:
					node.solution = self.mesh.westBC
					node.solved = True

				# Last column. East BC
				elif i == self.mesh.numOfxNodes-1:
					node.solution = self.mesh.eastBC
					node.solved = True
            
				else:
					pass

