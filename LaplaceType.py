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
	@Brief Sets up the A matrix for a 5 point grid
	"""
	def getAMatrix(self):
		self.__setACoefficients()
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
		A = np.zeros((numOfUnknowns,numOfUnknowns))

		# diagonal
		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(i)
			
			if node.boundary:
				A[i,i] = 1.0
			else:
				# main diagonal
				A[i,i] = self.CoeffC

				# the off off diagonals
				A[i,i+numOfRows] = self.CoeffB
				A[i,i-numOfRows] = self.CoeffD

				# the off diagonals
				A[i,i+1] = self.CoeffA
				
				A[i,i-1] = self.CoeffE

		return A

	"""
	@Brief Builds the B part of super matrix
	"""
	def getBMatrix(self):
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
		B = np.zeros((numOfUnknowns,numOfUnknowns))
		
		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(i)
			if node.boundary:
				B[i,i] = 0.0
			else:
				B[i,i] = 1.0
		return B
	"""
	@Brief Builds and returns the b vector 
	"""
	def getbVector(self):
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
		b = np.zeros((numOfUnknowns,1))

		for k in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(k)
			b[k] = 0.0
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
		

