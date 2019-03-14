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
		self.CoeffF = None
		self.CoeffG = None
		self.CoeffH = None
		self.CoeffI = None
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
		numOfColumns = self.mesh.numOfxiNodes
		numOfRows = self.mesh.numOfetaNodes
		numOfUnknowns = numOfColumns*numOfRows
		A = np.zeros((numOfUnknowns,numOfUnknowns))

		# diagonal
		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(i)
			
			if node.boundary:
				A[i,i] = 1.0
			else:
				self.__setACoefficients(node)
				# main diagonal k node
				A[i,i] = self.CoeffC

				# the off off diagonals
				# east node
				A[i,i+numOfRows] = self.CoeffB
				# west node
				A[i,i-numOfRows] = self.CoeffD
				# south west node
				A[i,i-numOfRows-1] = self.CoeffI
				# north west node
				A[i,i-numOfRows+1] = self.CoeffG
				# south east node
				A[i,i+numOfRows-1] = self.CoeffH
				# north east node
				A[i,i-numOfRows+1] = self.CoeffF

				# the off diagonals
				# north node
				A[i,i+1] = self.CoeffA
				# south node	
				A[i,i-1] = self.CoeffE

		return A

	"""
	@Brief Builds the B part of super matrix
	"""
	def getBMatrix(self):
		numOfColumns = self.mesh.numOfxiNodes
		numOfRows = self.mesh.numOfetaNodes
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
		numOfColumns = self.mesh.numOfxiNodes
		numOfRows = self.mesh.numOfetaNodes
		numOfUnknowns = numOfColumns*numOfRows
		b = np.zeros((numOfUnknowns,1))

		for k in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(k)
			b[k] = 0.0
		return b
			
	"""
	@brief Sets the coefficients for the A matrix
	
	@param node		node object
	"""
	def __setACoefficients(self, node):
		self.CoeffA = node.beta/self.mesh.deta**2. + node.Q/2./self.mesh.dxi
		self.CoeffB = node.alfa/self.mesh.dxi**2. + node.P/2./self.mesh.dxi
		self.CoeffC = -2.*node.alfa/self.mesh.dxi**2. - 2.*node.beta/self.mesh.deta**2.
		self.CoeffD = node.alfa/self.mesh.dxi**2. - node.P/2./self.mesh.dxi
		self.CoeffE = node.beta/self.mesh.deta**2. - node.Q/2./self.mesh.dxi
		self.CoeffF = 2.*node.gama/4./self.mesh.dxi/self.mesh.deta
		self.CoeffG = -2.*node.gama/4./self.mesh.dxi/self.mesh.deta
		self.CoeffH = -2.*node.gama/4./self.mesh.dxi/self.mesh.deta
		self.CoeffI = 2.*node.gama/4./self.mesh.dxi/self.mesh.deta
