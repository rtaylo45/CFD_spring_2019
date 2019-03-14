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
		self.CoeffF = None
		self.CoeffG = None
		self.CoeffH = None
		self.CoeffI = None
		self.northBC = None
		self.southBC = None
		self.eastBC = None
		self.westBC = None

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
	@ Brief Builds and return the C matrix
	"""
	def getCMatrix(self):
		numOfColumns = self.mesh.numOfxiNodes
		numOfRows = self.mesh.numOfetaNodes
		numOfUnknowns = numOfColumns*numOfRows
		C = np.zeros((numOfUnknowns,numOfUnknowns))

		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(i)

			coeffAy = -7.*node.beta/(2.*self.mesh.deta**2.)
			coeffBy = 4.*node.beta/self.mesh.deta**2.
			coeffCy = -1.*node.beta/(2.*self.mesh.deta**2.)
			coeffAx = -7.*node.alfa/(2.*self.mesh.dxi**2.)
			coeffBx = 4.*node.alfa/self.mesh.dxi**2.
			coeffCx = -1.*node.alfa/(2.*self.mesh.dxi**2.)
	
			if node.boundary:
				if node.boundaryLoc == 'north':
					# k node
					C[i,i] = coeffAy
					# south node
					C[i,i-1] = coeffBy
					# south south node
					C[i,i-2] = coeffCy	
				elif node.boundaryLoc == 'south':
					# k node
					C[i,i] = coeffAy
					# north node
					C[i,i+1] = coeffBy
					# north north node
					C[i,i+2] = coeffCy
				elif node.boundaryLoc =='west':
					# k node
					C[i,i] = coeffAx
					# east node
					C[i,i+numOfRows] = coeffBx
					# east east node
					C[i,i+2*numOfRows] = coeffCx
				elif node.boundaryLoc == 'east':
					# k node
					C[i,i] = coeffAx
					# west node
					C[i,i-numOfRows] = coeffBx
					# west west node
					C[i,i-2*numOfRows] = coeffCx
				else:
					pass
		return C

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

			if node.boundary:
				if node.boundaryLoc == 'north':
					b[k] = -(3.*node.beta/self.mesh.deta+node.Q)*1.0/node.detady
				else:
					b[k] = 0.0
			else:
				b[k] = node.VorticitySolution/self.dt
			
		return b

	"""
	@brief Sets the coefficients for the A matrix
	"""
	def __setACoefficients(self, node):
		Vg = node.vVelocity
		Ug = node.uVelocity
		deta = self.mesh.deta
		dxi = self.mesh.dxi
		alfa = node.alfa
		beta = node.beta
		gama = node.gama
		Q = node.Q
		P = node.P
		Re = self.Re
		dt = self.dt

		self.CoeffA = -Vg/2./deta - beta/Re/deta**2. - Q/2./Re/deta
		self.CoeffB = Ug/2./dxi - alfa/Re/dxi**2. - P/2./Re/dxi
		self.CoeffC = 1./dt + 2.*alfa/(Re*dxi**2.) + 2.*beta/(Re*deta**2.)
		self.CoeffD = -Ug/2./dxi - alfa/Re/dxi**2. + P/2./Re/dxi
		self.CoeffE = Vg/2./deta - beta/Re/deta**2. + Q/2./Re/deta
		self.CoeffF = -2.*gama/4./Re/dxi/deta
		self.CoeffG = 2.*gama/4./Re/dxi/deta
		self.CoeffH = 2.*gama/4./Re/dxi/deta
		self.CoeffI = -2.*gama/4./Re/dxi/deta

	"""
	@Brief Loops through the model to update the velocity
	"""
	def updateVelocities(self):
		for i in xrange(self.mesh.numOfxiNodes):
			for j in xrange(self.mesh.numOfetaNodes):
				node = self.mesh.getNodeByLoc(i,j)
				
				if not node.boundary:
					node.uVelocity = node.jac*((node.north.LaplaceSolution - 
					node.south.LaplaceSolution)/(2.*self.mesh.deta))

					node.vVelocity = node.jac*((node.west.LaplaceSolution -
					node.east.LaplaceSolution)/(2.*self.mesh.dxi))

