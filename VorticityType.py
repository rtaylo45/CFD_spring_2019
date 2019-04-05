"""
Author: Zack Taylor
Date: 2/7/19

The Vorticity package used by the physics driver
"""
import numpy as np
import matplotlib.pylab as plt

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

		for j in xrange(self.mesh.numOfxiNodes):
			for i in xrange(self.mesh.numOfetaNodes):
				node = self.mesh.getNodeByLoc(i,j)
			
				if node.boundary:
					A[node.absIndex,node.absIndex] = 1.0
				else:
					self.__setACoefficients(node)
					# main diagonal k node
					A[node.absIndex,node.absIndex] = self.CoeffC
					#print 'abs',node.absIndex, node.north.absIndex

					# the off off diagonals
					# east node
					A[node.absIndex,node.east.absIndex] = self.CoeffB
					# west node
					A[node.absIndex,node.west.absIndex] = self.CoeffD
					# south west node
					A[node.absIndex,node.south.west.absIndex] = self.CoeffI
					# north west node
					A[node.absIndex,node.north.west.absIndex] = self.CoeffG
					# south east node
					A[node.absIndex,node.south.east.absIndex] = self.CoeffH
					# north east node
					A[node.absIndex,node.north.east.absIndex] = self.CoeffF

					# the off diagonals
					# north node
					A[node.absIndex,node.north.absIndex] = self.CoeffA
					# south node	
					A[node.absIndex,node.south.absIndex] = self.CoeffE

		return A

	"""
	@ Brief Builds and return the C matrix
	"""
	def getCMatrix(self):
		numOfColumns = self.mesh.numOfxiNodes
		numOfRows = self.mesh.numOfetaNodes
		numOfUnknowns = numOfColumns*numOfRows
		C = np.zeros((numOfUnknowns,numOfUnknowns))

		for j in xrange(self.mesh.numOfxiNodes):
			for i in xrange(self.mesh.numOfetaNodes):
				node = self.mesh.getNodeByLoc(i,j)

				coeffAy = -7.*node.beta/(2.*self.mesh.deta**2.)
				coeffBy = 4.*node.beta/self.mesh.deta**2.
				coeffCy = -1.*node.beta/(2.*self.mesh.deta**2.)
	
				if node.boundary:
					if node.boundaryLoc == 'south':
						# k node
						C[node.absIndex,node.absIndex] = coeffAy
						# north node
						C[node.absIndex,node.north.absIndex] = coeffBy
						# north north node
						C[node.absIndex,node.north.north.absIndex] = coeffCy	
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

		self.CoeffA = Vg/2./deta - beta/Re/deta**2. - Q/2./Re/deta
		self.CoeffB = Ug/2./dxi - alfa/Re/dxi**2. - P/2./Re/dxi
		self.CoeffC = 1./dt + 2.*alfa/(Re*dxi**2.) + 2.*beta/(Re*deta**2.)
		self.CoeffD = -Ug/2./dxi - alfa/Re/dxi**2. + P/2./Re/dxi
		self.CoeffE = -Vg/2./deta - beta/Re/deta**2. + Q/2./Re/deta
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
					# Updates the computational velocity
					node.uVelocity = node.jac*((node.north.LaplaceSolution - 
					node.south.LaplaceSolution)/(2.*self.mesh.deta))

					node.vVelocity = -node.jac*((node.east.LaplaceSolution -
					node.west.LaplaceSolution)/(2.*self.mesh.dxi))

					# Updates the physical velocity
					node.uVelocityPhy = (((node.north.LaplaceSolution -
					node.south.LaplaceSolution)/(2.*self.mesh.deta)*node.detady) + 
					(node.east.LaplaceSolution - node.west.LaplaceSolution)/
					(2.*self.mesh.dxi)*node.dxidy)
	
					node.vVelocityPhy = -(((node.north.LaplaceSolution -
					node.south.LaplaceSolution)/(2.*self.mesh.deta)*node.detadx) +
					(node.east.LaplaceSolution - node.west.LaplaceSolution)/
					(2.*self.mesh.dxi)*node.dxidx)


