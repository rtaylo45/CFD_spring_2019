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
    @Brief Sets up the A matrix for a 5 point grid
	"""
	def getAMatrix(self):
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
				self.__setACoefficients(node)
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
	@ Brief Builds and return the C matrix
	"""
	def getCMatrix(self):
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
		C = np.zeros((numOfUnknowns,numOfUnknowns))

		coeffAy = -7./(2.*self.mesh.dy**2.)
		coeffBy = 4./self.mesh.dy**2.
		coeffCy = -1./(2.*self.mesh.dy**2.)
		coeffAx = -7./(2.*self.mesh.dx**2.)
		coeffBx = 4./self.mesh.dx**2.
		coeffCx = -1./(2.*self.mesh.dx**2.)

		for i in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(i)
	
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
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
		b = np.zeros((numOfUnknowns,1))

		for k in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(k)

			if node.boundary:
				if node.boundaryLoc == 'north':
					b[k] = -3.*1.0/self.mesh.dy
				else:
					b[k] = 0.0
			else:
				b[k] = node.VorticitySolution/self.dt
			
		return b

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
				
				if not node.boundary:
					node.uVelocity = ((node.north.LaplaceSolution - 
					node.south.LaplaceSolution)/(2.*self.mesh.dy))

					node.vVelocity = -((node.west.LaplaceSolution -
					node.east.LaplaceSolution)/(2.*self.mesh.dx))

