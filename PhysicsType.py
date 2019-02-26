"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
import matplotlib.pylab as plt
import numpy as np
import scipy.sparse.linalg as spla
from scipy import sparse as sp
import LaplaceType as LaType
import VorticityType as VortType
from math import log10

class Physics(object):

	"""
    @brief Initializes the Physics object

    @param Mesh 			The mesh object
    @param problem type		The type of problem we are trying to solve
    """
	def __init__(self, mesh, dt, Re, solveType=0, tol=1.0e-2):
		self.mesh = mesh
		self.solveType = solveType
		self.tol = tol
		self.iterations = 0
		self.LaplaceObj = LaType.Laplace(mesh)
		self.NavierObj = VortType.Vorticity(mesh, dt, Re)

	"""
	@Brief Solves the problem
    """
	def solve(self,solveType):

		if solveType==0:
			diffs, iterations = self.__gaussSeidel()
			return diffs, iterations

		elif solveType==1:
			timeSteps = []
			lapDiffs = []
			navDiffs = []
			timeStep = 0
			time = 0.0
			diff = 0.0

			while diff > -10.0:
				time = float(timeStep)*self.NavierObj.dt
				timeSteps.append(timeStep)
				
				# builds the super matrix
				S, BLap, ALap, DNav, CNav = self.__buildSuperMatrix()
				# builds the super b vector
				b = self.__buildSuperbVector()
				#print b
				# solves the system
				sol = self.mesh.solveLinalg(S,b)
			
				self.__unPackSolution(sol)

				print time, log10(self.__calcDiff('Lap') + self.__calcDiff('Vor'))
				diff =  log10(self.__calcDiff('Lap') + self.__calcDiff('Vor'))
				print 
				navDiffs.append(self.__calcDiff('Vor'))
				timeStep += 1

		return S,b, BLap, ALap, DNav, CNav
			
	def __unPackSolution(self, solutionVector):
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows

		for k in xrange(numOfUnknowns):
			node = self.mesh.getNodeByAbsIndex(k)
			# Update Laplace solution
			node.oldLaplaceSolution = node.LaplaceSolution
			node.LaplaceSolution = solutionVector[k]
			# Update Vorticity solution
			node.oldVorticitySolution = node.VorticitySolution
			node.VorticitySolution = solutionVector[k+numOfUnknowns]

	"""
	@Brief Loops over the mesh to set the error between exacpt and approx
	"""
	def __calcDiff(self,sol):
		diff = 0.0
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.boundary:
					if sol == 'Lap':
						diff += (node.oldLaplaceSolution-
						node.LaplaceSolution)**2.
					elif sol =='Vor':
						diff += (node.oldVorticitySolution -
						node.VorticitySolution)**2.
		diff = diff**(.5)/(self.mesh.numOfyNodes*self.mesh.numOfxNodes)
		return float(diff)

	"""
	@Brief Runs the Gauss-Seidel solver
	"""
	def __gaussSeidel(self):
		diff = 1.0
		diffs = []
		iterations = 0

		while diff > self.tol:
			# Loops through the model
			self.__sweepTopDown()

			diff = 0.0

			for k in xrange(self.mesh.maxSolIndex):
				node = self.mesh.getNodeBySolIndex(k)
				diff = diff + abs(node.oldSolution-node.solution)
			
			diffs.append(diff)
			iterations += 1
		return diffs, iterations

	"""
	@Brief Sweeps from the top of the domain to the bottom
	"""
	def __sweepTopDown(self):
		for i in xrange(self.mesh.numOfxNodes-1,-1,-1):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				node.oldSolution = copy.deepcopy(node.solution)
				if not node.boundary:
					self.__updateLaplace(node)

	"""
	@Brief Updates the solution for a sweeping pattern soltuion method

	@param node 	mesh node object
	"""
	def __updateLaplace(self,node):
		dx = self.mesh.dx
		dy = self.mesh.dy
		Te = node.east.solution
		Tw = node.west.solution
		Tn = node.north.solution
		Ts = node.south.solution

		temp1 = (2./dx**2. + 2./dy**2.)
		xcontribution = (Te+Tw)/(dx**2.)
		ycontribution = (Tn+Ts)/(dy**2.)

		node.solution = (1./temp1)*(xcontribution+ycontribution)


	"""
	@Brief Builds the big matrix used to for unsegregated solver
	|A B| = |e|
	|C D| = |f|
	"""
	def __buildSuperMatrix(self):
		numOfColumns = self.mesh.numOfxNodes
		numOfRows = self.mesh.numOfyNodes
		numOfUnknowns = numOfColumns*numOfRows
        
		# update velocities for D matrix
		self.NavierObj.updateVelocities()
	
		# builds the B matrix
		BLap = self.LaplaceObj.getBMatrix()

		# builds the A martix	
		ALap = self.LaplaceObj.getAMatrix()
		#plt.spy(ALap)
		#plt.show()

		# builds |A B| part of matrix
		topHalf = np.concatenate((ALap,BLap), axis=1)
		#plt.spy(topHalf)
		#plt.show()


		# builds the D matrix
		DNav = self.NavierObj.getAMatrix()
		#plt.spy(DNav)
		#plt.show()

		# builds the C matrix
		CNav = self.NavierObj.getCMatrix()
		#plt.spy(CNav)
		#plt.show()

		# builds |C D| part of matrix
		botHalf = np.concatenate((CNav,DNav),axis=1)
		#plt.spy(botHalf)
		#plt.show()

		# builds the super matrix
		SMatrix = np.concatenate((topHalf,botHalf))
		#plt.spy(SMatrix)
		#plt.show()

		return SMatrix, BLap, ALap, DNav, CNav
		

	"""
	@Brief builds the super b vector 
	|A B| = |e|
	|C D| = |f|
	"""
	def __buildSuperbVector(self):
		# builds the e vector
		e = self.LaplaceObj.getbVector()
		
		# builds the f vector
		f = self.NavierObj.getbVector()

		return np.concatenate((e,f))
	
