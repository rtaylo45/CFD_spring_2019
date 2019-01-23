"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
import numpy as np

class Physics(object):

	"""
    @brief Initializes the Physics object

    @param Mesh 			The mesh object
    @param problem type		The type of problem we are trying to solve

	"""
	def __init__(self, mesh, problemType, solveType=0, tol=1.0e-2):
		self.mesh = mesh
		self.problemType = problemType
		self.solveType = solveType
		self.tol = tol
		self.iterations = 0

		self.__runPreSolve()

	"""
	@brief Sets the initial guess for the solution
	"""
	def __runPreSolve(self):

		# Solution method is iterative using initial guesses for temperature
		if self.solveType==0:
			for i in xrange(self.mesh.numOfxNodes):
				for j in xrange(self.mesh.numOfyNodes):
					node = self.mesh.getNodeByLoc(i,j)
					initialGuess = 0.0
					if not node.solved:
						node.solution = initialGuess
	"""
	@Brief Solves the problem
	"""
	def solve(self,pattern):
		if pattern==0:
			self.__gaussSeidel()
		elif pattern==1:
			A = self.__getAMatrix()
			b = self.__getbVector()

			solutionVector = self.mesh.solveLinalg(A,b)

			self.__unPackSolution(solutionVector)
		self.__exactLaplace()

	def __unPackSolution(self, solutionVector):
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			node.solution = solutionVector[k]


	"""
	@Brief Loops over the mesh and sets the exact soltuion
	"""
	def __exactLaplace(self):
		h = self.mesh.yLength
		w = self.mesh.xLength
		n = 101

		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				x = node.x
				y = node.y
				exact = 0.0

				for k in xrange(1,n+1,2):
					temp = np.sinh(k*np.pi*y/h)
					temp1 = np.sin(k*np.pi*x/w)
					temp2 = 1./(k*np.sinh(k*np.pi*h/w))
					exact += temp*temp1*temp2
				exact = exact*400./np.pi
				node.exact = exact

	"""
	@Brief Loops over the mesh to set the error between exacpt and approx
	"""
	def __calcError(self):
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				try:
					absError = abs(node.exact-node.solution)/node.solution
					node.error = absError
				except:
					node.error = 0.0

	def __getAMatrix(self):
		numOfColumns = self.mesh.numOfxNodes-2
		numOfRows = self.mesh.numOfyNodes-2
		A = np.zeros((numOfColumns*numOfRows,numOfRows*numOfColumns))

		alpha = 1./(self.mesh.dx**2.)
		beta = -(2./(self.mesh.dx**2.) + 2./(self.mesh.dy**2.))
		gamma = 1./(self.mesh.dy**2.)

		# diagonal
		for i in xrange(numOfColumns*numOfRows):
			# main diagonal
			A[i,i] = beta
			# the off off diagonals
			if i-numOfRows<0:
				pass
			else:
				A[i-numOfRows,i] = alpha
				A[i,i-numOfRows] = alpha

			# the off diagonals
			if (i+1)%numOfRows:
				A[i+1,i] = gamma
				A[i,i+1] = gamma

		return A

	def __getbVector(self):
		B = np.zeros((self.mesh.maxSolIndex,1))
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			B[k] = node.source

		return B 













	def __gaussSeidel(self):
		diff = 1.0
		iterations = 0

		while diff > self.tol:
			# Loops through the model
			if iterations % 2 == 0:
				self.__sweepTopDown()
			else:
				self.__sweepBottomUp()

			diff = 0.0

			for node in self.mesh.nodes:
				if not node.solved:
					diff = diff + abs(node.solution-node.oldSolution)

			iterations += 1

	"""
	@Brief Sweeps from the top of the domain to the bottom
	"""
	def __sweepTopDown(self):
		for i in xrange(self.mesh.numOfxNodes-1,-1,-1):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				node.oldSolution = copy.deepcopy(node.solution)
				if not node.solved:
					self.__updateLaplace(node)

	def __sweepBottomUp(self):
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				node.oldSolution = copy.deepcopy(node.solution)
				if not node.solved:
					self.__updateLaplace(node)

	def __sweepAbsIndex(self):
		for node in self.mesh.nodes:
			if not node.solved:
				self.__updateSolution(node)
	def __sweepAbsIndexRev(self):
		pass

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









