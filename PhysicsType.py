"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
import numpy as np
import LaplaceType as LaType

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
		self.LaplaceObj = LaType.Laplace(mesh)

		self.__runPreSolve()

	"""
    @brief Sets the initial guess for the solution
    """
	def __runPreSolve(self):
	
		self.LaplaceObj.runPresolve()

        # Solution method is iterative using initial guesses for temperature
		if self.solveType==0:
			for i in xrange(self.mesh.numOfxNodes):
				for j in xrange(self.mesh.numOfyNodes):
					node = self.mesh.getNodeByLoc(i,j)
					initialGuess = 0.0
					if not node.solved:
						node.solution = initialGuess
						node.Vorticity = initialGuess
						node.StreamFunct = initialGuess
						node.xVelocity = initialGuess
						node.yVelocity = initialGuess
	"""
	@Brief Solves the problem
    """
	def solve(self,solveType):
		if solveType==0:
			diffs, iterations = self.__gaussSeidel()
			return diffs, iterations
		elif solveType==1:
			A = self.LaplaceObj.getAMatrix()
			b = self.LaplaceObj.getbVector()

			solutionVector = self.mesh.solveLinalg(A,b)

			self.__unPackSolution(solutionVector)
		self.__exactLaplace()
		self.__calcError()

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
					temp = np.sinh(k*np.pi*y/w)
					temp1 = np.sin(k*np.pi*x/w)
					temp2 = 1./(k*np.sinh(k*np.pi*h/w))
					exact += temp*temp1*temp2
				exact = exact*400./np.pi
				node.exact = exact

	"""
	@Brief Loops over the mesh to set the error between exacpt and approx
	"""
	def __calcError(self):
		diff = 0.0
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.solved:
					diff += (node.exact-node.solution)**2.
		diff = diff**(.5)/(self.mesh.numOfyNodes*self.mesh.numOfxNodes)
		self.mesh.globalError = float(diff)

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
				if not node.solved:
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









