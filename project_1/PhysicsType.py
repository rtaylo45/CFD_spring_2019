"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
from math import sin, sinh, pi 

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

	def solve(self):
		self.__sweep(0)
		self.__exactLaplace()
		self.__calcError()

	"""
	@Brief The part that actually solves the problem
	"""
	def __updateSolution(self, node):
		if self.problemType=="Laplace":
			realDiff = 1.0
			while realDiff > self.tol:
				oldSolution = copy.deepcopy(node.solution)
				# Updates the solution
				self.__updateLaplace(node)
				newSolution = node.solution
				realDiff = abs(newSolution-oldSolution)
				self.iterations += 1




	"""
	@Brief Loops through the mesh

	@pattern 	The pattern it uses to to sweep
	"""
	def __sweep(self,pattern):
		if pattern==0:
			self.__sweepTopDown()


	"""
	@Brief Sweeps from the top of the domain to the bottom
	"""
	def __sweepTopDown(self):
		for i in xrange(self.mesh.numOfxNodes-1,-1,-1):
			for j in xrange(self.mesh.numOfyNodes-1,-1,-1):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.solved:
					self.__updateSolution(node)



	def __updateLaplace(self,node):
		dx = self.mesh.dx
		dy = self.mesh.dy
		Te = node.east.solution
		Tw = node.west.solution
		Tn = node.north.solution
		Ts = node.south.solution

		temp1 = (dx**2.*dy**2.)/(2.*(dy**2.+dx**2.))
		xcontribution = (Te+Tw)/(dx**2.)
		ycontribution = (Tn+Tw)/(dy**2.)

		node.solution = temp1*(xcontribution+ycontribution)

	def __exactLaplace(self):
		h = self.mesh.yLength
		w = self.mesh.xLength
		n = 21

		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				x = node.x
				y = node.y
				exact = 0.0
				for k in xrange(1,n+1,2):
					top = sinh(n*pi*y/w)*sin(n*pi*x/w)
					bottom = 1./(n*sinh(n*pi*h/w))
					exact = exact + top/bottom

				node.exact = (400.0/pi)*exact


	def __calcError(self):
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				try:
					absError = abs(node.exact-node.solution)/node.solution
					node.error = absError
				except:
					node.error = 0.0























