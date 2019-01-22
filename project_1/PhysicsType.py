"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
#from math import sin, sinh, pi 
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
		self.__sweep(pattern)
		self.__exactLaplace()
		self.__calcError()

	"""
	@Brief The part that actually solves the problem
	"""
	def __updateSolution(self, node):
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
			#self.__sweepTopDown()
			self.__sweepBottomUp()
		elif pattern==1:
			A = self.__getAMatrix()
			print A
			#b = self.__getbVector()

			#solutionVector = self.mesh.solveMatrix(A,b)

			#self.__unPackSolution(solutionVector)



	"""
	@Brief Sweeps from the top of the domain to the bottom
	"""
	def __sweepTopDown(self):
		for i in xrange(self.mesh.numOfxNodes-1,-1,-1):
			for j in xrange(self.mesh.numOfyNodes-1,-1,-1):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.solved:
					self.__updateSolution(node)

	def __sweepBottomUp(self):
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.solved:
					self.__updateSolution(node)


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
		ycontribution = (Tn+Tw)/(dy**2.)

		node.solution = (1./temp1)*(xcontribution+ycontribution)

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

	"""
	@Brief Creats the A matrix used in the solution by Ax=b method
	"""
	def __getAMatrix(self):
		numOfColumns = self.mesh.numOfxNodes-2
		numOfRows = self.mesh.numOfyNodes-2
		lastNode = self.mesh.maxSolIndex-1
		A = np.zeros((numOfColumns*numOfRows,numOfRows*numOfColumns))

		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)

			# the first column
			if k < numOfRows:
				# first node (left hand bottom corner)
				if k==0:
					#A[0,0] = node.myCoeff
					#A[0,1] = node.northCoeff
					#A[0,numOfRows] = node.eastCoeff
					A[0,0] = 1
					A[0,1] = 2
					A[0,numOfRows] = 3

				elif k==numOfRows-1:
					#A[k,k] = node.myCoeff
					#A[k,k-1] = node.southCoeff
					#A[k+1,k] = node.eastCoeff
					A[k,k] = 8
					A[k,k-1] = 9
					A[k,k+numOfRows] = 10
				# Nodes between corners
				else:
					#A[k,k] = node.myCoeff
					#A[k,k-1] = node.southCoeff
					#A[k,k+1] = node.northCoeff
					#A[k,k+numOfRows] = node.eastCoeff
					A[k,k] = 4
					A[k,k-1] = 5
					A[k,k+1] = 6
					A[k,k+numOfRows] = 7

			# The last column
			elif k >= (numOfColumns-1)*numOfRows:
				print "last column",k
				# right hand bottom corner
				if k==(numOfColumns-1)*numOfRows:
					#A[k,k] = node.myCoeff
					#A[k,k+1] = node.northCoeff
					#A[k-numOfRows,k] = node.westCoeff
					A[k,k] = 11
					A[k,k+1] = 12
					A[k,k-numOfRows] = 13
				# right and top corner
				elif k==lastNode:
					A[k,k] = 18
					A[k,k-1] = 19
					A[k,k-numOfRows] = 20
				else:
					A[k,k] = 14
					A[k,k+1] = 15
					A[k,k-1] = 16
					A[k,k-numOfRows] = 17
			else:
				print "mid ",k
		return A






















