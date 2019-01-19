"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""

class Physics(object):

	"""
    @brief Initializes the Physics object

    @param Mesh 			The mesh object
    @param problem type		The type of problem we are trying to solve

	"""
	def __init__(self, mesh, problemType, solveType=0):
		self.mesh = mesh
		self.problemType = problemType
		self.solveType = solveType

		self.__runPreSolve()

	def __runPreSolve(self):

		# Solution method is iterative using initial guesses for temperature
		if self.solveType==0:
			for i in xrange(self.mesh.numOfxNodes):
				for j in xrange(self.mesh.numOfyNodes):
					node = self.mesh.getNodeByLoc(i,j)
					initialGuess = 0.0
					if not node.solved:
						node.solution = initialGuess

	def __updateSolution(self):
		pass

	def __sweep(self):
		pass
