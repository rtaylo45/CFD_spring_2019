"""
Author: Zack Taylor
Date: 1/18/19

The physics type applied problem specific physics to the mesh domain. This 
sets the solution.

"""
import copy
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

			ALap = self.LaplaceObj.getAMatrix()
			ALapCSC = sp.csc_matrix(ALap)

			while diff > -6.0:
				time = float(timeStep)*self.NavierObj.dt
				timeSteps.append(timeStep)
				bLap = self.LaplaceObj.getbVector()
			
				solLap = self.mesh.solveLinalg(ALap,bLap,A_=ALapCSC)
			
				self.__unPackSolution(solLap,'Lap')

				self.NavierObj.upDateBC()
				self.NavierObj.updateVelocities()

				ANav = self.NavierObj.getAMatrix()
				bNav = self.NavierObj.getbVector()
				lapDiffs.append(self.__calcDiff('Lap'))

				solNav = self.mesh.solveLinalg(ANav,bNav)
				self.__unPackSolution(solNav, "Vor")

				print time, log10(self.__calcDiff('Lap') + self.__calcDiff('Vor'))
				diff =  log10(self.__calcDiff('Lap') + self.__calcDiff('Vor'))
				print 
				navDiffs.append(self.__calcDiff('Vor'))
				timeStep += 1

		return timeSteps, lapDiffs, navDiffs
			
		#self.__exactLaplace()
		#self.__calcError()

	def __unPackSolution(self, solutionVector,var):
		for k in xrange(self.mesh.maxSolIndex):
			node = self.mesh.getNodeBySolIndex(k)
			if var == "Lap":
				node.oldLaplaceSolution = node.LaplaceSolution
				node.LaplaceSolution = solutionVector[k]
			if var =="Vor":
				node.oldVorticitySolution = node.VorticitySolution
				node.VorticitySolution = solutionVector[k]

	"""
	@Brief Loops over the mesh to set the error between exacpt and approx
	"""
	def __calcDiff(self,sol):
		diff = 0.0
		for i in xrange(self.mesh.numOfxNodes):
			for j in xrange(self.mesh.numOfyNodes):
				node = self.mesh.getNodeByLoc(i,j)
				if not node.solved:
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









