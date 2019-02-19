"""
Author: Zack Taylor
Date: 1/15/19

Mesh class used to generate and house geometry data.
Node class used to house nodal information

Note on the Mesh Class: The general convenction is to loop from x = 0, y to ymax.
and so on. In the apply BC method its difference to get the 100 BC to include the
corners of the domain. 

"""
import numpy as np
import time
import scipy.sparse.linalg as spla
from scipy import sparse as sp
import matplotlib.pyplot as plt
from math import log10
from mpl_toolkits.mplot3d import Axes3D

class Mesh(object):

	"""
	@brief Initializes the Mesh object. Serves as a data structure.

	@param xLength  Total length in x direction
	@param yLength  Total length in y direction
	@param xNodes   Number of nodes in x direction
	@param yNodes   Number of nodes in y direction
	"""
	def __init__(self, xLength, yLength, xNodes, yNodes):

		# Total length in x direction
		self.xLength = xLength
        # Total length in y direction
		self.yLength = yLength
        # Total number of nodes in x direction
		self.numOfxNodes = xNodes
        # Total number of nodes in y direction
		self.numOfyNodes = yNodes
        # Total number of mesh nodes
		self.numOfTotalNodes = xNodes*yNodes
        # List of all the mesh nodes
		self.nodes = []
        # List of all the nodes in the solution domain. This doesn't include
        # the boundary nodes
		self.solNodes = []
        # The delta x term
		self.dx = None
        # The delta y term
		self.dy = None
        # BC on north side of mesh
		self.northBC = None
        # BC type on the north side of mesh
		self.northBCType = None
        # BC on south side of mesh
		self.southBC = None
        # BC type on the south side of mesh
		self.southBCType = None
        # BC on east side of mesh
		self.eastBC = None
        # BC type on the east side of mesh
		self.eastBCType = None
        # BC on west side of mesh
		self.westBC = None
        # BC type on the west side of mesh
		self.westBCType = None
        # Containter that stores the map from mesh to ploting matrix
		self.jMeshToMatrix = None
        # The number of solution nodes
		self.maxSolIndex = None
		# The Global error
		self.globalError = None 

		# builds the geometry
		self.__buildGeo()
		self.__setSolutionIndex()
        
	"""
    @Brief Builds the Gemoetry and sets dx and dy
    """
	def __buildGeo(self):

		self.dx = self.xLength/(self.numOfxNodes-1.0)
		self.dy = self.yLength/(self.numOfyNodes-1.0)

		self.__createNodes()
		self.__connectNodes()

	"""
	@brief Builds the nodes and puts them in mesh data continer
	"""
	def __createNodes(self):

		absoluteIndex = 0
       
		# Loops from bottom left hand corner of mesh to top of mesh. Then goes
		# to next x index. 
		for i in xrange(self.numOfxNodes):
			for j in xrange(self.numOfyNodes):
				x = i*self.dx
				y = j*self.dy
				nodeObj = Node(i, j, absoluteIndex, x, y)
				self.nodes.append(nodeObj)
				absoluteIndex += 1

	"""
	@Brief Makes the node connections
	"""
	def __connectNodes(self):
		for i in xrange(self.numOfxNodes):
			for j in xrange(self.numOfyNodes):

				node = self.getNodeByLoc(i,j)

				node.north = self.getNodeByLoc(i,j+1)
				node.south = self.getNodeByLoc(i,j-1)
				node.west = self.getNodeByLoc(i-1,j)
				node.east = self.getNodeByLoc(i+1,j)

	"""
	@Brief Checks the i,j index to make sure they are in range

	@param i    x index
	@param j    y index
	"""
	def __checkij(self,i,j):
        # Check to make sure i index is in range
		if i < 0:
			return False
		elif i >= self.numOfxNodes:
			return False
		# Check to see if j index is in range
		elif j < 0:
			return False
		elif j >= self.numOfyNodes:
			return False
		else:
			return True

	"""
	@Brief Sets the solution index
	"""
	def __setSolutionIndex(self):
		k = 0 
		for i in xrange(1,self.numOfxNodes-1):
			for j in xrange(1,self.numOfyNodes-1):
				node = self.getNodeByLoc(i,j)
				node.solIndex = k
				self.solNodes.append(node)	
				k+=1
		self.maxSolIndex = k

	"""
	@Brief Setter for boundary condition

	@param side     Face where BC is applied
	@param BC       The boundary condition value
	@param BCType   The type of boundary condition being applied
	"""
	def setBC(self, side, BC, BCType=0):
		if side == "north":
			self.northBC = BC
			self.northBCType = BCType
		elif side == "south":
			self.southBC = BC
			self.southBCType = BCType
		elif side == "east":
			self.eastBC = BC
			self.eastBCType = BCType
		elif side == "west":
			self.westBC = BC
			self.westBCType = BCType
		else:
			print "Invalid BC"

	"""
	@Brief Returns node object 

	@param i    x index
	@param j    y index
	"""
	def getNodeByLoc(self, i, j):

		if self.__checkij(i,j):
			k = j + i*self.numOfyNodes
			node = self.nodes[k]
			return node
		else:
			return None

	"""
	@Brief Gets the node using the absolute index

	@param k    The absolute index
	"""
	def getNodeByAbsIndex(self,k):
		if k<0:
			return None
		else:
			return self.nodes[k]

	"""
	@Brief Gets the node using the solution index. Solution index doesn't 
		   induce the boundary conditions

    @param k    The solution index
    """
	def getNodeBySolIndex(self,k):
		if k<0:
			return None
		elif k>self.maxSolIndex-1:
			return None
		else:
			return self.solNodes[k]


	"""
	@Brief Plots the solution on a contour plot

	@param solution         The soltuion you are trying to plot. 
                            exact or approx
	@param plotType         The plot type either 2d or 3d
	@param numOfPlotLines   Number of lines in contour plots
	"""
	def plot(self,solution="approx",plotType='2d', dt=0.0, Re=0.0):
		x = []
		y = []
		for i in xrange(self.numOfxNodes):
			node = self.getNodeByLoc(i,0)
			x.append(node.x)

		for j in xrange(self.numOfyNodes-1,-1,-1):
			node = self.getNodeByLoc(0,j)
			y.append(node.y)


		X, Y = np.meshgrid(x, y)

		Solution = np.zeros(X.shape)

		h = 0
		for j in xrange(self.numOfyNodes-1,-1,-1):
			k = 0
			for i in xrange(0,self.numOfxNodes,1):
				node = self.getNodeByLoc(i,j)
				if solution=="Phi":
					Solution[h,k] = node.LaplaceSolution
					solTitle = '2-D Stream Function. Resolution '
					saveTitle = 'StreamFunctionRe'+str(Re)+'dt'+str(dt) +'Resolution'
					levels = [-1.0e-10,-1.0e-7,-1.0e-5,-1.0e-4,-0.01,
					-0.03,-0.05,-0.07,-0.09,-0.1,-0.11,-0.115,-0.1175,
					1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-4,2.5e-4,5.0e-4,
					1.0e-3,1.5e-3,3.0e-3]
				elif solution=="exact":
					Solution[h,k] = node.exact
				elif solution=="w":
					Solution[h,k] = node.VorticitySolution
					solTitle = '2-D Vorticity. Resolution '
					saveTitle = 'VorticityRe'+str(Re)+'dt'+str(dt) +'Resolution'
					levels = [-3.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0]
				k+= 1
			h+=1
        
		if plotType=='2d':

			cp = plt.contour(X, Y, Solution, sorted(levels), cmap='tab20')
			plt.title(solTitle +str(self.numOfxNodes)+' x '+str(self.numOfyNodes))
			plt.ylabel('y')
			plt.xlabel('x')
			plt.colorbar(cp)
			plt.savefig(saveTitle +str(self.numOfxNodes)+'x'+str(self.numOfyNodes)+'.png')
			plt.close()
			#plt.show()

		elif plotType=='3d':

			fig = plt.figure()
			ax = fig.gca(projection='3d')
			ax.plot_surface(X,Y,Solution,cmap=plt.cm.plasma,linewidth=0, antialiased=False)
			ax.set_xlabel('x (cm)')
			ax.set_ylabel('y (cm)')
			ax.set_zlabel('Temperature (c)')
			plt.title('2-D Laplace heat equation. Resolution '+str(self.numOfxNodes)+' x '+str(self.numOfyNodes))

			plt.show()

	"""
	@Brief Linear alg solver. Solves Ax=b
	
	@param A	A matrix size nxn
	@param b	b vector size nx1
	"""
	def solveLinalg(self,A,b,A_=None):
		if A_==None:
			A_ = sp.csc_matrix(A)
		x, exitCode = spla.gcrotmk(A_,b)
		end = time.time()
		return x

class Node(Mesh):

	"""
	@Brief Initilizes the node class

	@param i        x index
	@param j        y index
	@param absIndex Abosolute index
	@param x        x location on domain
	@param y        y location on domain
	"""
	def __init__(self, i, j, absIndex, x, y):
		# Mesh i index, x index
		self.i = i
		# Mesh j index, y index
		self.j = j
		# Absolution index. Includes the boundary nodes
		self.absIndex = absIndex
		# Solution index. Does not include boundary nodes
		self.solIndex = None
		# x position 
		self.x = x
		# y position
		self.y = y
		# the approx solution
		self.solution = None
		# the vorticity solution
		self.VorticitySolution = 0.0
		# the stream function solution
		self.LaplaceSolution = 0.0
		# the vorticity solution
		self.oldVorticitySolution = 0.0
		# the stream function solution
		self.oldLaplaceSolution = 0.0
		# the velocity in x direction
		self.uVelocity = 0.0
		# the velocity in the y direction
		self.vVelocity = 0.0
		# the old solution
		self.oldSolution = None
		# the exact solution
		self.exact = None
		# error between approx and exact
		self.error = None
        # logical saying if the nodes has been solved or not
		self.solved = False
		# General Source 
		self.source = 0.0
		# Stream function source
		self.LaplacSource = 0.0
		# Vorticity source
		self.VorticitySource = 0.0

        # Node connection information
		self.east = None
		self.west = None
		self.north = None
		self.south = None

	"""
	@Brief Sets the node connection

	@param location Location of node
	@param node     Node object

	Location index:
                    north = 0
                    south = 1
                    east  = 2
                    west  = 3
	"""
	def connect(self,location,node):
		if location == 0:
			self.north = node
		elif location == 1:
			self.south = node
		elif location == 2:
			self.east = node
		elif location == 3:
			self.west = node
		else:
			print "Invalid location"

