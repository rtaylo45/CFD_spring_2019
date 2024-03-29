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

	@param xiLength  Total length in xi direction
	@param etaLength  Total length in eta direction
	@param xiNodes   Number of nodes in xi direction
	@param etaNodes   Number of nodes in eta direction
	"""
	def __init__(self,xLength, yLength, xiLength, etaLength, xiNodes, etaNodes):

		# Total length in xi direction
		self.xiLength = xiLength
        # Total length in eta direction
		self.etaLength = etaLength
		# Total length in x direction
		self.xLength = xLength
        # Total length in y direction
		self.yLength = yLength
        # Total number of nodes in xi direction
		self.numOfxiNodes = xiNodes
        # Total number of nodes in eta direction
		self.numOfetaNodes = etaNodes
        # Total number of mesh nodes
		self.numOfTotalNodes = xiNodes*etaNodes
        # List of all the mesh nodes
		self.nodes = []
        # List of all the nodes in the solution domain. This doesn't include
        # the boundary nodes
		self.solNodes = []
        # The delta x term
		self.dx = None
        # The delta y term
		self.dy = None
		# the delta xi term for generalized coordinate system
		self.dxi = None
		# the delta eta term for generalized coordinate system
		self.deta = None
        # The number of solution nodes
		self.maxSolIndex = None
		# The Global error
		self.globalError = None 

		# builds the geometry
		self.__buildGeo()
		self.__setSolutionIndex()
		self.__setBoundaries()
		self.__setDerivatives()
		self.__calcCrossDerivatives()
		self.__calcJacobian()
		self.__calcTransformMetrics()
        
	"""
    @Brief Builds the Gemoetry and sets dx and dy
    """
	def __buildGeo(self):

		# sets the comp domain
		self.dxi = self.xiLength/(self.numOfxiNodes-1.0)
		self.deta = self.etaLength/(self.numOfetaNodes-1.0)

		# sets the physical domain
		self.dx = self.xLength/(self.numOfxiNodes-1.0)
		self.dy = self.yLength/(self.numOfetaNodes-1.0)

		self.__createNodes()
		self.__connectNodes()

	"""
	@brief Builds the nodes and puts them in mesh data continer
	"""
	def __createNodes(self):

		absoluteIndex = 0
		alen = np.zeros((self.numOfxiNodes,1))
		for i in xrange(1,self.numOfxiNodes):
			alen[i] = alen[i-1] + min(float(i),self.numOfxiNodes-i)**1.0
		alen = alen/alen[self.numOfxiNodes-1]
       
		# Loops from bottom left hand corner of mesh to top of mesh. Then goes
		# to next x index. 
		for i in xrange(self.numOfxiNodes):
			for j in xrange(self.numOfetaNodes):
				xi = float(i)*self.dxi
				eta = float(j)*self.deta
				x = float(alen[i])
				y = float(alen[j])
				#x = float(i)*self.dx
				#y = float(j)*self.dy

				nodeObj = Node(iPhy=i, jPhy=j, iComp=i, jComp=j, absIndex=absoluteIndex, 
				xi=xi, eta=eta, x=x, y=y)
				self.nodes.append(nodeObj)
				absoluteIndex += 1

	"""
	@Brief Makes the node connections
	"""
	def __connectNodes(self):
		for i in xrange(self.numOfxiNodes):
			for j in xrange(self.numOfetaNodes):

				node = self.getNodeByLoc(i,j)

				node.north = self.getNodeByLoc(i,j+1)
				node.south = self.getNodeByLoc(i,j-1)
				node.west = self.getNodeByLoc(i-1,j)
				node.east = self.getNodeByLoc(i+1,j)

	"""
	@Brief Checks the i,j index to make sure they are in range

	@param i    xi index
	@param j    eta index
	"""
	def __checkij(self,i,j):
        # Check to make sure i index is in range
		if i < 0:
			return False
		elif i >= self.numOfxiNodes:
			return False
		# Check to see if j index is in range
		elif j < 0:
			return False
		elif j >= self.numOfetaNodes:
			return False
		else:
			return True

	"""
	@Brief Sets the solution index
	"""
	def __setSolutionIndex(self):
		k = 0 
		for i in xrange(1,self.numOfxiNodes-1):
			for j in xrange(1,self.numOfetaNodes-1):
				node = self.getNodeByLoc(i,j)
				node.solIndex = k
				self.solNodes.append(node)	
				k+=1
		self.maxSolIndex = k

	"""
	@Brief Sets the boundaries
	"""
	def __setBoundaries(self):
		for j in xrange(self.numOfetaNodes):
			# west BC	
			i = 0
			node = self.getNodeByLoc(i,j)
			node.boundary = True
			node.boundaryLoc = 'west'

			# east BC
			i = (self.numOfxiNodes-1)
			node = self.getNodeByLoc(i,j)
			node.boundary = True
			node.boundaryLoc = 'east'
						
		for i in xrange(self.numOfxiNodes):
			# south BC	
			j = 0
			node = self.getNodeByLoc(i,j)
			node.boundary = True
			node.boundaryLoc = 'south'

			# north BC
			j = (self.numOfetaNodes-1)
			node = self.getNodeByLoc(i,j)
			node.boundary = True
			node.boundaryLoc = 'north'

	"""
	@Brief Calculates the first and second derivatives for physical to 
		computational mesh
	"""
	def __setDerivatives(self):
		
		# The first and second derivatives for inner points and boundaries
		for i in xrange(self.numOfxiNodes):
			for j in xrange(self.numOfetaNodes):
				node = self.getNodeByLoc(i,j)

				if not node.boundary:
					# first derivatives
					node.dxdxi = (node.east.x - node.west.x)/2./self.dxi
					node.dydxi = (node.east.y - node.west.y)/2./self.dxi
					node.dxdeta = (node.north.x - node.south.x)/2./self.deta
					node.dydeta = (node.north.y - node.south.y)/2./self.deta
				
					# second derivatives
					node.ddxdxidxi = (node.east.x - 2.*node.x + 
						node.west.x)/self.dxi**2.
					node.ddydxidxi = (node.east.y - 2.*node.y + 
						node.west.y)/self.dxi**2.
					node.ddxdetadeta = (node.north.x - 2.*node.x + 
						node.south.x)/self.deta**2.
					node.ddydetadeta = (node.north.y - 2.*node.y + 
						node.south.y)/self.deta**2.
					node.ddxdxideta = (node.north.east.x - node.south.east.x - 
						node.north.west.x + node.south.west.x)/4./self.dxi/self.deta 
					node.ddydxideta = (node.north.east.y - node.south.east.y -
						node.north.west.y + node.south.west.y)/4./self.dxi/self.deta

				else:
					# Left (west) boundary 
					if i == 0:
						# xi derivatives
						node.dxdxi = (-3.*node.x + 4.*node.east.x -
							node.east.east.x)/2./self.dxi
						node.dydxi = (-3.*node.y + 4.*node.east.y -
							node.east.east.y)/2./self.dxi
						node.ddxdxidxi = (2.*node.x - 5.*node.east.x +
							4.*node.east.east.x - node.east.east.east.x)/self.dxi**2.
						node.ddydxidxi = (2.*node.y - 5.*node.east.y +
							4.*node.east.east.y - node.east.east.east.y)/self.dxi**2.
						# eta derivatives
						if j > 0 and j < self.numOfetaNodes-1:
							node.dxdeta = (node.north.x - node.south.x)/2./self.deta	
							node.dydeta = (node.north.y - node.south.y)/2./self.deta	
							node.ddxdetadeta = (node.north.x - 2.*node.x + 
								node.south.x)/self.deta**2.
							node.ddydetadeta = (node.north.y - 2.*node.y + 
								node.south.y)/self.deta**2.
					# right (east) boundary
					if i == self.numOfxiNodes-1:
						# xi derivatives
						node.dxdxi = (3.*node.x - 4.*node.west.x +
							node.west.west.x)/2./self.dxi
						node.dydxi = (3.*node.y - 4.*node.west.y +
							node.west.west.y)/2./self.dxi
						node.ddxdxidxi = (2.*node.x - 5.*node.west.x +
							4.*node.west.west.x - node.west.west.west.x)/self.dxi**2.
						node.ddydxidxi = (2.*node.y - 5.*node.west.y +
							4.*node.west.west.y - node.west.west.west.y)/self.dxi**2.
						# eta derivatives
						if j > 0 and j < self.numOfetaNodes-1:
							node.dxdeta = (node.north.x - node.south.x)/2./self.deta	
							node.dydeta = (node.north.y - node.south.y)/2./self.deta	
							node.ddxdetadeta = (node.north.x - 2.*node.x + 
								node.south.x)/self.deta**2.
							node.ddydetadeta = (node.north.y - 2.*node.y + 
								node.south.y)/self.deta**2.
					# the bottom (south) boundary	
					if j == 0:
						# xi derivatives
						if i > 0 and i < self.numOfxiNodes-1:
							node.dxdxi = (node.east.x - node.west.x)/2./self.dxi	
							node.dydxi = (node.east.y - node.west.y)/2./self.dxi
							node.ddxdxidxi = (node.east.x - 2.*node.x + 
								node.west.x)/self.dxi**2.
							node.ddydxidxi = (node.east.y - 2.*node.y + 
								node.west.y)/self.dxi**2.
						# eta derivatives
						node.dxdeta = (-3.*node.x + 4.*node.north.x - 
							node.north.north.x)/2./self.deta
						node.dydeta = (-3.*node.y + 4.*node.north.y - 
							node.north.north.y)/2./self.deta
						node.ddxdetadeta = (2.*node.x - 5.*node.north.x + 
							4.*node.north.north.x - node.north.north.north.x)/self.deta**2.
						node.ddydetadeta = (2.*node.y - 5.*node.north.y + 
							4.*node.north.north.y - node.north.north.north.y)/self.deta**2.
					# top (north) boundary
					if j == self.numOfetaNodes-1:
						# eta derivatives
						node.dxdeta = (3.*node.x - 4.*node.south.x + 
							node.south.south.x)/2./self.deta
						node.dydeta = (3.*node.y - 4.*node.south.y + 
							node.south.south.y)/2./self.deta
						node.ddxdetadeta = (2.*node.x - 5.*node.south.x + 
							4.*node.south.south.x - node.south.south.south.x)/self.deta**2.
						node.ddydetadeta = (2.*node.y - 5.*node.south.y + 
							4.*node.south.south.y - node.south.south.south.y)/self.deta**2.
						# xi derivatives
						if i > 0 and i < self.numOfxiNodes-1:
							node.dxdxi = (node.east.x - node.west.x)/2./self.dxi	
							node.dydxi = (node.east.y - node.west.y)/2./self.dxi
							node.ddxdxidxi = (node.east.x - 2.*node.x + 
								node.west.x)/self.dxi**2.
							node.ddydxidxi = (node.east.y - 2.*node.y + 
								node.west.y)/self.dxi**2.
		
	"""
	@Brief computes the cross derivatives 
	"""
	def __calcCrossDerivatives(self):
		# cross derivates on the top and bottom excluding i = 0 and i = imax	
		for i in xrange(1, self.numOfxiNodes-1):
			for j in xrange(0, self.numOfetaNodes, self.numOfetaNodes-1):
				node = self.getNodeByLoc(i,j)

				node.ddxdxideta = (node.east.dxdeta - node.west.dxdeta)/2./self.dxi
				node.ddydxideta = (node.east.dydeta - node.west.dydeta)/2./self.dxi

		# cross derivatives at i = 0 and i = imax excluding the top and bottom 
		for i in xrange(0, self.numOfxiNodes, self.numOfxiNodes-1):
			for j in xrange(1, self.numOfetaNodes-1):
				node = self.getNodeByLoc(i,j)

				node.ddxdxideta = (node.north.dxdxi - node.south.dxdxi)/2./self.deta
				node.ddydxideta = (node.north.dydxi - node.south.dydxi)/2./self.deta

		# cross derivatives at the for corners
		node = self.getNodeByLoc(0,0)
		node.ddxdxideta = (-3.*node.dxdxi + 4.*node.north.dxdxi - 
			node.north.north.dxdxi)/2./self.deta
		node.ddydxideta = (-3.*node.dydxi + 4.*node.north.dydxi - 
			node.north.north.dydxi)/2./self.deta

		node = self.getNodeByLoc(0,self.numOfetaNodes-1)
		node.ddxdxideta = (3.*node.dxdxi - 4.*node.south.dxdxi + 
			node.south.south.dxdxi)/2./self.deta
		node.ddydxideta = (3.*node.dydxi - 4.*node.south.dydxi + 
			node.south.south.dydxi)/2./self.deta

		node = self.getNodeByLoc(self.numOfxiNodes-1,0)
		node.ddxdxideta = (-3.*node.dxdxi + 4.*node.north.dxdxi - 
			node.north.north.dxdxi)/2./self.deta
		node.ddydxideta = (-3.*node.dydxi + 4.*node.north.dydxi - 
			node.north.north.dydxi)/2./self.deta

		node = self.getNodeByLoc(self.numOfxiNodes-1, self.numOfetaNodes-1)
		node.ddxdxideta = (3.*node.dxdxi - 4.*node.south.dxdxi + 
			node.south.south.dxdxi)/2./self.deta
		node.ddydxideta = (3.*node.dydxi - 4.*node.south.dydxi + 
			node.south.south.dydxi)/2./self.deta
		
				
	"""
	@Brief computes the jacobian for the transformation
	"""
	def __calcJacobian(self):
		for i in xrange(self.numOfxiNodes):
			for j in xrange(self.numOfetaNodes):
				node = self.getNodeByLoc(i,j)
				
				node.jac = 1./(node.dxdxi*node.dydeta - node.dxdeta*node.dydxi)			

	"""
	@Brief computes the transform metric for the computational coor.
	"""
	def __calcTransformMetrics(self):
		for i in xrange(self.numOfxiNodes):
			for j in xrange(self.numOfetaNodes):
				node = self.getNodeByLoc(i,j)
				
				node.dxidx = node.jac*node.dydeta
				node.dxidy = -node.jac*node.dxdeta
				node.detadx = -node.jac*node.dydxi
				node.detady = node.jac*node.dxdxi

				node.alfa = node.dxidx**2. + node.dxidy**2.
				node.beta = node.detadx**2. + node.detady**2.
				node.gama = node.dxidx*node.detadx + node.dxidy*node.detady

				node.P = -(1.*node.alfa*(node.ddxdxidxi*node.dxidx + 
					node.ddydxidxi*node.dxidy) + 2.*node.gama*(node.ddxdxideta*
					node.dxidx + node.ddydxideta*node.dxidy) + 1.*node.beta*(
					node.ddxdetadeta*node.dxidx + node.ddydetadeta*node.dxidy))

				node.Q = -(1.*node.alfa*(node.ddxdxidxi*node.detadx + 
					node.ddydxidxi*node.detady) + 2.*node.gama*(node.ddxdxideta*
					node.detadx + node.ddydxideta*node.detady) + 1.*node.beta*(
					node.ddxdetadeta*node.detadx + node.ddydetadeta*node.detady))

	"""
	@Brief Returns node object 

	@param i    x index
	@param j    y index
	"""
	def getNodeByLoc(self, i, j):

		if self.__checkij(i,j):
			k = j + i*self.numOfetaNodes
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
		for i in xrange(self.numOfxiNodes):
			node = self.getNodeByLoc(i,0)
			x.append(node.x)

		for j in xrange(self.numOfetaNodes-1,-1,-1):
			node = self.getNodeByLoc(0,j)
			y.append(node.y)


		X, Y = np.meshgrid(x, y)

		Solution = np.zeros(X.shape)
		vVelocity = np.zeros(X.shape)
		uVelocity = np.zeros(X.shape)

		h = 0
		for j in xrange(self.numOfetaNodes-1,-1,-1):
			k = 0
			for i in xrange(0,self.numOfxiNodes,1):
				node = self.getNodeByLoc(i,j)
				if solution=="Phi":
					Solution[h,k] = node.LaplaceSolution
					solTitle = '2-D Stream Function. Resolution '
					saveTitle = 'StreamFunctionRe'+str(int(Re))+'Resolution'
					levels = [-1.0e-10,-1.0e-7,-1.0e-5,-1.0e-4,-0.01,
					-0.03,-0.05,-0.07,-0.09,-0.1,-0.11,-0.115,-0.1175,
					1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-4,2.5e-4,5.0e-4,
					1.0e-3,1.5e-3,3.0e-3]
				elif solution=="exact":
					Solution[h,k] = node.exact
				elif solution=="w":
					Solution[h,k] = node.VorticitySolution
					solTitle = '2-D Vorticity. Resolution '
					saveTitle = 'VorticityRe'+str(int(Re))+'Resolution'
					levels = [-3.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0]
				elif solution=="velocity":
					solTitle = '2-D Velocity vector. Resolution '
					saveTitle = 'VelocityRe'+str(int(Re))+'Resolution'
					vVelocity[h,k] = node.vVelocity
					uVelocity[h,k] = node.uVelocity	
				k+= 1
			h+=1
        
		if plotType=='2d':

			#cp = plt.contour(X, Y, Solution, sorted(levels), cmap='tab20')
			cp = plt.contour(X, Y, Solution,levels=sorted(levels), cmap='tab20')
			plt.title(solTitle + 'Re = '+ str(Re)+' '  +str(self.numOfxiNodes)+' x '+str(self.numOfetaNodes))
			plt.ylabel('eta')
			plt.xlabel('xi')
			plt.colorbar(cp)
			plt.savefig(saveTitle +str(self.numOfxiNodes)+'x'+str(self.numOfetaNodes)+'.png')
			#plt.show()
			plt.close()

		elif plotType=='3d':

			fig = plt.figure()
			ax = fig.gca(projection='3d')
			ax.plot_surface(X,Y,Solution,cmap=plt.cm.plasma,linewidth=0, antialiased=False)
			ax.set_xlabel('x (cm)')
			ax.set_ylabel('y (cm)')
			ax.set_zlabel('Temperature (c)')
			plt.title('2-D Laplace heat equation. Resolution '+str(self.numOfxiNodes)+' x '+str(self.numOfetaNodes))

			plt.show()
		

		elif plotType=='vector':
			fig = plt.figure()
			plt.title(solTitle + 'Re = '+ str(Re)+' '  +str(self.numOfxiNodes)+' x '+str(self.numOfetaNodes))
			Q = plt.quiver(X, Y, uVelocity, vVelocity,cmap='tab20')
			qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$1$', labelpos='E',
                   coordinates='figure')
			plt.ylabel('eta')
			plt.xlabel('xi')
			plt.savefig(saveTitle +str(self.numOfxiNodes)+'x'+str(self.numOfetaNodes)+'.png')
			#plt.show()
			plt.close()

	def plotCenterLineUVelocity(self):
		# center i node
		ic = (self.numOfxiNodes-1)/2
		plotVel = []
		ploty = []
		goalVel = [0.0,-0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
			-0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151,
			0.68717, 0.73722, 0.78871, 0.84123, 1.0]
		goaly = [0.0,0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 
			0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 
			0.9766, 1.0]
		for j in xrange(self.numOfetaNodes):
			node = self.getNodeByLoc(ic,j)
			plotVel.append(node.uVelocityPhy)
			ploty.append(node.y)
		plt.plot(plotVel,ploty, label="Solution")
		plt.plot(goalVel, goaly, 'o', label="Ghia")
		plt.grid()
		plt.legend()
		plt.ylabel("y")
		plt.ylim(0,1)
		plt.xlabel("U Velocity")
		plt.show()
		#plt.savefig("uVelocity.png")
		#plt.close()

	def plotCenterLineVVelocity(self):
		# center i node
		jc = (self.numOfetaNodes-1)/2
		plotVel = []
		plotx = []
		goalVel = [0.0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 
			0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914,
			-0.10313, -0.08864, -0.07391, -0.05906, 0.0]
		goalx = [0.0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
			0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 
			0.9688, 1.0]
		for i in xrange(self.numOfxiNodes):
			node = self.getNodeByLoc(i,jc)
			plotVel.append(-node.vVelocityPhy)
			plotx.append(node.x)
		plt.plot(plotx,plotVel, label="Solution")
		plt.plot(goalx, goalVel, 'o',label="Ghia")
		plt.legend()
		plt.grid()
		plt.ylabel("V Velocity")
		plt.xlabel("x")
		plt.xlim(0,1)
		#plt.savefig("vVelocity.png")
		plt.show()
		#plt.close()


	"""
	@Brief Linear alg solver. Solves Ax=b
	
	@param A	A matrix size nxn
	@param b	b vector size nx1
	"""
	def solveLinalg(self,A,b,A_=None):
		start = time.time()
		#if A_==None:
			#A_ = sp.csc_matrix(A)
		A = sp.csc_matrix(A)
		x = spla.spsolve(A,b)
		#x = np.linalg.solve(A,b)
		end = time.time()
		print end-start
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
	def __init__(self, iPhy=None, jPhy=None, absIndex=None, 
			iComp=None, jComp=None, x=None, y=None, xi=None, eta=None):
		# Physical Mesh i index, x index
		self.iPhy = iPhy
		# Physical Mesh j index, y index
		self.jPhy = jPhy
		# Computational Mesh i index, xi index
		self.iComp = iComp
		# Computational Mesh j index, eta index
		self.jComp = jComp
		# Absolution index. Includes the boundary nodes
		self.absIndex = absIndex
		# Solution index. Does not include boundary nodes
		self.solIndex = None
		# x position 
		self.x = x
		# y position
		self.y = y
		# xi position
		self.xi = xi
		# eta position
		self.eta = eta
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
		# the velocity in x direction computational
		self.uVelocity = 0.0
		# the velocity in the y direction computational
		self.vVelocity = 0.0
		# the velocity in x direction physical
		self.uVelocityPhy = 0.0
		# the velocity in the y direction physical
		self.vVelocityPhy = 0.0
		# the old solution
		self.oldSolution = None
		# the exact solution
		self.exact = None
		# error between approx and exact
		self.error = None
        # logical saying if the nodes has been solved or not
		self.boundary = False
		# Boundary location
		self.boundaryLoc = None
		# first derivatives of generalized coordinate system
		self.dxdxi = None
		self.dxdeta = None
		self.dydxi = None
		self.dydeta = None
		self.dxidx = None
		self.dxidy = None
		self.detadx = None
		self.detady = None
		# second derivatives of generalized coordinate system
		self.ddxdxidxi = None
		self.ddxdetadeta = None
		self.ddxdxideta = None
		self.ddydxidxi = None
		self.ddydetadeta = None
		self.ddydxideta = None
		# jacobian used  to transform physical domain to comp domain
		self.jac = None
		# transform metrics as defined by Daniel J. Garmann
		self.alfa = None
		self.beta = None
		self.gama = None
		self.P = None
		self.Q = None

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

