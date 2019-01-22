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
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
plt.style.use('seaborn-white')

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
        # BC on south side of mesh
        self.southBC = None
        # BC on east side of mesh
        self.eastBC = None
        # BC on west side of mesh
        self.westBC = None
        # Containter that stores the map from mesh to ploting matrix
        self.jMeshToMatrix = None
        # The number of solution nodes
        self.maxSolIndex = None

        # builds the geometry
        self.__buildGeo()
        
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
    @Brief Sets the boundary conditions for each face
    """
    def __applyBC(self):
        for j in xrange(self.numOfyNodes):
            for i in xrange(self.numOfxNodes):

                node = self.getNodeByLoc(i,j)
                # First row. South BC
                if j == 0:
                    node.solution = self.southBC 
                    node.solved = True 

                # Last Row. North BC
                elif j == self.numOfyNodes-1:
                    node.solution = self.northBC
                    node.solved = True

                # First column. West BC
                elif i == 0:
                    node.solution = self.westBC
                    node.solved = True

                # Last column. East BC
                elif i == self.numOfxNodes-1:
                    node.solution = self.eastBC
                    node.solved = True
            
                else:
                    pass

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
    @Brief Takes in the mesh i j index and transforms it to the plotting matrix 
           i,j (h,k).

    @param i    The i mesh index
    @param j    The j mesh index
    """
    def __mapMeshijToMatrixij(self,i,j):
        k = i
        h = self.jMeshToMatrix[j]
        return h,k

    """
    @Brief Generates a map for connecting the solution to the plotting matrix.
           This is done because Array for plot starts from top. The way the 
           mesh is created is from the bottom up
    """
    def __generateMapToMatrix(self):

        # j index vector
        jVect = [j for j in xrange(self.numOfyNodes)]
        self.jMeshToMatrix = jVect[::-1]

    """
    @Brief Sets the solution index for each node
    """
    def __setSolIndex(self):
        k = 0
        for i in xrange(self.numOfxNodes):
            for j in xrange(self.numOfyNodes):
                node = self.getNodeByLoc(i,j)
                if not node.solved:
                    node.solIndex = k
                    self.solNodes.append(node)
                    k += 1
        self.maxSolIndex = k

    """
    @Brief Sets the node coefficients for the A matrix
    """
    def __setAMatrixCoeff(self):
        alpha = 1./self.dx**2.
        beta = -(2./self.dx**2. + 2./self.dy**2.)
        gamma = 1./self.dy**2.
        for i in xrange(self.numOfxNodes):
            for j in xrange(self.numOfxNodes):
                node = self.getNodeByLoc(i,j)

                node.myCoeff = beta
                node.eastCoeff = alpha
                node.westCoeff = alpha
                node.northCoeff = gamma
                node.southCoeff = gamma

    """
    @Brief Runs functions prier to solve
    """
    def finalize(self):
        self.__applyBC()
        self.__generateMapToMatrix()
        self.__setSolIndex()
        self.__setAMatrixCoeff()

    """
    @Brief Setter for boundary condition

    @param side     Face where BC is applied
    @param BC       The boundary condition value
    """
    def setBC(self, side, BC):
        if side == "north":
            self.northBC = BC
        elif side == "south":
            self.southBC = BC
        elif side == "east":
            self.eastBC = BC
        elif side == "west":
            self.westBC = BC
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

    @param solution     The soltuion you are trying to plot. 
                        exact or approx
    """
    def plot2D(self,solution):
        x = []
        y = []
        for i in xrange(self.numOfxNodes):
            node = self.getNodeByLoc(i,0)
            x.append(node.x)

        for j in xrange(self.numOfyNodes):
            node = self.getNodeByLoc(0,j)
            y.append(node.y)


        X, Y = np.meshgrid(x, y)

        Solution = np.zeros(X.shape)

        for i in xrange(self.numOfxNodes):
            for j in xrange(self.numOfyNodes):
                h,k = self.__mapMeshijToMatrixij(i,j)
                node = self.getNodeByLoc(i,j)
                if solution=="approx":
                    Solution[h,k] = node.solution
                elif solution=="exact":
                    Solution[h,k] = node.exact

        plt.contour(X, Y, Solution,500, cmap='RdGy')
        #plt.title('Gibbs Energy at '+str(MolU)+' Mole U and '+str(MolTh)+' Mole Th')
        #plt.ylabel('Moles of O2')
        #plt.xlabel('Temperature (K)')
        plt.colorbar()
        #ax = sns.heatmap(Solution)
        #sns.plt.show()
        plt.show()

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
        # the exact solution
        self.exact = None
        # error between approx and exact
        self.error = None
        # logical saying if the nodes has been solved or not
        self.solved = False

        # Node connection information
        self.east = None
        self.west = None
        self.north = None
        self.south = None

        # Coefficients for A matrix
        self.eastCoeff = None
        self.westCoeff = None
        self.northCoeff = None
        self.southCoeff = None
        self.myCoeff = None


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





















