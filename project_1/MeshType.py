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
        self.xLength = xLength
        self.yLength = yLength
        self.numOfxNodes = xNodes
        self.numOfyNodes = yNodes
        self.numOfTotalNodes = xNodes*yNodes
        self.nodes = []
        self.dx = None
        self.dy = None
        self.northBC = None
        self.southBC = None
        self.eastBC = None
        self.westBC = None

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

    def __mapMeshijToMatrixij(self,i,j):
        k = i
        h = j+(self.numOfyNodes-1)-i
        return h,k

    """
    @Brief Runs functions prier to solve
    """
    def finalize(self):
        self.__applyBC()

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

    def plot2D(self):
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
        print X.shape
        print 
        print Solution

        for i in xrange(self.numOfxNodes):
            for j in xrange(self.numOfyNodes):
                h,k = self.__mapMeshijToMatrixij(i,j)
                print i,j, h,k
                node = self.getNodeByLoc(i,j)
                Solution[h,k] = node.solution

        #plt.contour(X, Y, Solution, 500, cmap='RdGy')
        #plt.title('Gibbs Energy at '+str(MolU)+' Mole U and '+str(MolTh)+' Mole Th')
        #plt.ylabel('Moles of O2')
        #plt.xlabel('Temperature (K)')
        #plt.colorbar()
        #plt.show()

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
        self.i = i
        self.j = j
        self.absIndex = absIndex
        self.x = x
        self.y = y
        self.solution = None
        self.exact = None
        self.error = None
        self.solved = False

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





















