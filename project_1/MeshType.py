"""
Author: Zack Taylor
Date: 1/15/19

Mesh class used to generate and house geometry data.
Node class used to house nodal information
"""

class Mesh(object):

    """
    @param xLength  Total length in x direction
    @param yLength  Total length in y direction
    @param xNodes   Number of nodes in x direction
    @param yNodes   Number of nodes in y direction
    """
    def __init__(self, xLength, yLength, xNodes, yNodes):
        self.xLength = xLength
        self.yLength = xLength
        self.numOfxNodes = xNodes
        self.numOfyNodes = yNodes
        self.numOfTotalNodes = xNodes*yNodes
        self.nodes = []
        self.dx = None
        self.dy = None

        # builds the geometry
        self.__buildGeo()
        

    def __buildGeo(self):

        self.dx = self.xLength/(self.numOfxNodes-1)
        self.dy = self.yLength/(self.numOfyNodes-1)

        self.__createNodes()
        self.__connectNodes()

    def __createNodes(self):

        absoluteIndex = 0
       
        # Loops from bottom left hand corner of mesh to top of mesh. Then goes
        # to next x index. 
        for i in xrange(self.numOfxNodes):
            for j in xrange(self.numOfyNodes):
                nodeObj = Node(i, j, absoluteIndex)
                self.nodes.append(nodeObj)
                absoluteIndex += 1

    def __connectNodes(self):
        for node in self.nodes:
            # First column. West BC
            if node.i == 0:
                node.approx = self.westBC
                node.solved = True

            # Last column. East BC
            elif node.i == self.numOfxNodes-1:
                node.approx = self.eastBC
                node.solved = True

            # First row. South BC
            elif node.j == 0:
                node.approx = self.southBC 
                node.solved = True 

            # Last Row. North BC
            elif node.j == self.numOfyNodes-1:
                node.approx = self.northBC
                node.solved = True

            else:
                print "Connection gone wrong"


    def setBC(side, BC):
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
                 
class Node(Mesh):

    """
    @param i    x index
    @param j    y index
    """
    def __init__(self, i, j, absIndex):
        self.i = i
        self.j = j
        self.absIndex = absIndex
        self.approx = 0.0
        self.exact = 0.0
        self.error = 0.0
        self.solved = False

        # Node connection information
        self.east = None
        self.west = None
        self.north = None
        self.south = None

    """
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





















