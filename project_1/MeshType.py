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
        self.Nodes = []
        self.dx = None
        self.dy = None

        # builds the geometry
        self.__buildGeo()
        

    def __self.buildGeo(self):

        self.dx = self.xLength/(self.numOfxNodes-1)
        self.dy = self.yLength/(self.numOfyNodes-1)

        self.__createNodes()
        self.__connectNodes()

    def __createNodes(self):
        pass

    def __connectNodes(self):
        pass


class Node(Mesh):

    """
    @param i    x index
    @param j    y index
    """
    def __init__(self, i, j)
        self.i = i
        self.j = j
        self.approx = 0.0
        self.exact = 0.0
        self.error = 0.0

        # Node connection information
        self.east = None
        self.west = None
        self.north = None
        self.south = None
