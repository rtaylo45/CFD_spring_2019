from MeshType import Node, Mesh
import pytest

def test_initialization():
	xLength = 10.0
	yLength = 7.0
	xNodes = 11
	yNodes = 8
	mesh = Mesh(xLength,yLength,xNodes,yNodes)
	assert mesh.xLength == xLength
	assert mesh.yLength == yLength
	assert mesh.numOfxNodes == xNodes
	assert mesh.numOfyNodes == yNodes
	assert mesh.dx == xLength/(xNodes-1.0)
	assert mesh.dy == yLength/(yNodes-1.0)

def test_BC():
	mesh = Mesh(10.0,10.0,5,5)
	mesh.setBC(side="south", BC=0.)
	mesh.setBC(side="east", BC=25.0)
	mesh.setBC(side="west", BC=50.0)
	mesh.setBC(side="north", BC=100.)
	mesh.finalize()

 	for j in xrange(mesh.numOfyNodes):
		for i in xrange(mesh.numOfxNodes):

			node = mesh.getNodeByLoc(i,j)
			# First row. South BC
			if j == 0:
				assert node.solution == 0.0
				assert node.solved == True

			# Last Row. North BC
			elif j == mesh.numOfyNodes-1:
				assert node.solution == 100.0
				assert node.solved == True

			# First column. West BC
			elif i == 0:
				assert node.solution == 50.0
 				assert node.solved == True

			# Last column. East BC
 			elif i == mesh.numOfxNodes-1:
				assert node.solution == 25.0
				assert node.solved == True
            
			else:
				assert node.solution == None
				assert node.solved == False