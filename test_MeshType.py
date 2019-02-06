from MeshType import Node, Mesh
import PhysicsType as Phy
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

def test_HeadDiffusionSmallGrid():
	
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=4, yNodes=5)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Finalize the mesh
	mesh.finalize()
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=1)

	assert pytest.approx(mesh.globalError) == 0.08223148

def test_HeadDiffusionMidGrid():
	
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Finalize the mesh
	mesh.finalize()
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=1)

	assert pytest.approx(mesh.globalError) == 0.00167193766785
