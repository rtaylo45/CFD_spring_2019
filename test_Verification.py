from MeshType import Node, Mesh
import PhysicsType as Phy
import pytest

def test_HeatDiffusionSmallGrid():
	
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

def test_HeatDiffusionMidGrid():
	
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
