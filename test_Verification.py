from MeshType import Node, Mesh
import PhysicsType as Phy
import pytest

def test_HeatDiffusionSmallGrid():

	# Test direct solve	
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=4, yNodes=5)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=1)

	assert mesh.globalError < 1.e-1

	# Test Gauss Seidel solve
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=4, yNodes=5)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=0)

	assert mesh.globalError < 1.e-1

def test_HeatDiffusionMidGrid():
	
	# Test direct solve	
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=1)

	assert mesh.globalError < 1.e-2

	# Test Gauss Seidel solve	
	# Builds mesh
	mesh = Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)
	
	# sets the boundary conditions
	mesh.setBC(side="south", BC=0., BCType=0)
	mesh.setBC(side="east", BC=0.0, BCType=0)
	mesh.setBC(side="west", BC=0.0, BCType=0)
	mesh.setBC(side="north", BC=100., BCType=0)
	
	# Generates the problem around the mesh
	problem = Phy.Physics(mesh=mesh, problemType="Laplace")
	
	# solve the system
	problem.solve(solveType=0)

	assert mesh.globalError < 1.e-2
