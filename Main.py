import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=1.0, yLength=1.0, xNodes=5, yNodes=5)

print "Setting BC"
# sets the boundary conditions
mesh.setBC(side="south", BC=0., BCType=0)
mesh.setBC(side="east", BC=0.0, BCType=0)
mesh.setBC(side="west", BC=0.0, BCType=0)
mesh.setBC(side="north", BC=1., BCType=0)

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

print "Solving problem"
# solve the system
problem.solve(solveType=1)

#print mesh.globalError

print "plotting solution"
mesh.plot(solution="Phi",plotType='3d', numOfPlotLines=900)
mesh.plot(solution="w",plotType='3d', numOfPlotLines=900)








