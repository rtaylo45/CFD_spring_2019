import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)

print "Setting BC"
# sets the boundary conditions
mesh.setBC(side="south", BC=0., BCType=0)
mesh.setBC(side="east", BC=0.0, BCType=0)
mesh.setBC(side="west", BC=0.0, BCType=0)
mesh.setBC(side="north", BC=100., BCType=0)

print "Finalizing mesh"
# Finalize the mesh
mesh.finalize()

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

print "Solving problem"
# solve the system
problem.solve(solveType=1)

print mesh.globalError

print "plotting solution"
#mesh.plot(solution="approx",plotType='2d', numOfPlotLines=900)








