import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=1.0, yLength=1.0, xNodes=20, yNodes=20)


print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

print "Setting BC"
# sets the boundary conditions
problem.LaplaceObj.setBC(side="north", BC=0., BCType=0)
problem.LaplaceObj.setBC(side="south", BC=0., BCType=0)
problem.LaplaceObj.setBC(side="east", BC=0., BCType=0)
problem.LaplaceObj.setBC(side="west", BC=0., BCType=0)

problem.NavierObj.setBC(side="north", BC=1., BCType=0)
problem.NavierObj.setBC(side="south", BC=0., BCType=0)
problem.NavierObj.setBC(side="east", BC=0., BCType=0)
problem.NavierObj.setBC(side="west", BC=0., BCType=0)

print "Solving problem"
# solve the system
problem.solve(solveType=1)

#print mesh.globalError

print "plotting solution"
mesh.plot(solution="Phi",plotType='2d', numOfPlotLines=30)
mesh.plot(solution="w",plotType='2d', numOfPlotLines=30)








