import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log
import time
#start = time.time()

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=15.0, yLength=20.0, xNodes=141, yNodes=181)

print "Setting BC"
# sets the boundary conditions
mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=0.0)
mesh.setBC(side="west", BC=50.0)
mesh.setBC(side="north", BC=100.)

print "Finalizing mesh"
# Finalize the mesh
mesh.finalize()

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

print "Solving problem"
# solve the system
problem.solve(solveType=1)
#iterats = [x for x in xrange(1,iterations+1,1)]
#diffs = [log(dif) for dif in diffs]

print "plotting solution"
#mesh.plot2D(solution="exact")
#mesh.plot(solution="exact",plotType='2d', numOfPlotLines=300)
mesh.plot(solution="approx",plotType='2d', numOfPlotLines=900)
#end = time.time()
#print "Program runtime: ",(end - start)

#plt.title('Global Residual vs. iteration')
#plt.ylabel('log(Residual)')
#plt.xlabel('Interations')
#plt.plot(iterats, diffs,'o')
#plt.show()








