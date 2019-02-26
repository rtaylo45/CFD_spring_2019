import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log10

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=1.0, yLength=1.0, xNodes=51, yNodes=51)
dt = 0.1
Re = 500.0

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, dt=dt,Re=Re)

print "Solving problem"
# solve the system
timeSteps, Diffs = problem.solve(solveType=1)


plt.plot(timeSteps,Diffs, 'o')
plt.title('Global Residual '+ str(mesh.numOfxNodes)+ ' x '+str(mesh.numOfyNodes))
plt.ylabel('Log Global Residual')
plt.xlabel('Time Step')
plt.grid()
plt.savefig('Re'+str(Re)+'dt'+str(dt)+'Residual'+str(mesh.numOfxNodes)+ 'x'+str(mesh.numOfyNodes)+'.png')
plt.close()
print "plotting solution"
mesh.plot(solution="Phi",plotType='2d', dt=dt, Re=Re)
mesh.plot(solution="w",plotType='2d', dt=dt, Re=Re)




