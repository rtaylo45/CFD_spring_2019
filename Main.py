import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log10

dt = 1.0
Re = 100.0
print "Time steps [s]", dt
print "Reynolds Number", Re
print 

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(yLength=1.0,xLength=1.0,xiLength=1.0, 
	etaLength=1.0, xiNodes=101, etaNodes=101)

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, dt=dt,Re=Re)
print "Solving problem"
# solve the system
timeSteps, Diffs = problem.solve(solveType=1)
plt.plot(timeSteps,Diffs, 'o')
plt.title('Global Residual '+ str(mesh.numOfxiNodes)+ ' x '+str(mesh.numOfetaNodes))
plt.ylabel('Log Global Residual')
plt.xlabel('Time Step')
plt.grid()
plt.savefig('Re'+str(int(Re))+'Residual'+str(mesh.numOfxiNodes)+ 'x'+str(mesh.numOfetaNodes)+'.png')
plt.close()

print "plotting solution"
mesh.plot(solution="Phi",plotType='2d', dt=dt, Re=Re)
mesh.plot(solution="w",plotType='2d', dt=dt, Re=Re)
mesh.plot(solution="velocity",plotType='vector', dt=dt, Re=Re)
mesh.plotCenterLineUVelocity()
mesh.plotCenterLineVVelocity()



