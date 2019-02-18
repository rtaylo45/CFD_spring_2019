import MeshType as Me
import PhysicsType as Phy
import matplotlib.pyplot as plt
from math import log10

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=1.0, yLength=1.0, xNodes=101, yNodes=101)
dt = 0.001
Re = 100.0

print "Setting up problem"
# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, dt=dt,Re=Re)

print "Setting BC"
problem.NavierObj.setBC(side="north", BC=1., BCType=0)
problem.NavierObj.setBC(side="south", BC=0., BCType=0)
problem.NavierObj.setBC(side="east", BC=0., BCType=0)
problem.NavierObj.setBC(side="west", BC=0., BCType=0)

print "Solving problem"
# solve the system
timeSteps, lapDiffsTemp, navDiffs = problem.solve(solveType=1)
lapDiffs = []
for val in lapDiffsTemp:
	try:
		lapDiffs.append(log10(val))
	except:
		lapDiffs.append(0.0)
navDiffs = [log10(val) for val in navDiffs]

#print mesh.globalError
plt.plot(timeSteps,lapDiffs,'o')
plt.title('Stream Function Residual '+ str(mesh.numOfxNodes)+ ' x '+str(mesh.numOfyNodes))
plt.ylabel('Log Global Residual')
plt.xlabel('Time Step')
plt.grid()
plt.savefig('StreamRe'+str(Re)+'dt'+str(dt)+'Residual'+str(mesh.numOfxNodes)+ 'x'+str(mesh.numOfyNodes)+'.png')
plt.close()

plt.plot(timeSteps,navDiffs, 'o')
plt.title('Vorticity Residual '+ str(mesh.numOfxNodes)+ ' x '+str(mesh.numOfyNodes))
plt.ylabel('Log Global Residual')
plt.xlabel('Time Step')
plt.grid()
plt.savefig('VorticityRe'+str(Re)+'dt'+str(dt)+'Residual'+str(mesh.numOfxNodes)+ 'x'+str(mesh.numOfyNodes)+'.png')
plt.close()
print "plotting solution"
mesh.plot(solution="Phi",plotType='2d', dt=dt, Re=Re)
mesh.plot(solution="w",plotType='2d', dt=dt, Re=Re)








