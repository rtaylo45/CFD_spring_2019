import MeshType as Me
import PhysicsType as Phy

# Builds mesh
print "Building Mesh"
mesh = Me.Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)

print "Setting BC"
# sets the boundary conditions
mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=25.0)
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
problem.solve(pattern=1)

print "plotting solution"
#mesh.plot2D(solution="exact")
mesh.plot(solution="approx",plotType='3d')








