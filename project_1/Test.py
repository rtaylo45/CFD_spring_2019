import MeshType as Me
import PhysicsType as Phy

# Builds mesh
mesh = Me.Mesh(xLength=15.0, yLength=20.0, xNodes=31, yNodes=41)

# sets the boundary conditions
mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=0.)
mesh.setBC(side="west", BC=0.)
mesh.setBC(side="north", BC=100.)

# Finalize the mesh
mesh.finalize()

# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

# solve the system
problem.solve(pattern=0)
"""
for i in xrange(mesh.numOfxNodes):
	for j in xrange(mesh.numOfyNodes):
		node = mesh.getNodeByLoc(i,j)
		print i,j,node.error

print "Iterations ", problem.iterations
"""
mesh.plot2D(solution="exact")
mesh.plot2D(solution="approx")

