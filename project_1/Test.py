import MeshType as Me
import PhysicsType as Phy

mesh = Me.Mesh(xLength=15, yLength=20, xNodes=4, yNodes=5)


mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=0.)
mesh.setBC(side="west", BC=0.)
mesh.setBC(side="north", BC=100.)

mesh.finalize()

problem = Phy.Physics(mesh=mesh, problemType="Laplace")
problem.solve()

"""
for i in xrange(mesh.numOfxNodes):
	for j in xrange(mesh.numOfyNodes):
		node = mesh.getNodeByLoc(i,j)
		print i,j,node.error

print "Iterations ", problem.iterations
"""

mesh.plot2D()

