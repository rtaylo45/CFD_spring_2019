import MeshType as Me
import PhysicsType as Phy

# Builds mesh
mesh = Me.Mesh(xLength=15.0, yLength=20.0, xNodes=5, yNodes=10)


# sets the boundary conditions
mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=0.)
mesh.setBC(side="west", BC=0.0)
mesh.setBC(side="north", BC=100.)

# Finalize the mesh
mesh.finalize()

# Generates the problem around the mesh
problem = Phy.Physics(mesh=mesh, problemType="Laplace")

# solve the system
problem.solve(pattern=1)


"""
for i in xrange(mesh.numOfxNodes):
	for j in xrange(mesh.numOfyNodes):
		node = mesh.getNodeByLoc(i,j)
		print "node",i,j
		print 
		try:
			print "north",node.north.i,node.north.j
		except:
			print "None"
		try:
			print "south", node.south.i, node.south.j
		except:
			print "None"
		try:
			print "east",node.east.i, node.east.j
		except:
			print "None"
		try:
			print "west",node.west.i, node.west.j
		except:
			print "None"
		print
		print 
"""
"""
for i in xrange(mesh.numOfxNodes):
	for j in xrange(mesh.numOfyNodes):
		node = mesh.getNodeByLoc(i,j)
		if node.solved:
			print "node",i,j, node.solution
"""
#print "Iterations ", problem.iterations

#mesh.plot2D(solution="exact")
mesh.plot2D(solution="approx")








