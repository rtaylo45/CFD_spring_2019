import MeshType as Me
import PhysicsType as Phy

mesh = Me.Mesh(xLength=15, yLength=20, xNodes=4, yNodes=5)


mesh.setBC(side="south", BC=0.)
mesh.setBC(side="east", BC=0.)
mesh.setBC(side="west", BC=0.)
mesh.setBC(side="north", BC=100.)

mesh.finalize()

problem = Phy.Physics(mesh=mesh, problemType="Laplace")
"""
for i in xrange(mesh.numOfxNodes):
	for j in xrange(mesh.numOfyNodes):
		node = mesh.getNodeByLoc(i,j)
		print 
		print "Self ",i,j
		try:
			print "north",node.north.i,node.north.j
		except:
			print
			pass
		try:
			print "south", node.south.i,node.south.j
		except:
			print 
			pass
		try:
			print "east", node.east.i,node.east.j
		except:
			print 
			pass
		try:
			print "west", node.west.i,node.west.j
		except:
			print
			pass
"""