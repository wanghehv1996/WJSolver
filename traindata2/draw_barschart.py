import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *
from mpl_toolkits.mplot3d.axes3d import Axes3D

def barschart2d():
	xres = 32
	yres = 32
	mat = numpy.fromfile('./wjacobi_data.dat', dtype = numpy.float64)
	mat = mat.reshape((-1, xres*yres*3+1))
	bars = numpy.zeros([xres,yres])
	# mat = mat[0:6]
	for index,line in enumerate(mat):
		error    = line[0          :xres*yres]#error
		residual = line[xres*yres  :xres*yres*2]#residual
		ck       = line[xres*yres*2:xres*yres*3]#ck
		omega    = line[xres*yres*3]#omega
		
		ck = ck.reshape((32,32))
		# k = numpy.argmax(numpy.abs(ck))
		k = numpy.unravel_index(numpy.argmax(ck, axis=None), ck.shape)
		bars[k] = bars[k]+1
		print k, omega

	X,Y = meshgrid(range(0,32) ,range(0,32))
	fig = plt.figure(figsize=(8,6))
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	# p = ax.plot_surface(X, Y, bars, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7)#,norm=norm)
	
	xpos = numpy.arange(0,xres,1)    # Set up a mesh of positions
	ypos = numpy.arange(0,yres,1)
	xpos, ypos = numpy.meshgrid(xpos+0.5, ypos+0.5)

	xpos = xpos.flatten()   # Convert positions to 1D array
	ypos = ypos.flatten()
	zpos = numpy.zeros(xres*yres)

	dx = 1 * numpy.ones_like(zpos)
	dy = dx.copy()
	dz = bars.flatten()
	ax.bar3d(xpos,ypos,zpos, dx, dy, dz, shade = True)

	plt.show()

	# 
	# plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.1, hspace=0.1)
	
	# ax.set_title('result')

	# print mat
	
	# cb = fig.colorbar(p, shrink=0.5)
	# #plt.savefig('lambdawkm.png',dpi=100)
	# plt.show()

barschart2d()