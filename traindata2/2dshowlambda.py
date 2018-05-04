import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *
from math import pi
from mpl_toolkits.mplot3d.axes3d import Axes3D


def ShowLambdaWkm():
	xres=32
	yres=32
	omega=1

	mat = numpy.zeros([xres,yres])
	for y in range(1,yres+1):
		for x in range(1,xres+1):

			mat[x-1,y-1] = 1.-omega*( pow(sin(pi*x/2/(xres+1)),2) + pow(sin(pi*y/2/(yres+1)),2) )
			print x,y,mat[x-1,y-1]

	xarr = linspace(1,xres,xres)
	yarr = linspace(1,yres,yres)
	print xarr, yarr

	X,Y = meshgrid(xarr ,yarr)

	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.1, hspace=0.1)
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	ax.set_title('result')

	print mat
	p = ax.plot_surface(X, Y, mat, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7)#,norm=norm)
	cb = fig.colorbar(p, shrink=0.5)
	#plt.savefig('lambdawkm.png',dpi=100)
	plt.show()

ShowLambdaWkm()

def ShowMat(mat, xres, yres):
	mat = mat.reshape((xres, yres))
	xarr = linspace(0, xres-1, xres)
	xaxis = linspace(0, xres-1, xres)
	yarr = linspace(0, yres-1, yres)
	yaxis = linspace(0, yres-1, xres)

	X,Y = meshgrid(xarr ,yarr)

	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.1, hspace=0.1)
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	ax.set_title('result')
	ax.set_xticks(xaxis)
	ax.set_yticks(yaxis)

	p = ax.plot_surface(X, Y, mat, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7)#,norm=norm)
	cb = fig.colorbar(p, shrink=0.5)
	plt.show()
	
