import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *
from math import pi
from mpl_toolkits.mplot3d.axes3d import Axes3D

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
	
def RenderMat(prefix, ind, xres, yres):
	indstr = str(ind)
	if len(indstr) == 1:
		indstr='000'+indstr
	elif len(indstr) == 2:
		indstr='00'+indstr
	elif len(indstr) == 3:
		indstr='0' + indstr
	filename = prefix+indstr+'.dat'
	print 'render iter '+indstr
	mat = numpy.fromfile(filename, dtype = numpy.float64)
	mat = mat.reshape((xres, yres))
	xarr = linspace(0, xres-1, xres)
	xaxis = linspace(0, xres, xres/4+1)
	yarr = linspace(0, yres-1, yres)
	yaxis = linspace(0, yres, yres/4+1)
	norm = colors.Normalize(vmin=0, vmax=1)

	X,Y = meshgrid(xarr ,yarr)

	fig = plt.figure(figsize=(8,6))
	plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.1, hspace=0.1)
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	ax.set_title('Iteration '+indstr)
	ax.set_xticks(xaxis)
	ax.set_yticks(yaxis)

	p = ax.plot_surface(X, Y, mat, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7,norm=norm)
	ax.set_zlim3d(-0.1,1.1)
	cb = fig.colorbar(p, shrink=0.5,norm=norm)
	plt.savefig(prefix+indstr+'.png',dpi=100)
	plt.close()

def RenderSerial(prefix, xres, yres):
	for ind in range(0,100):
		RenderMat(prefix,ind, xres, yres)

def RenderIter(filename, xres, yres):
	mat = numpy.fromfile(filename, dtype = numpy.float64)
	mat = mat.reshape((-1, xres*yres*3+1))
	# mat = mat[0:6]
	for index,line in enumerate(mat):
		error    = line[0          :xres*yres]#error
		residual = line[xres*yres  :xres*yres*2]#residual
		ck       = line[xres*yres*2:xres*yres*3]#ck
		omega    = line[xres*yres*3]#omega

		errorsum = numpy.sum(numpy.square(error))
		residualsum = numpy.sum(numpy.square(residual))

		indstr = str(index)
		if len(indstr) == 1:
			indstr='000'+indstr
		elif len(indstr) == 2:
			indstr='00'+indstr
		elif len(indstr) == 3:
			indstr='0' + indstr

		fig = plt.figure(figsize=(14,6))

		# axDown = fig.add_subplot(2,1,2)
		axDown = plt.axes([0., 0., 1.0, 0.1])
		axDown.text(0.5, 0, 'iter = '+indstr+
							' omg = '+str(omega)+
							'\ne^2  = '+str(errorsum)+
							' r^2 = '+str(residualsum),
							horizontalalignment='center')

		
		# plt.text(0.5, 0.1, 'iter = '+indstr+,
			# fontdict={'size': 16, 'color': 'r'})
		# axCk = fig.add_subplot(2, 2, 1, projection='3d')
		axCk = plt.axes([0., 0.1, 0.5, 0.9], projection='3d')
		# plt.legend(loc='best')
		plt.title('Spectrum')
		plt.ylabel('ck')
		plt.xlabel('k')

		# axWave = fig.add_subplot(2, 2, 2, projection='3d')
		axWave = plt.axes([0.5, 0.1, 0.5, 0.9], projection='3d')
		plt.legend(loc='best')
		plt.title('Wave')
		plt.ylabel('y')
		plt.xlabel('x')


		xarr = linspace(1,xres,xres)
		yarr = linspace(1,yres,yres)

		X,Y = meshgrid(xarr ,yarr)

		error = error.reshape((xres, yres))
		ck = ck.reshape((xres, yres))

		p = axCk.plot_surface(X, Y, numpy.abs(ck), rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7)#,norm=norm)
		axCk.set_zlim3d(0,1)
		cb = fig.colorbar(p, shrink=0.5)
		p = axWave.plot_surface(X, Y, error, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False,alpha=0.7)#,norm=norm)
		
		plt.savefig(indstr+'.png',dpi=100)
		plt.close()

RenderIter('wjacobi_data.dat',32,32)

