import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *

fig = plt.figure(figsize=(16,6))
axOmega = fig.add_subplot(1, 3, 1)
plt.legend(loc='upper right')
plt.title('Omega')
plt.ylabel('omega')
plt.xlabel('iteration')

axErr = fig.add_subplot(1, 3, 2)
plt.legend(loc='upper right')
plt.title('Error')
plt.ylabel('|error|^2')
plt.xlabel('iteration')

axR = fig.add_subplot(1, 3, 3)
plt.legend(loc='upper right')
plt.title('Residual')
plt.ylabel('|r|^2')
plt.xlabel('iteration')

plt.subplots_adjust(left=0.05, right=0.95)

res = 32*32
mat1 = numpy.fromfile('./wjacobi_data_rck.dat', dtype = numpy.float64)
print mat1.shape
print 25784320/32/32
print 25180*32*32
mat1 = mat1.reshape((-1, res*2))
print mat1.shape

error_line1=[]
residual_line1=[]
omega_line1=[]
for index,line in enumerate(mat1):
	# error    = line[0    :res]#error
	# residual = line[res  :res*2]#residual
	# ck       = line[res*2:res*3]#ck
	# omega    = line[res*3]#omega
	residual = line[0  :res]#residual
	ck       = line[res:res*2]#ck
	

	#find maxck
	maxck = argmax(abs(ck))
	#divide into k for x-axis, m for y-axis
	#need +1 op, because we use k+1 when calculating omega
	maxk = maxck%32 +1
	maxm = maxck/32 +1
	omega_km = 1./( sin(maxk*numpy.pi/2/(32+1))*sin(maxk*numpy.pi/2/(32+1)) + sin(maxm*numpy.pi/2/(32+1))*sin(maxm*numpy.pi/2/(32+1)));

	# error_line1.append(numpy.sum(numpy.square(error)))
	residual_line1.append(numpy.sum(numpy.square(residual)))
	# omega_line1.append(omega)
	# print index,numpy.sum(numpy.square(residual)),omega,error,residual,ck


# axErr.semilogy(range(0,len(error_line1)), error_line1)

axR.semilogy(range(0,len(residual_line1)), residual_line1)

axOmega.semilogy(range(0,len(omega_line1)), omega_line1)

axErr.legend()
axR.legend()
axOmega.legend()
# plt.savefig('wave.png',dpi=100)
plt.show()
plt.close()

