import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *

fig = plt.figure(figsize=(16,6))
axOmega = fig.add_subplot(1, 3, 1)
plt.title('Omega')
plt.ylabel('omega')
plt.xlabel('iteration')

axErr = fig.add_subplot(1, 3, 2)
plt.title('Relative Error')
plt.ylabel('error')
plt.xlabel('iteration')

axR = fig.add_subplot(1, 3, 3)
plt.title('Relative Residual')
plt.ylabel('residual')
plt.xlabel('iteration')

plt.subplots_adjust(left=0.05, right=0.95)

res = 32*32
mat1 = numpy.fromfile('./wjacobi_data_1.dat', dtype = numpy.float64)
mat1 = mat1.reshape((-1, res*3+1))

mat2 = numpy.fromfile('./wjacobi_data_2div3.dat', dtype = numpy.float64)
mat2 = mat2.reshape((-1, res*3+1))

mat3 = numpy.fromfile('./wjacobi_data_largest.dat', dtype = numpy.float64)
mat3 = mat3.reshape((-1, res*3+1))

error_line1=[]
residual_line1=[]
omega_line1=[]
for index,line in enumerate(mat1):
	error    = line[0    :res]#error
	residual = line[res  :res*2]#residual
	ck       = line[res*2:res*3]#ck
	omega    = line[res*3]#omega

	error_line1.append(numpy.sum(numpy.square(error)))
	residual_line1.append(numpy.sum(numpy.square(residual)))
	omega_line1.append(omega)
error_line1/=error_line1[0]
residual_line1/=residual_line1[0]

error_line2=[]
residual_line2=[]
omega_line2=[]
for index,line in enumerate(mat2):
	error    = line[0    :res]#error
	residual = line[res  :res*2]#residual
	ck       = line[res*2:res*3]#ck
	omega    = line[res*3]#omega

	error_line2.append(numpy.sum(numpy.square(error)))
	residual_line2.append(numpy.sum(numpy.square(residual)))
	omega_line2.append(omega)
error_line2/=error_line2[0]
residual_line2/=residual_line2[0]

error_line3=[]
residual_line3=[]
omega_line3=[]
for index,line in enumerate(mat3):
	error    = line[0    :res]#error
	residual = line[res  :res*2]#residual
	ck       = line[res*2:res*3]#ck
	omega    = line[res*3]#omega

	error_line3.append(numpy.sum(numpy.square(error)))
	residual_line3.append(numpy.sum(numpy.square(residual)))
	omega_line3.append(omega)
error_line3/=error_line3[0]
residual_line3/=residual_line3[0]


axErr.semilogy(range(0,len(error_line1)), error_line1,label='omega=1')
axErr.semilogy(range(0,len(error_line2)), error_line2,label='omega=0.666')
axErr.semilogy(range(0,len(error_line3)), error_line3,label='our method')

axR.semilogy(range(0,len(residual_line1)), residual_line1,label='omega=1')
axR.semilogy(range(0,len(residual_line2)), residual_line2,label='omega=0.666')
axR.semilogy(range(0,len(residual_line3)), residual_line3,label='our method')

axOmega.semilogy(range(0,len(omega_line1)), omega_line1,label='omega=1')
axOmega.semilogy(range(0,len(omega_line2)), omega_line2,label='omega=0.666')
axOmega.semilogy(range(0,len(omega_line3)), omega_line3,label='our method')

axErr.legend(loc=4)
axR.legend(loc=4)
axOmega.legend(loc=4)
plt.savefig('2dcompareomega.png',dpi=100)
# plt.show()
plt.close()

