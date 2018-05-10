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
mat1 = numpy.fromfile('./wjacobi_data_neigh.dat', dtype = numpy.float64)
mat1 = mat1.reshape((-1, res*3+1))
print mat1.shape
mat1 = mat1[0:-1]

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
	# print index,numpy.sum(numpy.square(residual)),omega,error,residual,ck

error_line1/=error_line1[0]
residual_line1/=residual_line1[0]
axErr.semilogy(range(0,len(error_line1)), error_line1)

axR.semilogy(range(0,len(residual_line1)), residual_line1)


omega_mat = [2,res]
omega_mat[0] = omega_line1
omega_mat[1] = range(0,len(omega_line1))
# print omega_mat
omega_mat = numpy.array(omega_mat)
omega_mat=omega_mat[:,omega_mat[0].argsort()]
# print omega_mat

omegaidxs=[]
omegax=[]
omegay=[]
oldomega = 0
for idx in range(0,len(omega_line1)):
# for omega in omega_line1:
	omega = omega_mat[0,idx]
	
	if abs(omega-oldomega<0.0000001):
		# omega_numline[-1]+=1
		omegaidxs.append(omega_mat[1,idx])
	else:
		print  oldomega, omegaidxs
		if oldomega!=0:
			omegax.append(1./oldomega)
			omegay.append(len(omegaidxs))
		omegaidxs=[omega_mat[1,idx]]
		oldomega = omega
print  1/oldomega, omegaidxs
omegax.append(1./oldomega)
omegay.append(len(omegaidxs))
axOmega.semilogy(range(0,len(omega_line1)), omega_line1)
# axOmega.plot(omegax, omegay)


axErr.legend()
axR.legend()
axOmega.legend()
# plt.savefig('wavebrn.png',dpi=100)
plt.show()
plt.close()

