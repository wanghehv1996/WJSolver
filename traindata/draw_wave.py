import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *

fig = plt.figure(figsize=(14,6))
axCk = fig.add_subplot(1, 2, 1)
# plt.legend(loc='best')
plt.title('Spectrum')
plt.ylabel('ck')
plt.xlabel('k')

axWave = fig.add_subplot(1, 2, 2)
plt.legend(loc='best')
plt.title('Wave')
plt.ylabel('y')
plt.xlabel('x')

res = 32
mat = numpy.fromfile('./wjacobi_data.dat', dtype = numpy.float64)
mat = mat.reshape((-1, res*3+1))
# mat = mat[0:5]
# indoff = 27;
indoff = 0
# mat = mat[27:32]
for index,line in enumerate(mat):
	# if (index % 4!=0):
	# 	continue
	error    = line[0   :res]#error
	residual = line[res  :res*2]#residual
	ck       = line[res*2:res*3]#ck
	omega    = line[res*3]#omega
	# plt.plot(range(0,res),ck,label='iter='+str(index)+' o'+str(omega))

	xlist=range(1,len(ck)+1,1)
	axCk.plot(xlist, ck,label='iter='+str(index+indoff))
	
	xlist=range(0,len(error),1)
	axWave.plot(xlist, error,label='iter='+str(index+indoff))

axWave.legend()
axWave.grid()
axCk.grid()
# axCk.legend()
plt.savefig('wave.png',dpi=100)
# plt.show()
plt.close()

