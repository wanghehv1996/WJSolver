import matplotlib.pyplot as plt
from matplotlib import *
from numpy import *
import math

fig = plt.figure(figsize=(8,6))
plt.legend(loc='upper right')
plt.title('lambda-k')
plt.ylabel('lambda')
plt.xlabel('k')

res = 32

omega1 = []
omega2d3 = []
omega2 = []
for k in range(0,res+1):
	omega1.append(1-2*pow( math.sin(k*math.pi/2/(res+1)) ,2))
	omega2d3.append(1-2*2/3*pow( math.sin(k*math.pi/2/(res+1)) ,2))
	omega2.append(1-2*2*pow( math.sin(k*math.pi/2/(res+1)) ,2))
	# print index,numpy.sum(numpy.square(residual)),omega,error,residual,ck

plt.plot(range(0,res+1), omega1, label = 'omega=1')
plt.plot(range(0,res+1), omega2d3, label = 'omega=2/3')
plt.plot(range(0,res+1), omega2, label = 'omega=2')
plt.axhline(0,color='black')
plt.ylim(-3, 3)
plt.legend()
# plt.savefig('1dlambda_k.png',dpi=100)
plt.show()
plt.close()

