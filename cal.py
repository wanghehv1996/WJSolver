from math import *

res = 32
def getLambda(omega, k):
	return 1.-2.*omega*sin(k*pi/2/(res+1))*sin(k*pi/2/(res+1))

def getOmega(lamb, k):
	return (1.-lamb)/(2*sin(k*pi/2/(res+1))*sin(k*pi/2/(res+1)))
def getK(lamb, omega):
	sinv = sqrt((1.-lamb)/omega)
	
	if sinv>1:
		return res
	return asin(sinv)*2*(res+1)/pi#boom here

lamb=0
k = 1
while k<32:
	omega = getOmega(lamb,k)
	print omega#last omega is useless
	k = getK(-lamb, omega)

	# k = int(k)+1
	print omega,k
	lamb = 0.1


