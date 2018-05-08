import math
import numpy

xres = 32
yres = 32
for maxk in range(1,33):
	for maxm in range(1,maxk):
		sink2 = math.sin(maxk*math.pi/2/(xres+1)) * math.sin(maxk*math.pi/2/(xres+1))
		sinm2 = math.sin(maxm*math.pi/2/(yres+1)) * math.sin(maxm*math.pi/2/(yres+1))
		print sink2+sinm2, maxk, maxm