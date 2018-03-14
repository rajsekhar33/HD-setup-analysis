from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
import struct as st

#Declare all parameters and filenames, file location
nx=256
ny=256
for filenumber in range(0,34,3):
	fileno=str(filenumber).rjust(4,'0')
	filedir="/mnt/lustre/ug4/ugrajs/cooling/thermal_heating/256/tabulated_cooling/F5e-1/k0-2/"
	filename=filedir+'sb'+fileno+'.dbl'
	print filename
	#Initialise all arrays as required into arrays of zeroes of required dimensions

	data = np.fromfile(filename, dtype= 'double', count=-1, sep='')
	data=data.reshape(256,256)

	x=np.zeros(nx)
	y=np.zeros(ny)
	x1=np.zeros((nx,ny))
	y1=np.zeros((nx,ny))
	dx,dy=1/nx,1/ny

	# generate 2 2d grids for the x & y bounds
	y1, x1 = np.mgrid[slice(-0.5, 0.5 + dy, dy),
			slice(-0.5, 0.5 + dx, dx)]


	#The following formulae should be changed slightly if we are plotting other velocity components

	plt.figure()
	plt.pcolormesh(x1,y1,np.log(data),cmap=plt.cm.PuBuGn,vmin=-60,vmax=-52)
	plt.colorbar()
	plt.title('Surface brightness at t='+str(0.2*filenumber))
	plt.savefig(filedir+'sbz'+fileno+'.png',dpi=250)
	plt.close()
