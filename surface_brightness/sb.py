from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
import struct as st
from matplotlib.colors import LogNorm
import pyPLUTO as pp

#Declare all parameters and filenames, file location
nx=256
ny=256
step_size=0.2
CONST_mp=1.67262171e-24
CONST_kB=1.3806505e-16
CONST_mu=0.5
CONST_pc=3.0856775807e18
UNIT_VELOCITY= (1.e8)
UNIT_DENSITY=  (CONST_mp*.1)
UNIT_LENGTH =  (CONST_pc*40.e3)
gamma  = 5./3.
def lam(T):
	return 2e-27*np.sqrt(T)

midx, midy=nx/2, ny/2
rmin,rmax=2.,nx/2.
no_bins=50
ratio=(rmax/rmin)**(1./no_bins)
nx,ny,nz=256,256,256

filedir="/mnt/lustre/ug4/ugrajs/cooling/thermal_heating/256/tabulated_cooling/F1e-3/k12/"
for filenumber in range(3,157,3):


	#Sort the data arrays into new 1D arrays, in increasing order of distance from the centre
	sx=(np.array(xrange(nx))-nx/2.)**2
	sy=(np.array(xrange(ny))-ny/2.)**2
	
	#This stores distance of each array point from the centre
	dist_array=(np.outer(sx,np.ones(ny))+np.outer(np.ones(nx),sy)).reshape(1,nx*ny)
	fileno=str(filenumber).rjust(4,'0')

	filename1=filedir+'sb'+fileno+'.dbl'
	filename2=filedir+'vlos'+fileno+'.dbl'
	filename3=filedir+'sigvlos'+fileno+'.dbl'
	print filename1
	#Initialise all arrays as required into arrays of zeroes of required dimensions
	store=np.zeros((nx*ny,4))
	data1 = np.fromfile(filename1, dtype= 'double', count=-1, sep='')
	store[:,0]=dist_array
	store[:,1]=data1
	data1=np.transpose(data1.reshape(256,256))

	data2 = np.fromfile(filename2, dtype= 'double', count=-1, sep='')
	store[:,2]=data2
	data2=np.transpose(data2.reshape(256,256))

	data3 = np.fromfile(filename3, dtype= 'double', count=-1, sep='')
	store[:,3]=data3
	data3=np.transpose(data3.reshape(256,256))

	#Sort the array on the basis of distancefrom the centre
	store=store[store[:,0].argsort()]
	
	dist=np.zeros((no_bins,4))
	dist[:,0]=np.array(xrange(0,no_bins))
	dist[:,0]=rmin*ratio**dist[:,0]*UNIT_LENGTH

	r_max=rmin#Set initial bin size
	index=0
	for i in xrange(no_bins):
		r_min=r_max
		r_max=r_min*ratio
		if(r_max>rmax):break
		delta_r=r_max-r_min
		sbr,vlosr,sigvlosr=0.,0.,0.	
		for j in xrange(index,nx*ny):
			if(store[:,0][j]<r_max):
				sbr+=store[:,1][j]
				vlosr+=store[:,2][j]
				sigvlosr+=store[:,3][j]
			else:
				dist[:,1][i]=sbr/delta_r
				dist[:,2][i]=vlosr/delta_r
				dist[:,3][i]=sigvlosr/delta_r
				if(sbr!=0): index=j
				break

	x1=np.zeros((nx,ny))
	y1=np.zeros((nx,ny))
	dx,dy=1/nx,1/ny

	# generate 2 2d grids for the x & y bounds
	y1, x1 = np.mgrid[slice(-0.5, 0.5 + dy, dy),
			slice(-0.5, 0.5 + dx, dx)]

	#Make all relevant plots
	
	plt.figure()
	plt.pcolormesh(x1,y1,data1, norm=LogNorm(vmin=0.001, vmax=1),cmap=plt.cm.PuBuGn)
	plt.colorbar()
	plt.title('Surface brightness at t='+str(step_size*filenumber))
	plt.savefig('sbz'+fileno+'.png',dpi=250)
	plt.close()

	plt.figure()
	plt.pcolormesh(x1,y1,data2,cmap=plt.cm.PuBuGn)
	plt.colorbar()
	plt.title('Emission weighted LOS velocity at t='+str(step_size*filenumber))
	plt.savefig(filedir+'vlosz'+fileno+'.png',dpi=250)
	plt.close()

	plt.figure()
	plt.pcolormesh(x1,y1,data3,norm=LogNorm(vmin=data3.min(), vmax=data3.max()),cmap=plt.cm.PuBuGn)
	plt.colorbar()
	plt.title('Emission weighted LOS velocity fluctuations at t='+str(step_size*filenumber))
	plt.savefig(filedir+'sigvlosz'+fileno+'.png',dpi=250)
	plt.close()

	plt.figure()
	plt.plot(dist[:,0],dist[:,1])
	plt.title('Surface brightness vs r at t='+str(step_size*filenumber))
	plt.xlabel('r')
	plt.ylabel('Surface brightness')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(filedir+'sbzr'+fileno+'.png',dpi=250)
	plt.close()

	plt.figure()
	plt.plot(dist[:,0],dist[:,2])
	plt.title('Emission weighted LOS velocity vs r at t='+str(step_size*filenumber))
	plt.xlabel('r')
	plt.ylabel('$v_{LOS}(r)$')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(filedir+'vloszr'+fileno+'.png',dpi=250)
	plt.close()

	plt.figure()
	plt.plot(dist[:,0],dist[:,3])
	plt.title('Emission weighted LOS velocity fluctuations vs r at t='+str(step_size*filenumber))
	plt.xlabel('r')
	plt.ylabel('$\sigma(v_{LOS})(r)$')
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(filedir+'sigvloszr'+fileno+'.png',dpi=250)
	plt.close()

