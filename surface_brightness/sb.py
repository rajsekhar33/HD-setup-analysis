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
diff=(rmax-rmin)*(1./no_bins)
nx,ny,nz=256,256,256

filedir="/mnt/lustre/ug4/ugrajs/cooling/turb_perturb/DkHC/"
for filenumber in range(2, 40, 2):


	#Sort the data arrays into new 1D arrays, in increasing order of distance from the centre
	sx=(np.array(xrange(nx))-nx/2.)**2
	sy=(np.array(xrange(ny))-ny/2.)**2
	
	#This stores distance of each array point from the centre
	dist_array=np.sqrt((np.outer(sx,np.ones(ny))+np.outer(np.ones(nx),sy)).reshape(1,nx*ny))
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
	dist[:,0]=(rmin+diff*dist[:,0])*UNIT_LENGTH

	r_max=rmin#Set initial bin size
	index=0
	for i in xrange(no_bins):
		r_min=r_max
		r_max=r_min+diff
		if(r_max>rmax):break
		sbr,vlosr,sigvlosr=0.,0.,0.	
		for j in xrange(index,nx*ny):
			if(store[:,0][j]<r_max):
				sbr+=store[:,1][j]
				vlosr+=store[:,2][j]
				sigvlosr+=store[:,3][j]
			else:
				dist[:,1][i]=sbr/diff/r_min
				dist[:,2][i]=vlosr/diff/r_min
				dist[:,3][i]=sigvlosr/diff/r_min
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
	plt.savefig(filedir+'sbz'+fileno+'.png',dpi=250)
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
	plt.plot(dist[:,0][1:],dist[:,1][1:])
	plt.title('Surface brightness vs r at t='+str(step_size*filenumber))
	plt.xlabel('r')
	plt.ylabel('Surface brightness')
	plt.ylim(0.05,1.5)
	#plt.xscale('log')
	#plt.yscale('log')
	plt.savefig(filedir+'sbzr'+fileno+'.png',dpi=250)
	plt.close()

	fig, ax = plt.subplots()
	ax.plot(dist[:,0][1:],dist[:,3][1:],label='$\sigma(v_{LOS})(r)$')
	#ax.plot(dist[:,0],dist[:,0]*1e-16,label='$r$')
	ax.set_title('Emission weighted LOS velocity fluctuations vs r at t='+str(step_size*filenumber))
	ax.set_xlabel('r')
	ax.set_ylabel('$\sigma(v_{LOS})(r)$')
	ax.set_ylim(5e6,3e8)
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	ax.legend()
	plt.savefig(filedir+'sigvloszr'+fileno+'.png',dpi=250)
	plt.close()
'''
plt.figure()
plt.plot(dist[:,0],dist[:,2])
plt.title('Emission weighted LOS velocity vs r at t='+str(step_size*filenumber))
plt.xlabel('r')
plt.ylabel('$v_{LOS}(r)$')
plt.xscale('log')
plt.yscale('log')
plt.savefig(filedir+'vloszr'+fileno+'.png',dpi=250)
plt.close()
'''

