from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
import struct as st

#Declare all parameters and filenames, file location
data='vort'
order=2
difftype='c'
filetype='.dbl'
filenumber=20
nx=512
ny=512
nz=512
nv=3

fileno=str(filenumber).rjust(4,'0')
filedir="/mnt/lustre/ug4/ugrajs/fiducial_runs/512/"
filename=filedir+data+'.'+difftype+str(order)+fileno+filetype
print filename
#Initialise all arrays as required into arrays of zeroes of required dimensions

data = np.zeros(nv*nx*ny*nz)
v=np.zeros((nv,nx,ny,nz))
vort=np.zeros((nx,ny))#To be changed if other components of vorticity are to be plotted

x=np.zeros(nx)
y=np.zeros(ny)
x1=np.zeros((nx,ny))
y1=np.zeros((nx,ny))
dx,dy=1/nx,1/ny

# generate 2 2d grids for the x & y bounds
y1, x1 = np.mgrid[slice(-0.5, 0.5 + dy, dy),
                slice(-0.5, 0.5 + dx, dx)]
#Open file to be read in binary format

with open(filename, 'rb') as f:
    data = f.read()

#Loop over the data in the file and store them in a 4D array, after unpacking the data from binary format

for dir in range(0,nv):
    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                l=8*(dir*nx*ny*nz+i*ny*nz+j*nz+k)
                v[dir,i,j,k]=st.unpack('d',data[l:l+8])[0]

#Loop again to make a 2D array, keeping one of x,y or z co-ordinates constant

for i in range(0,nx):
    for j in range(0,ny):
        x[i]=(i-nx/2)/nx
        y[i]=(j-ny/2)/ny
        vort[i,j]=v[0,i,j,64]
#The following formulae should be changed slightly if we are plotting other velocity components

plt.figure()
plt.pcolormesh(x1,y1,vort,cmap=plt.cm.PuBuGn,vmin=-300,vmax=300)
plt.colorbar()
plt.title('Vorticity at different grid points at t='+str(0.2*filenumber))
plt.savefig('vort'+'z'+difftype+str(order)+fileno+'.png',dpi=1000)
plt.show()
