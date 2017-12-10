from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os, sys
import struct as st

#Declare all parameters and filenames, file location
data='data'
filetype='.dbl'
filenumber=20
nx=1024
ny=1024
nz=1024

fileno=str(filenumber).rjust(4,'0')
filedir="/mnt/lustre/ug4/ugrajs/fiducial_runs/1024/"
filename=filedir+data+'.'+fileno+filetype
print filename
#Initialise all arrays as required into arrays of zeroes of required dimensions

data = np.zeros(nx*ny*nz)
rho=np.zeros((nx,ny,nz))
rho_c=np.zeros((nx,ny))#To be changed if other components of vorticity are to be plotted

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

for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nz):
      l=8*(i*ny*nz+j*nz+k)
      rho[i,j,k]=st.unpack('d',data[l:l+8])[0]
      rho_c[i,j]+=rho[i,j,k]
rho_m=np.mean(rho_c)
#The following formulae should be changed slightly if we are plotting other velocity components

plt.figure()
plt.pcolormesh(x1,y1,np.log(rho_c/rho_m),cmap=plt.cm.PuBuGn)
plt.colorbar()
plt.title(r'$log(\frac{\rho_c}{\rho_0})$ at t='+str(0.2*filenumber))
plt.savefig('column_density'+fileno+'_'+str(nx)+'.png',dpi=1000)
