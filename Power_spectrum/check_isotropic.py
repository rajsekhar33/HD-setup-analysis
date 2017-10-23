from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
sys.path.append("/home/rajsekhar/PLUTO/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()

#average is a quantity that denotes the number of k^2 data points we have averaged over to get E(k) and k at the desired point

average=20

#Declare all parameters and filenames, file location

filedir="/home/rajsekhar/MHD-TURBULE-01/HD-setup/Data/"

#n is an array that stores the size of the simulation domain

n=np.array([256,256,256])
z=1.0
solver="hllc"
filedir+="Z"+str(z)+'/'+str(n[0])+'/'+solver+'/'

filenumber=13

#Load data files using pp.pload

D=pp.pload(filenumber,w_dir=filedir)

#Perform fourier transform from complex to complex(I didn't do real to complex, because the later steps get less complicated if I do it this way

Vk1=np.fft.fftn(D.vx1,s=n)
Vk2=np.fft.fftn(D.vx2,s=n)
Vk3=np.fft.fftn(D.vx3,s=n)

Vk1=np.fft.fftshift(Vk1)
Vk2=np.fft.fftshift(Vk2)
Vk3=np.fft.fftshift(Vk3)

print("--- %s seconds ---" % (time.time() - start_time))

#Find the absolute values of Vk1, Vk2, Vk3

Vk1=np.square(np.absolute(Vk1))
Vk2=np.square(np.absolute(Vk2))
Vk3=np.square(np.absolute(Vk3))

#Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(nz/2-k)^2

k_sq=np.fromfunction(lambda i,j,k:((n[0]/2)-i)**2+((n[1]/2)-j)**2+((n[2]/2)-k)**2,n)

#Flatten the V_ks and k now, and store them in a single  2D array

K=np.transpose(np.vstack((np.ndarray.flatten(k_sq),np.ndarray.flatten(Vk1),np.ndarray.flatten(Vk2),np.ndarray.flatten(Vk3))))
#K[][0] stores k, K[][1] stores Vk1, K[][2] stores Vk2, K[][3] stores Vk3 
#V_sq stores value of velocity squared at each point
V_sq = np.add(np.add(np.square(np.absolute(Vk1)),np.square(np.absolute(Vk2))),np.square(np.absolute(Vk3)))
#E stores the values of energy at all grid points
E=np.multiply(k_sq,V_sq)
#const_x is the x-axis slice that we take in the k-space for our 2d plot
const_x=n[0]/8

#Take mean of a few kxs
E=np.mean(E[int(const_x):2*int(const_x)],axis=0)
#Also make buckets along ky and kz
bucket=2
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)
E=rebin(E,(int(n[1]/bucket),int(n[2]/bucket)))
#Create the grid by slicing required number of times
dy,dz=float(bucket)/n[1],float(bucket)/n[2]
z, y = np.mgrid[slice(-0.5, 0.5 , dz),slice(-0.5, 0.5 , dy)]

plt.figure()

plt.pcolormesh(y,z,E,norm=colors.LogNorm(vmin=None, vmax=None, clip=False),cmap='RdBu_r')
plt.colorbar()
plt.title('E(ky,kz) for $k_x$ averaged from $2\pi(-48)$  to $2\pi(-32)$')
plt.xlabel('$k_y$')
plt.ylabel('$k_z$')
plt.savefig(filedir+'E_isotropic_'+str(filenumber)+'.png')


print("--- %s seconds ---" % (time.time() - start_time))
#The following code adds up E(k) values for the same k^2

c1=Counter()
c2=Counter()
c3=Counter()
for k, v1, v2, v3 in K:
#c[k] stores all E(k) values in corresponding k^2 index of the array
    c1[k] += v1
    c2[k] += v2
    c3[k] += v3


#K_rad takes non-empty elements of c, and stores the corresponding k and E(k) values

VK_rad=[[np.sqrt(var),c1[var],c2[var],c3[var]] for var in c1 if var]
VK_rad=np.array(VK_rad)
#the previous procedure leaves out the value corresponding to k=0, so we insert it
#VK_rad=np.insert((VK_rad,0,c1[0],c2[0],c3[0]),axis=0)    
#a denotes the size of the truncated array so that k and E(k)  can be converted to 2-D arrays for averaging

a=np.shape(VK_rad)[0]-np.shape(VK_rad)[0]%average

#Truncate the array so that it can be reshaped

VK_rad=VK_rad[0:a]
VK_rad_k=VK_rad[:,0].reshape((-1,average))
VK_rad_V1=VK_rad[:,1].reshape((-1,average))
VK_rad_V2=VK_rad[:,2].reshape((-1,average))
VK_rad_V3=VK_rad[:,3].reshape((-1,average))

#Divide into bins of size average, sum up values inside each bin

VK_rad_k=np.mean(VK_rad_k,1)
VK1=np.sum(VK_rad_V1,1)
VK2=np.sum(VK_rad_V2,1)
VK3=np.sum(VK_rad_V3,1)


print("--- %s seconds ---" % (time.time() - start_time))

#Plot the data 

fig, ax = plt.subplots()
ax.loglog(VK_rad_k,VK1,label='$V_x^2$')
ax.loglog(VK_rad_k,VK2,label='$V_y^2$')
ax.loglog(VK_rad_k,VK3,label='$V_z^2$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.7, 1.0))
plt.xlabel('ln(k)')
plt.ylabel('ln(V(k)$^2$)')
plt.title('V(k)$^2$ along different components vs ln(k) for t='+str(filenumber/10))
plt.savefig(filedir+'log_Vk'+str(filenumber)+'.png')
print("--- %s seconds ---" % (time.time() - start_time))
