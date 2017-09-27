from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
from scipy.optimize import curve_fit
sys.path.append("/home/rajsekhar/PLUTO/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()
'''
#average is a quantity that denotes the number of k^2 data points we have averaged over to get E(k) and k at the desired point

average=20
'''
#bin_size denotes the division of k space into bins of size 2*pi*bin_size
bin_size=4

#n is an array that stores the size of the simulation domain
n=np.array([128,128,128])

#Declare all parameters and filenames, file location

filedir="/home/rajsekhar/MHD-TURBULE-01/HD-setup/Data/"

z=1.0
solver="hllc"
filedir+="Z"+str(z)+'/'+str(n[0])+'/'+solver+'/'
for filenumber in xrange(13,14):

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

    #Square and add the absolute values of Vk1, Vk2, Vk3

    Vk_sq=np.add(np.add(np.square(np.absolute(Vk1)),np.square(np.absolute(Vk2))),np.square(np.absolute(Vk3)))

    #Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(nz/2-k)^2
    k_sq=np.fromfunction(lambda i,j,k:((n[0]/2)-i)**2+((n[1]/2)-j)**2+((n[2]/2)-k)**2,n)

    #Energy density is given by k^2 * Vk^2
    E_k=np.multiply(Vk_sq,k_sq)


    #Flatten E_k and k now, and store them in a single  2D array

    K=np.transpose(np.vstack((np.ndarray.flatten(k_sq),np.ndarray.flatten(E_k))))

    #K[][0] stores k, K[][1] stores E(k). 

    print("--- %s seconds ---" % (time.time() - start_time))

    #The following code adds up E(k) values for the same k^2

    c=Counter()
    for k, v in K:
    #c[k] stores all E(k) values in corresponding k^2 index of the array
        c[k] += v

    #K_rad takes non-empty elements of c, and stores the corresponding k(taking square root of k_sq) and E(k) values
    K_rad=[[2*np.pi*np.sqrt(var),c[var]] for var in c if var]
    K_rad=np.array(K_rad)
    #the previous procedure leaves out the value corresponding to k=0, so we insert it
    K_rad=np.insert(K_rad,0,c[0],axis=0)
    K_E_comp=np.multiply(K_rad[:,1],np.power(K_rad[:,0],5/3))
    
    #Take bins of a certain size, and add up values corresponding to thse bins

    index=0
    K_avg=np.empty([0,3],dtype=float)
    for i in xrange(0,int((1/bin_size)*0.5*np.sqrt(np.sum(np.square(n))))):
        k_min=2*bin_size*np.pi*i
        k_max=2*bin_size*np.pi*(i+1)
        E=0
        E_comp=0
        for j in xrange(index,np.size(K_rad[:,0])):
            if K_rad[:,0][j]<k_max:
                E+=K_rad[:,1][j]
                E_comp+=K_E_comp[j]
            else:
                K_avg=np.append(K_avg,[[k_min,E,E_comp]],axis=0)
                index=j+1
                break
        
    """
    #a denotes the size of the truncated array so that k and E(k)  can be converted to 2-D arrays for averaging

    a=np.shape(K_rad)[0]-np.shape(K_rad)[0]%average

    #Truncate the array so that it can be reshaped

    K_rad=K_rad[0:a]
    K_rad_k=K_rad[:,0].reshape((-1,average))
    K_rad_E=K_rad[:,1].reshape((-1,average))
    
    #K_E_comp stores the values of energy compensated with the kolmogorov scaling, i.e. it stores E(k)*k^(5/3)
    
    K_E_comp=np.multiply(K_rad_E,np.power(K_rad_k,5/3))

    #divide into bins of size average, sum up values inside each bin

    K_rad_k=np.mean(K_rad_k,1)
    K_rad_E=np.sum(K_rad_E,1)
    K_E_comp=np.sum(K_E_comp,1) 
    
    K=np.transpose(np.vstack((K_rad_k,K_rad_E,K_E_comp)))
    print("--- %s seconds ---" % (time.time() - start_time))
    """    
    np.savetxt(filedir+'power_spectrum_'+str(filenumber)+'.txt',K_avg)
    #Curve fitting
    """
    def func(x, a1,a2):
        y1 = c1*((x-x1)**(a1) )[x<=xc]
        y2 = c2*((x-x2)**(a2) )[x>xc and x<=xd]
        y3 = c3*((x-x3)**(a3) )[x>xd]
        y = np.concatenate((y1,y2,y3))
        return y
    popt,pconv=curve_fit(func,K_rad_k,K_rad_E,bounds=(0, [2.,-2.,-3.])
    """
    #Plot the data 

    plt.figure()
    #Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
    plt.loglog(K_avg[:,0],K_avg[:,2])
    plt.xlabel('k')
    plt.ylabel('E(k)*$k^{5/3}$')
    #plt.ylim(10**11,10**15)
    plt.title('Compensated E(k) vs k for t='+str(filenumber/10))

    plt.savefig(filedir+'log_E_k_compensated'+str(filenumber)+'.png')
    #This is to plot the original power spectrum, without any compensation
    plt.figure()
    plt.loglog(K_avg[:,0],K_avg[:,1])
    plt.xlabel('k')
    plt.ylabel('E(k)')
    #plt.ylim(10**11,10**13)

    plt.title('E(k) vs k for t='+str(filenumber/10))
    plt.savefig(filedir+'log_E_k'+str(filenumber)+'.png')

    print("--- %s seconds ---" % (time.time() - start_time))
