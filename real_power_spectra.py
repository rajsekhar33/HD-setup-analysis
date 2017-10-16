from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from collections import Counter
import time
from scipy.optimize import curve_fit
sys.path.append("/home/rajsekhar/PLUTO41_old/Tools/pyPLUTO")
import pyPLUTO as pp

#Compute how long the simulation takes
start_time = time.time()
#bin_size denotes the division of k space into bins of size 2*pi*bin_size
bin_size=1

#n is an array that stores the size of the simulation domain
n=np.array([128,128,128])

#Declare all parameters and filenames, file location

filedir="/home/rajsekhar/PLUTO41_old/3D_turb/Tau_c_2/"
nu=0
if nu!=0:
    filedir=filedir+'with_viscosity/'
filedir=filedir+str(n[0])+'/'
if nu!=0:
    filedir=filedir+str(nu)+'/'

for filenumber in xrange(3,4):

    #Load data files using pp.pload

    D=pp.pload(filenumber,w_dir=filedir)

    #Perform fourier transform from real to complex

    Vk1=(float(1)/np.prod(n))*np.fft.rfftn(D.vx1,s=n)
    Vk2=(float(1)/np.prod(n))*np.fft.rfftn(D.vx2,s=n)
    Vk3=(float(1)/np.prod(n))*np.fft.rfftn(D.vx3,s=n)

    Vk1=np.fft.fftshift(Vk1,axes=(0,1))
    Vk2=np.fft.fftshift(Vk2,axes=(0,1))
    Vk3=np.fft.fftshift(Vk3,axes=(0,1))

    print("--- %s seconds ---" % (time.time() - start_time))

    #Square and add the absolute values of Vk1, Vk2, Vk3

    Vk_sq=np.add(np.add(np.square(np.absolute(Vk1)),np.square(np.absolute(Vk2))),np.square(np.absolute(Vk3)))

    #Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(k)^2
    k_sq=np.fromfunction(lambda i,j,k:(i-n[0]/2)**2+(j-n[1]/2)**2+k**2,np.shape(Vk1))

    #Energy density is given by Vk^2
    E_k=Vk_sq

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
    K_rad=[[2*np.pi*np.sqrt(var),c[var]/(float(2)*np.pi*bin_size)] for var in c if var]
    K_rad=np.array(K_rad)
    
    #Take bins of a certain size, and add up values corresponding to these bins
    #Some bins are scaled in stretched manner, the rest are in logarithmic, once k values are large enough
    index=0
    #no_bins sets the total no. of bins
    no_bins=100
    #ratio sets the no. of points to be scaled using stretched bin method and the no. of points to be stretched using logarithmic binning method
    ratio=10
    max_k=np.ndarray.max(K_rad[:,0])
    #r_bin1 sets the ratio between cnsecutive bin sizes for stretched binning method
    r_bin1=np.power(max_k/ratio,float(1)/(9*no_bins/10))
    #r_bin2 sets the ratio between consecutive k values for logarithmic binning method
    r_bin2=ratio**(10/no_bins)
    K_avg=np.empty([0,2],dtype=float)
    k_max=2*bin_size*np.pi-bin_size*np.pi
    for i in xrange(1,no_bins):
        k_min=k_max
        #for stretched binning method, at lower k values
        if k_min<=max_k/ratio:
            k_max=k_max+2*bin_size*np.pi*r_bin1**i
        #for logarithmic binning method, at higher k values
        else:
            k_max=k_max*r_bin2
        #leave the loop if k_max exceeds maximum possible k value
        if k_max>np.ndarray.max(K_rad[:,0]):
            break
        delta_k=k_max-k_min
        E=0
        for j in xrange(index,np.size(K_rad[:,0])):
            if K_rad[:,0][j]<k_max:
                #add up values in each bin
                E+=K_rad[:,1][j]
            else:
                #add the entry to the array
                K_avg=np.append(K_avg,[[k_min,E/(delta_k)]],axis=0)
                if(E!=0):
		    index=j
                break
    #calculate the compensated power spectrum by multiplying with k^(5/3)
    K_E_comp=np.multiply(K_avg[:,1],np.power(K_avg[:,0],5/3))
    #add this as the third column to K_avg array
    K_avg=np.transpose(np.vstack((K_avg[:,0],K_avg[:,1],K_E_comp)))
    np.savetxt(filedir+'power_spectrum_'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.txt',K_avg)
    
    np.savetxt('./'+str(n[0])+'/'+str(int(bin_size*100))+'/power_spectrum_'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.txt',K_avg)
    
#Curve fitting
    #Plot the data 

    plt.figure()
    #Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
    plt.plot(K_avg[:,0],K_avg[:,2],'o-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('k')
    plt.ylabel('E(k)*$k^{5/3}$')
    #plt.ylim(10**11,10**15)
    plt.title('Compensated E(k) vs k for t='+str(filenumber)+' '+'bin size = '+str(bin_size)+', size =' +str(n[0])+'*'+str(n[1])+'*'+str(n[2]))

    plt.savefig(filedir+'log_E_k_compensated'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.png')
    plt.savefig('./'+str(n[0])+'/'+str(int(bin_size*100))+'/log_E_k_compensated'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.png')
    
    #This is to plot the original power spectrum, without any compensation
    plt.figure()
    plt.plot(K_avg[:,0],K_avg[:,1],'o-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('k')
    plt.ylabel('E(k)')
    #plt.ylim(10**11,10**13)
    plt.title('E(k) vs k for t='+str(filenumber)+' '+'bin size = '+str(bin_size)+', size =' +str(n[0])+'*'+str(n[1])+'*'+str(n[2]))
    plt.savefig(filedir+'log_E_k'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.png')
    plt.savefig('./'+str(n[0])+'/'+str(int(bin_size*100))+'/log_E_k'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'_nu_'+str(nu)+'.png')

    print("--- %s seconds ---" % (time.time() - start_time))
