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

filedir="/home/rajsekhar/PLUTO41_old/3D_turb/Tau_c_20/128/"

for filenumber in xrange(3,4):

    #Load data files using pp.pload

    D=pp.pload(filenumber,w_dir=filedir)

    #Perform fourier transform from real to complex

    rho_k=(float(1)/np.prod(n))*np.fft.rfftn(D.rho,s=n)

    rho_k=np.fft.fftshift(rho_k,axes=(0,1))

    print("--- %s seconds ---" % (time.time() - start_time))

    #Square and add the absolute values of Vk1, Vk2, Vk3

    rhok_sq=np.square(np.absolute(rho_k))

    #Write k_sq as a functional 3d array, with value at each element given by (nx/2-i)^2+(ny/2-j)^2+(k)^2
    k_sq=np.fromfunction(lambda i,j,k:(i-n[0]/2)**2+(j-n[1]/2)**2+k**2,np.shape(rho_k))

    #Power spectrum for density is given by 2 * mod(rho_k)^2
    E_k=2*rhok_sq

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
    #the previous procedure leaves out the value corresponding to k=0, so we insert it
    K_rad=np.insert(K_rad,0,c[0]/(float(2)*np.pi*bin_size),axis=0)
    
    #Take bins of a certain size, and add up values corresponding to thse bins

    index=0
    K_avg=np.empty([0,2],dtype=float)
    for i in xrange(1,int((1/bin_size)*np.sqrt(np.sum(np.square(n/2))))):
        k_min=2*bin_size*np.pi*1.01*i-bin_size*np.pi
        k_max=2*bin_size*np.pi*1.01*i+bin_size*np.pi
        E=0
        for j in xrange(index,np.size(K_rad[:,0])):
            if K_rad[:,0][j]<k_max:
                E+=K_rad[:,1][j]
            else:
                K_avg=np.append(K_avg,[[0.5*(k_min+k_max),E]],axis=0)
                if(E!=0):
		    index=j
                break
    K_E_comp=np.multiply(K_avg[:,1],np.power(K_avg[:,0],5/3))
    K_avg=np.transpose(np.vstack((K_avg[:,0],K_avg[:,1],K_E_comp)))
    np.savetxt(filedir+'density_power_spectrum_'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.txt',K_avg)
    np.savetxt('./'+str(n[0])+'/density_power_spectrum_'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.txt',K_avg)
    
#Curve fitting
    #Plot the data 

    plt.figure()
    #Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
    plt.plot(K_avg[:,0],K_avg[:,2],'o-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('k')
    plt.ylabel(r'$\rho_k^2*k^{5/3}$')
    plt.title(r'Compensated $\rho_k^2$ vs k for t='+str(filenumber)+' '+'bin size = '+str(bin_size)+', size =' +str(n[0])+'*'+str(n[1])+'*'+str(n[2]))

    plt.savefig(filedir+'log_rho_k_compensated'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.png')
    plt.savefig('./'+str(n[0])+'/log_rho_k_compensated'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.png')
    
    #This is to plot the original power spectrum, without any compensation
    plt.figure()
    plt.plot(K_avg[:,0],K_avg[:,1],'o-')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('k')
    plt.ylabel(r'$\rho_k^2$')
    plt.title(r'$\rho_k^2$ vs k for t='+str(filenumber)+' '+'bin size = '+str(bin_size)+', size =' +str(n[0])+'*'+str(n[1])+'*'+str(n[2]))
    plt.savefig(filedir+'log_rho_k'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.png')
    plt.savefig('./'+str(n[0])+'/log_rho_k'+str(filenumber)+'_'+'bin_'+str(int(bin_size*100))+'_'+str(n[0])+'.png')

    print("--- %s seconds ---" % (time.time() - start_time))
