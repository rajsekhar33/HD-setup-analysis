import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats

#Compute how long the simulation takes
start_time = time.time()

n1=np.array([128,128,128])
n2=np.array([256,256,256])
n3=np.array([512,512,512])

#Declare all parameters and filenames, file location

filedir="/home/rajsekhar/PLUTO41_old/3D_turb/Tau_c_2/"
filenumber=3
bin_size=1
filedir1=filedir+str(n1[0])+'/'
filedir2=filedir+str(n2[0])+'/'
filedir3=filedir+str(n3[0])+'/'
#Load data files
file1=filedir1+"density_power_spectrum_"+str(filenumber)+'_bin_'+str(bin_size*100)+'_'+str(n1[0])+".txt"
file2=filedir2+"density_power_spectrum_"+str(filenumber)+'_bin_'+str(bin_size*100)+'_'+str(n2[0])+".txt"
file3=filedir3+"density_power_spectrum_"+str(filenumber)+'_bin_'+str(bin_size*100)+'_'+str(n3[0])+".txt"

fname1 = open(file1,'rt')
fname2 = open(file2,'rt')
fname3 = open(file3,'rt')

data1 = np.loadtxt(file1, delimiter=' ', usecols=(0,1,2))
data2 = np.loadtxt(file2, delimiter=' ', usecols=(0,1,2))
data3 = np.loadtxt(file3, delimiter=' ', usecols=(0,1,2))

k1=data1[:,0]
Ek1=data1[:,1]
Ek_comp1=data1[:,2]

k2=data2[:,0]
Ek2=data2[:,1]
Ek_comp2=data2[:,2]

k3=data3[:,0]
Ek3=data3[:,1]
Ek_comp3=data3[:,2]

#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
ax.plot(k1,Ek_comp1,'o-',label=str(n1[0])+'$^3$')
ax.plot(k2,Ek_comp2,'*-',label=str(n2[0])+'$^3$')
ax.plot(k3,Ek_comp3,'d-',label=str(n3[0])+'$^3$')

ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(r'$\rho_k^2*k^{5/3}$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
#plt.ylim(10**11,10**15)
plt.title(r'Compensated $\rho_k^2$ vs k for t='+str(float(filenumber))+' '+'bin size = '+str(bin_size))

plt.savefig('rho_k_compensated'+str(filenumber)+'bin_size'+str(int(bin_size*100))+'.png')

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
ax.plot(k1,Ek1,'o-',label=str(n1[0])+'$^3$')
ax.plot(k2,Ek2,'*-',label=str(n2[0])+'$^3$')
ax.plot(k3,Ek3,'d-',label=str(n3[0])+'$^3$')
ax.plot(k3,k3**(-5/3),'d-',label='$k^{-5/3}$')
ax.set_yscale('log')
ax.set_xscale('log')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
plt.xlabel('k')
plt.ylabel(r'$\rho(k)$')
#plt.ylim(10**11,10**15)
plt.title(r'$\rho_k^2$ vs k for t='+str(float(filenumber))+' '+'bin size = '+str(bin_size))

plt.savefig('rho_k'+str(filenumber)+'bin_size'+str(int(bin_size*100))+'.png')

print("--- %s seconds ---" % (time.time() - start_time))
