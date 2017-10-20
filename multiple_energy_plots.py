import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats

#Compute how long the simulation takes
start_time = time.time()
#Declare all parameters and filenames, file location

#Load data files
nu1=0
nu2=1.e-3
filenumber=3;
n=np.array([512,512,512])
file1='/home/rajsekhar/PLUTO41_old/3D_turb/Tau_c_2/512/Ek0003.txt'
file2='/home/rajsekhar/PLUTO41_old/3D_turb/Tau_c_2/with_viscosity/512/3/Ek0003.txt'

fname1 = open(file1,'rt')
fname2 = open(file2,'rt')

data1 = np.loadtxt(file1, delimiter=' ', usecols=(0,1,2))
data2 = np.loadtxt(file2, delimiter=' ', usecols=(0,1,2))
#data3 = np.loadtxt(file3, delimiter=' ', usecols=(0,1,2))

k1=data1[:,0][0:150]
Ek1=data1[:,1][0:150]
Ek_comp1=data1[:,2][0:150]

k2=data2[:,0][0:150]
Ek2=data2[:,1][0:150]
Ek_comp2=data2[:,2][0:150]
'''
k3=data3[:,0]
Ek3=data3[:,1]
Ek_comp3=data3[:,2]
'''
#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
ax.plot(k1,Ek_comp1,'o-',label=r'$\nu = $'+str(nu1))
ax.plot(k2,Ek_comp2,'*-',label=r'$\nu = $'+str(nu2))
#ax.plot(k3,Ek_comp3,'d-',label=str(n3[0])+'$^3$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel('E(k)*$k^{5/3}$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
#plt.ylim(10**11,10**15)
plt.title('Compensated E(k) vs k for t='+str(float(filenumber))+' size '+str(n[0])+'$^3$' )

plt.savefig('E_k_compensated'+str(filenumber)+'_size_'+str(n[0])+'.png')

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
ax.plot(k1,Ek1,'o-',label=r'$\nu = $'+str(nu1))
ax.plot(k2,Ek2,'*-',label=r'$\nu = $'+str(nu2))
#ax.plot(k3,Ek3,'d-',label=str(n3[0])+'$^3$')
ax.plot(k1,k1**(-5/3),'d-',label='$k^{-5/3}$')

ax.set_yscale('log')
ax.set_xscale('log')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
plt.xlabel('k')
plt.ylabel('E(k)')
#plt.ylim(10**11,10**15)
plt.title('E(k) vs k for t='+str(float(filenumber))+' size '+str(n[0])+'$^3$')

plt.savefig('E_k'+str(filenumber)+'_size_'+str(n[0])+'.png')

print("--- %s seconds ---" % (time.time() - start_time))
