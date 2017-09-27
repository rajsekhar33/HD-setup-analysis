import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats

#Compute how long the simulation takes
start_time = time.time()

n1=np.array([256,256,256])
n2=np.array([128,128,128])

#Declare all parameters and filenames, file location

filedir="/home/rajsekhar/MHD-TURBULE-01/HD-setup/Data/"
filenumber=13

z=1.0
solver="hllc"
filedir1=filedir+"Z"+str(z)+'/'+str(n1[0])+'/'+solver+'/'
filedir2=filedir+"Z"+str(z)+'/'+str(n2[0])+'/'+solver+'/'

#Load data files
file1=filedir1+"power_spectrum_"+str(filenumber)+".txt"
file2=filedir2+"power_spectrum_"+str(filenumber)+".txt"

fname1 = open(file1,'rt')
fname2 = open(file2,'rt')

data1 = np.loadtxt(file1, delimiter=' ', usecols=(0,1,2))
data2 = np.loadtxt(file2, delimiter=' ', usecols=(0,1,2))
k1=data1[:,0]
Ek1=data1[:,1]
Ek_comp1=data1[:,2]

k2=data2[:,0]
Ek2=data2[:,1]
Ek_comp2=data2[:,2]

#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
ax.loglog(k1,Ek_comp1,label=str(n1[0])+'$^3$')
ax.loglog(k2,Ek_comp2,label=str(n2[0])+'$^3$')
plt.xlabel('k')
plt.ylabel('E(k)*$k^{5/3}$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
#plt.ylim(10**11,10**15)
plt.title('Compensated E(k) vs k for t='+str(filenumber/10))

plt.savefig('E_k_compensated'+str(filenumber)+'.png')

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
ax.loglog(k1,Ek1,label=str(n1[0])+'$^3$')
ax.loglog(k2,Ek2,label=str(n2[0])+'$^3$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
plt.xlabel('k')
plt.ylabel('E(k)')
#plt.ylim(10**11,10**15)
plt.title('E(k) vs k for t='+str(filenumber/10))

plt.savefig('E_k'+str(filenumber)+'.png')

print("--- %s seconds ---" % (time.time() - start_time))