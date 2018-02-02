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

#Load data files
file1="test_power_spectrum_"+str(n1[0])+".txt"
file2="test_power_spectrum_"+str(n2[0])+".txt"

fname1 = open(file1,'rt')
fname2 = open(file2,'rt')

data1 = np.loadtxt(file1, delimiter=' ', usecols=(0,1))
data2 = np.loadtxt(file2, delimiter=' ', usecols=(0,1))
k1=data1[:,0]
Ek1=data1[:,1]

k2=data2[:,0]
Ek2=data2[:,1]

#Plot the data 

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
ax.loglog(k1,Ek1,label=str(n1[0])+'$^3$')
ax.loglog(k2,Ek2,label=str(n2[0])+'$^3$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.0))
plt.xlabel('k')
plt.ylabel('E(k)')
#plt.ylim(10**11,10**15)
plt.title('E(k) vs k')

plt.savefig('E_k_test.png')

print("--- %s seconds ---" % (time.time() - start_time))
