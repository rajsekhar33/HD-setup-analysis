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
start=40
end=50
time_step=0.2
k1=np.zeros((end-start,150))
Ek1=np.zeros((end-start,150))
Ekcomp1=np.zeros((end-start,150))
filedir1='/mnt/lustre/ug4/ugrajs/cooling/512/'
file=filedir1+'pluto_hst.out'
fname = open(file,'rt')
data1 = np.loadtxt(file, skiprows=1, usecols=(0,10))
#i=0 for density, i=1 for pressure
i=1
if (i==0):string='rho'
else:string='P'
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir1+'del'+string+'k'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k1[filenumber-start]=data[:,0][0:150]
   Ek1[filenumber-start]=data[:,1][0:150]
   epsilon1=np.average(data1[(data1[:,0]>time_step*filenumber)*(data1[:,0]<time_step*filenumber+1)])
   if(i==0):Ekcomp1[filenumber-start]=epsilon1**(-2.0/3.0)*data[:,2][0:150]
   else:Ekcomp1[filenumber-start]=epsilon1**(-4.0/3.0)*data[:,2][0:150]
k1=k1[0]
del_Ek1=np.std(Ek1,0)
Ek1=np.average(Ek1,0)
Ek_comp1=np.average(Ekcomp1,0)
del_Ek_comp1=np.std(Ekcomp1,0)
#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:],Ek_comp1[1:],yerr=del_Ek_comp1[1:],fmt='*-',label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(string+' Compensated E(k) '+string)
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 0.95))
plt.title(string+' Compensated E(k) '+string+' vs k for $512^3$' )
plt.savefig(string+'k_compensated_512.png',dpi=250)

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:],Ek1[1:],yerr=del_Ek1[1:],fmt='*-',label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
ax.plot(k1[1:],k1[1:]**(-7.0/3.0)*epsilon1**(-4.0/3.0),'-',label='$k^{-7/3}$')
ax.plot(k1[1:],k1[1:]**(-5.0/3.0)*epsilon1**(-2.0/3.0),'-',label='$k^{-5/3}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(string+' E(k)')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 0.95))
plt.title(string+r' E(k) vs k for $512^3$' )
plt.savefig(string+'k_512.png',dpi=250)
print("--- %s seconds ---" % (time.time() - start_time))
