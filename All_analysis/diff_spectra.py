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
start=20
end=30
time_step=1.0
no_bins=200

k1=np.zeros((end-start,no_bins))
Ek1=np.zeros((end-start,no_bins))
Ekcomp1=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp0020/'
file=filedir+'pluto_hst.out'
fname = open(file,'rt')
data1 = np.loadtxt(file, skiprows=1, usecols=(0,11))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filenamex1=filedir+'Ekx'+str(fileno)+'.txt'
   filenamex2=filedir+'Eky'+str(fileno)+'.txt'
   filenamex3=filedir+'Ekz'+str(fileno)+'.txt'
   fname = open(filenamex1,'rt')
   datax1 = np.loadtxt(filenamex1, usecols=(0,1,2))
   fname = open(filenamex2,'rt')
   datax2 = np.loadtxt(filenamex2, usecols=(0,1,2))
   fname = open(filenamex3,'rt')
   datax3 = np.loadtxt(filenamex3, usecols=(0,1,2))
   k1[filenumber-start]=datax1[:,0]
   Ek1[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]
   epsilon=np.average(data1[(data1[:,0]>time_step*filenumber)*(data1[:,0]<time_step*filenumber+1)])
   Ekcomp1[filenumber-start]=epsilon**(-2/3)*(datax1[:,2]+datax2[:,2]+datax3[:,2])

k1=k1[0]
del_Ek1=np.std(Ek1,0)
Ek1=np.average(Ek1,0)
Ek_comp1=np.average(Ekcomp1,0)
del_Ek_comp1=np.std(Ekcomp1,0)

k2=np.zeros((end-start,no_bins))
Ek2=np.zeros((end-start,no_bins))
Ekcomp2=np.zeros((end-start,no_bins))
#filedir2='/mnt/lustre/ug4/ugrajs/higher_k/512/k4-6/'
filedir2='/mnt/lustre/ug4/ugrajs/cooling/power_law/256/k0-2/'

#filedir2='/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k4-6/'
file=filedir2+'pluto_hst.out'
fname = open(file,'rt')
data2 = np.loadtxt(file, skiprows=1, usecols=(0,11))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir2+'Rhok'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k2[filenumber-start]=data[:,0]
   Ek2[filenumber-start]=data[:,1]
   Ekcomp2[filenumber-start]=data[:,2]

k2=k2[0]
del_Ek2=np.std(Ek2,0)
Ek2=np.average(Ek2,0)
Ek_comp2=np.average(Ekcomp2,0)
del_Ek_comp2=np.std(Ekcomp2,0)


#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:-10],Ek_comp1[1:-10],yerr=del_Ek_comp1[1:-10],fmt='*-',label=r'$V_k^2$')
ax.errorbar(k2[1:-10],Ek_comp2[1:-10],yerr=del_Ek_comp2[1:-10],fmt='d-',label=r'$\rho_k^2$')
'''
ax.errorbar(k3[1:-10],Ek_comp3[1:-10],yerr=del_Ek_comp3[1:-10],fmt='.-',label=r'$6 \leq |k_{driving}| \leq 8$')
ax.errorbar(k4[1:-10],Ek_comp4[1:-10],yerr=del_Ek_comp4[1:-10],fmt='o-',label=r'$8 \leq |k_{driving}| \leq 10$')
ax.errorbar(k5[1:-10],Ek_comp5[1:-10],yerr=del_Ek_comp5[1:-10],fmt='*-',label=r'$|k_{driving}|=12$')
'''
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel('E(k)*$k^{5/3}*\epsilon_V^{-2/3}$')
# Shrink current axis by 20%
box = ax.get_position()
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.8))
plt.title(r'Compensated $E_{\rho}(k)$ and $E_V(k)$ vs k for $256^3$' )
plt.savefig('E_k_compensated_rho-vel_lowk_256.png',dpi=250)


#Plot ratio
fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
plt.plot(k1[1:-10],Ek1[1:-10]/Ek2[1:-10])
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(r'$V_k^2/\rho_k^2$')
plt.ylim()
plt.title(r'Ratio of density and velocity power spectra' )
plt.savefig('Ratio_rho-vel_lowk_256.png',dpi=250)

print("--- %s seconds ---" % (time.time() - start_time))
