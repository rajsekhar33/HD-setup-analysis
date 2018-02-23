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
end=45
time_step=0.2
no_bins=200
n=np.array((512,512,512))

k1=np.zeros((no_bins))
Ek1=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/cooling/power_law/512/k0-2/'
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
   k1=datax1[:,0]
   Ek1[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]

del_k1=k1[1:]-k1[:-1]
k1=k1[:-1]
Ek1=np.average(Ek1,0)
vk1=np.sqrt(2.0*Ek1[:-1]*del_k1)
tmix1=1.0/(k1*vk1)



k2=np.zeros((no_bins))
Ek2=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/cooling/power_law/512/k4-6/'
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
   k2=datax1[:,0]
   Ek2[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]

del_k2=k2[1:]-k2[:-1]
k2=k2[:-1]
Ek2=np.average(Ek2,0)
vk2=np.sqrt(2.0*Ek2[:-1]*del_k2)
tmix2=1.0/(k2*vk2)



k3=np.zeros((no_bins))
Ek3=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/cooling/power_law/512/k0-2/k6-8/'
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
   k3=datax1[:,0]
   Ek3[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]
del_k3=k3[1:]-k3[:-1]
k3=k3[:-1]
Ek3=np.average(Ek3,0)
vk3=np.sqrt(2.0*Ek3[:-1]*del_k3)
tmix3=1.0/(k3*vk3)



k4=np.zeros((end-start,no_bins))
Ek4=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/cooling/power_law/512/k8-10/'
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
   k4=datax1[:,0]
   Ek4[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]
del_k4=k4[1:]-k4[:-1]
k4=k4[:-1]
Ek4=np.average(Ek4,0)
vk4=np.sqrt(2.0*Ek4[:-1]*del_k4)
tmix4=1.0/(k4*vk4)




k5=np.zeros((end-start,no_bins))
Ek5=np.zeros((end-start,no_bins))
filedir='/mnt/lustre/ug4/ugrajs/cooling/power_law/512/k12/'
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
   k5=datax1[:,0]
   Ek5[filenumber-start]=datax1[:,1]+datax2[:,1]+datax3[:,1]
del_k5=k5[1:]-k5[:-1]
Ek5=np.average(Ek5,0)
vk5=np.sqrt(2.0*Ek5[:-1]*del_k5)
k5=k5[:-1]
tmix5=1.0/(k5*vk5)


#Plot the data 

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:-10],tmix1[1:-10],yerr=np.std(tmix1[1:-10]),fmt='*-',label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
ax.errorbar(k2[1:-10],tmix2[1:-10],yerr=np.std(tmix2[1:-10]),fmt='d-',label=r'$4 \leq |k_{driving}| \leq 6$')
ax.errorbar(k3[1:-10],tmix3[1:-10],yerr=np.std(tmix3[1:-10]),fmt='.-',label=r'$6 \leq |k_{driving}| \leq 8$')
ax.errorbar(k4[1:-10],tmix4[1:-10],yerr=np.std(tmix4[1:-10]),fmt='o-',label=r'$8 \leq |k_{driving}| \leq 10$')
ax.errorbar(k5[1:-10],tmix5[1:-10],yerr=np.std(tmix5[1:-10]),fmt='*-',label=r'$|k_{driving}|=12$')
ax.plot(k1[1:100],10*k1[1:100]**(-2./3.),'-',label='$k^{-2/3}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel('$t_{mix}(k)$')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('$t_{mix}(k)$ vs k for $512^3$, different $k_{driving}$' )
plt.savefig('tmix_k_diff_k_512.png',dpi=250)
print("--- %s seconds ---" % (time.time() - start_time))
