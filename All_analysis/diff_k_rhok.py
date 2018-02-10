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

k1=np.zeros((end-start,no_bins))
Ek1=np.zeros((end-start,no_bins))
Ekcomp1=np.zeros((end-start,no_bins))
#filedir1='/mnt/lustre/ug4/ugrajs/fiducial_runs/512/'
filedir1='/mnt/lustre/ug4/ugrajs/cooling/512/'
file=filedir1+'pluto_hst.out'
fname = open(file,'rt')
data1 = np.loadtxt(file, skiprows=1, usecols=(0,10))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir1+'Rhok'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k1[filenumber-start]=data[:,0]
   Ek1[filenumber-start]=data[:,1]
   Ekcomp1[filenumber-start]=data[:,2]

k1=k1[0]
del_Ek1=np.std(Ek1,0)
Ek1=np.average(Ek1,0)
Ek_comp1=np.average(Ekcomp1,0)
del_Ek_comp1=np.std(Ekcomp1,0)
k2=np.zeros((end-start,no_bins))
Ek2=np.zeros((end-start,no_bins))
Ekcomp2=np.zeros((end-start,no_bins))
#filedir2='/mnt/lustre/ug4/ugrajs/higher_k/512/k4-6/'
filedir2='/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k4-6/'
file=filedir2+'pluto_hst.out'
fname = open(file,'rt')
data2 = np.loadtxt(file, skiprows=1, usecols=(0,10))
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

k3=np.zeros((end-start,no_bins))
Ek3=np.zeros((end-start,no_bins))
Ekcomp3=np.zeros((end-start,no_bins))
#filedir3='/mnt/lustre/ug4/ugrajs/higher_k/512/k6-8/'
filedir3='/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k6-8/'
file=filedir3+'pluto_hst.out'
fname = open(file,'rt')
data3 = np.loadtxt(file, skiprows=1, usecols=(0,10))
epsilon3=np.average(data3[(data3[:,0]>start*time_step)*(data3[:,0]<end*time_step)])
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir3+'Rhok'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k3[filenumber-start]=data[:,0]
   Ek3[filenumber-start]=data[:,1]
   Ekcomp3[filenumber-start]=data[:,2]
k3=k3[0]
del_Ek3=np.std(Ek3,0)
Ek3=np.average(Ek3,0)
Ek_comp3=np.average(Ekcomp3,0)
del_Ek_comp3=np.std(Ekcomp3,0)

k4=np.zeros((end-start,no_bins))
Ek4=np.zeros((end-start,no_bins))
Ekcomp4=np.zeros((end-start,no_bins))
filedir4='/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k8-10/'
file=filedir4+'pluto_hst.out'
fname = open(file,'rt')
data4 = np.loadtxt(file, skiprows=1, usecols=(0,10))
epsilon4=np.average(data4[(data4[:,0]>start*time_step)*(data4[:,0]<end*time_step)])
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir4+'Rhok'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k4[filenumber-start]=data[:,0]
   Ek4[filenumber-start]=data[:,1]
   Ekcomp4[filenumber-start]=data[:,2]
k4=k4[0]
del_Ek4=np.std(Ek4,0)
Ek4=np.average(Ek4,0)
Ek_comp4=np.average(Ekcomp4,0)
del_Ek_comp4=np.std(Ekcomp4,0)

k5=np.zeros((end-start,no_bins))
Ek5=np.zeros((end-start,no_bins))
Ekcomp5=np.zeros((end-start,no_bins))
filedir5='/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k12/'
file=filedir5+'pluto_hst.out'
fname = open(file,'rt')
data5 = np.loadtxt(file, skiprows=1, usecols=(0,10))
epsilon5=np.average(data5[(data5[:,0]>start*time_step)*(data5[:,0]<end*time_step)])
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir5+'Rhok'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k5[filenumber-start]=data[:,0]
   Ek5[filenumber-start]=data[:,1]
   Ekcomp5[filenumber-start]=data[:,2]
k5=k5[0]
del_Ek5=np.std(Ek5,0)
Ek5=np.average(Ek5,0)
Ek_comp5=np.average(Ekcomp5,0)
del_Ek_comp5=np.std(Ekcomp5,0)



#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:-10],Ek_comp1[1:-10],yerr=del_Ek_comp1[1:-10],fmt='*-',label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
ax.errorbar(k2[1:-10],Ek_comp2[1:-10],yerr=del_Ek_comp2[1:-10],fmt='d-',label=r'$4 \leq |k_{driving}| \leq 6$')
ax.errorbar(k3[1:-10],Ek_comp3[1:-10],yerr=del_Ek_comp3[1:-10],fmt='.-',label=r'$6 \leq |k_{driving}| \leq 8$')
ax.errorbar(k4[1:-10],Ek_comp4[1:-10],yerr=del_Ek_comp4[1:-10],fmt='o-',label=r'$8 \leq |k_{driving}| \leq 10$')
ax.errorbar(k5[1:-10],Ek_comp5[1:-10],yerr=del_Ek_comp5[1:-10],fmt='*-',label=r'$|k_{driving}|=12$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(r'$\rho$(k)*$k^{5/3}$')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.title(r'Compensated $\rho$(k) vs k for $512^3$, different $k_{driving}$' )
plt.savefig('Rho_k_compensated_diff_k_512.png',dpi=250)

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1[1:-10],Ek1[1:-10],yerr=del_Ek1[1:-10],fmt='*-',label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
ax.errorbar(k2[1:-10],Ek2[1:-10],yerr=del_Ek2[1:-10],fmt='d-',label=r'$4 \leq |k_{driving}| \leq 6$')
ax.errorbar(k3[1:-10],Ek3[1:-10],yerr=del_Ek3[1:-10],fmt='.-',label=r'$6 \leq |k_{driving}| \leq 8$')
ax.errorbar(k4[1:-10],Ek4[1:-10],yerr=del_Ek4[1:-10],fmt='o-',label=r'$8 \leq |k_{driving}| \leq 10$')
ax.errorbar(k5[1:-10],Ek5[1:-10],yerr=del_Ek5[1:-10],fmt='*-',label=r'$|k_{driving}|=12$')
ax.plot(k1[1:-10],k1[1:-10]**(-5/3),'-',label='$k^{-5/3}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel(r'$\rho$(k)')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title(r'$\rho$(k) vs k for $512^3$, different $k_{driving}$' )
plt.savefig('Rho_k_diff_k_512.png',dpi=250)
print("--- %s seconds ---" % (time.time() - start_time))
