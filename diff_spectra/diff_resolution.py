import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
plt.style.use('classic')
params = {'legend.fontsize':8.25,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)


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
filedir1='/mnt/lustre/ug4/ugrajs/fiducial_runs/128/'
file=filedir1+'pluto_hst.out'
fname = open(file,'rt')
data1 = np.loadtxt(file, skiprows=1, usecols=(0,10))
#Take the average value of epsilon for the time we are considering
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir1+'Ek'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
#scan from the binned Ek values
   k1[filenumber-start]=data[:,0][0:150]
   Ek1[filenumber-start]=data[:,1][0:150]
   #Take average value of epsioln between two such  time steps
   epsilon1=np.average(data1[(data1[:,0]>time_step*filenumber)*(data1[:,0]<time_step*filenumber+1)])
   Ekcomp1[filenumber-start]=epsilon1**(-2/3)*data[:,2][0:150]

k1=k1[0]
#Take average over different file numbers
del_Ek1=np.std(Ek1,0)
Ek1=np.average(Ek1,0)
Ek_comp1=np.average(Ekcomp1,0)
del_Ek_comp1=np.std(Ekcomp1,0)
k2=np.zeros((end-start,150))
Ek2=np.zeros((end-start,150))
Ekcomp2=np.zeros((end-start,150))
filedir2='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/'
file=filedir2+'pluto_hst.out'
fname = open(file,'rt')
data2 = np.loadtxt(file, skiprows=1, usecols=(0,10))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir2+'Ek'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k2[filenumber-start]=data[:,0][0:150]
   Ek2[filenumber-start]=data[:,1][0:150]
   epsilon2=np.average(data2[(data2[:,0]>time_step*filenumber)*(data2[:,0]<time_step*filenumber+1)])
   Ekcomp2[filenumber-start]=epsilon2**(-2/3)*data[:,2][0:150]

k2=k2[0]
del_Ek2=np.std(Ek2,0)
Ek2=np.average(Ek2,0)
Ek_comp2=np.average(Ekcomp2,0)
del_Ek_comp2=np.std(Ekcomp2,0)

k3=np.zeros((end-start,150))
Ek3=np.zeros((end-start,150))
Ekcomp3=np.zeros((end-start,150))
filedir3='/mnt/lustre/ug4/ugrajs/fiducial_runs/512/'
file=filedir3+'pluto_hst.out'
fname = open(file,'rt')
data3 = np.loadtxt(file, skiprows=1, usecols=(0,10))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir3+'Ek'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k3[filenumber-start]=data[:,0][0:150]
   Ek3[filenumber-start]=data[:,1][0:150]
   epsilon3=np.average(data3[(data3[:,0]>time_step*filenumber)*(data3[:,0]<time_step*filenumber+1)])
   Ekcomp3[filenumber-start]=epsilon3**(-2/3)*data[:,2][0:150]

k3=k3[0]
del_Ek3=np.std(Ek3,0)
Ek3=np.average(Ek3,0)
Ek_comp3=np.average(Ekcomp3,0)
del_Ek_comp3=np.std(Ekcomp3,0)
"""
k4=np.zeros((end-start,150))
Ek4=np.zeros((end-start,150))
Ekcomp4=np.zeros((end-start,150))
filedir4='/mnt/lustre/ug4/ugrajs/fiducial_runs/1024/'
file=filedir4+'pluto_hst.out'
fname = open(file,'rt')
data4 = np.loadtxt(file, skiprows=1, usecols=(0,10))
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir4+'Ek'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1,2))
   k4[filenumber-start]=data[:,0][0:150]
   Ek4[filenumber-start]=data[:,1][0:150]
   epsilon4=np.average(data4[(data4[:,0]>time_step*filenumber)*(data4[:,0]<time_step*filenumber+1)])
   Ekcomp4[filenumber-start]=epsilon4**(-2/3)*data[:,2][0:150]

k4=k4[0]
Ek4=np.average(Ek4,0)
del_Ek4=np.std(Ek4,0)
Ek_comp4=np.average(Ekcomp4,0)
del_Ek_comp4=np.std(Ekcomp4,0)
"""
#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
yerr=del_Ek_comp3
ax.errorbar(k1,Ek_comp1,yerr=del_Ek_comp1,fmt='*-',label=r'$n=128^3$', markeredgecolor='None')
ax.errorbar(k2,Ek_comp2,yerr=del_Ek_comp2,fmt='d-',label=r'$n=256^3$', markeredgecolor='None')
ax.errorbar(k3,Ek_comp3,yerr=del_Ek_comp3,fmt='.-',label=r'$n=512^3$', markeredgecolor='None')
#ax.errorbar(k4,Ek_comp4,yerr=del_Ek_comp4'.-',label=r'$n=1024^3$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
ax.set_ylim(1e-5,)
plt.ylabel('$E(k) k^{5/3} \epsilon^{-2/3}$')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, .0), ncol=3)
plt.title('Compensated E(k) vs k' )
plt.savefig('E_k_compensated_diff_resolution_fiducial.png',dpi=150)

#This is to plot the original power spectrum, without any compensation

fig, ax = plt.subplots()
fig.set_size_inches(6, 5)
ax.errorbar(k1,Ek1,yerr=del_Ek1,fmt='*-',label=r'$n=128^3$', markeredgecolor='None')
ax.errorbar(k2,Ek2,yerr=del_Ek2,fmt='d-',label=r'$n=256^3$', markeredgecolor='None')
ax.errorbar(k3,Ek3,yerr=del_Ek3,fmt='.-',label=r'$n=512^3$', markeredgecolor='None')
#ax.errorbar(k4,Ek4,yerr=del_Ek_comp4,fmt='.-',label=r'$n=1024^3$')
ax.plot(k1,10.*k1**(-5/3)*epsilon1**(-2/3),'-',label='$k^{-5/3}$')
ax.set_ylim(1e-7,)
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('k')
plt.ylabel('E(k)')
leg = ax.legend(loc='lower center', bbox_to_anchor=(0.5, .0), ncol=4)
plt.title('E(k) vs k' )
plt.savefig('E_k_diff_resolution_fiducial.png',dpi=150)
print("--- %s seconds ---" % (time.time() - start_time))
