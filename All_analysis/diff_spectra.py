import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
params = {'legend.fontsize':10,
          'legend.handlelength': 0.5}
plot.rcParams.update(params)



#Compute how long the simulation takes
start_time = time.time()
#Declare all parameters and filenames, file location

#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(6, 8)
for i in xrange(0, 4):
	#Load data files
	start=np.array((10, 2, 1, 1))
	amp=np.array((0.020,0.1,0.9,2.5))
	time_step=1.0
	no_bins=200
	no_files=5

	k1=np.zeros((no_files,no_bins))
	Ek1=np.zeros((no_files,no_bins))
	Ekcomp1=np.zeros((no_files,no_bins))
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*1000)).rjust(4,'0')+'/'
	file=filedir+'pluto_hst.out'
	fname = open(file,'rt')
	data1 = np.loadtxt(file, skiprows=1, usecols=(0,10))
	for filenumber in xrange(start[i], start[i]+no_files):
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
	   k1[filenumber-start[i]]=datax1[:,0]
	   Ek1[filenumber-start[i]]=datax1[:,1]+datax2[:,1]+datax3[:,1]
	   epsilon=np.average(data1[(data1[:,0]>time_step*filenumber)*(data1[:,0]<time_step*filenumber+1)])
	   Ekcomp1[filenumber-start[i]]=epsilon**(-2/3)*(datax1[:,2]+datax2[:,2]+datax3[:,2])

	k1=k1[0]
	del_Ek1=np.std(Ek1,0)
	Ek1=np.average(Ek1,0)
	Ek_comp1=np.average(Ekcomp1,0)
	del_Ek_comp1=np.std(Ekcomp1,0)

	k2=np.zeros((no_files,no_bins))
	Ek2=np.zeros((no_files,no_bins))
	Ekcomp2=np.zeros((no_files,no_bins))

	file=filedir+'pluto_hst.out'
	fname = open(file,'rt')
	data2 = np.loadtxt(file, skiprows=1, usecols=(0,11))
	for filenumber in xrange(start[i], start[i]+no_files):
	   fileno=str(filenumber).rjust(4,'0')
	   filename=filedir+'Rhok'+str(fileno)+'.txt'
	   fname = open(filename,'rt')
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k2[filenumber-start[i]]=data[:,0]
	   Ek2[filenumber-start[i]]=data[:,1]
	   Ekcomp2[filenumber-start[i]]=data[:,2]

	k2=k2[0]
	del_Ek2=np.std(Ek2,0)
	Ek2=np.average(Ek2,0)
	Ek_comp2=np.average(Ekcomp2,0)
	del_Ek_comp2=np.std(Ekcomp2,0)
	if (i==0): 
		ax1.errorbar(k1[1:-10],Ek_comp1[1:-10],yerr=del_Ek_comp1[1:-10],fmt='*-',label=r'$k^{5/3}V_k^2$, $A_{turb}=$'+str(amp[i]))
		ax1.errorbar(k2[1:-10],Ek_comp2[1:-10],yerr=del_Ek_comp2[1:-10],fmt='d-',label=r'$k^{5/3}\rho_k^2$, $A_{turb}=$'+str(amp[i]))
	#Plot ratio
	ax2.plot(k1[1:-10],Ek1[1:-10]/Ek2[1:-10], label='$A_{turb}=$'+str(amp[i]))
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('k')
ax1.set_ylabel('E(k)*$k^{5/3}*\epsilon_V^{-2/3}$')
# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(0.2, 0.95), ncol=2)
ax1.set_title(r'Compensated $E_{\rho}(k)$ and $E_V(k)$ vs k for $256^3$' )

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('k')
ax2.set_ylabel(r'$V_k^2/\rho_k^2$')
ax2.legend(loc='center left', bbox_to_anchor=(0, 0.95), ncol=4)
ax2.set_title(r'Ratio of velocity and density power spectra' )
plt.savefig('Ratio_rho-vel_lowk_256.png',dpi=250)

print("--- %s seconds ---" % (time.time() - start_time))
