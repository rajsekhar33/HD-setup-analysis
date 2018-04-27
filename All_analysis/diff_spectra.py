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

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(7, 7)

start=np.array((100, 75, 20, 14, 5, 2))
amp=np.array((0.005, 0.02, 0.1, 0.1, 0.9,2.5))
mach=((0.24, 0.43, 0.76, 0.90, 1.25, 2.1))
epsilon=((0.001255, 0.0088, 0.142, 0.133, 6.16, 25.45))
cs=((0.59, 0.69, 0.87, 0.75, 1.982, 1.97))
time_step=0.2
no_bins=200
no_files=1

for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	k1=np.zeros((no_files,no_bins))
	Ek1=np.zeros((no_files,no_bins))
	Ekcomp1=np.zeros((no_files,no_bins))
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*1000)).rjust(4,'0')+'/'
	file=filedir+'pluto_hst.out'
	data1 = np.loadtxt(file, skiprows=1, usecols=(0,10,14))
	for filenumber in xrange(start[i], start[i]+no_files):
	   fileno=str(filenumber).rjust(4,'0')
	   filenamex1=filedir+'Ekx'+str(fileno)+'.txt'
	   filenamex2=filedir+'Eky'+str(fileno)+'.txt'
	   filenamex3=filedir+'Ekz'+str(fileno)+'.txt'
	   datax1 = np.loadtxt(filenamex1, usecols=(0,1,2))
	   datax2 = np.loadtxt(filenamex2, usecols=(0,1,2))
	   datax3 = np.loadtxt(filenamex3, usecols=(0,1,2))
	   k1[filenumber-start[i]]=datax1[:,0]
	   Ek1[filenumber-start[i]]=datax1[:,1]+datax2[:,1]+datax3[:,1]
#Calculate average energy injection
	   epsilon= np.average(data[(data[:,0]>time_step*filenumber)*(data[:,0]<time_step*filenumber+1e-2)][:,11])
	   Ekcomp1[filenumber-start[i]]=epsilon[i]**(-2/3)*(datax1[:,2]+datax2[:,2]+datax3[:,2])

	k1=k1[0]
	del_Ek1=np.std(Ek1,0)
	Ek1=np.average(Ek1,0)
	Ek_comp1=np.average(Ekcomp1,0)
	del_Ek_comp1=np.std(Ekcomp1,0)

	k2=np.zeros((no_files,no_bins))
	Ek2=np.zeros((no_files,no_bins))
	Ekcomp2=np.zeros((no_files,no_bins))

	file=filedir+'pluto_hst.out'
	data2 = np.loadtxt(file, skiprows=1, usecols=(0,10))
	for filenumber in xrange(start[i], start[i]+no_files):
	   fileno=str(filenumber).rjust(4,'0')
	   filename=filedir+'Rhok'+str(fileno)+'.txt'
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k2[filenumber-start[i]]=data[:,0]
	   Ek2[filenumber-start[i]]=data[:,1]
#Calculate average energy injection
	   Ekcomp2[filenumber-start[i]]=epsilon[i]**(-2/3)*data[:,2]
	k2=k2[0]
	del_Ek2=np.std(Ek2,0)
	Ek2=np.average(Ek2,0)
	Ek_comp2=np.average(Ekcomp2,0)
	del_Ek_comp2=np.std(Ekcomp2,0)
#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
	if (i==2 or i==0): 
		ax1.errorbar(k1[1:-10],Ek_comp1[1:-10],yerr=del_Ek_comp1[1:-10],fmt='*', label=r'$\epsilon_v^{-2/3}k^{5/3}E_v(k)$, $\mathcal{M}=$'+str(mach[i]), markeredgecolor='None')
		ax1.errorbar(k2[1:-10],Ek_comp2[1:-10],yerr=del_Ek_comp2[1:-10],fmt='d', label=r'$\epsilon_v^{-2/3}k^{5/3}E_{\rho}(k)$, $\mathcal{M}= $'+str(mach[i]), markeredgecolor='None')
	#Plot ratio
	ax2.plot(k1[1:-10],Ek2[1:-10]*cs[i]**2./Ek1[1:-10], label='$\mathcal{M}=$'+str(mach[i]))
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_xlabel('k')
ax1.set_ylabel('$E(k)k^{5/3}\epsilon_V^{-2/3}$')
# Put a legend to the right of the current axis
ax1.legend(loc='lower right', bbox_to_anchor=(0.88, 0.0), ncol=2)
#ax1.set_title(r'Compensated $E_{\rho}(k)$ and $E_V(k)$ vs k for $256^3$' )
ax1.set_xlim(1e1, 2e3)

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('k')
ax2.set_ylabel(r'$\frac{\rho_k^2}{\left<\rho\right>^2}/\frac{V_k^2}{c_s^2}$')
#ax2.set_ylim(1.e-2,3e-1)
ax2.legend(loc='upper left', bbox_to_anchor=(0., 1.0), ncol=6)
#ax2.set_title(r'Ratio of velocity and density power spectra' )
ax2.grid(color='grey', linestyle='-', linewidth=0.2)
plt.savefig('Ratio_rho-vel_lowk_256.png',dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))
