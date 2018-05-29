import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
#plt.style.use('classic')
params = {'legend.fontsize':16,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.labelsize'] = 14 
plt.rcParams['ytick.labelsize'] = 14 
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(2.8, 3)

amp=np.array((0.005, 0.02, 0.1, 0.1, 0.9,2.5))
mach=((0.25, 0.45, 0.75, 0.90, 1.2, 2.1))
cs=((0.59, 0.69, 0.87, 0.75, 1.982, 1.97))
time_step=0.2
no_bins=200
no_files=5

start_array1=np.array(((1,27.0), (2,30.8), (3,20.4), (4,25.2), (5,17.0), (5,25.6)))
start_array2=np.array(((1,10.8), (1,12.8), (1,17.6), (1,20.4), (2,12.8)))
start_array3=np.array(((1,4.4), (1,5.4), (2,5.4), (3,5.6), (4,5.2)))
start_array4=np.array(((1,2.8), (2,3.2), (3,3.0), (4,4.2), (5,3.4)))
start_array5=np.array(((1,1.0), (2,1.4), (3,1.0), (4,1.2), (5,1.6)))
start_array6=np.array(((1,0.4), (2,0.4), (3,0.4), (4,0.4), (5,0.4)))

start=np.array((start_array1, start_array2, start_array3, start_array4, start_array5, start_array6))

colors=((230, 25, 75), (250, 190, 190) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (210, 245, 60), (145, 30, 180), (0, 128, 128), (240, 50, 230))
colors=np.array(colors)/255.

num_plots=3
spectra = [None] * (num_plots*2)

for i in xrange(0, amp.size):
	#Load data files
	#Declare all parameters and filenames, file location

	k1=np.zeros((no_files,no_bins))
	Ek1=np.zeros((no_files,no_bins))
	for j in xrange(no_files):
	   filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run'+str(int(start[i][:,0][j]))+'/'
	   file=filedir+'pluto_hst.out'
	   data1 = np.loadtxt(file, skiprows=1, usecols=(0,10,14))
	   fileno=str(int(start[i][:,1][j]/time_step)).rjust(4,'0')
	   filenamex1=filedir+'Ekx'+str(fileno)+'.txt'
	   filenamex2=filedir+'Eky'+str(fileno)+'.txt'
	   filenamex3=filedir+'Ekz'+str(fileno)+'.txt'
	   datax1 = np.loadtxt(filenamex1, usecols=(0,1,2))
	   datax = np.loadtxt(filenamex2, usecols=(0,1,2))
	   datax3 = np.loadtxt(filenamex3, usecols=(0,1,2))
	   k1[j]=datax1[:,0]
	   Ek1[j]=datax1[:,1]+datax[:,1]+datax3[:,1]

	k1=k1[0]

	k2=np.zeros((no_files,no_bins))
	Ek2=np.zeros((no_files,no_bins))

	for j in xrange(no_files):
	   filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run'+str(int(start[i][:,0][j]))+'/'
	   fileno=str(int(start[i][:,1][j]/time_step)).rjust(4,'0')
	   filename=filedir+'Rhoks'+str(fileno)+'.txt'
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k2[j]=data[:,0]
	   Ek2[j]=data[:,1]
#Calculate average energy injection
	k2=k2[0]
	
	ratiok=Ek2*cs[i]**2./Ek1
        del_ratiok=np.std(ratiok,0)
        ratio_k=np.average(ratiok,0)


#Calculate standard deviation and mean
	del_Ek1=np.std(Ek1,0)
	Ek1=np.average(Ek1,0)
	del_Ek2=np.std(Ek2,0)
	Ek2=np.average(Ek2,0)

	#Plot original power spectra
	if (i<num_plots):
	#first velocity spectra
		spectra[2*i]=ax1.errorbar(k1[1:-10],Ek1[1:-10]/cs[i]**2., yerr=del_Ek1[1:-10]/cs[i]**2., color=colors[i], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=r'$\frac{V_k^2}{c_s^2}$, $\mathcal{M}=$'+str(mach[i]), markersize=4., elinewidth=0.8)

	#pressure spectra
		spectra[2*i+1]=ax1.errorbar(k1[1:-10],Ek2[1:-10], yerr=del_Ek2[1:-10], color=colors[i], fmt='v', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=r'$\frac{\rho_k^2}{\left<\rho\right>^2}$, $\mathcal{M}=$'+str(mach[i]), markersize=4., elinewidth=0.8)

	#Plot ratio
	if (i<5):
		ax2.errorbar(k1[1:-10], ratio_k[1:-10], yerr=del_ratiok[1:-10], color=colors[i], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label='$\mathcal{M}=$'+str(mach[i]), markersize=2.5, elinewidth=0.8)

x=np.arange(10., 10**3., 1.)
fit, =ax1.plot(x, 1e-1*x**(-5./3.), label=r'$k^{-5/3}$', marker="d", markeredgecolor='none', markersize=0.5, linewidth=1.0)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_xlabel('k', fontsize=14)
ax1.set_ylabel(r'$\frac{V_k^2}{c_s^2}$, $\frac{\rho_k^2}{\left<\rho\right>^2}$', fontsize=18)
ax1.set_ylim(1.e-11,1e-1)
ax1.set_xlim(1.e1,1e3)
ax1.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.2)
spectra_legend=ax1.legend(handles=spectra, loc='lower left', bbox_to_anchor=(0., 0.0), ncol=3, fancybox=True, framealpha=0., fontsize=16.)
ax1.add_artist(spectra_legend)
ax1.legend(handles=[fit], loc='upper right', bbox_to_anchor=(1., 1.0), ncol=1, fancybox=True, framealpha=0., fontsize=16.)

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('k', fontsize=18)
ax2.set_ylabel(r'$\frac{V_k^2}{c_s^2}/\frac{\rho_k^2}{\left<\rho\right>^2}$', fontsize=18)
#ax2.set_ylim(1.e-10,1.)
ax2.set_xlim(1.e1,1e3)
ax2.legend(loc='lower right',  ncol=3, fancybox=True, framealpha=0.)
#ax.set_title(r'Ratio of velocity and density power spectra' )
ax2.grid(color='grey', linestyle='-', linewidth=0.2)
ax2.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
plt.savefig('ratio-rho-vel-spectra.png',dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))
