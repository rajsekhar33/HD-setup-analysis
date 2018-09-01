import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
import pyPLUTO as pp
#plt.style.use('classic')

params = {'legend.fontsize':16,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 16
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.major.size'] = 14
plt.rcParams['ytick.minor.size'] = 7
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(8.5, 8.)

time_step=0.2
no_bins=200
no_files=5

start=np.array((10., 51., 10., 42., 40., 15.))
wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'T-runs/TDh/', 'B-runs/BDh/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh')

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.

num_plots=4
spectra = [None] * (2)
ratio_spectra=[None]*start.size
for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	k1=np.zeros((no_files,no_bins))
	Ek1=np.zeros((no_files,no_bins))
	k2=np.zeros((no_files,no_bins))
	Ek2=np.zeros((no_files,no_bins))

	filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i]
	for j in xrange(no_files):
	   fileno=str(int(start[i]/time_step+0.5+j)).rjust(4,'0')
	   filename=filedir+'Rhoks'+str(fileno)+'.txt'
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k1[j]=data[:,0]
	   Ek1[j]=data[:,1]
	   filename=filedir+'Prsk'+str(fileno)+'.txt'
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k2[j]=data[:,0]
	   Ek2[j]=data[:,1]
	   all_data=pp.pload(int(start[i]/time_step+0.5+j), w_dir=filedir)
	   Ek2_mean=np.average(all_data.prs)
	   Ek2[j]=Ek2[j]/(Ek2_mean**2.)
#Calculate average energy injection
	k1=k1[0]
	k2=k2[0]
	
	ratiok=Ek2/Ek1
        del_ratiok=np.std(ratiok,0)
        ratio_k=np.average(ratiok,0)


#Calculate standard deviation and mean
	del_Ek1=np.std(Ek1,0)
	Ek1=np.average(Ek1,0)
	del_Ek2=np.std(Ek2,0)
	Ek2=np.average(Ek2,0)

	if (i==0):
                leg1=r'$A_k=\frac{|\rho_k|^2}{\left<\rho\right>^2}$'
                leg2=r'$A_k=\frac{|P_k|^2}{\left<P\right>^2}$'
	#Plot original power spectra
	#first velocity spectra
		spectra[2*i]=ax1.errorbar(k1[1:-10],Ek1[1:-10], yerr=del_Ek1[1:-10], color=colors[i], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=leg1, markersize=4., elinewidth=0.8)

		spectra[2*i+1]=ax1.errorbar(k1[1:-10],Ek2[1:-10], yerr=del_Ek2[1:-10], color=colors[i], fmt='v', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=leg2, markersize=4., elinewidth=0.8)
	elif (i<num_plots):
		ax1.errorbar(k1[1:-10],Ek1[1:-10], yerr=del_Ek1[1:-10], color=colors[i], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, markersize=4., elinewidth=0.8)

		ax1.errorbar(k1[1:-10],Ek2[1:-10], yerr=del_Ek2[1:-10], color=colors[i], fmt='v', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, markersize=4., elinewidth=0.8)

	#Plot ratio
	if (i<start.size):
		ratio_spectra[i]=ax2.errorbar(k1[1:-10], ratio_k[1:-10], yerr=del_ratiok[1:-10], color=colors[i], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=labels[i], markersize=4., elinewidth=1.)

x=np.arange(10., 10**3., 1.)
fit, =ax1.plot(x, 1e1*x**(-5./3.), label=r'$k^{-5/3}$', marker="d", markeredgecolor='none', markersize=0.5, linewidth=2.0, color=colors[9])
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_xlabel('k', fontsize=14)
ax1.set_ylabel('$A_k$', fontsize=20.)
#ax1.set_ylabel(r'$\frac{V_k^2}{c_s^2}$, $\frac{\rho_k^2}{\left<\rho\right>^2}$', fontsize=18)
ax1.set_ylim(1.e-6,0.2)
ax1.set_xlim(1.e1,1e3)
ax1.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.2)
spectra_legend=ax1.legend(handles=spectra, loc='lower left', bbox_to_anchor=(-0.05, -0.08), ncol=2, fancybox=True, framealpha=0., fontsize=28.)
ax1.add_artist(spectra_legend)
ax1.legend(handles=[fit], loc='upper right', bbox_to_anchor=(1., 1.0), ncol=1, fancybox=True, framealpha=0., fontsize=25.)

#fit, =ax2.plot(x, 4.*x**(-1./6.), label=r'$k^{-1/6}$', marker="d", markeredgecolor='none', markersize=0.5, linewidth=2.0, color=colors[9])
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('$k$', fontsize=20)
ax2.set_ylabel(r'$R_k$', fontsize=20)
ax2.set_ylim(1.e-3,3.)
ax2.set_xlim(1.e1,1e3)

ratio_legend=ax2.legend(handles=ratio_spectra, loc='lower left', bbox_to_anchor=(-0.02, -0.08), ncol=3, fancybox=True, framealpha=0., fontsize=20.)
ax2.add_artist(ratio_legend)
#ax2.legend(handles=[fit], loc='upper right', bbox_to_anchor=(.8, 1.0), ncol=1, fancybox=True, framealpha=0., fontsize=25.)

ax2.grid(color='grey', linestyle='-', linewidth=0.2)
ax2.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
plt.savefig('ratio-rho-prs-spectra-cool.png',dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))
