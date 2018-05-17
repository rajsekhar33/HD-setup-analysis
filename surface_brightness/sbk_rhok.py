import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
import pyPLUTO as pp

#plt.style.use('classic')
params = {'legend.fontsize':10.,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

colors=((230, 25, 75), (250, 190, 190) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (210, 245, 60), (145, 30, 180), (0, 128, 128), (240, 50, 230))
colors=np.array(colors)/255.

num_plots=4
spectra = [None] * (num_plots*2)
fit = [None] * (2)


#Initialise the figure


start=np.array((83, 100, 75, 35, 20, 18, 14, 5, 4, 2))
amp=np.array((0.005, 0.005, 0.02, 0.02, 0.10, 0.1, 0.1, 0.9, 0.9, 2.5))
mach=np.array((0.21, 0.24, 0.43, 0.50, 0.75, 0.8, 0.9, 1.2, 1.5, 2.0))
time_step=0.2
no_bins=200
no_files=2
sb_dev=np.zeros((start.size))
sb_mean=np.zeros((start.size))
for i in xrange(start.size):
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run1/'
	data=np.fromfile(filedir+'sbs'+str(start[i]).rjust(4,'0')+'.dbl')
	sb_dev[i]=np.std(data)
	sb_mean[i]=np.average(data)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(9, 7)

start=np.array((100, 75, 20, 14, 5, 2))
amp=np.array((0.005, 0.02, 0.1, 0.1, 0.9,2.5))
mach=((0.24, 0.43, 0.76, 0.90, 1.25, 2.1))

for i in xrange(start.size):
	#if (i%3<1e-1):
	if (i<5):
		#j=int(i/3)
		j=i
		filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run1/'
		rhok=np.loadtxt(filedir+'Rhoks'+str(start[i]).rjust(4,'0')+'.txt')
		sbk =np.loadtxt(filedir+'sbks'+str(start[i]).rjust(4,'0')+'.txt')

		if (j<4):
			spectra[2*j]=ax1.errorbar(rhok[:,0][1:-20], rhok[:,1][1:-20], fmt='d', color=colors[j], markeredgecolor=None, markersize=5.0, ecolor=None, capsize=None, barsabove=False, label=r'$\frac{|\rho_k|^2}{\left<\rho\right>^2}$, $\mathcal{M}=$'+str(mach[i]))
			spectra[2*j+1]=ax1.errorbar(sbk[:,0][1:-20], sbk[:,1][1:-20]/sb_mean[i]**2., fmt='*', color=colors[j], markeredgecolor=None, markersize=5.0, ecolor=None, capsize=None, barsabove=False, label=r'$\frac{|SB_k|^2}{\left<SB\right>^2}$, $\mathcal{M}=$'+str(mach[i]))

		ax2.errorbar(sbk[:,0][1:-20], rhok[:,1][1:-20]/sbk[:,1][1:-20]/sb_mean[i]**2./sbk[:,0][1:-20], fmt='o', color=colors[j], markeredgecolor=None, markersize=5.0, ecolor=None, capsize=None, barsabove=False, label= '$\mathcal{M}=$'+str(mach[i]))


ax1.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.2)
spectra_legend=ax1.legend(handles=spectra, loc='lower left', bbox_to_anchor=(0.0, 0.0), ncol=4)
ax1.add_artist(spectra_legend)
spectra_legend.get_frame().set_alpha(0.)

x=np.arange(10., 10**3., 1.)
fit[0], = ax1.plot(x, 1e-3*x**(-5./3.), label=r'$k^{-5/3}$')
fit[1], = ax1.plot(x, 8e-8*x**(-8./3.), label=r'$k^{-8/3}$')
fit_legend=ax1.legend(handles=fit, loc='upper right', bbox_to_anchor=(1., 1.0), ncol=2)
fit_legend.get_frame().set_alpha(0.)

ax1.add_artist(fit_legend)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel(r'$\frac{|\rho_k|^2}{\left<\rho\right>^2}$, $\frac{|SB_k|^2}{\left<SB\right>^2}$', fontsize=12)
ax1.set_xlim(1e1, 1e3)
ax1.set_ylim(1e-20,1e-2)

ax2.set_ylabel(r'$\left(\frac{|\rho_k|^2}{\left<\rho\right>^2}\right)/\left(k\frac{|SB_k|^2}{\left<SB\right>^2}\right)$', fontsize=12)
ax2.set_xlabel('$k$', fontsize=12)
ax2.set_ylim(1e8,5e10)
ax2.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax2.grid(color='grey', linestyle='-', linewidth=0.2)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.legend(loc='upper right', bbox_to_anchor=(1., 1.0), ncol=3, fancybox=True, framealpha=0.)

#plt.show()
plt.savefig('SBk_rhok.png', dpi=250)
plt.close()
