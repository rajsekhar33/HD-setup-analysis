import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
plt.style.use('classic')
params = {'legend.fontsize':10.,
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

start=np.array((83, 100, 75, 35, 20, 18, 14, 5, 4, 2))
amp=np.array((0.005, 0.005, 0.02, 0.02, 0.10, 0.1, 0.1, 0.9, 0.9, 2.5))
mach=np.array((0.21, 0.24, 0.43, 0.50, 0.75, 0.8, 0.9, 1.2, 1.5, 2.0))
time_step=0.2
no_bins=200
no_files=1
sb_dev=np.zeros((start.size))
sb_mean=np.zeros((start.size))

#fig, ax = plt.subplots()
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(7, 7)
for i in xrange(start.size):
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/'
        rhok=np.loadtxt(filedir+'Rhoks'+str(start[i]).rjust(4,'0')+'.txt')
        sbk =np.loadtxt(filedir+'sbks'+str(start[i]).rjust(4,'0')+'.txt')
	if (i%3<1e-1):
		ax1.errorbar(rhok[:,0][1:-10], rhok[:,1][1:-10], fmt='d', label=r'$\mathcal{M}=$'+str(mach[i]), markeredgecolor='None')
		ax2.errorbar(sbk[:,0][1:-10], sbk[:,1][1:-10], fmt='*', label=r'$\mathcal{M}=$'+str(mach[i]), markeredgecolor='None')

x=np.arange(10., 10**3., 1.)

ax1.plot(x, 5*x**(-5./3.), label=r'$k^{-5/3}$')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel(r'$|\rho_k|^2$', fontsize=12)
ax1.set_xlim(1e1, 1e3)
ax1.set_ylim(1e-11,)
ax1.legend(loc='lower left', ncol=5)

ax2.plot(x, 8e-8*x**(-8./3.), label=r'$k^{-8/3}$')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel(r'$|SB_k|^2$',  fontsize=12)
ax2.set_xlim(1e1,1e3)
ax2.legend(loc='lower left', ncol=5)
ax2.set_xlabel('k',  fontsize=12)
plt.savefig('SBk_rhok.png', dpi=250)
