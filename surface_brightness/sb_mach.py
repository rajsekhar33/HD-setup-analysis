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
fig, ax = plt.subplots()
fig.set_size_inches(7, 4)

start=np.array((100, 75, 35, 20, 18, 14, 5, 4, 2))
amp=np.array((0.005, 0.02, 0.02, 0.10, 0.1, 0.1, 0.9, 0.9, 2.5))
mach=np.array((0.25, 0.42, 0.50, 0.75, 0.8, 0.9, 1.2, 1.5, 2.0))
time_step=0.2
no_bins=200
no_files=1
sb_dev=np.zeros((start.size))
sb_mean=np.zeros((start.size))
for i in xrange(start.size):
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*1000)).rjust(4,'0')+'/'
	data=np.fromfile(filedir+'sbs'+str(start[i]).rjust(4,'0')+'.dbl')
	sb_dev[i]=np.std(data)
	sb_mean[i]=np.average(data)
ax.scatter(mach, sb_dev/sb_mean, label=r'$\frac{\delta (SB)}{\left<SB\right>}$ vs $\left<\mathcal{M}\right>_{rms}$')
m1=np.arange(0.2, 0.8, 0.01)
m2=np.arange(0.8, 2.0, 0.01)
ax.plot(m1, 0.13*m1**2., label=r'$\left<\mathcal{M}\right>_{rms}^2$')
ax.plot(m2, .12*m2, label=r'$\left<\mathcal{M}\right>_{rms}$')
ax.legend(loc='lower center', ncol=3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.2,4.)
ax.set_ylim(0.002,.3)
ax.set_xlabel('$\mathcal{M}$')
ax.set_ylabel(r'$\frac{\delta SB}{SB}$')
plt.savefig('sb_mach.png', dpi=250)
plt.close()

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(7, 7)
for i in xrange(start.size):
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*1000)).rjust(4,'0')+'/'
        rhok=np.loadtxt(filedir+'Rhoks'+str(start[i]).rjust(4,'0')+'.txt')
        sbk =np.loadtxt(filedir+'sbks'+str(start[i]).rjust(4,'0')+'.txt')
	if (i%3<1e-1):
		ax1.errorbar(rhok[:,0][1:-10], rhok[:,1][1:-10], fmt='d', label=r'$\mathcal{M}=$'+str(mach[i]), markeredgecolor='None')
		ax2.errorbar(sbk[:,0][1:-10], sbk[:,1][1:-10], fmt='*', label=r'$\mathcal{M}=$'+str(mach[i]), markeredgecolor='None')

x=np.arange(10., 10**3., 1.)

ax1.plot(x, x**(-5./3.), label=r'$k^{-5/3}$')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel(r'$|\rho_k|^2$')
ax1.set_xlim(1e1, 1e3)
ax1.set_ylim(1e-11,)
ax1.legend(loc='lower left', ncol=4)

ax2.plot(x, 1e-8*x**(-8./3.), label=r'$k^{-8/3}$')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel(r'$|SB_k|^2$')
ax2.set_xlim(1e1,1e3)
ax2.legend(loc='lower left', ncol=4)
ax2.set_xlabel('k')
plt.savefig('SBk_rhok.png', dpi=250)
