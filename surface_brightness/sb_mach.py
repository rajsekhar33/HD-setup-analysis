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

fig, ax = plt.subplots()
fig.set_size_inches(7., 7.)

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
ax.scatter(mach, sb_dev/sb_mean, label=r'$\frac{\delta (SB)}{\left<SB\right>}$ vs $\left<\mathcal{M}\right>_{rms}$')
m1=np.arange(0.21, 0.8, 0.01)
m2=np.arange(0.8, 2.0, 0.01)
ax.plot(m1, 0.24*m1**2., label=r'$\left<\mathcal{M}\right>_{rms}^2$')
ax.plot(m2, .24*m2, label=r'$\left<\mathcal{M}\right>_{rms}$')
ax.legend(loc='upper center', ncol=3, fancybox=True, framealpha=0., fontsize=15)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,3.)
ax.set_ylim(0.009,1.)
ax.set_xlabel('$\mathcal{M}$',  fontsize=12)
ax.set_ylabel(r'$\frac{\delta (SB)}{SB}$',  fontsize=14)

ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')

plt.savefig('sb_mach.png', dpi=250)
plt.close()

plt.close()
