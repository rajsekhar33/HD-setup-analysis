import numpy as np
import matplotlib.pyplot as plt
import pylab as plot

plt.style.use('classic')
params = {'legend.fontsize':9.0,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)


#Declare all parameters and filenames, file location

CONST_pc=3.0856775807e18
UNIT_VELOCITY= (1.e8)
UNIT_LENGTH =  (CONST_pc*40.e3)
UNIT_TIME=UNIT_LENGTH/UNIT_VELOCITY/(3.15e13)

no_files=7
#Load data files
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'no_turb/2e-1/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/')

labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')

fig, ax = plt.subplots()
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1, usecols=(0,13))
	ax.plot(data[:,0]*UNIT_TIME,data[:,1],label=labels[i1])
fig.set_size_inches(6, 5)
ax.set_xlabel(r'time (Myr)')
ax.set_ylabel(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$', fontsize=15)
ax.set_xlim(0,2200)
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ vs time')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), ncol=7)
plt.savefig('del-rho-time.png',dpi=250)
