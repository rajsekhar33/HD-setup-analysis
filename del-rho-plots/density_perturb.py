import numpy as np
import matplotlib.pyplot as plt
import pylab as plot

#plt.style.use('classic')

params = {'legend.fontsize':20,
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

#Declare all parameters and filenames, file location

CONST_pc=3.0856775807e18
UNIT_VELOCITY= (1.e8)
UNIT_LENGTH =  (CONST_pc*40.e3)
UNIT_TIME=UNIT_LENGTH/UNIT_VELOCITY/(3.15e13)

no_files=7
#Load data files
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/', 'no_turb/2e-1/')

labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh', 'QD')

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.

fig, ax = plt.subplots()
fig.set_size_inches(7, 7)
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1, usecols=(0,13))
	ax.plot(data[:,0]*UNIT_TIME,data[:,1],label=labels[i1], color=colors[i1])
fig.set_size_inches(8, 6)
ax.set_xlabel(r'time (Myr)', fontsize=18)
ax.set_ylabel(r'$\left<\delta\rho\right>_{rms}/\left<\rho\right>$', fontsize=20)
#ax.set_ylabel(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$', fontsize=24)
ax.set_xlim(0,2200)
ax.set_ylim(0.,3.3)

ax.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)
ax.grid(color='grey', linestyle='-', linewidth=0.2)

#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ vs time')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='upper right', ncol=4, bbox_to_anchor=(1.0, 0.95), fancybox=True, framealpha=0.)
plt.savefig('del-rho-time.png',dpi=250)
