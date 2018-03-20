import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
params = {'legend.fontsize':7.5,
          'legend.handlelength': 0.5}
plot.rcParams.update(params)

#Declare all parameters and filenames, file location

CONST_pc=3.0856775807e18
UNIT_VELOCITY= (1.e8)
UNIT_LENGTH =  (CONST_pc*40.e3)
UNIT_TIME=UNIT_LENGTH/UNIT_VELOCITY/(3.15e13)

no_files=10
#Load data files
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F1e-3/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-3/k12/','thermal_heating/256/tabulated_cooling/F1e-2/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-2/k12/', 'thermal_heating/256/tabulated_cooling/F1e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-1/k12/','thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'turb_perturb/F1e-3/k12/', 'turb_perturb/F1e-2/k12/', 'no_turb/')

labels=('Ckl', 'Ckh', 'HCkl4', 'HCkh4',  'HCkl3', 'HCkh3',  'HCkl2', 'HCkh2', 'HCkl1', 'DHCkh2', 'DHCkh1', 'DHC')

fig, ax = plt.subplots()
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1, usecols=(0,13))
	ax.plot(data[:,0]*UNIT_TIME,data[:,1],label=labels[i1])
fig.set_size_inches(6, 4)
ax.set_xlabel(r'time (Myr)')
ax.set_ylabel(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$')
ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ vs time')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='center left', bbox_to_anchor=(0.9, 0.4), ncol=1)
plt.savefig('del-rho-time.png',dpi=250)
