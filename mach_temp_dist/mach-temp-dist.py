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
nx,ny,nz=256,256,256
#no_files=13
#Load data files
#wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-2/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-3/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'thermal_heating/256/tabulated_cooling/F1e-1/k12/', 'thermal_heating/256/tabulated_cooling/F1e-2/k12/', 'thermal_heating/256/tabulated_cooling/F1e-3/k12/', 'no_turb/', 'turb_perturb/F1e-2/k12/', 'turb_perturb/F1e-3/k12/')
#labels=('Ckl', 'Ckh', 'HCkl1', 'HCkl2', 'HCkl3', 'HCkl4', 'HCkh1', 'HCkh2', 'HCkh3', 'HCkh4', 'DHC', 'DHCkh1', 'DHCkh2')
#t_mp=(2, 45, 4, 10, 14, 21, 30, 33, 51, 30, 79, 35, 36)
#t_turb=(1.2, 4.8, 2, 4, 6, 10.2, 6, 6, 12, 12, 2, 10.2, 10.2)
#step_size=(0.2, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 0.2)

no_files=2
wdir=('thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/')
labels=('HCkl', 'HCkh')
t_mp=(6., 30 )
t_turb=(2., 4.8)
step_size=(0.2, 0.2)

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1)
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'mach'+str(int(t_turb[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax1.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+' t='+str(int(t_turb[i1]*UNIT_TIME))+'Myr')
	file=filedir+'mach'+str(int(t_mp[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax1.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+' t='+str(int(t_mp[i1]*UNIT_TIME))+'Myr')
fig.set_size_inches(6, 7)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e-3, 5e1)
ax1.set_ylim(1e-6, )
ax1.set_xlabel(r'$\mathcal{M}$')
ax1.set_ylabel(r'Number fraction')
# Shrink current axis by 20%
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax1.legend(loc='lower left', bbox_to_anchor=(0.0, 0.0), ncol=3)
#plt.savefig('mach-dist.png',dpi=250)

for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'temp'+str(int(t_turb[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax2.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+', t='+str(int(t_turb[i1]*UNIT_TIME))+'Myr')
	file=filedir+'temp'+str(int(t_mp[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax2.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+', t='+str(int(t_mp[i1]*UNIT_TIME))+'Myr')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(1e4, 3e8)
ax2.set_ylim(1e-6, )
ax2.set_xlabel(r'$T$')
ax2.set_ylabel(r'Number fraction')
box = ax2.get_position()
ax2.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax2.legend(loc='lower left', bbox_to_anchor=(0.0, 0.0), ncol=3)
#plt.savefig('temp-dist.png',dpi=250)
plt.savefig('mach-temp-dist-HC.png',dpi=250, mode="expand", borderaxespad=0.)
plt.close()

