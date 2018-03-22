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
nx,ny,nz=256,256,256
#no_files=13
#Load data files
#wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-2/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-3/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'thermal_heating/256/tabulated_cooling/F1e-1/k12/', 'thermal_heating/256/tabulated_cooling/F1e-2/k12/', 'thermal_heating/256/tabulated_cooling/F1e-3/k12/', 'no_turb/', 'turb_perturb/F1e-2/k12/', 'turb_perturb/F1e-3/k12/')
#labels=('Ckl', 'Ckh', 'HCkl1', 'HCkl2', 'HCkl3', 'HCkl4', 'HCkh1', 'HCkh2', 'HCkh3', 'HCkh4', 'DHC', 'DHCkh1', 'DHCkh2')
#t_mp=(2, 45, 4, 10, 14, 21, 30, 33, 51, 30, 79, 35, 36)
#t_turb=(1.2, 4.8, 2, 4, 6, 10.2, 6, 6, 12, 12, 2, 10.2, 10.2)
#step_size=(0.2, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 0.2)

no_files=3
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'no_turb/')
labels=('Ckl', 'Ckh', 'DHC')
t_mp=(19, 45, 79)
t_turb=(1.2, 4.8, 2)
step_size=(0.2, 0.2, 0.2)

fig, ax = plt.subplots()
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'mach'+str(int(t_turb[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+' t='+str(int(t_turb[i1]*UNIT_TIME))+'Myr')
	file=filedir+'mach'+str(int(t_mp[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+' t='+str(int(t_mp[i1]*UNIT_TIME))+'Myr')
fig.set_size_inches(6, 4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\left<\mathcal{M}\right>_{rms}$')
ax.set_ylabel(r'Number density')
ax.set_title(r'Number density vs $\left<\mathcal{M}\right>_{rms}$')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.4), ncol=1)
#plt.savefig('mach-dist.png',dpi=250)
plt.savefig('mach-dist-extremes.png',dpi=250)
plt.close()

fig, ax = plt.subplots()
for i1 in xrange(0,no_files):
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	file=filedir+'temp'+str(int(t_turb[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+', t='+str(int(t_turb[i1]*UNIT_TIME))+'Myr')
	file=filedir+'temp'+str(int(t_mp[i1]/step_size[i1]+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax.plot(data[:,0],data[:,1]/(float(nx*ny*nz)),label=labels[i1]+', t='+str(int(t_mp[i1]*UNIT_TIME))+'Myr')
fig.set_size_inches(6, 4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\left<T\right>$')
ax.set_ylabel(r'Number density')
ax.set_title(r'Number density vs $\left<T\right>$')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.4), ncol=1)
#plt.savefig('temp-dist.png',dpi=250)
plt.savefig('temp-dist-extremes.png',dpi=250)
plt.close()

