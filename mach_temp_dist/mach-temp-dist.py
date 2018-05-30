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
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plot.rcParams.update(params)

plt.rc('text', usetex=True)

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

#wdir=('thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/')
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/')
#wdir=('no_turb/2e-1/', 'turb_perturb//DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')

#labels=('Bl', 'Bh')
labels=('Tl', 'Th')
#labels=('QD', 'TDh', 'BDh')

#num_plots=2
num_plots=2
#num_plots=3

distM = [None] * (num_plots)
distT = [None] * (num_plots)

#t_mp=(10., 42.)
t_mp=(10., 51.)
#t_mp=(10.4, 42., 15.)

#t_turb=(1., 4.8)
t_turb=(0.8, 4.8)
#t_turb=(0.2, 10.2, 3.)

step_size=np.array((0.2))#, 0.2, 0.2))

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1)
fig.set_size_inches(8, 9)

for i in xrange(num_plots):
	j=i
	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i]
	file=filedir+'mach'+str(int(t_turb[i]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	distM[j], =ax1.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), label=labels[i]+', t='+str(int(t_turb[i]*UNIT_TIME))+' Myr')
	file=filedir+'mach'+str(int(t_mp[i]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax1.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), dashes=[30, 5, 10, 5], )
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlim(1e-2, 1e1)
	ax1.set_ylim(1e-4, 4e-1)
	ax1.set_xlabel(r'$\mathcal{M}$', fontsize=20)
	ax1.set_ylabel(r'Volume fraction', fontsize=20)
	# Shrink current axis by 20%
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height])
	#plt.savefig('mach-dist.png',dpi=250)
	ax1.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
	ax1.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
	ax1.grid(color='grey', linestyle='-', linewidth=0.2)
	


	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i]
	file=filedir+'temp'+str(int(t_turb[i]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax2.plot(data[:,0], data[:,1]/(float(nx*ny*nz)))
	file=filedir+'temp'+str(int(t_mp[i]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	distT[j], =ax2.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), dashes=[30, 5, 10, 5], label=labels[i]+', t='+str(int(t_mp[i]*UNIT_TIME))+' Myr')
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.set_xlim(1e5, 1e8)
	ax2.set_ylim(1e-4, 2e-1)
	ax2.set_xlabel(r'$T(K)$', fontsize=20)
	ax2.set_ylabel(r'Volume fraction', fontsize=20)
	box = ax2.get_position()
	#plt.savefig('temp-dist.png',dpi=250)
	ax2.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
	ax2.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
	ax2.grid(color='grey', linestyle='-', linewidth=0.2)

distM_legend=ax1.legend(handles=distM, loc='upper left', ncol=2, fancybox=True, framealpha=0., bbox_to_anchor=(0., 1.05))
distT_legend=ax2.legend(handles=distT, loc='upper left', ncol=1, fancybox=True, framealpha=0., bbox_to_anchor=(0., 1.05))

#plt.savefig('mach-temp-dist-B.png',dpi=250, mode="expand", borderaxespad=0.)
plt.savefig('mach-temp-dist-T.png',dpi=250, mode="expand", borderaxespad=0.)
#plt.savefig('mach-temp-dist-D.png',dpi=250, mode="expand", borderaxespad=0.)

plt.close()

