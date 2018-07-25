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
wdir='F1e-1/'#, 'F5e-1/', 'F9e-1/')

labels='0.1'#, '0.5', '0.9')

num_plots=1

distM = [None] * (num_plots)
distT = [None] * (num_plots)

t_mp=10.
t_turb=1.
step_size=np.array((0.2))#, 0.2, 0.2))

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1)
fig.set_size_inches(8, 10.)

for i in xrange(num_plots):
	j=i
	filedir='/mnt/lustre/ug4/ugrajs/cooling/turb_perturb/DBh2e-1/'+wdir
	file=filedir+'mach'+str(int(t_turb/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	distM[j], =ax1.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), label='$f=$'+labels+', t='+str(int(t_turb*UNIT_TIME))+' Myr')
	file=filedir+'mach'+str(int(t_mp/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax1.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), dashes=[30, 5, 10, 5], )

	file=filedir+'temp'+str(int(t_turb/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	ax2.plot(data[:,0], data[:,1]/(float(nx*ny*nz)))
	file=filedir+'temp'+str(int(t_mp/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1, usecols=(0,1))
	distT[j], =ax2.plot(data[:,0], data[:,1]/(float(nx*ny*nz)), dashes=[30, 5, 10, 5], label='$f=$'+labels+', t='+str(int(t_mp*UNIT_TIME))+' Myr')

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

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(1e5, 1e8)
ax2.set_ylim(1e-4, 2e-1)
ax2.set_xlabel(r'$T(K)$', fontsize=20)
ax2.set_ylabel(r'Volume fraction', fontsize=20)
box = ax2.get_position()
#plt.savefig('temp-dist.png',dpi=250)

ax1.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.3)

ax2.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)
ax2.grid(color='grey', linestyle='-', linewidth=.3)

distM_legend=ax1.legend(handles=distM, loc='upper left', ncol=2, fancybox=True, framealpha=0., bbox_to_anchor=(0., 1.05))
distT_legend=ax2.legend(handles=distT, loc='upper left', ncol=1, fancybox=True, framealpha=0., bbox_to_anchor=(0., 1.05))

plt.savefig('mach-temp-dist-BDh.png',dpi=250, mode="expand", borderaxespad=0.)

plt.close()

