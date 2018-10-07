import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
from scipy.optimize import curve_fit
import pyPLUTO as pp

#plt.style.use('classic')
params = {'legend.fontsize':16,
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

UNIT_VELOCITY= (1.e8)
CONST_pc=3.0856775807e18
CONST_mp=1.67262171e-24
UNIT_LENGTH =  (CONST_pc*40.e3)
UNIT_TIME=UNIT_LENGTH/UNIT_VELOCITY/(3.15e13)
CONST_mu=0.5
CONST_kB=1.3806505e-16
#wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'no_turb/2e-1/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'B-runs/BDh2/')
#labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')
#labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh')
labels=('Tl', 'Th', 'Bl', 'Bh', 'BDh2')

step_size=0.2
end=np.array((24.8, 53.2, 22., 42.4, 31.))
start_time=0.2*np.ones(end.shape)
num_snap=(end-start_time)/0.2
#start_time sets time at which statistical equilibrium has been reached

temp_cutoff=1.225e6
leg_plot = [None] * (start_time.size)


colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
fig.set_size_inches(8., 6.)
ax.set_xlabel(r'$t$ $(Myr)$', fontsize=20)
ax.set_ylabel(r'Fraction of cold gas mass', fontsize=20)

ax.grid(color='black', linestyle='dashed', linewidth=.5)
ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)

for i1 in xrange(0, start_time.size):
#       fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i1]
	dist_data=np.zeros((int(num_snap[i1]), 2))
	for j in xrange(int(num_snap[i1])):
		snap=int(start_time[i1]/step_size+0.1)+j
		d1=pp.pload(snap, w_dir=filedir)
		temp=(d1.prs/d1.rho)*(UNIT_VELOCITY*UNIT_VELOCITY)*(CONST_mp*CONST_mu/CONST_kB)
		dist_data[j][0]=snap*0.2*UNIT_TIME
		dist_data[j][1]=np.sum(d1.rho[(temp<temp_cutoff)])/256.**3.
		print np.sum(d1.rho[(temp<temp_cutoff)])/256.**3.
	j1=i1
	if(i1==4):j1=i1+3
	leg_plot[i1],=ax.plot(dist_data[:,0], dist_data[:,1], color=colors[j1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=2., label=labels[i1])
        #Plot the data
legend1=ax.legend(handles=leg_plot[::2], loc='upper right', bbox_to_anchor=(1.0, 1.04), ncol=1, fontsize=18.)
ax.add_artist(legend1)
legend1.get_frame().set_alpha(0.)

legend2=ax.legend(handles=leg_plot[1::2], loc='upper left', bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=18.)
ax.add_artist(legend2)
legend2.get_frame().set_alpha(0.)
#ax.legend()

ax.set_xlim(0., 2e3)
#ax.set_ylim(0., 3e-3)
plt.savefig('cold-gas-time.png',dpi=250)

