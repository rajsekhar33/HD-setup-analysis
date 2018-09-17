import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

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
wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'T-runs/TDh/', 'B-runs/BDh/', 'Q-run/')

labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh', 'QD')

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.

fig, ax = plt.subplots()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
axins = plt.axes([0,0,1,1])
ip = InsetPosition(ax, [0.1,0.1,0.5,0.3])
axins.set_axes_locator(ip)
fig.set_size_inches(6, 5)
for i1 in xrange(0,no_files):
        filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i1]
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1, usecols=(0,13))
	ax.plot(data[:,0]*UNIT_TIME,data[:,1],label=labels[i1], color=colors[i1])
	axins.plot(data[:,0]*UNIT_TIME,data[:,1],label=labels[i1], color=colors[i1])
ax.set_xlabel(r'time (Myr)', fontsize=18)
ax.set_ylabel(r'$\left<\delta\rho\right>_{rms}/\left<\rho\right>$', fontsize=20)
#ax.set_ylabel(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$', fontsize=24)
ax.set_xlim(0,2200)
ax.set_ylim(0.,3.3)

ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
ax.grid(color='grey', linestyle='-', linewidth=0.2)

axins.set_xlim(0.,80.)
axins.set_ylim(0.,0.6)

plt.xticks(visible=True)
plt.yticks(visible=True)
mark_inset(ax, axins, loc1=2, loc2=4, fc="gray", ec="0.2")
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ vs time')
# Shrink current axis by 20%
# Put a legend to the bottom of the current axis
ax.legend(loc='upper right', ncol=4, bbox_to_anchor=(1.02, 1.02), fancybox=True, framealpha=0.5, fontsize=17.)
plt.savefig('del-rho-time.png',dpi=250)
