import numpy as np
import matplotlib.pyplot as plt
import pylab as plot

#plt.style.use('classic')
params = {'legend.fontsize':16,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 16
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.major.size'] = 14
plt.rcParams['ytick.minor.size'] = 7
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Declare all parameters and filenames, file location

#wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'T-runs/TDh/', 'B-runs/BDh/')
#labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')
wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'B-runs/BDh2/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'BDh2')

step_size=0.2
end=np.array((24.8, 53.2, 22.2, 53.2, 33.6))
start=np.ones(end.size)*step_size
#start=np.array((10., 51., 10., 42., 10.4, 42., 15.))
#Load data files
perturb1 = [None] * (start.size)
perturb2 = [None] * (start.size)
fit = [None] * 2
#t_start sets time at which statistical equilibrium has been reached

NUM_COLORS = 10
colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
for i1 in xrange(0, start.size):
#	fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i1]   
	file=filedir+'hot_mach'+str(int(start[i1]/step_size+0.1)).rjust(4,'0')+'-'+str(int(end[i1]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1)
	#Load data
	#Ignore data before statistical equilibrium state
	#Plot the data
	perturb1[i1], =ax.plot(data[:,1],data[:,0], label=labels[i1], color=colors[i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1)

fig.set_size_inches(8, 6.5)
ax.set_xlabel(r'$\left< \mathcal{M}\right>_{rms}$', fontsize=18)
ax.set_ylabel(r'$\delta R$', fontsize=18)
#ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.3,0.8,0.002)
y=np.arange(0.8,1.6,0.002)


#fit[0],= ax.plot(x1,0.6*x1**2,label=r'$\left< \mathcal{M}\right>_{rms}^2$', color=colors[9], linewidth=2.)
#fit[1],= ax.plot(y,0.48*y**1,label=r'$\left< \mathcal{M}\right>_{rms}$', color=colors[6], linewidth=2.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')

ax.set_xlim(0.2,5.)
#ax.set_ylim(0.5,4)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis

perturb1_leg=ax.legend(handles=perturb1, loc='upper left', bbox_to_anchor=(-0.05, 1.0), ncol=1, fontsize=25)
ax.add_artist(perturb1_leg)
perturb1_leg.get_frame().set_alpha(0.)
#fit_leg=ax.legend(handles=fit, loc='lower left', bbox_to_anchor=(0., 0.0), ncol=2, fontsize=25)
#fit_leg.get_frame().set_alpha(0.)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)
plt.savefig('rho-mach-cool.png',dpi=250)

