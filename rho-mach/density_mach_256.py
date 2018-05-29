import numpy as np
import matplotlib.pyplot as plt
import pylab as plot

#plt.style.use('classic')
params = {'legend.fontsize':14,
          'legend.handlelength': 2}
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3


plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.1,0.9,2.5))
t_start=np.array((2.0,1.0,0.5))
perturb = [None] * (amp.size*2)
fit = [None] * 2
#t_start sets time at which statistical equilibrium has been reached

NUM_COLORS = 10
colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i1]*10000)).rjust(5,'0')+'/run1/'
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1)
	#Load data
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	#Ignore data before statistical equilibrium state
	#Plot the data 
	perturb[i1], =ax.plot(data[:,11],1.5*data[:,12],label=r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $A_{turb}=$'+str(amp[i1]), color=colors[2*i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1)
	perturb[i1+3], =ax.plot(data[:,11],data[:,13],label=r'$\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$, $A_{turb}=$'+str(amp[i1]), color=colors[2*i1+1], marker="d", markeredgecolor='none', markersize=0.1, linewidth=1.0)

fig.set_size_inches(8, 6)
ax.set_xlabel(r'$\left< \mathcal{M}\right>_{rms}$', fontsize=16)
ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.3,0.8,0.002)
y=np.arange(0.8,1.6,0.002)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

fit[0],= ax.plot(x1,0.6*x1**2,label=r'$\left< \mathcal{M}\right>_{rms}^2$', color=colors[9], linewidth=2.)
fit[1],= ax.plot(y,0.48*y**1,label=r'$\left< \mathcal{M}\right>_{rms}$', color=colors[6], linewidth=2.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(0.3,2.)
ax.set_ylim(0.03,1)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis

perturb_leg=ax.legend(handles=perturb, loc='lower right', bbox_to_anchor=(1., 0.0), ncol=2, fontsize=16)
ax.add_artist(perturb_leg)
perturb_leg.get_frame().set_alpha(0.)
fit_leg=ax.legend(handles=fit, loc='upper left', bbox_to_anchor=(0., 1.0), ncol=2, fontsize=16)
fit_leg.get_frame().set_alpha(0.)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.75, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.5, top=True, right=True)
plt.savefig('rho-mach-256.png',dpi=250)

