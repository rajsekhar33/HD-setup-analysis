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
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plot.rcParams.update(params)

plt.rc('text', usetex=True)

wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'no_turb/2e-1/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')

step_size=0.2
end=np.array((24.8, 53.2, 24., 42.4, 19., 55.4, 31.))
start=np.ones(end.size)*step_size
start_time=np.array((10., 51., 10., 31., 6., 40., 10.))
#Load data files
perturb = [None] * (start.size)
#t_start sets time at which statistical equilibrium has been reached


#Load data files
perturb1 = [None] * (start.size)
perturb2 = [None] * (start.size)
fit = [None] * 2
#t_start sets time at which statistical equilibrium has been reached

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
fig.set_size_inches(8, 6.5)
ax.set_xlabel(r'$\mathcal{M}_{rms}$', fontsize=20)
ax.set_ylabel(r'$\delta R$', fontsize=20)
#ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.3,0.8,0.002)
y=np.arange(0.8,1.6,0.002)

fit[0],= ax.plot(x1,0.6*x1**2,label=r'$\mathcal{M}_{rms}^2$', color=colors[9], linewidth=2.)
fit[1],= ax.plot(y,0.48*y**1,label=r'$\mathcal{M}_{rms}$', color=colors[6], linewidth=2.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(4e-2,2.)
ax.set_ylim(0.03,1.)


fit_leg=ax.legend(handles=fit, loc='upper left', bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=25)
fit_leg.get_frame().set_alpha(0.)
ax.add_artist(fit_leg)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)

for i1 in xrange(0, start.size):
#       fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
        file=filedir+'hot_mach'+str(int(start[i1]/step_size+0.1)).rjust(4,'0')+'-'+str(int(end[i1]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file, skiprows=1)
        #Load data
        #Plot the data
        perturb1[i1], =ax.plot(data[:,1][:int(start_time[i1]/step_size)],data[:,0][:int(start_time[i1]/step_size)], label=labels[i1], color=colors[i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1.)
        perturb2[i1], =ax.plot(data[:,1][int(start_time[i1]/step_size):],data[:,0][int(start_time[i1]/step_size):], label=labels[i1], color=colors[i1], marker=".", markeredgecolor='none', markersize=0.2, linewidth=2.5)


perturb2_leg=ax.legend(handles=perturb2, loc='lower left', bbox_to_anchor=(0., 0.), ncol=2, fontsize=25)
ax.add_artist(perturb2_leg)
perturb2_leg.get_frame().set_alpha(0.)

plt.savefig('rho-mach-all.png',dpi=250)

