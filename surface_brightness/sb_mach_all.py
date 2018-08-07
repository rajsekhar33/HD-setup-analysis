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
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plot.rcParams.update(params)

plt.rc('text', usetex=True)

wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh')

step_size=0.2
end=np.array((24.8, 53.2, 24., 42.4, 55.4, 31.))
start=np.ones(end.size)*step_size
start_time=np.array((10., 51., 10., 31., 40., 10.))
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
fig.set_size_inches(6.5, 5.)
ax.set_xlabel(r'$\mathcal{M}_{rms}$', fontsize=20)
ax.set_ylabel(r'$\delta R$', fontsize=20)
#ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.3, 0.8, 0.002)
y=np.arange(0.8, 4.0, 0.002)

fit[0],= ax.plot(x1,0.25*x1**2,label=r'$\mathcal{M}_{rms}^2$', color=colors[6], linewidth=2.)
fit[1],= ax.plot(y,0.2*y**1,label=r'$\mathcal{M}_{rms}$', color=colors[9], linewidth=2.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(4e-1,3.)
ax.set_ylim(1e-2,4.)


fit_leg=ax.legend(handles=fit, loc='upper left', bbox_to_anchor=(-0.02, 1.02), ncol=2, fontsize=22)
fit_leg.get_frame().set_alpha(0.)
ax.add_artist(fit_leg)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)

for i1 in xrange(0, start.size):
#       fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	snaps=int((end[i1]-start[i1])/step_size-1)
	del_R=np.zeros((snaps))
	mach_file=filedir+'hot_mach'+str(int(start[i1]/step_size+0.1)).rjust(4,'0')+'-'+str(int(end[i1]/step_size+0.1)).rjust(4,'0')+'.txt'
	data=np.loadtxt(mach_file, skiprows=1)
	#print data.shape	
	for i in xrange(1, snaps):
		fileno=str(i).rjust(4, '0')
		sb=np.fromfile(filedir+'sbs'+fileno+'.dbl', dtype= 'double')
		del_R[i]=np.std(sb)/np.average(sb)
        #Load data
        #Plot the data
	#print del_R.shape
        perturb1[i1], =ax.plot(data[:,1][1:int(start_time[i1]/step_size)], del_R[1:int(start_time[i1]/step_size)], label=labels[i1], color=colors[i1], marker=".", markeredgecolor='none', markersize=4., linewidth=1.)
        perturb2[i1], =ax.plot(data[:,1][int(start_time[i1]/step_size):del_R.size], del_R[int(start_time[i1]/step_size):], label=labels[i1], color=colors[i1], marker=".", markeredgecolor='none', markersize=6., linewidth=1.)


perturb2_leg=ax.legend(handles=perturb2, loc='lower right', bbox_to_anchor=(1.02, -0.02), ncol=2, fontsize=20)
ax.add_artist(perturb2_leg)
perturb2_leg.get_frame().set_alpha(0.5)

plt.savefig('sb-mach-all.png',dpi=250)

