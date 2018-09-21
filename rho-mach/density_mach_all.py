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
plt.rcParams['xtick.labelsize'] = 17.
plt.rcParams['ytick.labelsize'] = 17.
plot.rcParams.update(params)

plt.rc('text', usetex=True)

wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'B-runs/BDh2/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'BDh2')

step_size=0.2
#end=np.array((25., 53.2, 22.2, 53.2))#, 55.4, 31.))
end=np.array((25., 53.2, 22.2, 53.2, 33.6))
start=np.ones(end.size)*step_size
start_time=np.array((10., 51., 10., 31., 10.))
#Load data files
perturb = [None] * (start.size)
#t_start sets time at which statistical equilibrium has been reached


#Load data files
perturb1 = [None] * (start.size)
perturb2 = [None] * 2
fit = [None] * 2
#t_start sets time at which statistical equilibrium has been reached

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
fig.set_size_inches(6.5, 5.)
ax.set_xlabel(r'$\mathcal{M}_{rms}$', fontsize=20)
ax.set_ylabel(r'$\delta R$', fontsize=20)
#ax.set_ylabel(r'$\frac{5}{3}\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.25,0.8,0.002)
y=np.arange(0.8,3.,0.002)

fit[0],= ax.plot(x1,0.6*x1**2,label=r'$\mathcal{M}_{rms}^2$', color=colors[6], linewidth=2.)
ax.plot(y,0.6*y**2, dashes=[2, 2],color=colors[6], linewidth=2.)
fit[1],= ax.plot(y,0.48*y**1,label=r'$\mathcal{M}_{rms}$', color=colors[9], linewidth=2.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(25e-2,3.)
ax.set_ylim(0.03,2.)


fit_leg=ax.legend(handles=fit, loc='upper left', bbox_to_anchor=(-0.05, 1.05), ncol=2, fontsize=22)
fit_leg.get_frame().set_alpha(0.)
ax.add_artist(fit_leg)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)

for i1 in xrange(0, start.size):
#       fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i1]
        file=filedir+'hot_mach'+str(int(start[i1]/step_size+0.1)).rjust(4,'0')+'-'+str(int(end[i1]/step_size+0.1)).rjust(4,'0')+'.txt'
	data = np.loadtxt(file)
        #Load data
	print file 
        #Plot the data
        ax.plot(data[:,2][:int(start_time[i1]/step_size)],5./3.*data[:,1][:int(start_time[i1]/step_size)], label=labels[i1], color=colors[i1], marker=".", markersize=2., linewidth=.5)
        perturb2[0], =ax.plot(data[:,2][int(start_time[i1]/step_size):],5./3.*data[:,1][int(start_time[i1]/step_size):], label=r'$\frac{5}{3}\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$', color=colors[i1], marker=".", markersize=3., linewidth=0.5)

        perturb1[i1],=ax.plot(data[:,2][:int(start_time[i1]/step_size)],data[:,0][:int(start_time[i1]/step_size)], label=labels[i1], color=colors[i1], marker="v", markersize=2., linewidth=.5)
        perturb2[1], =ax.plot(data[:,2][int(start_time[i1]/step_size):],data[:,0][int(start_time[i1]/step_size):], label=r'$\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', color=colors[i1], marker="v", markersize=3., linewidth=0.5)


perturb1_leg=ax.legend(handles=perturb1, loc='lower right', bbox_to_anchor=(1.02, -0.02), ncol=2, fontsize=22)
ax.add_artist(perturb1_leg)
perturb1_leg.get_frame().set_alpha(0.5)
perturb2_leg=ax.legend(handles=perturb2, loc='lower right', bbox_to_anchor=(1.02, 0.3), ncol=1, fontsize=22)
ax.add_artist(perturb2_leg)
perturb2_leg.get_frame().set_alpha(0.5)

plt.savefig('rho-prs-mach-all.png',dpi=250)

