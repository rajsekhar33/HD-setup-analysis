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

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.1,0.9,2.5,0.1))
start=np.array((2.0,1.0,0.8,0.6))
end=np.array((50., 50., 50., 12.))
step_size=0.2
perturb = [None] * (amp.size)
fit = [None] * 2
#t_start sets time at which statistical equilibrium has been reached

NUM_COLORS = 10
colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i1]*10000)).rjust(5,'0')+'/run2/'
	if (i1==3):
		filedir='/mnt/lustre/ug4/ugrajs/higher_k/256/k12/amp00100/'
	snaps=int((end[i1]-start[i1])/step_size-1)
	mach_rms=np.zeros((snaps))
	del_R=np.zeros((snaps))
	file1=filedir+'pluto_hst.out'
	data = np.loadtxt(file1, skiprows=1)
	for i in xrange(1, snaps):
                fileno=str(i).rjust(4, '0')
                sb=np.fromfile(filedir+'sbs'+fileno+'.dbl', dtype= 'double')
                del_R[i]=np.std(sb)/np.average(sb)
                mach_rms[i]=np.average(data[(data[:,0]>step_size*i)*(data[:,0]<step_size*i+1e-2)][:,11])
	
	#Ignore data before statistical equilibrium state
	#Plot the data
	if (i1==0): 
		perturb[i1], =ax.plot(mach_rms, del_R,label=r'$\delta R=\frac{\left<\delta SB \right>_{rms}}{\left< SB \right>}$, $A_{turb}=$'+str(amp[i1]), color=colors[2*i1+1], marker=".", markeredgecolor='none', markersize=1., linewidth=1.)
	elif(i1>0 and i1<3): 
		perturb[i1], =ax.plot(mach_rms, del_R, label=r'$A_{turb}=$'+str(amp[i1]), color=colors[2*i1+1], marker=".", markeredgecolor='none', markersize=1., linewidth=1.)
 	else:
		perturb[i1], =ax.plot(mach_rms, del_R, label=r'$A_{turb}=$'+str(amp[i1])+', $k_d=12$', color=colors[7], marker=".", markeredgecolor='none', markersize=1., linewidth=1.)

fig.set_size_inches(8, 6.5)
ax.set_xlabel(r'$\mathcal{M}_{rms}$', fontsize=18)
ax.set_ylabel(r'$\delta R$', fontsize=18)
#ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.1,0.8,0.002)
y=np.arange(0.8,4.0,0.002)

fit[0],= ax.plot(x1,0.25*x1**2,label=r'$\mathcal{M}_{rms}^2$', color=colors[6], linewidth=3.)
fit[1],= ax.plot(y,0.2*y**1,label=r'$\mathcal{M}_{rms}$', color=colors[9], linewidth=3.)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(0.1,4.)
ax.set_ylim(2e-3,1)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis

perturb_leg=ax.legend(handles=perturb, loc='upper left', bbox_to_anchor=(-0.05, 1.05), ncol=1, fontsize=22)
ax.add_artist(perturb_leg)
perturb_leg.get_frame().set_alpha(0.)
fit_leg=ax.legend(handles=fit, loc='lower right', bbox_to_anchor=(1.0, 0.0), ncol=2, fontsize=22)
fit_leg.get_frame().set_alpha(0.)

ax.grid(color='black', linestyle='dashed', linewidth=1.5, axis='x')
ax.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)
plt.savefig('sb-mach-fiducial.png',dpi=250)

