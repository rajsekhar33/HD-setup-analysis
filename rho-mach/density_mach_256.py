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
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.1,0.9,2.5))
t_start=np.array((2.,0.4,0.2))
t_start_kh=2.
perturb1 = [None] * (amp.size+1)
perturb2 = [None] * (amp.size+1)
fit = [None] * 2
#t_start sets time at which TE becomes increases by an order of magnitude
NUM_COLORS = 10
colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/fiducial_data/amp'+str(int(amp[i1]*10000)).rjust(5,'0')+'/run1/'
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1)
	#Load data
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	#Ignore data before statistical equilibrium state
	#Plot the data
	if (i1==0): 
		perturb1[i1], =ax.plot(data[:,11],5./3.*data[:,12],label=r'$\delta R=\frac{5}{3}\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $A_{turb}=$'+str(amp[i1]), color=colors[2*i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1.5)
		perturb2[i1], =ax.plot(data[:,11],data[:,13],label=r'$\delta R=\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$, $A_{turb}=$'+str(amp[i1]), color=colors[2*i1+1], marker="d", markeredgecolor='none', markersize=0.1, linewidth=1.5)
	else: 
		perturb1[i1], =ax.plot(data[:,11],5./3.*data[:,12],label=r'$A_{turb}=$'+str(amp[i1]), color=colors[2*i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1.5)
 
		perturb2[i1], =ax.plot(data[:,11],data[:,13],label=r'$A_{turb}=$'+str(amp[i1]), color=colors[2*i1+1], marker="d", markeredgecolor='none', markersize=0.1, linewidth=1.5)

fig.set_size_inches(6.5, 5.)
ax.set_xlabel(r'$\mathcal{M}_{rms}$', fontsize=18)
ax.set_ylabel(r'$\delta R$', fontsize=18)
#ax.set_ylabel(r'$\frac{5}{3}\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$', fontsize=16)
x1=np.arange(0.2,0.8,0.002)
y=np.arange(0.8,1.6,0.002)

dir_kh='/mnt/lustre/phy/phyprtek/RAJ_RUNS/fiducial_data/k12/amp00100/'
file_kh=dir_kh+'pluto_hst.out'
data_kh=np.loadtxt(file_kh, skiprows=1)
for i in xrange(0,data_kh.shape[0]):
  if data_kh[:,0][i]>t_start_kh:
    break
data_kh=data_kh[:][i:]
#Ignore data before statistical equilibrium state
perturb1[i1+1], =ax.plot(data_kh[:,11],5./3.*data_kh[:,12],label=r'$A_{turb}=$'+str(amp[0])+', $K_d=12$', color=colors[8], marker=".", markeredgecolor='none', markersize=0.1, linewidth=1.5)
perturb2[i1+1], =ax.plot(data_kh[:,11],data_kh[:,13],label=r'$A_{turb}=$'+str(amp[0])+', $K_d=12$', color=colors[7], marker="d", markeredgecolor='none', markersize=0.1, linewidth=1.5)


fit[0],= ax.plot(x1,0.6*x1**2,label=r'$\mathcal{M}_{rms}^2$', color=colors[6], linewidth=2.5)
ax.plot(y,0.6*y**2, color=colors[6], dashes=[2, 2],linewidth=2.5)
fit[1],= ax.plot(y,0.48*y**1,label=r'$\mathcal{M}_{rms}$', color=colors[9], linewidth=2.5)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(0.2,2.)
ax.set_ylim(0.03,1)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis

perturb1_leg=ax.legend(handles=perturb1, loc='upper left', bbox_to_anchor=(-0.05, 1.05), ncol=1, fontsize=20)
ax.add_artist(perturb1_leg)
perturb1_leg.get_frame().set_alpha(0.)
perturb2_leg=ax.legend(handles=perturb2, loc='lower right', bbox_to_anchor=(1.05, -0.02), ncol=1, fontsize=20)
ax.add_artist(perturb2_leg)
perturb2_leg.get_frame().set_alpha(0.)
fit_leg=ax.legend(handles=fit, loc='lower left', bbox_to_anchor=(.02, 0.3), ncol=1, fontsize=20.)
fit_leg.get_frame().set_alpha(0.)

ax.grid(color='black', linestyle='dashed', linewidth=.5, axis='x')
ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
plt.savefig('rho-mach-256.png',dpi=250)
