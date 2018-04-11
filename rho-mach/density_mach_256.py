import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
import seaborn as sns

#plt.style.use('classic')
params = {'legend.fontsize':10,
          'legend.handlelength': 2}
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.020,0.1,0.9,2.5))
t_start=np.array((5.0,2.0,1.0,0.5))
#t_start sets time at which statistical equilibrium has been reached

NUM_COLORS = 10
colors=((230, 25, 75), (250, 190, 190) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (210, 245, 60), (145, 30, 180), (0, 128, 128), (240, 50, 230))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i1]*1000)).rjust(4,'0')+'/'
	file=filedir+'pluto_hst.out'
	data = np.loadtxt(file, skiprows=1)
	#Load data
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	#Ignore data before statistical equilibrium state
	#Plot the data 
	line1 =ax.plot(data[:,11],1.5*data[:,12],label=r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $A_{turb}=$'+str(amp[i1]))
        line1[0].set_linewidth(.5)
	line1[0].set_color(colors[i1])
	line2=ax.plot(data[:,11],data[:,13],label=r'$\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$, $A_{turb}=$'+str(amp[i1]))
	line1[0].set_color(colors[i1+5])
        line2[0].set_linewidth(.5)
fig.set_size_inches(6, 4.5)
ax.set_xlabel(r'$\left< \mathcal{M}\right>_{rms}$')
ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$')
x1=np.arange(0.3,0.8,0.002)
y=np.arange(0.8,1.6,0.002)
ax.plot(x1,0.6*x1**2,label=r'$\left< \mathcal{M}\right>_{rms}^2$')
ax.plot(y,0.5*y**1,label=r'$\left< \mathcal{M}\right>_{rms}$')
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(0.3,2.)
ax.set_ylim(0.02,1)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.), ncol=2)
plt.savefig('rho-mach-256.png',dpi=250)

