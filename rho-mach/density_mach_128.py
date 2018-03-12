import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
params = {'legend.fontsize':7.5,
          'legend.handlelength': 0.5}
plot.rcParams.update(params)

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.020,0.1,0.9,2.5))
t_start=np.array((5.0,2.0,1.0,0.5))
#t_start sets time at which statistical equilibrium has been reached
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/128/amp'+str(int(amp[i1]*1000)).rjust(4,'0')+'/'
	file=filedir+'pluto_hst.out'
	fname = open(file,'rt')
	data = np.loadtxt(file, skiprows=1)
	#Load data
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	#Ignore data before statistical equilibrium state
	#Plot the data 
	ax.plot(data[:,11],1.5*data[:,12],label=r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $A_{turb}=$'+str(amp[i1]))
	ax.plot(data[:,11],data[:,13],label=r'$\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$, $A_{turb}=$'+str(amp[i1]))
fig.set_size_inches(6, 4)
ax.set_xlabel(r'$\left< \mathcal{M}\right>_{rms}$')
ax.set_ylabel(r'$\frac{1.5\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$, $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$')
x1=np.arange(0.2,0.8,0.002)
y=np.arange(0.7,1.6,0.002)
ax.plot(x1,0.6*x1**2,label=r'$\left< \mathcal{M}\right>_{rms}^2$')
ax.plot(y,0.5*y**1,label=r'$\left< \mathcal{M}\right>_{rms}$')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$ and $\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$  vs $\left< \mathcal{M}\right>_{rms}$')
ax.set_xlim(0.1,)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
# Put a legend to the bottom of the current axis
ax.legend(loc='center left', bbox_to_anchor=(0.7, 0.4), ncol=1)
plt.savefig('rho-mach-128.png',dpi=250)
