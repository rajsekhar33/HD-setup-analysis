import numpy as np
import matplotlib.pyplot as plt

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.005,0.020,0.1,0.9,2.5))
t_start=np.array((10.0,5.0,2.0,1.0,0.5))
#t_start sets time at which statistical equilibrium has been reached
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True)
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
	ax1.plot(data[:,11],data[:,12],label='amp='+str(amp[i1]))
	ax2.plot(data[:,11],data[:,13],label='amp='+str(amp[i1]))
fig.set_size_inches(10, 6)
ax1.set_xlabel('$\mathcal{M}$')
ax1.set_ylabel(r'$\frac{\delta\rho}{\rho}$')
ax2.set_xlabel('$\mathcal{M}$')
ax2.set_ylabel(r'$\frac{\delta P}{P}$')
x1=np.arange(0.2,0.8,0.002)
y=np.arange(0.8,2.2,0.002)
ax1.plot(x1,0.35*x1**2,label='$\mathcal{M}^2$')
ax1.plot(y,0.4*y-0.1,label='$\mathcal{M}$')
x2=np.arange(0.2,2.,0.002)
ax2.plot(x2,0.5*x2**2.,label='$\mathcal{M}^2$')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax1.set_title(r'$\frac{\delta\rho}{\rho}$ vs $\mathcal{M}$')
ax2.set_title(r'$\frac{\delta P}{P}$ vs $\mathcal{M}$')
ax1.set_xlim(0.1,)
ax2.set_xlim(0.1,)
# Shrink current axis by 20%
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height])

box = ax2.get_position()
ax2.set_position([box.x0+box.x0/5., box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('rho-mach-128.png',dpi=250)
