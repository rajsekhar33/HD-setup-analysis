import numpy as np
import matplotlib.pyplot as plt

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.005,0.020,0.1,0.9,2.5))
t_start=np.array((10.0,5.0,2.0,1.0,0.5))
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
	ax.plot(data[:,11],1.5*data[:,12],label=r'$1.5\frac{\delta\rho}{\rho}$, amp='+str(amp[i1]))
	ax.plot(data[:,11],data[:,13],label=r'$\frac{\delta P}{P}$,amp='+str(amp[i1]))
fig.set_size_inches(10, 6)
ax.set_xlabel('$\mathcal{M}$')
ax.set_ylabel(r'$\frac{\delta\rho}{\rho}$')
x1=np.arange(0.2,0.8,0.002)
y=np.arange(0.7,1.6,0.002)
ax.plot(x1,0.6*x1**2,label='$\mathcal{M}^2$')
ax.plot(y,0.5*y**1,label='$\mathcal{M}$')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title(r'$\frac{\delta\rho}{\rho}$ and $\frac{\delta P}{P}$  vs $\mathcal{M}$')
ax.set_xlim(0.1,)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('rho-mach-128.png',dpi=250)
