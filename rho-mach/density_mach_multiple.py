import numpy as np
import matplotlib.pyplot as plt

#Declare all parameters and filenames, file location

#Load data files
#amp=np.array((0.005,0.010,0.015,0.020,0.025,0.035,0.05))#,0.1,0.15,0.2,0.25,0.3))
#t_start=np.array((10.0,7.0,6.0,5.0,5.0,5.0,2.0))#,2.0,2.0,2.0,2.0,2.0))
#amp=np.array((0.35,0.6,0.9,1.5,2.5))
#t_start=np.array((1.5,1.0,0.8,0.5,0.5))
amp=np.array((0.005,0.010,0.015,0.020,0.025,0.035,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.6,0.9,1.5,2.5))
t_start=np.array((10.0,7.0,6.0,5.0,5.0,5.0,2.0,2.0,2.0,2.0,2.0,2.0,1.5,1.0,0.8,0.5,0.5))
fig, ax = plt.subplots(1)
fig.set_size_inches(6, 5)
for i1 in xrange(0,np.shape(amp)[0]):
#	fig, ax = plt.subplots(1)
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i1]*1000)).rjust(4,'0')+'/'
	file=filedir+'pluto_hst.out'
	fname = open(file,'rt')
	data = np.loadtxt(file, skiprows=1)
	#num_bins denotes the number of bins into which mach number distribution is divided into
	num_bins=50
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	#Plot the data 
	ax.plot(data[:,11],data[:,12],label='amp='+str(amp[i1]))
	fig.set_size_inches(6, 5)
	plt.xlabel('$\mathit{M}$')
	plt.ylabel(r'$\frac{\delta\rho}{\rho}$')
"""	plt.title(r'$\frac{\delta\rho}{\rho}$ vs $\mathit{M}$ for amp='+str(amp[i1]))
	plt.savefig('rho-mach-amp'+str(int(amp[i1]*1000)).rjust(3,'0')+'.png')
	plt.savefig(filedir+'rho-mach-amp'+str(int(amp[i1]*1000)).rjust(3,'0')+'.png')
	plt.clf()
"""	
x=np.arange(0.08,0.8,0.002)
#x=np.arange(0.4,0.8,0.002)
y=np.arange(0.8,2.2,0.002)
#ax.plot(x,0.35*x**2,label='$\mathit{M}^2$')
ax.plot(y,0.4*y-0.1,label='$\mathit{M}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.title(r'$\frac{\delta\rho}{\rho}$ vs $\mathit{M}$')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.savefig('rho-small-mach.png',dpi=1000)
#plt.savefig('rho-large-mach.png',dpi=1000)
plt.savefig('rho-mach.png',dpi=250)

