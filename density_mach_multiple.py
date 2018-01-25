import numpy as np
import matplotlib.pyplot as plt

#Declare all parameters and filenames, file location

#Load data files
amp=np.array((0.005,0.020,0.05,0.1,0.15,0.2,0.25,0.3))
t_start=np.array((10.0,5.0,2.0,2.0,2.0,2.0,2.0,2.0))
fig, ax = plt.subplots()
for i1 in xrange(0,np.shape(amp)[0]):
	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i1]*1000)).rjust(3,'0')+'/'
	file=filedir+'pluto_hst.out'
	fname = open(file,'rt')
	data = np.loadtxt(file, skiprows=1)
	#num_bins denotes the number of bins into which mach number distribution is divided into
	num_bins=500
	for i in xrange(0,data.shape[0]):
	  if data[:,0][i]>t_start[i1]:
	    break
	data=data[:][i:]
	data_sort=np.zeros((np.shape(data)[0],2))
	data_sort[:,0]=data[:,11]
	data_sort[:,1]=data[:,12]

	data_sort=np.sort(data_sort,axis=0,kind='quicksort')
	min_mach=data_sort[0][0]
	max_mach=np.max(data_sort[:,0])

	bin_size=(max_mach-min_mach)/num_bins
	#Take bins of a certain size, and add up values corresponding to thse bins
	data_binned=np.fromfunction(lambda i,j:(min_mach+i*bin_size)*(1-j),(num_bins,2))
	index=0
	for i in xrange(0,num_bins-1):
	  del_log_rho=0.
	  for j in xrange(index,np.size(data_sort[:,0])):
	     if data_sort[:,0][j]<data_binned[:,0][i+1]:
		del_log_rho+=data_sort[:,1][j]
	     else:
		data_binned[:,1][i]=del_log_rho
		index=j+1
		break

	#Plot the data 
	ax.plot(data_binned[:,0],data_binned[:,1],label='amp='+str(amp[i1]))
	fig.set_size_inches(6, 5)
	plt.xlabel('$\mathit{M}$')
	plt.ylabel(r'$\frac{\delta\rho}{\rho}$')
	leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 0.95))

	#plt.xlim(0.235,0.245)
plt.title(r'$\frac{\delta\rho}{\rho}$ vs $\mathit{M}$')
plt.savefig('rho-mach.png')

