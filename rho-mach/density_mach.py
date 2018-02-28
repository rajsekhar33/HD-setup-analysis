import numpy as np
import matplotlib.pyplot as plt

#Declare all parameters and filenames, file location

#Load data files
amp=0.01
filedir1='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp*1000)).rjust(3,'0')+'/'
file=filedir1+'pluto_hst.out'
fname = open(file,'rt')
data1 = np.loadtxt(file, skiprows=1)
#num_bins denotes the number of bins into which mach number distribution is divided into
num_bins=500
t_start1=7.0

for i in xrange(0,data1.shape[0]):
  if data1[:,0][i]>t_start1:
    break
data1=data1[:][i:]
data_sort1=np.zeros((np.shape(data1)[0],2))
data_sort1[:,0]=data1[:,11]
data_sort1[:,1]=data1[:,12]

data_sort1=np.sort(data_sort1,axis=0,kind='quicksort')
min_mach1=data_sort1[0][0]
max_mach1=np.max(data_sort1[:,0])

bin_size1=(max_mach1-min_mach1)/num_bins
#Take bins of a certain size, and add up values corresponding to thse bins
data_binned1=np.fromfunction(lambda i,j:(min_mach1+i*bin_size1)*(1-j),(num_bins,2))
index=0
for i in xrange(0,num_bins-1):
  del_log_rho=0.
  for j in xrange(index,np.size(data_sort1[:,0])):
     if data_sort1[:,0][j]<data_binned1[:,0][i+1]:
	del_log_rho+=data_sort1[:,1][j]
     else:
	data_binned1[:,1][i]=del_log_rho
	index=j+1
	break

#Plot the data 

plt.figure()
plt.plot(data_binned1[:,0],data_binned1[:,1])
plt.xlabel('$\mathcal{M}$')
plt.ylabel(r'$\frac{\delta\rho}{\rho}$')
#plt.xlim(0.235,0.245)
plt.title(r'$\frac{\delta\rho}{\rho}$ vs $\mathcal{M}$ for amp='+str(amp))
plt.savefig('rho-mach-amp'+str(int(amp*1000)).rjust(3,'0')+'.png')
plt.savefig(filedir1+'rho-mach-amp'+str(int(amp*1000)).rjust(3,'0')+'.png')

