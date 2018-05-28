import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
#plt.style.use('classic')
params = {'legend.fontsize':7.5,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

#Initialise the figure
fig, ax = plt.subplots()
fig.set_size_inches(7, 5)

amp=np.array((0.005, 0.02, 0.1, 0.1, 0.9,2.5))
mach=((0.25, 0.45, 0.75, 0.90, 1.2, 2.1))
cs=((0.59, 0.69, 0.87, 0.75, 1.982, 1.97))
time_step=0.2
no_bins=200
no_files=5

start_array1=np.array(((1,27.0), (2,30.8), (3,20.4), (4,25.2), (5,17.0), (5,25.6)))
start_array2=np.array(((1,10.8), (1,12.8), (1,17.6), (1,20.4), (2,12.8)))
start_array3=np.array(((1,4.4), (1,5.4), (2,5.4), (3,5.6), (4,5.2)))
start_array4=np.array(((1,2.8), (2,3.2), (3,3.0), (4,4.2), (5,3.4)))
start_array5=np.array(((1,1.0), (2,1.4), (3,1.0), (4,1.2), (5,1.6)))
start_array6=np.array(((1,0.4), (2,0.4), (3,0.4), (4,0.4), (5,0.4)))

start=np.array((start_array1, start_array2, start_array3, start_array4, start_array5, start_array6))

colors=((230, 25, 75), (250, 190, 190) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (210, 245, 60), (145, 30, 180), (0, 128, 128), (240, 50, 230))
colors=np.array(colors)/255.

num_plots=3
spectra = [None] * (num_plots*2)

for i in xrange(0, 1):
	#Load data files
	#Declare all parameters and filenames, file location

	Ek1=np.zeros((no_files,no_bins))
	Ekcomp1=np.zeros((no_files,no_bins))
	
	for j in xrange(no_files):
	   filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run'+str(int(start[i][:,0][j]))+'/'
	   file=filedir+'pluto_hst.out'
	   data1 = np.loadtxt(file, skiprows=1, usecols=(0,10,14))
	   fileno=str(int(start[i][:,1][j]/time_step)).rjust(4,'0')
	   filenamex1=filedir+'Ekx'+str(fileno)+'.txt'
	   filenamex2=filedir+'Eky'+str(fileno)+'.txt'
	   filenamex3=filedir+'Ekz'+str(fileno)+'.txt'
	   datax1 = np.loadtxt(filenamex1, usecols=(0,1,2))
	   datax2 = np.loadtxt(filenamex2, usecols=(0,1,2))
	   datax3 = np.loadtxt(filenamex3, usecols=(0,1,2))
	   k1=datax1[:,0]
	   Ek1[j]=datax1[:,1]+datax2[:,1]+datax3[:,1]
	   epsilon= np.average(data1[(data1[:,0]>start[i][:,1][j])*(data1[:,0]<start[i][:,1][j]+1e-1)][:,2])
	   print epsilon
	   Ekcomp1[j]=epsilon**(-2/3)*k1**(5./3.)*(datax1[:,1]+datax2[:,1]+datax3[:,1])

	for j in xrange(no_files-1):
	   ax.errorbar(k1[1:-10], (Ekcomp1[j]/Ekcomp1[j+1])[1:-10], color=colors[j], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label='$t=$'+str(start[i][:,1][j])+', run'+str(int(start[i][:,0][j])), markersize=1.5, elinewidth=0.4)


ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('k')
ax.set_ylabel(r'Ratio of compensated spectra')
#ax.set_ylim(1.e-10,1.)
ax.set_xlim(1.e1,1e3)
ax.legend(loc='lower right',  ncol=3, fancybox=True, framealpha=0.)
#ax.set_title(r'Ratio of velocity and density power spectra' )
ax.grid(color='grey', linestyle='-', linewidth=0.2)
ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
plt.savefig('ratio-vel-spectra.png',dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))
