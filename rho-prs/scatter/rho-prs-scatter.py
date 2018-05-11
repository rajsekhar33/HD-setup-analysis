import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import scipy.optimize as so
import pylab as plot
import pyPLUTO as pp
from matplotlib.colors import LogNorm

params = {'legend.fontsize':7.5,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plot.rcParams.update(params)

#Compute how long the simulation takes
start_time = time.time()


start=np.array((100, 75, 20, 14, 5, 2))
amp=np.array((0.005, 0.02, 0.1, 0.1, 0.9,2.5))
mach=((0.24, 0.43, 0.76, 0.90, 1.25, 2.1))
epsilon=((0.001255, 0.0088, 0.142, 0.133, 6.16, 25.45))
cs=((0.59, 0.69, 0.87, 0.75, 1.982, 1.97))
time_step=0.2
no_bins=200
no_files=2

colors=((230, 25, 75), (250, 190, 190) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (210, 245, 60), (145, 30, 180), (0, 128, 128), (240, 50, 230))
colors=np.array(colors)/255.

for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/'
	data = pp.pload(start[i], w_dir=filedir)
	Rho=np.ndarray.flatten(data.rho)
	Prs=np.ndarray.flatten(data.prs)
	H=np.zeros((no_bins,no_bins))
	H, xedges, yedges = np.histogram2d(Prs, Rho, bins=(no_bins,no_bins), normed=True)
	x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1, no_bins))
	y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((no_bins, 1))

	pdf = (H*(x_bin_sizes*y_bin_sizes))

	X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
	XX, YY = np.meshgrid(X,Y)
	Z = pdf#.T

	#Initialise the figure
	fig, ax = plt.subplots()
	fig.set_size_inches(7, 5)
	cax=ax.contourf(XX, YY, Z, norm=LogNorm(vmin=1e-8, vmax=.1))
	ax.set_ylabel(r'$\frac{\left<\delta\rho\right>_{rms}}{\left<\rho\right>}$')
	ax.set_xlabel(r'$\frac{\left<\delta P\right>_{rms}}{\left< P\right>}$')

	ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
	ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
	if(i==0): fig.colorbar(cax)
	plt.savefig('rho-prs-pdf-mach-'+str(int(mach[i]*1000)).rjust(3, '0')+'.png', dpi=200)
	plt.clf()
