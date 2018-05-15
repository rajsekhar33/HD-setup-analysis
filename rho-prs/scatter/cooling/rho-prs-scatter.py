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


start=np.array((10., 51., 10., 42., 10.4, 42., 15.))
no_bins=200
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'no_turb/2e-1/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')
time_step=0.2

for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i]
	
	data = pp.pload(int(start[i]/time_step), w_dir=filedir)
	Rho_mean=np.mean(data.rho)
	Prs_mean=np.mean(data.prs)
	Rhon=np.log(np.ndarray.flatten(data.rho)/Rho_mean)
	#print np.average(Rhon)
	Prsn=np.log(np.ndarray.flatten(data.prs)/Prs_mean)
	H=np.zeros((no_bins,no_bins))
	H, xedges, yedges = np.histogram2d(Prsn, Rhon, bins=(no_bins,no_bins), normed=True)
	x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1, no_bins))
	y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((no_bins, 1))

	pdf = (H*(x_bin_sizes*y_bin_sizes))

	X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
	XX, YY = np.meshgrid(X,Y)
	Z = pdf#.T

	#Initialise the figure
	fig, ax = plt.subplots()
	fig.set_size_inches(7, 5)
	cax=ax.pcolormesh(XX, YY, Z, norm=LogNorm(vmin=1e-6, vmax=.1), cmap='plasma')
	ax.set_ylabel(r'$log\left(\frac{\rho}{\left<\rho\right>}\right)$')
	ax.set_xlabel(r'$log\left(\frac{P}{\left< P\right>}\right)$')
	#ax.set_xlim(-25., 5.)
	#ax.set_ylim(-2., 4.)

	ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
	ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
	fig.colorbar(cax)
	#plt.savefig('rho-prs-pdf-scaled-'+labels[i]+'.png', dpi=200)
	plt.savefig('rho-prs-pdf-unscaled-'+labels[i]+'.png', dpi=200)
	plt.clf()
	plt.close()
