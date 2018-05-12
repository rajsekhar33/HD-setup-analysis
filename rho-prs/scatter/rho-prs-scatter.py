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
no_bins=200

for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/'
	data = pp.pload(start[i], w_dir=filedir)
	Rho_mean=np.mean(data.rho)
	Prs_mean=np.mean(data.prs)
	Rhon=np.log(np.ndarray.flatten(data.rho)/Rho_mean)
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
	cax=ax.contourf(XX, YY, Z, norm=LogNorm(vmin=1e-8, vmax=.1), )
	ax.set_ylabel(r'$log\left(\frac{\rho}{\left<\rho\right>}\right)$')
	ax.set_xlabel(r'$log\left(\frac{P}{\left< P\right>}\right)$')
	#ax.set_xlim(-4., 2.)
	#ax.set_ylim(-2., 2.)

	ax.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
	ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
	fig.colorbar(cax)
	#plt.savefig('rho-prs-pdf-scaled-mach-'+str(int(mach[i]*100)).rjust(3, '0')+'.png', dpi=200)
	plt.savefig('rho-prs-pdf-unscaled-mach-'+str(int(mach[i]*100)).rjust(3, '0')+'.png', dpi=200)
	#plt.savefig('rho-prs-pdf-scaled-mach-'+str(int(mach[i]*100)).rjust(3, '0')+'.png', dpi=200)
	plt.clf()
