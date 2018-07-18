import numpy as np
import matplotlib.pyplot as plt
import pylab as plot
from scipy.optimize import curve_fit

#plt.style.use('classic')
params = {'legend.fontsize':16,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 16
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.major.size'] = 14
plt.rcParams['ytick.minor.size'] = 7
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'no_turb/2e-1/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
wdir=('tabulated_cooling/256/k0-2/', 'tabulated_cooling/256/k12/', 'thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-1/k12/', 'turb_perturb/DkHC2e-1/', 'turb_perturb/DBh2e-1/F5e-1/')
#labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')
labels=('Tl', 'Th', 'Bl', 'Bh', 'TDh', 'BDh')

no_bins=200
UNIT_VELOCITY = 1.e8/1e5
#"""
sig=165.
sig1=sig+10.
sig2=sig-10.
"""
sig=150.
sig1=sig+70.
sig2=sig-70.
"""
step_size=0.2
end=np.array((24.8, 53.2, 24., 42.4, 55.4, 31.))
start_time=np.array((10., 51., 10., 31., 40., 10.))
num_snap=(end-start_time)/0.2
#start_time sets time at which statistical equilibrium has been reached

sb_vlos = [None] * (start_time.size+1)
fit_label = [None] * 2


colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
fig.set_size_inches(9.3, 6.2)
ax.set_xlabel(r'$V_{los}$ $(km/s)$', fontsize=20)
ax.set_ylabel(r'$SB$ (normalised)', fontsize=20)

ax.grid(color='black', linestyle='dashed', linewidth=.5)
ax.tick_params(axis='both', which='major', direction='out', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='out', length=5, width=0.5, top=True, right=True)

#define the gaussian gunction with all parameters
def f(x, mu, sigma, height): 
	return height * 1/(np.sqrt(2.*np.pi)*sigma)*np.exp(-((x-mu)/sigma)**2.)
x=np.arange(-1000., 1000., 0.1)


y=1/(np.sqrt(2.*np.pi)*sig)*np.exp(-x**2./(2.*sig**2.))
y1=1/(np.sqrt(2.*np.pi)*sig1)*np.exp(-x**2./(2.*sig**2.))
y2=1/(np.sqrt(2.*np.pi)*sig2)*np.exp(-x**2./(2.*sig**2.))

for i1 in xrange(0, start_time.size):
#       fig, ax = plt.subplots(1)
        filedir='/mnt/lustre/ug4/ugrajs/cooling/'+wdir[i1]
	dist_data=np.zeros((int(num_snap[i1]), no_bins))
	for j in xrange(int(num_snap[i1])):
		snap=int(start_time[i1]/step_size+0.1)+j
		file=filedir+'radiat_vel'+str(snap).rjust(4,'0')+'.txt'
		#Load data
		data = np.loadtxt(file)
		data[:,1]=data[:,1]/(np.sum(data[:,1])*no_bins*1.e1)
		dist_data[j]=data[:,1]
	del_dist_data=np.std(dist_data,0)
	dist_data=np.average(dist_data,0)
	x1=UNIT_VELOCITY*data[:,0]
	# curve fit [with only y-error]
	popt, pcov = curve_fit(f, x1, dist_data, sigma=del_dist_data, bounds=([-100., 0., 0.], [100., 1000., 5e2]))
	#popt, pcov = curve_fit(f, x1, dist_data, sigma=del_dist_data)
	perr = np.sqrt(np.diag(pcov))

	# prepare confidence level curves
	nstd = 1. # to draw 1-sigma intervals
	popt_up = popt+nstd*perr
	popt_dw = popt-nstd*perr

	fit = f(x1, *popt)
	height=popt[2]
	sb_vlos[i1],=ax.plot(x1, fit/height, color=colors[i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=2., label=labels[i1]+' $\sigma_v=$'+str(int(popt[1]))+'km/s')
        #Plot the data
        ax.fill_between(x1, (dist_data+0.5*del_dist_data)/height, (dist_data-0.5*del_dist_data)/height, color=colors[i1], alpha=.25)
        #ax.errorbar(x1, dist_data/height, yerr=del_dist_data/height, color=colors[i1], fmt=".", markeredgecolor='none', markersize=.1, linewidth=1.)
sb_vlos[i1+1],=ax.plot(x, y, label='Hitomi'+' $\sigma_v=$'+str(int(sig))+'km/s', color=colors[9], marker=".", markeredgecolor='none', markersize=0.1, linewidth=2.)
ax.fill_between(x, y1, y2, color=colors[9], alpha=.25)
sb_legend1=ax.legend(handles=sb_vlos[::2], loc='upper right', bbox_to_anchor=(1.05, 1.0), ncol=1, fontsize=24.)
ax.add_artist(sb_legend1)
sb_legend1.get_frame().set_alpha(0.)

sb_legend2=ax.legend(handles=sb_vlos[1::2], loc='upper left', bbox_to_anchor=(0.0, 1.0), ncol=1, fontsize=24.)
ax.add_artist(sb_legend2)
sb_legend2.get_frame().set_alpha(0.)
#ax.legend()

ax.set_xlim(-1e3, 1e3)
ax.set_ylim(0., 2.8e-3)
plt.savefig('sb-vlos-pdf.png',dpi=250)
