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
wdir=('thermal_heating/256/tabulated_cooling/F5e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F1e-1/k0-2/', 'thermal_heating/256/tabulated_cooling/F5e-2/k0-2/', 'thermal_heating/256/tabulated_cooling/F3e-2/k0-2/')
#labels=('Tl', 'Th', 'Bl', 'Bh', 'QD', 'TDh', 'BDh')
labels=('0.5', '0.1', '0.05', '0.03')

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
end=np.array((24., 55.8, 60.6, 41.6))
start_time=np.array((6.4, 8.6, 14.6, 11.2))
num_snap=(end-start_time)/0.2
#start_time sets time at which statistical equilibrium has been reached

sb_vlos = [None] * (start_time.size+1)
fit_label = [None] * 2


colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.
fig, ax = plt.subplots()
fig.set_size_inches(8., 6.)
ax.set_xlabel(r'$v_{LOS}$ $(km/s)$', fontsize=20)
ax.set_ylabel(r'$SB$ (normalised)', fontsize=20)

ax.grid(color='black', linestyle='dashed', linewidth=.5)
ax.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)

#define the gaussian gunction with all parameters
def f(x, mu, sigma, height): 
	return height * 1/(np.sqrt(2.*np.pi)*sigma)*np.exp(-0.5*((x-mu)/sigma)**2.)
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
	print popt, perr
	# prepare confidence level curves
	nstd = 1. # to draw 1-sigma intervals
	popt_up = popt+nstd*perr
	popt_dw = popt-nstd*perr

	fit = f(x1, *popt)
	height=popt[2]
        #Plot the data
        sb_vlos[i1],=ax.plot(x1, dist_data/height, color=colors[i1], marker=".", markeredgecolor='none', markersize=0.1, linewidth=2., label='$f=$'+labels[i1]+', $\sigma_v=$'+str(int(popt[1]))+'km/s')
	ax.fill_between(x1, (dist_data+del_dist_data)/height, (dist_data-del_dist_data)/height, color=colors[i1], alpha=.25)
        #ax.errorbar(x1, dist_data/height, yerr=del_dist_data/height, color=colors[i1], fmt=".", markeredgecolor='none', markersize=.1, linewidth=1.)
sb_vlos[i1+1],=ax.plot(x, y, label='Hitomi'+' $\sigma_v=$'+str(int(sig))+'km/s', color=colors[9], marker=".", markeredgecolor='none', markersize=0.1, linewidth=2.)

ax.fill_between(x, y1, y2, color=colors[9], alpha=.25)
sb_legend1=ax.legend(handles=sb_vlos[::2], loc='upper left', bbox_to_anchor=(0.0, 1.02), ncol=1, fontsize=18.)
ax.add_artist(sb_legend1)
sb_legend1.get_frame().set_alpha(0.)

sb_legend2=ax.legend(handles=sb_vlos[1::2], loc='upper right', bbox_to_anchor=(1.0, 1.02), ncol=1, fontsize=18.)
ax.add_artist(sb_legend2)
sb_legend2.get_frame().set_alpha(0.)
#ax.legend()

ax.set_xlim(-1e3, 1e3)
ax.set_ylim(0., 3e-3)
plt.savefig('sb-vlos-pdf-Bl.png',dpi=250)
ax.set_yscale('log')
ax.set_ylim(1e-5, 1e-2)
plt.savefig('sb-vlos-pdf-log-Bl.png',dpi=250)

