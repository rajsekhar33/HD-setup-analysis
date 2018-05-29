import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
import pyPLUTO as pp

#plt.style.use('classic')
params = {'legend.fontsize':5.5,
          'legend.handlelength': 0.5}
plt.rcParams['axes.linewidth'] = .5
plt.rcParams['xtick.major.size'] = 1.
plt.rcParams['xtick.minor.size'] = .5
plt.rcParams['ytick.major.size'] = 1.5
plt.rcParams['ytick.minor.size'] = .5
plt.rcParams['xtick.labelsize'] = 6.
plt.rcParams['ytick.labelsize'] = 6.
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(3.6, 3.1)

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

spectra = [None] * (2)
fit = [None] * (2)

sb_mean=np.zeros((start.size))
for i in xrange(amp.size):
	#Load data files
        #Declare all parameters and filenames, file location

        k1=np.zeros((no_files,no_bins))
        Ek1=np.zeros((no_files,no_bins))
        Ek2=np.zeros((no_files,no_bins))
	for j1 in xrange(no_files):
		filedir='/mnt/lustre/ug4/ugrajs/fiducial_runs/256/amp'+str(int(amp[i]*10000)).rjust(5,'0')+'/run'+str(int(start[i][:,0][j1]))+'/'
		fileno=str(int(start[i][:,1][j1]/time_step)).rjust(4,'0')
		data=np.fromfile(filedir+'sbs'+fileno+'.dbl')
		sb_mean=np.average(data)
		rhok=np.loadtxt(filedir+'Rhoks'+fileno+'.txt')
		sbk =np.loadtxt(filedir+'sbks'+fileno+'.txt')/sb_mean**2.
		k1[j1]=rhok[:,0]
		Ek1[j1]=rhok[:,1]

		Ek2[j1]=sbk[:,1]
	k1=k1[0]
	k2=k1

#Calcullate ratio of power spectra of density and surface brightness
	ratiok=Ek1/Ek2/k1
        del_ratiok=np.std(ratiok,0)
        ratio_k=np.average(ratiok,0)

#Calculate standard deviation and mean
        del_Ek1=np.std(Ek1,0)
        Ek1=np.average(Ek1,0)
        del_Ek2=np.std(Ek2,0)
        Ek2=np.average(Ek2,0)

	j=i
	if (j==0): 
		leg1=r'$A_k=\frac{|\rho_k|^2}{\left<\rho\right>^2}$'
		leg2=r'$A_k=\frac{|SB_k|^2}{\left<SB\right>^2}$'
		spectra[2*j]=ax1.errorbar(k1[1:-20], Ek1[1:-20], fmt='d', yerr=del_Ek1[1:-20], color=colors[j], markeredgecolor=None, markersize=2.0, ecolor=None, capsize=None, barsabove=False, label=leg1, elinewidth=0.5)
		spectra[2*j+1]=ax1.errorbar(k2[1:-20], Ek2[1:-20], yerr=del_Ek2[1:-20], fmt='*', color=colors[j], markeredgecolor=None, markersize=2.0, ecolor=None, capsize=None, barsabove=False, label=leg2, elinewidth=0.5)
	elif (j<4):
		ax1.errorbar(k1[1:-20], Ek1[1:-20], fmt='d', yerr=del_Ek1[1:-20], color=colors[j], markeredgecolor=None, markersize=2.0, ecolor=None, capsize=None, barsabove=False, elinewidth=0.5)
		ax1.errorbar(k2[1:-20], Ek2[1:-20], yerr=del_Ek2[1:-20], fmt='*', color=colors[j], markeredgecolor=None, markersize=2.0, ecolor=None, capsize=None, barsabove=False, elinewidth=0.5)

	ax2.errorbar(k2[1:-20], ratio_k[1:-20], yerr=del_ratiok[1:-20], fmt='o', color=colors[j], markeredgecolor=None, markersize=1.5, ecolor=None, capsize=None, barsabove=False, label= '$\mathcal{M}=$'+str(mach[i]), elinewidth=0.5)


ax1.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.2)
spectra_legend=ax1.legend(handles=spectra, loc='lower left', ncol=4, fontsize=10.)
ax1.add_artist(spectra_legend)
spectra_legend.get_frame().set_alpha(0.)

x=np.arange(10., 10**3., 1.)
fit[0], = ax1.plot(x, 1e-3*x**(-5./3.), label=r'$k^{-5/3}$', linewidth=0.8)
fit[1], = ax1.plot(x, 8e-8*x**(-8./3.), label=r'$k^{-8/3}$', linewidth=0.8)
fit_legend=ax1.legend(handles=fit, loc='upper right', bbox_to_anchor=(1., 1.0), ncol=2, fontsize=6.5)
fit_legend.get_frame().set_alpha(0.)

ax1.add_artist(fit_legend)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_ylabel(r'$A_k$', fontsize=5.)
ax1.set_title(r'$A_k$ vs $k$', fontsize=6.)
ax1.set_xlim(1e1, 1e3)
ax1.set_ylim(1e-20,1e-2)
ax2.set_ylabel(r'$R_k$', fontsize=6.)

#ax2.set_ylabel(r'$\left(\frac{|\rho_k|^2}{\left<\rho\right>^2}\right)/\left(k\frac{|SB_k|^2}{\left<SB\right>^2}\right)$', fontsize=5.)
ax2.set_xlabel('$k$', fontsize=5.)
ax2.set_ylim(1e3,1e5)
ax2.tick_params(axis='both', which='major', direction='out', length=6, width=0.5, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='out', length=3, width=0.25, top=True, right=True)
ax2.grid(color='grey', linestyle='-', linewidth=0.2)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.legend(loc='upper right', bbox_to_anchor=(1., 1.0), ncol=2, fancybox=True, framealpha=0., fontsize=7.)

#plt.show()
plt.savefig('SBk_rhok.png', dpi=400)
plt.close()
