import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats
import pylab as plot
import pyPLUTO as pp
#plt.style.use('classic')

params = {'legend.fontsize':16,
          'legend.handlelength': 1.0}
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.size'] = 16
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['ytick.major.size'] = 14
plt.rcParams['ytick.minor.size'] = 7
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plot.rcParams.update(params)

plt.rc('text', usetex=True)

#Compute how long the simulation takes
start_time = time.time()

#Initialise the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.set_size_inches(8., 8.)

wdir=('T-runs/Tl/', 'T-runs/Th/', 'B-runs/Bl/', 'B-runs/Bh/', 'B-runs/BDh2/')
labels=('Tl', 'Th', 'Bl', 'Bh', 'BDh2')
gamma=5./3.
time_step=0.2
no_bins=200
no_files=5

start=np.array((10., 51., 10., 42., 15.))

colors=((230, 25, 75) , (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), (183, 58, 12), (240, 50, 230), (250, 190, 190), (6, 71, 24))
colors=np.array(colors)/255.

num_plots=4
ratio_spectra=[None]*start.size
spectra = [None] * (2)

for i in xrange(0, start.size):
	#Load data files
	#Declare all parameters and filenames, file location

	k1=np.zeros((no_files,no_bins))
	Ek1=np.zeros((no_files,no_bins))
	filedir='/mnt/lustre/phy/phyprtek/RAJ_RUNS/cooling_data/'+wdir[i]   
	for j in xrange(no_files):
	   fileno=str(int(start[i]/time_step+0.5+j)).rjust(4,'0')
	   filenamex1=filedir+'Ekx'+str(fileno)+'.txt'
	   filenamex2=filedir+'Eky'+str(fileno)+'.txt'
	   filenamex3=filedir+'Ekz'+str(fileno)+'.txt'
	   datax1 = np.loadtxt(filenamex1, usecols=(0,1,2))
	   datax2 = np.loadtxt(filenamex2, usecols=(0,1,2))
	   datax3 = np.loadtxt(filenamex3, usecols=(0,1,2))
	   k1[j]=datax1[:,0]
	   Ek1[j]=datax1[:,1]+datax2[:,1]+datax3[:,1]
	   all_data=pp.pload(int(start[i]/time_step+0.5+j), w_dir=filedir)
	   cs=np.sqrt(np.mean(gamma*all_data.prs/all_data.rho))
	   Ek1[j]=Ek1[j]/cs**2.
	k1=k1[0]

	k2=np.zeros((no_files,no_bins))
	Ek2=np.zeros((no_files,no_bins))

	for j in xrange(no_files):
	   fileno=str(int(start[i]/time_step+0.5+j)).rjust(4,'0')
	   filename=filedir+'Rhoks'+str(fileno)+'.txt'
	   data = np.loadtxt(filename, usecols=(0,1,2))
	   k2[j]=data[:,0]
	   Ek2[j]=data[:,1]
#Calculate average energy injection
	k2=k2[0]
	
	ratiok=Ek2/Ek1
        del_ratiok=np.std(ratiok,0)
        ratio_k=np.average(ratiok,0)


#Calculate standard deviation and mean
	del_Ek1=np.std(Ek1,0)
	Ek1=np.average(Ek1,0)
	del_Ek2=np.std(Ek2,0)
	Ek2=np.average(Ek2,0)
	j1=i
	if (i==0):
                leg1=r'$A_k=\frac{|V_k|^2}{c_s^2}$'
                leg2=r'$A_k=\frac{|\rho_k|^2}{\left<\rho\right>^2}$'
	#Plot original power spectra
	#first velocity spectra
		spectra[2*i]=ax1.errorbar(k1[1:-10],Ek1[1:-10], yerr=del_Ek1[1:-10], color=colors[j1], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=leg1, markersize=4., elinewidth=0.8)

		spectra[2*i+1]=ax1.errorbar(k1[1:-10],Ek2[1:-10], yerr=del_Ek2[1:-10], color=colors[j1], fmt='v', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=leg2, markersize=4., elinewidth=0.8)
	elif (i<num_plots):
		if(i==4):j1=i+3
		ax1.errorbar(k1[1:-10],Ek1[1:-10], yerr=del_Ek1[1:-10], color=colors[j1], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, markersize=4., elinewidth=0.8)

		ax1.errorbar(k1[1:-10],Ek2[1:-10], yerr=del_Ek2[1:-10], color=colors[j1], fmt='v', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, markersize=4., elinewidth=0.8)

	#Plot ratio
	if (i<start.size):
		if(i==4):j1=i+3
		ratio_spectra[i]=ax2.errorbar(k1[1:-10], ratio_k[1:-10], yerr=del_ratiok[1:-10], color=colors[j1], fmt='o', markeredgecolor=None, ecolor=None, capsize=None, barsabove=False, label=labels[i], markersize=4., elinewidth=1.)
	print j1
x=np.arange(10., 10**3., 1.)
fit, =ax1.plot(x, 5.*x**(-5./3.), label=r'$k^{-5/3}$', marker="d", markeredgecolor='none', markersize=0.5, linewidth=2.0, color=colors[9])
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_xlabel('k', fontsize=14)
ax1.set_ylabel('$A_k$', fontsize=20.)
#ax1.set_ylabel(r'$\frac{V_k^2}{c_s^2}$, $\frac{\rho_k^2}{\left<\rho\right>^2}$', fontsize=18)
ax1.set_ylim(1e-5,2.)
ax1.set_xlim(1.e1,1e3)
ax1.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax1.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
ax1.grid(color='grey', linestyle='-', linewidth=0.2)
spectra_legend=ax1.legend(handles=spectra, loc='upper right', bbox_to_anchor=(1.05, 1.08), ncol=2, fancybox=True, framealpha=0., fontsize=25.)
ax1.add_artist(spectra_legend)
ax1.legend(handles=[fit], loc='lower left', bbox_to_anchor=(0.3, -0.05), ncol=1, fancybox=True, framealpha=0., fontsize=25.)

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel('$k$', fontsize=20)
ax2.set_ylabel('$\eta_k^2$', fontsize=20)
ax2.set_ylim(5e-2,1.e3)
ax2.set_xlim(1.e1,1e3)
ratio_legend=ax2.legend(loc='lower right',  ncol=5, fancybox=True, bbox_to_anchor=(1.12, -0.02), frameon=True, framealpha=0.5, fontsize=20)
ax2.add_artist(ratio_legend)
fit, =ax2.plot(x, 5.e-3*x**(3./2.), label=r'$k^{3/2}$', marker="d", markeredgecolor='none', markersize=0.5, linewidth=2.0, color=colors[9])
ax2.legend(handles=[fit], loc='upper right', bbox_to_anchor=(.6, 1.0), ncol=1, fancybox=True, framealpha=0., fontsize=25.)
ax2.grid(color='grey', linestyle='-', linewidth=0.2)
ax2.tick_params(axis='both', which='major', direction='in', length=10, width=1.0, top=True, right=True)
ax2.tick_params(axis='both', which='minor', direction='in', length=5, width=0.5, top=True, right=True)
plt.savefig('ratio-rho-vel-spectra-cool.png',dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))
