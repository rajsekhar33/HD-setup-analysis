import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
from scipy.optimize import curve_fit
from scipy import stats

#Compute how long the simulation takes
start_time = time.time()
#Declare all parameters and filenames, file location

#Load data files
file='/home/rajsekhar/Final_year_project/HD_Module/Analysis/diff_forcing/pluto_hst.out'
n=np.array([32,32,32])
fname = open(file,'rt')

data = np.loadtxt(file, skiprows=1, usecols=(0,1,2,3,4,5,6,7,8,9,10))

t=data[:,0]
dt=data[:,1]
mass=data[:,2]
TE=data[:,3]
KE1=data[:,4]
KE2=data[:,5]
KE3=data[:,6]
MOM1=data[:,7]
MOM2=data[:,8]
MOM3=data[:,9]
epsilon=data[:,10]
TOTE=TE+KE1+KE2+KE3
TOTKE=KE1+KE2+KE3

rate1=(TOTE[10:]-TOTE[:-10])/(t[10:]-t[:-10])
rate2=(TOTKE[10:]-TOTKE[:-10])/(t[10:]-t[:-10])

#Plot the data 

#Here we plot the compensated power spectrum, multiplying E(k) with k^(5/3)
fig, ax = plt.subplots()
ax.plot(t,epsilon,'-',label='$\epsilon$')
ax.plot(t[0:-10],rate1,'-',label=r'$\frac{d(Tot E)}{dt}$')
plt.xlabel('t')
plt.ylabel('Energy dissipation rate, $\epsilon$')
leg = ax.legend(loc=2, bbox_to_anchor=(0.75, 1.15))
plt.title('Comparing $\epsilon$ and Energy dissipation rates ' )

plt.savefig('epsilon.png')

print("--- %s seconds ---" % (time.time() - start_time))
