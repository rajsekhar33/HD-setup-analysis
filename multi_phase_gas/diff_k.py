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
start=40
end=50
time_step=0.2
num_bins=1000

temp1=np.zeros((end-start,num_bins*0.9))
count1=np.zeros((end-start,num_bins*0.9))
filedir1='/mnt/lustre/ug4/ugrajs/cooling/256/'
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir1+'temp'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1))
   temp1[filenumber-start]=data[:,0][:num_bins*0.9]
   count1[filenumber-start]=data[:,1][:num_bins*0.9]
temp1=temp1[0]
del_count1=np.std(count1,0)
#count1=np.average(count1,0)


temp2=np.zeros((end-start,num_bins*0.9))
count2=np.zeros((end-start,num_bins*0.9))
filedir2='/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k4-6/'
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir2+'temp'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1))
   temp2[filenumber-start]=data[:,0][:num_bins*0.9]
   count2[filenumber-start]=data[:,1][:num_bins*0.9]
temp2=temp2[0]
del_count2=np.std(count2,0)
#count2=np.average(count2,0)


temp3=np.zeros((end-start,num_bins*0.9))
count3=np.zeros((end-start,num_bins*0.9))
filedir3='/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k6-8/'
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir3+'temp'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1))
   temp3[filenumber-start]=data[:,0][:num_bins*0.9]
   count3[filenumber-start]=data[:,1][:num_bins*0.9]
temp3=temp3[0]
del_count3=np.std(count3,0)
#count3=np.average(count3,0)



temp4=np.zeros((end-start,num_bins*0.9))
count4=np.zeros((end-start,num_bins*0.9))
filedir4='/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k8-10/'
for filenumber in xrange(start,end):
   fileno=str(filenumber).rjust(4,'0')
   filename=filedir4+'temp'+str(fileno)+'.txt'
   fname = open(filename,'rt')
   data = np.loadtxt(filename, usecols=(0,1))
   temp4[filenumber-start]=data[:,0][:num_bins*0.9]
   count4[filenumber-start]=data[:,1][:num_bins*0.9]
temp4=temp4[0]
del_count4=np.std(count4,0)
#count4=np.average(count4,0)




#Plot the data 


#This is to plot the Temperature distribution function 

for i in xrange(0,end-start):
	time1=(start+i)*(time_step)
	fig, ax = plt.subplots()
	fig.set_size_inches(6, 5)
	ax.plot(temp1, count1[i],label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
	ax.plot(temp2, count2[i],label=r'$4 \leq |k_{driving}| \leq 6$')
	ax.plot(temp3, count3[i],label=r'$6 \leq |k_{driving}| \leq 8$')
	ax.plot(temp4, count4[i],label=r'$8 \leq |k_{driving}| \leq 10$')
	#ax.errorbar(temp1, count1, yerr=del_count1,label=r'$0 < |k_{driving}| \leq \sqrt{2}$')
	#ax.errorbar(temp2, count2, yerr=del_count2,label=r'$4 \leq |k_{driving}| \leq 6$')
	#ax.errorbar(temp3, count3, yerr=del_count3,label=r'$6 \leq |k_{driving}| \leq 8$')
	#ax.errorbar(temp4, count4, yerr=del_count4,label=r'$8 \leq |k_{driving}| \leq 10$')
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlabel('Temperature (K)')
	plt.ylabel('Number density')
	plt.title('Temperature distribution of gas, different $k_{driving}$' )
	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig('Temp_dist_256_diff_k_t-'+str(time1)+'.png',dpi=250)
	print("--- %s seconds ---" % (time.time() - start_time))
