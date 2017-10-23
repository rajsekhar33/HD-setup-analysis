import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift

np.random.seed(1001)
xmin = 0.0; xmax = 1.0; n = 10000 #choose to ve even
x = np.linspace(xmin,xmax,n); y = np.zeros(n)
nmin = 2; nmax = 8
kmin = nmin*2*math.pi/(xmax-xmin); kmax = nmax*2*math.pi/(xmax-xmin)

for nk in range(nmin, nmax+1):
    k = nk*2*math.pi/(xmax-xmin)
    A_k = np.abs(1.0/k)
    phi_k = np.random.uniform(0.0,2.*math.pi)
    y =  y + A_k*np.cos(k*x+phi_k)
    
print ( sum(abs(y)*abs(y))/n )

#plt.plot(x,y,'+-')
#plt.grid()
#plt.show()

fft_y = fft(y)
print ( sum(abs(fft_y)*abs(fft_y))/n/n, sum(2*abs(fft_y[1:n/2+1])*abs(fft_y[1:n/2+1])/n/n) )
#ya = ifft(fft_y)
#xax = np.linspace(0,n-1,n)
#plt.loglog(xax, np.abs(fft_y)/n,'-+',xax,xax**-1.)

kx = 2*math.pi*np.linspace(1,n/2,n/2)

plt.loglog(kx, 2*abs(fft_y[1:n/2+1])*abs(fft_y[1:n/2+1])/n/n, '-+', kx, 0.5*kx**-2)

plt.ylabel(r'|FFT|$^2$')
plt.xlabel('k')

plt.grid()
plt.show()