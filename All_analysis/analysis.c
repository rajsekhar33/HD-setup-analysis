#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>
#include "params.h"
#include "arrays.h"
#include "io.h"
#include "binning.h"
#include "binning2.h"
int main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; //corresponds to three velocity components
  double *store;//1d data array to store all scanned data
  Ek *E_k; //1d data array to store Spectral energy and corresponding k values in k-space
  Ek *E_k_added;//1d data array to store all E_k corresponding to same |k|
  Ek *E_k_binned;//1d data array to bin data
  Ek *E_k_comp;//1d data array to store compensated power spectrum
  double *prs; //1d data array to store pressure
  double *rho; //1d data array to store density
  double *temp; //1d array to store temperature data at each grid point in r-space 
  double *trc; //1d array to store tracer data at each grid point in r-space 
  double *sb; //1d array to store surface brightness data at each grid point in r-space 
  double *vlos; //1d array to store line of sight velocity data at each grid point in r-space 
  double *sigvlos; //1d array to store line of sight velocity dispersion data at each grid point in r-space 
  double *vx1, *vx2, *vx3;  //1d data arrays to store velocity components
  double *mach; //1d data array to store mach number data
  double cs; //speed of sound
  double no_density;//No density of particels in a particular grid cell
  double radiat_rate; //Rate of bremsstrahlung radiation
 
  int i, j, k, i1, j1, k1;
  fftw_plan p;
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz_r);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      for(k=0;k<nz_r;k++){
        if (out==NULL) printf("Allocation failure in Output fft array.\n");
      }
    }
  }
// create arrays to store various quantities
//Store stores all the real data, in the order:density, vx1, vx2, vx3 and pressure.
//The same order is followed for E_k, which stores the Fourier transformed data.
  store = (double *)array1d((nv+ntrc)*nx*ny*nz,sizeof(double));
  rho   = (double *)array1d(nx*ny*nz,sizeof(double));
  prs   = (double *)array1d(nx*ny*nz,sizeof(double));
  temp  = (double *)array1d(nx*ny*nz,sizeof(double));
  trc   = (double *)array1d(nx*ny*nz,sizeof(double));
  sb    = (double *)array1d(nx*ny,sizeof(double));
  vlos    = (double *)array1d(nx*ny,sizeof(double));
  sigvlos    = (double *)array1d(nx*ny,sizeof(double));
  vx1   = (double *)array1d(nx*ny*nz,sizeof(double));
  vx2   = (double *)array1d(nx*ny*nz,sizeof(double));
  vx3   = (double *)array1d(nx*ny*nz,sizeof(double));
  mach  = (double *)array1d(nx*ny*nz,sizeof(double));
  E_k   = (Ek *)array1d((nv+ntrc)*nx*ny*nz_r,sizeof(Ek));
  E_k_added = (Ek *)array1d((nx*nx+ny*ny+nz*nz)/2+1,sizeof(Ek));
  E_k_binned = (Ek *)array1d(no_bins,sizeof(Ek));
  E_k_comp = (Ek *)array1d(no_bins,sizeof(Ek));

  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
  printf("Memory allocation completed in %f seconds.\n", time_taken);
  double *in;

  for(i=f1;i<f2;i+=fstep){
//   read data into the array
     read_dbl(i,store);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("Reading competed in %f seconds.\n",time_taken);

//First temperature and Mach no. calculation    
     for(i1=0;i1<nx*ny*nz;i1++){
//assign value to each variable from the stored data
       rho[i1]=store[i1];
       vx1[i1]=store[i1+nx*ny*nz];
       vx2[i1]=store[i1+2*nx*ny*nz];
       vx3[i1]=store[i1+3*nx*ny*nz];
       prs[i1]=store[i1+4*nx*ny*nz];
       if(ntrc){
         trc[i1]=store[i1+5*nx*ny*nz];
         //if(fabs(trc[i1])>1.e-8 && fabs(trc[i1]-1.)>1.e-8) printf ("i1=%d, trc= %20.10e\n",i1, trc[i1]);
       }
       temp[i1]=(prs[i1]/rho[i1])*(UNIT_VELOCITY*UNIT_VELOCITY)*(CONST_mp*CONST_mu/CONST_kB);
//rho is density, prs is pressure, we bring temp to CGS units, but we calculate cs and v in code units, since mach number is dimensionless
       cs=sqrt(gamma*prs[i1]/rho[i1]);
       mach[i1]=sqrt(vx1[i1]*vx1[i1]+vx2[i1]*vx2[i1]+vx3[i1]*vx3[i1])/cs;
       if(temp[i1]<0.) {printf("Error here %lf, %lf %d\n",rho[i1], prs[i1], i1); exit(0);}
     }
     temp_count(temp,i);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Temperature calculation completed in %f seconds.\n", time_taken);

     mach_count(mach,i);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Mach no. calculation completed in %f seconds.\n", time_taken);
//Now surface brightness and vlos calculation
     for(i1=0;i1<nx;i1++){
       for(j1=0;j1<ny;j1++){
         sb[i1*ny+j1]=0.0;
         vlos[i1*ny+j1]=0.0;
         for(k1=0;k1<nz;k1++){
           no_density = rho[i1*ny*nz+j1*nz+k1]*UNIT_DENSITY/CONST_mu/CONST_mp;
           radiat_rate = no_density*no_density*lambda(temp[i1*ny*nz+j1*nz+k1]);
           sb[i1*ny+j1] += radiat_rate*(1.0/(double)nz)*UNIT_LENGTH;
//Mean velocity along line of sight
           vlos[i1*ny+j1] += vx3[i1*ny*nz+j1*nz+k1]*radiat_rate*(1.0/(double)nz)*UNIT_LENGTH*UNIT_VELOCITY;
         }
         vlos[i1*ny+j1]=vlos[i1*ny+j1]/sb[i1*ny+j1];
//Line of sight velocity dispersion calculation         
         sigvlos[i1*ny+j1]=0.0;
         for(k1=0;k1<nz;k1++){
           no_density = rho[i1*ny*nz+j1*nz+k1]*UNIT_DENSITY/CONST_mu/CONST_mp;
           radiat_rate = no_density*no_density*lambda(temp[i1*ny*nz+j1*nz+k1]);
           sigvlos[i1*ny+j1] += (vx3[i1*ny*nz+j1*nz+k1]*UNIT_VELOCITY-vlos[i1*ny+j1])*(vx3[i1*ny*nz+j1*nz+k1]*UNIT_VELOCITY-vlos[i1*ny+j1])*radiat_rate*(1.0/(double)nz)*UNIT_LENGTH;
         }
         sigvlos[i1*ny+j1]=sqrt(sigvlos[i1*ny+j1]/sb[i1*ny+j1]);
       }
     }
//Write this data into a file
     write_sb(i, sb);
     write_vlos(i, vlos, sigvlos);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Surface brightness and vlos calculation completed in %f seconds.\n", time_taken);

//Now do spectral analysis
//   do fft on each component seperately
     for(dir=0;dir<(nv+ntrc);dir++){
        in = &store[dir*nx*ny*nz];
        
// create fftw plan
        p=fftw_plan_dft_r2c_3d(nx,ny,nz,in,out,FFTW_ESTIMATE);
        fftw_execute(p);
        printf("FFT %d completed.\n",dir);
        fftshift(out);
        //The following function calculates |V_k_i|^2 values
	write_E_k(dir,&E_k[dir*nx*ny*nz_r],out);
     	counter(&E_k[dir*nx*ny*nz_r], E_k_added);//This adds up values at points having same |k|
        bin(E_k_added, E_k_binned);//This bins the data.
        calc_comp(E_k_comp, E_k_binned, dir);//This function calculates the compensated power spectrum
        printf("Computing Ek_%d completed.\n",dir);
        //write data to output
        write_file_Ek_binned(i,dir,E_k_binned,E_k_comp);//This prints the final data into a txt file
        current_time=clock()-start_time;
        time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
        printf("FFT and computing V_k competed in %f seconds.\n",time_taken);
    }
  }
  freearray1d((void *) store);
  freearray1d((void *) E_k);
  freearray1d((void *) E_k_added);
  freearray1d((void *) E_k_binned);
  freearray1d((void *) E_k_comp);
  freearray1d((void *) rho);
  freearray1d((void *) prs);
  freearray1d((void *) temp);
  freearray1d((void *) trc);
  freearray1d((void *) sb);
  freearray1d((void *) vx1);
  freearray1d((void *) vx2);
  freearray1d((void *) vx3);
  freearray1d((void *) mach);
  fftw_destroy_plan(p);
  if(out!=NULL) fftw_free(out);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
