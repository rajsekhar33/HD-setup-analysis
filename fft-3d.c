#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>
#include "params.h"
#include "binning.h"
#include "arrays.h"
#include "io.h"

void main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; //corresponds to three velocity components
  double ****velr; //4d data array to store velocity components in r-space
  Ek ****E_k; //4d data array to store Energy and corresponding k values in k-space
  Ek *E_k_added;//1d data array to store all E_k corresponding to same |k|
  Ek *E_k_binned;//1d data array to bin data
  Ek *E_k_comp;//1d data array to store compensated power spectrum
  int i,j,k;
  fftw_plan p;
  fftw_complex *out;
  double nz_r=0.5*nz+1;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz_r);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      for(k=0;k<nz_r;k++){
        if (out==NULL) printf("Allocation failure in Output fft array.\n");
      }
    }
  }
// create arrays to store various quantities
  velr = (double ****)array4d(nv,nx,ny,nz,sizeof(double));
  E_k = (Ek ****)array4d(nv,nx,ny,(nz/2+1),sizeof(Ek));
  E_k_added = (Ek *)array1d((nx*nx+ny*ny+nz*nz)/2+1,sizeof(Ek));
  E_k_binned = (Ek *)array1d(no_bins,sizeof(Ek));
  E_k_comp = (Ek *)array1d(no_bins,sizeof(Ek));
//  if(verbose) printarray4d(nv,nx,ny,nz,velr);

// create fftw plan
  double *in;
  in = &velr[0][0][0][0];
  p=fftw_plan_dft_r2c_3d(nx,ny,nz,in,out,FFTW_ESTIMATE);

  for(i=f1;i<f2;i++){
//   read data into the array
     read_dbl(i,velr);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("Reading competed in %f seconds.\n",time_taken);
//   if(verbose) printarray4d(nv,nx,ny,nz,velr);
//   do fft on each component seperately
     printf("Reading completed.\n");
     for(dir=0;dir<nv;dir++){
        in = &velr[dir][0][0][0];
        fftw_execute(p);
        printf("FFT %d completed.\n",dir);
        //fftshift(out,nz_r);
        //The following function calculates |V_k_i|^2 values
	write_E_k(dir,E_k,out);
     }
     counter(&E_k[0][0][0][0], E_k_added);//This adds up values at points having same |k|
     write_file_Ek(i,E_k_added);//This writes the above generated values into a txt file
     bin(E_k_added, E_k_binned);//This bins the data, and this function isn't entirely bug free yet. I have been first trying to get the fftshift right.
     calc_comp(E_k_added, E_k_binned);//This function calculates the compensated power spectrum
     printf("Computing Ek_%d completed.\n",dir);
     //write data to output
     write_file_Ek_binned(i,E_k_binned,E_k_comp);//This prints the final data into a txt file
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("FFT and computing V_k competed in %f seconds.\n",time_taken);
  }
  freearray4d((void ****) velr);
  freearray4d((void *) E_k);
  freearray1d((void *) E_k_added);
  freearray1d((void *) E_k_binned);
  freearray1d((void *) E_k_comp);
  fftw_destroy_plan(p);
  if(out!=NULL) fftw_free(out);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
