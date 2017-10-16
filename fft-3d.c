#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>
#include "arrays.h"
#include "params.h"
#include "io.h"

void main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; //corresponds to three velocity components
  double ****velr; //4d data array to store velocity components in r-space
  double ****velk; //4d data array to store velocity components in k-space
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
// create a 4-d array to store velocity components
  velr = (double ****)array4d(nv,nx,ny,nz,sizeof(double));
  velk = (double ****)array4d(nv,nx,ny,(nz/2+1),sizeof(double));
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
//     if(verbose) printarray4d(nv,nx,ny,nz,velr);
//   do fft on each component seperately
     printf("Reading completed.\n");
     for(dir=0;dir<nv;dir++){
        in = &velr[dir][0][0][0];
        fftw_execute(p);
        printf("FFT %d completed.\n",dir);
	write_velk(dir,velk,out);
	printf("Computing vk_%d completed.\n",dir);
     }
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("FFT and computing V_k competed in %f seconds.\n",time_taken);

     if(verbose) printarray4d(nv,nx,ny,nz/2+1,velk);
// create energy spectrum

// write data to output
// 0: dbl
// 1: vtk
   write_output(0,f1,nv,nx,ny,nz/2+1,velk);
   write_output(1,f1,nv,nx,ny,nz/2+1,velk);
  }
  freearray4d((void ****) velr);
  freearray4d((void ****) velk);
  fftw_destroy_plan(p);
  fftw_free(out);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
