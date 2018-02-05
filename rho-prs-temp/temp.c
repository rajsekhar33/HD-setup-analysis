#include<math.h> 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<time.h>
#include "params.h"
#include "arrays.h"
#include "io_temp.h"
#include "binning.h"

void main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; //corresponds to three velocity components
  double ****velr; //4d data array to store velocity components in r-space
  double *temp; //3d array to store temperature data at each grid point in r-space 
  Ek ****E_k; //4d data array to store Energy and corresponding k values in k-space 
  Ek *E_k_added;//1d data array to store all E_k corresponding to same |k| 
  Ek *E_k_binned;//1d data array to bin data 
  Ek *E_k_comp;//1d data array to store compensated power spectrum 
  int i,j,k,i1,j1,k1; 
  double nz_r=0.5*nz+1; 
// create arrays to store various quantities 
  printf("Memory allocation started\n"); 
  velr = (double ****)array4d(nv,nx,ny,nz,sizeof(double)); 
  temp=  (double *)array1d(nx*ny*nz,sizeof(double)); 
  current_time=clock()-start_time; 
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
  printf("Memory allocation completed in %f seconds.\n", time_taken); 
  double *in; 
  for(i=f1;i<f2;i++){ 
//   read data into the array 
     read_dbl(i,velr); 
     current_time=clock()-start_time; 
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Reading completed in %f seconds.\n",time_taken); 
 
    for(i1=0;i1<nx;i1++){ 
      for(j1=0;j1<ny;j1++){ 
        for(k1=0;k1<nz;k1++){ 
          temp[ny*nz*i1+nz*j1+k1]=(velr[1][i1][j1][k1]/velr[0][i1][j1][k1])*(UNIT_VELOCITY*UNIT_VELOCITY)*(CONST_mp*CONST_mu/CONST_kB);  
//velr[0] is density, velr[1] is pressure 
        } 
      } 
    } 
    temp_count(&temp[0],i); 
    current_time=clock()-start_time; 
    time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
    printf("Temperature calculation completed in %f seconds.\n", time_taken); 
  }
  freearray4d((void ****) velr);
  freearray1d((void ****) temp);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
