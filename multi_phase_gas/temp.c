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
  int dir; 
  double *prs; //1d data array to store pressure
  double *rho; //1d data array to store density
  double *temp; //3d array to store temperature data at each grid point in r-space 
  int i,i1; 
// create arrays to store various quantities 
  printf("Memory allocation started\n"); 
  rho = (double *)array1d(nx*ny*nz,sizeof(double)); 
  prs = (double *)array1d(nx*ny*nz,sizeof(double)); 
  temp=  (double *)array1d(nx*ny*nz,sizeof(double)); 
  current_time=clock()-start_time; 
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
  printf("Memory allocation completed in %f seconds.\n", time_taken); 
  double *in; 
  for(i=f1;i<f2;i++){ 
//   read data into the array 
     read_dbl(i,0,rho);
     read_dbl(i,4,prs); 
     current_time=clock()-start_time; 
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Reading completed in %f seconds.\n",time_taken); 
 
    for(i1=0;i1<nx*ny*nz;i1++){ 
      temp[i1]=(prs[i1]/rho[i1])*(UNIT_VELOCITY*UNIT_VELOCITY)*(CONST_mp*CONST_mu/CONST_kB);  
//rho is density, prs is pressure
      if(temp[i1]<0.) {printf("Error here %lf, %lf %d\n",rho[i1], prs[i1], i1); exit(0);}
    } 
    temp_count(temp,i); 
    current_time=clock()-start_time; 
    time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
    printf("Temperature calculation completed in %f seconds.\n", time_taken); 
  }
  freearray1d((void *) rho);
  freearray1d((void *) prs);
  freearray1d((void *) temp);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
