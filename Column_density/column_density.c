#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<time.h>
#include "params.h"
#include "arrays.h"
#include "io.h"

void main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; //corresponds to three velocity components
  double ***rho; //3d data array to store density
  double **rho_col; //2d data array to store column density
  int i,j,k;
// create arrays to store various quantities
  printf("Memory allocation started.\n");
  rho = (double ***)array3d(nx,ny,nz,sizeof(double));
  rho_col = (double **)array2d(nx,ny,sizeof(double));
//  if(verbose) printarray4d(nv,nx,ny,nz,velr);
  printf("Memory allocation completed.\n");
  for(i=f1;i<f2;i++){
//   read data into the array
     read_dbl(i,rho);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("Reading competed in %f seconds.\n",time_taken);
     for(i=0;i<nx;i++){
       for(j=0;j<ny;j++){
         for(k=0;k<nz;k++){
           if(k==0) rho_col[i][j]=0;
           rho_col[i][j]+=rho[i][j][k];
         }
       }
     }
     write_output(i,nx,ny,rho_col);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  }
  freearray4d((void ****) rho);
  freearray2d((void **) rho_col);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
