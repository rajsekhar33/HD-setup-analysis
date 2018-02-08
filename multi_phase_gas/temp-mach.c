#include<math.h> 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<time.h>
#include "params.h"
#include "arrays.h"
#include "io_temp.h"
//#include "binning.h"
#include "binning2.h"

void main()
{
  clock_t start_time, current_time;
  double time_taken;
  start_time=clock();
  int dir; 
  double *store; //1d array to store all scanned data
  double *prs; //1d data array to store pressure
  double *rho; //1d data array to store density
  double *temp; //1d array to store temperature data at each grid point in r-space 
  double *vx1, *vx2, *vx3;  //1d data arrays to store velocity components
  double *mach; //1d data array to store mach number data
  double cs; //speed of sound
  int i,i1; 
// create arrays to store various quantities 
  printf("Memory allocation started\n");
  store=(double *)array1d(5*nx*ny*nz,sizeof(double)); 
  rho = (double *)array1d(nx*ny*nz,sizeof(double)); 
  prs = (double *)array1d(nx*ny*nz,sizeof(double)); 
  temp=  (double *)array1d(nx*ny*nz,sizeof(double));
  vx1=  (double *)array1d(nx*ny*nz,sizeof(double));
  vx2=  (double *)array1d(nx*ny*nz,sizeof(double));
  vx3=  (double *)array1d(nx*ny*nz,sizeof(double));
  mach=  (double *)array1d(nx*ny*nz,sizeof(double));
   
  current_time=clock()-start_time; 
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
  printf("Memory allocation completed in %f seconds.\n", time_taken); 
  for(i=f1;i<f2;i++){ 
//   read data into the store array
     read_dbl(i, 5, store);//We are scanning 5 different variables into the storage array 
     current_time=clock()-start_time; 
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
     printf("Reading completed in %f seconds.\n",time_taken); 
 
    for(i1=0;i1<nx*ny*nz;i1++){
//assign value to each variable from the stored data
      rho[i1]=store[i1];
      vx1[i1]=store[i1+nx*ny*nz];
      vx2[i1]=store[i1+2*nx*ny*nz];
      vx3[i1]=store[i1+3*nx*ny*nz];
      prs[i1]=store[i1+4*nx*ny*nz];
 
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
 
  }
  freearray1d((void *) store);
  freearray1d((void *) rho);
  freearray1d((void *) prs);
  freearray1d((void *) temp);
  freearray1d((void *) vx1);
  freearray1d((void *) vx2);
  freearray1d((void *) vx3);
  freearray1d((void *) mach);
  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
  printf("Everything completed in %f seconds.\n",time_taken);
  exit(0);
}
