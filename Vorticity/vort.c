//This code calculates the vorticity of a velocity field
//Input and output are in form of 4d data array

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "arrays.h"
#include "params.h"
#include "vort_calculator.h"
#include "io.h"

void main()
{
  int dir; //corresponds to three velocity components
  double ****velr; //4d data array to store velocity components in r-space
  double ****vort; //4d data array to store vorticity components in r-space
  int i,j,k;
  
// create a 4-d array to store velocity and vorticity components
  velr = (double ****)array4d(nv,nx,ny,nz,sizeof(double));
  vort = (double ****)array4d(nv,nx,ny,nz,sizeof(double));
//  if(verbose) printarray4d(nv,nx,ny,nz,velr);

// verbose is a parameter, if set to a non-zero value, prints the data to the screen
  double *in;
  in = &velr[0][0][0][0];

  for(i=f1;i<f2;i++){
//   read data into the array
     read_dbl(i,velr);
//   Find all vorticity components
     vort_calculator(in,vort);  
     printf("Vorticity calculation completed.\n");
  if(verbose){
  	printarray4d(nv,nx,ny,nz,vort);
	}
// write data to output
// 0: dbl
// 1: vtk
   write_output(0,i,nv,nx,ny,nz,vort);
//   write_output(1,f1,nv,nx,ny,nz/2+1,vort);
  }
  
  freearray4d((void ****) velr);
  freearray4d((void ****) vort);
  printf("Everything completed\n");
//  print_output(0,f1,nv,nx,ny,nz);
  exit(0);
}
