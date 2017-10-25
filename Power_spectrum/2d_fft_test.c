#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>
#include "arrays.h"

#define CONST_PI 3.14159265358979 

void main (){
  
  srand(1000);
  double xmin=0.0;
  double xmax=1.0;
  double ymin=0.0;
  double ymax=1.0;
  int nx=300;
  int ny=300;
  double pi= CONST_PI;
  int phi_kmax=2*pi;
  double x[nx], y[ny];
  double diff_x=(xmax-xmin)/(nx-1);
  double diff_y=(ymax-ymin)/(ny-1);
  int i,j;
  for(i=0;i<nx;i++) x[i]=xmin+i*diff_x;
  for(i=0;i<ny;i++) y[i]=ymin+i*diff_y;
  int nmin_x=2;
  int nmax_x=8;
  int nmin_y=2;
  int nmax_y=8;
  double kmin_x=nmin_x*2*pi/(xmax-xmin);
  double kmax_x=nmax_x*2*pi/(xmax-xmin);
  double kmin_y=nmin_y*2*pi/(ymax-ymin);
  double kmax_y=nmax_y*2*pi/(ymax-ymin);
  double k_x,k_y,A_k,phi_k;
  int nkx, nky;
  double **f;
  f = (double **)array2d(nx,ny,sizeof(double));
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++) f[i][j]=0;
  }
  for(nkx=nmin_x;nkx<=nmax_x;nkx++){
    for(nky=nmin_y;nky<=nmax_y;nky++){
      phi_k = ((double)rand()/(double)(RAND_MAX)) * phi_kmax;
      k_x=nkx*2*pi/(xmax-xmin);
      k_y=nky*2*pi/(ymax-ymin);
      A_k=fabs(1/sqrt(k_x*k_x+k_y*k_y));
      for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
          f[i][j]+=A_k*cos(k_x*x[i]+k_y*y[i]+phi_k);
        }
      }
    }
  }
  for(nkx=nmin_x;nkx<=nmax_x;nkx++){
    for(nky=nmin_y;nky<=nmax_y;nky++){
      phi_k = ((double)rand()/(double)(RAND_MAX)) * phi_kmax;
      k_x=-nkx*2*pi/(xmax-xmin);
      k_y=nky*2*pi/(ymax-ymin);
      A_k=fabs(1/sqrt(k_x*k_x+k_y*k_y));
      for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
          f[i][j]+=A_k*cos(k_x*x[i]+k_y*y[i]+phi_k);
        }
      }
    }
  }

  double sum=0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++) sum+=fabs(f[i][j])*fabs(f[i][j])/(nx*ny);
  }
  printf("Sum %16.20lf\n",sum);
  fftw_plan p;
  fftw_complex *out;
  int ny_r=0.5*ny+1;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny_r);//Output is in form of a complex double 2D array
  for(i=0;i<nx;i++){
    for(j=0;j<ny_r;j++){
      if (out==NULL) printf("Allocation failure in Output fft array.\n");
    }
  }
  double *in;
  in=&f[0][0];//This is input double array
  p=fftw_plan_dft_r2c_2d(nx,ny,in,out,FFTW_ESTIMATE);
  fftw_execute(p);
  double sum_fft=0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny_r;j++) sum_fft+=2*cabs(out[ny_r*i+j]/(nx*ny))*cabs(out[ny_r*i+j]/(nx*ny));
  }
  printf("Sum_fft %16.20lf\n",sum_fft);
  /*double kx[(int)n_r];
  double kx_min=2*pi; 
  double kdiff=2*pi;
  for(i=0;i<n_r;i++) kx[i]=kx_min+i*kdiff;
  FILE *fp;
  char filename[100];
  strcpy(filename,"1d_test.txt");
  printf("%s\n",filename);
  fp = fopen(filename,"w");
  for(i=1;i<n_r;i++)  fprintf(fp,"%16.20lf %16.20lf %16.20lf\n",kx[i-1],2*cabs(out[i]/n)*cabs(out[i]/n),0.5*pow(kx[i-1],-2));
  fclose(fp);
*/

}
