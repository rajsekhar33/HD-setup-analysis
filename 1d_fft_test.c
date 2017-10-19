#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>

#define CONST_PI 3.14159265358979 

void main (){
  
  srand(1000);
  double xmin=0.0;
  double xmax=1.0;
  int n=1000;
  int phi_kmax=100;
  double pi= CONST_PI;
  double x[n];
  double diff=(xmax-xmin)/(n-1);
  int i;
  for(i=0;i<n;i++) x[i]=xmin+i*diff;
  int nmin=2;
  int nmax=8;
  double kmin=nmin*2*pi/(xmax-xmin);
  double kmax=nmax*2*pi/(xmax-xmin);
  double k,A_k,phi_k;
  int nk;
  double y[n];
  for(i=0;i<n;i++) y[i]=0;
  for(i=0;i<n;i++){
    for(nk=nmin;nk<=nmax;nk++){
      k=nk*2*pi/(xmax-xmin);
      A_k=fabs(1/k);
      phi_k = ((double)rand()/(double)(RAND_MAX)) * phi_kmax;
      y[i]+=A_k*cos(k*x[i]+phi_k);
    }
  }
  double sum=0;
  for(i=0;i<n;i++) sum+=fabs(y[i])*fabs(y[i])/n;
  printf("Sum %16.20lf\n",sum);
  fftw_plan p;
  fftw_complex *out;
  double n_r=0.5*n+1;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n_r);//Output is in form of a complex double 1D array
  for(i=0;i<n_r;i++){
    if (out==NULL) printf("Allocation failure in Output fft array.\n");
  }
  double *in;
  in=&y[0];//This is input double array
  p=fftw_plan_dft_r2c_1d(n,in,out,FFTW_ESTIMATE);
  fftw_execute(p);
  double sum_fft=0;
  for(i=0;i<n_r;i++) sum_fft+=2*cabs(out[i]/n)*cabs(out[i]/n);
  printf("Sum_fft %16.20lf\n",sum_fft);
  double kx[(int)n_r];
  double kx_min=2*pi; 
  double kdiff=2*pi;
  for(i=0;i<n_r;i++) kx[i]=kx_min+i*kdiff;
  FILE *fp;
  char filename[100];
  strcpy(filename,"1d_test.txt");
  printf("%s\n",filename);
  fp = fopen(filename,"w");
  for(i=1;i<n_r;i++)  fprintf(fp,"%16.20lf %16.20lf\n",kx[i-1],2*cabs(out[i]/n)*cabs(out[i]/n));
  fclose(fp);


}
