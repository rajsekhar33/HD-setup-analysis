#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<fftw3.h>
#include<time.h>
#include "params-spherical.h"
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
  Ek *sbk; //1d data array to store FFT data of surface brightness
  Ek *sbk_added;//1d data array to store all sbk corresponding to same |k|
  Ek *sbk_binned;//1d data array to bin data
  double **hot_emission;//2d data array to bin data
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

  FILE *hot_file; //file to store hot phase gas data 
  char filenumb1[5];
  char filenumb2[5];
  char filename[100];

  //Variables for naming hot gas data file 
  sprintf(filenumb1,"%04d",f1);
  sprintf(filenumb2,"%04d",f2);

  strcpy(filename,datdir);
  strcat(filename, "hot_mach");
  strcat(filename,filenumb1);
  strcat(filename,"-");
  strcat(filename,filenumb2);
  strcat(filename,".txt");
  printf("%s\n",filename);
  hot_file = fopen(filename,"w");

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
  sbk   = (Ek *)array1d(nx*ny_r,sizeof(Ek));
  sbk_added = (Ek *)array1d((nx*nx+ny*ny)/2+1,sizeof(Ek));
  sbk_binned = (Ek *)array1d(no_bins,sizeof(Ek));
//Although hot emission stores velocity and emission, they are both doubles 
  hot_emission = (double **)array2d(no_hot_bins, 2, sizeof(Ek));

  current_time=clock()-start_time;
  time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds 
  printf("Memory allocation completed in %f seconds.\n", time_taken);
  double *in;
  double r;//r stores the distance of the grid point from the centre
  int i2, j2, k2;
  double del_rho; //del_rho stores the density perturbation at a grid point
  double del_prs; //del_prs stores the pressure perturbation at a grid point
  double hot_rho0, hot_delrho, hot_prs0, hot_delprs, hot_delv, hot_cs, hot_count;
  for(i=f1;i<f2;i+=fstep){
//Initialise values to zero
     hot_rho0=0; hot_delrho=0; hot_prs0=0; hot_delprs=0; hot_delv=0; hot_cs=0; hot_count=0;

     for(i1=0;i1<no_hot_bins;i1++){
       hot_emission[i1][0]=vel_min+(double)(i1*del_v);
       hot_emission[i1][1]=0.;
     }
//   read data into the array
     read_dbl(i,store);
     current_time=clock()-start_time;
     time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
     printf("Reading competed in %f seconds.\n",time_taken);

//First temperature and Mach no. calculation    
     for(i1=0;i1<nx*ny*nz;i1++){
//assign value to each variable from the stored data
       k2=i1%nz;
       j2=(i1-k2)%(ny*nz)/nz;
       i2=(i1-k2-j2*nz)/ny/nz;
       r=sqrt((i2-nx/2.)*(i2-nx/2.)+(j2-ny/2.)*(j2-ny/2.)+(k2-nz/2.)*(k2-nz/2.))/(double)nx;
       del_rho=store[i1]-rho0;
       rho[i1]=rho0+del_rho/2.*(1-tanh((r-r0)/sigma));
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
//Calculations separate for hot gas
       if(temp[i1]>1e5){
         hot_count++; hot_cs+=cs; 
         hot_rho0+=rho[i1]; hot_delrho+=rho[i1]*rho[i1];
         hot_prs0+=prs[i1]; hot_delprs+=prs[i1]*prs[i1];
         hot_delv+=(vx1[i1])*(vx1[i1])+(vx2[i1])*(vx2[i1])+(vx3[i1])*(vx3[i1]);
         no_density = rho[i1]*UNIT_DENSITY/CONST_mu/CONST_mp;
         radiat_rate = no_density*no_density*lambda(temp[i1])*(UNIT_LENGTH/nx)*(UNIT_LENGTH/ny)*(UNIT_LENGTH/nz);
         //hot_bin(hot_emission, vx3[i1], radiat_rate);
         if(vx3[i1]<vel_max && vx3[i1]>vel_min){
           hot_emission[find_my_i(vx3[i1], hot_emission, 0, no_hot_bins)][1]+=radiat_rate;
         }
       }
       if(temp[i1]<0.) {printf("Error here %lf, %lf %d\n",rho[i1], prs[i1], i1); exit(0);}
     }
     hot_cs/=hot_count; 
     hot_rho0/=hot_count; hot_delrho=sqrt(hot_delrho/hot_count-hot_rho0*hot_rho0);
     hot_prs0/=hot_count; hot_delprs=sqrt(hot_delprs/hot_count-hot_prs0*hot_prs0);
     hot_delv=sqrt(hot_delv/hot_count);

     //Write to file
     write_file_hot_radiat_binned(i, hot_emission);
     fprintf(hot_file, "%16.20lf %16.20lf %16.20lf\n", hot_delprs/hot_prs0, hot_delrho/hot_rho0, hot_delv/hot_cs);

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
           if(temp[i1*ny*nz+j1*nz+k1]>5e6){
             no_density = rho[i1*ny*nz+j1*nz+k1]*UNIT_DENSITY/CONST_mu/CONST_mp;
             radiat_rate = no_density*no_density*lambda(temp[i1*ny*nz+j1*nz+k1]);
             sb[i1*ny+j1] += radiat_rate*(1.0/(double)nz)*UNIT_LENGTH;
//Mean velocity along line of sight
           vlos[i1*ny+j1] += vx3[i1*ny*nz+j1*nz+k1]*radiat_rate*(1.0/(double)nz)*UNIT_LENGTH*UNIT_VELOCITY;
           }
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
	write_E_k(&E_k[dir*nx*ny*nz_r],out);
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
//do fft for surface brightness
//
    p=fftw_plan_dft_r2c_2d(nx,ny,sb,out,FFTW_ESTIMATE);
    fftw_execute(p);
    printf("FFT of sb completed.\n");
    fftshift2(out);
//The following function calculates |V_k_i|^2 values
    write_sbk(sbk,out);
    counter2(sbk, sbk_added);//This adds up values at points having same |k|
//    bin2(sbk_added, sbk_binned);//This bins the data.
    bin(sbk_added, sbk_binned);//This bins the data. Follow same binning as 3D spectra, this helps in plotting ratios
    printf("Computing SBk completed.\n");
//write data to output
    write_file_sbk_binned(i,sbk_binned);//This prints the final data into a txt file
    current_time=clock()-start_time;
    time_taken=((double)current_time)/CLOCKS_PER_SEC; // in seconds
    printf("FFT and computing SB_k competed in %f seconds.\n",time_taken);
  }
  fclose(hot_file);
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
  freearray1d((void *) sbk);
  freearray1d((void *) sbk_added);
  freearray1d((void *) sbk_binned);
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
