void read_dbl(int f, double *arr)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];
   int i;
   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,prefix);
   strcat(filename,filenumb);
   strcat(filename,suffix);

   printf("%s\n",filename);
   fp = fopen(filename,"rb");
   if(fp==NULL){
     printf("fopen error (1)");
     exit(0);
   }
   for(i=0;i<(nv+ntrc)*nx*ny*nz;i++){
     fread(&arr[i],sizeof(double),1,fp);
   }
   fclose(fp);
   return;
}
//We use this function to shift the complex output data before we calculate the |Vk_i|^2 values
void fftshift(fftw_complex *out){
double complex temp;
int nz_r=nz/2+1;
int i,j,k,i_shifted,j_shifted;
  if (out==NULL) printf("Allocation failure in Output fft array.\n");
  for(i=0;i<nx;i++){
    for(j=0;j<ny/2;j++){
      for(k=0;k<nz_r;k++){
        i_shifted=(i+nx/2)%nx;
        j_shifted=(j+ny/2)%ny;
        temp=out[i*ny*nz_r+j*nz_r+k];
        out[i*ny*nz_r+j*nz_r+k]=out[i_shifted*ny*nz_r+j_shifted*nz_r+k];
        out[i_shifted*ny*nz_r+j_shifted*nz_r+k]=temp;
      }
    }
  }
}
//This function calculate Ek^2, from the FFT output data
void write_E_k(int dir, Ek *E_k, fftw_complex *out)
{
   int i,j,k;
   int nz_r=0.5*nz+1;
   double real, imag;
   double i_sq, j_sq;
   for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
         for(k=0;k<nz_r;k++){
            if (out==NULL) printf("Allocation failure in Output fft array.\n");
	    real = creal(out[i*ny*nz_r+j*nz_r+k]);
	    imag = cimag(out[i*ny*nz_r+j*nz_r+k]);
	    E_k[i*ny*nz_r+j*nz_r+k].energy=2*(1/pow(nx*ny*nz,2))*(real*real+imag*imag);
            E_k[i*ny*nz_r+j*nz_r+k].k_sq=(i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+k*k;
         }
      }
   }
   return;
}

void write_file_Ek_binned(int f, int dir, Ek *E_k_binned, Ek *E_k_comp)
{
   FILE *fp;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,dataname[dir]);
   strcat(filename,filenumb);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<no_bins;i++){
     fprintf(fp,"%16.20lf %16.20lf %16.20lf\n", E_k_binned[i].k_sq,E_k_binned[i].energy,E_k_comp[i].energy);
   }
   fclose(fp);
   return;
}
void write_sb(int f, double *sb)
{
   FILE *fp;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,dataname[6]);
   strcat(filename,filenumb);
   strcat(filename,".dbl");
   printf("%s\n",filename);
   fp = fopen(filename,"w");
   fwrite(sb, sizeof(double), nx*ny, fp);
   fclose(fp);
   return;
}
void write_vlos(int f, double *vlos, double *sigvlos)
{
   FILE *fp;
   long int offset;
   char filenumb[5];
   char filename[100];
   int i;
   for(i=0;i<2;i++){
     sprintf(filenumb,"%04d",f);
     strcpy(filename,datdir);
     strcat(filename,dataname[7+i]);
     strcat(filename,filenumb);
     strcat(filename,".dbl");
     printf("%s\n",filename);
     fp = fopen(filename,"w");
     if(i==0) fwrite(vlos, sizeof(double), nx*ny, fp);
     if(i==1) fwrite(sigvlos, sizeof(double), nx*ny, fp);
     fclose(fp);
   }
   return;
}

void write_file_binned(int f, double**arr, int type)
{
   FILE *fp;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   if(type==0) strcat(filename,"temp");
   else strcat(filename,"mach");
   strcat(filename,filenumb);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<no_bins;i++) fprintf(fp,"%16.20lf %16.20lf\n", arr[i][0], arr[i][1]);

   fclose(fp);
   return;
}
