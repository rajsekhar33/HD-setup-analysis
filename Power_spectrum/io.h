void read_dbl(int f, double ****vel)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];
   int i,j,k,dir;
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
   for(dir=0;dir<nv;dir++){
     for(i=0;i<nx;i++){
       for(j=0;j<ny;j++){
         for(k=0;k<nz;k++){
	   offset = i*ny*nz+j*ny+k + (dir+2)*nx*ny*nz;
	   fseek(fp,sizeof(double)*offset,SEEK_SET);
	   fread(&d,sizeof(double),1,fp);
	   vel[dir][i][j][k] = d;
	   }
         }
      }
   }
   fclose(fp);
   return;
}

void printvk(fftw_complex *out, int nz_r, int f, int j, int k){
  int i;
  FILE *fp;
  double d;
  long int offset;
  char filenumb[5];
  char filename[100];
  sprintf(filenumb,"%04d",f);
  strcpy(filename,"Vk_shifted");
  strcat(filename,filenumb);
  strcat(filename,".txt");
  printf("%s\n",filename);
  fp = fopen(filename,"w");
  offset=ny*nz_r*j+nz_r*k;
  for(i=0;i<nz_r;i++) fprintf(fp,"%16.20lf\n",cabs(out[offset+i])/(nx*ny*nz_r));
  fclose(fp);
  return;
}

void printarray4d(int n0, int n1, int n2, int n3, double ****vel)
{
   int n,i,j,k;
   for(n=0;n<n0;n++){
       for(i=0;i<n1;i++){
          for(j=0;j<n2;j++){
             for(k=0;k<n3;k++){
                printf("%d,%2d,%2d,%2d,%16.12f\n",n,i,j,k,vel[n][i][j][k]);
             }
          }
       }
    }
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
//This function calculate Vk_i^2 along the three directions
void write_E_k(int dir, Ek ****E_k, fftw_complex *out)
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
	    E_k[dir][i][j][k].energy=(1/pow(nx*ny*nz,2))*(real*real+imag*imag);
            E_k[dir][i][j][k].k_sq=(i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+k*k;
         }
      }
   }
   return;
}

void write_file_Ek_binned(int f, Ek *E_k_binned, Ek *E_k_comp)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,outprefix);
   strcat(filename,filenumb);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<no_bins;i++) fprintf(fp,"%16.20lf %16.20lf %16.20lf\n",E_k_binned[i].k_sq,E_k_binned[i].energy,E_k_comp[i].energy);

   fclose(fp);
   return;
}
void write_file_Ek(int f, Ek *E_k_added)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,"E_k_added");
   strcat(filename,filenumb);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<(nx*nx+ny*ny+nz*nz)/2+1;i++)  fprintf(fp,"%16.20lf %16.20lf\n",E_k_added[i].k_sq,E_k_added[i].energy);
   fclose(fp);
   return;
}
