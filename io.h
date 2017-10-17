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

   for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
         for(k=0;k<nz;k++){
            for(dir=0;dir<nv;dir++){
			   offset = i*ny*nz+j*ny+k + (dir+1)*nx*ny*nz;
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

void write_E_k(int dir, Ek ****E_k, fftw_complex *out)
{
   int i,j,k;
   double nz_r=0.5*nz+1;
   double real, imag;
   for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
         for(k=0;k<nz_r;k++){
            if (out==NULL) printf("Allocation failure in Output fft array.\n");
	    real = creal(out[i*nx*ny+j*ny+k]);
	    imag = cimag(out[i*nx*ny+j*ny+k]);
	    E_k[dir][i][j][k].energy=(1/pow(nx*ny*nz,2))*(real*real+imag*imag);
            if(i<=nx/2 && j<=ny/2)   E_k[dir][i][j][k].k_sq=i*i+j*j+k*k;
            else E_k[dir][i][j][k].k_sq=(i-nx)*(i-nx)+(j-ny)*(j-ny)+k*k;
         }
      }
   }
   return;
}

void write_output(int f, Ek *E_k_binned, Ek *E_k_comp)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,outprefix);
   strcat(filenumb,outprefix);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<no_bins;i++)  fprintf(fp,"%lf %lf %lf\n",E_k_binned[i].k_sq,E_k_binned[i].energy,E_k_comp[i].energy);
   fclose(fp);
   return;
}

