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
	   offset = i*ny*nz+j*ny+k + dir*nx*ny*nz;
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

void write_file_temp_binned(int f, double**temp)
{
   FILE *fp;
   double d;
   char filenumb[5];
   char filename[100];

   sprintf(filenumb,"%04d",f);
   strcpy(filename,datdir);
   strcat(filename,"temp");
   strcat(filename,filenumb);
   strcat(filename,".txt");
   printf("%s\n",filename);
   int i;
   fp = fopen(filename,"w");
   for(i=0;i<no_bins;i++) fprintf(fp,"%16.20lf %16.20lf\n", temp[i][0], temp[i][1]);

   fclose(fp);
   return;
}
