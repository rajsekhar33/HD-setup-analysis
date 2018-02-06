void read_dbl(int f, int dir, double *v)
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
   for(i=0;i<nx*ny*nz;i++){
     offset = i + dir*nx*ny*nz;
     fseek(fp,sizeof(double)*offset,SEEK_SET);
     fread(&d,sizeof(double),1,fp);
     //if(dir==1 && d>0.) {printf("No Error here %lf %d \n",d, offset);}
     if(dir==1 && d<0.) {printf("Error here 1 %lf %d \n",d, offset);exit(0);}
     v[i] = d;
   }
   fclose(fp);
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
