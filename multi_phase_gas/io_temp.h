void read_dbl(int f, int size, double *v)
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
   for(i=0;i<size*nx*ny*nz;i++){
     fread(&d,sizeof(double),1,fp);
     v[i] = d;
   }
   fclose(fp);
   return;
}
void write_file_binned(int f, double**arr, int type)
{
   FILE *fp;
   double d;
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
