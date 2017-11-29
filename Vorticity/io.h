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
	double n[3]={nx,ny,nz};
   for(dir=0;dir<nv;dir++){
      for(i=0;i<nx;i++){
         for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){
			   offset = i*ny*nz+j*ny+k + (dir+1)*nx*ny*nz;
			   fseek(fp,sizeof(double)*offset,SEEK_SET);
			   fread(&d,sizeof(double),1,fp);
			   /*if(dir==0) vel[dir][i][j][k] = ((j/n[1])*(j/n[1])*(j/n[1])*i/n[0])/128;
			   if(dir==1) vel[dir][i][j][k] = (k/n[2])*(k/n[2])/128;
			   if(dir==2) vel[dir][i][j][k] = ((j/n[1])*(j/n[1]))/128;*/		   
			   vel[dir][i][j][k] = d;
			   //printf("%16.12f\n",vel[dir][i][j][k]);
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



void write_output(int type, int f, int n0, int n1, int n2, int n3, double ****vort)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];
   char ordnumb[5];
   sprintf(filenumb,"%04d",f);
   sprintf(ordnumb,"%d",order);
   strcpy(filename,datdir);
   strcat(filename,outprefix);
   strcat(filename,typed);
   strcat(filename,ordnumb);
   strcat(filename,filenumb);
   
   if(type == 0){
     strcat(filename,".dbl");
   } 
   if(type == 1){
     strcat(filename,".vtk");
   }
   printf("\n");
   printf("%s\n",filename);

   int n,i,j,k;
   if(type == 0){
     fp = fopen(filename,"wb");

	for(n=0;n<n0;n++){
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){     
					fwrite(&vort[n][i][j][k],sizeof(double),1,fp);
					//fprintf(fp,"%lf\n",vort[n][i][j][k]);
				}
			}
		}
	}
     fclose(fp);
   }

   return;
}

void print_output(int type, int f, int n0, int n1, int n2, int n3)
{
   FILE *fp;
   double d;
   long int offset;
   char filenumb[5];
   char filename[100];
   char ordnumb[5];
   sprintf(filenumb,"%04d",f);
   sprintf(ordnumb,"%d",order);
   strcpy(filename,datdir);
   strcat(filename,outprefix);
   strcat(filename,typed);
   strcat(filename,ordnumb);
   strcat(filename,filenumb);
   
   if(type == 0){
     strcat(filename,".dbl");
   } 
   if(type == 1){
     strcat(filename,".vtk");
   }
   printf("\n");
   printf("%s\n",filename);

   int n,i,j,k;
   if(type == 0){
     fp = fopen(filename,"rb");

	for(n=0;n<n0;n++){
		for(i=0;i<n1;i++){
			for(j=0;j<n2;j++){
				for(k=0;k<n3;k++){     
					fread(&d,sizeof(double),1,fp);
					printf("%lf\n",d);
				}
			}
		}
	}
     fclose(fp);
   }

   return;
}

