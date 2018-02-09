int find_my_i(double Value,  double **binned, int min, int max);

void counter(Ek *Data, int dir, int f){
//This function bins fourier transformed daat into bins, and 
//then calls the function to write them into data files 
  int i, j, index;
  double ratio;
  double **binned;
  double K_mid=5e1;
  double diff=(K_mid-K_min)*4/no_bins;
  ratio=powl(K_max/K_mid,(double)(4.0/(3.0*no_bins))); 
  binned=  (double **)array2d(no_bins,3,sizeof(double));
  
  for(i=0;i<no_bins;i++){
    if(i<no_bins/4) binned[i][0]=K_min+i*diff;
    else binned[i][0]=K_mid*powl(ratio,(double)(i-no_bins/4));//Create bins uniform in the logarithmic space for last 3/4th, and uniform in linear space for the first 1/4th points
    binned[i][1]=0.;
    binned[i][2]=0.;
  }
  for(i=0;i<nx*ny*nz_r;i++){
    if(Data[i].k>K_min && Data[i].k<K_max){ 
      index=find_my_i(Data[i].k, binned, 0, no_bins);//Function to find the bin corresponding to a particular k
      if(index>=no_bins){
        printf("Error here.\n");
        exit(0);
      }
      binned[index][1]+=Data[i].energy;//Add energy to specific k bins
      
      if(dir!=4) binned[index][2]+=Data[i].energy*powl(Data[i].k, 5.0/3.0);//This is for compensation
      else binned[index][2]+=Data[i].energy*powl(Data[i].k, 7.0/3.0);
    }
  }
  write_file_binned(f, binned, dir);//This writes the above generated values into a txt file
}

int find_my_i(double Value, double **binned, int min, int max){
//This function divides the domain into halves till the required index is found.
  int pivot=(max+min)/2.0;
  if (min==pivot) return pivot;
  if (Value>binned[pivot][0]) return find_my_i(Value, binned, pivot, max);
  else return find_my_i(Value, binned, min, pivot);
}
