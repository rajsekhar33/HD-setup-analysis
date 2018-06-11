int find_my_i(double Value,  double **binned, int min, int max);

void temp_count(double *Temp, int f){
  
  int i, j, index;
  double ratio2;
  double **temp_binned;
  ratio2=powl(TEMP_max/TEMP_min,(double)1.0/no_bins); 
  temp_binned=  (double **)array2d(no_bins,2,sizeof(double));
  
  for(i=0;i<no_bins;i++){
    temp_binned[i][0]=TEMP_min*powl(ratio2,(double)i);
    temp_binned[i][1]=0.;
  }
  for(i=0;i<nx*ny*nz;i++){
    if(Temp[i]>TEMP_min && Temp[i]<TEMP_max) temp_binned[find_my_i(Temp[i], temp_binned, 0, no_bins)][1]+=1;
  } 
  write_file_binned(f, temp_binned, 0);
}

void mach_count(double *Mach, int f){
  
  int i, j, index;
  double ratio2;
  double **mach_binned;
  ratio2=powl(MACH_max/MACH_min,(double)1.0/no_bins); 
  mach_binned=  (double **)array2d(no_bins,2,sizeof(double));
  for(i=0;i<no_bins;i++){
    mach_binned[i][0]=MACH_min*powl(ratio2,(double)i);
    mach_binned[i][1]=0.;
  }
  for(i=0;i<nx*ny*nz;i++){
    if(Mach[i]>MACH_min && Mach[i]<MACH_max) mach_binned[find_my_i(Mach[i], mach_binned, 0, no_bins)][1]+=1;
  } 
  write_file_binned(f, mach_binned, 1);
}

int find_my_i(double Value, double **binned, int min, int max){
  int pivot=(max+min)/2.0;
  if (min==pivot) return pivot;
  if (Value>binned[pivot][0]) return find_my_i(Value, binned, pivot, max);
  else return find_my_i(Value, binned, min, pivot);
}

void hot_bin(Ek *emission, double vel, double radiat_rate){

  int my_i=(int)((vel-vel_min)*no_hot_bins/2.);
  //printf("%d\n", my_i);
  //printf("%16.20lf, %d, %16.20lf, %16.20lf\n", vel, my_i, emission[my_i].energy, radiat_rate);
  if (my_i>=0 && my_i<no_hot_bins){
    emission[my_i].energy+=radiat_rate;
  }
  //printf("%16.20lf\n", radiat_rate);
}
