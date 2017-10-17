typedef struct EK{
  double energy;//Stores the value of energy at this particular point in k space
  double k_sq;  //Stores the value of k_sq at this point in k space
} Ek;
const double pi=CONST_PI;

void counter(Ek *E_k, Ek *E_k_added){
  //This function adds up values of E_k for the same value of k
  int i;
  for(i=0;i<(nx*nx+ny*ny+nz*nz)/2+1;i++){
    E_k_added[i].energy=0;
    E_k_added[i].k_sq=2*pi*sqrt((double)i);//From here onwards k_sq stores the actual k value, not k^2
  }
  for(i=0;i<(nx*nx+ny*ny+nz*nz)/2+1;i++){
    if(E_k[i].energy!=0){  
      E_k_added[(int)E_k[i].k_sq].energy+=E_k[i].energy; //Here we add up E(k) for same value of k
    }
  }
}
//The following function bins E(k) vs k data
void bin(Ek *E_k_added, Ek *E_k_binned){
  int index=0;
  int i,j;
  int bin_no=0;//bin_no indicates the bin in which we are storing our data
  double max_k, r_bin1, r_bin2, k_min, k_max, delta_k, E;
  //max_k stores the maximum value of k in the entire k space
  max_k=2*pi*sqrt(((double)nx*nx+ny*ny+nz*nz)/2);
  //r_bin1 sets the ratio between cnsecutive bin sizes for stretched binning method
  r_bin1=pow(max_k/ratio,1/(9*no_bins/10));
  //r_bin2 sets the ratio between consecutive k values for logarithmic binning method
  r_bin2=pow(ratio,(10/no_bins));
  //each bin is bounded by k_min and k_max
  k_max=2*bin_size*pi-bin_size*pi;
  for (i=1;i<no_bins;i++){
    k_min=k_max;
    //for stretched binning method, at lower k values
    if (k_min<=max_k/ratio) k_max=k_max+2*bin_size*pi*pow(r_bin1,(double)i);
    //for logarithmic binning method, at higher k values
    else  k_max=k_max*r_bin2;
    //leave the loop if k_max exceeds maximum possible k value
    if (k_max>max_k) break;
    delta_k=k_max-k_min;
    E=0;
    for (j=index;j<(nx*nx+ny*ny+nz*nz)/2+1 ;j++){
      if (E_k_added[j].k_sq<k_max) E+=E_k_added[j].energy;
        //add up values in each bin
      else{
      //add the entry to the array
        E_k_binned[bin_no].k_sq=k_min;
        E_k_binned[bin_no].energy=E/(delta_k);
        bin_no++;
        if(E!=0) index=j;
        break;
      }
    }
  }
}
//The following function calculates the compensated power spectrum
void calc_comp(Ek *E_k_comp, Ek *E_k_binned){
  int tot_bins=no_bins;
  int i;
  for(i=0;i<tot_bins;i++){
    E_k_comp[i].k_sq=E_k_binned[i].k_sq;
    E_k_comp[i].energy=E_k_binned[i].energy * pow(E_k_binned[i].k_sq,5/3);
  }
}
