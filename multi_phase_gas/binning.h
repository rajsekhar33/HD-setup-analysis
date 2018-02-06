void swap(double* a, double* b);
int partition (double *arr, int low, int high);
void quickSort(double *arr, int low, int high);

void temp_count(double *temp, int f){
  int i, j, index;
  double ratio2;
  double **temp_binned;
  ratio2=powl(TEMP_max/TEMP_min,(double)1.0/no_bins); 
  temp_binned=  (double **)array2d(no_bins,2,sizeof(double));
  quickSort(temp, 0, nx*ny*nz);
  printf("%f %f\n",temp[0],temp[1]);
  for(i=0;i<nx*ny*nz;i++){
    if(temp[i]>TEMP_min){
      index=i;
      break;
    }
  }
  for(i=0;i<no_bins;i++){
    temp_binned[i][0]=TEMP_min*powl(ratio2,(double)i);
    temp_binned[i][1]=0.;
    for(j=index;j<nx*ny*nz;j++){
      if(temp[j]>TEMP_max) break;
      else if(temp[j]<temp_binned[i][0]) temp_binned[i][1]+=1;
      else{
        index=j+1;
        break;
      }
    }
  }
  
  write_file_temp_binned(f, temp_binned);
}
// A utility function to swap two elements
void swap(double* a, double* b){
    double t = *a;
    *a = *b;
    *b = t;
}
 
/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (double *arr, int low, int high){
    double pivot = arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
    int j;
    for (j = low; j <= high- 1; j++){
        // If current element is smaller than or
        // equal to pivot
        if (arr[j] <= pivot){
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
 
/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quickSort(double *arr, int low, int high){
    if(low < high){
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

