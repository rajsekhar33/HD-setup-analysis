#define nx 256
#define ny 256
#define nz 256
#define nv 5
#define prefix "data."
#define suffix ".dbl"
#define datdir "/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k6-8/"
#define f1 40 
#define f2 45
#define CONST_PI 3.14159265358979 
#define no_bins 200 //number of bins in k space that we divide our fourier transformed data among
#define  UNIT_VELOCITY (1.e8)
#define CONST_mp      1.67262171e-24     /**<  Proton mass.            */
#define CONST_kB      1.3806505e-16      /**<  Boltzmann constant.     */
#define CONST_mu      0.5      /**<  average mass in CONST_mp units    */
#define TEMP_max      1e8      /**<  Maximum temperature that we consider in our distribution        */
#define TEMP_min      1e4      /**<  Minimum temperature that we consider in our distribution        */
#define gamma         5./3.    //Value of gamma
#define MACH_max      1e2      /**<  Maximum mach no. that we consider in our distribution        */
#define MACH_min      1e-2     /**<  Minimum mach no. that we consider in our distribution        */
#define ratio 10.0
#define bin_size 1



int nz_r=0.5*nz+1;

char dataname[5][100]={"Rhok", "Ekx", "Eky", "Ekz", "Prsk"};


//The following structure stores both |k| value and the corresponding energy
typedef struct EK{
  double energy;//Stores the value of energy at this particular point in k space
  double k_sq;  //Stores the value of k_sq at this point in k space
} Ek;

