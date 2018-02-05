#define nx 512
#define ny 512
#define nz 512
#define nv 2
#define prefix "data."
#define suffix ".dbl"
#define datdir "/mnt/lustre/ug4/ugrajs/cooling/higher_k/512/k8-10/"
#define outprefix "delPk"
#define f1 20
#define f2 30
#define verbose 0
#define fn 3 
#define nu 2
#define bin_size 1
#define CONST_PI 3.14159265358979 
#define ratio 10.0
#define no_bins 200
#define  UNIT_VELOCITY (1.e8)
#define CONST_mp      1.67262171e-24     /**<  Proton mass.               */
#define CONST_kB      1.3806505e-16      /**<  Boltzmann constant.        */
#define CONST_mu      0.5      /**<  average mass in CONST_mp units        */
#define TEMP_max      1e8      /**<  Minimum temperature that we consider in our distribution        */
#define TEMP_min      1e6      /**<  Maximum temperature that we consider in our distribution        */


//The following structure stores both |k| value and the corresponding energy
typedef struct EK{
  double energy;//Stores the value of energy at this particular point in k space
  double k_sq;  //Stores the value of k_sq at this point in k space
} Ek;

