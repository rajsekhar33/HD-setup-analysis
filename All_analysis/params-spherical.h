#define nx 256
#define ny 256
#define nz 256
#define nv 5
#define prefix "data."
#define suffix ".dbl"
#define datdir "/mnt/lustre/ug4/ugrajs/cooling/turb_perturb/DBh2e-1/F7e-1/"
#define f1 1
#define f2 142 
#define fstep 1 
#define ntrc 0       //Number of tracers

#define CONST_PI      3.14159265358979   /**<  \f$ \pi \f$.               */
#define CONST_mp      1.67262171e-24     /**<  Proton mass.            */
#define CONST_kB      1.3806505e-16      /**<  Boltzmann constant.     */
#define CONST_mu      0.5      /**<  average mass in CONST_mp units    */
#define CONST_pc      3.0856775807e18    /**<  Parsec.                    */

#define UNIT_VELOCITY (1.e8)
#define UNIT_DENSITY  (CONST_mp*.1)
#define UNIT_LENGTH   (CONST_pc*40.e3)
#define gamma         5./3.    //Value of gamma

#define TEMP_max      1e9      /**<  Maximum temperature that we consider in our distribution        */
#define TEMP_min      1e4      /**<  Minimum temperature that we consider in our distribution        */
#define MACH_max      1e2      /**<  Maximum mach no. that we consider in our distribution        */
#define MACH_min      1e-3     /**<  Minimum mach no. that we consider in our distribution        */

#define ratio         10.0
#define no_bins       200 //number of bins in k space that we divide our fourier transformed data among
#define no_hot_bins   200 //number of bins into which we divide vlos-emission data 
#define lambda(T)     2e-27*sqrt(T)

//For introducing spherical density profile
#define sigma         0.2 //Sigma denotes decay length of tanh function
#define r0            0.5 //r0 denotes radius of the sphere beyond which we set
                          //density perturbations to zero
#define rho0          1.0 //Unperturbed density at each grid point
#define vel_min       -1.0 //minimum los_velocity taken for calculation
#define vel_max       1.0 //minimum los_velocity taken for calculation
double del_v=(vel_max-vel_min)/(double)no_hot_bins; //bin size for hot gas vlos

int nz_r=0.5*nz+1;
int ny_r=0.5*ny+1;

char dataname[10][100]={"Rhoks", "Ekx", "Eky", "Ekz", "Prsk", "Trck", "sbs", "vlos", "sigvlos", "radiat_vel"};
char sbname[100]={"sbks"};

//The following structure stores both |k| value and the corresponding energy
typedef struct EK{
  double energy;//Stores the value of energy at this particular point in k space
  double k_sq;  //Stores the value of k_sq at this point in k space
} Ek;
