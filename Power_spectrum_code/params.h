#define nx 256
#define ny 256
#define nz 256
#define nv 3
#define prefix "data."
#define suffix ".dbl"
#define datdir "/mnt/lustre/ug4/ugrajs/cooling/higher_k/256/k12/"
#define outprefix "Ek"
#define f1 40
#define f2 50
#define verbose 0
#define fn 3 
#define nu 2
#define bin_size 1
#define CONST_PI 3.14159265358979 
#define ratio 10.0
#define no_bins 200

//The following structure stores both |k| value and the corresponding energy
typedef struct EK{
  double energy;//Stores the value of energy at this particular point in k space
  double k_sq;  //Stores the value of k_sq at this point in k space
} Ek;
