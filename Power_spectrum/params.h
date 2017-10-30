#define nx 512
#define ny 512
#define nz 512
#define nv 3
#define prefix "data."
#define suffix ".dbl"
#define datdir "/home/rajsekhar/Final_year_project/HD_Module/Data/higher_k/512/"
#define outprefix "Ek"
#define f1 3
#define f2 4
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

