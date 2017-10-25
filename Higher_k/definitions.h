#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     0
#define  USER_DEF_CONSTANTS      4

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */


/* -- user-defined symbolic constants -- */

#define  KMAX                    sqrt(36.0)
#define  KMIN                    sqrt(16.0)
#define  TAU_C                   2.0
#define  TURB_AMP                0.1

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       NO
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      YES
#define  SHOCK_FLATTENING       NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          NO
#define  LIMITER                DEFAULT
