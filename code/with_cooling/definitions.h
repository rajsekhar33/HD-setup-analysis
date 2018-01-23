#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 POWER_LAW
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     0
#define  USER_DEF_CONSTANTS      9

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */


/* -- user-defined symbolic constants -- */

#define  UNIT_DENSITY            (CONST_mp*.1)
#define  UNIT_LENGTH             (CONST_pc*40.e3)
#define  UNIT_VELOCITY           (1.e8)
#define  KMAX                    sqrt(2.0)
#define  KMIN                    sqrt(0.0)
#define  TAU_C                   2.0
#define  TURB_AMP                0.1
#define  TAU_F                   TAU_C/100.0
#define  FSEED                   -2001

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       NO
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      YES
#define  SHOCK_FLATTENING       NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          NO
#define  LIMITER                DEFAULT
