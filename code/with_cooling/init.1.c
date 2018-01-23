/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sepy 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
   v[PRS] = 1.0;
  #endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;

  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int k, j, i;
  double g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3;
  double *dx1, *dx2, *dx3;
  double sendarray[8], recvarray[8], dvol;
  FILE *hist_file;

  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;
  #ifdef PARALLEL
  if (prank==0)
  #endif
  if (g_stepNumber==0) {
    hist_file = fopen ("pluto_hst.out", "w");
    fprintf(hist_file,"#time g_dt mass TE KE1 KE2 KE3 MOM1 MOM2 MOM3\n ");
  }
  else hist_file = fopen ("pluto_hst.out", "a");

  g_mass=0.0; g_TE=0.0; g_KE1=0.0; g_KE2=0.0; g_KE3=0.0;
  g_mom1=0.0; g_mom2=0.0; g_mom3=0.0;

  DOM_LOOP(k,j,i){
    dvol = dx1[i]*dx2[j]*dx3[k];
    g_mass += d->Vc[RHO][k][j][i]*dvol;
    g_TE += d->Vc[PRS][k][j][i]*dvol/(g_gamma-1.0);
    g_KE1 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]*dvol;
    g_KE2 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]*dvol;
    g_KE3 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]*dvol;
    g_mom1 += d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*dvol;
    g_mom2 += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*dvol;
    g_mom3 += d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*dvol;
  }

  #ifdef PARALLEL
   sendarray[0]=g_mass; sendarray[1]=g_TE; sendarray[2]=g_KE1; sendarray[3]=g_KE2;
   sendarray[4]=g_KE3; sendarray[5]=g_mom1; sendarray[6]=g_mom2; sendarray[7]=g_mom3;
   MPI_Reduce (sendarray, recvarray, 8, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
   if (prank == 0){
     g_mass=recvarray[0]; g_TE=recvarray[1]; g_KE1=recvarray[2]; g_KE2=recvarray[3];
     g_KE3=recvarray[4]; g_mom1=recvarray[5]; g_mom2=recvarray[6]; g_mom3=recvarray[7];
  #endif
  fprintf(hist_file,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n ", g_time
  , g_dt, g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3);
  fclose(hist_file);
  #ifdef PARALLEL
  }
  #endif
}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, k1, k2, k3, dirn;
  const int  KMAX_INT = (int) ( ((int)KMAX) +1 );
  double  *x1, *x2, *x3, *dx1, *dx2, *dx3;
  double kw1, kw2, kw3, phase, momx1, momx2, momx3, mass, dvol;

//  printf("%d %f %f %f\n",KMAX_INT, KMAX, TAU_C, TURB_AMP);

  x1 = grid[IDIR].x; x2 = grid[JDIR].x; x3 = grid[KDIR].x;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;

//  if (g_intStage==1) GetAcc(d); //Get amplitude and phases of acceleration field 
  
  GetAcc(d);

  if (side == 0) {    /* -- check solution inside domain -- */
  
//    printf("%d\n",g_intStage);
    mass=momx1=momx2=momx3=0.0;
    DOM_LOOP(k,j,i){
     dvol = dx1[i]*dx2[j]*dx3[k];
     for (k1=0; k1<=2*KMAX_INT; k1++)
     for (k2=0; k2<=2*KMAX_INT; k2++)
     for (k3=0; k3<=2*KMAX_INT; k3++) {
      kw1 = 2.*CONST_PI*(k1-KMAX_INT); kw2 = 2.*CONST_PI*(k2-KMAX_INT);
      kw3 = 2.*CONST_PI*(k3-KMAX_INT);
      if ( 0. < kw1*kw1+kw2*kw2+kw3*kw3 &&
         kw1*kw1+kw2*kw2+kw3*kw3 <= 4.*CONST_PI*CONST_PI*KMAX*KMAX ) {
//    printf("%d %d %d \n", k1-KMAX_INT, k2-KMAX_INT, k3-KMAX_INT); 
// apply turbulent forcing; only driving large scale modes. Follow Eswaran & Pope 1987
       phase = kw1*x1[i] + kw2*x2[j] + kw3*x3[k];
       momx1 += g_dt*d->Vc[RHO][k][j][i]*( d->Vacc[k3][k2][k1][0]*sin(phase)
                +d->Vacc[k3][k2][k1][3]*cos(phase) )*dvol;
       d->Vc[VX1][k][j][i] += g_dt*( d->Vacc[k3][k2][k1][0]*sin(phase)+d->Vacc[k3][k2][k1][3]*cos(phase) );
       momx2 += g_dt*d->Vc[RHO][k][j][i]*( d->Vacc[k3][k2][k1][1]*sin(phase)
                +d->Vacc[k3][k2][k1][4]*cos(phase) )*dvol;
       d->Vc[VX2][k][j][i] += g_dt*( d->Vacc[k3][k2][k1][1]*sin(phase)+d->Vacc[k3][k2][k1][4]*cos(phase) );
       momx3 += g_dt*d->Vc[RHO][k][j][i]*( d->Vacc[k3][k2][k1][2]*sin(phase)
                +d->Vacc[k3][k2][k1][5]*cos(phase) )*dvol;
       d->Vc[VX3][k][j][i] += g_dt*( d->Vacc[k3][k2][k1][2]*sin(phase)+d->Vacc[k3][k2][k1][5]*cos(phase) );
//       printf("%d %e\n",g_stepNumber,d->Vc[VX1][k][j][i]);
      }
     }
     mass += d->Vc[RHO][k][j][i]*dvol;
    }
    DOM_LOOP(k,j,i){
     d->Vc[VX1][k][j][i] -= momx1/mass;
     d->Vc[VX2][k][j][i] -= momx2/mass;
     d->Vc[VX3][k][j][i] -= momx3/mass;
    }

  }
}
/* ********************************************************************* */
void GetAcc(const Data *d)
{
  int dirn, k1, k2, k3;
  const int  KMAX_INT = (int) ( ((int)KMAX) +1 );
  double ran2(), q1, q2, q3, q4, fac, khdota, modk;
/*  FILE *turb_file;

  if (g_stepNumber==0) {
    turb_file = fopen ("turb.out", "w");
    fprintf(turb_file,"#g_stepNumber,dirn,k1,k2,k3,d->Vacc[k3][k2][k1][dirn],d->Vacc[k3][k2][k1][dirn+3]\n ");
  }
  else turb_file = fopen ("turb.out", "a"); */ 

    for (k3=0; k3<=2*KMAX_INT; k3++)
    for (k2=0; k2<=2*KMAX_INT; k2++)
    for (k1=0; k1<=2*KMAX_INT; k1++) {
      for (dirn=0; dirn<3; dirn++) {
         if (g_stepNumber==0){
//       fprintf(turb_file,"%5i %5i %5i %5i %5i %5i %10.5e %10.5e %10.5e %10.5e %10.5e\n",g_stepNumber,g_intStage,dirn,k1,k2,k3,g_time,d->Vacc[k3][k2][k1][dirn],d->Vacc[k3][k2][k1][dirn+3],q3,sqrt(1.-fac*fac));
           //printf("%d\n",g_iseed);
           q1 = ran2(&g_iseed); //printf("%d\n",g_iseed); 
           q2 = ran2(&g_iseed); 
           q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*CONST_PI*q2);
           q4 = sqrt(-2.0*log(q2+1.0e-20))*cos(2.0*CONST_PI*q1);
//           printf("%d %d %d\n",k1,k2,k3);
           d->Vacc[k3][k2][k1][dirn] = TURB_AMP*q3;
           d->Vacc[k3][k2][k1][dirn+3] = TURB_AMP*q4;
         } else {
           fac = exp(-g_dt/TAU_C);
           //printf("%d\n",g_iseed);
           q1 = ran2(&g_iseed); 
           //printf("%d\n",g_iseed);
           q2 = ran2(&g_iseed); 
           q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*CONST_PI*q2);
           q4 = sqrt(-2.0*log(q2+1.0e-20))*cos(2.0*CONST_PI*q1);
/*      if (k1==2 && k2==2 && k3==2)
       fprintf(turb_file,"%5i %5i %5i %5i %5i %5i %10.5e %10.5e %10.5e %10.5e %10.5e\n",g_stepNumber,g_intStage,dirn,k1,k2,k3,g_time,d->Vacc[k3][k2][k1][dirn],d->Vacc[k3][k2][k1][dirn+3],q3,sqrt(1.-fac*fac)); */
           d->Vacc[k3][k2][k1][dirn] =  TURB_AMP*q3*sqrt(1.-fac*fac)
                                     + d->Vacc[k3][k2][k1][dirn]*fac;
           d->Vacc[k3][k2][k1][dirn+3] = TURB_AMP*q4*sqrt(1.-fac*fac)
                                     + d->Vacc[k3][k2][k1][dirn+3]*fac;
         }
//       printf("here\n");
//       printf("%5i %5i\n",k1, KMAX_INT);
//      if (k1==2 && k2==2 && k3==2) 
//       fprintf(turb_file,"%5i %5i %5i %5i %5i %5i %10.5e %10.5e %10.5e %10.5e %10.5e\n",g_stepNumber,g_intStage,dirn,k1,k2,k3,g_time,d->Vacc[k3][k2][k1][dirn],d->Vacc[k3][k2][k1][dirn+3],q3,sqrt(1.-fac*fac)); 
      }
//make acceleration divergenceless
      khdota = (k1-KMAX_INT)*d->Vacc[k3][k2][k1][0] + (k2-KMAX_INT)*d->Vacc[k3][k2][k1][1]
             + (k3-KMAX_INT)*d->Vacc[k3][k2][k1][2];
      modk = sqrt( (k1-KMAX_INT)*(k1-KMAX_INT) + (k2-KMAX_INT)*(k2-KMAX_INT)
             + (k3-KMAX_INT)*(k3-KMAX_INT) );
      khdota /= (modk+1.e-20);
      d->Vacc[k3][k2][k1][0] -= khdota*(k1-KMAX_INT)/(modk+1.e-20);
      d->Vacc[k3][k2][k1][1] -= khdota*(k2-KMAX_INT)/(modk+1.e-20);
      d->Vacc[k3][k2][k1][2] -= khdota*(k3-KMAX_INT)/(modk+1.e-20);

//      printf("%e %e %e\n",khdota*(k1-KMAX_INT)/(modk+1.e-10),khdota,modk+1.e-10);


      khdota = (k1-KMAX_INT)*d->Vacc[k3][k2][k1][3] + (k2-KMAX_INT)*d->Vacc[k3][k2][k1][4]
             + (k3-KMAX_INT)*d->Vacc[k3][k2][k1][5];
      khdota /= (modk+1.e-20);
      d->Vacc[k3][k2][k1][3] -= khdota*(k1-KMAX_INT)/(modk+1.e-20);
      d->Vacc[k3][k2][k1][4] -= khdota*(k2-KMAX_INT)/(modk+1.e-20);
      d->Vacc[k3][k2][k1][5] -= khdota*(k3-KMAX_INT)/(modk+1.e-20);

    }
//    fclose(turb_file);

}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

/* RAN2 *************************************************************** */
/* Numerical Recipes routine ran2 to generate random number sequence 
 * with a seed */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
