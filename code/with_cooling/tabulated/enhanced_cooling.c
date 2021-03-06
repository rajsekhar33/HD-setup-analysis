/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Take a source step using power-law cooling.

  Integrate the ODE 
  \f[
       dp_{\rm cgs}/dt_{\rm cgs}
      = -(\Gamma - 1) \Lambda(\rho_{\rm cgs}, T_{\rm cgs})
       \qquad{\rm where}\qquad
       \Lambda = \frac{a_{br}}{(\mu m_H)^2} \rho_{\rm cgs}^2 \sqrt{T_{\rm cgs}}
      \,,\qquad
         a_{br} = 2.e-27   c.g.s
  \f]
  which accounts for bremmstrahlung cooling.
  Here the subscript 'cgs' means that the corresponding quantity is given
  in c.g.s units. 
  We denote with \c mu the molecular weight, \c mH the hydrogen mass (in c.g.s).
  
  The previous integration is carried out analytically since the density
  does not change during this step.
  In non-dimensional form:
  \f[
      \frac{dp}{dt} = -{\rm cost} \;\rho  (\rho p)^{\HALF}
  \f]
 
   [notice that since  \c p/rho=T/KELVIN this is equivalent to:
    \c dp/dt=-cost rho^2 (T/KELVIN)^(1/2) ]
 
  The quantity \c cost is determined by transforming the dimensional
  equation into the non-dimensional one.
  If p, rho and t are in code (non-dimensional) units and if \c L_0,
  \c rho_0, and \c V_0 are the unit length, density and velocity,
  then \c cost is found to be:
  \verbatim
                   a_br * (gamma - 1) * L_0 * rho_0
        cost = -------------------------------------------
                sqrt(kB * mu * mH) * kB * mu * mH * V_0^2
  \endverbatim
   where a_{br} = 2.e-27 (in c.g.s), kB is the Boltmann constant
  
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 28, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PowerLawCooling (Data_Arr VV, double dt, Time_Step *Dts, Grid *grid)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  cost, dE;
  double  rho, p, T, p_f, T_f;
  double sendarray, recvarray;
  double *dx1, *dx2, *dx3;
  double dvol;
  double mu=0.5, mue=1.0, mui=1.0;//mu is mean molecular weight

  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;

  cost  = UNIT_LENGTH*UNIT_DENSITY/(UNIT_VELOCITY*UNIT_VELOCITY);
  cost *= 2.e-27*(g_gamma-1.0)*sqrt(mu)/(mue*mui*CONST_mH*sqrt(CONST_mH*CONST_kB));

/*  -------------------------------------------------------------
                Integrate analytically
    -------------------------------------------------------------  */

  dE = 1.e-18;
  g_tot_cool = 0.0;
  DOM_LOOP(k,j,i){

    dvol = dx1[i]*dx2[j]*dx3[k];

/*  ----  Find initial temperature in Kelvin  ----  */

    rho = VV[RHO][k][j][i];
    p   = VV[PRS][k][j][i];

    T   = (mu*p/rho*KELVIN);

    if (T < g_minCoolingTemp || T > g_maxCoolingTemp) continue;

/*  ----  Find final energy  ----  */
    p_f = sqrt(p) - 3.0*0.5*cost*rho*sqrt(rho)*dt;
    //The factor of 3 is introduced here to enhance the cooling
    p_f = MAX(p_f, 0.0);
    p_f = p_f*p_f;

    T_f = mu*p_f/rho*KELVIN;
    T_f = MAX (T_f, g_minCoolingTemp);

/*  ----  Update Energy  ----  */

    p_f = T_f*rho/mu/KELVIN+1.e-8;

    VV[PRS][k][j][i] = p_f;
    dE = fabs(1.0 - p_f/p) + 1.e-18;

    Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE);

    g_tot_cool += (p-p_f)/(g_gamma-1.)/dt*dvol;  
  }

#ifdef PARALLEL
   sendarray = g_tot_cool;
   MPI_Allreduce (&sendarray, &recvarray, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   g_tot_cool = recvarray;
#endif
   print1("T_f = %20.10e, g_tot_cool=%20.30e\n",T_f,  g_tot_cool);
}

double MeanMolecularWeight (double *v)
{
  return(0.5);
}
