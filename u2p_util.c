
/*! \file u2p_util.c
    \brief U->P inversion functions used by all utoprim files
    
*/


//#include "u2p_defs.h"



/* 
 * reference implementation of transformation from
 * primitive to conserved variables.
 *
 * cfg 7-6-04
 *
 * input: 
 * 
 * primitive variables in 8 element array 
 * metric in contravariant and covariant form
 *
 * output:
 * 
 * conserved variables in 8 element array 
 * 
 */

///  primitive variables 
///  covariant (index dn) form of metric 
///  contravariant (index up) form of metric 
///  matrix of derivatives 


#include "global.grmhd.h"
//#include "decs.h"


static void primtoU_g(struct of_geom *ptrgeom, FTYPE *prim,FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE *U);
static void ucon_calc_g(FTYPE prim[],FTYPE *gcov,FTYPE *gcon,FTYPE ucon[]);
static void raise_g(FTYPE vcov[], FTYPE *gcon, FTYPE vcon[]);
static void lower_g(FTYPE vcon[], FTYPE *gcov, FTYPE vcov[]);
static void ncov_calc(FTYPE *gcon,FTYPE ncov[]) ;
static void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]); 


/// Converts primitive variables to conservative variables 
/// as they are defined with respect to the new primitive 
/// variable solver and not how they are used in the rest of HARM.
/// shouldn't use this function since not optimized to avoid catastrophic cancellations
static void primtoU_g(struct of_geom *ptrgeom, FTYPE *prim,FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE *U)
{
  int i,j ;
  FTYPE rho0 ;
  FTYPE ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],ncov[NDIM] ;
  FTYPE gamma,n_dot_b,bsq,u,p,w ;

  // preliminaries
  ucon_calc_g(prim,gcov,gcon,ucon) ;
  lower_g(ucon,gcov,ucov) ;
  ncov_calc(gcon,ncov) ;

  gamma = -ncov[0]*ucon[0] ;

  bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,gcov,bcov) ;

  n_dot_b = 0. ;
  for(i=0;i<NDIM;i++) n_dot_b += ncov[i]*bcon[i] ;
  bsq = 0. ;
  for(i=0;i<NDIM;i++) bsq += bcov[i]*bcon[i] ;

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho0,u) ;
  w = rho0 + u + p ;

  U[RHO] = gamma*rho0 ;

  for(i=0;i<NDIM;i++) 
    U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
      - (p + bsq/2.)*ncov[i] 
      + n_dot_b*bcov[i] ;

  U[BCON1] = prim[BCON1] ;
  U[BCON2] = prim[BCON2] ;
  U[BCON3] = prim[BCON3] ;

  return ;
}

//// find the contravariant fluid four-velocity from primitive variables plus the metric
static void ucon_calc_g(FTYPE prim[8],FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE ucon[NDIM])
{
  FTYPE u_tilde_con[NDIM] ;
  FTYPE u_tilde_sq ;
  FTYPE gamma,lapse ;
  int i,j ;
 
  u_tilde_con[0] = 0. ;
  u_tilde_con[1] = prim[UTCON1] ;
  u_tilde_con[2] = prim[UTCON2] ;
  u_tilde_con[3] = prim[UTCON3] ;

  u_tilde_sq = 0. ;
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++)
      u_tilde_sq += gcov[GIND(i,j)]*u_tilde_con[i]*u_tilde_con[j] ;
  u_tilde_sq = fabs(u_tilde_sq) ;

  gamma = sqrt(1. + u_tilde_sq) ;

  lapse = sqrt(-1./gcon[GIND(0,0)]) ;

  for(i=0;i<NDIM;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[GIND(0,i)] ;

  return ;
}

///  raise covariant vector vcov using gcon, place result in vcon 
static void raise_g(FTYPE vcov[NDIM], FTYPE gcon[SYMMATRIXNDIM], FTYPE vcon[NDIM])
{
  int i,j;



  for(i=0;i<NDIM;i++) {
    vcon[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcon[i] += gcon[GIND(i,j)]*vcov[j] ;
  }

  return ;
}
///  lower contravariant vector vcon using gcov, place result in vcov 
static void lower_g(FTYPE vcon[NDIM], FTYPE gcov[SYMMATRIXNDIM], FTYPE vcov[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcov[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcov[i] += gcov[GIND(i,j)]*vcon[j] ;
  }

  return ;
}

///  set covariant normal observer four-velocity 
static void ncov_calc(FTYPE gcon[SYMMATRIXNDIM],FTYPE ncov[NDIM])
{
  FTYPE lapse ;

  lapse = sqrt(-1./gcon[GIND(0,0)]) ;

  ncov[0] = -lapse ;
  ncov[1] = 0. ;
  ncov[2] = 0. ;
  ncov[3] = 0. ;

  return ;
}



///  set covariant normal observer four-velocity 
static void ncov_calc_fromlapse(FTYPE lapse,FTYPE ncov[NDIM]) 
{

  ncov[0] = -lapse ;
  ncov[1] = 0. ;
  ncov[2] = 0. ;
  ncov[3] = 0. ;

  return ;
}

///  calculate contravariant magnetic field four-vector b 
static void bcon_calc_g(FTYPE prim[8],FTYPE ucon[NDIM],FTYPE ucov[NDIM],FTYPE ncov[NDIM],FTYPE bcon[NDIM]) 
{
  FTYPE Bcon[NDIM] ;
  FTYPE u_dot_B ;
  FTYPE gamma ;
  int i ;

  Bcon[0] = 0. ;
  for(i=1;i<NDIM;i++) Bcon[i] = prim[BCON1+i-1] ;

  u_dot_B = 0. ;
  for(i=0;i<NDIM;i++) u_dot_B += ucov[i]*Bcon[i] ;

  gamma = -ucon[0]*ncov[0] ;
  for(i=0;i<NDIM;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


