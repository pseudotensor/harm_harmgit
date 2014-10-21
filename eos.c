#include "decs.h"

/*! \file eos.c
  \brief Equation of State general functions

  ////////////////////////////
  //
  // IMPLEMENT EOS
  //
  // u_rho0_p : used by initial conditions
  //
  // u_rho0_T : used by initial conditions
  //
  // pressure_rho0_u : used by inversion for initial guess and by rest of code to set pressure as functions of rho0 and u
  // 
  // dpdu_rho0_u  dpdrho0_rho0_u : used by sources or other such derivatives
  //
  // cs2_compute : used by vchar.c for characteristics
  //
  // pressure_wmrho0 : used by inversion
  //
  // compute_idwmrho0dp : used by 2D inversion
  //
  //
  // compute_entropy and compute_u_from_entropy : used by entropy evolution and inversion
  //
  ////////////////////////////
  */


#define OLDCALC 0

// below not used anymore since want compiler to know all files are included so don't have to repeat "make prep"
// now #if is inside each file
// this chooses the equation of state
//#if(WHICHEOS==GRBPWF99)
//#include "grbpwf99eos.c"
//#elif(WHICHEOS==IDEALGAS)
//#include "idealgaseos.c"
//#elif(WHICHEOS==MIGNONE)
//#include "mignoneeos.c"
//#elif(WHICHEOS==KAZFULL)
//#include "kazfulleos.c"
//#endif

// always include all of the ever-wanted EOSs so included in dependencies when compiling
#include "coldeos.c"
#include "idealgaseos.c"
#include "mignoneeos.c"
#include "grbpwf99eos.c"
#include "kazfulleos.c"






//////////////////////////////////////////////////
//
// wrappers for any EOS and any EOMTYPE
// Currently 20 wrapper functions
//
//////////////////////////////////////////////////

/// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_pressure_rho0_u[whicheos]))(EOSextra,rho0,u) );
}

/// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return( (*(ptr_u_rho0_p[whicheos]))(EOSextra,rho0,p) );
}

/// u(rho0, T) (used for initial conditions)
FTYPE u_rho0_T(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE T)
{
  return( (*(ptr_u_rho0_T[whicheos]))(EOSextra,rho0,T) );
}

/// dp(rho0, u)/du
FTYPE dpdu_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_dpdu_rho0_u[whicheos]))(EOSextra,rho0,u) );
}

/// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_dpdrho0_rho0_u[whicheos]))(EOSextra,rho0,u) );
}

/// sound speed squared (for vchar.c)
FTYPE cs2_compute(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_cs2_compute[whicheos]))(EOSextra,rho0,u) );
}

/// entropy as function of rho0 and internal energy (u)
/// S(rho0,u)
FTYPE compute_entropy(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_compute_entropy[whicheos]))(EOSextra,rho0,u) );
}

/// u(rho0,S)
FTYPE compute_u_from_entropy(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  return( (*(ptr_compute_u_from_entropy[whicheos]))(EOSextra,rho0,entropy) );
}

/// used for dudp_calc
FTYPE compute_dSdrho(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_compute_dSdrho[whicheos]))(EOSextra,rho0,u) );
}


/// used for dudp_calc
FTYPE compute_dSdu(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_compute_dSdu[whicheos]))(EOSextra,rho0,u) );
}




/// specific entropy as function of rho0 and internal energy (wmrho0)
/// specificS(rho0,wmrho0)
FTYPE compute_specificentropy_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_compute_specificentropy_wmrho0[whicheos]))(EOSextra,rho0,wmrho0) );
}

/// used for utoprim_jon entropy inversion
FTYPE compute_dspecificSdrho_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_compute_dspecificSdrho_wmrho0[whicheos]))(EOSextra,rho0,wmrho0) );
}

/// used for utoprim_jon entropy inversion
FTYPE compute_dspecificSdwmrho0_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_compute_dspecificSdwmrho0_wmrho0[whicheos]))(EOSextra,rho0,wmrho0) );
}




/// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_pressure_wmrho0[whicheos]))(EOSextra,rho0,wmrho0) );
}


/// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_compute_idwmrho0dp[whicheos]))(EOSextra,rho0,wmrho0) );
}


/// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*(ptr_compute_idrho0dp[whicheos]))(EOSextra,rho0,wmrho0) );
}


/// radiation rate
FTYPE compute_qdot(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_compute_qdot[whicheos]))(EOSextra,rho0,u) );
}

/// compute source terms for EOS
int compute_sources_EOS(int whicheos, FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  return( (*(ptr_compute_sources_EOS[whicheos]))(EOSextra,pr, geom, q, Ui, dUother, dUcomp));
}


/// compute any and all extra quantities for the EOS
void compute_allextras(int whicheos, int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  (*(ptr_compute_allextras[whicheos]))(justnum,EOSextra,rho0,u,numextrasreturn,extras);
  return;
}

/// compute all extras *and* processed quantities for this EOS
int get_extrasprocessed(int whicheos, int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  return((*(ptr_get_extrasprocessed[whicheos]))(doall,EOSextra,pr,extras,processed));
}


/// compute temperature
FTYPE compute_temp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*(ptr_compute_temp[whicheos]))(EOSextra,rho0,u) );
}

/// compute EOS parameters
void compute_EOS_parms(int whicheos, FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  (*(ptr_compute_EOS_parms[whicheos]))(EOSextra,prim);
}

/// compute EOS all parameters 
void compute_EOS_parms_full(int whicheos, FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  (*(ptr_compute_EOS_parms_full[whicheos]))(EOSextra,prim);
}

/// store EOS parameters
void store_EOS_parms(int whicheos, int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  (*(ptr_store_EOS_parms[whicheos]))(numparms,EOSextra,parlist);
}

/// get EOS parameters
void get_EOS_parms(int whicheos, int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  (*(ptr_get_EOS_parms[whicheos]))(numparms, EOSextra, parlist);
}

/// fix primitives for this eos
void fix_primitive_eos_scalars(int whicheos, FTYPE *EOSextra, FTYPE *pr)
{
  (*(ptr_fix_primitive_eos_scalars[whicheos]))(EOSextra, pr);
}

/// get all quantities needed for inversion of conserved -> primitive
void getall_forinversion(int whicheos, int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{
  return( (*(ptr_getall_forinversion[whicheos]))(eomtype, whichd, EOSextra,quant1,quant2,fun,dfunofrho,dfunofu) );
}



//////////////////////////////////////////
//
// ANY EOS
//


/// old function
/// p(rho0,w)
FTYPE pressure_rho0_w(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE w)
{
  FTYPE wmrho0=w-rho0;

#if(OLDCALC)
  return((GAMMA-1.)*(w - rho0)/GAMMA) ;
#else
  return(pressure_wmrho0(whicheos,EOSextra,rho0,wmrho0)) ;
#endif
}


/// pick EOMTYPE for EOS
int pickeos_eomtype(int whicheosinput, int whicheom, int *whicheosoutput)
{
  //COLDEOS
  //IDEALGAS
  //MIGNONE
  //GRBPWF99
  //KAZFULL



  // EOS functions used during inversion and other places
  if(whicheom==EOMCOLDGRMHD||whicheom==EOMFFDE||whicheom==EOMFFDE2){
    *whicheosoutput=COLDEOS;
  }
  else{
    *whicheosoutput=whicheosinput;
  }




  if(*whicheosoutput==IDEALGAS){
    // then need to set gamideal
    gamideal=gam; // OPENMPMARK: Since global, would need to make this threadprivate .. figure out how to avoid.
    // OPENMPMARK: Expensive for OpenMP to have to allow all those functions for EOS to be different for each thread (which generally be allowed).  So force pickeos_eomtype() to be called outside parallel region so can avoid having EOS functions ptr's as threadprivate
  }
  else{
    // otherwise assume gam and gamideal could be different, where gam is used to true EOS for some purpose while gamideal is particular always for ideal EOS
  }

  return(0);

}




/// setup array of pointers to point to different EOSs so can access later without use of global variables.  For OpenMP global vars would have to be threadprivate since each thread may use a different EOS, and using threadprivate is very slow.
int initeos_eomtype(void)
{
  int whicheos;

  //////////////////////////////////////////////////////
  whicheos=IDEALGAS;
  ptr_pressure_rho0_u[whicheos] = &pressure_rho0_u_idealgas;
  ptr_compute_u_from_entropy[whicheos] = &compute_u_from_entropy_idealgas;
  ptr_u_rho0_p[whicheos] = &u_rho0_p_idealgas;
  ptr_u_rho0_T[whicheos] = &u_rho0_T_idealgas;
  ptr_dpdu_rho0_u[whicheos] = &dpdu_rho0_u_idealgas;
  ptr_dpdrho0_rho0_u[whicheos] = &dpdrho0_rho0_u_idealgas;
  ptr_cs2_compute[whicheos] = &cs2_compute_idealgas;

  ptr_compute_dSdrho[whicheos] = &compute_dSdrho_idealgas;
  ptr_compute_dSdu[whicheos] = &compute_dSdu_idealgas;
  ptr_compute_entropy[whicheos] = &compute_entropy_idealgas;

  ptr_compute_dspecificSdrho_wmrho0[whicheos] = &compute_dspecificSdrho_wmrho0_idealgas;
  ptr_compute_dspecificSdwmrho0_wmrho0[whicheos] = &compute_dspecificSdwmrho0_wmrho0_idealgas;
  ptr_compute_specificentropy_wmrho0[whicheos] = &compute_specificentropy_wmrho0_idealgas;
      
  ptr_pressure_wmrho0[whicheos] = &pressure_wmrho0_idealgas;
  ptr_compute_idwmrho0dp[whicheos] = &compute_idwmrho0dp_idealgas;
  ptr_compute_idrho0dp[whicheos] = &compute_idrho0dp_idealgas;
      
  ptr_compute_qdot[whicheos] = &compute_qdot_idealgas;
  ptr_compute_sources_EOS[whicheos] = &compute_sources_EOS_idealgas;
  ptr_compute_allextras[whicheos] = &compute_allextras_idealgas;
  ptr_get_extrasprocessed[whicheos] = &get_extrasprocessed_idealgas;
      
  ptr_compute_temp[whicheos] = &compute_temp_idealgas;
      
  ptr_compute_EOS_parms[whicheos] = &compute_EOS_parms_idealgas;
  ptr_compute_EOS_parms_full[whicheos] = &compute_EOS_parms_idealgas; // same as non-full
  ptr_store_EOS_parms[whicheos] = &store_EOS_parms_idealgas;
  ptr_get_EOS_parms[whicheos] = &get_EOS_parms_idealgas;
  ptr_fix_primitive_eos_scalars[whicheos] = &fix_primitive_eos_scalars_idealgas;
  ptr_getall_forinversion[whicheos] = &getall_forinversion_idealgas;


  //////////////////////////////////////////////////////
  whicheos=MIGNONE;
  ptr_pressure_rho0_u[whicheos] = &pressure_rho0_u_mignone;
  ptr_compute_u_from_entropy[whicheos] = &compute_u_from_entropy_mignone;
  ptr_u_rho0_T[whicheos] = &u_rho0_T_mignone;
  ptr_dpdu_rho0_u[whicheos] = &dpdu_rho0_u_mignone;
  ptr_dpdrho0_rho0_u[whicheos] = &dpdrho0_rho0_u_mignone;
  ptr_cs2_compute[whicheos] = &cs2_compute_mignone;

  ptr_compute_dSdrho[whicheos] = &compute_dSdrho_mignone;
  ptr_compute_dSdu[whicheos] = &compute_dSdu_mignone;
  ptr_compute_entropy[whicheos] = &compute_entropy_mignone;

  ptr_compute_dspecificSdrho_wmrho0[whicheos] = &compute_dspecificSdrho_wmrho0_mignone;
  ptr_compute_dspecificSdwmrho0_wmrho0[whicheos] = &compute_dspecificSdwmrho0_wmrho0_mignone;
  ptr_compute_specificentropy_wmrho0[whicheos] = &compute_specificentropy_wmrho0_mignone;
      
  ptr_pressure_wmrho0[whicheos] = &pressure_wmrho0_mignone;
  ptr_compute_idwmrho0dp[whicheos] = &compute_idwmrho0dp_mignone;
  ptr_compute_idrho0dp[whicheos] = &compute_idrho0dp_mignone;
      
  ptr_compute_qdot[whicheos] = &compute_qdot_mignone;
  ptr_compute_sources_EOS[whicheos] = &compute_sources_EOS_mignone;
  ptr_compute_allextras[whicheos] = &compute_allextras_mignone;
  ptr_get_extrasprocessed[whicheos] = &get_extrasprocessed_mignone;

  ptr_compute_temp[whicheos] = &compute_temp_mignone;
      
  ptr_compute_EOS_parms[whicheos] = &compute_EOS_parms_mignone;
  ptr_compute_EOS_parms_full[whicheos] = &compute_EOS_parms_mignone; // same as non-full
  ptr_store_EOS_parms[whicheos] = &store_EOS_parms_mignone;
  ptr_get_EOS_parms[whicheos] = &get_EOS_parms_mignone;
  ptr_fix_primitive_eos_scalars[whicheos] = &fix_primitive_eos_scalars_mignone;
  ptr_getall_forinversion[whicheos] = &getall_forinversion_mignone;

  //////////////////////////////////////////////////////
  whicheos=GRBPWF99;
  ptr_pressure_rho0_u[whicheos] = &pressure_rho0_u_grbpwf99;
  ptr_compute_u_from_entropy[whicheos] = &compute_u_from_entropy_grbpwf99;
  ptr_u_rho0_T[whicheos] = &u_rho0_T_grbpwf99;
  ptr_dpdu_rho0_u[whicheos] = &dpdu_rho0_u_grbpwf99;
  ptr_dpdrho0_rho0_u[whicheos] = &dpdrho0_rho0_u_grbpwf99;
  ptr_cs2_compute[whicheos] = &cs2_compute_grbpwf99;

  ptr_compute_dSdrho[whicheos] = &compute_dSdrho_grbpwf99;
  ptr_compute_dSdu[whicheos] = &compute_dSdu_grbpwf99;
  ptr_compute_entropy[whicheos] = &compute_entropy_grbpwf99;

  ptr_compute_dspecificSdrho_wmrho0[whicheos] = &compute_dspecificSdrho_wmrho0_grbpwf99;
  ptr_compute_dspecificSdwmrho0_wmrho0[whicheos] = &compute_dspecificSdwmrho0_wmrho0_grbpwf99;
  ptr_compute_specificentropy_wmrho0[whicheos] = &compute_specificentropy_wmrho0_grbpwf99;
      
  ptr_pressure_wmrho0[whicheos] = &pressure_wmrho0_grbpwf99;
  ptr_compute_idwmrho0dp[whicheos] = &compute_idwmrho0dp_grbpwf99;
  ptr_compute_idrho0dp[whicheos] = &compute_idrho0dp_grbpwf99;
      
  ptr_compute_qdot[whicheos] = &compute_qdot_grbpwf99;
  ptr_compute_sources_EOS[whicheos] = &compute_sources_EOS_grbpwf99;
  ptr_compute_allextras[whicheos] = &compute_allextras_grbpwf99;
  ptr_get_extrasprocessed[whicheos] = &get_extrasprocessed_grbpwf99;

      
  ptr_compute_temp[whicheos] = &compute_temp_grbpwf99;
      
  ptr_compute_EOS_parms[whicheos] = &compute_EOS_parms_grbpwf99;
  ptr_compute_EOS_parms_full[whicheos] = &compute_EOS_parms_grbpwf99; // same as non-full
  ptr_store_EOS_parms[whicheos] = &store_EOS_parms_grbpwf99;
  ptr_get_EOS_parms[whicheos] = &get_EOS_parms_grbpwf99;
  ptr_fix_primitive_eos_scalars[whicheos] = &fix_primitive_eos_scalars_grbpwf99;
  ptr_getall_forinversion[whicheos] = &getall_forinversion_grbpwf99;


  //////////////////////////////////////////////////////
  whicheos=KAZFULL;
  ptr_pressure_rho0_u[whicheos] = &pressure_rho0_u_kazfull;
  ptr_compute_u_from_entropy[whicheos] = &compute_u_from_entropy_kazfull;
  ptr_u_rho0_T[whicheos] = &u_rho0_T_kazfull;
  ptr_dpdu_rho0_u[whicheos] = &dpdu_rho0_u_kazfull;
  ptr_dpdrho0_rho0_u[whicheos] = &dpdrho0_rho0_u_kazfull;
  ptr_cs2_compute[whicheos] = &cs2_compute_kazfull;

  ptr_compute_dSdrho[whicheos] = &compute_dSdrho_kazfull;
  ptr_compute_dSdu[whicheos] = &compute_dSdu_kazfull;
  ptr_compute_entropy[whicheos] = &compute_entropy_kazfull;

  ptr_compute_dspecificSdrho_wmrho0[whicheos] = &compute_dspecificSdrho_wmrho0_kazfull;
  ptr_compute_dspecificSdwmrho0_wmrho0[whicheos] = &compute_dspecificSdwmrho0_wmrho0_kazfull;
  ptr_compute_specificentropy_wmrho0[whicheos] = &compute_specificentropy_wmrho0_kazfull;
            
  ptr_pressure_wmrho0[whicheos] = &pressure_wmrho0_kazfull;
  ptr_compute_idwmrho0dp[whicheos] = &compute_idwmrho0dp_kazfull;
  ptr_compute_idrho0dp[whicheos] = &compute_idrho0dp_kazfull;
      
  ptr_compute_qdot[whicheos] = &compute_qdot_kazfull;
  ptr_compute_sources_EOS[whicheos] = &compute_sources_EOS_kazfull;
  ptr_compute_allextras[whicheos] = &compute_allextras_kazfull;
  ptr_get_extrasprocessed[whicheos] = &get_extrasprocessed_kazfull;

      
  ptr_compute_temp[whicheos] = &compute_temp_kazfull;
      
  ptr_compute_EOS_parms[whicheos] = &compute_EOS_parms_kazfull;
  ptr_compute_EOS_parms_full[whicheos] = &compute_EOS_parms_full_kazfull; // different than non-full
  ptr_store_EOS_parms[whicheos] = &store_EOS_parms_kazfull;
  ptr_get_EOS_parms[whicheos] = &get_EOS_parms_kazfull;
  ptr_fix_primitive_eos_scalars[whicheos] = &fix_primitive_eos_scalars_kazfull;
  ptr_getall_forinversion[whicheos] = &getall_forinversion_kazfull;


  //////////////////////////////////////////////////////
  whicheos=COLDEOS;
  ptr_pressure_rho0_u[whicheos] = &pressure_rho0_u_coldgrmhd;
  ptr_compute_u_from_entropy[whicheos] = &compute_u_from_entropy_coldgrmhd;
  ptr_u_rho0_T[whicheos] = &u_rho0_T_coldgrmhd;
  ptr_dpdu_rho0_u[whicheos] = &dpdu_rho0_u_coldgrmhd;
  ptr_dpdrho0_rho0_u[whicheos] = &dpdrho0_rho0_u_coldgrmhd;
  ptr_cs2_compute[whicheos] = &cs2_compute_coldgrmhd;

  ptr_compute_dSdrho[whicheos] = &compute_dSdrho_coldgrmhd;
  ptr_compute_dSdu[whicheos] = &compute_dSdu_coldgrmhd;
  ptr_compute_entropy[whicheos] = &compute_entropy_coldgrmhd;

  ptr_compute_dspecificSdrho_wmrho0[whicheos] = &compute_dspecificSdrho_wmrho0_coldgrmhd;
  ptr_compute_dspecificSdwmrho0_wmrho0[whicheos] = &compute_dspecificSdwmrho0_wmrho0_coldgrmhd;
  ptr_compute_specificentropy_wmrho0[whicheos] = &compute_specificentropy_wmrho0_coldgrmhd;
            
  ptr_pressure_wmrho0[whicheos] = &pressure_wmrho0_coldgrmhd;
  ptr_compute_idwmrho0dp[whicheos] = &compute_idwmrho0dp_coldgrmhd;
  ptr_compute_idrho0dp[whicheos] = &compute_idrho0dp_coldgrmhd;

  ptr_compute_qdot[whicheos] = &compute_qdot_coldgrmhd;
  ptr_compute_sources_EOS[whicheos] = &compute_sources_EOS_coldgrmhd;
  ptr_compute_allextras[whicheos] = &compute_allextras_coldgrmhd;
  ptr_get_extrasprocessed[whicheos] = &get_extrasprocessed_coldgrmhd;


  ptr_compute_temp[whicheos] = &compute_temp_coldgrmhd;

  ptr_compute_EOS_parms[whicheos] = &compute_EOS_parms_coldgrmhd;
  ptr_compute_EOS_parms_full[whicheos] = &compute_EOS_parms_coldgrmhd; // same as non-full
  ptr_store_EOS_parms[whicheos] = &store_EOS_parms_coldgrmhd;
  ptr_get_EOS_parms[whicheos] = &get_EOS_parms_coldgrmhd;
  ptr_fix_primitive_eos_scalars[whicheos] = &fix_primitive_eos_scalars_coldgrmhd;
  ptr_getall_forinversion[whicheos] = &getall_forinversion_coldgrmhd;


  //////////////////////////////////////////////////////



  return(0);

}





