#include "decs.h"



////////////////////////////
//
// IMPLEMENT EOS
//
// u_rho0_p : used by initial conditions
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
#include "grbpwf99eos.c"
#include "idealgaseos.c"
#include "mignoneeos.c"
#include "kazfulleos.c"








//////////////////////////////////////////////////////
//
// COLD EOS or COLD EOMTYPE
//
// then explicitly set pressure/ie to 0 so can keep some parts of code similar
//
//////////////////////////////////////////////////////

// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return(0.0) ;
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// used for dudp_calc
FTYPE compute_dSdrho_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// used for dudp_calc
FTYPE compute_dSdu_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

// u(rho0,S)
FTYPE compute_u_from_entropy_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  return(0.0) ;
}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0) ;
}


// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}


FTYPE compute_qdot_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0);
}

int compute_sources_EOS_coldgrmhd(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  return(0);
}


void compute_allextras_coldgrmhd(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  return;
}

int get_extrasprocessed_coldgrmhd(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  return(0);
}


FTYPE compute_temp_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0);
}

void compute_EOS_parms_coldgrmhd(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  return; // do nothing
}

void store_EOS_parms_coldgrmhd(int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  return; // do nothing
}
void get_EOS_parms_coldgrmhd(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  return; // do nothing
}





//////////////////////////////////////////////////
//
// wrappers for any EOS and any EOMTYPE
//
//////////////////////////////////////////////////

// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_pressure_rho0_u)(EOSextra,rho0,u) );
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return( (*ptr_u_rho0_p)(EOSextra,rho0,p) );
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_dpdu_rho0_u)(EOSextra,rho0,u) );
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_dpdrho0_rho0_u)(EOSextra,rho0,u) );
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_cs2_compute)(EOSextra,rho0,u) );
}


// used for dudp_calc
FTYPE compute_dSdrho(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_dSdrho)(EOSextra,rho0,u) );
}


// used for dudp_calc
FTYPE compute_dSdu(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_dSdu)(EOSextra,rho0,u) );
}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_entropy)(EOSextra,rho0,u) );
}

// u(rho0,S)
FTYPE compute_u_from_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  return( (*ptr_compute_u_from_entropy)(EOSextra,rho0,entropy) );
}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_pressure_wmrho0)(EOSextra,rho0,wmrho0) );
}


// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_compute_idwmrho0dp)(EOSextra,rho0,wmrho0) );
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return( (*ptr_compute_idrho0dp) (EOSextra,rho0,wmrho0) );
}


// radiation rate
FTYPE compute_qdot(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_qdot) (EOSextra,rho0,u) );
}


int compute_sources_EOS(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  return( (*ptr_compute_sources_EOS) (EOSextra,pr, geom, q, Ui, dUother, dUcomp));
}



void compute_allextras(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  (*ptr_compute_allextras) (justnum,EOSextra,rho0,u,numextrasreturn,extras);
  return;
}

int get_extrasprocessed(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  return((*ptr_get_extrasprocessed) (doall,EOSextra,pr,extras,processed));
}



FTYPE compute_temp(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return( (*ptr_compute_temp) (EOSextra,rho0,u) );
}


void compute_EOS_parms(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  (*ptr_compute_EOS_parms) (EOSextra,prim);
}

void store_EOS_parms(int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  (*ptr_store_EOS_parms) (numparms,EOSextra,parlist);
}
void get_EOS_parms(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  (*ptr_get_EOS_parms) (numparms, EOSextra, parlist);
}




//////////////////////////////////////////
//
// ANY EOS
//

// old function
// p(rho0,w)
FTYPE pressure_rho0_w(FTYPE *EOSextra, FTYPE rho0, FTYPE w)
{
  FTYPE wmrho0=w-rho0;

#if(OLDCALC)
  return((GAMMA-1.)*(w - rho0)/GAMMA) ;
#else
  return(pressure_wmrho0(EOSextra,rho0,wmrho0)) ;
#endif
}


// pick EOMTYPE for EOS
int pickeos_eomtype(int whicheos, int whicheom)
{
  //IDEALGAS
  //MIGNONE
  //GRBPWF99
  //KAZFULL


  if(whicheos==IDEALGAS){
    // then need to set gamideal
    gamideal=gam;
  }
  else{
    // otherwise assume gam and gamideal could be different, where gam is used to true EOS for some purpose while gamideal is particular always for ideal EOS
  }


  // EOS functions used during inversion and other places
  if(whicheom==EOMGRMHD){
    if(whicheos==IDEALGAS){
      ptr_pressure_rho0_u = &pressure_rho0_u_idealgas;
      ptr_compute_u_from_entropy = &compute_u_from_entropy_idealgas;
      ptr_u_rho0_p = &u_rho0_p_idealgas;
      ptr_dpdu_rho0_u = &dpdu_rho0_u_idealgas;
      ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_idealgas;
      ptr_cs2_compute = &cs2_compute_idealgas;
      ptr_compute_dSdrho = &compute_dSdrho_idealgas;
      ptr_compute_dSdu = &compute_dSdu_idealgas;
      ptr_compute_entropy = &compute_entropy_idealgas;
      
      ptr_pressure_wmrho0=&pressure_wmrho0_idealgas;
      ptr_compute_idwmrho0dp=&compute_idwmrho0dp_idealgas;
      ptr_compute_idrho0dp=&compute_idrho0dp_idealgas;
      
      ptr_compute_qdot=&compute_qdot_idealgas;
      ptr_compute_sources_EOS=&compute_sources_EOS_idealgas;
      ptr_compute_allextras=&compute_allextras_idealgas;
      ptr_get_extrasprocessed=&get_extrasprocessed_idealgas;
      
      ptr_compute_temp=&compute_temp_idealgas;
      
      ptr_compute_EOS_parms=&compute_EOS_parms_idealgas;
      ptr_store_EOS_parms=&store_EOS_parms_idealgas;
      ptr_get_EOS_parms=&get_EOS_parms_idealgas;
    }
    else if(whicheos==MIGNONE){
      ptr_pressure_rho0_u = &pressure_rho0_u_mignone;
      ptr_compute_u_from_entropy = &compute_u_from_entropy_mignone;
      ptr_u_rho0_p = &u_rho0_p_mignone;
      ptr_dpdu_rho0_u = &dpdu_rho0_u_mignone;
      ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_mignone;
      ptr_cs2_compute = &cs2_compute_mignone;
      ptr_compute_dSdrho = &compute_dSdrho_mignone;
      ptr_compute_dSdu = &compute_dSdu_mignone;
      ptr_compute_entropy = &compute_entropy_mignone;
      
      ptr_pressure_wmrho0=&pressure_wmrho0_mignone;
      ptr_compute_idwmrho0dp=&compute_idwmrho0dp_mignone;
      ptr_compute_idrho0dp=&compute_idrho0dp_mignone;
      
      ptr_compute_qdot=&compute_qdot_mignone;
      ptr_compute_sources_EOS=&compute_sources_EOS_mignone;
      ptr_compute_allextras=&compute_allextras_mignone;
      ptr_get_extrasprocessed=&get_extrasprocessed_mignone;


      ptr_compute_temp=&compute_temp_mignone;
      
      ptr_compute_EOS_parms=&compute_EOS_parms_mignone;
      ptr_store_EOS_parms=&store_EOS_parms_mignone;
      ptr_get_EOS_parms=&get_EOS_parms_mignone;
    }
    else if(whicheos==GRBPWF99){
      ptr_pressure_rho0_u = &pressure_rho0_u_grbpwf99;
      ptr_compute_u_from_entropy = &compute_u_from_entropy_grbpwf99;
      ptr_u_rho0_p = &u_rho0_p_grbpwf99;
      ptr_dpdu_rho0_u = &dpdu_rho0_u_grbpwf99;
      ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_grbpwf99;
      ptr_cs2_compute = &cs2_compute_grbpwf99;
      ptr_compute_dSdrho = &compute_dSdrho_grbpwf99;
      ptr_compute_dSdu = &compute_dSdu_grbpwf99;
      ptr_compute_entropy = &compute_entropy_grbpwf99;
      
      ptr_pressure_wmrho0=&pressure_wmrho0_grbpwf99;
      ptr_compute_idwmrho0dp=&compute_idwmrho0dp_grbpwf99;
      ptr_compute_idrho0dp=&compute_idrho0dp_grbpwf99;
      
      ptr_compute_qdot=&compute_qdot_grbpwf99;
      ptr_compute_sources_EOS=&compute_sources_EOS_grbpwf99;
      ptr_compute_allextras=&compute_allextras_grbpwf99;
      ptr_get_extrasprocessed=&get_extrasprocessed_grbpwf99;

      
      ptr_compute_temp=&compute_temp_grbpwf99;
      
      ptr_compute_EOS_parms=&compute_EOS_parms_grbpwf99;
      ptr_store_EOS_parms=&store_EOS_parms_grbpwf99;
      ptr_get_EOS_parms=&get_EOS_parms_grbpwf99;
    }
    else if(whicheos==KAZFULL){
      ptr_pressure_rho0_u = &pressure_rho0_u_kazfull;
      ptr_compute_u_from_entropy = &compute_u_from_entropy_kazfull;
      ptr_u_rho0_p = &u_rho0_p_kazfull;
      ptr_dpdu_rho0_u = &dpdu_rho0_u_kazfull;
      ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_kazfull;
      ptr_cs2_compute = &cs2_compute_kazfull;
      ptr_compute_dSdrho = &compute_dSdrho_kazfull;
      ptr_compute_dSdu = &compute_dSdu_kazfull;
      ptr_compute_entropy = &compute_entropy_kazfull;
      
      ptr_pressure_wmrho0=&pressure_wmrho0_kazfull;
      ptr_compute_idwmrho0dp=&compute_idwmrho0dp_kazfull;
      ptr_compute_idrho0dp=&compute_idrho0dp_kazfull;
      
      ptr_compute_qdot=&compute_qdot_kazfull;
      ptr_compute_sources_EOS=&compute_sources_EOS_kazfull;
      ptr_compute_allextras=&compute_allextras_kazfull;
      ptr_get_extrasprocessed=&get_extrasprocessed_kazfull;

      
      ptr_compute_temp=&compute_temp_kazfull;
      
      ptr_compute_EOS_parms=&compute_EOS_parms_kazfull;
      ptr_store_EOS_parms=&store_EOS_parms_kazfull;
      ptr_get_EOS_parms=&get_EOS_parms_kazfull;
    }

  }
  else{
    ptr_pressure_rho0_u = &pressure_rho0_u_coldgrmhd;
    ptr_compute_u_from_entropy = &compute_u_from_entropy_coldgrmhd;
    ptr_u_rho0_p = &u_rho0_p_coldgrmhd;
    ptr_dpdu_rho0_u = &dpdu_rho0_u_coldgrmhd;
    ptr_dpdrho0_rho0_u = &dpdrho0_rho0_u_coldgrmhd;
    ptr_cs2_compute = &cs2_compute_coldgrmhd;
    ptr_compute_dSdrho = &compute_dSdrho_coldgrmhd;
    ptr_compute_dSdu = &compute_dSdu_coldgrmhd;
    ptr_compute_entropy = &compute_entropy_coldgrmhd;

    ptr_pressure_wmrho0=&pressure_wmrho0_coldgrmhd;
    ptr_compute_idwmrho0dp=&compute_idwmrho0dp_coldgrmhd;
    ptr_compute_idrho0dp=&compute_idrho0dp_coldgrmhd;

    ptr_compute_qdot=&compute_qdot_coldgrmhd;
    ptr_compute_sources_EOS=&compute_sources_EOS_coldgrmhd;
    ptr_compute_allextras=&compute_allextras_coldgrmhd;
    ptr_get_extrasprocessed=&get_extrasprocessed_coldgrmhd;


    ptr_compute_temp=&compute_temp_coldgrmhd;

    ptr_compute_EOS_parms=&compute_EOS_parms_coldgrmhd;
    ptr_store_EOS_parms=&store_EOS_parms_coldgrmhd;
    ptr_get_EOS_parms=&get_EOS_parms_coldgrmhd;

  }

  return(0);

}


