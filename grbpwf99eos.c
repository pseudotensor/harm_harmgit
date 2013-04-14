// EOS from PWF99 (GODMARK: mignone really, not done yet)



// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{

  return(pressure_rho0_u_idealgas(EOSextra,rho0,u));
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return(u_rho0_p_idealgas(EOSextra,rho0,p));
}

// u(rho0, T) (used for initial conditions)
FTYPE u_rho0_T_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE T)
{
  return(u_rho0_T_idealgas(EOSextra,rho0,T));
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(dpdu_rho0_u_idealgas(EOSextra,rho0,u));
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(dpdrho0_rho0_u_idealgas(EOSextra,rho0,u));
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(cs2_compute_idealgas(EOSextra,rho0,u));
}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(compute_entropy_idealgas(EOSextra,rho0,u));
}


// u(rho0,S)
FTYPE compute_u_from_entropy_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  return(compute_u_from_entropy_idealgas(EOSextra,rho0,entropy));
}



// used for dudp_calc
FTYPE compute_dSdrho_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(compute_dSdrho_idealgas(EOSextra,rho0,u));
}


// used for dudp_calc
FTYPE compute_dSdu_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(compute_dSdu_idealgas(EOSextra,rho0,u));
}

FTYPE compute_specificentropy_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(compute_specificentropy_wmrho0_idealgas(EOSextra,rho0,wmrho0));
}



FTYPE compute_dspecificSdrho_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(compute_dspecificSdrho_wmrho0_idealgas(EOSextra,  rho0,  wmrho0));
}

FTYPE compute_dspecificSdwmrho0_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(compute_dspecificSdwmrho0_wmrho0_idealgas(EOSextra,  rho0,  wmrho0));
}



// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(pressure_wmrho0_idealgas(EOSextra,  rho0,  wmrho0));
}



// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(compute_idwmrho0dp_idealgas(EOSextra,  rho0,  wmrho0));
}

// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(compute_idrho0dp_idealgas(EOSextra,  rho0,  wmrho0));
}


// cooling or heating rate
FTYPE compute_qdot_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)        
{
  return(0.0);
}

int compute_sources_EOS_grbpwf99(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  int sc,pl,pliter;

  SCPHYSICSLOOP(sc) PLOOP(pliter,pl) dUcomp[sc][pl]=0.0;
  return(0);
}


void compute_allextras_grbpwf99(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  int i;

  if(justnum==0){

    // set rest to 0
    for(i=0;i<MAXNUMEXTRAS;i++){
      extras[i] = 0.0;
    }
  }

  *numextrasreturn=0;
}

int get_extrasprocessed_grbpwf99(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  int pi;
  int numextrasreturn;
  FTYPE rho0,u;
  rho0=pr[RHO];
  u=pr[UU];
  compute_allextras_grbpwf99(0, EOSextra, rho0, u, &numextrasreturn, extras);
  for(pi=0;pi<MAXPROCESSEDEXTRAS;pi++){
    processed[pi] = 0.0;
  }
  return(0);
}




// temperature
FTYPE compute_temp_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE temp;
  FTYPE pressure_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u);

  temp = pressure_rho0_u_grbpwf99(EOSextra,rho0,u)/rho0;
  if(temp<SMALL) temp=SMALL;

  return(temp);

}


void compute_EOS_parms_grbpwf99(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS],  FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  return; // do nothing                                                                                           
}


void store_EOS_parms_grbpwf99(int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  return; // do nothing
}


void get_EOS_parms_grbpwf99(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{

  *numparms=0;
  return; // do nothing
}

void fix_primitive_eos_scalars_grbpwf99(FTYPE *EOSextra, FTYPE *pr)
{
  return; // do nothing
}

void getall_forinversion_grbpwf99(int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{
  getall_forinversion_idealgas(eomtype,whichd,EOSextra,quant1,quant2,fun,dfunofrho,dfunofu);
  return;
}
