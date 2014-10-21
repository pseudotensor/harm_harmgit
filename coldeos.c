//////////////////////////////////////////////////////
//
// COLD EOS or COLD EOMTYPE==EOMCOLDGRMHD
//
// Also used for EOMTYPE==EOMFFDE || EOMTYPE==EOMFFDE2
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

// u(rho0, T) (used for initial conditions)
FTYPE u_rho0_T_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE T)
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




FTYPE compute_specificentropy_wmrho0_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);

}

FTYPE compute_dspecificSdrho_wmrho0_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);

}


FTYPE compute_dspecificSdwmrho0_wmrho0_coldgrmhd(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{

  return(0.0);

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

void fix_primitive_eos_scalars_coldgrmhd(FTYPE *EOSextra, FTYPE *pr)
{
  return; // do nothing
}

void getall_forinversion_coldgrmhd(int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{

  *fun=0.0;
  *dfunofrho=0.0;
  *dfunofu=0.0;

}
