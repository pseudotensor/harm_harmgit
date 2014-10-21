// TM EOS

// See entropy_new_inversion_mignone.nb

static FTYPE compute_inside_entropy_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
static FTYPE u_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);

// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;

  pressure = u*(2.0*rho0+u)/(3.0*(rho0+u));

  return(pressure);
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return( 1.5*(p + 3.0*p*p/(2.0*rho0+sqrt(9.0*p*p+4.0*rho0*rho0))) );
}

// u(rho0, T) (used for initial conditions)
FTYPE u_rho0_T_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE T)
{
  // solve from u(rho0,p) and T(rho0,p)
  return(0.5*(-2.*rho0 + 3.*rho0*T + sqrt(4.*rho0*rho0 + 9.*rho0*rho0*T*T)));
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dpdu;

  dpdu = 1.0/3.0*(1.0 + rho0*rho0/( (rho0+u)*(rho0+u)));

  return(dpdu);
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dpdrho0;

  dpdrho0 = u*u/(3.0*(rho0+u)*(rho0+u));

  return(dpdrho0) ;
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u_mignone(EOSextra,rho0,u);
  h=rho0+u+pressure; // not specific h

  cs2=pressure*(5.0*h - 8.0*pressure) / (3.0*h*(h-pressure));


  //  dualfprintf(fail_file,"cs2=%21.15g pressure=%21.15g\n",cs2,pressure);

  return(cs2);

}

// entropy density (per unit volume) as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE entropy,insideentropy;

  insideentropy = compute_inside_entropy_mignone(EOSextra,rho0,u);

  entropy = rho0*log(insideentropy);

  return(entropy);

}

//local aux function
static FTYPE compute_inside_entropy_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure,insideentropy;
  FTYPE pressure_rho0_u_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  

  pressure=pressure_rho0_u_mignone(EOSextra,rho0,u);

  // Entropy will be wrong when rho0 or pressure are non-positive
  if(rho0<SMALL) rho0=SMALL;
  if(u<SMALL) u=SMALL;
  if(pressure<SMALL) pressure=SMALL;
  
  insideentropy=pressure*(rho0+u)/pow(rho0,8.0/3.0);

  return(insideentropy);

}


// u(rho0,Sden)
FTYPE compute_u_from_entropy_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  FTYPE u;
  FTYPE expfactor;

  // u will be wrong when density is non-positive
  if(rho0<SMALL) rho0=SMALL;

  expfactor = exp(entropy/rho0)*pow(rho0,8.0/3.0);
  u = 3.0*expfactor/(rho0+sqrt(rho0*rho0+3.0*expfactor));

  return(u);

}



// used for dudp_calc
FTYPE compute_dSdrho_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dSdrho;

  // Limiting the below values will likely confuse utoprim.orig.c inversion if reached
  if(rho0<SMALL) rho0=SMALL;
  if(u<SMALL) u=SMALL;

  dSdrho = (-8.0/3.0) + 2.0*rho0/(2.0*rho0+u) + log(u*(2.0*rho0+u)/(3.0*pow(rho0,8.0/3.0)));

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dSdu;

  // Limiting the below values will likely confuse utoprim.orig.c inversion if reached
  if(rho0<SMALL) rho0=SMALL;
  if(u<SMALL) u=SMALL;

  dSdu = rho0*(1.0/u + 1.0/(2.0*rho0+u));

  return(dSdu);

}

FTYPE compute_entropy_wmrho0_mignone_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE u,P,entropy;
  FTYPE pressure_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);

  u=u_wmrho0_mignone(EOSextra, rho0, wmrho0);
  P=pressure_wmrho0_mignone(EOSextra, rho0, wmrho0);

  // Don't limit rho0 and pressure since this is used for iterative scheme that requires to know if beyond valid domain or not.  Nan will be terminated during inversion.
  //  if(rho0<SMALL) rho0=SMALL;
  //  if(u<SMALL) u=SMALL;
  //  if(P<SMALL) P=SMALL;

  entropy = rho0*log(P*(rho0+u)/pow(rho0,8.0/3.0));

  return(entropy);

}

// used for utoprim_jon.c entropy inversion
FTYPE compute_specificentropy_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE u,P,specificentropy;
  FTYPE pressure_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);

  u=u_wmrho0_mignone(EOSextra, rho0, wmrho0);
  P=pressure_wmrho0_mignone(EOSextra, rho0, wmrho0);

  // Don't limit rho0 and pressure since this is used for iterative scheme that requires to know if beyond valid domain or not.  Nan will be terminated during inversion.
  //  if(rho0<SMALL) rho0=SMALL;
  //  if(u<SMALL) u=SMALL;
  //  if(P<SMALL) P=SMALL;

  specificentropy = log(P*(rho0+u)/pow(rho0,8.0/3.0));

  return(specificentropy);

}




FTYPE compute_dSdrho_wmrho0_mignone_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE compute_entropy_wmrho0_mignone_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  FTYPE dSdrho;
  FTYPE entropy;

  entropy=compute_entropy_wmrho0_mignone_unused(EOSextra, rho0, wmrho0);

  dSdrho = (-5.0/3.0) + entropy/rho0 + rho0/(2.0*rho0+wmrho0) + (-3.0*wmrho0*wmrho0-6.0*wmrho0*rho0-5.0*rho0*rho0)/((2.0*rho0+wmrho0)*sqrt(9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+25.0*rho0*rho0));

  return(dSdrho);
}

// used for utoprim_jon.c entropy inversion
FTYPE compute_dspecificSdrho_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE sqrtthing,sqrtinside,dSdrho;

  sqrtinside=9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+25.0*rho0*rho0;
  sqrtthing=sqrt(sqrtinside);
  dSdrho = -(9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+15.0*rho0*rho0+5.0*wmrho0*sqrtthing + 7.0*rho0*sqrtthing)/(3.0*rho0*(wmrho0+2.0*rho0)*sqrtthing);

  return(dSdrho);
}

FTYPE compute_dSdwmrho0_wmrho0_mignone_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;
  FTYPE sqrtthing;

  sqrtthing = sqrt(9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+25.0*rho0*rho0);
  dSdchi = (rho0*(3.0*wmrho0*wmrho0+rho0*(5.0*rho0+sqrtthing)+wmrho0*(6.0*rho0+sqrtthing)))/(wmrho0*(2.0*rho0+wmrho0)*sqrtthing);

  return(dSdchi);

}

// used for utoprim_jon.c entropy inversion
FTYPE compute_dspecificSdwmrho0_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;
  FTYPE sqrtthing,sqrtinside;

  sqrtinside=9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+25.0*rho0*rho0;
  sqrtthing=sqrt(sqrtinside);

  dSdchi = (wmrho0+rho0-10.0*rho0*rho0/(3.0*sqrtthing)+sqrtthing/3.0)/(wmrho0*(wmrho0+2.0*rho0));

  return(dSdchi);

}






// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE QQ,delta,delta2;
  FTYPE pressure;

  QQ=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+QQ);
  delta2=delta/rho0;

  pressure=(5.0/8.0)*(wmrho0 - delta/(1.0+sqrt(1.0+delta2)));

  return(pressure);
}

// u(rho0, w-rho0 = u+p)
static FTYPE u_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE u;
  FTYPE sqrtthing;

  sqrtthing=sqrt(9.0*wmrho0*wmrho0+18.0*wmrho0*rho0+25.0*rho0*rho0);

  u = 3.0*wmrho0*(3.0*wmrho0+3.0*rho0+sqrtthing)/(4.0*(3.0*wmrho0+5.0*rho0+sqrtthing));

  return(u);
}


// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_mignone_old(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE QQ,delta,delta2;
  FTYPE ddeltadwmrho0,idwmrho0dp;

  QQ=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+QQ);
  delta2=delta/rho0;
  
  ddeltadwmrho0=18.0/25.0*(1.0+QQ);

  idwmrho0dp = 5.0/16.0*(2.0-ddeltadwmrho0/sqrt(1.0+delta2));

  return(idwmrho0dp);

}

// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE p;

  p = pressure_wmrho0_mignone(EOSextra, rho0, wmrho0);

  idwmrho0dp = (2.0*wmrho0 + 2.0*rho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idwmrho0dp);

}

// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  FTYPE idrho0dp;
  FTYPE p;

  p = pressure_wmrho0_mignone(EOSextra, rho0, wmrho0);

  idrho0dp = (2.0*wmrho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idrho0dp);
}


// cooling or heating rate
FTYPE compute_qdot_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)        
{
  return(0.0);
}

int compute_sources_EOS_mignone(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  int sc,pl,pliter;

  SCPHYSICSLOOP(sc) PLOOP(pliter,pl) dUcomp[sc][pl]=0.0;
  return(0);
}


void compute_allextras_mignone(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
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

int get_extrasprocessed_mignone(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  int pi;
  int numextrasreturn;
  FTYPE rho0,u;
  rho0=pr[RHO];
  u=pr[UU];
  compute_allextras_mignone(0, EOSextra, rho0, u, &numextrasreturn, extras);
  for(pi=0;pi<MAXPROCESSEDEXTRAS;pi++){
    processed[pi] = 0.0;
  }
  return(0);
}




// temperature
FTYPE compute_temp_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE temp;
  FTYPE pressure_rho0_u_mignone(FTYPE *EOSextra, FTYPE rho0, FTYPE u);

  temp = pressure_rho0_u_mignone(EOSextra,rho0,u)/rho0;
  if(temp<SMALL) temp=SMALL;

  return(temp);

}


void compute_EOS_parms_mignone(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS],  FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  return; // do nothing                                                                                           
}


void store_EOS_parms_mignone(int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  return; // do nothing
}


void get_EOS_parms_mignone(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{

  *numparms=0;
  return; // do nothing
}

void fix_primitive_eos_scalars_mignone(FTYPE *EOSextra, FTYPE *pr)
{
  return; // do nothing
}


// this could be optimized since redundant calculations
// worry about missed branch predictions?
void getall_forinversion_mignone(int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{

  if(eomtype==EOMGRMHD){
    if(whichd==CHIDIFF){
      // PofRHOCHI
      *fun=pressure_wmrho0_mignone(EOSextra, quant1, quant2);
      *dfunofrho=compute_idrho0dp_mignone(EOSextra, quant1, quant2);
      *dfunofu=compute_idwmrho0dp_mignone(EOSextra, quant1, quant2);
    }
    else if(whichd==UTOTDIFF){
      *fun=pressure_rho0_u_mignone(EOSextra, quant1, quant2);
      *dfunofrho=dpdrho0_rho0_u_mignone(EOSextra, quant1, quant2);
      *dfunofu=dpdu_rho0_u_mignone(EOSextra, quant1, quant2);
    }
  }
  else if(eomtype==EOMENTROPYGRMHD){
    // NOTE: for \chi version it's specific entropy and for u version it's entropy density
    if(whichd==CHIDIFF){
      *fun=compute_specificentropy_wmrho0_mignone(EOSextra, quant1, quant2);
      *dfunofrho=compute_dspecificSdrho_wmrho0_mignone(EOSextra, quant1, quant2);
      *dfunofu=compute_dspecificSdwmrho0_wmrho0_mignone(EOSextra, quant1, quant2);
    }
    else if(whichd==UTOTDIFF){
      *fun=compute_entropy_mignone(EOSextra, quant1, quant2);
      *dfunofrho=compute_dSdrho_mignone(EOSextra, quant1, quant2);
      *dfunofu=compute_dSdu_mignone(EOSextra, quant1, quant2);
    }
  }
  else if(eomtype==EOMCOLDGRMHD || eomtype==EOMFFDE || eomtype==EOMFFDE2){
    *fun=*dfunofrho=*dfunofu=0.0;
  }

  return; // done!
}
