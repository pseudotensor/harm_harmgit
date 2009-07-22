// EOS from PWF99 (GODMARK: mignone really, not done yet)



// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;

  pressure = u*(2.0*rho0+u)/(3.0*(rho0+u));

  return(pressure);
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return( 1.5*(p + 3.0*p*p/(2.0*rho0+sqrt(9.0*p*p+4.0*rho0*rho0))) );
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dpdu;

  dpdu = 1.0/3.0*(1.0 + rho0*rho0/( (rho0+u)*(rho0+u)));

  return(dpdu);
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dpdrho0;

  dpdrho0 = u*u/(3.0*(rho0+u)*(rho0+u));

  return(dpdrho0) ;
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u_grbpwf99(EOSextra,rho0,u);
  h=rho0+u+pressure; // not specific h

  cs2=pressure*(5.0*h - 8.0*pressure) / (3.0*h*(h-pressure));


  //  dualfprintf(fail_file,"cs2=%21.15g pressure=%21.15g\n",cs2,pressure);

  return(cs2);

}


// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE entropy;

  // not setup yet for TM EOS
  entropy=compute_entropy_idealgas(EOSextra,rho0,u);

  return(entropy);

}


// u(rho0,S)
FTYPE compute_u_from_entropy_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  FTYPE u;

  // not setup yet for TM EOS
  u = compute_u_from_entropy_idealgas(EOSextra,rho0,entropy);

  return(u);

}



// used for dudp_calc
FTYPE compute_dSdrho_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dSdrho;

  // not setup yet for TM EOS
  dSdrho=compute_dSdrho_idealgas(EOSextra,rho0,u);

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dSdu;

  // not setup yet for TM EOS
  dSdu = compute_dSdu_idealgas(EOSextra,rho0,u);

  return(dSdu);

}

FTYPE compute_specificentropy_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE specificentropy;

  // not setup yet for TM EOS
  specificentropy=compute_specificentropy_wmrho0_idealgas(EOSextra,rho0,wmrho0);

  return(specificentropy);

}



FTYPE compute_dspecificSdrho_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;

  // not setup yet for TM EOS
  dSdchi=compute_dspecificSdrho_wmrho0_idealgas(EOSextra,  rho0,  wmrho0);

  return(dSdchi);
}

FTYPE compute_dspecificSdwmrho0_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;

  // not setup yet for TM EOS
  dSdchi = compute_dspecificSdwmrho0_wmrho0_idealgas(EOSextra,  rho0,  wmrho0);

  return(dSdchi);

}



// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE pressure;

  Q=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+Q);
  delta2=delta/rho0;

  pressure=(5.0/8.0)*(wmrho0 - delta/(1.0+sqrt(1.0+delta2)));

  return(pressure);
}


// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_grbpwf99_old(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE ddeltadwmrho0,idwmrho0dp;

  Q=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+Q);
  delta2=delta/rho0;
  
  ddeltadwmrho0=18.0/25.0*(1.0+Q);

  idwmrho0dp = 5.0/16.0*(2.0-ddeltadwmrho0/sqrt(1.0+delta2));

  return(idwmrho0dp);

}

// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE p;

  p = pressure_wmrho0_grbpwf99(EOSextra, rho0, wmrho0);

  idwmrho0dp = (2.0*wmrho0 + 2.0*rho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idwmrho0dp);

}

// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0_grbpwf99(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  FTYPE idrho0dp;
  FTYPE p;

  p = pressure_wmrho0_grbpwf99(EOSextra, rho0, wmrho0);

  idrho0dp = (2.0*wmrho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idrho0dp);
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

