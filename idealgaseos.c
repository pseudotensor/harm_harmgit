// IDEAL GAS EOS

// P = (\GAMMA -1) u
// h = 1+\Gamma_r \Theta  ; \Theta=p/\rho_0  \Gamma_r = (\GAMMA/(\GAMMA-1))

#define GAMMA (gamideal)

// 1/\Gamma_r
#define GAMMAM1 (GAMMA-1.0)
#define IGAMMAR (GAMMAM1/GAMMA)


// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
#if(OLDCALC)
  return(GAMMAM1*u) ;
#else
  return((GAMMA - 1.)*u) ;
#endif

}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
	return(p/GAMMAM1) ;
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
	return(GAMMAM1) ;
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
	return(0.0) ;
}

// sound speed squared (for vchar.c)
FTYPE cs2_compute_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u(EOSextra,rho0,u);
  h=rho0+u+pressure;

  cs2=GAMMA*pressure/h;

  return(cs2);

}

// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
// entropy = \rho\ln( p^n/\rho^{n+1} )
FTYPE compute_entropy_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
  FTYPE pressure,indexn,entropy;


  pressure=pressure_rho0_u(EOSextra,rho0,u);
  indexn=1.0/GAMMAM1;

  if(rho0<SMALL) rho0=SMALL;
  if(pressure<SMALL) pressure=SMALL;
  
  entropy=rho0*log(pow(pressure,indexn)/pow(rho0,indexn+1.0));

  // OPTMARK: Note that above and below ideal gas pow() shows up as expensive in pfmon (didn't appear in gprof!)
  // 2nd most expensive function is these two pow() calls after get_state_uconucovonly()
  // OPTMARK: log is bit less of a concern apparently, even though shows up in pfmon it's small %.

  return(entropy);

}

// u(rho0,S)
FTYPE compute_u_from_entropy_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  FTYPE rho,ie,pressure;
  FTYPE indexn;
  FTYPE u;

  indexn=1.0/GAMMAM1;

  if(rho0<SMALL) rho0=SMALL;
  
  // entropy version of ie
  u=pow(pow(rho0,indexn+1.0)*exp(entropy/rho0),1.0/indexn)/GAMMAM1;


  return(u);

}

// used for dudp_calc
FTYPE compute_dSdrho_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE entropy;
  FTYPE dSdrho;
  FTYPE compute_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE u);

  entropy=compute_entropy(EOSextra,rho0,u);

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdrho=entropy/rho0-(indexn+1.0);

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE dSdu;

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdu=indexn*rho0/u;

  return(dSdu);

}



// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(IGAMMAR*wmrho0) ;
}

// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(GAMMAM1/GAMMA);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}


// cooling or heating rate
FTYPE compute_qdot_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0);
}

int compute_sources_EOS_idealgas(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  int sc,pl,pliter;

  SCPHYSICSLOOP(sc) PLOOP(pliter,pl) dUcomp[sc][pl]=0.0;
  return(0);
}


void compute_allextras_idealgas(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
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

int get_extrasprocessed_idealgas(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  int pi;
  int numextrasreturn;
  FTYPE rho0,u;
  rho0=pr[RHO];
  u=pr[UU];
  compute_allextras(0, EOSextra, rho0, u, &numextrasreturn, extras);
  for(pi=0;pi<MAXPROCESSEDEXTRAS;pi++){
    processed[pi] = 0.0;
  }
  return(0);
}



// temperature
FTYPE compute_temp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE temp;

  temp = GAMMAM1*u/rho0;

  return(temp);

}

void compute_EOS_parms_idealgas(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  return; // do nothing
}



void store_EOS_parms_idealgas(int numparms, FTYPE *EOSextra, FTYPE *parlist)
{
  return; // do nothing
}


void get_EOS_parms_idealgas(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
{

  *numparms=0;
  return; // do nothing
}

#undef GAMMA
#undef GAMMAM1
#undef IGAMMAR
