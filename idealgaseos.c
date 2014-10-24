/*! \file idealgaseos.c
  \brief IDEAL GAS EOS functions

  // P = (\GAMMA -1) u
  // h = 1+\Gamma_r \Theta  ; \Theta=p/\rho_0  \Gamma_r = (\GAMMA/(\GAMMA-1))
  // see entropy_new_inversion.nb
*/



#define GAMMA (gamideal)

// 1/\Gamma_r
#define GAMMAM1 (GAMMA-1.0)
#define IGAMMAR (GAMMAM1/GAMMA)


static FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);


/// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
#if(OLDCALC)
  return(GAMMAM1*u) ;
#else
  return((GAMMA - 1.)*u) ;
#endif

}

/// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  return(p/GAMMAM1) ;
}

// anywhere T or entropy are involved, really want number density, not mass density, so involves rho0/MUMEAN

/// u(rho0, T) (used for initial conditions)
FTYPE u_rho0_T_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE T)
{
  return(((rho0/MUMEAN)*T)/GAMMAM1) ;
}

/// dp(rho0, u)/du
FTYPE dpdu_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(GAMMAM1) ;
}

/// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(0.0) ;
}

/// sound speed squared (for vchar.c)
FTYPE cs2_compute_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u_idealgas(EOSextra,rho0,u);
  h=rho0+u+pressure+SMALL;

  cs2=GAMMA*pressure/h;
  if(cs2<SMALL) cs2=SMALL;

  return(cs2);

}

/// The resulting entropy is an entropy per unit volume.
/// The dimensionless specific entropy per baryon (or per gram really) would be this entropy divided by \rho_0
/// In general, entropy density often comes in units of erg/K/cc.  If one divides by k_b one gets 1/cc.  If one divides by n_b, one gets a dimensionless entropy per baryon.  If one divides by instead \rho_0=m_b n_b, one gets a entropy per baryon in per baryon rest-mass.
/// The prefactor constant doesn't matter for entropy evolution/inversion or sound speed evaluations as long as everything is self-consistent.
///
/// entropy (in energy/volume) as function of rho0 and internal energy (u)
/// Sden(rho0,u)
/// entropy density = \rho\ln( p^n/\rho^{n+1} )
FTYPE compute_entropy_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
  FTYPE pressure,indexn,entropy;


  pressure=pressure_rho0_u_idealgas(EOSextra,rho0,u);
  indexn=1.0/GAMMAM1;

  // Entropy will be wrong when rho0 or pressure are non-positive
  if(rho0<SMALL && pressure<SMALL*1E-10){
    rho0=SMALL;
    pressure=SMALL*1E-10;
  }
  else if(rho0<SMALL){
    rho0=SMALL;
  }
  else if(pressure<SMALL*1E-10){
    pressure=SMALL*1E-10;
    //    entropy=-BIG; // so not -inf
    // Too small entropy can cause precision/Newton inversion step issues, so implement floor on specific entropy
  }

  // normal:
  entropy=(rho0/MUMEAN)*log(pow(pressure,indexn)/pow(rho0/MUMEAN,indexn+1.0));

  // if p>>1 and rho0~0, then argument to log() could be out of range (i.e. inf).
  //  if(!isfinite(entropy)) entropy=0.0;

  // OPTMARK: Note that above and below ideal gas pow() shows up as expensive in pfmon (didn't appear in gprof!)
  // 2nd most expensive function is these two pow() calls after get_state_uconucovonly()
  // OPTMARK: log is bit less of a concern apparently, even though shows up in pfmon it's small %.

  return(entropy);

}

/// u(rho0,Sden = U[ENTROPY]/U[RHO]*pr[RHO])
FTYPE compute_u_from_entropy_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  FTYPE rho,ie,pressure;
  FTYPE indexn;
  FTYPE u;

  indexn=1.0/GAMMAM1;

  // u will be wrong when density is non-positive
  if(rho0<SMALL) rho0=SMALL;
  
  // entropy version of ie
  u=pow(pow(rho0/MUMEAN,indexn+1.0)*exp(entropy/(rho0/MUMEAN)),1.0/indexn)/GAMMAM1;


  return(u);

}



/// local aux function
static FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure,indexn,insideentropy;
  FTYPE pressure_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
  

  pressure=pressure_wmrho0_idealgas(EOSextra,rho0,wmrho0);
  indexn=1.0/GAMMAM1;

  // Don't limit rho0 and pressure since this is used for iterative scheme that requires to know if beyond valid domain or not.  Nan will be terminated during inversion.
  //  if(rho0<SMALL) rho0=SMALL;
  //  if(pressure<SMALL) pressure=SMALL;
  
  insideentropy=pow(pressure,indexn)/pow(rho0/MUMEAN,indexn+1.0);

  return(insideentropy);
}



/// used for dudp_calc
/// dSden/drho0
FTYPE compute_dSdrho_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE entropy;
  FTYPE dSdrho;
  FTYPE compute_entropy_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u);

  entropy=compute_entropy_idealgas(EOSextra,rho0,u);

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdrho=(entropy/(rho0/MUMEAN)-(indexn+1.0))/MUMEAN;

  return(dSdrho);

}


/// used for dudp_calc
/// dSden/du
FTYPE compute_dSdu_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE dSdu;

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdu=indexn*(rho0/MUMEAN)/u;

  return(dSdu);

}

/// entropy as function of rho0 and internal energy (u)
/// Sden(rho0,\chi=u+p)
/// entropy density = \rho\ln( p^n/\rho^{n+1} )
FTYPE compute_entropy_wmrho0_idealgas_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE insideentropy,entropy;

  insideentropy=compute_inside_entropy_wmrho0_idealgas(EOSextra, rho0, wmrho0);
  
  entropy=(rho0/MUMEAN)*log(insideentropy);

  //  if(!isfinite(entropy)) entropy=0.0;


  return(entropy);

}

/// specific entropy as function of rho0 and internal energy (u)
/// Ss(rho0,\chi=u+p)
/// specific entropy = \ln( p^n/\rho^{n+1} )
FTYPE compute_specificentropy_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE insideentropy,specificentropy;

  insideentropy=compute_inside_entropy_wmrho0_idealgas(EOSextra, rho0, wmrho0);
  
  specificentropy=log(insideentropy)/MUMEAN; // because Sspecific = S/rho alone by choice/definition.

  //  if(!isfinite(specificentropy)) specificentropy=0.0;


  return(specificentropy);

}

/// used for utoprim_jon when doing entropy evolution
/// Because P=(\gamma-1)u, then holding \chi=w-\rho_0 constant is the same as holding u constant
/// dSden/drho0
FTYPE compute_dSdrho_wmrho0_idealgas_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdrho;
  FTYPE insideentropy;


  insideentropy=compute_inside_entropy_wmrho0_idealgas(EOSextra, rho0, wmrho0);
  
  dSdrho=(GAMMA/(1.0-GAMMA) + log(insideentropy))/MUMEAN;

  //  if(!isfinite(dSdrho)) dSdrho=0.0;

  // Note that it makes no sense to speak of entropy changes with isothermal (GAMMA=1.0) gas since in the limit GAMMA->1, dSdrho->-\infty

  return(dSdrho);

}

/// dSspecific/drho0
FTYPE compute_dspecificSdrho_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdrho;

  // Don't limit rho0 and pressure since this is used for iterative scheme that requires to know if beyond valid domain or not.  Nan will be terminated during inversion.
  //  if(rho0<SMALL) rho0=SMALL;
  
  dSdrho=GAMMA/((1.0-GAMMA)*(rho0*MUMEAN)); // yes, really *MUMEAN because we define specific as S/rho and not S/(mu*rho)

  // Note that it makes no sense to speak of entropy changes with isothermal (GAMMA=1.0) gas since in the limit GAMMA->1, dSdrho->-\infty

  return(dSdrho);

}


/// used for utoprim_jon when doing entropy evolution
/// dSden/d\chi
FTYPE compute_dSdwmrho0_wmrho0_idealgas_unused(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;

  dSdchi = (rho0/MUMEAN)/(GAMMAM1*wmrho0);

  // Again, GAMMA->1 means dSdchi->\infty unless \chi->0 or rho0->0

  return(dSdchi);

}

/// used for utoprim_jon when doing entropy evolution
/// dSspecific/d\chi
FTYPE compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;

  dSdchi = 1.0/(GAMMAM1*wmrho0)/MUMEAN;

  // Again, GAMMA->1 means dSdchi->\infty unless \chi->0 or rho0->0

  return(dSdchi);

}



/// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(IGAMMAR*wmrho0) ;
}

/// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(GAMMAM1/GAMMA);
}


/// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(0.0);
}


/// cooling or heating rate
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
  compute_allextras_idealgas(0, EOSextra, rho0, u, &numextrasreturn, extras);
  for(pi=0;pi<MAXPROCESSEDEXTRAS;pi++){
    processed[pi] = 0.0;
  }
  return(0);
}



/// temperature
FTYPE compute_temp_idealgas(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE temp;

  temp = GAMMAM1*u/(SMALL+fabs(rho0/MUMEAN));
  // can't let temp go negative because have T^2 or T^4 terms in lambda or other functions
  if(temp<SMALL) temp=SMALL;

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

void fix_primitive_eos_scalars_idealgas(FTYPE *EOSextra, FTYPE *pr)
{
  return; // do nothing
}


/// this could be optimized since redundant calculations
/// worry about missed branch predictions?
void getall_forinversion_idealgas(int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{

  if(eomtype==EOMGRMHD){
    if(whichd==CHIDIFF){
      // PofRHOCHI
      *fun=pressure_wmrho0_idealgas(EOSextra, quant1, quant2);
      *dfunofrho=compute_idrho0dp_idealgas(EOSextra, quant1, quant2);
      *dfunofu=compute_idwmrho0dp_idealgas(EOSextra, quant1, quant2);
    }
    else if(whichd==UTOTDIFF){
      *fun=pressure_rho0_u_idealgas(EOSextra, quant1, quant2);
      *dfunofrho=dpdrho0_rho0_u_idealgas(EOSextra, quant1, quant2);
      *dfunofu=dpdu_rho0_u_idealgas(EOSextra, quant1, quant2);
    }
  }
  else if(eomtype==EOMENTROPYGRMHD){
    // NOTE: for \chi version it's specific entropy and for u version it's entropy density
    if(whichd==CHIDIFF){
      *fun=compute_specificentropy_wmrho0_idealgas(EOSextra, quant1, quant2);
      *dfunofrho=compute_dspecificSdrho_wmrho0_idealgas(EOSextra, quant1, quant2);
      *dfunofu=compute_dspecificSdwmrho0_wmrho0_idealgas(EOSextra, quant1, quant2);
    }
    else if(whichd==UTOTDIFF){
      *fun=compute_entropy_idealgas(EOSextra, quant1, quant2);
      *dfunofrho=compute_dSdrho_idealgas(EOSextra, quant1, quant2);
      *dfunofu=compute_dSdu_idealgas(EOSextra, quant1, quant2);
    }
  }
  else if(eomtype==EOMCOLDGRMHD || eomtype==EOMFFDE || eomtype==EOMFFDE2){
    *fun=*dfunofrho=*dfunofu=0.0;
  }

  return; // done!
}

#undef GAMMA
#undef GAMMAM1
#undef IGAMMAR
