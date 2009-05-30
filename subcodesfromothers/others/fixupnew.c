#include "decs.h"



int inflow_check_4vel(int dir, double *pr, struct of_geom *ptrgeom)
{
  int ii,jj;
  int iin,iout;

  if(dir==2) return(0);

  ii=ptrgeom->i;
  jj=ptrgeom->j;

  // for dir==1
  if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ iin=-1; iout=totalsize[1]; }
  else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ iin=0; iout=totalsize[1]; }
  
  if( 
     ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(pr[dir] > 0.)) 
     ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(pr[dir] < 0.)) 
     ) {
    pr[U1]=0;
  }
  return(0);
}

int inflow_check_3vel(int dir, double *pr, struct of_geom *ptrgeom)
{
  int ii,jj;
  int iin,iout;

  if(dir==2) return(0);


  ii=ptrgeom->i;
  jj=ptrgeom->j;

  // for dir==1
  if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ iin=-1; iout=totalsize[1]; }
  else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ iin=0; iout=totalsize[1]; }
  
  if( 
     ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(pr[dir] > 0.)) 
     ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(pr[dir] < 0.)) 
     ) {
    pr[U1]=0;
  }
  return(0);
}

int inflow_check_rel4vel(int dir, double *pr, struct of_geom *ptrgeom)
{
  double ucon[NDIM] ;
  int j,k ;
  double alpha,betacon1,gamma,vsq ;
  int ii,jj;
  int iin,iout;

  if(dir==2) return(0);


  ii=ptrgeom->i;
  jj=ptrgeom->j;

  ucon_calc(pr, ptrgeom, ucon) ;

  // for dir==1
  if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ iin=-1; iout=totalsize[1]; }
  else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ iin=0; iout=totalsize[1]; }
  

  if( 
     ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(ucon[dir] > 0.)) 
     ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(ucon[dir] < 0.)) 
     ) {
    /* find gamma and remove it from primitives */
    if(gamma_calc(pr,ptrgeom,&gamma)>=1){
      dualfprintf(fail_file,"gamma calc failed: inflow_check_rel4vel\n");
      if (fail(FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    pr[U1] /= gamma ;
    pr[U2] /= gamma ;
    pr[U3] /= gamma ;
    alpha = 1./sqrt(-ptrgeom->gcon[0][0]) ;
    betacon1 = ptrgeom->gcon[0][1]*alpha*alpha ;

    /* reset radial velocity so radial 4-velocity
     * is zero */
    pr[U1] = betacon1/alpha ; // gives 3-vel contravariant

    /* now find new gamma and put it back in */
    vsq = 0. ;
    //    SLOOP vsq += ptrgeom->gcon[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
    SLOOP vsq += ptrgeom->gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
    if(fabs(vsq)<1E-13) vsq=1E-13; // machine precision thing
    if(vsq>=1.0){
      dualfprintf(fail_file,"gamma calc failed: vsq=%g\n",vsq);
      if (fail(FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }

    gamma = 1./sqrt(1. - vsq) ;
    pr[U1] *= gamma ;
    pr[U2] *= gamma ;
    pr[U3] *= gamma ;

    /* done */
    return(0);
  }
  else

    return(0);

}


/* apply floors to density, internal energy */

// currently called before bound, which assumes bound sets boundary
// values exactly as wanted without any fixing.

#define JONFIXUP 1 // 0=gammie 1=jon's


#if 0
void fixup(FTYPE (*pv)[NH + 4][NPR],SFTYPE Dt)
{
  int i, j;

  ZSLOOP(-2, NR + 1, -2, NH + 1) pfixup(pv[i][j], i, j);

}
#endif

// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
int postbc_fixup(FTYPE (*pv)[N2 + 4][NPR],SFTYPE Dt)
{
  get_bsqflags(pv);
}

#if(JONFIXUP==1)

int fixup(FTYPE (*pv)[N2 + 4][NPR],SFTYPE Dt)
{
  int i, j;
  struct of_geom geom;

  // check for bad solutions
  fixup_utoprim(pv);

  //  ZSLOOP(-2,N1+1,-2,N2+1) {
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    get_geometry(i,j,CENT,&geom);
    if(fixup1zone(pv[i][j],&geom,Dt)>=1)
      FAILSTATEMENT("fixup.c:fixup()", "fixup1zone()", 1);
  }
  return(0);
}
#else

int fixup(FTYPE (*pv)[N2 + 4][NPR],SFTYPE Dt)
{
  int i, j, k;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE uuscal, uurscal, uutscal, rhoscal, rhorscal, rhotscal;
  FTYPE ftempA,ftempB;
  struct of_geom geom ;

  // check for bad solutions
  fixup_utoprim(pv);

  //  ZSLOOP(-2,N1+1,-2,N2+1) {
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    if(1){ // choice
      rhorscal = pow(r, -1.5);
    }
    else{
      rhorscal = pow(r, -2.0);
    }
    uurscal = rhorscal / r;


    if(0){ // choice
      ftempA=(RHOMAX+RHOMIN)*0.5/RHOMIN;
      ftempB=(RHOMAX-RHOMIN)*0.5/RHOMIN;
      // choice
      if(1) rhotscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  rhotscal = ftempA+ftempB*cos(2.0*M_PI*th);
      ftempA=(UUMAX+UUMIN)*0.5/UUMIN;
      // choice (make same as above)
      ftempB=(UUMAX-UUMIN)*0.5/UUMIN;
      if(1) uutscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  uutscal = ftempA+ftempB*cos(2.0*M_PI*th);
    }
    else{
      rhotscal = 1.0;
      uutscal = 1.0;
    }


    uuscal = UUMIN*uurscal*uutscal;
    rhoscal = RHOMIN*rhorscal*rhotscal;

    /* floor on density (momentum *not* conserved) */
    if (pv[i][j][RHO] < rhoscal) {
      fladd[RHO] +=
	  dVF * gdet[i][j][CENT] * (RHOMIN * rhoscal - pv[i][j][RHO]);
      pv[i][j][RHO] = rhoscal;
    }

    /* floor on internal energy */
    if (pv[i][j][UU] < uuscal) {
      fladd[UU] +=
	  dVF * gdet[i][j][CENT] * (uuscal - pv[i][j][UU]);
      pv[i][j][UU] = uuscal;
    }

    /* limit gamma wrt normal observer */
    if(WHICHVEL==2){
      get_geometry(i,j,CENT,&geom) ;
      if(limit_gamma(GAMMAMAX,pv[i][j],&geom)>=1)
	FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
    }
    /* 
       if(bsq_calc(pv[i][j],&geom,&bsq)>=1){
       fprintf(fail_file,"1init:bsq_calc: failure\n");
       fflush(fail_file); return(1); } if(bsq >
       MAXBSQOVERRHO*pv[i][j][RHO]){
       fladd[RHO]+=dVF*gdet[i][j][CENT]*(rhoscal-pv[i][j][RHO]);
       pv[i][j][RHO] = bsq/MAXBSQOVERRHO ; } if(bsq >
       MAXBSQOVERUU*pv[i][j][UU]){
       fladd[UU]+=dVF*gdet[i][j][CENT]*(uuscal-pv[i][j][UU]);
       pv[i][j][UU] = bsq/MAXBSQOVERUU ; } */
  }

  return(0);
}


#endif



// Dt==0 is non-accounting, Dt==dt is accounting
int fixup1zone(FTYPE *pr, struct of_geom *ptrgeom, SFTYPE Dt)
{
  int k;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE scaleh[NPR],scaler[NPR],scalemin[NPR],scalefactor[NPR];
  FTYPE ftempA,ftempB;
  struct of_state q;
  FTYPE prfloor[NPR];
  FTYPE Ui[NPR],Uf[NPR];
  int checkfl[NPR];
  int failreturn;


  coord(ptrgeom->i, ptrgeom->j, ptrgeom->p, X);
  bl_coord(X, &r, &th);


  // assign general floor variables
  // whether to check floor condition
  PLOOP{ checkfl[k]=0;}
  checkfl[RHO]=1;
  checkfl[UU]=1;
  
  ////////////////////
  // scaling functions
  ////
  if(1){ // choice
    scaler[RHO] = pow(r, -1.5);
  }
  else{
    scaler[RHO] = pow(r, -2.0);
  }
  scaler[UU] = scaler[RHO] / r;
  
  
  if(0){ // choice
    ftempA=(RHOMAX+RHOMIN)*0.5/RHOMIN;
    ftempB=(RHOMAX-RHOMIN)*0.5/RHOMIN;
    // choice
    if(1) scaleh[RHO] = ftempA+ftempB*cos(2.0*M_PI*X[2]);
    else  scaleh[RHO] = ftempA+ftempB*cos(2.0*M_PI*th);
    ftempA=(UUMAX+UUMIN)*0.5/UUMIN;
    // choice (make same as above)
    ftempB=(UUMAX-UUMIN)*0.5/UUMIN;
    if(1) scaleh[UU] = ftempA+ftempB*cos(2.0*M_PI*X[2]);
    else  scaleh[UU] = ftempA+ftempB*cos(2.0*M_PI*th);
  }
  else{
    scaleh[RHO] = 1.0;
    scaleh[UU] = 1.0;
  }

  scalemin[RHO]=RHOMINLIMIT;
  scalemin[UU]=UUMINLIMIT;
  scalefactor[RHO]=RHOMIN;
  scalefactor[UU]=UUMIN;

  PLOOP prfloor[k]=scalefactor[k]*scaler[k]*scaleh[k];

  PLOOP{
    if(checkfl[k]){
      if(prfloor[k]<scalemin[k]) prfloor[k]=scalemin[k];
    }
  }
  
  if(Dt==dt){
    // before any floor changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
    failreturn=primtoU(pr,&q,ptrgeom,Ui);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
  }
  // shouldn't fail since before and after states should be ok, as
  // long as would have changed the value.  Check will occur if
  // simulation continues ok.  Could place check inside if below.
  PLOOP{
    if ( checkfl[k]&&(prfloor[k] > pr[k]) ){
      //dualfprintf(fail_file,"%d : %d %d : %d : %d : %21.15g - %21.15g\n",k,ptrgeom->i,ptrgeom->j,ptrgeom->p,checkfl[k],prfloor[k],pr[k]); 
      // only add on full step since middle step is not really updating primitive variables
      pr[k]=prfloor[k];
    }
  }
  
  if(Dt==dt){
    
    // after any floor changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(pr,&q,ptrgeom,Uf);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");

    PLOOP{
      fladd[k] += dVF * (Uf[k]-Ui[k]);
    }
  }
  /* limit gamma wrt normal observer */
  if(WHICHVEL==2)  if(limit_gamma(GAMMAMAX,pr,ptrgeom)>=1)
    FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);

  /* 
     if(bsq_calc(pr,ptrgeom,&bsq)>=1){
     fprintf(fail_file,"1init:bsq_calc: failure\n");
     fflush(fail_file); return(1); } if(bsq >
     MAXBSQOVERRHO*pr[RHO]){
     fladd[RHO]+=dVF*gdet[i][j][gpos]*(rhoscal-pr[RHO]);
     pr[RHO] = bsq/MAXBSQOVERRHO ; } if(bsq >
     MAXBSQOVERUU*pr[UU]){
     fladd[UU]+=dVF*gdet[i][j][gpos]*(uuscal-pr[UU]);
     pr[UU] = bsq/MAXBSQOVERUU ; }
  */
  return(0);
}


#define AVG4_1(pr,i,j,k) (0.25*(pr[i][j+1][k]+pr[i][j-1][k]+pr[i-1][j][k]+pr[i+1][j][k]))
#define AVG4_2(pr,i,j,k) (0.25*(pr[i+1][j+1][k]+pr[i+1][j-1][k]+pr[i-1][j+1][k]+pr[i-1][j-1][k]))

int fixup_utoprim(FTYPE (*pv)[N2 + 4][NPR])
{
  int i,j,k;
  struct of_geom geom;

  // this average only works if using 4 velocity since only then guarenteed solution is good after interpolation
  if(WHICHVEL==1) return; // just stick with static, best can do

  // first bound failure flag
  bound_pflag(pflag);


  // first check for bad solutions and try to fix based upon average surrounding solution
  if(UTOPRIMADJUST==1){
    ZSLOOP(0 , N1-1, 0, N2-1) {
      if(pflag[i][j][FLAGUTOPRIMFAIL]!=0){
	if( // but if surrounded by good values
	   (pflag[i][j+1][FLAGUTOPRIMFAIL]==0)&&
	   (pflag[i][j-1][FLAGUTOPRIMFAIL]==0)&&
	   (pflag[i+1][j][FLAGUTOPRIMFAIL]==0)&&
	   (pflag[i-1][j][FLAGUTOPRIMFAIL]==0)
	   ){
	  if(debugfail>=1) dualfprintf(fail_file,"t=%g : i=%d j=%d : utoprim corrected1\n",t,startpos[1]+i,startpos[2]+j);
	  // then average
	  // don't mess with B field, it's correct already
	  for(k=0;k<=U3;k++) pv[i][j][k]=AVG4_1(pv,i,j,k);
	}
	else if( // but if surrounded by good values
		(pflag[i+1][j+1][FLAGUTOPRIMFAIL]==0)&&
		(pflag[i+1][j-1][FLAGUTOPRIMFAIL]==0)&&
		(pflag[i-1][j+1][FLAGUTOPRIMFAIL]==0)&&
		(pflag[i-1][j-1][FLAGUTOPRIMFAIL]==0)
		){
	  if(debugfail>=1) dualfprintf(fail_file,"t=%g : i=%d j=%d : utoprim corrected2\n",t,startpos[1]+i,startpos[2]+j);
	  // then average
	  for(k=0;k<=U3;k++) pv[i][j][k]=AVG4_2(pv,i,j,k);
	}
	else{// if(0){ // don't do for now
	  if(debugfail>=1) dualfprintf(fail_file,"t=%g : i=%d j=%d : utoprim corrected3\n",t,startpos[1]+i,startpos[2]+j);
	  // otherwise some surrounding solutions are also bad, so just keep static
	  // keeping static can lead to large regions of high gamma that are static and unstable to interactions.  Thus, reset static to normal observer flow for safety.
	  // limit_gamma?  Could also reset v completely to normal observer
	  // if(limit_gamma(1.0,pv[i][j],&geom)>=1)
	  //	  FAILSTATEMENT("fixup.c:fixup_utoprim()", "limit_gamma()", 1);

	  // average out densities
	  for(k=0;k<=UU;k++) pv[i][j][k]=0.5*(AVG4_1(pv,i,j,k)+AVG4_2(pv,i,j,k));
	  // make velocity the normal observer
#if(WHICHVEL==2)
	  pv[i][j][U1] = pv[i][j][U2] = pv[i][j][U3] = 0.0;
#elif(WHICHVEL==1)
	  get_geometry(i,j,CENT,&geom);
	  pv[i][j][U1] = geom.gcon[0][1]/geom.gcon[0][0] ;
	  pv[i][j][U2] = geom.gcon[0][2]/geom.gcon[0][0] ;
	  pv[i][j][U3] = geom.gcon[0][3]/geom.gcon[0][0] ;
#elif(WHICHVEL==0)
	  get_geometry(i,j,CENT,&geom);
	  pv[i][j][U1] = -geom.gcon[0][1]/sqrt(-geom.gcon[0][0]) ;
	  pv[i][j][U2] = -geom.gcon[0][2]/sqrt(-geom.gcon[0][0]) ;
	  pv[i][j][U3] = -geom.gcon[0][3]/sqrt(-geom.gcon[0][0]) ;
#endif
	}
	pflag[i][j][FLAGUTOPRIMFAIL]=0; // reset
      }
      // otherwise good solution out of utoprim
    }
  }
  return(0);
}

int get_bsqflags(FTYPE (*pv)[N2 + 4][NPR])
{
  int i,j;
  FTYPE bsq;
  struct of_geom geom;
  struct of_state q;

  ZSLOOP(-1 , N1, -1, N2) { // over range dq is computed

    get_geometry(i, j, CENT, &geom);

    if (get_state(pv[i][j], &geom, &q) >= 1)
      FAILSTATEMENT("fixup.c:get_bsqflags()", "get_state()", 1);

    bsq = dot(q.bcon, q.bcov);
    
    if(bsq/pv[i][j][RHO]>BSQORHOLIMIT) pflag[i][j][FLAGBSQORHO]=1;
    else  pflag[i][j][FLAGBSQORHO]=0;
    
    if(bsq/pv[i][j][UU]>BSQOULIMIT) pflag[i][j][FLAGBSQOU]=1;
    else  pflag[i][j][FLAGBSQOU]=0;

    if(
       (((LIMADJUST==1)||(LIMADJUST==3))&&(pflag[i][j][FLAGBSQORHO]==1)) ||
       (((LIMADJUST==2)||(LIMADJUST==3))&&(pflag[i][j][FLAGBSQOU]==1))
       ){
      pflag[i][j][FLAGREALLIM]=MINM;
    }
    else pflag[i][j][FLAGREALLIM]=lim;
  }
  return(0);
}


int limit_gamma(FTYPE gammamax, FTYPE*pr,struct of_geom *geom)
{
  FTYPE f,gamma,pref;

  //PLOOP{ dualfprintf(fail_file,"i=%d j=%d pr[%d]=%21.15g\n",startpos[1]+i,startpos[2]+j,k,pv[i][j][k]);}
  if(gamma_calc(pr,geom,&gamma)>=1){
    dualfprintf(fail_file,"limit_gamma: gamma calc failed\n");
    dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,gamma);
    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }
    
  if(gamma > gammamax) {
    if(debugfail>=2) dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,gamma);
    /* rescale velocities to reduce gamma 
     * to gammamax */
    pref=(gammamax*gammamax - 1.)/
      (gamma*gamma - 1.);
    if(pref<0.0){
      dualfprintf(fail_file,"limit_gamma: pref calc failed pref=%g\n",pref);
      dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,gamma);
      if (fail(FAIL_UTCALC_DISCR) >= 1)
	return (1);
    }

    f = sqrt(pref);
    pr[U1] *= f ;	
    pr[U2] *= f ;	
    pr[U3] *= f ;	
  }
  return(0);
}

void fix_flux(double F1[][N2+4][NPR], double F2[][N2+4][NPR])
{
        int i,j,k ;
        double sth ;
        double X[NDIM],r,th ;

	if(mycpupos[2]==0){	 
	  for(i=-1;i<N1;i++) {
	    // emf
	    F1[i][-1][B2] = -F1[i][0][B2] ;
	    PLOOP F2[i][0][k] = 0. ;
	  }
	}
	if(mycpupos[2]==ncpux2-1){
	  for(i=-1;i<N1;i++) {
	    // emf
	    F1[i][N2][B2] = -F1[i][N2-1][B2] ;	    
	    PLOOP F2[i][N2][k] = 0. ;
	  }
	}

}



#define UTLIMIT (50.0) // limit of u^t given no other info and need to fix u
#define UTFIX (utobs) //  if failed u^t, then fix to this level to keep normal, or use local value as desired value
#define GRADIENTFACTOR (.001) // how to set?

// pr is pr to be fixed
// prmodel is pr that has a model ucon[TT], since we want to interpolate based upon ucon[TT], not just a random one when inside step_ch.c in fluxcalc()
// prmodel is just pr if in utoprim.c since we want a real change, not just an interpolation, thus will not be used as basis
int check_pr(FTYPE *pr,FTYPE *prmodel, struct of_geom *ptrgeom,int modelpos)
{
  int failedcheck;
  FTYPE ucon[NPR],uconmodel[NPR];
  FTYPE prold[NPR],probs[NPR];
  FTYPE gradient[NDIM],normsq;
  int trialcount,ntrials;
  int i,j,k;
  FTYPE lastuttdiscr,uttdiscr0,utobs;
  FTYPE dampfactor,dampfactorchange;
  FTYPE newerr,olderr;
  int method;
  struct of_geom modelgeom;
  struct of_geom *ptrmodelgeom;
  int idel,jdel;
  FTYPE realutlimit,realdiscrlimit;
  FTYPE modeldiscr;
  int utinterp;
  FTYPE toldiscr;
  int modeli,modelj;
  int bctype;


  if(modelpos==-100){ // then utoprim called us from usrfun
    modelpos=0;
    ntrials=30.0; // don't try so hard since failure is more likely, leads to static solution
  }
  else{
    ntrials=200.0; // try really hard since using observer as solution is not fun
  }
  toldiscr=1.E-2;
  dampfactorchange=0.5;
  dampfactor=1.0;
  method=1; // 0=GRADIENTFACTOR 1=damp newton's method
  utinterp=0; // whether to fix interpolation (if bad) based upon a model pr
  

  if(method==1){
    // save old pr
    PLOOP probs[k]=prold[k]=pr[k];
  }

  whocalleducon=1; // turn off failures
  if(ucon_calc(pr, ptrgeom, ucon) >= 1) ucon[TT]=1E30; // bad, so keep going
  uttdiscr0=lastuttdiscr=uttdiscr; // ucon_calc() set uttdiscr
  if(ucon[TT]<UTLIMIT) return(0); // good, so just continue calculations

  if(modelpos==-2) return(1); // don't try to correct, just fail since ucon[TT] wasn't less than the limit (otherwise wouldn't reach here).  Used to just check, no fix.
  
  if(modelpos==-3){
    modelpos=-1;
    bctype=1; // so bound_prim called us
  }
  else bctype=0; // non-bound call
  
  // if we are here, original velocities need to be corrected

  // first find normal observer to get relevant u^t
  for(i=1;i<=3;i++){
    probs[U1+i-1] = ptrgeom->gcon[TT][i]/ptrgeom->gcon[TT][TT] ;
  }
  // ucon[TT] must be good for normal observer (?) hard to solve for u^t for normal observer in general
  // change of whocalleducon forces kill level check on normal observer
  whocalleducon=0;
  if(ucon_calc(probs, ptrgeom, ucon) >= 1){
    dualfprintf(fail_file,"Thought normal observer had to have good u^t!\n");
    return(1);
  }
  else utobs=ucon[TT];
  whocalleducon=1;


  if(utinterp&&(modelpos>=0)){ // use modelpos==-1 for no use of model
    // below for no interpolation but use of model
    if(modelpos==0){
      modeli=ptrgeom->i;
      modelj=ptrgeom->j;
      ptrmodelgeom=ptrgeom;
    }
    // modelpos>=1 for interpolations in fluxcalc() (model value on center)
    else if(modelpos==1){
      if(ptrgeom->p==FACE1){ idel=1; jdel=0; }
      else if(ptrgeom->p==FACE2){ idel=0; jdel=1; }
      modeli=(ptrgeom->i) -idel;
      modelj=(ptrgeom->j) -jdel;
      // then i-1,j is center position
      get_geometry(modeli,modelj,CENT,&modelgeom);
      ptrmodelgeom=&modelgeom;
    }
    else if(modelpos==2){
      modeli=ptrgeom->i;
      modelj=ptrgeom->j;
      // then i,j is center position
      get_geometry(modeli,modelj,CENT,&modelgeom);
      ptrmodelgeom=&modelgeom;
    }
    // determine model u^t
    if(ucon_calc(prmodel, ptrmodelgeom, uconmodel) >= 1){
      // model no good
      if(bctype==0){
	dualfprintf(fail_file,"serious failure.  On-grid values and fixed bc values used have u^t imaginary: modeli: %d modelj: %d uttdiscr: %21.15g\n",startpos[1]+modeli,startpos[2]+modelj,uttdiscr);
	whocalleducon=0; // turn on failures
	if (fail(FAIL_UTCALC_DISCR) >= 1)
	  return (1);
      }
      else uconmodel[TT]=realutlimit=1E30;
      // otherwise normal to sometimes encounter failure if using model in bc (which isn't currently)
    }
    else realutlimit=uconmodel[TT]; // model ut
    modeldiscr=uttdiscr;

    // upper limit (if model has large u^t, assume ok to have one)
    if(realutlimit>UTLIMIT) realutlimit=UTLIMIT;
  }
  else realutlimit=UTFIX; // no model, just fix based upon normal observer since no idea what is "ok" to have

  realdiscrlimit=1.0/(realutlimit*realutlimit);
  newerr=olderr=fabs(lastuttdiscr-realdiscrlimit)/realdiscrlimit;


  // LOOP


  // otherwise need to fix
  failedcheck=0;

  trialcount=0;
  while( ((newerr>toldiscr)&&(method==1)) ||((ucon[TT]>realutlimit)&&(method==0)) ){
    // see if we can fix things since outside limits

    // determine gradient
    normsq=0.0;
    for(i=1;i<=3;i++){
      gradient[i]=2.0*(ptrgeom->gcov[0][i]);
      for(j=1;j<=3;j++){
	// note that ucon is the same as pr here since ucon_calc sets spatial terms to pr
	gradient[i]+=2.0*ucon[j]*ptrgeom->gcov[i][j];
      }
      normsq+=gradient[i]*gradient[i];
    }
    // normalize gradient
    for(i=1;i<=3;i++){
      gradient[i]/=sqrt(normsq);
    }
    // save old pr and change new one    
    if(method==0){
      for(i=1;i<=3;i++){
	
	pr[U1+i-1]-=gradient[i]*GRADIENTFACTOR*((FTYPE)(ntrials)+1.0)/((FTYPE)(ntrials)-(FTYPE)(trialcount)+1.0);
	//pr[U1+i-1]-=gradient[i]*GRADIENTFACTOR;
      }
    }
    else if(method==1){
      for(i=1;i<=3;i++){
	prold[U1+i-1]=pr[U1+i-1];
	if(realdiscrlimit-uttdiscr>0)	pr[U1+i-1]-=gradient[i]*dampfactor;
	else 	pr[U1+i-1]+=gradient[i]*dampfactor;
      }
    }
    // get new ucon[TT], is it ok now?
    if(ucon_calc(pr, ptrgeom, ucon) >= 1) ucon[TT]=1E30; // bad bad, keep going
    newerr=fabs(uttdiscr-realdiscrlimit)/realdiscrlimit;
    if((method==1)&&(newerr>=olderr)){
      // then went too far (if going in right direction at all)
      dampfactor*=dampfactorchange;
      if(dampfactor<1E-10){
	if((fabs(ucon[TT]-realutlimit)/realutlimit)<0.5) break; // just be happy you got out alive
	else{
	  failedcheck=1;
	  if(debugfail>=1) dualfprintf(fail_file,"dumpfactor reached min\n");
	  break;
	}
      }
      // revert to old pr and start again
      for(i=1;i<=3;i++){
	pr[U1+i-1]=prold[U1+i-1];
      }
    }
    else{
      trialcount++; // only iterate if good step
      olderr=newerr;
    }
    if(debugfail>=2) {
      if((myid==0)&&(ptrgeom->i==117)&&(ptrgeom->j==-1)){
	dualfprintf(fail_file,"trialcount=%d uttdiscr0=%21.15g uttdiscr=%21.15g newerr: %21.15g dampfactor=%21.15g\n",trialcount,uttdiscr0,uttdiscr,newerr,dampfactor);
      }
    }
    // even if not bad, could still be too large, so check
    if(trialcount==ntrials){
      if((fabs(ucon[TT]-realutlimit)/realutlimit)<0.5) break; // just be happy you got out alive
      else{
	failedcheck=1;
	if(debugfail>=1) dualfprintf(fail_file,"number of trials reached max\n");
	break;
      }
    }
  }
  whocalleducon=0; // turn back on failures
  if(debugfail>=2){
    if((myid==0)&&(icurr==117)&&(jcurr==-1)){
      dualfprintf(fail_file,"ucon[TT]=%21.15g\n",ucon[TT]);
    }
  }

  if(failedcheck){
    if(debugfail>=1) dualfprintf(fail_file,"couldn't fix ucon[TT]=%21.15g newerr=%21.15g uttdiscr=%21.15g\n",ucon[TT],newerr,uttdiscr);
    if(debugfail>=1) dualfprintf(fail_file,"check_pr failure: t=%g , couldn't fix ucon: i=%d j=%d p=%d ucon[TT]=%g realutlimit=%g uconmodel[TT]=%g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,ucon[TT],realutlimit,uconmodel[TT]);
    if(debugfail>=1) dualfprintf(fail_file,"uttdiscr0=%21.15g uttdiscr=%21.15g realdiscrlimit=%21.15g modeldiscr=%21.15g dampfactor=%21.15g\n",uttdiscr0,uttdiscr,realdiscrlimit,modeldiscr,dampfactor);
    // will still run perhaps with UT>realutlimit, could stop it, but won't for now      
    if(debugfail>=1){
      PLOOP{
	dualfprintf(fail_file,"pr[%d]=%21.15g prmodel[%d]=%21.15g\n",k,pr[k],k,prmodel[k]);
      }
      if(debugfail>=1) dualfprintf(fail_file,"need better algorithm: check_pr failure: couldn't fix ucon: i=%d j=%d p=%d ucon[TT]=%g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,ucon[TT]);
    }

    // force to normal observer solution
    PLOOP pr[k]=probs[k];
  }

  return(0);
}
