
/* 
 *
 * invert U (conserved variables) to obtain
 * p (primitive variables).  
 * NB: pr must contain initial guess
 *
 *
 * new version using Del Zanna's scheme for reducing
 * complexity of solution.
 *
 * cfg, ldz 16 jan 03
 *
 */

#include "decs.h"

FTYPE Dc, Ec, Qsq, Sc, Bsq, inormc, Tsq;
FTYPE Wc;

/* pr *MUST* contain initial guess */
int Utoprim_ldz(FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr)
{

  FTYPE alpha, beta[NDIM], Tcov[NDIM], Tcon[NDIM];
  FTYPE B[NDIM], Q[NDIM], W, vsq, ut, gamma, v[NDIM];
  int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W);
  FTYPE ucon[NDIM];
  int i = 0, j = 0, k = 0;
  char testtext[NPR][50];
  // debug stuff
  FTYPE pr0[NPR];
  struct of_state q;
  FTYPE Ustart[NPR];
  FTYPE X[NDIM],r,th;
  int gotwack;


  // assume guess has limited gamma
  alpha = 1. / sqrt(-ptrgeom->gcon[0][0]);
  beta[0] = 0.;
  SLOOPA beta[j] = ptrgeom->gcon[0][j] * alpha * alpha;

  /* undo density addition, renormalized conserved variables */
  U[UU] -= U[RHO];
  PLOOP U[k] /= ptrgeom->g;

  // go ahead and evolve the field
  pr[B1] = U[B1];
  pr[B2] = U[B2];
  pr[B3] = U[B3];


  PLOOP pr0[k]=pr[k];
  // for debug
  /*
  if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
  if (primtoU(pr0, &q, ptrgeom, Ustart) >= 1)
    FAILSTATEMENT("utoprim_ldz.c:utoprim()", "primtoU()", 1);
  

  if(t>413.84){
    if((ptrgeom->i+startpos[1]==32)&&(ptrgeom->j+startpos[2]==0)){
      PLOOP{
	dualfprintf(fail_file,"32: pr[%d]=%21.15g\n",k,pr[k]);
      }
    }
    if((ptrgeom->i+startpos[1]==33)&&(ptrgeom->j+startpos[2]==0)){
      PLOOP{
	dualfprintf(fail_file,"33: pr[%d]=%21.15g\n",k,pr[k]);
      }
    }
  }
  // end for debug
  */



  /* raise index on mhd stress tensor */
  DLOOPA Tcov[j] = U[j + UU];
  raise(Tcov, ptrgeom, Tcon);

  Dc = alpha * U[0];
  SLOOPA Q[j] = alpha * (Tcon[j] + beta[j] * Tcon[0]);
  Ec = alpha * alpha * Tcon[0];
  SLOOPA B[j] = alpha * U[B1 + j - 1];

  /* calculate Q^2 */
  Qsq = 0.;
  SLOOP Qsq += Q[j] * Q[k] * ptrgeom->gcov[j][k];

  /* calc. S */
  Sc = 0.;
  SLOOP Sc += ptrgeom->gcov[j][k] * Q[j] * B[k];

  /* calc. Bsq */
  Bsq = 0.;
  SLOOP Bsq += ptrgeom->gcov[j][k] * B[j] * B[k];

  Tsq = Bsq * Qsq - Sc * Sc;

  /* calc. guess for W */
  // get_state hides ucon which hides type of pr we got
  if (ucon_calc(pr, ptrgeom, ucon) >= 1)
    FAILSTATEMENT("utoprim_ldz.c:Utoprim_ldz()", "ucon_calc()", 1);


  vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);
  if(fabs(vsq)<1E-13) vsq=1E-13; // just a machine precision fix (say when ucon==0)
  //  if(fabs(vsq-1.0)<1E-13) vsq=1.0-1E-13; // should never correct if vsq>1.0
  if((vsq>=1.0)||(vsq<0.0)){
    dualfprintf(fail_file,"t=%g : tofail in ldz: !! i=%d j=%d initial vsq=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,vsq);
  }

  inormc = 1.0/U[RHO];

  /* now solve */
  if(wvsq_solv_ldz(&vsq, &W)>=1){
    if(debugfail>=1) dualfprintf(fail_file,"t=%g : i=%d j=%d : failed to solve in ldz\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j);
    /*
    if (fail(FAIL_LDZ) >= 1)
    return(1);
    */
    // Flag bad solution to be treated later by fixup()
    pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]= 1;

    // note, we keep other variables static as reference for fixing
    // continue to evolve B (must conserve divb==0)
    //    if(debugfail>=1) PLOOP dualfprintf(fail_file,"%d %d U[%d]=%21.15g pr[%d]=%21.15g\n",ptrgeom->i+startpos[1],ptrgeom->j+startpos[2],k,U[k],k,pr[k]);
    return(0);
  }
  if((vsq>=1.0)||(vsq<=0.0)){
    if(fabs(vsq)<1E-10) vsq=1E-10; // machine precision thing
    else{
      dualfprintf(fail_file,"t=%g : unexpected result from utoprim_ldz: vsq=%21.15g passed\n",t,vsq);
      if (fail(FAIL_LDZ) >= 1)
	return(1);
    }
  }
  gamma = 1. / sqrt(1. - vsq);
  pr[RHO] = Dc / gamma;
  pr[UU] = ((1. - vsq) * W - pr[RHO]) / gam;

  SLOOPA v[j] = (Q[j] + (Sc / W) * B[j]) / (W + Bsq);

  if(WHICHVEL==VEL4){
    pr[U1] = gamma * (v[1] - beta[1]/alpha);
    pr[U2] = gamma * (v[2] - beta[2]/alpha);
    pr[U3] = gamma * (v[3] - beta[3]/alpha);
  }
  else if(WHICHVEL==VEL3){
    pr[U1] = alpha * v[1] - beta[1];
    pr[U2] = alpha * v[2] - beta[2];
    pr[U3] = alpha * v[3] - beta[3];
  }
  else if(WHICHVEL==VELREL4){
    pr[U1] = gamma*v[1] ;
    pr[U2] = gamma*v[2] ;
    pr[U3] = gamma*v[3] ;
  }

  // completely good solution
  pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]= 0;



  // check to make sure we get real solution, not just nan
  PLOOP{
    sprintf(testtext[k],"%g",pr[k]);
    if((!strcmp("nan",testtext[k]))||(!strcmp("inf",testtext[k]))||(!strcmp("-inf",testtext[k]))){
      if(debugfail>=1){
	dualfprintf(fail_file,"(SHOULD NEVER REACH THIS ANYMORE!) nan found as solution for ldz: i=%d j=%d p=%d k=%d\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,k);
	dualfprintf(fail_file,"failed to find solution in ldz\n");
	//PLOOP 	dualfprintf(fail_file,"ldz nan: pr[%d]=%21.15g U=%21.15g\n",k,pr[k],U[k]);
	//dualfprintf(fail_file,"vsq=%21.15g gamma=%21.15g\n",vsq,gamma);
      }
      if (fail(FAIL_LDZ) >= 1)
	return(1);
    }
  }
  // check for wacked solution (i.e. crazy but not nan)
  gotwack=0;
  PLOOP{
    if(fabs(pr[k])>prMAX[k]){
      dualfprintf(fail_file,"t=%21.15g wack solution found as solution for ldz: i=%d j=%d p=%d k=%d pr=%21.15g prMAX=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,k,pr[k],prMAX[k]);
      /*
      dualfprintf(fail_file,"wacked value: pr[%d]: %21.15g > %21.15g\n",k,pr[k],prMAX[k]);
      coord(ptrgeom->i,ptrgeom->j,ptrgeom->p,X);
      bl_coord(X,&r,&th);
      dualfprintf(fail_file,"i=%d j=%d \nx1=%21.15g x2=%21.15g \nr=%21.15g th=%21.15g \ng=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,X[1],X[2],r,th,ptrgeom->g);
      dualfprintf(fail_file,"ldz wack: pr0=\n");
      PLOOP 	dualfprintf(fail_file,"%21.15g\n",pr0[k]);
      dualfprintf(fail_file,"ldz wack: U0=\n");
      PLOOP 	dualfprintf(fail_file,"%21.15g\n",Ustart[k]);
      dualfprintf(fail_file,"ldz wack: Utarget=\n");
      PLOOP 	dualfprintf(fail_file,"%21.15g\n",U[k]);
      dualfprintf(fail_file,"ldz wack: pr=\n");
      PLOOP 	dualfprintf(fail_file,"%21.15g\n",pr[k]);
      PLOOP 	dualfprintf(fail_file,"vsq=%21.15g gamma=%21.15g\n",vsq,gamma);
      
      if (fail(FAIL_LDZ) >= 1)
	return(1);
      */
      gotwack=1;
    }
  }
  if(gotwack){
    // only got here because ldz algorithm believed it converged but it didn't really
    pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]= 1;

    // note, we keep other variables static as reference for fixing
    // continue to evolve B (must conserve divb==0)
    PLOOP pr[k]=pr0[k];// revert
    // but update field as necessary
    return(0);
  }
  return (0);

}

int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W)
{
  FTYPE tol, x1, x2;
  FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2, FTYPE xacc);
  extern void func(FTYPE *x, FTYPE *f, FTYPE *df);
  int nrunsafe(void (*funcd) (FTYPE*,FTYPE*,FTYPE*), FTYPE *guess);
  int nrunsafeorig(void (*funcd) (FTYPE*,FTYPE*,FTYPE*), FTYPE *guess);

  x1 = 0. - SMALL;
  x2 = 1. - 1.e-6;

  // x1 = 0.5*(*vsq) ;
  // x2 = 0.5 + 0.5*(*vsq) ;

  //  if(nrunsafe(func, vsq)>=1){
    if(nrunsafeorig(func, vsq)>=1){
    if(debugfail>=1) dualfprintf(fail_file,"t=%g : nrsafe failed\n",t);
    return(1);
  }
  // sometimes, however unlikely, it can (and has) happen that vsq is identically 1.0 to machine precision.  This result can, sometimes, have good enough error to pass nrunsafe.  We however cannot allow this, so fail on this case.
  //  if(*vsq>=1.0){
  // don't fail this since not so bad as a real failure unless actually >=1.0, we only want to rescale velocities, not screw up solution
  //  if(*vsq>=1.0-1.0/(GAMMAMAX*GAMMAMAX)){
  if( (*vsq>=1.0-1.0/(GAMMAFAIL*GAMMAFAIL))||(*vsq<0.0)){
    if(fabs(*vsq)<1E-10) *vsq=1E-10; // machine precision thing
    else{    
      if(debugfail>=1) dualfprintf(fail_file,"t=%g : nrunsafe passed vsq>~1.0  : vsq=%21.15g passed\n",t,*vsq);
      return(1); // i.e. failure, even if passed tolerance requirements
    }
  }


  // tol = 1.e-4 ;
  // *vsq = rtsafe(func,x1,x2,tol) ;

  *W = Wc;

  return(0);

}

/*
// low error mode
// choice
#define TOL 		1.e-15
#define MINTOL          1.e-5
#define NITERGOODMAX    50
#define NITERMAX  	50
#define NITERMIN  	2
*/

// middle error mode
// choice
#define TOL             1.e-10
#define MINTOL          1.e-5
#define NITERGOODMAX    20
#define NITERMAX        50 // let the solution really try
#define NITERMIN        3

/*
// high error mode
#define TOL 		1.e-5
#define MINTOL          1.e-3
#define NITERGOODMAX    5
#define NITERMAX  	10
#define NITERMIN  	0
*/

int nrunsafeorig(void (*funcd)(FTYPE*,FTYPE*,FTYPE*), FTYPE* guess)
{
  double f,df ;
  int n_iter ;

  (*funcd)(guess,&f,&df);

  n_iter = 0 ;
  while((fabs(f) > TOL && n_iter < NITERMAX)||(n_iter<NITERMIN)) {
    *guess -= f/df ;
    (*funcd)(guess,&f,&df);
    n_iter++ ;
    nstroke++;
  }

  if(n_iter == NITERMAX) {
    fprintf(stderr,"max iterations exceeded\n") ;
    return(1);
  }
  return(0);
}


int nrunsafe(void (*funcd) (FTYPE*, FTYPE*,FTYPE*), FTYPE *guess)
{
  FTYPE f, df;
  int n_iter,n_itergood;
  FTYPE lastf,lastdf;
  FTYPE dampfactor,dampfactorchange;
  FTYPE oldguess;
  int dodamp;
  FTYPE alpha;


  dodamp=1; // whether to use damped method
  dampfactor=1.0;
  dampfactorchange=0.8;
  alpha=1.0E-4; // used to base damp factor use (9.7 numerical recipies)
  lastf=f=1E30; // starting error
  oldguess=*guess;

  //  (*funcd) (guess, &f, &df);


  // f here is the error in the function being 0.
  // df here is the gradient of f w.r.t. x, our guess
  // guess is x, which df=df/dx.

  n_iter = n_itergood=0;
  while ((fabs(f) > TOL && n_iter < NITERMAX && n_itergood < NITERGOODMAX) || (n_iter < NITERMIN)) {
    nstroke++;
    (*funcd) (guess, &f, &df);
    if(n_iter==0){
      if(fabs(f)==1E30){
	// then absolute failure since initial guess should always be good!
	dualfprintf(fail_file,"initial guess failure: n_iter=%d guess=%g lastf=%g f=%g lastdf=%g df=%g Wc=%g\n",n_iter,oldguess,*guess,lastf,f,lastdf,df,Wc);
	return(1);
      }
    }
    if(dodamp){
      if(fabs(f)>=fabs(lastf)){
	//if(fabs(f)>=fabs(lastf)+alpha*df*(*guess-oldguess)){
	if(debugfail>=2){
	  if(fabs(f)==1E30){// then physics failure, output data
	    dualfprintf(fail_file,"reduced dampfactor: n_iter=%d oldguess=%g guess=%g lastf=%g f=%g lastdf=%g df=%g Wc=%g\n",n_iter,oldguess,*guess,lastf,f,lastdf,df,Wc);
	  }
	}
        // then we need to damp
	dampfactor*=dampfactorchange;
	*guess=oldguess-lastf/lastdf*dampfactor; // revert to previous solution with 1 damped
      }
      else{
	// otherwise good, just continue
	oldguess=*guess;
	*guess -= f / df * dampfactor;
	lastf=f;
	lastdf=df;
	n_itergood++;
      }
    }
    else{
      *guess -= f / df;
    }
    n_iter++;
  }

  if (n_itergood == NITERMAX) {
    if(fabs(f)>MINTOL){
      if(debugfail>=1){
	dualfprintf(fail_file, "t=%g : max iterations exceeded\n",t);
	dualfprintf(fail_file,"t=%g : ldz method unable to converge: guess=%21.15g f=%21.15g dampfactor=%g\n",t,*guess,f,dampfactor);
      }
      return(1);
    }
    else return(0);
  } else
    return(0);

}

void func(FTYPE *x, FTYPE *f, FTYPE *df)
{
  FTYPE vsq, glf1, cc, igam1, ee, wsq, dw;
  FTYPE a2, a1, a0, q, r, th, w1, vp2, cosa;

  vsq = *x;
  if(vsq>=1.0){
    if(vsq<1.0+1E-10){
      if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: vsqldz=%21.15g set to 1.0-1E-10",vsq);
      vsq=1.0-1E-10;
    }
    else{
      if(debugfail>=2) dualfprintf(fail_file,"tofail in ldz: vsq=%21.15g\n",vsq);
      *f = 1E30;
      *df = 1E30;
      Wc=1E30;
      return;
    }
  }
  if(vsq<0.0){
    if(vsq>0.0-1E-10){
      if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: vsqldz=%21.15g set to 0.0+1E-10",vsq);
      vsq=1.0+1E-10;
    }
    else{
      if(debugfail>=2) dualfprintf(fail_file,"tofail in ldz: vsq=%21.15g\n",vsq);
      *f = 1E30;
      *df = 1E30;
      Wc=1E30;
      return;
    }
  }
  glf1 = sqrt(1. - vsq);
  /*
  if(glf1==0.0){
    // should never reach this if above is caught
    if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: glf1=%21.15g\n",glf1);
    *f = 1E30;
    *df = 1E30;
    Wc=1E30;
    return;
  }
  */
  igam1 = (gam - 1.) / gam;
  cc = 1. / (1. - igam1 * (1. - vsq));
  ee = (Ec - igam1 * glf1 * Dc - 0.5 * Bsq) * cc;

  /* 
     if(ee < 0.) { fprintf(stderr,"ee < 0 in func\n") ; myexit(40) ; } */

  if (Tsq > 0.) {
    a2 = (2. * Bsq - ee) / 3.;
    a1 = (Bsq - 2. * ee) * Bsq;
    a0 = 0.5 * Tsq * cc - ee * Bsq * Bsq;
    q = a1 / 3. - a2 * a2;

    if(q>=0.0){
      if(q<1E-10){
	if(debugfail>=1) dualfprintf(fail_file,"wastofail in ldz: qldz=%21.15g set to -1E-10n",q);
	q=-1E-10; // minimum, failure usually at 1E-21 level or 0
      }
      else{
	if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: qldz=%21.15g\n",q);
	*f = 1E30;
	*df = 1E30;
	Wc=1E30;
	return;
      }
    }
    
    r = 0.5 * (a1 * a2 - a0) - a2 * a2 * a2;
    cosa = r / sqrt(-q * q * q);

    if(cosa>1.0){
      if(cosa < 1.0+1E-10){
	// occurs exceedingly alot (some machine precision issue)
	if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: cosaldz=%21.15g set to 1\n",cosa);
	cosa = 1.0;
      }
      else{
	if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: cosaldz=%21.15g\n",cosa);
	*f = 1E30;
	*df = 1E30;
	Wc=1E30;
	return;
      }
    }
    else if(cosa<-1.){
      if(cosa > -1.0-1E-10){
	if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: cosa=%21.15g set to -1\n",cosa);
	cosa = -1.0;
      }
      else{
	if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: cosa=%21.15g\n",cosa);
	*f = 1E30;
	*df = 1E30;
	Wc=1E30;
	return;
      }
    }
    th = acos(cosa);

    // fprintf(stderr,">>> %g %g %g\n",q,r/sqrt(-q*q*q),th) ;

    Wc = 2. * sqrt(-q) * cos(th / 3.) - a2;
    wsq = Wc * Wc;
    
    if(Wc+Bsq==0.0){
      // no good way to correct
      if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: Wc+Bsq=%21.15g Wc=%21.15g Bsq=%21.15g\n",Wc+Bsq,Wc,Bsq);
      *f = 1E30;
      *df = 1E30;
      Wc=1E30;
      return;
    }
    if(inormc>=1E30){
      // doesn't seem to ever happen
      if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: inormc=%21.15g\n",inormc);
      *f = 1E30;
      *df = 1E30;
      Wc=1E30;
      return;
    }
    
    w1 = 1. / (Wc + Bsq);
    
    if(1.+ 2. * w1 * (Wc - ee)==0.0){
      // doesn't seem to ever happen
      if(debugfail>=1) dualfprintf(fail_file,"tofail in ldz: 1. + 2. * w1 * (Wc - ee)=%21.15g\n",1. + 2. * w1 * (Wc - ee));
      *f = 1E30;
      *df = 1E30;
      Wc=1E30;
      return;
    }
    

    vp2 = Tsq * w1 * w1;
    dw = -igam1 * cc * (Wc - 0.5 * Dc / glf1) / (1. +
						 2. * w1 * (Wc - ee));
    *f = (wsq * vsq + vp2 * (2. * Wc + Bsq) - Qsq) * (inormc * inormc);
    *df = (wsq + 2. * (Wc * dw * (vsq - vp2 * w1))) * (inormc * inormc);
  } else {
    Wc = ee;
    wsq = Wc * Wc;
    dw = -igam1 * cc * (Wc - 0.5 * Dc / glf1);
    *f = wsq * vsq - Qsq;
    *df = wsq + 2. * Wc * dw * vsq;
  }
}

#define MAXIT 100

FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2, FTYPE xacc)
{
  void nrerror();
  int j;
  FTYPE df, dx, dxold, f, fh, fl;
  FTYPE temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    for (j = 0; j < 100; j++) {
      (*funcd) (0.01 * j, &fl, &df);
      if(debugfail>=1) dualfprintf(fail_file, "%d %g %g\n", j, fl, df);
    }
    /* 
     */
    nrerror("Root must be bracketed in rtsafe");
  }
  if (fabs(fl) < SMALL)
    return x1;
  if (fabs(fh) < SMALL)
    return x2;

  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    xh = x1;
    xl = x2;
  }
  rts = 0.5 * (x1 + x2);
  dxold = fabs(x2 - x1);
  dx = dxold;
  (*funcd) (rts, &f, &df);
  for (j = 1; j <= MAXIT; j++) {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0)
	|| (fabs(2.0 * f) > fabs(dxold * df))) {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
	return rts;
    } else {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts)
	return rts;
    }
    if (fabs(dx) < xacc)
      return rts;
    (*funcd) (rts, &f, &df);
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  nrerror("Maximum number of iterations exceeded in rtsafe");
  return 0.0;
}

#undef MAXIT
