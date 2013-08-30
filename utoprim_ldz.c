
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

#define JONCORRECTIONS 0
// 0: don't use, original
// 1: use jon corrections and checks

#include "decs.h"


static int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W);
static int nrunsafe(void (*funcd) (FTYPE*, FTYPE*,FTYPE*), FTYPE *guess);
static int nrunsafeorig(void (*funcd) (FTYPE*, FTYPE*,FTYPE*), FTYPE *guess);
static void func(FTYPE *x, FTYPE *f, FTYPE *df);
static FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2,
                    FTYPE xacc);
static void funcorig(FTYPE *x, FTYPE *f, FTYPE *df);


static FTYPE Dc, Ec, Qsq, Sc, Bsq, inormc, Tsq;
static FTYPE Wc;

static PFTYPE *glpflag; // global pflag for local file



/* pr *MUST* contain initial guess */
int Utoprim_ldz(FTYPE *U, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats)
{

  FTYPE alpha, beta[NDIM], Tcov[NDIM], Tcon[NDIM];
  FTYPE B[NDIM], Q[NDIM], W, vsq, ut, gamma, v[NDIM];
  int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W);
  FTYPE ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  int i = 0, j = 0, k = 0;
  int pl;
  char testtext[NPR][50];
  // debug stuff
  FTYPE pr0[NPR];
  struct of_state q;
  FTYPE Ustart[NPR];
  FTYPE X[NDIM],r,th;
  int gotwack;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"Utoprim_ldz() called in parallel region\n");
    myexit(9366983);
  }
#endif


  // assign global int pointer to lpflag pointer
  glpflag=lpflag;

#if( WHICHVEL != VELREL4 )
  stderrfprintf("Utoprim_ldz() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1);
#endif

  // assume guess has limited gamma
  alpha = 1. / sqrt(-ptrgeom->gcon[GIND(0,0)]);
  beta[0] = 0.;
  SLOOPA(j) beta[j] = ptrgeom->gcon[GIND(0,j)] * alpha * alpha;



  // go ahead and evolve the field
  pr[B1] = U[B1];
  pr[B2] = U[B2];
  pr[B3] = U[B3];

#if(JONCORRECTIONS)
  PALLLOOP(pl) pr0[pl]=pr[pl];
  // for debug
  /*
    if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    if (primtoU(NOTHING,pr0, &q, ptrgeom, Ustart, NULL) >= 1)
    FAILSTATEMENT("utoprim_ldz.c:utoprim()", "primtoU()", 1);
  

    if(t>413.84){
    if((ptrgeom->i+startpos[1]==32)&&(ptrgeom->j+startpos[2]==0)){
    PALLLOOP(pl){
    dualfprintf(fail_file,"32: pr[%d]=%21.15g\n",pl,pr[pl]);
    }
    }
    if((ptrgeom->i+startpos[1]==33)&&(ptrgeom->j+startpos[2]==0)){
    PALLLOOP(pl){
    dualfprintf(fail_file,"33: pr[%d]=%21.15g\n",pl,pr[pl]);
    }
    }
    }
    // end for debug
    */
#endif


  /* raise index on mhd stress tensor */
  DLOOPA(j) Tcov[j] = U[j + UU];
  raise_vec(Tcov, ptrgeom, Tcon);

  Dc = alpha * U[0];
  SLOOPA(j) Q[j] = alpha * (Tcon[j] + beta[j] * Tcon[0]);
  Ec = alpha * alpha * Tcon[0];
  SLOOPA(j) B[j] = alpha * U[B1 + j - 1];

  /* calculate Q^2 */
  Qsq = 0.;
  SLOOP(j,k) Qsq += Q[j] * Q[k] * ptrgeom->gcov[GIND(j,k)];

  /* calc. S */
  Sc = 0.;
  SLOOP(j,k) Sc += ptrgeom->gcov[GIND(j,k)] * Q[j] * B[k];

  /* calc. Bsq */
  Bsq = 0.;
  SLOOP(j,k) Bsq += ptrgeom->gcov[GIND(j,k)] * B[j] * B[k];

  Tsq = Bsq * Qsq - Sc * Sc;

  /* calc. guess for W */
  // get_state hides ucon which hides type of pr we got

#if(JONCORRECTIONS)
  if (ucon_calc(pr, ptrgeom, ucon, others) >= 1)
    FAILSTATEMENT("utoprim_ldz.c:Utoprim_ldz()", "ucon_calc()", 1);
#else
  ucon_calc(pr, ptrgeom, ucon, others);
#endif

  vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);
  inormc = 1.0/U[RHO];

#if(JONCORRECTIONS)
  if(fabs(vsq)<1E-13) vsq=1E-13; // just a machine precision fix (say when ucon==0)
  //  if(fabs(vsq-1.0)<1E-13) vsq=1.0-1E-13; // should never correct if vsq>1.0
  if((vsq>=1.0)||(vsq<0.0)){
    dualfprintf(fail_file,"t=%21.15g : tofail in ldz: !! i=%d j=%d initial vsq=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,vsq);
  }
  /* now solve */
  if(wvsq_solv_ldz(&vsq, &W)>=1){
    if(debugfail>=1) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : failed to solve in ldz\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j);
    /*
      if (fail(i,j,k,FAIL_LDZ) >= 1)
      return(1);
    */
    // Flag bad solution to be treated later by fixup()
    *glpflag = UTOPRIMFAILCONVUTSQ;

    // note, we keep other variables static as reference for fixing
    // continue to evolve B (must conserve divb==0)
    //    if(debugfail>=1) PALLLOOP(pl) dualfprintf(fail_file,"%d %d U[%d]=%21.15g pr[%d]=%21.15g\n",ptrgeom->i+startpos[1],ptrgeom->j+startpos[2],pl,U[pl],pl,pr[pl]);
    return(0);
  }
  if((vsq>=1.0)||(vsq<=0.0)){
    if(fabs(vsq)<1E-10) vsq=1E-10; // machine precision thing
    else{
      dualfprintf(fail_file,"t=%21.15g : unexpected result from utoprim_ldz: vsq=%21.15g passed\n",t,vsq);
      if (fail(i,j,k,FAIL_LDZ) >= 1)
        return(1);
    }
  }

#else
  wvsq_solv_ldz(&vsq, &W);
#endif

  gamma = 1. / sqrt(1. - vsq);
  pr[RHO] = Dc / gamma;
  pr[UU] = ((1. - vsq) * W - pr[RHO]) / gam;


  *pressure = (gam-1.0)*pr[UU]; // ideal gas only allowed

  // check densities for positivity
  if((pr[RHO]<0.0)||(pr[UU]<0.0)){
    //*glpflag = 1;
    if(debugfail>=1) dualfprintf(fail_file,"utoprim_ldz found negative density: %21.15g %21.15g\n",pr[RHO],pr[UU]);
    //    return(0);
  }


  SLOOPA(j) v[j] = (Q[j] + (Sc / W) * B[j]) / (W + Bsq);

#if(WHICHVEL==VEL4)
  pr[U1] = gamma * (v[1] - beta[1]/alpha);
  pr[U2] = gamma * (v[2] - beta[2]/alpha);
  pr[U3] = gamma * (v[3] - beta[3]/alpha);
#elif(WHICHVEL==VEL3)
  pr[U1] = alpha * v[1] - beta[1];
  pr[U2] = alpha * v[2] - beta[2];
  pr[U3] = alpha * v[3] - beta[3];
#elif(WHICHVEL==VELREL4)
  pr[U1] = gamma*v[1] ;
  pr[U2] = gamma*v[2] ;
  pr[U3] = gamma*v[3] ;
#endif

  // completely good solution
  *glpflag= UTOPRIMNOFAIL;


#if(JONCORRECTIONS)
  // check to make sure we get real solution, not just nan
  PALLLOOP(pl){
    sprintf(testtext[pl],"%21.15g",pr[pl]);
    if((!strcmp("nan",testtext[pl]))||(!strcmp("inf",testtext[pl]))||(!strcmp("-inf",testtext[pl]))){
      if(debugfail>=1){
        dualfprintf(fail_file,"(SHOULD NEVER REACH THIS ANYMORE!) nan found as solution for ldz: i=%d j=%d p=%d pl=%d\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,pl);
        dualfprintf(fail_file,"failed to find solution in ldz\n");
        //PALLLOOP(pl)  dualfprintf(fail_file,"ldz nan: pr[%d]=%21.15g U=%21.15g\n",pl,pr[pl],U[pl]);
        //dualfprintf(fail_file,"vsq=%21.15g gamma=%21.15g\n",vsq,gamma);
      }
      if (fail(i,j,k,FAIL_LDZ) >= 1)
        return(1);
    }
  }
#else
  // check to make sure we get real solution, not just nan
  PALLLOOP(pl){
    sprintf(testtext[pl],"%21.15g",pr[pl]);
    if((!strcmp("nan",testtext[pl]))||(!strcmp("inf",testtext[pl]))||(!strcmp("-inf",testtext[pl]))){
      dualfprintf(fail_file,"ldz nan/inf @ i=%d j=%d, pl=%d\n",ptrgeom->i,ptrgeom->j,pl); fflush(fail_file);
      *glpflag = UTOPRIMFAILCONVW;
      PALLLOOP(pl) pr[pl]=pr0[pl];// revert
      return(0);
      // fail, but correct
    }
  }
#endif

#if(1||JONCORRECTIONS)
  // check for wacked solution (i.e. crazy but not nan)
  gotwack=0;
  PALLLOOP(pl){
    if(fabs(pr[pl])>prMAX[pl]){
      dualfprintf(fail_file,"t=%21.15g wack solution found as solution for ldz: i=%d j=%d p=%d pl=%d pr=%21.15g prMAX=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,pl,pr[pl],prMAX[pl]);
      /*
        dualfprintf(fail_file,"wacked value: pr[%d]: %21.15g > %21.15g\n",pl,pr[pl],prMAX[pl]);
        bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->p,X, V);
        r=V[1]; th=V[2];
        dualfprintf(fail_file,"i=%d j=%d \nx1=%21.15g x2=%21.15g \nr=%21.15g th=%21.15g \ng=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,X[1],X[2],r,th,ptrgeom->g);
        dualfprintf(fail_file,"ldz wack: pr0=\n");
        PALLLOOP(pl)  dualfprintf(fail_file,"%21.15g\n",pr0[pl]);
        dualfprintf(fail_file,"ldz wack: U0=\n");
        PALLLOOP(pl)  dualfprintf(fail_file,"%21.15g\n",Ustart[pl]);
        dualfprintf(fail_file,"ldz wack: Utarget=\n");
        PALLLOOP(pl)  dualfprintf(fail_file,"%21.15g\n",U[pl]);
        dualfprintf(fail_file,"ldz wack: pr=\n");
        PALLLOOP(pl)  dualfprintf(fail_file,"%21.15g\n",pr[pl]);
        PALLLOOP(pl)  dualfprintf(fail_file,"vsq=%21.15g gamma=%21.15g\n",vsq,gamma);
      
        if (fail(i,j,k,FAIL_LDZ) >= 1)
        return(1);
      */
      gotwack=1;
    }
  }
  if(gotwack){
    // only got here because ldz algorithm believed it converged but it didn't really
    *glpflag= UTOPRIMFAILCONVW;

    // note, we keep other variables static as reference for fixing
    // continue to evolve B (must conserve divb==0)
    PALLLOOP(pl) pr[pl]=pr0[pl];// revert
    // but update field as necessary
    return(0);
  }
#endif
  return (0);

}

int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W)
{
  FTYPE tol, x1, x2;
  FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2, FTYPE xacc);
  int nrunsafe(void (*funcd) (FTYPE*,FTYPE*,FTYPE*), FTYPE *guess);
  int nrunsafeorig(void (*funcd) (FTYPE*,FTYPE*,FTYPE*), FTYPE *guess);

  x1 = 0. - SMALL;
  x2 = 1. - 1.e-6;

  // x1 = 0.5*(*vsq) ;
  // x2 = 0.5 + 0.5*(*vsq) ;

#if(JONCORRECTIONS)
  if(nrunsafe(func, vsq)>=1){
    if(debugfail>=1) dualfprintf(fail_file,"t=%21.15g : nrsafe failed\n",t);
    return(1);
  }
  // sometimes, however unlikely, it can (and has) happen that vsq is identically 1.0 to machine precision.  This result can, sometimes, have good enough error to pass nrunsafe.  We however cannot allow this, so fail on this case.
  //  if(*vsq>=1.0){
  // don't fail this since not so bad as a real failure unless actually >=1.0, we only want to rescale velocities, not screw up solution
  //  if(*vsq>=1.0-1.0/(GAMMAMAX*GAMMAMAX)){
  if( (*vsq>=1.0-1.0/(GAMMAFAIL*GAMMAFAIL))||(*vsq<0.0)){
    if(fabs(*vsq)<1E-10) *vsq=1E-10; // machine precision thing
    else{    
      if(debugfail>=1) dualfprintf(fail_file,"t=%21.15g : nrunsafe passed vsq>~1.0  : vsq=%21.15g passed\n",t,*vsq);
      return(1); // i.e. failure, even if passed tolerance requirements
    }
  }


  // tol = 1.e-4 ;
  // *vsq = rtsafe(func,x1,x2,tol) ;
#else
  nrunsafeorig(funcorig, vsq);
#endif

  *W = Wc;

  return(0);

}


/*
// low error mode
// choice
#define TOL   1.e-15
#define MINTOL          1.e-5
#define NITERGOODMAX    50
#define NITERMAX   50
#define NITERMIN   4
*/

// middle error mode
// choice
#define TOL   1.e-10
#define MINTOL          1.e-5
#define NITERGOODMAX    20
#define NITERMAX   20
#define NITERMIN   3


// seems to fail actually, (causes evolution problems and sometimes failure)
/*
// high error mode
#define TOL   1.e-5
#define MINTOL          1.e-3
#define NITERGOODMAX    5
#define NITERMAX   10
#define NITERMIN   0
*/

int nrunsafeorig(void (*funcd)(FTYPE*,FTYPE*,FTYPE*), FTYPE* guess)
{
  FTYPE f,df ;
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
    stderrfprintf("max iterations exceeded\n") ;
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
  alpha=1.0E-4; // used to base damp factor use (9.7 numerical recipes)
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
        dualfprintf(fail_file,"initial guess failure: n_iter=%d oldguess=%21.15g guess=%21.15g lastf=%21.15g f=%21.15g lastdf=%21.15g df=%21.15g Wc=%21.15g\n",n_iter,oldguess,*guess,lastf,f,lastdf,df,Wc);
        return(1);
      }
    }
    if(dodamp){
      if(fabs(f)>=fabs(lastf)){
        //if(fabs(f)>=fabs(lastf)+alpha*df*(*guess-oldguess)){
        if(debugfail>=2){
          if(fabs(f)==1E30){// then physics failure, output data
            dualfprintf(fail_file,"reduced dampfactor: n_iter=%d oldguess=%21.15g guess=%21.15g lastf=%21.15g f=%21.15g lastdf=%21.15g df=%21.15g Wc=%21.15g\n",n_iter,oldguess,*guess,lastf,f,lastdf,df,Wc);
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
        dualfprintf(fail_file, "t=%21.15g : max iterations exceeded\n",t);
        dualfprintf(fail_file,"t=%21.15g : ldz method unable to converge: guess=%21.15g f=%21.15g dampfactor=%21.15g\n",t,*guess,f,dampfactor);
      }
      return(1);
    }
    else return(0);
  } else
    return(0);

}

static void func(FTYPE *x, FTYPE *f, FTYPE *df)
{
  FTYPE vsq, glf1, cc, igam1, ee, wsq, dw;
  FTYPE a2, a1, a0, q, r, th, w1, vp2, cosa;

  vsq = *x;

  if(vsq>=1.0){
    if(vsq<1.0+NUMEPSILON*100.0){
      if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: vsqldz=%21.15g set to 1.0-1E-10",vsq);
      vsq=1.0-NUMEPSILON*100.0;
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
    if(vsq>0.0-NUMEPSILON*100.0){
      if(debugfail>=2) dualfprintf(fail_file,"wastofail in ldz: vsqldz=%21.15g set to 0.0+1E-10",vsq);
      vsq=0.0;
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
     if(ee < 0.) { stderrfprintf("ee < 0 in func\n") ; myexit(40) ; } */

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

    // stderrfprintf(">>> %21.15g %21.15g %21.15g\n",q,r/sqrt(-q*q*q),th) ;

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


static void funcorig(FTYPE *x,FTYPE *f,FTYPE *df)
{
  FTYPE vsq,glf1,cc,igam1,ee,wsq,dw ;
  FTYPE a2,a1,a0,q,r,th,w1,vp2,cosa ;

  vsq = *x ;
  glf1 = sqrt(1. - vsq) ;
  igam1 = (gam - 1.)/gam ;
  cc = 1./(1. - igam1*(1. - vsq)) ;
  ee = (Ec - igam1*glf1*Dc - 0.5*Bsq)*cc ;

  if(Tsq > 0.) {
    a2 = (2.*Bsq - ee)/3. ;
    a1 = (Bsq - 2.*ee)*Bsq ;
    a0 = 0.5*Tsq*cc - ee*Bsq*Bsq ;
    q = a1/3. - a2*a2 ;
    r = 0.5*(a1*a2 - a0) - a2*a2*a2 ;
    cosa = r/sqrt(-q*q*q) ;
    if(cosa > 1.) cosa = 1. ;
    th = acos(cosa) ;

    Wc = 2.*sqrt(-q)*cos(th/3.) - a2 ;
    wsq = Wc*Wc ;
    w1 = 1./(Wc + Bsq) ;
    vp2 = Tsq*w1*w1 ;
    dw = -igam1*cc*(Wc - 0.5*Dc/glf1)/(1. + 2.*w1*(Wc - ee)) ;
    *f =  ( wsq*vsq+vp2*(2.*Wc+Bsq)-Qsq    )*(inormc*inormc) ;
    *df = ( wsq + 2.*(Wc*dw*(vsq - vp2*w1)))*(inormc*inormc) ;
  }
  else {
    Wc = ee ;
    wsq = Wc*Wc ;
    dw = -igam1*cc*(Wc - 0.5*Dc/glf1) ;
    *f = wsq*vsq - Qsq ;
    *df = wsq + 2.*Wc*dw*vsq ;
  }
}




#define MAXIT 100

FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2, FTYPE xacc)
{
  void nrerror(char error_text[]);
  int j;
  FTYPE df, dx, dxold, f, fh, fl;
  FTYPE temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    for (j = 0; j < 100; j++) {
      (*funcd) (0.01 * j, &fl, &df);
      if(debugfail>=1) dualfprintf(fail_file, "%d %21.15g %21.15g\n", j, fl, df);
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
