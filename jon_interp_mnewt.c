#include "decs.h"


#define DEBUG (1)
#define MAXTRIAL 1000

int mnewt(int ntrial, int mintrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf, FTYPE tolxallowed, FTYPE tolfallowed, FTYPE tolxreport, FTYPE tolfreport, FTYPE *parms, int (*usrfun)(int n, FTYPE *, FTYPE *, FTYPE *, FTYPE **, FTYPE*)) // usrfun has extra FTYPE *norm to return
{
  int i = 0, j = 0, k = 0;
  FTYPE errx=1E30, errf=1E30, d;
  FTYPE lasterrx,lasterrf;
  FTYPE dampfactor,dampfactorchange;
  static int firstc = 1;
  static int *indx;
  static FTYPE **fjac, *fvec, *pp,*xold;
  // debug stuff
  static long count = 0;
  static long lastnstep = 0;
  static long calls = 0;
  // debug stuff
  static FTYPE *startx,*startfvec;
  int usrfunreturn;
  int lowtol[2]={0,0}; // 0=errx, 1=errf
  FTYPE X[NDIM],r,th;
#if(DEBUG==2)
  FTYPE trialvalue[MAXTRIAL][NPR];
  FTYPE trialerr[MAXTRIAL][2];
  FILE  * out;
#endif
  int truetrialnum=0;
  int donesincechange;
  static int DODAMP;
  FTYPE abstol=(NUMEPSILON*50.0); // near machine precision, typically best mnewt() can do
  FTYPE abstolhard=(NUMEPSILON); // stop trying to iterate below this
  int allownewdamp,numstabletot,countstable,numdampedtot,countdamped;
  int dampdeath;
  long int mnewtstepfail;
  int mnewtifail,mnewtjfail,mnewtpartialstepfail;
  FTYPE normf;
  
#if(DEBUG)
  calls++;
#endif

  if(tolx<abstol) tolx=abstol;
  if(tolf<abstol) tolf=abstol;


  // for debug purposes
  mnewtstepfail=11;
  mnewtifail=0;
  mnewtjfail=0;
  mnewtpartialstepfail=0;

#if(ROOTMETHOD==0)
  DODAMP=0; // ==1 is bad for some reason
#elif(ROOTMETHOD==1)
  // settings
  DODAMP=2;
#endif
  dampfactor=1.0;
  dampfactorchange=0.82;
  donesincechange=0;
  dampdeath=0;
  // for DODAMP==2
  allownewdamp=0; // start fresh
  numstabletot=8;
  numdampedtot=8;
  countstable=0;
  countdamped=0;


  if (firstc) {
    firstc = 0;

    indx = ivector(1, n);
    pp = dvector(1, n);
    xold = dvector(1, n);
    startx = dvector(1, n);
    startfvec = dvector(1, n);
    fvec = dvector(1, n);
    fjac = dmatrix(1, n, 1, n);
  }


  // save guess for debugging
#if(DEBUG)
  for (i = 1; i <= n; i++){
    startx[i]=x[i];
#if(DEBUG==2)
    trialvalue[truetrialnum][i-1]=x[i];
#endif
  }
#endif
  
  // initially
  for (i = 1; i <= n; i++) xold[i]=x[i];




  ///////////////////////
  //
  // DO MNEWT LOOP
  //

  for (k = 1; k <= ntrial; k++) {
    // count how many iterations
    truetrialnum++;
    nstroke++;

    // determine error in f, conserved quantities
    usrfunreturn=(*usrfun)(n, parms, x, fvec, fjac,&normf);
#if(DEBUG)
    for(i=1;i<=n;i++){ startfvec[i]=fvec[i]; }
#endif

    // check error of conserved quantity
    if (usrfunreturn>=1){
      // fix up the damping if we get a psychotic solution
      errf=1.E30;
      errx=1.E30;
      failed=0; // force no failure condition
    }
    else{ // estimate error normally
      errf = 0.0;
      for (i = 1; i <= n; i++) errf += fabs(fvec[i]);
      errf/=normf; // renormalize since truly U (conserved quantity) has no better significance than this error

      if (errf <= abstolhard) lowtol[1]=2;
      else if (errf <= tolf) lowtol[1]=1;
      else lowtol[1]=0;
    }

    // now fix x (primitive variables)

    if((DODAMP==1)||((DODAMP==2)&&(allownewdamp)&&(k>1)&&(lasterrf<errf))){ // only consider if damped failed when using damping odd/even trial
      // if we are here, the error actually increased
      // if increased, lower damp factor and back up
      countdamped++; // coming here counts as a general damped run
      if(countdamped>=numdampedtot){ allownewdamp=0; countdamped=0; countstable=0;}
      donesincechange=0;
      dampfactor*=dampfactorchange;

      dampdeath=0;
      // how hard do we want to make the code try for the highest precision (damping)?
      if((dampfactor<1E-5)&&(errf<tolf*100.0)){
        // then we'll assume this is good enough and any smaller dampfactor won't get us much less error
        dampdeath=1;
      }
      else if(dampfactor<1E-7){
        if(debugfail>=1){
          dualfprintf(fail_file,"mnewt: ok, we really need something better\n: dampfactor: %g errf: %g\n",dampfactor,errf); 
        }
        dampdeath=1; // let the bottom non-convergence criterea handle things
      }
      if(dampdeath==0){
        for (i = 1; i <= n; i++) x[i]=xold[i];
        k--; // want absolute number of trials to be fixed
      }
    }
    else{
      // then error is decreasing or first trial, good!  So let's continue changing x
      //if(donesincechange==5) dampfactor/=dampfactorchange;
      donesincechange++;
      if(allownewdamp==0) countstable++;
      if(allownewdamp) countdamped++; // coming here counts as a general damped run
      //else       countstable++; // then counts as stable run

      if(countstable>=numstabletot){ allownewdamp=1; countdamped=0; countstable=0;}
      if(countdamped>=numdampedtot){ allownewdamp=0; countdamped=0; countstable=0;}
      lasterrf=errf;
      // save old x to go back to it
      for (i = 1; i <= n; i++) xold[i]=x[i]; 
      
      for (i = 1; i <= n; i++) pp[i] = -fvec[i];
      ludcmp(fjac, n, indx, &d);
      lubksb(fjac, n, indx, pp);
      // DAMP (faster to damp every other one)
      // 
      if(DODAMP) for (i = 1; i <= n; i++) pp[i] = dampfactor*pp[i];

      ///////
      // evaluate x error (already properly normalized since directly what we seek)
      ///////
      errx = 0.0;
      for (i = 1; i <= n; i++) {
        errx += fabs(pp[i]);
        x[i] += pp[i];
      }

      if (errx <= abstolhard) lowtol[0]=2;
      else if (errx <= tolx) lowtol[0]=1;
      else lowtol[0]=0;

    }// else if error decreased
#if(DEBUG==2)
    // some debug stuff, done every trial type
    for (i = 1; i <= n; i++) {
      trialvalue[truetrialnum][i-1]=x[i];
    }
    if(truetrialnum>MAXTRIAL){ dualfprintf(fail_file,"oops! %d %d\n",truetrialnum,MAXTRIAL); fflush(fail_file); myexit(1);}
    trialerr[truetrialnum-1][0]=errf;
    trialerr[truetrialnum-1][1]=errx;
#endif

    if(dampdeath) break; // problem with damping, too strong, etc.

    // tolerance conditions
    // change below to && since sometimes error estimate can be near 0
    if(k>mintrial){// only quit if done more than mintrial trials
      //      if((lowtol[0]==2)&&(lowtol[1]==2)) break; // end immediately since we are unable to go further anyways below machine precision
      if((lowtol[0]>=1)&&(lowtol[1]>=1)) break; // then exactly what we wanted      
      // if 0 0 or 1 0 or 0 1, then continue to try to find better solution
    }

  }// over trials


  ///////////////////////////////
  //
  // done with MNEWT, now debug stuff comes
  //
  ///////////////////////////////

  // see what's going on
  // some counting on this run of mnewt, failed or not
#if(DEBUG)
  if (lastnstep < nstep) {
    logfprintf("#1 count/zone: %g calls: %g\n",
               ((FTYPE) count) / ((FTYPE) (N1 * N2)),
               ((FTYPE)calls) / ((FTYPE)(N1 * N2)));
    logfprintf("count: %ld zones: %d calls: %ld\n",
               count,N1 * N2,calls);
    /*
      myfprintf(stderr,"count: %ld zones: %d calls: %d\n",
      count,totalzones,calls);
    */
    myfprintf(logfull_file,"#1 count/zone: %g calls: %g\n",
              ((FTYPE) count) / ((FTYPE) (totalzones)),((FTYPE)
                                                        calls) / ((FTYPE)(totalzones)));
    myfprintf(stderr,"#1 count/zone: %g calls: %g\n",
              ((FTYPE) count) / ((FTYPE) (totalzones)),((FTYPE)
                                                        calls) / ((FTYPE)(totalzones)));
    count = (long)(k - 1);
    lastnstep = nstep;
    calls = 0;
  } else {
    count += (long)(k - 1);
  }
#endif


  // determine if can leave or not
  if((lowtol[0]>=1)&&(lowtol[1]>=1)) return(0);
  // otherwise we didn't reach tolerances we wanted note that both
  // tolerances are asked for, not just one or the other, since one
  // variable may have low error but the other high error, and that's
  // not what we want.  We want both to be low.



  ////////////////////////////////////
  //
  // rest of this is error analysis (i.e. failure, or acceptable failure)
  //
  //


  // if we got here, we never converged to desired tolerance in the specified maximum number of trials

  // only report if pseudo-bad convergence i.e. not near limits since that produces too much data
  if ((errf >= tolfreport)||(errx >= tolxreport)) {
    stderrfprintf("Couldn't find solution: %g %g\n",errf,errx);
  }
  // assume won't fail if not too bad convergence
  if ((errf <= tolfallowed)&&(errx<=tolxallowed)) {
    return (0);
  }
  else{ // then something is wrong
    return(1);
  }
}
