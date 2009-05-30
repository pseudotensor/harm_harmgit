#include "decs.h"


#define DEBUG (1)

#define NONFAILMODE (1)

int mnewt(int ntrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf)
{
  FTYPE normf;
  int i = 0, j = 0, k = 0;
  FTYPE errx, errf, d;
  FTYPE lasterrx,lasterrf;
  FTYPE dampfactor,dampfactorchange;
  static int firstc = 1;
  static int *indx;
  static FTYPE **fjac, *fvec, *pp,*xold;
  // debug stuff
  static int count = 0;
  static long lastnstep = 0;
  static int calls = 0;
  // debug stuff
  static FTYPE *startx,*startfvec;
  int usrfunreturn;
  struct of_geom geom;
  FTYPE X[NDIM],r,th;

#if(DEBUG)
  calls++;
#endif


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
  for (i = 1; i <= n; i++){
    startx[i]=x[i];
  }
  
  // initially
  for (i = 1; i <= n; i++) xold[i]=x[i];
  dampfactor=1.0;
  dampfactorchange=0.5;

  for (k = 1; k <= ntrial; k++) {
    nstroke++;
    usrfunreturn=usrfun(x, n, fvec, fjac,&normf);

    // got fvec to use:
    for(i=1;i<=n;i++){ startfvec[i]=fvec[i]; }

    /*
      if((icurr+startpos[1]==402)&&(jcurr+startpos[2]==144)){
      for(i=1;i<=n;i++) { dualfprintf(fail_file, "startfvec[%d]=%21.15g startx[%d]=%21.15g x[%d]=%21.15g\n",i,startfvec[i],i,startx[i],i,x[i]); }
      for(i=1;i<=n;i++){
      for(j=1;j<=n;j++){ dualfprintf(fail_file, "fjac[%d][%d]=%21.15g ",i,j,fjac[i][j]); }
      dualfprintf(fail_file,"\n");
      }
      }
    */
    if (usrfunreturn>=1){
            
      // fix up the damping if we get a psychotic solution
      errf=1.E30;
      errx=1.E30;
      failed=0; // force no failure condition
    }
    else{ // estimate error normally
      errf = 0.0;
      for (i = 1; i <= n; i++)
	errf += fabs(fvec[i]);
      if (errf <= tolf) {
#if(DEBUG)
	if (lastnstep < nstep) {
	  fprintf(log_file,"#1 count/zone: %g calls: %d\n",
		     (FTYPE) count / ((FTYPE) (N1 * N2)),
		     calls / (N1 * N2)); fflush(log_file);
	  mpiisum0(&count,0);
	  mpiisum0(&calls,0);
	  myfprintf(logfull_file,"#1 count/zone: %g calls: %d\n",
		     (FTYPE) count / ((FTYPE) (totalzones)),
		     calls / (totalzones));
	  myfprintf(stderr,"#1 count/zone: %g calls: %d\n",
		     (FTYPE) count / ((FTYPE) (totalzones)),
		     calls / (totalzones));
	  count = k - 1;
	  lastnstep = nstep;
	  calls = 0;
	} else {
	  count += k - 1;
	}
#endif
	// return supposed to be here, fixed in another code
      }
    }
    // if we are here, the error is either not low enough or the error actually increased
    // if increased, lower damp factor and back up
    if((k>1)&&(lasterrf<errf)){
      dampfactor*=dampfactorchange;
      // how hard do we want to make the code try for the highest precision?
      if((dampfactor<1E-5)&&(errf<1E-8)){
	// then we'll assume this is good enough and any smaller dampfactor won't get us much less error
	break;
      }
      else if(dampfactor<1E-7){
	if(debugfail>=1){
	  dualfprintf(fail_file,"mnewt: ok, we really need something better\n: dampfactor: %g errf: %g\n",dampfactor,errf);	
	}
	break; // let the bottom non-convergence criterea handle things
      }
      for (i = 1; i <= n; i++) x[i]=xold[i];
      k--; // want absolute number of trials to be fixed
    }
    else{// then error is decreasing or first trial, good!  So let's continue changing x
      lasterrf=errf;
      // save old x to go back to it
      for (i = 1; i <= n; i++) xold[i]=x[i]; 
      
      for (i = 1; i <= n; i++)
	pp[i] = -fvec[i];
      ludcmp(fjac, n, indx, &d);
      lubksb(fjac, n, indx, pp);
      // DAMP
      for (i = 1; i <= n; i++)
	pp[i] = dampfactor*pp[i];
      
      /*
	if((icurr+startpos[1]==402)&&(jcurr+startpos[2]==144)){
	for(i=1;i<=n;i++){
	dualfprintf(fail_file,"pp[%d]: %21.15g\n",i,pp[i]);
	}
	}
      */
      errx = 0.0;
      for (i = 1; i <= n; i++) {
	errx += fabs(pp[i]);
	x[i] += pp[i];
      }
      // don't care how close p is, just care how close U is.
      if (0&&(errx <= tolx)) {
#if(DEBUG)
	if (lastnstep < nstep) {
	  trifprintf("#2 count/zone: %g calls: %d\n",
		     (FTYPE) count / ((FTYPE) (N1 * N2)),
		     calls / (N1 * N2));
	  fflush(log_file);
	  count = k;
	  lastnstep = nstep;
	  calls = 0;
	} else {
	  count += k;
	}
#endif
	return (0);


      }


    }
  }
  // if we got here, we never converged to desired tolerance in the specified maximum number of trials

  // only report if pseudo-bad convergence i.e. not near limits since that produces too much data
  if ((errf >= 1E-8)||(errx >= 1E-8)) {
    if(debugfail>=2){
      dualfprintf(fail_file,"proc: %d, mnewt didn't converge: i=%d j=%d, errf: %g errx: %g dampfactor: %g\n",myid, startpos[1]+icurr,startpos[2]+jcurr,errf,errx,dampfactor);
    }
  }
  // assume won't fail if not too bad convergence if <=1E-4
  if ((errf <= 1E-4)&&(errx<=1E-4)) {
    return (0); // for now
  }
  else{ // if >1E-4, then something is wrong
    if(debugfail>=1){
      dualfprintf(fail_file, "mnewt:usrfun: (k=%d) failure\n", k);
    }
    if(debugfail>=2){
      // COMMENT
      // if failed, probably went outside allowed solution given constraints
      // placed on p.  How can this occur when U and p are related?  Is U somehow
      // more flexible?
      coord(icurr,jcurr,CENT,X);
      bl_coord(X,&r,&th);
      get_geometry(icurr,jcurr,CENT,&geom);
      dualfprintf(fail_file,"i=%d j=%d \nx1=%21.15g x2=%21.15g \nr=%21.15g th=%21.15g \ng=%21.15g\n",startpos[1]+icurr,startpos[2]+jcurr,X[1],X[2],r,th,geom.g);
      for(i=1;i<=n;i++) { dualfprintf(fail_file, "startfvec[%d]=%21.15g startx[%d]=%21.15g x[%d]=%21.15g\n",i,startfvec[i],i,startx[i],i,x[i]); }
      for(i=1;i<=n;i++){
	for(j=1;j<=n;j++){ dualfprintf(fail_file, "fjac[%d][%d]=%21.15g ",i,j,fjac[i][j]); }
	dualfprintf(fail_file,"\n");
      }
      //failed = 3;		// source of failure (nonconvergence)
    }
    FAILSTATEMENT("mnewt.c", "convergence", 1);
  }
}
