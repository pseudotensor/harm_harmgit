#include "decs.h"


#define DEBUG (2)
#define MAXTRIAL 1000

#define NONFAILMODE (1)

int mnewt(int ntrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf)
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
  struct of_geom geom;
  FTYPE X[NDIM],r,th;
#if(DEBUG==2)
  FTYPE trialvalue[MAXTRIAL][NPR-3];
  FTYPE trialerr[MAXTRIAL][2];
#endif
  int truetrialnum=0;
  int donesincechange;
  static int DODAMP;
  // 0=no
  // 1=yes globally (kinda bad)
  // 2=yes, flipping on/off damping every so number of iterations, and letting new damped value settle in for a number of iterations
  int allownewdamp,numstabletot,countstable,numdampedtot,countdamped;
  
#if(DEBUG)
  calls++;
#endif

  // settings
  DODAMP=1;
  dampfactor=1.0;
  dampfactorchange=0.5;
  donesincechange=0;
  if((myid==2)&&(icurr+startpos[1]==0)&&(jcurr+startpos[2]==230)&&(realnstep==198858)){
    
    DODAMP=2;
    allownewdamp=0; // start fresh
    numstabletot=5;
    numdampedtot=5;
    countstable=0;
    countdamped=0;
    
    // testing, works well, turn on the if just inside the ntrial loop
    //    DODAMP=0;
  }

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
#if(DEBUG==2)
    trialvalue[truetrialnum][i-1]=x[i];
#endif
  }
  
  // initially
  for (i = 1; i <= n; i++) xold[i]=x[i];

  for (k = 1; k <= ntrial; k++) {
    /*
      if(DODAMP==0){
      if((k==10)&&(myid==2)&&(icurr+startpos[1]==0)&&(jcurr+startpos[2]==230)&&(realnstep==198858)){
      DODAMP=1;
      dualfprintf(fail_file,"true=%d\n",truetrialnum); fflush(fail_file);
      }
      }
    */
    truetrialnum++;
    nstroke++;
    usrfunreturn=usrfun(x, n, fvec, fjac);

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
      for (i = 1; i <= n; i++)	errf += fabs(fvec[i]);
      if (errf <= tolf) {
#if(DEBUG)
	if (lastnstep < nstep) {
	  fprintf(log_file,"#1 count/zone: %g calls: %g\n",
		     ((FTYPE) count) / ((FTYPE) (N1 * N2)),
		     ((FTYPE)calls) / ((FTYPE)(N1 * N2))); fflush(log_file);
	  fprintf(log_file,"count: %ld zones: %d calls: %d\n",
		     count,N1 * N2,calls); fflush(log_file);
	  mpildsum0(&count,0);
	  mpildsum0(&calls,0);
	  /*
	  myfprintf(stderr,"count: %ld zones: %d calls: %d\n",
		     count,totalzones,calls); fflush(log_file);
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
	// see what's going on
        if((myid==2)&&(icurr+startpos[1]==0)&&(jcurr+startpos[2]==230)&&(realnstep==198858)){
	  dualfprintf(fail_file,"1trueerrf=%21.15g\n",errf);
          errf=1E30; // pretend
          break;
        }

	return(0);// we are good!
      }
    }
  

    if(((DODAMP==1)||(DODAMP==2)&&(allownewdamp))&&(k>1)&&(lasterrf<errf)){ // only consider if damped failed when using damping odd/even trial
      // if we are here, the error actually increased
      // if increased, lower damp factor and back up
      countdamped++; // coming here counts as a general damped run
      if(countdamped>=numdampedtot){ allownewdamp=0; countdamped=0; countstable=0;}
      donesincechange=0;
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
      for (i = 1; i <= n; i++){
	x[i]=xold[i];
#if(DEBUG==2)
	trialvalue[truetrialnum][i-1]=x[i];
#endif
      }
#if(DEBUG==2)
      if(truetrialnum>MAXTRIAL){ dualfprintf(fail_file,"oops! %d %d\n",truetrialnum,MAXTRIAL); fflush(fail_file); myexit(1);}
#endif

      k--; // want absolute number of trials to be fixed
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
      
      for (i = 1; i <= n; i++)	pp[i] = -fvec[i];
      ludcmp(fjac, n, indx, &d);
      lubksb(fjac, n, indx, pp);
      // DAMP (faster to damp every other one)
      // 
      if(DODAMP) for (i = 1; i <= n; i++)	pp[i] = dampfactor*pp[i];
      
      errx = 0.0;
      for (i = 1; i <= n; i++) {
	errx += fabs(pp[i]);
	x[i] += pp[i];
#if(DEBUG==2)
	trialvalue[truetrialnum][i-1]=x[i];
#endif
      }
#if(DEBUG==2)
      if(truetrialnum>MAXTRIAL){ dualfprintf(fail_file,"oops! %d %d\n",truetrialnum,MAXTRIAL); fflush(fail_file); myexit(1);}
#endif
      if((myid==2)&&(icurr+startpos[1]==0)&&(jcurr+startpos[2]==230)&&(realnstep==198858)){
	dualfprintf(fail_file,"#3: %d %d %d %d %d %21.15g %21.15g\n",allownewdamp,countstable,countdamped,k,truetrialnum,errf,errx);
      }
      // don't care how close p is, just care how close U is.  well,
      // if near machine precision, never going to do better, just
      // flipping primitive variables around at machine level
      // shouldn't influence conserved quantities that much
#if(DEBUG==2)
      trialerr[truetrialnum-1][0]=errf;
      trialerr[truetrialnum-1][1]=errx;
#endif
      if (0&&(errx <= tolx)) {
#if(DEBUG)
	if (lastnstep < nstep) {
	  trifprintf("#2 count/zone: %g calls: %g\n",
		     ((FTYPE) count) / ((FTYPE) (N1 * N2)),
		     ((FTYPE)calls) /((FTYPE) (N1 * N2)));
	  fflush(log_file);
	  count = (long)(k);
	  lastnstep = nstep;
	  calls = 0;
	} else {
	  count += (long)(k);
	}
#endif
	return (0);
      }// tolx
    }// else if error decreased
  }// over trials

  // see what's going on
  if((myid==2)&&(icurr+startpos[1]==0)&&(jcurr+startpos[2]==230)&&(realnstep==198858)){
    dualfprintf(fail_file,"2trueerrf=%21.15g\n",errf);
    errf=1E30; // pretend
  }


  ////////////////////////////////////
  //
  // rest of this is error analysis (i.e. failure)
  //
  //


  // if we got here, we never converged to desired tolerance in the specified maximum number of trials

  // only report if pseudo-bad convergence i.e. not near limits since that produces too much data
  if ((errf >= 1E-8)||(errx >= 1E-8)) {
    if(debugfail>=2){
      dualfprintf(fail_file,"proc: %d, t=%21.15g realnstep=%ld mnewt didn't converge (k=%d true=%d): i=%d j=%d, errf: %g errx: %g dampfactor: %g\n",myid,t,realnstep, k, truetrialnum, startpos[1]+icurr,startpos[2]+jcurr,errf,errx,dampfactor);
    }
  }
  // assume won't fail if not too bad convergence if <=1E-4
  if ((errf <= 1E-4)&&(errx<=1E-4)) {
    return (0); // for now
  }
  else{ // if >1E-4, then something is wrong
    if(debugfail>=1){
      dualfprintf(fail_file, "mnewt:usrfun: (k=%d truetrialnum=%d) failure\n", k,truetrialnum);
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
      dualfprintf(fail_file,"x->%21.15g``20, y->%21.15g``20\n",X[1],X[2]);
      for(i=1;i<=n;i++) { dualfprintf(fail_file, "startfvec[%d]=%21.15g startx[%d]=%21.15g x[%d]=%21.15g\n",i,startfvec[i],i,startx[i],i,x[i]); }
      for(i=1;i<=n;i++){
	for(j=1;j<=n;j++){
	  dualfprintf(fail_file, "fjac[%d][%d]=%21.15g ",i,j,fjac[i][j]);
	}
	dualfprintf(fail_file,"\n");
      }
      
#if(DEBUG==2)
      dualfprintf(fail_file,"mnewt={");
      for(i=0;i<n;i++){ // seperately for each primitive variable
	dualfprintf(fail_file,"{");
	for(j=0;j<truetrialnum;j++){	  
	  dualfprintf(fail_file,"%21.15g``20 ",trialvalue[j][i]);
	  if(j<truetrialnum-1) 	    dualfprintf(fail_file,",");	    
	}
	dualfprintf(fail_file,"}\n");
	if(i<n-1) 	    dualfprintf(fail_file,",");	    
      }
      dualfprintf(fail_file,"}\n");
      
      dualfprintf(fail_file,"mnewterr={");
      for(i=0;i<=1;i++){ // seperately for each primitive variable
	dualfprintf(fail_file,"{");
	for(j=0;j<truetrialnum;j++){	  
	  dualfprintf(fail_file,"%21.15g``20 ",trialerr[j][i]);
	  if(j<truetrialnum-1) 	    dualfprintf(fail_file,",");	    
	}
	dualfprintf(fail_file,"}\n");
	if(i<1) 	    dualfprintf(fail_file,",");	    
      }
      dualfprintf(fail_file,"}\n");
      
#endif
      //failed = 3;		// source of failure (nonconvergence)
    }
    FAILSTATEMENT("mnewt.c", "convergence", 1);
  }
}
