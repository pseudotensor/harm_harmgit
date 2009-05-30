#include "decs.h"


#define DEBUG (1)
#define DEBUGPOINT (0)
#define MAXTRIAL 1000

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
  FTYPE abstol=(NUMEPSILON*50.0); // near machine precision
  int allownewdamp,numstabletot,countstable,numdampedtot,countdamped;
  int dampdeath;
  long int mnewtstepfail;
  int mnewtifail,mnewtjfail,mnewtpartialstepfail;
  FTYPE tolfallowed,tolxallowed,tolfreport,tolxreport;
  FTYPE normf;
#if(DEBUG)
  calls++;
#endif

  tolfallowed=tolxallowed=1E-4;
  tolfreport=tolxreport=tolf*1.E3;

  // for debug purposes
  mnewtstepfail=9198;
  mnewtifail=0;
  mnewtjfail=63;
  mnewtpartialstepfail=0;

  // settings
  DODAMP=2;
  dampfactor=1.0;
  dampfactorchange=0.5;
  donesincechange=0;
  dampdeath=0;
  // for DODAMP==2
  allownewdamp=0; // start fresh
  numstabletot=5;
  numdampedtot=5;
  countstable=0;
  countdamped=0;

#if(DEBUGPOINT)
  if((myid==2)&&(icurr+startpos[1]==mnewtifail)&&(jcurr+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
    
    DODAMP=2;
    
    // testing, works well, turn on the if just inside the ntrial loop
    //        DODAMP=1;
  }
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

  for (k = 1; k <= ntrial; k++) {
    truetrialnum++;
    nstroke++;

    // determine error in f, conserved quantities
    usrfunreturn=usrfun(x, n, fvec, fjac,&normf);
#if(DEBUG)
    for(i=1;i<=n;i++){ startfvec[i]=fvec[i]; }
#endif

    if (usrfunreturn>=1){
            
      // fix up the damping if we get a psychotic solution
      errf=1.E30;
      errx=1.E30;
      failed=0; // force no failure condition
    }
    else{ // estimate error normally
      errf = 0.0;
      for (i = 1; i <= n; i++)	errf += fabs(fvec[i]);
      errf/=normf; // renormalize since truly U (conserved quantity) has no better significance than this error

      if (errf <= abstol) lowtol[1]=2;
      else if (errf <= tolf) lowtol[1]=1;
      else lowtol[1]=0;
    }

    // now fix x (primitive variables)

    if(((DODAMP==1)||(DODAMP==2)&&(allownewdamp))&&(k>1)&&(lasterrf<errf)){ // only consider if damped failed when using damping odd/even trial
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
	for (i = 1; i <= n; i++)	x[i]=xold[i];
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
      
      for (i = 1; i <= n; i++)	pp[i] = -fvec[i];
      ludcmp(fjac, n, indx, &d);
      lubksb(fjac, n, indx, pp);
      // DAMP (faster to damp every other one)
      // 
      if(DODAMP) for (i = 1; i <= n; i++)	pp[i] = dampfactor*pp[i];

      ///////
      // evaluate x error (already properly normalized since directly what we seek)
      ///////
      errx = 0.0;
      for (i = 1; i <= n; i++) {
	errx += fabs(pp[i]);
	x[i] += pp[i];
      }

      if (errx <= abstol) lowtol[0]=2;
      else if (errx <= tolx) lowtol[0]=1;
      else lowtol[0]=0;

    }// else if error decreased
#if(DEBUGPOINT)
    //if((myid==2)&&(icurr+startpos[1]==mnewtifail)&&(jcurr+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
      //	if (errf <=1E-20){
      //for (i = 1; i <= n; i++) dualfprintf(fail_file,"wtf: true=%d k=%d errx=%21.15g pp[%d]=%25.17g\n",truetrialnum, k, errx,i,pp[i]);
      //	}
    //}

#endif
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
    if((lowtol[0]==2)||(lowtol[1]==2)) break; // end immediately since we are unable to go further anyways below machine precision
    else if((lowtol[0]==1)&&(lowtol[1]==1)) break; // then exactly what we wanted      
    // if 0 0 or 1 0 or 0 1, then continue to try to find better solution

  }// over trials


  ///////////////////////////////
  //
  // done with MNEWT, now debug stuff comes
  //
  ///////////////////////////////

  // see what's going on
#if(DEBUGPOINT)
  if((myid==2)&&(icurr+startpos[1]==mnewtifail)&&(jcurr+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
    dualfprintf(fail_file,"trueerrx=%25.17g trueerrf=%25.17g\n",errx,errf);
    errf=1E30; errx=1E30; lowtol[1]=0; lowtol[0]=0; // pretend    
  }
#endif
  // some counting on this run of mnewt, failed or not
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
    if(debugfail>=2){
      dualfprintf(fail_file,"proc: %d, t=%25.17g realnstep=%ld\nmnewt didn't converge (k=%d true=%d): i=%d j=%d, errf: %g errx: %g dampfactor: %g\n",myid,t,realnstep, k, truetrialnum, startpos[1]+icurr,startpos[2]+jcurr,errf,errx,dampfactor);
    }
  }
  // assume won't fail if not too bad convergence if <=1E-4
  if ((errf <= tolfallowed)&&(errx<=tolxallowed)) {
    return (0); // for now
  }
  else{ // if >1E-4, then something is wrong
#if(DEBUG)
    if(debugfail>=1){
      dualfprintf(fail_file, "mnewt: (k=%d truetrialnum=%d) failure\n", k,truetrialnum);
    }
    if(debugfail>=2){
      // COMMENT
      // if failed, probably went outside allowed solution given constraints
      // placed on p.  How can this occur when U and p are related?  Is U somehow
      // more flexible?
      coord(icurr,jcurr,CENT,X);
      bl_coord(X,&r,&th);
      get_geometry(icurr,jcurr,CENT,&geom);
      dualfprintf(fail_file,"i=%d j=%d \nx1=%25.17g x2=%25.17g \nr=%25.17g th=%25.17g \ng=%25.17g\n",startpos[1]+icurr,startpos[2]+jcurr,X[1],X[2],r,th,geom.g);
      dualfprintf(fail_file,"x->%25.17g``20, y->%25.17g``20\n",X[1],X[2]);
      for(i=1;i<=n;i++) { dualfprintf(fail_file, "startfvec[%d]=%25.17g startx[%d]=%25.17g x[%d]=%25.17g\n",i,startfvec[i],i,startx[i],i,x[i]); }
      for(i=1;i<=n;i++){
	for(j=1;j<=n;j++){
	  dualfprintf(fail_file, "fjac[%d][%d]=%25.17g ",i,j,fjac[i][j]);
	}
	dualfprintf(fail_file,"\n");
      }
#endif      
#if(DEBUG==2)
#if(0)
      dualfprintf(fail_file,"mnewt={");
      for(i=0;i<n;i++){ // seperately for each primitive variable
	dualfprintf(fail_file,"{");
	for(j=0;j<truetrialnum+1;j++){	  
	  dualfprintf(fail_file,"%25.17g``20 ",trialvalue[j][i]);
	  if(j<truetrialnum-1) 	    dualfprintf(fail_file,",");	    
	}
	dualfprintf(fail_file,"}\n");
	if(i<n-1) 	    dualfprintf(fail_file,",");	    
      }
      dualfprintf(fail_file,"};\n");
      dualfprintf(fail_file,"mnewterr={");
      for(i=0;i<=1;i++){ // seperately for each primitive variable
	dualfprintf(fail_file,"{");
	for(j=0;j<truetrialnum;j++){	  
	  dualfprintf(fail_file,"%25.17g``20 ",trialerr[j][i]);
	  if(j<truetrialnum-1) 	    dualfprintf(fail_file,",");	    
	}
	dualfprintf(fail_file,"}\n");
	if(i<1) 	    dualfprintf(fail_file,",");	    
      }
      dualfprintf(fail_file,"};\n");
#else
      out=fopen("mnewtvaluelist.txt","wt");
      if(out==NULL){ fprintf(stderr,"cannot open mnewtvaluelist.txt\n"); exit(1);}
      for(j=0;j<truetrialnum+1;j++){
	for(i=0;i<n;i++) {
	  fprintf(out,"%25.17g ",trialvalue[j][i]);
	}
	fprintf(out,"\n");
      }
      fclose(out);

      out=fopen("mnewterrlist.txt","wt");
      if(out==NULL){ fprintf(stderr,"cannot open mnewterrlist.txt\n"); exit(1);}
      for(j=0;j<truetrialnum;j++){
	for(i=0;i<2;i++) {
	  fprintf(out,"%25.17g ",trialerr[j][i]);
	}
	fprintf(out,"\n");
      }
      fclose(out);
#endif      
#endif
      //failed = 3;		// source of failure (nonconvergence)
    }
    FAILSTATEMENT("mnewt.c", "convergence", 1);
  }
}
