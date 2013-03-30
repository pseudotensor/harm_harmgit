#include "decs.h"


#define DEBUG (0)
#define DEBUGPOINT (0)
#define MAXTRIAL 1000

// whether to check if UU or RHO go negative in damping routine and damp if so.
// problems when they are negative at start?
#define CHECKUURHO 0

// only p (x) is allowed to be offset.  Rest are forced to be 1..n
#define PITERLOOP(i) for(i = startp+1; i <= startp+n; i++)

#define P0ITERLOOP(i) for(i = 1; i <= n; i++)


int mnewt(FTYPE *U_target,FTYPE *pr0,int numnormterms,int whichcons, int primtoUcons, struct of_geom *ptrgeom, int whethertoreport,int ntrial, int mintrial, FTYPE x[], int n, int startp, FTYPE tolx, FTYPE tolf, FTYPE tolxallowed, FTYPE tolfallowed, FTYPE tolxreport, FTYPE tolfreport, struct of_newtonstats *newtonstats)
{
  extern int usrfun(FTYPE *U_target,FTYPE *pr0,int numnormterms,int whichcons, int primtoUcons, struct of_geom *ptrgeom, FTYPE *prguess, int n, FTYPE *beta, FTYPE **alpha,FTYPE*norm);
  int i = 0, j = 0, k = 0;
  FTYPE errx=1E30, errf=1E30, d;
  FTYPE lasterrx,lasterrf;
  FTYPE dampfactor,dampfactorchange;
  FTYPE redampfactor,redampfactorrho,redampfactoru;
  int redamp;
  // debug stuff
  static long count = 0;
  static long lastnstep = 0;
  static long calls = 0;
  static FTYPE *startx,*startfvec;
  //
  int usrfunreturn;
  struct of_geom geom;
  int lowtol[2]={0,0}; // 0=errx, 1=errf
  FTYPE X[NDIM],V[NDIM];
#if(DEBUG==2)
  FTYPE trialvalue[MAXTRIAL][NPR];
  FTYPE trialerr[MAXTRIAL][2];
  FILE  * out;
#endif
  int truetrialnum=0;
  int donesincechange;
  int DODAMP;
  FTYPE abstol=(NUMEPSILON*50.0); // near machine precision, typically best mnewt() can do
  FTYPE abstolhard=(NUMEPSILON); // stop trying to iterate below this
  int allownewdamp,numstabletot,countstable,numdampedtot,countdamped;
  int dampdeath;
  long int mnewtstepfail;
  int mnewtifail,mnewtjfail,mnewtpartialstepfail;
  int shouldleave;
  FTYPE norm;
#if(DEBUG)
  calls++;
#endif




#if(USEOPENMP)
  // maintain thread safety
  int *indx;
  FTYPE **fjac, *fvec, *pp,*xold,*testx,*pptest;
  indx = ivector(1, n);
  pp = dvector(1, n);
  pptest = dvector(1, n);
  xold = dvector(1, n);
  testx = dvector(1, n);
  startx = dvector(1, n);
  startfvec = dvector(1, n);
  fvec = dvector(1, n);
  fjac = dmatrix(1, n, 1, n);
#else
  static int firstc = 1;
  static int *indx;
  static FTYPE **fjac, *fvec, *pp,*xold,*testx,*pptest;

  if (firstc) {
    firstc = 0;

    indx = ivector(1, n);
    pp = dvector(1, n);
    pptest = dvector(1, n);
    xold = dvector(1, n);
    testx = dvector(1, n);
    startx = dvector(1, n);
    startfvec = dvector(1, n);
    fvec = dvector(1, n);
    fjac = dmatrix(1, n, 1, n);
  }
#endif


  if(tolx<abstol) tolx=abstol;
  if(tolf<abstol) tolf=abstol;


  // for debug purposes
  mnewtstepfail=11;
  mnewtifail=0;
  mnewtjfail=0;
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
  //if((myid==2)&&(ptrgeom->i+startpos[1]==mnewtifail)&&(ptrgeom->j+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
    
  DODAMP=2;
    
  // testing, works well, turn on the if just inside the ntrial loop
  //        DODAMP=1;
}
#endif


// save guess for debugging
#if(DEBUG)

PITERLOOP(i){
  startx[i]=x[i];
#if(DEBUG==2)
  trialvalue[truetrialnum][i-1]=x[i];
#endif
}

#endif
  
// initially
PITERLOOP(i) xold[i]=x[i];




             ///////////////////////
             //
             // DO MNEWT LOOP
             //
             newtonstats->lntries=0;

             //
             for (k = 1; k <= ntrial; k++) {
               // count how many iterations
               truetrialnum++;
               (newtonstats->nstroke)++;
               (newtonstats->lntries)++;

    
               usrfunreturn=usrfun(U_target,pr0,numnormterms,whichcons, primtoUcons,ptrgeom, x, n, fvec, fjac,&norm);


               // DEBUG:
               //    P0ITERLOOP(i) dualfprintf(fail_file,"usrfunreturn=%d k=%d x[%d]=%21.15g\n",usrfunreturn,k,i,x[i]);
               //    P0ITERLOOP(i)P0ITERLOOP(j) dualfprintf(fail_file,"fjac[%d][%d]=%21.15g norm=%21.15g\n",i,j,fjac[i][j],norm);




#if(0)
               // DEBUG: OPENMPMARK: Not thread safe
               // determine error in f, conserved quantities
               P0ITERLOOP(i){
                 if(!finite(x[i])) return(1); // just fail if nan'ed or inf'ed out
                 if(!finite(fvec[i])) return(1); // just fail if nan'ed or inf'ed out
                 P0ITERLOOP(j) if(!finite(fjac[i][j])) return(1); // just fail if nan'ed or inf'ed out
               }
#endif


#if(DEBUG)
               P0ITERLOOP(i){ startfvec[i]=fvec[i]; }
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
                 P0ITERLOOP(i) errf += fabs(fvec[i]);
                 errf/=norm; // renormalize since truly U (conserved quantity) has no better significance than this error

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
                   if(debugfail>=1 && whethertoreport){
                     dualfprintf(fail_file,"mnewt: ok, we really need something better\n: dampfactor: %21.15g errf: %21.15g\n",dampfactor,errf); 
                   }
                   dampdeath=1; // let the bottom non-convergence criterea handle things
                 }
                 if(dampdeath==0){
                   PITERLOOP(i) x[i]=xold[i];
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
                 PITERLOOP(i) xold[i]=x[i]; 
      

                 // here we evaluate J^i_j \delta x^j = -F^i where J is the Jacobian dF^i/dx_j , x^j are the primitives, and F^i are the difference between the target and current conserved quantities (i.e. we are finding the zeroes)
                 // fjac and fvec are already normalized, so pp result is correctly normalized
                 P0ITERLOOP(i) pp[i] = -fvec[i];
                 // fjac here determines pp
                 ludcmp(fjac, n, indx, &d);
                 lubksb(fjac, n, indx, pp);


                 // DEBUG:
                 //      P0ITERLOOP(i)P0ITERLOOP(j) dualfprintf(fail_file,"ifjac[%d][%d]=%21.15g pp[%d]=%21.15g fvec[%d]=%21.15g\n",i,j,fjac[i][j],i,pp[i],i,fvec[i]);



                 // DAMP (faster to damp every other one)
                 // 
                 if(DODAMP){
                   P0ITERLOOP(i) pptest[i] = dampfactor*pp[i];

#if(CHECKUURHO)
                   // check if constraints on some quantities are violated.
                   // particular to inversion doing
                   testx[RHO+1] = x[RHO+1]+pptest[RHO+1];
                   testx[UU+1] = x[UU+1]+pptest[UU+1];
                   redamp=0;
                   redampfactor=redampfactorrho=redampfactoru=dampfactor;

                   if(testx[RHO+1]<0){ redampfactorrho=-x[RHO+1]/pp[RHO+1]*0.1; redamp=1;}
                   if(testx[UU+1]<0){ redampfactoru=-x[UU+1]/pp[UU+1]*0.1; redamp=1;}
                   if(redamp){
                     if(redampfactorrho<redampfactoru) redampfactor=redampfactorrho;
                     else  redampfactor=redampfactoru;

                     P0ITERLOOP(i) pptest[i] = redampfactor*pp[i];
                     if(redampfactor<1E-7){
                       if(debugfail>=1){
                         dualfprintf(fail_file,"mnewt: redampfactor death : ok, we really need something better\n: redampfactor: %21.15g errf: %21.15g\n",redampfactor,errf);
                       }
                       dampdeath=1; // let the bottom non-convergence criterea handle things
                       break;
                     }
                   }
#endif
                   P0ITERLOOP(i) pp[i] = pptest[i];
                 }


                 ///////
                 // evaluate x error (already properly normalized since directly what we seek)
                 ///////
                 errx = 0.0;
                 P0ITERLOOP(i){
                   errx += fabs(pp[i]);
                   // one location where x and another thing are coupled, so shift x manually
                   x[startp+i] += pp[i]; 

                   // DEBUG:
                   // dualfprintf(fail_file,"pp[%d]=%21.15g\n",i,pp[i]);


                 }

                 if (errx <= abstolhard) lowtol[0]=2;
                 else if (errx <= tolx) lowtol[0]=1;
                 else lowtol[0]=0;

               }// else if error decreased



#if(DEBUGPOINT)
               //if((myid==2)&&(ptrgeom->i+startpos[1]==mnewtifail)&&(ptrgeom->j+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
               // if (errf <=1E-20){
               //P0ITERLOOP(i) dualfprintf(fail_file,"wtf: true=%d k=%d errx=%21.15g pp[%d]=%21.15g\n",truetrialnum, k, errx,i,pp[i]);
               // }
               //}

#endif
#if(DEBUG==2)
               // some debug stuff, done every trial type
               PITERLOOP(i) {
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
#if(DEBUGPOINT)
//  if((myid==2)&&(ptrgeom->i+startpos[1]==mnewtifail)&&(ptrgeom->j+startpos[2]==mnewtjfail)&&(realnstep==mnewtstepfail)&&(partialstep==mnewtpartialstepfail)){
//    dualfprintf(fail_file,"trueerrx=%21.15g trueerrf=%21.15g\n",errx,errf);
//    errf=1E30; errx=1E30; lowtol[1]=0; lowtol[0]=0; // pretend    
//  }
#endif


// some counting on this run of mnewt, failed or not
#if(DEBUG)
if (lastnstep < nstep) {
  logfprintf("#1 count/zone: %21.15g calls: %21.15g\n",
             ((FTYPE) count) / ((FTYPE) (N1 * N2)),
             ((FTYPE)calls) / ((FTYPE)(N1 * N2)));
  logfprintf("count: %ld zones: %d calls: %ld\n",
             count,N1 * N2,calls);

  if(0){// can't assume all CPUs get here, since may use different inversions
    mpildsum0(&count,0);
    mpildsum0(&calls,0);
    /*
      myfprintf(stderr,"count: %ld zones: %d calls: %d\n",
      count,totalzones,calls);
    */
    myfprintf(logfull_file,"#1 count/zone: %21.15g calls: %21.15g\n",
              ((FTYPE) count) / ((FTYPE) (totalzones)),((FTYPE)
                                                        calls) / ((FTYPE)(totalzones)));
    myfprintf(stderr,"#1 count/zone: %21.15g calls: %21.15g\n",
              ((FTYPE) count) / ((FTYPE) (totalzones)),((FTYPE)
                                                        calls) / ((FTYPE)(totalzones)));
  }

  count = (long)(k - 1);
  lastnstep = nstep;
  calls = 0;
 } else {
  count += (long)(k - 1);
 }
#endif

// determine if can leave or not
if((lowtol[0]>=1)&&(lowtol[1]>=1)) shouldleave=1;
// otherwise we didn't reach tolerances we wanted note that both
// tolerances are asked for, not just one or the other, since one
// variable may have low error but the other high error, and that's
// not what we want.  We want both to be low.





if(!shouldleave){
  ////////////////////////////////////
  //
  // rest of this is error analysis (i.e. failure, or acceptable failure)
  //
  //


  // if we got here, we never converged to desired tolerance in the specified maximum number of trials

  // only report if pseudo-bad convergence i.e. not near limits since that produces too much data
  if ((errf >= tolfreport)||(errx >= tolxreport)) {
    if(debugfail>=2){
      dualfprintf(fail_file,"proc: %d, t=%21.15g realnstep=%ld\nmnewt didn't converge (k=%d true=%d): i=%d j=%d k=%d, errf: %21.15g errx: %21.15g dampfactor: %21.15g\n",myid,t,realnstep, k, truetrialnum, startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,errf,errx,dampfactor);
    }
  }
  // assume won't fail if not too bad convergence
  if ((errf <= tolfallowed)&&(errx<=tolxallowed)) {
    shouldleave=1; // for now
  }
  else{



    /////////////
    //
    // if >1E-4, then something is very wrong
    //
    ///////////



#if(DEBUG)
    if(debugfail>=1){
      dualfprintf(fail_file, "mnewt: (k=%d truetrialnum=%d) failure\n", k,truetrialnum);
    }
    if(debugfail>=2){
      // COMMENT
      // if failed, probably went outside allowed solution given constraints
      // placed on p.  How can this occur when U and p are related?  Is U somehow
      // more flexible?
      coord(ptrgeom->i,ptrgeom->j,ptrgeom->k,CENT,X);
      bl_coord(X,V);
      dualfprintf(fail_file,"i=%d j=%d k=%d \nx1=%21.15g x2=%21.15g x3=%21.15g \nV1=%21.15g V2=%21.15g V3=%21.15g \ng=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,X[1],X[2],X[3],V[1],V[2],V[3],ptrgeom->g);
      dualfprintf(fail_file,"X1->%21.15g``20, X2->%21.15g``20, X3->%21.15g``20\n",X[1],X[2],X[3]);
      P0ITERLOOP(i) { dualfprintf(fail_file, "startfvec[%d]=%21.15g startx[%d]=%21.15g x[%d]=%21.15g\n",i,startfvec[i],i,startx[startp+i],i,x[startp+i]); }
      P0ITERLOOP(i){
        P0ITERLOOP(j){
          dualfprintf(fail_file, "fjac[%d][%d]=%21.15g ",i,j,fjac[i][j]);
        }
        dualfprintf(fail_file,"\n");
      }
#if(DEBUG==2)
#if(0)
      dualfprintf(fail_file,"mnewt={");
      for(i=0;i<n;i++){ // seperately for each primitive variable
        dualfprintf(fail_file,"{");
        for(j=0;j<truetrialnum+1;j++){   
          dualfprintf(fail_file,"%21.15g``20 ",trialvalue[j][i]);
          if(j<truetrialnum-1)      dualfprintf(fail_file,",");     
        }
        dualfprintf(fail_file,"}\n");
        if(i<n-1)      dualfprintf(fail_file,",");     
      }
      dualfprintf(fail_file,"};\n");
      dualfprintf(fail_file,"mnewterr={");
      for(i=0;i<=1;i++){ // seperately for each primitive variable
        dualfprintf(fail_file,"{");
        for(j=0;j<truetrialnum;j++){   
          dualfprintf(fail_file,"%21.15g``20 ",trialerr[j][i]);
          if(j<truetrialnum-1)      dualfprintf(fail_file,",");     
        }
        dualfprintf(fail_file,"}\n");
        if(i<1)      dualfprintf(fail_file,",");     
      }
      dualfprintf(fail_file,"};\n");
#else
      out=fopen("mnewtvaluelist.txt","wt");
      if(out==NULL){ stderrfprintf("cannot open mnewtvaluelist.txt\n"); exit(1);}
      for(j=0;j<truetrialnum+1;j++){
        for(i=0;i<n;i++) {
          fprintf(out,"%21.15g ",trialvalue[j][i]);
        }
        fprintf(out,"\n");
      }
      fclose(out);

      out=fopen("mnewterrlist.txt","wt");
      if(out==NULL){ stderrfprintf("cannot open mnewterrlist.txt\n"); exit(1);}
      for(j=0;j<truetrialnum;j++){
        for(i=0;i<2;i++) {
          fprintf(out,"%21.15g ",trialerr[j][i]);
        }
        fprintf(out,"\n");
      }
      fclose(out);
#endif      
#endif
      //failed = 3;  // source of failure (nonconvergence)
    }
#endif

    if(whethertoreport){
      FAILSTATEMENT("mnewt.c", "convergence", 1);
    }
    else shouldleave=2;
  }
 }



#if(USEOPENMP)
// maintain thread safety
free_ivector(indx,1, n);
free_dvector(pp,1, n);
free_dvector(pptest,1, n);
free_dvector(xold,1, n);
free_dvector(testx,1, n);
free_dvector(startx,1, n);
free_dvector(startfvec,1, n);
free_dvector(fvec,1, n);
free_dmatrix(fjac,1, n, 1, n);
#endif

if(shouldleave==1) return(0);
 else return(1);

}
