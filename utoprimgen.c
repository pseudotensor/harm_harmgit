


#include "decs.h"

// Note that if using primtoU() or primtoflux(), these respect nprlist for PLOOP, hence why these tests use PLOOP while others use PALLLOOP



// fractional error above which inversion is reported on as having made a signficant error
#if(PRODUCTION)
#define CHECKONINVFRAC (1E-1)
#else
#define CHECKONINVFRAC (1E-2)
#endif

// whether to fail if check on inversion fails
#define FAILIFBADCHECK 1

// fraction upon if greater will treat as failure if FAILIFBADCHECK is triggered
#define CHECKONINVFRACFAIL (1E-1)

static int negdensitycheck(FTYPE *prim, PFTYPE *pflag);
static int check_on_inversion(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats);
static int debug_utoprimgen(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew);
static int compare_ffde_inversions(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats);

static int trycoldinversion(PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats);


int Utoprimgen(int finalstep, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  int otherfail;
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  extern int Utoprim_coldgrmhd(FTYPE *U, struct of_geom *geom, FTYPE *pr, int *positivityproblem);
  int positivityproblem;
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
  int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
  int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats);
  int Utoprimgen_tryagain2(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats);
  void convert_U_removerestmassfromuu(int utoprimverison, int removerestmassfromuu, FTYPE *U);
  int invert_scalars(FTYPE *Uold, FTYPE *Ugeomfree0,FTYPE *Ugeomfree,FTYPE *pr0,FTYPE *pr);
  int pl,pliter;
  PFTYPE lpflag;
  
  
  ////////////////////////
  //
  // INVERT Magnetic field or restore magnetic field for SPLITNPR case
  //
  ////////////////////////

#if(SPLITNPR)
  if(nprlist[nprstart]==B1 && nprlist[nprend]==B3){ // then assume doing B1,B2,B3
    // invert magnetic field only
    PLOOPBONLY(pl) pr[pl]=U[pl]/ptrgeom->e[pl]; // direct inversion from geometry-free conserved quantity
    return(0);
  }
  else{
    // otherwise assume want to invert ONLY non-B components and assume using same memory space
    // can achieve this by setting Ugeomfree to *existing* pr (assume memory space for SPLITNPR step 0 and 1 didn't change
    // overwrites update (that is no field update done on second pass) to magnetic field
    // Note that field U not updated in this pass, so need to fill with something.
    // Assume pr here is pbb that contains field-update's Bnew
    PLOOPBONLY(pl) U[pl]=pr[pl]*ptrgeom->e[pl];
  }
#else
  // invert magnetic field
  // direct inversion from geometry-free conserved quantity
  // assume normal inversion takes care of this
  //  PLOOPBONLY(pl) pr[pl]=U[pl]/ptrgeom->e[pl];
#endif






  //  if(U[ENTROPY]<0.0){
  //    dualfprintf(fail_file,"BADENTROPY: i=%d j=%d nstep=%ld steppart=%d : bad U[ENTROPY]=%21.15g\n",ptrgeom->i,ptrgeom->j,nstep,steppart,U[ENTROPY]);
  //  }

  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);
  // e.g. if REMOVERESTMASSFROMUU==1 and inputtype==UEVOLVE, then Ugeomfree has energy T^t_t that includes rest-mass
#if(SPLITNPR)
  if(!(nprlist[nprstart]==B1 && nprlist[nprend]==B3)){
    // then need to update field quantities since UtoU uses PLOOP not PALLLOOP -- GODMARK: could use PALLLOOP in UtoU, but then adds excessive code to phys.c related to SPLITNPR
    PLOOPBONLY(pl) Ugeomfree[pl]=U[pl]/ptrgeom->e[pl];
  }
#endif









  /////////////////////////////////////////////////////////
  //
  // backup
  //
  ////////////////////////////////////////////////////////
  PALLLOOP(k){
    Uold[k]=Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }
  

  /////////////////////////////////////////////////////////////
  //
  // convert U into form appropriate for inversion routine
  //
  /////////////////////////////////////////////////////////////
  convert_U_removerestmassfromuu(UTOPRIMVERSION,REMOVERESTMASSFROMUU,Ugeomfree);
  convert_U_removerestmassfromuu(UTOPRIMVERSION,REMOVERESTMASSFROMUU,Ugeomfree0);





  ////////////////////////
  //
  // INVERT pseudo-passive scalars
  //
  // all pseudo-passive scalars are trivially inverted
  //
  ////////////////////////

  invert_scalars(Uold,Ugeomfree0,Ugeomfree,pr0,pr);






  ////////////////////////
  //
  // DO INVERSION
  //
  ////////////////////////

  // assume solution is good unless flagged as bad
  GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;



  if(EOMTYPE==EOMGRMHD){
    ///////////////////////////////////////////////////
    //
    ///////////// HOT GRMHD
    //
    ///////////////////////////////////////////////////


    if(DOENTROPY!=DOEVOLVEDIRECTENTROPY){
      // then do energy equation version
    
      if(UTOPRIMVERSION!=UTOPRIMCOMPARE) Utoprimgen_pick(UTOPRIMVERSION, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,newtonstats);
      else Utoprimgen_compare(EVOLVENOENTROPY,Ugeomfree,ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,newtonstats);

      // try other methods (assumes all methods can handle WHICHVEL, etc. used for primary model)
      // right now, all can handle WHICHVEL==VELREL4 and energy equation evolution and REMOVERESTMASSFROMUU=0,1
#if((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU<=1)&&(UTOPRIMTRYAGAIN))
      Utoprimgen_tryagain(EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr,newtonstats);
#elif((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU==2)&&(UTOPRIMTRYAGAIN))
      // Can only try again using same type of U since tryagain code doesn't convert U 
      Utoprimgen_tryagain2(EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr,newtonstats);
#endif
    }
    else{ // direct entropy evolution

      // do inversion for entropy version of EOMs
      // only one inversion is setup to handle this
      PALLLOOP(k) prother[k]=pr0[k];
      whichentropy=EVOLVEFULLENTROPY;

      ////////////////////////
      // get entropy evolution (don't use failure -- otherfail)
      MYFUN(Utoprim(whichentropy,Uold, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim", 1);

    } // end if doing direct entropy evolution




    ////////////////////
    //  If hot GRMHD failed or gets suspicious solution, revert to cold GRMHD if solution is cold
    ///////////////////
    if(HOT2COLD){
      int hotpflag;

      // get failure flag
      hotpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);

      if(hotpflag){
	trycoldinversion(hotpflag, pr0, pr, Ugeomfree, Ugeomfree0, ptrgeom,newtonstats);
      }// end if hotpflag

    }




  }
  else if(EOMTYPE==EOMCOLDGRMHD){
    ///////////////////////////////////////////////////
    //
    ///////////// COLDGRMHD
    //
    ///////////////////////////////////////////////////



    if(1){
      // Jon's inversion
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,newtonstats);
    }
    else if(0){ // not working yet
      MYFUN(Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem),"step_ch.c:Utoprimgen()", "Utoprim_coldgrmhd", 1);
      //Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem);
    }



  }
  else if(EOMTYPE==EOMFFDE){

    ///////////////////////////////////////////////////
    //
    ///////////// FORCE FREE
    //
    ///////////////////////////////////////////////////



    // GODMARK: inversions lead to different behavior (start with torus with rho=u=0 but loop of field)!
    
    if(0){ // Jon's old inversion
      MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);
    }
    else if(1){
      // Jon's inversion
      Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,newtonstats);
    }
    else if(0){
      compare_ffde_inversions(&lpflag, pr0, pr, ptrgeom, Ugeomfree0, Ugeomfree, Uold, Unew,newtonstats);
    }


  }



  ///////////////////////////////////////////////////
  //
  ///////////// Report failure
  //
  // complain if pflag set
  // only output if failed
  //
  ///////////////////////////////////////////////////

  ////////////
  //
  // modify failure flag if necessary
  //
  ////////////

#if(EOMTYPE == EOMGRMHD)
  if(finalstep){
    negdensitycheck(pr, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL));
  }
#endif

  // for now only report if not just negative density failure
  lpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);


  if(IFUTOPRIMFAILSOFT(lpflag)){
    // then don't report info for now SUPERGODMARK
  }
  else if(lpflag &&(debugfail>=1)){
    dualfprintf(fail_file, "Failed to find a prim. var. solution!! t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d : errx=%21.15g\n",t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag,newtonstats->lerrx);
  }





#if(PRODUCTION==0)
  ///////////////////////////////////////////////////
  //
  ///////////// LOTS OF DEBUG STUFF
  //
  ///////////////////////////////////////////////////
  debug_utoprimgen(&lpflag, pr0, pr, ptrgeom, Uold, Unew);



  //////////////////////////////////////
  //
  // check up on solution from inversion if PRODUCTION==0
  // Assume if PRODUCTION==1 that inversion always does good job
  //
  ////////////////////////////////////////
  check_on_inversion(&lpflag, pr0, pr, ptrgeom, Uold, Unew,newtonstats);
#endif


     
  return(0);

}






// use cold grmhd if hot grmhd gives u<0 since then cold is valid
#define USECOLDIFHOTUNEG 1
// use cold grmhd if hot grmhd leads to convergence or other "no solution" type failure
#define USECOLDIFHOTFAILCONV 1

// try cold inversion of hot one fails
int trycoldinversion(PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats)
{
  int k;
  FTYPE prhot[NPR],prcold[NPR];
  PFTYPE coldpflag;


  //  dualfprintf(fail_file,"Got here in trycoldinversion\n");

  if(IFUTOPRIMFAILSOFTRHORELATED(hotpflag) || (USECOLDIFHOTUNEG==0 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) ){
    // then maybe not so bad failure
    
  }
  else if( (USECOLDIFHOTUNEG==1 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) || (USECOLDIFHOTFAILCONV==1 && hotpflag!=0) ){
    
    
    
    // then bad failure, so try to use cold grmhd
    // restore backup in case previous inversion changed things
    PALLLOOP(k){
      prhot[k]=pr[k];
      Ugeomfree[k]=Ugeomfree0[k];
      pr[k]=pr0[k];
    }
    
    
    // get cold inversion
    MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(EOMCOLDGRMHD,GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),Ugeomfree, ptrgeom, &coldpflag, prcold,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);
    
    
    ///////////////////////////////////
    //
    // check if cold solution is good
    if(coldpflag==UTOPRIMNOFAIL){

      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL; // default then is that no failure

      // then keep cold solution
      PALLLOOP(k) pr[k]=prcold[k];


      // but set internal energy to previous value (should really evolve with entropy equation, but if negligible and no strong shocks, then ok )
      // GODMARK: another ok way is to set special failure that only averages internal energy!  Then can evolve at least -- in some diffusive way
      
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and good! hotpflag=%d coldpflag=%d\n",hotpflag,coldpflag);

      ///////////////////////////////
      //
      // decide how to use cold inversion solution    
      if(IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)){
	// since cold approximation is very good, then use cold solution and just set internal energy to 0
	// if internal energy is actually small, then just set it to 0
	// works for Hubble flow!
	pr[UU]=0.0;
      }
      else{
	//////////////
	//  if internal energy is not negligible or unknown, then should average or evolve!
	pr[UU]=pr0[UU];
	// then very bad failure, so try cold inversion and average internal energy for now
	GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMFAILU2AVG1; // assume only internal energy needs correcting by averaging
	
      }// end if (else) trying cold inversion

    }// else if coldpflag is bad
    else{
      // then both hot and cold are bad, so keep hot
      // GODMARK: Could try going to force-free and "failing" the parallel velocity so it gets averaged like internal energy in cold case!'
      // only can go to force-free if b^2/rho>>1 as well
      // keep hotpflag and keep hot solution
      PALLLOOP(k) pr[k]=prhot[k];
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and bad! hotpflag=%d coldpflag=%d\n",hotpflag,coldpflag);
    }
        
  }// if bad failure of some kind

  return(0);
}




// function to check on success of inversion by comparing original conserved quantities (Uold) with new conserved quantities (Unew(p(Uold))) as computed from primitives (p(Uold)) obtained from the inversion.
static int check_on_inversion(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats)
{
  int badinversion,badinversionfail;
  FTYPE Unormalnew[NPR],Unormalold[NPR];
  FTYPE fdiff[NPR];
  struct of_state q;
  int pl,pliter;
  int j,k,jj;
  FTYPE errornorm;




#if(CHECKONINVERSION && (EOMTYPE!=EOMFFDE)) // in force-free projection occurs that makes momenta terms not the same

  //  if(1 || !(*lpflag)){ // only report if not failure
  if(!(*lpflag)){ // only report if not failure
    MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included


    // assume not bad
    badinversion=0;
    badinversionfail=0;


    /////////////
    //
    // first create a normalized (geometrically) conserved quantity
    //
    /////////////
    // assume default value assignment for normalized version of U
    PLOOP(pliter,pl){
      Unormalold[pl] = Uold[pl];
      Unormalnew[pl] = Unew[pl];
    }      

#if(REMOVERESTMASSFROMUU==2)
    // now make U_0 and U_i comparable
    // like \rho\gamma^2 v^2/(1+\gamma) + (u+p)\gamma^2 \sim \rho \gamma v^2 + (u+p)\gamma^2
    // In the limit as v->0, this can be compared against the normalized version of the U_i terms
    Unormalold[UU] = Uold[UU];
    Unormalnew[UU] = Unew[UU];

    // now deal with momentum terms and correct for geometry so all U_i's are comparable
    // (\rho+u+p)\gamma u_i -> (\rho+u+p)\gamma^2 v^2
    SLOOPA(jj) Unormalold[UU+jj] = Uold[UU+jj]*(q.ucon[jj]/q.ucon[TT]);
    SLOOPA(jj) Unormalnew[UU+jj] = Unew[UU+jj]*(q.ucon[jj]/q.ucon[TT]);
#endif
    //////////////////////////////
    //
    // now check for errors
    //
    //////////////////////////////
    PLOOP(pliter,pl){

      // default is to assume nothing wrong
      fdiff[pl]=0.0;

      // no point checking if inversion doesn't handle or is inconsistent with conservation of that quantity
      if(EOMTYPE==EOMFFDE && (pl==RHO || pl==UU) ) continue;
      if(EOMTYPE==EOMCOLDGRMHD  && (pl==UU) ) continue;
      // leave geometry out of it
      //      Unormalnew[pl]*=ptrgeom->gdet;
      //      Unormalold[pl]*=ptrgeom->gdet;
      if(pl==RHO || pl==UU) fdiff[pl] = fabs(Unormalnew[pl]-Unormalold[pl])/(fabs(Unormalnew[pl])+fabs(Unormalold[pl])+SMALL);
      else if(pl==U1 || pl==U2 || pl==U3){

	errornorm  = THIRD*(fabs(Unormalnew[U2])+fabs(Unormalold[U2])+fabs(Unormalnew[U2])+fabs(Unormalold[U2])+fabs(Unormalnew[U3])+fabs(Unormalold[U3]));
#if(REMOVERESTMASSFROMUU==2)
	// only valid comparison if rest-mass taken out of energy term and modify U_i term to be comparable with U_t term
	errornorm = MAX(errornorm,0.5*(fabs(Unormalold[UU])+fabs(Unormalnew[UU])));
#endif
			  
	fdiff[pl] = fabs(Unormalnew[pl]-Unormalold[pl]) / (errornorm+SMALL);

      }
      else if(pl==B1 || pl==B2 || pl==B3){
	fdiff[pl] = fabs(Unormalnew[pl]-Unormalold[pl])/(THIRD*(fabs(Unormalnew[B1])+fabs(Unormalold[B1])+fabs(Unormalnew[B2])+fabs(Unormalold[B2])+fabs(Unormalnew[B3])+fabs(Unormalold[B3]) )+SMALL);
      }
    }
    
    // broke loop to check multiple directions

    PLOOP(pliter,pl){
      if((*lpflag)!=0 || fdiff[pl]>CHECKONINVFRAC){
	if(
	   ( (pl>=RHO)&&(pl<=B3)&&((fabs(Unormalold[pl])>SMALL)&&(fabs(Unormalnew[pl])>SMALL)) )
	    ){
	  badinversion++;
	  dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",pl,fdiff[pl],Unormalold[pl],Unormalnew[pl]);
	}
      }
      if(fdiff[pl]>CHECKONINVFRACFAIL){
	if((pl>=RHO)&&(pl<=B3)&&((fabs(Unormalold[pl])>SMALL)&&(fabs(Unormalnew[pl])>SMALL))){
	  badinversionfail++;
	}
      }
    }


    
    if(FAILIFBADCHECK && badinversionfail){
      (*lpflag)=UTOPRIMFAILCONVBADINVERTCOMPARE;
    }
    if(badinversion){
      dualfprintf(fail_file,"Bad inversion (or possibly Bad U(p) calculation):\n");
      dualfprintf(fail_file,"t=%21.15g nstep=%ld stepart=%d :: i=%d j=%d :: lntries=%d lerrx=%21.15g\n",t,nstep,steppart,ptrgeom->i,ptrgeom->j,newtonstats->lntries,newtonstats->lerrx);
      PLOOP(pliter,pl) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Uold[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g)\n",pl,Uold[pl],pl,Uold[pl]*ptrgeom->gdet,pl,pr[pl],pl,pr0[pl]);
      
      dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->gdet);
      
      DLOOP(j,k) dualfprintf(fail_file,"gcon=%21.15g gcov=%21.15g\n",ptrgeom->gcov[GIND(j,k)],ptrgeom->gcon[GIND(j,k)]);

      for(k=0;k<4;k++){
	dualfprintf(fail_file,"q.ucon[%d]=%21.15g q.ucov[%d]=%21.15g q.bcon[%d]=%21.15g q.bcov[%d]=%21.15g\n",k,q.ucon[k],k,q.ucov[k],k,q.bcon[k],k,q.bcov[k]);
      }

      // only really need the below quantities to check on inversion in mathematica
      // Use Eprime_inversion.nb to check on utoprim_jon.c inversion
      for(j=0;j<NUMINVPROPERTY;j++) dualfprintf(fail_file,"%sstr=\"%21.15g\";\n",newtonstats->invpropertytext[j],newtonstats->invproperty[j]);
    }

  }


#endif


  return(0);


}
			   






// compare ffde inversions
static int compare_ffde_inversions(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats)
{
  int j,k;
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  struct of_state q;
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);


  MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);

  PALLLOOP(k) Ugeomfree[k]=Ugeomfree0[k]; // make sure good conserved quantity
      
  Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), prother, newtonstats);

  // now get conserved quantity from pr to check
  // find U(p)
  MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Upr),"step_ch.c:advance()", "primtoU()", 1);

  MYFUN(get_state(prother, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,prother, &q, ptrgeom, Uprother),"step_ch.c:advance()", "primtoU()", 1);


  // now compare

  //PALLLOOP(k)
  dualfprintf(fail_file,"nstep=%ld stepart=%d i=%d j=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j);
  for(k=U1;k<B3;k++){
    /*
      if(
      ((fabs(Ugeomfree[k]-Upr[k])/(fabs(Ugeomfree[k])+fabs(Upr[k])+SMALL)) > 1E-12)||
      ((fabs(Ugeomfree[k]-Uprother[k])/(fabs(Ugeomfree[k])+fabs(Uprother[k])+SMALL)) > 1E-12)
      ){
      dualfprintf(fail_file,"DIFF: %21.15g %21.15g :: k=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(Ugeomfree[k]-Upr[k])/(fabs(Ugeomfree[k])+fabs(Upr[k])+SMALL)),(fabs(Ugeomfree[k]-Uprother[k])/(fabs(Ugeomfree[k])+fabs(Uprother[k])+SMALL)),k,Ugeomfree[k],pr[k],prother[k],Upr[k],Uprother[k]);
      }
    */
    if(
       ((fabs(pr[k]-prother[k])/(fabs(pr[k])+fabs(prother[k])+SMALL)) > 1E-12)&&( (fabs(pr[k])>1E-15)||(fabs(prother[k])>1E-15) )
       ){
      dualfprintf(fail_file,"DIFF: %21.15g :: k=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(pr[k]-prother[k])/(fabs(pr[k])+fabs(prother[k])+SMALL)),k,Ugeomfree[k],pr[k],prother[k],Upr[k],Uprother[k]);
    }
  }


  return(0);

}




// random debugging code
static int debug_utoprimgen(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew)
{
  int j,k;
  FTYPE fdiff[NPR];






#define CRAZYDEBUG 0
#if(CRAZYDEBUG) // begin crazy debug stuff
  if(1|| newtonstats->lntries>3){
    if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){
    //    if(nstep==0 && steppart==0 && ptrgeom->i==4 && ptrgeom->j==36){
    //  if(nstep==0 && steppart==0 && ptrgeom->i == 5 && ptrgeom->j == 31 ){
    dualfprintf(fail_file,"nstep=%ld stepart=%d :: i=%d j=%d :: lntries=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j,newtonstats->lntries);
    
    //    PALLLOOP(k) dualfprintf(fail_file,"Ugeomfree[%d]=%21.15g pr[%d]=%21.15g\n",k,Ugeomfree[k],k,pr[k]);
    PALLLOOP(k) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Uold[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g)\n",k,Uold[k],k,Uold[k]*ptrgeom->gdet,k,pr[k],k,pr0[k]);

    dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->gdet);

    DLOOP(j,k) dualfprintf(fail_file,"gcon=%21.15g gcov=%21.15g\n",ptrgeom->gcov[GIND(j,k)],ptrgeom->gcon[GIND(j,k)]);
    //    myexit(777);
    //  }
    }
  }


  // test inversion for accuracy in conservative space
  //  if(nstep==0 && steppart==0 && ptrgeom->i == 4 && ptrgeom->j == 36 ){
  if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){

#if(0)
pr[0]= 7.22714301361038e-06 ;
pr[1]=-2.55753797775927e-07 ;
pr[2]= 0.000892168434512681 ;
pr[3]=  -0.0349621334882251 ;
pr[4]=                   -0 ;
pr[5]= -0.00101880685475446 ;
pr[6]=   0.0399035382308458 ;
pr[7]=                    0 ;
#elif(0)
pr[0]= 7.46157819677347e-06;
pr[1]=  7.3407120498644e-06;
pr[2]= 0.000367464781319951;
pr[3]=  -0.0144105143008605;
pr[4]=                    0;
pr[5]= -0.00101880685475446;
pr[6]=   0.0399035382308458;
 pr[7]=                    0;
#elif(0)
pr[0]=  7.5124289176258e-06 ;
pr[1]= 1.33752037209996e-08 ;
pr[2]= 1.33529579432262e-07 ;
pr[3]=-5.23276639757274e-06 ;
pr[4]=                    0 ;
pr[5]= -0.00101880685475444 ;
pr[6]=   0.0399035382308459 ;
 pr[7]=                    0 ;

#endif

  MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included

  for(k=0;k<4;k++){
    dualfprintf(fail_file,"q.ucon[%d]=%21.15g q.ucov[%d]=%21.15g q.bcon[%d]=%21.15g q.bcov[%d]=%21.15g\n",k,q.ucon[k],k,q.ucov[k],k,q.bcon[k],k,q.bcov[k]);
  }

  PLOOP(pliter,k){
    Unew[k]*=ptrgeom->gdet;
    Uold[k]*=ptrgeom->gdet;
    fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k]+Uold[k])+1E-30);
    //    if(fdiff[k]>1E-10){
      if((k>=RHO)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
	dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",k,fdiff[k],Uold[k],Unew[k]);
      }
      //    }
  }

  dualfprintf(fail_file,"dt=%21.15g\n",dt);


  myexit(124);
  }
#endif// end crazy debug stuff
#undef CRAZYDEBUG


  return(0);


}



static int negdensitycheck(FTYPE *prim, PFTYPE *pflag)
{

  //  return(0);

  //Inversion from the average value succeeded or has a negative density or internal energy
  if(IFUTOPRIMFAILSOFT(*pflag)) {
      
    if( prim[UU] < 0.0 ) {
      *pflag = UTOPRIMFAILU2AVG2;
    }
  }

  return(0);
}





// perform U->p inversion on advected scalars
int invert_scalars(FTYPE *Uold, FTYPE *Ugeomfree0,FTYPE *Ugeomfree,FTYPE *pr0,FTYPE *pr)
{
  FTYPE ye;
  FTYPE myrhouu0,oneOmyrhouu0;
  
  // Note that Y_e = Y_l - Y_\nu such that Y_e>0 and Y_e<1 as handled in EOS

  // avoid division by 0 but allow sign and 1/0 ->0 assumed 0 geometry  means value can be taken to be 0
  myrhouu0=Ugeomfree[RHO];
  oneOmyrhouu0=sign(Ugeomfree[RHO])/(fabs(Ugeomfree[RHO])+SMALL);

#if(DOYL!=DONOYL)
  pr0[YL] = pr[YL] = Ugeomfree[YL]*oneOmyrhouu0;
#endif
#if(DOYNU!=DONOYNU)
  pr0[YNU] = pr[YNU] = Ugeomfree[YNU]*oneOmyrhouu0;
#endif

  // Change the primitives to be limited by constraints on Y_e
  // Do this so if out of bounds, Y_e can still recover eventually if fluid goes out of region where Y_\nu is large
#if(DOYL!=DONOYL && DOYNU!=DONOYNU)
  ye = pr[YL] - pr[YNU];

  // have pr[YL] and pr[YNU] share the blame
  if(ye<0.0){
    pr0[YL] = pr[YL]  += -0.5*ye;
    pr0[YNU]= pr[YNU] += +0.5*ye;
    // now effective ye is 0.0
    ye=0.0;
    // fix conserved quantities too
    Uold[YL]=Ugeomfree0[YL]=Ugeomfree[YL] = pr[YL]*myrhouu0;
    Uold[YNU]=Ugeomfree0[YNU]=Ugeomfree[YNU] = pr[YNU]*myrhouu0;
  }

  if(ye>1.0){
    pr0[YL] = pr[YL]  += -0.5*ye+0.5;
    pr0[YNU]= pr[YNU] += +0.5*ye-0.5;
    // Now effective ye is 1.0
    ye=1.0;
    Uold[YL]=Ugeomfree0[YL]=Ugeomfree[YL] = pr[YL]*myrhouu0;
    Uold[YNU]=Ugeomfree0[YNU]=Ugeomfree[YNU] = pr[YNU]*myrhouu0;
  }

  if(pr[YL]<0.0){
    pr0[YL] = pr[YL]  = 0.0;
    pr0[YNU]= pr[YNU] = 0.0;
    // Now effective ye is 0.0
    ye=0.0;
    Uold[YL]=Ugeomfree0[YL]=Ugeomfree[YL] = pr[YL]*myrhouu0;
    Uold[YNU]=Ugeomfree0[YNU]=Ugeomfree[YNU] = pr[YNU]*myrhouu0;
    
  }

  if(pr[YNU]<0.0){
    pr0[YL] = pr[YL]  = ye; // Y_l = Y_e in this case
    pr0[YNU]= pr[YNU] = 0.0;
    Uold[YL]=Ugeomfree0[YL]=Ugeomfree[YL] = pr[YL]*myrhouu0;
    Uold[YNU]=Ugeomfree0[YNU]=Ugeomfree[YNU] = pr[YNU]*myrhouu0;
  }

#endif


  return(0);
}





// convert U (conserved quantities) to form acceptable to utoprimversion for a given removerestmassfromuu
void convert_U_removerestmassfromuu(int utoprimversion, int removerestmassfromuu, FTYPE *U)
{

  ///////////////////////////////////////////////////
  //
  ///////////// Setup U[UU] with or without rest-mass for different inversions and setting of REMOVERESTMASSFROMUU
  //
  // check if inversion is chosen consistent with REMOVERESTMASSFROMUU
  // At this point, Ugeomfree[UU] in "standard form with rest-mass" for REMOVERESTMASSFROMUU=0,1.  If ==2, then no rest-mass
  //
  ///////////////////////////////////////////////////

  if(removerestmassfromuu==2){
    // 5D1 handles removerestmassfromuu internally
    // NONRELCOMPAT only accepts UU without rest-mass, so good to go
    if( utoprimversion!=UTOPRIM5D1 && (utoprimversion!=UTOPRIMJONNONRELCOMPAT) ){
      //    dualfprintf(fail_file,"This code does not handle removerestmassfromuu==2: other Utoprim's besides 5D1 and 2DFINAL\n");
      //    myexit(1);
      // then change conserved quantity so can work
      U[UU]-=U[RHO]; // "add in rest-mass" so in correct form
    }
  }
  else if(removerestmassfromuu==0 || removerestmassfromuu==1){
    if(utoprimversion==UTOPRIMJONNONRELCOMPAT){
      // then change conserved quantity so can work for UTOPRIMJONNONRELCOMPAT
      U[UU]+=U[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
    }
  }


}










// Used for dissipation calculation only
int Utoprimdiss(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, PFTYPE *otherpflag, struct of_newtonstats *newtonstats)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  FTYPE fdiff[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int pl,pliter;

  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);


  /////////////////////////////////////////////////////////
  //
  // backup
  //
  ////////////////////////////////////////////////////////
  PALLLOOP(k){
    Uold[k]=Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }
  

  /////////////////////////////////////////////////////////////
  //
  // convert U into form appropriate for inversion routine (probably does nothing since assume UTOPRIM5D1 used) and feed Utprim Uold anyways
  //
  /////////////////////////////////////////////////////////////
  convert_U_removerestmassfromuu(UTOPRIM5D1,REMOVERESTMASSFROMUU,Ugeomfree);
  convert_U_removerestmassfromuu(UTOPRIM5D1,REMOVERESTMASSFROMUU,Ugeomfree0);


  ////////////////////////
  //
  // DO INVERSION

  // assume solution is good unless flagged as bad
  *otherpflag=UTOPRIMNOFAIL;

  // do separate inversion for entropy version of EOMs
  // only one inversion is setup to handle this
  whichentropy=WHICHENTROPYEVOLVE;

  ////////////////////////
  // get entropy evolution (don't use failure -- otherpflag)
  // this inversion knows about global settings for DOENTROPY and increases tolerance for doing comparison
  MYFUN(Utoprim(whichentropy, Uold, ptrgeom, otherpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim", 1);

  // DEBUG:
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"otherpflag=%d Uold[%d]=%21.15g pr0[%d]=%21.15g pr[%d]=%21.15g\n",*otherpflag,pl,Uold[pl],pl,pr0[pl], pl,pr[pl]);

  
  // now get conserved quantity from pr to check
  // find U(p)
  //  MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  //  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1);
  
  
  //  if(0&& *otherpflag){
  //    PLOOP(pliter,k){
  //  fdiff[k] = fabs(Unew[k]-Uold[k])/(fabs(Unew[k]+Uold[k])+1E-30);
  //      if(fdiff[k]>1E-10){
  //	if((k>=U1)&&(k<=B3)&&((fabs(Uold[k])>1E-20)||(fabs(Unew[k])>1E-20))){
  //	  dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %g %g\n",k,fdiff[k],Uold[k],Unew[k]);
  //	}
  //      }
  //	}
  //      }

     
  return(0);
}






//////////////////////////////
//
// COMPARISON of 2 (currently only 2) methods
int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
  int pl,pliter;
  FTYPE Uo1[NPR], Uo2[NPR];
  FTYPE ptest1[NPR],ptest2[NPR];
  FTYPE Uo[NPR], po[NPR];
  int test;
  PFTYPE lpflag1,lpflag2;
  int lntries1,lntries2;
  FTYPE lerrx1,lerrx2;


  // backup
  PALLLOOP(pl){
    ptest1[pl]=ptest2[pl]=pr[pl];
    Uo1[pl]=Uo2[pl]=Ugeomfree[pl];
  }

  // currently comparing 5D1 and LDZ
  //  Utoprimgen_pick(UTOPRIM5D1, EVOLVENOENTROPY, Uo1, ptrgeom, lpflag, ptest1,newtonstats);
  //  Utoprimgen_pick(UTOPRIMLDZ, EVOLVENOENTROPY, Uo2, ptrgeom, lpflag, ptest2,newtonstats);

  // remove rest-mass
  Uo1[UU]+=Uo1[RHO];
  Utoprimgen_pick(UTOPRIMJONNONRELCOMPAT, EVOLVENOENTROPY, Uo1, ptrgeom, &lpflag1, ptest1,newtonstats);
  lntries1=newtonstats->lntries; newtonstats->lntries=0; lerrx1=newtonstats->lerrx;
  Utoprimgen_pick(UTOPRIM2DFINAL, EVOLVENOENTROPY, Uo2, ptrgeom, &lpflag2, ptest2,newtonstats);
  lntries2=newtonstats->lntries; lerrx2=newtonstats->lerrx;

#define ERRORCONST (1E-10)

  PALLLOOP(pl){
    if(ptest1[pl]!=0.0) test=fabs(ptest1[pl]-ptest2[pl])/ptest1[pl]>ERRORCONST;
    else test=(ptest1[pl]-ptest2[pl])>ERRORCONST;
    if(test){
      dualfprintf(fail_file,"utoprimdiff: %d %ld :: %d %d %d :: %d ::  %21.15g   %21.15g   %21.15g %21.15g :: %d %d :: %21.15g %21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pl,ptest1[pl]-ptest2[pl],(ptest1[pl]!=0.0) ? fabs(ptest1[pl]-ptest2[pl])/ptest1[pl] : ptest1[pl]-ptest2[pl],ptest1[pl],ptest2[pl],lntries1,lntries2,lerrx1,lerrx2);
    }
  }
  if(lpflag1!=0 || lpflag2!=0){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2);
  }
  if(lntries1>10 || lntries2>10){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g lntries1=%d lntries2=%d\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2,lntries1,lntries2);
  }

  
  // always do (use old utoprim)
  if(lpflag1==0) PALLLOOP(pl){
    pr[pl]=ptest1[pl];
    *lpflag=lpflag1;
  }
  else{
    PALLLOOP(pl) pr[pl]=ptest2[pl];
    *lpflag=lpflag2;
  }

  

  return(0);
}




// just picks the algorithm to invert
int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  if(which==UTOPRIM5D1){      MYFUN(Utoprim(parameter,Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim", 1);}
  else if(which==UTOPRIMLDZ){ MYFUN(Utoprim_ldz(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_ldz", 1);}
  else if(which==UTOPRIM2D){ MYFUN(Utoprim_2d(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d", 1);}
  else if(which==UTOPRIM1D){ MYFUN(Utoprim_1d(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d", 1);}
  else if(which==UTOPRIM1DOPT){ MYFUN(Utoprim_1d_opt(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d_opt", 1);}
  else if(which==UTOPRIM1DFINAL){ MYFUN(Utoprim_1d_final(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d_final", 1);}
  else if(which==UTOPRIM2DFINAL){ MYFUN(Utoprim_2d_final(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d_final", 1);}
  else if(which==UTOPRIMJONNONRELCOMPAT){ MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(EOMTYPE,GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);}
  else if(which==UTOPRIM5D2){ MYFUN(Utoprim_5d2_final(Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_5d2_final", 1);}
  
  return(0);
}




// Utoprimgen try again.  Can be whatever sequence or number of inversions
int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats);

  Utoprimgen_tryagain_substep(UTOPRIMJONNONRELCOMPAT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM2DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  // can't trust below really (causes some kind of superslow-down when doing torus problem with  MCSTEEP)
  //  Utoprimgen_tryagain_substep(UTOPRIM5D1, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  // LDZ as normally used generates nan's
  //  Utoprimgen_tryagain_substep(UTOPRIMLDZ, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM1DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM1DOPT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM5D2, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM2D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM1D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);

  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}


// Utoprimgen try again.  Can be whatever sequence or number of inversions
// only those algorithms designed directly for REMOVERESTMASSFROMUU==2
int Utoprimgen_tryagain2(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats);
  void convert_U_removerestmassfromuu(int utoprimverison, int removerestmassfromuu, FTYPE *U);
  FTYPE Uorig0[NPR],Uorig[NPR];
  int pl,pliter;


  Utoprimgen_tryagain_substep(UTOPRIMJONNONRELCOMPAT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);
  Utoprimgen_tryagain_substep(UTOPRIM5D1, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,newtonstats);

  PALLLOOP(pl){
    Uorig0[pl]=Ugeomfree0[pl];
    Uorig[pl]=Ugeomfree[pl];
  }
  // assume at this point that U is such that setup for REMOVERESTMASSFROMUU==2 so that UU has rest-mass subtracted.
  // so need to add it back in for 2D FINAL method
  // convert U into form appropriate for inversion routine 
  convert_U_removerestmassfromuu(UTOPRIM2DFINAL,REMOVERESTMASSFROMUU,Uorig);
  convert_U_removerestmassfromuu(UTOPRIM2DFINAL,REMOVERESTMASSFROMUU,Uorig0);
  // now in form usable by UTOPRIM2DFINAL
  Utoprimgen_tryagain_substep(UTOPRIM2DFINAL, parameter, Uorig0, Uorig, ptrgeom, lpflag, pr0, pr,newtonstats);


  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}





// need to keep this function of to date if going to add other "non-critical" failures.  Point is that some inversions can't handle non-rel case and might get positive results when shouldn't have.
int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int k;

  // if really a bad failure that don't want / can't handle, then try again
  if(
     (*lpflag!=0)
     &&(!((IFUTOPRIMFAILSOFTNOTRHORELATED(*lpflag))&&(STEPOVERNEGU)))
     &&(!((*lpflag==UTOPRIMFAILRHONEG)&&(STEPOVERNEGRHO)))
     &&(!((*lpflag==UTOPRIMFAILRHOUNEG)&&(STEPOVERNEGRHOU)))
     ){
    // restore
    PALLLOOP(k){
      Ugeomfree[k]=Ugeomfree0[k];
      pr[k]=pr0[k];
    }
    // try again
    MYFUN(Utoprimgen_pick(which,parameter,Ugeomfree, ptrgeom, lpflag, pr,newtonstats),"step_ch.c:Utoprimgen_tryagain_substep()", "Utoprimgen_pick", 1);
  }

  return(0);
}




// there may be something wrong with this function -- didn't work in TIMEORDER==4, had to do standard method
// could have just been that I wasn't bounding after using this
int Utoprimloop(FTYPE (*U)[NSTORE2][NSTORE3][NPR],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], struct of_newtonstats *newtonstats)
{
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int i,j,k;

  ZLOOP{
    get_geometry(i, j, k, CENT, ptrgeom);
    // invert True U->p
    MYFUN(Utoprimgen(0,EVOLVEUTOPRIM,UEVOLVE,MAC(U,i,j,k), ptrgeom, MAC(prim,i,j,k),newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
  }
  return(0);
}


int primtoUloop(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q;
  int i,j,k;

  ZLOOP{
    MYFUN(get_state(MAC(prim,i,j,k), ptrgeom, &q),"step_ch.c:primtoUloop()", "get_state()", 1);
    get_geometry(i, j, k, CENT, ptrgeom);
    // forward calculate U(p)
    MYFUN(primtoU(UEVOLVE,MAC(prim,i,j,k), &q, ptrgeom, MAC(U,i,j,k)),"step_ch.c:primtoUloop()", "primtoU()", 1);
  }
  return(0);
}



// filter out velocity along field line
void filterffde(int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  struct of_newtonstats newtonstats;

  // now filter velocity
  get_geometry(i,j,k,CENT,ptrgeom);
  get_state(pr,ptrgeom,&q);
  primtoU(UNOTHING,pr,&q,ptrgeom,U);

  //  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion
  Utoprimgen(0,EVOLVEUTOPRIM,UNOTHING,U,ptrgeom,prout,&newtonstats);

  PALLLOOP(pl) pr[pl]=prout[pl];
  // kill densities
  pr[RHO]=pr[UU]=0.0;



}
