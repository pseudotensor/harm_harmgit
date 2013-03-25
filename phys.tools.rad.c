#include "decs.h"

static int f_implicit_lab(int whichcall, int showmessages, int allowlocalfailurefixandnoreport, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f);


static void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G, FTYPE *chireturn);
static void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu, FTYPE *chireturn);
void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM]);
static void koral_explicit_source_rad(FTYPE *pr, FTYPE *Unew, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *qnew, FTYPE (*dUcomp)[NPR]);
static int koral_explicit_source_rad_prepare(FTYPE *pr, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE *prnew, FTYPE *Unew, struct of_state *qnew);

static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gd, FTYPE chi, struct of_geom *ptrgeom, FTYPE *dtsub);

static int Utoprimgen_failwrapper(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats);

static int simplefast_rad(int dir, struct of_geom *geom,struct of_state *q, FTYPE vrad2,FTYPE *vmin, FTYPE *vmax);


static int opacity_interpolated_urfconrel(FTYPE tautotmax, FTYPE *pp,struct of_geom *ptrgeom,FTYPE *Av, FTYPE Erf,FTYPE gammarel2,FTYPE *Erfnew, FTYPE *urfconrel);

static FTYPE compute_dt(FTYPE *CUf, FTYPE dtin);

static int get_m1closure_gammarel2(int showmessages, struct of_geom *ptrgeom, FTYPE *Av, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn);
static int get_m1closure_gammarel2_cold(int showmessages, struct of_geom *ptrgeom, FTYPE *Av, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel);


static int get_m1closure_Erf(struct of_geom *ptrgeom, FTYPE *Av, FTYPE gammarel2, FTYPE *Erfreturn);

static int get_m1closure_urfconrel(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Av, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad);
static int get_m1closure_urfconrel_olek(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Av, FTYPE gammarel2, FTYPE delta, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad);



// mnemonics for return modes so schemes know how failed and what to do.
#define UTOPRIMGENWRAPPERRETURNFAILRAD 1
#define UTOPRIMGENWRAPPERRETURNFAILMHD 2

// wrapper for Utoprimgen() that returns non-zero if failed in some way so know can't continue with that method
static int Utoprimgen_failwrapper(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats)
{

  //calculating primitives  
  // OPTMARK: Should optimize this to  not try to get down to machine precision
  // KORALTODO: NOTEMARK: If failure, then need to really fix-up or abort this implicit solver!
  // Need to check if Utoprimgen fails in sense of stored inversion failure and revert to explicit or other if happens -- maybe CASE dependent!
  MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, U, ptrgeom, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);

  // check how inversion did.  If didn't succeed, then check if soft failure and pass.  Else if hard failure have to return didn't work.
  int lpflag,lpflagrad;
  lpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
  lpflagrad=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);
  if(IFUTOPRIMFAILSOFT(lpflag)){
    // assume soft failure ok, but reset
	GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
	if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got soft MHD failure inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(IFUTOPRIMRADFAIL(lpflagrad)){
    // have to reduce Newton step if getting failure.
	GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
	if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got some radiation inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    return(UTOPRIMGENWRAPPERRETURNFAILRAD);
  }
  else if( IFUTOPRIMFAIL(lpflag) || IFUTOPRIMRADFAIL(lpflagrad) ){
	// these need to get fixed-up, but can't, so return failure
	if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got hard failure of inversion (MHD part only considered as hard) in f_implicit_lab(): ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
	return(UTOPRIMGENWRAPPERRETURNFAILMHD);
  }
  else{
	// no failure
    //	dualfprintf(fail_file,"No failure in Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  

 

  if(showmessages && debugfail>=2){
	static int maxlntries=0,maxnstroke=0;
	int diff;
	diff=0;
	// For RADSHADOW, gets up to 5
	if(newtonstats->lntries>maxlntries){ maxlntries=newtonstats->lntries; diff=1;}
	if(newtonstats->nstroke>maxnstroke){ maxnstroke=newtonstats->nstroke; diff=1;}
	// only report if grew beyond prior maximum
	if(diff) dualfprintf(fail_file,"newtonsteps: lntries=%d (max=%d) nstroke=%d (max=%d) logerror=%g\n",newtonstats->lntries,maxlntries,newtonstats->nstroke,maxnstroke,newtonstats->lerrx);
  }


  return(0);
}


//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame, returns ultimate deltas ***********************
//******* the fiducial approach *****************************************
//**********************************************************************


// sign of G that goes between Koral determination of G and HARM source term.
#define SIGNGD (-1.0)
//#define SIGNGD 1.0

#define SIGNGD2 (1.0) // sign that goes into implicit differencer

#define SIGNGD3 (1.0) // sign that goes into implicit solver

//uu0 - original cons. qty
//uu -- current iteration
//f - (returned) errors

//NOTE: uu WILL be changed from inside this fcn
static int f_implicit_lab(int whichcall, int showmessages, int allowlocalfailurefixandnoreport, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f)
{
  struct of_state q;
  FTYPE pp[NPR];
  int pliter, pl;
  int iv;
  struct of_newtonstats newtonstats;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.

  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;

  // get primitive (don't change pp0 or uu0)
  PLOOP(pliter,pl) pp[pl] = pp0[pl];

  // get change in conserved quantity between fluid and radiation (equal and opposite 4-force)
  // required for inversion to get P(U) for MHD and RAD variables
  DLOOPA(iv) uu[UU+iv] = uu0[UU+iv] - (uu[URAD0+iv]-uu0[URAD0+iv]);

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"f_implicit_lab: wc=%d i=%d j=%d pl=%d uu=%g\n",whichcall, ptrgeom->i,ptrgeom->j,pl,uu[pl]);
  
  // Get P(U)
  int failreturn;
  failreturn=Utoprimgen_failwrapper(showmessages,allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
  if(failreturn){
	if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of f_implicit_lab()\n");
	return(failreturn);
  }

  // re-get needed q's
  get_state_uconucovonly(pp, ptrgeom, &q);
  get_state_uradconuradcovonly(pp, ptrgeom, &q);

  //radiative covariant four-force
  FTYPE Gd[NDIM];
  FTYPE chireturn;
  calc_Gd(pp, ptrgeom, &q, Gd, &chireturn);

  // compute difference vector between original and new 4-force's effect on conserved radiative quantities
  // i.e. f->0 as change in conserved quantity approaches current value of 4-force
  DLOOPA(iv) f[iv] = (uu[URAD0+iv] - uu0[URAD0+iv]) + (SIGNGD2 * localdt * Gd[iv]);


  //  dualfprintf(fail_file,"i=%d Gd=%g %g %g %g uuG=%g chi=%g\n",ptrgeom->i,Gd[0],Gd[1],Gd[2],Gd[3],(SIGNGD2 * localdt * Gd[iv]),chireturn);


  return 0;
} 


static FTYPE compute_dt(FTYPE *CUf, FTYPE dtin)
{
  // what's applied to source and flux terms to get update (see global.stepch.h and step_ch.c:get_truetime_fluxdt() and step_ch.c:setup_rktimestep()) to get Uf
  // We don't use the ucum update version of dt.  As part of the RK method, the ucum update is separate from the substeps used to get information on updates(?). GODMARK.
  return(CUf[2]*dtin);
}

// compute changes to U (both T and R) using implicit method
// KORALTODO: If doing implicit, should also add geometry source term that can sometimes be stiff.
static int koral_implicit_source_rad(FTYPE *pin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int i1,i2,i3,iv,ii,jj,pliter,sc;
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uu[NPR]; 
  FTYPE uuporig[NPR],uu0orig[NPR];
  FTYPE f1[NDIM],f2[NDIM],f3[NDIM],x[NDIM];
  FTYPE realdt;
  FTYPE radsource[NPR], deltas[NDIM]; 
  int pl;
  static long long int numimplicits=0;
  static long long int numoff1iter=0,numofiter=0;



#if(0)
  // DEBUG
#define RHO_AMB (MPERSUN*MSUN/(LBAR*LBAR*LBAR)) // in grams per cm^3 to match koral's units with rho=1
  dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d\n",nstep,steppart,dt,ptrgeom->i);
  pl=RHO; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]/RHO_AMB,Uiin[pl]/RHO_AMB);
  pl=UU;    jj=0; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]/RHO_AMB,Uiin[pl]/RHO_AMB);
  pl=U1;    jj=1; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
  pl=U2;    jj=2; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
  pl=U3;    jj=3; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
  pl=URAD0; jj=0; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]/RHO_AMB,Uiin[pl]/RHO_AMB);
  pl=URAD1; jj=1; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
  pl=URAD2; jj=2; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
  pl=URAD3; jj=3; dualfprintf(fail_file,"%d %g %g\n",pl,pin[pl]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]))  ,  Uiin[pl]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])));
#endif



  // static counter for diagnosing issues
  numimplicits++;



  // setup implicit loops
  realdt = compute_dt(CUf,dt);
  FTYPE DAMPFACTOR=1.0; // factor by which step Newton's method.
  FTYPE fracdtuu0=1.0,fracdtG=1.0,fracuup=1.0; // initially try full realstep step

 
  if(USEDUINRADUPDATE){
    // uu0 will hold original vector of conserved
    // here original means U[before fluxes, geometry, etc.] + dU[due to fluxes, geometry, etc. already applied and included in dUother]
    // This is required for stiff source term so immediately have balance between fluxes+geometry+radiation.  Otherwise, radiation diffuses.
    // I'm guessing that even though one uses RK2 or RK3, the first step generates large radiative velocities without any balanced source term since U isn't updated yet.  One would hope RK2 would recover on the final substep, but it doesn't!  In RK2, upon the final substep, that velocity is present for the radiation source term.  But it's also present for the fluxes!  That is, if there were no flux update on the final substep, then the source would have balanced the previous flux, but yet another flux is done, so there can be no balance.  This leads to a run-away velocity that would be similar to the \tau\sim 1 case.
    PLOOP(pliter,pl) uu[pl]=uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
    // Note that "q" isn't used in this function or used in function call, so don't have to update it here.
  }
  else{
    PLOOP(pliter,iv) uu[iv] = uu0[iv] = Uiin[iv];
  }  
  


  ////////////////////////////////
  // START IMPLICIT ITERATIONS
  ////////////////////////////////
  int iter=0;
  int failreturn;
  int f1iter;
  int checkconv,changeotherdt;
  int showmessages=0; // by default 0, don't show any messages for inversion stuff during implicit solver, unless debugging.  Assume any moment of inversion failure is corrected for now unless failure of final inversion done outside implicit solver.
  int showmessagesheavy=0;  // very detailed for common debugging
  int allowlocalfailurefixandnoreport=0; // must be 0 so implicit method knows when really failure

#define MAXF1TRIES 100 // 100 might sound like alot, but Jacobian does 4*4=16 inversions each iteration, and that 100 is only typically needed for very first iteration.
  // goes to f1iter=10 for RADPULSE KAPPAES=1E3 case.  Might want to scale maximum iterations with \tau, although doubling of damping means exponential w.r.t. f1iter, so probably 100 is certainly enough since 2^(-100) is beyond machine precision.

#define DAMPDELTA (0.5) // choose, but best to choose 1/Integer so no machine precision issues.
#define DAMPUNDELTA (1.0/DAMPDELTA)

  // initialize previous 'good inversion' based uu's
  PLOOP(pliter,pl)  uupp[pl]=uuporig[pl]=uup[pl]=uu0orig[pl]=uu[pl];
  
  do{
    iter++;
    
    //vector of conserved at the previous two iterations
    PLOOP(pliter,pl)  uupp[pl]=uup[pl]; // uupp will have solution for inversion: P(uupp)
    PLOOP(pliter,pl)  uup[pl]=uu[pl]; // uup will not necessarily have P(uup) because uu used Newton step.
    PLOOP(pliter,pl)  uuporig[pl]=uu[pl];
    
    //values at zero state
    for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
      int whichcall=1;
      failreturn=f_implicit_lab(whichcall,showmessages, allowlocalfailurefixandnoreport, pin, uu0, uu, fracdtG*realdt, ptrgeom, f1); // modifies uu
      if(failreturn){
        // if initial uu failed, then should take smaller jump from Uiin->uu until settled between fluid and radiation.
        // If here, know original Uiin is good, so take baby steps in case uu needs heavy raditive changes.
        //      return(1);

        // if f1 fails, try going back to Uiin a bit
        fracdtuu0*=DAMPDELTA; // DAMP Uiin->uu0 step that may be too large and generated too large G

        PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0

        if(iter==1){
          // modifies uup and uu as if starting over, and then another call to f_implicit_lab(f1) will change uu by generally smaller amount.
          PLOOP(pliter,pl) uup[pl]=uu[pl]=uu0[pl];
        }
        else{
          // if here, then assume prior uup was good in sense that no failures like P(uup) is good inversion.  And uu=uup before f_implicit_lab(f1) is called.
          // No, not necessarily true, because uu updated with some-sized Newton step without checking if inversion is good for that uu.
          // So need to damp between original uup and uupp .  This essentially damps Newton step now that have knowledge the Newton step was too large as based upon P(U) failure.

          // Avoid G-damping because not needed so far.  If added, competes in non-trivial way with fracuup damping that's required separately for large forces in some problems beyond iter=1 (e.g. NTUBE=31).
          //fracdtG*=DAMPDELTA; // DAMP give only fraction of 4-force to let uu to catch-up

          fracuup*=DAMPDELTA; // DAMP in case Newton step is too large after iter>1 and stuck with certain uu from Newton step that gets stored in uup above.

          PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uuporig[pl];
          //          PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uup[pl];
          PLOOP(pliter,pl) uup[pl]=uu[pl]; // store new version of prior Newton step
        }
        // keep below so can count inversion failures against retry successes in the failure file.
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f1 failed: iter=%d  Backing-up both uu0 and G.: f1iter=%d fracdtuu0=%g fracdtG=%g fracuup=%g\n",iter,f1iter,fracdtuu0,fracdtG,fracuup);
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
      }
      else{
        // then success, so was able to do inversion P(U) with change in radiation on fluid: P(uu[fluid] = uu0[fluid] - (uu[rad]-uu0[rad]))
        // This doesn't necessarily mean could do P(uu0) except for iter=1.
        // Corresponds to success for a certain P(uu0,uu) pair.
        // || ptrgeom->i==23 && ptrgeom->j==24
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"SUCCESS: pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
        break;
      }
    }// end loop over f1iter
    if(f1iter==MAXF1TRIES){
      if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES\n");
      // KORALTODO: If only soft fluid failure (i.e. rho<0 or u<0, then actually allow somehow)  UTOPRIMFAILFIXEDENTROPY or UTOPRIMFAILFIXEDCOLD (but right now those are already always allowed -- but watch for them to see if that's a problem).
      return(1);
    }

    
    // diagnose
    numoff1iter += f1iter;


    if(showmessagesheavy) dualfprintf(fail_file,"i=%d f1: %g %g %g %g\n",ptrgeom->i,f1[0],f1[1],f1[2],f1[3]);


    // see if pre-convergence (happens if force is small or no force at all.  Can't necessarily continue since Jacobian can require arbitrarily large dU on fluid and fail to invert even if no fluid-radiation interaction!
    //test pre-convergence using initial |dU/U|
    // KORALTODO: This isn't a completely general error check since force might be large for fluid that needs itself to have more accuracy, but if using ~NUMEPSILON, won't resolve 4-force of radiation on fluid to better than that.
    FTYPE LOCALPREIMPCONV=(10.0*NUMEPSILON); // more strict than later tolerance
    DLOOPA(ii){
      f3[ii]=f1[ii];
      f3[ii]=fabs(f3[ii]/uu0[URAD0]);
    }
    if(f3[0]<LOCALPREIMPCONV && f3[1]<LOCALPREIMPCONV && f3[2]<LOCALPREIMPCONV && f3[3]<LOCALPREIMPCONV){
      if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter PREDONE1=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3[0],f3[1],f3[2],f3[3]);
      break;
    }




    
    /////////////
    //
    //calculating approximate Jacobian: dUresid(dUrad,G(Urad))/dUrad = dy(x)/dx
    //
    /////////////

    FTYPE localIMPEPS=IMPEPS;
    FTYPE FRACIMPEPSCHANGE=0.1;
    FTYPE del;
    DLOOPA(ii){
	  DLOOPA(jj){

        while(1){

          if(uup[jj+URAD0]==0.) del=localIMPEPS*uup[URAD0];
          else del=localIMPEPS*uup[jj+URAD0];

          // KORALTODO: offset uu (KORALTODO: How to ensure this doesn't have machine precision problems or is good enough difference?)
          // KORALTODO: If ultrarel, then even this small "del" might be too large change in uu and might have bad P(U).  Need to control this "del" better.
          uu[jj+URAD0]=uup[jj+URAD0]-del;
	
          // get dUresid for this offset uu
          int whichcall=2;
          failreturn=f_implicit_lab(whichcall,showmessages,allowlocalfailurefixandnoreport, pin,uu0,uu,fracdtG*realdt,ptrgeom,f2);
          if(failreturn){
            if(showmessages&& debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed: ii=%d jj=%d.  Trying smaller localIMPEPS=%g (giving del=%g) to %g\n",ii,jj,localIMPEPS,del,localIMPEPS*FRACIMPEPSCHANGE);
            localIMPEPS*=FRACIMPEPSCHANGE;
            if(localIMPEPS<NUMEPSILON*10.0){
              if(debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed: ii=%d jj=%d with localIMPEPS=%g (giving del=%g)\n",ii,jj,localIMPEPS,del);
              return(1); // can't go below machine precision for difference else will be 0-0
            }
          }
          else{
            // didn't fail
            break;
          }
        }

	
        if(showmessagesheavy) dualfprintf(fail_file,"i=%d f2: %g %g %g %g\n",ptrgeom->i,f2[0],f2[1],f2[2],f2[3]);

        // get Jacobian (uncentered, ok?  Probably actually best.  Don't want to go back along unknown trajectory in U that might lead to bad P(U))
		J[ii][jj]=(f2[ii] - f1[ii])/(uu[jj+URAD0]-uup[jj+URAD0]);

        // restore uu after getting changed by f_implicit_lab(f2)
		uu[jj+URAD0]=uup[jj+URAD0];
      }
    }
    

    /////////////////////
    //
	//invert Jacobian
    //
    /////////////////////
	failreturn=inverse_44matrix(J,iJ);
    if(failreturn){
      if(debugfail>=2) dualfprintf(fail_file,"inverse_44matrix(J,iJ) failed\n");
      return(1);
    }
   
    /////////
    //
    //updating x, start with previous uu = uup
    //
    /////////
    DLOOPA(ii) x[ii]=uup[ii+URAD0];

    /////////
    //
    // step forward uu=x
    // DAMPFACTOR unused so far because don't know a priori whether to damp.  fracuup does post-inversion effective damping of this Newton step.
    //
    /////////
    DLOOPA(ii){
      DLOOPA(jj){
		x[ii] -= DAMPFACTOR*iJ[ii][jj]*f1[jj];
      }
	}

    /////////
    //
    // assign new uu
    //
    /////////
    DLOOPA(ii) uu[ii+URAD0]=x[ii];


    /////////
    //
    // try to increase DAMPs if was damping or check convergence if no damping active.
    //
    /////////
    checkconv=1;
    changeotherdt=1;
    //    if(fracuup!=1.0){
    if(fabs(fracuup-1.0)>10.0*NUMEPSILON){
      // try increasing amount of uu used
      fracuup*=DAMPUNDELTA;
      checkconv=0;
      changeotherdt=0; // ensure fracuup back to 1.0 first before reverting others.
    }
    //    if(fracdtuu0!=1.0){
    if(fabs(fracdtuu0-1.0)>10.0*NUMEPSILON && changeotherdt){
      // try increasing uu0 away from Uiin to account for full dUother
      fracdtuu0*=DAMPUNDELTA;
      PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0
      checkconv=0;
    }
    //    if(fracdtG!=1.0){
    if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
      // try increasing amount of G applied
      fracdtG*=DAMPUNDELTA;
      checkconv=0;
    }      

    /////////
    //
    // test convergence after Newton step
    // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
    // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
    /////////
    if(checkconv){
      //test convergence using |dU/U|
      DLOOPA(ii){
        f3[ii]=(uu[ii+URAD0]-uup[ii+URAD0]);
        f3[ii]=fabs(f3[ii]/uup[URAD0]);
      }
      
      // see if |dU/U|<tolerance for all components (KORALTODO: What if one component very small and sub-dominant?)
      if(f3[0]<IMPTRYCONV && f3[1]<IMPTRYCONV && f3[2]<IMPTRYCONV && f3[3]<IMPTRYCONV){
        if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iterDONE1=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3[0],f3[1],f3[2],f3[3]);
        break;
      }
    }
 

    /////////
    //
    // see if took too many Newton steps
    //
    /////////
    if(iter>IMPMAXITER){
	  // KORALTODO: Need backup that won't fail.

      //test convergence using |dU/U|
      DLOOPA(ii){
        f3[ii]=(uu[ii+URAD0]-uup[ii+URAD0]);
        f3[ii]=fabs(f3[ii]/uup[URAD0]);
      }

      if(f3[0]<IMPALLOWCONV && f3[1]<IMPALLOWCONV && f3[2]<IMPALLOWCONV && f3[3]<IMPALLOWCONV){
        if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iterDONE1=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3[0],f3[1],f3[2],f3[3]);
        if(debugfail>=2) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab().  But error was ok at %g %g %g %g : checkconv=%d (if checkconv=0, could be issue!)\n",IMPMAXITER,f3[0],f3[1],f3[2],f3[3],checkconv);
        // NOTE: If checkconv=0, then wasn't ready to check convergence and smallness of f3 might only mean smallness of fracuup.  So look for "checkconv=0" cases in fail output.
        break;
      }
      else{
        if(debugfail>=2) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab().  Bad error.\n",IMPMAXITER);
        return(1);
      }
    }

    if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3[0],f3[1],f3[2],f3[3]);
  }// end do
  while(1);
  

  // diagnose
  numofiter+=iter;


  if(showmessagesheavy) dualfprintf(fail_file,"numimplicits=%lld averagef1iter=%g averageiter=%g\n",numimplicits,(FTYPE)numoff1iter/(FTYPE)numimplicits,(FTYPE)numofiter/(FTYPE)numimplicits);


  // get source update as "dU" = dU/dt using real dt that used during implicit iterations, and will eventually use to update U in advance.c.
  DLOOPA(jj) deltas[jj]=(uu[URAD0+jj]-uu0[URAD0+jj])/realdt;

  // apply source update as force
  PLOOP(pliter,pl) radsource[pl] = 0;
  DLOOPA(jj) radsource[UU+jj]    = -SIGNGD3*deltas[jj];
  DLOOPA(jj) radsource[URAD0+jj] = +SIGNGD3*deltas[jj];

  // DEBUG:
  //  DLOOPA(jj) dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d implicitGd[%d]=%g %g\n",nstep,steppart,ptrgeom->i,jj,radsource[UU+jj],radsource[URAD0+jj]);

  // store source update in dUcomp for return.
  sc = RADSOURCE;
  PLOOP(pliter,pl) dUcomp[sc][pl] += radsource[pl];



  return(0);
  
}


// get dt for explicit sub-cyclings
static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother,  FTYPE *CUf, FTYPE *Gd, FTYPE chi, struct of_geom *ptrgeom, FTYPE *dtsub)
{
  int jj;
  int pliter,pl;
  FTYPE idtsub;
  //
  FTYPE Umhd,Urad,Gtot,iUmhd,iUrad;
  //
  FTYPE idtsubs,idtsubt;
  FTYPE idtsubmhd,idtsubrad;
  FTYPE Usmhd,Usrad,Gstot,iUsmhd,iUsrad;
  FTYPE Utmhd,Utrad,Gttot,iUtmhd,iUtrad;
  FTYPE Gddt[NDIM];

  // NOTE: The timestep does not fllow NR1992 S19.2 on diffusion equation step size.  That applies if diffusion was part of flux calculation, not source term.  And it applies to the radiative velocity limiter for the advection's effective wave speed.
  // The relevant NR1992 is S16.6 on stiff source terms.

  // get G*dt that can be compared to Umhd or Urad
  FTYPE realdt=compute_dt(CUf,dt);
  DLOOPA(jj) Gddt[jj] = Gd[jj]*realdt;

  // get updated uu from dUother in case leads to different result
  FTYPE U0[NPR];
  PLOOP(pliter,pl) U0[pl]=UFSET(CUf,dt,Ui[pl],Uf[pl],dUother[pl],0.0);

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d U0=%g realdt=%g dt=%g Ui=%g Uf=%g dUother=%g\n",pl,U0[pl],realdt,dt,Ui[pl],Uf[pl],dUother[pl]);

  // get original U
  FTYPE U[NPR];
  PLOOP(pliter,pl) U[pl]=Ui[pl];


  // get smallest timestep for stiff source terms of 8 equations with a single source term vector.
  // Based upon NR 16.6.6 with removal of factor of two
  if(method==TAUSUPPRESS){
    // KORALTODO: \tau suppression not general enough because G\propto R\tau + \gamma R \tau + \gamma B, and further, this assumes T\sim R.  If T<<R, then suppression by \tau won't be enough to treat stiffness effect on fluid's T.
    // dR^t_t/R^t_t \sim c \gamma_{fluid} \chi dt = \sim c \gamma_{fluid} \tau dt/dx , but might as well just use first version
    
    if(0){
      // use approximate dt along each spatial direction.  chi is based in orthonormal basis
      // get maximum \tau for all relevant directions
      FTYPE dxortho[NDIM],tautotdir[NDIM];
      SLOOPA(jj) dxortho[jj] = (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])));
      // KORALTODO: Should this only depend upon kappa and not kappaes?  Stiffness in source term comes from full chi, so seems to be full chi.
      SLOOPA(jj) tautotdir[jj] = chi * dxortho[jj];
      // only include relevant directions
      FTYPE taumax=SMALL+MAX(MAX(tautotdir[1]*N1NOT1,tautotdir[2]*N2NOT1),tautotdir[3]*N3NOT1);
      
      idtsub=taumax/realdt;
    }
    else{
      FTYPE uu0=q->ucon[TT]; // what enters G for dR^t_t/R^t_t from time part
      //      FTYPE ratchangeRtt=chi * uu0 * realdt * 1.0; // 1.0 = c (1st term)
      FTYPE ratchangeRtt=SMALL+fabs(chi * uu0 * uu0 * realdt * 1.0); // 1.0 = c (2nd term with chi instead of kappaes to be sever and account for B-based term)
      idtsub = ratchangeRtt/realdt; // if ratchange=1, then in principle right at edge of big change.
      // this is like having a "source speed" of v_s\sim \tau \gamma^2 c and limiting the timestep so the source wave only reaches across a cell dxortho in time dt.

      //      dualfprintf(fail_file,"uu0=%g ratchangeRtt=%g idtsub=%g\n",uu0,ratchangeRtt,idtsub);

    }


    //    dualfprintf(fail_file,"i=%d dtsub0=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);

    // account for case where effect on fluid is more than on radiation (where above would only account for effect on radiation)

    // first compare to original U
    jj=TT; Umhd=SMALL+fabs(U[UU+jj]);
    jj=TT; Urad=fabs(U[URAD0+jj]);
    idtsub=MAX(idtsub,idtsub*Urad/Umhd);

    //    dualfprintf(fail_file,"i=%d dtsub1=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);
       
    // also compare against changed U=U0
    jj=TT; Umhd=SMALL+fabs(U0[UU+jj]);
    jj=TT; Urad=fabs(U0[URAD0+jj]);
    idtsub=MAX(idtsub,idtsub*Urad/Umhd);

    //    dualfprintf(fail_file,"i=%d dtsub2=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);
       
  }
  // below is if Diffusion is part of source term as in Koral
  // source term should lead to small (<1/2) change in conserved quantities
  else if(method==SPACETIMESUBSPLITNONE){
	// merged space-time to avoid negligible total momentum with large update needing to be resolved.
    // !! Doesn't make sense if Umhd>>Urad, then in optically thick case the radiation needs itself to be stiffly evolved despite much larger Umhd !!
	Umhd=Urad=Gtot=0.0;
	DLOOPA(jj) Umhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	DLOOPA(jj) Urad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	DLOOPA(jj) Gtot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	iUmhd=1.0/(fabs(Umhd)+SMALL);
	iUrad=1.0/(fabs(Urad)+SMALL);
	idtsub=SMALL+fabs(Gtot*MIN(iUmhd,iUrad));

    //    dualfprintf(fail_file,"UMHD: %g %g %g %g %g %g\n",Umhd,Urad,Gtot,iUmhd,iUrad);

  }
  else if(method==SPACETIMESUBSPLITTIME){
	// won't be efficient if v~0
	// if v<<1 and G is mid-range but still negligible, then dt will be incredibly small and code will halt.
	Usmhd=Usrad=Gstot=0.0;
	Utmhd=Utrad=Gttot=0.0;
	SLOOPA(jj) Usmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Utmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	SLOOPA(jj) Usrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Utrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	SLOOPA(jj) Gstot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Gttot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	iUsmhd=1.0/(fabs(Usmhd)+SMALL);
	iUtmhd=1.0/(fabs(Utmhd)+SMALL);
	iUsrad=1.0/(fabs(Usrad)+SMALL);
	iUtrad=1.0/(fabs(Utrad)+SMALL);
	idtsubs=SMALL+fabs(Gstot*MIN(iUsmhd,iUsrad));
	idtsubt=SMALL+fabs(Gttot*MIN(iUtmhd,iUtrad));
	idtsub=MAX(idtsubs,idtsubt);
  }
  else if(method==SPACETIMESUBSPLITALL){
	// won't be efficient if flow becomes grid-aligned or if v~0
	Usmhd=Usrad=Gstot=0.0;
	idtsub=0.0;
	DLOOPA(jj){
	  Umhd = fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  Urad = fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  Gtot = fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  iUmhd=1.0/(fabs(Umhd)+SMALL);
	  iUrad=1.0/(fabs(Urad)+SMALL);
	  idtsub=MAX(idtsub,SMALL+fabs(Gtot*MIN(iUmhd,iUrad)));
	}
  }
  else if(method==SPACETIMESUBSPLITSUPERALL){
	// won't be efficient if flow becomes grid-aligned or if v~0 or if radiation neglibile contribution to fluid dynamics
    // but for stable evolution of the radiation independent of the fluid, need separate mhd and rad conditionals.
	Usmhd=Usrad=Gstot=0.0;
	idtsub=0.0;
	DLOOPA(jj){
	  Umhd = fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  Urad = fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  Gtot = fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	  iUmhd=1.0/(fabs(Umhd)+SMALL);
	  iUrad=1.0/(fabs(Urad)+SMALL);
	  idtsub=MAX(idtsub,SMALL+fabs(Gtot*iUmhd));
	  idtsub=MAX(idtsub,SMALL+fabs(Gtot*iUrad));
	}
  }
  else if(method==SPACETIMESUBSPLITMHDRAD){
	// merged space-time to avoid negligible total momentum with large update needing to be resolved.
	Umhd=Urad=Gtot=0.0;
    idtsub=0.0;
	DLOOPA(jj) Umhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	DLOOPA(jj) Urad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	DLOOPA(jj) Gtot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	iUmhd=1.0/(fabs(Umhd)+SMALL);
	iUrad=1.0/(fabs(Urad)+SMALL);
	idtsub=MAX(idtsub,SMALL+fabs(Gtot*iUmhd));
	idtsub=MAX(idtsub,SMALL+fabs(Gtot*iUrad));

    //    dualfprintf(fail_file,"UMHD: %g %g %g %g %g %g\n",Umhd,Urad,Gtot,iUmhd,iUrad);

  }
  else if(method==SPACETIMESUBSPLITTIMEMHDRAD){
	// won't be efficient if v~0
	// if v<<1 and G is mid-range but still negligible, then dt will be incredibly small and code will halt.
	Usmhd=Usrad=Gstot=0.0;
	Utmhd=Utrad=Gttot=0.0;
	SLOOPA(jj) Usmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Utmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	SLOOPA(jj) Usrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Utrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
	SLOOPA(jj) Gstot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	jj=TT;     Gttot += fabs(Gddt[jj]*Gddt[jj]*ptrgeom->gcon[GIND(jj,jj)]);
	iUsmhd=1.0/(fabs(Usmhd)+SMALL);
	iUtmhd=1.0/(fabs(Utmhd)+SMALL);
	iUsrad=1.0/(fabs(Usrad)+SMALL);
	iUtrad=1.0/(fabs(Utrad)+SMALL);
    idtsub=SMALL;
	idtsub=MAX(idtsub,SMALL+fabs(Gstot*iUsmhd));
	idtsub=MAX(idtsub,SMALL+fabs(Gstot*iUsrad));
	idtsub=MAX(idtsub,SMALL+fabs(Gttot*iUtmhd));
	idtsub=MAX(idtsub,SMALL+fabs(Gttot*iUtrad));
  }


  
  // what to return
  *dtsub=COURRADEXPLICIT/idtsub;


  //  dualfprintf(fail_file,"*dtsub=%g idtsub=%g method=%d\n",*dtsub,idtsub,method);
 
}


// prepare for explicit sub-cycling
static int koral_explicit_source_rad_prepare(FTYPE *pr, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother, FTYPE *prnew, FTYPE *Unew, struct of_state *qnew)
{
  int pliter,pl;
  struct of_newtonstats newtonstats;
  int showmessages=0;  // don't show during prepare unless debugging
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int allowlocalfailurefixandnoreport=0; // need to know if failed

  int nochange=1;
  if(USEDUINRADUPDATE){
    // backup pr and Uiin
    PLOOP(pliter,pl) prnew[pl]=pr[pl];
    PLOOP(pliter,pl) Unew[pl]=UFSET(CUf,dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);

    // Get P(U0) in case dUother is non-zero
    if(Utoprimgen_failwrapper(showmessages,allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, Unew, ptrgeom, prnew, &newtonstats)){
      // failed to get updated primitive
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed in koral_explcit_source_rad(), so will use primitive without flux+geom update.\n");
      // KORALTODO: If only soft fluid failure (i.e. rho<0 or u<0, then actually allow somehow)
    }
    else{
      nochange=0;
      // then got good new primitive that is stored in prnew already
      // also need to update qnew (doesn't use q)
      get_state_uconucovonly(prnew, ptrgeom, qnew);
      get_state_uradconuradcovonly(prnew, ptrgeom, qnew);
    }
  }

  if(nochange){
    // revert P and U and copy over q
    PLOOP(pliter,pl) prnew[pl]=pr[pl];
    PLOOP(pliter,pl) Unew[pl]=Uiin[pl];
    *qnew=*q;
  }

  //  dualfprintf(fail_file,"prepare: nochange=%d\n",nochange);
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Uiin=%g Unew=%g\n",pl,Uiin[pl],Unew[pl]);


  return(nochange);

}



// compute changes to U (both T and R) using implicit method
// NOTEMARK: The explicit scheme is only stable if the fluid speed is order the speed of light.  Or, one can force the explicit scheme to use vrad=c.
// Only change dUcomp, and can overwrite prnew, Unew, and qnew since "prepare" function isolated original values already
static void koral_explicit_source_rad(FTYPE *prnew, FTYPE *Unew, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *qnew, FTYPE (*dUcomp)[NPR])
{
  FTYPE Gd[NDIM], radsource[NPR];
  int pliter, pl, jj, sc;
  FTYPE chi;
  FTYPE pr0[NPR],U0[NPR];
  int numcycles=0;

  struct of_newtonstats newtonstats;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.

  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;


  // backup pr and Uiin (Uiin should already be changed with dUother and new q should be computed in q)
  PLOOP(pliter,pl) pr0[pl]=prnew[pl];
  PLOOP(pliter,pl) U0[pl]=Unew[pl];



  // initialize source update
  sc = RADSOURCE;
  PLOOP(pliter,pl) radsource[pl] = 0;


  // Based upon size of Gd, sub-cycle this force.
  // 1) calc_Gd()
  // 2) locally set dtsub~dt/\tau or whatever it should be
  // 3) update T^t_\nu and R^t_\nu
  // 4) U->P locally
  // 5) repeat.


  FTYPE dttrue=0.0,dtcum=0.0;  // cumulative sub-cycle time
  FTYPE dtdiff;
  FTYPE dtsub,dtsubold,dtsubuse;
  FTYPE realdt=compute_dt(CUf,dt);
  int method;
  int showmessages=0;
  int allowlocalfailurefixandnoreport=0; // need to see if any failures.

  int itersub=0;
  while(1){

	// get 4-force
	calc_Gd(prnew, ptrgeom, qnew, Gd, &chi);


	// see if need to sub-cycle
	// dynamically change dt to allow chi to change during sub-cycling.
    // use subcycling if directly requested or if reversion from implicit and that reversion was requested
	if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODEXPLICITSUBCYCLE || (WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICIT || WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICITEXPLICITCHECK) && IMPLICITREVERTEXPLICIT){

      // SPACETIMESUBSPLITMHDRAD doesn't work -- generates tons of noise in prad1 with COURRADEXPLICIT=0.2, and was asymmetric in x.
      method=TAUSUPPRESS; // forced -- only method that is efficient and effective and noise free at moderate optical depths.

	  // get dt for explicit stiff source term sub-cycling
      FTYPE Ualt[NPR];
      FTYPE dUotheralt[NPR]={0.0};
      PLOOP(pliter,pl) Ualt[pl]=Unew[pl]; // at most doubles value
      // assume Unew already accounts for Uf and dUother
	  get_dtsub(method,prnew,qnew,Unew, Ualt,dUotheralt, CUf, Gd, chi, ptrgeom, &dtsub);
      if(realdt/dtsub>MAXSUBCYCLES){
        if(debugfail>=2) dualfprintf(fail_file,"dtsub very small: %g with realdt=%g and only allowing MAXSUBCYCLES=%d subcycles, so limit dtsub\n",dtsub,realdt,MAXSUBCYCLES);
        dtsub=realdt/(FTYPE)MAXSUBCYCLES;
      }
      if(itersub==0) dtsubold=dtsubuse=dtsub;
      else{
        // ensure don't change step too fast
        if(dtsub>dtsubold*(1.0+MAXEXPLICITSUBSTEPCHANGE)) dtsubuse=dtsubold*(1.0+MAXEXPLICITSUBSTEPCHANGE);
        else dtsubuse=dtsub;
        // override if dtsub is larger than realdt, indicating really done with iterations and reached some equilibrium, so no longer necessary to check vs. dtsubold
        if(dtsub>realdt) dtsubuse=dtsub;
        // need to compare with previous actual dt
        dtsubold=dtsubuse;
      }

	  //	  dualfprintf(fail_file,"i=%d chi=%g realdt=%g fakedt=%g dtsub=%g dtcum=%g dttrue=%g itersub=%d\n",ptrgeom->i,chi,realdt,fakedt,dtsub,dtcum,dttrue,itersub);

	  // time left to go in sub-cycling
	  dtdiff=MAX(realdt-dtcum,0.0);
      dttrue=MIN(dtsubuse,dtdiff);

	  if(realdt>dtsubuse && itersub==0 || dtdiff>0.0 && itersub>0){ // then sub-cycling

        
        if(showmessages&&debugfail>=2) dualfprintf(fail_file,"DoingSUBCYCLE: iter=%d : dtsub=%g dtsubuse=%g dtdiff=%g dttrue=%g dtcum=%g realdt=%g\n",itersub,dtsub,dtsubuse,dtdiff,dttrue,dtcum,realdt);


		// initialize counters
		newtonstats.nstroke=newtonstats.lntries=0;

		// get change in conserved quantity between fluid and radiation
		// take step, but only go up to exactly realdt in total time
		DLOOPA(jj) radsource[UU+jj]    += -SIGNGD*Gd[jj]*dttrue/realdt;
		DLOOPA(jj) radsource[URAD0+jj] += +SIGNGD*Gd[jj]*dttrue/realdt;
		dtcum+=dttrue; // cumulative dttrue
	  
		// get Unew[U0,dUrad]
		PLOOP(pliter,pl) Unew[pl]       = U0[pl];
		DLOOPA(jj)       Unew[UU+jj]    = U0[UU+jj]    + radsource[UU+jj]*realdt;
		DLOOPA(jj)       Unew[URAD0+jj] = U0[URAD0+jj] + radsource[URAD0+jj]*realdt;
	  
		// Get U->P
		// OPTMARK: Should optimize this to  not try to get down to machine precision
		// KORALTODO: NOTEMARK: If failure, then need to really fix-up or abort this implicit solver!
		// Get P(U)
		if(Utoprimgen_failwrapper(showmessages, allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, Unew, ptrgeom, prnew, &newtonstats)){
		  dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of koral_explicit_source_rad()\n");
		  //myexit(347366436); // KORALTODO: if explicit dies, nothing to revert to, except to assume force not important for overall dynamics.
          // KORALTODO: If only soft fluid failure (i.e. rho<0 or u<0, then actually allow somehow)
		  return;
		}

		// re-get needed q's
		get_state_uconucovonly(prnew, ptrgeom, qnew);
		get_state_uradconuradcovonly(prnew, ptrgeom, qnew);

		itersub++;
	  }
	  else if(itersub==0){
		// equal and opposite forces
		DLOOPA(jj) radsource[UU+jj]    = -SIGNGD*Gd[jj];
		DLOOPA(jj) radsource[URAD0+jj] = +SIGNGD*Gd[jj];
		break;
	  }
	  else break;
	}
	else{
	  // not allowing sub-cycling
	  // equal and opposite forces
	  DLOOPA(jj) radsource[UU+jj]    = -SIGNGD*Gd[jj];
	  DLOOPA(jj) radsource[URAD0+jj] = +SIGNGD*Gd[jj];
	  break;
	}

  }// done looping



#if(0)
#define RHO_AMB (MPERSUN*MSUN/(LBAR*LBAR*LBAR)) // in grams per cm^3 to match koral's units with rho=1
  // DEBUG
  FTYPE dxdxp[NDIM][NDIM];
  dxdxprim_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,dxdxp);
  DLOOPA(jj) dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d explicitGd[%d]=%g vs=%g\n",nstep,steppart,ptrgeom->i,jj,Gd[jj]/RHO_AMB*sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])),dx[1]/(dt/cour)*dxdxp[1][1]);
#endif


  // apply 4-force as update in dUcomp[][]
  // only changed this quantity, none of other among function arguments
  PLOOP(pliter,pl) dUcomp[sc][pl] += radsource[pl];

  
}

void koral_source_rad(FTYPE *pin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *dUother, FTYPE (*dUcomp)[NPR])
{
  int pliter,pl;
  int showmessages=0; // 0 ok if not debugging and think everything works.
  int showmessagesheavy=0;


  // make code use "orig" values that below can modify, so return preserves no changes to these no matter what internal functions do.
#if(1)
  // save pin, Uiin, Ufin, and q to avoid being modified upon return
  FTYPE pinorig[NPR],Uiinorig[NPR],Ufinorig[NPR];
  struct of_state qorigmem;
  struct of_state *qorig=&qorigmem;
  PLOOP(pliter,pl){
    pinorig[pl]=pin[pl];
    Uiinorig[pl]=Uiin[pl];
    Ufinorig[pl]=Ufin[pl];
  }
  *qorig=*q;
#else
  FTYPE *pinorig,*Uiinorig,*Ufinorig;
  struct of_state *qorig;

  pinorig=pin;
  Uiinorig=Uiin;
  Ufinorig=Ufin;
  qorig=q;
#endif


  //  if(ptrgeom->i==10 && ptrgeom->k==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"BEFORE KORALSOURCE: pl=%d U=%g\n",pl,Uiinorig[pl]);
  //  }

  // KORALTODO: SUPERGODMARK: Need to add NR 2007 page940 17.5.2 StepperSie method here as higher-order alternative if 1st order Newton breaks

  if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODEXPLICIT || WHICHRADSOURCEMETHOD==RADSOURCEMETHODEXPLICITSUBCYCLE){

    FTYPE prnew[NPR],Unew[NPR];
    struct of_state qnew;

    // prepare for explicit sub-cycling, possibly trying to use updated U[dUgeom+dUriemann] to start with so can have balanced forces
    int failprepare=koral_explicit_source_rad_prepare(pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother , prnew, Unew, &qnew);
    if(failprepare){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"explicit failed to prepare, which means subcycle won't be correctly determining Gd or dtsub: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
      // still do explicit anyways, since best can do with the choice of method -- will fail possibly to work if stiff regime, but ok in non-stiff.
    }
	koral_explicit_source_rad( prnew, Unew, CUf, ptrgeom, &qnew, dUcomp);
  }
  else if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICIT){
	
	int failimplicit=koral_implicit_source_rad( pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);
	
	if(failimplicit){
	  if(IMPLICITREVERTEXPLICIT){
        FTYPE prnew[NPR],Unew[NPR];
        struct of_state qnew;
        int failprepare=koral_explicit_source_rad_prepare(pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother , prnew, Unew, &qnew);
        if(failprepare){
          if(showmessages && debugfail>=2) dualfprintf(fail_file,"implicit failed, and then explicit failed to prepare, which means subcycle won't be correctly determining Gd or dtsub: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
          // still do explicit anyways, since best can do with the choice of method -- will fail possibly to work if stiff regime, but ok in non-stiff.
        }
        koral_explicit_source_rad(prnew, Unew, CUf, ptrgeom, &qnew, dUcomp);
      }
	  else{
		// then can only fail
		myexit(39475251);
	  }
	}

  }
  else if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICITEXPLICITCHECK){



    // Check if can completely skip optical depth effects
    FTYPE Gd[NDIM],chi,dtsub;
    calc_Gd(pinorig, ptrgeom, qorig, Gd, &chi);
    int method=TAUSUPPRESS; //WHICHSPACETIMESUBSPLIT; // chosen in global.nondepmnemonics.rad.h or by user in init.h (assume quite conservative since can just do implicit step)
    get_dtsub(method,pinorig,qorig,Uiinorig, Ufinorig, dUother, CUf, Gd, chi, ptrgeom, &dtsub);
    FTYPE realdt = compute_dt(CUf,dt);
    if(NUMEPSILON*dtsub>=realdt){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"Nothing in source step to do because effect on radiation or fluid is below machine precision\n");
      return;
    }


#if(0)
    // NO: Need to compare R\tau vs. E in sense of G/U.
    // first check if \tau is exceedingly small, so fluid-rad interaction can be skipped.  Avoids inversion for explicit prepare.
    FTYPE tautot[NDIM],tautotmax;
    calc_tautot(pinorig, ptrgeom, tautot, &tautotmax);
    if(tautotmax<MINTAUSOURCE){
      if(debugfail>=2) dualfprintf(fail_file,"maximum optical depth for ijk=%d %d %d was tautotmax=%g , which is negligible enough to completely avoid source term explicit or implicit.\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,tautotmax);
      //return; // then nothing to do.
    }
    else{
      if(debugfail>=2) dualfprintf(fail_file,"maximum optical depth for ijk=%d %d %d was tautotmax=%g.\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,tautotmax);
    }
#endif




    FTYPE prnew[NPR],Unew[NPR];
    struct of_state qnew;
    int failprepare=0;
    int doimplicit=1;


    // first check if can just step with explicit scheme without sub-cycling
    if(dtsub>=realdt){
      // then just take explicit step!
      failprepare=koral_explicit_source_rad_prepare(pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother , prnew, Unew, &qnew);
      if(failprepare){
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"was testing if can do explicit instead of implicit, but explicit failed to prepare, which means subcycle won't be correctly determining Gd or dtsub: ijk=%d %d %d.\n So in this case, just avoid explicit and do implicit.\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
        doimplicit=1;
      }
      else{
        doimplicit=0;
        // assumes below doesn't modify pinorig,Uiinorig, Ufinorig, CUf,ptrgeom, or qorig
        koral_explicit_source_rad(prnew, Unew, CUf, ptrgeom, &qnew, dUcomp);
        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Was able to take explicit step: ijk=%d %d %d : realdt=%g dtsub=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,realdt,dtsub);
      }
	}
    else{
      doimplicit=1;
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Prepared for explicit step, but would require subcycling, so do implicit directly: ijk=%d %d %d : realdt=%g dtsub=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,realdt,dtsub);
    }



    if(doimplicit){
	  // else just continue with implicit, with only cost extra calc_Gd() evaluation
	
	  // DEBUG (to comment out!):
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Had to take implicit step\n");

      //      PLOOP(pliter,pl) dualfprintf(fail_file,"BEFORE: i=%d j=%d pl=%d pin=%g pinorig=%g Uiin=%g Uiinorig=%g Ufin=%g Ufinorig=%g\n",ptrgeom->i,ptrgeom->j,pl,pin[pl],pinorig[pl],Uiin[pl],Uiinorig[pl],Ufin[pl],Ufinorig[pl]);
	
	  int failimplicit=koral_implicit_source_rad( pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);
	  if(showmessages && debugfail>=2 && failimplicit) dualfprintf(fail_file,"Implicit failed in RADSOURCEMETHODIMPLICITEXPLICITCHECK.\n");

      //      PLOOP(pliter,pl) dualfprintf(fail_file,"AFTER: i=%d j=%d pl=%d pin=%g pinorig=%g Uiin=%g Uiinorig=%g Ufin=%g Ufinorig=%g\n",ptrgeom->i,ptrgeom->j,pl,pin[pl],pinorig[pl],Uiin[pl],Uiinorig[pl],Ufin[pl],Ufinorig[pl]);


	  if(failimplicit){
		if(IMPLICITREVERTEXPLICIT){
          failprepare=koral_explicit_source_rad_prepare(pinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother , prnew, Unew, &qnew);
          if(failprepare){
            if(debugfail>=2) dualfprintf(fail_file,"implicit failed, and failed to prepare for explicit, which means subcycle won't be correctly determining Gd or dtsub: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
          }
		  koral_explicit_source_rad(prnew, Unew, CUf, ptrgeom, &qnew, dUcomp);
		}
		else{
		  // then can only fail
		  myexit(39475252);
		}
	  }
	}

  }
  else if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODNONE){

  }
  else{

	dualfprintf(fail_file,"3 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
	myexit(18754363);

  }

 
  // if(ptrgeom->i==10 && ptrgeom->k==0){
  //    dualfprintf(fail_file,"dUcomp=%g %g %g %g\n",dUcomp[RADSOURCE][UU],dUcomp[RADSOURCE][U1],dUcomp[RADSOURCE][U2],dUcomp[RADSOURCE][U3]);
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"AFTER KORALSOURCE: pl=%d U=%g\n",pl,Uiinorig[pl]);
  //  }

}



//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************
//absorption
void calc_kappa(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa)
{

  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappa()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  *kappa = calc_kappa_user(rho,T,xx,yy,zz);
  //  dualfprintf(fail_file,"kappaabs=%g\n",*kappa);
}

//scattering
void calc_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa)
{  
  extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappaes()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  *kappa = calc_kappaes_user(rho,T,xx,yy,zz);
  //  dualfprintf(fail_file,"kappaes=%g\n",*kappa);
}

void calc_chi(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *chi)
{
  FTYPE kappa,kappaes;
  calc_kappa(pr,ptrgeom,&kappa);
  calc_kappaes(pr,ptrgeom,&kappaes);
  
  *chi=kappa+kappaes;
}


static void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G, FTYPE* chireturn)
{
  calc_Gu(pp, ptrgeom, q, G, chireturn);
  indices_21(G, G, ptrgeom);
}


//**********************************************************************
//****** takes radiative stress tensor and gas primitives **************
//****** and calculates contravariant four-force ***********************
//**********************************************************************
static void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu, FTYPE* chireturn) 
{
  int i,j,k;
  
  //radiative stress tensor in the lab frame
  FTYPE Rij[NDIM][NDIM];

  //this call returns R^i_j, i.e., the first index is contra-variant and the last index is co-variant
  mhdfull_calc_rad(pp, ptrgeom, q, Rij);

  //  FTYPE Rijkoral[NDIM][NDIM];
  //  DLOOP(i,j) Rijkoral[i][j] =0;
  //  DLOOP(i,j) DLOOPA(k) Rijkoral[i][j] += Rij[i][k]*ptrgeom->gcon[GIND(k,j)];
  //  DLOOP(i,j) dualfprintf(fail_file,"i=%d j=%d Rij=%g Rijkoral=%g\n",i,j,Rij[i][j],Rijkoral[i][j]);

  //the four-velocity of fluid in lab frame
  FTYPE *ucon,*ucov;

  ucon = q->ucon;
  ucov = q->ucov;
  
  
  //gas properties
#if(0)
  FTYPE rho=pp[RHO];
  FTYPE u=pp[UU];
  FTYPE p= pressure_rho0_u_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho,u);
  FTYPE T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  FTYPE B = SIGMA_RAD*pow(T,4.)/Pi;
  FTYPE Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
#else
  FTYPE rho=pp[RHO];
  FTYPE u=pp[UU];
  FTYPE T=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho,u);
  FTYPE B=0.25*ARAD_CODE*pow(T,4.)/Pi; // KORALTODO: Why Pi here?
  // KORALTODO: Also, doesn't this only allow for thermal black body emission?
  //  dualfprintf(fail_file,"rho=%g u=%g T=%g B=%g Erad=%g\n",rho,u/rho,T*TEMPBAR,B/rho,pp[PRAD0]/rho);
#endif
  FTYPE kappa,kappaes;
  calc_kappa(pp,ptrgeom,&kappa);
  calc_kappaes(pp,ptrgeom,&kappaes);
  FTYPE chi=kappa+kappaes;
  
  //contravariant four-force in the lab frame
  
  //R^a_b u_a u^b
  FTYPE Ruu=0.; DLOOP(i,j) Ruu+=Rij[i][j]*ucov[i]*ucon[j];

#if(0) // DEBUG
  DLOOPA(j) dualfprintf(fail_file,"ucov[%d]=%g ucon[%d]=%g\n",j,ucov[j],j,ucon[j]);
  FTYPE *uradcon,*uradcov;
  uradcon=q->uradcon;
  uradcov=q->uradcov;
  DLOOPA(j) dualfprintf(fail_file,"uradcov[%d]=%g uradcon[%d]=%g\n",j,uradcov[j]*sqrt(fabs(ptrgeom->gcon[GIND(j,j)])),j,uradcon[j]*sqrt(fabs(ptrgeom->gcov[GIND(j,j)])));
  FTYPE Rijuu[NDIM][NDIM];
  DLOOP(i,j) Rijuu[i][j]=0.0;
  int b;
  DLOOP(i,j) DLOOPA(b) Rijuu[i][j]+=Rij[i][b]*ptrgeom->gcon[GIND(b,j)];

  DLOOP(i,j) dualfprintf(fail_file,"i=%d j=%d Rij=%g Rijuu=%g\n",i,j,Rij[i][j]/rho,Rijuu[i][j]/rho);
#endif
  
  FTYPE Ru;
  DLOOPA(i){

    Ru=0.; DLOOPA(j) Ru+=Rij[i][j]*ucon[j];

    Gu[i]=-chi*Ru - (kappaes*Ruu + kappa*4.*Pi*B)*ucon[i]; // KORALTODO: Is 4*Pi here ok?


    //    dualfprintf(fail_file,"i=%d Ruu=%g Ru=%g kappa=%g kappaes=%g chi=%g B=%g ucon=%g %g %g %g\n",i,Ruu,Ru,kappa,kappaes,chi,B,ucon[0],ucon[1],ucon[2],ucon[3]);
#if(0)
	dualfprintf(fail_file,"i=%d : Ruu=%g Ru=%g chi=%g Gu=%g\n",i,Ruu/rho,Ru/RHO_AMB*sqrt(fabs(ptrgeom->gcov[GIND(i,i)])),chi,Gu[i]/RHO_AMB*sqrt(fabs(ptrgeom->gcov[GIND(i,i)])));
#endif

  }


  *chireturn=chi; // if needed

}

int vchar_all(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxall, FTYPE *vminall,int *ignorecourant)
{
  FTYPE vminmhd,vmaxmhd;
  FTYPE vminrad,vmaxrad;
  FTYPE vminrad2,vmaxrad2;
  
  vchar_each(pr, q, dir, geom, &vmaxmhd, &vminmhd, &vmaxrad, &vminrad, &vmaxrad2, &vminrad2, ignorecourant);
  // below correct even if EOMRADTYPE==EOMRADNONE because vchar_each() sets both mhd and rad to be mhd and so below always chooses the mhd values.
  *vminall=MIN(vminmhd,vminrad2); // vminrad2 for dt
  *vmaxall=MAX(vmaxmhd,vmaxrad2); // vmaxrad2 for dt
  

  return(0);
}

int vchar_each(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxmhd, FTYPE *vminmhd, FTYPE *vmaxrad, FTYPE *vminrad, FTYPE *vmaxrad2, FTYPE *vminrad2,int *ignorecourant)
{
  
  vchar(pr, q, dir, geom, vmaxmhd, vminmhd,ignorecourant);
  if(EOMRADTYPE!=EOMRADNONE){
    vchar_rad(pr, q, dir, geom, vmaxrad, vminrad, vmaxrad2, vminrad2, ignorecourant);
  }
  else{// default as if no other values for wave speeds
    *vmaxrad2=*vmaxrad=*vmaxmhd;
    *vminrad2=*vminrad=*vminmhd;
  }

  return(0);
}

int vchar_rad(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin, FTYPE *vmax2, FTYPE *vmin2,int *ignorecourant)
{

  
  // compute chi
  // Assume computed as d\tau/dorthonormallength as defined by user.
  // Assume \chi defined in fluid frame (i.e. not radiation frame).
  
  // KORALTODO: below 2 lines were how harm = koralinsert was before
  // KORALTODO: Need to still use total wave speed for timestep limit!
  FTYPE kappa,chi;
  calc_chi(pr,geom,&chi);
  //  calc_kappa(pr,geom,&kappa);
  //  chi=kappa;


  //characterisitic wavespeed in the radiation rest frame
  FTYPE vrad2=THIRD;
  FTYPE vrad2limited;

  if(chi>0.0){// && WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICIT){
    // NOT DOING THIS:
    // compute tautot assuming chi is optical depth per unit grid dx[1-3].  I.e. calc_chi() computes grid-based opacity
    // tautot is the total optical depth of the cell in dim dimension
    //  tautot = chi * dx[dir];

    // DOING THIS:
    // KORALTODO: Approximation to any true path, but approximation is sufficient for approximate wave speeds.
    // \tau_{\rm tot}^2 \approx \chi^2 [dx^{dir} \sqrt{g_{dirdir}}]^2 
    FTYPE tautotsq,vrad2tau;
    // Note that tautot is frame independent once multiple \chi by the cell length.  I.e. it's a Lorentz invariant.
    tautotsq = chi*chi * dx[dir]*dx[dir]*(geom->gcov[GIND(dir,dir)]);

    //    dualfprintf(fail_file,"chi=%g dx=%g dir=%d tautot=%g\n",chi,dx[dir],dir,sqrt(tautotsq));
  
    vrad2tau=(4.0/3.0)*(4.0/3.0)/tautotsq; // KORALTODO: Why 4.0/3.0 ?  Seems like it should be 2.0/3.0 according to NR1992 S19.2.6 or NR2007 S20.2 with D=1/(3\chi), but twice higher speed is only more robust.
    vrad2limited=MIN(vrad2,vrad2tau);

	// NOTEMARK: For explicit method, this will lead to very large dt relative to step desired by explicit method, leading to ever more sub-cycles for WHICHRADSOURCEMETHOD==RADSOURCEMETHODEXPLICITSUBCYCLE method.

    // TODOMARK: I wonder if another possibility is to use a speed limiter in the advection equation.  With my pseudo-Newtonian code is has a limiter on the sound and Alfven speeds following the idea of limiting the Alfven speed by Miller & Stone (2000, http://adsabs.harvard.edu/abs/2000ApJ...534..398M).  That is, there must be a way to insert a term into the radiation advection equations to limit the velocity to ~c/\tau that only becomes effective at and beyond that speed.  Then the Jacobian would be modified (Or thinking of how the Jacobian could get modified, one gets a different equation of motion).

  }
  else{
    vrad2limited=vrad2;
  }

  // for setting flux so diffusive term is not exaggerated in high \tau regions
  simplefast_rad(dir,geom,q,vrad2limited,vmin,vmax);
#if(1)
  *vmin2=*vmin;
  *vmax2=*vmax;
#elif(0)
  // for setting timestep since advective part has no knowledge of \tau-limited velocity
  simplefast_rad(dir,geom,q,vrad2,vmin2,vmax2);
#else
  // works even if not using damping of implicit solver
  simplefast_rad(dir,geom,q,2.0/3.0,vmin2,vmax2);
#endif
  // KORALTODO: It seems to not always be required (shock tests don't need it even for non-rel optically thick problem), so maybe there's a way to have vrad2 -> c to keep the optically thick regime stable by interpolating vrad2 to become -> c when \tau>>1.

  
  return(0);
}

// get lab-frame 3-velocity for radiative emission in radiative frame
static int simplefast_rad(int dir, struct of_geom *geom,struct of_state *q, FTYPE vrad2,FTYPE *vmin, FTYPE *vmax)
{
  extern int simplefast(int dir, struct of_geom *geom,struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax);

  //need to substitute ucon,ucov with uradcon,uradcov to fool simplefast
  FTYPE ucon[NDIM],ucov[NDIM];
  int ii;
  DLOOPA(ii){
    ucon[ii]=q->ucon[ii];
    ucov[ii]=q->ucov[ii];
    q->ucon[ii]=q->uradcon[ii];
    q->ucov[ii]=q->uradcov[ii];
  }

  //calculating vmin, vmax
  simplefast(dir,geom,q,vrad2,vmin,vmax);

  //restoring gas 4-velocities
  DLOOPA(ii){
    q->ucon[ii]=ucon[ii];
    q->ucov[ii]=ucov[ii];
  }


#if(0)
  // Cartesian-Minkowski speed-of-light limit of radiation velocity
 FTYPE dxdxp[NDIM][NDIM];
  dxdxprim_ijk(geom->i, geom->j, geom->k, geom->p, dxdxp);
  // characeristic wavespeeds are 3-velocity in lab-frame
  *vmin=-1.0/dxdxp[dir][dir]; // e.g. dxdxp = dr/dx1
  *vmax=+1.0/dxdxp[dir][dir];
#endif


  return(0);
}

// convert primitives to conserved radiation quantities
int p2u_rad(FTYPE *pr, FTYPE *Urad, struct of_geom *ptrgeom, struct of_state *q)
{
  FTYPE Uraddiag;
  
  mhd_calc_rad(pr, TT, ptrgeom, q, Urad);

  // SUPERGODMARK: ensure \detg applied when filling conserved quantity

  return(0);
}

// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  void compute_1plusud0(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])

  // urad^\mu
  // ucon_calc() assumes primitive velocities are in U1 through U3, but otherwise calculation is identical for radiation velocity, so just shift effective list of primitives so ucon_calc() operates on U1RAD through U3RAD
  MYFUN(ucon_calc(&pr[URAD1-U1], ptrgeom, q->uradcon,q->othersrad) ,"phys.c:get_state()", "ucon_calc()", 1);
  // urad_\mu
  lower_vec(q->uradcon, ptrgeom, q->uradcov);


  return (0);
}


void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM])
{
  int jj,kk;
  
  if(EOMRADTYPE!=EOMRADNONE){
	DLOOPA(jj){
	  mhd_calc_rad( pr, jj, ptrgeom, q, &(radstressdir[jj][0]) );
	}
  }
  else DLOOP(jj,kk) radstressdir[jj][kk]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.
}

// compute radiation stres-energy tensor
void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *ptrgeom, struct of_state *q, FTYPE *radstressdir)
{
  int jj;

  // R^{dir}_{jj} radiation stress-energy tensor
  if(EOMRADTYPE==EOMRADEDD){
	// force radiation frame to be fluid frame
    DLOOPA(jj) radstressdir[jj]=THIRD*pr[PRAD0]*(4.0*q->ucon[dir]*q->ucov[jj] + delta(dir,jj));
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){
    DLOOPA(jj) radstressdir[jj]=THIRD*pr[PRAD0]*(4.0*q->uradcon[dir]*q->uradcov[jj] + delta(dir,jj));
  }
  else DLOOPA(jj) radstressdir[jj]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.


}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* takes E and F^i from primitives and calculates radiation stress ****
//******* tensor R^ij using M1 closure scheme *****************************
//**********************************************************************
// Use: For initial conditions and dumping
// pp : fluid frame orthonormal basis radiation conserved quantities
// Rij : fluid frame orthonormal basis radiation stress-energy tensor
int calc_Rij_ff(FTYPE *pp, FTYPE Rij[][NDIM])
{
  FTYPE E=pp[PRAD0];
  FTYPE F[3]={pp[PRAD1],pp[PRAD2],pp[PRAD3]};

  FTYPE nx,ny,nz,nlen,f;

  nx=F[0]/E;
  ny=F[1]/E;
  nz=F[2]/E;
  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
  if(EOMRADTYPE==EOMRADEDD){
	f=1./3.; // f and Rij are both as if nx=ny=nz=0
	//  f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){

	if(nlen>=1.) f=1.; // KORALTODO: limiter, but only used so far for IC
	else  f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  //M1
  }
  else if(EOMRADTYPE==EOMRADNONE){

  }
  else{
	dualfprintf(fail_file,"1 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
	myexit(837453242);
  }
  
  ////////// Get R^{ij} in orthonormal fluid frame 
  Rij[0][0]=E;

  if(EOMRADTYPE==EOMRADEDD){
	// KORALTODO: Below 3 should be zero for Eddington approximation!  Why didn't original koral have that?
	Rij[0][1]=Rij[1][0]=0.0;
	Rij[0][2]=Rij[2][0]=0.0;
	Rij[0][3]=Rij[3][0]=0.0;
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){
	Rij[0][1]=Rij[1][0]=F[0];
	Rij[0][2]=Rij[2][0]=F[1];
	Rij[0][3]=Rij[3][0]=F[2];
  }
  else if(EOMRADTYPE==EOMRADNONE){

  }
  else{
	dualfprintf(fail_file,"2 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
	myexit(837453243);
  }


  // normalize n^i for Rij calculation
  if(nlen>0){
	nx/=nlen;
	ny/=nlen;
	nz/=nlen;
  }
  else{
	;
  }

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}



FTYPE my_min(FTYPE a, FTYPE b)
{
  if(a<b) return a;
  else return b;
}

FTYPE my_sign(FTYPE x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 0.;
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 4by4 matrix
int inverse_44matrix(FTYPE a[][NDIM], FTYPE ia[][NDIM])
{
  FTYPE mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  FTYPE	tmp[12]; FTYPE	src[16]; FTYPE det;
  /* transpose matrix */
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

 
  /* calculate matrix inverse */
  det = 1.0/det; 

  if(isnan(det)){
    dualfprintf(fail_file,"det in inverse 4x4 zero\n");
	return(1); // indicates failure
	//    myexit(13235);
  }

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      ia[i][j]= dst[i*4+j];

  return 0;
}






/*********************************************************************************/
/****** radiative ortonormal ff primitives (E,F^i) <-> primitives in lab frame  *******/
// Used only for initial conditions
/*********************************************************************************/
// whichvel: input vel type for U1-U3
// whichcoord: input coord type for both U1-U3 and URAD1-URAD3
// whichdir: LAB2FF or FF2LAB  . In addition, here lab means HARM-lab different by alpha factor from true lab.
// i,j,k,loc = standard grid location
// ptrgeom: any input geometry for the lab frame (ptrgeom could be from MCOORD, PRIMECOORDS, etc.) (same for pin's velocity as well as orthonormal basis)
//          If ptrgeom==NULL, then use i,j,k,loc to get geometry in whichcoord coordinates
// pradffortho: radiation primitives (PRAD0-3) should be fluid-frame orthonormal basis values (i.e. E,F in fluid frame orthonormal basis)
// pin: inputs for primitives (i.e. whichvel for U1-U3 and whichcoord for U1-U3,URAD1-URAD3)
// pout: outputs for normal primitives (i.e. WHICHVEL for U1-U3,URAD1-URAD3 but still in whichcoord coordinates)
int prad_fforlab(int *whichvel, int *whichcoord, int whichdir, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{
  FTYPE Rijff[NDIM][NDIM],Rijlab[NDIM][NDIM],U[NPR]={0};
  int pliter,pl;
  int primcoord;
  int jj,kk;
  struct of_geom geomtousedontuse;
  struct of_geom *ptrgeomtouse=&geomtousedontuse;


  if(ptrgeom==NULL){
    if(*whichcoord!=PRIMECOORDS){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_genloc(getprim,*whichcoord,i,j,k,loc,ptrgeomtouse);
    }
    else{
      get_geometry(i, j, k, loc, ptrgeomtouse);
    }
  }
  else{
    // then assumes ptrgeom is in *whichcoord coordinates
    ptrgeomtouse=ptrgeom;
  }


  // set primitive that can use as pre-existing fluid velocity if need to use for reduction
  // also use pout instead of pin so preserves pin no matter what (unless user set pin=pout)
  PLOOP(pliter,pl) pout[pl]=pin[pl];


  // radiative stress tensor in the fluid frame orthonormal basis
  // assuming input pradffortho for radiation is in fluid frame orthonormal basis, but in "primitive" format so using pradffortho[PRAD0-PRAD3]
  // gets R^{ij} in fluid frame orthonormal basis from primitive quantities in fluid frame orthonormal basis
  calc_Rij_ff(pradffortho,Rijff);
  
  //  PLOOPRADONLY(pl) dualfprintf(fail_file,"pl=%d pout=%g\n",pl,pout[pl]);
  //  DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijff=%g\n",jj,kk,Rijff[jj][kk]);


  // get ucon (assumed primitive velocity in ptrgeomtouse coordinates)
  FTYPE ucon[NDIM],others[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(*whichvel,pout,ptrgeomtouse,ucon,others);


  // also convert whichvel ucon to WHICHVEL primitive velocity for use by u2p_rad() and as needed for consistent final output from this function and as possible backup value
  if(*whichvel!=WHICHVEL) ucon2pr(WHICHVEL,ucon,ptrgeomtouse,pout);

  
  // transform and boost (ultimately converts pradffortho -> Rijff -> Rijlab -> U)
  if(whichdir==FF2LAB){
	int tconcovtypeA=TYPEUCON;
	int tconcovtypeB=TYPEUCON;
	if(*whichcoord==PRIMECOORDS) primcoord=1;
	else primcoord=0;
	tensor_lab2orthofluidorback(primcoord, FF2LAB, ptrgeomtouse, TYPEUCON, ucon, tconcovtypeA, tconcovtypeB, Rijff, Rijlab);

	//R^munu -> R^mu_nu so in standard form to extract conserved quantity R^t_\nu
	indices_2221(Rijlab,Rijlab,ptrgeomtouse);

    // Store radiation conserved quantity from R^t_\nu .  u2p_rad() below only uses radiation U's.
    // for true lab to fake-harm lab, end up dividing by alpha (see vector_harm2orthofluidorback() in tetrad.c)
    FTYPE alpha=ptrgeomtouse->alphalapse;
    U[URAD0]=Rijlab[TT][TT]/alpha;
    U[URAD1]=Rijlab[TT][RR]/alpha;
    U[URAD2]=Rijlab[TT][TH]/alpha;
    U[URAD3]=Rijlab[TT][PH]/alpha;

    // 

  }
  else{
	dualfprintf(fail_file,"prad_fforlab() not yet setup for lab2ff since not needed.");
	myexit(652526624);
  }



  //convert to WHICHVEL PRIMECOORDS velocity type primitives
  PFTYPE lpflag=UTOPRIMNOFAIL,lpflagrad=UTOPRIMRADNOFAIL;
  int showmessages=1; // LEAVE on (not normal debugging)
  int allowlocalfailurefixandnoreport=1;
  // NOTEMARK: lpflag=UTOPRIMNOFAIL means accept input pout for velocity to maybe be used in local reductions to fluid frame.
  // u2p_rad() only uses U[URAD0-URAD3]
  // generally u2p_rad() could use all of pout[] except only assigns pout[PRAD0-PRAD3] and doesn't use that for anything except as "static" solution (i.e. uses pin effectively)
  u2p_rad(showmessages, allowlocalfailurefixandnoreport, U, pout, ptrgeomtouse, &lpflag, &lpflagrad);


  // so now both fluid and radiation velocities are in WHICHVEL whichcoord format
  // inversion returns WHICHVEL velocity type, so pass that back
  //  *whichvel=WHICHVEL;



  // no longer force return of WHICHVEL, just change both to whichvel.
  FTYPE uconback[NDIM],othersback[NUMOTHERSTATERESULTS];
  // for fluid
  ucon_calc_whichvel(WHICHVEL,pout,ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,pout);
  // for radiation
  ucon_calc_whichvel(WHICHVEL,&pout[URAD1-U1],ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,&pout[URAD1-U1]);



  // DEBUG:
  if(lpflag!=UTOPRIMNOFAIL || lpflagrad!=UTOPRIMRADNOFAIL){ // DEBUG with 1||
	dualfprintf(fail_file,"Failed to invert during prad_fforlab() with whichdir=%d.  Assuming fixups won't be applied: %d %d\n",whichdir,lpflag,lpflagrad);
	dualfprintf(fail_file,"ijk=%d %d %d : %d\n",ptrgeomtouse->i,ptrgeomtouse->j,ptrgeomtouse->k,ptrgeomtouse->p);
	PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pin=%g U=%g\n",pl,pin[pl],U[pl]);
	DLOOPA(jj) dualfprintf(fail_file,"jj=%d ucon=%g\n",jj,ucon[jj]);
	DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijff=%g Rijlab=%g\n",jj,kk,Rijff[jj][kk],Rijlab[jj][kk]);
	DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d gcov=%g gcon=%g\n",jj,kk,ptrgeomtouse->gcov[GIND(jj,kk)],ptrgeomtouse->gcon[GIND(jj,kk)]);
	PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pout=%g\n",pl,pout[pl]);
	myexit(189235);
	//if(ptrgeomtouse->i==700) myexit(189235);
	// KORALTODO: Check whether really succeeded?  Need to call fixups?  Probably, but need per-cell fixup.  Hard to do if other cells not even set yet as in ICs.  Should probably include fixup process during initbase.c stuff.
  }



  return 0;
} 


// for BCs, to take E[radiation frame] and u^i as radiation primitives in whichvel/whichcoord
// obtains WHICHVEL/PRIMECOORD primitives
int primefluid_EVrad_to_primeall(int *whichvel, int *whichcoord, struct of_geom *ptrgeom, FTYPE *pin, FTYPE *pout)
{
  int pliter,pl;
  int i=ptrgeom->i;
  int j=ptrgeom->j;
  int k=ptrgeom->k;
  int loc=ptrgeom->p;

  // copy over
  PLOOP(pliter,pl) pout[pl]=pin[pl];

  // get metric grid geometry for these ICs
  int getprim=0;
  struct of_geom geomrealdontuse;
  struct of_geom *ptrgeomreal=&geomrealdontuse;
  gset_genloc(getprim,*whichcoord,i,j,k,loc,ptrgeomreal);

  FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(*whichvel,&pout[URAD1-U1],ptrgeomreal,uradcon,othersrad);

  // now convert velocity so in PRIMECOORDS assuming whichcoord=MCOORD
  mettometp_genloc(i,j,k,loc,uradcon);

  if(*whichcoord!=MCOORD){
	dualfprintf(fail_file,"primefluid_EVrad_to_primeall() needs whichcoord (%d) to be MCOORD (%d)\n",whichcoord,MCOORD);
	myexit(87345246);
  }

  // assumed already inputted PRIMECOORDS WHICHVEL for fluid velocity, so no conversion for the fluid velocity

  // now go from ucon[PRIMECOORDS] -> primitive[PRIMECOORDS] for radiation velocity and get WHICHVEL version
  ucon2pr(WHICHVEL,uradcon,ptrgeom,&pout[URAD1-U1]);

  // now all PRIMECOORDS WHICHVEL type assuming ptrgeom inputted PRIMECOORDS version as expected
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  
  return(0);
}


// Input: start with pin [with fluid in whichvel velocity and whichcoordfluid coordinates (PRIMECOORDS or MCOORD) and radiation as E,F in fluid frame orthonormal basis in whichcoordrad coordinates]
// Output: pout [with all WHICHVEL PRIMECOORDS and radiation using velocity primitive]
//
// Useful for BCs when have (say) VEL3,MCOORD for fluid velocity as well as E,F in ff for radiation and need normal WHICHVEL PRIMECOORDS fluid velocity as well as normal velocity primitive for radiation
int whichfluid_ffrad_to_primeall(int *whichvel, int *whichcoordfluid, int *whichcoordrad, struct of_geom *ptrgeomprimecoords, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{
  int pliter,pl;
  int i=ptrgeomprimecoords->i;
  int j=ptrgeomprimecoords->j;
  int k=ptrgeomprimecoords->k;
  int loc=ptrgeomprimecoords->p;


  //  PLOOP(pliter,pl) dualfprintf(fail_file,"ijk=%d %d %d pl=%d pin0=%g pout0=%g\n",i,j,k,pl,pin[pl],pout[pl]);

  // prad_fforlab() should only use radiation primitives, but copy all primitives so can form ucon for transformation
  PLOOP(pliter,pl) pout[pl]=pin[pl];

  // 4 cases:
  // rad   fluid
  // PRIME PRIME
  // PRIME other
  // other PRIME
  // other other

  // get real geometry if needed
  struct of_geom geomrealraddontuse;
  struct of_geom *ptrgeomrealrad=&geomrealraddontuse;
  struct of_geom geomrealfluiddontuse;
  struct of_geom *ptrgeomrealfluid=&geomrealfluiddontuse;
  if(*whichcoordrad!=PRIMECOORDS){
    int getprim=0;
    gset_genloc(getprim,*whichcoordrad,i,j,k,loc,ptrgeomrealrad);
  }
  else ptrgeomrealrad=ptrgeomprimecoords;
  if(*whichcoordfluid!=PRIMECOORDS){
    int getprim=0;
    gset_genloc(getprim,*whichcoordfluid,i,j,k,loc,ptrgeomrealfluid);
  }
  else ptrgeomrealfluid=ptrgeomprimecoords;



  // make whichcoord for fluid same as for rad before continuing (make fluid same as radiation by only changing fluid whichcoord)
  if(*whichcoordfluid!=*whichcoordrad){
    FTYPE ucon[NDIM];
    FTYPE others[NUMOTHERSTATERESULTS];

    // pr->ucon for fluid
    if (pr2ucon(*whichvel,pout, ptrgeomrealfluid, ucon) >= 1) FAILSTATEMENT("bounds.koral.c:bl2met2metp2v_genloc() for radiation", "pr2ucon()", 2);

    if(*whichcoordfluid==PRIMECOORDS){ // then radiation is not PRIMECOORDS, so go to radiation coords
      metptomet_genloc(i,j,k,loc,ucon); // now ucon is in MCOORD
      // convert MCOORD->whichcoordrad
      coordtrans(MCOORD,*whichcoordrad,i,j,k,loc,ucon);
    }
    else if(*whichcoordrad==PRIMECOORDS){ // then go to PRIMECOORDS for fluid
      // convert whichcoordfluid->MCOORD
      coordtrans(*whichcoordfluid,MCOORD,i,j,k,loc,ucon);
      mettometp_genloc(i,j,k,loc,ucon); // now ucon is in PRIMECOORDS
    }
    else{ // then neither is PRIMECOORDS, so just transform and skip mettometp() or metptomet()
      // convert whichcoordfluid->whichcoordrad
      coordtrans(*whichcoordfluid,*whichcoordrad,i,j,k,loc,ucon);
    }

    // get whichvel primitive
    ucon2pr(*whichvel,ucon,ptrgeomrealrad,pout);

    // changed fluid to have same whichcoord as radiation (set by radiation), so set that
    *whichcoordfluid=*whichcoordrad;

  }
  else{
    // otherwise, whichcoord same for fluid and radiation, so can continue
  }



  // get WHICHVEL primitives (still will be in whichcoord coordinates)
  // whichvel here is for fluid velocity (prad_fforlab() converts velocity to WHICHVEL for consistency with only currently allowed output of radiation velocity)
  // pradffortho assumed as in orthonormal fluid frame, but coordinates of whichcoordrad
  int whichframedir=FF2LAB; // fluid frame orthonormal to lab-frame
  prad_fforlab(whichvel, whichcoordrad, whichframedir, i, j, k, loc, ptrgeomrealrad, pradffortho, pout, pout);

  // output from prad_fforlab() is always WHICHVEL for both fluid and radiation primitives
  // changed whichvel's, so report that back if needed
  //  *whichvel=WHICHVEL;
  // above change of whichvel no longer true (and anyways, whichvel was changed in prad_fforlab() directly)
 
  // PLOOP(pliter,pl) dualfprintf(fail_file,"ijk=%d %d %d pl=%d pout=%g\n",i,j,k,pl,pout[pl]);

 
  // output from prad_fforlab() not yet necessarily PRIMECOORDS.
  if(*whichcoordrad==MCOORD){
    // Get all primitives in WHICHVEL/PRIMECOORDS (no conversion for WHICHVEL since prad_fforlab() already put quantities in WHICHVEL due to u2p_rad() only setup for WHICHVEL)
    if (bl2met2metp2v_genloc(*whichvel, *whichcoordrad, pout, i,j,k,loc) >= 1){
      FAILSTATEMENT("bounds.koral.c:bound_radatmbeaminflow()", "bl2ks2ksp2v()", 1);
    }
  }

  // changed coordinates to PRIMECOORDS, so set that as the case
  *whichcoordfluid=*whichcoordrad=PRIMECOORDS;

  return(0);

}



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//calculates general Lorenz matrix for lab <-> ff
//whichdir: [LAB2FF, FF2LAB]
// SUPERGODMARK: What is lab frame velocity?
// UNUSED -- KORALTODO :REMOVE?
int calc_Lorentz_laborff(int whichdir,FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE L[][NDIM])
{
  int ii,jj,kk;
  int verbose=0;

  //"to" frame 4-velocity
  FTYPE ucon[NDIM],ucov[NDIM];
  //"from" frame 4-velocity
  FTYPE wcon[NDIM],wcov[NDIM];

  if(whichdir==LAB2FF)
    {
      //fluid frame
      DLOOPA(ii) 
      {
	ucon[ii]=q->ucon[ii];
	ucov[ii]=q->ucov[ii];
      }
      
      //lab frame
      //TODO: KerrShild???
      wcon[0]=1./sqrt(-ptrgeom->gcov[GIND(0,0)]);
      wcon[1]=wcon[2]=wcon[3]=0.;
      indices_21(wcon,wcov,ptrgeom);
    }
  else if(whichdir==FF2LAB)
    {
      //fluid frame
      DLOOPA(ii) 
      {
	wcon[ii]=q->ucon[ii];
	wcov[ii]=q->ucov[ii];
      }

      //lab frame
      ucon[0]=1./sqrt(-ptrgeom->gcov[GIND(0,0)]);
      ucon[1]=ucon[2]=ucon[3]=0.;
      indices_21(ucon,ucov,ptrgeom);
    }

  //temporary Om matrix
  FTYPE Om[NDIM][NDIM];

  DLOOP(ii,jj)
      Om[ii][jj]=ucon[ii]*wcov[jj]-wcon[ii]*ucov[jj];
  
  //Lorentz factor = -w^mu u_mu
  FTYPE gam=0.;
  DLOOPA(ii)
    gam+=wcon[ii]*ucov[ii];

  FTYPE Omsum;
  //Lorentz matrix components
  DLOOP(ii,jj)
    {
      Omsum=0.;
      DLOOPA(kk)
	Omsum+=Om[ii][kk]*Om[kk][jj];
	
      L[ii][jj]=delta(ii,jj)+1./(1.+gam)*Omsum+Om[ii][jj];
    }

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost lab <-> fluid frame
//whichdir: [LAB2FF, FF2LAB]
int boost22_laborff(int whichdir, FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom)
{ 
  int ii,jj,kk,ll;
  FTYPE Tt[NDIM][NDIM];

  //general Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  calc_Lorentz_laborff(whichdir,pp,q,ptrgeom,L);

  //copying
  DLOOP(ii,jj)
    Tt[ii][jj]=T1[ii][jj];
  
  //boosting
  DLOOP(ii,jj)
    {
      T2[ii][jj]=0.;
      DLOOP(kk,ll)
	T2[ii][jj]+=L[ii][kk]*L[jj][ll]*Tt[kk][ll];
    }

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost from lab to fluid frame
// UNUSED -- KORALTODO :REMOVE?
int boost2_laborff(int whichdir, FTYPE A1[NDIM],FTYPE A2[NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom)
{ 
  int ii,jj,kk,ll;
  FTYPE At[NDIM];

  //general Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  calc_Lorentz_laborff(whichdir,pp,q,ptrgeom,L);

  //copying
  DLOOPA(ii)
    At[ii]=A1[ii];
  
  //boosting
  DLOOPA(ii)
    {
      A2[ii]=0.;
      DLOOPA(kk)
	A2[ii]+=L[ii][kk]*At[kk];
    }

  return 0; 
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//simple Lorentz boosts between ortonormal frames
// UNUSED -- KORALTODO :REMOVE?
int boost2_zamo2ff(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(ZAMO2FF, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

// UNUSED -- KORALTODO :REMOVE?
int boost2_ff2zamo(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(FF2ZAMO, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost fluid frame -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
int boost2_fforzamo(int whichdir, FTYPE A1[NDIM],FTYPE A2[NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{
  // ptrgeom not used for now
  int i,j,k,l;
  FTYPE At[NDIM];
  FTYPE *ulab = q->ucon;
  FTYPE uzamo[NDIM];

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  //proper velocity for ZAMO
  FTYPE vpr[NDIM];
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  if(whichdir==ZAMO2FF){
    vpr[0]*=-1.0;
    vpr[1]*=-1.0;
    vpr[2]*=-1.0;
  }

  //Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  
  //Lorentz factor
  FTYPE vpr2=dot3(vpr,vpr); 
  FTYPE gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2; // NOTEMARK: maybe catestrophic cancellation issue
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }


  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }
  

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }


  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost ZAMO -> fluid frame
// UNUSED -- KORALTODO :REMOVE?
int boost22_zamo2ff(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(ZAMO2FF, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

// UNUSED -- KORALTODO :REMOVE?
int boost22_ff2zamo(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(FF2ZAMO, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost fluid frame -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
int boost22_fforzamo(int whichdir, FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  // ptrgeom not used for now
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];
  FTYPE *ulab = q->ucon;
  FTYPE uzamo[NDIM];

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  //proper velocity for ZAMO
  FTYPE vpr[NDIM];
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  if(whichdir==ZAMO2FF){
    vpr[0]*=-1.0;
    vpr[1]*=-1.0;
    vpr[2]*=-1.0;
  }

  //Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  
  //Lorentz factor
  FTYPE vpr2=dot3(vpr,vpr); 
  FTYPE gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2; // NOTEMARK: maybe catestrophic cancellation issue
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }


  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation ZAMO -> lab
// UNUSED -- KORALTODO :REMOVE?
int trans22_zamo2lab(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE elo[][NDIM])
{
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=elo[i][k]*elo[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation lab -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
int trans22_lab2zamo(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE eup[][NDIM])
{
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=eup[k][i]*eup[l][j]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation lab -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
int trans2_lab2zamo(FTYPE *u1,FTYPE *u2,FTYPE eup[][NDIM])
{
  int i,j,k;
  FTYPE ut[NDIM];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*eup[j][i];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation ZAMO -> lab
// UNUSED -- KORALTODO :REMOVE?
int trans2_zamo2lab(FTYPE *u1,FTYPE *u2,FTYPE elo[][NDIM])
{
  int i,j,k;
  FTYPE ut[NDIM];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*elo[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// T^ij -> T^i_j
int indices_2221(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<NDIM;k++)
	    {
	      Tt[i][j]+=T1[i][k]*ptrgeom->gcov[GIND(k,j)];
	    }	  
	}
    }

   for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A^i -> A^_j
int indices_21(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
	{
	  At[i]+=A1[k]*ptrgeom->gcov[GIND(i,k)];
	}	  
    }

   for(i=0;i<NDIM;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A_i -> A^_j
int indices_12(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
	{
	  At[i]+=A1[k]*ptrgeom->gcon[GIND(i,k)];
	}	  
    }

   for(i=0;i<NDIM;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}





//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//**********************************************************************
//**********************************************************************
//
///////////////
//
// Invert U->direct Primitive for radiation
// OLD (i.e. no longer true): (must come after HD or MHD or whatever sets velocity of fluid, because radiation needs to have updated velocity so that can define fluid frame)
// old code inside utoprimgen.c was:
//    struct of_state qrad;
// this uses new pr to get only ucon and ucov
//get_state_uconucovonly(pr, ptrgeom, &qrad); // OLD
// get new radiation primitives
//
// NEW (currently true): fluid frame no longer needed because go directly from lab-frame conserved quantities to lab-frame primitive quantities.
//
// ERADLIMIT applied here, but (SUPERGODMARK KORALTODO) should in general treat Erf<something as failure in inversion and try to average-out first before reducing to ERf=something and zero relative velocity.  And not sure what that something should be except as relative to other energy densities like how handle UU and BSQ relative to RHO.  Should do that instead and just have very small floor. 
//
// uu: Conserved quantities with URAD0,1,2,3 as radiation conserved quantities
// pp: primitives with PRAD0,1,2,3 as radiation primitive quantities
// ptrgeom: Standard pointer to geometry
// lpflag: see gobal.nondepmnemonics.h .  Tells u2p_rad() if can use/trust fluid velocity.
// lpflagrad: Should be set to indicate success of u2p_rad() inversion
//
// NOTES:
//
// Using *lpflag<=UTOPRIMNOFAIL to check for fluid inversion success rather than a SOFTer condition (e.g. no fail or IFUTOPRIMFAILSOFT==1) because only want to trust fluid as reduction of M1 in case where velocity is accurate with non-negative densities.
//
///////////////
int u2p_rad(int showmessages, int allowlocalfailurefixandnoreport, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  int jj,kk;
  FTYPE pp[NPR];
  int pliter,pl;


  if(WHICHVEL!=VELREL4){
    dualfprintf(fail_file,"u2p_rad() only setup for relative 4-velocity, currently.\n");
    myexit(137432636);
  }


  // copy over pin so pin isn't modified until end
  PLOOP(pliter,pl) pp[pl]=pin[pl];

  //////////////////////
  //
  // Prepare inversion from U->p for radiation assuming M1 closure
  //
  //////////////////////

  *lpflagrad=UTOPRIMRADNOFAIL;


  //conserved - R^t_mu
  FTYPE Av[NDIM]={uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3]};
  //indices up - R^tmu
  indices_12(Av,Av,ptrgeom);


  FTYPE gammarel2,delta,numerator,divisor;
  FTYPE Erf;
  FTYPE urfconrel[NDIM];


  if(EOMRADTYPE==EOMRADEDD){
	// NOTEMARK: Can't use normal inversion that assumes R^t_i are independently evolved because they will generally lead to different velocity than fluid.

	// radiation is same as fluid gamma (assume fluid has already been inverted)
	urfconrel[1]=pp[PRAD1]=pp[U1];
	urfconrel[2]=pp[PRAD2]=pp[U2];
	urfconrel[3]=pp[PRAD3]=pp[U3];
	
	// get gammarel2
	FTYPE gammarel,qsq;
	gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammarel,&qsq);
	gammarel2=gammarel*gammarel;
	
    FTYPE alpha=ptrgeom->alphalapse; //sqrtl(-1./ptrgeom->gcon[GIND(0,0)]);
	// get energy density in fluid frame from lab-frame
	Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){

    // get \gamma^2 for relative 4-velocity
    get_m1closure_gammarel2(showmessages,ptrgeom,Av,&gammarel2,&delta,&numerator,&divisor);

    if(0){
      // testing
      FTYPE Avnew[NDIM]={Av[0],Av[1],Av[2],Av[3]};
      FTYPE urfconrelnew[NDIM];
      FTYPE gammarel2new,deltanew,numeratornew,divisornew,Erfnew;
      //      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avnew,&gammarel2new,&deltanew,&numeratornew,&divisornew,&Erfnew,urfconrelnew);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avnew,NULL,&deltanew,&numeratornew,&divisornew,&Erfnew,urfconrelnew);
    }



    // get E in radiation frame
    get_m1closure_Erf(ptrgeom,Av,gammarel2,&Erf);

    // get relative 4-velocity
    if(CASECHOICE==JONCHOICE) get_m1closure_urfconrel(showmessages,allowlocalfailurefixandnoreport,ptrgeom,pp,Av,gammarel2,delta,numerator,divisor,&Erf,urfconrel,lpflag,lpflagrad);
    else if(CASECHOICE==OLEKCHOICE) get_m1closure_urfconrel_olek(showmessages,allowlocalfailurefixandnoreport,ptrgeom,pp,Av,gammarel2,delta,&Erf,urfconrel,lpflag,lpflagrad);


  }// end if M1
  else{
	dualfprintf(fail_file,"No such EOMRADTYPE=%d in u2p_rad()\n",EOMRADTYPE);
	myexit(368322162);
  }


  //new primitives (only uses urfcon[1-3])
  pin[PRAD0]=Erf;
  pin[PRAD1]=urfconrel[1];
  pin[PRAD2]=urfconrel[2];
  pin[PRAD3]=urfconrel[3];

  //  DLOOPA(jj){
  //    if(!isfinite(pin[PRAD0+jj])){
  //      dualfprintf(fail_file,"caughtnan: jj=%d : ijk=%d %d %d\n",jj,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  //    }
  //  }

  if(DORADFIXUPS==1 || allowlocalfailurefixandnoreport==0){
	// KORALTODO: Problem is fixups can average across shock or place where (e.g.) velocity changes alot, and averaging diffuses shock and can leak-out more failures.
  }
  else{
	// CASE reductions (so set as no failure so fixups don't operate -- but might also want to turn off CHECKINVERSIONRAD else that routine won't know when to ignore bad U->P->U cases.)
	*lpflagrad=UTOPRIMRADNOFAIL;
  }

  return 0;
}





// interpolate between optically thick and thin limits when no u2p_rad() inversion solution
// KORALTODO SUPERGODMARK: Maybe should reduce to fluid frame or gammamax frame as linearly determined by opacity.
// I.e. if optically thick, then reduce to fluid frame.  If optically thin, reduce to gammamax.
//      else if(M1REDUCE==TOOPACITYDEPENDENTFRAME) opacity_interpolated_urfconrel(tautotmax,pp,ptrgeom,Av,Erf,gammarel2,&Erf,urfconrel);
static int opacity_interpolated_urfconrel(FTYPE tautotmax, FTYPE *pp,struct of_geom *ptrgeom,FTYPE *Av, FTYPE Erf,FTYPE gammarel2,  FTYPE *Erfnew, FTYPE *urfconrel)
{
  int jj;
  FTYPE alpha=ptrgeom->alphalapse; //sqrtl(-1./ptrgeom->gcon[GIND(0,0)]);

  //  dualfprintf(fail_file,"Erf=%g gammarel2=%g\n",Erf,gammarel2);

  FTYPE gammafluid,gammarel2fluid,qsqfluid,Erffluid;
  gamma_calc_fromuconrel(&pp[U1-1],ptrgeom,&gammafluid,&qsqfluid);
  gammarel2fluid=gammafluid*gammafluid;
  get_m1closure_Erf(ptrgeom, Av, gammarel2fluid, &Erffluid);
  if(Erffluid<ERADLIMIT) Erffluid=ERADLIMIT;

  FTYPE gammarad,gammarel2rad,qsqrad,Erfrad;
  gamma_calc_fromuconrel(&pp[URAD1-1],ptrgeom,&gammarad,&qsqrad);
  gammarel2rad=gammarad*gammarad;
  get_m1closure_Erf(ptrgeom, Av, gammarel2rad, &Erfrad);
  if(Erfrad<ERADLIMIT) Erfrad=ERADLIMIT;

  // now set urfconrel.  Choose fluid if tautotmax>=2/3 (updated fluid value), while choose previous radiation value (i.e. static!)
  FTYPE tautotmaxlim=MIN(fabs(tautotmax),1.0); // limit for interpolation below
  // done with Erf, so get Erfnew (Erf and *Erfnew might be same variable, but Erf passed by value so changing *Erfnew won't change Erf anyways)
  *Erfnew = (1.0-tautotmaxlim)*Erfrad + tautotmaxlim*Erffluid;
  SLOOPA(jj) urfconrel[jj] = (1.0-tautotmaxlim)*pp[URAD1+jj-1] + tautotmaxlim*pp[U1+jj-1];

  //  dualfprintf(fail_file,"i=%d tautotmax=%g tautotmaxlim=%g\n",ptrgeom->i,tautotmax,tautotmaxlim);
  //  SLOOPA(jj) dualfprintf(fail_file,"jj=%d Erfrad=%g Erffluid=%g gammarad=%g gammafluid=%g Erfnew=%g urfconrel=%g\n",jj,Erfrad,Erffluid,gammarad,gammafluid,*Erfnew,urfconrel[jj]);

  return(0);
}



// get's gamma^2 for lab-frame gamma
static int get_m1closure_gammarel2(int showmessages, struct of_geom *ptrgeom, FTYPE *Av, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn)
{
  FTYPE gamma2,gammarel2,delta,numerator,divisor;

  if(1){
    // has some catastrophic cancellation issue for non-moving velocity at very low E\sim 1E-92 (as in RADPULSE test if no temperature conversion)

    //g_munu R^tmu R^tnu
    int jj,kk;
    FTYPE gRR=0.0;
    DLOOP(jj,kk) gRR += ptrgeom->gcov[GIND(jj,kk)]*Av[jj]*Av[kk];

    //the quadratic equation for u^t of the radiation rest frame (urf[0])
    // Formed as solution for solving two equations (R^{t\nu} R^t_\nu(E,ut) and R^{tt}(E,ut)) for ut
    //supposed to provide two roots for (u^t)^2 of opposite signs
    FTYPE a,b,c;
    a=16.*gRR;
    b=8.*(gRR*ptrgeom->gcon[GIND(0,0)]+Av[0]*Av[0]);
    c=ptrgeom->gcon[GIND(0,0)]*(gRR*ptrgeom->gcon[GIND(0,0)]-Av[0]*Av[0]);
    delta=b*b-4.*a*c;

    numerator=0.5*(-b-sqrt(delta));
    divisor=a;

    gamma2=numerator/divisor; // lab-frame gamma^2
    //if unphysical try the other root
    if(gamma2<=0.){
      numerator=0.5*(-b+sqrt(delta));
      divisor=a;
      gamma2=  numerator/divisor; 
    }
    
    *numeratorreturn=numerator;
    *divisorreturn=divisor;
  }
  //    dualfprintf(fail_file,"GAMMA2CHECK: ijk=%d %d %d : %g %g : a=%g b=%g c=%g : delta=%g gRR=%g Av0123=%g %g %g %g : gamma2=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,0.5*(-b-sqrt(delta))/a,0.5*(-b+sqrt(delta))/a,a,b,c,delta,gRR,Av[0],Av[1],Av[2],Av[3],gamma2);


  else{
    // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
    // well, then had problems for R~1E-14 for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
    // see u2p_inversion.nb
    static FTYPE gctt, gv11, gv12,  gv13,  gv14,  gv22,  gv23,  gv24,  gv33,  gv34,  gv44,  Rtt,  Rtx,  Rty,  Rtz;
    gv11=ptrgeom->gcov[GIND(0,0)];
    gv12=ptrgeom->gcov[GIND(0,1)];
    gv13=ptrgeom->gcov[GIND(0,2)];
    gv14=ptrgeom->gcov[GIND(0,3)];
    gv22=ptrgeom->gcov[GIND(1,1)];
    gv23=ptrgeom->gcov[GIND(1,2)];
    gv24=ptrgeom->gcov[GIND(1,3)];
    gv33=ptrgeom->gcov[GIND(2,2)];
    gv34=ptrgeom->gcov[GIND(2,3)];
    gv44=ptrgeom->gcov[GIND(3,3)];
    Rtt=Av[0];
    Rtx=Av[1];
    Rty=Av[2];
    Rtz=Av[3];
    gctt=ptrgeom->gcon[GIND(0,0)];

    delta = (1. + 3.*gctt*gv11)*((Rtt)*(Rtt)) + 
      6.*gctt*Rtt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 
      3.*gctt*(gv22*((Rtx)*(Rtx)) + 2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 
               2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*((Rtz)*(Rtz)));

    divisor=(gv11*((Rtt)*(Rtt)) + 2.*gv12*Rtt*Rtx + gv22*((Rtx)*(Rtx)) + 2.*gv13*Rtt*Rty + 
       2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 2.*(gv14*Rtt + gv24*Rtx + gv34*Rty)*Rtz + 
             gv44*((Rtz)*(Rtz)));

    numerator=(-0.25*((1. + gctt*gv11)*((Rtt)*(Rtt)) + 
                     gctt*(gv22*((Rtx)*(Rtx)) + 2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 
                           2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*((Rtz)*(Rtz))) + 
                     Rtt*(2.*gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 
                          Sqrt(delta))));

    gamma2 = numerator/divisor;
  }


  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  FTYPE alpha=ptrgeom->alphalapse;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  // check for machine error away from 1.0 that happens sometimes
  if(gammarel2>GAMMASMALLLIMIT && gammarel2<1.0){
    //	if(debugfail>=2) dualfprintf(fail_file,"Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);
    gammarel2=1.0;
  }

  //  dualfprintf(fail_file,"gammarel2=%g gamma2=%g delta=%26.20g\n",gammarel2,gamma2,delta);

  *gammarel2return=gammarel2;
  *deltareturn=delta;
  *numeratorreturn=numerator;
  *divisorreturn=divisor;
  return(0);

}




// get Erf
static int get_m1closure_Erf(struct of_geom *ptrgeom, FTYPE *Av, FTYPE gammarel2, FTYPE *Erfreturn)
{
  FTYPE alpha=ptrgeom->alphalapse;

  ////////////
  //
  // get initial attempt for Erf
  // If delta<0, then gammarel2=nan and Erf<RADLIMIT check below will fail as good.
  //
  ////////////
  *Erfreturn = 3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  return(0);
}


// 0 or 1
#define TRYCOLD 0
// COLD leads to some bombs/issues near boundaries.  Not sure why.

// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Av, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  FTYPE Erf=*Erfreturn; // get initial Erf
  FTYPE gammamax=GAMMAMAXRAD;
  int jj,kk;


  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  //////////////////////
  //
  // First case is if gammarel>gammamax, then set gammarel=gammamax unless Erf<ERADLIMIT (~0) in which case set Erf=ERADLIMIT and gammarel=1.
  // Note, can't set urfcon[0]=gammamax in case gammamax still remains space-like, e.g. inside horizon if gammamax isn't big enough.
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  int nonfailure=gammarel2>=1.0 && Erf>ERADLIMIT && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADLIMIT (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Av[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25 && delta>=0.0 && divisor!=0.0) || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADLIMIT;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=gammarel2>gammamax*gammamax && Erf>=ERADLIMIT || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;



  if(nonfailure){
    // get good relative velocity
    FTYPE gammarel=sqrt(gammarel2);
    FTYPE alpha=ptrgeom->alphalapse;

    SLOOPA(jj) urfconrel[jj] = alpha * (Av[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);

    *Erfreturn=Erf; // pass back new Erf to pointer
    return(0);

    //        dualfprintf(fail_file,"NO failure: %g %g ijk=%d %d %d\n",Erf,gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(failure1){
    if(TRYCOLD){
      gammarel2=pow(1.0+10.0*NUMEPSILON,2.0);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Av,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
    }
    else{
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;
     
      SLOOPA(jj) urfconrel[jj] = 0.0; // consistent with gammarel2=1
    }
    if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE3A;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE3A: normal gamma, but Erf<ERADLIMIT. ijk=%d %d %d : %ld %d %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);

  }    
  else if(failure2){
    if(TRYCOLD){
      gammarel2=pow(1.0+10.0*NUMEPSILON,2.0);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Av,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
    }
    else{
      FTYPE gammarel2orig=gammarel2;
      // override
      gammarel2=1.0;
      FTYPE gammarel=1.0;  // use this below

      // get new Erf(gammarel)
      get_m1closure_Erf(ptrgeom, Av, gammarel2, &Erf);
      if(Erf<ERADLIMIT)  Erf=ERADLIMIT;
    
      SLOOPA(jj) urfconrel[jj] = 0.0;

    }
    if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2A;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2A: normal gamma, but Erf<ERADLIMIT. ijk=%d %d %d : %ld %d %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
  }
  else{
    if(TRYCOLD){
      gammarel2=gammamax*gammamax;
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Av,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
      if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
    }
    else{
      FTYPE gammarel=gammamax;
      gammarel2=gammamax*gammamax;

      // get new Erf(gammarel)
      get_m1closure_Erf(ptrgeom, Av, gammarel2, &Erf);


      // Check if Erf is too small with gamma->gammamax
      if(Erf<ERADLIMIT){
        if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1A;
        // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
        Erf=ERADLIMIT;

        // can't use normal velocity with small Erf -- fails with inf or nan
        SLOOPA(jj) urfconrel[jj] = 0.0;

        if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
      }
      else{
        // if Erf normal, assume ok to have gammamax for radiation.  This avoids fixups, which can generate more oscillations.
        if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammarel2=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);

        // regardless of Erf value, now that have some Erf, ensure gamma=gammamax
        // lab-frame radiation relative 4-velocity
        FTYPE alpha=ptrgeom->alphalapse;
        SLOOPA(jj) urfconrel[jj] = alpha * (Av[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
        
        // compute \gammarel using this (gammatemp can be inf if Erf=ERADLIMIT, and then rescaling below will give urfconrel=0 and gammarel=1
        FTYPE gammatemp,qsqtemp;
        int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);
        MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);

        if(!isfinite(gammatemp)){
          SLOOPA(jj) urfconrel[jj] =0.0;
        }
        else if(0&&gammatemp<=gammamax){
          // do nothing, don't make gamma larger just to get consistency
        }
        else{
          // now rescale urfconrel[i] so will give desired \gammamax
          SLOOPA(jj) urfconrel[jj] *= (gammamax/gammatemp);
        }
	
#if(PRODUCTION==0)
        // check that gamma really correctly gammamax
        FTYPE gammatemp2,qsqtemp2;
        MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp2,&qsqtemp2),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
        if(showmessages) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammamax=%g gammatemp=%g gammatemp2=%g ijk=%d %d %d : %ld %d %g\n",gammamax,gammatemp,gammatemp2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
#endif
      }
      //		  if(showmessages && debugfail>=2) DLOOPA(jj) dualfprintf(fail_file,"CASE1B: urfconrel[%d]=%g uu[%d]=%g\n",jj,urfconrel[jj],jj,uu[URAD0+jj]);

      //    SLOOPA(jj) urfconrel[jj] = 0.0; // consistent with gammarel2=1
    }
  }


  // if here, then one failure mode.  See if optically thick or thin and reduce to (e.g.) fluid frame if thick

  // can't use normal velocity with small Erf -- fails with inf or nan
  // setup "old" pp in case used
  pp[PRAD0] = Erf;
  SLOOPA(jj) pp[PRAD1+jj-1] = urfconrel[jj];


  FTYPE tautot[NDIM],tautotmax;
  if(M1REDUCE==TOOPACITYDEPENDENTFRAME){
    // then will possibly need tautotmax
    // get tautot based upon previous pp in order to determine what to do in case of failure
    calc_tautot(pp, ptrgeom, tautot, &tautotmax);
  }

  
  if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
  else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;
  else if(M1REDUCE==TOOPACITYDEPENDENTFRAME) opacity_interpolated_urfconrel(tautotmax,pp,ptrgeom,Av,Erf,gammarel2,&Erf,urfconrel);
     

  *Erfreturn=Erf; // pass back new Erf to pointer
  return(0);
}




// get contravariant relative 4-velocity in lab frame using Olek's koral choices
static int get_m1closure_urfconrel_olek(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Av, FTYPE gammarel2, FTYPE delta, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  FTYPE Erf=*Erfreturn; // get initial Erf
  FTYPE gammamax=GAMMAMAXRAD;
  int jj,kk;


  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  int failure1=gammarel2>1.01*gammamax*gammamax || gammarel2<0. || delta<0.;
  int failure2=gammarel2<1. || delta<0.; // NOTE: first failure1 already catches delta<0.



  if(failure1){
    //    if(failure1 && tautotmax<TAUFAILLIMIT){ // works for DBLSHADOW

    FTYPE gammarel=gammamax;
    gammarel2=gammamax*gammamax;

    // get new Erf(gammarel)
    get_m1closure_Erf(ptrgeom, Av, gammarel2, &Erf);


    // Check if Erf is too small with gamma->gammamax
    if(Erf<ERADLIMIT){
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1A;
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;

      // can't use normal velocity with small Erf -- fails with inf or nan
      SLOOPA(jj) urfconrel[jj] = 0.0;

      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }
    else{
      // if Erf normal, assume ok to have gammamax for radiation.  This avoids fixups, which can generate more oscillations.
      if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);

      // regardless of Erf value, now that have some Erf, ensure gamma=gammamax
      // lab-frame radiation relative 4-velocity
      FTYPE alpha=ptrgeom->alphalapse;
      SLOOPA(jj) urfconrel[jj] = alpha * (Av[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
        
      // compute \gammarel using this (gammatemp can be inf if Erf=ERADLIMIT, and then rescaling below will give urfconrel=0 and gammarel=1
      FTYPE gammatemp,qsqtemp;
      int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);
      MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
        
      // now rescale urfconrel[i] so will give desired \gammamax
      SLOOPA(jj) urfconrel[jj] *= (gammamax/gammatemp);
	
#if(PRODUCTION==0)
      // check that gamma really correctly gammamax
      FTYPE gammatemp2,qsqtemp2;
      MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp2,&qsqtemp2),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
      if(showmessages) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammamax=%g gammatemp=%g gammatemp2=%g ijk=%d %d %d\n",gammamax,gammatemp,gammatemp2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif
    }

  }
  //////////////////////
  //
  // Second case is if gammarel<1 or delta<0, then set gammarel=1.  If Erf<ERADLIMIT (~0), then set Erf=ERADLIMIT and gammarel=1.
  // Can't assume this condition is equivalent to large gamma, because if not, then leads to crazy boost of energy.
  //
  //////////////////////
  else if(failure2){


    FTYPE gammarel2orig=gammarel2;
    // override
    gammarel2=1.0;
    FTYPE gammarel=1.0;  // use this below

    // get new Erf(gammarel)
    get_m1closure_Erf(ptrgeom, Av, gammarel2, &Erf);
    SLOOPA(jj) urfconrel[jj] = 0.0;


    if(Erf<ERADLIMIT){ // JCM
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2A;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2A: gamma<1 or delta<0 and Erf<ERADLIMIT : gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }
    else{
      // normal Erf
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2B: gamma<1 or delta<0 and Erf normal : gammamax=%g gammarel2orig=%21.15g gammarel2=%21.15g delta=%g : i=%d j=%d k=%d\n",gammamax,gammarel2orig,gammarel2,delta,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }


      
  }
  //////////////////////
  //
  // Third case is if no bad conditions, then try regular calculation.  If Erf<ERADLIMIT, then already caught with first condition
  //
  //////////////////////
  else{

    if(Erf<ERADLIMIT){
      Erf=ERADLIMIT;
      SLOOPA(jj) urfconrel[jj] = 0.0;
    }
    
    // get good relative velocity
    FTYPE gammarel=sqrt(gammarel2);
    FTYPE alpha=ptrgeom->alphalapse;

    // NOTEMARK: This overwrites choice above for urfconrel when Erf<ERADLIMIT.
    SLOOPA(jj) urfconrel[jj] = alpha * (Av[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);


    //        dualfprintf(fail_file,"NO failure: %g %g ijk=%d %d %d\n",Erf,gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }

      
  *Erfreturn=Erf; // pass back new Erf to pointer
  return(0);
}









// get's gamma assuming fixed E rather than using original R^{tt} that we assume is flawed near floor regions.  We want to preserve R^{ti} (i.e momentum)
static int get_m1closure_gammarel2_cold(int showmessages, struct of_geom *ptrgeom, FTYPE *Av, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel)
{
  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=ptrgeom->alphalapse;

  static FTYPE gctt, gv11, gv12,  gv13,  gv14,  gv22,  gv23,  gv24,  gv33,  gv34,  gv44,  Rtt,  Rtx,  Rty,  Rtz;
  gv11=ptrgeom->gcov[GIND(0,0)];
  gv12=ptrgeom->gcov[GIND(0,1)];
  gv13=ptrgeom->gcov[GIND(0,2)];
  gv14=ptrgeom->gcov[GIND(0,3)];
  gv22=ptrgeom->gcov[GIND(1,1)];
  gv23=ptrgeom->gcov[GIND(1,2)];
  gv24=ptrgeom->gcov[GIND(1,3)];
  gv33=ptrgeom->gcov[GIND(2,2)];
  gv34=ptrgeom->gcov[GIND(2,3)];
  gv44=ptrgeom->gcov[GIND(3,3)];
  //  Rtt=Av[0];
  Rtx=Av[1];
  Rty=Av[2];
  Rtz=Av[3];
  gctt=ptrgeom->gcon[GIND(0,0)];


  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamax=GAMMAMAXRAD;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);



  // but check if utsq is too small
  FTYPE utsqmina=(0.5*(-8.*gctt*Power(gv12,2)*Power(Rtx,2) + 8.*gv22*Power(Rtx,2) + 
       8.*gctt*gv11*gv22*Power(Rtx,2) - 16.*gctt*gv12*gv13*Rtx*Rty + 
       16.*gv23*Rtx*Rty + 16.*gctt*gv11*gv23*Rtx*Rty - 
       8.*gctt*Power(gv13,2)*Power(Rty,2) + 8.*gv33*Power(Rty,2) + 
       8.*gctt*gv11*gv33*Power(Rty,2) - 16.*gctt*gv12*gv14*Rtx*Rtz + 
       16.*gv24*Rtx*Rtz + 16.*gctt*gv11*gv24*Rtx*Rtz - 
       16.*gctt*gv13*gv14*Rty*Rtz + 16.*gv34*Rty*Rtz + 
       16.*gctt*gv11*gv34*Rty*Rtz - 8.*gctt*Power(gv14,2)*Power(Rtz,2) + 
       8.*gv44*Power(Rtz,2) + 8.*gctt*gv11*gv44*Power(Rtz,2) - 
       1.*Sqrt(Power(8.*gctt*Power(gv12,2)*Power(Rtx,2) - 8.*gv22*Power(Rtx,2) - 
            8.*gctt*gv11*gv22*Power(Rtx,2) + 16.*gctt*gv12*gv13*Rtx*Rty - 
            16.*gv23*Rtx*Rty - 16.*gctt*gv11*gv23*Rtx*Rty + 
            8.*gctt*Power(gv13,2)*Power(Rty,2) - 8.*gv33*Power(Rty,2) - 
            8.*gctt*gv11*gv33*Power(Rty,2) + 16.*gctt*gv12*gv14*Rtx*Rtz - 
            16.*gv24*Rtx*Rtz - 16.*gctt*gv11*gv24*Rtx*Rtz + 
            16.*gctt*gv13*gv14*Rty*Rtz - 16.*gv34*Rty*Rtz - 
            16.*gctt*gv11*gv34*Rty*Rtz + 8.*gctt*Power(gv14,2)*Power(Rtz,2) - 
            8.*gv44*Power(Rtz,2) - 8.*gctt*gv11*gv44*Power(Rtz,2),2) - 
          4.*(16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
             32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
             16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
             32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 
             32.*gv13*gv14*Rty*Rtz - 32.*gv11*gv34*Rty*Rtz + 
             16.*Power(gv14,2)*Power(Rtz,2) - 16.*gv11*gv44*Power(Rtz,2))*
           (Power(gctt,2)*Power(gv12,2)*Power(Rtx,2) + gctt*gv22*Power(Rtx,2) - 
             1.*Power(gctt,2)*gv11*gv22*Power(Rtx,2) + 
             2.*Power(gctt,2)*gv12*gv13*Rtx*Rty + 2.*gctt*gv23*Rtx*Rty - 
             2.*Power(gctt,2)*gv11*gv23*Rtx*Rty + 
             Power(gctt,2)*Power(gv13,2)*Power(Rty,2) + gctt*gv33*Power(Rty,2) - 
             1.*Power(gctt,2)*gv11*gv33*Power(Rty,2) + 
             2.*Power(gctt,2)*gv12*gv14*Rtx*Rtz + 2.*gctt*gv24*Rtx*Rtz - 
             2.*Power(gctt,2)*gv11*gv24*Rtx*Rtz + 
             2.*Power(gctt,2)*gv13*gv14*Rty*Rtz + 2.*gctt*gv34*Rty*Rtz - 
             2.*Power(gctt,2)*gv11*gv34*Rty*Rtz + 
             Power(gctt,2)*Power(gv14,2)*Power(Rtz,2) + gctt*gv44*Power(Rtz,2) - 
             1.*Power(gctt,2)*gv11*gv44*Power(Rtz,2)))))/
   (16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
     32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
     16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
     32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 32.*gv13*gv14*Rty*Rtz - 
     32.*gv11*gv34*Rty*Rtz + 16.*Power(gv14,2)*Power(Rtz,2) - 
    16.*gv11*gv44*Power(Rtz,2));

  FTYPE utsqminb=(0.5*(-8.*gctt*Power(gv12,2)*Power(Rtx,2) + 8.*gv22*Power(Rtx,2) + 
       8.*gctt*gv11*gv22*Power(Rtx,2) - 16.*gctt*gv12*gv13*Rtx*Rty + 
       16.*gv23*Rtx*Rty + 16.*gctt*gv11*gv23*Rtx*Rty - 
       8.*gctt*Power(gv13,2)*Power(Rty,2) + 8.*gv33*Power(Rty,2) + 
       8.*gctt*gv11*gv33*Power(Rty,2) - 16.*gctt*gv12*gv14*Rtx*Rtz + 
       16.*gv24*Rtx*Rtz + 16.*gctt*gv11*gv24*Rtx*Rtz - 
       16.*gctt*gv13*gv14*Rty*Rtz + 16.*gv34*Rty*Rtz + 
       16.*gctt*gv11*gv34*Rty*Rtz - 8.*gctt*Power(gv14,2)*Power(Rtz,2) + 
       8.*gv44*Power(Rtz,2) + 8.*gctt*gv11*gv44*Power(Rtz,2) + 
       Sqrt(Power(8.*gctt*Power(gv12,2)*Power(Rtx,2) - 8.*gv22*Power(Rtx,2) - 
           8.*gctt*gv11*gv22*Power(Rtx,2) + 16.*gctt*gv12*gv13*Rtx*Rty - 
           16.*gv23*Rtx*Rty - 16.*gctt*gv11*gv23*Rtx*Rty + 
           8.*gctt*Power(gv13,2)*Power(Rty,2) - 8.*gv33*Power(Rty,2) - 
           8.*gctt*gv11*gv33*Power(Rty,2) + 16.*gctt*gv12*gv14*Rtx*Rtz - 
           16.*gv24*Rtx*Rtz - 16.*gctt*gv11*gv24*Rtx*Rtz + 
           16.*gctt*gv13*gv14*Rty*Rtz - 16.*gv34*Rty*Rtz - 
           16.*gctt*gv11*gv34*Rty*Rtz + 8.*gctt*Power(gv14,2)*Power(Rtz,2) - 
           8.*gv44*Power(Rtz,2) - 8.*gctt*gv11*gv44*Power(Rtz,2),2) - 
         4.*(16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
            32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
            16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
            32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 
            32.*gv13*gv14*Rty*Rtz - 32.*gv11*gv34*Rty*Rtz + 
            16.*Power(gv14,2)*Power(Rtz,2) - 16.*gv11*gv44*Power(Rtz,2))*
          (Power(gctt,2)*Power(gv12,2)*Power(Rtx,2) + gctt*gv22*Power(Rtx,2) - 
            1.*Power(gctt,2)*gv11*gv22*Power(Rtx,2) + 
            2.*Power(gctt,2)*gv12*gv13*Rtx*Rty + 2.*gctt*gv23*Rtx*Rty - 
            2.*Power(gctt,2)*gv11*gv23*Rtx*Rty + 
            Power(gctt,2)*Power(gv13,2)*Power(Rty,2) + gctt*gv33*Power(Rty,2) - 
            1.*Power(gctt,2)*gv11*gv33*Power(Rty,2) + 
            2.*Power(gctt,2)*gv12*gv14*Rtx*Rtz + 2.*gctt*gv24*Rtx*Rtz - 
            2.*Power(gctt,2)*gv11*gv24*Rtx*Rtz + 
            2.*Power(gctt,2)*gv13*gv14*Rty*Rtz + 2.*gctt*gv34*Rty*Rtz - 
            2.*Power(gctt,2)*gv11*gv34*Rty*Rtz + 
            Power(gctt,2)*Power(gv14,2)*Power(Rtz,2) + gctt*gv44*Power(Rtz,2) - 
            1.*Power(gctt,2)*gv11*gv44*Power(Rtz,2)))))/
   (16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
     32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
     16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
     32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 32.*gv13*gv14*Rty*Rtz - 
     32.*gv11*gv34*Rty*Rtz + 16.*Power(gv14,2)*Power(Rtz,2) - 
    16.*gv11*gv44*Power(Rtz,2));

  dualfprintf(fail_file,"utsq=%g utsqmina=%g utsqminb=%g\n",utsq,utsqmina,utsqminb);
  // KORALTODO: override (only applicable for first root)  Unsure if 2nd root used for GR in ergosphere.  e.g. gv11 switches sign!
  if(utsq<utsqmina && utsqmina>utsqminb) utsq=utsqmina;
  if(utsq<utsqminb && utsqminb>utsqmina) utsq=utsqminb;


  

  // get new Av[0]=R^{tt}
  
  Av[0]=(-1.*(gctt + 4.*utsq)*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
       0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
          36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
           (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
    (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));

  Erf=(-3.*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
       0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
          36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
           (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
    (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));

  dualfprintf(fail_file,"NOR SOL: %g %g :: %g %g %g\n",Av[0],Erf,Rtx,Rty,Rtz);
  
  delta=0; // not yet

  if(1){
    // alt solution
    FTYPE Avalt=(-1.*(gctt + 4.*utsq)*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq - 
       0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
          36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
           (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
      (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));


    FTYPE Erfalt=(-3.*gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) - 12.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
     0.5*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
        36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
         (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq))))/
      (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));

    dualfprintf(fail_file,"ALT SOL: %g %g : %g %g %g\n",Avalt,Erfalt,Rtx,Rty,Rtz);
  }


  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  int jj;
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Av[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  dualfprintf(fail_file,"NORM ROOT 4-vel: %g %g %g\n",urfconrel[1],urfconrel[2],urfconrel[3]);

  
  *Erfreturn=Erf; // pass back new Erf to pointer


  return(0);
}




























//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
// UNUSED -- KORALTODO :REMOVE?
int f_u2prad_num(FTYPE *uu,FTYPE *pp, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM], FTYPE elo[][NDIM],FTYPE *f)
{
  FTYPE Rij[NDIM][NDIM];
  FTYPE ppp[NPR];

  calc_Rij_ff(pp,Rij);
  boost22_ff2zamo(Rij,Rij,pp,q,ptrgeom,eup);
  trans22_zamo2lab(Rij,Rij,elo);
  indices_2221(Rij,Rij,ptrgeom);

  FTYPE gdet=ptrgeom->gdet;

  // f = error
  f[0]=-Rij[0][0]+uu[URAD0];
  f[1]=-Rij[0][1]+uu[URAD1];
  f[2]=-Rij[0][2]+uu[URAD2];
  f[3]=-Rij[0][3]+uu[URAD3];

  return 0;
} 


// U->P inversion for Eddington approximation using Newton method.
// UNUSED -- KORALTODO :REMOVE?
int u2p_rad_num(FTYPE *uu, FTYPE *pp, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM], FTYPE elo[][NDIM])
{
  FTYPE pp0[NPR],pporg[NPR];
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE x[NDIM],f1[NDIM],f2[NDIM],f3[NDIM];
  int i,j,k,iter=0;



  for(i=PRAD0;i<=PRAD3;i++)
    {
      pporg[i]=pp[i];
    }
  
  do
    {
      iter++;
      for(i=PRAD0;i<NPR;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,q,ptrgeom,eup,elo,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+PRAD0]=pp[j+PRAD0]+PRADEPS*pp[PRAD0];
	    
	      f_u2prad_num(uu,pp,q,ptrgeom,eup,elo,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[PRAD0]);

	      pp[j+PRAD0]=pp0[j+PRAD0];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+PRAD0];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}

      for(i=0;i<4;i++)
	{
	  pp[i+PRAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+PRAD0]-pp0[i+PRAD0]);
	  f3[i]=fabs(f3[i]/pp0[PRAD0]);
	}

      if(f3[0]<RADCONV && f3[1]<RADCONV && f3[2]<RADCONV && f3[3]<RADCONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in u2prad_num()\n");
	  
	  for(i=PRAD0;i<NPR;i++)
	    {
	      pp[i]=pporg[i];
	    }
	  
	  return -1;

	  break; // WHY HERE IF UNREACHABLE?
	}
     
    }
  while(1);
  
  if(pp[PRAD0]<ERADLIMIT) 
    {
      printf("%g=PRAD0<ERADLIMIT=%g u2prad()\n",pp[PRAD0],ERADLIMIT);
      pp[PRAD0]=ERADLIMIT;
    }
  

  return 0;
}

 



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives fluid frame -> ZAMO
// Use: Maybe dumping
// pp1 : Full set of primitives with radiation primitives replaced by fluid frame radiation \hat{E} and \hat{F}
// pp2 : Full set of primitives with ZAMO frame for radiation conserved quantities
// UNUSED -- KORALTODO :REMOVE?
int prad_ff2zamo(FTYPE *pp1, FTYPE *pp2, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  FTYPE Rij[NDIM][NDIM];
  int i,j;

  // set all (only non-rad needed) primitives
  int pliter,pl;
  PLOOP(pliter,pl) pp2[pl]=pp1[pl];

  calc_Rij_ff(pp1,Rij);
  boost22_ff2zamo(Rij,Rij,pp1,q,ptrgeom,eup);

  // overwrite pp2 with new radiation primitives
  pp2[PRAD0]=Rij[0][0];
  pp2[PRAD1]=Rij[0][1];
  pp2[PRAD2]=Rij[0][2];
  pp2[PRAD3]=Rij[0][3];

  return 0;
} 

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives ZAMO -> fluid frame - numerical solver
// Use: Error (f) for iteration
// ppff : HD primitives (i.e. WHICHVEL velocity) and radiation primitives of \hat{E} and \hat{F} (i.e. fluid frame conserved quantities)
// ppzamo : \hat{E} & \hat{F} -> E & F in ZAMO
// UNUSED -- KORALTODO :REMOVE?
int f_prad_zamo2ff(FTYPE *ppff, FTYPE *ppzamo, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM],FTYPE *f)
{
  FTYPE Rij[NDIM][NDIM];
  calc_Rij_ff(ppff,Rij);
  boost22_ff2zamo(Rij,Rij,ppff,q,ptrgeom,eup);

  f[0]=-Rij[0][0]+ppzamo[PRAD0];
  f[1]=-Rij[0][1]+ppzamo[PRAD1];
  f[2]=-Rij[0][2]+ppzamo[PRAD2];
  f[3]=-Rij[0][3]+ppzamo[PRAD3];


  return 0;
} 



// SUPERGODMARK: Is q constant ok?
// numerical iteration to find ff from zamo
// Use: Initial conditions
// ppzamo : E & F in ZAMO -> \hat{E} & \hat{F}
// UNUSED -- KORALTODO :REMOVE?
int prad_zamo2ff(FTYPE *ppzamo, FTYPE *ppff, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  FTYPE pp0[NPR],pp[NPR];
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE x[NDIM],f1[NDIM],f2[NDIM],f3[NDIM];
  int i,j,k,iter=0;


  
  //initial guess
  for(i=0;i<NPR;i++)
    {
      pp[i]=ppzamo[i];
    }


  //debug
  f_prad_zamo2ff(ppzamo,pp,q,ptrgeom,eup,f1);
  for(i=0;i<4;i++)
    {
      x[i]=pp[i+PRAD0];
    }  
 
  do
    {
      iter++;
      for(i=PRAD0;i<NPR;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_prad_zamo2ff(pp,ppzamo,q,ptrgeom,eup,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+PRAD0]=pp[j+PRAD0]+PRADEPS*pp[PRAD0];
	    
	      f_prad_zamo2ff(pp,ppzamo,q,ptrgeom,eup,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[PRAD0]);

	      pp[j+PRAD0]=pp0[j+PRAD0];
	    }
	}
  
      //inversion
      inverse_44matrix(J,iJ);

      //updating unknowns
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+PRAD0];
	}      

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}

      for(i=0;i<4;i++)
	{
	  pp[i+PRAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+PRAD0]-pp0[i+PRAD0]);
	  f3[i]=fabs(f3[i]/pp0[PRAD0]);
	}

      if(f3[0]<PRADCONV && f3[1]<PRADCONV && f3[2]<PRADCONV && f3[3]<PRADCONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in prad_zamo2ff()\n");
	  break;
	}
    }
  while(1);



  //returning prad
  for(i=0;i<NPR;i++)
    {
      ppzamo[i]=pp[i];
    }

  return 0;
}











//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
int calc_tautot(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *tautot, FTYPE *tautotmax)
{
#if(0)
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  FTYPE rho=pp[0];
  FTYPE u=pp[1];  
  FTYPE pr=(gamideal-1.)*(u);
  FTYPE T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  FTYPE kappa=calc_kappa_user(rho,T,xx[1],xx[2],xx[3]);
  FTYPE chi=kappa+calc_kappaes_user(rho,T,xx[1],xx[2],xx[3]);
#endif
  //xx[0] holds time
  FTYPE kappa,kappaes,chi;
  calc_kappa(pp,ptrgeom,&kappa);
  calc_kappaes(pp,ptrgeom,&kappaes);
  chi=kappa+kappaes;
  int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1}; // want to ignore non-used dimensions

  int jj;
  *tautotmax=0.0;
  SLOOPA(jj){
    tautot[jj]=chi * (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])))*NxNOT1[jj];
    *tautotmax=MAX(*tautotmax,tautot[jj]);
  }

  return 0;
}

//**********************************************************************
//******* calculates abs opacity over dx[] ***************************
//**********************************************************************
int calc_tauabs(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *tauabs, FTYPE *tauabsmax)
{
#if(0)
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  FTYPE rho=pp[0];
  FTYPE u=pp[1];  
  FTYPE pr=(gamideal-1.)*(u);
  FTYPE T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  FTYPE kappa=calc_kappa_user(rho,T,xx[1],xx[2],xx[3]);
#endif
  FTYPE kappa;
  calc_kappa(pp,ptrgeom,&kappa);

  int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1}; // want to ignore non-used dimensions

  int jj;
  *tauabsmax=0.0;
  SLOOPA(jj){
    tauabs[jj]=kappa * (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])))*NxNOT1[jj];
    *tauabsmax=MAX(*tauabsmax,tauabs[jj]);
  }

  return 0;
}






//**********************************************************************
//suplementary routines for conversions
//**********************************************************************
FTYPE calc_PEQ_ufromTrho(FTYPE T,FTYPE rho)
{
#if(0)
  //  FTYPE p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;
  FTYPE p=rho*T; // /MU_GAS;
  FTYPE u=p/(gamideal-1.);
#else
  // if use local function function instead of below directly,
  // then assume user doesn't care about position for EOS.
  FTYPE u=u_rho0_T_simple(0, 0, 0, CENT, rho, T);
#endif
  return u;
}

FTYPE calc_PEQ_Tfromurho(FTYPE u,FTYPE rho)
{
#if(0)
  FTYPE p=u*(gamideal-1.);
  //  FTYPE T=p/(K_BOLTZ*rho/MU_GAS/M_PROTON);
  //  FTYPE T=p/(rho/MU_GAS);
  FTYPE T=p/(rho);
#else
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
#endif
  return T;
}

// E=urad=arad T^4
FTYPE calc_LTE_EfromT(FTYPE T)
{
  //  return 4.*SIGMA_RAD*T*T*T*T;
  return (ARAD_CODE*T*T*T*T);
}

// E=urad=arad T^4 and just solve for T
FTYPE calc_LTE_TfromE(FTYPE E )
{
  //  return sqrt(sqrt((E/4./SIGMA_RAD)));
  return (sqrt(sqrt((E/ARAD_CODE))));
}


FTYPE calc_LTE_Efromurho(FTYPE u,FTYPE rho)
{
#if(0)
  FTYPE p=(gamideal-1.)*(u);
  //  FTYPE T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  //  FTYPE T=p*MU_GAS/rho;
  FTYPE T=p/rho;
#else
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
#endif

  return (calc_LTE_EfromT(T));
}



// set velocity based upon ncon and gammamax and return in whichvel format for the ptrgeom geometry/coords
int set_ncon_velocity(int whichvel, FTYPE gammamax, FTYPE *ncon, struct of_geom *ptrgeom, FTYPE *uconwhichvel)
{
  int ii,jj;
  FTYPE ncondefault[NDIM]={0.0,-1.0, 0.0, 0.0}; //for radially flowing photons
  FTYPE *nconuse;

  // default is radial motion in zamo frame
  if(ncon==NULL) nconuse=ncondefault;
  else nconuse=ncon;

  // compute \gammarel using nconuse
  FTYPE gammatemp,qsq;
  gamma_calc_fromuconrel(nconuse,ptrgeom,&gammatemp,&qsq);

  // now rescale nconuse[i] so will give desired \gammamax
  FTYPE prtemp[NPR];
  SLOOPA(ii) prtemp[U1+ii-1] = nconuse[ii]*(gammamax/gammatemp);

  // now get u^\mu[lab]
  FTYPE uconlab[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(WHICHVEL,prtemp,ptrgeom,uconlab,others);

  // now get other whichvel type
  FTYPE prtemp2[NPR];
  ucon2pr(whichvel,uconlab,ptrgeom,prtemp2);

  SLOOPA(jj) uconwhichvel[jj] = prtemp2[U1+jj-1];

  return(0);

}
