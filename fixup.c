// all fixup stuff only called for non-B advance


#include "decs.h"


static int simple_average(int startpl, int endpl, int i, int j, int k,PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR]);
static int general_average(int startpl, int endpl, int i, int j, int k, PFTYPE mypflag, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *geom);
static int fixup_nogood(int startpl, int endpl, int i, int j, int k, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom);
static int fixuputoprim_accounting(int i, int j, int k, PFTYPE mypflag, PFTYPE (*lpflag)[NSTORE2][NSTORE3][NUMPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *geom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
static int fixup_negdensities(int *fixed, int startpl, int endpl, int i, int j, int k, PFTYPE mypflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *geom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

static int superdebug_utoprim(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);




/* apply floors to density, internal energy */

// currently called before bound, which assumes bound sets boundary
// values exactly as wanted without any fixing.

#define JONFIXUP 1 // 0=gammie 1=jon's

#if 0
void fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR],int finalstep)
{
  int i, j, k;



  COMPZLOOP{ pfixup(MAC(pv,i,j,k), i, j, k);}

}
#endif

// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
// operations that only need to be done inside computational loop
int pre_fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR])
{

  // grab b^2 flags (above fixup may change u or rho, so must do this after)
  get_bsqflags(stage,pv);


  return(0);
}



// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
// this function actually changes primitives
int post_fixup(int stageit,int finalstep, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{
  int stage,stagei,stagef;
  int boundstage;

#if(UTOPRIMADJUST!=0)
  ////////////////////////////////////
  //
  // utoprim fixup of primitive solution



  if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
  else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
  else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}

  if(SIMULBCCALC>=1) boundstage=STAGE0;
  else boundstage=STAGEM1;


  for(stage=stagei;stage<=stagef;stage++){



    // first bound failure flag
    // OPTMARK: could optimize bound of pflag since often failures don't occur (just ask if any failures first), although probably negligible performance hit
    if(stage<STAGE2){
      bound_pflag(boundstage, finalstep, boundtime, GLOBALPOINT(pflag), USEMPI);
      if(stage!=STAGEM1) boundstage++;
    }


    // check for bad solutions and set as failure if good is reasonably bad
#if(CHECKSOLUTION)
    fixup_checksolution(stage,pv,finalstep);

    // check solution changed pflag, so have to bound again
    if(stage<STAGE2){
      bound_pflag(boundstage, finalstep, boundtime, pflag, USEMPI);
      if(stage!=STAGEM1) boundstage++;
    }

#endif


    // fixup before new solution (has to be here since need previous stage's failure flag)
    fixup_utoprim(stage,pv,pbackup,ucons,finalstep);

    
#if(0)
    // GODMARK: I don't see why need to bound pflag since already done with using pflag
    if(stage<STAGE2){
      if(stage!=STAGEM1){
        bound_pflag(boundstage, finalstep, boundtime, pflag, USEMPI);
        boundstage++;
      }
    }
#endif

  }
#endif

  ////////////////////////////////////
  //
  // standard fixup of floor
  // in this case stageit=-1, so does all stages
  //  fixup(stageit,pv,finalstep);




  return(0);
}


// this function just reports problems, but doesn't fix them
int post_fixup_nofixup(int stageit, int finalstep, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{

  fixup_utoprim_nofixup(STAGEM1,pv,pbackup,ucons,finalstep);

  return(0);
}





#if(JONFIXUP==1)

int fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  int i, j, k;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;



  COMPZLOOP{
    //    dualfprintf(fail_file,"i=%d j=%d k=%d\n",i,j,k); fflush(fail_file);
    get_geometry(i,j,k,CENT,ptrgeom);
    if(fixup1zone(MAC(pv,i,j,k),MAC(ucons,i,j,k), ptrgeom,finalstep)>=1)
      FAILSTATEMENT("fixup.c:fixup()", "fixup1zone()", 1);
  }
  return(0);
}
#else

// GAMMIE OLD FIXUP (not kept up to date)
int fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int finalstep)
{
  int i,j,k,pl,pliter;
  int ip, jp, kp, im, jm, km;
  FTYPE bsq, del;
  FTYPE ftempA,ftempB;
  FTYPE prfloor[NPR];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;




  COMPZLOOP {
    get_geometry(i,j,k, CENT,ptrgeom) ;
    // densities
    if(DOEVOLVERHO||DOEVOLVEUU) set_density_floors(ptrgeom,MAC(pv,i,j,k),prfloor);
  

    if(DOEVOLVERHO ){
      /* floor on density (momentum *not* conserved) */
      if (MACP0A1(pv,i,j,k,RHO) < prfloor[RHO]) {
#if(FLOORDIAGS)
        fladd[RHO] +=
          dVF * ptrgeom->gdet * (prfloor[RHO] - MACP0A1(pv,i,j,k,RHO));
#endif
        MACP0A1(pv,i,j,k,RHO) = prfloor[RHO];
      }
    }
    
    if(DOEVOLVEUU){
      /* floor on internal energy */
      if (MACP0A1(pv,i,j,k,UU) < prfloor[UU]) {
#if(FLOORDIAGS)
        fladd[UU] +=
          dVF * ptrgeom->gdet * (prfloor[UU] - MACP0A1(pv,i,j,k,UU));
#endif
        MACP0A1(pv,i,j,k,UU) = prfloor[UU]; // REBECCAMARK
      }
    }

    /* limit gamma wrt normal observer */
#if(WHICHVEL==VELREL4)
    if(limit_gamma(GAMMAMAX,MAC(pv,i,j,k),MAC(ucons,i,j,k), ptrgeom,-1)>=1)  // need general accounting for entire routine.
      FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
#endif
  }

  return(0);}


#endif



// choose whether within correctable/diagnosticable region
int diag_fixup_correctablecheck(int docorrectucons, struct of_geom *ptrgeom)
{
  int is_within_correctable_region;
  int docorrectuconslocal;

  if(DOONESTEPDUACCOUNTING){
    docorrectuconslocal=docorrectucons;
  }
  else{
    ///////////
    //
    // determine if within correctable region
    //
    ///////////
    if( DOENOFLUX != NOENOFLUX ) {
      is_within_correctable_region=((ptrgeom->i)>=Uconsevolveloop[FIS])&&((ptrgeom->i)<=Uconsevolveloop[FIE])&&((ptrgeom->j)>=Uconsevolveloop[FJS])&&((ptrgeom->j)<=Uconsevolveloop[FJE])&&((ptrgeom->k)>=Uconsevolveloop[FKS])&&((ptrgeom->k)<=Uconsevolveloop[FKE]);
    }
    else{
      is_within_correctable_region=1; // assume diag_fixup() only called where ok to do change to ucons!
    }


    ///////////
    //
    // determine if should do correction to ucons
    // only correct once -- should really put correction somewhere else.
    //
    ///////////
    docorrectuconslocal=docorrectucons  && is_within_correctable_region;

  }

  return(docorrectuconslocal);


}



// record who called the diag_fixup routine
int count_whocalled(struct of_geom *ptrgeom, int finalstep, int whocalled)
{
  int tscale;

  /////////////////////
  //
  // Count times in diag_fixup() and who called diag_fixup()
  // count every time corrects, not just on conserved quantity tracking time
  //
  /////////////////////
  if(DODEBUG){

    if(whocalled>=NUMFAILFLOORFLAGS || whocalled<COUNTNOTHING || ptrgeom->i<-N1BND || ptrgeom->i>N1-1+N1BND ||ptrgeom->j<-N2BND || ptrgeom->j>N2-1+N2BND ||ptrgeom->k<-N3BND || ptrgeom->k>N3-1+N3BND ){
      dualfprintf(fail_file,"In diag_fixup() whocalled=%d for i=%d j=%d k=%d\n",whocalled,ptrgeom->i,ptrgeom->j,ptrgeom->k);
      myexit(24683463);
    }

    if(whocalled!=COUNTNOTHING){
      int indexfinalstep;
      indexfinalstep=0;
      TSCALELOOP(tscale) GLOBALMACP0A3(failfloorcount,ptrgeom->i,ptrgeom->j,ptrgeom->k,indexfinalstep,tscale,whocalled)++;
      if(finalstep){
        indexfinalstep=1;
        // iterate finalstep version
        TSCALELOOP(tscale) GLOBALMACP0A3(failfloorcount,ptrgeom->i,ptrgeom->j,ptrgeom->k,indexfinalstep,tscale,whocalled)++;
      }
    }// end if counting something

  }// end if DODEBUG


  return(0);

}



// compute dU and account
// Assumes Ui,Uf are UDIAG form
// Assumes ucons is UEVOLVE form
int diag_fixup_dUandaccount(FTYPE *Ui, FTYPE *Uf, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep, int whocalled, int docorrectuconslocal)
{
  FTYPE dUincell[NPR];
  int is_within_diagnostic_region;
  FTYPE deltaUavg[NPR],Uiavg[NPR];
  FTYPE Uprefixup[NPR],Upostfixup[NPR];
  int pliter,pl;
  int enerregion;



  /////////////
  //
  // Get deltaUavg[] and also modify ucons if required and should
  //
  /////////////

  if(DOENOFLUX != NOENOFLUX) {  //SASMARKx: adjust the conserved quantity to correspond to the adjusted primitive quanitities
    // Correction to conserved quantities not exactly accurate because using point values where should use averaged values
    // notice that geometry comes after subtractions/additions of EOMs
    UtoU(UDIAG,UEVOLVE,ptrgeom,Ui,Uprefixup);  // convert from UDIAG -> UEVOLVE
    UtoU(UDIAG,UEVOLVE,ptrgeom,Uf,Upostfixup); // convert from UDIAG -> UEVOLVE
          
    PALLLOOP(pl) deltaUavg[pl] = Uf[pl]-Ui[pl];
          
    if(docorrectuconslocal){
      // correct ucons if requested
      //adjust the averaged conserved quantity by the same amt. as the point conserved quantity
      PALLLOOP(pl) ucons[pl] += Upostfixup[pl] - Uprefixup[pl];  

      // old code: UtoU(UDIAG,UEVOLVE,ptrgeom,Uf,ucons); // convert from UNOTHING->returntype (jon's comment)
      // the above line actually converts fixed up U from diagnostic form of U (with gdet) 
      // to evolution form of U (maybe withnogdet) and replaces the avg. conserved quantity (ADT)
    }

  }
  else if(0){

    // this method doesn't work:
    UtoU(UEVOLVE,UDIAG,ptrgeom,ucons,Uiavg); // convert from UNOTHING->returntype

    if(docorrectuconslocal){
      // notice that geometry comes after subtractions/additions of EOMs
      UtoU(UDIAG,UEVOLVE,ptrgeom,Uf,ucons); // convert from UNOTHING->returntype
    }
          
    PALLLOOP(pl) deltaUavg[pl] = Uf[pl]-Uiavg[pl];
  }
  else{
    // original HARM method
    // don't modify ucons

    PALLLOOP(pl) deltaUavg[pl] = Uf[pl]-Ui[pl];
  }



  //only do aggregate accounting, after the fact (just before taking the new time step)
  if(DOONESTEPDUACCOUNTING && docorrectuconslocal < 0 ||  DOONESTEPDUACCOUNTING==0){

    ///////////////////
    //
    // get correction
    //
    //////////////////
    PALLLOOP(pl){
            
      // dUincell means already (e.g.) (dU0)*(\detg')*(dV') = integral of energy in cell = dUint0 in SM
      // So compare this to (e.g.) (U0)*(\detg')*(dV') = U0*gdet*dV in SM
      dUincell[pl]=dVF * deltaUavg[pl];

      if(DOFLOORDIAG){
        // only store this diagnostic once (not for each enerregion)
        // Note that unlike failfloorcount[], failfloordu[] is independent of fladd and fladdterms that are integrated simultaneously rather than in dump_ener.c
        // Also note that failfloordu not stored in restart file, so like spatial debug info it is lost upon restart.
        GLOBALMACP0A1(failfloordu,ptrgeom->i,ptrgeom->j,ptrgeom->k,pl)+=dUincell[pl];
      }

    }// end over pl's




    //////////////
    //
    // Loop over ENERREGIONs
    //
    //////////////
    ENERREGIONLOOP(enerregion){

      // setup pointers to enerregion diagnostics
      enerpos=enerposreg[enerregion];
      fladd=fladdreg[enerregion];
      fladdterms=fladdtermsreg[enerregion];


      ///////////
      //
      // determine if within diagnostic region
      //
      ///////////
      is_within_diagnostic_region=WITHINENERREGION(enerpos,ptrgeom->i,ptrgeom->j,ptrgeom->k);


      /////////////////////////
      //
      // diagnostics (both for enerregion and single-region types)
      //
      /////////////////////////
      if(is_within_diagnostic_region){

        PALLLOOP(pl){
            
          // dUincell means already (e.g.) (dU0)*(\detg')*(dV') = integral of energy in cell = dUint0 in SM
          // So compare this to (e.g.) (U0)*(\detg')*(dV') = U0*gdet*dV in SM
          dUincell[pl]=dVF * deltaUavg[pl];

          fladdterms[whocalled][pl] += (SFTYPE)dUincell[pl];
          fladd[pl] += dUincell[pl];
            

        }// end over pl's
      }// end if within diagnostic region

    }// end over enerregions

  }// end if doing accounting


  return(0);
}


// single call in step_ch.c:post_advance() to do all diag_fixup() diagnostic dU stores.  Still allows counts by other diag_fixup calls.
// for DOONESTEPDUACCOUNTING==1
int diag_fixup_allzones(int truestep, int finalstep, FTYPE (*pf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{

  if(truestep && finalstep){
    int i, j, k, pliter,pl;
    struct of_geom *ptrgeom;
    struct of_geom geomdontuse;

    if(DOENOFLUX!=NOENOFLUX){
      dualfprintf(fail_file,"Cannot use diag_fixup_allzones() with DOENOFLUX==NOENOFLUX\n");
      myexit(3487622211);
    }
    
    ptrgeom=&(geomdontuse);
    
    COMPZLOOP{
      get_geometry(i, j, k, CENT, ptrgeom);
      // account for change of conserved quantities
      // Primitives have been modified by fixup1zone() in advance.c (floors).
      // During this call, called from post_advance(), pf also modified by bounds (poledeath,gammadeath) and also post_fixup() (failures, checks, limits).
      //
      // ucons=unewglobal that stores last steps final substep version of ucum that is full U[]
      // ucons has yet to be modified at all, so is true conserved quantity without corrections (as long as avoided corrections to ucons during other diag_fixup calls).
      // GODMARK: So this method only works if NOENOFLUX==1, since otherwise *need* to modify U[] during modification of p[] since know how much to modify.

      int docorrectucons=-1; // -1 is like 1, but is used to tell if coming from this function (if -1) or not (if 1)
      diag_fixup_Ui_pf(docorrectucons,MAC(ucons,i,j,k),MAC(pf,i,j,k),ptrgeom,finalstep,COUNTONESTEP);
    }
  }    


  return(0);

}



// account for changes by tracking conserved quantities
// accounts for both failures and floor recoveries
// this modifies unew if on finalstep to be consistent with floor-limited primitive
// diagnostics only for actions on conservative quantities
// assume COUNT types are of PFTYPE
int diag_fixup(int docorrectucons, FTYPE *pr0, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep, int whocalled)
{
  struct of_state q;
  FTYPE Uicent[NPR],Ufcent[NPR];
  int failreturn;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  FTYPE deltaUavg[NPR],Uiavg[NPR];
  FTYPE Uprefixup[NPR],Upostfixup[NPR];
  int docorrectuconslocal;



#if(DOSUPERDEBUG)
  superdebug_utoprim(pr0,pr,ptrgeom,whocalled);
  // collect values for non-failed and failed zones
#endif



  // count whocalled diag_fixup()
  count_whocalled(ptrgeom, finalstep, whocalled);



  /////////////////////////////////////////
  //
  // Account for changes in primitives or conserved quantities due to fixups (floor or failures or any other thing that can call diag_fixup()
  //
  // only account if on full timestep
  // ucum (unew) only inverted to primitives on final substep.  Any other conserved or primitive corrections do not matter since they only affected fluxes that go into true conserved quantity that is ucum (unew)
  //
  /////////////////////////////////////////

  if(finalstep > 0){

    ///////////
    // determine if within correctable region
    ///////////
    docorrectuconslocal=diag_fixup_correctablecheck(docorrectucons,ptrgeom);



    ////////////////////////
    //      
    // Get Uicent and Ufcent.  Don't do this inside enerregion because no point since assume diag_fixup() called in limited regions of i,j,k anyways.
    //
    // only account if within active zones for that region
    //
    // Only valid if not higher order method or if MERGED method where conserved (except field) are at points as also the primitives are
    //
    ////////////////////////


    // before any changes
    failreturn=get_state(pr0,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Uicent);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
        

    // after any changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Ufcent);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");

    // if Uicent and Ufcent are both from pi and pf at CENT, then B1,B2,B3 entries are agreeably located even for FLUXB==FLUXCTSTAG


    // Get deltaUavg[] and also modify ucons if required and should
    diag_fixup_dUandaccount(Uicent, Ufcent, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);


  }// end if finalstep>0



  return(0);
}



// like diag_fixup(), but input initial conserved quantity as Ui and final primitive as pf
// Must use this when pi[Ui] doesn't exist and had to use non-hot-MHD inversion.
// Assumes Ui is like unewglobal, so UEVOLVE type
// Assume ultimately hot MHD equations are used, so need to get new Uf that'll differ from Ui
// Also don't know Uf quite yet.
int diag_fixup_Ui_pf(int docorrectucons, FTYPE *Uievolve, FTYPE *pf, struct of_geom *ptrgeom, int finalstep, int whocalled)
{
  struct of_state q;
  FTYPE Ufcent[NPR],Uicent[NPR],ucons[NPR];
  int failreturn;
  int pliter,pl,enerregion;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int docorrectuconslocal;



  // count whocalled diag_fixup()
  count_whocalled(ptrgeom, finalstep, whocalled);


  if(finalstep > 0){

    // determine if within correctable region
    docorrectuconslocal=diag_fixup_correctablecheck(docorrectucons,ptrgeom);



    //////////////////
    //
    // Get ucons
    //
    //////////////////

    // GODMARK: NOENOFLUX==0 not accounted for (have to change unewglobal or something like that)
    PLOOP(pliter,pl) ucons[pl]=Uievolve[pl];


    //////////////////
    //
    // Get Ufcent(pf[cent])
    //
    //////////////////

    failreturn=get_state(pf,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pf,&q,ptrgeom,Ufcent);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");


    //////////////////
    //
    // Get Uicent
    //
    //////////////////
    
    // UEVOLVE -> UDIAG
    // Assumes that Uievolve is like unewglobal and is at general U[] position (i.e. U[B1..B3] staggered and otherwise centered for FLUXB==FLUXCTSTAG)
    UtoU(UEVOLVE,UDIAG,ptrgeom,Uievolve,Uicent);

    // Override B1..B3 with correct centered versions (correct both for value and geometry)
    // ensure B^i is really at center even if FLUXB==FLUXCTSTAG (that would have Ui[B1..B3] at staggered)
    // Assumes, as very generally true, that U[B1..B3] never change and can never be adjusted.
    PLOOPBONLY(pl) Uicent[pl]=Ufcent[pl];


    ////////////////
    //
    // Get deltaUavg[] and also modify ucons if required and should
    //
    ////////////////
    diag_fixup_dUandaccount(Uicent, Ufcent, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);



  }// end if finalstep>0



  return(0);
}




// account for changes by tracking conserved quantities
// accounts for both failures and floor recoveries
// only called on final step of RK once unew is defined since only on final step is unew modified if floor encountered
// ONLY used by phys.ffde.c inversion routine when E^2>B^2
// Assume Ui and Uf in UDIAG form
int diag_fixup_U(int docorrectucons, FTYPE *Ui, FTYPE *Uf, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep,int whocalled)
{
  FTYPE Uicent[NPR],Ufcent[NPR];
  struct of_state q;
  int failreturn;
  int pliter,pl,enerregion, tscale;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int docorrectuconslocal;


  // count whocalled diag_fixup()
  count_whocalled(ptrgeom, finalstep, whocalled);
  


  if(finalstep>0){ // only account if on full timestep (assume only called if finalstep==1


    ///////////
    // determine if within correctable region
    ///////////
    docorrectuconslocal=diag_fixup_correctablecheck(docorrectucons,ptrgeom);


    ///////////
    //
    // First get correction (don't do inside enerregion loop since would be overly expensive and assume will need correction for at least one enerregion)
    //
    ///////////

    ///////////
    //
    // Change ucons (GODMARK: assumes Uf at CENT since uses single ptrgeom -- even though ucons is normally capable of being staggered for B1..B3)
    //
    ///////////
    if(DOENOFLUX != NOENOFLUX){ // JONMARK
      // notice that geometry comes after subtractions/additions of EOMs
      UtoU(UDIAG,UEVOLVE,ptrgeom,Uf,ucons); // convert from UDIAG->UEVOLVE 
    }


    // get Uicent and Ufcent
    // Assumes that Ui is like unewglobal and is at general U[] position (i.e. U[B1..B3] staggered and otherwise centered for FLUXB==FLUXCTSTAG)
    PLOOP(pliter,pl){
      Uicent[pl]=Ui[pl];
      Ufcent[pl]=Uf[pl];
    }

    // ensure B^i is really at center even if FLUXB==FLUXCTSTAG (that would have Ui[B1..B3] at staggered)
    // Assumes, as very generally true, that U[B1..B3] never change and can never be adjusted.
    PLOOPBONLY(pl) Uicent[pl]=Ufcent[pl];


    // Get deltaUavg[] and also modify ucons if required and should
    diag_fixup_dUandaccount(Uicent, Ufcent, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);


  }



  return(0);
}



// 0 = primitive (adds rho,u in comoving frame)
// 1 = conserved but rho,u added in ZAMO frame
// 2 = conserved but ignore strict rho,u change for ZAMO frame and instead conserved momentum (doesn't keep desired u/rho, b^2/rho, or b^2/u and so that itself can cause problems
#define FIXUPTYPE 1

// finalstep==0 is non-accounting, finalstep==1 is accounting
int fixup1zone(FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{
  int pliter,pl;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE ftempA,ftempB;
  struct of_state q;
  struct of_state dq;
  FTYPE prfloor[NPR];
  FTYPE prdiag[NPR];
  FTYPE pr0[NPR];
  FTYPE prnew[NPR];
  FTYPE U[NPR];
  int checkfl[NPR];
  int failreturn;
  int didchangeprim;
  FTYPE scalemin[NPR];
  //  FTYPE ucovzamo[NDIM];
  //  FTYPE uconzamo[NDIM];
  FTYPE dpr[NPR];
  FTYPE dU[NPR];
  //  FTYPE P,Pnew;
  int jj;
  int badinversion;



  
  // assign general floor variables
  // whether to check floor condition
  PALLLOOP(pl){
    checkfl[pl]=0;
    pr0[pl]=pr[pl];
    prdiag[pl]=pr0[pl];
  }

  // shouldn't fail since before and after states should be ok, as
  // long as would have changed the value.  Check will occur if
  // simulation continues ok.  Could place check inside if below.
  didchangeprim=0;



  ////////////
  //
  // Set which quantities to check
  //
  ////////////
  if(DOEVOLVERHO){
    checkfl[RHO]=1;
  }
  if(DOEVOLVEUU){
    checkfl[UU]=1;
  }


    
  ////////////
  //
  // Only apply floor if cold or hot GRMHD
  //
  ////////////
  if(DOEVOLVERHO||DOEVOLVEUU){


    //////////////
    //
    // get floor value
    //
    //////////////
    set_density_floors(ptrgeom,pr,prfloor);
    scalemin[RHO]=RHOMINLIMIT;
    scalemin[UU]=UUMINLIMIT;
    
    

    //////////////
    //
    // Set super low floor
    //
    //////////////
    PALLLOOP(pl){
      if(checkfl[pl]){
        if(prfloor[pl]<scalemin[pl]) prfloor[pl]=scalemin[pl];
      }
    }
    


    /////////////////////////////
    //
    // Get new primitive if went beyond floow
    //
    /////////////////////////////

    PALLLOOP(pl){
      if ( checkfl[pl]&&(prfloor[pl] > pr[pl]) ){
        didchangeprim=1;
        //dualfprintf(fail_file,"%d : %d %d %d : %d : %d : %21.15g - %21.15g\n",pl,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,checkfl[pl],prfloor[pl],pr[pl]); 
        // only add on full step since middle step is not really updating primitive variables
        prnew[pl]=prfloor[pl];
      }
      else prnew[pl]=pr[pl];
    }




    //////////////////////////
    //
    // ONLY do something if want lower than floor
    //
    //////////////////////////
    if(didchangeprim){

#if(FIXUPTYPE==0)
      // effectively adds mass/internal energy in comoving frame, which can lead to instabilities as momentum is added
      // For example, occurs on poles where u^r\sim 0 (stagnation surface) which launches artificially high u^t stuff only because goes below floor for a range of radii and so adds momentum to low density material
      // 
      PALLLOOP(pl){
        pr[pl]=prnew[pl];
      }


#elif(FIXUPTYPE==1 || FIXUPTYPE==2)
      // mass and internal energy added in frame not necessarily the comoving frame
      // using a frame not directly associated with comoving frame avoids arbitrary energy-momentum  growth
      // GODMARK: FIXUPTYPE==1 doesn't exactly match between when b^2/\rho_0>BSQORHOULIMIT such that amount of mass added will force equality of b^2/\rho_0==BSQORHOLIMIT, so this may lead to problems.
      // physically FIXUPTYPE==1 models some non-local transport of baryons and energy to a location that supposedly occurs when b^2/rho_0 is too large.
      // FIXUPTYPE==2 models a local injection of baryons/energy with momentum conserved and energy-momentum conserved if mass injected.  Injection will slow flow.  Essentially there is an ad hoc conversion of kinetic/thermal energy into mass energy.

      // compute original conserved quantities
      failreturn=get_state(pr,ptrgeom,&q);
      if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
      failreturn=primtoU(UNOTHING,pr,&q,ptrgeom,U);
      if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");

      // get change in primitive quantities
      PALLLOOP(pl) dpr[pl]=0.0; // default
      // use ZAMO velocity as velocity of inserted fluid
      for(pl=RHO;pl<=UU;pl++) dpr[pl]=prnew[pl]-pr[pl];
      set_zamo_velocity(WHICHVEL,ptrgeom,dpr);

      // get change in conserved quantities
      failreturn=get_state(dpr,ptrgeom,&dq);
      failreturn=primtoU(UNOTHING,dpr,&dq,ptrgeom,dU);
      if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");


      if(FIXUPTYPE==1){
        // then done, dU is right
      }
      else if(FIXUPTYPE==2){
        // then don't allow momentum to change regardless of meaning for implied rho,u
        dU[U1]=dU[U2]=dU[U3]=0.0;

        pl=UU;
        if ( checkfl[pl]&&(prfloor[pl] > pr[pl]) ){
          // then must change dU[UU]
        }
        else dU[UU]=0.0; // if only mass added, then no change needed to energy-momentum


      }

      // get final new conserved quantity
      PALLLOOP(pl) U[pl]+=dU[pl];



      // pr finally changes here
      // get primitive associated with new conserved quantities
      struct of_newtonstats newtonstats;
      failreturn=Utoprimgen(finalstep,OTHERUTOPRIM,UNOTHING,U,ptrgeom,pr,&newtonstats);
      badinversion = (failreturn>=1 || IFUTOPRIMFAIL(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)));

      if(badinversion){
        if(debugfail>=2) dualfprintf(fail_file,"Utoprimgen failed in fixup.c");
        // if problem with Utoprim, then just modify primitive quantities as normal without any special constraints
        PALLLOOP(pl){
          pr[pl]=prnew[pl];
        }
      }
#endif


    }// end if didchangeprim

  }// end if cold or hot GRMHD




  ///////////
  //
  // since inflow check is on boundary values, no need for inflow check here
  //
  ///////////


  ///////////////////////////////
  //
  // account for primitive changes
  //
  ///////////////////////////////
  if(didchangeprim&&FLOORDIAGS){// FLOORDIAGS includes fail diags
    int docorrectucons=1;
    diag_fixup(docorrectucons,prdiag, pr, ucons, ptrgeom, finalstep,COUNTFLOORACT);
    // now prdiag=pr as far as diag_fixup() is concerned, so next changes are new changes (i.e. don't cumulative w.r.t. pr0 multiple times, since that (relative to pr each time) would add each prior change to each next change)
    PALLLOOP(pl) prdiag[pl]=pr[pl];
  }



  ////////////////////
  //
  // limit gamma wrt normal observer
  //
  ////////////////////
#if(WHICHVEL==VELREL4)
  int docorrectucons=1;
  didchangeprim=0;

  failreturn=limit_gamma(GAMMAMAX,pr,ucons,ptrgeom,-1);
  if(failreturn>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
  if(failreturn==-1) didchangeprim=1;

  if(didchangeprim&&FLOORDIAGS){// FLOORDIAGS includes fail diags
    diag_fixup(docorrectucons,prdiag, pr, ucons, ptrgeom, finalstep,COUNTLIMITGAMMAACT);
    PALLLOOP(pl) prdiag[pl]=pr[pl];
  }

#endif// end if WHICHVEL==VEL4REL





  //////////////////////////
  // now keep track of modified primitives via conserved quantities

  //  if(didchangeprim){
  // assume once we go below floor, all hell will break loose unless we calm the storm by shutting down this zone's relative velocity
  // normal observer velocity
  // i.e. consider this a failure
  //GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)= 1;
  //  }

  


  return(0);
}


// number of vote checks per zone
#define MAXVOTES 8

#define NUMCHECKS 2
#define ISGAMMACHECK 0
#define ISUUCHECK 0

// GODMARK: function is 2D right now, but works in 3D, it just uses only x-y plane for checking

// check whether solution seems reasonable
// useful if b^2/rho\gg 1 or try to approach stationary model where variations in space shouldn't be large zone to zone
// checks whether u or gamma is much different than surrounding zones.
// if true, then flag as failure, else reasonable solution and keep
// can't assume failed zones are reasonably set
// fixup_checksolution() currently only uses pflag[FLAGUTOPRIMFAIL]
int fixup_checksolution(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],int finalstep)
{
  //  int inboundloop[NDIM];
  //  int outboundloop[NDIM];
  //  int innormalloop[NDIM];
  //  int outnormalloop[NDIM];
  //  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  //  int riin,riout,rjin,rjout,rkin,rkout;
  //  int dosetbc[COMPDIM*2];
  //  int ri;
  //  int boundvartype=BOUNDINTTYPE;
  // extra memory
  FTYPE (*gammacheck)[NSTORE2][NSTORE3][NPR];
  extern void get_advance_startendindices(int *is,int *ie,int *js,int *je,int *ks,int *ke);
  int is,ie,js,je,ks,ke;

  ///////////////////////////
  //
  // setup memory space for gammacheck
  //
  //////////////////////////

  // use UU space for holding \gamma
  gammacheck = GLOBALPOINT(ptemparray); // should be ok to use since fixup_checksolution() not called inside higher order or EOS or in fixup_utoprim() where also used


  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  //  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);
  // +-1 extra so can do check
  is=Uconsevolveloop[FIS]-SHIFT1;
  ie=Uconsevolveloop[FIE]+SHIFT1;
  js=Uconsevolveloop[FJS]-SHIFT2;
  je=Uconsevolveloop[FJE]+SHIFT2;
  ks=Uconsevolveloop[FKS]-SHIFT3;
  ke=Uconsevolveloop[FKE]+SHIFT3;




  ///////////////////////////////////
  //
  // determine gamma to be used
  //
  ///////////////////////////////////
  //  LOOPXalldir{
  //////  COMPZSLOOP(is,ie,js,je,ks,ke){
#pragma omp parallel  // don't need full (i.e. don't need EOS here)
  {
    int i,j,k;
    FTYPE ucon[NDIM],others[NUMOTHERSTATERESULTS];
    FTYPE qsq;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //    if(1|| (GLOBALMACP0A1(pflag,i,j,k,FLAGBSQORHO)||GLOBALMACP0A1(pflag,i,j,k,FLAGBSQOU))&&(IFUTOPRIMFAILORFIXED(GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)))){
      if(1){
        get_geometry(i,j,k,CENT,ptrgeom);
      
#if(WHICHVEL==VELREL4)
        MYFUN(gamma_calc(MAC(pv,i,j,k),ptrgeom,&MACP0A1(gammacheck,i,j,k,UU),&qsq),"fixup_checksolution: gamma calc failed\n","fixup.c",1);
#else
        if (ucon_calc(MAC(pv,i,j,k), ptrgeom, ucon, others) >= 1)  FAILSTATEMENT("fixup.c:fixup_checksolution()", "ucon_calc()", 1);
        MACP0A1(gammacheck,i,j,k,UU)=ucon[TT];
#endif
      }// end if 1
    }// end 3D LOOP
  }// end parallel region






  ///////////////////////////////////
  //
  // see if percent differences (actually factors) are too large
  //
  ///////////////////////////////////

  ////////  COMPZLOOP {
#pragma omp parallel  // don't need copyin, doesn't currently use global non-array vars
  {
    int i,j,k;
    int l;
    FTYPE percdiff[NUMCHECKS][MAXVOTES];
    int vote[NUMCHECKS];
    int numvotes[NUMCHECKS];
    int checkcondition[NUMCHECKS];
    int checki;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPZLOOP;


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);



      // doesn't need bound of pflag since don't check pflag of surrounding values (assumes if comparing with failure point, failure point is reasonable afterwards if before)
      //    if(1 || (GLOBALMACP0A1(pflag,i,j,k,FLAGBSQORHO)||GLOBALMACP0A1(pflag,i,j,k,FLAGBSQOU))&&(IFUTOPRIMFAILORFIXED(GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)))){// if b^2/{rho,u}\gg 1 and not failure already, check if solution is reasonable

      checkcondition[ISGAMMACHECK]=(MACP0A1(gammacheck,i,j,k,UU)>=2.0);
      checkcondition[ISUUCHECK]=1;

      // use fabs in case gamma<0 or especially if u<zerouuperbaryon*prim[RHO] that can easily happen
      if(checkcondition[ISGAMMACHECK]){
        percdiff[ISGAMMACHECK][0]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,i,jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,i,jp1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][1]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,i,jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,i,jm1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][2]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),j,k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,ip1mac(i),j,k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][3]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,im1mac(i),j,k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,im1mac(i),j,k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;

        percdiff[ISGAMMACHECK][4]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,ip1mac(i),jp1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][5]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,ip1mac(i),jm1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][6]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,ip1mac(i),jp1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
        percdiff[ISGAMMACHECK][7]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,im1mac(i),jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(gammacheck,im1mac(i),jm1mac(j),k,UU)/MACP0A1(gammacheck,i,j,k,UU)) : -1;
      }
      
      if(checkcondition[ISUUCHECK]){
        percdiff[ISUUCHECK][0]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,i,jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,i,jp1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][1]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,i,jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,i,jm1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][2]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),j,k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,ip1mac(i),j,k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][3]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,im1mac(i),j,k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,im1mac(i),j,k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
      
        percdiff[ISUUCHECK][4]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,ip1mac(i),jp1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][5]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,ip1mac(i),jm1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][6]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,ip1mac(i),jp1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,ip1mac(i),jp1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
        percdiff[ISUUCHECK][7]=(IFUTOPRIMNOFAILORFIXED(GLOBALMACP0A1(pflag,im1mac(i),jm1mac(j),k,FLAGUTOPRIMFAIL))) ? fabs(MACP0A1(pv,im1mac(i),jm1mac(j),k,UU)/MACP0A1(pv,i,j,k,UU)) : -1;
      }

      //////////////////////////
      //
      // determine how many votes
      //
      //////////////////////////
      for(checki=0;checki<NUMCHECKS;checki++){
        vote[checki]=0;
        numvotes[checki]=0;
      }
      for(l=0;l<MAXVOTES;l++){
        // No vote for failed zones
        if(checkcondition[ISGAMMACHECK] && percdiff[ISGAMMACHECK][l]>=0.0){
          if( (fabs(percdiff[ISGAMMACHECK][l])>GAMMAPERCDIFFMAX)||(fabs(percdiff[ISGAMMACHECK][l])<1.0/GAMMAPERCDIFFMAX) ) vote[ISGAMMACHECK]++;
          numvotes[ISGAMMACHECK]++;
        }
        if(checkcondition[ISUUCHECK] && percdiff[ISUUCHECK][l]>=0.0){
          if((DOEVOLVEUU)&& ((fabs(percdiff[ISUUCHECK][l])>UPERCDIFFMAX)||(fabs(percdiff[ISUUCHECK][l])<1.0/UPERCDIFFMAX)) ) vote[ISUUCHECK]++;
          numvotes[ISUUCHECK]++;
        }
      }

      /////////////////////
      //    
      // Use majority rule:  This allows shock fronts, but no faster members.
      // Don't include failures in voting process (either count or total voting)
      // > is used in case degenerate condition 0>0 so if no votes and no total, then do nothing
      //
      /////////////////////
      checki=ISGAMMACHECK;
      if(checkcondition[checki] && (vote[checki]>numvotes[checki]*0.5)){
        // then majority rules
        //      stderrfprintf("caught one-0: %d %d\n",i,j);
        GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILGAMMAPERC;
      }
      checki=ISUUCHECK;
      if(checkcondition[checki] && (vote[checki]>numvotes[checki]*0.5)){
        // then majority rules
        //      stderrfprintf("caught one-1: %d %d\n",i,j);
        GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILUPERC;
      }

    }// end COMPZLOOP
  }// end parallel region

  return(0);
}













// GODMARK: fixup_utoprim() is 2D function that works in 3D, but doesn't consider 3rd direction (i.e. it doesn't average in 3rd direction)


#if(MPIEQUALNONMPI==1)
#define ORDERINDEPENDENT 1 // no choice
// whether fixup_utoprim should be order-independent or not
#else
#define ORDERINDEPENDENT 1 // NO choice anymore, with fixup done before assigning initial D0, gamma, etc.
#endif

// whether to do specific chosen points for averages or do maximal average
#define GENERALAVERAGE 1

// whether to conserve D when averaging (uses original D to keep D constant -- less problematic compared to how used in limit_gamma() where original \gamma might be quite large)
// This is probably not useful
// More useful to use conserved D from unew as reference D0 so particle mass really conserved -- and this is ok compared to limit_gamma since newly averaged 4-velocity will not imply a large \gamma (unless averaging nearly failed regions!)
// So disable for now
#define DO_CONSERVE_D_INFAILFIXUPS 0 

#define HANDLEUNEG 0
// seems to keep failing with this, so probably treating u<zerouuperbaryon*prim[RHO] like failure is better idea
// 0: treat as full failure
// 1: treat as floor issue

#define HANDLERHONEG 0
// seems to keep failing with this, so probably treating rho<0 like failure is better idea
// 0: treat as full failure
// 1: treat as floor issue

#define HANDLERHOUNEG 0
// seems to keep failing with this, so probably treating rho<0 like failure is better idea
// 0: treat as full failure
// 1: treat as floor issue


// fix the bad solution as determined by utoprim() and fixup_checksolution()
// needs fail flag over -1..N, but uses p at 0..N-1
int fixup_utoprim(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR];
  extern void get_advance_startendindices(int *is,int *ie,int *js,int *je,int *ks,int *ke);
  int is,ie,js,je,ks,ke;




  // this average only works if using 4 velocity since only then guaranteed solution is good after interpolation
  if(WHICHVEL==VEL3) return(0); // just stick with static, best can do
  if(EOMTYPE==EOMFFDE) return(0); // nothing to do




  if(ORDERINDEPENDENT){
    // then on the right-side of equations should appear ptoavg and on left side pv
    // put another way, ptoavg contains constant thing never to be modified.  Usually then, pv is never used except to have something assigned to it, but could be conditions and such post-assignment that modifies the assignment but NOT the behavior of the routine otherwise.
    ptoavg=GLOBALPOINT(ptemparray);
    // ptemparray is used since in step_ch.c not needed after each advance()
    // ptemparray is template for values used to fixup failures.
    // In fixup-loop below, only change/fixup pv, not prc.  Otherwise order-dependent and not MPI friendly
    // ptoavg is only used for surrounding zones used to fixup local zone
    //LOOPXalldir

    get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);
    
    // +-1 extra so can do check
    is += -SHIFT1;
    ie += +SHIFT1;
    js += -SHIFT2;
    je += +SHIFT2;
    ks += -SHIFT3;
    ke += +SHIFT3;

    //////    COMPZSLOOP(is,ie,js,je,ks,ke) PALLLOOP(pl) MACP0A1(ptoavg,i,j,k,pl)=MACP0A1(pv,i,j,k,pl);
    copy_3dnpr(is,ie,js,je,ks,ke,pv,ptoavg); // GODMARK: Won't PLOOP be PALLLOOP here, so ok?
  }
  else ptoavg=pv;

  ///////////////////////////////////////////
  //
  // Save original pflag before changed.  This ensures thread-safety
  //
  ///////////////////////////////////////////

  // shouldn't modify pflagfailorig once inside parallel region below
  // Should only modify final (true) pflag value
  copy_3dpftype_special_fullloop(GLOBALPOINT(pflag),GLOBALPOINT(pflagfailorig));





  ///////////////////////////////////
  //
  // first check for bad solutions and try to fix based upon average surrounding solution
  //
  //////////////////////////////////

  /////////  COMPZLOOP
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // accounting requires state stuff
  {
    int i,j,k,pl,pliter;
    FTYPE gamma,alpha,vsq,ucon[NDIM],others[NUMOTHERSTATERESULTS];
    FTYPE qsq;
    FTYPE pr0[NPR];
    PFTYPE mypflag;
    int fixed;
    int startpl,endpl;
    FTYPE ftemp;
    FTYPE D0;
    //
    int limitedgamma;
    int nogood;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;


    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPZLOOP;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      // get failure flag
      mypflag=GLOBALMACP0A0(pflagfailorig,i,j,k);



      if(IFUTOPRIMFAILFIXED(mypflag)){
        /////////////////
        //
        // see if utoprim() previously fixed so can do diagnostics
        // only do accounting
        //
        //////////////////

        // set pre-fixed primitives
        PALLLOOP(pl)    pr0[pl]=MACP0A1(ptoavg,i,j,k,pl);
        get_geometry(i,j,k,CENT,ptrgeom);

        /////////////////////////////////
        //
        // ACCOUNTING (static or average)
        //
        /////////////////////////////////
        fixuputoprim_accounting(i, j, k, mypflag, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);


      }
      else if( IFUTOPRIMFAIL(mypflag)){
        /////////////////
        //
        // see if utoprim() failed
        //
        //////////////////

        fixed=0; // assume not fixed yet

        // set pre-fixed primitives
        // put back inside "if" when not superdebugging since wasteful of cpu
        PALLLOOP(pl)    pr0[pl]=MACP0A1(ptoavg,i,j,k,pl);
        get_geometry(i,j,k,CENT,ptrgeom);

        
        /////////////////////////////
        //
        // Check if want to average
        //
        // only modified if doing some kind of averaging, else static since already kept old value in inversion method
        //
        //////////////////////////////
        if(UTOPRIMADJUST==UTOPRIMAVG){

          /////////////////
          //
          // choose which range of quantities to average
          //
          //////////////////
          // field is evolved fine, so only average non-field
          if(mypflag==UTOPRIMFAILRHOUNEG && HANDLERHOUNEG==1){
            startpl=RHO;
            endpl=UU;
          }
          else if(mypflag==UTOPRIMFAILU2AVG1 || mypflag==UTOPRIMFAILU2AVG2 || mypflag==UTOPRIMFAILU2AVG1FROMCOLD || mypflag==UTOPRIMFAILU2AVG2FROMCOLD || mypflag==UTOPRIMFAILUPERC || mypflag==UTOPRIMFAILUNEG && (HANDLEUNEG==1) ){
            startpl=UU;
            endpl=UU;
          }
          else if(mypflag==UTOPRIMFAILRHONEG && HANDLERHONEG==1){
            startpl=RHO;
            endpl=RHO;
          }
          else{
            // then presume inversion failure with no solution or assuming rho<=0 or u<=zerouuperbaryon*prim[RHO] is bad inversion if HANDLE?NEG==0
            startpl=RHO;
            endpl=U3;
          }





          //////////////////////////////
          //
          // fixup negative densities
          //
          //////////////////////////////
          fixup_negdensities(&fixed, startpl, endpl, i, j, k, mypflag, pv,ptoavg, ptrgeom, pr0, ucons, finalstep);

          if(fixed==1 && (startpl<=RHO && endpl>=U1)){
            // then fixup but only changed densities, so still need to process non-densities
            startpl=U1; // start at U1 (first velocity) and finish at same ending if was ending on some velocity
            fixed=0; // then reset fixed->0 so can still modify these remaining quantities -- otherwise v^i would be unchanged even if wanted to average that out for (e.g.) negative density results for inversions.
          }


          //////////////////////////////
          //
          // other kinds of failures not caught by above (inversion convergence type failures)
          //
          //////////////////////////////
          if(fixed==0){


            /////////////////////
            //
            // fix using average of surrounding good values, if they exist
            //
            //////////////////////
            nogood=0;
#if(GENERALAVERAGE==1)
            nogood=general_average(startpl, endpl, i, j, k, mypflag, GLOBALPOINT(pflagfailorig) ,pv,ptoavg,ptrgeom);
#else
            nogood=simple_average(startpl,endpl,i,j,k, GLOBALPOINT(pflagfailorig) ,pv,ptoavg);
#endif

            /////////////////////
            //
            // If no good surrounding value found, then average bad values in some way
            //
            //////////////////////
            if(nogood){
              fixup_nogood(startpl, endpl, i, j, k, pv,ptoavg, pbackup, ptrgeom);
            }


            /////////////////////
            //
            // Things to do only if modifying density
            //
            //////////////////////
            if(startpl<=RHO && endpl>=RHO){


#if(DO_CONSERVE_D_INFAILFIXUPS)

              if(mypflag==UTOPRIMFAILGAMMAPERC || 1){ // GODMARK: always doing it
                // Use D0 to constrain how changing u^t changes rho
                // GODMARK: Why not used evolved D=\rho_0 u^t  from conserved quantity?
                // GODMARK: See fixup.c's limit_gamma() notes on why using conserved version of D not good to use
                // Here we ignore all conserved quantities and just ensure that D0 is conserved (close) to original value after averaging that assumes original value was reasonable
                // This is probably not necessary or useful
                D0 = MACP0A1(ptoavg,i,j,k,RHO)*ucon[TT];

                ///////////////////////////////////////////
                //
                // constrain change in density so conserve particle number
                // always do it?
                //
                //////////////////////////////////////////
                if (ucon_calc(MAC(pv,i,j,k), ptrgeom, ucon,others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);
                MACP0A1(pv,i,j,k,RHO) = D0/ucon[TT];
              }
#endif

            }// end over density





            /////////////////////
            //
            // Things to do only if modifying velocity
            //
            //////////////////////
            if(startpl<=U1 && endpl>=U3){

#if(WHICHVEL==VELREL4)
              //////////////
              //
              // check gamma to so calibrate new gamma to no larger than previous gamma
              /////////////
              MYFUN(gamma_calc(MAC(ptoavg,i,j,k),ptrgeom,&gamma,&qsq),"fixup_utoprim: gamma calc failed\n","fixup.c",1);
              if (ucon_calc(MAC(ptoavg,i,j,k), ptrgeom, ucon, others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);
              //          alpha = 1. / sqrt(-ptrgeom->gcon[GIND(0,0)]);
              alpha = ptrgeom->alphalapse;
              vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);

              limitedgamma=0;
              if(gamma>GAMMAMAX){
                limitedgamma=1;
                if(debugfail>=2) dualfprintf(fail_file,"initial gamma: %21.15g,  max: %21.15g, initial vsq: %21.15g\n",gamma,GAMMAMAX,vsq);
                // limit:
                gamma=GAMMAMAX;
              }
#else
              limitedgamma=0;
              if (ucon_calc(MAC(ptoavg,i,j,k), ptrgeom, ucon,others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);
#endif
          

              /////////////
              //
              // check new gamma to make sure smaller than original (i.e. for pv, not original ptoavg)
              //
              /////////////
            
#if(WHICHVEL==VELREL4)
              if(limit_gamma(gamma,MAC(pv,i,j,k),MAC(ucons,i,j,k),ptrgeom,-1)>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 2);

              if(debugfail>=3){
                if(limitedgamma){
                  // check gamma
                  MYFUN(gamma_calc(MAC(pv,i,j,k),ptrgeom,&gamma,&qsq),"fixup_utoprim: gamma calc failed\n","fixup.c",2);
                  dualfprintf(fail_file,"final gamma: %21.15g\n",gamma);
                }
              }
#endif

            }// end if dealing with velocity




#if(0)
            // DEBUG problem of launch with pressureless stellar model collapse
            PALLLOOP(pl) dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d j=%d k=%d pl=%d pv=%21.15g ptoavg=%21.15g\n",nstep,steppart,i,j,k,pl,MACP0A1(pv,i,j,k,pl),MACP0A1(ptoavg,i,j,k,pl));
#endif



          } // end if fixed==0

        }// end if not keeping static
        // else kept static
        
        
        
        /////////////////////////////////
        //
        // ACCOUNTING (static or average)
        //
        /////////////////////////////////
        fixuputoprim_accounting(i, j, k, mypflag, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);
          
          
          
      }// end if failure
    }// end over COMPZLOOP loop
  }// end over parallel region

  return(0);
}









// fix the bad solution as determined by utoprim() and fixup_checksolution()
// needs fail flag over -1..N, but uses p at 0..N-1
int fixup_utoprim_nofixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR];

  // this average only works if using 4 velocity since only then guaranteed solution is good after interpolation
  if(WHICHVEL==VEL3) return(0); // just stick with static, best can do
  if(EOMTYPE==EOMFFDE) return(0); // nothing to do


  ///////////////////////////////////
  //
  // Just loop over and only report problems
  //
  //////////////////////////////////
  
  ptoavg=pv;


  /////////  COMPZLOOP
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // accounting requires state stuff
  {
    int i,j,k,pl,pliter;
    FTYPE pr0[NPR];
    PFTYPE mypflag;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;


    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPZLOOP;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      // get failure flag
      mypflag=GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL);


      if( IFUTOPRIMFAIL(mypflag)){
        
        PALLLOOP(pl)    pr0[pl]=MACP0A1(ptoavg,i,j,k,pl);
        get_geometry(i,j,k,CENT,ptrgeom);
        
        /////////////////////////////////
        //
        // ACCOUNTING (static or average)
        //
        /////////////////////////////////
        fixuputoprim_accounting(i, j, k, mypflag, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);
          
      }// end if failure
    }// end over COMPZLOOP loop
  }// end over parallel region

  return(0);
}









// fixup negative densities
static int fixup_negdensities(int *fixed, int startpl, int endpl, int i, int j, int k, PFTYPE mypflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE prguess[NPR];

  ////////////////////////////
  //
  // deal with negative internal energy as special case of failure
  //

  if(*fixed!=0){
    if(mypflag==UTOPRIMFAILUNEG){
          
      if(STEPOVERNEGU==NEGDENSITY_NEVERFIXUP){ if(HANDLEUNEG==1) *fixed=1; }
      else if((STEPOVERNEGU==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP && finalstep)){

        if(HANDLEUNEG==1){
          // set back to floor level
          set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess);
          // GODMARK -- maybe too agressive, maybe allow more negative?
                
          if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED){
            // then pv is previous timestep value and can use to make fix
            if(-MACP0A1(pv,i,j,k,UU)<prguess[UU]){ *fixed=1; MACP0A1(pv,i,j,k,UU)=prguess[UU];} // otherwise assume really so bad that failure
          }
          else{
            // just treat as floor for all failures since do not know what updated quantity is
            MACP0A1(pv,i,j,k,UU)=prguess[UU];
            *fixed=1;
          }
        }// end if handling u<zerouuperbaryon*prim[RHO] in special way
      }// end if not allowing negative u or if allowing but not yet final step
      else if((STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP)&&(!finalstep)){
        if(HANDLEUNEG==1) *fixed=1; // tells rest of routine to leave alone and say ok solution, but don't use it to fix convergence failures for other zones
      }
    }// end if u<zerouuperbaryon*prim[RHO]
  }// end if not fixed


  ////////////////////////////
  //
  // deal with negative density as special case of failure (as for internal energy above)
  //

  if(*fixed!=0){
    if(mypflag==UTOPRIMFAILRHONEG){
          
      if(STEPOVERNEGRHO==NEGDENSITY_NEVERFIXUP){ if(HANDLERHONEG) *fixed=1; }
      else if((STEPOVERNEGRHO==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGRHO==NEGDENSITY_FIXONFULLSTEP && finalstep)){

        if(HANDLERHONEG==1){
          // set back to floor level
          set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess);
          // GODMARK -- maybe too agressive, maybe allow more negative?
                
          if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED){
            // then pv is previous timestep value and can use to make fix
            if(-MACP0A1(pv,i,j,k,RHO)<prguess[RHO]){ *fixed=1; MACP0A1(pv,i,j,k,RHO)=prguess[RHO];} // otherwise assume really so bad that failure
          }
          else{
            // just treat as floor for all failures since do not know what updated quantity is
            MACP0A1(pv,i,j,k,RHO)=prguess[RHO];
            *fixed=1;
          }
        }// end if handling rho<0 in special way
      }// end if not allowing negative rho or if allowing but not yet final step
      else if((STEPOVERNEGRHO==NEGDENSITY_FIXONFULLSTEP)&&(!finalstep)){
        if(HANDLERHONEG) *fixed=1; // tells rest of routine to leave alone and say ok solution, but don't use it to fix convergence failures for other zones
      }
    }// end if rho<0
  }// end if not fixed



  ////////////////////////////
  //
  // deal with negative density and negative internal energy as special case of failure (as for internal energy above)
  //

  if(*fixed!=0){
    if(mypflag==UTOPRIMFAILRHOUNEG){

      if(STEPOVERNEGRHOU==NEGDENSITY_NEVERFIXUP){ if(HANDLERHOUNEG) *fixed=1; }
      // GODMARK: Why use STEPOVERNEGU and STEPOVERNEGRH instead of STEPOVERNEGRHOU below?
      else if( (STEPOVERNEGRHOU==NEGDENSITY_ALWAYSFIXUP)  ||(STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP && STEPOVERNEGRHO==NEGDENSITY_FIXONFULLSTEP && finalstep)){

        if(HANDLERHOUNEG==1){
          // set back to floor level
          set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess);
          // GODMARK -- maybe too agressive, maybe allow more negative?
                
          if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED){
            // then pv is previous timestep value and can use to make fix
            if(-MACP0A1(pv,i,j,k,UU)<prguess[UU]){ *fixed=1; MACP0A1(pv,i,j,k,UU)=prguess[UU];} // otherwise assume really so bad that failure
            if(-MACP0A1(pv,i,j,k,RHO)<prguess[RHO]){ *fixed=1; MACP0A1(pv,i,j,k,RHO)=prguess[RHO];} // otherwise assume really so bad that failure
          }
          else{
            // just treat as floor for all failures since do not know what updated quantity is
            MACP0A1(pv,i,j,k,UU)=prguess[UU];
            MACP0A1(pv,i,j,k,RHO)=prguess[RHO];
            *fixed=1;
          }
        }// end if handling rho<0 and u<zerouuperbaryon*prim[RHO] in special way
      }// end if not allowing negative rho or if allowing but not yet final step
      else if(STEPOVERNEGRHOU==NEGDENSITY_FIXONFULLSTEP &&(!finalstep)){
        if(HANDLERHOUNEG) *fixed=1; // tells rest of routine to leave alone and say ok solution, but don't use it to fix convergence failures for other zones
      }
    }// end if rho<0 and u<zerouuperbaryon*prim[RHO]
  }// end if not fixed
        
  return(0);
}




// DOCOUNTNEG???? only applies for STEPOVERNEG???==NEGDENSITY_NEVERFIXUP

// whether to count any substep u<zerouuperbaryon*prim[RHO] as failure in debug data
// 2: always counted
// 1: only counted on final substep
// 0: never counted
#define DOCOUNTNEGU 1

// whether to count any substep rho<0 as failure in debug data
#define DOCOUNTNEGRHO 1

// whether to count any substep rho<0 && u<zerouuperbaryon*prim[RHO] case as failure in debug data
#define DOCOUNTNEGRHOU 1


// whether to make conserved quantity consistent with the associated fixed primitive quantity.
// required at least at the finalstep.
// 0: don't do it
// 1: adjust only if finalstep==1
// 2: adjust for all substeps
#define ADJUSTCONSERVEDQUANTITY 0

// ACCOUNTING (under any circumstance, static or average)
static int fixuputoprim_accounting(int i, int j, int k, PFTYPE mypflag, PFTYPE (*lpflag)[NSTORE2][NSTORE3][NUMPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  PFTYPE utoprimfailtype;
  int doadjustcons;
  struct of_state q;
  FTYPE (*utoinvert)[NSTORE2][NSTORE3][NPR];
  int docorrectucons;
  int pliter,pl;


  // account for changes by tracking conserved quantities
  // note that if new utoprim solution was impossible, we are using here prior solution as new state, which means the diagnostics won't be acccurate.  There's no way to fix this unless the fluxes always give a well-defined new primitive variable.
  // specific value of the mypflag>0 communicates the type of failure
  if(IFUTOPRIMNOFAIL(mypflag)){
    utoprimfailtype=-1;
    docorrectucons=0;
  }
  else if(
          (mypflag==UTOPRIMFAILCONV)|| // only used by 5D method currently
          (mypflag==UTOPRIMFAILCONVGUESSUTSQ)|| // rest are only used by 1D/2D method currently
          (mypflag>=UTOPRIMFAILCONVRET)||
          (mypflag==UTOPRIMFAILCONVW)||
          (mypflag==UTOPRIMFAILCONVUTSQVERYBAD)||
          (mypflag==UTOPRIMFAILNANGUESS)||
          (mypflag==UTOPRIMFAILNANRESULT)||
          (mypflag==UTOPRIMFAILCONVBADINVERTCOMPARE)||
          (mypflag==UTOPRIMFAILCONVUTSQ)||
          (mypflag==UTOPRIMFAILFAKEVALUE) // fake value for avoiding MPI boundary call and so MPI boundary values for fixup
          ){
    utoprimfailtype=COUNTUTOPRIMFAILCONV;
    docorrectucons=1;
  }
  else if(mypflag==UTOPRIMFAILRHONEG){
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.
    if(STEPOVERNEGRHO==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGRHO==1){
        if(finalstep){
          utoprimfailtype=COUNTUTOPRIMFAILRHONEG;
          docorrectucons=1;
        }
        else{
          utoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGRHO==2){
        utoprimfailtype=COUNTUTOPRIMFAILRHONEG;
        docorrectucons=1;
      }
      else{
        utoprimfailtype=-1;
        docorrectucons=0;
      }
    }   
    else if((STEPOVERNEGRHO)&&(!finalstep)){
      utoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      utoprimfailtype=COUNTUTOPRIMFAILRHONEG;
      docorrectucons=1;
    }
  }
  else if(mypflag==UTOPRIMFAILUNEG || mypflag==UTOPRIMFAILU2AVG1|| mypflag==UTOPRIMFAILU2AVG2 || mypflag==UTOPRIMFAILU2AVG1FROMCOLD|| mypflag==UTOPRIMFAILU2AVG2FROMCOLD){ // GODMARK: maybe want separate accounting
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.

    // default counting:
    if(mypflag==UTOPRIMFAILU2AVG1FROMCOLD|| mypflag==UTOPRIMFAILU2AVG2FROMCOLD) utoprimfailtype=COUNTCOLD;
    else utoprimfailtype=COUNTUTOPRIMFAILUNEG;

    // now set whether to ucons correction or override counting
    if(STEPOVERNEGU==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGU==1){
        if(finalstep){
          docorrectucons=1;
        }
        else{
          utoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGU==2){
        docorrectucons=1;
      }
      else{
        utoprimfailtype=-1;
        docorrectucons=0;
      }
    }   
    else if((STEPOVERNEGU)&&(!finalstep)){
      utoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      docorrectucons=1;
    }
  }
  else if(mypflag==UTOPRIMFAILRHOUNEG){
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.
    if(STEPOVERNEGRHOU==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGRHOU==1){
        if(finalstep){
          utoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
          docorrectucons=1;
        }
        else{
          utoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGRHOU==2){
        utoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
        docorrectucons=1;
      }
      else{
        utoprimfailtype=-1;
        docorrectucons=0;
      }
    }   
    else if((STEPOVERNEGRHOU)&&(!finalstep)){
      utoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      utoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
      docorrectucons=1;
    }
  }
  else if(mypflag==UTOPRIMFAILGAMMAPERC){
    utoprimfailtype=COUNTGAMMAPERC;
    docorrectucons=1;
  }
  else if(mypflag==UTOPRIMFAILUPERC){
    utoprimfailtype=COUNTUPERC;
    docorrectucons=1;
  }
  else if(mypflag==UTOPRIMFAILFIXEDENTROPY){
    utoprimfailtype=COUNTENTROPY;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mypflag==UTOPRIMFAILFIXEDCOLD){
    utoprimfailtype=COUNTCOLD;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mypflag==UTOPRIMFAILFIXEDBOUND1){
    utoprimfailtype=COUNTBOUND1;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mypflag==UTOPRIMFAILFIXEDBOUND2){
    utoprimfailtype=COUNTBOUND2;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mypflag==UTOPRIMFAILFIXEDONESTEP){
    utoprimfailtype=COUNTONESTEP;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mypflag==UTOPRIMFAILFIXEDUTOPRIM){
    dualfprintf(fail_file,"prior pflag not cleared: nstep=%ld steppart=%d t=%21.15g i=%d j=%d k=%d \n",nstep,steppart,t,i,j,k);
    utoprimfailtype=-1;
    docorrectucons=0;
  }
  else{
    dualfprintf(fail_file,"No such pflag failure type: %d for t=%21.15g nstep=%ld steppart=%d i=%d j=%d k=%d\n",mypflag,t,nstep,steppart,i,j,k);
    myexit(1);
  }



  if(utoprimfailtype!=-1){
    // diagnostics
    FTYPE prdiag[NPR],pr[NPR];
    PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
    diag_fixup(docorrectucons,prdiag, MAC(pv,i,j,k), MAC(ucons,i,j,k), ptrgeom, finalstep,(int)utoprimfailtype);
    PLOOP(pliter,pl) prdiag[pl]=MACP0A1(pv,i,j,k,pl);

    ////////////////
    //
    // reset true pflag counter to "no" (fixed) failure
    //
    ////////////////
    MACP0A1(lpflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILFIXEDUTOPRIM;


    // now adjust uf or ucum so agrees [already done in diag_fixup() but only if final step.  This allows one to do it on any step in case want to control behavior of conserved quantities or primitive quantities used to create fluxes.
    if(ADJUSTCONSERVEDQUANTITY&&docorrectucons){

      if((ADJUSTCONSERVEDQUANTITY==1)&&(finalstep)) doadjustcons=1;
      else if(ADJUSTCONSERVEDQUANTITY==2) doadjustcons=1;
      else doadjustcons=0;

      if(doadjustcons){
        //utoinvert=ucons;
        //          if(finalstep){ // last call, so ucum is cooked and ready to eat!
        //          }
        //          else{ // otherwise still iterating on primitives
        //            utoinvert=ulast;
        //          }
            
        MYFUN(get_state(MAC(pv,i,j,k), ptrgeom, &q),"fixup.c:fixup_utoprim()", "get_state()", 1);
        MYFUN(primtoU(UEVOLVE,MAC(pv,i,j,k), &q, ptrgeom, MAC(ucons,i,j,k)),"fixup.c:fixup_utoprim()", "primtoU()", 1);
      }
    }


  }// end if utoprimfailtype!=-1



#if(DOSUPERDEBUG)
  // only used to keep track of normal non-failure points
  // some non-utoprim-failure points will be picked up by this.  I.E. those with limitgammaact/flooract/inflowact
  if(utoprimfailtype==-1) superdebug_utoprim(pr0,MAC(pv,i,j,k),ptrgeom,utoprimfailtype);
  // collect values for non-failed and failed zones
#endif


  return(0);
}






// sums need to be symmerized in order to preserve exact symmetry
#define AVG4_1(pr,i,j,k,pl) (0.25*((MACP0A1(pr,i,jp1mac(j),k,pl)+MACP0A1(pr,i,jm1mac(j),k,pl))+(MACP0A1(pr,im1mac(i),j,k,pl)+MACP0A1(pr,ip1mac(i),j,k,pl))))  // symmetrized sum
#define AVG4_2(pr,i,j,k,pl) (0.25*((MACP0A1(pr,ip1mac(i),jp1mac(j),k,pl)+MACP0A1(pr,ip1mac(i),jm1mac(j),k,pl))+(MACP0A1(pr,im1mac(i),jp1mac(j),k,pl)+MACP0A1(pr,im1mac(i),jm1mac(j),k,pl))))  //symmetrized sum

#define AVG2_1(pr,i,j,k,pl) (0.5*(MACP0A1(pr,i,jp1mac(j),k,pl)+MACP0A1(pr,i,jm1mac(j),k,pl)))  
#define AVG2_2(pr,i,j,k,pl) (0.5*(MACP0A1(pr,ip1mac(i),j,k,pl)+MACP0A1(pr,im1mac(i),j,k,pl)))
#define AVG2_3(pr,i,j,k,pl) (0.5*(MACP0A1(pr,ip1mac(i),jp1mac(j),k,pl)+MACP0A1(pr,im1mac(i),jm1mac(j),k,pl)))
#define AVG2_4(pr,i,j,k,pl) (0.5*(MACP0A1(pr,ip1mac(i),jm1mac(j),k,pl)+MACP0A1(pr,im1mac(i),jp1mac(j),k,pl)))

#define AVG4_1(pr,i,j,k,pl) (0.25*((MACP0A1(pr,i,jp1mac(j),k,pl)+MACP0A1(pr,i,jm1mac(j),k,pl))+(MACP0A1(pr,im1mac(i),j,k,pl)+MACP0A1(pr,ip1mac(i),j,k,pl))))  // symmetrized sum
#define AVG4_2(pr,i,j,k,pl) (0.25*((MACP0A1(pr,ip1mac(i),jp1mac(j),k,pl)+MACP0A1(pr,ip1mac(i),jm1mac(j),k,pl))+(MACP0A1(pr,im1mac(i),jp1mac(j),k,pl)+MACP0A1(pr,im1mac(i),jm1mac(j),k,pl))))  //symmetrized sum


// whether doing causal averaging for non-U2AVG failures
#define CAUSALAVG 0

#define CAUSALAVG_WHEN_U2AVG 0 // causal way -- uses same normal loops
#define SIMPLEAVG_WHEN_U2AVG 1 // normal average
#define MAXUPOSAVG_WHEN_U2AVG 2 // Sasha way -- take min of positive internal energies
#define CAUSAL_THENMIN_WHEN_U2AVG 3 // Jon way -- use min of ANY (pos or neg) answer between (1) causal and (2) min of all positive
              
//#define HOWTOAVG_WHEN_U2AVG SIMPLEAVG_WHEN_U2AVG // can artificially pump up internal energy as in caustic test
//#define HOWTOAVG_WHEN_U2AVG CAUSALAVG_WHEN_U2AVG
#define HOWTOAVG_WHEN_U2AVG MAXUPOSAVG_WHEN_U2AVG
//#define HOWTOAVG_WHEN_U2AVG CAUSAL_THENMIN_WHEN_U2AVG

// more general averaging procedure that tries to pick best from surrounding good values
static int general_average(int startpl, int endpl, int i, int j, int k, PFTYPE mypflag, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom)
{
  int doavgcausal,failavglooptype;
  int pl,pliter;
  int ii,jj,kk;
  int jjj;
  struct of_state state_pv;
  FTYPE cminmax[NDIM][NUMCS];
  FTYPE superfast[NDIM];
  int numavg;
  int numavg0,numavg1;
  FTYPE mysum[2][NPR];
  int numupairs,qq,thisnotfail,thatnotfail,rnx,rny,rnz;
  FTYPE ref;
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE lastmin[NPR],avganswer0[NPR],avganswer1[NPR];
  int ignorecourant=0;
  int factor;



  if(debugfail>=2) dualfprintf(fail_file,"uc2general: mypflag=%d startpl=%d endpl=%d :: i=%d j=%d k=%d\n",mypflag,startpl,endpl,i,j,k); // not too much


  if(mypflag==UTOPRIMFAILU2AVG1 || mypflag==UTOPRIMFAILU2AVG2 || mypflag==UTOPRIMFAILU2AVG1FROMCOLD || mypflag==UTOPRIMFAILU2AVG2FROMCOLD){
    if(HOWTOAVG_WHEN_U2AVG==CAUSALAVG_WHEN_U2AVG){ // causal type
      doavgcausal=1;
      failavglooptype=0;
    }
    else if(HOWTOAVG_WHEN_U2AVG==MAXUPOSAVG_WHEN_U2AVG){ // Sasha type
      doavgcausal=0; // doesn't hurt, it's just wasted process
      failavglooptype=1;
    }
    else if(HOWTOAVG_WHEN_U2AVG==SIMPLEAVG_WHEN_U2AVG){ // simple type
      doavgcausal=0;
      failavglooptype=0;
    }
    else if(HOWTOAVG_WHEN_U2AVG==CAUSAL_THENMIN_WHEN_U2AVG){ // Jon's mixed type
      doavgcausal=1;
      failavglooptype=2; // do both loops
    }
  }
  else if(CAUSALAVG){ // other types of failures besides U2AVG
    failavglooptype=0;
    doavgcausal=1;
  }
  else{
    doavgcausal=0;
    failavglooptype=0;
  }  


 
  if(doavgcausal){
    //////////
    // first get wave speed to check if superfast so need causal averaging
    MYFUN(get_state(MAC(ptoavg,i,j,k), ptrgeom, &state_pv),"fixup.c:fixup_utporim()", "get_state()", 1);
    if(N1NOT1){
      MYFUN(vchar(MAC(ptoavg,i,j,k), &state_pv, 1, ptrgeom, &cminmax[1][CMAX], &cminmax[1][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar() dir=1or2", 1);
    }
    if(N2NOT1){
      MYFUN(vchar(MAC(ptoavg,i,j,k), &state_pv, 2, ptrgeom, &cminmax[2][CMAX], &cminmax[2][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar() dir=1or2", 2);
    }
    if(N3NOT1){
      MYFUN(vchar(MAC(ptoavg,i,j,k), &state_pv, 3, ptrgeom, &cminmax[3][CMAX], &cminmax[3][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar() dir=1or2", 3);
    }
    SLOOPA(jjj){
      superfast[jjj]=0; // assume subfast

#if(1)
      if(cminmax[jjj][CMIN]>0.0){ superfast[jjj]=1.0;} // superfast to right
      else if(cminmax[jjj][CMAX]<0.0){ superfast[jjj]=-1.0;} // superfast to left
#else
      if(MACP0A1(ptoavg,i,j,k,U1-1+jjj)>0.0){ superfast[jjj]=1.0;} // wind blows to right
      else if(MACP0A1(ptoavg,i,j,k,U1-1+jjj)<0.0){ superfast[jjj]=-1.0;} // wind blows to left
#endif
    }

  } // end if doing causal loop -- end getting wave speeds


          


        
  ///////////////////////////////////////////////////////////////
  //
  // average all surrounding good values (keeps symmetry)
  //
  /////////////////////////////////////////////////////////////
  numavg=0;
  numavg0=numavg1=0;
  for(pl=startpl;pl<=endpl;pl++){
    mysum[0][pl]=0.0;
    mysum[1][pl]=0.0;
    lastmin[pl]=1E50;
  }
  // size of little averaging box
  rnx=(SHIFT1==1) ? 3*SHIFT1 : 1;
  rny=(SHIFT2==1) ? 3*SHIFT2 : 1;
  rnz=(SHIFT3==1) ? 3*SHIFT3 : 1;
            
  // number of unique pairs
  factor=SHIFT1+SHIFT2+SHIFT3;
  numupairs=(int)(pow(3,factor)-1)/2;
            
  for(qq=0;qq<numupairs;qq++){
              
    // 1-d to 3D index
    ii=(int)(qq%rnx)-SHIFT1;
    jj=(int)((qq%(rnx*rny))/rnx)-SHIFT2;
    kk=(int)(qq/(rnx*rny))-SHIFT3;
              
    thisnotfail=IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i+ii,j+jj,k+kk));
    thatnotfail=IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i-ii,j-jj,k-kk));

    if(doavgcausal){
#if(1)
      if( ii<0 && superfast[1]>0 && thisnotfail ) thatnotfail=0; // then don't need this one
      if( ii>0 && superfast[1]<0 && thatnotfail ) thisnotfail=0; // then don't need this one
      if( jj<0 && superfast[2]>0 && thisnotfail ) thatnotfail=0; // then don't need this one
      if( jj>0 && superfast[2]<0 && thatnotfail ) thisnotfail=0; // then don't need this one
      if( kk<0 && superfast[3]>0 && thisnotfail ) thatnotfail=0; // then don't need this one
      if( kk>0 && superfast[3]<0 && thatnotfail ) thisnotfail=0; // then don't need this one
#elif(0)
      if(i<totalsize[1]/2 && ii<0) thatnotfail=0;
      if(i<totalsize[1]/2 && ii>0) thisnotfail=0;
      if(i>=totalsize[1]/2 && ii<0) thisnotfail=0;
      if(i>=totalsize[1]/2 && ii>0) thatnotfail=0;
#endif

    }

#if(1)
    // bigger than myself
    ref=MACP0A1(ptoavg,i,j,k,UU);
#elif(0)
    // bigger than 0
    ref=0.0;
#endif

    //      if(failavglooptype==1  || failavglooptype==2){  // only use if positive

    // number of quantities one summed
    numavg0+=thisnotfail+thatnotfail;
            
    if(failavglooptype==0 || failavglooptype==2){
      if(debugfail>=3) dualfprintf(fail_file,"uc2: i=%d i+ii=%d j=%d j+jj=%d k=%d k+kk=%d (e.g.) pl=%d pv=%21.15g\n",i,i+ii,j,j+jj,k,k+kk,RHO,MACP0A1(pv,i,j,k,RHO)); // still bit much
      for(pl=startpl;pl<=endpl;pl++){
        mysum[qq%2][pl]+=MACP0A1(ptoavg,i+ii,j+jj,k+kk,pl)*thisnotfail + MACP0A1(ptoavg,i-ii,j-jj,k-kk,pl)*thatnotfail;
        if(debugfail>=4) dualfprintf(fail_file,"uc2: i=%d i+ii=%d j=%d j+jj=%d k=%d k+kk=%d pl=%d pv=%21.15g\n",pl,i,i+ii,j,j+jj,k,k+kk,MACP0A1(pv,i,j,k,pl)); // bit much, just do for all pl above

#if(0) // DEBUG
        // DEBUG problem of launch with pressureless stellar model collapse
        dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d j=%d k=%d pl=%d pv=%21.15g thisnotfail=%d ptoavg1=%21.15g thatnotfail=%d ptoavg2=%21.15g :: ii=%d jj=%d kk=%d numavg0=%d\n",nstep,steppart,i,j,k,pl,MACP0A1(pv,i,j,k,pl),thisnotfail,MACP0A1(ptoavg,i+ii,j+jj,k+kk,pl),thatnotfail,MACP0A1(ptoavg,i-ii,j-jj,k-kk,pl),ii,jj,kk,numavg0);
#endif
      }
    }
    if(failavglooptype==1 || failavglooptype==2){ // only for U2AVG
      for(pl=startpl;pl<=endpl;pl++){
#if(0)
        ftemp=MACP0A1(ptoavg,i+ii,j+jj,k+kk,pl);
        if(ftemp>=ref){
          lastmin[pl]=MIN(lastmin[pl],ftemp); // smallest positive number
          numavg1++;
        }
                
        ftemp=MACP0A1(ptoavg,i-ii,j-jj,k-kk,pl);
        if(ftemp>=ref){
          lastmin[pl]=MIN(lastmin[pl],ftemp);
          numavg1++;
        }
#else
        ftemp1=MACP0A1(ptoavg,i+ii,j+jj,k+kk,pl);
        ftemp2=MACP0A1(ptoavg,i-ii,j-jj,k-kk,pl);
        if(ftemp1>=ref && ftemp2>=ref){
          lastmin[pl]=MIN(MIN(lastmin[pl],ftemp1),ftemp2); // smallest positive number if both of pair are larger than my value
          numavg1++;
        }
#endif
      }
    }
  } //end loop over pairs

        



  ///////////////////////
  //
  // all loops over surrounding points is done, now get average answer
  //
  ////////////////////////
  if(failavglooptype==0 || failavglooptype==2){       
    if(numavg0!=0) for(pl=startpl;pl<=endpl;pl++){
        avganswer0[pl]=(mysum[0][pl]+mysum[1][pl])/((FTYPE)(numavg0));
      }
  }
  if(failavglooptype==1 || failavglooptype==2){
    for(pl=startpl;pl<=endpl;pl++){
      avganswer1[pl]=lastmin[pl];
    }
  }
          

  ///////////////
  //
  // choose final answer depending upon loop type
  //
  ///////////////
  if(failavglooptype==0 || ((failavglooptype==2)&&(numavg1==0)) ){  
    if(numavg0!=0) for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer0[pl];
    numavg=numavg0;
  }
  if(failavglooptype==1 || ((failavglooptype==2)&&(numavg0==0)) ){  
    //      if(numavg1!=0 && (MACP0A1(pv,i,j,k,pl)<avganswer1[pl]) ) for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl]; // else keep same as original answer
    //      if(numavg1!=0 && (MACP0A1(ptoavg,i,j,k,pl)<avganswer1[pl]) ) for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl]; // else keep same as original answer
    if(numavg1!=0) for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl]; // else keep same as original answer
    //      if(numavg1==2) for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl]; // else keep same as original answer
    numavg=numavg1;
  }
  if(failavglooptype==2 && (numavg0!=0) && (numavg1!=0) ){ // here if both numavg0!=0 and numavg1!=0
    //for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=MIN(MIN(avganswer1[pl],avganswer0[pl]),MACP0A1(pv,i,j,k,pl));
    if(MACP0A1(pv,i,j,k,pl)<avganswer1[pl] && MACP0A1(pv,i,j,k,pl)<avganswer0[pl]){
      for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=MIN(avganswer1[pl],avganswer0[pl]);
    }
    else if(MACP0A1(pv,i,j,k,pl)<avganswer1[pl] ){ // Sasha scheme
      for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl];
    }
    else if(MACP0A1(pv,i,j,k,pl)<avganswer0[pl] ){ // Causal scheme
      for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer0[pl];
    }
    //for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer0[pl];
    //for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=avganswer1[pl];
    numavg=MAX(numavg0,numavg1); // only matters now that this is nonzero
  }


  if(debugfail>=2) dualfprintf(fail_file,"uc2general: mypflag=%d numavg=%d startpl=%d endpl=%d :: i=%d j=%d k=%d\n",mypflag,numavg,startpl,endpl,i,j,k);

          
  if( mypflag==UTOPRIMFAILU2AVG2){
    numavg++; // assume at least always one good one so don't treat as real failure if no good values surrounding
  }
  // else use real numavg
          
          
  if(numavg==0){
    return(1);
  }
          
  // good value exist
  return(0);
}




// simple average of good surrounding zones
static int simple_average(int startpl, int endpl, int i, int j, int k,PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;

  if(debugfail>=2) dualfprintf(fail_file,"uc2simple: i=%d j=%d k=%d\n",i,j,k); // not too much



  /////////////
  //
  // 4 VALUES
  //
  /////////////
  if( // but if surrounded by good values
     (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jp1mac(j),k)))&&
     (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jm1mac(j),k)))&&
     (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),j,k)))&&  
     (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),j,k)))    
      ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected1\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    // don't mess with B field, it's correct already
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG4_1(ptoavg,i,j,k,pl);  
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if surrounded by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jp1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jm1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jp1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected2\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG4_2(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  /////////////
  //
  // 2 VALUES
  //
  /////////////
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jp1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected3\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG2_1(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),j,k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),j,k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected4\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG2_2(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jp1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected5\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG2_3(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jm1mac(j),k)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jp1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected6\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=AVG2_4(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  /////////////
  //
  // SINGLE VALUES
  //
  /////////////
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jp1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),jp1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),j,k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,ip1mac(i),jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,i,jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jm1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),j,k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,im1mac(i),jp1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),jp1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A0(lpflagfailorig,i,jp1mac(j),k)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    for(pl=startpl;pl<=endpl;pl++){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,i,jp1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else{
    return(1);
  }

  // reaches here if ok
  return(0);
}





#define NOGOODSTATIC 0
#define NOGOODRESET 1
#define NOGOODAVERAGE 2

// no good values.  Choices: STATIC, RESET to something, and AVERAGE failed (t-dt surrounding) values.
#define TODONOGOOD NOGOODAVERAGE
// 0=STATIC // model-independent, but really bad idea -- causes death of computation for accretion disks in jet region.
// 1=RESET // model-dependent
// 2=AVERAGE // model-independent -- probably best idea -- and then hope for the best.


// as stated
// 0 : pv
// 1: pbackup
#define WHICHPTOUSEWHENNOGOOD 0


// how to fixup when no good surrounding values to use
static int fixup_nogood(int startpl, int endpl, int i, int j, int k, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom)
{
  int pl,pliter;
  FTYPE prguess[NPR];
  FTYPE (*ptoavgwhennogood)[NSTORE2][NSTORE3][NPR];
  int ti,tj,tk,resetregion;



  if(debugfail>=2) dualfprintf(fail_file,"uc2nogood: i=%d j=%d k=%d\n",i,j,k); // not too much


  // Determine which primitive to use to average when no good solution
  // exists around to average
  if(WHICHPTOUSEWHENNOGOOD==0){
    ptoavgwhennogood=ptoavg;
  }
  else if(WHICHPTOUSEWHENNOGOOD==1){
    ptoavgwhennogood=pbackup;
  }


  ti=startpos[1]+i;
  tj=startpos[2]+j;
  tk=startpos[3]+k;


  /////////////
  //
  // no good values.  Choices: STATIC, RESET to something, and AVERAGE failed (t-dt surrounding) values.
  //
  /////////////


  /////////////////////////////////////////////////
  //
  /////// BELOW STUFF IS MODEL DEPENDENT
  //
  /////////////////////////////////////////////////

#if(USERRESETREGION==1)
  // tj = 0,1,totalsize[2]-2,totalsize[2]-1 are reset region
  // ti<10 near BH is reset region
  //  resetregion=(tj>=-2 && tj<=1 || tj>=totalsize[2]-2 && tj<=totalsize[2]+1) || (ti<10);
  resetregion=(tj<=1 || tj>=totalsize[2]-2) || (ti<10);
#elif(USERRESETREGION==0)
  resetregion=0;
#endif



  if(resetregion || TODONOGOOD==NOGOODRESET){// if no solution, revert to normal observer and averaged densities

    // model-dependent: assumes failures occur mostly in jet region near where matter is mostly freefalling near black hole where b^2/rho_0>>1

    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected9\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // otherwise some surrounding solutions are also bad, so just keep static
    // keeping static can lead to large regions of high gamma that are static and unstable to interactions.  Thus, reset static to normal observer flow for safety.
    // setting to normal observer for surrounding high gamma flows leads to large shocks and high internal energy.
    //
    // thus not sure what to do. (should have gotten correct U in first place! -- that's what you should do)
    //
    // limit_gamma?  Could also reset v completely to normal observer
    // if(limit_gamma(1.0,MAC(pv,i,j,k),ptrgeom,-1)>=1)
    //    FAILSTATEMENT("fixup.c:fixup_utoprim()", "limit_gamma()", 1);

    // average out densities
    for(pl=0;pl<=UU;pl++) MACP0A1(pv,i,j,k,pl)=0.5*(AVG4_1(ptoavgwhennogood,i,j,k,pl)+AVG4_2(ptoavgwhennogood,i,j,k,pl));
    // make velocity the normal observer
    set_atmosphere(0,WHICHVEL,ptrgeom,prguess);
    for(pl=U1;pl<=U3;pl++) MACP0A1(pv,i,j,k,pl)=prguess[pl];
    // average out things above field "pl"
    for(pl=B3+1;pl<NPR;pl++) MACP0A1(pv,i,j,k,pl)=0.5*(AVG4_1(ptoavgwhennogood,i,j,k,pl)+AVG4_2(ptoavgwhennogood,i,j,k,pl));

  }
  /////////////////////////////////////////////////
  //
  /////// BELOW STUFF IS NOT MODEL DEPENDENT
  //
  /////////////////////////////////////////////////
  else if(TODONOGOOD==NOGOODAVERAGE){// if no solution, revert to normal observer and averaged densities
    // average out densities+velocity
    for(pl=startpl;pl<=endpl;pl++) MACP0A1(pv,i,j,k,pl)=0.5*(AVG4_1(ptoavgwhennogood,i,j,k,pl)+AVG4_2(ptoavgwhennogood,i,j,k,pl));
  }
  else if(TODONOGOOD==NOGOODSTATIC){// if no solution, revert to normal observer and averaged densities
    // don't change -- stay at previous timestep value (or whatever utoprim() left it as).
  }
  else{
    dualfprintf(fail_file,"No condition for failure in fixup_nogood()\n");
    myexit(348576346);
  }

  return(0);
}




// OPENMPNOTE: Assume superdebug_utoprim() not called by multiple threads (i.e. when debugging, not using multiple threads), so static's are ok (including firsttime)
static int superdebug_utoprim(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled)
{
  int pl,pliter;
  struct of_state q;
  FTYPE Ui[NPR],Uf[NPR];
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;
  int failreturn;
  static FILE * super_fail_file;
  static int firsttime=1;
  int output;
  static int countoutput=0;

#if(USEOPENMP)
  dualfprintf(fail_file,"Cannot Superdebug with OpenMP\n");
  myexit(45968546);
#endif


  if(firsttime){
    super_fail_file=fopen("superdebug.out","wt");
    if(super_fail_file==NULL){
      dualfprintf(fail_file,"Cannot open superdebug.out\n");
      myexit(1);
    }
    countoutput=0;
  }

  if(whocalled==-1){
    // then only output if every so often, randomly in step
    if(ranc(0,0)>0.99) output=1;
    else output=0;
  }
  else output=1;

  if(output){
    // before any changes
    failreturn=get_state(pr0,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Ui);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
    
    // after any changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Uf);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");


    coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
    bl_coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, V);
    
    
    // used to determine nature of pre and post failure quantities
    fprintf(super_fail_file,"%21.15g %ld %d %d %d ",t,realnstep,steppart,numstepparts,whocalled);
    fprintf(super_fail_file,"%21.15g %21.15g %21.15g %d %d %d ",V[1],V[2],V[3],startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k);
    PLOOP(pliter,pl) fprintf(super_fail_file,"%21.15g %21.15g %21.15g %21.15g ",pr0[pl],pr[pl],Ui[pl],Uf[pl]);
    fprintf(super_fail_file,"\n");
    if(!(countoutput%1000)) fflush(super_fail_file);
    countoutput++;
  }



  firsttime=0;

  return(0);
}



int set_density_floors_default(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  struct of_state q;
  FTYPE U[NPR];
  FTYPE Upbinf[NPR];
  FTYPE r,th,X[NDIM],V[NDIM];
  FTYPE bsq;
  int pl,pliter;

  coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
  bl_coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, V);
  r=V[1];
  th=V[2];



  ////////////////////
  // scaling functions
  ////
  if(1){ // choice
    if(rescaletype==0){
      // to be like Bondi

      // even if d/dt->0 and d/dphi->0, this doesn't conserve E and L along magnetic flux lines in the funnel region.  This is bad!
      prfloor[RHO] = prfloorcoef[RHO]*pow(r, -1.5); 
      prfloor[UU] = prfloorcoef[UU]*pow(r,-2.5);
    }
    else if(rescaletype==1){
      // to conserve E and L along magnetic lines, have b^2/\rho\sim const.
      prfloor[RHO] = prfloorcoef[RHO]*pow(r, -2.7);
      prfloor[UU] = prfloorcoef[UU]*pow(r,-2.7) ;
    }
    else if(rescaletype==2){
      // only set absolute limits and let density be set later

      // absolute limits determined by what density contrast can be handled by flux calculation.  This is at least limited by machine precision, but more likely an instability is generated for some  ratio of densities that's much smaller than machine precision.

      // we should fix absolute floor so that this point is not too much different than surrounding points.
      // this is difficult because process depends on how done, so just set reasonable minimum

      // Bondi like
      // coefficient should be small enough so that at outer radius the energy per baryon at infinity limit of (say) 100 can be maintained (b^2 \sim 0.01 r^{-2.7} for \rho=1 density maximum)
      // these will do for rout=10^4
      //      prfloor[RHO] = 1E-12*pow(r, -1.5);
      //prfloor[UU] = 1E-14*prfloor[RHO]/r ;
      


      // to best conserve E and L along magnetic field lines
      if (get_state(pr, ptrgeom, &q) >= 1)
        FAILSTATEMENT("fixup.c:set_density_floors()", "get_state() dir=0", 1);

      MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, U),"fixup.c:set_density_floors()", "primtoU()", 1);

      // now have U[UU] and U[PH]
      // note that U[UU]/(gdet*rho*u^t) is conserved along field lines, even in non-ideal MHD, where A_{\phi} is still a good stream function in non-ideal MHD.
      // same in terms of U[PH]
      // this is really a coordinate dependent quantity, but if field lines connect to infinity, then this is the same at infinity.

      // Conserved quantity per baryon (or rest-mass unit)
      //      PALLLOOP(k) Upbinf=U[k]/(ptrgeom->gdet * pr[RHO]*(q.ucon[TT]));
      //Upbinf[UU]*=-1.0; // time component has - sign in this -+++ signature code

      // could inject mass or enthalpy, we choose mass for now since rho>h typically (except at very edge of torus)
      prfloor[RHO]=-U[UU]/(ptrgeom->gdet * GAMMAMAX * (q.ucon[TT]));
      prfloor[UU]=prfloor[RHO]*0.01; // cold injection

    }
    else if(rescaletype==3){
      // for dipole
      prfloor[RHO] = prfloorcoef[RHO]*pow(r/Rin, -5);
      prfloor[UU] = prfloorcoef[UU]*pow(r/Rin,-6) ;
    }
    else if(rescaletype==4){
      // for jet injection with maximum b^2/rho and b^2/u and maximum u/rho
      
      if(bsq_calc(pr,ptrgeom,&bsq)>=1){
        dualfprintf(fail_file,"bsq_calc:bsq_calc: failure\n");
        return(1);
      }
      prfloor[UU]=MAX(bsq/BSQOULIMIT,zerouuperbaryon*MAX(pr[RHO],SMALL));
      // below uses max of present u and floor u since present u may be too small (or negative!) and then density comparison isn't consistent with final floor between u and rho
      prfloor[RHO]=MAX(MAX(bsq/BSQORHOLIMIT,max(pr[UU],prfloor[UU])/UORHOLIMIT),SMALL);
    }
  }
  else{
    prfloor[RHO] = RHOMIN*pow(r, -2.0);
    prfloor[UU] = UUMIN*prfloor[RHO] ;
  }






  // some old attempts
  /*
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

    //    else{
    if(1){ 
    rhotscal = 1.0;
    uutscal = 1.0;
    }


    prfloor[UU] = UUMIN*uurscal*uutscal;
    prfloor[RHO] = RHOMIN*rhorscal*rhotscal;
  */


  return(0);
}



#define FLOORDAMPFRAC (0.1)
#define NUMBSQFLAGS 5

int get_bsqflags(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE bsq;
  struct of_state q;
  FTYPE gamma,X[NDIM],V[NDIM];
  FTYPE qsq;
  FTYPE prfloor[NPR];
  int flags[NUMBSQFLAGS]={0};
  int limgen;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int loc=CENT;




  //  if((CHECKSOLUTION==0)&&(LIMADJUST==0)&&(FLUXADJUST==0)) return(0); // otherwise do it
  // fixup_checksolution() only uses failure flag right now
  if((LIMADJUST==0)&&(FLUXADJUST==0)) return(0); // otherwise do it



  // temporary hack
  limgen=MAX(MAX(lim[1],lim[2]),lim[3]); // not reached if((LIMADJUST==0)&&(FLUXADJUST==0))



  COMPDQZLOOP { // over range wspeed and dq are computed // OPENMPOPTMARK: Could optimize this, but not using it currently

    get_geometry(i, j, k, loc, ptrgeom);


#if(DOEVOLVERHO)
    // b^2
    if (get_state(MAC(pv,i,j,k), ptrgeom, &q) >= 1)
      FAILSTATEMENT("fixup.c:get_bsqflags()", "get_state()", 1);
    bsq = dot(q.bcon, q.bcov);
    // initial
    GLOBALMACP0A1(pflag,i,j,k,FLAGBSQORHO)=0;

    // now do checks 10/31/2003 : constants fixed based upon accretion
    // models and locations of disk, corona, and funnel (b^2/rho>1)

    // b^2/\rho
    if(bsq/MACP0A1(pv,i,j,k,RHO)>BSQORHOLIMIT) flags[0]=2;
    else if(bsq/MACP0A1(pv,i,j,k,RHO)>0.9*BSQORHOLIMIT) flags[0]=1;
    else flags[0]=0;

#endif

#if(DOEVOLVEUU)
    // b^2/u

    if(bsq/MACP0A1(pv,i,j,k,UU)>BSQOULIMIT) flags[3]=2;
    else if(bsq/MACP0A1(pv,i,j,k,UU)>0.9*BSQOULIMIT) flags[3]=1;
    else flags[3]=0;
#endif



#if(WHICHVEL==VELREL4)
    // gamma
    if(gamma_calc(MAC(pv,i,j,k),ptrgeom,&gamma,&qsq)>=1){
      dualfprintf(fail_file,"gamma calc failed: get_bsqflags\n");
      if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    // \gamma relative to 4-velocity
    if(gamma>2.0*GAMMADAMP) flags[1]=2;
    else if(gamma>GAMMADAMP) flags[1]=1;
    else flags[1]=0;
#endif
    // these density mod on limiters may be too much diffusivity in otherwise ok regions


#if(0)
    // density floors
    set_density_floors(ptrgeom,MAC(pv,i,j,k),prfloor);
    // rho/rho_{fl}    
    // GODMARK, 0.1 here    
    if(prfloor[RHO]/MACP0A1(pv,i,j,k,RHO)>FLOORDAMPFRAC) flags[2]=2;
    else if(prfloor[RHO]/MACP0A1(pv,i,j,k,RHO)>0.1*FLOORDAMPFRAC) flags[2]=1;
    else flags[2]=0;

    // u/u_{fl}

    if(prfloor[UU]/MACP0A1(pv,i,j,k,UU)>FLOORDAMPFRAC) flags[4]=2;
    else if(prfloor[UU]/MACP0A1(pv,i,j,k,UU)>0.1*FLOORDAMPFRAC) flags[4]=1;
    else flags[4]=0;
#endif


    // now check our flag state
    GLOBALMACP0A1(pflag,i,j,k,FLAGBSQORHO)=flags[0];
    GLOBALMACP0A1(pflag,i,j,k,FLAGBSQOU)=flags[3];


#if(LIMADJUST>0)
    // now check our flag state
    if((flags[0]==2)||(flags[1]==2)||(flags[2]==2)||(flags[3]==2)||(flags[4]==2)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALLIM)=limgen-1; }
    else if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALLIM)=limgen-1; }
    else{GLOBALMACP0A1(pflag,i,j,k,FLAGREALLIM)=limgen;}
    if(GLOBALMACP0A1(pflag,i,j,k,FLAGREALLIM)<0) GLOBALMACP0A1(pflag,i,j,k,FLAGREALLIM)=0;
#endif

#if(FLUXADJUST>0)
    // now check our flag state
    if((flags[0]==2)||(flags[1]==2)||(flags[2]==2)||(flags[3]==2)||(flags[4]==2)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod-1;}
    else if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod; }
    else{GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod;}
    if(GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)<0) GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=0;
#endif


  }
  return(0);
}



// whether to limit gamma inside ergosphere
#define GAMMAERGOLIMIT 0

// whether to conserve (at least) D=\rho_0 u^t when modifying \gamma
// risky to assume D=rho_0 u^t conserved when limiting \gamma because \gamma may be large due to error simply in T^t_\nu evolution alone
// So D evolution may be normal and more accurate, and so effectively error in T^t_\nu leads to large u^t but D changes little
// so this fix would add a great deal of rest-mass
// So for now just limit velocity ignoring all conservation laws
#define DO_CONSERVE_D 0

// limit \gamma=\alpha u^t and u^t
int limit_gamma(FTYPE gammamax, FTYPE*pr, FTYPE *ucons, struct of_geom *ptrgeom,int finalstep)
{
  FTYPE f,gamma,pref;
  FTYPE qsq;
  FTYPE alpha;
  FTYPE pr0[NPR];
  int pl,pliter;
  FTYPE realgammamax;
  FTYPE r,X[NDIM],V[NDIM];
  int didchange;
  FTYPE uu0,uu0max;
  int i=ptrgeom->i;
  int j=ptrgeom->j;
  int k=ptrgeom->k;
  int loc=ptrgeom->p;




  // assume didn't change primitives
  didchange=0;



  /////////////////////////////////
  //
  // Ad-hoc ergo fix for force-free black hole problem (enabled very rarely)
  //
  /////////////////////////////////
#if(GAMMAERGOLIMIT)
  coord_ijk(i,j,k,loc,X);
  bl_coord_ijk(i,j,k,loc,V);
  r=V[1];

  // force flow to not move too fast inside ergosphere
  if(r<2) realgammamax=3;
  else realgammamax=GAMMAMAX;
#else
  realgammamax=gammamax;
  if(realgammamax<=1.0) return(0); // nothing to do
#endif






  /////////////////////////////////
  //
  // Get \gamma
  //
  /////////////////////////////////
  //PALLLOOP(pl){ dualfprintf(fail_file,"i=%d j=%d k=%d pr[%d]=%21.15g\n",startpos[1]+i,startpos[2]+j,startpos[3]+k,pl,MAC(pv,i,j,k));}
  if(gamma_calc(pr,ptrgeom,&gamma,&qsq)>=1){
    dualfprintf(fail_file,"limit_gamma: gamma calc failed\n");
    dualfprintf(fail_file,"i=%d j=%d k=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,gamma);
    if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  


  ////////////////////////////////
  //
  // set pre-changed primitive
  //
  ////////////////////////////////
  PALLLOOP(pl)    pr0[pl]=pr[pl];


    
  ////////////////////////////////
  //
  // If \gamma>\gammamax, then force \gamma=\gammamax
  //
  ////////////////////////////////
  if((gamma > realgammamax && (gamma!=1.0))) {    

    // rescale velocities to reduce gamma to realgammamax
    pref=(realgammamax*realgammamax - 1.)/(gamma*gamma - 1.);

    if(debugfail>=3){ // special >=3 since often called when floor or fixup
      dualfprintf(fail_file,"nstep=%ld steppart=%d t=%21.15g :: i=%d j=%d k=%d :: pref=%21.15g oldgamma=%21.15g realgammamax=%21.15g\n",nstep,steppart,t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pref,gamma,realgammamax);
    }

    if(pref<0.0){
      dualfprintf(fail_file,"limit_gamma: pref calc failed pref=%21.15g\n",pref);
      dualfprintf(fail_file,"i=%d j=%d k=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,gamma);
      if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }

    f = sqrt(pref);
    pr[U1] *= f ;       
    pr[U2] *= f ;       
    pr[U3] *= f ;


#if(DO_CONSERVE_D)
    //    alpha = alpha = 1./sqrt(-ptrgeom->gcon[GIND(TT,TT)]) ;
    //uu0old=gamma/alpha;
    // force conservation of particle number
    pr[RHO] = pr0[RHO]*gamma/realgammamax; // can do this since alpha is constant and so cancels
#endif

    gamma=gammamax; // reset gamma for next check
    didchange=1; // indicate did change primitive



  }




  ///////////////////
  //
  // limiting u^t makes no sense except as an indicator of when in difficult regime
  //
  ///////////////////
#if(0)
  // limit u^t too since maybe alpha very small
  //alpha = 1./sqrt(-ptrgeom->gcon[GIND(TT,TT)]) ;
  alpha = ptrgeom->alphalapse;
  uu0=gamma/alpha;
  uu0max=realgammamax;
  // since alpha is always >=1, then limit on uu0 always overrides limit on gamma?
  // here realgammamax is u^t not u^t\alpha
  if((uu0 > uu0max && (uu0!=1.0))) {    

    //dualfprintf(fail_file,"gamma=%21.15g realgammamax=%21.15g\n",gamma,gammamax);

    if(debugfail>=2) dualfprintf(fail_file,"i=%d j=%d k=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,gamma);
    // rescale velocities to reduce gamma to realgammamax
    pref=(uu0max*uu0max*alpha*alpha - 1.)/(gamma*gamma - 1.);
    if(pref<0.0){
      dualfprintf(fail_file,"limit_gamma: pref calc failed pref=%21.15g\n",pref);
      dualfprintf(fail_file,"i=%d j=%d k=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,gamma);
      if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }

    f = sqrt(pref);
    pr[U1] *= f ;       
    pr[U2] *= f ;       
    pr[U3] *= f ;

#if(DO_CONSERVE_D)
    //    alpha = alpha = 1./sqrt(-ptrgeom->gcon[GIND(TT,TT)]) ;
    //uu0old=gamma/alpha;
    // force conservation of particle number
    pr[RHO] = pr0[RHO]*uu0/uu0max;
#endif

    didchange=1;
  }
#endif




  ///////////////////
  //
  // Account for changes in conserved quantities via changes in: \rho_0 and pr[U1..U3]
  //
  ///////////////////
  if(didchange){
    FTYPE prdiag[NPR];
    PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
    diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTLIMITGAMMAACT);
    PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    return(-1);// indicates did change primitive
  }


  return(0); // indicates didn't change anything
}







// check_pr() is purely 2D function

#define UTLIMIT (50.0) // limit of u^t given no other info and need to fix u
#define UTFIX (utobs) //  if failed u^t, then fix to this level to keep normal, or use local value as desired value
#define GRADIENTFACTOR (.001) // how to set?

// pr is pr to be fixed
// prmodel is pr that has a model ucon[TT], since we want to interpolate based upon ucon[TT], not just a random one when inside step_ch.c in fluxcalc()
// prmodel is just pr if in utoprim.c since we want a real change, not just an interpolation, thus will not be used as basis
int check_pr(FTYPE *pr,FTYPE *prmodel, FTYPE *ucons, struct of_geom *ptrgeom,int modelpos, int finalstep)
{

#if(WHICHVEL==VEL3)
  int failedcheck;
  FTYPE ucon[NPR],uconmodel[NPR];
  FTYPE others[NUMOTHERSTATERESULTS],othersmodel[NUMOTHERSTATERESULTS];
  FTYPE prold[NPR],probs[NPR];
  FTYPE gradient[NDIM],normsq;
  int trialcount,ntrials;
  int i,j,k;
  int pl,pliter;
  FTYPE lastuttdiscr,uttdiscr0,utobs;
  FTYPE dampfactor,dampfactorchange;
  FTYPE newerr,olderr;
  int method;
  struct of_geom modelgeomdontuse;
  struct of_geom *ptrmodelgeom=&modelgeomdontuse;
  int idel,jdel,kdel;
  FTYPE realutlimit,realdiscrlimit;
  FTYPE modeldiscr;
  int utinterp;
  FTYPE toldiscr;
  int modeli,modelj,modelk;
  int bctype;
  FTYPE pr0[NPR];
  int loc=CENT;


  if(modelpos==-100){ // then utoprim called us from usrfun
    modelpos=0;
    ntrials=30; // don't try so hard since failure is more likely, leads to static solution
  }
  else{
    ntrials=200; // try really hard since using observer as solution is not fun
  }
  toldiscr=1.E-2;
  dampfactorchange=0.5;
  dampfactor=1.0;
  method=1; // 0=GRADIENTFACTOR 1=damp newton's method
  utinterp=0; // whether to fix interpolation (if bad) based upon a model pr
  

  if(method==1){
    // save old pr
    PALLLOOP(pl) probs[pl]=prold[pl]=pr[pl];
  }

  whocalleducon=1; // turn off failures
  if(ucon_calc(pr, ptrgeom, ucon,others) >= 1) ucon[TT]=1E30; // bad, so keep going
  uttdiscr0=lastuttdiscr=uttdiscr; // ucon_calc() set uttdiscr
  if(ucon[TT]<UTLIMIT) return(0); // good, so just continue calculations


  /////////////////////////////////////////////////////
  //
  // if here, need to fix
  //
  // set pre-primitive
  PALLLOOP(pl)    pr0[pl]=pr[pl];



  if(modelpos==-2) return(1); // don't try to correct, just fail since ucon[TT] wasn't less than the limit (otherwise wouldn't reach here).  Used to just check, no fix.
  
  if(modelpos==-3){
    modelpos=-1;
    bctype=1; // so bound_prim called us
  }
  else bctype=0; // non-bound call
  

  //////////////////
  //
  // if we are here, original velocities need to be corrected

  // first find normal observer to get relevant u^t
  for(i=1;i<=3;i++){
    probs[U1+i-1] = ptrgeom->gcon[GIND(TT,i)]/ptrgeom->gcon[GIND(TT,TT)] ;
  }
  // ucon[TT] must be good for normal observer (?) hard to solve for u^t for normal observer in general
  // change of whocalleducon forces kill level check on normal observer
  whocalleducon=0;
  if(ucon_calc(probs, ptrgeom, ucon, others) >= 1){
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
      modelk=ptrgeom->k;
      ptrmodelgeom=ptrgeom;
    }
    // modelpos>=1 for interpolations in fluxcalc() (model value on center)
    else if(modelpos==1){
      if(ptrgeom->p==FACE1){ idel=1; jdel=0; kdel=0; }
      else if(ptrgeom->p==FACE2){ idel=0; jdel=1; kdel=0; }
      else if(ptrgeom->p==FACE3){ idel=0; jdel=0; kdel=1; }
      modeli=(ptrgeom->i) -idel;
      modelj=(ptrgeom->j) -jdel;
      modelk=(ptrgeom->k) -kdel;
      // then i-1,j is center position
      get_geometry(modeli,modelj,modelk,loc,ptrmodelgeom);
    }
    else if(modelpos==2){
      modeli=ptrgeom->i;
      modelj=ptrgeom->j;
      modelk=ptrgeom->k;
      // then i,j is center position
      get_geometry(modeli,modelj,modelk,loc,ptrmodelgeom);
    }
    // determine model u^t
    if(ucon_calc(prmodel, ptrmodelgeom, uconmodel, othersmodel) >= 1){
      // model no good
      if(bctype==0){
        dualfprintf(fail_file,"serious failure.  On-grid values and fixed bc values used have u^t imaginary: modeli: %d modelj: %d uttdiscr: %21.15g\n",startpos[1]+modeli,startpos[2]+modelj,uttdiscr);
        whocalleducon=0; // turn on failures
        if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
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
      gradient[i]=2.0*(ptrgeom->gcov[GIND(0,i)]);
      for(j=1;j<=3;j++){
        // note that ucon is the same as pr here since ucon_calc sets spatial terms to pr
        gradient[i]+=2.0*ucon[j]*ptrgeom->gcov[GIND(i,j)];
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
        if(realdiscrlimit-uttdiscr>0)   pr[U1+i-1]-=gradient[i]*dampfactor;
        else    pr[U1+i-1]+=gradient[i]*dampfactor;
      }
    }
    // get new ucon[TT], is it ok now?
    if(ucon_calc(pr, ptrgeom, ucon,others) >= 1) ucon[TT]=1E30; // bad bad, keep going
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

  if(failedcheck){
    if(debugfail>=1) dualfprintf(fail_file,"couldn't fix ucon[TT]=%21.15g newerr=%21.15g uttdiscr=%21.15g\n",ucon[TT],newerr,uttdiscr);
    if(debugfail>=1) dualfprintf(fail_file,"check_pr failure: t=%21.15g , couldn't fix ucon: i=%d j=%d k=%d p=%d ucon[TT]=%21.15g realutlimit=%21.15g uconmodel[TT]=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,ptrgeom->p,ucon[TT],realutlimit,uconmodel[TT]);
    if(debugfail>=1) dualfprintf(fail_file,"uttdiscr0=%21.15g uttdiscr=%21.15g realdiscrlimit=%21.15g modeldiscr=%21.15g dampfactor=%21.15g\n",uttdiscr0,uttdiscr,realdiscrlimit,modeldiscr,dampfactor);
    // will still run perhaps with UT>realutlimit, could stop it, but won't for now      
    if(debugfail>=1){
      PALLLOOP(pl){
        dualfprintf(fail_file,"pr[%d]=%21.15g prmodel[%d]=%21.15g\n",pl,pr[pl],pl,prmodel[pl]);
      }
      if(debugfail>=1) dualfprintf(fail_file,"need better algorithm: check_pr failure: couldn't fix ucon: i=%d j=%d k=%d p=%d ucon[TT]=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,ptrgeom->p,ucon[TT]);
    }

    // force to normal observer solution
    PALLLOOP(pl) pr[pl]=probs[pl];
  }

  // account for changes
  FTYPE prdiag[NPR];
  PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
  diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTLIMITGAMMAACT);
  PLOOP(pliter,pl) prdiag[pl]=pr[pl];

#endif

  return(0);
}


// GODMARK: check this function for correctness
int inflow_check_4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{
  int ii,jj,kk;
  int iin,iout;
  int jjn,jout;
  int kkn,kout;
  FTYPE pr0[NPR],prdiag[NPR];
  int pl,pliter;


  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;


  if(dir==1){
    // for dir==1
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN1) ){ iin=-1; iout=totalsize[1]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN3)||(ptrgeom->p==CORN2) ){ iin=0; iout=totalsize[1]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[U1]=0;

      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==FIXEDOUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==FIXEDOUTFLOW)&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      // then inflow according to Bondi-like atmosphere
      set_atmosphere(1,WHICHVEL,ptrgeom,pr);

      // below never really accounted for since on boundary zones
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
  }
  else if(dir==2){
    // for dir==2
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE1)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN2) ){ jjn=-1; jout=totalsize[2]; }
    else if((ptrgeom->p==FACE2)||(ptrgeom->p==CORN1)||(ptrgeom->p==CORN3) ){ jjn=0; jout=totalsize[2]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW)&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[U2]=0;

      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
  }
  else if(dir==3){
    // for dir==3
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN3) ){ kkn=-1; kout=totalsize[3]; }
    else if((ptrgeom->p==FACE3)||(ptrgeom->p==CORN1)||(ptrgeom->p==CORN2) ){ kkn=0; kout=totalsize[3]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[3]+kk<=kkn)&&(BCtype[X3DN]==OUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[3]+kk>=kout)&&(BCtype[X3UP]==OUTFLOW)&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[U3]=0;

      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
  }
  else return(1); // uh

  return(0);
}


int inflow_check_3vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{

  return(inflow_check_4vel(dir,pr,ucons, ptrgeom,-1));

}

// GODMARK: check for correctness
int inflow_check_rel4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{
  FTYPE ucon[NDIM] ;
  FTYPE others[NUMOTHERSTATERESULTS];
  int ii,jj,kk,loc;
  int j,k ;
  int pl,pliter;
  FTYPE alpha,betacon,gamma,vsq ;
  FTYPE qsq;
  int iin,iout;
  int jjn,jout;
  int kkn,kout;
  int dofix=0;
  FTYPE pr0[NPR];


  //  return(0);
  ii=ptrgeom->i; 
  jj=ptrgeom->j;
  kk=ptrgeom->k;
  loc=ptrgeom->p;


  //ucon_calc(pr, ptrgeom, ucon,others) ;

  //  PALLLOOP(pl) dualfprintf(fail_file,"nstep=%ld steppart=%d t=%21.15g :: pl=%d %21.15g\n",nstep,steppart,t,pl,pr[pl]);

  MYFUN(ucon_calc(pr, ptrgeom, ucon, others),"fixup.c:inflow_check_rel4vel()", "ucon_calc() dir=0", 1);

  //  DLOOPA(pl) dualfprintf(fail_file,"ucon[%d]=%21.15g\n",pl,ucon[pl]);

  dofix=0;
  if(dir==1){
    // for dir==1
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN1) ){ iin=-1; iout=totalsize[1]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN2)||(ptrgeom->p==CORN3) ){ iin=0; iout=totalsize[1]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==OUTFLOWNOINFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==OUTFLOWNOINFLOW)&&(ucon[dir] < 0.)) 
        ) {
      dofix=1;
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==FIXEDOUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==FIXEDOUTFLOW)&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      // then inflow according to Bondi-like atmosphere
      // GODMARK: need to ensure this gives well-defined answers during init.c processing
      set_atmosphere(1,WHICHVEL,ptrgeom,pr); // can change all pr's
      dofix=0; // assume in boundary
    }
  }
  else if(dir==2){
    // for dir==2
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE1)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN2) ){ jjn=-1; jout=totalsize[2]; }
    else if((ptrgeom->p==FACE2)||(ptrgeom->p==CORN1)||(ptrgeom->p==CORN3) ){ jjn=0; jout=totalsize[2]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW || BCtype[X2DN]==OUTFLOWNOINFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW || BCtype[X2UP]==OUTFLOWNOINFLOW)&&(ucon[dir] < 0.)) 
        ) {
      dofix=2;
    }
  }
  else if(dir==3){
    // for dir==3
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE1)||(ptrgeom->p==FACE2)||(ptrgeom->p==CORN3) ){ kkn=-1; kout=totalsize[3]; }
    else if((ptrgeom->p==FACE3)||(ptrgeom->p==CORN1)||(ptrgeom->p==CORN2) ){ kkn=0; kout=totalsize[3]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[3]+kk<=kkn)&&(BCtype[X3DN]==OUTFLOW || BCtype[X3DN]==OUTFLOWNOINFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[3]+kk>=kout)&&(BCtype[X3UP]==OUTFLOW || BCtype[X3UP]==OUTFLOWNOINFLOW)&&(ucon[dir] < 0.)) 
        ) {
      dofix=3;
    }
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach dir=%d\n",dir);
    return(1); // uh
  }



  if(dofix){
    // set pre-primitive
    PALLLOOP(pl)    pr0[pl]=pr[pl];
    FTYPE prdiag[NPR];
    PLOOP(pliter,pl) prdiag[pl]=pr0[pl];


    /* find gamma and remove it from primitives */
    if(gamma_calc(pr,ptrgeom,&gamma,&qsq)>=1){
      dualfprintf(fail_file,"gamma calc failed: inflow_check_rel4vel\n");
      if (fail(ii,jj,kk,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    pr[U1] /= gamma ;
    pr[U2] /= gamma ;
    pr[U3] /= gamma ;
    alpha = 1./sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    
    /* reset radial velocity so radial 4-velocity
     * is zero */
    if(dofix==1){
      betacon = ptrgeom->gcon[GIND(0,1)]*alpha*alpha ;
      pr[U1] = betacon/alpha ; // gives 3-vel contravariant
    }
    else if(dofix==2){
      betacon = ptrgeom->gcon[GIND(0,2)]*alpha*alpha ;
      pr[U2] = betacon/alpha ; // gives 3-vel contravariant
    }
    else if(dofix==3){
      betacon = ptrgeom->gcon[GIND(0,3)]*alpha*alpha ;
      pr[U3] = betacon/alpha ; // gives 3-vel contravariant
    }
    /* now find new gamma and put it back in */
    vsq = 0. ;
    SLOOP(j,k) vsq += ptrgeom->gcov[GIND(j,k)]*pr[U1+j-1]*pr[U1+k-1] ;

    if(vsq<0.0){
      if(vsq>-NUMEPSILON*100.0) vsq=0.0; // machine precision thing
      else if (fail(ii,jj,kk,loc,FAIL_VSQ_NEG) >= 1){
        trifprintf("vsq=%21.15g\n",vsq);
        return (1);
      }
    }
 
    // it's possible that setting ucon(bc comp)->0 leads to v>c, so just reduce gamma to GAMMAMAX if that's the case
    if(vsq>=1.0){
      if(debugfail>=1) dualfprintf(fail_file,"i=%d j=%d k=%d inflow limit required change in gamma (after dofix==%d): vsq=%21.15g newvsq=%21.15g\n",ii+startpos[1],jj+startpos[2],kk+startpos[3],dofix,vsq,1.0-1.0/(GAMMAMAX*GAMMAMAX));

      // new vsq
      vsq = 1.0-1.0/(GAMMAMAX*GAMMAMAX);

    }

    gamma = 1./sqrt(1. - vsq) ;
    pr[U1] *= gamma ;
    pr[U2] *= gamma ;
    pr[U3] *= gamma ;

    // only for boundary conditions, not active zones, hence -1.0 instead of finalstep
    diag_fixup(1,prdiag, pr, ucons, ptrgeom, finalstep,COUNTINFLOWACT);
    PLOOP(pliter,pl) prdiag[pl]=pr[pl];

    /* done */
    return(-1);
  }
  else return(0);

}
 

// only correct for polar axis at both inner and outer x2/theta grid edge.
// OPENMPOPTMARK: Could optimize this, but not using it currently
void fix_flux(FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl ;
  FTYPE sth ;
  FTYPE X[NDIM],V[NDIM],r,th ;
  int inboundloop[NDIM];
  int outboundloop[NDIM];
  int innormalloop[NDIM];
  int outnormalloop[NDIM];
  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  int riin,riout,rjin,rjout,rkin,rkout;
  int dosetbc[COMPDIM*2];
  int ri;
  int boundvartype=BOUNDFLUXTYPE;

  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);
  //enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
  //  localenerpos=enerposreg[enerregion];


  // this has nothing to deal with MPI-boundaries, so ok as is
  // only applies for polar axis
  if(mycpupos[2]==0){    
    if(BCtype[X2DN]==POLARAXIS){
      LOOPX2dir{
        // emf should be antisymmetric around polar axes? // GODMARK: how does this mix with metric?
        LOOPBOUND2IN{
          MACP0A1(F1,i,j,k,B1) = 0;
          MACP0A1(F1,i,j,k,B2) = -MACP0A1(F1,i,-(jp1mac(j)),k,B2) ; // symmetric positions around polar axis, but antisymmetric value
          MACP0A1(F3,i,j,k,B2) = -MACP0A1(F3,i,-(jp1mac(j)),k,B2) ; // symmetric positions around polar axis, but antisymmetric value
          MACP0A1(F3,i,j,k,B3) = 0;
        }
        // all should be 0 except kinetic energy flux
        // GODMARK: I'm unsure if emf is not unlike, say, \Omega, which is a well-defined thing on the axis.
        PALLLOOP(pl) if(pl!=U2) MACP0A1(F2,i,0,k,pl) = 0. ;
      }
    }
  }
  if(mycpupos[2]==ncpux2-1){
    if(BCtype[X2UP]==POLARAXIS){
      LOOPX2dir{
        // emf
        LOOPBOUND2OUT{
          MACP0A1(F1,i,j,k,B2) = -MACP0A1(F1,i,jrefshiftmac(j),k,B2) ;
          MACP0A1(F3,i,j,k,B2) = -MACP0A1(F3,i,jrefshiftmac(j),k,B2) ;
        }
        // GODMARK: unsure
        PALLLOOP(pl) if(pl!=U2) MACP0A1(F2,i,N2,k,pl) = 0. ;
      }
    }
  }

  // avoid mass flux in wrong direction, so consistent with velocity fix
  // how to treat other fluxes?
  if(mycpupos[1]==0){    
    if(BCtype[X1DN]==OUTFLOW){
      LOOPX1dir{
        ri=riin;
        LOOPBOUND1IN{
          if(MACP0A1(F1,i,j,k,RHO)>0) MACP0A1(F1,i,j,k,RHO)=0;
        }
      }
    }
  }
  if(mycpupos[1]==ncpux1-1){
    if(BCtype[X1UP]==OUTFLOW){
      LOOPX1dir{
        LOOPBOUND1OUT{
          if(MACP0A1(F1,i,j,k,RHO)<0) MACP0A1(F1,i,j,k,RHO)=0;
        }
      }
    }
  }

}


// counter for EOS table lookup failures
void diag_eosfaillookup(int i, int j, int k)
{
  int whocalled = COUNTEOSLOOKUPFAIL;
  int tscale;

  // count every time corrects
  if(DODEBUG){
    if(i==AVOIDI && j==AVOIDJ && k==AVOIDK){
      // then do nothing since don't know where failure occurred
    }
    else if(i<-N1BND || i>N1-1+N1BND ||j<-N2BND || j>N2-1+N2BND ||k<-N3BND || k>N3-1+N3BND ){
      dualfprintf(fail_file,"In diag_eosfaillookup() whocalled=%d for i=%d j=%d k=%d\n",whocalled,i,j,k);
      myexit(3472762356);
    }
    int indexfinalstep;
    FINALSTEPLOOP(indexfinalstep) TSCALELOOP(tscale) GLOBALMACP0A3(failfloorcount,i,j,k,indexfinalstep,tscale,whocalled)++;
  }
}
