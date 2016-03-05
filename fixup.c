
/*! \file fixup.c
    \brief All routines related to fixing up solution when inversion fails or floors on densities are hit
    // all fixup stuff only called for non-B advance

*/



#include "decs.h"


static int simple_average(int startpl, int endpl, int i, int j, int k, int doingmhd, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3][NUMFAILPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR]);
static int general_average(int useonlynonfailed, int numbndtotry, int maxnumbndtotry, int startpl, int endpl, int i, int j, int k, int doingmhd, PFTYPE mhdlpflag, PFTYPE radlpflag, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3][NUMFAILPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom);
static int fixup_nogood(int startpl, int endpl, int i, int j, int k, int doingmhd, PFTYPE mhdlpflag, PFTYPE radlpflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom);
static int fixuputoprim_accounting(int i, int j, int k, PFTYPE mhdlpflag, PFTYPE radlpflag, int limitgammamhd, int limitgammarad, PFTYPE (*lpflag)[NSTORE2][NSTORE3][NUMPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *geom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
static int fixup_negdensities(int whichtofix, int *fixed, int startpl, int endpl, int i, int j, int k, PFTYPE mhdlpflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *geom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

static int superdebug_utoprim(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);



// whether to add radiation YFL4,YFL5 to gas YFL2,YFL3
#define YFLADDRADTOGAS 0



/// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
/// operations that only need to be done inside computational loop
int pre_fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR])
{

  // grab b^2 flags (above fixup may change u or rho, so must do this after)
  get_bsqflags(stage,pv);


  return(0);
}



/// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
/// this function actually changes primitives
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
    
    //    int long long nstepfake;
    //    nstepfake=nstep;
    //    nstep=0;
    //    image_dump(-4);
    //    nstep=nstepfake;

    // fixup before new solution (has to be here since need previous stage's failure flag)
    fixup_utoprim(stage,pv,pbackup,ucons,finalstep);

    // after fixup_utoprim(), the floors/ceilings won't be satisfied anymore, especially near sharp jumps in densities, so repeat.
    fixup(stage, pv, ucons, finalstep);

    
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


/// this function just reports problems, but doesn't fix them
int post_fixup_nofixup(int stageit, int finalstep, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{

  fixup_utoprim_nofixup(STAGEM1,pv,pbackup,ucons,finalstep);

  return(0);
}






/// apply floors to density, internal energy
/// currently called before bound, which assumes bound sets boundary
/// values exactly as wanted without any fixing.
int fixup(int stage,FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  int i, j, k;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;



  COMPZLOOP{
    //    dualfprintf(fail_file,"i=%d j=%d k=%d\n",i,j,k); fflush(fail_file);
    get_geometry(i,j,k,CENT,ptrgeom);
    if(fixup1zone(0,MAC(pv,i,j,k),MAC(ucons,i,j,k), ptrgeom,finalstep)>=1)
      FAILSTATEMENT("fixup.c:fixup()", "fixup1zone()", 1);
  }
  return(0);
}



/// choose whether within correctable/diagnosticable region
int diag_fixup_correctablecheck(int docorrectucons, struct of_geom *ptrgeom)
{
  int is_within_correctable_region;
  int docorrectuconslocal;

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


  return(docorrectuconslocal);


}



/// record who called the diag_fixup routine
int count_whocalled(int i, int j, int k, int finalstep, int whocalled, CTYPE toadd)
{
  int tscale;

  /////////////////////
  //
  // Count times in diag_fixup() and who called diag_fixup()
  // count every time corrects, not just on conserved quantity tracking time
  //
  /////////////////////
  if(DODEBUG){

    if(whocalled>=NUMFAILFLOORFLAGS || whocalled<COUNTNOTHING || i<-N1BND || i>N1-1+N1BND ||j<-N2BND || j>N2-1+N2BND ||k<-N3BND || k>N3-1+N3BND ){
      dualfprintf(fail_file,"In diag_fixup() whocalled=%d for i=%d j=%d k=%d\n",whocalled,i,j,k);
      myexit(24683463);
    }

    if(whocalled>=COUNTREALSTART){
      int indexfinalstep;
      indexfinalstep=0;
      TSCALELOOP(tscale) GLOBALMACP0A3(failfloorcount,i,j,k,indexfinalstep,tscale,whocalled)+=toadd;
      if(finalstep){
        indexfinalstep=1;
        // iterate finalstep version
        TSCALELOOP(tscale) GLOBALMACP0A3(failfloorcount,i,j,k,indexfinalstep,tscale,whocalled)+=toadd;
      }
    }// end if counting something

  }// end if DODEBUG


  return(0);

}



/// compute dU and account
/// Assumes Ui,Uf are UDIAG form
/// Assumes ucons is UEVOLVE form
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
    // don't modify ucons.  That is, we ignore docorrectuconslocal.

    PALLLOOP(pl) deltaUavg[pl] = Uf[pl]-Ui[pl];
  }



  //only do aggregate accounting, after the fact (just before taking the new time step)
  if(DOONESTEPDUACCOUNTING && whocalled==COUNTONESTEP ||  DOONESTEPDUACCOUNTING==0){

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

          fladd[pl] += dUincell[pl];

        }// end over pl's
      }// end if within diagnostic region

    }// end over enerregions

  }// end if doing accounting

  if(whocalled>=0 && whocalled!=COUNTONESTEP){

    //////////////
    //
    // Loop over ENERREGIONs
    //
    //////////////
    ENERREGIONLOOP(enerregion){

      // setup pointers to enerregion diagnostics
      enerpos=enerposreg[enerregion];
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
     

        }// end over pl's
      }// end if within diagnostic region

    }// end over enerregions

  }// end if doing accounting


  return(0);
}



int consfixup_allzones(int finaluu, FTYPE (*pf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{


  int i, j, k, pliter,pl;
  struct of_geom *ptrgeom;
  struct of_geom geomdontuse;

   
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

    consfixup_1zone(finaluu, i,j,k, ptrgeom, &MACP0A1(pf,i,j,k,0), &MACP0A1(ucons,i,j,k,0));


  }


  return(0);

}


int consfixup_1zone(int finaluu, int i, int j, int k, struct of_geom *ptrgeom, FTYPE *pf, FTYPE *ucons)
{
  if(ENSURECONS==0) return(0);

  int pliter,pl;
  // now that have averaged values, try to enforce energy conservation by borrowing from radiation
  FTYPE *pp=pf;
  struct of_state qp;
  get_state(pp,ptrgeom,&qp);
  int uutype=UNOTHING;
  FTYPE uu[NPR],uuabs[NPR];
  primtoU(uutype,pp,&qp,ptrgeom,uu,uuabs);
  FTYPE uu0[NPR];
  if(finaluu==1){
#if(WHICHEOM==WITHGDET)
    PLOOP(pliter,pl) uu0[pl]=ucons[pl]*ptrgeom->igdetnosing; // put in UNOTHING form
#else
    PLOOP(pliter,pl) uu0[pl]=ucons[pl]*ptrgeom->ieomfuncnosing[pl]; // put in UNOTHING form
#endif
  }
  else{
    //    PLOOP(pliter,pl) uu0[pl] = GLOBALMACP0A1(uu0old,i,j,k,pl); // USE OF GLOBALS from phys.tools.rad.c // already in UNOTHING form
  }
  

  
  FTYPE uunew1[NPR];
  PLOOP(pliter,pl) uunew1[pl]=uu[pl];
  //
  FTYPE dugas[4];
  int jj;
  DLOOPA(jj) uunew1[UU+jj] = uu[UU+jj];
  DLOOPA(jj) dugas[jj] = uu[UU+jj]-uu0[UU+jj];
  DLOOPA(jj) uunew1[URAD0+jj] = uu0[URAD0+jj] - dugas[jj];


  
#if(0)
  FTYPE uunew2[NPR];
  PLOOP(pliter,pl) uunew2[pl]=uu[pl];
  //
  FTYPE durad[4];
  DLOOPA(jj) uunew2[URAD0+jj] = uu[URAD0+jj];
  DLOOPA(jj) durad[jj] = uu[URAD0+jj]-uu0[URAD0+jj];
  DLOOPA(jj) uunew2[UU+jj] = uu0[UU+jj] - durad[jj];
  //
  // choose which or how much of each uunew and uunew2
  FTYPE uunewfinal[NPR];
  PLOOP(pliter,pl) uunewfinal[pl] = uu[pl];
  //
  DLOOPA(jj){
    uunewfinal[UU+jj]   = 0.5*(uunew1[UU+jj] + uunew2[UU+jj]);
    uunewfinal[URAD0+jj] = 0.5*(uunew1[URAD0+jj] + uunew2[URAD0+jj]);
  }
#endif

  // now invert just radiation
  int showmessages=0;
  int allowlocalfailurefixandnoreport=1;
  int whichcap=CAPTYPEBASIC; // finalcheck type

  int didrad=0;
  if(-uunew1[URAD0]<-uu[URAD0]){ // only modify radiation if forcing less radiation, not more as that can run-away
    FTYPE GAMMAMAXRADIMPLICITSOLVER=GAMMAMAXRAD;
    PFTYPE lpflag=0;
    PFTYPE lpflagrad=0;
    FTYPE pprad[NPR];
    PLOOP(pliter,pl) pprad[pl]=pp[pl];
    // u2p_rad takes UNOTHING form
    u2p_rad(showmessages, allowlocalfailurefixandnoreport,GAMMAMAXRADIMPLICITSOLVER,whichcap,uunew1,pprad,ptrgeom,&lpflag,&lpflagrad);
    int radinvmod=(int)(lpflagrad);
    if(RADINVOK(radinvmod) && isfinite(pprad[URAD0])==1 && isfinite(pprad[URAD1])==1 && isfinite(pprad[URAD2])==1 && isfinite(pprad[URAD3])==1 && pprad[URAD0] < pp[URAD0]){
      PLOOP(pliter,pl) if(RADFULLPL(pl)) pp[pl] = pprad[pl];
      didrad=1;
    }
    else{
      didrad=0;
      //    dualfprintf(fail_file,"issue ijk=%d %d %d radinvmod=%d\n",i,j,k,radinvmod);
      //    PLOOP(pliter,pl) if(RADFULLPL(pl)) dualfprintf(fail_file,"pprad[%d]=%g\n",pl,pprad[pl]);
    }
  }


#if(0)
  if(0&&didrad==0){
    int finalstep=finaluu;
    int eomtypelocal=EOMTYPE;
    int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
    int modprim=0;
    int checkoninversiongas=CHECKONINVERSION;
    int checkoninversionrad=CHECKONINVERSIONRAD;
    FTYPE dissmeasure=0.0;
  

    struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    // KORALNOTE: If make newtonstats.tryconv~convabs, then if convabs~1E-12, then MHD inversion may return error~1E-10 in terms of how measured with f_error_check(), so must try harder than expected.
    newtonstats.tryconv=1E-3;
    newtonstats.tryconvultrarel=1E-5;
    newtonstats.extra_newt_iter=1;
    newtonstats.extra_newt_iter_ultrarel=1;
#define MINTRYCONVFORMHDINVERSION (1E-4) // assume not failure if got down to this much. -- don't have to be related to implicit allowance.

    newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION;
    newtonstats.maxiter=100;

    FTYPE pptry[NPR];
    PLOOP(pliter,pl) pptry[pl] = pf[pl];
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtypelocal,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UNOTHING,uunewfinal, &qp, ptrgeom, dissmeasure, pptry, pptry, &newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);

    PFTYPE *lpflag,*lpflagrad;
    lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
    lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);
  

    if(IFUTOPRIMFAILSOFT(*lpflag)){
      // PLOOP(pliter,pl) pf[pl] = pptry[pl]; 
    }
    else if(IFUTOPRIMRADFAIL(*lpflagrad)){
    }
    else if( IFUTOPRIMFAIL(*lpflag) || IFUTOPRIMRADFAIL(*lpflagrad) ){
    }
    else{
      PLOOP(pliter,pl) pf[pl] = pptry[pl]; 
    }
  }  
#endif

  struct of_state qf;
  get_state(pf,ptrgeom,&qf);
  FTYPE uf[NPR],ufabs[NPR];
  primtoU(uutype,pf,&qf,ptrgeom,uf,ufabs); // should be closer to having uf[UU]+uf[URAD0] = ucons[UU]+ucons[URAD0]

  FTYPE uudiag[NPR],ufdiag[NPR];
  UtoU(UNOTHING,UDIAG,ptrgeom,uu,uudiag);  // convert from UNOTHING -> UDIAG
  UtoU(UNOTHING,UDIAG,ptrgeom,uf,ufdiag); // convert from UNOTHING -> UDIAG


  // Get deltaUavg[] and also modify ucons if required and should
  int whocalled=COUNTUCONSFIXUP;
  int docorrectuconslocal=0;
  int finalstep=finaluu;
  diag_fixup_dUandaccount(uudiag, ufdiag, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);




  // so exact energy conservation at cost of unknown effect on radiation
  // so one-step accounting below will include a floor for both gas and radiation, but sum of floors will always be zero net gain of energy-momentum as long as radiation has radinvmod==0

  return(0);
}



/// single call in step_ch.c:post_advance() to do all diag_fixup() diagnostic dU stores.  Still allows counts by other diag_fixup calls.
/// for DOONESTEPDUACCOUNTING==1
/// ucons is in UEVOLVE form (originates from unewglobal in step_ch_full() called in main.c)
/// only makes sense if do this with truestep==1 and finalstep==1
/// uses diag_fixup(modcons=-1,...) argument to ensure all other diag calls do not get counted in failfloordu accounting.
int diag_fixup_allzones(FTYPE (*pf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{

  int i, j, k, pliter,pl;
  struct of_geom *ptrgeom;
  struct of_geom geomdontuse;

  if(DOENOFLUX!=NOENOFLUX){
    dualfprintf(fail_file,"Cannot use diag_fixup_allzones() with DOENOFLUX!=NOENOFLUX\n");
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

    int docorrectucons=(DOENOFLUX != NOENOFLUX); // make any needed corrections if doing corrections
    int finalstep=1; // if here, always on finalstep=1
    FTYPE Uf[NPR]; // for returning back pf->Uf result so don't have to repeat if needed later
    diag_fixup_Ui_pf(docorrectucons,MAC(ucons,i,j,k),MAC(pf,i,j,k),ptrgeom,finalstep,COUNTONESTEP, Uf); // Uf in UDIAG form

    FTYPE *pr=MAC(pf,i,j,k);

    int map[NPR];
    FTYPE uconmap[NPR];
    PLOOP(pliter,pl){
      uconmap[pl]=1.0; // default

      if(pl==RHO) map[pl]=YFL1;
      else if(pl==UU) map[pl]=YFL2;
      else if(pl==U3) map[pl]=YFL3;
      else if(pl==URAD0) map[pl]=YFL4;
      else if(pl==URAD3) map[pl]=YFL5;
      else map[pl]=VARNOTDEFINED;
    }


    FTYPE ucon[NDIM],others[NUMOTHERSTATERESULTS];
    ucon_calc(pr,ptrgeom,ucon,others);
    if(YFL1>=0) uconmap[YFL1]=ucon[TT];
    if(YFL2>=0) uconmap[YFL2]=-ucon[TT]; // NOTEMARK: Same sign as in other places like utoprimgen.c
    if(YFL3>=0) uconmap[YFL3]=ucon[TT];
    
    if(YFL4>=0 || YFL5>=0){
      FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];
      ucon_calc(&pr[URAD1-U1],ptrgeom,uradcon,othersrad);
      if(YFL4>=0) uconmap[YFL4]=-uradcon[TT]; // NOTEMARK: Same sign as in other places like utoprimgen.c
      if(YFL5>=0) uconmap[YFL5]=uradcon[TT];
    }

    // assume fixed floor for YFLx like at t=0
    FTYPE offset[NPR];

    FTYPE rhofloor=pr[RHO]*NUMEPSILON*10.0;
    FTYPE vfloor=NUMEPSILON*10.0;
    FTYPE enfloor=ERADLIMIT + (pr[URAD0])*NUMEPSILON*10.0;
    PLOOP(pliter,pl) offset[pl]=SMALL;
    if(YFL1>=0) offset[YFL1] = SMALL + rhofloor; // rho floor
    if(YFL2>=0) offset[YFL2] = SMALL + rhofloor*vfloor*vfloor + UULIMIT; // -T^t_t-rho u^r floor
    if(YFL3>=0) offset[YFL3] = SMALL + rhofloor*vfloor; // T^t_phi floor
    if(YFL4>=0) offset[YFL4] = SMALL + enfloor; // -R^t_t floor
    if(YFL5>=0) offset[YFL5] = SMALL + enfloor*vfloor; // R^t_\phi floor
  

    FTYPE uconsnothing[NPR];
    UtoU(UEVOLVE,UNOTHING,ptrgeom,MAC(ucons,i,j,k),uconsnothing);

    FTYPE Ufnothing[NPR];
    UtoU(UDIAG,UNOTHING,ptrgeom,Uf,Ufnothing);

    int mapvar;
    PLOOP(pliter,pl){
      mapvar=map[pl];
      if(mapvar>=0){
      
        // effective source term for floor scalar accounting for all changes to pl not accounted for by conserved quantity additions (that includes no floors or failures or anything else except fluxes)
        // Assumes floor and mass have same velocity
        // dpl can actually be negative or positive, but only negative when failure
        //      FTYPE pl=MACP0A1(pf,i,j,k,PL); // new final total
        FTYPE pltotal=Ufnothing[pl]/uconmap[mapvar];
        //FTYPE dpl=pl - uconsnothing[pl]/uconmap[mapvar];
        FTYPE dpl=(Ufnothing[pl] - uconsnothing[pl])/uconmap[mapvar];

        if(YFLADDRADTOGAS){
          int plalt;
          if(pl==URAD0){
            plalt=UU;
            dpl+=(Ufnothing[plalt] - uconsnothing[plalt])/uconmap[map[plalt]];
          }
          if(pl==URAD3){
            plalt=U3;
            dpl+=(Ufnothing[plalt] - uconsnothing[plalt])/uconmap[map[plalt]];
          }
        }

        FTYPE plfl;
        //      FTYPE plfl = MACP0A1(pf,i,j,k,mapvar)*pltotal); // final plfl without source term (if failure, then yflx and pl come from averaging, then plfl and pltotal will not be related by conserved fluxes and (e.g.) yflx can become >1
        //        FTYPE yflx = uconsnothing[mapvar]/uconsnothing[pl]; // yflx expected from ucons that accounts for fluxes but no sources, but uconsnothing[pl] could be negative or zero and would have led to failure
        plfl = uconsnothing[mapvar]/uconmap[mapvar]; // plfl from ucons, not consistent with final yflx,pltotal if failure, but averaged yflx probably worse than yflx linked to pltotal
        //      if(plfl<SMALL) plfl=SMALL;

        FTYPE plflfinal = plfl + dpl;
        if(1&&(mapvar==YFL2 || mapvar==YFL4)){ // just energy densities.  Can't apply to angular momenta that can be any sign.  Could apply to density, but doesn't seem to need it.
          // at least for non-densities, especially energy densities, having near 0 or negative values leads to an instability and crazy run-away in the values due to fluxes.
          // So won't be able to track losses of energy, only gains, unless split gains and losses.
          // Or maybe need floor at t=0 at least so that not dealing with crazy small values?
          FTYPE offsetfull=0.0;
          if(pl==UU&&0){ // need to compare to zero when rest-mass added
            offsetfull=offset[mapvar]-uconsnothing[RHO]/ucon[TT];
          }
          else offsetfull=offset[mapvar];
          if(plflfinal<offsetfull) plflfinal=offsetfull; // allow flux and source to compensate to give final >0 quantity
        }

#if(0)
        // plfl could have been negative, as if failure or other process pulled away total rest-mass.  Since force not have floor on value, then that means plfl can go greater than pltotal sometimes.  So avoid.
        if(plflfinal>pltotal) plflfinal=pltotal; // limit so effective true Y_fl<=1
#endif

        // set actual total change in effective floor
        if(DOYFL==1){ // here get true Y_fl=plflfinal/pltotal
          MACP0A1(pf,i,j,k,mapvar) = plflfinal/(SMALL+fabs(pltotal)); // newYflx = newplfl / newpltotal
        }
        else if(DOYFL==2){ // source term to density scalar, assuming pf was evolved via fluxes already for YFLx itself, so here we just add error term from normal evolved quantities
          //          if(mapvar==YFL4 && plflfinal<ERADLIMIT) plflfinal=ERADLIMIT;
          MACP0A1(pf,i,j,k,mapvar) = plflfinal;
        }
        //      dualfprintf(fail_file,"pltotal=%g dpl=%g plfl=%g plflfinal=%g yfl=%g\n",pltotal,dpl,plfl,plflfinal,MACP0A1(pf,i,j,k,YFL));


      } // end if doing this Yflx
    }// end loop over Yflx
  }// end spatial loop


  return(0);

}



/// account for changes by tracking conserved quantities
/// accounts for both failures and floor recoveries
/// this modifies unew if on finalstep to be consistent with floor-limited primitive
/// diagnostics only for actions on conservative quantities
/// assume COUNT types are of PFTYPE
int diag_fixup(int docorrectucons, FTYPE *pr0, FTYPE *pr, FTYPE *uconsinput, struct of_geom *ptrgeom, int finalstep, int doingmhdfixup, int whocalled)
{
  struct of_state q;
  FTYPE Uicent[NPR],Ufcent[NPR];
  int failreturn;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  FTYPE deltaUavg[NPR],Uiavg[NPR];
  FTYPE Uprefixup[NPR],Upostfixup[NPR];
  int docorrectuconslocal;
  int pliter,pl;

  // store ucons and only change if needed (and handle avoidance of B1,B2,B3 in case ucons points to unewglobal staggered field, while here we operate on centers)
  FTYPE ucons[NPR];
  if(uconsinput!=NULL){
    PLOOP(pliter,pl) ucons[pl]=uconsinput[pl];
  }

#if(DOSUPERDEBUG)
  superdebug_utoprim(pr0,pr,ptrgeom,whocalled);
  // collect values for non-failed and failed zones
#endif



  // count whocalled diag_fixup()
  if(doingmhdfixup) count_whocalled(ptrgeom->i,ptrgeom->j,ptrgeom->k, finalstep, whocalled,1);



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
    failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Uicent, NULL);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
 

    // after any changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Ufcent, NULL);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");

    // if Uicent and Ufcent are both from pi and pf at CENT, then B1,B2,B3 entries are agreeably located even for FLUXB==FLUXCTSTAG

    // Get deltaUavg[] and also modify ucons if required and should
    diag_fixup_dUandaccount(Uicent, Ufcent, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);


    if(uconsinput!=NULL && docorrectuconslocal){
      
      // copy over ucons result in case changed.
      PLOOP(pliter,pl){
        // if staggered field, avoid modifying field since at FACEs, not CENT where pr lives.
        if(BPL(pl)==0 && FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){
          uconsinput[pl] = ucons[pl];
        }
      }
    }


  }// end if finalstep>0



  return(0);
}



/// like diag_fixup(), but input initial conserved quantity as Ui and final primitive as pf
/// Must use this when pi[Ui] doesn't exist and had to use non-hot-MHD inversion.
/// Assumes Ui is like unewglobal, so UEVOLVE type
/// Assume ultimately hot MHD equations are used, so need to get new Uf that'll differ from Ui
/// Also don't know Uf quite yet.
int diag_fixup_Ui_pf(int docorrectucons, FTYPE *Uievolve, FTYPE *pf, struct of_geom *ptrgeom, int finalstep, int whocalled, FTYPE *Uf)
{
  struct of_state q;
  FTYPE Ufcent[NPR],Uicent[NPR],uconsinput[NPR],ucons[NPR];
  int failreturn;
  int pliter,pl,enerregion;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int docorrectuconslocal;





  // count whocalled diag_fixup()
  count_whocalled(ptrgeom->i,ptrgeom->j,ptrgeom->k, finalstep, whocalled,1);


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
    failreturn=primtoU(UDIAG,pf,&q,ptrgeom,Ufcent, NULL);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");

    if(Uf!=NULL) PLOOP(pliter,pl) Uf[pl]=Ufcent[pl];// store result for returning before gets modified

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

    PFTYPE *lpflag,*lpflagrad;
    lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
    lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);

    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"lpflag=%d lpflagrad=%d\n",*lpflag,*lpflagrad);
    //    }

    if(*lpflag>UTOPRIMNOFAIL){
      // then assume fixup_utoprim() needs to operate and will also handle accounting
      PLOOP(pliter,pl) if(RADPL(pl)==0) Ufcent[pl] = Uicent[pl];
    }
    if(IFUTOPRIMRADFAIL(*lpflagrad)){
      // then assume fixup_utoprim() needs to operate and will also handle accounting
      PLOOP(pliter,pl) if(RADPL(pl)==1) Ufcent[pl] = Uicent[pl];
    }



    ////////////////
    //
    // Get deltaUavg[] and also modify ucons if required and should
    //
    ////////////////
    diag_fixup_dUandaccount(Uicent, Ufcent, ucons, ptrgeom, finalstep, whocalled, docorrectuconslocal);


    if(docorrectuconslocal){
      // copy over ucons result in case changed.
      PLOOP(pliter,pl){
        // if staggered field, avoid modifying field since at FACEs, not CENT where pr lives.
        if(BPL(pl)==0 && FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){
          uconsinput[pl] = ucons[pl];
        }
      }
    }


  }// end if finalstep>0



  return(0);
}




/// account for changes by tracking conserved quantities
/// accounts for both failures and floor recoveries
/// only called on final step of RK once unew is defined since only on final step is unew modified if floor encountered
/// ONLY used by phys.ffde.c inversion routine when E^2>B^2
/// Assume Ui and Uf in UDIAG form
int diag_fixup_U(int docorrectucons, FTYPE *Ui, FTYPE *Uf, FTYPE *uconsinput, struct of_geom *ptrgeom, int finalstep,int whocalled)
{
  FTYPE Uicent[NPR],Ufcent[NPR];
  struct of_state q;
  int failreturn;
  int pliter,pl,enerregion, tscale;
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int docorrectuconslocal;


  // store ucons and only change if needed (and handle avoidance of B1,B2,B3 in case ucons points to unewglobal staggered field, while here we operate on centers)
  FTYPE ucons[NPR];
  if(uconsinput!=NULL){
    PLOOP(pliter,pl) ucons[pl]=uconsinput[pl];
  }

  // count whocalled diag_fixup()
  count_whocalled(ptrgeom->i,ptrgeom->j,ptrgeom->k, finalstep, whocalled,1);
  


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

    if(uconsinput!=NULL && docorrectuconslocal){
      // copy over ucons result in case changed.
      PLOOP(pliter,pl){
        // if staggered field, avoid modifying field since at FACEs, not CENT where pr lives.
        if(BPL(pl)==0 && FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){
          uconsinput[pl] = ucons[pl];
        }
      }
    }


  }




  return(0);
}




#if(VARTOINTERP==PRIMTOINTERP_GDETFULLVERSION_WALD)
#define FIXUPTYPE 0 // or else would create spurious Poynting flux
#else
#define FIXUPTYPE 0 // too expensive if inversion fails alot as can happen near floors, near poles, or with radiation.
#endif
/// 0 = primitive (adds rho,u in comoving frame)
/// 1 = conserved but rho,u added in ZAMO frame
/// 2 = conserved but ignore strict rho,u change for ZAMO frame and instead conserved momentum (doesn't keep desired u/rho, b^2/rho, or b^2/u and so that itself can cause problems
///
/// finalstep==0 is non-accounting, finalstep==1 is accounting
int fixup1zone(int docorrectucons, FTYPE *pr, FTYPE *uconsinput, struct of_geom *ptrgeom, int finalstep)
{
  int pliter,pl;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE ftempA,ftempB;
  struct of_state q;
  struct of_state dq;
  FTYPE prfloor[NPR],prceiling[NPR];
  FTYPE prdiag[NPR];
  FTYPE pr0[NPR];
  FTYPE prmhdnew[NPR];
  FTYPE U[NPR],U0[NPR];
  int checkfl[NPR];
  int failreturn;
  int didchangeprim;
  FTYPE scalemin[NPR];
  //  FTYPE ucovzamo[NDIM];
  //  FTYPE uconzamo[NDIM];
  FTYPE dprmhd[NPR];
  FTYPE prmhd[NPR];
  FTYPE dU[NPR];
  //  FTYPE P,Pnew;
  int jj;
  int badinversion;
  PFTYPE oldmhdpflag,oldradpflag;



  // store ucons and only change if needed (and handle avoidance of B1,B2,B3 in case ucons points to unewglobal staggered field, while here we operate on centers)
  FTYPE ucons[NPR];
  if(uconsinput!=NULL){
    PLOOP(pliter,pl) ucons[pl]=uconsinput[pl];
  }


  // assign general floor variables
  // whether to check floor condition
  PALLLOOP(pl){
    checkfl[pl]=0;
    pr0[pl]=pr[pl];
    prmhd[pl]=pr0[pl];
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
  // KORALTODO: Keep as mhd only for now.
  int DOEVOLVEURAD0=PRAD0>=0 && EOMRADTYPE!=EOMRADNONE;
  if(DOEVOLVEURAD0) checkfl[PRAD0]=1;


  ////////////////////
  //
  // limit gamma wrt normal observer -- this can change b^2, so do this first.
  //
  ////////////////////
#if(WHICHVEL==VELREL4)
  //  int docorrectucons=(DOENOFLUX != NOENOFLUX);
  //  didchangeprim=0;

  //  if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
  //    dualfprintf(fail_file,"BEFORE IN FIXUP1ZONE LIMITGAMMA: finalstep=%d\n",finalstep);
  //  }
  failreturn=limit_gamma(docorrectucons,GAMMAMAX,GAMMAMAXRAD,prmhd,ucons,ptrgeom,-1);
  if(failreturn>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
  if(failreturn==-1) didchangeprim=1;
  //  PALLLOOP(pl) prdiag[pl]=prmhd[pl];
  //  if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
  //    dualfprintf(fail_file,"AFTER IN FIXUP1ZONE LIMITGAMMA: didchangeprim=%d\n",didchangeprim);
  //  }

  //  if(didchangeprim&&FLOORDIAGS){// FLOORDIAGS includes fail diags
  //    int doingmhdfixup=1;
  //    diag_fixup(docorrectucons,prdiag, prmhd, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTLIMITGAMMAACT);
  //    PALLLOOP(pl) prdiag[pl]=prmhd[pl];
  //  }

#endif// end if WHICHVEL==VEL4REL


    
  ////////////
  //
  // Only apply floor if cold or hot GRMHD
  //
  ////////////
  if(DOEVOLVERHO||DOEVOLVEUU||DOEVOLVEURAD0){


    //////////////
    //
    // get floor value (computes, e.g., bsq)
    //
    //////////////
    set_density_floors(ptrgeom,prmhd,prfloor,prceiling);
    scalemin[RHO]=RHOMINLIMIT;
    scalemin[UU]=UUMINLIMIT;
    if(URAD0>=0) scalemin[URAD0]=ERADLIMIT;
    
    

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
    // Get new primitive if went beyond floor
    //
    /////////////////////////////

    /// default
    PALLLOOP(pl){
      prmhdnew[pl]=prmhd[pl];
    }

    PALLLOOP(pl){
      if ( checkfl[pl]&&(prmhd[pl] < prfloor[pl] && prmhd[pl]<prceiling[pl]) ){
        didchangeprim=1;
        //dualfprintf(fail_file,"%d : %d %d %d : %d : %d : %21.15g - %21.15g\n",pl,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,checkfl[pl],prfloor[pl],prmhd[pl]); 
        // only add on full step since middle step is not really updating primitive variables
        prmhdnew[pl]=prfloor[pl];
      }
    }

    PALLLOOP(pl){
      if ( checkfl[pl]&&(prmhd[pl] > prfloor[pl] && prmhd[pl]>prceiling[pl]) ){
        didchangeprim=1;
        //dualfprintf(fail_file,"%d : %d %d %d : %d : %d : %21.15g - %21.15g\n",pl,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,checkfl[pl],prfloor[pl],prmhd[pl]); 
        // only add on full step since middle step is not really updating primitive variables
        prmhdnew[pl]=prceiling[pl];
      }
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
        prmhd[pl]=prmhdnew[pl];
      }


#elif(FIXUPTYPE==1 || FIXUPTYPE==2)
      // mass and internal energy added in frame not necessarily the comoving frame
      // using a frame not directly associated with comoving frame avoids arbitrary energy-momentum  growth
      // GODMARK: FIXUPTYPE==1 doesn't exactly match between when b^2/\rho_0>BSQORHOULIMIT such that amount of mass added will force equality of b^2/\rho_0==BSQORHOLIMIT, so this may lead to problems.
      // physically FIXUPTYPE==1 models some non-local transport of baryons and energy to a location that supposedly occurs when b^2/rho_0 is too large.
      // FIXUPTYPE==2 models a local injection of baryons/energy with momentum conserved and energy-momentum conserved if mass injected.  Injection will slow flow.  Essentially there is an ad hoc conversion of kinetic/thermal energy into mass energy.

      // compute original conserved quantities
      failreturn=get_state(prmhd,ptrgeom,&q);
      if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
      failreturn=primtoU(UNOTHING,prmhd,&q,ptrgeom,U, NULL);
      if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");

      // store original U
      PLOOP(pliter,pl) U0[pl]=U[pl];

      // get change in primitive quantities
      PALLLOOP(pl) dprmhd[pl]=0.0; // default
      // use ZAMO velocity as velocity of inserted fluid
      PALLLOOP(pl) if(pl==RHO || pl==UU || pl==URAD0){ dprmhd[pl]=prmhdnew[pl]-prmhd[pl];}
      set_zamo_velocity(WHICHVEL,ptrgeom,dprmhd);

      // get change in conserved quantities
      failreturn=get_state(dprmhd,ptrgeom,&dq);
      failreturn=primtoU(UNOTHING,dprmhd,&dq,ptrgeom,dU, NULL);
      if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");


      if(FIXUPTYPE==1){
        // then done, dU is right
      }
      else if(FIXUPTYPE==2){
        // then don't allow momentum to change regardless of meaning for implied rho,u
        dU[U1]=dU[U2]=dU[U3]=0.0;
        if(URAD0>=0) dU[URAD0]=dU[URAD1]=dU[URAD2]=dU[URAD3]=0.0;

        pl=UU;
        if ( checkfl[pl]&&(prfloor[pl] > prmhd[pl] || prceiling[pl] < prmhd[pl]) ){
          // then must change dU[UU]
        }
        else dU[UU]=0.0; // if only mass added, then no change needed to energy-momentum


      }

      // get final new conserved quantity
      PALLLOOP(pl) U[pl]+=dU[pl];
      // except, if fixed-up u_g or rho because <0, then just set entropy.
      // If u_g>0 and rho>0, then assume entropy adjustment also ok.
      //      if((prmhd[UU]<=0.0 || prmhd[RHO]<=0.0) && ENTROPY>0) U[ENTROPY] = U0[ENTROPY];
      // assume this adjustment is energy-only based.
      // must do this because otherwise if u_g<0 or rho<0, then entropy is ill-defined, and here assume floor always activated related to too small u_g or rho so that entropy would be badly defined or ill-defined.
      if(ENTROPY>=0) U[ENTROPY] = U0[ENTROPY];


      // pr finally changes here
      // get primitive associated with new conserved quantities
      oldmhdpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
      oldradpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);

      int showmessages=0; // messages not important if fixup doens't work, unless debugging.
      int allowlocalfailurefixandnoreport=1; 
      int eomtype=EOMDEFAULT;
      FTYPE dissmeasure=-1.0; // assume energy try ok
      int whichcap=CAPTYPEBASIC;
      int whichmethod=MODEDEFAULT;
      int modprim=0;
      int checkoninversiongas=CHECKONINVERSION;
      int checkoninversionrad=CHECKONINVERSIONRAD;

      struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
      if(0){
        newtonstats.nstroke=newtonstats.lntries=0;
      }
      else{
        ////////
        // be quick about checking this
        newtonstats.nstroke=newtonstats.lntries=0;
        // set inputs for errors, maxiters, etc.
#define IMPTRYCONVCONSTFORFX (1E-6)
#define IMPOKCONVCONSTFORFX (1E-4)
#define IMPALLOWCONVCONSTFORFX (1E-3)
#define IMPMAXITERFORFX (10)
        newtonstats.tryconv=1E-2*MAX(IMPTRYCONVCONSTFORFX,IMPOKCONVCONSTFORFX);
        newtonstats.tryconvultrarel=1E-2*MAX(IMPTRYCONVCONSTFORFX,IMPOKCONVCONSTFORFX);
        newtonstats.extra_newt_iter=1;
        newtonstats.extra_newt_iter_ultrarel=2;
        newtonstats.mintryconv=IMPALLOWCONVCONSTFORFX;
        newtonstats.maxiter=IMPMAXITERFORFX;
        //
      }
      //      dualfprintf(fail_file,"BEFORE FIXUPUTOPRIMGEN\n");
      failreturn=Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,OTHERUTOPRIM,UNOTHING,U,&q, ptrgeom,dissmeasure,prmhd,prmhd,&newtonstats);
      //      dualfprintf(fail_file,"AFTER FIXUPUTOPRIMGEN\n");
      // have to add since takes effort.s
      nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;

      // KORALNOTEMARK: Only changing floor related to MHD fluid so far, so no check on failure of radiation inversion.
      badinversion = (failreturn>=1 || IFUTOPRIMFAIL(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)));

      static long long int utoprimgenfixup=0,utoprimgenfixupbad=0;
      utoprimgenfixup++;
      if(badinversion) utoprimgenfixupbad++;
      if(debugfail>=2) if(utoprimgenfixup%totalzones==0) dualfprintf(fail_file,"UTOPRIMGENFIXUP: %lld (bad=%lld) : %ld %d\n",utoprimgenfixup,utoprimgenfixupbad,nstep,steppart);


      if(badinversion){
        if(debugfail>=2) dualfprintf(fail_file,"Utoprimgen failed in fixup.c");
        // if problem with Utoprim, then just modify primitive quantities as normal without any special constraints
        PALLLOOP(pl){
          prmhd[pl]=prmhdnew[pl];
        }
      }

      // in any case, this is not normal inversion procedure, so clear the flags
      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=oldmhdpflag;
      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=oldradpflag;

#endif


      // get min or max of comoving and ZAMO frame versions of densities
      // generally trying to keep atmosphere non-relativistic
      //      pl=RHO;   prmhd[pl]=MIN(prmhd[pl],prmhdnew[pl]);
      // get larger of rho
      pl=RHO;   prmhd[pl]=MAX(prmhd[pl],prmhdnew[pl]);
      // get smallest of two
      pl=UU;    prmhd[pl]=MIN(prmhd[pl],prmhdnew[pl]);
      pl=URAD0; prmhd[pl]=MIN(prmhd[pl],prmhdnew[pl]);


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
    int docorrectucons2=(DOENOFLUX != NOENOFLUX);
    int doingmhdfixup=1;
    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"FLOORACTBEFORE: steppart=%d nstep=%ld\n",steppart,nstep);
    //    }

    PFTYPE *lpflag,*lpflagrad;
    lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
    lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);

    if(*lpflag>UTOPRIMNOFAIL || IFUTOPRIMRADFAIL(*lpflagrad)){
      // by this point, pr0 or prdiag are not necessarily consistent with ucons, which if finalstep=1 is the final conserved quantity.  This occurs when no implicit solution and no explicit inversion.  Happens when (e.g.) U[RHO]<0.
      FTYPE Uievolve[NPR];
      PLOOP(pliter,pl) Uievolve[pl] = ucons[pl];
      diag_fixup_Ui_pf(docorrectucons2, Uievolve, prmhd, ptrgeom, finalstep,COUNTFLOORACT,NULL);
    }
    else{
      // if no failure, then fixup_utoprim won't be called or do diagnostics, so floor will be preserved and so do accounting here.
      diag_fixup(docorrectucons2,prdiag, prmhd, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTFLOORACT);
    }

    


    // now prdiag=prmhd as far as diag_fixup() is concerned, so next changes are new changes (i.e. don't cumulative w.r.t. pr0 multiple times, since that (relative to prmhd each time) would add each prior change to each next change)
    //    PALLLOOP(pl) prdiag[pl]=prmhd[pl];
    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"FLOORACTAFTER: steppart=%d nstep=%ld\n",steppart,nstep);
    //    }
  }






    // assume once we go below floor, all hell will break loose unless we calm the storm by shutting down this zone's relative velocity
    // normal observer velocity
    // i.e. consider this a failure
    //GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)= 1;


  //////////////////////////
  //
  // now keep track of modified primitives via conserved quantities
  //
  //////////////////////////
  if(didchangeprim){

    ///////////////
    //
    // Assign final primitives
    //
    ////////////////

    PALLLOOP(pl) pr[pl]=prmhd[pl];

    // SUPERGODMARK: only doing rest-mass density here for now.  Still account for YFL2-5 via finalstep==1 source injection, but sub-steps only then handle fluxes for YFL2-5 currently unless add how conserved quantities changed here.  In any case, this is estimate for finalstep actual ucons vs. Uf(pf) changes.
    if(DOYFL){
      FTYPE drho=(pr[RHO]-pr0[RHO]);
      if(DOYFL==1){
        pr[YFL1] = (pr0[YFL1]*pr0[RHO] + drho)/pr[RHO];// have floor set new floor scalar fraction.  Adds mass at same velocity for FIXUPTYPE==0
        if(pr[YFL1]<0.0) pr[YFL1]=0.0;
        if(pr[YFL1]>1.0) pr[YFL1]=1.0;
      }
      else if(DOYFL==2){
        pr[YFL1] = pr0[YFL1] + drho;
        if(pr[YFL1]<SMALL) pr[YFL1]=SMALL;
        if(pr[YFL1]>pr[RHO]) pr[YFL1]=pr[RHO];
      }


      if(DOYFL==2&&0){
        // also iterate on other YFLx's
        struct of_state qnew;
        get_state(pr,ptrgeom,&qnew);
        FTYPE Unew[NPR];
        primtoU(UNOTHING,pr,&qnew,ptrgeom,Unew, NULL);

        struct of_state q0;
        get_state(pr0,ptrgeom,&q0);
        FTYPE U0[NPR];
        primtoU(UNOTHING,pr0,&q0,ptrgeom,U0, NULL);

        if(YFL2>=0) pr[YFL2]+= -(Unew[UU]-U0[UU])/qnew.ucon[TT];
        if(YFL3>=0) pr[YFL3]+= +(Unew[U3]-U0[U3])/qnew.ucon[TT];
        if(YFLADDRADTOGAS){
          if(YFL2>=0) pr[YFL2]+= -(Unew[URAD0]-U0[URAD0])/qnew.uradcon[TT];
          if(YFL3>=0) pr[YFL3]+= +(Unew[URAD3]-U0[URAD3])/qnew.uradcon[TT];
        }

        if(YFL4>=0){
          pr[YFL4]+= -(Unew[URAD0]-U0[URAD0])/qnew.uradcon[TT];
          if(pr[YFL4]<ERADLIMIT) pr[YFL4]=ERADLIMIT;
        }
        if(YFL5>=0) pr[YFL5]+= (Unew[URAD3]-U0[URAD3])/qnew.uradcon[TT];
      }
    }


#if(0)
    ///////////////
    //
    // now ensure primitive and conserved consistent for RK3 and other methods that use Uf
    // assumes "ucons" is utoinvert in advance.c so it's uf or utempcum when necessary as associated with inversion (internal implicit or external).
    //
    //////////////
    // here's where we correct ucons
    docorrectucons=1;
    struct of_state qnew;
    get_state(pr,ptrgeom,&qnew);

    primtoU(UEVOLVE,pr,&qnew,ptrgeom,ucons, NULL);
#endif

    return(-1); // -1 means made changes

  }


  if(uconsinput!=NULL && docorrectucons){
    // copy over ucons result in case changed.
    PLOOP(pliter,pl){
      // if staggered field, avoid modifying field since at FACEs, not CENT where pr lives.
      if(BPL(pl)==0 && FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){
        uconsinput[pl] = ucons[pl];
      }
    }
  }


  //  if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
  //    struct of_state qb;
  //    get_state(pr,ptrgeom,&qb);
  //    FTYPE ub[NPR],ubabs[NPR];
  //    int uutype=UDIAG;
  //    primtoU(uutype,pr,&qb,ptrgeom,ub,ubabs);
  //    dualfprintf(fail_file,"URHOINFIXUPU1ZONE=%21.15g\n",ub[RHO]*dx[1]*dx[2]*dx[3]);
  //  }
  


  return(0);
}


// number of vote checks per zone
#define MAXVOTES 8

#define NUMCHECKS 2
#define ISGAMMACHECK 0
#define ISUUCHECK 0

/// GODMARK: function is 2D right now, but works in 3D, it just uses only x-y plane for checking
/// check whether solution seems reasonable
/// useful if b^2/rho\gg 1 or try to approach stationary model where variations in space shouldn't be large zone to zone
/// checks whether u or gamma is much different than surrounding zones.
/// if true, then flag as failure, else reasonable solution and keep
/// can't assume failed zones are reasonably set
/// fixup_checksolution() currently only uses pflag[FLAGUTOPRIMFAIL]
/// KORALNOTEMARK: So far fixup_checksolution() only setup for MHD fluid, but ok since don't use this function anymore (i.e. CHECKSOLUTION set to 0 usually because this function can cause problems).
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
        // stderrfprintf("caught one-0: %d %d\n",i,j);
        GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILGAMMAPERC;
      }
      checki=ISUUCHECK;
      if(checkcondition[checki] && (vote[checki]>numvotes[checki]*0.5)){
        // then majority rules
        // stderrfprintf("caught one-1: %d %d\n",i,j);
        GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILUPERC;
      }

    }// end COMPZLOOP
  }// end parallel region

  return(0);
}













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
// Also, during failure, can't assume mass can be so easily conserved and D might be bad
// So disable for now
// enable for finalstep for YFL, but for RHO causes high gamma fixed in regions, so disable for that
#define DO_CONSERVE_D_INFAILFIXUPS (finalstep==1 && (SCALARPL(pl))) // not pl==RHO

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


/// fix the bad solution as determined by utoprim() and fixup_checksolution()
/// needs fail flag over -NUMPFLAGBNDx .. ..N-1+NUMPFLAGBNDx , and uses p at 0..N-1
int fixup_utoprim(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR];
  extern void get_advance_startendindices(int *is,int *ie,int *js,int *je,int *ks,int *ke);
  int is,ie,js,je,ks,ke;




  // this average only works if using 4 velocity since only then guaranteed solution is good after interpolation
  if(WHICHVEL==VEL3) return(0); // just stick with static, best can do
  if(EOMTYPE==EOMFFDE || EOMTYPE==EOMFFDE2) return(0); // nothing to do




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
    
    // +-NUMPFLAGBNDx extra so can do check
    is += -NUMPFLAGBND1;
    ie += +NUMPFLAGBND1;
    js += -NUMPFLAGBND2;
    je += +NUMPFLAGBND2;
    ks += -NUMPFLAGBND3;
    ke += +NUMPFLAGBND3;

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
  // only copies fail flags
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
    PFTYPE mhdlpflag,radlpflag;
    int fixedmhd,fixedrad;
    int startpl,endpl;
    FTYPE ftemp;
    FTYPE D0;
    //
    int limitedgamma;
    int nogood;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    int fixingmhd,fixingrad;
    int limitgammamhd,limitgammarad;


    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPZLOOP;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      

     
      // get failure flag
      mhdlpflag=GLOBALMACP0A1(pflagfailorig,i,j,k,FLAGUTOPRIMFAIL);
      radlpflag=GLOBALMACP0A1(pflagfailorig,i,j,k,FLAGUTOPRIMRADFAIL);

      // assume not fixing anything
      fixingmhd=IFUTOPRIMFAIL(mhdlpflag);
      fixingrad=IFUTOPRIMRADFAIL(radlpflag); // not just a hard failure, but some softish ones too

      
      
      if(IFUTOPRIMFAILFIXED(mhdlpflag) && IFUTOPRIMFAILFIXED(radlpflag)){ // IF BOTH MHD AND RAD FIXED-UP ELSEWHERE
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
        fixuputoprim_accounting(i, j, k, mhdlpflag, radlpflag, 0, 0, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);


      }
      else if(fixingmhd || fixingrad){





        if(fixingmhd){

          ///////////////////////////////////////////////////////////
          /////////////////
          //
          // see if MHD utoprim() failed
          //
          ///////////////////////////////////////////////////////////
          //////////////////

          fixedmhd=0; // assume not fixed yet

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
            // Doesn't include any floor scalars (e.g. YFLx) presumed to be forced as fully conservative post-adjusted ucon[TT] unlike (normally) evolved density
            //
            //////////////////
            // field is evolved fine, so only average non-field
            if(mhdlpflag==UTOPRIMFAILRHOUNEG && HANDLERHOUNEG==1){
              startpl=RHO;
              endpl=UU;
            }
            else if(mhdlpflag==UTOPRIMFAILU2AVG1 || mhdlpflag==UTOPRIMFAILU2AVG2 || mhdlpflag==UTOPRIMFAILU2AVG1FROMCOLD || mhdlpflag==UTOPRIMFAILU2AVG2FROMCOLD || mhdlpflag==UTOPRIMFAILUPERC || mhdlpflag==UTOPRIMFAILUNEG && (HANDLEUNEG==1) ){
              startpl=UU;
              endpl=UU;
            }
            else if(mhdlpflag==UTOPRIMFAILRHONEG && HANDLERHONEG==1){
              startpl=RHO;
              endpl=RHO;
            }
            else if(mhdlpflag==UTOPRIMFAILURHO2AVG1FROMFFDE){
              startpl=RHO;
              endpl=UU;
              // SUPERGODMARK: should also average v along v.B
            }
            else{
              // then presume inversion failure with no solution or assuming rho<=0 or u<=zerouuperbaryon*prim[RHO] is bad inversion if HANDLE?NEG==0
              startpl=RHO;
              endpl=NPR-1;
            }







            //////////////////////////////
            //
            // other kinds of failures not caught by above (inversion convergence type failures)
            //
            //////////////////////////////
            if(fixedmhd==0){


              /////////////////////
              //
              // fix using average of surrounding good values, if they exist
              //
              //////////////////////
              nogood=0;
              int doingmhd=1;
#if(GENERALAVERAGE==1)
              int useonlynonfailed=1;
              int maxnumbndtotry=MAXNUMPFLAGBND; // choice smaller or equal to MAXNUMPFLAGBND
              int numbndtotry;
              // TODOMARK: could also do multi-pass fixup_utoprim() with BCs done as squeeze in on failure region from outside, in case region is not caught by this.
              for(numbndtotry=1;numbndtotry<=maxnumbndtotry;numbndtotry++){
                nogood=general_average(useonlynonfailed,numbndtotry,maxnumbndtotry,startpl, endpl, i, j, k, doingmhd, mhdlpflag, radlpflag, GLOBALPOINT(pflagfailorig) ,pv,ptoavg,ptrgeom);
                if(nogood==0) break; // then success before using any more points
              }
#else
              nogood=simple_average(startpl,endpl,i,j,k, doingmhd, GLOBALPOINT(pflagfailorig) ,pv,ptoavg);
#endif


              if(nogood){
                //////////////////////////////
                //
                // fixup negative densities (if no good case)
                //
                //////////////////////////////
                fixup_negdensities(EOMSETMHD, &fixedmhd, startpl, endpl, i, j, k, mhdlpflag, pv,ptoavg, ptrgeom, pr0, ucons, finalstep);

                if(fixedmhd==1 && (startpl<=RHO && endpl<=UU)){
                  // assume success for densities, so no longer nogood.
                  nogood=0;
                }

                if(fixedmhd==1 && (startpl<=RHO && endpl>=U1)){
                  // then fixup but only changed densities, so still need to process non-densities
                  startpl=U1; // start at U1 (first velocity) and finish at same ending if was ending on some velocity
                  fixedmhd=0; // then reset fixedmhd->0 so can still modify these remaining quantities -- otherwise v^i would be unchanged even if wanted to average that out for (e.g.) negative density results for inversions.
                }
              }

              /////////////////////
              //
              // If no good surrounding value found, then average bad values in some way
              //
              //////////////////////
              if(nogood){
                //                dualfprintf(fail_file,"fixupnogood: mhd %d %d %d : %d %d\n",i,j,k,startpl,endpl);
                fixup_nogood(startpl, endpl, i, j, k, doingmhd, mhdlpflag, radlpflag, pv,ptoavg, pbackup, ptrgeom);
              }



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
                //   alpha = 1. / sqrt(-ptrgeom->gcon[GIND(0,0)]);
                alpha = ptrgeom->alphalapse;
                vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);

                limitedgamma=0;
                if(gamma>GAMMAMAX){
                  limitedgamma=1;
                  if(debugfail>=2) dualfprintf(fail_file,"MHD: initial gamma: %21.15g,  max: %21.15g, initial vsq: %21.15g\n",gamma,GAMMAMAX,vsq);
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
                //                if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
                //                  dualfprintf(fail_file,"BEFORE IN FIXUPUTOPRIM LIMITGAMMA: finalstep=%d flag=%d\n",finalstep,mhdlpflag);
                //                }
                if(limitgammamhd=limit_gamma(0,gamma,GAMMAMAXRAD,MAC(pv,i,j,k),MAC(ucons,i,j,k),ptrgeom,-1)>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 2);
                //                if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
                //                  dualfprintf(fail_file,"AFTER IN FIXUPUTOPRIM LIMITGAMMA: finalstep=%d\n",finalstep);
                //                }

                if(debugfail>=3){
                  if(limitedgamma){
                    // check gamma
                    MYFUN(gamma_calc(MAC(pv,i,j,k),ptrgeom,&gamma,&qsq),"fixup_utoprim: gamma calc failed\n","fixup.c",2);
                    dualfprintf(fail_file,"MHD: final gamma: %21.15g\n",gamma);
                  }
                }
#endif

              }// end if dealing with velocity



              /////////////////////
              //
              // Things to do only if modifying density (must come after modifying velocity that enters ucon[TT])
              //
              //////////////////////
              PLOOP(pliter,pl){
                if(pl==RHO || pl==YFL1 || pl==YFL2 || pl==YFL3){
                  if(startpl<=pl && endpl>=pl){


                    if(DO_CONSERVE_D_INFAILFIXUPS){
                      
                      if(mhdlpflag==UTOPRIMFAILGAMMAPERC || 1){ // GODMARK: always doing it

                        if (ucon_calc(MAC(ptoavg,i,j,k), ptrgeom, ucon, others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);

                        // Use D0 to constrain how changing u^t changes rho
                        // GODMARK: Why not used evolved D=\rho_0 u^t  from conserved quantity?
                        // GODMARK: See fixup.c's limit_gamma() notes on why using conserved version of D not good to use
                        // Here we ignore all conserved quantities and just ensure that D0 is conserved (close) to original value after averaging that assumes original value was reasonable
                        if(finalstep==1){
                          // this is not reasonable, because ignores any mass injection required during failure to make sense of solution
                          FTYPE Unothing[NPR];
                          UtoU(UEVOLVE,UNOTHING,ptrgeom,MAC(ucons,i,j,k),Unothing);
                          D0 = Unothing[pl] ;
                        }
                        else{
                          // This is probably not necessary or useful
                          D0 = MACP0A1(ptoavg,i,j,k,pl)*ucon[TT];
                        }
                        
                        ///////////////////////////////////////////
                        //
                        // constrain change in density so conserve particle number
                        //
                        //////////////////////////////////////////
                        if(D0>0.0 && pl==RHO || pl==YFL1 || pl==YFL2 || pl==YFL3){ // only makes sense if D0 implies positive density at least
                          MACP0A1(pv,i,j,k,pl) = D0/ucon[TT];
                        }
                      }
                    }
                  }
                }
              }// end over density or scalars



#if(0)
              // DEBUG problem of launch with pressureless stellar model collapse
              PALLLOOP(pl) dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d j=%d k=%d pl=%d pv=%21.15g ptoavg=%21.15g\n",nstep,steppart,i,j,k,pl,MACP0A1(pv,i,j,k,pl),MACP0A1(ptoavg,i,j,k,pl));
#endif



            } // end if fixedmhd==0

          }// end if not keeping static
        }// end if mhd failed






      
        if(fixingrad){ // RAD FAILED

          ///////////////////////////////////////////////////////////
          ////////////////////
          //
          // see if u2p_rad() failed in some way
          //
          ////////////////////
          ///////////////////////////////////////////////////////////

          fixedrad=0; // assume not fixedrad yet

          // set pre-fixedrad primitives
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
            if(radlpflag==UTOPRIMFAILU2AVG1 || radlpflag==UTOPRIMFAILU2AVG2 || radlpflag==UTOPRIMFAILU2AVG1FROMCOLD || radlpflag==UTOPRIMFAILU2AVG2FROMCOLD || radlpflag==UTOPRIMFAILUPERC || radlpflag==UTOPRIMFAILUNEG && (HANDLEUNEG==1) || radlpflag==UTOPRIMRADFAILERFNEG){
              startpl=PRAD0; // use as minimum for PRAD0 and average for PRAD1-PRAD3
              endpl=PRAD3;
            }
            else if(radlpflag==UTOPRIMRADFAILGAMMAHIGH){ // fixing gammarad so no constant high value from hitting ceiling on gammaradmax
              startpl=PRAD1;
              endpl=PRAD3;
            }
            else{
              // then presume inversion failure with no solution
              startpl=RHO;
              endpl=NPR-1;
            }





            //////////////////////////////
            //
            // other kinds of failures not caught by above (inversion convergence type failures)
            //
            //////////////////////////////
            if(fixedrad==0){


              /////////////////////
              //
              // fix using average of surrounding good values, if they exist
              //
              //////////////////////
              nogood=0;
              int doingmhd=0;
#if(GENERALAVERAGE==1)
              int useonlynonfailed=1;
              int maxnumbndtotry=MAXNUMPFLAGBND; // choice smaller or equal to MAXNUMPFLAGBND
              int numbndtotry;
              for(numbndtotry=1;numbndtotry<=maxnumbndtotry;numbndtotry++){
                nogood=general_average(useonlynonfailed,numbndtotry,maxnumbndtotry,startpl, endpl, i, j, k, doingmhd, mhdlpflag, radlpflag, GLOBALPOINT(pflagfailorig) ,pv,ptoavg,ptrgeom);
                if(nogood==0) break; // then success before using any more points
              }
#else
              nogood=simple_average(startpl,endpl,i,j,k, doingmhd, GLOBALPOINT(pflagfailorig) ,pv,ptoavg);
#endif


              if(nogood && 0){ // avoid fixup_negdensities() because resets to ERADLIMIT
                //////////////////////////////
                //
                // fixup negative densities
                //
                //////////////////////////////
                fixup_negdensities(EOMSETRAD,&fixedrad, startpl, endpl, i, j, k, radlpflag, pv,ptoavg, ptrgeom, pr0, ucons, finalstep);

                if(fixedrad==1 && (startpl==URAD0 && endpl==URAD0)){
                  // assume success for densities, so no longer nogood.
                  nogood=0;
                }
                if(fixedrad==1 && (startpl<=URAD0 && endpl>=URAD1)){
                  // then fixup but only changed densities, so still need to process non-densities
                  startpl=URAD1; // start at U1 (first velocity) and finish at same ending if was ending on some velocity
                  fixedrad=0; // then reset fixedmhd->0 so can still modify these remaining quantities -- otherwise v^i would be unchanged even if wanted to average that out for (e.g.) negative density results for inversions.
                }
              }

              //              dualfprintf(fail_file,"Trying no good? nogood=%d\n",nogood);

              /////////////////////
              //
              // If no good surrounding value found, then average bad values in some way
              //
              //////////////////////
              if(1){ // should stay 1, see below.
                if(nogood){
                  //                  dualfprintf(fail_file,"fixupnogood: rad %d %d %d : %d %d\n",i,j,k,startpl,endpl);
                  fixup_nogood(startpl, endpl, i, j, k, doingmhd, mhdlpflag, radlpflag, pv,ptoavg, pbackup, ptrgeom);
                }
              }
              else{
                // normally assume u2p_rad() did ok local fixup if non-local fixup doesn't work out.
                // NO, reset failure to non-failure if ok local fixup.  If here, want fixup even of radiation quantity, such as when implicit solver fails.
              }

              //              if(MACP0A1(pv,i,j,k,URAD0)<1.5*ERADLIMIT){
              //                dualfprintf(fail_file,"Why not fixed? %d %d %d: %21.15g err=%d\n",i,j,k,MACP0A1(pv,i,j,k,URAD0),radlpflag);
              //                myexit(0);
              //              }


              /////////////////////
              //
              // Things to do only if modifying radiation velocity
              //
              //////////////////////
              if(startpl<=PRAD1 && endpl>=PRAD3){

#if(WHICHVEL==VELREL4)
                //////////////
                //
                // check gamma to so calibrate new gamma to no larger than previous gamma
                /////////////
                MYFUN(gamma_calc(&MACP0A1(ptoavg,i,j,k,URAD1-U1),ptrgeom,&gamma,&qsq),"fixup_utoprim: gamma calc failed\n","fixup.c",1);
                if (ucon_calc(&MACP0A1(ptoavg,i,j,k,URAD1-U1), ptrgeom, ucon, others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);
                //   alpha = 1. / sqrt(-ptrgeom->gcon[GIND(0,0)]);
                alpha = ptrgeom->alphalapse;
                vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);

                limitedgamma=0;
                if(gamma>GAMMAMAXRAD){
                  limitedgamma=1;
                  if(debugfail>=2) dualfprintf(fail_file,"RAD: initial gamma: %21.15g,  max: %21.15g, initial vsq: %21.15g\n",gamma,GAMMAMAXRAD,vsq);
                  // limit:
                  gamma=GAMMAMAXRAD;
                }
#else
                limitedgamma=0;
                if (ucon_calc(&MACP0A1(ptoavg,i,j,k,URAD1-U1), ptrgeom, ucon,others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);
#endif
   

                /////////////
                //
                // check new gamma to make sure smaller than original (i.e. for pv, not original ptoavg)
                //
                /////////////
     
#if(WHICHVEL==VELREL4)
                //                if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
                //                  dualfprintf(fail_file,"BEFORE IN FIXUPUTOPRIM LIMITGAMMA2: finalstep=%d radlpflag=%d\n",finalstep,radlpflag);
                //                }
                if(limitgammarad=limit_gamma(0,GAMMAMAX,gamma,MAC(pv,i,j,k),MAC(ucons,i,j,k),ptrgeom,-1)>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 2);
                //                if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
                //                  dualfprintf(fail_file,"AFTER IN FIXUPUTOPRIM LIMITGAMMA2: finalstep=%d\n",finalstep);
                //                }

                if(debugfail>=3){
                  if(limitedgamma){
                    // check gamma
                    MYFUN(gamma_calc(&MACP0A1(pv,i,j,k,URAD1-U1),ptrgeom,&gamma,&qsq),"fixup_utoprim: gamma calc failed\n","fixup.c",2);
                    dualfprintf(fail_file,"RAD: final gamma: %21.15g\n",gamma);
                  }
                }
#endif

              }// end if dealing with velocity



              /////////////////////
              //
              // Things to do only if modifying density (must come after modifying velocity that enters ucon[TT])
              //
              //////////////////////
              PLOOP(pliter,pl){
                if(pl==YFL4 || pl==YFL5){
                  if(startpl<=pl && endpl>=pl){


                    if(DO_CONSERVE_D_INFAILFIXUPS){
                      
                      if(mhdlpflag==UTOPRIMFAILGAMMAPERC || 1){ // GODMARK: always doing it

                        if (ucon_calc(&MACP0A1(pv,i,j,k,URAD1-U1), ptrgeom, ucon,others) >= 1) FAILSTATEMENT("fixup.c:utoprimfail_fixup()", "ucon_calc()", 1);

                        if(finalstep==1){
                          // this is not reasonable, because ignores any mass injection required during failure to make sense of solution
                          FTYPE Unothing[NPR];
                          UtoU(UEVOLVE,UNOTHING,ptrgeom,MAC(ucons,i,j,k),Unothing);
                          D0 = Unothing[pl] ;
                        }
                        else{
                          // This is probably not necessary or useful
                          D0 = MACP0A1(ptoavg,i,j,k,pl)*ucon[TT];
                        }
                        
                        ///////////////////////////////////////////
                        //
                        // constrain change in density so conserve particle number
                        //
                        //////////////////////////////////////////
                        if(pl==YFL4 || pl==YFL5){
                          MACP0A1(pv,i,j,k,pl) = D0/ucon[TT];
                        }
                      }
                    }
                  }
                }
              }// end over density or scalars





#if(0)
              // DEBUG problem of launch with pressureless stellar model collapse
              PALLLOOP(pl) dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d j=%d k=%d pl=%d pv=%21.15g ptoavg=%21.15g\n",nstep,steppart,i,j,k,pl,MACP0A1(pv,i,j,k,pl),MACP0A1(ptoavg,i,j,k,pl));
#endif


              //              int finaluu=0;
              //              consfixup_1zone(finaluu,i,j,k,ptrgeom, &MACP0A1(pv,i,j,k,0),NULL);


            } // end if fixedrad==0
          }// end if not keeping static
        }// end if radiation failed 
 

        /////////////////////////////////
        //
        // ACCOUNTING (static or average)
        //
        /////////////////////////////////
        fixuputoprim_accounting(i, j, k, mhdlpflag, radlpflag, limitgammamhd, limitgammarad, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);


      }// if fixing mhd or fixing radiation  
   
    }// end over COMPZLOOP loop
  }// end over parallel region

  return(0);
}








/// report problems, but doesn't fix them
int fixup_utoprim_nofixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR];

  // this average only works if using 4 velocity since only then guaranteed solution is good after interpolation
  if(WHICHVEL==VEL3) return(0); // just stick with static, best can do
  if(EOMTYPE==EOMFFDE || EOMTYPE==EOMFFDE2) return(0); // nothing to do


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
    PFTYPE mhdlpflag,radlpflag;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;


    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPZLOOP;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      // get failure flag
      mhdlpflag=GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL);
      radlpflag=GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMRADFAIL);


      if(IFUTOPRIMFAIL(mhdlpflag)||IFUTOPRIMFAIL(radlpflag)){
 
        PALLLOOP(pl)    pr0[pl]=MACP0A1(ptoavg,i,j,k,pl);
        get_geometry(i,j,k,CENT,ptrgeom);
 
        /////////////////////////////////
        //
        // ACCOUNTING (static or average)
        //
        /////////////////////////////////
        fixuputoprim_accounting(i, j, k, mhdlpflag, radlpflag, 0, 0, GLOBALPOINT(pflag),pv,ptoavg, ptrgeom, pr0, ucons, finalstep);
   
      }// end if failure
    }// end over COMPZLOOP loop
  }// end over parallel region

  return(0);
}









/// fixup negative densities
static int fixup_negdensities(int whicheomset, int *fixed, int startpl, int endpl, int i, int j, int k, PFTYPE lpflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom, FTYPE *pr0, FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE prguess[NPR],prceiling[NPR];


  if(whicheomset==EOMSETMHD){
    ////////////////////////////
    //
    // deal with negative internal energy as special case of failure
    //

    if(*fixed!=0){
      if(lpflag==UTOPRIMFAILUNEG){
   
        if(STEPOVERNEGU==NEGDENSITY_NEVERFIXUP){ if(HANDLEUNEG==1) *fixed=1; }
        else if((STEPOVERNEGU==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP && finalstep)){

          if(HANDLEUNEG==1){
            // set back to floor level
            set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess,prceiling);
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
      if(lpflag==UTOPRIMFAILRHONEG){
   
        if(STEPOVERNEGRHO==NEGDENSITY_NEVERFIXUP){ if(HANDLERHONEG) *fixed=1; }
        else if((STEPOVERNEGRHO==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGRHO==NEGDENSITY_FIXONFULLSTEP && finalstep)){

          if(HANDLERHONEG==1){
            // set back to floor level
            set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess,prceiling);
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
      if(lpflag==UTOPRIMFAILRHOUNEG){

        if(STEPOVERNEGRHOU==NEGDENSITY_NEVERFIXUP){ if(HANDLERHOUNEG) *fixed=1; }
        // GODMARK: Why use STEPOVERNEGU and STEPOVERNEGRH instead of STEPOVERNEGRHOU below?
        else if( (STEPOVERNEGRHOU==NEGDENSITY_ALWAYSFIXUP)  ||(STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP && STEPOVERNEGRHO==NEGDENSITY_FIXONFULLSTEP && finalstep)){

          if(HANDLERHOUNEG==1){
            // set back to floor level
            set_density_floors(ptrgeom,MAC(pv,i,j,k),prguess,prceiling);
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
  }
  else if(whicheomset==EOMSETRAD){

    ////////////////////////////
    //
    // deal with negative internal energy as special case of failure
    //

    if(*fixed!=0){
      if(lpflag==UTOPRIMRADFAILERFNEG){
   
        if(STEPOVERNEGU==NEGDENSITY_NEVERFIXUP){ if(HANDLEUNEG==1) *fixed=1; }
        else if((STEPOVERNEGU==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP && finalstep)){

          if(HANDLEUNEG==1){
            // set back to floor level
            prguess[URAD0]=ERADLIMIT;
  
            if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED){
              // then pv is previous timestep value and can use to make fix
              if(-MACP0A1(pv,i,j,k,URAD0)<prguess[URAD0]){ *fixed=1; MACP0A1(pv,i,j,k,URAD0)=prguess[URAD0];} // otherwise assume really so bad that failure
            }
            else{
              // just treat as floor for all failures since do not know what updated quantity is
              MACP0A1(pv,i,j,k,UU)=prguess[URAD0];
              *fixed=1;
            }
          }// end if handling u<zerouuperbaryon*prim[RHO] in special way
        }// end if not allowing negative u or if allowing but not yet final step
        else if((STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP)&&(!finalstep)){
          if(HANDLEUNEG==1) *fixed=1; // tells rest of routine to leave alone and say ok solution, but don't use it to fix convergence failures for other zones
        }
      }// end if u<zerouuperbaryon*prim[RHO]
    }// end if not fixed


  }
  

 
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

/// ACCOUNTING (under any circumstance, static or average) when modified quantities due to inversion issue
static int fixuputoprim_accounting(int i, int j, int k, PFTYPE mhdlpflag, PFTYPE radlpflag, int limitgammamhd, int limitgammarad, PFTYPE (*lpflag)[NSTORE2][NSTORE3][NUMPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom, FTYPE *pr0, FTYPE (*uconsinput)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  PFTYPE mhdutoprimfailtype,radutoprimfailtype;
  int doadjustcons;
  struct of_state q;
  FTYPE (*utoinvert)[NSTORE2][NSTORE3][NPR];
  int docorrectucons;
  int pliter,pl;



  // store ucons and only change if needed (and handle avoidance of B1,B2,B3 in case ucons points to unewglobal staggered field, while here we operate on centers)
  FTYPE ucons[NPR];
  if(uconsinput!=NULL){
    PLOOP(pliter,pl) ucons[pl]=MACP0A1(uconsinput,i,j,k,pl);
  }

  //////////////////////
  //
  // MHD
  //
  //////////////////////

  // account for changes by tracking conserved quantities
  // note that if new utoprim solution was impossible, we are using here prior solution as new state, which means the diagnostics won't be acccurate.  There's no way to fix this unless the fluxes always give a well-defined new primitive variable.
  // specific value of the mhdlpflag>0 communicates the type of failure
  if(IFUTOPRIMNOFAIL(mhdlpflag)){
    mhdutoprimfailtype=-1;
    docorrectucons=0;
  }
  else if(
          (mhdlpflag==UTOPRIMFAILCONV)|| // only used by 5D method currently
          (mhdlpflag==UTOPRIMFAILCONVGUESSUTSQ)|| // rest are only used by 1D/2D method currently
          (mhdlpflag>=UTOPRIMFAILCONVRET)||
          (mhdlpflag==UTOPRIMFAILCONVW)||
          (mhdlpflag==UTOPRIMFAILCONVUTSQVERYBAD)||
          (mhdlpflag==UTOPRIMFAILNANGUESS)||
          (mhdlpflag==UTOPRIMFAILNANRESULT)||
          (mhdlpflag==UTOPRIMFAILCONVBADINVERTCOMPARE)||
          (mhdlpflag==UTOPRIMFAILCONVUTSQ)||
          (mhdlpflag==UTOPRIMFAILFAKEVALUE) // fake value for avoiding MPI boundary call and so MPI boundary values for fixup
          ){
    mhdutoprimfailtype=COUNTUTOPRIMFAILCONV;
    docorrectucons=1;
  }
  else if(mhdlpflag==UTOPRIMFAILRHONEG){
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.
    if(STEPOVERNEGRHO==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGRHO==1){
        if(finalstep){
          mhdutoprimfailtype=COUNTUTOPRIMFAILRHONEG;
          docorrectucons=1;
        }
        else{
          mhdutoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGRHO==2){
        mhdutoprimfailtype=COUNTUTOPRIMFAILRHONEG;
        docorrectucons=1;
      }
      else{
        mhdutoprimfailtype=-1;
        docorrectucons=0;
      }
    } 
    else if((STEPOVERNEGRHO)&&(!finalstep)){
      mhdutoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      mhdutoprimfailtype=COUNTUTOPRIMFAILRHONEG;
      docorrectucons=1;
    }
  }
  else if(mhdlpflag==UTOPRIMFAILUNEG || mhdlpflag==UTOPRIMFAILU2AVG1|| mhdlpflag==UTOPRIMFAILU2AVG2 || mhdlpflag==UTOPRIMFAILU2AVG1FROMCOLD|| mhdlpflag==UTOPRIMFAILU2AVG2FROMCOLD || mhdlpflag==UTOPRIMFAILURHO2AVG1FROMFFDE){ // GODMARK: maybe want separate accounting
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.

    // default counting:
    if(mhdlpflag==UTOPRIMFAILU2AVG1FROMCOLD|| mhdlpflag==UTOPRIMFAILU2AVG2FROMCOLD) mhdutoprimfailtype=COUNTCOLD;
    else mhdutoprimfailtype=COUNTUTOPRIMFAILUNEG;

    // default counting:
    if(mhdlpflag==UTOPRIMFAILURHO2AVG1FROMFFDE) mhdutoprimfailtype=COUNTFFDE;
    else mhdutoprimfailtype=COUNTUTOPRIMFAILUNEG;

    // now set whether to ucons correction or override counting
    if(STEPOVERNEGU==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGU==1){
        if(finalstep){
          docorrectucons=1;
        }
        else{
          mhdutoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGU==2){
        docorrectucons=1;
      }
      else{
        mhdutoprimfailtype=-1;
        docorrectucons=0;
      }
    } 
    else if((STEPOVERNEGU)&&(!finalstep)){
      mhdutoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      docorrectucons=1;
    }
  }
  else if(mhdlpflag==UTOPRIMFAILRHOUNEG){
    // whether to count uneg as failure in diagnostic reporting or not
    // should really have a new diagnostic for substep u<zerouuperbaryon*prim[RHO] 's.
    if(STEPOVERNEGRHOU==NEGDENSITY_NEVERFIXUP){
      if(DOCOUNTNEGRHOU==1){
        if(finalstep){
          mhdutoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
          docorrectucons=1;
        }
        else{
          mhdutoprimfailtype=-1;
          docorrectucons=0;
        }
      }
      else if(DOCOUNTNEGRHOU==2){
        mhdutoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
        docorrectucons=1;
      }
      else{
        mhdutoprimfailtype=-1;
        docorrectucons=0;
      }
    } 
    else if((STEPOVERNEGRHOU)&&(!finalstep)){
      mhdutoprimfailtype=-1;
      docorrectucons=0;
    }
    else{
      mhdutoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
      docorrectucons=1;
    }
  }
  else if(mhdlpflag==UTOPRIMFAILGAMMAPERC){
    mhdutoprimfailtype=COUNTGAMMAPERC;
    docorrectucons=1;
  }
  else if(mhdlpflag==UTOPRIMFAILUPERC){
    mhdutoprimfailtype=COUNTUPERC;
    docorrectucons=1;
  }
  else if(mhdlpflag==UTOPRIMFAILFIXEDENTROPY){
    mhdutoprimfailtype=COUNTENTROPY;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mhdlpflag==UTOPRIMFAILFIXEDCOLD){
    mhdutoprimfailtype=COUNTCOLD;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mhdlpflag==UTOPRIMFAILFIXEDBOUND1){
    mhdutoprimfailtype=COUNTBOUND1;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mhdlpflag==UTOPRIMFAILFIXEDBOUND2){
    mhdutoprimfailtype=COUNTBOUND2;
    docorrectucons=0; // account, but don't change conserved quantities
  }
  else if(mhdlpflag==UTOPRIMFAILFIXEDUTOPRIM){
    dualfprintf(fail_file,"prior pflag not cleared: nstep=%ld steppart=%d t=%21.15g i=%d j=%d k=%d \n",nstep,steppart,t,i,j,k);
    mhdutoprimfailtype=-1;
    docorrectucons=0;
  }
  else{
    dualfprintf(fail_file,"No such pflag failure type: %d for t=%21.15g nstep=%ld steppart=%d i=%d j=%d k=%d\n",mhdlpflag,t,nstep,steppart,i,j,k);
    myexit(1);
  }


  //////////////////////
  //
  // RADIATION
  //
  //////////////////////

  if(IFUTOPRIMNOFAIL(radlpflag)){
    radutoprimfailtype=-1;
    docorrectucons=0;
  }
  else if(radlpflag>=UTOPRIMRADFAILCASE1A || radlpflag<=UTOPRIMRADFAILCASE3B || radlpflag==UTOPRIMRADFAILERFNEG){
    radutoprimfailtype=COUNTRADLOCAL;
    docorrectucons=1;
  }
  else if(radlpflag>=UTOPRIMRADFAILBAD1 || radlpflag==UTOPRIMRADFAILGAMMAHIGH || radlpflag==UTOPRIMRADFAILERFNEG){
    radutoprimfailtype=COUNTRADNONLOCAL;
    docorrectucons=1;
  }
  else if(radlpflag==UTOPRIMRADFAILFIXEDUTOPRIMRAD){
    dualfprintf(fail_file,"prior rad pflag not cleared: nstep=%ld steppart=%d t=%21.15g i=%d j=%d k=%d \n",nstep,steppart,t,i,j,k);
    radutoprimfailtype=-1;
    docorrectucons=0;
  }



  //////////////////////
  //
  // Check if either MHD or RADIATION was corrected
  //
  //////////////////////

  if(mhdutoprimfailtype!=-1 || radutoprimfailtype!=-1 || limitgammamhd==-1 || limitgammarad==-1){ 
    if(docorrectucons==1){
      docorrectucons=(DOENOFLUX != NOENOFLUX);
    }

    ////////////////
    //
    // reset true pflag counter to "no" (fixed) failure so diagnostics account for change
    //
    ////////////////
    MACP0A1(lpflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMFAILFIXEDUTOPRIM;
    MACP0A1(lpflag,i,j,k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILFIXEDUTOPRIMRAD;



    // diagnostics
    //    FTYPE prdiag[NPR],pr[NPR];
    //    PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
    //    int doingmhdfixup=(mhdutoprimfailtype!=-1); // diag_fixup() accounting diagnostics part doesn't have radiation yet

    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"BEFORE IN FIXUPUTOPRIM: finalstep=%d\n",finalstep);
    //    }


    //    diag_fixup(docorrectucons,prdiag, MAC(pv,i,j,k), ucons, ptrgeom, finalstep, doingmhdfixup, (int)mhdutoprimfailtype); // do corrections in general, but only include accounting for MHD case so far (KORALTODO, not too crucial, just diags).
    //    PLOOP(pliter,pl) prdiag[pl]=MACP0A1(pv,i,j,k,pl);

    // pr0=prdiag only correct as reference primitive if was no failure of any kind, but some kind of failure if here, so need to use ucons as reference to what should be.
    FTYPE Uievolve[NPR];
    PLOOP(pliter,pl) Uievolve[pl] = MACP0A1(uconsinput,i,j,k,pl);
    diag_fixup_Ui_pf(docorrectucons, Uievolve, MAC(pv,i,j,k), ptrgeom, finalstep,(int)mhdutoprimfailtype, NULL);

    



    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"AFTER IN FIXUPUTOPRIM: UievolveRHO=%21.15g\n",Uievolve[RHO]*dx[1]*dx[2]*dx[3]);
    //    }


    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      struct of_state qb;
    //      get_state(MAC(pv,i,j,k),ptrgeom,&qb);
    //      FTYPE ub[NPR],ubabs[NPR];
    //      int uutype=UDIAG;
    //      primtoU(uutype,MAC(pv,i,j,k),&qb,ptrgeom,ub,ubabs);
    //      dualfprintf(fail_file,"URHOINFIXUPUTOPRIMACCTT=%21.15g\n",ub[RHO]*dx[1]*dx[2]*dx[3]);
    //    }



    // now adjust uf or ucum so agrees [already done in diag_fixup() but only if final step.  This allows one to do it on any step in case want to control behavior of conserved quantities or primitive quantities used to create fluxes.
    if(ADJUSTCONSERVEDQUANTITY&&docorrectucons){

      if((ADJUSTCONSERVEDQUANTITY==1)&&(finalstep)) doadjustcons=1;
      else if(ADJUSTCONSERVEDQUANTITY==2) doadjustcons=1;
      else doadjustcons=0;

      if(doadjustcons){
        //utoinvert=ucons;
        //     if(finalstep){ // last call, so ucum is cooked and ready to eat!
        //     }
        //     else{ // otherwise still iterating on primitives
        //       utoinvert=ulast;
        //     }
     
        MYFUN(get_state(MAC(pv,i,j,k), ptrgeom, &q),"fixup.c:fixup_utoprim()", "get_state()", 1);
        MYFUN(primtoU(UEVOLVE,MAC(pv,i,j,k), &q, ptrgeom, ucons, NULL),"fixup.c:fixup_utoprim()", "primtoU()", 1);
      }
    }

    // copy over ucons result in case changed.
    if(uconsinput!=NULL && docorrectucons){
      PLOOP(pliter,pl){
        // if staggered field, avoid modifying field since at FACEs, not CENT where pr lives.
        if(BPL(pl)==0 && FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){
          MACP0A1(uconsinput,i,j,k,pl) = ucons[pl];
        }
      }
    }

  }// end if mhdutoprimfailtype!=-1 || radutoprimfailtype!=-1



#if(DOSUPERDEBUG)
  // only used to keep track of normal non-failure points
  // some non-utoprim-failure points will be picked up by this.  I.E. those with limitgammaact/flooract/inflowact
  if(mhdutoprimfailtype==-1) superdebug_utoprim(pr0,MAC(pv,i,j,k),ptrgeom,mhdutoprimfailtype);
  // collect values for non-failed and failed zones
#endif


  return(0);
}


// what quantities to average or treat in fixup_utoprim()
// only touch mhd quantities for mhd case, rad quantities for rad case, and never touch field
#define PLOOPSTARTEND(pl) for(pl=startpl;pl<=endpl;pl++) if((NONRADFULLPL(pl) && doingmhd || RADFULLPL(pl) && doingmhd==0) && BPL(pl)==0)




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

/// more general averaging procedure that tries to pick best from surrounding good values
static int general_average(int useonlynonfailed, int numbndtotry, int maxnumbndtotry, int startpl, int endpl, int i, int j, int k, int doingmhd, PFTYPE mhdlpflag, PFTYPE radlpflag, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3][NUMFAILPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom)
{
  int doavgcausal,failavglooptype;
  int pl,pliter;
  int ii,jj,kk;
  int jjj;
  struct of_state state_pv;
  FTYPE cminmax[NDIM][NUMCS];
  FTYPE superfast[NDIM];
  int numavg[NPR];
  int numavg0,numavg1[NPR];
  FTYPE mysum[2][NPR];
  int numupairs,qq,thisnotfail,thatnotfail,rnx,rny,rnz;
  FTYPE ref[NPR];
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE lastmin[NPR],avganswer0[NPR],avganswer1[NPR];
  int ignorecourant=0;
  int factor;
  PFTYPE lpflag;
  int whichpflag;

  
  prod0dualfprintf(debugfail>=2,fail_file,"uc2generalA: doingmhd=%d mhdlpflag=%d radlpflag=%d startpl=%d endpl=%d :: i=%d j=%d k=%d\n",doingmhd,mhdlpflag,radlpflag,startpl,endpl,i,j,k); // not too much
  //  if(debugfail>=2) for(pl=startpl; pl<=endpl; pl++) dualfprintf(fail_file,"uc2generalSTART: pl=%d pv=%g\n",pl,MACP0A1(pv,i,j,k,pl));


  if(doingmhd==1){ // if MHD
    whichpflag=FLAGUTOPRIMFAIL;
    lpflag=mhdlpflag;

    if(mhdlpflag==UTOPRIMFAILU2AVG1 || mhdlpflag==UTOPRIMFAILU2AVG2 || mhdlpflag==UTOPRIMFAILU2AVG1FROMCOLD || mhdlpflag==UTOPRIMFAILU2AVG2FROMCOLD){
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
      failavglooptype=0; // must stay 0 since can't assume a density of some kind with min value
    }  
  }
  else{ // if RAD
    whichpflag=FLAGUTOPRIMRADFAIL;
    lpflag=radlpflag;
    doavgcausal=0;
    failavglooptype=0; // must stay 0
  }




 
  if(doavgcausal){
    //////////
    // first get wave speed to check if superfast so need causal averaging
    MYFUN(get_state(MAC(ptoavg,i,j,k), ptrgeom, &state_pv),"fixup.c:fixup_utporim()", "get_state()", 1);
    if(N1NOT1){
      MYFUN(vchar_all(MAC(ptoavg,i,j,k), &state_pv, 1, ptrgeom, &cminmax[1][CMAX], &cminmax[1][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar_all() dir=1or2", 1);
    }
    if(N2NOT1){
      MYFUN(vchar_all(MAC(ptoavg,i,j,k), &state_pv, 2, ptrgeom, &cminmax[2][CMAX], &cminmax[2][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar_all() dir=1or2", 2);
    }
    if(N3NOT1){
      MYFUN(vchar_all(MAC(ptoavg,i,j,k), &state_pv, 3, ptrgeom, &cminmax[3][CMAX], &cminmax[3][CMIN],&ignorecourant),"fixup.c:fixup_utoprim()", "vchar_all() dir=1or2", 3);
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
  numavg0=0;
  PLOOPSTARTEND(pl){
    numavg1[pl]=0;
    numavg[pl]=0;
    mysum[0][pl]=0.0;
    mysum[1][pl]=0.0;
    lastmin[pl]=VERYBIG;
  }
  // size of little averaging box, from -NUMPFLAGBNDx to +NUMPFLAGBNDx
  int rbndx=MIN(numbndtotry,NUMPFLAGBND1);
  int rbndy=MIN(numbndtotry,NUMPFLAGBND2);
  int rbndz=MIN(numbndtotry,NUMPFLAGBND3);
  rnx=(SHIFT1==1) ? (rbndx*2+1) : 1;
  rny=(SHIFT2==1) ? (rbndy*2+1) : 1;
  rnz=(SHIFT3==1) ? (rbndz*2+1) : 1;
     
  // number of unique pairs ( NUMPFLAGBND1*NUMPFLAGBND2*NUMPFLAGBND3 + NUMPFLAGBND1
  numupairs=(int)(rnx*rny*rnz-1)/2;
     
  for(qq=0;qq<numupairs;qq++){
       
    // 1-d to 3D index
    ii=(int)(qq%rnx)-rbndx;
    jj=(int)((qq%(rnx*rny))/rnx)-rbndy;
    kk=(int)(qq/(rnx*rny))-rbndz;

    if(useonlynonfailed==1){
      if(doingmhd){
        thisnotfail=IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i+ii,j+jj,k+kk,whichpflag));
        thatnotfail=IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i-ii,j-jj,k-kk,whichpflag));
      }
      else{
        thisnotfail=IFUTOPRIMRADNOFAILORFIXED(MACP0A1(lpflagfailorig,i+ii,j+jj,k+kk,whichpflag));
        thatnotfail=IFUTOPRIMRADNOFAILORFIXED(MACP0A1(lpflagfailorig,i-ii,j-jj,k-kk,whichpflag));
      }
    }
    else{
      // here lpflagfailorig can be just NULL since not used
      thisnotfail=thatnotfail=1; // just ignores whether failed
    }




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
    // bigger than myself and absolute limit
    PLOOPSTARTEND(pl){
      ref[pl]=0.0; // default
      if(pl==RHO) ref[pl]=MAX(1.001*RHOMINLIMIT,MACP0A1(ptoavg,i,j,k,pl)); // MHD
      else if(pl==UU) ref[pl]=MAX(1.001*UUMINLIMIT,MACP0A1(ptoavg,i,j,k,pl)); // MHD
      else if(pl==URAD0) ref[pl]=MAX(1.001*ERADLIMIT,MACP0A1(ptoavg,i,j,k,URAD0)); // RAD
    }
#elif(0)
    // bigger than 0
    PLOOPSTARTEND(pl) ref[pl]=0.0;
#endif

    //     if(failavglooptype==1  || failavglooptype==2){  // only use if positive

    // number of quantities one summed
    numavg0+=thisnotfail+thatnotfail;
     
    // do sum or min of these 2 pairs
    PLOOPSTARTEND(pl){
      ftemp1=MACP0A1(ptoavg,i+ii,j+jj,k+kk,pl);
      ftemp2=MACP0A1(ptoavg,i-ii,j-jj,k-kk,pl);

      /////////
      // AVG type
      /////////
      mysum[qq%2][pl]+=ftemp1*thisnotfail + ftemp2*thatnotfail;

      /////////
      // MIN type
      /////////
#if(0)
      if(ftemp1>ref[pl] && thisnotfail){
        lastmin[pl]=MIN(lastmin[pl],ftemp1); // smallest positive number
        numavg1[pl]++;
      }
      
      if(ftemp2>ref[pl] && thatnotfail){
        lastmin[pl]=MIN(lastmin[pl],ftemp2);
        numavg1[pl]++;
      }
#else
      //      if(pl==8) dualfprintf(fail_file,"ii=%d jj=%d kk=%d i+ii=%d j+jj=%d k+kk=%d pl=%d testing: %g>%g && %g>%g : %d %d\n",ii,jj,kk,i+ii,j+jj,k+kk,pl,ftemp1,ref[pl],ftemp2,ref[pl],thisnotfail,thatnotfail);
      if(ftemp1>ref[pl] && thisnotfail && ftemp2>ref[pl] && thatnotfail){
        lastmin[pl]=MIN(MIN(lastmin[pl],ftemp1),ftemp2); // smallest positive number if both of pair are larger than my value (else keep my small value)
        numavg1[pl]+=2;
      }
      else if(ftemp1>ref[pl] && thisnotfail){
        lastmin[pl]=MIN(lastmin[pl],ftemp1); // smallest positive number if both of pair are larger than my value (else keep my small value)
        numavg1[pl]++;
      }
      else if(ftemp2>ref[pl] && thatnotfail){
        lastmin[pl]=MIN(lastmin[pl],ftemp2); // smallest positive number if both of pair are larger than my value (else keep my small value)
        numavg1[pl]++;
      }
#endif
    }// pl loop
  } //end loop over pairs

 



  ///////////////////////
  //
  // all loops over surrounding points is done, now get average answer
  //
  ////////////////////////
  //////////
  // NORMAL:
  //////////
  if(numavg0!=0){
    PLOOPSTARTEND(pl){
      avganswer0[pl]=(mysum[0][pl]+mysum[1][pl])/((FTYPE)(numavg0));
    }
  }
  PLOOPSTARTEND(pl){
    avganswer1[pl]=lastmin[pl];
  }
  
   

  ///////////////
  //
  // choose final answer (avg or min)
  //
  ///////////////

  int doingavgtype[NPR];
  PLOOPSTARTEND(pl){
    doingavgtype[pl]=-1; // default avoidance of pl

    // default
    if(failavglooptype==0 || ((failavglooptype==2)&&(numavg1[pl]==0)) ){
      doingavgtype[pl]=1;
    }
    if(failavglooptype==1 || ((failavglooptype==2)&&(numavg0==0)) ){  
      doingavgtype[pl]=0;
    }


    // override
    if(doingmhd==0 && radlpflag==UTOPRIMRADFAILERFNEG && pl==URAD0 || doingmhd==1 && lpflag==UTOPRIMFAILRHONEG && pl==RHO || doingmhd==1 && lpflag==UTOPRIMFAILUNEG && pl==UU || doingmhd==1 && lpflag==UTOPRIMFAILRHOUNEG && (pl==RHO || pl==UU)){
      doingavgtype[pl]=0; // min type
    }
    else doingavgtype[pl]=1; // avg type
  }



  int numavgfinal=256; // just large number (at least (2*NxBND)**3
  PLOOPSTARTEND(pl){// should only be over specific density(s)
    // choose min
    if(numavg1[pl]!=0 && doingavgtype[pl]==0){
      MACP0A1(pv,i,j,k,pl)=avganswer1[pl]; // else keep same as original answer
      //        dualfprintf(fail_file,"HERE: pl=%d %21.15g %d\n",pl,MACP0A1(pv,i,j,k,pl),numavg1[pl]);
      numavg[pl]=numavg1[pl];
    }
    // choose avg
    if(numavg0!=0 && doingavgtype[pl]==1){
      MACP0A1(pv,i,j,k,pl)=avganswer0[pl];
      numavg[pl]=numavg0;
    }
    numavgfinal=MIN(numavg[pl],numavgfinal); // get minimum number of good points over those wanted "average".  If one was bad, have to do nogood version
  }// end pl loop

  

  if( lpflag==UTOPRIMFAILU2AVG2){ // lpflag here is either for MHD or RAD depending upon doingmhd=1 or 0
    if(numbndtotry==maxnumbndtotry) numavgfinal++; // assume at least always one good one so don't treat as real failure if no good values surrounding.  But only do so if last chance for no good average.
  }
  // else use real numavgfinal
   
   
  if(numavgfinal==0){
    return(1);
  }
   
  // good value exist
  return(0);
}




/// simple average of good surrounding zones
static int simple_average(int startpl, int endpl, int i, int j, int k,int doingmhd, PFTYPE (*lpflagfailorig)[NSTORE2][NSTORE3][NUMFAILPFLAGS],FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  int whichpflag;

  prod0dualfprintf(debugfail>=2,fail_file,"uc2simple: i=%d j=%d k=%d\n",i,j,k); // not too much



  if(doingmhd==1){ // if MHD
    whichpflag=FLAGUTOPRIMFAIL;
  }
  else{ // if RAD
    whichpflag=FLAGUTOPRIMRADFAIL;
  }


  /////////////
  //
  // 4 VALUES
  //
  /////////////
  if( // but if surrounded by good values
     (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jp1mac(j),k,whichpflag)))&&
     (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jm1mac(j),k,whichpflag)))&&
     (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),j,k,whichpflag)))&&  
     (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),j,k,whichpflag)))    
      ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected1\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    // don't mess with B field, it's correct already
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=AVG4_1(ptoavg,i,j,k,pl);  
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if surrounded by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jp1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jm1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jp1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected2\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    PLOOPSTARTEND(pl){
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
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jp1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected3\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=AVG2_1(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),j,k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),j,k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected4\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=AVG2_2(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jp1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected5\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=AVG2_3(ptoavg,i,j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good values
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jm1mac(j),k,whichpflag)))&&
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jp1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected6\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then average
    PLOOPSTARTEND(pl){
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
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jp1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),jp1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),j,k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,ip1mac(i),jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,ip1mac(i),jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,i,jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jm1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),jm1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),j,k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),j,k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,im1mac(i),jp1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
      MACP0A1(pv,i,j,k,pl)=MACP0A1(ptoavg,im1mac(i),jp1mac(j),k,pl);
      if(debugfail>=2) dualfprintf(fail_file,"uc2: pl=%d pv=%21.15g\n",pl,MACP0A1(pv,i,j,k,pl));
    }
  }
  else if( // but if "surrounded" by good value
          (IFUTOPRIMNOFAILORFIXED(MACP0A1(lpflagfailorig,i,jp1mac(j),k,whichpflag)))
           ){
    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d k=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j,startpos[3]+k);
    // then ASSIGN
    PLOOPSTARTEND(pl){
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


/// how to fixup when no good surrounding values to use
//static int fixup_nogood(int startpl, int endpl, int i, int j, int k, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom)
static int fixup_nogood(int startpl, int endpl, int i, int j, int k, int doingmhd, PFTYPE mhdlpflag, PFTYPE radlpflag, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*ptoavg)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], struct of_geom *ptrgeom)
{
  int pl,pliter;
  FTYPE prguess[NPR];
  FTYPE (*ptoavgwhennogood)[NSTORE2][NSTORE3][NPR];
  int ti,tj,tk,resetregion;


  
  prod0dualfprintf(debugfail>=2,fail_file,"uc2nogood: i=%d j=%d k=%d\n",i,j,k); // not too much


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
    // if(limit_gamma(0,1.0,MAC(pv,i,j,k),ptrgeom,-1)>=1)
    //   FAILSTATEMENT("fixup.c:fixup_utoprim()", "limit_gamma()", 1);

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
#if(GENERALAVERAGE==1)
    // use general average without regard to failure
    int useonlynonfailed=0; // any type of cell, failed or not, will be used in the average
    int numbndtotry=1;
    int maxnumbndtotry=1; //MAXNUMPFLAGBND; // choice smaller or equal to MAXNUMPFLAGBND
    int nogood=general_average(useonlynonfailed,numbndtotry,maxnumbndtotry,startpl, endpl, i, j, k, doingmhd, mhdlpflag, radlpflag, NULL ,pv,ptoavg,ptrgeom);
    // nogood will be 1
#else
    PLOOPSTARTEND(pl) MACP0A1(pv,i,j,k,pl)=0.5*(AVG4_1(ptoavgwhennogood,i,j,k,pl)+AVG4_2(ptoavgwhennogood,i,j,k,pl)); // like simple average
#endif
    
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




/// OPENMPNOTE: Assume superdebug_utoprim() not called by multiple threads (i.e. when debugging, not using multiple threads), so static's are ok (including firsttime)
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
    failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Ui, NULL);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
    
    // after any changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Uf, NULL);
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

/// default way to set density floors
int set_density_floors_default(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor, FTYPE *prceiling)
{
  struct of_state q;
  FTYPE U[NPR];
  FTYPE bsq;

  if(rescaletype==2){
    // to best conserve E and L along magnetic field lines
    if (get_state(pr, ptrgeom, &q) >= 1)
      FAILSTATEMENT("fixup.c:set_density_floors()", "get_state() dir=0", 1);
    
    MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, U, NULL),"fixup.c:set_density_floors()", "primtoU()", 1);
  }

  if(rescaletype==4||rescaletype==5){
    if(bsq_calc(pr,ptrgeom,&bsq)>=1){
      dualfprintf(fail_file,"bsq_calc:bsq_calc: failure\n");
      return(1);
    }
  }

  set_density_floors_default_alt(ptrgeom, &q, pr, U, bsq, prfloor, prceiling);

  return(0);
}

/// alternative default way to set density floors
int set_density_floors_default_alt(struct of_geom *ptrgeom, struct of_state *q, FTYPE *pr, FTYPE *U, FTYPE bsq, FTYPE *prfloor, FTYPE *prceiling)
{
  FTYPE Upbinf[NPR];
  FTYPE r,th,X[NDIM],V[NDIM];
  int pl,pliter;
  int i=ptrgeom->i,j=ptrgeom->j,k=ptrgeom->k,p=ptrgeom->p;

  coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
  bl_coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, V);
  r=V[1];
  th=V[2];

  // default
  prceiling[RHO]=BIG;
  prceiling[UU]=BIG;


  if(PRAD0>=0.0) prfloor[PRAD0]=ERADLIMIT;

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
      



      // now have U[UU] and U[PH]
      // note that U[UU]/(gdet*rho*u^t) is conserved along field lines, even in non-ideal MHD, where A_{\phi} is still a good stream function in non-ideal MHD.
      // same in terms of U[PH]
      // this is really a coordinate dependent quantity, but if field lines connect to infinity, then this is the same at infinity.

      // Conserved quantity per baryon (or rest-mass unit)
      //      PALLLOOP(k) Upbinf=U[k]/(ptrgeom->gdet * pr[RHO]*(q.ucon[TT]));
      //Upbinf[UU]*=-1.0; // time component has - sign in this -+++ signature code

      // could inject mass or enthalpy, we choose mass for now since rho>h typically (except at very edge of torus)
      prfloor[RHO]=-U[UU]/(ptrgeom->gdet * GAMMAMAX * (q->ucon[TT]));
      prfloor[UU]=prfloor[RHO]*0.01; // cold injection

    }
    else if(rescaletype==3){
      // for dipole
      prfloor[RHO] = prfloorcoef[RHO]*pow(r/Rin, -5);
      prfloor[UU] = prfloorcoef[UU]*pow(r/Rin,-6) ;
    }
    else if(rescaletype==4||rescaletype==5){
      // for jet injection with maximum b^2/rho and b^2/u and maximum u/rho
      
      // below uses max of present u and floor u since present u may be too small (or negative!) and then density comparison isn't consistent with final floor between u and rho

      if(0 && RESTARTMODE==1){ // not on by default  // choose
        // temporarily used to restart if cleaned field or increased resolution of field, because in high b^2/rho or b^2/u regions, the field dissipation due to small-scale structure in the field leads to heating that can become relativistic and disrupt the solution
        // Then also temporarily set UORHOLIMIT to 1.0 or some limiting thing
        FTYPE UORHOLIMITTEMP=0.5; // choose
        FTYPE TRESTART=5930; // choose
        FTYPE TRESTARTE=5921; // choose
        FTYPE BSQORHOLIMITTEMP=BSQORHOLIMIT*2.0;
        int POLESIZEE=8;
        int POLESIZE=3;
        if(t<TRESTARTE && ISSPCMCOORD(MCOORD) && (startpos[2]+j<POLESIZEE || startpos[2]+j>totalsize[2]-1-POLESIZEE)){
          // temporarily allow much higher field due to very early adjustments in field near pole.  This avoids dumping mass into the pole at learly times.
          // But don't let u_g respond to field at all
          prfloor[RHO]=MAX(bsq/BSQORHOLIMITTEMP,SMALL);
          prfloor[UU]=MAX(MIN(pr[UU],UORHOLIMITTEMP*prfloor[RHO]),SMALL);
        }
        else if(t<TRESTART){
          //          if(t<TRESTARTE && ISSPCMCOORD(MCOORD) && (startpos[2]+j<=1 || startpos[2]+j>=totalsize[2]-2)){
          //            // don't allow floor injection of mass density related to u/rho during very early phase
          //            prfloor[RHO]=MAX(bsq/BSQORHOLIMIT);
          //            prfloor[UU]=MIN(pr[UU],UORHOLIMITTEMP*prfloor[RHO]);
          //          }
          //          else{
          //          }
          // don't allow floor injection of mass density related to u/rho during very early phase
          prfloor[RHO]=MAX(bsq/BSQORHOLIMIT,SMALL);
          prfloor[UU]=MAX(MIN(pr[UU],UORHOLIMITTEMP*prfloor[RHO]),SMALL);
        }
        else if(ISSPCMCOORD(MCOORD) && (startpos[2]+j<POLESIZE || startpos[2]+j>totalsize[2]-1-POLESIZE)){
          prfloor[RHO]=MAX(bsq/BSQORHOLIMIT,SMALL);
          // still limit u near pole
          prfloor[UU]=MAX(MIN(pr[UU],UORHOLIMITTEMP*prfloor[RHO]),SMALL);
          //prceiling[UU]=MAX(prfloor[UU],MIN(prceiling[UU],UORHOLIMIT*MAX(pr[RHO],prfloor[RHO])));
        }
        else{ // NORMAL CASE when restarting,etc.
          prfloor[RHO]=MAX(bsq/BSQORHOLIMIT,SMALL);
          prfloor[UU]=MAX(bsq/BSQOULIMIT,zerouuperbaryon*MAX(prfloor[RHO],SMALL));
          //          prfloor[UU]=MAX(MIN(MAX(pr[UU],prfloor[UU]),UORHOLIMIT*prfloor[RHO]),SMALL);
          prceiling[UU]=MAX(prfloor[UU],MIN(prceiling[UU],UORHOLIMIT*MAX(pr[RHO],prfloor[RHO])));

        }// end NORMAL else conditional
      }
      else{ // FULLY NORMAL CASE
        prfloor[RHO]=MAX(bsq/BSQORHOLIMIT,SMALL);
        prfloor[UU]=MAX(bsq/BSQOULIMIT,zerouuperbaryon*MAX(prfloor[RHO],SMALL));
        //        prfloor[UU]=MAX(MIN(MAX(pr[UU],prfloor[UU]),UORHOLIMIT*prfloor[RHO]),SMALL);
        prceiling[UU]=MAX(prfloor[UU],MIN(prceiling[UU],UORHOLIMIT*MAX(pr[RHO],prfloor[RHO])));

        if(rescaletype==5){//WALD
          if(1){
            // do nothing, so constant bsq/rho
          }
          if(0){
            if(r>500.0){
              prfloor[RHO]*=pow(500,-1.5);
              prfloor[UU]*=pow(500,-2.5);
            }
            else{
              prfloor[RHO]*=(pow(r,-1.5)+pow(500,-1.5));
              prfloor[UU]*=(pow(r,-2.5)+pow(500,-2.5));
            }
          }
          if(0){
            FTYPE R=r*sin(th);
            FTYPE Z=fabs(r*cos(th));
            if(Z>500.0){
              prfloor[RHO]*=pow(500,-1.5);
              prfloor[UU]*=pow(500,-2.5);
            }
            else{
              prfloor[RHO]*=pow(Z,-1.5);
              prfloor[UU]*=pow(Z,-2.5);
            }
          }
          if(0){ // WALD TEST
            FTYPE R=r*sin(th);
            FTYPE Z=fabs(r*cos(th));
            prfloor[RHO]*=MIN(1.0/(pow(Z,1.5)*pow(R,-1.5)),1.0);
            prfloor[UU]*=MIN(1.0/(pow(Z,2.5)*pow(R,-2.5)),1.0);
          }

        }// end rescapetype==5
      }// else NORMAL case



    }// end rescaletype==4 || rescaletype==5
  }
  else{
    prfloor[RHO] = RHOMIN*pow(r, -2.0);
    prfloor[UU] = UUMIN*prfloor[RHO] ;
  }
  
  // KORALTODO: Maybe implement floor to ERAD relative to other energy densities like how handle UU and BSQ relative to RHO.






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

/// get flags that indicate properties of flow, like bsq/rho
int get_bsqflags(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE bsq;
  struct of_state q;
  FTYPE gamma,X[NDIM],V[NDIM];
  FTYPE qsq;
  FTYPE prfloor[NPR];
  PFTYPE flags[NUMBSQFLAGS]={0};
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

    // KORALTODO: Maybe implement floor to ERAD relative to other energy densities like how handle UU and BSQ relative to RHO.


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
    set_density_floors(ptrgeom,MAC(pv,i,j,k),prfloor,prceiling);
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
    if((flags[0]==2)||(flags[1]==2)||(flags[2]==2)||(flags[3]==2)||(flags[4]==2)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod[RHO]-1;}
    else if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){ GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod[RHO]; }
    else{GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)=fluxmethod[RHO];}
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

/// limit \gamma=\alpha u^t and u^t
int limit_gamma(int docorrectucons, FTYPE gammamax, FTYPE gammamaxrad, FTYPE*pr, FTYPE *ucons, struct of_geom *ptrgeom,int finalstep)
{
  FTYPE f,gamma,pref;
  FTYPE qsq;
  FTYPE alpha;
  FTYPE radgamma,radqsq;
  FTYPE pr0[NPR];
  int pl,pliter;
  FTYPE realgammamax,realgammamaxrad;
  FTYPE r,X[NDIM],V[NDIM];
  int didchange;
  FTYPE uu0,uu0max;
  int i=ptrgeom->i;
  int j=ptrgeom->j;
  int k=ptrgeom->k;
  int loc=ptrgeom->p;

  //  if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
  //    dualfprintf(fail_file,"inside limit_gamma(): finalstep=%d nstep=%ld steppart=%d\n",finalstep,nstep,steppart);
  //  }

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
  realgammamaxrad=GAMMAMAXRAD;
#else
  realgammamax=gammamax;
  realgammamaxrad=gammamaxrad;
  if(realgammamax<=1.0 && realgammamaxrad<=1.0) return(0); // nothing to do
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

  if(EOMRADTYPE!=EOMRADNONE){
    if(gamma_calc(&pr[URAD1-U1],ptrgeom,&radgamma,&radqsq)>=1){
      dualfprintf(fail_file,"limit_gamma: gamma calc failed: rad\n");
      dualfprintf(fail_file,"i=%d j=%d k=%d : rad : oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,gamma);
      if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
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



  // KORAL:
  if(EOMRADTYPE!=EOMRADNONE){

    if((radgamma > realgammamaxrad && (radgamma!=1.0))) {    

      // rescale velocities to reduce radgamma to realgammamaxrad
      pref=(realgammamaxrad*realgammamaxrad - 1.)/(radgamma*radgamma - 1.);

      if(debugfail>=3){ // special >=3 since often called when floor or fixup
        dualfprintf(fail_file,"nstep=%ld steppart=%d t=%21.15g :: i=%d j=%d k=%d :: pref=%21.15g oldradgamma=%21.15g realgammamaxrad=%21.15g\n",nstep,steppart,t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pref,radgamma,realgammamaxrad);
      }

      if(pref<0.0){
        dualfprintf(fail_file,"limit_gamma: pref calc failed pref=%21.15g\n",pref);
        dualfprintf(fail_file,"i=%d j=%d k=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,radgamma);
        if (fail(i,j,k,loc,FAIL_UTCALC_DISCR) >= 1)
          return (1);
      }

      f = sqrt(pref);
      pr[URAD1] *= f ; 
      pr[URAD2] *= f ; 
      pr[URAD3] *= f ;


      radgamma=gammamaxrad; // reset radgamma for next check
      didchange=1; // indicate did change primitive
    }
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
    int doingmhdfixup=1;
    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"LIMITGAMMABEFORE: steppart=%d nstep=%ld\n",steppart,nstep);
    //    }
    diag_fixup(docorrectucons,prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTLIMITGAMMAACT);
    PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    //    if((startpos[1]+ptrgeom->i==17) && (startpos[2]+ptrgeom->j)==0){
    //      dualfprintf(fail_file,"LIMITGAMMAAFTER: steppart=%d nstep=%ld\n",steppart,nstep);
    //    }

    return(-1);// indicates did change primitive
  }


  return(0); // indicates didn't change anything
}







// check_pr() is purely 2D function

#define UTLIMIT (50.0) // limit of u^t given no other info and need to fix u
#define UTFIX (utobs) //  if failed u^t, then fix to this level to keep normal, or use local value as desired value
#define GRADIENTFACTOR (.001) // how to set?

/// pr is pr to be fixed
/// prmodel is pr that has a model ucon[TT], since we want to interpolate based upon ucon[TT], not just a random one when inside step_ch.c in fluxcalc()
/// prmodel is just pr if in utoprim.c since we want a real change, not just an interpolation, thus will not be used as basis
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
        if(realdiscrlimit-uttdiscr>0) pr[U1+i-1]-=gradient[i]*dampfactor;
        else  pr[U1+i-1]+=gradient[i]*dampfactor;
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
  int doingmhdfixup=1;
  int docorrectucons=(DOENOFLUX != NOENOFLUX);
  diag_fixup(docorrectucons,prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTLIMITGAMMAACT);
  PLOOP(pliter,pl) prdiag[pl]=pr[pl];

#endif

  return(0);
}


/// Check for inflow and modify flow if didn't want inflow (for 4vel)
/// GODMARK: check this function for correctness
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
       ((startpos[1]+ii<=iin)&&((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW))&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW))&&(pr[U1+dir-1] < 0.)) 
        ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[U1]=0;
      
      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[1]+ii<=iin)&&((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW))&&(pr[URAD1+dir-1] > 0.)) 
                                ||((startpos[1]+ii>=iout)&&((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW))&&(pr[URAD1+dir-1] < 0.)) 
                                )
       ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[URAD1]=0;
      
      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
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
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
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
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW)&&(pr[URAD1+dir-1] > 0.)) 
                                ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW)&&(pr[URAD1+dir-1] < 0.)) 
                                )
       ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[URAD2]=0;

      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
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
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[3]+kk<=kkn)&&(BCtype[X3DN]==OUTFLOW)&&(pr[URAD1+dir-1] > 0.)) 
                                ||((startpos[3]+kk>=kout)&&(BCtype[X3UP]==OUTFLOW)&&(pr[URAD1+dir-1] < 0.)) 
                                )
       ) {
      // set pre-primitive
      PALLLOOP(pl)    pr0[pl]=pr[pl];
      pr[URAD3]=0;

      // account for changes
      PLOOP(pliter,pl) prdiag[pl]=pr0[pl];
      int doingmhdfixup=1;
      diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
      PLOOP(pliter,pl) prdiag[pl]=pr[pl];
    }
  }
  else return(1); // uh

  return(0);
}


/// Check for inflow and modify flow if didn't want inflow (for 3vel)
int inflow_check_3vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{

  return(inflow_check_4vel(dir,pr,ucons, ptrgeom,-1));

}

/// Check for inflow and modify flow if didn't want inflow (for rel4vel)
// GODMARK: check for correctness
int inflow_check_rel4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep)
{
  FTYPE ucon[NDIM],uradcon[NDIM] ;
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE othersrad[NUMOTHERSTATERESULTS];
  int ii,jj,kk,loc;
  int j,k ;
  int pl,pliter;
  FTYPE alpha,betacon,gamma,vsq ;
  FTYPE qsq;
  int iin,iout;
  int jjn,jout;
  int kkn,kout;
  int dofix=0, dofixrad=0;
  FTYPE pr0[NPR];


  //  return(0);
  ii=ptrgeom->i; 
  jj=ptrgeom->j;
  kk=ptrgeom->k;
  loc=ptrgeom->p;


  //ucon_calc(pr, ptrgeom, ucon,others) ;

  //  PALLLOOP(pl) dualfprintf(fail_file,"nstep=%ld steppart=%d t=%21.15g :: pl=%d %21.15g\n",nstep,steppart,t,pl,pr[pl]);

  MYFUN(ucon_calc(pr, ptrgeom, ucon, others),"fixup.c:inflow_check_rel4vel()", "ucon_calc() dir=0", 1);
  MYFUN(ucon_calc(&pr[URAD1-U1], ptrgeom, uradcon, othersrad),"fixup.c:inflow_check_rel4vel()", "ucon_calc() dir=0", 1);

  //  DLOOPA(pl) dualfprintf(fail_file,"ucon[%d]=%21.15g uconrad[%d]=%21.15g\n",pl,ucon[pl],pl,uradcon[pl]);

  dofix=dofixrad=0;
 
  if(dir==1){
    // for dir==1
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)||(ptrgeom->p==FACE3)||(ptrgeom->p==CORN1) ){ iin=-1; iout=totalsize[1]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN2)||(ptrgeom->p==CORN3) ){ iin=0; iout=totalsize[1]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[1]+ii<=iin)&&((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW) || BCtype[X1DN]==OUTFLOWNOINFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[1]+ii>=iout)&&((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW) || BCtype[X1UP]==OUTFLOWNOINFLOW)&&(ucon[dir] < 0.)) 
        ) {
      dofix=1;
    }
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[1]+ii<=iin)&&((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW) || BCtype[X1DN]==OUTFLOWNOINFLOW)&&(uradcon[dir] > 0.)) 
                                ||((startpos[1]+ii>=iout)&&((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW) || BCtype[X1UP]==OUTFLOWNOINFLOW)&&(uradcon[dir] < 0.)) 
                                )
       ) {
      dofixrad=1;
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
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW || BCtype[X2DN]==OUTFLOWNOINFLOW)&&(uradcon[dir] > 0.)) 
                                ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW || BCtype[X2UP]==OUTFLOWNOINFLOW)&&(uradcon[dir] < 0.)) 
                                )
       ) {
      dofixrad=2;
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
    if(EOMRADTYPE!=EOMRADNONE&&(
                                ((startpos[3]+kk<=kkn)&&(BCtype[X3DN]==OUTFLOW || BCtype[X3DN]==OUTFLOWNOINFLOW)&&(uradcon[dir] > 0.)) 
                                ||((startpos[3]+kk>=kout)&&(BCtype[X3UP]==OUTFLOW || BCtype[X3UP]==OUTFLOWNOINFLOW)&&(uradcon[dir] < 0.)) 
                                )
       ) {
      dofixrad=3;
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
    int doingmhdfixup=1;
    diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
    PLOOP(pliter,pl) prdiag[pl]=pr[pl];

  }



  if(dofixrad){

    // set pre-primitive
    PALLLOOP(pl)    pr0[pl]=pr[pl];
    FTYPE prdiag[NPR];
    PLOOP(pliter,pl) prdiag[pl]=pr0[pl];


    /* find gamma and remove it from primitives */
    if(gamma_calc(&pr[URAD1-U1],ptrgeom,&gamma,&qsq)>=1){
      dualfprintf(fail_file,"gamma calc failed: inflow_check_rel4vel\n");
      if (fail(ii,jj,kk,loc,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    pr[URAD1] /= gamma ;
    pr[URAD2] /= gamma ;
    pr[URAD3] /= gamma ;
    alpha = 1./sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    
    /* reset radial velocity so radial 4-velocity
     * is zero */
    if(dofixrad==1){
      betacon = ptrgeom->gcon[GIND(0,1)]*alpha*alpha ;
      pr[URAD1] = betacon/alpha ; // gives 3-vel contravariant
    }
    else if(dofixrad==2){
      betacon = ptrgeom->gcon[GIND(0,2)]*alpha*alpha ;
      pr[URAD2] = betacon/alpha ; // gives 3-vel contravariant
    }
    else if(dofixrad==3){
      betacon = ptrgeom->gcon[GIND(0,3)]*alpha*alpha ;
      pr[URAD3] = betacon/alpha ; // gives 3-vel contravariant
    }
    /* now find new gamma and put it back in */
    vsq = 0. ;
    SLOOP(j,k) vsq += ptrgeom->gcov[GIND(j,k)]*pr[URAD1+j-1]*pr[URAD1+k-1] ;

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
    pr[URAD1] *= gamma ;
    pr[URAD2] *= gamma ;
    pr[URAD3] *= gamma ;

    // only for boundary conditions, not active zones, hence -1.0 instead of finalstep
    int doingmhdfixup=1;
    diag_fixup((DOENOFLUX != NOENOFLUX),prdiag, pr, ucons, ptrgeom, finalstep,doingmhdfixup,COUNTINFLOWACT);
    PLOOP(pliter,pl) prdiag[pl]=pr[pl];

  }
  

  if(dofix||dofixrad){
    return(-1);
  }

  return(0);

}
 

// fixup fluxes
/// only correct for polar axis at both inner and outer x2/theta grid edge.
/// OPENMPOPTMARK: Could optimize this, but not using it currently
void fix_flux(FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
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
    if((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW)){
      LOOPX1dir{
        ri=riin;
        LOOPBOUND1IN{
          if(MACP0A1(F1,i,j,k,RHO)>0) MACP0A1(F1,i,j,k,RHO)=0;
        }
      }
    }
  }
  if(mycpupos[1]==ncpux1-1){
    if((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW)){
      LOOPX1dir{
        LOOPBOUND1OUT{
          if(MACP0A1(F1,i,j,k,RHO)<0) MACP0A1(F1,i,j,k,RHO)=0;
        }
      }
    }
  }

}


/// counter for EOS table lookup failures
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
