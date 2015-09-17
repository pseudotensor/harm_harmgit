#include "decs.h"

/*! \file advance.c
  \brief Takes RK sub-step

  Takes RK substep by doing 1) flux and source 2) flux+source->dU 3) dU->U 4) U->P
*/



/// static declarations
static int compute_dt_fromsource(struct of_geom *ptrgeom, struct of_state *state, FTYPE *U, FTYPE *pr, FTYPE *dUevolve, FTYPE *dUgeomevolveUU, FTYPE *dtij, FTYPE *gravitydt);
static int dUtodt(struct of_geom *ptrgeom, struct of_state *state, FTYPE *pr, FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *dUgeomgravity, FTYPE *accdt, FTYPE *gravitydt);
static int check_point_vs_average(int timeorder, int numtimeorders, PFTYPE *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats);



static int prepare_globaldt(
                            int truestep, 
                            FTYPE ndt1,FTYPE ndt2,FTYPE ndt3,
                            FTYPE accdt,int accdti,int accdtj,int accdtk,
                            FTYPE gravitydt,int gravitydti,int gravitydtj,int gravitydtk,
                            FTYPE *ndt);


//FTYPE globalgeom[NPR];
static void flux2dUavg(int whichpl, int i, int j, int k, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE *dUavg1,FTYPE *dUavg2,FTYPE *dUavg3);
static void dUtoU(int timeorder, int whichpl, int i, int j, int k, int loc, FTYPE *dUgeom, FTYPE (*dUcomp)[NPR], FTYPE *dUriemann, FTYPE *CUf, FTYPE *CUnew, FTYPE *Ui,  FTYPE *uf, FTYPE *ucum);
static void dUtoU_check(int timeorder, int i, int j, int k, int loc, int pl, FTYPE *dUgeom, FTYPE (*dUcomp)[NPR], FTYPE *dUriemann, FTYPE *CUf, FTYPE *CUnew, FTYPE *Ui,  FTYPE *Uf, FTYPE *ucum);
static int asym_compute_1(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
static int asym_compute_2(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);


static FTYPE fractional_diff( FTYPE a, FTYPE b );

static FTYPE compute_dissmeasure(int timeorder, int i, int j, int k, int loc, FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *CUf, FTYPE *CUnew, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE *Ui,  FTYPE *Uf, FTYPE *tempucum);


// AVG_2_POINT functions:
static void debugeno_compute(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE (*U)[NSTORE2][NSTORE3][NPR],FTYPE (*debugvars)[NSTORE2][NSTORE3][NUMENODEBUGS]);

static void show_fluxes(int i, int j, int k, int loc, int pl,FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);


static int advance_standard(int truestep,int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                            FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                            FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                            FTYPE (*ui)[NSTORE2][NSTORE3][NPR], FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
                            FTYPE *CUf,FTYPE *CUnew,SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime, int stagenow, int numstages, FTYPE *ndt);
static int advance_standard_orig(int truestep,int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                            FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                            FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                            FTYPE (*ui)[NSTORE2][NSTORE3][NPR], FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
                            FTYPE *CUf,FTYPE *CUnew,SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime, int stagenow, int numstages, FTYPE *ndt);
static int advance_finitevolume(int truestep,int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                                FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                                FTYPE (*ui)[NSTORE2][NSTORE3][NPR],FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
                                FTYPE *CUf,FTYPE *CUnew, SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime,  int stagenow, int numstages, FTYPE *ndt);




/// things to do before any interpolation or advance step
/// includes pre-computed things for interpolation and advance that (e.g.) aren't required to perform for each interpolation or advance call or portion of a call
void pre_interpolate_and_advance(FTYPE (*pb)[NSTORE2][NSTORE3][NPR])
{

  /////////////////////////////////////
  //
  // Compute and Store (globally) the get_state() data for the CENT position to avoid computing later and for merged-higher-order method
  //
  /////////////////////////////////////
#if(STOREFLUXSTATE||STORESHOCKINDICATOR)
  // NOTE: This is done before advance since always needed, and only needed once for all dimensions, and don't here instead of inside advance() since needed during fluxcalc() that is called first before any use of get_geometry() that we would use to put this call with
  compute_and_store_fluxstatecent(pb);
  // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif


  // indicate to any needed piece of code that computed pre-interpolates if needed.
  global_preinterpolate=1;

}


/// pi: initial values at t=t0 to compute Ui
/// pb: values used to compute flux/source
/// pf: solution using flux(pb) from pi's Ui -> Uf
/// pi, pb, and pf can all be the same since
/// 1) pb used first on a stencil, not modified, to compute fluxes
/// 2) pf=pi is assigned by value at each zone
/// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
///
/// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance(int truestep, int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
            FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
            FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
            FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
            FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
            FTYPE (*ui)[NSTORE2][NSTORE3][NPR],FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
            FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime, int timeorder, int numtimeorders, FTYPE *ndt)
{


  ////////////////
  //
  // setup state and interpolation stuff for interpolation and advance calls
  //
  ////////////////
  pre_interpolate_and_advance(pb);



  ////////////////
  //
  // advance
  //
  ////////////////
  if(DOENOFLUX==ENOFINITEVOLUME){
    MYFUN(advance_finitevolume(truestep,stage,pi,pb,pf,pstag,pl_ct, pr_ct, F1, F2, F3, vpot,ui,uf,ucum,CUf,CUnew,fluxdt,boundtime,fluxtime,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_finitevolume()", 1);
  }
  else if((DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFLUXSPLIT)){
    // new standard preserves conserved quantities even when metric changes
    //    MYFUN(advance_standard_orig(truestep,stage,pi,pb,pf,pstag,pl_ct, pr_ct, F1, F2, F3, vpot,ui,uf,ucum,CUf,CUnew,fluxdt,boundtime,fluxtime,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_standard()", 1);
    MYFUN(advance_standard(truestep,stage,pi,pb,pf,pstag,pl_ct, pr_ct, F1, F2, F3, vpot,ui,uf,ucum,CUf,CUnew,fluxdt,boundtime,fluxtime,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_standard()", 1);
  }
  else{
    dualfprintf(fail_file,"No such DOENOFLUX=%d\n",DOENOFLUX);
    myexit(1);
  }


  return(0);


}



// Notes:
// loop range defined with +SHIFT? so consistent with requirement by IF3DSPCTHENMPITRANSFERATPOLE or consistent with setting upper face flux and using it, which is only used for FLUXBSTAG for evolving field on faces.  This forces field on upper face to be evolved as required in some cases.
// This is bit excessive for non-face quantities, but fake partial evolution of centered value at "N" is ok as long as don't hit NaN's that slow down code.





/// this method guarantees conservation of non-sourced conserved quantities when metric is time-dependent
/// this method has updated field staggered method
/// Note that when dt==0.0, assume no fluxing, just take ucum -> ui -> {uf,ucum} and invert.  Used with metric update.
///
/// NEW: like advance_standard_orig(), but removed debug info and set field "inversion" first so have centered value for source() so have it for any point-use of values like in implicit solver for radiation-fluid interaction.
static int advance_standard(
                            int truestep,
                            int stage,
                            FTYPE (*pi)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pb)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                            FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                            FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                            FTYPE (*ui)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*uf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*ucum)[NSTORE2][NSTORE3][NPR], 
                            FTYPE *CUf, 
                            FTYPE *CUnew, 
                            SFTYPE fluxdt,
                            SFTYPE boundtime,
                            SFTYPE fluxtime,
                            int timeorder, int numtimeorders,
                            FTYPE *ndt)
{
  FTYPE ndt1, ndt2, ndt3;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  static SFTYPE dt4diag_willbe=0;
  int finalstep,initialstep;
  FTYPE accdt, accdt_ij;
  int accdti,accdtj,accdtk;
  FTYPE gravitydt, gravitydt_ij;
  int gravitydti,gravitydtj,gravitydtk;
  //  FTYPE (*dUriemannarray)[NSTORE2][NSTORE3][NPR];
  FTYPE (*ucumformetric)[NSTORE2][NSTORE3][NPR];
  int enerregion;
  int *localenerpos;
  int jj;
  FTYPE (*utoinvert)[NSTORE2][NSTORE3][NPR];
  FTYPE *utoinvert1;
  FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR];
  FTYPE (*useducum)[NSTORE2][NSTORE3][NPR];
  FTYPE *useducum1;
  FTYPE (*preupoint)[NSTORE2][NSTORE3][NPR];
  FTYPE (*myupoint)[NSTORE2][NSTORE3][NPR];
  FTYPE *myupoint1;
  FTYPE (*myupointuf)[NSTORE2][NSTORE3][NPR];
  FTYPE (*olduf)[NSTORE2][NSTORE3][NPR];
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int is,ie,js,je,ks,ke;
  int doingextrashiftforstag;




  if((GLOBALBCMOVEDWITHACTIVESECTION==0 || USEMPI) && truestep && DOGRIDSECTIONING){
    // then need to fill pf with pi everywhere in global region, so that non-ACTIVE region will be set for boundary conditions within the domain.
    int i,j,k,pliter,pl;
    LOOPF{
      PLOOP(pliter,pl) MACP0A1(pf,i,j,k,pl) = MACP0A1(pi,i,j,k,pl);
    }
  }


  if(FLUXB==FLUXCTSTAG){
    // fill in tempucum with changes that are not controlled well in space, but later fill in ucum from this in controlled way
    tempucum=GLOBALPOINT(utemparray);
    useducum=tempucum; // unless changed
  }
  else{
    // no special shifting of indices occurs that requires loop shifting
    tempucum=ucum;
    useducum=ucum;
  }
  olduf=GLOBALPOINT(oldufstore);


  ucumformetric=GLOBALPOINT(ucumformetric);// temporary space for ucum for metric that is same time as "pb", so not updated yet or is ui



  /////////////////////////////////////////////
  //
  // Setup function tasks
  //
  ////////////////////////////////////////////


  // for ZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];


  accdt=BIG; // initially no limit to dt due to acceleration
  accdti=accdtj=accdtk=-100;
  gravitydt=BIG; // initially no limit to dt due to time derivatives in gravity
  gravitydti=gravitydtj=gravitydtk=-100;





  /////////
  //
  // set initialstep and finalstep to tell some procedures and diagnostic functions if should be accounting or not
  //
  /////////
  if(timeorder==0){
    initialstep=1;
  }
  else{
    initialstep=0;
  }

  if(timeorder==numtimeorders-1){
    finalstep=1;
  }
  else{
    finalstep=0;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  // whether need to compute flux
  // only compute flux if flux is added to get Uf or Ucum
  ////////////////////////////////////////////////////////////////////////////////
  int doflux = (CUf[2]!=0.0 || CUnew[1]!=0.0); // in reality, only do flux if CUf[2]!=0 as CUnew is just to accumulate to get U^{n+1}
  // current implicit dt factor
  // Since we don't pass CUnew[NUMPREDTCUFS+timeorder], we assume never try to add M^i into U^{n+1} unless computed for current U^i or previous U^i
  FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];


  /////////
  //
  // set dt4diag for source diagnostics
  //
  /////////
  if(timeorder==numtimeorders-1 && (nstep%DIAGSOURCECOMPSTEP==0) ){
    // every 4 full steps as well as on final step of full step (otherwise diag_source_comp() too expensive)
    dt4diag=dt4diag_willbe;
    dt4diag_willbe=0;
  }
  else{
    dt4diag_willbe+=dt;
    dt4diag=-1.0;
  }



  /////////////////////////////////////////////
  //
  // Setup loops [+1 extra compared to normal comp region if doing FLUXCTSTAG]
  //
  ////////////////////////////////////////////
  get_flux_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);



  /////////////////////////////////////////////
  //
  // Set initial ui, temporary space, so ok that COMPZLOOP() goes over +1 in FLUXB==FLUXCTSTAG case
  //
  ////////////////////////////////////////////
  if(timeorder==0){
    // last timestep's final ucum is stored into ucumformetric and ui as initial term in ucum
    // copy ucum -> {ucumformetric,ui}
    if(doingmetricsubstep()) copy_3dnpr_2ptrs(is,ie,js,je,ks,ke,ucum,ucumformetric,ui);
    else copy_3dnpr(is,ie,js,je,ks,ke,ucum,ui); // only need to setup ui then
  }
  else{
    // preserve this time's value of ucum for the metric (for timeorder==0 ucumformetric is assigned from ui)
    // copy ucum -> ucumformetric
    if(doingmetricsubstep()) copy_3dnpr(is,ie,js,je,ks,ke,ucum,ucumformetric);
  }



  /////////////////////////////////////////////
  //
  // Compute flux
  //
  ////////////////////////////////////////////

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif




  ndt1=ndt2=ndt3=BIG; // if doflux==0 or truestep==0, then setting ndt=BIG means doesn't affect any other sub-steps or steps that are trying to set dt as minimum over sub-step dt's.
  if(doflux && truestep){ // only do if not just passing through
    if(1){
      // NORMAL:
      // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
      MYFUN(fluxcalc(stage,initialstep,finalstep,pb,pstag,pl_ct, pr_ct, vpot,F1,F2,F3,CUf,CUnew,fluxdt,fluxtime,&ndt1,&ndt2,&ndt3),"advance.c:advance_standard()", "fluxcalcall", 1);
    }
  }// end if not just passing through
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil



#if(PRODUCTION==0)
  trifprintf( "1f");
#endif

 




  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //
  // now update get flux update [loop should go over normal computational region +1 on outer edges so automatically set field if staggered.  Since we only set tempucum so far, ucum in that +1 region is unaffected]
  //
  /////////////////////////////////////////////////////


  int didreturnpf;
  if(truestep){


    // initialize uf and ucum if very first time here since ucum is cumulative (+=) [now tempucum is cumulative]
    // copy 0 -> {uf,tempucum,olduf}
    if(timeorder==0) init_3dnpr_3ptrs(is, ie, js, je, ks, ke,0.0, uf,tempucum,olduf);



#if(WHICHEOM==WITHNOGDET && (NOGDETB1==1 || NOGDETB2==1 || NOGDETB3==1) )
    // for FLUXB==FLUXCTSTAG, assume no source term for field
    if(FLUXB==FLUXCTSTAG){
      dualfprintf(fail_file,"Not setup for field source term if staggered field\n");
      myexit(176023);
    }
#endif

  
   

    ////////////////////////
    //
    // FIRST UPDATE FIELD using its flux
    //
    ////////////////////////
    int doother=DOALLPL; // default
    if(FLUXB==FLUXCTSTAG){
      doother=DONONBPL;
      
      // then field version of ui[] is stored as "conserved" value at FACE, not CENT
      //PLOOPBONLY(pl) MACP0A1(ui,i,j,k,pl) is itself // FACE (see step_ch.c's setup_rktimestep and know that ui=unew for before first substep)

      if(NOGDETB1 ||NOGDETB2 ||NOGDETB3){
        dualfprintf(fail_file,"This approach requires B have no source terms.\n");
        myexit(34983463);
      }

#pragma omp parallel // OPENMPGLOBALPRIVATEFORINVERSION // don't need EOS stuff here since not getting state or doing inversion yet.
      {
        int pl,pliter,i,j,k;
        // zero-out dUgeom in case non-B pl's used.
        FTYPE dUgeom[NPR]={0.0},dUriemann[NPR+NSPECIAL]={0.0},dUriemann1[NPR+NSPECIAL]={0.0},dUriemann2[NPR+NSPECIAL]={0.0},dUriemann3[NPR+NSPECIAL]={0.0};
        FTYPE dUcomp[NUMSOURCES][NPR];
        int sc;
        PLOOP(pliter,pl) SCLOOP(sc) dUcomp[sc][pl]=0.0;


        /////////////////////////////////////////////////////
        // [Loop goes over up to +1 normal computational region for getting new staggered U if doing FLUXCTSTAG] so can do interpolation on it to get centered field
        /////////////////////////////////////////////////////
        OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);

          // dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
          if(doflux){
            flux2dUavg(DOBPL,i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
            PLOOPBONLY(pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is one type of avg->point mistake
          }
          else PLOOPBONLY(pl) dUriemann[pl]=0.0;

          // Get update
          // ui is itself at FACE as already set
          // this overwrites uf[B1,B2,B3], while need original uf[B1,B2,B3] for some RK methods, so store as olduf outside this loop, already.
          dUtoU(timeorder,DOBPL,i,j,k,CENT,dUgeom, dUcomp, dUriemann, CUf, CUnew, MAC(ui,i,j,k), MAC(uf,i,j,k), MAC(tempucum,i,j,k));
          if(finalstep==0){
            PLOOP(pliter,pl) if(BPL(pl)) MACP0A1(olduf,i,j,k,pl) = MACP0A1(uf,i,j,k,pl);
          }

        }//end loop
      }// end parallel



      // first copy over all quantities as point, which is true except for fields if FLUXRECON active
      // copy utoinvert -> myupoint
      // only copy magnetic field pl's -- later copy rest when needed for inversion
      // KORALNOTE: For general RK method, need to feed-into implicit solver field from Uf, not tempucum.  So seems I need both fields to be processed for such RK methods.   One for actual final field, and one for "initial+flux" form into implicit solver that just contributes to tempucum.

      // if using staggered grid for magnetic field, then need to convert ucum to pstag to new pb/pf
      // GODMARK: Use of globals
      myupoint=GLOBALPOINT(upointglobal);
      
      if(1){ // normal switching case
      
        if(finalstep==1) preupoint=tempucum;
        else preupoint=uf;
      
        copy_tempucum_finalucum(DOBPL,Uconsevolveloop,preupoint,myupoint);

        // uses tempucum and gets reaveraged field into myupoint
        if(extrazones4emf && dofluxreconevolvepointfield==0) field_integrate_fluxrecon(stage, pb, preupoint, myupoint);

        // first pb entry is used for shock indicator, second is filled with new field
        // myupoint goes in as staggered point value of magnetic flux and returned as centered point value of magnetic flux
        interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstag, myupoint, NULL);
      }

      // special higher-time-order uf-always calculation of field for when final ucum is different from final uf
      // if CUnew[2]==1 on final step, then means Uf we compute is same as ucum.
      if(finalstep==1 && CUnew[2]!=1.0){

        // still have to do uf calculation
        myupointuf=GLOBALPOINT(upointglobaluf);
        // get other non-field things in case used
        // NO, not used, so can avoid.
        //        int pliter;
        //        // NOT RIGHT, shoud full copy: PLOOP(pliter,pl) MACP0A1(myupointuf,i,j,k,pl) = MACP0A1(myupoint,i,j,k,pl);

        preupoint=uf;
      
        copy_tempucum_finalucum(DOBPL,Uconsevolveloop,preupoint,myupointuf);

        // uses tempucum and gets reaveraged field into myupointuf
        if(extrazones4emf && dofluxreconevolvepointfield==0) field_integrate_fluxrecon(stage, pb, preupoint, myupointuf);

        // first pb entry is used for shock indicator, second is filled with new field
        // myupointuf goes in as staggered point value of magnetic flux and returned as centered point value of magnetic flux
        interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstag, myupointuf, NULL);

      }
      else{
        // then no need for separate myupointuf, so just point to same space
        myupointuf=myupoint;
      }


      ////////////////////    
      // now myupoint contain CENTered point conserved (to be converted to primitive quantity later) ready for MHD or RAD inversion procedures (or implicit use of such inversions)
      // myupointuf contains always uf version of field
      // myupoint constains uf version except for finalstep=1, when it contains tempucum version
      ////////////////////

    }// end if staggered field method
  




    // get loop range (only needs to be over same location as final centered primitive to invert.
    get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);


    ////////////////////////
    //
    // UPDATE NON_FIELD using flux and source
    // PERFORM INVERSION
    // DO FIXUP1ZONE
    //
    ////////////////////////
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM  // <-- only includes more than OPENMPGLOBALPRIVATEFORINVERSION needed for inversion
    {
      int pl,pliter,i,j,k;
      struct of_geom geomdontuse;
      struct of_geom *ptrgeom=&geomdontuse;

      // for source()
      FTYPE Uitemp[NPR];
      // set all to zero in case doother=DONONBPL in which case need rest of calculations to know no change to field
      FTYPE dUgeom[NPR]={0},dUriemann[NPR]={0},dUriemann1[NPR+NSPECIAL]={0},dUriemann2[NPR+NSPECIAL]={0},dUriemann3[NPR+NSPECIAL]={0},dUcomp[NUMSOURCES][NPR]={{0}};
      struct of_state qdontuse;
      struct of_state *qptr=&qdontuse;
      struct of_state qdontuse2;
      struct of_state *qptr2=&qdontuse2; // different qptr since call normal and special get_state()


      // for inversion
      FTYPE prbefore[NPR];
      struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
      int showmessages=1;
      int allowlocalfailurefixandnoreport=1; // allow local fixups
      // initialize counters
      newtonstats.nstroke=newtonstats.lntries=0;
      int eomtype;



      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
      
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // setup default eomtype
        eomtype=EOMDEFAULT;


        // set geometry for centered zone to be updated
        get_geometry(i, j, k, CENT, ptrgeom);
      
        // find Ui(pi)
        // force use of primitive to set Ui since otherwise wherever corrected/changed primitive (in fixup, etc.) then would have to change conserved quantity, while same since both are point values
        // only field for staggered method is special point value at faces that needs to come from conserved quantity
        MYFUN(get_state(MAC(pi,i,j,k), ptrgeom, qptr),"step_ch.c:advance()", "get_state()", 1);
        MYFUN(primtoU(UEVOLVE,MAC(pi,i,j,k), qptr, ptrgeom, Uitemp, NULL),"step_ch.c:advance()", "primtoU()", 1);

        if(FLUXB==FLUXCTSTAG || DOENOFLUX != NOENOFLUX ){
          // then field version of ui[] is stored as "conserved" value at FACE, not CENT
          PLOOPNOB1(pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // CENT
          //PLOOPBONLY(pl) MACP0A1(ui,i,j,k,pl) is itself // FACE (see step_ch.c's setup_rktimestep and know that ui=unew for before first substep)
          PLOOPNOB2(pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // CENT
        }
        else{
          PLOOP(pliter,pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // all at CENT
        }

        //        if(myid==5 && nstep==1 && steppart==0 && ptrgeom->i==19 && ptrgeom->j==15){
        //          PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pi=%21.15g ui=%21.15g\n",MACP0A1(pi,i,j,k,pl),MACP0A1(ui,i,j,k,pl)*ptrgeom->igdetnosing);
        //        }
 
        // dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
        if(doflux){
          flux2dUavg(doother,i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
          PLOOP(pliter,pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is one type of avg->point mistake
        }
        else{
          PLOOP(pliter,pl) dUriemann[pl]=0.0;
        }

        // save pi before it gets modified in case pf=pi or pf=pb as a pointer.
        FTYPE piorig[NPR],pborig[NPR];
        PALLLOOP(pl){
          piorig[pl] = MACP0A1(pi,i,j,k,pl);
          pborig[pl] = MACP0A1(pb,i,j,k,pl);
        }


        /////////////
        // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
        // MAJORNOTE: final step often sets pointers of pi=pf, in order to use arbitrary guess, so must set pf ONLY once pi is done being used.
        // pb is probably closer than pi for inversion
        ////////////
        PALLLOOP(pl) MACP0A1(pf,i,j,k,pl) = MACP0A1(pb,i,j,k,pl);
        if(FLUXB==FLUXCTSTAG){
          // then upoint actually contains centered updated field used for source() and starting point for more readily getting inversion
          // myupointuf contains always uf version, not tempucum version, of updated field.
          PLOOPBONLY(pl) MACP0A1(pf,i,j,k,pl)=MACP0A1(myupointuf,i,j,k,pl)*ptrgeom->igdetnosing; // valid for this setup of NOGDETB1/B2/B3==0
        }


        /////////////////////////////////////////////////////
        // get state since both source() and dUtodt() need same state
        // From pb, so different than state for Ui(pi)
        // Inside source() might use pb,pf,etc. to get new state.
        MYFUN(get_stateforsource(pborig, ptrgeom, &qptr2) ,"advance.c:()", "get_state() dir=0", 1);
      

        // note that uf and ucum are initialized inside setup_rktimestep() before first substep

        // get dissmeasure
        FTYPE dissmeasure=compute_dissmeasure(timeorder,i,j,k,ptrgeom->p,MAC(pf,i,j,k),ptrgeom,qptr2,CUf, CUnew, F1, F2, F3, MAC(ui,i,j,k),MAC(olduf,i,j,k), MAC(tempucum,i,j,k)); // GODMARK:  If doflux==0, then this actually uses old F's.

        // find dU(pb)
        // so pf contains updated field at cell center for use in (e.g.) implicit solver that uses inversion P(U)
        // Note that uf[B1,B2,B3] is already updated, but need to pass old uf for RK3/RK4, so use olduf.
        MYFUN(source(piorig, pborig, MAC(pf,i,j,k), &didreturnpf, &eomtype, ptrgeom, qptr2, MAC(ui,i,j,k), MAC(olduf,i,j,k), CUf, CUimp,dissmeasure, dUriemann, dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
        // assumes final dUcomp is nonzero and representative of source term over this timestep
        


#if(SPLITNPR)
        // don't update metric if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
            if(DODIAGS){
#if(ACCURATESOURCEDIAG>=1)
              diag_source_all(ptrgeom,dUgeom,fluxdt);
#else
              diag_source_all(ptrgeom,dUgeom,dt4diag);
#endif
#if(ACCURATESOURCEDIAG>=2)
              diag_source_comp(ptrgeom,dUcomp,fluxdt);
#else
              diag_source_comp(ptrgeom,dUcomp,dt4diag);
#endif
            }
          }

      
        /////////////
        //
        // Get update: tempucum
        //
        /////////////
        //        PLOOP(pliter,pl) globalgeom[pl]=ptrgeom->gdet;
        // ok to use "uf" here instead of "olduf" since field is not done here.
        FTYPE ufconsider[NPR],tempucumconsider[NPR];
        PLOOP(pliter,pl){
          ufconsider[pl]=MACP0A1(olduf,i,j,k,pl);
          tempucumconsider[pl]=MACP0A1(tempucum,i,j,k,pl);
        }
        dUtoU(timeorder,doother,i,j,k,ptrgeom->p,dUgeom, dUcomp, dUriemann, CUf, CUnew, MAC(ui,i,j,k), ufconsider, tempucumconsider);

        // KORALNOTE: field dUtoU loses previous time uf needed for RK3 and RK4, so need to store it for safe keeping
        // this oldud is used for B1,B2,B3 and required in cases when CUf[1]!=0.0, as even TVD RK2 method needs because previous Uf used even for updating new Uf.
        ///////
        // Save uf here because below modify uf to account for floors.
        // The later modification of uf is required for >=RK3 methods that use olduf
        FTYPE origolduf[NPR];
        if(doother==DONONBPL){
          PLOOP(pliter,pl){
            if(BPL(pl)==0){
              origolduf[pl] = MACP0A1(olduf,i,j,k,pl);
              MACP0A1(olduf,i,j,k,pl) = ufconsider[pl];
            }
          }
        }
        else if(doother==DOALLPL){
          PLOOP(pliter,pl){
            origolduf[pl] = MACP0A1(olduf,i,j,k,pl);
            MACP0A1(olduf,i,j,k,pl) = ufconsider[pl];
          }
        }


        ////////////////////////////
        // Choose what to invert
        ////////////////////////////
        if(finalstep){
          // invert ucum on final step
          utoinvert = tempucum;
          useducum=tempucum;

          utoinvert1 = tempucumconsider;
          useducum1 = tempucumconsider;

        }
        else{
          // invert uf on substeps
          utoinvert = uf; // should always be uf, not olduf.
          // tempucum just cumulates for now
          useducum=tempucum;

          utoinvert1 = ufconsider;
          useducum1 = tempucumconsider;
        }
        if(FLUXB==FLUXCTSTAG){
          // copy over evolved field.  If finalstep=1, then myupoint contains ucum field as required.  Else has uf field as required.
          PLOOPBONLY(pl) utoinvert1[pl] = MACP0A1(myupoint,i,j,k,pl);
        } // else already there as point centered quantity
        
        

#if(0)
        ////////////////////////////
        // setup myupoint to invert
        ////////////////////////////
        if(FLUXB==FLUXCTSTAG){
          // already have field in myupoint, just copy the others over
          PLOOP(pliter,pl) if(doother==DOALLPL || doother==DONONBPL && BPL(pl)==0 || doother==DOBPL && BPL(pl)==1) MACP0A1(myupoint,i,j,k,pl)=MACP0A1(utoinvert,i,j,k,pl);
        }
        else{
          // utoinvert never reassigned from global a_utoinvert assignment since if here not doing FLUXCTSTAG
          myupoint=utoinvert;
        }
        myupoint1=utoinvert1;
#endif

        ////////////////////////////
        // INVERT [loop only over "centered" cells]
        //
        // Utoprimgen() properly  uses eomtype as potentially changed in source() call
        //
        ////////////////////////////

        /////////////////
        //
        // if doing inversion, then check if should use entropy or energy if originally assuming should use energy
        //
        /////////////////
        int eomtypelocal;
        eomtypelocal=eomtype;



        /////////////////
        //
        // Setup and Do inversion
        //
        /////////////////
        if(finalstep){ // last call, so ucum is cooked and ready to eat!
          // store guess for diss_compute before changed by normal inversion
          PALLLOOP(pl) prbefore[pl]=MACP0A1(pf,i,j,k,pl);

          if(CUnew[2]==1.0){
            // then use original eomtypelocal
          }
          else{
            // then Uf is not ucum on finalstep=1, so any prior implicit inversion was not final inversion.
            // So change any do nothing to do something
            if(eomtypelocal==EOMDIDGRMHD) eomtypelocal=EOMGRMHD;
            else if(eomtypelocal==EOMDIDENTROPYGRMHD) eomtypelocal=EOMENTROPYGRMHD;
            else if(eomtypelocal==EOMDIDCOLDGRMHD) eomtypelocal=EOMCOLDGRMHD;
            else if(eomtypelocal==EOMDIDFFDE) eomtypelocal=EOMFFDE;
            else if(eomtypelocal==EOMDIDFFDE2) eomtypelocal=EOMFFDE2;
          }
        }



        // actual inversion
        int whichcap=CAPTYPEBASIC;
        int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
        int modprim=0;
        int checkoninversiongas=CHECKONINVERSION;
        int checkoninversionrad=CHECKONINVERSIONRAD;
        MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtypelocal,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE,utoinvert1, qptr2, ptrgeom, dissmeasure, piorig, MAC(pf,i,j,k),&newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
        nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;
          
          
#if(DODISS||DODISSVSR)
        if(finalstep){
          // then see what entropy inversion would give
          diss_compute(EVOLVEUTOPRIM,UEVOLVE,utoinvert1,ptrgeom,prbefore,MAC(pf,i,j,k),&newtonstats);
        }
#endif
          
        
        
        ////////////////////////////
        //
        // Do fixup1zone and adjust dUriemann, uf, and tempucum
        //
        ////////////////////////////
#if(SPLITNPR)
        // don't update metric if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
            // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
            // SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
            if((STEPOVERNEGU==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGRHO==NEGDENSITY_ALWAYSFIXUP)||(STEPOVERNEGRHOU==NEGDENSITY_ALWAYSFIXUP)||(finalstep)){
              FTYPE utoinvertlocal[NPR];
              PLOOP(pliter,pl) utoinvertlocal[pl] = utoinvert1[pl];
              int docorrectucons=0;
              int didfixup=fixup1zone(docorrectucons,MAC(pf,i,j,k),utoinvertlocal, ptrgeom,finalstep);
            }// if might do fixups
#endif
          }// end doing single-point fixups


        
        
        /////////////////
        //
        // Save results to uf and tempucum (whether fixup or not)
        //
        //////////////////
        PLOOP(pliter,pl){
          // only save if not already updated uf and tempucum separately
          if(doother==DOALLPL || doother==DONONBPL && BPL(pl)==0 || doother==DOBPL && BPL(pl)==1){
            MACP0A1(uf,i,j,k,pl)=ufconsider[pl];
            MACP0A1(tempucum,i,j,k,pl)=tempucumconsider[pl];
          }
        }

      
            
      
        /////////////////////////////////////
        //
        // get timestep limit from acceleration
        //
        /////////////////////////////////////
#if(LIMITDTWITHSOURCETERM)
#if(SPLITNPR)
        // don't do dUtodt if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
            // geometry is post-metric update, but should still give good estimate of future dt
            dUtodt(ptrgeom, qptr2, MAC(pb,i,j,k), dUgeom, dUriemann, dUcomp[GEOMSOURCE], &accdt_ij, &gravitydt_ij);

#pragma omp critical
            {
              if(accdt_ij<accdt){
                accdt=accdt_ij;
                accdti=i;
                accdtj=j;
                accdtk=k;
              }
              if(gravitydt_ij<gravitydt){
                gravitydt=gravitydt_ij;
                gravitydti=i;
                gravitydtj=j;
                gravitydtk=k;
              }
            }// end critical region
          }// end block that may mean: if not only doing B1-B3

#endif // end if doing LIMITDTWITHSOURCETERM
          



      } // end COMPZLOOP :: end looping to obtain dUriemann and full unew update
    }// end parallel block
  } // end if truestep
  else{
    // then nothing to do since nothing changed
    // previously updated dU and got new ucum as fed into metric, but now metric has its own ucummetric so all that is not needed
    // SUPERGODMARK: no, I guess I don't recall what was being done for metric and why when passing through with dt==0.0

    // just need to copy ui -> {uf,tempucum} to duplicate behavior of dUtoU()
    copy_3dnpr_2ptrs(is,ie,js,je,ks,ke,ui,uf,tempucum);
  }





#if(PRODUCTION==0)
  trifprintf( "#0m");
#endif



 
    
  /////////////////////////////////////////////
  //
  // EVOLVE (update/compute) METRIC HERE
  // In general computes stress-energy tensor (T) from pb and T is then used to compute metric
  // Done here instead of after flux since horizon_flux() updates flux through boundary that would change metric
  // want metric to be at same time as everything else done with pb so RK method is consistent
  //
  // uses unew that's NOT updated yet
  /////////////////////////////////////////////
  if(truestep){
    if(doingmetricsubstep()){
#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
        {
          compute_new_metric_substep(CUf,CUnew,pb,ucumformetric); // CHANGINGMARK : Not sure if CUnew here is correct
        }
    }// end if doing metric substepping
  }





  ////////////////////////////////
  //
  // compute flux diagnostics (accurately using all substeps)
  //
  ///////////////////////////////
  if(doflux && truestep){
    // must come after metric changes that can change where flux surface is since otherwise when flux surface changes, we won't track this substep's flux through the new surface but the old surface (which could even be at r=0 and have no flux)
    // if using unew, then since metric update above uses old unew, need present flux at new horizon surface
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
        diag_flux_pureflux(pb,F1, F2, F3, fluxdt); // fluxdt is true dt for this flux as added in dUtoU() as part of ucum
#endif
      }
  }


  
#if(PRODUCTION==0)
  trifprintf( "#0s");
#endif



  ///////////////////////////////////////
  //
  // Copy over tempucum -> ucum per pl to account for staggered field
  //
  ///////////////////////////////////////
  if(finalstep && FLUXB==FLUXCTSTAG){
    // copy over new ucum in only desired locations irrespective of where tempucum was updated
    // copy tempucum -> ucum
    copy_tempucum_finalucum(DOALLPL,Uconsevolveloop,tempucum,ucum); // fill-in all pl's for storage and for next step.
  }




  /////////////////////////////////
  //
  // If not fixing up primitives after inversion immediately, then fix up all zones at once afterwards
  //
  /////////////////////////////////

#if(SPLITNPR)
  // don't update metric if only doing B1-B3
  if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
#if(FIXUPZONES==FIXUPALLZONES)
      //      fixup(stage,pf,useducum,finalstep); // GODMARK: if want to correct useducum, then have to be more careful like when doing fixup1zone()
      fixup(stage,pf,utoinvert,finalstep);
#endif  
    }



  /////////////////////////////////
  //
  // Determine next timestep from waves, fluxes, and source updates
  //
  /////////////////////////////////


  prepare_globaldt(truestep,ndt1,ndt2,ndt3,accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,ndt);



#if(PRODUCTION==0)
  trifprintf( "2f");
#endif


  return (0);
}










/// compute dissipation measure for determining if can use entropy equations of motion or must use energy equations of motion
static FTYPE compute_dissmeasure(int timeorder, int i, int j, int k, int loc, FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *CUf, FTYPE *CUnew, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE *ui,  FTYPE *uf, FTYPE *tempucum)
{
  FTYPE dissmeasure;
  int pliter,pl;

  if(DODISSMEASURE==0) return(-1.0);

  if(NSPECIAL>0){
    FTYPE dUgeomtemp[NPR+NSPECIAL];
    FTYPE uitemp[NPR+NSPECIAL];
    FTYPE uftemp[NPR+NSPECIAL];
    FTYPE tempucumtemp[NPR+NSPECIAL];
    PLOOP(pliter,pl){
      dUgeomtemp[pl]=0.0; // geometry not relevant for this calculation
      uitemp[pl]=ui[pl];
      uftemp[pl]=uf[pl];
      tempucumtemp[pl]=tempucum[pl];
    }
    //
    int plsp;
    int specialfrom;
    int truenspecial=0;
#if(NSPECIAL==6)

#if(EOMRADTYPE!=EOMRADNONE)
    int plspeciallist[NSPECIAL]={SPECIALPL1,SPECIALPL2,SPECIALPL3,SPECIALPL4,SPECIALPL5,SPECIALPL6};
    truenspecial=NSPECIAL;
#else
    int plspeciallist[NSPECIAL]={SPECIALPL1,SPECIALPL2,SPECIALPL3,SPECIALPL4,SPECIALPL5};
    truenspecial=5;
#endif

#elif(NSPECIAL==2)
    int plspeciallist[NSPECIAL]={SPECIALPL1,SPECIALPL2};
    truenspecial=NSPECIAL;
#elif(NSPECIAL==1)
    int plspeciallist[NSPECIAL]={SPECIALPL1};
    truenspecial=NSPECIAL;
#endif
    PLOOPSPECIALONLY(plsp,truenspecial){
      specialfrom=plspeciallist[plsp-NPR];
      dUgeomtemp[plsp]=0.0;
      uitemp[plsp] = uitemp[specialfrom];
      uftemp[plsp] = uftemp[specialfrom];
      tempucumtemp[plsp] = tempucumtemp[specialfrom];
    }
    //
    FTYPE dUriemanntemp[NPR+NSPECIAL]={0},dUriemann1temp[NPR+NSPECIAL]={0},dUriemann2temp[NPR+NSPECIAL]={0},dUriemann3temp[NPR+NSPECIAL]={0};
    //
    // get flux update for all NPR and NSEPCIAL quantities
    flux2dUavg(DOSPECIALPL,i,j,k,F1,F2,F3,dUriemann1temp,dUriemann2temp,dUriemann3temp);
    //
    // get dUriemann for NPR and truenspecial quantities
    PLOOP(pliter,pl) dUriemanntemp[pl]=dUriemann1temp[pl]+dUriemann2temp[pl]+dUriemann3temp[pl];
    PLOOPSPECIALONLY(plsp,truenspecial){
      dUriemanntemp[plsp]=dUriemann1temp[plsp]+dUriemann2temp[plsp]+dUriemann3temp[plsp];
    }
    // Get final U for NPR and NSPECIAL quantities
    //    FTYPE dUcomptemp[NUMSOURCES][NPR+NSPECIAL];
    FTYPE dUcomptemp[NUMSOURCES][NPR];
    int sc;
    PLOOP(pliter,pl) SCLOOP(sc) dUcomptemp[sc][pl]=0.0;
    //    PLOOPSPECIALONLY(plsp,truenspecial){
    //      SCLOOP(sc) dUcomptemp[sc][plsp]=0.0;
    //    }
    dUtoU(timeorder,DOSPECIALPL,i,j,k,loc,dUgeomtemp, dUcomptemp, dUriemanntemp, CUf, CUnew, uitemp, uftemp, tempucumtemp);



    //////////////
    //
    // now get dissipation measure.  dissmeasure<0.0 means dissipation occuring (e.g. in density field)
    //
    /////////////
    
    // First, get dUnondiss and dUdiss and norm
    FTYPE dUdissplusnondiss[NPR];
    FTYPE dUnondiss[NPR];
    FTYPE dUdiss[NPR];
    FTYPE dissmeasurepl[NPR];
    FTYPE signdiss;
    FTYPE norm[NPR];
    PLOOPSPECIALONLY(plsp,truenspecial){
      specialfrom=plspeciallist[plsp-NPR];
      dUdissplusnondiss[specialfrom]=(uftemp[specialfrom] - uitemp[specialfrom]);
      dUnondiss[specialfrom]=(uftemp[plsp] - uitemp[plsp]);
      dUdiss[specialfrom]=dUdissplusnondiss[specialfrom] - dUnondiss[specialfrom];
            
      norm[specialfrom]=SMALL+(fabs(dUdiss[specialfrom])+fabs(dUnondiss[specialfrom]));
    }

    // get actual norm that's over multiple components for the magentic field
    FTYPE actualnorm[NPR];
    PLOOPSPECIALONLY(plsp,truenspecial){
      specialfrom=plspeciallist[plsp-NPR];
      if(specialfrom==RHO || specialfrom==UU){
        actualnorm[specialfrom] = SMALL + norm[specialfrom];
      }
      else if(specialfrom==B1 || specialfrom==B2 || specialfrom==B3){
        // for magnetic field, actualnorm is computed as:
        // |dUdiss(B1)|^2 |dUnondiss(B1)|^2 + |dUdiss(B2)|^2 + |dUnondiss(B2)|^2 + |dUdiss(B3)|^2 + |dUnondiss(B3)|^2 + |dUdiss(UU)| + |dUnondiss(UU)|
        // So we account for any magnetic energy change and help with normalization by also using total energy change
        actualnorm[specialfrom]=0.0;
        int jj;
        SLOOPA(jj){
          pl=B1+jj-1;
          actualnorm[specialfrom] += SMALL+(fabs(dUdiss[pl]*dUdiss[pl])+fabs(dUnondiss[pl]*dUnondiss[pl]));
        }
        actualnorm[specialfrom] += fabs(dUdiss[UU]) + fabs(dUnondiss[UU]); // add energy as reference for normalization to avoid weak magnetic fields suggesting strong dissipation.  GODMARK: But, kinda risky, because nominally want to use energy to capture any field dissipation, not just when the field is important to total energy dissipation.
        //  But much better than (say) using total energy density as reference.  We  use actual energy dissipation as reference.  So if reconnection is at all important to total dissipation, it will be accounted for.  When total energy dissipation is dominated by non-reconnection, that missed dissipation is a small correction.
        // But, because total energy dissipation is less trustable due to average vs. point issue in general creating artificial heating or cooling, this is a bit less trustable as even a normalization.
      }
      else if(specialfrom==URAD0){
        actualnorm[specialfrom] = SMALL + norm[specialfrom];
      }
    }
    
    /////////////
    //
    // compute dissmeasure for each quantity using actualnorm
    //
    /////////////
    PLOOPSPECIALONLY(plsp,truenspecial){
      specialfrom=plspeciallist[plsp-NPR];

      // set sign in front of dissmeasure
      if(specialfrom==RHO) signdiss=+1.0; // dUdiss>0 means adding density due to dissipation (i.e. convergence, not divergence)
      else if(specialfrom==UU) signdiss=-1.0; // dUdiss<0 means adding energy due to dissipation (i.e. dissipation)
      else if(specialfrom==B1 || specialfrom==B2 || specialfrom==B3) signdiss=-1.0; // dUdiss<0 means lost magentic energy due to dissipation (i.e. dissipation)
      else if(specialfrom==URAD0) signdiss=-1.0; // dUdiss<0 means adding energy due to dissipation (i.e. dissipation)

      // get dissmeasure, which is normalized dUdiss with correct sign
      if(specialfrom==RHO || specialfrom==UU || specialfrom==URAD0){
        // standard dUdiss/actualnorm
        dissmeasurepl[specialfrom] = -signdiss*(dUdiss[specialfrom])/actualnorm[specialfrom];
      }
      else if(specialfrom==B1 || specialfrom==B2 || specialfrom==B3){
        // special use of square of dUdiss: dUdiss^2(B1,B2,B3)/actualnorm
        dissmeasurepl[specialfrom] = -signdiss*(sign(dUdiss[specialfrom])*dUdiss[specialfrom]*dUdiss[specialfrom])/actualnorm[specialfrom];
      }

      //      dualfprintf(fail_file,"plsp=%d specialfrom=%d dissmeasurepl=%g\n",plsp,specialfrom,dissmeasurepl[specialfrom]);

      //    dualfprintf(fail_file,"DISS (ijk=%d %d %d): (%d %d) ui=%g %g : uf=%g %g : dUdiss=%g dUnondiss=%g : dU1=%g %g : dU2=%g %g : dissmeasure=%g\n",i,j,k,specialfrom,plsp,uitemp[specialfrom],uitemp[plsp],uftemp[specialfrom],uftemp[plsp],dUdiss,dUnondiss,dUriemann1temp[specialfrom]*CUf[2]*dt,dUriemann1temp[plsp]*CUf[2]*dt,dUriemann2temp[specialfrom]*CUf[2]*dt,dUriemann2temp[plsp]*CUf[2]*dt,dissmeasure);

      // DEBUG
      if(DODISSMEASURE){
        GLOBALMACP0A1(dissmeasurearray,i,j,k,plsp-NPR)=dissmeasurepl[specialfrom];
      }

    }

    ///////////
    //
    // Choose which *total* dissipation measurement to use.  Must somehow combine different original equations of motion into a single measure since energy vs. entropy changes involve energy equation that combines both shocks and reconnection.
    //
    ///////////
    if(truenspecial==1){
      dissmeasure=dissmeasurepl[SPECIALPL1];
    }
    // Can't trust energy, since average vs. point issue itself leads to erreoneous apparent dissipation (which when using energy equation alone is what leads to inappropriate heating in divergent flows, as well as inappropriate cooling in divergent flows, often in stripes in u_g due to sensitivity to higher-order interpolation schemes like PPM).
    //          dissmeasure=dissmeasurepl[SPECIALPL2];
    else if(truenspecial>=5){
      // can't trust energy, but density doesn't account for magnetic energy.  So use density for shocks and magnetic energy flux for reconnection.
      dissmeasure=dissmeasurepl[SPECIALPL1];

      FTYPE dissmeasurefield=dissmeasurepl[B1];
      dissmeasurefield=MIN(dissmeasurefield,dissmeasurepl[B2]);
      dissmeasurefield=MIN(dissmeasurefield,dissmeasurepl[B3]);
      
      // add field dissipation measure, but only account if above machine error (i.e. relative to total energy as defined above)
      dissmeasure=MIN(dissmeasure,10.0*NUMEPSILON+dissmeasurefield);

      if(truenspecial>=6 && SPECIALPL6>=0){
        // add radiation pressure to total pressure if optically thick
        FTYPE tautot[NDIM],tautotmax;
        calc_tautot(pr, ptrgeom, NULL, tautot, &tautotmax); // high accuracy not required

        FTYPE dissswitch=MIN(tautotmax/TAUTOTMAXSWITCH,1.0);

        dissmeasure = MIN(dissmeasure,dissmeasurepl[SPECIALPL6])*dissswitch + dissmeasure*(1.0-dissswitch);
      }
    }

    if(truenspecial>0 && DODISSMEASURE){
      // DEBUG
      GLOBALMACP0A1(dissmeasurearray,ptrgeom->i,ptrgeom->j,ptrgeom->k,NSPECIAL)=dissmeasure;
    }


  }

  dissmeasure=-1.0; // force use of energy unless energy fails.

  return(dissmeasure);
}

















/// this method guarantees conservation of non-sourced conserved quantities when metric is time-dependent
/// this method has updated field staggered method
/// Note that when dt==0.0, assume no fluxing, just take ucum -> ui -> {uf,ucum} and invert.  Used with metric update.
static int advance_standard_orig(
                            int truestep,
                            int stage,
                            FTYPE (*pi)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pb)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                            FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                            FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                            FTYPE (*ui)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*uf)[NSTORE2][NSTORE3][NPR],
                            FTYPE (*ucum)[NSTORE2][NSTORE3][NPR], 
                            FTYPE *CUf, 
                            FTYPE *CUnew, 
                            SFTYPE fluxdt,
                            SFTYPE boundtime,
                            SFTYPE fluxtime,
                            int timeorder, int numtimeorders,
                            FTYPE *ndt)
{
  FTYPE ndt1, ndt2, ndt3;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  static SFTYPE dt4diag_willbe=0;
  int finalstep,initialstep;
  FTYPE accdt, accdt_ij;
  int accdti,accdtj,accdtk;
  FTYPE gravitydt, gravitydt_ij;
  int gravitydti,gravitydtj,gravitydtk;
  //  FTYPE (*dUriemannarray)[NSTORE2][NSTORE3][NPR];
  FTYPE (*ucumformetric)[NSTORE2][NSTORE3][NPR];
  int enerregion;
  int *localenerpos;
  int jj;
  FTYPE (*utoinvert)[NSTORE2][NSTORE3][NPR];
  FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR];
  FTYPE (*useducum)[NSTORE2][NSTORE3][NPR];
  FTYPE (*myupoint)[NSTORE2][NSTORE3][NPR];
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int is,ie,js,je,ks,ke;
  int doingextrashiftforstag;






  if(FLUXB==FLUXCTSTAG){
    // fill in tempucum with changes that are not controlled well in space, but later fill in ucum from this in controlled way
    tempucum=GLOBALPOINT(utemparray);
    useducum=tempucum; // unless changed
  }
  else{
    // no special shifting of indices occurs that requires loop shifting
    tempucum=ucum;
    useducum=ucum;
  }


  ucumformetric=GLOBALPOINT(ucumformetric);// temporary space for ucum for metric that is same time as "pb", so not updated yet or is ui



  /////////////////////////////////////////////
  //
  // Setup function tasks
  //
  ////////////////////////////////////////////


  // for ZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];


  accdt=BIG; // initially no limit to dt due to acceleration
  accdti=accdtj=accdtk=-100;
  gravitydt=BIG; // initially no limit to dt due to time derivatives in gravity
  gravitydti=gravitydtj=gravitydtk=-100;




#if(ASYMDIAGCHECK)
  dualfprintf(fail_file,"BEGINNING steppart=%d nstep=%ld\n",steppart,nstep);

  dualfprintf(fail_file,"pi in advance\n");
  asym_compute_1(pi);

  dualfprintf(fail_file,"pb in advance\n");
  asym_compute_1(pb);
#endif

  /////////
  //
  // set initialstep and finalstep to tell some procedures and diagnostic functions if should be accounting or not
  //
  /////////
  if(timeorder==0){
    initialstep=1;
  }
  else{
    initialstep=0;
  }

  if(timeorder==numtimeorders-1){
    finalstep=1;
  }
  else{
    finalstep=0;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // whether need to compute flux
  // only compute flux if flux is added to get Uf or Ucum
  ////////////////////////////////////////////////////////////////////////////////
  int doflux = (CUf[2]!=0.0 || CUnew[1]!=0.0);
  // current implicit dt factor
  // Since we don't pass CUnew[NUMPREDTCUFS+timeorder], we assume never try to add M^i into U^{n+1} unless computed for current U^i or previous U^i
  FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];



  /////////
  //
  // set dt4diag for source diagnostics
  //
  /////////
  if(timeorder==numtimeorders-1 && (nstep%DIAGSOURCECOMPSTEP==0) ){
    // every 4 full steps as well as on final step of full step (otherwise diag_source_comp() too expensive)
    dt4diag=dt4diag_willbe;
    dt4diag_willbe=0;
  }
  else{
    dt4diag_willbe+=dt;
    dt4diag=-1.0;
  }



  /////////////////////////////////////////////
  //
  // Setup loops [+1 extra compared to normal comp region if doing FLUXCTSTAG]
  //
  ////////////////////////////////////////////
  get_flux_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);



  /////////////////////////////////////////////
  //
  // Set initial ui, temporary space, so ok that COMPZLOOP() goes over +1 in FLUXB==FLUXCTSTAG case
  //
  ////////////////////////////////////////////
  if(timeorder==0){
    // last timestep's final ucum is stored into ucumformetric and ui as initial term in ucum
    // copy ucum -> {ucumformetric,ui}
    if(doingmetricsubstep()) copy_3dnpr_2ptrs(is,ie,js,je,ks,ke,ucum,ucumformetric,ui);
    else copy_3dnpr(is,ie,js,je,ks,ke,ucum,ui); // only need to setup ui then
  }
  else{
    // preserve this time's value of ucum for the metric (for timeorder==0 ucumformetric is assigned from ui)
    // copy ucum -> ucumformetric
    if(doingmetricsubstep()) copy_3dnpr(is,ie,js,je,ks,ke,ucum,ucumformetric);
  }


  /////////////////////////////////////////////
  //
  // Compute flux
  //
  ////////////////////////////////////////////

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif




  ndt1=ndt2=ndt3=BIG;
  if(doflux && truestep){ // only do if not just passing through

    if(1){
      // NORMAL:
      // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
      MYFUN(fluxcalc(stage,initialstep,finalstep,pb,pstag,pl_ct, pr_ct, vpot,F1,F2,F3,CUf,CUnew,fluxdt,fluxtime,&ndt1,&ndt2,&ndt3),"advance.c:advance_standard_orig()", "fluxcalcall", 1);
    }

  
#if(0)// DEBUG:
    if(1){
      ndt1donor=ndt2donor=ndt3donor=BIG;
      MYFUN(fluxcalc_donor(stage,pb,pstag,pl_ct, pr_ct, vpot,GLOBALPOINT(F1EM),GLOBALPOINT(F2EM),GLOBALPOINT(F3EM),CUf,CUnew,fluxdt,fluxtime,&ndt1donor,&ndt2donor,&ndt3donor),"advance.c:advance_standard_orig()", "fluxcalcall", 1);
    }
    // DEBUG:
    if(1){
      int i,j,k,pl,pliter;
      FULLLOOP{
        PLOOP(pliter,pl){
          if(N1>1) MACP0A1(F1,i,j,k,pl)=GLOBALMACP0A1(F1EM,i,j,k,pl);
          if(N2>1) MACP0A1(F2,i,j,k,pl)=GLOBALMACP0A1(F2EM,i,j,k,pl);
          if(N3>1) MACP0A1(F3,i,j,k,pl)=GLOBALMACP0A1(F3EM,i,j,k,pl);
        }
      }
      ndt1=ndt1donor;
      ndt2=ndt2donor;
      ndt3=ndt3donor;
    }
#endif


#if(0)// DEBUG:
    if(1){
      int i,j,k,pliter,pl;
      dualfprintf(fail_file,                 "BADLOOPCOMPARE: nstep=%ld steppart=%d\n",nstep,steppart);
      dualfprintf(fail_file,                 "ndt1orig=%21.15g ndt1new=%21.15g ndt2orig=%21.15g ndt2new=%21.15g\n",ndt1,ndt1donor,ndt2,ndt2donor);
      COMPFULLLOOP{
        if(i==390 && j==1 && k==0){
          dualfprintf(fail_file,                 "i=%d j=%d k=%d\n",i,j,k);
          PLOOP(pliter,pl) dualfprintf(fail_file,"          pl=%d F1orig=%21.15g F1new=%21.15g  :: F2orig=%21.15g F2new=%21.15g \n",pl,MACP0A1(F1,i,j,k,pl),GLOBALMACP0A1(F1EM,i,j,k,pl),MACP0A1(F2,i,j,k,pl),GLOBALMACP0A1(F2EM,i,j,k,pl));
        }
      }
    }
#endif
  
#if(0)// DEBUG:
    if(1){
      int i,j,k,pliter,pl;
      FTYPE diff1,diff2;
      FTYPE sum1,sum2;
      dualfprintf(fail_file,                 "BADLOOPCOMPARE: nstep=%ld steppart=%d\n",nstep,steppart);
      dualfprintf(fail_file,                 "ndt1orig=%21.15g ndt1new=%21.15g ndt2orig=%21.15g ndt2new=%21.15g\n",ndt1,ndt1donor,ndt2,ndt2donor);
      COMPFULLLOOP{
        PLOOP(pliter,pl){
          diff1=fabs(MACP0A1(F1,i,j,k,pl)-GLOBALMACP0A1(F1EM,i,j,k,pl));
          sum1=fabs(MACP0A1(F1,i,j,k,pl))+fabs(GLOBALMACP0A1(F1EM,i,j,k,pl))+SMALL;
          diff2=fabs(MACP0A1(F2,i,j,k,pl)-GLOBALMACP0A1(F2EM,i,j,k,pl));
          sum2=fabs(MACP0A1(F2,i,j,k,pl))+fabs(GLOBALMACP0A1(F2EM,i,j,k,pl))+SMALL;
          if(diff1/sum1>100.0*NUMEPSILON || diff2/sum2>100.0*NUMEPSILON){
            dualfprintf(fail_file,                 "i=%d j=%d k=%d\n",i,j,k);
            dualfprintf(fail_file,"          pl=%d diff1/sum1=%21.15g F1orig=%21.15g F1new=%21.15g  :: diff2/sum2=%21.15g F2orig=%21.15g F2new=%21.15g \n",pl,diff1/sum1,MACP0A1(F1,i,j,k,pl),GLOBALMACP0A1(F1EM,i,j,k,pl),diff2/sum2,MACP0A1(F2,i,j,k,pl),GLOBALMACP0A1(F2EM,i,j,k,pl));
          }
        }
      }
    }
#endif

  }// end if not just passing through






#if(PRODUCTION==0)
  trifprintf( "1f");
#endif




  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil


  




  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //
  // now update get flux update [loop should go over normal computational region +1 on outer edges so automatically set field if staggered.  Since we only set tempucum so far, ucum in that +1 region is unaffected]
  //
  /////////////////////////////////////////////////////


  int didreturnpf;
  if(truestep){


    // initialize uf and ucum if very first time here since ucum is cumulative (+=) [now tempucum is cumulative]
    // copy 0 -> {uf,tempucum}
    if(timeorder==0) init_3dnpr_2ptrs(is, ie, js, je, ks, ke,0.0, uf,tempucum);



#if(WHICHEOM==WITHNOGDET && (NOGDETB1==1 || NOGDETB2==1 || NOGDETB3==1) )
    // for FLUXB==FLUXCTSTAG, assume no source term for field
    if(FLUXB==FLUXCTSTAG){
      dualfprintf(fail_file,"Not setup for field source term if staggered field\n");
      myexit(176023);
    }
#endif


#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM
    {
      int pl,pliter,i,j,k;
      struct of_geom geomdontuse;
      struct of_geom *ptrgeom=&geomdontuse;
      FTYPE Uitemp[NPR];
      FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR+NSPECIAL],dUriemann2[NPR+NSPECIAL],dUriemann3[NPR+NSPECIAL],dUcomp[NUMSOURCES][NPR];
      struct of_state qdontuse;
      struct of_state *qptr=&qdontuse;
      struct of_state qdontuse2;
      struct of_state *qptr2=&qdontuse2; // different qptr since call normal and special get_state()
      // setup default eomtype
      int eomtype;


      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

      
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // setup default eomtype
        eomtype=EOMDEFAULT;


        // set geometry for centered zone to be updated
        get_geometry(i, j, k, CENT, ptrgeom);
      
        // find Ui(pi)
        // force use of primitive to set Ui since otherwise wherever corrected/changed primitive (in fixup, etc.) then would have to change conserved quantity, while same since both are point values
        // only field for staggered method is special point value at faces that needs to come from conserved quantity
        MYFUN(get_state(MAC(pi,i,j,k), ptrgeom, qptr),"step_ch.c:advance()", "get_state()", 1);
        MYFUN(primtoU(UEVOLVE,MAC(pi,i,j,k), qptr, ptrgeom, Uitemp, NULL),"step_ch.c:advance()", "primtoU()", 1);

        if(FLUXB==FLUXCTSTAG || DOENOFLUX != NOENOFLUX ){
          // then field version of ui[] is stored as "conserved" value at FACE, not CENT
          PLOOPNOB1(pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // CENT
          //PLOOPBONLY(pl) MACP0A1(ui,i,j,k,pl) is itself // FACE (see step_ch.c's setup_rktimestep and know that ui=unew for before first substep)
          PLOOPNOB2(pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // CENT
        }
        else{
          PLOOP(pliter,pl) MACP0A1(ui,i,j,k,pl)=Uitemp[pl]; // all at CENT
        }
 
        // dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
        if(doflux){
          flux2dUavg(DOALLPL,i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
          PLOOP(pliter,pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is one type of avg->point mistake
        }
        else{
          PLOOP(pliter,pl) dUriemann[pl]=0.0;
        }


        /////////////
        //
        // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
        //
        // MAJORNOTE: Since final step often sets pointers of pi=pf, in order to use arbitrary guess, must set here once done with pi,pb,pf.
        //
        ////////////
        // copy pb->pf as initial pf
        // can do this now that don't need pi anymore
        PALLLOOP(pl) MACP0A1(pf,i,j,k,pl) = MACP0A1(pb,i,j,k,pl);


        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        //
        // now update get source update (only affects stress-energy tensor in general)
        // [Loop goes over up to +1 normal computational region for getting new staggered U if doing FLUXCTSTAG]
        //
        /////////////////////////////////////////////////////
 
        // get state since both source() and dUtodt() need same state
        // From pb, so different than state for Ui(pi)
        MYFUN(get_stateforsource(MAC(pb,i,j,k), ptrgeom, &qptr2) ,"advance.c:()", "get_state() dir=0", 1);
      

        // note that uf and ucum are initialized inside setup_rktimestep() before first substep


        // find dU(pb)
        MYFUN(source(MAC(pi,i,j,k), MAC(pb,i,j,k), MAC(pf,i,j,k), &didreturnpf, &eomtype, ptrgeom, qptr2, MAC(ui,i,j,k), MAC(uf,i,j,k), CUf, CUimp, 0.0, dUriemann, dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
        // assumes final dUcomp is nonzero and representative of source term over this timestep
        // KORALTODO: Not using eomtype for this method because outside loop.  Would have to store eomtype in array or something!
 
        // DEBUG (to check if source updated pf guess)
        //        if(didreturnpf==1) dualfprintf(fail_file,"didreturnpf=%d\n",didreturnpf);


#if(SPLITNPR)
        // don't update metric if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
            if(DODIAGS){
#if(ACCURATESOURCEDIAG>=1)
              diag_source_all(ptrgeom,dUgeom,fluxdt);
#else
              diag_source_all(ptrgeom,dUgeom,dt4diag);
#endif
#if(ACCURATESOURCEDIAG>=2)
              diag_source_comp(ptrgeom,dUcomp,fluxdt);
#else
              diag_source_comp(ptrgeom,dUcomp,dt4diag);
#endif
            }
          }

      

        // Get update
        dUtoU(timeorder,DOALLPL, i,j,k,ptrgeom->p,dUgeom, dUcomp, dUriemann, CUf, CUnew, MAC(ui,i,j,k), MAC(uf,i,j,k), MAC(tempucum,i,j,k));
      
      
      
        // get timestep limit from acceleration
#if(LIMITDTWITHSOURCETERM)
#if(SPLITNPR)
        // don't do dUtodt if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
            // geometry is post-metric update, but should still give good estimate of future dt
            dUtodt(ptrgeom, qptr2, MAC(pb,i,j,k), dUgeom, dUriemann, dUcomp[GEOMSOURCE], &accdt_ij, &gravitydt_ij);

#pragma omp critical
            {
              if(accdt_ij<accdt){
                accdt=accdt_ij;
                accdti=i;
                accdtj=j;
                accdtk=k;
              }
              if(gravitydt_ij<gravitydt){
                gravitydt=gravitydt_ij;
                gravitydti=i;
                gravitydtj=j;
                gravitydtk=k;
              }
            }// end critical region
          }// end block that may mean: if not only doing B1-B3

#endif // end if doing LIMITDTWITHSOURCETERM
          





#if(FLUXDUMP==1)
        fluxdumpdt=dt; // store current dt so when dump fluxdump use that dt instead of updated dt
        fluxdumprealnstep=realnstep;
        // DEBUG - DIAG:
        PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,0*NPR + pl)=dUgeom[pl];

        if(N1>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,1*NPR + pl)=dUriemann1[pl];
        else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,1*NPR + pl)=0.0;
  
        if(N2>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,2*NPR + pl)=dUriemann2[pl];
        else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,2*NPR + pl)=0.0;

        if(N3>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,3*NPR + pl)=dUriemann3[pl];
        else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,3*NPR + pl)=0.0;
#endif

      } // end COMPZLOOP :: end looping to obtain dUriemann and full unew update
    }// end parallel block
  } // end if truestep
  else{
    // then nothing to do since nothing changed
    // previously updated dU and got new ucum as fed into metric, but now metric has its own ucummetric so all that is not needed
    // SUPERGODMARK: no, I guess I don't recall what was being done for metric and why when passing through with dt==0.0

    // just need to copy ui -> {uf,tempucum} to duplicate behavior of dUtoU()
    copy_3dnpr_2ptrs(is,ie,js,je,ks,ke,ui,uf,tempucum);
  }







#if(PRODUCTION==0)
  trifprintf( "#0m");
#endif

  
    
  /////////////////////////////////////////////
  //
  // EVOLVE (update/compute) METRIC HERE
  // In general computes stress-energy tensor (T) from pb and T is then used to compute metric
  // Done here instead of after flux since horizon_flux() updates flux through boundary that would change metric
  // want metric to be at same time as everything else done with pb so RK method is consistent
  //
  // uses unew that's NOT updated yet
  /////////////////////////////////////////////
  if(truestep){
    if(doingmetricsubstep()){
#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
        {
          compute_new_metric_substep(CUf,CUnew,pb,ucumformetric); // CHANGINGMARK : Not sure if CUnew here is correct
        }
    }// end if doing metric substepping
  }




  ////////////////////////////////
  //
  // compute flux diagnostics (accurately using all substeps)
  //
  ///////////////////////////////
  if(doflux && truestep){
    // must come after metric changes that can change where flux surface is since otherwise when flux surface changes, we won't track this substep's flux through the new surface but the old surface (which could even be at r=0 and have no flux)
    // if using unew, then since metric update above uses old unew, need present flux at new horizon surface
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
        diag_flux_pureflux(pb,F1, F2, F3, fluxdt); // fluxdt is true dt for this flux as added in dUtoU() as part of ucum
#endif
      }
  }






  
#if(PRODUCTION==0)
  trifprintf( "#0s");
#endif





  ///////////////////////////////////////
  //
  // Copy over tempucum -> ucum per pl to account for staggered field
  //
  // And choose which RK-quantity to invert
  //
  ///////////////////////////////////////
  if(finalstep){
    if(FLUXB==FLUXCTSTAG){
      // copy over new ucum in only desired locations irrespective of where tempucum was updated
      // copy tempucum -> ucum
      copy_tempucum_finalucum(DOALLPL,Uconsevolveloop,tempucum,ucum);
    }
    utoinvert = ucum;
    useducum=ucum;
  }
  else{
    // tempucum cumulates for now
    utoinvert = uf;
    // tempucum cumulates for now
    useducum=tempucum;
  }





  ////////////////////////////
  //
  // split ZLOOP above and below to allow staggered field method
  // [loop only over "centered" cells]
  //
  ////////////////////////////
  if(FLUXB==FLUXCTSTAG){
    // if using staggered grid for magnetic field, then need to convert ucum to pstag to new pb/pf


    // GODMARK: Use of globals
    myupoint=GLOBALPOINT(upointglobal);

    // first copy over all quantities as point, which is true except for fields if FLUXRECON active
    // copy utoinvert -> myupoint
    // copy all pl's since myupoint eventually used to invert rest of non-field quantities
    copy_tempucum_finalucum(DOALLPL,Uconsevolveloop,utoinvert,myupoint);


    if(extrazones4emf && dofluxreconevolvepointfield==0){
      //bound_uavg(STAGEM1,utoinvert); // DEBUG
      // uses utoinvert and gets reaveraged field into myupoint
      field_integrate_fluxrecon(stage, pb, utoinvert, myupoint);
    }


    // first pb entry is used for shock indicator, second is filled with new field
    // myupoint goes in as staggered point value of magnetic flux and returned as centered point value of magnetic flux
    // first pb entry is used for shock indicator, second is filled with new field
    // myupoint goes in as staggered point value of magnetic flux and returned as centered point value of magnetic flux
    if(1){ // must manually override to use below donor version
      interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstag, myupoint, pf);
    }
    else{
      //      interpolate_ustag2fieldcent_donor(stage, boundtime, timeorder, numtimeorders, pb, pstag, myupoint, pf);
    }

    ////////////////////    
    // now myupoint contains centered point conserved quantities ready for inversion
    ////////////////////

  }
  else{
    // utoinvert never reassigned from global a_utoinvert assignment since if here not doing FLUXCTSTAG
    myupoint=utoinvert;
  }








  ////////////////////////////
  //
  // INVERT [loop only over "centered" cells]
  //
  ////////////////////////////


  // get loop range
  get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);

#pragma omp parallel OPENMPGLOBALPRIVATEFORINVERSION
  {
    int pl,pliter,i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE prbefore[NPR];
    struct of_newtonstats newtonstats;  setnewtonstatsdefault(&newtonstats);
    int showmessages=1;
    int allowlocalfailurefixandnoreport=1; // allow local fixups

    // setup default eomtype
    int eomtype;
    
    OPENMP3DLOOPVARSDEFINE;  OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;


#pragma omp for schedule(OPENMPVARYENDTIMESCHEDULE(),OPENMPCHUNKSIZE(blocksize)) reduction(+: nstroke)
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
 
      // setup default eomtype
      // KORALTODO: Note that eomtype should have been set from source(), but different loop.  But not using this "orig"  method for koral evolution, so ok.
      eomtype=EOMDEFAULT;

      // set geometry for centered zone to be updated
      get_geometry(i, j, k, CENT, ptrgeom);

      
      // invert U->p
      if(finalstep){ // last call, so ucum is cooked and ready to eat!
        // store guess for diss_compute before changed by normal inversion
        PALLLOOP(pl) prbefore[pl]=MACP0A1(pf,i,j,k,pl);

        FTYPE dissmeasure=-1.0; // assume energy try ok
        int whichcap=CAPTYPEBASIC;
        int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
        int modprim=0;
        int checkoninversiongas=CHECKONINVERSION;
        int checkoninversionrad=CHECKONINVERSIONRAD;
        MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE,MAC(myupoint,i,j,k), NULL, ptrgeom, dissmeasure, MAC(pi,i,j,k), MAC(pf,i,j,k),&newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
        nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;


#if(DODISS||DODISSVSR)
        // then see what entropy inversion would give
        diss_compute(EVOLVEUTOPRIM,UEVOLVE,MAC(myupoint,i,j,k),ptrgeom,prbefore,MAC(pf,i,j,k),&newtonstats);
#endif
 
      }
      else{ // otherwise still iterating on primitives
        FTYPE dissmeasure=-1.0; // assume energy try ok
        int whichcap=CAPTYPEBASIC;
        int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
        int modprim=0;
        int checkoninversiongas=CHECKONINVERSION;
        int checkoninversionrad=CHECKONINVERSIONRAD;
        MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE,MAC(myupoint,i,j,k), NULL, ptrgeom, dissmeasure, MAC(pi,i,j,k), MAC(pf,i,j,k),&newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
        nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;
      }




      
#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
        {
          // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
          // SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
          if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(finalstep)){
            int docorrectucons=0;
            MYFUN(fixup1zone(docorrectucons,MAC(pf,i,j,k),MAC(useducum,i,j,k), ptrgeom,finalstep),"fixup.c:fixup()", "fixup1zone()", 1);
          }
#endif
        }


    }// end COMPZLOOP
  }// end parallel section









#if(ASYMDIAGCHECK)
  dualfprintf(fail_file,"useducum in advance\n");
  asym_compute_2(useducum);

  dualfprintf(fail_file,"ENDING steppart=%d nstep=%ld\n",steppart,nstep);
#endif




  /////////////////////////////////
  //
  // If not fixing up primitives after inversion immediately, then fix up all zones at once afterwards
  //
  /////////////////////////////////

#if(SPLITNPR)
  // don't update metric if only doing B1-B3
  if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
#if(FIXUPZONES==FIXUPALLZONES)
      fixup(stage,pf,useducum,finalstep);
#endif  
    }


  /////////////////////////////////
  //
  // Determine next timestep from waves, fluxes, and source updates
  //
  /////////////////////////////////


  prepare_globaldt(truestep,ndt1,ndt2,ndt3,accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,ndt);



#if(PRODUCTION==0)
  trifprintf( "2f");
#endif

  return (0);
}

















/// finite volume method NOT SETUP FOR CONSISTENT METRIC EVOLUTION YET -- EASY, JUST NOT DOING IT YET -- FOLLOW ABOVE AS EXAMPLE OF WHAT TO DO
/// also not setup for staggered field method
static int advance_finitevolume(
                                int truestep,
                                int stage,
                                FTYPE (*pi)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pb)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                                FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                                FTYPE (*ui)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*uf)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*ucum)[NSTORE2][NSTORE3][NPR], 
                                FTYPE *CUf, 
                                FTYPE *CUnew, 
                                SFTYPE fluxdt,
                                SFTYPE boundtime,
                                SFTYPE fluxtime,
                                int timeorder, int numtimeorders,
                                FTYPE *ndt)
{
  int sc;
  FTYPE ndt1, ndt2, ndt3;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  static SFTYPE dt4diag_willbe=0;
  int finalstep,initialstep;
  FTYPE accdt, accdt_ij;
  int accdti,accdtj,accdtk;
  FTYPE gravitydt, gravitydt_ij;
  int gravitydti,gravitydtj,gravitydtk;
  int enerregion;
  int *localenerpos;
  int jj;
  FTYPE (*utoinvert)[NSTORE2][NSTORE3][NPR];
  FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR];
  FTYPE (*useducum)[NSTORE2][NSTORE3][NPR];
  FTYPE (*myupoint)[NSTORE2][NSTORE3][NPR];
  FTYPE (*ucumformetric)[NSTORE2][NSTORE3][NPR];
  FTYPE (*dUgeomarray)[NSTORE2][NSTORE3][NPR];
  // staggered field function:
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int docons,dosource;
  int locpl[NPR];
  int is,ie,js,je,ks,ke;
  int doingextrashiftforstag;


  if(FLUXB==FLUXCTSTAG){
    tempucum=GLOBALPOINT(utemparray);
    useducum=tempucum; // unless changed
  }
  else{
    // no special shifting of indices occurs that requires loop shifting
    tempucum=ucum;
    useducum=ucum;
  }


  ucumformetric=GLOBALPOINT(ucumformetric);// temporary space for ucum for metric that is same time as "pb", so not updated yet or is ui
  dUgeomarray=GLOBALPOINT(dUgeomarray);// temporary space for dUgeomarray


  
  {// begin block
    int pl,pliter;
    // any cons
    docons=0;
    PLOOP(pliter,pl) docons+=do_conserved_integration[pl];
    docons*=(DOENOFLUX == ENOFINITEVOLUME);

    // any source
    dosource=0;
    PLOOP(pliter,pl) dosource+=do_source_integration[pl];
    dosource*=(DOENOFLUX == ENOFINITEVOLUME);
  }// end block


  /////////////////////////////////////////////
  //
  // Setup loops and dt's
  //
  ////////////////////////////////////////////


  // for CZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];


  accdt=BIG; // initially no limit to dt due to acceleration
  accdti=accdtj=accdtk=-100;
  gravitydt=BIG; // initially no limit to dt due to time-changes in gravity
  gravitydti=gravitydtj=gravitydtk=-100;
 

  /////////
  //
  // set initialstep and finalstep to tell some procedures and diagnostic functions if should be accounting or not
  //
  /////////
  if(timeorder==0){
    initialstep=1;
  }
  else{
    initialstep=0;
  }

  if(timeorder==numtimeorders-1){
    finalstep=1;
  }
  else{
    finalstep=0;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // whether need to compute flux
  // only compute flux if flux is added to get Uf or Ucum
  ////////////////////////////////////////////////////////////////////////////////
  int doflux = (CUf[2]!=0.0 || CUnew[1]!=0.0);
  // current implicit dt factor
  // Since we don't pass CUnew[NUMPREDTCUFS+timeorder], we assume never try to add M^i into U^{n+1} unless computed for current U^i or previous U^i
  FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];

  /////////
  //
  // set dt4diag for source diagnostics
  //
  /////////
  if(timeorder==numtimeorders-1 && (nstep%DIAGSOURCECOMPSTEP==0) ){
    // every 4 full steps as well as on final step of full step (otherwise diag_source_comp() too expensive)
    dt4diag=dt4diag_willbe;
    dt4diag_willbe=0;
  }
  else{
    dt4diag_willbe+=dt;
    dt4diag=-1.0;
  }


  /////////////////////////////////////////////
  //
  // Setup loops
  //
  ////////////////////////////////////////////
  get_flux_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);



  /////////////////////////////////////////////
  //
  // Setup RK stuff
  //
  ////////////////////////////////////////////

  // setup RK's uinitial (needed since sometimes set ui=0 inside advance())
  // unew should be read in, now assign to uinitial for finite volume method or normal method when dt=0.0 or moving grid
  // below can't be CZLOOP since need uinitial in ghost+active region
  // GODMARK: since ui=ucum (average quantities) then if change primitive somehow (fixup, etc.) then must change corresponding average conserved quantity somehow (this is one reason why using average values is annoying, although in some cases we use Sasha's method to treat average like point for *change* in average conserved quantity)
  // otherwise assume ui didn't change.  Present RK schemes assume this.  Otherwise have to keep track of pf/Uf pairs in RK stepping
  /////////////////////////////////////////////
  //
  // Set initial ui, temporary space, so ok that COMPZLOOP() goes over +1 in FLUXB==FLUXCTSTAG case
  //
  ////////////////////////////////////////////
  if(timeorder==0){
    // last timestep's final ucum is stored into ucumformetric and ui as initial term in ucum
    // copy ucum -> {ucumformetric,ui}
    if(doingmetricsubstep()) copy_3dnpr_2ptrs(is,ie,js,je,ks,ke,ucum,ucumformetric,ui);
    else copy_3dnpr(is,ie,js,je,ks,ke,ucum,ui); // only need to setup ui then
  }
  else{
    // preserve this time's value of ucum for the metric (for timeorder==0 ucumformetric is assigned from ui)
    // copy ucum -> ucumformetric
    if(doingmetricsubstep()) copy_3dnpr(is,ie,js,je,ks,ke,ucum,ucumformetric);
  }



  /////////////////////////////////////////////
  //
  // Compute flux
  //
  ////////////////////////////////////////////

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif

  ndt1=ndt2=ndt3=BIG;
  if(doflux && truestep){
    // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
    MYFUN(fluxcalc(stage,initialstep,finalstep,pb,pstag,pl_ct, pr_ct, vpot,F1,F2,F3,CUf,CUnew,fluxdt,fluxtime,&ndt1,&ndt2,&ndt3),"advance.c:advance_standard_orig()", "fluxcalcall", 1);
  }

#if(PRODUCTION==0)
  trifprintf( "1f");
#endif
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil







#if(PRODUCTION==0)
  trifprintf( "#0m");
#endif


  /////////////
  //
  // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
  //
  // Since final step often sets pointers of pi=pf, in order to use arbitrary guess, must set here once done with pi,pb,pf.
  //
  ////////////
  // setup initial guess for inversion
  // use pb since should be closer to pf
  // copy pb->pf
  // OK to do here, because never really use pi
  copy_3dnpr_fullloop(pb,pf);


 
  /////////////////////////
  //
  // SOURCE TERM
  //
  ////////////////////////
  int didreturnpf; // used to have source() pass back better guess for Utoprimgen() than pb.

  if(truestep){
    // GODMARK: other/more special cases?
#if(WHICHEOM==WITHNOGDET && (NOGDETB1==1 || NOGDETB2==1 || NOGDETB3==1) )
    // for FLUXB==FLUXCTSTAG, assume no source term for field
    if(FLUXB==FLUXCTSTAG){
      dualfprintf(fail_file,"Not setup for field source term if staggered field\n");
      myexit(176023);
    }
#endif

    //    COMPZSLOOP(is,ie,js,je,ks,ke){
#if(STOREFLUXSTATE)
    // call get_stateforsource, which just reads existing data instead of generating from EOS, etc.
#pragma omp parallel 
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM
#endif
    {
      int i,j,k,pl,pliter;
      FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR+NSPECIAL],dUriemann2[NPR+NSPECIAL],dUriemann3[NPR+NSPECIAL],dUcomp[NUMSOURCES][NPR];
      struct of_geom geomdontuse;
      struct of_geom *ptrgeom=&geomdontuse;
      struct of_state qdontuse;
      struct of_state *qptr=&qdontuse;
      // setup default eomtype
      int eomtype;

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

      
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // setup default eomtype
        eomtype=EOMDEFAULT;

        // find dU(pb)
        // only find source term if non-Minkowski and non-Cartesian
        // set geometry for centered zone to be updated
        get_geometry(i, j, k, CENT, ptrgeom);

        // get state
        MYFUN(get_stateforsource(MAC(pb,i,j,k), ptrgeom, &qptr) ,"advance.c:()", "get_state() dir=0", 1);


        // dUriemann is volume averaged quantity (here this calcuation is done in case want to limit sources)
        if(doflux){
          flux2dUavg(DOALLPL,i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
          PLOOP(pliter,pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is entirely consistent with point->averages
        }
        else{
          PLOOP(pliter,pl) dUriemann[pl]=0.0;
        }

     
        // get source term (point source, don't use to update diagnostics)
        MYFUN(source(MAC(pi,i,j,k), MAC(pb,i,j,k), MAC(pf,i,j,k), &didreturnpf, &eomtype, ptrgeom, qptr, MAC(ui,i,j,k), MAC(uf,i,j,k), CUf, CUimp, 0.0, dUriemann, dUcomp, MAC(dUgeomarray,i,j,k)),"step_ch.c:advance()", "source", 1);
      }// end COMPZLOOP


      // Store diagnostics related to component form of sources. Done here since don't integrate-up compnents of source.  So only point accurate
#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
        {
#pragma omp critical // since diagnostics store in same global cumulative variables
          {
            if(DODIAGS){
#if(ACCURATESOURCEDIAG>=2)
              diag_source_comp(ptrgeom,dUcomp,fluxdt);
#else
              diag_source_comp(ptrgeom,dUcomp,dt4diag);
#endif
            }
          }
        }
      

      // volume integrate dUgeom
      // c2a_1 c2a_2 c2a_3
      if(dosource){
        // need to limit source averaging -- GODMARK
        PALLLOOP(pl) locpl[pl]=CENT;
        PALLLOOP(pl) whichpltoavg[pl]=do_source_integration[pl];// default
        PALLLOOP(pl) ifnotavgthencopy[pl]=1-do_source_integration[pl];// default
        avg2cen_interp(locpl,whichpltoavg, ifnotavgthencopy, ENOSOURCETERM, ENOCENT2AVGTYPE, pb, POINT(dUgeomarray), POINT(dUgeomarray));
      }
    }// end parallel region
  }// end if truestep










  ///////////////////////////////////////////////////////
  //
  // update Ui to Uf (ultimately to ucum)
  //
  ///////////////////////////////////////////////////////

  //  COMPZSLOOP(is,ie,js,je,ks,ke){
#if(STOREFLUXSTATE)
  // call get_stateforsource, which just reads existing data instead of generating from EOS, etc.
#pragma omp parallel 
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM
#endif
  {
    int i,j,k,pl,pliter;
    FTYPE fdummy;
    FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR+NSPECIAL],dUriemann2[NPR+NSPECIAL],dUriemann3[NPR+NSPECIAL],dUcomp[NUMSOURCES][NPR];
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    struct of_state qdontuse;
    struct of_state *qptr=&qdontuse;
    // GODMARK: KORALTODO: setup default eomtype (not using source() from above!)
    int eomtype;

    OPENMP3DLOOPVARSDEFINE;  OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

      
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      // setup default eomtype
      eomtype=EOMDEFAULT;


      // get geometry at center where source is located
      get_geometry(i, j, k, CENT, ptrgeom);

      // get state
      MYFUN(get_stateforsource(MAC(pb,i,j,k), ptrgeom, &qptr) ,"advance.c:()", "get_state() dir=0", 1);
  
      /////////////
      //
      // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
      //
      ////////////
  
      // initialize uf and ucum if very first time here since ucum is cumulative (+=)
      //    if(timeorder==0) PLOOP(pliter,pl) MACP0A1(uf,i,j,k,pl)=MACP0A1(tempucum,i,j,k,pl)=MACP0A1(ucum,i,j,k,pl)=0.0;
      if(timeorder==0) PLOOP(pliter,pl) MACP0A1(uf,i,j,k,pl)=MACP0A1(tempucum,i,j,k,pl)=0.0;
  
      // NEED TO DEFINE Ui on other substeps besides the 0th one
      // find Ui(pi)



      if(truestep){
        // get source term (volume integrated)
        PLOOP(pliter,pl) dUgeom[pl]=MACP0A1(dUgeomarray,i,j,k,pl);


        // store diagnostics related to source.  Totals, so use volume-integrated source.
#if(SPLITNPR)
        // don't update metric if only doing B1-B3
        if(advancepassnumber==-1 || advancepassnumber==1)
#endif
          {
#pragma omp critical // since diagnostics store in same global cumulative variables
            {
            if(DODIAGS){
#if(ACCURATESOURCEDIAG>=1)
              diag_source_all(ptrgeom,dUgeom,fluxdt);
#else
              diag_source_all(ptrgeom,dUgeom,dt4diag);
#endif
            }
          }   
        }
    

        // dUriemann is volume averaged quantity
        if(doflux){
          flux2dUavg(DOALLPL,i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
          PLOOP(pliter,pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is entirely consistent with point->averages
        }
        else{
          PLOOP(pliter,pl) dUriemann[pl]=0.0;
        }
      }
      else{
        PLOOP(pliter,pl){
          dUriemann[pl]=0.0;
          dUgeom[pl]=0.0;
        }
      }

      // find uf==Uf and additional terms to ucum
      dUtoU(timeorder,DOALLPL, i,j,k,ptrgeom->p,dUgeom, dUcomp, dUriemann, CUf, CUnew, MAC(ui,i,j,k), MAC(uf,i,j,k), MAC(tempucum,i,j,k));



#if(LIMITDTWITHSOURCETERM)
      // SUPERGODMARK: no longer have access to dUcomp : NEED TO FIX
      // below is correct, but excessive
      // get source term again in order to have dUcomp (NEED TO FIX)
      MYFUN(source(MAC(pi,i,j,k), MAC(pb,i,j,k), MAC(pf,i,j,k), &didreturnpf, &eomtype, ptrgeom, qptr, MAC(ui,i,j,k), MAC(uf,i,j,k), CUf, CUimp, 0.0, dUriemann, dUcomp, &fdummy),"step_ch.c:advance()", "source", 2);


      dUtodt(ptrgeom, qptr, MAC(pb,i,j,k), dUgeom, dUriemann, dUcomp[GEOMSOURCE], &accdt_ij, &gravitydt_ij);
      // below is old before grid sectioning
      //#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON>=1) )
      //    // avoid limiting dt if inside horizon
      //    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
      //#endif
#pragma omp critical
      {
        if(accdt_ij<accdt){
          accdt=accdt_ij;
          accdti=i;
          accdtj=j;
          accdtk=k;
        }
        if(gravitydt_ij<gravitydt){
          gravitydt=gravitydt_ij;
          gravitydti=i;
          gravitydtj=j;
          gravitydtk=k;
        }
      }// end critical region
#endif




#if(FLUXDUMP==1)
      PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,0*NPR + pl)=dUgeom[pl];

      if(N1>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,1*NPR + pl)=dUriemann1[pl];
      else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,1*NPR + pl)=0.0;
  
      if(N2>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,2*NPR + pl)=dUriemann2[pl];
      else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,2*NPR + pl)=0.0;

      if(N3>1) PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,3*NPR + pl)=dUriemann3[pl];
      else PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,3*NPR + pl)=0.0;
#endif
  

    }// end 3D LOOP
  }//end parallel region





  /////////////////////////////////////////////
  //
  // EVOLVE (update/compute) METRIC HERE
  // In general computes stress-energy tensor (T) from pb and T is then used to compute metric
  // Done here instead of after flux since horizon_flux() updates flux through boundary that would change metric
  // want metric to be at same time as everythin else done with pb so RK method is consistent
  //
  // uses unew that's NOT updated yet
  /////////////////////////////////////////////
  if(truestep){
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
        compute_new_metric_substep(CUf,CUnew,pb,ucumformetric); // CHANGINGMARK : Not sure if CUnew here is correct
      }
  }



  ////////////////////////////////
  //
  // compute flux diagnostics (accurately using all substeps)
  //
  ///////////////////////////////
  if(doflux && truestep){
    // must come after metric changes that can change where flux surface is since otherwise when flux surface changes, we won't track this substep's flux through the new surface but the old surface (which could even be at r=0 and have no flux)
    // if using unew, then since metric update above uses old unew, need present flux at new horizon surface
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
        diag_flux_pureflux(pb,F1, F2, F3, fluxdt); // fluxdt is true dt for this flux as added in dUtoU() as part of ucum update
#endif
      }
  }







  ///////////////////////////////////////
  //
  // Copy over tempucum -> ucum per pl to account for staggered field
  //
  // and volume differentiate the conserved quantity
  //
  ///////////////////////////////////////
  if(finalstep){
    // last call, so ucum is cooked and ready to eat!
    if(FLUXB==FLUXCTSTAG){
      copy_tempucum_finalucum(DOALLPL,Uconsevolveloop,tempucum,ucum);
    }
    // useducum only used after ucum assigned to ui above
    utoinvert=ucum;
    useducum=ucum;
    //    ubound=utoinvert;
  }
  else{ // otherwise still iterating on primitives
    // don't need to copy over tempucum->ucum yet
    utoinvert=uf;
    useducum=tempucum;
    //    ubound=utoinvert;
  }

  // GODMARK: use of globals
  myupoint=GLOBALPOINT(upointglobal);



  // to debug ENO
  //#if(DOENODEBUG)
  //debugeno_compute(pb,utoinvert,enodebugarray); //SASDEBUG -- OLD USE: now assign debugs inside reconstructeno code
  //#endif





  // a2c_1 a2c_2 a2c_3
  if(docons){
    int pl,pliter;
    // conserved quantity is limited later after primitive is obtained
    PALLLOOP(pl) locpl[pl]=CENT;
    PALLLOOP(pl) whichpltoavg[pl]=do_conserved_integration[pl];// default
    PALLLOOP(pl) ifnotavgthencopy[pl]=1-do_conserved_integration[pl];// default
    avg2cen_interp(locpl,whichpltoavg, ifnotavgthencopy, ENOCONSERVED, ENOAVG2CENTTYPE, pb, utoinvert, myupoint);  //SASMARK:  pb's for shock indicators should be defined on ALL grid, not only on ghost+active.  Maybe should use pi instead because define everywhere?
  }
  else{
    //    myupoint=utoinvert;
    copy_tempucum_finalucum(DOALLPL,Uconsevolveloop,utoinvert,myupoint);
  }



  ////////////////////////////
  //
  // split CZLOOP above and below to allow staggered field method
  //
  ////////////////////////////
  if(FLUXB==FLUXCTSTAG){
    // if using staggered grid for magnetic field, then need to convert ucum to pstag to new pb/pf

    // GODMARK: If had c2a/a2c with 3-point outputs, then could do avg2cen_interp and below at once

    // first pb entry is used for shock indicator, second is filled with new field
    interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstag, myupoint, pf);

    ////////////////////    
    // now utoinvert contains centered point conserved quantities ready for inversion
    ////////////////////

  }




  /////////////////////////////////
  //
  // Invert U -> p
  //
  // and fixup p if inversion bad
  // and compute dissipation rate if requested
  //
  /////////////////////////////////


  // get loop range
  get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);

#pragma omp parallel OPENMPGLOBALPRIVATEFORINVERSION
  {
    int i,j,k,pl,pliter;
    FTYPE prbefore[NPR];
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
    int showmessages=1;
    int allowlocalfailurefixandnoreport=1; // allow local fixups
    
    // setup default eomtype
    int eomtype;


    OPENMP3DLOOPVARSDEFINE;  OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    // ptr different when in parallel region
    ptrgeom=&geomdontuse;

#pragma omp for schedule(OPENMPVARYENDTIMESCHEDULE(),OPENMPCHUNKSIZE(blocksize)) reduction(+: nstroke)
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      // setup default eomtype (KORALTODO: GODMARK: Doesn't use source() version of eomtype from above)
      eomtype=EOMDEFAULT;
 
      // set geometry for centered zone to be updated
      get_geometry(i, j, k, CENT, ptrgeom);


      // store guess for diss_compute before changed by normal inversion
      PALLLOOP(pl) prbefore[pl]=MACP0A1(pf,i,j,k,pl);
  

      // invert point U-> point p
      int dissmeasure=-1.0; // assume ok to try energy
      int whichcap=CAPTYPEBASIC;
      int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
      int modprim=0;
      int checkoninversiongas=CHECKONINVERSION;
      int checkoninversionrad=CHECKONINVERSIONRAD;
      MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM, UEVOLVE, MAC(myupoint,i,j,k), NULL, ptrgeom, dissmeasure, MAC(pi,i,j,k), MAC(pf,i,j,k),&newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
      nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;

      //If using a high order scheme, need to choose whether to trust the point value
      if(docons){
        MYFUN(check_point_vs_average(timeorder, numtimeorders, GLOBALMAC(pflag,i,j,k),MAC(pb,i,j,k),MAC(pf,i,j,k),MAC(myupoint,i,j,k),MAC(utoinvert,i,j,k),ptrgeom,&newtonstats),"advance.c:advance_finitevolume()", "check_point_vs_average()", 1);
      }


#if(DODISS||DODISSVSR)
      if(finalstep){
        // then see what entropy inversion would give
        diss_compute(EVOLVEUTOPRIM,UEVOLVE,MAC(useducum,i,j,k),ptrgeom,prbefore, MAC(pf,i,j,k),&newtonstats);
      }
#endif


#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
        {
          // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
          // SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
          if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(finalstep)){
            int docorrectucons=0;
            MYFUN(fixup1zone(docorrectucons,MAC(pf,i,j,k), MAC(useducum,i,j,k), ptrgeom, finalstep),"advance.c:advance_finitevolume()", "fixup1zone()", 1);
          }
#endif
        }
    }// end COMPZLOOP
  }// end parallel region





  /////////////////////////////////
  //
  // If not fixing up primitives after inversion immediately, then fix up all zones at once afterwards
  //
  /////////////////////////////////

#if(SPLITNPR)
  // don't update metric if only doing B1-B3
  if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
#if(FIXUPZONES==FIXUPALLZONES)
      fixup(stage,pf,useducum,finalstep);
#endif  
    }



  /////////////////////////////////
  //
  // Determine next timestep from waves, fluxes, and source updates
  //
  /////////////////////////////////

  prepare_globaldt(truestep,ndt1,ndt2,ndt3,accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,ndt);
 

#if(PRODUCTION==0)
  trifprintf( "2f");
#endif

  return (0);
}





/// some dt calculations done at end of each substep
static int prepare_globaldt(
                            int truestep,
                            FTYPE ndt1,FTYPE ndt2,FTYPE ndt3,
                            FTYPE accdt,int accdti,int accdtj,int accdtk,
                            FTYPE gravitydt,int gravitydti,int gravitydtj,int gravitydtk,
                            FTYPE *ndt)
{
  FTYPE wavedt;
  int jj;



  ////////////////
  //
  // unsplit multidimensional Courant condition
  //
  ////////////////
  if(PERCELLDT==0){
    wavedt = MINDTSET(ndt1,ndt2,ndt3);
  }
  else{
    wavedt = ndt1; // full result stored in *any* of ndt1,2,3 regardless of whether dimension is relevant
  }



  ///////////////
  //
  // minimize dt over all "operators"
  //
  ///////////////
  *ndt = defcon * MIN(wavedt , accdt);

#if(USEGRAVITYDTINDTLIMIT)
  // use gravitydt (often too small, but sometimes accdt/ndt not small enough)
  *ndt = MIN(*ndt,gravitydt);
#endif


  /////////
  //
  // Store some global variables for timestep
  //
  /////////

  gravitydtglobal = gravitydt;
  sourcedtglobal  = accdt; // accdt includes gravitydtglobal
  wavedtglobal    = wavedt;

#if(PRODUCTION==0)
  if(truestep){
    // report per-CPU time-step limited every 100 time steps

    // GODMARK: 1 : do always
    if(1|| nstep%DTr==0){
      logdtfprintf("nstep=%ld steppart=%d :: dt=%g ndt=%g ndt1=%g ndt2=%g ndt3=%g\n",nstep,steppart,dt,*ndt,ndt1,ndt2,ndt3);
      SLOOPA(jj) logdtfprintf("dir=%d wavedti=%d wavedtj=%d wavedtk=%d\n",jj,waveglobaldti[jj],waveglobaldtj[jj],waveglobaldtk[jj]);
      logdtfprintf("accdt=%g (accdti=%d accdtj=%d accdtk=%d) :: gravitydt=%g (gravitydti=%d gravitydtj=%d gravitydtk=%d) :: gravityskipstep=%d\n",accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,gravityskipstep);
    }
  }
#endif


  return(0);
}











/// check whether point conserved quantity inverted successfully to point primitive.
///   if unsuccessful, then see if want to revert to average conserved quantity and invert that
///   if Uavg->p unsuccessful, then leave as failure
/// if Upoint->p is good, then check if p from Upoint is much different than p from Uavg.  If so, limit change
/// upoint only needed for diagnostics
static int check_point_vs_average(int timeorder, int numtimeorders, PFTYPE *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats)
{
  FTYPE pavg[NPR];  //atch for temporary storage of primitives obtained from inverting the averaged conserved quantities
  PFTYPE invert_from_point_flag, invert_from_average_flag;
  FTYPE frac_avg_used;  //this is be used for flux interpolation limiting
  int pl,pliter;
  int is_convergence_failure;
  int avgschemeatall;
  int finalstep;
  FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout );
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; // allow local fixups

  // setup default eomtype
  int eomtype=EOMDEFAULT;


  finalstep=timeorder == numtimeorders-1;



  avgschemeatall=(interporder[avgscheme[1]]>3) ||  (interporder[avgscheme[2]]>3) ||  (interporder[avgscheme[3]]>3);
  if(avgschemeatall==0) return(0); // since nothing to do


  invert_from_point_flag = lpflag[FLAGUTOPRIMFAIL];


  if( 0 && debugfail >= 1 && (IFUTOPRIMFAILSOFT(invert_from_point_flag)) ) {
    dualfprintf( fail_file, "t = %g, nstep = %ld, steppart = %d, i = %d, j = %d, rho = %21.15g, u = %21.15g, fracneg = %21.15g\n",
                 t, realnstep, steppart, ptrgeom->i + startpos[1], ptrgeom->j + startpos[2],
                 pf[RHO], pf[UU], (pf[RHO]>0)?(-pf[UU]/(pf[RHO]+DBL_MIN)):(-pf[RHO]/(pf[UU]+DBL_MIN)) );
  }


  //WHAT IF INTERNAL ENERGY BECOMES SLIGHTLY NEGATIVE?  WE STILL CAN DO THE LIMITING IN PRIM QUANTITIES! -- coorrected but check! -- SUPERSASMARK TODO atch
  if( LIMIT_AC_PRIM_FRAC_CHANGE &&
      (
       IFUTOPRIMNOFAILORFIXED(invert_from_point_flag) || //atch added the below to still do the pt. vs. avg. check on primitives if the internal energy goes neg.
       ( (IFUTOPRIMFAILSOFTNOTRHORELATED(invert_from_point_flag)) && (0 != STEPOVERNEGU) ) || //intermediate substep with stepping over u < 0
       ( (IFUTOPRIMFAILSOFTRHORELATED(invert_from_point_flag)) && (0 != STEPOVERNEGRHO) ) //intermediate substep with stepping over rho < 0
       )
      ) {


    //make a copy of the initial guess so that not to modify the original pb's
    PLOOP(pliter,pl) pavg[pl] = pb[pl];
    //invert the average U -> "average" p
    int dissmeasure=-1.0; // assume ok to try energy
    int whichcap=CAPTYPEBASIC;
    int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
    int modprim=0;
    int checkoninversiongas=CHECKONINVERSION;
    int checkoninversionrad=CHECKONINVERSIONRAD;
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM, UEVOLVE, uavg, NULL, ptrgeom, dissmeasure, pavg, pavg,newtonstats),"step_ch.c:advance()", "Utoprimgen", 3);

    invert_from_average_flag = GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);

    //Inversion from the average value succeeded or has a negative density or internal energy
    if( IFUTOPRIMNOFAILORFIXED(invert_from_average_flag) || IFUTOPRIMFAILSOFT(invert_from_average_flag) ) {
      //Inversion from both from the point and the average values succeeded
      //checks if the states' gamma factors and densities are different by more than a certain fraction
      //and if different, modify the point values such that they are not further than by MAX_AC_PRIM_FRAC_CHANGE
      //from the average ones
      frac_avg_used = limit_prim_correction(MAX_AC_PRIM_FRAC_CHANGE, ptrgeom, pavg, pf);

#if( DOENODEBUG )  //atch debug; note that use the location with pl =0 , interporflux = 0, & dir = 1 since limiting the change in prim quantities is not a per-component operation
      if(  frac_avg_used > 0.01 ) {
        MACP0A4(enodebugarray,ptrgeom->i,ptrgeom->j,ptrgeom->k,1-1,0,0,ENODEBUGPARAM_LIMCORR_PRIM)++;
      }
#endif
      
      if( pf[UU] < 0.0 && timeorder == numtimeorders-1 ) {
        lpflag[FLAGUTOPRIMFAIL] = UTOPRIMFAILU2AVG2;
      }

      //lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;  //unneeded since it is alrady == UTOPRIMNOFAIL

    } // end if both point and average did NOT fail
    else {
      //Inversion from the point values succeeded but that from the average failed:
      //retain the point value

      //set the inversion error flag to correspond to the inversion from the average value
      lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;
      frac_avg_used = 0.0;  //used point value, i.e. zero fracion of the average value
    }
  }
  else if( INVERTFROMAVERAGEIFFAILED && (IFUTOPRIMFAIL(invert_from_point_flag)) ) {  //failure  //atch correct
    //inversion from the point value failed

    // if last substep -> revert to the average value, else if only negative densities then allow for substep.  If other type of failures, then never allow and revert to Uavg


    // GODMARK:
    // CHECK IF INVERSION FAILED.  IF FAILED, TRY TO USE uavg:
    // invert U->p

    // GODMARK: decide whether best to revert to average or not

    //=1 if it is a non-convergence failure; =0 if it is an occurrence of a negative density
    is_convergence_failure = !IFUTOPRIMFAILSOFT(invert_from_point_flag);

    if((timeorder==numtimeorders-1 /*&& (1 == is_nonconvergence_failure || FIXUPZONES == FIXUPNOZONES)*/ ) ||   // last substep, then DO invert from average IF no fixups or non-convergence failure (ADT) SASMARKx

       ( 1 == is_convergence_failure ) || //non-convergence (no solution in primitive quantities) error, so have to fix it up

       ( (timeorder<numtimeorders-1) && (
                                         ((IFUTOPRIMFAILSOFTNOTRHORELATED(invert_from_point_flag)) && (0 == STEPOVERNEGU))||
                                         ((invert_from_point_flag==UTOPRIMFAILRHOUNEG) && (0 == STEPOVERNEGRHOU)) ||
                                         ((IFUTOPRIMFAILSOFTRHORELATED(invert_from_point_flag)) && (0 == STEPOVERNEGRHO))
                                         )
         )  //intermediate substep with no stepping over u < 0, rho<0, or both <0
       ) {
      if(debugfail >= 1) {
        dualfprintf( fail_file, "Inversion from the point value failed.  Using the inversion from the average value.\n" );
      }

      //make a copy of the initial guess so that not to modify the original pb's
      PLOOP(pliter,pl) pf[pl] = pb[pl];
      //invert the average U -> "average" p
      int dissmeasure=-1.0; // assume ok to try energy
      int whichcap=CAPTYPEBASIC;
      int whichmethod=MODEPICKBEST; // try to choose best option for this "external" inversion
      int modprim=0;
      int checkoninversiongas=CHECKONINVERSION;
      int checkoninversionrad=CHECKONINVERSIONRAD;
      MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM, UEVOLVE, uavg, NULL, ptrgeom, dissmeasure, pb, pf,newtonstats),"step_ch.c:advance()", "Utoprimgen", 3);
      //      invert_from_average_flag = lpflag[FLAGUTOPRIMFAIL];


      //Have the results from the inversion from the average value -- copy the result over
      //      PLOOP(pliter,pl) pf[pl] = pavg[pl];
      //      lpflag[FLAGUTOPRIMFAIL] = invert_from_average_flag;

      //old code:
      //MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,EVOLVEUTOPRIM, UEVOLVE, avg, NULL, ptrgeom, dissmeasure, pb, pf,&newtonstats),"step_ch.c:advance()", "Utoprimgen", 2);

      frac_avg_used = 1.0; //reverted to the average value

    }
    else {
      frac_avg_used = 0.0; //used the point value
    }
  }

  return(0);
}





#define COMPARE_GAMMA 0

///If density or gamma-factors are different by more than fractional_difference_threshold for states pin & pout, 
///if different -- correct pout such that it is not more than fractional_difference_threshold away from pin.
/// externally referenced
FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout )
{
  FTYPE gammain = 0.0, gammaout = 0.0;
  FTYPE frac_start, frac_end, frac_diff;
  FTYPE fraction_input_value;
  FTYPE comovingenergyin, comovingenergyout;
  int pl,pliter;
  FTYPE bsqin,bsqout;
  struct of_state qin, qout;
  int jj;
  FTYPE bdotuin, bdotuout;


#if( COMPARE_GAMMA ) 
  gamma_calc( pin, geom, &gammain );
  gamma_calc( pout, geom, &gammaout );
#endif

  get_state(pin, geom, &qin);
  get_state(pout, geom, &qout);

  bsqout = dot(qout.bcon, qout.bcov);
  bsqin = dot(qin.bcon, qin.bcov);

  bdotuin = 0.0;   DLOOPA(jj) bdotuin+=(qin.ucov[jj])*(qin.bcon[jj]);
  bdotuout = 0.0;  DLOOPA(jj) bdotuout+=(qout.ucov[jj])*(qout.bcon[jj]);

  // u.T.u comoving energy density
  comovingenergyin = pin[RHO] + pin[UU] + bsqin*0.5 - bdotuin*bdotuin;
  comovingenergyout = pout[RHO] + pout[UU] + bsqout*0.5 - bdotuout*bdotuout;


 
#if( COMPARE_GAMMA ) 
  frac_diff = MAX( fractional_diff(gammain, gammaout), 
                   fractional_diff( comovingenergyin, comovingenergyout ) );
#else
  frac_diff = fractional_diff( comovingenergyin, comovingenergyout );
#endif

  //fractional difference at which the reduction to the input value starts
  frac_start = 0.5 * fractional_difference_threshold;

  //fractional difference after which only the input value is used
  frac_end = fractional_difference_threshold;

  //the fraction of the input value used in the output; increases linearly from 0 to 1 for fractions going from frac_start to frac_end
  fraction_input_value = MAX( 0., MIN(1., (frac_diff - frac_start)/(frac_end - frac_start) ) );

  if( 0.0 != fraction_input_value ){
    //states are too different: reverted to primitives that correspond to average conserved quantities because trust them more than point values
    dualfprintf( fail_file, "States are too different, using %3d%% of the average values: i = %d, j = %d, k = %d, nstep = %ld, steppart = %d, t = %21.15g\n", 
                 (int)(100. * fraction_input_value), geom->i, geom->j, geom->k, nstep, steppart, t );
    if( debugfail >= 2 ){
      dualfprintf( fail_file, "Prim. pt. value (gamma, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n",  gammaout, pout[RHO], pout[UU] );
      dualfprintf( fail_file, "Prim. avg value (gamma, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", gammain, pin[RHO], pin[UU] );
      dualfprintf( fail_file, "Frac. difference(ganna, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", 
                   fractional_diff(gammain, gammaout),
                   fractional_diff(pin[RHO], pout[RHO]),  
                   fractional_diff(pin[UU], pout[UU])
                   );
    }
  }

  PLOOP(pliter,pl) {
    pout[pl] = fraction_input_value * pin[pl] + (1. - fraction_input_value) * pout[pl];
  }

  return( fraction_input_value );
}



///Returns the fractional difference between a & b
static FTYPE fractional_diff( FTYPE a, FTYPE b )
{
  FTYPE frac_diff;

  frac_diff = 2. * fabs( a - b ) / ( fabs(a) + fabs(b) + DBL_MIN );

  return( frac_diff );

}

























/// get dUavg
static void flux2dUavg(int whichpl, int i, int j, int k, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE *dU1avg,FTYPE *dU2avg,FTYPE *dU3avg)
{
  FTYPE idx1,idx2,idx3;
  int pl,pliter;

#if(VOLUMEDIFF==0)
  idx1=1.0/dx[RR];
  idx2=1.0/dx[TH];
  idx3=1.0/dx[PH];
#else
  idx1=GLOBALMETMACP0A1(idxvol,i,j,k,RR);
  idx2=GLOBALMETMACP0A1(idxvol,i,j,k,TH);
  idx3=GLOBALMETMACP0A1(idxvol,i,j,k,PH);
#endif

  int special;
  if(whichpl==DOSPECIALPL){
    special=NSPECIAL;
    whichpl=DOALLPL; // override
  }
  else special=0;

  // initialize for simplicity so don't have to do it conditionally on N>1
  PALLLOOPSPECIAL(pl,special){
    dU1avg[pl]=0;
    dU2avg[pl]=0;
    dU3avg[pl]=0;
  }

  if(FLUXB==FLUXCD){ // don't use volume reg. since differencing is large
    
    if(whichpl==DOALLPL || whichpl==DONONBPL){

      PLOOPNOB1(pl){
#if(N1>1)
        dU1avg[pl]=(
                    - (MACP0A1(F1,ip1mac(i),j,k,pl) - MACP0A1(F1,i,j,k,pl)) *idx1
                    );
#endif
#if(N2>1)
        dU2avg[pl]=(
                    - (MACP0A1(F2,i,jp1mac(j),k,pl) - MACP0A1(F2,i,j,k,pl)) *idx2
                    );
#endif

#if(N3>1)
        dU3avg[pl]=(
                    - (MACP0A1(F3,i,j,kp1mac(k),pl) - MACP0A1(F3,i,j,k,pl)) *idx3
                    );
#endif

      }
    
      // rest of variables (if any) are normal
      PLOOPNOB2SPECIAL(pl,special){
#if(N1>1)
        dU1avg[pl]=(
                    - (MACP0A1(F1,ip1mac(i),j,k,pl) - MACP0A1(F1,i,j,k,pl)) *idx1
                    );
#endif
#if(N2>1)
        dU2avg[pl]=(
                    - (MACP0A1(F2,i,jp1mac(j),k,pl) - MACP0A1(F2,i,j,k,pl)) *idx2
                    );
#endif
#if(N3>1)
        dU3avg[pl]=(
                    - (MACP0A1(F3,i,j,kp1mac(k),pl) - MACP0A1(F3,i,j,k,pl)) *idx3
                    );
#endif

      }
    }// end if doing non-B U

    if(whichpl==DOALLPL || whichpl==DOBPL){

      // simple version that assumes Fi[Bi] is set to 0 in flux.c for FLUXCD (which it is currently)
      PLOOPBONLY(pl){
#if(N1>1)
        dU1avg[pl]=(
                    - (MACP0A1(F1,ip1mac(i),j,k,pl) - MACP0A1(F1,im1mac(i),j,k,pl)) *idx1
                    );
#endif
#if(N2>1)
        dU2avg[pl]=(
                    - (MACP0A1(F2,i,jp1mac(j),k,pl) - MACP0A1(F2,i,jm1mac(j),k,pl)) *idx2
                    );
#endif
#if(N3>1)
        dU3avg[pl]=(
                    - (MACP0A1(F3,i,j,kp1mac(k),pl) - MACP0A1(F3,i,j,km1mac(k),pl)) *idx3
                    );
#endif
      }
    }// end if doing B U

  }// end FLUXCD

  else{

    
    FTYPE lowerF2,upperF2;
#if(IF3DSPCTHENMPITRANSFERATPOLE==1)
    // Ensure that spherical polar axis flux only contributes to B2 for special3dspc=1 with polar cut-out.
    // Allows more arbitrary interpolations across axis (including with gdet factors)
    // "lower" means we will zero-out lower F2, while "upper" means we will zero-out upper F2 in terms of j in Du2avg expression
    int preconditionF2lower=(ISSPCMCOORD(MCOORD) && FLUXB==FLUXCTSTAG && special3dspc && (startpos[2]+j==0 || startpos[2]+j==totalsize[2]));
    int preconditionF2upper=(ISSPCMCOORD(MCOORD) && FLUXB==FLUXCTSTAG && special3dspc && (startpos[2]+j==-1 || startpos[2]+j==totalsize[2]-1));
#endif

    // other (normal) FLUXB methods, including FLUXCTSTAG
    PALLLOOPSPECIAL(pl,special) {
      if(whichpl==DONONBPL && BPL(pl)==1 || whichpl==DOBPL && BPL(pl)==0) continue;


#if(IF3DSPCTHENMPITRANSFERATPOLE==1)
      if(pl!=B2 && preconditionF2lower) lowerF2=0.0;
      else lowerF2=1.0;
      if(pl!=B2 && preconditionF2upper) upperF2=0.0;
      else upperF2=1.0;
#else
      lowerF2=1.0;
      upperF2=1.0;
#endif


#if(N1>1)
      dU1avg[pl]=(
                  - (MACP0A1(F1,ip1mac(i),j,k,pl) - MACP0A1(F1,i,j,k,pl)) *idx1
                  );
#endif
#if(N2>1)
      dU2avg[pl]=(
                  - (upperF2*MACP0A1(F2,i,jp1mac(j),k,pl) - lowerF2*MACP0A1(F2,i,j,k,pl)) *idx2
                  );
#endif
#if(N3>1)
      dU3avg[pl]=(
                  - (MACP0A1(F3,i,j,kp1mac(k),pl) - MACP0A1(F3,i,j,k,pl)) *idx3
                  );
#endif
    }




  }




}





/// convert point versions of U_i^{n} and dU -> U_i^{n+1} and other versions
static void dUtoU(int timeorder, int whichpl, int i, int j, int k, int loc, FTYPE *dUgeom, FTYPE (*dUcomp)[NPR], FTYPE *dUriemann, FTYPE *CUf, FTYPE *CUnew, FTYPE *Ui,  FTYPE *Uf, FTYPE *ucum)
{
  int pl,pliter;
  void dUtoU_check(int timeorder, int i, int j, int k, int loc, int pl, FTYPE *dUgeom, FTYPE (*dUcomp)[NPR], FTYPE *dUriemann, FTYPE *CUf, FTYPE *CUnew, FTYPE *Ui,  FTYPE *Uf, FTYPE *ucum);

  int special;
  if(whichpl==DOSPECIALPL){
    special=NSPECIAL;
    whichpl=DOALLPL; // override
  }
  else special=0;



  // get physics non-geometry source terms
  FTYPE dUrad[NPR+NSPECIAL],dUnonrad[NPR+NSPECIAL]; // NSPECIAL part only used if whichpl==DOSPECIALPL
  int sc;
  PALLLOOPSPECIAL(pl,special){
    dUrad[pl]=0.0; // init as zero
    // only sc=RADSOURCE is in implicit part, so only separate that out
    sc=RADSOURCE;
    if(pl<NPR) dUrad[pl] = dUcomp[sc][pl]; // dUcomp always NPR, but dUrad over full range possible

    // all terms except RADSOURCE
    dUnonrad[pl] = dUgeom[pl] - dUrad[pl];
  }

  // initialize dUradall as zero.  Ensure initialized for any possible pl's
  FTYPE dUradall[NPR+NSPECIAL][MAXTIMEORDER];
  int ii,jj;
  for(ii=0;ii<MAXTIMEORDER;ii++){
    PLOOP(pliter,pl) dUradall[pl][ii]=0.0;
    PALLLOOPSPECIAL(pl,special) dUradall[pl][ii]=0.0;
    PLOOPBONLY(pl) dUradall[pl][ii]=0.0;
  }


  // store physics dU's that could need an implicit treatment
  if(EOMRADTYPE!=EOMRADNONE && special==0){ // don't access Mradk if from compute_dissmeasure that uses special!=0
    // only non-field case
    if(whichpl==DOALLPL || whichpl==DONONBPL){
      PLOOP(pliter,pl){
        GLOBALMACP1A1(Mradk,timeorder,i,j,k,pl) = dUrad[pl]; // USE OF GLOBALS // only NPR pl's exist for Mradk
      }
      
      // now assign dUradall for radiation
      for(ii=0;ii<=timeorder;ii++){ // timeorder goes from 0..TIMEORDER-1 inclusive
        PLOOP(pliter,pl) dUradall[pl][ii] = GLOBALMACP1A1(Mradk,ii,i,j,k,pl); // radiation could affect any pl's, but only NPR pl's exist for Mradk GODMARK -- issue with dissmeasure?
      }
    }
  }
  else{
    // then assume not related to implicit radiation but still done per timeorder as if explicit
    ii=timeorder;
    PALLLOOPSPECIAL(pl,special) dUradall[pl][ii] = dUrad[pl];
  }


  if(whichpl==DOALLPL){
    // finally assign new Uf and ucum
    // store uf to avoid recomputing U(pf) used later as pb for advance()
    PALLLOOPSPECIAL(pl,special) Uf[pl] = UFSET(CUf,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);
    
    
    // how much of Ui, dU, and Uf to keep for final solution
    // ultimately ucum is actual solution used to find final pf
    PALLLOOPSPECIAL(pl,special) ucum[pl] += UCUMUPDATE(CUnew,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);

#if(PRODUCTION==0)
    if(FLUXB!=FLUXCTSTAG){// turned off by default for FLUXB==FLUXCTSTAG since even with PRODUCTION==0, FLUXB==FLUXCTSTAG's extended loop causes output at edges.
      PLOOP(pliter,pl) dUtoU_check(timeorder,i,j,k,loc,pl, dUgeom, dUcomp, dUriemann, CUf, CUnew, Ui,  Uf, ucum);
    }
#endif

  }
  else if(whichpl==DOBPL){
    PLOOPBONLY(pl) Uf[pl] = UFSET(CUf,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);
    PLOOPBONLY(pl) ucum[pl] += UCUMUPDATE(CUnew,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);
  }
  else if(whichpl==DONONBPL){
    PALLLOOPSPECIAL(pl,special) if(!BPL(pl)) Uf[pl] = UFSET(CUf,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);
    PALLLOOPSPECIAL(pl,special) if(!BPL(pl)) ucum[pl] += UCUMUPDATE(CUnew,dt,Ui[pl],Uf[pl],dUriemann[pl],dUnonrad[pl],dUradall[pl]);
  }

  //  if(nstep==4 && steppart==0 && whichpl==DONONBPL){
  //    PALLLOOPSPECIAL(pl,special) dualfprintf(fail_file,"UtoU: %21.15g %21.15g %21.15g %21.15g\n",(CUf[0])*(Ui[pl]/globalgeom[pl]),(CUf[1])*(Uf[pl]/globalgeom[pl]),(CUf[2])*(dt)*((dUriemann[pl]/globalgeom[pl])+(dUgeom[pl]/globalgeom[pl])),CUf[2]*dt);
  //  }


}



/// Check result of dUtoU()
static void dUtoU_check(int timeorder, int i, int j, int k, int loc, int pl, FTYPE *dUgeom, FTYPE (*dUcomp)[NPR], FTYPE *dUriemann, FTYPE *CUf, FTYPE *CUnew, FTYPE *Ui,  FTYPE *Uf, FTYPE *ucum)
{
  int showfluxes;
  void show_fluxes(int i, int j, int k, int loc, int pl,FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);


  // default
  showfluxes=0;

  if(!isfinite(Uf[pl])){
    dualfprintf(fail_file,"dUtoU after: nan found for Uf[%d]=%21.15g\n",pl,Uf[pl]);
    dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
    showfluxes=1;
  }
  if(!isfinite(ucum[pl])){
    dualfprintf(fail_file,"dUtoU after: nan found for ucum[%d]=%21.15g\n",pl,ucum[pl]);
    dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
    showfluxes=1;
  }

  if(showfluxes){
    show_fluxes(i,j,k,loc,pl,GLOBALPOINT(F1),GLOBALPOINT(F2),GLOBALPOINT(F3));
  }


}


/// Check ucum
/// externally referenced.
void ucum_check(int i, int j, int k, int loc, int pl, FTYPE *ucum)
{
  int showfluxes;


  // default
  showfluxes=0;

  if(!isfinite(ucum[pl])){
    dualfprintf(fail_file,"ucum_check: nan found at i=%d j=%d k=%d loc=%d for ucum[pl=%d]=%21.15g\n",i,j,k,loc,pl,ucum[pl]);
    showfluxes=1;
  }

  if(showfluxes){
    show_fluxes(i,j,k,loc,pl,GLOBALPOINT(F1),GLOBALPOINT(F2),GLOBALPOINT(F3));
  }


}

/// debug show fluxes
static void show_fluxes(int i, int j, int k, int loc, int pl,FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  if(N1>1){
    dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d :: F1[i]=%21.15g F1[i+1]=%21.15g dF1/dx1=%21.15g \n",pl,i,j,k,MACP0A1(F1,i,j,k,pl),MACP0A1(F1,i+SHIFT1,j,k,pl),(MACP0A1(F1,i+SHIFT1,j,k,pl)-MACP0A1(F1,i,j,k,pl))/dx[1]);
  }
  if(N2>1){
    dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d :: F2[j]=%21.15g F2[j+1]=%21.15g dF2/dx2=%21.15g\n",pl,i,j,k,MACP0A1(F2,i,j,k,pl),MACP0A1(F2,i,j+SHIFT2,k,pl),(MACP0A1(F2,i,j+SHIFT2,k,pl)-MACP0A1(F2,i,j,k,pl))/dx[2]);
  }
  if(N3>1){
    dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d :: F3[k]=%21.15g F3[k+1]=%21.15g dF3/dx3=%21.15g \n",pl,i,j,k,MACP0A1(F3,i,j,k,pl),MACP0A1(F3,i,j,k+SHIFT3,pl),(MACP0A1(F3,i,j,k+SHIFT3,pl)-MACP0A1(F3,i,j,k,pl))/dx[3]);
  }

}



/// find global dt.  Useful if needed not during evolution, such as at t=0 or for restarting the run if restarting finished run that has a generally smaller dt than should use (including possibly dt=0)
/// externally referenced
int set_dt(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE *dt)
{
  struct of_state state;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int i,j,k;
  int pl,pliter;
  int jj,kk,ll;
  int dir,ignorecourant;
  FTYPE cmax1,cmin1;
  FTYPE cmax2,cmin2;
  FTYPE cmax3,cmin3;
  int wavendti[NDIM],wavendtj[NDIM],wavendtk[NDIM];
  int accdti,accdtj,accdtk;
  int gravitydti,gravitydtj,gravitydtk;
  FTYPE tempwavedt,tempaccdt,tempgravitydt;
  FTYPE dtij[NDIM]={BIG}, wavedt, accdt, gravitydt;
  FTYPE wavendt[NDIM];
  FTYPE wavedt_1,wavedt_2,wavedttemp;
  FTYPE ndtfinal;
  FTYPE dUgeom[NPR],dUcomp[NUMSOURCES][NPR];
  int enerregion;
  int *waveenerpos, *sourceenerpos;
  FTYPE X[NDIM],V[NDIM],Vp1[NDIM];
  FTYPE dUriemann[NPR];
  FTYPE Ugeomfree[NPR],U[NPR];




  wavedt_1=wavedt_2=wavedt=accdt=gravitydt=ndtfinal=BIG;
  wavendt[1]=wavendt[2]=wavendt[3]=BIG;
  wavendti[1]=wavendtj[1]=wavendtk[1]=-100;
  wavendti[2]=wavendtj[2]=wavendtk[2]=-100;
  wavendti[3]=wavendtj[3]=wavendtk[3]=-100;
  accdti=accdtj=accdtk=-100;
  gravitydti=gravitydtj=gravitydtk=-100;

  enerregion=OUTSIDEHORIZONENERREGION; // consistent with flux update (except when doing WHAM)
  sourceenerpos=enerposreg[enerregion];

  int didreturnpf=0; // used to have source() pass back better guess for Utoprimgen() than pb.

  //  COMPFULLLOOP{ // want to use boundary cells as well to limit dt (otherwise boundary-induced changes not treated)
  LOOPWITHINACTIVESECTIONEXPAND1(i,j,k){ // only 1 grid cell within boundary.  Don't need to use extrapolated cells deep inside boundary because those cells are not evolved.  More consistent with how flux.c sets dt.

    // includes "ghost" zones in case boundary drives solution
    get_geometry(i, j, k, CENT, ptrgeom);
    
    // need full state for vchar()
    MYFUN(get_state(MAC(prim,i,j,k), ptrgeom, &state),"advance.c:set_dt()", "get_state()", 1);
    
    
#if(N1>1)
    dir=1;
    MYFUN(vchar_all(MAC(prim,i,j,k), &state, dir, ptrgeom, &cmax1, &cmin1,&ignorecourant),"advance.c:set_dt()", "vchar_all() dir=1", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax1),fabs(cmin1));
    if(FORCESOLVEL){
      dtij[dir] = cour*dx[dir]*sqrt(ptrgeom->gcov[GIND(dir,dir)]);
    }
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif
    
#if(N2>1)
    dir=2;
    MYFUN(vchar_all(MAC(prim,i,j,k), &state, dir, ptrgeom, &cmax2, &cmin2,&ignorecourant),"advance.c:set_dt()", "vchar_all() dir=2", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax2),fabs(cmin2));
    if(FORCESOLVEL){
      dtij[dir] = cour*dx[dir]*sqrt(ptrgeom->gcov[GIND(dir,dir)]);
    }
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif

#if(N3>1)
    dir=3;
    MYFUN(vchar_all(MAC(prim,i,j,k), &state, dir, ptrgeom, &cmax3, &cmin3,&ignorecourant),"restart.c:set_dt()", "vchar_all() dir=3", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax3),fabs(cmin3));
    if(FORCESOLVEL){
      dtij[dir] = cour*dx[dir]*sqrt(ptrgeom->gcov[GIND(dir,dir)]);
    }
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif



    // also get PERCELLDT==1 version
    wavedttemp = MINDTSET(dtij[1],dtij[2],dtij[3]);
    if(wavedttemp<wavedt_2) wavedt_2=wavedttemp;



    ////////////////////
    //
    // deal with source terms
    //
    //////////////////////
    if(WITHINENERREGION(sourceenerpos,i,j,k)){


#if(LIMITDTWITHSOURCETERM)

      // conserved quantity without geometry
      MYFUN(primtoU(UEVOLVE, MAC(prim,i,j,k), &state, ptrgeom, U, NULL),"step_ch.c:advance()", "primtoU()", 1);
      PLOOP(pliter,pl) Ugeomfree[pl] = U[pl]*(ptrgeom->IEOMFUNCNOSINGMAC(pl));

      // get source term
      // GODMARK: here dUriemann=0, although in reality this setting of dt is related to the constraint trying to make
      PLOOP(pliter,pl) dUriemann[pl]=0.0;
      FTYPE CUf[NUMDTCUFS]={0.0}; // no update yet!
      int timeorder=0; // fake timeorder
      FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];
      // modifies prim() to be closer to final, which is ok here.
      // setup default eomtype
      int eomtype=EOMDEFAULT;
      MYFUN(source(MAC(prim,i,j,k), MAC(prim,i,j,k), MAC(prim,i,j,k), &didreturnpf, &eomtype, ptrgeom, &state, U, U, CUf, CUimp, 0.0, dUriemann, dUcomp, dUgeom),"advance.c:set_dt()", "source", 1);

      // get dt limit
      compute_dt_fromsource(ptrgeom,&state,MAC(prim,i,j,k), Ugeomfree, dUgeom, dUcomp[GEOMSOURCE], &tempaccdt, &tempgravitydt);
      if(accdt>tempaccdt){
        accdt=tempaccdt;
        accdti=i;
        accdtj=j;
        accdtk=k;
      }
      if(gravitydt>tempgravitydt){
        gravitydt=tempgravitydt;
        gravitydti=i;
        gravitydtj=j;
        gravitydtk=k;
      }

#if(0)// DEBUG
      if(i==-1 || i==0){
        dualfprintf(fail_file,"BANG: i=%d\n",i);
        PLOOP(pliter,pl) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl,MACP0A1(prim,i,j,k,pl));
        PLOOP(pliter,pl) dualfprintf(fail_file,"dUgeom[%d]=%21.15g dUcompgeomsource=%21.15g\n",pl,dUgeom[pl],dUcomp[GEOMSOURCE][pl]);
        if(i==-1){
          // side-by-side
          coord(i,j,k,CENT,X);
          bl_coord(X,V);
          coord(ip1mac(i),j,k,CENT,X);
          bl_coord(X,Vp1);
          dualfprintf(fail_file,"r(i=-1)=%21.15g r(i=0)=%21.15g\n",V[1],Vp1[1]);
          dualfprintf(fail_file,"gdet(i=-1)=%21.15g gdet(i=0)=%21.15g\n",GLOBALMETMACP1A0(gdet,CENT,i,j,k),GLOBALMETMACP1A0(gdet,CENT,ip1mac(i),j,k));
          DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcov(i=-1)=%21.15g gcov(i=0)=%21.15g\n",jj,kk,GLOBALMETMACP1A2(gcov,CENT,i,j,k,jj,kk),GLOBALMETMACP1A2(gcov,CENT,ip1mac(i),j,k,jj,kk));
          DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcovlast(i=-1)=%21.15g gcovlast(i=0)=%21.15g\n",jj,kk,GLOBALMETMACP1A2(gcovlast,CENT,i,j,k,jj,kk),GLOBALMETMACP1A2(gcovlast,CENT,ip1mac(i),j,k,jj,kk));
          DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcon(i=-1)=%21.15g gcon(i=0)=%21.15g\n",jj,kk,GLOBALMETMACP1A2(gcon,CENT,i,j,k,jj,kk),GLOBALMETMACP1A2(gcon,CENT,ip1mac(i),j,k,jj,kk));
          DLOOP(jj,kk) DLOOPA(ll) dualfprintf(fail_file,"%d %d %d conn(i=-1)=%21.15g conn(i=0)=%21.15g\n",jj,kk,ll,GLOBALMETMACP0A3(conn,i,j,k,jj,kk,ll),GLOBALMETMACP0A3(conn,ip1mac(i),j,k,jj,kk,ll));
        }
      }
#endif


#endif


    } // end if within source enerregion

  } // end of loop






    // GODMARK: note that in normal advance, wavendt[i] is over each CPU region and wavedt computed for each CPU and then minimized over all CPUs -- so not perfectly consistent with MPI
    // here we preserve perfect MPI domain decomposition
  mpifmin(&wavendt[1]);
  mpifmin(&wavendt[2]);
  mpifmin(&wavendt[3]);
  // single all-CPU wavedt for PERCELLDT==0 version
  wavedt_1 = MINDTSET(wavendt[1],wavendt[2],wavendt[3]); // wavendt[i] is over entire region for each i

  // minimize per-cell dt over all CPUs for PERCELLDT==1 version
  mpifmin(&wavedt_2);

  // get actual version to use
  if(PERCELLDT==0) wavedt=wavedt_1;
  else wavedt=wavedt_2;

  // single all-CPU accdt and gravitydt
  mpifmin(&accdt);
  mpifmin(&gravitydt);

  wavedtglobal=wavedt;
  sourcedtglobal=accdt;
  gravitydtglobal=gravitydt;


  // find global minimum value of wavendt over all cpus
  ndtfinal=MIN(wavedt,MIN(accdt,gravitydt));

#if(1)
  // below single line only right if 1-CPU
  SLOOPA(jj) dualfprintf(log_file,"dtij[%d]=%21.15g wavendti=%d wavendtj=%d wavendtk=%d\n",jj,wavendt[jj],wavendti[jj],wavendtj[jj],wavendtk[jj]);
  dualfprintf(log_file,"wavedt_1=%21.15g wavedt_2=%21.15g\n",wavedt_1,wavedt_2); // report this so can compare PERCELLDT==0 vs. 1 at least when starting or restarting runs
  dualfprintf(log_file,"ndtfinal=%21.15g wavedt=%21.15g accdt=%21.15g gravitydt=%21.15g\n",ndtfinal,wavedt,accdt,gravitydt); 
  dualfprintf(log_file,"accdti=%d accdtj=%d accdtk=%d :: gravitydti=%d  gravitydtj=%d  gravitydtk=%d\n",accdti,accdtj,accdtk,gravitydti,gravitydtj,gravitydtk);
#endif

  *dt = ndtfinal;


  // don't step beyond end of run
  if (t + *dt > tf) *dt = tf - t;
  
  return(0);
}


// 0.5 not good enough for pressureless collapse
// normal cour=0.8/4 works for presureless collapse for longsteps, so use 0.1 to be safe since rarely gravity conquers timestep
// but cour=0.8 and GRAVITYCOUR = 0.1 doesn't even work for longsteps!
#define GRAVITYCOUR (0.1)

/// note that dUevolve and dUgeomevolve are really dU/dt (i.e. per unit dt)
static int compute_dt_fromsource(struct of_geom *ptrgeom, struct of_state *state, FTYPE *pr, FTYPE *U, FTYPE *dUevolvedt, FTYPE *dUgeomevolvedt, FTYPE *dtij, FTYPE *gravitydt)
{
  FTYPE dUevolve[NDIM],dUgeomevolve[NDIM];
  FTYPE dUd[NDIM],dUu[NDIM];
  int jj,kk;
  FTYPE rhoprime[NDIM];
  FTYPE ag[NDIM],dtsource[NDIM];
  FTYPE rho,u,P,bsq,w,eta;
  FTYPE mydU[NDIM];
  FTYPE mydUdtgravity, rhoprimegravity, aggravity;
  FTYPE frac;
  int i,j,k,loc;
  extern void compute_dr(int i, int j, int k, FTYPE *dr);
  FTYPE dr,dphidt,phi,tempdt;
  FTYPE veleff;
  FTYPE v1max,v1min;
  FTYPE mydx[NDIM];


  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;

  // default is no limit on dt due to flux differential or source terms
  *dtij=BIG;
  // default is no limit due to time-dependence of gravity
  *gravitydt=BIG;


  // convert from dU/dt to dU
  DLOOPA(jj){
    dUevolve[UU+jj] = dUevolvedt[UU+jj]*dt;
    dUgeomevolve[UU+jj] = dUgeomevolvedt[UU+jj]*dt;
  }


  DLOOPA(jj){
    dUd[jj]=dUevolve[UU+jj]*(ptrgeom->IEOMFUNCNOSINGMAC(UU+jj)); // remove geometry
  }
  raise_vec(dUd,ptrgeom,dUu);

  // comparing time update with time value, so keep lower as conserved quantity
  mydU[TT]=dUd[TT];
  mydUdtgravity = dUgeomevolvedt[UU]*(ptrgeom->IEOMFUNCNOSINGMAC(UU)); // pure gravity piece

  mydx[TT]=dt; // so in the end dt_time <=C*dt/(dU/rhoprime)
  SLOOPA(jj){
    // treating momentum update as giving dv^i, so for getting Courant condition of moving across grid, need upper
    mydU[jj] = dUu[jj];
    mydx[jj] = dx[jj];
  }
    

  bsq = dot(state->bcon, state->bcov);
  rho=pr[RHO];
  u=pr[UU];
  P=pressure_rho0_u_simple(i,j,k,loc,rho,u);
  w = rho+u+P;
  eta = w+bsq;





  /////////////////////////
  //
  // New method for dealing with source terms
  //
  /////////////////////////

  DLOOPA(jj){
    // U[UU] is like eta but with \gamma^2 factors that are present in the momentum terms
    // account for geometry prefactor of conserved quantities and put in geometry consistent with source term
    // U[UU] \sim eta \gamma^2 when neglecting stress terms
    // GODMARK: Might want to be more careful like in utoprim_jon.c in how normalize errors
    // GODMARK: Consider REMOVERESTMASSFROMUU==2 effects
    rhoprime[jj]=MAX(fabs(eta),fabs(U[UU]));

    // update to 3-velocity v^i is approximately due to this acceleration
    // ag^j \sim (\rho\gamma^2 v^i)/(\rho\gamma^2) \sim dv^i/dt
    // below is really = acceleration * dt
    ag[jj]=SMALL+fabs(mydU[jj]/rhoprime[jj]); // acceleration.  SMALL is so doesn't identically vanish and we get nan below

    // for time, idea is to keep d\rho/\rho\lesssim cour, so dtnew = dtold *\rho/d\rho = dtold/ag[TT]
    if(jj==TT) dtsource[jj]=cour*(mydx[jj]/ag[jj]);
    // dt = dx/ag since ag is really  = dv^i = dx^i/dt
    else dtsource[jj]=cour*(mydx[jj]/ag[jj]); // characteristic time-step for increasing velocity so mass would cross a grid
      
  }




  /////////////////////////
  //
  // Self-gravity
  //
  /////////////////////////

#if(DOSELFGRAVVSR)
  // make sure metric is not varying too fast
  // just look at perturbed part of g_{tt} -- an invariant in stationary metrics, so good indicator of gauge-invariant time-dependence of metric
  // only look at loc=CENT since others should be similar to order unity -- also CENT won't diverge at r=0
  //  frac = (METMACP1A1(gcovpertlast,CENT,i,j,k,TT) - METMACP1A1(gcovpert,CENT,i,j,k,TT))/(fabs(METMACP1A1(gcovpertlast,CENT,i,j,k,TT)) + fabs(METMACP1A1(gcovpert,CENT,i,j,k,TT)));
  // above frac has no dt, so no measure of what dt should be

  // \Gamma^t_{tt} measures dg_{tt}/dt
  // \Gamma^t_{tt} \approx g_{tt},t g^{tt}  so that g_{tt},t = \Gamma^t_{tt}/g^{tt}
  // now form same construct as with (dU/dt)/U as above
  // g^{tt} can't go to 0 unless in bad coordinate system (BL-coords on horizon)
  // frac is approximately g_{tt},t/g^{tt} \aprox 1/dt
  // GODMARK: below might be problem in non-relativistic case



  // get \Gamma_{ttt} = dphi/dt ~ 1/2 g_{tt,t}
  dphidt=0.0;
  DLOOPA(jj) dphidt += GLOBALMETMACP0A3(conn,i,j,k,jj,TT,TT)*(ptrgeom->gcov[GIND(jj,TT)]);
  dphidt = fabs(dphidt); // don't care about sign, just magnitude

  // get \phi ~ -(g_{tt} +1)/2
  // phi by itself has no meaning, but as a reference for changes in phi in time it does
  phi = -(1.0+ptrgeom->gcov[GIND(TT,TT)])*0.5;
  phi = fabs(phi); // sign not important

#if(0)

  // treat dt ~ \phi / (d\phi/dt)
  //  frac = fabs(GLOBALMETMACP0A3(conn,i,j,k,TT,TT,TT)/((1+GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT))*(1+GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT))));
  frac = fabs(dphidt/phi);
  *gravitydt = cour*(GRAVITYCOUR/frac);
  // this dt keeps frac~cour
  //  *gravitydt = cour*(GRAVITYCOUR/frac); // GRAVITYCOUR is additional courant factor on gravitational term


  // treat d\phi/dt as v^2/dt
  compute_dr( i,  j,  k, &dr);
  //dgttdt = fabs(GLOBALMETMACP0A3(conn,i,j,k,TT,TT,TT)/(GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT)));
  //tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dgttdt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dphidt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term

  if(tempdt<*gravitydt) *gravitydt=tempdt;
#elif(0)

  frac = fabs(ptrgeom->gcon[GIND(TT,TT)]*dphidt);
  *gravitydt = cour*(GRAVITYCOUR/frac);




#elif(1)

  ///////////////////////////////
  //
  // 3 methods to limit dt
  //
  //////////////////////////////

  // LOCAL DPHI/DT LIMIT
  // treat dt ~ \phi / (d\phi/dt)
  //  frac = fabs(GLOBALMETMACP0A3(conn,i,j,k,TT,TT,TT)/((1+GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT))*(1+GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT))));
  frac = fabs(dphidt/phi);
  *gravitydt = cour*(GRAVITYCOUR/frac);
  // this dt keeps frac~cour
  //  *gravitydt = cour*(GRAVITYCOUR/frac); // GRAVITYCOUR is additional courant factor on gravitational term



  compute_dr( i,  j,  k, &dr);


  // LOCAL DU/DT LIMIT BUT TREAT AS VELOCITY LIMITED BY SOL

#if(REMOVERESTMASSFROMUU==2)
  rhoprimegravity=fabs(U[UU]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#else
  // remove rest-mass
  rhoprimegravity=fabs(U[UU]+U[RHO]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#endif
  aggravity=SMALL+fabs(mydUdtgravity/rhoprimegravity);
  veleff = aggravity*mydx[1];

  // get speed of light in 1-direction (dx^1/dt)
  sol(pr,state,1,ptrgeom,&v1max,&v1min);
  // limit "velocity" to speed of light
  if(veleff>fabs(v1min)) veleff=fabs(v1min);
  if(veleff>fabs(v1max)) veleff=fabs(v1max);

  tempdt=GRAVITYCOUR*cour*mydx[1]/veleff;
  if(tempdt<*gravitydt) *gravitydt=tempdt;


  // LOCAL LIMIT ON DPHI/DT WHERE PHI~V^2
  //dgttdt = fabs(GLOBALMETMACP0A3(conn,i,j,k,TT,TT,TT)/(GLOBALMETMACP1A2(gcon,CENT,i,j,k,TT,TT)));
  //tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dgttdt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dphidt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  if(tempdt<*gravitydt) *gravitydt=tempdt;


#endif




#endif // end if DOSELFGRAVVSR==1










  /////////////////////////
  //
  // Finally store source term's version of limited dt to be used later
  //
  /////////////////////////

  //  if(ptrgeom->i==30 && ptrgeom->j==31){
  //  DLOOPA(jj) dualfprintf(fail_file,"dtsource[%d]=%21.15g : %21.15g : %21.15g %21.15g : %21.15g\n",jj,dtsource[jj],cour,mydx[jj],ag[jj],dUevolve[UU+jj]);
  //} 

  //  dualfprintf(fail_file,"i=%d mydUdtgravity=%21.15g rhoprimegravity=%21.15g rhoprime[TT]=%21.15g mydU[TT]=%21.15g\n",ptrgeom->i,mydUdtgravity,rhoprimegravity,rhoprime[TT],mydU[TT]);

  // always do time-component
  // accounts for thermal changes if cooling function or geometry changes if metric changing
  if (dtsource[TT] < *dtij) *dtij = dtsource[TT];

#if(N1>1)
  if (dtsource[RR] < *dtij) *dtij = dtsource[RR];
#endif
#if(N2>1)
  if (dtsource[TH] < *dtij) *dtij = dtsource[TH];
#endif
#if(N3>1)
  if (dtsource[PH] < *dtij) *dtij = dtsource[PH];
#endif


  return(0);

}



/// compute dt from full coordinate acceleration
static int dUtodt(struct of_geom *ptrgeom, struct of_state *qaddr, FTYPE *pr, FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *dUgeomgravity, FTYPE *accdt, FTYPE *gravitydt)
{
  int pl,pliter;
  FTYPE dUtotal[NPR];
  FTYPE Ugeomfree[NPR],U[NPR];



  PLOOP(pliter,pl) {
    // while each piece may be "large", when summed if small then final change to conserved quantity is small, so that's all that's relevant
    dUtotal[pl]=dUriemann[pl]+dUgeom[pl];
    // GODMARK:
    //    dUtotal[pl]=fabs(dUriemann[pl])+fabs(dUgeom[pl]);
  }

  // conserved quantity without geometry
  MYFUN(primtoU(UEVOLVE, pr, qaddr, ptrgeom, U, NULL),"step_ch.c:advance()", "primtoU()", 1);
  PLOOP(pliter,pl) Ugeomfree[pl] = U[pl]*(ptrgeom->IEOMFUNCNOSINGMAC(pl));


  compute_dt_fromsource(ptrgeom, qaddr, pr, Ugeomfree, dUtotal, dUgeomgravity, accdt, gravitydt);
    


  return(0);

}

