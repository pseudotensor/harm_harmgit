
/*! \file flux.c
    \brief Routines for computing fluxes except those special routines directly related to FLUXB==FLUXCTSTAG or FLUXCTTOTH

    // OPTMARK: Should redo flux's so that fluxes are accessed by MAC(F1,j,k,i) MAC(F2,k,i,j) MAC(F3,i,j,k) for faster differencing in advance.c
    // Maybe not important

*/


#include "decs.h"




static void leftrightcompute(int i, int j, int k, int dir, int is, int ie, int js, int je, int ks, int ke, int *computewithleft, int *computewithright);





/// get whether should compute using left/right states
static void leftrightcompute(int i, int j, int k, int dir, int is, int ie, int js, int je, int ks, int ke, int *computewithleft, int *computewithright)
{

#if(MERGEDC2EA2CMETHOD)
  // if doing merged method, then expanded "flux" calculation to get state in full (i.e. left,cent,right) in the "-1" and "N" cells
  if(dir==1 && (i==is) || dir==2 && (j==js) || dir==3 && (k==ks)){ *computewithleft=0; *computewithright=1; }
  else if(dir==1 && (i==ie) || dir==2 && (j==je) || dir==3 && (k==ke)){ *computewithleft=1; *computewithright=0; }
  else{ *computewithleft=1; *computewithright=1; }
#else
  // normal routine when always get final dissipative flux
  *computewithleft=1; *computewithright=1;
#endif


}





/// see fluxcompute.c for non-computer science, real physics calculations of flux
int fluxcalc(int stage,
             int initialstep, int finalstep,
             FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
             FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
             FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
             FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
             FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
             FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
             FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
             FTYPE *CUf,
             FTYPE *CUnew,
             SFTYPE fluxdt,
             SFTYPE fluxtime,
             FTYPE *ndt1,
             FTYPE *ndt2,
             FTYPE *ndt3
             )
{
  int fluxcalc_flux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, SFTYPE time, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop);
  void fix_flux(FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]) ;
  FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (**ptrfluxvec)[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE *ndtvec[NDIM];
  int Nvec[NDIM];
  int flux_point2avg(int stage, int whichmaorem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  void preinterp_flux_point2avg(void);
  int i,j,k;
  int pl,pliter;
  int dir;
  int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int cleanup_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int zero_out_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvecem[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);

  // face2faceloop[dir=facedir] and face2cornloop[EMFdir=edgedir][EMFodir1][EMFodir2]
  // face2centloop separately used during separate part of advance()
  struct of_loop cent2faceloop[NDIM],face2cornloop[NDIM][NDIM][NDIM];


  /////////////////////////
  //
  // SETUP dimensionality
  //
  /////////////////////////

  fluxvec[1]=F1;
  fluxvec[2]=F2;
  fluxvec[3]=F3;

  fluxvecEM[1]=GLOBALPOINT(F1EM); // more temporary than F1,F2,F3, so don't pass to here and just keep global
  fluxvecEM[2]=GLOBALPOINT(F2EM);
  fluxvecEM[3]=GLOBALPOINT(F3EM);

  dqvec[1]=GLOBALPOINT(dq1);
  dqvec[2]=GLOBALPOINT(dq2);
  dqvec[3]=GLOBALPOINT(dq3);

  ndtvec[1]=ndt1;
  ndtvec[2]=ndt2;
  ndtvec[3]=ndt3;

  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  if(splitmaem) ptrfluxvec=fluxvecEM;
  else ptrfluxvec=fluxvec;


  ///////////////////////////////////////////////
  //
  // Zero-out fluxes so can easily update full region in advance.c without special conditions
  //
  ///////////////////////////////////////////////
  //  zero_out_fluxes(Nvec,ptrfluxvec); // NOT FOR NOW



  // Below no longer needed since always zero-out all fluxes (reverted)
  // zero-out EMFs if evolving/tracking vector potential so boundary regions appear simple in DUMPS when showing boundary cells
  if(TRACKVPOT && FULLOUTPUT!=0){
    zero_out_emf_fluxes(Nvec,ptrfluxvec);
  }

  ///////////////////////////////////////////////
  //
  // some pre-interplation flux averaging setups
  //
  ///////////////////////////////////////////////
  preinterp_flux_point2avg();



  ///////////////////////////////////////////////
  //
  // Compute normal flux at face
  // Involves interpolating CENT -> FACE1,2,3
  // Final p_l p_r results of interpolation are stored in pl_ct pr_ct if needed by SPLITNPR or FLUXB==FLUXCTSTAG
  // Wavespeeds may also be stored globally if STOREWAVESPEEDS>0
  // In all cases wavespeed constraint on timestep is set here
  //
  // assume fluxvec is MA only if splitmaem==1
  //
  ///////////////////////////////////////////////
  
  //  CUf[2]=current dt to be used on flux
  fluxcalc_flux(stage, pr, pstag, pl_ct, pr_ct, Nvec, dqvec, fluxvec, fluxvecEM, CUf[2], fluxtime, ndtvec, cent2faceloop);







#if(0)
  // DEBUG:
  if(Nvec[1]>1) FULLLOOP PLOOP(pliter,pl) MACP1A1(fluxvec,1,i,j,k,pl)+=MACP1A1(fluxvecEM,1,i,j,k,pl);
  if(Nvec[2]>1) FULLLOOP PLOOP(pliter,pl) MACP1A1(fluxvec,2,i,j,k,pl)+=MACP1A1(fluxvecEM,2,i,j,k,pl);
  if(Nvec[3]>1) FULLLOOP PLOOP(pliter,pl) MACP1A1(fluxvec,3,i,j,k,pl)+=MACP1A1(fluxvecEM,3,i,j,k,pl);
#endif



  //////////////////////////////
  //
  // FIXFLUX
  //
  /////////////////////////////

  if(FIXUPFLUX && splitmaem==0){ // for now can't use fix_flux with splitmaem since it resets field part to 0 that I'm using as diagonal gas pressure term
    fix_flux(pr, fluxvec[1], fluxvec[2], fluxvec[3]);
    if(splitmaem) fix_flux(pr, fluxvecEM[1], fluxvecEM[2], fluxvecEM[3]);
#if(PRODUCTION==0)
    trifprintf( "x");
#endif
  }


  //////////////////////////////
  //
  // ADJUST FLUXES
  //
  /////////////////////////////

  if(ADJUSTFLUX){
    extern void adjust_flux(SFTYPE fluxtime,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);

    adjust_flux(fluxtime, pr, fluxvec[1], fluxvec[2], fluxvec[3]);
#if(PRODUCTION==0)
    trifprintf( "q");
#endif
  }




  //////////////////////////////
  //
  // FLUXCTSTAG -- overwrite field fluxes (emf's) with correct CORN1,2,3 points values
  //
  // assumes fluxcalc_flux() above called fluxcalc_standard_4fluxctstag() so pl_ct and pr_ct are set with FACE1,2,3 values of all quantities (including field along dir that comes from pstagscratch[] in this case)
  //
  /////////////////////////////
  if(FLUXB==FLUXCTSTAG){
#if(STOREWAVESPEEDS==0 || USESTOREDSPEEDSFORFLUX==0)
    dualfprintf(fail_file,"STOREWAVESPEEDS,USESTOREDSPEEDSFORFLUX must be >0 when FLUXB==FLUXCTSTAG\n");
    // must store because vchar() cannot be computed at CORN since only have partial velocities and partial fields and no densities
    // really only STOREWAVESPEEDS must be >0, but for now assume both
    myexit(175106);
#endif


    MYFUN(fluxcalc_fluxctstag(stage, initialstep, finalstep, pr, pstag, pl_ct, pr_ct, GLOBALPOINT(pvbcorninterp), GLOBALPOINT(wspeed), GLOBALPOINT(prc), GLOBALPOINT(pleft), GLOBALPOINT(pright), GLOBALPOINT(fluxstatecent), GLOBALPOINT(fluxstate), GLOBALPOINT(geomcornglobal), Nvec, dqvec, ptrfluxvec, vpot, CUf, CUnew, fluxdt, fluxtime, cent2faceloop, face2cornloop),"flux.c:fluxcalc()", "fluxcalc_fluxctstag", 0);

  }// end if staggered method where can update A_i directly



#if(PRODUCTION==0)
  trifprintf( "c");
#endif



#if(0)
  // DEBUG:
  FULLLOOP{
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
        PLOOP(pliter,pl){
          if(pl>=B1 && pl<=B2){
            if(!isfinite(MACP1A1(fluxvec,dir,i,j,k,pl))){
              dualfprintf(fail_file,"GOD FLUXA: dir=%d pl=%d i=%d j=%d k=%d\n",dir,pl,i,j,k,MACP1A1(fluxvec,dir,i,j,k,pl));
            }
          }
        }
      }
    }
  }
#endif





  //////////////////////////////
  //
  // convert point FLUXES to surface-averaged fluxes (for field, if FLUXB==FLUXCTSTAG, then should do line integral of emf ENOMARK)
  //
  /////////////////////////////
  if((interporder[avgscheme[1]]>3) ||  (interporder[avgscheme[2]]>3) ||  (interporder[avgscheme[3]]>3)){
    if(splitmaem){
      // assume fluxvec is MA only if splitmaem==1
      //flux_point2avg(stage, ISMAONLY, pr, Nvec, fluxvec,NULL);
      // below version forces summation of MA+EM for stencil
      flux_point2avg(stage, ISMAONLY, pr, Nvec, fluxvec,fluxvecEM);
      flux_point2avg(stage, ISEMONLY, pr, Nvec, fluxvecEM, NULL);
    }
    else{
      flux_point2avg(stage, ISMAANDEM, pr, Nvec, fluxvec,NULL);
    }
  }





#if(0)
  // DEBUG: // also helped after flux_point2avg() but not before! -- so higher order code is problem or fed bad data
  // bound_flux(-1,F1,F2,F3);
#endif



  //////////////////////////////
  //
  // FLUXCT : must come after modifying fluxes so that divb=0 is preserved to machine error
  // i.e. flux_ct modifies fluxes in special way just before updating field so that divb=0 is conserved
  // Method is only 2nd order at best since assumes volume average is same as surface average
  //
  /////////////////////////////
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)||(FLUXB==FLUXCD)){

    MYFUN(flux_ct(stage, initialstep, finalstep, pr, GLOBALPOINT(emf), GLOBALPOINT(vconemf), GLOBALPOINT(dq1), GLOBALPOINT(dq2), GLOBALPOINT(dq3), ptrfluxvec[1], ptrfluxvec[2], ptrfluxvec[3], vpot, Nvec, CUf, CUnew, fluxdt, fluxtime),"step_ch.c:advance()", "flux_ct",1);

  }


  //////////////////////
  //
  // sum up MA+EM
  // if splitmaem==1, then MA and EM split up to this point
  //
  ///////////////////////

  if(splitmaem) fluxsum(Nvec, fluxvec,fluxvecEM);
  //////////////////
  //
  // now fluxvec has MA+EM and fluxvecEM is no longer needed
  //
  //////////////////



  //////////////////////////////
  //
  // Merged c2e-c2a method
  //
  /////////////////////////////
  if(MERGEDC2EA2CMETHOD){
    mergedc2ea2cmethod_compute(Nvec,fluxvec);
  }





  //////////////////////////////
  //
  // Clean-up fluxes so zero outside well-defined computational box
  // Turns out to be easier to clean-up rather than ensure computed in right places
  //
  /////////////////////////////
  cleanup_fluxes(Nvec,ptrfluxvec);


  //////////////////////////////
  //
  // DEBUG:
  //
  /////////////////////////////
#if(FLUXDUMP==1)
  // this accounts for final flux
  FULLLOOP{
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
        PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=MACP1A1(fluxvec,dir,i,j,k,pl);
      }
      else{
        PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=0.0;
      }
    }
  }    
#endif



  // DEBUG: // helped!
  //  bound_prim(-1,F1,NULL,NULL,1);
  //  bound_prim(-1,F2,NULL,NULL,1);

  // DEBUG: // also helped!
  //bound_flux(-1,F1,F2,F3);


  return(0);
  
  
}



/// ensure fluxes only exist on well-defined computational box
/// makes advance.c, GRIDSECTIONING, and adaptive time-stepping easier to code for
/// Note that the branch prediction penalty for these conditions inside the LOOP is low since just zero out nearby memory elements
int cleanup_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{


  
  // for nowait to matter, has to be on inner loop but outer (dir) loop must be inside the parallel construct since a parallel construct has an impassible barrier.  If one put the parallel region inside the outer (dir) loop, the nowait would be useless.  Then interior things that depend upon "dir" must be made private, including "dir" itself
#pragma omp parallel 
  {
    int i,j,k;
    int dir;
    int pl,pliter;
    int is,ie,js,je,ks,ke;
    int B1is,B1ie,B1js,B1je,B1ks,B1ke;
    int B2is,B2ie,B2js,B2je,B2ks,B2ke;
    int B3is,B3ie,B3js,B3je,B3ks,B3ke;
    int loop[NUMFLUXLOOPNUMBERS];
    //  int odir1,odir2;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL; // doesn't depend upon dir, so can be outside DIMENLOOP.  However, blockijk must be private, so go ahead and put inside parallel region

    // pick reference loop as where evolve conserved quantities
    // must shift since indices not used as arguments to LOOP
    loop[FIS]=Uconsevolveloop[FIS]+SHIFTX1DN;
    loop[FIE]=Uconsevolveloop[FIE]+SHIFTX1UP;
    loop[FJS]=Uconsevolveloop[FJS]+SHIFTX2DN;
    loop[FJE]=Uconsevolveloop[FJE]+SHIFTX2UP;
    loop[FKS]=Uconsevolveloop[FKS]+SHIFTX3DN;
    loop[FKE]=Uconsevolveloop[FKE]+SHIFTX3UP;



    DIMENLOOP(dir){
      if(Nvec[dir]>1){

        //      get_odirs(dir,&odir1,&odir2);

        is=loop[FIS];
        ie=loop[FIE]+(dir==1)*SHIFT1;
        js=loop[FJS];
        je=loop[FJE]+(dir==2)*SHIFT2;
        ks=loop[FKS];
        ke=loop[FKE]+(dir==3)*SHIFT3;


        if(FLUXB==FLUXCTSTAG){
          // E1: Needed: i==0..N1-1 j=0..N2   k=0..N3
          // E2: Needed: i==0..N1   j=0..N2-1 k=0..N3
          // E3: Needed: i==0..N1   j=0..N2   k=0..N3-1
          // Note Fi[Bi]=0 always and should already be set everywhere

          // B1
          B1is=loop[FIS];
          B1ie=loop[FIE]+(dir==1 || dir==2 || dir==3)*SHIFT1; // for F1[B1],F2[B1]=E3,F3[B1]=E2
          B1js=loop[FJS];
          B1je=loop[FJE]+(dir==1 || dir==2)*SHIFT2; // for F1[B1],F2[B1]=E3
          B1ks=loop[FKS];
          B1ke=loop[FKE]+(dir==1 || dir==3)*SHIFT3; // for F1[B1],F3[B1]=E2

          // B2
          B2is=loop[FIS];
          B2ie=loop[FIE]+(dir==1 || dir==2)*SHIFT1; // for F1[B2]=E3,F2[B2],F3[B2]=E1
          B2js=loop[FJS];
          B2je=loop[FJE]+(dir==1 || dir==2 || dir==3)*SHIFT2; // for F1[B2]=E3,F2[B2],F3[B2]=E1
          B2ks=loop[FKS];
          B2ke=loop[FKE]+(dir==2 || dir==3)*SHIFT3; // for F2[B2],F3[B2]=E1

          // B3
          B3is=loop[FIS];
          B3ie=loop[FIE]+(dir==1 || dir==3)*SHIFT1; // for F1[B3]=E2,F3[B3]
          B3js=loop[FJS];
          B3je=loop[FJE]+(dir==2 || dir==3)*SHIFT2; // for F1[B3]=E2,F2[B3]=E1,F3[B3]
          B3ks=loop[FKS];
          B3ke=loop[FKE]+(dir==1 || dir==2 || dir==3)*SHIFT3; // for F1[B3]=E2,F2[B3]=E1,F3[B3]
        }
        else{
          B1is=B2is=B3is=is;
          B1ie=B2ie=B3ie=ie;
          B1js=B2js=B3js=js;
          B1je=B2je=B3je=je;
          B1ks=B2ks=B3ks=ks;
          B1ke=B2ke=B3ke=ke;
        }


        // dualfprintf(fail_file,"dir=%d Clean1: is=%d ie=%d js=%d je=%d ks=%d ke=%d\n",dir,is,ie,js,je,ks,ke);
        // dualfprintf(fail_file,"CleanB1: B1is=%d B1ie=%d B1js=%d B1je=%d B1ks=%d B1ke=%d\n",B1is,B1ie,B1js,B1je,B1ks,B1ke);
        // dualfprintf(fail_file,"CleanB2: B2is=%d B2ie=%d B2js=%d B2je=%d B2ks=%d B2ke=%d\n",B2is,B2ie,B2js,B2je,B2ks,B2ke);
        // dualfprintf(fail_file,"CleanB3: B3is=%d B3ie=%d B3js=%d B3je=%d B3ks=%d B3ke=%d\n",B3is,B3ie,B3js,B3je,B3ks,B3ke);


        ////      COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);

          // below means we are not within the computational block
          if(! (i>=is && i<=ie && j>=js  && j<=je && k>=ks && k<=ke) ){
            PLOOPNOB1(pl) MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
            PLOOPNOB2(pl) MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
          }
          if(! (i>=B1is && i<=B1ie && j>=B1js && j<=B1je && k>=B1ks && k<=B1ke) ){
            pl=B1; MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
          }
          if(! (i>=B2is && i<=B2ie && j>=B2js && j<=B2je && k>=B2ks && k<=B2ke) ){
            pl=B2; MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
          }
          if(! (i>=B3is && i<=B3ie && j>=B3js && j<=B3je && k>=B3ks && k<=B3ke) ){
            pl=B3; MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
          }
        }// end 3D loop
      }// end if doing this dimen
    }// end over dimens
  }// end parallel region (no need for added barrier since parallel region has impassible barrier)




  return(0);

}



/// zero-out fluxes in COMPFULLLOOP so advance.c updates all possible quantities and does so with simple non-conditional code
/// it's ok to zero-out beyond standard+-1 box a bit.  Just ensure zero out enough
int zero_out_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    int dir;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

    DIMENLOOP(dir){
      if(Nvec[dir]>1){
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);
          ////COMPFULLLOOP
          PLOOP(pliter,pl) MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
        }// end over 3D loop
      }// if dir
    }// over dirs
  }// end parallel region (no need for added barrier since parallel region has impassible barrier)



  return(0);
}

/// zero-out fluxes that correspond to EMFs if evolving/tracking A_i since want boundary values to be plottable in SM or whatever
/// This is  necessary since also use fluxes for other purposes in-between steps
/// GODMARK: If in BCs have partial EMF and not full EMF for updating A_i, then may also cause A_i to appear funny
int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    int dir;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

    DIMENLOOP(dir){
      if(Nvec[dir]>1){
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);
          ////COMPFULLLOOP
          PLOOPBONLY(pl) MACP1A1(fluxvec,dir,i,j,k,pl)=0.0;
        }// end over 3D loop
      }// if dir
    }// over dirs
  }// end parallel region (no need for added barrier since parallel region has impassible barrier)


  return(0);
}




// fill EM version if necessary when having new CT EMFs
int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    int dir;

    if(splitmaem){

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
      DIMENLOOP(dir){
        if(Nvec[dir]>1){
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
          OPENMP3DLOOPBLOCK{
            OPENMP3DLOOPBLOCK2IJK(i,j,k);
            //// COMPFULLLOOP
            PLOOPBONLY(pl) MACP1A1(fluxvecEM,dir,i,j,k,pl)=MACP1A1(fluxvec,dir,i,j,k,pl);
          }// end over 3D loop
        }// if dir
      }// over dirs
    }// end if splitmaem
  }// end parallel region


  return(0);
}

/// sums MA and EM fluxes
int fluxsum_old(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    int dir;
    
    if(splitmaem){

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
      DIMENLOOP(dir){
        if(Nvec[dir]>1){
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
          OPENMP3DLOOPBLOCK{
            OPENMP3DLOOPBLOCK2IJK(i,j,k);
            ////COMPFULLLOOP
            PLOOP(pliter,pl) MACP1A1(fluxvec,dir,i,j,k,pl)+=MACP1A1(fluxvecEM,dir,i,j,k,pl);
          }// end over 3D loop
        }// end if dir
      }// end over dirs
    }// end if splitmaem
  }// end parallel region


  return(0);
}


/// sums MA (with FLUXSPLITMA(dir) containing flux[UU+dir] component) and EM fluxes
int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    int dir;
    
    if(splitmaem){

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

      DIMENLOOP(dir){
        if(Nvec[dir]>1){
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
          OPENMP3DLOOPBLOCK{
            OPENMP3DLOOPBLOCK2IJK(i,j,k);
            //// COMPFULLLOOP{
#if(SPLITPRESSURETERMINFLUXMA)
            // add diagonal pressure term back to normal FLUX term
            MACP1A1(fluxvec,dir,i,j,k,UU+dir)+=MACP1A1(fluxvec,dir,i,j,k,FLUXSPLITPMA(dir));
            // reset temporary storage
            MACP1A1(fluxvec,dir,i,j,k,FLUXSPLITPMA(dir))=0.0;
#endif
#if(SPLITPRESSURETERMINFLUXEM)
            // add diagonal pressure term back to normal FLUX term
            MACP1A1(fluxvec,dir,i,j,k,UU+dir)+=MACP1A1(fluxvecEM,dir,i,j,k,FLUXSPLITPEM(dir));
            // reset temporary storage
            MACP1A1(fluxvecEM,dir,i,j,k,FLUXSPLITPEM(dir))=0.0;
#endif
            // now do normal addition
            PLOOP(pliter,pl) MACP1A1(fluxvec,dir,i,j,k,pl)+=MACP1A1(fluxvecEM,dir,i,j,k,pl);
          }// over 3D LOOP
        }// if dir
      }// over dirs
    }// end if splitmaem
  }// end parallel region


  return(0);
}






/// wrapper for CENT_to_FACE1,2,3 used to compute flux at face
int fluxcalc_flux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, SFTYPE time, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop)
{
  int fluxcalc_flux_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata );
  int dir;
  int idel, jdel, kdel, face;
  int is, ie, js, je, ks, ke;
  int didassigngetstatecentdata;


  ///////////////////////////////////////////////
  //
  // LOOP OVER DIMENSIONS interpolating within each dimension separately for that flux -- assumes flux determined by 1-D Riemann problem
  //
  ///////////////////////////////////////////////
  didassigngetstatecentdata=0;
  
  // OPENMPOPTMARK: Can actually do each direction independently, so able to parallelize this! (But would required, e.g., an array per dir for (e.g.) pleft/pright,wspeedtemp (global variable used inside for temp space for rescale().  *ndt is done per-dir, so that's ok (it would be fine already since using "critical" that applies to *all* threads not just one team))).
  // Note also ok to do each dir separately with FLUXB==FLUXCTSTAG.
  // OPENMPOPTMARK: Otherwise, looks doable.

  DIMENLOOP(dir){

    // set dimension having no influence on dt by default
    *(ndtvec[dir])=BIG;
    
    // don't skip dimension if doesn't exist, will be taken care of inside fluxcalc_flux_1d()


    // get loop details
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    face = fluxloop[dir][FFACE];

    //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
    is=fluxloop[dir][FIS];
    ie=fluxloop[dir][FIE];
    js=fluxloop[dir][FJS];
    je=fluxloop[dir][FJE];
    ks=fluxloop[dir][FKS];
    ke=fluxloop[dir][FKE];


    MYFUN(fluxcalc_flux_1d(stage, pr, pstag, pl_ct, pr_ct, dir, time, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dqvec[dir], fluxvec[dir], fluxvecEM[dir], CUf, ndtvec[dir], &cent2faceloop[dir], &didassigngetstatecentdata),"flux.c:fluxcalc()", "fluxcalc_flux_1d", dir);

#if(PRODUCTION==0)
    trifprintf("%d",dir);
#endif
  }// end DIMENLOOP(dir)






  ///////////////////
  //
  // set per-cell minimum for dt
  //
  ///////////////////
  if(PERCELLDT){
    FTYPE ndtveclocal[NDIM];

    {  
      int dimen;
      // set dimension having no influence on dt by default
      DIMENLOOP(dimen){
        *(ndtvec[dimen])=BIG;
        ndtveclocal[dimen]=BIG;
      }
    }

    // whether dimension is relevant
    int doingdimen[NDIM];
    doingdimen[1]=N1NOT1;
    doingdimen[2]=N2NOT1;
    doingdimen[3]=N3NOT1;

    int dimenorig=1; // choose one dimension to stick things into

#pragma omp parallel
    {
      int i,j,k;
      FTYPE wavedt;
      int dimen;

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(WITHINACTIVESECTIONEXPAND1IS,WITHINACTIVESECTIONEXPAND1IE,WITHINACTIVESECTIONEXPAND1JS,WITHINACTIVESECTIONEXPAND1JE,WITHINACTIVESECTIONEXPAND1KS,WITHINACTIVESECTIONEXPAND1KE);
    

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);
 
        // get dt for each dimension at each grid point -- but only if dimension is relevant (otherwise stays as BIG and doesn't affect wavedt)
        DIMENLOOP(dimen) if(doingdimen[dimen]) ndtveclocal[dimen]=GLOBALMACP0A1(dtijk,i,j,k,dimen);

        if(i<-1 || j<-1 || k<-1 || i>N1 || j>N2 || k>N3 || i<0 && j<0 || i<0 && k<0 || j<0 && k<0 || i>N1-1 && j>N2-1 || i>N1-1 && k>N3-1 || j>N2-1 && k>N3-1 || ndtveclocal[1] <0.0 || ndtveclocal[2] <0.0 || ndtveclocal[3] <0.0) continue; // avoid too deep into boundary or if never set. // SUPERGODMARK:  KINDA HACK, can improve.

        // set local per-cell dt
        // sum of inverses is proper for unsplit scheme based upon split interpolations/fluxes.
        wavedt = MINDTSET(ndtveclocal[1],ndtveclocal[2],ndtveclocal[3]);
 
        // use dimen=1 to store result
        dimen=dimenorig;
#pragma omp critical
        {
          if (wavedt < *(ndtvec[dimen]) ){
            *ndtvec[dimen] = wavedt;
            // below are global so can report when other dt's are reported in advance.c
            waveglobaldti[dimen]=i;
            waveglobaldtj[dimen]=j;
            waveglobaldtk[dimen]=k;
          }
        }// end critical region
      }//end over 3dloopblock
    }// end over parallel region

    // store result in all of ndt1,ndt2,ndt3 so have result no matter how many dimensions working on, and just use correctly later.
    {
      int dimen;
      DIMENLOOP(dimen){
        if(doingdimen[dimen]){
          *ndtvec[dimen]=*ndtvec[dimenorig];
          waveglobaldti[dimen]=waveglobaldti[dimenorig];
          waveglobaldtj[dimen]=waveglobaldtj[dimenorig];
          waveglobaldtk[dimen]=waveglobaldtk[dimenorig];
        }
      }
    }


  }// end over doing PERCELLDT



  return(0);

}



/// wrapper for different standard 1-D flux calculators
/// 1-D interpolate and get flux for that direction (assumes purely 1-D Riemann problem)
int fluxcalc_flux_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata )
{
  int fluxcalc_standard(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);
  int fluxcalc_standard_4fluxctstag(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);

  //  int fluxcalc_fluxspliteno(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE *ndt);
  int Nvec[NDIM];
  int i,j,k,pl,pliter;
  int odir1,odir2;


  
  if(DOENOFLUX!=ENOFLUXSPLIT){

    if(FLUXB==FLUXCTSTAG){
      MYFUN(fluxcalc_standard_4fluxctstag(stage,pr,pstag,pl_ct, pr_ct, dir,time, is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard_4fluxctstag()", 1);
    }
    else{
      // use older code that doesn't store into pl_ct and pr_ct since not needed and then waste of memory
      MYFUN(fluxcalc_standard(stage,pr,pstag,pl_ct, pr_ct, dir,time,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard()", 1);
    }
  }
  else{
    //MYFUN(fluxcalc_fluxspliteno(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt),"flux.c:fluxcalc_flux_1d()", "fluxcalc_fluxspliteno()", 1);
  }



  ////////////////////////////////////
  //
  // ensure flux of field set correctly to 0 in 1D
  //
  // if 1D in z and x, then no Ey so no F3[B1]
  // if 1D in z and y, then no Ex so no F3[B2]
  //
  // Stupid: This is automatically done since if no z-dir, then F3 not even used
  //
  ////////////////////////////////////

  /*
    Nvec[1]=N1;
    Nvec[2]=N2;
    Nvec[3]=N3;
  

    odir1=dir%3+1;
    odir2=(dir+1)%3+1;
    if(splitmaem==0){
    if(Nvec[dir]==1 && Nvec[odir1]==1) COMPFULLLOOP MACP0A1(F,i,j,k,B1-1+odir1)=0.0;
    if(Nvec[dir]==1 && Nvec[odir2]==1) COMPFULLLOOP MACP0A1(F,i,j,k,B1-1+odir2)=0.0;
    COMPFULLLOOP MACP0A1(F,i,j,k,B1-1+dir)=0.0; // flux along field is always 0 due to antisymmetry of Faraday/Maxwell
    }
    else{
    // only need to change FEM since above F doesn't contain this EMF in this case
    if(Nvec[dir]==1 && Nvec[odir1]==1) COMPFULLLOOP MACP0A1(FEM,i,j,k,B1-1+odir1)=0.0;
    if(Nvec[dir]==1 && Nvec[odir2]==1) COMPFULLLOOP MACP0A1(FEM,i,j,k,B1-1+odir2)=0.0;
    COMPFULLLOOP MACP0A1(FEM,i,j,k,B1-1+dir)=0.0; // flux along field is always 0 due to antisymmetry of Faraday/Maxwell
    // don't reset field part of MA flux since using that for other purposes (diagonal pressure term for now)
    }
  */



  return(0);
  
}

/// determine whether interpolating real primitives or substep
void set_normal_realisinterp(int *realisinterp)
{


  // real means normal primitive list = {rho,u,v1,v2,v3,B1,B2,B3} with v^i as WHICHVEL velocity
  // check:

  if(npr2interpstart==0 && npr2interpend==7 &&
     npr2interplist[RHO]==RHO &&
     npr2interplist[UU]==UU &&
     npr2interplist[U1]==U1 &&
     npr2interplist[U2]==U2 &&
     npr2interplist[U3]==U3 &&
     npr2interplist[B1]==B1 &&
     npr2interplist[B2]==B2 &&
     npr2interplist[B3]==B3
     ){

    // 1 here stands for NPR2INTERP type of loop
    // only do if some type of interpolation to be done
    *realisinterp=(RESCALEINTERP==0 || VARTOINTERP==PRIMTOINTERP);
  }
  else{
    *realisinterp=0; // must be forced to be zero
  }
  


}




/// wrapper for rescale()
void rescale_calc_full(int dir,FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR])
{

  ////COMPFULLLOOP{
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP // generally requires full information
  {
    int i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

    // generally ptr's are different inside parallel block
    ptrgeom=&geomdontuse;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
      // get geometry for center pre-interpolated values
      get_geometry(i, j, k, CENT, ptrgeom); 
      // assume no need for a guess to p2interp to get pr (consistent with no unrescale being done after interpolation)
      if(npr2interpstart<=npr2interpend) rescale(DORESCALE,dir,MAC(pr,i,j,k),ptrgeom,MAC(p2interp,i,j,k));
    }// end COMPFULLLOOP
  }// end parallel region


}





/// original flux calculator that gets F in "dir".  At end global pleft,pright,dq also set and if STOREWAVESPEEDS>0 then wavespeeds stored globally
int fluxcalc_standard(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time,int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  int getplpr(int dir, SFTYPE time, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
  FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP];
  //  int odir1,odir2;
  int Nvec[NDIM];
  int realisinterp;
  int dointerpolation;
  void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop, int *didassigngetstatecentdata);
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);




  ///////////////////////////////////////////////
  //
  // setup 1D reduction of EMF,flux calculation
  // here dir is face dir
  //  odir1=dir%3+1;
  //  odir2=(dir+1)%3+1;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  set_normal_realisinterp(&realisinterp);



  ////////////////////////////////////////////
  //
  // skip to next dir if no such dimension
  //
  /////////////////////////////////////////////

  // if nothing to interpolate, then quit
  if(npr2interpstart<=npr2interpend && Nvec[dir]>1){
    dointerpolation=1;
  }
  else{
    // just get loops and copy over to gp_?
    dointerpolation=0;

    //do limited number of things when not interpolating that dimension
    do_noninterpolation_dimension(ORIGINALFLUXCALC, dointerpolation,  realisinterp, dir, idel, jdel, kdel, pr, pl_ct, pr_ct, cent2faceloop, didassigngetstatecentdata);

    // skip real interpolation since no variation in that direction
    return(0);
  }




  /////////////////////////////
  //
  // setup wavedt
  //
  /////////////////////////////
  waveglobaldti[dir]=-100;
  waveglobaldtj[dir]=-100;
  waveglobaldtk[dir]=-100;



  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  if(npr2interpstart<=npr2interpend){
    // assume if DOEXTRAINTERP==1, then must get here
    p2interp=GLOBALPOINT(prc); // it's different
  }
  rescale_calc_full(dir,pr,p2interp);
#else
  p2interp=pr; // it's itself
#endif



  // STOREWAVESPEEDS==0 : use very local estimate at fluxes FACE?
  // STOREWAVESPEEDS==1 : use original pre-computed (just below) wavespeed at CENT that is "max-averaged" to FACE?
  // STOREWAVESPEEDS==2 : use very local estimate at fluxes FACE?, but compute and store during loop 
  if(STOREWAVESPEEDS==1){

#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
        MYFUN(get_global_wavespeeds_full(dir,is,ie,js,je,ks,ke,idel,jdel,kdel,POINT(pr),GLOBALPOINT(wspeed)),"flux.c:fluxcalc_standard()", "get_global_wavespeeds_full()", 0);
      }
  } // end if storing wavespeeds








  /////////////////////////////////////
  //
  // evaluate slopes (dq) or get pleft/pright of (possibly rescaled) primitive variables
  // c2e reconstruction: p2interp -> pleft & pright (indexed by grid cell #) -- atch comment
  //
  // pleft,pright also considered very temporary and so ok to be global (unless some computational overlapping is desired)
  //
  /////////////////////////////////////
  slope_lim(dointerpolation,realisinterp,dir,idel,jdel,kdel,pr,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),cent2faceloop);







  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

  //#if((SIMULBCCALC==2)&&(TYPE2==1))
  //  COMPFZLOOP(is,js,ks){
  //#else
  //  COMPZSLOOP( is, ie, js, je, ks, ke ) {{
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // requires full information
  {
    int i,j,k,pl,pliter;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
    FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
    FTYPE *p2interp_l,*p2interp_r;
    FTYPE dtij;
    FTYPE ctop,ctoprad;
    int computewithleft,computewithright;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
    // Setup rescale pointer reference
#if(RESCALEINTERP)
    p2interp_l=pstore_l; // then need separate storage
    p2interp_r=pstore_r;
#else
    p2interp_l=p_l; // p2interp_? is final p_?
    p2interp_r=p_r;
#endif


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
  


      ////////////////////////
      //
      // get the geometry for the flux face
      //
      ///////////////////////
      get_geometry(i, j, k, face, ptrgeom);


      //////////////////////////////////
      //
      // use p2interp,dq,pleft,pright to get p_l and p_r
      //
      /////////////////////////////////

      if(npr2interpstart<=npr2interpend){
        MYFUN(getplpr(dir,time,idel,jdel,kdel,i,j,k,ptrgeom,pr,pstag,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
#if(SPLITNPR || FIELDSTAGMEM)
        // then means there is going to be a second pass, so store into memory
        PINTERPLOOP(pliter,pl){
          MACP1A1(pl_ct,dir,i,j,k,pl)=p_l[pl];
          MACP1A1(pr_ct,dir,i,j,k,pl)=p_r[pl];
        }
        PNOTINTERPLOOP(pliter,pl){ // restore those things didn't interpolate
          p_l[pl]=MACP1A1(pl_ct,dir,i,j,k,pl);
          p_r[pl]=MACP1A1(pr_ct,dir,i,j,k,pl);
        }
#endif
      }
      else{
#if(SPLITNPR || FIELDSTAGMEM)
        // GODMARK: for now assume interpolations either all done or none done, else need exclusion list for interpolations
        // then just get from memory
        PLOOPALLINTERP(pl){
          p_l[pl]=MACP1A1(pl_ct,dir,i,j,k,pl);
          p_r[pl]=MACP1A1(pr_ct,dir,i,j,k,pl);
        }
#endif
#if(SPLITNPR==0 || FIELDSTAGMEM==0)
        dualfprintf(fail_file,"Should not be using gp_? in flux.c when SPLITNPR==0 or FIELDSTAGMEM==0\n");
        myexit(16760276);
#endif
      }



      // determine whether should compute things using left/right states
      leftrightcompute(i, j, k, dir, is, ie, js, je, ks, ke, &computewithleft, &computewithright);



      /////////////////////////////////////
      //
      // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
      //
      /////////////////////////////////////
#if(STOREFLUXSTATE)
      if(computewithleft) compute_and_store_fluxstate(dir, ISLEFT, i, j, k, ptrgeom, p_l);
      if(computewithright) compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, ptrgeom, p_r);
      // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif



      if(computewithleft&&computewithright){

        //////////////////////////////////
        //
        // actually compute the dissipative flux
        //
        /////////////////////////////////

        if(splitmaem==0){
          MYFUN(flux_compute_general(i, j, k, dir, ptrgeom, CUf,  MAC(pr,i,j,k), p_l, p_r, MAC(F,i,j,k), &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
        }
        else{
          MYFUN(flux_compute_splitmaem(i, j, k, dir, ptrgeom, CUf,  MAC(pr,i,j,k), p_l, p_r, MAC(F,i,j,k), MAC(FEM,i,j,k), &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
        }


        /////////////////////////////
        //
        // evaluate restriction on timestep
        //
        ///////////////////////////////
        // below is old before new grid sectioning
        //#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON==2))
        // avoid limiting dt if inside horizon
        // GODMARK: can only do this if boundary condition does not drive solution's dt behavior
        //    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
        //#endif


        ////////////////////////
        // set dt for this cell in this direction
        dtij = cour * dx[dir] / ctop;


        /////////////////////////////////
        //
        // save minimum dt
        //
        /////////////////////////////////
        if(PERCELLDT==0){
          // only set timestep if in computational domain or just +-1 cell beyond.  Don't go further since end up not really using that flux or rely on the stability of fluxes beyond that point.
          // Need +-1 in case flow is driven by injection boundary conditions rather than what's on grid
          if(WITHINACTIVESECTIONEXPAND1(i,j,k)){
#pragma omp critical
            {// *ndt and waveglobaldt's have to have blocked write access for OpenMP
              if (dtij < *ndt){
                *ndt = dtij;
                // below are global so can report when other dt's are reported in advance.c
                waveglobaldti[dir]=i;
                waveglobaldtj[dir]=j;
                waveglobaldtk[dir]=k;
              }
            }// end critical region
          }// end if within dt-setting section
        }
        else{
          GLOBALMACP0A1(dtijk,i,j,k,dir) = dtij;
        }




      }// end if doing computing using both left+right states (required also for ctop to exist)


    }// end 3D loop
  }// end parallel region

  


  return (0);
}






/// assign global variables needed when not doing certain dimensions
/// For example, to keep algorithm simple, always assume all directions of pl_ct and pr_ct exist, but if doing 2D simulation then we don't fill "left-right" states, so instead we need to compute or assign it so that data exists
/// also for the CTSTAG method, need loop information about going from center to face even if that dimension doesn't exist
void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  int i,j,k;
  int pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);
  void compute_and_store_fluxstate_assign(int dimeninput, int dimenoutput, int isleftrightinput, int isleftrightoutput, int i, int j, int k);


  if(FLUXB==FLUXCTSTAG){
    // if storing pl_ct and pr_ct, then copy 1D result in case used in some way not associated with whether dimension exists or not (not expensive)
    // For example, this is used for FLUXCTSTAG method in case dimension doesn't exist (Nvec[dir]==1) but still access gp_{l,r} instead of special conditions
    ////////COMPFULLLOOP PINTERPLOOP(pliter,pl) MACP1A1(pl_ct,dir,i,j,k,pl)=MACP1A1(pr_ct,dir,i,j,k,pl)=MACP0A1(pr,i,j,k,pl);
    copy_3dnpr2interp_2ptrs_fullloop(pr,pl_ct[dir],pr_ct[dir]);
  }
 
  
  if(SPLITNPR||FLUXB==FLUXCTSTAG){
    // just get loop information
    if(whichfluxcalc==ORIGINALFLUXCALC){
      slope_lim(dointerpolation, realisinterp,dir,idel,jdel,kdel,NULL,NULL,NULL,NULL,NULL,cent2faceloop);
    }
    else{
      slope_lim_cent2face(dointerpolation, realisinterp,dir,idel,jdel,kdel,NULL,NULL,NULL,NULL,NULL,cent2faceloop);
    }
  }
  

}



/// standard (non-field flux) calculation but setup to store results of interpolation so can be used for fluxctstag calculation
/// set pl_ct and pr_ct with FACE interpolations from CENT (including field face from pstagscratch[])
/// At end global pleft,pright,dq also set and if STOREWAVESPEEDS>0 then wavespeeds stored globally
int fluxcalc_standard_4fluxctstag(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  //  int odir1,odir2;
  int Nvec[NDIM];
  int realisinterp;
  int dointerpolation;
  void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop, int *didassigngetstatecentdata);





  ///////////////////////////////////////////////
  //
  // setup 1D reduction of EMF,flux calculation
  // here dir is face dir
  //  odir1=dir%3+1;
  //  odir2=(dir+1)%3+1;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  // only do interpolation if some type of interpolation to be done
  set_normal_realisinterp(&realisinterp);
  if(FLUXB==FLUXCTSTAG){
    realisinterp=0; // override since not interpolating p[B1+dir-1]
  }
  

  ////////////////////////////////////////////
  //
  // skip to next dir if no such dimension
  //
  /////////////////////////////////////////////
  if(npr2interpstart<=npr2interpend && Nvec[dir]>1){
    dointerpolation=1;
  }
  else{
    // just get loops, nothing to copy
    dointerpolation=0;

    //do limited number of things when not interpolating that dimension
    do_noninterpolation_dimension(NEWFLUXCALC, dointerpolation,  realisinterp, dir, idel, jdel, kdel, pr, pl_ct, pr_ct, cent2faceloop,didassigngetstatecentdata);

    // nothing else to do
    return(0);
  }






  /////////////////////////////
  //
  // setup wavedt
  //
  /////////////////////////////
  waveglobaldti[dir]=-100;
  waveglobaldtj[dir]=-100;
  waveglobaldtk[dir]=-100;


  

  //////////////////////////
  //
  // get wavespeeds if storing
  //
  ////////////////////////////

  if(STOREWAVESPEEDS==1){
#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
        MYFUN(get_global_wavespeeds_full(dir,is,ie,js,je,ks,ke,idel,jdel,kdel,POINT(pr),GLOBALPOINT(wspeed)),"flux.c:fluxcalc_standard()", "get_global_wavespeeds_full()", 0);
      }
  } // end if storing wavespeeds






  //////////////////////////
  //
  // obtain pl_ct and pr_ct (point face quantities) from pr (point centered quantity)
  //
  ////////////////////////////
  interpolate_prim_cent2face(stage, realisinterp, pr, pstag, pl_ct, pr_ct, dir, time, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dq, cent2faceloop);





  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

  //#if((SIMULBCCALC==2)&&(TYPE2==1))
  //  COMPFZLOOP(is,js,ks)
  //#else
  //#endif
  ////  COMPZSLOOP( is, ie, js, je, ks, ke ){
#if(STOREFLUXSTATE)
#pragma omp parallel  // then flux_compute() below uses *stored* state already
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // requires full information
#endif
  {
    int i, j, k;
    FTYPE dtij;
    FTYPE ctop;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    int computewithleft,computewithright;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
    // generally ptr's are different inside parallel block
    ptrgeom=&geomdontuse;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);



      // determine whether should compute things using left/right states
      leftrightcompute(i, j, k, dir, is, ie, js, je, ks, ke, &computewithleft, &computewithright);




      if(computewithleft&&computewithright){

        ////////////////////////
        //
        // get the geometry for the flux face
        //
        ///////////////////////
        get_geometry(i, j, k, face, ptrgeom); // OPTMARK: Seems should put back together interpolate_prim_cent2face() with flux_compute stuff below so only 1 call to get_geometry @ face....not sure why this way?!


        //////////////////////////////////
        //
        // actually compute the dissipative flux
        //
        /////////////////////////////////

        if(splitmaem==0){
          MYFUN(flux_compute_general(i, j, k, dir, ptrgeom, CUf,  MAC(pr,i,j,k), MACP1A0(pl_ct,dir,i,j,k), MACP1A0(pr_ct,dir,i,j,k), MAC(F,i,j,k), &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
        }
        else{
          MYFUN(flux_compute_splitmaem(i, j, k, dir, ptrgeom, CUf,  MAC(pr,i,j,k), MACP1A0(pl_ct,dir,i,j,k), MACP1A0(pr_ct,dir,i,j,k), MAC(F,i,j,k), MAC(FEM,i,j,k), &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
        }


        /////////////////////////////
        //
        // evaluate restriction on timestep
        //
        ///////////////////////////////
        // below is old before new grid sectioning
        //#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON==2))
        //    // avoid limiting dt if inside horizon
        //    // GODMARK: can only do this if boundary condition does not drive solution's dt behavior
        //    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
        //#endif


        ////////////////////////
        // set dt for this cell in this direction
        dtij = cour * dx[dir] / ctop;


        /////////////////////////////////
        //
        // save minimum dt
        //
        /////////////////////////////////
        if(PERCELLDT==0){

          // only set timestep if in computational domain or just +-1 cell beyond.  Don't go further since end up not really using that flux or rely on the stability of fluxes beyond that point.
          // Need +-1 in case flow is driven by injection boundary conditions rather than what's on grid
          if(WITHINACTIVESECTIONEXPAND1(i,j,k)){
#pragma omp critical
            {// *ndt and waveglobaldt's have to have blocked write access for OpenMP
              if (dtij < *ndt){
                *ndt = dtij;
                // below are global so can report when other dt's are reported in advance.c
                waveglobaldti[dir]=i;
                waveglobaldtj[dir]=j;
                waveglobaldtk[dir]=k;
              }
            }// end critical region
          }// end if within dt-setting section
        }
        else{
          GLOBALMACP0A1(dtijk,i,j,k,dir) = dtij;
        }



      }// end if doing computing using both left+right states (required also for ctop to exist)


    }// end FLUX LOOP
  }// end parallel region





  return (0);
}










/// normal interpolation of CENT quantities to FACE quantities
/// sets global variables pl_ct and pr_ct to p_l and p_r from interpolations
int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, SFTYPE time, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop)
{
  void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  int getplpr(int dir, SFTYPE time, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);
  FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP];
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int realis,realie,realjs,realje,realks,realke;
  int dointerpolation;





  // if inside this function, then doing interpolation
  dointerpolation=1;


  /////////////////////////////////////
  //
  // setup interpolation so avoids staggered field for field along "dir" direction
  // avoid magnetic field interpolation along "dir" direction using slope_lim() that is for cent->edges only
  // inside getplpr() staggered field is assigned to final p_l p_r before flux is computed
  //
  /////////////////////////////////////
  if(FLUXB==FLUXCTSTAG){
    int pl,pl2;
    ////////////////////////////////////////////
    //
    // save choice for interpolations
    nprlocalstart=npr2interpstart;
    nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];


#pragma omp parallel private(pl,pl2)
    { // must set npr2interp stuff inside parallel region since threadprivate
      // choice for range of PLOOPINTERP
      // check if wanted to interpolate B along dir, and if so remove
      for(pl=npr2interpstart;pl<=npr2interpend;pl++){
        if(npr2interplist[pl]==B1+dir-1){
          for(pl2=pl+1;pl2<=npr2interpend;pl2++) npr2interplist[pl2-1]=npr2interplist[pl2]; // moving upper to lower index
          npr2interpend--; // lost field along dir, so one less thing to do
          break;
        }
      }
    }

  }



  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  // assume if DOEXTRAINTERP==1, then must get here
  p2interp=GLOBALPOINT(prc); // it's different
  rescale_calc_full(dir,pr,p2interp);
#else
  p2interp=pr; // it's itself
#endif






  /////////////////////////////////////
  //
  // evaluate slopes (dq) or get pleft/pright of (possibly rescaled) primitive variables
  // c2e reconstruction: p2interp -> pleft & pright (indexed by grid cell #) -- atch comment
  //
  /////////////////////////////////////
  slope_lim_cent2face(dointerpolation,realisinterp,dir,idel,jdel,kdel,pr,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),cent2faceloop);


  // override normal fluxloop
  if(extrazones4emf){
    // don't really need fluxes in this domain, but do need interpolated face values to be transferred from pleft/pright to pl_ct pr_ct
    realis=emfUconsloop[FIS]-SHIFT1;
    realie=emfUconsloop[FIE]+SHIFT1;
    realjs=emfUconsloop[FJS]-SHIFT2;
    realje=emfUconsloop[FJE]+SHIFT2;
    realks=emfUconsloop[FKS]-SHIFT3;
    realke=emfUconsloop[FKE]+SHIFT3;
  }
  else{ // otherwise normal
    realis=is;
    realie=ie;
    realjs=js;
    realje=je;
    realks=ks;
    realke=ke;
  }



  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

  //#if((SIMULBCCALC==2)&&(TYPE2==1))
  //  COMPFZLOOP(realis,realjs,realks)
  //#else
  //#endif
  ///// COMPZSLOOP( realis, realie, realjs, realje, realks, realke ){

#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // storage of state requires full copyin()
  {
    FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
    FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
    FTYPE *p2interp_l,*p2interp_r;
    int pl,pliter;
    int i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    int computewithleft,computewithright;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP( realis, realie, realjs, realje, realks, realke );
    
    // Setup rescale pointer reference
#if(RESCALEINTERP)
    p2interp_l=pstore_l; // then need separate storage
    p2interp_r=pstore_r;
#else
    p2interp_l=p_l; // p2interp_? is final p_?
    p2interp_r=p_r;
#endif


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      ////////////////////////
      //
      // get the geometry for the flux face
      //
      ///////////////////////
      get_geometry(i, j, k, face, ptrgeom);


      //////////////////////////////////
      //
      // use p2interp,dq,pleft,pright to get p_l and p_r
      //
      /////////////////////////////////

      MYFUN(getplpr(dir,time,idel,jdel,kdel,i,j,k,ptrgeom,pr,pstag,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
      if(SPLITNPR || FLUXB==FLUXCTSTAG){
        // then means there is going to be a second pass, so store into memory
        PINTERPLOOP(pliter,pl){
          MACP1A1(pl_ct,dir,i,j,k,pl)=p_l[pl];
          MACP1A1(pr_ct,dir,i,j,k,pl)=p_r[pl];
        }
      }
      if(FLUXB==FLUXCTSTAG){
        // then also get B in dir direction not included in interpolation in slope_lim() above but included in getplpr() from pstagscratch[]
        pl = B1+dir-1;
        MACP1A1(pl_ct,dir,i,j,k,pl)=p_l[pl];
        MACP1A1(pr_ct,dir,i,j,k,pl)=p_r[pl];
      }


      // determine whether should compute things using left/right states
      // note that should only do non-left non-right state stuff (like computing dt) if both left+right states being done
      leftrightcompute(i, j, k, dir, is, ie, js, je, ks, ke, &computewithleft, &computewithright);


      /////////////////////////////////////
      //
      // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
      //
      /////////////////////////////////////
#if(STOREFLUXSTATE)
      if(computewithleft) compute_and_store_fluxstate(dir, ISLEFT, i, j, k, ptrgeom, p_l);
      if(computewithright) compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, ptrgeom, p_r);
      // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif




    }// end COMPZLOOP
  }// end parallel region







#if(0)
  // DEBUG: (didn't matter)
  bound_prim(STAGEM1,pl_ct[dir]);
  bound_prim(STAGEM1,pr_ct[dir]);
#endif




  if(FLUXB==FLUXCTSTAG){
    int pl,pliter;
    ////////////////////////////////////////////
    //
    // restore choice for interpolations
#pragma omp parallel private(pl)
    { // must set npr2interp stuff inside parallel region since threadprivate
      npr2interpstart=nprlocalstart;
      npr2interpend=nprlocalend;
      PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
    }
  }
  

  return(0);

}



/// assign one fluxstate to another fluxstate (usually used to assign left->right or visa versa when dimension doesn't exist
/// Note that isleftright(input/output) can only be ISLEFT,ISRIGHT and not ISMIDDLE unless rewrote function to use fluxstatecent, but function not needed anymore
/// GODMARK: no longer using this since store cent and copy it to ISLEFT,ISRIGHT
void compute_and_store_fluxstate_assign(int dimeninput, int dimenoutput, int isleftrightinput, int isleftrightoutput, int i, int j, int k)
{

  // copy entire structure
  GLOBALMACP1A1(fluxstate,dimenoutput,i,j,k,isleftrightoutput)=GLOBALMACP1A1(fluxstate,dimeninput,i,j,k,isleftrightinput);

}

/// compute and store get_state() data for p_l and p_r type objects
/// need to call this whenever compute left-right primitives used for computing flux
void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr)
{
  int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);

  // get state
  pureget_stateforfluxcalc(pr,geom,&GLOBALMACP1A1(fluxstate,dimen,i,j,k,isleftright));

}



/// compute and store get_state() data for centered state
void compute_and_store_fluxstatecent(FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  int pureget_stateforfluxcalcorsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int is,ie,js,je,ks,ke,di,dj,dk;
  //  FTYPE (*shocktemparray)[NSTORE2][NSTORE3][NPR];
  FTYPE (*shocktemparray)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  int startorderi,endorderi;
  extern FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P);
  extern FTYPE  Divcalc(int dir, FTYPE Fi, FTYPE *V, FTYPE *P);



  const int Nvec[NDIM]={0,N1,N2,N3};
  const int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1};

  // setup temporary space
  shocktemparray = GLOBALPOINT(emf);

  // define +-1 in every direction loop range
  int loc=CENT;
  int continuous=0;
  set_interppoint_loop_ranges_3Dextended(ENOINTERPTYPE, loc, continuous, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);


#if(STOREFLUXSTATE||STORESHOCKINDICATOR)

  /////////////////////////
  //
  // loop over all i,j,k
  //
  /////////////////////////
  //  COMPZSLOOP(is,ie,js,je,ks,ke){
  //  COMPFULLLOOP{ // NOTE: Using COMPFULLLOOP now since using this CENT for many things, including global wave speed calculation that currently uses COMPFULLOOP
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // requires full copyin()
  {
    int i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    int dir;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
    

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

 

      // set geometry for centered zone to be updated
      get_geometry(i, j, k, CENT, ptrgeom);

      // get state
      pureget_stateforfluxcalcorsource(MAC(pr,i,j,k),ptrgeom,&GLOBALMAC(fluxstatecent,i,j,k));

      //      FTYPE flux[NDIM];
      //      mhd_calc(MAC(pr,i,j,k),TT,ptrgeom,&GLOBALMAC(fluxstatecent,i,j,k),flux,NULL);

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
      // get true directional energy flux instead of fake flux
      //      FTYPE fluxrad[NDIM];
      //      mhd_calc_rad(MAC(pr,i,j,k),TT,ptrgeom,&GLOBALMAC(fluxstatecent,i,j,k),fluxrad,NULL);
#endif


#if(STORESHOCKINDICATOR)
      // see more details in interpline.c:get_V_and_P()


      DIMENLOOP(dir){
        if(NxNOT1[dir]){
#if(VLINEWITHGDETRHO==0)
          MACP1A0(shocktemparray,SHOCKPLSTOREVEL1+dir-1,i,j,k)=MACP0A1(pr,i,j,k,UU+dir);
#else
          //MACP1A0(shocktemparray,SHOCKPLSTOREVEL1+dir-1,i,j,k)=(ptrgeom->gdet)*MACP0A1(pr,i,j,k,RHO)*(GLOBALMAC(fluxstatecent,i,j,k).ucon[dir]);
          MACP1A0(shocktemparray,SHOCKPLSTOREVEL1+dir-1,i,j,k)=(ptrgeom->gdet)*(GLOBALMAC(fluxstatecent,i,j,k).ucon[dir]);
          //          MACP1A0(shocktemparray,SHOCKPLSTOREVEL1+dir-1,i,j,k) = (ptrgeom->gdet)*(-flux[dir]);
#endif
          
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
#if(VLINEWITHGDETRHO==0)
          MACP1A0(shocktemparray,SHOCKRADPLSTOREVEL1+dir-1,i,j,k) = MACP0A1(pr,i,j,k,URAD0+dir);
#else
          //MACP1A0(shocktemparray,SHOCKRADPLSTOREVEL1+dir-1,i,j,k) = (ptrgeom->gdet)*MACP0A1(pr,i,j,k,URAD0)*(GLOBALMAC(fluxstatecent,i,j,k).uradcon[dir]);
          MACP1A0(shocktemparray,SHOCKRADPLSTOREVEL1+dir-1,i,j,k) = (ptrgeom->gdet)*(GLOBALMAC(fluxstatecent,i,j,k).uradcon[dir]);
          //MACP1A0(shocktemparray,SHOCKRADPLSTOREVEL1+dir-1,i,j,k) = (ptrgeom->gdet)*(-fluxrad[dir]);
#endif
#endif
        }
      }


      // get total MHD pressure
      FTYPE pmhd=GLOBALMAC(fluxstatecent,i,j,k).pressure + 0.5*GLOBALMAC(fluxstatecent,i,j,k).bsq;
      MACP1A0(shocktemparray,SHOCKPLSTOREPTOT,i,j,k)=pmhd;

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
      FTYPE prad=(4.0/3.0-1.0)*MACP0A1(pr,i,j,k,PRAD0);  // KORALNOTE: recall pressure just along diagonal and no velocity in R^\mu_\nu
      MACP1A0(shocktemparray,SHOCKRADPLSTOREPTOT,i,j,k) = prad;

      // add radiation pressure to total pressure if optically thick
      FTYPE tautot[NDIM],tautotmax;
      // &GLOBALMAC(fluxstatecent,i,j,k)
      calc_tautot(&MACP0A1(pr,i,j,k,0), ptrgeom, NULL, tautot, &tautotmax); // high accuracy not required

      MACP1A0(shocktemparray,SHOCKPLSTOREPTOT,i,j,k) += MIN(tautotmax/TAUTOTMAXSWITCH,1.0)*prad;
      MACP1A0(shocktemparray,SHOCKRADPLSTOREPTOT,i,j,k) += MIN(tautotmax/TAUTOTMAXSWITCH,1.0)*pmhd;
#endif

#endif // end if STORESHOCKINDICATOR



#if(STOREFLUXSTATE)

      // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen][ISLEFT,ISRIGHT]
      // If doing staggered field method, then also assumes left,right stored like pl_ct and pr_ct are stored
      // note we are copying entire structure here

      DIMENLOOP(dir){

        if(NxNOT1[dir]==0 || FIELDSTAGMEM){
          // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen][ISLEFT,ISRIGHT]
          GLOBALMACP1A1(fluxstate,dir,i,j,k,ISLEFT)=GLOBALMAC(fluxstatecent,i,j,k);
          GLOBALMACP1A1(fluxstate,dir,i,j,k,ISRIGHT)=GLOBALMAC(fluxstatecent,i,j,k);
        }
      }

#endif// end if STOREFLUXSTATE

    }// end 3D loop
  }// end parallel region  
  
#endif// end if STOREFLUXSTATE|| STORESHOCKINDICATOR






  {// begin block



#if(STORESHOCKINDICATOR)
    // OPTMARK: Consider storing temporary velvec aligned with extraction direction (dir) .  Although ptot scalar has no preferred choice and final pr should have standard choice.
    // Ficalc() requires vel and p exist +-2 beyond where shock indicator is computed.
    // Assume shock indicator required at largest possible domain, which is then 2-inwards from outer boundaries only in each direction.
    // shockindicatorarray[] computed here is used in interpline.c when calling get_1d_line_shockarray() to get shockindicator[]
      
#if(MAXBND<4)
    dualfprintf(fail_file,"MAXBND should be 4 for shockindicator???\n");
    myexit(8465684);
#endif
  
    startorderi = - (7)/2; // order=7 fixed for shock detector
    endorderi   = - startorderi;

  

#pragma omp parallel
    {
      int dir;
      int imod,jmod,kmod;
      int gotgeometry=0;
      int i,j,k,l;
      int pl;
      OPENMP3DLOOPVARSDEFINE;
      //    FTYPE *primptr[MAXNPR]; // number of pointers
      FTYPE *velptr,*ptotptr;
      //    FTYPE a_primstencil[MAXNPR][MAXSPACEORDER];
      FTYPE a_velstencil[MAXSPACEORDER],a_ptotstencil[MAXSPACEORDER]; // Shouldn't need anymore than 10      


      // shift pointer
      //    PALLREALLOOP(pl){
      //      primptr[pl] = a_primstencil[pl] - startorderi;
      //    }
      velptr = a_velstencil - startorderi;
      ptotptr = a_ptotstencil - startorderi;


      //            if(dir==1) OPENMP3DLOOPSETUPFULLINOUT2DIR1;
      //            if(dir==2) OPENMP3DLOOPSETUPFULLINOUT2DIR2;
      //            if(dir==3) OPENMP3DLOOPSETUPFULLINOUT2DIR3;
      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        FTYPE radshock={0},mhdshock={0}; // KORALTODO: Not used yet, but can be roughly used to pull-over multiple dimensions for shock indicator
        
        DIMENLOOP(dir){
          if(NxNOT1[dir]){
            
            
            // extract stencil of data for Ficalc()
            for(l=startorderi;l<=endorderi;l++){
              imod=i+(dir==1)*l;
              jmod=j+(dir==2)*l;
              kmod=k+(dir==3)*l;

              imod=MAX(-N1BND,imod);
              imod=MIN(N1-1+N1BND,imod);

              jmod=MAX(-N2BND,jmod);
              jmod=MIN(N2-1+N2BND,jmod);

              kmod=MAX(-N3BND,kmod);
              kmod=MIN(N3-1+N3BND,kmod);


              // PALLREALLOOP(pl) primptr[pl][l] = MACP0A1(pr,imod,jmod,kmod,pl);
              velptr[l] = MACP1A0(shocktemparray,SHOCKPLSTOREVEL1+dir-1,imod,jmod,kmod);
              ptotptr[l] = MACP1A0(shocktemparray,SHOCKPLSTOREPTOT,imod,jmod,kmod);
            }
      
            FTYPE Fi;
            //      GLOBALMACP0A1(shockindicatorarray,i,j,k,SHOCKPLDIR1+dir-1)=Ficalc(dir,&velptr[0],&ptotptr[0],&primptr[0]);
            Fi=Ficalc(dir,&velptr[0],&ptotptr[0]);
            //            Fi=1.0;
            if(NSPECIAL>=1&&0 && DODISSMEASURE){
              FTYPE dissit=fabs(GLOBALMACP0A1(dissmeasurearray,i,j,k,1));
              dissit = MAX(0.0,MIN(1,2*(dissit-0.3)));
              Fi=MIN(1.0,MAX(Fi,dissit));
            }
            if(steppart>0) Fi=MAX(Fi,GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k));
            GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k)=Fi;

            // DEBUG
            if(DODISSMEASURE){
              GLOBALMACP0A1(dissmeasurearray,i,j,k,NSPECIAL+1+dir-1)=GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k);
            }
            //            if(i==399){
            //              for(l=startorderi;l<=endorderi;l++){
            //                dualfprintf(fail_file,"GAS: l=%d vel=%g ptot=%g\n",l,velptr[l],ptotptr[l]);
            //              }                
            //            }


            if(DIVERGENCEMETHOD==DIVMETHODPREFLUX){
              GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k)=Divcalc(dir,Fi,&velptr[0],&ptotptr[0]);
            }

            if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE){
              // extract stencil of data for Ficalc()
              for(l=startorderi;l<=endorderi;l++){
                imod=i+(dir==1)*l;
                jmod=j+(dir==2)*l;
                kmod=k+(dir==3)*l;

                imod=MAX(-N1BND,imod);
                imod=MIN(N1-1+N1BND,imod);

                jmod=MAX(-N2BND,jmod);
                jmod=MIN(N2-1+N2BND,jmod);

                kmod=MAX(-N3BND,kmod);
                kmod=MIN(N3-1+N3BND,kmod);

                // PALLREALLOOP(pl) primptr[pl][l] = MACP0A1(pr,imod,jmod,kmod,pl);
                velptr[l] = MACP1A0(shocktemparray,SHOCKRADPLSTOREVEL1+dir-1,imod,jmod,kmod);
                ptotptr[l] = MACP1A0(shocktemparray,SHOCKRADPLSTOREPTOT,imod,jmod,kmod);
              }
                
              //      GLOBALMACP0A1(shockindicatorarray,i,j,k,SHOCKRADPLDIR1+dir-1)=Ficalc(dir,&velptr[0],&ptotptr[0],&primptr[0]);
              FTYPE Firad;
              Firad=Ficalc(dir,&velptr[0],&ptotptr[0]);
              //              Firad*=2.0;
              //              Firad=MIN(1.0,Firad);
              if(NSPECIAL>=6&&0 && DODISSMEASURE){
                FTYPE dissit=fabs(GLOBALMACP0A1(dissmeasurearray,i,j,k,5));
                dissit = MAX(0.0,MIN(1,2*(dissit-0.3)));
                Firad=MIN(1.0,MAX(Firad,dissit));
              }
              if(steppart>0) Firad=MAX(Firad,GLOBALMACP1A0(shockindicatorarray,SHOCKRADPLDIR1+dir-1,i,j,k));
              GLOBALMACP1A0(shockindicatorarray,SHOCKRADPLDIR1+dir-1,i,j,k)=Firad;

              //DEBUG
              //              dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d dir=%d Fi=%g Firad=%g\n",nstep,steppart,i,dir,Fi,Firad);
              //              if(i==399){
              //                for(l=startorderi;l<=endorderi;l++){
              //                  dualfprintf(fail_file,"RAD: l=%d vel=%g ptot=%g\n",l,velptr[l],ptotptr[l]);
              //                }                
              //              }

              // DEBUG
              if(DODISSMEASURE){
                GLOBALMACP0A1(dissmeasurearray,i,j,k,NSPECIAL+1+3+dir-1)=GLOBALMACP1A0(shockindicatorarray,SHOCKRADPLDIR1+dir-1,i,j,k);
              }


              if(DIVERGENCEMETHOD==DIVMETHODPREFLUX){
                GLOBALMACP1A0(shockindicatorarray,DIVRADPLDIR1+dir-1,i,j,k)=Divcalc(dir,Firad,&velptr[0],&ptotptr[0]);
              }
            }// end if doing radiation
          }// end if doing dir
        }//end over dir loop

        // set geometry for centered zone to be updated
        struct of_geom geomdontuse;
        struct of_geom *ptrgeom=&geomdontuse;
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
        get_geometry(i, j, k, CENT, ptrgeom);
        gotgeometry=1;

        // use maximum shock indicator if optically thick
        // KORALTODO: using tau here, but tau is only over a cell.  If resolution changes, different tau, yet physics is the same.  How to manage?
        // KORALNOTE: Need this because radiation can be escaping from shock, so radiation itself is not converging, while still sufficient for MHD shock conditions.
        FTYPE tautot[NDIM],tautotmax;
        // &GLOBALMAC(fluxstatecent,i,j,k)
        calc_tautot(&MACP0A1(pr,i,j,k,0), ptrgeom, NULL, tautot, &tautotmax); // high accuracy not required

        DIMENLOOP(dir){
          if(NxNOT1[dir]){
            FTYPE Firad=GLOBALMACP1A0(shockindicatorarray,SHOCKRADPLDIR1+dir-1,i,j,k);
            FTYPE Fi=GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k);
            FTYPE Fiswitch =MIN(tautotmax/TAUTOTMAXSWITCH,1.0);

            GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k) = Fiswitch*MIN(1.0,MAX(Firad,Fi)) + (1.0-Fiswitch)*Fi;
            GLOBALMACP1A0(shockindicatorarray,SHOCKRADPLDIR1+dir-1,i,j,k) = Fiswitch*MIN(1.0,MAX(Firad,Fi)) + (1.0-Fiswitch)*Firad;



            if(DIVERGENCEMETHOD==DIVMETHODPREFLUX){
              // don't consider radiation divergence, and certainly  not as if shock indiator value.
              //            GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k) = MAX(MIN(tautotmax/TAUTOTMAXSWITCH,1.0)*GLOBALMACP1A0(shockindicatorarray,DIVRADPLDIR1+dir-1,i,j,k),GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k));
              //            GLOBALMACP1A0(shockindicatorarray,DIVRADPLDIR1+dir-1,i,j,k) = MAX(MIN(tautotmax/TAUTOTMAXSWITCH,1.0)*GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k),GLOBALMACP1A0(shockindicatorarray,DIVRADPLDIR1+dir-1,i,j,k));
            }
          }
        }
#endif


        ////////
        //
        // compute best indication that entropy equation should be used based upon smoothness of flow
        //
        /////////

        if(DIVERGENCEMETHOD==DIVMETHODPREFLUX){

          // set geometry for centered zone to be updated
          if(gotgeometry==0) get_geometry(i, j, k, CENT, ptrgeom);

          // first see which divergence is largest in absolute terms accounting for dimensions
          FTYPE divcond=-1.0,absdivcond=0.0; // 0.0 here is crucial as reference for divergent flow
          FTYPE divcondtest[NDIM],absdivcondtest;
          int largestdir=-1;
          DIMENLOOP(dir){
            if(NxNOT1[dir]){
              if(VLINEWITHGDETRHO==0){
                divcondtest[dir]=GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k)/(sqrt(fabs(ptrgeom->gcon[GIND(dir,dir)])));
              }
              else{
                divcondtest[dir]=GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k)/(sqrt(fabs(ptrgeom->gcon[GIND(dir,dir)]))*ptrgeom->gdet);
              }
              absdivcondtest=fabs(divcondtest[dir]);
              if(absdivcondtest>absdivcond){
                largestdir=dir;
                absdivcond=absdivcondtest;
                divcond=divcondtest[dir];
              }
            }
          }
          // ensure largest divergence dominates
          DIMENLOOP(dir){
            if(NxNOT1[dir]){
#define MAXSHOCK (0.1)
#define FACTORBIGDIV (10.0)
              // only consider fluid, since if optically thick already accounted for radiation in this indicator
              if(GLOBALMACP1A0(shockindicatorarray,SHOCKPLDIR1+dir-1,i,j,k)>MAXSHOCK){
                divcond=-BIG;
              }
              // check that divergence is some factor more than convergence in any other direction
              if(dir!=largestdir && fabs(divcond)<FACTORBIGDIV*fabs(divcondtest[dir]) && divcondtest[dir]<0.0){
                divcond=-BIG;
              }
            }
          }
          // divcond now holds divergence condition for largest divergent part of flow.  Replace all dimensions with this uni-dimensional condition
          DIMENLOOP(dir){
            if(NxNOT1[dir]){
              GLOBALMACP1A0(shockindicatorarray,DIVPLDIR1+dir-1,i,j,k)=divcond;
            }
          }
        }// end if(DIVERGENCEMETHOD==DIVMETHODPREFLUX)

      }// end 3D loop
    }// end parallel region







    // if want to use contact steepener methods, then store that too.
    //  if(CONTACTINDICATOR){
    ///////////////////
    //
    // 1D GET CONTACT INDICATOR
    //
    ////////////////////
    //use primitive values that correspond to the quantities being interpolated
    //  }



    // GODMARK: At end, could combine into multi-dimensional shock indicator.  But should have measure to neglect other minor dimensions where (e.g.) velocity or pressure jumps not important)


#endif // end if STORESHOCKINDICATOR
  }// end block



}








/// use dq,pleft,pright to obtain p_l and p_r for CENT to FACE
int getplpr(int dir, SFTYPE time, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *ptrgeom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r)
{
  void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r);
  int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r);
  int pl,pliter;


  
  //////////////////////////////////////
  //
  // interpolate primitive using slope (dq) or directly from pleft and pright
  // For FV: p_left, p_right (indexed by grid cell #) -> p2interp_l, p2interp_r (indexed by interface #) -- atch comment
  //
  // always do since p2interp is good and dq/pleft/pright should have stored quantities or just-computed quantities
  /////////////////////////////////////
  getp2interplr(dir,idel,jdel,kdel,i,j,k,p2interp,dq,pleft,pright,p2interp_l,p2interp_r);
   


  /////////////////////////////////////
  //      
  // after interpolation, unrescale from p2interp to normal primitive 
  // p2interp_? is p_? if not rescaling, so no need to assign p2interp_? to p_? if not rescaling
  //
  ///////////////////////////////////
  if(RESCALEINTERP && npr2interpstart<=npr2interpend){
    // only do if some interpolation done (consistent with no rescale() done in fluxcalc_standard() above
    // setup plausible p_l/p_r in case used for some reason (e.g. inversion starting guess)
    // this is sufficient for utoprim_1D...a better guess does no better (see interpU code)
    PINTERPLOOP(pliter,pl){
      p_l[pl]=MACP0A1(pr,i,j,k,pl);
      p_r[pl]=MACP0A1(pr,i,j,k,pl);
    }
    rescale(UNRESCALE,dir,p_l,ptrgeom,p2interp_l);
    rescale(UNRESCALE,dir,p_r,ptrgeom,p2interp_r);
  }




  //////////////////////////////
  //
  // Must preserve divb in 1D Riemann problem, so B^{dir} must be continuous
  //
  // GODMARK: Should really interpolate SUCH THAT this is automatically satisfied?
  // Yes, should use larger stencil and interpolate such that constant. Essentially choose
  // large stencil and effectively choosing which points we trust (not necessarily local points)
  //
  // GODMARK: This does not enforce E_\perp to be continuous for stationary flow!
  //
  ///////////////////////////////


  if(FLUXB==FLUXCTSTAG){
    // exactly correct that there is only 1 value
    pl = B1+dir-1;
    p_l[pl]=p_r[pl]=MACP0A1(pstag,i,j,k,pl);
  }
  else{
#if(BDIRCONT)
    // should really interpolate such that p_l=p_r
    pl = B1+dir-1;
    p_l[pl]=p_r[pl]=0.5*(p_l[pl]+p_r[pl]);
#endif
  }


#if(BOUNDPLPR)
  set_plpr(dir,time,i,j,k,pr,p_l,p_r);
#endif

  ///////////////////////
  //
  // correct interpolated quantities
  // no fixup accounting for these intermediate quantities
  //
  //////////////////////
  MYFUN(check_plpr(dir, i, j, k, idel, jdel, kdel, ptrgeom, pr, p_l, p_r),"step_ch.c:fluxcalc()", "check_plpr()", 1);


  return(0);

}





/// check that left and right primitives are valid
int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r)
{
  int pl,pliter;

#if(EVOLVECHECKS)
#if(WHICHVEL==VEL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_4vel(dir,p_l,NULL,-1);
  inflow_check_4vel(dir,p_r,NULL,geom,-1);
#endif
#elif(WHICHVEL==VEL3)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_3vel(dir,p_l,NULL,geom,-1);
  inflow_check_3vel(dir,p_r,NULL,geom,-1);
#endif
#if(JONCHECKS2)
  // must verify if this p makes sense (u^t sense)
  MYFUN(check_pr(p_l,MAC(pr,i-idel,j-jdel,k-kdel),NULL,geom,1,-1.0),"flux.c:check_plpr()", "check_pr()", 1); 
  
  MYFUN(check_pr(p_r,MAC(pr,i,j,k),NULL,geom,2,-1),"flux.c:check_plpr()", "check_pr()", 2);
#endif
#elif(WHICHVEL==VELREL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_rel4vel(dir,p_l,NULL,geom,-1);
  inflow_check_rel4vel(dir,p_r,NULL,geom,-1);
#endif
  // need to limit gamma since gamma may be large for interpolated value and would lead to bad fluxes
  MYFUN(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,p_l,NULL,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 1);  //jon corr, see email re: SPINT warnings from 4/24/2006 10:54 p.m.
  MYFUN(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,p_r,NULL,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 2);  //jon corr
#endif// end if WHICHVEL==VEL4REL      
#endif


#if(PRODUCTION==0 && MERGEDC2EA2CMETHOD==0)
  // NOTE: with merged method use expanded loops that will use bogus data, but done for simplicity.  So no nan check possible there -- assume fill bad regions with zeroes
  PINTERPLOOP(pliter,pl) if(!isfinite(p_l[pl])){
    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d j=%d k=%d :: p_l is not finite pl=%d p_l=%21.15g\n",nstep,steppart,i,j,k,pl,p_l[pl]);
  }
  PINTERPLOOP(pliter,pl) if(!isfinite(p_r[pl])){
    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d j=%d k=%d :: p_r is not finite pl=%d p_r=%21.15g\n",nstep,steppart,i,j,k,pl,p_r[pl]);
  }
#endif

  // DEBUG:
  //  dualfprintf(fail_file,"CHECK dir=%d :: i=%d j=%d k=%d\n",dir,i,j,k);


  return(0);
}







///
///
/// interpolate primitive using slope or just copy pleft/pright into the correct place
///
/// |=interface
/// i=zone center of ith zone
///
/// |              |     p2interp(i)    |
/// |              |       dq(i)        |
/// |        p_l(i)|p_r(i)   i          |
/// |              |pleft(i)   pright(i)|
///
///
///
//////////////////////////////////////
/// GODMARK: as Sasha mentions, for shifting stencil shouldn't just extrapolate value from non-local values, but use other slope and most LOCAL value to obtain extrapolation.
///
void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r)
{
  FTYPE Xcent[NDIM],Xleft[NDIM];
  FTYPE Vcent[NDIM],Vleft[NDIM];
  FTYPE rleft,rcenter,thleft,thcent;
  int pl,pliter;
  int locallim;
  int choose_limiter(int dir, int i, int j, int k, int pl);


  ////////////
  //
  // reshuffle dq's such that reconstruction within active zones does not use the boundary values
  //
  ////////////
  remapdq( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);


  ////////////
  //
  // Do interpolation
  //
  ////////////
  
  if((HORIZONSUPERFAST)&&((lim[dir]<PARA)&&(LIMADJUST==0))&&(dir==1)){
    // since this uses dq's to shift, must always have dq's everywhere.  Hence LIMADJUST==0 and lim<PARA must be true
    bl_coord_ijk_2(im1mac(i), j, k, CENT, Xleft,Vleft);
    bl_coord_ijk_2(i, j, k, CENT, Xcent,Vcent);
    rleft=Vleft[1];
    thleft=Vleft[2];
    rcenter=Vcent[1];
    thcent=Vcent[2];
    
    PINTERPLOOP(pliter,pl){
      locallim=choose_limiter(dir, i,j,k,pl);
      int usedq = usedqarray[locallim];
      // get interpolated quantity
      if(usedq){
        if(rleft>Rhor) p2interp_l[pl] = MACP0A1(p2interp,i - idel,j - jdel,k - kdel,pl) + 0.5 * MACP0A1(dq,i - idel,j - jdel,k - kdel,pl);
        else p2interp_l[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP0A1(dq,i,j,k,pl);
        if(rcenter>Rhor) p2interp_r[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP0A1(dq,i,j,k,pl);
        else p2interp_r[pl] = MACP0A1(p2interp,i+idel,j+jdel,k+kdel,pl) - 1.5 * MACP0A1(dq,i+idel,j+jdel,k+kdel,pl);
      }
      else{
        p2interp_l[pl] = MACP0A1(pright,i-idel,j-jdel,k-kdel,pl);
        p2interp_r[pl] = MACP0A1(pleft,i,j,k,pl);
      }
    } // end PLOOPINTERP
  }// if horizonsuperfast
  else{
    /////////////////////////////
    //
    // standard routine
    //
    /////////////////////////////
    PINTERPLOOP(pliter,pl){
      locallim=choose_limiter(dir, i,j,k,pl);
      int usedq = usedqarray[locallim];
      // get interpolated quantity
      if(usedq){
        p2interp_l[pl] = MACP0A1(p2interp,i - idel,j - jdel,k - kdel,pl) + 0.5 * MACP0A1(dq,i - idel,j - jdel,k - kdel,pl);
        p2interp_r[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP0A1(dq,i,j,k,pl);
      }
      else{
        //p_l & p_r for the current interface -- atch comment
        p2interp_l[pl] = MACP0A1(pright,i-idel,j-jdel,k-kdel,pl);
        p2interp_r[pl] = MACP0A1(pleft,i,j,k,pl);
      }
    }
  }



  //////////////
  //
  // for boundaries where grid cells are avoided: set outer interface primitives at the boundary equal to inner interface primitive at that boundary
  // This function is what's used to force primitives on flux surfaces to be chosen as desired by user rather than interpolated by code
  // This is useful (required) if want to properly and exactly define surface of an object where (e.g.) perfect conductivity and other specifications are desired.
  //
  //////////////
  remapplpr( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);


}





/// choose limiter
int choose_limiter(int dir, int i, int j, int k, int pl)
{
#if(LIMADJUST==LIMITERFIXED)
  // TESTMARK
  // if(i>=4) return(lim[dir]);
  //else return(DONOR);
  return(lim[dir]);
#else

#if(HYDROLIMADJUSTONLY)
  if(pl<B1) return(MACP0A1(pflag,i,j,k,FLAGREALLIM));
  else return(lim[dir]);
#else
  return(MACP0A1(pflag,i,j,k,FLAGREALLIM));
#endif

#endif

}









/// slope_lim() is provided p2interp and returns pleft/pright
///
/// |=interface
/// i=zone center of ith zone
///
/// |              |      p2interp(i)   |
/// |         pl(i)|pr(i)    i          |
/// |         Fl(i)|Fr(i)    i          |
/// |         Ul(i)|Ur(i)    i          |
/// |              |pleft(i)   pright(i)|
/// |              |F(i)                |
///
void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop)
{
  int pl,pliter;




  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation) slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    int loc=CENT;
    int continuous=0;
    get_loop(INTERPPOINTTYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation){
      PINTERPLOOP(pliter,pl){
        slope_lim_pointtype(ENOINTERPTYPE, realisinterp, pl, dir, loc, continuous, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}



/// slope_lim_cent2face() is provided p2interp and returns pleft/pright
/// gets interpolations in expanded region for FLUXRECON && FLUXCTSTAG method if updating quasi-deaveraged field instead of point value
///
/// |=interface
/// i=zone center of ith zone
///
/// |              |      p2interp(i)   |
/// |         pl(i)|pr(i)    i          |
/// |         Fl(i)|Fr(i)    i          |
/// |         Ul(i)|Ur(i)    i          |
/// |              |pleft(i)   pright(i)|
/// |              |F(i)                |
///
void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop)
{
  int pl,pliter;
  int interporflux;

  
  if(extrazones4emf){
    interporflux=ENOINTERPTYPE4EMF;
    
  }
  else{
    interporflux=ENOINTERPTYPE;
  }
  


  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, interporflux, dir, cent2faceloop);
    if(dointerpolation) slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, interporflux, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    int loc=CENT;
    int continuous=0;
    get_loop(INTERPPOINTTYPE, interporflux, dir, cent2faceloop);
    if(dointerpolation){
      PINTERPLOOP(pliter,pl){
        slope_lim_pointtype(interporflux, realisinterp, pl, dir, loc, continuous, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}










// local corresponding value to global STOREWAVESPEEDS
#define STOREWAVESPEEDMETHOD 1

/// compute donor-based flux, bypassing normal fluxcalc().  Allows one to check consistency.  Or if really want DONOR, then much faster.
/// To use this, remove "_donor" on function name and rename normal function to be with (e.g.) "_normal" on end of function name.
/// compute donor-based flux, bypassing normal fluxcalc().  Allows one to check consistency.  Or if really want DONOR, then much faster.
//DEBUGREPLACEFUNNAME//int fluxcalc
int fluxcalc_donor
(int stage,
 FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
 FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
 FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
 FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
 FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
 FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
 FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
 FTYPE CUf,
 FTYPE fluxdt,
 FTYPE *ndt1,
 FTYPE *ndt2,
 FTYPE *ndt3
 )
{
  int dir,i,j,k,pl,pliter;
  struct of_geom geomcdontuse;
  struct of_geom *ptrgeomc=&geomcdontuse;

  struct of_geom geomf1dontuse;
  struct of_geom *ptrgeomf1=&geomf1dontuse;
  struct of_geom geomf2dontuse;
  struct of_geom *ptrgeomf2=&geomf2dontuse;

  struct of_geom geomf1udontuse;
  struct of_geom *ptrgeomf1u=&geomf1udontuse;
  struct of_geom geomf2udontuse;
  struct of_geom *ptrgeomf2u=&geomf2udontuse;

  struct of_geom geomcorn3dontuse;
  struct of_geom *ptrgeomcorn3=&geomcorn3dontuse;
  //    FTYPE Flux[3][3][NDIM][NPR];
  //    FTYPE ctop[3][3][NDIM];
  FTYPE primcent[3][3][NPR];
  FTYPE prim[3][3][NPR];
  FTYPE uconf1[3][3][NDIM];
  FTYPE uconf2[3][3][NDIM];
  FTYPE uconcent[3][3][NDIM];
  int ii,jj,kk;
  int shifti,shiftj;
  FTYPE dtij[NDIM];
  //GLOBALMAC(F1,ii,jj,k)
  FTYPE emf2d[3][3],c2d[NUMCS][COMPDIM];
  FTYPE others[100];
  struct of_state state;
  FTYPE cminmax[3][3][NUMCS][NDIM];
  int ignorecourant=0;
  FTYPE dB[2],ctop[2],ctoporig[2];
  FTYPE ctopold[NDIM];
  FTYPE emffinal;
  struct of_state state_c, state_l, state_r;
  struct of_state *ptrstate_c, *ptrstate_l, *ptrstate_r;
  FTYPE F_c[NPR], F_l[NPR], F_r[NPR];
  FTYPE U_c[NPR], U_l[NPR], U_r[NPR];
  FTYPE ocminmax_l[NUMCS], ocminmax_r[NUMCS], ocminmax[NUMCS];
  FTYPE ocminmaxrad_l[NUMCS], ocminmaxrad_r[NUMCS], ocminmaxrad[NUMCS];
  extern   int flux_compute(int i, int j, int k, int dir, struct of_geom *geom, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE ctopmhd, FTYPE *cminmaxrad_l, FTYPE *cminmaxrad_r, FTYPE *cminmaxrad, FTYPE ctoprad, FTYPE CUf, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, FTYPE *F);
  extern int p2SFUevolve(int dir, int isleftright, FTYPE *p, struct of_geom *geom, struct of_state **ptrstate, FTYPE *F, FTYPE *U);
  FTYPE ctopother[NDIM];
  FTYPE *ctopptr=&ctopother[0];
  FTYPE ctopother2[NDIM];
  FTYPE *ctopptr2=&ctopother2[0];

  FTYPE ctopradother[NDIM];
  FTYPE *ctopradptr=&ctopradother[0];
  FTYPE ctopradother2[NDIM];
  FTYPE *ctopradptr2=&ctopradother2[0];

  FTYPE F[NPR];
  FTYPE p_l[NPR],p_r[NPR];
  FTYPE pvec_l[NDIM][NPR],pvec_r[NDIM][NPR];
  extern int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  int is,ie,js,je,ks,ke;


#if(EOMRADTYPE!=EOMRADNONE)
  dualfprintf(fail_file,"fluxcalc_donor not setup for radiation.");
  myexit(1);
#endif






  // default
  ptrstate_c = &state_c;
  ptrstate_l = &state_l;
  ptrstate_r = &state_r;

  //  COMPFULLLOOP{
  //  COMPZLOOP{
  // ODDITY: ADDED EXTRA SHIFT
  // ODDITY: ADDED EXTRA LOWER SHIFT
#define MYCOMPLOOPF3 for(k=-SHIFT3+SHIFTX3DN;k<=N3-1+SHIFT3+SHIFT3+SHIFTX3UP;k++)
#define MYCOMPLOOPF2 for(j=-SHIFT2+SHIFTX2DN;j<=N2-1+SHIFT2+SHIFT2+SHIFTX2UP;j++)
#define MYCOMPLOOPF1 for(i=-SHIFT1+SHIFTX1DN;i<=N1-1+SHIFT1+SHIFT1+SHIFTX1UP;i++)
  //#define MYCOMPLOOPF3 for(k=0+SHIFTX3DN;k<=N3-1+SHIFT3+SHIFTX3UP;k++)
  //#define MYCOMPLOOPF2 for(j=0+SHIFTX2DN;j<=N2-1+SHIFT2+SHIFTX2UP;j++)
  //#define MYCOMPLOOPF1 for(i=0+SHIFTX1DN;i<=N1-1+SHIFT1+SHIFTX1UP;i++)
  
  
  //#define MYCOMPLOOPF3 for(k=ks+SHIFTX3DN;k<=ke+SHIFTX3UP;k++)
  //#define MYCOMPLOOPF2 for(j=js+SHIFTX2DN;j<=je+SHIFTX2UP;j++)
  //#define MYCOMPLOOPF1 for(i=is+SHIFTX1DN;i<=ie+SHIFTX1UP;i++)

#define MYCOMPLOOPF MYCOMPLOOPF3 MYCOMPLOOPF2 MYCOMPLOOPF1
  
  MYCOMPLOOPF{
    // set primitives (effective interpolation)
    kk=0;
    for(ii=im1mac(i);ii<=ip1mac(i);ii++){
      for(jj=jm1mac(j);jj<=jp1mac(j);jj++){
        shifti=1+ii-i;
        shiftj=1+jj-j;

        PLOOP(pliter,pl) primcent[shifti][shiftj][pl] = MACP0A1(pr,ii,jj,kk,pl);
        PLOOP(pliter,pl) prim[shifti][shiftj][pl] = MACP0A1(pr,ii,jj,kk,pl);
        PLOOPBONLY(pl) prim[shifti][shiftj][pl] = MACP0A1(pstag,ii,jj,kk,pl); // replace with face values (should be same really)

        // if(i==390 && j==1 && k==0){
        //   dualfprintf(fail_file,"i=%d j=%d k=%d prim[B2]=%21.15g\n",ii,jj,k,prim[shifti][shiftj][B2]);
        // }

        get_geometry(ii,jj,kk, CENT, ptrgeomc);
        ucon_calc(prim[shifti][shiftj], ptrgeomc, uconcent[shifti][shiftj],others);

        // only geometrical difference compared to above
        get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
        ucon_calc(prim[shifti][shiftj], ptrgeomf1, uconf1[shifti][shiftj],others);

        // only geometrical difference compared to above
        get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
        ucon_calc(prim[shifti][shiftj], ptrgeomf2, uconf2[shifti][shiftj],others);

      }
    }



    ////////////////////////////////////////////////
    // set flux at faces (effective flux_compute)
    dir=1;

    ii=i;
    jj=j;
    kk=k;
    get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
    PLOOP(pliter,pl) p_l[pl] = prim[0][1][pl];
    PLOOP(pliter,pl) p_r[pl] = prim[1][1][pl];
    p_l[B1]=p_r[B1]=prim[1][1][B1]; // wasn't enforced at first


    get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
    get_geometry(ii,jp1mac(jj),kk, FACE2, ptrgeomf2u);
    get_geometry(ii,jj,kk, CENT, ptrgeomc);
    p_r[B2]=0.5*(prim[1][1][B2]*ptrgeomf2->gdet+prim[1][2][B2]*ptrgeomf2u->gdet)/(ptrgeomc->gdet);
    //      if(i==390 && j==1 && k==0){
    // dualfprintf(fail_file,"ONE: ii=%d jj=%d kk=%d :: %21.15g %21.15g : %21.15g %21.15g : %21.15g\n",ii,jj,kk,prim[1][1][B2],ptrgeomf2->gdet,prim[1][2][B2],ptrgeomf2u->gdet,ptrgeomc->gdet);
    // dualfprintf(fail_file,"ONEB: %d %d %d :: %d %d %d :: %d %d %d\n",ii,jj,kk,ii,jp1mac(jj),kk,ii,jj,kk);
    //      }
    get_geometry(im1mac(ii),jj,kk, FACE2, ptrgeomf2);
    get_geometry(im1mac(ii),jp1mac(jj),kk, FACE2, ptrgeomf2u);
    get_geometry(im1mac(ii),jj,kk, CENT, ptrgeomc);
    p_l[B2]=0.5*(prim[0][1][B2]*ptrgeomf2->gdet+prim[0][2][B2]*ptrgeomf2u->gdet)/(ptrgeomc->gdet);
    //      if(i==390 && j==1 && k==0){
    // dualfprintf(fail_file,"TWO: ii=%d jj=%d kk=%d :: %21.15g %21.15g : %21.15g %21.15g : %21.15g\n",ii,jj,kk,prim[0][1][B2],ptrgeomf2->gdet,prim[0][2][B2],ptrgeomf2u->gdet,ptrgeomc->gdet);
    // dualfprintf(fail_file,"ONEB: %d %d %d :: %d %d %d :: %d %d %d\n",im1mac(ii),jj,kk,im1mac(ii),jp1mac(jj),kk,im1mac(ii),jj,kk);
    //      }

    PLOOP(pliter,pl) pvec_l[dir][pl]=p_l[pl];
    PLOOP(pliter,pl) pvec_r[dir][pl]=p_r[pl];

    if(0){
      // issue with wavespeed?
      flux_compute_general(i, j, k, dir, ptrgeomf1, CUf, NULL, p_l, p_r, MAC(F1,i,j,k), &ctopptr[dir]);
    }
    else{

      MYFUN(p2SFUevolve(dir, ISLEFT, p_l, ptrgeomf1, &ptrstate_l, F_l, U_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
      MYFUN(p2SFUevolve(dir, ISRIGHT, p_r, ptrgeomf1, &ptrstate_r, F_r, U_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);
      if(STOREWAVESPEEDMETHOD==0){
        // == 0 method
        // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
        MYFUN(vchar(p_l, ptrstate_l, dir, ptrgeomf1, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        MYFUN(vchar(p_r, ptrstate_r, dir, ptrgeomf1, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr[1]);
      }
      else{
        // ==1 method
        get_geometry(im1mac(ii),jj,kk, CENT, ptrgeomc);
        get_state(primcent[0][1],ptrgeomc,&state);
        MYFUN(vchar(primcent[0][1], &state, dir, ptrgeomc, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        get_geometry(ii,jj,kk, CENT, ptrgeomc);
        get_state(primcent[1][1],ptrgeomc,&state);
        MYFUN(vchar(primcent[1][1], &state, dir, ptrgeomc, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr[1]); 
      }



      //      c2d[CMIN][0] = MAX(0.,-ocminmax[CMIN]);
      //      c2d[CMAX][0] = MAX(0.,+ocminmax[CMAX]);
      // cminmax_calc() already sets the below with positive sign and MAX'ed against zero
      c2d[CMIN][0] = ocminmax[CMIN];
      c2d[CMAX][0] = ocminmax[CMAX];
 
      // now get flux
      MYFUN(flux_compute(i, j, k, dir, ptrgeomf1, ocminmax_l,ocminmax_r, ocminmax, ctopptr[1], ocminmaxrad_l,ocminmaxrad_r, ocminmaxrad, ctopradptr[1], CUf, p_l, p_r, U_l, U_r, F_l, F_r, F),"step_ch.c:fluxcalc()", "flux_compute", 1);
      PLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl)=F[pl];

#if(0)
      // DEBUG:
      if(i==390 && j==1 && k==0){
        dualfprintf(fail_file,"\n");
        PLOOP(pliter,pl) dualfprintf(fail_file,"F1: pl=%d p_l=%21.15g p_r=%21.15g U_l=%21.15g U_r=%21.15g F_l=%21.15g F_r=%21.15g F=%21.15g\n",pl,p_l[pl],p_r[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl],F[pl]);
        dualfprintf(fail_file,"NEW: ocminmax_l[CMIN]=%21.15g ocminmax_r[CMIN]=%21.15g ocminmax_l[CMAX]=%21.15g ocminmax_r[CMAX]=%21.15g ocminmax[CMIN]=%21.15g ocminmax[CMAX]=%21.15g ctopptr[1]=%21.15g\n",ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],ocminmax[CMIN],ocminmax[CMAX],ctopptr[1]);
      }
#endif

      // get lower in j for wave speeds for EMF
      ii=i;
      jj=jm1mac(j);
      kk=k;
      get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
      PLOOP(pliter,pl) p_l[pl] = prim[0][0][pl];
      PLOOP(pliter,pl) p_r[pl] = prim[1][0][pl];
      p_l[B1]=p_r[B1]=prim[1][0][B1]; // wasn't enforced at first

      get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
      get_geometry(ii,jp1mac(jj),kk, FACE2, ptrgeomf2u);
      get_geometry(ii,jj,kk, CENT, ptrgeomc);
      p_r[B2]=0.5*(prim[1][0][B2]*ptrgeomf2->gdet+prim[1][1][B2]*ptrgeomf2u->gdet)/(ptrgeomc->gdet);

      get_geometry(im1mac(ii),jj,kk, FACE2, ptrgeomf2);
      get_geometry(im1mac(ii),jp1mac(jj),kk, FACE2, ptrgeomf2u);
      get_geometry(im1mac(ii),jj,kk, CENT, ptrgeomc);
      p_l[B2]=0.5*(prim[0][0][B2]*ptrgeomf2->gdet+prim[0][1][B2]*ptrgeomf2u->gdet)/(ptrgeomc->gdet);

      MYFUN(p2SFUevolve(dir, ISLEFT, p_l, ptrgeomf1, &ptrstate_l, F_l, U_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
      MYFUN(p2SFUevolve(dir, ISRIGHT, p_r, ptrgeomf1, &ptrstate_r, F_r, U_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);

      if(STOREWAVESPEEDMETHOD==0){
        // ==0 method
        // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
        MYFUN(vchar(p_l, ptrstate_l, dir, ptrgeomf1, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        MYFUN(vchar(p_r, ptrstate_r, dir, ptrgeomf1, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr2[1]);
        // have cmin,cmax,ctop here
      }
      else{
        // ==1 method
        get_geometry(im1mac(ii),jj,kk, CENT, ptrgeomc);
        get_state(primcent[0][0],ptrgeomc,&state);
        MYFUN(vchar(primcent[0][0], &state, dir, ptrgeomc, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        get_geometry(ii,jj,kk, CENT, ptrgeomc);
        get_state(primcent[1][0],ptrgeomc,&state);
        MYFUN(vchar(primcent[1][0], &state, dir, ptrgeomc, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr2[1]);
      }

      // cminmax_calc() already sets the below with positive sign and MAX'ed against zero
      c2d[CMIN][0] = fabs(MAX(c2d[CMIN][0],ocminmax[CMIN]));
      c2d[CMAX][0] = fabs(MAX(c2d[CMAX][0],ocminmax[CMAX]));
      ctoporig[0]  = MAX(c2d[CMIN][0],c2d[CMAX][0]);

    }


    /////////////////////////////////////////
    //
    dir=2;

    ii=i;
    jj=j;
    kk=k;
    get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
    PLOOP(pliter,pl) p_l[pl] = prim[1][0][pl];
    PLOOP(pliter,pl) p_r[pl] = prim[1][1][pl];
    p_l[B2]=p_r[B2]=prim[1][1][B2]; // wasn't enforced at first

    get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
    get_geometry(ip1mac(ii),jj,kk, FACE1, ptrgeomf1u);
    get_geometry(ii,jj,kk, CENT, ptrgeomc);
    p_r[B1]=0.5*(prim[1][1][B1]*ptrgeomf1->gdet+prim[2][1][B1]*ptrgeomf1u->gdet)/(ptrgeomc->gdet);

    get_geometry(ii,jm1mac(jj),kk, FACE1, ptrgeomf1);
    get_geometry(ip1mac(ii),jm1mac(jj),kk, FACE1, ptrgeomf1u);
    get_geometry(ii,jm1mac(jj),kk, CENT, ptrgeomc);
    p_l[B1]=0.5*(prim[1][0][B1]*ptrgeomf1->gdet+prim[2][0][B1]*ptrgeomf1u->gdet)/(ptrgeomc->gdet);

    PLOOP(pliter,pl) pvec_l[dir][pl]=p_l[pl];
    PLOOP(pliter,pl) pvec_r[dir][pl]=p_r[pl];

    if(0){
      // issue with wavespeed?
      flux_compute_general(i, j, k, dir, ptrgeomf2, CUf, NULL, p_l, p_r, MAC(F2,i,j,k), &ctopptr[dir]);
    }
    else{

      MYFUN(p2SFUevolve(dir, ISLEFT, p_l, ptrgeomf2, &ptrstate_l, F_l, U_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
      MYFUN(p2SFUevolve(dir, ISRIGHT, p_r, ptrgeomf2, &ptrstate_r, F_r, U_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);
      
      if(STOREWAVESPEEDMETHOD==0){
        // == 0 method
        // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
        MYFUN(vchar(p_l, ptrstate_l, dir, ptrgeomf2, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        MYFUN(vchar(p_r, ptrstate_r, dir, ptrgeomf2, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr[2]);
        // have cmin,cmax,ctop here
      }
      else{
        // == 1 method
        get_geometry(ii,jm1mac(jj),kk,CENT, ptrgeomc);
        get_state(primcent[1][0],ptrgeomc,&state);
        MYFUN(vchar(primcent[1][0], &state, dir, ptrgeomc, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        get_geometry(ii,jj,kk, CENT, ptrgeomc);
        get_state(primcent[1][1],ptrgeomc,&state);
        MYFUN(vchar(primcent[1][1], &state, dir, ptrgeomc, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr[2]);
      }

      // cminmax_calc() already sets the below with positive sign and MAX'ed against zero
      c2d[CMIN][1] = ocminmax[CMIN];
      c2d[CMAX][1] = ocminmax[CMAX];


      // now get flux
      MYFUN(flux_compute(i, j, k, dir, ptrgeomf2, ocminmax_l,ocminmax_r, ocminmax, ctopptr[2], ocminmaxrad_l,ocminmaxrad_r, ocminmaxrad, ctopradptr[2], CUf, p_l, p_r, U_l, U_r, F_l, F_r, F),"step_ch.c:fluxcalc()", "flux_compute", 1);
      PLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl)=F[pl];

#if(0)
      // DEBUG:
      if(i==390 && j==1 && k==0){
        dualfprintf(fail_file,"\n");
        // PLOOP(pliter,pl) dualfprintf(fail_file,"prim: pl=%d prim[1][0]=%21.15g prim[1][1]=%21.15g\n",prim[1][0][pl],prim[1][1][pl]);
        PLOOP(pliter,pl) dualfprintf(fail_file,"F2: pl=%d p_l=%21.15g p_r=%21.15g U_l=%21.15g U_r=%21.15g F_l=%21.15g F_r=%21.15g F=%21.15g\n",pl,p_l[pl],p_r[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl],F[pl]);
        dualfprintf(fail_file,"NEW: ocminmax_l[CMIN]=%21.15g ocminmax_r[CMIN]=%21.15g ocminmax_l[CMAX]=%21.15g ocminmax_r[CMAX]=%21.15g ocminmax[CMIN]=%21.15g ocminmax[CMAX]=%21.15g ctopptr[1]=%21.15g\n",ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],ocminmax[CMIN],ocminmax[CMAX],ctopptr[2]);
      }
#endif

      // get lower in i for wave speeds for EMF
      ii=im1mac(i);
      jj=j;
      kk=k;
      get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
      PLOOP(pliter,pl) p_l[pl] = prim[0][0][pl];
      PLOOP(pliter,pl) p_r[pl] = prim[0][1][pl];
      p_l[B2]=p_r[B2]=prim[0][1][B2]; // wasn't enforced at first

      get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
      get_geometry(ip1mac(ii),jj,kk, FACE1, ptrgeomf1u);
      get_geometry(ii,jj,kk, CENT, ptrgeomc);
      p_r[B1]=0.5*(prim[0][1][B1]*ptrgeomf1->gdet+prim[1][1][B1]*ptrgeomf1u->gdet)/(ptrgeomc->gdet);

      get_geometry(ii,jm1mac(jj),kk, FACE1, ptrgeomf1);
      get_geometry(ip1mac(ii),jm1mac(jj),kk, FACE1, ptrgeomf1u);
      get_geometry(ii,jm1mac(jj),kk, CENT, ptrgeomc);
      p_l[B1]=0.5*(prim[0][0][B1]*ptrgeomf1->gdet+prim[1][0][B1]*ptrgeomf1u->gdet)/(ptrgeomc->gdet);

      MYFUN(p2SFUevolve(dir, ISLEFT, p_l, ptrgeomf2, &ptrstate_l, F_l, U_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
      MYFUN(p2SFUevolve(dir, ISRIGHT, p_r, ptrgeomf2, &ptrstate_r, F_r, U_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);
      if(STOREWAVESPEEDMETHOD==0){
        // ==0 method
        // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
        MYFUN(vchar(p_l, ptrstate_l, dir, ptrgeomf2, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        MYFUN(vchar(p_r, ptrstate_r, dir, ptrgeomf2, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr2[2]);
        // have cmin,cmax,ctop here
      }
      else{
        // ==1 method
        get_geometry(ii,jm1mac(jj),kk,CENT, ptrgeomc);
        get_state(primcent[0][0],ptrgeomc,&state);
        MYFUN(vchar(primcent[0][0], &state, dir, ptrgeomc, &ocminmax_l[CMAX], &ocminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
        get_geometry(ii,jj,kk,CENT, ptrgeomc);
        get_state(primcent[0][1],ptrgeomc,&state);
        MYFUN(vchar(primcent[0][1], &state, dir, ptrgeomc, &ocminmax_r[CMAX], &ocminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
        cminmax_calc(ocminmax_l[CMIN],ocminmax_r[CMIN],ocminmax_l[CMAX],ocminmax_r[CMAX],&ocminmax[CMIN],&ocminmax[CMAX],&ctopptr2[2]);
      }

      // cminmax_calc() already sets the below with positive sign and MAX'ed against zero
      c2d[CMIN][1] = fabs(MAX(c2d[CMIN][1],ocminmax[CMIN]));
      c2d[CMAX][1] = fabs(MAX(c2d[CMAX][1],ocminmax[CMAX]));
      ctoporig[1]  = MAX(c2d[CMIN][1],c2d[CMAX][1]);

    }
    

    // zero-out as test comparison with normal code
    if(1){
      dir=1; if(i==389) PLOOP(pliter,pl) if(pl!=U1) MACP0A1(F1,i,j,k,pl)=0.0;
      dir=1; if(i==419) PLOOP(pliter,pl) if(pl!=U1) MACP0A1(F1,i,j,k,pl)=0.0;
      dir=2; if(i==0) PLOOP(pliter,pl) if(pl!=U2) MACP0A1(F2,i,j,k,pl)=0.0;
      dir=2; if(i==N2) PLOOP(pliter,pl) if(pl!=U2) MACP0A1(F2,i,j,k,pl)=0.0;
    }


    ///////////////////
    // now compute correct staggered EMF at CORN3
    get_geometry(i,j,k, CORN3, ptrgeomcorn3);


#define WHICHEMF 1

#if(WHICHEMF==0)
    // not really spatially correct for emf2d and vchar not like normal code

    // EMF3
    kk=0;
    for(ii=im1mac(i);ii<=ip1mac(i);ii++){
      for(jj=jm1mac(j);jj<=jp1mac(j);jj++){
        shifti=1+ii-i;
        shiftj=1+jj-j;
   
        // use uconf1 with B1 since that's how packed for interpolation in normal routine.  Can try alternatives.
        // emf_3 = gdet(B^2 v^1 - B^1 v^2) = F1[B2] or -F2[B1]
        emf2d[shifti][shiftj] = (
                                 //     + prim[shifti][shiftj][B1]*uconf1[shifti][shiftj][2]/uconf1[shifti][shiftj][TT]
                                 //     - prim[shifti][shiftj][B2]*uconf2[shifti][shiftj][1]/uconf2[shifti][shiftj][TT]
                                 + prim[shifti][shiftj][B2]*uconf2[shifti][shiftj][1]/uconf2[shifti][shiftj][TT]
                                 - prim[shifti][shiftj][B1]*uconf1[shifti][shiftj][2]/uconf1[shifti][shiftj][TT]
                                 );

        // GODMARK: vchar directly at CORN3, which is different than averaging procedure done in normal code
        get_state(prim[shifti][shiftj], ptrgeomcorn3, &state);
        dir=1; MYFUN(vchar(prim[shifti][shiftj], &state, dir, ptrgeomcorn3, &cminmax[shifti][shiftj][CMAX][dir], &cminmax[shifti][shiftj][CMIN][dir],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 10);
        dir=2; MYFUN(vchar(prim[shifti][shiftj], &state, dir, ptrgeomcorn3, &cminmax[shifti][shiftj][CMAX][dir], &cminmax[shifti][shiftj][CMIN][dir],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 11);


        // dualfprintf(fail_file,"i=%d j=%d k=%d ii=%d jj=%d kk=%d ",i,j,k,ii,jj,kk);
        // PLOOP(pliter,pl) dualfprintf(fail_file,"prim[%d]=%21.15g ",pl,prim[shifti][shiftj][pl]);
        // dualfprintf(fail_file,"\n");
      }
    }

    // use vchar generated at CORN3
    //c[CMIN,CMAX][0=odir1,1=odir2]
    // dir=1:
    ctop[0]    = MAX((-cminmax[0][0][CMIN][1]),(-cminmax[1][0][CMIN][1]));
    ctop[0]    = MAX(ctop[0],(-cminmax[0][1][CMIN][1]));
    ctop[0]    = MAX(ctop[0],(-cminmax[1][0][CMIN][1]));
    ctop[0]    = MAX(ctop[0],(-cminmax[1][1][CMIN][1]));
    ctop[0]    = MAX(ctop[0],0.0);
    
    ctop[0]    = MAX(ctop[0],(cminmax[0][0][CMAX][1]));
    ctop[0]    = MAX(ctop[0],(cminmax[1][0][CMAX][1]));
    ctop[0]    = MAX(ctop[0],(cminmax[0][1][CMAX][1]));
    ctop[0]    = MAX(ctop[0],(cminmax[1][1][CMAX][1]));

    // dir=2:
    ctop[1]    = MAX((-cminmax[0][0][CMIN][2]),(-cminmax[1][0][CMIN][2]));
    ctop[1]    = MAX(ctop[1],(-cminmax[0][1][CMIN][2]));
    ctop[1]    = MAX(ctop[1],(-cminmax[1][0][CMIN][2]));
    ctop[1]    = MAX(ctop[1],(-cminmax[1][1][CMIN][2]));
    ctop[1]    = MAX(ctop[1],0.0);

    ctop[1]    = MAX(ctop[1],(cminmax[0][0][CMAX][2]));
    ctop[1]    = MAX(ctop[1],(cminmax[1][0][CMAX][2]));
    ctop[1]    = MAX(ctop[1],(cminmax[0][1][CMAX][2]));
    ctop[1]    = MAX(ctop[1],(cminmax[1][1][CMAX][2]));

    // use vchar generated and "max-averaged" as in normal code
    //    ctop[0] = ctoporig[0];
    //    ctop[1] = ctoporig[1];

#else

    // 1 1
    ii=i;
    jj=j;
    kk=k;
    get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
    ucon_calc(primcent[1][1], ptrgeomf1, uconf1[1][1],others);

    get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
    ucon_calc(primcent[1][1], ptrgeomf2, uconf2[1][1],others);

    emf2d[1][1] = (
                   + prim[1][1][B2]*uconf2[1][1][1]/uconf2[1][1][TT]
                   - prim[1][1][B1]*uconf1[1][1][2]/uconf1[1][1][TT]
                   );

    // 0 1
    ii=im1mac(i);
    jj=j;
    kk=k;
    get_geometry(ii+1,jj,kk, FACE1, ptrgeomf1);
    ucon_calc(primcent[0][1], ptrgeomf1, uconf1[0][1],others);

    get_geometry(ii,jj,kk, FACE2, ptrgeomf2);
    ucon_calc(primcent[0][1], ptrgeomf2, uconf2[0][1],others);

    emf2d[0][1] = (
                   + prim[0][1][B2]*uconf2[0][1][1]/uconf2[0][1][TT]
                   - prim[1][1][B1]*uconf1[0][1][2]/uconf1[0][1][TT]
                   );

    // 1 0
    ii=i;
    jj=jm1mac(j);
    kk=k;
    get_geometry(ii,jj,kk, FACE1, ptrgeomf1);
    ucon_calc(primcent[1][0], ptrgeomf1, uconf1[1][0],others);

    get_geometry(ii,jj+1,kk, FACE2, ptrgeomf2);
    ucon_calc(primcent[1][0], ptrgeomf2, uconf2[1][0],others);

    emf2d[1][0] = (
                   + prim[1][1][B2]*uconf2[1][0][1]/uconf2[1][0][TT]
                   - prim[1][0][B1]*uconf1[1][0][2]/uconf1[1][0][TT]
                   );

    // 0 0
    ii=im1mac(i);
    jj=jm1mac(j);
    kk=k;
    get_geometry(ii+1,jj,kk, FACE1, ptrgeomf1);
    ucon_calc(primcent[0][0], ptrgeomf1, uconf1[0][0],others);

    get_geometry(ii,jj+1,kk, FACE2, ptrgeomf2);
    ucon_calc(primcent[0][0], ptrgeomf2, uconf2[0][0],others);

    emf2d[0][0] = (
                   + prim[0][1][B2]*uconf2[0][0][1]/uconf2[0][0][TT]
                   - prim[1][0][B1]*uconf1[0][0][2]/uconf1[0][0][TT]
                   );

    // don't compute cminmax[][][CMIN/CMAX][dir] since not needed

    // use vchar generated and "max-averaged" as in normal code
    ctop[0] = ctoporig[0];
    ctop[1] = ctoporig[1];
#endif



    //edgedir=3 odir1=1 odir2=2
    // dB in perp to B-dir at CORN3
    // dir=1:
    dB[0] = prim[1][1][B1]-prim[1][0][B1];
    // dir=2:
    dB[1] = prim[1][1][B2]-prim[0][1][B2];


    emffinal=(
              + 0.25*(emf2d[0][0]+emf2d[0][1]+emf2d[1][0]+emf2d[1][1])
              - 0.50*(ctop[0]*dB[1] - ctop[1]*dB[0])
              );

#if(0)
    // DEBUG:
    if(i==390 && j==1 && k==0){
      dualfprintf(fail_file,"NEW: emf2d[1][1]=%21.15g emf2d[1][0]=%21.15g emf2d[0][1]=%21.15g emf2d[0][0]=%21.15g ctop[0]=%21.15g  ctop[1]=%21.15g dB[0]=%21.15g  dB[1]=%21.15g emffinal=%21.15g gdetcorn3=%21.15g\n",emf2d[1][1],emf2d[1][0],emf2d[0][1],emf2d[0][0],ctop[0],ctop[1],dB[0],dB[1],emffinal,ptrgeomcorn3->gdet);
      dualfprintf(fail_file,"NEW: c2d[CMIN][0]=%21.15g c2d[CMAX][0]=%21.15g c2d[CMIN][1]=%21.15g c2d[CMAX][1]=%21.15g\n",c2d[CMIN][0],c2d[CMAX][0],c2d[CMIN][1],c2d[CMAX][1]);
    }
#endif
    
    emffinal *= (ptrgeomcorn3->gdet);
    
    // zero-out as test comparison with normal code
    if(1){
      if(i==389) emffinal=0.0;
      if(i==419) emffinal=0.0;
      // ODDITY: Nan grabbed at 419 with i+1 by ucum_check() [only for field]
      if(i==419+1) emffinal=0.0;
      if(j==0) emffinal=0.0;
      if(j==N2) emffinal=0.0;
      //      if(i==N2+1) emffinal=0.0;
    }

    MACP0A1(F1,i,j,k,B2) = emffinal;
    MACP0A1(F2,i,j,k,B1) = -MACP0A1(F1,i,j,k,B2);

    //    dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d j=%d k=%d :: %21.15g %21.15g\n",nstep,steppart,i,j,k,MACP0A1(F1,i,j,k,B1),MACP0A1(F2,i,j,k,B2));
    // ODDLY: required below to preserve divb=0
    // NOT ODD, diffusive term will be c*(dB1) and c*(dB2) despite zero non-diffusive flux.
    // Well, diffusive term is zero once set staggered field -> face flux correctly
    MACP0A1(F1,i,j,k,B1) = 0.0;
    MACP0A1(F2,i,j,k,B2) = 0.0;


    // set dt
    if(WITHINACTIVESECTIONEXPAND1(i,j,k)){
      //      dualfprintf(fail_file,"INSIDEDTCHECK: nstep=%ld steppart=%d i=%d j=%d k=%d\n",nstep,steppart,i,j,k);
      dir=1;
      dtij[dir] = cour * dx[dir] / ctopptr[dir]; // use wave speed from face flux calculation
      //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
      is=fluxloop[dir][FIS]+SHIFTX1DN;
      ie=fluxloop[dir][FIE]+SHIFTX1UP;
      js=fluxloop[dir][FJS]+SHIFTX2DN;
      je=fluxloop[dir][FJE]+SHIFTX2UP;
      ks=fluxloop[dir][FKS]+SHIFTX3DN;
      ke=fluxloop[dir][FKE]+SHIFTX3UP;
      if(dtij[dir]<*ndt1 && i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke){ // mimics normal code
        *ndt1=dtij[dir];
        dualfprintf(fail_file,"NEW: Got lower dir=%d i=%d j=%d k=%d ndt=%21.15g\n",dir,i,j,k,*ndt1);
      }

      dir=2;
      dtij[dir] = cour * dx[dir] / ctopptr[dir]; // use wave speed from face flux calculation
      is=fluxloop[dir][FIS]+SHIFTX1DN;
      ie=fluxloop[dir][FIE]+SHIFTX1UP;
      js=fluxloop[dir][FJS]+SHIFTX2DN;
      je=fluxloop[dir][FJE]+SHIFTX2UP;
      ks=fluxloop[dir][FKS]+SHIFTX3DN;
      ke=fluxloop[dir][FKE]+SHIFTX3UP;
      if(dtij[dir]<*ndt2 && i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke){ // mimics normal code
        *ndt2=dtij[dir];
        dualfprintf(fail_file,"NEW: Got lower dir=%d i=%d j=%d k=%d ndt=%21.15g\n",dir,i,j,k,*ndt2);
      }
    }
  }

  FTYPE (**ptrfluxvec)[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int Nvec[NDIM];
  int cleanup_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);

  fluxvec[1]=F1;
  fluxvec[2]=F2;
  fluxvec[3]=F3;

  ptrfluxvec=fluxvec;

  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  dualfprintf(fail_file,"BEFORE CLEAN\n");
  cleanup_fluxes(Nvec,ptrfluxvec);
  dualfprintf(fail_file,"AFTER CLEAN\n");



  return(0);


}
