
// OPTMARK: Should redo flux's so that fluxes are accessed by MAC(F1,j,k,i) MAC(F2,k,i,j) MAC(F3,i,j,k) for faster differencing in advance.c
// Maybe not important

#include "decs.h"



// see fluxcompute.c for non-computer science, real physics calculations of flux
int fluxcalc(int stage,
	     FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
	     FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
	     FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
	     FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
	     FTYPE (*F1)[NSTORE2][NSTORE3][NPR], 
	     FTYPE (*F2)[NSTORE2][NSTORE3][NPR], 
	     FTYPE (*F3)[NSTORE2][NSTORE3][NPR], 
 	     FTYPE CUf,
	     FTYPE fluxdt,
	     FTYPE *ndt1,
	     FTYPE *ndt2,
	     FTYPE *ndt3
	     )
{
  int fluxcalc_flux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop);
  void fix_flux(FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]) ;
  FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR];
  FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR];
  FTYPE (**ptrfluxvec)[NSTORE2][NSTORE3][NPR];
  FTYPE *ndtvec[NDIM];
  int Nvec[NDIM];
  int flux_point2avg(int stage, int whichmaorem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR]);
  void preinterp_flux_point2avg(void);
  int i,j,k;
  int pl,pliter;
  int dir;
  int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR]);
  int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR]);
  int cleanup_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR]);
  int zero_out_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR]);
  int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR]);

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
  // Wavespeeds may also be stored globally if STOREWAVESPEEDS==1
  // In all cases wavespeed constraint on timestep is set here
  //
  // assume fluxvec is MA only if splitmaem==1
  //
  ///////////////////////////////////////////////
  
  fluxcalc_flux(stage, pr, pstag, pl_ct, pr_ct, Nvec, dqvec, fluxvec, fluxvecEM, CUf, ndtvec, cent2faceloop);





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
  // FLUXCTSTAG -- overwrite field fluxes (emf's) with correct CORN1,2,3 points values
  //
  // assumes fluxcalc_flux() above called fluxcalc_standard_4fluxctstag() so pl_ct and pr_ct are set with FACE1,2,3 values of all quantities (including field along dir that comes from pstagscratch[] in this case)
  //
  /////////////////////////////
  if(FLUXB==FLUXCTSTAG){
#if(STOREWAVESPEEDS==0 || USESTOREDSPEEDSFORFLUX==0)
    dualfprintf(fail_file,"STOREWAVESPEEDS,USESTOREDSPEEDSFORFLUX must be 1 when FLUXB==FLUXCTSTAG\n");
    // must store because vchar() cannot be computed at CORN since only have partial velocities and partial fields and no densities
    // really only STOREWAVESPEEDS must be 1, but for now assume both
    myexit(175106);
#endif

    MYFUN(fluxcalc_fluxctstag(stage, pr, pstag, pl_ct, pr_ct, GLOBALPOINT(pbcorninterp), GLOBALPOINT(pvcorninterp), GLOBALPOINT(wspeed), GLOBALPOINT(prc), GLOBALPOINT(pleft), GLOBALPOINT(pright), GLOBALPOINT(fluxstatecent), GLOBALPOINT(fluxstate), GLOBALPOINT(geomcornglobal), Nvec, dqvec, ptrfluxvec, CUf, cent2faceloop, face2cornloop),"flux.c:fluxcalc()", "fluxcalc_fluxctstag", 0);

    //////////////////////////////
    //
    // User "boundary conditions" to modify EMFs before used
    //
    /////////////////////////////
    if(DOGRIDSECTIONING){
      adjust_fluxctstag_emfs(pr,Nvec,ptrfluxvec);
    }


    ////////////////
    // Before higher-order operations on flux, track vector potential update
    // so updating point value of A_i
    ////////////////
    update_vpot(stage, pr, ptrfluxvec, fluxdt, vpot);

  }

  


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

    MYFUN(flux_ct(stage, pr, GLOBALPOINT(emf), GLOBALPOINT(vconemf), GLOBALPOINT(dq1), GLOBALPOINT(dq2), GLOBALPOINT(dq3), ptrfluxvec[1], ptrfluxvec[2], ptrfluxvec[3]),"step_ch.c:advance()", "flux_ct",1);

    // Note that user modifications to EMFs done inside flux_ct() before offset Flux computed

    ////////////////
    // TOTH CT method doesn't cleanly differentiate between point update and average update of A_i, so just stick to TOTH CT EMF itself
    ////////////////
    update_vpot(stage, pr, ptrfluxvec, fluxdt, vpot);

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
#if(FLUXDUMP)
  // this accounts for final flux
  FULLLOOP{
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
	PLOOP(pliter,pl) MACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=MACP1A1(fluxvec,dir,i,j,k,pl);
      }
      else{
	PLOOP(pliter,pl) MACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=0.0;
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



// ensure fluxes only exist on well-defined computational box
// makes advance.c, GRIDSECTIONING, and adaptive time-stepping easier to code for
// Note that the branch prediction penalty for these conditions inside the LOOP is low since just zero out nearby memory elements
int cleanup_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
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



// zero-out fluxes in COMPFULLLOOP so advance.c updates all possible quantities and does so with simple non-conditional code
// it's ok to zero-out beyond standard+-1 box a bit.  Just ensure zero out enough
int zero_out_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
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

// zero-out fluxes that correspond to EMFs if evolving/tracking A_i since want boundary values to be plottable in SM or whatever
// This is  necessary since also use fluxes for other purposes in-between steps
// GODMARK: If in BCs have partial EMF and not full EMF for updating A_i, then may also cause A_i to appear funny
int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
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
int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR])
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

// sums MA and EM fluxes
int fluxsum_old(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR])
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


// sums MA (with FLUXSPLITMA(dir) containing flux[UU+dir] component) and EM fluxes
int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR])
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
	    ////	COMPFULLLOOP{
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






// wrapper for CENT_to_FACE1,2,3 used to compute flux at face
int fluxcalc_flux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvecEM[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop)
{
  int fluxcalc_flux_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata );
  int i,j,k,pl,pliter;
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


    MYFUN(fluxcalc_flux_1d(stage, pr, pstag, pl_ct, pr_ct, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dqvec[dir], fluxvec[dir], fluxvecEM[dir], CUf, ndtvec[dir], &cent2faceloop[dir], &didassigngetstatecentdata),"flux.c:fluxcalc()", "fluxcalc_flux_1d", dir);

#if(PRODUCTION==0)
    trifprintf("%d",dir);
#endif
  }// end DIMENLOOP(dir)

  return(0);

}



// wrapper for different standard 1-D flux calculators
// 1-D interpolate and get flux for that direction (assumes purely 1-D Riemann problem)
int fluxcalc_flux_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata )
{
  int fluxcalc_standard(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);
  int fluxcalc_standard_4fluxctstag(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);

  //  int fluxcalc_fluxspliteno(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE *ndt);
  int Nvec[NDIM];
  int i,j,k,pl,pliter;
  int odir1,odir2;


  
  if(DOENOFLUX!=ENOFLUXSPLIT){

    if(FLUXB==FLUXCTSTAG){
      MYFUN(fluxcalc_standard_4fluxctstag(stage,pr,pstag,pl_ct, pr_ct, dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard_4fluxctstag()", 1);
    }
    else{
      // use older code that doesn't store into pl_ct and pr_ct since not needed and then waste of memory
      MYFUN(fluxcalc_standard(stage,pr,pstag,pl_ct, pr_ct, dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard()", 1);
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




// wrapper for rescale()
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
      if(npr2interpstart<=npr2interpend) rescale(1,dir,MAC(pr,i,j,k),ptrgeom,MAC(p2interp,i,j,k));
    }// end COMPFULLLOOP
  }// end parallel region


}





// original flux calculator that gets F in "dir".  At end global pleft,pright,dq also set and if STOREWAVESPEEDS==1 then wavespeeds stored globally
int fluxcalc_standard(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
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
    p2interp=prc; // it's different
  }
  rescale_calc_full(dir,pr,p2interp);
#else
  p2interp=pr; // it's itself
#endif




#if(STOREWAVESPEEDS) // otherwise use very local estimate
#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
	MYFUN(get_global_wavespeeds_full(dir,is,ie,js,je,ks,ke,idel,jdel,kdel,POINT(pr),GLOBALPOINT(wspeed)),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
      }
#endif // end if storing wavespeeds









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
    FTYPE ctop;

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
	MYFUN(getplpr(dir,idel,jdel,kdel,i,j,k,ptrgeom,pr,pstag,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
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



      /////////////////////////////////////
      //
      // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
      //
      /////////////////////////////////////
#if(STOREFLUXSTATE)
      compute_and_store_fluxstate(dir, ISLEFT, i, j, k, ptrgeom, p_l);
      compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, ptrgeom, p_r);
      // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif




      //////////////////////////////////
      //
      // actually compute the flux
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

      // only set timestep if in computational domain or just +-1 cell beyond.  Don't go further since end up not really using that flux or rely on the stability of fluxes beyond that point.
      // Need +-1 in case flow is driven by injection boundary conditions rather than what's on grid
      if(WITHINACTIVESECTIONEXPAND1(i,j,k)){
	dtij = cour * dx[dir] / ctop;

#pragma omp critical
	{
	  if (dtij < *ndt){
	    *ndt = dtij;
	    // below are global so can report when other dt's are reported in advance.c
	    waveglobaldti[dir]=i;
	    waveglobaldtj[dir]=j;
	    waveglobaldtk[dir]=k;
	  }
	}// end critical region
      }// end if within dt-setting section

    }// end 3D loop
  }// end parallel region

  


  return (0);
}






// assign global variables needed when not doing certain dimensions
// For example, to keep algorithm simple, always assume all directions of pl_ct and pr_ct exist, but if doing 2D simulation then we don't fill "left-right" states, so instead we need to compute or assign it so that data exists
// also for the CTSTAG method, need loop information about going from center to face even if that dimension doesn't exist
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



// standard (non-field flux) calculation but setup to store results of interpolation so can be used for fluxctstag calculation
// set pl_ct and pr_ct with FACE interpolations from CENT (including field face from pstagscratch[])
// At end global pleft,pright,dq also set and if STOREWAVESPEEDS==1 then wavespeeds stored globally
int fluxcalc_standard_4fluxctstag(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*F)[NSTORE2][NSTORE3][NPR], FTYPE (*FEM)[NSTORE2][NSTORE3][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
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

#if(STOREWAVESPEEDS) // otherwise use very local estimate
#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
	MYFUN(get_global_wavespeeds_full(dir,is,ie,js,je,ks,ke,idel,jdel,kdel,POINT(pr),GLOBALPOINT(wspeed)),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
      }
#endif // end if storing wavespeeds






  //////////////////////////
  //
  // obtain pl_ct and pr_ct (point face quantities) from pr (point centered quantity)
  //
  ////////////////////////////
  interpolate_prim_cent2face(stage, realisinterp, pr, pstag, pl_ct, pr_ct, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dq, cent2faceloop);





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

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
    // generally ptr's are different inside parallel block
    ptrgeom=&geomdontuse;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      ////////////////////////
      //
      // get the geometry for the flux face
      //
      ///////////////////////
      get_geometry(i, j, k, face, ptrgeom); // OPTMARK: Seems should put back together interpolate_prim_cent2face() with flux_compute stuff below so only 1 call to get_geometry @ face....not sure why this way?!


      //////////////////////////////////
      //
      // actually compute the flux
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
      // only set timestep if in computational domain or just +-1 cell beyond.  Don't go further since end up not really using that flux or rely on the stability of fluxes beyond that point.
      // Need +-1 in case flow is driven by injection boundary conditions rather than what's on grid
      if(WITHINACTIVESECTIONEXPAND1(i,j,k)){
	dtij = cour * dx[dir] / ctop;

#pragma omp critical
	{
	  if (dtij < *ndt){
	    *ndt = dtij;
	    // below are global so can report when other dt's are reported in advance.c
	    waveglobaldti[dir]=i;
	    waveglobaldtj[dir]=j;
	    waveglobaldtk[dir]=k;
	  }
	}// end critical region
      }// end if within dt-setting section

    }// end FLUX LOOP
  }// end parallel region





  return (0);
}










// normal interpolation of CENT quantities to FACE quantities
// sets global variables pl_ct and pr_ct to p_l and p_r from interpolations
int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop)
{
  void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop);
  int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
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



  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  // assume if DOEXTRAINTERP==1, then must get here
  p2interp=prc; // it's different
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

      MYFUN(getplpr(dir,idel,jdel,kdel,i,j,k,ptrgeom,pr,pstag,p2interp,dq,GLOBALPOINT(pleft),GLOBALPOINT(pright),p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
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


      /////////////////////////////////////
      //
      // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
      //
      /////////////////////////////////////
#if(STOREFLUXSTATE)
      compute_and_store_fluxstate(dir, ISLEFT, i, j, k, ptrgeom, p_l);
      compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, ptrgeom, p_r);
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
    npr2interpstart=nprlocalstart;
    npr2interpend=nprlocalend;
    PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
  }
  

  return(0);

}



// assign one fluxstate to another fluxstate (usually used to assign left->right or visa versa when dimension doesn't exist
// Note that isleftright(input/output) can only be ISLEFT,ISRIGHT and not ISMIDDLE unless rewrote function to use fluxstatecent, but function not needed anymore
// GODMARK: no longer using this since store cent and copy it to ISLEFT,ISRIGHT
void compute_and_store_fluxstate_assign(int dimeninput, int dimenoutput, int isleftrightinput, int isleftrightoutput, int i, int j, int k)
{

  // copy entire structure
  GLOBALMACP2A0(fluxstate,dimenoutput,isleftrightoutput,i,j,k)=GLOBALMACP2A0(fluxstate,dimeninput,isleftrightinput,i,j,k);

}

// compute and store get_state() data for p_l and p_r type objects
// need to call this whenever compute left-right primitives used for computing flux
void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr)
{
  int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);

  // get state
  pureget_stateforfluxcalc(pr,geom,&GLOBALMACP2A0(fluxstate,dimen,isleftright,i,j,k));

}

// compute and store get_state() data for centered state
void compute_and_store_fluxstatecent(FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  int pureget_stateforfluxcalcorsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int is,ie,js,je,ks,ke,di,dj,dk;





  // define +-1 in every direction loop range
  set_interppoint_loop_ranges_3Dextended(ENOINTERPTYPE, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);



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

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
    

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

 

      // set geometry for centered zone to be updated
      get_geometry(i, j, k, CENT, ptrgeom);

      // get state
      pureget_stateforfluxcalcorsource(MAC(pr,i,j,k),ptrgeom,&GLOBALMAC(fluxstatecent,i,j,k));

      // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen][ISLEFT,ISRIGHT]
      // If doing staggered field method, then also assumes left,right stored like pl_ct and pr_ct are stored
      // note we are copying entire structure here
#if(N1==1 || FIELDSTAGMEM)
      // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen][ISLEFT,ISRIGHT]
      GLOBALMACP2A0(fluxstate,1,ISLEFT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
      GLOBALMACP2A0(fluxstate,1,ISRIGHT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
#endif
#if(N2==1 || FIELDSTAGMEM)
      // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen,ISLEFT,ISRIGHT]
      GLOBALMACP2A0(fluxstate,2,ISLEFT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
      GLOBALMACP2A0(fluxstate,2,ISRIGHT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
#endif
#if(N3==1 || FIELDSTAGMEM)
      // if dimension doesn't exist, then copy over fluxstatecent to fluxstate[dimen,ISLEFT,ISRIGHT]
      GLOBALMACP2A0(fluxstate,3,ISLEFT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
      GLOBALMACP2A0(fluxstate,3,ISRIGHT,i,j,k)=GLOBALMAC(fluxstatecent,i,j,k);
#endif

    

    }// end 3D loop
  }// end parallel region  

}








// use dq,pleft,pright to obtain p_l and p_r for CENT to FACE
int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *ptrgeom, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r)
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
    rescale(-1,dir,p_l,ptrgeom,p2interp_l);
    rescale(-1,dir,p_r,ptrgeom,p2interp_r);
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
  set_plpr(dir,i,j,k,pr,p_l,p_r);
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





// check that left and right primitives are valid
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
  MYFUN(limit_gamma(GAMMAMAX,p_l,NULL,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 1);  //jon corr, see email re: SPINT warnings from 4/24/2006 10:54 p.m.
  MYFUN(limit_gamma(GAMMAMAX,p_r,NULL,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 2);  //jon corr
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







//////////////////////////////////////
//
// interpolate primitive using slope or just copy pleft/pright into the correct place
//
// |=interface
// i=zone center of ith zone
//
// |              |     p2interp(i)    |
// |              |       dq(i)        |
// |        p_l(i)|p_r(i)   i          |
// |              |pleft(i)   pright(i)|
//
//
//
//////////////////////////////////////

// GODMARK: as Sasha mentions, for shifting stencil shouldn't just extrapolte value from non-local values, but use other slope and most LOCAL value to obtain extrapolation.

void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r)
{
  FTYPE Xcent[NDIM],Xleft[NDIM];
  FTYPE Vcent[NDIM],Vleft[NDIM];
  FTYPE rleft,rcent,thleft,thcent;
  int pl,pliter;
  int locallim;
  int choose_limiter(int dir, int i, int j, int k, int pl);

  //reshuffle dq's such that reconstruction within active zones does not use the boundary values
  remapdq( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);
  
  if((HORIZONSUPERFAST)&&((lim[dir]<PARA)&&(LIMADJUST==0))&&(dir==1)){
    // since this uses dq's to shift, must always have dq's everywhere.  Hence LIMADJUST==0 and lim<PARA must be true
    bl_coord_ijk_2(im1mac(i), j, k, CENT, Xleft,Vleft);
    bl_coord_ijk_2(i, j, k, CENT, Xcent,Vcent);
    rleft=Vleft[1];
    thleft=Vleft[2];
    rcent=Vcent[1];
    thcent=Vcent[2];

    PINTERPLOOP(pliter,pl){
      locallim=choose_limiter(dir, i,j,k,pl);
      // get interpolated quantity
      if((locallim<PARA)&&(LIMADJUST==0)){
	if(rleft>Rhor) p2interp_l[pl] = MACP0A1(p2interp,i - idel,j - jdel,k - kdel,pl) + 0.5 * MACP0A1(dq,i - idel,j - jdel,k - kdel,pl);
	else p2interp_l[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP0A1(dq,i,j,k,pl);
	if(rcent>Rhor) p2interp_r[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP0A1(dq,i,j,k,pl);
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
      // get interpolated quantity
      if((locallim<PARA)&&(LIMADJUST==0)){
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

  //for boundaries where grid cells are avoided: set outer interface primitives at the boundary equal to inner interface primitive at that boundary
  remapplpr( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);
}





// choose limiter
int choose_limiter(int dir, int i, int j, int k, int pl)
{
#if(LIMADJUST==LIMITERFIXED)
  // TESTMARK
  // if(i>=4)	return(lim[dir]);
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









// slope_lim() is provided p2interp and returns pleft/pright
//
// |=interface
// i=zone center of ith zone
//
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//
void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *cent2faceloop)
{
  int pl,pliter;




  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation) slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    get_loop(INTERPPOINTTYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation){
      PINTERPLOOP(pliter,pl){
	slope_lim_pointtype(ENOINTERPTYPE, realisinterp, pl, dir, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}



// slope_lim_cent2face() is provided p2interp and returns pleft/pright
// gets interpolations in expanded region for FLUXRECON && FLUXCTSTAG method if updating quasi-deaveraged field instead of point value
//
// |=interface
// i=zone center of ith zone
//
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//
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
    get_loop(INTERPPOINTTYPE, interporflux, dir, cent2faceloop);
    if(dointerpolation){
      PINTERPLOOP(pliter,pl){
	slope_lim_pointtype(interporflux, realisinterp, pl, dir, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}












