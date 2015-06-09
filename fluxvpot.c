#include "decs.h"

/*! \file fluxvpot.c
  \brief All things related to vector potential: A_mu
 */

/// get directions for cyclic variables
/// signflux determines signature of relationship between flux and A_i and likewise the EMF_i = -\detg E_i
/// As in fluxct.c:
/// For evolution:
/// d_t A1 = EMF1 = -\detg E1 = F2[B3] = -F3[B2]
/// d_t A2 = EMF2 = -\detg E2 = F3[B1] = -F1[B3]
/// d_t A3 = EMF3 = -\detg E3 = F1[B2] = -F2[B1]
/// 
/// For initialization:
/// A1 = F2[B3] = -F3[B2]
/// A2 = F3[B1] = -F1[B3]
/// A3 = F1[B2] = -F2[B1]
/// \detg B1 = d_2 A3 - d_3 A2
/// \detg B2 = d_3 A1 - d_1 A3
/// \detg B3 = d_1 A2 - d_2 A1
/// 
/// opposite ordering required to be used when F[odir1] doesn't exist because N[odir1]==1
/// Should use signflux when changing between A_i <-> flux
int get_fluxpldirs(int *Nvec, int dir, int *fluxdir, int* pldir, int *plforflux, FTYPE *signflux)
{
  int odir1,odir2;


  // get cyclic other dimensions
  get_odirs(dir,&odir1,&odir2);

  if(Nvec[odir1]>1)     { *fluxdir=odir1; *pldir=odir2; *plforflux = B1+ *pldir-1; *signflux=1.0; }   //  A[dir] = F[odir1][B1+odir2-1] // normal ordering
  else if(Nvec[odir2]>1){ *fluxdir=odir2; *pldir=odir1; *plforflux = B1+ *pldir-1; *signflux=-1.0; }  //  A[dir] = F[odir2][B1+odir1-1] // opposite ordering
  else{ *fluxdir=0; *pldir=0; *plforflux=0; }
  //  else{
  //    dualfprintf(fail_file,"Shouldn't reach EMF line integration with Nvec=%d %d %d\n",Nvec[1],Nvec[2],Nvec[3]);
  //    myexit(917515);
  //  }

  return(0);
}


/// cyclic other dimensions
/// dir=1 odir1=2 odir2=3 so signflux=1 for this ordering
/// dir=2 odir1=3 odir2=1 so signflux=1 for this ordering
/// dir=3 odir1=1 odir2=2 so signflux=1 for this ordering
void get_odirs(int dir,int *odir1,int *odir2)
{
  // m%3+1 gives next 1->2,2->3,3->1
  // 3-(4-m)%3 = (dir+1)%3+1 gives previous 1->3,2->1,3->2
  *odir1=dir%3+1; // if dir==3, then odir1=1
  *odir2=(dir+1)%3+1; // if dir==3, then odir2=2
}


/// set locations for vpot or emf and number of independent memory locations for each A_i or E_i to consider
/// GODMARK: Notice the positions loc[][] are only for initializing field, while FLUXCTTOTH evolves A_i really at FACE
int set_location_fluxasemforvpot(int dir, int *numdirs, int *odir1, int *odir2, int *loc)
{
  int fieldloc[NDIM];

  get_numdirs_fluxasemforvpot(numdirs,fieldloc); // field loc has locations of each field component (not used here)


  get_odirs(dir,odir1,odir2);


  if(FLUXB==FLUXCTHLL){
    // for this method we use F1/F2/F3 and locations are different

    loc[*odir1]=FACE1-1 + *odir1;
    loc[*odir2]=FACE1-1 + *odir2;

  }
  else{

    loc[*odir1]=CORN1-1 + dir;
    loc[*odir2]=CORN1-1 + dir;

  }

  return(0);



}

/// set locations for vpot or emf and number of independent memory locations for each A_i or E_i to consider
int get_numdirs_fluxasemforvpot(int *numdirs, int *fieldloc)
{

  ///////////
  //
  // These choices must be made consistent with use of vpot2field()
  //
  ///////////


  if(! (FLUXB==FLUXCTHLL || FLUXB==FLUXCTSTAG) ){
    *numdirs = 1;
    fieldloc[1]=CENT;
    fieldloc[2]=CENT;
    fieldloc[3]=CENT;
  }
  else if(FLUXB==FLUXCTHLL){
    // for this method we use F1/F2/F3 and locations are different
    *numdirs=2;
    fieldloc[1]=CENT;
    fieldloc[2]=CENT;
    fieldloc[3]=CENT;
  }
  else{

    if(extrazones4emf==0){

      if(DOENOFLUX==ENOFLUXRECON && dofluxreconevolvepointfield==1){
        // for this method we use F1/F2/F3 but location is same corner
        *numdirs=2;
      }
      else{
        // unless emfextrazones4emf==0 and fluxrecon and doing point field, then don't need both directions
        *numdirs=1;
      }
    }
    else{
      // then just using A not F1/F2/F3
      *numdirs=1;
    }

    if(FLUXB==FLUXCTSTAG){
      fieldloc[1]=FACE1;
      fieldloc[2]=FACE2;
      fieldloc[3]=FACE3;
    }
    else{
      // non-HLL version
      fieldloc[1]=CENT;
      fieldloc[2]=CENT;
      fieldloc[3]=CENT;
    }

  }

  return(0);



}



/// assumes normal field p
/// pleft used as temp var if FLUXB==FLUXCTSTAG
/// assigns conserved field in UEVOLVE form (i.e. with gdet always)
/// implicitly Flux F1,F2,F3 are inputted into function
int vpot2field(SFTYPE time, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;
  int numdirs, fieldloc[NDIM];


  ////////////////
  //
  // get whether using A or F1/F2/F3 and where final field located per direction of field
  //
  ////////////////
  get_numdirs_fluxasemforvpot(&numdirs,fieldloc);


  ////////////////
  //
  // One of 3 things:
  // 1) single line-integrate vector potential (A_i) for FV method
  // 2) double transverse quasi-deaverage (per A_i) for FLUXRECON method if evolving Bhat
  // 3) single quasi-deaverage (per flux direction) for FLUXRECON if evolving point field
  //
  /////////////////
  if(DOENOFLUX==NOENOFLUX){
    // nothing to do
  }
  else if(DOENOFLUX==ENOFLUXRECON || DOENOFLUX == ENOFINITEVOLUME){

    if(numdirs==2){
      // difference between CTHLL and CTSTAG is position of flux
      // then treat flux itself as flux and assume flux set as vector potential
      // don't actually use A here
      // if doing FLUXRECON and evolving points, then 
      vectorpot_useflux(STAGEM1, pfield, F1, F2, F3);
    }
    else{
      // SUPERGODMARK: should tell where field is located
      vectorpot_fluxreconorfvavg(STAGEM1, pfield, A, F1, F2, F3);
    }
  }
  else{
    dualfprintf(fail_file,"No vectorpot for this case of DOENOFLUX=%d\n",DOENOFLUX);
    myexit(83465765);
  }




  ////////////////
  //
  // Now use higher-order vector potential to obtain (primitive AND conserved) field that satisfies divb=0 in certain form
  //
  /////////////////

  if(numdirs==1){
    // use of A
    if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)){
      vpot2field_centeredfield(A,pfield,ucons);
    }
    else if(FLUXB==FLUXCD){
      vpot2field_centeredfield(A,pfield,ucons); // probably wrong for FLUXCD? GODMARK
      dualfprintf(fail_file,"FLUXCD probably not setup right in vpot2field()\n");
      myexit(196726);
    }
    else if(FLUXB==FLUXCTSTAG){
      // overwrites any old choice for unew for field only
      // pstag is in general not a point value of field!
      vpot2field_staggeredfield(A,pfield,ucons);
    }
    else{
      dualfprintf(fail_file,"Undefined numdirs==1 with FLUXB==%d\n",FLUXB);
      myexit(2752632);
    }
  }
  else if(numdirs==2){

    // use of F1/F2/F3
    vpot2field_useflux(fieldloc,pfield,ucons,F1,F2,F3);

    // FLUXCTHLL is just for testing how behaves without divb=0 control
    // with fluxcthll always numdirs==2 so uses flux instead of A
    
    // overwrites any old choice for unew for field only
    // pstag is in general not a point value of field!
    if(! (FLUXB==FLUXCTHLL || FLUXB==FLUXCTSTAG) ){
      dualfprintf(fail_file,"Undefined numdirs==2 with FLUXB==%d\n",FLUXB);
      myexit(23578234);
    }

    if(FLUXB==FLUXCTSTAG && extrazones4emf==0 && DOENOFLUX == ENOFLUXRECON){
      // DEBUG:
      //      bound_pstag(STAGEM1,t,pfield, pstag, ucons, 1, USEMPI);

      // then need to obtian Bhat so can compute divb
      // In this case unew is inputted point conserved field
      field_Bhat_fluxrecon(pfield,ucons,Bhat);
    }
    
  }


  ///////////////
  //
  // copy result of vector potential calculation to staggered field if doing FLUXCTSTAG
  //
  ///////////////
  if(FLUXB==FLUXCTSTAG){
    copy_3d_fieldonly_fullloop(pfield,pstag);
    // from now on, pfield is assumed to be at CENT
  }




  ////////////////
  //
  // Before ucons->{upoint,ppoint} that involves interpolation, need to make sure setup pre-interpolates if haven't already.
  // Uses pglobal like initbase.c consistent with first calculation of pre_interpolation
  //
  /////////////////

  if(STORESHOCKINDICATOR==1 && global_preinterpolate==0){
    pre_interpolate_and_advance(GLOBALPOINT(pglobal));
  }

  ////////////////
  //
  // convert average or quasi-deaveraged field to point field (ulast and pfield)
  // Ok to use ulastglobal as temporary space
  //
  /////////////////

  ucons2upointppoint(time,pfield,pstag,ucons,uconstemp,pfield);



  // Since above procedures changed pfield that is probably pcent that is p, we need to rebound p since pfield was reset to undefined values in ghost cells since A_i isn't determined everywhere
  // alternatively for evolve_withvpot() could have inputted not the true p or some copy of it so wouldn't have to bound (except up to machine error difference when recomputed field using A_i)
  int finalstep=1; // assume user wants to know that conserved quants changed
  bound_prim(STAGEM1,finalstep, time,pfield,pstag,ucons, USEMPI);


  return(0);

}






/// convert conservative U to stag point P and CENT point U and CENT point primitive
int ucons2upointppoint(SFTYPE boundtime, FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*ulast)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;
  int numdirs, fieldloc[NDIM];


  ////////////////
  //
  // get whether using A or F1/F2/F3 and where final field located per direction of field
  //
  ////////////////
  get_numdirs_fluxasemforvpot(&numdirs,fieldloc);



  // SUPERDEBUG:
  //bound_pstag(STAGEM1, t, pfield, pfield, unew);
  //bound_pstag(STAGEM1, t, pfield, pfield, unew);
  //    dualfprintf(fail_file,"DEBUG Got here\n");


  // if restarting, then unew was read-in and then later assigned to unitial, so by the time init is done we have staggered fields in unew/unitial for any case
  // from staggered field still need to fill pfield with centered version
  // for first call to interpolation need real primitive, but don't yet have CENT field, so use staggered version as "guess"-- GODMARK not strictly grid-correct



  ////////////////
  //
  // convert average or quasi-deaveraged field to point field (ulast)
  //
  /////////////////

  if(numdirs==1){
    if(DOENOFLUX == ENOFINITEVOLUME){
      // unew and ulast both inputted so function compares results with original
      deaverage_fields_fv(pfield,unew,ulast);
    }
    else if(DOENOFLUX == ENOFLUXRECON){
      // unew and ulast both inputted so function compares results with original
      field_integrate_fluxrecon(STAGEM1,pfield,unew,ulast);
    }
    else{
      // second order methods have same average and point value
      copy_3d_fieldonly_fullloop(unew,ulast);
    }
  }
  else{
    // FLUXB==FLUXCTHLL has point value as conserved field value like ENOFLUXRECON has for non-field values
    copy_3d_fieldonly_fullloop(unew,ulast);
  }


  ////////////////
  //
  // Convert point staggered field to point CENT field for staggered field method
  //
  /////////////////

  if(FLUXB==FLUXCTSTAG){

    // ok to insert pfield for real primitive since any update to pfield is ok to be used since only used for shock indicator
    interpolate_ustag2fieldcent(STAGEM1,boundtime, -1,-1,pfield,pstag,ulast,pcent);
    // now pstagscratch will have correct staggered point value and ulast (if needed) will be replaced with center conserved value and pfield will contain centered field primitive value

    // now field sits in both centered and staggered positions with initial data

  }


  return(0);
}



/// initialize vector potential given user function
/// assumes normal non-staggered field in pr
/// Notice that if using F (flux), then *location* can be different for (e.g.) F1[B2] and F2[B1] while if using A_3 then no choice in varying positions
int init_vpot(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];


#if(TRACKVPOT)
  A=vpot; // real t=0 value put here and used to track A_i for diagnostics
#else
  A=Atemp; // dummy memory space not used till computation so safe.
#endif



  // prim -> A using user function init_vpot_user()
  init_vpot_justAcov(t, prim, A);  // t is ok since always t=tstart

  // A->Flux if required
  init_vpot_toF(A, F1, F2, F3);

  // assigns value to p and pstagscratch and unew (uses F1/F2/F3 if doing FLUXCTHLL)
  // often user just uses init_vpot2field() unless not always using vector potential
  init_vpot2field_user(t, A,prim,pstag,ucons,Bhat);  // t is ok since always t=tstart

  // in case user perturbs vpot, need to ensure ghost cells updated.
  if(1){
    int boundvartype=BOUNDVPOTTYPE;
    int finalstep=0; // assume user wants to know if initial conserved quants changed
    int doboundmpi=1;
    int doboundnonmpi=0;
    bound_vpot(STAGEM1, finalstep, t, boundvartype, vpot, doboundmpi, doboundnonmpi);
  }
  


  return(0);

}



/// get A_i in PRIMECOORDS for FLUXB==FLUXCTSTAG at CORNi for each A_i, the natural staggered field locations for A_i (i.e. not all same location for all A_i)
int get_vpot_fluxctstag_primecoords(SFTYPE time, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *vpot)
{
  int dir;
  FTYPE X[NDIM],V[NDIM];
  int loc;
  FTYPE vpotuser[NDIM];
  int ii,jj,kk;
  int whichcoord;
  int userdir;


  // copied from init_vpot_justAcov()
  // loop over A_l's
  DIMENLOOP(dir){
    
    if(dir==1) loc=CORN1;
    else if(dir==2) loc=CORN2;
    else if(dir==3) loc=CORN3;
    
    
    bl_coord_ijk_2(i, j, k, loc, X, V); // face[odir2] and flux[odir2]
    // dxdxprim_ijk(i, j, k, loc, dxdxp);
    
    // get user vpot in user coordinates (assume same coordinates for all A_{userdir}
    DLOOPA(userdir){
      init_vpot_user(&whichcoord, userdir, time, i,j,k, loc, prim, V, &vpotuser[userdir]);
    }
    
    // convert from user coordinate to PRIMECOORDS
    ucov_whichcoord2primecoords(whichcoord, i, j, k, loc, vpotuser);
    
    // now copy the only component required
    vpot[dir]=vpotuser[dir];
    
  }// loop over A_i


  return(0);

}



/// initialize vector potential given user function
/// assumes normal non-staggered field in pr
int init_vpot_justAcov(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int Nvec[NDIM];
  int odir1[NDIM],odir2[NDIM];
  int numdirs;
  int locvpot[NDIM][NDIM];
  int dir;



  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;



  // get location of vpot and number of locations to set
  DIMENLOOP(dir){
    set_location_fluxasemforvpot(dir, &numdirs, &odir1[dir], &odir2[dir], locvpot[dir]);

    // these quantities are associated with a single choice for relationship between A_i and F_j[B_k]
    //    get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir], &signforfluxlist[dir]);
  }


  if(numdirs==2){
    // will be done in another function
  }
  else{
    trifprintf("Initialize field from vector potential (really A)\n");
    //    DIMENLOOP(dir) COMPFULLLOOPP1 MACP1A0(A,dir,i,j,k) = 0.;
    DIMENLOOP(dir) init_3dvpot_fullloopp1(0.0,A[dir]);
#if(TRACKVPOT)
    init_3dvpot_fullloopp1(0.0,A[TT]);
#endif
  }
    

  trifprintf("Initialize field from vector potential assign: init_vpot_justAcov()\n");



  // GODMARK: Caution: Possible to use quantity off grid
  // (e.g. density) to define lower corner value of A, which then
  // defines B at center for lower cells.
  // Do have *grid* quantities for everywhre A is.
  
  ////////////  COMPFULLLOOPP1{
#pragma omp parallel private(dir) OPENMPGLOBALPRIVATEFULL // user could do anything
  {
    FTYPE X[NDIM],V[NDIM];
    FTYPE dxdxp[NDIM][NDIM];
    int otherdir;
    FTYPE vpotuser[NDIM];
    int userdir;
    int whichcoord;
    int loc;
    int i,j,k,pl,pliter;
    //  int fluxdir,pldir,plforflux;
    //  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
    //  FTYPE signforfluxlist[NDIM];
    //  FTYPE signforflux;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULLP1;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

    
      // can't define F1/F2/F3 in outermost part ofCOMPFULLLOOPP1 since don't exist there
      // can only define as if doing COMPFULLLOOP
      if(numdirs==2 && (i>OUTFULL1 || j>OUTFULL2 || k>OUTFULL3) ) continue;

      // loop over A_l's
      DIMENLOOP(dir){

        // these quantities are associated with a single choice for relationship between A_i and F_j[B_k]
        //      fluxdir=fluxdirlist[dir];
        //      pldir=pldirlist[dir];
        //      plforflux=plforfluxlist[dir];
        //      signforflux=signforfluxlist[dir];


        // loop over positions to get vector potential
        for(otherdir=1;otherdir<=numdirs;otherdir++){
 
          if(otherdir==1){
            loc=locvpot[dir][odir1[dir]];
          }
          else if(otherdir==2){
            loc=locvpot[dir][odir2[dir]];
          }

          bl_coord_ijk_2(i, j, k, loc, X, V); // face[odir2] and flux[odir2]
          // dxdxprim_ijk(i, j, k, loc, dxdxp);


          // get user vpot in user coordinates (assume same coordinates for all A_{userdir}
          DLOOPA(userdir){
            init_vpot_user(&whichcoord, userdir, time, i,j,k, loc, prim, V, &vpotuser[userdir]);
          }

          // convert from user coordinate to PRIMECOORDS
          ucov_whichcoord2primecoords(whichcoord, i, j, k, loc, vpotuser);

          // for numdirs=anything (used for input for 2Flux and for 2Flux required if doing TRACKVPOT
          MACP1A0(A,dir,i,j,k) = vpotuser[dir];

        }// end for loop over otherdirs
      }// end over A_i


    }// end over i,j,k

  }// end parallel region





  return(0);

}



/// initialize vector potential given user function
/// assumes normal non-staggered field in pr
/// Notice that if using F (flux), then *location* can be different for (e.g.) F1[B2] and F2[B1] while if using A_3 then no choice in varying positions
int init_vpot_toF(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int Nvec[NDIM];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int odir1[NDIM],odir2[NDIM];
  int numdirs;
  int locvpot[NDIM][NDIM];
  int dir;



  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  fluxvec[1] = F1;
  fluxvec[2] = F2;
  fluxvec[3] = F3;


  // get location of vpot and number of locations to set
  DIMENLOOP(dir){
    set_location_fluxasemforvpot(dir, &numdirs, &odir1[dir], &odir2[dir], locvpot[dir]);

    // these quantities are associated with a single choice for relationship between A_i and F_j[B_k]
    //    get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir], &signforfluxlist[dir]);
  }


  if(numdirs==2){
    // then fill F1/F2/F3 instead of A[]
    trifprintf("Initialize field from vector potential (really flux)\n");
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
        init_3dnpr_fullloop_flux(0.0,fluxvec[dir]);
      }
    }

  }
  else{
    // then A already initialized
  }
    

  trifprintf("Initialize field from vector potential assign: init_vpot_toF()\n");



  // GODMARK: Caution: Possible to use quantity off grid
  // (e.g. density) to define lower corner value of A, which then
  // defines B at center for lower cells.
  // Do have *grid* quantities for everywhre A is.
  
  ////////////  COMPFULLLOOPP1{
#pragma omp parallel private(dir) OPENMPGLOBALPRIVATEFULL // user could do anything
  {
    int otherdir;
    int userdir;
    int whichcoord;
    int loc;
    int i,j,k,pl,pliter;
    //  int fluxdir,pldir,plforflux;
    //  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
    //  FTYPE signforfluxlist[NDIM];
    //  FTYPE signforflux;

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULLP1;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

    
      // can't define F1/F2/F3 in outermost part ofCOMPFULLLOOPP1 since don't exist there
      // can only define as if doing COMPFULLLOOP
      if(numdirs==2 && (i>OUTFULL1 || j>OUTFULL2 || k>OUTFULL3) ) continue;

      // loop over A_l's
      DIMENLOOP(dir){

        // these quantities are associated with a single choice for relationship between A_i and F_j[B_k]
        //      fluxdir=fluxdirlist[dir];
        //      pldir=pldirlist[dir];
        //      plforflux=plforfluxlist[dir];
        //      signforflux=signforfluxlist[dir];


        // loop over positions to get vector potential
        for(otherdir=1;otherdir<=numdirs;otherdir++){
 
          if(otherdir==1){
            loc=locvpot[dir][odir1[dir]];
          }
          else if(otherdir==2){
            loc=locvpot[dir][odir2[dir]];
          }


          // normal ordering is A[dir] = fluxvec[fluxdir][plforflux] with fluxdir=odir1 and plforflux from odir2
 
          if(numdirs==2){
            // then using F1/F2/F3
            // Notice that fluxvec here has no \detg in it, as consistent with taking B = curlA ($\detg B^i = d_j A_k \epsilon^{ijk}$) when used
            // Notice that both fluxes are assigned a value in some cases, as required
            // e.g. F1[B2]=A3 and F2[B1]=-A3
            if(otherdir==1 && Nvec[odir1[dir]]>1) MACP1A1(fluxvec,odir1[dir],i,j,k,B1-1+odir2[dir])=MACP1A0(A,dir,i,j,k);
            if(otherdir==2 && Nvec[odir2[dir]]>1) MACP1A1(fluxvec,odir2[dir],i,j,k,B1-1+odir1[dir])=-MACP1A0(A,dir,i,j,k); // opposite ordering


          }
          else{
            // then already assigned to A
          }
        }// end for loop over otherdirs
      }// end over A_i


    }// end over i,j,k

  }// end parallel region





  return(0);

}












/// copy A_i to F_i[B1..B3] -- for use with some field evolution methods during time evolution
/// Note can't just use this in init_vpot() since there F's normally associated with a single A_i might have different positions (e.g. FLUXCTHLL method)
int copy_vpot2flux(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int Nvec[NDIM];
  int odir1[NDIM],odir2[NDIM];
  int numdirs;
  int locvpot[NDIM][NDIM];
  int dir;




  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  fluxvec[1] = F1;
  fluxvec[2] = F2;
  fluxvec[3] = F3;


  // get location of vpot and number of locations to set
  DIMENLOOP(dir){
    set_location_fluxasemforvpot(dir, &numdirs, &odir1[dir], &odir2[dir], locvpot[dir]);
  }


  if(numdirs==2){

    // 0-out flux (F1/F2/F3)
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
        init_3dnpr_fullloop_flux(0.0,fluxvec[dir]);
      }
    }



    // can't define F1/F2/F3 in outermost part of COMPFULLLOOPP1 since don't exist there
    // can only define as if doing FULLLOOP
    //////    COMPFULLLOOP{
#pragma omp parallel private(dir) // no copyin since no use of non-array globals (yet)
    {
      int otherdir;
      int i,j,k;
      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        DIMENLOOP(dir){

          // loop over positions to get vector potential
          for(otherdir=1;otherdir<=numdirs;otherdir++){
 
            // then using F1/F2/F3
            // Notice that fluxvec here has no \detg in it, as consistent with taking B = curlA ($\detg B^i = d_j A_k \epsilon^{ijk}$) when used
            // Notice that both fluxes are assigned a value in some cases, as required
            if(otherdir==1 && Nvec[odir1[dir]]>1) MACP1A1(fluxvec,odir1[dir],i,j,k,B1-1+odir2[dir])=MACP1A0(A,dir,i,j,k);
            if(otherdir==2 && Nvec[odir2[dir]]>1) MACP1A1(fluxvec,odir2[dir],i,j,k,B1-1+odir1[dir])=-MACP1A0(A,dir,i,j,k); // opposite ordering

          }// end for loop over otherdirs
        }// end over A_i


      }// end over i,j,k
    }//end if parallel region
  }// end if numdirs==2 (and so actually need to copy vpot to flux)

  return(0);

}





/// once this function is done, ensure obtain higher-order with test=1102, etc.
/// GODMARK: what about fields in direction with no dimension, like B3 in 2D?  Can't obtain B3 from A_i??
///          Seems fine B3=A2,1 +-A1,2
///          What about B1,B2?  Only needs A_3 since don't need B1=A_2,3 and B2=A_1,3 ?  Can I assume those are zero?
/// In 1D only have 0=A1,1 B3=A2,1 and B2=A3,1, so need to be able to obtain B1 somehow -- or just don't change B1 since must be constant in time
/// Note if FLUXCTHLL or FLUXCTTOTH, then no single A_i is evolved, so can't evolve with A_i in that case since field is determined by 2 different versions of A_i that we don't track (nor should!)
int evolve_withvpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{

  if(EVOLVEWITHVPOT==0){
    


  }
  else{

    // 1) fill F or vpot(A,prim);
    // 2) call vpot2field()

    // perhaps split init_vpot so can be used for this function too
    // **TODO** essentially need a separate function that takes point A_i and fills F if needed rather than all at once as in init_vpot()

    // only copies if necessary for how method is setup (e.g. staggered point evolution method)
    // i.e. F is not really flux here, just place-holder for A_i at different positions as required by some other methods, such as FLUXCTTOTH or FLUXCTSTAG w/ higher-order corrections.
    copy_vpot2flux(vpot,F1,F2,F3);


    // Note that vpot2field bounds pstag and input prim as required since A_i doesn't exist everywhere
    // GODMARK: use of many globals: ok for now since evolve_withvpot() not used for any other purpose
    vpot2field(time, vpot,prim, pstag, unew, Bhat,F1,F2,F3,Atemp,uconstemp);



#if(0)
    // DEBUG: Trying to get test=1102 to work with EVOLVEWITHVPOT (also turned on bound_flux_fluxrecon() in update_vpot() so E_i fully periodic in BZs before A_i iterated
    bound_pstag(STAGEM1,t,prim, pstag, unew);
    bound_prim(STAGEM1,t,prim, pstag, unew);
    bound_uavg(STAGEM1,t,prim, pstag, unew);
#endif
  }

  return(0);

}





/// find divb
/// if higher order method, then must use conserved value U[] assumed to then exist and be used for field
void setfdivb(FTYPE *divb, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], int i, int j, int k)
{



  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)){
    if(DOENOFLUX==NOENOFLUX){ SETFDIVBFLUXCTTOTHPRIM((*divb),p,i,j,k);}
    else{ SETFDIVBFLUXCTTOTH((*divb),U,i,j,k);}
  }
  else if(FLUXB==FLUXCD){
    if(DOENOFLUX==NOENOFLUX){ SETFDIVBFLUXCDPRIM((*divb),p,i,j,k);}
    else{ SETFDIVBFLUXCD((*divb),U,i,j,k);}
  }
  else if(FLUXB==FLUXCTHLL){
    if(DOENOFLUX==NOENOFLUX){ SETFDIVBFLUXCTTOTHPRIM((*divb),p,i,j,k);} // fake
    else{  SETFDIVBFLUXCTTOTH((*divb),U,i,j,k);} // fake
  }
  else if(FLUXB==FLUXCTSTAG){
    // for STAG, always can use U since always used, but can use pstag too for lower order method
    //if(DOENOFLUX==NOENOFLUX) SETFDIVBFLUXCTSTAG((*divb),pstag,i,j,k);

    if(DOENOFLUX==ENOFLUXRECON && extrazones4emf==0){
      // evolving point field but need Bhat to check divb=0
      // apparently must use fixed stencil so could compute when divb requested instead of always
      // now compute divBhat that should be 0 to machine accuracy
      SETFDIVBFLUXCTSTAG((*divb),Bhat,i,j,k);
    }
    else{
      if(DOENOFLUX==NOENOFLUX){
        // then can use pstag itself that is always bounded properly
        SETFDIVBFLUXCTSTAGPRIM((*divb),pstag,i,j,k);
      }
      else{
        // then must use higher-order version
        // for diagnostics in MPI can bound U before dumping, but should only bound MPI since no treatment of U for real boundaries
        SETFDIVBFLUXCTSTAG((*divb),U,i,j,k);
      }
    }

    //    *divb = 
    //      +(MACP0A1(U,ip1mac(i),j,k,B1)-MACP0A1(U,i,j,k,B1))/dx[1]
    //      +(MACP0A1(U,i,jp1mac(j),k,B2)-MACP0A1(U,i,j,k,B2))/dx[2]
    //      +(MACP0A1(U,i,j,kp1mac(k),B3)-MACP0A1(U,i,j,k,B3))/dx[3]
    //      ;
  }
  else{
    dualfprintf(fail_file,"No such FLUXB==%d in setfdivb()\n",FLUXB);
    myexit(817163);
  }

}







/// Used with FLUXB==FLUXCTHLL
/// actually uses F1/F2/F3 instead of inputted A[]
/// When using FLUXCTHLL, doesn't preserve divb, but gives correct CENT position of field given face vector potential
/// If flux is really vector potential at corners and vector potential is direction-dependent quasi-deaveraged version
int vpot2field_useflux(int *fieldloc,FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int Nvec[NDIM];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];






  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  fluxvec[1] = F1;
  fluxvec[2] = F2;
  fluxvec[3] = F3;



  // Assumes defined vector potential with signature such that
  // B = +curlA
  // \detg B^i = +d_j ([tijk] A_k)
  // and HARM conventions for what F1/F2/F3 mean in terms of EMF-like object
  // Note that dB/dt = -curlE has opposite sign prefactor that is folded into E when doing flux calculation, so flux is really -E
  // So when letting dA/dt = -E = +emf = +flux with consistent cyclic indices



#pragma omp parallel 
  {

    int i,j,k;
    FTYPE igdetgnosing[NDIM];
    int dir;
    int pl,pliter;
    int odir1,odir2;
    int jj;
    struct of_geom geomfdontuse[NDIM];
    struct of_geom *ptrgeomf[NDIM];
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL; // doesn't depend upon dir, but blockijk must be private, so put in parallel loop

    // assign memory
    DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);


    /////////////////
    //
    // A_1
    //
    /////////////////

    dir=1;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ////////      COMPFULLLOOP{ // COMPFULLLOOP allows since A_i exists atCOMPFULLLOOPP1 and so always accessing valid A_i
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // no wait allowed as in fluxct.c
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // ufield doesn't require geometry

        //    myA[3]=-F2[B1]
        //    myA[2]=F3[B1];
        MACP0A1(ufield,i,j,k,B1)  = 0.0;
#if(N2>1)
        MACP0A1(ufield,i,j,k,B1) += -(MACP0A1(F2,i,jp1mac(j),k,B1)-MACP0A1(F2,i,j,k,B1))/(dx[2]);
#endif
#if(N3>1)
        MACP0A1(ufield,i,j,k,B1) += -(MACP0A1(F3,i,j,kp1mac(k),B1)-MACP0A1(F3,i,j,k,B1))/(dx[3]);
#endif

        get_geometry(i, j, k, fieldloc[dir], ptrgeomf[dir]);
        igdetgnosing[dir] = sign(ptrgeomf[dir]->gdet)/(fabs(ptrgeomf[dir]->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];
      }// end 3D LOOP
    }// end if doing A1

    /////////////////
    //
    // A_2
    //
    /////////////////
 
    dir=2;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ///////COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // no wait allowed as in fluxct.c
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        //    myA[1]=-F3[B2];
        //    myA[3]=F1[B2];
      
        MACP0A1(ufield,i,j,k,B2)  = 0.0;
#if(N3>1)
        MACP0A1(ufield,i,j,k,B2) += -(MACP0A1(F3,i,j,kp1mac(k),B2)-MACP0A1(F3,i,j,k,B2))/(dx[3]);
#endif
#if(N1>1)
        MACP0A1(ufield,i,j,k,B2) += -(MACP0A1(F1,ip1mac(i),j,k,B2)-MACP0A1(F1,i,j,k,B2))/(dx[1]);
#endif

        get_geometry(i, j, k, fieldloc[dir], ptrgeomf[dir]);
        igdetgnosing[dir] = sign(ptrgeomf[dir]->gdet)/(fabs(ptrgeomf[dir]->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];
      }
    }

    /////////////////
    //
    // A_3
    //
    /////////////////

    dir=3;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ////////COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // no wait allowed as in fluxct.c
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        //    myA[2]=-F1[B3];
        //    myA[1]=F2[B3];
      
        MACP0A1(ufield,i,j,k,B3)  = 0.0;
#if(N1>1)
        MACP0A1(ufield,i,j,k,B3) += -(MACP0A1(F1,ip1mac(i),j,k,B3)-MACP0A1(F1,i,j,k,B3))/(dx[1]);
#endif
#if(N2>1)
        MACP0A1(ufield,i,j,k,B3) += -(MACP0A1(F2,i,jp1mac(j),k,B3)-MACP0A1(F2,i,j,k,B3))/(dx[2]);
#endif

        get_geometry(i, j, k, fieldloc[dir], ptrgeomf[dir]);
        igdetgnosing[dir] = sign(ptrgeomf[dir]->gdet)/(fabs(ptrgeomf[dir]->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];
      }
    }


  }// end parallel region (and implied barrier)









#if(FLUXDUMP==1)
  {
    int i,j,k,dir,pl,pliter;
    // this accounts for final flux
    FULLLOOP{
      DIMENLOOP(dir){
        if(Nvec[dir]>1){
          PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=MACP1A1(fluxvec,dir,i,j,k,pl);
        }
        else{
          PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=0.0L;
        }
      }
    }
  }
#endif

  

  return(0);
}







////////////////////////////////////
///
/// Evolve A_i (either just updating A_i to be consistent with fluxes/emf, or manipulate A_i/EMFs/Fluxes and then update A_i)
///
/// Like advance() but for A_i with allowed modifications to either EMF or A_i directly
///
/// Can evolve with vpot (i.e. could have modified A_i or may want A_i to be machine accurate instead of drifting from B^i)
/// this used to be in post_advance(), but was excessively doing full evolve_withvpot().  But we already adjusted A_i and can now get corrected Fluxes compared to how computed in fluxcalc_fluxctstag.  The flux hasn't yet been used to compute anything, so this correction is in the right location.
///    if(EVOLVEWITHVPOT){
/// Normally do:
/// a) EMF->A_i
/// b) EMF->F->B in all B's forms (prim,pstag,unew,Bhat).
/// But then if want A_i to be machine error consistent with these B's, need to do A_i->B.  Otherwise, A_i evolution will drift from B's evolution even with same EMFs.
/// In order to avoid repeating A_i->B, can take new and old A_i to create EMF that's used along normal (EMF->F->B) path.  This way one really has (dA_i = A_iold-A_inew)->EMF->F->B.
/// So then don't need to separately adjust EMFs using adjust_fluxctstag_emf(), only have to adjust A_i.  Also solves any inconsistency in how EMF and A_i are modified using adjust_'s.
/// Could decide to only adjust EMFs and not A_i (but then A_i will diverge due to machine error, so haven't solved the original issue -- only solved a redundancy issue).
/// Problem is still there if use dA_i.  B's and A_i will diverge from one another.  Can only have 1 fundamental, or shift fundamental every so often.
///
/// But if want internal BC solution to be smooth extention (not some fake accurate solution), then need to manipulate A_i directly.  Else staggered B^i interpolation will also be jumpy.
/// So let's shift to A_i every step in post_advance() so A_i controls B^i, but let's reset EMF_i from dA_i so don't have to adjust EMF_i directly and separately from A_i.  Otherwise, unless perfectly  match EMF and A_i, may cause inconsistencies -- for example, in B_CENT computed from B_STAG that's first computed from EMFs, while shifted to A_i only after.  So B_CENT would be inconsistent with new B_STAG from A_i
///  }
/// SUPERGODMARK: Note: By relying upon A_i, one eventually cumulates catastrophic cancellation as (e.g.) dA_2+=EMF_2 dt such that after double precision (10^{14} or so) steps the calculation is totally corrupted.  Long before this, machine errors become larger and larger.
int evolve_vpotgeneral(int whichmethod, int stage,
                       int initialstep, int finalstep,
                       FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                       int *Nvec,
                       FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],
                       FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                       FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime,
                       FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]
                       )
{
  int dimen;
  FTYPE (*vpot0)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  FTYPE (*vpotlast)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  FTYPE (*vpotcum)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  vpot0   =vpot+1*NDIM;
  vpotlast=vpot+2*NDIM;
  vpotcum =vpot+3*NDIM;



  if(EVOLVEWITHVPOT==0 && TRACKVPOT==0){
    dualfprintf(fail_file,"Shouldn't be in evolve_vpotgeneral() with EVOLVEWITHVPOT==0 && TRACKVPOT==0\n");
    myexit(45986252);
  }


  //////////////////
  //
  // Can either:
  // 1) user adjust only EMF
  // 2) user adjust both EMF and VPOT.  Reason both: EMF reset on surfaces used to set new VPOT.  VPOT adjusted to deal with interior so Bstag^i and Bcent^i are smooth through surfaces.
  // In case VPOT is adjusted, EMF is fixed-up to agree with VPOT by updating temporary VPOT that is then modified that is then used to back-create the EMF that is used to update to a final VPOT
  // This ensures consistency between user EMF and user VPOT modifications
  //
  //////////////////



  // user adjust of EMFs
  if(MODIFYEMFORVPOT==MODIFYEMF || MODIFYEMFORVPOT==MODIFYVPOT) adjust_emfs(fluxtime, whichmethod, pr, Nvec, fluxvec, emf);




  if(MODIFYEMFORVPOT==MODIFYVPOT || TRACKVPOT>0 || EVOLVEWITHVPOT>0){
    if(initialstep==1){
      // For all current TIMEORDERs, Ui is always from pi so always from last ucum or vpotcum
      for(dimen=0;dimen<NDIM;dimen++){
        // copy over vpotcum -> vpot0
        copy_3dvpot_fullloopp1(vpot[dimen],vpot0[dimen]); // really copying vpotcum -> vpot0 since already copied vpotcum->vpot (below in this function, or at t=0 vpot is correct to use)
        // must first set 0->{vpotlast,vpotcum} if first step (i.e. vpotcum (like ucum) will add vpot0 as required, so starts at 0, not vpot0, here)
        init_3dvpot_fullloopp1(0.0,vpotlast[dimen]);
        init_3dvpot_fullloopp1(0.0,vpotcum[dimen]);
      }
    }
    else{
      // then no change in vpot0 for current TIMEORDER's
    }
  }


  // EVOLVEWITHVPOT==1 means use A_i as ultimate basis for B^i. So can modify A_i and affect resulting B^i while keeping divB=0
  // MODIFYEMFORVPOT==MODIFYVPOT means allow modification of fluxes/emfs using A_i (no point unless A_i is used to compute B^i)


  if(MODIFYEMFORVPOT==MODIFYVPOT && EVOLVEWITHVPOT>0){  // overall, get new emf/fluxes


    // flux->A_i : get new vpot using fluxes
    // like getting Uf from fluxes
    // vpotcum argument is set to NULL since don't want to update ucum-like quantity yet
    update_vpot(whichmethod, stage, pr, fluxvec, emf, CUf, CUnew, fluxdt, vpot, vpot0, vpotlast, NULL);




    //////////////////////////////
    //
    // User "boundary conditions" to modify A_i's before used
    // They will be used in step_ch.c calling evolve_withvpot() that resets B^i using A_i
    //
    /////////////////////////////

    // A_i -> Amod_i (-> A_i)
    adjust_vpot(fluxtime, whichmethod, pr, Nvec, vpot);



    // now obtain EMF_i from A_i and Aold_i
    // this sets EMFs so that one only has to adjust A_i and not also EMF_i
    // A_i,Aold_i -> EMF_i (in F_i)
    // like deriving flux from Uf and Ui (NOTE: not derived from Ucum)
    set_emfflux(whichmethod, stage,pr,fluxvec,emf,CUf, CUnew, fluxdt,vpot,vpot0,vpotlast);


  }


  if(MODIFYEMFORVPOT==MODIFYVPOT || TRACKVPOT>0 || EVOLVEWITHVPOT>0){ // overall, update A_i using emf/fluxes (either original or new from modified A_i)

    ////////////////
    //
    // Update A_i from EMF_i (directly computed, and then either modified directly or through A_i)
    // Before higher-order operations on flux, track vector potential update
    // so updating point value of A_i
    //
    ////////////////

    // EMF_i (in F_i) -> A_i
    // like regetting Uf and Ucum from fluxes
    // only cumulate at this "real" update_vpot() call
    update_vpot(whichmethod, stage, pr, fluxvec, emf, CUf, CUnew, fluxdt, vpot, vpot0, vpotlast, vpotcum);



    // update vpotlast
    for(dimen=0;dimen<NDIM;dimen++) copy_3dvpot_fullloopp1(vpot[dimen],vpotlast[dimen]);


    if(finalstep==1){
      // Doing this for 3 reasons:
      // 1) Copy to vpot so (at top of this function) the (vpotcum->vpot)->vpot0 sets up vpot0
      // 2) In post_advance() will use vpot to recompute all B's, so copy over vpotcum -> vpot so vpot ready for that computation (otherwise could reference vpotcum there)
      // 3) vpotcum will be what's reported in diagnostics with this assignment
      for(dimen=0;dimen<NDIM;dimen++) copy_3dvpot_fullloopp1(vpotcum[dimen],vpot[dimen]);
    }


  }



  return(0);
}






/// wrapper for adjust_fluxctstag_emfs() and adjust_fluxcttoth_emfs()
/// calls user functions
void adjust_emfs(SFTYPE fluxtime, int whichmethod, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{


  if((whichmethod==FLUXCTTOTH)||(whichmethod==ATHENA2)||(whichmethod==ATHENA1)||(whichmethod==FLUXCD)){
    adjust_fluxcttoth_emfs(fluxtime,pr,emf);
  }
  else{
    adjust_fluxctstag_emfs(fluxtime,pr,Nvec,fluxvec);
  }


}




/// wrapper for adjust_fluxctstag_vpot() and adjust_fluxcttoth_vpot()
/// calls user functions
void adjust_vpot(SFTYPE fluxtime, int whichmethod, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
 
  if((whichmethod==FLUXCTTOTH)||(whichmethod==ATHENA2)||(whichmethod==ATHENA1)||(whichmethod==FLUXCD)){
    adjust_fluxcttoth_vpot(fluxtime,pr,Nvec,vpot);
  }
  else{
    adjust_fluxctstag_vpot(fluxtime,pr,Nvec,vpot);
  }


}




/// update A_i
int update_vpot(int whichmethod, int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt,FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpot0)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpotlast)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpotcum)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  extern int bound_flux_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int dir;
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  FTYPE signforfluxlist[NDIM];
  int odir1[NDIM],odir2[NDIM];
  int locvpot[NDIM][NDIM];
  int Nvec[NDIM];
  int numdirs;







#if(0)
  // DEBUG ONLY:  Bound before vpot updated so EMF is bounded so vpot looks normal
  bound_flux_fluxrecon(stage,pr,Nvec,fluxvec);
#endif



  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;
  

  DIMENLOOP(dir){
    // GODMARK: Should synchronize how these functions operate/generate results
    // initialize directions and variables for cyclic access  
    get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir], &signforfluxlist[dir]);
    // get location of vpot and number of locations to set
    set_location_fluxasemforvpot(dir, &numdirs, &odir1[dir], &odir2[dir], locvpot[dir]);
  }





#pragma omp parallel private(dir) //nothing needed to copyin()
  {
    int is,ie,js,je,ks,ke;
    int i,j,k,pl,pliter;
    int fluxdir,pldir,plforflux;
    FTYPE signforflux;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE igdetvpot;
    FTYPE myemf;

    OPENMP3DLOOPVARSDEFINE;



    // GODMARK: time-part not updated and assumed to be 0
    DIMENLOOP(dir){

      //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
      // since looping over edges (emfs) and flux loop different than emf loop, then expand loops so consistent with both fluxes corresponding to that emf
      is=emffluxloop[dir][FIS];
      ie=emffluxloop[dir][FIE];
      js=emffluxloop[dir][FJS];
      je=emffluxloop[dir][FJE];
      ks=emffluxloop[dir][FKS];
      ke=emffluxloop[dir][FKE];

      fluxdir=fluxdirlist[dir];
      pldir=pldirlist[dir];
      plforflux=plforfluxlist[dir];
      signforflux=signforfluxlist[dir];

#if(DEBUGNSBH)
      dualfprintf(fail_file,"nstep=%ld steppart=%d dir=%d fluxdir=%d plforflux=%d locvpot=%d geomeodir=%d pldir=%d signforflux=%21.15g\n",nstep,steppart,dir,fluxdir,plforflux,locvpot[dir][odir2[dir]],B1-1+odir1[dir],pldir,signforflux);
#endif


      if(Nvec[fluxdir]>1){
 
        // \partial_t A_i = -\detg E_i = EMF
        // e.g. A_1 += -\detg E_1
        // GODMARK: vpot exists atCOMPFULLLOOPP1, but fluxvec does not

        //      if(dir==3){
        //      }


        ///#if((SIMULBCCALC==2)&&(TYPE2==1))
        ///      COMPFZLOOP(is,js,ks)
        ///#else
        //#endif
        ///////      COMPZSLOOP( is, ie, js, je, ks, ke ){ // slightly expanded compared to normal flux calculation due to needing emf that is for 2 different fluxes
      
        OPENMP3DLOOPSETUP( is, ie, js, je, ks, ke );
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // Can use "nowait" since each vpot[dir] setting is independent from prior loops
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);


          // GODMARK: Assume doesn't matter what order odir1/odir2 come in here
          // get_geometry(i, j, k, locvpot[dir][odir2[dir]], ptrgeom); // get geometry at CORN[dir] where emf is located
          // which field ptrgeom->e doesn't matter (see fluxctstag.c)
          // igdetvpot=sign(ptrgeom->e[B1-1+odir1[dir]])/(fabs(ptrgeom->e[B1-1+odir1[dir]])+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet

          // \partial_t A_i = EMF = -\detg E_i where E_i=flux/\detg
          // Note that fluxvec as created in flux.c (etc.) has \detg in front of EMF, so remove for A_i as required
          // Note sign depends upon conventions of how associated flux in either direction with single EMF = -\detg E = dA/dt (hence + sign for cyclic indices)
          // MACP1A0(vpot,dir,i,j,k) += fluxdt*MACP1A1(fluxvec,fluxdir,i,j,k,plforflux)*igdetvpot;

          // e.g. from get_fluxpldirs:
          // dA_1 += dt F3[B2] if N3>1 else dA_1 += dt F2[B3] if N2>1 else don't do


          if((whichmethod==FLUXCTTOTH)||(whichmethod==ATHENA2)||(whichmethod==ATHENA1)||(whichmethod==FLUXCD)){
            // then need to use EMF not fluxes since only EMF is filled -- to update A_i
            // loop setup over fluxes so only reach this once (not twice) so only use EMF once as required
            //     MACP1A0(vpot,dir,i,j,k) += fluxdt*MACP1A0(emf,dir,i,j,k);
            myemf=MACP1A0(emf,dir,i,j,k);
          }
          else{
            // then use fluxes directly to update A_i
            // note that only reach this once (not twice) so only use one flux, not its redundant partner
            //     MACP1A0(vpot,dir,i,j,k) += signforflux*fluxdt*MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);
            myemf=signforflux*MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);
          }

          // update like in advance using flux2dUavg() & dUtoU().  But only 1 flux position rather than difference of fluxes.
          FTYPE ftemp[MAXTIMEORDER]={0.0};
          if(vpot!=NULL) MACP1A0(vpot,dir,i,j,k)        = UFSET      (CUf  ,dt,MACP1A0(vpot0,dir,i,j,k),MACP1A0(vpotlast,dir,i,j,k),myemf,0.0,ftemp); // notice vpotlast[] used here not vpot
          if(vpotcum!=NULL) MACP1A0(vpotcum,dir,i,j,k) += UCUMUPDATE(CUnew,dt,MACP1A0(vpot0,dir,i,j,k),MACP1A0(vpot,dir,i,j,k)    ,myemf,0.0,ftemp); // notice vpot[] used here, just like in dUtoU()


#if(DEBUGNSBH)
          if(dir==2 && i==26 && j==40){
            dualfprintf(fail_file,"inside update_vpot: %21.15g (%21.15g %21.15g %21.15g) dt=%21.15g vpot0=%21.15g vpotlast=%21.15g myemf=%21.15g\n",MACP1A0(vpot,dir,i,j,k),CUf[0],CUf[1],CUf[2],dt,MACP1A0(vpot0,dir,i,j,k),MACP1A0(vpotlast,dir,i,j,k),myemf);
          }
#endif


          // if(dir==3){
          //   dualfprintf(fail_file,"i=%d j=%d k=%d fluxdt=%21.15g vpot=%21.15g fluxvec=%21.15g geomvpot=%21.15g\n",i,j,k,fluxdt,MACP1A0(vpot,dir,i,j,k),MACP1A1(fluxvec,fluxdir,i,j,k,plforflux),geomvpot);
          // }

          // after each full step (for example) could compute magnetic field from vpot so no roundoff error progression and not wasteful for higher-order timestepping


        }// end 3D LOOP
      }// end if fluxdir exists
    }// end over dirs
  }// end over parallel region (and implied barrier)
  
  

  return(0);

}




/// fix EMF_i (in Flux space so ready to be used) from A_i and Aold_i
/// pretty similar to update_vpot(), but inverted assignment.  So removed detailed comments/debug from this function.
int set_emfflux(int whichmethod, int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt,FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpot0)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpotlast)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int dir;
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  FTYPE signforfluxlist[NDIM];
  int odir1[NDIM],odir2[NDIM];
  int locvpot[NDIM][NDIM];
  int Nvec[NDIM];
  int numdirs;





  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;
  

  DIMENLOOP(dir){
    // GODMARK: Should synchronize how these functions operate/generate results
    // initialize directions and variables for cyclic access  
    get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir], &signforfluxlist[dir]);
    // get location of vpot and number of locations to set
    set_location_fluxasemforvpot(dir, &numdirs, &odir1[dir], &odir2[dir], locvpot[dir]);
  }





#pragma omp parallel private(dir) //nothing needed to copyin()
  {
    int is,ie,js,je,ks,ke;
    int i,j,k,pl,pliter;
    int fluxdir,pldir,plforflux;
    FTYPE signforflux;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE igdetvpot;
    FTYPE myemf;

    OPENMP3DLOOPVARSDEFINE;



    // GODMARK: time-part not updated and assumed to be 0
    DIMENLOOP(dir){

      //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
      // since looping over edges (emfs) and flux loop different than emf loop, then expand loops so consistent with both fluxes corresponding to that emf
      is=emffluxloop[dir][FIS];
      ie=emffluxloop[dir][FIE];
      js=emffluxloop[dir][FJS];
      je=emffluxloop[dir][FJE];
      ks=emffluxloop[dir][FKS];
      ke=emffluxloop[dir][FKE];

      fluxdir=fluxdirlist[dir];
      pldir=pldirlist[dir];
      plforflux=plforfluxlist[dir];
      signforflux=signforfluxlist[dir];


      if(Nvec[fluxdir]>1){
 
      
        OPENMP3DLOOPSETUP( is, ie, js, je, ks, ke );
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // Can use "nowait" since each vpot[dir] setting is independent from prior loops
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);


          // \partial_t A_i = EMF = -\detg E_i where E_i=flux/\detg
          // e.g. from get_fluxpldirs:
          // dA_1 += dt F3[B2] if N3>1 else dA_1 += dt F2[B3] if N2>1 else don't do


          myemf = dUfromUFSET(CUf,dt,MACP1A0(vpot0,dir,i,j,k),MACP1A0(vpotlast,dir,i,j,k),MACP1A0(vpot,dir,i,j,k));
   

          if((whichmethod==FLUXCTTOTH)||(whichmethod==ATHENA2)||(whichmethod==ATHENA1)||(whichmethod==FLUXCD)){
            // then need to fill EMF since final flux is non-trivial
            MACP1A0(emf,dir,i,j,k)=myemf;
          }
          else{
            // So just invert from update_vpot(): roughly: dA =  Anew - Aold = dt F
            // below memory flux is always there -- how chosen
            MACP1A1(fluxvec,fluxdir,i,j,k,plforflux)=myemf/signforflux;
            // below memory is not always there, so must check
            if(Nvec[plforflux-B1+1]>1){
              MACP1A1(fluxvec,plforflux-B1+1,i,j,k,B1+fluxdir-1)=-MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);
            }
     
          }


        }// end 3D LOOP
      }// end if fluxdir exists
    }// end over dirs
  }// end over parallel region (and implied barrier)
  
  




  return(0);

}



/// renormalize field using user norm
/// GODMARK: in reality should renormalize A_i and then recompute and can renormalize each A_i independently
/// this type of function assumes renormalizing field energy density in lab-frame
int normalize_field_withnorm(FTYPE norm, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel
  {
    int i,j,k,pl,pliter;
    
    OPENMP3DLOOPVARSDEFINE;

    ///////  COMPZLOOP {
    OPENMP3DLOOPSETUPZLOOP;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // next vpot assignment is independent
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

    
      // normalize primitive
      MACP0A1(prim,i,j,k,B1) *= norm;
      MACP0A1(prim,i,j,k,B2) *= norm;
      MACP0A1(prim,i,j,k,B3) *= norm;

      // normalize conserved quantity
      MACP0A1(ucons,i,j,k,B1) *= norm;
      MACP0A1(ucons,i,j,k,B2) *= norm;
      MACP0A1(ucons,i,j,k,B3) *= norm;

      // normalize staggered field primitive
      if(FLUXB==FLUXCTSTAG){
        MACP0A1(pstag,i,j,k,B1) *= norm;
        MACP0A1(pstag,i,j,k,B2) *= norm;
        MACP0A1(pstag,i,j,k,B3) *= norm;
      }

      // normalize higher-order field
      if(HIGHERORDERMEM){
        MACP0A1(Bhat,i,j,k,B1) *= norm;
        MACP0A1(Bhat,i,j,k,B2) *= norm;
        MACP0A1(Bhat,i,j,k,B3) *= norm;      
      }
    }// end 3D LOOP



    // normalize vector potential if tracking
    if(TRACKVPOT){

      /////////      COMPFULLLOOPP1{
      OPENMP3DLOOPSETUPFULLP1;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // next vpot assignment is independent
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);
 
        MACP1A0(vpot,1,i,j,k) *= norm;
        MACP1A0(vpot,2,i,j,k) *= norm;
        MACP1A0(vpot,3,i,j,k) *= norm;
      }// end 3D LOOP
    } 

  }// end parallel region (and implied barrier)

  return(0);

}






/// used when not using vector potential and just assigning conserved quantities as point values from primitives for field
/// uses global p and pstagscratch
int assign_fieldconservatives_pointvalues(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{



  //////  COMPFULLLOOP{
  ////////////  COMPLOOPF{
#pragma omp parallel
  {
    int i,j,k,pl,pliter;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    struct of_geom geomfdontuse;
    struct of_geom *ptrgeomf=&geomfdontuse;
    
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      if(FLUXB==FLUXCTSTAG){
        PLOOPBONLY(pl){
          // get point value
          get_geometry(i, j, k, FACE1+(pl-B1), ptrgeomf);
          MACP0A1(ucons,i,j,k,pl)=MACP0A1(pstag,i,j,k,pl)*(ptrgeomf->gdet);
        }
      }
      else{
        PLOOPBONLY(pl){
          // assume field is centered
          get_geometry(i, j, k, CENT, ptrgeom);
          MACP0A1(ucons,i,j,k,pl)=MACP0A1(prim,i,j,k,pl)*(ptrgeom->gdet);
        }
      }
    }
  }// end parallel region


  return(0);
}





/// this assigns rough pstag value from p
/// in case not using vector potential
int assignrough_primitive_pstag(int i,int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{

  if(FLUXB==FLUXCTSTAG && pstag!=NULL && p!=NULL){

    int myim1,myjm1,mykm1;
    int myip1,myjp1,mykp1;
    int dir,jj,pl;
    struct of_gdetgeom geomfdontuse[NDIM];
    struct of_gdetgeom *ptrgeomf[NDIM];

    // generally ptr's are different inside parallel block
    DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);

    // limit to avoid where have no data
    myim1=MAX(-N1BND,im1mac(i));
    myjm1=MAX(-N2BND,jm1mac(j));
    mykm1=MAX(-N3BND,km1mac(k));

    myip1=MIN(N1+N1BND-1,ip1mac(i));
    myjp1=MIN(N2+N2BND-1,jp1mac(j));
    mykp1=MIN(N3+N3BND-1,kp1mac(k));

   
    // do FACE1 for B1
    dir=1; pl=B1;
    MACP0A1(pstag,i,j,k,pl)=0.5*(MACP0A1(p,i,j,k,pl)+MACP0A1(p,myim1,j,k,pl));
    get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
    MACP0A1(ucons,i,j,k,B1-1+dir)=MACP0A1(pstag,i,j,k,pl)*ptrgeomf[dir]->IEOMFUNCNOSINGMAC(pl);

    // do FACE2 for B2
    dir=2; pl=B2;
    MACP0A1(pstag,i,j,k,pl)=0.5*(MACP0A1(p,i,j,k,pl)+MACP0A1(p,i,myjm1,k,pl));
    get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
    MACP0A1(ucons,i,j,k,B1-1+dir)=MACP0A1(pstag,i,j,k,pl)*ptrgeomf[dir]->IEOMFUNCNOSINGMAC(pl);

    // do FACE3 for B3
    dir=3; pl=B3;
    MACP0A1(pstag,i,j,k,pl)=0.5*(MACP0A1(p,i,j,k,pl)+MACP0A1(p,i,j,mykm1,pl));
    get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
    MACP0A1(ucons,i,j,k,B1-1+dir)=MACP0A1(pstag,i,j,k,pl)*ptrgeomf[dir]->IEOMFUNCNOSINGMAC(pl);

  }


  

  return(0);
}



/// zero-out field for primitives (and pstag too)
int init_zero_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{



  ////////////  COMPLOOPF{
#pragma omp parallel
  {
    int i,j,k;
    int jj;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      MACP0A1(prim,i,j,k,B1)=MACP0A1(prim,i,j,k,B2)=MACP0A1(prim,i,j,k,B3)=0.0;
      MACP0A1(ucons,i,j,k,B1)=MACP0A1(ucons,i,j,k,B2)=MACP0A1(ucons,i,j,k,B3)=0.0;
      if(HIGHERORDERMEM){
        MACP0A1(Bhat,i,j,k,B1)=MACP0A1(Bhat,i,j,k,B2)=MACP0A1(Bhat,i,j,k,B3)=0.0;      
      }
      if(FLUXB==FLUXCTSTAG){
        MACP0A1(pstag,i,j,k,B1)=MACP0A1(pstag,i,j,k,B2)=MACP0A1(pstag,i,j,k,B3)=0.0;
      }
    }
  }// end parallel region


  if(TRACKVPOT){
    /////   COMPFULLLOOPP1{
#pragma omp parallel
    {
      int i,j,k;
      int jj;
      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULLP1;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        DLOOPA(jj) MACP1A0(vpot,jj,i,j,k) = 0.0;
      }
    } 
  }// end parallel region


  
  return(0);
}
