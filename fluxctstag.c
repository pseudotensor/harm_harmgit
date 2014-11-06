#include "decs.h"


/*! \file fluxctstag.c
  \brief all things related to FLUXB==FLUXCTSTAG

  // Except fluxcalc_standard_4fluxctstag() (and its interpolate_prim_cent2face()) that could be generally used to replace standard method
  // and except more general functions with simple FLUXB==FLUXCTSTAG conditions like getplpr()
  // and except functions that call these things, of course

  ///////////////////////////////////////////////////////
  //
  // OUTLINE of entire staggered field procedure:
  //
  //
  // 1) init.c has user set A_i at CORN_i (and user can control whether A_i set or B^i set via fieldfrompotential[])
  // 2) initbase.c computes B^i @ FACE_i and unew=\detg B^i at FACE_i using vpot2field() and B^i and \detg B^i at CENT through interpolation:
  //   a) Creates higher-order A_i or flux using:
  //     i) vectorpot_useflux() when using flux
  //     ii) vectorpot_fluxreconorfvavg() for FV and FLUXRECON methods otherwise
  //   b) Obtains point field and conserved field using:
  //     i)  vpot2field_staggeredfield() if higher-order method involves higher-order field evolution
  //     ii) vpot2field_useflux() if higher-order method involves evolving point fields
  //   c) Obtains de-averaged point conserved staggered conserved field using deaverage_ustag2pstag()
  //   d) Interpolates staggered FACE1,2,3 field to CENT using interpolate_pfield_face2cent()
  // 3) advance.c: Uses Ui=unew, and flux-updates unew and uf, which are used to obtain B^i and \detg B^i at cent:
  //   a) Computes \detg B^i [ui] at FACE1,2,3 for ui at FACE1,2,3 for field in advance_standard()  (Ui at CENT from pi for non-field)
  //   b) Computes pstag at FACE1,2,3 from updated unew/uf at FACE1,2,3 using deaverage_ustag2pstag() in advance_standard()
  //   c) Computes \detg B^i [upoint] at CENT from pstag at FACE1,2,3 using  interpolate_pfield_face2cent() in advance_standard()
  //
  // So updating conserved \detg B^i at FACE1,2,3 and keeping pstag consistent with this using de-averaging function
  // Inversion is peformed on interpolated B^i at CENT obtained from face pstag


  // NOTES:

  // 1) Presently bound pstag because corner regions need EMF and interpolate face to edge so need boundary face values
  //    For example, for EMF3 we need B1 along dir=2 and B2 along dir=1 so can reconstruct field to corner. This requires (a limited) bounding of pstag.  I don't see a way around this without using CENT to get everything and then violating locality of pstag with the EMFs.
  // In non-staggered scheme those face values could come from centered values.  Would be inconsistent to do that for stag method, so need to find another way if want to only bound CENT primitives.
  // 

  // 2) Note that boundary conditions are applied to staggered field in logical way similar to CENT field (with one exception).
  //    For MPI this is correct.  For periodic this is correct.  For reflecting BC this is correct as long as boundary conditions supplied and (for polar axis) \gdetg=0.0 exactly.  For analytical setting of BCs this is correct.  For outflow this is *sufficient* since don't really need to evolve that last cell
  //    This makes code simpler with this assumption so don't have to extra-evolve the last upper face field value
  // For FLUXRECON, this means bound_flux() must set upper flux so staggered field is set by boundary conditions
  // For IF3DSPCTHENMPITRANSFERATPOLE, this is no longer true.  In this case must evolve that outer cell to evolve B2 at outer pole for any-all i,k,myid.

  // 3) Note that for some problems the staggered method reaches nearly 0 field so that the normalized divb is erroneously large due to hitting machine errors relative to the dominate field strength evolved.  Maybe should output unnormalized divb too
  //    TOTH doesn't have this problem because diffusion is high enough to keep field value large.  In 0-field regions TOTH generates checkerboard pattern with much higher field values and so divb appears to stay small.



  // 4) If A_\phi\propto r^{-p} \theta^0 near the pole, then B2\propto 1/\theta near the pole.  But that would imply infinite energy density because B2hat \propto 1/\theta as well.

  */




static void rescale_calc_stagfield_full(int *Nvec, FTYPE (*pstag)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP]);







/// This vpot2field function is for staggered method to evolve quasi-deaveraged fields (i.e. not point fields) when doing higher order
/// compute field at FACE1,2,3 from vector potential A at CORN1,2,3
/// assumes normal field p
/// assume if 1D then along the 1D direction the field doesn't change and input pfield and ufield are already correct and set
int vpot2field_staggeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR])
{
  int Nvec[NDIM];


#if(FIELDSTAGMEM==0)
  dualfprintf(fail_file,"Set FIELDSTAGMEM==1 to do FLUXCTSTAG\n");
  myexit(7158915);
#endif

  /* flux-staggered */

  // note A_i is always at CORN1,2,3 (edges) for any method, so init.c remains the same
  // all comments same as for centered field -- only change is no averaging of A_i to faces from edges (CORN1,2,3)


  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  ////////////////////////
  //
  // now with double-quasi-deaveraged A_i compute B^i
  //
  // use Nvec in case 1-D then if Nvec[odir1]==1 && Nvec[odir2]==1 then don't assign dir field since constant in time and not zero as would be determined here
  // pfield should be from de-(laterally)-averaged ufield at FACE1,2,3
  //
  ////////////////////////

  // since nowait'ed above (because each loop sets ufield,pfield for each dir that don't depend upon eachother), need barrier here at end to stop before continuing (implied in parallel reigon)
#pragma omp parallel 
  {

    int i,j,k;
    int dir;
    FTYPE igdetgnosing[NDIM];
    int odir1,odir2;
    int jj;
    struct of_gdetgeom geomfdontuse[NDIM];
    struct of_gdetgeom *ptrgeomf[NDIM];
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL; // doesn't depend upon dir, but blockijk must be private


    // generally ptr's are different inside parallel block
    DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);



    dir=1;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      //////////      COMPFULLLOOP{ // COMPFULLLOOP allows since A_i exists at COMPFULLLOOPP1 and so always accessing valid A_i
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // ufield doesn't require geometry
        MACP0A1(ufield,i,j,k,B1)  = +(NOAVGCORN_1(A[3],i,jp1mac(j),k)-NOAVGCORN_1(A[3],i,j,k))/(dx[2]);
        MACP0A1(ufield,i,j,k,B1) += -(NOAVGCORN_1(A[2],i,j,kp1mac(k))-NOAVGCORN_1(A[2],i,j,k))/(dx[3]);

        get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
        set_igdetsimple(ptrgeomf[dir]);
        igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];

      }/// end 3D LOOP
    }// end if
 
    dir=2;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ///////COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        MACP0A1(ufield,i,j,k,B2)  = +(NOAVGCORN_2(A[1],i,j,kp1mac(k))-NOAVGCORN_2(A[1],i,j,k))/(dx[3]);
        MACP0A1(ufield,i,j,k,B2) += -(NOAVGCORN_2(A[3],ip1mac(i),j,k)-NOAVGCORN_2(A[3],i,j,k))/(dx[1]);

        get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
        set_igdetsimple(ptrgeomf[dir]);
        igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];

      }
    }

    dir=3;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      //////////      COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        MACP0A1(ufield,i,j,k,B3)  = +(NOAVGCORN_3(A[2],ip1mac(i),j,k)-NOAVGCORN_3(A[2],i,j,k))/(dx[1]);
        MACP0A1(ufield,i,j,k,B3) += -(NOAVGCORN_3(A[1],i,jp1mac(j),k)-NOAVGCORN_3(A[1],i,j,k))/(dx[2]);

        get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
        set_igdetsimple(ptrgeomf[dir]);
        igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing[dir];

      }
    }
  }// end full parallel region (with implied barrier)



#if(0)
  bound_prim(STAGEM1,t,ufield, 0);
  bound_prim(STAGEM1,t,pfield, 0);
#endif

  

  return(0);
}








/// wrapper for:
/// 1) interpolate FACE 2 CORN
/// 2) loop over dimensions setting field flux dimension-by-dimension using multi-D interpolated CORN quantities
/// At present, original flux as emf is computed like normal flux even if overwritten here, and shouldn't be much more expensive doing that there since primary cost is interpolation whose results are required and used here
int fluxcalc_fluxctstag(int stage,
                        int initialstep, int finalstep,
                        FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                        //   FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                        FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                        FTYPE (*wspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3],
                        FTYPE (*prc)[NSTORE2][NSTORE3][NPR2INTERP],
                        FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP],
                        FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP],
                        struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
                        struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
                        FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                        int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],
                        FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                        FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime, struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM]
                        )
{
  int interpolate_prim_face2corn(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*primface_l)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*primface_r)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                                 //     FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                                 FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                                 struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM],
                                 FTYPE (*prc)[NSTORE2][NSTORE3][NPR2INTERP],
                                 FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP],
                                 FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP],
                                 int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP]);
  int dir;
  int idel, jdel, kdel, face;
  int is, ie, js, je, ks, ke;
  int fluxcalc_fluxctstag_emf_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int dir, int odir1, int odir2, int is, int ie, int js, int je, int ks, int ke, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, struct of_loop (*face2cornloop)[NDIM],
                                 //     FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                                 FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                                 FTYPE (*wspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3]);
  int edgedir, odir1,odir2;
  //  static int firsttime=1;
  int i,j,k;
  int enerregion;
  int *localenerpos;



  // forCOMPZSLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];


  //////////////////////////
  //
  // first obtain pbinterp and pvinterp (point CORN1,2,3 quantities) from pl_ct and pr_ct (FACE1,2,3 quantities)
  // This can't be done per dimension since flux needs multi-D quantities for 2D Riemann problem...hence why stored and outside dimenloop
  //
  ////////////////////////////
  interpolate_prim_face2corn(pr, pl_ct, pr_ct, pvbcorn, cent2faceloop, face2cornloop,prc,pleft,pright,Nvec,dqvec);




  ///////////////////////////////////////////////
  //
  // LOOP OVER EDGES (corresponds to EMF,edge directions, NOT face directions)
  //
  ///////////////////////////////////////////////
  DIMENLOOP(dir){

    edgedir=dir;

    ///////////////////////////////////////////////
    //
    // other dimensions
    // will be setting flux associated with d_t(Bodir1) and d_t(Bodir2) and setting flux(Bdir)->0

    get_odirs(dir,&odir1,&odir2);
    // skip to next dir if 1D such that EMF[dir] not needed since always cancels with itself
    // assumes set to 0 and if set to 0 once then always 0 (assume set to 0 in 
    if(Nvec[odir1]==1 && Nvec[odir2]==1){
      //      if(firsttime){
      // then ensure that really 0 and should remain 0 for entire evolution
      // don't need to set this except once assuming no other code sets to non-zero
      // No!  If those directions are 0 then flux isn't defined nor used
      //COMPZSLOOP( -N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND ){
      //   MACP1A1(fluxvec,odir1,i,j,k,B1-1+odir2)=MACP1A1(fluxvec,odir2,i,j,k,B1-1+odir1)=MACP1A1(fluxvec,dir,i,j,k,B1-1+dir)=0.0;
      // }
      //      }
      continue;
    }

    //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
    // since looping over edges (emfs) and flux loop different than emf loop, then expand loops so consistent with both fluxes corresponding to that emf
    is=emffluxloop[dir][FIS];
    ie=emffluxloop[dir][FIE];
    js=emffluxloop[dir][FJS];
    je=emffluxloop[dir][FJE];
    ks=emffluxloop[dir][FKS];
    ke=emffluxloop[dir][FKE];


    // dir corrsponds to *edge,emf* NOT face
    // CUf[2] is for flux updates
    MYFUN(fluxcalc_fluxctstag_emf_1d(stage, pr, dir, odir1, odir2, is, ie, js, je, ks, ke, fluxvec, CUf[2], face2cornloop[edgedir],pvbcorn,wspeed),"flux.c:fluxcalc()", "fluxcalc_fluxctstag_1d", dir);

#if(PRODUCTION==0)
    trifprintf("%d",dir);
#endif
  }// end DIMENLOOP(dir)


#if(0)
  bound_flux(STAGEM1,t,fluxvec[1],fluxvec[2],fluxvec[3], 0);
  if(N1>1) bound_prim(STAGEM1,t,fluxvec[1], 0);
  if(N2>1) bound_prim(STAGEM1,t,fluxvec[2], 0);
  if(N3>1) bound_prim(STAGEM1,t,fluxvec[3], 0);
#endif


#if(0)
  // SUPERDEBUG: shouldn't be needed -- test

  dualfprintf(fail_file,"nstep=%ld steppart=%d\n",nstep,steppart);
  LOOPF_23{ // debug only
    dualfprintf(fail_file,"j=%d emf[0]=%21.15g emf[N1]=%21.15g diff=%21.15g\n",j,MACP1A1(fluxvec,1,0,j,k,B2),MACP1A1(fluxvec,1,N1,j,k,B2),MACP1A1(fluxvec,1,0,j,k,B2)-MACP1A1(fluxvec,1,N1,j,k,B2));
    //MACP1A1(fluxvec,1,N1,j,k,B2)=MACP1A1(fluxvec,1,0,j,k,B2);
    //    MACP1A1(fluxvec,2,N1,j,k,B1)=-MACP1A1(fluxvec,1,0,j,k,B2);

    //    MACP1A1(fluxvec,1,N1,j,k,B2)=-MACP1A1(fluxvec,2,N1,j,k,B1);
    //=-MACP1A1(fluxvec,2,N1,j,k,B1)=
  }
#endif




  if(EVOLVEWITHVPOT>0 ||  TRACKVPOT>0){
    // Evolve A_i
    evolve_vpotgeneral(FLUXB, stage, initialstep, finalstep, pr, Nvec, fluxvec, NULL, CUf, CUnew, fluxdt, fluxtime, vpot);
  }


  int fluxvpot_modifyemfsuser=0;
  fluxvpot_modifyemfsuser=(EVOLVEWITHVPOT>0 ||  TRACKVPOT>0)&&(MODIFYEMFORVPOT==MODIFYEMF || MODIFYEMFORVPOT==MODIFYVPOT);

  if(fluxvpot_modifyemfsuser==0){// if didn't already call adjust_emfs() in fluxvpot above, have to allow user to be able to still modify emfs calling function directly
    // User "boundary conditions" to modify EMFs/FLUXes
    adjust_fluxctstag_emfs(fluxtime,pr,Nvec,fluxvec);
  }




  //  firsttime=0;
  return(0);

}







/// use global pbcorninterp and pvcorninterp at CORN1,2,3 to obtain EMFs at CORN1,2,3
/// NO interpolation of quantities (except kinda wavespeeds) done here
///
/// assumes 2-D Riemann problem
///
/// When time for flux calculation:
///
/// 1) Wavespeed calculation:
///
///    MACP3A0(wspeed,EOMSETMHD,dir,CMIN/CMAX,i,j,k) : Use global wspeed located at FACE (as computed by flux_standard() and put at FACE by global_vchar())
///
///    CMIN,CMAX correspond to left,right going waves at FACEdir
///    In Del Zanna et al. (2003) equation 43-45, we have
///    alpha_x^+ = max(0,wspeed[dir corresponding to x][CMAX][index],wspeed[dir corresponding to x][CMAX][index-1 in y direction])
///    That is, wspeed located at FACE means only need to take max over 2 remaining speeds since already got wspeed by max,min of CENT speeds -- so realy are going over all 4 states for max,min each for each direction
///
/// 2) Flux calculation
///
///    F1[B1]=F2[B2]=F3[B3]=0
///    F2[B3]=-F3[B2] (+-E1)
///    F1[B3]=-F3[B1] (+-E2)
///    F1[B2]=-F2[B1] (+-E3)
///
///    So only really 3 fluxes to be set corresponding to EMFs E1,E2,E3.
///    If 3D, then all fields use E in orthogonal plane (e.g. Bx uses d_2(E3) and d_3(E2), etc.)
///    If 2D in (say) x-y plane, then B1 only evolves by d_2(E3), B2 only evolves by d_1(E3), and B3 evolves by d_1(E2) and d_2(E1)
///    If 1D in (say) x-dir, then d_1(E3) evolves B2 and d_1(E2) evolves B3 and B1 doesn't evolve.  Hence if 1-D in x-dir don't need to compute E1
///    

/// pvcorninterp and pbcorninterp are defined to be used like below initializaiton in set_arrays.c
///   for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=1;pl<=COMPDIM;pl++) for(m=0;m<NUMCS+1;m++) for(l=0;l<NUMCS;l++)  GLOBALMACP1A3(pvbcorninterp,pl2,i,j,k,pl,m,l)=valueinit;
int fluxcalc_fluxctstag_emf_1d(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int dir, int odir1, int odir2, int is, int ie, int js, int je, int ks, int ke, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE CUf, struct of_loop (*face2cornloop)[NDIM],
                               //          FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                               FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                               FTYPE (*wspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3])
{
  int idel1,jdel1,kdel1;
  int idel2,jdel2,kdel2;
  int Nvec[NDIM];



  ///////////////////////////////////////////////
  //
  // get direction offsets for array accesses
  //
  ///////////////////////////////////////////////
  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  idel1 = fluxloop[odir1][FIDEL];
  jdel1 = fluxloop[odir1][FJDEL];
  kdel1 = fluxloop[odir1][FKDEL];
  
  idel2 = fluxloop[odir2][FIDEL];
  jdel2 = fluxloop[odir2][FJDEL];
  kdel2 = fluxloop[odir2][FKDEL];

  

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
  /////   COMPZSLOOP( is, ie, js, je, ks, ke ){ // slightly expanded compared to normal flux calculation due to needing emf that is for 2 different fluxes

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    struct of_gdetgeom gdetgeomdontuse;
    struct of_gdetgeom *ptrgdetgeom=&gdetgeomdontuse;
    int m,l;
    FTYPE emf2d[COMPDIM-1][COMPDIM-1],c2d[NUMCS][COMPDIM-1];
    FTYPE dB[COMPDIM-1],ctop[COMPDIM-1];
    FTYPE emffinal;
    FTYPE geomcornodir1,geomcornodir2;
    FTYPE topwave[COMPDIM-1],bottomwave[COMPDIM-1];
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
    int fluxmethodlocal;
    FTYPE cmaxfactorodir1;
    FTYPE cmaxfactorodir2;
    FTYPE cminfactorodir1;
    FTYPE cminfactorodir2;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      // defaults
      fluxmethodlocal=fluxmethod[B1]; // assume all B's are same.
      cmaxfactorodir1=1.0;
      cmaxfactorodir2=1.0;
      cminfactorodir1=1.0;
      cminfactorodir2=1.0;
#if(OUTERRADIALSUPERFAST)
      // force effective superfast condition on outer radial boundaries
      if(ISBLACKHOLEMCOORD(MCOORD)){
        if(dir==3 && odir1==1 && startpos[1]+i==totalsize[1]){ fluxmethodlocal=HLLFLUX; cminfactorodir1=0.0; }
        if(dir==3 && odir2==1 && startpos[1]+i==totalsize[1]){ fluxmethodlocal=HLLFLUX; cminfactorodir2=0.0; }
        if(dir==3 && odir1==1 && startpos[1]+i==0){ fluxmethodlocal=HLLFLUX; cmaxfactorodir1=0.0; }
        if(dir==3 && odir2==1 && startpos[1]+i==0){ fluxmethodlocal=HLLFLUX; cmaxfactorodir2=0.0; }
      }
#endif

      // pvcorninterp and pbcorninterp are defined to be used like below initializaiton in set_arrays.c
      //    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=1;pl<=COMPDIM;pl++) for(m=0;m<NUMCS+1;m++) for(l=0;l<NUMCS;l++)  GLOBALMACP1A3(pvbcorninterp,pl2,i,j,k,pl,m,l)=valueinit;


      /////////////////////////////////
      //
      // pvcorn must contain velocity in same frame as field primitive in order to avoid interpolating yet another set of quantities (even can't do just 1 extra since 2 directions to interpolate from and no way to keep symmetry of interpolations without doing 2 extras and (e.g.) averaging)
      // assume final EMF being too large for the local field ($B^2-E^2>0$ or $\gamma^2\sim B^2/(B^2-E^2)$) will not be a problem since extra EMF at CORN would only increase field
      // not formally inconsistent since didn't interpolate other velocity or field to same location
      //
      //  for each edge/EMF (dir):
      //     1) There are 2 values for fields in each directions (odir1,odir2)
      //     2) There are 2x2 values for velocities in directions (odir1,odir2)
      //
      //     3) Need to work out where interpolated quantities go
      // 


      // unsure if signature and overall sign will be right -- GODMARK
      // i,j,k should already be taken into account (i.e. no offsets to add)
      // emf in dir direction
      for(m=0;m<NUMCS;m++) for(l=0;l<NUMCS;l++){
          // emf[+- in odir1][+- in odir2]
          // velocity in same positions as emf
          // for example, emf3[+-x][+-y] = By[+-x]*vx[+-x][+-y] - Bx[+-y]*vy[+-x][+-y]
          // below requires velocity to be lab-frame 3-velocity consistent with lab-frame 3-field primitive

          // see fluxct.c for signature definition of EMF and flux

          // m%3+1 gives next 1->2,2->3,3->1
          // 3-(4-m)%3 = (dir+1)%3+1 gives previous 1->3,2->1,3->2
          // so odir1 is forward cyclic
          // so odir2 is backward cyclic
          // notice that the emf in total uses 4 fields and 4*2 velocities for 12 total quantities
          // not all pbcorn[][?] positions are used (have 3 only need 2 per dir)
          // pvcorn[which corner][which component in pl form][+-odir1][+-odir2]
          // pbcorn[which corner][which component in pl form][+-remaining direction that is not corn nor pl-dir]

          emf2d[m][l] = 
            //     + MACP3A0(pbcorn,dir,B1-1+odir2,m,i,j,k)*MACP4A0(pvcorn,dir,U1-1+odir1,m,l,i,j,k) 
            //     - MACP3A0(pbcorn,dir,B1-1+odir1,l,i,j,k)*MACP4A0(pvcorn,dir,U1-1+odir2,m,l,i,j,k);
            + MACP1A3(pvbcorn,dir,i,j,k,odir2,NUMCS,m)*MACP1A3(pvbcorn,dir,i,j,k,odir1,m,l) 
            - MACP1A3(pvbcorn,dir,i,j,k,odir1,NUMCS,l)*MACP1A3(pvbcorn,dir,i,j,k,odir2,m,l);
        }

      
      ///////////////////////////
      //
      // get wave speeds (these are not interpolated yet to CORNER, they start at FACE regardless of STOREWAVESPEEDS==1,2)
      // note wspeed has NO sign information (for any case, including as set by global_vchar())
      // c[CMIN,CMAX][0=odir1,1=odir2]
      // need to determine i,j,k to choose based upon odir value
      // -?del? since going from FACE to CORN
      // del2 for c2d[][0] since wspeed[odir1] is wave going in odir1-direction, whereas other wavespeed to MAX with is in other (odir2) direction
      if(Nvec[odir1]>1){
        c2d[CMIN][0] = fabs(MAX(0.,MAX(+MACP3A0(wspeed,EOMSETMHD,odir1,CMIN,i,j,k),+MACP3A0(wspeed,EOMSETMHD,odir1,CMIN,i-idel2,j-jdel2,k-kdel2))))*cminfactorodir1;
        c2d[CMAX][0] = fabs(MAX(0.,MAX(+MACP3A0(wspeed,EOMSETMHD,odir1,CMAX,i,j,k),+MACP3A0(wspeed,EOMSETMHD,odir1,CMAX,i-idel2,j-jdel2,k-kdel2))))*cmaxfactorodir1;
      }
      else{
        // GODMARK: shoud just set the offending wavespeed (associated with non-existing dimensional direction) to 1.0 so don't have this conditional
        // then speed doesn't matter, not set, so scale-out
        c2d[CMIN][0] = 1.0;
        c2d[CMAX][0] = 1.0;
      }
      
      if(Nvec[odir2]>1){
        c2d[CMIN][1] = fabs(MAX(0.,MAX(+MACP3A0(wspeed,EOMSETMHD,odir2,CMIN,i,j,k),+MACP3A0(wspeed,EOMSETMHD,odir2,CMIN,i-idel1,j-jdel1,k-kdel1))))*cminfactorodir2;
        c2d[CMAX][1] = fabs(MAX(0.,MAX(+MACP3A0(wspeed,EOMSETMHD,odir2,CMAX,i,j,k),+MACP3A0(wspeed,EOMSETMHD,odir2,CMAX,i-idel1,j-jdel1,k-kdel1))))*cmaxfactorodir2;
      }
      else{
        // then speed doesn't matter, not set, so scale-out
        c2d[CMIN][1] = 1.0;
        c2d[CMAX][1] = 1.0;
      }

      ctop[0]    = MAX(c2d[CMIN][0],c2d[CMAX][0]);
      ctop[1]    = MAX(c2d[CMIN][1],c2d[CMAX][1]);


      //////////////////////////////////
      //
      // compute conserved dissipation term 
      //
      /////////////////////////////////


      // "upper" minus "lower" fields
      // dB[?] corresponds to ? meaning field direction, not interpolation or any other direction
      // dB[0] corresponds to dB["odir1"]
      //      dB[0] = MACP3A0(pbcorn,dir,B1-1+odir1,1,i,j,k) - MACP3A0(pbcorn,dir,B1-1+odir1,0,i,j,k) ;
      dB[0] = MACP1A3(pvbcorn,dir,i,j,k,odir1,NUMCS,1) - MACP1A3(pvbcorn,dir,i,j,k,odir1,NUMCS,0) ;
      // dB[1] corresponds to dB["odir2"]
      //      dB[1] = MACP3A0(pbcorn,dir,B1-1+odir2,1,i,j,k) - MACP3A0(pbcorn,dir,B1-1+odir2,0,i,j,k) ;
      dB[1] = MACP1A3(pvbcorn,dir,i,j,k,odir2,NUMCS,1) - MACP1A3(pvbcorn,dir,i,j,k,odir2,NUMCS,0) ;



      //////////////////////////////////
      //
      // compute emf
      //
      /////////////////////////////////

      // Del Zanna et al. (2003) has opposite sign in equations 44,45 since they define the electric field (E_i) whereas we define EMF=-gdet E_i
      // see fluxct.c comments
      // Sign of dissipative term is as for HLL flux, which since EMF and flux have different sign relationships for each direction/field, we have:
      // 
      // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]  : Dissipative terms: -(dB3) and +(dB2)
      // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]  : Dissipative terms: -(dB1) and +(dB3)
      // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]  : Dissipative terms: -(dB2) and +(dB1)
      //
      // So since dissipative term order is EMF[edgedir] += dB[odir1]  + dB[odir2]
      // then we have:
      // edgedir=1 odir1=2 odir2=3 , so need -dB[1] and +dB[0]
      // edgedir=2 odir1=3 odir2=1 , so need -dB[1] and +dB[0]
      // edgedir=3 odir1=1 odir2=2 , so need -dB[1] and +dB[0]
      // so same formula for all EMFs


      // below topwave[] are positive definite but could be 0.0 and then HLLFLUX formula not right
      topwave[0]=c2d[CMIN][0] + c2d[CMAX][0];
      topwave[1]=c2d[CMIN][1] + c2d[CMAX][1];


      if( (fluxmethodlocal==HLLFLUX) && (topwave[0]>SMALL) && (topwave[1]>SMALL) ){

        bottomwave[0]=1.0/topwave[0];
        bottomwave[1]=1.0/topwave[1];

        // HLL
        emffinal = 
          // non-dissipative term
          + (
             +c2d[CMAX][0]*c2d[CMAX][1]*emf2d[0][0]  // emf has -odir1 -odir2, so wavespeed has +odir1 +odir2
             +c2d[CMAX][0]*c2d[CMIN][1]*emf2d[0][1]  // emf has -odir1 +odir2, so wavespeed has +odir1 -odir2
             +c2d[CMIN][0]*c2d[CMAX][1]*emf2d[1][0]  // emf has +odir1 -odir2, so wavespeed has -odir1 +odir2
             +c2d[CMIN][0]*c2d[CMIN][1]*emf2d[1][1]  // emf has +odir1 +odir2, so wavespeed has -odir1 -odir2
             )*bottomwave[0]*bottomwave[1]
          // dissipative terms
          - EMFDISSIPATION*(
             // dB has d(B[odir2]) so wavespeed has +-odir1  (note d(B[odir1]) for +-odir1 is 0 due to divb=0) (i.e. otherwise would be 4 dissipation terms for 2D Riemann problem)
             c2d[CMIN][0]*c2d[CMAX][0]*bottomwave[0]*dB[1]
             )
          + EMFDISSIPATION*(
             // dB has d(B[odir1]) so wavespeed has +-odir2  (note d(B[odir2]) for +-odir2 is 0 due to divb=0) (i.e. otherwise would be 4 dissipation terms for 2D Riemann problem)
             c2d[CMIN][1]*c2d[CMAX][1]*bottomwave[1]*dB[0]
             )
          ;
      }
      else{ // assume if not HLL then LAXF since only 2 methods for this routine
        // LAXF
        // dB and ctop flip odir's due to divb=0(i.e. otherwise would be 4 dissipation terms for 2D Riemann problem)
        emffinal =
          + 0.25*(emf2d[0][0]+emf2d[0][1]+emf2d[1][0]+emf2d[1][1])
          - EMFDISSIPATION*0.50*(ctop[0]*dB[1] - ctop[1]*dB[0]);
      }





      //////////////////////////////////
      //
      // set the fluxes (final flux is expected to have geometry on it)
      //
      /////////////////////////////////

      // B inside pbcorn was setup with geometry
      // assume that if on singularity where gdet=0 then interpolation gave good value on singularity
      // GODMARK: seems like issue with choice of how to treat gdet is related to presence (or not) of line current on singularity
      //          For example, for POLARAXIS E_z will not be 0 if first used gdet and then interpolated, while if used gdet afterwards then will be 0
      //          Same ambiguity exists in fluxct.c
#if(CORNGDETVERSION==1)
      // even though only want geom.e, geom.e calculation could depend upon many things
      get_geometry_geomeonly(i, j, k, CORN1-1+dir, ptrgdetgeom); // get geometry at CORN[dir] where emf is located

      // GODMARK: why use ptrgeom->e here and ptrgeom->g in most other places?
      geomcornodir1=ptrgdetgeom->EOMFUNCMAC(B1-1+odir1); // which field ptrgeom->e doesn't matter as mentioned below
      geomcornodir2=ptrgdetgeom->EOMFUNCMAC(B1-1+odir2);
#else
      geomcornodir1=1.0;
      geomcornodir2=1.0;
#endif



#if(0) // DEBUG:
      if(i==26 && j==40 && k==0){
        dualfprintf(fail_file,"ORIG: emf2d[1][1]=%21.15g emf2d[1][0]=%21.15g emf2d[0][1]=%21.15g emf2d[0][0]=%21.15g ctop[0]=%21.15g  ctop[1]=%21.15g dB[0]=%21.15g  dB[1]=%21.15g emffinal=%21.15g gdetcorn3=%21.15g\n",emf2d[1][1],emf2d[1][0],emf2d[0][1],emf2d[0][0],ctop[0],ctop[1],dB[0],dB[1],emffinal,ptrgdetgeom->gdet);
        dualfprintf(fail_file,"ORIG: c2d[CMIN][0]=%21.15g c2d[CMAX][0]=%21.15g c2d[CMIN][1]=%21.15g c2d[CMAX][1]=%21.15g\n",c2d[CMIN][0],c2d[CMAX][0],c2d[CMIN][1],c2d[CMAX][1]);
      }
#endif

    


      // see fluxct.c for definitions of signature
      // e.g. edgedir=3 gives F[1][B2]=E3 and F[2][B1]=-E3  F[3][B3]=0, which is correct.
      // Checked that correct for all edgedir's
      // notice that geometries for field-type primitive must all be same for field since otherwise emf ptrgeom->e factor has a mixed association
      if(Nvec[odir1]>1) MACP1A1(fluxvec,odir1,i,j,k,B1-1+odir2) = + emffinal*geomcornodir2;
      if(Nvec[odir2]>1) MACP1A1(fluxvec,odir2,i,j,k,B1-1+odir1) = - emffinal*geomcornodir1;
      if(Nvec[dir]>1) MACP1A1(fluxvec,dir,i,j,k,B1-1+dir)     =   0.0;

#if((WHICHEOM==WITHNOGDET)&&(NOGDETB1>0)||(NOGDETB2>0)||(NOGDETB3>0))
      dualfprintf(fail_file,"Makes no sense to have field with different geometry factor since flux,emf mixes those factors\n");
      myexit(196826);
#endif

      // done!

    }// end 3D LOOP
  }// end parallel region





#if(0)
  bound_flux(STAGEM1,t,fluxvec[1],fluxvec[2],fluxvec[3], 0);
  if(N1>1) bound_prim(STAGEM1,t,fluxvec[1], 0);
  if(N2>1) bound_prim(STAGEM1,t,fluxvec[2], 0);
  if(N3>1) bound_prim(STAGEM1,t,fluxvec[3], 0);
#endif




  return (0);
}






/// wrapper for *continuous* interpolation for FACE_to_CENT
void slope_lim_continuous_e2z(int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *face2centloop)
{
  extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
  extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int loc, int continuous, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
  int pl,pliter;




  // TODOMARK: Can take discrete derivative so edge quantity at center, then c2e fully correct, and correctly interpolates as continuous.
  // Like:
  //       0      1      2      3
  //   B1     B1     B1     B1     B1 
  //->   dB1    dB1    dB1    dB1
  //-> interpolate as center to face giving:
  //->  dl dr  dl dr  dl dr  dl dr  dl dr
  //-> Now have left-right face values as well as center value of dB
  //

  
  PINTERPLOOP(pliter,pl){ // probably only pl=B1+dir-1;
    // difference oriented on difference in dir, otherwise full
  }

  // SUPERGODMARK: preal is located at CENT always, but p2interp will be at FACE, so have to average preal inside to get correct location


  // SUPERGODMARK: Use continuous version, overrides use of (e.g.) PARALINE in favor of MC continuous.
  if(OLDNONCONT==1 && LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    get_loop(INTERPLINETYPE, ENOINTERPTYPE, dir, face2centloop);

    // 1 below means primitives instead of conserved quantities (used for loops)
    // GODMARK: ENOMARK: pstag below may only contain field so can't be used for pressure indicator
    // GODMARK: set_interppoint() inside this function sets starting and ending position for loops, and as set c2e always needs more than e2c for obtaining flux, so leaving as for c2e is fine for now
    slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
      
    // ENOMARK: for ENO should really  use special _e2c_cont() function that assumes continuity is required GODMARK
    //    slope_lim_linetype_e2c_cont(realisinterp, ENOPRIMITIVE, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, NULL, p2interp, pcent);
  }
  else{
    // Should really interpolate such that continuous GODMARK
    // GODMARK: set_interppoint() inside this function sets starting and ending position for loops, and as set c2e always needs more than e2c for obtaining flux, so leaving as for c2e is fine for now
    int loc;
    int continuous;
    if(OLDNONCONT){
      loc=CENT;
      continuous=0;
    }
    else{
      loc=FACE1+dir-1;
      continuous=1;
    }
    get_loop(INTERPPOINTTYPE, ENOINTERPTYPE, dir, face2centloop);
    PINTERPLOOP(pliter,pl){
      slope_lim_pointtype(ENOINTERPTYPE, realisinterp, pl, dir, loc, continuous, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright); // GODMARK: overwritting dq from other type of interpolation
    }
  }



  // TODOMARK: Then after interpolation, just sum-up from left-most-edge-value using derivative at edge to get centered quantities
  // Get final values by doing:
  //
  // for linear interpolations, is continuous linear function within cell
  //
  // Bc = Bfl + (3.0/8.0)*dl + (1.0/8.0)*dr
  //
  // for parabolic interpolations, is continous parabolic function within cell
  //
  // Bc = Bfl + (1.0/3.0)*dB + (5.0/24.0)*dl - (1.0/24.0)*dr
  //
  // assumes original interpolation of dB properly reduces in shocks, so in shock one reduces down to possibly DONOR so dl=dr=0








#if(0)
  bound_prim(STAGEM1,t,pleft, 0);
  bound_prim(STAGEM1,t,pright, 0);
#endif

}





/// wrapper for taking staggered conserved field quantity and obtaining centered field quantity
/// upoint enters as stag and leaves as CENT
/// once this function is done we have:
/// 1) pstag will have correct staggered point value
/// 2) upoint (if needed) will be replaced with center conserved value
/// 3) pfield will contain centered field primitive value
/// Note that no averaging or deaveraging occurs in this function -- everything is points
int interpolate_ustag2fieldcent(int stage, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*upoint)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR])
{
  int ustagpoint2pstag(FTYPE (*upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
  int interpolate_pfield_face2cent(FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucent)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR], struct of_loop *face2centloop, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*prc)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], int *Nvec);
  struct of_loop face2cent[NDIM];
  int finalstep;
  FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP];
  int Nvec[NDIM];



#if(0)
  bound_prim(STAGEM1,boundtime,upoint, 0);
#endif



  // setup spatial sizes
  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  // setup dq's if needed
  // called by advance() that doesn't have dq1,dq2,dq3 passed to it (yet), so access global dq1,dq2,dq3
  dqvec[1]=GLOBALPOINT(dq1);
  dqvec[2]=GLOBALPOINT(dq2);
  dqvec[3]=GLOBALPOINT(dq3);



  // Actually updates field using conserved update (and inverts field)
  // next use of pstagscratch by fluxcalc() will be correct (see  setup_rktimestep() commends for RK4)
  ustagpoint2pstag(upoint,pstag);


  if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;




  /////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Bound pstag
  //
  // must bound pstag at FACE1,2,3 so enough boundary zones to interpolate to get field at CENT
  // note that unlike other conserved quantities, since field trivially inverted there is no issue with possible failure modes or fixups, so only need to bound once
  // That is, we always update field as required for divb=0
  // bound_pstag() takes care of which quantities to bound (only bounding B1,B2,B3)

  bound_pstag(stage, finalstep, boundtime, preal, pstag, upoint, USEMPI);


  // note that ustag isn't bounded, but is used for divb calculation, which is thus only valid at active CENT cells -- but that's all that's in normal dumps unless FULLOUTPUT is used


  // pstagescratch should contain result of deaverage_ustag2pstag() -- gets utoinvert ready for inversion using guess pb
  // other quantities in utoinvert are unchanged by this function (and if needed de-averaging this didn't do it!)
  interpolate_pfield_face2cent(preal,pstag,upoint,pcent,face2cent,dqvec,GLOBALPOINT(prc),GLOBALPOINT(pleft),GLOBALPOINT(pright),Nvec);


#if(0)
  bound_prim(STAGEM1,boundtime,pcent);
#endif

  // pstagscratch is at FACE1,2,3 and is primary field variable
  // pfield is at CENT and is dependent field variable
  //
  // This means pstagscratch must be bounded as a flux1,2,3
  // In MPI use pstag with boundflux()
  // In user bound, user must treat field differently knowing where field is located

  // pfield at CENT does NOT need to be bounded itself since pstag is bounded and then interpolated in extended off-dir range in interpolate_pfield_face2cent() that uses slope_lim_continuous_e2z()

  // So only need to bound pstag like a flux -- need to bound before FACE_2_CENT interpolation


  // now have field at CENT wherever needed for both EMF at CORN1,2,3  AND  have enough for normal fluxes at faces
  // Note don't use CENT field to get field in direction of flux of any quantity -- revert to existing staggered value
  // This is why don't need to bound the field at CENT




  return(0);

}







#define USTAG2PSTAGCONSTRAINEDLOOP 1 // should be 1

/// this function called in initbase.c and during advance()
/// here ustag is point value
/// if dimension doesn't exist, then ustag and pstag are really at CENT effectively
int ustagpoint2pstag(FTYPE (*ustag)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR])
{
  int Nvec[NDIM];
  

  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  //////////////////////////
  //
  // get pstag from unew
  // p should be de-averaged field at FACE1,2,3
  //
  ////////////////////////////



#pragma omp parallel 
  {
    void ustag2pstag(int dir, int i, int j, int k, FTYPE (*ustag)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
    extern void get_stag_startendindices(int *loop, int dir, int *is,int *ie,int *js,int *je,int *ks,int *ke);
    int i,j,k,pl,pliter;
    int dir;
    FTYPE igdetgnosing;
    int is,ie,js,je,ks,ke;
    struct of_gdetgeom gdetgeomfdontuse[NDIM];
    struct of_gdetgeom *ptrgdetgeomf[NDIM];
    int jj;

    OPENMP3DLOOPVARSDEFINE;


    // generally ptr's are different inside parallel block
    DLOOPA(jj) ptrgdetgeomf[jj]=&(gdetgeomfdontuse[jj]);


#if(USTAG2PSTAGCONSTRAINEDLOOP==0)


    // Not well-constrained
    DIMENLOOP(dir){
      pl=B1-1+dir;

      OPENMP3DLOOPSETUPFULL;
      /////////    COMPFULLLOOP{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // "nowait" ok because each pstag[B1,B2,B3] set independently
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
        // get geometry for face pre-interpolated values
        get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgdetgeomf[dir]); // FACE1,FACE2,FACE3 each
        set_igdetsimple(ptrgdetgeomf[dir]);
        igdetgnosing = ptrgdetgeomf[dir]->igdetnosing;      
        MACP0A1(pstag,i,j,k,pl)  = MACP0A1(ustag,i,j,k,pl)*igdetgnosing;
      }
    }



#else // else if constrained loops



    // well-constrained to only update pstag where really want update
    // This ensures don't modify pstag outside well-defined box that defines current computational space
    // do pl==B1
    pl=B1;
    dir=1;
    get_stag_startendindices(Uconsevolveloop, dir, &is,&ie,&js,&je,&ks,&ke);
    //  dualfprintf(fail_file,"dir=%d is=%d ie=%d js=%d je=%d ks=%d ke=%d\n",dir,is,ie,js,je,ks,ke);
    //////  COMPZSLOOP(is,ie,js,je,ks,ke){
    OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // "nowait" ok because each pstag[B1,B2,B3] set independently
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
    
      ustag2pstag(dir, i, j, k, ustag,pstag);
    }



    // do pl==B2
    pl=B2;
    dir=2;
    get_stag_startendindices(Uconsevolveloop, dir, &is,&ie,&js,&je,&ks,&ke);
    //  dualfprintf(fail_file,"dir=%d is=%d ie=%d js=%d je=%d ks=%d ke=%d\n",dir,is,ie,js,je,ks,ke);
    //////  COMPZSLOOP(is,ie,js,je,ks,ke){
    OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      ustag2pstag(dir, i, j, k, ustag,pstag);
    }



    // do pl==B3
    pl=B3;
    dir=3;
    get_stag_startendindices(Uconsevolveloop, dir, &is,&ie,&js,&je,&ks,&ke);
    //  dualfprintf(fail_file,"dir=%d is=%d ie=%d js=%d je=%d ks=%d ke=%d\n",dir,is,ie,js,je,ks,ke);
    //////  COMPZSLOOP(is,ie,js,je,ks,ke){
    OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      ustag2pstag(dir, i, j, k, ustag,pstag);
    }


#endif
  }// end parallel region (with implicit barrier)






#if(0)
  bound_prim(STAGEM1,t,pstag);
#endif


  return(0);
}




/// U staggered field to P staggered field
/// OPTMARK: presume this function gets inlined
void ustag2pstag(int dir, int i, int j, int k, FTYPE (*ustag)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR])
{
  FTYPE igdetgnosing;
  int pl,pliter;
  struct of_gdetgeom gdetgeomfdontuse[NDIM];
  struct of_gdetgeom *ptrgdetgeomf[NDIM];
  int jj;

  // assign memory
  DLOOPA(jj) ptrgdetgeomf[jj]=&(gdetgeomfdontuse[jj]);


  pl=B1-1+dir;

  // get geometry for face pre-interpolated values
  get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgdetgeomf[dir]); // FACE1,FACE2,FACE3 each
  set_igdetsimple(ptrgdetgeomf[dir]);
  igdetgnosing = ptrgdetgeomf[dir]->igdetnosing;


  MACP0A1(pstag,i,j,k,pl)  = MACP0A1(ustag,i,j,k,pl)*igdetgnosing;


}






// if not rescaling, then default is to interpolate \detg B^i rather than B^i -- more accurate for field-aligned flows (e.g. monopole)
// 0 : don't use gdet rescale.  Use normal rescale or no rescale.
// 1 : use gdet rescale
// 2 : use gdet rescale dependent on SPC coordinates.  \detg B1 alone 1-dir, B2 along 2-dir, and B3 along 3-dir (but \detg B3 just as fine)
//
// \detg B1 alone 1-dir for typical split-monopolar type field.
// B2 along 2-dir because B2 and B2hat are regular near pole and division by near-zero would be inaccurate at pole. Assumes B2 flips sign across pole in correct sense in boundary cells as viewed by active cells.  Also for EMF_1 or EMF_3 wouldn't make sense unless also interpolated \detg v2 that makes no sense either.
// B3 or \detg B3 along 3-dir is same.
//
// 1 might be fine because this is only used for face2cent, and using "1" instead of "2" allows even jumps in B2 (e.g. due to sign change) across the pole so CENT version of B2 higher order.
// However, experimentally found 1 leads to MANY more hot->entropy reductions.  So should use "2"
//
///////
#define IFNOTRESCALETHENUSEGDET 1

#define IFNOTRESCALETHENUSEGDETswitch(dir) (IFNOTRESCALETHENUSEGDET==1 || (IFNOTRESCALETHENUSEGDET==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (dir==1 || dir==3))))



/// wrapper for rescale() used for staggered field
static void rescale_calc_stagfield_full(int *Nvec, FTYPE (*pstag)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP])
{
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];


  ////////////////////////////////////////////
  //
  // save choice for interpolations so can restore global variables after done
  {
    int pl,pliter;
    nprlocalstart=npr2interpstart;
    nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];
  }


#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERPFULLNPR2INTERP // requires full copyin() changes npr2interp stuff
  {
    int i,j,k;
    int jj;
    int pl,dir;
    struct of_geom geomfdontuse[NDIM];
    struct of_geom *ptrgeomf[NDIM];
    struct of_gdetgeom gdetgeomfdontuse[NDIM];
    struct of_gdetgeom *ptrgdetgeomf[NDIM];
    extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;


    // generally ptr's are different inside parallel block
    // assign memory
    DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);
    DLOOPA(jj) ptrgdetgeomf[jj]=&(gdetgeomfdontuse[jj]);


    // LOOP over faces (dimensions)
    DIMENLOOP(dir){ // dimen loop because rescale() is at different locations
      //    odir1=dir%3+1;
      //    odir2=(dir+1)%3+1;
    
      //    if(Nvec[dir]==1 || (Nvec[dir]!=1 && (Nvec[odir1]==1 && Nvec[odir2]==1) ) ) continue; // later will copy if no such dimension or such dimension but other dimensions don't exist so field must be constant
      if(Nvec[dir]==1) continue; // later will copy if no such dimension


      pl=B1-1+dir;
      npr2interpstart=npr2interpend=0; npr2interplist[0]=pl; // B1,B2,B3 each


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) // NO, rescale() may set arbitrary p2interp's //nowait // p2interp[B1,B2,B3] set independently
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


#if(RESCALEINTERPFLUXCTSTAG && RESCALEINTERP)
        // get geometry for face pre-interpolated values
        get_geometry(i, j, k, FACE1-1+dir, ptrgeomf[dir]); // FACE1,FACE2,FACE3 each
        rescale(DORESCALE,dir,MAC(pstag,i,j,k),ptrgeomf[dir],MAC(p2interp,i,j,k));
#else
        MACP0A1(p2interp,i,j,k,pl) = MACP0A1(pstag,i,j,k,pl);

        if(IFNOTRESCALETHENUSEGDETswitch(dir)){
          // get geometry for face pre-interpolated values
          get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgdetgeomf[dir]); // FACE1,FACE2,FACE3 each
          MACP0A1(p2interp,i,j,k,pl) *= (ptrgdetgeomf[dir]->gdet);
        }   
#endif
 
      }// end COMPFULLLOOP
    }// end over dirs
  }// end parallel region (with implicit barrier)



  ////////////////////////////////////////////
  //
  // restore choice for interpolations from global variables
#pragma omp parallel
  { // must set npr2interp stuff inside parallel region since threadprivate
    int pl;
    npr2interpstart=nprlocalstart;
    npr2interpend=nprlocalend;
    PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
  }


}










/// interpolates field at FACE to field at CENT assuming field along itself is *continuous*
///
/// |=interface
/// i=zone center of ith zone
///
/// |            |         i                 |
/// |         pstag(i)     i                 |
/// |            |         i                 |
/// |         p2interp(i)  i                 |
/// |            |         i                 |
/// |          dq(i)       i                 |
/// |            |         i                 |
/// |pleft(i)    |     pright(i) pleft(i+1)  |    //  notice odd pleft(i) position
/// |            |         i                 |
/// |            |      p_l p_r              |
/// |            |         i                 |
/// |            |        pcent(i)           |
/// |            |         i                 |
///
/// exactly correct for arbitrary order since input is pstag (point quantity) and output is ucent,pcent (point quantities)
/// ustag only used to get non-field ucent at end of function... otherwise use pstag
/// ucent has ONLY field values updated since ucent may already contain correct values for other quantities (means if not, then outside this function must set ucent to something before using full quantity set!)
/// this function called in initbase.c and during advance()
///
/// INPUTS: Nvec, preal, pstag[B1,B2,B3]
/// OUTPUTS: ucent[B1,B2,B3] (if not NULL), pcent[B1,B2,B3], face2centloop
/// TEMPVARS: dqvec(dir,{B1,B2,B3}), prc[B1,B2,B3], pleft[B1,B2,B3], pright[B1,B2,B3]
/// OPENMPOPTMARK: pleft,pright only exist per dir but used for all dirs, so can't use  nowait when involving pleft/pright
/// We don't use wavespeeds here, so don't worry about wspeedtemp being only for one dir
/// 
///
int interpolate_pfield_face2cent(FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucent)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR], struct of_loop *face2centloop, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*prc)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], int *Nvec)
{
  FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP];
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];



  //////////////////////////////////////
  // PROCEDURE:
  // 1) save interp list
  // 2) setup which quantities to interpolate
  // 3) rescale pstag (no need to unrescale since have separate memory for p2interp)
  // 4) interpolate rescaled pstag (p2interp)
  // 5) obtain pcent from dq or pleft/pright for simple interpolation/averaging method
  // 6) copy fields if not interpolated/rescaled
  // 7) restore interp list

  // GODMARK: in 1D don't have to interpolate face2cent for field in dir-direction since has to be constant

  ////////////////////////////////////////////
  //
  // save choice for interpolations so can restore global variables after done
  {
    int pl,pliter;
    nprlocalstart=npr2interpstart;
    nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];
  }


  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
  if(RESCALEINTERPFLUXCTSTAG && RESCALEINTERP || IFNOTRESCALETHENUSEGDET){
    p2interp=prc; // it's different
    // rescale or multiply by \sqrt{-g}
    rescale_calc_stagfield_full(Nvec, pstag,p2interp);
  }
  else{
    p2interp=pstag; // it's itself
  }




  //////////////////////////
  //
  // LOOP OVER FACES for each dimension
  //
  // Peform interpolation on point quantity (even for ENO/FV method a2c_perps would be used before here)
  //
  ////////////////////////////

#pragma omp parallel OPENMPGLOBALPRIVATEFORGEOMNPR2INTERP
  {
    extern int choose_limiter(int dir, int i, int j, int k, int pl);
    void slope_lim_continuous_e2z(int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *face2centloop);
    int pl,pliter;
    int idel,jdel,kdel;
    int i,j,k;
    int dir;
    int locallim;
    FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
    FTYPE *p2interp_l,*p2interp_r;
    FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
    int usedq;
    int realisinterp;
    int odir1,odir2;
    FTYPE igdetgnosing;
    int is,ie,js,je,ks,ke;
    struct of_geom geomcdontuse;
    struct of_geom *ptrgeomc=&geomcdontuse;
    struct of_geom geomfdontuse[NDIM];
    struct of_geom *ptrgeomf[NDIM];
    struct of_gdetgeom gdetgeomcdontuse;
    struct of_gdetgeom *ptrgdetgeomc=&gdetgeomcdontuse;
    struct of_gdetgeom gdetgeomfdontuse[NDIM];
    struct of_gdetgeom *ptrgdetgeomf[NDIM];
    int jj;
    FTYPE ucentgdet;

    OPENMP3DLOOPVARSDEFINE;

    // assign memory
    DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);
    DLOOPA(jj) ptrgdetgeomf[jj]=&(gdetgeomfdontuse[jj]);

    // Setup rescale pointer reference
#if(RESCALEINTERPFLUXCTSTAG && RESCALEINTERP)
    p2interp_l=pstore_l; // then need separate storage
    p2interp_r=pstore_r;
#else
    p2interp_l=p_l; // p2interp_? is final p_?
    p2interp_r=p_r;
#endif




    // GODMARK:  for now use standard interpolation and just force continuity by averaging (probably leads to unstable algorithm)
    DIMENLOOP(dir){

      //   odir1=dir%3+1;
      //   odir2=(dir+1)%3+1;
    
      ////////////
      //
      // for TESTNUMBER 1002, stag, recon, below causes peak of B2/B3 to be smaller -- NO IDEA WHY SUPERGODMARK
      // for same test, stag, fv method has small peak even with Nvec[dir]==1 only
      //
      //if(Nvec[dir]==1 || (Nvec[dir]!=1 && (Nvec[odir1]==1 && Nvec[odir2]==1) ) ) continue; // later will copy if no such dimension or such dimension but other dimensions don't exist so field must be constant
      //
      ////////////


      // setup which quantity to interpolate
      pl = B1-1+dir;
      npr2interpstart=npr2interpend=0; npr2interplist[0]=pl; // B1,B2,B3 each for each dir=1,2,3




      if(Nvec[dir]==1){


        //////////////////////////
        //
        // Copy if no such dimension or such dimension but other dimensions don't exist so field must be constant
        // copy over those things not interpolating (and hence not rescalig either)
        // unlike FACE_to_EDGE, these degenerate cases result in a trivial copy from the fake FACE to CENT
        // what's constant is \detg B^i, not B^i, so need geometry if field along 1D dir with no other dimensions
        // no idel,jdel,kdel for simplicity
        OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // don't wait for each dir since independent and all memory written to (and used in later loops) is independent for each dir.
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);

          if(pcent!=NULL){
            MACP0A1(pcent,i,j,k,pl)=MACP0A1(pstag,i,j,k,pl);
          }

          // get ucent if required
          if(ucent!=NULL){// OPTMARK: Is this expensive?
            // Note: If WHICHEOM==WITHNOGDET and turned off \detg for fields, then staggered method doesn't work, so ok to assume gdet below and assume standard primitive field such that \detg B^i = conserved quantity
            get_geometry_gdetonly(i, j, k, CENT, ptrgdetgeomc); // final quantity is at CENT
            MACP0A1(ucent,i,j,k,pl)  = MACP0A1(pstag,i,j,k,pl)*(ptrgdetgeomc->gdet); // exactly correct (even for ENO/FV)
          }// end if ucent!=NULL

        }// end 3D LOOP




      }
      else{ // Then Nvec[dir]!=1, so should treat field generally along this direction


   
        // get loop details
        idel = fluxloop[dir][FIDEL];
        jdel = fluxloop[dir][FJDEL];
        kdel = fluxloop[dir][FKDEL];

        //////////////////////
        // interpolate    
        //////////////////////
        realisinterp=0; // since only dealing with fields
        slope_lim_continuous_e2z(realisinterp, dir, idel, jdel, kdel, preal, p2interp, dqvec[dir], pleft, pright, &(face2centloop[dir]));

  
        ///////////////////
        // get p_l p_r
        //////////////////////////////////////
        //
        // interpolate primitive using slope (dq) or directly from pleft and pright
        // For FV: p_left, p_right (indexed by grid cell #) -> p2interp_l, p2interp_r (indexed by interface #) -- atch comment
        //
        // always do since p2interp is good and dq/pleft/pright should have stored quantities or just-computed quantities
        /////////////////////////////////////

        // Assume for now that limiter is not per i,j,k but only per dir (unlike normal interpolation)
        // no HORIZONSUPERFAST here
        locallim=choose_limiter(dir, 0,0,0,pl);
        usedq = usedqarray[locallim];
    
        // using COMPZSLOOP since final CENT quantity is used wherever centered primitive is needed for flux (which is sometimes transverse direction)
        // do maximal loop but avoid going out of bounds when accessing dqvec,pleft,pright
        // NOW do constrained loop
        get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);
        ////////      COMPZSLOOP(is,ie,js,je,ks,ke){
        OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) ///nowait // don't wait since each dir is independent (NO: pleft,pright not functions of dir, so each dir not independent)
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);



#if(OLDNONCONT)
          // note that interpolation from FACE_to_CENT is different than from FACE_to_EDGE or CENT_to_FACE
          // if not rescaling or auto-gdet rescaling, then p2interp_? is already really p_?, otherwise p2interp_? is separate memory from p_?
          if(usedq){
            p2interp_l[pl] = MACP0A1(p2interp,i,j,k,pl) + 0.5 * MACP1A1(dqvec,dir,i,j,k,pl);
            p2interp_r[pl] = MACP0A1(p2interp,i + idel,j + jdel,k + kdel,pl) - 0.5 * MACP1A1(dqvec,dir,i + idel,j + jdel,k + kdel,pl);
          }
          else{
            p2interp_l[pl] = MACP0A1(pright,i,j,k,pl);
            p2interp_r[pl] = MACP0A1(pleft,i+idel,j+jdel,k+kdel,pl);
          }
#else
          // note that interpolation from FACE_to_CENT is different than from FACE_to_EDGE or CENT_to_FACE
          // if not rescaling or auto-gdet rescaling, then p2interp_? is already really p_?, otherwise p2interp_? is separate memory from p_?
          if(usedq){
            p2interp_r[pl] = p2interp_l[pl] = MACP0A1(p2interp,i,j,k,pl) + 0.5 * MACP1A1(dqvec,dir,i,j,k,pl);
          }
          else{
            p2interp_r[pl] = p2interp_l[pl] = MACP0A1(pright,i,j,k,pl); // left at same i,j,k just set equal to right, so can ignore left.
          }
#endif

          // TEST ddq=0 equiv or behavior like Athena that just averages the faces so Bi along i is assumed to be piece-wise continuous linear
          //   p2interp_r[pl] = p2interp_l[pl] = 0.5*(MACP0A1(p2interp,i,j,k,pl) + MACP0A1(p2interp,i+idel,j+jdel,k+kdel,pl));

 


          /////////////////////////////////////
          //
          // after interpolation, unrescale from p2interp to normal primitive 
          //
          /////////////////////////////////////
#if(RESCALEINTERPFLUXCTSTAG && RESCALEINTERP)
          get_geometry(i, j, k, CENT, ptrgeomc); // final quantity is at CENT
          rescale(UNRESCALE,dir,p_l,ptrgeomc,p2interp_l);
          rescale(UNRESCALE,dir,p_r,ptrgeomc,p2interp_r);

          if(ucent!=NULL) MACP0A1(ucent,i,j,k,pl)  = 0.5*(p_l[pl]+p_r[pl])*(ptrgeomc->gdet); // exactly correct (even for ENO/FV)

#elif(IFNOTRESCALETHENUSEGDET)
          get_geometry_gdetonly(i, j, k, CENT, ptrgdetgeomc); // final quantity is at CENT

   
          if(IFNOTRESCALETHENUSEGDETswitch(dir)){
            set_igdetsimple(ptrgdetgeomc);
            igdetgnosing=ptrgdetgeomc->igdetnosing;
            ucentgdet=1.0;
          }
          else{
            igdetgnosing=1.0;
            ucentgdet=(ptrgdetgeomc->gdet);
          }

          // Assign \detg B^i (must come before p_{lr} mod since p2interp_{lr} pointer = p_{lr} pointer
          if(ucent!=NULL) MACP0A1(ucent,i,j,k,pl)=0.5*(p2interp_l[pl]+p2interp_r[pl])*ucentgdet; // go ahead and assign ucent if this method

          // now get B^i from \detg B^i
          // remove \detg used during interpolation
          p_l[pl] = p2interp_l[pl]*igdetgnosing;
          p_r[pl] = p2interp_r[pl]*igdetgnosing;


#else
          // otherwise p2interp_{l,r} are really just pointing to p_l and p_r
          if(ucent!=NULL){
            get_geometry_gdetonly(i, j, k, CENT, ptrgdetgeomc); // final quantity is at CENT
            // Note: If WHICHEOM==WITHNOGDET and turned off \detg for fields, then staggered method doesn't work, so ok to assume gdet below and assume standard primitive field such that \detg B^i = conserved quantity
            MACP0A1(ucent,i,j,k,pl)  = 0.5*(p2interp_l[pl]+p2interp_r[pl])*(ptrgdetgeomc->gdet); // exactly correct (even for ENO/FV)
          }
   
#endif


          if(pcent!=NULL){
            // now set pcent -- GODMARK -- should interpolate such that only 1 continuous value
            // Must preserve divb in 1D Riemann problem, so B^{dir} must be continuous
            // Yes, should use larger stencil and interpolate such that constant. Essentially choose
            // large stencil and effectively choosing which points we trust (not necessarily local points)
            MACP0A1(pcent,i,j,k,pl)=0.5*(p_l[pl]+p_r[pl]);
          }

        }// endCOMPZSLOOP

      }// end if N>1

    }// end DIMENLOOP
  }// end parallel region (with implicit barrier)  






  ////////////////////////////////////////////
  //
  // restore choice for interpolations from global variables
#pragma omp parallel
  { // must set npr2interp stuff inside parallel region since threadprivate
    int pl;
    npr2interpstart=nprlocalstart;
    npr2interpend=nprlocalend;
    PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
  }


#if(0)
  if(pcent!=NULL) bound_prim(STAGEM1,t,pcent);
  if(ucent!=NULL) bound_prim(STAGEM1,t,ucent);
#endif

  return(0);

}






/// slope_lim_face2corn() is provided p2interp and returns pleft/pright
/// differs from slope_lim() by range over which quantities required
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
void slope_lim_face2corn(int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *face2cornloop)
{
  extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
  extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int loc, int continuous, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
  int pl,pliter;
  int interporflux;


  if(extrazones4emf){
    interporflux=ENOINTERPTYPE4EMF;
    
  }
  else{
    interporflux=ENOINTERPTYPE;
  }

  // SUPERGODMARK: primreal is located at CENT always, but p2interp can be at FACE in orthogonal direction, so have to average preal inside to get correct location


  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, interporflux, dir, face2cornloop);
    slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, interporflux, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    int loc=FACE1+dir-1; // CENT relative to direction of interpolation, but FACE1+dir-1 relative to primreal
    int continuous=0;
    get_loop(INTERPPOINTTYPE, interporflux, dir, face2cornloop);
    PINTERPLOOP(pliter,pl){
      slope_lim_pointtype(interporflux,realisinterp, pl, dir, loc, continuous, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
    }
  }


#if(0)
  bound_prim(STAGEM1,t,pleft);
  bound_prim(STAGEM1,t,pright);
#endif


}








/// interpolates certain B and v at FACE1,2,3 to CORN1,2,3
///
///
///
/// MACP1A1(primface_l,dir,i,j,k,pl), etc.
/// pbinterp[l,r][4 things]
/// pvinterp[l,r][u,d][4 things]
///
/// Interpolation list:
///
/// FACE1: B1 -> +-C2 and +-C3   l/r V2 -> lr/+-C3   l/r V3 -> lr/+-C2
/// FACE2: B2 -> +-C1 and +-C3   l/r V1 -> lr/+-C3   l/r V3 -> lr/+-C1
/// FACE3: B3 -> +-C1 and +-C2   l/r V1 -> lr/+-C2   l/r V2 -> lr/+-C1
///
/// Hence taking primface_{l or r} for field and primface_{l and r} for velocity and obtaining
/// pfield[l,r][4]
/// pvelocity[l,r][u,d][4]  where l,r,u,d are determine in some cyclic way for a given quantity
///
/// |=interface
/// i=zone center of ith zone
/// |            |                     |
/// |            |                     |
/// |            |                     |
/// |            |                     |
/// |            |                     |
/// |   i-1,j+1  |        i,j+1        |   i+1,j+1
/// |   i-1,j    |        i,j          |   i+1,j
/// |   i-1,j-1  |        i,j-1        |   i+1,j-1
/// |            |                     |
/// |            |                     |  ^ dir=2
///  -> dir=1                              |
///
/// All inputs are points and all outputs are points
///
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// below from fluxcalc type function:
///
///      for(m=0;m<NUMCS;m++) for(l=0;l<NUMCS;l++){
/// emf[+- in odir1][+- in odir2]
/// velocity in same positions as emf
/// for example, emf3[+-x][+-y] = By[+-x]*vx[+-x][+-y] - Bx[+-y]*vy[+-x][+-y]
/// below requires velocity to be lab-frame 3-velocity consistent with lab-frame 3-field primitive

/// see fluxct.c for signature definition of EMF and flux

/// m%3+1 gives next 1->2,2->3,3->1
/// 3-(4-m)%3 = (dir+1)%3+1 gives previous 1->3,2->1,3->2
/// so odir1 is forward cyclic
/// so odir2 is backward cyclic
///   emf2d[m][l] = 
///     + MACP1A3(pvbcorn,dir,i,j,k,odir2,NUMCS,m)*MACP1A3(pvbcorn,dir,i,j,k,odir1,m,l) 
///     - MACP1A3(pvbcorn,dir,i,j,k,odir1,NUMCS,l)*MACP1A3(pvbcorn,dir,i,j,k,odir2,m,l);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// pvcorn and pbcorn are defined to be used like below initializaiton in set_arrays.c
///    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=1;pl<=COMPDIM;pl++) for(m=0;m<NUMCS+1;m++) for(l=0;l<NUMCS;l++)  GLOBALMACP1A3(pvbcorninterp,pl2,i,j,k,pl,m,l)=valueinit;
///

/// JCM: I originally thought the below as 1 would make little sense, because then one would interpolate (e.g.) \detg B1 through the axis when B1 for a non-tilted dipolar field would be constant.  One would be making B1(\theta)\propto \theta near the pole, which requires higher-order interpolation.
/// However, for tilted dipolar fields, B1 and B3 across the pole will necessarily depend upon resolution (see tilteddipole.nb), and B1 can be flat across the pole looking like a bug even if correct.  The interpolation for B1stag -> B1 across the pole will then be poor and reduce near the flat part.
/// But, interpolating \detg B1 avoids this because it relies on A_\phi,\theta that has no glitch.  The \detg smooths it out.
/// This isn't just a numerical issue, it's about how the field is represented.  The flat B1 across the pole is innaccurate to the continuous solution, but one can still have a robust/accurate discrete method.  The flat part can be just fine in general if viewed correctly and treated optimally.
/// 0 : just use direct Bi in transverse interpolation
/// 1 : use \detg Bi in transverse directions
/// 2 : if SPC use B1 in 2-dir and 3-dir .  use \detg B2 in 1-dir and 3-dir  .  use \detg B3 in 1-dir and 2-dir (because B3 generally blows-up at pole while \detg B3 is flat)
#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD 2
/// 1 and 2 don't work for unknown reasons (code eventually crashes due to polar region) for nsdipole problem

///#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitch(facedir,interpdir) (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==1 || (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (facedir==2&&interpdir==1) || (facedir==3&&(interpdir==1||interpdir==2))))))

/// since only store 1 thing to interpolate in any direction, has to be uniform use of \detg per Bi, which works-out to be ok.
///#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuni(facedir) (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==1 || (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (facedir==2||facedir==3))))

/// B1 along 2 and 3-dir.  3-dir isn't issue as \detg B1 same.  But B1 regular at pole and no sign change.  So shouldn't introduce \detg that would make \detg B1\propto \theta and require obtaining B1 at pole by division by small value.  It would also make cancellation for EMF_3: -v1 B2 + v2 B1 terms harder unless also interpolated \detg v1, but that also wouldn't make sense.

/// Interpolate \detg B2 along 1-dir and 3-dir.  No singularity issue, but may be more accurate to include \detg.

/// Interpolate B3 along 1-dir and 2-dir.  Probably some other interpolation better for 1-dir, but for 2-dir, if interpolate \detg B3 across pole, then would have to interpolate \detg v3 for EMF_1: -v2 B3 + v3 B2 to have consistent cancellation on the polar cut-out.  Otherwise B2 at the pole (which must be regular) would result from wrongly cancelling terms in EMF_1.

/// All uses gdet, except for 2-dir only B3 uses gdet.
///#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuni(facedir) (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==1 || (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (facedir==2))))
#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(facedir,interpdir) (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==1 || (INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (facedir==1 && (interpdir==1 || interpdir==3) || facedir==2 && (interpdir==1 || interpdir==3) ||  facedir==3 && (interpdir==1 || interpdir==3)   )))) // GODMARK: Want to allow interpdir==2 with facedir==3, but not working yet.




////////////////////////////////////////////////////////////////
///
/// whether to include gdet in velocity
/// NOTEMARK: Don't have full control since precompute with or without gdet factor for *both* interpdir's.  So can't choose differently for each interpdir.
/// Choosing 0 means prior default of nothing done to velocity when interpolating it.

/// Choose 1 to always use gdet for all velocity components for all interpolation directions.

/// Choosing 2 currently avoids application on v1.  Assumes also interpolate gdet*B2 and gdet*B3 in fluxctstag and in rescale_interp and bounds extrapolations.
///  Note that odir is always interpdir for velocity, so if only want to ever interpolate uu3 in 2-dir, then that'll never be done here.  Just avoid odir=2 to avoid interpolating v2 across 2-dir.  So that's why have odir==1 || odir==2 for this "2" option.
///  2 works fine.  Not contentious with any dir==2 since never want to interpolate v2 in 2-dir.
#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITY 0

/// e.g. dir=facedir=1 means odir12=23 so v2 and v3 would have gdet applied (v3 crucial for non-axisymmetry near poles).  gdet on v2 is crucial if interpolating gdet*B2 in dir=2 in IFNOTRESCALETHENUSEGDET in interpolate_pfield_face2cent() for choices for field-directed interpolation.
/// e.g. dir=facedir=2 means odir12=13 so v1 and v3 would have gdet applied (v3 crucial for non-axisymmetry near poles).  v1 not crucial.
/// e.g. dir=facedir=3 means odir12=12 so v1 and v2 would have gdet applied (application to v1 and v2 not crucial for non-axisymmetry near poles)
#define INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(facedir,odir) (INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITY==1 || (INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITY==2 && (ISSPCMCOORD(MCOORD)==0 || ISSPCMCOORD(MCOORD)==1 && (odir==1 || odir==3) )))





/// INPUTS: Nvec, pr, primface_l[dir], primface_r[dir]
/// OUTPUTS: pbcorn[dir][side], pvcorn[dir][side][side], cent2faceloop, face2cornloop
/// TEMPVARS: prc pleft pright dqvec[dir]
/// OPENMPOPTMARK: prc (used for p2interp in funny way), pleft,pright only exist per dir but used for all dirs, so can't use  nowait when involving pleft/pright/p2interp
/// Note: wavespeeds already stored, so don't worry about wspeedtmp being only 1 dir at a time.  Further, don't use wavespeeds here, only used directly in EMF calculation
///
int interpolate_prim_face2corn(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*primface_l)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*primface_r)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                               //FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                               FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                               struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM],
                               FTYPE (*prc)[NSTORE2][NSTORE3][NPR2INTERP],
                               FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP],
                               FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP],
                               int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP])
{
  FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP];
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];



  // GODMARK: face2corn duplicates cent2face if in 1D
  // This causes FLUXCTSTAG method to be slightly slower than FLUXCTTOTH method even in 1D, and somewhat still in 2D since comparing simple averaging with full interpolation

  ////////////////////////////////////////////
  //
  // save choice for interpolations
  {
    int pl,pliter;
    nprlocalstart=npr2interpstart;
    nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];
  }


  ///////////////////////////////////////
  //
  // Procedure: (follow interpolate_pfield_face2cent() above)

  //
  // 1) Take primface_l and primface_2 and convert velocity term (have all components) to lab-frame 3-velocity (can't modify primface!!) -- so need another space? (vconemf?) FTYPE BASEMACP0A1(vconemf,N?M,...,NDIM-1);
  // 2) rescale() can be setup, but only for field since rescale() assumes input is primitive and after interpolation wouldn't be able to recover lab-frame 3-velocity cmponents if were using primitive velocity (unless WHICHVEL=VEL3, rarely used so no code for that special case when rescale() could be used)
  //
  // 3) For each of 3 edges, interpolate up and down (into pleft,pright) from each 1 FACE to 2 different CORNs
  //
  // 4) take interpolated pleft,pright and form p_l and p_r type quantities inside pbcorn[] and pvcorn[]
  //
  // 5) 
  //
  //
  //
  //
  //
  //
  //

  // Interpolation list:
  //
  // FACE1: B1 -> +-C2 and +-C3   l/r V2 -> lr/+-C3   l/r V3 -> lr/+-C2
  // FACE2: B2 -> +-C1 and +-C3   l/r V1 -> lr/+-C3   l/r V3 -> lr/+-C1
  // FACE3: B3 -> +-C1 and +-C2   l/r V1 -> lr/+-C2   l/r V2 -> lr/+-C1

  // loop over each edge(CORN1,2,3), interpolating in both odir1,odir2 directions and if that odir1,2 direction has N[odir1,2]==1, then just copy instead of interpolate
  // only have 1 dq,pleft,pright for NPR2INTERP quantities

  // sort above list by edge(C1,C2,C3):

  //  +-C1: FACE2 B2 interpolated in dir=3 ...........


  ///////////
  // interpolate dimension-by-dimension to allow future optimizations based upon memory localization
  // also this allows feeding in primface_l and primface_r directly
  ///////////

  // total of 6 sets of interpolations
  // note have chosen to mix velocity,field interpolations with certain symmetry.  In the end this means velocity is interpolated along its own direction only -- maybe not best?
  // note primface_l[i][Bi]=primface_r[i][Bi]
  // note not using many of the primface_l and primface_r (e.g. primface_lr[1][B2] -- cross fields since less local compared to directly evolved staggered field), but others are simply not used
  //
  // +-dir=1 : FACE2_to_CORN3 1 B2 & 2 V1 (primface_l[2][B2] and primface_lr[2][U1]) , FACE3_to_CORN2 1 B3, 2 V1 (primface_l[3][B3] and primface_lr[3][U1])
  // +-dir=2 : FACE1_to_CORN3 1 B1 & 2 V2 (primface_l[1][B1] and primface_lr[1][U2]) , FACE3_to_CORN1 1 B3, 2 V2 (primface_l[3][B3] and primface_lr[3][U2])
  // +-dir=3 : FACE1_to_CORN2 1 B1 & 2 V3 (primface_l[1][B1] and primface_lr[1][U3]) , FACE2_to_CORN1 1 B2, 2 V3 (primface_l[2][B2] and primface_lr[2][U3])


  /////////////////////////////////////////////////
  //
  // before computing EMF, need to convert velocity primitive to something consistent with field primitive frame
  // i.e. B^i (lab-frame 3-field) is normal field primitive, so use v^i (lab-frame 3-velocity) for velocity
  // otherwise number of interpolations increases to 2*4 extra per edge(CORN) since need symmetric value of the other velocity or uu0 if sending ucon
  // so just send lab-frame 3-velocity v^i as consistent with lab-frame field B^i in Faraday,Maxwell tensors

  // instead of using rescale(), just assume always better to inteprolate gdet B^i since the flux is conserved
#define NUMSTAGINTERP 6
#define BFACEINTERP 0
#define BFACEINTERPODIR1 0
#define BLFACEINTERP 0 // same as BFACEINTERP
#define BRFACEINTERP 1 // only used if Nvec[dir]==1 [No, but used to hold non-gdet applied version if required]
#define BFACEINTERPODIR2 1 // used to hold standard B field, not gdet*B, in case required.
#define VLODIR1INTERP 2
#define VRODIR1INTERP 3
#define VLODIR2INTERP 4
#define VRODIR2INTERP 5


#if(NUMSTAGINTERP>NPR2INTERP)
#error "Cannot have NUMSTAGINTERP>NPR2INTERP.  Create new memory space if have to."
#endif


  // holds quantities prepared for interpolation with space used as listed above
  p2interp=prc;
 


  // Loop over faces that contains the original interpolation so only have to compute \detg B^i and v^i once instead of multiple times for each +-dir=1,2,3
  // FACE1: B1 -> +-C2 and +-C3   l/r V2 -> lr/+-C3   l/r V3 -> lr/+-C2
  // FACE2: B2 -> +-C1 and +-C3   l/r V1 -> lr/+-C3   l/r V3 -> lr/+-C1
  // FACE3: B3 -> +-C1 and +-C2   l/r V1 -> lr/+-C2   l/r V2 -> lr/+-C1






#pragma omp parallel OPENMPGLOBALPRIVATEFORUCONANDGEOMNPR2INTERP
  {
    void slope_lim_face2corn(int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], struct of_loop *face2cornloop);
    int pl,pliter;
    int idel,jdel,kdel;
    int idel1,jdel1,kdel1;
    int idel2,jdel2,kdel2;
    int i,j,k;

    struct of_geom geomfdontuse;
    struct of_geom *ptrgeomf=&geomfdontuse;

    struct of_gdetgeom gdetgeomfdontuse;
    struct of_gdetgeom *ptrgdetgeomf=&gdetgeomfdontuse;

    struct of_gdetgeom gdetgeomcdontuse;
    struct of_gdetgeom *ptrgdetgeomc=&gdetgeomcdontuse;

    struct of_gdetgeom gdetgeomcorndontuse;
    struct of_gdetgeom *ptrgdetgeomcorn=&gdetgeomcorndontuse;

    int dir,interpdir,edgedir;
    int locallim;
    FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
    FTYPE *p2interp_l,*p2interp_r;
    FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
    int enerregion;
    int *localenerpos;
    int odir1,odir2;
    int EMFodir1,EMFodir2;
    struct of_state ql,qr;
    struct of_state *ptrql,*ptrqr;
    FTYPE igdetgnosing,igdetgnosingvel;
    int usedq;
    int whichodir;
    int BFACEINTERPCURRENT;
    int Aodir1,Aodir2;
    int Bodir1,Bodir2;
    int Codir1,Codir2;
    int Dodir1,Dodir2;
    extern int choose_limiter(int dir, int i, int j, int k, int pl);
    int realisinterp;
    int is,ie,js,je,ks,ke;
    FTYPE *prface_l,*prface_r;



    OPENMP3DLOOPVARSDEFINE;


    // default state pointer (from now on below code should only use ptrql,ptrqr
    ptrql=&ql;
    ptrqr=&qr;

    p2interp_l=pstore_l;
    p2interp_r=pstore_r;


    ////////////////////////////////////////////
    //
    // Loop over faces
    //
    ///////////////////////////////////////////
    DIMENLOOP(dir){




      ///////////////////////////////////////////////
      //
      // other dimensions
      get_odirs(dir,&odir1,&odir2);

      //  if(Nvec[odir1]==1 && Nvec[odir2]==1) then EMF[dir]==0, but still copy over results in accordance with this routine since results here are for other EMFs
      // for example, if 2D in dirs=1,2 and 1-D in dir=3, then FACE3 still gives EMF1,EMF2 whose differences in dirs=1,2 is used for evolution
      // later explicit interpolations are avoided if interpdir is 1-D

      ///////////////////////////////////////////////
      //
      // get direction offsets for array accesses
      //
      ///////////////////////////////////////////////
      idel1 = fluxloop[odir1][FIDEL];
      jdel1 = fluxloop[odir1][FJDEL];
      kdel1 = fluxloop[odir1][FKDEL];
  
      idel2 = fluxloop[odir2][FIDEL];
      jdel2 = fluxloop[odir2][FJDEL];
      kdel2 = fluxloop[odir2][FKDEL];




      ////////////////////////////////////////////
      //
      // Convert and store primitive into p2interp[] -- these are *all* the quantities needed from this face-dir
      //
      ///////////////////////////////////////////

      // DEBUG:
      //    dualfprintf(fail_file,"%d %d : %d %d : %d %d\n",cent2faceloop[dir].is, cent2faceloop[dir].ie, cent2faceloop[dir].js, cent2faceloop[dir].je, cent2faceloop[dir].ks, cent2faceloop[dir].ke);

      // Was using:
      //     COMPZSLOOP( -N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND ){
      // But, need loops to be controlled to don't access beyond where left-right faces set since ucon can fail with arbitrary input or nan's unlike other types of calculations
      // cent2faceloop is over CENT positions, which now since accessing the faces needs to be transcribed to FACE values interior to those CENT values
      // One shifts up in interpdir (dir) flux direction because interpolated -1 extra CENT cell and N CENT cell to get flux at 0 and N.
      // Note that +SHIFT term just translates pleft/pright type locations to p_l/p_r type locations
      is=cent2faceloop[dir].is + SHIFT1*(dir==1);
      ie=cent2faceloop[dir].ie;
      js=cent2faceloop[dir].js + SHIFT2*(dir==2);
      je=cent2faceloop[dir].je;
      ks=cent2faceloop[dir].ks + SHIFT3*(dir==3);
      ke=cent2faceloop[dir].ke;

      /////////COMPZSLOOP( is, ie, js, je, ks, ke ){
      OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) ///nowait // can't use "nowait" when using p2interp for each dir in next loop
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        // setup which primeface_l and primeface_r to use
        // if Nvec[dir]==1, handled specially later
        prface_l=MACP1A0(primface_l,dir,i,j,k);
        prface_r=MACP1A0(primface_r,dir,i,j,k);

 
        // get geometry for face pre-interpolated values

 
#if(STOREFLUXSTATE==0)
        get_geometry(i, j, k, FACE1-1+dir, ptrgeomf); // at face[dir]

        // VELs: note don't use velocity in "dir" direction
        // LEFTVEL: compute and store v^i at face (input to ucon_calc() is primitive list as correct in primface_l,r)
        MYFUN(ucon_calc(prface_l, ptrgeomf, ptrql->ucon, ptrql->others) ,"flux.c:interpolate_face2corn()", "ucon_calc()", 1);
        // RIGHTVEL: compute and store v^i at face
        MYFUN(ucon_calc(prface_r, ptrgeomf, ptrqr->ucon, ptrqr->others) ,"flux.c:interpolate_face2corn()", "ucon_calc()", 2);

#else

        // really only need i,j,k in geomf for get_stateforfluxcalc(), unless doing INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELD!=0
        ptrgeomf->i=i;
        ptrgeomf->j=j;
        ptrgeomf->k=k;
        ptrgeomf->p=FACE1-1+dir; // "p" not used by get_stteforfluxcalc(), but ISLEFT/ISRIGHT is w.r.t. just offset from face but located at the position of the face

        get_stateforfluxcalc(dir, ISLEFT, prface_l, ptrgeomf, &ptrql);
        get_stateforfluxcalc(dir, ISRIGHT, prface_r, ptrgeomf, &ptrqr);

        // Just always do to avoid complicated conditional
        get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgdetgeomf); // at face[dir]


#endif


        MACP0A1(p2interp,i,j,k,BFACEINTERPODIR2) = MACP0A1(p2interp,i,j,k,BFACEINTERPODIR1) = prface_l[B1-1+dir]; // note that prface_l[dir]=prface_r[dir] at face[dir]
        // p2interp located at face[dir] before interpolation, so gdet should be there
        // BFACE: compute and store \detg B^i and prepare for 2-way interpolation of that single field (notice that primface_l,r same for face field
        // note that since interpolating \detg B^i, don't have to unrescale because can just use this to obtain EMF w/ gdet
        // odir1 and odir2 here refer to possible choices for interpdir, which is always different than [and transverse to] dir.
        if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(dir,odir1)){
          MACP0A1(p2interp,i,j,k,BFACEINTERPODIR1) *= (ptrgdetgeomf->gdet);
        }
        if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(dir,odir2)){
          MACP0A1(p2interp,i,j,k,BFACEINTERPODIR2) *= (ptrgdetgeomf->gdet);
        }


        // VELs: note don't use velocity in "dir" direction
        // GODMARK: Could/Should allow application of rescale on v^i, such as gdet.  This would allow (e.g.) SPC to handle uu3 and B3 the same with gdet, so that interpolation doesn't operate on divergent quantity and operates consistently on these so EMF1 is correctly cancels terms.
        // LEFTVEL: compute and store v^i at face (input to ucon_calc() is primitive list as correct in primface_l,r)
        MACP0A1(p2interp,i,j,k,VLODIR1INTERP) = (ptrql->ucon[odir1])/(ptrql->ucon[TT]);
        MACP0A1(p2interp,i,j,k,VLODIR2INTERP) = (ptrql->ucon[odir2])/(ptrql->ucon[TT]);
        // RIGHTVEL: compute and store v^i at face
        MACP0A1(p2interp,i,j,k,VRODIR1INTERP) = (ptrqr->ucon[odir1])/(ptrqr->ucon[TT]);
        MACP0A1(p2interp,i,j,k,VRODIR2INTERP) = (ptrqr->ucon[odir2])/(ptrqr->ucon[TT]);

        // See comments related to switch prior to this whole function
        if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir1)){
          MACP0A1(p2interp,i,j,k,VLODIR1INTERP) *= (ptrgdetgeomf->gdet);
          MACP0A1(p2interp,i,j,k,VRODIR1INTERP) *= (ptrgdetgeomf->gdet);
        }
        if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir2)){
          MACP0A1(p2interp,i,j,k,VLODIR2INTERP) *= (ptrgdetgeomf->gdet);
          MACP0A1(p2interp,i,j,k,VRODIR2INTERP) *= (ptrgdetgeomf->gdet);
        }
      



      }// end COMPZSLOOP






      ///////////////
      //
      // now send relevant terms in p2interp to interpolator for each odir1,odir2 directions
      //
      ///////////////








      // MACP4A0(pvcorn,corner/emf/edge dir,i,j,k,which velocity,l/r,u/d)
      // note that emf[+- in EMFodir1][+- in EMFodir2] implies pvcorn[l/r][u/d] =  pvcorn[l/r in EMFodir1][u/d in EMFodir2] since have matching [l][m] when accessing emf[] and pvcorn[] where in this comment edgedir=dir and odir1 and odir2 are interpdir and dir (order?)
      //
      // translation:
      //
      // in EMF calculation function:
      //   edgedir=1 odir1=2 odir2=3
      //   edgedir=2 odir1=3 odir2=1
      //   edgedir=3 odir1=1 odir2=2
      //
      // in this function for interpdir=odir1 interpolation:
      //   face dir=1 odir1=interpdir=2 odir2=edgedir=3  EMFodir1=1  EMFodir2=2
      //   face dir=2 odir1=interpdir=3 odir2=edgedir=1  EMFodir1=2  EMFodir2=3
      //   face dir=3 odir1=interpdir=1 odir2=edgedir=2  EMFodir1=3  EMFodir2=1
      //
      //
      // Hence, EMFodir1 is facedir and EMFodir2 is interpdir
      //
      // in this function for interpdir=odir2 interpolation:
      //   face dir=1 odir2=interpdir=3 odir1=edgedir=2  EMFodir1=3  EMFodir2=1
      //   face dir=2 odir2=interpdir=1 odir1=edgedir=3  EMFodir1=1  EMFodir2=2
      //   face dir=3 odir2=interpdir=2 odir1=edgedir=1  EMFodir1=2  EMFodir2=3
      //
      //
      // Hence, EMFodir1 is interpdir and EMFodir2 is facedir
      //
      // So need to define EMFodir1 and EMFodir2 from edgedir



      //////////////////////////////////////////
      //
      // Loop over other directions not in face-dir
      //
      //////////////////////////////////////////

      for(whichodir=0;whichodir<=1;whichodir++){

        if(whichodir==0){ // whichodir==0 arbitrarily corresponds to interpdir=odir1
          ///////////////////////////
          // interpolate in odir1 direction (places quantities at edge/emf/corner odir2)

          npr2interpstart=0;
          npr2interpend=2; // 3 things
          npr2interplist[0]=BFACEINTERPODIR1; // always interpdir field
          BFACEINTERPCURRENT=npr2interplist[0];
          // ODIR? should be same as interpdir
          npr2interplist[1]=VLODIR1INTERP; // hence [1] is for previous p_l
          npr2interplist[2]=VRODIR1INTERP; // hence [2] is for previous p_r

          // face is dir, and interpolation direction and edge direction are orthogonal to that
          interpdir=odir1;
          edgedir=odir2;

          // del's associated with interpolation direction
          idel=idel1;
          jdel=jdel1;
          kdel=kdel1;
        }
        else{

          ///////////////////////////
          // interpolate in odir2 direction (places quantities at edge/emf/corner odir1)

          npr2interpstart=0;
          npr2interpend=2; // 3 things
          npr2interplist[0]=BFACEINTERPODIR2; // always interpdir field
          BFACEINTERPCURRENT=npr2interplist[0];
          // ODIR? should be same as interpdir
          npr2interplist[1]=VLODIR2INTERP; // hence [1] is for previous p_l
          npr2interplist[2]=VRODIR2INTERP; // hence [2] is for previous p_r

          interpdir=odir2;
          edgedir=odir1;

          // del's associated with interpolation direction
          idel=idel2;
          jdel=jdel2;
          kdel=kdel2;
        }



        //////////////////////
        //
        // set EMFodir's
        //
        //////////////////////
        get_odirs(edgedir,&EMFodir1,&EMFodir2);



        //////////////////////
        //
        //    Setup access to emf-like quantities that have 3-dimensions and going to 2-dimensions
        //
        //////////////////////

        if(EMFodir1==interpdir){
          // then first entry should contain 0/1 for *current* interpolation
          // then EMFodir1 is interpdir corresponding to direction for current interpolation
          // then EMFodir2 is face-dir corresponding to direction for previous interpolation
          Aodir1=0; Aodir2=0;
          Bodir1=1; Bodir2=0;
          Codir1=0; Codir2=1;
          Dodir1=1; Dodir2=1;
        }
        else{
          // then first entry should contain 0/1 for *previous* interpolation
          // then EMFodir1 is face-dir corresponding to direction for previous interpolation
          // then EMFodir2 is interpdir corresponding to direction for current interpolation
          Aodir1=0; Aodir2=0;
          Bodir1=0; Bodir2=1;
          Codir1=1; Codir2=0;
          Dodir1=1; Dodir2=1;
        }




        //////////////////////
        // interpolate    
        //
        // GODMARK: Note put primface_l[dir] (yes, face-dir) into slope_lim as a "good primitive" to be used for shock or other indicators -- should really use average of primface_l and primface_r for symmetry considerations -- otherwise not used -- ENOMARK
        //
        // only really do interpolation if dimension exists ...
        //////////////////////
        //  (Nvec[interpdir]>1 && (! (Nvec[edgedir]==1 && Nvec[dir]==1) ))      // ... or even if dimension exists but orthogonal dimensions do not then just copy over result -- need to modify more things to do this
        //      if(Nvec[interpdir]>1){

        // if Nvec[dir]==1, then that face quantity is really a centered quantity in some CORN plane.  See below for how assignment is done.

        if(!(Nvec[interpdir]==1 || Nvec[dir]==1&&Nvec[interpdir]!=1   )){
          realisinterp=0; // since only ever limited set of quantities
          slope_lim_face2corn(realisinterp, interpdir,idel,jdel,kdel,pr,p2interp,dqvec[interpdir],pleft,pright, &(face2cornloop[edgedir][EMFodir1][EMFodir2]));
        }



  

        ///////////////////
        // get p_l p_r
        //////////////////////////////////////
        //
        // interpolate primitive using slope (dq) or directly from pleft and pright
        //
        /////////////////////////////////////
 
        // Assume for now that limiter is not per i,j,k but only per dir (unlike normal interpolation) (also no pl dependence)
        // no HORIZONSUPERFAST here
        locallim=choose_limiter(interpdir, 0,0,0,B1);
        usedq = usedqarray[locallim];


        // face2corn is at effective-CENT relative to edgedir
        // loop below will be at effective-FACEs, so extended in interpdir direction
        // One shifts up in interpdir direction because interpolated -1 extra "CENT" cell and N "CENT" cell to get corner at 0 and N.
        if(!(Nvec[interpdir]==1|| Nvec[dir]==1&&Nvec[interpdir]!=1  )){
          is=face2cornloop[edgedir][EMFodir1][EMFodir2].is + SHIFT1*(interpdir==1);
          ie=face2cornloop[edgedir][EMFodir1][EMFodir2].ie;
          js=face2cornloop[edgedir][EMFodir1][EMFodir2].js + SHIFT2*(interpdir==2);
          je=face2cornloop[edgedir][EMFodir1][EMFodir2].je;
          ks=face2cornloop[edgedir][EMFodir1][EMFodir2].ks + SHIFT3*(interpdir==3);
          ke=face2cornloop[edgedir][EMFodir1][EMFodir2].ke;
        }
        else{
          // since just copying, revert to locations where interpolated CENT -> FACE
          // Note that +SHIFT term just translates pleft/pright type locations to p_l/p_r type locations
          is=cent2faceloop[interpdir].is + SHIFT1*(interpdir==1);
          ie=cent2faceloop[interpdir].ie;
          js=cent2faceloop[interpdir].js + SHIFT2*(interpdir==2);
          je=cent2faceloop[interpdir].je;
          ks=cent2faceloop[interpdir].ks + SHIFT3*(interpdir==3);
          ke=cent2faceloop[interpdir].ke;
          //   is=-N1BND+idel;
          //   ie=N1-1+N1BND;
          //   js=-N2BND+jdel;
          //   je=N2-1+N2BND;
          //   ks=-N3BND+kdel;
          //   ke=N3-1+N3BND;
        }

        // dualfprintf(fail_file,"edgedir=%d EMFodir1=%d EMFodir2=%d :: is=%d ie=%d js=%d je=%d ks=%d ke=%d\n",edgedir,EMFodir1,EMFodir2,is,ie,js,je,ks,ke);
 
        /////// COMPZSLOOP( -N1BND+idel, N1-1+N1BND, -N2BND+jdel, N2-1+N2BND, -N3BND+kdel, N3-1+N3BND ){
        // OPENMP3DLOOPSETUP( -N1BND+idel, N1-1+N1BND, -N2BND+jdel, N2-1+N2BND, -N3BND+kdel, N3-1+N3BND );
        OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) //nowait // Don't wait since each direction is independent // NO: use of pleft, pright, and some use of p2interp (that depend on each dir) from loop above makes not possible to use "nowait"
        OPENMP3DLOOPBLOCK{
          OPENMP3DLOOPBLOCK2IJK(i,j,k);

   
          // if(Nvec[interpdir]>1 && (! (Nvec[edgedir]==1 && Nvec[dir]==1) ))
          // if(Nvec[interpdir]>1){
          if(!(Nvec[interpdir]==1|| Nvec[dir]==1&&Nvec[interpdir]!=1  )){

            if(usedq){
              PINTERPLOOP(pliter,pl){
                // FACE_to_CORN interpolation is same as if doing CENT_to_EDGE from point of view of indicies to use and pleft,pright assignments
                p2interp_l[pl] = MACP0A1(p2interp,i - idel,j - jdel,k - kdel,pl) + 0.5 * MACP1A1(dqvec,interpdir,i - idel,j - jdel,k - kdel,pl);
                p2interp_r[pl] = MACP0A1(p2interp,i,j,k,pl) - 0.5 * MACP1A1(dqvec,interpdir,i,j,k,pl);
              }
            }
            else{
              PINTERPLOOP(pliter,pl){
                p2interp_l[pl] = MACP0A1(pright,i-idel,j-jdel,k-kdel,pl);
                p2interp_r[pl] = MACP0A1(pleft,i,j,k,pl);
              }
            }

          }
          else if(Nvec[interpdir]==1){

            // if no interpolation, just copy result from pre-interpolated p2interp[] avoiding,bypassing dq,pleft,pright
            PINTERPLOOP(pliter,pl){
              p2interp_r[pl] = p2interp_l[pl] = MACP0A1(p2interp,i,j,k,pl);
            }



          }
          else if(Nvec[dir]==1&&Nvec[interpdir]!=1){

            // if interpolation already done, just copy over left-right values from cent2face operation
            // E.g. If N1==1 and dir==1, then B1,v2 are already interpolated with interpdir==2 to CORN3=FACE2 *and* B1,v3 are already interpolated with interpdir==3 to CORN2=FACE3.
            // E.g. If N2==1 and dir==2, then B2,v1 are already interpolated with interpdir==1 to CORN3=FACE1 *and* B2,v3 are already interpolated with interpdir==3 to CORN1=FACE3.
            // E.g. If N3==1 and dir==3, then B3,v1 are already interpolated with interpdir==1 to CORN2=FACE1 *and* B3,v2 are already interpolated with interpdir==2 to CORN1=FACE2.
            // 
            // Since face is the corner we want, so same ptrgeomf is good.  That is, FACE1-1+dir = desired CORN.
            // But needed to change l/r dir from "dir" to "interpdir"

            // Have to compute v^i=u^i/u^t here since didn't know interpdir in earlier part of this whole function
            prface_l=MACP1A0(primface_l,interpdir,i,j,k);
            prface_r=MACP1A0(primface_r,interpdir,i,j,k);


            // get geometry for face pre-interpolated values
#if(STOREFLUXSTATE==0)
            get_geometry(i, j, k, FACE1-1+dir, ptrgeomf); // at face[dir]=CORN(perp to dir and interpdir)
            MYFUN(ucon_calc(prface_l, ptrgeomf, ptrql->ucon, ptrql->others) ,"flux.c:interpolate_face2corn()", "ucon_calc()", 1);
            MYFUN(ucon_calc(prface_r, ptrgeomf, ptrqr->ucon, ptrqr->others) ,"flux.c:interpolate_face2corn()", "ucon_calc()", 2);
#else

            ptrgeomf->i=i;
            ptrgeomf->j=j;
            ptrgeomf->k=k;
            ptrgeomf->p=FACE1-1+dir; // "p" not used by get_stateforfluxcalc()

            get_stateforfluxcalc(interpdir, ISLEFT, prface_l, ptrgeomf, &ptrql);
            get_stateforfluxcalc(interpdir, ISRIGHT, prface_r, ptrgeomf, &ptrqr);

            // Just always do to avoid complicated conditional
            // as when storing, still at original face locations already interpolated/computed during normal flux calculation
            get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgdetgeomf); // at face[dir]


#endif

            // now copy over values
            p2interp_l[BFACEINTERPCURRENT] = prface_l[B1-1+dir];
            p2interp_r[BFACEINTERPCURRENT] = prface_r[B1-1+dir];
            // applying gdet if just copying result so same application of factors as when setting up p2interp[]
            if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(dir,interpdir)){
              p2interp_l[BFACEINTERPCURRENT] *= (ptrgdetgeomf->gdet);
              p2interp_r[BFACEINTERPCURRENT] *= (ptrgdetgeomf->gdet);
            }

            // [1,2] are from previous (i.e. face) interpolation.  But if Nvec[dir=facedir]=1, then those are same values
            // as above setup has, always dealing with v^{interpdir}
            // here, the l,r indicate across interpdir (not dir as otherwise when requiring interpolation)
            p2interp_l[npr2interplist[1]] = p2interp_l[npr2interplist[2]] = (ptrql->ucon[interpdir])/(ptrql->ucon[TT]);
            p2interp_r[npr2interplist[1]] = p2interp_r[npr2interplist[2]] = (ptrqr->ucon[interpdir])/(ptrqr->ucon[TT]);
            // applying gdet if just copying result so same application of factors as when setting up p2interp[]
            if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir1) && npr2interplist[1]==VLODIR1INTERP || INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir2) && npr2interplist[1]==VLODIR2INTERP){
              p2interp_l[npr2interplist[1]] *= (ptrgdetgeomf->gdet);
              p2interp_l[npr2interplist[2]] *= (ptrgdetgeomf->gdet);
              p2interp_r[npr2interplist[1]] *= (ptrgdetgeomf->gdet);
              p2interp_r[npr2interplist[2]] *= (ptrgdetgeomf->gdet);
            }
          }
          else{
            dualfprintf(fail_file,"Shouldn't reach here in fluxctstag.c\n");
            myexit(837434873);
          }



          //////////////////////////////////////
          // DONE with interpolations or copies
          //////////////////////////////////////



          ////////////////////
          //
          // Set lab-frame 3-magnetic field at CORN
          //
          // MACP1A3(pbcorn,corner/emf/edge dir,i,j,k,which field,+-present interpdir)

          int didsetigdet;
          didsetigdet=0;
          if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(dir,interpdir)==1 && CORNGDETVERSION==1){     
            get_geometry_gdetonly(i, j, k, CORN1-1+edgedir, ptrgdetgeomcorn); // at CORN[dir]
            // then unrescale field since will multiply geometry once have final EMF (avoids line currents)
            set_igdetsimple(ptrgdetgeomcorn);
            igdetgnosing = ptrgdetgeomcorn->igdetnosing;
            didsetigdet=1;
          }
          else if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFFIELDswitchuniboth(dir,interpdir)==0 && CORNGDETVERSION==0){
            get_geometry_gdetonly(i, j, k, CORN1-1+edgedir, ptrgdetgeomcorn); // at CORN[dir]
            // then add gdet now since will not multiply geometry once have final EMF
            igdetgnosing = ptrgdetgeomcorn->gdet; // here igdet really is just geom
          }
          else{
            // nothing to do
            igdetgnosing = 1.0;
          }


          //   MACP3A0(pbcorn,edgedir,B1-1+dir,0,i,j,k) = p2interp_l[npr2interplist[0]]*igdetgnosing;
          //   MACP3A0(pbcorn,edgedir,B1-1+dir,1,i,j,k) = p2interp_r[npr2interplist[0]]*igdetgnosing;
          MACP1A3(pvbcorn,edgedir,i,j,k,dir,NUMCS,0) = p2interp_l[npr2interplist[0]]*igdetgnosing;
          MACP1A3(pvbcorn,edgedir,i,j,k,dir,NUMCS,1) = p2interp_r[npr2interplist[0]]*igdetgnosing;




          ////////////////////
          //
          // Set lab-frame 3-velocity at CORN
          //
          //      pvcorn[EMFdir][i][j][k][whichvel][l/r in EMFodir1][u/d in EMFodir2]
          //
          // for example:
          // edgedir=1: EMFodir1=2 EMFodir2=3 then emf[0][0] is emf[left for EMFodir1][left for EMFodir2]
          // so if interpolated in 2-dir and EMFodir1==2, then emf[0/1 filled with current p_l and p_r][0/1 filled with previous p_l p_r]


          // de-applying any gdet factor
          if(INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir1) && npr2interplist[1]==VLODIR1INTERP || INCLUDEGDETINTRANSVERSEINTERPLATIONOFVELOCITYswitchuniboth(dir,odir2) && npr2interplist[1]==VLODIR2INTERP){
            if(didsetigdet){
              igdetgnosingvel=igdetgnosing;
            }
            else{
              get_geometry_gdetonly(i, j, k, CORN1-1+edgedir, ptrgdetgeomcorn); // at CORN[dir]
              set_igdetsimple(ptrgdetgeomcorn);
              igdetgnosingvel = ptrgdetgeomcorn->igdetnosing;
            }
          }
          else igdetgnosingvel=1.0;


          // npr2interplist[1,2] constains [u,d for velocity in interpdir direction] 
          MACP1A3(pvbcorn,edgedir,i,j,k,interpdir,Aodir1,Aodir2) = p2interp_l[npr2interplist[1]] *= igdetgnosingvel; // current p_l for previous p_l 
          MACP1A3(pvbcorn,edgedir,i,j,k,interpdir,Bodir1,Bodir2) = p2interp_r[npr2interplist[1]] *= igdetgnosingvel; // current p_r for previous p_l
          MACP1A3(pvbcorn,edgedir,i,j,k,interpdir,Codir1,Codir2) = p2interp_l[npr2interplist[2]] *= igdetgnosingvel; // current p_l for previous p_r
          MACP1A3(pvbcorn,edgedir,i,j,k,interpdir,Dodir1,Dodir2) = p2interp_r[npr2interplist[2]] *= igdetgnosingvel; // current p_r for previous p_r



        }// endCOMPZSLOOP
  
      }// end loop over (whichodir) other directions // at end of loop, have pbcorn,pvcorn for this 1 face interpolated to 2 corners
      


    }// end DIMENLOOP // at end of loop, have pbcorn,pvcorn for 3 edges


  }// end over parallel region (and implied barrier)



  ////////////////////////////////////////////
  //
  // restore choice for interpolations
  //
  ///////////////////////////////////////////
#pragma omp parallel
  { // must set npr2interp stuff inside parallel region since threadprivate
    int pl;
    npr2interpstart=nprlocalstart;
    npr2interpend=nprlocalend;
    PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
  }






  return(0);

}











/// pure DONOR version:
/// upoint is input staggered updated field that needs to be converted to pstag and then "interpolated" to pcent and upoint that is upon output at CENT
/// To use this, just rename this as without "_donor" and rename normal code to have (e.g.) "_normal" on end of function name.
int interpolate_ustag2fieldcent_donor
//int interpolate_ustag2fieldcent
(int stage, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*upoint)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl;
  int ii,jj,kk;
  struct of_geom geomcdontuse;
  struct of_geom *ptrgeomc=&geomcdontuse;
  struct of_geom geomfdontuse;
  struct of_geom *ptrgeomf=&geomfdontuse;
  struct of_geom geomfudontuse;
  struct of_geom *ptrgeomfu=&geomfudontuse;
  int finalstep;

  //  COMPFULLLOOP{

  // must carefully replace pstag only on specific locations
#define MYOCOMPLOOPF3 for(k=SHIFTX3DN;k<=N3-1+SHIFT3+SHIFTX3UP;k++)
#define MYOCOMPLOOPF2 for(j=SHIFTX2DN;j<=N2-1+SHIFT2+SHIFTX2UP;j++)
#define MYOCOMPLOOPF1 for(i=SHIFTX1DN;i<=N1-1+SHIFT1+SHIFTX1UP;i++)
#define MYOCOMPLOOPF MYOCOMPLOOPF3 MYOCOMPLOOPF2 MYOCOMPLOOPF1

  MYOCOMPLOOPF{
    PLOOPBONLY(pl){
      get_geometry(i, j, k, FACE1+pl-B1, ptrgeomf);

      // first convert staggered upoint to pstag
      MACP0A1(pstag,i,j,k,pl)=MACP0A1(upoint,i,j,k,pl)*sign(ptrgeomf->gdet)/(fabs(ptrgeomf->gdet)+SMALL);
      // GODMARK: Some problem with gdet not being 0 on axis!!  Leads to pstag[B2]=-2E50, but apparently overwritten?
    }
  }

  //  COMPZLOOP{
  //    dualfprintf(fail_file,"DEATH1: i=%d j=%d k=%d :: %21.15g %21.15g %21.15g\n",i,j,k,MACP0A1(pstag,i,j,k,B1),MACP0A1(pstag,i,j,k,B2),MACP0A1(pstag,i,j,k,B3));
  //  }

  // bound new pstag
  if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;
  bound_pstag(stage, finalstep, boundtime, preal, pstag, upoint, USEMPI);


  //  COMPFULLLOOP{
  //    dualfprintf(fail_file,"DEATH2: i=%d j=%d k=%d :: %21.15g %21.15g %21.15g\n",i,j,k,MACP0A1(pstag,i,j,k,B1),MACP0A1(pstag,i,j,k,B2),MACP0A1(pstag,i,j,k,B3));
  //  }

  // now can define everything from pstag
  // GODMARK: this defines pcent and upoint on larger domain than normal code, so may mask problem.
  //  COMPFULLLOOP{

  // correctly similar constrained loop as in normal code:
  int is,ie,js,je,ks,ke;
  get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);
  COMPZSLOOP(is,ie,js,je,ks,ke){

    get_geometry(i, j, k, CENT, ptrgeomc);
    PLOOPBONLY(pl){
      get_geometry(i, j, k, FACE1+pl-B1, ptrgeomf);
      if(pl==B1){
        ii=ip1mac(i);
        jj=j;
        kk=k;
      }
      else if(pl==B2){
        ii=i;
        jj=jp1mac(j);
        kk=k;
      }
      else if(pl==B3){
        ii=i;
        jj=j;
        kk=kp1mac(k);
      }

      if(ii==i && jj==j && kk==k){
        // then just copy
        MACP0A1(pcent,i,j,k,pl)=MACP0A1(pstag,i,j,k,pl);
      }
      else{
        // now "interpolate" pstag -> pcent
        get_geometry(ii, jj, kk, FACE1+pl-B1, ptrgeomfu);
        // below consistent with interpolate_ustag2fieldcent(), but can try different way 
        MACP0A1(pcent,i,j,k,pl)=0.5*( MACP0A1(pstag,i,j,k,pl)*(ptrgeomf->gdet/ptrgeomc->gdet) + MACP0A1(pstag,ii,jj,kk,pl)*(ptrgeomfu->gdet/ptrgeomc->gdet));
      }

      // finally get conserved quantity at CENT for inversion
      MACP0A1(upoint,i,j,k,pl)=MACP0A1(pcent,i,j,k,pl)*(ptrgeomc->gdet);
    }// end over pl
  }// end loop

  return(0);
}
