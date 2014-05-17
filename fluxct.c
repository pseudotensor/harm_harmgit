
#include "decs.h"

/*! \file fluxct.c
  \brief TOTH CT method for preserving divb=0 (FLUXB==FLUXCTTOTH)

  ////////////////////////////////
  //
  // Notes on sign conventions:
  //
  ///////////////////////////////

  // flux part is just average of same emf term at 4 different edge locations of (B^2 v^1 - B^1 v^2)
  //    COMPEMFZLOOP{
  //COMPCOMPLOOPINFP1{ // constrain or control better? GODMARK
  // B^i = \dF^{it}
  // E_i = - [ijk] v^j B^k  , such that (\detg B^i),t = - (\detg(B^i v^j - B^j v^i)),j = - (\detg [ijk] E_k),j = ([ijk] emf[k]),j
      
  // -> E_1 = v^3 B^2 - v^2 B^3
  // -> E_2 = v^1 B^3 - v^3 B^1
  // -> E_3 = v^2 B^1 - v^1 B^2

  // emf[i] = - \detg E_i

  // And notice that Fj[Bi] = \dF^{ij} = B^i v^j - B^j v^i , where j=dir

  // so:
  // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
  // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
  // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]

  // Notice only 6 independent ways.  The diagonal terms vanish (e.g. Fi[Bi]=0).
  */


/// compute field at CENT from vector potential A at CORN1,2,3
/// assumes normal field p
int vpot2field_centeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR])
{
  int Nvec[NDIM];

  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  /* flux-ct */

  // A[1] located at CORN1
  // A[2] located at CORN2
  // A[3] located at CORN3

  // F_{\mu\nu} \equiv A_{\nu,\mu} - A_{\mu,\nu}

  // B^i \equiv \dF^{it}

  // F_{ij} = \detg B^k [ijk]

  // F_{\theta\phi} = \detg B^r

  // F_{\phi r} = \detg B^\theta

  // F_{r\theta} = \detg B^\phi


  // \detg B^x = A_{z,y} - A_{y,z}
  // \detg B^y = A_{x,z} - A_{z,x}
  // \detg B^z = A_{y,x} - A_{x,y}

  // loop optimized for filling cells rather than for speed

      
  // can do full loop since A[1,2,3] are defined even on their plus parts (redundant but simpler than messing with B's post-hoc by linear extrapolation or something)
  // puts burden on computing A[1,2,3] (such as having radius or anything at plus part)
  // this burden is easier since coord() and bl_coord() have no limits on the i,j,k they are given.
  // so in principle one can give an analytic expression
  // however, depends on metric and such things if converting expression from one basis to another.
  // however, this way is more correct than post-hoc extrapolation.
  // only an issue for fixed boundary conditions, where one could just specify the flux directly.
  // ufield explicitly needed for FV method


  // since nowait'ed above (because each loop sets ufield,pfield for each dir that don't depend upon eachother), need barrier here at end to stop before continuing (implied in parallel reigon)
#pragma omp parallel 
  {
    int i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE igdetgnosing;
    int dir;
    int odir1,odir2;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL; // doesn't depend upon dir, but blockijk must be private, so put in parallel loop


    // generally ptr's are different inside parallel block
    ptrgeom=&geomdontuse;


    /////////////////
    //
    // B^1
    //
    /////////////////

    dir=1;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ////    COMPFULLLOOP{ // COMPFULLLOOP allows since A_i exists at COMPFULLLOOPP1 and so always accessing valid A_i

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) nowait // nowait valid because each next loop writes to independent memory regions (or to local temporary variables like igdetgnosing that is overwritten for each iteration and so doesn't matter).  And also don't require result of one loop for next loop (this is often not true!)
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // ufield doesn't require geometry
        MACP0A1(ufield,i,j,k,B1)  = +(AVGCORN_1(A[3],i,jp1mac(j),k)-AVGCORN_1(A[3],i,j,k))/(dx[2]);
        MACP0A1(ufield,i,j,k,B1) += -(AVGCORN_1(A[2],i,j,kp1mac(k))-AVGCORN_1(A[2],i,j,k))/(dx[3]);

        get_geometry(i, j, k, CENT, ptrgeom);
        igdetgnosing = sign(ptrgeom->gdet)/(fabs(ptrgeom->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing;
      }
    }// end if doing this dir


    /////////////////
    //
    // B^2
    //
    /////////////////
    
    dir=2;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ////    COMPFULLLOOP{

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) nowait
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        MACP0A1(ufield,i,j,k,B2)  = +(AVGCORN_2(A[1],i,j,kp1mac(k))-AVGCORN_2(A[1],i,j,k))/(dx[3]);
        MACP0A1(ufield,i,j,k,B2) += -(AVGCORN_2(A[3],ip1mac(i),j,k)-AVGCORN_2(A[3],i,j,k))/(dx[1]);

        get_geometry(i, j, k, CENT, ptrgeom);
        igdetgnosing = sign(ptrgeom->gdet)/(fabs(ptrgeom->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing;
      }
    } // end if dir
  

    /////////////////
    //
    // B^3
    //
    /////////////////

    dir=3;
    get_odirs(dir,&odir1,&odir2);
    if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
      ////    COMPFULLLOOP{

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) nowait
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        MACP0A1(ufield,i,j,k,B3)  = +(AVGCORN_3(A[2],ip1mac(i),j,k)-AVGCORN_3(A[2],i,j,k))/(dx[1]);
        MACP0A1(ufield,i,j,k,B3) += -(AVGCORN_3(A[1],i,jp1mac(j),k)-AVGCORN_3(A[1],i,j,k))/(dx[2]);

        get_geometry(i, j, k, CENT, ptrgeom);
        igdetgnosing = sign(ptrgeom->gdet)/(fabs(ptrgeom->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
        MACP0A1(pfield,i,j,k,B1-1+dir)  = MACP0A1(ufield,i,j,k,B1-1+dir)*igdetgnosing;
      }
    }// end if dir
  
  }// end full parallel region (with implied barrier)





  

  return(0);
}





/// compute flux for FLUXCTTOTH method
int flux_ct(int stage,
            int initialstep, int finalstep,
            FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int *Nvec, FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime)
{
  int flux_ct_computeemf(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int flux_ct_diffusivecorrections(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int flux_ct_emf2flux(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);






  MYFUN(flux_ct_computeemf(stage, pb, emf, vconemf, dq1, dq2, dq3, F1, F2, F3),"step_ch.c:advance()", "flux_ct",1);

  MYFUN(flux_ct_diffusivecorrections(stage, pb, emf, vconemf, dq1, dq2, dq3, F1, F2, F3),"step_ch.c:advance()", "flux_ct",1);


  if(EVOLVEWITHVPOT>0 ||  TRACKVPOT>0){
    // Evolve A_i
    // Had to be here where EMFs are at standard CORN1,2,3 positions and before final F assigned
    // TOTH CT method doesn't cleanly differentiate between point update and average update of A_i, so just stick to TOTH CT EMF itself
    evolve_vpotgeneral(FLUXB, stage, initialstep, finalstep, pb, Nvec, NULL, emf, CUf, CUnew, fluxdt, fluxtime, vpot);
  }

  int fluxvpot_modifyemfsuser=0;
  fluxvpot_modifyemfsuser=(EVOLVEWITHVPOT>0 ||  TRACKVPOT>0)&&(MODIFYEMFORVPOT==MODIFYEMF || MODIFYEMFORVPOT==MODIFYVPOT);

  if(fluxvpot_modifyemfsuser==0){// if didn't already call adjust_emfs() in fluxvpot above, have to allow user to be able to still modify emfs calling function directly
    // User "boundary conditions" to modify EMFs before used to get fluxes
    adjust_fluxcttoth_emfs(fluxtime,pb,emf);
  }

  //////////////
  //
  // must come last for FLUXCTTOTH to produce correct fluxes from point corner EMFs
  //
  //////////////
  MYFUN(flux_ct_emf2flux(stage, pb, emf, vconemf, dq1, dq2, dq3, F1, F2, F3),"step_ch.c:advance()", "flux_ct",1);



  return(0);
}






/// Compute TOTH EMF
/// OPENMPMARK: Apparently the below loops are expensive due to OpenMP overhead.  Using static schedule helps overhead a bit
int flux_ct_computeemf(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  // full-type geometry below
  FTYPE coefemf[NDIM];
  extern int choose_limiter(int dir, int i, int j, int k, int pl);


  // Note that Toth method needs off-direction flux beyond where normal-flux needed.  For example, for dir=1 setting F1(i) needs F23[i-1],F23[i].  So fluxloop needs to define F2/F3 as if cell centered quantity.
  // Also, if modify flux[B1,B2,B3] after set by flux calculation, then that information can propagate from outer boundaries to active region
  // When FLUXCTTOTH method used, must not modify fluxes within ghost region (say setting F1[B2]=-F2[B1]) since info will move to active region
  // e.g. if put NaN in EMF3 at very outer edge, then reaches active domain by this method




  if(FLUXB==FLUXCTHLL){
    dualfprintf(fail_file,"Makes no sense to call flux_ct with FLUXB==FLUXCTHLL\n");
    myexit(9176325);
  }




  ///////////////////////
  //
  //  COMPUTE PRE-FUNCTIONS used by Athena method
  //
  ///////////////////////
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)){
    // compute v^i

    // loop must go over COMPEMFZLOOP's range minus 1 for i and j and k
    //    COMPPREEMFZLOOP{
    //    COMPFULLLOOP{ // GODMARK: could try to be more constrained
    // use ucon_calc() to get v^i 
#pragma omp parallel 
    {
      int i, j, k, l;
      struct of_geom geomcfulldontuse;
      struct of_geom *ptrgeomcfull=&geomcfulldontuse;
      FTYPE ucon[NDIM];
      FTYPE others[NUMOTHERSTATERESULTS];

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
    

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        // EMF below is based upon averaging of zone-centered quantities, so use CENT here (i.e. not CORN)
        get_geometry(i, j, k, CENT, ptrgeomcfull);
        MYFUN(ucon_calc(MAC(pb,i,j,k), ptrgeomcfull, ucon, others),"fluxct.c:flux_ct()", "ucon_calc() dir=0", 1);


        // ptrgeom->gdet is \detg for EMF flux (ptrgeom->e is EOM factor for flux equation)
#if(CORNGDETVERSION)
        for(l=U1;l<=U3;l++) MACP0A1(vconemf,i,j,k,l)=(ucon[l-U1+1]/ucon[TT]); // put in at end
#else
        for(l=U1;l<=U3;l++) MACP0A1(vconemf,i,j,k,l)=(ucon[l-U1+1]/ucon[TT])*(ptrgeomcfull->EOMFUNCMAC(l));
#endif
      }// end 3D LOOP
    }// end parallel region
  }// end if ATHENA1||ATHENA2










#if(CORNGDETVERSION)

  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)||(FLUXB==FLUXCD)){

    /////////////
    //
    // strip off geometry factor, which is added when computing the EMF
    //
    ////////////
    //    COMPFULLLOOP{ // GODMARK: could try to be more constrained
#pragma omp parallel 
    {
      int i,j,k;
      // gdet-type geometry below
      struct of_gdetgeom geomf1dontuse,geomf2dontuse,geomf3dontuse;
      struct of_gdetgeom *ptrgeomf1=&geomf1dontuse,*ptrgeomf2=&geomf2dontuse,*ptrgeomf3=&geomf3dontuse;

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
    

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        ////////////////////
        // F1
        ////////////////////
#if(N1>1)
        get_geometry_gdetmix(i,j,k,FACE1,ptrgeomf1);
        MACP0A1(F1,i,j,k,B1)*=(ptrgeomf1->IEOMFUNCNOSINGMAC(B1));
        MACP0A1(F1,i,j,k,B2)*=(ptrgeomf1->IEOMFUNCNOSINGMAC(B2));
        MACP0A1(F1,i,j,k,B3)*=(ptrgeomf1->IEOMFUNCNOSINGMAC(B3));
#endif

        ////////////////////
        // F2
        ////////////////////
#if(N2>1)
        get_geometry_gdetmix(i,j,k,FACE2,ptrgeomf2);
        MACP0A1(F2,i,j,k,B1)*=(ptrgeomf2->IEOMFUNCNOSINGMAC(B1));
        MACP0A1(F2,i,j,k,B2)*=(ptrgeomf2->IEOMFUNCNOSINGMAC(B2));
        MACP0A1(F2,i,j,k,B3)*=(ptrgeomf2->IEOMFUNCNOSINGMAC(B3));
#endif
    
        ////////////////////
        // F3
        ////////////////////
#if(N3>1)
        get_geometry_gdetmix(i,j,k,FACE3,ptrgeomf3);
        MACP0A1(F3,i,j,k,B1)*=(ptrgeomf3->IEOMFUNCNOSINGMAC(B1));
        MACP0A1(F3,i,j,k,B2)*=(ptrgeomf3->IEOMFUNCNOSINGMAC(B2));
        MACP0A1(F3,i,j,k,B3)*=(ptrgeomf3->IEOMFUNCNOSINGMAC(B3));
#endif
      }// end 3D loop
    }// end parallel region

    // don't need to put back geometry factor since in the end the EMF defines the Flux completely (except for FLUXB==FLUXCTHLL)
  }// end if ATHENA1||ATHENA2||FLUXCTTOTH||FLUXCD
#endif // end if CORNGDETVERSION






  // GODMARK: strange one has to do this.  Related to FLUXCT not reducing correctly for plane-parallel grid-aligned flows?
  if((N2>1)&&(N3>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; // or 0 if neither as controlled by below
  if((N1>1)&&(N3>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; // or 0 if neither as controlled by below
  if((N1>1)&&(N2>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; // or 0 if neither as controlled by below





  /////////////////////////
  //
  // COMPUTE EMF
  //
  //////////////////////////  

  if((FLUXB==FLUXCTTOTH)||(FLUXB==ATHENA1)||(FLUXB==ATHENA2)){
    /* calculate EMFs */
    /* Toth approach: just average */


    // emf_i centered on corner of plane with normal direction i

    // Note, do need F1[-N1BND] but don't need F2[-N2BND] or F3[-N3BND] in averaging process, which makes normally need to be complicated.
    // Note: As discussed in superdefs.h, superdefs.pointers.h, and set_arrays_multidimen.c, F1,F2,F3,F1EM,F2EM,F3EM pointers are normal arrays, but their BASE arrays have extra memory at *bottom* so can access F[-NBND-1] without segfaulting.  Resulting data not used, but makes loops simpler.


#pragma omp parallel  // globalprivate stuff needed for final CORNGDET setting if doing that method
    {
      int i,j,k;
      // gdet-type geometry below
      struct of_gdetgeom geomf1dontuse,geomf2dontuse,geomf3dontuse;
      struct of_gdetgeom *ptrgeomf1=&geomf1dontuse,*ptrgeomf2=&geomf2dontuse,*ptrgeomf3=&geomf3dontuse;
      OPENMP3DLOOPVARSDEFINE;
      

      //    //    COMPLOOPINFP1dir1full
      //      OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,INFULLP12,OUTFULL2,INFULLP13,OUTFULL3);
      ////    COMPLOOPINFP1dir2full
      //// OPENMP3DLOOPSETUP(INFULLP11,OUTFULL1,-N2BND,N2-1+N2BND,INFULLP13,OUTFULL3);
      ////    COMPLOOPINFP1dir3full
      //      OPENMP3DLOOPSETUP(INFULLP11,OUTFULL1,INFULLP12,OUTFULL2,-N3BND,N3-1+N3BND);
      
      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) ///nowait // Can use "nowait" because each emf[{1,2,3}] set independently in each successive loop
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);
 

        ////////////////////
        // EMF1
        ////////////////////
#if((N2>1)||(N3>1))
        MACP1A0(emf,1,i,j,k) =
          coefemf[1] * (
#if(N2>1)
                        + MACP0A1(F2,i,j,k,B3) + MACP0A1(F2,i,j,km1mac(k),B3)
#endif
#if(N3>1)
                        - MACP0A1(F3,i,j,k,B2) - MACP0A1(F3,i,jm1mac(j),k,B2)
#endif
                        );
#else // end if doing EMF1
        MACP1A0(emf,1,i,j,k)=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif // end if not doing EMF1

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CORN1,ptrgeomf1);
        // obviously ptrgeom->EOMFUNCMAC(B2) has to be equal to ptrgeom->EOMFUNCMAC(B3) for this method
        MACP1A0(emf,1,i,j,k) *= (ptrgeomf1->EOMFUNCMAC(B2));
#endif // end if CORNGDETVERSION



        ////////////////////
        // EMF2
        ////////////////////
#if((N1>1)||(N3>1))
        MACP1A0(emf,2,i,j,k) =
          coefemf[2] * (
#if(N3>1)
                        + MACP0A1(F3,i,j,k,B1) + MACP0A1(F3,im1mac(i),j,k,B1)
#endif
#if(N1>1)
                        - MACP0A1(F1,i,j,k,B3) - MACP0A1(F1,i,j,km1mac(k),B3)
#endif
                        );
#else // end if doing EMF2
        MACP1A0(emf,2,i,j,k)=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif // end if not doing EMF2

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CORN2,ptrgeomf2);
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        MACP1A0(emf,2,i,j,k) *=(ptrgeomf2->EOMFUNCMAC(B1));
#endif // end if CORNGDETVERSION




        ////////////////////
        // EMF3
        ////////////////////
#if((N1>1)||(N2>1))
        MACP1A0(emf,3,i,j,k) =
          coefemf[3] * (
#if(N1>1)
                        + MACP0A1(F1,i,j,k,B2) + MACP0A1(F1,i,jm1mac(j),k,B2)
#endif
#if(N2>1)
                        - MACP0A1(F2,i,j,k,B1) - MACP0A1(F2,im1mac(i),j,k,B1)
#endif
                        );
#else // end if doing EMF3
        MACP1A0(emf,3,i,j,k)=0.0; // not really 0, but differences in emf will be 0, and that's all that matters
#endif // end if not doing EMF3

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CORN3,ptrgeomf3);
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        MACP1A0(emf,3,i,j,k) *=(ptrgeomf3->EOMFUNCMAC(B1));
#endif // end if CORNGDETVERSION

      }// end 3D LOOP
    }// end parallel region (and implied barrier)









  }// end if FLUXCT TOTH
  else if(FLUXB==FLUXCD){




    // Centered Difference (CD)
    // here emf is actually electric field
    // in Toth 2000, the F1=f^x F2=f^y
    // see Toth 2000 equation 28/29
    // 0.125 here takes care of larger differencing used in advance() in advance.c
    // all emf's end up at CENT
    //    COMPEMFZLOOP{
    //    COMPLOOPOUTFM1{ // constrain or control better? GODMARK



    // GODMARK: Why did Charles change sign in front of F's (see old code) ?  He changed it as if real sign such that emf[i] = \detg E_i but then later flux was wrong sign since he set flux=emf

    // GODMARK: should just compute F with diffusive term at CENT directly instead of averaging flux

#pragma omp parallel  // globalprivate stuff needed for final CORNGDET setting if doing that method
    {
      int i,j,k;
      struct of_gdetgeom geomcdontuse;
      struct of_gdetgeom *ptrgeomc=&geomcdontuse;
      OPENMP3DLOOPVARSDEFINE; 


      ////    COMPLOOPOUTFM1dir1full
      //OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,INFULL2,OUTFULLM12,INFULL3,OUTFULLM13);
      ////    COMPLOOPOUTFM1dir2full
      ////      OPENMP3DLOOPSETUP(INFULL1,OUTFULLM11,-N2BND,N2-1+N2BND,INFULL3,OUTFULLM13);
      ////    COMPLOOPOUTFM1dir3full
      ////OPENMP3DLOOPSETUP(INFULL1,OUTFULLM11,INFULL2,OUTFULLM12,-N3BND,N3-1+N3BND);

      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) /////nowait // Can use "nowait" since each emf[{1,2,3}] set independently in each successive loop
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        ////////////////////
        // EMF1
        ////////////////////
#if((N2>1)||(N3>1))
        MACP1A0(emf,1,i,j,k) =
          0.125 * (
                   + MACP0A1(F2,i,j,k,B3) + MACP0A1(F2,i,jp1mac(j),k,B3)
                   - MACP0A1(F3,i,j,k,B2) - MACP0A1(F3,i,j,kp1mac(k),B2)
                   );

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CENT,ptrgeomc);
        // obviously ptrgeom->EOMFUNCMAC(B2) has to be equal to ptrgeom->EOMFUNCMAC(B3) for this method
        MACP1A0(emf,1,i,j,k) *=(ptrgeomc->EOMFUNCMAC(B2));
#endif // end if CORNGDETVERSION
 
#endif // end if doing EMF1



        ////////////////////
        // EMF2
        ////////////////////
#if((N1>1)||(N3>1))

        MACP1A0(emf,2,i,j,k) =
          0.125 * (
                   + MACP0A1(F3,i,j,k,B1) + MACP0A1(F3,i,j,kp1mac(k),B1)
                   - MACP0A1(F1,i,j,k,B3) - MACP0A1(F1,ip1mac(i),j,k,B3)
                   );

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CENT,ptrgeomc);
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        MACP1A0(emf,2,i,j,k) *=(ptrgeomc->EOMFUNCMAC(B1));
#endif // end if CORNGDETVERSION

#endif // end if doing EMF2



        ////////////////////
        // EMF3
        ////////////////////
#if((N1>1)||(N2>1))
        MACP1A0(emf,3,i,j,k) =
          0.125 * (
                   + MACP0A1(F1,i,j,k,B2) + MACP0A1(F1,ip1mac(i),j,k,B2)
                   - MACP0A1(F2,i,j,k,B1) - MACP0A1(F2,i,jp1mac(j),k,B1)
                   );

#if(CORNGDETVERSION)// then tack on geometry
        get_geometry_gdetmix(i,j,k,CENT,ptrgeomc);
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        MACP1A0(emf,3,i,j,k) *=(ptrgeomc->EOMFUNCMAC(B1));
#endif // end if CORNGDETVERSION

#endif // end if doing EMF3

      }// end 3D LOOP

    }// end parallel region (and implied barrier)




  }// end if FLUXCD



  return(0);
}

















/// Compute "diffusive corrections" to avoid field loop boost type issue
/// OPENMPMARK: Apparently the below loops are expensive due to OpenMP overhead.  Using static schedule helps overhead a bit
int flux_ct_diffusivecorrections(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  // full-type geometry below
  FTYPE coefemf[NDIM];
  extern int choose_limiter(int dir, int i, int j, int k, int pl);







  ///////////////////////////////////////////////////////
  //
  // ADD DIFFUSIVE CORRECTIONS
  //
  ///////////////////////////////////////////////////////




  // add diffusive term
  // must come after FLUXCT
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)){

    // Stone & Gardiner point out that Toth FLUXCT and FLUXCD not consistent with underlying integration algorithm for plane-parallel, grid-aligned flows.
    // fix is to change 0.25 to 0.5 and use a diffusive term as in ATHENA1

    /* Stone & Gardiner eq. 39 */
    // Charles Gammie (8/17/05) says this is simple, but repairs HARM defect that 2D does not reduce to 1D when waves are along coordinate lines


    ////    COMPLOOPINFP1 // constrain or control better? GODMARK
    //    COMPEMFZLOOP{
#pragma omp parallel 
    {
      int i,j,k,l;
      // gdet-type geometry below
      struct of_gdetgeom geomf1dontuse,geomf2dontuse,geomf3dontuse;
      struct of_gdetgeom *ptrgeomf1=&geomf1dontuse,*ptrgeomf2=&geomf2dontuse,*ptrgeomf3=&geomf3dontuse;
      FTYPE diffusiveterm[NDIM];
      FTYPE emfmm[NDIM],emfpm[NDIM],emfmp[NDIM],emfpp[NDIM],alpha[NDIM] ;


      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(INFULLP11,OUTFULL1,INFULLP12,OUTFULL2,INFULLP13,OUTFULL3);

      // generally ptr's are different inside parallel block
      ptrgeomf1=&geomf1dontuse;
      ptrgeomf2=&geomf2dontuse;
      ptrgeomf3=&geomf3dontuse;


#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // {emf}_i=-\epsilon_{ijk} v^i B^k


        // average of the below results in averaged emf located at corner (CORN)
        // same sign as F's (B^2 v^1 - B^1 v^2) with gdet built into vconemf

        // could remove gdet from vconemf and F1/F2 and put gdet at CORN for final result!  Avoids axis problems?
        // applies to above Toth version as well!
        // GODMARK

#if((N2>1)||(N3>1))
        // emf_1
        emfmp[1] = 
          MACP0A1(pb,i  ,jm1mac(j),k  ,B3)*MACP0A1(vconemf,i  ,jm1mac(j),k  ,U2) -
          MACP0A1(pb,i  ,jm1mac(j),k  ,B2)*MACP0A1(vconemf,i  ,jm1mac(j),k  ,U3) ;
        emfmm[1] = 
          MACP0A1(pb,i  ,jm1mac(j),km1mac(k),B3)*MACP0A1(vconemf,i  ,jm1mac(j),km1mac(k),U2) -
          MACP0A1(pb,i  ,jm1mac(j),km1mac(k),B2)*MACP0A1(vconemf,i  ,jm1mac(j),km1mac(k),U3) ;
        emfpm[1] = 
          MACP0A1(pb,i  ,j  ,km1mac(k),B3)*MACP0A1(vconemf,i  ,j  ,km1mac(k),U2) -
          MACP0A1(pb,i  ,j  ,km1mac(k),B2)*MACP0A1(vconemf,i  ,j  ,km1mac(k),U3) ;
        emfpp[1] = 
          MACP0A1(pb,i  ,j  ,k  ,B3)*MACP0A1(vconemf,i  ,j  ,k  ,U2) -
          MACP0A1(pb,i  ,j  ,k  ,B2)*MACP0A1(vconemf,i  ,j  ,k  ,U3) ;
#endif

#if((N1>1)||(N3>1))
        // emf_2
        emfmp[2] = 
          MACP0A1(pb,im1mac(i),j  ,k  ,B1)*MACP0A1(vconemf,im1mac(i),j  ,k  ,U3) -
          MACP0A1(pb,im1mac(i),j  ,k  ,B3)*MACP0A1(vconemf,im1mac(i),j  ,k  ,U1) ;
        emfmm[2] = 
          MACP0A1(pb,im1mac(i),j  ,km1mac(k),B1)*MACP0A1(vconemf,im1mac(i),j  ,km1mac(k),U3) -
          MACP0A1(pb,im1mac(i),j  ,km1mac(k),B3)*MACP0A1(vconemf,im1mac(i),j  ,km1mac(k),U1) ;
        emfpm[2] = 
          MACP0A1(pb,i  ,j  ,km1mac(k),B1)*MACP0A1(vconemf,i  ,j  ,km1mac(k),U3) -
          MACP0A1(pb,i  ,j  ,km1mac(k),B3)*MACP0A1(vconemf,i  ,j  ,km1mac(k),U1) ;
        emfpp[2] = 
          MACP0A1(pb,i  ,j  ,k  ,B1)*MACP0A1(vconemf,i  ,j  ,k  ,U3) -
          MACP0A1(pb,i  ,j  ,k  ,B3)*MACP0A1(vconemf,i  ,j  ,k  ,U1) ;
#endif

#if((N1>1)||(N2>1))
        // emf_3 (same sign as FLUXCT method as emf[3])
        emfmp[3] = 
          MACP0A1(pb,im1mac(i),j  ,k  ,B2)*MACP0A1(vconemf,im1mac(i),j  ,k  ,U1) -
          MACP0A1(pb,im1mac(i),j  ,k  ,B1)*MACP0A1(vconemf,im1mac(i),j  ,k  ,U2) ;
        emfmm[3] = 
          MACP0A1(pb,im1mac(i),jm1mac(j),k  ,B2)*MACP0A1(vconemf,im1mac(i),jm1mac(j),k  ,U1) -
          MACP0A1(pb,im1mac(i),jm1mac(j),k  ,B1)*MACP0A1(vconemf,im1mac(i),jm1mac(j),k  ,U2) ;
        emfpm[3] = 
          MACP0A1(pb,i  ,jm1mac(j),k  ,B2)*MACP0A1(vconemf,i  ,jm1mac(j),k  ,U1) -
          MACP0A1(pb,i  ,jm1mac(j),k  ,B1)*MACP0A1(vconemf,i  ,jm1mac(j),k  ,U2) ;
        emfpp[3] = 
          MACP0A1(pb,i  ,j  ,k  ,B2)*MACP0A1(vconemf,i  ,j  ,k  ,U1) -
          MACP0A1(pb,i  ,j  ,k  ,B1)*MACP0A1(vconemf,i  ,j  ,k  ,U2) ;
#endif

        for(l=1;l<=3;l++){
          diffusiveterm[l]= 0.25*(emfmp[l] + emfmm[l] + emfpm[l] + emfpp[l]);
        }

#if(CORNGDETVERSION)// then tack on geometry
      
        get_geometry_gdetmix(i,j,k,CORN1,ptrgeomf1);
        get_geometry_gdetmix(i,j,k,CORN2,ptrgeomf2);
        get_geometry_gdetmix(i,j,k,CORN3,ptrgeomf3);
      
        // obviously ptrgeom->EOMFUNCMAC(B2) has to be equal to ptrgeom->EOMFUNCMAC(B3) for this method
        diffusiveterm[1] *=(ptrgeomf1->EOMFUNCMAC(B2));
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        diffusiveterm[2] *=(ptrgeomf2->EOMFUNCMAC(B1));
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        diffusiveterm[3] *=(ptrgeomf3->EOMFUNCMAC(B1));
#endif

        // now add diffusive term to emf
        // notice original emf multiplied by 2 to account for diffusive term being subtracted, so result is consistent
        for(l=1;l<=3;l++){
          MACP1A0(emf,l,i,j,k) = 2.0*MACP1A0(emf,l,i,j,k) - diffusiveterm[l];
        }


      }// end COMPEMFZLOOP
    }// end parallel region

  }// end ATHENA1
      









  // add another diffusive flux correction
  // must come after ATHENA1
  if(FLUXB==ATHENA2){



    if(LIMADJUST>0 || DODQMEMORY==0 || choose_limiter(1, 0,0,0,RHO)>=PARA || choose_limiter(2, 0,0,0,RHO)>=PARA || choose_limiter(3, 0,0,0,RHO)>=PARA){ // note only look at one point and one variable here for test
      dualfprintf(fail_file,"Cannot use Athena2 with limadjust since dq's not defined\n");
      myexit(11);
    }

    // GODMARK
    // should just use unrescaled p2interp stored, but need one for each direction.
    // should also use wave speeds from edges?




    /* Stone & Gardiner eq. 48 */
    // Charles Gammie (8/17/05) says does much better on flux loop advection test than ordinary HARM
    //    COMPLOOPINFP1{ // constrain or control better? GODMARK
    //    COMPEMFZLOOP{
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // requires full copyin()
    {
      // Below stuff is for Athena 1 and Athena 2
      int dir;
      int i,j,k,l,pl,pliter;
      // gdet-type geometry below
      struct of_gdetgeom geomf1dontuse,geomf2dontuse,geomf3dontuse;
      struct of_gdetgeom *ptrgeomf1=&geomf1dontuse,*ptrgeomf2=&geomf2dontuse,*ptrgeomf3=&geomf3dontuse;
      FTYPE diffusiveterm[NDIM];
      FTYPE emfmm[NDIM],emfpm[NDIM],emfmp[NDIM],emfpp[NDIM],alpha[NDIM] ;
      // Gammie stuff
      FTYPE B1pp,B1pm,B1mp,B1mm;
      FTYPE B2pp,B2pm,B2mp,B2mm ;
      FTYPE B3pp,B3pm,B3mp,B3mm ;
      FTYPE U1pp,U1pm,U1mp,U1mm;
      FTYPE U2pp,U2pm,U2mp,U2mm ;
      FTYPE U3pp,U3pm,U3mp,U3mm ;
      //  FTYPE cms_func(FTYPE *prim_var) ;
      FTYPE B1d,B1u,B1l,B1r;
      FTYPE B2d,B2u,B2l,B2r;
      FTYPE B3d,B3u,B3l,B3r;
      FTYPE pbavg[NPR];
      struct of_state state;
      int ignorecourant;
      FTYPE cmax1,cmin1,cmax2,cmin2,ctop1,ctop2;
      struct of_geom geomco1dontuse,geomco2dontuse,geomco3dontuse;
      struct of_geom *ptrgeomco1=&geomco1dontuse,*ptrgeomco2=&geomco2dontuse,*ptrgeomco3=&geomco3dontuse;


      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(INFULLP11,OUTFULL1,INFULLP12,OUTFULL2,INFULLP13,OUTFULL3);


      // generally ptr's are different inside parallel block
      ptrgeomf1=&geomf1dontuse;
      ptrgeomf2=&geomf2dontuse;
      ptrgeomf3=&geomf3dontuse;
      ptrgeomco1=&geomco1dontuse;
      ptrgeomco2=&geomco2dontuse;
      ptrgeomco3=&geomco3dontuse;

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        // dq1 and dq2 and dq3 are well-defined from fluxcalc() (both directions)
     

        // simple average of 1-D linear extrapolation here, where could use p2interp data saved from step_ch.c?  still available by this function call?


#if((N2>1)||(N3>1))
        // for emf_1

        // B2 located at FACE2 @ (k-1)
        B2d = 0.5*(
                   MACP0A1(pb,i,jm1mac(j),km1mac(k),B2) + 0.5*MACP0A1(dq2,i,jm1mac(j),km1mac(k),B2) +
                   MACP0A1(pb,i,j  ,km1mac(k),B2) - 0.5*MACP0A1(dq2,i,j  ,km1mac(k),B2)
                   ) ;

        // B2 located at FACE2 @ (k)
        B2u = 0.5*(
                   MACP0A1(pb,i,jm1mac(j),k,B2) + 0.5*MACP0A1(dq2,i,jm1mac(j),k,B2) +
                   MACP0A1(pb,i,j  ,k,B2) - 0.5*MACP0A1(dq2,i,j  ,k,B2)
                   ) ;

        // B3 located at FACE3 @ (j-1)
        B3l = 0.5*(
                   MACP0A1(pb,i,jm1mac(j),km1mac(k),B3) + 0.5*MACP0A1(dq3,i,jm1mac(j),km1mac(k),B3) +
                   MACP0A1(pb,i,jm1mac(j),k  ,B3) - 0.5*MACP0A1(dq3,i,jm1mac(j),k  ,B3)
                   ) ;
        // B3 located at FACE3 @ j
        B3r = 0.5*(
                   MACP0A1(pb,i,j,km1mac(k),B3) + 0.5*MACP0A1(dq3,i,j,km1mac(k),B3) +
                   MACP0A1(pb,i,j,k  ,B3) - 0.5*MACP0A1(dq3,i,j,k  ,B3)
                   ) ;

        // B2 for all centers around CORN1
        // B2[j - mp][k - mp]
        B2mm = MACP0A1(pb,i,jm1mac(j),km1mac(k),B2) ;
        B2mp = MACP0A1(pb,i,jm1mac(j),k  ,B2) ;
        B2pm = MACP0A1(pb,i,j  ,km1mac(k),B2) ;
        B2pp = MACP0A1(pb,i,j  ,k  ,B2) ;

        // B3 for all centers around CORN1
        // B3[j - mp][k - mp]
        B3mm = MACP0A1(pb,i,jm1mac(j),km1mac(k),B3) ;
        B3mp = MACP0A1(pb,i,jm1mac(j),k  ,B3) ;
        B3pm = MACP0A1(pb,i,j  ,km1mac(k),B3) ;
        B3pp = MACP0A1(pb,i,j  ,k  ,B3) ;

        // compute characteristic velocity -- only for Athena2 method


        // average pb to CORN1 for average phase speed there
        PLOOP(pliter,pl) pbavg[pl]=0.25*(MACP0A1(pb,i,j,k,pl)+MACP0A1(pb,i,jm1mac(j),k,pl)+MACP0A1(pb,i,j,km1mac(k),pl)+MACP0A1(pb,i,jm1mac(j),km1mac(k),pl));
        get_geometry(i, j, k, CORN1, ptrgeomco1); // used here and below emf's
        MYFUN(get_state(pbavg, ptrgeomco1, &state),"step_ch.c:flux_ct()", "get_state()", 1);
        // KORALNOTE: Correct because this applies to the diffusive term only for the magnetic field evolution, so needs to stay vchar()
        dir=2; MYFUN(vchar(pbavg, &state, dir, ptrgeomco1, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
        dir=3; MYFUN(vchar(pbavg, &state, dir, ptrgeomco1, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
        ctop1 = max(fabs(cmax1), fabs(cmin1));
        ctop2 = max(fabs(cmax2), fabs(cmin2));
        //      alpha=0.5*(ctop1+ctop2); // use average?
        //      alpha=max(ctop1,ctop2); // use maximum?
        // seems alpha can be arbitrary since 0 is ATHENA1

        //      alpha = dx1/dt ; /* crude approx */


        // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
        // if fabs were around each sutracted term, then would be ok (e.g. fabs(B1d-B1u))

        // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
        diffusiveterm[1] =  0.125*(
                                   +ctop1*(
                                           + B2d - B2mm - B2u + B2mp
                                           + B2d - B2pm - B2u + B2pp
                                           )
                                   +ctop2*(
                                           + B3r - B3pm - B3l + B3mm
                                           + B3r - B3pp - B3l + B3mp
                                           )
                                   ) ;
#endif

#if((N1>1)||(N3>1))
        // for emf_2

        // B3 located at FACE3 @ (i-1)
        B3d = 0.5*(
                   MACP0A1(pb,im1mac(i),j,km1mac(k),B3) + 0.5*MACP0A1(dq3,im1mac(i),j,km1mac(k),B3) +
                   MACP0A1(pb,im1mac(i),j,k  ,B3) - 0.5*MACP0A1(dq3,im1mac(i),j,k  ,B3)
                   ) ;

        // B3 located at FACE3 @ (i)
        B3u = 0.5*(
                   MACP0A1(pb,i,j,km1mac(k),B3) + 0.5*MACP0A1(dq3,i,j,km1mac(k),B3) +
                   MACP0A1(pb,i,j,k  ,B3) - 0.5*MACP0A1(dq3,i,j,k  ,B3)
                   ) ;

        // B3 located at FACE1 @ (k-1)
        B1l = 0.5*(
                   MACP0A1(pb,im1mac(i),j,km1mac(k),B1) + 0.5*MACP0A1(dq1,im1mac(i),j,km1mac(k),B1) +
                   MACP0A1(pb,i  ,j,km1mac(k),B1) - 0.5*MACP0A1(dq1,i  ,j,km1mac(k),B1)
                   ) ;
        // B1 located at FACE1 @ k
        B1r = 0.5*(
                   MACP0A1(pb,im1mac(i),j,k,B1) + 0.5*MACP0A1(dq1,im1mac(i),j,k,B1) +
                   MACP0A1(pb,i  ,j,k,B1) - 0.5*MACP0A1(dq1,i  ,j,k,B1)
                   ) ;

        // B3 for all centers around CORN1
        // B3[k - mp][i - mp]
        B3mm = MACP0A1(pb,im1mac(i),j,km1mac(k),B3) ;
        B3mp = MACP0A1(pb,i  ,j,km1mac(k),B3) ;
        B3pm = MACP0A1(pb,im1mac(i),j,k  ,B3) ;
        B3pp = MACP0A1(pb,i  ,j,k  ,B3) ;

        // B1 for all centers around CORN1
        // B1[k - mp][i - mp]
        B1mm = MACP0A1(pb,im1mac(i),j,km1mac(k),B1) ;
        B1mp = MACP0A1(pb,i  ,j,km1mac(k),B1) ;
        B1pm = MACP0A1(pb,im1mac(i),j,k  ,B1) ;
        B1pp = MACP0A1(pb,i  ,j,k  ,B1) ;

        // compute characteristic velocity -- only for Athena2 method


        // average pb to CORN1 for average phase speed there
        PLOOP(pliter,pl) pbavg[pl]=0.25*(MACP0A1(pb,i,j,k,pl)+MACP0A1(pb,im1mac(i),j,k,pl)+MACP0A1(pb,i,j,km1mac(k),pl)+MACP0A1(pb,im1mac(i),j,km1mac(k),pl));
        get_geometry(i, j, k, CORN2, ptrgeomco2); // used here and below emf's
        MYFUN(get_state(pbavg, ptrgeomco2, &state),"step_ch.c:flux_ct()", "get_state()", 1);
        dir=3; MYFUN(vchar(pbavg, &state, dir, ptrgeomco2, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
        dir=1; MYFUN(vchar(pbavg, &state, dir, ptrgeomco2, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
        ctop1 = max(fabs(cmax1), fabs(cmin1));
        ctop2 = max(fabs(cmax2), fabs(cmin2));
        //      alpha=0.5*(ctop1+ctop2); // use average?
        //     alpha=max(ctop1,ctop2); // use maximum?
        // seems alpha can be arbitrary since 0 is ATHENA1

        //      alpha = dx1/dt ; /* crude approx */


        // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
        // if fabs were around each sutracted term, then would be ok (e.g. fabs(B3d-B3u))

        // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
        diffusiveterm[2] =  0.125*(
                                   +ctop1*(
                                           + B3d - B3mm - B3u + B3mp
                                           + B3d - B3pm - B3u + B3pp
                                           )
                                   +ctop2*(
                                           + B1r - B1pm - B1l + B1mm
                                           + B1r - B1pp - B1l + B1mp
                                           )
                                   ) ;
#endif

#if((N1>1)||(N2>1))
        // for emf_3

        // B1 located at FACE1 @ (j-1)
        B1d = 0.5*(
                   MACP0A1(pb,im1mac(i),jm1mac(j),k,B1) + 0.5*MACP0A1(dq1,im1mac(i),jm1mac(j),k,B1) +
                   MACP0A1(pb,i  ,jm1mac(j),k,B1) - 0.5*MACP0A1(dq1,i  ,jm1mac(j),k,B1)
                   ) ;

        // B1 located at FACE1 @ j
        B1u = 0.5*(
                   MACP0A1(pb,im1mac(i),j,k,B1) + 0.5*MACP0A1(dq1,im1mac(i),j,k,B1) +
                   MACP0A1(pb,i  ,j,k,B1) - 0.5*MACP0A1(dq1,i  ,j,k,B1)
                   ) ;

        // B2 located at FACE2 @ (i-1)
        B2l = 0.5*(
                   MACP0A1(pb,im1mac(i),jm1mac(j),k,B2) + 0.5*MACP0A1(dq2,im1mac(i),jm1mac(j),k,B2) +
                   MACP0A1(pb,im1mac(i),j  ,k,B2) - 0.5*MACP0A1(dq2,im1mac(i),j  ,k,B2)
                   ) ;
        // B2 located at FACE2 @ i
        B2r = 0.5*(
                   MACP0A1(pb,i,jm1mac(j),k,B2) + 0.5*MACP0A1(dq2,i,jm1mac(j),k,B2) +
                   MACP0A1(pb,i,j  ,k,B2) - 0.5*MACP0A1(dq2,i,j  ,k,B2)
                   ) ;

        // B1 for all centers around CORN3
        // B1[i - mp][j - mp]
        B1mm = MACP0A1(pb,im1mac(i),jm1mac(j),k,B1) ;
        B1mp = MACP0A1(pb,im1mac(i),j  ,k,B1) ;
        B1pm = MACP0A1(pb,i  ,jm1mac(j),k,B1) ;
        B1pp = MACP0A1(pb,i  ,j  ,k,B1) ;

        // B2 for all centers around CORN3
        // B2[i - mp][j - mp]
        B2mm = MACP0A1(pb,im1mac(i),jm1mac(j),k,B2) ;
        B2mp = MACP0A1(pb,im1mac(i),j  ,k,B2) ;
        B2pm = MACP0A1(pb,i  ,jm1mac(j),k,B2) ;
        B2pp = MACP0A1(pb,i  ,j  ,k,B2) ;

        // compute characteristic velocity -- only for Athena2 method


        // average pb to CORN3 for average phase speed there
        PLOOP(pliter,pl) pbavg[pl]=0.25*(MACP0A1(pb,i,j,k,pl)+MACP0A1(pb,i,jm1mac(j),k,pl)+MACP0A1(pb,im1mac(i),j,k,pl)+MACP0A1(pb,im1mac(i),jm1mac(j),k,pl));
        get_geometry(i, j, k, CORN3, ptrgeomco3); // used here and below emf's
        MYFUN(get_state(pbavg, ptrgeomco3, &state),"step_ch.c:flux_ct()", "get_state()", 1);
        dir=1; MYFUN(vchar(pbavg, &state, dir, ptrgeomco3, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
        dir=2; MYFUN(vchar(pbavg, &state, dir, ptrgeomco3, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
        ctop1 = max(fabs(cmax1), fabs(cmin1));
        ctop2 = max(fabs(cmax2), fabs(cmin2));
        //      alpha=0.5*(ctop1+ctop2); // use average?
        //      alpha=max(ctop1,ctop2); // use maximum?
        // seems alpha can be arbitrary since 0 is ATHENA1

        //      alpha = dx1/dt ; /* crude approx */


        // GODMARK: seems to have left/right and up/down asymmetry due to subtraction
        // if fabs were around each sutracted term, then would be ok (e.g. fabs(B1d-B1u))

        // notice that ctop1 and ctop2 have different "units", so cannot use with B2/B3 arbitrarily, must be consistent.
        diffusiveterm[3] =  0.125*(
                                   +ctop1*(
                                           + B1d - B1mm - B1u + B1mp
                                           + B1d - B1pm - B1u + B1pp
                                           )
                                   +ctop2*(
                                           + B2r - B2pm - B2l + B2mm
                                           + B2r - B2pp - B2l + B2mp
                                           )
                                   ) ;
#endif




        /////////////////////////
        //
        // add geometry
        //
        ////////////////////////

        // must add geometry (no choice since above diffusive correction never had geometry)
        get_geometry_gdetmix(i,j,k,CORN1,ptrgeomf1);
        get_geometry_gdetmix(i,j,k,CORN2,ptrgeomf2);
        get_geometry_gdetmix(i,j,k,CORN3,ptrgeomf3);
      
        // obviously ptrgeom->EOMFUNCMAC(B2) has to be equal to ptrgeom->EOMFUNCMAC(B3) for this method
        diffusiveterm[1] *=(ptrgeomf1->EOMFUNCMAC(B2));
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        diffusiveterm[2] *=(ptrgeomf2->EOMFUNCMAC(B1));
        // obviously ptrgeom->EOMFUNCMAC(B1) has to be equal to ptrgeom->EOMFUNCMAC(B2) for this method
        diffusiveterm[3] *=(ptrgeomf3->EOMFUNCMAC(B1));

        //////////////
        //
        // now add diffusive term to emf
        //
        //////////////
        for(l=1;l<=3;l++){
          MACP1A0(emf,l,i,j,k)+= - diffusiveterm[l];
        }


      }// end EMF loop
    }// end parallel region

  }// end if athena2



  return(0);
}

















/// TOTH: EMF->FLUX
/// OPENMPMARK: Apparently the below loops are expensive due to OpenMP overhead.  Using static schedule helps overhead a bit
int flux_ct_emf2flux(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  // full-type geometry below
  FTYPE coefemf[NDIM];
  extern int choose_limiter(int dir, int i, int j, int k, int pl);





  ////////////////////////////////////
  ///////////////////////////////////
  //
  // compute flux from EMF (signs must be consistent so that similar signed quantity in end)
  // e.g. flux terms originally involving B^2 v^1 - B^1 v^2 must come out with same sign at end
  //
  // This calculation takes care of case when EMF is not needed, so above need not worry about resetting values of emf to 0 or something
  //
  // we set flux to 0 when that dimension doesn't matter, but note that flux really is not 0.  Rather the differences in the flux are 0 across that direction.  However, this difference never enters since this flux is only used as flux differences.
  // These flux differences are *assumed* to vanish if that dimension has N=1.  However, a problem could be setup such that those differences would NOT be 0.
  // For example, a wedge out of an axisymmetric slice that has boundaries that are not symmetric around the equator.  This would have flux differences that would lead to changes in the quantity.  This would be a quasi-2D problem modelled in 1D.
  //
  //////////////////////////////////
  //////////////////////////////////

  // Presume emf[N+SHIFT] so can combine loops, which just ignores values at emf[N+SHIFT-1].  Increases cache miss for F1,F2,F3, but less for emf(1,2,3) and less loop overhead, especially for OpenMP

  if((FLUXB==FLUXCTTOTH)||(FLUXB==ATHENA2)||(FLUXB==ATHENA1)){
    /* rewrite EMFs as fluxes, after Toth */


#pragma omp parallel  // no copyin() needed since not calling any functions that use global vars
    {
      int i,j,k;
      OPENMP3DLOOPVARSDEFINE;

      ////    COMPLOOPOUTFM1dir1full{ // constrain or control better? GODMARK
      //////      OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,INFULL2,OUTFULLM12,INFULL3,OUTFULLM13);
      ////    COMPLOOPOUTFM1dir2full{ // constrain or control better? GODMARK
      //      OPENMP3DLOOPSETUP(INFULL1,OUTFULLM11,-N2BND,N2-1+N2BND,INFULL3,OUTFULLM13);
      //    COMPF2CTZLOOP {
      ////    COMPLOOPOUTFM1dir3full{ // constrain or control better? GODMARK
      //OPENMP3DLOOPSETUP(INFULL1,OUTFULLM11,INFULL2,OUTFULLM12,-N3BND,N3-1+N3BND);
      //    COMPF3CTZLOOP {

      
      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) // not needed with just 1 loop: nowait // Can "nowait" since each F1,F2,F3 set independently on successive loops
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        /////////////////////////////////////
        // F1
        ////////////////////////////////////
#if(N1>1)
        // below line always true
        MACP0A1(F1,i,j,k,B1) = 0.;
        MACP0A1(F1,i,j,k,B2) = 0.5 * (MACP1A0(emf,3,i,j,k) + MACP1A0(emf,3,i,jp1mac(j),k)); // put emf3 back to FACE1
        MACP0A1(F1,i,j,k,B3) = - 0.5 * (MACP1A0(emf,2,i,j,k) + MACP1A0(emf,2,i,j,kp1mac(k))); // put emf2 back to FACE1
#endif



        /////////////////////////////////////
        // F2
        ////////////////////////////////////
#if(N2>1)
        MACP0A1(F2,i,j,k,B1) = - 0.5 * (MACP1A0(emf,3,i,j,k) + MACP1A0(emf,3,ip1mac(i),j,k));
        // below line always true
        MACP0A1(F2,i,j,k,B2) = 0.;
        MACP0A1(F2,i,j,k,B3) = 0.5 * (MACP1A0(emf,1,i,j,k) + MACP1A0(emf,1,i,j,kp1mac(k)));
#endif

        /////////////////////////////////////
        // F3
        ////////////////////////////////////
#if(N3>1)
        MACP0A1(F3,i,j,k,B1) = 0.5 * (MACP1A0(emf,2,i,j,k) + MACP1A0(emf,2,ip1mac(i),j,k));
        MACP0A1(F3,i,j,k,B2) = - 0.5 * (MACP1A0(emf,1,i,j,k) + MACP1A0(emf,1,i,jp1mac(j),k));
        // below line always true
        MACP0A1(F3,i,j,k,B3) = 0.;
#endif
      }// end loop block [using outer loop since even though 3X more expensive cache-wise for F1, less so for emf.  Reduces OpenMP overhead alot.
    }// end parallel region (and implied barrier)





  } // end if FLUXCT
  else if(FLUXB==FLUXCD){





    // F's are emf's, where dF/dx is (F(i+1)-F(i-1))/dx
    // fluxes remain at center!

#pragma omp parallel  // no copyin() required
    {
      int i,j,k;
      OPENMP3DLOOPVARSDEFINE; 



      //    COMPF1CTZLOOP {
      ////    COMPFULLLOOP{ // SHOULD try to be more constrained GODMARK


      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) // nowait // Can "nowait" since each successive loop sets F1,F2,F3 independently
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        /////////////////////////////////////
        // F1
        ////////////////////////////////////
#if(N1>1)
        // always below line
        MACP0A1(F1,i,j,k,B1) = 0.0;
        MACP0A1(F1,i,j,k,B2) = MACP1A0(emf,3,i,j,k);
        MACP0A1(F1,i,j,k,B3) = -MACP1A0(emf,2,i,j,k);
#endif // end if doing F1



        /////////////////////////////////////
        // F2
        ////////////////////////////////////
#if(N2>1)
        MACP0A1(F2,i,j,k,B1) = -MACP1A0(emf,3,i,j,k);
        // always below line
        MACP0A1(F2,i,j,k,B2) = 0.;
        MACP0A1(F2,i,j,k,B3) = MACP1A0(emf,1,i,j,k);
#endif // end if doing F2



        /////////////////////////////////////
        // F3
        ////////////////////////////////////
#if(N3>1)
        MACP0A1(F3,i,j,k,B1) = MACP1A0(emf,2,i,j,k);
        MACP0A1(F3,i,j,k,B2) = -MACP1A0(emf,1,i,j,k);
        // always below line
        MACP0A1(F3,i,j,k,B3) = 0.;
#endif // end if doing F3
      }// end loop block [using outer loop since even though 3X more expensive cache-wise for F1, less so for emf.  Reduces OpenMP overhead alot.
    }// end parallel region (and implied barrier)



  } // end if FLUXCD



  return(0);
}
