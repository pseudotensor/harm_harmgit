#include "decs.h"

/*! \file copyandinit_functions.c&
    \brief copy/init arrays
    //
    // Various mostly general copy and initialization functions that operate on 3D Loop
    //
    /////////////////////////

*/






/// loop range for inversion or final centered field primitive
void get_inversion_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke)
{


#if(0)
  // define loop range
  // +SHIFT? is for IF3DSPCTHENMPITRANSFERATPOLE
  // with cleanup_fluxes() in flux.c so that fluxes are zeroed-out outside well-defined computational box, then *always* update +-1 from "normal" conservatives so that all possible fluxes are accounted in changes due to interior and some exterior cells
  // This way all advance.c is simple and just expanded by 1 cell effectively
  // This allows use of adaptive time-stepping such that surrounding cells are properly updated and keep pace with the effective time of the RK-stepping
  // Also ensures divb=0 and flux conservation under any case since always take into account fluxes through the well-defined computational box
  // must constrain result since this is final centered value of field and only have enough information to be well-defined on the current computational box over a finite range
  // Only expand if on outer edge of not-evolved region in order to (primarily) preserve divb=0
  
  // must constrain result since this is final centered value of field and only have enough information to be well-defined on the current computational box over a finite range
  // Only expand if on outer edge of not-evolved region in order to (primarily) preserve divb=0

  // 1|| because realized don't want to adjust "boundary cells".  Want to keep them fixed.  Violates conservation (and divb=0) unless separately evolve that other region.  This is ok since just approximating evolution in non-evolved region when moving full boundary.
  // if(subgrid inner boundary>global active grid inner boundary)
  if(AVOIDADVANCESHIFTX1DN||enerposreg[ACTIVEREGION][X1DN]>enerposreg[ACTIVEREGION][X1DN]) *is=Uconsevolveloop[FIS];
  else *is=Uconsevolveloop[FIS]-SHIFT1;

  if(AVOIDADVANCESHIFTX1UP||enerposreg[ACTIVEREGION][X1UP]<enerposreg[ACTIVEREGION][X1UP]) *ie=Uconsevolveloop[FIE];
  else *ie=Uconsevolveloop[FIE]+SHIFT1;

  if(AVOIDADVANCESHIFTX2DN||enerposreg[ACTIVEREGION][X2DN]>enerposreg[ACTIVEREGION][X2DN]) *js=Uconsevolveloop[FJS];
  else *js=Uconsevolveloop[FJS]-SHIFT2;

  if(AVOIDADVANCESHIFTX2UP||enerposreg[ACTIVEREGION][X2UP]<enerposreg[ACTIVEREGION][X2UP]) *je=Uconsevolveloop[FJE];
  else *je=Uconsevolveloop[FJE]+SHIFT2;

  if(AVOIDADVANCESHIFTX3DN||enerposreg[ACTIVEREGION][X3DN]>enerposreg[ACTIVEREGION][X3DN]) *ks=Uconsevolveloop[FKS];
  else *ks=Uconsevolveloop[FKS]-SHIFT3;

  if(AVOIDADVANCESHIFTX3UP||enerposreg[ACTIVEREGION][X3UP]<enerposreg[ACTIVEREGION][X3UP]) *ke=Uconsevolveloop[FKE];
  else *ke=Uconsevolveloop[FKE]+SHIFT3;
#else

  // this loop range must be equal to that used in copy_tempucum_finalucum() for centered quantities there (i.e. ignore FLUXB==FLUXCTSTAG in that function as compared to here).
  *is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
  *ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
  *js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
  *je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
  *ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
  *ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);

  // new method is always centered inversion or primitive for field
  //*is=loop[FIS];
  //  *ie=loop[FIE];
  //  *js=loop[FJS];
  //  *je=loop[FJE];
  //  *ks=loop[FKS];
  //  *ke=loop[FKE];

#endif

}

/// determine loop range for ustagpoint2pstag() in fluxctstag.c
void get_stag_startendindices(int *loop, int dir, int *is,int *ie,int *js,int *je,int *ks,int *ke)
{
  
  // must constrain result since this is final centered value of field and only have enough information to be well-defined on the current computational box over a finite range
  // Only expand if on outer edge of not-evolved region in order to (primarily) preserve divb=0

  // 1|| for same reason as for centered quantities as discussed above.
  // if(subgrid inner boundary>global active grid inner boundary)
  if(AVOIDADVANCESHIFTX1DN||enerposreg[ACTIVEREGION][X1DN]>enerposreg[ACTIVEREGION][X1DN]) *is=loop[FIS];
  else *is=loop[FIS]-SHIFT1;

  if(AVOIDADVANCESHIFTX1UP||enerposreg[ACTIVEREGION][X1UP]<enerposreg[ACTIVEREGION][X1UP]) *ie=loop[FIE]+SHIFT1*(dir==1);
  else *ie=loop[FIE]+SHIFT1;

  if(AVOIDADVANCESHIFTX2DN||enerposreg[ACTIVEREGION][X2DN]>enerposreg[ACTIVEREGION][X2DN]) *js=loop[FJS];
  else *js=loop[FJS]-SHIFT2;

  if(AVOIDADVANCESHIFTX2UP||enerposreg[ACTIVEREGION][X2UP]<enerposreg[ACTIVEREGION][X2UP]) *je=loop[FJE]+SHIFT2*(dir==2);
  else *je=loop[FJE]+SHIFT2;

  if(AVOIDADVANCESHIFTX3DN||enerposreg[ACTIVEREGION][X3DN]>enerposreg[ACTIVEREGION][X3DN]) *ks=loop[FKS];
  else *ks=loop[FKS]-SHIFT3;

  if(AVOIDADVANCESHIFTX3UP||enerposreg[ACTIVEREGION][X3UP]<enerposreg[ACTIVEREGION][X3UP]) *ke=loop[FKE]+SHIFT3*(dir==3);
  else *ke=loop[FKE]+SHIFT3;

}



/// determine loop range for *using* fluxes
/// shifts extra +1 to get field update.
/// Assumes applied to temporary ucum and that final ucum more controlled
/// This forces flux update so face fields included, but final ucum/pf still only updated at center.
/// +SHIFT1/2/3 required in particular for IF3DSPCTHENMPITRANSFERATPOLE.  In all other cases turns out wouldn't have needed this, but still ok to do it in general -- especially for AMR
/// Notes from above:
/// then always use +-1 expanded loop even for inversion
/// This overdoes corner cells for inversion for centered U, but just cycles through existing values so ok
/// Even ok if extends to true boundary cell
/// with cleanup_fluxes() in flux.c so that fluxes are zeroed-out outside well-defined computational box, then *always* update +-1 from "normal" conservatives so that all possible fluxes are accounted in changes due to interior and some exterior cells
/// This way all advance.c is simple and just expanded by 1 cell effectively
/// This allows use of adaptive time-stepping such that surrounding cells are properly updated and keep pace with the effective time of the RK-stepping
/// Also ensures divb=0 and flux conservation under any case since always take into account fluxes through the well-defined computational box
/// must constrain result since this is final centered value of field and only have enough information to be well-defined on the current computational box over a finite range
/// Only expand if on outer edge of not-evolved region in order to (primarily) preserve divb=0
void get_flux_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke)
{

  // this loop range must be equal or larger than that used in copy_tempucum_finalucum()
  *is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
  *ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
  *js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
  *je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
  *ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
  *ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);


  if(FLUXB==FLUXCTSTAG){
    // generic shift upwards to compute necessary things for staggered field.  In end presume however shifted here, ucum and final primitives are computed on highly controlled locations only

    // overrides:
    *ie=loop[FIE]+SHIFT1;
    *je=loop[FJE]+SHIFT2;
    *ke=loop[FKE]+SHIFT3;
  }
  else{
    // default is good
  }
}








/// copy tempucum -> ucum so ucum only has updates where wanted for each pl.
/// avoids NaN or out of bounds assignments into ucum (which shows up in, e.g., dump diagnostics)
/// Also ensures that final unewglobal values are only updated within the well-defined computational box
/// Note that inversion U->p is also within well-defined computational box.
/// Note that ustag->pstag occurs in fluxctstag.c in well-defined computational box [with extra face for each B1,B2,B3 as required]
/// So overall prim,pstag,unew are only updated where desired -- no leakage unlike fluxes, temp primitives, and other things.
void copy_tempucum_finalucum(int whichpl, int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel  // just copying, only need PLOOP even for ucum_check()
  {
    int i,j,k;
    int pl,pliter;
    int is,ie,js,je,ks,ke;
    extern void ucum_check(int i, int j, int k, int loc, int pl, FTYPE *ucum);

    // loop range where final ucum is set from scratch ucum that may have been set outside desired region for simplicity of loop structures
    is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
    ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
    js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
    je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
    ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
    ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);


    if(FLUXB==FLUXCTSTAG){
      if(whichpl==DOALLPL || whichpl==DONONBPL){
        // do non-field quantities
        copy_3d_nofield_nowait(is, ie, js, je, ks, ke, tempucum,ucum);

#if(PRODUCTION==0)
#pragma omp barrier // force barrier since otherwise nowait will leak into here with undefined values in general
        COMPZSLOOP(is,ie,js,je,ks,ke){
          PLOOPNOB1(pl) ucum_check(i,j,k,CENT,pl, MAC(ucum,i,j,k));
          PLOOPNOB2(pl) ucum_check(i,j,k,CENT,pl, MAC(ucum,i,j,k));
        }
#endif
      }

      // do field quantities
      if(whichpl==DOALLPL || whichpl==DOBPL){

        // do pl==B1
        pl=B1;
        is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
        ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
        js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
        je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
        ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
        ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);

        ie=loop[FIE]+SHIFT1; // always shift - override
        copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );

#if(PRODUCTION==0)
#pragma omp barrier // force barrier since otherwise nowait will leak into here with undefined values in general
        COMPZSLOOP(is,ie,js,je,ks,ke){
          ucum_check(i,j,k,FACE1,pl, MAC(ucum,i,j,k));
        }
#endif


        // do pl==B2
        pl=B2;
        is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
        ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
        js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
        je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
        ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
        ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);

        je=loop[FJE]+SHIFT2;
        copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );


#if(PRODUCTION==0)
#pragma omp barrier // force barrier since otherwise nowait will leak into here with undefined values in general
        COMPZSLOOP(is,ie,js,je,ks,ke){
          ucum_check(i,j,k,FACE2,pl, MAC(ucum,i,j,k));
        }
#endif


        // do pl==B3
        pl=B3;
        is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
        ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
        js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
        je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
        ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
        ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);

        ke=loop[FKE]+SHIFT3;
        copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );


#if(PRODUCTION==0)
#pragma omp barrier // force barrier since otherwise nowait will leak into here with undefined values in general
        COMPZSLOOP(is,ie,js,je,ks,ke){
          ucum_check(i,j,k,FACE3,pl, MAC(ucum,i,j,k));
        }
#endif
      }// end if need to to B type pl

      // now ucum is assigned only where should be changed

    }
    else{
      // nothing to do since tempucum is ucum and all at CENT
      // just check
#if(PRODUCTION==0)
#pragma omp barrier // force barrier since otherwise nowait will leak into here with undefined values in general
      COMPZSLOOP(is,ie,js,je,ks,ke){
        PLOOP(pliter,pl){
          if(whichpl==DONONBPL && BPL(pl)==1 || whichpl==DOBPL && BPL(pl)==0) continue;
          ucum_check(i,j,k,CENT,pl, MAC(ucum,i,j,k));
        }
      }
#endif
    }

  }// end parallel region (with implicit barrier)

}


/// like copy_tempucum_finalucum() but for field only
/// Used for setting up point value of field in advance.c
void copy_tempucum_finalucum_fieldonly(int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel  // just copying, only need PLOOP even for ucum_check()
  {
    int i,j,k;
    int pl,pliter;
    int is,ie,js,je,ks,ke;
    extern void ucum_check(int i, int j, int k, int loc, int pl, FTYPE *ucum);

    // loop range where final ucum is set from scratch ucum that may have been set outside desired region for simplicity of loop structures
    is=loop[FIS]-SHIFT1*(AVOIDADVANCESHIFTX1DN==0);
    ie=loop[FIE]+SHIFT1*(AVOIDADVANCESHIFTX1UP==0);
    js=loop[FJS]-SHIFT2*(AVOIDADVANCESHIFTX2DN==0);
    je=loop[FJE]+SHIFT2*(AVOIDADVANCESHIFTX2UP==0);
    ks=loop[FKS]-SHIFT3*(AVOIDADVANCESHIFTX3DN==0);
    ke=loop[FKE]+SHIFT3*(AVOIDADVANCESHIFTX3UP==0);


    if(FLUXB==FLUXCTSTAG){

      // do pl==B1
      pl=B1;
      ie=loop[FIE]+SHIFT1; // always shift - override
      copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );

#if(PRODUCTION==0)
      COMPZSLOOP(is,ie,js,je,ks,ke){
        ucum_check(i,j,k,FACE1,pl, MAC(ucum,i,j,k));
      }
#endif


      // do pl==B2
      pl=B2;
      je=loop[FJE]+SHIFT2;
      copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );


#if(PRODUCTION==0)
      COMPZSLOOP(is,ie,js,je,ks,ke){
        ucum_check(i,j,k,FACE2,pl, MAC(ucum,i,j,k));
      }
#endif


      // do pl==B3
      pl=B3;
      ke=loop[FKE]+SHIFT3;
      copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, tempucum, ucum );


#if(PRODUCTION==0)
      COMPZSLOOP(is,ie,js,je,ks,ke){
        ucum_check(i,j,k,FACE3,pl, MAC(ucum,i,j,k));
      }
#endif


      // now ucum is assigned only where should be changed

    }
    else{
      // nothing to do since tempucum is ucum and all at CENT
      // just check
#if(PRODUCTION==0)
      COMPZSLOOP(is,ie,js,je,ks,ke){
        PLOOPBONLY(pl) ucum_check(i,j,k,CENT,pl, MAC(ucum,i,j,k));
      }
#endif
    }

  }// end parallel region (with implicit barrier)

}








/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      ////      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);
      }
    }// end 3D loop


  }// end parallel region

}

// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
// put as function because then wrap-up OpenMP stuff
void copy_3dnpr_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  

  copy_3dnpr(is, ie, js, je, ks, ke,source, dest);

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3d_nofield(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOPNOB1(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);
      PLOOPNOB2(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

    }// end 3D loop


  }// end parallel region

}


/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOPBONLY(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

    }// end 3D loop


  }// end parallel region

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3d_fieldonly_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  

  copy_3d_fieldonly(is, ie, js, je, ks, ke, source, dest);


}


/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
/// Presumes parallel region is outside function
void copy_3d_nofield_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;

  // already inside parallel region
  OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) nowait
  OPENMP3DLOOPBLOCK{
    OPENMP3DLOOPBLOCK2IJK(i,j,k);

    //      COMPZSLOOP(is,ie,js,je,ks,ke){
    PLOOPNOB1(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);
    PLOOPNOB2(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

  }// end 3D loop


}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
/// Presumes parallel region is outside function
void copy_3d_fieldonly_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;

  // already inside parallel region
  OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE()) nowait
  OPENMP3DLOOPBLOCK{
    OPENMP3DLOOPBLOCK2IJK(i,j,k);

    //      COMPZSLOOP(is,ie,js,je,ks,ke){
    PLOOPBONLY(pl) MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

  }// end 3D loop


}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3d_onepl(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel
  {
    int i,j,k;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

    }// end 3D loop


  }// end parallel region

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
/// Presumes parallel region is outside function
void copy_3d_onepl_nowait(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;


  // already inside parallel region
  OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
  OPENMP3DLOOPBLOCK{
    OPENMP3DLOOPBLOCK2IJK(i,j,k);

    //      COMPZSLOOP(is,ie,js,je,ks,ke){
    MACP0A1(dest,i,j,k,pl)=MACP0A1(source,i,j,k,pl);

  }// end 3D loop

}


/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3d_onepl_fullloop(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  

  copy_3d_onepl(is, ie, js, je, ks, ke, pl, source, dest);

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
/// Presumes parallel region is outside function
void copy_3d_onepl_fullloop_nowait(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  

  copy_3d_onepl_nowait(is, ie, js, je, ks, ke, pl, source, dest);

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void init_3dnpr_flux(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);


#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest,i,j,k,pl)=initvalue;
      }
    }// end 3D loop


  }// end parallel region

}

/// init 3d npr type array
void init_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{


#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);


#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest,i,j,k,pl)=initvalue;
      }
    }// end 3D loop


  }// end parallel region

}

/// init 3d npr type array over full loop
void init_3dnpr_fullloop(FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  
  init_3dnpr(is,ie,js,je,ks,ke,initvalue,dest);

}

/// init 3d npr type array over full loop over flux positions
void init_3dnpr_fullloop_flux(FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  
  init_3dnpr_flux(is,ie,js,je,ks,ke,initvalue,dest);

}


/// initialize single pre-component of the vpot type array
void init_3dvpot(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{


#pragma omp parallel
  {
    int i,j,k;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
      MAC(dest,i,j,k)=initvalue;
    }// end 3D loop
  }// end parallel region

}


/// initialize single pre-component of the vpot type array
void init_3dvpot_fullloopp1(FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int is=-N1BND;
  int ie=N1-1+N1BND+SHIFT1;
  int js=-N2BND;
  int je=N2-1+N2BND+SHIFT2;
  int ks=-N3BND;
  int ke=N3-1+N3BND+SHIFT3;
  
  init_3dvpot(is,ie,js,je,ks,ke,initvalue,dest);

}



/// initialize single pre-component of the vpot type array
void copy_3dvpot(int is, int ie, int js, int je, int ks, int ke, FTYPE (*source)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{


#pragma omp parallel
  {
    int i,j,k;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
      MAC(dest,i,j,k)=MAC(source,i,j,k);
    }// end 3D loop
  }// end parallel region

}


/// initialize single pre-component of the vpot type array
void copy_3dvpot_fullloopp1(FTYPE (*source)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int is=-N1BND;
  int ie=N1-1+N1BND+SHIFT1;
  int js=-N2BND;
  int je=N2-1+N2BND+SHIFT2;
  int ks=-N3BND;
  int ke=N3-1+N3BND+SHIFT3;
  
  copy_3dvpot(is,ie,js,je,ks,ke,source,dest);

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void init_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest1,i,j,k,pl)=MACP0A1(dest2,i,j,k,pl)=initvalue;
      }
    }// end 3D loop


  }// end parallel region

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void init_3dnpr_3ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR],FTYPE (*dest3)[NSTORE2][NSTORE3][NPR])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest1,i,j,k,pl)=MACP0A1(dest2,i,j,k,pl)=MACP0A1(dest3,i,j,k,pl)=initvalue;
      }
    }// end 3D loop


  }// end parallel region

}




/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR])
{

#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pliter,pl){
        MACP0A1(dest1,i,j,k,pl)=MACP0A1(dest2,i,j,k,pl)=MACP0A1(source,i,j,k,pl);
      }
    }// end 3D loop


  }// end parallel region

}







/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dpftype_special(int is, int ie, int js, int je, int ks, int ke,PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3][NUMFAILPFLAGS])
{


#pragma omp parallel
  {
    int pf;
    int i,j,k;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      FAILPFLAGLOOP(pf) MACP0A1(destspecial,i,j,k,pf)=MACP0A1(source,i,j,k,pf);

    }// end 3D loop

  }// end parallel region

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dpftype_special_fullloop(PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3][NUMFAILPFLAGS])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  
  // override
  get_inversion_startendindices(Uconsevolveloop,&is,&ie,&js,&je,&ks,&ke);
  
  // +-NUMPFLAGBNDx extra so can do check
  is += -NUMPFLAGBND1;
  ie += +NUMPFLAGBND1;
  js += -NUMPFLAGBND2;
  je += +NUMPFLAGBND2;
  ks += -NUMPFLAGBND3;
  ke += +NUMPFLAGBND3;


  copy_3dpftype_special(is, ie, js, je, ks, ke,source, destspecial);

}



/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dnpr2interp_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP])
{


#pragma omp parallel OPENMPGLOBALPRIVATEPLOOPINTERPONLY
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PINTERPLOOP(pliter,pl){
        MACP0A1(dest1,i,j,k,pl)=MACP0A1(dest2,i,j,k,pl)=MACP0A1(source,i,j,k,pl);
      }
    }// end 3D loop


  }// end parallel region

}

/// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
/// put as function because then wrap-up OpenMP stuff
void copy_3dnpr2interp_2ptrs_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP])
{
  int i,j,k,pl,pliter;
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;

  copy_3dnpr2interp_2ptrs(is, ie, js, je, ks, ke,source, dest1, dest2);

}

