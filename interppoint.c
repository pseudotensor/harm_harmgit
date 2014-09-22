
/*! \file interppoint.c
     \brief Spatial Interpolation for fluxes based upon providing each point
     // Inefficient compared to interpline.c, but simpler.
*/

#include "decs.h"





/// whether to send all pl's for access by interpolator (such as for steepening and flattening)
/// whether to use contact steepener for reallim=MCSTEEP
#define DOUSEPPMCONTACTSTEEP 1 // used alone is unstable
#define DOUSEPARAMONOCHECK 1 // above must be used with this
/// using the 2 above, contacts are as steep as with PARA!

#define DOUSEPARAFLAT 1 // avoids oscillations in strong shocks

#define DOINGALLMCSTEEP (DOUSEPPMCONTACTSTEEP&&DOUSEPARAFLAT&&DOUSEPARAMONOCHECK)

#define WHICH3POINTLIMT MC

/// below not accounted for unless DOINGALLMCSTEEP==1
/// for blast wave (sasha TESTNUMBER==8) actually helps contact but hurts only peak a bit -- contact much better relatively speaking)
#define DQALLOWEXTREMUM 1 // assume USEPARAMONOCHECK enabled if using this!!! (otherwise crazy)


/// whether to interpolate and enforce steepening on u and rho in consistent way
/// see comments where used first as to why turned off
#define CONNECTUANDRHO (0 && (DOEVOLVERHO&&DOEVOLVEUU))

/// whether to interpolate and enforce steepening on B and rho in consistent way
/// Should I turn this off?  Maybe for same reason of U version turned off -- may cause current sheet to spontaneously dissipate -- just use separate steepening
/// Seems to cause lots of problems with MCSTEEP and MHD torus problem -- too much steepening
#define CONNECTBANDRHO (0 && (DOEVOLVERHO) )


#if(DQALLOWEXTREMUM==0)
#define MCSTEEPMINMOD(a,b) MINMOD(a,b)
#else
#define MCSTEEPMINMOD(a,b) MINMODB(a,b)
#endif


/// whether primreal is p2interp
/// can't assume this is true in general (e.g. for FIELDCTSTAG input is not even complete range)
//#define YREALISY (RESCALEINTERP==0 || VARTOINTERP==PRIMTOINTERP)




/// interpolation with loop over POINTS
void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int loc, int continuous, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP])
{
  void slope_lim_point_c2e(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right);
  void slope_lim_point_e2c_continuous(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right);
  void slope_lim_point_allpl(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int startorderi, int endorderi, FTYPE **yreal, FTYPE **y, FTYPE *dq,FTYPE *left,FTYPE *right);
  extern int choose_limiter(int dir, int i, int j, int k, int pl);
  FTYPE interplist[MAXSPACEORDER];
  FTYPE realinterplist[MAXSPACEORDER];
  FTYPE interplistpl[NPR2INTERP][MAXSPACEORDER];
  FTYPE realinterplistpl[NPR2INTERP][MAXSPACEORDER];
  FTYPE *yin;
  int ijkshift;
  int startorderi,endorderi;
  int reallim;
  int is,ie,js,je,ks,ke,di,dj,dk;




  // make sure dq's exist if going to use them
  if( (DODQMEMORY==0)&&(LIMADJUST==LIMITERFIXED) ){
    dualfprintf(fail_file,"Must have dq's when using MC or lower second order methods\n");
    myexit(17);
  }

  // set range of positions interpolated to
  set_interppoint_loop_ranges(interporflux, dir, loc, continuous, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);
  
  //  loc=CENT;
  //  continuous=0;

  {
    int i,j,k;
    // GODMARK: for now assume not doing per-point choice of limiter (note originally pl didn't actually have correct choice of limiter since pl is fake and pl set to final pl at end of loop
    i=j=k=0;
    reallim=choose_limiter(dir, i,j,k,pl);
    if(continuous==0){
      // get starting point for stencil, assumed quasi-symmetric (including even/odd size of stencil)
      // for 2nd order: i-1, i, i+1
      startorderi = - interporder[reallim]/2;
      endorderi   = - startorderi;
    }
    else{
      // get starting and ending point for stencil centered around face
      // for 2nd order: i-1, i, i+1, i+2
      // NOTE: Overall loop range is more limited, so extra i+2 not an issue
      startorderi = - interporder[reallim]/2;
      endorderi   = - startorderi + 1;
    }
  }


  ////////////////////
  //
  // For limiters that use specific information about the equations and need all variables, use slope_lim_point_allpl()
  //
  ///////////////////
  if(reallim==PARAFLAT || reallim == MCSTEEP){


    //////////////////////////////
    //
    // Loop over points and quantities both
    //
    //////////////////////////////
#if( DOUSEPARAFLAT || DOUSEPPMCONTACTSTEEP)
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP
#else
    // don't need full copyin() unless pressure needs to be computed as computed if prior conditional holds
#pragma omp parallel
#endif
    {
      int i,j,k,l;
      int plpl,pliter;
      FTYPE *ypl[NPR2INTERP];
      FTYPE *yrealpl[NPR2INTERP];

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

      //////////////////////////////
      //
      // pointer shift
      //
      //////////////////////////////

      PINTERPLOOP(pliter,plpl){
        //      dualfprintf(fail_file,"PLOOPINTERP: plpl=%d\n",plpl);
        ypl[plpl] = interplistpl[plpl] - startorderi;
      }
      PALLREALLOOP(plpl){ // need all variables for real quantity
        // in general need space to be separate in case modify ypl's in some way during interpolation, so can't have pointer yrealpl->ypl
        //      dualfprintf(fail_file,"PALLREALLOOP: plpl=%d\n",plpl);
        yrealpl[plpl] = realinterplistpl[plpl] - startorderi;
      }


#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        PINTERPLOOP(pliter,plpl) for(l=startorderi;l<=endorderi;l++){
          // get interpolation points, where y[0] is point of interest for which interpolation is found.
          ypl[plpl][l]=MACP0A1(p2interp,i + l*idel,j + l*jdel,k + l*kdel,plpl);
        }

        if(realisinterp){
          // need all quantities for real var
          // PLOOPINTERP is used because PALLREALLOOP can be different (extra things at end not part of "real" set
          PINTERPLOOP(pliter,plpl) for(l=startorderi;l<=endorderi;l++){
            yrealpl[plpl][l]=ypl[plpl][l];// faster to copy from ypl than duplicating primreal access
          }
        }
        else{
          // need all quantities for real var
          PALLREALLOOP(plpl) for(l=startorderi;l<=endorderi;l++){
            yrealpl[plpl][l]=MACP0A1(primreal,i + l*idel,j + l*jdel,k + l*kdel,plpl);
          }
        }

        slope_lim_point_allpl(i, j, k, loc, realisinterp, dir,reallim,startorderi,endorderi,yrealpl,ypl,
                              &MACP0A1(dq,i,j,k,0),&MACP0A1(pleft,i,j,k,0),&MACP0A1(pright,i,j,k,0)
                              );
      }// end over loop
    }// end parallel region
  } // end if doing all plpl at once
  else{

    //////////////////////////////
    //
    // Loop over points only
    //
    //////////////////////////////

    // For limiters that are general, use slope_lim_point()
    // get interpolation points, where y[0] is point of interest for which interpolation is found.
#if( DOUSEPARAFLAT || DOUSEPPMCONTACTSTEEP)
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP
#else
    // don't need full copyin() unless pressure needs to be computed as computed if prior conditional holds
#pragma omp parallel
#endif
    {
      int i,j,k,l;
      FTYPE *y;
      FTYPE *yreal;

      OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);

      //////////////////////
      //
      // shift for easy use and clarity of loop and easy of use by slope_lim_point()
      //
      ////////////////////////

      y = interplist - startorderi;
      yreal = realinterplist - startorderi;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        for(l=startorderi;l<=endorderi;l++){
          y[l]=MACP0A1(p2interp,i + l*idel,j + l*jdel,k + l*kdel,pl);
        }
        if(realisinterp){
          for(l=startorderi;l<=endorderi;l++){
            yreal[l]=y[l];
          }
        }
        else{
          // if loc!=CENT, then because primreal always at CENT, assume use of yreal will be treated as being CENT even if inputted quantity to interpolate does not start (or even end at) CENT
          for(l=startorderi;l<=endorderi;l++){
            yreal[l]=MACP0A1(primreal,i + l*idel,j + l*jdel,k + l*kdel,pl);
          }
        }

        //        dualfprintf(fail_file,"godpl=%d\n",pl);

        if(continuous==0){
          slope_lim_point_c2e(i, j, k, loc, realisinterp, dir, reallim,pl,startorderi,endorderi,yreal,y,
                              &MACP0A1(dq,i,j,k,pl),&MACP0A1(pleft,i,j,k,pl),&MACP0A1(pright,i,j,k,pl)
                              );
        }
        else{
          slope_lim_point_e2c_continuous(i, j, k, loc, realisinterp, dir, reallim,pl,startorderi,endorderi,yreal,y,
                              &MACP0A1(dq,i,j,k,pl),&MACP0A1(pleft,i,j,k,pl),&MACP0A1(pright,i,j,k,pl)
                              );
        }

      }// end over loop
    }// end parallel region
  }// end if doing per pl


  if(reallim==PARAFLAT){ // assumes last chosen reallim is constant // SUPERGODMARK
    pl=NPR2INTERP; // finish the loop outside this function // SPLITNPR ok
    if(LIMADJUST!=LIMITERFIXED){
      dualfprintf(fail_file,"Can't use PARAFLAT with LIMADJUST!=LIMITERFIXED\n");
      myexit(1);
    }
  }

}








/// given dir=fluxdir, return ranges for loop over which obtain those interpolated quantities
/// these is,ie,js,je,ks,ke are used inside a loop that has grid section SHIFTS already, so don't add shift here
/// Note that if direction doesn't exist should still return reasonable result that is not out of range (e.g. dir==3 if N3==1 should still give ks=ke=0)
int set_interppoint_loop_ranges(int interporflux, int dir, int loc, int continuous, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{

  if(continuous==0){
    //  if(useghostplusactive){
    set_interppoint_loop_expanded(interporflux, dir, loc, is, ie, js, je, ks, ke, di, dj, dk);
    //  }
    //  else{
    //    set_interppoint_loop(interporflux, dir, loc, is, ie, js, je, ks, ke, di, dj, dk);
    //  }
  }
  else{
    set_interppoint_loop_expanded_face2cent(interporflux, dir, loc, is, ie, js, je, ks, ke, di, dj, dk);
  }


  return(0);
}



int set_interppoint_loop_ranges_3Dextended(int interporflux, int loc, int continuous, int *maxis, int *maxie, int *maxjs, int *maxje, int *maxks, int *maxke, int *di, int *dj, int *dk)
{
  int firsttime;
  int dimen;
  int is,ie,js,je,ks,ke;

  ////////////////////////
  //
  // define loop range
  // loop range should be same as where computed fluxes, or where formed primitive interpolation, as in interppoint.c:
  // since dimension-dependent, and need CENT states in all such cells, expand loop over all dimensions
  //
  ////////////////////////
  firsttime=1;
  DIMENLOOP(dimen){

    if(firsttime){
      firsttime=0;
      // set range of positions interpolated to
      set_interppoint_loop_ranges(interporflux, dimen, loc, continuous, maxis, maxie, maxjs, maxje, maxks, maxke, di, dj, dk);
    }
    else{
      // max here means largest extent in the direction of that boundary as from interior of domain
      set_interppoint_loop_ranges(interporflux, dimen, loc, continuous, &is, &ie, &js, &je, &ks, &ke, di, dj, dk);
      if(is<*maxis) *maxis=is;
      if(ie>*maxie) *maxie=ie;
      if(js<*maxjs) *maxjs=js;
      if(je>*maxje) *maxje=je;
      if(ks<*maxks) *maxks=ks;
      if(ke>*maxke) *maxke=ke;
    }
  }

  return(0);

}


/// range over which need CENT cells defined for EMF in plane defined by corner
/// Note that is not like range needed for interpolation itself, which requires for FLUXCTSTAG all off-dir cells for Flux^{dir}
/// Note that corner ranges from i=0..N so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
void set_interppoint_loop_ranges_2D_EMF_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{
  int *myUconsloop;


  if(interporflux==ENOINTERPTYPE4EMF){
    myUconsloop=emfUconsloop;
  }
  else{
    myUconsloop=Uconsloop;
  }


  if(corner==1){

    *is=0;
    //    *ie=N1-1;
    *ie=N1;

    *js=myUconsloop[FJS]-SHIFT2;
    *je=myUconsloop[FJE]+SHIFT2;

    *ks=myUconsloop[FKS]-SHIFT3;
    *ke=myUconsloop[FKE]+SHIFT3;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(corner==2){

    *is=myUconsloop[FIS]-SHIFT1;
    *ie=myUconsloop[FIE]+SHIFT1;

    *js=0;
    //    *je=N2-1;
    *je=N2;

    *ks=myUconsloop[FKS]-SHIFT3;
    *ke=myUconsloop[FKE]+SHIFT3;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(corner==3){

    *is=myUconsloop[FIS]-SHIFT1;
    *ie=myUconsloop[FIE]+SHIFT1;

    *js=myUconsloop[FJS]-SHIFT2;
    *je=myUconsloop[FJE]+SHIFT2;

    *ks=0;
    //    *ke=N3-1;
    *ke=N3;

    *di=1;
    *dj=1;
    *dk=1;

  }
  else{
    dualfprintf(fail_file,"No such dir=%d in set_interppoint_loop_ranges_2D_EMF_formerged()\n",corner);
    myexit(9894387);
  }



}

/// range over which need CORN cells have geomcorn defined
/// Note that corner ranges from i=0..N so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
void set_interppoint_loop_ranges_geomcorn_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{
  int *myUconsloop;


  if(interporflux==ENOINTERPTYPE4EMF){
    myUconsloop=emfUconsloop;
  }
  else{
    myUconsloop=Uconsloop;
  }


  if(corner==1){

    *is=0;
    //    *ie=N1-1;
    *ie=N1;

    *js=myUconsloop[FJS]-SHIFT2;
    *je=myUconsloop[FJE]+SHIFT2*2;

    *ks=myUconsloop[FKS]-SHIFT3;
    *ke=myUconsloop[FKE]+SHIFT3*2;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(corner==2){

    *is=myUconsloop[FIS]-SHIFT1;
    *ie=myUconsloop[FIE]+SHIFT1*2;

    *js=0;
    //    *je=N2-1;
    *je=N2;

    *ks=myUconsloop[FKS]-SHIFT3;
    *ke=myUconsloop[FKE]+SHIFT3*2;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(corner==3){

    *is=myUconsloop[FIS]-SHIFT1;
    *ie=myUconsloop[FIE]+SHIFT1*2;

    *js=myUconsloop[FJS]-SHIFT2;
    *je=myUconsloop[FJE]+SHIFT2*2;

    *ks=0;
    //    *ke=N3-1;
    *ke=N3;

    *di=1;
    *dj=1;
    *dk=1;

  }
  else{
    dualfprintf(fail_file,"No such dir=%d in set_interppoint_loop_ranges_geomcorn_formerged()\n",corner);
    myexit(9894387);
  }



}


////////////
/// OBSOLETE
////////////
/// Using set_interppoint_loop_expanded() now for either case
/// Setup loop over region of points
/// very similar to set_interp_loop() in interpline.c
/// off-dir directions have been expanded to account for EMF type calculations
void set_interppoint_loop(int interporflux, int dir, int loc, int continuous, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{

  // interporflux never matters if ghost+active not used

  if(dir==1){

    *is=-SHIFT1;
    *ie=N1-1+SHIFT1;

    *js=INFULL2; // 0
    *je=OUTFULL2; // N2-1;

    *ks=INFULL3; //0;
    *ke=OUTFULL3; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==2){
    *is=INFULL1; //0;
    *ie=OUTFULL1; //N1-1;

    *js=-SHIFT2;
    *je=N2-1+SHIFT2;

    *ks=INFULL3; // 0;
    *ke=OUTFULL3; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==3){
    *is=INFULL1; // 0;
    *ie=OUTFULL1; // N1-1;

    *js=INFULL2; // 0;
    *je=OUTFULL2; // N2-1;

    *ks=-SHIFT3;
    *ke=N3-1+SHIFT3;

    *di=1;
    *dj=1;
    *dk=1;

  }
  else{
    dualfprintf(fail_file,"No such dir=%d in set_interppoint_loop()\n",dir);
    myexit(9894386);
  }


}







/// Setup loop over region of points for finite volume method (or any method that uses extended ghost+active grid)
/// Note that +-SHIFT? gives interpolation at maximal face (i.e. for dir=1 and ncpux1=1, i=0 and i=N1), so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
void set_interppoint_loop_expanded(int interporflux, int dir, int loc, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{
  int *myUconsloop;


  if(interporflux==ENOINTERPTYPE4EMF){
    myUconsloop=emfUconsloop;
  }
  else{
    myUconsloop=Uconsloop;
  }


  if(dir==1){

    // -1 and +1 so can get flux at face from the obtained primitive
    *is=myUconsloop[FIS]-SHIFT1;
    *ie=myUconsloop[FIE]+SHIFT1;

    *js=fluxloop[dir][FJS];
    *je=fluxloop[dir][FJE];

    *ks=fluxloop[dir][FKS];
    *ke=fluxloop[dir][FKE];

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==2){

    *is=fluxloop[dir][FIS];
    *ie=fluxloop[dir][FIE];

    *js=myUconsloop[FJS]-SHIFT2;
    *je=myUconsloop[FJE]+SHIFT2;

    *ks=fluxloop[dir][FKS];
    *ke=fluxloop[dir][FKE];

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==3){

    *is=fluxloop[dir][FIS];
    *ie=fluxloop[dir][FIE];

    *js=fluxloop[dir][FJS];
    *je=fluxloop[dir][FJE];

    *ks=myUconsloop[FKS]-SHIFT3;
    *ke=myUconsloop[FKE]+SHIFT3;

    *di=1;
    *dj=1;
    *dk=1;

  }
  else{
    dualfprintf(fail_file,"No such dir=%d in set_interppoint_loop_expanded()\n",dir);
    myexit(9894387);
  }

}




/// Setup loop over region of points for finite volume method (or any method that uses extended ghost+active grid)
/// Note that +-SHIFT? gives interpolation at maximal face (i.e. for dir=1 and ncpux1=1, i=0 and i=N1), so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
void set_interppoint_loop_expanded_face2cent(int interporflux, int dir, int loc, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{
  int *myUconsloop;


  if(interporflux==ENOINTERPTYPE4EMF){
    myUconsloop=emfUconsloop;
  }
  else{
    myUconsloop=Uconsloop;
  }

  // 0 to N-1 for actual values because others are center boundary cells
  if(dir==1){
    *is=myUconsloop[FIS];
    *ie=myUconsloop[FIE];

    *js=fluxloop[dir][FJS];
    *je=fluxloop[dir][FJE];

    *ks=fluxloop[dir][FKS];
    *ke=fluxloop[dir][FKE];

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==2){

    *is=fluxloop[dir][FIS];
    *ie=fluxloop[dir][FIE];

    *js=myUconsloop[FJS];
    *je=myUconsloop[FJE];

    *ks=fluxloop[dir][FKS];
    *ke=fluxloop[dir][FKE];

    *di=1;
    *dj=1;
    *dk=1;
  }
  else if(dir==3){

    *is=fluxloop[dir][FIS];
    *ie=fluxloop[dir][FIE];

    *js=fluxloop[dir][FJS];
    *je=fluxloop[dir][FJE];

    *ks=myUconsloop[FKS];
    *ke=myUconsloop[FKE];

    *di=1;
    *dj=1;
    *dk=1;

  }
  else{
    dualfprintf(fail_file,"No such dir=%d in set_interppoint_loop_expanded_face2cent()\n",dir);
    myexit(9894387);
  }

}









/// interpolate from a center point to a left/right interface
void slope_lim_point_c2e(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  extern void para(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para2(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para3(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para4(int realisinterp, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void parajon(int i, int j, int k, int loc, int realisinterp, int dir, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void csslope(int pl, FTYPE *y, FTYPE *dq);
  void mp5(FTYPE *y, FTYPE *lout, FTYPE *rout);
  void eppm(FTYPE *y, FTYPE *lout, FTYPE *rout);
  FTYPE mydq = 0.0; //to avoid copying of nan's


  // parabolic reconstruction
  switch(reallim){

  case PARA:
    // para
#if(WHICHPARA==PARA1)
    para(y,left,right);
#elif(WHICHPARA==PARA2)
    para2(y,left,right);
#elif(WHICHPARA==PARA3)
    para3(y,left,right);
#elif(WHICHPARA==PARA4)
    para4(realisinterp, pl,y,left,right);
#elif(WHICHPARA==PARAJON)
    parajon(i,j,k,loc,realisinterp, dir, pl,y,left,right);
#endif
    break;

  case CSSLOPE:
    csslope(pl, y, &mydq);
    *left =y[0] - 0.5* mydq;
    *right=y[0] + 0.5* mydq;
#if(DODQMEMORY)
    // for now force both ways always so have information
    *dq=mydq; // only set to dq if necessary since DODQMEMORY may be 0
#endif
    break;

  case MP5:
    mp5(y,left,right);
    break;

  case EPPM:
    eppm(y,left,right);
    break;

  default:

    //    if(pl==RHO || pl==U1+dir-1 || pl==UU) slope_lim_3points(MINM, y[-1], y[0], y[1], &mydq);
    //    else slope_lim_3points(reallim, y[-1], y[0], y[1], &mydq);
    slope_lim_3points(reallim, y[-1], y[0], y[1], &mydq);

    *left =y[0] - 0.5* mydq;
    *right=y[0] + 0.5* mydq;
#if(DODQMEMORY)
    // for now force both ways always so have information
    *dq=mydq; // only set to dq if necessary since DODQMEMORY may be 0
#endif
    break;
  }


}

/// DEVELOPING RIGHT NOW
/// interpolate from a center point to a left/right interface
void slope_lim_point_e2c_continuous(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  void slope_lim_4points_e2c_continuous(int reallim, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr, FTYPE *dq);
  FTYPE ddq = 0.0; //to avoid copying of nan's
  FTYPE mydq = 0.0; //to avoid copying of nan's

  FTYPE yl,yc,yr,yrr;
  yl=y[-1];
  yc=y[0];
  yr=y[1];
  yrr=y[2];
  if(reallim>MC) reallim=MC; // SUPERGODMARK TODOMARK (just changes local value)
  slope_lim_4points_e2c_continuous(reallim, yl, yc, yr, yrr, &ddq);

  if(reallim<=MC || 1){ // || 1 for now SUPERGODMARK

    // Obtain lower mydq by just setting true Bcenter(x=1/2)=left or =right immediately below and solving for mydq.  So fake slope just to get the final value correct for *right
    // *left can be set, but not actual left value at next cell center.  So really only should use *right
    mydq=2.*(0.25*ddq - 1.*yc + 1.*yl - 0.5*(-1.*yc + yr)); // only for fake left that will just return same Bcenter(x=1/2) as *right
    *left =yc - 0.5* mydq; // fake, just gives back B(x=1/2) as well

    // now overwrite mydq with correct value assuming will only use *right as result for centered final value
    mydq=-2.*(0.25*ddq - 0.5*(-1.*yc + yr));// only for right
    *right=yc + 0.5* mydq; // B(x=1/2)

    //    dualfprintf(fail_file,"POINT4: ijk=%d %d %d : %g %g %g %g : %g %g\n",i,j,k,yl,yc,yr,yrr,mydq,*right);
  }
  else{
    dualfprintf(fail_file,"NOT SETUP YET FOR slope_lim_point_e2c_continuous\n");
    myexit(972372252);
  }

#if(DODQMEMORY)
  // for now force both ways always so have information
  *dq=mydq; // only set to dq if necessary since DODQMEMORY may be 0
#endif
  

}

/// interpolate from a center point to a left/right interface
void slope_lim_point_allpl(int i, int j, int k, int loc, int realisinterp, int dir,int reallim, int startorderi, int endorderi, FTYPE **yreal, FTYPE **y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  extern void parapl(int i, int j, int k, int loc, int realisinterp, int dir, FTYPE **yreal, FTYPE **y, FTYPE *lout, FTYPE *rout);
  void mcsteeppl(int i, int j, int k, int loc, int realisinterp, int dir, FTYPE **yreal, FTYPE **y, FTYPE *lout, FTYPE *rout);

  switch(reallim){
  case PARAFLAT:
    parapl(i, j, k, loc, realisinterp, dir,yreal,y,left,right);
    break;
  case MCSTEEP:
    mcsteeppl(i, j, k, loc, realisinterp, dir,yreal,y,left,right);
    break;
  default:
    dualfprintf(fail_file,"No such allpl reallim=%d\n",reallim);
    myexit(189676);
    break;
  }



}



/// limited slopes using 3 points
/// Gives slope that will be assumed to be located at same location as yc -- center of cell
void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE dq1l,dq1r,dq2;
  extern void get_limit_slopes(int reallim, int extremum, FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dqout);


  dq1l = 1.0*(yc - yl);
  dq1r = 1.0*(yr - yc);
  dq2  = 0.5*(yr - yl);

  get_limit_slopes(reallim, 0, &dq1l, &dq1r, &dq2, dq);

}


/// limited 2nd derivatives using 4 points
/// see e2c_continuous.nb
void slope_lim_4points_e2c_continuous(int reallim, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr, FTYPE *ddq)
{
  FTYPE ddq1l,ddq1r,ddq2;
  extern void get_limit_slopes(int reallim, int extremum, FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dqout);


  ddq1l=0.5*(-2.0*yc+yl+yr);
  ddq1r=0.5*(yc-2.0*yr+yrr);
  ddq2=SIXTH*(-3.*yc + yl + 3.*yr - 1.*yrr) + 
    SIXTH*(3.*yc - 1.*yl - 3.*yr + yrr) + 
    0.5*(-1.*yc + yl - 1.*yr + yrr);
  

  // use same procedure as slope limiter
  get_limit_slopes(reallim, 0, &ddq1l, &ddq1r, &ddq2, ddq);

}


/// dq1 = 1-sided dq2 = centered
void get_limit_slopes(int reallim, int extremum, FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dq)
{
  FTYPE Dqc,Dqp,Dqm,s,theta;


  switch(reallim){

  case MC:
    // monotonized central (Woodward) (Barth-Jespersen)
    Dqm = 2.0 * dq1l[0];
    Dqp = 2.0 * dq1r[0];
    Dqc = dq2[0];
    *dq=MINMODGEN(extremum,Dqc,MINMODGEN(extremum,Dqm,Dqp));
    break;

  case VANL:
    /* van leer slope limiter */
    Dqm = dq1l[0];
    Dqp = dq1r[0];
    s = Dqm * Dqp;
    if (s <= 0.)  *dq= 0.;
    else *dq= (2.0 * s / (Dqm + Dqp));
    break;

  case MINM:
    /* minmod slope limiter (crude but robust) */
    theta = 1.0;
    Dqm = theta * dq1l[0];
    Dqp = theta * dq1r[0];
    Dqc = dq2[0];
    *dq = MINMODGEN(extremum,MINMODGEN(extremum,Dqm,Dqc),Dqp);
    break;

  case NLIMCENT:
    // Centered slope (Fromm)
    *dq = dq2[0];
    break;

  case NLIMUP:
    // Upwind slope (Beam-Warming)
    *dq = dq1l[0];
    break;

  case NLIMDOWN:
    // Downwind slope (Lax-Wendroff)
    *dq = dq1r[0];
    break;

  case DONOR:
    // no slope
    *dq=(0.0);
    break;

  default:
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
    break;
  }

}


void slope_lim_3points_old(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE Dqm, Dqp, Dqc, s;

  if (reallim == MC) {
    Dqm = 2.0 * (yc - yl);
    Dqp = 2.0 * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    s = Dqm * Dqp;
    if (s <= 0.)  *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
        *dq= (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
        *dq= (Dqp);
      else
        *dq= (Dqc);
    }
  }
  /* van leer slope limiter */
  else if (reallim == VANL) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.)
      *dq= 0.;
    else
      *dq= (2.0 * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  else if (reallim == MINM) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.) *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp)) *dq= Dqm;
      else *dq= Dqp;
    }
  }
  else if (reallim == NLIM) {
    Dqc = 0.5 * (yr - yl);
    *dq= (Dqc);
  }
  else if (reallim == DONOR) {
    *dq=(0.0);
  }
  else {
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
  }



}



/// see Mignone & Bodo (2005) astro-ph/0506414 equations 27-31
/// see Colella (1985), Saltzman (1994)
/// second order with 4th order steepeners
void csslope(int pl, FTYPE *y, FTYPE *dq)
{
  FTYPE s, sm, sp;
  FTYPE Dql, Dqlm, Dqlp;
  FTYPE Dqp, Dqpp;
  FTYPE Dqm, Dqmm;
  FTYPE Dqc, Dqcm, Dqcp;
  FTYPE Dqbp,Dqbm;
  FTYPE alpha;


  Dqp=(y[1]-y[0]);
  Dqpp=(y[2]-y[1]);

  Dqm=(y[0]-y[-1]);
  Dqmm=(y[-1]-y[-2]);

  Dqc=0.5*(y[1]-y[-1]);
  Dqcm=0.5*(y[0]-y[-2]);
  Dqcp=0.5*(y[2]-y[0]);

  // Dqmm  Dqm Dqp  Dqpp
  //    Dqcm Dqc Dqcp
  //     sm   q   sp
  //    Dqlm Dql Dqlp
  //    Dqbm     Dqbp

  s=0.5*(sign(Dqp)+sign(Dqm));
  sp=0.5*(sign(Dqpp)+sign(Dqp));
  sm=0.5*(sign(Dqm)+sign(Dqmm));

  // MB05 use alpha=2 for 1-D and alpha=2,1.25,1 for rho,v,p (respectively) for 2D.
  // alpha=[1-2].  alpha=2 more compressive, alpha=1 less compressive.
  if(pl==RHO) alpha=2.0;
  else if((pl>=U1)&&(pl<=U3)) alpha=1.25;
  else if(pl==UU) alpha=1.0;
  else if((pl>=B1)&&(pl<=B3)) alpha=1.25;
  else alpha=2.0;

  Dql=alpha*min(fabs(Dqp),fabs(Dqm));
  Dqlp=alpha*min(fabs(Dqpp),fabs(Dqp));
  Dqlm=alpha*min(fabs(Dqm),fabs(Dqmm));

  Dqbp=sp*min(Dqlp,fabs(Dqcp));
  Dqbm=sm*min(Dqlm,fabs(Dqcm));

  *dq=s*min(fabs(FOURTHIRD*Dqc-SIXTH*(Dqbp+Dqbm)),Dql);


}




/// only test that MCSTEEP fails on is 105 (high relativistic single shock)
/// used when lim=MCSTEEP
void mcsteeppl(int i, int j, int k, int loc, int realisinterp, int dir, FTYPE **yrealpl, FTYPE **ypl, FTYPE *loutpl, FTYPE *routpl)
{
  FTYPE dqpl0[NPR2INTERP][8];
  FTYPE *dqpl[NPR2INTERP];
  FTYPE *y;
  //,*yreal;
  FTYPE *dq;
  FTYPE a_P[10];
  FTYPE *V,*P;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  extern void getPressure(int whicheom, int i, int j, int k, int loc, FTYPE **yrealpl, FTYPE *P);
  extern void parasteep(int dir, int pl, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  extern void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r);
  extern FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P);
  extern void checkparamonotonicity(int smooth, int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout);
  int pl,pliter;
  int dqrange;
  FTYPE Fi;
  void get_mcsteep_dqs(int dqrange, FTYPE *y, FTYPE *dq);
  int mm;
  int odir1,odir2;
  FTYPE a_ddq[7];
  FTYPE *ddq;

  if(EOMRADTYPE!=EOMRADNONE){
    dualfprintf(fail_file,"mcsteppl not setup for koral\n");
    myexit(3752352523);
  }

  // orthogonal to "dir" direction
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;


  ddq=a_ddq+3; // shifted sufficiently


  // consistent with MCSTEEP using 7 points
  dqrange = 5; // dq's will exist from -2,-1,0,1,2 @ CENT and ddq computed from -1,0,1,2 @ FACE

  // shift P
  P=a_P + 4; // P accessed from -3..3 ( shifted sufficiently)


  // assume velocity is istelf
  V = yrealpl[U1+dir-1];

  // shift dq
  PINTERPLOOP(pliter,pl){
    dqpl[pl]=dqpl0[pl]+4; // shifted sufficiently
  }





#if(1 && (CONNECTUANDRHO||CONNECTBANDRHO) )
  if(realisinterp){
    // convert u to (u/rho)
    // When steepening density, if small internal energy, then internal energy will leak out and lead to large u/rho
    // so generically interpolate u/rho if steepening to avoid this
    // if YREALISY==0 or 1, then assume user knows what they are doing in general since YREALISY is not actually generally correct
    // Note that this doesn't change yrealpl since that is stored in separate memory space

    // Problem with this steepening of U is that it can build up a shock and strong dissipation can spontaneously occur in energy equation

    for(mm=-interporder[MCSTEEP]/2;mm<=interporder[MCSTEEP]/2;mm++){
#if(CONNECTUANDRHO)
      ypl[UU][mm] /= (fabs(ypl[RHO][mm])+SMALL);
#endif
#if(CONNECTBANDRHO)
      ypl[B1+odir1-1][mm] /= (fabs(ypl[RHO][mm])+SMALL);
      ypl[B1+odir2-1][mm] /= (fabs(ypl[RHO][mm])+SMALL);
#endif
    }
    // convert input l,r values
#if(CONNECTUANDRHO)
    routpl[UU] /= (fabs(routpl[RHO])+SMALL);
    loutpl[UU] /= (fabs(loutpl[RHO])+SMALL);
#endif
#if(CONNECTBANDRHO)
    routpl[B1+odir1-1] /= (fabs(routpl[RHO])+SMALL);
    loutpl[B1+odir1-1] /= (fabs(loutpl[RHO])+SMALL);
    routpl[B1+odir2-1] /= (fabs(routpl[RHO])+SMALL);
    loutpl[B1+odir2-1] /= (fabs(loutpl[RHO])+SMALL);
#endif

    // GODMARK: Seems like also want to do this for the perpendicular field components since B/rho is advected
    // i.e. Induction equation can be written as:
    // D/Dt(B^j/\rho_0) = (B^i/\rho_0 \nabla_i)v^j
    // Mass must move with orthogonal field, so need to steepen these consistently
    // However, if B^2>>\rho_0, then noise may exist in \rho_0 while want evolution of B^i to be still accurate

  }
#endif

  ///////////////
  //
  // Loop over variables and get interpolated left/right values within cell
  //
  //////////////




  // get pressures for all points since needed for reduction or steepening
#if( DOUSEPARAFLAT || DOUSEPPMCONTACTSTEEP)
  getPressure(EOMSETMHD,i,j,k,loc,yrealpl, P);
#endif

  // KORALTODO: Not using EOMSETRAD yet.



  // computed only once for all variables
#if( DOUSEPARAFLAT )
  Fi = Ficalc(dir,V,P);
#else
  Fi = 0.0;
#endif


  PINTERPLOOP(pliter,pl){

    y=ypl[pl];
    //    yreal=yrealpl[pl];
    dq=dqpl[pl];


    /////////////////////////
    //
    // Get needed slopes
    //
    ////////////////////////

#if(DOINGALLMCSTEEP)
    get_mcsteep_dqs(dqrange, y, dq);
#else

    slope_lim_3points(WHICH3POINTLIMT, y[-1], y[0], y[1],&dq[0]);

#if(DOUSEPPMCONTACTSTEEP)
    slope_lim_3points(WHICH3POINTLIMT, y[-2], y[-1], y[0],&dq[-1]);
    slope_lim_3points(WHICH3POINTLIMT, y[0], y[1], y[2],&dq[1]);
#endif

#if(DOUSEPARAMONOCHECK)
    slope_lim_3points(WHICH3POINTLIMT, y[-3], y[-2], y[-1],&dq[-2]);
    slope_lim_3points(WHICH3POINTLIMT, y[1], y[2], y[3],&dq[2]);
#endif

#endif



    /////////////////////////////
    //
    // get 3-point solution (only need dq[-1,0,1] for parasteep, but need more for paramonocheck with Duez allowance of non-monotonicity
    // get normal centered 3-point result:
    //
    /////////////////////////////

    loutpl[pl] = y[0] - 0.5* dq[0];
    routpl[pl] = y[0] + 0.5* dq[0];


    /////////////////////////////
    //
    // Steepen
    //
    /////////////////////////////


#if(DOUSEPPMCONTACTSTEEP)
    parasteep(dir,pl,V,P,ypl[pl],dqpl[pl],&loutpl[pl],&routpl[pl]);
#endif

    /////////////////////////////
    //
    // Flatten
    //
    /////////////////////////////


#if(DOUSEPARAFLAT)
    paraflatten(dir,pl,ypl[pl],Fi,&loutpl[pl],&routpl[pl]);
#endif


    /////////////////////////////
    //
    // Monotonicity check
    //
    /////////////////////////////

#if(DOUSEPARAMONOCHECK)
    for(mm=-dqrange/2+1;mm<=dqrange/2;mm++){
      ddq[mm] = dq[mm] - dq[mm-1];
    }
    int smooth=0;
    checkparamonotonicity(smooth, dqrange, pl, ypl[pl], ddq, dqpl[pl], &loutpl[pl], &routpl[pl], &loutpl[pl], &routpl[pl]);
#endif


  }



#if(1 && (CONNECTUANDRHO||CONNECTBANDRHO) )
  if(realisinterp){
    // convert (u/rho) back to u
    for(mm=-interporder[MCSTEEP]/2;mm<=interporder[MCSTEEP]/2;mm++){
      // convert back ypl in case use it later for something
#if(CONNECTUANDRHO)
      ypl[UU][mm] *= fabs(ypl[RHO][mm]);
#endif
#if(CONNECTBANDRHO)
      ypl[B1+odir1-1][mm] *= fabs(ypl[RHO][mm]);
      ypl[B1+odir2-1][mm] *= fabs(ypl[RHO][mm]);
#endif
    }
#if(CONNECTUANDRHO)
    // convert back l,r values
    routpl[UU] *= fabs(routpl[RHO]);
    loutpl[UU] *= fabs(loutpl[RHO]);
#endif
#if(CONNECTBANDRHO)
    routpl[B1+odir1-1] *= fabs(routpl[RHO]);
    loutpl[B1+odir1-1] *= fabs(loutpl[RHO]);
    routpl[B1+odir2-1] *= fabs(routpl[RHO]);
    loutpl[B1+odir2-1] *= fabs(loutpl[RHO]);
#endif

  }
#endif





}



/// use this to avoid calling slope_lim_3points() so many times
void get_mcsteep_dqs(int dqrange, FTYPE *y, FTYPE *dq)
{
  int mm;
  FTYPE a_dq1[10],a_dq2[10];
  FTYPE *dq1,*dq2;
  FTYPE Dqp,Dqm,Dqc;
  FTYPE Dqvanl;


  // shifted dq (sufficiently shifted)
  dq1=a_dq1+dqrange;
  dq2=a_dq2+dqrange;


  // get slopes
  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {
    dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face
    dq2[mm] = 0.5 *(y[mm+1]-y[mm-1]); // slope centered at cell center
  }
  mm=dqrange/2+1; // get last dq1 (can't do in loop above since +1 would mean dq2 beyond data range
  dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face

  // Dqm(0) = 2.0*(dq1[0])
  // Dqp(0) = 2.0*(dq1[1]) // so need to go 1 farther to get all needed Dqp's

  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {
    Dqc = dq2[mm];        // normal

    // MC:
#if(WHICH3POINTLIMT==MC)

    Dqm = 2.0 * dq1[mm];   // steepened
    Dqp = 2.0 * dq1[mm+1]; // steepened
    dq[mm]=MCSTEEPMINMOD(Dqc,MCSTEEPMINMOD(Dqm,Dqp));

#elif(WHICH3POINTLIMT==MINM)

    Dqm = dq1[mm];
    Dqp = dq1[mm+1];
    dq[mm] = MCSTEEPMINMOD(Dqm,Dqp);

#elif(WHICH3POINTLIMT==VANL)

    Dqm = dq1[mm];
    Dqp = dq1[mm+1];

    s = Dqm * Dqp;
    if (s <= 0.)  dq[mm] = 0.;
    else dq[mm] = (2.0 * s / (Dqm + Dqp));

#endif
  }


}


#define MIN3(x,y,z) (min(min(x,y),z))
#define MAX3(x,y,z) (max(max(x,y),z))

#define MIN4(w,x,y,z) (min(w,MIN3(x,y,z)))
#define MAX4(w,x,y,z) (max(w,MAX3(x,y,z)))

#define MINMODMP5(x,y) (0.5*(sign(x)+sign(y))*min(fabs(x),fabs(y)))

#define MINMOD4(w,x,y,z) (0.125*(sign(w)+sign(x))*fabs( (sign(w)+sign(y))*(sign(w)+sign(z)) ) * MIN4(fabs(w),fabs(x),fabs(y),fabs(z)) )

#define MEDIAN(x,y,z) ((x) + MINMODMP5((y)-(x),(z)-(x)));


/// http://mesa.sourceforge.net/pdfs/suresh+huynh_97.pdf
/// http://iopscience.iop.org/0264-9381/31/1/015005/pdf/0264-9381_31_1_015005.pdf
/// http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1304.5544
/// see Suresh & Huynh (1997) and Mosta et al. (2014)
/// Use of alpha=4 below may require cour=0.2 but cour=0.4 said to be ok in practice
void mp5(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  void mp5face(FTYPE yll, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr, FTYPE *out);

  FTYPE yll=y[-2];
  FTYPE yl=y[-1];
  FTYPE yc=y[0];
  FTYPE yr=y[1];
  FTYPE yrr=y[2];
  
  mp5face(yll,yl,yc,yr,yrr,rout); // right face relative to center (U^L_{i+1/2})
  mp5face(yrr,yr,yc,yl,yll,lout); // left face relative to center (U^R_{i-1/2})

}

/// if yll,yl,yc,yr,yrr are really in that order, then get lout.
/// reverse order gives rout
/// Takes 5 *point* positions and gives back *point* interface values, while original takes averages and gets points.
void mp5face(FTYPE yll, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr, FTYPE *out)
{

  //  dualfprintf(fail_file,"yll=%g yl=%g yc=%g yr=%g yrr=%g\n",yll,yl,yc,yr,yrr);

  // original UL eq2.1
  FTYPE UL = 0.0078125*(90.*yc - 20.*yl + 3.*yll + 60.*yr - 5.*yrr); // see poly.nb fpr
  //  FTYPE UL = (1.0/60.0)*(2.0*yll-13.0*yl+47.0*yc+27.0*yr-3.0*yrr);

  // Mosta et al. (2014) uses alpha and \tilde{alpha} for same thing and doesn't define \alpha.
  FTYPE alpha=4.0; // see entire paper.
  FTYPE UMP = yc + MINMODMP5(yr-yc,alpha*(yc-yl)); // eq 2.12

  FTYPE epsilonmp5=1E-10; // used by Mosta et al. (2014)
  FTYPE fabsU=(1.0/5.0)*sqrt(yll*yll + yl*yl + yc*yc + yr*yr + yrr*yrr); // L2 norm
  int nolimit=(UL-yc)*(UL-UMP)<=epsilonmp5*fabsU; // eq 2.30


  // unlimited
  *out = UL;

  if(nolimit==0){ // then get limited contribution
    // 2nd derivatives based on points or averages gives same 2nd derivative for point or average reconstruction
    FTYPE Dm = yll - 2.0*yl + yc; // eq 2.19
    FTYPE D0 = yl - 2.0*yc + yr; // eq 2.19
    FTYPE Dp = yc - 2.0*yr + yrr; // eq 2.19

    FTYPE DM4p=MINMOD4(4.0*D0-Dp,4.0*Dp-D0,D0,Dp); // eq 2.27
    FTYPE DM4m=MINMOD4(4.0*D0-Dm,4.0*Dm-D0,D0,Dm); // eq 2.27 also

    FTYPE UUL = yc + alpha*(yc-yr); // eq 2.8
    FTYPE UAV = 0.5*(yc+yr); // eq 2.16
    //    FTYPE UFL = yc + 0.5*(yc-yl);
    //    FTYPE UFR = yr + 0.5*(yr-yrr);
    //FTYPE UMD = median(UAV,UFL,UFR); // replaced by Eq2.27's DM4
    FTYPE UMD = UAV - 0.5*DM4p; // eq 2.28 (correct for points or averages)
    FTYPE dd = 4.0*DM4m;
    FTYPE ULC = 0.75*yc + 0.375*(dd + 2.*yc - 1.*yl) - 0.125*yl; // see poly3-2.nb fnewright (original eq2.29 for averages->points)
    //FTYPE ULC = yc + 0.5*(yc-yl) + (4.0/3.0)*DM4m;

    FTYPE UMIN = max(MIN3(yc,yr,UMD),MIN3(yc,UUL,ULC)); // eq 2.24a
    FTYPE UMAX = min(MAX3(yc,yr,UMD),MAX3(yc,UUL,ULC)); // eq 2.24b

    // new limited value
    *out += MINMOD(UMIN-UL,UMAX-UL); // eq 2.26 using eq 2.5

    //    dualfprintf(fail_file,"UL=%g UMP=%g fabsU=%g nolimit=%d UUL=%g UAV=%g UMD=%g ULC=%g UMIN=%g UMAX=%g\n",UL,UMP,fabsU,nolimit,UUL,UAV,UMD,ULC,UMIN,UMAX);
  }
  else{
    //    dualfprintf(fail_file,"UL=%g UMP=%g fabsU=%g nolimit=%d\n",UL,UMP,fabsU,nolimit);
  }


}


/// Enhanced PPM
///The Einstein Toolkit source code says that the relevant references are:
///Colella & Sekora 2008 (http://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/ColellaSekora.pdf)
///                       McCorquodale & Colella 2011 (http://msp.org/camcos/2011/6-1/camcos-v6-n1-p01-s.pdf)
void eppm(FTYPE *y, FTYPE *lout, FTYPE *rout)
{

}
