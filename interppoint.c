
#include "decs.h"




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// POINTS METHODS
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// whether to send all pl's for access by interpolator (such as for steepening and flattening)




// whether to use contact steepener for reallim=MCSTEEP
#define DOUSEPPMCONTACTSTEEP 1 // used alone is unstable
#define DOUSEPARAMONOCHECK 1 // above must be used with this
// using the 2 above, contacts are as steep as with PARA!

#define DOUSEPARAFLAT 1 // avoids oscillations in strong shocks

#define DOINGALLMCSTEEP (DOUSEPPMCONTACTSTEEP&&DOUSEPARAFLAT&&DOUSEPARAMONOCHECK)

#define WHICH3POINTLIMT MC

// below not accounted for unless DOINGALLMCSTEEP==1
// for blast wave (sasha TESTNUMBER==8) actually helps contact but hurts only peak a bit -- contact much better relatively speaking)
#define DQALLOWEXTREMUM 1 // assume USEPARAMONOCHECK enabled if using this!!! (otherwise crazy)


// whether to interpolate and enforce steepening on u and rho in consistent way
// see comments where used first as to why turned off
#define CONNECTUANDRHO (0 && (DOEVOLVERHO&&DOEVOLVEUU))

// whether to interpolate and enforce steepening on B and rho in consistent way
// Should I turn this off?  Maybe for same reason of U version turned off -- may cause current sheet to spontaneously dissipate -- just use separate steepening
// Seems to cause lots of problems with MCSTEEP and MHD torus problem -- too much steepening
#define CONNECTBANDRHO (0 && (DOEVOLVERHO) )


#if(DQALLOWEXTREMUM==0)
#define MCSTEEPMINMOD(a,b) MINMOD(a,b)
#else
#define MCSTEEPMINMOD(a,b) MINMODB(a,b)
#endif


// whether primreal is p2interp
// can't assume this is true in general (e.g. for FIELDCTSTAG input is not even complete range)
//#define YREALISY (RESCALEINTERP==0 || VARTOINTERP==PRIMTOINTERP)




// interpolation with loop over POINTS
void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP])
{
  void slope_lim_point(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right);
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
  int loc;


  loc=CENT; // SUPERGODMARK: Need to pass loc in


  // make sure dq's exist if going to use them
  if( (DODQMEMORY==0)&&(LIMADJUST==LIMITERFIXED) ){
    dualfprintf(fail_file,"Must have dq's when using MC or lower second order methods\n");
    myexit(17);
  }


  // set range of positions interpolated to
  set_interppoint_loop_ranges(interporflux, dir, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);



  {
    int i,j,k;
    // GODMARK: for now assume not doing per-point choice of limiter (note originally pl didn't actually have correct choice of limiter since pl is fake and pl set to final pl at end of loop
    i=j=k=0;
    reallim=choose_limiter(dir, i,j,k,pl);
    // get starting point for stencil, assumed quasi-symmetric (including even/odd size of stencil)
    startorderi = - interporder[reallim]/2;
    endorderi   = - startorderi;
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
          for(l=startorderi;l<=endorderi;l++){
            yreal[l]=MACP0A1(primreal,i + l*idel,j + l*jdel,k + l*kdel,pl);
          }
        }
      
        slope_lim_point(i, j, k, loc, realisinterp, dir, reallim,pl,startorderi,endorderi,yreal,y,
                        &MACP0A1(dq,i,j,k,pl),&MACP0A1(pleft,i,j,k,pl),&MACP0A1(pright,i,j,k,pl)
                        );
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








// given dir=fluxdir, return ranges for loop over which obtain those interpolated quantities
// these is,ie,js,je,ks,ke are used inside a loop that has grid section SHIFTS already, so don't add shift here
// Note that if direction doesn't exist should still return reasonable result that is not out of range (e.g. dir==3 if N3==1 should still give ks=ke=0)
int set_interppoint_loop_ranges(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
{
  
  //  if(useghostplusactive){
  set_interppoint_loop_expanded(interporflux, dir, is, ie, js, je, ks, ke, di, dj, dk);
  //  }
  //  else{
  //    set_interppoint_loop(interporflux, dir, is, ie, js, je, ks, ke, di, dj, dk);
  //  }


  return(0);
}



int set_interppoint_loop_ranges_3Dextended(int interporflux, int *maxis, int *maxie, int *maxjs, int *maxje, int *maxks, int *maxke, int *di, int *dj, int *dk)
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
      set_interppoint_loop_ranges(interporflux, dimen, maxis, maxie, maxjs, maxje, maxks, maxke, di, dj, dk);
    }
    else{
      // max here means largest extent in the direction of that boundary as from interior of domain
      set_interppoint_loop_ranges(interporflux, dimen, &is, &ie, &js, &je, &ks, &ke, di, dj, dk);
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


// range over which need CENT cells defined for EMF in plane defined by corner
// Note that is not like range needed for interpolation itself, which requires for FLUXCTSTAG all off-dir cells for Flux^{dir}
// Note that corner ranges from i=0..N so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
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

// range over which need CORN cells have geomcorn defined
// Note that corner ranges from i=0..N so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
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


///////////
// OBSOLETE
///////////
// Using set_interppoint_loop_expanded() now for either case
// Setup loop over region of points
// very similar to set_interp_loop() in interpline.c
// off-dir directions have been expanded to account for EMF type calculations
void set_interppoint_loop(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
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







// Setup loop over region of points for finite volume method (or any method that uses extended ghost+active grid)
// Note that +-SHIFT? gives interpolation at maximal face (i.e. for dir=1 and ncpux1=1, i=0 and i=N1), so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
void set_interppoint_loop_expanded(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk)
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










// interpolate from a center point to a left/right interface
void slope_lim_point(int i, int j, int k, int loc, int realisinterp, int dir, int reallim, int pl, int startorderi, int endorderi, FTYPE *yreal, FTYPE*y, FTYPE *dq,FTYPE *left,FTYPE *right)
{
  extern void para(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para2(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para3(FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void para4(int realisinterp, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  extern void parajon(int i, int j, int k, int loc, int realisinterp, int dir, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void csslope(int pl, FTYPE *y, FTYPE *dq);
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

// interpolate from a center point to a left/right interface
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


void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE dq1l,dq1r,dq2;
  extern void get_limit_slopes(int reallim, int extremum, FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dqout);


  dq1l = 1.0*(yc - yl);
  dq1r = 1.0*(yr - yc);
  dq2  = 0.5*(yr - yl);

  get_limit_slopes(reallim, 0, &dq1l, &dq1r, &dq2, dq);

}


// dq1 = 1-sided dq2 = centered
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





// see Mignone & Bodo (2005) astro-ph/0506414 equations 27-31
// see Colella (1985), Saltzman (1994)
// second order with 4th order steepeners
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




// only test that MCSTEEP fails on is 105 (high relativistic single shock)

// used when lim=MCSTEEP
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
  extern void getPressure(int i, int j, int k, int loc, FTYPE **yrealpl, FTYPE *P);
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
  getPressure(i,j,k,loc,yrealpl, P);
#endif




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



// use this to avoid calling slope_lim_3points() so many times
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
