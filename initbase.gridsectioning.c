
/*! \file initbase.gridsectioning.c
     \brief General initialization of code related to grid sectioning

// USAGEMARK: on using grid sectioning:
// 1) findandsetactivesection() calls two problem dependent functions.  They specify limits of integration range in terms of absolute grid indicies and the update period.
// 2) Big jumps issue when exposing new cells. If jump is more than MAXBND, then cells are defined by (e.g.) t=0 data instead of being defined by bounds.c because jumped more than MAXBND and bounds only does MAXBND cells: Solution: User should ensure that jumps are <MAXBND or bound region outside standard boundary region in a reasonable way.
// 3) Jon version of OUTFLOW -- field divb problem when *adding* grid cells since outflowed field doesn't satisfy divb=0.  Assume not adding and if do then bounds were written specially so added grid cells will have divb=0.
*/


#include "decs.h"







static int check_limitsinbox(FTYPE rlo, FTYPE rhi, int iloabs, int ihiabs, int doreport);


/// Sets enerregion for ACTIVEREGION
/// Similar to setflux() and other functions, but for ACTIVEREGION/COMPLOOP enerregion.
/// No user/problem-dependent code here
int setgridsectioning(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{

  if(DOGRIDSECTIONING==0){
    dualfprintf(fail_file,"Shouldn't be here in setgridsectioning() with DOGRIDSECTIONING=%d\n",DOGRIDSECTIONING);
    myexit(3139686);
  }

  if(initialcall==1){
    // force update to ACTIVEREGION
    // initialize grid sectioning to full grid at first
    init_gridsectioning();
  }
  else{
    findandsetactivesection(initialcall, timeorder, numtimeorders, thenstep, thetime );
  }

  return(0);
}



/// initialize grid sectioning
/// if not restarting, then just sets to full grid (at end of init.c real first section set)
/// if restarting, then go ahead and set to section given by global_enerregiondef
/// No user/problem-dependent code here
int init_gridsectioning(void)
{
  int doprintout;
  int dimen;
  int badglobal_sectiondef;
  int faketimeorder,fakenumtimeorders;





  if( DOGRIDSECTIONING ){
    trifprintf("Initializing grid sectioning: BEGIN\n");

    faketimeorder=-1;
    fakenumtimeorders=-1;

    if( RESTARTMODE == 0 ) {
      findandsetactivesection(1,faketimeorder,fakenumtimeorders,nstep, t ); //SASMARK SECTIONMARK
    }
    else {
      //set up active section arrays either by using read-in parameters of the active section OR, if this
      //info is absent from restart file, from the current time
      doprintout = 1;
      badglobal_sectiondef=0;
      DIMENLOOP(dimen){
        if(global_enerregiondef[ACTIVEREGION][POINTDOWN][dimen] < - MAXBND || global_enerregiondef[ACTIVEREGION][POINTUP][dimen] > totalsize[dimen] + MAXBND-1 || global_enerregiondef[ACTIVEREGION][POINTDOWN][dimen] >= global_enerregiondef[ACTIVEREGION][POINTUP][dimen]){
          badglobal_sectiondef=1;
        }
      }
      if(badglobal_sectiondef){
        trifprintf( "Sectioning info mangled; regenerating it for current time t = %21.15g\n", t );
        findandsetactivesection(1,faketimeorder,fakenumtimeorders,nstep, t );
      }
      else {
        setactivesection( global_enerregiondef[ACTIVEREGION], doprintout );
      }
    }
    trifprintf("Initializing grid sectioning: END\n");
  }



  return(0);
}


/// bound non-active region fully
/// needed so RK-stepping has defined primitive at intermediate calculations
/// user could use if(WITHINACTIVESECTION(ri,rj,rk)) to avoid boundary cells being set using undefined non-active cells
/// use FULLLOOP to set original boundary cells to correct value so boundary conditions operate nominally
/// No user/problem-dependent code here
/// OPENMPOPTMARK: Code not used (deprecated), so don't optimize
int bound_gridsectioning(int primtype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep)
{
  FTYPE (*primsource)[NSTORE2][NSTORE3][NPR];
  FTYPE (*primdest)[NSTORE2][NSTORE3][NPR];
  int i,j,k,pl,pliter;
  int dimen;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q;
  FTYPE Unew[NPR];





  return(0); // no longer need to bound due to grid sectioning // carefully bound and iterate only those things needed with constrained loops.



  dualfprintf(fail_file,"Should not be in bound_gridsectioning EVER: deprecated code\n");
  myexit(248967233);


  if(DOGRIDSECTIONING==0){
    dualfprintf(fail_file,"Should not be in bound_gridsectioning\n");
    myexit(248967234);
  }

  // global variables to point to in order to get defined primitives for intermediate RK primitives
  // GLOBALMARK
  if(primtype==CENTEREDPRIM){
    primsource=GLOBALPOINT(pglobal);
    primdest=prim;
  }
  else{
    primsource=GLOBALPOINT(pstagglobal);
    primdest=pstag;
  }





  //only for grid sectioning:  to fill in quasi-ghost zones (zones not already filled by bound() call)
  if(primdest != primsource) {

    // Need to set real boundary conditions to something unless user uses WITHINACTIVESECTION(ri,rj,rk) inside their bounds.c
    // This sets non-active region
    if(primtype==STAGGEREDPRIM){
      // Never reach here since only 1 version of pstag=pstagglobal for all RK substeps
      FULLLOOP{
        if(!WITHINACTIVESTAGWITHBNDSECTIONX1(i,j,k)){
          pl=B1; MACP0A1(primdest,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl);  //SASMARK SECTIONING
        }
        if(!WITHINACTIVESTAGWITHBNDSECTIONX2(i,j,k)){
          pl=B2; MACP0A1(primdest,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl);  //SASMARK SECTIONING
        }
        if(!WITHINACTIVESTAGWITHBNDSECTIONX3(i,j,k)){
          pl=B3; MACP0A1(primdest,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl);  //SASMARK SECTIONING
        }
        if(!WITHINACTIVEWITHBNDSECTION(i,j,k)){
          PLOOP(pliter,pl){
            if(pl!=B1 && pl!=B2 && pl!=B3) MACP0A1(primdest,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl);  //SASMARK SECTIONING
          }
        }
      }
    }
    else{
      FULLLOOP{
        if(!WITHINACTIVEWITHBNDSECTION(i,j,k)){
          PLOOP(pliter,pl) MACP0A1(primdest,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl);  //SASMARK SECTIONING
        }
      }
    }
  }



  if(finalstep){
    ///////////
    // GODMARK: not sure if below needed  
    if(FLUXB==FLUXCTSTAG || DOENOFLUX != NOENOFLUX ){
      // then need to deal with unew
      // force unew to be as prim (so second order)
      
      // Note that unlike primitive, conserved quantity is multi-positional so only one thing to modify

      if(FLUXB==FLUXCTSTAG && primtype==STAGGEREDPRIM){

        FULLLOOP{
          if(!WITHINACTIVESTAGWITHBNDSECTIONX1(i,j,k)){
            dimen=1;
            get_geometry(i, j, k, FACE1-1+dimen, ptrgeom);
            pl=B1-1+dimen;
            MACP0A1(ucons,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl)*(ptrgeom->EOMFUNCMAC(pl)); // UEVOLVE
          }
          if(!WITHINACTIVESTAGWITHBNDSECTIONX2(i,j,k)){
            dimen=2;
            get_geometry(i, j, k, FACE1-1+dimen, ptrgeom);
            pl=B1-1+dimen;
            MACP0A1(ucons,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl)*(ptrgeom->EOMFUNCMAC(pl)); // UEVOLVE
          }
          if(!WITHINACTIVESTAGWITHBNDSECTIONX3(i,j,k)){
            dimen=3;
            get_geometry(i, j, k, FACE1-1+dimen, ptrgeom);
            pl=B1-1+dimen;
            MACP0A1(ucons,i,j,k,pl) = MACP0A1(primsource,i,j,k,pl)*(ptrgeom->EOMFUNCMAC(pl)); // UEVOLVE
          }
          // no need to fill ucons[!=B1,B2,B3]
        }
      }// end if prim was staggered

      if(FLUXB==FLUXCTSTAG && primtype==CENTEREDPRIM){
        // unew for field is not set
        FULLLOOP{
          if(!WITHINACTIVEWITHBNDSECTION(i,j,k)){
            get_geometry(i, j, k, CENT, ptrgeom);
            MYFUN(get_state(MAC(primsource,i,j,k), ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
            MYFUN(primtoU(UEVOLVE,MAC(primsource,i,j,k), &q, ptrgeom, Unew, NULL),"initbase.gridsectioning.c:bound_gridsectioning()", "primtoU()", 1);
            PLOOPNOB1(pl) MACP0A1(ucons,i,j,k,pl)=Unew[pl];
            PLOOPNOB2(pl) MACP0A1(ucons,i,j,k,pl)=Unew[pl];
          }
        }
      }// end if prim was centered

      if(FLUXB!=FLUXCTSTAG && primtype==CENTEREDPRIM){
        // unew for field is set
        FULLLOOP{
          if(!WITHINACTIVEWITHBNDSECTION(i,j,k)){
            get_geometry(i, j, k, CENT, ptrgeom);
            MYFUN(get_state(MAC(primsource,i,j,k), ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
            MYFUN(primtoU(UEVOLVE,MAC(primsource,i,j,k), &q, ptrgeom, Unew, NULL),"initbase.gridsectioning.c:bound_gridsectioning()", "primtoU()", 1);
            PLOOP(pliter,pl) MACP0A1(ucons,i,j,k,pl)=Unew[pl];
          }
        }
      }// end if prim was centered

    }
  }// end if finalstep, upon which then only should change conserved quantity that (unew) that is summed over all RK substeps  

  return(0);

}



//// Finds the index of a grid cell that conains a given radius.
//// Should be called after grid is setup for all cpus.
//// \param xr (i) radius
//// \param xi (o) index of a grid cell that contains that radius (relative to current cpu)
//// \param xcpupos1 (o) index of the cpu that contains that radius
//// \return 0 on success
/// Somewhat user-dependent since in general might want to find an arbitrary grid index from an arbitrary physical position
/// Notes:
/// 1) BCtype could change in time, so don't assume fixed -- ok so far
int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi)
{
  int i, j, k, ii;
  FTYPE r1, r2;
  FTYPE X[NDIM],V[NDIM];
  int gotit;
  int fromwhere;


  // need to find horizon and place horizoni on right-hand-side of location

  // definition of horizoni must be consistent so fluxes are consistent and have conservation
  fromwhere=1; // force to be on upside unless Rhor=0, which is caught first

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon
  // notice that only 1 CPU will get horizon since stop process once found
  // notice that radius(horizoni) is below or equal to actual horizon radius


  *xi = FLUXNOTONGRID;
  *xcpupos1 = -1;
  gotit = 0;
  for (ii = numprocs - 1; ii >= 0; ii--) { // should get done by first row
    if (ii == myid) {
      for (i = N1-1; i >= 0; i--) {
        if( BCtype[X2UP] == POLARAXIS ) {
          j = totalsize[2]-1-startpos[2]; //on the polar axis
          k = totalsize[3]-1-startpos[3]; //on the polar axis
        }
        else if( BCtype[X2DN] == POLARAXIS ) {
          j = 0-startpos[2]; //on the polar axis
          k = 0-startpos[3]; //on the polar axis
        }
        else {
          j = N2 / 2;             // doesn't matter (spherical polar assumed)
          k = N3 / 2;             // doesn't matter (spherical polar assumed)
        }
        coord(i, j, k, FACE1, X);
        bl_coord(X, V);
        r1=V[1];
        coord(ip1mac(i), j, k, FACE1, X);
        bl_coord(X, V);
        r2=V[1];
        // looking between FACE1's r value and upper FACE1's r value, so loop is from i=N1-1..i=0

        if(ii==myid && myid==0 && i==0){ //radius to the left of the grid
          // special check in case radius inside inner-most radial grid
          if(xr<=r1){ 
            // then horizon off grid or right on edge, but still ok
            // treat as if horizon is off grid if right on edge
            *xi = 0;
            *xcpupos1=mycpupos[1];
            break;
          }
        }

        if(ii==myid && myid==numprocs-1 && i==N1-1){ //radius to the right of the grid
          // special check in case radius outside outer-most radial grid
          if(xr>=r2){ 
            // then radius off grid or right on edge, but still ok
            // treat as if horizon is off grid if right on edge
            //            *xi = N1-1;
            *xi = N1; // Become N1 so that code knows off grid from FACE
            *xcpupos1=mycpupos[1];
            break;
          }
        }

        if (fromwhere!=2){
          if(xr >= r1 && xr < r2){ // note that if strictly on r2, then next CPU should pick it up
            *xi = i;
            *xcpupos1 = mycpupos[1];
            break;
          }
        }
        else if (fromwhere==2){
          if(xr >= r1 && xr < r2){
            *xi = ip1mac(i);
            *xcpupos1 = mycpupos[1];
            if(*xi>=N1){
              *xi=0;
              ++(*xcpupos1);
            }
            else{
              // then on original CPU
              *xcpupos1 = mycpupos[1];
            }
            break;
          }
        }
      }
    }

    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(xi, 1, MPI_INT, MPIid[ii], MPI_COMM_GRMHD);
      MPI_Bcast(xcpupos1, 1, MPI_INT, MPIid[ii], MPI_COMM_GRMHD);
#endif
    }
    if (*xi >= 0) gotit = 1;                // can stop entire process

    // keep horizoni as relative to CPU with horizon so any CPU knows where horizon is
    //    if (mycpupos[1] != horizoncpupos1) {
    //  horizoni = FLUXNOTONGRID;
    //}                           // reset if not right cpu group
    if (gotit) break;
  }




  if(gotit==0){
    dualfprintf(fail_file,"Never found grid cell corresponding to radius %21.15g : fromwhere=%d\n",xr,fromwhere);
    myexit(6246);
  }


  /////////////////////////////////
  //
  // report some information
#if(PRODUCTION==0)
  if(1||fromwhere==0) {
    trifprintf("xi: %d xcpupos1: %d\n", *xi, *xcpupos1);
    // just a check
    dualfprintf(log_file,"xi: %d mycpupos[1]: %d xcpupos1: %d\n", *xi, mycpupos[1], *xcpupos1);

    trifprintf("end: find_horizon\n");
  }
#endif

  return(0);
}



/// Sets up enerregion when not just initializing.
/// Sets up ACTIVEREGION for grid sectioning
/// Call's user function to get active region boundaries
/// Finally sets-up active region loops
/// No user/problem-dependent code here
/// OPENMPMARK: Assume findandsetactivesection() not called by multiple threads, so ok to  have static and firsttime.
int findandsetactivesection(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  int iloabs, ihiabs;
  int jloabs, jhiabs;
  int kloabs, khiabs;
  int updateeverynumsteps;
  int everynumsteps;
  int doprintout,doset;
  int sectiondef[NUMUPDOWN][NDIM];
  static int firsttimerealgridset=1;


#if( DOGRIDSECTIONING == 0 )
  dualfprintf(fail_file,"Got inside findandsetactivesection() with DOGRIDSECTIONING == 0\n");
  myexit(246983462);
#endif



#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"setgeneral_enerregion() called in parallel region\n");
    myexit(853736);
  }
#endif


  if(initialcall==1){
    // this indicates very first setup call that should (in general) set full grid right now

    // nothing interesting to report for initialization since not even bl_coord() defined yet
    doprintout = 0;

    // normal full total grid
    sectiondef[POINTDOWN][1]=0;
    sectiondef[POINTUP][1]=totalsize[1]-1;
    sectiondef[POINTDOWN][2]=0;
    sectiondef[POINTUP][2]=totalsize[2]-1;
    sectiondef[POINTDOWN][3]=0;
    sectiondef[POINTUP][3]=totalsize[3]-1;

    // compute numcompzones
    compute_numcompzones(sectiondef,&realtotalcompzones);

    setactivesection( sectiondef, doprintout );
    
    return(0); // done!!
  }



  //////////////////////////
  //
  // Update period choices
  //
  //////////////////////////

  // PROBLEM DEPENDENT FUNCTION:
  theproblem_set_enerregionupdate(initialcall, timeorder, numtimeorders, thenstep, thetime, &updateeverynumsteps, &everynumsteps);

  // see if time for update
  doset = (initialcall==-1 || (nstep % updateeverynumsteps==0)  && (timeorder==numtimeorders-1) );
  if(!doset) return(0); // nothing to do

  // see if should print out section information
  // print out if first time in the real setting section (here)
  doprintout = (nstep % everynumsteps==0) && (timeorder==numtimeorders-1) || (firsttimerealgridset==1);
  firsttimerealgridset=0;

  /////////////////////////////////
  //
  // BELOW CALL's PROBLEM DEPENDENT FUNCTION for setting indices of active grid
  //
  /////////////////////////////////

  // Below function should exist in user's init.c that can call existing functions
  theproblem_set_enerregiondef(initialcall, timeorder, numtimeorders, thenstep, thetime, sectiondef);


  /////////////////
  //
  // Set active region (block shape in i,j,k)
  //
  /////////////////



  setactivesection( sectiondef, doprintout );

  // compute numcompzones
  compute_numcompzones(sectiondef,&realtotalcompzones);


  return( 0 );
}


int compute_numcompzones(int (*sectiondef)[NDIM], long long int *localnumcompzones)
{

  *localnumcompzones=(long long int)(sectiondef[POINTUP][1]-sectiondef[POINTDOWN][1]+1)*(long long int)(sectiondef[POINTUP][2]-sectiondef[POINTDOWN][2]+1)*(long long int)(sectiondef[POINTUP][3]-sectiondef[POINTDOWN][3]+1);

  return(0);
}
        



/// Given absolute integer grid positions, determine active loop ranges
/// No user/problem-dependent code here
int setactivesection(int (*sectiondef)[NDIM], int doprintout)
{
  int dimen;
  int updowniter;


  /////////////
  //
  // Check hi/lo indicies
  //
  /////////////
  assert( DOGRIDSECTIONING != 1, "setactivesection(): grid sectioning should be enabled\n" );
  DIMENLOOP(dimen){
    assert( sectiondef[POINTDOWN][dimen] >= sectiondef[POINTUP][dimen], "setactivesection(): hi/lo indices out of order: dimen=%d losectiondef = %d, hisectiondef = %d\n", dimen, sectiondef[POINTDOWN][dimen], sectiondef[POINTUP][dimen] );
  }

  // get region and its additional boundary region
  setgeneral_enerregion(sectiondef, doprintout, ACTIVEREGION, ACTIVEWITHBNDREGION);


  return( 0 );
}












/////////////////////////////////////
///
/// Example user-dependent code
///
/////////////////////////////////////




/// specify grid sectioning for Sasha's wind problem
/// example user-dependent code
int setsashawind_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{
  int cpupos1lo,cpupos1hi;
  int ilo,ihi;
  int iloabs,ihiabs;
  FTYPE rlo,rhi;
  int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
  int extra_safe_cells;
  FTYPE t0,t1;



  extra_safe_cells = MAXBND;
  t_transition_in=1; // was in bounds.c
  
  // Sasha's Wind
  rlo = 0.02 * MAX(0,thetime - t_transition_in);  //SASMARK SECTIONMARK -- hardcoded value of mm (or v^p)
  //this may cause problems if the actual grid inner radius is smaller than 0.1
  if( rlo < 0.1 ) {
    rlo = 0.1;
  }

  //discrete changes in lower radius of active section
  //rlo = pow( 2., floor(log(rlo)/M_LN2l) );

  rhi = 3. * thetime + 3. * (2 * M_PI / 0.25);  //SASMARK SECTIONMARK -- hardcoded value of mm (or v^p)
  if( rhi > Rout ) {
    rhi = Rout;
  }

  //rlo = Rin;
  rhi = Rout;  //on bg don't care where the outer boundary is

  //X1DN boundary of active region
  findindexfromradius( rlo, &cpupos1lo, &ilo );

  //X1UP boundary of active region
  findindexfromradius( rhi, &cpupos1hi, &ihi );

  //find absolute index of the X1DN boundary
  iloabs = cpupos1lo * N1 + ilo;

  //find absolute index of the X1UP boundary
  ihiabs = cpupos1hi * N1 + ihi;

  //add extra cells for safety to ensure shock does not come too close numerically
  //to the X1UP boundary
  ihiabs += extra_safe_cells;

  //make sure the X1UP boundary is on the grid
  if( ihiabs >= totalsize[1] ) {
    ihiabs = totalsize[1] - 1;
  }
  

  // problem Sasha refers to is boundary conditions set worse than on-grid set and kink and behavior at true boundary leads to wave going back
  ihiabs = (totalsize[1] - 1) - MAXBND;  //to avoid problems at the upper boundary

  //iloabs = MAX( iloabs, ihiabs - totalsize[1]/2 );  //don't care for BG (was done on mako so that no more than half grid is in active section)


  enerregiondef[POINTDOWN][1]=iloabs;
  enerregiondef[POINTUP][1]=ihiabs;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
  

  return(0);
}



/// example Sasha's wind setting of update periods
int sashawind_set_enerregionupdate(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps)
{

  ////////////////////
  //
  // Setup update period
  //
  ///////////////////
  *updateeverynumsteps=100;
  //number of steps after which position/size of active section is updated
  *everynumsteps = *updateeverynumsteps;

  return(0);
}



/// specify grid sectioning for Jon's torus problem
/// example user-dependent code
int torus_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{
  int cpupos1lo,cpupos1hi;
  int ilo,ihi;
  int iloabs,ihiabs;
  FTYPE rlo,rhi;
  int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
  int extra_safe_cells;
  FTYPE t0,t1;



  // see advance.c: Whether to allow shift in evolved quantities to preserve conservation and divb=0.  Set to zero if exposing that surface in time.  Set to 1 if absorbing that surface in time and relying on it to inject a solution.
  AVOIDADVANCESHIFTX1DN= 0;
  AVOIDADVANCESHIFTX1UP= 0;
  AVOIDADVANCESHIFTX2DN= 1;
  AVOIDADVANCESHIFTX2UP= 1;
  AVOIDADVANCESHIFTX3DN= 1;
  AVOIDADVANCESHIFTX3UP= 1;
  GLOBALBCMOVEDWITHACTIVESECTION= 1;


  extra_safe_cells = MAXBND;

  // then switch to 6..Rout for active section with decreasing to original at late time
  t0=50;
  t1=150;
  rlo = MAX(Rin-SMALL,6.0*MIN(1.0,1.0-(thetime-t0)/(t1-t0)));


  
  //rlo = Rin;
  rhi = Rout;  //on bg don't care where the outer boundary is

  //X1DN boundary of active region
  findindexfromradius( rlo, &cpupos1lo, &ilo );

  //X1UP boundary of active region
  findindexfromradius( rhi, &cpupos1hi, &ihi );

  //find absolute index of the X1DN boundary
  iloabs = cpupos1lo * N1 + ilo;

  //find absolute index of the X1UP boundary
  ihiabs = cpupos1hi * N1 + ihi;

  //add extra cells for safety to ensure shock does not come too close numerically
  //to the X1UP boundary
  ihiabs += extra_safe_cells;

  //make sure the X1UP boundary is on the grid
  if( ihiabs >= totalsize[1] ) {
    ihiabs = totalsize[1] - 1;
  }


  enerregiondef[POINTDOWN][1]=iloabs;
  enerregiondef[POINTUP][1]=ihiabs;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;


  return(0);

}




/// specify grid sectioning for Jon's torus problem
/// example user-dependent code
int jet_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{
  int cpupos1lo,cpupos1hi;
  int ilo,ihi;
  int iloabs,ihiabs;
  FTYPE rlo,rhi;
  int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
  int extra_safe_cells;
  FTYPE t0,t1;
  FTYPE velgrid;
  int pl,pliter;
  int i,j,k;
  FTYPE newDTref,ftemp;
  FTYPE Rint0,Routt0;
  int reportbox;





#if(0)
  // DEBUG:
  enerregiondef[POINTDOWN][1]=0;
  enerregiondef[POINTUP][1]=88;

  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
  return(0);
#endif




  // see advance.c: Whether to allow shift in evolved quantities to preserve conservation and divb=0.  Set to zero if exposing that surface in time.  Set to 1 if absorbing that surface in time and relying on it to inject a solution.
  // We don't preserve divb inside boundary being absorbed since missing inside-boundary EMF could have cancelled active emf at edge of grid leading to very different (say) B2 just inside boundary that determines inside-boundary divb value.
  AVOIDADVANCESHIFTX1DN= 1;
  AVOIDADVANCESHIFTX1UP= 0;
  AVOIDADVANCESHIFTX2DN= 1;
  AVOIDADVANCESHIFTX2UP= 1;
  AVOIDADVANCESHIFTX3DN= 1;
  AVOIDADVANCESHIFTX3UP= 1;
  GLOBALBCMOVEDWITHACTIVESECTION= 0;








  // setup
  extra_safe_cells = MAXBND;
  t_transition_in=5E3;
  t_transition_out=1.0; // start moving almost immediately
  Rint0=Rin*0.99;
  Routt0=200.0;



  /////////////////////////
  //
  // X1UP
  //
  /////////////////////////
  if(thetime<t_transition_out){

    rhi=Routt0;

    // have to always fix or at least always set EMF's to zero so divb stays zero when exposing solution
    //    BCtype[X1UP]=OUTFLOW;

  }
  else{
    velgrid=1.1; // move slight faster than speed of light so that relativistic flow doesn't ram into outer partially reflecting wall
    rhi = Routt0 + velgrid*MAX(0,thetime - t_transition_out);

    // GODMARK: If were to change rhi (Rout), worry about exposing monopoles.  Is this taken care of?
    
  }


  //  rhi = Rout;  //on bg don't care where the outer boundary is


  /////////////////////////
  //
  // X1DN -- controls dumping too
  //
  /////////////////////////
  if(thetime<t_transition_in){
    rlo=Rint0; // ensure gets ilo=0

    // CHANGINGMARK:
    ftemp=50.0;
    //    ftemp=5.0;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=ftemp;
    DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= DTdumpgen[MAINDUMPTYPE];
    DTdumpgen[ENERDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
    DTdumpgen[IMAGEDUMPTYPE] = DTdumpgen[ENERDUMPTYPE];
    DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];
    DTdumpgen[DEBUGDUMPTYPE] = DTdumpgen[MAINDUMPTYPE];

    BCtype[X1DN]=OUTFLOW;
    
  }
  else{

    // Can just set newDTref=thetime since diagnostics pre-determine next dump time at each dump rather than at every moment, but the above is more robust and keeps simpler times for dumps
    // CHANGINGMARK:
    ftemp=50.0;
    //ftemp=5.0;



    velgrid=0.2;
    rlo = Rint0 + Routt0 + velgrid*MAX(0,thetime - t_transition_in);

    //rlo = Rin;


    //discrete changes in lower radius of active section
    //rlo = pow( 2., floor(log(rlo)/M_LN2l) );

    // steps DT in powers of 10X
    //    newDTref = MAX(thetime,pow((int)log10(fabs(thetime)),10.0));

    //    newDTref = DTdumpgen[MAINDUMPTYPE]*10.0;

    newDTref = round(rlo*ftemp);


    DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = newDTref;
    DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= DTdumpgen[MAINDUMPTYPE];
    DTdumpgen[ENERDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
    DTdumpgen[IMAGEDUMPTYPE] = DTdumpgen[ENERDUMPTYPE];
    DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];
    DTdumpgen[DEBUGDUMPTYPE] = DTdumpgen[MAINDUMPTYPE];

    // change to free outflow?  Allows inflow and outflow.  Always extrapolates.
    //    BCtype[X1DN]=FREEOUTFLOW;
  }




  ////////////////////
  //
  // Get global index from global positions
  //
  ////////////////////


  //X1DN boundary of active region
  findindexfromradius( rlo, &cpupos1lo, &ilo );


  //X1UP boundary of active region
  findindexfromradius( rhi, &cpupos1hi, &ihi );


  //find absolute index of the X1DN boundary
  iloabs = cpupos1lo * N1 + ilo;


  //find absolute index of the X1UP boundary
  ihiabs = cpupos1hi * N1 + ihi;



  /////////////////////
  //
  // check box limits
  //
  /////////////////////
  reportbox=1;
  if(check_limitsinbox(rlo, rhi, iloabs,ihiabs,reportbox)){
    //////////////////////////////////
    // Then need to end computation
    //
    // equivalent to reaching final time
    t=tf;// force to end
  }
  else{
    //////////////////////////////////
    //
    // then still doing computation

    //add extra cells for safety to ensure shock does not come too close numerically
    //to the X1UP boundary
    //  ihiabs += extra_safe_cells;

    //make sure the X1UP boundary is on the grid
    if( ihiabs >= totalsize[1] ) {
      ihiabs = totalsize[1] - 1;
    }

    ////////////////////
    //
    // Enforce relationship between ilo and ihi so that minimum number of cells operated on regardless of how (say) slow outer boundary moves
    //
    ////////////////////
    // Don't use for now, not necessarily needed if smart enough about moving those boundaries
    //  ihiabs = MAX(iloabs+MINGRIDSECTIONCELLSX1,ihiabs);

    if(thetime<t_transition_in){
      iloabs=0; // enforce
      // SUPERGODMARK: inconsistent with setting of horizon and how that affects computations?
    }
  }


  /////////////////////
  //
  // Set final section/region limits in terms of global index
  //
  /////////////////////
  enerregiondef[POINTDOWN][1]=iloabs;
  enerregiondef[POINTUP][1]=ihiabs;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;


  return(0);

}






/// check whether chosen radii and computed indices suggest should stop calculation
int check_limitsinbox(FTYPE rlo, FTYPE rhi, int iloabs, int ihiabs, int doreport)
{


  if(iloabs>=totalsize[1]){
    // >= , since checking if "radius" was within faces.  Then if i=totalsize[1], then beyond last face
    if(doreport) trifprintf("Lower grid section moved beyond outer box face\n");
    return(1);
  }

  if(ihiabs<0){
    // <0 since ==0 is between face at 0 and 1
    if(doreport) trifprintf("Upper grid section moved inside inner box face\n");
    return(1);
  }

  if(iloabs>ihiabs){
    // can be equal, but not beyond eachother
    // when equal, plausible that "radius" sits between both faces so need to compute that single cell
    // need to check "radius" itself to ensure radii haven't moved beyond eachother
    if(doreport) trifprintf("Lower grid section moved beyond upper grid section\n");
    return(1);
  }

  if(rlo>=rhi){
    // then no region to compute
    if(doreport) trifprintf("Lower radius of grid section moved beyond radius of upper grid section\n");
    return(1);
  }


  // if reached here, then ok to continue computing
  return(0);
}













/// specify MPI task rank ordering
/// example user-dependent code
int jet_set_myid(void)
{
  
  if(USEMPI==0){
    return(0); // nothing to do ever
  }

  if(DOGRIDSECTIONING==0){
    return(0); // nothing to do so far as I can tell
  }
  else{

    // Can choose to rearrange MPI tasks

    // Assume TACC Ranger affinity order for MPI tasks by rank.
    // http://services.tacc.utexas.edu/index.php/ranger-user-guide

    // Assume we use "8way" option on Ranger so that tasks are originally:
    // rank0 -> Socket0
    // rank1 -> Socket0
    // rank2 -> Socket1
    // rank3 -> Socket1
    // rank4 -> Socket2
    // rank5 -> Socket2
    // rank6 -> Socket3
    // rank7 -> Socket4

    // We want ranks ordered such that upper-half of physical grid has even ranks and lower half has odd ranks.

    // reorder MPI ranks so that each socket has 2 MPI tasks


    int ranki,rankj,rankk,origid,newid;
    for(rankk=0;rankk<ncpux3;rankk++){
      for(rankj=0;rankj<ncpux2;rankj++){
        for(ranki=0;ranki<ncpux1;ranki++){
          origid=ranki + rankj*ncpux1 + rankk*ncpux1*ncpux2;
          if(ranki<ncpux1/2){
            newid=ranki + rankj*(ncpux1/2) + rankk*(ncpux1/2)*ncpux2;
            // converts 0,1,2,3,4,... -> 0,2,4,6,8,...
            MPIid[origid]=newid*2;
          }
          else if(ranki>=ncpux1/2){
            newid=(ranki-ncpux1/2) + rankj*(ncpux1/2) + rankk*(ncpux1/2)*ncpux2;
            // converts 0,1,2,3,4,... -> 1,3,5,7,9,...
            MPIid[origid]=newid*2+1;
          }
        }
      }
    }// end over rankk


    // Note that unlike TACC Ranger, TACC Lonestar has 2 sockets with 2 cores per socket, but sockets *share* main memory.  So no special socket association is required for optimal memory use.
    // http://services.tacc.utexas.edu/index.php/lonestar-user-guide

 

  }



  return(0);
}
