
/*! \file initbase.enerregions.c
     \brief General initialization of code related to loops over active cells
*/


#include "decs.h"


// initialize enerregion LOOP: ENERREGIONLOOP()
// No user/problem-dependent code here
void reset_dothisenerregion(int initialcall)
{
  int enerregion;

  // these set some "enerregions"
  // first set which ener regions to do
  ENERREGIONALLLOOP(enerregion) dothisenerreg[enerregion]=1;
  // turn off some enerregions
  if(DOGRIDSECTIONING==0){
    dothisenerreg[ACTIVEREGION]=0;
    dothisenerreg[ACTIVEWITHBNDREGION]=0;
  }
  if(DOJETDIAG==0){
    dothisenerreg[INNERJETREGION]=0;
    dothisenerreg[OUTERJETREGION]=0;
  }
}





// universal function to recompute boundaries for enerregions and ACTIVEREGION
// initialcall: 0 : not initial call, just normal evolution call
// initialcall: 1 : true initial call -- probably first time function called
// intiialcall: 2 : probably not true initial call, but for GRIDSECTION still treat as first call since still need full COMP loops
// initialcall: 3 : not true initial call, but force update to GRIDSECTION
// No user/problem-dependent code here
// Can't assume recompute_fluxpositions() called *every* timestep or substep.  Some steps might be skipped for various reasons to avoid grid changes
int recompute_fluxpositions(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  int normalinitialcall;
  int compinitialcall;

  // put here anything when changing horizonti (apart from recomputing metric, which is done outside)                     
  // if got here then horizonti changed, so need new positions of fluxes

  if(initialcall==2){
    compinitialcall=1;
    normalinitialcall=initialcall;
  }
  else if(initialcall==3){
    compinitialcall=-1;
    normalinitialcall=initialcall;
  }
  else{
    compinitialcall=initialcall;
    normalinitialcall=initialcall;
  }

  if(normalinitialcall==1){
    // (equivalent to horizon not existing -- just normal gridding) 
    horizoni = 0;
    horizoncpupos1 = 0;
  }
  else{
    // otherwise assume set from outside already
  }



  // set enerregions
  // must be done before any ENERREGIONLOOP() call
  reset_dothisenerregion(normalinitialcall);

  // set enerregions
  setflux(normalinitialcall,timeorder,numtimeorders,thenstep,thetime);
  sethorizonflux(normalinitialcall,timeorder,numtimeorders,thenstep,thetime);
  settrueglobalregion(normalinitialcall,timeorder,numtimeorders,thenstep,thetime);

  if(DOJETDIAG)  setjetflux(normalinitialcall,timeorder,numtimeorders,thenstep,thetime);



  // set ACTIVEREGION
  // should come after other definitions of enerregions in case they are used to define activeregion!
  if(DOGRIDSECTIONING){
    setgridsectioning(compinitialcall,timeorder, numtimeorders, nstep,t);
  }


  return(0);

}




// Universal function, where one provides absolute integer grid positions, then it determine ener loop ranges and flux boundaries for provided region and +boundary version if desired
// No user/problem-dependent code here
// OPENMPMARK: Assume setgeneral_enerregion() not called by multiple threads, so ok to  have static and firsttime.
int setgeneral_enerregion(int (*enerregiondef)[NDIM], int doprintout, int whichregion, int whichbndregion)
{
  int dimen;
  int dir;
  int dirsign;
  int enerregion,enerregionglobal;
  int local_enerregiondef[NUMUPDOWN][NDIM];
  FTYPE X[NDIM],V[NDIM];
  int Nvec[NDIM];
  int Nbndvec[NDIM];
  int ti,tj,tk;
  int updowniter;
  int updowniteri,updowniterj,updowniterk;
  int rind,rindglobal;
  int *localenerpos,*localenerposglobal;
  int *localdoflux,*localdofluxglobal;
  static int firsttime=1;
  static int previouslysetglobal_enerregiondef[NUMENERREGIONS];
  int totaldiff;
  int SHIFTdimen[NDIM];


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"setgeneral_enerregion() called in parallel region\n");
    myexit(853736);
  }
#endif

  //////////////////////////////////////////////
  //
  // initialize check if previously set global_enerregiondef[]
  //
  //////////////////////////////////////////////
  if(firsttime){
    ENERREGIONLOOP(enerregion) previouslysetglobal_enerregiondef[enerregion]=0;
    firsttime=0;
  }


  

  //////////////////////
  //
  // setup pointers and perform some checks
  //
  //////////////////////
  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  Nbndvec[0]=0;
  Nbndvec[1]=N1BND;
  Nbndvec[2]=N2BND;
  Nbndvec[3]=N3BND;


  if(whichregion!=NULLENERREGIONS){
    enerregion=whichregion;
    localenerpos=enerposreg[enerregion];  //activesection
    localdoflux=dofluxreg[enerregion];    //activefluxsection

    if(previouslysetglobal_enerregiondef[whichregion]==1){

      // see if anything actually changed
      totaldiff=0;
      DIMENLOOP(dimen){
        for(updowniter=0;updowniter<NUMUPDOWN;updowniter++){
          // GLOBALMARK
          totaldiff += abs(global_enerregiondef[enerregion][updowniter][dimen] - enerregiondef[updowniter][dimen]);
        }
      }

      if(totaldiff>0){
        // check if exposing more than N?BND cells
        DIMENLOOP(dimen){
          updowniter=POINTDOWN;
          if(global_enerregiondef[enerregion][updowniter][dimen] - enerregiondef[updowniter][dimen]>Nbndvec[dimen]){
            // then exposing and too large an exposure for normal bound call
            dualfprintf(fail_file,"SUPERWARNING!!! : Exposing more cells than bounded on inner boundary: enerregion=%d updowniter=%d dimen=%d\n",enerregion,updowniter,dimen);
          }
          updowniter=POINTUP;
          if(enerregiondef[updowniter][dimen]-global_enerregiondef[enerregion][updowniter][dimen]>Nbndvec[dimen]){
            // then exposing and too large an exposure for normal bound call
            dualfprintf(fail_file,"SUPERWARNING!!! : Exposing more cells than bounded on outer boundary: enerregion=%d updowniter=%d dimen=%d\n",enerregion,updowniter,dimen);
          }
        }// dimen loop
      }

    }// if previously set global_enerregiondef[][][]
    else{
      totaldiff=NUMUPDOWN*COMPDIM; // note that there at least should be a change since undefined at first
    }


    // now that checks are complete, set global section definition for restart purposes and future check purposes
    DIMENLOOP(dimen){
      for(updowniter=0;updowniter<NUMUPDOWN;updowniter++){
        // GLOBALMARK
        global_enerregiondef[enerregion][updowniter][dimen] = enerregiondef[updowniter][dimen];
      }
    }

    // now note that previously set global_enerregiondef[][][] so can operate checks
    previouslysetglobal_enerregiondef[whichregion]=1;

  }



  //////////////////////////////////////
  //
  // check if anything to do/change
  //
  //////////////////////////////////////
  if(totaldiff==0){
    // then nothing to do/change, so just quit
    return(0);
  }
  


  //////////////////////
  //
  // setup pointers for +bnd version of grid (only relevant if normal enerregion done)
  //
  //////////////////////
  if(whichbndregion!=NULLENERREGIONS){

    if(whichregion==NULLENERREGIONS){
      dualfprintf(fail_file,"For now assume normal enerregion must exist for bnd region to exist.  This is so enerregiondef[][] remains well-defined, which applies to both versions.\n");
      myexit(896746644);
    }


    enerregionglobal=whichbndregion;
    localenerposglobal=enerposreg[enerregionglobal];  //activewithbndsection
    localdofluxglobal=dofluxreg[enerregionglobal];    //activewithbndfluxsection

    // for restart, recompute_fluxpositions() will redo bnd-dependent regions
  }



  // for enerpos: loop range for this CPU

  // for doflux: flux boundary for this CPU
  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.


  SHIFTdimen[0]=0;
  SHIFTdimen[1]=SHIFT1;
  SHIFTdimen[2]=SHIFT2;
  SHIFTdimen[3]=SHIFT3;


  /////////////////////
  //
  //define relative indices for active section boundaries
  //
  /////////////////////
  DIMENLOOP(dimen){

    // define this CPUs relative section position
    for(updowniter=0;updowniter<NUMUPDOWN;updowniter++){
      local_enerregiondef[updowniter][dimen] = enerregiondef[updowniter][dimen] - startpos[dimen];
    }


    dirsign=-1;
    
    if(whichregion!=NULLENERREGIONS){
      //active section interacts with the current processor
      rind=localenerpos[DIRFROMDIMEN(dimen,dirsign)] = MAX( 0, local_enerregiondef[POINTDOWN][dimen] );
      int rindflux=local_enerregiondef[POINTDOWN][dimen];
      if( Nvec[dimen]>1 && rindflux >= 0 && rindflux < Nvec[dimen] ){
        localdoflux[DIRFROMDIMEN(dimen,dirsign)]=rindflux;
        if(doprintout) logfprintf("proc: %d doing enerregion=%d flux %d rind=%d\n",myid,enerregion,DIRFROMDIMEN(dimen,dirsign),rind);
      }
      else localdoflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
    }
    

    if(whichbndregion!=NULLENERREGIONS){
      rindglobal=localenerposglobal[DIRFROMDIMEN(dimen,dirsign)] = MAX( -Nbndvec[dimen], local_enerregiondef[POINTDOWN][dimen]-Nbndvec[dimen] );
      int rindglobalflux=local_enerregiondef[POINTDOWN][dimen]-Nbndvec[dimen];
      if( Nvec[dimen]>1 && rindglobalflux >= -Nbndvec[dimen] && rindglobalflux < Nvec[dimen]+Nbndvec[dimen] ){
        localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=rindglobalflux;
        if(doprintout) logfprintf("proc: %d doing enerregion=%d sectionflux %d rindglobal=%d\n",myid,enerregionglobal,DIRFROMDIMEN(dimen,dirsign),rindglobal);
      }
      else localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
    }
    
    dirsign=1;

    if(whichregion!=NULLENERREGIONS){
      rind = localenerpos[DIRFROMDIMEN(dimen,dirsign)] = MIN( Nvec[dimen]-1, local_enerregiondef[POINTUP][dimen] );
      int rindflux = local_enerregiondef[POINTUP][dimen];
      //(local_enerregiondef[POINTUP][dimen]+1) is the location of face index
      if( Nvec[dimen]>1 && rindflux >= 0 && rindflux < Nvec[dimen] ){
        localdoflux[DIRFROMDIMEN(dimen,dirsign)]=rindflux + SHIFTdimen[dimen];  //need to add 1 to get upper edge index for the face given rind is cell center index
        if(doprintout) logfprintf("proc: %d doing enerregion=%d flux %d rind=%d\n",myid,enerregion,DIRFROMDIMEN(dimen,dirsign),rind);
      }
      else localdoflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
    }

    if(whichbndregion!=NULLENERREGIONS){
      rindglobal = localenerposglobal[DIRFROMDIMEN(dimen,dirsign)] = MIN( Nvec[dimen]-1+Nbndvec[dimen], local_enerregiondef[POINTUP][dimen]+Nbndvec[dimen] );
      int rindglobalflux = local_enerregiondef[POINTUP][dimen]+Nbndvec[dimen];
      if( Nvec[dimen]>1 && rindglobalflux >= -Nbndvec[dimen] && rindglobalflux < Nvec[dimen]+Nbndvec[dimen] ){
        // localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=rindglobal + SHIFTdimen[dimen];  //need to add 1 to get upper edge index for the face given rindglobal is cell center index
        localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=rindglobalflux; // +1 value isn't set generally by BC's.  flux for this type of quantity not useful once beyond box anyways
        if(doprintout) logfprintf("proc: %d doing enerregion=%d flux %d using rindglobal=%d\n",myid,enerregionglobal,DIRFROMDIMEN(dimen,dirsign),rindglobal);
      }
      else localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
    }

  }




  if(whichregion!=NULLENERREGIONS){
    // check if really on grid
    if(localenerpos[X1UP]<localenerpos[X1DN] && Nvec[1]>1 || localenerpos[X2UP]<localenerpos[X2DN] && Nvec[2]>1 || localenerpos[X3UP]<localenerpos[X3DN] && Nvec[3]>1){
      DIMENLOOP(dimen){
        dirsign=-1;
        localdoflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
        dirsign=1;
        localdoflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
      }
    }
  }
  if(whichbndregion!=NULLENERREGIONS){
    // check if really on grid
    if(localenerposglobal[X1UP]<localenerposglobal[X1DN] && Nvec[1]>1 || localenerposglobal[X2UP]<localenerposglobal[X2DN] && Nvec[2]>1 || localenerposglobal[X3UP]<localenerposglobal[X3DN] && Nvec[3]>1){
      DIMENLOOP(dimen){
        dirsign=-1;
        localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
        dirsign=1;
        localdofluxglobal[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;
      }
    }
  }

  /////////////////////////////////////
  //
  // Print out some diagnostic information about this enerregion
  //
  //////////////////////////////////////

  // only print if user desires and there really was a change
  if(doprintout && totaldiff>0){
    // fluxes are on edges of zone, so 0 and N are on edge fluxes
    if(!specialstep){
      if(whichregion!=NULLENERREGIONS){
        DIRLOOP(dir) logfprintf("proc: %d enerregion=%d: doflux[%d]=%d enerpos[%d]=%d\n",myid,enerregion,dir,localdoflux[dir],dir,localenerpos[dir]);
      }
      if(whichbndregion!=NULLENERREGIONS){
        DIRLOOP(dir) logfprintf("proc: %d enerregion=%d: doflux[%d]=%d enerpos[%d]=%d\n",myid,enerregionglobal,dir,localdofluxglobal[dir],dir,localenerposglobal[dir]);
      }
    }
  }

  ///////////////
  //
  // fluxes are on edges of zone, so 0 and N are on edge fluxes
  //
  ///////////////
  // only print if user desires and there really was a change
  if(doprintout && totaldiff>0){
    if(whichregion!=NULLENERREGIONS){ // only print out non-bnd region
      DIRLOOP(dir) logfprintf("proc: myid=%d  :: t=%21.15g nstep=%ld enerregion=%d section: totaldiff=%d localdoflux[%d]=%d localenerpos[%d]=%d\n",myid,t,nstep,enerregion,totaldiff,dir,localdoflux[dir],dir,localenerpos[dir]);
      
      // full 3D cube outputted (2^3=8 3D points)
      for(updowniteri=NUMUPDOWN-1;updowniteri>=0;updowniteri--) for(updowniterj=NUMUPDOWN-1;updowniterj>=0;updowniterj--) for(updowniterk=NUMUPDOWN-1;updowniterk>=0;updowniterk--){
            ti=enerregiondef[updowniteri][1] + (updowniteri==POINTUP);
            tj=enerregiondef[updowniterj][2] + (updowniterj==POINTUP);
            tk=enerregiondef[updowniterk][3] + (updowniterk==POINTUP);
            // Can't use bl_coord_ijk() or bl_coord_ijk2() below since generally know that i,j,k requested can be beyond stored grid
            bl_coord_coord( ti, tj, tk, CORNT, X, V );
            logfprintf( "t = %21.15g, ud_{i,j,k} = %d %d %d :: CORNT_enerregiondef_{i,j,k} = %d %d %d :: V_{1,2,3} = %21.15g %21.15g %21.15g \n", t, updowniteri, updowniterj, updowniterk, ti, tj, tk, V[1], V[2], V[3] );
     
          }
    }
  }

  return( 0 );
}






// enerregion that is ignorant of horizon -- just total computational domain
int setflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  int enerregion;
  int enerregiondef[NUMUPDOWN][NDIM];
  int doprintout;


  // which enerregion to create
  enerregion=GLOBALENERREGION;

  // set enerregiondef[][]
  setflux_set_enerregiondef(initialcall, timeorder, numtimeorders, thenstep, thetime, enerregiondef );

  // get region and its additional boundary region
  doprintout=!specialstep && initialcall==1;
  setgeneral_enerregion(enerregiondef, doprintout, enerregion, NULLENERREGIONS);

  return(0);
}


// actual definition of enerregion for setflux()
int setflux_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // normal full total grid
  enerregiondef[POINTDOWN][1]=0;
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
  
  return(0);
}





// enerregion that includes horizon
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int sethorizonflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  int enerregion;
  int enerregiondef[NUMUPDOWN][NDIM];
  int doprintout;


  // which enerregion to create
  enerregion=OUTSIDEHORIZONENERREGION;

  // set enerregiondef[][]
  sethorizonflux_set_enerregiondef(initialcall, timeorder, numtimeorders, thenstep, thetime, enerregiondef );

  // get region and its additional boundary region
  doprintout=!specialstep && initialcall==1;
  setgeneral_enerregion(enerregiondef, doprintout, enerregion, NULLENERREGIONS);

  return(0);
}


// determine if this cpu is doing what flux through horizon
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int sethorizonflux_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // entire region outside horizon
  enerregiondef[POINTDOWN][1]=horizoni+horizoncpupos1*N1;
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;

  return(0);
}



// enerregion where quantities are evolved (true global region)
// this can be used to determine where quantities are evolved and so active domain only
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int settrueglobalregion(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  int enerregion,enerbndregion;
  int enerregiondef[NUMUPDOWN][NDIM];
  int doprintout;


  // which enerregion to create
  enerregion=TRUEGLOBALENERREGION;
  // this ranges over entire domain where ANYTHING is done at all
  // This is used to avoid any unnecessary access to cells that will never be used
  enerbndregion=TRUEGLOBALWITHBNDENERREGION;

  // set enerregiondef[][]
  settrueglobalregion_set_enerregiondef(initialcall, timeorder, numtimeorders, thenstep, thetime, enerregiondef );

  // get region and its additional boundary region
  doprintout=!specialstep && initialcall==1;
  setgeneral_enerregion(enerregiondef, doprintout, enerregion, enerbndregion);

  return(0);
}



int settrueglobalregion_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // entire region outside horizon
  // below tries to say that inside horizoni-1-N1BND the evolution doesn't matter so ignore it.
  // used, e.g., in advance.c.  However, so that rest of code doesn't have to change, must ensure that
  //   the other zones are set to something so no Nan's and need to be set such that time-like (e.g. for boundprim())
  // otherwise can go through code and base loops on this enerregion
  //  if(horizoni>0 && horizoncpupos1==0){
  //    // horizoni-1 accounts for fact that horizoni is ceil() of horizon position, not floor()
  //    enerregiondef[POINTDOWN][1]=MAX(0,horizoni-1+horizoncpupos1*N1-N1BND);
  //  }
  //  else{
  // standard case:
  // should be no smaller than 0
  enerregiondef[POINTDOWN][1]=MAX(0,horizoni+horizoncpupos1*N1-N1BND);
  //}

  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;

  return(0);
}





// enerregion for jet parts of grid
// determines the flux positions for each CPU for the jet region (jetedge) AND the range of volume integration for cons. variables in jet (jetpos).
// ISSUEMARK: the below only works for a grid with full Pi.  Crashes code at runtime otherwise!  should fix.
// ISSUEMARK: assume this is an intrinsically axisymmetric function, so k (\phi) is just carried along -- no truncation in \phi
int setjetflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  FTYPE X[NDIM],V[NDIM],r,th;
  int i,j,k,pl,pliter;
  int dir;
  FTYPE startth,endth,thetajet;
  int jetedge[NUMJETS];
  int doprintout;


  doprintout=(!specialstep && initialcall==1);

  if(defcoord==EQMIRROR){
    dualfprintf(fail_file,"setjetflux() not setup to work for non-full-Pi grids\n");
    myexit(1);
  }


  // jet region is assumed to be within a constant theta slice
  // this is theta w.r.t. polar axis
  thetajet=M_PI*0.5-h_over_r_jet;
  // find j for which theta is at our prespecified point

  i=0;j=0;k=0;
  bl_coord_coord(i, j, k, FACE2, X, V);
  startth=V[2];
  i=0;j=N2;k=0;
  bl_coord_coord(i, j, k, FACE2, X, V);
  endth=V[2];

  // assumes 0<thetajet<Pi/2
  if((fabs(startth-thetajet)/thetajet<1E-8)||
     (fabs(endth-thetajet)/thetajet<1E-8)||
     (fabs(startth-(M_PI-thetajet))/thetajet<1E-8)||
     (fabs(endth-(M_PI-thetajet))/thetajet<1E-8)
     ){
    dualfprintf(fail_file,"thetajet is on top of grid, move h_over_r_jet a bit\n");
    myexit(1);
  }


  ////////////////////
  //
  // INNERJET
  //

  // setup pointers
  enerpos=enerposreg[INNERJETREGION];
  doflux=dofluxreg[INNERJETREGION];


  // see if jet edge is related to this CPU
  // assumes increasing j is increasing th
  if((startth<=thetajet)&&(endth<=thetajet)){
    // if cpu entirely within inner theta jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[INNERJET]=FLUXNOTONGRID;
  }
  else if((startth<thetajet)&&(endth>thetajet)){
    // if inner jet edge is on this CPU but not on boundary
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;

    // default:
    enerpos[X2UP]=0;
    jetedge[INNERJET]=FLUXNOTONGRID;

    i=0;
    for(j=0;j<=OUTM2;j++){
      bl_coord_coord(i, j, k, FACE2, X, V);
      r=V[1];
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>thetajet){
        enerpos[X2UP]=jm1mac(j);
        jetedge[INNERJET]=j;
        break;
      }
    }
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if((startth>=thetajet)&&(endth>=thetajet)){
    // if cpu is entirely not contained in inner jet
    enerpos[X1DN]=FLUXNOTONGRID;
    enerpos[X1UP]=FLUXNOTONGRID;
    enerpos[X2DN]=FLUXNOTONGRID;
    enerpos[X2UP]=FLUXNOTONGRID;
    enerpos[X3DN]=FLUXNOTONGRID;
    enerpos[X3UP]=FLUXNOTONGRID;
    jetedge[INNERJET]=FLUXNOTONGRID;
  }
  else{
    trifprintf("problem with INNERJET setjetflux()\n");
    myexit(1);
  }



  // left edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1DN]!=FLUXNOTONGRID)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    if(doprintout) trifprintf("proc: %d doing inner jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=FLUXNOTONGRID;

  // right edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1UP]!=FLUXNOTONGRID)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(doprintout) trifprintf("proc: %d doing inner jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=FLUXNOTONGRID;

  // lower theta boundary
  if((N2>1)&&(enerpos[X2DN]!=FLUXNOTONGRID)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(doprintout) trifprintf("proc: %d doing inner jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=FLUXNOTONGRID;
  
  // upper theta boundary
  if((N2>1)&&(enerpos[X2UP]!=FLUXNOTONGRID)&&(jetedge[INNERJET]!=FLUXNOTONGRID)){ // only get flux if CPU has edge
    doflux[X2UP]=jetedge[INNERJET];
    if(doprintout) trifprintf("proc: %d doing inner jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=FLUXNOTONGRID;

  if((N3>1)&&(enerpos[X3DN]!=FLUXNOTONGRID)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(doprintout) trifprintf("proc: %d doing inner jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=FLUXNOTONGRID;

  // right edge (any directional condition would do)
  if((N3>1)&&(enerpos[X3UP]!=FLUXNOTONGRID)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(doprintout) trifprintf("proc: %d doing inner jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=FLUXNOTONGRID;


  if(doprintout) {
    DIRLOOP(dir) trifprintf("proc: %d %d innerjet: doflux[%d]=%d enerpos[%d]=%d\n",myid,INNERJETREGION,dir,doflux[dir],dir,enerpos[dir]);
  }

  /////////////////////
  //
  // OUTERJET
  //

  // setup pointers
  enerpos=enerposreg[OUTERJETREGION];
  doflux=dofluxreg[OUTERJETREGION];


  // see if outer jet edge is related to this CPU
  if((startth<=M_PI-thetajet)&&(endth<=M_PI-thetajet)){
    // if cpu entirely not within outer jet region
    enerpos[X1DN]=FLUXNOTONGRID;
    enerpos[X1UP]=FLUXNOTONGRID;
    enerpos[X2DN]=FLUXNOTONGRID;
    enerpos[X2UP]=FLUXNOTONGRID;
    enerpos[X3DN]=FLUXNOTONGRID;
    enerpos[X3UP]=FLUXNOTONGRID;
    jetedge[OUTERJET]=FLUXNOTONGRID;
  }
  else if((startth<M_PI-thetajet)&&(endth>M_PI-thetajet)){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;

    // default:
    enerpos[X2DN]=0;
    jetedge[OUTERJET]=FLUXNOTONGRID;
    // if outer jet edge is on this CPU but not on boundary
    i=0;k=0;
    for(j=0;j<=OUTM2;j++){
      bl_coord_coord(i, j, k, FACE2, X, V);
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>M_PI-thetajet){
        enerpos[X2DN]=jm1mac(j);
        jetedge[OUTERJET]=jm1mac(j);
        break;
      }
    }
    enerpos[X2UP]=N2-1;

    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;

    //    dualfprintf(fail_file,"HIT2: %d %d %d %d %d %d\n",enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]);
  }
  else if((startth>=M_PI-thetajet)&&(endth>=M_PI-thetajet)){
    // if cpu is entirely containe within outer jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[OUTERJET]=FLUXNOTONGRID;
  }
  else{
    trifprintf("problem with OUTERJET setjetflux()\n");
    myexit(1);
  }

  if((N1>1)&&(enerpos[X1DN]!=FLUXNOTONGRID)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    if(doprintout) trifprintf("proc: %d doing outer jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=FLUXNOTONGRID;

  if((N1>1)&&(enerpos[X1UP]!=FLUXNOTONGRID)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(doprintout) trifprintf("proc: %d doing outer jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=FLUXNOTONGRID;

  if((N2>1)&&(enerpos[X2DN]!=FLUXNOTONGRID)&&(jetedge[OUTERJET]!=FLUXNOTONGRID)){
    doflux[X2DN]=jetedge[OUTERJET];
    if(doprintout) trifprintf("proc: %d doing outer jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=FLUXNOTONGRID;

  if((N2>1)&&(enerpos[X2UP]!=FLUXNOTONGRID)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(doprintout) trifprintf("proc: %d doing outer jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=FLUXNOTONGRID;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  if((N3>1)&&(enerpos[X3DN]!=FLUXNOTONGRID)&&(mycpupos[3]==0)){
    doflux[X3DN]=0; 
    if(doprintout) trifprintf("proc: %d doing outer jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=FLUXNOTONGRID;

  if((N3>1)&&(enerpos[X3UP]!=FLUXNOTONGRID)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(doprintout) trifprintf("proc: %d doing outer jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=FLUXNOTONGRID;


  if(1||doprintout){
    DIRLOOP(dir) trifprintf("proc: %d %d outerjet: doflux[%d]=%d enerpos[%d]=%d\n",myid,OUTERJETREGION,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}
