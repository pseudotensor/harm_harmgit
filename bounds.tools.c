
#include "decs.h"



#define ADJUSTFLUXCT 0 // whether to adjust fluxCT

// GODMARK: something seriously wrong with OUTEREXTRAP=1 (EOMFFDE)

#define DEBUGINOUTLOOPS 0

#define OUTEREXTRAP 3
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())
// 3: Treat same as HORIZONEXTRAP==3

// how to bound near horizon
#define HORIZONEXTRAP 3
// as above and also:
// 3: jon version

// number of iterations to get v^\phi\sim constant
#define NUMITERVPHI 5


// number of zones to use pole crushing regularizations

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole
#if(PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP==VARTOINTERP)
#define POLEDEATH0 (N2BND==0 ? 0 : MIN(2,N2BND)) 
#else
#define POLEDEATH0 (N2BND==0 ? 0 : MIN(1,N2BND)) // with expansion by 1 point if detects jumps in densities or Lorentz factor (see poldeath())
#endif
//#define MAXPOLEDEATH N2BND // can't be larger than N2BND
#define MAXPOLEDEATH (N2BND==0 ? 0 : MIN(2,N2BND)) // can't be larger than N2BND
#define DEATHEXPANDAMOUNT 0

#define POLEINTERPTYPE 3 // 0=set uu2=bu2=0, 1=linearly interpolate uu2,bu2  2=interpolate B_\phi into pole  3 =linearly for uu2 unless sucking on pole

// number of zones to enforce Lorentz factor to be small
// notice that at pole uu1 and uu2 are artificially large and these regions can lead to runaway low densities and so even higher uu1,uu3
// problem with POLEGAMMADEATH is that at large radius as fluid converges toward pole the fluid stagnates and can fall back at larger angles for no reason -- even for simple torus problem this happens when GAMMAPOLE=1.001
#if(PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP==VARTOINTERP)
#define POLEGAMMADEATH0 0
#else
#define POLEGAMMADEATH0 1
#endif
// maximum allowed Lorentz factor near the pole (set to something large that should be allowed by solution -- problem and grid dependent)
//#define GAMMAPOLE (2.0)


#define GAMMAPOLEOUTGOING 1.1 // keep low
#define GAMMAPOLEOUTGOINGPOWER 1.0
#define GAMMAPOLEOUTGOINGRADIUS 10.0 // very model dependent
#define GAMMAPOLEINGOING GAMMAMAX


// factor by which to allow quantities to jump near pole
#define POLEDENSITYDROPFACTOR 5.0
#define POLEGAMMAJUMPFACTOR 2.0

#if(DOPOLEDEATH)
// avoid such things if N2==1
#define POLEDEATH (N2==1 ? 0 : POLEDEATH0)
#else
#define POLEDEATH 0
#endif

#if(DOPOLEGAMMADEATH)
#define POLEGAMMADEATH (N2==1 ? 0 : POLEGAMMADEATH0)
#else
#define POLEGAMMADEATH 0
#endif

#define POLESMOOTH DOPOLESMOOTH

// whether to average in radius for poledeath
#define AVERAGEINRADIUS 0 // not correct  across MPI boundaries since have to shift near boundary yet need that last cell to be consistent with as if no MPI boundary // OPENNPMARK: Also not correct for OpenMP
#define RADIUSTOSTARTAVERAGING 7 // should be beyond horizon so doesn't diffuse across horizon
#define RADIUSTOAVOIDRADIALSUCK (2.0*Rhor)

#if( POLEDEATH > N2BND )
#error POLEDEATH should be <= N2BND
#endif

#if( MAXPOLEDEATH > N2BND )
#error MAXPOLEDEATH should be <= N2BND
#endif


// X1DN FIXEDUSINGPANALYTIC
//
///////////////////////////
// Currently assume completely general situation where 
// only triggers on BCs, but across all CPUs.  Use grid sectioning to enforce per-CPU dependence if desired.  Any other CPUs that have BCs set will have BCs overwritten by MPI routines
// SUPERGODMARK: Should be able to use set_boundloop()'s result if included FIXED version, but currently it only handles OUTFLOW types for grid sectioning
int bound_x1dn_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int pl,pliter;

#if(ANALYTICMEMORY==0)
  dualfprintf(fail_file,"ANALYTICMEMORY==0 but called bound_x1dn_analytic\n");
  myexit(786752);
#endif


  // then no need to set pglobal or pstagglobal since evolution sets appropriately in constrained way on a well-defined computational box
  dualfprintf(fail_file,"No use for this currently\n");
  myexit(34269834);



  if(ispstag){
    COMPFULLLOOP{
      // note we assume "i=0" is boundary cell to be fixed
      // This ensures divb=0, but may be inconsistent with code's treatement of true i=0 if doing OUTFLOW
      if(WITHINACTIVESTAGBNDSECTIONX1DN(i,j,k)){
	pl=B1; MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
      }
      if(WITHINACTIVEBNDSECTIONX1DN(i,j,k)){
	PLOOP(pliter,pl) if(pl!=B1) MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
      }
    }
  }
  else{
    COMPFULLLOOP{
      if(WITHINACTIVEBNDSECTIONX1DN(i,j,k)){
	PLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(panalytic,i,j,k,pl);
      }
    }
	
  }


  return(0);
}


// X1 inner NSSURFACE (created from OUTFLOW/FIXEDOUTFLOW)
int bound_x1dn_nssurface(
		       int finalstep, int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{
  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);
  FTYPE get_omegaf_code(FTYPE t, FTYPE dt, FTYPE steppart);
  extern FTYPE global_vpar0;
  
  
  
#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    struct of_geom geom[NPR],rgeom[NPR];
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];
    FTYPE X[NDIM], rX[NDIM];
    FTYPE V[NDIM], rV[NDIM];
    FTYPE rV2[NDIM], rV3[NDIM];
    FTYPE bsq;
    FTYPE dxdxp[NDIM][NDIM];
    FTYPE rdxdxp[NDIM][NDIM];
    FTYPE rBur, Bur;
    FTYPE vpar;
    FTYPE rucon[NDIM];
    FTYPE thetapc;
    //FTYPE bval;
    FTYPE BSQORHOBND = 0.5*BSQORHOLIMIT;
    FTYPE BSQOUBND = 0.5*BSQORHOLIMIT;
    
    
    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }
    
    
    if(BCtype[X1DN]==NSSURFACE){
      
      
      if ( (totalsize[1]>1) && (mycpupos[1] == 0) ) { 
	/* inner r boundary condition: u, just copy */
	
	OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	OPENMPBCLOOPBLOCK{
	  OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);
	  
	  
	  
	  ri=riin;
	  rj=j;
	  rk=k;
	  
	  PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
	  pl=B2; bl_coord_ijk(ri,rj,rk,dirprim[pl],rV2);
	  pl=B3; bl_coord_ijk(ri,rj,rk,dirprim[pl],rV3);
	  
	  //outflow everything first
	  LOOPBOUND1INSPECIAL{ // bound entire region inside non-evolved portion of grid
	    PBOUNDLOOP(pliter,pl) {
	      //treat B2 and B3 symmetrically since for a tilted dipole in 3D they become intermixed
	      if (pl==B2) {
		//SASMARK: works accurately only for radial grid
		bl_coord_ijk(i,j,k,dirprim[pl],V);
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (rV2[1]*rV2[1])/(V[1]*V[1]);
	      }
	      else if (pl==B3) {
		bl_coord_ijk(i,j,k,dirprim[pl],V);
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (rV3[1]*rV3[1])/(V[1]*V[1]);
	      }
	      else {
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	      }
	    }
	  }
	  
	  //location of reference staggered radial magnetic field
	  dxdxprim_ijk(ri, rj, rk, FACE1, rdxdxp);
	  bl_coord_ijk(ri, rj, rk, FACE1, rV);
	  rBur = GLOBALMACP0A1(pstagglobal,ri,rj,rk,B1)*rdxdxp[1][1];
	  
	  //fix things: rho, u, B^1, Omega_F
	  if(ispstag){
	    //note that this does not set B^r at r = R_in -- this is enforced through setting emf = 0
	    LOOPBOUND1INSPECIAL{ // bound entire region inside non-evolved portion of grid
	      PBOUNDLOOP(pliter,pl) {
		if(pl==B1) {
#if( NSBC_ASSUME_DIPOLE_FIELD )
		  //This way of setting ghost B[1] relies on the innermost active Bstag[1], which is
		  //evolved according to prescribed emfs (that correspond to rotation with omegaf).
		  //In dipole, radial component scales exactly as 1/r^3, hence
		  //scale radial field as Bur ~ 1/r^3 using the innermost active Bstag[1] as a reference value
		  dxdxprim_ijk(i, j, k, dirprim[pl], dxdxp);
		  bl_coord_ijk(i, j, k, dirprim[pl], V);
		  Bur = rBur*rV[1]*rV[1]*rV[1]/(V[1]*V[1]*V[1]);
		  MACP0A1(prim,i,j,k,pl) = Bur/dxdxp[1][1];
#else
		  //THIS ONLY WORKS IN AXISYMMETRY or if pstaganalytic is updated on every substep
		  MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
#endif
		}
	      }
	    }
	  }//else ispstag
	  else{
	    LOOPBOUND1INSPECIAL{ // bound entire region inside non-evolved portion of grid
	      PBOUNDLOOP(pliter,pl) {
		if( pl==RHO 
		    || pl==UU || pl==U1 
		    || pl==U2 || pl==U3) { //XXX for now outflow everything, worry about velocity later
		  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl); //GLOBALMACP0A1(panalytic,i,j,k,pl);
		}
		else if( pl == B1 ){
#if( NSBC_ASSUME_DIPOLE_FIELD )
		  //scale radial field as Bur ~ 1/r^3
		  dxdxprim_ijk(i, j, k, dirprim[pl], dxdxp);
		  bl_coord_ijk(i, j, k, dirprim[pl], V);
		  Bur = rBur*rV[1]*rV[1]*rV[1]/(V[1]*V[1]*V[1]);
		  MACP0A1(prim,i,j,k,pl) = Bur/dxdxp[1][1];
#else
		  //THIS ONLY WORKS IN AXISYMMETRY or if pstaganalytic is updated on every substep
		  //prim[B1] is centered field; need to interpolate from pstaganalytic[B1]:
		  MACP0A1(prim,i,j,k,pl) = 0.5*(GLOBALMACP0A1(pstaganalytic,i,j,k,pl)+GLOBALMACP0A1(pstaganalytic,i+1,j,k,pl)); //GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
#endif
		}
	      }
	      pl=U1; get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	      bl_coord_ijk(i, j, k, dirprim[pl], V);
	      
	      //bval = MACP0A1(prim,i,j,k,B1);
	      
	      //dualfprintf(fail_file, "bval = %21.15g\n", bval);
	      
	      
	      //compute radial 4-velocity in reference cell
	      pr2ucon(WHICHVEL, MAC(prim,ri,rj,rk), ptrrgeom[pl], rucon);
	      
	      //always do this to avoid noise on switching
	      if( 1 || rucon[1] > 0 ){
		//compute parallel velocity in reference cell
		compute_vpar(MAC(prim,i,j,k), ptrgeom[pl], &vpar);
		
		//if u^r > 0, fix vpar in ghost cells to preselected value
		
		thetapc = sqrt(Rin * get_omegaf_phys(t,dt,steppart));
		
		if( fabs(V[2]) < thetapc || fabs(V[2]-M_PI) < thetapc ){
		  vpar = 0.9;
		}
		if(rucon[1] > 0){
		  vpar = global_vpar0;
		} //else leave it as it is

		//set velocity due to magnetic field rotation + sliding down the field with vpar = global_vpar0
		set_vel_stataxi( ptrgeom[pl], get_omegaf_code(t,dt,steppart), vpar, MACP0A0(prim,i,j,k) );
		
		//set vpar to what we want -- this is to ensure that the correct version of parallel velocity (w.r.t. drift velocity) is used
		//SASMARK: this sometimes leads to superluminal veloicities and hence code crashes
		set_vpar(vpar, ptrgeom[pl], MAC(prim,i,j,k));
		
		if( rucon[1] > 0 ) {
		  //now that velocity is set, can compute bsq
		  if( bsq_calc(MAC(prim,i,j,k), ptrgeom[pl], &bsq) >= 1 ) { 
		    FAILSTATEMENT("bounds.tools.c:bound_x1nd_nssurface()", "bsq_calc()", 2);
		  }
		  
		  //set rho and u consistent with density floor
		  MACP0A1(prim,i,j,k,RHO) = bsq/BSQORHOBND;
		  MACP0A1(prim,i,j,k,UU) = bsq/BSQOUBND;
		}
		//*** the below already done ***
		//else {
		//  //if velocity points into the star, outflow density and internal energy
		//  MACP0A1(prim,i,j,k,RHO) = MACP0A1(prim,ri,rj,rk,RHO);
		//  MACP0A1(prim,i,j,k,UU) = MACP0A1(prim,ri,rj,rk,UU);
		//}
	      }
	    }
	  }//end else ispstag
	  
	}//openmpblock
      }//end if mycpupos[1]==0
    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946);
    }
    
    
  }// end parallel region
  
  
  return(0);
}


// X1UP FIXEDUSINGPANALYTIC
//
///////////////////////////
// Currently assume completely general situation where 
// only triggers on BCs, but across all CPUs.  Use grid sectioning to enforce per-CPU dependence if desired.  Any other CPUs that have BCs set will have BCs overwritten by MPI routines
// SUPERGODMARK: Should be able to use set_boundloop()'s result if included FIXED version, but currently it only handles OUTFLOW types for grid sectioning
int bound_x1up_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int pl,pliter;


#if(ANALYTICMEMORY==0)
  dualfprintf(fail_file,"ANALYTICMEMORY==0 but called bound_x1dn_analytic\n");
  myexit(786753);
#endif



  // then no need to set pglobal or pstagglobal since evolution sets appropriately in constrained way on a well-defined computational box
  dualfprintf(fail_file,"No use for this currently\n");
  myexit(34269834);


  if(ispstag){
    COMPFULLLOOP{
      if(WITHINACTIVESTAGBNDSECTIONX1UP(i,j,k)){
	pl=B1; MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
      }
      if(WITHINACTIVEBNDSECTIONX1UP(i,j,k)){
	PLOOP(pliter,pl) if(pl!=B1) MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(pstaganalytic,i,j,k,pl);
      }
    }
  }
  else{
    COMPFULLLOOP{
      if(WITHINACTIVEBNDSECTIONX1UP(i,j,k)){
	PLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = GLOBALMACP0A1(panalytic,i,j,k,pl);
      }
    }
	
  }


  return(0);
}


// X1 inner OUTFLOW/FIXEDOUTFLOW
int bound_x1dn_outflow(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{



#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    struct of_geom geom[NPR],rgeom[NPR];
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];


    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==FREEOUTFLOW)){


      if ( (totalsize[1]>1) && (mycpupos[1] <= horizoncpupos1)) { // now all CPUs inside CPU with horizon will be using this (GODMARK: reference value needs to be chosen somehow for CPUs not on active grid)
	/* inner r boundary condition: u, just copy */

	OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	OPENMPBCLOOPBLOCK{
	  OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



	  ri=riin;
	  rj=j;
	  rk=k;


#if(HORIZONEXTRAP==0)
	  PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  LOOPBOUND1INSPECIAL{ // bound entire region inside non-evolved portion of grid
	    PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	  }
#elif(HORIZONEXTRAP==1)
	  PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  LOOPBOUND1INSPECIAL{ // bound entire region inside non-evolved portion of grid

	    PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

	    for(pl=RHO;pl<=UU;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]) ;
	    }
	    pl=U1;	    // treat U1 as special
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
	    for(pl=U2;pl<=U3;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. + (i-ri)*dx[1]) ;
	    }
	    pl=B1; // treat B1 special
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]);
	    for(pl=B2;pl<=B3;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. + (i-ri)*dx[1]) ;
	    }
	    for(pl=B3+1;pl<NPRBOUND;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]);
	    }
	  }
#elif(HORIZONEXTRAP==2)
	  get_geometry(ri, rj, rk, dirprim[0], ptrrgeom[0]);
	  rescale(1,1,MAC(prim,ri,rj,rk),ptrrgeom[0],prescale);
	  LOOPBOUND1INSPECIAL{
	    // set guess
	    PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,k,pl);
	    get_geometry(i, j, k, dirprim[0], ptrgeom[0]);	    
	    rescale(-1,1,MAC(prim,i,j,k),ptrgeom[0],prescale);
	  }
#elif(HORIZONEXTRAP==3)
	  extrapfunc(X1DN,j,k,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
#endif




	  if(ispstag==0){

	    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
	      // GODMARK: assume all velocities at same location when doing inflow check
	      LOOPBOUND1INSPECIAL{
#if(WHICHVEL==VEL4)
		get_geometry(i, j, k, dirprim[U1], ptrgeom[U1]);
		inflow_check_4vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1], 0) ;
#elif(WHICHVEL==VEL3)
		get_geometry(i, j, k, dirprim[U1], ptrgeom[U1]);
		inflow_check_3vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1], 0) ;
		// projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
		if(jonchecks){
		  //fixup1zone(MAC(prim,i,j,k),ptrgeom[U1],0);
		  failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),ptrgeom[U1],-3);
		  if(failreturn){
		    dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
		    if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
		  }
		}
#endif
#elif(WHICHVEL==VELREL4)
		get_geometry(i,j,k,dirprim[U1],ptrgeom[U1]) ;
		inflow_check_rel4vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1],0) ;
		if(limit_gamma(GAMMAMAX,MAC(prim,i,j,k),NULL,ptrgeom[U1],0)>=1)
		  FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	      }
	    } // end if not allowing inflow
	  }


	}// end 2 3

      }// end if mycpupos[1]==0



    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946836);
    }


  }// end parallel region


  return(0);
}




// X1 outer OUTFLOW/FIXEDOUTFLOW
int bound_x1up_outflow(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{




#pragma omp parallel  // assume don't require EOS
  {

    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==FREEOUTFLOW)){

      // outer r BC:
      if ( (totalsize[1]>1) && (mycpupos[1] == ncpux1 - 1) ) {


	OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	OPENMPBCLOOPBLOCK{
	  OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


	  ri=riout;
	  rj=j;
	  rk=k;

#if(OUTEREXTRAP==0)
	  PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#elif(OUTEREXTRAP==1)
	  PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  LOOPBOUND1OUT{
	    PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	    for(pl=RHO;pl<=UU;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]) ;
	    }
	    pl=U1; // treat U1 as special
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - 2*(i-ri)*dx[1]) ;
	    for(pl=U2;pl<=U3;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
	    }
	    pl=B1; // treat B1 special
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]) ;
	    for(pl=B2;pl<=B3;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
	    }
	    for(pl=B3+1;pl<NPRBOUND;pl++){
	      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (ptrrgeom[pl]->gdet/ptrgeom[pl]) ;
	    }
	  }
#elif(OUTEREXTRAP==2)
	  get_geometry(ri, rj, rk, dirprim[0], ptrrgeom[0]);
	  rescale(1,1,MAC(prim,ri,rj,rk),ptrrgeom[0],prescale);
	  LOOPBOUND1OUT{
	    // set guess
	    PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,rk,pl);
	    get_geometry(i, j, k, dirprim[0], ptrgeom[0]);
	    rescale(-1,1,MAC(prim,i,j,k),ptrgeom[0],prescale);
	  }
#elif(OUTEREXTRAP==3)
	  extrapfunc(X1UP,j,k,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
#endif





	  if(ispstag==0){

	    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){

	      LOOPBOUND1OUT{
#if(WHICHVEL==VEL4)
		get_geometry(i, j, k, dirprim[U1], ptrgeom[U1]);
		inflow_check_4vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1],0) ;
#elif(WHICHVEL==VEL3)
		get_geometry(i, j, k, dirprim[U1], ptrgeom[U1]);
		inflow_check_3vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1],0) ;
		// projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
		if(jonchecks){
		  //fixup1zone(MAC(prim,i,j,k),ptrgeom[U1],0);
		  failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),ptrgeom[U1],-3);
		  if(failreturn){
		    dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
		    if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
		  }
		}
#endif
#elif(WHICHVEL==VELREL4)
		get_geometry(i,j,k,dirprim[U1],ptrgeom[U1]) ;
		// 	      dualfprintf(fail_file,"JUST BEFORE INFLOWCHECK: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1) *sqrt(geom[U1].gcov[GIND(1,1)]),MACP0A1(prim,i,j,k,U2) *sqrt(geom[U1].gcov[GIND(2,2)]),MACP0A1(prim,i,j,k,U3) *sqrt(geom[U1].gcov[GIND(3,3)]));
		// 	      dualfprintf(fail_file,"JUST BEFORE INFLOWCHECK: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1),MACP0A1(prim,i,j,k,U2),MACP0A1(prim,i,j,k,U3));
		//	      DLOOP(jj,kk) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",jj,kk,geom[U1].gcov[GIND(jj,kk)]);
		inflow_check_rel4vel(1,MAC(prim,i,j,k),NULL,ptrgeom[U1],0) ;
		// 	      dualfprintf(fail_file,"JUST BEFORE LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1) *sqrt(geom[U1].gcov[GIND(1,1)]),MACP0A1(prim,i,j,k,U2) *sqrt(geom[U1].gcov[GIND(2,2)]),MACP0A1(prim,i,j,k,U3) *sqrt(geom[U1].gcov[GIND(3,3)]));
		// 	      dualfprintf(fail_file,"JUST BEFORE LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1),MACP0A1(prim,i,j,k,U2),MACP0A1(prim,i,j,k,U3));
		if(limit_gamma(GAMMAMAX,MAC(prim,i,j,k),NULL,ptrgeom[U1], 0)>=1)
		  FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
		//	      dualfprintf(fail_file,"JUST AFTER LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1) *sqrt(geom[U1].gcov[GIND(1,1)]),MACP0A1(prim,i,j,k,U2) *sqrt(geom[U1].gcov[GIND(2,2)]),MACP0A1(prim,i,j,k,U3) *sqrt(geom[U1].gcov[GIND(3,3)]));
		//	      dualfprintf(fail_file,"JUST AFTER LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,MACP0A1(prim,i,j,k,U1),MACP0A1(prim,i,j,k,U2),MACP0A1(prim,i,j,k,U3));
#endif	
	      }// end over i
	    }// end if not allowing inflow

	  }



	}// end 2 3
      }// end if mycpu is correct

    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946838);
    }

  }// end parallel region


  return(0);
}




// X1 inner R0SING
int bound_x1dn_sym(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{


#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;


  
    if((BCtype[X1DN]==R0SING)||(BCtype[X1DN]==SYMM)||(BCtype[X1DN]==ASYMM) ){



      /* inner radial BC (preserves u^t rho and u) */
      if ( (totalsize[1]>1) && (mycpupos[1] == 0)) {
	////////	LOOPX1dir{

	{ // start block
	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);

	    rj=j;
	    rk=k;
	    LOOPBOUND1IN{
	      PBOUNDLOOP(pliter,pl){
		// SECTIONMARK: assume r=0 singularity can't move
		if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ) ri = -i; // FACE1 values
		else ri=-i-1; // "CENT" values for purposes of this BC
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	      }// over pl
	    }// over boundary zones
	  }
	}// end block



	if((BCtype[X1DN]==R0SING)||(BCtype[X1DN]==ASYMM) ){

	  /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	  ////	LOOPX1dir{

	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


	    // SECTIONMARK: assume r=0 singularity can't move
	    i=0;
	    PBOUNDLOOP(pliter,pl){
	      if(pl==U1 || pl==B1){
		if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ){
		  MACP0A1(prim,i,j,k,pl) = 0.0;
		}
	      }// else don't do this pl
	    } // end over pl
	
	    LOOPBOUND1IN {
	      PBOUNDLOOP(pliter,pl){
		if(pl==U1 || pl==B1){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		}// end if right pl
	      } // end over pl
	    } // end over boundary zones
	  }// end loop 23
	}
      } //end if inner CPU wall
    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946839);
    }

  } // end parallel region

  if(BCtype[X1DN]==R0SING){
    bound_x1dn_r0singfixinterior(boundstage,finalstep, boundtime, whichdir, boundvartype, dirprim,ispstag, prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
  }



  return(0);

}



// X2 inner POLARAXIS
int bound_x2dn_polaraxis_full3d(
				int whichcall,
				int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
				int *inboundloop,
				int *outboundloop,
				int *innormalloop,
				int *outnormalloop,
				int (*inoutlohi)[NUMUPDOWN][NDIM],
				int riin, int riout, int rjin, int rjout, int rkin, int rkout,
				int *dosetbc,
				int enerregion,
				int *localenerpos
		       )
{



#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;



    if(BCtype[X2DN]==POLARAXIS && special3dspc){

      /* inner polar BC (preserves u^t rho and u) */
      if ( (totalsize[2]>1) && (mycpupos[2] == 0) ) {

	if(ncpux3==1){
	  // if ncpux3>1 then this is handled by MPI
	  ////////	  LOOPX2dir{

	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);

	    ri=i;
	    rk=(k+(int)N3/2)%N3;
	    LOOPBOUND2IN{
	      PBOUNDLOOP(pliter,pl){
		// assume polar axis can't move SECTIONMARK
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = -j; // FACE2 values
		else rj=-j-1; // "CENT" values for purposes of this BC
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		// flip sign
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		}

#if(DEBUGINOUTLOOPS)		
		dualfprintf(fail_file,"INNER pole1: ispstag=%d  pl=%d :: ri=%d rj=%d rk=%d i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,i,j,k);
		if(!isfinite(MACP0A1(prim,ri,rj,rk,pl))){
		  dualfprintf(fail_file,"INNER pole1: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
		}
#endif


	      }// end over pl
	    }// end over j

	    // also do j=0
	    j=0;
	    PBOUNDLOOP(pliter,pl){
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		// only do j=0 if at pole exactly
		rj = -j; // FACE2 values
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	  
		// flip sign
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		}

#if(DEBUGINOUTLOOPS)		
		dualfprintf(fail_file,"INNER pole2: ispstag=%d  pl=%d :: ri=%d rj=%d rk=%d i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,i,j,k);
		if(!isfinite(MACP0A1(prim,ri,rj,rk,pl))){
		  dualfprintf(fail_file,"INNER pole2: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
		}
#endif


	      }
	    }

	  }// end over i,k
	}// end if ncpux3==1

	// SUPERGODMARK: continue to use for now
	// only do poledeath() after MPI call (i.e. whichcall==2)
	if(BCtype[X2DN]==POLARAXIS && (whichcall==2 && ncpux3>1 || whichcall==1 && ncpux3==1) ){
	  if(POLEDEATH) poledeath(X2DN,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	  else if(POLESMOOTH) polesmooth(X2DN,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	}


      }//end if inner CPU wall





    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946840);
    }



  } // end parallel region




  //end if special transmissive BCs at poles
  return(0);
}





// X2 inner POLARAXIS
int bound_x2dn_polaraxis(
			 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
			 int *inboundloop,
			 int *outboundloop,
			 int *innormalloop,
			 int *outnormalloop,
			 int (*inoutlohi)[NUMUPDOWN][NDIM],
			 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
			 int *dosetbc,
			 int enerregion,
			 int *localenerpos
		       )
{



#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;


    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){


      if (  (totalsize[2]>1) && (mycpupos[2] == 0) ){

	/////      LOOPX2dir{

	{// block region
	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);
	
	    ri=i;
	    rk=k;
	    LOOPBOUND2IN{
	      PBOUNDLOOP(pliter,pl){
		// assume polar axis can't move SECTIONMARK
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = -j; // FACE2 values
		else rj=-j-1; // "CENT" values for purposes of this BC
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	      }
	    }
	  }
	}// end blocked region


	if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

	  /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	  ////	LOOPX2dir{

	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


	    // assume polar axis can't move SECTIONMARK
	    j=0;
	    PBOUNDLOOP(pliter,pl){
	      if(pl==U2 || pl==B2){
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		  MACP0A1(prim,i,j,k,pl) = 0.0;
		}
	      }// else don't do this pl
	    } // end over pl
	
	
	    LOOPBOUND2IN {
	      PBOUNDLOOP(pliter,pl){
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		} // end if right pl
	      } // end over pl
	    } // end over boundary cells
	  }// end loop 13


	}// end if polar or asym condition

	if(BCtype[X2DN]==POLARAXIS){
	  if(POLEDEATH) poledeath(X2DN,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	  else if(POLESMOOTH) polesmooth(X2DN,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	}
	
      }// end if inner CPU wall



    }// end if polar asym or asym
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946841);
    }

  }// end parallel region

    

  return(0);
}





// X2 outer POLARAXIS full3d
int bound_x2up_polaraxis_full3d(
				int whichcall,
				int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
				int *inboundloop,
				int *outboundloop,
				int *innormalloop,
				int *outnormalloop,
				int (*inoutlohi)[NUMUPDOWN][NDIM],
				int riin, int riout, int rjin, int rjout, int rkin, int rkout,
				int *dosetbc,
				int enerregion,
				int *localenerpos
		       )
{





#pragma omp parallel  // assume don't require EOSs
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;



    if(BCtype[X2UP]==POLARAXIS && special3dspc){
	

      if (  (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ) {


	if(ncpux3==1){
	  // if ncpux3>1 then this is handled by MPI

	  //////	  LOOPX2dir{
	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


	    ri=i;
	    rk=(k+(int)N3/2)%N3;
	    LOOPBOUND2OUT{
	      PBOUNDLOOP(pliter,pl){
		// assume polar axis can't move SECTIONMARK
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = 2*N2-j; // FACE2 values
		else rj=2*(N2-1)-j+1; // "CENT" values for purposes of this BC
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		// flip sign
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		}

#if(DEBUGINOUTLOOPS)		
		dualfprintf(fail_file,"OUTER pole1: ispstag=%d  pl=%d :: ri=%d rj=%d rk=%d i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,i,j,k);
		if(!isfinite(MACP0A1(prim,ri,rj,rk,pl))){
		  dualfprintf(fail_file,"OUTER pole1: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
		}
#endif


	      }// end over pl
	    }// end over j

	    // also do j=N2
	    j=N2;
	    PBOUNDLOOP(pliter,pl){
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		// only do j=0 if at pole exactly
		rj = 2*N2-j; // FACE2 values
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	  
		// flip sign
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		}

#if(DEBUGINOUTLOOPS)		
		dualfprintf(fail_file,"OUTER pole2: ispstag=%d  pl=%d :: ri=%d rj=%d rk=%d i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,i,j,k);
		if(!isfinite(MACP0A1(prim,ri,rj,rk,pl))){
		  dualfprintf(fail_file,"OUTER pole2: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
		}
#endif

	      }
	    }

	  }// end over i,k
	}// end if ncpux3==1

	// SUPERGODMARK: continue to use for now
	// only do poledeath() after MPI call (i.e. whichcall==2)
	if(BCtype[X2UP]==POLARAXIS && (whichcall==2 && ncpux3>1 || whichcall==1 && ncpux3==1) ){
	  if(POLEDEATH) poledeath(X2UP,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	  else if(POLESMOOTH) polesmooth(X2UP,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	}

      }// end if outer CPU wall

    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946842);
    }

  }// end parallel region
    

  return(0);
}





// X2 outer POLARAXIS
int bound_x2up_polaraxis(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{





#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;




  
    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){


      if (  (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ) {

	//////      LOOPX2dir{
	{// block region
	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


	    ri=i;
	    rk=k;
	    LOOPBOUND2OUT{
	      PBOUNDLOOP(pliter,pl){
		// assume polar axis can't move SECTIONMARK
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = 2*N2-j; // FACE2 values
		else rj=2*(N2-1)-j+1; // "CENT" values for purposes of this BC
		MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	      }
	    }
	  }
	}
    
	if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

	  /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	  //////	LOOPX2dir{
	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
	    OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


	    j=N2;
	    PBOUNDLOOP(pliter,pl){
	      if(pl==U2 || pl==B2){
		if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		  MACP0A1(prim,i,j,k,pl) = 0.0;
		}
	      }// else don't do this pl
	    } // end over pl

	    LOOPBOUND2OUT {
	      PBOUNDLOOP(pliter,pl){
		if(pl==U2 || pl==B2){
		  MACP0A1(prim,i,j,k,pl) *= -1.;
		} // end if right pl
	      } // end over pl
	    } // end over bondary cells
	  }// end loop 13

	}// end if reflecting type conditions


	if(BCtype[X2UP]==POLARAXIS){
	  if(POLEDEATH) poledeath(X2UP,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	  else if(POLESMOOTH) polesmooth(X2UP,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
	}
      
      }// end if mycpupos[2]==ncpux2-1
      
  

    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946843);
    }

  }// end parallel region


  return(0);
}




// X3 inner periodic
int bound_x3_periodic(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{



#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    int locali1,globali1;
    int locali2,globali2;
    int ri1,ri2;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int horizonti;
    int jj,kk;


    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){


      // periodic x3
      if ( (totalsize[3]>1) && (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
	// just copy from one side to another

	////	LOOPX3dir{

	OPENMPBCLOOPVARSDEFINELOOPX3DIR; OPENMPBCLOOPSETUPLOOPX3DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	OPENMPBCLOOPBLOCK{
	  OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j);


	  // copy from upper side to lower boundary zones
	  ri=i;
	  rj=j;
	  rk=rkout;
	  LOOPBOUND3IN PBOUNDLOOP(pliter,pl){
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+1+k,pl);

	    
#if(DEBUGINOUTLOOPS)		
	    dualfprintf(fail_file,"INNER X3: ispstag=%d  pl=%d :: ri=%d rj=%d rk=%d (%d) i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,rk+SHIFT3+k,i,j,k);
	    if(!isfinite(MACP0A1(prim,ri,rj,rk+SHIFT3+k,pl))){
	      dualfprintf(fail_file,"INNER X3: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
	    }
#endif

	  }

	  // copy from lower side to upper boundary zones
	  ri=i;
	  rj=j;
	  rk=rkin;
	  LOOPBOUND3OUT PBOUNDLOOP(pliter,pl){
	    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(k-N3),pl);

#if(DEBUGINOUTLOOPS)		
	    dualfprintf(fail_file,"OUTER X3: ispstag=%d pl=%d :: ri=%d rj=%d rk=%d (%d) i=%d j=%d k=%d\n",ispstag,pl,ri,rj,rk,rk+(k-N3),i,j,k);
	    if(!isfinite(MACP0A1(prim,ri,rj,rk+(k-N3),pl))){
	      dualfprintf(fail_file,"INNER X3: caught copying nan ri=%d rj=%d rk=%d pl=%d\n",ri,rj,rk,pl);
	    }
#endif

	  }
	}
      }


    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946844);
    }


  }// end parallel region

  return(0);
}




// X1 inner R0SING
int bound_x1dn_r0singfixinterior(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{
  int i,j,k,pl,pliter;
  int locali1,globali1;
  int locali2,globali2;
  int ri1,ri2;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  int horizonti;
  int jj,kk;



  if(  (totalsize[1]>1) && BCtype[X1DN]==R0SING){ // global condition that all CPUs know about

    /////////////////
    //
    // now deal with region interior to the horizon
    // This just makes values reasonable for diagnostics to be ignorant of loop ranges -- results are NOT used! -- no FLUXCTSTAG issues
    // similar to bound_spacetime_inside_horizon() in metric.c

  
    horizonti=N1*horizoncpupos1+horizoni;

    if(horizonti>0){ // global condition that all CPUs know about
      // for self-gravity, only loop over region outside horizon
      enerregion=OUTSIDEHORIZONENERREGION;
      enerpos=enerposreg[enerregion];


      // copy from very outer boundary or from location on grid defined by where last ghost cell is inside horizon

      // this doesn't do anything to cells inside horizon but outside fake region
      // apparently in this case uu1>0 even though inside the horizon
      globali1=N1*horizoncpupos1+horizoni-N1BND;
      if(globali1<startpos[1]-N1BND) locali1=-N1BND;
      else if(globali1>endpos[1]+N1BND-1) locali1=N1+N1BND-1;
      else locali1=horizoni-N1BND;
      ri1=locali1;



      // include all zones inside horizon as outflow so remove any peak that was accreted as black hole
      // flux through horizon of mass should be consistent with either method, but seems more reasonable to do below
      globali2=N1*horizoncpupos1+horizoni;
      if(globali2<startpos[1]) locali2=0;
      else if(globali2>endpos[1]-1) locali2=N1-1;
      else locali2=horizoni;
      ri2=locali2;
      // GODMARK: somehow leads to density sticking to value despite much lower active density


      LOOPF{

	rj=j;
	rk=k;
      
	//      if(WITHINENERREGION(enerpos,i,j,k)){
	// then do nothing
	//      }
	//else{
	// assume horizon on negative side of "i" so don't modify right side of "i" that happens to be in the outer radial boundary

	if(i < globali1-startpos[1]){
	  // then inside horizon
	  
	  // set primitives to be trivialized
 
	  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri1,rj,rk,pl);

	  // reset pflag in case old failure inside horizon
#if(UTOPRIMADJUST==UTOPRIMAVG)
	  GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
#endif


	  if(ispstag==0){
	    // making it small artificially forces timestep to be small since timestep still limited within fake boundary region
	    //	  MACP0A1(prim,i,j,k,RHO) = SMALL;
	    // assumed acting on relative 4-velocity that stays physical no matter what
	    if(MACP0A1(prim,i,j,k,U1)>0) MACP0A1(prim,i,j,k,U1)=-0.5*(fabs(MACP0A1(prim,ri1,rj,rk,U1))+fabs(MACP0A1(prim,ri1+1,rj,rk,U1)));
	  }
	}

	if(0 && i < globali2-startpos[1]){
	  // then on or inside horizon and inside N1BND more grid cells for interpolation
	  
	  if(ispstag==0){
	    // assumed acting on relative 4-velocity that stays physical no matter what
	    if(MACP0A1(prim,i,j,k,U1)>0) MACP0A1(prim,i,j,k,U1)=-0.5*(fabs(MACP0A1(prim,ri1,rj,rk,U1))+fabs(MACP0A1(prim,ri1+1,rj,rk,U1)));
	  }

	}
	//}
      }


    }
  }
  else{
    dualfprintf(fail_file,"Shouldn't be here in bounds\n");
    myexit(3946845);
  }



  return(0);
}







// Check that boundary conditions were set properly
int bound_checks1(
		       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
		       int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       )
{
  int i,j,k,pl,pliter;
  int locali1,globali1;
  int locali2,globali2;
  int ri1,ri2;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  int horizonti;
  int jj,kk;



#if(PRODUCTION==0)
  // can only check after all directions done, and we assume whichdir==3 last to be done
  int trigger=0;
  if(ispstag){
    ZSLOOP(inoutlohi[POINTDOWN][POINTDOWN][1], inoutlohi[POINTUP][POINTUP][1],inoutlohi[POINTDOWN][POINTDOWN][2], inoutlohi[POINTUP][POINTUP][2],inoutlohi[POINTDOWN][POINTDOWN][3], inoutlohi[POINTUP][POINTUP][3]){
      // Can't use FULLLOOP since different boundary types go different depths into boundary cells
      //      PALLLOOP(pl){
      PBOUNDLOOP(pliter,pl){
	// for staggered field, avoid inner-most cell since isn't needed or computed and so can be NaN
	if(pl==B1 && i>=-N1BND+SHIFT1 || pl==B2 && j>=-N2BND+SHIFT2 || pl==B3 && k>=-N3BND+SHIFT3){
	  if(!finite(MACP0A1(prim,i,j,k,pl))){
	    trigger++;
	    dualfprintf(fail_file,"whichdir=%d ispstag=%d trigger=%d :: BC didn't set properly: #1: i=%d j=%d k=%d ti=%d tj=%d tk=%d pl=%d\n",whichdir,ispstag,trigger,i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
	  }
	}
      }
    }
  }
  else{
    ZSLOOP(inoutlohi[POINTDOWN][POINTDOWN][1], inoutlohi[POINTUP][POINTUP][1],inoutlohi[POINTDOWN][POINTDOWN][2], inoutlohi[POINTUP][POINTUP][2],inoutlohi[POINTDOWN][POINTDOWN][3], inoutlohi[POINTUP][POINTUP][3]){
      //      PALLLOOP(pl){
      PBOUNDLOOP(pliter,pl){
	if(!finite(MACP0A1(prim,i,j,k,pl))){
	  trigger++;
	  dualfprintf(fail_file,"whichdir=%d ispstag=%d trigger=%d :: BC didn't set properly: #1: i=%d j=%d k=%d ti=%d tj=%d tk=%d pl=%d\n",whichdir,ispstag,trigger,i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
	}
      }
    }
  }
  //  if(trigger>0) myexit(92489346); // first calls to bound will have no field, for example, so can't exit
#endif



  return (0);
}





// OPTMARK: all get_geometries, etc. could be simplified if take advantage of fact that only really 2 cases: all CENT or CENT for all except B's.
// otherwise extrapfunc() and poledeath()'s get_geometry() dominates CPU time when few active cells or many boundary cells

// boundary = X1DN or X1UP
// assume if ispstag==1 then only doing field part -- otherwise logic would get more complicated
int extrapfunc(int boundary, int j,int k,
	       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
	       int *inboundloop,
	       int *outboundloop,
	       int *innormalloop,
	       int *outnormalloop,
	       int (*inoutlohi)[NUMUPDOWN][NDIM],
	       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
	       int *dosetbc,
	       int enerregion,
	       int *localenerpos
	       )
{
  int i,pl,pliter;
  int ri,rj,rk;
  int iloopstart,iloopend,iloopstep;
  int ri2,ri3;
  // for Bd3 copy
  FTYPE Bd3,Bd3ri2,Bd3ri3;
  int jj;
  FTYPE Bu1,Bu2,gcon03,gcon13,gcon23,gcon33;
  FTYPE gcov01,gcov02,gcov11,gcov12,gcov21,gcov22,gcov03,gcov13,gcov23;
  // Mdot copy
  FTYPE Mdot,ucon[NDIM],uconref[NDIM];
  FTYPE Xr[NPR][NDIM],Vr[NPR][NDIM],X[NPR][NDIM],V[NPR][NDIM];
  FTYPE primtemp[NPR];
  FTYPE *ucontouse, ucon2[NDIM];
  int itervphi;
  FTYPE uconrefri2[NDIM];
  FTYPE Mdotri2;
  FTYPE uconrefri3[NDIM];
  FTYPE othersref[NUMOTHERSTATERESULTS],othersrefri2[NUMOTHERSTATERESULTS],othersrefri3[NUMOTHERSTATERESULTS],others[NUMOTHERSTATERESULTS],othersrefnew[NUMOTHERSTATERESULTS];
  FTYPE Mdotri3;
  FTYPE Dqp,Dqm,Dqc,dqMdot;
  FTYPE myMdot;
  FTYPE dq[NPR];
  FTYPE dqucon[NDIM];
  FTYPE dqBd3,myBd3;
  FTYPE dqgdetpl[NPR];
  FTYPE dqlogdensity[NPR];
  FTYPE uconreftouse[NDIM];
  FTYPE xfrac,yfrac;
  FTYPE ftemp;
  FTYPE signdq;
  FTYPE xdqfrac,ydqfrac;
  FTYPE ftemp2,linearvalue,expvalue;
  FTYPE D0,uconrefnew[NDIM];
  FTYPE mydqlog;
  FTYPE igdetnosing;
  struct of_geom geomdontuse[NPR];
  struct of_geom *ptrgeom[NPR];
  struct of_geom rgeomdontuse[NPR];
  struct of_geom *ptrrgeom[NPR];
  struct of_geom ri2geomdontuse[NPR];
  struct of_geom *ptrri2geom[NPR];
  struct of_geom ri3geomdontuse[NPR];
  struct of_geom *ptrri3geom[NPR];

  // assign memory
  PALLLOOP(pl){
    ptrgeom[pl]=&(geomdontuse[pl]);
    ptrrgeom[pl]=&(rgeomdontuse[pl]);
    ptrri2geom[pl]=&(ri2geomdontuse[pl]);
    ptrri3geom[pl]=&(ri3geomdontuse[pl]);
  }


  rj=j;
  rk=k;


  // setup multiple reference points
  if(boundary==X1DN){
    ri = riin;
    ri2=ri+1;
    ri3=ri+2;
    // iterate up (due to use of single loop form)
    iloopstart=LOWERBOUND1;
    iloopend=ri-1;
    iloopstep=1;
    signdq=1.0;
  }
  else if(boundary==X1UP){
    ri = riout;
    ri2=ri-1;
    ri3=ri-2;
    // iterate up (due to use of single loop form)
    iloopstart=ri+1;
    iloopend=UPPERBOUND1;
    iloopstep=1;
    signdq=-1.0;
  }
  else{
    dualfprintf(fail_file,"extrapfunc not setup for this boundary = %d\n",boundary);
    myexit(811751);
  }



  /////////////
  //
  // get coordinates of reference location
  //
  /////////////
  PBOUNDLOOP(pliter,pl){
    get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
    get_geometry(ri2, rj, rk, dirprim[pl], ptrri2geom[pl]);
    get_geometry(ri3, rj, rk, dirprim[pl], ptrri3geom[pl]);

    bl_coord_ijk_2(ri,rj,rk,dirprim[pl],Xr[pl],Vr[pl]); // reference locations for B2/U2
  }



  /////////////
  //
  // obtain reference B_\phi
  //
  /////////////
  Bd3=0.0; SLOOPA(jj) Bd3 += MACP0A1(prim,ri,rj,rk,B1+jj-1)*(ptrrgeom[B1+jj-1]->gcov[GIND(3,jj)]);
  Bd3ri2=0.0; SLOOPA(jj) Bd3ri2 += MACP0A1(prim,ri2,rj,rk,B1+jj-1)*(ptrri2geom[B1+jj-1]->gcov[GIND(3,jj)]);
  Bd3ri3=0.0; SLOOPA(jj) Bd3ri3 += MACP0A1(prim,ri3,rj,rk,B1+jj-1)*(ptrri3geom[B1+jj-1]->gcov[GIND(3,jj)]);
  // B_\phi slope
  Dqp=Bd3ri3-Bd3ri2;
  Dqm=Bd3ri2-Bd3;
  Dqc=Bd3ri3-Bd3;
  dqBd3 = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);



  ///////////////////////////////////////////
  //
  // reference value for \detg rho u^{x1}
  //
  ///////////////////////////////////////////


  if(ispstag==0){
    ucon_calc(MAC(prim,ri,rj,rk),ptrrgeom[U1],uconref,othersref);
    Mdot = (ptrrgeom[U1]->gdet)*MACP0A1(prim,ri,rj,rk,RHO)*uconref[1];

    // ri2
    ucon_calc(MAC(prim,ri2,rj,rk),ptrri2geom[U1],uconrefri2,othersrefri2);
    Mdotri2 = (ptrri2geom[U1]->gdet)*MACP0A1(prim,ri2,rj,rk,RHO)*uconrefri2[1];

    // ri3
    ucon_calc(MAC(prim,ri3,rj,rk),ptrri3geom[U1],uconrefri3,othersrefri3);
    Mdotri3 = (ptrri3geom[U1]->gdet)*MACP0A1(prim,ri3,rj,rk,RHO)*uconrefri3[1];

    // Mdot slope
    Dqp=Mdotri3-Mdotri2;
    Dqm=Mdotri2-Mdot;
    Dqc=Mdotri3-Mdot;
    dqMdot = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }



  //////////////////////
  //
  // prim[pl] slope
  //
  //////////////////////
  PBOUNDLOOP(pliter,pl){
    // determine MINM slope for extrapolation
    Dqp=MACP0A1(prim,ri3,rj,rk,pl)-MACP0A1(prim,ri2,rj,rk,pl);
    Dqm=MACP0A1(prim,ri2,rj,rk,pl)-MACP0A1(prim,ri,rj,rk,pl);
    Dqc=MACP0A1(prim,ri3,rj,rk,pl)-MACP0A1(prim,ri,rj,rk,pl);
    dq[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }


  //////////////////////////
  //
  // log(prim[pl]) slope
  //
  //////////////////////////
  PBOUNDLOOP(pliter,pl){
    // determine MINM slope for extrapolation
    Dqp=log(SMALL+fabs(MACP0A1(prim,ri3,rj,rk,pl)))-log(SMALL+fabs(MACP0A1(prim,ri2,rj,rk,pl)));
    Dqm=log(SMALL+fabs(MACP0A1(prim,ri2,rj,rk,pl)))-log(SMALL+fabs(MACP0A1(prim,ri,rj,rk,pl)));
    Dqc=log(SMALL+fabs(MACP0A1(prim,ri3,rj,rk,pl)))-log(SMALL+fabs(MACP0A1(prim,ri,rj,rk,pl)));
    dqlogdensity[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }

  ////////////
  //
  // gdet * prim[pl] slope
  //
  ////////////
  PBOUNDLOOP(pliter,pl){
    // determine MINM slope for extrapolation
    Dqp=(ptrri3geom[pl]->gdet)*MACP0A1(prim,ri3,rj,rk,pl)-(ptrri2geom[pl]->gdet)*MACP0A1(prim,ri2,rj,rk,pl);
    Dqm=(ptrri2geom[pl]->gdet)*MACP0A1(prim,ri2,rj,rk,pl)-(ptrrgeom[pl]->gdet)*MACP0A1(prim,ri,rj,rk,pl);
    Dqc=(ptrri3geom[pl]->gdet)*MACP0A1(prim,ri3,rj,rk,pl)-(ptrrgeom[pl]->gdet)*MACP0A1(prim,ri,rj,rk,pl);
    dqgdetpl[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }


  //////////////////////
  //
  // u^i slope
  //
  //////////////////////
  if(ispstag==0){
    DLOOPA(jj){
      // determine MINM slope for extrapolation
      Dqp=uconrefri3[jj]-uconrefri2[jj];
      Dqm=uconrefri2[jj]-uconref[jj];
      Dqc=uconrefri3[jj]-uconref[jj];
      dqucon[jj] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
    }

    // get fraction of farther away reference value to use
    // Fraction of normal reference value is determined by how similar u^t is for (e.g. for inner boundary) ri and ri+1
    xfrac = fabs(uconref[TT]-uconrefri2[TT])/uconref[TT];
#define XFRAC1 (0.1) // 0.1 was chosen from experience with BH+torus and what the u^t value was doing at the pole near the horizon -- JCM found once d(u^t)>0.1 then started to grow unbound.  Certainly related to resolution so assume resolving of horizon is generally similar to the 64^2 case with Rout=40 and exp grid
#define YFRAC1 (0.0)
#define XFRAC2 (0.2)
#define YFRAC2 (1.0)
#define YFRAC12(frac) (YFRAC1 + (YFRAC2-YFRAC1)/(XFRAC2-XFRAC1)*(frac-XFRAC1))

    if(xfrac<XFRAC1) yfrac=YFRAC1;
    else if(xfrac<XFRAC2) yfrac=YFRAC12(xfrac);
    else yfrac=YFRAC2;


#define XDQFRAC1 (0.1)
#define YDQFRAC1 (0.0)
#define XDQFRAC2 (0.2)
#define YDQFRAC2 (1.0)
#define YDQFRAC12(frac) (YDQFRAC1 + (YDQFRAC2-YDQFRAC1)/(XDQFRAC2-XDQFRAC1)*(frac-XDQFRAC1))

    // limit slope in extrapolation if excessive slope
    xdqfrac = fabs(uconref[TT]-uconrefri2[TT])/uconref[TT];
    if(xdqfrac<XDQFRAC1) ydqfrac=YDQFRAC1;
    else if(xdqfrac<XDQFRAC2) ydqfrac= YDQFRAC12(xdqfrac);
    else ydqfrac=YDQFRAC2;
    
    // modify dqucon based upon ydqfrac
    // Uses original MINM slope unless slope is too large, then start to decrease to no slope
    // only apply to velocity related slopes
    DLOOPA(jj){
      dqucon[jj] = dqucon[jj]*(1.0-ydqfrac);
      dq[UU+jj] = dq[UU+jj]*(1.0-ydqfrac);
      dqlogdensity[UU+jj] = dqlogdensity[UU+jj]*(1.0-ydqfrac);
      dqMdot = dqMdot*(1.0-ydqfrac);
      dqgdetpl[UU+jj] = dqgdetpl[UU+jj]*(1.0-ydqfrac);
    }



  }// end if ispstag==0
  else{
    yfrac=ydqfrac=0.0;
  }







  /////////////////////
  //
  // LOOP over i and pl
  //
  /////////////////////

  LOOPBOUNDGENMORE(i,iloopstart,iloopend,iloopstep){
  


    //////////////////
    //
    // default copy
    //
    //////////////////
    PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,rk,pl);

    //////////////////
    //
    // get coordinates of local location for all quantities
    //
    //////////////////
    PBOUNDLOOP(pliter,pl){
      get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
      bl_coord_ijk_2(i,j,k,dirprim[pl],X[pl],V[pl]); // reference locations for B2/U2
    }





    ///////////////////////////////
    //
    // NON-magnetic field parts (\rho,u,u^i, etc.)
    //
    ///////////////////////////////
    if(ispstag==0){


    ///////////////////////////////
    //
    // Densities and variables at higher than B3
    //
    ///////////////////////////////

    /// did assume densities roughly constant (good for disk, not as good for polar region of mostly radial Bondi-like flow)
    // now focus on accuracy for polar regions to maintain stability (avoid dissipation terms doing funny things)

    // using below Bondi-like scaling with Mdot scaling for ucon leads to space-like answers  Near BH densities should become roughly constant, and this doesn't lead to problems, so use instead
    //	  pl=RHO; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*pow(V[pl][1]/Vr[pl][1],-3.0/2.0);
    //	  pl=UU;  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*pow(V[pl][1]/Vr[pl][1],-5.0/2.0);

#define POWEREXTRAFRAC (0.1)
#define POWERRATIO (1.0+POWEREXTRAFRAC) // ratio of densities not allow to go bigger than this when using dqlogdensity


#if(1)
      for(pl=0;pl<NPRBOUND;pl++){
	if(pl<=UU || pl>B3){

	  // only use linear if exponentiation causes growth of value, not decreasion
	  ftemp = exp(-signdq*dqlogdensity[pl]) - POWERRATIO;

	  // limit interpolation strength to fixed log-slope
	  if(ftemp>=0.0) mydqlog=log(POWERRATIO)/(-signdq);
	  else mydqlog = dqlogdensity[pl];

	  // log extrap (very speculative and can cause problems if used alone when (say) density is super low on ri+1 and relatively high on ri, then i will be super huge
	  // should use this when values that go into slope are much different, or equally when dqlogdensity is large
	  MACP0A1(prim,i,j,k,pl) = exp(log(SMALL+fabs(MACP0A1(prim,ri,rj,rk,pl))) + mydqlog*(i-ri));

	  // DEBUG:
	  //	  dualfprintf(fail_file,"i=%d j=%d pl=%d ftemp=%21.15g linearvalue=%21.15g expvalue=%21.15g final=%21.15g\n",i,j,pl,ftemp,linearvalue,expvalue,MACP0A1(prim,i,j,k,pl));

	}
      }

#else

      // override using OLD WAY for densities
#if(0)
    // linear extrap (leads to negative values, which can cause problems with flux inversions)
    pl=RHO; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) + dq[pl]*(i-ri);
    pl=UU;  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) + dq[pl]*(i-ri);
#else
    // log extrap
    pl=RHO; MACP0A1(prim,i,j,k,pl) = exp(log(SMALL+fabs(MACP0A1(prim,ri,rj,rk,pl))) + dqlogdensity[pl]*(i-ri));
    pl=UU;  MACP0A1(prim,i,j,k,pl) = exp(log(SMALL+fabs(MACP0A1(prim,ri,rj,rk,pl))) + dqlogdensity[pl]*(i-ri));
#endif

#endif



	  
      ///////////////////////////////
      //
      // Velocity
      //
      ///////////////////////////////

#if(0)
      // appears too speculative and use of v^3 causes ucon2pr to lead to no solution within boundary conditions -- suggestive of some inconsistency
      //  (primitive can be anything)

      // combined with \rho_0=constant, this implies \detg \rho_0 u^x1 == constant as required for Mdot=constant for radial flow
      //	  pl=U1; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*fabs((ptrrgeom[pl]->gdet)/(ptrgeom[pl]->gdet));
      myMdot = Mdot + dqMdot*(i-ri);
      ucon[1] = myMdot/((ptrgeom[U1]->gdet)*MACP0A1(prim,i,j,k,RHO));
      // assume ucon[2] and ucon[3] are just copied directly
      // assume mostly radial flow, but visibly appears like B2 so treat similarly so no large kink in behavior
      ucon[2] = uconref[2] + dqucon[2]*(i-ri);
      // assume v^\phi \sim \Omega doesn't change much
      ucon[3] = uconref[3] + dqucon[3]*(i-ri);

      // using non-computed u^t in ucon2pr means slightly inconsistent with desired u^t[u^i], but assume approximately correct and is sufficient for extrapolation to get primitive
      //	  ucon[TT] = uconref[TT];

      // fill so can use standard functions
      primtemp[U1] = ucon[1];
      primtemp[U2] = ucon[2];
      primtemp[U3] = ucon[3];
      // get real ucon with real u^t
      ucon_calc_4vel_bothut(primtemp, ptrgeom[U1], ucon, ucon2, others);
      // now get primitive velocity (choosing 4-velocity branch closest to reference velocity)
      // using the below causes flipping of choice in accretion flow region of boundary conditions due to being too far away from original uconref[TT]
      //	  if(fabs(ucon[TT]-uconref[TT])<fabs(ucon2[TT]-uconref[TT])) ucontouse=ucon;
      //	  else ucontouse=ucon2;
      ucontouse=ucon;

      /////////////////
      //
      // iterate to get closer to v^\phi\sim constant (safer than setting v^i directly and using vcon2pr())
      //
      /////////////////

      for(itervphi=0;itervphi<=NUMITERVPHI;itervphi++){
	// assume v^\phi \sim \Omega doesn't change much
	primtemp[U3] = uconref[3]/uconref[TT]*ucontouse[TT];
	// get real ucon with real u^t
	ucon_calc_4vel_bothut(primtemp, ptrgeom[U1], ucon, ucon2, others);
	// now get primitive velocity (choosing 4-velocity branch closest to reference velocity)
	//	    if(fabs(ucon[TT]-uconref[TT])<fabs(ucon2[TT]-uconref[TT])) ucontouse=ucon;
	//	    else ucontouse=ucon2;
	ucontouse=ucon;
      }

      /////////////////////////
      // get final primitive from ucon
      ucon2pr(WHICHVEL,ucontouse,ptrgeom[U1],MAC(prim,i,j,k)); // fills only U1,U2,U3


#elif(0)

      // GODMARK: Was causing problems for Rebecca


      // BH-pole corner angular sucking drive u^t>>0 and then failures occur
      // avoid exaggerating via unwanted radial sucking by using reference value with

      // whether to conserve (at least) D=\rho_0 u^t when modifying \gamma
      // see notes in fixup.c for limit_gamma() regarding why this is not a good idea
#define DO_CONSERVE_D_INBOUNDS 0 // CHANGINGMARK



#if(DO_CONSERVE_D_INBOUNDS)
      // GODMARK: here we actually modify active values of quantities (should we preserve D=\rho_0 u^t ? -- i.e. change \rho_0?)
      // This correction comes after any other velocity modifications
      // of all conservation laws, at least conserve particle number
      D0 = MACP0A1(prim,ri,rj,rk,RHO)*uconref[TT];
#endif

      for(pl=U1;pl<=U3;pl++){

	// switch away from using near-BC value as reference if going bad since don't want to exaggerate it
	ftemp = MACP0A1(prim,ri,rj,rk,pl)*(1.0-yfrac) + MACP0A1(prim,ri2,rj,rk,pl)*yfrac;
	// interpolate
	MACP0A1(prim,i,j,k,pl) = (ftemp + dq[pl]*(i-ri));


	//	dualfprintf(fail_file,"i=%d j=%d k=%d pl=%d ftemp=%21.15g yfrac=%21.15g dq[pl]=%21.15g i=%d ri=%d ri2=%d ri3=%d rj=%d rk=%d :: prim=%21.15g pri2=%21.15g pri3=%21.15g\n",i,j,k,pl,ftemp,yfrac,dq[pl],i,ri,ri2,ri3,rj,rk,MACP0A1(prim,i,j,k,pl),MACP0A1(prim,ri2,rj,rk,pl),MACP0A1(prim,ri3,rj,rk,pl)); // CHANGINGMARK



	if(V[pl][RR]<RADIUSTOAVOIDRADIALSUCK){
	  // only do this if near a BH
	  // also use yfrac to reset on-grid value since apparently it's going bad
	  // this also keeps interpolation scheme passing through to boundary and so avoiding unwanted dissipation that can evolve flow away from sanity
	  // linear extrap
	  // GODMARK: Should probably make sure that u^r doesn't change sign?!
	  MACP0A1(prim,ri,rj,rk,pl) = MACP0A1(prim,ri,rj,rk,pl)*(1.0-yfrac) + (MACP0A1(prim,ri2,rj,rk,pl) + dq[pl]*(ri-(ri2)))*yfrac;
	}

      }// end over velocities

#if(DO_CONSERVE_D_INBOUNDS)
      // get new u^t
      ucon_calc(MAC(prim,ri,rj,rk),ptrrgeom[U1],uconrefnew, othersrefnew);
      // recompute \rho_0 to conserve particle number (generally will increase \rho_0)
      MACP0A1(prim,ri,rj,rk,RHO) = D0/uconrefnew[TT];
#endif




#elif(0)
      // linear extrap
      for(pl=U1;pl<=U3;pl++){

	// switch away from using near-BC value as reference if going bad since don't want to exaggerate it
	ftemp = MACP0A1(prim,ri,rj,rk,pl)*(1.0-yfrac) + MACP0A1(prim,ri3,rj,rk,pl)*yfrac;
	// interpolate
	MACP0A1(prim,i,j,k,pl) = ftemp + dq[pl]*(i-ri);

      }


#elif(1)
      // linear extrap
      for(pl=U1;pl<=U3;pl++){

	ftemp = MACP0A1(prim,ri,rj,rk,pl);
	// interpolate
	MACP0A1(prim,i,j,k,pl) = ftemp + dq[pl]*(i-ri);

      }


#endif	  


    } // end if ispstag==0



    ///////////////////////////////
    //
    // MAGNETIC FIELD
    //
    ///////////////////////////////

    ////////
    //
    // assume \detg B^x1 == constant (well, extrapolated now)
    //
    // GODMARK: Should probably make sure that B^r and B^\phi don't change sign! (well, maybe at least B^r)
    // Should probably preserve signature of B^\phi/B^r = -\Omega_F/c -- suggests interpolating B^\phi/B^r instead of B_\phi -- well, can't divide by B^r
    //
    ////////
    for(pl=B1;pl<=B2;pl++){
      igdetnosing=sign(ptrgeom[pl]->gdet)/(fabs(ptrgeom[pl]->gdet)+SMALL); // avoids 0.0 for any sign of ptrgeom->gdet
      MACP0A1(prim,i,j,k,pl) =  (MACP0A1(prim,ri,rj,rk,pl)*fabs((ptrrgeom[pl]->gdet) + dqgdetpl[pl]*(i-ri))*igdetnosing);
    }

    ///////////////
    //
    // obtian generally correct B^\phi[B^r,B^\theta,B_\phi]
    //
    ///////////////
    if(dirprim[B3]==FACE3){
      // then assume staggered method
      // need to average fields
      //      Bu1=0.25*(MACP0A1(prim,i,j,k,B1)+MACP0A1(prim,ip1,j,k,B1)+MACP0A1(prim,i,j,km1,B1)+MACP0A1(prim,ip1,j,km1,B1));
      //Bu2=0.25*(MACP0A1(prim,i,j,k,B2)+MACP0A1(prim,i,jp1,k,B2)+MACP0A1(prim,i,j,km1,B2)+MACP0A1(prim,i,jp1,km1,B2));
      // assume average just results in the same value since can't average non-set values
      Bu1=MACP0A1(prim,i,j,k,B1);
      Bu2=MACP0A1(prim,i,j,k,B2);
    }
    else{
      // then assume all fields at CENT
      Bu1=MACP0A1(prim,i,j,k,B1);
      Bu2=MACP0A1(prim,i,j,k,B2);
    }
    gcon03=ptrgeom[B3]->gcon[GIND(0,3)];
    gcon13=ptrgeom[B3]->gcon[GIND(1,3)];
    gcon23=ptrgeom[B3]->gcon[GIND(2,3)];
    gcon33=ptrgeom[B3]->gcon[GIND(3,3)];

    gcov01=ptrgeom[B3]->gcov[GIND(0,1)];
    gcov02=ptrgeom[B3]->gcov[GIND(0,2)];
    gcov11=ptrgeom[B3]->gcov[GIND(1,1)];
    gcov12=gcov21=ptrgeom[B3]->gcov[GIND(1,2)];
    gcov22=ptrgeom[B3]->gcov[GIND(2,2)];
    gcov03=ptrgeom[B3]->gcov[GIND(0,3)];
    gcov13=ptrgeom[B3]->gcov[GIND(1,3)];
    gcov23=ptrgeom[B3]->gcov[GIND(2,3)];
	  
    // Bd3fromBu.nb (just moved signs)
    myBd3=Bd3 + dqBd3*(i-ri);
    ftemp=(1.0 - gcon03*gcov03 - gcon13*gcov13 - gcon23*gcov23);
    igdetnosing=sign(ftemp)/(fabs(ftemp)+SMALL);
    pl=B3; MACP0A1(prim,i,j,k,pl) = (myBd3*gcon33 + Bu1*gcon03*gcov01 + Bu2*gcon03*gcov02 + Bu1*gcon13*gcov11 + Bu2*gcon13*gcov12 + Bu1*gcon23*gcov21 + Bu2*gcon23*gcov22)*igdetnosing;

    // old B_\phi imprecise copy (need to avoid singularity anywways)
    //	  pl=B3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*fabs((ptrrgeom[pl]->gcov[GIND(3,3)])/(ptrgeom[pl]->gcov[GIND(3,3)]));






  } // end LOOP over i and pl




  return(0);



}




#define UTHETAPOLEDEATH ((PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP==VARTOINTERP)?(1):(0)) //only interpolate u^\theta instead of u^2 here when doing the same in interp

#if(0 == UTHETAPOLEDEATH)
#   define MACP0A1mod(prim,ri,rj,rk,pl) MACP0A1(prim,ri,rj,rk,pl)
#else
//average u^\theta=u^2*dxdxp22 as opposed to u^2
static FTYPE MACP0A1mod(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int ri, int rj, int rk, int pl);

static FTYPE MACP0A1mod(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int ri, int rj, int rk, int pl)
{
    FTYPE dxdxp[NDIM][NDIM];
    int dir;

#if(PRODUCTION==0)
    if ( dir < U2 || dir > B3 ) {
	dualfprintf( fail_file, "dir cannot be < U2 or > B3 in MACP0A1mod()\n" );
	exit(3232);
    }
#endif
    
    dir = 1 + (pl-U1)%3; //direction that corresponds to pl (e.g., 2 for U2 or B2); assumes contiguous ordering of pl: <..>,U1,U2,U3,B1,B2,B3,<..>
    
    dxdxprim_ijk(ri, rj, rk, CENT, dxdxp);
    
    return( MACP0A1(prim,ri,rj,rk,pl) * dxdxp[dir][dir] );   
}
#endif

#define POLEOFFSET (POLESMOOTH)

int polesmooth(int whichx2,
	      int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
	      int *inboundloop,
	      int *outboundloop,
	      int *innormalloop,
	      int *outnormalloop,
	      int (*inoutlohi)[NUMUPDOWN][NDIM],
	      int riin, int riout, int rjin, int rjout, int rkin, int rkout,
	      int *dosetbc,
	      int enerregion,
	      int *localenerpos)
{
  
  int i, j, k, pl, j0, dj;
  int ri, rj, rk;
  FTYPE ucon[NDIM], uconspc[NDIM], uconxyz[NDIM], uconspcavg[NDIM], uconxyzavg[NDIM];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom = &geomdontuse;
  FTYPE X[NDIM], V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE *pr, newpr[NPR], r, th, ph;
  FTYPE rho, ug;

  if( ncpux3 != 1 ){
    dualfprintf(fail_file, "polesmooth() is currently designed to work only for 1 cpu in the phi-direction\n" );
    myexit(2323);
  }
  
  if(ispstag) return(0);
  
  /////////////////////
  //
  // Loop over i
  if (whichx2==X2DN) {
    j0 = -POLEOFFSET; //starting from this j
    rj = POLEOFFSET;  //until this j
    dj = 1;
    if(mycpupos[2]!=0) return(0); //not on pole cpu
  }
  else{
    j0 = N2-1+POLEOFFSET;  //starting from this j
    rj = N2-1-POLEOFFSET; //until this j
    dj = -1;
    if(mycpupos[2]!=ncpux2-1) return(0); //not on pole cpu
  }
  for (i=-N1BND; i<N1+N1BND; i++) {
    //zero out velocities
    DLOOPA(pl) {
      uconspcavg[pl] = 0;
      uconxyzavg[pl] = 0;
    }
    rho = 0;
    ug = 0;
    ri = i;
    for (k=0; k<N3; k++) {
      rk = k;
      bl_coord_ijk_2(ri,rj,rk,CENT,X,V);
      r = V[1];
      th = V[2];
      ph = V[3];
      get_geometry(ri, rj, rk, CENT, ptrgeom);
      dxdxprim_ijk(ri, rj, rk, CENT, dxdxp);

      pr = MAC(prim,ri,rj,rk);
      //obtain coordinate 4-velocity
      pr2ucon(WHICHVEL, pr, ptrgeom, ucon);
      //4-vel in spherical polar
      uconspc[TT] = 0;
      uconspc[RR] = dxdxp[RR][RR]*ucon[RR];
      uconspc[TH] = dxdxp[TH][RR]*ucon[RR] + dxdxp[TH][TH]*ucon[TH];
      uconspc[PH] = dxdxp[PH][PH]*ucon[PH];
      //4-vel in Cartesian, x, y, z (see tiltedAphi.m)
      uconxyz[0]  = 0;
      uconxyz[1]  = -(r*sin(ph)*sin(th)*uconspc[PH]) + cos(ph)*sin(th)*uconspc[RR] + r*cos(ph)*cos(th)*uconspc[TH];
      uconxyz[2]  = r*cos(ph)*sin(th)*uconspc[PH] + sin(ph)*sin(th)*uconspc[RR] + r*cos(th)*sin(ph)*uconspc[TH];
      uconxyz[3]  = cos(th)*uconspc[RR] - r*sin(th)*uconspc[TH];
      DLOOPA(pl) {
	uconspcavg[pl] += uconspc[pl];
	uconxyzavg[pl] += uconxyz[pl];
      }
      rho += pr[RHO];
      ug += pr[UU];
    }
    //divide by numer of terms (assume single CPU)
    SLOOPA(pl) {
      uconspcavg[pl] /= (FTYPE)N3;
      uconxyzavg[pl] /= (FTYPE)N3;
    }
    rho /= (FTYPE)N3;
    ug /= (FTYPE)N3;
    
    //now populate the affected cells with averaged rho, ug, and velocities
    for (k=-N3BND; k<N3+N3BND; k++) {
      for (j=j0; j != rj; j+=dj) {
	//NOTE below uses j (as opposed to rj)
	bl_coord_ijk_2(i,j,k,CENT,X,V);
	r = V[1];
	th = V[2];
	ph = V[3];
	get_geometry(i, j, k, CENT, ptrgeom);
	dxdxprim_ijk(i, j, k, CENT, dxdxp);
	
	//obtain the spherical polar velocities from Cartesian (x,y,z) velocities
	uconspc[RR] = cos(ph)*sin(th)*uconxyzavg[1] + sin(ph)*sin(th)*uconxyzavg[2] + cos(th)*uconxyzavg[3];
	uconspc[TH] = pow(r,-1)*(cos(ph)*cos(th)*uconxyzavg[1] + cos(th)*sin(ph)*uconxyzavg[2] - sin(th)*uconxyzavg[3]);
	uconspc[PH] = cos(ph)*1./sin(th)*pow(r,-1)*(-(tan(ph)*uconxyzavg[1]) + uconxyzavg[2]);
	
	//add in average rotational velocity to account for any rotation around the polar axis, if any
	uconspc[PH] += uconspcavg[PH];
	
	//convert back to code coordinates
	//first, get coordinate 4-velocity
	ucon[RR] = uconspc[RR]/dxdxp[RR][RR];
	ucon[TH] = (uconspc[TH] - dxdxp[TH][RR]*ucon[RR]) / dxdxp[TH][TH];
	ucon[PH] = uconspc[PH]/dxdxp[PH][PH];
	
	//convert coordinate 4-velocity into into internal coordinate velocity (usually VELREL4 velocity)
	ucon2pr(WHICHVEL,ucon,ptrgeom,newpr);
	//put the updated velocity back into pr
	pr = MAC(prim,i,j,k);
	SLOOPA(pl) pr[UU+pl] = newpr[UU+pl];
	pr[RHO] = rho;
	pr[UU] = ug;
      }
    }      
  }
  
  return(0);
  
}

// interpolate across pole to remove numerical errors there
// Note that this function is done over all zones
int poledeath(int whichx2,
	      int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
	      int *inboundloop,
	      int *outboundloop,
	      int *innormalloop,
	      int *outnormalloop,
	      int (*inoutlohi)[NUMUPDOWN][NDIM],
	      int riin, int riout, int rjin, int rjout, int rkin, int rkout,
	      int *dosetbc,
	      int enerregion,
	      int *localenerpos)
{




  // when poledeath() called, already in parallel region
  //#pragma omp parallel  // assume don't require EOS
  //  {

  int jj;
  int i,j,k;
  int pl,pliter;
  int rim1,rip1;
  int ri,rj,rk;
  int rj0;
  int rjstag;
  int rjtest,rjstag0,deathstagjs0,deathstagje0,rjstagtest,deathstagjstest,deathstagjetest;
  int poleloc;
  int deathjs0,deathje0;
  int deathjs,deathje;
  int deathstagjs,deathstagje;
  int gammadeathjs,gammadeathje;
  FTYPE X[NPR][NDIM];
  FTYPE V[NPR][NDIM];
  FTYPE Xr[NPR][NDIM];
  FTYPE X0[NDIM];
  int lowrho,lowuu,highu1;
  int deathjstest,deathjetest;
  FTYPE gamma,gammarj0,gammarjtest;
  FTYPE qsq,qsqrj0,qsqrjtest;
  FTYPE ftemp;
  FTYPE ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE gammavaluelimit;
  int doavginradius[NPR];
  int pl2;
  struct of_geom geomdontuse[NPR];
  struct of_geom *ptrgeom[NPR];
  struct of_geom rgeomdontuse[NPR];
  struct of_geom *ptrrgeom[NPR];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE prdiag[NPR],pr[NPR];
  int jstep;



  // OPENMPMARK: Can't do this check because if reduce to 1 thread (even in OpenMP) then omp_in_parallel() is false!
  //#if(USEOPENMP)
  //  // ensure within parallel region
  //  if(!omp_in_parallel()){
  //    dualfprintf(fail_file,"poledeath() called outside parallel region\n");
  //    myexit(84968346);
  //  }
  //#endif


  // assign memory
  PALLLOOP(pl){
    ptrgeom[pl]=&(geomdontuse[pl]);
    ptrrgeom[pl]=&(rgeomdontuse[pl]);
  }


  // note that doesn't matter the order of the j-loop since always using reference value (so for loop doesn't need change in <= to >=)
  if(whichx2==X2DN){

    jstep=-1; // direction of j loop, so starts with active cells in case modify boundary cells as dependent upon active cells.

    rj0 = POLEDEATH;
    rjtest = rj0+DEATHEXPANDAMOUNT; // used to ensure near pole the density doesn't drop suddenly
    poleloc = 0;
    // for POLEDEATH==2, then deathjs,je=-2,-1,0,1 as required for CENT quantities rj=2
    deathjs0 = 0-POLEDEATH;
    deathje0 = 0+POLEDEATH-1;

    deathjstest = deathjs0-DEATHEXPANDAMOUNT;
    deathjetest = deathje0+DEATHEXPANDAMOUNT;
    if(deathjstest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathjstest=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(deathjetest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathjetest=inoutlohi[POINTDOWN][POINTDOWN][2];

    // assume for POLEDEATH==1 that B2 set correctly as 0 on pole and only do something if POLEDEATH>1
    // if POLEDEATH==2 then B2 set at  -1,0,1 and will correctly set B2 to 0 at pole rj=2
    rjstag0 = rj0;
    deathstagjs0 = 0-POLEDEATH+1;
    deathstagje0 = 0+POLEDEATH-1;

    rjstagtest = rjtest;
    deathstagjstest = deathstagjs0-DEATHEXPANDAMOUNT;
    deathstagjetest = deathstagje0+DEATHEXPANDAMOUNT;
    if(deathstagjstest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathstagjstest=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(deathstagjetest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathstagjetest=inoutlohi[POINTDOWN][POINTDOWN][2];

    // assumes velocity is always CENT
    gammadeathjs=0-POLEGAMMADEATH;
    gammadeathje=0+POLEGAMMADEATH-1;
    if(gammadeathjs<inoutlohi[POINTDOWN][POINTDOWN][2]) gammadeathjs=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(gammadeathje<inoutlohi[POINTDOWN][POINTDOWN][2]) gammadeathje=inoutlohi[POINTDOWN][POINTDOWN][2];

    // NO, don't do below.  Since poledeath called after MPI, need to set ghost cells as consistent with how active cells would have been set by other part of grid or other CPUs
    //    if(special3dspc){
    //      // then assume poledeath called *after* MPI (so have full and correct information across pole), so only should modify active cells and not boundary cells
    //      deathjs0 = 0;
    //      deathjstest = 0;
    //      deathstagjs0 = 0;
    //      deathstagjstest = 0;
    //      gammadeathjs=0;
    //    }

  }
  else if(whichx2==X2UP){
    rj0=N2-1-POLEDEATH;
    rjtest = rj0-DEATHEXPANDAMOUNT;
    poleloc=N2;
    // if POLEDEATH==2 then CENTs set at N2-2,N2-1,N2,N2+1 rj=N2-3
    deathjs0 = N2-1+1-POLEDEATH;
    deathje0 = N2-1+POLEDEATH;

    deathjstest = deathjs0-DEATHEXPANDAMOUNT;
    deathjetest = deathje0+DEATHEXPANDAMOUNT;
    if(deathjstest>inoutlohi[POINTUP][POINTUP][2]) deathjstest=inoutlohi[POINTUP][POINTUP][2];
    if(deathjetest>inoutlohi[POINTUP][POINTUP][2]) deathjetest=inoutlohi[POINTUP][POINTUP][2];

    // if POLEDEATH==2, then B2 is set at N2-1,N2,N2+1 rj=N2-3
    if(dirprim[B2]==FACE2) rjstag0=N2-POLEDEATH;
    else if(dirprim[B2]==CENT) rjstag0=rj0;
    deathstagjs0 = N2+1-POLEDEATH;
    deathstagje0 = N2-1+POLEDEATH;

    rjstagtest = rjtest;
    deathstagjstest = deathstagjs0-DEATHEXPANDAMOUNT;
    deathstagjetest = deathstagje0+DEATHEXPANDAMOUNT;
    if(deathstagjstest>inoutlohi[POINTUP][POINTUP][2]) deathstagjstest=inoutlohi[POINTUP][POINTUP][2];
    if(deathstagjetest>inoutlohi[POINTUP][POINTUP][2]) deathstagjetest=inoutlohi[POINTUP][POINTUP][2];


    // assumes velocity is always CENT .  If POLEDEATH==2, N2-2,N2-1,N2,N2+1
    gammadeathjs = N2-1+1-POLEGAMMADEATH;
    gammadeathje = N2-1+POLEGAMMADEATH;
    if(gammadeathjs>inoutlohi[POINTUP][POINTUP][2]) gammadeathjs=inoutlohi[POINTUP][POINTUP][2];
    if(gammadeathje>inoutlohi[POINTUP][POINTUP][2]) gammadeathje=inoutlohi[POINTUP][POINTUP][2];


    //    if(special3dspc){
    //      // then assume poledeath called *after* MPI (so have full and correct information across pole), so only should modify active cells and not boundary cells
    //      deathje0 = N2-1;
    //      deathjetest = N2-1;
    //      deathstagje0 = N2-1;
    //      deathstagjetest = N2-1;
    //      gammadeathje=N2-1;
    //    }

  }
  



  /////////////////////
  //
  // Loop over i,k


  if(POLEDEATH){

    ///////    LOOPX2dir{
    OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMPBCLOOPBLOCK{
      OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


      // set reference positions in i,k
      ri = i;
      rk = k;


      //////////////////
      //
      // set radial positions to average over
      //
      //////////////////
      if(ri==inoutlohi[POINTDOWN][POINTDOWN][1]){
	// shift up, but not including ri
	rim1=ri+1;
	rip1=ri+2;
      }
      else if(ri==inoutlohi[POINTUP][POINTUP][1]){
	// shift down, but not including ri
	rim1=ri-2;
	rip1=ri-1;
      }
      else{
	// around ri
	rim1=ri-1;
	rip1=ri+1;
      }
    




      /////////////
      // first check for significant density drops near the axis or Lorentz factor spikes up near axis
      lowrho=lowuu=0;
      highu1=0;


      if(0){ // disabled for now

	get_geometry(ri, rj0, rk, dirprim[U1], ptrrgeom[U1]);
	gamma_calc(MAC(prim,ri,rj0,rk), ptrrgeom[U1],&gammarj0,&qsqrj0);
    
	get_geometry(ri, rjtest, rk, dirprim[U1], ptrrgeom[U1]);
	gamma_calc(MAC(prim,ri,rjtest,rk), ptrrgeom[U1],&gammarjtest,&qsqrjtest);



	for (j = deathjs0; j <= deathje0; j++) { // (no prims modified here, so no need for special j loop or diag_fixup)

	  //POLEDENSITYDROPFACTOR
	  //POLEGAMMAJUMPFACTOR


	  pl=RHO;
	  if(fabs(MACP0A1(prim,i,j,k,pl))< fabs(MACP0A1(prim,ri,rj0,rk,pl))/POLEDENSITYDROPFACTOR || fabs(MACP0A1(prim,i,j,k,pl))< fabs(MACP0A1(prim,ri,rjtest,rk,pl))/POLEDENSITYDROPFACTOR){
	    lowrho=1;
	  }
	  pl=UU;
	  if(fabs(MACP0A1(prim,i,j,k,pl))< fabs(MACP0A1(prim,ri,rj0,rk,pl)/POLEDENSITYDROPFACTOR) || fabs(MACP0A1(prim,i,j,k,pl))< fabs(MACP0A1(prim,ri,rjtest,rk,pl))/POLEDENSITYDROPFACTOR){
	    lowuu=1;
	  }


	  // get Lorentz factor
	  get_geometry(i, j, k, dirprim[U1], ptrrgeom[U1]);
	  gamma_calc(MAC(prim,i,j,k), ptrrgeom[U1],&gamma,&qsq);


	  //      pl=U1;
	  //      if(fabs(MACP0A1(prim,i,j,k,pl))> fabs(MACP0A1(prim,ri,rj0,rk,pl))/POLEJUMPFACTOR || fabs(MACP0A1(prim,i,j,k,pl))> fabs(MACP0A1(prim,ri,rjtest,rk,pl))/POLEJUMPFACTOR ){
	  //	highu1=1;
	  //      }

	  // check variation of Lorentz factor
	  if(gamma> gammarj0*POLEGAMMAJUMPFACTOR || gamma> gammarjtest*POLEGAMMAJUMPFACTOR ){
	    highu1=1;
	  }


	}
    
	if(lowuu || lowrho || highu1){
	  rj=rjtest; // expand copy procedure to use better reference value
	  deathjs=deathjstest;
	  deathje=deathjetest;

	  rjstag=rjstagtest; // expand copy procedure to use better reference value
	  deathstagjs=deathstagjstest;
	  deathstagje=deathstagjetest;
	}
	else{
	  // then use normal reference and range
	  rj=rj0;
	  deathjs=deathjs0;
	  deathje=deathje0;

	  rjstag=rjstag0;
	  deathstagjs=deathstagjs0;
	  deathstagje=deathstagje0;
	}
      }
      else{
	rj=rjtest; // expand copy procedure to use better reference value
	deathjs=deathjstest;
	deathje=deathjetest;
      
	rjstag=rjstagtest; // expand copy procedure to use better reference value
	deathstagjs=deathstagjstest;
	deathstagje=deathstagjetest;
      }




      ////////////////////////////////
      //
      // reference location geometry and coordinates
      //
      ////////////////////////////////

      PBOUNDLOOP(pliter,pl) {
	coord_ijk(ri,rj,rk,dirprim[pl],Xr[pl]); // reference locations for B2/U2
	if(pl==B2 && dirprim[B2]==FACE2){
	  get_geometry(ri, rjstag, rk, dirprim[pl], ptrrgeom[pl]);
	}
	else{
	  get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
	}
      }






      //#define DEATHLOOP2 for (j = MIN(deathstagjs,deathjs); j <= MAX(deathstagje,deathje); j++) {
      //      for (j = MIN(deathstagjs,deathjs); j <= MAX(deathstagje,deathje); j++) {
      // always iterate from active to ghost cell, so can modify active and have ready for ghost
      // didn't redefine js and je, they always go from small to large values, which is why ordering below
#define DEATHLOOP2(whichx2) for (j = (whichx2==X2UP ? MIN(deathstagjs,deathjs) : MAX(deathstagje,deathje)) ; (whichx2==X2UP ? (j <= MAX(deathstagje,deathje)) : (j >= MIN(deathstagjs,deathjs)) )   ; (whichx2==X2UP ? j++ : j--) )

      ////////////////////////////////
      //
      // Loop over points to be modified (prim's are only modifed here)
      //
      ////////////////////////////////

      DEATHLOOP2(whichx2){



	//////////////
	//
	// setup initial pr
	//
	//////////////
	PLOOP(pliter,pl) prdiag[pl]=MACP0A1(prim,i,j,k,pl);
	int madechange=0;




	PBOUNDLOOP(pliter,pl) {
	  // pole location
	  coord_ijk(i,poleloc,k,FACE2,X0); // pole locations for B2/U2

	  // current location
	  bl_coord_ijk_2(i,j,k,dirprim[pl],X[pl],V[pl]);
	  get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

	  if(AVERAGEINRADIUS && (V[pl][RR]>RADIUSTOSTARTAVERAGING)){
	    doavginradius[pl]=1;
	  }
	  else doavginradius[pl]=0;

	}




	if(ispstag==0){
	  // symmetric (if reflecting BC at pole) quantities
	  if(j>=deathjs && j<=deathje){


	    //////////
	    // u1, and u3, average in radius too!
	    // copying this means copying \Omega_F in magnetically-dominated regime beyond LC
	    for(pl=U1;pl<=U3;pl++){
	      if(pl==U2) continue;
		
	      if(doavginradius[pl]) MACP0A1(prim,i,j,k,pl) = THIRD*(MACP0A1(prim,rim1,rj,rk,pl)+MACP0A1(prim,ri,rj,rk,pl)+MACP0A1(prim,rip1,rj,rk,pl));
	      else MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
	      madechange++;
	    }



	    if(EOMTYPE!=EOMFFDE){
	      //////////
	      // for densities
	      // this helps remove drop-outs in density at high b^2/\rho_0 and high b^2/u
	      for(pl=RHO;pl<=UU;pl++){
		
		if(doavginradius[pl]) MACP0A1(prim,i,j,k,pl) = THIRD*(MACP0A1(prim,rim1,rj,rk,pl)+MACP0A1(prim,ri,rj,rk,pl)+MACP0A1(prim,rip1,rj,rk,pl));
		else MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
		madechange++;
	      }


	      if(BCSIGMACONSTATPOLE==1){
		
		pl=RHO;
		FTYPE bsqr,sigmar;
		
		if (bsq_calc(MAC(prim,ri,rj,rk), ptrrgeom[pl], &bsqr) >= 1) FAILSTATEMENT("bounds.tools.c:poledeath()", "bsq_calc()", 1);
		sigmar=bsqr/(2.0*fabs(MACP0A1(prim,ri,rj,rk,pl)));
		
		// ensure constant sigma at pole
		FTYPE bsq,sigma;
		if (bsq_calc(MAC(prim,i,j,k), ptrgeom[pl], &bsq) >= 1) FAILSTATEMENT("bounds.tools.c:poledeath()", "bsq_calc()", 2);
		sigma=bsq/(2.0*fabs(MACP0A1(prim,i,j,k,pl)));
		
		if(sigma>BSQORHOLIMIT*0.5 && sigmar<BSQORHOLIMIT*0.5){
		  MACP0A1(prim,i,j,k,pl) = fabs(MACP0A1(prim,i,j,k,pl)*(sigma/sigmar));
		  madechange++;
		}
		// do nothing different than simple copy in any other cases.  NOTE: If not high sigma, else if near-pole is low sigma, then this feeds in mass crazily
	      }// end if BCSIGMACONSTATPOLE==1
	    }// end if EOMTYPE!=EOMFFDE


	  }// end if correct j
	}// end if ispstag==0



	// symmetric (if reflecting BC at pole) quantities
	if(j>=deathjs && j<=deathje){
	  // B1 left alone
	}



	if(
	   dirprim[B2]==FACE2 && j>=deathstagjs && j<=deathstagje || 
	   dirprim[B2]==CENT && j>=deathjs && j<=deathje
	   ){
	  
	  if(N2==1){
	    /////////////////////////////
	    // B2:
#if(POLEINTERPTYPE==0)
	    // if flow converges toward pole, then this loses information about the velocity and field approaching the pole
	    // anti-symmetric (if reflecting BC at pole) quantities:
	    MACP0A1(prim,i,j,k,B2) = 0.;
	    madechange++;

#elif(POLEINTERPTYPE==1 || POLEINTERPTYPE==2)
	    // anti-symmetric (if reflecting BC at pole):
	    // assume X[2] goes through 0 at the pole and isn't positive definite
	    pl=B2;
	    if(doavginradius[pl]) ftemp=THIRD*(MACP0A1(prim,rim1,rjstag,rk,pl) + MACP0A1(prim,ri,rjstag,rk,pl) + MACP0A1(prim,rip1,rjstag,rk,pl));
	    else  ftemp=MACP0A1(prim,ri,rjstag,rk,pl);
	    MACP0A1(prim,i,j,k,pl) = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
	    madechange++;

#elif(POLEINTERPTYPE==3)
	    // then don't modify B2 -- trying to avoid instabilities related to divb=0 violation.  And seems B2 behaves ok
#endif
	  }
	  else{
	    // then don't modify B2 to ensure divB=0
	  }

	}





	// symmetric quantity:      
	if(j>=deathjs && j<=deathje){

	  ////////////////
	  // B3:
	  pl=B3;

	  if(N3==1){

#if(POLEINTERPTYPE==0 || POLEINTERPTYPE==1)
	    // symmetric:
	    if(doavginradius[pl]) ftemp=THIRD*(MACP0A1(prim,rim1,rj,rk,pl) + MACP0A1(prim,ri,rj,rk,pl) + MACP0A1(prim,rip1,rj,rk,pl));
	    else ftemp=MACP0A1(prim,ri,rj,rk,pl);
	    MACP0A1(prim,i,j,k,B3) = ftemp;
	    madechange++;

#elif(POLEINTERPTYPE==2)
	    // symmetric:
	    // approximate B_\phi copy, which (unlike copying B3) can resolve singular currents on axis
	    if(doavginradius[pl]) ftemp=THIRD*(MACP0A1(prim,rim1,rj,rk,pl) + MACP0A1(prim,ri,rj,rk,pl) + MACP0A1(prim,rip1,rj,rk,pl));
	    else ftemp=MACP0A1(prim,ri,rj,rk,pl);
	    MACP0A1(prim,i,j,k,pl) =  ftemp*fabs((ptrrgeom[pl]->gcov[GIND(3,3)])/(ptrgeom[pl]->gcov[GIND(3,3)]));
	    madechange++;

#elif(POLEINTERPTYPE==3)
	    // symmetric:
	    if(doavginradius[pl]) ftemp=THIRD*(MACP0A1(prim,rim1,rj,rk,pl) + MACP0A1(prim,ri,rj,rk,pl) + MACP0A1(prim,rip1,rj,rk,pl));
	    else ftemp=MACP0A1(prim,ri,rj,rk,pl);
	    MACP0A1(prim,i,j,k,B3) = ftemp;
	    madechange++;
#endif
	  }// end if N3==1
	  else{
	    // then do nothing if in 3D
	  }


	  if(ispstag==0){

	    ///////////////////////////////////
	    //
	    // do rest if any -- assumed at CENT
	    //
	    ///////////////////////////////////
	    for(pl=B3+1;pl<NPRBOUND;pl++){
	      if(doavginradius[pl]) ftemp=THIRD*(MACP0A1(prim,rim1,rj,rk,pl) + MACP0A1(prim,ri,rj,rk,pl) + MACP0A1(prim,rip1,rj,rk,pl));
	      else ftemp=MACP0A1(prim,ri,rj,rk,pl);
	      MACP0A1(prim,i,j,k,pl)=ftemp;
	      madechange++;
	    }
	  }

	}// end over CENT-type quantities



	if(madechange){
	  ///////////
	  //
	  // accounting for on-grid changes
	  //
	  ////////////
	  int modcons=1;
	  FTYPE *ucons;
	  ucons=GLOBALMAC(unewglobal,i,j,k); // GODMARK -- need to pass ucons to bounds if really want !=NOENOFLUX method to work
	  // GODMARK: assume all quantities at the same location since only ispstag==0 modifies relevant primitves, so ptrgeom[pl]->ptrgoem
	  // in general, not sure which pl really exists at this point, so pick first in PBOUNDLOOP loop
	  struct of_geom *fixupptrgeom;
	  PBOUNDLOOP(pliter,pl) {
	    fixupptrgeom=ptrgeom[pl];
	    break;
	  }
	  PLOOP(pliter,pl) pr[pl]=MACP0A1(prim,i,j,k,pl);
	  diag_fixup(modcons,prdiag,pr,ucons,fixupptrgeom,finalstep,COUNTBOUND1);
	  PLOOP(pliter,pl) prdiag[pl]=pr[pl];
	}

      
      }// end loop over j




      ////////////////////////////////
      //
      // Loop over points to be modified (prim's are only modifed here)
      //
      // independent loop over j from above density loop since U2 changes may required RHO and those must already all be changed for all j
      //
      ////////////////////////////////

      DEATHLOOP2(whichx2){



	//////////////
	//
	// setup initial pr
	//
	//////////////
	PLOOP(pliter,pl) prdiag[pl]=MACP0A1(prim,i,j,k,pl);
	int madechange=0;




	PBOUNDLOOP(pliter,pl) {
	  // pole location
	  coord_ijk(i,poleloc,k,FACE2,X0); // pole locations for B2/U2

	  // current location
	  bl_coord_ijk_2(i,j,k,dirprim[pl],X[pl],V[pl]);
	  get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

	  if(AVERAGEINRADIUS && (V[pl][RR]>RADIUSTOSTARTAVERAGING)){
	    doavginradius[pl]=1;
	  }
	  else doavginradius[pl]=0;

	}

      
      
	if(ispstag==0){
	
	  if(j>=deathjs && j<=deathje){
	    //////////////////////////
	    // U2:
	    pl=U2;


	    if(special3dspc){
	      // U2 not necessarily anti-symmetric in this case.
	      // This will be kinda odd if POLEDEATH>1 due to comparing non-local regions.
	      // But, for now, POLEDEATH<=1 is set so make sense.

	      int jother;
	      FTYPE signD;
	      if(j<N2/2){
		jother=-1-j;
	      }
	      else{
		jother=N2-1+(N2-j);
	      }

	      if(j>=jother) signD=+1.0;
	      else signD=-1.0;


	      // see if sucking on pole
	      FTYPE rhovjhere,rhovjother,rhovDiff;
	      rhovjhere =MACP0A1(prim,i,j,k,RHO)*MACP0A1mod(prim,i,j,k,U2);
	      rhovjother=MACP0A1(prim,i,jother,k,RHO)*MACP0A1mod(prim,i,jother,k,U2);

	      rhovDiff=signD*(rhovjhere-rhovjother);
	      // same gdet, so no need to multiply both by same factor for below test
	      // make change to both simultaneously so that when other j is hit, rhovDiff condition is no longer hit and all is consistent as if separate memory field for entire poledeath before final copy-over to primitive memory space.
	      // But do active grid cells first (determined by DEATHLOOPJ) so diag_fixup() occurs on active region for accounting.
	      if(rhovDiff>0.0){
		// then sucking on pole
		// must change active and ghost cells consistently for any number of CPUs
		// so average-out the suck to zero suck by changing the values equally (as weighted by mass)
		FTYPE dU2     =-signD*rhovDiff*0.5/MACP0A1(prim,i,j,k,RHO);
		FTYPE U2jhere = MACP0A1mod(prim,i,j,k,U2);
		FTYPE U2jother= MACP0A1mod(prim,i,jother,k,U2);

		if( (fabs(U2jhere)>fabs(dU2))&&(fabs(U2jother)>fabs(dU2)) ){
		  // only change if in both cases we lower the velocity, not increase.
		  MACP0A1(prim,i,j,k,U2)      += dU2;
		  MACP0A1(prim,i,jother,k,U2) -= dU2;
		  madechange++;
		}
		else{ // then jhere or jother is changed by an absolute magnitude more than its value, which we want to avoid
		  // crush a bit only as much as leaves smaller value changed as much as possible without increasing its magnitude
		  if(fabs(U2jhere)<fabs(U2jother)){
		    // then drop jhere as closest to fixed value (100% change), and reset jother with same value
		    MACP0A1(prim,i,j,k,U2) *= -1.0;
		    MACP0A1(prim,i,jother,k,U2) = MACP0A1(prim,i,j,k,U2);
		    madechange++;
		  }
		  else{
		    // then drops down to value matching other side so D=0 in the end still, so still no sucking.
		    MACP0A1(prim,i,jother,k,U2) *= -1.0;
		    MACP0A1(prim,i,j,k,U2) = MACP0A1(prim,i,jother,k,U2);
		    madechange++;
		  }
		}// end else abs mag of change is larger than 100% for one of the values

		//		  MACP0A1(prim,i,j,k,U2)      =0.0;
		//		  madechange++;
		//		MACP0A1(prim,i,jother,k,U2) += +rhovDiff*0.5/MACP0A1(prim,i,jother,k,RHO); // this taken care of by other j in ghost region

		// Note that for anti-symmetric U2 (i.e. reflective BCs around pole) this is same as crushing regularization leading to U2->0
	      }
	      else{
		// then just enforce linear behavior near pole
		// NOT YET
	      }


	    }
	    else{


#if(POLEINTERPTYPE==0)
	      // if flow converges toward pole, then this loses information about the velocity and field approaching the pole
	      // anti-symmetric (if reflecting BC at pole) quantities:
	      MACP0A1(prim,i,j,k,pl) = 0.;
	      madechange++;

#elif(POLEINTERPTYPE==1 || POLEINTERPTYPE==2)
	      // anti-symmetric (if reflecting BC at pole):
	      // assume X[2] goes through 0 at the pole and isn't positive definite
	      if(doavginradius[pl]) ftemp=THIRD*(MACP0A1mod(prim,rim1,rj,rk,pl) + MACP0A1mod(prim,ri,rj,rk,pl) + MACP0A1mod(prim,rip1,rj,rk,pl));
	      else ftemp=MACP0A1mod(prim,ri,rj,rk,pl);
	      MACP0A1(prim,i,j,k,pl) = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
	      madechange++;

#elif(POLEINTERPTYPE==3)

	      // anti-symmetric (if reflecting BC at pole):

	      // assume X[2] goes through 0 at the pole and isn't positive definite

	      // choose reference value
	      // 	ftemp=THIRD*(MACP0A1mod(prim,rim1,rj,rk,pl) + MACP0A1mod(prim,ri,rj,rk,pl) + MACP0A1mod(prim,rip1,rj,rk,pl));
	      ftemp=MACP0A1mod(prim,ri,rj,rk,pl);

	      if(whichx2==X2DN && ftemp>0.0){
		// then sucking on \theta=0 pole
		// try to minimize sucking on pole by finding minimum U2 around
		for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,MACP0A1mod(prim,ri,jj,rk,pl));
		if(doavginradius[pl]){
		  for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,MACP0A1mod(prim,rip1,jj,rk,pl));
		  for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,MACP0A1mod(prim,rim1,jj,rk,pl));
		}

		ftemp=0.0; // try crushing sucking GODMARK

		// assume ftemp is at reference location
		MACP0A1(prim,i,j,k,pl) = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
		madechange++;
	      }
	      else if(whichx2==X2UP && ftemp<0.0){
		// then sucking on \theta=\pi pole
		for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,MACP0A1mod(prim,ri,jj,rk,pl));
		if(doavginradius[pl]){
		  for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,MACP0A1mod(prim,rip1,jj,rk,pl));
		  for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,MACP0A1mod(prim,rim1,jj,rk,pl));
		}

		ftemp=0.0; // try crushing sucking GODMARK

		// assume ftemp is at reference location (same formula for both poles)
		MACP0A1(prim,i,j,k,pl) = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
		madechange++;
	      }
	      else{
		// otherwise enforce natural regular linear behavior on U2
		if(doavginradius[pl]) ftemp=THIRD*(MACP0A1mod(prim,rim1,rj,rk,pl) + MACP0A1mod(prim,ri,rj,rk,pl) + MACP0A1mod(prim,rip1,rj,rk,pl));
		else ftemp=MACP0A1mod(prim,ri,rj,rk,pl);

		MACP0A1(prim,i,j,k,pl) = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
		madechange++;
	      }

#endif // endif POLEINTERPTYPE==3
	    }// end if special3dspc==0


#if( UTHETAPOLEDEATH )
	    //if interpolated u^\theta, now convert back to u^2
	    dxdxprim_ijk(i, j, k, CENT, dxdxp);
	    MACP0A1(prim,i,j,k,pl) /= dxdxp[1 + (pl-U1)%3][1 + (pl-U1)%3];
	    madechange++;
#endif


	  }// end if correct j range
	}// end if ispstag==0



	if(madechange){ // only need accounting if actually made a change
	  ///////////
	  //
	  // accounting for on-grid changes
	  //
	  ////////////
	  int modcons=1;
	  FTYPE *ucons;
	  ucons=GLOBALMAC(unewglobal,i,j,k); // GODMARK -- need to pass ucons to bounds if really want !=NOENOFLUX method to work
	  // GODMARK: assume all quantities at the same location since only ispstag==0 modifies relevant primitves, so ptrgeom[pl]->ptrgoem
	  // in general, not sure which pl really exists at this point, so pick first in PBOUNDLOOP loop
	  struct of_geom *fixupptrgeom;
	  PBOUNDLOOP(pliter,pl) {
	    fixupptrgeom=ptrgeom[pl];
	    break;
	  }

	  PLOOP(pliter,pl) pr[pl]=MACP0A1(prim,i,j,k,pl);
	  diag_fixup(modcons,prdiag,pr,ucons,fixupptrgeom,finalstep,COUNTBOUND1);
	  PLOOP(pliter,pl) prdiag[pl]=pr[pl];
	}


      }// end loop over j
    }// end loop 13
  }// end if POLEDEATH









  if(POLEGAMMADEATH){ // should come after POLEDEATH to be useful
    if(ispstag==0){
      // fixup


      // assume only dealing with velocities at same location (loop assumes CENT)
      pl=U1;


      ///////    LOOPX2dir{
      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
	OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


	for (j = gammadeathjs; j <= gammadeathje; j++) { // currently not multiple-point dependent, so normal j loop is fine


	  //////////////
	  //
	  // setup initial pr (already using pl, so use pl2 for prdiag assignment/setup)
	  //
	  //////////////
	  PLOOP(pliter,pl2) prdiag[pl2]=MACP0A1(prim,i,j,k,pl2);
	  int madechange=0;


	  ///////////
	  //
	  // forever-kept debug
	  //
	  ///////////
	  PBOUNDLOOP(pliter,pl2){
	    if( !isfinite( MACP0A1(prim,i,j,k,pl2))){
	      dualfprintf(fail_file,"BNDNAN: ispstag=%d t=%21.15g nstep=%ld steppart=%d :: i=%d j=%d k=%d pl2=%d prim=%21.15g\n",ispstag,t,nstep,steppart,i,j,k,pl2,MACP0A1(prim,i,j,k,pl2));
	    }
	  }



	  ///////////
	  //
	  // get u^r and limit if necessary
	  //
	  ///////////
	  get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

	  ucon_calc(MAC(prim,i,j,k),ptrgeom[pl],ucon, others);
	  // only limit velocity if outgoing relative to grid (GODMARK: only valid in KS or BL-like coordinates such that u^r>0 means outgoing w.r.t. an observer at infinity)
	  if(ucon[RR]>0.0){

	    bl_coord_ijk_2(i,j,k,dirprim[pl],X[pl],V[pl]);


	    if(V[pl][RR]<=GAMMAPOLEOUTGOINGRADIUS){
	      // flat \gamma limit up to GAMMAPOLEOUTGOINGRADIUS
	      gammavaluelimit = GAMMAPOLEOUTGOING;
	    }
	    else{
	      gammavaluelimit = GAMMAPOLEOUTGOING*pow(V[pl][RR]/GAMMAPOLEOUTGOINGRADIUS,GAMMAPOLEOUTGOINGPOWER);
	    }

	    limit_gamma(gammavaluelimit,MAC(prim,i,j,k),NULL,ptrgeom[pl],-1);
	    madechange++;
	  }
	  else{
	    limit_gamma(GAMMAPOLEINGOING,MAC(prim,i,j,k),NULL,ptrgeom[pl],-1);
	    madechange++;
	  }




	  if(madechange){
	    ///////////
	    //
	    // accounting for on-grid changes
	    //
	    ////////////
	    int modcons=1;
	    FTYPE *ucons;
	    ucons=GLOBALMAC(unewglobal,i,j,k); // GODMARK -- need to pass ucons to bounds if really want !=NOENOFLUX method to work
	    // GODMARK: assume all quantities at the same location since ispstag==0 is assumed in this section, so ptrgeom[pl]->ptrgoem
	    // in general, not sure which pl2 really exists at this point, so pick first in PBOUNDLOOP loop
	    struct of_geom *fixupptrgeom;
	    PBOUNDLOOP(pliter,pl2) {
	      fixupptrgeom=ptrgeom[pl2];
	      break;
	    }
	    
	    PLOOP(pliter,pl2) pr[pl2]=MACP0A1(prim,i,j,k,pl2);
	    diag_fixup(modcons,prdiag,pr,ucons,fixupptrgeom,finalstep,COUNTBOUND2);
	    PLOOP(pliter,pl2) prdiag[pl2]=pr[pl2];

	  }


	}// end over j
      }// end loop 13
    }// end if ispstag==0
  }// end if POLEGAMMADEATH



  return(0);


}




// Toth version (different because emf vs. flux, where emf has extra zone -- otherwise same!
//reset emf's (E_3) at the boundaries in x1-x2 plane to zero
void user1_adjust_fluxcttoth_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{
  int i, j, k;
  int dir;
  int dimen;
  int dirsign;


  if(ADJUSTFLUXCT==0){
    // then no need to set EMF's to zero to maintain divb=0, although can choose to still set EMF's to constant along appropriate directions to maintain stationarity all the way up to including the boundary. Not sure why that would be important to maintain.
    return;
  }
  

  if( DOGRIDSECTIONING==0) {  //don't do anything if sectioning not done
    return;
  }

  
  if(BCtype[X1DN]==FIXEDUSEPANALYTIC || BCtype[X1UP]==FIXEDUSEPANALYTIC){
    // only do if fixing BCs, not (e.g.) outflowing

    //x1-dimension
    dimen = 1;
    DIRSIGNLOOP(dirsign){
      dir = DIRFROMDIMEN( dimen, dirsign );

      if(dir==X1DN && BCtype[X1DN]==FIXEDUSEPANALYTIC || dir==X1UP && BCtype[X1UP]==FIXEDUSEPANALYTIC){ // otherwise don't do

	i = dofluxreg[ACTIVEREGION][dir];

	//if boundary is not on this processor, do not modify emf's
	if( i < -MAXBND ) continue;
	if( i > N1-1+MAXBND) continue;

	//the boundary is on the processor, so reset emf's to zero at the boundary
	if(dir==X1UP && BCtype[dir]==FIXEDUSEPANALYTIC){ // i.e. don't reset EMF2=0 for lower boundary since that corresponds to rotation
	  COMPFULLLOOPP1_23{ 
	    MACP1A0(emf,2,i,j,k) = 0.0;
	  }
	}
	COMPFULLLOOPP1_23{ 
	  MACP1A0(emf,3,i,j,k) = 0.0;
	}
      }
    }
  }


#if(0) // below not currently used

  //x2-dimension
  dimen = 2;
  

  //SASMARKXXX -- DO demand stationarity at the polar axis
  DIRSIGNLOOP(dirsign){
    dir = DIRFROMDIMEN( dimen, dirsign );

    j = dofluxreg[ACTIVEREGION][dir];

    //if boundary is not on this processor, do not modify emf's
    if( j < -MAXBND ) continue;
    if( j > N2-1+MAXBND ) continue;

    //the boundary is on the processor, so reset emf's to zero at the boundary
    COMPFULLLOOPP1_13{
      MACP1A0(emf,1,i,j,k) = 0.0;
      MACP1A0(emf,3,i,j,k) = 0.0;
    }
  }

  //x3-dimension
  dimen = 3;


  DIRSIGNLOOP(dirsign){
    dir = DIRFROMDIMEN( dimen, dirsign );

    k = dofluxreg[ACTIVEREGION][dir];

    //if boundary is not on this processor, do not modify emf's
    if( k < -MAXBND ) continue;
    if( k > N3-1+MAXBND ) continue;

    //the boundary is on the processor, so reset emf's to zero at the boundary
    COMPFULLLOOPP1_23{
      MACP1A0(emf,1,i,j,k) = 0.0;
      MACP1A0(emf,2,i,j,k) = 0.0;
    }
  }
#endif
  
  
}

// STAG version (different because emf vs. flux, where emf has extra zone -- otherwise same!
//reset emf's (E_3) at the boundaries in x1-x2 plane to zero
void user1_adjust_fluxctstag_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
{
  int i, j, k;
  int dir;
  int dimen;
  int dirsign;
  struct of_geom geom;
  struct of_geom *ptrgeom = &geom;
#if( N3 > 1 && DONSEMFS )
  FTYPE V_ph1[NDIM];
  FTYPE V_ph2[NDIM];
  FTYPE V_th1[NDIM];
  FTYPE V_th2[NDIM];
  FTYPE omega;
  FTYPE dflux;
  FTYPE aflux;
  FTYPE nflux;
#endif
  
  if(ADJUSTFLUXCT==0 && DOADJUSTEMFS==0){
    // then no need to set EMF's to zero to maintain divb=0, although can choose to still set EMF's to constant along appropriate directions to maintain stationarity all the way up to including the boundary. Not sure why that would be important to maintain.
    return;
  }

  if( DOGRIDSECTIONING==0 && DOADJUSTEMFS == 0 ) {  //don't do anything if sectioning not done
    return; //commented in order to fix the fields at inner radius
  }

  if(BCtype[X1DN]==FIXEDUSEPANALYTIC || BCtype[X1UP]==FIXEDUSEPANALYTIC 
  || BCtype[X1DN]==NSSURFACE){
    // only do if fixing BCs, not (e.g.) outflowing

    //x1-dimension
    dimen = 1;
    DIRSIGNLOOP(dirsign){
      dir = DIRFROMDIMEN( dimen, dirsign );


      if(dir==X1DN && BCtype[X1DN]==FIXEDUSEPANALYTIC || dir==X1UP && BCtype[X1UP]==FIXEDUSEPANALYTIC 
      || dir==X1DN && BCtype[X1DN]==NSSURFACE ){ // otherwise don't do

	//if boundary is not on this processor, do not modify emf's
	if( DOGRIDSECTIONING ) {
	  i = dofluxreg[ACTIVEREGION][dir];
	  if( i < -MAXBND ) continue;
	  if( i > N1-1+MAXBND) continue;
	}
	else if( totalsize[1] > 0 && mycpupos[1] == 0 && dir == X1DN ) {
	  i = 0;
	}
	else {
	  break;
	}

	
	//the boundary is on the processor, so reset emf's to zero at the boundary
	COMPFULLLOOPP1_23{
	  // EMF[3]:
#if(N1>1)
	  MACP1A1(fluxvec,1,i,j,k,B2) = 0.0;
	  MACP1A1(fluxvec,1,i,j,k,B1) = 0.0; // always zero
#endif
#if(N2>1)
	  MACP1A1(fluxvec,2,i,j,k,B1) = 0.0;
	  MACP1A1(fluxvec,2,i,j,k,B2) = 0.0; // always zero
#endif

	}
      }


      if( dir==X1DN && BCtype[X1DN]==NSSURFACE ){

	//if boundary is not on this processor, do not modify emf's
	if( DOGRIDSECTIONING ) {
	  i = dofluxreg[ACTIVEREGION][dir];
	  if( i < -MAXBND ) continue;
	  if( i > N1-1+MAXBND) continue;
	}
	else if( totalsize[1] > 0 && mycpupos[1] == 0 ) {
	  //work with the inner radial boundary
	  i = 0;
	}
	else {
	  break;
	}

	
	//the boundary is on the processor, so reset emf's to zero at the boundary
	COMPFULLLOOPP1_23{
	  // EMF[2]:
#if(N3>1 && DONSEMFS)
	  if( j >= 0 && j <= N2 && k >= 0 && k <= N3 ){
	    //get_geometry(i, j, k+1, CORN2, ptrgeom_ph2); //CORN2 -- "corner" in 1-3 plane, think this is where E_theta is located
	    bl_coord_ijk(i, j, k  , CORN2, V_ph1);
	    bl_coord_ijk(i, j, k+1, CORN2, V_ph2);
	    //get_geometry(i, j  , k, CORN3, ptrgeom_th1);
            //get_geometry(i, j+1, k, CORN3, ptrgeom_th2);
	    bl_coord_ijk(i, j,   k, CORN3, V_th1);
	    bl_coord_ijk(i, j+1, k, CORN3, V_th2);
	    //km1 = km1mac(k);
	    //km1 = max(km1, INFULL3);
	    omega = get_omegaf_phys(t, dt, steppart);
	    dflux = dfluxns(V_ph1[1], omega, V_ph1[3], V_th1[2], V_th2[2], t, dt);
	    //myB1 = 0.5 * ( GLOBALMACP0A1(pstagglobal,i,j,k,B1)+GLOBALMACP0A1(pstagglobal,i,j,km1,B1) );
	    //d(gdet*B1)/dt = -dF3(B1)/dx3
	    //dflux = d(gdet*B1*dx2*dx3) = -dF3(B1)*dx2*dt  
	    //F3(B1) = dflux / (dx2*dt) <-- make sure sign correct
#if(0)
	    if( 1 || j == 1 && k == 1 ){  
	      get_geometry(i, j, k  , FACE1, ptrgeom); //CORN2 -- "corner" in 1-3 	    //get_geometry(i, j, k  , CORN2, ptrgeom_ph1); //CORN2 -- "corner" in 1-3 plane, think this is where E_theta  is   located
	      aflux = vpotns_flux(V_ph1[1],V_th1[2],V_th2[2],V_ph1[3]-omega*t,V_ph2[3]-omega*t);
	      nflux = ptrgeom->gdet * GLOBALMACP0A1(pstagglobal,i,j,k,B1)*dx[2]*dx[3]; //MACP0A1(prim,i,j,k,B1)
	      dualfprintf(fail_file, "t = %21.15g, nstep = %ld, sp = %d, fluxvec[3][%d][%d][%d][B1] = %g, emf = %g, aflux = %g, nflux = %g\n",
			  t, nstep, steppart, i, j, k, MACP1A1(fluxvec,3,i,j,k,B1), dflux / (dx[2] * dt),
			  aflux, nflux);
	      //dualfprintf(fail_file,"r=%g, th1=%g, th2=%g, ph1=%g, ph2=%g\n", V_ph1[1], V_th1[2], V_th2[2], V_ph1[3], V_ph2[3]);
	    
	    }
#endif
	    MACP1A1(fluxvec,3,i,j,k,B1) = dflux / (dx[2] * dt);   // rotation, E_2 = (-[v x B])_2 = - v^3 B^1
#if(N1>1)
	    //adjust "symmetric" flux F1[B3] = - F3[B1] (see fluxvpot.c:set_emfflux())
	    MACP1A1(fluxvec,1,i,j,k,B3) = -MACP1A1(fluxvec,3,i,j,k,B1);
#endif
	    MACP1A1(fluxvec,3,i,j,k,B3) = 0.0; // always zero
	  }
#endif
	}
      }

      if(dir==X1UP && BCtype[dir]==FIXEDUSEPANALYTIC){ // i.e. don't reset EMF2=0 for lower boundary since that corresponds to rotation

	//if boundary is not on this processor, do not modify emf's
	i = dofluxreg[ACTIVEREGION][dir];
	if( i < -MAXBND ) continue;
	if( i > N1-1+MAXBND) continue;


	COMPFULLLOOPP1_23{
	  // EMF[2]:
#if(N1>1)
	  MACP1A1(fluxvec,1,i,j,k,B3) = 0.0;
	  MACP1A1(fluxvec,1,i,j,k,B1) = 0.0; // always zero
#endif
#if(N3>1)
	  MACP1A1(fluxvec,3,i,j,k,B1) = 0.0;
	  MACP1A1(fluxvec,3,i,j,k,B3) = 0.0; // always zero
#endif
	}
      }


    }
  }

  
  
}


FTYPE get_omegaf_prefactor( FTYPE t_transition, FTYPE t, FTYPE dt, FTYPE steppart )
{
  return(1.);

#if(0)  //too hard to keep track of varying omega, so keep it constant
  
  FTYPE sin1st, sin2nd, sin4th;
  FTYPE tsteppart;
  FTYPE prefact = 1.;     //no prefactor by default
  
  if(TIMEORDER==2){ 
    tsteppart = t + 0.5 * dt; //tsteppartf;  //time of the end of the current substep so that the boundaries are OK
    
    if( 1 == steppart ) {
      tsteppart += 0.5 * dt;
    }
  }
  else{
    tsteppart=t;
  }
  
  // trace back in time to get the advected value of Omega (assume moves at speed of light)
  //if( 1 ) { 
  //    coord( i, j, k, geom->p, X );
  //   bl_coord( X, V );
  //    r = V[1];
  //    tsteppart -= (r-Rin);  //trace back in time to the star surface r = Rin
  //}
  
  tsteppart -= 0.1;  //wait for deltat = 0.1 to avoid startup effects; assume dt < 0.1
  
  
  if( tsteppart < 0 ) {
    tsteppart = 0.0;
  }
  
  
  if( tsteppart < t_transition ) {  //smoothly switch on the rotation with time
    sin1st = sin( (tsteppart/t_transition-1) * M_PI_2l );
    sin2nd = sin1st * sin1st;
    sin4th = sin2nd * sin2nd;
    prefact = 1 - sin4th;
  }
  
  return( prefact );
#endif
  
}




// adjust final fluxes
// NOT USED
void user1_adjust_final_fluxes(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
{
  int i, j, k;
  int dir;
  int dimen;
  int dirsign;
  
  if( DOGRIDSECTIONING==0) {  //don't do anything if sectioning not done
    return;
  }

  
  //x1-dimension
  dimen = 1;
  DIRSIGNLOOP(dirsign){
    dir = DIRFROMDIMEN( dimen, dirsign );


    //if boundary is not on this processor, do not modify emf's
    i = dofluxreg[ACTIVEREGION][dir];
    if( i < -MAXBND ) continue;
    if( i > N1-1+MAXBND) continue;
    
    
    //the boundary is on the processor, so reset emf's to zero at the boundary
    COMPFULLLOOPP1_23{
      // EMF[3]:
#if(N1>1)
      MACP1A1(fluxvec,1,i,j,k,B2) = 0.0;
      MACP1A1(fluxvec,1,i,j,k,B1) = 0.0; // always zero
#endif
#if(N2>1)
      MACP1A1(fluxvec,2,i,j,k,B1) = 0.0;
      MACP1A1(fluxvec,2,i,j,k,B2) = 0.0; // always zero
#endif
      // EMF[2]:
#if(N1>1)
      MACP1A1(fluxvec,1,i,j,k,B3) = 0.0;
      MACP1A1(fluxvec,1,i,j,k,B1) = 0.0; // always zero
#endif
#if(N3>1)
      MACP1A1(fluxvec,3,i,j,k,B1) = 0.0;
      MACP1A1(fluxvec,3,i,j,k,B3) = 0.0; // always zero
#endif
    }
  }


  
  
}





// check that \detg=0 on singularities
// Don't use |\theta-0|<small number since can be quite small yet not on singularity
// use integer-based grid position based detection as consistent with bondary conditions
// Called for fresh and restart run
void check_spc_singularities_user(void)
{
  int i,j,k;
  int poleloop;
  int needzero;
  int singfound;
  int nonsingfound;
  int numlocs,indloc,loc;
  int pl,pliter;
  LOCALMETRICTEMPVARS;
  int doprintout;

  if(nstep==0) doprintout=1;
  else doprintout=0; // avoid print out when evolving metric and recalling this -- assume faily similar situation as at nstep==0

  ///////////////////////////
  //
  // X1 inner r=0 singularity
  //
  ///////////////////////////


  /* inner radial BC (preserves u^t rho and u) */
  if(mycpupos[1] == 0 && BCtype[X1DN]==R0SING){
    i=0;
      
    //    LOOPX1dir{
    // No, should really be full loops
    LOOPF_23{
    
      for(indloc=0;indloc<numlocs;indloc++){

	if(indloc==0) loc=FACE1;
	else if(indloc==1) loc=CORN3;
	else if(indloc==2) loc=CORN2;
	else if(indloc==3) loc=CORNT;

	// get local metric quantities for this loc,i,j,k
	GETLOCALMETRIC(loc,i,j,k);

	singfound=0;
	singfound+=(localgdet[0]!=0.0);
#if(WHICHEOM!=WITHGDET)
	PLOOP(pliter,pl) singfound+=(LOCALEOMFUNCMAC(pl)!=0.0);
#endif

	if(singfound){
	  if(doprintout) dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d loc=%d with gdet=%21.15g so resetting it to 0.0\n",i,j,k,loc,localgdet[0]);
	  localgdet[0]=0.0;
#if(WHICHEOM!=WITHGDET)
	  PLOOP(pliter,pl) LOCALEOMFUNCMAC(pl)=0.0;
#endif
	}
      }// end over indloc
    }// end loop 23
  }



  ///////////////////////////
  //
  // X2 POLARAXES
  //
  ///////////////////////////



  if(special3dspc && COORDSINGFIX>0){
    // Then presume really doing 3D problem with \phi direction
    // Then presume that (for simplicity of the method) offsetting pole by \epsilon from \theta=0,\pi so that \detg is generally non-zero
    // Then flux through pole, but transmissive BCs at pole ensure flux into "\epsilon hole" are same as out of the hole on other side
    // In any case that flux is very small.  Key point is that for d_t(\detg B2) = -d_1(\detg E3) + d_3(\detg E1), the small angular term in \detg cancels leaving regular B2 on the pole that enters interpolations of B2.

    // Check that \detg is NOT zero now
    needzero=0;
  }
  else{
    // Then presume no \epsilon offset at pole
    needzero=1;
  }



  for(poleloop=0;poleloop<=1;poleloop++){
    if(poleloop==0 && mycpupos[2] == 0 && BCtype[X2DN]==POLARAXIS){
      j=0;
    }
    else if(poleloop==1 && mycpupos[2] == ncpux2-1 && BCtype[X2UP]==POLARAXIS){
      j=N2-1+SHIFT2;
    }
    else continue;
    
    
    //    LOOPX2dir{
    // No, should really be full loops:
    LOOPF_13{

      for(indloc=0;indloc<numlocs;indloc++){

	if(indloc==0) loc=FACE2;
	else if(indloc==1) loc=CORN3;
	else if(indloc==2) loc=CORN1;
	else if(indloc==3) loc=CORNT;

	// get local metric quantities for this loc,i,j,k
	GETLOCALMETRIC(loc,i,j,k);

	singfound=0;
	singfound+=(localgdet[0]!=0.0);
#if(WHICHEOM!=WITHGDET)
	PLOOP(pliter,pl) singfound+=(LOCALEOMFUNCMAC(pl)!=0.0);
#endif

	nonsingfound=0;
	nonsingfound+=(fabs(localgdet[0])<=NUMEPSILON);
#if(WHICHEOM!=WITHGDET)
	PLOOP(pliter,pl) nonsingfound+=(fabs(LOCALEOMFUNCMAC(pl))<=NUMEPSILON);
#endif

	if(needzero && singfound){
	  if(doprintout) dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d loc=%d with gdet=%21.15g so resetting it to 0.0\n",i,j,k,loc,localgdet[0]);
	  localgdet[0]=0.0;
#if(WHICHEOM!=WITHGDET)
	  PLOOP(pliter,pl) LOCALEOMFUNCMAC(pl)=0.0;
#endif
	}
	else if(needzero==0 && nonsingfound){
	  dualfprintf(fail_file,"Detected |\\detg|<=NUMEPSILON at i=%d j=%d k=%d loc=%d with gdet=%21.15g: Must be non-zero\n",i,j,k,loc,localgdet[0]);
	  myexit(9813523);
	}
      }// end over indloc
    }// end loop 13
  }// end over poles

}


int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr)
{
  int i,j,k;
  int pliter, pliter2, pl, pl2;
  FTYPE vcon[NDIM]; // coordinate basis vcon
  extern FTYPE Omegastar;
  extern FTYPE Bpole;
  extern FTYPE Vpar;
  FTYPE prnew[NPR];
  FTYPE Bcon[NDIM];
  extern int OBtopr_general(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr_general2(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr_general3(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr_general3p(FTYPE omegaf, FTYPE v0, FTYPE *Bccon, struct of_geom *geom, FTYPE *pr);
  //FTYPE X[NDIM], V[NDIM], r;
  extern FTYPE t_transition;
  
  PBOUNDLOOP(pliter,pl){
    if(!isfinite(pr[pl]) && (pl!=U1 && pl!=U2 && pl!=U3)){ //check if any of the primitives other than velocity is a nan
      //don't do anything
      return(0);
    }
  }
  
  i=geom->i;
  j=geom->j;
  k=geom->k;
  
  
  
  Bcon[0]=0;
  Bcon[1]=pr[B1];
  Bcon[2]=pr[B2];
  Bcon[3]=pr[B3];
  
  
#if(EOMTYPE==EOMFFDE) // what would be done in force-free with no parallel velocity
  
  //	  dualfprintf(fail_file,"Omegastar=%21.15g dxdxp[3][3]=%21.15g B1=%21.15g B2=%21.15g B3=%21.15g\n",Omegastar,dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);
  
  ///////////////////////////////
  //
  // new way to get velocity
  // surface rotates with angular frequency Omegastar to observer at infinity
  if(OBtopr_general(omegaf,Bcon,geom,prnew)>=1){
    dualfprintf(fail_file, "OBtopr(bounds): space-like error in init_postfield()\n");
    dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
    failed=1;
    return(1);
  }
  // assign answer
  pr[U1]=prnew[U1];
  pr[U2]=prnew[U2];
  pr[U3]=prnew[U3];
  
  //	  if(t>1.9 && t<2.1){
  //dualfprintf(fail_file,"t=%21.15g i=%d j=%d\n",t,i,j);
  //  dualfprintf(fail_file,"newus: %21.15g %21.15g %21.15g\n",prnew[U1],prnew[U2],prnew[U3]);
  //  dualfprintf(fail_file,"Omegastar'=%21.15g Bcon1=%21.15g Bcon2=%21.15g Bcon3=%21.15g\n",Omegastar/dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);
  // }
  
#endif
  
#if(0)
  // set  	    ucon[TT,etc.]
  vcon[RR]=0; // surface that completely dissipates normal direction momentum
  vcon[TH]=0; // "" for this component
  // below assumes no phi mixing with other directions in grid
  vcon[PH]=omegaf; // surface rotates with angular frequency Omegastar to observer at infinity
  
  // get in terms of primitive velocity
  MYFUN(vcon2pr(WHICHVEL,vcon,geom,pr),"bounds.ns.c:bound_prim_user()", "vcon2pr()", 1);
#endif
  
#if(EOMTYPE!=EOMFFDE) 
  
  ///////////////////////////////
  //
  // new way to get velocity
  // surface rotates with angular frequency Omegastar to observer at infinity
  if(EOMTYPE==EOMCOLDGRMHD && MCOORD==KSCOORDS) {
    if(OBtopr_general3(omegaf,vpar,Bcon,geom,prnew)>=1){
      dualfprintf(fail_file, "OBtopr(bounds): space-like error in init_postfield()\n");
      dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
      failed=1;
      return(1);
    }
  }
  else{
    //if(OBtopr_general3p(omegaf,vpar,Bcon,geom,prnew)>=1){  //set poloidal velocity
    if(OBtopr_general(omegaf,Bcon,geom,prnew)>=1){ //pure force-free
    //if(OBtopr_general2(omegaf,vpar,Bcon,geom,prnew)>=1){  //pure force-free + vparallel
      dualfprintf(fail_file, "OBtopr(bounds): space-like error in init_postfield()\n");
      dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
      failed=1;
      return(1);
    }
  }
  // assign answer
  pr[U1]=prnew[U1];
  pr[U2]=prnew[U2];
  pr[U3]=prnew[U3];
  
  //dualfprintf( fail_file, "Actual vpar = %21.15g\n", pr[U1] );
  
#endif
  
  
  PBOUNDLOOP(pliter,pl){
    if(!isfinite(pr[pl])){
      dualfprintf(fail_file,"not finite in set_vel_stataxi: i=%d j=%d steppart=%d nstep=%ld p = %d omegaf = %g vpar = %g\n",i+startpos[1],j+startpos[2],steppart,nstep, geom->p, omegaf, vpar);
      PBOUNDLOOP(pliter2,pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,pr[pl2]);
      fflush( stderr );
      fflush( stdout );
      myexit(2525);
    }
  }
  
  
  
  return(0);
  
}



// DEBUG special3dspc
// check matching CPUs have correct information across pole
// stick in debugspecial3dspc(0, whichdir, ispstag, prim); calls before and after various stages defined through (e.g.) which=0,1,2,3,4,...
int debugspecial3dspc(int which, int whichdir, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pliter,pl;


  dualfprintf(fail_file,"nstep=%ld steppart=%d which=%d whichdir=%d ispstag=%d\n",nstep,steppart,which,whichdir,ispstag);
  for(pl=0;pl<=7;pl++){
    if(ispstag==1 && (pl<5 || pl>7)) continue;
    dualfprintf(fail_file,"pl=%d\n",pl);
    for(k=-N3BND;k<=N3BND;k++){ // was looking at just k=0, but seeing now if all other k are good.
      dualfprintf(fail_file,"k=%d\n",k);
      dualfprintf(fail_file,"%02s | i=\n   ","j");
      for(i=-N1BND;i<=N1BND;i++){
	if(ispstag==1 && (pl==5 && i==-N1BND)) dualfprintf(fail_file,"%021s ","NA");
	else dualfprintf(fail_file,"%19s%02d "," ",i);
      }
      dualfprintf(fail_file,"\n");
      for(j=-N2BND;j<=N2BND;j++){
	dualfprintf(fail_file,"%02d ",j);
	for(i=-N1BND;i<=N1BND;i++){
	  if(ispstag==1 && (pl==5 && i==-N1BND || pl==6 && j==-N2BND)) dualfprintf(fail_file,  "%21s ","NA");
	  else  dualfprintf(fail_file,  "%21.15g ",MACP0A1(prim,i,j,k,pl));
	}
	dualfprintf(fail_file,"\n");
      }
      dualfprintf(fail_file,"\n");
    }
  }// end over pl


  return(0);

}
