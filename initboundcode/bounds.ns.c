
#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: something seriously wrong with EXTRAP=1 (EOMFFDE)

// FILE SIMILAR TO for fishmon but with NSSURFACE for inner radial boundary


// 0: original fix B^r, \Omega_F in ghost zones only (no flux bound using plr), outflow rest, and clean B^\theta in ghost zones
// 1: as 0, but also set fluxes through pl and pr using quasi-analytical face
// 2: as 1, but also choose boundary analytical values to be extrapolated through analytical boundary value at surface
// 3: as 2, but fix v^i not parallel inside star instead of extrapolating it


// whether to set flux of NS using face values (analytical + outflowed)
#define SETNSFLUXRHO 1
#define SETNSFLUXUU 1
#define SETNSFLUXV1 1 // screws up B2 near equator
#define SETNSFLUXV2 0 // causes problems for B2 deep inside star near poles
#define SETNSFLUXV3 1
#define SETNSFLUXB1 0 // causes problems for B2 near star near poles
#define SETNSFLUXB2 1
#define SETNSFLUXB3 1

//int setnsflux[NPR]={SETNSFLUXRHO,SETNSFLUXUU,SETNSFLUXV1,SETNSFLUXV2,SETNSFLUXV3,SETNSFLUXB1,SETNSFLUXB2,SETNSFLUXB3};
int setnsflux[NPR]={0,0,0,0,0,0,0,0};
//int setnsflux[NPR]={1,1,  1,1,1,   1,1,1};
//int setnsflux[NPR]={0,0,0,1,0,0,0,0}; // V2 problem
//int setnsflux[NPR]={1,1,1,0,0,0,0,0}; // V1 problem
//int setnsflux[NPR]={0,0,0,0,0,1,0,0}; // B1 problem
//int setnsflux[NPR]={1,1,1,1,1,0,0,0};

// whether to
// 2: set via stationarity condition (\rho u^p / B^p = const[face value])
// 1: set ghost zones analytically
// 0: extrapolate through face
#define SETGHOSTRHO 1

// whether to
// 1: set ghost zones analytically
// 0: extrapolate through face
#define SETGHOSTUU 1

// whether to
// 1: set ghost zones analytically (i.e. stationarity conditions with fixed vpar)
// 0: extrapolate through face (no constraint using stationarity -- just extrapolate through -- vpar not constrained)
#define SETGHOSTV 1

// whether to
// 1: set ghost zones analytically
// 0: extrapolate through face
#define SETGHOSTB1 1

// whether to clean B^\theta for star
#define CLEANBH 0

#define EXTRAP 1
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
//#define POLEDEATH N2BND
#define POLEDEATH 0
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole


// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPN2 LOOPN3 LOOPBOUNDIN1 & LOOPBOUNDOUT1
// X2 in&out: LOOPF1 LOOPN3 LOOPBOUNDIN2 & LOOPBOUNDOUT2  // LOOPF1 ok if X1 dimension not there, then LOOPF1->LOOPN1
// X3 in&out: LOOPF1 LOOPF2 LOOPBOUNDIN3 & LOOPBOUNDOUT3  // as above

int bound_prim_user(int boundstage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  void basic_outflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  //
  void bound_field_outflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void compute_aphi_fromoutflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*aphifromoutflow)[NSTORE2][NSTORE3]);
  void bound_x2_aphi(FTYPE (*aphifromoutflow)[NSTORE2][NSTORE3]);
  void compute_btheta_fromaphi(FTYPE (*aphifromoutflow)[NSTORE2][NSTORE3], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void clean_btheta_x1inner(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void clean_btheta_x1outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

  int bound_vel_from_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

  int bound_x1_outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x2_inner_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x2_outer_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x3_periodic(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  FTYPE BASEMAC(aphifromoutflow,N1BND+1,N2M,N3M); // only inner boundary ghost zones plus 1 active zone right on boundary
  FTYPE PTRDEFGLOBALMAC(aphifromoutflow,N1BND+1,N2M,N3M);
  int i,j,k;


  // pointer shift
  GLOBALPOINT(aphifromoutflow) = PTRMAC(aphifromoutflow,N1BND+1,N2M,N3M) (&(BASEMAC(aphifromoutflow,N1BND,N2BND,N3BND)));

  //  dualfprintf(fail_file,"primmem: %ld %d\n",prim,prim);

#if(0)
  // SUPERGODMARK
  if(t>1.0 && t<3.0){
    FULLLOOP{
      MACP0A1(prim,i,j,k,U1)=0.0;
      MACP0A1(prim,i,j,k,U2)=0.0;
      MACP0A1(prim,i,j,k,U3)=0.0;
    }
  }
#endif

  // polar axis bound so field is bounded fully
  x2_inner_polar(prim);
  x2_outer_polar(prim);

  // basic outflow
  basic_outflow(prim);

  // polar axis bound so field is bounded fully
  x2_inner_polar(prim);
  x2_outer_polar(prim);

   
      
  // first outflow B^\theta and B^\phi however approximately wanted
  //  dualfprintf(fail_file,"B1\n");
  bound_field_outflow(prim);

  // polar axis bound so field is bounded fully
  //      dualfprintf(fail_file,"B2\n");
  x2_inner_polar(prim);
  //      dualfprintf(fail_file,"B3\n");
  x2_outer_polar(prim);

      
#if(0)
  // next compute vector potential (A_\phi) at corners from analytical B^r at FACE1 and outflowed B^\theta at center
  // gives all aphi needed, so don't need to bound aphi
  //      dualfprintf(fail_file,"B4\n");
  compute_aphi_fromoutflow(prim,GLOBALPOINT(aphifromoutflow));

  // compute B^\theta from A_\phi
  //      dualfprintf(fail_file,"B5\n");
  compute_btheta_fromaphi(GLOBALPOINT(aphifromoutflow),prim);
#endif
#if(CLEANBH)
  clean_btheta_x1inner(prim);
#endif

  // polar axis bound after B^\theta adjusted
  //  dualfprintf(fail_file,"B6\n");
  x2_inner_polar(prim);
  //      dualfprintf(fail_file,"B7\n");
  x2_outer_polar(prim);

  // compute v^i from B^i and \Omega_F
  //      dualfprintf(fail_file,"B8\n");
#if(0)
  //  Maybe only want to extrapolate through B1 and still fix v^i?

  // if setting velocity by extrapolating \tilde{u}^i from FACE1 value, then don't do this
  // should see how makes difference
  MYFUN(bound_vel_from_field(prim),"bounds.ns.c:bound_prim_user()", "bound_vel_from_field()", 1);
#endif



  /////////// x1

  //      dualfprintf(fail_file,"B9\n");
  MYFUN(bound_x1_outer(prim),"bounds.ns.c:bound_prim_user()", "bound_x1_outer()", 1);

  // bound x2 again since needed for outer radial boundary to be consistent
  //      dualfprintf(fail_file,"B10\n");
  x2_inner_polar(prim);
  //      dualfprintf(fail_file,"B11\n");
  x2_outer_polar(prim);

#if(CLEANBH)
  clean_btheta_x1outer(prim);
#endif
  // polar bound after B^\theta adjusted
  x2_inner_polar(prim);
  x2_outer_polar(prim);



  // x3
  //      dualfprintf(fail_file,"B12\n");
  x3_periodic(prim);
  //  dualfprintf(fail_file,"B13\n");


  return (0);
}





void basic_outflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;

  LOOPF2 LOOPF3{
    ri=0;
    rj=j;
    rk=k;


    
    // loop over ghost zones inside star
    LOOPBOUND1IN{
      PLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl);
    }
  }

  LOOPF2 LOOPF3{
    ri=N1-1;
    rj=j;
    rk=k;


    
    // loop over outer zones
    LOOPBOUND1OUT{
      PLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl);
    }
  }


}


///////////////////////////
//
// X1 inner OUTFLOW/FIXEDOUTFLOW
//
///////////////////////////
// first outflow B^\theta and B^\phi however approximately wanted
void bound_field_outflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;
  struct of_geom geom,rgeom,rrgeom;
  FTYPE X[NDIM],V[NDIM],rX[NDIM],rV[NDIM],dxdxp[NDIM][NDIM],rdxdxp[NDIM][NDIM];
  struct of_geom fgeom;
  FTYPE fX[NDIM],fV[NDIM],fdxdxp[NDIM][NDIM];
  FTYPE rrX[NDIM],rrV[NDIM],rrdxdxp[NDIM][NDIM];
  FTYPE Bcon[NDIM],prnew[NPR];
  FTYPE Bconnew[NDIM];
  extern void getnewucon(FTYPE *uconmetpin, FTYPE *rV, struct of_geom *rptrgeom, FTYPE (*rdxdxp)[NDIM], FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetpout);
  FTYPE v1,v2,v3,myA,myB,slope,newv;
  extern int OBtopr_general2(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  FTYPE prface[NPR];
  void set_face1(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface);
  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);
  FTYPE up2face,Bp2face,up2cent,Bp2cent;
  int jj,kk;




  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==NSSURFACE)){




      LOOPF2 LOOPF3{
        ri=0;
        rj=j;
        rk=k;

        // get reference grid parameters
        get_geometry(ri, rj, rk, CENT, &rgeom);


#if(1) // usually need face values -- don't HAVE to use them
        get_geometry(ri, rj, rk, FACE1, &fgeom);
        coord(ri, rj, rk, FACE1, fX);
        bl_coord( fX, fV );
        dxdxprim(fX, fV, fdxdxp);

        // get analytical solution at NS surface
        set_face1(ri,rj,rk,prim,prface);

        //PLOOP(pliter,pl) dualfprintf(fail_file,"ri=%d rj=%d rk=%d pl=%d prface=%21.15g\n",ri,rj,rk,pl,prface[pl]);
#endif




        // loop over ghost zones inside star
        LOOPBOUND1IN{


          // reference geometry
          coord(ri, rj, rk, CENT, rX);
          bl_coord( rX, rV );
          dxdxprim(rX, rV, rdxdxp);
          get_geometry(ri, rj, rk, CENT, &rgeom);

          // ghost cell geometry
          coord(i, j, k, CENT, X);
          bl_coord( X, V );
          dxdxprim(X, V, dxdxp);
          get_geometry(i, j, k, CENT, &geom);





          // below not used yet
          ///   pl=UU; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          pl=UU; MACP0A1(prim,i,j,k,pl) = 0.0;


          // outflow field, NS field completely reconnects through surface
          //   for(pl=B1;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl);

          // outflow B^\theta B^\phi
          // SUPERGODMARK : removed B2 -- worried about divb ... probably just normalization is crazy
          //   MACP0A1(prim,i,j,k,B2) = GLOBALMACP0A1(panalytic,i,j,k,B2);

          //   for(pl=B2;pl<=B3;pl++){
          //  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          // }
          //pl = B2;
          // for dipole, constant is B^\theta \propto 1/(\detg r^2)
          // below seems best
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)*rV[1]*rV[1]/(GLOBALMACP1A0(gdet,CENT,i,j,k)*V[1]*V[1]) ;
          // seems to do ok
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/(GLOBALMACP1A0(gdet,CENT,i,j,k)) ;
          // does worst as for polar artifacts
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) ;

#if(0)


          // reference field
          Bcon[0]=0;
          Bcon[1]=MACP0A1(prim,ri,rj,rk,B1);
          Bcon[2]=MACP0A1(prim,ri,rj,rk,B2);
          Bcon[3]=MACP0A1(prim,ri,rj,rk,B3);

          getnewucon(Bcon, rV, &rgeom, rdxdxp, V, &geom, dxdxp, Bconnew);

          MACP0A1(prim,i,j,k,B2)=Bconnew[TH];
          MACP0A1(prim,i,j,k,B3)=Bconnew[PH];
#endif
#if(0)
          pl = B3;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/(GLOBALMACP1A0(gdet,CENT,i,j,k)) ;

          pl=B2;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;

          //pl=B3;
          //MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;


          //   for(pl=B2;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#endif
#if(1)
          // extrapolation must be consistent with face method


          // MOST ROBUST (choice for B2 doesn't matter if using cleaning and using FLIPGDETAXIS 0
          pl = B3;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          pl=B2;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          //MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          // try supressing rise near pole/star of B^\theta (odd that makes one side (pi-pole) large + and then - value with cleaning?
          //      MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k)*fabs(sin(V[2])) ;
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl); // if cleaning, doesn't matter what this is except the startingn value that is solvable if FLIPGDETAXIS=0 and not otherwise
          // below only matters if FLIPGDETAXIS==1 since otherwise one solves for these values in the end anyways
          //   if(startpos[2]+j==0) MACP0A1(prim,i,j,k,pl)=0.0; // crushing regularization
          //   if(startpos[2]+j==totalsize[2]-1) MACP0A1(prim,i,j,k,pl)=0.0; // crushing regularization

#endif
#if(0)

          if(sign(MACP0A1(prim,ri,rj,rk,B3))==sign(-prface[B1]*GLOBALMACP1A0(pother,OMEGAFFACE1,ri,rj,rk))){
            // then normal situation so can extrapolate 
            MACP0A1(prim,i,j,k,B3) = MACP0A1(prim,ri,rj,rk,B3)*pow(rV[1]/V[1],3);
          }
          else{ // then assume out of equilibrium so don't extrapolate
            MACP0A1(prim,i,j,k,B3) = MACP0A1(prim,ri,rj,rk,B3);
          }

          pl=B2;
          //    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          // pretty good
          //    MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          MACP0A1(prim,i,j,k,pl) = prface[pl] * (fgeom.g)/(geom.g);
          // for dipole mostly constant

          if(VARTOINTERPFIELD==PULSARFIELD2){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * dxdxp[2][2]*pow(V[1],4)/(rdxdxp[2][2]*pow(rV[1],4));
          }
          else if(VARTOINTERPFIELD==PULSARFIELD3){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(sqrt(fabs(geom.gcov[GIND(2,2)]))*pow(V[1],3.0))/(sqrt(fabs(rgeom.gcov[GIND(2,2)]))*pow(rV[1],3.0));
          }


#endif


          ////////////////////
          // fix quantities within entire surface to analytical values
          for(pl=UU;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          //   for(pl=UU;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
    



          ///////////////////////////////////////
          //
          // set B^r
          //
          ///////////////////////////////////////


#if(SETGHOSTB1==1)
          ///////////////////////////////////////////////////////////////////
          //
          /////////////// FIX SOME QUANTITIES to analytic solution
          //
          //////////////////////////////////////////////////////////////////
          // fix B^r
          // "analytical" value from FLUXCT
          //MACP0A1(prim,i,j,k,B1) = GLOBALMACP0A1(panalytic,i,j,k,B1);
          // analytical value from perfect differencing
          MACP0A1(prim,i,j,k,B1) = GLOBALMACP1A0(pother,B1CENT,i,j,k); // SUPERGODMARK
          //MACP0A1(prim,i,j,k,B1) = GLOBALMACP1A0(pother,B1FACE1,i,j,k); // SUPERGODMARK

#elif(SETGHOSTB1==0)

          //////////////////////////////
          //
          // extrapolate through surface analytical value
          //
          ////////////////////////////////


          // choose to extrapolate the same as in interpolation routines (flux.c) or something close


          // B1
          // should really be fixing orthonormal B^r and not B1 !  Otherwise problem changes when changing grid/resolution
          // if going to use other components, then really need to do this after B^\theta B^\phi are done changing (i.e. after B^\theta (B2?) cleaned)
          // or perhaps that's not so important and can clean B^\theta anyways and use origianl B2 as estimate to constrain B1 from B^r estimated B\theta
#if(0)
          v1=prface[B1]*sqrt(fabs(fgeom.gcov[GIND(1,1)]))*pow(fV[1],3); // close to constant for dipole
          v2=MACP0A1(prim,ri,rj,rk,B1)*sqrt(fabs(rgeom.gcov[GIND(1,1)]))*pow(rV[1],3);
          slope=(v2-v1)/(rV[1]-fV[1]);

          newv=slope*(V[1]-fV[1])+v1;
          MACP0A1(prim,i,j,k,B1) = newv/(sqrt(fabs(geom.gcov[GIND(1,1)]))*pow(V[1],3));
#elif(1)
          v1=prface[B1];
          v2=MACP0A1(prim,ri,rj,rk,B1);
          slope=(v2-v1)/(rV[1]-fV[1]);

          newv=slope*(V[1]-fV[1])+v1;
          MACP0A1(prim,i,j,k,B1) = newv;
#elif(0)   
          MACP0A1(prim,i,j,k,B1) = prface[B1];
#elif(0)
          // preserves flux if face has preserved flux
          // causes uu0>>1 on entire surface
          MACP0A1(prim,i,j,k,B1) = prface[B1]*(fgeom.g)/(geom.g);
#elif(0)

          // second reference geometry
          coord(ri+1, rj, rk, CENT, rrX);
          bl_coord( rrX, rrV );
          dxdxprim(rrX, rrV, rrdxdxp);
          get_geometry(ri+1, rj, rk, CENT, &rrgeom);

          v1=prface[B1]*sqrt(fabs(fgeom.gcov[GIND(1,1)]))*pow(fV[1],3); // close to constant for dipole
          v2=MACP0A1(prim,ri,rj,rk,B1)*sqrt(fabs(rgeom.gcov[GIND(1,1)]))*pow(rV[1],3);
          v3=MACP0A1(prim,ri+1,rj,rk,B1)*sqrt(fabs(rrgeom.gcov[GIND(1,1)]))*pow(rrV[1],3);

          myA=( (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) )/(fV[1]-rV[1]);
          myB=(v1-v2)/(fV[1]-rV[1]) + (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) ;

          newv=v1+myA*(V[1]-fV[1])*(V[1]-fV[1]) + myB*(V[1]-fV[1]);
          MACP0A1(prim,i,j,k,B1) = newv/(sqrt(fabs(geom.gcov[GIND(1,1)]))*pow(V[1],3));
#elif(0)

          // second reference geometry
          coord(ri+1, rj, rk, CENT, rrX);
          bl_coord( rrX, rrV );
          dxdxprim(rrX, rrV, rrdxdxp);
          get_geometry(ri+1, rj, rk, CENT, &rrgeom);

          v1=prface[B1];
          v2=MACP0A1(prim,ri,rj,rk,B1);
          v3=MACP0A1(prim,ri+1,rj,rk,B1);

          myA=( (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) )/(fV[1]-rV[1]);
          myB=(v1-v2)/(fV[1]-rV[1]) + (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) ;

          newv=v1+myA*(V[1]-fV[1])*(V[1]-fV[1]) + myB*(V[1]-fV[1]);
          MACP0A1(prim,i,j,k,B1) = newv;
#endif




#endif



          ///////////////////////////////////////
          //
          // set velocity (field should be set by now)
          //
          ///////////////////////////////////////



#if(SETGHOSTV==1)// set via stationary conditions

          // assume vparface1 is same as vparcent
          set_vel_stataxi(&geom, GLOBALMACP1A0(pother,OMEGAFCENT,i,j,k),GLOBALMACP1A0(pother,VPARFACE1,i,j,k),MAC(prim,i,j,k));

#elif(SETGHOSTV==0)// otherwise will set via stationarity conditions

          //   for(pl=U1;pl<=B3;pl++) dualfprintf(fail_file,"i=%d j=%d pl=%d prface=%21.15g :: %21.15g\n",i,j,pl,prface[pl],GLOBALMACP1A0(pother,VPARFACE1,i,j,k));

          // now can interpolate through primitive velocities and then RESET boundary values to satisfy stationary/axisymmetric conditions

          ///////////////////////////////
          // NEW: Don't worry about proxy values inside star being stationary.  Just make sure boundary velocity is stationary where flux is set
          // compare with setting velocity via field, omegaf, and vpar in standard vel-set function
          // \tilde{u}^{jj}
#if(0)
          SLOOPA(jj){
            v1=prface[U1+jj-1];
            v2=MACP0A1(prim,ri,rj,rk,U1+jj-1);
            slope=(v2-v1)/(rV[1]-fV[1]);
            newv=slope*(V[1]-fV[1])+v1;
            MACP0A1(prim,i,j,k,U1+jj-1) = newv;
          }
#elif(1)
          SLOOPA(jj){
            MACP0A1(prim,i,j,k,U1+jj-1) = prface[U1+jj-1];
          }
#endif
          //   PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d j=%d pl=%d pother=%21.15g\n",i,j,pl,GLOBALMACP1A0(pother,RHOFACE1+pl,ri,rj,rk));

#endif




          ///////////////////////////////////////
          //
          // set density
          //
          ///////////////////////////////////////



#if(SETGHOSTRHO==1)
   
          // fix \rho_0
          MACP0A1(prim,i,j,k,RHO) = GLOBALMACP0A1(panalytic,i,j,k,RHO);

#elif(SETGHOSTRHO==0)

          // rho
#if(0)
          v1=prface[RHO]*pow(fV[1],2); // assume density goes like 1/r^2
          v2=MACP0A1(prim,ri,rj,rk,RHO)*pow(rV[1],2);
          slope=(v2-v1)/(rV[1]-fV[1]);

          newv=slope*(V[1]-fV[1])+v1;
          MACP0A1(prim,i,j,k,RHO) = newv/(pow(V[1],2));
#elif(0)
          v1=prface[RHO];
          v2=MACP0A1(prim,ri,rj,rk,RHO);
          slope=(v2-v1)/(rV[1]-fV[1]);

          newv=slope*(V[1]-fV[1])+v1;
          MACP0A1(prim,i,j,k,RHO) = newv;
#elif(1)
          MACP0A1(prim,i,j,k,RHO) = prface[RHO];
#elif(0)

          // second reference geometry
          coord(ri+1, rj, rk, CENT, rrX);
          bl_coord( rrX, rrV );
          dxdxprim(rrX, rrV, rrdxdxp);
          get_geometry(ri+1, rj, rk, CENT, &rrgeom);

          v1=prface[RHO];
          v2=MACP0A1(prim,ri,rj,rk,RHO);
          v3=MACP0A1(prim,ri+1,rj,rk,RHO);

          myA=( (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) )/(fV[1]-rV[1]);
          myB=(v1-v2)/(fV[1]-rV[1]) + (v1-v3)/(fV[1]-rrV[1]) + (v3-v2)/(rV[1]-rrV[1]) ;

          newv=v1+myA*(V[1]-fV[1])*(V[1]-fV[1]) + myB*(V[1]-fV[1]);
          MACP0A1(prim,i,j,k,RHO) = newv;
#endif




#elif(SETGHOSTRHO==2)

          // rho
#if(1) // donor using stationary/axisymmetric conditions for density

          up2face=0.0;
          SLOOP(jj,kk){
            up2face+=prface[U1+jj-1]*prface[U1+kk-1]*fgeom.gcov[GIND(jj,kk)];
            //     dualfprintf(fail_file,"prface[%d]=%21.15g\n",jj,prface[U1+jj-1]);
          }
          Bp2face=0.0;
          SLOOP(jj,kk){
            Bp2face+=prface[B1+jj-1]*prface[B1+kk-1]*fgeom.gcov[GIND(jj,kk)];
            //     dualfprintf(fail_file,"prface[%d]=%21.15g\n",jj,prface[B1+jj-1]);
          }

          up2cent=0.0;
          SLOOP(jj,kk){
            up2cent+=MACP0A1(prim,i,j,k,U1+jj-1)*MACP0A1(prim,i,j,k,U1+kk-1)*geom.gcov[GIND(jj,kk)];
          }
          Bp2cent=0.0;
          SLOOP(jj,kk){
            Bp2cent+=MACP0A1(prim,i,j,k,B1+jj-1)*MACP0A1(prim,i,j,k,B1+kk-1)*geom.gcov[GIND(jj,kk)];
          }
   
          myA=sqrt(fabs(up2face)/fabs(Bp2face));
          myB=sqrt(fabs(up2cent)/fabs(Bp2cent));
          if(myB!=0.0 && Bp2face!=0.0 && Bp2cent!=0.0){
            MACP0A1(prim,i,j,k,RHO) = prface[RHO]*myA/myB;
          }
          else MACP0A1(prim,i,j,k,RHO)=prface[RHO];

          //   dualfprintf(fail_file,"i=%d j=%d myA=%21.15g myB=%21.15g up2face=%21.15g Bp2face=%21.15g up2cent=%21.51g Bp2cent=%21.15g\n",i,j,myA,myB,up2face,Bp2face,up2cent,Bp2cent);

#endif




#endif








        }
      }
      // now RHO, B^r, B^\theta, B^\phi, and U1 are set according to boundary conditions at surface
      // now can set v^i [other than v^r] and clean B^\theta if desired

    }
  }

}




// set prface for fluxes at boundary and interpolation through boundary value
// Analytical set of B1, omegaf, and vpar
void set_face1(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface)
{
  int pl,pliter;
  FTYPE Bcon[NDIM];
  extern int OBtopr_general2(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  FTYPE Xc[NDIM],Vc[NDIM];
  struct of_geom geomc;
  FTYPE dxdxpc[NDIM][NDIM];
  //
  FTYPE Xf[NDIM],Vf[NDIM];
  struct of_geom geomf;
  FTYPE dxdxpf[NDIM][NDIM];




  if(startpos[1]+i==0){

    // face geometry
    get_geometry(i, j, k, FACE1, &geomf);
    coord(i, j, k, FACE1, Xf);
    bl_coord( Xf, Vf );
    dxdxprim(Xf, Vf, dxdxpf);


    // center geometry
    get_geometry(i, j, k, CENT, &geomc);
    coord(i, j, k, CENT, Xc);
    bl_coord( Xc, Vc );
    dxdxprim(Xc, Vc, dxdxpc);


    prface[RHO] = GLOBALMACP1A0(pother,RHOFACE1,i,j,k);
    prface[UU] = GLOBALMACP1A0(pother,UUFACE1,i,j,k);

    ///////////////////////////////////
    //
    /////////// B1

    // set field at face1
    //    prface[B1] = GLOBALMACP1A0(pother,B1FACE1,i,j,k); // set analytically
    //prface[B1] = GLOBALMACP1A0(pother,B1FLUXCTFACE1,i,j,k); // set from FLUXCT-type definition at FACE1
    

    // best
    prface[B1] = GLOBALMACP0A1(panalytic,i,j,k,B1)*sqrt(fabs(geomc.gcov[GIND(1,1)]))*pow(Vc[1],3)/(sqrt(fabs(geomf.gcov[GIND(1,1)]))*pow(Vf[1],3)); // set from offset of analytical solution

    // leads to large uu0 on entire surface
    //    prface[B1] = MACP0A1(prim,i,j,k,B1);


    // this preserves divb since flux at every radii preserved
    // leads to large uu0 on entire surface
    //    prface[B1] = MACP0A1(prim,i,j,k,B1)*(geomc.g)/(geomf.g);

    // ok, but B2 has relatively large oscillatory value near pole-star interface
    //    prface[B1] = MACP0A1(prim,i,j,k,B1)*sqrt(fabs(geomc.gcov[GIND(1,1)]))*pow(Vc[1],3)/(sqrt(fabs(geomf.gcov[GIND(1,1)]))*pow(Vf[1],3));


    //    prface[B1] = 0.5*(GLOBALMACP0A1(panalytic,im1,j,k,B1)+GLOBALMACP0A1(panalytic,i,j,k,B1)); // set from offset of analytical solution


    ///////////////////////////////////
    //
    /////////// B2
    //    prface[B2] = MACP0A1(prim,i,j,k,B2);
    // outflow using gdet
    //    if(MACP0A1(prim,i,j,k,B2)>0.0){
    prface[B2] = MACP0A1(prim,i,j,k,B2)*(geomc.g)/(geomf.g);
    // }
    //else{
    //  prface[B2] = 0.0; // limit
    // }
    
    ///////////////////////////////////
    //
    /////////// B3
    if(VARTOINTERPFIELD==PULSARFIELD2){
      prface[B2] = MACP0A1(prim,i,j,k,B2)*(dxdxpc[2][2]*pow(Vc[1],4.0))/(dxdxpf[2][2]*pow(Vf[1],4.0));
    }
    else if(VARTOINTERPFIELD==PULSARFIELD3){
      prface[B2] = MACP0A1(prim,i,j,k,B2)*(sqrt(fabs(geomc.gcov[GIND(2,2)]))*pow(Vc[1],3.0))/(sqrt(fabs(geomf.gcov[GIND(2,2)]))*pow(Vf[1],3.0));
    }

    // outflow as copy
    if(sign(MACP0A1(prim,i,j,k,B3))==sign(-prface[B1]*GLOBALMACP1A0(pother,OMEGAFFACE1,i,j,k))){
      // then normal situation so can extrapolate 
      prface[B3] = MACP0A1(prim,i,j,k,B3)*pow(Vc[1]/Vf[1],3);
    }
    else{ // then assume out of equilibrium so don't extrapolate
      prface[B3] = MACP0A1(prim,i,j,k,B3);
    }
    
    // get \tilde{u}^i @ FACE1 from field at face
    Bcon[0]=0;
    Bcon[1]=prface[B1];
    Bcon[2]=prface[B2];
    Bcon[3]=prface[B3];
    OBtopr_general2(GLOBALMACP1A0(pother,OMEGAFFACE1,i,j,k), GLOBALMACP1A0(pother,VPARFACE1,i,j,k), Bcon, &geomf, prface);
      

    // now set pother so can be used by rest of routines in any other functions
    PLOOP(pliter,pl){
      GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k)=prface[pl];
    }

  }




}




void set_plpr(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r)
{
  int pl,pliter;


  // need a new function that sets quantities at the face1 to be used when setting p_l and p_r
  // need to set v^i at the FACE1
  if(startpos[1]+i==0  && dir == 1){
    // then at stellar surface for this flux
    //    dualfprintf(fail_file,"nstep=%ld steppart=%d j=%d\n",nstep,steppart,j);
    PLOOP(pliter,pl){
 
      if(setnsflux[pl]){
        if( 0 && pl == U1 ) {
          dualfprintf(fail_file,"n = %ld sp = %d j = %d pl=%d p_l=%21.15g p_r=%21.15g pother=%21.15g dp/p = %21.15g, dplpr/p = %21.15g\n",
                      nstep,steppart,startpos[2]+j,pl,p_l[pl],p_r[pl],GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k), 
                      fabs((p_l[pl]-GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k))/GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k)),
                      2*fabs((p_l[pl]-p_r[pl])/(p_r[pl]+p_l[pl])) 
                      ); 
        }

        p_l[pl]=p_r[pl]=GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k); // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
 
        //p_l[pl]=p_r[pl]; // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
        //p_r[pl]=p_l[pl];
      }

    }
  }

}



// start from poles since error accumulated at equator will be smaller
// seems to work!
// && FLIPGDETAXIS is because if gdet flips and B2 flips, then term solving for cancels across polar axis.  Otherwise can solve.  SUPERGODMARK -- sure?  can always solve even if would have cancelled!
void clean_btheta_x1inner(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE divbr;
  struct of_geom geompp,geompm,geommp,geommm;
  int ip,im,jp,jm;



  k=0;
  for(i=-1;i>=-N1BND;i--) for(j=0;j<=OUTM2;j++){ // going from down to up

      if(startpos[2]+j>=totalsize[2]/2){
      }
      else if(startpos[2]+j == 0 && FLIPGDETAXIS==0){
        // if FLIPGDETAXIS, then can't solve for near-axial value of B2
        // if FLIPGDETAXIS=0, then can solve but must solve such that near-axial values are same instead of just original polar boundary value
        // here pair is MACP0A1(prim,im,jp,k,B2) and MACP0A1(prim,im,jm,k,B2)
        // and other pair is MACP0A1(prim,ip,jp,k,B2) and MACP0A1(prim,ip,jm,k,B2) (must treat as pair since iterating and haven't yet updated polar boundary values

        // then on lower \theta
        ip=ip1;
        im=i;
        jp=j; // otherwise similar to above section
        jm=jm1; // otherwise similar to above section

        get_geometry(im, jp, k, CENT, &geommp); // local point
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm); // anti-sym to local point

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jp,k,B2) = (-geompp.g*MACP0A1(prim,ip,jp,k,B2)*2.0 - 2.0*dx[2]*divbr)/(2.0*geommp.g);
      }
      else if(startpos[2]+j == 0 && FLIPGDETAXIS){
      }
      else if(startpos[2]+j<=totalsize[2]/2-1){
        // then on lower \theta
        ip=ip1;
        im=i;
        jp=j; // otherwise similar to above section
        jm=jm1; // otherwise similar to above section

        get_geometry(im, jp, k, CENT, &geommp); // local point
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jp,k,B2) = (-geompp.g*MACP0A1(prim,ip,jp,k,B2)+geompm.g*MACP0A1(prim,ip,jm,k,B2)+geommm.g*MACP0A1(prim,im,jm,k,B2) - 2.0*dx[2]*divbr)/(geommp.g);
      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  

  k=0;
  for(i=-1;i>=-N1BND;i--) for(j=OUTM2;j>=0;j--){ // then going from up to down

      if(startpos[2]+j<=totalsize[2]/2-1){
      }
      else if(startpos[2]+j == totalsize[2]-1 && FLIPGDETAXIS){
      }
      else if(startpos[2]+j == totalsize[2]-1 && FLIPGDETAXIS==0){
        // pairs are MACP0A1(prim,im,jm,k,B2) and MACP0A1(prim,im,jp,k,B2)
        // then on upper \theta
        ip=ip1;
        im=i;
        jp=jp1;
        jm=j;

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm); // local point

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jm,k,B2) = (-geompm.g*MACP0A1(prim,ip,jm,k,B2)*2.0 + 2.0*dx[2]*divbr)/(2.0*geommm.g);
      }
      else  if(startpos[2]+j>=totalsize[2]/2){
        // then on upper \theta
        ip=ip1;
        im=i;
        jp=jp1;
        jm=j;

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm); // local point

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jm,k,B2) = (geommp.g*MACP0A1(prim,im,jp,k,B2)-geompm.g*MACP0A1(prim,ip,jm,k,B2)+geompp.g*MACP0A1(prim,ip,jp,k,B2) + 2.0*dx[2]*divbr)/(geommm.g);
      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  



}


// assume equatorial value to be as outflowed and then clean rest of B^2 to satisfy FLUXCT version of divb=0
// too much accumulated error near poles
void clean_btheta_x1inner_old(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE divbr;
  struct of_geom geompp,geompm,geommp,geommm;
  int ip,im,jp,jm;

  k=0;
  for(i=-N1BND;i<0;i++) for(j=0;j<=OUTM2;j++){
      //for(i=-N1BND;i<0;i++) for(j=-1;j<=OUTM2+1;j++){

      if(startpos[2]+j == totalsize[2]/2 ) {
        // then leave as outflow but symmetrize
        MACP0A1(prim,i,j,k,B2) = 0.5*(MACP0A1(prim,i,j-1,k,B2)+MACP0A1(prim,i,j,k,B2));
      }
      else if(startpos[2]+j == totalsize[2]/2-1){
        // then leave as outflow but symmetrize
        MACP0A1(prim,i,j,k,B2) = 0.5*(MACP0A1(prim,i,j,k,B2)+MACP0A1(prim,i,j+1,k,B2));
      }
      else if(startpos[2]+j<totalsize[2]/2-1){
      }
      else  if(startpos[2]+j>totalsize[2]/2){
      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }
    }



  k=0;
  for(i=-1;i>=-N1BND;i--) for(j=0;j<=OUTM2;j++){ // going from down to up
      //for(i=-1;i>=-N1BND;i--) for(j=-1;j<=OUTM2+1;j++){ // going from down to up

      if(startpos[2]+j == totalsize[2]/2 ) {
      }
      else if(startpos[2]+j == totalsize[2]/2-1){
      }
      else if(startpos[2]+j<totalsize[2]/2-1){
      }
      else  if(startpos[2]+j>totalsize[2]/2){
        // then on upper \theta
        ip=ip1;
        im=i;
        jp=j; // otherwise similar to above section
        jm=jm1; // otherwise similar to above section

        get_geometry(im, jp, k, CENT, &geommp); // local point
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jp,k,B2) = (-geompp.g*MACP0A1(prim,ip,jp,k,B2)+geompm.g*MACP0A1(prim,ip,jm,k,B2)+geommm.g*MACP0A1(prim,im,jm,k,B2) - 2.0*dx[2]*divbr)/(geommp.g);
      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  

  k=0;
  for(i=-1;i>=-N1BND;i--) for(j=OUTM2;j>=0;j--){ // then going from up to down
      //for(i=-1;i>=-N1BND;i--) for(j=OUTM2+1;j>=-1;j--){ // then going from up to down

      if(startpos[2]+j == totalsize[2]/2 ) {
      }
      else if(startpos[2]+j == totalsize[2]/2-1){
      }
      else if(startpos[2]+j<totalsize[2]/2-1){
        // then on lower \theta
        ip=ip1;
        im=i;
        jp=jp1;
        jm=j;

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm); // local point

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,im,jm,k,B2) = (geommp.g*MACP0A1(prim,im,jp,k,B2)-geompm.g*MACP0A1(prim,ip,jm,k,B2)+geompp.g*MACP0A1(prim,ip,jp,k,B2) + 2.0*dx[2]*divbr)/(geommm.g);
      }
      else  if(startpos[2]+j>totalsize[2]/2){
      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  



}




void clean_btheta_x1outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE divbr;
  struct of_geom geompp,geompm,geommp,geommm;
  int ip,im,jp,jm;



  k=0;
  for(i=N1;i<=N1+N1BND-1;i++) for(j=0;j<=OUTM2;j++){ // going from down to up

      if(startpos[2]+j>=totalsize[2]/2){
      }
      else if(startpos[2]+j == 0 && FLIPGDETAXIS==0){
        // pairs are MACP0A1(prim,ip,jp,k,B2) and MACP0A1(prim,ip,jm,k,B2)
        // then on lower \theta
        ip=i;
        im=im1;
        jp=j; // otherwise similar to above section
        jm=jm1; // otherwise similar to above section

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp); // local point
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,ip,jp,k,B2) = (-geommp.g*MACP0A1(prim,im,jp,k,B2)*2.0 - 2.0*dx[2]*divbr)/(2.0*geompp.g);
      }
      else if(startpos[2]+j == 0 && FLIPGDETAXIS){
      }
      else if(startpos[2]+j<=totalsize[2]/2-1){
        // then on lower \theta
        ip=i;
        im=im1;
        jp=j; // otherwise similar to above section
        jm=jm1; // otherwise similar to above section

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp); // local point
        get_geometry(ip, jm, k, CENT, &geompm);
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,ip,jp,k,B2) = (-geommp.g*MACP0A1(prim,im,jp,k,B2)+geompm.g*MACP0A1(prim,ip,jm,k,B2)+geommm.g*MACP0A1(prim,im,jm,k,B2) - 2.0*dx[2]*divbr)/(geompp.g);

      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  

  k=0;
  for(i=N1;i<=N1+N1BND-1;i++) for(j=OUTM2;j>=0;j--){ // then going from up to down

      if(startpos[2]+j<=totalsize[2]/2-1){
      }
      else if(startpos[2]+j ==totalsize[2]-1 && FLIPGDETAXIS){
      }
      else if(startpos[2]+j ==totalsize[2]-1 && FLIPGDETAXIS==0){
        // pairs are MACP0A1(prim,ip,jm,k,B2) and MACP0A1(prim,ip,jp,k,B2)
        // then on upper \theta
        ip=i;
        im=im1;
        jp=jp1;
        jm=j;

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm); // local point
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,ip,jm,k,B2) = (-geommm.g*MACP0A1(prim,im,jm,k,B2)*2.0 + 2.0*dx[2]*divbr)/(2.0*geompm.g);
      }
      else  if(startpos[2]+j>=totalsize[2]/2) {
        // then on upper \theta
        ip=i;
        im=im1;
        jp=jp1;
        jm=j;

        get_geometry(im, jp, k, CENT, &geommp);
        get_geometry(ip, jp, k, CENT, &geompp);
        get_geometry(ip, jm, k, CENT, &geompm); // local point
        get_geometry(im, jm, k, CENT, &geommm);

        divbr = ( (geompp.g * MACP0A1(prim,ip,jp,k,B1)+ geompm.g*MACP0A1(prim,ip,jm,k,B1)) - (geommp.g*MACP0A1(prim,im,jp,k,B1)+geommm.g*MACP0A1(prim,im,jm,k,B1)))/(2.0*dx[1]) ;
        MACP0A1(prim,ip,jm,k,B2) = (geompp.g*MACP0A1(prim,ip,jp,k,B2)+geommp.g*MACP0A1(prim,im,jp,k,B2)-geommm.g*MACP0A1(prim,im,jm,k,B2) + 2.0*dx[2]*divbr)/(geompm.g);


      }
      else{
        dualfprintf(fail_file,"problem 117\n");
        myexit(117);
      }

    }  



}



void compute_aphi_fromoutflow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*aphifromoutflow)[NSTORE2][NSTORE3])
{
  int i,j,k;
  int ii;
  struct of_geom geom,rgeom;
  struct of_geom geom1,geom2;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE rdxdxp[NDIM][NDIM];
  FTYPE V[NDIM],X[NDIM];
  FTYPE rV[NDIM],rX[NDIM];
  FTYPE Bcon[NDIM];
  FTYPE B2face2;
  extern FTYPE (*brface)[NSTORE2][NSTORE3];
  extern FTYPE (*aphicorn)[NSTORE2][NSTORE3];

  ////////////////////////////////////////////////////////////////////////////
  //
  // Compute A_\phi using outflowed B^\theta and analytic B^r
  //
  // Need A_\phi at corners from j = 0 .. N2 and and i=-N1BND..0 and k=0..N3
  //
  ////////////////////////////////////////////////////////////////////////////
      
  //  dualfprintf(fail_file,"asdf=%d\n",aphicorn);

  for(i=-N1BND;i<=0;i++) for(j=0;j<=OUTM2;j++) for(k=0;k<=OUTM3;k++){


        MACP0A0(aphifromoutflow,i,j,k)=MAC(aphicorn,i,j,k);
        for(ii=-N1NOT1;ii>=i;ii--){  // from -1 down to i

          //      dualfprintf(fail_file,"i=%d j=%d ii=%d\n",i,j,ii);
          get_geometry(ii, j, k, FACE2, &geom);
          //      coord(ii, j, k, FACE2, X);
          //      bl_coord( X, V );
          //      dxdxprim(X, V, dxdxp);

          // FLUXCT Toth centered method
          get_geometry(ii, j  , k, CENT, &geom1);
          get_geometry(ii, j-1, k, CENT, &geom2);
          B2face2=0.5*(geom1.g*MACP0A1(prim,ii,j,k,B2)+geom2.g*MACP0A1(prim,ii,j-1,k,B2))/(geom.g);

          //      B2face2=0.5*(MACP0A1(prim,ii,j,k,B2)+MACP0A1(prim,ii,j-1,k,B2));

          MACP0A0(aphifromoutflow,i,j,k)-= (geom.g) * B2face2 * dx[1]; // centered integral
          // gives aphi from corner to corner
        }

        //    dualfprintf(fail_file,"aphifromoutflow[%d][%d][%d]=%21.15g %21.15g %21.15g\n",i,j,k,MACP0A0(aphifromoutflow,i,j,k),B2face2,MAC(aphicorn,i,j,k));
   
      }
}










void compute_btheta_fromaphi(FTYPE (*aphifromoutflow)[NSTORE2][NSTORE3], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  struct of_geom geom,rgeom;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE rdxdxp[NDIM][NDIM];
  FTYPE V[NDIM],X[NDIM];
  FTYPE rV[NDIM],rX[NDIM];
  FTYPE Bcon[NDIM];
  extern FTYPE (*brface)[NSTORE2][NSTORE3];
  extern FTYPE (*aphicorn)[NSTORE2][NSTORE3];


  ////////////////////////////////////////////////////////////////////////////
  //
  ///////////////////////////////////////////////   REDO B^\theta to be consistent with divb=0 in FLUXCT form
  //
  ////////////////////////////////////////////////////////////////////////////

  LOOPN2 LOOPN3{
    ri=0;
    rj=j;
    rk=k;
    LOOPBOUND1IN{

      //      dualfprintf(fail_file,"bgod %d %d\n",i,j);

      get_geometry(i, j, k, CENT, &geom);
      //      coord(i, j, k, CENT, X);
      //      bl_coord( X, V );
      //      dxdxprim(X, V, dxdxp);

      MACP0A1(prim,i,j,k,B2)=0.0;
      //      MACP0A1(prim,i,j,k,B2)  = +(AVGCORN_2(A[1],i,j,kp1mac(k))-AVGCORN_2(A[1],i,j,k))/(geom.g*dx[3]); // not 3D
      // FLUXCT Toth centered method
      MACP0A1(prim,i,j,k,B2) += -(AVGCORN_2(aphifromoutflow,ip1mac(i),j,k)-AVGCORN_2(aphifromoutflow,i,j,k))/(geom.g*dx[1]);
      //      MACP0A1(prim,i,j,k,B2) += -(AVGCORN_2(aphicorn,ip1mac(i),j,k)-AVGCORN_2(aphicorn,i,j,k))/(geom.g*dx[1]);
      //      dualfprintf(fail_file,"agod %d %d\n",i,j);


      MACP0A1(prim,i,j,k,UU) = MACP0A0(aphifromoutflow,i,j,k);

    }
  }
}



int bound_vel_from_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  struct of_geom geom;
  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);


  ////////////////////////////////////////////////////////////////////////////
  //
  ///////////////////////////////////////////////   SET VELOCITY
  //
  ////////////////////////////////////////////////////////////////////////////

  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==NSSURFACE)){


      LOOPN2 LOOPN3{
        LOOPBOUND1IN{

          get_geometry(i, j, k, CENT, &geom);
          //   coord(i, j, k, CENT, X);
          //bl_coord( X, V );
          //dxdxprim(X, V, dxdxp);

   
          set_vel_stataxi(&geom, GLOBALMACP1A0(pother,OMEGAFCENT,i,j,k), GLOBALMACP1A0(pother,VPARFACE1,i,j,k), MAC(prim,i,j,k));


        }
      }
    }
  }

  return(0);

}




int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr)
{
  int i,j,k;
  int pl, pl2;
  FTYPE vcon[NDIM]; // coordinate basis vcon
  extern FTYPE Omegastar;
  extern FTYPE Bpole;
  extern FTYPE Vpar;
  FTYPE prnew[NPR];
  FTYPE Bcon[NDIM];
  extern int OBtopr_general(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr_general2(FTYPE omegaf, FTYPE vr, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);



  i=geom->i;
  j=geom->j;
  k=geom->k;



  Bcon[0]=0;
  Bcon[1]=pr[B1];
  Bcon[2]=pr[B2];
  Bcon[3]=pr[B3];


#if(0) // what would be done in force-free with no parallel velocity

  //   dualfprintf(fail_file,"Omegastar=%21.15g dxdxp[3][3]=%21.15g B1=%21.15g B2=%21.15g B3=%21.15g\n",Omegastar,dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);

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

  //   if(t>1.9 && t<2.1){
  //dualfprintf(fail_file,"t=%21.15g i=%d j=%d\n",t,i,j);
  //  dualfprintf(fail_file,"newus: %21.15g %21.15g %21.15g\n",prnew[U1],prnew[U2],prnew[U3]);
  //  dualfprintf(fail_file,"Omegastar'=%21.15g Bcon1=%21.15g Bcon2=%21.15g Bcon3=%21.15g\n",Omegastar/dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);
  // }

#endif

#if(0)
  // set       ucon[TT,etc.]
  vcon[RR]=0; // surface that completely dissipates normal direction momentum
  vcon[TH]=0; // "" for this component
  // below assumes no phi mixing with other directions in grid
  vcon[PH]=omegaf; // surface rotates with angular frequency Omegastar to observer at infinity

  // get in terms of primitive velocity
  MYFUN(vcon2pr(WHICHVEL,vcon,geom,pr),"bounds.ns.c:bound_prim_user()", "vcon2pr()", 1);
#endif

#if(1) 

  ///////////////////////////////
  //
  // new way to get velocity
  // surface rotates with angular frequency Omegastar to observer at infinity
  if(OBtopr_general2(omegaf,vpar,Bcon,geom,prnew)>=1){
    dualfprintf(fail_file, "OBtopr(bounds): space-like error in init_postfield()\n");
    dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
    failed=1;
    return(1);
  }
  // assign answer
  pr[U1]=prnew[U1];
  pr[U2]=prnew[U2];
  pr[U3]=prnew[U3];

#endif


  PLOOP(pliter,pl){
    if(!isfinite(pr[pl])){
      dualfprintf(fail_file,"i=%d j=%d steppart=%d nstep=%ld\n",i,j,steppart,nstep);
      PLOOP(pliter,pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,pr[pl2]);
      myexit(2525);
    }
  }



  return(0);

}





void x2_inner_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;



  ///////////////////////////
  //
  // X2 inner POLARAXIS
  //
  ///////////////////////////


  /* inner polar BC (preserves u^t rho and u) */
  if (mycpupos[2] == 0) {
    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
      LOOPF1 LOOPF3{
        ri=i;
        rj=0;
        rk=k;
        LOOPBOUND2IN{
          PBOUNDLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j-1),rk,pl);
          // dualfprintf(fail_file,"i=%d j=%d ri=%d rj=%d  :: rj+(rj-j-1)=%d\n",i,j,ri,rj,rj+(rj-j-1));
        }
      }
    }

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
        LOOPBOUND2IN {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= -1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= -1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= -1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= -1.;
          }
        }
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPF3 {
        for (j = 0-POLEDEATH; j < 0+POLEDEATH; j++) {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            MACP0A1(prim,i,j,k,U2) *= 0;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= 0.;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
        }
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==0


}



void x2_outer_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;

  ///////////////////////////
  //
  // X2 outer POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == ncpux2-1) {
    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
      LOOPF1 LOOPN3{
        ri=i;
        rj=N2-1;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j+1),rk,pl);
      }
    }

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
        LOOPBOUND2OUT {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= -1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= -1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= -1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= -1.;
          }
        }
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPF3 {
        for (j = N2-1+1-POLEDEATH; j <= N2-1+POLEDEATH; j++) {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            MACP0A1(prim,i,j,k,U2) *= 0;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= 0.;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
        }
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==ncpux2-1

}





int bound_x1_outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;
  struct of_geom geom,rgeom;
  FTYPE prescale[NPR];
#if(WHICHVEL==VEL3)
  int failreturn;
#endif

  ///////////////////////////
  //
  // X1 outer OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////


  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
      /* outer r BC: outflow */

      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#elif(EXTRAP==1)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT{
          for(pl=RHO;pl<=UU;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          }
          pl=U1; // treat U1 as special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - 2*(i-ri)*dx[1]) ;
          for(pl=U2;pl<=U3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
          }
          pl=B1; // treat B1 special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          for(pl=B2;pl<=B3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
          }
        }
#elif(EXTRAP==2)
        ri=N1-1;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,MAC(prim,ri,rj,rk),&rgeom,prescale);
        LOOPBOUND1OUT{
          // set guess
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,rk,pl);
          get_geometry(i, j, k, CENT, &geom);
          rescale(-1,1,MAC(prim,i,j,k),&geom,prescale);
        }
#endif

        LOOPBOUND1OUT{
#if(WHICHVEL==VEL4)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_4vel(1,MAC(prim,i,j,k),&geom,0) ;
#elif(WHICHVEL==VEL3)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_3vel(1,MAC(prim,i,j,k),&geom,0) ;
          // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
          if(jonchecks){
            //fixup1zone(MAC(prim,i,j,k),&geom,0);
            failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),&geom,-3);
            if(failreturn){
              dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
              if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
            }
          }
#endif
#elif(WHICHVEL==VELREL4)
          get_geometry(i,j,k,CENT,&geom) ;
          inflow_check_rel4vel(1,MAC(prim,i,j,k),&geom,0) ;
          if(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,MAC(prim,i,j,k),&geom, 0)>=1)
            FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif 
        }
      }// end 2 3
    }// end if correct bound type
  }// end if mycpu is correct


  return(0);
}





void x3_periodic(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;

  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another

      LOOPF1 LOOPF2{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=j;
        rk=N3;
        LOOPBOUND3IN PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+k,pl);

        // copy from lower side to upper boundary zones
        ri=i;
        rj=j;
        rk=0;
        LOOPBOUND3OUT PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(k-N3),pl);
      }
    }
  }
}








// see interpline.c
int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])
{
  int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM]);

  if(doinverse==0){
    flip_y(iterglobal, recontype, bs, be, yin);
  }
  else{
    flip_y(iterglobal, recontype, bs, be, yin);
    flip_y(iterglobal, recontype, bs, be, yout);
  }

  return(0);

}


#include "reconstructeno.h"

int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM])
{
  int pl,myi;


#if( WENO_DIR_FLIP_CONS_SIGN_DN )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_DN && (recontype == CVT_C2A || recontype == CVT_A2C) && mycpupos[iterglobal] == 0 ) { 
    PLOOP(pliter,pl) 
      for( myi = bs; myi < 0; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif
 
#if( WENO_DIR_FLIP_CONS_SIGN_UP )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_UP && (recontype == CVT_C2A || recontype == CVT_A2C)  && mycpupos[iterglobal] == numbercpu[iterglobal] - 1 ) { 
    PLOOP(pliter,pl) 
      for( myi = N1*(iterglobal==1) + N2*(iterglobal==2) + N3*(iterglobal==3); myi <= be; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif


  return(0);

}


///Called after the MPI boundary routines
int bound_prim_user_after_mpi(int boundstage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  return(0);
}


void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE *p2interp_l, FTYPE *p2interp_r )
{

}



void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, 
                FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE *p2interp_l, FTYPE *p2interp_r )
{
}
