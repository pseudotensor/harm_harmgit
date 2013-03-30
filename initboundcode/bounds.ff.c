
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

extern int DISKDIR;

#if(0)
//int setnsflux[NPR]={SETNSFLUXRHO,SETNSFLUXUU,SETNSFLUXV1,SETNSFLUXV2,SETNSFLUXV3,SETNSFLUXB1,SETNSFLUXB2,SETNSFLUXB3};
int setnsflux[NPR]  ={1,1,1,1,1,1,1,1};
int setdiskflux[NPR]={1,1,1,1,1,1,1,1};
int setpolarflux[NPR]={1,1,1,1,1,1,1,1};

//whether to avoid using boundary values in the reconstruction; reconavoidbcxx[0] controls dir == 1
int a_reconavoidbcdn[NDIM-1] = { 1, 1, 0 };  //avoid using ghost zones for reconstructions on the star surface (X1DN) 
                                             //and on the disk surface (X2DN); !!! assumes DISKDIR = X1DN (SASMARK)

int a_reconavoidbcup[NDIM-1] = { 0, 1, 0 };  //do not avoid using ghost zones for reconstructions 
#endif

#if(1)
//int setnsflux[NPR]={SETNSFLUXRHO,SETNSFLUXUU,SETNSFLUXV1,SETNSFLUXV2,SETNSFLUXV3,SETNSFLUXB1,SETNSFLUXB2,SETNSFLUXB3};
int setnsflux[NPR]  ={0,0,0,0,0,0,0,0};
int setdiskflux[NPR]={0,0,0,0,0,0,0,0};
int setpolarflux[NPR]={1,1,1,1,1,1,1,1};

//whether to avoid using boundary values in the reconstruction; reconavoidbcxx[0] controls dir == 1
int a_reconavoidbcdn[NDIM-1] = { 0, 0, 0 };  //avoid using ghost zones for reconstructions on the star surface (X1DN) 
                                             //and on the disk surface (X2DN); !!! assumes DISKDIR = X1DN (SASMARK)

int a_reconavoidbcup[NDIM-1] = { 0, 1, 0 };  //do not avoid using ghost zones for reconstructions 
#endif



//shift the pointer such that reconavoidbc[dir] can be used
int *reconavoidbcdn = &a_reconavoidbcdn[-1];
int *reconavoidbcup = &a_reconavoidbcup[-1];

#define NOSETFACE 0
#define BCSETFACE 1
#define RECONSETFACE 2

//whether to set face values in boundary condition routines or by reconstruction
#define WHICHSETFACE RECONSETFACE
//#define WHICHSETFACE BCSETFACE
//#define WHICHSETFACE NOSETFACE


//int setnsflux[NPR]={1,1,  1,1,1,   1,1,1};
//int setnsflux[NPR]={0,0,0,1,0,0,0,0}; // V2 problem
//int setnsflux[NPR]={1,1,1,0,0,0,0,0}; // V1 problem
//int setnsflux[NPR]={0,0,0,0,0,1,0,0}; // B1 problem
//int setnsflux[NPR]={1,1,1,1,1,0,0,0};

// whether to
// 2: set via stationarity condition (\rho u^p / B^p = const[face value])
// 1: set ghost zones analytically
#define SETGHOSTRHO 1


// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
//#define POLEDEATH N2BND
#define POLEDEATH 2
//#define POLEDEATH 0
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
  void x1_star(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x1_outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x2_inner_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x2_outer_polar(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x2_disk(int dir, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  void x3_periodic(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  static int firsttime=1;



  //////////////////////////////////////////////////////////////////////////
  //
  /// FILL IN GHOST ZONES with DUMMY VALUES
  //
  //////////////////////////////////////////////////////////////////////////

  // polar axis bound so field is bounded fully
  x2_inner_polar(prim);
  x2_outer_polar(prim);

  // basic outflow
  basic_outflow(prim);

  // polar axis bound so field is bounded fully
  x2_inner_polar(prim);
  x2_outer_polar(prim);



  //////////////////////////////////////////////////////////////////////////
  //
  /// FILL IN GHOST ZONES with DUMMY VALUES
  //
  //////////////////////////////////////////////////////////////////////////


  //outflow into the star (lower-r)
  x1_star(prim);
  x1_outer(prim);

  // polar axis bound and disk bound so field is bounded fully

  if(DISKDIR==X2DN){
    if(firsttime) trifprintf("x2_outer\n");
    x2_outer_polar(prim);
  }
  else if(DISKDIR==X2UP){
    if(firsttime) trifprintf("x2_inner\n");
    x2_inner_polar(prim);
  }
  if(firsttime) trifprintf("x2_disk\n");
  x2_disk(DISKDIR,prim);

    
  // x3
  x3_periodic(prim);


  firsttime=0;
  return (0);
}


void x1_star(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);
  void bound_face1dn(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface);


  int ri, rj, rk;
  int fi, fj, fk;
  int i, j, k;
  int pl,pliter;
  struct of_geom geom, rgeom;
  FTYPE X[NDIM], V[NDIM], rX[NDIM], rV[NDIM];
  FTYPE omegaf;
  FTYPE prface[NPR];




  if( BCtype[X1DN] != NSSURFACE ) {
    dualfprintf( fail_file, "Unknown BCtype[X1DN] = %d\n", BCtype[X1DN] );
    myexit(1);
  }

  if( mycpupos[1] == 0 ){ 


    ri=0;
    fi=0;

    LOOPF2 LOOPF3 {
      rj=j;
      rk=k;

      fj=j;
      fk=k;

      // set face value
#if(WHICHSETFACE == BCSETFACE)
      bound_face1dn(fi,fj,fk,prim,prface);
#endif

      get_geometry( ri, rj, rk, CENT, &rgeom );
      coord( ri, rj, rk, CENT, rX );
      bl_coord( rX, rV );

      // loop over ghost zones at inner edge (star)
      LOOPBOUND1IN{

        get_geometry( i, j, k, CENT, &geom );
        coord( i, j, k, CENT, X );
        bl_coord( X, V );


        // outflow densities
        pl = RHO; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
        pl = UU; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

        //reset the normal magnetic field to the initial value
        pl = B1; MACP0A1(prim,i,j,k,pl) = MACP0A1(panalytic,i,j,k,pl);


        if(DISKDIR==X2DN && (startpos[2]+j==totalsize[2]-1 || startpos[2]+j==totalsize[2]-2 )){
          pl = B2; MACP0A1(prim,i,j,k,pl) = 0.0;
        }
        else if(DISKDIR==X2UP && (startpos[2]+j==0 || startpos[2]+j==1 )){
          pl = B2; MACP0A1(prim,i,j,k,pl) = 0.0;
        }
        else{
          // outflow \detg B^2
          pl = B2; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.g)/(geom.g);
        }

        // outflow B_\phi
        pl = B3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.gcov[GIND(3,3)])/(geom.gcov[GIND(3,3)]);

        //reset velocities to stationary ones
        set_vel_stataxi(&geom,GLOBALMACP1A0(pother,OMEGAFCENT,i,j,k),GLOBALMACP1A0(pother,VPARCENT,i,j,k),MAC(prim,i,j,k));


      }
    }
  }
}

//outflow some quantities from the active grid to the interface and set the quantities to be fixed on the interface
//saves the face values in GLOBALMACP1A0(pother,xxxFACE1,i,j,k)
void bound_face1dn(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface)
{
  void outflow_face1dn(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface);
  void set_face1dn(int i, int j, int k, FTYPE *prface);

  int pl,pliter;

  //outflow from active grid to face
  outflow_face1dn( i, j, k, prim, prface );

  //set quantities that are to be fixed on the face
  set_face1dn( i, j, k, prface);

  //save face values in a global variable
  PLOOP(pliter,pl) GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k) = prface[pl];
}


// set prface for fluxes at boundary and interpolation through boundary value
// Analytical set of B1, omegaf, and vpar
//assumes that prface is populated with outflown values for pl = 0 (RHO), 1 (UU), 6 (B2), 7 (B3)
void set_face1dn(int i, int j, int k, FTYPE *prface)
{
  int pl,pliter;
  FTYPE Bcon[NDIM];
  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);
  FTYPE Xf[NDIM],Vf[NDIM];
  struct of_geom geomf;
  FTYPE dxdxpf[NDIM][NDIM];


  if(is_on_surface(X1DN,i,j,k,FACE1)){

    // face geometry
    get_geometry(i, j, k, FACE1, &geomf);
    coord(i, j, k, FACE1, Xf);
    bl_coord( Xf, Vf );
    dxdxprim(Xf, Vf, dxdxpf);


    prface[RHO] = GLOBALMACP1A0(pother,RHOFACE1,i,j,k);
    prface[UU] = GLOBALMACP1A0(pother,UUFACE1,i,j,k);
    prface[B1] = GLOBALMACP1A0(pother,B1FACE1,i,j,k);

    //reset velocities to stationary ones (VPAR analytical)
    set_vel_stataxi(&geomf,GLOBALMACP1A0(pother,OMEGAFFACE1,i,j,k),GLOBALMACP1A0(pother,VPARFACE1,i,j,k),prface);
  }
}


void outflow_face1dn(int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface)
{
  int pl,pliter;
  FTYPE Bcon[NDIM];
  FTYPE Xc[NDIM],Vc[NDIM];
  struct of_geom geomc;
  FTYPE dxdxpc[NDIM][NDIM];
  //
  FTYPE Xf[NDIM],Vf[NDIM];
  struct of_geom geomf;
  FTYPE dxdxpf[NDIM][NDIM];
  int ri,rj,rk;
  FTYPE *rprim;


  if(is_on_surface(X1DN,i,j,k,FACE1)){

    // face geometry
    get_geometry(i, j, k, FACE1, &geomf);
    coord(i, j, k, FACE1, Xf);
    bl_coord( Xf, Vf );
    dxdxprim(Xf, Vf, dxdxpf);


    // center geometry
    // for X1DN reference is i,j,k
    ri=i;
    rj=j;
    rk=k;

    rprim=MAC(prim,ri,rj,rk);
    get_geometry(ri, rj, rk, CENT, &geomc);
    coord(ri, rj, rk, CENT, Xc);
    bl_coord( Xc, Vc );
    dxdxprim(Xc, Vc, dxdxpc);

    //first, outflow everything
    PLOOP(pliter,pl) prface[pl] = rprim[pl];

    ///////////////////////////////////
    //
    /////////// B2
    //    prface[B2] = rprim[B2];
    // outflow using gdet
    if(DISKDIR==X2DN && (startpos[2]+j==totalsize[2]-1 || startpos[2]+j==totalsize[2]-2 )){
      pl = B2; prface[pl] = 0.0;
    }
    else if(DISKDIR==X2UP && (startpos[2]+j==0 || startpos[2]+j==1 )){
      pl = B2; prface[pl] = 0.0;
    }
    else{
      // outflow \detg B^2
      pl = B2; prface[pl] = rprim[pl]*(geomc.g)/(geomf.g);
    }
    
    ///////////////////////////////////
    //
    /////////// B3

    // outflow as copy
    //    if(sign(rprim[B3])==sign(-prface[B1]*GLOBALMACP1A0(pother,OMEGAFFACE1,i,j,k))){
    //      // then normal situation so can extrapolate 
    //      prface[B3] = rprim[B3]*pow(Vc[1]/Vf[1],3);
    //    }
    //    else{ // then assume out of equilibrium so don't extrapolate
    //      prface[B3] = rprim[B3];
    //    }

    // B_\phi is constant along field lines, so this is closer
    prface[B3] = rprim[B3]*(geomc.gcov[GIND(3,3)])/(geomf.gcov[GIND(3,3)]);
  }
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






void set_plpr(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r)
{
  void set_face1dn(int i, int j, int k, FTYPE *prface);
  void set_face2(int dir, int i, int j, int k, FTYPE *prface);

  int pl,pliter;
  FTYPE prface[NPR];


  PLOOP(pliter,pl) prface[pl] = 0.5 * ( p_l[pl] + p_r[pl] );  //initialize prface with boundary values

  // need a new function that sets quantities at the face1 to be used when setting p_l and p_r
  // need to set v^i at the FACE1
  if(is_on_surface(X1DN,i,j,k,FACE1) && dir == 1){
    // then at stellar surface for this flux
    //    dualfprintf(fail_file,"nstep=%ld steppart=%d j=%d\n",nstep,steppart,j);

    PLOOP(pliter,pl){
#if( WHICHSETFACE == RECONSETFACE ) 
      set_face1dn( i, j, k, prface );
#elif( WHICHSETFACE == BCSETFACE )
      PLOOP(pliter,pl) prface[pl] = GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k); // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
#endif

      if(setnsflux[pl]){
        if( 0 && pl == U1 ) {
          dualfprintf(fail_file,"n = %ld sp = %d j = %d pl=%d p_l=%21.15g p_r=%21.15g pother=%21.15g dp/p = %21.15g, dplpr/p = %21.15g\n",
                      nstep,steppart,startpos[2]+j,pl,p_l[pl],p_r[pl],GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k), 
                      fabs((p_l[pl]-GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k))/GLOBALMACP1A0(pother,RHOFACE1+pl,i,j,k)),
                      2*fabs((p_l[pl]-p_r[pl])/(p_r[pl]+p_l[pl])) 
                      ); 
        }

        p_l[pl]=p_r[pl]=prface[pl];

        //p_l[pl]=p_r[pl]; // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
        //p_r[pl]=p_l[pl];
      }

    }
  }
  else if(is_on_surface(DISKDIR,i,j,k,FACE2) && dir == 2){
    //  else if( startpos[2] + j == totalsize[2] && dir == 2 ){

#if( WHICHSETFACE == RECONSETFACE ) 
    set_face2( DISKDIR, i, j, k, prface );
#elif( WHICHSETFACE == BCSETFACE )
    PLOOP(pliter,pl) prface[pl] = GLOBALMACP1A0(pother,RHOFACE2+pl,i,j,k);  // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
#endif

    PLOOP(pliter,pl) {
      if( setdiskflux[pl] ) {
        p_l[pl]=p_r[pl]=prface[pl];
      }
    }
  }
  else if(BCtype[X2UP] == POLARAXIS && is_on_surface(X2UP,i,j,k,FACE2) && dir == 2){
    //  else if( startpos[2] + j == totalsize[2] && dir == 2 ){

#if( WHICHSETFACE == RECONSETFACE ) 
    set_face2( X2UP, i, j, k, prface );
#elif( WHICHSETFACE == BCSETFACE )
    set_face2( X2UP, i, j, k, prface );
    //PLOOP(pliter,pl) prface[pl] = GLOBALMACP1A0(pother,RHOFACE2+pl,i,j,k);  // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
#endif

    PLOOP(pliter,pl) {
      if( setpolarflux[pl] ) {
        p_l[pl]=p_r[pl]=prface[pl];
      }
    }
  }
  else if(BCtype[X2DN] == POLARAXIS && is_on_surface(X2DN,i,j,k,FACE2) && dir == 2){
    //  else if( startpos[2] + j == totalsize[2] && dir == 2 ){

#if( WHICHSETFACE == RECONSETFACE ) 
    set_face2( X2DN, i, j, k, prface );
#elif( WHICHSETFACE == BCSETFACE )
    set_face2( X2DN, i, j, k, prface );
    //PLOOP(pliter,pl) prface[pl] = GLOBALMACP1A0(pother,RHOFACE2+pl,i,j,k);  // pl=0..NPR-1  (assumes ordering of pother[RHOFACE1->B3FACE1] as for standard primitives)
#endif

    PLOOP(pliter,pl) {
      if( setpolarflux[pl] ) {
        p_l[pl]=p_r[pl]=prface[pl];
      }
    }
  }
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
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
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




///////////////////////////
//
// X2 DISK
//
///////////////////////////
void x2_disk(int dir, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  // time while static boundary conditions
  //#define TSTATIC (15.0)

  // whether to use new method to limit v^i to be time-like
#define LIMIT3VELBOUND 1


  // 0 = 0th order prims
  // 1 = linear prims
  // 2 = 0th order ramesh self-sim
  // 3 = linear ramesh self-sim
#define BOUNDARYINTERPTYPE 3

#define TSTATIC (0.0)



  int set_vel_stataxi(struct of_geom *geom, FTYPE omegaf, FTYPE vpar, FTYPE *pr);
  void bound_face2disk(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface);
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;
  FTYPE X[NDIM], V[NDIM];
  struct of_geom rgeom;
  FTYPE rX[NDIM], rV[NDIM];
  struct of_geom geom;
  FTYPE omegaf;
  FTYPE prface[NPR];
  int fi,fj,fk;
  int dobounds;
  int sj,ej;

  if( BCtype[dir] != DISKSURFACE ) {
    dualfprintf( fail_file, "Unknown BCtype[%d] = %d\n", dir,BCtype[dir] );
    myexit(1);
  }

  
  dobounds=0;

  if(dir==X2UP && mycpupos[2] == ncpux2-1){
    sj=OUTBOUNDLO2;
    ej=OUTBOUNDHI2;
    rj = N2 - 1;
    fj = N2;
    dobounds=1;
    
  }
  else if(dir==X2DN && mycpupos[2] == 0){
    sj=INBOUNDLO2;
    ej=INBOUNDHI2;
    rj = 0;
    fj = 0;
    dobounds=1;
  }


  if(dobounds) {
    
    LOOPF1 LOOPF3{ // loop over entire domain looking for equator

      ri = i;
      rk = k;
      
      get_geometry( ri, rj, rk, CENT, &rgeom );
      coord( ri, rj, rk, CENT, rX );
      bl_coord( rX, rV );

      // set face value
      fi=ri;
      fk=rk;

#if(WHICHSETFACE == BCSETFACE)
      bound_face2disk(dir, fi, fj, fk, prim, prface);
#endif

      // loop over ghost zones inside disk
      for(j=sj;j<=ej;j++){

        get_geometry( i, j, k, CENT, &geom );
        coord( i, j, k, CENT, X );
        bl_coord( X, V );


        // set densities analytically
        pl = RHO; MACP0A1(prim,i,j,k,pl) = MACP0A1(panalytic,i,j,k,pl);
        pl = UU; MACP0A1(prim,i,j,k,pl) = MACP0A1(panalytic,i,j,k,pl);

        // outflow \detg B^1
        pl = B1; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.g)/(geom.g);

        // set the normal magnetic field to the initial value
        pl = B2; MACP0A1(prim,i,j,k,pl) = MACP0A1(panalytic,i,j,k,pl);

        // outflow B_\phi
        pl = B3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.gcov[GIND(3,3)])/(geom.gcov[GIND(3,3)]);

        //reset velocities to stationary ones (VPAR analytical)
        set_vel_stataxi(&geom,GLOBALMACP1A0(pother,OMEGAFCENT,i,j,k),GLOBALMACP1A0(pother,VPARCENT,i,j,k),MAC(prim,i,j,k));

      }// end loop over ghost cells
    }// end loop over domain of x1 and x3
  }// end if mycpupos[2]==ncpux2-1
}



///////////////////////////
//
// set analytical values on disk surface
//
///////////////////////////
void outflow_face2disk(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface)
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
  int ri,rj,rk;



  if(is_on_surface(dir,i,j,k,FACE2)){


    // face geometry
    get_geometry(i, j, k, FACE2, &geomf);
    coord(i, j, k, FACE2, Xf);
    bl_coord( Xf, Vf );
    dxdxprim(Xf, Vf, dxdxpf);


    // center geometry
    ri=i;
    rj=j - (DIRSIGN(dir)==1);
    rk=k;
    // for X2DN, center reference cell is at ri,rj,rk
    get_geometry(ri, rj, rk, CENT, &geomc);
    coord(ri, rj, rk, CENT, Xc);
    bl_coord( Xc, Vc );
    dxdxprim(Xc, Vc, dxdxpc);

    //first, outflow everything
    PLOOP(pliter,pl) prface[pl] = MACP0A1(prim,ri,rj,rk,pl);

    //////
    //
    // Outflow tangential fields
    //
    //////

    // outflow B1
    prface[B1]=MACP0A1(prim,ri,rj,rk,B1)*(geomc.g)/(geomf.g);

    // outflow B_\phi
    prface[B3]=MACP0A1(prim,ri,rj,rk,B3)*fabs(geomc.gcov[GIND(3,3)])/(geomf.gcov[GIND(3,3)]);
  }// end if on disk surface

}


//sets face values on the disk surface whose location is given by dir = X2DN or X2UP.
//assumes that prface is populated with outflown values for pl = 0 (RHO), 1 (UU), 5 (B1), 7 (B3)
void set_face2(int dir, int i, int j, int k, FTYPE *prface)
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
  int ri,rj,rk;



  if(is_on_surface(dir,i,j,k,FACE2) && dir == DISKDIR){
    // face geometry
    get_geometry(i, j, k, FACE2, &geomf);
    coord(i, j, k, FACE2, Xf);
    bl_coord( Xf, Vf );
    dxdxprim(Xf, Vf, dxdxpf);


    // densities
    prface[RHO] = GLOBALMACP1A0(pother,RHOFACE2,i,j,k);
    prface[UU] = GLOBALMACP1A0(pother,UUFACE2,i,j,k);

    // set B2
    prface[B2]=GLOBALMACP1A0(pother,B2FACE2,i,j,k);

    //get velocity
    set_vel_stataxi(&geomf,GLOBALMACP1A0(pother,OMEGAFFACE2,i,j,k),GLOBALMACP1A0(pother,VPARFACE2,i,j,k),prface);

  }// end if on disk surface
  else if(is_on_surface(dir,i,j,k,FACE2) && dir != DISKDIR && BCtype[dir] == POLARAXIS ) { //polar axis
    //reset asymmetric quantities to zero at the axis
    prface[B2] = 0.;
    prface[U2] = 0.;
  }


}

//outflows, sets, and saves interface values at the disk surface whose location is given by dir = X2DN or X2UP
void bound_face2disk(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface)
{
  void outflow_face2disk(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *prface);
  void set_face2(int dir, int i, int j, int k, FTYPE *prface);

  int pl,pliter;

  outflow_face2disk( dir, i, j, k, prim, prface);
  set_face2( dir, i, j, k, prface);

  // now set pother so can be used by rest of routines in any other functions
  PLOOP(pliter,pl){
    GLOBALMACP1A0(pother,RHOFACE2+pl,i,j,k)=prface[pl];
  }
}

///////////////////////////
//
// X1 outer radial edge
//
///////////////////////////
void x1_outer(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ri,rj,rk;
  int pl,pliter;
  FTYPE X[NDIM], V[NDIM];
  struct of_geom rgeom;
  FTYPE rX[NDIM], rV[NDIM];
  struct of_geom geom;




  if( mycpupos[1] == ncpux1-1 ) {


    LOOPF2 LOOPF3{ // loop over entire domain looking for equator

      ri = N1 - 1;
      rj = j;
      rk = k;
      
      get_geometry( ri, rj, rk, CENT, &rgeom );
      coord( ri, rj, rk, CENT, rX );
      bl_coord( rX, rV );

      // loop over ghost zones at outer edge
      LOOPBOUND1OUT{

        get_geometry( i, j, k, CENT, &geom );
        coord( i, j, k, CENT, X );
        bl_coord( X, V );

        // outflow densities
        pl = RHO; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
        pl = UU; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

        // outflow B^1
        pl = B1; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

        // outflow \detg B^2
        pl = B2; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.g)/(geom.g);

        // outflow B_\phi
        pl = B3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*(rgeom.gcov[GIND(3,3)])/(geom.gcov[GIND(3,3)]);

        // outflow velocities (copy)
        pl = U1; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
        pl = U2; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
        pl = U3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

      }// end loop over ghost cells
    }// end loop over domain of x2 and x3
  }// end if mycpupos[1]==ncpux1-1
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
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
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
int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM])
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



void bl_coord_2d(FTYPE *X, FTYPE *V1, FTYPE *V2)
{
  FTYPE V[NDIM];

  bl_coord( X, V );
  *V1 = V[1];
  *V2 = V[2];

}

// Remaps the slopes (dq) such that one avoids using the boundary (ghost) zones 
//
// |=interface
// i=zone center of ith zone
//
// |              |       dq(i)        |
// |         pl(i)|pr(i)    i          |
// |              |pleft(i)   pright(i)|

void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE *p2interp_l, FTYPE *p2interp_r )
{
  extern int choose_limiter(int dir, int i, int j, int k, int pl);


  int ijk = idel * i + jdel * j + kdel * k;
  int ijkabs = startpos[dir] + ijk;
  int locallim;
  int pl,pliter;

  PLOOP(pliter,pl) {
    locallim = choose_limiter(dir, i,j,k,pl);
    if( (locallim<PARA) && (LIMADJUST==0) ) {
      //borrow the slope for reconstruction near the boundary from the neigbouring grid cell further from the boundary
      //lower boundary
      if( reconavoidbcdn[dir] == 1 && ijkabs == 0 ) {
        MACP0A1(dq,i,j,k,pl) = MACP0A1(dq,i + idel,j + jdel,k + kdel,pl);
      }
      //upper boundary
      if( reconavoidbcup[dir] == 1 && ijkabs == totalsize[dir] - 1 ) {
        MACP0A1(dq,i,j,k,pl) = MACP0A1(dq,i - idel,j - jdel,k - kdel,pl);
      }
    }
  }


  if(dir==2){
    if(ijkabs==totalsize[dir]-2){
      MACP0A1(dq,i,j,k,B2)=(0.0-MACP0A1(p2interp,i,j,k,B2))/(1.5);
      MACP0A1(dq,i,j,k,U2)=(0.0-MACP0A1(p2interp,i,j,k,U2))/(1.5);
    }
    else if(ijkabs==totalsize[dir]-1){
      MACP0A1(dq,i,j,k,B2)=(0.0-MACP0A1(p2interp,i,j,k,B2))/(0.5);
      MACP0A1(dq,i,j,k,U2)=(0.0-MACP0A1(p2interp,i,j,k,U2))/(0.5);
    }
  }

}

// Adjusts p_l & p_r at the boundary such that they do not depend on ghost zones;  
// NOTE: this fcn assumes remapdq() has been run before reconstruction to adjust the slopes (dq) or reconstruction
//
// |=interface
// i=zone center of ith zone
//
// |              |       dq(i)        |
// |         pl(i)|pr(i)    i          |
// |              |pleft(i)   pright(i)|
void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, 
                FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE *p2interp_l, FTYPE *p2interp_r )
{
  extern int choose_limiter(int dir, int i, int j, int k, int pl);

  int ijk = idel * i + jdel * j + kdel * k;
  int ijkabs = startpos[dir] + ijk;
  int locallim;
  int pl,pliter;

  PLOOP(pliter,pl) {
    locallim = choose_limiter(dir, i,j,k,pl);
    if( (locallim<PARA) && (LIMADJUST==0) ) {
      //borrow the slope for reconstruction near the boundary from the neigbouring grid cell further from the boundary
      //lower boundary
      if( reconavoidbcdn[dir] == 1 && ijkabs == 0 ) {
        p2interp_l[pl] = p2interp_r[pl];
      }
      //upper boundary
      if( reconavoidbcup[dir] == 1 && ijkabs == totalsize[dir] ) {
        p2interp_r[pl] = p2interp_l[pl];
      }
    }
  }
}


