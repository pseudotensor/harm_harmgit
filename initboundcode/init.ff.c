
/* 
 *
 NS disk or not and various NS field setups
 *
 *
 */

#include "decs.h"


#define KEPDISK 0
#define DONUTDISK 1
#define NODISK 2

#define DISKTYPE NODISK




#define SLOWFAC 1.0  /* reduce u_phi by this amount */

SFTYPE rhomaxold,umaxold,rhomax=0,umax=0,bsq_max=0,beta,rin;


// 0 = A_\phi in prime coordinate basis at CORN
// 1 = B^r in prime coordinate basis at FACE1
// 2 = B^\theta in prime coordinate basis at FACE2
// pother[]
FTYPE (*brface)[NSTORE2][NSTORE3];
FTYPE (*aphicorn)[NSTORE2][NSTORE3];
FTYPE Vpar;
FTYPE Omegastar;
FTYPE Bpole;

int DISKDIR;
// -1 : disk special wedge or no disk
// else dir as in global.h

#define SLOWFAC 1.0  /* reduce u_phi by this amount */

static FTYPE mm;

#define NTHETAMAX 10005

static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac,mydThetac;
static int NTHETA;

static FTYPE  INDEXN;  //SASMARKX:  changed to FTYPE from int, otherwise INDEXN = (1.0/4.0) would be truncated to zero


// 0 : split monopole
// 1 : monopole
// 2 : Wald
// 3 : BZ paraboloidal
// 4 : GRMHD nearly-paraboloidal
// 5 : BZ para but monopole near horizon
// 6 : GRMHD nearly-paraboloidal but monopole near horizon
// 7 : constant velocity Ramesh disk // can do cylindrical coordinatse also
#define PROBLEMTYPE 7

#define RAMESHTYPE 5  // 1 -- nu = 1, M = 0.25, s = 0.5 solution





int pre_init_specific_init(void)
{
  void get_nearlypara_data(void);
  void get_ramesh_data(void);

  // used in boundary conditions
  aphicorn = pother[APHICORN];
  brface = pother[B1FACE1];

  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  Omegastar=0;

  get_ramesh_data();

  Vpar=0.0; // parallel-to-field velocity at surface of NS in fraction of speed of light (orthonormal basis)

  Bpole = 0.0;

  return(0);
}


int init_conservatives(FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*Utemp)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  extern int pi2Uavg(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*Upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*Uavg)[NSTORE2][NSTORE3][NPR]);

  trifprintf("begin init_conservatives\n");
  pi2Uavg(p, Utemp, U);
  trifprintf("end init_conservatives\n");

  return(0);

}


int post_init_specific_init(void)
{
  FTYPE X[NDIM],V[NDIM];
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;
  //UTOPRIMVERSION = UTOPRIMCOMPARE;

  //  UTOPRIMVERSION =   UTOPRIM5D2;
  //UTOPRIMVERSION =   UTOPRIM5D1;
  //UTOPRIMVERSION =   UTOPRIM2D;
  //  UTOPRIMVERSION =   UTOPRIM1D;
  //UTOPRIMVERSION =   UTOPRIM1DOPT;
  //UTOPRIMVERSION =   UTOPRIM1DFINAL;
  //UTOPRIMVERSION =   UTOPRIM2DFINAL;

  // to ease into failures :) :(
  //dt=1E-5; // SUPERGODMARK

  coord(-N1BND, 0, 0, FACE1, X);
  bl_coord(X, V);
  if(V[1]<0.0){
    dualfprintf(fail_file,"Cannot have negative inner radius (r=%21.15g)\n",V[1]);
    myexit(235);
  }



  return(0);
}

int init_grid(void)
{

#if(MCOORD==SPCMINKMETRIC)


  if(0){
    R0 = -2.0;
    Rout = 1E3;
  }
  else if(1){
    R0 = -3.0;
    Rout = 1E4;
  }


  //  Rin=setRin(setihor());
  Rin = 1.0;

  hslope=1.0;

  // define coordinate type
  defcoord = RAMESHCOORDS_HALFDISK;

  /* output choices */
  tf = 10.0*Rout;

  DTd = 5;   /* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;   /* logfile frequency, in units of M */
  DTi = tf/100.0;   /* image file frequ., in units of M */
  DTdebug = DTd; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;   /* restart file period in steps */



#else
#error MCOORD != SPCMINKMETRIC
#endif

  return(0);
}

int init_global(void)
{

  if(1){// pulsar ffde
    GAMMAMAX=1000.0;
    GAMMAFAIL=100.0*GAMMAMAX;
    BSQORHOLIMIT=20.0;
    BSQOULIMIT=100.0;
    GAMMADAMP=5.0;

    rescaletype=1;
  }

  POSDEFMETRIC=0;


  SAFE=1.3;
  //  cour = 0.9;
  cour=0.8;
  //  cour = 0.5;

  ///////////////////////
  //
  // ENO-RELATED STUFF
  //
  ///////////////////////
  //  avgscheme=WENO5BND;
  avgscheme[1]=avgscheme[2]=avgscheme[3]=DONOR;
  
  do_transverse_flux_integration = 0;
  do_source_integration = 0;
  do_conserved_integration = 0;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;


#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  //lim = WENO5BND;
  //lim = PARA;
  TIMEORDER=4;
  //TIMEORDER=2;
  //  TIMEORDER=1;
  //lim[1]=lim[2]=lim[3] = DONOR;
  //lim[1]=lim[2]=lim[3] = WENO5BND;
  lim[1]=lim[2]=MC;
  //  lim[1]=lim[2]=MINM;
  //lim = MC;
  //  TIMEORDER=2;
  //  DOENOFLUX = ENOFINITEVOLUME;
  DOENOFLUX = NOENOFLUX;
  fluxmethod=HLLFLUX;
  //fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  //FLUXB = FLUXCTHLL;
  //  UTOPRIMVERSION=UTOPRIM5D1;

#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim[1]=lim[2]=lim[3] = MC;
  TIMEORDER=2;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif



  ranc(1,7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=NOCOOLING;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=NSSURFACE;

  DISKDIR=X2DN;
  //DISKDIR=X2UP;

  if(DISKDIR==X2DN){
    BCtype[X2UP]=POLARAXIS;
  }
  else if(DISKDIR==X2UP){
    BCtype[X2DN]=POLARAXIS;
  }
  else{
    dualfprintf(fail_file,"Please define DISKDIR properly you bitch\n");
    myexit(666);
  }

  BCtype[DISKDIR]=DISKSURFACE;

  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;




  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }





  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  //get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  //set_atmosphere(0,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  //if(pr[RHO]<pratm[RHO]){
  //  PLOOP(pliter,pl) pr[pl]=pratm[pl];
  //}
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}

int init_grid_post_set_grid(void)
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=sqrt(-1);
  Risco=sqrt(-1); //rmso_calc(PROGRADERISCO);

  return(0);

}



int init_primitives(FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{
  int whichvel, whichcoord;
  int initreturn;
  int i = 0, j = 0, k = 0, l;
  int pl,pliter;
  struct of_geom geom;
  FTYPE r,th,X[NDIM],V[NDIM];
  FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  int normalize_densities(FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
  int init_vpot(int l, int i, int j, int k, FTYPE *A);
  int normalize_field(FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
  int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *p);
  int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *pr);
  int init_vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pr)[NSTORE2][NSTORE3][NPR]);
  void set_analytical_face(void);


  ///////////////////////////////////
  //
  // Set some analytical solutions to be used by bounds()
  //
  ///////////////////////////////////
  // get analytical field for surface
  trifprintf("Set analytical face\n");
  set_analytical_face();


  ///////////////////////////////////
  //
  // Assign primitive variables
  //
  ///////////////////////////////////
  trifprintf("Assign primitives\n");

  // assume we start in bl coords and convert to KSprim
  FULLLOOP{
    initreturn=init_dsandvels(&whichvel, &whichcoord,i,j,k,MAC(p,i,j,k)); // request densities for all computational centers
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel,whichcoord,MAC(p,i,j,k), i,j,k) >= 1)
        FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
  }

  /////////////////////////////
  //
  // normalize density if wanted
  //
  /////////////////////////////// 
  // at this point densities are still standard, so just send "p"
  trifprintf("Normalize densities\n");
  normalize_densities(p);


  /////////////////////////////
  //
  // Define an atmosphere if wanted
  //
  /////////////////////////////// 

#if(DISKTYPE!=NODISK)

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // normalized atmosphere
  trifprintf("Add atmosphere\n");
  ZLOOP{
    initreturn=init_atmosphere(&whichvel, &whichcoord,i,j,k,MAC(p,i,j,k));
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel, whichcoord,MAC(p,i,j,k), i,j,k) >= 1)
        FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
  }
#endif
#endif



  // copy over initial solution as analytic solution
  // SET ANALYTIC SOLUTION FROM vector potential-based solution
  // NEEDED FOR BOUND in case uses panalytic
  COMPZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND){
    PLOOP(pliter,pl) MACP0A1(panalytic,i,j,k,pl)=MACP0A1(p,i,j,k,pl);

    // 0 out these things so dump files are readable by SM
    PLOOP(pliter,pl) MACP0A1(udump,i,j,k,pl)=0.0;
#if(CALCFARADAYANDCURRENTS)
    DLOOPA(pl) MACP0A1(jcon,i,j,k,pl)=0.0;
    for(pl=0;pl<NUMFARADAY;pl++) MACP0A1(fcon,i,j,k,pl)=0.0;
#endif
  }

  pdump=panalytic;
  // dump analytic solution
  if (dump(9995) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }
  pdump=p;




  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #1\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  trifprintf("Fixup#1\n");
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  trifprintf("Bound#1\n");
  MYFUN(bound_prim(STAGEM1,p),"init.c:init()", "bound_prim()", 1);

  trifprintf("pre_fixup#1\n");
  MYFUN(pre_fixup(STAGEM1,p),"init.c:init()", "pre_fixup()", 1);
#endif


  /////////////////////////////
  //
  // Initialize field from vector potential
  //
  /////////////////////////////// 
  A=emf; // dummy memory space not used till computation so safe.


  trifprintf("Initialize field from vector potential\n");
  COMPFULLLOOPP1{
    for(l=1;l<=3;l++) MACP1A0(A,l,i,j,k) = 0.;
  }

  COMPFULLLOOPP1{
    // GODMARK: Caution: Possible to use quantity off grid
    // (e.g. density) to define lower corner value of A, which then
    // defines B at center for lower cells.
    // Do have *grid* quantities for everywhre A is.
    for(l=1;l<=3;l++) init_vpot(l,i,j,k,&MACP1A0(A,l,i,j,k)); // request vector potential for all computational corners
  }
  trifprintf("Initialize field from vector potential assign\n");

  init_vpot2field(A,p);

  normalize_field(p);

  // copy over initial solution as analytic solution
  // SET ANALYTIC SOLUTION FROM vector potential-based solution
  COMPZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND) PLOOP(pliter,pl){
    MACP0A1(panalytic,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
  }


  pdump=panalytic;
  // dump analytic solution
  if (dump(9996) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }
  pdump=p;


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #2\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  trifprintf("Fixup#1\n");
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  trifprintf("Bound#1\n");
  MYFUN(bound_prim(STAGEM1,p),"init.c:init()", "bound_prim()", 1);

  trifprintf("pre_fixup#1\n");
  MYFUN(pre_fixup(STAGEM1,p),"init.c:init()", "pre_fixup()", 1);
#endif

  pdump=panalytic;
  // dump analytic solution
  if (dump(9997) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }
  pdump=p;


  // dump other
  if (dumpother(0) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }



  return(0);


}


// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int ii, int jj, int kk, FTYPE *pr)
{
  FTYPE X[NDIM], V[NDIM], r, th;
  struct of_geom geom;

  FTYPE Fcov[NDIM][NDIM];
  FTYPE Mcon[NDIM][NDIM];
  FTYPE Mcov[NDIM][NDIM];
  FTYPE etacov[NDIM],etacon[NDIM];
  FTYPE Ecov[NDIM],Econ[NDIM],Bcov[NDIM],Bcon[NDIM];
  FTYPE  alpha;

  void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM]);
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);

  struct of_state q;
  FTYPE faradaytest[NDIM][NDIM];


  if(EOMTYPE!=EOMFFDE){
    dualfprintf(fail_file,"Are you sure?\n");
    myexit(1);
  }


  coord(ii, jj, kk, CENT, X);
  bl_coord(X, V);
  r = V[1];
  th = V[2];
  get_geometry(ii, jj, kk, CENT, &geom); 

#define B0 1.0

  pr[RHO]=pr[UU]=0.0;
  pr[U1]=pr[U2]=pr[U3]=0.0;
  pr[B2]=pr[B3]=0;

  pr[B1]=0.0; // not used   
  pr[B2]=0.0;
  pr[B3]=0.0;

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;

  return(0); // shouldn't reach here
}








// assumes we are fed the true densities
int normalize_densities(FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;

#if(DISKTYPE!=NODISK && 0)
  rhomax=0;
  umax=0;
  ZLOOP{
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];

    if (MACP0A1(p,i,j,k,RHO) > rhomax)   rhomax = MACP0A1(p,i,j,k,RHO);
    if (MACP0A1(p,i,j,k,UU) > umax && r > rin)    umax = MACP0A1(p,i,j,k,UU);
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %21.15g umax: %21.15g\n", rhomax, umax);


  ZLOOP{
    MACP0A1(p,i,j,k,RHO) /= rhomax;
    MACP0A1(p,i,j,k,UU) /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

#endif

  return(0);
}














#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define DIPOLEFIELD 3 // i.e. no disk field with dipole field (like pulsar model)
#define DISKDIPOLE 4


#define FIELDTYPE DIPOLEFIELD


// assumes normal field in pr
int init_vpot(int l, int i, int j, int k, FTYPE *A)
{
  FTYPE X[NDIM],V[NDIM];
  struct of_geom geom;
  FTYPE dxdxp[NDIM][NDIM];
  void set_ramesh_solution( int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE *A);
  int pp;




  if(l==3){// A_\phi

    pp=CORN3;
    get_geometry(i, j, k, pp, &geom);
    coord(i, j, k, pp, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    //set vector potential and save the initial (analytic value) of omegaf
    set_ramesh_solution( l, &geom, V, p, A ); 

    // assume we were setting A_\phi and not A_{x^3}, so need to correct
    // correct for 2\pi factor
    *A *=dxdxp[PH][PH];

  }


  

  return(0);

}

// Vector potential as function of real coordinates
// assumes normal field in pr
int init_vpot_V(int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE *A)
{
  int get_analytic_vpot( int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *A );

  *A=0;

  if(l==3){// A_\phi
    get_analytic_vpot( l, geom, V, prim, A );
  }

  return(0);

}


void set_analytical_face(void)
{
  FTYPE X[NDIM],V[NDIM];
  struct of_geom geom;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE dxdxph[NDIM][NDIM];
  FTYPE dxdxpl[NDIM][NDIM];
  FTYPE Xh[NDIM],Vh[NDIM];
  struct of_geom geomh,geoml;
  FTYPE Xl[NDIM],Vl[NDIM];
  int init_vpot_V(int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE *A);
  FTYPE Ah,Al;
  FTYPE Ahh,Ahl,Alh,All;
  FTYPE dxdxphh[NDIM][NDIM];
  FTYPE dxdxpll[NDIM][NDIM];
  FTYPE dxdxphl[NDIM][NDIM];
  FTYPE dxdxplh[NDIM][NDIM];
  struct of_geom geomhh,geomll,geomhl,geomlh;
  FTYPE Xhh[NDIM], Xll[NDIM], Xhl[NDIM], Xlh[NDIM];
  FTYPE Vhh[NDIM], Vll[NDIM], Vhl[NDIM], Vlh[NDIM];
  int i,j,k;
  int jj;
  int l;
  FTYPE prtemp[NPR];
  int init_nodisk(int pp, int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr);
  int whichcoord,whichvel;
  int get_analytic_omegaf( struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *omegaf );
  FTYPE omegaf;



  // get analytical A_\phi at corner
  FULLLOOP{

    get_geometry(i, j, k, CORN3, &geom);
    coord(i, j, k, CORN3, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    l=3;
    init_vpot_V(l,&geom,V,p,&GLOBALMACP1A0(pother,APHICORN,i,j,k));
    
    // assume we were setting A_\phi and not A_{x^3}, so need to correct
    // correct for 2\pi factor
    GLOBALMACP1A0(pother,APHICORN,i,j,k) *=dxdxp[PH][PH];

  }

  // get analytical A_\phi at CENT
  FULLLOOP{

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    l=3;
    init_vpot_V(l,&geom,V,p,&GLOBALMACP1A0(pother,APHICENT,i,j,k));
    
    // assume we were setting A_\phi and not A_{x^3}, so need to correct
    // correct for 2\pi factor
    GLOBALMACP1A0(pother,APHICENT,i,j,k) *=dxdxp[PH][PH];

  }

  // get analytical B^r in prime coordinate basis at FACE1
  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);

    DLOOPA(jj) Xh[jj]=Xl[jj]=X[jj];

    Xh[2]+=1E-5;
    Xl[2]-=1E-5;

    bl_coord(Xh, Vh);
    bl_coord(Xl, Vl);

    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);


    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah); // geom not used
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al); // geom not used
    Al *=dxdxpl[PH][PH];

    // use analytical definition of B^{x^1}
    GLOBALMACP1A0(pother,B1FACE1,i,j,k)=(Ah-Al)/(geom.g*(Xh[2]-Xl[2]));
    
  }

  // get analytical B^r in prime coordinate basis at CENT
  FULLLOOP{

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);

    DLOOPA(jj) Xh[jj]=Xl[jj]=X[jj];

    Xh[2]+=1E-5;
    Xl[2]-=1E-5;

    bl_coord(Xh, Vh);
    bl_coord(Xl, Vl);


    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);

    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah);
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al);
    Al *=dxdxpl[PH][PH];

    // use analytical definition of B^{x^1}
    GLOBALMACP1A0(pother,B1CENT,i,j,k)=(Ah-Al)/(geom.g*(Xh[2]-Xl[2]));


  }

  // get ZEUS-like B^r in prime coordinate basis at FACE1
  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);

    get_geometry(i, j, k, CORN3, &geoml);
    coord(i, j, k, CORN3, Xl);
    bl_coord(Xl, Vl);

    get_geometry(i, jp1, k, CORN3, &geomh);
    coord(i, jp1, k, CORN3, Xh);
    bl_coord(Xh, Vh);

    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);

    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah);
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al);
    Al *=dxdxpl[PH][PH];

    // use ZEUS definition of B^{x^1}
    GLOBALMACP1A0(pother,B1ZEUS,i,j,k)=(Ah-Al)/(geom.g*(Xh[2]-Xl[2]));


  }

  // get FLUXCT-like B^r in prime coordinate basis at FACE1
  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);

    get_geometry(im1, j, k, FACE2, &geomll);
    coord(im1, j, k, FACE2, Xll);
    bl_coord(Xll, Vll);

    get_geometry(i, jp1, k, FACE2, &geomhh);
    coord(i, jp1, k, FACE2, Xhh);
    bl_coord(Xhh, Vhh);

    get_geometry(im1, jp1, k, FACE2, &geomlh);
    coord(im1, jp1, k, FACE2, Xlh);
    bl_coord(Xlh, Vlh);

    get_geometry(i, j, k, FACE2, &geomhl);
    coord(i, j, k, FACE2, Xhl);
    bl_coord(Xhl, Vhl);


    dxdxprim(Xhh, Vhh, dxdxphh);
    dxdxprim(Xll, Vll, dxdxpll);
    dxdxprim(Xhl, Vhl, dxdxphl);
    dxdxprim(Xlh, Vlh, dxdxplh);

    l=3;
    init_vpot_V(l,&geom,Vhh,p,&Ahh);
    Ahh *=dxdxphh[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vll,p,&All);
    All *=dxdxpll[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vhl,p,&Ahl);
    Ahl *=dxdxphl[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vlh,p,&Alh);
    Alh *=dxdxplh[PH][PH];

    GLOBALMACP1A0(pother,B1FLUXCTFACE1,i,j,k)=((Ahh+Alh)*0.5 - (Ahl+All)*0.5)/(geom.g*dx[2]);


  }

  
  // get analytical B^\theta in prime coordinate basis at FACE2
  FULLLOOP{

    get_geometry(i, j, k, FACE2, &geom);
    coord(i, j, k, FACE2, X);
    bl_coord(X, V);

    DLOOPA(jj) Xh[jj]=Xl[jj]=X[jj];

    Xh[1]+=1E-5;
    Xl[1]-=1E-5;

    bl_coord(Xh, Vh);
    bl_coord(Xl, Vl);

    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);

    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah);
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al);
    Al *=dxdxpl[PH][PH];

    // use analytical definition of B^{x^2}
    if(geom.g!=0.0){
      GLOBALMACP1A0(pother,B2FACE2,i,j,k)=-(Ah-Al)/(geom.g*(Xh[1]-Xl[1]));
    }
    else{
      GLOBALMACP1A0(pother,B2FACE2,i,j,k)=0.0;
    }

  }

  // get analytical B^\theta in prime coordinate basis at CENT
  FULLLOOP{

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);

    DLOOPA(jj) Xh[jj]=Xl[jj]=X[jj];

    Xh[1]+=1E-5;
    Xl[1]-=1E-5;

    bl_coord(Xh, Vh);
    bl_coord(Xl, Vl);

    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);

    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah);
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al);
    Al *=dxdxpl[PH][PH];

    // use analytical definition of B^{x^2}
    GLOBALMACP1A0(pother,B2CENT,i,j,k)=-(Ah-Al)/(geom.g*(Xh[1]-Xl[1]));

  }

  // get ZEUS-like B^\theta in prime coordinate basis at FACE2
  FULLLOOP{

    get_geometry(i, j, k, FACE2, &geom);
    coord(i, j, k, FACE2, X);
    bl_coord(X, V);

    get_geometry(i, j, k, CORN3, &geoml);
    coord(i, j, k, CORN3, Xl);
    bl_coord(Xl, Vl);

    get_geometry(ip1, j, k, CORN3, &geomh);
    coord(ip1, j, k, CORN3, Xh);
    bl_coord(Xh, Vh);

    dxdxprim(Xh, Vh, dxdxph);
    dxdxprim(Xl, Vl, dxdxpl);

    l=3;
    init_vpot_V(l,&geom,Vh,p,&Ah);
    Ah *=dxdxph[PH][PH];

    l=3;
    init_vpot_V(l,&geom,Vl,p,&Al);
    Al *=dxdxpl[PH][PH];

    // use ZEUS-like definition of B^{x^2}
    if(geom.g!=0.0){
      GLOBALMACP1A0(pother,B2ZEUS,i,j,k)=-(Ah-Al)/(geom.g*(Xh[1]-Xl[1]));
    }
    else{
      GLOBALMACP1A0(pother,B2ZEUS,i,j,k)=0.0;
    }

  }


  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    // in coordinate basis
    GLOBALMACP1A0(pother,OMEGAFFACE1,i,j,k)=mm*pow(B0*pow(Rin,nu),-1.0/nu)/dxdxp[3][3];

  }



  FULLLOOP{

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    //get field angular velocity
    get_analytic_omegaf( &geom, V, p, &omegaf);

    // in coordinate basis
    GLOBALMACP1A0(pother,OMEGAFCENT,i,j,k)=omegaf / dxdxp[3][3];

  }

  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    prtemp[RHO] = 0; //SASMARK
    //init_nodisk(FACE1, &whichvel, &whichcoord, i, j, k, prtemp);
    //    set_density_floors(&geom,prtemp,prtemp);
    GLOBALMACP1A0(pother,RHOFACE1,i,j,k)=prtemp[RHO];

  }

  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    // assume Vpar in orthonormal basis
    if(GLOBALMACP1A0(pother,B1FACE1,i,j,k)>0){
      GLOBALMACP1A0(pother,VPARFACE1,i,j,k)=Vpar;
    }
    else{
      GLOBALMACP1A0(pother,VPARFACE1,i,j,k)=-Vpar;
    }

  }

  FULLLOOP{

    get_geometry(i, j, k, CENT, &geom);
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    // assume Vpar in orthonormal basis
    GLOBALMACP1A0(pother,VPARCENT,i,j,k)=Vpar;


  }

  FULLLOOP{

    get_geometry(i, j, k, FACE1, &geom);
    coord(i, j, k, FACE1, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    // not used
    GLOBALMACP1A0(pother,UUFACE1,i,j,k)=0.0;

  }
  



  //disk SASMARK
  FULLLOOP{ 
    get_geometry(i, j, k, FACE2, &geom);
    coord(i, j, k, FACE2, X);
    bl_coord(X, V);
    dxdxprim(X, V, dxdxp);

    // in coordinate basis 
    //this expression is true in midplane (used only there), since there V[1] = R
    GLOBALMACP1A0(pother,OMEGAFFACE2,i,j,k)= (mm/V[1]) / dxdxp[3][3];
    GLOBALMACP1A0(pother,B2FACE2SIMPLE,i,j,k) = -(B0*nu*pow(V[1],nu-2)) * (geom.gcov[GIND(2,2)]);

    prtemp[RHO] = 0; //SASMARK
    //init_nodisk(FACE2, &whichvel, &whichcoord, i, j, k, prtemp);
    //    set_density_floors(&geom,prtemp,prtemp);
    GLOBALMACP1A0(pother,RHOFACE2,i,j,k)=prtemp[RHO];

    // assume Vpar in orthonormal basis
    GLOBALMACP1A0(pother,VPARFACE2,i,j,k)=Vpar;

    // not used
    GLOBALMACP1A0(pother,UUFACE2,i,j,k)=0.0;

  }



}



int init_vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  extern int vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
  int i,j,k;

  if(FLUXB==FLUXCTHLL){
    // then assume analytically setting B^i
    FULLLOOP{
      MACP0A1(pr,i,j,k,B1)=GLOBALMACP1A0(pother,B1CENT,i,j,k);
      MACP0A1(pr,i,j,k,B2)=GLOBALMACP1A0(pother,B2CENT,i,j,k);
      MACP0A1(pr,i,j,k,B3)=0.0;
    }
  }

  return(vpot2field(A,pr));

}



// assumes normal field definition
int normalize_field(FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;

#if(0) // assume field already setup to be normalized relative to NS field

  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, k, CENT, &geom);    

    if(FIELDTYPE==VERTFIELD){
      coord(i, j, k, CENT, X);
      bl_coord(X, V);
      r=V[1];
      th=V[2];
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
        if (bsq_calc(MAC(p,i,j,k), &geom, &bsq_ij) >= 1)
          FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
 
        if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
      }
    }
    else{
      if (bsq_calc(MAC(p,i,j,k), &geom, &bsq_ij) >= 1)
        FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  trifprintf("initial beta: %21.15g (should be %21.15g)\n", beta_act,beta);
  norm = sqrt(beta_act / beta);
  
  bsq_max = 0.;
  ZLOOP {
    MACP0A1(p,i,j,k,B1) *= norm;
    MACP0A1(p,i,j,k,B2) *= norm;
    MACP0A1(p,i,j,k,B3) *= norm;

    get_geometry(i, j, k, CENT, &geom);
    if (bsq_calc(MAC(p,i,j,k), &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
    if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    
  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);
  trifprintf("final beta: %21.15g (should be %21.15g)\n", beta_act,beta);
#endif

  return(0);
}



#undef SLOWFAC

SFTYPE lfish_calc(SFTYPE r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
           ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
            sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
            ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) +
                                          pow(a,
                                              2) * (2. + r))) / sqrt(1 +
                                                                     (2.
                                                                      *
                                                                      a)
                                                                     /
                                                                     pow
                                                                     (r,
                                                                      1.5)
                                                                     -
                                                                     3.
                                                                     /
                                                                     r)))
          / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
             (pow(a, 2) + (-2. + r) * r))
          );
}


FTYPE nz_func(FTYPE R)
{
  return(
         sqrt(
              (3.*a*a - 4.*a*sqrt(R) + R*R)/
              pow(R*(a + pow(R,1.5)),2)
              )
         ) ;


}



// for axisymmetric case
// converts metric' to metric
void metptomet_axisym(FTYPE *uconmetp, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmet)
{
  int j,k;
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM];
  /* transform ucon */
  // this is u^j = T^j_k u^k, as in above

  DLOOPA(j) tmp[j] = 0.;
  DLOOP(j,k) tmp[j] += dxdxp[j][k] * uconmetp[k];
  DLOOPA(j) uconmet[j] = tmp[j];

}
// converts metric to metric'
void mettometp_axisym(FTYPE *uconmet, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp)
{
  int j,k;
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM];
  FTYPE bottom;

  /* transform ucon */
  // this is u^j = (iT)^j_k u^k, as in mettobl() above

  idxdxp[0][0]=1.0/dxdxp[1][1];
  idxdxp[0][1]=0.0;
  idxdxp[0][2]=0.0;
  idxdxp[0][3]=0.0;

  bottom=(dxdxp[2][2]*dxdxp[1][1]-dxdxp[2][1]*dxdxp[1][2]);

  idxdxp[1][0]=0.0;
  idxdxp[1][1]=dxdxp[2][2]/bottom;
  idxdxp[1][2]=dxdxp[1][2]/(-bottom);
  idxdxp[1][3]=0.0;

  idxdxp[2][0]=0.0;
  idxdxp[2][1]=dxdxp[2][1]/(-bottom);
  idxdxp[2][2]=dxdxp[1][1]/bottom;
  idxdxp[2][3]=0.0;

  idxdxp[3][0]=0.0;
  idxdxp[3][1]=0.0;
  idxdxp[3][2]=0.0;
  idxdxp[3][3]=1.0/dxdxp[3][3];


  DLOOPA(j) tmp[j] = 0.;
  DLOOP(j,k) tmp[j] += idxdxp[j][k] * uconmet[k];
  DLOOPA(j) uconmetp[j] = tmp[j];


}


void getconsts(FTYPE *uconmetp, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconconsts)
{
  int j,k;
  FTYPE uconortho[NDIM],uconmet[NDIM];
  void metptomet_axisym(FTYPE *uconmetp, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmet);
  FTYPE sinh1,cosh1;
  FTYPE sinh2,cosh2;
  FTYPE sinh3,cosh3;

  // simplest case
  //  sinh1=sinh2=sin(V[2]);
  //cosh1=cosh2=cos(V[2]);
  sinh1=sinh2=1.0;
  cosh1=cosh2=1.0;
  sinh3=1.0;
  cosh3=1.0;

  

  // get 4-vector in normal metric coordinates
  metptomet_axisym(uconmetp, ptrgeom, dxdxp, uconmet);

  // assume spherical polar and get orthonormal versions

  //  uconortho[0]=uconmet[0]*sqrt(fabs(ptrgeom->gcov[GIND(0,0)]));
  uconortho[1]=uconmet[1];
  uconortho[2]=uconmet[2]*V[1];
  uconortho[3]=uconmet[3]*V[1]*sinh3;

  // should be ok for radial interpolation

  uconconsts[1]=uconortho[1]*V[1]*V[1]*V[1]/(SMALL+cosh1); // may be problem at equator
  uconconsts[2]=uconortho[2]*V[1]*V[1]*V[1]/(SMALL+sinh2); // may be problem at pole
  uconconsts[3]=uconortho[3]*V[1]; // assume B^{\hat{\phi}}\propto 1/r

  // ignore time and phi components since not interpolating them
}

void undoconsts(FTYPE *uconconst, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp)
{
  int j,k;
  FTYPE uconortho[NDIM],uconmet[NDIM];
  void mettometp_axisym(FTYPE *uconmet, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp);
  FTYPE sinh1,cosh1;
  FTYPE sinh2,cosh2;
  FTYPE sinh3,cosh3;

  // simplest case
  //  sinh1=sinh2=sin(V[2]);
  //cosh1=cosh2=cos(V[2]);
  sinh1=sinh2=1.0;
  cosh1=cosh2=1.0;
  sinh3=1.0;
  cosh3=1.0;


  // now assign new function given constant interpolation from reference value
  uconortho[1]=uconconst[1]/(V[1]*V[1]*V[1])*cosh1; // might cause problems if not exactly symmetric solution
  uconortho[2]=uconconst[2]/(V[1]*V[1]*V[1])*sinh2;
  uconortho[3]=uconconst[3]/V[1];

  //  uconmet[0]=uconortho[0]/sqrt(fabs(ptrgeom->gcov[GIND(0,0)]));
  uconmet[1]=uconortho[1];
  uconmet[2]=uconortho[2]/V[1];
  uconmet[3]=uconortho[3]/(V[1]*(SMALL+sinh3)); // might be problem at poles

  mettometp_axisym(uconmet, ptrgeom, dxdxp,uconmetp);
  

}


void getnewucon(FTYPE *uconmetpin, FTYPE *rV, struct of_geom *rptrgeom, FTYPE (*rdxdxp)[NDIM], FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetpout)
{
  void getconsts(FTYPE *uconmetp, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconconst);
  void undoconsts(FTYPE *uconconst, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp);
  FTYPE uconconst[NDIM];

  // get constants at reference location
  getconsts(uconmetpin, rV, rptrgeom, rdxdxp, uconconst);

  undoconsts(uconconst, V, ptrgeom, dxdxp, uconmetpout);

}


void get_ramesh_data(void)
{
  int i;
  FILE * inTheta;
  FTYPE ftemp,ftemp1,ftemp2;


  ////////////////////////
  // remove empty line at top of file and convert D to E
  // and set NTHETA
  //
  // CHOOSE which Ramesh problem
  //
  ///////////////////////////////

  if(myid==0){
#if(RAMESHTYPE==0)
    // for nu.75_m.5.txt
    Ttpow = 0; // not known right now
    NTHETA=1314;
    if( (inTheta=fopen("nu.75_m.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu.75_m.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==1)
    // for nu1.0_m.25.txt
    // T(\theta)\propto \theta^{Ttpow} -> \theta_j \propto r^{-1/(1+Ttpow/nu)}
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("nu1.0_m.25.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu1.0_m.25.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==125)
    // or nu1.0_m.25_hres.txt
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("nu1.0_m.25_hres.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu1.0_m.25_hres.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==2)
    // for nu.75_m.0005.txt
    Ttpow = 0.0; // not known right now
    NTHETA=1206;
    if( (inTheta=fopen("nu.75_m.0005.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open nu.75_m.0005.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==3)
    // for ramesh_nu1_m.25_s1.txt
    Ttpow = 1.0;
    NTHETA=2030;
    if( (inTheta=fopen("ramesh_nu1_m.25_s1.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1_m.25_s1.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==4)
    // for ramesh_nu.75_m.1_s.5.txt
    Ttpow = 3.0;
    NTHETA=1505;
    if( (inTheta=fopen("ramesh_nu.75_m.1_s.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu.75_m.1_s.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==5)
    // for ramesh_nu.75_m.1_s.28146.txt
    Ttpow = 1.361;
    NTHETA=2032;
    if( (inTheta=fopen("ramesh_nu.75_m.1_s.28146.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu.75_m.1_s.28146.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==6)
    // for ramesh_nu1.25_m.4_soM2.5.txt
    Ttpow = 0.0; // special solution
    NTHETA=1502;
    if( (inTheta=fopen("ramesh_nu1.25_m.4_soM2.5.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1.25_m.4_soM2.5.txt\n");
      myexit(100);
    }
#elif(RAMESHTYPE==7)
    // for ramesh_nu1.25_m.4_soM1.6.txt
    Ttpow = 2.0-nu; // for \nu>1
    NTHETA=1865;
    if( (inTheta=fopen("ramesh_nu1.25_m.4_soM1.6.txt","rt"))==NULL){
      dualfprintf(fail_file,"Cannot open ramesh_nu1.25_m.4_soM1.6.txt\n");
      myexit(100);
    }
#endif 

    if(NTHETA>=NTHETAMAX){
      dualfprintf(fail_file,"NTHETA=%d and NTHETAMAX=%d\n",NTHETA,NTHETAMAX);
      myexit(7);
    }


    ///////////////////////////////////////////
    //
    // READ IN DATA FILE
    //
    ////////////////////////////////////////////

    // skip first blank line
    //while(fgetc(inTheta)!='\n'); // skip first line, a comment

    // read in data file

    // read in first row (parameters)
    fscanf(inTheta,"%lf %lf %lf %lf",&nu,&mm,&ss,&ucrit);
    //dualfprintf(fail_file,"got0: %g %g %g %g\n",nu,mm,ss,ucrit);

    // NTHETA does not include first line



    for(i=0;i<NTHETA;i++){
      // while(!feof(inTheta)){
      fscanf(inTheta,"%lf %lf %lf",&ftemp,&ftemp1,&ftemp2);
      // Z/R itself
      theta[i]=ftemp; // ftemp=u= Z/R=tan(\theta) //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
      Thetavstheta[i]=ftemp1;
      dThetadthetavstheta[i]=ftemp2;
#if(0)
      // theta
      theta[NTHETA-1-i]=atan(1.0/ftemp); // true theta since file has ftemp=u= Z/R=tan(\theta)
      Thetavstheta[NTHETA-1-i]=ftemp1;
      dThetadthetavstheta[NTHETA-1-i]=ftemp2;
#endif
      //   dualfprintf(fail_file,"got1: %d %g %g %g\n",i,theta[NTHETA-1-i],Thetavstheta[NTHETA-1-i],dThetadthetavstheta[NTHETA-1-i]);
      //   i++;
    }
    // NTHETA=i;

    // \theta_j \propto r^{-jetalpha}
    jetalpha = 1.0/(1.0+(Ttpow/nu));

    trifprintf("Ramesh parameters: %g %g %g %g %g %g\n",nu,mm,ss,ucrit,Ttpow,jetalpha);


  }





  // send data to all CPUs
#if(USEMPI)
  MPI_Bcast(&nu,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&mm,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&ss,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&ucrit,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&Ttpow,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&jetalpha,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);


  MPI_Bcast(&NTHETA,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&theta,NTHETA,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&Thetavstheta,NTHETA,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dThetadthetavstheta,NTHETA,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
#endif


}



////////////////////////////////////////////////////
//
// Ramesh constant velocity disk
//
///////////////////////////////////////////////////
void set_ramesh_solution(int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  FTYPE Xc[NDIM],Vc[NDIM],rc,thc;
  int ii, jj, kk;
  struct of_geom geomc;
  FTYPE th2;
  FILE * inTheta;
  //  static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac;
  int dumi;
  int i,k;
  static int firsttime=1;
  FTYPE rprime,zprime;
  FTYPE B0para;
  FTYPE rtrans;
  //  int NTHETA;
  FTYPE Ac;
  FTYPE theu;
  int get_analytic_vpot( int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *A );
  int get_analytic_omegaf( struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *omegaf );
  void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer);
  FTYPE R_smooth;



  ///////////////////////////////////////////
  //
  // Already have data, now interpolate data to HARM grid
  //
  ////////////////////////////////////////////


  r = V[1];
  th = V[2];

  ii = geom->i;
  jj = geom->j;
  kk = geom->k;

  coord(ii, jj, kk, CENT, Xc);
  bl_coord(Xc, Vc);
  rc = Vc[1];
  thc = Vc[2];
  get_geometry(ii, jj, kk, CENT, &geomc);
  
  ////////
  //
  //done with setting geometries and r, th
  //
  ////////

  get_analytic_vpot( l, geom, V, prim, A );

  // for testing purposes (dump999x), place omegafanalytic into internal energy.  This is set back to 0 before evolution.
  if( ii < N1 + N1BND && jj < N2 + N2BND && kk < N3 + N3BND ) {  //SASMARKx: to avoid getting out of bounds
    //    MACP0A1(prim,ii,jj,kk,RHO)=*A;
    MACP0A1(prim,ii,jj,kk,B1)=0.0;
    MACP0A1(prim,ii,jj,kk,B2)=0.0;
    MACP0A1(prim,ii,jj,kk,B3)=0.0;
  } 

  ////////////////
  //
  // need omegaf analytic! (only used at center)
  //
  ////////////////

  //  get_analytic_omegaf( &geomc, Vc, prim, &MACP0A1(omegafanalytic,ii,jj,kk,0) );
  //MACP0A1(omegafanalytic,ii,jj,kk,0) /= dxdxp[3][3];

  // for testing purposes (dump999x), place omegafanalytic into internal energy.  This is set back to 0 before evolution.
  //  MACP0A1(prim,ii,jj,kk,UU)=MACP0A1(omegafanalytic,ii,jj,kk,0);

}

int get_analytic_omegaf( struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *omegaf )
{
  int get_analytic_vpot( int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *A );

  FTYPE r, th;
  FTYPE A;
  FTYPE Xc,Vc;
  FTYPE rc;
  
  int ii, jj, kk, pp;
  int phi_component = 3;



  r = V[1];
  th = V[2];

  ii = geom->i;
  jj = geom->j;
  kk = geom->k;
  pp = geom->p;

  get_analytic_vpot( phi_component, geom, V, prim, &A );

  //  if( r <= Rin * (1.0 + NUMEPSILON) ) {
  //  if( r < 0.5(Rin+rc) ) {
  //  if( r <= Rin || (startpos[1]+ii==0 && geom->p==FACE1) ) {
  if(is_inside_surface(X1DN,ii,jj,kk,pp)){
    // fixed \omega inside star (surface)
    *omegaf = mm*pow(B0*pow(Rin,nu),-1.0/nu); // equatorial value for all of star
  }
  else{
    *omegaf = mm*pow(A/Thetavstheta[0],-1.0/nu);
  }

  //  if(th<0.0 || th>M_PI ){ // assuming consistently using boundary condition where u^\phi is antisymmetric
  //    *omegaf *= -1.0;
  //  }

  return( 0 );
}




// Returns Ramesh's analytic expression for vector potential A_\phi
int get_analytic_vpot( int l, struct of_geom *geom, FTYPE *V, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *A )
{
  void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer);
  
  FTYPE r, th;
  FTYPE rprime, zprime;
  FTYPE theu;
  int i;
  FTYPE myTheta;
  FTYPE th2;
  FTYPE signt;
  FTYPE A3sign;



  if( l != 3 ) {
    *A = 0;
    return( 0 );
  }

  r = V[1];
  th = V[2];


  if(th<0.0){
    // function is even
    th2=-th;
  }
  else if(th<=M_PI*0.5){
    // theta is itself
    th2=th;
  }
  else if(th<=M_PI){
    //    th2=M_PI_2l; //M_PI-th;  //outflow aphi along r = const
    th2=M_PI-th;  // JONMARK
  }
  else{
    // function is even
    th2=-(M_PI-th);
  }



  // R
  rprime=fabs(r)*fabs(sin(th2));
  // z
  zprime=fabs(r)*fabs(cos(th2));

  theu = fabs(cos(th2)/sin(th2));

  // \theta is not uniform, so loop to find first theta
  for(i=0;i<NTHETA;i++){
    if(theta[i]>theu) break;   //SASMARK:  actually, theta[i] is not \theta but u = Z/R = tan(\theta)
  }

  if(i==NTHETA) i=NTHETA-1;
  
  interpfun(i, theu, theta, Thetavstheta, &myTheta);


  if(DISKDIR==-1){
  }
  else if(DISKDIR==X2DN){
    if(th<M_PI*0.5){
      myTheta = 2.0 - myTheta;
    }
  }
  else if(DISKDIR==X2UP){
    if(th>M_PI*0.5){
      myTheta = 2.0 - myTheta;
    }
  }


  // A_\phi = P(R,Z) = R^\nu T(u)
  *A = B0*pow(rprime,nu)*myTheta;

#if(FLIPGDETAXIS==0)
  // causes kink at pole in B^r and B^\theta with flip of gdet sign
  if(th<0.0) A3sign=-1.0;
  else if(th>M_PI) A3sign=-1.0;
  else A3sign=1.0;
#else
  A3sign=1.0; // to be used with gdet flip of sign
#endif

  *A *=A3sign;   //SASMARK:  B^r = dA_\phi/d\theta has to be symmetric, so flipping the sign is correct here

  
  return(0);
}


#define LINEARTYPE 0
#define LOGTYPE 1
#define QUADRATICTYPE 2

#define INTERPTYPE LOGTYPE
//#define INTERPTYPE LINEARTYPE
//#define INTERPTYPE QUADRATICTYPE

void interpfun(int i, FTYPE th2, FTYPE *theta, FTYPE *fun, FTYPE *answer)
{
  FTYPE slope,intercept;
  FTYPE slope1,slope2,xminusx0;
  FTYPE f0,f1,f2,x0,x1,x2;



  if(INTERPTYPE==LINEARTYPE){
    // linearly interpolate \Theta using \theta
    //      *answer = fun[i-1] + (fun[i]-fun[i-1])/(theta[i]-theta[i-1])*(th2-theta[i-1]);    
    slope = (fun[i]-fun[i-1])/(theta[i]-theta[i-1]);
    intercept = fun[i-1];
    *answer = slope*(th2-theta[i-1]) + intercept;

  }
  else if(INTERPTYPE==QUADRATICTYPE){
    // quadratically interpolate \Theta using \theta
    if(i-1<0){
      f0=fun[i];
      f1=fun[i+1];
      f2=fun[i+2];
      x0=theta[i];
      x1=theta[i+1];
      x2=theta[i+2]; 
    }
    else if(i+1>=NTHETA){
      f0=fun[i-2];
      f1=fun[i-1];
      f2=fun[i];
      x0=theta[i-2];
      x1=theta[i-1];
      x2=theta[i]; 
    }
    else{
      f0=fun[i-1];
      f1=fun[i];
      f2=fun[i+1];
      x0=theta[i-1];
      x1=theta[i];
      x2=theta[i+1];
    }

    slope2 = ((f0-f2)/(x0-x2) - (f2-f1)/(x2-x1))/(x0-x1);
    slope1 = (f0-f1)/(x0-x1) + (f0-f2)/(x0-x2) - (f2-f1)/(x2-x1);
    xminusx0 = (th2-x0);

    *answer = slope2*pow(xminusx0,2.0) + slope1*xminusx0 + f0;
  }
  else if(INTERPTYPE==LOGTYPE){
    // log interpolate \Theta using \theta
    slope = log(fun[i]/fun[i-1])/log(theta[i]/theta[i-1]);
    if(fabs(slope)<1E-10 || !isfinite(slope)) *answer=fun[0];  //SASMARKx
    else if(fun[i-1]<0.0){
      // assume bi-log
      *answer=-exp( slope*log(th2/theta[i-1])+log(-fun[i-1]) );
    }
    else *answer=exp( slope*log(th2/theta[i-1])+log(fun[i-1]) );

    //dualfprintf(fail_file,"ii=%d jj=%d slope=%g myTheta=%g\n",ii,jj,slope,myTheta);
  }

}

