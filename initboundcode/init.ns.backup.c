
/* 
 *
 NS disk or not and various NS field setups
 *
 *
 */

#include "decs.h"





#define SLOWFAC 1.0  /* reduce u_phi by this amount */

SFTYPE rhomaxold,umaxold,rhomax=0,umax=0,bsq_max=0,beta,rin;

int pre_init_specific_init(void)
{

  if(0){ // pulsar ffde case

    // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
    h_over_r=1.0;
    // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
    h_over_r_jet=h_over_r;


    // see avery_NS_GRB.nb
    // polar field for NS in units of $\sqrt{\rho_{disk} c^2}$
    // weak field
    //Bpole=2.1E-6;
    // strong field
    Bpole=1.0; 
    
    Omegastar=0.0216; // in units of per GM/c^3

  }
  else if(1){
    h_over_r=0.3;
    h_over_r_jet=h_over_r;

    Bpole=1.0; 
    //    Omegastar=0.0216; // in units of per GM/c^3
    Omegastar=0.0013;
    //    Omegastar=0.0433248; // in units of per GM/c^3

  }

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
  dt=1E-5;

  return(0);
}

int init_grid(void)
{
  SFTYPE rh;
  
  

  if(0){ // pulsar ffde
    defcoord = 10;
    a = 0.20 ;
    Rin=4.84; // in GM/c^2, NS surface

    //R0 = 0.0;
    //  R0=0.0;
    // see ns-boundarylayer.nb
    R0 = 0.903*Rin; // accounts for boundary layer with 20 zones above surface for N1=256 Rout=100GM/c^2

    Rhor=rhor_calc(0);
    Rout = 400.;
  }
  else if(1){ // ns MHD 1

    // metric stuff first
    a = 0.20 ;

    defcoord = 10;
    Rin=4.84; // in GM/c^2, NS surface

    //R0 = 0.0;
    //  R0=0.0;
    // see ns-boundarylayer.nb
    //  R0 = 0.903*Rin; // accounts for boundary layer with 20 zones above surface for N1=256 Rout=100GM/c^2
    //  R0 = 0.786*Rin;
    //    R0 = 0.5*Rin; // so midpoint of grid is lightcylinder
    R0 = 0.3*Rin;

    Rhor=rhor_calc(0);
    Rout = 1E4;
 
  }
  else if(0){ // ns MHD 2

    // metric stuff first
    a = 0.20 ;

    defcoord = 10;
    Rin=4.84; // in GM/c^2, NS surface

    //R0 = 0.0;
    //  R0=0.0;
    // see ns-boundarylayer.nb
    //  R0 = 0.903*Rin; // accounts for boundary layer with 20 zones above surface for N1=256 Rout=100GM/c^2
    //  R0 = 0.786*Rin;
    R0 = 0.5*Rin; // so midpoint of grid is lightcylinder

    Rhor=rhor_calc(0);
    Rout = 1E6;
 
  }


  hslope = 1.04*pow(h_over_r,2.0/3.0);


  return(0);
}

int init_global(void)
{

  if(0){// pulsar ffde
    GAMMAMAX=1000.0;
    GAMMAFAIL=100.0*GAMMAMAX;
    BSQORHOLIMIT=20.0;
    BSQOULIMIT=100.0;
    GAMMADAMP=5.0;

    BCtype[X1UP]=OUTFLOW;


  }
  else if(1){ // NS MHD 1
    GAMMAMAX=50.0;
    GAMMAFAIL=100.0*GAMMAMAX;
    BSQORHOLIMIT=20.0;
    BSQOULIMIT=100.0;
    GAMMADAMP=5.0;

    rescaletype=1;
    //    RHOMIN=1E-10;
    //    UUMIN=1E-11;

    RHOMIN=1E-2;
    UUMIN=1E-4;

    tf = 5E4;

    //    BCtype[X1UP]=FIXEDOUTFLOW;
    BCtype[X1UP]=OUTFLOW;



  }
  else if(0){ // NS MHD 2
  
    GAMMAMAX=25.0;
    GAMMAFAIL=100.0*GAMMAMAX;
    BSQORHOLIMIT=20.0;
    BSQOULIMIT=100.0;
    GAMMADAMP=5.0;

    // set to get sigma=10^6 like suspected in crab near NS
    // ( i.e. b^2/rho=10^6 with b=1 means RHOMIN=1E-6)
    RHOMIN=1.e-6;
    UUMIN =1.e-8;
    RHOMINLIMIT=1.e-20;
    UUMINLIMIT =1.e-20;

    rescaletype=0;

    tf = 50000.0;

    BCtype[X1UP]=FIXEDOUTFLOW;

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
  avgscheme=DONOR;
  
  do_transverse_flux_integration = 1;
  do_source_integration = 1;
  do_conserved_integration = 1;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;


#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  //lim = WENO5BND;
  //lim = PARA;
  TIMEORDER=4;
  //  lim = DONOR;
  lim = MC;
  //  TIMEORDER=2;
  //  DOENOFLUX = ENOFINITEVOLUME;
  DOENOFLUX = NOENOFLUX;
  //fluxmethod=HLLFLUX;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  //FLUXB = FLUXCTHLL;
  //  UTOPRIMVERSION=UTOPRIM5D1;

#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim = MC;
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

  //BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=NSSURFACE;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
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



  /* output choices */

  DTd = 50.0;   /* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2.0;   /* logfile frequency, in units of M */
  DTi = 2.0;   /* image file frequ., in units of M */
  DTdebug = DTd; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;   /* restart file period in steps */

  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  set_atmosphere(0,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  if(pr[RHO]<pratm[RHO]){
    PLOOP(pliter,pl) pr[pl]=pratm[pl];
  }
  

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
  Risco=rmso_calc(PROGRADERISCO);

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




  // copy over initial solution as analytic solution
  // SET ANALYTIC SOLUTION FROM vector potential-based solution
  // NEEDED FOR BOUND in case uses panalytic
  COMPZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND, -N3BND, N3-1+N3BND) PLOOP(pliter,pl){
    MACP0A1(panalytic,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
  }


  pdump=panalytic;
  // dump analytic solution
  if (dump(9995) >= 1){
    dualfprintf(fail_file,"unable to print dump file\n");
    return (1);
  }




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




  return(0);


}


#define KEPDISK 0
#define DONUTDISK 1
#define NODISK 2

#define DISKTYPE NODISK

// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr)
{
  int init_donut(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr);
  int init_kepdisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr);
  int init_nodisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr);


  if(DISKTYPE==KEPDISK){
    return(init_kepdisk(whichvel,whichcoord,i,j,k,pr));
  }
  else if(DISKTYPE==DONUTDISK){
    return(init_donut(whichvel,whichcoord,i,j,k,pr));
  }
  else if(DISKTYPE==NODISK){
    return(init_nodisk(whichvel,whichcoord,i,j,k,pr));
  }

  return(1); // shouldn't reach here
}



// as in fishmon
int init_donut(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeom,geom;
  

  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];


  rin = 6. ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  kappa = 1.e-3 ;
  beta = 1.e2 ;
  randfact = 4.e-2;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];



  sth = sin(th);
  cth = cos(th);

  /* calculate lnh */
  DD = r * r - 2. * r + a * a;
  AA = (r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth;
  SS = r * r + a * a * cth * cth;
  
  thin = M_PI / 2.;
  sthin = sin(thin);
  cthin = cos(thin);
  DDin = rin * rin - 2. * rin + a * a;
  AAin = (rin * rin + a * a) * (rin * rin + a * a)
    - DDin * a * a * sthin * sthin;
  SSin = rin * rin + a * a * cthin * cthin;
  
  if (r >= rin) {
    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
                               (AA * sth * AA * sth))) / (SS * DD /
                                                          AA))
      - 0.5 * sqrt(1. +
                   4. * (l * l * SS * SS) * DD / (AA * AA * sth *
                                                  sth))
      - 2. * a * r * l / AA -
      (0.5 *
       log((1. +
            sqrt(1. +
                 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin *
                                                      sthin *
                                                      sthin))) /
           (SSin * DDin / AAin))
       - 0.5 * sqrt(1. +
                    4. * (l * l * SSin * SSin) * DDin / (AAin *
                                                         AAin *
                                                         sthin *
                                                         sthin))
       - 2. * a * rin * l / AAin);
  } else
    lnh = 1.;
  

  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {


    get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {
    hm1 = exp(lnh) - 1.;
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
    u = kappa * pow(rho, gam) / (gam - 1.);
    ur = 0.;
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}




int init_kepdisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  FTYPE V[NDIM];
  FTYPE nz_func(FTYPE R) ;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];

  beta=1.e2;
  randfact=1E-3;
  //rin = (1. + h_over_r)*Risco;
  //rin=Risco;
  rin = 6;



  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    rho = 1.e-7*RHOMIN ;
    u = 1.e-7*UUMIN ;

    get_geometry(i, j, k, CENT, &geom); // true coordinate system

    // normal observer velocity
    if(*whichvel==VEL4){
      ur = -geom.gcon[GIND(0,1)]/sqrt(-geom.gcon[GIND(0,0)]) ;
      uh = -geom.gcon[GIND(0,2)]/sqrt(-geom.gcon[GIND(0,0)]) ;
      up = -geom.gcon[GIND(0,3)]/sqrt(-geom.gcon[GIND(0,0)]) ;
    }
    else if(*whichvel==VEL3){
      ur = geom.gcon[GIND(0,1)]/geom.gcon[GIND(0,0)] ;
      uh = geom.gcon[GIND(0,2)]/geom.gcon[GIND(0,0)] ;
      up = geom.gcon[GIND(0,3)]/geom.gcon[GIND(0,0)] ;
    }
    else if(*whichvel==VELREL4){
      ur = 0.0;
      uh = 0.0;
      up = 0.0;
    }

    pr[RHO] = rho;
    pr[UU] = u;
    pr[U1] = ur;
    pr[U2] = uh;
    pr[U3] = up;

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {

    H = h_over_r*R ;
    nz = nz_func(R) ;
    z = r*cos(th) ;
    S = 1./(H*H*nz) ;
    cs = H*nz ;

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H))
      * taper_func(R,rin) ;
    u = rho*cs*cs/(gam - 1.) ;
    ur = 0. ;
    uh = 0. ;
    up = 1./(pow(r,1.5) + a) ;
    // solution for 3-vel


    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;
    
    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}



int init_nodisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom geom;


  get_geometry(i, j, k, CENT, &geom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,&geom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system


  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}


// assumes we are fed the true densities
int normalize_densities(FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;

#if(DISKTYPE!=NODISK)
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
  SFTYPE rho_av, q;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom geom;
  FTYPE fieldhor;
  FTYPE nz_func(FTYPE R) ;
  FTYPE nz,z,S,cs,rho,u_ref,u_av;
  FTYPE dxdxp[NDIM][NDIM];



  fieldhor=3.0*h_over_r;




  *A=0;

  if(l==3){// A_\phi
    
    get_geometry(i, j, k, CORN3, &geom);
    coord(i, j, k, CORN3, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    dxdxprim(X, V, dxdxp);



    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      *A += 0.5*r*sin(th) ;
    }
    // dipole field version
    // see ns_dipole.nb and Rezzolla, Ahmedov, Miller (2001) eqs 54-56
    if((FIELDTYPE==DIPOLEFIELD)||(FIELDTYPE==DISKDIPOLE)){
      *A += 0.5*Bpole*Rin*Rin*Rin*sin(th)*sin(th)/r ; // assumes Rin is surface of NS
    }
    /* field-in-disk version */
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)||(FIELDTYPE==DISKDIPOLE)){
      if(DISKTYPE==DONUTDISK){ // single loop
        // average of density that lives on CORN3
 
 
        // since init_vpot() is called for all i,j,k, can't use
        // non-existence values, so limit averaging:
        if((i==-N1BND)&&(j==-N2BND)){
          rho_av = MACP0A1(p,i,j,k,RHO);
          u_av = MACP0A1(p,i,j,k,UU);
        }
        else if(i==-N1BND){
          rho_av = AVGN_2(p,i,j,k,RHO);
          u_av = AVGN_2(p,i,j,k,UU);
        }
        else if(j==-N2BND){
          rho_av = AVGN_1(p,i,j,k,RHO);
          u_av = AVGN_1(p,i,j,k,UU);
        }
        else{ // normal cells
          rho_av = AVGN_for3(p,i,j,k,RHO);
          u_av = AVGN_for3(p,i,j,k,UU);
        }
 
        q = rho_av / rhomax - 0.2;
 
        // if (q > 0.)      *A += q;
 
        //if (q > 0.)      *A += q*sqrt(2.0*(gam-1.0)*u_av/beta);
        if (q > 0.)      *A += q*sqrt(2.0*10.0/beta);


      }
      else if(DISKTYPE==KEPDISK){ // multiple loops
        // in MPI, don't have "equatorial value" for all cpus so easily, so just compute it
        // just replaced R with r and remaining th with Pi/2
        H = h_over_r*r ;
        nz = nz_func(r) ;
        z = 0.0 ;
        S = 1./(H*H*nz) ;
        cs = H*nz ;
 
        rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H))* taper_func(r,rin) ;
        u_ref = rho*cs*cs/(gam - 1.) ;
        // u_ref/=rhomaxold; // SUPERGODMARK?
        /*
          u_ref = 0.25*(
          MAC(p,i,N2/2,UU) +
          MAC(p,i-1,N2/2,UU) +
          MAC(p,i,N2/2-1,UU) +
          MAC(p,i-1,N2/2-1,UU)) ;
        */
#define STARTFIELD (1.1*rin)
 
        if(r > STARTFIELD) q = ((u_av/u_ref) - 0.2)*pow(r,0.25) ;
        else q = 0. ;
 
        if(q > 0.){
          *A += q*q*sin(log(r/STARTFIELD)/fieldhor)*taper_func(r,STARTFIELD);
          //trifprintf("%d %d u_ref=%g A=%g\n",i,j,u_ref,A[i][j]);
        }
      }
    }

    // assume we were setting A_\phi and not A_{x^3}, so need to correct
    // correct for 2\pi factor
    *A *=dxdxp[PH][PH];

  }


  

  return(0);

}


int init_vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  extern int vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*p)[NSTORE2][NSTORE3][NPR]);


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
