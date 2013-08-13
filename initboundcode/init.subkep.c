
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0  /* reduce u_phi by this amount */
#define MAXPASSPARMS 10


static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;


static FTYPE nz_func(FTYPE R) ;



int prepre_init_specific_init(void)
{
  int funreturn;
  
  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);

}


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.3;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  rhodisk=1.0;

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  return(0);
}

int set_fieldfrompotential(int *fieldfrompotential)
{
  int pl,pliter;

  // default (assume all fields are from potential)
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;


  return(0);
}


int init_conservatives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*Utemp)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  int fieldfrompotential[NDIM];

  set_fieldfrompotential(fieldfrompotential);

  funreturn=user1_init_conservatives(fieldfrompotential, prim,pstag, Utemp, U);
  if(funreturn!=0) return(funreturn);


  return(0);
 
}


int post_init_specific_init(void)
{
  int funreturn;

  funreturn=user1_post_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}


#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define SUBKEP 3

// For blandford problem also need to set:
// 0) WHICHPROBLEM 2
// 1) a=0.92
// 2) Rout=1E3
// 3) hslope=0.3
// 4) BSQORHOLIMIT=1E2;
// 5) BSQOULIMIT=1E3;
// 6) UORHOLIMIT=1E3;
// 7) tf=1E4; // and maybe DTd's
// 8) randfact = .2;
// 8.5) Choose fieldtype: FIELDTYPE BLANDFORDQUAD or DISKFIELD (SS dipole)
// 9) lim=PARALINE; FLUXB=FLUXCTSTAG; TIMEORDER=4;
// 10) N1,N2,N3 in init.h
// 11) MAXBND 4
// 12) PRODUCTION 1
// 13) USE<systemtype>=1
// 14) setup batch

/*

  Models to run:

  Constant parameters:

  1) Rout=1E3 and run for tf=1E4 (so will take 5X longer than compared to Orange run at 128x128x32)

  2) BSQORHOLIMIT=1E3, etc.

  3) PARALINE, FLUXCTSTAG, TO=4

  4) Form of A_\phi fixed

  Field parameter studies in 2D axisymmetry at 256^2:

  1) H/R=0.3, a=0.9: LS quadrapole,  LS dipole, SS quadrapole, SS dipole

 

  Spin parameter study in 2D axisymmetry at 256^2:

 

  1) H/R=0.3, LS quadrapole: a=-.999,-.99,-.9,-.5,-0.2,0,.2,.5,.9,.99,.999

  H/R parameter study in 2D axisymmetry at 256^2:

  1) a=0.9 LS quadrapole with H/R=0.1,0.3,0.9,1.5

  2D Fiducial Models:

  1) Using a=0.9, H/R=0.3, LS quad and LS dipole, do two 2D fudicial models at: 1024^2

  3D Fiducial Models:

  1) Using a=0.9, H/R=0.3, LS quadrapole and LS dipole, do two 3D fiducial models at 2 different resolutions: 128x128x32 and 256x256x64

  Questions for Roger:

  1) Choice for disk thickness?
  2) Choice for field shape -- specifically?
  3) Choice for flux threading disk vs. BH initially?
  4) Ask about BZ77 and residual A_\phi at pole
  5) 

*/


#define WHICHPROBLEM SUBKEP // choice


int init_defcoord(void)
{
  
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // define coordinate type
  defcoord = JET3COORDS;
#elif(WHICHPROBLEM==GRBJET)
  // define coordinate type
  defcoord = JET4COORDS;
#elif(WHICHPROBLEM==SUBKEP)
  // define coordinate type
  defcoord = LOGRSINTH;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first
  //  a = 0.9375 ;
  a = 0.0;
  

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==SUBKEP)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  Rout = 40.0;
#elif(WHICHPROBLEM==GRBJET)
  R0 = -3.0;
  Rout = 1E5;
#endif

 
  Rhor=rhor_calc(0);

  //  hslope = 0.3;
  // hslope = 1.04*pow(h_over_r,2.0/3.0);
  hslope = 1.0;

  setRin_withchecks(&Rin);





  return(0);
}



int init_global(void)
{
  int pl,pliter;
  int funreturn;


  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);


  //////////////////
  // overrides for more detailed problem dependence


  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  FLUXB=FLUXCTTOTH;

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=1E3;
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
#elif(WHICHPROBLEM==GRBJET)
  BCtype[X1UP]=FIXEDOUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E3;
  BSQOULIMIT=1E4;
  RHOMIN = 23.0;
  UUMIN = 1.7;
#elif(WHICHPROBLEM==SUBKEP)
  BCtype[X1UP]=FIXEDOUTFLOW;
  //BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E2;
  BSQOULIMIT=1E3;
  UORHOLIMIT=1E3;
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
#endif






#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==SUBKEP)
  /* output choices */
  tf = 2000.0;

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 50.;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 50.0;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 2.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 2.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 50.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 3000;

#elif(WHICHPROBLEM==GRBJET)
  /* output choices */
  tf = 5E5;
  
  DTd = 250.;                 /* dumping frequency, in units of M */
  DTavg = 250.0;
  DTener = 2.0;                       /* logfile frequency, in units of M */
  DTi = 10.0;                 /* image file frequ., in units of M */
  DTdebug = 250.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 1000;                  /* restart file period in steps */
#endif



  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int funreturn;

  //  funreturn=user1_init_atmosphere(whichvel, whichcoord,i, j, k, pr);
  //  if(funreturn!=0) return(funreturn);

  return(0);

}

int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);



  beta = 1.e2 ;
  randfact = 4.e-2;

#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET)
#elif(WHICHPROBLEM==SUBKEP)
#endif

  



  //SASMARK restart: need to populate panalytic with IC's
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
    MYFUN(init_primitives(panalytic,pstaganalytic,GLOBALPOINT(utemparray),vpotanalytic,Bhatanalytic,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }




  
  // diagnostic
  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    i=INFULL1;
    j=k=0;
    coord(i,j,k, FACE1, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    trifprintf("rmin(i=%d,X=%21.15g): %21.15g\n", i,X[1],r);
    trifprintf("rmin/rh: %21.15g\n", r / Rhor );
    //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
    if(r/Rhor<=1.0){
      trifprintf("inner grid is inside horizon\n");
    }
    else{
      trifprintf("inner grid is outside horizon\n");
    }
  }


  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;

  funreturn=user1_init_primitives(prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

  return(0);


}



int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_subkep(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==SUBKEP)
  return(init_dsandvels_subkep(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#endif

}


// unnormalized density
int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeomdontuse;
  struct of_geom *ptrrealgeom=&realgeomdontuse;
  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];
  int pl,pliter;






  kappa = 1.e-3 ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  

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


    get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrrealgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

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

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}



int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th,ph;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  int pl,pliter;


  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  ph=V[3];




  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    get_geometry(i, j, k, CENT, ptrgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

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

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H)) * taper_func(R,rin) ;
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

    if(FLUXB==FLUXCTSTAG){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}

int init_dsandvels_subkep(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th,ph;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  int pl,pliter;


  dt=1E-5;

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  ph=V[3];




  /* region outside disk */
  R = r*sin(th) ;

  if(1||r < 6.0) {

    // assumes VELREL4
    pr[U1]=0;
    pr[U2]=0;
    pr[U3]=0;

    if(FLUXB==FLUXCTSTAG){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
  }
  else {

    // this gets BL-geometry in native BL spc coordinates (not PRIMECOORDS) 
    struct of_geom geombl; 
    gset(0,BLCOORDS,i,j,k,&geombl);

    pr[U1] = (geombl.gcon[GIND(0,1)])/(geombl.gcon[GIND(0,0)]) ;
    pr[U2] = (geombl.gcon[GIND(0,2)])/(geombl.gcon[GIND(0,0)]) ;
    pr[U3] = (geombl.gcon[GIND(0,3)])/(geombl.gcon[GIND(0,0)]) ;


    //  pr[U1] = -1e-2*pow(r,-1./2.) ;
    // pr[U2] = 0.0 ;    
    // pr[U3] = 0.0;

    if(FLUXB==FLUXCTSTAG){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    // DEBUG:
    //    PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d j=%d k=%d pl=%d pr=%21.15g\n",i,j,k,pl,pr[pl]);

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;

  }


  // set densities
  pr[RHO] = RHOMIN*pow(r,-2.0) ;
  pr[UU]  = UUMIN*pow(r,-5./2.) ;


  return(0);
}








#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3
#define SUBKEP1 4

#define FIELDTYPE SUBKEP1


FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower)
{
  FTYPE fneg,fpos;
  FTYPE gpara;

  fneg=1.0-pow(cos(th),thpower);
  fpos=1.0+pow(cos(th),thpower);
  gpara=0.5*(myr*fneg + 2.0*fpos*(1.0-log(fpos)));
  // remove BZ77 Paraboloidal divb!=0 at pole
  gpara=gpara-2.0*(1.0-log(2.0));

  return(gpara);


}

FTYPE setblandfordfield(FTYPE r, FTYPE th)
{
  FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower);
  FTYPE rshift,myr,rpower,myz,myR,myvert;
  FTYPE thother,thpower,gparalow,gparahigh,mygpara;
  FTYPE aphi;


  rshift=4.0;
  rpower=0.75;
  thpower=4.0;
  

  myr=pow(r+rshift,rpower);
  myz=myr*cos(th);
  myR=myr*sin(th);
  myvert = (th>M_PI*0.5) ? (myr*sin(th)) : (myr*sin(-th));

  thother=M_PI-th;
  gparalow=setgpara(myr,th,thpower);
  gparahigh=setgpara(myr,thother,thpower);
  mygpara=(th<0.5*M_PI) ? gparalow : gparahigh;

  // GOOD:
  // aphi=mygpara;
  // aphi=mygpara*cos(th); // B1 diverges at pole
  //  aphi=mygpara*cos(th)*sin(th); // doesn't diverge as much
  //  aphi=mygpara*cos(th)*sin(th)*sin(th); // old choice before subtracted original BZ77 problem
  aphi=mygpara*cos(th); // latest choice
  //aphi=myvert*cos(th); // vert with quad
  //aphi=myR*cos(th);

  // BAD:
  // aphi=myvert;
  
  

  return(aphi);


}

// assumes normal field in pr
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE r,th;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);




  vpot=0.0;


  if(l==3){// A_\phi

    r=V[1];
    th=V[2];


    if(FIELDTYPE==SUBKEP1){
      vpot += 0.0;
    }

    // Blandford quadrapole field version
    if(FIELDTYPE==BLANDFORDQUAD){
      vpot += setblandfordfield(r,th);
    }

    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      FTYPE rpow;
      rpow=3.0/4.0; // Using rpow=1 leads to quite strong field at large radius, and for standard atmosphere will lead to \sigma large at all radii, which is very difficult to deal with -- especially with grid sectioning where outer moving wall keeps opening up highly magnetized region
      vpot += 0.5*pow(r,rpow)*sin(th)*sin(th) ;
    }


    /* field-in-disk version */
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){
      // average of density that lives on CORN3


      // since init_vpot() is called for all i,j,k, can't use
      // non-existence values, so limit averaging:
      if((i==-N1BND)&&(j==-N2BND)){
        rho_av = MACP0A1(prim,i,j,k,RHO);
      }
      else if(i==-N1BND){
        rho_av = AVGN_2(prim,i,j,k,RHO);
      }
      else if(j==-N2BND){
        rho_av = AVGN_1(prim,i,j,k,RHO);
      }
      else{ // normal cells
        rho_av = AVGN_for3(prim,i,j,k,RHO);
      }

      q = rho_av / rhomax - 0.2;

      if (q > 0.)      vpot += q;
    }
  }

  //////////////////////////////////
  //
  // finally assign what's returned
  //
  //////////////////////////////////
  *A = vpot;
  *whichcoord = MCOORD;



  return(0);

}



int init_vpot2field_user(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  int fieldfrompotential[NDIM];

  funreturn=user1_init_vpot2field_user(fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
 

  return(0);


}



// assumes we are fed the true densities
int normalize_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  FTYPE parms[MAXPASSPARMS];
  int eqline;

 

  return(0);
}



// assumes we are fed the true densities
int getmax_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax)
{
  int funreturn;

  funreturn=user1_getmax_densities(prim,rhomax, umax);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}



// get maximum b^2 and p_g
int get_maxes(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max)
{
  int funreturn;
  int eqslice;
  FTYPE parms[MAXPASSPARMS];
  
  if(FIELDTYPE==VERTFIELD || FIELDTYPE==BLANDFORDQUAD){
    eqslice=1;
  }
  else{
    eqslice=0;
  }

  parms[0]=rin;

  funreturn=user1_get_maxes(eqslice, parms,prim, bsq_max, pg_max);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;

 
  funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
  if(funreturn!=0) return(funreturn);
 
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
                                                                     pow(r,
                                                                         1.5)
                                                                     -
                                                                     3.
                                                                     /
                                                                     r)))
          / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
             (pow(a, 2) + (-2. + r) * r))
          );
}


// UUMIN/RHOMIN used for atmosphere

// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int set_atmosphere(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int pl,pliter;
  int i,j,k;



  int atmospheretype=1;



  bl_coord_ijk_2(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, V);
  r=V[1];
  th=V[2];


  if(r<6.0){
    int funreturn;
    funreturn=user1_set_atmosphere(atmospheretype, whichcond, whichvel, ptrgeom, pr);
    if(funreturn!=0) return(funreturn);
  }
  else{

    // Bondi like initial atmosphere
    //    rho = RHOMIN * 1.E-14;
    //    u = UUMIN * 1.E-14;

    // default
    PLOOP(pliter,pl) prlocal[pl]=pr[pl];

    if((EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD)){
      // Bondi-like atmosphere
      if(rescaletype==4){
        if(atmospheretype==1){
          // couple rescaletype to atmosphere type
          prlocal[RHO] = RHOMIN*pow(r,-2.0);
        }
        else if(atmospheretype==2){
          // couple rescaletype to atmosphere type
          if(r>40.0) prlocal[RHO] = RHOMIN*pow(r,-2.0);
          else prlocal[RHO] = RHOMIN*pow(40.0,-2.0);
        }
      }
      else{
        prlocal[RHO] = RHOMIN*pow(r,-1.5);
      }
    }
    else{
      prlocal[RHO] = 0;
    }


    if(EOMTYPE==EOMGRMHD){
      // Bondi-like atmosphere
      prlocal[UU]  = UUMIN*pow(r,-2.5);
    }
    else{
      prlocal[UU]  = 0;
    }

    // this gets BL-geometry in native BL spc coordinates (not PRIMECOORDS)                                                                    

    i=ptrgeom->i;
    j=ptrgeom->j;
    k=ptrgeom->k;
    struct of_geom geombl;
    FTYPE prblcoord[NPR];
    //  FTYPE pstag[NPR];

    gset(0,BLCOORDS,i,j,k,&geombl);

    //    prblcoord[U1] = (geombl.gcon[GIND(0,1)])/(geombl.gcon[GIND(0,0)]) ;
    prblcoord[U1] = UUMIN/RHOMIN/10.0*pow(r,-0.5);
    prblcoord[U2] = (geombl.gcon[GIND(0,2)])/(geombl.gcon[GIND(0,0)]) ;
    prblcoord[U3] = (geombl.gcon[GIND(0,3)])/(geombl.gcon[GIND(0,0)]) ;

    prblcoord[B1] = 0.0;
    prblcoord[B2] = 0.0;
    prblcoord[B3] = 0.0;

    // PLOOP(pliter,pl) pstag[pl]=prblcoord[pl];  
    // transform BLCOORDS to PRIMECOORDS:  
    //  MYFUN(transform_primitive_vB(VEL3, BLCOORDS, i,j,k, prblcoord, pstag),"init.c:init_primitives","transform_primitive_vB()",0);
    if (bl2met2metp2v(VEL3,BLCOORDS,prblcoord, i,j,k) >= 1) FAILSTATEMENT("initbase.c:transform_primitive_vB()", "bl2ks2ksp2v()", 1);

  
    prlocal[U1]=prblcoord[U1];
    prlocal[U2]=prblcoord[U2];
    prlocal[U3]=prblcoord[U3];

    if(whichcond==1){
      if(100.0*dt>TAUADJUSTATM){
        dualfprintf(fail_file,"dt=%21.15g and TAUADJUSTATM=%21.15g\n",dt,TAUADJUSTATM);
        myexit(1);
      }
      // TAUADJUSTATM must be >> dt always in order for this to make sense (i.e. critical damping to fixed solution)
      PLOOP(pliter,pl) pr[pl] = pr[pl]+(prlocal[pl]-pr[pl])*dt/TAUADJUSTATM;
    }
    else if(whichcond==0){ 
      PLOOP(pliter,pl) pr[pl] = prlocal[pl];
      // very specific
      // always outflow field
      //    pr[B1] = pr[B2] = pr[B3] = 0;
    }
    else if(whichcond==-1){ 
      // t=0, just sets as "floor"
      PLOOP(pliter,pl) pr[pl] = prlocal[pl];
      pr[B1] = pr[B2] = pr[B3] = 0;
    }

    return(0);
  }// end if beyond 6M
 


  return(0);
}




int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  int funreturn;
  
  funreturn=set_density_floors_default(ptrgeom, pr, prfloor);
  if(funreturn!=0) return(funreturn);

  return(0);
}




static FTYPE nz_func(FTYPE R)
{
  return(
         sqrt(
              (3.*a*a - 4.*a*sqrt(R) + R*R)/
              pow(R*(a + pow(R,1.5)),2)
              )
         ) ;


}



// Setup problem-dependent grid sectioning
int theproblem_set_enerregiondef(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // Torus problem
  //  torus_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);
  //jet_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);

#if(1)
  enerregiondef[POINTDOWN][1]=0;
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
#endif

  return(0);
}


int theproblem_set_enerregionupdate(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps)
{

  ////////////////////
  //
  // Setup update period
  //
  ///////////////////

  ////////
  //
  // number of steps after which position/size of active section is updated
  //
  ////////
#if(N3==1)
  //  *updateeverynumsteps=100;
  *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
#else
  //  *updateeverynumsteps=10;
  *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
#endif

  ////////
  //
  //number of steps after which position/size of active section is reported to file
  //
  ////////
  *everynumsteps = *updateeverynumsteps*100;

  return(0);
}


// specify MPI task rank ordering
// example user-dependent code
int theproblem_set_myid(void)
{
  int retval;
 
  // default is to do nothing
  //  retval=jet_set_myid();
  retval=0;

  // do other things?

  return(retval);

}


