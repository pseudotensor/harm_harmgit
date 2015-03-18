// uses init.fishmon.h and bounds.fishmon.c
// See also coord.c, init.tools.c, and restart.c if changes to header/dumpfile
// For thickdisk7->fishmon(new):
// a) must change TRACKVPOT->0.  If doing zakamskabig,then use some large-box versions of parameters given by WHICHPROBLEM==THICKDISK case.
// b) Also, if restarting from old zakamskabig, must truncate read restart of normal dump file, (upper pole file should zero out values as desired), NUMFAILFLOORFLAGS-> (NUMFAILFLOORFLAGS-2) upon read, NUMDUMPTYPES->NUMDUMPTYPES-1
// c) Also, set makehead.inc USETACCLONESTAR4->1


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

#define USER_THETAROTMETRIC (0.0)
#define USER_THETAROTPRIMITIVES (0.0) // probably want to choose 0, so initial conditions are as if no tilt

#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define THICKDISK 3 // very similar to NORMALTORUS


// For blandford problem also need to set:
// 0) WHICHPROBLEM
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


#define WHICHPROBLEM NORMALTORUS //THICKDISK // NORMALTORUS // choice


static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;


static FTYPE nz_func(FTYPE R) ;

FTYPE normglobal;
int inittypeglobal; // for bounds to communicate detail of what doing



int prepre_init_specific_init(void)
{
  int funreturn;

  
  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);

}


int pre_init_specific_init(void)
{

  if(WHICHPROBLEM==THICKDISK){
    // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
    h_over_r=1.0;
    // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
    h_over_r_jet=M_PI*0.4;
  }
  else{
    // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
    h_over_r=0.3;
    // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
    h_over_r_jet=2.0*h_over_r;
  }

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


  if(WHICHPROBLEM==THICKDISK){
    cour=0.8;
    //  fluxmethod= HLLFLUX;
  }
  else{
    // leave as default
  }

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}



int init_defcoord(void)
{
  
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // define coordinate type
  defcoord = JET3COORDS;
#elif(WHICHPROBLEM==THICKDISK)
  defcoord = JET6COORDS;
#elif(WHICHPROBLEM==GRBJET)
  // define coordinate type
  defcoord = JET4COORDS;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first
  a = 0.9375 ;


  

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  Rout = 40.0;
#elif(WHICHPROBLEM==THICKDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  //  Rout = 1E3;
  Rout = 1.3*2*1E4;
#elif(WHICHPROBLEM==GRBJET)
  R0 = -3.0;
  Rout = 1E5;
#endif

 
  Rhor=rhor_calc(0);

  //  hslope = 0.3;
  hslope = 1.04*pow(h_over_r,2.0/3.0);


  setRin_withchecks(&Rin);



  // set global THETAROTPRIMITIVES
  if(ALLOWMETRICROT){
    THETAROTPRIMITIVES=USER_THETAROTPRIMITIVES; // 0 to M_PI : what thetarot to use when primitives are set
  }
  else{
    THETAROTPRIMITIVES=0.0; // DO NOT CHANGE
  }

  if(ALLOWMETRICROT){
    THETAROTMETRIC = USER_THETAROTMETRIC; // defines metric generally
  }
  else{
    THETAROTMETRIC = 0.0;
  }


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
  //  FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E5; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=1E10;
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
#elif(WHICHPROBLEM==THICKDISK)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=50.0; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=1E10;
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
#endif





  // default dumping period
  int idt;
  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=50.0;

  // ener period
  DTdumpgen[ENERDUMPTYPE] = 2.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 2.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 1000;
  DTfake=MAX(1,DTr/10);


#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  /* output choices */
  tf = 2000.0;
#elif(WHICHPROBLEM==THICKDISK)
  /* output choices */
  tf = 1.3E4*2.0;


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
  DTfake=MAX(1,DTr/10);
#endif



  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int funreturn;

  funreturn=user1_init_atmosphere(whichvel, whichcoord,i, j, k, pr);
  if(funreturn!=0) return(funreturn);

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


  // defaults
  beta = 1.e2 ;
  randfact = 4.e-2;



#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
#elif(WHICHPROBLEM==THICKDISK)
  //  beta = 1.e2 ;
  //  beta = 20.0;
  beta = 10.0*4.0; // 256x128x256
  //  randfact = 4.e-2;
  randfact = 0.1;
  //rin = Risco;
  //  rin = 15.0 ;
  rin = 10.0 ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET)
#endif



  if(RESTARTMODE==1){
    dualfprintf(fail_file,"WARNING: On lonestar4 with zakamskabig restart (first restart attempt), reaching inside below conditional eventually leads to core dump on several tasks.  Some core dumps are good, and show i index inside transform_primitive_pstag() called by transform_primitive_vB called in user1_init_primitives() for the first exposed loop -- leads to crazy i index (very negative or positive).  So memory segfault occurs.  No idea where that i is being changed.  Maybe memory leak, but all cores that had good core dump looked ok.  So maybe something odd going on.  Note pglobal and panalytic and other arrays have funny or NULL-like address, but aren't used so probably just an optimization thing.\n");
    

    /*
      (gdb) bt
      #0  transform_primitive_pstag (whichvel=0, whichcoord=1, i=-1534771426, j=-4, k=11,
      p=0x130b79f0, pstag=0x130ba300) at fluxvpot.c:1196
      #1  0x0000000000408989 in transform_primitive_vB (whichvel=0, whichcoord=1, i=-1534771426,
      j=-4, k=11, p=0x130b79f0, pstag=0x130ba300) at initbase.c:2718
      #2  0x0000000000405826 in user1_init_primitives (prim=0x0, pstag=0x1,
      ucons=0x3f385eb2a4853f1e, vpot=0xfffffffffffffffc, Bhat=0xb, panalytic=0x130b79f0,
      pstaganalytic=0x130ba300, vpotanalytic=0x0, Bhatanalytic=0x0, F1=0x12af5968,
      F2=0x1b22a28, F3=0x1fe8028, Atemp=0x797e380) at init.tools.c:340
      #3  0x000000000040466e in init_grid_post_set_grid (prim=0x0, pstag=0x1,
      ucons=0x3f385eb2a4853f1e, vpot=0xfffffffffffffffc, Bhat=0xb, panalytic=0x130b79f0,
      pstaganalytic=0x130ba300, vpotanalytic=0x0, Bhatanalytic=0x0, F1=0x12af5968,
      F2=0x1b22a28, F3=0x1fe8028, Atemp=0x797e380) at init.c:367
      #4  0x000000000040b8bd in init (argc=0x0, argv=0x1) at initbase.c:143
      #5  0x00000000004bd96a in main (argc=6, argv=0x7fffd5b2ec08) at main.c:30
      (gdb) print myid
      $4 = 1380
      (gdb) print startpos
      $5 = {0, 136, 64, 112}
      (gdb) print mycpux1
      No symbol "mycpux1" in current context.
      (gdb) print numprocs
      $6 = 1536


    */

  }


  //SASMARK restart: need to populate panalytic with IC's
  if(0&& RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
    //    MYFUN(init_primitives(prim,pstag,ucons,vpot,Bhat,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);

    // utemparray only used otherwise in advance.c
    MYFUN(init_primitives(panalytic,pstaganalytic,GLOBALPOINT(utemparray),vpotanalytic,Bhatanalytic,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }

  trifprintf("BEGIN check_rmin\n");
  // check rmin
  check_rmin();
  trifprintf("END check_rmin\n");


  // check that singularities are properly represented by code
  trifprintf("BEGIN check_spc_singularities_user\n");
  // SUPERGODMARK: Goes very slowly sometimes randomly for unknown reasons.
  dualfprintf(fail_file,"WARNING: check_spc_singularities_user() oddly stalls sometimes...\n");
  check_spc_singularities_user();
  trifprintf("END check_spc_singularities_user\n");
  dualfprintf(fail_file,"WARNUNG: done with check_spc_singularities_user(), but it sometimes stalls or goes very very slow for no apparently good reason.  E.g., on NAUTILUS with -O0, very slow checks.  But just putting dualfprintf before and after the above call leads to quick finish.");

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;
  int inittype;
  FTYPE thetarotorig;


  thetarotorig=THETAROT;
  THETAROT = THETAROTPRIMITIVES; // define rho,u,v,B as if no rotation (but metric might still be used, so still use set_grid_all() in initbase.c)


  inittype=1;

  funreturn=user1_init_primitives(inittype, prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

  THETAROT = thetarotorig; // back to previous version

  return(0);


}



int init_dsandvels(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

  // assume inittype not used, pos==CENT, and time doesn't matter (e.g. only used at t=0)

#if(WHICHPROBLEM==NORMALTORUS||WHICHPROBLEM==THICKDISK)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
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



  // default
  kappa = 1.e-3 ;
  rmax = 12. ;

#if(WHICHPROBLEM==THICKDISK)
  kappa = 1.e-3 ;
  //  rmax = 1E2 ;
  rmax = 1E2 ;
  //  rmax = 60 ;
#endif


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








#define DISK1FIELD 0
#define DISK2FIELD 1
#define VERTFIELD 2
#define DISK1VERT 3
#define DISK2VERT 4
#define BLANDFORDQUAD 5
#define TOROIDALFIELD 6


#if(WHICHPROBLEM==THICKDISK)
//#define FIELDTYPE TOROIDALFIELD
#define FIELDTYPE DISK2FIELD
#else
#define FIELDTYPE DISK1FIELD
#endif




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
// SUPERNOTE: A_i must be computed consistently across all CPUs.  So, for example, cannot use randomization of vector potential here.
int init_vpot_user(int *whichcoord, int l, SFTYPE time, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av,q;
  FTYPE r,th,ph;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);


#define FRACAPHICUT 0.2
  //#define FRACAPHICUT 0.1



  vpot=0.0;



  // since init_vpot() is called for all i,j,k, can't use
  // non-existence values, so limit averaging:
  if((i==-N1BND)&&(j==-N2BND)){
    rho_av = MACP0A1(prim,i,j,k,RHO);
    u_av = MACP0A1(prim,i,j,k,UU);
  }
  else if(i==-N1BND){
    rho_av = AVGN_2(prim,i,j,k,RHO);
    u_av = AVGN_2(prim,i,j,k,UU);
  }
  else if(j==-N2BND){
    rho_av = AVGN_1(prim,i,j,k,RHO);
    u_av = AVGN_1(prim,i,j,k,UU);
  }
  else{ // normal cells
    rho_av = AVGN_for3(prim,i,j,k,RHO);
    u_av = AVGN_for3(prim,i,j,k,UU);
  }





  if(FIELDTYPE==TOROIDALFIELD){
    
    if(l==2){// A_\theta (MCOORD)
      
      r=V[1];
      th=V[2];

      //      q = r*r*r*fabs(sin(th)) * 1.0 ; // constant B^\phi
      q = r*r*r;

      q=q/(r); // makes more uniform in radius

      q = q*(u_av / umax - FRACAPHICUT); // weight by internal energy density
      //      q = (rho_av / rhomax - FRACAPHICUT);

      if(q<0.0) q=0.0;

      vpot += q;
      
    }
  }


  FTYPE rpow;
  rpow=3.0/4.0; // Using rpow=1 leads to quite strong field at large radius, and for standard atmosphere will lead to \sigma large at all radii, which is very difficult to deal with -- especially with grid sectioning where outer moving wall keeps opening up highly magnetized region
  //  FTYPE FIELDROT=M_PI*0.5;
  FTYPE FIELDROT=0.0;
  FTYPE hpow=2.0;


  if(l==2){// A_\theta

    r=V[1];
    th=V[2];
    ph=V[3];


    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISK1VERT)||(FIELDTYPE==DISK2VERT)){
      vpot += -(pow(r,rpow)*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph));
    }


  }

  if(l==3){// A_\phi

    r=V[1];
    th=V[2];
    ph=V[3];


    // Blandford quadrapole field version
    if(FIELDTYPE==BLANDFORDQUAD){
      vpot += setblandfordfield(r,th);
    }

    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISK1VERT)||(FIELDTYPE==DISK2VERT)){
      //vpot += 0.5*pow(r,rpow)*sin(th)*sin(th) ;
      vpot += pow(r,rpow)*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
    }


    /* field-in-disk version */
    if(FIELDTYPE==DISK1FIELD || FIELDTYPE==DISK1VERT){
      q = rho_av / rhomax - 0.2;
      if (q > 0.)      vpot += q;
    }


    if(FIELDTYPE==DISK2FIELD || FIELDTYPE==DISK2VERT){
      // average of density that lives on CORN3


#define FRACAPHICUT 0.2
      //#define FRACAPHICUT 0.1

      //      q = (rho_av / rhomax - FRACAPHICUT);
      q = (u_av / umax - FRACAPHICUT);

      //#define QPOWER 0.5
#define QPOWER (1.0)

#define POWERNU (2.0)
      //#define POWERNU (4.0)

      //      if (q > 0.)      vpot += q*q*pow(r*fabs(sin(th)),POWERNU);
      FTYPE fact1,fact2,SSS,TTT;
      fact1=pow(fabs(q),QPOWER)*pow(r*fabs(sin(th)),POWERNU);
      SSS=rin*0.5;
      TTT=0.28;
      fact2=sin(log(r/SSS)/TTT);

      if (q > 0.)      vpot += fact1*fact2;
      //      if (q > 0.)      vpot += q*q;

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



int init_vpot2field_user(SFTYPE time, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  int fieldfrompotential[NDIM];

  funreturn=user1_init_vpot2field_user(time, fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
 

  return(0);


}



// assumes we are fed the true densities
int normalize_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  FTYPE parms[MAXPASSPARMS];
  int eqline;

  eqline=1;
  parms[0]=rin;
  parms[1]=rhodisk;

  funreturn=user1_normalize_densities(eqline, parms, prim, &rhomax, &umax);
  if(funreturn!=0) return(funreturn);
 

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
int get_maxes(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min)
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
  parms[1]=100.0;
  parms[2]=0.3; // THETAEQ

  funreturn=user1_get_maxes(eqslice, parms,prim, bsq_max, pg_max, beta_min);
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
  int funreturn;
  int atmospheretype;

  if(WHICHPROBLEM==NORMALTORUS ||WHICHPROBLEM==THICKDISK || WHICHPROBLEM==KEPDISK){
    atmospheretype=1;
  }
  else if(WHICHPROBLEM==GRBJET){
    atmospheretype=2;
  }
  else {
    atmospheretype=1; // default
  }
 
  funreturn=user1_set_atmosphere(atmospheretype, whichcond, whichvel, ptrgeom, pr);
  if(funreturn!=0) return(funreturn);
 
  return(0);

}



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor, FTYPE *prceiling)
{
  int funreturn;
  
  funreturn=set_density_floors_default(ptrgeom, pr, prfloor, prceiling);
  if(funreturn!=0) return(funreturn);

  return(0);
}

int set_density_floors_alt(struct of_geom *ptrgeom, struct of_state *q, FTYPE *pr, FTYPE *U, FTYPE bsq, FTYPE *prfloor, FTYPE *prceiling)
{
  int funreturn;
  
  funreturn=set_density_floors_default_alt(ptrgeom, q, pr, U, bsq, prfloor, prceiling);
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



void adjust_fluxctstag_vpot(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  // not used
}

void adjust_fluxcttoth_vpot(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  // not used
}


void adjust_fluxcttoth_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{
  // not used
}

void adjust_fluxctstag_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  // not used
}



// not used
FTYPE calc_kappa_user(FTYPE rho, FTYPE B, FTYPE Tg, FTYPE Tr,FTYPE x,FTYPE y,FTYPE z)
{
  return(0.0);
}

FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z)
{  
  return(0.0);
}


// User's cooling function:

#define USERTHETACOOL       (h_over_r)  /* should be same as h_over_r */
#define USERTAUCOOL         (2.0*M_PI)          /* cooling time in number of rotational times : really USERTAUCOOL=2*M_PI would be 1 rotational time */
#define USERNOCOOLTHETAFACT     (1.0)           /* this times h_over_r and no more cooling there*/


int coolfunc_user(FTYPE h_over_r, FTYPE *pr, struct of_geom *geom, struct of_state *q,FTYPE (*dUcomp)[NPR])
{
  FTYPE X[NDIM],V[NDIM],r,th,R,Wcirc,cs_circ,rho,u,P,w,wcirc,dUcool;
  FTYPE taper0;
  int ii,jj, kk, pp;
  FTYPE pressure;
  FTYPE enk, enk0;
        
  FTYPE rpho;      
  FTYPE photoncapture;
  FTYPE rincool;
  FTYPE nocoolthetafactor,thetacool,taucool;



  // setup for macros
  nocoolthetafactor=USERNOCOOLTHETAFACT;
  thetacool=USERTHETACOOL;
  taucool=USERTAUCOOL;

  ii=geom->i;
  jj=geom->j;
  kk=geom->k;
  pp=geom->p;

  /* cooling function for maintaining fixed H/R */
  rho = pr[RHO] ;
  u = pr[UU] ;
  P=pressure_rho0_u_simple(ii,jj,kk,pp,rho,u);
  w = rho + u + P ;

  bl_coord_ijk_2(ii,jj,kk,CENT,X, V) ;
  r=V[1];
  th=V[2];
 
  rpho=2.0*(1.0+cos(2.0/3.0*(acos(-a))));

  //    trifprintf("rphoton=%lf\n", rpho);
  if(1 || r>rpho){ //SASMARK: cool always, including inside photon orbit
    photoncapture=1.0 ;
    //  trifprintf("r=%lf, photoncapture=%lf, rph=%lf \ n", r, photoncapture, rpho); 
  }
  else{
    photoncapture=0.0 ;
  }


  R = r*sin(th) ;
  enk=u*(gam-1.)/(pow(rho, gam));
  //enk0 = 1.e-3; //same as kappa in init.fishmon.c -- somehow wrong, is it because of wrong gam?!
  enk0 = 0.0043; //as read from ic's for thick torus
  //enk0=0.00016; //User's version
  //    enk0=0.00161;
  //    rin = (1. + h_over_r)*Risco;
  rincool=10.;
  /* crude approximation */
  Wcirc = 1./(a + pow(R,1.5)) ;
  cs_circ = thetacool/sqrt(R) ;
  //        wcirc = rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;

  wcirc =   rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;



  //    trifprintf("photoncapture=%lf, r=%lf, rpho=%lf \n", photoncapture, r, rpho);

  //  if(t > 0.){
  // dUcool = -(Wcirc/taucool)*( (w - wcirc)*(q->ucon[TT])) ;
  //     if(t > 0.){


  if(t > 0. && dt < taucool/Wcirc  && log(enk/enk0) > 0.) {

    //            dUcool = -(Wcirc/taucool)*( (w - wcirc)*(q->ucon[TT])*(q->ucov[TT])) ;


         

    // dUcool=-u*(Wcirc/taucool)*log(enk/enk0)*(q->ucon[TT])*photoncapture;



    dUcool=-u*(Wcirc/taucool)*log(enk/enk0)*photoncapture;

    //    dUcool*=COOLTAPER1(th);
    //  dUcool*=taper_func(R,Rhor);

    // shape function to avoid problems near pole
    //taper0=COOLTAPER(0);
    //dUcool*=1./(1.-1./taper0)+1./(1.-taper0)*COOLTAPER(th);
    //dUcool*=COOLTAPER2(th);
    // dUcool*=COOLTAPER3(th);
    // dUcool*=taper_func(R,Rhor); // don't cool inside horizon
  }
  else{
    dUcool = 0. ;
    // dUcool = (-u*log(enk/enk0)/dt)*(q->ucon[TT]) ;

  }

  //    dUcomp[RADSOURCE][UU]=dUcool;


  dUcomp[RADSOURCE][UU]=dUcool*(q->ucov[TT]);
  dUcomp[RADSOURCE][U1]=dUcool*(q->ucov[RR]);
  dUcomp[RADSOURCE][U2]=dUcool*(q->ucov[TH]);
  dUcomp[RADSOURCE][U3]=dUcool*(q->ucov[PH]);

  //                    trifprintf("ducomps are %g %g %g %g \n", dUcomp[RADSOURCE][UU], dUcomp[RADSOURCE][U1], dUcomp[RADSOURCE][U2],   dUcomp[RADSOURCE][U3]); 
  return(0) ;
}


int init_postvpot(int i, int j, int k, FTYPE *pr, FTYPE *pstag, FTYPE *ucons){

  return(0);
}



void set_coord_parms_nodeps_user(int defcoordlocal)
{

}


void set_coord_parms_deps_user(int defcoordlocal)
{
}

void write_coord_parms_user(int defcoordlocal, FILE *out)
{
}
void read_coord_parms_user(int defcoordlocal, FILE *in)
{
}

void read_coord_parms_mpi_user(int defcoordlocal)
{
}


void blcoord_user(FTYPE *X, FTYPE *V)
{
}


void dxdxp_analytic_user(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{

}
void set_points_user(void)
{

}
FTYPE setRin_user(int ihor, FTYPE ihoradjust)
{
  return(0.0);
}
 


int setihor_user(void)
{
  return(0);
}
