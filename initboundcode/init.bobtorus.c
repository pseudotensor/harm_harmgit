
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define MAXPASSPARMS 10

#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define THINDISKFROMMATHEMATICA 3
#define THINTORUS 4
#define THICKDISKFROMMATHEMATICA 5


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

// AKMARK: which problem
//#define WHICHPROBLEM THINDISKFROMMATHEMATICA // choice
//#define WHICHPROBLEM THICKDISKFROMMATHEMATICA // choice
#define WHICHPROBLEM THINTORUS
//#define WHICHPROBLEM NORMALTORUS

#define TORUSHASBREAKS 0   // AKMARK: 0 for usual torus, 1 for 3-region torus (constant l in regions 1 and 3)

#define DO_REMAP_MPI_TASKS (0)  //remap cores for performance (currently only on 8-core-per-node machines)

static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;

static FTYPE toruskappa;   // AKMARK: entropy constant KK from mathematica file
static FTYPE torusn;   // AKMARK: n from mathematica file (power of lambda in DHK03)
static FTYPE torusrmax;   // AKMARK: torus pressure max

static int read_data(FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR]);

#define SLOWFAC 1.0  /* reduce u_phi by this amount */

#include "init.torus.h"

#undef SLOWFAC
             
int prepre_init_specific_init(void)
{
  int funreturn;
  
  /////////////////////
  //PHI GRID SETUP
  /////////////////////

  dofull2pi = 0;   // AKMARK: do full phi
  
  global_fracphi = 0.5;   //phi-extent measured in units of 2*PI, i.e. 0.25 means PI/2; only used if dofull2pi == 0
  
  //binaryoutput=MIXEDOUTPUT;  //uncomment to have dumps, rdumps, etc. output in binary form with text header
   
  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);

}


// AKMARK: h/r
int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
#if( WHICHPROBLEM == THINDISKFROMMATHEMATICA )
  h_over_r=0.1;
#elif( WHICHPROBLEM == THINTORUS )
  h_over_r=0.1;
#elif( WHICHPROBLEM == THICKDISKFROMMATHEMATICA )
  h_over_r=0.3;
#else
  h_over_r=0.3;
#endif
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

  TIMEORDER = 2;
  DTr = 2000;

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



// AKMARK: grid coordinates
int init_defcoord(void)
{
  
  // define coordinate type
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  defcoord = SJETCOORDS;
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  defcoord = REBECCAGRID ;
#elif(WHICHPROBLEM==THINTORUS)
  defcoord = SJETCOORDS;
  //defcoord = REBECCAGRID ;
#elif(WHICHPROBLEM==GRBJET)
  // define coordinate type
  defcoord = JET4COORDS;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first


  // AKMARK: spin
#if(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  a = 0.;
#elif(WHICHPROBLEM==THINTORUS)
  a = 0.9;
#elif(WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  a = 0.;
#else
  a = 0.95;   //so that Risco ~ 2
#endif

 
  Rhor=rhor_calc(0);

  // AKMARK: hslope
  hslope = 0.13;  //sas: use a constant slope as Jon suggests in the comments
  //hslope = 1.04*pow(h_over_r,2.0/3.0);


  // AKMARK: inner (Rin) and outer (Rout) radii of simulation domain, R0
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.8 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3 * Rin;
  Rout = 1.e4;
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.92 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 50.;
#elif(WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.9 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 200.;
#elif(WHICHPROBLEM==THINTORUS)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.88 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 5.e4;
#elif(WHICHPROBLEM==GRBJET)
  setRin_withchecks(&Rin);
  R0 = -3.0;
  Rout = 1E5;
#endif

  /////////////////////
  // RADIAL GRID SETUP
  /////////////////////
  global_npow=1.0;  //don't change it, essentially equivalent to changing cpow2

  //radial hyperexponential grid
  global_npow2=4.0; //power exponent
  global_cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
  global_rbr = 100.;  //radius at which hyperexponentiation kicks in
  
  /////////////////////
  //ANGULAR GRID SETUP
  /////////////////////

  //transverse resolution fraction devoted to different components
  //(sum should be <1)
  global_fracdisk = 0.2;
  global_fracjet = 0.5;

  global_jetnu = 0.75;  //the nu-parameter that determines jet shape

  //subtractor, controls the size of the last few cells close to axis:
  //if rsjet = 0, then no modification <- *** default for use with grid cylindrification
  //if rsjet ~ 0.5, the grid is nearly vertical rather than monopolar,
  //                which makes the timestep larger
  global_rsjet = 0.0; 

  //distance at which theta-resolution is *exactly* uniform in the jet grid -- want to have this at BH horizon;
  //otherwise, near-uniform near jet axis but less resolution (much) further from it
  //the larger r0grid, the larger the thickness of the jet 
  //to resolve
  global_r0grid = Rin;    

  //distance at which jet part of the grid becomes monopolar
  //should be the same as r0disk to avoid cell crowding at the interface of jet and disk grids
  global_r0jet = 3;
    
  //distance after which the jet grid collimates according to the usual jet formula
  //the larger this distance, the wider is the jet region of the grid
  global_rjetend = 15;
    
  //distance at which disk part of the grid becomes monopolar
  //the larger r0disk, the larger the thickness of the disk 
  //to resolve
  global_r0disk = global_r0jet;

  //distance after which the disk grid collimates to merge with the jet grid
  //should be roughly outer edge of the disk
  global_rdiskend = 80;
  
  global_x10 = 3.3;  //radial distance in MCOORD until which the innermost angular cell is cylinrdical
  global_x20 = -1. + 1./totalsize[2];     //This restricts grid cylindrification to the one 
  //single grid closest to the pole (other cells virtually unaffeced, so there evolution is accurate).  
  //This trick minimizes the resulting pole deresolution and relaxes the time step.
  //The innermost grid cell is evolved inaccurately whether you resolve it or not, and it will be fixed
  //by POLEDEATH (see bounds.tools.c).

  return(0);
}



int init_global(void)
{
  int pl,pliter;
  int funreturn;


  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);

  cour = 0.8; //increase courant factor

  //////////////////
  // overrides for more detailed problem dependence


  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  lim[1] = lim[2] = lim[3] = PARALINE; //sas: it's already set in init.tools.c but reset it here just to make sure
  //also need to ensure that in para_and_paraenohybrid.h JONPARASMOOTH is set to 0 (resolves disk best) or 1 (resolves jet best)

  // AKMARK: cooling
#if(  WHICHPROBLEM==THINDISKFROMMATHEMATICA )
  cooling = COOLREBECCATHINDISK; //do Rebecca-type cooling; make sure enk0 is set to the same value as p/rho^\Gamma in the initial conditions (as found in dump0000).
#elif( WHICHPROBLEM==THINTORUS )
  //cooling = COOLREBECCATHINDISK; //do Rebecca-type cooling; make sure enk0 is set to the same value as p/rho^\Gamma in the initial conditions (as found in dump0000).
  cooling = NOCOOLING; //no cooling
#else
  cooling = NOCOOLING; //no cooling
#endif

  // AKMARK: Toth vs stag
  //FLUXB = FLUXCTTOTH;
  FLUXB = FLUXCTSTAG;

  // AKMARK: floors
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  //SASMARK: decrease magnetization by 2x to make it easier (still is around ~45>>1)
  BSQORHOLIMIT=0.5*1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=0.5*1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=0.5*1E3;
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






  // AKMARK: dumping frequencies, final time
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  /* output choices */
  tf = 1e4;

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
  DTr = 2000;

#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS)
  /* output choices */
  tf = 20000.;

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 100.;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 100.0;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 10.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 10.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 100.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 20000;

#elif(WHICHPROBLEM==GRBJET)
  /* output choices */
  tf = 5E5;
  
  DTd = 250.;                 /* dumping frequency, in units of M */
  DTavg = 250.0;
  DTener = 2.0;                       /* logfile frequency, in units of M */
  DTi = 10.0;                 /* image file frequ., in units of M */
  DTdebug = 250.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;                  /* restart file period in steps */
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



  beta = 1.e2 ;   // AKMARK: plasma beta (pgas/pmag)
  randfact = 4.e-2; //sas: as Jon used for 3D runs but use it for 2D as well

  toruskappa = 0.01;   // AKMARK: entropy constant KK from mathematica file
  torusn = 2. - 1.97;   // AKMARK: n from mathematica file (power of lambda in DHK03)
  torusrmax = 20.;   // AKMARK: torus pressure max

  // AKMARK: torus inner radius
#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
  toruskappa = 1e-3; 
  torusrmax = 12.; 
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  rin = 20. ;
#elif(WHICHPROBLEM==THINTORUS)
  rin = 10. ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET)
  rin = Risco;
#endif

  


#if( ANALYTICMEMORY == 1 && WHICHPROBLEM != THINDISKFROMMATHEMATICA && WHICHPROBLEM != THICKDISKFROMMATHEMATICA )
  //SASMARK restart: need to populate panalytic with IC's; DO NOT do this 
  //when reading the ICs in from a file since then need to carry the file around
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
    MYFUN(init_primitives(panalytic,pstaganalytic,GLOBALPOINT(utemparray),vpotanalytic,Bhatanalytic,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }
#endif
  // check rmin
  check_rmin();


  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA ) 
  //read initial conditions from input file for the Mathematica-generated thin-disk ICs
#if(ANALYTICMEMORY==0)
#error Cannot do THINDISKFROMMATHEMATICA problem with ANALYTICMEMORY==0.  Please set ANALYTICMEMORY = 1.
#endif
  read_data(panalytic);
#endif

  funreturn=user1_init_primitives(prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

  return(0);


}









int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindiskfrommathematica(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thintorus(int *whichvel, int*whichcoord, int ti, int tj, int tk, FTYPE *pr, FTYPE *pstag);

  // AKMARK: check which function is called for your WHICHPROBLEM, and change parameters in it below
#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  return(init_dsandvels_thindiskfrommathematica(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINTORUS)
  return(init_dsandvels_thintorus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#endif

}



#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3

//#define FIELDTYPE BLANDFORDQUAD
#define FIELDTYPE DISKFIELD



// assumes normal field in pr
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av, q;
  FTYPE r,th;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);
#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
  FTYPE fieldhor;
#endif



  vpot=0.0;


  if(l==3){// A_\phi

    r=V[1];
    th=V[2];



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

      // AKMARK: magnetic loop radial wavelength
#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
#define STARTFIELD (1.1*rin)
      fieldhor=0.28;
#elif(WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
#define STARTFIELD (1.1*rin)
      fieldhor=0.28;
#endif
      // average of density that lives on CORN3
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
        u_av=AVGN_for3(prim,i,j,k,UU);
      }

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
      //SASMARK: since u was randomly perturbed, may need to sync the u across tiles to avoid monopoles
      if(r > STARTFIELD) q = ((u_av/umax) - 0.2)*pow(r,0.75) ;
      else q = 0. ;
      //trifprintf("rhoav=%g q=%g\n", rho_av, q);

      if(q > 0.){
        //       vpot += q*q*sin(log(r/STARTFIELD)/fieldhor)* (1. + 0.02 * (ranc(0,0) - 0.5))  ;
        vpot += q*q*sin(log(r/STARTFIELD)/fieldhor) ;
      }
#else
      q = rho_av / rhomax - 0.2;
      if (q > 0.)      vpot += q;
#endif
    
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

  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM==THINTORUS){
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



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  int funreturn;
  
  funreturn=set_density_floors_default(ptrgeom, pr, prfloor);
  if(funreturn!=0) return(funreturn);

  return(0);
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
  int group_by_node_set_myid(int n1tile, int n2tile, int n3tile);
  int retval;
 
  // default is to do nothing
  //  retval=jet_set_myid();
#if(DO_REMAP_MPI_TASKS)
  retval=group_by_node_set_myid( 2, 2, 2 );  //2x2x2 for Odyssey/QueenBee
  //retval=group_by_node_set_myid( 3, 2, 2 );  //3x2x2 for Kraken
#else
  retval=0;
#endif  
  // do other things?

  return(retval);

}

