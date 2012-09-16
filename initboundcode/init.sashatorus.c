
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
#define FLOORFACTOR (2.)
#define BSQOUPREFACT (5.)

#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define THINDISKFROMMATHEMATICA 3
#define THINTORUS 4
#define THICKDISKFROMMATHEMATICA 5
#define NSTAR 6

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
//#define WHICHPROBLEM NSTAR

#define TORUSHASBREAKS 0   // AKMARK: 0 for usual torus, 1 for 3-region torus (constant l in regions 1 and 3)

#define DO_REMAP_MPI_TASKS (0)  //remap cores for performance (currently only on 8-core-per-node machines)

static const FTYPE aphipow = 2.5;
static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;

static FTYPE toruskappa;   // AKMARK: entropy constant KK from mathematica file
static FTYPE torusn;   // AKMARK: n from mathematica file (power of lambda in DHK03)
extern FTYPE torusrmax;   // AKMARK: torus pressure max
FTYPE t_transition;
FTYPE global_vpar0;
FTYPE global_dipole_alpha;
FTYPE global_OmegaNS;

static int read_data(FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR]);
FTYPE is_inside_torus_freeze_region( FTYPE r, FTYPE th );
void add_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
void add_3d_fieldonly_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
void null_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
void null_3d_fieldonly_fullloop(FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

#if( DOFREEZETORUS )
void add_torus_magnetic_fields(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );
void save_torus_allvars(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );
FTYPE BASEMACP0A1(global_ptorus,N1M,N2M,N3M,NPR);       /* space for primitive vars */
FTYPE BASEMACP0A1(global_pstagtorus,N1M,N2M,N3M,NPR);       /* space for primitive vars */
FTYPE BASEMACP0A1(global_uconstorus,N1M,N2M,N3M,NPR);       /* space for conserved vars */
FTYPE PTRDEFGLOBALMACP0A1(global_ptorus,N1M,N2M,N3M,NPR);
FTYPE PTRDEFGLOBALMACP0A1(global_pstagtorus,N1M,N2M,N3M,NPR);
FTYPE PTRDEFGLOBALMACP0A1(global_uconstorus,N1M,N2M,N3M,NPR);
#endif

#define SLOWFAC 1.0		/* reduce u_phi by this amount */

#include "init.torus.h"

#undef SLOWFAC

FTYPE normglobal;
int inittypeglobal; // for bounds to communicate detail of what doing

int prepre_init_specific_init(void)
{
  int funreturn;
#if( DOFREEZETORUS )
  FTYPE valueinit = 0.0;
  int i, j, k, pliter, pl;
#endif
  
#if( DOFREEZETORUS )
  GLOBALPOINT(global_ptorus) = (FTYPE PTRMACP0A1(global_ptorus,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(global_ptorus,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(global_pstagtorus) = (FTYPE PTRMACP0A1(global_pstagtorus,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(global_pstagtorus,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(global_uconstorus) = (FTYPE PTRMACP0A1(global_uconstorus,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(global_uconstorus,N1BND,N2BND,N3BND,0)));
  
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(global_ptorus,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(global_pstagtorus,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(global_uconstorus,i,j,k,pl) = valueinit;
  }
#endif
  
  /////////////////////
  //PHI GRID SETUP
  /////////////////////

  dofull2pi = 1;   // AKMARK: do full phi
  
  global_fracphi = 1.0;   //phi-extent measured in units of 2*PI, i.e. 0.25 means PI/2; only used if dofull2pi == 0
  
  binaryoutput=MIXEDOUTPUT;  //uncomment to have dumps, rdumps, etc. output in binary form with text header
   
  t_transition = 1.;
  global_vpar0 = 0.0;
  
#if(WHICHPROBLEM==NSTAR)
  global_dipole_alpha = 60. * M_PI / 180.;
  global_OmegaNS = 0.2;
#endif
  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);

}

#if( DOFREEZETORUS )
void add_torus_magnetic_fields(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{
  add_3d_fieldonly_fullloop(GLOBALPOINT(global_ptorus),prim);
  add_3d_fieldonly_fullloop(GLOBALPOINT(global_pstagtorus),pstag);
  add_3d_fieldonly_fullloop(GLOBALPOINT(global_uconstorus),ucons);
}

void save_torus_allvars(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{
  copy_3dnpr_fullloop(prim,GLOBALPOINT(global_ptorus));
  copy_3dnpr_fullloop(pstag,GLOBALPOINT(global_pstagtorus));
  copy_3dnpr_fullloop(ucons,GLOBALPOINT(global_uconstorus));
}
#endif

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
  //UTOPRIMVERSION = UTOPRIM2DFINAL;

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
  DTr = 6000;
  tf = 1e5;
  DOENERDIAG=0;
  //DOAVGDIAG=0;  //set here to override after restart
  //DODUMPDIAG=0; //=0 switches off all dumps (including floor dumps)
  
  t_transition = 1.;
  global_vpar0 = 0.0;
  
  BSQORHOLIMIT=FLOORFACTOR*1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=FLOORFACTOR*BSQOUPREFACT*1E2; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=FLOORFACTOR*BSQOUPREFACT*1E2;

  
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
#elif(WHICHPROBLEM==NSTAR)
  defcoord = SNSCOORDS;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first
#if(WHICHPROBLEM==NSTAR)
  MBH=0.2;
#endif  
  

// AKMARK: spin
#if(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  a = 0.;
#elif(WHICHPROBLEM==THINTORUS)
  a = 0.99;
#elif(WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  a = 0.;
#elif(WHICHPROBLEM==NSTAR)
  //flat metric so use this instead of Omega_F
  a = 0.0;  //Omega_F = a; phi-velocity: v^\phi = a
#else
  a = 0.95;   //so that Risco ~ 2
#endif

#if( DOAUTOCOMPUTEENK0 )
  //this should be set to the final entropy constant
  //this is done *automatically* for the case of WHICHPROBLEM == THINTORUS
  //and THINTORUS_NORMALIZE_DENSITY == 1
  global_toruskappafinal = 0.01;  
#endif
  
  toruskappa = 0.01;   // AKMARK: entropy constant KK from mathematica file
  torusn = 2. - 1.75;   // AKMARK: n from mathematica file (power of lambda in DHK03)
  torusrmax = 34.; //22.7; //37.1; //22.82; //34.1;   // AKMARK: torus pressure max
  
  beta = 1.e2 ;   // AKMARK: plasma beta (pgas/pmag)
  randfact = 4.e-2; //sas: as Jon used for 3D runs but use it for 2D as well
  
#if(MCOORD==KSCOORDS)
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);
  
  // AKMARK: torus inner radius
#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
  toruskappa = 1e-3; 
  torusrmax = 12.; 
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  rin = 20. ;
#elif(WHICHPROBLEM==THINTORUS)
  rin = 15. ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET)
  rin = Risco;
#endif

#endif

#if(WHICHPROBLEM==NSTAR)
  rin = 1.;
#endif

  Rhor = rin;
  Risco = rin;
 
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
  Rin = 0.83 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 1.e5;
#elif(WHICHPROBLEM==GRBJET)
	setRin_withchecks(&Rin);
	R0 = -3.0;
  Rout = 1E5;
#elif(WHICHPROBLEM==NSTAR)
  Rin = rin;
  Rout = 1e3;
  R0 = 0.5;
#endif

  /////////////////////
  // RADIAL GRID SETUP
  /////////////////////
  global_npow=1.0;  //don't change it, essentially equivalent to changing cpow2

  //radial hyperexponential grid
  global_npow2=4.0; //power exponent
  global_cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
  global_rbr = 1000.;  //radius at which hyperexponentiation kicks in
  if(WHICHPROBLEM==NSTAR){
    global_rbr = 100.;
  }
  
  /////////////////////
  //ANGULAR GRID SETUP (so far irrelevant for WHICHPROBLEM==NSTAR)
  /////////////////////

  //transverse resolution fraction devoted to different components
  //(sum should be <1)
  global_fracdisk = 0.36;
  global_fracjet = 0.3;

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
  global_r0grid = 5.0*Rin;    

  //distance at which jet part of the grid becomes monopolar
  //should be the same as r0disk to avoid cell crowding at the interface of jet and disk grids
  global_r0jet = Rin;
    
  //distance after which the jet grid collimates according to the usual jet formula
  //the larger this distance, the wider is the jet region of the grid
  global_rjetend = 5;
    
  //distance at which disk part of the grid becomes monopolar
  //the larger r0disk, the larger the thickness of the disk 
  //to resolve
  global_r0disk = Rin+0*global_r0jet;

  //distance after which the disk grid collimates to merge with the jet grid
  //should be roughly outer edge of the disk
  global_rdiskend = 300.;
  
  global_x10 = 3.0;  //radial distance in MCOORD until which the innermost angular cell is cylinrdical
  global_x20 = -1. + 1./totalsize[2];     //This restricts grid cylindrification to the one 
    //single grid closest to the pole (other cells virtually unaffeced, so there evolution is accurate).  
    //This trick minimizes the resulting pole deresolution and relaxes the time step.
    //The innermost grid cell is evolved inaccurately whether you resolve it or not, and it will be fixed
    //by POLEDEATH (see bounds.tools.c).
  
  if (0&&defcoord==SNSCOORDS) {
    global_x20 = 0;  //value of |x2| at which to start transition from cylindric to spherical coords
    hslope = 2.; //hslope*global_x20 is the value of |x2| at which this transition finishes
  }

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
  //lim[1] = lim[2] = lim[3] = PARALINE; //sas: it's already set in init.tools.c but reset it here just to make sure
  lim[1] = lim[2] = lim[3] = MC; //sas: it's already set in init.tools.c but reset it here just to make sure
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
#if(THINTORUS_NORMALIZE_DENSITY && WHICHPROBLEM == THINTORUS)
//scale up to match the usual values after the corresponding density normalization in init_thintorus()
  RHOMIN *= 70;
  UUMIN *= 70;
#endif
#elif(WHICHPROBLEM==NSTAR)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=NSSURFACE;
  //  rescaletype=1;
  rescaletype=4;
  //SASMARK: decrease magnetization by 2x to make it easier (still is around ~45>>1)
  BSQORHOLIMIT=FLOORFACTOR*1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=FLOORFACTOR*BSQOUPREFACT*1E2; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=FLOORFACTOR*BSQOUPREFACT*1E2;
  RHOMIN = 1E-4;
  UUMIN = 1E-4;
  GAMMADAMP=50.0;
  
  if(DOEVOLVERHO){
    // GODMARK -- unstable beyond about 25, but can sometimes get away with 1000
    GAMMAMAX=100.0; // when we think gamma is just too high and may cause unstable flow, but solution is probably accurate.
  }
  else{
    GAMMAMAX=2000.0;
  }
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
  tf = 1e5; //also check post_init_specific_init()

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 10.;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 10.0;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 100.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 5.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 100.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 6000;  //also see post_init_specific_init()

#elif(WHICHPROBLEM == NSTAR)
  /* output choices */
  tf = 10.; //also check post_init_specific_init()
  
  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 2*M_PI/global_OmegaNS;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 2*M_PI/global_OmegaNS;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 2*M_PI/global_OmegaNS;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = (1./32.)*2*M_PI/global_OmegaNS;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];
  
  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 2*M_PI/global_OmegaNS;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 6000;  //also see post_init_specific_init()
  
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

int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);

//  BSQORHOLIMIT=FLOORFACTOR*1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
//  BSQOULIMIT=FLOORFACTOR*BSQOUPREFACT*1E2; // was 1E3 but latest BC test had 1E4
//  UORHOLIMIT=FLOORFACTOR*BSQOUPREFACT*1E2;

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
#if(MCOORD==KSCOORDS)
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);
#else
  Rhor = rin;
  Risco = rin;
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



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;
  int inittype;

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA ) 
  //read initial conditions from input file for the Mathematica-generated thin-disk ICs
#if(ANALYTICMEMORY==0)
#error Cannot do THINDISKFROMMATHEMATICA problem with ANALYTICMEMORY==0.  Please set ANALYTICMEMORY = 1.
#endif
  read_data(panalytic);
#endif

  inittype=1;

  funreturn=user1_init_primitives(inittype, prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

  return(0);


}



int init_dsandvels(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindiskfrommathematica(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thintorus(int *whichvel, int*whichcoord, int ti, int tj, int tk, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_nstar(int *whichvel, int*whichcoord, int ti, int tj, int tk, FTYPE *pr, FTYPE *pstag);

// AKMARK: check which function is called for your WHICHPROBLEM, and change parameters in it below
#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA)
  return(init_dsandvels_thindiskfrommathematica(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINTORUS)
  return(init_dsandvels_thintorus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==NSTAR)
  return(init_dsandvels_nstar(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#endif

}



#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3
#define DISKBHFIELD 4
#define NSFIELD 5

//Options for DISKBHFIELD
#define BHFIELDVAL (1.0) //(3.4142135623730950488)
#define BHFIELDNU (-1.0) //negative means hyperbolic field lines
#define BHFIELDALPHA (1.5)
#define BHFIELDR (1.)    //field stops at r = BHFIELDR * rin
#define BHFIELDDZDR (.9) //asymptotic slope of field boundary (=dz/dr)
#define BHFIELDLOGPOW (2.) //power to which the log prefactor is taken

#define NSFIELDVAL (1.5*3.162277660168379332*2*3*3*0.5)

//#define FIELDTYPE DISKBHFIELD
//#define FIELDTYPE DISKFIELD
#define FIELDTYPE NSFIELD

FTYPE vpotbh_normalized( FTYPE r, FTYPE th )
{
  FTYPE rh;
  FTYPE vpotbh;
  rh = rhor_calc(0);
  if(BHFIELDNU>=0) {
    //normalized vector potential: total vpot through BH equals some constant order unity
    vpotbh = pow(r/rh,BHFIELDNU)*(1 - fabs(cos(th)));
    if( vpotbh > (1 - fabs(cos(M_PI/4.))) ) vpotbh = (1 - fabs(cos(M_PI/4.)));
    //rescale the flux to amplitude given by BHFLUX and add it up to vector potential
  }
  else if(BHFIELDNU<0) {
    //roughly uniform Bz at constant slices of z = r*cos(th) nearly all the way to the edges of the torus
    vpotbh = pow( r*sin(th)/(BHFIELDR*rin),2 ) 
      / pow( 1+pow(
		    fabs(r*cos(th))/(BHFIELDR*rin*BHFIELDDZDR*(1+pow(log10(1+r/rin),BHFIELDLOGPOW))),
	       BHFIELDALPHA),
	2./BHFIELDALPHA );
    if( vpotbh > 1 ) vpotbh = 1;
  }
  return(vpotbh);
}

FTYPE vpotns_normalized( int i, int j, int k, int loc, FTYPE *V, int l )
{
  FTYPE get_ns_alpha();
  FTYPE alpha = get_ns_alpha();
  FTYPE vpot;
  FTYPE r = V[1], th = V[2], ph = V[3], phi = V[3], theta = V[2];  //the latter two to avoid typos
#if(0)
  //normalized vector potential: total vpot through NS equals some constant order unity
  //vpot = 1 - fabs(cos(th));  //split-monopole
  //vpot = 1 - cos(th);        //monopole
  vpot = sin(th)*sin(th)/r;        //dipole
  return(vpot);
#elif(1)
  //tilted dipole
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE Ad1, Ad2, Ad3;
  FTYPE Adr, Adtheta, Adphi;

  dxdxprim_ijk( i, j, k, loc, dxdxp );

//incorrect:
//  Adr = pow(r,-2)*(-(r*cos(ph)*pow(cos(th),2)*sin(alpha)) + 
//		  pow(sin(th),2)*sin(ph)*(-(cos(alpha)*cos(ph)) + r*sin(alpha)*sin(ph)) + 
//		  r*cos(alpha)*cos(th)*pow(cos(ph),2)*sin(th));
//  
//  Adtheta = -(pow(r,-2)*sin(ph)*(r*pow(cos(th),2)*sin(alpha) + 
//			     sin(th)*(-(r*cos(ph)*cos(alpha + th)) + cos(alpha)*sin(ph)*sin(th))));
//  
//  Adphi = -(pow(r,-2)*sin(th)*(cos(th)*(-(r*sin(alpha)) + cos(alpha)*sin(ph)) + r*cos(alpha)*cos(ph)*sin(th)));

#if(0)
  //mu-Omega plane at phi = pi/2 (along y-axis) at t = 0
  Adr = 0;
  Adtheta = cos(ph)*pow(r,-1)*sin(alpha);
  Adphi = pow(r,-1)*sin(th)*(-(cos(th)*sin(alpha)*sin(ph)) + 
			     cos(alpha)*sin(th));
#elif(1)  
  //mu-Omega plane at phi = 0 (along x-axis) at t = 0
  Adr = 0;
  Adtheta = -sin(alpha)*sin(ph) / r;
  Adphi = sin(th) * ( -cos(th)*cos(ph)*sin(alpha) + cos(alpha)*sin(th) ) / r; 
#endif
  
  if( 1 == l ){
    //Ad1 = dxdxp[1][1] * Adr + dxdxp[2][1] * Adtheta;
    Ad1 = Adr;
    return(NSFIELDVAL*Ad1);
  }
  else if( 2 == l ){
    //Ad2 = dxdxp[1][2] * Adr + dxdxp[2][2] * Adtheta;
    Ad2 = Adtheta;
    return(NSFIELDVAL*Ad2);
  }
  else if( 3 == l ){
    //Ad3 = dxdxp[3][3] * Adphi;
    Ad3 = Adphi;
    return(NSFIELDVAL*Ad3);
  }
#endif
}

FTYPE get_ns_alpha()
{
  return(global_dipole_alpha);
}
//returns enclosed flux between ph1, ph2 and th1, th2
//flux = \int_phi1^phi2 A_phi dphi|_th1^th2 - \int_th1^th2 A_th dth|_ph1^ph2
FTYPE vpotns_flux( FTYPE r, FTYPE th1, FTYPE th2, FTYPE ph1, FTYPE ph2)
{
  FTYPE get_ns_alpha();
  FTYPE alpha = get_ns_alpha();
  FTYPE sinth1 = sin(th1);
  FTYPE sinth2 = sin(th2);
  FTYPE sinth1sq = sinth1*sinth1;
  FTYPE sinth2sq = sinth2*sinth2;
  FTYPE int_Ath_dth = - ( (th2-th1)*(sin(ph2)-sin(ph1))*sin(alpha) );
  FTYPE int_Aph_dph =   ( (ph2-ph1)*(sinth2sq-sinth1sq)*cos(alpha) 
			  -(sin(2*th2)-sin(2*th1))*(sin(ph2)-sin(ph1))*sin(alpha)*0.5 );
  return( NSFIELDVAL*(-int_Ath_dth + int_Aph_dph)/r );
}

FTYPE dfluxns( FTYPE r, FTYPE Omega, FTYPE phi, FTYPE th1, FTYPE th2, FTYPE t, FTYPE dt )
{
  FTYPE vpotns_flux( FTYPE r, FTYPE th1, FTYPE th2, FTYPE ph1, FTYPE ph2);
  FTYPE phi2 = phi - Omega*t;
  FTYPE phi1 = phi2 - Omega*dt;
  return( vpotns_flux(r, th1, th2, phi1, phi2) );
}

FTYPE is_inside_torus_freeze_region( FTYPE r, FTYPE th )
{
  FTYPE vpotbh_normalized( FTYPE r, FTYPE th );
  int is_inside;
  FTYPE vpotbh;
  
  vpotbh = vpotbh_normalized(r,th);
  is_inside = (vpotbh>=1);
  
  return(is_inside);
}

// assumes normal field in pr
// SUPERNOTE: A_i must be computed consistently across all CPUs.  So, for example, cannot use randomization of vector potential here.
int init_vpot_user(int *whichcoord, int l, SFTYPE time, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  FTYPE vpotbh_normalized( FTYPE r, FTYPE th );
  FTYPE vpotns_normalized( int i, int j, int k, int loc, FTYPE *V, int l );
  SFTYPE rho_av, u_av, q;
  FTYPE r,th,ph;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);
#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM==THICKDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
  FTYPE fieldhor;
#endif
  FTYPE rh;
  FTYPE vpotbh, vpotns;
  
#if(MCOORD==KSCOORDS)
  rh = rhor_calc(0);
#else
  rh = rin;
#endif


  vpot=0.0;
  r=V[1];
  th=V[2];
  ph=V[3];
  
  
  
  if(l==-3){// A_\phi for bh field
    r=V[1];
    th=V[2];
    if( FIELDTYPE==DISKBHFIELD ) {
      vpotbh = vpotbh_normalized(r, th);
      vpot += BHFIELDVAL * vpotbh;
    }
  }
  
  //NS field
  if( FIELDTYPE==NSFIELD && l >= 1 && l <= 3 ) {
    vpotns = vpotns_normalized(i, j, k, loc, V, l);
    vpot += vpotns;
  }

  if(l==3){// A_\phi for disk


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
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKBHFIELD)||(FIELDTYPE==DISKVERT)){

// AKMARK: magnetic loop radial wavelength
#if( WHICHPROBLEM == THINTORUS ) 
#define STARTFIELD (1.*rin)
      fieldhor=0.194;
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
      if(r > STARTFIELD) q = (pow(r,aphipow)*(rho_av/rhomax));
      else q = 0. ;
      //trifprintf("rhoav=%g q=%g\n", rho_av, q);

      if(q > 0.){
	//       vpot += q*q*sin(log(r/STARTFIELD)/fieldhor)* (1. + 0.02 * (ranc(0,0) - 0.5))  ;
	vpot += q*q; //*sin(log(r/STARTFIELD)/fieldhor) ;
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



int init_vpot2field_user(SFTYPE time, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int getmax_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax);
  int normalize_field_local_nodivb(FTYPE targbeta, FTYPE rhomax, FTYPE amax, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], 
				   FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], 
				   FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], 
				   FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
  FTYPE get_maxval(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir );
  int compute_vpot_from_gdetB1( FTYPE (*prim)[NSTORE2][NSTORE3][NPR], 
			       FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], 
			       FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], 
			       FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
  
  int normalize_field_diskonly(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);

  int add_vpot_bhfield_user_allgrid( FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  FTYPE get_maxprimvalrpow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE rpow, int pl );
  
  int funreturn;
  int fieldfrompotential[NDIM];
  FTYPE rhomax, umax, amax;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int i, j, k;
  FTYPE vpar;
  FTYPE bsq;

  //convert A to staggered pstag, centered prim and ucons, unsure about Bhat  
  funreturn=user1_init_vpot2field_user(time, fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
  
#if(DO_OPTIMIZE_DISK_FLUX) //if optimizing disk flux
  getmax_densities(prim, &rhomax, &umax);
  //amax = get_maxval( A, 3 );
  amax = get_maxprimvalrpow( prim, aphipow, RHO );
  trifprintf("amax = %g\n", amax);
  
  //by now have the fields (centered and stag) computed from vector potential
    
  //here need to 
  //1) compute bsq
  //2) rescale field components such that beta = p_g/p_mag is what I want (constant in the main disk body and tapered off to zero near torus edges)
  normalize_field_local_nodivb( beta, rhomax, amax, prim, pstag, ucons, A, Bhat );
  
  //3) re-compute vector potential by integrating up ucons (requires playing with MPI)
  compute_vpot_from_gdetB1( prim, pstag, ucons, A, Bhat );
	
  //4) call user1_init_vpot2field_user() again to recompute the fields
  //convert A to staggered pstag, centered prim and ucons, unsure about Bhat  
  funreturn=user1_init_vpot2field_user(fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
#endif

  //normalize disk field -- the usual normalize, just call it from here instead of the usual place in init.tools.c:user1_init_primitives()
  normalize_field_diskonly(prim, pstag, ucons, A, Bhat);

#if( DOFREEZETORUS )
  //save the torus field: it will be added on later
  save_torus_allvars(prim, pstag, ucons, A );
  //null out the field everywhere
  null_3d_fieldonly_fullloop(prim);
  null_3d_fieldonly_fullloop(pstag);
  null_3d_fieldonly_fullloop(ucons);
#endif
  
#if(FIELDTYPE==DISKBHFIELD) //if doing bh field
  //compute vpot again after field was normalized
  compute_vpot_from_gdetB1( prim, pstag, ucons, A, Bhat );

  //add vpot describing BH field 
  add_vpot_bhfield_user_allgrid( A, prim );

  funreturn=user1_init_vpot2field_user(fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
#endif
  
  
#if( WHICHPROBLEM == NSTAR && 0)  //USE THIS WITH GRAVITY TO RESET 4-vel to zero
  //ensure that VEL4 field velocity is zero (so that the only motion is
  //parallel to field lines)
  //compute parallel velocity (to the poloidal field) of a ZAMO
  compute_vpar(pr, ptrgeom, &vpar);
  //set_vpar(global_vpar0, GAMMA_MAX, ptrgeom, pr);
  
  //set field velocity to zero
  //for this, first, reset full 4-velocity (VEL4) to zero
  DLOOPA(pl) ucon[pl] = 0.0;
  ucon2pr(WHICHVEL, ucon, ptrgeom, pr);  //this does not use t-component of ucon, so no need to set it
  
  //then reinstate the ZAMO velocity along field lines
  set_vpar(vpar, GAMMAMAX, ptrgeom, pr);
#endif

#if( WHICHPROBLEM == NSTAR )  //set parallel velocity to global_vpar0 AND set rho and u at 1/FRACBSQORHO and 1/FRACBSQOU times the floor
  FULLLOOP {
    get_geometry(i, j, k, CENT, ptrgeom);
    //then reinstate the ZAMO velocity along field lines
    set_vpar(global_vpar0, GAMMAMAX, ptrgeom, MAC(prim,i,j,k));
    if(bsq_calc(MAC(prim,i,j,k),ptrgeom,&bsq)>=1){
      dualfprintf(fail_file,"bsq_calc:bsq_calc: failure\n");
      return(1);
    }
    MAC(prim,i,j,k)[RHO] = bsq/(FRACBSQORHO*BSQORHOLIMIT);
    MAC(prim,i,j,k)[UU]  = bsq/(FRACBSQOU*BSQOULIMIT);
  }    
#endif
  
  return(0);


}

int add_vpot_bhfield_user_allgrid( FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i, j, k, loc, userdir, dir, whichcoord;
  FTYPE vpotuser[NDIM];
  FTYPE X[NDIM], V[NDIM];
  
  //add in bh field vpot -- call:
  ZSLOOP(0,N1-1+SHIFT1,0,N2-1+SHIFT2,0,N3-1+SHIFT3) {
    // get user vpot in user coordinates (assume same coordinates for all A_{userdir})
    DLOOPA(userdir){
      loc = CORN1 - 1 + userdir; //CORRECT?
      bl_coord_ijk_2(i, j, k, loc, X, V); 
      //note minus sign: -userdir; this indicates to only init bh field
      init_vpot_user(&whichcoord, -userdir, t, i,j,k, loc, prim, V, &vpotuser[userdir]);
    }
    // convert from user coordinate to PRIMECOORDS
    ucov_whichcoord2primecoords(whichcoord, i, j, k, loc, vpotuser);
    
    DIMENLOOP(dir){
      NOAVGCORN_1(A[dir],i,j,k) += vpotuser[dir];
    }
  }
  return(0);
}

//compute vector potential assuming B_\phi = 0 and zero flux at poles
//(not tested in non-axisymmetric field distribution but in principle should work)
int compute_vpot_from_gdetB1( FTYPE (*prim)[NSTORE2][NSTORE3][NPR], 
				 FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], 
				 FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], 
				 FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int i, j, k;
  int jj;
  int cj;
  int dj, js, je, jsb, jeb;
  int dointegration;
  struct of_gdetgeom geomfdontuse[NDIM];
  struct of_gdetgeom *ptrgeomf[NDIM];
  FTYPE gdetnosing;
  FTYPE igdetgnosing[NDIM];
  int dir;
  int finalstep;
  
  if(steppart==TIMEORDER-1) finalstep=1; else finalstep=0;
  
  if(FLUXB!=FLUXCTSTAG) {
    dualfprintf( fail_file, "compute_vpot_from_gdetB1 only works for FLUXB = FLUXCTSTAG\n" );
    myexit(121);
  }

  DLOOPA(jj) ptrgeomf[jj]=&(geomfdontuse[jj]);

  //first, bound to ensure consistency of magnetic fields across tiles
  bound_allprim(STAGEM1,finalstep,t,prim,pstag,ucons, USEMPI);

  if( ncpux2 == 1 ) {
    //1-cpu version
    for (i=0; i<N1+SHIFT1; i++) {
      for (k=0; k<N3+SHIFT3; k++) {
	//zero out starting element of vpot
	NOAVGCORN_1(A[3],i,0,k) = 0.0;
	//integrate vpot along the theta line
	for (j=0; j<N2/2; j++) {
	  dir=1;
	  get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
	  set_igdetsimple(ptrgeomf[dir]);
	  igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
	  gdetnosing = 1.0/igdetgnosing[dir];
	  //take a loop along j-line at a fixed i,k and integrate up vpot
	  NOAVGCORN_1(A[3],i,jp1mac(j),k) = NOAVGCORN_1(A[3],i,j,k) + MACP0A1(pstag,i,j,k,B1)*gdetnosing*dx[2];
	  NOAVGCORN_1(A[1],i,j,k) = 0;
	  NOAVGCORN_1(A[2],i,j,k) = 0;
	}
	NOAVGCORN_1(A[3],i,N2,k) = 0.0;
	//integrate vpot along the theta line
	for (j=N2; j>N2/2; j--) {
	  dir=1;
	  get_geometry_gdetonly(i, jm1mac(j), k, FACE1-1+dir, ptrgeomf[dir]);
	  set_igdetsimple(ptrgeomf[dir]);
	  igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
	  gdetnosing = 1.0/igdetgnosing[dir];
	  //take a loop along j-line at a fixed i,k and integrate up vpot
	  NOAVGCORN_1(A[3],i,jm1mac(j),k) = NOAVGCORN_1(A[3],i,j,k) - MACP0A1(pstag,i,jm1mac(j),k,B1)*gdetnosing*dx[2];
	  NOAVGCORN_1(A[1],i,j,k) = 0;
	  NOAVGCORN_1(A[2],i,j,k) = 0;
	}
      }
    }
  }
  else if( 0&&ncpux2 == 2 ) {
    //dualfprintf( fail_file, "Got into 2-cpu version, mycpupox[2] = %d\n", mycpupos[2] );
    //2-cpu version
    for (i=0; i<N1+SHIFT1; i++) {
      for (k=0; k<N3+SHIFT3; k++) {
	if( mycpupos[2] == 0 ) {
	  //zero out starting element of vpot
	  NOAVGCORN_1(A[3],i,0,k) = 0.0;
	  //integrate vpot along the theta line
	  for (j=0; j<N2; j++) {
	    dir=1;
	    get_geometry_gdetonly(i, j, k, FACE1-1+dir, ptrgeomf[dir]);
	    set_igdetsimple(ptrgeomf[dir]);
	    igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
	    gdetnosing = 1.0/igdetgnosing[dir];
	    //take a loop along j-line at a fixed i,k and integrate up vpot
	    NOAVGCORN_1(A[3],i,jp1mac(j),k) = NOAVGCORN_1(A[3],i,j,k) + MACP0A1(pstag,i,j,k,B1)*gdetnosing*dx[2];
	    NOAVGCORN_1(A[1],i,j,k) = 0;
	    NOAVGCORN_1(A[2],i,j,k) = 0;
	  }
	  //copy A[3] -> B[3] before bounding
	  MACP0A1(pstag,i,N2-1,k,B3) = NOAVGCORN_1(A[3],i,N2,k);
	}
	else {
	  NOAVGCORN_1(A[3],i,N2,k) = 0.0;
	  //integrate vpot along the theta line
	  for (j=N2; j>0; j--) {
	    dir=1;
	    get_geometry_gdetonly(i, jm1mac(j), k, FACE1-1+dir, ptrgeomf[dir]);
	    set_igdetsimple(ptrgeomf[dir]);
	    igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
	    gdetnosing = 1.0/igdetgnosing[dir];
	    //take a loop along j-line at a fixed i,k and integrate up vpot
	    NOAVGCORN_1(A[3],i,jm1mac(j),k) = NOAVGCORN_1(A[3],i,j,k) - MACP0A1(pstag,i,jm1mac(j),k,B1)*gdetnosing*dx[2];
	    NOAVGCORN_1(A[1],i,j,k) = 0;
	    NOAVGCORN_1(A[2],i,j,k) = 0;
	  }
	  //copy A[3] -> B[3] before bounding
	  MACP0A1(pstag,i,0,k,B3) = NOAVGCORN_1(A[3],i,0,k);
	}
      }
    }
    //just in case, wait until all CPUs get here
#if(USEMPI)
    MPI_Barrier(MPI_COMM_GRMHD);
#endif
    //bound here
    bound_allprim(STAGEM1,finalstep,t,prim,pstag,ucons, USEMPI);
    //ensure consistency of vpot across the midplane
    if( mycpupos[2] == ncpux2/2 ) {
      for (i=0; i<N1+SHIFT1; i++) {
	for (k=0; k<N3+SHIFT3; k++) {
	  NOAVGCORN_1(A[3],i,0,k) = MACP0A1(pstag,i,-1,k,B3);
	}
      }
    }
  }
  else {
#if( USEMPI )
    for( cj = 0; cj < ncpux2/2; cj++ ) {
      if( mycpupos[2] == cj ){ 
	dj = +1;
	js = 0;
	jsb = 0;
	je = N2;
	jeb = N2 - 1;
	dointegration = 1;
      }
      else if( ncpux2 - mycpupos[2] - 1 == cj ){ 
	dj = -1;
	js = N2;
	jsb = N2-1;
	je = 0;
	jeb = 0;
	dointegration = 1;
      }
      else {
	dointegration = 0;  //skip directly to bounding
      }
      
      if( 1 == dointegration ) {
	//then it's the turn of the current row of CPUs to pick up where the previous row has left it off
	//since pstag is bounded unlike A, use pstag[B3] as temporary space to trasnfer values of A[3] between CPUs
	//initialize lowest row of A[3]
	for (i=0; i<N1+SHIFT1; i++) {
	  for (k=0; k<N3+SHIFT3; k++) {
	    //zero out or copy starting element of vpot
	    if( 0 == cj ) {
	      //if CPU is at physical boundary, initialize (zero out) A[3]
	      NOAVGCORN_1(A[3],i,js,k) = 0.0;
	    }
	    else {
	      //else copy B[3] (which was bounded below) -> A[3]
	      NOAVGCORN_1(A[3],i,js,k) = MACP0A1(pstag,i,jsb-dj,k,B3);
	    }
	    //integrate vpot along the theta line
	    for (j=js; j!=je; j+=dj) {
	      dir=1;
	      get_geometry_gdetonly(i, j-js+jsb, k, FACE1-1+dir, ptrgeomf[dir]);
	      set_igdetsimple(ptrgeomf[dir]);
	      igdetgnosing[dir] = ptrgeomf[dir]->igdetnosing;
	      gdetnosing = 1.0/igdetgnosing[dir];
	      //take a loop along j-line at a fixed i,k and integrate up vpot
	      NOAVGCORN_1(A[3],i,j+dj,k) = NOAVGCORN_1(A[3],i,j,k) + dj * MACP0A1(pstag,i,j-js+jsb,k,B1)*gdetnosing*dx[2];
	    }
	    //copy A[3] -> B[3] before bounding
	    MACP0A1(pstag,i,jeb,k,B3) = NOAVGCORN_1(A[3],i,je,k);
	  }
	}
      }
      //just in case, wait until all CPUs get here
      MPI_Barrier(MPI_COMM_GRMHD);
      //bound here
      bound_allprim(STAGEM1,finalstep,t,prim,pstag,ucons, USEMPI); 
    }
    //ensure consistency of vpot across the midplane
    if( mycpupos[2] == ncpux2/2 ) {
      for (i=0; i<N1+SHIFT1; i++) {
	for (k=0; k<N3+SHIFT3; k++) {
	  NOAVGCORN_1(A[3],i,0,k) = MACP0A1(pstag,i,-1,k,B3);
	}
      }
    }
#endif
  }  
  //need to zero out pstag[B3] everywhere
  FULLLOOP{
    MACP0A1(pstag,i,j,k,B3)=0.0;
    NOAVGCORN_1(A[1],i,j,k) = 0;
    NOAVGCORN_1(A[2],i,j,k) = 0;
  }
  return(0);
}

int normalize_field_local_nodivb(FTYPE targbeta, FTYPE rhomax, FTYPE amax, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], 
				 FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], 
				 FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], 
				 FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  FTYPE compute_rat(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], 
		    FTYPE rhomax, FTYPE amax, FTYPE targbeta, int loc, int i, int j, int k);
  int i,j,k;
  FTYPE ratc_ij, ratc_im1j, ratc_ijm1;
  FTYPE ratf1_ij, ratf2_ij;
  int finalstep;
  
  if(steppart==TIMEORDER-1) finalstep=1; else finalstep=0;
  
  //FULLLOOP{
    //THIS CHANGES IC's!
    //firstly, decrease B2 by 2x -- this leads to a more uniform final beta distribution
    //MACP0A1(prim,i,j,k,B2) *= 0.5;
    //only need to do so on centered fields that bsq depends on
  //}
  
  bound_allprim(STAGEM1,finalstep,t,prim,pstag,ucons, USEMPI);
  
  //Now rescale staggered field components (ucons)
  //Only need to rescale B1cons since B2 will be reconstructed only from B1cons by integrating up vector potential
  //ZLOOP{
  ZSLOOP(0,N1-1+SHIFT1,0,N2-1,0,N3-1) {
    //cell centered ratio in this cell
    ratc_ij   = compute_rat(prim, A, rhomax, amax, targbeta, CENT, i, j, k);
    
    //and in two neighboring cells
    ratc_im1j = compute_rat(prim, A, rhomax, amax, targbeta, CENT, im1mac(i), j, k);
    ratc_ijm1 = compute_rat(prim, A, rhomax, amax, targbeta, CENT, i, jm1mac(j), k);
	
    //ratios centered at FACE1 and FACE2, respectively:
    ratf1_ij = 0.5*(ratc_im1j + ratc_ij);
    ratf2_ij = 0.5*(ratc_ijm1 + ratc_ij);

    // normalize primitive
    // DON'T DO THIS SO AS NOT TO INTRODUCE ORDER DEPENDENCE IN LOOP
    //MACP0A1(prim,i,j,k,B1) *= ratc_ij;
    //MACP0A1(prim,i,j,k,B2) *= ratc_ij;
    //MACP0A1(prim,i,j,k,B3) *= ratc_ij;
    
    // normalize conserved quantity
    //MACP0A1(ucons,i,j,k,B1) *= ratf1_ij;
    //MACP0A1(ucons,i,j,k,B2) *= 0*ratf2_ij;
    //MACP0A1(ucons,i,j,k,B3) *= 0*ratc_ij;
    
    // normalize staggered field primitive
    if(FLUXB==FLUXCTSTAG){
      MACP0A1(pstag,i,j,k,B1) *= ratf1_ij;
    //  MACP0A1(pstag,i,j,k,B2) *= 0*ratf2_ij;
    //  MACP0A1(pstag,i,j,k,B3) *= 0*ratc_ij;  //since assuming axisymmetry, ratc_ij is good enough
    }
    
    // normalize higher-order field
    // SASMARK: have no clue what Bhat is and where in the cell it is located
    // SASMARK: therefore, will normalize it as if it were at the cell center
    //if(HIGHERORDERMEM){
    //  MACP0A1(Bhat,i,j,k,B1) *= ratc_ij;
    //  MACP0A1(Bhat,i,j,k,B2) *= 0*ratc_ij;
    //  MACP0A1(Bhat,i,j,k,B3) *= 0*ratc_ij;      
    //}
  }
  
  //bound_allprim(STAGEM1,t,prim,pstag,ucons, 1, USEMPI);

  return(0);
}

//Returns: factor to multiply field components by to get the desired
//value of beta: targbeta = p_g/p_mag
FTYPE compute_rat(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
		  FTYPE rhomax, FTYPE amax, FTYPE targbeta, int loc, int i, int j, int k)
{
  FTYPE bsq_ij,pg_ij,beta_ij,rat_ij;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM], V[NDIM];
  FTYPE rat, ratc;
  FTYPE profile;
  FTYPE r;

  get_geometry(i, j, k, loc, ptrgeom);
  coord(i,j,k,loc,X);
  bl_coord(X,V);
  r = V[1];
  
  if( bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1 ){
    FAILSTATEMENT("init.sashatorus.c:compute_rat()", "bsq_calc()", 1);
  }
  pg_ij=pressure_rho0_u_simple(i,j,k,loc,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
  beta_ij=2*pg_ij/(bsq_ij+SMALL);
  rat_ij = sqrt(beta_ij / targbeta); //ratio at CENT
  
  //rescale rat_ij so that:
  // rat_ij = 1 inside the main body of torus
  // rat_ij = 0 outside of main body of torus
  // rat_ij ~ rho in between
  //ASSUMING DENSITY HAS ALREADY BEEN NORMALIZED -- SASMARK
  //profile = ( r*r*MACP0A1(prim,i,j,k,RHO)/578.36728604946757 - 0.05 ) / 0.1;
  //profile = ( pow(r,aphipow)*MACP0A1(prim,i,j,k,RHO)/amax - 0.001 ) / 0.05;
  profile = ( log10(pow(r,aphipow)*MACP0A1(prim,i,j,k,RHO)/amax+SMALL) + 3 ) / 1.0;
  //profile = ( MACP0A1(prim,i,j,k,RHO)/rhomax - 0.05 ) / 0.1;
  //profile = ( sqrt(NOAVGCORN_1(A[3],i,j,k)/amax) - 0.05 ) / 0.1;
  if(profile<0) profile = 0;
  if(profile>1) profile = 1;  
  rat_ij *= profile;
  
  return(rat_ij);
}

FTYPE get_maxprimvalrpow(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE rpow, int pl )
{
  int i,j,k;
  FTYPE X[NDIM], V[NDIM];
  FTYPE r;
  
  FTYPE val;
  FTYPE maxval = -VERYBIG;
  
  ZLOOP {
    coord(i,j,k,CENT,X);
    bl_coord(X,V);
  
    r = V[1];
    val = pow(r,rpow)*MACP0A1(prim,i,j,k,pl);
    if ( val > maxval)      maxval = val;
  }
  
  mpimax(&maxval);
  
  
  return(maxval);
  
}


FTYPE get_maxval(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir )
{
  int i,j,k;
  
  FTYPE maxval = -VERYBIG;
  
  ZLOOP {
    if (NOAVGCORN_1(A[dir],i,j,k) > maxval)      maxval = NOAVGCORN_1(A[dir],i,j,k);
  }
  
  mpimax(&maxval);
 
  
  return(maxval);
  
}

FTYPE get_minval(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir )
{
  int i,j,k;
  
  FTYPE minval = VERYBIG;
  
  ZLOOP {
    if (NOAVGCORN_1(A[dir],i,j,k) < minval)      minval = NOAVGCORN_1(A[dir],i,j,k);
  }
  
  mpimin(&minval);
  
  return(minval);
  
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

  //nothing to do here: normalization moved to init_vpot2field_user()
  
  //funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
  //if(funreturn!=0) return(funreturn);
 
  return(0);

}

// assumes normal field definition
int normalize_field_diskonly(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  
#if(DO_NORMALIZE_FIELD)
  funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
  if(funreturn!=0) return(funreturn);
#endif
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
  else if(WHICHPROBLEM==NSTAR){
    atmospheretype=4;
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


// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
// put as function because then wrap-up OpenMP stuff
void add_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  
  
#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOPBONLY(pl) MACP0A1(dest,i,j,k,pl)+=MACP0A1(source,i,j,k,pl);
      
    }// end 3D loop
    
    
  }// end parallel region
  
}

// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
// put as function because then wrap-up OpenMP stuff
void null_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  
  
#pragma omp parallel 
  {
    int i,j,k,pl,pliter;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke);
    
#pragma omp for schedule(OPENMPFULLNOVARYSCHEDULE())
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      
      //      COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOPBONLY(pl) MACP0A1(dest,i,j,k,pl)=0;
      
    }// end 3D loop
    
    
  }// end parallel region
  
}

// general purpose copy machine for 3D arrays with only size NPR appended onto the end of array
// put as function because then wrap-up OpenMP stuff
void add_3d_fieldonly_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  
  
  add_3d_fieldonly(is, ie, js, je, ks, ke, source, dest);
  
  
}

void null_3d_fieldonly_fullloop(FTYPE (*dest)[NSTORE2][NSTORE3][NPR])
{
  int is=-N1BND;
  int ie=N1-1+N1BND;
  int js=-N2BND;
  int je=N2-1+N2BND;
  int ks=-N3BND;
  int ke=N3-1+N3BND;
  
  
  null_3d_fieldonly(is, ie, js, je, ks, ke, dest);
  
  
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

void adjust_fluxctstag_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
{
  
  user1_adjust_fluxctstag_emfs(fluxtime, prim, Nvec, fluxvec);

}

FTYPE get_omegaf_phys(FTYPE t, FTYPE dt, FTYPE steppart)
{
  extern FTYPE t_transition;
  FTYPE get_omegaf_prefactor( FTYPE t_transition, FTYPE t, FTYPE dt, FTYPE steppart );
  
  return( global_OmegaNS * get_omegaf_prefactor( t_transition, t, dt, steppart ) );
  
}

FTYPE get_omegaf_code(FTYPE t, FTYPE dt, FTYPE steppart)
{
  extern FTYPE t_transition;
  FTYPE get_omegaf_prefactor( FTYPE t_transition, FTYPE t, FTYPE dt, FTYPE steppart );
  FTYPE dxdxp[NDIM][NDIM];
  
  //assume dxdxp[3][3] is independent of location
  dxdxprim_ijk(0, 0, 0, CENT, dxdxp);
  
  return( get_omegaf_phys( t, dt, steppart ) / dxdxp[3][3] );
  
}



