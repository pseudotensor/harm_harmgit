// uses init.fishmon.h and bounds.fishmon.c
// See also coord.c, init.tools.c, and restart.c if changes to header/dumpfile
// For thickdisk7->fishmon(new):
// a) must change TRACKVPOT->0.  If doing zakamskabig,then use some large-box versions of parameters given by WHICHPROBLEM==THICKDISK case.
// b) Also, if restarting from old zakamskabig, must truncate read restart of normal dump file, (upper pole file should zero out values as desired), NUMFAILFLOORFLAGS-> (NUMFAILFLOORFLAGS-2) upon read, NUMDUMPTYPES->NUMDUMPTYPES-1
// c) Also, set makehead.inc USETACCLONESTAR4->1


/* 
 *
 * Mark: Tilted BP problem 
 * 
 *
 * 
 *
 */

#include "decs.h"

#define SLOWFAC 1.0		/* reduce u_phi by this amount */
#define MAXPASSPARMS 10

//#define THETAROTMETRIC (0.5*0.7)
#define USER_THETAROTMETRIC (0.05) // arctan(0.2) = 0.19739556
#define USER_THETAROTPRIMITIVES (0.0) // probably want to choose 0, so initial conditions are as if no tilt

#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define THICKDISK 3 // very similar to NORMALTORUS
#define THINBP 4

#define WHICHPROBLEM THINBP

static FTYPE RHOMINEVOLVE,UUMINEVOLVE;
static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin,rinfield; // OPENMPMARK: Ok file global since set as constant before used

static FTYPE rhodisk;


static FTYPE taper_func_exp(FTYPE R,FTYPE rin, FTYPE POTENTIALorPRESSURE ) ;  // MAVARA added June 3 2013

static FTYPE nz_func(FTYPE R) ;   // MARKNOTE torus calculation
static FTYPE taper_func(FTYPE R,FTYPE rin, FTYPE rpow) ;
static FTYPE taper_funcB(FTYPE R,FTYPE rin, FTYPE rpow) ;
static FTYPE integrate_vpot_r(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int i, int j, int k) ;

int calc_da3vsr(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]); // For setting beta const as a function of r at initial conditions
static SFTYPE *da3vsr;
static SFTYPE *da3vsr_tot;
static SFTYPE *tempstore;
static SFTYPE *tempstore_tot;
static SFTYPE *da3vsr_integrated;

FTYPE normglobal;
int inittypeglobal; // for bounds to communicate detail of what doing



int prepre_init_specific_init(void)
{
  int funreturn;
  

  /////////////////////
  //PHI GRID SETUP
  /////////////////////

  dofull2pi = 1;   // MAVARANOTE: do less than full phi    ; This line and the next taken from init.sashatorus.c
  
  global_fracphi = 1.0; //0.25;   //phi-extent measured in units of 2*PI, i.e. 0.25 means PI/2; only used if dofull2pi == 0

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

  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  if(PRODUCTION==0){
    binaryoutput=TEXTOUTPUT;
  }

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
  else if(WHICHPROBLEM==THINBP){
    // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
    h_over_r=0.1;
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
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;  // MAVARA pl = "primative list"


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
  else if(WHICHPROBLEM==THINBP){
    cour=0.8;
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
#elif(WHICHPROBLEM==THINBP)
  // define coordinate type
  defcoord = BPTHIN1;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first

  a = 0.5; //9375 ;

  

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
#elif(WHICHPROBLEM==THINBP)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.; // LR -.3 sometimes
  Rout = 40.0;
  Rin=1.4; //LR .75
  if(totalsize[1]<32) Rout=50.0;
  else if(totalsize[1]<=64) Rout=1.E3;
  else Rout=1.E5;
  Rout=1.E4; // MAVARA tilttest used this
#endif

 
  Rhor=rhor_calc(0);

  //  hslope = 0.3;
  //  hslope = 1.04*pow(h_over_r,2.0/3.0);
  hslope = h_over_r;


  //setRin_withchecks(&Rin);


  
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
  //FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;

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
#elif(WHICHPROBLEM==THINBP)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=1E2; // may have to make smaller if problems
  BSQOULIMIT=1E30;
  UORHOLIMIT=1E2;
  //  UORHOLIMIT=10.0;
  // JCM: Have to choose below so that Mdot from atmosphere is not important compared to true Mdot for thin disk.
  RHOMIN = 1.E-4;
  UUMIN = 1E-6;
  RHOMINEVOLVE = 1.E-4;
  UUMINEVOLVE = 1E-6;

  cooling=COOLUSER; // MARKTODO should override these values set in initbase, right?
  // cooling=NOCOOLING;
  gam=4./3.;

#elif(WHICHPROBLEM==THICKDISK)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=50.0; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=50.0;
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
  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=2000.0;

  // ener period
  DTdumpgen[ENERDUMPTYPE] = 2000.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 4.0; // was 5 after 2.0
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 16000; // was 1000
  DTfake=MAX(1,DTr/10);


#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
/* output choices */
  tf = 4000.0;
#elif(WHICHPROBLEM==THINBP)
/* output choices */
  tf = 1000000.0;
#elif(WHICHPROBLEM==THICKDISK)
  /* output choices */
  tf = 1.3E4*2.0;
#elif(WHICHPROBLEM==GRBJET)
  /* output choices */
  tf = 5E5;
  
  DTd = 250.;                 /* dumping frequency, in units of M */
  DTavg = 250.0;
  DTener = 20.0;                       /* logfile frequency, in units of M */
  DTi = 50.0;                 /* image file frequ., in units of M */
  DTdebug = 250.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 1000;                  /* restart file period in steps */
  DTfake=MAX(1,DTr/10);
#endif



  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *pr)
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


  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);


  // defaults
  beta = 1.e2 ;
  randfact = 4.e-2;



#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
#elif(WHICHPROBLEM==THINBP)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
  rinfield = 10.0;
  beta = 200.;
  randfact = 10.e-2; //4.e-2;
  //  fieldnormalizemin = 3. * Risco;
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
  //check_spc_singularities_user();
  trifprintf("END check_spc_singularities_user\n");
  dualfprintf(fail_file,"WARNUNG: done with check_spc_singularities_user(), but it sometimes stalls or goes very very slow for no apparently good reason.  E.g., on NAUTILUS with -O0, very slow checks.  But just putting dualfprintf before and after the above call leads to quick finish.");

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
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
  int init_dsandvels_bpthin(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

  // assume inittype not used, pos==CENT, and time doesn't matter (e.g. only used at t=0)

#if(WHICHPROBLEM==NORMALTORUS||WHICHPROBLEM==THICKDISK)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINBP)
  return(init_dsandvels_bpthin(whichvel, whichcoord,  i,  j,  k, pr, pstag));
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
    pr[UU] = u * (1. + randfact * (ranc(0,0) - 0.5));
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

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H)) * taper_func(R,rin,-1.0) ;
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



int init_dsandvels_bpthin(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
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
  //FTYPE dxdxp[NDIM][NDIM];
#define UGPOW_MACRO (-.6)
  FTYPE UGPOW=UGPOW_MACRO ;



  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  ph=V[3];


  /* if(j == totalsize[2]/ 2){

    FTYPE dr, dth, dph, dR, dTH, dPH;
    FTYPE dxdxp[NDIM][NDIM];
    dxdxprim_ijk(i, j, k, CENT, dxdxp);
    dr=dx[1]*dxdxp[1][1]; // delta(r) = dx[1]<this is width in cartesian grid> * <scaling to physical boyer-lindquist>dxdp[][]
    dth=dx[2]*dxdxp[2][2];                    // just chain rule:  dx' = dv/dx
    dph=dx[3]*dxdxp[3][3];

    dR=dr;
    dTH=dth*r;
    dPH=dph*sin(th)*r;

    trifprintf("At r = %g  Ratios: %g (dR) %g %g",r,dR, dTH/dR, dPH/dR);
    }
  */

  /* region outside disk */
  R = r*sin(th) ;

  /*if(k==0 && j == totalsize[2]/2){
    trifprintf("RADIUS cents: i=%g , r = %g \n",i,r);
    }*/

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
    cs = H*nz*sqrt(gam) ; // To understand this factor of gam, and how this equation is approximate, consult Eqn. 7.43 of Melia and definition of T~cs^2/gam

    rho = (S/sqrt(2.*M_PI*H*H)) * (pow(R/rin,3./2+UGPOW)) * exp(-z*z/(2.*H*H))  * taper_func(r,rin, 3.0) ; //3 is best i think //* ( 1.0 - R*R/((pow(R,1.5) + a)*(pow(R,1.5) + a)) )  ;// taper_func(R,rin,-1.0) ;
    u = rho*cs*cs/(gam - 1.)/gam;
    ur = 0. ;
    uh = 0. ;
    up = 1./(pow(r,1.5) + a) ;     // MARKNOTE  angular, not linear
    // solution for 3-vel

    
    pr[RHO] = rho ;
    pr[UU] = u * (1. +randfact * (ranc(0,0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up ; //* (1. +0.10 * (ranc(0,0) - 0.5));

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
#define DISKVERTBP 7


#if(WHICHPROBLEM==THICKDISK)
//#define FIELDTYPE TOROIDALFIELD
#define FIELDTYPE DISK2FIELD
#elif(WHICHPROBLEM==THINBP)
#define FIELDTYPE DISKVERTBP
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
  FTYPE idxdxp[NDIM][NDIM];

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
#define RPOW2_MACRO (0.1)  
  FTYPE rpow2=RPOW2_MACRO; //.1 abandoned on april 20, 2014 due to too much flux on BH when taper function in rho removed; //3.0/4.0; // has two meanings. was originally used for, as below, some power outside the break radius. now hijacked for field radial dependence inside rtransition.
  FTYPE FIELDROT=0.0;
  FTYPE hpow=2.0; // MAVARANOTE originally 2.0
#define RBREAK_MACRO (10000000.0)
  FTYPE RBREAK=RBREAK_MACRO;
#define RTRANSITION_MACRO (rin) 
  FTYPE RTRANSITION=RTRANSITION_MACRO; 
  FTYPE UGPOW=UGPOW_MACRO ;
  FTYPE normalize;

  if(l==2){// A_\theta

    r=V[1];
    th=V[2];
    ph=V[3];


    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISK1VERT)||(FIELDTYPE==DISK2VERT)){
      vpot += -(pow(r,rpow)*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph));
    }
    if(FIELDTYPE==DISKVERTBP){
      if(0){
	if(r<RBREAK && r<1.1*rin){
	  vpot += -taper_funcB(1.1*rin,rin,rpow) * (pow(1.1*rin,rpow)*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph));
	}
	else if(r<RBREAK && r>=1.1*rin){
	  vpot += -taper_funcB(r,rin,rpow) * (pow(r,rpow)*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph));
	}
	else if(r>=RBREAK) vpot += -taper_funcB(r,rin,rpow) * ((pow(RBREAK,rpow)/pow(RBREAK,rpow2)*pow(r,rpow2))*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph)); 
      }
      else if(0){
	//dxdxprim_ijk(i, j, k, loc, dxdxp); //FACE1???
	vpot += da3vsr_integrated[startpos[1]+i]; // * dxdxp[2][2]; //should be set to corner! //integrate_vpot_r(prim,i,j,k);
	vpot = vpot * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph);
	//	if(r<RBREAK) vpot += -(pow(r,rpow)*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph));
	//	else if(r>=RBREAK) vpot += -((pow(RBREAK,rpow)/pow(RBREAK,rpow2)*pow(r,rpow2))*pow(sin(th),hpow)*sin(FIELDROT)*sin(ph)); 
      }
      else if(1){
	idxdxprim_ijk(i, j, k, CORN2, idxdxp); // CORN2 for l==2     
	if(r<RTRANSITION) vpot += - (1.+0.*sin(ph*8.))*(1./1.) * pow(tempstore_tot[0]/tempstore_tot[0],(UGPOW/2.+1.5)/1.0) * da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph); //6 is what bptf5 was, etc. //4 works pretty well//MAVARATEMP
	else if(r>=RTRANSITION && r<RBREAK) vpot += -  (1.+0.*sin(ph*8.))*da3vsr_integrated[startpos[1]+i]  * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph);
	else if(r>=RBREAK) vpot += - (1.+0.*sin(ph*8.))*da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph); 
	else vpot += 0.0 ; //-  pow(1.5/tempstore_tot[0],rpow2) * da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph); //MAVARATEMP was 0.0 normally
      }
      else if(0){
	if(r<RTRANSITION) vpot += RTRANSITION/10.;
	else vpot +=  pow(r/500.,4.2) * pow(sin(th),hpow)*sin(FIELDROT)*sin(ph);
      }
      else {
	vpot += 0.0;
      }
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
 
    if(FIELDTYPE==DISKVERTBP){
      //vpot += 0.5*pow(r,rpow)*sin(th)*sin(th) ;
      if(0){
	if(r<RBREAK && r < 1.1*rin){
	  vpot += 0.0 ; //taper_funcB(1.1*rin,rin,rpow) * pow(1.1*rin,rpow)*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
	}
	else if(r<RBREAK && r >= 1.1*rin){
	  vpot += taper_funcB(r,rin,rpow) * pow(r,rpow)*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT)) - taper_funcB(1.1*rin,rin,rpow) * pow(1.1*rin,rpow)*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
	}
	else if(r>=RBREAK) vpot += taper_funcB(r,rin,rpow) * (pow(RBREAK,rpow)/pow(RBREAK,rpow2)*pow(r,rpow2))*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT)) - taper_funcB(1.1*rin,rin,rpow) * pow(1.1*rin,rpow)*pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
      }
      else if(1){
	idxdxprim_ijk(i, j, k, CORN3, idxdxp);
	// vpot += da3vsr_integrated[startpos[1]+i] * idxdxp[3][3] * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT)); // * dxdxp[2][2]; //should be set to corner! //integrate_vpot_r(prim,i,j,k);
	//if(r<=0)	vpot = 0.0;
	//else    
	//vpot = vpot 
	if(r<RTRANSITION) vpot += (1.+0.*sin(ph*8.))*(1./1.) * pow(tempstore_tot[0]/tempstore_tot[0],(UGPOW/2.+1.5)/1.0) * da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT)); //MAVARATEMP
	else if(r>=RTRANSITION && r<RBREAK) vpot += (1.+0.*sin(ph*8.))*da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
	else if(r>=RBREAK) vpot += (1.+0.*sin(ph*8.))*da3vsr_integrated[startpos[1]+i]  * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT)); // was using a .2 multiplyer for sin(ph*8.) term
	else vpot += 0.0; //pow(1.5/tempstore_tot[0],rpow2) * da3vsr_integrated[startpos[1]+i] * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
      }
      else if(0){
	if(r<RTRANSITION) vpot += -RTRANSITION/10.;
	else vpot += - pow(r/500.,4.2) * pow(sin(th),hpow)*(cos(FIELDROT) - cos(ph)*cot(th)*sin(FIELDROT));
      }
      else {
	vpot += 0.0;
      }
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
  else if(FIELDTYPE==DISKVERTBP){
    eqslice=1;
  }
  else{
    eqslice=0;
  }


  parms[0]=rinfield;

  int user2_get_maxes(int eqslice, FTYPE *parms, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min);
  funreturn=user2_get_maxes(eqslice, parms,prim, bsq_max, pg_max, beta_min);

  if(funreturn!=0) return(funreturn);
 
  return(0);
}

/*
// get maximum b^2 and p_g and minimum of beta    MAVARA This is done just as user1_get_maxes except things are normalized only considering values outside a certain minimum radius.
// OPENMPOPTMARK: No benefit from OpenMP due to critical region required and too user specific
int user1_get_maxes2(int eqslice, FTYPE *parms, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min)
{
  int i,j,k;
  FTYPE bsq_ij,pg_ij,beta_ij;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  int gotnormal;
  FTYPE rin;
  int loc;

  
  

  rin=parms[0];


  bsq_max[0] = SMALL;
  pg_max[0]= 0.;
  beta_min[0]=VERYBIG;
  gotnormal=0; // to check if ever was in location where wanted to normalize
  loc=CENT;

  ZLOOP {
    get_geometry(i, j, k, loc, ptrgeom);

    if(eqslice){
      bl_coord_ijk_2(i, j, k, loc, X, V);
      r=V[1];
      th=V[2];
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
	gotnormal=1;
	if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
	if (bsq_ij > bsq_max[0])      bsq_max[0] = bsq_ij;

	pg_ij=pressure_rho0_u_simple(i,j,k,loc,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
	if (pg_ij > pg_max[0])      pg_max[0] = pg_ij;

	beta_ij=pg_ij/(bsq_ij*0.5);

	if (beta_ij < beta_min[0])      beta_min[0] = beta_ij;

      }
    }
    else{
      gotnormal=1;
      if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      if (bsq_ij > bsq_max[0])      bsq_max[0] = bsq_ij;

      pg_ij=pressure_rho0_u_simple(i,j,k,loc,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
      if (pg_ij > pg_max[0])      pg_max[0] = pg_ij;

      beta_ij=pg_ij/(bsq_ij*0.5);
      
      if (beta_ij < beta_min[0])      beta_min[0] = beta_ij;

    }
  }


  mpiisum(&gotnormal);

  if(gotnormal==0){
    dualfprintf(fail_file,"Never found place to normalize field\n");
    if(N2==1 && N3==1){
      ZLOOP {
	bl_coord_ijk_2(i, j, k, loc, X, V);
	dxdxprim_ijk(i, j, k, loc, dxdxp);
	r=V[1];
	th=V[2];

	dualfprintf(fail_file,"i=%d j=%d k=%d V[1]=%21.15g dxdxp[1][1]*dx[1]=%21.15g dx[1]=%21.15g\n",i,j,k,V[1],dxdxp[1][1]*dx[1],dx[1]);
      }
    }
    myexit(111);
  }

  mpimax(bsq_max);
  mpimax(pg_max);
  mpimin(beta_min);


  return(0);

}
*/


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

  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==THICKDISK || WHICHPROBLEM==KEPDISK){
    atmospheretype=1;
  }
  else if(WHICHPROBLEM==GRBJET){
    atmospheretype=2;
  }
  else if(WHICHPROBLEM==THINBP){
    atmospheretype=1;
  }
  else{
    atmospheretype=1; // default
  }

  int user2_set_atmosphere(int atmospheretype, int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
 
  funreturn=user2_set_atmosphere(atmospheretype, whichcond, whichvel, ptrgeom, pr);
  if(funreturn!=0) return(funreturn);
 
  return(0);

}

// UUMIN/RHOMIN used for atmosphere

// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int user2_set_atmosphere(int atmospheretype, int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int pl,pliter;

  // Bondi like initial atmosphere
  //    rho = RHOMIN * 1.E-14;
  //    u = UUMIN * 1.E-14;
  bl_coord_ijk_2(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, V);
  r=V[1];
  th=V[2];

  // default
  PLOOP(pliter,pl) prlocal[pl]=pr[pl];

  if(DOEVOLVERHO){
    // Bondi-like atmosphere
    if(rescaletype==4){
      if(atmospheretype==1){
        // couple rescaletype to atmosphere type
        prlocal[RHO] = RHOMIN*pow(r,-1.5);
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


  if(DOEVOLVEUU){
    // Bondi-like atmosphere
    prlocal[UU]  = UUMIN*pow(r,-2.5);
  }
  else{
    prlocal[UU]  = 0;
  }

    
  // bl-normal observer (4-vel components)
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    prlocal[U1] = -ptrgeom->gcon[GIND(0,1)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    prlocal[U2] = -ptrgeom->gcon[GIND(0,2)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    prlocal[U3] = -ptrgeom->gcon[GIND(0,3)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
  }
  else if(whichvel==VEL3){
    prlocal[U1] = ptrgeom->gcon[GIND(0,1)]/ptrgeom->gcon[GIND(0,0)] ;
    prlocal[U2] = ptrgeom->gcon[GIND(0,2)]/ptrgeom->gcon[GIND(0,0)] ;
    prlocal[U3] = ptrgeom->gcon[GIND(0,3)]/ptrgeom->gcon[GIND(0,0)] ;
    // GAMMIE
    //ur = -1./(r*r);
    //uh=up=0.0;
  }
  else if(whichvel==VELREL4){
    prlocal[U1] = 0.0;
    prlocal[U2] = 0.0;
    prlocal[U3] = 0.0;
  }
  
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

}




int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor, FTYPE *prceiling)
{
  int funreturn;
  
  funreturn=set_density_floors_default(ptrgeom, pr, prfloor, prceiling);

  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  bl_coord_ijk_2(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, V);
  r=V[1];
  th=V[2];

  prlocal[RHO] = RHOMINEVOLVE*pow(r,-1.5);
  prlocal[UU]  = UUMINEVOLVE*pow(r,-2.5);
  

  if(pr[RHO]<prlocal[RHO]) pr[RHO]=prlocal[RHO];
  if(pr[UU]<prlocal[UU]) pr[UU]=prlocal[UU];

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

static FTYPE taper_func_exp(FTYPE R,FTYPE rin, FTYPE POTENTIALorPRESSURE)  // MAVARA added June 3 2013
{
  // POTENTIALorPRESSURE = 1. if used in setting pressure, 2. if used for potential
  FTYPE softer = .3; // 2.0 works ok for resolution of 64 cells for inner 40 M radius
  if(R <= rin)
    return(0.) ;
  else
    return(1. - POTENTIALorPRESSURE * exp((0.95*rin - R)*softer)) ;

}

static FTYPE taper_func(FTYPE R,FTYPE rin, FTYPE rpow)
{
  if(R < rin)
    return(0.) ;
  else
    return(pow( 1.-pow(0.95*rin/R,rpow) , 1.)) ; // was  3 and 2 for powers before 11/10/2012 MAVARA
}

static FTYPE taper_funcB(FTYPE R,FTYPE rin, FTYPE rpow)
{
  if(R < 1.1*rin)
    return(0.) ;
  else
    return(1. - (rpow/(rpow-3.0))*rin*rin*rin / pow(R-0.1*rin,3.0)) ;
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

void adjust_fluxctstag_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR])
{
  // not used
}




// User's cooling function:

#define USERTHETACOOL       (0.1)	/* should be same as h_over_r */
#define USERTAUCOOL         (2.0*M_PI*0.1)	        /* cooling time in number of rotational times : really USERTAUCOOL=2*M_PI would be 1 rotational time */
#define USERNOCOOLTHETAFACT     (1.0)           /* this times h_over_r and no more cooling there*/

// This implementation of cooling is as in Noble et. al. 2009, to simulate a radiative cooling source term which keeps the disk thin to a target H/r

int coolfunc_user(FTYPE h_over_r, FTYPE *pr, struct of_geom *geom, struct of_state *q,FTYPE (*dUcomp)[NPR])
{
        FTYPE X[NDIM],V[NDIM],r,th,R,Wcirc,cs_circ,rho,u,e,P,w,wcirc,dUcool;
	FTYPE taper0;
	int ii,jj, kk, pp;
	FTYPE pressure;
	FTYPE enk, enk0;
        FTYPE Tfix;
	FTYPE Yscaling;

        FTYPE rpho;      
	FTYPE photoncapture;
	FTYPE rincool;
	FTYPE nocoolthetafactor,thetacool,taucool;
	FTYPE bsq_ijcool;


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
	e = u/rho ;   // specific internal energy density
	P=pressure_rho0_u_simple(ii,jj,kk,pp,rho,u);
        w = rho + u + P ; // enthalpy?

	if (bsq_calc(pr, geom, &bsq_ijcool) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
        bl_coord_ijk_2(ii,jj,kk,CENT,X, V) ;
	r=V[1];
	th=V[2];
 
      	rpho=2.0*(1.0+cos(2.0/3.0*(acos(-a))));
	/*
	if(ii==0 && jj==0 && kk==0 && 1) printf("in cooling function: ");
	if(bsq_ijcool/rho > BSQORHOLIMIT || bsq_ijcool/u > BSQOULIMIT || u/rho > UORHOLIMIT){
	  printf("in cooling: bsq/rho=%21.15g, at ijk= %d,%d,%d \n",bsq_ijcool/rho,ii,jj,kk);
	  printf("in cooling: bsq/u=%21.15g \n",bsq_ijcool/u);
	  printf("in cooling: u/rho=%21.15g \n",u/rho);
	}
	*/

	//	trifprintf("rphoton=%lf\n", rpho);
	if(1 || r>rpho){ //SASMARK: cool always, including inside photon orbit
	  photoncapture=1.0 ;
	  //  trifprintf("r=%lf, photoncapture=%lf, rph=%lf \ n", r, photoncapture, rpho); 
	}
	else{
	  photoncapture=0.0 ;
	}


	//if(r>=Rhor)
	//  R = r*sin(th) ;     // r in Noble paper
	//else
	//  R = Rhor; // Don't want target temperature to be overly high in column above and below BH where R<<Rhor    
        R = r*sin(th) ;     // r in Noble paper

	rincool=10.;
        /* crude approximation */
        Wcirc = 1./(a + pow(r,1.5)) ;   // Omega in Noble paper
        cs_circ = thetacool/sqrt(R) ;
	//        wcirc = rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;

        wcirc =   rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;

	Tfix = (thetacool * R * Wcirc) * (thetacool * R * Wcirc);
	Yscaling = (gam-1.)*e/(Tfix);


	//R=r*sin(th); //revert to regular form for everything other than where Tfix/Wcirc has previous adjusted version internally

	if(t > 0. && Yscaling > 1.0 ) { //&& bsq_ijcool*bsq_ijcool*.5/(gam-1)/u >= 0.005) { 
	  if(R*R*h_over_r*h_over_r/1./(bsq_ijcool/rho) > taucool*taucool/Wcirc/Wcirc ) {    
	    dUcool = - (Wcirc/taucool) * rho*Tfix/(gam-1.) * pow( Yscaling - 1.,1.) * photoncapture * q->ucon[TT]  ; 
	  }
	  else if( (h_over_r*R/1.)/sqrt(bsq_ijcool/rho) > 2.*dt){
	    dUcool = - (sqrt(bsq_ijcool/rho)/(h_over_r*R/1.)) * rho*Tfix/(gam-1.) * pow( Yscaling - 1.,1.) * photoncapture * q->ucon[TT]  ;
	  }
	  else{
	    dUcool = - (1./(2.*dt)) * rho*Tfix/(gam-1.) * pow( Yscaling - 1.,1.) * photoncapture * q->ucon[TT]  ;
	  }
	}
	else if(0 && t > 0. && Yscaling > 1.0 && bsq_ijcool*bsq_ijcool*.5/(gam-1)/u < 0.005) { // MAVARANOTE This ceiling has been moved to fixup.c so more immediate 
	  dUcool = - rho*(e - 1.*Tfix/(gam-1)) / dt * photoncapture * q->ucon[TT]  ;
	} 
	else{
	  dUcool = 0. ;
	}

	

	/*
	if(t > 0. && Yscaling > 1.0 ) {

	  dUcool = - (Wcirc/taucool) * u * sqrt( Yscaling - 1.) * photoncapture * q->ucon[TT]  ; // MAVARA temporarily added 0.1 factor to slow cooling to see if it makes a difference on 11/10/2013
	  //	  dUcool=-u*(Wcirc/taucool)*log(enk/enk0)*photoncapture;
	}
        else{
	    dUcool = 0. ;
	}
	*/

	dUcomp[RADSOURCE][UU]=dUcool*(q->ucov[TT]);
	dUcomp[RADSOURCE][U1]=dUcool*(q->ucov[RR]);
	dUcomp[RADSOURCE][U2]=dUcool*(q->ucov[TH]);
	dUcomp[RADSOURCE][U3]=dUcool*(q->ucov[PH]);

	//			trifprintf("ducomps are %g %g %g %g \n", dUcomp[RADSOURCE][UU], dUcomp[RADSOURCE][U1], dUcomp[RADSOURCE][U2],	dUcomp[RADSOURCE][U3]); 
        return(0) ;
}




// get maximum b^2 and p_g and minimum of beta
// OPENMPOPTMARK: No benefit from OpenMP due to critical region required and too user specific
int user2_get_maxes(int eqslice, FTYPE *parms, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min)
{
  int i,j,k;
  FTYPE bsq_ij,pg_ij,beta_ij;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  int gotnormal;
  FTYPE rinlocal;
  int loc;

  
  

  rinlocal=parms[0];


  bsq_max[0] = SMALL;
  pg_max[0]= 0.;
  beta_min[0]=VERYBIG;
  gotnormal=0; // to check if ever was in location where wanted to normalize
  loc=CENT;

  ZLOOP {
    get_geometry(i, j, k, loc, ptrgeom);

    if(eqslice){
      bl_coord_ijk_2(i, j, k, loc, X, V);
      r=V[1];
      th=V[2];
      //      if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      //trifprintf("betaforcheck = %e \n",bsq_ij);

      if(j == 0 && (r<0.9*120.)&&(r>100./2.0)&&(fabs(th-M_PI*0.5)<10.0*M_PI*dx[2]*hslope)){
        gotnormal=1;
        if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
        if (bsq_ij > bsq_max[0])      bsq_max[0] = bsq_ij;

        pg_ij=pressure_rho0_u_simple(i,j,k,loc,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
        if (pg_ij > pg_max[0])      pg_max[0] = pg_ij;

        beta_ij=pg_ij/(bsq_ij*0.5);
	printf("betaforcheck = %f \n",beta_ij);
        if (beta_ij < beta_min[0])      beta_min[0] = beta_ij;

	}
    }
    else{
      gotnormal=1;
      if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      if (bsq_ij > bsq_max[0])      bsq_max[0] = bsq_ij;

      pg_ij=pressure_rho0_u_simple(i,j,k,loc,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
      if (pg_ij > pg_max[0])      pg_max[0] = pg_ij;

      beta_ij=pg_ij/(bsq_ij*0.5);
      
      if (beta_ij < beta_min[0])      beta_min[0] = beta_ij;
	 
    }
  }


  mpiisum(&gotnormal);

  if(gotnormal==0){
    dualfprintf(fail_file,"Never found place to normalize field\n");
    if(N2==1 && N3==1){
      ZLOOP {
        bl_coord_ijk_2(i, j, k, loc, X, V);
        dxdxprim_ijk(i, j, k, loc, dxdxp);
        r=V[1];
        th=V[2];

        dualfprintf(fail_file,"i=%d j=%d k=%d V[1]=%21.15g dxdxp[1][1]*dx[1]=%21.15g dx[1]=%21.15g\n",i,j,k,V[1],dxdxp[1][1]*dx[1],dx[1]);
      }
    }
    myexit(111);
  }

  mpimax(bsq_max);
  mpimax(pg_max);
  mpimin(beta_min);


  return(0);

}



FTYPE integrate_vpot_r(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int i, int j, int k)
{
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM],V[NDIM],r,th,ph;
  FTYPE dxdxp[NDIM][NDIM];
  int loc;
  FTYPE dr, dth, dph, dR, dTH, dPH;
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE ucon[NDIM];
  int ir;

  loc = CENT;

  FTYPE vpottemp = 0.0;

  for(ir=0; ir<=i; ir++){
    coord(ir, j, k, CENT, X);
    bl_coord(X, V);
    //bl_coord_ijk_2(ir, j, k, loc, X, V);
    get_geometry(ir, j, k, loc, ptrgeom);
    r=V[1];
    th=V[2];
    ph=V[3];
    dxdxprim_ijk(ir, j, k, CENT, dxdxp);
    dr=dx[1]*dxdxp[1][1]; // delta(r) = dx[1]<this is width in cartesian grid> * <scaling to physical boyer-lindquist>dxdp[][]
    dth=dx[2]*dxdxp[2][2];                    // just chain rule:  dx' = dv/dx
    dph=dx[3]*dxdxp[3][3];

    dR=dr;
    dTH=dth*r;
    dPH=dph*sin(th)*r;
    
    ucon_calc(MAC(prim,ir,j,k),ptrgeom,ucon, others);

    if(r<=rin) vpottemp += 0.0;
    else { 
      if(prim[ir-2][j][k][RHO] < prim[ir][j][k][RHO] / 100. ) vpottemp += 0.0 ;
      else vpottemp += dR * pow(.1,.5) * ucon[TT] * sqrt( (gam - 1.0)*prim[ir][j][k][UU] ) ; 
    }  
  }

  return( vpottemp );
}




int calc_da3vsr(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int ii,jj=0,kk=0;
  FTYPE X[NDIM],V[NDIM],XX[NDIM],VV[NDIM],rr,r,th;
  FTYPE idxdxp[NDIM][NDIM],dxdxp[NDIM][NDIM];
  FTYPE ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE Rtransition = -1.;
  FTYPE UGPOW = UGPOW_MACRO;
  int trackingticker = 0;

  trifprintf("Starting calc_da3vsr \n");
  da3vsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE)); // add test to see if these work                                                                        
  da3vsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
  da3vsr_integrated=(SFTYPE*)malloc((ncpux1*N1)*sizeof(SFTYPE));
  tempstore=(SFTYPE*)malloc(3*sizeof(SFTYPE)); // add test to see if these work                                                                        
  tempstore_tot=(SFTYPE*)malloc(3*sizeof(SFTYPE));
  tempstore[0]=0.0;
  tempstore[1]=0.0;
  tempstore[2]=0.0;
  FTYPE startval;
  FTYPE rpow2=RPOW2_MACRO; 
  FTYPE switchmad1;
  FTYPE switchmad2;

  if(da3vsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open lumvsr_tot memory \n");
    myexit(1);
  }

  printf("Check 1\n");
  //trifprintf("startpos[2]==totalsize[2]/2: %g ",startpos[2]==totalsize[2]/2 );

  for(ii=0; ii<N1*ncpux1; ii++) da3vsr[ii]=0.;

  for(ii=0; ii<N1; ii++) {                                                                                                                                
                                                                                                                                                          
    if(mycpupos[3] == 0 && (mycpupos[2]==ncpux2/2 && ncpux2>1 || ncpux2==1)){//|| startpos[2]==totalsize[2]/2 ) {  //&& !(startpos[1]==0 && ii==0)                                                                                       
      printf("Getting called 1");

      if(ncpux2==1) jj = N2/2; // MAVARACHANGE?
      //bl_coord_ijk_2(ii, jj, kk, FACE2, X, V);                                
      coord(ii, jj, kk, CORN3, XX); // USING THIS TO MAKE SURE i have the same if statement below as I do when I set the pot.
      coord(ii, jj, kk, FACE2, X);
      bl_coord(X, V);
      bl_coord(XX, VV);
      get_geometry(ii, jj, kk, FACE2, ptrgeom);                                                                             
      //gset_genloc(1,BLCOORDS,ii,jj,kk,FACE2,ptrgeom);
      //gset_genloc(int getprim, int whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom)

      r=V[1];   
      rr=VV[1];
      th=V[2];  

      if(ii==0) startval=rr;

      ucon_calc(MAC(prim,ii,jj,kk),ptrgeom,ucon, others);                                                                                                 
      idxdxprim_ijk(ii, jj, kk, CORN3, idxdxp);                                                                                                           
      dxdxprim_ijk(ii, jj, kk, CORN3, dxdxp);                                                                          
                                 
      if(rr>=RTRANSITION_MACRO && rr < RBREAK_MACRO) {

	switchmad2 = 0.0; //0.5+1.0/M_PI*atan((rr-3000000000000.)/4.); 
	switchmad1 = 1.0; //0.5-1.0/M_PI*atan((rr-30000000000000.)/4.);

	if(tempstore[0] > 0.1 && tempstore[1] < 0.1) {
	  tempstore[1] = rr; 
	  printf("rr tempstore[1]: %21.15g , at tracking= %d \n",rr,trackingticker);
	  trackingticker++;
	}
	if(tempstore[0] <0.1 && startval<=RTRANSITION_MACRO) { //this if is true if r of ii=o for this core was <= RTRANSITION and rr and reached > RTRANSITION as well...only happens on the right single core. //trackingticker == 0 && rr< RBREAK_MACRO/2. && rr<(RTRANSITION_MACRO*1.2) ) {
	  Rtransition=rr;
	  tempstore[0]=rr; //MAVARACHANGE from rr on April 28th so the power law multiple in the A_phi setting section works for r<rtransition
	  tempstore[2] = - (switchmad1*sqrt(1.)+switchmad2*1.)* ptrgeom->gdet * idxdxp[3][3] * dx[1] * ucon[TT] / sqrt(ptrgeom->gcov[GIND(2,2)]) * sqrt( (2./beta) * (gam - 1.0) * MACP0A1(prim,ii,jj,kk,UU) ) ;   // r/(UGPOW + 4.5)^2 was the multiplier
	  printf("rr tempstore[0]: %21.15g , at tracking= %d \n",rr,trackingticker);
	  trackingticker++;
	}

	da3vsr[startpos[1]+ii] = - (switchmad1*sqrt(1.)+switchmad2*1.)* ptrgeom->gdet * idxdxp[3][3] * dx[1] * ucon[TT] / sqrt(ptrgeom->gcov[GIND(2,2)]) * sqrt( (2./beta) * (gam - 1.0) * MACP0A1(prim,ii,jj,kk,UU) ) ; //- ptrgeom->gdet * idxdxp[3][3] * dx[1] / sqrt(ptrgeom->gcov[GIND(2,2)]) * ucon[TT] * sqrt( (2./beta) * (gam - 1.0) * MACP0A1(prim,ii,jj,kk,UU) ) ; //- ptrgeom->gdet * dxdxp[1][1] * dx[1] / sqrt(ptrgeom->gcov[GIND(2,2)]) /r * sqrt( (2./beta) * (gam - 1.0) * MACP0A1(prim,ii,jj,kk,UU) ) ; //- ptrgeom->gdet * dx[1] / sqrt(ptrgeom->gcov[GIND(2,2)]) * pow(.1,.5) * ucon[TT] * sqrt( (gam - 1.0) * MACP0A1(prim,ii,jj,kk,UU) )  ; 

	printf("idxdxp[3][3] checkval: %21.15g , at r= %21.15g \n",idxdxp[3][3],rr);
	printf("idxdxp[2][2] checkval: %21.15g , at r= %21.15g \n",idxdxp[2][2],rr);
	printf("idxdxp[1][1] checkval: %21.15g , at r= %21.15g \n",idxdxp[1][1],rr);
	printf("sqrt(ptrgeom->gcov[GIND(2,2)]) checkval: %21.15g , at r= %21.15g \n",sqrt(ptrgeom->gcov[GIND(2,2)]),rr);
	printf("ucon[TT] checkval: %21.15g , at r= %21.15g \n",ucon[TT],rr);
	printf("ucon[TT] checkval: %21.15g , at r= %21.15g \n",ucon[TT],rr);
	printf("ptrgeom->gdet checkval: %21.15g , at r= %21.15g \n",ptrgeom->gdet,rr);
      }
      else {//if(startpos[2]==totalsize[2]/2) 
	da3vsr[startpos[1]+ii] = 0.0 ;
	printf("Getting called  %d \n", ii);
      } // the if(r < RBREAK_MACRO) above means that when I integrate through to get all the A_3 values at the end of this$
    
      //if(rr>=RTRANSITION_MACRO && da3vsr[0] == 0.) da3vsr[0] = idxdxp[3][3]; // *7*

    }
    else{ 
      da3vsr[startpos[1]+ii] = 0.0;
    }
                                                                                                                                                   
    
                                                                                                                                                          
    //   da3vsr_tot[startpos[1]+ii] = 0.0;                                                                                                                
    // if(integrate(ncpux1*N1,&lumvsr[0],&lumvsr_tot[0],CUMULATIVETYPE,enerregion)>=1) return(1);                                                         
  }                                                                                                                                                    
  
  sleep(10);
  

  //for(ii=0; ii<ncpux1*N1; ii++) da3vsr_tot[ii] = 0.0;
  printf("Check 1\n core id: %d",myid);
  printf("Check 2\n");
  for(ii=0; ii<ncpux1*N1; ii++) printf("valuein %d : %21.15g \n", ii, da3vsr[ii]);
  
  if(integrate(ncpux1*N1,&da3vsr[0],&da3vsr_tot[0], CUMULATIVETYPE3, 0) > 0) trifprintf("Failed to perform integrate() across cpus. \n") ; // 0 is just a filler for an integer$
  sleep(10);
  printf("Check 3\n core id: %d",myid);
  for(ii=0; ii<ncpux1*N1; ii++) printf("valueout %d : %21.15g \n", ii, da3vsr_tot[ii]);
  for(ii=0; ii<ncpux1*N1; ii++) trifprintf("da3vsr_tot %d : %21.15g \n", ii, da3vsr_tot[ii]);

  sleep(10);
  if(integrate(3,&tempstore[0],&tempstore_tot[0], CUMULATIVETYPE3, 0) > 0) trifprintf("Failed to perform integrate() across cpus. \n") ;
  trifprintf("Ending values of tempstore0 and tempstore1 , and tempstore[2] = %21.15g   %21.15g  %21.15g \n", tempstore_tot[0], tempstore_tot[1],tempstore_tot[2]);


  da3vsr_integrated[0] = 1.0 * tempstore_tot[2] / (pow(tempstore_tot[1]/tempstore_tot[0],UGPOW/2+1.5)-1.0) ;

  for(ii=1; ii<N1*ncpux1; ii++) da3vsr_integrated[ii] = da3vsr_integrated[ii-1] + da3vsr_tot[ii-1] ;

  for(ii=0; ii<ncpux1*N1; ii++) trifprintf("value %d : %21.15g \n", ii, da3vsr_integrated[ii]);

  

  // /*became outdated on april 25th 2014 when I added the above sections to exactly match A_phi form on each side of RTRANSITION by exactly matching them at tempstore[0]
  //////////////////////////////////////
  /////////// Finds value of potential at transition radius and sets all values in da3vsr_integrated previous to that to that value  
  ii=0;
  while(da3vsr_integrated[ii]<da3vsr_integrated[0]+.000001 && da3vsr_integrated[ii]>da3vsr_integrated[0]-0.000001) ii++;
  int  ioffirstpast_rtransition=ii-1;
  printf("ioffirstpast_rtransition= %d \n",ioffirstpast_rtransition);
  FTYPE temp1 = da3vsr_integrated[ioffirstpast_rtransition];
  FTYPE temp2 = da3vsr_integrated[ioffirstpast_rtransition+1];
  
// /for(ii=0; ii<=ioffirstpast_rtransition; ii++){ da3vsr_integrated[ii] = da3vsr_integrated[ii] + (1./(pow(tempstore_tot[1]/tempstore_tot[0],.3)-1.0)) * (temp2-temp1) ;} // I don't need to divide by dr since I'm setting the delta A of Ain and Aout equal to each other, so the dr's or dx's would cancel anyway   ;(da3vsr_integrated[ncpux1*N1+1] - da3vsr_integrated[ncpux1*N1]); 
  
  // /for(ii=ioffirstpast_rtransition+1; ii<ncpux1*N1; ii++){ da3vsr_integrated[ii] = da3vsr_integrated[ii] * (1./(pow(tempstore_tot[1]/tempstore_tot[0],.3)-1.0)) ;} // was paired with above line // / etc.*****
  
  printf("value added to pot = %21.15g \n",(1./(pow(tempstore_tot[1]/tempstore_tot[0],rpow2)-1.0)) * (da3vsr_integrated[ioffirstpast_rtransition+1]-da3vsr_integrated[ioffirstpast_rtransition])) ;
  /////*/
  //for(ii=ioffirstpast_rtransition; ii<N1*ncpux1; ii++) da3vsr_integrated[ii] = da3vsr_integrated[ii] + (temp1/.2+temp2);  // MAVARANOTE here I multiple by rpow2 = .1 so that when I divide by it later when using the functional form of vpot_<rtransition the edges at the transition of vpot forms line up
  
  for(ii=0; ii<ncpux1*N1; ii++) if(mycpupos[1]==0 && mycpupos[2]==0) trifprintf("value %d : %21.15g \n", ii, da3vsr_integrated[ii]);
  ///////////

  trifprintf("Ending calc_da3vsr Totalsize[2]=%d \n", totalsize[2]);

  return(0);
  


}
