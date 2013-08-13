#include "decs.h"


//////////////
//
// initial/boundary condition routines for NS-BH problem
//
//////////////





// TODO:
// 1) U7/gdet and B3 don't agree.  Well, they agree outside surface.  But why not inside? A: Because bound B3 but not U7 that just comes from update of fluxes.
// 2) Check why fluxdump out of synch
// 3) Why values for dU are different than B3 for first dB3?
// 4) Check regarding when update_vpot() is called.  Show workflow.  Done: and fixed some optimizaiton issues with fluxctstag
// 5) dU7/gdet is correct sign -- apparent chunkiness at large angle is just plotting (can plot negative and see more symmetric or use symmetric switch in twod.m and careful with image cursor since flipped -- only color is correct).
// 6) F17/gdet makes sense to generate B3 sign everywhere.
// 7) B3 large radius near NS near equator is not right. : More correct at small time.  Maybe already cancelling wave?  Kinda chunky at 45degs at large radius.
// ***8) B3 small radius is mostly wrong sign and very chunky.  : More correct at small time. : 
// 9) uu3 wrong sign for larger radius.  omegaf (both omegaf1 and omegaf2) also wrong sign at larger radius.  Also, both are changing along along NS surface despite fixing PLPR : Was metric terms at small time.
// 10) uu2 correct sign everywhere except right around polar axis of NS -- but sign follows switches of B1,B2 for desired rotation.
// 11) uu1<0 everywhere even at larger radius where should be uu1>0.  Caused by wrong sign for uu3?  Or maybe transient due to field moving near t=0?
// ***12) Figure out why compact for A_3 leads to kinks.  Improve A_{dir} interpolation so better near pole.
// 13) Why does F->U3/omegaf get updated before (1 step off) of EMF->B3?  : Was metric stuff.
//       SUBSTEP0: t=0 so nothing happens
//       SUBSTEP1: t=bit more : EMF/F non-zero on surface of NS, and EMF/F used to update B3stag&B3/U3, so both should be updated

// ***14) uu0 dies on surface with new additions for bound_ (mono and usecompact)
// Seems to be usecompact, death in uu0 on surface at checkaphi 2 (normal DT for dumps) even if mono3=mono4=1 are set.

// 15) Why is B3 in ghost cells being set to 0?  (see EMF @ i=26 j=40 debug code for i=25 and i=26) : Seems was due to usecompact==0
//       Also, Why is p[B3] changing sign with equal magnitude right across the NS surface?  Seems A_i not set well since pstag update leads to jumpiness? .. usecompact=1 in vpot doesn't help.  Better with general usecompact that actually includes A_2 (duh).

// ***16) Use globally fixed curvature sign around NS and convert fixed corners to boundary corners even if jumping by more than 1 cell.  Flux between will be able to move around, but only within bounds of fixed points.
//        But unsure about Ad1,Ad2.  Can skip them (or looks like would skip all of them!) since they carry rotation.
//        Unsure about what can skip or if can skip anything at all.
//     First, try to convert only a few points to BC type, where have problems.  See if even works.  Only modify for A_3 for now.  Some general sense of relaxing A_i when convexifies in perp1i and perp2i plane.
//        Point is that can allow slight velocity *across* surface till that settles down, but need fixed points to nail-down the field so field doesn't find path to slip.

// 17) Evolve A_i (i.e. A_1 and A_2) directly instead of using EMFs for points on NS surface.  More accurate and no dissipation issue.  Like fixing A_3, but fixing all by updating via EMF_i computed on surface from known values.  Actually, with A_1 and A_2 at FACE2 and FACE1 (respectively), since BOUNDPLPR fixes exactly using v,B right there, same as fixing.  In any case, don't really know B^2 at FACE1 and B^1 at FACE2.  Interpolation would be as good as doing now.


// NEW:
// 1) B3 wrong sign at quasi-poles because B1 wrong sign?
// 2) B3 wrong signn for entire inner-radial region at late time. Field lines bent wrong way???!!!  All other v^i and B^p make sense.  At very early times, sign of B3 makes sense (swept back).  All those closed lines (into BH) face forward?  Odd, makes no sense.
// 3) Removal of kink and lower uu0 leads to field reversal.  Also, actually more kinky elsewhere on NS surface at slightly large radii nearer to equator.  And that's where uu0 is slightly high (\gamma\sim 3), whereas without kink-fix (freely floating A3) that part is not slighlty high (gamma\sim 1.4)



/////////////////
//
// NOTE the start-up sequence:
//
// 0) pre-dump) t=0, switches->0, so omegaf=0.  Initially have B3=uu3=0.  Note that vpot,B,etc. are renormalized so BCs will pass through unnormalized fields at first.
// 1) dump0) t=0, switches->0, so omegaf=0.  Initially have B3=uu3=0
// 2) dump1) t=completed first substep, fluxtime=0 during substep so switches still off.  boundtime=some dt
// 3) dump2) t=completed second substep, fluxtime=some dt so switches->on (so omegaf=omegak if hard to full on switch).  boundtime=some next dt.  In BL, uu3>0 on NS surface -> EMF1/2 non-zero -> B3 non-zero.  But also, in KS, uu3>0 means u_1 non-zero even if B3=0 such that u^1=u^2=0 (i.e. stationary condition on  NS surface).  So u.B non-zero, so T^p_\phi non-zero, so uu3 *on grid* becomes non-zero.
// ...
//
//
/////////////////







#define SLOWFAC 1.0  /* reduce u_phi by this amount */
#define MAXPASSPARMS 10


static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;


static FTYPE nz_func(FTYPE R) ;


// from phys.tools.c
extern int OBtopr_general3(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);


static int init_dsandvels_nsbh(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag);



static int check_nsdepth(SFTYPE time);


static int pos_NS(SFTYPE time, FTYPE *V, FTYPE *Vns, FTYPE *Vcartns, FTYPE *absrdiff);

static int setNSparms(SFTYPE time, FTYPE *rns, FTYPE *omegak, FTYPE *omegaf, FTYPE *Rns, FTYPE *Rsoft, FTYPE *v0);

static int is_dir_insideNS(int dir,int i, int j, int k, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom);
static int is_pos_insideNS(int pos,int i, int j, int k, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom, int *faces);
static int is_dir_onactivegrid(int dir, int i, int j, int k);
static int is_ongrid(int dir, int i, int j, int k);


static int get_del(int i, int j, int k, int ii, int jj, int kk, int *delorig, int *del);


static int rescale_pl(SFTYPE time, int dir, struct of_geom *ptrgeom, FTYPE *prin, FTYPE *prout);
static int unrescale_pl(SFTYPE time, int dir, struct of_geom *ptrgeom, FTYPE *prin, FTYPE *prout);
static FTYPE rescale_A(SFTYPE time, int dir, int i, int j, int k);


static int get_fixedvalues(int dir, int i, int j, int k, int ii, int jj, int kk);

static int checkmono4(int firstpos, int lastpos, FTYPE *xpos, FTYPE y0, FTYPE y1, FTYPE y2, FTYPE y3, FTYPE y4);
static int checkmono3(int firstpos, int lastpos, FTYPE *xpos, FTYPE y0, FTYPE y1, FTYPE y2, FTYPE y3);


static void bound_gen_deeppara(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);


static int get_compactvalue(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, int *del, int *usecompactpl, FTYPE *compactvalue, FTYPE *retiref, FTYPE *retjref, FTYPE *retkref, int *retposcorn, int *reticorn, int *retjcorn, int *retkcorn);

static int get_fulldel_interpolation(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, FTYPE iref, FTYPE jref, FTYPE kref, int poscorn, int icorn, int jcorn, int kcorn, int *del, int usecompact, FTYPE *compactvalue, int *isongrid, int *isinside, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom, int *numpointsused, FTYPE *xpos, int *compactpoint, FTYPE *prfulldel, FTYPE *prfulldelreduced);

static int checkucon_modifyuu(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, int iii, int jjj, int kkk, int iiii, int jjjj, int kkkk, FTYPE *xpos, FTYPE localsingle[5][NPR], struct of_geom *(ptrgeom[5]), int usecompact);




#define switcht(t,t1,t2) (t>t1 && t<t2 ? (t-t1)/(t2-t1) : (t<=t1 ? 0.0 : 1.0) )


// loop modelled after COMPLOOPN in global.comploops.h
#define COMPLOOPN3kk(kk) for(kk=0+SHIFTX3DN;kk<=N3-1+SHIFTX3UP;kk++)
#define COMPLOOPN2jj(jj) for(jj=0+SHIFTX2DN;jj<=N2-1+SHIFTX2UP;jj++)
#define COMPLOOPN1ii(ii) for(ii=0+SHIFTX1DN;ii<=N1-1+SHIFTX1UP;ii++)
#define COMPLOOPNiijjkk(ii,jj,kk) COMPLOOPN3kk(kk) COMPLOOPN2jj(jj) COMPLOOPN1ii(ii)




FTYPE normglobal;
int inittypeglobal; // for bounds to communicate detail of what doing


#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define NSBH 3

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


#define WHICHPROBLEM NSBH // choice





int prepre_init_specific_init(void)
{
  int funreturn;
  
  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  // initial value for overall A_i normalization
  normglobal=1.0;
  inittypeglobal=-100;

  if(BOUNDPLPR==0){
    dualfprintf(fail_file,"Must have BOUNDPLPR=1 for NSBH problem\n");
    myexit(869346313);
  }


  return(0);

}


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=1.0;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=M_PI*0.4;

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

  cour=0.8;
  //  fluxmethod= HLLFLUX;

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
  //  defcoord = JET3COORDS;
  defcoord = JET6COORDS;
#elif(WHICHPROBLEM==GRBJET || WHICHPROBLEM==NSBH)
  // define coordinate type
  defcoord = JET4COORDS;
  //  defcoord = LOGRSINTH;
  //  defcoord = JET6COORDS;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first
  // TEST
  //  a = 0.9 ;
  a = 0.0;
  

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  //  Rout = 1E3;
  Rout = 1.3*2*1E4;
  hslope = 1.04*pow(h_over_r,2.0/3.0);
  //  hslope = 0.3;
#elif(WHICHPROBLEM==GRBJET)
  R0 = -3.0;
  Rout = 1E5;
  hslope = 1.04*pow(h_over_r,2.0/3.0);
  //  hslope = 0.3;
#elif(WHICHPROBLEM==NSBH)
  Rout = 2E2;
  // defcoord==JET6COORDS
  // R0=0.0;
  //hslope = 1.04*pow(h_over_r,2.0/3.0);
  // defcoord==JET4COORDS
  R0 = 0.0;
  hslope = 0.3;
#endif

 
  Rhor=rhor_calc(0);



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
  //  FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;

  // TEST
  //  lim[1] = lim[2] = lim[3] = PARALINE;
  lim[1] = lim[2] = lim[3] = MC;

  DOENOFLUX=NOENOFLUX;
  

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==NSBH)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=50.0; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=50.0;

  //  TEST
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
  //RHOMIN = 1E+3;
  //UUMIN = 1E+2;

#elif(WHICHPROBLEM==GRBJET)
  BCtype[X1UP]=FIXEDOUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E3;
  BSQOULIMIT=1E4;
  RHOMIN = 23.0;
  UUMIN = 1.7;
#endif


  if(EOMTYPE==EOMFFDE){
    //GAMMAMAX=2000.0;
    GAMMAMAX=200.0; // problem with catastrohpic cancellation not always allowing gamma=2000
    GAMMAFAIL=10.0*GAMMAMAX; // when we think gamma is rediculous as to mean failure and solution is not accurate.
  }
  




#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==NSBH)
  /* output choices */
  tf = 1.3E4*2.0;

  //#define BASEDT 0.005
  //#define SHORTDT 0.005

  // GODMARK DEBUG TEST TODO
  //#define BASEDT 0.00001
  //#define SHORTDT 0.00001
#define BASEDT 1.0
#define SHORTDT 1.0

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = BASEDT;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= BASEDT;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = SHORTDT;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = SHORTDT;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = BASEDT;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 1000;

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



  // GODMARK TODO DEBUG TEST
  DODIAGEVERYSUBSTEP=0;


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



  ///////////////////
  //
  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  //
  ///////////////////
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);



  ////////////////
  //
  // things to do BEFORE primitives are assigned but AFTER grid is setup
  //
  /////////////////

  if(WHICHPROBLEM==NSBH){
    //////////////////
    //
    // check if NS has sufficient grid depth for boundary conditions to be used for interpolations on smoothish functions
    //
    //////////////////
    check_nsdepth(t); // t is ok here

  }

  

  //////////////////////
  //
  // Set primitives
  //
  //////////////////////


  //  beta = 1.e2 ;
  //  beta = 20.0;
  beta = 10.0*4.0; // 256x128x256
  //  randfact = 4.e-2;
  randfact = 0.1;

#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  //  rin = 15.0 ;
  rin = 10.0 ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET || WHICHPROBLEM==NSBH)
#endif



  //SASMARK restart: need to populate panalytic with IC's
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
    //    MYFUN(init_primitives(prim,pstag,ucons,vpot,Bhat,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);

    // utemparray only used otherwise in advance.c
    MYFUN(init_primitives(panalytic,pstaganalytic,GLOBALPOINT(utemparray),vpotanalytic,Bhatanalytic,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }

  // check rmin
  check_rmin();


  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;
  int inittype;


  inittype=1;
  inittypeglobal=inittype; // required for bound to know that field not set yet
  funreturn=user1_init_primitives(inittype, prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);

  
  if(WHICHPROBLEM==NSBH){
    // need to use field obtained above to re-get behavior inside NS at t=0 so solution stationary one
    // this doesn't require iteration, just that A_\phi and B^i set after density,u, v^i set, even though A_i doesn't depend upon them.
    inittype=2;
    inittypeglobal=inittype; // tells bound that field set already, so can play with density as dependent upon field
    funreturn+=user1_init_primitives(inittype, prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  }


  if(funreturn!=0) return(funreturn);

  return(0);


}




int init_dsandvels(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);


#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==NSBH)
  return(init_dsandvels_nsbh(inittype, pos, whichvel, whichcoord, time, i, j, k, pr, pstag));
#endif

  // WHICHPROBLEM==GRBJET doesn't have initial dsandvels to set beyond atmosphere

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
  //  rmax = 1E2 ;
  rmax = 1E2 ;
  //  rmax = 60 ;
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

    if(FLUXB==FLUXCTSTAG && pstag!=NULL){
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

    if(FLUXB==FLUXCTSTAG && pstag!=NULL){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}


















// NS is axisymmetric or fully non-axisymmetric
#define DONUT 0
#define SPHERICAL 1


#define NSSHAPE DONUT
#define RSEP 6.0
#define BPOLE 1.0

// TEST
// barely can do SIGMA0=50.0 at 128x64 with Rout=2E2 with omegak=0.  Problems with omegak!=0
//#define SIGMA0 50.0
#define SIGMA0 10.0

// sigma to start with for NS pole
// seems as long as only fix EMF (not rest of fluxes) that behavior is good even at high sigma, so don't have to vary this so much
// That's good, so that don't have mass in the way!
//#define SIGMA0 (0.01) // not real sigma when allowing for diffusion out of NS if don't bound plpr at surface
//#define SIGMAT0 (0.01)
#define SIGMAT0 (SIGMA0) // sigma to start with for NS pole

#define HOTSIGMA(sigma) (sigma*10.0)


// get whether inside NS or not  (used by bounds.nsbh.c as well)
int get_insideNS(SFTYPE time, int i, int j, int k, FTYPE *V, FTYPE *Bcon, FTYPE *vcon, FTYPE *inout)
{
  int insideNS;
  FTYPE r=V[RR],th=V[TH],ph=V[PH]; // time so more general
  FTYPE R=r*sin(th);
  FTYPE x,y,z;
  FTYPE xns,yns,zns;
  FTYPE absrdiff;
  FTYPE Vns[NDIM],Vcartns[NDIM];



  // Cartesian position
  // assumes V is spherical polar coordinates
  x = r*cos(ph)*sin(th);
  y = r*sin(ph)*sin(th);
  z = r*cos(th);


  // default is outside
  insideNS=0;

  // get NS position
  pos_NS(time,V, Vns, Vcartns, &absrdiff);

  // get NS parameters
  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  setNSparms(time,&rns, &omegak, &omegaf,  &Rns, &Rsoft, &v0);

  if(absrdiff<=Rns){
    insideNS=1;
  }



  if(Bcon!=NULL){


    // get difference vector for spherical polar coordinates
    //    ncon[TT]=0.0;
    //    ncon[RR]=rns-r;
    //    ncon[TH]=thns-th;
    //    ncon[PH]=phns-ph;

    // Cart orthonormal difference vector from center of NS to V[Cart]
    FTYPE ncart[NDIM];
    ncart[TT]=0.0;
    ncart[1]=(Vcartns[1]-x); // x
    ncart[2]=(Vcartns[2]-y); // y
    ncart[3]=(Vcartns[3]-z); // z

    // get BL geometry in spherical polar coordinates
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    gset(0,BLCOORDS,i,j,k,ptrgeom);


    /////////////////////////////
    // Get Bcart orthonormal
    FTYPE Bcov[NDIM];
    lower_vec(Bcon,ptrgeom,Bcov);
    FTYPE Bsq,absB;
    Bsq=0.0;
    int jj;
    DLOOPA(jj) Bsq+=Bcon[jj]*Bcov[jj];
    absB=sqrt(Bsq);

    // orthonormal Bo^i, such that Bo^i Bo^i = 1
    FTYPE Bo[NDIM];
    Bo[0]=0.0;
    Bo[1]=Bcon[1]/absB*(1.0);
    Bo[2]=(Bcon[2]/absB)*(r);
    Bo[3]=(Bcon[3]/absB)*(r*sin(th));

    // orthonormal Cart Bcart^i
    FTYPE Bcart[NDIM];
    Bcart[TT]=0.0;
    Bcart[1]=Bo[1]*sin(Bo[2])*cos(Bo[3]);
    Bcart[2]=Bo[1]*sin(Bo[2])*sin(Bo[3]);
    Bcart[3]=Bo[1]*cos(Bo[2]);

    /////////////////////////////
    // Get vcart orthonormal
    FTYPE vcov[NDIM];
    lower_vec(vcon,ptrgeom,vcov);
    FTYPE vsq,absv;
    vsq=0.0;
    DLOOPA(jj) vsq+=vcon[jj]*vcov[jj];
    absv=sqrt(vsq);

    // orthonormal vo^i, such that vo^i vo^i = 1
    FTYPE vo[NDIM];
    vo[0]=0.0;
    vo[1]=vcon[1]/absv*(1.0);
    vo[2]=(vcon[2]/absv)*(r);
    vo[3]=(vcon[3]/absv)*(r*sin(th));

    // orthonormal Cart vcart^i
    FTYPE vcart[NDIM];
    vcart[TT]=0.0;
    vcart[1]=vo[1]*sin(vo[2])*cos(vo[3]);
    vcart[2]=vo[1]*sin(vo[2])*sin(vo[3]);
    vcart[3]=vo[1]*cos(vo[2]);


    // get PRIMECOORDS version
    //mettometp(i,j,k,CENT,ncon);
    // get geometry
    //      get_geometry(i, j, k, CENT, ptrgeom);
    //    gset(0,BLCOORDS,i,j,k,ptrgeom);
    // lower
    //    lower_vec(ncon,ptrgeom,ncov);

    FTYPE Binout=0.0;
    DLOOPA(jj) Binout+=Bcart[jj]*ncart[jj];

    FTYPE Vinout=0.0;
    DLOOPA(jj) Vinout+=vcart[jj]*ncart[jj];

    // so now can compare sign of inout with sign of B^i v_i
    if(sign(Binout)==sign(Vinout)) *inout=1.0;
    else *inout=-1.0;

    //    dualfprintf(fail_file,"FOO: %d %d %d : %21.15g\n",i,j,k,*inout);

    // override:
    if(th<M_PI*0.5){
      *inout=-1.0;
    }
    else{
      *inout=+1.0;
    }

  }


  return(insideNS);


}


















// get position of NS
int pos_NS(SFTYPE time, FTYPE *V, FTYPE *Vns, FTYPE *Vcartns, FTYPE *absrdiff)
{
  int insideNS;
  FTYPE r=V[RR],th=V[TH],ph=V[PH]; // time so more general
  FTYPE R=r*sin(th);
  FTYPE x,y,z;
  FTYPE thns,phns;
  FTYPE xns,yns,zns;


  // get NS parameters
  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  setNSparms(time,&rns, &omegak, &omegaf,  &Rns, &Rsoft, &v0);


  // Cartesian position
  // assumes V is spherical polar coordinates
  x = r*cos(ph)*sin(th);
  y = r*sin(ph)*sin(th);
  z = r*cos(th);


  // depends upon shape
  if(NSSHAPE==SPHERICAL){

    /////////////////
    //
    // spherical NS
    //
    /////////////////
    

    thns=0.5*M_PI;

    xns=rns*cos(time*omegak); // and orbits
    yns=rns*sin(time*omegak); // and orbits
    zns=0.0; // no change with time

    Rns=sqrt(xns*xns+zns*zns);
    if(zns!=0.0) thns=atan2(Rns,zns);
    else thns=0.5*M_PI;
    if(xns!=0.0) phns=atan2(yns,xns);
    else if(yns>0.0) phns=0.5*M_PI;
    else if(yns<0.0) phns=1.5*M_PI;
    

    Vcartns[TT]=time;
    Vcartns[1]=xns;
    Vcartns[2]=yns;
    Vcartns[3]=zns;

    // get |rvec-rnsvec|
    *absrdiff=sqrt( (xns-x)*(xns-x) + (zns-z)*(zns-z) + (yns-y)*(yns-y) );

    
  }
  else if(NSSHAPE==DONUT){

    /////////////////
    //
    // donut NS (no change of shape with time)
    //
    /////////////////
    zns=0.0; // no change with time

    thns=0.5*M_PI;
    phns=ph; // never introduce difference in \phi


    xns=R;
    yns=y; // force as if no difference in \phi
    
    Vcartns[TT]=time;
    Vcartns[1]=rns; // quasi-position
    Vcartns[2]=yns;
    Vcartns[3]=zns;

    // get |rvec-rnsvec|  quasi-version for axisymmetry
    *absrdiff=sqrt( (rns-R)*(rns-R) + (zns-z)*(zns-z) + (yns-y)*(yns-y) );

  }


  Vns[TT]=time;
  Vns[1]=rns;
  Vns[2]=thns;
  Vns[3]=phns;



  return(0);

}











#define NUMMASKFLAGS 18

#define NSMASKINSIDE 0
#define NSMASKSHELL 1
#define NSMASKDISTTOSHELL 2
#define NSMASKDISTTOSHELLCORN1 (NSMASKDISTTOSHELL+1) // should be 2+1=3
#define NSMASKDISTTOSHELLCORN2 (NSMASKDISTTOSHELL+2) // should be 2+2=4
#define NSMASKDISTTOSHELLCORN3 (NSMASKDISTTOSHELL+3) // should be 2+3=5
#define NSMASKCLOSEI 6
#define NSMASKCLOSEJ 7
#define NSMASKCLOSEK 8
#define NSMASKCLOSEICORN1 9  //  NSMASKCLOSEICORN1 + 3*(dir-1) + [ijk as CORN_dir-1=0,1,2]
#define NSMASKCLOSEJCORN1 10
#define NSMASKCLOSEKCORN1 11
#define NSMASKCLOSEICORN2 12
#define NSMASKCLOSEJCORN2 13
#define NSMASKCLOSEKCORN2 14
#define NSMASKCLOSEICORN3 15
#define NSMASKCLOSEJCORN3 16
#define NSMASKCLOSEKCORN3 17

int PTRDEFGLOBALMACP0A1(nsmask,N1M,N2M,N3M,NUMMASKFLAGS);
int BASEMACP0A1(nsmask,N1M,N2M,N3M,NUMMASKFLAGS); 




// see if this grid pos (CENT, FACE1, etc.) are on NS surface or not
int is_pos_insideNS(int pos,int i, int j, int k, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom, int *faces)
{
  int dir,odir1,odir2;
  int isinside;




  if(pos==CENT || pos==CORN1 || pos==CORN2 || pos==CORN3){
    // dir is CORN_{dir}
    if(pos==CENT){
      dir=0;
    }
    else if(pos==CORN1){
      dir=1;
    }
    else if(pos==CORN2){
      dir=2;
    }
    else if(pos==CORN3){
      dir=3;
    }

    isinside=is_dir_insideNS(dir,i, j, k, hasmask, hasinside, reallyonsurface, cancopyfrom);

  }
  else if(pos==FACE1 || pos==FACE2 || pos==FACE3){

    dir=pos-FACE1+1; // face dir

    // see if FACEdir is NS boundary
    int im,jm,km;
    im=i-(dir==1)*N1NOT1;
    jm=j-(dir==2)*N2NOT1;
    km=k-(dir==3)*N3NOT1;

    int ip,jp,kp;
    ip=i+(dir==1)*N1NOT1;
    jp=j+(dir==2)*N2NOT1;
    kp=k+(dir==3)*N3NOT1;

    int faceinside,faceshell,faceboth;
    faceinside=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1    && GLOBALMACP0A1(nsmask,im,jm,km,NSMASKSHELL)==1); // |    shell    |face_i    i inNS    |
    faceshell =(GLOBALMACP0A1(nsmask,im,jm,km,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k,NSMASKSHELL)==1);    // |    inNS     |face_i    i shell   |
    faceboth  =(GLOBALMACP0A1(nsmask,im,jm,km,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);   // |    inNS     |face_i    i inNS    |

    // need faces to determine whether should use left or right interpolated primitive
    faces[0]=faceinside;
    faces[1]=faceshell;
    faces[2]=faceboth;

    // translate
    //    isinside=faceboth;
    isinside=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1    && GLOBALMACP0A1(nsmask,im,jm,km,NSMASKINSIDE)==1);
    //    *hasmask=faceinside||faceshell;
    *hasmask=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKSHELL)==1    || GLOBALMACP0A1(nsmask,im,jm,km,NSMASKSHELL)==1);
    //    *hasinside=faceinside||faceshell||faceboth;
    *hasinside=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1    || GLOBALMACP0A1(nsmask,im,jm,km,NSMASKINSIDE)==1);
    *reallyonsurface=(isinside==0 && *hasmask==1 && *hasinside==1);
    //    *reallyonsurface=faceinside||faceshell; // both versions should give same result


    // get odir's
    get_odirs(dir,&odir1,&odir2);

    // for hasinside and cancopyfrom: is extended to include active cells compared to reallyonsurface because otherwise no points to copy from would exist (or compactly exist) in some cases
    int ii,jj,kk;
    for(ii=i-N1NOT1*(odir1==1)-N1NOT1*(odir2==1);ii<=i;ii++){
      for(jj=j-N2NOT1*(odir1==2)-N2NOT1*(odir2==2);jj<=j;jj++){
        for(kk=k-N3NOT1*(odir1==3)-N3NOT1*(odir2==3);kk<=k;kk++){
          if(*hasinside==0 && *hasmask==1 &&
             (GLOBALMACP0A1(nsmask,ii-(dir==1)*N1NOT1                  ,jj-(dir==2)*N2NOT1                  ,kk-(dir==3)*N3NOT1,NSMASKINSIDE)==1)
             //             || (GLOBALMACP0A1(nsmask,ii+(dir==1)*N1NOT1                  ,jj+(dir==2)*N2NOT1                  ,kk+(dir==3)*N3NOT1,NSMASKINSIDE)==1)
             ){
            *hasinside=1;
          }
        }
      }
    }


    *cancopyfrom=(isinside==0 && *hasmask==1 && *hasinside==1);


    // GODMARK TODO SUPERMARK: TO FIX FOR CURVATURE VERSION OF ON SURFACE


  }
  else{
    dualfprintf(fail_file,"is_pos_insideNS() not yet setup for pos=%d ijk=%d %d %d\n",pos,i,j,k);
    myexit(87346623);
  }



  return(isinside);


}



// check if A_i is fully inside NS (else should be fixed or evolved and not extrapolated/interpolated)
// also report if shell is a neighbor
// if A_i is just inside by half-cell (at CENT relative to actual surface), allow hasinside==1 if is mask cell and inside cells are next to it.  So if need to know if A_dir or EMF_dir are actually on the surface, should ask if isinside==0 && hasmask==1 && hasinside==1 && (not case that all 4 cells around corner are INSIDE and not case that all 4 cells around corner are notINSIDE)
int is_dir_insideNS(int dir,int i, int j, int k, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom)
{
  int odir1,odir2;
  int isinside;


  // GODMARK TODO SUPERMARK: TO FIX FOR CURVATURE VERSION OF ON SURFACE


  if(dir==0){
    isinside=GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE);
    *hasmask=GLOBALMACP0A1(nsmask,i,j,k,NSMASKSHELL);
    if(*hasmask==1) *hasinside=1; // for loc=CENT, if on mask, then must be inside type nearby.
    else *hasinside=0;

    *reallyonsurface=0; // never so
    *cancopyfrom=(isinside==0 && *hasmask==1 && *hasinside==1);
  }
  else{

    // get odir's
    get_odirs(dir,&odir1,&odir2);

 
    // Only bound if A_i is such that a boundary cell, which occurs if *all* cells directly *touching* A_i are boundary cells.  Otherwise, should be set as surface value.
    isinside=0;
    if(
       (i>=-N1BND+N1NOT1 && j>=-N2BND+N2NOT1 && k>=-N3BND+N3NOT1)         && // avoid out of bounds
       (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKINSIDE)==1)                                     && 
       (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)                   && 
       (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKINSIDE)==1)                   && 
       (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)
       ){
      isinside=1;
    }


    // see if one of the neighbors is a mask.
    *hasmask=0;
    if(
       (i>=-N1BND+N1NOT1 && j>=-N2BND+N2NOT1 && k>=-N3BND+N3NOT1)         && // avoid out of bounds
       (
        (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKSHELL)==1)                                     ||
        (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKSHELL)==1)                   || 
        (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKSHELL)==1)                   || 
        (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKSHELL)==1)
        )
       ){
      *hasmask=1;

      // DEBUG:
      //      dualfprintf(fail_file,"GOTHASMASK: dir=%d : %d %d %d %d\n",dir,
      //    (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKSHELL)==1),                                     
      //    (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKSHELL)==1)                   ,
      //    (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKSHELL)==1)                   ,
      //    (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKSHELL)==1)
      //    );

    }


    // see if one of the neighbors is inside.
    // for dimensions that don't exist, treated as FACE values
    *hasinside=0;
    if(
       (i>=-N1BND+N1NOT1 && j>=-N2BND+N2NOT1 && k>=-N3BND+N3NOT1)         && // avoid out of bounds
       (
        (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKINSIDE)==1)                                     ||
        (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)                   || 
        (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKINSIDE)==1)                   || 
        (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)
        )
       ){
      *hasinside=1;

      // DEBUG:
      //      dualfprintf(fail_file,"GOTHASINSIDE: dir=%d : %d %d %d %d\n",dir,
      //    (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKINSIDE)==1)                                     ,
      //    (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)                   ,
      //    (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKINSIDE)==1)                   ,
      //    (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)
      //    );


    }


    // indicates really on surface (note: only true for dir==3 with also N3==1
    *reallyonsurface=(isinside==0 && *hasmask==1 && *hasinside==1);



    // check if dir-direction corresponds to surface normal direction (so dir value is treated as loc=CENT quantity)
    // This ensures if mask cell nearby but no inside cell that (e.g. if dir==2 and N3NOT1==0) that know that this is near surface by half-cell
    // This A_dir or EMF_dir can then be used to copy into the ghost cells since it's either fixed (on surface exactly) or active (for some values, no surface values and only active ones can be source)
    int ii,jj,kk;
    for(ii=i-N1NOT1*(odir1==1)-N1NOT1*(odir2==1);ii<=i;ii++){
      for(jj=j-N2NOT1*(odir1==2)-N2NOT1*(odir2==2);jj<=j;jj++){
        for(kk=k-N3NOT1*(odir1==3)-N3NOT1*(odir2==3);kk<=k;kk++){
          if(*hasinside==0 && *hasmask==1 &&
             (GLOBALMACP0A1(nsmask,ii-(dir==1)*N1NOT1                  ,jj-(dir==2)*N2NOT1                  ,kk-(dir==3)*N3NOT1,NSMASKINSIDE)==1)                   || 
             (GLOBALMACP0A1(nsmask,ii+(dir==1)*N1NOT1                  ,jj+(dir==2)*N2NOT1                  ,kk+(dir==3)*N3NOT1,NSMASKINSIDE)==1)
             ){
            *hasinside=1;
          }
        }
      }
    }
       

    // in cases where no exact on-surface value, indicate that can copy from active value on active grid
    *cancopyfrom=(isinside==0 && *hasmask==1 && *hasinside==1);
    

    //////
    //
    //    So if isinside==0 and hasmask==1 and hasinside==1, then this point is right on (near) the NS surface, so can use that cell to copy from for BCs
    //
    /////






    // GODMARK TODO SUPERMARK: TO FIX FOR CURVATURE VERSION OF ON SURFACE

    // TEST
    if(dir==3){
      if(
         (i==36 && j==23 || i==37 && j==23 || i==38 & j==24 || i==39 && j==25 || i==40 && j==26 || i==40 && j==27)
         || (i==36 && j==41 || i==37 && j==41 || i==38 & j==40 || i==39 && j==39 || i==40 && j==38 || i==40 && j==37)
         ){
        isinside=1;
        *hasmask=1;
        *hasinside=1;
        *reallyonsurface=0;
        *cancopyfrom=0;
      }
    }
    


  }


  return(isinside);
}








// get number of fixedvalues for implied interpolation
int get_fixedvalues(int dir, int i, int j, int k, int ii, int jj, int kk)
{
  int fixedvalues=0;

  // avoid nearest neighbor if it means extrapolation will require passing through more than one fixed point.  That would limit freedom too much and generate kinks for A_i that is present on surface in some cases.
  // Note that one could reduce to lower order if outer-most points are the ones that are fixed, but in general it's always the nearest ones that are fixed since that's the surface.
  // potentially 1-2 (or more?) values for extrapolation may be fixed, so check, and don't consider this as new nearest neighbor if >1 fixed values
  // dir=0 (CENT) would always have only 0 fixed values for concave NS surface.
  if(dir>=1 && dir<=3){

    // get del
    int deloriglocal[NDIM],dellocal[NDIM];
    get_del(i,j,k,ii,jj,kk,deloriglocal,dellocal);

    int isinsidelocal[NDIM],hasmasklocal[NDIM],hasinsidelocal[NDIM],reallyonsurfacelocal[NDIM],cancopyfromlocal[NDIM];
    int fixedvalue[NDIM];

    // get ii,jj,kk result
    isinsidelocal[1]=is_dir_insideNS(dir, ii, jj, kk, &hasmasklocal[1], &hasinsidelocal[1], &reallyonsurfacelocal[1], &cancopyfromlocal[1]);
    fixedvalue[1]=reallyonsurfacelocal[1];  //(isinsidelocal[1]==0 && hasmasklocal[1]==1 && hasinsidelocal[1]==1);

    
    // get iii,jjj,kkk result
    int iii,jjj,kkk;
    iii=ii+dellocal[1];
    jjj=jj+dellocal[2];
    kkk=kk+dellocal[3];
    isinsidelocal[2]=is_dir_insideNS(dir, iii, jjj, kkk, &hasmasklocal[2], &hasinsidelocal[2], &reallyonsurfacelocal[2], &cancopyfromlocal[2]);
    fixedvalue[2]=reallyonsurfacelocal[2];  //=(isinsidelocal[2]==0 && hasmasklocal[2]==1 && hasinsidelocal[2]==1);

    // get iiii,jjjj,kkkk result
    int iiii,jjjj,kkkk;
    iiii=ii+2*dellocal[1];
    jjjj=jj+2*dellocal[2];
    kkkk=kk+2*dellocal[3];
    isinsidelocal[3]=is_dir_insideNS(dir, iiii, jjjj, kkkk, &hasmasklocal[3], &hasinsidelocal[3], &reallyonsurfacelocal[3], &cancopyfromlocal[3]);
    fixedvalue[3]=reallyonsurfacelocal[3];  //=(isinsidelocal[3]==0 && hasmasklocal[3]==1 && hasinsidelocal[3]==1);

    // now have all 3 points.  Is ok if 1 is local, but try to avoid >1
    int fiter;
    SLOOPA(fiter) fixedvalues+=fixedvalue[fiter];

  }
  else fixedvalues=0;


  return(fixedvalues);

}



// check if NS has sufficient grid depth for boundary conditions to be used for interpolations on smoothish functions
// ensure NS is at least 2*N1BND in size in i, 2*N2BND in size in j, 2*N3BND in size in k, so that BCs exist
// Can't have just one point if using more than one point for BCs since then interpolations will see jumps.
// Need to avoid seeing jumps, which requires 2*N?BND in each dimension.
// Near edges of NS, can have less than 2*N?BND since (say) interpolation across (not into) NS will not have big jumps.
// And closer to edge interpolation is, less strong those jumps will be.
int check_nsdepth(SFTYPE time)
{

  int i,j,k;
  int pl,pliter;
  int imask[N1M]={0};
  int jmask[N2M]={0};
  int kmask[N3M]={0};
  int *imaskptr=imask+N1BND;
  int *jmaskptr=jmask+N2BND;
  int *kmaskptr=kmask+N3BND;
  int icount=0,jcount=0,kcount=0;
  int pos=CENT; // CENT forms basis of mask
  FTYPE X[NDIM],V[NDIM];
  int isinside;
  FTYPE ncon[NDIM];
  //
  int ishell[N1M]={0};
  int jshell[N2M]={0};
  int kshell[N3M]={0};
  int *ishellptr=ishell+N1BND;
  int *jshellptr=jshell+N2BND;
  int *kshellptr=kshell+N3BND;
  int isinside2,isshell;
  int icountshell=0,jcountshell=0,kcountshell=0;
  int countinsidezones=0,countshellzones=0;
  FTYPE XX[NDIM],VV[NDIM];
  //
  int ii,jj,kk;
  static int firsttime=1;
  FTYPE nconcon[NDIM];
  FILE *nscheck;





  // setup global pointer for NS masks to be used in boundary condition calls to speed up sub-looping
  if(firsttime==1){
    firsttime=0;
    GLOBALPOINT(nsmask) = (int PTRMACP0A1(nsmask,N1M,N2M,N3M,NUMMASKFLAGS)) (&(BASEMACP0A1(nsmask,N1BND,N2BND,N3BND,0)));
  }

  // always reset mask and shell to zero since only change to 1 if is mask or is shell
  int pf;
  FULLLOOP  for(pf=0;pf<NUMMASKFLAGS;pf++) GLOBALMACP0A1(nsmask,i,j,k,pf) = 0;



  // initialize
  countinsidezones=0;
  // don't include actual or MPI boundaries in count since they don't form real domain of NS that must be resolved
  //  COMPFULLLOOP{
  COMPLOOPN{

    bl_coord_ijk(i,j,k,pos,V);
    isinside=get_insideNS(time,i, j, k, V, NULL,NULL,NULL);
    if(isinside){
      imaskptr[i]=1;
      jmaskptr[j]=1;
      kmaskptr[k]=1;

      GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)=1;
      countinsidezones++;

      for(kk=k-N3NOT1;kk<=k+N3NOT1;kk++){ // TODO: avoid boundary cells?  No need, since setting and extrapolation process can operate inside boundary cells just fine.  Assumes MPI is always last to bound so no machine inconsistencies creep in.
        for(jj=j-N2NOT1;jj<=j+N2NOT1;jj++){
          for(ii=i-N1NOT1;ii<=i+N1NOT1;ii++){

            bl_coord_ijk(ii,jj,kk,pos,VV);
            isinside2=get_insideNS(time,ii, jj, kk, VV, NULL, NULL, NULL);
     

            // if original cell inside, but shifting 1 away is not, then that's the shell bounding the inside region!
            if(isinside2==0){
              ishellptr[ii]=1;
              jshellptr[jj]=1;
              kshellptr[kk]=1;

              GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)=1;
              // can't count "countshellzones" here since will have duplicate offset cells

            }
          } // end over ii
        } // end over jj
      }// end over kk

    } // end if inside


  }// end over COMPLOOPN




  /////////////
  //
  // count linear extents of mask and shell zones
  //
  /////////////
  COMPLOOPN1 if(imaskptr[i]==1) icount++;
  COMPLOOPN2 if(jmaskptr[j]==1) jcount++;
  COMPLOOPN3 if(kmaskptr[k]==1) kcount++;

  COMPLOOPN1 if(ishellptr[i]==1) icountshell++;
  COMPLOOPN2 if(jshellptr[j]==1) jcountshell++;
  COMPLOOPN3 if(kshellptr[k]==1) kcountshell++;


  //////////////
  //
  // Count shell zones
  //
  //////////////
  // don't include actual or MPI boundaries in count since they don't form real domain of NS that must be resolved
  // initialize
  countshellzones=0;
  //  COMPFULLLOOP{
  COMPLOOPN{
    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKSHELL)==1){
      countshellzones++;
    }
  }



  //////////////
  //
  // Collect over multiple processes
  //
  //////////////

#if(USEMPI)
  int isend;
  int jsend;
  int ksend;

  isend=icount;  MPI_Allreduce(&isend, &icount, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);
  jsend=jcount;  MPI_Allreduce(&jsend, &jcount, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);
  ksend=kcount;  MPI_Allreduce(&ksend, &kcount, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);

  isend=icountshell;  MPI_Allreduce(&isend, &icountshell, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);
  jsend=jcountshell;  MPI_Allreduce(&jsend, &jcountshell, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);
  ksend=kcountshell;  MPI_Allreduce(&ksend, &kcountshell, 1, MPI_INT, MPI_SUM, MPI_COMM_GRMHD);
#endif


  //////////////
  //
  // Report results
  //
  //////////////

  trifprintf("NSicount=%d NSjcount=%d NSkcount=%d\n",icount,jcount,kcount);
  trifprintf("NS needs count of: %d %d %d\n",2*N1BND,2*N2BND,3*N3BND);
  trifprintf("NSicountshell=%d NSjcountshell=%d NSkcountshell=%d\n",icountshell,jcountshell,kcountshell);
  trifprintf("NSinsidezones=%d NSshellzones=%d\n",countinsidezones,countshellzones);



  //////////////
  //
  // Find distance to nearest shell grid
  //
  //////////////

  // over i,j,k for which need distance
  int lll;
  int dir;

  int localisinside[NDIM], localhasmask[NDIM], localhasinside[NDIM],localreallyonsurface[NDIM],localcancopyfrom[NDIM];
  int hasmask,hasinside,reallyonsurface,cancopyfrom;


  //////////////////////////
  //
  // main loop over i,j,k
  //
  //////////////////////////
  COMPFULLLOOP{
    

    // initialize nsmask
    // only relevant for cells that act as boundary cells (i.e. inside NS in relevant sense)
    // CENT and A_{dir} @ CORN_{dir}
    for(dir=0;dir<=3;dir++){
      localisinside[dir]=is_dir_insideNS(dir, i, j, k, &localhasmask[dir], &localhasinside[dir], &localreallyonsurface[dir], &localcancopyfrom[dir]);

      // set up default
      GLOBALMACP0A1(nsmask,i,j,k,NSMASKDISTTOSHELL+dir)=0;
      GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 0)=0;
      GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 1)=0;
      GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 2)=0;

    }




    // if at least inside for some cell position, can enter ii,jj,kk loop
    if(localisinside[0]==1 || localisinside[1]==1 || localisinside[2]==1 || localisinside[3]==1){


#define MINDISTLOCAL 1.01
#define MAXNEIGH 10 // maximum number of last few neighbors (so many because tracking DUP's as well

      FTYPE dist; // for CENT,CORN1,CORN2,CORN3
      // start assuming large and find minimum
      FTYPE mindist[NDIM][MAXNEIGH]={{BIG}}; // [which dir][which neighbor]=[this neighbor's distance, corresponding to some smallest(smallish) distance relative to other neighbors]

      int lastneigh[NDIM][MAXNEIGH][NDIM]={{{-100}}}; //[which dir][which last neighbor][which i,j,k = 1,2,3]=[ii or jj or kk]
      int lastfixed[NDIM][MAXNEIGH]={{NBIGM}}; // [which dir][which neighbor]=[number of fixed points for implied 3-point interpolation (prfulldel)]
      int lastcount[NDIM][MAXNEIGH]={{0}}; // [which dir][which neighbor]=[number of DUPs so far for this entry]
      int numneigh[NDIM]={0}; // [which dir]=[number of neighbors]
      int count[NDIM][MAXNEIGH]={{0}}; // [which dir][which neighbor]=[number of neighbors that have same dist (duplicates)]
      int neighbase=0; // neighbase is fixed to be 0 as per below algorithm's method of shifting
      int newfixed[NDIM];

      // GODMARK: Why do above nested curly braces not set initial assignment?
      // Seems I have to do it myself for mindist at least
      for(dir=0;dir<=3;dir++){
        numneigh[dir]=0;
        int neighiter;
        for(neighiter=0;neighiter<MAXNEIGH;neighiter++){
          mindist[dir][neighiter]=BIG;
          count[dir][neighiter]=0;

          int spaceiter;
          DLOOPA(spaceiter) lastneigh[dir][neighiter][spaceiter]=-100;
          lastfixed[dir][neighiter]=NBIGM;
          lastcount[dir][neighiter]=0;
        }
      }

      // http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp?topic=%2Fcom.ibm.xlcpp8a.doc%2Flanguage%2Fref%2Faryin.htm
      // static int number[3] = { [0] = 5, [2] = 7 }; // skip over elements in init of array

      ///////////////////
      //
      // Get neighbor(s) (loop over ii,jj,kk)
      //
      ///////////////////
      COMPLOOPNiijjkk(ii,jj,kk){
 
        // CENT and CORN_dir (GODMARK: TODO: This code has become very expensive in bound if doing 3D, so maybe ignore maxdist and just use whatever find?)
        for(dir=0;dir<=3;dir++){ // loop over dir's

          if(localisinside[dir]==1){ // only do this if inside NS in relevant sense
            isinside=is_dir_insideNS(dir, ii, jj, kk, &hasmask, &hasinside, &reallyonsurface, &cancopyfrom);

    
            //     if(isinside==0 && hasmask==1 && hasinside==1){ // if get inside here, then ii,jj,kk for dir should be cell that can be copied FROM and is as close as possible to surface for such a source of data
            if(cancopyfrom){ // if get inside here, then ii,jj,kk for dir should be cell that can be copied FROM and is as close as possible to surface for such a source of data

              // DEBUG:
              //       dualfprintf(fail_file,"WTF: dir=%d ijk=%d %d %d : ijk2=%d %d %d : mindist=%21.15g\n",dir,i,j,k,ii,jj,kk,mindist[dir][0]);


              // get dist for this ii,jj,kk
              dist=sqrt( (ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k) );

              // get fixed values count
              newfixed[dir]=get_fixedvalues(dir, i, j, k, ii, jj, kk);

              // DEBUG:
              //       if(dir==3 && i==36 && j==24){
              //  dualfprintf(fail_file,"FU: dir=%d ijk=%d %d %d : ijk2=%d %d %d : dist=%21.15g : numneigh=%d mindist=%21.15g\n",dir,i,j,k,ii,jj,kk,dist,numneigh[dir],mindist[dir][0]);
              //       }


              int founddup=0;
              // first loop over any existing neighbors and see if dist is similar to their dist that was chosen as some prior mindist
              {
                int neighiter;
                for(neighiter=0;neighiter<MIN(MAXNEIGH,numneigh[dir]);neighiter++){ // only up to numneigh[dir] (unless larger than MAXNEIGH) since DUPs only exist if array filled-in already

                  if(newfixed[dir]==lastfixed[dir][neighiter] && fabs(dist-mindist[dir][neighiter])<1E-3){ // if within same cell distance and fixed interpolation count, then treat as duplicate since can't distinguish
                    // DEBUG:
                    dualfprintf(fail_file,"DUP: count=%d : dir=%d ijk=%d %d %d : ijk2=%d %d %d : dist=%21.15g fixed=%d : lorig123=%d %d %d\n",count[dir][neighiter],dir,i,j,k,ii,jj,kk,dist,newfixed[dir],i+lastneigh[dir][neighiter][1],j+lastneigh[dir][neighiter][2],k+lastneigh[dir][neighiter][3]);

                    // then just keep adding to same neighbor
                    count[dir][neighiter]++;
                    // accumulate ii,jj,kk as vector sum
                    lastneigh[dir][neighiter][1]+=(ii-i);
                    lastneigh[dir][neighiter][2]+=(jj-j);
                    lastneigh[dir][neighiter][3]+=(kk-k);
                    // add up number of fixed points
                    //      lastfixed[dir][neighiter]+=newfixed[dir]; // note we are summing number of fixed, so in the end question is to obtain lowest total or equally number of fixed per duplicate
                    lastfixed[dir][neighiter]=newfixed[dir]; // same fixed, but ok to assign again
                    mindist[dir][neighiter]=dist; // same if dup, but ok to assign again

                    founddup=1; // indicate that found duplicate

                    // DEBUG:
                    dualfprintf(fail_file,"DUP2: count=%d : dir=%d ijk=%d %d %d : ijk2=%d %d %d : dist=%21.15g fixed=%d : l123=%d %d %d\n",count[dir][neighiter],dir,i,j,k,ii,jj,kk,dist,newfixed[dir],i+lastneigh[dir][neighiter][1],j+lastneigh[dir][neighiter][2],k+lastneigh[dir][neighiter][3]);


                  }
                }
              }

              int gotmin=0;
              int whichmin=0;
              if(founddup==0){
                {
                  // second, loop over and see if really closer or fewer fixed points than any others
                  int neighiter;

                  // then really smaller distance (used to avoid machine precision issues) or a lower fixed count
                  //    if(dist<mindist[dir][neighiter]-1E-3){ // then really smaller distance (used to avoid machine precision issues)
                  // (newfixed[dir]==lastfixed[dir][neighiter] && dist<mindist[dir][neighiter]-1E-3)
                  // (dist<mindist[dir][neighiter]-1E-3) || (newfixed[dir]<=lastfixed[dir][neighiter] && fabs(dist-mindist[dir][neighiter])<1E-3)


                  // gotmin=0: same or smaller distance and same or fewer fixed points (most preferred to trigger on)
                  // this should be good since if ever both were same (so undetermined duplicate), above duplicate check should merge the points already
                  for(neighiter=0;neighiter<MAXNEIGH;neighiter++){ // not up to numneigh[dir] because it's 0 at first -- and in general want to trigger on being smaller than existing large values if filling-in new point
                    if( newfixed[dir]<=lastfixed[dir][neighiter] && (dist<mindist[dir][neighiter]-1E-3 || fabs(dist-mindist[dir][neighiter])<1E-3) ){
                      gotmin=1;
                      whichmin=neighiter;
                      break; // will use whichmin to insert this data into storage arrays
                    }
                  }


                  if(gotmin==0){
                    dualfprintf(fail_file,"Still gotmin==0 (This is ok, but check that really wanted to ignore this point): dir=%d ijk=%d %d %d iijjkk=%d %d %d : dist=%21.15g\n",dir,i,j,k,ii,jj,kk,dist);
                  }

                }
              }

              // DEBUG:
              //       if(dir==3 && i==36 && j==24){
              //  dualfprintf(fail_file,"FU2: dir=%d ijk=%d %d %d : ijk2=%d %d %d : gotmin=%d (%21.15g %21.15g)\n",dir,i,j,k,ii,jj,kk,gotmin,dist,mindist[dir][whichmin]);
              //       }


              // see if min over all existing stored neighbors
              //if(gotmin>=MIN(MAXNEIGH,numneigh[dir])){ // note that this gets triggered if never any other neighbors yet, as required.
              if(gotmin!=0){ // note that this gets triggered if never any other neighbors yet, as required.
                // DEBUG:
                dualfprintf(fail_file,"NOTDUP : gotmin=%d count=%d : dir=%d ijk=%d %d %d : ijk2=%d %d %d : dist: %21.15g mindist: %21.15g fixed=%d\n",gotmin,count[dir][neighbase],dir,i,j,k,ii,jj,kk,dist,mindist[dir][neighbase],newfixed[dir]);
                // then got closer neighbor, so shift neighbor list and assign
  
                // shift neighbors (lose neighbors if there are more than MAXNEIGH neighbors, which is fine since don't expect too many equally valid neighbors)
                int neighiter;
                for(neighiter=MAXNEIGH-1;neighiter>=whichmin+1;neighiter--){// only down to whichmin+1, since [whichmin] is old point for which new point is just closer. So [whichmin] is moved up to [whichmin+1] to make room for new point at [whichmin].  This sorts by distance, else could lose good small distance out the right side of array!
                  count[dir][neighiter]=count[dir][neighiter-1];
                  int siter; SLOOPA(siter) lastneigh[dir][neighiter][siter]=lastneigh[dir][neighiter-1][siter];
                  lastfixed[dir][neighiter]=lastfixed[dir][neighiter-1];
                  mindist[dir][neighiter]=mindist[dir][neighiter-1];
                }
  
                // now assign values for (not add values to) new neighbor
                numneigh[dir]++; // really new neighbor at new distance
                // new values for new neighbor
                count[dir][whichmin]=1; // sets as not duplicate count
                lastneigh[dir][whichmin][1]=(ii-i);
                lastneigh[dir][whichmin][2]=(jj-j);
                lastneigh[dir][whichmin][3]=(kk-k);
                lastfixed[dir][whichmin]=newfixed[dir];
                mindist[dir][whichmin]=dist;

                dualfprintf(fail_file,"NOTDUP2 : whichmin=%d count=%d : dir=%d ijk=%d %d %d : ijk2=%d %d %d : dist: %21.15g mindist: %21.15g fixed=%d l123=%d %d %d\n",whichmin,count[dir][whichmin],dir,i,j,k,ii,jj,kk,dist,mindist[dir][whichmin],newfixed[dir],i+lastneigh[dir][whichmin][1],j+lastneigh[dir][whichmin][2],k+lastneigh[dir][whichmin][3]);

              }// end if really got new min over all prior neighbors

            }// end if not inside (and if A_i then has mask)
          }
        }// over dir
      
      }// over ii,jj,kk




      
#define DOSTUCKFIX 1

      // check if i,j,k==ii,jj,kk and shift over such that removes such a useless point
      if(DOSTUCKFIX){
        for(dir=0;dir<=3;dir++){

          numneigh[dir]=MIN(numneigh[dir],MAXNEIGH); // set number of neighbors as real possible upper limit since may reduce below MAXNEIGH now
          for(int neighiter=0;neighiter<MIN(numneigh[dir],MAXNEIGH);neighiter++){

            if(lastneigh[dir][neighiter][1]==0 && lastneigh[dir][neighiter][2]==0 && lastneigh[dir][neighiter][3]==0){
              // then shift on top of this useless point (even if formed from duplicates, since couldn't choose among them unless one would break symmetry (or find another indicator))

              for(int neighiter2=neighiter;neighiter2<MIN(numneigh[dir],MAXNEIGH)-1;neighiter2++){ // -1 since lost data from removal
                count[dir][neighiter2]=count[dir][neighiter2+1];
                int siter; SLOOPA(siter) lastneigh[dir][neighiter2][siter]=lastneigh[dir][neighiter2+1][siter];
                lastfixed[dir][neighiter2]=lastfixed[dir][neighiter2+1];
                mindist[dir][neighiter2]=mindist[dir][neighiter2+1];
              }


              // reduce count (will go below MAXNEIGH now so don't show (in debug) the upper-end that was used to copy to lower end)
              numneigh[dir]--;
              // and need to drop back one, since will next be doing neighiter++ and want to evaluate the same neighiter as here because it's actually new data not checked yet
              neighiter--;

              // DEBUG:
              dualfprintf(fail_file,"DOSTUCKFIX: dir=%d ijk=%d %d %d : (new)neighiter=%d (new)numneigh=%d\n",dir,i,j,k,neighiter,numneigh[dir]);


            }// end if stuck point
          }// end over neighiter
        }// end over dir
      }// end if fixing stuck points








#define DONEIGHBORFIX 1


      // TODO: If not last, but second to last, was part of a DUP group, then choosing just the last one of the DUP's will produce asymmetries.  But other DUP members may be later than once mindist triggered on non-DUP.
      // Maybe re-find all similar dist neighbors?  Expensive?
      // Or, keep track of not just mindist but all prior stored neighbors?  probably good.  Problem: DUP versions need to be fixed below, which comes after
      
      // choose neighbor that has lowest number of fixed values used for extrapolation.  This gives extrapolation more freedom/flexibility for field moving on surface to avoid kinks through the surface.
      if(DONEIGHBORFIX){
        for(dir=0;dir<=3;dir++){ // only A_i (dir=1,2,3) could have multiple fixed points, but no harm in keeping this general

          // count[dir][neighbase] contains final count of dups for final ii,jj,kk
          //   if(count[dir][neighbase]==1){ // GODMARK TODO OPTMARK: for now only do for count==1 since count>1 has to be modified to be a closer neighbor (DODUPFIX) and count>1 (expensively) handled in interpolation code so already can reduce to lower fixedvalues.  If also dealt with count>1, would make interpolations less expensive

          // loop over neighbors in order from latest to oldest:  in general because last ones would have lower dist[] values and so closer as desirable.
          int neighiter;
          //   FTYPE lowestfixed=(FTYPE)lastfixed[dir][neighbase]/((FTYPE)count[dir][neighbase]); // start with existing neighbor's fixed point count per duplicate
          FTYPE lowestfixed=(FTYPE)lastfixed[dir][neighbase]; // now grouping by fixed counts too and not adding-up fixed counts for dups, so only look at fixed count itself
          for(neighiter=1;neighiter<MIN(numneigh[dir],MAXNEIGH);neighiter++){ // only iterate for existing number of neighbors actually stored up to maximum number stored.
            //     FTYPE otherfixed=(FTYPE)lastfixed[dir][neighiter]/((FTYPE)(count[dir][neighiter]));
            FTYPE otherfixed=(FTYPE)lastfixed[dir][neighiter];
            if(otherfixed<lowestfixed-1E-3){ // only choose if strictly lower, not just equal, since if only equal then might as well use closer neighbor.  1E-3 avoids machine precision issues so really must be actually lower, not just lower by machine error

              // DEBUG:
              dualfprintf(fail_file,"CHOSENEWNEIGHBOR: dir=%d ijk=%d %d %d oldijk2=%d %d %d newijk2=%d %d %d : oldfixed=%21.15g otherfixed=%21.15g : oldcount=%d newcount=%d : oldmindist=%21.15g newmindist=%21.15g\n",dir,i,j,k,lastneigh[dir][neighbase][1],lastneigh[dir][neighbase][2],lastneigh[dir][neighbase][3],lastneigh[dir][neighiter][1],lastneigh[dir][neighiter][2],lastneigh[dir][neighiter][3],lowestfixed,otherfixed,count[dir][neighbase],count[dir][neighiter],mindist[dir][neighbase],mindist[dir][neighiter]);

              lowestfixed=otherfixed;
              // overwrite [neighbase] with this version
              count[dir][neighbase]=count[dir][neighiter];
              lastneigh[dir][neighbase][1]=lastneigh[dir][neighiter][1];
              lastneigh[dir][neighbase][2]=lastneigh[dir][neighiter][2];
              lastneigh[dir][neighbase][3]=lastneigh[dir][neighiter][3];
              lastfixed[dir][neighbase]=lastfixed[dir][neighiter];
              mindist[dir][neighbase]=mindist[dir][neighiter];
       
            }// end if strictly lower fixed count
          }// end over neighiter
          //}// end if final ii,jj,kk was not a duplicate
        }// end over dir

      }// end if fixing neighbors



      
      // DEBUG:
      for(dir=0;dir<=3;dir++){ // only A_i (dir=1,2,3) could have multiple fixed points, but no harm in keeping this general
        for(int neighiter=0;neighiter<MIN(numneigh[dir],MAXNEIGH);neighiter++){
          dualfprintf(fail_file,"REPORTBEFOREDUPFIX: dir=%d ijk=%d %d %d : neighiter=%d :: count=%d initial iijjkk=%d %d %d fixed=%d mindist=%21.15g\n",dir,i,j,k,neighiter,count[dir][neighiter],i+lastneigh[dir][neighiter][1],j+lastneigh[dir][neighiter][2],k+lastneigh[dir][neighiter][3],lastfixed[dir][neighiter],mindist[dir][neighiter]);
        }
      }




#define DODUPFIX 1

      // NOTEMARK: When DUPFIX ends up using corner cell(s) instead of shorter direct cells, then on substeps the copy will be less up-to-date than the copies from face directions since flux update is along face directions and not corners.  Should be updated within another substep.
      // This causes, e.g., pr[U3] to be zero in ghost cells because that ghost copies from corner to itself that hasn't been updated in a single substep once rotation is turned on.


      if(DODUPFIX){
        ///////////////////
        //
        // divide out counts for multiple same-distance neighbors if exist
        //
        ///////////////////
        int interpproblem;
        int delsuperorig[4],delorig[4],del[4];
        int cancopyfromold;
        int whilecheckiter;
        int isongrid;
        int trialshift;
        for(dir=0;dir<=3;dir++){
          if(localisinside[dir]==1){ // only do this if inside NS in relevant sense

            if(count[dir][neighbase]==0){
              dualfprintf(fail_file,"Never found nearest neighbor! dir=%d : %d %d %d\n",dir,i,j,k);
            }
            else if(count[dir][neighbase]==1){
              // then nothing else to do except worry that too many fixed values may be used in the interpolation if along grid lines when dealing with A_i that can be on surface of NS)
            }
            else if(count[dir][neighbase]>1){
              // if multiple cases, then iterate back towards i,j,k until find more closest neighbor along that path.


              // set initial ii,jj,kk
              ii=i+lastneigh[dir][neighbase][1];
              jj=j+lastneigh[dir][neighbase][2];
              kk=k+lastneigh[dir][neighbase][3];


              // DEBUG:
              dualfprintf(fail_file,"COUNT: dir=%d ijk=%d %d %d : count=%d initial iijjkk=%d %d %d\n",dir,i,j,k,count[dir][neighbase],ii,jj,kk);


              if(i==ii && j==jj & k==kk){
                dualfprintf(fail_file,"DUP got ijk==ijk2: %d %d %d : %d %d %d\n",i,j,k,ii,jj,kk);
                myexit(18492352);
              }

              // check if already a mask cell
              isongrid=(ii>=-N1BND && ii<=N1M-1 &&  jj>=-N2BND && jj<=N2M-1 &&  kk>=-N3BND && kk<=N3M-1);
              if(isongrid==0){ isinside=hasmask=hasinside=reallyonsurface=cancopyfrom=0; }
              else isinside=is_dir_insideNS(dir, ii, jj, kk, &hasmask, &hasinside, &reallyonsurface, &cancopyfrom);


              //       if(! (hasmask==1 && hasinside==1 )){ // only need to do something if not already on boundary (boundary that is values that can copy from)
              if(cancopyfrom==0){ // only need to do something if not already on boundary (boundary that is values that can copy from)

                // DEBUG:
                dualfprintf(fail_file,"INITIALDUPFIX: dir=%d count=%d : %d %d %d : %d %d %d : ihhrc=%d %d %d %d %d\n",dir,count[dir][neighbase],i,j,k,ii,jj,kk, isinside, hasmask, hasinside, reallyonsurface, cancopyfrom);

                // get delsuperorig
                interpproblem=get_del(i,j,k,ii,jj,kk,delsuperorig,del);


                whilecheckiter=0;
                trialshift=0;
                int ii0,jj0,kk0;
                int ii0true,jj0true,kk0true;
                while(1){
                  // get lowest order offset for each ii,jj,kk that is tried
                  interpproblem=get_del(i,j,k,ii,jj,kk,delorig,del);

                  // original ii,jj,kk before changed for test, so is always outside NS
                  ii0=ii;
                  jj0=jj;
                  kk0=kk;

                  // in some rare cases, may still be inside NS even if dup, so then need to move out instead of in.
                  int signdir=(isinside==1 ? 1 : -1);

                  if(trialshift==0){

                    // original using actual path
                    ii0true=ii;
                    jj0true=jj;
                    kk0true=kk;

                    // shift ii,jj,kk back to i,j,k
                    ii=ii+signdir*del[1];
                    jj=jj+signdir*del[2];
                    kk=kk+signdir*del[3];
                  }
                  else if(trialshift>0){
                    // shift by direction with largest absolute offset and then only by 1 unit each step in that direction, unless that direction has no offset to begin with
                    ii=ii+signdir*ROUND2INT(sign(del[1])) * MIN(1,abs(del[1])) * (abs(del[1])>=abs(del[2]) && abs(del[1])>=abs(del[3]));
                    jj=jj+signdir*ROUND2INT(sign(del[2])) * MIN(1,abs(del[2])) * (abs(del[2])>=abs(del[1]) && abs(del[2])>=abs(del[3]));
                    kk=kk+signdir*ROUND2INT(sign(del[3])) * MIN(1,abs(del[3])) * (abs(del[3])>=abs(del[1]) && abs(del[3])>=abs(del[2]));

                    trialshift++;

                    // DEBUG:
                    dualfprintf(fail_file,"2Doing trialshift: dir=%d : %d %d %d : %d %d %d\n",dir,i,j,k,ii,jj,kk);
                  }


                  // check if external to NS (or on shell)
                  // only check if on the grid
                  isongrid=(ii>=-N1BND && ii<=N1M-1 &&  jj>=-N2BND && jj<=N2M-1 &&  kk>=-N3BND && kk<=N3M-1);
                  if(isongrid==0){ isinside=hasmask=hasinside=reallyonsurface=cancopyfrom=0; }
                  else isinside=is_dir_insideNS(dir, ii, jj, kk, &hasmask, &hasinside, &reallyonsurface, &cancopyfrom);

                  // DEBUG:
                  dualfprintf(fail_file,"DUPFIXING: dir=%d count=%d : %d %d %d : %d %d %d : %d %d %d %d %d\n",dir,count[dir][neighbase],i,j,k,ii,jj,kk,isinside,hasmask,hasinside,reallyonsurface,&cancopyfrom);
                  int delloop;
                  SLOOPA(delloop) dualfprintf(fail_file,"DUPFIXING2: del=%d : %d %d %d\n",delloop,del[delloop],delorig[delloop],delsuperorig[delloop]);


                  if(interpproblem==0 && isongrid==1 && cancopyfrom==1){
                    // then this ii,jj,kk the one we want to keep
                    break;
                  }

                  if(interpproblem==0 && isongrid==1 && isinside==1){// then inside (went too far)
  
                    // then went too far, so back to original before inside
                    ii=ii0;
                    jj=jj0;
                    kk=kk0;

                    if(cancopyfromold==1){
                      // then done!  back-up cell is a mask cell, so that's good
                      break;
                    }
                    else{
                      // try shifting by 1 unit (trying to get a mask cell)

                      // DEBUG:
                      dualfprintf(fail_file,"Doing trialshift: dir=%d : ijk=%d %d %d : ts=%d : ijk2=%d %d %d\n",dir,i,j,k,trialshift,ii,jj,kk);

                      trialshift++;
                    }
                  }// end if inside NS


                  // check how trial shift is going
                  if(trialshift>MAX(3,MAX( MAX( abs(delsuperorig[1]),abs(delsuperorig[2]) ) , abs(delsuperorig[3]) )  ) ){
                    // then never found mask cell when back shifting, so revert to latest cell that's outside NS
                    // can't trust odd path, so revert to true original
                    ii=ii0true;
                    jj=jj0true;
                    kk=kk0true;

                    dualfprintf(fail_file,"Got duplicate, found external that's close, but it's not a mask value: dir=%d : ijk=%d %d %d : ts=%d : ijk2=%d %d %d\n",dir,i,j,k,trialshift,ii,jj,kk);
                    dualfprintf(fail_file,"count [dir=%d]=%d\n",dir,count[dir]);
                    int testtest;
                    for(testtest=1;testtest<=3;testtest++) dualfprintf(fail_file,"dir=%d delorig=%d del=%d\n",testtest,delorig[testtest],del[testtest]);

                    break;
                  }

                  // if here, then still going
                  cancopyfromold=cancopyfrom;

                  if(whilecheckiter<100) whilecheckiter++;
                  else{
                    dualfprintf(fail_file,"whilecheckiter=%d caught at dir=%d : %d %d %d : %d %d %d\n",whilecheckiter,dir,i,j,k,ii,jj,kk);
                  }
                }// end while loop

     
                // DEBUG:
                dualfprintf(fail_file,"DUPDONE: dir=%d count=%d : %d %d %d : %d %d %d\n",dir,count[dir][neighbase],i,j,k,ii,jj,kk);


                // get final modified ii,jj,kk
                lastneigh[dir][neighbase][1]=(ii-i);
                lastneigh[dir][neighbase][2]=(jj-j);
                lastneigh[dir][neighbase][3]=(kk-k);
  
              }// end if original ii,jj,kk is not a mask cell already
            }// end if multiple counts for this dir
          }// end if inside NS in relevant sense
     
        }// over all dirs=0,1,2,3
 
      }// end if DUPFIX

      
      
      /////////////////
      //
      // Make final assignments to nsmask array
      //
      ////////////////
      for(dir=0;dir<=3;dir++){
        if(localisinside[dir]==1){ // only do this if inside NS in relevant sense
   
          // store ii,jj,kk that is nearest neighbor to i,j,k 
          // uses [neighbase] as reference
          GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 0)=lastneigh[dir][neighbase][1];
          GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 1)=lastneigh[dir][neighbase][2];
          GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 2)=lastneigh[dir][neighbase][3];
  
          // force lowest min to be 1.01
          // assume never itself (i.e. 0 distance) or would use that fact
          if(mindist[dir][neighbase]<MINDISTLOCAL){
            mindist[dir][neighbase]=MINDISTLOCAL; // if i,j,k is not on shell, and found shell, then must be 1 cell away, so avoid (int)(0.9999)->0 and set as 1.01 so (int)1.01 -> 1
          }

          // now have mindist for this i,j,k
          GLOBALMACP0A1(nsmask,i,j,k,NSMASKDISTTOSHELL+dir)=(int)mindist[dir][neighbase];

        }// end if inside NS in relevant sense
      }// over all dirs=0,1,2,3

    }// end if cell is inside NS in sense needed for BCs
  }// over each i,j,k

  




  //////////////////////////
  //
  // report/debug
  //
  //////////////////////////
  if(numprocs==1){

    nscheck=fopen("nscheck.dat","wt");
    if(nscheck==NULL){
      dualfprintf(fail_file,"Couldn't open nscheck.dat\n");
      myexit(350683463);
    }


    int maski;
    COMPLOOPN{
      fprintf(nscheck,"%d %d %d",i,j,k);
      for(maski=0;maski<NUMMASKFLAGS;maski++){
        fprintf(nscheck," %d",GLOBALMACP0A1(nsmask,i,j,k,maski));
      }
      fprintf(nscheck,"\n");
    }

    fclose(nscheck);


    if(N3==1){ // for debug only


      if(1){
        // just plot all types separately and merge manually by shifting graphic or paper
        // in Konsole shell, can use smaller font, print to file as .ps (all 4 dirplots), open in photoshop, overlay all 4 appropriately, save, flatten, copy section, paste into new image, set size as 8.5x?" size, and print.
        int dirplot,dira;
        for(dirplot=0;dirplot<=4;dirplot++){
          dualfprintf(fail_file,"dirplot=%d\n",dirplot);
          k=0;
          COMPLOOPN2{
            COMPLOOPN1{
              int isinsideplot[NDIM],hasmaskplot[NDIM],hasinsideplot[NDIM],reallyonsurfaceplot[NDIM],cancopyfromplot[NDIM];
              int isonsurfaceplot[NDIM];
    
              if(dirplot<4) dira=dirplot;
              else if(dirplot==4) dira=0;
 
              isinsideplot[dira]=is_dir_insideNS(dira, i, j, k, &hasmaskplot[dira], &hasinsideplot[dira], &reallyonsurfaceplot[dira], &cancopyfromplot[dira]);
              isonsurfaceplot[dira]=cancopyfromplot[dira]; //(isinsideplot[dira]==0 && hasmaskplot[dira]==1 && hasinsideplot[dira]==1);
   
              if(dira==0){ // CENT
                if(dirplot==0){
                  if(isinsideplot[dira]) dualfprintf(fail_file,"I");
                  else  dualfprintf(fail_file,"O");
                }
                else if(dirplot==4){
                  if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKSHELL)) dualfprintf(fail_file,"S");
                  else  dualfprintf(fail_file," ");
                }  
              }
              else if(dira==1){ // CORN1==FACE2
                if(isonsurfaceplot[dira]) dualfprintf(fail_file,"-");
                else  dualfprintf(fail_file," ");
              }
              else if(dira==2){ // CORN2==FACE1
                if(isonsurfaceplot[dira]) dualfprintf(fail_file,"|");
                else  dualfprintf(fail_file," ");
              }
              else if(dira==3){ // CORN3
                if(isonsurfaceplot[dira]) dualfprintf(fail_file,".");
                else  dualfprintf(fail_file," ");
              }
            }// end COMPLOOPN1
            dualfprintf(fail_file,"\n");
          }// end COMPLOOPN2
        }// end over dirs
      }


    }// end if N3==1


  }//end if numprocs==1


  //////////////
  //
  // Fail if not enough zones to deal with boundary conditions in each dimension
  //
  //////////////

  if(icount<2*N1BND || jcount<2*N2BND || kcount<2*N3BND){
    myexit(2498346112);
  }
  



  return(0);

}


// Take B^i and get B_\phi
int rescale_Bcon_to_Bcovphi(FTYPE *pr,struct of_geom *ptrrgeom, FTYPE *Bd3)
{
  int jj;

  /////////////
  //
  // obtain reference B_\phi (not true B_\phi = *F_{\phi t} but close enough for NS where g_{ti} is small)
  //
  /////////////
  *Bd3=0.0;
  SLOOPA(jj) *Bd3 += pr[B1+jj-1]*(ptrrgeom->gcov[GIND(3,jj)]);

  return(0);

}


// Take B^1 B^2 and B_3 and get back B^3
int rescale_Bconp_and_Bcovphi_to_Bconphi(FTYPE *pr, struct of_geom *ptrgeom, FTYPE Bd3)
{
  FTYPE Bu1,Bu2,gcon03,gcon13,gcon23,gcon33;
  FTYPE gcov01,gcov02,gcov11,gcov12,gcov21,gcov22,gcov03,gcov13,gcov23;
  FTYPE myBd3;
  FTYPE ftemp,iftempnosing;
  int pl;


  Bu1=pr[B1];
  Bu2=pr[B2];

  gcon03=ptrgeom->gcon[GIND(0,3)];
  gcon13=ptrgeom->gcon[GIND(1,3)];
  gcon23=ptrgeom->gcon[GIND(2,3)];
  gcon33=ptrgeom->gcon[GIND(3,3)];

  gcov01=ptrgeom->gcov[GIND(0,1)];
  gcov02=ptrgeom->gcov[GIND(0,2)];
  gcov11=ptrgeom->gcov[GIND(1,1)];
  gcov12=gcov21=ptrgeom->gcov[GIND(1,2)];
  gcov22=ptrgeom->gcov[GIND(2,2)];
  gcov03=ptrgeom->gcov[GIND(0,3)];
  gcov13=ptrgeom->gcov[GIND(1,3)];
  gcov23=ptrgeom->gcov[GIND(2,3)];
   
  // Bd3fromBu.nb (just moved signs)
  myBd3=Bd3; // + dqBd3*(i-ri);
  ftemp=(1.0 - gcon03*gcov03 - gcon13*gcov13 - gcon23*gcov23);
  iftempnosing=sign(ftemp)/(fabs(ftemp)+SMALL);
  pl=B3; pr[pl] = (myBd3*gcon33 + Bu1*gcon03*gcov01 + Bu2*gcon03*gcov02 + Bu1*gcon13*gcov11 + Bu2*gcon13*gcov12 + Bu1*gcon23*gcov21 + Bu2*gcon23*gcov22)*iftempnosing;

  // old B_\phi imprecise copy (need to avoid singularity anywways)
  //   pl=B3; MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl)*fabs((ptrrgeom[pl]->gcov[GIND(3,3)])/(ptrgeom[pl]->gcov[GIND(3,3)]));

  return(0);

}




// take prin[pl] -> prout[all same except pl=B3 has Bd3 in it]
int rescale_pl(SFTYPE time, int dir, struct of_geom *ptrgeom, FTYPE *prin, FTYPE *prout)
{
  int pliter,pl;
  int rescale_Bcon_to_Bcovphi(FTYPE *pr,struct of_geom *ptrrgeom, FTYPE *Bd3);

  if(dir==0){ // then CENT prim
    // default is all same in and out
    PLOOP(pliter,pl) prout[pl]=prin[pl];
    
    
    // get B_\phi and use it as thing to interpolate
#if(1) // GODMARK TODO DEBUG
    rescale_Bcon_to_Bcovphi(prin,ptrgeom, &prout[B3]);
#endif
    // so prout[B3] has Bd3 in it
  }
  else{ // then A_{dir}
    
    FTYPE rescale;
    rescale=rescale_A(time, dir,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    prout[0]=prin[0]*rescale;

  }


  return(0);

}



// take prin[all same except pl=B3 has Bd3 in it] -> prout[normal pl list with pl=B3 having B^3 in it]
int unrescale_pl(SFTYPE time, int dir, struct of_geom *ptrgeom, FTYPE *prin, FTYPE *prout)
{
  int pliter,pl;
  int rescale_Bconp_and_Bcovphi_to_Bconphi(FTYPE *pr, struct of_geom *ptrgeom, FTYPE Bd3);

  if(dir==0){ // then CENT prim
    
    // default is all same in and out
    PLOOP(pliter,pl) prout[pl]=prin[pl];
    
    
    // puts Bconphi into prim[B3]
    // so prin has B^1 and B^2 and prin[B3] contains Bd3
#if(1) // GODMARK TODO DEBUG
    rescale_Bconp_and_Bcovphi_to_Bconphi(prout, ptrgeom, prin[B3]);
#endif
    // so prout[B3] will be filled with B^3 as desired
  }
  else{ // then A_{dir}
    
    FTYPE rescale;
    rescale=rescale_A(time,dir,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    prout[0]=prin[0]/rescale;

  }
  

  return(0);

}






// check for monotonicity of using 4 point *interpolation*
int checkmono4(int firstpos, int lastpos, FTYPE *xpos, FTYPE y0, FTYPE y1, FTYPE y2, FTYPE y3, FTYPE y4)
{
  FTYPE x0=xpos[0]; // x0 in middle (y0 not set, so not used)
  FTYPE x1=xpos[1]; // x1 on edge
  FTYPE x2=xpos[2];
  FTYPE x3=xpos[3];
  FTYPE x4=xpos[4]; // x4 on other edge

  FTYPE xextreme[2];

  int mono;

  // /home/jon/testing_nsbh.nb

  xextreme[0]=0.16666666666666666*pow(x1*(x1 - 1.*x4)*x4*(y2 - 1.*y3) + 
                                      (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,2.) + 
                                      (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,2.) + 
                                      x3*((-1.*y2 + y4)*pow(x1,2.) + (-1.*y1 + y2)*pow(x4,2.)) + 
                                      x2*((y3 - 1.*y4)*pow(x1,2.) + (-1.*y1 + y4)*pow(x3,2.) + (y1 - 1.*y3)*pow(x4,2.)),
                                      -1.)*(2.*(x1*(x1 - 1.*x4)*x4*(x1 + x4)*(y2 - 1.*y3) + 
                                                (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,3.) + 
                                                (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,3.) + 
                                                x3*((-1.*y2 + y4)*pow(x1,3.) + (-1.*y1 + y2)*pow(x4,3.)) + 
                                                x2*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                    (y1 - 1.*y3)*pow(x4,3.))) - 
                                            1.*pow(-12.*(x1*(x1 - 1.*x4)*x4*(y2 - 1.*y3) + 
                                                         (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,2.) + 
                                                         (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,2.) + 
                                                         x3*((-1.*y2 + y4)*pow(x1,2.) + (-1.*y1 + y2)*pow(x4,2.)) + 
                                                         x2*((y3 - 1.*y4)*pow(x1,2.) + (-1.*y1 + y4)*pow(x3,2.) + 
                                                             (y1 - 1.*y3)*pow(x4,2.)))*
                                                   ((x1 - 1.*x3)*y4*pow(x1,2.)*pow(x3,2.) - 1.*y2*pow(x1,3.)*pow(x3,2.) + 
                                                    y2*pow(x1,2.)*pow(x3,3.) + y2*pow(x1,3.)*pow(x4,2.) - 
                                                    1.*y3*pow(x1,3.)*pow(x4,2.) + y1*pow(x3,3.)*pow(x4,2.) - 
                                                    1.*y2*pow(x3,3.)*pow(x4,2.) + 
                                                    pow(x2,3.)*((-1.*y3 + y4)*pow(x1,2.) + (y1 - 1.*y4)*pow(x3,2.) + 
                                                                (-1.*y1 + y3)*pow(x4,2.)) - 1.*y2*pow(x1,2.)*pow(x4,3.) + 
                                                    y3*pow(x1,2.)*pow(x4,3.) - 1.*y1*pow(x3,2.)*pow(x4,3.) + 
                                                    y2*pow(x3,2.)*pow(x4,3.) + 
                                                    pow(x2,2.)*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                                (y1 - 1.*y3)*pow(x4,3.))) + 
                                                   4.*pow(x1*(x1 - 1.*x4)*x4*(x1 + x4)*(y2 - 1.*y3) + 
                                                          (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,3.) + 
                                                          (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,3.) + 
                                                          x3*((-1.*y2 + y4)*pow(x1,3.) + (-1.*y1 + y2)*pow(x4,3.)) + 
                                                          x2*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                              (y1 - 1.*y3)*pow(x4,3.)),2.),0.5));


  xextreme[1]=0.16666666666666666*pow(x1*(x1 - 1.*x4)*x4*(y2 - 1.*y3) + 
                                      (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,2.) + 
                                      (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,2.) + 
                                      x3*((-1.*y2 + y4)*pow(x1,2.) + (-1.*y1 + y2)*pow(x4,2.)) + 
                                      x2*((y3 - 1.*y4)*pow(x1,2.) + (-1.*y1 + y4)*pow(x3,2.) + (y1 - 1.*y3)*pow(x4,2.)),
                                      -1.)*(2.*(x1*(x1 - 1.*x4)*x4*(x1 + x4)*(y2 - 1.*y3) + 
                                                (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,3.) + 
                                                (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,3.) + 
                                                x3*((-1.*y2 + y4)*pow(x1,3.) + (-1.*y1 + y2)*pow(x4,3.)) + 
                                                x2*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                    (y1 - 1.*y3)*pow(x4,3.))) + 
                                            pow(-12.*(x1*(x1 - 1.*x4)*x4*(y2 - 1.*y3) + 
                                                      (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,2.) + 
                                                      (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,2.) + 
                                                      x3*((-1.*y2 + y4)*pow(x1,2.) + (-1.*y1 + y2)*pow(x4,2.)) + 
                                                      x2*((y3 - 1.*y4)*pow(x1,2.) + (-1.*y1 + y4)*pow(x3,2.) + 
                                                          (y1 - 1.*y3)*pow(x4,2.)))*
                                                ((x1 - 1.*x3)*y4*pow(x1,2.)*pow(x3,2.) - 1.*y2*pow(x1,3.)*pow(x3,2.) + 
                                                 y2*pow(x1,2.)*pow(x3,3.) + y2*pow(x1,3.)*pow(x4,2.) - 
                                                 1.*y3*pow(x1,3.)*pow(x4,2.) + y1*pow(x3,3.)*pow(x4,2.) - 
                                                 1.*y2*pow(x3,3.)*pow(x4,2.) + 
                                                 pow(x2,3.)*((-1.*y3 + y4)*pow(x1,2.) + (y1 - 1.*y4)*pow(x3,2.) + 
                                                             (-1.*y1 + y3)*pow(x4,2.)) - 1.*y2*pow(x1,2.)*pow(x4,3.) + 
                                                 y3*pow(x1,2.)*pow(x4,3.) - 1.*y1*pow(x3,2.)*pow(x4,3.) + 
                                                 y2*pow(x3,2.)*pow(x4,3.) + 
                                                 pow(x2,2.)*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                             (y1 - 1.*y3)*pow(x4,3.))) + 
                                                4.*pow(x1*(x1 - 1.*x4)*x4*(x1 + x4)*(y2 - 1.*y3) + 
                                                       (-1.*x4*y1 - 1.*x1*y3 + x4*y3 + x3*(y1 - 1.*y4) + x1*y4)*pow(x2,3.) + 
                                                       (x4*y1 + x1*y2 - 1.*x4*y2 - 1.*x1*y4)*pow(x3,3.) + 
                                                       x3*((-1.*y2 + y4)*pow(x1,3.) + (-1.*y1 + y2)*pow(x4,3.)) + 
                                                       x2*((y3 - 1.*y4)*pow(x1,3.) + (-1.*y1 + y4)*pow(x3,3.) + 
                                                           (y1 - 1.*y3)*pow(x4,3.)),2.),0.5));


  // ensure extremum is outside the range of the values used for interpolation
  if(xextreme[0]<=xpos[firstpos] && xextreme[0]>=xpos[lastpos] && xextreme[1]<=xpos[firstpos] && xextreme[1]>=xpos[lastpos]){
    mono=1;
  }
  else mono=0;

  // GODMARK TODO TEST DEBUG:
  // Forcing mono for (at least) field will lead to more kinks because of overlapping fixed values along interpolation that will require non-monotonicity to compress field together leading to smoother interpolation
  // mono=1;

  return(mono);

}



// check for monotonicity of using 3 point *extrapolation*
int checkmono3(int firstpos, int lastpos, FTYPE *xpos, FTYPE y0, FTYPE y1, FTYPE y2, FTYPE y3)
{
  FTYPE x0=xpos[0]; // x0 on edge (y0 not set, so not used)
  FTYPE x1=xpos[1];
  FTYPE x2=xpos[2];
  FTYPE x3=xpos[3]; // x3 on other edge

  FTYPE xextreme;

  int mono;

  // /home/jon/testing_nsbh.nb

  xextreme=0.5*(x1 - 1.*x2)*(x1 - 1.*x3)*(x2 - 1.*x3)*
    (-1.*(y2 - 1.*y3)*pow(x1,2.)*pow(x1 - 1.*x2,-1.)*pow(x1 - 1.*x3,-1.)*
     pow(x2 - 1.*x3,-1.) - 1.*(-1.*y1 + y3)*pow(x1 - 1.*x2,-1.)*pow(x2,2.)*
     pow(x1 - 1.*x3,-1.)*pow(x2 - 1.*x3,-1.) - 
     1.*(y1 - 1.*y2)*pow(x1 - 1.*x2,-1.)*pow(x1 - 1.*x3,-1.)*pow(x2 - 1.*x3,-1.)*
     pow(x3,2.))*pow(-1.*x3*y1 - 1.*x1*y2 + x3*y2 + x2*(y1 - 1.*y3) + x1*y3,-1.);

  
  // ensure extremum is outside the range of the values used for interpolation
  if(xextreme<=xpos[firstpos] && xextreme>=xpos[lastpos]){
    mono=1;
  }
  else mono=0;
  

  // GODMARK TODO TEST DEBUG:
  // Forcing mono for (at least) field will lead to more kinks because of overlapping fixed values along interpolation that will require non-monotonicity to compress field together leading to smoother interpolation
  //  mono=1;


  return(mono);
}








// bound CENTered quantities
int bound_prim_user_dir_nsbh(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int bound_prim_user_dir_nsbh_new(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  //  int bound_prim_user_dir_nsbh_old1(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  //  int bound_prim_user_dir_nsbh_old2(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  int boundreturn;




  boundreturn=0;
  // first get old type BCs
  // TEST DEBUG:  With _new() and removed avoidance of Bfield setting, still showed constant but wrong omegaf -- as if calling _old1() mattered while shouldn't.
  //  boundreturn+=bound_prim_user_dir_nsbh_old1(boundstage, finalstep, boundtime, whichdir, boundvartype, prim);
  // now overwrite with new parabolic interpolation type, but can turn on/off which prim[pl] this is done for
  boundreturn+=bound_prim_user_dir_nsbh_new(boundstage, finalstep, boundtime, whichdir, boundvartype, prim);


  return(boundreturn);
}



// special NS boundary code for centerered primitives
// averages all points within NS using parabolic interpolation
// TODO GODMARK: Note, could use face or corner analytical values at surface to pass interpolation through --  might be more stable by fixing extrapolation to surface value.  That would be more consistent with how A_i is treated.  Potentially how extrapolate into NS shouldn't matter too much, but apparently does alot (oscillating B_phi and OmegaF with this method vs. fixed wrong omegaf and B_\phi for older bound method)
// When DUP copy ends up with corner-like active cell to be copied from, then value can be quite off compared to nearer cells leading to lots of pointy behavior near surface.  Do like with A_dir and usecompact==1 so DUP result is more compact and doesn't miss changes near to surface.  That is, pointy behavior doesn't smooth-off and stays there for long time.  E.g. shear in U3 generates B3, but shear of U3 along surface from pointy copy means B3 that's generated preserved shear.
// One would think EMF that fixes U3 to NS surface rate would avoid that problem in long-term sense.
int bound_prim_user_dir_nsbh_new(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  static int firstwhichdir;
  static int firsttime=1;



  dualfprintf(fail_file,"BOUNDNEW: nstep=%ld steppart=%d boundtime=%21.15g : whichdir=%d\n",nstep,steppart,boundtime,whichdir);


  //////////////////
  //
  // Only need to do this once, not per-dimension, so setup which dimension can always trigger on (assumes dimensionality of problem doesn't change)
  //
  //////////////////
  if(firsttime==1){
    firsttime=0;
    firstwhichdir=whichdir;
  }



  if(WHICHPROBLEM==NSBH){
    //////////////////
    //
    // check if NS has sufficient grid depth for boundary conditions to be used for interpolations on smoothish functions
    //
    // if NS is moving, then must do this every substep and before boundary call since boundary call corresponds to at new time.
    //
    //////////////////
    //    check_nsdepth(boundtime);    // GODMARK: TODO: In 3D will need this to be time-dependent, but now that code for getting distances is too slow.
  }



  //////////////////////
  //
  // Set (time-dependent) primitives except pstag and field centered
  // copied from init.tools.c from init_primitives user1 function
  //
  //////////////////////


  if(whichdir==firstwhichdir){ // only need to set NS boundary conditions once -- not per dimension

    // DEBUG:
    dualfprintf(fail_file,"BOUNDTYPE: inittypeglobal=%d\n",inittypeglobal);

    
    int Nvec[NDIM];
    Nvec[1]=N1;
    Nvec[2]=N2;
    Nvec[3]=N3;

    bound_gen_deeppara(boundtime, prim, Nvec, NULL); // NULL = vpot since only bounding prim here


  }// end if doing first whichdir


  return(0);

 
}// end function







#define FULLINTERP 0
#define REDUCEDINTERP 1


// generalized boundary conditions function for deeppara method
void bound_gen_deeppara(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{


  //#pragma omp parallel //private(i,j,k,initreturn,whichvel,whichcoord) OPENMPGLOBALPRIVATEFULL
#pragma omp parallel OPENMPGLOBALPRIVATEFULL
  {
    int whichvel,whichcoord,i,j,k;
    int ii,jj,kk;
    int dir,dirstart,dirend;
    int odir1,odir2;
    int pos;
    int isinsideNSijk,hasmaskijk,hasinsideijk,reallyonsurfaceijk,cancopyfromijk;
    int deliter;
    FTYPE dumbinf;


    if(PRODUCTION==0){
      dumbinf=1.0;
      dumbinf/=0.0;
    }
    else dumbinf=0.0;



    // choose which thing interpolating
    if(vpot!=NULL){
      dirstart=1;
      dirend=3;
    }
    else{
      dirstart=0;
      dirend=0;
    }


    //////////////////
    // LOOP
    OPENMP3DLOOPVARSDEFINE;
    ////////  COMPFULLLOOP{
    OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      // loop over dirs.  dir=0 for CENT and dir=1,2,3 for CORN's.
      for(dir=dirstart;dir<=dirend;dir++){

        // get odir's
        get_odirs(dir,&odir1,&odir2);

        // set pos for computing geometry and primitive
        if(dir==0) pos=CENT;
        else if(dir==1) pos=CORN1;
        else if(dir==2) pos=CORN2;
        else if(dir==3) pos=CORN3;


 
        // Only bound if A_i  is such that a boundary cell, which occurs if *all* cells directly *touching* A_i are boundary cells.  Otherwise, should be set as surface value.
        isinsideNSijk=is_dir_insideNS(dir,i,j,k, &hasmaskijk, &hasinsideijk, &reallyonsurfaceijk, &cancopyfromijk);
        // isinsideNSijk=GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE);
        // if inside NS, then use shell values to form average (distance weighted) answer for boundary values

        if(isinsideNSijk){

          int interpproblem;
          int delorig[5];
          int del[5];
          int pliter,pl;
          int countonsurface[NDIM]={0};
          int countinsidesurface[NDIM]={0};
          int countiter;


          // get nearest neighbor for this dir=0 (CENT) or A_{dir}
          ii=i+GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 0);
          jj=j+GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 1);
          kk=k+GLOBALMACP0A1(nsmask,i,j,k,NSMASKCLOSEICORN1 + 3*(dir-1) + 2);


          interpproblem=0;
          // get lowest order offsets in each dimension
          interpproblem+=get_del(i,j,k,ii,jj,kk,delorig,del);


          // DEBUG:
          int spaceiter;
          SLOOPA(spaceiter) dualfprintf(fail_file,"PDFC2b: siter=%d : delorig=%d del=%d\n",spaceiter,delorig[spaceiter],del[spaceiter]);


          // now have offset for each cell to grab from starting ii,jj,kk position
          if(interpproblem==0){

            ///////////////
            //
            // get compactvalue corresponding to perpendicular interpolation across nearest points to i,j,k
            //
            ///////////////
            int usecompact,usecompactpl[NPR];
            FTYPE compactvalue[NPR];
            PLOOP(pliter,pl) compactvalue[pl]=dumbinf;
            FTYPE iref,jref,kref;
            int poscorn,icorn,jcorn,kcorn;
            usecompact=get_compactvalue(time, prim, Nvec, vpot,  dir,  pos, i,  j,  k,  ii,  jj,  kk, del, usecompactpl, compactvalue,&iref,&jref,&kref,&poscorn,&icorn,&jcorn,&kcorn);


            dualfprintf(fail_file,"BEFORE PRFULLDEL\n");


            ///////////////
            //
            // get interpolation along full del[1,2,3] direction
            // assumes optimal to interpolate along diagonal.  Sometimes, however, leads to kink along grid directions and code notices the kinks and dissipation is inevitable.
            //
            ///////////////
            int gotvaluedel[NDIM];
            int numpointsuseddel[NDIM];
            int gotfullvalue;
            FTYPE prfulldel[2][NPR];
            for(int rediter=0;rediter<2;rediter++) PLOOP(pliter,pl) prfulldel[rediter][pl]=dumbinf;
            int isongridfull[5],isinsidefull[5],hasmaskfull[5],hasinsidefull[5],reallyonsurfacefull[5],cancopyfromfull[5];
            FTYPE xposfull[5];
            int numpointsusedfulldel;
            int compactpointfull;
  
            gotvaluedel[0]=gotfullvalue=!get_fulldel_interpolation(time, prim, Nvec, vpot, dir, pos, i, j, k, ii, jj, kk, iref, jref, kref, poscorn, icorn, jcorn, kcorn, del, usecompact, compactvalue, isongridfull, isinsidefull, hasmaskfull, hasinsidefull, reallyonsurfacefull, cancopyfromfull, &numpointsusedfulldel, xposfull, &compactpointfull, prfulldel[0], prfulldel[1]);
            numpointsuseddel[0]=numpointsusedfulldel; // for later

            for(countiter=1;countiter<=numpointsusedfulldel;countiter++){// avoid 0 since that's i,j,k itself, and go to =numpoints since that's 1,2,3,...numpoints all used to get value at i,j,k
              //       if(isongridfull[countiter]==1 && isinsidefull[countiter]==0 && hasmaskfull[countiter]==1 && hasinsidefull[countiter]==1) countonsurface[0]++;
              if(isongridfull[countiter]==1 && reallyonsurfacefull[countiter]==1) countonsurface[0]++; // actually only restrict if really on the surface where value is fixed, else still flexible point.
              if(isinsidefull[countiter]==1) countinsidesurface[0]++; // gotvaluedel==0 if any points inside surface, so not needed yet
            }


     
            dualfprintf(fail_file,"BEFORE PRDEL123\n");
     
            ////////////////
            //
            // get interpolation along individual del directions if del!=0
            //
            ////////////////
            FTYPE prdel[2][NDIM][NPR];
            int isongriddel[NDIM][5],isinsidedel[NDIM][5],hasmaskdel[NDIM][5],hasinsidedel[NDIM][5],reallyonsurfacedel[NDIM][5],cancopyfromdel[NDIM][5];
            FTYPE xposdel[NDIM][5];
            int compactpointdel[NDIM];

            for(int rediter=0;rediter<2;rediter++) SLOOPA(deliter) PLOOP(pliter,pl) prdel[rediter][deliter][pl]=dumbinf;
     

            int countdeldimen=0;

            int ijk[NDIM]; ijk[1]=i; ijk[2]=j; ijk[3]=k;
     
            SLOOPA(deliter){ if(del[deliter]!=0) countdeldimen++; }
            if(countdeldimen>=2){ // only relevant (different than get_fulldel_interpolation() above that obtained prfulldel) if at least two del's are non-zero

              int delalt[NDIM][NDIM];
              int iijjkkalt[NDIM][NDIM];

              SLOOPA(deliter){// over grid-aligned directions (deliter refers to directions 1,2, or 3)

                // no need to do interpolation if no del in that direction
                if(del[deliter]==0) continue;

                // get odir's
                int odir1del,odir2del;
                get_odirs(deliter,&odir1del,&odir2del);

  
                delalt[deliter][deliter]=del[deliter]; // delalt[which direction][amount of del1,del2,del3 (1,2,3) for that direction] = [del along del-direction]
                delalt[deliter][odir1del]=delalt[deliter][odir2del]=0; // del=0 for orthogonal directions for interpolation along del-direction 

                // Any del's added will only be farther from NS, so always outside -- because NS has regular constant closed curvature.
                // iref,jref,kref,poscorn,icorn,jcorn,kcorn won't be used since usecompact=0
                int localusecompact=0; // compactvalue won't be used
                int iidel,jjdel,kkdel;
                int isinsidetemp,hasmasktemp,hasinsidetemp,reallyonsurfacetemp,cancopyfromtemp;

                // starting guess for iijjkkdel
                iidel = i + (ii-i)*(deliter==1);
                jjdel = j + (jj-j)*(deliter==2);
                kkdel = k + (kk-k)*(deliter==3);

                // new reference point for nearest neighbor, which if here, should be on surface or outside NS.  Since using single index ii or jj or kk generally won't be far enough along grid lines, search for valid point
                int numtries=0;
#define MAXNUMTRIES (Nvec[deliter])
                while(1){
                  isinsidetemp=is_dir_insideNS(dir, iidel, jjdel, kkdel, &hasmasktemp, &hasinsidetemp, &reallyonsurfacetemp, &cancopyfromtemp);
                  if(isinsidetemp==1){
                    iidel += (deliter==1)*(delalt[deliter][1]>=1 ? 1 : (delalt[deliter][1]<=-1 ? -1 : 0));
                    jjdel += (deliter==2)*(delalt[deliter][2]>=1 ? 1 : (delalt[deliter][2]<=-1 ? -1 : 0));
                    kkdel += (deliter==3)*(delalt[deliter][3]>=1 ? 1 : (delalt[deliter][3]<=-1 ? -1 : 0));
                  }
                  else{
                    // then done finding near-surface (strictly outside) point
                    break;
                  }
                  numtries++;
                  if(numtries>MAXNUMTRIES){
                    dualfprintf(fail_file,"Reached numtries=%d for dir=%d ijk=%d %d %d : iijjkkdel=%d %d %d\n",numtries,dir,i,j,k,iidel,jjdel,kkdel);
                    myexit(633254236);
                  }
                }

                // DEBUG:
                dualfprintf(fail_file,"dir=%d ijk=%d %d %d : deliter=%d odir12del=%d %d iijjkkdel=%d %d %d\n",dir,i,j,k,deliter,odir1del,odir2del,iidel,jjdel,kkdel);

  
                gotvaluedel[deliter]=!get_fulldel_interpolation(time, prim, Nvec, vpot, dir, pos, i, j, k, iidel, jjdel, kkdel, iref, jref, kref, poscorn, icorn, jcorn, kcorn, delalt[deliter], localusecompact, compactvalue, isongriddel[deliter], isinsidedel[deliter], hasmaskdel[deliter], hasinsidedel[deliter], reallyonsurfacedel[deliter], cancopyfromdel[deliter], &numpointsuseddel[deliter],  xposdel[deliter], &compactpointdel[deliter], prdel[0][deliter],prdel[1][deliter]);

                // iterate if other points used are actually on surface (only occurs for A_i)
                for(countiter=1;countiter<=numpointsuseddel[deliter];countiter++){ // avoid 0 point that is i,j,k itself, and go out to =numpoints since includes one of points to get value at i,j,k
                  //    if(isongriddel[deliter][countiter]==1 && isinsidedel[deliter][countiter]==0 && hasmaskdel[deliter][countiter]==1 && hasinsidedel[deliter][countiter]==1) countonsurface[deliter]++;
                  if(isongriddel[deliter][countiter]==1 && reallyonsurfacedel[deliter][countiter]==1) countonsurface[deliter]++; // only count of really on surface; trying to avoid fixed values
                  if(isinsidedel[deliter][countiter]==1) countinsidesurface[deliter]++; // gotvaluedel==0 if any points inside surface, so not needed yet
                }
              }


              // DEBUG:
              // compare various versions
              PLOOP(pliter,pl){
                if(vpot!=NULL) pl=0;

                dualfprintf(fail_file,"dir=%d ijk=%d %d %d : pl=%d : prfulldel=%21.15g prdel1=%21.15g prdel2=%21.15g prdel3=%21.15g\n",dir,i,j,k,pl,prfulldel[0][pl],prdel[0][1][pl],prdel[0][2][pl],prdel[0][3][pl]);
                dualfprintf(fail_file,"dir=%d ijk=%d %d %d : pl=%d : prfulldelreduced=%21.15g prdelreduced1=%21.15g prdelreduced2=%21.15g prdelreduced3=%21.15g\n",dir,i,j,k,pl,prfulldel[1][pl],prdel[1][1][pl],prdel[1][2][pl],prdel[1][3][pl]);

                if(vpot!=NULL) pliter=nprend+1;
              }


            }// end if at least 2 dimensions have del   


            dualfprintf(fail_file,"AFTER PRDEL123\n");






            ///////////////////////
            //
            // Determine which path to use for extrapolation
            //
            ///////////////////////     
     
     
            FTYPE prnew[2][NPR];

            // default
            PLOOP(pliter,pl) prnew[FULLINTERP][pl]=prfulldel[FULLINTERP][pl];



            if(countdeldimen>=2){
#if(0)
              // just average all versions (prfulldel, and any prdel's that are valid)
              PLOOP(pliter,pl) prnew[0][pl]=prfulldel[pl];
              if(del[1]!=0) PLOOP(pliter,pl) prnew[0][pl]+=prdel[1][pl];
              if(del[2]!=0) PLOOP(pliter,pl) prnew[0][pl]+=prdel[2][pl];
              if(del[3]!=0) PLOOP(pliter,pl) prnew[0][pl]+=prdel[3][pl];
              // normalize
              PLOOP(pliter,pl) prnew[0][pl]/=(1.0+(FTYPE)countdeldimen);

#elif(1)
              // give preference to those cases that have fewest fixed values
              for(int rediter=0;rediter<2;rediter++) PLOOP(pliter,pl) prnew[rediter][pl]=0.0;

              FTYPE totalweight[2]={0.0,0.0};
              int ignore[NDIM],ignoreiter;
              /////////////////////
              //
              // translate some prfull results into del[0] form for simplicity of routines below
              //
              /////////////////////
              if(usecompact==1){
                xposdel[0][0]=xposfull[0];
                xposdel[0][1]=xposfull[compactpointfull];
                xposdel[0][2]=xposfull[1];
                xposdel[0][3]=xposfull[2];
              }
              else{
                xposdel[0][0]=xposfull[0];
                xposdel[0][1]=xposfull[1];
                xposdel[0][2]=xposfull[2];
                xposdel[0][3]=xposfull[3];
              }
              ignore[0]=(gotvaluedel[0]==0);
              SLOOPA(ignoreiter) ignore[ignoreiter]=(del[ignoreiter]==0 || gotvaluedel[ignoreiter]==0);
              for(int rediter=0;rediter<2;rediter++) PLOOP(pliter,pl) prdel[rediter][0][pl]=prfulldel[rediter][pl];


              /////////////////////////////
              //
              // Get conditional averaged value of both rediter's
              //
              /////////////////////////////
              for(int rediter=0;rediter<2;rediter++){ // over normal and reduced value

                // also weight by distance to first point, so closer points are used if all equal counts for number of points used that are on surface (i.e. countonsurface)
                if(gotvaluedel[0]==1 && (countonsurface[0]<=countonsurface[1] || ignore[1]) && (countonsurface[0]<=countonsurface[2] || ignore[2]) && (countonsurface[0]<=countonsurface[3] || ignore[3])){
                  totalweight[rediter]+=1.0/(xposdel[0][1]*xposdel[0][1]);
                  PLOOP(pliter,pl) prnew[rediter][pl]+=prdel[rediter][0][pl]/(xposdel[0][1]*xposdel[0][1]);
                }
                if(gotvaluedel[1]==1 && del[1]!=0 && (countonsurface[1]<=countonsurface[0] || ignore[0]) && (countonsurface[1]<=countonsurface[2] || ignore[2]) && (countonsurface[1]<=countonsurface[3] || ignore[3])){
                  totalweight[rediter]+=1.0/(xposdel[1][1]*xposdel[1][1]);
                  PLOOP(pliter,pl) prnew[rediter][pl]+=prdel[rediter][1][pl]/(xposdel[1][1]*xposdel[1][1]);
                }
                if(gotvaluedel[2]==1 && del[2]!=0 && (countonsurface[2]<=countonsurface[0] || ignore[0]) && (countonsurface[2]<=countonsurface[1] || ignore[1]) && (countonsurface[2]<=countonsurface[3] || ignore[3])){
                  totalweight[rediter]+=1.0/(xposdel[2][1]*xposdel[2][1]);
                  PLOOP(pliter,pl) prnew[rediter][pl]+=prdel[rediter][2][pl]/(xposdel[2][1]*xposdel[2][1]);
                }
                if(gotvaluedel[3]==1 && del[3]!=0 && (countonsurface[3]<=countonsurface[0] || ignore[0]) && (countonsurface[3]<=countonsurface[1] || ignore[1]) && (countonsurface[3]<=countonsurface[2] || ignore[2])){
                  totalweight[rediter]+=1.0/(xposdel[3][1]*xposdel[3][1]);
                  PLOOP(pliter,pl) prnew[rediter][pl]+=prdel[rediter][3][pl]/(xposdel[3][1]*xposdel[3][1]);
                }
 
       
                if(totalweight[rediter]!=0.0) PLOOP(pliter,pl) prnew[rediter][pl]/=totalweight[rediter];
                else{
                  dualfprintf(fail_file,"Never got preference for countonsurface.  Shouldn't happen: counts=%d %d %d %d : gots=%d %d %d %d\n",countonsurface[0],countonsurface[1],countonsurface[2],countonsurface[3],gotfullvalue,gotvaluedel[1],gotvaluedel[2],gotvaluedel[3]);
                  myexit(48225252);
                }



                //  DEBUG:
                DLOOPA(deliter){
                  if(deliter==0 || del[deliter]!=0){
                    dualfprintf(fail_file,"rediter=%d deliter=%d countonsurface=%d ignore=%d gotvaluedel=%d totalweight=%21.15g\n",rediter,deliter,countonsurface[deliter],ignore[deliter],gotvaluedel[deliter],totalweight[rediter]);
                    for(spaceiter=0;spaceiter<=numpointsuseddel[deliter];spaceiter++) dualfprintf(fail_file,"point=%d xpos=%21.15g\n",spaceiter,xposdel[deliter][spaceiter]);
                  }
                }
       
              }



              if(1){
                //////////////
                //
                // Try to avoid kinks along all paths to i,j,k by manipulating result at i,j,k
                //
                ///////////////

                // Try using value that is closest to the conditional weighted averaged linear interpolation
                // the reduced average is the most conservative value, but least accurate if higher-order terms are required

                FTYPE error[NDIM][NPR];
                DLOOPA(deliter){
                  PLOOP(pliter,pl){
                    if(vpot!=NULL) pl=0;

                    if(deliter==0 || del[deliter]!=0) error[deliter][pl]=fabs(prdel[FULLINTERP][deliter][pl]-prnew[REDUCEDINTERP][pl]); // compare individual fully interpolated value with conditional weighted averaged linear interpolation value
                    else error[deliter][pl]=BIG;

                    // DEBUG:
                    dualfprintf(fail_file,"dir=%d ijk=%d %d %d : pl=%d deliter=%d error=%21.15g\n",dir,i,j,k,pl,deliter,error[deliter][pl]);

                    if(vpot!=NULL) pliter=nprend+1;
                  }
                }
  
                // now overwrite prnew[0] that is the fully interpolated value with value closest to "conditional weighted averaged linear" value
                int choose=0;
                PLOOP(pliter,pl){
                  if(vpot!=NULL) pl=0;

                  if(error[0][pl]<=error[1][pl] && error[0][pl]<=error[2][pl] && error[0][pl]<=error[3][pl]){ // <= so catches case when all errors are zero when reduced case used as full case
                    choose=0;
                  }
                  else if(error[1][pl]<=error[0][pl] && error[1][pl]<=error[2][pl] && error[1][pl]<=error[3][pl]){
                    choose=1;
                  }
                  else if(error[2][pl]<=error[0][pl] && error[2][pl]<=error[1][pl] && error[2][pl]<=error[3][pl]){
                    choose=2;
                  }
                  else if(error[3][pl]<=error[0][pl] && error[3][pl]<=error[1][pl] && error[3][pl]<=error[2][pl]){
                    choose=3;
                  }
                  else{
                    dualfprintf(fail_file,"Never got small error: dir=%d ijk=%d %d %d\n",dir,i,j,k);
                    myexit(34765634);
                  }

                  prnew[FULLINTERP][pl]=prdel[FULLINTERP][choose][pl];

                  // DEBUG:
                  dualfprintf(fail_file,"dir=%d ijk=%d %d %d : pl=%d choose=%d\n",dir,i,j,k,pl,choose);

                  if(vpot!=NULL) pliter=nprend+1;
                }// end over pl
              }  


#endif

            }







            //////////////
            //
            // store old values and assign final vpot/prim
            //
            /////////////
            FTYPE prold[NPR];
            if(vpot!=NULL){
              prold[0]=MACP1A0(vpot,dir,i,j,k);
              MACP1A0(vpot,dir,i,j,k) = prnew[FULLINTERP][0];
            }
            else{
              PLOOP(pliter,pl){
                prold[pl]=MACP0A1(prim,i,j,k,pl);
                MACP0A1(prim,i,j,k,pl) = prnew[FULLINTERP][pl];
              }
            }


            // DEBUG:
            if(vpot!=NULL){
              dualfprintf(fail_file,"dir=%d ijk=%d %d %d prold=%21.15g prnew=%21.5g\n",dir,i,j,k,prold[0],prnew[FULLINTERP][0]);
            }
            else{
              PLOOP(pliter,pl){
                dualfprintf(fail_file,"dir=%d ijk=%d %d %d pl=%d : prold=%21.15g prnew=%21.5g\n",dir,i,j,k,pl,prold[pl],prnew[FULLINTERP][pl]);
              }
            }

   

          }// end if interpproblem==0 from get_del()
          else{

            dualfprintf(fail_file,"Never found shell to use for: %d %d %d.  Means inside NS, but no nearby shell (can occur for MPI).  Assume if so, then value doesn't matter (not used), so revert to fixed initial value for NS.\n",i,j,k);


            if(vpot!=NULL){
              FTYPE vpotlocal[NDIM];
              get_vpot_fluxctstag_primecoords(time,i,j,k,prim,vpotlocal);
              MACP1A0(vpot,dir,i,j,k)=vpotlocal[dir];
            }
            else{
    
              int inittype;
              int initreturn;

              // evolve type (don't set inittypeglobal, set by init_primitives)
              inittype=0;

              initreturn=init_dsandvels(inittype, CENT, &whichvel, &whichcoord,time,i,j,k,MAC(prim,i,j,k),NULL);

              if(initreturn>0){
                FAILSTATEMENT("init.c:init_primitives()", "init_dsandvels()", 1);
              }
              else if(initreturn==0) MYFUN(transform_primitive_vB(whichvel, whichcoord, i,j,k, prim, NULL),"init.c:init_primitives","transform_primitive_vB()",0);
            }

          }//end if interpproblem==1

 
        }// end if inside NS
      }// end over dirs
    }// end loop block
  }// end parallel region


  return;
}




// whether to (by default) allow use of usecompact==1
// NOTEMARK: The problem with usecompact=1 with A_\phi (or anything on surface) is that the compact point will necessarily use all fixed points.  So in the end, interpolation will be through 2 fixed points instead of only 1 or 0, which constrains the interpolation to allow more kinky or jumpy behavior near surface.
//#define DEFAULTUSECOMPACT 1
#define DEFAULTUSECOMPACT 0


// Check if can use compact point instead of outer point and compute *compactvalue*
int get_compactvalue(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, int *del, int *usecompactpl, FTYPE *compactvalue, FTYPE *retiref, FTYPE *retjref, FTYPE *retkref, int *retposcorn, int *reticorn, int *retjcorn, int *retkcorn)
{
  int usecompact;
  int ialt[4],jalt[4],kalt[4]; // doesn't include ijk, so only 4 needed when going from 0..3
  struct of_geom geomdontusecorn;
  struct of_geom *ptrgeomcorn;
  struct of_geom geomdontuse[5];
  // must be array of pointers so each pointer can change
  struct of_geom *(ptrgeom[5]);
  int ptriter;
  int firstpos,lastpos;
  int pl,pliter;
  int isongridall,isinsideall;
  int isongrid[5],isinside[5],hasmask[5],hasinside[5],reallyonsurface[5],cancopyfrom[5];
  int yyy;
  int yyystart,yyyend;
  FTYPE xpos[5];
  FTYPE prpass[NPR];
  FTYPE localsingle[5][NPR];
  int constrained;
  int mono4[NPR],mono3left[NPR],mono3right[NPR];



  // must assign these pointers inside parallel region since new pointer per parallel thread
  for(ptriter=0;ptriter<5;ptriter++){
    ptrgeom[ptriter]=&(geomdontuse[ptriter]);
  }
  ptrgeomcorn=&geomdontusecorn;



  //default
  usecompact=0;
  PLOOP(pliter,pl) usecompactpl[pl]=0;


  // compact value position (all dir same since a difference)
  // Note that ijkref are based upon pos position
  FTYPE iref,jref,kref;
  iref=(FTYPE)ii-0.5*sign((FTYPE)del[1])*N1NOT1;
  jref=(FTYPE)jj-0.5*sign((FTYPE)del[2])*N2NOT1;
  kref=(FTYPE)kk-0.5*sign((FTYPE)del[3])*N3NOT1;

  

  // GODMARK TODO: special case for axisymmetry only
  // see if can use more compact (interpolated) point for parabolic extrapolation
  if(abs(del[1])==abs(del[2]) && abs(del[1])>0){ // GODMARK TODO: Special case


    // get geometry for iref,jref,kref so can unrescale below
    int icorn,jcorn,kcorn,poscorn;


    // only for axisym del case: GODMARK
    if(dir==0) poscorn=CORN3;
    else if(dir==1) poscorn=FACE1;
    else if(dir==2) poscorn=FACE2;
    else if(dir==3) poscorn=CENT;


    // convert ijkref in pos form to ijkcorn in poscorn form
    // all dir=0,1,2,3 are same for ijkcorn since ii,jj,kk formed for dir's pos and then any differences (0 or +1 or -1) are based upon differences
    //  if(iref>i) icorn=ii;
    //  else icorn=ii-1*ROUND2INT(sign((FTYPE)del[1]))*N1NOT1;
    //  if(jref>j) jcorn=jj;
    //  else jcorn=jj-1*ROUND2INT(sign((FTYPE)del[2]))*N2NOT1;
    //  if(kref>k) kcorn=kk;
    //  else kcorn=kk-1*ROUND2INT(sign((FTYPE)del[3]))*N3NOT1;


    // convert i,j,k for a given pos conversion
    // NOTE: not relying upon round do get to integer -- that's just to get to correct integer value without machine error issues
    if(pos==CENT && poscorn==CORN3){
      icorn=ROUND2INT(iref+0.5);
      jcorn=ROUND2INT(jref+0.5);
      kcorn=ROUND2INT(kref);
    }
    else if(pos==CORN3 && poscorn==CENT){
      icorn=ROUND2INT(iref-0.5);
      jcorn=ROUND2INT(jref-0.5);
      kcorn=ROUND2INT(kref);
    }
    else if(pos==CORN1 && poscorn==FACE1){
      icorn=ROUND2INT(iref+0.5);
      jcorn=ROUND2INT(jref-0.5);
      kcorn=ROUND2INT(kref);
    }
    else if(pos==CORN2 && poscorn==FACE2){
      icorn=ROUND2INT(iref-0.5);
      jcorn=ROUND2INT(jref+0.5);
      kcorn=ROUND2INT(kref);
    }
    else{
      dualfprintf(fail_file,"No such pos/poscorn combo: pos=%d poscorn=%d\n",pos,poscorn);
      myexit(81721352);
    }


    // seems all dir's are same for ijkalt since formed from differences
#define FUNCI(tj) (ROUND2INT(iref - ((jj-jref)/(ii-iref))*(tj-jref)) ) // rotated around iref,jref
#define FUNCJ(ti) (ROUND2INT(jref - ((ii-iref)/(jj-jref))*(ti-iref)) ) // rotated around iref,jref
 
    jalt[0]=jref-0.5-1.0 + 0;
    jalt[1]=jref-0.5-1.0 + 1;
    jalt[2]=jref-0.5-1.0 + 2;
    jalt[3]=jref-0.5-1.0 + 3;
       
    ialt[0]=FUNCI(jalt[0]);
    ialt[1]=FUNCI(jalt[1]);
    ialt[2]=FUNCI(jalt[2]);
    ialt[3]=FUNCI(jalt[3]);
     
    kalt[0]=kalt[1]=kalt[2]=kalt[3]=k;


    // DEBUG:
    dualfprintf(fail_file,"ijk=%d %d %d : iijjkk=%d %d %d : ijkref=%21.15g %21.15g %21.15g : ijkcorn=%d %d %d poscorn=%d ijkalt0=%d %d %d ijkalt1=%d %d %d ijkalt2=%d %d %d ijkalt3=%d %d %d\n",i,j,k,ii,jj,kk,iref,jref,kref,icorn,jcorn,kcorn,poscorn,ialt[0],jalt[0],kalt[0],ialt[1],jalt[1],kalt[1],ialt[2],jalt[2],kalt[2],ialt[3],jalt[3],kalt[3]);

     

    // ensure points are on grid while not inside NS
    // Since dealing with A_i before bounded, must avoid edges of grid where there can be nan's at this point
       
    isongridall=isinsideall=0; // init
    for(yyy=0;yyy<=3;yyy++){
      isongrid[yyy+1]=is_dir_onactivegrid(dir,ialt[yyy],jalt[yyy],kalt[yyy]);
      isinside[yyy+1]=is_dir_insideNS(dir,ialt[yyy],jalt[yyy],kalt[yyy], &hasmask[yyy+1],&hasinside[yyy+1],&reallyonsurface[yyy+1],&cancopyfrom[yyy+1]);
  
      isongridall+=isongrid[yyy+1]; // should add up to 4
      isinsideall+=isinside[yyy+1]; // should stay 0
    }


    int cando4=0;
    int cando3left=0;
    int cando3right=0;
    int cando2=0;
    int donecompact=0;

    // default in case all cando's will be zero, so by default avoid loop that sets values up
    yyystart=0;
    yyyend=-1;
      
    if(isongrid[2]==1 && isongrid[3]==1 && isinside[2]==0 && isinside[3]==0){ // then assume 2 points are on surface at corner of NS (index 2,3 for isongrid[] and isinside[])
      yyystart=1;
      yyyend=2;
      cando2=1;
    }
    if(isongrid[1]==1 && isongrid[2]==1 && isongrid[3]==1 && isinside[1]==0 && isinside[2]==0 && isinside[3]==0){// then can do 3 point interpolation
      yyystart=0;
      yyyend=2;
      cando3left=1;
    }
    if(isongrid[2]==1 && isongrid[3]==1 && isongrid[4]==1 && isinside[2]==0 && isinside[3]==0 && isinside[4]==0){ // then can do 3 point interpolation
      yyystart=1;
      yyyend=3;
      cando3right=1;
    }
    if(isongridall==4 && isinsideall==0){// then all points are good to use
  
      // below yyystart,yyyend overrides smaller stencils above
      yyystart=0;
      yyyend=3;
      cando4=1;
    }



    // positions and values along special path
    xpos[0]=0.0;


    for(yyy=yyystart;yyy<=yyyend;yyy++){
      // get value to interpolate
      if(vpot!=NULL) prpass[0]=MACP1A0(vpot,dir,ialt[yyy],jalt[yyy],kalt[yyy]);
      else PLOOP(pliter,pl) prpass[pl]=MACP0A1(prim,ialt[yyy],jalt[yyy],kalt[yyy],pl);


      // rescale (ptrgeom[0] never use to get compactvalue)
      get_geometry(ialt[yyy], jalt[yyy], kalt[yyy], pos, ptrgeom[yyy+1]);
      rescale_pl(time,dir,ptrgeom[yyy+1], prpass ,localsingle[yyy+1]);

      // reference point is between i,j,k and ii,jj,kk
      xpos[yyy+1]=sqrt((iref-ialt[yyy])*(iref-ialt[yyy]) + (jref-jalt[yyy])*(jref-jalt[yyy]) + (kref-kalt[yyy])*(kref-kalt[yyy]));
      if(yyy<=1) xpos[yyy+1]*=-1.0; // really negative relative to ijkref
    }



    ////////////////////////
    //
    // Attempt 4 point interpolation
    //
    ////////////////////////
    if(donecompact==0 && cando4==1){
      PLOOP(pliter,pl){
    
        if(vpot!=NULL) pl=0;
    
        // check if multiple points are constrained (fixed in time) to surface values, which will restrict the freedom of the field lines too much and require non-monotonic fit
        if(vpot!=NULL){
          constrained=0;
          //   if(isinside[1]==0 && hasmask[1]==1 && hasinside[1]==1) constrained++;
          //   if(isinside[2]==0 && hasmask[2]==1 && hasinside[2]==1) constrained++;
          //   if(isinside[3]==0 && hasmask[3]==1 && hasinside[3]==1) constrained++;
          //   if(isinside[4]==0 && hasmask[4]==1 && hasinside[4]==1) constrained++;
          if(reallyonsurface[1]==1) constrained++;
          if(reallyonsurface[2]==1) constrained++;
          if(reallyonsurface[3]==1) constrained++;
          if(reallyonsurface[4]==1) constrained++;
        }
        else constrained=0; // ignore constraint question for non-field

        if(constrained>=2) mono4[pl]=1; // avoid overconstraint for field
        else{
          // check for monotonicity of using 4 points
          // [0] is in middle, so order is really 1,2,0,3,4.  And only need to be monotonic between 2 & 3 since other points can actually be correctly non-monotonic.  Just don't want to introduce extra non-monotonicity between points
          firstpos=2; lastpos=3; // 0 is really at "2.5"
          mono4[pl]=checkmono4(firstpos,lastpos,xpos,localsingle[0][pl],localsingle[1][pl],localsingle[2][pl],localsingle[3][pl],localsingle[4][pl]);
        }

        if(mono4[pl]==1){

          // perform interpolation for these points if monotonic using 4 points
          localsingle[0][pl] =
            +localsingle[1][pl]*(xpos[0]-xpos[2])*(xpos[0]-xpos[3])*(xpos[0]-xpos[4])/((xpos[1]-xpos[2])*(xpos[1]-xpos[3])*(xpos[1]-xpos[4]))
            +localsingle[2][pl]*(xpos[0]-xpos[1])*(xpos[0]-xpos[3])*(xpos[0]-xpos[4])/((xpos[2]-xpos[1])*(xpos[2]-xpos[3])*(xpos[2]-xpos[4]))
            +localsingle[3][pl]*(xpos[0]-xpos[1])*(xpos[0]-xpos[2])*(xpos[0]-xpos[4])/((xpos[3]-xpos[1])*(xpos[3]-xpos[2])*(xpos[3]-xpos[4]))
            +localsingle[4][pl]*(xpos[0]-xpos[1])*(xpos[0]-xpos[2])*(xpos[0]-xpos[3])/((xpos[4]-xpos[1])*(xpos[4]-xpos[2])*(xpos[4]-xpos[3]))
            ;
      
          usecompactpl[pl]=1; // says compact value is good to go
          donecompact=1; // says no attempt any other weaker interpolations
        }// end if monotonic


        if(vpot!=NULL) pliter=nprend+1;

      }// end pl
   
    }
    else{
      PLOOP(pliter,pl) mono4[pl]=0; // so diags appear reasonable looking
    }


    ////////////////////////
    //
    // Attempt one of two 3 point interpolations
    //
    ////////////////////////
    if(donecompact==0 && (cando3left==1 || cando3right==1)){
      FTYPE xposnew[5],yfun[5];

      PLOOP(pliter,pl){
    
        if(vpot!=NULL) pl=0;
    
        // check for monotonicity of using 3 points
        if(cando3left==1){
          xposnew[0]=xpos[0]; yfun[0]=localsingle[0][pl]; // yfun[0] not used
          xposnew[1]=xpos[1]; yfun[1]=localsingle[1][pl];
          xposnew[2]=xpos[2]; yfun[2]=localsingle[2][pl];
          xposnew[3]=xpos[3]; yfun[3]=localsingle[3][pl];
          // check if multiple points are constrained (fixed in time) to surface values, which will restrict the freedom of the field lines too much and require non-monotonic fit
          if(vpot!=NULL){
            constrained=0;
            //     if(isinside[1]==0 && hasmask[1]==1 && hasinside[1]==1) constrained++;
            //     if(isinside[2]==0 && hasmask[2]==1 && hasinside[2]==1) constrained++;
            //     if(isinside[3]==0 && hasmask[3]==1 && hasinside[3]==1) constrained++;
            if(reallyonsurface[1]==1) constrained++;
            if(reallyonsurface[2]==1) constrained++;
            if(reallyonsurface[3]==1) constrained++;
          }
          else constrained=0; // ignore constraint question for non-field

          if(constrained>=2) mono3left[pl]=1; // avoid overconstraint for field
          else{
            // xpos[0] in middle.  Points are: 1,2,0,3 so check between 2 and 3
            firstpos=2; lastpos=3;
            if(dir==0) mono3left[pl]=checkmono3(firstpos,lastpos,xposnew,yfun[0],yfun[1],yfun[2],yfun[3]);
            else mono3left[pl]=1; // assume want field to always use mono (never any discontinuities)
          }
        }
        else mono3left[pl]=0;

        if(cando3right==1){
          xposnew[0]=xpos[0]; yfun[0]=localsingle[0][pl]; // yfun[0] not used
          xposnew[1]=xpos[2]; yfun[1]=localsingle[2][pl];
          xposnew[2]=xpos[3]; yfun[2]=localsingle[3][pl];
          xposnew[3]=xpos[4]; yfun[3]=localsingle[4][pl];
          // check if multiple points are constrained (fixed in time) to surface values, which will restrict the freedom of the field lines too much and require non-monotonic fit
          if(vpot!=NULL){
            constrained=0;
            //     if(isinside[2]==0 && hasmask[2]==1 && hasinside[2]==1) constrained++;
            //     if(isinside[3]==0 && hasmask[3]==1 && hasinside[3]==1) constrained++;
            //     if(isinside[4]==0 && hasmask[4]==1 && hasinside[4]==1) constrained++;
            if(reallyonsurface[2]==1) constrained++;
            if(reallyonsurface[3]==1) constrained++;
            if(reallyonsurface[4]==1) constrained++;
          }
          else constrained=0; // ignore constraint question for non-field

          if(constrained>=2) mono3right[pl]=1; // avoid overconstraint for field
          else{
            // xpos[0] in middle.  Points are: 2,0,3,4 so check between 2 and 3
            firstpos=1; lastpos=3;
            if(dir==0) mono3right[pl]=checkmono3(firstpos,lastpos,xposnew,yfun[0],yfun[1],yfun[2],yfun[3]);
            else mono3right[pl]=1; // assume want field to always use mono (never any discontinuities)
          }
        }
        else mono3right[pl]=0;

        if(mono3left[pl]==1 || mono3right[pl]==1){

          // perform interpolation for these points if monotonic using 3 points
          localsingle[0][pl] =
            +yfun[1]*(xposnew[0]-xposnew[2])*(xposnew[0]-xposnew[3])/((xposnew[1]-xposnew[2])*(xposnew[1]-xposnew[3]))
            +yfun[2]*(xposnew[0]-xposnew[1])*(xposnew[0]-xposnew[3])/((xposnew[2]-xposnew[1])*(xposnew[2]-xposnew[3]))
            +yfun[3]*(xposnew[0]-xposnew[1])*(xposnew[0]-xposnew[2])/((xposnew[3]-xposnew[1])*(xposnew[3]-xposnew[2]))
            ;

          usecompactpl[pl]=1; // says compact value is good to go
          donecompact=1; // says no attempt any other weaker interpolations
        }// end if monotonic


        if(vpot!=NULL) pliter=nprend+1;

      }// end pl
   
    }
    else{
      PLOOP(pliter,pl){
        mono3left[pl]=0; // so diags appear reasonable looking
        mono3right[pl]=0; // so diags appear reasonable looking
      }
    }
       


    ////////////////////////
    //
    // Attempt 2 point interpolation
    //
    ////////////////////////
    if(donecompact==0 && cando2==1){

      // perform (linear) interpolation for these points (only yyy=1,2 corresponding to localsingle[2,3] and xpos[2,3])
      PLOOP(pliter,pl){

        if(vpot!=NULL) pl=0;

        localsingle[0][pl] =
          +localsingle[2][pl]*(xpos[0]-xpos[3])/((xpos[2]-xpos[3]))
          +localsingle[3][pl]*(xpos[0]-xpos[2])/((xpos[3]-xpos[2]))
          ;
    
        usecompactpl[pl]=1;
        donecompact=1;
    
        if(vpot!=NULL) pliter=nprend+1;

      }// end pl
  
    }// end for cando2

     

    // check if ever did compact interpolation
    if(donecompact==0){
      usecompact=0;
      dualfprintf(fail_file,"WARNING: Reached donecompact=0 for all attempts, so no compact possible: ijk=%d %d %d dir=%d\n",i,j,k,dir);
    }
    else{


      // choose default
      usecompact=DEFAULTUSECOMPACT; // default is compactvalue was computed




      PLOOP(pliter,pl){
        if(vpot!=NULL) pl=0;

        if(usecompactpl[pl]==0) usecompact=0; // disable usecompact if couldn't do it for some reason for any pl

        if(vpot!=NULL) pliter=nprend+1;

      }
    }


    // transfer over compact interpolated value
    if(usecompact==1){
      // get ijkref geometry
      get_geometry(icorn, jcorn, kcorn, poscorn, ptrgeomcorn);
      // unrescale
      unrescale_pl(time,dir,ptrgeomcorn, localsingle[0],compactvalue); // so compactvalue is back to normal prim type (i.e. not rescaled)
    }


    // DEBUG:
    dualfprintf(fail_file,"dir=%d usecompactC=%d: ijk=%d %d %d : ijkref=%21.15g %21.15g %21.15g : can4332=%d %d %d %d\n",dir,usecompact,i,j,k,iref,jref,kref,cando4,cando3left,cando3right,cando2);
    for(yyy=0;yyy<=3;yyy++){
      dualfprintf(fail_file,"COMPACTMOREC: %d %d (%d %d) : ialt=%d jalt=%d kalt=%d\n",isongrid[yyy+1],isinside[yyy+1],isongridall,isinsideall,ialt[yyy],jalt[yyy],kalt[yyy]);
    }
    PLOOP(pliter,pl){
      if(vpot!=NULL) pl=0;

      int spaceiter;
      for(spaceiter=yyystart+1;spaceiter<=yyyend+1;spaceiter++){
        dualfprintf(fail_file,"COMPACTPDF2C: siter=%d : ison=%d isins=%d los[pl=%d]=%21.15g %21.15g\n",spaceiter,isongrid[spaceiter],isinside[spaceiter],pl,localsingle[spaceiter][pl],xpos[spaceiter]);
      }
      dualfprintf(fail_file,"COMPACTPDF3C: %ld %21.15g : pl=%d : los0=%21.15g : compactvalue=%21.15g : mono4=%d : mono3left=%d mono3right=%d :  usecompact=%d\n",nstep,t,pl,localsingle[0][pl],compactvalue[pl],mono4[pl],mono3left[pl],mono3right[pl],usecompactpl[pl]);

      if(vpot!=NULL) pliter=nprend+1;
    }
     



    // other things to return
    *retposcorn=poscorn;
    *reticorn=icorn;
    *retjcorn=jcorn;
    *retkcorn=kcorn;


  }// end if possible to use compact


  // other things to return
  *retiref=iref;
  *retjref=jref;
  *retkref=kref;



  return(usecompact);
}





#define CHECKUCON 1
// interpolate along del[1,2,3] with usecompact option
int get_fulldel_interpolation(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, FTYPE iref, FTYPE jref, FTYPE kref, int poscorn, int icorn, int jcorn, int kcorn, int *del, int usecompact, FTYPE *compactvalue, int *isongrid, int *isinside, int *hasmask, int *hasinside, int *reallyonsurface, int *cancopyfrom, int *numpointsused, FTYPE *xpos, int *compactpoint, FTYPE *prfulldel, FTYPE *prfulldelreduced)
{
  int iii,jjj,kkk;
  int iiii,jjjj,kkkk;
  int spaceiter;
  FTYPE localsingle[5][NPR];
  //  FTYPE xpos[5];
  int testi;
  int mono3[NPR];
  int constrained;
  struct of_geom geomdontuse[5];
  // must be array of pointers so each pointer can change
  struct of_geom *(ptrgeom[5]);
  FTYPE prpass[NPR];
  int firstpos,lastpos;
  int ptriter;
  int pliter,pl;



  // must assign these pointers inside parallel region since new pointer per parallel thread
  for(ptriter=0;ptriter<5;ptriter++){
    ptrgeom[ptriter]=&(geomdontuse[ptriter]);
  }

  // get geometry for i,j,k so can unrescale below
  get_geometry(i, j, k, pos, ptrgeom[0]);
 

  iii=ii+del[1];
  jjj=jj+del[2];
  kkk=kk+del[3];
   
  iiii=ii+2*del[1];
  jjjj=jj+2*del[2];
  kkkk=kk+2*del[3];


  // DEBUG:
  //     dualfprintf(fail_file,"iiijjjkkk=%d %d %d i4=%d %d %d\n",iii,jjj,kkk,iiii,jjjj,kkkk);


  // compactpoint tells (upon return) which point was used for the compact value (xpos, etc.)
  *compactpoint=3;
  // numpointsused is over all pl, so take max case below
  *numpointsused=0;


  // ensure points are on grid while not inside NS
  // Since dealing with A_i before bounded, must avoid edges of grid where there can be nan's at this point
  // For dir==0, not quite necessary, but ok since NS shouldn't be so close to edge of grid.
  isongrid[0]=is_dir_onactivegrid(dir,i,j,k);
  isongrid[1]=is_dir_onactivegrid(dir,ii,jj,kk);
  isongrid[2]=is_dir_onactivegrid(dir,iii,jjj,kkk);
  if(usecompact==0) isongrid[3]=is_dir_onactivegrid(dir,iiii,jjjj,kkkk);
  else isongrid[3]=MAX(is_dir_onactivegrid(dir,i,j,k),is_dir_onactivegrid(dir,ii,jj,kk));

  isinside[0]=is_dir_insideNS(dir,i,j,k, &hasmask[0],&hasinside[0],&reallyonsurface[0],&cancopyfrom[0]);
  isinside[1]=is_dir_insideNS(dir,ii,jj,kk, &hasmask[1],&hasinside[1],&reallyonsurface[1],&cancopyfrom[1]);
  isinside[2]=is_dir_insideNS(dir,iii,jjj,kkk, &hasmask[2],&hasinside[2],&reallyonsurface[2],&cancopyfrom[2]);
  if(usecompact==0) isinside[3]=is_dir_insideNS(dir,iiii,jjjj,kkkk, &hasmask[3],&hasinside[3],&reallyonsurface[3],&cancopyfrom[3]);
  else  isinside[3]=is_dir_insideNS(dir,ii,jj,kk, &hasmask[3],&hasinside[3],&reallyonsurface[3],&cancopyfrom[3]); // assume this already taken care of by usecompact==0,1



  if(isongrid[1]==1 && isongrid[2]==1 && isongrid[3]==1 && isinside[1]==0 && isinside[2]==0 && isinside[3]==0){
       
       
    // position along special path
    xpos[0]=0.0;
    xpos[1]=sqrt((i-ii)*(i-ii) + (j-jj)*(j-jj) + (k-kk)*(k-kk));
    xpos[2]=sqrt((i-iii)*(i-iii) + (j-jjj)*(j-jjj) + (k-kkk)*(k-kkk));
    if(usecompact==0) xpos[3]=sqrt((i-iiii)*(i-iiii) + (j-jjjj)*(j-jjjj) + (k-kkkk)*(k-kkkk));
    else xpos[3]=sqrt((i-iref)*(i-iref) + (j-jref)*(j-jref) + (k-kref)*(k-kref));
       
       
    // get geometry for ii,jj,kk
    get_geometry(ii, jj, kk, pos, ptrgeom[1]);
    // get geometry for iii,jjj,kkk
    get_geometry(iii, jjj, kkk, pos, ptrgeom[2]);
    // get geometry for iiii,jjjj,kkkk
    if(usecompact==0){
      get_geometry(iiii, jjjj, kkkk, pos, ptrgeom[3]);
    }
    else{
      get_geometry(icorn, jcorn, kcorn, poscorn, ptrgeom[3]);
    }




    // get rescaled prim to interpolate
    // get value to rescale
    if(vpot!=NULL) prpass[0]=MACP1A0(vpot,dir,ii,jj,kk);
    else PLOOP(pliter,pl) prpass[pl]=MACP0A1(prim,ii,jj,kk,pl);
    rescale_pl(time,dir,ptrgeom[1], prpass,localsingle[1]);

    if(vpot!=NULL) prpass[0]=MACP1A0(vpot,dir,iii,jjj,kkk);
    else PLOOP(pliter,pl) prpass[pl]=MACP0A1(prim,iii,jjj,kkk,pl);
    rescale_pl(time,dir,ptrgeom[2], prpass,localsingle[2]);

    if(usecompact==0){
      if(vpot!=NULL) prpass[0]=MACP1A0(vpot,dir,iiii,jjjj,kkkk);
      else PLOOP(pliter,pl) prpass[pl]=MACP0A1(prim,iiii,jjjj,kkkk,pl);
      rescale_pl(time,dir,ptrgeom[3], prpass,localsingle[3]);
    }
    else rescale_pl(time,dir,ptrgeom[3], compactvalue,localsingle[3]); // using compactvalue[pl] at ijkcorn and poscorn


       
       
    FTYPE reducedcase[NPR];
    int numpointsreduced;
    FTYPE error;
    error=0;
    // perform parabolic interpolation for these points
    PLOOP(pliter,pl){
      if(vpot!=NULL) pl=0;

      // check if multiple points are constrained (fixed in time) to surface values, which will restrict the freedom of the field lines too much and require non-monotonic fit
      if(vpot!=NULL){
        constrained=0;
        // if(isinside[1]==0 && hasmask[1]==1 && hasinside[1]==1) constrained++;
        // if(isinside[2]==0 && hasmask[2]==1 && hasinside[2]==1) constrained++;
        // if(isinside[3]==0 && hasmask[3]==1 && hasinside[3]==1) constrained++;
        if(reallyonsurface[1]==1) constrained++;
        if(reallyonsurface[2]==1) constrained++;
        if(reallyonsurface[3]==1) constrained++;
      }
      else constrained=0; // ignore constraint question for non-field

      if(constrained>=2){
        // constrained==1 would be normal extrapolation through 1 NS surface value for field or 1
        // then force assumption that mono3=1 so even if checkmono3 would give mono3=0, we allow for some non-monotonicity to allow freedom in field to avoid it being overconstrained and generating kinks
        mono3[pl]=1;
      }
      else{
        // check for monotonicity of using 3 point *extrapolation*
        // points must be in order
        // only need to avoid adding non-monotonicity in extrapolation region (so range from 0 -> 1 for usecompact==0 and 0->3 (really compact point) for usecompact==1)
        if(usecompact==0){ firstpos=0; lastpos=1; } // points are 0,1,2,3 so check between 0 and 1
        else { firstpos=0; lastpos=3; } // points are 0,3,1,2 so check between 0 and 3
        if(dir==0) mono3[pl]=checkmono3(firstpos,lastpos,xpos,localsingle[0][pl],localsingle[1][pl],localsingle[2][pl],localsingle[3][pl]);
        else mono3[pl]=1; // assume want field to always use mono (never any discontinuities)
      }


      if(mono3[pl]==1){
        *numpointsused=MAX(*numpointsused,3);

        localsingle[0][pl] =
          +localsingle[1][pl]*(xpos[0]-xpos[2])*(xpos[0]-xpos[3])/((xpos[1]-xpos[2])*(xpos[1]-xpos[3]))
          +localsingle[2][pl]*(xpos[0]-xpos[1])*(xpos[0]-xpos[3])/((xpos[2]-xpos[1])*(xpos[2]-xpos[3]))
          +localsingle[3][pl]*(xpos[0]-xpos[1])*(xpos[0]-xpos[2])/((xpos[3]-xpos[1])*(xpos[3]-xpos[2]))
          ;

      }// end if mono




      {// get linear and constant versions
        // GODMARK TODO: won't be continuous when mono switches!
        int other1,other2;
        if(usecompact==0) { other1=1; other2=2; }
        else { other1=3; other2=1; }

        // For interpolation of A_i in a dir\neq i, constant value is fine and may even be preferred
        if(
           ((dir==1 &&(!(del[2]==0&&del[3]==0)))&&(vpot!=NULL) ||
            (dir==2 &&(!(del[1]==0&&del[3]==0)))&&(vpot!=NULL) ||
            (dir==3 &&(!(del[1]==0&&del[2]==0)))&&(vpot!=NULL)) ||
           (vpot==NULL)
           )
          error=fabs(localsingle[other1][pl]-localsingle[other2][pl])/(fabs(localsingle[other1][pl])+fabs(localsingle[other2][pl]));
        else{
          error=0;
        }



        FTYPE caseA,caseB;
    
        // linear extrapolation
        caseA=
          +localsingle[other1][pl]*(xpos[0]-xpos[other2])/((xpos[other1]-xpos[other2]))
          +localsingle[other2][pl]*(xpos[0]-xpos[other1])/((xpos[other2]-xpos[other1]))
          ;
    
        caseB= +localsingle[other1][pl]; // nearest neighbor extrapolation
    
    
        // GODMARK TODO: free parameters.
        FTYPE switcherror = switcht(error,0.5,0.8);

        reducedcase[pl]=caseB*switcherror + (1.0-switcherror)*caseA;
        numpointsreduced=ROUND2INT(1*switcherror + (1.0-switcherror)*2); // estimate

      }// end getting linear and constant versions


      
      if(mono3[pl]==0){// then revert to using 2 points and so must be mono
        localsingle[0][pl]=reducedcase[pl];
        *numpointsused=MAX(*numpointsused,numpointsreduced);
      }


      if(vpot!=NULL) pliter=nprend+1;

    }// end over pl




    if(CHECKUCON && vpot==NULL){ // only makes sense to change U1-U3 if modifying prim
      // ignore returned value
      checkucon_modifyuu(time, prim, Nvec, vpot, dir, pos, i, j, k, ii, jj, kk, iii, jjj, kkk, iiii, jjjj, kkkk, xpos, localsingle, ptrgeom, usecompact);
    }


    //////////////
    //
    // unrescale final prim/vpot
    //
    //////////////
    unrescale_pl(time,dir,ptrgeom[0],localsingle[0], prfulldel);

    //get reduced case
    unrescale_pl(time,dir,ptrgeom[0],reducedcase, prfulldelreduced);


       
  }// end if ongrid and is not inside NS
  else{
    dualfprintf(fail_file,"UNEXPECTED: near boundary: dir=%d : ijk=%d %d %d : ijk2=%d %d %d ijk3=%d %d %d ijk4=%d %d %d\n",dir,i,j,k,ii,jj,kk,iii,jjj,kkk,iiii,jjjj,kkkk);
    *numpointsused=0;
    return(1);
  }


     
     
  // DEBUG:
  dualfprintf(fail_file,"PDFC: dir=%d ijk=%d %d %d : ijk2=%d %d %d : ijk3=%d %d %d ijk4=%d %d %d\n",dir,i,j,k,ii,jj,kk,iii,jjj,kkk,iiii,jjjj,kkkk);
  SLOOPA(spaceiter) dualfprintf(fail_file,"PDFC2: siter=%d : ison=%d isins=%d xpos=%21.15g\n",spaceiter,isongrid[spaceiter],isinside[spaceiter],xpos[spaceiter]);
  if(vpot!=NULL){
    pl=0;
    SLOOPA(spaceiter) dualfprintf(fail_file,"si=%d pl=%d : los=%21.15g\n",spaceiter,pl,localsingle[spaceiter][pl]);
    dualfprintf(fail_file,"PDFC3v: %ld t=%21.15g : pl=%d : los0=%21.15g : vpot=%21.15g : vpotreduced=%21.15g : mono3=%d\n",nstep,t,pl,localsingle[0][pl],prfulldel[pl],prfulldelreduced[pl],mono3[pl]);
  }
  else{
    SLOOPA(spaceiter) PLOOP(pliter,pl) dualfprintf(fail_file,"si=%d pl=%d : %21.15g\n",spaceiter,pl,localsingle[spaceiter][pl]);
    PLOOP(pliter,pl) dualfprintf(fail_file,"PDFC3p: %ld t=%21.15g : pl=%d : los0=%21.15g : prim=%21.15g : primreduced=%21.15g: mono3=%d\n",nstep,t,pl,localsingle[0][pl],prfulldel[pl],prfulldelreduced[pl],mono3[pl]);
  }

  // DEBUG:
  for(spaceiter=0;spaceiter<=*numpointsused;spaceiter++){
    dualfprintf(fail_file,"point=%d of <=%d : dir=%d ijk=%d %d %d : isinside=%d hasmask=%d hasinside=%d reallyonsurface=%d cancopyfrom=%d\n",spaceiter,*numpointsused,dir,i,j,k,isinside[spaceiter],hasmask[spaceiter],hasinside[spaceiter],reallyonsurface[spaceiter],cancopyfrom[spaceiter]);
  }



  return(0); // 0 indicates got value (i.e. no error)
}




// check that uu0 is reasonable value and if not, then modify how UU1-UU3 interpolated
int checkucon_modifyuu(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int dir, int pos, int i, int j, int k, int ii, int jj, int kk, int iii, int jjj, int kkk, int iiii, int jjjj, int kkkk, FTYPE *xpos, FTYPE localsingle[5][NPR], struct of_geom *(ptrgeom[5]), int usecompact)
{
  int baducononsurface;
  int ptriter;
  struct of_geom geomdontuse4ucon[5];
  struct of_geom *(ptrgeom4ucon[5]);
  FTYPE ucon[5][NDIM],others[5][NUMOTHERSTATERESULTS];
  int pliter,pl;



  // pos used for ptrgeom is correct for prim since not doing this for vpot
  for(ptriter=0;ptriter<5;ptriter++) ptrgeom4ucon[ptriter]=ptrgeom[ptriter];


  // get \gamma to ensure not going to copy-in bad (i.e. limited due to bsq low and gamma\sim GAMMAMAX) values.
  // both vpot and prim BC versions use prim for this!
  ucon_calc(MAC(prim,ii,jj,kk),ptrgeom4ucon[1],ucon[1],others[1]);
  ucon_calc(MAC(prim,iii,jjj,kkk),ptrgeom4ucon[2],ucon[2],others[2]);
  if(usecompact==0) ucon_calc(MAC(prim,iiii,jjjj,kkkk),ptrgeom4ucon[3],ucon[3],others[3]);
  else{
    FTYPE ucontemp1[NDIM],ucontemp2[NDIM];
    FTYPE otherstemp1[NUMOTHERSTATERESULTS],otherstemp2[NUMOTHERSTATERESULTS];
    
    ucon_calc(MAC(prim,i,j,k),ptrgeom4ucon[0],ucontemp1,otherstemp1);
    ucon_calc(MAC(prim,ii,jj,kk),ptrgeom4ucon[1],ucontemp1,otherstemp2);
    int jjiter;
    DLOOPA(jjiter) ucon[3][jjiter]=0.5*(ucontemp1[jjiter]+ucontemp2[jjiter]); // only need ucon[TT] below, so assume rest are only roughly accurate (if use others, would be bad near BH where u^\mu interpolations can lead to undefined velocities
    for(jjiter=0;jjiter<NUMOTHERSTATERESULTS;jjiter++) others[3][jjiter]=0.5*(otherstemp1[jjiter]+otherstemp2[jjiter]);
    
  }

  //////////
  //     
  // only use point (localsingle[0]) if not out of wack when doing FFDE (or lower GAMMAMAX if doing MHD)
  //
  // NOTE: If usecompact==1, then ucon[3][TT] is actually nearest to i,j,k.  So should treat ucon[1,2,3] as if could be at any position.  So 2^3=8 choices.
  //
  //////////
  baducononsurface=0;
  if(ucon[1][TT]<0.5*GAMMAMAX && ucon[2][TT]<0.5*GAMMAMAX && ucon[3][TT]<0.5*GAMMAMAX){ // <<<
    // then assume velocities are good as interpolated
  }
  else if(ucon[1][TT]<0.5*GAMMAMAX && ucon[2][TT]<0.5*GAMMAMAX  && ucon[3][TT]>=0.5*GAMMAMAX){ // <<>

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[1][pl]*(xpos[0]-xpos[2])/((xpos[1]-xpos[2]))
          +localsingle[2][pl]*(xpos[0]-xpos[1])/((xpos[2]-xpos[1]))
          ;
      }
    }
  }
  else if(ucon[1][TT]<0.5*GAMMAMAX && ucon[2][TT]>=0.5*GAMMAMAX && ucon[3][TT]>=0.5*GAMMAMAX){ // <>>

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[1][pl];
        ;
      }
    }
  }
  else if(ucon[1][TT]<0.5*GAMMAMAX && ucon[2][TT]>=0.5*GAMMAMAX && ucon[3][TT]<0.5*GAMMAMAX){ // <><

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[1][pl]*(xpos[0]-xpos[1])/((xpos[3]-xpos[1]))
          +localsingle[3][pl]*(xpos[0]-xpos[3])/((xpos[1]-xpos[3]))
          ;
      }
    }
  }
  else if(ucon[1][TT]>=0.5*GAMMAMAX && ucon[2][TT]<0.5*GAMMAMAX && ucon[3][TT]<0.5*GAMMAMAX){ // ><<

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[2][pl]*(xpos[0]-xpos[2])/((xpos[3]-xpos[2]))
          +localsingle[3][pl]*(xpos[0]-xpos[3])/((xpos[2]-xpos[3]))
          ;
      }
    }
  }
  else if(ucon[1][TT]>=0.5*GAMMAMAX && ucon[2][TT]<0.5*GAMMAMAX && ucon[3][TT]>=0.5*GAMMAMAX){ // ><>

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[2][pl];
        ;
      }
    }
  }
  else if(ucon[1][TT]>=0.5*GAMMAMAX && ucon[2][TT]>=0.5*GAMMAMAX  && ucon[3][TT]<0.5*GAMMAMAX){ // >><

    PLOOP(pliter,pl){
      if(pl>=U1 || pl<=U3){
        localsingle[0][pl] =
          +localsingle[3][pl];
        ;
      }
    }
  }
  else{ // >>>
    // just use bad copy -- nothing else to do
    baducononsurface=1;
  }


  // DEBUG:
  int spaceiter;
  SLOOPA(spaceiter) dualfprintf(fail_file,"siter=%d : ucon=%21.15g\n",spaceiter,ucon[spaceiter][TT]);


  return(baducononsurface);
}









// special NS boundary code for staggered fields
int bound_pstag_user_dir_nsbh(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  // v1)
  // 1) So first bound pstag fields such that B^n is fixed and Bperp1 & Bperp2 are copied from exterior
  // 2) Then maybe directly bound A_i that produces B^n , so then no monopoles created as NS moves around (for 3D case)
  //    also, otherwise don't know how to modify B^i.  If just refix B^n to desired solution, will mess up divB=0 just outside star!
  //    also, then less sensitive to how treat edges of boundary for NS since divb=0 for sure no mater what I do.
  // So bound full Aperp1 & Aperp2 on surface, but copy Apar (or leave alone and just copy Bperp1 and Bperp2 as above)
  // Requires EVOLVEWITHVPOT==1 (NO, not modifying pstag)
  //
  // old way requires setting B^i and A_i


  // v2)
  // For B^n = B^x : Copy A_x into BC from active grid
  //                 Linearly interpolate A_z into BC
  //                 Linearly interpolate A_y into BC
  // etc.
  //
  // This forces B^y and B^z to be close to active values so field is smooth through surface as good for boundary conditions to avoid dissipation at surface.
  // Note that A_perp1n & A_perp2n must be at least linearly extrapolated to obtain non-zero field at face.
  //
  // Otherwise, E_i even in active region just beyond NS surface will have dissipative term that don't want.
  // So A_i is dealt with at same time as when fix some components of A_i, rather than here.


  // v3)
  // Difficult to get A_i (e.g. A_1 and A_2) to be smooth into NS.  Ok to set Bstag^i inside NS for interpolation of Bstag^i -> Bcent^i as long as Bstag^i is later overwritten to use A_i versions so divB=0 is preserved in case NS moves.


  // v3 requires EVOLVEWITHVPOT>0 because change Bstag^i here as boundary condition inside NS so smooth creation of Bcent^i and need to recover Bstag^i that gives divB=0 and this requires reverting to Bstag^i from A_i.  Just in case NS moves.
  if(EVOLVEWITHVPOT==0){
    dualfprintf(fail_file,"Need EVOLVEWITHVPOT>0\n");
    myexit(869346311);
  }

  if(TRACKVPOT==0){
    dualfprintf(fail_file,"Good to have TRACKVPOT>0 since otherwise derived fieldcalc(smcalc) A_i looks bad even if B^i is good.\n");
    myexit(869346312);
  }









  return(0);
}





// DOSETFIX:
// 0 = don't set A_i to fixed values inside (and on surface of) NS
// 1 = do
// always should do this

// DOSETEXTRAPDIRECT
// 0 = don't extrapolate relevant A_i into NS (temporarily overwrites some fixed values)
// 1 = do
// like DOSETEXTRAP, but directly assigns values for each A_i at each grid point in such a way that is perfectly assigned even if multiple directions can interp/extrap to that point

// DOSETEXTRAP:
// 0 = don't extrapolate relevant A_i into NS (temporarily overwrites some fixed values)
// 1 = do
// restricted assignments so doesn't set A_i that might need to be fixed (which within same loop would be used as extrapolation rather than as fixed), but multiple assignments for interp/extrap possible, so not perfectly symmetric process.

// DOSETFINALFORCE:
// 0 = don't fix A_i on NS surface
// 1 = do
// should be done with DOSETEXTRAP==1 to fix-up NS surface values

// DONSBOUNDPLPR:
// 0 = don't bound p_l and p_r at flux faces
// 1 = only use one-sided interpolations from p_l,p_r
// 2 = like 1, but also force rho,u,v^i,B^i on NS surface
// problems with 1 || 2 if shocks/discontinuities near NS surface

// DOSETFIXEMF:
// 0 = don't fix E_i on NS surface
// 1 = do
// always do this, seems works well -- even better than withh DONSBOUNDPLPR or all things not fixed
// DOSETFIXEMF is like DONSBOUNDPLPR, but for EMFs since separately computed from normal face fluxes


// GODMARK: TODO: Solution appears to stabilize if only DOSETFIX==1, but at high density due to diffusion term on density flux.
// Can turn on things slowly
// 1) DOSETFIX==1 only up to t\sim 40M or roughly 1 orbital period at NS
// 2) DONSBOUNDPLPR=1 next, up to another 40M, so density becomes correct.  Somehow slowly turn this on over ~40M.
// 3) DOSETEXTRAP + DOSETFINALFORCE=1 next, up to another 40M so interpolation is really dissipationless for field.  Slowly turn on over ~40M
// So final solution should appear by around 40M+40M+40M = 120M.
// So need outer boundary of >120M


// GODMARK: TODO:
/////////////// 1) Problem with rotation being met: Try FFDE (done: immediately turns crusty in uu0 near NS when rotation is on, \Omega_F crusty and probably not matched.  Already problems at low omega.)
//                 But, also crusty in uu0 on NS surface (and in current sheet to BH) just before rotaion turned on.  Same spots where kinks in field are in both FFDE and MHD solutions.
/////////////// 2) Turn off BH gravity (just spherical polar Mink) (very much doubt problem since BH spin works)
/////////////// 3) Look at NS init/bounds code (done, got some ideas)
// 4) Problem with crustiness in uu0 near NS: Try lower resolution to see if crustier.  Try higher res to see if goes away (NOTE: NO early-time crustiness with FFDE!, so probably due to jumps in rho,u,v,B near NS as saw with NS by itself.)
// 4.5) So try MHD with new super-average for CENT BCs.
/////////////// 5) Try MC (done, field kink at surface on top-right of NS remains.  Rotation: no better at same time, but by t=154, kinda has some *small* regions outside NS that match.  Bad still though.)
/////////////// 6) Try RK4 (doubt problem)
/////////////// 7) Ensure LAXF (done)
/////////////// 8) Get NS init/bounds code working and try -- ensure STAG works (don't think that's problem)
/////////////// 9) Try spinning BH (done, works fine)
/////////////// 10) Estimate which field lines open up (draw, think: R_L\sim 14 -> r_L\sim 6+14\sim 20).  But indeed, interior may not open up!  But doesn't it have differential rotation?
// 11) Try sigma<1 with rotation
// 12) Try with only emf fixed -- i.e. what is origin of asymmetry?
// 13) Inspect EMF's in file with rotation.  too low?
/////////////// 14) Review parameters (numerical and physical) (nothing odd)
/////////////// 15) Maybe problem with TRACKVPOT or UPDATEVPOT (recall had problems before eventually in accretion  simulations).  Maybe update field in wrong place? (well, BH spin works)
// 16) Have to extrapolate deeper into NS so interpolation leads to no dissipation near NS.  Can't only set flux at NS surface. (done for CENT, still kinky and crusty at surface -- checking on rotation....  Still need to try for A_i as well.)
////////////// 17) Don't worry about proxy values inside NS being stationary, just ensure flux exactly at NS boundary agrees with stationarity.  interior are just proxy for interpolation for exterior. (done)
/////////////// 18) In bounds.ns.c: Copying B1 led to large uu0 on entire surface, while fitted extrapolation did well. (that was copying to face, don't do that and shouldn't -- and super-average now doing should lead to semi-decent extrapolation.).
// 19) Only linearly extrapolate if see that steady. B^r \Omega_F\sim B^\phi
/////////////// 20) Too low rotation to see on timescales simulated? (have to go to ~100M for 1 rotation, but Alfven wave moves at ~c, so should reach new solution over light crossing time as Alfven wave goes out (or cancels)
// 20) New super-average slow: precompute? (and limit to, say, 10 nearest cells?)
// 21) Use super-average for A_i, but need to average dA_i ?
// 22) Omega_F *inside* star is not constant for fixed BCs inside star -- problem?  Maybe surface of NS is set such that omegaf is not contsant?! (for FFDE)  Check MHD.  (Maybe just that copied ghost cells -- but what about before did that for testnsbh_ffde?)
//
//
// IMPORTANT:
// A) OMEGA: shows very different (and non-constant and even opposite signs) for \Omega_F computed from emf/B type calculations, while doing EMF!!  That's bad! (fixed -- was new ffde problem --- with OMEGADIFF=0 there are still same problems!?!?)
//////////// B) in init_dsandvels_nsbh(), averaging to CORN_i may be bad since doesn't avoid out-of-NS values, but still should give correct \Omega_F (CORN_i (i.e. emf) now uses correct shifted field)
//////////// C) So check stat/axi version of v^i ... error or wrong?  (no, was use of inversion bug with pr->U as fix)
//
//
// 23) uu0 large on pole and outer equator and pulls that into BCs.  Maybe avoid high gamma values for extrapolation (e.g. if near GAMMAMAX as if limited).  testnsbh_ffde didn't have bad consequences or as bad uu0 -- so shows that how interpolates CENT quantities into NS does matter!
// 24) Since how interpolate CENT quantities seems to matter still, use same types of tricks used for NS (and BH) when it was at r=0.  This includes avoiding cells that are bad, avoiding spikes in quantities, etc.
//     1) POLEDEATH to avoid suck near NS pole?  Odd that needed?  What does it mean?  Why treated like current sheet (i.e. why does uu0>>>1)?
///       minimize sucking by finding min "U2" around or crush "U2"=0.  B3 symmetric if 2D axi.
//     2) POLEGAMMADEATH0 -- force uu0 back to normal value even on grid?  Use limit_gamma()
//     3) extrapolate B_\phi (what I do: Bd3fromBu.nb) , gdet rho u^x1 (no)
//     4) Only use other value if gamma's match to 10% (from bounds.tools.c, see YFRAC12 and YDQFRAC12)
//     5) limit excessive slopes?  use MINM?
//     6) for densities: only use linear if exponentiation causes growth of value, not decreasion
//     7) linear for u^i
//     8) linearly interpolate gdet B^1 and gdet B^2 (kinda like what's already done by interpolating A_i, but would be equivalent to parabolic on A_i instead of linear)
// 25) Need analytical B^n for CORN_i positions rather than interpolating so EMF with fixed velocities is more accurate.  Note that B^n interpolation only uses fixed analytical values assuming A_i's on surface are fixed (and they are).


//  NOTES:
// 1) Ad3 as in vpotdump's is at loc=CORN3, so will appear non-symmetric around centered NS position if Ad3 is plotted as a centered quantity without any interpolation to CENT
// 2) Kink in Ad3 near NS on right side is reasonable given gravity pulls flux in and no allowed adjustment of internal field.


#define DOSETFIX 1
#define DOSETEXTRAPDIRECT 1
#define DOSETEXTRAP 0
#define DOSETFINALFORCE 0
#define DONSBOUNDPLPR 2
#define DOSETFIXEMF 1




// TIME SWITCHES
#if(0) //  DEBUG GODMARK TODO TEST

#define SWITCHT0 0.0
#define SWITCHT1 20.0
#define SWITCHT2 40.0 // by \sim 40, NS magnetosphere is (with or without rotation) settled to new geometry without drastic changes near NS
#define SWITCHT3 60.0
#define SWITCHT4 80.0
#define SWITCHT5 100.0
#define SWITCHT6 120.0

#else

#define SWITCHT0 0.0
#define SWITCHT1 0.0
#define SWITCHT2 0.0 // by \sim 40, NS magnetosphere is (with or without rotation) settled to new geometry without drastic changes near NS
#define SWITCHT3 0.0
#define SWITCHT4 0.0
#define SWITCHT5 0.0
#define SWITCHT6 0.0

#endif


// GODMARK TODO:
// 1) Need to control flux onto domain so v0~0 for closed field lines, corresponding to |J.B|/(|J||B|) \lesssim 0.1
// 2) Ramp-up rotation.
// 3) Ramp-up mass-loading (reasonable looking dipolar-like solution at high mass loading)
// 4) Need B_\phi\sim constant to avoid dissipative solution.  Force B_\phi\sim 0 until t\sim 40M after which extrapolate.  This allows damping and destruction of colliding Alfven waves, and then later non-zero B_\phi for "jet" regions (although why isn't setting flux at boundary enough?  maybe other cells farther outside NS see jump near NS and that leads to problems?)
// 5) Check that *F_{t\phi} and \Omega_F are \sim constant.


// special NS boundary code for A_i for staggered fields
// adjust A_i during evolution as a boundary condition
// After A_i is updated:
// force Aperp1=Aperp2 to be NS dipole field so B^n is good.  Forces only relaxed constraint to be Bperp1 & Bperp2 derived from differences on Aperp1 & Aperp2 or from Apar in full.
// Also where A_i is extrapolated (instead of only fixed), which takes care of need to set staggered fields and centered fields.  Both will be correct with evolved and then fixed-up A_i from this function
// NOTE: vpot[1-3][i][j][k]  not vpot[i][j][k][1-3]
// NOTE: get_vpot_fluxctstag_primecoords() gets A_i at CORNi, not at same location for staggered field.
void adjust_fluxctstag_vpot(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;
  //  void adjust_fluxctstag_vpot_dosetfix(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  //  void adjust_fluxctstag_vpot_dosetfix_new(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  void adjust_fluxctstag_vpot_dosetfix_newest(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  //  void adjust_fluxctstag_vpot_dosetextrapdirect(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  //  void adjust_fluxctstag_vpot_dosetextrapdirect_deep(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  void adjust_fluxctstag_vpot_dosetextrapdirect_deeppara(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  //  void adjust_fluxctstag_vpot_dosetfinalforce(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  //  void adjust_fluxctstag_vpot_dosetextrap(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);



  dualfprintf(fail_file,"adjustvpot: nstep=%ld steppart=%d fluxtime=%21.15g\n",nstep,steppart,fluxtime);


  if(DOSETFIX){
    //      adjust_fluxctstag_vpot_dosetfix_new(fluxtime,prim,Nvec,vpot); // ignorant of details about surface
    adjust_fluxctstag_vpot_dosetfix_newest(fluxtime,prim,Nvec,vpot); // directly uses is_dir_insideNS()
  }


  if(DOSETEXTRAPDIRECT){
    //      adjust_fluxctstag_vpot_dosetextrapdirect_deep(fluxtime,prim,Nvec,vpot);
    adjust_fluxctstag_vpot_dosetextrapdirect_deeppara(fluxtime,prim,Nvec,vpot);
  }
  


  /*   else{ */
  /*     // NOT REALLY CORRECT since should avoid surface A_i ! */

  /*     if(DOSETFIX){ */
  /*       adjust_fluxctstag_vpot_dosetfix(fluxtime,prim,Nvec,vpot); */
  /*     } */


  /*     //  if(DOSETEXTRAPDIRECT && t>SWITCHT1){// hard switch in time */
  /*     if(DOSETEXTRAPDIRECT){ */
  /*       adjust_fluxctstag_vpot_dosetextrapdirect(fluxtime,prim,Nvec,vpot); */
  /*     } */





  /*     if(DOSETEXTRAP && t>SWITCHT1){// hard switch in time */
  /*       //if(DOSETEXTRAP){ */
  /*       adjust_fluxctstag_vpot_dosetextrap(fluxtime,prim,Nvec,vpot); */
  /*     } */
  /*     if(DOSETFINALFORCE && t>SWITCHT1){// hard switch in time */
  /*       //  if(DOSETFINALFORCE){ */
  /*       adjust_fluxctstag_vpot_dosetfinalforce(fluxtime,prim,Nvec,vpot); */
  /*     } */
  /*   } */


  return;

}


// For those A_i components that are constrained to be fixed in time, fix A_i to the init_vpot_user() version fo A_i within the NS and on the NS surface
void adjust_fluxctstag_vpot_dosetfix_newest(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE vpotlocal[NDIM];
  int isinsidelocal[NDIM],hasmasklocal[NDIM],hasinsidelocal[NDIM],reallyonsurfacelocal[NDIM],cancopyfromlocal[NDIM];
  int fixedvalue[NDIM];


  //////////////
  //
  // Set surface of NS to analytical solution if possible.  Only carefully set surface values, not interior values that will be use deep (i.e. full interior) extrapolation
  //
  //////////////


  COMPFULLLOOP{

    for(int dir=1;dir<=3;dir++){

      if(Nvec[dir]==1){ // only directly fix to constant in time from t=0 if that dimension doesn't exist.  Otherwise, generally have to assume EMF_i determimes evolution (e.g. can then be stationary in time, but not fixed in time)
        // get i,j,k data
        isinsidelocal[dir]=is_dir_insideNS(dir, i, j, k, &hasmasklocal[dir], &hasinsidelocal[dir], &reallyonsurfacelocal[dir], &cancopyfromlocal[dir]);
        fixedvalue[dir]=reallyonsurfacelocal[dir];  //(isinsidelocal[dir]==0 && hasmasklocal[dir]==1 && hasinsidelocal[dir]==1);
      }
    }


    if(reallyonsurfacelocal[1]||reallyonsurfacelocal[2]||reallyonsurfacelocal[3]){ // if one of dir's is really on surface, then compute A_i at pos=CORN_i (i.e. natural position for each A_i)

      get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);

      for(int dir=1;dir<=3;dir++){ // at least one Nvec[dir] will be not 1
        // only fix if really on surface and dimension is ignorable so can fix A_i for all time to constant value
        if(Nvec[dir]==1 && reallyonsurfacelocal[dir]==1)   MACP1A0(vpot,dir,i,j,k)       =vpotlocal[dir];
      }
    }


  }// end i,j,k loop

}






// rescale factor to obtain more constant A*rescaleA dependence.
// Not subtracting off original A_i since that would be too speculative since some regions shift alot in angle (e.g. polar region where A_\phi goes through zero)
//
// This rescale factor makes difference near "equatorial" plane where otherwise (using linear extrapolation) field flips sign once just inside NS
FTYPE rescale_A(SFTYPE time, int dir, int i, int j, int k)
{
  int pos;
  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  FTYPE X[NDIM],V[NDIM];
  FTYPE X2[NDIM],V2[NDIM];

  // get NS parameters
  setNSparms(time,&rns,&omegak, &omegaf, &Rns,&Rsoft,&v0);


  // set pos for computing geometry and primitive (not done so far)
  if(dir==1) pos=CORN1;
  else if(dir==2) pos=CORN2;
  else if(dir==3) pos=CORN3;


  // determine rescale function
  bl_coord_ijk(i,j,k,pos,V);
  V[0]=time; FTYPE r=V[RR],th=V[TH],ph=V[PH],R=r*sin(th); // use time instead of V[0] or t since could be used at fluxtime instead of t or V[0]
  FTYPE x,y,z;



  FTYPE xns,yns,zns;
  FTYPE Vns[NDIM],Vcartns[NDIM];
  // Cartesian position
  //   x = r*cos(ph)*sin(th);
  //   y = r*sin(ph)*sin(th);
  //   z = r*cos(th);


  // get NS position
  FTYPE absrdiff;
  pos_NS(time, V, Vns, Vcartns, &absrdiff);
  //   xns=Vcartns[1];
  //   yns=Vcartns[2];
  //   zns=Vcartns[3];



  FTYPE rescaleA;

  if(dir==1){
    rescaleA=1.0;
  }
  else if(dir==2){
    rescaleA=1.0;
  }
  if(dir==3){
    // only rescale by radial (away from NS center) part of A_\phi (same for either NSSHAPE)
    
    // i.e. A_\phi \propto 1/rescaleA at t=0
    //    rescaleA=pow(absrdiff+Rsoft,1.0) ;
    // causes problems near poles where things change alot from dipolar

    rescaleA=1.0;
  }

  
  // return thing that one should multiply by in order to remove dependence
  // devide by it in order to recover dependence
  return(rescaleA);
}







// get lowest order offset spacing in each dimension between 2 grid points
int get_del(int i, int j, int k, int ii, int jj, int kk, int *delorig, int *del)
{
  int interpproblem;


  delorig[1]=ii-i;
  delorig[2]=jj-j;
  delorig[3]=kk-k;

  interpproblem=0;
  // rescale aspect ratio so that one of the del's is exactly 1 in absolute value
  // only valid for dimensions that exist
  // Using <= instead of = so catches equality case.  Doesn't create asymmetry since still same equal steps away from original ii,jj,kk point in the end
  if(N1NOT1==1 && delorig[1]!=0 && (delorig[2]==0 || fabs(delorig[1])<=fabs(delorig[2])) && (delorig[3]==0 || fabs(delorig[1])<=fabs(delorig[3]) ) ){ // then delorig[1] smallest in size
    del[1]=delorig[1]/abs(delorig[1]);
    del[2]=ROUND2INT((FTYPE)(delorig[2])/((FTYPE)fabs(delorig[1]))); // round to nearest integer
    del[3]=ROUND2INT((FTYPE)(delorig[3])/((FTYPE)fabs(delorig[1]))); // round to nearest integer
  }
  else if(N2NOT1==1 && delorig[2]!=0 && (delorig[1]==0 || fabs(delorig[2])<=fabs(delorig[1])) && (delorig[3]==0 || fabs(delorig[2])<=fabs(delorig[3]) ) ){ // then delorig[2] smallest non-zero in size
    del[2]=delorig[2]/abs(delorig[2]);
    del[1]=ROUND2INT((FTYPE)(delorig[1])/((FTYPE)fabs(delorig[2]))); // round to nearest integer
    del[3]=ROUND2INT((FTYPE)(delorig[3])/((FTYPE)fabs(delorig[2]))); // round to nearest integer
  }
  else if(N3NOT1==1 && delorig[3]!=0 && (delorig[1]==0 || fabs(delorig[3])<=fabs(delorig[1])) && (delorig[2]==0 || fabs(delorig[3])<=fabs(delorig[2]) ) ){ // then delorig[3] smallest in size
    del[3]=delorig[3]/abs(delorig[3]);
    del[1]=ROUND2INT((FTYPE)(delorig[1])/((FTYPE)fabs(delorig[3]))); // round to nearest integer
    del[2]=ROUND2INT((FTYPE)(delorig[2])/((FTYPE)fabs(delorig[3]))); // round to nearest integer
  }
  else interpproblem=1;


  if(del[1]==0 && del[2]==0 && del[3]==0){
    dualfprintf(fail_file,"All dels are zero: %d %d %d : %d %d %d\n",i,j,k,ii,jj,kk);
    myexit(363468346);
  }


  return(interpproblem);

}




// whether A_{dir} or EMF_{dir} is on *active* grid (i.e. instead of on boundary cells)
int is_dir_onactivegrid(int dir, int i, int j, int k)
{
  int isongrid;

  if(dir==0){ // CENT
    isongrid=(i>=0 && i<=N1-1 &&  j>=0 && j<=N2-1 &&  k>=0 && k<=N3-1);
  }
  else if(dir==1){ // CORN_1 (sits lower j-k)
    isongrid=(i>=0 && i<=N1-1 &&  j>=N2NOT1 && j<=N2-1 &&  k>=N3NOT1 && k<=N3-1);
  }
  else if(dir==2){ // CORN_2 (sits lower i-k)
    isongrid=(i>=N1NOT1 && i<=N1-1 &&  j>=0 && j<=N2-1 &&  k>=N3NOT1 && k<=N3-1);
  }
  else if(dir==3){ // CORN_3 (sits lower i-j)
    isongrid=(i>=N1NOT1 && i<=N1-1 &&  j>=N2NOT1 && j<=N2-1 &&  k>=0 && k<=N3-1);
  }

  return(isongrid);
}

// whether CENT or A_{dir} or EMF_{dir} is on full grid
// CENT and A_i same since A_i extends up to +1 in each direction compared to CENT
int is_ongrid(int dir, int i, int j, int k)
{
  int isongrid;

  if(dir==0){ // CENT
    isongrid=(i>=-N1BND && i<=N1M-1 &&  j>=-N2BND && j<=N2M-1 &&  k>=-N3BND && k<=N3M-1);
  }
  else if(dir==1){ // CORN_1
    isongrid=(i>=-N1BND && i<=N1M-1 &&  j>=-N2BND && j<=N2M-1 &&  k>=-N3BND && k<=N3M-1);
  }
  else if(dir==2){ // CORN_2
    isongrid=(i>=-N1BND && i<=N1M-1 &&  j>=-N2BND && j<=N2M-1 &&  k>=-N3BND && k<=N3M-1);
  }
  else if(dir==3){ // CORN_3
    isongrid=(i>=-N1BND && i<=N1M-1 &&  j>=-N2BND && j<=N2M-1 &&  k>=-N3BND && k<=N3M-1);
  }

  return(isongrid);
}



// similar to adjust_fluxctstag_vpot_dosetextrapdirect_deep(), but instead of super average over many nearby zones, choose nearest neighbor and it's path-like partners for parabolic interpolation along that path
// Only adjusts *interior* of NS, *not* surface since those can be updated from the EMFs.
// Or put another way, cannot predict A_i
void adjust_fluxctstag_vpot_dosetextrapdirect_deeppara(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{

  dualfprintf(fail_file,"adjustvpotdeep: nstep=%ld steppart=%d fluxtime=%21.15g\n",nstep,steppart,fluxtime);

  bound_gen_deeppara(fluxtime, prim, Nvec, vpot); // with vpot !=NULL, boundary conditions will be applied to only vpot and prim will be used if necessary to compute quantities


  return;
}














void adjust_fluxcttoth_vpot(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{

  // not used
}


void adjust_fluxcttoth_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] )
{

  // do nothing since assume using staggered fields

}









// special NS boundary code for EMFs for staggered fields
//void adjust_fluxctstag_emfs_new(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
// GODMARK TODO DEBUG
// NOTE: BOUNDPLPR will set EMF1,EMF2 for axisymmetric case since those are at faces.  And EMF3 just zero.  So will be exactly correct except zeroing-out EMF3 on boundary.  That's a good check that both functions are doing the same thing -- so more general 3D case is more likely bug free.
void adjust_fluxctstag_emfs(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

  // v1) force EMFperp1=EMFperp2=0 after computed EMF so that B^n stays fixed.
  // only good if object is stationary since EMF at midpoint in time, but still ok to do since A_i BC will fix-up the final field on the NS surface and interior

  // v2) No, don't want to set EMF to zero since EMF is at midpoint in time, so would force no update to that grid cell where NS was, while NS will be in new position invalidating exact EMF=0.  In essence, the presence of the NS creates a modified EMF (due to the conductive surface) that forces A_i on surface to be what NS wants.  So no need to force EMF itself, just force A_i and it will correspond to some effectively modified EMF (compared to flux update) for a given cell.


  // v3) No, erroneous EMF on NS boundary due to dissipative terms can lead to large fields just above surface.
  // Next timestep will take care of proper field evolution near NS given (e.g.) motion of NS through grid.



  int insideNS;
  int i,j,k,pl,pliter;
  struct of_gdetgeom gdetgeomdontuse;
  struct of_gdetgeom *ptrgdetgeom=&gdetgeomdontuse;
  int inittype;
  int dir,pos,odir1,odir2;
  int initreturn;
  int whichvel,whichcoord;
  FTYPE emf,gdet;
  FTYPE pr[NPR];
  int lll;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE Bcon[NDIM],ucon[NDIM],vcon[NDIM],others[NUMOTHERSTATERESULTS];
  FTYPE Bconl[NDIM],Bconr[NDIM];
  FTYPE Bcontouse[NDIM];
  int localisinside,localhasmask,localhasinside,localreallyonsurface,localcancopyfrom;
  int localin,odir1left,odir2left,odir12left;
  int m,l;


  dualfprintf(fail_file,"adjustemfs: nstep=%ld steppart=%d fluxtime=%21.15g\n",nstep,steppart,fluxtime);

  
  if(DOSETFIXEMF==0) return;


#if(0)
  // DEBUG:
  COMPLOOPP1{
    dualfprintf(fail_file,"IJKBAD1: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,1,i,j,k,B2));
    dualfprintf(fail_file,"IJKBAD2: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,2,i,j,k,B1));
    dualfprintf(fail_file,"IJKBAD3: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,1,i,j,k,B3));
    dualfprintf(fail_file,"IJKBAD4: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,2,i,j,k,B3));
  }   
#endif 




  // -1 means evolve type and set conditions as if inside NS, but don't have good interpolation for corner currently, so don't assume input is at loc==pos
  //  inittype=-1;

  // -3 means evolve type and set conditions as if inside NS, but only assume input primitive is good for B orthogonal to EMF_{dir} located at loc==pos.
  inittype=-3;



  // LOOP over i,j,k
  COMPFULLLOOP{


    // over dir (i.e. EMF's corners)
    for(dir=1;dir<=3;dir++){

      // get odir's
      get_odirs(dir,&odir1,&odir2);
      
      // skip if no such dimension relevant for non-zero EMF (assumes other dimension has all homogeneous terms (including geometry!)
      if(Nvec[odir1]==1 && Nvec[odir2]==1) continue;


      // set pos for computing geometry and primitive
      pos=CORN1+dir-1;
      
      // get where this pos is relative to inside and shell
      localisinside=is_dir_insideNS(dir, i, j, k, &localhasmask, &localhasinside, &localreallyonsurface, &localcancopyfrom);


      // EMF really on the NS surface
      if(localreallyonsurface==1){


        // get whether inside NS surrounding CORN_{dir} position
        localin=   GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE);
        odir1left= GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1,j-(odir1==2)*N2NOT1,k-(odir1==3)*N3NOT1,NSMASKINSIDE);
        odir2left= GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1,j-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1,NSMASKINSIDE);
        odir12left=GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir1==3)*N3NOT1-(odir2==3)*N3NOT1,NSMASKINSIDE);



        ////////////////////////
        //
        // Setup guess for primitive at CORN_{dir}
        //
        // prim contains pb for *CENT* value of primitive, which is a poor guess in terms of preserving symmetry.  But it'll be overwritten by init_dsandvels() below
        // But will be different B[dir] than from PLPR that uses directly interpolated value that is exactly in right place if one of odir's is zero.
        //
        ////////////////////////
        PLOOP(pliter,pl) pr[pl]=MACP0A1(prim,i,j,k,pl);

        // defaults
        Bcontouse[dir]=pr[B1+dir-1];
        Bcontouse[odir1]=pr[B1+odir1-1];
        Bcontouse[odir2]=pr[B1+odir2-1];


        ////////////////////////
        //
        // Now get B and v from actual interpolations
        //
        ////////////////////////
      

        ////////////////////////////////
        //
        // use actual field that was interpolated to CORN_i
        // Would have just computed these if doing evolution, which is only time EMF modified as in this function
        // don't worry about VEL since is forced within NS and on surface -- i.e. it's fixed while only originally have field at FACEs.
        // See fluxctstag.c for how pvbcorn accessed:
        // pvbcorninterp[dir = EMF_dir][i,j,k][which v or B component: odir1 or odir2][0,1 for v +- in odir1, NUMCS for B][0,1 for v +- in odir2 or B jump that is perp to both dir and component chosen);
        // 

        Bconl[odir1]=GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir1,NUMCS,0); // jump in odir2
        Bconr[odir1]=GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir1,NUMCS,1); // jump in odir2
        Bconl[odir2]=GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir2,NUMCS,0); // jump in odir1
        Bconr[odir2]=GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir2,NUMCS,1); // jump in odir1

 
        // pvcorn[which corner][which component in pl form][+-odir1][+-odir2]
        //  GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir2,m,l); // vel for jump in 
        //  GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir1,m,l); // vel
        // v^{odir2}  with jump in +-odir1 and +-odir2
        // v^{odir1}  with jump in +-odir1 and +-odir2
        //   pr[U1-1+odir2]=0.25*(GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir2,m,l);
        //   pr[U1-1+odir1]=GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir1,m,l);
 

        // Also deal with field that doesn't directly come into EMF, but indirectly can through velocity chosen.
        // Instead of defaulting to init_dsandvels_nsbh averaged B[dir]:
        // Better to use p_l and p_r from flux interpolation.  So will be consistent with flux interpolation for EMF.
        // CORN2 uses dir=FACE=1 and averages (or selects) in 3 (k,k-1) and averages (or selects) in 1 (l/r)
        // So in general, B[dir] is *from* CENT with respect to odir1,odir2 coordinates, but *located* at FACE[odir1,odir2,lr]
        // Determination of which B[dir] to use for averaging or selecting (of 8 choices or 4 l/r pairs)


        // TODO: Could compute normal field analytically at CORN_i and use that -- most correct.
        // GODMARK DEBUG: 
        if(0){
          Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
          Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]);
   
          for(m=0;m<NUMCS;m++){
            for(l=0;l<NUMCS;l++){
              pr[U1-1+odir2]+=0.25*GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir2,m,l);
              pr[U1-1+odir1]+=0.25*GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,odir1,m,l);
            }
          }

          // recall pl/pr are located at same location (same i,j,k) for a given interface.
          Bcontouse[dir]=(1.0/8.0)*(
                                    +GLOBALMACP1A1(gp_l,odir1,i,j,k,B1+dir-1)+GLOBALMACP1A1(gp_r,odir1,i,j,k,B1+dir-1)
                                    +GLOBALMACP1A1(gp_l,odir2,i,j,k,B1+dir-1)+GLOBALMACP1A1(gp_r,odir2,i,j,k,B1+dir-1)
                                    +GLOBALMACP1A1(gp_l,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1)+GLOBALMACP1A1(gp_r,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1)
                                    +GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1)+GLOBALMACP1A1(gp_r,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1)
                                    );

        }
        else{

          // Decide which corner is more inside NS (so more representative of desired surface field).  Might be better (more stable) than just averaging.
          // For field, want to choose surface value if possible.  If not possible, then must choose completely external field.  This defines what must be copied from exterior.  Average if more than one.
          // For velocity, interpolation comes from CENT, which will be overwritten by init_dsandvels().

          // one case can be analytically overridden, but other two cases really does mix internal and external A_{dir}
          // For A_i, always some mixture of l+r unless override with analytical solution in some cases
          // But when dealing with Bperp1 and Bperp2, copying from exterior makes sense since those B^i have no fixed values from surface of NS.
          // That is, could imagine A_dir from exterior and A_dir at surface and forming derivative at surface to get B^i, but ultimately the value of B^i is free and not fixed by NS.
          // 2^4-2=14 cases for each odir(1,2)
          // can't be 1 1 1 1 or 0 0 0 0 -- at least one cell has to inside and at least one has to be not inside


          // skip 1 1 1 1
          FTYPE totalweight;
          if     (localin==1 && odir1left==1 && odir2left==1 && odir12left==0){
            Bcontouse[odir1]=Bconl[odir1]; // one inside NS, choose surface one
            Bcontouse[odir2]=Bconl[odir2]; // one inside NS, choose surface one
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            // This will then behave like PLPR for setting flux when doing field-corrected fix for using prext[]
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==1 && odir2left==0 && odir12left==1){
            Bcontouse[odir1]=Bconl[odir1]; // one inside NS, choose surface one
            Bcontouse[odir2]=Bconr[odir2]; // one inside
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==1 && odir2left==0 && odir12left==0){
            Bcontouse[odir1]=Bconl[odir1]; // can't obtain from surface alone, so use exterior
            Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]); // both are on surface -- could choose analytical solution, but just average surface solution
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==0 && odir2left==1 && odir12left==1){
            Bcontouse[odir1]=Bconr[odir1]; // one is inside NS, so then choose surface one
            Bcontouse[odir2]=Bconl[odir2]; // one is inside NS, so then choose surface one
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==0 && odir2left==1 && odir12left==0){
            Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
            Bcontouse[odir2]=Bconl[odir2]; // can't obtain from surface alone, so use exterior
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==0 && odir2left==0 && odir12left==1){
            Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
            Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]);
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==1 && odir1left==0 && odir2left==0 && odir12left==0){
            Bcontouse[odir1]=Bconr[odir1]; // on one surface
            Bcontouse[odir2]=Bconr[odir2]; // on one surface
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          // recover 0 1 1 1
          else if     (localin==0 && odir1left==1 && odir2left==1 && odir12left==1){
            Bcontouse[odir1]=Bconr[odir1];
            Bcontouse[odir2]=Bconr[odir2];
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if     (localin==0 && odir1left==1 && odir2left==1 && odir12left==0){
            Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
            Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]);
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==0 && odir1left==1 && odir2left==0 && odir12left==1){
            Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
            Bcontouse[odir2]=Bconr[odir2]; // one inside, so use exterior one
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==0 && odir1left==1 && odir2left==0 && odir12left==0){
            Bcontouse[odir1]=Bconr[odir1]; // can get from one surface
            Bcontouse[odir2]=Bconl[odir2]; // can get from one surface
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==0 && odir1left==0 && odir2left==1 && odir12left==1){
            Bcontouse[odir1]=Bconr[odir1]; // can't get from surface, so use exterior
            Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]);
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i,j,k,B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==0 && odir1left==0 && odir2left==1 && odir12left==0){
            Bcontouse[odir1]=Bconl[odir1];
            Bcontouse[odir2]=Bconr[odir2];
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_l,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i,j,k,B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else if(localin==0 && odir1left==0 && odir2left==0 && odir12left==1){
            Bcontouse[odir1]=Bconl[odir1];
            Bcontouse[odir2]=Bconl[odir2];
            // B[dir] is from CENT (on odir1-odir2 plane that is CORNdir), so not naturally at pos, so must come from external since not directly settable on surface
            Bcontouse[dir]=0.0;  totalweight=0.0;
            if(Nvec[odir1]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir1,i-N1NOT1*(odir2==1),j-N2NOT1*(odir2==2),k-N3NOT1*(odir2==3),B1+dir-1); }
            if(Nvec[odir2]>1){ totalweight+=1.0; Bcontouse[dir]+=GLOBALMACP1A1(gp_r,odir2,i-N1NOT1*(odir1==1),j-N2NOT1*(odir1==2),k-N3NOT1*(odir1==3),B1+dir-1); }
            if(totalweight!=0.0) Bcontouse[dir]/=totalweight;
          }
          else{
            dualfprintf(fail_file,"bad set: dir=%d ijk=%d %d %d : %d %d %d %d\n",dir,i,j,k,localin,odir1left,odir2left,odir12left);
            dualfprintf(fail_file,"bad set2: %d %d %d %d %d\n",localisinside,localhasmask,localhasinside,localreallyonsurface,localcancopyfrom);
            myexit(198352546);
          }
        }// end else 1

        ////////////////////////////////
        //
        // NOTE: Bcontouse[dir] from face interpolation already striped of gdet
        //
        // Bcontouse[odir1/odir2] from corner interpolation may have gdet or not
        //

        // get "gdet" factor (really, whatever factor also used in fluxctstag.c:fluxcalc_fluxctstag_emf_1d() )
        // used to get true B^i since needed, and used later to get flux from emf
        get_geometry_geomeonly(i, j, k, pos, ptrgdetgeom); // get geometry at CORN[dir] where emf is located
        gdet=ptrgdetgeom->EOMFUNCMAC(B1-1+odir1); // which field ptrgeom->e doesn't matter as mentioned below


#if(CORNGDETVERSION==1)
        // then pvbcorninterp does *not* contain any extra gdet factors, just B^i
        pr[B1+odir1-1]=Bcontouse[odir1];
        pr[B1+odir2-1]=Bcontouse[odir2];
#else
        // then pvbcorninterp contains gdet factor, so \detg B^i  needs to be converted to B^i
        FTYPE igdetnosing=sign(gdet)/(fabs(gdet)+SMALL);
        pr[B1+odir1-1]=Bcontouse[odir1]*igdetnosing;
        pr[B1+odir2-1]=Bcontouse[odir2]*igdetnosing;
#endif
        // still need to get pr[B1+dir-1], which will be from interpolation in init_dsandvels_nsbh()  


        dualfprintf(fail_file,"EMFCOMPARE1: dir=%d ijk=%d %d %d : %d %d %d %d\n",dir,i,j,k,localin,odir1left,odir2left,odir12left);
        dualfprintf(fail_file,"EMFCOMPARE2: %d %d %d %d %d\n",localisinside,localhasmask,localhasinside,localreallyonsurface,localcancopyfrom);
        PLOOP(pl,pliter) dualfprintf(fail_file,"EMFCOMPARE3: pl=%d prim=%21.15g\n",pl,pr[pl]);

    
        ///////////////////////
        //
        // now constrain velocity (and possibly also field)
        //
        // get raw primitive in whichvel/whichcoord coordinates
        // Note: Don't need inputted (to function) prim[] since fixing values
        // prim is default field if interpolation-extrapolation to loc=pos fails to be contained entirely within NS
        // ok that don't fill "pstag" entry with pext[] type quantity, since not setting densities here and don't care about v||B
        ///////////////////////
        initreturn=init_dsandvels(inittype, pos, &whichvel, &whichcoord, fluxtime, i, j, k, pr, NULL);

        // if successfully got raw primitive, then transform raw primitive to WHICHVEL/PRIMECOORD primitive and compute EMF
        if(initreturn>0){
          FAILSTATEMENT("init.c:set_plpr_nsbh()", "init_dsandvels()", 1);
        }
        else if(initreturn==0){   
          // general location transformation for v and B from whichvel/whichcoord to WHICHVEL/(MCOORD->PRIMECOORDS)
          bl2met2metp2v_genloc(whichvel, whichcoord, pr, i, j, k, pos);


          // DEBUG:
          PLOOP(pliter,pl){
            dualfprintf(fail_file,"COMPAREEMF: dir=%d : ijk=%d %d %d : pl=%d pr=%21.15g\n",dir,i,j,k,pl,pr[pl]);
          }


          // set B^\mu
          Bcon[TT]=0.0;
          SLOOPA(lll) Bcon[lll] = pr[B1+lll-1];


#if(0)
          // get u^\mu
          get_geometry(i, j, k, pos, ptrgeom);
          ucon_calc(pr, ptrgeom, ucon, others);

          // set v^i
          vcon[TT]=0.0;
          SLOOPA(lll) vcon[lll] = ucon[lll]/ucon[TT];


          // now set EMF_i
          // see fluxct.c for signature of emf[v,B] as related to fluxes
          //
          // so:
          // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
          // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
          // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]
          //
          // get_odirs(): dir=3 gives odir1=1 and odir2=2
          // so cyclic in: dir,odir1,odir2
          emf = Bcon[odir2] * vcon[odir1] - Bcon[odir1] * vcon[odir2];
#else

          // bit excessive for flux since only need 1 of them not all 3 or whole matrix that's computed internally
          FTYPE flux[NPR];
          struct of_state q;
          get_geometry(i, j, k, pos, ptrgeom);
          get_state(pr, ptrgeom, &q) ;
          dualfprintf(fail_file,"gdet=%21.15g : B1=%21.15g B2=%21.15g B3=%21.15g uu0=%21.15g uu1=%21.15g uu2=%21.15g uu3=%21.15g\n",gdet,pr[B1],pr[B2],pr[B3],q.ucon[TT],q.ucon[1],q.ucon[2],q.ucon[3]);
          if(Nvec[odir1]>1){
            dualfaradayspatial_calc(pr,odir1,&q,&flux[B1]); // fills B1->B3
            emf=flux[B1+odir2-1];
          }
          else{
            dualfaradayspatial_calc(pr,odir2,&q,&flux[B1]); // fills B1->B3
            emf=-flux[B1+odir1-1];
          }
       
#endif


          // DEBUG:
          if(Nvec[odir1]>1) dualfprintf(fail_file,"1EMF(%d): ijk=%d %d %d : old: %21.15g  new: %21.15g\n",dir,i,j,k,MACP1A1(fluxvec,odir1,i,j,k,B1-1+odir2)/gdet,emf);
          if(Nvec[odir2]>1) dualfprintf(fail_file,"2EMF(%d): ijk=%d %d %d : old: %21.15g  new: %21.15g\n",dir,i,j,k,MACP1A1(fluxvec,odir2,i,j,k,B1-1+odir1)/gdet,-emf);


#if(1)
          if(dir==1 || dir==2){
            // DEBUG: CHECK that \Omega_F is as expected
            // get NS parameters
            FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
            setNSparms(fluxtime, &rns, &omegak, &omegaf,  &Rns, &Rsoft, &v0);

            FTYPE dxdxp[NDIM][NDIM];
            dxdxprim_ijk(i,j,k,pos,dxdxp);

            FTYPE omegaftest;
            if(dir==1){
              omegaftest=(-emf/Bcon[2])*dxdxp[PH][PH];
            }
            if(dir==2){
              omegaftest=(emf/Bcon[1])*dxdxp[PH][PH];
            }

            dualfprintf(fail_file,"OMEGA: %d %d %d : %21.15g : %21.15g =?= %21.15g\n",i,j,k,omegak,omegaf,omegaftest);
            if(fabs(omegaftest-omegaf)>1E-5){
              dualfprintf(fail_file,"OMEGADIFF: %d %d %d : %21.15g : %21.15g =?= %21.15g\n",i,j,k,omegak,omegaf,omegaftest);
            }
          }
#endif

    
          // assign
          if(Nvec[odir1]>1) MACP1A1(fluxvec,odir1,i,j,k,B1-1+odir2) = + emf*gdet;
          if(Nvec[odir2]>1) MACP1A1(fluxvec,odir2,i,j,k,B1-1+odir1) = - emf*gdet;
          if(Nvec[dir]>1) MACP1A1(fluxvec,dir,i,j,k,B1-1+dir)     =   0.0;
        }//end if initreturn==0
        // else initreturn<0 then no change (i.e. nothing to set for boundary condition)
      }// end if EMF_{dir} is on inside of NS or its surface
    }// end over dir
  }// end i,j,k LOOP



  return;


}






// special NS boundary code for modifying/remapping primitives at flux interfaces
// Note that this remap or setting of p2interp_l and p2interp_r assumes RESCAPEINTERP==0 *or* that are setting rescaled quantities directly.
void remapplpr_nsbh( int dir, int idel, int jdel, int kdel, int i, int j, int k, 
                     FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], 
                     FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
                     FTYPE *p2interp_l, FTYPE *p2interp_r )
{


  // Not used


}


// special NS boundary code for primitives at interfaces used to compute fluxes
// Needs BOUNDPLPR==1
// p_l and p_r are normal set of primitives (not rescaled)
// as shown in flux.c, the positions are such that:  |   i-1  p_l|p_r  i   |
// So p_l and p_r are on lower edge always
void set_plpr_nsbh(int dir, SFTYPE fluxtime, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r)
{
  int initreturn,inittype;
  int whichvel,whichcoord;
  int pliter,pl;
  FTYPE pr[NPR], prother[NPR], pstag[NPR];
  FTYPE pnew_l[NPR],pnew_r[NPR];
  FTYPE prext[NPR];



  dualfprintf(fail_file,"setplpr: nstep=%ld steppart=%d dir=%d fluxtime=%21.15g ijk=%d %d %d\n",nstep,steppart,dir,fluxtime,i,j,k);



  if(DONSBOUNDPLPR==0) return;


  // -2 means evolve type and set conditions as if inside NS, and input values are at inputted pos (so can avoid excessive interpolation inside init_dsandvels_nsbh())
  inittype=-2;


  int isinside,hasmask,hasinside,reallyonsurface,cancopyfrom;
  // get surface position
  int pos=FACE1+dir-1;
  int faces[3];
  isinside=is_pos_insideNS(pos,i, j, k, &hasmask, &hasinside, &reallyonsurface, &cancopyfrom,faces);
  int faceinside,faceshell,faceboth;
  faceinside=faces[0];
  faceshell=faces[1];
  faceboth=faces[2];


  
  //  if(faceinside || faceshell || faceboth){    // if here, then this boundary is a NS surface boundary, so set primitives here
  if(reallyonsurface){    // if here, then this pos (dir) is a on NS surface, so set primitives so flux is set (previously also set fluxes inside NS, but not necessary and wasteful)


    ///////////////////
    //
    // Use one-sided interpolations to set NS boundary values
    //
    ///////////////////
    if(DONSBOUNDPLPR>=1){
      // set primitives
      // not using inputted (into function) prim[] since pr is to be set
      // internally uses pglobal and pstagglobal for field
      // p_l or p_r is default field if interpolation-extrapolation to loc=pos fails to be contained entirely within NS
      // decide if should use interpolated face value (say one-side) or use interpolation provided by init_dsandvels_nsbh that tries to exclude outside NS values of field
      // default first (used only when no other possibility.  Corresponds to pb @ CENT -- so bad for symmetry)
      PLOOP(pliter,pl) prext[pl]=pr[pl]=MACP0A1(prim,i,j,k,pl);
      // used code's one-sided interpolation
      if(faceinside){
        PLOOP(pliter,pl){
          pr[pl]=p_r[pl]; // assume want values inside NS for fixed quantities
          prext[pl]=p_l[pl]; // assume want values outside NS for outflow quantities
        }
      }
      if(faceshell){
        PLOOP(pliter,pl){
          pr[pl]=p_l[pl]; // assume want values inside NS
          prext[pl]=p_r[pl];
        }
      }
      if(faceboth){
        PLOOP(pliter,pl){
          pr[pl]=0.5*(p_l[pl]+p_r[pl]);
          prext[pl]=pr[pl]; // force prext as same as pr
        }
      }


      // fix pr so that prext is used for quantities that should be effectively copied from exterior to NS
      // for dir==1, B2,B3 should come from exterior, etc.
      // Actually, more generally, Bperp1 and Bperp2 come from mixter of values from the surface and exterior for A_i, so this may be overdoing it.
      // But seems ok, since (at least B3) has no fixed value and comes from derivatives in A_i.  So when dealing with B^i (not A_i), this copying from exterior alone seems reasonable.
      if(1){
        // This seems required to have stable (albeit inaccurate) omegaf.  Probably related to B3 and dir==1 and dir==2.
        // But still not correct value of omegaf -- jumps at surface.  B3 has jump?
        // Note, I wasn't doing this before with old1 bound and oscillated as well.  Maybe oscillations related to para for B3 but old bound fixes B3 bit more so doesn't oscillate when use those CENT values here.

        // So don't have to do this inside init_dsandvels_nsbh() for DONSBOUNDPLPR>=2 below
        if(dir==1){
          pr[B2]=prext[B2];
          pr[B3]=prext[B3];
        }
        else if(dir==2){
          pr[B1]=prext[B1];
          pr[B3]=prext[B3];
        }
        else if(dir==3){
          pr[B1]=prext[B1];
          pr[B2]=prext[B2];
        }
      }


      // DEBUG:
      PLOOP(pliter,pl){
        dualfprintf(fail_file,"COMPAREPLPR: dir=%d : ijk=%d %d %d : fff=%d %d %d : pl=%d prim=%21.15g pr=%21.15g prext=%21.15g (p_l=%21.15g p_r=%21.15g)\n",dir,i,j,k,faceinside,faceshell,faceboth,pl,MACP0A1(prim,i,j,k,pl),pr[pl],prext[pl],p_l[pl],p_r[pl]);
      }


    }// end if doing BOUNDPLPR>=1

    


    ///////////////////
    //
    // Use strict NS boundary values (required to have (e.g.) \Omega_F fixed for non-EMF fluxes)
    //
    // pr -> prother -> strict NS prother
    //
    ///////////////////
    if(DONSBOUNDPLPR>=2){


      PLOOP(pliter,pl) prother[pl]=pr[pl];

      initreturn=init_dsandvels(inittype, pos, &whichvel, &whichcoord, fluxtime, i, j, k, prother, prext); // only time really input prext, to be used as value that is external to NS if required. -- so far not required.

      if(initreturn>0){
        FAILSTATEMENT("init.c:set_plpr_nsbh()", "init_dsandvels()", 1);
      }
      else if(initreturn==0){

        // DEBUG:
        PLOOP(pliter,pl) dualfprintf(fail_file,"COMPAREPLPR2.0[%d]=%21.15g\n",pl,prother[pl]);

        // general location transformation for v and B from whichvel/whichcoord to WHICHVEL/(MCOORD->PRIMECOORDS)
        bl2met2metp2v_genloc(whichvel, whichcoord, prother, i, j, k, pos);

      }

      PLOOP(pliter,pl) pr[pl]=prother[pl];

      // DEBUG:
      PLOOP(pliter,pl) dualfprintf(fail_file,"COMPAREPLPR2.5[%d]=%21.15g\n",pl,pr[pl]);


      // else initreturn<0 then no change (i.e. nothing to set for boundary condition)
    }// end if using strict NS boundary values
   



    ///////////////
    //
    // pr -> pnew_l,pnew_r at same physical position
    //
    ///////////////
    PLOOP(pliter,pl) pnew_l[pl]=pnew_r[pl]=pr[pl];
    //      PLOOPNOB1(pl) p_l[pl]=p_r[pl]=pr[pl];
    //      PLOOPNOB2(pl) p_l[pl]=p_r[pl]=pr[pl];


    //////////////////
    //
    // pnew_l -> p_l and pnew_r -> p_r
    //
    // Can use soft switch in time to ramp-up to less dissipative solution
    //
    //////////////////
    // only turn this on once sure no shocks or discontinuities formed near NS surface, else leads to instabilities near surface since not enough dissipation for shocks/discontinuities to spread-out
    FTYPE switchvalue=switcht(fluxtime,SWITCHT2,SWITCHT4);
    //    switchvalue=1.0; // override TEST
    PLOOP(pliter,pl){

      // DEBUG:
      //      if(faceinside||faceshell) dualfprintf(fail_file,"PLPR: %d %d %d :: pl=%d : new: %21.15g oldl=%21.15g oldr=%21.15g ns=%21.15g\n",i,j,k,pl,pr[pl],p_l[pl],p_r[pl],prother[pl]);

      p_l[pl] = p_l[pl] * (1.0-switchvalue) + pnew_l[pl]*switchvalue ;
      p_r[pl] = p_r[pl] * (1.0-switchvalue) + pnew_r[pl]*switchvalue ;
    }



  }// done if flux is on boundary of NS

}




// initialize AND evolve NS-BH setup
// inittype == -3 : evolve and set regardless of whether inside NS or not, but use input B^i orthogonal to EMF_{dir} located at loc=pos=CORN_{dir} and otherwise use internal interpolation.
// inittype == -2 : evolve and set regardless of whether inside NS or not, but use input prim for field for interpolation instead of poor interpolation inside init_dsandvels_nsbh().    Use pstag as external to NS values -- but note that already assume copied over Bperp1 and Bperp2 from exterior to NS into pr from what is inside pstag.
// inittype == -1 : evolve and set regardless of whether inside NS or not
// inittype == 0 : evolve (currently only centered BC's)
// inittype == 1 : field not set yet
// inittype == 2 : field set, so can setup stationary conditions
// also currently uses pglobal and pstagglobal because need multiple points to average to correct position
//
// returns: -1 : no change
//           0 : good and keep change
//           1 : bad - stop code
//
// pr is filled with loc=pos if possible (or loc=CENT for inittype not -2 or -3) but for values internal to NS (if could choose between internal and external -- p_l and p_r values)
// pstag is filled with loc=pos value of quantities for values that are external to the NS
//
// Both pr and pstag are PRIMECOORDS/WHICHVEL values
//
int init_dsandvels_nsbh(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE V[NDIM],r,th,ph;
  FTYPE dxdxp[NDIM][NDIM];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_geom geomdontusebl;
  struct of_geom *ptrgeombl=&geomdontusebl;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  int pl,pliter;
  FTYPE localpr[NPR];
  int insideNS;
  int jj;


  dualfprintf(fail_file,"init_dsandvels_nsbh: nstep=%ld steppart=%d time=%21.15g : ijk=%d %d %d\n",nstep,steppart,time,i,j,k);


  if(FLUXB!=FLUXCTSTAG){
    dualfprintf(fail_file,"nsbh not setup for Toth\n");
    myexit(21359283);
  }
  

  // get coordinate stuff
  bl_coord_ijk(i,j,k,pos,V);
  V[0]=time; // override with desired time in case later things use V[0]
  dxdxprim_ijk(i,j,k,pos,dxdxp);
  r=V[1];
  th=V[2];
  ph=V[3];






  // check if inside or outside NS (based upon pos=CENT of cell)
  insideNS=GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE);
  //  dualfprintf(fail_file,"Doing: inittype=%d : insideNS=%d : %d %d %d\n",inittype,insideNS,i,j,k);




  //  dualfprintf(fail_file,"init_dsandvels_nsbh(): %d %d : %d : %d %d %d\n",insideNS,inittype,pos,i,j,k);



  if(insideNS==0 && (inittype>0) || inittype==1){ // outside NS and setting up atmospheree
    // if inittype==1, then field not set yet, so don't set anything (i.e. leave as atmosphere)

    get_geometry(i, j, k, pos, ptrgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system
  

    // initialize field, assume field set later
    if(FLUXB==FLUXCTSTAG && pstag!=NULL){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    //    dualfprintf(fail_file,"Do 1: inittype=%d : insideNS=%d : %d %d %d\n",inittype,insideNS,i,j,k);
    return(0);
  }
  else if(insideNS==1 && inittype!=1 || inittype<=-1){ // inside NS (or forced assumed inside NS), so need to set rho,u,v^i as required for stationary B^n, etc.  Also used to evolve these quantities


    dualfprintf(fail_file,"INSET: nstep=%ld steppart=%d t=%21.15g dt=%21.15g : time=%21.15g\n",nstep,steppart,t,dt,time);

    FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
    setNSparms(time,&rns,&omegak, &omegaf, &Rns,&Rsoft,&v0);
    FTYPE v0touse=v0; // default


    // modify v0 so that don't mass-load lines that don't need the mass, such as closed field lines.
    // assumes those field lines accumulate mass, so sigma drops, so just detect changes in sigma across NS surface
    // or just vary density injected, and keep v0 fixed.
    //    v0touse=v0*(sigmaout/sigmain);


    //    dualfprintf(fail_file,"GOD: %21.15g %21.15g %21.15g\n",rns,omegak,omegaf,Rns);

    ///////////////
    //
    // BCs:
    //
    // if inside NS, then:
    //
    // Bpar = B^n = Bdipole(vec{rns}) (using A_\mu using init_vpot_user() )
    // Bperp1 = copy from exterior
    // Bperp2 = copy from exterior
    //
    // vpar (along poloidal line?) \sim 0.5 (not relevant for FF)

    // GODMARK: really only want v0>0 if open field lines, so can use (say) presence of non-zero J.B to determine that.
    // so v0\sim 0.2 J.B/(|J|*|B|+SMALL)
    // too large, and how set v^i can be bad since not boosting properly

    // vperp1 & vperp2 : stat-axi using Omegaorbit and OmegaNS
    //
    // rho0 from (B^n)^2(\theta=0 for dipole)/(rho0 c^2) = 100 (not relevant in FF)
    // u : copy from exterior (not relevant in FF)
    //
    // EMF:
    //
    // EMFperp1 = EMFperp2 = 0 (so d_t(B^n)=0 )  *ON* boundary AND inside (must force since v interpolate to boundary where EMFperp1/2 exist)
    // EMFpar (set via v and B automatically from BCs if set correctly here)



    ///////////////
    //
    // init localpr -- assumes loc=CENT (or close to loc=pos) at this point, changed (by averaging) below OR by choosing pstag that (based upon inittype==-3 or -2) can be at loc=pos
    //
    ///////////////
    PLOOP(pliter,pl) localpr[pl]=pr[pl];
    


    ////////////////////
    //
    // GET FIELD at correct location (USES GLOBAL pglobal or pstagglobal!).  Don't have locally passed pstag as BCs setup.
    //
    ////////////////////

    FTYPE count;
    FTYPE weight,totalweight;
    FTYPE interppr[NPR]; // temporary pr


    PLOOP(pliter,pl) interppr[pl]=0.0; // so can cumulate it below

    //////////////////////////////
    //
    // use staggered or centered fields to get other positions for fields
    // Fill localpr with  new field values
    // These fields are PRIMECOORDS values
    //
    //////////////////////////////
    if(pos==CENT){
      // no change
    }
    else if(pos==FACE1){

      if(inittype==-2){
        // then localpr[B1,B2,B3] are already interpolated well by code as given by pr[B1,B2,B3] and already copied
        // But for B2,B3 assume already filled pr (and so localpr) with values external to NS

        // DEBUG:
        // dualfprintf(fail_file,"Got 0: %21.15g %21.15g %21.15g\n",localpr[B1],localpr[B2],localpr[B3]);
      }
      else{
        localpr[B1]=GLOBALMACP0A1(pstagglobal,i,j,k,B1);

        // average B2 center to get FACE1
        // average B3 center to get FACE1
        if(i<=N1M-1-N1NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i,j,k,B2)*weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i,j,k,B3)*weight;


          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j,k,B2)*weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j,k,B3)*weight;

          if(totalweight>=0.5){
            localpr[B2]=interppr[B2]/totalweight;
            localpr[B3]=interppr[B3]/totalweight;
          }
        }
      }// done if inittype!=-2
    }
    else if(pos==FACE2){

      if(inittype==-2){
        // then interppr[B1,B2,B3] are already interpolated well by code
        // But for B1,B3 assume already filled pr (and so localpr) with values external to NS
      }
      else{
        localpr[B2]=GLOBALMACP0A1(pstagglobal,i,j,k,B2);

        // average B1 center to get FACE2
        // average B3 center to get FACE2
        if(j<=N2M-1-N2NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j,k,B1)*weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i,j,k,B3)*weight;


          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j-N2NOT1,k,B1)*weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i,j-N2NOT1,k,B3)*weight;

          if(totalweight>=0.5){
            localpr[B1]=interppr[B1]/totalweight;
            localpr[B3]=interppr[B3]/totalweight;
          }
        }
      }// done if inittype!=-2



    }
    else if(pos==FACE3){

      if(inittype==-2){
        // then interppr[B1,B2,B3] are already interpolated well by code
        // But for B1,B2 assume already filled pr (and so localpr) with values external to NS
      }
      else{
        localpr[B3]=GLOBALMACP0A1(pstagglobal,i,j,k,B3);

        // average B1 center to get FACE3
        // average B2 center to get FACE3
        if(k<=N3M-1-N3NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j,k,B1)*weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i,j,k,B2)*weight;


          weight=(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j,k-N3NOT1,B1)*weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i,j,k-N3NOT1,B2)*weight;

          if(totalweight>=0.5){
            localpr[B1]=interppr[B1]/totalweight;
            localpr[B2]=interppr[B2]/totalweight;
          }
        }
      }// done if inittype!=-2

    }
    else if(pos==CORN1){

      // if inittype==-3, B1,B2,B3 are already set properly by _emf function (B1 using averaging unless Nvec[1]==1)

      if(inittype==-2 || inittype==-3){
        // then interppr[B1,B2,B3] are already interpolated well by code
      }
      else{

        // average B1 center to get CORN1
        if(k>=-N3BND+N3NOT1 && j>=-N2BND+N2NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j,k,B1)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j,k-N3NOT1,B1)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j-N2NOT1,k,B1)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pglobal,i,j-N2NOT1,k-N3NOT1,B1)*weight;
   
          if(totalweight>=0.5) localpr[B1]=interppr[B1]/totalweight;
        }

        // average B2 face2 in interpdir=3 to get CORN1
        if(k>=-N3BND+N3NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pstagglobal,i,j,k,B2)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pstagglobal,i,j,k-N3NOT1,B2)*weight;
   
          if(totalweight>=0.5){
            interppr[B2]/=totalweight;
            localpr[B2]=interppr[B2];
          }
        }
        // average B3 face3 in interpdir=2 to get CORN1
        if(j>=-N2BND+N2NOT1 && inittype!=-3){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pstagglobal,i,j,k,B3)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pstagglobal,i,j-N2NOT1,k,B3)*weight;

          if(totalweight>=0.5){
            interppr[B3]/=totalweight;
            localpr[B3]=interppr[B3];
          }
        }
      

        // DEBUG:
        if(fabs(pr[B2]-localpr[B2])/(fabs(pr[B2])+fabs(localpr[B2]))>1E-1){
          dualfprintf(fail_file,"CORN1DIFFB2: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B2],localpr[B2]);
        }
        if(fabs(pr[B3]-localpr[B3])/(fabs(pr[B3])+fabs(localpr[B3]))>1E-1){
          dualfprintf(fail_file,"CORN1DIFFB3: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B3],localpr[B3]);
        }

      }// end if interptype!=-2 && !=-3

    }
    else if(pos==CORN2){

      if(inittype==-2 || inittype==-3){
        // then interppr[B1,B2,B3] are already interpolated well by code (or as best possible)
      }
      else{

        // average B2 center to get CORN2
        if(k>=-N3BND+N3NOT1 && i>=-N1BND+N1NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i,j,k,B2)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i,j,k-N3NOT1,B2)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j,k,B2)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j,k-N3NOT1,B2)*weight;
   
          if(totalweight>=0.5) localpr[B2]=interppr[B2]/totalweight;
        }

        // average B1 face1 in interpdir=3 to get CORN2
        if(k>=-N3BND+N3NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pstagglobal,i,j,k,B1)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pstagglobal,i,j,k-N3NOT1,B1)*weight;

          if(totalweight>=0.5){
            interppr[B1]/=totalweight;
            localpr[B1]=interppr[B1];
          }
        }
        // average B3 face3 in interpdir=1 to get CORN2
        if(i>=-N1BND+N1NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pstagglobal,i,j,k,B3)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pstagglobal,i-N1NOT1,j,k,B3)*weight;

          if(totalweight>=0.5){
            interppr[B3]/=totalweight;
            localpr[B3]=interppr[B3];
          }
        }


        // DEBUG:
        if(fabs(pr[B1]-localpr[B1])/(fabs(pr[B1])+fabs(localpr[B1]))>1E-1){
          dualfprintf(fail_file,"CORN2DIFFB1: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B1],localpr[B1]);
        }
        if(fabs(pr[B3]-localpr[B3])/(fabs(pr[B3])+fabs(localpr[B3]))>1E-1){
          dualfprintf(fail_file,"CORN2DIFFB3: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B3],localpr[B3]);
        }
 
      }// end if interptype!=-2 && !=-3

    }
    else if(pos==CORN3){

      if(inittype==-2 || inittype==-3){
        // then interppr[B1,B2,B3] are already interpolated well by code
      }
      else{

        // average B3 center to get CORN3
        if(j>=-N2BND+N2NOT1 && i>=-N1BND+N1NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i,j,k,B3)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i,j-N2NOT1,k,B3)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j,k,B3)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B3]+=GLOBALMACP0A1(pglobal,i-N1NOT1,j-N2NOT1,k,B3)*weight;
   
          if(totalweight>=0.5) localpr[B3]=interppr[B3]/totalweight;
        }

        // average B1 face1 in interpdir=2 to get CORN3
        if(j>=-N2BND+N2NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pstagglobal,i,j,k,B1)*weight;
          weight=(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B1]+=GLOBALMACP0A1(pstagglobal,i,j-N2NOT1,k,B1)*weight;

          if(totalweight>=0.5){
            interppr[B1]/=totalweight;
            localpr[B1]=interppr[B1];
          }
        }
        // average B2 face2 in interpdir=1 to get CORN3
        if(i>=-N1BND+N1NOT1){
          totalweight=0.0;    
          weight=(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pstagglobal,i,j,k,B2)*weight;
          weight=(GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1);
          totalweight+=weight;
          interppr[B2]+=GLOBALMACP0A1(pstagglobal,i-N1NOT1,j,k,B2)*weight;
   
          if(totalweight>=0.5){
            interppr[B2]/=totalweight;
            localpr[B2]=interppr[B2];
          }
        }

        // DEBUG:
        if(fabs(pr[B1]-localpr[B1])/(fabs(pr[B1])+fabs(localpr[B1]))>1E-1){
          dualfprintf(fail_file,"CORN3DIFFB1: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B1],localpr[B1]);
        }
        if(fabs(pr[B2]-localpr[B2])/(fabs(pr[B2])+fabs(localpr[B2]))>1E-1){
          dualfprintf(fail_file,"CORN3DIFFB2: %d %d %d : %21.15g %21.15g\n",i,j,k,pr[B2],localpr[B2]);
        }


      }// end if interptype!=-2 && !=-3

    }
    else if(pos==CORNT){ // GODMARK: not controlled with mask yet, but not needed yet

      if(inittype==-2){
        // then interppr[B1,B2,B3] are already interpolated well by code
      }
      else{
        if(k>=-N3BND+N3NOT1 && j>=-N2BND+N2NOT1){
          interppr[B1]=0.25*(GLOBALMACP0A1(pstagglobal,i,j,k,B1)+GLOBALMACP0A1(pstagglobal,i,j-N2NOT1,k,B1)+GLOBALMACP0A1(pstagglobal,i,j,k-N3NOT1,B1)+GLOBALMACP0A1(pstagglobal,i,j-N2NOT1,k-N3NOT1,B1)); // average B1 face1 to get CORNT
        }
        if(i>=-N1BND+N1NOT1 && k>=-N3BND+N3NOT1){
          interppr[B2]=0.25*(GLOBALMACP0A1(pstagglobal,i,j,k,B2)+GLOBALMACP0A1(pstagglobal,i-N1NOT1,j,k,B2)+GLOBALMACP0A1(pstagglobal,i,j,k-N3NOT1,B2)+GLOBALMACP0A1(pstagglobal,i-N1NOT1,j,k-N3NOT1,B2)); // average B2 face2 to get CORNT
        }
        if(i>=-N1BND+N1NOT1 && j>=-N2BND+N2NOT1){
          interppr[B3]=0.25*(GLOBALMACP0A1(pstagglobal,i,j,k,B3)+GLOBALMACP0A1(pstagglobal,i,j-N2NOT1,k,B3)+GLOBALMACP0A1(pstagglobal,i-N1NOT1,j,k,B3)+GLOBALMACP0A1(pstagglobal,i-N1NOT1,j-N2NOT1,k,B3)); // average B3 face3 to get CORNT
        }
      }
    }



    ////////////////////
    //
    // copy field back in (in case changed)
    //
    ///////////////////
    PLOOP(pliter,pl) pr[pl]=localpr[pl];



    ////////////////////
    //
    // set density
    //
    ///////////////////
    FTYPE sigmaswitch;
    FTYPE sigma;

    // 0->1
    sigmaswitch=switcht(time,SWITCHT0,SWITCHT2);
    // set sigma with (in general) a temporal switch
    sigma=SIGMA0*(sigmaswitch) + SIGMAT0*(1.0-sigmaswitch);
    // GODMARK TODO:
    // can also compare local interior (within NS) sigma (i.e. using local bsq instead of BPOLE*BPOLE) and exterior (to NS) sigma
    // Then increase sigma if exterior is too low, and decrease sigma if exterior is too high.
    // All the while,  keeping v0 fixed
      

    //    dualfprintf(fail_file,"omegafswitch=%21.15g omegaf=%21.15g : sigmaswitch=%21.15g sigma=%21.15g\n",omegafswitch,omegaf,sigmaswitch,sigma);


    ////////////////
    // things to set
    ////////////////
    // set \mu or \sigma_0 approximately
    pr[RHO]=BPOLE*BPOLE/(2.0*sigma); // set density for all inittype==0,1,2


    if((inittype==-2 || inittype==0 && inittypeglobal==2) && pstag!=NULL){
      // when inittype==-2, pstag holds other side of primitive value used for flux calculation, which would correspond to the exterior of the NS

      // accurate (i.e. using correct pos for B^i) local b^2
      FTYPE bsq_int;
      get_geometry(i, j, k, pos, ptrgeom);
      if (bsq_calc(pr, ptrgeom, &bsq_int) >= 1) FAILSTATEMENT("init.c:init_dsandvels_nsbh()", "bsq_calc()", 1);
      FTYPE sigma_int;
      sigma_int=bsq_int/(2.0*fabs(pr[RHO])+SMALL);

      // external b^2
      FTYPE bsq_ext;
      if (bsq_calc(pstag, ptrgeom, &bsq_ext) >= 1) FAILSTATEMENT("init.c:init_dsandvels_nsbh()", "bsq_calc()", 2);
      FTYPE sigma_ext;
      sigma_ext=bsq_ext/(2.0*fabs(pstag[RHO])+SMALL);


      FTYPE modinflow = (fabs(sigma_ext)/(fabs(sigma_int)+SMALL)); // avoid nan-out when on inittype==0 && inittypeglobal==2 but still just bounding before 


      if(EOMTYPE!=EOMFFDE){
#if(0)
        // if external sigma is low/high, feed in less/more material so environment enters an equilibrium
        if(sigma_ext<sigma_int){
          pr[RHO]*=modinflow;
          // still won't be stationary until infinite density beyond NS
        }
#else
        // instead, mod v0
        // force no injection of mass from NS if sigma low enough exterior to NS
        // drop down v0 if exterior sigma is much lower than desired
        v0touse=v0*switcht(modinflow,0.2,1.0);
#endif
      }
      else{
        v0touse=0.0;
      }


      // GODMARK: TODO: WHY THE HECK does this get reached with B=0 (RHO=0 too?) even though inittypeglobal==2?  Need to figure this out!  Different variable than global common one used?


      // DEBUG:
      //      dualfprintf(fail_file,"sigmafix: inittype=%d itg=%d : ijk = %d %d %d :: %21.15g %21.15g\n",inittype,inittypeglobal,i,j,k,sigma_int,sigma_ext);

    }

    // other things to set
    // EMFperp1=EMFperp2=0 on boundary and inside (so B^normal remains constant).  But why B^normal for funny surface?  Cause that's what I can control to machine error to avoid drift of field lines.


    ////////////////
    // other things also copied in
    ////////////////
    //    if(inittype==0 || inittype==2){
    // set set them initially
    // pr[UU] is FIXED if super-slow flow, and generally expect NS to be cold so maybe that's good.
    // otherwise, have to extrapolate from active domain
    pr[UU]=BPOLE*BPOLE/(2.0*HOTSIGMA(sigma)); // so cold -- will copy u into ghost cells and overwrite this in end


    // assume Bperp1 and Bperp2 don't have to be set yet

    // asssume EMFpar copied in

    // GODMARK: TODO: eventually should outflow/copy/extrapolate pr[UU]







    // DEBUG
    //    SLOOPA(jj) dualfprintf(fail_file,"1Bcon[KSP,%d]=%21.15g\n",jj,localpr[B1+jj-1]);


    // convert PRIMECOORDS B^i to BL-coord version (should be ok since outside horizon inside NS)
    // using localpr (oringally at CENT unless inittype==-2,-3 and then field is at loc=pos) might be problem at loc==pos -- generally should only transform field!  GODMARK
    metp2met2bl_genloc(WHICHVEL,BLCOORDS,localpr,i,j,k,pos);

    // get BL geom so can lower B^i[BL]
    gset_genloc(0,BLCOORDS,i,j,k,pos,ptrgeombl);


    // Bcon [BL]
    FTYPE Bcon[NDIM];
    // Fill Bcon[BL]
    Bcon[TT]=0;
    SLOOPA(jj) Bcon[jj]=localpr[B1+jj-1];

    // DEBUG
    //    SLOOPA(jj) dualfprintf(fail_file,"1Bcon[BL,%d]=%21.15g\n",jj,Bcon[jj]);

    // set this using B^n, but what about time component?  Bcon.normalveccov=|Bcon||normalveccov|?  orthonormal or what?
    // want Btotal^i n_i = (Bn^i Bn_i)/|Bn|  so that v^i = v0 B^i/|Bn| = ??
    //      FTYPE Bccov[NDIM];
    //      // get Bccov
    //      lower_vec(Bccon,ptrgeom,Bccov);
    //      FTYPE normalvec[NDIM];
    // fill normalvec
    //      normalvec[TT]=0;
    //      SLOOPA(jj) normalvec[jj]=Bcon[jj]/
    // get stationary v^i (hard to decide what normal B is, and doesn't make much sense for racketting surface -- B_p or total B make more sense)
    //      OBtopr_general3n(omegaf,v0,Bcon,normalvec,ptrgeom,localpr);


    // get stationary v^i (just along full B)
    // obtains pr[U1,U2,U3] in PRIMECOORDS
    //      OBtopr_general3(omegaf,v0,Bcon,ptrgeom,localpr);

    // get Bcov[BL]
    FTYPE Bcov[NDIM];
    lower_vec(Bcon,ptrgeombl,Bcov);
    FTYPE Bsq=0.0+SMALL;
    // get B^2 and |B| in BL
    DLOOPA(jj) Bsq+=Bcon[jj]*Bcov[jj];
    FTYPE absB=sqrt(fabs(Bsq));



    // get fluid (not EM) v^i for stationary solution in 3-velocity BL-coords
    FTYPE vconperpB[NDIM],vconparB[NDIM],vcon[NDIM];



    /////////////
    //
    // force corotation of field with NS
    //
    /////////////



    // how fluid moves on the NS *surface*
    FTYPE vconperpBfluid[NDIM];
    vconperpBfluid[1] = 0.0;
    vconperpBfluid[2] = 0.0;
    vconperpBfluid[3] = omegaf;


    if(EOMTYPE==EOMGRMHD||EOMTYPE==EOMENTROPYGRMHD||EOMTYPE==EOMCOLDGRMHD){

      vconperpB[0] = 0.0;
      vconperpB[1] = vconperpBfluid[1];
      vconperpB[2] = vconperpBfluid[2];
      vconperpB[3] = vconperpBfluid[3];

      // motion along field line
      vconparB[0] = 0.0;
      vconparB[1] = v0touse*Bcon[1]/absB;
      vconparB[2] = v0touse*Bcon[2]/absB;
      vconparB[3] = v0touse*Bcon[3]/absB;

      // setting vtot = vfluid
      DLOOPA(jj) vcon[jj] = vconparB[jj] + vconperpB[jj];
    }
    else{// EOMTYPE==EOMFFDE
      //      vconperpB[1] =        - Bcon[1]/(absB*absB)*(Bcov[TT] + omegaf * Bcov[PH]);
      //      vconperpB[2] =        - Bcon[2]/(absB*absB)*(Bcov[TT] + omegaf * Bcov[PH]);
      //      vconperpB[3] = omegaf - Bcon[3]/(absB*absB)*(Bcov[TT] + omegaf * Bcov[PH]);

      // written more generally and such that no catastrophic cancellation if (e.g.) toroidally-dominated (i.e. removed 1-B_\phi B^\phi/B^2)
      vconperpB[0] = 0.0;
      vconperpB[1] = (vconperpBfluid[1]*(Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3]) - vconperpBfluid[2]*Bcov[2]*Bcon[1] - vconperpBfluid[3]*Bcov[3]*Bcon[1] - Bcon[1]*Bcov[TT])/Bsq;
      vconperpB[2] = (vconperpBfluid[2]*(Bcon[1]*Bcov[1] + Bcon[3]*Bcov[3]) - vconperpBfluid[1]*Bcov[1]*Bcon[2] - vconperpBfluid[3]*Bcov[3]*Bcon[2] - Bcon[2]*Bcov[TT])/Bsq;
      vconperpB[3] = (vconperpBfluid[3]*(Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2]) - vconperpBfluid[1]*Bcov[1]*Bcon[3] - vconperpBfluid[2]*Bcov[2]*Bcon[3] - Bcon[3]*Bcov[TT])/Bsq;

      // then parallel component v0touse doesn't enter the equations, and need to set v to vEM instead of vtot=vEM+vparB
      vconparB[0] = 0.0;
      vconparB[1] = 0.0;
      vconparB[2] = 0.0;
      vconparB[3] = 0.0;

      // setting vtot = vperpB
      DLOOPA(jj) vcon[jj] = vconparB[jj] + vconperpB[jj];

    }



    // old
    // pick sign of v^i so that pointing OUT of NS and doesn't depend upon sign of B^i
    // assumes spherical polar coordinates for V[]
    //    FTYPE t=V[TT],r=V[RR],th=V[TH],ph=V[PH];
    //    FTYPE R=r*sin(th);
    //    FTYPE x,y,z;
    // Cartesian position
    // assumes V is spherical polar coordinates
    //    x = r*cos(ph)*sin(th);
    //    y = r*sin(ph)*sin(th);
    //    z = r*cos(th);


    // Check signature of vcon[BL] using Bcon[BL]
    FTYPE inout;
    get_insideNS(time, i, j, k, V, Bcon, vconparB, &inout);

    // correct signature of vcon along B so fluid flow is always going *out* of NS
    DLOOPA(jj) vcon[jj]=vconparB[jj]*inout + vconperpB[jj];



    // DEBUG:
    FTYPE vphieff=omegaf*(r*sin(th));
    dualfprintf(fail_file,"i=%d j=%d k=%d : vphieff=%21.15g v0touse=%21.15g\n",i,j,k,vphieff,v0touse);
    SLOOPA(jj) dualfprintf(fail_file,"vconbl[%d]=%21.15g\n",jj,vcon[jj]);



    ////////////
    //
    // Before final conversion back to PRIMECOORDS/WHICHVEL, must ensure 3-velocity is good.  For force-free, easiest to pass through inversion as filter to limit velocity
    //
    ////////////
    if(EOMTYPE==EOMFFDE){

      // pr will contain WHICHVEL velocity, but velocity and field are in BL coords
      limit_3vel_ffde(Bcon,ptrgeombl,vcon,pr);


      // DEBUG:
      SLOOPA(jj) dualfprintf(fail_file,"u[WHICHVEL,BL][%d]=%21.15g\n",jj,pr[U1+jj-1]);


      // pr will obtain correct pr[RHO,UU] below
      // convert WHICHVEL utilde^i and 3-B B^i from BL coords to WHICHVEL/PRIMECOORDS at pos position
      bl2met2metp2v_genloc(WHICHVEL,BLCOORDS,pr,i,j,k,pos);

    }
    else{
      // fill pr with v^i [BL]
      SLOOPA(jj) pr[U1+jj-1]=vcon[jj];
      // fill pr with B[BL] since assume vel and B both use same coordinates
      SLOOPA(jj) pr[B1+jj-1]=Bcon[jj];

      // pr will obtain correct pr[RHO,UU] below
      // convert 3-vel v^i and 3-B B^i from BL coords to WHICHVEL/PRIMECOORDS at pos position
      bl2met2metp2v_genloc(VEL3,BLCOORDS,pr,i,j,k,pos);
    }



    // DEBUG:
    //    SLOOPA(jj) dualfprintf(fail_file,"2Bcon[BL,%d]=%21.15g\n",jj,Bcon[jj]);



    // DEBUG
    //    SLOOPA(jj) dualfprintf(fail_file,"2Bcon[KSP,%d]=%21.15g\n",jj,pr[B1+jj-1]);

    // now all is in WHICHVEL/PRIMECOORDS, so no need to convert anything.  Do the above conversion here so don't convert pstag.


    // DEBUG:
    //    SLOOPA(jj) dualfprintf(fail_file,"vconprimecoords[%d]=%21.15g\n",jj,pr[U1+jj-1]);




    //      *whichvel=VEL3;
    //      *whichcoord=BLCOORDS;
    //      return(0);
      















    // if reach here, then didn't  use BL-coords and velocity and field are same as before
    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    //    dualfprintf(fail_file,"Do 2: inittype=%d : insideNS=%d : %d %d %d\n",inittype,insideNS,i,j,k);
    return(0);
  }
  else{

    // no change
    // Don't modify prim or pstag so keeps original value (such as for outside NS during evolution, since then no boundary condition to apply)


    // default is do nothing with pr
    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    //    dualfprintf(fail_file,"Do nothing: inittype=%d : insideNS=%d : %d %d %d\n",insideNS,inittype,i,j,k);
    return(-1);
  }


}










#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3
#define TOROIDALFIELD 4
#define NSFIELDOFFSET1 5


//#define FIELDTYPE TOROIDALFIELD
//#define FIELDTYPE DISKFIELD
#define FIELDTYPE NSFIELDOFFSET1


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



int setNSparms(SFTYPE time, FTYPE *rns, FTYPE *omegak, FTYPE *omegaf, FTYPE *Rns, FTYPE *Rsoft, FTYPE *v0)
{

  *rns=RSEP; // orbital separation
  *omegak=1.0/(pow( (*rns) ,3.0/2.0)+a);
  *Rns=2.0; // as for similar size of non-rotating BH
  *Rsoft=(*Rns) *0.1; // soften potential so like interior solution for dipole
  *v0=0.2;


  // TEST
  //  *omegak=0;

  // DEBUG:
  //  dualfprintf(fail_file,"setNSparms: setNSparms: time=%21.15g\n",time);

  
  // 0->1 after sigma has increased somewhat
  FTYPE omegafswitch=switcht(time,SWITCHT4,SWITCHT6);
  // TEST
  //  omegafswitch=1.0;
  
  *omegaf=(*omegak)*(omegafswitch) + 0.0*(1.0-omegafswitch);



  return(0);
}










// offset NS dipole
FTYPE setnsoffset1(int l, SFTYPE time, FTYPE *V)
{
  FTYPE A3sign;
  FTYPE aphi;
  V[0]=time; // override
  FTYPE r=V[RR],th=V[TH],ph=V[PH],R=r*sin(th); // time is for time variable
  FTYPE x,y,z;
  FTYPE xns,yns,zns;
  FTYPE Vns[NDIM],Vcartns[NDIM];


  // GODMARK TODO: Need to boost A_\mu in general, unless motion not very relativistic


  // Cartesian position
  x = r*cos(ph)*sin(th);
  y = r*sin(ph)*sin(th);
  z = r*cos(th);


  // get NS position
  FTYPE absrdiff;
  pos_NS(time, V, Vns, Vcartns, &absrdiff);
  xns=Vcartns[1];
  yns=Vcartns[2];
  zns=Vcartns[3];

  
  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  setNSparms(time, &rns,&omegak, &omegaf, &Rns,&Rsoft,&v0);



  // initialize aphi
  aphi=0.0;


  if(l==3){
    
#if(FLIPGDETAXIS==0)
    // causes kink at pole in B^r and B^\theta with flip of gdet sign
    if(th<0.0) A3sign=-1.0;
    else if(th>M_PI) A3sign=-1.0;
    else A3sign=1.0;
#else
    A3sign=1.0; // to be used with gdet flip of sign
#endif



    if(NSSHAPE==SPHERICAL){
      /////////////////
      //
      // spherical NS
      //
      /////////////////

      // Use Rin so same size as BH (as if MBH=10Msun and MNS=1.4Msun)
      aphi += normglobal*A3sign*0.5*BPOLE*Rns*Rns*Rns*(xns-x)*pow(xns-x,1.0)/pow(absrdiff+Rsoft,3.0) ;

      // GODMARK: Have to ensure sign makes sense across pole, etc.

    }// end if spherical NS
    else if(NSSHAPE==DONUT){
      /////////////////
      //
      // donut NS
      //
      /////////////////
      FTYPE aphidipole;

      aphidipole=normglobal*A3sign*0.5*BPOLE*Rns*Rns*Rns*(rns-R)*pow(fabs(rns-R),1.0)/pow(absrdiff+Rsoft,3.0) ;

      if(1){ // too out of balance at t=0
        // Damp potential near R=0 since otherwise (in axisymmetry) the field must have a monopole if B^\theta finite near \theta=0
        // That is, if field points into pole, it also points into pole on other side of pole.  But then don't just have kink in B, have monopole.
 
        // So aphi->0 on order of Rin or as \theta->0
        // aphidipole*=1.0/(1.0+1.0/(R/Rns));

        // So aphi->0 on order of rns or as \theta->0
        aphidipole*=1.0/(1.0+1.0/(R/rns));
      }
      if(0){ // makes no sense in axisymmetry
        // force current sheet at pole so kink in field is not a monopole
        if(th<0 || th>M_PI) aphidipole*=-1.0;
      }

      aphi += aphidipole;





    }// end if spherical NS





  } // end l==3




  return(aphi);


}










// assumes normal field in pr
// SUPERNOTE: A_i must be computed consistently across all CPUs.  So, for example, cannot use randomization of vector potential here.
int init_vpot_user(int *whichcoord, int l, SFTYPE time, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av,q;
  FTYPE r,th;
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




  if(FIELDTYPE==NSFIELDOFFSET1){
    vpot += setnsoffset1(l,time,V);
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
  void adjust_fluxctstag_vpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  int Nvec[NDIM];


  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  if(FLUXB==FLUXCTSTAG){
    // before ccompute anything, fix-up A_i as done during evolution after A_i is updated
    adjust_fluxctstag_vpot(time, prim, Nvec, A);
  }
  else{
    dualfprintf(fail_file,"SUPERWARNING: Should use staggered field for NSBH problem\n");
    myexit(568234642);
  }


  // now compute (prim/pstag/ucons/Bhat)[A_i]
  funreturn=user1_init_vpot2field_user(time, fieldfrompotential, A, prim, pstag, ucons, Bhat); // time ultimately needed for a boundary call
  if(funreturn!=0) return(funreturn);
 

  return(0);


}



// assumes we are fed the true densities
int normalize_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  FTYPE parms[MAXPASSPARMS];
  int eqline;


  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK){
    eqline=1;
    parms[0]=rin;
    parms[1]=rhodisk;
    
    funreturn=user1_normalize_densities(eqline, parms, prim, &rhomax, &umax);
    if(funreturn!=0) return(funreturn);
  }
  else if(WHICHPROBLEM==NSBH){

    // done when normalize field
  }

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
  int normalize_field_nsbh(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
  int normalize_densities_postnormalizefield(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

 
  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK){
    funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
    if(funreturn!=0) return(funreturn);
  }
  else if(WHICHPROBLEM==NSBH){
    funreturn=normalize_field_nsbh(t, prim, pstag, ucons, vpot, Bhat); // t is ok here
    if(EOMTYPE!=EOMFFDE) normalize_densities_postnormalizefield(t, prim); // t is ok here

    if(funreturn!=0) return(funreturn);    
  }
 
  return(0);

}


#if(EOMTYPE==EOMGRMHD||EOMTYPE==EOMCOLDGRMHD||EOMTYPE==EOMENTROPYGRMHD)
#define NORMALIZEFIELDMETHOD 1 // choice
// 0 : by sigma (not best for varying density with time)
// 1 : by b^2\sim B^2 (best for varying density with time)
#else
#define NORMALIZEFIELDMETHOD 1 // no choice
#endif

// assumes normal field definition for NSBH problem
int normalize_field_nsbh(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{


  FTYPE sigma_pole, bsq_pole, norm;
  int get_sigmabsq_atpole(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *sigma_pole, FTYPE *bsq_pole);
  
  // get initial maximum
  get_sigmabsq_atpole(time, prim, &sigma_pole,&bsq_pole);
  trifprintf("initial sigma_pole: %21.15g bsq_pole: %21.15g\n", sigma_pole,bsq_pole);

  if(NORMALIZEFIELDMETHOD==0){
    // get normalization parameter (uses final sigma to set field strength)
    norm = sqrt(SIGMA0/sigma_pole);
    trifprintf("initial sigma_pole: %21.15g (should be %21.15g) norm=%21.15g\n", sigma_pole,SIGMA0,norm);
  }
  else if(NORMALIZEFIELDMETHOD==1){
    // get normalization parameter (uses final sigma to set field strength)
    norm = sqrt(BPOLE*BPOLE/bsq_pole);
    trifprintf("initial bsq_pole: %21.15g (should be %21.15g) norm=%21.15g\n", bsq_pole,BPOLE*BPOLE,norm);
  }
  
  
  // not quite right since only correct static field energy, not moving field energy
  normalize_field_withnorm(norm, prim, pstag, ucons, vpot, Bhat);


  // get new maxes to check if normalization is correct
  get_sigmabsq_atpole(time, prim, &sigma_pole,&bsq_pole);
  trifprintf("new initial sigma_pole: %21.15g bsq_pole: %21.15g\n", sigma_pole,bsq_pole);


  if(NORMALIZEFIELDMETHOD==0){
    trifprintf("initial sigma_pole: %21.15g (should be %21.15g) norm=%21.15g\n", sigma_pole,SIGMA0,norm);
  }
  else if(NORMALIZEFIELDMETHOD==1){
    trifprintf("initial bsq_pole: %21.15g (should be %21.15g) norm=%21.15g\n", bsq_pole,BPOLE*BPOLE,norm);
  }



  // for use when calling init_vpot_user later
  // if called this twice, need both renormalizations as multiplied by each other
  normglobal*=norm;


  return(0);
}



// normalize densities after field has been normalized
int normalize_densities_postnormalizefield(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  int loc;

  
  

  loc=CENT;
  // get NS position
  FTYPE Vns[NDIM],Vcartns[NDIM],absrdiff,xns,yns,zns;
   
  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  setNSparms(time, &rns,&omegak, &omegaf, &Rns,&Rsoft,&v0);

  FTYPE R;
  int isinside;

  FTYPE bsq_ij,sigma_ij;
  FTYPE sigmahot_ij;
  FTYPE uorho;


  // LOOP
  FULLLOOP{
    get_geometry(i, j, k, loc, ptrgeom);

    bl_coord_ijk_2(i, j, k, loc, X, V);
    r=V[1];
    th=V[2];
    R=r*sin(th);

    isinside=get_insideNS(time, i, j, k, V, NULL,NULL,NULL);

    if(!isinside){

      if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

      // \sigma \sim \mu \sim b^2/(2\rho_0)
      sigma_ij=bsq_ij/(2.0*fabs(MACP0A1(prim,i,j,k,RHO)+SMALL));

      if(sigma_ij>BSQORHOLIMIT*0.5){
        dualfprintf(fail_file,"RENORMDEN1: %d %d %d %21.15g %21.15g\n",i,j,k,sigma_ij,MACP0A1(prim,i,j,k,RHO));
        MACP0A1(prim,i,j,k,RHO)*=sigma_ij/(BSQORHOLIMIT*0.5);
      }

      // \sigma \sim \mu \sim b^2/(2u)
      sigmahot_ij=bsq_ij/(2.0*fabs(MACP0A1(prim,i,j,k,UU))+SMALL);
        
      if(sigmahot_ij>BSQOULIMIT*0.5){
        dualfprintf(fail_file,"RENORMDEN2: %d %d %d %21.15g %21.15g\n",i,j,k,sigmahot_ij,MACP0A1(prim,i,j,k,UU));
        MACP0A1(prim,i,j,k,UU)*=sigmahot_ij/(BSQOULIMIT*0.5);
      }

      // u/\rho_0
      uorho=MACP0A1(prim,i,j,k,UU)/(fabs(MACP0A1(prim,i,j,k,RHO))+SMALL);

      if(uorho>UORHOLIMIT){
        dualfprintf(fail_file,"RENORMDEN3: %d %d %d %21.15g %21.15g %21.15g\n",i,j,k,uorho,MACP0A1(prim,i,j,k,RHO),MACP0A1(prim,i,j,k,UU));
        MACP0A1(prim,i,j,k,RHO)*=uorho/UORHOLIMIT;
      }

      // old, stupid, maybe:
      //      MACP0A1(prim,i,j,k,RHO)*=normglobal*normglobal;
      //      MACP0A1(prim,i,j,k,UU)*=normglobal*normglobal;

    }
  }

  return(0);


}






// get \sigma at pole of NS
int get_sigmabsq_atpole(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *sigma_pole, FTYPE *bsq_pole)
{
  int i,j,k;
  FTYPE bsq_ij,sigma_ij;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  int gotnormal;
  FTYPE rin;
  int loc;

  
  

  sigma_pole[0] = SMALL;
  bsq_pole[0] = SMALL;
  gotnormal=0; // to check if ever was in location where wanted to normalize
  loc=CENT;


  // get NS position
  FTYPE Vns[NDIM],Vcartns[NDIM],absrdiff,xns,yns,zns;
    

  FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
  setNSparms(time, &rns,&omegak, &omegaf, &Rns,&Rsoft,&v0);

  FTYPE R;
  int isinside;


  // LOOP
  ZLOOP {
    get_geometry(i, j, k, loc, ptrgeom);

    bl_coord_ijk_2(i, j, k, loc, X, V);
    r=V[1];
    th=V[2];
    R=r*sin(th);

    pos_NS(time, V, Vns, Vcartns, &absrdiff);
    xns=Vcartns[1];
    yns=Vcartns[2];
    zns=Vcartns[3];

    isinside=get_insideNS(time, i, j, k, V, NULL,NULL,NULL);

    //    dualfprintf(fail_file,"Got out: %d %d %d : %21.15g %21.15g %21.15g %21.15g\n",i,j,k, absrdiff,Rns,R,Vcartns[1]);
    
    if(isinside==1 && fabs(absrdiff-Rns)<Rns*0.2 && fabs(R-Vcartns[1])<Rns*0.2){ // then assume near pole of NS

      gotnormal=1;
      if (bsq_calc(MAC(prim,i,j,k), ptrgeom, &bsq_ij) >= 1) FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

      // \sigma \sim \mu \sim b^2/(2\rho_0)
      sigma_ij=bsq_ij/(2.0*fabs(MACP0A1(prim,i,j,k,RHO))+SMALL);

      //      dualfprintf(fail_file,"SIGMA: %d %d %d : sigma_ij=%21.15g\n",i,j,k,sigma_ij);

      // if multiple positions found, keep getting maximum
      if (sigma_ij > sigma_pole[0])      sigma_pole[0] = sigma_ij;
      if (bsq_ij > bsq_pole[0])      bsq_pole[0] = bsq_ij;
    }
  }


  mpiisum(&gotnormal);

  if(gotnormal==0){
    dualfprintf(fail_file,"Never found place to normalize field for NS\n");
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

  mpimax(sigma_pole);
  mpimax(bsq_pole);

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

  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK){
    atmospheretype=1;
  }
  else if(WHICHPROBLEM==GRBJET || WHICHPROBLEM==NSBH){
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


