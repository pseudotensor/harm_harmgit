

/*! \file init.koral.c
     \brief USER initial conditions for Koral/RAD simulations
*/

// issues:

// runnew8: RADTUBE: oscillations and slow.  averageiter=6.7
// runnew9: RADSHADOW : FAILINFOs without mathematica solution and CASE and slow.  prad1 fills-in shadow.  maybe radfixups worked to help that.
// runnew10: RADDBLSHADOW: Odd spot in prad0 in top-right in im8p0s0l0001.r8 .  More fake substructure near left wall, but more narrow shadow beam.
// runnew11: ATMSTATIC: prad0,1 evolve majorly, but that's correct.
// runnew12: RADATM: fine.
// runnew13: RADBEAM2D: ok
// runnew14: RADWALL: Bit noisy.  Probably need to use stronger shock condition for radiation.
// runnew15: RADWAVE: Unsure if right.  too fast.
// runnew16: RADBONDI  : goes crazy, bad in prad.
// runnew17: RADDOT: CASEGEN failures!!! and "Bad inv" problems.  Has stripes along certain directions unlike with MINM and older code.  Maybe PARA shock reduction will help here too. But expanding flow?  Or try 2Dify shock flattener.
// RADNT: IC look odd, but need to look at next step.
// runnew23: raddonut still was running and very bad behavior
// runnew18: RADCYLBEAM? was still running.
// runnew19: RADBEAM2DKSVERT: Goes crazy.
// runnew20: RADCYLBEAMCART : Maybe ok.



#include "decs.h"

static SFTYPE rhomax=0,umax=0,uradmax=0,utotmax=0,pmax=0,pradmax=0,ptotmax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin,rinfield,routfield; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;
static FTYPE nz_func(FTYPE R) ;
static FTYPE taper_func2(FTYPE R,FTYPE rin, FTYPE rpow) ;
static int fieldprim(int whichmethod, int whichinversion, int *whichvel, int*whichcoord, int ii, int jj, int kk, FTYPE *pr);

///////
// WALD STUFF
FTYPE B0WALD; // set later
// DOWALDDEN=0 : Turn off doing WALD
// DOWALDDEN=1 : among other things, set densities as floor-like with below b^2/rho at horizon.  Should also choose FIELDTYPE==FIELDWALD.
// DOWALDDEN=2 : monopole case against disk equator
int DOWALDDEN=0;
int WALDWHICHACOV=3; // which component : -1 is all of them
//FTYPE BSQORHOWALD=50.0; // leads to too large b^2/rho and uu0 cylindrical shock forms at r\sim 2r_g and remains forever (at least at 128x64)
FTYPE BSQORHOWALD=100.0;
FTYPE aforwald;
FTYPE TILTWALD;



//FTYPE thindiskrhopow=-3.0/2.0; // can make steeper like -0.7
//FTYPE thindiskrhopow=-0.2; // closer to NT73
FTYPE thindiskrhopow=-0.5; // closer to thick disk // SUPERMADNEW # -1 so Qmri constant vs. radius.

FTYPE normglobal;
int inittypeglobal; // for bounds to communicate detail of what doing

#define SLOWFAC 1.0  /* reduce u_phi by this amount */
#define MAXPASSPARMS 10

//#define THETAROTMETRIC (0.5*0.7)
//#define USER_THETAROTMETRIC (M_PI*0.25)
#define USER_THETAROTMETRIC (M_PI*0.25) // WALD
#define USER_THETAROTPRIMITIVES (0.0) // probably want to choose 0, so initial conditions are as if no tilt // WALD -> make same as USER_THETAROTMETRIC
  


#define NOFIELD -1
#define DISK1FIELD 0
#define DISK2FIELD 1
#define VERTFIELD 2
#define DISK1VERT 3
#define DISK2VERT 4
#define BLANDFORDQUAD 5
#define TOROIDALFIELD 6
#define OHSUGAFIELD 7
#define MONOPOLAR 8
#define OLEKFIELD 9
#define FIELDJONMAD 10
#define FIELDWALD 11
#define MONOPOLE 12
#define SPLITMONOPOLE 13




// NOTE on units:
// Many things below are so far in code units, not physical units, so they don't need conversion.  This includes:
//
// Rin, Rout, tf, DTdumpgen



int RADBEAM2D_BEAMNO;
int RADBEAM2D_FLATBACKGROUND;
FTYPE RADBEAM2D_RHOAMB;
FTYPE RADBEAM2D_TAMB;
int RADBEAM2D_BLOB;
FTYPE RADBEAM2D_BLOBW;
FTYPE RADBEAM2D_BLOBP;
FTYPE RADBEAM2D_BLOBX;
FTYPE RADBEAM2D_BLOBZ;
FTYPE RADBEAM2D_PAR_D;
FTYPE RADBEAM2D_PAR_E;
int RADBEAM2D_IFBEAM;
FTYPE RADBEAM2D_TLEFT;
FTYPE RADBEAM2D_NLEFT;
FTYPE RADBEAM2D_BEAML;
FTYPE RADBEAM2D_BEAMR;



int RADBEAM2DKSVERT_BEAMNO; // global for bounds.koral.c

FTYPE RADBEAMFLAT_FRATIO;
FTYPE RADBEAMFLAT_ERAD;
FTYPE RADBEAMFLAT_RHO;
FTYPE RADBEAMFLAT_UU;

FTYPE RADATM_MDOTEDD;
FTYPE RADATM_LUMEDD;
int RADATM_THINRADATM;
FTYPE RADATM_FERATIO;
FTYPE RADATM_FRATIO;
FTYPE RADATM_RHOAMB;
FTYPE RADATM_TAMB;

FTYPE RADSHADOW_NLEFT;
FTYPE RADSHADOW_ANGLE;
FTYPE RADSHADOW_TLEFTOTAMB;
FTYPE RADSHADOW_BEAMY;

FTYPE RADDBLSHADOW_NLEFT;
FTYPE RADDBLSHADOW_ANGLE;
FTYPE RADDBLSHADOW_TLEFTOTAMB;
FTYPE RADDBLSHADOW_BEAMY;




int RADWAVE_NWAVE;
int RADWAVE_NUMERO;
int RADWAVE_WAVETYPE;
FTYPE RADWAVE_PP;
FTYPE RADWAVE_CC;
FTYPE RADWAVE_KAPPA;
FTYPE RADWAVE_RHOFAC;
FTYPE RADWAVE_B0;
FTYPE RADWAVE_DRRE;
FTYPE RADWAVE_DRIM;
FTYPE RADWAVE_DVRE;
FTYPE RADWAVE_DVIM;
FTYPE RADWAVE_DV2RE;
FTYPE RADWAVE_DV2IM;
FTYPE RADWAVE_DB2RE;
FTYPE RADWAVE_DB2IM;
FTYPE RADWAVE_DURE;
FTYPE RADWAVE_DUIM;
FTYPE RADWAVE_DERE;
FTYPE RADWAVE_DEIM;
FTYPE RADWAVE_DFRE;
FTYPE RADWAVE_DFIM;
FTYPE RADWAVE_DF2RE;
FTYPE RADWAVE_DF2IM;
FTYPE RADWAVE_OMRE;
FTYPE RADWAVE_OMIM;
FTYPE RADWAVE_DTOUT1;
FTYPE RADWAVE_RHOZERO;
FTYPE RADWAVE_AAA;
FTYPE RADWAVE_KK;
FTYPE RADWAVE_UINT;
FTYPE RADWAVE_TEMP;
FTYPE RADWAVE_ERAD;
FTYPE RADWAVE_VX;
FTYPE RADWAVE_KAPPAES;
FTYPE RADWAVE_ERADFACTOR;
FTYPE RADWAVE_GASFACTOR;



FTYPE RADBONDI_TESTNO;
FTYPE RADBONDI_PRADGAS;
FTYPE RADBONDI_TGAS0;
FTYPE RADBONDI_MDOTPEREDD;
FTYPE RADBONDI_MDOTEDD;
FTYPE RADBONDI_MINX;
FTYPE RADBONDI_MAXX;


FTYPE RADDOT_XDOT;
FTYPE RADDOT_YDOT;
FTYPE RADDOT_ZDOT;
int RADDOT_IDOT;
int RADDOT_JDOT;
int RADDOT_KDOT;
FTYPE RADDOT_FYDOT;
FTYPE RADDOT_LTEFACTOR;
FTYPE RADDOT_URFX;
FTYPE RADDOT_F1;
FTYPE RADDOT_F2;



FTYPE RADNT_MINX;
FTYPE RADNT_MAXX;
FTYPE RADNT_KKK;
FTYPE RADNT_ELL;
FTYPE RADNT_UTPOT;
FTYPE RADNT_RHOATMMIN;
FTYPE RADNT_RHODONUT;
FTYPE RADNT_UINTATMMIN;
FTYPE RADNT_ERADATMMIN;
FTYPE RADNT_DONUTTYPE;
FTYPE RADNT_INFLOWING;
FTYPE RADNT_TGASATMMIN;
FTYPE RADNT_TRADATMMIN;
FTYPE RADNT_ROUT;
FTYPE RADNT_OMSCALE;
FTYPE RADNT_FULLPHI;
FTYPE RADNT_DONUTRADPMAX;
FTYPE RADNT_HOVERR;
FTYPE RADNT_LPOW;

int RADCYLJET_TYPE;
FTYPE RADCYLJET_VRSTAR;
FTYPE RADCYLJET_EHATJET;
FTYPE RADCYLJET_RHOJET;
FTYPE RADCYLJET_TEMPJET;
FTYPE RADCYLJET_UJET;

int RADDONUT_OPTICALLYTHICKTORUS;

static int get_full_rtsolution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom);
static int make_nonrt2rt_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom);
static int donut_analytical_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom, FTYPE *ptptr);
static int process_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom, FTYPE *ptptr);




FTYPE normglobal;

int prepre_init_specific_init(void)
{
  int funreturn;

  funreturn=user1_prepre_init_specific_init();


  // set periodicity in x1,x2,x3 directions

  // fully periodic problems
  if(WHICHPROBLEM==FLATNESS){
    periodicx1=periodicx2=periodicx3=1;
  }
  // if ever only 1D problems
  else if(WHICHPROBLEM==RADTUBE || WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR){
    periodicx1=0;
    periodicx2=periodicx3=1;
  }
  // if ever only 1D problems
  else if(WHICHPROBLEM==RADWAVE){
    periodicx1=1;
  }
  // problems with no necessary symmetry
  else if(WHICHPROBLEM==RADBEAMFLAT || WHICHPROBLEM==RADPULSE3D || WHICHPROBLEM==RADDBLSHADOW || WHICHPROBLEM==RADWALL || WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADDOT || WHICHPROBLEM==RADBEAM2DKSVERT || WHICHPROBLEM==RADCYLBEAMCART){
    periodicx1=periodicx2=periodicx3=0;
  }
  // periodic in x3
  else if(WHICHPROBLEM==RADCYLBEAM){
    periodicx1=periodicx2=0;
    periodicx3=1;
  }
  // periodic in x2
  else if(WHICHPROBLEM==RADSHADOW){
    periodicx1=periodicx3=0;
    periodicx2=1;
  }
  // spherical polar problems:
  else if(WHICHPROBLEM==RADDONUT || WHICHPROBLEM==ATMSTATIC || WHICHPROBLEM==RADATM || WHICHPROBLEM==RADBONDI || WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK){
    periodicx1=periodicx2=0;
    periodicx3=1;
  }
  else if(WHICHPROBLEM==KOMIPROBLEM){
    periodicx1=0;
    periodicx2=periodicx3=1;
  }
  // periodic in x3
  else if(WHICHPROBLEM==RADCYLJET){
    periodicx1=periodicx2=0;
    periodicx3=1;
  }
  // assume spherical polar problems:
  else{
    periodicx1=periodicx2=0;
    periodicx3=1;
  }

  // Also: SET USEROMIO to 0 or 1 in mympi.definit.h (needs to be 0 for TEXTOUTPUT)
  if(PRODUCTION==0||DOWALDDEN!=0){ // for now DOWALDDEN!=0
    binaryoutput=TEXTOUTPUT; // WALDPRODUCTION
    // KRAKEN: comment out above.  And change mympi.definit.h's USEROMIO 0 to 1 for the "choice" version.
  }

  //  binaryoutput=TEXTOUTPUT;



  return(0);

}


int pre_init_specific_init(void)
{
  // defaults
  //  h_over_r=0.3;
  h_over_r=0.1;
  h_over_r_jet=2.0*h_over_r;

  if(WHICHPROBLEM==RADDONUT){
    if(RADNT_DONUTTYPE==DONUTTHINDISK||RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      //      h_over_r=0.02;
      h_over_r=0.1; // SUPERMADNEW
    }
    else{
      h_over_r=0.2;
      h_over_r_jet=2.0*h_over_r;
    }
  }
  



  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;


  if(WHICHPROBLEM==RADCYLJET){
    static int firsttime=1;
    if(firsttime==1){
      firsttime=0;
      int itid;
      for(itid=0;itid<numprocs;itid++){
        if(itid==myid){
          FILE *fstar;
          if((fstar=fopen("star.txt","rt"))==NULL){
            dualfprintf(fail_file,"Couldn't open star.txt, assume values not used.\n");
          }
          else{
            logfprintf("opened star.txt and got contents\n");
            fscanf(fstar,"%d %lf %lf %lf %lf %lf",&RADCYLJET_TYPE,&RADCYLJET_RHOJET,&RADCYLJET_UJET,&RADCYLJET_EHATJET,&RADCYLJET_TEMPJET,&RADCYLJET_VRSTAR);
            fclose(fstar);
          }
        }
#if(USEMPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
      }
    }
  }


  return(0);
}

int set_fieldfrompotential(int *fieldfrompotential)
{
  int pl,pliter;

  // default (assume all fields are from potential)
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;

  // force B3=0 so only using poloidal part of Wald solution.
  //int pl=B3; fieldfrompotential[pl-B1+1]=0;
  
  //In the case of Komissarov's tests, set up the field directly in all cases
  if( WHICHPROBLEM==KOMIPROBLEM || WHICHPROBLEM == RADWAVE){
    PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=0;
  }
 



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



  trifprintf("WHICHPROBLEM: %d\n",WHICHPROBLEM);
  // print out units and some constants
  trifprintf("Constants\n");
  trifprintf("LBAR=%g TBAR=%g VBAR=%g RHOBAR=%g MBAR=%g UBAR=%g TEMPBAR=%g\n",LBAR,TBAR,VBAR,RHOBAR,MBAR,UBAR,TEMPBAR); 
  trifprintf("ARAD_CODE=%26.20g OPACITYBAR=%g KAPPA_ES_CODE(1,1)=%g KAPPA_FF_CODE(1,1,1)=%g KAPPA_BF_CODE(1,1,1)=%g KAPPA_GENFF_CODE(1,1)=%g\n",ARAD_CODE,OPACITYBAR,KAPPA_ES_CODE(1,1),KAPPA_FF_CODE(1,1,1),KAPPA_BF_CODE(1,1,1),KAPPA_GENFF_CODE(1,1,1));
  trifprintf("ARAD_CODE_DEF=%g\n",ARAD_CODE_DEF);
  trifprintf("GAMMAMAXRAD=%g\n",GAMMAMAXRAD);

  trifprintf("MASSCM=%g 1 koral unit = %g harm units (g/cm^3)\n",MASSCM,KORAL2HARMRHO(1.0));

  if(myid==0){
    // 22 things
#define DIMVARLIST GGG,CCCTRUE,MSUNCM,MPERSUN,LBAR,TBAR,VBAR,RHOBAR,MBAR,ENBAR,UBAR,TEMPBAR,ARAD_CODE_DEF,XFACT,YFACT,ZFACT,MUMEAN,MUMEAN,OPACITYBAR,MASSCM,KORAL2HARMRHO(1.0),TEMPMIN
#if(REALTYPE==FLOATYPE || REALTYPE==DOUBLETYPE)
#define DIMTYPELIST "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"
#elif(REALTYPE==LONGDOUBLETYPE)
#define DIMTYPELIST "%26.20Lg %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g\n"
#else
#error "no such type."
#endif

    FILE *dimfile;
    dimfile=fopen("dimensions.txt","wt");
    if(dimfile!=NULL){
      fprintf(dimfile,DIMTYPELIST,DIMVARLIST);
      fclose(dimfile);
    }
    else{
      dualfprintf(fail_file,"Could not open dimensions.txt\n");
      myexit(394582526);
    }
  }


  // override anything set by restart file.
#if(PRODUCTION==0)
  debugfail=2;
#else
  debugfail=0;
#endif

  tf = 200000;



  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  MBH=1.0;

  return(0);

}





int init_global(void)
{
  int pl,pliter;
  int funreturn;

  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);

  init_defcoord(); // just avoids splitting function call, here sets a

  // default
  ARAD_CODE=ARAD_CODE_DEF;


  //  ERADLIMIT=UUMIN; // set same for now
  ERADLIMIT=UUMINLIMIT; // seems fine.

  if(WHICHPROBLEM==RADDONUT){
    ERADLIMIT=1E-20; // choose so smallest value PRAD0 likely to obtain, but no smaller else problems with YFL4 (well, still issues)
  }



  // maximum radiation frame lorentz factor
  //GAMMAMAXRAD=10000.0; // problems with PARA or TIMEORDER=3 for NLEFT=0.99999 with RADBEAM2D, so stick to gammamax=100 in general unless for test.
  // NOTE: Even PARA,TO=3 can handle cartesian beams like RADBEAM2D or RADSHADOW without problems, but if injection gamma > GAMMAMAXRAD, then that limiting process causes problems currently.  Looking into it.
  GAMMAMAXRAD=100.0;



  //////////////////
  // overrides for more detailed problem dependence
  //  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  //  TIMEORDER=4;
  TIMEORDER=3; // more smooth accurate solution than TIMEORDER=4 or 2 (midpoint or TVD)

  if(EOMTYPE==EOMFFDE||EOMTYPE==EOMFFDE2){
    //    lim[1]=lim[2]=lim[3]=DONOR;
    //    lim[1]=lim[2]=lim[3]=MINM;
    lim[1]=lim[2]=lim[3]=MC;
    //  lim[1]=lim[2]=lim[3]=WENO5BND;
    //lim[1]=lim[2]=lim[3]=PARAFLAT;
    //lim[1]=lim[2]=lim[3]=PARALINE;
  }
  else{
    //    lim[1]=lim[2]=lim[3]=DONOR;
    //    lim[1]=lim[2]=lim[3]=MINM;
    //lim[1]=lim[2]=lim[3]=MC;
    //  lim[1]=lim[2]=lim[3]=WENO5BND;
    //lim[1]=lim[2]=lim[3]=PARAFLAT;
    lim[1]=lim[2]=lim[3]=PARALINE;
  }

  //    cour=0.1;
  //    cour=0.5;
  //    cour=0.9; // works fine, but 0.8 more generally good.  Although sometimes cour=0.9 actually gives a bit smoother solution.
  //  cour=0.8;
  cour=0.49999; // 0.8 is too unstable for RADBEAM2D with curved flow relative to grid.

  if(DOWALDDEN){
    PALLLOOP(pl) fluxmethod[pl]=HLLFLUX; // lower errors in unresolved regions
    lim[1]=lim[2]=lim[3]=MC; // to preserve symmetry better
  }
  else{
    PALLLOOP(pl) fluxmethod[pl]=LAXFFLUX; //HLLFLUX;
    // HLL leads to problems with radiation and realistic opacities.
    PALLLOOP(pl) if(RADFULLPL(pl)) fluxmethod[pl]=LAXFFLUX;
  }

  //FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;
  

  //  rescaletype=1;
  rescaletype=4;

  //    if(DOWALDDEN) rescaletype=4;
  if(DOWALDDEN) rescaletype=5; // like 4, but b^2/rho scales as 1/r away from horizon

  BSQORHOLIMIT=1E3; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK // was 2E2 but 
  BSQOULIMIT=1E9; // was 1E3 but latest BC test had 1E4.  was 1E5 but needed like 1E7 to 1E8 to avoid gastemperature in funnel being repeatedly forced up even when Compton and other processes keep low.  Also makes next solution guess for implicit solver very different, and takes longer to converge. // Up to 1E9 to allow T same for higher BSQORHOLIMIT=1E3
  UORHOLIMIT=1E10; // has to be quite high, else hit floor in high optical depth cases and run-away injection of u and then rho.
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
  //OSMARK: where is DTr1 defined? what is DTfake?
  DTfake=MAX(1,DTr/10); 




  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==FLATNESS){
 
    //  lim[1]=lim[2]=lim[3]=MINM;
    //    cour=0.5;
    gam=gamideal=5.0/3.0;
    cooling=KORAL;

    BCtype[X1UP]=PERIODIC; // OUTFLOW;
    BCtype[X1DN]=PERIODIC;
    BCtype[X2UP]=PERIODIC; // OUTFLOW;
    BCtype[X2DN]=PERIODIC;
    BCtype[X3UP]=PERIODIC; // OUTFLOW;
    BCtype[X3DN]=PERIODIC;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.05;   // default dumping period

    DTr = 100; //number of time steps for restart dumps
    tf = 10.0; //final time
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/


  if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D){

    //    TIMEORDER=3;
    //    lim[1]=lim[2]=lim[3]=PARALINE;

    //    cour=0.1;
    //    cour=0.5;
    //    cour=0.9; // works fine, but 0.8 more generally good.
    cour=0.8;
    gam=gamideal=5.0/3.0;
    cooling=KORAL;

    BCtype[X1UP]=OUTFLOW;
    BCtype[X1DN]=OUTFLOW;
    BCtype[X2UP]=OUTFLOW;
    BCtype[X2DN]=OUTFLOW;
    BCtype[X3UP]=OUTFLOW; 
    BCtype[X3DN]=OUTFLOW;

    // sigmarad = 1.56E-64
    // arad=4*sigmarad/c
    // NOTE: Koral code has different values than paper
    ARAD_CODE=4.0*1.56E-64*(TEMPBAR*TEMPBAR*TEMPBAR*TEMPBAR); // to match koral and avoiding real units


    int idt;
    //  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E3;

    DTr = 100; //number of time steps for restart dumps
    if(WHICHPROBLEM==RADPULSEPLANAR){
      tf=1E5;
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=100.0; // Koral output steps
      //for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1; // testing
    }
    else if(WHICHPROBLEM==RADPULSE){
      tf = 35; //final time
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    }
    else if(WHICHPROBLEM==RADPULSE3D){
      tf = 70; //final time
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    }

    //    DODIAGEVERYSUBSTEP = 1;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  // TOTRY: SEe if koral is going Erf<0
  // TOTRY: matching F/E and seeing koral fails.

  if(WHICHPROBLEM==RADBEAMFLAT){
    //cour=0.8; // this or with old MINDTSET, causes Erf<0 for default koral test
    //    cour=0.5;
    //    lim[1]=lim[2]=lim[3]=MINM;
    //    lim[1]=lim[2]=lim[3]=PARALINE;
    // gam=gamideal=5.0/3.0;
    gam=gamideal=4.0/3.0; // koral now
    cooling=KORAL;

    //    RADBEAMFLAT_FRATIO=0.995; // koral at some point.
    RADBEAMFLAT_FRATIO=0.99995;
    RADBEAMFLAT_ERAD=1./RHOBAR; // 1g/cm^3 worth of energy density in radiation
    RADBEAMFLAT_RHO=1./RHOBAR; // 1g/cm^3
    RADBEAMFLAT_UU=0.1/RHOBAR; // 0.1g/cm^3 worth of energy density in fluid

    // avoid hitting gamma ceiling
    GAMMAMAXRAD=MAX(GAMMAMAXRAD,2.0*1.0/sqrt(1.0-RADBEAMFLAT_FRATIO*RADBEAMFLAT_FRATIO));


    BCtype[X1UP]=OUTFLOW;
    BCtype[X1DN]=RADBEAMFLATINFLOW;
    BCtype[X2UP]=OUTFLOW;
    BCtype[X2DN]=OUTFLOW;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    //  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.05;   // default dumping period
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;   // like problem24 in koral

    DTr = 100; //number of time steps for restart dumps
    tf = 10.0; //final time
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADTUBE){

    // 1,2,3,31,4,41,5
    //#define NTUBE 1
    #define NTUBE 31 // harder near t=0 at discontinuity
    //#define NTUBE 5
    //#define NTUBE 3

    //    lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    // should have PARA(LINE) not oscillate so much at cusp
    // Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
    //cour=0.5;
    cooling=KORAL;

    // arad = 4*sigmarad/c (so removed /4. from koral sigma setup).
    // Note, sigmarad or arad is not arbitrary -- value chosen to IC in radiative-hydro balance for each separately the left and right states.
    if(NTUBE==1){
      gam=gamideal=5./3.;
      ARAD_CODE=(1e-8/pow(calc_PEQ_Tfromurho(3.e-5/(gam-1.),1.),4.));
    }
    else if(NTUBE==2){
      gam=gamideal=5./3.;
      ARAD_CODE=(2e-5/pow(calc_PEQ_Tfromurho(4.e-3/(gam-1.),1.),4.));
    }
    else if(NTUBE==3){
      gam=gamideal=2.;
      ARAD_CODE=(2./pow(calc_PEQ_Tfromurho(60./(gam-1.),1.),4.));
    }
    else if(NTUBE==31){
      gam=gamideal=2.;
      ARAD_CODE=(2./pow(calc_PEQ_Tfromurho(60./(gam-1.),1.),4.));
    }
    else if(NTUBE==4){
      gam=gamideal=5./3.;
      ARAD_CODE=(.18/pow(calc_PEQ_Tfromurho(6.e-3/(gam-1.),1.),4.));
    }
    else if(NTUBE==41){
      gam=gamideal=5./3.;
      ARAD_CODE=(.18/pow(calc_PEQ_Tfromurho(6.e-3/(gam-1.),1.),4.));
    }
    else if(NTUBE==5){
      gam=gamideal=2.;
      ARAD_CODE=(2./pow(calc_PEQ_Tfromurho(60./(gam-1.),1.),4.));
    }


    trifprintf("RADTUBE NTUBE=%d ARAD_CODE=%g SIGMARAD_CODE=%g\n",NTUBE,ARAD_CODE,ARAD_CODE/4.0);


    BCtype[X1UP]=FREEOUTFLOW;
    BCtype[X1DN]=FREEOUTFLOW;
    BCtype[X2UP]=OUTFLOW; // NOTEMARK: Koral sets fixed BCs.  We can do that following the IC choices, but not necessary.
    BCtype[X2DN]=OUTFLOW;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    if(NTUBE==5){
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1.0;
    }
    else{
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=2.0;
    }

    DTr = 100; //number of time steps for restart dumps
    //  tf = 100.0; //final time (seems almost good enough to get quasi-steady solution for these steady tube tests)
    if(NTUBE==5) tf=15.0; // koral paper shows about t=13, and code has no problem going much further than 15.
    else tf = 3E2; //final time (good enough to see any evolution and errored evolution)
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADSHADOW){

    //    lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    // should have PARA(LINE) not oscillate so much at cusp
    // Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
    //    cour=0.5;
    gam=gamideal=1.4;
    cooling=KORAL;
    ARAD_CODE=1E7*1E-5*(2.5E-9/7.115025791e-10); // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

    RADSHADOW_NLEFT=0.99999;
    RADSHADOW_ANGLE=0.0;
    RADSHADOW_TLEFTOTAMB=100.0;
    RADSHADOW_BEAMY=0.3;

    // avoid hitting gamma ceiling
    GAMMAMAXRAD=MAX(GAMMAMAXRAD,2.0*1.0/sqrt(1.0-RADSHADOW_NLEFT*RADSHADOW_NLEFT));


    BCtype[X1UP]=FREEOUTFLOW;
    BCtype[X1DN]=RADSHADOWINFLOW;
    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.5;

    DTr = 100; //number of time steps for restart dumps
    //  tf = 100.0; //final time (seems almost good enough to get quasi-steady solution for these steady tube tests)
    tf = 10.0; //final time
  }




  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADDBLSHADOW){

    //    lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    // lim[1]=lim[2]=lim[3]=MC; // MC gets totally bonkers answer with NLEFT=0.99999
    //lim[1]=lim[2]=lim[3]=PARALINE; // bonkers answer for NLEFT=0.99999, ok for NLEFT=0.99
    // should have PARA(LINE) not oscillate so much at cusp
    // Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
    //cour=0.5;
    // cour=0.49; // doesn't help oscillations for NLEFT=0.99999 with MINM
    gam=gamideal=1.4;
    cooling=KORAL;
    // ARAD_CODE=1E-30;
    //    ARAD_CODE=1E7*1E-5*(2.5E-9/7.115025791e-10); // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.
    ARAD_CODE=1E7*1E-5*1E-10*(6E-9/1.7E-25); // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.  And also similar value as in Figure 11 koral paper plot.  As long as prad0<<u and prad0<<rho, solution is independent of ARAD because 4-force off radiation on the fluid is negligible.  Then kappa just sets what rho becomes \tau\sim 1 and nothing about the fluid is affected.



    //    RADDBLSHADOW_NLEFT=0.99999; // Works well with MINM (only 49 total failures at relatively early time for otherwise default setup).  very hard on code -- only MINM with jon choice for CASES works.
    //    RADDBLSHADOW_NLEFT=0.99; // koral paper
    //    RADDBLSHADOW_NLEFT=0.999; // latest koral (ok to use, weak oscillations with LAXF)
    //PALLLOOP(pl) fluxmethod[pl]=HLLFLUX; // smaller oscillations even at 0.99999


    //  RADDBLSHADOW_NLEFT=0.7;
    //  RADDBLSHADOW_NLEFT=0.93;
   
    //    angle=0.4; // koral paper
    //    angle=0.3; // latest koral

    //    RADDBLSHADOW_NLEFT=0.99; // what's in HARMRAD
    //    RADDBLSHADOW_NLEFT=0.99999; // works but noisy
    RADDBLSHADOW_NLEFT=0.999;
    RADDBLSHADOW_ANGLE=0.4;
    RADDBLSHADOW_TLEFTOTAMB=100.0;
    RADDBLSHADOW_BEAMY=0.3;

    // avoid hitting gamma ceiling
    GAMMAMAXRAD=MAX(GAMMAMAXRAD,2.0*1.0/sqrt(1.0-RADDBLSHADOW_NLEFT*RADDBLSHADOW_NLEFT));


    BCtype[X1UP]=FREEOUTFLOW;
    BCtype[X1DN]=RADSHADOWINFLOW;
    BCtype[X2UP]=RADSHADOWINFLOWX2UP;
    BCtype[X2DN]=ASYMM;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=2E-1;

    DTr = 100; //number of time steps for restart dumps
    //  tf = 100.0; //final time (seems almost good enough to get quasi-steady solution for these steady tube tests)
    // tf = 200.0; //final time (far past plot so see evolves or stationary).
    tf = 20.0; //final time (far past plot so see evolves or stationary).
  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADBEAM2DKS){

    RADBEAM2D_BEAMNO=1; // 1-4
    // whether constant or radially varying background
    // ==0 doesn't make much sense for Minkowski without gravity, because flow reverses due to chosen high density
    RADBEAM2D_FLATBACKGROUND=1;


    //lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    //    cour=0.5;
    // cour=0.2; // doesn't seem to help avoid failures for this test

  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    if(!(ISSPCMCOORDNATIVE(MCOORD))){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADBEAM2D\n",MCOORD);
      myexit(3434628752);
    }

    gam=gamideal=1.4;
    cooling=KORAL;
    // ARAD_CODE=ARAD_CODE_DEF*1E5; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

    BCtype[X1UP]=RADBEAM2DFLOWINFLOW;
    if(MCOORD==KSCOORDS||BLCOORDS){
      BCtype[X1DN]=HORIZONOUTFLOW; // if SPCMINKMETRIC with no gravity, suckingin on boundary can leave prad0 very small and then pradi~c for no good reason
    }
    else{
      //    BCtype[X1DN]=OUTFLOW;
      BCtype[X1DN]=FREEOUTFLOW;
    }
    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    // BCtype[X3UP]=FREEOUTFLOW;
    BCtype[X3UP]=OUTFLOW;
    BCtype[X3DN]=RADBEAM2DBEAMINFLOW;



    FTYPE DTOUT1;
    if (RADBEAM2D_BEAMNO==1){
      tf = 10.0*(M_PI*0.5)*3.0; //final time
      DTOUT1=tf/100.0; //dt for basic output
    }
    else if (RADBEAM2D_BEAMNO==2){
      tf = 10.0*(M_PI*0.5)*6.0; //final time
      DTOUT1=tf/100.0; //dt for basic output
    }
    else if (RADBEAM2D_BEAMNO==3){
      tf = 10.0*(M_PI*0.5)*15.0; //final time
      DTOUT1=tf/100.0; //dt for basic output
    }
    else if (RADBEAM2D_BEAMNO==4){
      tf = 10.0*(M_PI*0.5)*40.0; //final time
      DTOUT1=tf/100.0; //dt for basic output
    }

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=DTOUT1;
    //    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.001; // testing

    DTr = 100; //number of time steps for restart dumps

    //    DODIAGEVERYSUBSTEP = 1;

  }
  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADBEAM2DKSVERT){

    RADBEAM2DKSVERT_BEAMNO=5; // 1-5
    // whether constant or radially varying background
    // ==0 doesn't make much sense for Minkowski without gravity, because flow reverses due to chosen high density
    RADBEAM2D_FLATBACKGROUND=1;


    //    lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    //    cour=0.5;


  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    if(!(ISSPCMCOORDNATIVE(MCOORD))){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADBEAM2DKSVERT\n",MCOORD);
      myexit(3434628752);
    }

    gam=gamideal=1.4;
    cooling=KORAL;
    // ARAD_CODE=ARAD_CODE_DEF*1E5; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

    BCtype[X1UP]=RADBEAM2DFLOWINFLOW;
    //BCtype[X1DN]=OUTFLOW;
    BCtype[X1DN]=HORIZONOUTFLOW;
    BCtype[X2UP]=RADBEAM2DKSVERTBEAMINFLOW;
    BCtype[X2DN]=OUTFLOW;
    // BCtype[X3UP]=FREEOUTFLOW;
    BCtype[X3UP]=OUTFLOW;
    BCtype[X3DN]=OUTFLOW;


    FTYPE DTOUT1;
    if (RADBEAM2DKSVERT_BEAMNO==1){
      DTOUT1=1; //dt for basic output
    }
    else if (RADBEAM2DKSVERT_BEAMNO==2){
      DTOUT1=.4; //dt for basic output
    }
    else if (RADBEAM2DKSVERT_BEAMNO==3){
      DTOUT1=1.; //dt for basic output
    }
    else if (RADBEAM2DKSVERT_BEAMNO==4){
      DTOUT1=.25; //dt for basic output
    }
    else if (RADBEAM2DKSVERT_BEAMNO==5){
      DTOUT1=.4; //dt for basic output
    }

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=DTOUT1;
    // for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.001;

    DTr = 100; //number of time steps for restart dumps
    tf = 20.0; //final time

    //    DODIAGEVERYSUBSTEP = 1;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==ATMSTATIC){

    //    lim[1]=lim[2]=lim[3]=MINM; 
    //    lim[1]=lim[2]=lim[3]=PARALINE; // actually more error in u^r than MINM for inner radial boundary points (~factor of two larger u^r).
    // NOTE: with FTYPE as double, not enough precision to have good convergence for u^r -- just noise. See makefile.notes for how to go to ldouble and then has same error as koral.
    //    cour=0.5;


    if(!(ISSPCMCOORDNATIVE(MCOORD))){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for ATMSTATIC\n",MCOORD);
      myexit(3434628752);
    }

 //   a=0.0; // no spin in case use MCOORD=KSCOORDS
    gam=gamideal=1.4;
    cooling=KORAL;
    ARAD_CODE=0.0;

    // HORIZONOUTFLOW or HORIZONOUTFLOWSTATIC leads to little bit more static solution near inner radial boundary due to higher-order interpolation.  Could also fix values as in Koral, but odd to fix values for incoming flow.
    //    BCtype[X1UP]=HORIZONOUTFLOWSTATIC;
    //    BCtype[X1DN]=HORIZONOUTFLOWSTATIC;

    // FIXEDUSEPANALYTIC gives solution just like koral where no boundary effects, but a bit odd to generally fix values for incoming flow.
    BCtype[X1UP]=FIXEDUSEPANALYTIC;
    BCtype[X1DN]=FIXEDUSEPANALYTIC;

    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    BCtype[X3UP]=PERIODIC;
    BCtype[X3DN]=PERIODIC;


    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E2;

    DTr = 100; //number of time steps for restart dumps
    tf = 1E7; //final time

    //    FLUXDISSIPATION=0.0;

    //    DODIAGEVERYSUBSTEP = 1;

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADATM){

    //lim[1]=lim[2]=lim[3]=MINM; // MINM gets larger error and jump in v1 at outer edge
    //    lim[1]=lim[2]=lim[3]=PARALINE;
    // Koral uses MINMOD_THETA2 (MC?)
    // koral paper uses MP5
    //    lim[1]=lim[2]=lim[3]=MC;
    //    cour=0.5;
    

    if(!(ISSPCMCOORDNATIVE(MCOORD))){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADATM\n",MCOORD);
      myexit(3434628753);
    }

  //  a=0.0; // no spin in case use MCOORD=KSCOORDS
    gam=gamideal=1.4;
    cooling=KORAL;
    //    ARAD_CODE=0.0;

    //    BCtype[X1UP]=RADATMBEAMINFLOW;
    //    BCtype[X1DN]=RADATMBEAMINFLOW;
    // really same as above, just simpler to avoid mistakes and can focus on init.c
    BCtype[X1UP]=FIXEDUSEPANALYTIC;
    BCtype[X1DN]=FIXEDUSEPANALYTIC;

    //    BCtype[X1UP]=HORIZONOUTFLOW;
    //    BCtype[X1DN]=HORIZONOUTFLOW;
    //BCtype[X1DN]=OUTFLOW;

    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    BCtype[X3UP]=PERIODIC;
    BCtype[X3DN]=PERIODIC;


    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E6;
    //for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E5; // DEBUG

    DTr = 100; //number of time steps for restart dumps
    tf = 2E9; //final time
    
    //    tf=1E8; // profiling // SUPERTODOMARK

    //    DODIAGEVERYSUBSTEP = 1;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADWALL){

    //    lim[1]=lim[2]=lim[3]=MINM; // Messy with PARALINE
    //    cour=0.5;


    gam=gamideal=5.0/3.0;
    cooling=KORAL;

    BCtype[X1UP]=OUTFLOW;
    BCtype[X1DN]=RADWALLINFLOW;
    BCtype[X2UP]=RADWALLINFLOW;
    BCtype[X2DN]=ASYMM;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.5;

    DTr = 100; //number of time steps for restart dumps
    tf = 50.0; //final time
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADWAVE){

    //    lim[1]=lim[2]=lim[3]=MINM; // generates glitch at extrema in prad0
    //lim[1]=lim[2]=lim[3]=MC; // less of a glitch near extrema in prad0
    //    cour=0.5;


    cooling=KORAL;
    gam=gamideal=5./3.;


    // KORALTODO: See Jiang, Stone, Davis (2012) for S6.1.2 for linear MHD-radiation compressible wave tests
    
    RADWAVE_NWAVE=5; // 1,2,3,4,5 .  And for 5 can choose NUMERO=41,11,1
    //    RADWAVE_NWAVE=1; // GOOD
    //    RADWAVE_NWAVE=2; // GOOD
    //    RADWAVE_NWAVE=3; // GOOD
    //    RADWAVE_NWAVE=4; // gets noisy in prad1 by t~30 with MINM or MC  -- check koral when Olek makes it work.  KORALTODO
    //    RADWAVE_NUMERO=11; // GOOD
    RADWAVE_NUMERO=104;
    //RADWAVE_NUMERO=41; // OK if don't use check if can do explicit.  So use this to show how should more generally improve the tau based suppression check!  But, DAMPS significantly! Smaller IMPCONV doesn't help.  Check with koral KORALTODO.  MC doesn't help/change much.
    //RADWAVE_NUMERO=1; // wierd jello oscillations in prad0, and no wave motion -- like in koral though.  KORALTODO.  With only implicit, jello is different (smaller IMPCONV doesn't help and larger IMPEPS doesn't help).

    // NUMERO=41 corresponds to Jiang et al. (2002) PP=100, sigma=10 (2nd row, 2nd column in Table B1) 11 to PP=0.01, sigma=0.01 (1st row, 1st column).

    //NUMERO 1 was supposed to be his original test (1st row, 1st column) but, as you mention, it turned out to be jelly. The reason is that the initial conditions from the table were not precise enough to hit the acoustic mode and much faster radiation mode quickly dominates causing the jelly behavior. I had difficult time with that and decided to derive the numbers by myself (PROBLEMS/RADWAVE/disp107.nb). They were a bit different and made the difference so that one sees only the acoustic mode. The reason is most likely the fact that my dispersion relation and the code are somewhat relativistic.




    // defaults
    RADWAVE_KAPPA=RADWAVE_KAPPAES=0.0;


    /* Pre-SASHA
    if(RADWAVE_NWAVE==5){ //sound wave with radiation set up according to Jiang+12

      if(RADWAVE_NUMERO==41){
        RADWAVE_PP=100.;
        RADWAVE_CC=1.e2;
        RADWAVE_KAPPA=10.;
        RADWAVE_RHOFAC=0.01;
        RADWAVE_DRRE=(1.e-3*RADWAVE_RHOFAC);
        RADWAVE_DRIM=0.;
        RADWAVE_DVRE=(4.06372e-6*RADWAVE_RHOFAC);
        RADWAVE_DVIM=(6.90937e-6*RADWAVE_RHOFAC);
        RADWAVE_DURE=(9.88671e-4*RADWAVE_RHOFAC);
        RADWAVE_DUIM=(6.97077e-6*RADWAVE_RHOFAC);
        RADWAVE_DERE=(-4.52724e-5*RADWAVE_RHOFAC);
        RADWAVE_DEIM=(2.78566e-5*RADWAVE_RHOFAC);
        RADWAVE_DFRE=(-5.83678e-6*RADWAVE_RHOFAC);
        RADWAVE_DFIM=(-9.48194e-6*RADWAVE_RHOFAC);
        RADWAVE_OMRE=0.0255331;
        RADWAVE_OMIM=0.0434128;
        RADWAVE_DTOUT1=1.e-0;
      }

      if(RADWAVE_NUMERO==11){
        RADWAVE_PP=0.01;
        RADWAVE_CC=1.e2;
        RADWAVE_KAPPA=0.01;
        RADWAVE_RHOFAC=0.01;
        RADWAVE_DRRE=(1.e-3*RADWAVE_RHOFAC);
        RADWAVE_DRIM=0.;
        RADWAVE_DVRE=(9.99998e-6*RADWAVE_RHOFAC);
        RADWAVE_DVIM=(8.48878e-9*RADWAVE_RHOFAC);
        RADWAVE_DURE=(1.66666e-3*RADWAVE_RHOFAC);
        RADWAVE_DUIM=(2.82938e-6*RADWAVE_RHOFAC);
        RADWAVE_DERE=(1.95853e-8*RADWAVE_RHOFAC);
        RADWAVE_DEIM=(1.91123e-7*RADWAVE_RHOFAC);
        RADWAVE_DFRE=(-1.33508e-5*RADWAVE_RHOFAC);
        RADWAVE_DFIM=(4.23463e-6*RADWAVE_RHOFAC);
        RADWAVE_OMRE=6.28317e-2;
        RADWAVE_OMIM=5.33366e-5;
        RADWAVE_DTOUT1=1.e-0;
      }


      if(RADWAVE_NUMERO==1){
        RADWAVE_PP=0.01;
        RADWAVE_CC=1.e4;
        RADWAVE_KAPPA=0.01;
        RADWAVE_DRRE=1.e-3;
        RADWAVE_DRIM=0.;
        RADWAVE_DVRE=9.7548e-8;
        RADWAVE_DVIM=7.92788e-9;
        //RADWAVE_DPRE=1.61075e-3;
        //RADWAVE_DPIM=2.07402e-4;
        RADWAVE_DURE=1.57546e-3;
        RADWAVE_DUIM=2.57783e-4;
        RADWAVE_DERE=1.6874e-8; // 1.79137e-8
        RADWAVE_DEIM=9.48966e-9; // 8.56498e-9
        RADWAVE_DFRE=-1.77115e-6; // -1.32035e-6
        RADWAVE_DFIM=3.65291e-6; // 3.88814e-6
        RADWAVE_OMRE=6.12912e-4; // 7.99077
        RADWAVE_OMIM=4.98123e-5; // 0.512336
        RADWAVE_DTOUT1=1.e-2;
      }

      RADWAVE_RHOZERO=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=((1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC)) ; // to get proper sound speed
      RADWAVE_TEMP=(calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO)) ; // temperature from rho and uint
      ARAD_CODE=((3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP)); //to get the proper radiation to gas pressure ratio, PP=4 sig T^4 / P
      RADWAVE_ERAD=(calc_LTE_EfromT(RADWAVE_TEMP)) ; // to get thermal equilibrium, E=4 sig T^4
    }


    if(RADWAVE_NWAVE==1){ //density wave advected with the gas
      // NO RADIATION
      RADWAVE_PP=0.1;
      RADWAVE_CC=1.e6;
      RADWAVE_VX=1.e-3;
      RADWAVE_DTOUT1=(.05/RADWAVE_VX);
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_ERAD=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
    }

    if(RADWAVE_NWAVE==2){ //hydro sound wave
      // NO RADIATION
      RADWAVE_PP=0.01;
      RADWAVE_CC=1.e6;
      RADWAVE_DTOUT1=(.05*RADWAVE_CC);
      RADWAVE_VX=0.;
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_ERAD=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
    }

    if(RADWAVE_NWAVE==3){ //radiative density wave advected with the gas
      //      FLUXDISSIPATION=(0.0);
      RADWAVE_PP=10.;
      RADWAVE_CC=1.e6;
      RADWAVE_VX=1.e-2;
      RADWAVE_DTOUT1=(.0005/RADWAVE_VX);
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
      RADWAVE_ERAD=calc_LTE_EfromT(RADWAVE_TEMP);
      RADWAVE_KAPPAES=10.;
    }


    if(RADWAVE_NWAVE==4){ //sound wave with radiation, set up without the phase shifts etc.
      //      FLUXDISSIPATION=(0.0);
      RADWAVE_PP=1.;
      RADWAVE_CC=1.e2;
      RADWAVE_DTOUT1=(.005*RADWAVE_CC);
      RADWAVE_VX=0.;
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-1;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
      RADWAVE_ERAD=calc_LTE_EfromT(RADWAVE_TEMP);
      RADWAVE_KAPPA=100.;
      RADWAVE_ERADFACTOR=.5;
      RADWAVE_GASFACTOR=.5;
    }
      */ // end pre-SASHA

    // start post-SASHA

    if(RADWAVE_NWAVE==5){ //sound wave with radiation set up according to Jiang+12

      if(RADWAVE_NUMERO==1) { // sound wave
        RADWAVE_WAVETYPE=0;  //sound wave
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.;
        RADWAVE_PP=0.01;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00010000000000000002*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000015228426395939093*RADWAVE_RHOFAC;
        RADWAVE_DUIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DERE=0*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-2.4365482233502554e-8*RADWAVE_RHOFAC;
        RADWAVE_DFIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.6283185307179587;
        RADWAVE_OMIM=0;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==10) { // fast magnetosonic wave
        RADWAVE_WAVETYPE=2;  //fast
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.01;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.0001602940583015828*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=-0.00009790871410382318*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.00001522842639593907*RADWAVE_RHOFAC;
        RADWAVE_DUIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0.00016230255678865884*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_DERE=0*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-3.905642029683238e-8*RADWAVE_RHOFAC;
        RADWAVE_DFIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=2.3855930340017838e-8*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_OMRE=1.007157271948693;
        RADWAVE_OMIM=0;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==11) { // slow magnetosonic wave
        RADWAVE_WAVETYPE=1;  //slow
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.01;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00006177069527516586*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0.00010011836806552759*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.00001522842639593909*RADWAVE_RHOFAC;
        RADWAVE_DUIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=-0.0000625515978604029*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_DERE=-4.0433761094240253e-25*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-1.5050727782781537e-8*RADWAVE_RHOFAC;
        RADWAVE_DFIM=0.*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=-2.4394323183478818e-8*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0.*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.3881167249671897;
        RADWAVE_OMIM=0;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==101) { // radiation-modified sound wave
        RADWAVE_WAVETYPE=0;  //sound
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.0000997992249118626*RADWAVE_RHOFAC;
        RADWAVE_DVIM=2.552072175928721e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000015155652908079845*RADWAVE_RHOFAC;
        RADWAVE_DUIM=7.696929719530536e-7*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DERE=1.3314776991134588e-10*RADWAVE_RHOFAC;
        RADWAVE_DEIM=3.6001746512388956e-8*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-2.5247126486226934e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=7.400407810152034e-8*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.627057023634126;
        RADWAVE_OMIM=0.016035142398657175;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }

      if(RADWAVE_NUMERO==102) { // radiation-modified sound wave
        RADWAVE_WAVETYPE=0;  //sound
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.0002662507979198814*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.0000633514446524509*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000011706978034894262*RADWAVE_RHOFAC;
        RADWAVE_DUIM=1.881532710186292e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DERE=0.00020541918444857084*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0.00014985861843019722*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.000020730815318727886*RADWAVE_RHOFAC;
        RADWAVE_DFIM=0.000037755564579364684*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0*RADWAVE_RHOFAC;
        RADWAVE_OMRE=1.67290310151504;
        RADWAVE_OMIM=0.39804886622888025;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==103) { // radiation-modified sound wave
        RADWAVE_WAVETYPE=0;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00009290985655754982*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.000014438241113392203*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.00001179766900341804*RADWAVE_RHOFAC;
        RADWAVE_DUIM=3.0429210480592464e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DERE=1.9819771721015766e-6*RADWAVE_RHOFAC;
        RADWAVE_DEIM=2.206454751998895e-6*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-4.3677706188561827e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=4.316214437243016e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.5837698456145599;
        RADWAVE_OMIM=0.09071814442518213;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==104) { // radiation-modified sound wave
        RADWAVE_WAVETYPE=0;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00007748046325109175*RADWAVE_RHOFAC;
        RADWAVE_DVIM=3.583193353788357e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DURE=9.141343124149056e-6*RADWAVE_RHOFAC;
        RADWAVE_DUIM=3.83221342644378e-7*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=0*RADWAVE_RHOFAC;
        RADWAVE_DERE=-1.5219307455357163e-7*RADWAVE_RHOFAC;
        RADWAVE_DEIM=9.380406060306363e-7*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.000019366644866579624*RADWAVE_RHOFAC;
        RADWAVE_DFIM=-7.930468856866346e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=0*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.4868241082927276;
        RADWAVE_OMIM=0.02251386783330655;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==1001){ //radiation-modified fast magnetosonic
        RADWAVE_WAVETYPE=2;  //fast
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00016025131429328265*RADWAVE_RHOFAC;
        RADWAVE_DVIM=7.238312005077197e-7*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=-0.00009795442630848571*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=9.836789501779977e-7*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000015198360895974991*RADWAVE_RHOFAC;
        RADWAVE_DUIM=4.815752909936621e-7*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0.00016234366410161697*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-8.96662164240542e-7*RADWAVE_RHOFAC;
        RADWAVE_DERE=1.4842118188293356e-9*RADWAVE_RHOFAC;
        RADWAVE_DEIM=6.063223162955078e-8*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-3.9543271084234507e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=8.51051304663626e-8*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=2.3667952154599258e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=2.1118238693659835e-8*RADWAVE_RHOFAC;
        RADWAVE_OMRE=1.0068887034237715;
        RADWAVE_OMIM=0.004547965563908265;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==1101){ //radiation-modified slow magnetosonic
        RADWAVE_WAVETYPE=1;  //slow
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.0000615332754996702*RADWAVE_RHOFAC;
        RADWAVE_DVIM=1.8313980164851912e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0.00009897721183016215*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=6.541857910619384e-6*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000015017423510649482*RADWAVE_RHOFAC;
        RADWAVE_DUIM=1.2229894345580098e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=-0.00006148820918323782*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-5.883153382953974e-6*RADWAVE_RHOFAC;
        RADWAVE_DERE=1.9170283012363001e-10*RADWAVE_RHOFAC;
        RADWAVE_DEIM=2.187214582100525e-8*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-1.6518053269343755e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=7.175204181969396e-8*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=-2.2367888827266106e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=-7.43141463935117e-8*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.38662497252216144;
        RADWAVE_OMIM=0.011507013108777591;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==1002){ //radiation-modified fast magnetosonic, opt THICK
        RADWAVE_WAVETYPE=2;  //fast
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.0002784991109850316*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.000052380393168228656*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=-0.000028109315354609682*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=6.2558750019767086e-6*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000011730539980454472*RADWAVE_RHOFAC;
        RADWAVE_DUIM=1.7129010401392157e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0.00011016964855189374*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-4.033370850235303e-6*RADWAVE_RHOFAC;
        RADWAVE_DERE=0.0002072941184085973*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0.00013636362726103654*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.00001833308481505651*RADWAVE_RHOFAC;
        RADWAVE_DFIM=0.00003636638174654657*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=2.6758144678721757e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=1.2427179588260004e-6*RADWAVE_RHOFAC;
        RADWAVE_OMRE=1.7498615222037273;
        RADWAVE_OMIM=0.3291157167389043;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }

      if(RADWAVE_NUMERO==1003){ //radiation-modified fast magnetosonic, opt THICK
        RADWAVE_WAVETYPE=2;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00015921491906000263*RADWAVE_RHOFAC;
        RADWAVE_DVIM=4.267309467171878e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=-0.00009865659642379489*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=6.116424112313706e-6*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000013105512199532442*RADWAVE_RHOFAC;
        RADWAVE_DUIM=2.2690789485828284e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0.00016304450066324336*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-5.540155699476499e-6*RADWAVE_RHOFAC;
        RADWAVE_DERE=2.9534637071306857e-6*RADWAVE_RHOFAC;
        RADWAVE_DEIM=1.5968078176265759e-6*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-2.7219589023767955e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=6.086535550665933e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=3.2386284213141576e-9*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=2.3827073246896776e-8*RADWAVE_RHOFAC;
        RADWAVE_OMRE=1.0003768401215956;
        RADWAVE_OMIM=0.02681229614532269;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==1004){ //radiation-modified fast magnetosonic, opt thin
        RADWAVE_WAVETYPE=2;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00015164754212433659*RADWAVE_RHOFAC;
        RADWAVE_DVIM=2.9243074320866887e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=-0.00011147156590526742*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=6.596175366942675e-7*RADWAVE_RHOFAC;
        RADWAVE_DURE=9.205934027867234e-6*RADWAVE_RHOFAC;
        RADWAVE_DUIM=7.437886816058989e-7*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=0.00017478715294162717*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-1.8658034881621974e-6*RADWAVE_RHOFAC;
        RADWAVE_DERE=-4.7021288964385735e-7*RADWAVE_RHOFAC;
        RADWAVE_DEIM=1.9823405974044627e-6*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.00003794222976803031*RADWAVE_RHOFAC;
        RADWAVE_DFIM=-3.18097469382448e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=0.00002693491298739615*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=2.670466935767298e-6*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.952829608545529;
        RADWAVE_OMIM=0.018373965490963148;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      
      if(RADWAVE_NUMERO==1102){ //radiation-modified slow magnetosonic, opt THICK
        RADWAVE_WAVETYPE=1;  //slow
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00008342687348874559*RADWAVE_RHOFAC;
        RADWAVE_DVIM=0.000012082877629545738*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0.00011363325579507992*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0.00027269723955801375*RADWAVE_RHOFAC;
        RADWAVE_DURE=9.461886706340277e-6*RADWAVE_RHOFAC;
        RADWAVE_DUIM=1.2137574760252086e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=-0.0000803822945939874*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-0.00030311425160368743*RADWAVE_RHOFAC;
        RADWAVE_DERE=0.000025966624942266354*RADWAVE_RHOFAC;
        RADWAVE_DEIM=0.00009678910913258438*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.000019826286725607242*RADWAVE_RHOFAC;
        RADWAVE_DFIM=5.476096510182763e-6*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=3.6607454416002e-6*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=-1.1474975240229893e-6*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.5241865057284164;
        RADWAVE_OMIM=0.0759189591904107;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }

      if(RADWAVE_NUMERO==1103){ //radiation-modified slow magnetosonic, opt THICK
        RADWAVE_WAVETYPE=1;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=0.1;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=10;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.000055107074739986735*RADWAVE_RHOFAC;
        RADWAVE_DVIM=6.81771795916397e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0.00007683477442788993*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=0.000018069585727045885*RADWAVE_RHOFAC;
        RADWAVE_DURE=0.000010353641210907202*RADWAVE_RHOFAC;
        RADWAVE_DUIM=2.5489931913014527e-6*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=-0.00004163521053362644*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-0.00001542206148990843*RADWAVE_RHOFAC;
        RADWAVE_DERE=9.058920389631598e-7*RADWAVE_RHOFAC;
        RADWAVE_DEIM=1.8594869909418343e-6*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-3.830409107982378e-7*RADWAVE_RHOFAC;
        RADWAVE_DFIM=1.992679542320905e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=2.114058251324126e-9*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=-6.394153931776251e-9*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.346247962327932;
        RADWAVE_OMIM=0.04283698530951345;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
        
      if(RADWAVE_NUMERO==1104){ //radiation-modified slow magnetosonic, opt thin
        RADWAVE_WAVETYPE=1;
        RADWAVE_RHOFAC=0.001;
        RADWAVE_B0=0.10075854437197568;
        RADWAVE_PP=10;
        RADWAVE_CC=10;
        RADWAVE_KAPPA=0.1;
        RADWAVE_DRRE=0.001*RADWAVE_RHOFAC;
        RADWAVE_DRIM=0*RADWAVE_RHOFAC;
        RADWAVE_DVRE=0.00005019053133722909*RADWAVE_RHOFAC;
        RADWAVE_DVIM=2.3866183515916313e-6*RADWAVE_RHOFAC;
        RADWAVE_DV2RE=0.00006765290591654829*RADWAVE_RHOFAC;
        RADWAVE_DV2IM=3.826394449423638e-6*RADWAVE_RHOFAC;
        RADWAVE_DURE=9.134497469805529e-6*RADWAVE_RHOFAC;
        RADWAVE_DUIM=2.481944791269412e-7*RADWAVE_RHOFAC;
        RADWAVE_DB2RE=-0.00003511412717901583*RADWAVE_RHOFAC;
        RADWAVE_DB2IM=-1.2206629792761417e-6*RADWAVE_RHOFAC;
        RADWAVE_DERE=-7.354752234301807e-8*RADWAVE_RHOFAC;
        RADWAVE_DEIM=6.007447328384492e-7*RADWAVE_RHOFAC;
        RADWAVE_DFRE=-0.000012540740009734634*RADWAVE_RHOFAC;
        RADWAVE_DFIM=-5.536217727919607e-7*RADWAVE_RHOFAC;
        RADWAVE_DF2RE=-0.000014894816213683286*RADWAVE_RHOFAC;
        RADWAVE_DF2IM=-5.731053967780647e-6*RADWAVE_RHOFAC;
        RADWAVE_OMRE=0.3153564090576144;
        RADWAVE_OMIM=0.0149955653605657;
        RADWAVE_DTOUT1=2*M_PI/RADWAVE_OMRE/10.;
      }
      RADWAVE_RHOZERO=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=((1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC)) ; // to get proper sound speed
      RADWAVE_TEMP=(calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO)) ; // temperature from rho and uint
      ARAD_CODE=((3.*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP)); //to get the proper radiation to gas pressure ratio, PP=4 sig T^4 / P
      RADWAVE_ERAD=(calc_LTE_EfromT(RADWAVE_TEMP)) ; // to get thermal equilibrium, E=4 sig T^4
      
      dualfprintf(fail_file,"RADWAVE_RHOZERO=%g, RADWAVE_KK=%g, RADWAVE_UINT=%g, RADWAVE_TEMP=%g, ARAD_CODE=%g, RADWAVE_ERAD=%21.15g\n",
                  RADWAVE_RHOZERO, RADWAVE_KK, RADWAVE_UINT, RADWAVE_TEMP, ARAD_CODE, RADWAVE_ERAD);
      if(RADWAVE_NWAVE==5){
        FILE *out;
        if((out=fopen("radtestparams.dat","wt"))==NULL){
          dualfprintf(fail_file,"Couldn't write radtestparams.dat file\n");
          myexit(1);
        }
        else{
          fprintf(out,"#%20s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n",
                  "RADWAVE_RHOZERO",
                  "RADWAVE_KK",
                  "RADWAVE_UINT",
                  "RADWAVE_ERAD",
                  "RADWAVE_DRRE",
                  "RADWAVE_RHOFAC",
                  "RADWAVE_B0",
                  "RADWAVE_PP",
                  "RADWAVE_CC",
                  "RADWAVE_KAPPA",
                  "RADWAVE_DRRE",
                  "RADWAVE_DRIM",
                  "RADWAVE_DVRE",
                  "RADWAVE_DVIM",
                  "RADWAVE_DV2RE",
                  "RADWAVE_DV2IM",
                  "RADWAVE_DURE",
                  "RADWAVE_DUIM",
                  "RADWAVE_DB2RE",
                  "RADWAVE_DB2IM",
                  "RADWAVE_DERE",
                  "RADWAVE_DEIM",
                  "RADWAVE_DFRE",
                  "RADWAVE_DFIM",
                  "RADWAVE_DF2RE",
                  "RADWAVE_DF2IM",
                  "RADWAVE_OMRE",
                  "RADWAVE_OMIM",
                  "RADWAVE_DTOUT1",
                  "RADWAVE_WAVETYPE");
          fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
                  RADWAVE_RHOZERO,
                  RADWAVE_KK,
                  RADWAVE_UINT,
                  RADWAVE_ERAD,
                  RADWAVE_DRRE,
                  RADWAVE_RHOFAC,
                  RADWAVE_B0,
                  RADWAVE_PP,
                  RADWAVE_CC,
                  RADWAVE_KAPPA,
                  RADWAVE_DRRE,
                  RADWAVE_DRIM,
                  RADWAVE_DVRE,
                  RADWAVE_DVIM,
                  RADWAVE_DV2RE,
                  RADWAVE_DV2IM,
                  RADWAVE_DURE,
                  RADWAVE_DUIM,
                  RADWAVE_DB2RE,
                  RADWAVE_DB2IM,
                  RADWAVE_DERE,
                  RADWAVE_DEIM,
                  RADWAVE_DFRE,
                  RADWAVE_DFIM,
                  RADWAVE_DF2RE,
                  RADWAVE_DF2IM,
                  RADWAVE_OMRE,
                  RADWAVE_OMIM,
                  RADWAVE_DTOUT1,
                  (double)RADWAVE_WAVETYPE
                  );
          fclose(out);
        }
      }

    }


    if(RADWAVE_NWAVE==1){ //density wave advected with the gas
      // NO RADIATION
      RADWAVE_PP=0.1;
      RADWAVE_CC=1.e6;
      RADWAVE_VX=1.e-3;
      RADWAVE_DTOUT1=(.05/RADWAVE_VX);
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_ERAD=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
    }

    if(RADWAVE_NWAVE==2){ //hydro sound wave
      // NO RADIATION
      RADWAVE_PP=0.01;
      RADWAVE_CC=1.e6;
      RADWAVE_DTOUT1=(.05*RADWAVE_CC);
      RADWAVE_VX=0.;
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_ERAD=1.;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
    }

    if(RADWAVE_NWAVE==3){ //radiative density wave advected with the gas
      //      FLUXDISSIPATION=(0.0);
      RADWAVE_PP=10.;
      RADWAVE_CC=1.e6;
      RADWAVE_VX=1.e-2;
      RADWAVE_DTOUT1=(.0005/RADWAVE_VX);
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-5;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
      RADWAVE_ERAD=calc_LTE_EfromT(RADWAVE_TEMP);
      RADWAVE_KAPPAES=10.;
    }


    if(RADWAVE_NWAVE==4){ //sound wave with radiation, set up without the phase shifts etc.
      //      FLUXDISSIPATION=(0.0);
      RADWAVE_PP=1.;
      RADWAVE_CC=1.e2;
      RADWAVE_DTOUT1=(.005*RADWAVE_CC);
      RADWAVE_VX=0.;
      RADWAVE_RHOZERO=1.;
      RADWAVE_AAA=1.e-1;
      RADWAVE_KK=2.*Pi;
      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
      RADWAVE_ERAD=calc_LTE_EfromT(RADWAVE_TEMP);
      RADWAVE_KAPPA=100.;
      RADWAVE_ERADFACTOR=.5;
      RADWAVE_GASFACTOR=.5;
    }

//    if(RADWAVE_NWAVE==5){ //fast radiation-modified magnetosonic waves #1001 //Sasha
//      FLUXDISSIPATION=(0.0);
//      RADWAVE_PP=0.1;
//      RADWAVE_CC=10.;
//      RADWAVE_DTOUT1=(.005*RADWAVE_CC);
//      RADWAVE_VX=0.;
//      RADWAVE_RHOZERO=1.;
//      RADWAVE_AAA=1.e-1;
//      RADWAVE_KK=2.*Pi;
//      RADWAVE_UINT=(1./RADWAVE_CC/RADWAVE_CC)*RADWAVE_RHOZERO/gam/(gam-1.-1./RADWAVE_CC/RADWAVE_CC);
//      RADWAVE_TEMP=calc_PEQ_Tfromurho(RADWAVE_UINT,RADWAVE_RHOZERO);
//      ARAD_CODE=(3.0*RADWAVE_PP*(gam-1.)*RADWAVE_UINT/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP/RADWAVE_TEMP);
//      RADWAVE_ERAD=calc_LTE_EfromT(RADWAVE_TEMP);
//      RADWAVE_KAPPA=0.1;
//      RADWAVE_ERADFACTOR=.5;
//      RADWAVE_GASFACTOR=.5;
//    }


    // end post-SASHA


    BCtype[X1UP]=PERIODIC;
    BCtype[X1DN]=PERIODIC;
    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    BCtype[X3UP]=PERIODIC; 
    BCtype[X3DN]=PERIODIC;

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=RADWAVE_DTOUT1;

    DTr = 100; //number of time steps for restart dumps
    if(RADWAVE_VX==0.0) tf = MAX(100.0*RADWAVE_DTOUT1,5.0/RADWAVE_CC);
    else tf = MAX(100.0*RADWAVE_DTOUT1,5.0/MIN(RADWAVE_VX,RADWAVE_CC));

    if(RADWAVE_NWAVE==5){
      tf = 1.1*2*M_PI/RADWAVE_OMRE;
    }



    //    DODIAGEVERYSUBSTEP = 1;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==KOMIPROBLEM){
    
    lim[1]=lim[2]=lim[3]=MINM; // ok with fast shock
    //lim[1]=lim[2]=lim[3]=PARALINE; // good with fast shock
    //    lim[1]=lim[2]=lim[3]=MC;  // not too good with fast shock
    //    lim[1]=lim[2]=lim[3]=MP5;  // kinda ok with fast shock
    cour=0.499;
    cooling=NOCOOLING;
    gam=gamideal=4./3.;
    GAMMAMAX=100.0;
    BSQORHOLIMIT=1E10;
    BSQOULIMIT=1E10;
    UORHOLIMIT=1E10;


    BCtype[X1UP]=FREEOUTFLOW;
    BCtype[X1DN]=FREEOUTFLOW;
    BCtype[X2UP]=OUTFLOW;
    BCtype[X2DN]=OUTFLOW;
    BCtype[X3UP]=PERIODIC;
    BCtype[X3DN]=PERIODIC;
    
    DTr = 100; //number of time steps for restart dumps
    
    //set final time
    
    //fast shock
    if(WHICHKOMI==1){
      tf = 2.5;
    }
    //slow shock
    else if(WHICHKOMI==2){
      tf = 2.0;
    }
    //fast switch-off rarefaction
    else if(WHICHKOMI==3){
      tf = 1.0;
    }
    //slow switch-on rarefaction
    else if(WHICHKOMI==4){
      tf = 2.0;
    }
    //alfven wave
    else if(WHICHKOMI==5){
      tf = 2.0;
    }
    //compound wave
    else if(WHICHKOMI==6){
      tf = 1.5;  //also 0.1 and 0.75 are other times
    }
    //Shock tube 1
    else if(WHICHKOMI==7){
      tf = 1.0;
    }
    //Shock tube 2
    else if(WHICHKOMI==8){
      tf = 1.0;
    }
    //Collision
    else if(WHICHKOMI==9){
      tf = 1.22;
    }

    if(WHICHKOMI>=101 && WHICHKOMI<=109){
      GAMMAMAX=2000.0;
    }


    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1*tf;
  }



  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADBONDI){

    // lim[1]=lim[2]=lim[3]=MINM; // too low order for ~100 points
    //    lim[1]=lim[2]=lim[3]=PARALINE;
    //    cour=0.5;


    if(!(ISSPCMCOORDNATIVE(MCOORD))){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADBONDI\n",MCOORD);
      myexit(3434628752);
    }

 //   a=0.0; // no spin in case use MCOORD=KSCOORDS
    cooling=KORAL;
    // ARAD_CODE=ARAD_CODE_DEF*1E5; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

    BCtype[X1UP]=RADBONDIINFLOW;
    //    BCtype[X1DN]=OUTFLOW;
    BCtype[X1DN]=HORIZONOUTFLOW; // although more specific extrapolation based upon solution might work better
    BCtype[X2UP]=PERIODIC;
    BCtype[X2DN]=PERIODIC;
    // BCtype[X3UP]=FREEOUTFLOW;
    BCtype[X3UP]=OUTFLOW;
    BCtype[X3DN]=OUTFLOW;



    RADBONDI_TESTNO=2;

    if(RADBONDI_TESTNO==0){ // E1T6 (#1 in koral Table 5)
      RADBONDI_PRADGAS=1.2e-4/RHOBAR;
      RADBONDI_TGAS0=1e6/TEMPBAR;
      RADBONDI_MDOTPEREDD=1.;
    }

    if(RADBONDI_TESTNO==1){ // E10T5 (#2 in koral Table 5 and #1 in Fragile paper)
      RADBONDI_PRADGAS=1.2e-7/RHOBAR;
      RADBONDI_TGAS0=1e5/TEMPBAR;
      RADBONDI_MDOTPEREDD=10.;
    }

    if(RADBONDI_TESTNO==2){ // E10T6 (#3 in koral Table 5 and #3 in Fragile paper)
      RADBONDI_PRADGAS=1.2e-4/RHOBAR;
      RADBONDI_TGAS0=1.e6/TEMPBAR;
      RADBONDI_MDOTPEREDD=10.;
    }

    if(RADBONDI_TESTNO==3){ // E10T7 (#4 in koral Table 5 and #5 in Fragile paper)
      RADBONDI_PRADGAS=1.2e-1/RHOBAR;
      RADBONDI_TGAS0=1e7/TEMPBAR;
      RADBONDI_MDOTPEREDD=10.;
    }

    // koral skips E30T6 that is #6 in Fragile paper

    if(RADBONDI_TESTNO==4){ // E100T6 (#5 in koral Table 5 and #7 in Fragile paper)
      RADBONDI_PRADGAS=1.2e-5/RHOBAR; // note koral paper has 1.2E-4 that is wrong.
      RADBONDI_TGAS0=1e6/TEMPBAR;
      RADBONDI_MDOTPEREDD=100.;
    } 

    // koral skips E300T6 that is #8 in Fragile paper

    RADBONDI_MDOTEDD=(2.23/16.*1e18*MPERSUN)/(MBAR/TBAR); //Mdot converted to code units
    gam=gamideal=(1.+1./3.*((1.+RADBONDI_PRADGAS)/(.5+RADBONDI_PRADGAS)));

    trifprintf("RADBONDI: %g %g %g %g %g\n",RADBONDI_PRADGAS,RADBONDI_TGAS0,RADBONDI_MDOTPEREDD,RADBONDI_MDOTEDD,gam);

    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=10.0;

    DTr = 100; //number of time steps for restart dumps
    tf = 100*DTdumpgen[0]; // 100 dumps(?)

    //    DODIAGEVERYSUBSTEP = 1;

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADDOT){

    //    lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
    //    cour=0.5;

 //   a=0.0; // no spin in case use MCOORD=KSCOORDS
    gam=gamideal=4.0/3.0;
    cooling=KORAL;
    ARAD_CODE=ARAD_CODE_DEF*1E-20; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

    BCtype[X1UP]=OUTFLOW;
    BCtype[X1DN]=OUTFLOW;
    BCtype[X2UP]=OUTFLOW;
    BCtype[X2DN]=OUTFLOW;
    BCtype[X3UP]=OUTFLOW;
    BCtype[X3DN]=OUTFLOW;

    
    int idt;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.05;

    DTr = 100; //number of time steps for restart dumps
    tf = 100.0*DTdumpgen[0];

    //    DODIAGEVERYSUBSTEP = 1;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK || WHICHPROBLEM==RADDONUT || WHICHPROBLEM==RADCYLBEAM || WHICHPROBLEM==RADCYLBEAMCART){


    // NOTENOTENOTE: Also do following before running with RADDONUT:
    // 1) coord.c: rbr=5E2 -> 1E2 OR use Sasha or Jon coordinate setup for all parameters.
    // 2) DORADFIXUPS 0
    // 3) BSQORHOLIMIT=1E2;
    // 4) #define CONNDERTYPE DIFFGAMMIE
    // 5) #define N1 32

    TIMEORDER=2; // faster and sufficient


    int set_fieldtype(void);
    int FIELDTYPE=set_fieldtype();

    if(FIELDTYPE==FIELDJONMAD){
      // then funnel becomes too optically thick and traps radiation and accelerates radiation down into BH, leading to bad physical energy conservation even if total energy-momentum conservation equations used very well.
      BSQORHOLIMIT=3E2; // back to 100 if using YFL1 and AVOIDTAUFLOOR==2 in phys.tools.rad.c
      BSQOULIMIT=1E9;
    }
    else{
      BSQORHOLIMIT=1E2;
      BSQOULIMIT=1E8;
    }



    if(DOWALDDEN){// override if doing wald
      BSQORHOLIMIT=BSQORHOWALD;
    }
 
    UULIMIT=1E-20; // something small


    // Torus setups:
    // 3) gam=gamideal=5.0/3.0; RADNT_ELL=4.5; RADNT_UTPOT=0.9999999; RADNT_ROUT=2.0;  RADNT_RHODONUT=3.0; RADDONUT_OPTICALLYTHICKTORUS=1;     RADNT_KKK=1.e-1 * (1.0/powl(RADNT_RHODONUT,gam-1.0)); // where KKK doesn't matter.
    // ATM setups:
    // 1,2,3)  RADNT_ROUT=2.0; RADNT_RHOATMMIN=RADNT_RHODONUT*1E-4;  RADNT_TGASATMMIN = 1.e9/TEMPBAR;  RADNT_UINTATMMIN= (calc_PEQ_ufromTrho(RADNT_TGASATMMIN,RADNT_RHOATMMIN));    RADNT_TRADATMMIN = 1.e7/TEMPBAR;    RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));




    /////////////////////////////////
    // DONUT selections
    RADNT_DONUTTYPE=DONUTTHINDISK3; // SUPERMADNEW
    //RADNT_DONUTTYPE=DONUTOLEK;
    //RADNT_DONUTTYPE=DONUTOHSUGA;
    RADDONUT_OPTICALLYTHICKTORUS=1; // otherwise, pressure only from gas.
    RADNT_INFLOWING=0;
    RADNT_OMSCALE=1.0;
    // gas (not radiation) EOS \gamma value:
    gam=gamideal=5.0/3.0; // Ohsuga choice, assumes pairs not important.

    FTYPE gamtorus;

  if(WHICHPROBLEM==RADDONUT){
    if(RADNT_DONUTTYPE==DONUTTHINDISK||RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      // NOTEMARK: SUPERNOTE: Should set njet=0.0 in coord.c for thin disk.
      // NOTEMARK: Also choose npow=4.0;
      // NOTEMARK: Choose rbr=1E2;
      // NOTEMARK: Choose Qjet=1.9
      // NOTEMARK: htheta=0.1;
      // NOTEMARK: h0=0.1 instead of h0=0.3
      // NOTEMARK: Force only theta1=th0;
      //      h_over_r=0.1;
      h_over_r=0.1;// SUPERMADNEW
      h_over_r_jet=2.0*h_over_r;
    }
    if(RADNT_DONUTTYPE==DONUTOLEK || RADNT_DONUTTYPE==DONUTOHSUGA){
      //  h_over_r=0.3;
      h_over_r=0.2;
      h_over_r_jet=2.0*h_over_r;
    }

  }





    /////////////////////////////////////
    // DONUT TYPE and PARAMETERS
    //////////////////////////////////

    if(RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      RADDONUT_OPTICALLYTHICKTORUS=1; // otherwise, pressure only from gas.
      if(RADDONUT_OPTICALLYTHICKTORUS==1) gamtorus=4.0/3.0; // then should be as if gam=4/3 so radiation supports torus properly at t=0
      else gamtorus=gam;

      if(1){// SUPERMADNEW
        RADNT_RHODONUT=1E-2;
        RADNT_RHODONUT*=40.0;
 
        int set_fieldtype(void);
        int FIELDTYPE=set_fieldtype();

 
        if(a==0.8 && FIELDTYPE==FIELDJONMAD){
          trifprintf("Condition=1%21.15g\n",1);
          RADNT_RHODONUT/=(2.0*138.0);
          RADNT_RHODONUT/=(2.8); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(4.4); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(3.75); // Mdot\sim 135Ledd/c^2
        }

        trifprintf("spin=%g",a);
        trifprintf("FIELD=%g",FIELDTYPE);
        trifprintf("Jon=%g",FIELDJONMAD);
        if(a==0.5 && FIELDTYPE==FIELDJONMAD){  //Danilo
          trifprintf("Condition=2%21.15g\n",2);
          RADNT_RHODONUT/=(2.0*138.0);
          RADNT_RHODONUT/=(2.8); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(4.4); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(3.75); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT/=(143.0); // Mdot\sim 0.3Ledd/c^2
          RADNT_RHODONUT/=(3.);
          RADNT_RHODONUT*=(195.12195); //RIAF-beta60-2
          RADNT_RHODONUT/=(2.82);
          RADNT_RHODONUT/=(1.02);
          //RADNT_RHODONUT/=(1.003);
          trifprintf("After Conditionals: RADNT_RHODONUT=%21.15g\n",RADNT_RHODONUT);
        }

        if(a==0.8 && FIELDTYPE!=FIELDJONMAD){
          trifprintf("Condition=3%21.15g\n",3);
          RADNT_RHODONUT/=(33.0);
          RADNT_RHODONUT/=(2.7); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(2.6); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(1.4); // Mdot\sim 135Ledd/c^2
        }

        if(a==0.0 && FIELDTYPE==FIELDJONMAD){
          trifprintf("Condition=4%21.15g\n",4);
          RADNT_RHODONUT/=(2.0*138.0);
          RADNT_RHODONUT/=(2.8); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(3.9); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(1.8); // Mdot\sim 135Ledd/c^2
        }

        if(a==0.0 && FIELDTYPE!=FIELDJONMAD){
          trifprintf("Condition=5%21.15g\n",5);
          RADNT_RHODONUT/=(33.0);
          RADNT_RHODONUT/=(12.5); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(1.4); // Mdot\sim 135Ledd/c^2
          RADNT_RHODONUT*=(0.71); // Mdot\sim 135Ledd/c^2
        }
        RADNT_TRADATMMIN = 1.e5/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));

      }
      if(0){
        RADNT_RHODONUT=1E-2; // NT73 with MBH=10msun, a=0, Mdot=5Ledd/c^2
        RADNT_TRADATMMIN = 1.e5/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));
      }
      if(0){
        RADNT_RHODONUT=1E-10;
        RADNT_TRADATMMIN = 1.e5/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));
      }
      if(0){
        RADNT_RHODONUT=1.0;
        RADNT_TRADATMMIN = 1.e7/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));
      }

    }
    else if(RADNT_DONUTTYPE==DONUTOLEK){
      if(1){
        RADDONUT_OPTICALLYTHICKTORUS=1; // otherwise, pressure only from gas.
        // Mdot~135Ledd/c^2
        //
        if(RADDONUT_OPTICALLYTHICKTORUS==1) gamtorus=4.0/3.0; // then should be as if gam=4/3 so radiation supports torus properly at t=0
        else gamtorus=gam;

        //    RADNT_RHODONUT=1E-5;
        RADNT_RHODONUT=3.0; // gives 0.26 final density peak if RAD_ELL=3.5
        //RADNT_RHODONUT = KORAL2HARMRHO(1.0); // equivalent to koral's non-normalization
        //    RADNT_ELL=4.5; // torus specific angular momentum
        RADNT_ELL=4.5; // torus specific angular momentum
        RADNT_UTPOT=0.9999999; // scales rin for donut
        RADNT_KKK=1.e-1 * (1.0/pow(RADNT_RHODONUT,gamtorus-1.0)); // no effect with the scaling with density put in.

        RADNT_TRADATMMIN = 1.e7/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));

      }


      if(0){
        // THIN DISK with Mdot~7Ledd/c^2
        RADDONUT_OPTICALLYTHICKTORUS=0; // otherwise, pressure only from gas.
        
        if(RADDONUT_OPTICALLYTHICKTORUS==1) gamtorus=4.0/3.0; // then should be as if gam=4/3 so radiation supports torus properly at t=0
        else gamtorus=gam;
        
        RADNT_RHODONUT=3.0/2E4; // gives 0.26 final density peak if RAD_ELL=3.5
        RADNT_ELL=4.5; // torus specific angular momentum
        RADNT_UTPOT=0.9999999; // scales rin for donut
        RADNT_KKK=1.e-1 * ((gamtorus-1.0)/(gamtorus)/pow(RADNT_RHODONUT,gamtorus-1.0)); // no effect with the scaling with density put in.

        RADNT_TRADATMMIN = 1.e7/TEMPBAR;
        RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));
      }

      
    }
    else{
      RADDONUT_OPTICALLYTHICKTORUS=1; // otherwise, pressure only from gas.
      
      if(RADDONUT_OPTICALLYTHICKTORUS==1) gamtorus=4.0/3.0; // then should be as if gam=4/3 so radiation supports torus properly at t=0
      else gamtorus=gam;
      
      RADNT_RHODONUT=1E-2; // actual torus maximum density
      RADNT_DONUTRADPMAX=20.0; // radius of pressure maximum
      RADNT_HOVERR=0.5; // H/R\sim c_s/v_K at torus pressure maximum
      RADNT_LPOW=0.1; // l\propto r^lpow , lpow=0.5 would be Keplerian, lpow=0 would be constant angular momentum.

      RADNT_TRADATMMIN = 1.e7/TEMPBAR;
      RADNT_ERADATMMIN= (calc_LTE_EfromT(RADNT_TRADATMMIN));
    }


    ///////////////////////////////////
    // DONUT atmosphere:
    RADNT_ROUT=2.0; // what radius ATMMIN things are defining
    //RADNT_RHOATMMIN=KORAL2HARMRHO(1.e-4);
    //    RADNT_RHOATMMIN= KORAL2HARMRHO(1.e-2); // current koral choice
    RADNT_RHOATMMIN=RADNT_RHODONUT*1E-6;
    //    RADNT_TGASATMMIN = 1.e11/TEMPBAR;
    RADNT_TGASATMMIN = 1.e9/TEMPBAR;

    if(1){ // SUPERMADNEW
      RADNT_RHOATMMIN=RADNT_RHODONUT*1E-5;
    }

    RADNT_UINTATMMIN= (calc_PEQ_ufromTrho(RADNT_TGASATMMIN,RADNT_RHOATMMIN));
    // need external radiation energy density to be lower than interior of torus, else drives photons into torus from overpressured atmosphere and is more difficult to evolve.
    //    RADNT_TRADATMMIN = 1.e9/TEMPBAR;

    trifprintf("RADNT_RHOATMMIN=%g RADNT_RHOATMMIN=%g RADNT_UINTATMMIN=%g RADNT_ERADATMMIN=%g RADNT_RHODONUT=%g \n",RADNT_RHOATMMIN,RADNT_RHOATMMIN,RADNT_UINTATMMIN,RADNT_ERADATMMIN,RADNT_RHODONUT);



    // TOTRY: Om not happening even if set!

    //    lim[1]=lim[2]=lim[3]=MINM; // too low order for ~100 points
    //    if(WHICHPROBLEM==RADDONUT) lim[1]=lim[2]=lim[3]=PARALINE; // try later
    //    cour=0.5;

    if(!ISSPCMCOORDNATIVE(MCOORD) && (WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK || WHICHPROBLEM==RADDONUT) ){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADNT,\n",MCOORD);
      myexit(3434628752);
    }
    else if(WHICHPROBLEM==RADCYLBEAM && MCOORD!=CYLMINKMETRIC){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be CYLMINKMETRIC for RADCYLBEAM.\n",MCOORD);
      myexit(2493434634);
    }
    else if(WHICHPROBLEM==RADCYLBEAMCART && MCOORD!=CARTMINKMETRIC2){
      dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be CARTMINKMETRIC2 for RADCYLBEAMCART.\n",MCOORD);
      myexit(2493434635);
    }

    if(DOWALDDEN) cooling=NOCOOLING;
    else cooling=KORAL;

    // Danilo tilted
    //    cooling=COOLREBECCATHINDISK;
    cooling=COOLUSER;


    // ARAD_CODE=ARAD_CODE_DEF*1E5; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.
    //    GAMMAMAXRAD=1000.0; // Koral limits for this problem.
    GAMMAMAXRAD=50.0L; // Koral limits for this problem.
    GAMMAMAXRADFAIL=50.0L;
    GAMMAMAX=15.0L; // MHD


    ////////////
    //
    // BOUNDARY CONDITIONS

    if(WHICHPROBLEM==RADCYLBEAM){
      //      BCtype[X1DN]=ASYMM;
      //      BCtype[X1DN]=SYMM;
      BCtype[X1DN]=CYLAXIS;
      BCtype[X1UP]=RADCYLBEAMBC;
      BCtype[X2DN]=OUTFLOW;
      BCtype[X2UP]=OUTFLOW;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    else if(WHICHPROBLEM==RADCYLBEAMCART){
      BCtype[X1DN]=RADCYLBEAMCARTBC;
      BCtype[X1UP]=RADCYLBEAMCARTBC;
      BCtype[X2DN]=RADCYLBEAMCARTBC;
      BCtype[X2UP]=RADCYLBEAMCARTBC;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    else{

      BCtype[X1DN]=HORIZONOUTFLOW; // although more specific extrapolation based upon solution might work better

      if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK) BCtype[X1UP]=RADNTBC; // inflow analytic
      else if(WHICHPROBLEM==RADDONUT){
        //        BCtype[X1UP]=RADNTBC; // inflow analytic
        //else BCtype[X1UP]=FIXEDUSEPANALYTIC; // fixed analytic // little silly for most of outer boundary, so avoid // KORALTODO: Also causes hellish problems with solution and implicit solver at the X1UP boundary surface (not just near torus)
        BCtype[X1UP]=HORIZONOUTFLOW;
        //      BCtype[X1DN]=OUTFLOW;
        //      BCtype[X1UP]=OUTFLOW;
        
        //        if(DOWALDDEN) BCtype[X1UP]=OUTFLOW;
      }
      
      if(WHICHPROBLEM==RADFLATDISK)  BCtype[X2DN]=ASYMM; // if non-zero Rin_array[2]
      else if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADDONUT){
        BCtype[X2DN]=POLARAXIS; // assumes Rin_array[2]=0
      }


      
      if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK) BCtype[X2UP]=RADNTBC; // disk condition (with ASYMM done first)
      else  if(WHICHPROBLEM==RADDONUT){
        //        BCtype[X2UP]=ASYMM; // with donut, let free, so ASYMM condition across equator (hemisphere)
        BCtype[X2UP]=POLARAXIS; // assumes Rin_array[2]=pi (full sphere)
      }

      if(DOWALDDEN==2) BCtype[X2UP]=WALDMONOBC;

      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }




    ////////////
    // DUMP PERIODS

    int idt;
    if(WHICHPROBLEM==RADCYLBEAMCART){
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=10.0;
    }
    else{
      //      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1.0;
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=4.0;
      //for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    }
    
    if(totalsize[1]<=64){
      DTr = 100; //number of time steps for restart dumps
    }
    else{
      DTr = 1000;
    }

    if(WHICHPROBLEM==RADDONUT && totalsize[1]>64){
      // then, not testing, so full production mode for dumps
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=50.0;
      DTdumpgen[FIELDLINEDUMPTYPE]=4.0;
      DTdumpgen[IMAGEDUMPTYPE]=4.0;
      DTr=5000;
      if(PRODUCTION>=2){
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=100.0;
        DTdumpgen[DEBUGDUMPTYPE]=400.0;
        DTdumpgen[FIELDLINEDUMPTYPE]=4.0;
        DTdumpgen[IMAGEDUMPTYPE]=4.0;
        DTdumpgen[ENERDUMPTYPE]=4.0;
        DTr=5000;
      }
    }



    // tf = 100*DTdumpgen[0]; // 100 dumps(?)
    //    tf = 2000*DTdumpgen[0]; // koral in default setup does 1000 dumps
    tf = 1e5; //Danilo-Time

    if(DOWALDDEN) tf=1e5;

    //    DODIAGEVERYSUBSTEP = 1;

    //    if(WHICHPROBLEM==RADDONUT) DODIAGEVERYSUBSTEP = 1;

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADCYLJET){
    
    TIMEORDER=2;
    //    lim[1]=lim[2]=lim[3]=MC;
    //    lim[1]=lim[2]=lim[3]=DONOR;
    lim[1]=lim[2]=lim[3]=MINM;
    //lim[1]=lim[2]=lim[3]=MC;
    //  lim[1]=lim[2]=lim[3]=WENO5BND;
    //lim[1]=lim[2]=lim[3]=PARAFLAT;
    //lim[1]=lim[2]=lim[3]=PARALINE;
    
    cour=0.1;

    //    int set_fieldtype(void);
    //    int FIELDTYPE=set_fieldtype();

    BSQORHOLIMIT=1E2;
    BSQOULIMIT=1E8;
    cooling=KORAL;

    //gam=gamideal=5.0/3.0; // 4/3 if pairs
    gam=gamideal=4.0/3.0; // 4/3 if pairs or to compare with radiation-dominated case

    GAMMAMAXRAD=50.0L; // increase for higher jet speeds
    GAMMAMAXRADFAIL=50.0L; // increase for higher jet speeds
    GAMMAMAX=15.0L; // increase for higher jet speeds


    // 1: original Jane version
    // 2: pure jet version with reflective edge, uniform within, and radiative flux at boundary
    // 3: have high density in boundary that is fixed.
    // 4: hot jet but otherwise like #1
    // 5: rigid wall with panalytic setting of boundary conditions as flux, but outflow effectively
    // 6: like #1, but 2D and with absorption opacity so T matters
    RADCYLJET_TYPE=6;


    if(RADCYLJET_TYPE==6){
      //lim[1]=lim[2]=lim[3]=MINM;
      lim[1]=lim[2]=lim[3]=PARALINE;
    }
    else{
      lim[1]=lim[2]=lim[3]=MINM;
    }

    ////////////
    //
    // BOUNDARY CONDITIONS

    if(RADCYLJET_TYPE==1 || RADCYLJET_TYPE==4){
      //      BCtype[X1DN]=ASYMM;
      //      BCtype[X1DN]=SYMM;
      BCtype[X1DN]=CYLAXIS;
      BCtype[X1UP]=OUTFLOW;
      BCtype[X2DN]=OUTFLOW;
      BCtype[X2UP]=OUTFLOW;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    else if(RADCYLJET_TYPE==2||RADCYLJET_TYPE==3){
      //      BCtype[X1DN]=ASYMM;
      //      BCtype[X1DN]=SYMM;
      BCtype[X1DN]=CYLAXIS;
      BCtype[X1UP]=RADCYLJETBC;
      BCtype[X2DN]=OUTFLOW;
      BCtype[X2UP]=OUTFLOW;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    else if(RADCYLJET_TYPE==5){
      BCtype[X1DN]=CYLAXIS;
      BCtype[X1UP]=FREEOUTFLOW; //FIXEDUSEPANALYTIC;
      BCtype[X2DN]=OUTFLOW;
      BCtype[X2UP]=OUTFLOW;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    else if(RADCYLJET_TYPE==6){
      //      BCtype[X1DN]=ASYMM;
      //      BCtype[X1DN]=SYMM;
      BCtype[X1DN]=CYLAXIS;
      BCtype[X1UP]=OUTFLOW;
      BCtype[X2DN]=RADCYLJETBC;
      BCtype[X2UP]=OUTFLOW;
      BCtype[X3UP]=PERIODIC;
      BCtype[X3DN]=PERIODIC;
    }
    
    ////////////
    // DUMP PERIODS

    int idt;
    //      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1.0;
    //    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.01;
    //for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1.0;
    //for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=10.0;
    
   // tf = 200000;
      tf = 1e5; //Danilo-Time
    //    DODIAGEVERYSUBSTEP = 1;

    //    if(WHICHPROBLEM==RADDONUT) DODIAGEVERYSUBSTEP = 1;

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  return(0);

}


int init_defcoord(void)
{

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


  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==FLATNESS){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;
  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR){
   // a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=-50.0; 
    Rin_array[2]=-1.0;
    Rin_array[3]=-1.0;

    Rout_array[1]=50.0; 
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADPULSE3D){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=-50.0; 
    Rin_array[2]=-50.0; 
    Rin_array[3]=-50.0; 

    Rout_array[1]=50.0; 
    Rout_array[2]=50.0; 
    Rout_array[3]=50.0; 
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAMFLAT){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADTUBE){
   // a=0.0; // no spin in case use MCOORD=KSCOORDS

    if(NTUBE==5){
      defcoord = UNIFORMCOORDS;
      Rin_array[1]=-20.0;
      Rin_array[2]=-1.0;
      Rin_array[3]=-1.0;

      Rout_array[1]=20.0;
      Rout_array[2]=1.0;
      Rout_array[3]=1.0;
    }
    else{
      defcoord = UNIFORMCOORDS;
      Rin_array[1]=-15.0;
      Rin_array[2]=-1.0;
      Rin_array[3]=-1.0;

      Rout_array[1]=15.0;
      Rout_array[2]=1.0;
      Rout_array[3]=1.0;
    }
  }
  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADSHADOW){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=-1.0;
    Rin_array[2]=-1.0;
    Rin_array[3]=-1.0;
 
    Rout_array[1]=3.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADDBLSHADOW){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=-6.0;
    Rin_array[2]=0.0;
    Rin_array[3]=-1.0;
 
    Rout_array[1]=3.0;
    Rout_array[2]=1.5;
    Rout_array[3]=1.0;
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADBEAM2DKS){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;
  


    if (RADBEAM2D_BEAMNO==1){
      Rin_array[1]=2.6;
      Rout_array[1]=3.5;
    }
    else if (RADBEAM2D_BEAMNO==2){
      Rin_array[1]=5.5;
      Rout_array[1]=7.5;
    }
    else if (RADBEAM2D_BEAMNO==3){
      Rin_array[1]=14.5;
      Rout_array[1]=20.5;
    }
    else if (RADBEAM2D_BEAMNO==4){
      Rin_array[1]=30;
      Rout_array[1]=50;
    }

    Rin_array[2]=0.99*M_PI*0.5;
    Rout_array[2]=1.01*M_PI*0.5;

    Rin_array[3]=0.0;
    //    Rout_array[3]=M_PI*0.25;
    Rout_array[3]=M_PI*0.5;

 
  }
  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAM2DKSVERT){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;
  


    if (RADBEAM2DKSVERT_BEAMNO==1){
      Rin_array[1]=2.3;
      Rout_array[1]=3.5;
    }
    else if (RADBEAM2DKSVERT_BEAMNO==2){
      Rin_array[1]=5.5;
      Rout_array[1]=12.5;
    }
    else if (RADBEAM2DKSVERT_BEAMNO==3){
      // for MKS
      //Rin_array[1]=2.5;
      //   Rout_array[1]=3.0;
      // for UNI
      Rin_array[1]=14.5;
      Rout_array[1]=20.5;
    }
    else if (RADBEAM2DKSVERT_BEAMNO==4){
      Rin_array[1]=30;
      Rout_array[1]=50;
    }
    else if (RADBEAM2DKSVERT_BEAMNO==5){
      Rin_array[1]=5.5;
      Rout_array[1]=12.5;
    }

    // NOTE: For testing radiation through polar axis region with axis *on* grid.
    Rin_array[2]= -0.25*Pi/2.;
    Rout_array[2]=+0.27*Pi/2.;

    Rin_array[3]=0.;
    Rout_array[3]=0.01*Pi/4.; // 8.*Pi/4.
      
 
  }
  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==ATMSTATIC){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
 
    Rin_array[1]=1E6;
    Rout_array[1]=2E6;

    Rin_array[2]=0.99*M_PI*0.5;
    Rout_array[2]=1.01*M_PI*0.5;

    Rin_array[3]=-1.0;
    Rout_array[3]=1.0;
 
  }
  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADATM){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
 
    Rin_array[1]=1E6;
    Rout_array[1]=1.4E6;

    Rin_array[2]=0.99*M_PI*0.5;
    Rout_array[2]=1.01*M_PI*0.5;

    Rin_array[3]=-1.0;
    Rout_array[3]=1.0;
 
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADWALL){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=-6;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=3.0;
    Rout_array[2]=1.5;
    Rout_array[3]=1.0;
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADWAVE){
    //a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
 #if(WHICHPROBLEM==KOMIPROBLEM)
    FTYPE xl, xc, xr;
   // a=0.0; // no spin in case use MCOORD=KSCOORDS
    
    defcoord = UNIFORMCOORDS;
    
    //fast shock
    if(WHICHKOMI==1){
      xl = -1.0;
      xc =  0.0;
      xr =  1.0;
    }
    //slow shock
    else if(WHICHKOMI==2){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    //fast switch-off rarefaction
    else if(WHICHKOMI==3){
      xl = -1.0;
      xc =  0.0;
      xr =  1.0;
    }
    //slow switch-on rarefaction
    else if(WHICHKOMI==4){
      xl = -0.75;
      xc =  0.0;
      xr =  1.2;
    }
    //alfven wave
    else if(WHICHKOMI==5){
      xl = -1.0;
      xc =  0.0;
      xr =  1.5;
    }
    //compound wave (plot is -0.5 to 1.5, but need 4 cells between -0.025 and 0.0 according to Figure 5)
    else if(WHICHKOMI==6){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    //Shock tube 1
    else if(WHICHKOMI==7){
      xl = -1.0;
      xc =  0.0;
      xr =  1.5;
    }
    //Shock tube 2
    else if(WHICHKOMI==8){
      xl = -1.2;
      xc =  0.0;
      xr =  1.2;
    }
    //Collision
    else if(WHICHKOMI==9){
      xl = -1.0;
      xc =  0.0;
      xr =  1.0;
    }


    if(WHICHKOMI==101){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==102){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==103){
      xl = -1.5;
      xc =  0.0;
      xr =  0.5;
    }
    if(WHICHKOMI==104){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==105){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==106){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==107){
      xl = -0.5;
      xc =  0.0;
      xr =  1.5;
    }
    if(WHICHKOMI==108){
      xl = -1.0;
      xc =  0.0;
      xr =  2.0;
    }
    if(WHICHKOMI==109){
      xl = -1.0;
      xc =  0.0;
      xr =  2.0;
    }

    Rin = Rin_array[1]=xl;
    Rin_array[2]=0;
    Rin_array[3]=0;
    
    Rout = Rout_array[1]=xr;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;

#endif


  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBONDI){
 //   a=0.0; // no spin in case use MCOORD=KSCOORDS

    // TOTRY: Change horizon interpolation to be like koral.
    // PURE HYDRO: entropy inversions don't lead to major problems.
    // TRY pure HD.
    // TRY moving the inner boundary outward a bit.

    // KORALTODO: 2.5 to 2E4 in paper
    //    RADBONDI_MINX=3.5;
    //    RADBONDI_MAXX=2e3;

    if(MCOORD==KSCOORDS){
      RADBONDI_MINX=1.9; // few cells inside horizon // 1.9 good for N1=512
      RADBONDI_MAXX=2e4;
    }
    else{
      RADBONDI_MINX=2.5;
      RADBONDI_MAXX=2e4;
    }

    //#define LOGXGRID
    //    FTYPE LOGPAR1=2.2;
    //    FTYPE LOGPAR2=2.;

    // defcoord = UNIFORMCOORDS;
    defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions


    if(RADBONDI_TESTNO==1){
      if(MCOORD==KSCOORDS){
        RADBONDI_MINX=1.9; // few cells inside horizon // 1.9 good for N1=512
      }
      else{
        RADBONDI_MINX=2.5;
      }
      R0=(1.85/1.9)*RADBONDI_MINX; // due to cold temp, energy equation evolution requires (with PARA) more cells near BH.  Something like MP5 (that does some avg->point conversions) as in Koral wouldn't need to do that.
      // also may require not so large tolerance, but so far not limiting solution like energy equation does.
    }
    else{
      // normal MINX ok.
      R0=0.9*RADBONDI_MINX;
      // can reduce noise in u_g in solutions if use larger R0.
    }


    Rin_array[2]=.99*Pi/2.;
    Rout_array[2]=1.01*Pi/2.;
    Rin_array[3]=-1.;
    Rout_array[3]=1.;


    if(0){
      ///////////
      // 2D MAGBONDI TEST
      RADBONDI_MINX=1.9;
      RADBONDI_MAXX=2e4;
      //    RADBONDI_MAXX=2e2;
      R0=1.88;
      Rin_array[2]=0.0;
      Rout_array[2]=Pi;
    }

    if(0){
      if(RADBONDI_TESTNO==3){
        // normal 1D non-MAG BONDI=3 test that ran for paper
        RADBONDI_MINX=2.5;
        RADBONDI_MAXX=2e4;
        R0=2.2;
        Rin_array[2]=.99*Pi/2.;
        Rout_array[2]=1.01*Pi/2.;
      }
    }


    //////////// nothing else to set for RADBONDI

    Rin=RADBONDI_MINX;
    Rout=RADBONDI_MAXX;


    if(R0>=RADBONDI_MINX){
      dualfprintf(fail_file,"Must have R0=%g < RADBONDI_MINX=%g\n",R0,RADBONDI_MINX);
      myexit(243532469);
    }

    Rhor=rhor_calc(0);
    if(RADBONDI_MINX>Rhor && MCOORD==KSCOORDS){
      dualfprintf(fail_file,"WARNING: Have boundary oustide horizon in KS coords\n");
    }


  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADDOT){
  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    defcoord = UNIFORMCOORDS;
    Rin_array[1]=0;
    Rin_array[2]=0;
    Rin_array[3]=0;

    Rout_array[1]=1.0;
    Rout_array[2]=1.0;
    Rout_array[3]=1.0;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADNT){
  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    if(1){
      //      RADNT_MINX=1.7; // allowed in KSCOORDS
      RADNT_MINX=1.9; // allowed in KSCOORDS (latest koral)
      RADNT_MAXX=50.0;
    }
    else{
      RADNT_MINX=1.5*Rhor;
      RADNT_MAXX=40.0; // 27.8
    }

    // KORALTODO: Why doesn't koral just use same log coords as used for RADBONDI?  Change of feature set.
    // defcoord = UNIFORMCOORDS;
    defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
    R0=0.0;
    Rin=RADNT_MINX;
    Rout=RADNT_MAXX;

    //    Rin_array[2]=0.2*Pi/4.;
    Rin_array[2]=0.0*Pi/4.;
    Rout_array[2]=Pi/2.;
    Rin_array[3]=-1.;
    Rout_array[3]=1.;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADFLATDISK){
  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    RADNT_MINX=4.0;
    RADNT_MAXX=100.0;

    if(0){
      // KORALTODO: The below use of log radial grid leads to opacity related failure of problem when using pre-recent koral kappa (i.e. not zero)
      defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
      R0=0.0;
      Rin=RADNT_MINX;
      Rout=RADNT_MAXX;
    }
    else{
      defcoord = UNIFORMCOORDS;
      Rin_array[1]=RADNT_MINX;
      Rout_array[1]=RADNT_MAXX;
    }


    //    Rin_array[2]=0.0*Pi/4.;
    //    Rin_array[2]=0.01*Pi/4.;
    Rin_array[2]=0.1*Pi/4.;
    Rout_array[2]=Pi/2.;
    Rin_array[3]=-1.;
    Rout_array[3]=1.;

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADDONUT){
    //    a=0.0; // no spin in case use MCOORD=KSCOORDS

    // metric stuff first
    a = 0.5;  // WALD
    

    if(1){
      RADNT_MINX=1.7; // allows in KSCOORDS
      //      RADNT_MAXX=50.0;
      //      RADNT_MAXX=60.0;
      //      RADNT_MAXX=400.0; // what was using before, but problems at outer radial edge develop after t\sim 5600
      RADNT_MAXX=1E4;
    }
    else{
      RADNT_MINX=1.8*Rhor;
      RADNT_MAXX=40.0; // 27.8
    }

    // KORALTODO: Why doesn't koral just use same log coords as used for RADBONDI?
    // defcoord = UNIFORMCOORDS;
    //    defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
    Rin=RADNT_MINX;
    Rout=RADNT_MAXX;

    if(DOWALDDEN==0) defcoord=JET6COORDS;
    else if(DOWALDDEN){
      //      Rout=2000.0; // new normal
      Rout=400.0; // newer normal
      defcoord=USERCOORD;
    }
    //    else if(DOWALDDEN==2) defcoord=JET6COORDS;


    if(0){
      Rin_array[1]=Rin;
      Rout_array[1]=Rout;
      Rin_array[2]=0.0*Pi/4.; // but koral currently uses 0.5*Pi/4
      Rout_array[2]=Pi; // KORALNOTE: Different from KORAL code test
      Rin_array[3]=-1.;
      Rout_array[3]=1.;
      defcoord=LOGRUNITH;
    }


    Rhor=rhor_calc(0);


    //  hslope = 0.3;
    hslope = 1.04*pow(h_over_r,2.0/3.0);
    // NOTEMARK: Should change h0 in coord.c from h0=0.3 to h0=0.1 or something for thin disks


    if(DOWALDDEN){ // to give more resolution near BH
      R0=0.7;
    }
    else{
      //    R0=0.0;
      R0=0.2;
    }



    if(Rout<1E3){
      Rin=1.05;
    }
    if(a==0.5){
      Rin=1.3; // a = 0.5 for totalsize[1]=128
      //      setRin_withchecks(&Rin);
    }
    if(a==0.8){
      Rin=1.2; // a = 0.8 for totalsize[1]=128
      //      setRin_withchecks(&Rin);
    }

    if(DOWALDDEN && a<0.9 && totalsize[1]>64) Rin=1.4;
    
    if(totalsize[1]<32*4&&DOWALDDEN==0){
      dualfprintf(fail_file,"RADDONUT setup for 128x64 with that grid\n");
      //myexit(28634693);
    }

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADCYLBEAM){
 //   a=0.0; // no spin in case use MCOORD=KSCOORDS

    if(WHICHPROBLEM==RADCYLBEAM) MBH=0.0; // because CYLMINKMETRIC really has gravity if choose MBH!=0

    // TOTRY: find azimuthal flux.

    RADNT_FULLPHI=(Pi/2.5);
    //    RADNT_FULLPHI=(2.0*Pi);

    RADNT_MINX=0.0; // all the way to the R=0 origin
    RADNT_MAXX=20.0;
 
    if(1){
      defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
      R0=-1.0;
      Rin=RADNT_MINX;
      Rout=RADNT_MAXX;
      
      Rin_array[2]=-1.0;
      Rout_array[2]=1.0;

      Rin_array[3]=0.0;
      Rout_array[3]=RADNT_FULLPHI;
    }
    else{
      defcoord = UNIFORMCOORDS;
     
      Rin_array[1]=RADNT_MINX;
      Rout_array[1]=RADNT_MAXX;

      Rin_array[2]=-1.0;
      Rout_array[2]=1.0;

      Rin_array[3]=0.0;
      Rout_array[3]=RADNT_FULLPHI;
    }

  }

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADCYLBEAMCART){
  //  a=0.0;

    RADNT_MINX=-20.0;
    RADNT_MAXX=+20.0;
 
    defcoord = UNIFORMCOORDS;
     
    Rin_array[1]=RADNT_MINX;
    Rout_array[1]=RADNT_MAXX;

    Rin_array[2]=RADNT_MINX;
    Rout_array[2]=RADNT_MAXX;

    Rin_array[3]=0.0;
    Rout_array[3]=1.0;

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADCYLJET){
  //  a=0.0; // no spin in case use MCOORD=KSCOORDS

    MBH=0.0; // because CYLMINKMETRIC really has gravity if choose MBH!=0

    // TOTRY: find azimuthal flux.

    //    RADNT_FULLPHI=(Pi/2.5);
    RADNT_FULLPHI=(2.0*Pi);

    if(RADCYLJET_TYPE==1 || RADCYLJET_TYPE==4){
      RADNT_MINX=1E-2; // all the way to the R=0 origin
      RADNT_MAXX=10.0;
    }
    else if(RADCYLJET_TYPE==2||RADCYLJET_TYPE==3){
      RADNT_MINX=1E-2; // all the way to the R=0 origin
      RADNT_MAXX=1.0;
    }
    else if(RADCYLJET_TYPE==5){
      RADNT_MINX=1E-2; // all the way to the R=0 origin
      RADNT_MAXX=1.0;
    }
    else if(RADCYLJET_TYPE==6){
      RADNT_MINX=1E-2; // all the way to the R=0 origin
      RADNT_MAXX=15.0;
    }
 
    if(RADCYLJET_TYPE==1 || RADCYLJET_TYPE==4 ||RADCYLJET_TYPE==2||RADCYLJET_TYPE==3){
      defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
      R0=-1.0;
      Rin=RADNT_MINX;
      Rout=RADNT_MAXX;
      Rin_array[1]=RADNT_MINX;
      Rout_array[1]=RADNT_MAXX;
      
      Rin_array[2]=-1.0;
      Rout_array[2]=1.0;

      Rin_array[3]=0.0;
      Rout_array[3]=RADNT_FULLPHI;
    }
    if(RADCYLJET_TYPE==5){
      defcoord = UNIFORMCOORDS;
     
      Rin_array[1]=RADNT_MINX;
      Rout_array[1]=RADNT_MAXX;

      Rin_array[2]=-1.0;
      Rout_array[2]=1.0;

      Rin_array[3]=0.0;
      Rout_array[3]=RADNT_FULLPHI;
    }
    if(RADCYLJET_TYPE==6){
      defcoord = LOGRUNITH; // Uses R0, Rin, Rout and Rin_array,Rout_array for 2,3 directions
      R0=-0.2;
      Rin=RADNT_MINX;
      Rout=RADNT_MAXX;
      Rin_array[1]=RADNT_MINX;
      Rout_array[1]=RADNT_MAXX;
      
      Rin_array[2]=0.0;
      Rout_array[2]=30.0*45.0; //Rout*10.0;

      Rin_array[3]=0.0;
      Rout_array[3]=RADNT_FULLPHI;
    }

  }


  /*************************************************/
  /*************************************************/
  /*************************************************/


  return(0);
}


int init_grid(void)
{
  
  init_defcoord(); // just avoids splitting function call, here sets R0,Rin,Rout

  return(0);
}



// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int funreturn;

  if(WHICHPROBLEM==RADDONUT){
    if(0&&(RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3)){ // keep as 0&& because not setting E properly in init_atmosphere()
      funreturn=user1_init_atmosphere(whichvel, whichcoord,i, j, k, pr);
      if(funreturn!=0) return(funreturn);
      return(0);
    }
    else{
      // NO atmosphere
      // tells to do no coordinate transformations
      *whichvel=WHICHVEL;
      *whichcoord=PRIMECOORDS;
      return(-1);
    }
  }
  else{
    // NO atmosphere
    // tells to do no coordinate transformations
    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(-1);
  }


  return(-1); // no atmosphere set, so don't do anything at all

}

int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);




  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);



  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();
  


  // defaults
  randfact = 0.1;
  beta=100.0;
  rin=6.0;
  rinfield=rin;
  routfield=12.0;

  if(WHICHPROBLEM==RADDONUT){

    if(FIELDTYPE==DISK2FIELD){
      //rin=6.0; // old setting
      rinfield=rin=9.0;
      routfield=13.0;
      //  beta=100.0; // was used for rada0.94 etc.
      beta = 10.0;
    }
    else if(FIELDTYPE==FIELDJONMAD){
      rin=9.0;
      rinfield=12;
      beta = 30.0;
    }
    else if(FIELDTYPE==FIELDWALD || FIELDTYPE==MONOPOLE){
      rin=9.0;
      rinfield=1.0;
      routfield=2.0;
      // so field at horizon is BSQORHOWALD
      beta = ((gam-1.0)*RADNT_UINTATMMIN)/(RADNT_RHOATMMIN*BSQORHOWALD*0.5);
    }

    if(RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      //      rinfield=1.1*Risco;
      rinfield=Risco;
      //      beta=1E30;
      beta=10.0;

      if(FIELDTYPE==FIELDJONMAD) rin=0.0;
      else rin=rinfield;
    }

    // SUPERMADNEW
    if(RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      rin=1.1*Risco;
      rinfield=1.1*Risco;
 
      rin=0.0;
 
      rinfield=20.0; // works for any a for these models
      routfield=40.0;

      //rinfield=Risco;
      //      beta=1E30;
      if(FIELDTYPE==FIELDJONMAD){
        beta=10.0; // so nearly MAD at t=0
      }
      else{
        beta=5.0; // so perturbations are at 2 but rest less and tends to be beta\sim 20
      }
    }


  }





  if(WHICHPROBLEM==RADDOT){
    RADDOT_XDOT=(20.0/41.0)*(Rout_array[1]-Rin_array[1]) + Rin_array[1];
    RADDOT_YDOT=(10.0/41.0)*(Rout_array[2]-Rin_array[2]) + Rin_array[2];
    RADDOT_ZDOT=(20.0/41.0)*(Rout_array[3]-Rin_array[3]) + Rin_array[3];

    // get X1,X2,X3 of dot assuming UNIFORMCOORDS
    FTYPE myX[NDIM]={0.0};
    FTYPE dxdxp[NDIM][NDIM];
    dxdxprim_ijk(0, 0, 0, CENT, dxdxp);
    myX[1]= startx[1] + (RADDOT_XDOT-Rin_array[1])/dxdxp[1][1];
    myX[2]= startx[2] + (RADDOT_YDOT-Rin_array[2])/dxdxp[2][2];
    myX[3]= startx[3] + (RADDOT_ZDOT-Rin_array[3])/dxdxp[3][3];

    trifprintf("RADDOT: myX: %g %g %g\n",myX[1],myX[2],myX[3]);

    // get nearest i,j,k
    extern void icoord_round(FTYPE *X,int loc, int *i, int *j, int *k);
    icoord_round(myX,CENT,&RADDOT_IDOT,&RADDOT_JDOT,&RADDOT_KDOT);

    RADDOT_FYDOT=0.3;
    RADDOT_LTEFACTOR=1.;
    RADDOT_URFX=0.;
    RADDOT_F1=100;
    RADDOT_F2=10000;
    
    trifprintf("RADDOT: %g %g %g : %d %d %d\n",RADDOT_XDOT,RADDOT_YDOT,RADDOT_ZDOT,RADDOT_IDOT,RADDOT_JDOT,RADDOT_KDOT);
  }

  if(ISBLACKHOLEMCOORD(MCOORD)){
    trifprintf("BEGIN check_rmin\n");
    // check rmin
    check_rmin();
    trifprintf("END check_rmin\n");
  }

  if(1){
    // check that singularities are properly represented by code
    trifprintf("BEGIN check_spc_singularities_user\n");
    // SUPERGODMARK: Goes very slowly sometimes randomly for unknown reasons.
    dualfprintf(fail_file,"WARNING: check_spc_singularities_user() oddly stalls sometimes...\n");
    check_spc_singularities_user();
    trifprintf("END check_spc_singularities_user\n");
    dualfprintf(fail_file,"WARNING: done with check_spc_singularities_user(), but it sometimes stalls or goes very very slow for no apparently good reason.  E.g., on NAUTILUS with -O0, very slow checks.  But just putting dualfprintf before and after the above call leads to quick finish.\n");
  }
  
  return(0);

}



//****************************************//
//****************************************//
//****************************************//
// Problem setup constants that only modifies things in init.c (not init.h)
//****************************************//
//****************************************//

#if(WHICHPROBLEM==FLATNESS)

#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif

//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D)




#if(WHICHPROBLEM==RADPULSEPLANAR)

#define KAPPA 0.
//#define KAPPAES (0.0)
//#define KAPPAES (1E-7)
//#define KAPPAES (1E-4*1.09713E-18*1E6)
//#define KAPPAES (1E-4*1.09713E-18*1E3)
//#define KAPPAES (1E-4*1.09713E-18*1E-0)
//#define KAPPAES (1E-4*1.09713E-18*0.2)
//#define KAPPAES (1E-4*1.09713E-18*1E-1)
//#define KAPPAES (1E-4*1.09713E-18*1E-2)
//#define KAPPAES (1E-4*1.09713E-18*1E-3)
//#define KAPPAES (1E-4*1.09713E-18*1E-3*1E-10)

// assume KAPPA defines fraction of FF opacity
//#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
//#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))

#define KAPPAES (1E3) // takes VERY long time with sub-cycling, but works.
//#define KAPPAES (1E2) // goes with sub-cycling at "ok" rate for this test.
//#define KAPPAES (10.0)
//#define KAPPAES (1E-1)
//#define KAPPAES (1.0)
//#define KAPPAES (1E-10)

#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA)
#define KAPPAESUSER(rho,T) (rho*KAPPAES)


#else // PULSE and PULSE3D

// KAPPAs are fraction of physical FF and ES opacities
#define KAPPA 0.
#define KAPPAES (SMALL)

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif


#endif


//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADBEAMFLAT)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif

//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADTUBE)

#define KAPPAESUSER(rho,T) (0.0)

#if(NTUBE==1)
#define KAPPAUSER(rho,B,Tg,Tr) (0.4*rho)
#elif(NTUBE==2)
#define KAPPAUSER(rho,B,Tg,Tr) (0.2*rho)
#elif(NTUBE==3)
#define KAPPAUSER(rho,B,Tg,Tr) (0.3*rho)
#elif(NTUBE==31)
#define KAPPAUSER(rho,B,Tg,Tr) (25*rho)
#elif(NTUBE==4)
#define KAPPAUSER(rho,B,Tg,Tr) (0.08*rho)
#elif(NTUBE==41)
#define KAPPAUSER(rho,B,Tg,Tr) (0.7*rho)
#elif(NTUBE==5)
#define KAPPAUSER(rho,B,Tg,Tr) (1000*rho)
#endif


#endif

//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADSHADOW || WHICHPROBLEM==RADDBLSHADOW)

//#define KAPPAUSER(rho,B,Tg,Tr) (rho*1E2)
#define KAPPAUSER(rho,B,Tg,Tr) (rho*1.0) // paper
#define KAPPAESUSER(rho,T) (rho*0.0)


#endif


//****************************************//
//****************************************//


//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADBEAM2DKS || WHICHPROBLEM==RADBEAM2DKSVERT)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif


#if(WHICHPROBLEM==ATMSTATIC)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif


#if(WHICHPROBLEM==RADATM)


#define KAPPA 0.
#define KAPPAES 1. // only scattering

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif


#if(WHICHPROBLEM==RADWALL)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif




#if(WHICHPROBLEM==RADWAVE)

#define KAPPAUSER(rho,B,Tg,Tr) (rho*RADWAVE_KAPPA)
#define KAPPAESUSER(rho,T) (rho*RADWAVE_KAPPAES)

#endif



#if(WHICHPROBLEM==KOMIPROBLEM)

#define KAPPAUSER(rho,B,Tg,Tr) (0.)
#define KAPPAESUSER(rho,T) (0.)

#endif


#if(WHICHPROBLEM==RADBONDI)


#define KAPPA 1.0
#define KAPPAES 1.0

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif


#if(WHICHPROBLEM==RADDOT)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg,Tr))
// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_BASIC_CODE(rho,T))


#endif

#if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK)

#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA_ES_CODE(rho,Tg)/1E14*0.1) // wierd use of kappa_{es} in koral
//#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA_ES_CODE(rho,T)/1E14*0.0)
#define KAPPAESUSER(rho,T) (0.0)

#endif

#if(WHICHPROBLEM==RADDONUT)
// kappa can't be zero or else flux will be nan
//#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA_ES_CODE(rho,T)/1E14*1.0) // wierd use of kappa_{es} in koral
//#define KAPPAESUSER(rho,T) (0.0)

// KORALNOTE: Different than koral code test, but as if full problem.
#define KAPPA 1.0
#define KAPPAES 1.0

// KORALTODO: Put a lower limit on T~1E4K so not overly wrongly opaque in spots where u_g->0 anomologously?
// assume KAPPA defines fraction of FF opacity
//#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*KAPPA_FF_CODE(rho,Tg+TEMPMIN))
// accounts for low temperatures so non-divergent and more physical

#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA*(KAPPA_GENFF_CODE(SMALL+rho,Tg+TEMPMIN,Tr+TEMPMIN)+KAPPA_SYN_CODE(SMALL+B,Tg+TEMPMIN,Tr+TEMPMIN)))

// TODOMARK: Compute FF number opacity
#define KAPPANUSER(rho,B,Tg,Tr) (rho*KAPPA*(KAPPA_GENFF_CODE(SMALL+rho,Tg+TEMPMIN,Tr+TEMPMIN)+KAPPAN_SYN_CODE(SMALL+B,Tg+TEMPMIN,Tr+TEMPMIN)))

// assume KAPPAES defines fraction of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


#endif

#if(WHICHPROBLEM==RADCYLBEAM || WHICHPROBLEM==RADCYLBEAMCART)

#define KAPPAUSER(rho,B,Tg,Tr) (rho*KAPPA_ES_BASIC_CODE(rho,Tg)/1E14*0.0) // note 0.0
#define KAPPAESUSER(rho,T) (0.0)

#endif

#if(WHICHPROBLEM==RADCYLJET)

#define KAPPAUSER(rho,B,Tg,Tr) (rho*(KAPPA_FF_CODE(SMALL+rho,Tg+TEMPMIN,Tr+TEMPMIN))) // SMALL
#define KAPPAESUSER(rho,Tg) (rho*KAPPA_ES_BASIC_CODE(rho,Tg)/100.0)

#endif


#ifndef KAPPANUSER
// in case haven't defined KAPPANUSER, just use energy opacity
#define KAPPANUSER(rho,B,Tg,Tr) KAPPAUSER(rho,B,Tg,Tr)
#endif

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

// When setting primitives, put conditionals around PRAD? or URAD? variables to if want to be able to set EOMRADTYPE to EOMRADNONE and have work as non-radiation problem
int init_dsandvels(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_koral(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

  // assume inittype not used, pos==CENT, and time doesn't matter (e.g. only used at t=0)

  init_dsandvels_koral(whichvel, whichcoord,  i,  j,  k, pr, pstag);


  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();

  if(FIELDTYPE==SPLITMONOPOLE || FIELDTYPE==MONOPOLE || FIELDTYPE==FIELDWALD){
    // generate field with same whichcoord as init_dsandvels_koral() set
    int whichmethod=0;
    int whichinversion=0;
    fieldprim(whichmethod, whichinversion, whichvel, whichcoord, i, j, k, pr); // assume pstag set as average of pr
  }


  // assume any floor set at t=0 is part of solution
  if(DOYFL){
    FTYPE rhofloor=pr[RHO]*NUMEPSILON*10.0;
    FTYPE vfloor=NUMEPSILON*10.0;
    FTYPE enfloor=ERADLIMIT+(pr[URAD0]+pr[UU])*NUMEPSILON*10.0;
    if(YFL1>=0) pr[YFL1] = SMALL + rhofloor; // rho floor
    if(YFL2>=0) pr[YFL2] = SMALL + rhofloor*vfloor*vfloor; // -T^t_t-rho u^r floor
    if(YFL3>=0) pr[YFL3] = SMALL + rhofloor*vfloor; // T^t_phi floor
    if(YFL4>=0) pr[YFL4] = SMALL + enfloor; // -R^t_t floor
    if(YFL5>=0) pr[YFL5] = SMALL + enfloor*vfloor; // R^t_\phi floor
  }

  return(0);
}


// unnormalized density
int init_dsandvels_koral(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  FTYPE X[NDIM],V[NDIM];
  int pl,pliter;

  //  coord(i, j, k, CENT, X);
  //  bl_coord(X, V);
  //  r=V[1];
  //  th=V[2];

  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==FLATNESS){

    
    pr[RHO] = 1./RHOBAR ; // i.e. 1g/cm^3
    pr[UU] = 0.1/RHOBAR; // i.e. c^2 * 1g/cm^3 of energy density
    pr[U1] = 0 ;
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    if(PRAD0>=0){
      pr[URAD0] = 1./RHOBAR; // i.e. c^2 * 1g/cm^3 of energy density
      pr[URAD1] = 0 ;
      pr[URAD2] = 0 ;    
      pr[URAD3] = 0 ;
    }

    *whichvel=WHICHVEL;
    *whichcoord=CARTMINKMETRIC2;

    // KORALTODO: no transformation for radiation.  Would give same result as assuming in fluid frame because vfluid=0 here and F=0 here.

    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAMFLAT){

    
    pr[RHO] = RADBEAMFLAT_RHO ;
    pr[UU] = RADBEAMFLAT_UU;
    pr[U1] = 0 ;
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }



    if(PRAD0>=0){
      // new way: correctly transform -- also how koral currently setup
      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = RADBEAMFLAT_ERAD;
      pradffortho[PRAD1] = 0;
      pradffortho[PRAD2] = 0;
      pradffortho[PRAD3] = 0;


      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      *whichvel=VEL4;
      *whichcoord=MCOORD;
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL,pradffortho,pr, pr);
    }
    else if(PRAD0>=0){
      // old way: don't transform, leave as radiation frame E.
      pr[PRAD0] = RADBEAMFLAT_ERAD;
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;
      *whichvel=WHICHVEL;
      *whichcoord=MCOORD;
    }
    else{
      *whichvel=WHICHVEL;
      *whichcoord=MCOORD;
    }

    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D){
    // now avoids use of real units as in Koral paper (unlike Koral code)

    FTYPE Trad,Tgas,ERAD,uint;
    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];


    if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSE3D){
      rsq=(xx)*(xx)+(yy)*(yy)+(zz)*(zz);
    }
    else if(WHICHPROBLEM==RADPULSEPLANAR){
      rsq=(xx)*(xx);

    }

    //FTYPE RHO_AMB (1.e0) // in grams per cm^3
    // FTYPE RHO_AMB=(MPERSUN*MSUN/(LBAR*LBAR*LBAR)); // in grams per cm^3 to match koral's units with rho=1
    FTYPE RHO_AMB=1.0;
    FTYPE T_AMB=1.0E6/TEMPBAR;

    FTYPE BLOBP=100.;
    FTYPE BLOBW=5.;

    // radiation temperature is distributed
    Trad=T_AMB*(1.+BLOBP*exp(-rsq/(BLOBW*BLOBW)));
    ERAD=calc_LTE_EfromT(Trad);

    //flat gas profiles
    Tgas=T_AMB;
    FTYPE rho;
    rho=RHO_AMB;
    uint=calc_PEQ_ufromTrho(Tgas,rho);

    // dualfprintf(fail_file,"IC i=%d Trad=%g ERAD=%g Tgas=%g rho=%g uint=%g\n",i,Trad,ERAD,Tgas,rho,uint);

   
    pr[RHO] = rho;
    pr[UU] = uint;
    pr[U1] = 0 ;
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    FTYPE Fx,Fy,Fz;
    Fx=Fy=Fz=0.0;

    if(PRAD0>=0){
      pr[URAD0] = ERAD ;
      pr[URAD1] = Fx ;
      pr[URAD2] = Fy ;    
      pr[URAD3] = Fz ;
    }

    // KORALTODO: no transformation, but only because tuned units to be like koral and so ERAD gives same value and also because no Flux.   Also, would give same result as assuming in fluid frame because vfluid=0 here and F=0 here.

    *whichvel=WHICHVEL;
    *whichcoord=CARTMINKMETRIC2;
    return(0);
  }



  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADTUBE){

    FTYPE rho,mx,my,mz,m,ERAD,uint,E0,Fx,Fy,Fz,pLTE;  
    FTYPE rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut,ux;
    FTYPE xx,yy,zz;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];
  

    // set fluid values.  Also set radiation ff primitives (E,F^i), which are set in fluid frame orthonormal basis
    if(xx<(Rout_array[1]+Rin_array[1])/2.){
      rho=1.;
      if(NTUBE==1) {uint = 3.e-5 / (gamideal - 1.); ERAD=1.e-8; Fx=1.e-2*ERAD;ux=0.015;}
      if(NTUBE==2) {uint = 4.e-3 / (gamideal - 1.);ERAD=2.e-5; Fx=1.e-2*ERAD;ux=0.25;}
      if(NTUBE==3 || NTUBE==31) {uint = 60. / (gamideal - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=10.;}
      if(NTUBE==4 || NTUBE==41) {uint = 6.e-3 / (gamideal - 1.);ERAD=0.18; Fx=1.e-2*ERAD;ux=0.69;}   
      if(NTUBE==5) {uint = 60. / (gamideal - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=1.25;}
    }
    else{
      if(NTUBE==1) {rho=2.4;uint = 1.61e-4/ (gamideal - 1.); ERAD=2.51e-7; Fx=1.e-2*ERAD;ux=6.25e-3;}
      if(NTUBE==2) {rho=3.11;uint = 0.04512 / (gamideal - 1.);ERAD=3.46e-3; Fx=1.e-2*ERAD;ux=0.0804;}
      if(NTUBE==3 || NTUBE==31) {rho=8.0;uint = 2.34e3 / (gamideal - 1.);ERAD=1.14e3; Fx=1.e-2*ERAD;ux=1.25;}
      if(NTUBE==4 || NTUBE==41) {rho=3.65;uint =3.59e-2 / (gamideal - 1.);ERAD=1.30; Fx=1.e-2*ERAD;ux=0.189;}   
      if(NTUBE==5) {rho=1.0;uint = 60. / (gamideal - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=1.10;}
    }
    Fy=Fz=0.0;

  
    pr[RHO] = rho;
    pr[UU] = uint;
    pr[U1] = ux ; // ux is 4-velocity
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
  
    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL4;
    *whichcoord=CARTMINKMETRIC2;

    if(PRAD0>=0){
      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERAD;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;

      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL,pradffortho,pr, pr);

      //  PLOOPRADONLY(pl) dualfprintf(fail_file,"FOO1: i=%d pl=%d pr=%g\n",ptrgeomreal->i,pl,pr[pl]);

    }
  
    return(0);

  }


  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADSHADOW || WHICHPROBLEM==RADDBLSHADOW){

    // FTYPE MASS=10.0;
    FTYPE TAMB=1.e7/TEMPBAR;
    FTYPE RHOAMB=1.e-4;
    FTYPE RHOBLOB=1.e3;
    FTYPE BLOBW=0.22;

    FTYPE Trad,Tgas,ERAD;
    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    FTYPE rho,uint,Fx,Fy,Fz,pLTE;  

    /*****************************/
 
    // FTYPE pamb=calc_PEQ_ufromTrho(TAMB,RHOAMB);
    rsq=xx*xx+yy*yy+zz*zz;
    rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;      
    Tgas=TAMB*RHOAMB/rho;
    // Paper says T = T0*rho/RHOAMB
    // Tgas = TAMB*rho/RHOAMB;
    // for constant gas pressure, P=\rho T implies rho T = constant so that T\propto 1/rho
    uint=calc_PEQ_ufromTrho(Tgas,rho);

    Trad=TAMB;
    ERAD=calc_LTE_EfromT(Trad);

    // dualfprintf(fail_file,"i=%d j=%d rho=%g Trad=%g uint=%g ERAD=%g\n",i,j,rho,Trad,uint,ERAD);

    Fx=0.;
    Fy=0.;
    Fz=0.;              

    FTYPE VV=0.; // assumed orthonormal 4-velocity

    pr[RHO] = rho;
    pr[UU] = uint;
    pr[U1] = -VV ;
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL4;
    *whichcoord=CARTMINKMETRIC2;
    

    if(PRAD0>=0){
      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERAD;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;

      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL,pradffortho,pr, pr);
      // PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pr=%g\n",pl,pr[pl]);
      
    }


    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADBEAM2DKS || WHICHPROBLEM==RADBEAM2DKSVERT){


    RADBEAM2D_RHOAMB=1.e0/RHOBAR;
    RADBEAM2D_TAMB=1e7/TEMPBAR;
    RADBEAM2D_BLOB=0; // whether to put blob in way of beam
    RADBEAM2D_BLOBW=.1;
    RADBEAM2D_BLOBP=100000.;
    RADBEAM2D_BLOBX=10.;
    RADBEAM2D_BLOBZ=Pi/20.;
    RADBEAM2D_PAR_D=1./RHOBAR;
    RADBEAM2D_PAR_E=1e-4/RHOBAR;

    // BEAM PROPERTIES
    RADBEAM2D_IFBEAM=1; // whether to have a beam
    RADBEAM2D_TLEFT=1e9/TEMPBAR;
    //    RADBEAM2D_NLEFT=0.95; // >~0.95 and code fails with SPCMINKMETRIC for BEAMNO=1
    //    RADBEAM2D_NLEFT=0.99; // testing GODMARK KORALTODO
    //RADBEAM2D_NLEFT=0.999; // code
    //    RADBEAM2D_NLEFT=0.99999; // paper  // major problems with SPCMINKMETRIC
    RADBEAM2D_NLEFT=0.9999; // works with harmrad paper setup with PPM

    // avoid hitting gamma ceiling
    GAMMAMAXRAD=MAX(GAMMAMAXRAD,4.0*1.0/sqrt(1.0-RADBEAM2D_NLEFT*RADBEAM2D_NLEFT));


    if (RADBEAM2D_BEAMNO==1){
      RADBEAM2D_BEAML=2.9;
      RADBEAM2D_BEAMR=3.1;
    }
    else if (RADBEAM2D_BEAMNO==2){
      RADBEAM2D_BEAML=5.8;
      RADBEAM2D_BEAMR=6.2;
    }
    else if (RADBEAM2D_BEAMNO==3){
      RADBEAM2D_BEAML=15.5;
      RADBEAM2D_BEAMR=16.5;
    }
    else if (RADBEAM2D_BEAMNO==4){
      RADBEAM2D_BEAML=37;
      RADBEAM2D_BEAMR=43;
    }


    FTYPE Fx,Fy,Fz;

    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    *whichvel=VEL4;
    *whichcoord=MCOORD;

    // if BLOB, override density with blob density
    FTYPE rhoblob;
    rhoblob=RADBEAM2D_RHOAMB*(1.+RADBEAM2D_BLOBP*exp(-((xx-RADBEAM2D_BLOBX)*(xx-RADBEAM2D_BLOBX)+(yy)*(yy)+(zz-RADBEAM2D_BLOBZ)*(zz-RADBEAM2D_BLOBZ))/RADBEAM2D_BLOBW/RADBEAM2D_BLOBW));

    //zaczynam jednak od profilu analitycznego:   
    FTYPE ERADAMB;
    FTYPE rho,uint,Vr;
    if(RADBEAM2D_FLATBACKGROUND){
      Vr=0.;
      rho=RADBEAM2D_RHOAMB;
      if(RADBEAM2D_BLOB) rho=rhoblob;
      uint=calc_PEQ_ufromTrho(RADBEAM2D_TAMB,rho);
      ERADAMB=calc_LTE_EfromT(RADBEAM2D_TAMB);
      //   ERADAMB=calc_LTE_Efromurho(uint,rho);
    }
    else{
      FTYPE r=V[1];
      FTYPE mD=RADBEAM2D_PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
      FTYPE mE=RADBEAM2D_PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
      Vr=sqrt(2./r)*(1.-2./r);

      // get metric grid geometry for these ICs
      int getprim=0;
      struct of_geom geomrealdontuse;
      struct of_geom *ptrgeomreal=&geomrealdontuse;
      gset(getprim,*whichcoord,i,j,k,ptrgeomreal);

      FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeomreal->gcov[GIND(1,1)]);
      rho=RADBEAM2D_PAR_D/(r*r*sqrt(2./r));
      if(RADBEAM2D_BLOB) rho += rhoblob;
      FTYPE T=RADBEAM2D_TAMB;
      //   FTYPE ERAD=calc_LTE_EfromT(T);
      uint=mE/W;
      ERADAMB=calc_LTE_Efromurho(uint,rho);
    }

    Fx=Fy=Fz=0;

    //test!
    //Vr=0.7;

    pr[RHO] = rho ;
    pr[UU]  = uint;
    pr[U1]  = -Vr; // VEL4
    pr[U2]  = 0 ;  // static in VEL4
    pr[U3]  = 0 ;  // static in VEL4

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    if(PRAD0>=0){
      pr[PRAD0] = ERADAMB;
      pr[PRAD1] = 0.0 ; // static in VEL4
      pr[PRAD2] = 0.0 ;    
      pr[PRAD3] = 0.0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERADAMB;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;

      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL,pradffortho,pr, pr);
    }

    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==ATMSTATIC){

    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    *whichvel=VEL4;
    *whichcoord=MCOORD;


    FTYPE rho0=1.;
    FTYPE r0=2.e6;
    FTYPE u0=0.0001;
    FTYPE r=xx;

    FTYPE rho,uint;
    rho=rho0*r0/r;
    uint=u0*r0*r0/r/r;
    
    //    FTYPE E=exp(1.);
    
    uint=(4*r*u0 - 4*r*u0*gamideal - 2*r*rho0 + 2*r0*rho0 + r*r0*rho0*Log(-2 + r) - r*r0*rho0*Log(r) - r*r0*rho0*Log(-2 + r0) + r*r0*rho0*Log(r0))/(4*r - 4*r*gamideal);

    FTYPE uradx,urady,uradz;
    uradx=urady=uradz=0.;
    
    pr[RHO] = rho ;
    pr[UU]  = uint;
    pr[U1]  = 0 ;
    pr[U2]  = 0 ;    
    pr[U3]  = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    if(PRAD0>=0){
      // pr[PRAD0] = ERADLIMIT;
      pr[PRAD0] = uint*1E-20;
      pr[PRAD1] = uradx ;
      pr[PRAD2] = urady ;    
      pr[PRAD3] = uradz ;
    }
    // no transformations required since only setting fluid-frame E that is PRAD0 itself. (i.e. urad(xyz)=0 and ufluid=0)

    // *whichvel=WHICHVEL;
    *whichvel=VEL4;
    *whichcoord=MCOORD;

    return(0);
  }


  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADATM){


    // KORALTODO: t=0 R^t_t [lab] R^t_x [lab] or prad0 (E rad frame] prad1 [rad 4-vel] are off by about 16% compared to koral  UNITS CONSTANTS.

    RADATM_MDOTEDD=(2.23/16.*1e18*MPERSUN)/(MBAR/TBAR);
    RADATM_LUMEDD=(1.25e38*MPERSUN)/(ENBAR/TBAR);
    RADATM_THINRADATM=1;
    //    RADATM_FERATIO=.99999; // koral code
    RADATM_FERATIO=.99; // koral paper

    // avoid hitting gamma ceiling
    GAMMAMAXRAD=MAX(GAMMAMAXRAD,2.0*1.0/sqrt(1.0-RADATM_FERATIO*RADATM_FERATIO));


#define WHICHRADATM 0 // 0,1,2,3

    if(WHICHRADATM==0){
      RADATM_FRATIO=1E-10; // 1 = edd limit.  They ran 1E-10, 0.1, 0.5, 1.0.
    }
    if(WHICHRADATM==1){
      RADATM_FRATIO=.1; // 1 = edd limit.  They ran 1E-10, 0.1, 0.5, 1.0.
    }
    if(WHICHRADATM==2){
      RADATM_FRATIO=.5; // 1 = edd limit.  They ran 1E-10, 0.1, 0.5, 1.0.
    }
    if(WHICHRADATM==3){
      RADATM_FRATIO=1.0; // 1 = edd limit.  They ran 1E-10, 0.1, 0.5, 1.0.
    }


    RADATM_RHOAMB=1E-15/RHOBAR; // 1E-15 is in cgs
    RADATM_TAMB=1.e6/TEMPBAR;


    FTYPE MINX=Rin_array[1];
    FTYPE kappaesperrho=calc_kappaes_user(1,0, 0,0,0);
    FTYPE FLUXLEFT=RADATM_FRATIO/kappaesperrho/pow(MINX,2.0);



    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    *whichvel=VEL4;
    *whichcoord=MCOORD;

    //at outern boundary
    FTYPE f;
    if(WHICHRADSOURCEMETHOD!=SOURCEMETHODNONE){
      f = (FTYPE)kappaesperrho*FLUXLEFT*MINX*MINX;
    }
    else f=0;


    FTYPE p0=RADATM_RHOAMB*RADATM_TAMB;
    FTYPE KKK=p0/pow(RADATM_RHOAMB,gamideal);
    FTYPE C3=gamideal*KKK/(gamideal-1.)*pow(RADATM_RHOAMB,gamideal-1.)-(1.-f)*(1./MINX+0.*1./MINX/MINX+0.*4./3./MINX/MINX/MINX);

    FTYPE rho=pow((gamideal-1.0)/gamideal/KKK*(C3+(1.-f)*(1./xx + 0.*1./xx/xx + 0.*4./3./xx/xx/xx)),1./(gamideal-1.0));

    FTYPE pre=KKK*pow(rho,gamideal);

    FTYPE uint=pre/(gamideal-1.0);

    FTYPE Fz=0;
    FTYPE Fy=0.;
    FTYPE Fx=FLUXLEFT*(MINX/xx)*(MINX/xx);

    FTYPE ERAD;
    if(RADATM_THINRADATM){
      ERAD=Fx/RADATM_FERATIO;
    }
    else{
      ERAD=calc_LTE_EfromT(calc_PEQ_Tfromurho(uint,rho));
    }

    if(0){//DEBUG:
      dualfprintf(fail_file,"i=%d f=%g p0=%g KKK=%g C3=%g rho=%g uint=%g Fx=%g ERAD=%g : kappaesperrho=%g \n",i,f,p0,KKK,C3,rho,uint,Fx,ERAD , kappaesperrho);
      
      dualfprintf(fail_file,"i=%d f=%g p0=%g KKK=%g C3=%g rho=%g uint=%g Fx=%g ERAD=%g : kappaesperrho=%g \n",i,f,p0*UBAR,KKK*UBAR/pow(RHOBAR,gamideal),C3,rho*RHOBAR,uint*UBAR,Fx*ENBAR/TBAR/LBAR/LBAR,ERAD*UBAR , kappaesperrho*OPACITYBAR);
    }

    

    pr[RHO] = rho ;
    pr[UU]  = uint;
    pr[U1]  = 0 ;
    pr[U2]  = 0 ;    
    pr[U3]  = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    *whichvel=VEL4;
    *whichcoord=MCOORD;

    if(PRAD0>=0){
      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERAD;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;

      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL,pradffortho, pr, pr);

      if(0){ // DEBUG

        // get metric grid geometry for these ICs
        int getprim=0;
        struct of_geom geomrealdontuse;
        struct of_geom *ptrgeomreal=&geomrealdontuse;
        gset(getprim,*whichcoord,i,j,k,ptrgeomreal);

        dualfprintf(fail_file,"AFTER: i=%d rho=%g uint=%g vx=%g ERAD=%g uradx=%g\n",i,pr[RHO]*RHOBAR,pr[UU]*UBAR,pr[U1]*sqrt(ptrgeomreal->gcov[GIND(1,1)])*VBAR,pr[URAD0]*UBAR,pr[URAD1]*sqrt(ptrgeomreal->gcov[GIND(1,1)])*VBAR);
      }

      // compared to koral, this is how koral would get CGS:
      // fprintf(stderr,"i=%d f=%g p0=%g KKK=%g C3=%g rho=%g uint=%g Fx=%g ERAD=%g : kappaesperrho=%g \n",ix,f,endenGU2CGS(p0),endenGU2CGS(KKK)/powl(rhoGU2CGS(1.0),GAMMA),C3,rhoGU2CGS(pp[0]),endenGU2CGS(pp[1]),fluxGU2CGS(Fx), endenGU2CGS(E) , kappaGU2CGS(KAPPAES));
      // fprintf(stderr,"AFTER: i=%d rho=%g uint=%g vx=%g ERAD=%g uradx=%g\n",ix,rhoGU2CGS(pp[0]),endenGU2CGS(pp[1]),velGU2CGS(pp[2]),endenGU2CGS(pp[6]),velGU2CGS(pp[7]));


    }

    return(0);
  }



  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADWALL){

    
    // direct assignments since simple
    pr[RHO] = 1.0 ;
    pr[UU] = 1.0;
    pr[U1] = 0 ;
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    if(PRAD0>=0){
      // direct assignments since simple
      pr[PRAD0] = 1.0;
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;
    }

    // no transformations required since only setting fluid-frame E that is PRAD0 itself since ufluid=F=0

    *whichvel=WHICHVEL;
    *whichcoord=CARTMINKMETRIC2;
    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADWAVE){

    FTYPE rho,ERAD,uint,Fx,Fy,Fz;
    FTYPE vx;
    FTYPE xx,yy,zz;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];


    // default
    Fx=Fy=Fz=0.0;


    FTYPE time=0.;

    //Jiang+12 waves
    if(RADWAVE_NWAVE==5){

      //printf("RHOZERO = %g\nUINT = %g\nT = %g\nERAD = %g\nARAD = %g\n",RADWAVE_RHOZERO,RADWAVE_UINT,RADWAVE_TEMP,RADWAVE_ERAD,ARAD_RAD_CODE);getchar();


      rho=RADWAVE_RHOZERO*(1+RADWAVE_DRRE*exp(-RADWAVE_OMIM*time)*(cos(RADWAVE_OMRE*time-RADWAVE_KK*xx)-RADWAVE_DRIM/RADWAVE_DRRE*sin(RADWAVE_OMRE*time-RADWAVE_KK*xx)));
      //FTYPE RADWAVE_DURE=RADWAVE_DPRE/(gam-1.); FTYPE RADWAVE_DUIM=RADWAVE_DPIM/(gam-1.);
      uint=RADWAVE_UINT*(1.+RADWAVE_DURE*exp(-RADWAVE_OMIM*time)*(cos(RADWAVE_OMRE*time-RADWAVE_KK*xx)-RADWAVE_DUIM/RADWAVE_DURE*sin(RADWAVE_OMRE*time-RADWAVE_KK*xx))) ;
      FTYPE cs=1/RADWAVE_CC;
      vx=0. + RADWAVE_DVRE*exp(-RADWAVE_OMIM*time)*(cos(RADWAVE_OMRE*time-RADWAVE_KK*xx)-RADWAVE_DVIM/RADWAVE_DVRE*sin(RADWAVE_OMRE*time-RADWAVE_KK*xx)) ; //RADWAVE_DVRE absolute!
      ERAD=RADWAVE_ERAD*(1+RADWAVE_DERE*exp(-RADWAVE_OMIM*time)*(cos(RADWAVE_OMRE*time-RADWAVE_KK*xx)-RADWAVE_DEIM/RADWAVE_DERE*sin(RADWAVE_OMRE*time-RADWAVE_KK*xx)));
      Fx=0. + RADWAVE_ERAD*RADWAVE_DFRE*exp(-RADWAVE_OMIM*time)*(cos(RADWAVE_OMRE*time-RADWAVE_KK*xx)-RADWAVE_DFIM/RADWAVE_DFRE*sin(RADWAVE_OMRE*time-RADWAVE_KK*xx));
      Fz=Fy=0.;

      //rho=RADWAVE_RHOZERO;
      //uint=RADWAVE_UINT;
      //ERAD=RADWAVE_ERAD;
      //vx=0;
      //Fx=0.;
    }

    //hydro density wave
    if(RADWAVE_NWAVE==1){
      rho=RADWAVE_RHOZERO*(1.+RADWAVE_AAA*cos(RADWAVE_KK*xx));
      uint=RADWAVE_UINT;
      vx=RADWAVE_VX;
      ERAD=1E-10*uint; // no radiation
      Fx=Fz=Fy=0.;
    }

    //hydro sound wave
    if(RADWAVE_NWAVE==2){
      rho=RADWAVE_RHOZERO*(1.+RADWAVE_AAA*cos(RADWAVE_KK*xx));
      uint=RADWAVE_UINT*(1.+gam*RADWAVE_AAA*cos(RADWAVE_KK*xx));
      FTYPE cs=1./RADWAVE_CC;
      vx=RADWAVE_AAA*cos(RADWAVE_KK*xx)*cs;
      ERAD=RADWAVE_ERAD; // KORALTODO: Why does koral not set #define RADIATION for this test?  avoid test.
      Fx=Fz=Fy=0.;
    }

    //radiative hydro density wave
    if(RADWAVE_NWAVE==3){
      rho=RADWAVE_RHOZERO*(1.+RADWAVE_AAA*cos(RADWAVE_KK*xx));
      uint=RADWAVE_UINT;
      vx=RADWAVE_VX;
      ERAD=RADWAVE_ERAD;
      Fx=Fz=Fy=0.;
    }

    //radiative sound wave
    if(RADWAVE_NWAVE==4){
      rho=RADWAVE_RHOZERO*(1.+RADWAVE_GASFACTOR*RADWAVE_AAA*cos(RADWAVE_KK*xx));
      uint=RADWAVE_UINT*(1.+RADWAVE_GASFACTOR*gam*RADWAVE_AAA*cos(RADWAVE_KK*xx));
      FTYPE cs=1./RADWAVE_CC;
      vx=RADWAVE_GASFACTOR*RADWAVE_AAA*cos(RADWAVE_KK*xx)*cs;
      ERAD=RADWAVE_ERAD*(1.+RADWAVE_ERADFACTOR*RADWAVE_AAA*cos(RADWAVE_KK*xx));
      Fx=Fz=Fy=0.;
    }



    pr[RHO] = rho;
    pr[UU] = uint;
    pr[U1] = vx ; // vx is 3-velocity
    pr[U2] = 0 ;    
    pr[U3] = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
  
    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    *whichvel=VEL3;
    *whichcoord=CARTMINKMETRIC2;

    if(PRAD0>=0){
      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERAD;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;


      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL, pradffortho, pr, pr);
      //  PLOOPRADONLY(pl) dualfprintf(fail_file,"FOO1: i=%d pl=%d pr=%g\n",ptrgeomreal->i,pl,pr[pl]);
    }

 
    return(0);
  

  }


  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==KOMIPROBLEM){

    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    FTYPE x = V[1];

    if(WHICHKOMI>=1 && WHICHKOMI<=9){
      *whichvel=VEL4;
      *whichcoord=CARTMINKMETRIC2;

      FTYPE pleft[NPR], pright[NPR], P;
      //zero out initial conditions
      PALLLOOP(pl) pleft[pl] = 0.;
      PALLLOOP(pl) pright[pl] = 0.;
    
    
      //fast shock
      if(WHICHKOMI==1){
        //left state
        pleft[U1] = 25.0;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 20.0;
        pleft[B2] = 25.02;
        pleft[B3] =  0.0;
        P = 1.0;
        pleft[RHO] = 1.;
        pleft[UU] = P/(gam-1);

        //right state
        pright[U1] = 1.091;
        pright[U2] =  0.3923;
        pright[U3] =  0.0;
        pright[B1] = 20.0;
        pright[B2] = 49.0;
        pright[B3] =  0.0;
        P = 367.5;
        pright[RHO] = 25.48;
        pright[UU] = P/(gam-1);
      }
      //slow shock
      else if(WHICHKOMI==2){
        //left state
        pleft[U1] = 1.53;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 10.0;
        pleft[B2] = 18.28;
        pleft[B3] = 0.0;
        P = 10.0;
        pleft[RHO] = 1.;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] =  .9571;
        pright[U2] = -0.6822;
        pright[U3] =  0.0;
        pright[B1] = 10.0;
        pright[B2] = 14.49;
        pright[B3] =  0.0;
        P = 55.36;
        pright[RHO] = 3.323;
        pright[UU] = P/(gam-1);
      }
      //fast switch-off rarefaction
      else if(WHICHKOMI==3){
        //left state
        pleft[U1] = -2.0;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 2.0;
        pleft[B2] = 0.0;
        pleft[B3] = 0.0;
        P = 1.0;
        pleft[RHO] = 0.1;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] = -0.212;
        pright[U2] = -0.590;
        pright[U3] =  0.0;
        pright[B1] =  2.0;
        pright[B2] =  4.710;
        pright[B3] =  0.0;
        P = 10.0;
        pright[RHO] = 0.562;
        pright[UU] = P/(gam-1);
      }
      //slow switch-on rarefaction
      else if(WHICHKOMI==4){
        //left state
        pleft[U1] = -0.765;
        pleft[U2] = -1.386;
        pleft[U3] =  0.0;
        pleft[B1] = 1.0;
        pleft[B2] = 1.022;
        pleft[B3] = 0.0;
        P = 0.1;
        pleft[RHO] = 1.78e-3;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] =  0.0;
        pright[U2] =  0.0;
        pright[U3] =  0.0;
        pright[B1] =  1.0;
        pright[B2] =  0.0;
        pright[B3] =  0.0;
        P = 1.0;
        pright[RHO] = 1.0;
        pright[UU] = P/(gam-1);
      }
      //alfven wave and compound wave
      else if(WHICHKOMI==5 || WHICHKOMI==6){
        //left state
        pleft[U1] = 0;
        pleft[U2] = 0;
        pleft[U3] = 0.0;
        pleft[B1] = 3.0;
        pleft[B2] = 3.0;
        pleft[B3] = 0.0;
        P = 1.0;
        pleft[RHO] = 1.0;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] =  3.70;
        pright[U2] =  5.76;
        pright[U3] =  0.0;
        pright[B1] =  3.0;
        pright[B2] =  -6.857;
        pright[B3] =  0.0;
        P = 1.0;
        pright[RHO] = 1.0;
        pright[UU] = P/(gam-1);
      }
      //Shock tube 1
      else if(WHICHKOMI==7){
        //left state
        pleft[U1] =  0.0;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 1.0;
        pleft[B2] = 0.0;
        pleft[B3] = 0.0;
        P = 1000.0;
        pleft[RHO] = 1.0;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] =  0.0;
        pright[U2] =  0.0;
        pright[U3] =  0.0;
        pright[B1] =  1.0;
        pright[B2] =  0.0;
        pright[B3] =  0.0;
        P = 1.0;
        pright[RHO] = 0.1;
        pright[UU] = P/(gam-1);
      }
      //Shock tube 2
      else if(WHICHKOMI==8){
        //left state
        pleft[U1] =  0.0;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 0.0;
        pleft[B2] = 20.0;
        pleft[B3] = 0.0;
        P = 30.;
        pleft[RHO] = 1.0;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] =  0.0;
        pright[U2] =  0.0;
        pright[U3] =  0.0;
        pright[B1] =  0.0;
        pright[B2] =  0.0;
        pright[B3] =  0.0;
        P = 1.0;
        pright[RHO] = 0.1;
        pright[UU] = P/(gam-1);
      }
      //Collision
      else if(WHICHKOMI==9){
        //left state
        pleft[U1] =  5.0;
        pleft[U2] =  0.0;
        pleft[U3] =  0.0;
        pleft[B1] = 10.0;
        pleft[B2] = 10.0;
        pleft[B3] = 0.0;
        P = 1.0;
        pleft[RHO] = 1.0;
        pleft[UU] = P/(gam-1);
      
        //right state
        pright[U1] = -5.0;
        pright[U2] =  0.0;
        pright[U3] =  0.0;
        pright[B1] =  10.0;
        pright[B2] = -10.0;
        pright[B3] =  0.0;
        P = 1.0;
        pright[RHO] = 1.0;
        pright[UU] = P/(gam-1);
      }

      if(WHICHKOMI==5 || WHICHKOMI==6){ // Alfven wave or compound wave
        FTYPE xcl,xcr,xint;
        //        xcl=xcr=0.0;
        if(WHICHKOMI==5){ xcl=-0.5; xcr=0.0; }
        if(WHICHKOMI==6){ xcl=-0.025; xcr=0.0; } // 4 cells according to figure 5 in Komi 1999, but seems to really be 2 cells unless expanded range of box relative to WHICHKOMI==5 for odd reason
        if(x<=xcl) PALLLOOP(pl) pr[pl] = pleft[pl];
        else if(x>xcr) PALLLOOP(pl) pr[pl] = pright[pl];
        else{
          //          FTYPE xint=(x- xcl)/(xcr - xcl);
          // default left-right state values
          //          PALLLOOP(pl) pr[pl] = pleft[pl] + (pright[pl]-pleft[pl])*xint;

          // but rotate field by \pi
          if(WHICHKOMI==5){ xcl=-0.5; xcr=0.0; }
          if(WHICHKOMI==6){ xcl=-0.025; xcr=0.0; } // 4 cells according to figure 5 in Komi 1999, but seems to really be 2 cells unless expanded range of box relative to WHICHKOMI==5 for odd reason
          FTYPE phi0;
          if(x<=xcl) phi0=0.0;
          else if((x>xcl)&&(x<xcr)) phi0=0.5*M_PI*(x-xcl)/(xcr-xcl);
          else if(x>=xcr) phi0=0.5*M_PI;

          FTYPE phi1;
          if(x<=xcl) phi1=0.0;
          else if((x>xcl)&&(x<xcr)) phi1=M_PI*(x-xcl)/(xcr-xcl);
          else if(x>=xcr) phi1=M_PI;

          //          xint=sin(phi0);
          // rotate all by pi
          PALLLOOP(pl) pr[pl] = pleft[pl]  + (pright[pl]-pleft[pl])*sin(phi0);
          //          PALLLOOP(pl) if(pl==B2) pr[pl] = pleft[B2]*sin(phi1) + pright[B2]*cos(phi1);
          // set B3 so bsq is constant
          
          FTYPE usq = pr[U1]*pr[U1]+pr[U2]*pr[U2]+pr[U3]*pr[U3];
          FTYPE gamma = sqrt(1.0 + fabs(usq));

          FTYPE usqleft = pleft[U1]*pleft[U1]+pleft[U2]*pleft[U2]+pleft[U3]*pleft[U3];
          FTYPE gammaleft = sqrt(1.0 + fabs(usqleft));

          // see komi_fake_alfven.nb
          FTYPE mybsqconst=Power(pleft[B3],2)/Power(gammaleft,2) + 
            (-1.*pleft[B1]*pleft[U1] - 1.*pleft[B2]*pleft[U2])*
            (pleft[B1]*pleft[U1] + pleft[B2]*pleft[U2]) + 
            Power(pleft[B1]/gammaleft + (pleft[U1]*
                                          (pleft[B1]*pleft[U1] + pleft[B2]*pleft[U2]))/gammaleft,2) + 
            Power(pleft[B2]/gammaleft + (pleft[U2]*
                                          (pleft[B1]*pleft[U1] + pleft[B2]*pleft[U2]))/gammaleft,2);

          dualfprintf(fail_file,"x=%g mybsqconst=%g\n",x,mybsqconst);

          FTYPE disc1=Power(gamma,2)*mybsqconst - 1.*Power(pr[B1],2) - 1.*Power(pr[B2],2) - 
            2.*Power(pr[B1],2)*Power(pr[U1],2) + 
            Power(gamma,2)*Power(pr[B1],2)*Power(pr[U1],2) - 
            1.*Power(pr[B1],2)*Power(pr[U1],4) - 4.*pr[B1]*pr[B2]*pr[U1]*pr[U2] + 
            2.*Power(gamma,2)*pr[B1]*pr[B2]*pr[U1]*pr[U2] - 
            2.*pr[B1]*pr[B2]*Power(pr[U1],3)*pr[U2] - 2.*Power(pr[B2],2)*Power(pr[U2],2) + 
            Power(gamma,2)*Power(pr[B2],2)*Power(pr[U2],2) - 
            1.*Power(pr[B1],2)*Power(pr[U1],2)*Power(pr[U2],2) - 
            1.*Power(pr[B2],2)*Power(pr[U1],2)*Power(pr[U2],2) - 
            2.*pr[B1]*pr[B2]*pr[U1]*Power(pr[U2],3) - 1.*Power(pr[B2],2)*Power(pr[U2],4);
          
          FTYPE disc2=Power(gamma,2)*mybsqconst - 1.*Power(pr[B1],2) - 1.*Power(pr[B2],2) - 
            2.*Power(pr[B1],2)*Power(pr[U1],2) + 
            Power(gamma,2)*Power(pr[B1],2)*Power(pr[U1],2) - 
            1.*Power(pr[B1],2)*Power(pr[U1],4) - 4.*pr[B1]*pr[B2]*pr[U1]*pr[U2] + 
            2.*Power(gamma,2)*pr[B1]*pr[B2]*pr[U1]*pr[U2] - 
            2.*pr[B1]*pr[B2]*Power(pr[U1],3)*pr[U2] - 2.*Power(pr[B2],2)*Power(pr[U2],2) + 
            Power(gamma,2)*Power(pr[B2],2)*Power(pr[U2],2) - 
            1.*Power(pr[B1],2)*Power(pr[U1],2)*Power(pr[U2],2) - 
            1.*Power(pr[B2],2)*Power(pr[U1],2)*Power(pr[U2],2) - 
            2.*pr[B1]*pr[B2]*pr[U1]*Power(pr[U2],3) - 1.*Power(pr[B2],2)*Power(pr[U2],4);
          
          if(disc1>=0.0){ PALLLOOP(pl) if(pl==B3) pr[pl] = +sqrt(disc1); }
          else if(disc2>=0){ PALLLOOP(pl) if(pl==B3) pr[pl] = -sqrt(disc2); }
          else PALLLOOP(pl){ if(pl==B3) pr[pl] = 0.0; }

          //          PALLLOOP(pl) if(pl==B2){ xint=cos(phi0); pr[pl] = pleft[pl] + (pright[pl]-pleft[pl])*xint; }
          //          PALLLOOP(pl) if(pl==B3){ xint=sin(phi0); pr[pl] = pleft[B2] + (pright[B2]-pleft[B2])*xint; }

        }
      }
      else{
        if (x<=0) {
          PALLLOOP(pl) pr[pl] = pleft[pl];
        }
        else if (x>0) {
          PALLLOOP(pl) pr[pl] = pright[pl];
        }
      }
    }

    if(WHICHKOMI>=101 && WHICHKOMI<=109){
      FTYPE E[NDIM],B[NDIM];
      FTYPE x0,dx0;
      FTYPE bcon[NDIM],vcon[NDIM],econ[NDIM];
      FTYPE phi0;
      FTYPE KK;
      FTYPE B0;
      FTYPE muf;

      // defaults
      if(EOMTYPE==EOMFFDE||EOMTYPE==EOMFFDE2){
        pr[RHO]=pr[UU]=0;
      }
      else{
        pr[RHO]=0.1; // all fields below are order unity overall, but can pass through zero.
        pr[UU]=pr[RHO]; // so relativistically hot
      }

      pr[U1]=pr[U2]=pr[U3]=0.0;
      pr[B2]=pr[B3]=0;
      pr[B1]=0;

      
      extern void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr);
      extern void computeKK(FTYPE *pr, struct of_geom *geom, FTYPE *KK);
      extern void EBvetatopr(FTYPE *Econ, FTYPE *Bcon, FTYPE *veta, struct of_geom *geom, FTYPE *pr);
      extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
      extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);


      int TESTNUMBER=WHICHKOMI-101; // to translate to init.komtests.c numbering

      *whichvel=WHICHVEL;
      //*whichvel=VEL4;
      *whichcoord=MCOORD;
      // get metric grid geometry for these ICs
      int getprim=0;
      struct of_geom geomdontuse;
      struct of_geom *ptrgeom=&geomdontuse;
      gset(getprim,*whichcoord,i,j,k,ptrgeom);

      //      *whichvel=WHICHVEL;
      //      *whichcoord=PRIMECOORDS;
      //      struct of_geom geomdontuse;
      //      struct of_geom *ptrgeom=&geomdontuse;
      //      int loc=CENT;
      //      get_geometry(i,j,k,loc,ptrgeom);


      if(TESTNUMBER==0){ // Fast wave
        tf = 1;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        //  tf = 1;
        //  DTd=1E-5;
  

        E[1]=0;
        E[2]=0;
        B[3]=0.0;
        B[1]=0.0;
        x0=0.0;
        dx0=0.1;
        if(x-x0<=-dx0) B[2]=1.0;
        else if((x-x0>-dx0)&&(x-x0<dx0)) B[2]=1.0-(0.3/(2.0*dx0))*((x-x0)+dx0);
        else if(x-x0>=dx0) B[2]=0.7;

        muf=1.0;

        E[3]=1.0-muf*B[2];

        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);
        
        //  pr[U1]=0.9;
        
        //dualfprintf(fail_file,"pr[U1]=%21.15g pr[U2]=%21.15g\n",pr[U1],pr[U2]);
        
        computeKK(pr,ptrgeom,&KK);
        
        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);
        
      }
      if(TESTNUMBER==1){ // comoving Fast wave (NOT a Komissarov test)
        //tf = 1;
        //  DTd=tf/10.0;

        tf = 1;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E-5;
        
        
        bcon[3]=0.0;
        bcon[1]=0.0;
        x0=0.0;
        dx0=0.1;
        if(x-x0<=-dx0) bcon[2]=1.0;
        else if((x-x0>-dx0)&&(x-x0<dx0)) bcon[2]=1.0-(0.3/(2.0*dx0))*((x-x0)+dx0);
        else if(x-x0>=dx0) bcon[2]=0.7;
        
        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention
        
        
        vcon[1]=0.9999;
        vcon[2]=vcon[3]=0;
        
        //        int jj;
        //        SLOOPA(jj) vcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) bcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        vbtopr(vcon,bcon,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

      }
      if(TESTNUMBER==2){ // (nondegenerate) Alfven wave (not going to work with HARM)
        tf = 2;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        bcon[1]=bcon[2]=1.0;
  
        x0=0.0;
        if(x-x0<=-0.1) bcon[3]=1.0;
        else if((x-x0>-0.1)&&(x-x0<0.1)) bcon[3]=1.0+3.0/2.0*((x-x0)+0.1);
        else if(x-x0>=0.1) bcon[3]=1.3;

        econ[2]=econ[3]=0.0;

        // can be + or -
        //FTYPE CONSTECON1=1.3;
        //  econ[1]=-sqrt(-CONSTECON1+bcon[1]*bcon[1]+bcon[2]*bcon[3]+bcon[3]*bcon[3]);
        //  econ[1]=1-0.5*bcon[3];
        FTYPE CONSTECON1=1.0;
        econ[1]=sqrt(-CONSTECON1 + bcon[3]*bcon[3]);
        //  econ[1]=0.0;
        //  econ[1]=-bcon[3];

        vcon[1]=-0.5;
        vcon[2]=vcon[3]=0;

        //        int jj;
        //        SLOOPA(jj) econ[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) bcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) vcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBvetatopr(econ, bcon, vcon, ptrgeom, pr);
        //  vbtopr(vcon,bcon,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

      }

      if(TESTNUMBER==3){ // Degenerate Alfven wave
        tf = 2;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;
        bcon[1]=0.0;

  
        x0=0.0;
        if(x-x0<=-0.1) phi0=0.0;
        else if((x-x0>-0.1)&&(x-x0<0.1)) phi0=5.0/2.0*M_PI*((x-x0)+0.1);
        else if(x-x0>=0.1) phi0=M_PI*0.5;

        bcon[2]=2.0*cos(phi0);
        bcon[3]=2.0*sin(phi0);


        vcon[1]=0.5;
        vcon[2]=vcon[3]=0;

        //        int jj;
        //        SLOOPA(jj) vcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) bcon[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        vbtopr(vcon,bcon,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }
      if(TESTNUMBER==4){ // Three-wave problem
        tf = .75;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        x0=0.5;
        if(x<x0){
          B[1]=1.0;
          B[2]=1.5;
          B[3]=3.5;
          E[1]=-1.0;
          E[2]=-0.5;
          E[3]=0.5;
        }
        else{
          B[1]=1.0;
          B[2]=2.0;
          B[3]=2.3;
          E[1]=-1.5;
          E[2]=1.3;
          E[3]=-0.5;
        }

        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }

      if(TESTNUMBER==5){ // B^2-E^2<0 problem
        tf = .02;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        x0=0.5;
        if(x<x0){
          B[0]=0.0;
          B[1]=1.0;
          B[2]=1.0;
          B[3]=1.0;
          E[0]=0.0;
          E[1]=0.0;
          E[2]=0.5;
          E[3]=-0.5;
        }
        else{
          B[0]=0.0;
          B[1]=1.0;
          B[2]=-1.0;
          B[3]=-1.0;
          E[0]=0.0;
          E[1]=0.0;
          E[2]=0.5;
          E[3]=-0.5;
        }

        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }
      if(TESTNUMBER==6){ // smoothed B^2-E^2<0 problem
        tf = .02;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        x0=0.5;
        if(x-x0<-0.1){
          B[1]=1.0;
          B[2]=1.0;
          B[3]=1.0;
          E[1]=0.0;
          E[2]=0.5;
          E[3]=-0.5;
        }
        else if(x-x0>0.1){
          B[1]=1.0;
          B[2]=-1.0;
          B[3]=-1.0;
          E[1]=0.0;
          E[2]=0.5;
          E[3]=-0.5;
        }
        else if((x-x0>=-0.1)&&(x-x0<=0.1)){
          B[1]=1.0;
          B[2]=1.0+(x-x0+0.1)*(-2.0/0.2);
          B[3]=1.0+(x-x0+0.1)*(-2.0/0.2);
          E[1]=0.0;
          E[2]=0.5;
          E[3]=-0.5;
        }

        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }

      if(TESTNUMBER==7){ // Komissarov 2004 C3.1 Alfven wave
        // PARA generates crap on left side, but wave doesn't move
        // MC does very well
        // no obvious difference between HLL and LAXF
        // Athena1/2 ok
        tf = 2.0;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        //  B[1]=B[2]=E[3]=E[2]=0;
        B[1]=B[2]=E[3]=1.0;
        E[2]=0;


        if(x<=0.5){
          B[3]=1.0;
        }
        else if(x>=0.2+0.5){
          B[3]=1.3;
        }
        else{
          B[3]=1.0+0.15*(1.0+sin(5.0*M_PI*(x-0.1-0.5)));
        }
        E[1]=-B[3];

        //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }
      if(TESTNUMBER==8){ // Komissarov 2004 C3.2 Current Sheet
        tf = 1.0;
        int idt;
        for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=tf/10.0;

        E[1]=E[2]=E[3]=0.0;
        B[3]=0.0;
        B[1]=1.0;

        //B0=0.5; // fine
        B0 = 2.0;

        if(x<=0.5){
          B[2]=B0;
        }
        else{
          B[2]=-B0;
        }


        //        int jj;
        //        SLOOPA(jj) E[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        //        SLOOPA(jj) B[jj]*=sqrt(ptrgeom->gcon[GIND(jj,jj)]);
        EBtopr(E,B,ptrgeom,pr);
        //EBtopr_2(E,B,ptrgeom,pr);

        computeKK(pr,ptrgeom,&KK);

        dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


      }

    }// end over init.komtests.c from FFDE code tests


    
  
    if(FLUXB==FLUXCTSTAG){
      //can ignore half a cell shift for B1: it does not change across the interface so does not matter
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }
    
    
    pr[PRAD0] = 10.0*ERADLIMIT ; // so doesn't hit floor and confuse debug info
    pr[PRAD1] = 0 ;
    pr[PRAD2] = 0 ;
    pr[PRAD3] = 0 ;
    
    
    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBONDI){


    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    *whichcoord=MCOORD;
    // get metric grid geometry for these ICs
    int getprim=0;
    struct of_geom geomrealdontuse;
    struct of_geom *ptrgeomreal=&geomrealdontuse;
    gset(getprim,*whichcoord,i,j,k,ptrgeomreal);

    
    FTYPE rho,ERAD,uint;
    FTYPE rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;

    FTYPE Fx,Fy,Fz;
    Fx=Fy=Fz=0;

    //at outern boundary
    r=RADBONDI_MAXX;
    ur=-sqrt(2./r);
    rho0=-RADBONDI_MDOTPEREDD*RADBONDI_MDOTEDD/(4.*Pi*r*r*ur);
    Tgas0=RADBONDI_TGAS0;
            
    //at given cell
    r=xx;
    ur=-sqrt(2./r);    
    ut=sqrt((-1.-ur*ur*ptrgeomreal->gcov[GIND(1,1)])/ptrgeomreal->gcov[GIND(0,0)]);
    vx=ur/ut;  
    rho=-RADBONDI_MDOTPEREDD*RADBONDI_MDOTEDD/(4.*Pi*r*r*ur);
    Tgas=Tgas0*pow(rho/rho0,gam-1.);      

    uint=calc_PEQ_ufromTrho(Tgas,rho);

    pgas=rho*Tgas;
    prad=RADBONDI_PRADGAS*pgas;
    ERAD=prad*3.;
    
    pr[RHO] = rho ;
    pr[UU]  = uint;
    // KORAL sets VEL3 type quantity vx, but then feeds that directly into koral's prad_ff2lab() that takes in relative 4-velocity.  So set ur instead as 4-velocity and assume not extremely close to hole.
    //    pr[U1]  = vx;
    pr[U1]  = ur;
    pr[U2]  = 0 ;    
    pr[U3]  = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    
    //    *whichvel=VEL3;
    *whichvel=VEL4;

    if(PRAD0>=0){
      pr[PRAD0] = ERAD;
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;


      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = ERAD;
      pradffortho[PRAD1] = Fx;
      pradffortho[PRAD2] = Fy;
      pradffortho[PRAD3] = Fz;


      if(0){//DEBUG:
        dualfprintf(fail_file,"i=%d rho=%g uint=%g ERAD=%g\n",i,rho,uint,ERAD);
      
        dualfprintf(fail_file,"i=%d rho=%g uint=%g ERAD=%g\n",i,rho*RHOBAR,uint*UBAR,ERAD*UBAR);
      }


      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,ptrgeomreal, pradffortho, pr, pr);

      if(0){ // DEBUG
        dualfprintf(fail_file,"AFTER: i=%d rho=%g uint=%g vx=%g ERAD=%g uradx=%g\n",i,pr[RHO]*RHOBAR,pr[UU]*UBAR,pr[U1]*sqrt(ptrgeomreal->gcov[GIND(1,1)])*VBAR,pr[URAD0]*UBAR,pr[URAD1]*sqrt(ptrgeomreal->gcov[GIND(1,1)])*VBAR);
      }
    }

    return(0);
  }


  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADDOT){

    FTYPE xx,yy,zz,rsq;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    xx=V[1];
    yy=V[2];
    zz=V[3];

    pr[RHO] = 1.0;
    pr[UU]  = 1.0;
    pr[U1]  = 0 ;
    pr[U2]  = 0 ;    
    pr[U3]  = 0 ;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    *whichvel=VEL4;
    *whichcoord=MCOORD;

    if(PRAD0>=0){
      pr[PRAD0] = 0 ; // so triggers failure if used
      pr[PRAD1] = 0 ;
      pr[PRAD2] = 0 ;    
      pr[PRAD3] = 0 ;

      //E, F^i in orthonormal fluid frame
      FTYPE pradffortho[NPR];
      pradffortho[PRAD0] = RADDOT_LTEFACTOR*calc_LTE_Efromurho(pr[RHO],pr[UU]);
      pradffortho[PRAD1] = 0;
      pradffortho[PRAD2] = 0;
      pradffortho[PRAD3] = 0;

      if(startpos[1]+i==RADDOT_IDOT && startpos[2]+j==RADDOT_JDOT && startpos[3]+k==RADDOT_KDOT){
        //      dualfprintf(fail_file,"GOT INITIAL DOT\n");
        if(N1==1) pradffortho[PRAD0] *= RADDOT_F1;
        else{
          pradffortho[PRAD0]*=RADDOT_F2;
          pradffortho[PRAD2]=RADDOT_FYDOT*pradffortho[PRAD0];
        }
      }// end if DOT

      
      // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
      prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,NULL, pradffortho, pr, pr);
    }

    return(0);
  }

  if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK || WHICHPROBLEM==RADDONUT){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not just BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3; // use 3-velocity since later set donut using VEL3
    // get metric grid geometry for these ICs
    int getprim=0;
    struct of_geom geomraddontuse;
    struct of_geom *ptrgeomrad=&geomraddontuse;
    gset(getprim,*whichcoord,i,j,k,ptrgeomrad);

    
    // KORAL:
    // atmtype=0 : -1.5 -2.5
    // atmtype=1 : -2.0 -2.5

    if(1){ // SUPERMADNEW
      int set_fieldtype(void);
      int FIELDTYPE=set_fieldtype();
      if(FIELDTYPE==FIELDJONMAD){
        // so bsq/rho isn't super high at large radii with current choice for magnetic field
        pr[RHO]=RADNT_RHOATMMIN*pow(r/RADNT_ROUT,-1.5); // SUPERMADNEW
      }
    }
    else{
      pr[RHO]=RADNT_RHOATMMIN*pow(r/RADNT_ROUT,-1.5);
    }

    pr[UU]=RADNT_UINTATMMIN*pow(r/RADNT_ROUT,-2.5);

    set_zamo_velocity(*whichvel,ptrgeomrad,pr); // only sets U1-U3 to zamo

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }


    // KORAL:
    // atmtype=0 : pr[URAD0]=ERADATMMIN and zamo vels
    // atmtype=1 : pr[URAD0]=ERADATMMIN and ncon={0,-1,0,0} with set_ncon_velocity(whichvel,1000.0,ncon,ptrgeomrad,uconreal);
    // atmtype=2 : pr[URAD0]=ERADATMMIN*pow(rout/r,4) pr[URAD1-URAD3]  with ncon={0,-gammamax*(pow(r/rout,1.)),0,0} again using set_non_velocity() with gammamax=10.

    int returndonut;
    // set as in fluid frame
    FTYPE pradffortho[NPR],pradfforthoatm[NPR];
    if(1){
      
      if(PRAD0>=0){
        // choose atmosphere level of radiation
        pradfforthoatm[PRAD0]=RADNT_ERADATMMIN*pow(r/RADNT_ROUT,-2.5);
        pradfforthoatm[PRAD1]=0;
        pradfforthoatm[PRAD2]=0;
        pradfforthoatm[PRAD3]=0;
        
        // copy as backup or atmosphere
        pradffortho[PRAD0]=pradfforthoatm[PRAD0];
        pradffortho[PRAD1]=pradfforthoatm[PRAD1];
        pradffortho[PRAD2]=pradfforthoatm[PRAD2];
        pradffortho[PRAD3]=pradfforthoatm[PRAD3];
      }


      // So donut already has ambient in whichvel whichcoord, now get donut
      if(WHICHPROBLEM==RADDONUT){

        // donut expects fluid frame values in pr.  JCM sets as output in fluid frame so use same conversion below.
        if(PRAD0>=0){
          pr[PRAD0]=pradffortho[PRAD0];
          pr[PRAD1]=pradffortho[PRAD1];
          pr[PRAD2]=pradffortho[PRAD2];
          pr[PRAD3]=pradffortho[PRAD3];
        }
        
        // ADD DONUT
        returndonut=get_full_rtsolution(whichvel,whichcoord,RADDONUT_OPTICALLYTHICKTORUS, pr,X,V,&ptrgeomrad);

        //        if(i==0 && j==0 && myid==0){
        //          dualfprintf(fail_file,"whichvel=%d whichcoord=%d return=%d\n",*whichvel,*whichcoord,returndonut);
        //          dualfprintf(fail_file,"rho=%g u=%g u1=%g u3=%g\n",pr[RHO],pr[UU],pr[U1],pr[U3]);
        //        }

        if(PRAD0>=0){
          // donut returns fluid frame orthonormal values for radiation in pp
          pradffortho[PRAD0]=pr[PRAD0];
          pradffortho[PRAD1]=pr[PRAD1];
          pradffortho[PRAD2]=pr[PRAD2];
          pradffortho[PRAD3]=pr[PRAD3];
          
          //          if(i==21 && j==3 && k==0){
          //            dualfprintf(fail_file,"ZOOM: CHECK: ijk=%d %d %d : %g %g %g %g\n",i,j,k,pradffortho[PRAD0],pradffortho[PRAD1],pradffortho[PRAD2],pradffortho[PRAD3]);
          //          }
        }

    
      }


      if(PRAD0>=0){

        //        if(i==21 && j==3 && k==0){
        //          PLOOP(pliter,pl) dualfprintf(fail_file,"ZOOM: CHECKPREpradfforlab: pl=%d pr=%21.15g\n",pl,pr[pl]);
        //        }

        // Transform these fluid frame E,F^i to lab frame coordinate basis primitives
        prad_fforlab(whichvel, whichcoord, FF2LAB, i,j,k,CENT,ptrgeomrad, pradffortho, pr, pr);

        //        if(i==21 && j==3 && k==0){
        //          PLOOP(pliter,pl) dualfprintf(fail_file,"ZOOM: CHECKPOSTpradfforlab: pl=%d pr=%21.15g\n",pl,pr[pl]);
        //        }
        
        if(debugfail>=2 && returndonut==0 && pradfforthoatm[PRAD0]>pradffortho[PRAD0]){
          dualfprintf(fail_file,"WARNING: Torus radiation pressure below atmosphere.\n");
          dualfprintf(fail_file,"CHECKPOST: ijk=%d %d %d : %g %g %g %g\n",i,j,k,pr[PRAD0],pr[PRAD1],pr[PRAD2],pr[PRAD3]);
        }
      }


#if(0)
      // report zamo
      FTYPE prreport[NPR];
      set_zamo_velocity(*whichvel,ptrgeomrad,prreport);
      dualfprintf(fail_file,"ZAMO: ijk=%d %d %d : %g %g %g : fluid: %g %g %g\n",i,j,k,prreport[U1],prreport[U2],prreport[U3],pr[U1],pr[U2],pr[U3]);
      int jj,kk;
      DLOOP(jj,kk) dualfprintf(fail_file,"gn%d%d=%26.20g\n",jj+1,kk+1,ptrgeomrad->gcon[GIND(jj,kk)]);
      DLOOP(jj,kk) dualfprintf(fail_file,"gv%d%d=%26.20g\n",jj+1,kk+1,ptrgeomrad->gcov[GIND(jj,kk)]);
#endif

    }
    else{


      if(PRAD0>=0){
        // like latest koral that assumes radiation frame is zamo and RADNT_ERADATMMIN is actually in that frame, so no transformation for E
        // So this is somewhat inconsistent with boundary conditions
        // KORAL:
        pr[PRAD0] = RADNT_ERADATMMIN; // assumed as lab-frame ZAMO frame value
        set_zamo_velocity(*whichvel,ptrgeomrad,&pr[URAD1-U1]); // only sets URAD1-URAD3 to zamo
      }
    }

    //    dualfprintf(fail_file,"returning: whichvel=%d whichcoord=%d\n",*whichvel,*whichcoord);

    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADCYLBEAM || WHICHPROBLEM==RADCYLBEAMCART){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3;

    
    // as in koral:
    pr[RHO]=1.0;
    pr[UU]=0.1;
    pr[U1]=0.0;
    pr[U2]=0.0;
    pr[U3]=0.0;

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
    
    // radiation primitives directly
    pr[URAD0]=0.0001;
    pr[URAD1]=0.0;
    pr[URAD2]=0.0;
    pr[URAD3]=0.0;


    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    return(0);
  }

  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADCYLJET && (RADCYLJET_TYPE==1 || RADCYLJET_TYPE==4)){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3;

#define fsx(x) ((x)>0 ? exp(-1.0/(x)) : 0.0)
#define gsx(x) (fsx(x)/(fsx(x) + fsx(1.0-(x))))
#define stepfunction(x) (1.0 - gsx((x)-0.5))
#define stepfunctionab(x,a,b) (((a)-(b))*stepfunction(x) + (b))   
#define stepfunctionab2(x,a,b) (((a)-(b))*stepfunction(x)*stepfunction(x) + (b))   
 
    FTYPE rhojet=1E-5;
    FTYPE E0perRho0;
    FTYPE Rho0;
    FTYPE Ehatjet;
    FTYPE Ehatstar;

    FTYPE vz0=0.5; // 0.99;
    FTYPE vr0=0;

    FTYPE L0;
    FTYPE E0;
    FTYPE Ehat0;
    FTYPE Fr0;

    if(RADCYLJET_TYPE==1){
      Rho0=1.0*rhojet;
      E0perRho0=1E-6;
      L0=0;
      E0=E0perRho0*Rho0*rhojet;
      Ehat0=E0+vr0*L0;
      Fr0=-L0-vr0*E0;
      Ehatjet=Ehat0*1E-2;
      Ehatstar=Ehat0;
    }
    else if(RADCYLJET_TYPE==4){
      Rho0=1000.0*rhojet;
      E0perRho0=1E2;
      L0=0;
      E0=E0perRho0*Rho0*rhojet;
      Ehat0=E0+vr0*L0;
      Fr0=-L0-vr0*E0;
      Ehatjet=Ehat0*0.01;
      Ehatstar=Ehat0;
    }


    FTYPE r0=RADNT_MINX;
    FTYPE r1=RADNT_MAXX;

    pr[RHO]=stepfunctionab(r,rhojet,Rho0);
    pr[U1]=vr0*stepfunctionab(r,1.0,0.001); // R
    pr[U2]=vz0*stepfunctionab(r,1,0.001); // z
    pr[U3]=0.0; // phi

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
    
    // radiation primitives directly
    pr[URAD0]=Ehat0*stepfunctionab(r,Ehatjet/Ehat0,Ehatstar/Ehat0); // not quite right.
    pr[URAD1]=pr[U1];
    pr[URAD2]=pr[U2];
    pr[URAD3]=pr[U3];

    // ensure total pressure is constant by varying gas temperature, but assume LTE at t=0
    // P = arad T^4 + rho*T  = arad T^4 + (gam-1)*u = Ehat0 + (gam-1)*u ->
    // u = (P - Ehat0)/(gam-1)
    // assume rad pressure in star balances 
    if(RADCYLJET_TYPE==1){
      FTYPE ptot=1.1*(4.0/3.0-1.0)*Ehatstar;
      pr[UU]=( ptot - (4.0/3.0-1.0)*pr[URAD0])/(gam-1.0);
    }
    else{
      pr[UU]=0.1*pr[URAD0];
    }

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    // THINGS TO TRY:

    return(0);
  }
  else if(WHICHPROBLEM==RADCYLJET && (RADCYLJET_TYPE==2||RADCYLJET_TYPE==3)){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3;

    FTYPE rhojet=0.1;
    FTYPE vz0=0.5; // 0.99;
    FTYPE vr0=0;
    FTYPE Ehat0=0.1*rhojet;

    FTYPE Ehatjet=Ehat0*1E-2;

    FTYPE r0=RADNT_MINX;
    FTYPE r1=RADNT_MAXX;

    pr[RHO]=rhojet;
    pr[U1]=vr0;
    pr[U2]=vz0;
    pr[U3]=0.0; // phi

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
    
    // radiation primitives directly
    pr[URAD0]=Ehatjet;
    pr[URAD1]=pr[U1];
    pr[URAD2]=pr[U2];
    pr[URAD3]=pr[U3];

    FTYPE ptot=1.5*(4.0/3.0-1.0)*Ehatjet;
    pr[UU]=( ptot - (4.0/3.0-1.0)*pr[URAD0])/(gam-1.0);

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    return(0);
  }
  else if(WHICHPROBLEM==RADCYLJET && (RADCYLJET_TYPE==5)){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3;

    FTYPE rhojet=0.1*100.0;
    FTYPE vz0=0.5; // 0.99;
    FTYPE vr0=0;
    FTYPE Ehat0=0.1*rhojet;

    FTYPE Ehatjet=Ehat0*1E-2;

    FTYPE r0=RADNT_MINX;
    FTYPE r1=RADNT_MAXX;

    pr[RHO]=rhojet;
    pr[U1]=vr0;
    pr[U2]=vz0;
    pr[U3]=0.0; // phi

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
    
    // radiation primitives directly
    pr[URAD0]=Ehatjet;
    pr[URAD1]=pr[U1];
    pr[URAD2]=pr[U2];
    pr[URAD3]=pr[U3];

    pr[UU]=0.1*pr[RHO];

    RADCYLJET_RHOJET=pr[RHO];
    RADCYLJET_UJET=pr[UU];
    RADCYLJET_EHATJET=pr[URAD0];
    RADCYLJET_VRSTAR=pr[U1];


    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    return(0);
  }
  /*************************************************/
  /*************************************************/
  else if(WHICHPROBLEM==RADCYLJET && (RADCYLJET_TYPE==6)){
    FTYPE r,th,ph;
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    ph=V[3];

    *whichcoord=MCOORD; // not BLCOORD, in case setting for inside horizon too when MCOORD=KSCOORDS or if using SPCMINKMETRIC
    *whichvel=VEL3;

    FTYPE Rho0=1E-5*100.0*100.0/79.4*40.0;
    pr[RHO]=Rho0;

    pr[U1]=0.0;
    pr[U2]=0.0;
    pr[U3]=0.0;

    
    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;
    
    // radiation primitives directly
    pr[URAD1]=pr[U1];
    pr[URAD2]=pr[U2];
    pr[URAD3]=pr[U3];

    // ensure thermal equilibrium in star
    FTYPE Tstar;
    if(WHICHRADSOURCEMETHOD==SOURCEMETHODNONE) Tstar=1510.0*1.0E7/TEMPBAR/0.8;
    else Tstar=6.0E7/TEMPBAR;
    // P = (arad/3)T^4 + rho T
    pr[UU]=u_rho0_T_simple(i,j,k,CENT,pr[RHO],Tstar);
    pr[URAD0]=calc_LTE_EfromT(Tstar);

    if(1){
      static int firsttime=1;
      if(firsttime){
        // cs2 ~ gam*Ptot/rho and need vz0>cs.
        // look at sqrt(cs2tot) in SM to check or obtain here.
        FTYPE ptot=((gam-1.0)*pr[UU] + (4.0/3.0-1.0)*pr[URAD0]);
        FTYPE ptrue;
        if(WHICHRADSOURCEMETHOD==SOURCEMETHODNONE) ptrue=(gam-1.0)*pr[UU];
        else ptrue=ptot;
        FTYPE gamptot=(gam*(gam-1.0)*pr[UU] + (4.0/3.0)*(4.0/3.0-1.0)*pr[URAD0]);
        FTYPE cs2 = gamptot/pr[RHO];
        dualfprintf(fail_file,"STAR: ptrue=%21.15g ptot=%21.15g ptrue/rho=%21.15g cs2=%21.15g cs=%21.15g\n",ptrue,ptot,ptrue/pr[RHO],cs2,sqrt(cs2));
        firsttime=0;
      }
    }

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    // really "star" values
    RADCYLJET_RHOJET=pr[RHO];
    RADCYLJET_UJET=pr[UU];
    RADCYLJET_EHATJET=pr[URAD0];
    RADCYLJET_TEMPJET=Tstar;
    RADCYLJET_VRSTAR=pr[U1];

    static int firsttime=1;
    if(myid==0&&firsttime==1&&i==0&&j==0&&k==0){
      firsttime=0;
      FILE *fstar;
      if((fstar=fopen("star.txt","wt"))==NULL){
        dualfprintf(fail_file,"Couldn't open star.txt\n");
        myexit(1);
      }
      else{
        fprintf(fstar,"%d %21.15g %21.15g %21.15g %21.15g %21.15g\n",RADCYLJET_TYPE,RADCYLJET_RHOJET,RADCYLJET_UJET,RADCYLJET_EHATJET,RADCYLJET_TEMPJET,RADCYLJET_VRSTAR);
        fclose(fstar);
      }
    }

    return(0);
  }

  return(0);
}

// analytical solution for RADDONUT donut
// get full radiative donut solution
// input pp has backup (e.g. atmosphere) values in *whichvel, *whichcoord with geometry ptrgeom
// input pp[PRAD0-PRAD3] are fluid frame orthonormal values
// RETURNS: pp and can change whichvel and whichcoord and ptrptrgeom (these must all 3 be consistent with the returned pp)
static int get_full_rtsolution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom)
{
  int jj,kk;
  int pliter,pl;
  int i=(*ptrptrgeom)->i;
  int j=(*ptrptrgeom)->j;
  int k=(*ptrptrgeom)->k;
  int loc=(*ptrptrgeom)->p;
  FTYPE r=V[1];
  FTYPE E,Fx,Fy,Fz;
  FTYPE ppback[NPR];
  PLOOP(pliter,pl) ppback[pl]=pp[pl]; // for initial backup, used to see if torus below atmosphere in donut
  // whichvelback and whichcoordback hold values before any make_nonrt2rt_solution call
  int whichvelback=*whichvel;
  int whichcoordback=*whichcoord;

  // get actual BL-coords in case ptrgeom not, because we need BL-coords to get sizes of cells in finite differences below
  int whichcoordreal=MCOORD;
  int getprim=0;
  struct of_geom geombldontuse;
  struct of_geom *ptrgeombl=&geombldontuse;
  gset(0,whichcoordreal,i,j,k,ptrgeombl); // for local differences using V.


  // get donut with whichvel / whichcoord / ptrptrgeom (or changed inside)
  int anret=make_nonrt2rt_solution(whichvel, whichcoord, opticallythick, pp,X,V,ptrptrgeom);
  if(DOWALDDEN){
    anret=1; // never inside donut for Wald
  }

  //  dualfprintf(fail_file,"anret=%d pp[RHO]=%g\n",anret,pp[RHO]);

  // see if inside donut solution or outside
  if(anret==0){ // then inside donut
    int anretlocal=0;

    ///////////////// STAGE1
    if(PRAD0>=0){
      E=pp[PRAD0];
      Fx=pp[PRAD1];
      Fy=pp[PRAD2];
      Fz=pp[PRAD3];
    }



    ////////////////// STAGE2

    //estimating F = -1/chi E,i
    //http://www.astro.wisc.edu/~townsend/resource/teaching/astro-310-F08/23-rad-diffusion.pdf
    // -(1/chi)dPrad/dx = Frad/c
    FTYPE kappa,kappaes,chi;
    // use of V assumes user knows which coordinates they are in (e.g. r,th,ph vs. x,y,z)
    FTYPE Tg=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
    FTYPE bsq,B;
    bsq_calc(pp,*ptrptrgeom,&bsq);
    B=sqrt(bsq);
    chi=
      calc_kappa_user(pp[RHO],B,Tg,Tg,V[1],V[2],V[3])
      +
      calc_kappaes_user(pp[RHO],calc_PEQ_Tfromurho(pp[UU],pp[RHO]),V[1],V[2],V[3]);

    FTYPE Vtemp1[NDIM],Vtemp2[NDIM];
    FTYPE Xtemp1[NDIM],Xtemp2[NDIM];
    FTYPE pptemp[NPR],E1=0,E2=0;
    getprim=(whichcoordback==PRIMECOORDS); // as consistent with pptemp being fed-in as backup
    int anretmin=0;
    struct of_geom geomtdontuse;
    struct of_geom *ptrgeomt=&geomtdontuse;

    ////////////////// STAGE2A

    //r dimension
    Xtemp1[0]=X[0];
    Xtemp1[1]=1.01*X[1];
    Xtemp1[2]=1.0*X[2];
    Xtemp1[3]=1.0*X[3];
    bl_coord(Xtemp1,Vtemp1); // only needed for Vtemp1[1] for radius in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp1,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp1,Vtemp1,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){
      E1=pptemp[PRAD0]; // fluid frame E
    }

    Xtemp2[0]=X[0];
    Xtemp2[1]=.99*X[1];
    Xtemp2[2]=1.0*X[2];
    Xtemp2[3]=1.0*X[3];
    bl_coord(Xtemp2,Vtemp2); // only needed for Vtemp2[1] for radius in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp2,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp2,Vtemp2,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){ 
      E2=pptemp[PRAD0]; // fluid frame E
    }

    // fluid frame Fx
    //    Fx=(E2-E1)/(.02*V[1]*(ptrgeombl->gcov[GIND(1,1)]))/chi/3.;
    Fx=-THIRD*(E2-E1)/((Vtemp2[1]-Vtemp1[1])*sqrt(fabs(ptrgeombl->gcov[GIND(1,1)])))/chi;
    //    dualfprintf(fail_file,"E2=%g E1=%g Fx=%g\n",E2,E1,Fx);

    ////////////////// STAGE2B

    //th dimension
    Xtemp1[0]=X[0];
    Xtemp1[1]=1.0*X[1];
    Xtemp1[2]=1.01*X[2];
    Xtemp1[3]=1.0*X[3];
    bl_coord(Xtemp1,Vtemp1); // only needed for Vtemp1[2] for theta in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp1,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp1,Vtemp1,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){ 
      E1=pptemp[PRAD0]; // fluid frame E1
    }

    Xtemp2[0]=X[0];
    Xtemp2[1]=1.0*X[1];
    Xtemp2[2]=0.99*X[2];
    Xtemp2[3]=1.0*X[3];
    bl_coord(Xtemp2,Vtemp2); // only needed for Vtemp2[2] for theta in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp2,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp2,Vtemp2,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){ 
      E2=pptemp[PRAD0]; // fluid frame E2
    }

    // fluid frame Fy
    //    Fy=(E2-E1)/(.02*V[2]*(ptrgeombl->gcov[GIND(2,2)]))/chi/3.;
    Fy=-THIRD*(E2-E1)/((Vtemp2[2]-Vtemp1[2])*sqrt(fabs(ptrgeombl->gcov[GIND(2,2)])))/chi;
    //dualfprintf(fail_file,"E2=%g E1=%g Fy=%g\n",E2,E1,Fy);

    ////////////////// STAGE2C

    //ph dimension
    //    Fz=0.; // fluid frame Fz
    Xtemp1[0]=X[0];
    Xtemp1[1]=1.0*X[1];
    Xtemp1[2]=1.0*X[2];
    Xtemp1[3]=1.01*X[3];
    bl_coord(Xtemp1,Vtemp1); // only needed for Vtemp1[3] for phi in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp1,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp1,Vtemp1,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){ 
      E1=pptemp[PRAD0]; // fluid frame E1
    }

    Xtemp2[0]=X[0];
    Xtemp2[1]=1.0*X[1];
    Xtemp2[2]=1.0*X[2];
    Xtemp2[3]=0.99*X[3];
    bl_coord(Xtemp2,Vtemp2); // only needed for Vtemp2[3] for phi in make_nonrt2rt_solution()
    gset_X(getprim,whichcoordback,i,j,k,NOWHERE,Xtemp2,ptrgeomt);
    PLOOP(pliter,pl) pptemp[pl]=ppback[pl]; // ppback that holds unchanged pp, not modified by previous make_nonrt2rt_solution() call.

    anretlocal=make_nonrt2rt_solution(&whichvelback, &whichcoordback, opticallythick, pptemp,Xtemp2,Vtemp2,&ptrgeomt);
    if(anretlocal<0) anretmin=-1;
    if(PRAD0>=0){ 
      E2=pptemp[PRAD0]; // fluid frame E2
    }

    // fluid frame Fz
    Fz=-THIRD*(E2-E1)/((Vtemp2[3]-Vtemp1[3])*sqrt(fabs(ptrgeombl->gcov[GIND(3,3)])))/chi;
    //dualfprintf(fail_file,"E2=%g E1=%g Fz=%g\n",E2,E1,Fz);



    ////////////////// STAGE3

    if(anretmin<0){ // then ended up with one point outside solution
      Fx=Fy=Fz=0.;
      //      dualfprintf(fail_file,"OUTSIDE: ijk=%d %d %d\n",i,j,k);
    }
    else{ // then inside solution
      //      dualfprintf(fail_file,"INSIDE: ijk=%d %d %d\n",i,j,k);
      FTYPE MAXFOE;
      
      MAXFOE=0.70; // KORALTODO: SUPERGODMARK: Kraken failed with this set to 0.99 where hit 3vel calculation with disc=+1E-16.  So just marginal beyond v=c, but with 0.99 why did that happen?

      //      if(r>1E3) MAXFOE=0.1;


      FTYPE Fl=sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
      if(Fl>MAXFOE*E){
        Fx=Fx/Fl*MAXFOE*E;
        Fy=Fy/Fl*MAXFOE*E;
        Fz=Fz/Fl*MAXFOE*E;
        FTYPE Flnew=sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
        if(PRODUCTION==0) dualfprintf(fail_file,"limited: Fl=%g E=%g Fx=%g Fy=%g Fz=%g : Flnew=%g\n",Fl,E,Fx,Fy,Fz,Flnew);
      }
    }


    if(PRAD0>=0){ 
      //saving ff values to pp[] (so any function using this function should know pp has fluid frame orthonormal values in pp[PRAD0-PRAD3] as was in the input as well.
      pp[PRAD1]=Fx;
      pp[PRAD2]=Fy;
      pp[PRAD3]=Fz;
    }



  }// end if adding donut


  if(anret==0) return(anret);
  else return(-1);

}


// convert non-radiative solution to radiative one by splitting total pressure up into gas + radiation pressure with optical depth corrections
// expects pp[PRAD0-PRAD3] to be fluid frame orthonormal, while pp[U1-U3] is ptrgeom whichvel whichcoord lab frame value (doesn't change U1-U3)
// returns: pp and can change whichvel and whichcoord
static int make_nonrt2rt_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V, struct of_geom **ptrptrgeom)
{
  // get location
  int i=(*ptrptrgeom)->i;
  int j=(*ptrptrgeom)->j;
  int k=(*ptrptrgeom)->k;
  int p=(*ptrptrgeom)->p;



  /////
  //
  // choose non-radiative solution to start from
  //
  /////
  // total torus pressure
  FTYPE pt;
  int usingback=donut_analytical_solution(whichvel, whichcoord, opticallythick, pp, X, V, ptrptrgeom, &pt);
  //int usingback=process_solution(whichvel, whichcoord, opticallythick, pp, X, V, ptrptrgeom, &pt);
  
  // assign to names
  FTYPE rho=pp[RHO];
  FTYPE uint=pp[UU];
  FTYPE E=pp[URAD0];
  FTYPE Fx=pp[URAD1];
  FTYPE Fy=pp[URAD2];
  FTYPE Fz=pp[URAD3];

  //  dualfprintf(fail_file,"rho=%g usingback=%d\n",rho,usingback);

  /////////////////////////////////////
  //
  // get actual temperature and separate internal energies for gas and radiation
  //
  /////////////////////////////////////
  FTYPE P,aaa,bbb;
  P=pt; // torus total pressure
  FTYPE Teq=calc_PEQ_Tfromurho(uint,rho);


  // initial equilibrium optically thick version of gas temperature
  FTYPE Tgas=Teq;

  // loop over to get correct gas temperature since kappa(Tgas) that affects Tgas.
#define NUMTgasITERS 4 // for now, was 4.
  int iter=0;
  for(iter=0;iter<NUMTgasITERS;iter++){

    // 2-stream approximation for pressure
    FTYPE kappaabs,kappaes,kappatot;
    // use of V assumes user knows which coordinates they are in (e.g. r,th,ph vs. x,y,z)
    FTYPE bsq,B;
    bsq_calc(pp,*ptrptrgeom,&bsq);
    B=sqrt(bsq);
    kappaabs=calc_kappa_user(rho,B,Tgas,Tgas,V[1],V[2],V[3]);
    kappaes=calc_kappaes_user(rho,Tgas,V[1],V[2],V[3]);
    kappatot=kappaabs+kappaes;

    // fake integral of opacity as if opacity roughly uniform and heading out away in angular distance at fixed angluar extent.  Could integrate to pole, but harder with MPI.
    FTYPE r=V[1];
    FTYPE th=V[2];
    FTYPE dl;
    dl = r*MAX(h_over_r-fabs(th-M_PI*0.5),0.0);

    FTYPE dxdxp[NDIM][NDIM];
    if(p==NOWHERE) dxdxprim_ijk(i, j, k, CENT, dxdxp); // don't use p, use CENT, because nowhere when differencing
    else dxdxprim_ijk(i, j, k, p, dxdxp); // don't use p, use CENT, because nowhere when differencing
    dl = r*MAX(dl,dx[2]*dxdxp[2][2]);

    FTYPE tautot,tauabs;
    tautot = kappatot*dl+SMALL;
    tauabs = kappaabs*dl+SMALL;
    FTYPE gtau = (tautot/2.0 + 1.0/sqrt(3.0))/(tautot/2.0 + 1.0/sqrt(3.0) + 1.0/(3.0*tauabs));
    
    // dualfprintf(fail_file,"dl=%g tautot=%g tauabs=%g\n",dl,tautot,tauabs);



    if(usingback==0){// && gtau>1E-10){
      gtau=1.0; // TESTING
      opticallythick=1;
      //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
      aaa=(4.0/3.0-1.0)*ARAD_CODE*gtau; // Prad=(gamma-1)*urad*g[tau] and gamma=4/3 and urad=arad*T^4*g[tau]
      bbb=rho;

      //      dualfprintf(fail_file,"aaa=%g bbb=%g P=%g gtau=%g\n",aaa,bbb,P,gtau);


      FTYPE naw1=cbrt(9.*aaa*Power(bbb,2.) - Sqrt(3.)*Sqrt(27.*Power(aaa,2.)*Power(bbb,4.) + 256.*Power(aaa,3.)*Power(P,3.)));
      Tgas=-Sqrt((-4.*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2.,0.3333333333333333)*Power(3.,0.6666666666666666)*aaa))/2. + Sqrt((4.*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2.,0.3333333333333333)*Power(3.,0.6666666666666666)*aaa) + (2.*bbb)/(aaa*Sqrt((-4.*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2.,0.3333333333333333)*Power(3.,0.6666666666666666)*aaa))))/2.;
      // cap Tgas for optically thin regions where density can be very low
      if(Tgas>1E12/TEMPBAR) Tgas=1E12/TEMPBAR;
      else{
        if(fabs(rho*Tgas + aaa*pow(Tgas,4.0)*gtau - P)>1E-5){
          dualfprintf(fail_file,"2: NOT right equation: rho=%26.21g Tgas=%26.21g aaa=%26.21g Ptried=%26.21g P=%26.21g\n",rho,Tgas,aaa,rho*Tgas + aaa*pow(Tgas,4.),P);
        }
      }

      if(!isfinite(Tgas)) Tgas=TEMPMIN;

      //dualfprintf(fail_file,"Tgas=%g\n",Tgas);


    }
    else{

      //      gtau=SMALL; // say atmosphere is always optically thin, to avoid overly radiative atmosphere for low rest-mass densities.
      opticallythick=0;
      Tgas = P/rho;
    }



    if(opticallythick){
      uint=calc_PEQ_ufromTrho(Tgas,rho);
      // overwrite uint only if actually optically thick, otherwise stick with original pressure that assumed pure gas-based pressure

      E=calc_LTE_EfromT(Tgas); // assume LTE, so thermal equilibrium between gas and radiation
      Fx=Fy=Fz=0.;
    }
    else{
      // else keep E that is background E.
      // stick with uint -- GODMARK: WHY?
      // GODMARK: Should set uint more generally, not just if(opticallythick) above
    }

    //    FTYPE rhoatm=RADNT_RHOATMMIN*pow(r/RADNT_ROUT,-1.5);
    //    dualfprintf(fail_file,"rhodonut4=%g uint=%g Tgas=%g gtau=%g tautot=%g tauabs=%g dl=%g E=%g usingback=%d ratiorho=%g\n",rho,uint,Tgas,gtau,tautot,tauabs,dl,E,usingback,rho/rhoatm);
  }


  // if doing radiation, then modify radiation using this result
  pp[RHO]=rho;
  pp[UU]=uint;
  if(PRAD0>=0){ 
    // fluid frame orthonormal values for radiation
    pp[PRAD0]=E;
    pp[PRAD1]=Fx;
    pp[PRAD2]=Fy;
    pp[PRAD3]=Fz;
  }

  if(usingback){
    return(-1); // tells to not form radiative flux for atmosphere
  }
  return 0;
}




// analytical solution for RADDONUT donut
// expects pp[PRAD0-PRAD3] to be fluid frame orthonormal, while pp[U1-U3] is ptrgeom whichvel whichcoord lab frame value (doesn't change U1-U3)
// RETURNS: pp, uTptr, ptptr and can return whichvel and whichcoord and ptrptrgeom
static int donut_analytical_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom, FTYPE *ptptr)
{
  int usingback=0;
  FTYPE Vphi=0.0,Vr=0.0,Vh=0.0;
  FTYPE D,W,uT=1.0,uphi,uPhi,rho,ucon[NDIM],uint,E,Fx,Fy,Fz;
  //  FTYPE rho,uint,uT=1.0,E,Fx,Fy,Fz;

  // set backup atmosphere value for primitives
  int pliter,pl;
  FTYPE ppback[NPR];
  PLOOP(pliter,pl) ppback[pl]=pp[pl];

  // get location (NOTE: not affected by whichcoord)
  FTYPE r=V[1];
  FTYPE th=V[2];
  FTYPE Rcyl=fabs(r*sin(th));
  int i=(*ptrptrgeom)->i;
  int j=(*ptrptrgeom)->j;
  int k=(*ptrptrgeom)->k;
  int p=(*ptrptrgeom)->p;


  // choose \gamma ideal gas constant for torus (KORALTODO: Could interpolate based upon some estimate of \tau)
  FTYPE gamtorus;
  if(opticallythick==1) gamtorus=4.0/3.0; // then should be as if gam=4/3 so radiation supports torus properly at t=0
  else gamtorus=gam;
  // total torus pressure that will be distributed among gas and radiation pressure
  FTYPE pt;




  // get density and velocity for torus
  if(RADNT_DONUTTYPE==DONUTTHINDISK3){

    int mycoords=BLCOORDS;
    if(*whichcoord!=mycoords){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_X(getprim,mycoords,i,j,k,NOWHERE,X,*ptrptrgeom);
    }
    *whichcoord=mycoords;
    *whichvel=VEL3;




    /* region outside disk */
    Rhor=rhor_calc(0);
    Risco=rmso_calc(PROGRADERISCO);
    R = MAX(Rhor,r*sin(th)) ;

#define NTROUT (70.0)
    // NT73
    FTYPE rp;
    FTYPE Rp=R;
    if(r>NTROUT){
      rp=NTROUT;
      Rp=rp*sin(th);
    }
    else if(r>=1.5*Risco && r<NTROUT){
      rp=r;
      Rp=R;
    }
    else{
      rp=1.5*Risco;
      Rp=MAX(Rhor,rp*sin(th)) ;
    }

    // alphavisc=0.1 leads to Mdot=0.16MdotEdd (subeddd)
    // alphavisc=1.0 leads to Mdot=2.8MdotEdd (subeddc)

    FTYPE alphavisc=0.4;
    FTYPE m=10.0; // Mbh in Msun
    FTYPE Mx=m;
    FTYPE mdot=7.33; // Mdot per Ledd/c^2
    FTYPE mdotperledd=mdot;
    FTYPE y=pow(rp,0.5);
    FTYPE AA=1.0+a*a*pow(y,-4.0)+2.0*a*a*pow(y,-6.0);
    FTYPE BB=1.0+a*pow(y,-3.0);
    FTYPE CC=1-3.0*pow(y,-2.0)+2.0*a*pow(y,-3.0);
    FTYPE DD=1-2.0*pow(y,-2.0)+a*a*pow(y,-4.0);
    FTYPE EE=1+4.0*a*a*pow(y,-4.0)-4.0*a*a*pow(y,-6.0)+3.0*a*a*a*a*pow(y,-8.0);


    FTYPE QQ=(1.0000000000000002*(1 + 0.5/Power(rp,1.5))*
     (-2.0574261905426465 + 1.*Sqrt(rp) - 
       0.5160444431189778*
        Log(1.9035389107262222*
          (-1.5320888862379562 + 1.*Sqrt(rp))) + 
       0.07635182233306995*
        Log(0.5847509232408147*
          (-0.34729635533386044 + 1.*Sqrt(rp))) + 
       1.1896926207859089*
        Log(0.25401267427810215*
          (1.8793852415718166 + 1.*Sqrt(rp))) - 
       0.75*Log(0.4860441675121526*Sqrt(rp))))/
   (Sqrt(1 + 1./Power(rp,1.5) - 2.999999999999999/rp)*Sqrt(rp));

    FTYPE Sigma1=(5.000000000000001*(1 + 0.1874999999999999/Power(rp,4) - 
       0.9999999999999999/Power(rp,3) + 
       0.9999999999999998/Power(rp,2))*
     Power(1 + 0.5/Power(rp,1.5),2)*
     (1 + 1./Power(rp,1.5) - 2.999999999999999/rp)*Power(rp,2))/
   (alphavisc*mdotperledd*
     Power(1 + 0.49999999999999994/Power(rp,3) + 
       0.24999999999999994/Power(rp,2),2)*
     (-2.0574261905426465 + 1.*Sqrt(rp) - 
       0.5160444431189778*
        Log(1.9035389107262222*
          (-1.5320888862379562 + 1.*Sqrt(rp))) + 
       0.07635182233306995*
        Log(0.5847509232408147*
          (-0.34729635533386044 + 1.*Sqrt(rp))) + 
       1.1896926207859089*
        Log(0.25401267427810215*
          (1.8793852415718166 + 1.*Sqrt(rp))) - 
       0.75*Log(0.4860441675121526*Sqrt(rp))));

    FTYPE H1=(100000.00000000001*mdotperledd*Mx*
     Power(1 + 0.49999999999999994/Power(rp,3) + 
       0.24999999999999994/Power(rp,2),2)*
     (-2.0574261905426465 + 1.*Sqrt(rp) - 
       0.5160444431189778*
        Log(1.9035389107262222*
          (-1.5320888862379562 + 1.*Sqrt(rp))) + 
       0.07635182233306995*
        Log(0.5847509232408147*
          (-0.34729635533386044 + 1.*Sqrt(rp))) + 
       1.1896926207859089*
        Log(0.25401267427810215*
          (1.8793852415718166 + 1.*Sqrt(rp))) - 
       0.75*Log(0.4860441675121526*Sqrt(rp))))/
   ((1 + 0.1874999999999999/Power(rp,4) - 
       0.9999999999999999/Power(rp,3) + 
       0.9999999999999998/Power(rp,2))*
     Power(1 + 0.5/Power(rp,1.5),2)*
     (1 + 0.24999999999999994/Power(rp,2) - 
       1.9999999999999993/rp)*Sqrt(rp));


    // rescale
    Sigma1=Sigma1/MBAR*LBAR*LBAR; // Sigma = g/cm^2
    H1=H1/LBAR; // H=cm

    


    ////////////////////////////
    // Set H : disk height
    FTYPE Hp=H1;
    FTYPE minhor=h_over_r*0.3;
    if(Hp/rp<minhor){
      Hp=minhor*rp; // go to constant H/R so disk not unresolved
    }


    FTYPE z = r*cos(th) ;
    FTYPE zp = rp*cos(th) ;


    //////////////////////////
    // SET DENSITY
    // Sigma = rho*H
    FTYPE rho01;
    rho01=Sigma1/(2.0*Hp); // so Sigma unchanged.
    
    //
    // simple power-law with some vertical distribution
    //
    //rho=RADNT_RHODONUT * exp(-z*z/(2.*H*H))*pow(r,thindiskrhopow);
    rho=rho01*exp(-zp*zp/(2.*Hp*Hp));

    // truncate disk
    if(r>NTROUT){
      rho*=exp(-r/NTROUT)/exp(-1.0);
    }

    //////////////////////////
    // ensure c_s/v_k = (H/R) , which should be violated in the ISCO such that K=P/rho^Gamma=constant but H/R still thins towards the horizon.
    // P = (gamma-1)*u and cs2 = \gamma P/(\rho+u+P)
    FTYPE omegakep=1./(pow(r,1.5) + a);
    FTYPE omegakeprp=1./(pow(rp,1.5) + a);

    // set uint
    FTYPE uintfact = (1.0 + randfact * (ranc(0,0) - 0.5));
    uint = rho*pow(Hp/Rp,2.0)*pow(rp*omegakeprp,2.0)/(gamtorus*(gamtorus-1.0))*uintfact;

    // Set pressure
    // K = P/rho^Gamma 
    FTYPE pint=(gamtorus-1.0)*uint;
    FTYPE KK= pint*pow(rho,-gamtorus);
    // get total pressure
    pt = uint * (gamtorus-1.0); // torus pressure


    //////////////////////////
    // SET
    // set radial velocity
    if(r<1.5*Rhor){
      Vr = 0.0;
      Vh =0.0;
      Vphi = 0.0;
      *whichvel=VELREL4;
      *whichcoord=KSCOORDS;
      usingback=1;
    }
    else{
      Vphi=omegakep; // now 3-vel
      Vr=-fabs(alphavisc*(Hp/rp)*(Hp/rp)*Vphi);
      Vh=0.0;
      *whichvel=VEL3;
      *whichcoord=BLCOORDS;
    }


    // see if should still use backup non-torus values
    if(rho<ppback[RHO]){
      usingback=1;
    }

    //    usingback=0;


    //    dualfprintf(fail_file,"rhodonut1=%g uint=%g\n",rho,uint);

    //    if(fabs(r-Risco)<0.1*Risco && fabs(th-0.5*M_PI)<0.2*h_over_r){
    //      dualfprintf(fail_file,"rhodonut5=%g uint=%g : r,th=%g %g usingback=%d RADNT_RHODONUT=%g rhoisco=%g\n",rho,uint,r,th,usingback,RADNT_RHODONUT,rhoisco);
    //    }


  }
  // get density and velocity for torus
  if(RADNT_DONUTTYPE==DONUTTHINDISK2){

    int mycoords=BLCOORDS;
    if(*whichcoord!=mycoords){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_X(getprim,mycoords,i,j,k,NOWHERE,X,*ptrptrgeom);
    }
    *whichcoord=mycoords;
    *whichvel=VEL3;




    /* region outside disk */
    Rhor=rhor_calc(0);
    if(1){
      Risco=rmso_calc(PROGRADERISCO);
    }
    else{ // SUPERMADNEW
      Risco=10.0;
    }
    R = MAX(Rhor,r*sin(th)) ;

    ////////////////////////////
    // Set H : disk height
    FTYPE H = h_over_r*R;
    FTYPE Hisco = h_over_r*Risco;

    // or try:
    //    FTYPE xrpow=0.5;
    FTYPE xrpow=1.0;
    FTYPE xr=pow(r,xrpow);
    FTYPE xr0=pow(Risco,xrpow);
    FTYPE xr1=pow(1.0,xrpow);
    //    H = 0.1*h_over_r*R + h_over_r*R * (MAX(0.0,xr-xr0)/xr) ;
    //    Hisco = 0.1*h_over_r*Risco + h_over_r*Risco * (MAX(0.0,xr0-xr0)/xr0) ;
    if(r>=Risco){
      H = h_over_r*R ;
    }
    else{
      //      H = h_over_r*R * pow(r/Risco,2.0);
      H = h_over_r*R * pow(r/Risco,.5); // SUPERMADNEW
    }

    FTYPE z = r*cos(th) ;
    FTYPE zisco = Risco*cos(th) ;

    //////////////////////////
    // SET DENSITY
    //
    // simple power-law with some vertical distribution
    //
    FTYPE rhoisco;
    if(1){
      // below Gaussian for isothermal gas, but works better for adiabatic EOS too.
      rhoisco=RADNT_RHODONUT * exp(-zisco*zisco/(2.*Hisco*Hisco))*pow(Risco,thindiskrhopow);
      rho=RADNT_RHODONUT * exp(-z*z/(2.*H*H))*pow(r,thindiskrhopow);
    }
    if(0){
      // for adiabatic EOS, need improved vertical distribution
      FTYPE NN=1.0/(gamtorus-1.0);
      rhoisco=RADNT_RHODONUT * pow(1.0-zisco*zisco/(Hisco*Hisco),NN)*pow(Risco/Risco,thindiskrhopow); if(rhoisco<0.0||!isfinite(rhoisco)) rhoisco=0.0;
      rho=RADNT_RHODONUT * pow(1.0-z*z/(H*H),NN)*pow(r/Risco,thindiskrhopow); if(rho<0.0||!isfinite(rho)) rho=0.0;
    }
    // density shouldn't peak at ISCO, but further out even for H/R=0.05 (see Penna et al. 2010 figure 11).
    FTYPE Rtrans;
    //    Rtrans=2.0*Risco;
    Rtrans=1.0*Risco; // SUPERMADNEW
    if(r<Rtrans){
      //      FTYPE rhopowisco=10.0;
      FTYPE rhopowisco=7.0; // this gives similar result to HD case, although not MHD-turbulence case, but want kinda HD equilibrium as much as possible.
      rho = rho*pow(r/Rtrans,rhopowisco);
    }


    // see if should still use backup non-torus values
    if(rho<ppback[RHO]){
      usingback=1;
    }

    if(DOWALDDEN){
      usingback=1;
    }

    //////////////////////////
    // ensure c_s/v_k = (H/R) , which should be violated in the ISCO such that K=P/rho^Gamma=constant but H/R still thins towards the horizon.
    // P = (gamma-1)*u and cs2 = \gamma P/(\rho+u+P)
    FTYPE omegakep=1./(pow(r,1.5) + a);
    FTYPE omegakepisco=1./(pow(Risco,1.5) + a);
    FTYPE omega,omegaisco=omegakepisco;

    FTYPE sigma=r*r+a*a*cos(th)*cos(th);
    FTYPE gppvsr=sigma + a*a*(1+2*r/sigma)*sin(th)*sin(th);
    FTYPE sigmaisco=Risco*Risco+a*a*cos(th)*cos(th);
    FTYPE gppvsrisco=sigmaisco + a*a*(1+2*Risco/sigmaisco)*sin(th)*sin(th);
    if(1||r>=Risco){ // final HD solution is Keplerian over all radii!
      omega=omegakep;
    }
    else{
      omega=omegaisco*gppvsrisco/gppvsr;
    }
    uPhi = omega  ; // Omega = v^\phi = 3-vel = d\phi/dt
    uPhi*=SLOWFAC;
    //    FTYPE grr=(*ptrptrgeom)->gcov[GIND(3,3)]; // BL
    FTYPE gtt=-(1.0-2.0*r/sigma);
    FTYPE gpp=sin(th)*sin(th)*(sigma+a*a*(1.0+2.0*r/sigma)*sin(th)*sin(th)); // KS
    uT = 1.0/pow(-gtt-uPhi*uPhi*gpp,0.5); // approximate \gamma = 1/sqrt(1-v^2) ~ 1/sqrt(1-vphi^2)
    FTYPE Delta=r*r-2*r+a*a;
    FTYPE AA=pow(r*r+a*a,2.0)-a*a*Delta*sin(th)*sin(th);
    FTYPE gttbl=-(1.0-2.0*r/sigma);
    FTYPE gppbl=AA*sin(th)*sin(th)/sigma;
    FTYPE uTbl = 1.0/pow(-gttbl-uPhi*uPhi*gppbl,0.5); // \gamma = 1/sqrt(1-v^2) ~ 1/sqrt(1-vphi^2) for BL-coords with u^r=u^\theta=0 in BL coords

    //////////////////////////
    // SET
    if(r>=Risco){
      Vphi=uPhi; // now 3-vel
    }
    else{
      //      Vphi=uPhi/pow(r/Risco,4.0)/0.6/0.80;
      Vphi=uPhi;
    }
    // set radial velocity
    //    if(r<Risco){
    if(r<1.1*Rhor){
      Vr = ppback[U1]; // zamo vel
      Vh = ppback[U2]; // zamo vel
    }
    else{
      Vr=Vh=0.0;
    }

    //////////////////////////
    // TEST if velocity gives reasonable u^t
    FTYPE prtest[NPR];
    FTYPE ucontest[NDIM];
    FTYPE others[NUMOTHERSTATERESULTS];
    prtest[U1]=Vr;
    prtest[U2]=Vh;
    prtest[U3]=Vphi;
    //    dualfprintf(fail_file,"BEFORE real ut\n");
    //    int badut=ucon_calc_3vel(prtest,*ptrptrgeom,ucon,others);
    failed=-1; // GLOBAL: indicate not really failure, so don't print out debug info
    int badut=ucon_calc_whichvel(*whichvel,prtest,*ptrptrgeom,ucon,others);
    failed=0; // GLOBAL: reset failure flag since just test
    //if(badut==0) dualfprintf(fail_file,"real ut=%g\n",ucon[TT]);

    // CHECK TEST
    if(usingback==0){
      // CATCH
      //    if(r<Risco || badut){// uT>10.0 || !isfinite(uT) || uTbl>10.0 || !isfinite(uTbl) || uT<1.0-1E-7 || uTbl<1.0-1E-7){
      if(badut==1){
        //      uT=10.0;
        //      uTbl=10.0;
        Vphi=0.0;//uPhi*uT; ///pow(r/Risco,4.0); // Vphi is 4-vel = d\phi/d\tau = u^t \Omega
        // set radial velocity
        Vr = 0. ; // zero zamo vel
        Vh = 0. ; // zero zamo vel
        *whichvel=VELREL4;
        if(*whichcoord==BLCOORDS) *whichcoord=KSCOORDS;
      }
      else if(r<1.1*Rhor){
        //      uT=10.0;
        //      uTbl=10.0;
        Vphi=0.0;//uPhi*uT; ///pow(r/Risco,4.0); // Vphi is 4-vel = d\phi/d\tau = u^t \Omega
        // set radial velocity
        Vr = 0. ; // zero zamo vel
        Vh = 0. ; // zero zamo vel
        *whichvel=VELREL4;
        if(*whichcoord==BLCOORDS) *whichcoord=KSCOORDS;
      }
      else{
        // i.e. keep whichvel, whichcoord the same
        *whichvel=VEL3;
        *whichcoord=BLCOORDS;
      }
      //dualfprintf(fail_file,"uT=%g uTbl=%g\n",uT,uTbl);
    }
    else{
      // i.e. keep whichvel, whichcoord the same
      *whichvel=VEL3;
      *whichcoord=BLCOORDS;
    }

    // OVERRIDE
    if(r>1.2*Rhor){
      *whichvel=VEL3;
      *whichcoord=BLCOORDS;
    }


    // set uint
    FTYPE uintfact = (1.0 + randfact * (ranc(0,0) - 0.5));
    uint = rho*pow(H/R,2.0)*pow(r*omega,2.0)/(gamtorus*(gamtorus-1.0))*uintfact;
    FTYPE uintisco = rhoisco*pow(Hisco/Risco,2.0)*pow(Risco*omegaisco,2.0)/(gamtorus*(gamtorus-1.0))*uintfact;

    // Set pressure
    // K = P/rho^Gamma 
    FTYPE pint=(gamtorus-1.0)*uint;
    FTYPE KK= pint*pow(rho,-gamtorus);
    FTYPE pintisco=(gamtorus-1.0)*uintisco;
    FTYPE KKisco= pintisco*pow(rhoisco,-gamtorus);
    if(r>=Risco){
      pint = pint;
    }
    else{
      // Have constant K inside isco, not uint
      pint = pintisco*pow(rhoisco,-gamtorus)*pow(rho,gamtorus);
    }
    // get total pressure
    pt = uint * (gamtorus-1.0); // torus pressure





    //    dualfprintf(fail_file,"rhodonut1=%g uint=%g\n",rho,uint);

    //    if(fabs(r-Risco)<0.1*Risco && fabs(th-0.5*M_PI)<0.2*h_over_r){
    //      dualfprintf(fail_file,"rhodonut5=%g uint=%g : r,th=%g %g usingback=%d RADNT_RHODONUT=%g rhoisco=%g\n",rho,uint,r,th,usingback,RADNT_RHODONUT,rhoisco);
    //    }


  }
  // get density and velocity for torus
  else if(RADNT_DONUTTYPE==DONUTTHINDISK){

    int mycoords=BLCOORDS;
    if(*whichcoord!=mycoords){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_X(getprim,mycoords,i,j,k,NOWHERE,X,*ptrptrgeom);
    }
    *whichcoord=mycoords;
    *whichvel=VEL3;


    SFTYPE sth, cth;
    //    struct of_geom geomdontuse;
    //    struct of_geom *ptrgeom=&geomdontuse;
    /* for disk interior */
    FTYPE R,H,nz,z,S,cs ;
    SFTYPE rh;
    //    int pl,pliter;


    /* region outside disk */
    Rhor=rhor_calc(0);
    Risco=rmso_calc(PROGRADERISCO);
    R = MAX(Rhor,r*sin(th)) ;
    
    //if(R < rin) { // assume already have atmosphere as background solution
    FTYPE SMALLESTH=1E-15;
      
    //    H = SMALLESTH + h_over_r*R ; // NOTEMARK: Could choose H=constant

    //    FTYPE Rtrans=1.5*Risco;
    FTYPE xrpow=0.5;
    FTYPE xr=pow(r,xrpow);
    FTYPE xr0=pow(Risco,xrpow);
    FTYPE xr1=pow(1.0,xrpow);

    FTYPE HS = h_over_r*R;
    FTYPE HSisco;
    HSisco = h_over_r*Risco;

    FTYPE Hisco;
    if(0){
      H = 0.1*h_over_r*R + h_over_r*R * (MAX(0.0,xr-xr0)/xr) ; // NOTEMARK: Could choose H=constant
      Hisco = 0.1*h_over_r*R;
    }
    else{
      H = HS;
      Hisco = HSisco;
    }

    // fix nz
    nz = nz_func(R) ;
    nz = pow(R,-1.5);
    FTYPE nzisco = nz_func(Risco) ;
    if(nz>1.0) nz=1.0;
    if(r<Risco) nz=nzisco;
    if(!isfinite(nz)) nz=1.0;

    z = r*cos(th) ;
    FTYPE zisco = Risco*cos(th) ;
    S = 1./(HS*HS*nz) ;
    FTYPE Sisco=1./(HSisco*HSisco*nzisco) ;
    cs = H*nz ;
    FTYPE csisco = Hisco*nzisco ;

    FTYPE rho0=RADNT_RHODONUT*(S/sqrt(2.*M_PI*HS*HS)) * exp(-z*z/(2.*H*H))*pow(r,3.0/2.0+thindiskrhopow);
    FTYPE rho0isco=RADNT_RHODONUT*(Sisco/sqrt(2.*M_PI*HSisco*HSisco)) * exp(-zisco*zisco/(2.*Hisco*Hisco))*pow(Risco,3.0/2.0+thindiskrhopow);
    if(r>Risco){
      rho = rho0 * (MAX(0.0,xr-xr1)/xr) ;
    }
    else{
      FTYPE rhopowisco=10.0;
      rho = (rho0isco  * (MAX(0.0,xr-xr1)/xr)) * pow(r/Risco,rhopowisco);
    }
    uint = rho*cs*cs/(gamtorus - 1.) ;
    Vr = 0. ;
    uPhi = 1./(pow(r,1.5) + a) ;
    FTYPE uPhiisco = 1./(pow(Risco,1.5) + a) ;
    // solution for 3-vel


    //    uint = uint * (1. + (1.+3.*(rin/R)*(rin/R))*randfact * (ranc(0,0) - 0.5));
    uint *= (1.0 + randfact * (ranc(0,0) - 0.5));
    uPhi*=SLOWFAC;
    Vphi=uPhi; // approximate
    FTYPE Vphiisco=uPhiisco;

    //////////
    // get total pressure
    pt = uint * (gamtorus-1.0); // torus pressure

    // see if should still use backup non-torus values
    if(rho<ppback[RHO]){
      usingback=1;
    }

    if(DOWALDDEN){
      usingback=1;
    }



    if(r<Risco){
      FTYPE sigma=r*r+a*a*cos(th)*cos(th);
      FTYPE gpp=sigma + a*a*(1+2*r/sigma)*sin(th)*sin(th);
      FTYPE sigmaisco=Risco*Risco+a*a*cos(th)*cos(th);
      FTYPE gppisco=sigmaisco + a*a*(1+2*Risco/sigmaisco)*sin(th)*sin(th);
      Vphi=Vphiisco*gppisco/gpp;
    }

    if(1||(r<1.3*Rhor || r<Risco) && usingback==0){
      *whichvel=VELREL4;
      if(*whichcoord==BLCOORDS) *whichcoord==KSCOORDS;
    }
    else{
      *whichvel=VEL3;
      //      *whichcoord=BLCOORDS;
    }

    //    if(r<Risco){ // TEMP
    //      return -1; //outside disk
    //    }



    //    dualfprintf(fail_file,"nz=%g z=%g S=%g cs=%g H=%g h_over_r=%g\n",nz,z,S,cs,H,h_over_r);

    //    dualfprintf(fail_file,"rhodonut1=%g uint=%g\n",rho,uint);


  }
  else if(RADNT_DONUTTYPE==DONUTOLEK){

    int mycoords=BLCOORDS;
    //    int mycoords=KSCOORDS;
    if(*whichcoord!=mycoords){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_X(getprim,mycoords,i,j,k,NOWHERE,X,*ptrptrgeom);
    }
    *whichcoord=mycoords;
    *whichvel=VEL3;

    FTYPE podpierd=-(((*ptrptrgeom)->gcon[GIND(0,0)])-2.*RADNT_ELL*((*ptrptrgeom)->gcon[GIND(0,3)])+RADNT_ELL*RADNT_ELL*((*ptrptrgeom)->gcon[GIND(3,3)]));
    FTYPE ut=-1./sqrt(podpierd);
    if(!isfinite(ut)) ut=-1.0; // so skips donut, but doesn't give false in condition below (i.e. condition controlled by only podpierd)

    ut/=RADNT_UTPOT; //rescales rin
    // below assumes no torus in unbound region
    if(ut<-1 || podpierd<0. || r<3. || RADNT_DONUTTYPE==NODONUT || RADNT_INFLOWING){
      // allow torus to be unbound with radiation -- not a problem
      //  if(podpierd<0. || r<3. || RADNT_DONUTTYPE==NODONUT || RADNT_INFLOWING)
      //return -1; //outside donut
      usingback=1;
    }

    FTYPE h=-1./ut;
    FTYPE eps=(h-1.)/gamtorus;
    if(usingback==0){
      // from P=K rho^gamma
      rho=pow(eps*(gamtorus-1.)/RADNT_KKK,1./(gamtorus-1.));
      uint=rho*eps;
      pt = uint * (gamtorus-1.0); // torus pressure
      uphi=-RADNT_ELL*ut; // i.e. constant angular momentum torus
      uT=((*ptrptrgeom)->gcon[GIND(0,0)])*ut+((*ptrptrgeom)->gcon[GIND(0,3)])*uphi;
      uPhi=((*ptrptrgeom)->gcon[GIND(3,3)])*uphi+((*ptrptrgeom)->gcon[GIND(0,3)])*ut;
      Vphi=uPhi/uT;
      Vr=0.;
    }



    // see if should still use backup non-torus values
    if(rho<ppback[RHO]){
      usingback=1;
    }
    if(DOWALDDEN){
      usingback=1;
    }

    


    //    dualfprintf(fail_file,"rhodonut1=%g eps=%g uint=%g usingback=%d\n",rho,eps,uint,usingback);

  }
  else if(RADNT_DONUTTYPE==DONUTOHSUGA){

    int mycoords=BLCOORDS;
    if(*whichcoord!=mycoords){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_X(getprim,mycoords,i,j,k,NOWHERE,X,*ptrptrgeom);
    }
    *whichcoord=mycoords;
    *whichvel=VEL3;

    FTYPE rs=2.0; // Schwarzschild radius
    FTYPE r0=RADNT_DONUTRADPMAX;  // pressure maximum radius // FREEPAR
    FTYPE l0=pow(r0,1.5)/(r0-rs); // angular momentum at pressure maximum [type in S3.3 in Ohsuga in text for l0]
    //    FTYPE aa=0.46; // FREEPAR 
    FTYPE aa=RADNT_LPOW; // FREEPAR 
    FTYPE ll=l0*pow(Rcyl/r0,aa); // l ~ r^aa (i.e. non-constant angular momentum)
    FTYPE rhoc=RADNT_RHODONUT; // density at center of torus
    FTYPE nn=1.0/(gamtorus-1.0);  //3.0; // -> 1 + 1/n = gamma -> n = 1/(gamma-1) : gamma=4/3 -> n=3
    // so p = P = K \rho^((n+1)/n) = K \rho^(1 + 1/n) = K \rho^\gamma  [typo in S3.3 in Ohsuga in text part]
    // Neutron stars are well modeled by polytropes with index about in the range between n=0.5 and n=1.
    // A polytrope with index n=1.5  is a good model for fully convective star cores (like those of red giants), brown dwarfs, giant gaseous planets (like Jupiter), or even for rocky planets.
    // Main sequence stars like our Sun and relativistic degenerate cores like those of white dwarfs are usually modeled by a polytrope with index n=3, corresponding to the Eddington standard model of stellar structure.
    //A polytrope with index n=5  has an infinite radius. It corresponds to the simplest plausible model of a self-consistent stellar system, first studied by A. Schuster in 1883.
    //A polytrope with index n=\infty corresponds to what is called isothermal sphere, that is an isothermal self-gravitating sphere of gas, whose structure is identical with the structure of a collisionless system of stars like a globular cluster.
    //Note that the higher the polytropic index, the more condensed at the centre is the density distribution.
    //    FTYPE phi = -(1.0+(*ptrptrgeom)->gcov[GIND(0,0)]); // rough potential, accurate for a=0
    FTYPE phi = -1.0/(r-rs);
    FTYPE phieq = -1.0/(Rcyl-rs);
    FTYPE phir0 = -1.0/(r0-rs); // estimate, good enough for that large radius
    FTYPE phieff=phi + 0.5*pow(ll/Rcyl,2.0)/(1.0-aa);
    FTYPE phieffr0=phir0 + 0.5*pow(l0/r0,2.0)/(1.0-aa);
    FTYPE phieffeq=phieq + 0.5*pow(ll/Rcyl,2.0)/(1.0-aa);
    FTYPE rhot,rhotarg,rhoteq,rhotargeq;
    if(0){
      // Ohsuga & Mineshige (2011):
      FTYPE epsilon0=1.45E-3; // FREEPAR // shifts how hot the flow is
      rhotarg=1.0 - (1.0/(gamtorus*epsilon0*epsilon0*phir0*phir0))*(phieff-phieffr0)/(nn+1.0);
      rhotargeq=1.0 - (1.0/(gamtorus*epsilon0*epsilon0*phir0*phir0))*(phieffeq-phieffr0)/(nn+1.0);
      rhot = rhoc*pow(rhotarg,nn); // [typo in paper Eq18]
      rhoteq = rhoc*pow(rhotargeq,nn); // [typo in paper Eq18]
      pt = rhoc*gamtorus*epsilon0*epsilon0*phir0*phir0*pow(rhot/rhoc,1.0+1.0/nn);
    }
    else{
      // Kato et al. (2004):
      //FTYPE vs0=5.6E-3; // FREEPAR // Table1
      //      FTYPE vs0=1E-1;
      FTYPE HoR0=RADNT_HOVERR;
      FTYPE vphi0=l0/r0;
      FTYPE vs0=HoR0*vphi0;
      rhotarg=1.0 - (gamtorus/(vs0*vs0))*(phieff-phieffr0)/(nn+1.0);
      rhotargeq=1.0 - (gamtorus/(vs0*vs0))*(phieffeq-phieffr0)/(nn+1.0);
      rhot = rhoc*pow(rhotarg,nn);
      rhoteq = rhoc*pow(rhotargeq,nn);
      pt = rhoc*vs0*vs0/(gamtorus)*pow(rhot/rhoc,1.0+1.0/nn);
      FTYPE Eth0=vs0*vs0/(gamtorus*fabs(phir0));
    }

    if(r<3.0 || rhotarg<0.0 || phieff<phieffr0 || rhot<0.0 || rhot>rhoc || phieffeq>0.0 || rhoteq>rhoc || rhoteq<0.0 ){
      usingback=1;
    }

    rho=rhot; // torus density
    uint=pt/(gamtorus-1.0); // torus internal energy density assuming ideal gas with gamtorus

    // see if should still use backup non-torus values
    if(rho<ppback[RHO]){
      usingback=1;
    }

    if(DOWALDDEN){
      usingback=1;
    }

    Vphi=ll/(r*r);
    Vr=0.;

    dualfprintf(fail_file,"donutohsuga: r(%d)=%g th(%d)=%g Rcyl=%g : l0=%g ll=%g phi=%g phieff=%g phir0=%g phieffr0=%g rhot=%g rhotarg=%g pt=%g Vphi=%g\n",r,(*ptrptrgeom)->i,th,(*ptrptrgeom)->j,Rcyl,l0,ll,phi,phieff,phir0,phieffr0,rhot,rhotarg,pt,Vphi);

  }



  //////////
  // see if should still use backup non-torus values
  // could use existing backup, or change.
  if(usingback==1){

    //use existing backup
    rho=ppback[RHO];
    uint=ppback[UU];

    // fix total pressure
    pt = pressure_rho0_u_simple(i,j,k,CENT,rho,uint);
    //      if(EOMRADTYPE!=EOMRADNONE) pt += pp[URAD0]*(4.0/3.0-1.0); // don't include if later say optically thin


    // get KSCOORDS ZAMO since original ppback[] might have been (e.g.) BLCOORDS/VEL3
    // overwrite ppback with MCOORD ZAMO
    *whichvel=VELREL4;
    *whichcoord=MCOORD;
    int getprim=0;
    gset(getprim,*whichcoord,i,j,k,*ptrptrgeom);
    set_zamo_velocity(*whichvel,*ptrptrgeom,ppback);

    // get velocities from new backup
    Vr=ppback[U1];
    Vh=ppback[U2];
    Vphi=ppback[U3];

    // use existing backup
    if(URAD0>=0){
      // still fluid-frame
      E=ppback[URAD0];
      Fx=ppback[URAD1];
      Fy=ppback[URAD2];
      Fz=ppback[URAD3];
    }

  }
  else{

    if(URAD0>=0){
      // not setting radiation here, so use backup for now
      E=pp[URAD0];
      Fx=pp[URAD1];
      Fy=pp[URAD2];
      Fz=pp[URAD3];
    }
  }

  ////
  //
  // assign torus density-velocity values for return
  //
  /////
  pp[RHO]=rho;
  pp[UU]=uint;
  *ptptr=pt;
  pp[U1]=Vr;
  pp[U2]=Vh;
  pp[U3]=Vphi;
  if(URAD0>=0){
    pp[URAD0]=E;
    pp[URAD1]=Fx;
    pp[URAD2]=Fy;
    pp[URAD3]=Fz;
  }

  return(usingback);
}




// return pp (if not already set) for gas primitives as well as assumed total pressure consistent with original simulation's pressure and ideal gas constant gamtorus.
// use instead of donut_analytical_solution() everywhere function called (currently just 1 location)
static int process_solution(int *whichvel, int *whichcoord, int opticallythick, FTYPE *pp,FTYPE *X, FTYPE *V,struct of_geom **ptrptrgeom, FTYPE *ptptr)
{

  FTYPE gamtorus=4.0/3.0; // whatever you had in original restart file.
  // total torus pressure that will be distributed among gas and radiation pressure
  *ptptr=(gamtorus-1.0)*pp[UU];

  int usingback=0; // never use backup, always a solution

  return(usingback);
}




// take global primitive data (say read in by restart_read()) and create radiation component.
// to be called externally after restart_init()
int process_restart_toget_radiation(void)
{
  int i,j,k;

  DUMPGENLOOP{ // could make this OpenMP'ed, but just done once.  Needs to be no more or less than what restart_init() and so restart_read() and so dumpgen() and so DUMPGENLOOP uses.
    int whichvel=WHICHVEL, whichcoord=PRIMECOORDS; // for the below to make sense and work (without further transformation calls), this should be WHICHVEL and PRIMECOORDS always.
    int loc=CENT;
    FTYPE X[NDIM],V[NDIM];
    bl_coord_ijk_2(i,j,k,loc,X,V);
    struct of_geom geomrealdontuse;
    struct of_geom *ptrgeomrad=&geomrealdontuse;
    get_geometry(i,j,k,loc,ptrgeomrad);
    int thick=1; // see how used
    FTYPE *pr=&GLOBALMACP0A1(pglobal,i,j,k,0);
    int returndonut=get_full_rtsolution(&whichvel,&whichcoord,thick, pr,X,V,&ptrgeomrad);
    FTYPE pradffortho[NPR];
    if(PRAD0>=0){
      // donut returns fluid frame orthonormal values for radiation in pr
      pradffortho[PRAD0]=pr[PRAD0];
      pradffortho[PRAD1]=pr[PRAD1];
      pradffortho[PRAD2]=pr[PRAD2];
      pradffortho[PRAD3]=pr[PRAD3];
    }
    prad_fforlab(&whichvel, &whichcoord, FF2LAB, i,j,k,loc,ptrgeomrad, pradffortho, pr, pr);
    // now all pr is PRIMECOORDS, WHICHVEL in lab-frame.
  }

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

static FTYPE taper_func2(FTYPE R,FTYPE rin, FTYPE rpow)
{
  if(R < rin)
    return(0.) ;
  else
    return(pow( 1.-pow(0.95*rin/R,0.5) , 1.)) ; // was  3 and 2 for powers before 11/10/2012 MAVARA
}


int set_fieldtype(void)
{
  int FIELDTYPE;

  if(WHICHPROBLEM==RADDONUT){

    if(RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3){
      //      FIELDTYPE=DISK2FIELD;
      FIELDTYPE=FIELDJONMAD; // SUPERMADNEW

    }
    else if(RADNT_DONUTTYPE==DONUTOLEK){
      //FIELDTYPE=VERTFIELD; // DISK2VERT//DISK2FIELD
      FIELDTYPE=DISK2FIELD;

      //FIELDTYPE=OLEKFIELD;
      //FIELDTYPE=FIELDJONMAD;

      if(DOWALDDEN==1){ // nothing to do with densities, just for simplicity
        FIELDTYPE=FIELDWALD; // WALD field
      }
      else if(DOWALDDEN==2) FIELDTYPE=MONOPOLE;


    }
    else if(RADNT_DONUTTYPE==DONUTOHSUGA){
      FIELDTYPE=OHSUGAFIELD;
    }
    else{
      //FIELDTYPE=VERTFIELD;
      FIELDTYPE=DISK2FIELD; // default
      //FIELDTYPE=DISK1FIELD;
    }
  }
  else if(WHICHPROBLEM==RADBONDI){
    //FIELDTYPE=MONOPOLAR; // for mag bondi
    FIELDTYPE=NOFIELD;
  }
  else{
    FIELDTYPE=NOFIELD;
  }


  return(FIELDTYPE);
}



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
  SFTYPE rho_av, p_av,u_av,q;
  FTYPE r,th,ph;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);

  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();



  FTYPE FRACAPHICUT;

  
  //#define FRACAPHICUT 0.1
  if(FIELDTYPE==DISK2FIELD){
    FRACAPHICUT=0.2; // for weak field
    //    FRACAPHICUT=0.001; // for disk-filling field that is more MAD like


    if(WHICHPROBLEM==RADDONUT && (RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3)){
      FRACAPHICUT=1E-6;
    }

  }
  



  vpot=0.0;

  FTYPE *pr = &MACP0A1(prim,i,j,k,0);

#if(1)
  int ibound=0,jbound=0,kbound=0;
  if(i==-N1BND || i==N1-1+N1BND) ibound=1 && N1NOT1;
  if(j==-N2BND || j==N2-1+N2BND) jbound=1 && N2NOT1;
  if(k==-N3BND || k==N3-1+N3BND) kbound=1 && N3NOT1;

  int iout=0,jout=0,kout=0;
  if(i<-N1BND || i>N1-1+N1BND) iout=1 && N1NOT1;
  if(j<-N2BND || j>N2-1+N2BND) jout=1 && N2NOT1;
  if(k<-N3BND || k>N3-1+N3BND) kout=1 && N3NOT1;
#endif

  // since init_vpot() is called for all i,j,k, can't use
  // non-existence values, so limit averaging:
  if(iout||jout||kout){
    rho_av=p_av=0.0;
  }
  else if(ibound||jbound||kbound){
    //    if(ibound||jbound||kbound){
    //  if(i==-N1BND && j==-N2BND){
    rho_av = pr[RHO];
    p_av = pressure_rho0_u_simple(i,j,k,CENT,pr[RHO],pr[UU]);
    if(EOMRADTYPE!=EOMRADNONE) p_av += pr[URAD0]*(4.0/3.0-1.0);
  }
  /* else if(i==-N1BND){ */
  /*   rho_av = AVGN_2(prim,i,j,k,RHO); */
  /*   u_av = AVGN_2(prim,i,j,k,UU); // simple cheat to avoid defining new AVGN macros */
  /*   p_av = pressure_rho0_u_simple(i,j,k,loc,rho_av,u_av); */
  /*   if(EOMRADTYPE!=EOMRADNONE) p_av += AVGN_2(prim,i,j,k,URAD0)*(4.0/3.0-1.0); */
  /* } */
  /* else if(j==-N2BND){ */
  /*   rho_av = AVGN_1(prim,i,j,k,RHO); */
  /*   u_av = AVGN_1(prim,i,j,k,UU); */
  /*   p_av = pressure_rho0_u_simple(i,j,k,loc,rho_av,u_av); */
  /*   if(EOMRADTYPE!=EOMRADNONE) p_av += AVGN_1(prim,i,j,k,URAD0)*(4.0/3.0-1.0); */
  /* } */
  else{ // normal cells
    if(loc==CORN3){
      rho_av = AVGN_for3(prim,i,j,k,RHO);
      u_av = AVGN_for3(prim,i,j,k,UU);
      p_av = pressure_rho0_u_simple(i,j,k,loc,rho_av,u_av);
      if(EOMRADTYPE!=EOMRADNONE) p_av += AVGN_for3(prim,i,j,k,URAD0)*(4.0/3.0-1.0);
    }
    else if(loc==CORN2){
      rho_av = AVGN_for2(prim,i,j,k,RHO);
      u_av = AVGN_for2(prim,i,j,k,UU);
      p_av = pressure_rho0_u_simple(i,j,k,loc,rho_av,u_av);
      if(EOMRADTYPE!=EOMRADNONE) p_av += AVGN_for2(prim,i,j,k,URAD0)*(4.0/3.0-1.0);
    }
    else if(loc==CORN1){
      rho_av = AVGN_for1(prim,i,j,k,RHO);
      u_av = AVGN_for1(prim,i,j,k,UU);
      p_av = pressure_rho0_u_simple(i,j,k,loc,rho_av,u_av);
      if(EOMRADTYPE!=EOMRADNONE) p_av += AVGN_for1(prim,i,j,k,URAD0)*(4.0/3.0-1.0);
    }
    else{
      dualfprintf(fail_file,"No such setup for loc=%d\n",loc);
      myexit(34782985);
    }
  }




  if(FIELDTYPE==TOROIDALFIELD){
    
    if(l==2){// A_\theta (MCOORD)
      
      r=V[1];
      th=V[2];

      //      q = r*r*r*fabs(sin(th)) * 1.0 ; // constant B^\phi
      q = r*r*r;

      q=q/(r); // makes more uniform in radius

      q = q*(p_av / ptotmax - FRACAPHICUT); // weight by pressure
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
  // FTYPE rpow=1.0; // previous SUPERMAD, now just use 3/4
  rpow=5.0/4.0; // 1.0 so libetatot constant vs. radius


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


    if(FIELDTYPE==MONOPOLAR){
      vpot += (1.0-cos(th));
    }

    if(FIELDTYPE==OLEKFIELD){
      vpot += MAX(pow(r,4.0)*pow(rho_av,2.0)*1E40-0.02,0.0)*pow(sin(th),4.0);
    }

#define JONMADHPOW (4.0)
#define JONMADR0 (0.0)
#define JONMADRIN (2.0)
#define JONMADROUT (NTROUT)

    FTYPE jonmadhpow;
    FTYPE interp;
    interp = (r-Rhor)/(r+1E-10);
    if(interp<0.0) interp=0.0;
    if(interp>1.0) interp=1.0;

    jonmadhpow = JONMADHPOW*interp + 0.0*( 1.0-interp );
    
    FTYPE thetafactout,thetafactin,thetafact1;
    thetafactin=1.0-cos(th);
    if(th>0.5*M_PI) thetafactin=1.0-cos(M_PI-th);
    thetafactout=pow(sin(th),1.0+jonmadhpow);
    thetafact1 = thetafactout*interp + thetafactin*( 1.0-interp );

    FTYPE interp2;
    interp2 = (r-JONMADROUT)/(r+1E-10);
    if(interp2<0.0) interp2=0.0;
    if(interp2>1.0) interp2=1.0;

    FTYPE thetafactout2,thetafactin2,thetafact2;
    thetafactin2=thetafact1;
    thetafactout2=1.0-cos(th);
    if(th>0.5*M_PI) thetafactout2=1.0-cos(M_PI-th);
    thetafact2 = thetafactout2*interp2 + thetafactin2*( 1.0-interp2 );

    FTYPE thetafact;
    thetafact=thetafact2;

    FTYPE rfact;
    rfact=MAX(pow(r-JONMADR0,rpow)*1E40-0.02,0.0);
    if(r>=JONMADROUT){
      rfact=MAX(pow(JONMADROUT-JONMADR0,rpow)*1E40-0.02,0.0);
    }
    if(r<=JONMADRIN){
      rfact=MAX(pow(JONMADRIN-JONMADR0,rpow)*1E40-0.02,0.0);
    }

    if(FIELDTYPE==FIELDJONMAD){
      
      vpot += rfact*thetafact;
        
      if(V[2]<1E-5 || V[2]>M_PI-1E-5){
        vpot=0;
      }


    }





    /* field-in-disk version */
    if(FIELDTYPE==DISK1FIELD || FIELDTYPE==DISK1VERT){
      q = rho_av / rhomax - 0.2;
      if(r<rin) q=0.0;
      if (q > 0.)      vpot += q;
    }

    if(FIELDTYPE==OHSUGAFIELD){
      q = rho_av / rhomax - 0.2;
      if(r<rin) q=0.0;
      if (q > 0.)      vpot += q;
    }


    if(FIELDTYPE==DISK2FIELD || FIELDTYPE==DISK2VERT){
      // average of density that lives on CORN3


      //#define FRACAPHICUT 0.1
      //#define FRACAPHICUT 0.1

      //      q = (rho_av / rhomax - FRACAPHICUT);
      //      q = (p_av / ptotmax - FRACAPHICUT); // was used for rada0.94, etc. models.
      if(WHICHPROBLEM==RADDONUT && (RADNT_DONUTTYPE==DONUTTHINDISK || RADNT_DONUTTYPE==DONUTTHINDISK2||RADNT_DONUTTYPE==DONUTTHINDISK3)){
        q = (p_av);
        //        dualfprintf(fail_file,"qorig=%g\n",q);
        if(rho_av/(rhomax*pow(r,thindiskrhopow)) - FRACAPHICUT<0) q=0;
        //        dualfprintf(fail_file,"qnew=%g : %g %g\n",q,rho_av/rhomax,FRACAPHICUT);
      }
      else{
        q = 1E-30*(p_av / ptotmax*1E30 - FRACAPHICUT);
        //        dualfprintf(fail_file,"qorig=%g\n",q);
        if(rho_av/rhomax-FRACAPHICUT<0) q=0;
        //        dualfprintf(fail_file,"qnew=%g : %g %g\n",q,rho_av/rhomax,FRACAPHICUT);
      }

      //#define QPOWER 0.5
#define QPOWER (1.0)

#define POWERNU (2.0) // 2.5 for toroidal field SUPERMAD paper
      //#define POWERNU (4.0)

      //      if (q > 0.)      vpot += q*q*pow(r*fabs(sin(th)),POWERNU);
      FTYPE fact1,fact2,SSS,TTT;
      fact1=pow(fabs(q),QPOWER)*pow(r*fabs(sin(th)),POWERNU);
      // fact1=pow(fabs(q),QPOWER)*pow(r,POWERNU); // for SUPERMAD paper
      if(r<rin) fact1=0.0;
      SSS=rin*0.5;
      TTT=0.28;
      fact2=sin(log(r/SSS)/TTT);
      fact2=1.0; // forces to avoid flipping.

      if (q > 0.)      vpot += fact1*fact2;
      //      if (q > 0.)      vpot += q*q;
      if(PRODUCTION==0) dualfprintf(fail_file,"ijk=%d %d %d : ptotmax=%g p_av=%g q=%g %g %g rin=%g vpot=%g : rho_av=%g rhomax=%g\n",i,j,k,ptotmax,p_av,q,fact1,fact2,rin,vpot,rho_av,rhomax);
    }



  }

  if(FIELDTYPE==SPLITMONOPOLE || FIELDTYPE==MONOPOLE || FIELDTYPE==FIELDWALD){

    r=V[1];
    th=V[2];
    ph=V[3];


    if(FIELDTYPE==SPLITMONOPOLE || FIELDTYPE==MONOPOLE) return(0); // otherwise setup poloidal components using vector potential
    else if(FIELDTYPE==FIELDWALD){

      FTYPE mcov[NDIM],mcon[NDIM],kcov[NDIM],kcon[NDIM];
      

      mcon[TT]=0;
      mcon[RR]=0;
      mcon[TH]=0;
      mcon[PH]=1.0;
      
      kcon[TT]=1.0;
      kcon[RR]=0;
      kcon[TH]=0;
      kcon[PH]=0;

      // get metric grid geometry for these ICs
      int getprim=0;
      struct of_geom geomrealdontuse;
      struct of_geom *ptrgeomreal=&geomrealdontuse;
      *whichcoord = MCOORD;
      gset_genloc(getprim,*whichcoord,i,j,k,loc,ptrgeomreal);

      lower_vec(mcon,ptrgeomreal,mcov);
      lower_vec(kcon,ptrgeomreal,kcov);
      

      // below so field is b^2/rho=BSQORHOWALD at horizon
      FTYPE Rhorlocal=rhor_calc(0);
      B0WALD=sqrt(BSQORHOWALD*RADNT_RHOATMMIN*pow(Rhorlocal/RADNT_ROUT,-1.5));
      TILTWALD=THETAROT;
      aforwald=a;
      
      if(1||fabs(TILTWALD-0.0)<1E-10){ // 1|| for now since issue with uu0 at maximum near outer radius when tilt=90deg-1E-5
        if(l==WALDWHICHACOV || WALDWHICHACOV==-1){
          vpot += -0.5*B0WALD*(mcov[l]+2.0*aforwald*kcov[l]);
        }


      }
      else{ // need all A_i's for tilted field

        // STEP1:
        // r,\theta,\phi are such that BH spin stays in z-hat and field is tilted.  So we need to provide r,\theta,\phi in coordinates wher z-hat is along spin axis.

        // Original V at this point is with tilted BH spin, but solution below is for spin along z-axis, so find Vmetric[corresponding to when BH spin is along z-axis]
        FTYPE Vmetric[NDIM];
        rotate_VtoVmetric(MCOORD,TILTWALD,V,Vmetric);

        FTYPE rV=Vmetric[1];
        FTYPE thV=Vmetric[2];
        FTYPE phV=Vmetric[3];

        // STEP2: Now get solution in setup where BH spin along z-axis always and it is the field that is rotated
        // B0y=0
        FTYPE B0z=1.0; // B0z
        FTYPE B0x=0.0; // B0x

        FTYPE delta,sigma,rp,rm,psi;

        delta= rV*rV - 2.0*MBH*rV + a*a;
        sigma = rV*rV + a*a*cos(thV)*cos(thV);
        rp = MBH + sqrt(MBH*MBH - a*a);
        rm = MBH - sqrt(MBH*MBH - a*a);
        
        psi = phV ; // + a/(rp-rm)*log((rV-rp)/(rV-rm)); // ph originally is phi[BL] and psi is phi[KS].

        FTYPE Acovblnonrot[NDIM];

        Acovblnonrot[0] = -1.*a*B0z + (a*B0z*MBH*rV*(1. + Power(Cos(thV),2)))/sigma + 
          (a*B0x*MBH*Cos(thV)*(rV*Cos(psi) - 1.*a*Sin(psi))*Sin(thV))/sigma;

        Acovblnonrot[1] = -1.*B0x*(-1.*MBH + rV)*Cos(thV)*Sin(psi)*Sin(thV);

        Acovblnonrot[2] = -1.*B0x*(Power(rV,2)*Power(Cos(thV),2) + Power(a,2)*Cos(2.*thV) - 
                                   1.*MBH*rV*Cos(2.*thV))*Sin(psi) - 
          1.*a*B0x*Cos(psi)*(MBH*Power(Cos(thV),2) + rV*Power(Sin(thV),2));
        
        
        Acovblnonrot[3] = -1.*B0x*Cos(thV)*(delta*Cos(psi) + 
                                            (MBH*(Power(a,2) + Power(rV,2))*(rV*Cos(psi) - 1.*a*Sin(psi)))/
                                            sigma)*Sin(thV) + B0z*(0.5*(Power(a,2) + Power(rV,2)) - 
                                                                   (1.*Power(a,2)*MBH*rV*(1. + Power(Cos(thV),2)))/sigma)*
          Power(Sin(thV),2);

        // STEP3: now take Acovblnonrot -> Acovksnonrot
        // NOTE: bl2ks transformation assumes BH spin is pointing z-hat, just like Acovblnonrot[] assumes, so have to apply this here before rotation vector.
        FTYPE Acovksnonrot[NDIM];
        int jj;
        DLOOPA(jj) Acovksnonrot[jj]=Acovblnonrot[jj];
        bltoks_ucov(i,j,k,loc,Acovksnonrot);

        // STEP4: now take Acovksnonrot -> Acovks.
        // That is, while coordinate position was correctly mapped in STEP1, coordinate vector still needs to be rotated
        FTYPE Acovks[NDIM];
        DLOOPA(jj) Acovks[jj]=Acovksnonrot[jj];
        transVmetrictoV_ucov(TILTWALD,Vmetric,Acovks);

        // now pluck out only the A_l one wanted
        vpot += Acovks[l];
 
      }
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






// setup primitive field (must use whichcoord inputted)
// sets CENT value of field primitive
static int fieldprim(int whichmethod, int whichinversion, int *whichvel, int*whichcoord, int ii, int jj, int kk, FTYPE *pr)
{



  /////
  //
  // Get FIELDTYPE
  //
  /////
  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();


  /////
  //
  // get X and V
  //
  ////
  FTYPE X[NDIM],V[NDIM];
  coord(ii, jj, kk, CENT, X);
  bl_coord(X, V);
  FTYPE r=V[1];
  FTYPE th=V[2];
  FTYPE ph=V[3];


  /////
  //
  // Get ptrgeom
  //
  /////
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE realX[NDIM];
  if(*whichcoord==PRIMECOORDS){
    get_geometry(ii,jj,kk,CENT,ptrgeom);
    int j;
    DLOOPA(j) realX[j]=X[j];
    //    dualfprintf(fail_file,"using PRIMECOORDS\n");
  }
  else{
    // get metric grid geometry for these ICs
    int getprim=0;
    gset(getprim,*whichcoord,ii,jj,kk,ptrgeom);
    int j;
    DLOOPA(j) realX[j]=V[j];
    //    dualfprintf(fail_file,"using WHICHCOORD=%d\n",*whichcoord);
  }



  /////
  //
  // set field
  //
  ////
  FTYPE Rhorlocal=rhor_calc(0);
  B0WALD=sqrt(BSQORHOWALD*RADNT_RHOATMMIN*pow(Rhorlocal/RADNT_ROUT,-1.5));
  //
  if(FIELDTYPE==SPLITMONOPOLE){
    if(th<M_PI*0.5)  pr[B1]=B0WALD/ptrgeom->gdet;
    else pr[B1]=-B0WALD/ptrgeom->gdet;
    pr[B2]=pr[B3]=0;
  }
  else if(FIELDTYPE==MONOPOLE){
    // Ruben's talk says they set $\dF^{tr} = C\sin{\theta}/\detg$.
    pr[B1]=B0WALD/ptrgeom->gdet;
    pr[B2]=pr[B3]=0;
  }
  else if(FIELDTYPE==FIELDWALD){

    int doit=ii==N1-1 && jj==N2/4-1 && kk==0;
    doit=0; // no debug normally


    FTYPE Fcov[NDIM][NDIM];
    FTYPE Mcon[NDIM][NDIM];
    FTYPE etacov[NDIM],etacon[NDIM];
    FTYPE Ecov[NDIM],Econ[NDIM],Bcov[NDIM],Bcon[NDIM];
    FTYPE alpha;
    FTYPE Jcon[NDIM];
    int j,k;

    void Fcov_numerical(int whichcoord, FTYPE *X, FTYPE (*Fcov)[NDIM]);
    void Jcon_numerical(int whichcoord, FTYPE *X, FTYPE *Jcon);
    extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
    extern void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);
    extern void raise_A(FTYPE (*Acov)[NDIM], struct of_geom *geom, FTYPE (*Acon)[NDIM]);//extern void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM]);
    extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
    extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);

    
    FTYPE prold[NPR];
    int pliter,pl;
    PLOOP(pliter,pl) prold[pl]=pr[pl];
    
    // check that J^\mu=0
    Jcon_numerical(*whichcoord, X, Jcon);

    FTYPE Fcon[NDIM][NDIM];
    FTYPE Fud[NDIM][NDIM];
    FTYPE Fdu[NDIM][NDIM];
    FTYPE Mud[NDIM][NDIM];
    FTYPE Mdu[NDIM][NDIM];
    FTYPE Mcov[NDIM][NDIM];
    if(whichmethod==0){
      //    first get F_{\mu\nu}
      Fcov_numerical(*whichcoord, X, Fcov);
      
      // get *F^{\mu\nu}
      MtoF(3,Fcov,ptrgeom,Mcon);
    }
    else{
      //      int j,k;
      DLOOP(j,k) Mcon[j][k]=0.0;
      SLOOPA(j) Mcon[j][0] = Bcon[j] = prold[B1+j-1]; // use original field from staggered field to centered from A_i
      SLOOPA(j) Mcon[0][j] = -Mcon[j][0];

      if(1){ // set *F_{\phi t}=0 so no angular momentum flux
        FTYPE myBd3=0.0;
        FTYPE Bu1=Bcon[1];
        FTYPE Bu2=Bcon[2];
        FTYPE Ed1=Mcon[2][3]; // *F^{23} = E_1 = F_{10}/detg
        FTYPE Ed2=Mcon[3][1]; // *F^{31} = E_2 = F_{20}/detg
        FTYPE Ed3=Mcon[1][2]; // *F^{12} = E_3 = F_{30}/detg
        FTYPE gv00=ptrgeom->gcov[GIND(0,0)];
        FTYPE gv01=ptrgeom->gcov[GIND(0,1)];
        FTYPE gv02=ptrgeom->gcov[GIND(0,2)];
        FTYPE gv03=ptrgeom->gcov[GIND(0,3)];
        FTYPE gv13=ptrgeom->gcov[GIND(1,3)];
        FTYPE gv23=ptrgeom->gcov[GIND(2,3)];
        FTYPE gv33=ptrgeom->gcov[GIND(3,3)];
        
        FTYPE denom=(-Power(gv03,2) + gv00*gv33);
        //        if(fabs(denom)<0.2 || fabs(r-Rhor)/Rhor<0.3){
        if(fabs(r-Rhor)/Rhor<0.3){
          Bcon[3] = Mcon[3][0] = 0.0;
        }
        else{
          Bcon[3] = Mcon[3][0] =(Bu1*gv01*gv03 + Bu2*gv02*gv03 - Bu1*gv00*gv13 - Ed3*gv02*gv13 + 
                                 Ed2*gv03*gv13 - Bu2*gv00*gv23 + Ed3*gv01*gv23 - Ed1*gv03*gv23 - 
                                 Ed2*gv01*gv33 + Ed1*gv02*gv33 + myBd3)/denom;
          //        Bcon[3] = Mcon[3][0] = (-(Bu1*gv01*gv03) - Bu2*gv02*gv03 + Bu1*gv00*gv13 + Bu2*gv00*gv23 + myBd3)/(pow(gv03,2) - gv00*gv33);
        }


        Mcon[0][3] = -Mcon[3][0];
      }



      MtoF(0,Mcon,ptrgeom,Fcov);
    }

    raise_A(Fcov,ptrgeom,Fcon);
    indices_2221(Fcon,Fud,ptrgeom);
    indices_2212(Fcon,Fdu,ptrgeom);

    indices_2221(Mcon,Mud,ptrgeom);
    indices_2212(Mcon,Mdu,ptrgeom);
    lower_A(Mcon,ptrgeom,Mcov);


    //    dualfprintf(fail_file,"MtoF\n");
    //    DLOOP dualfprintf(fail_file,"Mcon[%d][%d]=%21.15g\n",j,k,Mcon[j][k]);



    // T^\mu_\nu
    FTYPE Tud[NDIM][NDIM];
    int ll,mm;
    FTYPE Fsq=0.0;
    DLOOP(j,k) Fsq += Fcon[j][k]*Fcov[j][k];
    DLOOP(j,k) Tud[j][k] = 0.0;
    DLOOP(j,k){
      DLOOPA(ll) Tud[j][k] += Fud[j][ll]*Fdu[k][ll];
      Tud[j][k] += - 0.25*delta(j,k)*Fsq;
    }
    if(doit) DLOOP(j,k) dualfprintf(fail_file,"Tud[%d][%d]=%g\n",j,k,Tud[j][k]);


    // lapse
    // define \eta_\alpha
    // assume always has 0 value for space components
    alpha = 1./sqrt(-ptrgeom->gcon[GIND(0,0)]);
      
    etacov[TT]=-alpha; // any constant will work.
    SLOOPA(j) etacov[j]=0.0; // must be 0

    // shift
    // define \eta^\beta
    raise_vec(etacov,ptrgeom,etacon);
    //    dualfprintf(fail_file,"raise\n");

    //    DLOOPA dualfprintf(fail_file,"etacon[%d]=%21.15g etacov[%d]=%21.15g\n",j,etacon[j],j,etacov[j]);

    // Betacon
    FTYPE Betacon[NDIM],Betacov[NDIM];
    DLOOPA(j) Betacon[j]=0.0;
    DLOOP(j,k) Betacon[j]+=etacov[k]*Mcon[k][j];
    lower_vec(Betacon,ptrgeom,Betacov);

    //    dualfprintf(fail_file,"Fcov\n");
    //DLOOP dualfprintf(fail_file,"Fcov[%d][%d]=%21.15g\n",j,k,Fcov[j][k]);
    //    dualfprintf(fail_file,"%21.15g %21.15g\n",j,k,Fcov[0][3],Fcov[3][0]);

    // then get Eeta^\alpha
    FTYPE Eetacov[NDIM],Eetacon[NDIM];
    DLOOPA(j) Eetacov[j]=0.0;
    DLOOP(j,k) Eetacov[j]+=etacon[k]*Fcov[j][k];
    raise_vec(Eetacov,ptrgeom,Eetacon);
    if(doit) dualfprintf(fail_file,"etacon=%g %g %g %g Eetacov[3]=%21.15g Fcov03=%g\n",etacon[0],etacon[1],etacon[2],etacon[3],Eetacov[3],Fcov[0][3]);
    //DLOOPA dualfprintf(fail_file,"Eetacon[%d]=%2.15g\n",j,Eetacon[j]);

    // T^\mu_\nu from Eeta and Beta
    extern void EBtoT(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE (*T)[NDIM]);
    FTYPE TfromEB[NDIM][NDIM];
    EBtoT(Eetacon,Betacon,ptrgeom,TfromEB);
    if(doit) DLOOP(j,k) dualfprintf(fail_file,"TudfromEB[%d][%d]=%g\n",j,k,TfromEB[j][k]);



    if(whichinversion==0){

      FTYPE U[NPR];
      DLOOPA(j) U[UU+j] = Tud[0][j];
      SLOOPA(j) U[B1+j-1] = Bcon[j] = Mcon[j][0]; // B^i = \dF^{it}
      Bcon[TT]=0.0; // force


      if(doit) dualfprintf(fail_file,"CONS: %g %g %g %g : %g %g %g: F=%g %g\n",U[UU],U[U1],U[U2],U[U3],Bcon[1],Bcon[2],Bcon[3],Fcov[0][3],Fcov[3][0]);
      
      struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
      int eomtypelocal=EOMFFDE;
      //      int eomtypelocal=EOMTYPE;
      struct of_state qdontuse;
      struct of_state *qptr=&qdontuse;
      // assume if needed rest of pr already set
      SLOOPA(j) pr[B1+j-1]=Bcon[j];
      pr[RHO]=pr[UU]=0.0;
      //SLOOPA(j) pr[U1+j-1]=0.0; // only valid if WHICHVEL==VELREL4 // just assume use "old" versions
      if(doit) dualfprintf(fail_file,"BEFORE Utoprimgen()\n");
      if(doit) PLOOP(pliter,pl) dualfprintf(fail_file,"oldpr[%d]=%g\n",pl,pr[pl]);
      Utoprimgen(0,0,0,0,1,&eomtypelocal,CAPTYPEBASIC,0,1,EVOLVEUTOPRIM,UNOTHING,U,qptr,ptrgeom,0,pr,pr,&newtonstats);
      PLOOP(pliter,pl) if(pl>=URAD1 && pl<=URAD3) pr[pl]=prold[pl];
      if(doit) dualfprintf(fail_file,"AFTER Utoprim()\n");
      if(doit) PLOOP(pliter,pl) dualfprintf(fail_file,"newpr[%d]=%g\n",pl,pr[pl]);

      if(EOMTYPE!=EOMFFDE && EOMTYPE!=EOMFFDE2){
        PLOOP(pliter,pl){
          if(pl==RHO || pl==UU || pl>B3 || pl>=URAD1 && pl<=URAD3) pr[pl]=prold[pl];
        }
        
        struct of_state qdontuse4;
        struct of_state *qptr4=&qdontuse4;
        get_state(pr,ptrgeom,qptr4);
        FTYPE Ueomtype[NPR],Ueomtypeabs[NPR];
        primtoflux(UNOTHING,pr,qptr4,TT,ptrgeom,Ueomtype,Ueomtypeabs);
        eomtypelocal=EOMTYPE; // now do EOMTYPE
        Utoprimgen(0,0,0,0,1,&eomtypelocal,CAPTYPEBASIC,0,1,EVOLVEUTOPRIM,UNOTHING,Ueomtype,qptr4,ptrgeom,0,pr,pr,&newtonstats);
        //        PLOOP(pliter,pl) if(pl>=URAD1 && pl<=URAD3) pr[pl]=prold[pl];
        if(doit) PLOOP(pliter,pl) dualfprintf(fail_file,"newprMHD[%d]=%g\n",pl,pr[pl]);
        
      }
 
 

    }
    else if(whichinversion==1){
      

      //    DLOOPA dualfprintf(fail_file,"Econ[%d]=%21.15g Bcon[%d]=%21.15g\n",j,Econ[j],j,Bcon[j]);


      // Use Eeta,Beta to get primitives (v will be in WHICHVEL)
      // ASSUMES FORCE-FREE!
      if(doit) dualfprintf(fail_file,"BEFORE EBtopr()\n");
      if(doit) dualfprintf(fail_file,"%g %g %g %g : %g %g %g %g\n",Eetacon[0],Eetacon[1],Eetacon[2],Eetacon[3],Betacon[0],Betacon[1],Betacon[2],Betacon[3]);
      EBtopr(Eetacon,Betacon,ptrgeom,pr);
      if(doit) dualfprintf(fail_file,"AFTER EBtopr()\n");
      if(doit) PLOOP(pliter,pl) dualfprintf(fail_file,"newpr[%d]=%g\n",pl,pr[pl]);
    }

    if(EOMTYPE==EOMFFDE||EOMTYPE==EOMFFDE2){
      // revert radiation primitives to unmodified values (in whichvel, whichcoord)
      PLOOP(pliter,pl) if(pl>=URAD1 && pl<=URAD3) pr[pl]=prold[pl];
    }
    //    if(EOMTYPE!=EOMFFDE&& EOMTYPE!=EOMFFDE2){
    //      // revert rho and ug and anything beyond field
    //      PLOOP(pliter,pl) if((pl>=RHO && pl<=UU) || (pl>B3 && pl>URAD3)) pr[pl]=prold[pl];      
    //    }

  
    // some checks
    struct of_state qdontuse2;
    struct of_state *qptr2=&qdontuse2;
    get_state(pr,ptrgeom,qptr2);
    bcon_calc(pr,qptr2->ucon,qptr2->ucov,qptr2->bcon);
    FTYPE fdd03;
    fdd03 = ptrgeom->gdet * (qptr2->ucon[1]*qptr2->bcon[2] - qptr2->ucon[2]*qptr2->bcon[1]) ;
    if(doit) dualfprintf(fail_file,"fdd03=%g\n",fdd03);
    
    FTYPE mhdflux[NDIM][NDIM];
    DLOOPA(j)  mhd_calc(pr, j, ptrgeom, qptr2, mhdflux[j], NULL);
    if(doit) dualfprintf(fail_file,"mhdflux=%g %g %g\n",mhdflux[1][0],mhdflux[2][0],mhdflux[3][0]);
    FTYPE mhdfluxma[NDIM][NDIM];
    DLOOPA(j)  mhd_calc_ma(pr, j, ptrgeom, qptr2, mhdfluxma[j], NULL, NULL, NULL);
    if(doit) dualfprintf(fail_file,"mhdfluxma=%g %g %g\n",mhdfluxma[1][0],mhdfluxma[2][0],mhdfluxma[3][0]);
    FTYPE mhdfluxem[NDIM][NDIM];
    DLOOPA(j)  mhd_calc_em(pr, j, ptrgeom, qptr2, mhdfluxem[j], NULL);
    if(doit) dualfprintf(fail_file,"mhdfluxem=%g %g %g\n",mhdfluxem[1][0],mhdfluxem[2][0],mhdfluxem[3][0]);
    
    

    // transform WHICHVEL back to *whichvel (only for plasma, not radiation)
    FTYPE ucontemp[NDIM];
    //    dualfprintf(fail_file,"BEFORE pr2ucon\n");
    int return1=pr2ucon(WHICHVEL,pr,ptrgeom,ucontemp);
    //    dualfprintf(fail_file,"AFTER pr2ucon\n");
    if(failed || return1){
      for(pl=RHO;pl<=U3;pl++) pr[pl]=prold[pl];
      failed=0;
    }
    else{
      //      dualfprintf(fail_file,"BEFORE ucon2pr: %d\n",*whichvel);
      ucon2pr(*whichvel,ucontemp,ptrgeom,pr);
      //      dualfprintf(fail_file,"AFTER ucon2pr\n");
      //      PLOOP(pliter,pl) dualfprintf(fail_file,"from ucon2pr[%d]=%g\n",pl,pr[pl]);

      if(0&&r<Rhor){
        for(pl=RHO;pl<=U3;pl++) pr[pl]=prold[pl] + (pr[pl]-prold[pl])/(Rhor-0.0)*(r-0.0);
      }
      else{
        // keep
      }
    }


    if(EOMTYPE==EOMFFDE || EOMTYPE==EOMFFDE2){
      //      pr[U1]=pr[U2]=pr[U3]=0.0;
      //    dualfprintf(fail_file,"EBtopr\n");
    }
    else{ // just use atmosphere (normally a ZAMO) frame for fluid frame
      //      for(pl=RHO;pl<=U3;pl++) pr[pl]=prold[pl];
    }


    int BOOSTFIELD=0; // for moving BH problem // WALD

    if(BOOSTFIELD){
      // BOOST of field
      FTYPE xx=r*sin(th)*cos(ph);
      FTYPE yy=r*sin(th)*sin(ph);
      FTYPE zz=r*cos(th);
      FTYPE lambdatrans[NDIM][NDIM];
      FTYPE ilambdatrans[NDIM][NDIM];
      // assume time doesn't change or mix with space
      lambdatrans[TT][TT]=1.0;
      SLOOPA(j) lambdatrans[TT][j] = lambdatrans[j][TT] = 0.0;

      ilambdatrans[TT][TT]=1.0;
      SLOOPA(j) ilambdatrans[TT][j] = ilambdatrans[j][TT] = 0.0;

      // rest come from definitions of {x,y,z}(r,\theta,\phi)
      // assumes orthonormal to orhonormal!
      lambdatrans[1][RR] = sin(th)*cos(ph);
      lambdatrans[1][TH] = cos(th)*cos(ph);
      lambdatrans[1][PH] = -sin(ph);

      lambdatrans[2][RR] = sin(th)*sin(ph);
      lambdatrans[2][TH] = cos(th)*sin(ph);
      lambdatrans[2][PH] = cos(ph);

      lambdatrans[3][RR] = cos(th);
      lambdatrans[3][TH] = -sin(th);
      lambdatrans[3][PH] = 0.0;

      // Cart 2 SPC
      ilambdatrans[1][RR] = sin(th)*cos(ph);
      ilambdatrans[1][TH] = cos(th)*cos(ph);
      ilambdatrans[1][PH] = -sin(ph);

      ilambdatrans[2][RR] = sin(th)*sin(ph);
      ilambdatrans[2][TH] = cos(th)*sin(ph);
      ilambdatrans[2][PH] = cos(ph);

      ilambdatrans[3][RR] = cos(th);
      ilambdatrans[3][TH] = -sin(th);
      ilambdatrans[3][PH] = 0.0;

      // quasi-orthonormal
      FTYPE finalvec[NDIM];
      finalvec[TT]=0.0;
      finalvec[RR]=0.3; // x
      finalvec[TH]=0; // y
      finalvec[PH]=0; // z


      // transform from ortho Cart to ortho SPC
      FTYPE tempcomp[NDIM];
      DLOOPA(j) tempcomp[j]=0.0;
      DLOOP(j,k){
        tempcomp[k] += ilambdatrans[j][k]*finalvec[j];
      }
      DLOOPA(j) finalvec[j]=tempcomp[j]; // spc


      // add non-ortho coordinate basis velocity to ortho
      finalvec[TT]=0.0;
      finalvec[RR]=pr[U1] + finalvec[1]/sqrt(ptrgeom->gcov[GIND(1,1)]);
      finalvec[TH]=pr[U2] + finalvec[2]/sqrt(ptrgeom->gcov[GIND(2,2)]);
      finalvec[PH]=pr[U3] + finalvec[3]/sqrt(ptrgeom->gcov[GIND(3,3)]);

      pr[U1] = finalvec[RR];
      pr[U2] = finalvec[TH];
      pr[U3] = finalvec[PH];
    }    


    // stick J^\mu into dump file
#if(0)
    for(k=U1;k<=U1+3;k++){
      pr[k] = Jcon[k-U1];
    }
#endif
    //    dualfprintf(fail_file,"ii=%d jj=%d\n",ii,jj);

#if(0)
    // E.B
    k=U1;
    pr[k]=0;
    DLOOPA(j) pr[k]+=Ecov[j]*Bcon[j];

    // B^2-E^2
    k=UU;
    pr[k]=0;
    DLOOPA(j) pr[k]+=Bcov[j]*Bcon[j];
    //    DLOOPA pr[k]+=Bcon[j];
    DLOOPA(j) pr[k]-=Ecov[j]*Econ[j];

    pr[B1]=Econ[1]*geom.g/etacov[TT];
#endif


#if(0)
    struct of_state q;
    FTYPE faradaytest[NDIM][NDIM];

    // check where faraday changed
    get_state(pr,ptrgeom,&q);

    faraday_calc(0,q.bcon,q.ucon,&geom,faradaytest);
    //    DLOOP dualfprintf(fail_file,"%21.15g  %21.15g\n",faradaytest[j][k],Fcov[j][k]);
    DLOOP{
      if(fabs(faradaytest[j][k]-Fcov[j][k])>1E-10){
        dualfprintf(fail_file,"1 %d %d : %21.15g  %21.15g\n",ii,jj,faradaytest[j][k],Fcov[j][k]);
      }
    }
    if(fabs(faradaytest[0][3])>1E-10) dualfprintf(fail_file,"1 Fcov=%21.15g faraday=%21.15g\n",Fcov[0][3],faradaytest[0][3]);
#endif
 
  }// end if FIELDWALD


  return(0);
}

// compute after A_\mu -> Bstag^i -> Bcent^i to use proper centered field in inversion to get centered velocity for, e.g., Wald problem
int init_postvpot(int i, int j, int k, FTYPE *pr, FTYPE *pstag, FTYPE *ucons){

  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();

  if(FIELDTYPE==SPLITMONOPOLE || FIELDTYPE==MONOPOLE || FIELDTYPE==FIELDWALD){

    int whichvel=WHICHVEL;
    int whichcoord=PRIMECOORDS;
    int whichmethod=1;
    int whichinversion=0;
    fieldprim(whichmethod, whichinversion, &whichvel, &whichcoord, i, j, k, pr);

  }

  return(0);
}


// Also tried:

// For final conservation and noise
// 1) uniform grid near hole for half of grid (/home/jon/coordradial.nb)
// 2) Turned off POLEDEATH and GAMMAPOLEDEATH
// 3) advance_standard_orig() instead.

// For hot MHD part in cylindrical part of jet that gets in way
// 1) Cleaning/cooling off as floor held in wave front.  Still happens.

// For initial energy-momentum flux:
// 1) u^i at r>3
// 2) B_\phi (this function) called after the second set_grid_all() and before copy_prim2panalytic() in initbase.c
// 3) F_{it}=0 before used but after set
// 4) a=0 before first set_grid_all, a=0.8 before second
// 5) Wald A_\mu has zero a term.
// 6) Only use A_\phi
// 7) fieldfrompotential[B3]=0 so only use poloidal Wald.
// 8) Ecov=0
// 9) Econ=0
void shortout_Bd3(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pliter,pl;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FULLLOOP{
    get_geometry(i,j,k,CENT,ptrgeom);
    FTYPE Bcov[NDIM];
    FTYPE *pr = &MACP0A1(prim,i,j,k,0);

    FTYPE Bu1,Bu2,gcon03,gcon13,gcon23,gcon33;
    FTYPE gcov01,gcov02,gcov11,gcov12,gcov21,gcov22,gcov03,gcov13,gcov23;

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
    
    FTYPE myBd3=0.0;

    FTYPE ftemp=(1.0 - gcon03*gcov03 - gcon13*gcov13 - gcon23*gcov23);
    FTYPE igdetnosing=sign(ftemp)/(fabs(ftemp)+SMALL);
    pl=B3; pr[pl] = (myBd3*gcon33 + Bu1*gcon03*gcov01 + Bu2*gcon03*gcov02 + Bu1*gcon13*gcov11 + Bu2*gcon13*gcov12 + Bu1*gcon23*gcov21 + Bu2*gcon23*gcov22)*igdetnosing;

    FTYPE ucon[NDIM];
    pr2ucon(WHICHVEL,pr,ptrgeom,ucon);

    FTYPE X[NDIM],V[NDIM];
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    FTYPE r,th,ph;
    r=V[1];
    th=V[2];
    ph=V[3];

    if(r>3.0) ucon[1]=ucon[2]=ucon[3]=0.0;

    ucon2pr(WHICHVEL,ucon,ptrgeom,pr);
  }

  extern void filterffde(int i, int j, int k, FTYPE *pr);

  COMPFULLLOOP{
    filterffde(i,j,k,GLOBALMAC(pglobal,i,j,k));
  }

}


//#define FCOVDERTYPE DIFFNUMREC
#define FCOVDERTYPE DIFFGAMMIE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define FCOVDXDELTA 1E-5
#elif(REALTYPE==LONGDOUBLETYPE)
#define FCOVDXDELTA 1E-6
#endif

void Fcov_numerical(int whichcoord, FTYPE *X, FTYPE (*Fcov)[NDIM])
{
  int j,k,l;
  FTYPE Xhk[NDIM], Xlk[NDIM];
  FTYPE Xhj[NDIM], Xlj[NDIM];
  FTYPE Vhk[NDIM], Vlk[NDIM];
  FTYPE Vhj[NDIM], Vlj[NDIM];
  FTYPE mcovhj,mcovlj,kcovhj,kcovlj;
  FTYPE mcovhk,mcovlk,kcovhk,kcovlk;
  FTYPE mcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j); // i not used
  FTYPE kcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j); // i not used
  extern int dfridr(FTYPE (*func)(struct of_geom *,FTYPE*,int,int), struct of_geom *ptrgeom, FTYPE *X,int ii, int jj, int kk, FTYPE *ans);


  // setup dummy grid location since dxdxp doesn't need to know if on grid since don't store dxdxp (needed for dfridr())
  struct of_geom geom;
  struct of_geom *ptrgeom;
  ptrgeom=&geom;
  ptrgeom->i=0;
  ptrgeom->j=0;
  ptrgeom->k=0;
  ptrgeom->p=NOWHERE;

  // get V in case whichcoord!=PRIMECOORDS
  //  FTYPE V[NDIM];
  //  FTYPE dxdxp[NDIM][NDIM];
  //  FTYPE idxdxp[NDIM][NDIM];
  //  if(whichcoord!=PRIMECOORDS){
  //    bl_coord(X, V);
  //    dxdxprim(X, V, dxdxp);
  //    idxdxprim(dxdxp, idxdxp);
  //  }
  //  else{
  //    DLOOPA(j) V[j]=X[j];
  //    // PRIMECOORDS, then no transformation
  //    DLOOP(j,k) idxdxp[j][k]=0.0;
  //    DLOOPA(j) idxdxp[j][j]=1.0;
  //  }


  // GET Fcov
  if(FCOVDERTYPE==DIFFGAMMIE){

    for(k=0;k<NDIM;k++){
      for(j=0;j<NDIM;j++){

	  for(l=0;l<NDIM;l++) Xlk[l]=Xhk[l]=Xlj[l]=Xhj[l]=X[l]; // location of derivative
	  Xhk[k]+=FCOVDXDELTA; // shift up
	  Xlk[k]-=FCOVDXDELTA; // shift down

	  Xhj[j]+=FCOVDXDELTA; // shift up
	  Xlj[j]-=FCOVDXDELTA; // shift down

      if(whichcoord!=PRIMECOORDS){
        bl_coord(Xhk, Vhk);
        bl_coord(Xlk, Vlk);
        bl_coord(Xhj, Vhj);
        bl_coord(Xlj, Vlj);
      }
      else{
        int jj;
        DLOOPA(jj){
          Vhk[jj]=Xhk[jj];
          Vlk[jj]=Xlk[jj];
          Vhj[jj]=Xhj[jj];
          Vlj[jj]=Xlj[jj];
        }
      }


	  //	  dualfprintf(fail_file,"got here1: k=%d j=%d\n",k,j);

	  
	  mcovhj=mcov_func_mcoord(ptrgeom,Xhk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.1: k=%d j=%d\n",k,j);
	  mcovlj=mcov_func_mcoord(ptrgeom,Xlk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.2: k=%d j=%d\n",k,j);
	  mcovhk=mcov_func_mcoord(ptrgeom,Xhj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.3: k=%d j=%d\n",k,j);
	  mcovlk=mcov_func_mcoord(ptrgeom,Xlj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.4: k=%d j=%d\n",k,j);

	  kcovhj=kcov_func_mcoord(ptrgeom,Xhk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.5: k=%d j=%d\n",k,j);
	  kcovlj=kcov_func_mcoord(ptrgeom,Xlk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.6: k=%d j=%d\n",k,j);
	  kcovhk=kcov_func_mcoord(ptrgeom,Xhj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.7: k=%d j=%d\n",k,j);
	  kcovlk=kcov_func_mcoord(ptrgeom,Xlj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.8: k=%d j=%d\n",k,j);
      
      //      dualfprintf(fail_file,"got here2: j=%d k=%d idxdxp=%g %g\n",j,k,idxdxp[j][k],idxdxp[k][j]);

      aforwald=a;

      // F_{\mu\nu} = A_\nu,\mu - A_\mu,\nu
      // A_\phi only
      Fcov[j][k]=0.0;

      if(j==WALDWHICHACOV || WALDWHICHACOV==-1){
        Fcov[j][k] += 0.5*B0WALD*(
                                 +(mcovhj - mcovlj) / (Vhk[k] - Vlk[k])
                                 +2.0*aforwald*(
                                                +(kcovhj - kcovlj) / (Vhk[k] - Vlk[k])
                                                )
                                 );
      }
      if(k==WALDWHICHACOV || WALDWHICHACOV==-1){
        Fcov[j][k] += 0.5*B0WALD*(
                                 -(mcovhk - mcovlk) / (Vhj[j] - Vlj[j])
                                 +2.0*aforwald*(
                                                -(kcovhk - kcovlk) / (Vhj[j] - Vlj[j])
                                                )
                                 );
      }
      }// j
    }// k
  }
  else if(FCOVDERTYPE==DIFFNUMREC){

    aforwald=a;

    for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
        Fcov[j][k]=0.0;

        // 0 in dfridr not used
        FTYPE ans1; dfridr(mcov_func_mcoord,ptrgeom,X,0,j,k,&ans1);
        FTYPE ans2; dfridr(mcov_func_mcoord,ptrgeom,X,0,k,j,&ans2);
        FTYPE ans3; dfridr(kcov_func_mcoord,ptrgeom,X,0,j,k,&ans3);
        FTYPE ans4; dfridr(kcov_func_mcoord,ptrgeom,X,0,k,j,&ans4);
        if(j==WALDWHICHACOV || WALDWHICHACOV==-1){
          Fcov[j][k] += B0WALD*(ans1 +2.0*aforwald*(ans3));
        }
        if(k==WALDWHICHACOV || WALDWHICHACOV==-1){
          Fcov[j][k] += B0WALD*(ans2 +2.0*aforwald*(ans4));
        }
      }

    dualfprintf(fail_file,"NOT SETUP FOR NUMREC\n");
    myexit(23483466);
    
  }// end DIFFNUMREC

}// end Fcov_numerical()





// returns MCOORD m_\mu form of m^\mu={0,0,0,1} value for jth element
FTYPE mcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int ii, int jj) // i not used
{
  FTYPE mcon[NDIM];
  FTYPE mcov[NDIM];

  int getprim=0;
  int whichcoord=MCOORD;
  gset_X(getprim, whichcoord, ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, ptrgeom);

  //  dualfprintf(fail_file,"got here3.3\n");
  mcon[TT]=0.0;
  mcon[RR]=0.0;
  mcon[TH]=0.0;
  mcon[PH]=1.0;
  //  dualfprintf(fail_file,"got here3.4\n");

  // lower only needs geom->gcov
  lower_vec(mcon,ptrgeom,mcov);
  //  dualfprintf(fail_file,"got here3.5\n");

  return(mcov[jj]);
}

// returns MCOORD k_\mu form of k^\mu={1,0,0,0} value for jth element
FTYPE kcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int ii, int jj) // i not used
{
  FTYPE kcon[NDIM];
  FTYPE kcov[NDIM];

  int getprim=0;
  int whichcoord=MCOORD;
  gset_X(getprim, whichcoord, ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, ptrgeom);

  kcon[TT]=1.0;
  kcon[RR]=0.0;
  kcon[TH]=0.0;
  kcon[PH]=0.0;

  // lower only needs geom->gcov
  lower_vec(kcon,ptrgeom,kcov);

  return(kcov[jj]);
}

void Jcon_numerical(int whichcoord, FTYPE *X, FTYPE *Jcon)
{
  int j,k,l;
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE Vh[NDIM], Vl[NDIM];
  FTYPE Fconh,Fconl;
  FTYPE Fcon_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  extern int dfridr(FTYPE (*func)(struct of_geom *,FTYPE*,int,int), struct of_geom *ptrgeom, FTYPE *X,int ii, int jj, int kk, FTYPE *ans);


  // setup dummy grid location since dxdxp doesn't need to know if on grid since don't store dxdxp (needed for dfridr())
  struct of_geom geom;
  struct of_geom *ptrgeom;
  ptrgeom=&geom;
  ptrgeom->i=0;
  ptrgeom->j=0;
  ptrgeom->k=0;
  ptrgeom->p=NOWHERE;




  if(FCOVDERTYPE==DIFFGAMMIE){

    for(k=0;k<NDIM;k++){
      
      Jcon[k] = 0;
      for(j=0;j<NDIM;j++){

        for(l=0;l<NDIM;l++) Xl[l]=Xh[l]=X[l]; // location of derivative
        Xh[j]+=FCOVDXDELTA; // shift up
        Xl[j]-=FCOVDXDELTA; // shift down

        if(whichcoord!=PRIMECOORDS){
          bl_coord(Xh, Vh);
          bl_coord(Xl, Vl);
        }
        else{
          int jj;
          DLOOPA(jj){
            Vh[jj]=Xh[jj];
            Vl[jj]=Xl[jj];
          }
        }
        
        // F^{kj}
        Fconh=Fcon_func_mcoord(ptrgeom,Xh,k,j);
        Fconl=Fcon_func_mcoord(ptrgeom,Xl,k,j);
	
        // J^\mu = {F^{\mu\nu}}_{;\nu}
        // \detg J^k = F^{kj}_{;j} = (\detg F^{kj})_{,j} <--- thing actually computed and returned
        Jcon[k] += (Fconh - Fconl) / (Vh[j] - Vl[j]) ;
      } // j
    }// k

  }
  else if(FCOVDERTYPE==DIFFNUMREC){

    for(k=0;k<NDIM;k++){
      Jcon[k] = 0;
      for(j=0;j<NDIM;j++){
        FTYPE ans1; dfridr(Fcon_func_mcoord,ptrgeom,X,0,j,k,&ans1);
        Jcon[k]+=ans1;
      }
    }

    // get V in case whichcoord!=PRIMECOORDS
    FTYPE V[NDIM];
    FTYPE dxdxp[NDIM][NDIM];
    FTYPE idxdxp[NDIM][NDIM];
    if(whichcoord!=PRIMECOORDS){
      bl_coord(X, V);
      dxdxprim(X, V, dxdxp);
      idxdxprim(dxdxp, idxdxp);
    }
    else{
      DLOOPA(j) V[j]=X[j];
      // PRIMECOORDS, then no transformation
      DLOOP(j,k) idxdxp[j][k]=0.0;
      DLOOPA(j) idxdxp[j][j]=1.0;
    }

    // get whichcoord version -- assumes Fcon is in whichcoord
    FTYPE Jconorig[NDIM];
    DLOOPA(j) Jconorig[j]=Jcon[j];
    DLOOPA(j) Jcon[j]=0.0;
    DLOOP(j,k){
      Jcon[k] += Jconorig[j] * idxdxp[k][j];
    }

  }





}

#undef FCOVDERTYPE
#undef FCOVDXDELTA




// returns MCOORD F^{ii jj}
FTYPE Fcon_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int ii, int jj)
{
  int getprim=0;
  int whichcoord=MCOORD;
  gset_X(getprim, whichcoord, ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X, ptrgeom);

  FTYPE Fcov[NDIM][NDIM];
  void Fcov_numerical(int whichcoord, FTYPE *X, FTYPE (*Fcov)[NDIM]);
  Fcov_numerical(MCOORD, X, Fcov);

  // get covariant Maxwell from contravariant
  FTYPE Fcon[NDIM][NDIM];
  extern void raise_A(FTYPE (*Acov)[NDIM], struct of_geom *geom, FTYPE (*Acon)[NDIM]);//extern void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM]);
  raise_A(Fcov, ptrgeom, Fcon) ;

  // return gdet*F^{ii jj}
  return(ptrgeom->gdet*Fcon[ii][jj]);
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

  int getmax_densities_full(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax, SFTYPE *uradmax, SFTYPE *utotmax, SFTYPE *pmax, SFTYPE *pradmax, SFTYPE *ptotmax);
  int funreturn=getmax_densities_full(prim,&rhomax,&umax,&uradmax,&utotmax,&pmax,&pradmax,&ptotmax);


  if(0){
    FTYPE parms[MAXPASSPARMS];
    int eqline;

    eqline=1;
    parms[0]=rin;
    parms[1]=rhodisk;

    funreturn+=user1_normalize_densities(eqline, parms, prim, &rhomax, &umax);
  }

  return(funreturn);
}



// assumes we are fed the true densities
int getmax_densities_full(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax, SFTYPE *uradmax, SFTYPE *utotmax, SFTYPE *pmax, SFTYPE *pradmax, SFTYPE *ptotmax)
{
  int funreturn;

  funreturn=user1_getmax_densities_full(prim,rhomax, umax, uradmax, utotmax, pmax, pradmax, ptotmax);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}


// get maximum b^2 and p_tot
int get_maxes(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *ptot_max, FTYPE *beta_min)
{
  int funreturn;
  int eqslice=0;
  FTYPE parms[MAXPASSPARMS];

  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();
  
  if(FIELDTYPE==VERTFIELD || FIELDTYPE==BLANDFORDQUAD || FIELDTYPE==DISK2FIELD || FIELDTYPE==FIELDJONMAD || FIELDTYPE==FIELDWALD){
    eqslice=1;
  }
  else{
    eqslice=0;
  }

  if(FIELDTYPE==FIELDJONMAD || FIELDTYPE==FIELDWALD){
    parms[0]=rinfield;
    parms[1]=routfield;
  }
  else{
    parms[0]=rinfield;
    parms[1]=routfield;
  }

  parms[2]=0.05; // for THETAEQ for near equator only so easy to understand how beta enters.  Still should check vertical distribution.

  funreturn=user1_get_maxes(eqslice, parms,prim, bsq_max, ptot_max, beta_min);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;

  int set_fieldtype(void);
  int FIELDTYPE=set_fieldtype();

  if(FIELDTYPE!=NOFIELD && FIELDTYPE!=FIELDWALD && DOWALDDEN==0){
    dualfprintf(fail_file,"DID NORM FIELD\n");
    
    funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
    if(funreturn!=0) return(funreturn);
  }  
 
  return(0);

}



#undef SLOWFAC


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

  atmospheretype=1; // default
 
  funreturn=user1_set_atmosphere(atmospheretype, whichcond, whichvel, ptrgeom, pr);
  if(funreturn!=0) return(funreturn);
 
  return(0);

}



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor, FTYPE *prceiling)
{
  int funreturn;

  int pliter,pl;
  PLOOP(pliter,pl){
    prfloor[RHO]=RHOMINLIMIT;
    prfloor[UU]=UUMINLIMIT;

    prceiling[RHO]=BIG;
    prceiling[UU]=BIG;

    if(PRAD0>=0){ 
      prfloor[PRAD0]=ERADLIMIT;
      prceiling[PRAD0]=BIG;
    }
  }

  // default is for spherical flow near BH
  if(WHICHPROBLEM==RADDONUT){
    // KORALTODO: floor currently causes injection of hot matter and run-away problems with radiation.
    funreturn=set_density_floors_default(ptrgeom, pr, prfloor, prceiling);

    
    // absolute floor, because when magnetic field is zero in some region, then density can go to zero rather than being limited, and then radiation can radically push it around or density gradient can be too extreme for code.
    FTYPE V[NDIM];
    bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,V);
    //    FTYPE lowcoef=MIN(1E-7,RADNT_RHOATMMIN/RADNT_RHODONUT);
    FTYPE lowcoef=MIN(1E-9,RADNT_RHOATMMIN/RADNT_RHODONUT); // SUPERMADNEW
    FTYPE lowpow=-2.0;
    FTYPE rholimit=RADNT_RHODONUT*(lowcoef*pow(V[1],lowpow));
    if(pr[RHO]<rholimit) pr[RHO]=rholimit;

    if(funreturn!=0) return(funreturn);
  }

  return(0);
}

int set_density_floors_alt(struct of_geom *ptrgeom, struct of_state *q, FTYPE *pr, FTYPE *U, FTYPE bsq, FTYPE *prfloor, FTYPE *prceiling)
{
  int funreturn;

  int pliter,pl;
  PLOOP(pliter,pl){
    prfloor[RHO]=RHOMINLIMIT;
    prfloor[UU]=UUMINLIMIT;

    prceiling[RHO]=BIG;
    prceiling[UU]=BIG;

    if(PRAD0>=0){ 
      prfloor[PRAD0]=ERADLIMIT;
      prceiling[PRAD0]=BIG;
    }
  }

  // default is for spherical flow near BH
  if(WHICHPROBLEM==RADDONUT){
    // KORALTODO: floor currently causes injection of hot matter and run-away problems with radiation.
    funreturn=set_density_floors_default_alt(ptrgeom, q, pr, U, bsq, prfloor, prceiling);

    FTYPE V[NDIM];
    bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,V);
    //   FTYPE lowcoef=MIN(1E-7,RADNT_RHOATMMIN/RADNT_RHODONUT);
    FTYPE lowcoef=MIN(1E-9,RADNT_RHOATMMIN/RADNT_RHODONUT); // SUPERMADNEW
    FTYPE lowpow=-2.0;
    FTYPE rholimit=RADNT_RHODONUT*(lowcoef*pow(V[1],lowpow));
    if(pr[RHO]<rholimit) pr[RHO]=rholimit;



    if(funreturn!=0) return(funreturn);
  }

  return(0);
}






// Setup problem-dependent grid sectioning
int theproblem_set_enerregiondef(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // Torus problem
  //  torus_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);
  //jet_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);

  if(1){
    enerregiondef[POINTDOWN][1]=0;
    enerregiondef[POINTUP][1]=totalsize[1]-1;
    enerregiondef[POINTDOWN][2]=0;
    enerregiondef[POINTUP][2]=totalsize[2]-1;
    enerregiondef[POINTDOWN][3]=0;
    enerregiondef[POINTUP][3]=totalsize[3]-1;
  }

  if(DOGRIDSECTIONING && WHICHPROBLEM==RADCYLJET && RADCYLJET_TYPE==6){

    
    // see advance.c: Whether to allow shift in evolved quantities to preserve conservation and divb=0.  Set to zero if exposing that surface in time.  Set to 1 if absorbing that surface in time and relying on it to inject a solution.
    AVOIDADVANCESHIFTX1DN= 1;
    AVOIDADVANCESHIFTX1UP= 1;
    AVOIDADVANCESHIFTX2DN= 1;
    AVOIDADVANCESHIFTX2UP= 1;
    AVOIDADVANCESHIFTX3DN= 1;
    AVOIDADVANCESHIFTX3UP= 1;
    GLOBALBCMOVEDWITHACTIVESECTION= 1;



    enerregiondef[POINTDOWN][2]=0;
    
    enerregiondef[POINTUP][2]=N2BND*2 + round(((thetime+0.0)/45.0)*((FTYPE)(totalsize[2]-1))*(45.0/Rout_array[2]));
    
    if(enerregiondef[POINTUP][2]>totalsize[2]-1){
      enerregiondef[POINTUP][2]=totalsize[2]-1;
    }
  }

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
  if(N3==1){
    //  *updateeverynumsteps=100;
    *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
  }
  else{
    //  *updateeverynumsteps=10;
    *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
  }

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


void adjust_flux(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

  // X1UP
  if(WHICHPROBLEM==RADCYLJET && RADCYLJET_TYPE==5){
    int i,j,k,pl ;
    FTYPE sth ;
    FTYPE X[NDIM],V[NDIM],r,th ;
    int inboundloop[NDIM];
    int outboundloop[NDIM];
    int innormalloop[NDIM];
    int outnormalloop[NDIM];
    int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
    int riin,riout,rjin,rjout,rkin,rkout;
    int dosetbc[COMPDIM*2];
    int ri;
    int boundvartype=BOUNDFLUXTYPE;

    ////////////////////////
    //
    // set bound loop
    //
    ///////////////////////
    set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);
    //enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
    //  localenerpos=enerposreg[enerregion];



    // enforce true reflective condition on (near) polar axis
    if(mycpupos[1]==0){  
      LOOPX1dir{
        i=0;

        PALLLOOP(pl) if(pl!=U1 && pl!=URAD1) MACP0A1(F1,i,j,k,pl)=0;
      }
    }

    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;

    // prescribe outer boundary condition on wall
    if(mycpupos[1]==ncpux1-1){
      LOOPX1dir{
        i=N1;

        PALLLOOP(pl) if(pl!=U1 && pl!=URAD1) MACP0A1(F1,i,j,k,pl)=0;

        if(1){
          get_geometry(i,j,k,FACE1,ptrgeom) ;

          FTYPE t0=1.0;
          FTYPE interptime=(fluxtime<1.0 ? (fluxtime/t0) : 1.0);
          if(interptime>1.0) interptime=1.0;
          if(interptime<0.0) interptime=0.0;

          //        interptime=0;


          // ensure fake pressure equilibrium
          FTYPE Ehatstar=1.0*RADCYLJET_EHATJET; // 1.0 ensures pressure equlibrium, regardless of actual physics
          FTYPE ustar=1.0*RADCYLJET_UJET; // 1.0 ensures pressure equlibrium, regardless of actual physics
          FTYPE rhostar=1.0*RADCYLJET_RHOJET;

          // but drive radiative flux into jet as if vr0*Ehatstar
          FTYPE vr0; //=-1E-6; //RADCYLJET_VRSTAR;
          //        vr0=0.0;
          // should really set flux directly, or based upon flux-limited diffusion set vr
          // But jump is infinite, but can imagine as if not.

          // set flux of radiation: R^t_x1.  Just set Ehat=Ehat0 since non-rel injection, and compute literal flux
          FTYPE prflux[NPR];
          PALLLOOP(pl) prflux[pl]=0;
          prflux[RHO]=0.0; //rhostar;  // this has no effect
          //          prflux[UU]=ustar; // this term ensures pressure equilibrium
          prflux[UU]=0.0; //MACP0A1(prim,i,j,k,UU);
          prflux[U1]=0.0;
          prflux[U2]=0.0;
          prflux[U3]=0.0;
          prflux[URAD0]=Ehatstar*(1.0-interptime) + interptime*Ehatstar*10.0;

          //          prflux[URAD1]=interptime*vr0/sqrt(fabs(ptrgeom->gcov[GIND(1,1)])); // assume vr0 is orthonormal, so get coordinate out of it.
          FTYPE kappa=calc_kappaes_user(rhostar,1.0,0,0,0);
          vr0=sqrt(1.0/(3.0*kappa*t0));
          static int firsttime=1;
          if(firsttime) dualfprintf(fail_file,"vr0=%21.15g\n",vr0);
          firsttime=0;
          prflux[URAD1]=vr0/sqrt(fabs(ptrgeom->gcov[GIND(1,1)])); // assume vr0 is orthonormal, so get coordinate out of it.
          prflux[URAD2]=0.0;
          prflux[URAD3]=0.0;

          struct of_state q;
          get_state(prflux, ptrgeom, &q);
          FTYPE fluxrad[NPR];
          primtoflux(UEVOLVE, prflux, &q, 1, ptrgeom, fluxrad, NULL);

          //int pliter;
          //PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d fluxrad=%21.15g F1=%21.15g\n",pl,fluxrad[pl],MACP0A1(F1,i,j,k,pl));
          //myexit(0);

          // at least pressure term needs to be there
          //MACP0A1(F1,i,j,k,U1)=fluxrad[U1]; // use original outflowed, so pressure balance maintained regardless of vr
          //MACP0A1(F1,i,j,k,URAD1)=fluxrad[URAD1]; // use original outflowed, so pressure balance maintained regardless of vrad_r
          // will this lead, for large energy flux, inconsistently small momentum flux?

          // setup so exactly pressure balanced, so no momentum flux, implying no motion of boundary.
          //          MACP0A1(F1,i,j,k,U1)=MACP0A1(F1,im1mac(i),j,k,U1);
          //MACP0A1(F1,i,j,k,URAD1)=MACP0A1(F1,im1mac(i),j,k,URAD1);

          // But, even with pressure balance, *still* inject energy.  So while boundary cannot move, we still leak in energy.
          //          MACP0A1(F1,i,j,k,U1)=fluxrad[URAD1];
          MACP0A1(F1,i,j,k,U1)=MACP0A1(F1,im1mac(i),j,k,U1);

          MACP0A1(F1,i,j,k,URAD1)=fluxrad[URAD1];
          MACP0A1(F1,i,j,k,URAD0)=fluxrad[URAD0];

        }
        // IDEAS (for failures)

        // *) Turn on flux slowly so not shock.  Slower.
        // *) Why cylindrical boundary goofy and evolving? (not body forces change, but maybe non-linear conserves)
        //   e.g. gdet R^r_r = gdet (Ehat \gamma^2 vr^2 + Erad/3) = r Erad .  So linear, so should be ok.  Must be body force not right really?

      }
    }
  }


  // X2DN
  if(WHICHPROBLEM==RADCYLJET && RADCYLJET_TYPE==6){
    int i,j,k,pl ;
    FTYPE sth ;
    FTYPE X[NDIM],V[NDIM],r,th ;
    int inboundloop[NDIM];
    int outboundloop[NDIM];
    int innormalloop[NDIM];
    int outnormalloop[NDIM];
    int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
    int riin,riout,rjin,rjout,rkin,rkout;
    int dosetbc[COMPDIM*2];
    int ri;
    int boundvartype=BOUNDFLUXTYPE;

    ////////////////////////
    //
    // set bound loop
    //
    ///////////////////////
    set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);
    //enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
    //  localenerpos=enerposreg[enerregion];


    // X1DN
    // enforce true reflective condition on (near) polar axis
    if(mycpupos[1]==0){  
      LOOPX1dir{
        i=0;

        PALLLOOP(pl) if(pl!=U1 && pl!=URAD1) MACP0A1(F1,i,j,k,pl)=0;
      }
    }

    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;

    // prescribe outer boundary condition on wall
    if(mycpupos[2]==0){
      LOOPX2dir{
        j=0;

        //        PALLLOOP(pl) if(pl!=U2 && pl!=URAD2) MACP0A1(F2,i,j,k,pl)=0;

        if(1){
          get_geometry(i,j,k,FACE2,ptrgeom) ;
          bl_coord_ijk(i,j,k,FACE2,V);

          FTYPE t0=1.0;
          FTYPE interptime=(fluxtime<1.0 ? (fluxtime/t0) : 1.0);
          if(interptime>1.0) interptime=1.0;
          if(interptime<0.0) interptime=0.0;

          FTYPE prflux[NPR];
          extern int jetbound(int i, int j, int k, int loc, FTYPE *prin, FTYPE *prflux, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
          int insidejet=jetbound(i,j,k,FACE2,MAC(prim,i,j,k),prflux,prim);
          

          if(insidejet==0){
            //PALLLOOP(pl) if(pl!=U2 && pl!=URAD2) MACP0A1(F2,i,j,k,pl)=0;
            //            pl=U2; MACP0A1(F2,i,j,k,pl)=MACP0A1(F2,i,jp1mac(j),k,pl);
            //            pl=URAD2; MACP0A1(F2,i,j,k,pl)=MACP0A1(F2,i,jp1mac(j),k,pl);
          }
          else{
            struct of_state q;
            get_state(prflux, ptrgeom, &q);
            FTYPE flux[NPR];
            primtoflux(UEVOLVE, prflux, &q, 2, ptrgeom, flux, NULL);

            // assign flux
            PALLLOOP(pl) if(pl!=U2 && pl!=URAD2) MACP0A1(F2,i,j,k,pl)=0;
            PALLLOOP(pl) MACP0A1(F2,i,j,k,pl)=flux[pl];
          }

        }
        // IDEAS

      }
    }
  }

}

int jetbound(int i, int j, int k, int loc, FTYPE *prin, FTYPE *prflux, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  get_geometry(i,j,k,loc,ptrgeom) ;
  FTYPE X[NDIM],V[NDIM],r,th ;
  bl_coord_ijk(i,j,k,loc,V);



  int insidejet;
  FTYPE width=1.0;


  if(V[1]<=width){ // ||1 if want fixed conditions outside of jet as well.
    insidejet=1;
  }
  else{
    insidejet=0;
  }

  int pliter,pl;
  if(insidejet==0){
    PALLLOOP(pl) prflux[pl]=prin[pl];

  }
  else{


    // these are really star values
    FTYPE Ehatstar=RADCYLJET_EHATJET;
    FTYPE ustar=RADCYLJET_UJET;
    FTYPE Tstar=RADCYLJET_TEMPJET;
    FTYPE pradstar; pradstar=(4.0/3.0-1.0)*calc_LTE_EfromT(Tstar);
    FTYPE ujet;
    if(WHICHRADSOURCEMETHOD==SOURCEMETHODNONE) ujet=ustar*1E-2*99.99;
    else ujet=ustar*1E-2;

 
    FTYPE rhostar=RADCYLJET_RHOJET;
    FTYPE rhojet; rhojet=rhostar*1E-2;
    
    //    FTYPE vz0;
    

    //            FTYPE interpi=exp(-V[1]*V[1]/(2.0*width*width));

    // gamma^2 = 1+usq -> usq = gamma^2-1 -> |u| = sqrt(gamma^2-1)
    FTYPE gammacore; gammacore = 7.0;
    FTYPE uzcore; uzcore=sqrt(gammacore*gammacore-1.0);
    FTYPE rhocore; rhocore=rhojet;
    
    FTYPE efluxcore; efluxcore=rhocore*uzcore*uzcore;

    FTYPE eflux; eflux = efluxcore*stepfunctionab2(V[1],1,0);

    FTYPE uz=uzcore*stepfunctionab2(V[1],1,0); // orthonormal uz
    // don't use stepfunctionab() on vz because then drops in gamma very fast.
    // gamma^2 = 1/(1-v^2) -> 1-v^2 = 1/gamma^2 -> v^2 = 1-1/gamma^2 -> |v| = sqrt(1-1/gamma^2) or |v| = uz0/gamma
    //    FTYPE vz0=sqrt(1.0-1.0/(gamma*gamma));
    FTYPE rho; rho = eflux/(uz*uz);
    if(rho>rhostar) rho=rhostar;

    //stepfunctionab2(V[1],rhojet,rhostar); // so that rho drops faster so no density edge next to jet
    //FTYPE ug=stepfunctionab(V[1],ujet,ustar);


    FTYPE Tjet = 1.62*Tstar;
    FTYPE Temp;
    Temp = stepfunctionab2(V[1],Tjet,Tstar);
    FTYPE ug;
    ug = Temp*rho/(gam-1.0);
    //    FTYPE ug=stepfunctionab(V[1],ujet,ustar);
    //FTYPE Tstar=1.0E7/TEMPBAR;
    //          // P = (arad/3)T^4 + rho T
    //          pr[UU]=u_rho0_T_simple(i,j,k,CENT,pr[RHO],Tstar);
    //          pr[URAD0]=calc_LTE_EfromT(Tstar);

    // set flux of radiation: R^t_x1.  Just set Ehat=Ehat0 since non-rel injection, and compute literal flux
    PALLLOOP(pl) prflux[pl]=0;
    prflux[RHO]=rho;
    prflux[UU]=ug;
    FTYPE Tgas=compute_temp_simple(i,j,k,loc,prflux[RHO],prflux[UU]);

    prflux[U1]=0.0; //MACP0A1(prim,i,jp1mac(j),k,U1);
    prflux[U2]=uz/sqrt(fabs(ptrgeom->gcov[GIND(2,2)]));
    prflux[U3]=0.0; //MACP0A1(prim,i,jp1mac(j),k,U3);

    prflux[URAD0]=calc_LTE_EfromT(Tgas); // Thermal equilibrium  // Ehatstar;
    prflux[URAD1]=0.0; //MACP0A1(prim,i,jp1mac(j),k,URAD1);
    prflux[URAD2]=prflux[U2]; // optically thick
    prflux[URAD3]=0.0; //MACP0A1(prim,i,jp1mac(j),k,URAD3);

    FTYPE ucon[NDIM];
    FTYPE others[NUMOTHERSTATERESULTS];
    ucon_calc(prflux,ptrgeom,ucon,others);

    static int firsttime=1;
    if(firsttime==1){
      // cs2 ~ gam*Ptot/rho and need vz0>cs.
      // look at sqrt(cs2tot) in SM to check or obtain here.
      FTYPE ptot=((gam-1.0)*prflux[UU] + (4.0/3.0-1.0)*prflux[URAD0]);
      FTYPE ptrue;
      if(WHICHRADSOURCEMETHOD==SOURCEMETHODNONE) ptrue=(gam-1.0)*prflux[UU];
      else ptrue=ptot;
      FTYPE gamptot=(gam*(gam-1.0)*prflux[UU] + (4.0/3.0)*(4.0/3.0-1.0)*prflux[URAD0]);
      FTYPE cs2 = gamptot/prflux[RHO];
      dualfprintf(fail_file,"JET: ptrue=%21.15g ptot=%21.15g cs2=%21.15g cs=%21.15g\n",ptrue,ptot,cs2,sqrt(cs2));

      dualfprintf(fail_file,"Tgasinj=%21.15g u/rho=%21.15g prad0/rho=%21.15g\n",Tgas,prflux[UU]/prflux[RHO],prflux[URAD0]/prflux[RHO]);
      dualfprintf(fail_file,"TEST: %21.15g\n",0.99*stepfunctionab(1E-2,1.0,0.0));
      dualfprintf(fail_file,"TESTUCON0: %21.15g %21.15g %21.15g %21.15g\n",ucon[TT],ucon[1]*sqrt(ptrgeom->gcov[GIND(1,1)]),ucon[2]*sqrt(ptrgeom->gcov[GIND(2,2)]),ucon[3]*sqrt(ptrgeom->gcov[GIND(3,3)]));
      dualfprintf(fail_file,"Need: pradstar/rhojet=%21.15g > %21.15g and taujet=%21.15g >>1\n",pradstar/prflux[RHO],ucon[TT]*ucon[TT]*calc_kappaes_user(prflux[RHO],Temp,0,0,0)*width*width/(2.0*Rout_array[2]),calc_kappaes_user(prflux[RHO],Temp,0,0,0)*width);

      FTYPE ljet,ljet1,ljet2,ljet3;
      ljet1 = prflux[RHO]*ucon[TT]*ucon[TT]*M_PI*width*width;
      ljet2 = (prflux[UU] + (gam-1.0)*prflux[UU])*ucon[TT]*ucon[TT]*M_PI*width*width;
      ljet3 = (prflux[URAD0] + (4.0/3.0-1.0)*prflux[URAD0])*ucon[TT]*ucon[TT]*M_PI*width*width;
      ljet=ljet1+ljet2+ljet3;
      

      FTYPE lrad,lradalt1,lradalt2;
      lrad =      (pradstar/calc_kappaes_user(prflux[RHO],Temp,0,0,0))*(2.0*M_PI*Rout_array[2]); // best estimate, except rho and T vary inside jet so only applies to core of jet evaluated first.  Assumes diffusion scale is width, and it could be smaller in which case flux closer to lradalt2.
      lradalt1 =  (pradstar/calc_kappaes_user(rhostar,Tstar,0,0,0))*(2.0*M_PI*Rout_array[2]); // underestimate
      lradalt2 = (ARAD_CODE*pow(Tstar,4.0)/4.0)*(2.0*M_PI*width*Rout_array[2]); // overestimate

      dualfprintf(fail_file,"Need: ljet=%21.15g < lrad=%21.15g  lradalt1=%21.15g  lradalt2=%21.15g\n",ljet,lrad,lradalt1,lradalt2);
      dualfprintf(fail_file,"ljet1=%21.15g ljet2=%21.15g ljet3=%21.15g\n",ljet1,ljet2,ljet3);

      //sm: set lrad=(-Rud01)*dx2*dx3*gdet if(ti==51)
      //sm: set lradtot=SUM(lrad) print {lradtot}

      firsttime=0;
    }
  }
  return(0);
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
  int i,j,k;

  if(DOWALDDEN&&0){
    if(mycpupos[1]==ncpux1-1){
      dualfprintf(fail_file,"did modify\n");
      FULLLOOP{
        if(i==N1){
          MACP1A1(fluxvec,2,N1,j,k,B1)=MACP1A1(fluxvec,1,N1,j,k,B2)=0.0;
        }
      }
    }
  }
  else{
    // not used
  }
}


//**********************************************************************
//******* user opacities ****************************************************
//**********************************************************************

//////////////////////////////
// ff and s (but means something different for ff and s)
//listlog10xi = Range[-8, 8, 1];
//listxi = 10^listlog10xi;

//step1mexpf = 0.2
//listlog101mexpf = Range[-3, 0, step1mexpf];
//list1mexpf = 10^listlog101mexpf;
//list1mexpf = Sort[Join[{0}, list1mexpf]];
//listexpf = 1 - list1mexpf;
//////////////////////////////


//////////////////////////////
// dc:
//step1mexpf = 0.2*2
//listlog101mexpf = Range[-3, 0, step1mexpf];
//list1mexpf = 10^listlog101mexpf;
//list1mexpf = Sort[Join[{0}, list1mexpf]];
//listexpf = 1 - list1mexpf;

//stepthetag = 0.5*2
//listlogthetag = Range[-4, 1, stepthetag];
//listthetag = 10^listlogthetag;

//stepthetae = 0.5
//listlogthetae = Range[-5, 2, stepthetae];
//listthetae = 10^listlogthetae;
//////////////////////////////


// as long as scaled-out prefactor, then near Te=Tg and \mu=0 should be order unity and shouldn't be beyond float limit of 1E-30
// and accuracy of double not needed even if iterations in implicit solver required double, because calculations are not that accurate.
//#define TBLTYPE static const float
//#include "opacitytables.c"

// Steps:
//
// 1) For .nb's in #2, choose MBH, etc. so cut-off chosen and chosen same way for all .nb's.
// 1.5) Ensure size of tables at bottom is consistent as hardcoded size
// 2) run jon@physics-179:/data/jon/pseudotensor@gmail.com/harm_math/opacityfits\ doublecompton_opacity_mu_cutoff2.nb and freefree_opacity_mu_cutoff.nb and synch_opacity_mu_cutoff.nb (with picktemp=0,whichtemp=1 and picktemp=1,whichtemp=1 and picktemp=1,whichtemp=3)
// 3) Run "sed" at bottom of each .nb after everything done.
// 4) In opacityfits directory run: cat energyopacityffmu.dat4 numberopacityffmu.dat4 energyopacitydcmu.dat4 numberopacitydcmu.dat4 energyopacityedcmu.dat4 numberopacityedcmu.dat4 energyopacitysamu.dat4 numberopacitysamu.dat4  energyopacitysbmu.dat4 numberopacitysbmu.dat4  energyopacityscmu.dat4 numberopacityscmu.dat4  > opacitytables.temp.c ; sed 's/"//g' opacitytables.temp.c > opacitytables.c
// 5) copy opacitytables.c to local harmgit directory

// Opal table in GN93hz
//TABLE # 73  $G&N'93 Solar$    X=0.7000 Y=0.2800 Z=0.0200 dXc=0.0000 dXo=0.0000
//#include "opacityopaltables.c"




// KAPPAUSER is optical depth per unit length per unit rest-mass energy density
// calc_kappa_user and calc_kappan_user and calc_kappaes_user return optical depth per unit length.

// energy absorption
FTYPE calc_kappa_user(FTYPE rho, FTYPE B, FTYPE Tg,FTYPE Tr,FTYPE x,FTYPE y,FTYPE z)
{
  //  if(WHICHPROBLEM==RADDONUT && nstep>100){
  //    return(0.0);
  //  }
  //  else return(KAPPAUSER(rho,B,Tg,Tr));
  return(KAPPAUSER(rho,B,Tg,Tr));
}

// number absorption
FTYPE calc_kappan_user(FTYPE rho, FTYPE B, FTYPE Tg,FTYPE Tr,FTYPE x,FTYPE y,FTYPE z)
{
  //  if(WHICHPROBLEM==RADDONUT && nstep>100){
  //    return(0.0);
  //  }
  //  else return(KAPPANUSER(rho,B,Tg,Tr));
  return(KAPPANUSER(rho,B,Tg,Tr));
}

//scattering
FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z)
{  
  return(KAPPAESUSER(rho,T));

}

// User's cooling function:

#define USERTHETACOOL       (h_over_r)	/* should be same as h_over_r */
#define USERTAUCOOL         (2.0*M_PI*h_over_r)	        /* cooling time in number of rotational times : really USERTAUCOOL=2*M_PI would be 1 rotational time */
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









#define JET6LIKEUSERCOORD 0
#define UNIHALFUSERCOORD 1
#define ORIGWALD 2

#define WHICHUSERCOORD JET6LIKEUSERCOORD


// for defcoord=JET6COORDS like USERCOORDS
static FTYPE npow,r1jet,njet1,njet,r0jet,rsjet,Qjet, ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3, rs, r0,npow2,cpow2,rbr,x1br, h0,cpow3; 


void set_coord_parms_nodeps_user(int defcoordlocal)
{
  if(1){
    // see jet3coords_checknew.nb

    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    npow=1.0;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid


    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    //    cpow3=0.01;
    cpow3=1.0;
    //radius at which hyperexponentiation kicks in
    //    rbr = 1E3;
    if(DOWALDDEN){
      rbr = 5E7; // WALD 5E2->5E7
      //power exponent
      npow2=4.0; // WALD: 6.0->4.0
    }
    else{
      rbr = 5E2; // WALD 5E2->5E7
      //power exponent
      npow2=6.0; // WALD: 6.0->4.0
    }



    // must be same as in dxdxp()
    // GODMARK: Note njet here is overwritten by njet later, but could have been different values if setup variable names differently.
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(1){
      r1jet=2.8;
      njet=0.3;
      r0jet=15.0;
      rsjet=40.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
    }

    // for switches from normal theta to ramesh theta
    rs=40.0; // shift
    r0=20.0; // divisor
 
    // for theta1
    //    hslope=0.3 ; // resolve inner-radial region near equator
    r0jet3=20.0; // divisor
    rsjet3=0.0; // subtractor

    // for theta2
    h0=0.3; // inner-radial "hslope" for theta2
    //h0=0.1; // inner-radial "hslope" for theta2 // for thinner disks, change this.
    // GODMARK: Note that this overwrites above njet!
    // power \theta_j \propto r^{-njet}
    if(DOWALDDEN==1){
      njet=1.0;
    }
    else if(DOWALDDEN==2) njet=0.0;


    // see fix_3dpoledtissue.nb
#if(0)
    ntheta=21.0;
    htheta=0.15;
    rsjet2=5.0;
    r0jet2=2.0;
#else
    ntheta=5.0;
    htheta=0.15;
    rsjet2=5.0;
    r0jet2=2.0;
#endif

  }

}


void set_coord_parms_deps_user(int defcoordlocal)
{
  if(1){
    x1br = log( rbr - R0 ) / npow;  //the corresponding X[1] value
  }

}

void write_coord_parms_user(int defcoordlocal, FILE *out)
{
  if(1){
    fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,r0jet,rsjet,Qjet,ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3,rs,r0,npow2,cpow2,rbr,x1br,cpow3);
  }

}
void read_coord_parms_user(int defcoordlocal, FILE *in)
{

  if(1){
    fscanf(in,HEADER19IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet,&ntheta,&htheta,&rsjet2,&r0jet2,&rsjet3,&r0jet3,&rs,&r0,&npow2,&cpow2,&rbr,&x1br,&cpow3);

  }
}

void read_coord_parms_mpi_user(int defcoordlocal)
{
  if(1){
#if(USEMPI)
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&ntheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&htheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
#endif
  }

}


void blcoord_user(FTYPE *X, FTYPE *V)
{
  extern FTYPE mysin(FTYPE th);

  if(1){

#if(0) // no change in exponentiation
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],npow)) ;
#elif(WHICHUSERCOORD==UNIHALFUSERCOORD)

    Rout=2000.0;
    theexp = npow*X[1];
    npow=1.0;
    FTYPE gconst1=1.0;
    FTYPE gconst2=gconst1*.000001;
    V[1] = R0 + gconst1*X[1] + gconst2*exp(theexp);


#elif(WHICHUSERCOORD==JET6LIKEUSERCOORD)


#define cr(x) (exp(-1.0/(x)))
#define tr(x) (cr(x)/(cr(x) + cr(1.0-(x))))
#define trans(x,L,R) ((x)<=(L) ? 0.0 : ((x)>=(R) ? 1.0 : tr(((x)-(L))/((R)-(L)))) )
#define transR(x,center,width) ( 0.5*(1.0-tanh(+((x)-(center))/(width))))
#define transL(x,center,width) ( 0.5*(1.0-tanh(-((x)-(center))/(width))))
#define transM(x,center,width) ( exp(-pow(((x)-(center))/((width)*0.5),2.0) ) )

#define plateau(x,L,R,W) (trans(x,(L)-0.5*(W),(L)+0.5*(W))*(1.0-trans(x,(R)-0.5*(W),(R)+0.5*(W))))


    if(1){
      FTYPE theexp = npow*X[1];
      if( X[1] > x1br ) {
        theexp += cpow2 * pow(X[1]-x1br,npow2);
      }
      V[1] = R0+exp(theexp);
    }
    else{

#define line1r(x,w) (Rout)
#define line2r(x,w) (Routeq)
#define line3r(x,w) (Rout))
      //#define wparsam(x,r) (h0 + pow( ((r)-rsjet3)/r0jet3 , -njet))
      //#define wparsam(x,r) (h0 + pow( ((r)-0.0)/4.2 , -njet))
  //#define wparsam(x,r) (h0 + pow(0.15 + ((r)-0.0)/10.0 , -njet))
#define thetajon(x,w,xp1,xp2) (line1r(x,w)*(1.0-trans(x,xp1,xp2)) + line2r(x,w)*trans(x,xp1,xp2))

      //startx[1] = pow(log(Rin-R0),1.0/npow);
      //dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
      //startx[1] = log(Rin-R0)/npow;

      // 
      FTYPE mysx1=log(Rin-R0)/npow;

#define lineeqr(x,w) (R0 + exp(npow*(X[1]-mysx1)*0.70 + npow*mysx1  ) )
#define linepoler(x,w) (R0 + exp(npow*X[1]))
#define thetaLr(x,wp,weq,xp1,xp2) ( linepoler(x,wp)*(1.0-trans(x,xp1,xp2)) + lineeqr(x,weq)*trans(x,xp1,xp2) )
#define thetajon2(x,wp,weq,xp1,xp2) ( x<0.5 ? thetaLr(x,wp,weq,xp1,xp2) : thetaLr(1.0-x,wp,weq,xp1,xp2) )

      FTYPE Routvsx2,R0vsx2;
      FTYPE Routeq = 100.0;

      if(1){

        FTYPE xpole=0.25;
        FTYPE eqslope=0.5;
        FTYPE xeq=0.5;

        V[1] = thetajon2(X[2],0,0,xpole,xeq);

      }

      if(0){
        //      FTYPE poleslope=wparsam(X[2],V[1]);
        FTYPE xpole=0.25;
        //      FTYPE eqslope=0.1;
        FTYPE eqslope=0.5;
        FTYPE xeq=0.5;
        Routvsx2 = thetajon2(X[2], 0.0, 0.0, xpole, xeq);

        FTYPE R0eq=0.0;

        R0vsx2 = R0eq + R0*((Routvsx2-Routeq)/(Rout-Routeq));
      }


      if(0){
        //V[1] = R0+exp(npow*X[1]);
        // go startx[1] to totalsize[1]*dx[1]+startx[1] at pole and equator and everywhere.
        //      But go to reduced value of V[1] near equator related to Rout as will be used at pole via setting of X[1]'s first and last values.
        //FTYPE AA = 1.0 - Routeq/Rout;
        //      V[1] = Rin + (V[1]-Rin)*(1.0 - AA*sin(X[2]*M_PI));
        //      V[1] *= (1.0 - AA*sin(X[2]*M_PI));
        
        Routvsx2 = thetajon(X[2],0.0,0.25,0.75);

        R0vsx2=R0;
        
      }

      //      V[1] = R0 + exp(npow*( (X[1]-startx[1])*Routvsx2/Rout+startx[1])   );

      //      V[1] = R0 + exp(npow*X[1]);

      //      V[1] = R0 + (V[1]-R0)*Routvsx2/Rout + (1.0-Routvsx2/Rout)*(0.38+0.018);

      

      
      if(0){
        V[1] = R0+exp(npow*X[1]);
        // go startx[1] to totalsize[1]*dx[1]+startx[1] at pole and equator and everywhere.
        //      But go to reduced value of V[1] near equator related to Rout as will be used at pole via setting of X[1]'s first and last values.
        //        FTYPE Routeq = 100.0;
        FTYPE AA = 1.0 - Routeq/Rout;
        V[1] = Rin + (V[1]-Rin)*(1.0 - AA*sin(X[2]*M_PI));
      }
    }


    //    FTYPE npowtrue,npowlarger=10.0;
    //    FTYPE npowrs=1E3;
    //    FTYPE npowr0=2E2;
    //    npowtrue = npow + (npowlarger-npow)*(0.5+1.0/M_PI*atan((V[1]-npowrs)/npowr0));
    //    V[1] = R0+exp(pow(X[1],npowtrue)) ;
#elif(0)
    // avoid jump in grid at rbr
    // determine switches
    FTYPE r0rbr=rbr/2.0;
    FTYPE switch0 = 0.5+1.0/M_PI*atan((V[1]-rbr)/r0rbr); // 1 for outer r

    FTYPE V1 = R0+exp(npow*X[1]);
    FTYPE V2 = R0+exp(npow*X[1] + cpow2 * pow(cpow3*(X[1]-x1br*1.0),npow2));

    V[1] = V1*(1.0-switch0) + V2*switch0;

#endif








    FTYPE theta1,theta2,arctan2;


#if(0)
    // JET3COORDS-based:
    FTYPE myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    theta1 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
#else
    // RAMESH BASED
    // myhslope here is h2 in MCAF paper
    //    // h0 here is h3 in MCAF paper
    //FTYPE njetvsr;
    //if(V[1]<rbr) njetvsr=njet;
    //    else njetvsr=njet/(V[1])*rbr;
    //else njetvsr=
    //njetvsr=njet;

    FTYPE localrbr=Rout; //500.0; // rbr;
    //    FTYPE localrbr=rbr;
    FTYPE localrbrr0=100.0; //MAX(100.0,localrbr/2.0);

    FTYPE switch0 = 0.5+1.0/M_PI*atan((V[1]-localrbr)/localrbrr0); // large r
    FTYPE switch2 = 1.0-switch0; // small r

    //    switch0=0.0; switch2=1.0;

    FTYPE myhslope1=h0 + pow( (V[1]-rsjet3)/r0jet3 , njet);
    FTYPE myhslope2=h0 + pow( (localrbr-rsjet3)/r0jet3 , njet);
    FTYPE myhslope = myhslope1*switch2 + myhslope2*switch0;

    //    myhslope = pow(pow(myhslope1,-2.0) + pow(myhslope2,-2.0),-0.5);

    myhslope=myhslope1;

    // determine theta2
    FTYPE myx2;
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    FTYPE th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    // JET3COORDS-based:
    myhslope1=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    myhslope2=2.0-Qjet*pow(localrbr/r1jet,-njet*(0.5+1.0/M_PI*atan(localrbr/r0jet-rsjet/r0jet)));
    myhslope = myhslope1*switch2 + myhslope2*switch0;
    // myhslope here is h0 in MCAF paper
    FTYPE th0 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);


    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rs)/r0); // switch in .nb file // switch0->1 as r->infinity
    switch2 = 1.0-switch0; // for inner radial region

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.
    //    theta1=th0;
    theta1=th2;

#endif

    if(0){
      // fix_3dpoledtissue.nb based:
      theta2 = M_PI*0.5*(htheta*(2.0*X[2]-1.0)+(1.0-htheta)*pow(2.0*X[2]-1.0,ntheta)+1.0);
      
      // generate interpolation factor
      arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-rsjet2)/r0jet2) );
      
      // now interpolate between them
      V[2] = theta2 + arctan2*(theta1-theta2);
    }


    //V[2] = theta1;

    if(1){
      FTYPE fraceq=0.3;
      FTYPE fracpole=(1.0-fraceq)/2.0;
      FTYPE x2p1=0.0+fracpole;
      FTYPE x2p2=1.0-fracpole;
      FTYPE swide=0.04; //1E-1;

      // s(x) = 0.5 + 0.5 tanh((x-a)/b)

      //      FTYPE switchh0 = 0.5+1.0/M_PI*atan((X[2]-x2p1)/swide);
      //      FTYPE switchh2 = 0.5+1.0/M_PI*atan((X[2]-x2p2)/swide);

      FTYPE switchh0 = 0.5+0.5*tanh((X[2]-x2p1)/swide);
      FTYPE switchh2 = 0.5+0.5*tanh((X[2]-x2p2)/swide);

      FTYPE eqh=0.1;
      FTYPE theq = M_PI * X[2] + ((1. - eqh) * 0.5) * mysin(2. * M_PI * X[2]);
    
      FTYPE thup=switchh0*theq + (1.0-switchh0)*theta1;

      V[2]=switchh2*theta1 + (1.0-switchh2)*thup;

    }
    else{
      V[2]=theta1;
    }

    if(0){ // Sam Gralla
      FTYPE transwidth=0.06; //1E-1;
      FTYPE xcent=0.5; // fixed
      FTYPE transR=0.5*(1.0-tanh(+(X[2] - xcent)/transwidth));
      FTYPE transL=0.5*(1.0-tanh(-(X[2] - xcent)/transwidth));

      h0=0.0;
      //      FTYPE wpar=h0 + pow( (V[1]-rsjet3)/r0jet3 , -njet);
      FTYPE wpar=pow( (V[1])/1 , -njet);

      FTYPE line1 = wpar*X[2];
      FTYPE line2 = M_PI + wpar*(X[2]-1.0);

      V[2] = line1*transR + line2*transL;

      //      dualfprintf(fail_file,"X[2]=%g V[2]=%g\n",X[2],V[2]);

    }

    if(1){ // Sam Gralla 2

      h0=0.0;

#define line1(x,w) ((x)*(w))
#define line2(x,w) ((x)*(w)+M_PI-(w))
#define line3(x,w) ((x)*(w))
      //#define wparsam(x,r) (h0 + pow( ((r)-rsjet3)/r0jet3 , -njet))
      //#define wparsam(x,r) (h0 + pow( ((r)-0.0)/4.2 , -njet))
#define wparsam(x,r) (h0 + pow(0.15 + ((r)-0.0)/10.0 , -njet))
      //#define wparsam(x,r) (h0 + pow(0.19 + ((r)-0.0)/20.0 , -njet)) // widerjet, good for MHD or tilt=90deg
#define thetasam(x,r,w,xp1,xp2) (line1(x,w)*(1.0-trans(x,xp1,xp2)) + line2(x,w)*trans(x,xp1,xp2))

        V[2] = thetasam(X[2],V[1],wparsam(X[2],V[1]),0.25,0.75);
      //      V[2] = thetasam(X[2],V[1],1.0/V[1],0.2,0.8);

      //      dualfprintf(fail_file,"tr=%g %g %g %g\n",tr(0.5),line1(0.5,wparsam(0.5,V[1])),trans(X[2],0.2,0.8),line2(0.5,wparsam(0.5,V[1])));
      

    }

    if(0){ // Sam Gralla 3

      h0=0.0;

#define lineeq(x,w) ((x)*(w)+(0.5*M_PI)-(0.5*w))
#define linepole(x,w) (line1(x,w))
#define thetaL(x,wp,weq,xp1,xp2) ( linepole(x,wp)*(1.0-trans(x,xp1,xp2)) + lineeq(x,weq)*trans(x,xp1,xp2) )
#define thetasam2(x,wp,weq,xp1,xp2) ( x<0.5 ? thetaL(x,wp,weq,xp1,xp2) : -thetaL(1.0-x,wp,weq,xp1,xp2)+M_PI )

      FTYPE poleslope=wparsam(X[2],V[1]);
      FTYPE xpole=0.25;
      //      FTYPE eqslope=0.1;
      FTYPE eqslope=0.5;
      FTYPE xeq=0.5;
      V[2] = thetasam2(X[2], poleslope, eqslope, xpole, xeq);
      

    }





    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }


  if(WHICHUSERCOORD==ORIGWALD){
     V[1] = R0+exp(X[1]) ;
     V[2] = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
     V[3]=2.0*M_PI*X[3];
  }


}


void dxdxp_analytic_user(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  dualfprintf(fail_file,"Should not be computing USERCOORDS analytically\n");
  myexit(34698346);
  dxdxp[3][3] = 2.0*M_PI;    

}
void set_points_user(void)
{
  if(WHICHUSERCOORD==ORIGWALD){
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  }

  if(WHICHUSERCOORD==UNIHALFUSERCOORD){
    startx[1] = 0.3999985081775278946780799743777598329673;
    startx[2] = 0.;
    startx[3] = 0.;
    
    FTYPE endx1=21.40529883372801383045167738115556702610;
    dx[1] = (endx1-startx[1]) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }

  if(WHICHUSERCOORD==JET6LIKEUSERCOORD){
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    if(DOWALDDEN==1) dx[2] = 1. / totalsize[2];
    else if(DOWALDDEN==2) dx[2] = 0.5 / totalsize[2];
    dx[3] = 1.0/totalsize[3];

#if(1)
    startx[1] = log(Rin-R0)/npow;

    trifprintf( "ITERATIVE dx1: Rout=%21.15g R0=%21.15g npow=%21.15g cpow2=%21.15g cpow3=%21.15g npow2=%21.15g x1br=%21.15g rbr=%21.15g\n",Rout,R0,npow,cpow2,cpow3,npow2,x1br,rbr);

    FTYPE x1max0, x1max,dxmax;
    int iter;
    const FTYPE RELACC = NUMEPSILON*100.0;
    const int ITERMAX = 100;


    if( Rout < rbr ) {
      x1max = log(Rout-R0)/npow;
    }
    else {
      x1max0 = 1.1*x1br;
      x1max = 1.2*x1br;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {

        // trifprintf( "iter=%d x1max=%21.15g x2max0=%21.15g\n",iter,x1max0,x1max);

        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;

        if(1){
          dxmax= (pow( (log(Rout-R0) - npow*x1max0)/cpow2, 1./npow2 ) + x1br*1.0) - x1max0;
        }
        else{
          // f-f0 = (x-x0)*dfdx -> if f=Rout -> x = (Rout-f0)/dfdx+x0
          
          FTYPE dVdx1=(npow + cpow2*npow2*cpow3*pow(cpow3*(x1max0-x1br*1.0),npow2-1.0)) * exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          FTYPE V0 = R0 + exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          
          dxmax=(Rout-V0)/dVdx1; // x-x0

          dualfprintf(fail_file,"dVdx1=%g V0=%g dxmax=%g x1max=%g x1max0=%g\n",dVdx1,V0,dxmax,x1max,x1max0);
        }

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;


      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    dx[1] = ( x1max - startx[1] ) /totalsize[1];
#endif


  }

}
FTYPE setRin_user(int ihor, FTYPE ihoradjust)
{
  if(1){
    FTYPE ftemp;

     // see jet3coords_checknew.nb (and fix_3dpolestissue.nb) to have chosen Rin and ihor and compute required R0
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else if(npow2>0){
      return(1.2);
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%21.15g totalsize[1]=%d Rhor=%21.15g R0=%21.15g npow=%21.15g Rout=%21.15g\n",ihoradjust,totalsize[1],Rhor,R0,npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/npow) - ihoradjust*pow(log(Rout-R0),1.0/npow))/(totalsize[1]-ihoradjust),npow)));
    }
  }

  if(WHICHUSERCOORD==ORIGWALD){
    FTYPE ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));

  }
}
 

#define MAXIHOR 10
#define FRACN1 (0.1)
#define ADJUSTFRACT (0.25)

int setihor_user(void)
{
  // set to smaller of either totalsize[1]*0.1 or MAXIHOR
  if(totalsize[1]*FRACN1>MAXIHOR) return((int)MAXIHOR);
  else return((int)((FTYPE)totalsize[1]*(FTYPE)FRACN1));
}
#undef MAXIHOR
#undef FRACN1
#undef ADJUSTFRACT
