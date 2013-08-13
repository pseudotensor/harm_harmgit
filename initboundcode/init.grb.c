
// initial conditions for GRB models




//////////////
//
// General instructions for GRB model:
//
//
// 1) First create EOS table as described in helm/Makefile code
// 2) Copy table to location to be used by "dorungrb.sh"
//    Also ensure stellar models in location used by "dorungrb.sh"
// 3) Ensure lsoffset is set correctly from helm/jon_lsbox.f -> kazfulleos.c (TRUENUCLEAROFFSET and NUCLEAROFFSET)
// 4) cp init.grb.c init.c ; cp init.grb.h init.h ; cp bounds.grb.c bounds.c
// 5) Ensure kazfulleos.c has correctly set TABLE parameters
//    This includes "EOSN?" and "EOSSIMPLEN?" 
// 6) make superclean ; make ; make prep
// 7) 






#include "decs.h"



#define NRADIALMAX 10000
#define MAXPASSPARMS 10

static FTYPE rho[NRADIALMAX],ie[NRADIALMAX],vr[NRADIALMAX],radius[NRADIALMAX],omega3[NRADIALMAX],yl[NRADIALMAX],ynu[NRADIALMAX],ynu0[NRADIALMAX],hcm[NRADIALMAX];
static int NRADIAL;
static int INDEXN;

static FTYPE tauff; 

// changed SLOWFAC, beta,Rbeta, DTd

#define SLOWFAC 0.5  /* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,Rbeta=0.0,rin;
FTYPE rhodisk;

FTYPE Rb,Bpole;


////////////////////////////////////////////////////
//
// which problem to study
#define MASSONE 0 // stick in BH at t=0 like MW99/Proga03/etc.
#define VARYMASS 1 // allow mass and spin to vary (EVOLVEMETRIC must be turned on)
#define UNIFORMDENSITY 2 // Avery's test of uniform collapse
#define GRBPROBLEM 3

//#define WHICHPROBLEM MASSONE
//#define WHICHPROBLEM VARYMASS
//#define WHICHPROBLEM UNIFORMDENSITY
#define WHICHPROBLEM GRBPROBLEM




/////////////////////////////////////////////////////
//
// whether to read initial data and from what file
#define NODATA 0
#define GRBDATA1 1

#if(WHICHPROBLEM==UNIFORMDENSITY)
#define READINITIALDATA NODATA // no choice
#else
#define READINITIALDATA GRBDATA1
//#define READINITIALDATA NODATA
#endif




/////////////////////////////////////////////////////
//
// which magnetic field to use
#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define DIPOLEFIELD 3
#define GRBFIELD1 4
#define GRBFIELD2 5


// only applies for certain WHICHPROBLEMs
//#define FIELDTYPE DISKFIELD
//#define FIELDTYPE DISKVERT
//#define FIELDTYPE DIPOLEFIELD
//#define FIELDTYPE GRBFIELD1
#define FIELDTYPE GRBFIELD2



//////////////////////////////////////////////////
//
// which stellar type (see  init.readdata.c)
// only 0 for now
#define STELLARTYPE 0


/////////////////////////////////////////////////////
//
// which rotation profile in star to use
// the below only apply if READINITIALDATA!=NODATA
#define ROT_MW99 0
#define ROT_PROGA 1
#define ROT_FUJFAST 2
#define ROT_FUJSLOW 3

#define ROTATIONPROFILE ROT_MW99


int fieldfrompotential[NDIM];

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
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  Bpole=1.0;
  Rb=1.0;

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
  extern void get_stellar_data(void);
  FILE *units_file;

  ////////////////////////////
  //
  // GRB units for any WHICHPROBLEM
  //
  ///////////////////////////
  //  mb=mn;
  amu = 1.6605402E-24;
  // definition of baryon mass uses the below
  mb = amu;

  // units
  // We will input in cgs and scale by these units
  // below fixes overall scale of everything

  // fixed and scales r, a, QBH, and MBH such that r will be fixed and others vary
  Mcgs=1.7*msun; // true mass in grams
  Lunit=G*Mcgs/(C*C);
  Tunit=Lunit/C; // so velocity unit is C
  Vunit=Lunit/Tunit; // follows from above 2 things
  Ccode=C/Vunit; // speed of light in code units

  // specify type of rho unit
  // GODMARK: for now rho0 must have same units as Pressureunit
  // This is because often I have rho0+u or something and so must be equal units
  rho0unittype=1;


  if(rho0unittype==0){
    rhounit=1.0; // g/cc
    rhomassunit=rhounit;
    Munit=rhounit*pow(Lunit,3.0); 
    mbwithrhounit=mb/Munit;
    //Munit=Mcgs;
    ////rhounit=Munit/pow(Lunit,3.0);

    // below is only for scale-free simulation
    //  Mdot=0.1*msun;
    //Mdotc=0.35;// approx code units mass accretion rate
    //  rhounit=Mdot/(Mdotc*Lunit*Lunit*C);
    //Munit=rhounit*pow(Lunit,3.0);

    // other conversion factors
    mdotunit=Munit/Tunit;
    energyunit=Munit*Vunit*Vunit;
    mbcsq=mb*Vunit*Vunit/energyunit; // m_b c^2 / Mc^2
    edotunit=energyunit/Tunit;
    Pressureunit=rhounit*Vunit*Vunit;
    Tempunit=Pressureunit*mb/(rhounit*kb);
    Bunit=Vunit*sqrt(rhounit);


    // physics stuff that reqiures Mcgs,mb,edotunit
    ledd=4.*M_PI*C*G*Mcgs*mb/sigmamat;
    leddcode=ledd/edotunit;
  }
  else{
    // choose mass to always be energy unit
    // if so, then need to read-in EOS and convert units correctly

    rhounit=1.0*Vunit*Vunit; // 1g/cc -> 1erg/cc (energy density actualy)
    rhomassunit=rhounit/(Vunit*Vunit); // mass density actually

    Munit=rhounit*pow(Lunit,3.0);  // Mc^2 unit actually
    mbwithrhounit=mb*Vunit*Vunit/Munit; // m_b c^2 / Mc^2
    mdotunit=Munit/Tunit; // \dot{M}c^2 actually
    energyunit=Munit;
    mbcsq=mb*Vunit*Vunit/energyunit; // m_b c^2 / Mc^2
    edotunit=energyunit/Tunit;
    Pressureunit=rhounit;
    Tempunit=Pressureunit*mb*C*C/(rhounit*kb);
    Bunit=sqrt(rhounit);


    // physics stuff that reqiures Mcgs,mb,edotunit
    ledd=4.*M_PI*C*G*Mcgs*mb/sigmamat;
    leddcode=ledd/edotunit;
  }



  trifprintf("Mcgs=%21.15g Lunit=%21.15g Tunit=%21.15g rhounit=%21.15g Munit=%21.15g\n",Mcgs,Lunit,Tunit,rhounit,Munit);
  trifprintf("Vunit=%21.15g mdotunit=%21.15g energyunit=%21.15g Pressureunit=%21.15g Tempunit=%21.15g Bunit=%21.15g\n",Vunit,mdotunit,energyunit,Pressureunit,Tempunit,Bunit);

  


  ////////////////////////
  //
  // UNITS and initial metric parameters FOR WHICHPROBLEM
  //
  ////////////////////////

                                                                                                                     
  ////////////////////////////////////////////////////////////
  //
#if(WHICHPROBLEM==MASSONE)



  if(DOEVOLVEMETRIC==1){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOEVOLVEMETRIC=0\n",WHICHPROBLEM);
    myexit(26629);
  }
  if(DOSELFGRAVVSR==1){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOSELFGRAVVSR=0\n",WHICHPROBLEM);
    myexit(26630);
  }
  if(MCOORD!=KSCOORDS){
    dualfprintf(fail_file,"Expected MCOORD==KSCOORDS\n");
    myexit(26631);
  }

  // units first
  // metric stuff first
  MBH0=1.0;
  a0=0.9375/MBH0;
  QBH0=0.0;

#if(READINITIALDATA==NODATA)
  // can use below when not using stellar model
  Mfactor=1.0;
  Jfactor=1.0;
  rhodisk=1.0;
  rhofactor=1.0;
#else
  if(rho0unittype==0){
    // use units to set parameters
    Mfactor = G*Munit/(C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=Vunit*Vunit;
  }
  else{
    // Munit in mass c^2, so divide by another c^2 for Mfactor
    // use units to set parameters
    Mfactor = G*Munit/(C*C*C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=1.0;
  }
#endif






  ////////////////////////////////////////////////////////////
  //
#elif(WHICHPROBLEM==UNIFORMDENSITY)

  if(DOEVOLVEMETRIC==0){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOEVOLVEMETRIC=1\n",WHICHPROBLEM);
    myexit(26629);
  }
  if(DOSELFGRAVVSR==0){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOSELFGRAVVSR=1\n",WHICHPROBLEM);
    myexit(26630);
  }
  if(!  (MCOORD==KS_BH_TOV_COORDS || MCOORD==KS_TOV_COORDS || MCOORD==BL_TOV_COORDS)  ){
    dualfprintf(fail_file,"Expected MCOORD==KS_BH_TOV_COORDS, KS_TOV_COORDS, or BL_TOV_COORDS\n");
    myexit(26631);
  }


  // units first
  Mfactor=1.0;
  Jfactor=1.0;
  //rhodisk=1E-11; // so not already a black hole!
  //rhodisk=1E-10; // so not already a black hole!
  rhodisk=1E-9; // so forms a black hole!

  // metric stuff first
  Rin=0.0; // chosen here first, repeated below since overwritten
  MBH0=rhodisk*4.0*M_PI/3.0*pow(Rin,3.0); // assumes Rin with rhodisk uniform density
  a0=0.0;
  QBH0=0.0;

  tauff=sqrt(3.0*M_PI/(32.0*rhodisk)); // free-fall time in code units (G=c=1)

  trifprintf(" tauff = %21.15g\n",tauff);




  ////////////////////////////////////////////////////////////
  //
#elif(WHICHPROBLEM==VARYMASS)

  if(DOEVOLVEMETRIC==0){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOEVOLVEMETRIC=1\n",WHICHPROBLEM);
    myexit(26629);
  }
  if(DOSELFGRAVVSR==1){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOSELFGRAVVSR=0\n",WHICHPROBLEM);
    myexit(26630);
  }
  if(!  (MCOORD==KSCOORDS || MCOORD==BLCOORDS || MCOORD==KS_BH_TOV_COORDS || MCOORD==KS_TOV_COORDS || MCOORD==BL_TOV_COORDS)  ){
    dualfprintf(fail_file,"Expected MCOORD==KSCOORDS, BLCOORDS, KS_BH_TOV_COORDS, KS_TOV_COORDS, or BL_TOV_COORDS\n");
    myexit(26631);
  }

  if(rho0unittype==0){
    // use units to set parameters
    Mfactor = G*Munit/(C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=Vunit*Vunit;
  }
  else{
    // Munit in mass c^2, so divide by another c^2 for Mfactor
    // use units to set parameters
    Mfactor = G*Munit/(C*C*C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=1.0;
  }

#if(READINITIALDATA==NODATA)
  // below 2 come from stellar model and what exists within Lunit at t=0
  // MBH is in # of Mcgs since how used in metric.
  // That is, in metric, r= radius per unit Lunit and MBH is per unit Lunit assuming G=c=1
  MBH0=(0.001*msun)/(Lunit*C*C/G); // used in metric
  a0=0.1*MBH0; // j=0.1
  QBH0=0.0;


  // Below for normalizing initial conditions
  rhodisk=1E10*rhounit; // disk peak density.  Will scale internal energy density and field components (though beta) too.
  // velocity already scaled by speed of light
#elif(READINITIALDATA==GRBDATA1)

  // then assume will be read in since data file has MBH0,a0,QBH0

#endif






  ////////////////////////////////////////////////////////////
  //
#elif(WHICHPROBLEM==GRBPROBLEM)

  if(DOEVOLVEMETRIC==0){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOEVOLVEMETRIC=1\n",WHICHPROBLEM);
    myexit(26629);
  }
  if(DOSELFGRAVVSR==0){
    dualfprintf(fail_file,"WHICHPROBLEM=%d uses DOSELFGRAVVSR=1\n",WHICHPROBLEM);
    myexit(26630);
  }
  if(!  (MCOORD==KSCOORDS || MCOORD==BLCOORDS || MCOORD==KS_BH_TOV_COORDS || MCOORD==KS_TOV_COORDS || MCOORD==BL_TOV_COORDS)  ){
    dualfprintf(fail_file,"Expected MCOORD==KSCOORDS, BLCOORDS, KS_BH_TOV_COORDS, KS_TOV_COORDS, or BL_TOV_COORDS\n");
    myexit(26631);
  }


  if(rho0unittype==0){
    // use units to set parameters
    Mfactor = G*Munit/(C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=Vunit*Vunit;
  }
  else{
    // Munit in mass c^2, so divide by another c^2 for Mfactor
    // use units to set parameters
    Mfactor = G*Munit/(C*C*C*C*Lunit); // converts mass in code units to mass in metric length units (for evolution of BH mass)
    Jfactor = G*Munit/(C*C*C*C*Lunit);
    rhodisk=1.0; // just for initial conditions and RHOMIN/UUMIN
    // below converts mass density to energy density
    rhofactor=1.0;
  }





#endif



  trifprintf("Mfactor=%21.15g Jfactor=%21.15g rhodisk=%21.15g\n",Mfactor,Jfactor,rhodisk);









  if(READINITIALDATA==GRBDATA1){

    if(DOEVOLVERHO==0||DOEVOLVEUU==0||DOEVOLVEYL==0||DOEVOLVEYNU==0){
      dualfprintf(fail_file,"WARNING: Should be GRMHD model, normally with evolution of Y_l and Y_\\nu.\n");
    }
  }





  ///////////////////////////////////////////////
  //
  // READ DATA FROM FILE (or not)
  //
  ///////////////////////////////////////////////

#if(READINITIALDATA==NODATA)
  // the nothing to read in
#elif(READINITIALDATA==GRBDATA1)
  // read stellar model data into predefined arrays
  get_stellar_data();
#endif


  ///////////////////////////////////////////////
  //
  // Override metric parameters read from data file
  //
  ///////////////////////////////////////////////

#if(WHICHPROBLEM==MASSONE)
  // GODMARK: override stellar model's BH mass
  MBH0=1.0;
  a0=0.9375/MBH0;
  QBH0=0.0;
#endif



  ///////////////////////////////////////////////
  //
  // READ EOS from data file
  //
  ///////////////////////////////////////////////

#if(WHICHEOS==KAZFULL && EOMTYPE!=EOMCOLDGRMHD)
  read_setup_eostable();
#endif

  ///////////////////////////////////////////////
  //
  // Output units to file
  //
  ///////////////////////////////////////////////

  units_file=fopen("0_units.dat","wt");
  if(units_file==NULL){
    dualfprintf(fail_file,"Cannot open 0_units.dat\n");
    myexit(2863);
  }

  fprintf(units_file,"%d\n",rho0unittype);
  fprintf(units_file,"%21.15g %21.15g %21.15g\n",Munit,Lunit,Tunit);
  fprintf(units_file,"%21.15g %21.15g %21.15g %21.15g\n",Mfactor,Jfactor,Mcgs,rhodisk);
  fprintf(units_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",rhounit,rhomassunit,Vunit,mdotunit,energyunit,Pressureunit,Tempunit,Bunit);
  fclose(units_file);

  return(0);

}




int init_defcoord(void)
{

#if(WHICHPROBLEM==VARYMASS) // new method
  //////////////////////////////////
  //  new M anything G=c=1 method (could be normal torus model)
  // originally envisioned without self-gravity but evolution of metric
  // Uses relatively large Rin so that material just adds up to BH inside
  // Original idea was that Rin must be large to avoid horizon (and rminus inside horizon) reaching grid
  // Now if horizon passes Rin we will automagically move the evolved part of the grid outward (benefit of self-gravity tests where black hole grows substantially)
  // So this method is somewhat outdated, but still useable to avoid evolving inside Rin
  //

  // define coordinate type
  defcoord = JET4COORDS;
  




#elif(WHICHPROBLEM==GRBPROBLEM)
  //////////////////////////////////
  //  new M anything G=c=1 method with full GRB model with self-gravity in mind
  // starts with Rin=0 and MBH=0


#if(0)

  // define coordinate type
  //  defcoord = JET4COORDS;
  defcoord=LOGRSINTH; // best if following spherical collapse
#elif(1)

  // define coordinate type
  //  defcoord = JET4COORDS;
  //  defcoord=LOGRSINTH; // best if following spherical collapse

  defcoord=UNI2LOG; // Afactor, Nstar, and Rstar set in coord.c
  //defcoord = UNIRSINTH;
  //defcoord = UNIRSINTH2; // try with dxdxp11!=1.0
#else
  defcoord = UNIRSINTH;

#endif  







#elif(WHICHPROBLEM==MASSONE)
  /////////////////////////////
  /////////// old GM=c=1 method
  

  // define coordinate type
  defcoord = JET4COORDS;








#elif(WHICHPROBLEM==UNIFORMDENSITY)
  /////////////////////////////
  /////////// Studying self-gravitational collapse with formation of black hole
  //  defcoord = 0;
  defcoord = UNIRSINTH;


#endif

  return(0);
}





int init_grid(void)
{
  SFTYPE rh;



  a = a0 ;
  MBH=MBH0;
  QBH=QBH0;
  


#if(WHICHPROBLEM==VARYMASS) // new method
  //////////////////////////////////
  //  new M anything G=c=1 method (could be normal torus model)
  // originally envisioned without self-gravity but evolution of metric
  // Uses relatively large Rin so that material just adds up to BH inside
  // Original idea was that Rin must be large to avoid horizon (and rminus inside horizon) reaching grid
  // Now if horizon passes Rin we will automagically move the evolved part of the grid outward (benefit of self-gravity tests where black hole grows substantially)
  // So this method is somewhat outdated, but still useable to avoid evolving inside Rin
  //

  Rin=(0.95*2.0*Lunit)/Lunit; // fixed in time
  Rout=(2E3*Lunit)/Lunit; // fixed in time

  R0=-2.0;
  Rhor=rhor_calc(0);
  hslope = 0.3;

  




#elif(WHICHPROBLEM==GRBPROBLEM)
  //////////////////////////////////
  //  new M anything G=c=1 method with full GRB model with self-gravity in mind
  // starts with Rin=0 and MBH=0


#if(0)
  Rin=0.0; // at least initially
  Rout=(2E3*Lunit)/Lunit; // fixed in time
  R0=-2.0;
  Rhor=rhor_calc(0);
  hslope = 0.3;

#elif(1)
  Rin=0.0; // at least initially
  Rout=(1E4*Lunit)/Lunit; // fixed in time
  R0=-2.0;
  Rhor=rhor_calc(0);
  hslope = 0.3;

#else
  R0 = 0.0;
  Rout = 2E3;
  Rin = 0.0;
  hslope = 1.0;

#endif  







#elif(WHICHPROBLEM==MASSONE)
  /////////////////////////////
  /////////// old GM=c=1 method
  


  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
#if(0)
  R0 = -3.0;
  Rhor=rhor_calc(0);
  Rout = 1E5;
#endif
#if(1)
  R0 = -2.0;
  Rhor=rhor_calc(0);
  Rout = 3E3;
#endif

  //  setRin_withchecks(&Rin);
  Rin = 0.95 * Rhor;


  hslope = 0.3;








#elif(WHICHPROBLEM==UNIFORMDENSITY)
  /////////////////////////////
  /////////// Studying self-gravitational collapse with formation of black hole
  


  R0 = 0.0;
  Rout = 2E3;
  //  Rin = 1.0;
  Rin = 0.0;

  hslope = 1.0;



#endif
 
  

  trifprintf("R0=%21.15g Rin=%21.15g Rout=%21.15g hslope=%21.15g defcoord=%d\n",R0,Rin,Rout,hslope,defcoord);





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


  /* some physics parameters */
#if(WHICHEOS==IDEALGAS)
  //  gam = 4. / 3.; // for stellar model leads to early bounce
  //  gam = 1.1; // for stellar model leads to formation of 0.1M black hole
  gam = 1.2;
#else
  gam = -1.0; // <0 indicates not ideal gas EOS
#endif  
  gamideal=4.0/3.0; // ideal gas EOS index to use if off table

  // below assumes tabulated EOS has minimum internal energy of zero, which is approximately correct currently.
  zerouuperbaryon=-TRUENUCLEAROFFSETPRIMARY; // definition of zero-point energy per baryon (i.e. such that internal energy cannot be lower than this times baryon density)


  //  cooling=NOCOOLING;
  cooling=COOLEOSGENERAL; // neutrino cooling via tabulated Kaz EOS

  trifprintf("Cooling is set to: %d\n",cooling);
  if(cooling!=2) dualfprintf(fail_file,"SUPERWARNING: Neutrino cooling not enabled.\n");



  //  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  //  FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;


  /////////////////////////////////////////
  //
  // BOUNDARY CONDITIONS

  // R0SING is used when wanting to follow full collapse and may not have any black hole at t=0
  // then inside horizon (once formed) outflow will occur as in bounds.c and metric-extrapolation in metric_selfgravity_or_evolvemetric.c

#if(WHICHPROBLEM==UNIFORMDENSITY)
  //  BCtype[X1UP]=FREEOUTFLOW;
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=R0SING;
#elif(WHICHPROBLEM==GRBPROBLEM)
  //  BCtype[X1UP]=FREEOUTFLOW;
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=R0SING;
#else
  BCtype[X1UP]=FIXEDOUTFLOW;
  //BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
#endif

  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;


  ////////////////////
  //
  // ATMOSPHERE/FLOOR
  //
  ///////////////////

  // new type of density/internal energy injection model that limits b^2/rho and b^2/u
  // see collapsar_environment_density.nb
  rescaletype=4;
  BSQORHOLIMIT=1E2;
  BSQOULIMIT=1E3;
  UORHOLIMIT=1E3;

  // not dimensionless
  RHOMIN = 1E-3*rhodisk;
  UUMIN = 1E-4*rhodisk;

  trifprintf("RHOMIN=%g UUMIN=%g\n",RHOMIN,UUMIN);


  /////////////////////////////////////
  //
  // time in fixed units of GMcgs/c^3
  //
  /////////////////////////////////////

  /* output choices */
#if(WHICHPROBLEM==UNIFORMDENSITY)
  tf = tauff*.986153; // ZEUS time whereby density is 3 orders of magnitude higher
  //tf = tauff*.81831; // Broderick & Rathore (2006)


  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = tf/100.0;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= DTdumpgen[MAINDUMPTYPE];
  // ener period
  DTdumpgen[ENERDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = DTdumpgen[MAINDUMPTYPE];

  DTr = 100;   /* restart file period in steps */



#else


#if(0)
  tf = 5E5;

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = tf/1000.0;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= DTdumpgen[MAINDUMPTYPE];
  // ener period
  DTdumpgen[ENERDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/100.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = DTdumpgen[MAINDUMPTYPE];


  DTr = 100;   /* restart file period in steps */
#else
  // DEBUG:
  //  tf = 6483.63983628599; // dump~766
  //  tf = 4450.31668938777; // dump=100
  //tf = 4403.48368138165;
  tf = 1E5;


  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = tf/100.0;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= DTdumpgen[MAINDUMPTYPE];
  // ener period
  DTdumpgen[ENERDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/100.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/10.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = DTdumpgen[MAINDUMPTYPE];

  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;   /* restart file period in steps */
#endif



#endif



  trifprintf("tf=%g DTd=%g DTavg=%g DTener=%g DTi=%g DTdebug=%g DTr=%d\n",tf,DTdumpgen[MAINDUMPTYPE],DTdumpgen[AVG1DUMPTYPE],DTdumpgen[ENERDUMPTYPE],DTdumpgen[IMAGEDUMPTYPE],DTdumpgen[DEBUGDUMPTYPE],DTr);

  return(0);

}




// assumes normalized density
// only used if not reading in stellar model
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom realgeomdontuse;
  struct of_geom *ptrrealgeom=&realgeomdontuse;
  FTYPE pratm[NPR];


#if(READINITIALDATA==GRBDATA1)
  return(0);
#endif


  get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,ptrrealgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  //  if(pr[RHO]<pratm[RHO]){
  if(pr[RHO]<SMALL){
    PLOOP(pliter,pl) pr[pl]=pratm[pl];
  }
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);

}



// do after grid is set
// checks location of horizon
int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  extern void check_spc_singularities_user(void);

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);



  //SASMARK restart: need to populate panalytic with IC's
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
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







////////////////////////////
//
// set primitives at t=0
//
///////////////////////////
int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;
  int getmax_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax);

  funreturn=user1_init_primitives(prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

#if(READINITIALDATA==NODATA)

  // standard user1_init_primitives() result

#elif(READINITIALDATA==GRBDATA1)

  // Make sure normalize_densities() and init_atmopshere() do nothing in this case

  // then assume density already normalized properly and that atmosphere is how wanted
  // still get rhomax and umax
  getmax_densities(prim,&rhomax, &umax);

#endif


  return(0);


}




// set unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_disk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_star(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_uniformdensity(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);



#if(READINITIALDATA==NODATA)

#if(WHICHPROBLEM!=UNIFORMDENSITY)
  init_disk(whichvel,whichcoord,i,j,k,pr,pstag);
#else
  init_uniformdensity(whichvel,whichcoord,i,j,k,pr,pstag);
#endif


#elif(READINITIALDATA==GRBDATA1)

  init_star(whichvel,whichcoord,i,j,k,pr,pstag);

#endif


  return(0);

}



// set solution to be a disk with atmosphere
int init_disk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
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
  FTYPE jbh;
  int pl,pliter;


  if(EOMTYPE!=EOMGRMHD){
    dualfprintf(fail_file,"Should be GRMHD model\n");
    myexit(7724);
  }




  // rin : inner radius of torus
  // rmax : location of pressure maximum
  // notice that outer radius of torus is defined through torus needing to be an equilibrium

  if(WHICHPROBLEM==MASSONE){
    rin = 6.0 ;
    rmax = 12.0 ;
  }
  else if(WHICHPROBLEM==VARYMASS || WHICHPROBLEM==GRBPROBLEM){
    rin = 6.0/MBH; // notice that this is not at same radius away from BH as above
    rmax = 12/MBH;
  }

  l = lfish_calc(rmax) ;
  kappa = 1.e-3 ;
  beta = 100.0 ;
  Rbeta=0.0; // means not normalized at one location, but pgmax/pbmax
  // beta = 1E30 ;
  randfact = 4.e-2;
  jbh=a/MBH; // normalize so can use formulas without MBH
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);

  r=V[1]/MBH; // normalize so can use formulas without MBH
  th=V[2];



  sth = sin(th);
  cth = cos(th);

  /* calculate lnh  (r, jbh, etc. all dimensionless) */
  DD = r * r - 2. * r + jbh * jbh;
  AA = (r * r + jbh * jbh) * (r * r + jbh * jbh) - DD * jbh * jbh * sth * sth;
  SS = r * r + jbh * jbh * cth * cth;
  
  thin = M_PI / 2.;
  sthin = sin(thin);
  cthin = cos(thin);
  DDin = rin * rin - 2. * rin + jbh * jbh;
  AAin = (rin * rin + jbh * jbh) * (rin * rin + jbh * jbh)
    - DDin * jbh * jbh * sthin * sthin;
  SSin = rin * rin + jbh * jbh * cthin * cthin;
  
  if (r >= rin) {
    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
                               (AA * sth * AA * sth))) / (SS * DD /
                                                          AA))
      - 0.5 * sqrt(1. +
                   4. * (l * l * SS * SS) * DD / (AA * AA * sth *
                                                  sth))
      - 2. * jbh * r * l / AA -
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
       - 2. * jbh * rin * l / AAin);
  } else
    lnh = 1.;
  



  // just define some field
  pr[B1]=0.0;
  pr[B2]=0.0;
  pr[B3]=0.0;
  if(FLUXB==FLUXCTSTAG){
    // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
    PLOOPBONLY(pl) pstag[pl]=pr[pl];
  }


  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {
    

    //    get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
    //set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    pr[RHO]=0.0;
    pr[UU]=0.0;

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
    up = 2. * jbh * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    // put back in MBH since length and time units above was without MBH
    pr[RHO] *=rhodisk; // arbitrary density
    pr[UU]  *=rhodisk; // arbitrary density
    pr[U1]  *=1.0; // dr/dt in BL-coords
    pr[U2]  /=MBH; // d\theta/dt    
    pr[U3]  /=MBH; // d\phi/dt


    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}



// set solution to be uniform density (to follow collapse and BH formation)
int init_uniformdensity(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeom,geom;
  int pl,pliter;

  if(EOMTYPE!=EOMCOLDGRMHD){
    dualfprintf(fail_file,"Should be p=0 model\n");
    myexit(7723);
  }


  beta = 1E30;
  Rbeta=0.0; // means not normalized at one location, but pgmax/pbmax


  coord(i, j, k, CENT, X);
  bl_coord(X, V);

  r=V[1];
  th=V[2];

  // so that within volume GM/c^2 not too small (need to resolve gravity relativistically)
  // density is set 1 here and normalized later to be rhodisk
  pr[RHO]=1.0; 
  pr[UU]=0.0;

  // just define some field
  pr[B1]=0.0;
  pr[B2]=0.0;
  pr[B3]=0.0;
  if(FLUXB==FLUXCTSTAG){
    // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
    PLOOPBONLY(pl) pstag[pl]=pr[pl];
  }


  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}



// set vector potential
// assumes normal field in pr
// SUPERNOTE: A_i must be computed consistently across all CPUs.  So, for example, cannot use randomization of vector potential here.
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE r,th;
  FTYPE vpot;
  FTYPE Rtrans;




  vpot=0;

  if(l==3){// A_\phi

    r=V[1];
    th=V[2];


    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      vpot += 0.5*pow(r,3.0/4.0)*sin(th)*sin(th) ;
    }
    if(FIELDTYPE==GRBFIELD1){


      Rtrans=14341.8; // Nagataki 2007
      if(r>Rtrans){
        // dipole
        vpot += 0.5*Bpole*Rb*Rb*sin(th)*sin(th)/((r+SMALL)/Rtrans) ;
      }
      //      else if(r>
      // uniform field
      // vpot += 0.5*Bpole*Rb*Rb*(Rb/Rtrans)*sin(th)*sin(th)*(r/Rtrans)*(r/Rtrans) ;
      // uniform field makes no sense

      //}
      else{
        // split-monopole makes sense
        if(th<M_PI*0.5){
          vpot += -0.5*Bpole*Rb*Rb*cos(th);
        }
        else{
          vpot += 0.5*Bpole*Rb*Rb*cos(th);
        }
      }
    }


    if(FIELDTYPE==GRBFIELD2){
      // dipole with finite interior


      Rtrans=3.6E9/Lunit; // Nagataki 2007
      if(r>Rtrans){
        // dipole
        vpot += 0.5*Bpole*Rb*Rb*sin(th)*sin(th)/((r+SMALL)/Rtrans) ;
      }
      else{
        // uniform field
        vpot += 0.5*Bpole*Rb*Rb*sin(th)*sin(th)*(r/Rtrans)*(r/Rtrans) ;
      }
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
// different than getmax_densities() in that rin is used to limit u for beta and field
// also normalize density scale so density is unity
int normalize_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  FTYPE parms[MAXPASSPARMS];
  int eqline;

#if(READINITIALDATA==GRBDATA1)
  return(0);
#endif

  eqline=1;
  parms[0]=rin*MBH;
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
  
  eqslice=1;
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


#undef SLOWFAC


SFTYPE lfish_calc(SFTYPE oldr)
{
  FTYPE newr,jbh;
  FTYPE ud3;

  newr=oldr;// already normalized
  jbh=a/MBH;

  ud3= (((pow(jbh, 2) - 2. * jbh * sqrt(newr) + pow(newr, 2)) *
         ((-2. * jbh * newr * (pow(jbh, 2) - 2. * jbh * sqrt(newr) + pow(newr, 2))) /
          sqrt(2. * jbh * sqrt(newr) + (-3. + newr) * newr) +
          ((jbh + (-2. + newr) * sqrt(newr)) * (pow(newr, 3) +
                                                pow(jbh,
                                                    2) * (2. + newr))) / sqrt(1 +
                                                                              (2.
                                                                               *
                                                                               jbh)
                                                                              /
                                                                              pow
                                                                              (newr,
                                                                               1.5)
                                                                              -
                                                                              3.
                                                                              /
                                                                              newr)))
        / (pow(newr, 3) * sqrt(2. * jbh * sqrt(newr) + (-3. + newr) * newr) *
           (pow(jbh, 2) + (-2. + newr) * newr))
        );

  //  return(ud3*MBH); // puts "units" back in
  return(ud3); // assumes still normalized by MBH

}



int set_atmosphere(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  int set_atmosphere_grb1(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
  int set_atmosphere_grb2(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr);


#if(READINITIALDATA==NODATA)
  return(set_atmosphere_grb1(whichcond, whichvel, ptrgeom, pr));
#elif(READINITIALDATA==GRBDATA1)
  return(set_atmosphere_grb2(whichcond, whichvel, ptrgeom, pr));
#endif


}



// This shouldn't be called at t=0 until after panalytic is set
// assumes atmosphere or any kind of "FIXEDINFLOW" called by inflow_check or other functions uses original solution
int set_atmosphere_grb2(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  int ii,jj,kk;
  int pl,pliter;



  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;

  PLOOP(pliter,pl){
    pr[pl]=GLOBALMACP0A1(panalytic,ii,jj,kk,pl);
  }

  return(0);


}





// old function to set angular velocity of stellar material
int set_phi_velocity_grb1(FTYPE *X, FTYPE *V, FTYPE dxdxp[NDIM][NDIM], struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE omegak,omega,lspecific;
  FTYPE r,th;



  r=V[1];
  th=V[2];



  // HACK 
  if(r>5.0){
    // then setup rotating atmosphere with constant specific angular momentum (on cylinders) 
    // notice that the 4.0 is in terms of G=c=1 
    pr[U3] += 4.0*dxdxp[3][3]*ptrgeom->gcon[GIND(3,3)]*sin(th)*sin(th);
  }

#if(0)
  // setup as in Proga et al. (2003) 
  if(r<2.0){
    // else not rotation anymore than from ZAMO spin 
  }
  else if(r<840.59){ // Helium-Hydrogen transition 
    omegak=1.0/(pow(r,3.0/2.0)); // 1.0 here is in terms of Mcgs 
    omega=sqrt(0.02)*omegak; // 0.02 = centrifugal/(R-cylindrical gravity) 
    lspecific=pow(r*sin(th),2.0)*omega; // u_\phi 
    if(lspecific>13.29) lspecific=13.29;
    pr[U3] += lspecific*dxdxp[3][3]*ptrgeom->gcon[GIND(3,3)]; // u^\phi 
  }
  // else not rotating 
  //} 
#endif

  return(0);

}




/////////////////////////////////
//
// old function to set atmosphere (before stellar model)

// UUMIN/RHOMIN used for atmosphere
// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int set_atmosphere_grb1(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int pl,pliter;
  int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
  int set_phi_velocity_grb1(FTYPE *X, FTYPE *V, FTYPE dxdxp[NDIM][NDIM], struct of_geom *ptrgeom, FTYPE *pr);




  // Bondi like initial atmosphere
  //    rho = RHOMIN * 1.E-14;
  //    u = UUMIN * 1.E-14;
  coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
  bl_coord_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, V);
  dxdxprim_ijk(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, dxdxp);

  r=V[1];
  th=V[2];

  // default
  PLOOP(pliter,pl) prlocal[pl]=pr[pl];


  /////////////////////////////////////////////////////////////
  //
  // SET DENSITIES


#if((EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD))
  // Bondi-like atmosphere
  if(rescaletype==4){
    // couple rescaletype to atmosphere type
    if(r>40.0) prlocal[RHO] = RHOMIN*pow(r,-2.0);
    else prlocal[RHO] = RHOMIN*pow(40.0,-2.0);
    //prlocal[RHO] = RHOMIN*pow(r,-2.0);
  }
  else{
    prlocal[RHO] = RHOMIN*pow(r,-1.5);
  }
#else
  prlocal[RHO] = 0;
#endif

#if(EOMTYPE==EOMGRMHD)
  // Bondi-like atmosphere
  if(rescaletype==4){
    // couple rescaletype to atmosphere type
    if(r>40.0) prlocal[UU]  = UUMIN*pow(r,-2.5);
    else prlocal[UU]  = UUMIN*pow(40.0,-2.5);
    //prlocal[UU]  = UUMIN*pow(r,-2.5);
  }
  else{
    prlocal[UU]  = UUMIN*pow(r,-2.5);
  }
#else
  prlocal[UU]  = 0;
#endif


  /////////////////////////////////////////////////////////////
  //
  // SET Velocity



  set_zamo_velocity(whichvel, ptrgeom,prlocal);
  set_phi_velocity_grb1(X,V,dxdxp,ptrgeom,prlocal); // adds to zamo


  /////////////////////////////////////////////////////////////
  //
  // Set depending upon whichcond

  
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
  int trueglobalregion_enerdef[NUMUPDOWN][NDIM];



  // GODMARK:
  // TODO:
  // 1) With ACTIVEREGION used, need to include Mvsr(r), etc. within grid but not inside black hole.  How to grow black hole then?  Don't, but need to add to Mvsr(r) at least.  See SECTIONMARK in metric_selfgravity_or_evolvemetric.c
  // 2) get_new_metric_parms() accretion to BH only, need to accrete to normal non-BH mass if inner radius is far beyond horizon radius

  // if(Rin>rhor) then grow Mvsr[Rin] to model growth beyond domain.  if(Rin<rhor) then grow MBH.
  // But Mvsr integral can only be including r>rhor, not just i>=iactiveregion. Or ignore Rin vs. Rhor difference for accretion rate?  What did before, just do again.

  // apparently can use enerregiondef[POINTDOWN][1] for A=OUTSIDEHORIZONENERREGION and B=ACTIVEREGION.  If A>B, then add mass to MBH, else add mass to Mvsr(r<

  // GODMARK: Should also follow check like in Sasha's version of enerregiondef[][]

#if(1)
  // move active region with black hole horizon (but inside by some number of cells: MAXBND for now)
  enerregiondef[POINTDOWN][1]=MAX(horizoni+horizoncpupos1*N1-MAXBND,0);
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
#else
  // no limiting of loops to be beyond black hole
  enerregiondef[POINTDOWN][1]=0;
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
#endif


  
  // now GRIDSECTION can define ACTIVEREGION




  return(0);
}


int theproblem_set_enerregionupdate(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps)
{

  ////////////////////
  //
  // Setup update period
  //
  ///////////////////
  *updateeverynumsteps=1; // SUPERGODMARK
  //number of steps after which position/size of active section is updated
  *everynumsteps = *updateeverynumsteps;

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



// GRB read of initial data
#include "init.readdata.c"
