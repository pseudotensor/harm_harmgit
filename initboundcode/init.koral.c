
#include "decs.h"


// NOTE on units:
// Many things below are so far in code units, not physical units, so they don't need conversion.  This includes:
//
// Rin, Rout, tf, DTdumpgen


int BEAMNO,FLATBACKGROUND; // global for bounds.koral.c



FTYPE normglobal;

int prepre_init_specific_init(void)
{
  int funreturn;


  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  periodicx1=periodicx2=periodicx3=1;

  // Also: SET USEROMIO to 0 or 1 in mympi.definit.h (needs to be 0 for TEXTOUTPUT)
  binaryoutput=TEXTOUTPUT;

  return(0);

}


int pre_init_specific_init(void)
{

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

  // print out units and some constants
  trifprintf("Constants\n");
  trifprintf("LBAR=%g TBAR=%g VBAR=%g RHOBAR=%g MBAR=%g UBAR=%g TEMPBAR=%g\n",LBAR,TBAR,VBAR,RHOBAR,MBAR,UBAR,TEMPBAR); 
  trifprintf("ARAD_CODE=%26.20g OPACITYBAR=%g KAPPA_ES_CODE(1,1)=%g KAPPA_FF_CODE(1,1)=%g KAPPA_BF_CODE(1,1)=%g\n",ARAD_CODE,OPACITYBAR,KAPPA_ES_CODE(1,1),KAPPA_FF_CODE(1,1),KAPPA_BF_CODE(1,1));
  trifprintf("ARAD_CODE_DEF=%g\n",ARAD_CODE_DEF);


  //  cour=0.8;
  

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}





int init_global(void)
{
  int pl,pliter;
  int funreturn;

  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);


  // default
  ARAD_CODE=ARAD_CODE_DEF;


  //  ERADLIMIT=UUMIN; // set same for now
  ERADLIMIT=UUMINLIMIT; // seems fine.
  // maximum radiation frame lorentz factor
  GAMMAMAXRAD=10000.0;
  //  GAMMAMAXRAD=100.0;



  //////////////////
  // overrides for more detailed problem dependence
  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  //  FLUXB=FLUXCTTOTH;
  FLUXB=FLUXCTSTAG;
  

  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=1E3;
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
  //OSMARK: where is DTr1 defined? what is DTfake?
  DTfake=MAX(1,DTr/10); 




  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==FLATNESS){
	
	//  lim[1]=lim[2]=lim[3]=MINM;
	cour=0.8;
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

    lim[1]=lim[2]=lim[3]=MINM;
	//  cour=0.5;
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
      //      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1; // testing
    }
    else if(WHICHPROBLEM==RADPULSE){
      tf = 35; //final time
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    }
	else if(WHICHPROBLEM==RADPULSE3D){
      tf = 70; //final time
      for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;
    }
  }

  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADBEAMFLAT){
	cour=0.8;
	gam=gamideal=5.0/3.0;
	cooling=KORAL;

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
#define NTUBE 31
//#define NTUBE 5
//#define NTUBE 3

	lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
	// should have PARA(LINE) not oscillate so much at cusp
	// Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
	cour=0.8;
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

	lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
	// should have PARA(LINE) not oscillate so much at cusp
	// Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
	cour=0.8;
	gam=gamideal=1.4;
	cooling=KORAL;
	ARAD_CODE=1E7*1E-5*(2.5E-9/7.115025791e-10); // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

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

	lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
	//	lim[1]=lim[2]=lim[3]=MC; // MC gets totally bonkers answer with NLEFT=0.99999
	//lim[1]=lim[2]=lim[3]=PARALINE; // bonkers answer for NLEFT=0.99999, ok for NLEFT=0.99
	// should have PARA(LINE) not oscillate so much at cusp
	// Also should eliminate PARA's zig-zag steps in internal energy density in other tests.
	cour=0.8;
	//	cour=0.49; // doesn't help oscillations for NLEFT=0.99999 with MINM
	gam=gamideal=1.4;
	cooling=KORAL;
	ARAD_CODE=1E-30;
	//ARAD_CODE=1E7*1E-5*(2.5E-9/7.115025791e-10); // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

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
	tf = 200.0; //final time (far past plot so see evolves or stationary).
  }


  /*************************************************/
  /*************************************************/
  /*************************************************/

  if(WHICHPROBLEM==RADBEAM2D){

    BEAMNO=1; // 1-4
// whether constant or radially varying background
// ==0 doesn't make much sense for Minkowski without gravity, because flow reverses due to chosen high density
    FLATBACKGROUND=1;


	lim[1]=lim[2]=lim[3]=MINM; // NTUBE=1 has issues near cusp, so use MINM
	a=0.0; // no spin in case use MCOORD=KSCOORDS

	if(!(ISSPCMCOORDNATIVE(MCOORD))){
	  dualfprintf(fail_file,"Must choose MCOORD (currently %d) to be spherical polar grid type for RADBEAM2D\n",MCOORD);
	  myexit(3434628752);
	}

	cour=0.8;
	gam=gamideal=1.4;
	cooling=KORAL;
	//	ARAD_CODE=ARAD_CODE_DEF*1E5; // tuned so radiation energy flux puts in something much higher than ambient, while initial ambient radiation energy density lower than ambient gas internal energy.

	BCtype[X1UP]=RADBEAM2DFLOWINFLOW;
    BCtype[X1DN]=OUTFLOW;
    //	BCtype[X1DN]=HORIZONOUTFLOW;
	BCtype[X2UP]=PERIODIC;
	BCtype[X2DN]=PERIODIC;
	//	BCtype[X3UP]=FREEOUTFLOW;
	BCtype[X3UP]=OUTFLOW;
	BCtype[X3DN]=RADBEAM2DBEAMINFLOW;


	FTYPE DTOUT1;
	if (BEAMNO==1){
	  DTOUT1=.1; //dt for basic output
	}
	else if (BEAMNO==2){
	  DTOUT1=.4; //dt for basic output
	}
	else if (BEAMNO==3){
	  DTOUT1=1.; //dt for basic output
	}
	else if (BEAMNO==4){
	  DTOUT1=.25; //dt for basic output
	}

	int idt;
	for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=DTOUT1;

	DTr = 100; //number of time steps for restart dumps
	tf = 20.0; //final time

    //    DODIAGEVERYSUBSTEP = 1;

  }



  /*************************************************/
  /*************************************************/
  /*************************************************/

  return(0);

}


int init_defcoord(void)
{


  /*************************************************/
  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==FLATNESS){


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
  if(WHICHPROBLEM==RADBEAM2D){


	defcoord = UNIFORMCOORDS;
	Rin_array[1]=0;
	Rin_array[2]=0;
	Rin_array[3]=0;

	Rout_array[1]=1.0;
	Rout_array[2]=1.0;
	Rout_array[3]=1.0;
  


	if (BEAMNO==1){
	  Rin_array[1]=2.6;
	  Rout_array[1]=3.5;
	}
	else if (BEAMNO==2){
	  Rin_array[1]=5.5;
	  Rout_array[1]=7.5;
	}
	else if (BEAMNO==3){
	  Rin_array[1]=14.5;
	  Rout_array[1]=20.5;
	}
	else if (BEAMNO==4){
	  Rin_array[1]=30;
	  Rout_array[1]=50;
	}

	Rin_array[2]=0.99*M_PI*0.5;
	Rout_array[2]=1.01*M_PI*0.5;

	Rin_array[3]=0.0;
	Rout_array[3]=M_PI*0.25;

	
  }


  return(0);
}


int init_grid(void)
{
  

  return(0);
}



// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int funreturn;

  // NO atmosphere
  //  funreturn=user1_init_atmosphere(whichvel, whichcoord,i, j, k, pr);
  //  if(funreturn!=0) return(funreturn);
  //  return(0);

  // tells to do no coordinate transformations
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  //  return(0);

  return(-1); // no atmosphere set, so don't do anything at all

}

int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);



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
  dualfprintf(fail_file,"WARNING: done with check_spc_singularities_user(), but it sometimes stalls or goes very very slow for no apparently good reason.  E.g., on NAUTILUS with -O0, very slow checks.  But just putting dualfprintf before and after the above call leads to quick finish.\n");

  
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
#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


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
//#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
//#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))

#define KAPPAES (1E3)
//#define KAPPAES (10.0)
//#define KAPPAES (1E-1)

#define KAPPAUSER(rho,T) (rho*KAPPA)
#define KAPPAESUSER(rho,T) (rho*KAPPAES)


#else // PULSE and PULSE3D

// KAPPAs are fraction of physical FF and ES opacities
#define KAPPA 0.
#define KAPPAES 1.e-30

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


#endif


#endif


//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADBEAMFLAT)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


#endif

//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADTUBE)

#define KAPPAESUSER(rho,T) (0.0)

#if(NTUBE==1)
#define KAPPAUSER(rho,T) (0.4*rho)
#elif(NTUBE==2)
#define KAPPAUSER(rho,T) (0.2*rho)
#elif(NTUBE==3)
#define KAPPAUSER(rho,T) (0.3*rho)
#elif(NTUBE==31)
#define KAPPAUSER(rho,T) (25*rho)
#elif(NTUBE==4)
#define KAPPAUSER(rho,T) (0.08*rho)
#elif(NTUBE==41)
#define KAPPAUSER(rho,T) (0.7*rho)
#elif(NTUBE==5)
#define KAPPAUSER(rho,T) (1000*rho)
#endif


#endif

//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADSHADOW || WHICHPROBLEM==RADDBLSHADOW)

#define KAPPAUSER(rho,T) (rho*1E2) // seems to allow photon build-up at front edge of blob
//#define KAPPAUSER(rho,T) (rho*1E0) // seems radiation pentrates blob too much compared to koral paper
#define KAPPAESUSER(rho,T) (rho*0.0)


#endif


//****************************************//
//****************************************//


//****************************************//
//****************************************//


#if(WHICHPROBLEM==RADBEAM2D)


#define KAPPA 0.
#define KAPPAES 0.

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


#endif





int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;
  int inittype;


  THETAROT = 0.0; // define rho,u,v,B as if no rotation


  inittype=1;

  funreturn=user1_init_primitives(inittype, prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);



  return(0);


}



int init_dsandvels(int inittype, int pos, int *whichvel, int*whichcoord, SFTYPE time, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_flatness(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

  // assume inittype not used, pos==CENT, and time doesn't matter (e.g. only used at t=0)

  init_dsandvels_flatness(whichvel, whichcoord,  i,  j,  k, pr, pstag);

  return(0);
}


// unnormalized density
int init_dsandvels_flatness(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  FTYPE X[NDIM],V[NDIM];
  struct of_geom realgeomdontuse;
  struct of_geom *ptrrealgeom=&realgeomdontuse;
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


	pr[URAD0] = 1./RHOBAR; // i.e. c^2 * 1g/cm^3 of energy density
	pr[URAD1] = 0 ;
	pr[URAD2] = 0 ;    
	pr[URAD3] = 0 ;

	*whichvel=WHICHVEL;
	*whichcoord=CARTMINKMETRIC2;
	return(0);
  }





  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAMFLAT){

    
	pr[RHO] = RADBEAMFLAT_RHO/RHOBAR ;
	pr[UU] = RADBEAMFLAT_UU/RHOBAR; // RADBEAMFLAT_UU was set in per c^2 units
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


	pr[PRAD0] = RADBEAMFLAT_ERAD/RHOBAR; // RADBEAMFLAT_ERAD was set in per c^2 units
	pr[PRAD1] = 0 ;
	pr[PRAD2] = 0 ;    
	pr[PRAD3] = 0 ;

	// no transformations required since only setting fluid-frame E that is PRAD0 itself.

	*whichvel=WHICHVEL;
	*whichcoord=CARTMINKMETRIC2;
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
    //	FTYPE RHO_AMB=(MPERSUN*MSUN/(LBAR*LBAR*LBAR)); // in grams per cm^3 to match koral's units with rho=1
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

	//	dualfprintf(fail_file,"IC i=%d Trad=%g ERAD=%g Tgas=%g rho=%g uint=%g\n",i,Trad,ERAD,Tgas,rho,uint);

   
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

	pr[URAD0] = ERAD ;
	pr[URAD1] = Fx ;
	pr[URAD2] = Fy ;    
	pr[URAD3] = Fz ;

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

	Fy=Fz=0.0;
	pr[URAD0] = ERAD ;
	pr[URAD1] = Fx ;
	pr[URAD2] = Fy ;    
	pr[URAD3] = Fz ;

	// setup vel type and coord type for prad_fforlab() based upon input velocity type and coordinate/metric type from data above
	*whichvel=VEL4;
	*whichcoord=CARTMINKMETRIC2;

	// TODO: need to convert these radiation things from the fluid frame (as defined) to the lab-frame
	// make a prad_ff2lab() like in koral's frames.c using Jon's new tetrad conversion stuff.
	// But in that case, here, lab frame is Minkowski so need to use gset() to make ptrgeom using CARTMINKMETRIC2.
  
	// get metric grid geometry for these ICs
	int getprim=0;
	struct of_geom geomdontuse;
	struct of_geom *ptrgeom=&geomdontuse;
	gset(getprim,*whichcoord,i,j,k,ptrgeom);

	// transform radiation primitives to lab-frame
	FTYPE prrad[NPR],prradnew[NPR];
	PLOOP(pliter,pl) prrad[pl]=pr[pl]; // prad_fforlab() should only use radiation primitives, but copy all primitives so can form ucon for transformation
	int whichframedir=FF2LAB; // fluid frame orthonormal to lab-frame
	prad_fforlab(*whichvel, *whichcoord, whichframedir, prrad, prradnew, ptrgeom);
	// overwrite radiation primitives with new lab-frame values
	PLOOPRADONLY(pl) pr[pl]=prradnew[pl];

	//  PLOOPRADONLY(pl) dualfprintf(fail_file,"FOO1: i=%d pl=%d pr=%g\n",ptrgeom->i,pl,pr[pl]);

	// inversion returns WHICHVEL velocity type, so pass that back
	*whichvel=WHICHVEL;
	*whichcoord=CARTMINKMETRIC2;
  
	return(0);
  

  }


  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADSHADOW || WHICHPROBLEM==RADDBLSHADOW){

	//	FTYPE MASS=10.0;
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
	
	//	FTYPE pamb=calc_PEQ_ufromTrho(TAMB,RHOAMB);
	rsq=xx*xx+yy*yy+zz*zz;
	rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;      
	Tgas=TAMB*RHOAMB/rho;
	// Paper says T = T0*rho/RHOAMB
	//	Tgas = TAMB*rho/RHOAMB;
	// for constant gas pressure, P=\rho T implies rho T = constant so that T\propto 1/rho
	uint=calc_PEQ_ufromTrho(Tgas,rho);

	Trad=TAMB;
	ERAD=calc_LTE_EfromT(Trad);

	//	dualfprintf(fail_file,"i=%d j=%d rho=%g Trad=%g uint=%g ERAD=%g\n",i,j,rho,Trad,uint,ERAD);

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

	pr[URAD0] = ERAD ;
	pr[URAD1] = Fx ;
	pr[URAD2] = Fy ;    
	pr[URAD3] = Fz ;

	*whichvel=VEL4;
	*whichcoord=CARTMINKMETRIC2;

	// get metric grid geometry for these ICs
	int getprim=0;
	struct of_geom geomdontuse;
	struct of_geom *ptrgeom=&geomdontuse;
	gset(getprim,*whichcoord,i,j,k,ptrgeom);


	// now need to transform these fluid frame E,F^i to lab frame coordinate basis primitives
	FTYPE prrad[NPR],prradnew[NPR];
	PLOOP(pliter,pl) prrad[pl]=pr[pl]; // prad_fforlab() should only use radiation primitives, but copy all primitives so can form ucon for transformation
	int whichframedir=FF2LAB; // fluid frame orthonormal to lab-frame
	prad_fforlab(*whichvel, *whichcoord, whichframedir, prrad, prradnew, ptrgeom);
	// overwrite radiation primitives with new lab-frame values
	PLOOPRADONLY(pl) pr[pl]=prradnew[pl];
   
	//	PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pr=%g\n",pl,pr[pl]);

	*whichvel=WHICHVEL;
	*whichcoord=CARTMINKMETRIC2;
	return(0);
  }



  /*************************************************/
  /*************************************************/
  if(WHICHPROBLEM==RADBEAM2D){


	FTYPE RHOAMB=1.e0/RHOBAR;
	FTYPE TAMB=1e7/TEMPBAR;
	int BLOB=0; // whether to put blob in way of beam
	FTYPE BLOBW=.1;
	FTYPE BLOBP=100000.;
	FTYPE BLOBX=10.;
	FTYPE BLOBZ=Pi/20.;
	FTYPE PAR_D=1./RHOBAR;
	FTYPE PAR_E=1e-4/RHOBAR;


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
	rhoblob=RHOAMB*(1.+BLOBP*exp(-((xx-BLOBX)*(xx-BLOBX)+(yy)*(yy)+(zz-BLOBZ)*(zz-BLOBZ))/BLOBW/BLOBW));

	//zaczynam jednak od profilu analitycznego:   
	FTYPE ERADAMB;
	FTYPE rho,uint,Vr;
	if(FLATBACKGROUND){
	  Vr=0.;
	  rho=RHOAMB;
	  if(BLOB) rho=rhoblob;
	  uint=calc_PEQ_ufromTrho(TAMB,rho);
	  ERADAMB=calc_LTE_EfromT(TAMB);
	  //	  ERADAMB=calc_LTE_Efromurho(uint,rho);
	}
	else{
	  FTYPE r=V[1];
	  FTYPE mD=PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
	  FTYPE mE=PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
	  Vr=sqrt(2./r)*(1.-2./r);

	  // get metric grid geometry for these ICs
	  int getprim=0;
	  struct of_geom geomdontuse;
	  struct of_geom *ptrgeom=&geomdontuse;
	  gset(getprim,*whichcoord,i,j,k,ptrgeom);

	  FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeom->gcov[GIND(1,1)]);
	  rho=PAR_D/(r*r*sqrt(2./r));
	  if(BLOB) rho += rhoblob;
	  FTYPE T=TAMB;
	  //			FTYPE ERAD=calc_LTE_EfromT(T);
	  uint=mE/W;
	  ERADAMB=calc_LTE_Efromurho(uint,rho);
	}



	//test!
	//Vr=0.7;

	FTYPE uradx,urady,uradz;
	uradx=urady=uradz=0.;
    
	pr[RHO] = rho ;
	pr[UU]  = uint;
	pr[U1]  = -Vr;
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


	pr[PRAD0] = ERADAMB;
	pr[PRAD1] = uradx ;
	pr[PRAD2] = urady ;    
	pr[PRAD3] = uradx ;

	// no transformations required since only setting fluid-frame E that is PRAD0 itself. (i.e. urad(xyz)=0)

	//	*whichvel=WHICHVEL;
	*whichvel=VEL4;
	*whichcoord=MCOORD;

	return(0);
  }





}


#define NOFIELD -1
#define DISK1FIELD 0
#define DISK2FIELD 1
#define VERTFIELD 2
#define DISK1VERT 3
#define DISK2VERT 4
#define BLANDFORDQUAD 5
#define TOROIDALFIELD 6

#define FIELDTYPE NOFIELD


// assumes normal field in pr
// SUPERNOTE: A_i must be computed consistently across all CPUs.  So, for example, cannot use randomization of vector potential here.
int init_vpot_user(int *whichcoord, int l, SFTYPE time, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av,q;
  FTYPE r,th,ph;
  FTYPE vpot;

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


      vpot += 0;
      
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
   return(0);
}


// assumes normal field definition
int normalize_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;

 
  //  funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
  //  if(funreturn!=0) return(funreturn);
 
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



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  int funreturn;
  
  // default is for spherical flow near BH
  //  funreturn=set_density_floors_default(ptrgeom, pr, prfloor);

  int pliter,pl;
  PLOOP(pliter,pl){
    prfloor[RHO]=1E-10*RHOMIN;
    prfloor[UU]=1E-10*UUMIN;

    prfloor[PRAD0]=ERADLIMIT;
  }

  //  if(funreturn!=0) return(funreturn);

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


//**********************************************************************
//******* user opacities ****************************************************
//**********************************************************************

// \kappa is optical depth per unit length per unit rest-mass energy density

//absorption
FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z)
{
  return(KAPPAUSER(rho,T));
}

//scattering
FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z)
{  
  return(KAPPAESUSER(rho,T));

}

