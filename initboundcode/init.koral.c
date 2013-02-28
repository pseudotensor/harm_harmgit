
#include "decs.h"


// NOTE on units:
// Many things below are so far in code units, not physical units, so they don't need conversion.  This includes:
//
// Rin, Rout, tf, DTdumpgen






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

  // print out units and some constants
  trifprintf("Constants\n");
  trifprintf("LBAR=%g TBAR=%g VBAR=%g RHOBAR=%g MBAR=%g UBAR=%g TEMPBAR=%g\n",LBAR,TBAR,VBAR,RHOBAR,MBAR,UBAR,TEMPBAR); 
  trifprintf("ARAD_CODE=%g OPACITYBAR=%g KAPPA_ES_CODE(1,1)=%g KAPPA_FF_CODE(1,1)=%g KAPPA_BF_CODE(1,1)=%g\n",ARAD_CODE,OPACITYBAR,KAPPA_ES_CODE(1,1),KAPPA_FF_CODE(1,1),KAPPA_BF_CODE(1,1));


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



int init_defcoord(void)
{

  /*************************************************/
  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==FLATNESS)


  defcoord = UNIFORMCOORDS;
  Rin_array[1]=0;
  Rin_array[2]=0;
  Rin_array[3]=0;

  Rout_array[1]=1.0;
  Rout_array[2]=1.0;
  Rout_array[3]=1.0;

#endif


  /*************************************************/
  /*************************************************/
  /*************************************************/

#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR)

  defcoord = UNIFORMCOORDS;
  Rin_array[1]=-50.0; 
  Rin_array[2]=-1.0;
  Rin_array[3]=-1.0;

  Rout_array[1]=50.0; 
  Rout_array[2]=1.0;
  Rout_array[3]=1.0;

#endif

  /*************************************************/
  /*************************************************/
  /*************************************************/

#if(WHICHPROBLEM==RADPULSE3D)

  defcoord = UNIFORMCOORDS;
  Rin_array[1]=-50.0; 
  Rin_array[2]=-50.0; 
  Rin_array[3]=-50.0; 

  Rout_array[1]=50.0; 
  Rout_array[2]=50.0; 
  Rout_array[3]=50.0; 

#endif

  /*************************************************/
  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==RADBEAMFLAT)


  defcoord = UNIFORMCOORDS;
  Rin_array[1]=0;
  Rin_array[2]=0;
  Rin_array[3]=0;

  Rout_array[1]=1.0;
  Rout_array[2]=1.0;
  Rout_array[3]=1.0;

#endif
  /*************************************************/
  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==RADTUBE)


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

#endif


  return(0);
}


int init_grid(void)
{
  

  return(0);
}



int init_global(void)
{
  int pl,pliter;
  int funreturn;

  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);


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

#if(WHICHPROBLEM==FLATNESS)  

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
#endif

  /*************************************************/
  /*************************************************/
  /*************************************************/


#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D)

  //  lim[1]=lim[2]=lim[3]=MINM;
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

  int idt;
  //  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=1E3;
  for(idt=0;idt<NUMDUMPTYPES;idt++) DTdumpgen[idt]=0.1;

  DTr = 100; //number of time steps for restart dumps
  tf = 1E2; //final time

#endif

  /*************************************************/
  /*************************************************/
  /*************************************************/

#if(WHICHPROBLEM==RADBEAMFLAT)  
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

#endif


  /*************************************************/
  /*************************************************/
  /*************************************************/

#if(WHICHPROBLEM==RADTUBE)
  cour=0.8;
  gam=gamideal=GAMMATUBE;
  cooling=KORAL;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=OUTFLOW;
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
  tf = 10.0; //final time

#endif

  /*************************************************/
  /*************************************************/
  /*************************************************/

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

//#define RHO_AMB (1.e0) // in grams per cm^3
#define RHO_AMB (MPERSUN*MSUN/(LBAR*LBAR*LBAR)) // in grams per cm^3 to match koral's units with rho=1
#define T_AMB (1.0E6) // in Kelvin

#define BLOBP 100.
#define BLOBW 5.



#if(WHICHPROBLEM==RADPULSEPLANAR)

#define KAPPA 0.
//#define KAPPAES (0.0)
//#define KAPPAES (1E-7)
#define KAPPAES (1E-4*1.09713E-18*1E3)
//#define KAPPAES (1E-4*1.09713E-18*1E-0)
//#define KAPPAES (1E-4*1.09713E-18*0.2)
//#define KAPPAES (1E-4*1.09713E-18*1E-1)
//#define KAPPAES (1E-4*1.09713E-18*1E-2)
//#define KAPPAES (1E-4*1.09713E-18*1E-3)
//#define KAPPAES (1E-4*1.09713E-18*1E-3*1E-10)

// assume KAPPA defines fraction of FF opacity
#define KAPPAUSER(rho,T) (rho*KAPPA*KAPPA_FF_CODE(rho,T))
// assume KAPPAES defines fractoin of ES opacity
#define KAPPAESUSER(rho,T) (rho*KAPPAES*KAPPA_ES_CODE(rho,T))


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



#if(WHICHPROBLEM==RADBEAMFLAT)

//#define RADBEAMFLAT_FRATIO 0.95
#define RADBEAMFLAT_FRATIO 0.995 // making like problem24 in koral code
#define RADBEAMFLAT_RHO 1. // 1g/cm^3
#define RADBEAMFLAT_ERAD 1. // 1g/cm^3 worth of energy density in radiation
#define RADBEAMFLAT_UU 0.1 // 0.1g/cm^3 worth of energy density in fluid

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
#if(WHICHPROBLEM==FLATNESS)  

  // outsideness
  if (0) {

    get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrrealgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {
    
    
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
#endif 






  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==RADBEAMFLAT)  

  // outsideness
  if (0) {

    get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrrealgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {
    
    
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


    pr[URAD0] = RADBEAMFLAT_ERAD/RHOBAR; // RADBEAMFLAT_ERAD was set in per c^2 units
    pr[URAD1] = 0 ;
    pr[URAD2] = 0 ;    
    pr[URAD3] = 0 ;

    *whichvel=WHICHVEL;
    *whichcoord=CARTMINKMETRIC2;
    return(0);
  }
#endif




  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D)

  // outsideness
  if (0) {

    get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrrealgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {

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

    // radiation temperature is distributed
    Trad=(T_AMB/TEMPBAR)*(1.+BLOBP*exp(-rsq/(BLOBW*BLOBW)));
    ERAD=calc_LTE_EfromT(Trad);

    //flat gas profiles
    Tgas=(T_AMB/TEMPBAR);
	FTYPE rho;
    rho=(RHO_AMB/RHOBAR);
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
#endif




  /*************************************************/
  /*************************************************/
#if(WHICHPROBLEM==RADTUBE)

  FTYPE rho,mx,my,mz,m,ERAD,uint,E0,Fx,Fy,Fz,pLTE;  
  FTYPE rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut,ux;
  FTYPE xx,yy,zz;
  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  

  if(xx<(Rout_array[1]+Rin_array[1])/2.){
    rho=1.;
    if(NTUBE==1) {uint = 3.e-5 / (GAMMATUBE - 1.); ERAD=1.e-8; Fx=1.e-2*ERAD;ux=0.015;}
    if(NTUBE==2) {uint = 4.e-3 / (GAMMATUBE - 1.);ERAD=2.e-5; Fx=1.e-2*ERAD;ux=0.25;}
    if(NTUBE==3 || NTUBE==31) {uint = 60. / (GAMMATUBE - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=10.;}
    if(NTUBE==4 || NTUBE==41) {uint = 6.e-3 / (GAMMATUBE - 1.);ERAD=0.18; Fx=1.e-2*ERAD;ux=0.69;}	  
    if(NTUBE==5) {uint = 60. / (GAMMATUBE - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=1.25;}
  }
  else{
	if(NTUBE==1) {rho=2.4;uint = 1.61e-4/ (GAMMATUBE - 1.); ERAD=2.51e-7; Fx=1.e-2*ERAD;ux=6.25e-3;}
	if(NTUBE==2) {rho=3.11;uint = 0.04512 / (GAMMATUBE - 1.);ERAD=3.46e-3; Fx=1.e-2*ERAD;ux=0.0804;}
	if(NTUBE==3 || NTUBE==31) {rho=8.0;uint = 2.34e3 / (GAMMATUBE - 1.);ERAD=1.14e3; Fx=1.e-2*ERAD;ux=1.25;}
	if(NTUBE==4 || NTUBE==41) {rho=3.65;uint =3.59e-2 / (GAMMATUBE - 1.);ERAD=1.30; Fx=1.e-2*ERAD;ux=0.189;}	  
	if(NTUBE==5) {rho=1.0;uint = 60. / (GAMMATUBE - 1.);ERAD=2.; Fx=1.e-2*ERAD;ux=1.10;}
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

  pr[URAD0] = ERAD ;
  pr[URAD1] = Fx ;
  pr[URAD2] = Fy ;    
  pr[URAD3] = Fz ;

  // TODO: need to conver these radiation things from the fluid frame (as defined) to the lab-frame
  // make a prad_ff2lab() like in koral's frames.c using Jon's new tetrad conversion stuff.
  // But in that case, here, lab frame is Minkowski so need to use gset() to make ptrgeom using CARTMINKMETRIC2.
  
  *whichvel=WHICHVEL;
  *whichcoord=CARTMINKMETRIC2;
  return(0);
  


#endif 




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

