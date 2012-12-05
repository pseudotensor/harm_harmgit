#define MINKOWSKI //metric

#define RADIATION //whether to solve for radiation (or pure hydro)
#define EDDINGTON_APR //Eddington approximation (P=1/3 I)
//#define EXPLICIT_RAD_SOURCE //whether to impose explicit treatment of the radiative four force terms
//#define IMPLICIT_FF_RAD_SOURCE //whether to use the explicit-implicit approximate implicit method

#define MASS 1./MSUNCM //mass in Msun of BH used for unit conversion
#define U2PPREC 1.e-7 //precision of the numerical hydro solver for conserved to primitives solver
#define U2PRADPREC 1.e-7 //precision of the numerical radiation converter, used only for Eddington approximation
#define RHOFLOOR 1.e-50 //floors (are not strictly enforced)
#define UFLOOR 1.e-45
#define EFLOOR 1.e-45

#define TMAX 1.e10 //end time
#define NOUTSTOP 1e3 //number of outputs to stop

#define NX 800 //x-resoltution
#define NY 1 //y-resolution
#define NZ 1 //z=rezolution

#define MINX -15. //size of the domain
#define MAXX 15.
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.

//#define SPECIFIC_BC //whether the boundary conditions are more complicated and given in PROBLEMS/XXX/bc.c
//#define PERIODIC_XBC //periodic BC
#define COPY_XBC //outflow BC
#define COPY_YBC
#define COPY_ZBC

#define TSTEPLIM .5 //Courant limiter
#define INITTSTEPLIM (TSTEPLIM/100.) //Courant limiter for the first step
#define RK3STEPPING //type of time stepping
#define INT_ORDER 1 //order of reconstruction - 1 (linear), 2 (parabolic), 4 (MP5)
#define FLUXLIMITER 0 //type of flux limiter for linear reconstruction - 0 (generalized minmod)
#define MINMOD_THETA 1. //theta for the limiter - 1 (minmod) - 2 (MC)

#define DTOUT1 2.e-0 //time step for outputs
#define ALLSTEPSOUTPUT 0 //0 (print every DTOUT1), 1 (print every step)
//#define PRINTGC_LEFT //if x-left ghost cells are to be printed out, 
//#define PRINTGC_RIGHT //if x-right ghost cells are to be printed out 
#define VERBOSE0 0 //verbose level for some routines

#undef MUGAS
#define MUGAS 1. //mean molecular weight (default = 1)

//problem specific definitions
#define NTUBE 2

#undef SIGMA_RAD
#if (NTUBE==1)
#define SIGMA_RAD (1e-8/powl(calc_PEQ_Tfromurho(3.e-5/(GAMMA-1.),1.),4.)/4.)
#define GAMMA 5.l/3.l
#elif (NTUBE==2)
#define SIGMA_RAD (2e-5/powl(calc_PEQ_Tfromurho(4.e-3/(GAMMA-1.),1.),4.)/4.)
#define GAMMA 5.l/3.l
#elif (NTUBE==3)
#define SIGMA_RAD (2./powl(calc_PEQ_Tfromurho(60./(GAMMA-1.),1.),4.)/4.)
#define GAMMA 2.l
#elif (NTUBE==31)
#define SIGMA_RAD (2./powl(calc_PEQ_Tfromurho(60./(GAMMA-1.),1.),4.)/4.)
#define GAMMA 2.l
#elif (NTUBE==4)
#define SIGMA_RAD (.18/powl(calc_PEQ_Tfromurho(6.e-3/(GAMMA-1.),1.),4.)/4.)
#define GAMMA 5.l/3.l
#elif (NTUBE==41)
#define SIGMA_RAD (.18/powl(calc_PEQ_Tfromurho(6.e-3/(GAMMA-1.),1.),4.)/4.)
#define GAMMA 5.l/3.l
#elif (NTUBE==5)
#define SIGMA_RAD (2./powl(calc_PEQ_Tfromurho(60./(GAMMA-1.),1.),4.)/4.)
#define GAMMA 2.l
#undef DTOUT1
#undef MINX
#undef MAXX
#define DTOUT1 1.e-0
#define MINX -20.
#define MAXX 20.
#endif




