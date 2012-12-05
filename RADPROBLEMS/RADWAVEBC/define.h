#define U2PPREC 1.e-7
#define U2PRADPREC 1.e-7
#define VERBOSE0 0
#define MINKOWSKI

#define RK3STEPPING
#define INT_ORDER 1
#define NX 800
#define NY 1
#define NZ 1
#define PRINTGC_LEFT
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)

#define SPECIFIC_BC
#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define ALLSTEPSOUTPUT 0

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 3.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

//#define EDDINGTON_APR
#define RADIATION

#define KAPPAES 0.

#define NWAVE 1
//#define RADOUTPUTINZAMO

#undef SIGMA_RAD

#if (NWAVE==1)
#define PP 0.01
#define CC 1.e2
#define KK 2.*Pi
#define KAPPA 0.01
#define OMEGA KK/CC
#define RHO 1.
#define DRHO 1.e-3
#define ERAD 1.
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#define DTOUT1 10.

#endif


#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40
