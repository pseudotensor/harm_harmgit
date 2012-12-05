#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-7
#define U2PPREC 1.e-7

#define TMAX 1.e10
#define RADIATION
#define SCHWARZSCHILD
#define RK3STEPPING
//#define SPHERICAL
//#define MINKOWSKI
//#define SKIPLORENTZ
#define NX 112
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define INT_ORDER 1
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 1.e1
#define ALLSTEPSOUTPUT 0
#define VERBOSE0 0

#define TESTNO 2

#if (TESTNO==1)
#define PRADGAS 1.2e-7
#define TGAS0 1e5
#define MDOT 10.
#endif

#if (TESTNO==2)
#define PRADGAS 1.2e-4
#define TGAS0 1.e6
#define MDOT 10.
#endif

#if (TESTNO==3)
#define PRADGAS 1.2e-1
#define TGAS0 1e7
#define MDOT 10.
#endif

#if (TESTNO==4)
#define PRADGAS 1.2e-5
#define TGAS0 1e6
#define MDOT 100.
#endif 

#define MASS 3.
#define MDOTEDD 2.23/16.*1e18*MASS //cm/s
#define RHOAMB 1.e-25
#define TAMB 1.e5
#define GAMMA (long double)(1.+1./3.*((1.+PRADGAS)/(.5+PRADGAS)))
#undef MUGAS
#define MUGAS .5
//#define EXPLICITRADFORCE
//#define EDDINGTON_APR
#define MINX 3.5
#define MAXX 2e3
#define LOGXGRID
#define LOGPAR1 2.2
#define LOGPAR2 2.

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define SPECIFIC_BC
//#define COPY_XBC
//#define COPY_YBC
//#define COPY_ZBC
#define EFLOOR 1.e-40
#define CGSOUTPUT

//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT
//#define OUTPUTINZAMO
