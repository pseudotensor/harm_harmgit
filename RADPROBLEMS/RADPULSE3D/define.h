#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
//#define CYLINDRICAL
//#define KERR
//#define BHSPIN 0.
//#define SCHWARZSCHILD
//#define SPHERICAL
#define MINKOWSKI

#define RK3STEPPING
#define INT_ORDER 1
#define NX 51
#define NY 51
#define NZ 51
#define YSLICE 25
#define TSTEPLIM 1.
#define INITTSTEPLIM (TSTEPLIM/10.)
//#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 2.
#define DTOUT1 5.e0
#define ALLSTEPSOUTPUT 1
#define GAMMA (ldouble)(5./3.)
#define MINX -50.
#define MAXX 50.
#define MINY -50.
#define MAXY 50.
#define MINZ -50.
#define MAXZ 50.

#define RHO_AMB 1.e0
#define T_AMB 1.e6

#define EFLOOR 1.e-50
#define BLOBP 100.
#define BLOBW 5.
//#define YZSLICE

#define RADIATION
//#define GASRADOFF

#define KAPPA 0.
#define KAPPAES 1.e-6

//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

//#define EXPLICIT_RAD_SOURCE
#define IMPLICIT_FF_RAD_SOURCE

#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40
