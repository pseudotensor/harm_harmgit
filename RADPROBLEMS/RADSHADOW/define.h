#define TMAX 1.e10
#define RADIATION
//#define MINKOWSKI
//#define CYLINDRICAL
//#define SPHERICAL
//#define SCHWARZSCHILD
#define NX 100
#define NY 50
#define NZ 1
#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING
#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define GAMMA (1.4)
#define DTOUT1 1.e0
#define EDDINGTON_APR

//#define EXPLICIT_RAD_SOURCE

//#define GASRADOFF
//#define RADSOURCEOFF

#define MASS 10.

#define MINX -1
#define MAXX 3
#define MINY -1
#define MAXY 1
#define MINZ -1.
#define MAXZ 1.

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

#define TAMB 1.e7
#define TLEFT 1.e9

#define RHOAMB 1.e-4

#define RHOBLOB 1.e3

#define BLOBW 5.e-2
#define KAPPA 1.e0

#define NLEFT 0.99999

#define SPECIFIC_BC

#define RADOUTPUTINZAMO
#define PRINTGC_LEFT

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-6
#define RADFORCEPREC 1.e-5
