#define U2PPREC 1.e-7
#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-7
#define VERBOSE0 0
//#define CYLINDRICAL
//#define KERR
#define BHSPIN 0.
#define SCHWARZSCHILD
//#define SPHERICAL
//#define MINKOWSKI
#define RADIATION
#define ANAL_PROFILE
#define NX 212
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define SPECIFIC_BC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 1.e0
#define ALLSTEPSOUTPUT 0
#define GAMMA (ldouble)(5./3.)
#define KAPPAES kappaCGS2GU(0.4)
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
//#define LOGXGRID
#define MINX 10.
#define MAXX 40.
#define MINY 0*Pi/2.
#define MAXY 1.*Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define RHOFLOOR 1.e-20
#define UFLOOR 1.e-15
#define PAR_D 1.e0
#define PAR_E 1.e-2
#define RHO_AMB 1.e-3
#define U_AMB 1.e-7
#define BLOB_R 35
#define BLOB_D .2
#define BLOB_RHORATIO 0.
#define EFLOOR 1.e-40
//#define CGSOUTPUT
#define PRINTGC_LEFT
#define PRINTGC_RIGHT
#define MASS 3.
#define MDOTEDD 2.23/16.*1e18*MASS //cm/s
#define MDOT 0.001
#define LUMEDD 1.25e38*MASS //erg/s
#define LUM 10.
#define TAMB 1.e7
