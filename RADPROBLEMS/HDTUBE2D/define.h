#define CYLINDRICAL
//#define FLATINIT
#define NX 100
#define NY 1
#define NZ 100
#define TSTEPLIM .5//courant limiter
#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define KAPPADIFF .1
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define INT_ORDER 1
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 1.e-1 //dt for basic output
#define ALLSTEPSOUTPUT 0
#define GAMMA (5./3.)
//nonrel:
//#define GAMMA (1.4)
#define MINX 10
#define MAXX 14
#define MINY -1.
#define MAXY 1.
#define MINZ -.3
#define MAXZ .3
//#define MAXZ (.1/((MAXX+MINX)/2.))
#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define ST_P1 1.
#define ST_P5 .1
/*
//non-rel
#define ST_RHO1 1.e5
#define ST_RHO5 .125e5
#define ST_U1 2.5
#define ST_U5 .25
*/
//rel
#define ST_RHO1 10.
#define ST_RHO5 1.
#define ST_U1 20.
#define ST_U5 1.e-7

double ST_P3;
#define YZXDUMP
#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define EFLOOR 1.e-40
