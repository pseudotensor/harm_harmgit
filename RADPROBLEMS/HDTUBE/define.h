#define MINKOWSKI
#define RK3STEPPING
#define SOURCETERMS
#define NX 500
#define NY 1
#define NZ 1
#define TSTEPLIM .5//kind of courant limiter
#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define KAPPADIFF .1
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1
#define INT_ORDER 1
#define DTOUT1 5.e0 //dt for basic output
#define ALLSTEPSOUTPUT 0
#define GAMMA (5./3.)
//nonrel:
//#define GAMMA (1.4)
#define MINX -50
#define MAXX 50
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.
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
//#define ST_U5 1.e-7
#define ST_U5 (1.e-8/(GAMMA-1.))

double ST_P3;





#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define EFLOOR 1.e-40
