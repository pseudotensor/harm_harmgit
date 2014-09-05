
/*! \file u2p_defs.h
    \brief U->P inversion definitions for all utoprim files
    
*/


#include "decs.h"
#define GAMMA (gamideal)  /* Adiabatic index used for the state equation */

// jon commented those out:

/*Added these for the old version of HARM: */

/* #define FTYPE double */

/* FILE *fail_file, *log_file; */

/* //static int debugfail = 2; */


/* #define NUMEPSILON (2.2204460492503131e-16) */

#define CHANGEDTOOLDER 0

#define CRAZYDEBUG 0



#define OPTIMIZED 1

#define NEWCONVERGE 1

#define DOHISTOGRAM 0

/// at what \gamma^2 does the inversion change character to become very accurate?
#define GAMMASQCHECKRESID (1e2)


// trying to control repeated cycles
//#define CYCLESTOP 10 // 10 is too dangerous
//#define NUMCYCLES 3 // 3 is too dangerous
#define CYCLESTOP 200 // don't limit for now to avoid issues with bad solutions in highly magnetized cases.
#define NUMCYCLES 200 // ""

/// Note that if using -pc64 -mp that error in inversion seems to be limited for doubles to 1E-11 instead of 1E-15
#if(PRECISEINVERSION)
#define ITERDAMPSTART 10 // iteration to start using damp if haven't already
#define MAX_NEWT_ITER 100     /* Max. # of Newton-Raphson iterations for find_root_2D();  Hardly ever case where need 200 iterations */
#define NEWT_TOL   (1E3*NUMEPSILON)    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define NEWT_TOL_ULTRAREL   (5.0*NUMEPSILON)    /* Min. of tolerance allowed for Newton-Raphson iterations */
//#define MIN_NEWT_TOL  (1E5*NUMEPSILON)    /* Max. of tolerance allowed for Newton-Raphson iterations */

#if(REALTYPE==FLOATTYPE)
#define MIN_NEWT_TOL  (1E-2) 
#elif(REALTYPE==DOUBLETYPE)
#define MIN_NEWT_TOL  (1E-6) 
#elif(REALTYPE==LONGDOUBLETYPE)
#define MIN_NEWT_TOL  (1E-8)
#endif


#else
#define ITERDAMPSTART 10 // iteration to start using damp if haven't already
#define MAX_NEWT_ITER 30     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   (1E5*NUMEPSILON)    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define NEWT_TOL_ULTRAREL   (1E1*NUMEPSILON)    /* Min. of tolerance allowed for Newton-Raphson iterations */
//#define MIN_NEWT_TOL  (1E5*NUMEPSILON)    /* Max. of tolerance allowed for Newton-Raphson iterations */

#if(REALTYPE==FLOATTYPE)
#define MIN_NEWT_TOL  (1E-2) 
#elif(REALTYPE==DOUBLETYPE)
#define MIN_NEWT_TOL  (1E-2) 
#elif(REALTYPE==LONGDOUBLETYPE)
#define MIN_NEWT_TOL  (1E-2) 
#endif


#endif

/// minimum relative error expected in \gamma
/// used to allow ultrarelativistic results to have larger relative errors since limited by machine precision
#define MINERROREXPECTED (NUMEPSILON*10.0)


/* #define MAX_NEWT_ITER 30     /\* Max. # of Newton-Raphson iterations for find_root_2D(); *\/ */
/* #define NEWT_TOL   1.0e-10    /\* Min. of tolerance allowed for Newton-Raphson iterations *\/ */
/* #define MIN_NEWT_TOL  1.0e-10    /\* Max. of tolerance allowed for Newton-Raphson iterations *\/ */

#if(PRECISEINVERSION)
/// notice for simple 1-D waves in 2D slab along direction that inversion gave poor results when compiled with -mp -pc64 and 0 extra iterations
/// Seems extra noise in field (with vpot for field) caused problems in inversion error estimates
/// even 1 extra iteration fixes this problem for that test
/// let's do 2 extra for precise inversion
#define EXTRA_NEWT_ITER 2
#define EXTRA_NEWT_ITER_ULTRAREL 2
#else
/// always do 1 extra to avoid error estimate issues as above
#define EXTRA_NEWT_ITER 1
#define EXTRA_NEWT_ITER_ULTRAREL 2
#endif

#define CYCLE_BREAK_PERIOD 1000 /* change newton step by random factor every this number of newton iterations*/
#define CHECK_FOR_STALL 0     /* check for stationary newton stepping */


#define NEWT_TOL2     (NUMEPSILON)      /* TOL of new DBL gnr2 method */
#define MIN_NEWT_TOL2 (1E4*NUMEPSILON)  /* TOL of new DBL gnr2 method */


#define SCALEMAX      1.0e2    /* Max. value of the factor used to scale the Newton step */
#define TOL_LINE_STEP NEWT_TOL /* Minimum value of Max(dx/x) in line search algorithm */
#define GRADMIN       (1E4*NUMEPSILON)  /* Magnitude of gradient below which we say we are at a local min. */
#define NEWT_FUNC_TOL 1.0e-5  /* Max. ratio of the final and initial resid magnitudes to be considered converged */


//#define W_TOO_BIG (GAMMAFAIL*GAMMAFAIL)   
//  \gamma^2 (\rho_0 + u + p) is assumed to always be smaller than this.  This
//      is used to detect solver failures





#define GAMMASQ_TOO_BIG (GAMMAFAIL*GAMMAFAIL)

#define UTSQ_TOO_BIG    ((GAMMAFAIL-1.0)*(GAMMAFAIL-1.0))
#define UT_TOO_BIG (GAMMAFAIL-1.0) 
/* \tilde{u}^2 is assumed to be smaller
   than this.  Used to detect solver
   failures */

//#define VSQ_TOO_BIG (1.0-1.0/((UT_TOO_BIG+1.0)*(UT_TOO_BIG+1.0)))  //atch: corrected to be a square
#define VSQ_TOO_BIG ((GAMMASQ_TOO_BIG-1.0)/(GAMMASQ_TOO_BIG)) // avoids catastrophic cancellation unlike above


#define FAIL_VAL  BIG

/// seems code can recover gracefully when this happens, so don't worry about it
/// which check
/// 0: no check
/// 1: original check
/// 2: new check allowing somewhat negative result
/// 3: just set to 0 if negative
#define VSQNEGCHECK 2

/// For <=-1 the equations become complex
#define UTSQNEGLIMIT (-0.9*UTSQ_TOO_BIG)

/// positive value of max negative utsq
//#define MAXNEGUTSQ (1E-10) // greater than negative of this but <0 makes utsq=0
#define MAXNEGUTSQ (fabs(UTSQNEGLIMIT)) // greater than negative of this but <0 makes utsq=0
//#define MAXNEGVSQ (1.0-1.0/(MAXNEGUTSQ+1.0))
#define MAXNEGVSQ (0.9*VSQ_TOO_BIG)



/* some mnemonics */
/* for primitive variables */
#ifndef RHO
#define RHO  0 
#endif

#ifndef UU
#define UU  1 
#endif

#define UTCON1  2
#define UTCON2  3
#define UTCON3  4
#define BCON1 5
#define BCON2 6
#define BCON3 7

/* for conserved variables */
#define QCOV0 1
#define QCOV1 2
#define QCOV2 3
#define QCOV3 4


#define MYMAX(a,b) ( ((a) > (b)) ? (a) : (b) )


