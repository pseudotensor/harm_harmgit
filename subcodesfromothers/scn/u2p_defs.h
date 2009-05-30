
#include "decs.h"
#define GAMMA	gam  /* Adiabatic index used for the state equation */

/*Added these for the old version of HARM: */

#define FTYPE double

FILE *fail_file, *log_file;

//static int debugfail = 2;


#define NUMEPSILON (2.2204460492503131e-16)

/* loop over Primitive variables */
#define PLOOP for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
/* loop over all for j and Space for k; second rank loop */
#define DSLOOP for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)


#define OPTIMIZED 1

#define NEWCONVERGE 1

#define MAX_NEWT_ITER 30     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-10    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2

#define CYCLE_BREAK_PERIOD 1000 /* change newton step by random factor every this number of newton iterations*/
#define CHECK_FOR_STALL 0     /* check for stationary newton stepping */


#define NEWT_TOL2     1.0e-15      /* TOL of new DBL gnr2 method */
#define MIN_NEWT_TOL2 1.0e-10  /* TOL of new DBL gnr2 method */


#define SCALEMAX      1.0e2    /* Max. value of the factor used to scale the Newton step */
#define TOL_LINE_STEP NEWT_TOL /* Minimum value of Max(dx/x) in line search algorithm */
#define GRADMIN       1.0e-10  /* Magnitude of gradient below which we say we are at a local min. */
#define NEWT_FUNC_TOL 1.0e-5  /* Max. ratio of the final and initial resid magnitudes to be considered converged */


#define W_TOO_BIG	1.e9	/* \gamma^2 (\rho_0 + u + p) is assumed
                                  to always be smaller than this.  This
				  is used to detect solver failures */
#define UTSQ_TOO_BIG	1.e9    /* \tilde{u}^2 is assumed to be smaller
                                  than this.  Used to detect solver
				  failures */

#define FAIL_VAL  1.e30




/* some mnemonics */
/* for primitive variables */
#ifndef RHO
#define RHO 	0 
#endif

#ifndef UU
#define UU 	1 
#endif

#define UTCON1 	2
#define UTCON2 	3
#define UTCON3 	4
#define BCON1	5
#define BCON2	6
#define BCON3	7

/* for conserved variables */
#define QCOV0	1
#define QCOV1	2
#define QCOV2	3
#define QCOV3	4


#define MYMAX(a,b) ( ((a) > (b)) ? (a) : (b) )


