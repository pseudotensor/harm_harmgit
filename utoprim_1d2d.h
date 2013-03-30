#if(PRECISEINVERSION)
#define MAX_NEWT_ITER 200     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-15    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#else
#define MAX_NEWT_ITER 20     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-10    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-4    /* Max. of tolerance allowed for Newton-Raphson iterations */
#endif

// whether showing extra information and doing extra checks
// 0=do 1=don't
#define OPTIMIZED 1

#define CYCLE_BREAK_PERIOD 10 /* change newton step by random factor every this number of newton iterations*/  
#define CHECK_FOR_STALL 1     /* check for stationary newton stepping */                                       

#define SCALEMAX      1.0e2    /* Max. value of the factor used to scale the Newton step */
#define TOL_LINE_STEP NEWT_TOL /* Minimum value of Max(dx/x) in line search algorithm */
#define GRADMIN       1.0e-10  /* Magnitude of gradient below which we say we are at a local min. */
#define NEWT_FUNC_TOL 1.0e-5  /* Max. ratio of the final and initial resid magnitudes to be considered converged */


#define W_TOO_BIG 1.e9 /* \gamma^2 (\rho_0 + u + p) is assumed
                          to always be smaller than this.  This
                          is used to detect solver failures */
#define UTSQ_TOO_BIG ((GAMMAFAIL-1.0)*(GAMMAFAIL-1.0))
#define UT_TOO_BIG (GAMMAFAIL-1.0)      /* \tilde{u}^2 is assumed to be smaller
                                           than this.  Used to detect solver
                                           failures */
#define VSQ_TOO_BIG (1.0-1.0/UT_TOO_BIG)

#define MAXNEGUTSQ (1E-10) // greater than negative of this but <0 makes utsq=0
#define MAXNEGVSQ (1.0-1.0/(MAXNEGUTSQ+1.0))


#define FAIL_VAL (1.e30)








/* some mnemonics */
/* for primitive variables */
/* #define RHO  0  */
/* #define UU   1 */

#define UTCON1  2
#define UTCON2  3
#define UTCON3  4
#define BCON1   5
#define BCON2   6
#define BCON3   7

/* for conserved variables */
#define QCOV0   1
#define QCOV1   2
#define QCOV2   3
#define QCOV3   4





static FTYPE minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?      \
                   (minarg1) : (minarg2))

static FTYPE maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?      \
                   (maxarg1) : (maxarg2))


