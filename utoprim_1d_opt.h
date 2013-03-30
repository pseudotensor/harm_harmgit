
#define OPTIMIZED 1

#define NDIM 4
#define MYNPRINVERT 8

#if(PRECISEINVERSION)
#define MAX_NEWT_ITER 200     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-15    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#else
#define MAX_NEWT_ITER 30     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-10    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#endif

#define CYCLE_BREAK_PERIOD 1000 /* change newton step by random factor every this number of newton iterations*/  
#define CHECK_FOR_STALL 1     /* check for stationary newton stepping */                                       
#define EXTRA_NEWT_ITER 2

#define USE_LINE_SEARCH 0
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

static FTYPE minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?      \
                   (minarg1) : (minarg2))

static FTYPE maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?      \
                   (maxarg1) : (maxarg2))

/* some mnemonics */
/* for primitive variables */
/* #define RHO  0  */
/* #define UU  1 */

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


/* Dimension of or Number of Indep. variables used for the Newton-Raphson method of the primitive variable solver */
#define NEWT_DIM 1


/* functions used in grmhd */
static void bcon_calc_g(FTYPE prim[MYNPRINVERT],FTYPE ucon[NDIM],FTYPE ucov[NDIM],FTYPE ncov[NDIM],FTYPE bcon[NDIM])  ;
static void lower_g(FTYPE vcon[NDIM], FTYPE gcov[SYMMATRIXNDIM], FTYPE vcov[NDIM]) ;
static void ncov_calc(FTYPE g00,FTYPE ncov[NDIM]) ;
static void raise_g(FTYPE vcov[NDIM], FTYPE gcon[SYMMATRIXNDIM], FTYPE vcon[NDIM]) ;
static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR], FTYPE *pressure);
static void ucon_calc_g(FTYPE prim[MYNPRINVERT],FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE ucon[NDIM]) ;

static FTYPE pressure_rho0_u_1dopt(FTYPE rho0, FTYPE u) ;
static FTYPE pressure_rho0_w_1dopt(FTYPE rho0, FTYPE w) ;
static FTYPE utsq_calc(FTYPE W) ;
static FTYPE gammasq_calc(FTYPE W) ;

static FTYPE dpress_dW( FTYPE W, FTYPE gamma, FTYPE dg );
static FTYPE dgamma_dW(FTYPE W, FTYPE gamma);
static FTYPE res_sq_1d_orig(FTYPE x[]);

static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
                                   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), 
                                   FTYPE (*res_func) (FTYPE []) );  

//int    is_nan_inf( FTYPE x );

static void primtoU_g( FTYPE prim[MYNPRINVERT], FTYPE gcov[SYMMATRIXNDIM], FTYPE gcon[SYMMATRIXNDIM],FTYPE U[MYNPRINVERT] );


static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n);

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  );


static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) );
