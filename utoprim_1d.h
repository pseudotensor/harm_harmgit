#include "utoprim_1d2d.h"


#define NDIM 4
#define MYNPRINVERT 8

/* Dimension of or Number of Indep. variables used for the Newton-Raphson method of the primitive variable solver */
#define NEWT_DIM 1

#define SCALEMAX      1.0e2    /* Max. value of the factor used to scale the Newton step */
#define TOL_LINE_STEP NEWT_TOL /* Minimum value of Max(dx/x) in line search algorithm */
#define GRADMIN       1.0e-10  /* Magnitude of gradient below which we say we are at a local min. */
#define NEWT_FUNC_TOL 1.0e-5  /* Max. ratio of the final and initial resid magnitudes to be considered converged */

#define EXTRA_NEWT_ITER 2
#define USE_LINE_SEARCH_ALWAYS 0
/*

  USE_LINE_SEARCH = 0    line-searching is used only when normal NR fails to
  converge

  = 1     line-searching is always used
  < 0     line-searching is never used

*/








/* functions used in grmhd */
static void bcon_calc_g(FTYPE prim[MYNPRINVERT],FTYPE ucon[NDIM],FTYPE ucov[NDIM],FTYPE ncov[NDIM],FTYPE bcon[NDIM])  ;
static void lower_g(FTYPE vcon[NDIM], FTYPE gcov[SYMMATRIXNDIM], FTYPE vcov[NDIM]) ;
static void ncov_calc(FTYPE gcon[SYMMATRIXNDIM],FTYPE ncov[NDIM]) ;
static void raise_g(FTYPE vcov[NDIM], FTYPE gcon[SYMMATRIXNDIM], FTYPE vcon[NDIM]) ;
static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR], FTYPE *pressure);
static void ucon_calc_g(FTYPE prim[MYNPRINVERT],FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE ucon[NDIM]) ;

static FTYPE pressure_rho0_u_1d(FTYPE rho0, FTYPE u) ;
static FTYPE pressure_rho0_w_1d(FTYPE rho0, FTYPE w) ;
static FTYPE utsq_calc(FTYPE W) ;
static FTYPE gammasq_calc(FTYPE W) ;

static FTYPE dpress_dW( FTYPE W, FTYPE gamma, FTYPE dg );
static FTYPE dgamma_dW(FTYPE W, FTYPE gamma);
static FTYPE res_sq_1d_orig(FTYPE x[]);

static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
                                   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), 
                                   FTYPE (*res_func) (FTYPE []) );  

//int    is_nan_inf( FTYPE x );

static void primtoU_g(struct of_geom *ptrgeom, FTYPE prim[MYNPRINVERT], FTYPE gcov[SYMMATRIXNDIM], FTYPE gcon[SYMMATRIXNDIM],FTYPE U[MYNPRINVERT] );


static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n);

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  );


static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) );
