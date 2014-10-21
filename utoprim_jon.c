
/*! \file utoprim_jon.c
  \brief Jon's U->P inversion
     
  // OPENMPNOTE: This function should only have local variables in order to be parallelized, otherwise global shared variables can change while computing.
  // OPENMPNOTE: Should avoid use of static variables inside functions, for example when first time called to initialize things.  Only done right now in bin_newt_data() that is not used.
  */


#include "u2p_defs.h" // only includes #define's
#include "utoprim_jon.h" // only includes #define's


/// 0 = no debug indexes set
/// 1 = set global (file scope) index for debugging
///     OPENMPMARK: This is not thread safe!, so only set to 1 if not using more than one OpenMP thread.
#define DEBUGINDEX 0

#if(DEBUGINDEX)
// 
static int ifileglobal,jfileglobal,kfileglobal,pfileglobal;
#endif


/// Main function called:
static int Utoprim_new_body(int showmessages, int eomtype, PFTYPE *lpflag, int whicheos, FTYPE *EOSextra, FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats);

/// Metric-related things:
static void raise_g(FTYPE vcov[], FTYPE *gcon, FTYPE vcon[]);
static void lower_g(FTYPE vcon[], FTYPE *gcov, FTYPE vcov[]);
static void ncov_calc_fromlapse(FTYPE lapse,FTYPE ncov[]) ;
static void ncov_calc(FTYPE gcon[SYMMATRIXNDIM],FTYPE ncov[]) ;
static int coldgrhd(PFTYPE *lpflag, FTYPE Qtsq, FTYPE D, FTYPE *Wp);
static int vsqgtr1(FTYPE W,FTYPE Bsq,FTYPE QdotBsq, FTYPE Qtsq);

/// Kinetic-related things:
static FTYPE vsq_calc(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE dvsq_dW(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE Eprime_Wp(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int Psq_Wp(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, FTYPE *resid, FTYPE *norm);
static FTYPE dPsqdWp_Wp(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE utsq_calc(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE gammasq_calc_W(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE utsq_calc(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);

/// validate solution:
static void validate_Wp(FTYPE Wpold, FTYPE *Wpnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void check_utsq_fail(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);

/// EOS-related things for hot GRMHD evolution
static FTYPE pressure_W_vsq(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE pressure_Wp_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dpdW_calc_vsq(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dpdWp_calc_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dpdvsq_calc2_Wp(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dpdvsq_calc(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dpdWp_calc_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);

/// EOS-related things for entropy evolution
static FTYPE Ss_W_vsq(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE Ss_Wp_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dSsdW_calc_vsq(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dSsdWp_calc_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dSsdvsq_calc2_Wp(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);
static FTYPE dSsdvsq_calc(int whicheos, FTYPE *EOSextra, FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dSsdWp_calc_utsq(int whicheos, FTYPE *EOSextra, FTYPE Wp, FTYPE D, FTYPE utsq);

/// Find Root wrappers:
static int find_root_1D_gen(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_1D_gen_scn(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_1D_gen_Eprime(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_1D_gen_Psq(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_3D_gen_Palpha(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_2D_gen(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);
static int find_root_1D_gen_Sc(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);

/// Newton methods:
static int general_newton_raphson(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x[], int n, int do_line_search,
                                  void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), 
                                  FTYPE (*res_func) (FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE res_sq_1d_orig( FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void func_1d_orig_scn(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE res_sq_1d_orig_scn( FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE res_sq_1d_orig( FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE res_sq_1d_orig( FTYPE [] , FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static FTYPE res_sq_1d_orig( FTYPE [] , FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static  void func_vsq(   FTYPE [], FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static  void func_vsq2(   FTYPE [], FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static  FTYPE  res_sq_vsq( FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static  FTYPE  res_sq_vsq2( FTYPE [] , FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int x1_of_x0(FTYPE x0, FTYPE *x1, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void my_lnsrch(int eomtype, int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
                      FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);


/// Newton checks:
static int newt_repeatcheck(int n, FTYPE errx, FTYPE errx_old, FTYPE errx_oldest, FTYPE *dx, FTYPE *dx_old, FTYPE *x, FTYPE *x_old, FTYPE *x_older, FTYPE *x_olderer);
static int newt_cyclecheck(int n, FTYPE trueerror, FTYPE errx, FTYPE errx_old, FTYPE errx_oldest, FTYPE *dx, FTYPE *dx_old, FTYPE *x, FTYPE *x_old, FTYPE *x_older, FTYPE *x_olderer, int n_iter, FTYPE (*xlist)[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER], FTYPE (*dxlist)[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER], FTYPE *errxlist, FTYPE *trueerrorlist);
static int newt_errorcheck(int n, FTYPE NEWT_TOL_VAR, FTYPE NEWT_TOL_ULTRAREL_VAR, FTYPE errx, FTYPE trueerrortotal, FTYPE x0, FTYPE dx0, FTYPE *x, FTYPE *dx, FTYPE *wglobal);
static int newt_extracheck(int n, int EXTRA_NEWT_ITER_VAR, int EXTRA_NEWT_ITER_ULTRAREL_VAR, FTYPE errx, FTYPE x0, FTYPE dx0, FTYPE *x, FTYPE *dx, FTYPE *wglobal);
static void bin_newt_data(FTYPE MAX_NEWT_ITER_VAR, FTYPE errx, int niters, int conv_type, int print_now  ) ;

static void validate_x(void (*ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra),FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static void check_on_inversion(int showmessages, FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats);



/// Setup inversion:
static int compute_setup_quantities(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, FTYPE *Bsq,FTYPE *QdotB,FTYPE *QdotBsq,FTYPE *Qtsq,FTYPE *Qdotn,FTYPE *Qdotnp,FTYPE *D, FTYPE *Sc, int whicheos, FTYPE *EOSextra);
static int compute_Wp(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *Wpabs, FTYPE *etaabs, FTYPE *gamma, FTYPE *gammasq, FTYPE *rho0, FTYPE *u, FTYPE *p, FTYPE *wmrho0, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int set_guess_Wp(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int check_Wp(PFTYPE *lpflag, int eomtype, FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp_last, FTYPE Wp, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int Wp2prim(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, FTYPE *pressure, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
static int verify_Wlast(FTYPE u, FTYPE p, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);


/// other inversions:
static int forcefree_inversion(int cleaned, struct of_geom *ptrgeom, FTYPE W, FTYPE *Qtcon, FTYPE Bsq, FTYPE *Bcon, FTYPE *Bcov, FTYPE Qtsq, FTYPE *U, FTYPE *prim, FTYPE *gammaret);





#include "u2p_util.c" // only includes static functions
#include "utoprim_jon_eos.c" // only includes static functions


/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

/// The ONLY global function:
int Utoprim_jon_nonrelcompat_inputnorestmass(int showmessages, int eomtype, FTYPE *EOSextra, FTYPE *U, struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats)
{

  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, k,ret; 
  int pl,pl2;
  FTYPE gcov[SYMMATRIXNDIM], gcon[SYMMATRIXNDIM], alpha;
  const int ltrace = 0;
  /* these variables need to be shared between the functions
     Utoprim_1D, residual, and utsq */
  int real_dimen_newton;
  int whicheos;
  int pliter;

  //  showmessages=1; // FORCE

  // DEBUG:
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"%d %d %d : pl=%d U=%21.15g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pl,U[pl]);


#if(0)
  U[0]=      100.03640347219;
  U[1]=-2.21974020354085e-25;
  U[2]=-2.50158167405965e-12;
  U[3]=                    0;
  U[4]=                    0;
  U[5]=                    0;
  U[6]=                    0;
  U[7]=                    0;

  prim[0]=     99.8798655540244;
  prim[1]=  9.6414760309742e-26;
  prim[2]=-9.98706197921322e-14;
  prim[3]=                    0;
  prim[4]=                    0;
  prim[5]=                    0;
  prim[6]=                    0;
  prim[7]=                    0;
#endif



  ////////////////////////
  //
  // EOS functions used during inversion
  // GODMARK: WHICHEOS could be input instead of using global definition
  //
  ////////////////////////
  pickeos_eomtype(WHICHEOS,eomtype,&whicheos);




#if(!OPTIMIZED)
  // make sure nan and inf's are not in initial guess or U
  PLOOPALLINVERT(pl){
    if(!isfinite(prim[pl])){
      dualfprintf(fail_file,"Guess for p is NAN\n");
      PLOOPALLINVERT(pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,prim[pl2]);
      *lpflag=  UTOPRIMFAILNANGUESS;
      myexit(876);
    }
    if(!isfinite(U[pl])){
      dualfprintf(fail_file,"Value of U is NAN\n");
      PLOOPALLINVERT(pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,U[pl2]);
      *lpflag=  UTOPRIMFAILNANGUESS;
      myexit(877);
    }
  }
#endif



  /////////////////
  //
  // initialize iteration counter
  //
  /////////////////
  newtonstats->lntries=0;



  /////////////////
  //
  // Setup dimensionality of method for a given eomtype and choice if inverter
  //
  /////////////////
  if(eomtype==EOMGRMHD){
    if(WHICHHOTINVERTER==1) real_dimen_newton=1;
    else   if(WHICHHOTINVERTER==2) real_dimen_newton=2;
    else   if(WHICHHOTINVERTER==3) real_dimen_newton=1;
    else{
      dualfprintf(fail_file,"No such hot inverter %d\n",WHICHHOTINVERTER);
      myexit(89285323);
    }
  }
  else if(eomtype==EOMCOLDGRMHD){
    if(WHICHCOLDINVERTER==0) real_dimen_newton=1;
    else if(WHICHCOLDINVERTER==1) real_dimen_newton=1;
    else   if(WHICHCOLDINVERTER==2) real_dimen_newton=1;
    else   if(WHICHCOLDINVERTER==3) real_dimen_newton=1;
    else{
      dualfprintf(fail_file,"No such cold inverter %d\n",WHICHCOLDINVERTER);
      myexit(89285324);
    }
  }
  else if(eomtype==EOMENTROPYGRMHD){
    if(WHICHENTROPYINVERTER==0) real_dimen_newton=1;
    else{
      dualfprintf(fail_file,"No such entropy inverter %d\n",WHICHENTROPYINVERTER);
      myexit(89285325);
    }
  }
  else real_dimen_newton=0; // not using Newton's method



#if( WHICHVEL != VELREL4 ) 
  stderrfprintf("Utoprim_2d() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1); 
#endif



  ///////////
  //
  // First update the primitive B-fields
  //
  ///////////
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  
 

  ///////////
  //
  // Calculate the transform the CONSERVATIVE variables into the new system
  //
  ///////////

  // Set the geometry variables
  //  alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(0,0)]);
  alpha = ptrgeom->alphalapse;

  for( i = RHO; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] ;
  }

#if(DOENTROPY!=DONOENTROPY)
  i=ENTROPY;
  U_tmp[i] = alpha * U[i] ;
#endif

 
#if(CRAZYDEBUG)
  if(ptrgeom->i==0 && ptrgeom->j==31 && nstep==9 && steppart==2){
    PLOOPALLINVERT(k) dualfprintf(fail_file,"U_tmp[%d]=%21.15g\n",k,U_tmp[k]);
  }
#endif

 
  /////////////
  //
  // Calculate the transform the PRIMITIVE variables into the new system
  //
  /////////////

  // default is to treat all quantities as scalars
  // Note that need PALLLOOP here since while not setting all quantities, may use all quantities
  PALLLOOP(i) prim_tmp[i] = prim[i];
  // field
  for( i = BCON1; i <= BCON3; i++ ) prim_tmp[i] = alpha*prim[i];



  /////////////
  //
  // Perform the inversion from U -> P
  //
  /////////////
  ret = Utoprim_new_body(showmessages, eomtype, lpflag, whicheos,EOSextra, U_tmp, ptrgeom, prim_tmp,pressure,newtonstats);



  //////////
  //
  // Check conservative variable transformation
  //
  //////////
#if(!OPTIMIZED)
  if( ltrace ) {
    for(i = 0; i < NDIM; i++ ) {
      for(j = 0; j < NDIM; j++ ) {
        gcov[GIND(i,j)] = ptrgeom->gcov[GIND(i,j)];
        gcon[GIND(i,j)] = ptrgeom->gcon[GIND(i,j)];
        dualfprintf(fail_file,"gcov,gcon %d %d = %21.15g %21.15g \n", i, j, gcov[GIND(i,j)], gcon[GIND(i,j)]);
      }
    }
    dualfprintf(fail_file,"gdet = %21.15g \n", ptrgeom->g);
    //    primtoU_g(ptrgeom, prim_tmp, gcov, gcon, U_tmp2 ); 
    //    PALLLOOP(i){
    //      dualfprintf(fail_file, "Utoprim_1d(): Utmp1[%d] = %21.15g , Utmp2[%d] = %21.15g , dUtmp[%d] = %21.15g \n", 
    //        i, U_tmp[i], i, U_tmp2[i], i, fabs( (U_tmp[i]-U_tmp2[i]) / ( (U_tmp2[i]!=0.) ? U_tmp2[i] : 1. ) )  ); 
    //    }
  }
#endif



  //////////
  //
  // Transform new primitive variables back
  //
  ///////////
  if( ret == 0 ) {
    // only changed RHO through U3 (otherwise have to remove alpha from B^i)
    for( i = 0; i <=U3; i++ ) prim[i] = prim_tmp[i];
#if(DOENTROPY!=DONOENTROPY)
    // if entropy is also copy entropy primitive that is internal energy density
    i=ENTROPY;
    prim[i] = prim_tmp[i];
#endif

  }
  else{
    // ensure failure reported if ret!=0
    if(IFUTOPRIMNOFAILORFIXED(*lpflag)){
      //      dualfprintf(fail_file,"Internal ret=%d failure\n",ret);
      *lpflag= ret+UTOPRIMFAILCONVRET+1;// related to UTOPRIMFAILCONVRET
    }
  }


  //  PLOOP(pliter,pl) dualfprintf(fail_file,"internal pl=%d prim=%g\n",pl,prim[pl]);


  ////////
  //
  // revert EOS for WHICHEOS,EOMTYPE
  // not really  necessary unless whicheos used after the below call
  //
  /////////
  pickeos_eomtype(WHICHEOS,EOMTYPE,&whicheos); // see comments in eos.c for pickeos_eomtype_old()

  return( 0 ) ;

}


/**********************************************************************************

attempt an inversion from U to prim using the initial guess
prim.

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu  |
             \  alpha B^i        /
              \     ...        /
              \  specific entropy /


             /    rho        \
  P = |    uu         |
             | \tilde{u}^i   |
             \  alpha B^i   /
              \     ...        /
              \   uu from entropy (copied to uu) /


              *lpflag:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence (occurrence of "nan" or "+/-inf" ;
      
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
      3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: utsq<0 or utsq > UTSQ_TOO_BIG   with new  W;
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(int showmessages, int eomtype, PFTYPE *lpflag, int whicheos, FTYPE *EOSextra, FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats)
{
  FTYPE Qconp[NDIM],Qcovp[NDIM];
  FTYPE Qtcon[4],Bcon[4],Bcov[4];
  FTYPE primother[NPR];
  FTYPE W_last,W,Wp_last,Wp0,Wp,Wp1,Wp2,Wp3,Wpother;
  int i,j, retval0, retval, retval1,retval2, retval3 ;
  int retvalother,retvalother2;
  // below used to be static global in this file
  FTYPE Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc ;
  FTYPE wglobal[3]; // [0] is normalization for Wp and [1] is normalization for Wtest, and [2] is for rho0gammasqresidabs term
  int Wdependentsolution;


#if(0) // old globals used for debugging -- can't be used during OpenMP computations
  static FTYPE Qcovorig[NDIM],Qconorig[NDIM],Qdotnorig,Qsqorig,Qtsqorig;
  static FTYPE Qtsqnew;
#endif


   
#if(DEBUGINDEX)
  // initialize debug indices
  ifileglobal=ptrgeom->i;
  jfileglobal=ptrgeom->j;
  kfileglobal=ptrgeom->k;
  pfileglobal=ptrgeom->p;
#endif

  const int ltrace = 0;
  const int ltrace2 = 0;





#if(!OPTIMIZED)
  if( ltrace ) {
    for(i=0;i<NPR;i++) dualfprintf(fail_file,"%d %21.15g %21.15g\n",i,prim[i],U[i]) ;
    
  }
#endif




  //////////
  //
  // Assume no inversion failure initially
  //
  //////////
  *lpflag= UTOPRIMNOFAIL;
  retval=0;


  //////////////
  //
  // compute quantities necessary for inversion in a certain frame
  //
  ///////////////
  compute_setup_quantities(prim, U, ptrgeom, Qtcon, Bcon, Bcov, &Bsq,&QdotB,&QdotBsq,&Qtsq,&Qdotn,&Qdotnp,&D,&Sc,whicheos,EOSextra);




#if(0)
  // DEBUG:
  //  if(ptrgeom->i==1 && ptrgeom->j==18){
  DLOOPA(i) dualfprintf(fail_file,"Qtcon[%d]=%21.15g Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",i,Qtcon[i],i,Bcon[i],i,Bcov[i]);
  dualfprintf(fail_file,"Bsq=%21.15g QdotB=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Qdotn=%21.15g Qdotnp=%21.15g D=%21.15g Sc=%21.15g\n",Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc);
  //  }
#endif





  ////////////////////////////////
  //
  // CHECK WHICH EOM TYPE and proceed accordingly
  //
  ////////////////////////////////

  if(eomtype==EOMFFDE || eomtype==EOMFFDE2){
    Wdependentsolution=0;
    if(eomtype==EOMFFDE){
      FTYPE gammaret0;
      W_last=0.0;
      int cleaned=1;
      if(forcefree_inversion(cleaned, ptrgeom, W_last, Qtcon, Bsq, Bcon, Bcov, Qtsq, U, prim, &gammaret0)) retval++;
    }
    else if(eomtype==EOMFFDE2){ // does not so good.
      int cleaned=0;
      FTYPE gammaret0;

      // get Wp, wglobal, etaabs, gamma, gammaabs
      FTYPE gammasq,gamma,rho0,u,p,wmrho0;
      FTYPE Wpabs, etaabs;
      compute_Wp(showmessages, lpflag, eomtype, prim, ptrgeom, &W_last, &Wp_last, &Wpabs, &etaabs, &gamma, &gammasq, &rho0, &u, &p, &wmrho0, wglobal, Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

      if(forcefree_inversion(cleaned, ptrgeom, W_last, Qtcon, Bsq, Bcon, Bcov, Qtsq, U, prim, &gammaret0)) retval++;
    }
    else{ // doesn't do as well as EOMFFDE, so smallest instant value of gamma not only important thing, since eventually leads to larger gamma.

      // get Wp, wglobal, etaabs, gamma, gammaabs
      FTYPE gammasq,gamma,rho0,u,p,wmrho0;
      FTYPE Wpabs, etaabs;
      compute_Wp(showmessages, lpflag, eomtype, prim, ptrgeom, &W_last, &Wp_last, &Wpabs, &etaabs, &gamma, &gammasq, &rho0, &u, &p, &wmrho0, wglobal, Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

      // then perform the analytic inversion with Lorentz factor limit
      FTYPE gammaret1;
      FTYPE prim1[NPR];
      int pliter,pl;
      PLOOP(pliter,pl) prim1[pl] = prim[pl];
      int ret1=forcefree_inversion(0, ptrgeom, W_last, Qtcon, Bsq, Bcon, Bcov, Qtsq, U, prim1, &gammaret1);

      FTYPE gammaret2;
      FTYPE prim2[NPR];
      PLOOP(pliter,pl) prim2[pl] = prim[pl];
      int ret2=forcefree_inversion(1, ptrgeom, W_last, Qtcon, Bsq, Bcon, Bcov, Qtsq, U, prim2, &gammaret2);
      if(gammaret1>gammaret2){
        PLOOP(pliter,pl) prim[pl] = prim2[pl];
        if(ret2) retval++;
      }
      else{
        PLOOP(pliter,pl) prim[pl] = prim1[pl];
        if(ret1) retval++;
      }
    }
  }
  else{
    Wdependentsolution=1;

    /////////////
    //
    // SETUP ITERATIVE METHODS (good for GRMHD or entropy GRMHD or cold GRMHD)
    //
    /////////////
    retval+=set_guess_Wp(showmessages, lpflag, eomtype, prim, ptrgeom, &W_last, &Wp_last,wglobal, Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);


    /////////////////////////////////////////////////////////////////////////
    //
    //
    // HOT GRMHD INVERTERS (obtain Wp)
    //
    /////////////////////////////////////////////////////////////////////////
    if(eomtype==EOMGRMHD){

      // METHOD specific:
#if(WHICHHOTINVERTER==3)
      retval+= find_root_1D_gen_Eprime(showmessages, lpflag, eomtype, Wp_last, &Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
#elif(WHICHHOTINVERTER==2)
      retval+= find_root_2D_gen(showmessages, lpflag, eomtype, Wp_last, &Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
#elif(WHICHHOTINVERTER==1)
      retval+= find_root_1D_gen(showmessages, lpflag, eomtype, Wp_last, &Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
#else
      dualfprintf(fail_file,"No such WHICHHOTINVERTER=%d\n",WHICHHOTINVERTER);
      myexit(89285221);
#endif


#if(0)
      retval = find_root_2D_gen(showmessages, lpflag, eomtype, Wp_last, &Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
      retvalother = find_root_1D_gen(showmessages, lpflag, eomtype, Wp_last, &Wpother, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
      retvalother2=find_root_1D_gen_scn(showmessages, lpflag, eomtype, W_last, &W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);

      //  Wp=Wpother;
      if(retval!=retvalother || (fabs(Wp-Wpother)/(Wp+Wpother)>1E-13)){
        dualfprintf(fail_file,"2D/1D: i=%d j=%d :: %d %d %d :: Wp_last=%21.15g  Wp=%21.15g Wpother=%21.15g\n",ptrgeom->i,ptrgeom->j,retval,retvalother,retvalother2,Wp_last,Wp,Wpother);
        dualfprintf(fail_file,"scn W_last=%g W=%g D=%g  Wp_last=%g Wp=%g\n",W_last,W,D,W_last-D,W-D);
        dualfprintf(fail_file,"rho0=%g u=%g p=%g gamma=%g gammasq=%g w=%g utsq=%g\n",rho0,u,p,gamma,gammasq,w,utsq);
      }
#endif // end comparing SCN and Jon method



    } // end if hot GRMHD
  

    /////////////////////////////////////////////////////////////////////////
    //
    //
    // COLD GRMHD INVERTERS (obtain Wp)
    //
    /////////////////////////////////////////////////////////////////////////
    else if(eomtype==EOMCOLDGRMHD){



#if(WHICHCOLDINVERTER==0)
      retval0 = find_root_1D_gen(showmessages, lpflag, eomtype, Wp_last, &Wp0, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);


      //  if(IFUTOPRIMFAIL(*lpflag)){
      //    dualfprintf(fail_file,"E'[W,v^2] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g D=%21.15g\n",Wp0,Qdotnp,Qdotn,D);
      //    *lpflag=0; // let E' try again
      //  }
      //  else{
      //    dualfprintf(fail_file,"E'[W',v^2] equation gave good answer: Wp=%21.15g  Qdotnp=%21.15g\n",Wp0,Qdotnp);
      //  }

      retval+=retval0;
      Wp=Wp0;




#elif(WHICHCOLDINVERTER==1)



      // Do inversion using E'[W'] equation
      retval1 = find_root_1D_gen_Eprime(showmessages, lpflag, eomtype, Wp_last, &Wp1, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);





      //  if(IFUTOPRIMFAIL(*lpflag)){
      //    dualfprintf(fail_file,"E'[W'] equation gave bad answer: Wp=%21.15g  Bsq=%21.15g Qdotnp=%21.15g Qdotn=%21.15g D=%21.15g\n",Wp,Bsq,Qdotnp,Qdotn,D);

      //    for(i=0;i<4;i++) Qcovp[i] = U[QCOV0+i] ;
      //    raise_g(Qcovp,ptrgeom->gcon,Qconp) ;  

      //    for(i=0;i<4;i++) dualfprintf(fail_file,"%d Qcovp=%21.15g Qconp=%21.15g\n",i,Qcovp[i],Qconp[i]);


      //    *lpflag=0; // let momentum try again
      //  }
      //  else{
      //    dualfprintf(fail_file,"Wp=%21.15g Bsq=%21.15g QdotB=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Qdotnp=%21.15g D=%21.15g\n",Wp,Bsq,QdotB,QdotBsq,Qtsq,Qdotnp,D,Sc);
      //    dualfprintf(fail_file,"E'[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp,Qdotnp);
      //  }

  
      retval+=retval1;
      Wp=Wp1;



#elif(WHICHCOLDINVERTER==2)
      // Do inversion using P^2[W'] equation


#if(CRAZYDEBUG&&DEBUGINDEX)
      if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
        dualfprintf(fail_file,"Wp_last=%21.15g\n",Wp_last);
      }
#endif

      retval2 = find_root_1D_gen_Psq(showmessages, lpflag, eomtype, Wp_last, &Wp2, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
      //  if(IFUTOPRIMFAIL(*lpflag)){
      //    dualfprintf(fail_file,"P^2[W'] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g\n",Wp2,Qdotnp,Qdotn);
      //    myexit(0);
      //  }
      //  else{
      //    dualfprintf(fail_file,"P^2[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp2,Qdotnp);
      //  }


      //  Wp = (Wp>Wp2) ? Wp2 : Wp;
      retval+=retval2;
      Wp=Wp2;



#elif(WHICHCOLDINVERTER==3)
      // Do inversion using P^\alpha[W'] equation


      retval3 = find_root_3D_gen_Palpha(showmessages, lpflag, eomtype, Wp_last, &Wp3, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);



      //  if(IFUTOPRIMFAIL(*lpflag)){
      //    dualfprintf(fail_file,"P^2[W'] equation gave bad answer: Wp=%21.15g Qdotnp=%21.15g Qdotn=%21.15g\n",Wp2,Qdotnp,Qdotn);
      //    myexit(0);
      //  }
      //  else{
      //    dualfprintf(fail_file,"P^2[W'] equation gave good answer: Wp=%21.15g Qdotnp=%21.15g\n",Wp2,Qdotnp);
      //  }


      //  Wp = (Wp>Wp3) ? Wp3 : Wp;
      retval+=retval3;
      Wp=Wp3;

#else
      dualfprintf(fail_file,"No such WHICHCOLDINVERTER=%d\n",WHICHCOLDINVERTER);
      myexit(89285224);
#endif




    } // end if eomtype==EOMCOLDGRMHD



    /////////////////////////////////////////////////////////////////////////
    //
    //
    // ENTROPY GRMHD INVERTERS (obtain Wp)
    //
    /////////////////////////////////////////////////////////////////////////
    else if(eomtype==EOMENTROPYGRMHD){



      // METHOD specific:
#if(WHICHENTROPYINVERTER==0)
      retval+= find_root_1D_gen_Sc(showmessages, lpflag, eomtype, Wp_last, &Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
#else
      dualfprintf(fail_file,"No such WHICHENTROPYINVERTER=%d\n",WHICHENTROPYINVERTER);
      myexit(89285225);
#endif



    } // end if-else structure over eomtype
  }// end if doing non-analytic methods





  /////////////////////
  //
  // Use root if no failure in finding root
  //
  /////////////////////

  // check on result of inversion even if failure (i.e. even if retval>0)
  // Can be done for Wdependentsolution==0, but ignore Wp
  check_on_inversion(showmessages, prim, U, ptrgeom, Wp, Qtcon, Bcon, Bcov,retval, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);


  /// check on W' and convert W' -> prim[]
  if(retval==0 && Wdependentsolution==1){
    //  Check if solution was found
    retval+=check_Wp(lpflag, eomtype, prim, U, ptrgeom, Wp_last, Wp, retval, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // should add in case retval!=0 before this call
    if(retval) return(retval);


    // find solution
    retval+=Wp2prim(showmessages, lpflag, eomtype,prim, pressure, U, ptrgeom, Wp, Qtcon, Bcon, Bcov, retval, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);
    if(retval) return(retval);
  }


  /* done! */
  return(retval) ;

}






///returns global quanities Bsq,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc
/// returns other local quantities
static int compute_setup_quantities(FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, FTYPE *Bsq,FTYPE *QdotB,FTYPE *QdotBsq,FTYPE *Qtsq,FTYPE *Qdotn,FTYPE *Qdotnp,FTYPE *D, FTYPE *Sc, int whicheos, FTYPE *EOSextra)
{
  int i,j;
  FTYPE Qcov[4],Qcon[4],Qcovp[4],Qconp[4],ncov[4],ncon[4],Qsq,Qtcov[4];
  FTYPE alpha,alphasq,Dfactor;



  //  for(i=0;i<NPR;i++) dualfprintf(fail_file,"U[%d]=%g\n",i,U[i]);
  //  exit(0);

  *D = U[RHO] ;

#if(DOENTROPY!=DONOENTROPY)
  *Sc = U[ENTROPY];
#endif

  // get B^i
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  // get B_\mu
  lower_vec(Bcon,ptrgeom,Bcov) ;

  // get B^2
  *Bsq = 0. ;
  for(i=1;i<4;i++) *Bsq += Bcon[i]*Bcov[i] ;

  // get \eta^\mu and \eta_\mu
  ncov_calc_fromlapse(ptrgeom->alphalapse,ncov) ;
  raise_vec(ncov,ptrgeom,ncon);

  // -P_\alpha (P_t has D term inside it)
  // Q_\nu = \alpha \{T^t_\nu\}  = -T.\eta
  
  // Q'_\nu = Q_\nu + \delta^t_\nu D
  for(i=0;i<4;i++) Qcovp[i] = U[QCOV0+i] ;
  raise_vec(Qcovp,ptrgeom,Qconp) ;  



  // P^2, which can be written without use of T^t_t
  //  Qsq = 0. ;
  //for(i=0;i<4;i++) Qsq += Qcovp[i]*Qconp[i] ;

  // P^2 = Qtilde^2 can be written without T^t_t
  // Notice U !=Qcon, but only a problem for time term that is not used
  // T^t_t not really needed (cancels in the end)
  //  DLOOP(j) Qtcov[j]=0.0;
  //  DSLOOP(j,i) Qtcov[j] = U[QCOV0+i]*(delta(i,j)+ncon[i]*ncov[j]) ;
  //  raise_vec(Qtcov,ptrgeom,Qtcon) ;

  // P^\alpha = Qtcon^\alpha
  DLOOPA(j) Qtcon[j]=0.0;
  SSLOOP(j,i) Qtcon[j] += U[QCOV0+i]*(ptrgeom->gcon[GIND(i,j)] + ncon[i]*ncon[j]) ;
  //  for(i=0;i<4;i++){
  //    for(j=1;j<4;j++){
  //      Qtcon[i] = U[QCOV0+j]*(ptrgeom->gcon[GIND(j,i)]+ncon[j]*ncon[i]);
  //    }
  //  }
  lower_vec(Qtcon,ptrgeom,Qtcov);
  *Qtsq = 0. ;
  for(i=0;i<4;i++) *Qtsq += Qtcov[i]*Qtcon[i] ;

#if(0)
  DLOOPA(j) dualfprintf(fail_file,"Qtcon[%d]=%21.15g Qtcov[%d]=%21.15g\n",j,Qtcon[j],j,Qtcov[j]);
#endif


  // Qcon/Qcov  = -\eta.T = \alpha T^t_\nu :: where T^t_t *includes* rest-mass
  // avoid using this Q if possible to avoid large errors when v or u is small
  // put back in the D
  for(i=0;i<4;i++){
    if(i==0) Qcov[i]=Qcovp[i] - (*D);
    else Qcov[i]=Qcovp[i];
  }
  raise_vec(Qcov,ptrgeom,Qcon) ;  


  


  //  dualfprintf(fail_file,"D=%g\n",D);
  //  for(i=0;i<4;i++) dualfprintf(fail_file,"Qcon[%d]=%g Qcov[%d]=%g\n",i,Qcon[i],i,Qcov[i]);
  //  exit(0);



  // P.B (P_t not used because B^t=0, so ok to use Qcov)
  *QdotB = 0. ;
  for(i=1;i<4;i++) *QdotB += Qcov[i]*Bcon[i] ;
  *QdotBsq = (*QdotB)*(*QdotB) ;

  // Energy
  // Qdotnp=-E'=-E+D automatically
#if(0)
  // This does not resolve metric term to have gravity taken into account in non-rel limit
  // Note that gravity will be de-resolved once $r\sim 10^{14}M$
  *Qdotnp = Qconp[0]*ncov[0] + (*D)*(1.0-1.0/fabs(ncov[TT])) ; // -Qdotn-W = -Qdotnp-Wp
#else
  ////////////////
  //
  // Get 1-1/|\eta_t| in way that avoids catastrophic cancellation in non-rel gravity limit
  //
  ///////////////
  // 1-1/|\eta_t| = 1-1/\alpha = (\alpha-1)/\alpha = (\alpha^2-1)/(\alpha(\alpha+1)) = (-gcovpert[TT] + \alpha^2 g_{ti} g^{ti})/(\alpha(\alpha+1))
  //
  // \beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
  //  betasqoalphasq = 0.0;
  //  SLOOPA(j) betasqoalphasq += (ptrgeom->gcov[GIND(TT,j)])*(ptrgeom->gcon[GIND(TT,j)]);

  alpha = ptrgeom->alphalapse;
  alphasq = alpha*alpha;

  Dfactor = (-ptrgeom->gcovpert[TT] + alphasq*(ptrgeom->betasqoalphasq))/(alphasq+alpha);

  // now Compute Qdotnp
  *Qdotnp = Qconp[0]*ncov[0] + (*D)*(Dfactor) ; // -Qdotn-W = -Qdotnp-Wp
#endif


  *Qdotn  =  Qcon[0]*ncov[0] ;


#if(0)
  // deprecated DEBUG:
  // used to check on old vs. new method
  *Qtsqnew=*Qtsq;
  //  dualfprintf(fail_file,"1: Qtsq=%21.15g\n",Qtsq);
  for(i=0;i<4;i++) Qcovorig[i] = U[QCOV0+i] ;
  Qcovorig[TT] -=D;
  raise_vec(Qcovorig,ptrgeom,Qconorig) ;
  Qdotnorig = Qconorig[0]*ncov[0] ;
  *Qsqorig = 0. ;
  for(i=0;i<4;i++) *Qsqorig += Qcovorig[i]*Qconorig[i] ;
  *Qtsqorig = *Qsqorig + Qdotnorig*Qdotnorig ;
  //  dualfprintf(fail_file,"2: Qtsq=%21.15g\n",*Qtsq);
  *Qtsq=*Qtsqorig;
#endif


  // P^2 (since D already removed from E, then cancelling terms are small)
  // so no problem with catastrophic cancellation since involves energy term
  //  *Qtsq = *Qsq + (*Qdotnp)*(*Qdotnp) ; // P^2 in eta frame


  //  dualfprintf(fail_file,"Qdotn=%g Qtsq=%g\n",*Qdotn,*Qtsq);
  //exit(0);


#if(CRAZYDEBUG)
  // debug Qsq
  // deprecated DEBUG:
  if(ptrgeom->i==0 && ptrgeom->j==31 && nstep==9 && steppart==2){
    Qsq = 0. ;
    for(i=0;i<4;i++) *Qsq += Qcov[i]*Qcon[i] ;
    for(i=0;i<4;i++) dualfprintf(fail_file,"Qcon[%d]=%21.15g Qcov[%d]=%21.15g Qsq=%21.15g Qtsq=%21.15g\n",i,Qcon[i],i,Qcov[i],Qsq,Qtsq);
  }
#endif


  return(0);

}


/// compute Wp from primitives
/// also get normalizations wglobal
static int compute_Wp(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *Wpabs, FTYPE *etaabs, FTYPE *gamma, FTYPE *gammasq, FTYPE *rho0, FTYPE *u, FTYPE *p, FTYPE *wmrho0, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE utsq;
  FTYPE Ss;
  FTYPE w;
  int i,j;
  int jj;
  FTYPE Bcon[NDIM],Bcov[NDIM],bsq;
  FTYPE nuabs;

  /* calculate W from last timestep and use 
     for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += ptrgeom->gcov[GIND(i,j)]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < MAXNEGUTSQ) ) { 
    if(debugfail>=2) dualfprintf(fail_file,"utsq was %21.15g now 0.0 @ i=%d j=%d k=%d\n",utsq,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    utsq = 0.0;
  }
  if(utsq >= 0. && utsq < UTSQ_TOO_BIG) {
    // then good, if nan then fails with next else
  }
  else{
    if( debugfail>=2 ){
      dualfprintf(fail_file,"Utoprim_new(): utsq < 0 or utsq>UTSQ_TOO_BIG in utoprim_jon attempt.  Resetting to utsq=0.0 from utsq = %21.15g @ i=%d j=%d k=%d\n", utsq,ptrgeom->i,ptrgeom->j,ptrgeom->k) ;
      for(i=1;i<4;i++)
        for(j=1;j<4;j++) dualfprintf(fail_file,"i=%d j=%d :: gcov=%21.15g primi=%21.15g primj=%21.15g\n",i,j,ptrgeom->gcov[GIND(i,j)],prim[UTCON1+i-1],prim[UTCON1+j-1]) ;
    }

    // then set some kind of guess since guess may just be bad and may still be able to find a solution
    utsq = 0.0;
    
    //    *lpflag= UTOPRIMFAILCONVGUESSUTSQ; // guess failure actually
    // check on v^2>1 failure
    //check_utsq_fail(-1E30, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);
    //return(1) ;
  }

  *gammasq = 1. + utsq ;
  *gamma  = sqrt(*gammasq);
 
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  *rho0 = (D) / (*gamma) ;
  *u = prim[UU] ;
  if(whicheos==KAZFULL){
    *u = MAX(0.0,*u); // near degeneracy, allow somewhat "hot" solution so guess doesn't give immediately bad Wp and utsq(Wp)
  }
  *p = pressure_rho0_u(whicheos,EOSextra,*rho0,*u) ; // p(rho0,u) just used for guess, while p(rho0,\chi) used to get solution since assume don't know \chi yet.  Even if had initial \chi, wouldn't be final \chi anyways.
  *wmrho0=(*u)+(*p);
  w = (*rho0) + (*wmrho0) ;




  // W'=W-D (D removed from W)
  *Wp_last = (D)*utsq/(1.0+(*gamma)) + ((*u) + (*p))*(*gammasq);
  *W_last = *Wp_last + (D);


  // need b^2 to normalize W since error in W limited by b^2
  // GODMARK: Computing b^2 this way is too expensive, just use B^2 since just used to make an error estimate
  // Misses u.B term
  //  bsq_calc(prim, ptrgeom, &bsq);
  Bcon[TT]=0.0; SLOOPA(jj) Bcon[jj]=prim[B1-1+jj];
  lower_vec(Bcon,ptrgeom,Bcov);
  bsq = 0.0; SLOOPA(jj) bsq += Bcon[jj]*Bcov[jj];
  bsq *= ((ptrgeom->alphalapse * ptrgeom->alphalapse))/(*gammasq);

  // Get normalizations
  *Wpabs = fabs((D)*utsq)/(1.0+fabs(*gamma)) + (fabs(*u)+fabs(*p))*fabs(*gammasq) + SQRTMINNUMREPRESENT;
  nuabs = fabs(*u) + fabs(*p) + fabs(bsq);
  *etaabs = fabs(*rho0) + nuabs; 
  wglobal[2]=GAMMASQCHECKRESID * fabs(*rho0);


  return(0);
}


/////////////////////////////////
///
/// SETUP ITERATIVE METHODS (good for GRMHD or entropy evolution or cold GRMHD)
/// Doesn't use Sc since initial guess for W can be found from last primitives alone
///
/////////////////////////////////
static int set_guess_Wp(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Ss;
  FTYPE gammasq,gamma,rho0,u,p,wmrho0;
  FTYPE Wpabs, etaabs;
  int i,j;
  FTYPE Wpcoldhd;
  int numattemptstofixguess;
  int jj;


  // get Wp, wglobal, etaabs, gamma, gammaabs
  compute_Wp(showmessages, lpflag, eomtype, prim, ptrgeom, W_last, Wp_last, &Wpabs, &etaabs, &gamma, &gammasq, &rho0, &u, &p, &wmrho0, wglobal, Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);




  //  dualfprintf(fail_file,"Wlast=%g Wplast=%g\n",*W_last,*Wp_last);
  
  //  *wglobal=w+bsq; // used to validate normalized version of W
  
  //////////////// RIGHT
  //SUPERGODMARK: what if in nonrel case rho0 is identically zero? looks like rel. case but it's not (need rhorel for non-rel case)
  if( JONGENNEWTSTUFF &&   gammasq * etaabs >= (wglobal[2])  ){
    // RIGHT

  
    ////////////////// WRONG
    //  if( JONGENNEWTSTUFF && (w+bsq)*gammasq >= GAMMASQCHECKRESID ) { // WRONG
    ///////////////// WRONG

    int gotcold=0;
    if(coldgrhd(lpflag, Qtsq, D, &Wpcoldhd)){
      Wpcoldhd=*Wp_last; // if failed to get coldGRHD Wp, then just use Wp_last
      gotcold=1;
    }
    
    //  *wglobal=MAX(MAX(Wpcoldhd,*Wp_last),Wpabs); // used to validate normalized version of W -- accurate for cold GRHD
    *wglobal= MAX(MAX(Wpcoldhd,*Wp_last),Wpabs); // used to validate normalized version of W -- accurate for cold GRHD
    //*wglobal=Wpabs;
    
    if(gotcold) *Wp_last=MAX(2.0*Wpcoldhd,*Wp_last); // use cold HD to upgrade Wp if needed.  Don't upgrade from cold if just set to Wp_last from guess.

    // fudge
    //    *Wp_last = *W_last;
  }
  else{
    *wglobal=Wpabs; // non-relativistically correct
  }

  // include Wp-like term with rest-mass so becomes W-like term
  wglobal[1]=fabs(rho0)+wglobal[0]+SQRTMINNUMREPRESENT;


  //  if((nstep==4)&&(ptrgeom->i==0)&&(ptrgeom->j==47)){
  //  dualfprintf(fail_file,"utsq=%g p=%g w=%g W_last=%g Wp_last=%g\n",utsq,p,w,*W_last,*Wp_last);
  //    exit(0);
  //  }

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"pre verify\n");
    dualfprintf(fail_file,"rho0=%21.15g u=%21.15g p=%21.15g\n",rho0,u,p);
    
    dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g gammasq=%21.15g w=%21.15g\n",*W_last,*Wp_last,gammasq,w);
  }
#endif



  verify_Wlast(u, p, ptrgeom, W_last, Wp_last, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // W_last and Wp_last and rest of variables are already pointers

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"post verify\n");
    
    dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g gammasq=%21.15g w=%21.15g\n",*W_last,*Wp_last,gammasq,w);
  }
#endif




     
  //  dualfprintf(fail_file,"MID: Wlast=%g Wplast=%g\n",*W_last,*Wp_last);


#define MAXNUMGUESSCHANGES (1000)


  // sometimes above gives invalid guess (Wp=0 or utsq<0) so fix
  numattemptstofixguess=0;
  FTYPE utsq,Ss0;
  while(1){
    // under all eomtype's, ensure utsq reasonable for guess so Newton's method starts at reasonable value around which to work from
    // check  utsq from this guess for W
    utsq=utsq_calc(*W_last, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems

    if(eomtype==EOMENTROPYGRMHD){
      // for Entropy method, must also ensure Ss is reasonable
      Ss=Ss_Wp_utsq(whicheos,EOSextra,*Wp_last,D,utsq);
      // see what initial entropy is.  Even if Ss is not a nan or inf, can be so small that Newton stepping will be wrong, so keep increasing Wp until actually hotter than older entropy.  This ensures start hotter than required and not bumping-up against limit of method.
      // Noticed that entropy inversion will think it succeeds (low errx and dx) but it really fails very catastrophically! (conserved quantities don't match) for this reason.  Typically would give very low internal energy result, when should have been much higher.  Doesn't do that anymore.
      // But, once that very low entropy occurs, then noticed a secondary bad effect.  The updated conserved quantity would then be VERY large, leading to a VERY large entropy and huge internal energy.  While this is avoided with no internal energy drop-outs, still a concern that very low entropy regions can cause problematic fluxes and give high entropy.  Probably need larger floor on entropy than SMALL in EOS.
      Ss0=compute_specificentropy_wmrho0(whicheos,EOSextra,rho0,wmrho0);
    }
    else{
      // just some dummy value so passes below conditional without affect
      Ss=SMALL;
      Ss0=0.0;
    }

#define FRACSs0 (0.5) // want to have initial entropy not too far below previous entropy value
    if(utsq>=0.0 && utsq==utsq && isfinite(utsq) && (Ss==Ss && isfinite(Ss) && Ss>=Ss0-FRACSs0*fabs(Ss0))){
      // if utsq=nan or inf, will fail to reach here
      // if Ss=nan or inf, will fail to reach here
      if(numattemptstofixguess>0) if(showmessages && debugfail>=3) dualfprintf(fail_file,"GOOD Initial guess #%d/%d [i=%d j=%d k=%d] for W=%21.15g Wp=%21.15g Wp/D=%21.15g gives bad utsq=%21.15g Ss=%21.15g D=%21.15g u=%21.15g p=%21.15g gamma=%21.15g Ss0=%21.15g\n",numattemptstofixguess,MAXNUMGUESSCHANGES,ptrgeom->i,ptrgeom->j,ptrgeom->k,*W_last,*Wp_last,*Wp_last/D,utsq,Ss,D,u,p,gamma,Ss0);
      break;
    }
    else{
#if(PRODUCTION==0)
      if(showmessages && debugfail>=3) dualfprintf(fail_file,"Initial guess #%d/%d [i=%d j=%d k=%d] for W=%21.15g Wp=%21.15g Wp/D=%21.15g gives bad utsq=%21.15g Ss=%21.15g D=%21.15g u=%21.15g p=%21.15g gamma=%21.15g Ss0=%21.15g\n",numattemptstofixguess,MAXNUMGUESSCHANGES,ptrgeom->i,ptrgeom->j,ptrgeom->k,*W_last,*Wp_last,*Wp_last/D,utsq,Ss,D,u,p,gamma,Ss0);
#endif


      PFTYPE coldpflag;
      FTYPE coldWp;
      coldgrhd(&coldpflag, Qtsq, D, &coldWp);
      if(*Wp_last<coldWp) *Wp_last=coldWp;


      if(numattemptstofixguess>0){
        // then cold fix was not enough:
        //#define FACTORSHIFT (10.0) // bit much
#define FACTORSHIFT (2.0)
        *Wp_last = MAX(MAX(fabs(*Wp_last)*FACTORSHIFT,NUMEPSILON*100.0*fabs(D)),KINDASMALL);
      }

      // set new W
      *W_last = *Wp_last + (D);
    }


    if(numattemptstofixguess>MAXNUMGUESSCHANGES){
      *lpflag= UTOPRIMFAILCONVGUESSUTSQ; // guess failure actually
      // check on v^2>1 failure
      return(1) ;
    }


    numattemptstofixguess++;

  }


  //  dualfprintf(fail_file,"NEW: Wlast=%g Wplast=%g\n",*W_last,*Wp_last);


  // DEBUG:
  //  dualfprintf(fail_file,"DEBUG: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",utsq,rho0,u,p,bsq,*Wp_last,wglobal[2]);
  //int pl; PALLLOOP(pl) dualfprintf(fail_file,"DEBUG2: prim[%d]=%21.15g\n",pl,prim[pl]);
  //  int jjm,kkm; DLOOP(jjm,kkm) dualfprintf(fail_file,"DEBUG3: gcov=%21.15g gcon=%21.15g\n",ptrgeom->gcov[GIND(jjm,kkm)],ptrgeom->gcon[GIND(jjm,kkm)]);
  //for(pl=FIRSTEOSGLOBAL;pl<=LASTEOSGLOBAL;pl++) dualfprintf(fail_file,"DEBUG3: EOSextra[%d]=%21.15g\n",pl,EOSextra[pl]);





  return(0);


}


////////////////////////////////////////
///
/// Good for GRMHD or entropy GRMHD or cold GRMHD
///
/// Make sure that W is large enough so that v^2 < 1 : 
///
/// GODMARK: apparently this is necessary for 1D method, otherwise it can't find solution (2D method seems to be ok)
///
/////////////////////////////////////// 
static int verify_Wlast(FTYPE u, FTYPE p, struct of_geom *ptrgeom, FTYPE *W_last, FTYPE *Wp_last, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  int i_increase;


  //  FTYPE Wpnewsmall=-3.0*MAX(MAX(MAX(MAX(MAX(fabs(D),fabs(Qdotn)),fabs(Bsq)),fabs(Qdotnp)),fabs(wglobal[0])),fabs(wglobal[1]));
  // allow this constraint on W_last and Wp_last so starts hotter than could be.

  if(*W_last<-(D) || *Wp_last<0){
    *W_last = -(D);
    *Wp_last = 0;
#if(DEBUGINDEX)
    if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"W_last=%21.15g Wp_last=%21.15g\n",*W_last,*Wp_last);
    }
#endif
  }


  if(1){
    i_increase = 0;
    while( vsqgtr1(*W_last,Bsq,QdotBsq,Qtsq) && (i_increase < 10) ) {
    
      //    dualfprintf(fail_file,"bumped up %d %d %g %g\n",ptrgeom->i,ptrgeom->j,*W_last,*Wp_last);


      *Wp_last*= 10.;
      *W_last = *Wp_last+(D);

      //    dualfprintf(fail_file,"bumped up %d %d\n",ptrgeom->i,ptrgeom->j);


      i_increase++;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"badval :  W = %21.15g, i_increase = %d \n", *W_last, i_increase); 
#endif
    }
#if(!OPTIMIZED)
    if( i_increase >= 10 ) { 
      dualfprintf(fail_file,"i_increase is too large, i_increase = %d , W = %21.15g \n", i_increase, *W_last);
    }
#endif
   

#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"u = %21.15g,  p = %21.15g,  Bsq = %21.15g, Qtsq = %21.15g \n",u,p,Bsq,Qtsq);
      dualfprintf(fail_file,"Bcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
      dualfprintf(fail_file,"Bcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcov[0],Bcov[1],Bcov[2],Bcov[3]);
      dualfprintf(fail_file,"Qcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcon[0],Qcon[1],Qcon[2],Qcon[3]);
      dualfprintf(fail_file,"Qcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcov[0],Qcov[1],Qcov[2],Qcov[3]);
      dualfprintf(fail_file,"call find_root\n") ;  
    
    }
#endif
  }


  return(0);

}




////////////////////////////////////////
///
///  Check if solution was found
///
/////////////////////////////////////// 
static int check_Wp(PFTYPE *lpflag, int eomtype, FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp_last, FTYPE Wp, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  FTYPE Wtest;

   


  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (Wp == FAIL_VAL) ) {
    if( debugfail>=2 ) {
      dualfprintf(fail_file, "Failed to find a prim. var. solution!! Wp_last=%21.15g Bsq=%21.15g QdotBsq=%21.15g Qdotn=%21.15g D=%21.15g Qtsq=%21.15g \n",Wp_last,Bsq,QdotBsq,Qdotn,D,Qtsq);
      dualfprintf(fail_file, "Utoprim_new_body(): bad newt failure, nstep,steppart :: t,i,j,k, p[0-%d], U[0-%d] = %ld %d :: %21.15g %d %d %d ", NPR, NPR, nstep,steppart,t, ptrgeom->i, ptrgeom->j,ptrgeom->k );  
      int pliter,pl;
      PINVERTLOOP(pliter,pl){
        dualfprintf(fail_file, "%21.15g ", prim[pl]);
      }
      PINVERTLOOP(pliter,pl){
        dualfprintf(fail_file, "%21.15g ", U[pl]);
      }
      dualfprintf(fail_file, "\n");
    }      
    *lpflag= retval+UTOPRIMFAILCONVRET+1;// related to UTOPRIMFAILCONVRET
    if(debugfail>=1) dualfprintf(fail_file,"Got here Wp_last=%21.15g Wp=%21.15g retval=%d\n",Wp_last,Wp,retval);
    return(retval);
  }
  else{
    // Note: Wp comparable to D in scale/dimension
    // Note: Wp/wglobal[1] comparable to \gamma^2
    Wtest=Wp; ///wglobal[1]; // normalize to old densities

    // GODMARK
    //      if(Wtest<=0. && fabs(Wtest)>1E-4){
    // Complex result if W+D<0 or Wp-D+D<0 or Wp<0
    FTYPE Wpnewsmall=-3.0*MAX(MAX(MAX(MAX(MAX(fabs(D),fabs(Qdotn)),fabs(Bsq)),fabs(Qdotnp)),fabs(wglobal[0])),fabs(wglobal[1]));
    //    if(Wtest<=-(D)){
    if(Wtest<Wpnewsmall){
      if(debugfail>=1){
        dualfprintf(fail_file,"Wtest1 failed :: Wp=%21.15g D=%21.15g wglobal=%21.15g\n",Wp,(D),wglobal[1]);
        dualfprintf(fail_file,"nstep=%ld steppart=%d :: t=%21.15g :: i=%d j=%d k=%d p=%d\n",nstep,steppart,t,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);
      }
      *lpflag= retval+UTOPRIMFAILCONVRET+1;// related to UTOPRIMFAILCONVRET
      // reset Wp
      //      Wp = -(D)+fabs((D))*0.01;
      Wp = Wpnewsmall;
      //
      return(retval) ;
    }
    else if(Wtest < Wpnewsmall || Wtest/wglobal[1] > GAMMASQ_TOO_BIG) {
      //      dualfprintf(fail_file,"Wtest2 failed Wp=%21.15g\n",Wp);



      //if(Wtest <= 0.) {
      if( debugfail>=2 ) {
        dualfprintf(fail_file,"Wtest2 failure %21.15g %21.15g %21.15g %21.15g\n",Wtest,D,wglobal[1],GAMMASQ_TOO_BIG) ;
        dualfprintf(fail_file, "Utoprim_new_body(): Wtest<Wpnewsmall or Wtest=toobig failure, t,i,j,k, p[0-7], U[0-7] = %21.15g %d %d %d ", t, ptrgeom->i, ptrgeom->j,ptrgeom->k );  
        int pliter,pl;
        PINVERTLOOP(pliter,pl){
          dualfprintf(fail_file, "%21.15g ", prim[pl]);
        }
        PINVERTLOOP(pliter,pl){
          if(isfinite(U[pl])) dualfprintf(fail_file, "%21.15g ", U[pl]);
          else dualfprintf(fail_file, "nan ");
        }
        dualfprintf(fail_file, "\n");
      }      

      *lpflag=  UTOPRIMFAILCONVW;
      return(retval) ;
    }
  }


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"(Wp,Wp_last,Bsq,Qtsq,QdotB,gammasq,Qdotn) %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
                Wp,Wp_last,
                Bsq,Qtsq,QdotB,gammasq,Qdotn) ;
    dualfprintf(fail_file,"done find_root\n") ; 
  }
  if( ltrace2 ) {
    dualfprintf(fail_file, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    
  }
#endif

  /*
    if( nstep == 3 ) { 
    if( ((ptrgeom->i==89) && (ptrgeom->j==88||ptrgeom->j==111)) || (ptrgeom->i==92&&(ptrgeom->j==86||ptrgeom->j==113)) 
    || (ptrgeom->i==92&&(ptrgeom->j==86||ptrgeom->j==113) ) || (ptrgeom->i==94&&(ptrgeom->j==84||ptrgeom->j==115) ) 
    || (ptrgeom->i==105&&(ptrgeom->j==84||ptrgeom->j==115) ) || (ptrgeom->i==107&&(ptrgeom->j==86||ptrgeom->j==113) ) 
    || (ptrgeom->i==110&&(ptrgeom->j==88||ptrgeom->j==111) ) ) {
    dualfprintf(fail_file, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    }
    }
  */

  return(0);


}


/// check on v^2>1 failure
static void check_utsq_fail(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

#if(DEBUGINDEX)
  if(steppart==1 && nstep==32 && ifileglobal==0 &&  jfileglobal==19){
    dualfprintf(fail_file,"Qtsq=%21.15g Bsq=%21.15g D=%21.15g QdotB=%21.15g Wp=%21.15g\n",Qtsq,Bsq,D,QdotB,Wp);
    myexit(12245);
  }
#endif

}


static void check_on_inversion(int showmessages, FTYPE *prim, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  int checki;

#if(CHECKONINVERSION)
  checki=0;
  strcpy(newtonstats->invpropertytext[checki++],"Qdotnp");
  strcpy(newtonstats->invpropertytext[checki++],"Qtsq");
  //  strcpy(newtonstats->invpropertytext[checki++],"Qtsqorig");
  strcpy(newtonstats->invpropertytext[checki++],"Bsq");
  strcpy(newtonstats->invpropertytext[checki++],"DD");
  strcpy(newtonstats->invpropertytext[checki++],"Sc");
  strcpy(newtonstats->invpropertytext[checki++],"QdotB");
  strcpy(newtonstats->invpropertytext[checki++],"WWp");
  strcpy(newtonstats->invpropertytext[checki++],"Qtcon1");
  strcpy(newtonstats->invpropertytext[checki++],"Qtcon2");
  strcpy(newtonstats->invpropertytext[checki++],"Qtcon3");
  strcpy(newtonstats->invpropertytext[checki++],"Bcon1");
  strcpy(newtonstats->invpropertytext[checki++],"Bcon2");
  strcpy(newtonstats->invpropertytext[checki++],"Bcon3");
  // then store quantities required to check up on inversion in mathematica
  checki=0;
  newtonstats->invproperty[checki++]=Qdotnp;
  newtonstats->invproperty[checki++]=Qtsq;
  //  newtonstats->invproperty[checki++]=Qtsqorig;
  //  newtonstats->invproperty[checki++]=Qtsqnew;
  newtonstats->invproperty[checki++]=Bsq;
  newtonstats->invproperty[checki++]=D;
  newtonstats->invproperty[checki++]=Sc;
  newtonstats->invproperty[checki++]=QdotB;
  newtonstats->invproperty[checki++]=Wp;
  newtonstats->invproperty[checki++]=Qtcon[1];
  newtonstats->invproperty[checki++]=Qtcon[2];
  newtonstats->invproperty[checki++]=Qtcon[3];
  newtonstats->invproperty[checki++]=Bcon[1];
  newtonstats->invproperty[checki++]=Bcon[2];
  newtonstats->invproperty[checki++]=Bcon[3];
  if(checki>NUMINVPROPERTY){
    dualfprintf(fail_file,"Not enough memory for invproperty: checki=%d NUMINVPROPERTY=%d\n",checki,NUMINVPROPERTY);
  }
#endif
}



static int Wp2prim(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE *prim, FTYPE *pressure, FTYPE *U, struct of_geom *ptrgeom, FTYPE Wp, FTYPE *Qtcon, FTYPE *Bcon, FTYPE *Bcov, int retval, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE W,vsq;
  FTYPE utsq;
  int i;
  FTYPE rho0;
  FTYPE gtmp,gamma,gammasq,w,wmrho0,p,u,tmpdiff;
  FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);
  FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq);
  int retvsq;
  FTYPE pr0[NPR];
  int pl,pl2;
  int foundnan;

#if(DEBUGINDEX)
  //  if(ifileglobal==1 && jfileglobal==10 && nstep==137 && steppart==3){
  //Wp = -0.513243850397847;
  //    Wp = -0.4043038500443759;
  //Wp = -0.0010560020881852734;
  //Wp =  0.001055708263718976;
  //dualfprintf(fail_file,"Testing\n");
  //}
#endif

#if(!OPTIMIZED)
  if(!isfinite(Wp)){
    PLOOPALLINVERT(pl) prim[pl]=pr0[pl]; // since might contaminate solution
    *lpflag=  UTOPRIMFAILCONVW;
    return(1);
  }
#endif





  ////////////
  //
  // Invert (trivially) field components
  //
  ////////////
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  ///////////
  //
  // set backup of primitive
  //
  ///////////
  PLOOPALLINVERT(pl) pr0[pl]=prim[pl];

  
  ////////////
  //
  // Set dependent quantity W[W']
  //
  ////////////
  W=Wp+D;




  ////////////
  //
  // Compute:
  // 1) Either v^2 or \tilde{u}^2
  // 2) \gamma
  // 3) \chi = u+p
  // 4) \rho_0
  // 5) u
  //
  ////////////

#if(PRIMFROMVSQ==1)
  // Calculate vsq
  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;

  retvsq=validate_vsq(vsq,&vsq, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);


  if( retvsq==1 || retvsq==2) {
    //(vsq >= 1.0) || (vsq<0.0) ) {
    if( debugfail>=2 ) { 
      dualfprintf(fail_file,"vsq failure:  vsq = %21.15g , W = %21.15g \n",vsq, W) ;
      dualfprintf(fail_file, "Utoprim_new_body(): utsq==bad failure, t,i,j,k, p[0-7], U[0-7] = %21.15g %d %d %d ", t, ptrgeom->i, ptrgeom->j,ptrgeom->k );  
      PINVERTLOOP(pliter,pl){
        dualfprintf(fail_file, "%21.15g ", prim[pl]);
      }
      PINVERTLOOP(pliter,pl){
        dualfprintf(fail_file, "%21.15g ", U[pl]);
      }
      dualfprintf(fail_file, "\n");
    }      

    // check on v^2>1 failure
    check_utsq_fail(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

    *lpflag=  UTOPRIMFAILCONVUTSQVERYBAD;
    return(retval) ;
  }
  else if(retvsq==3){
    // if( vsq > VSQ_TOO_BIG ) {
    
    // check on v^2>VSQ_TOO_BIG failure (not really a failure of convergence, but of confidence in the solution)
    check_utsq_fail(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);
    
    *lpflag=  UTOPRIMFAILCONVUTSQ;
    return(retval) ;
  }

  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  gammasq=gamma*gamma;
  // can compute utsq or gamma from Wp alone instead of using vsq
  utsq = gammasq*vsq;

  // GODMARK: can compute u+p directly from utsq, Wp, and D.
  wmrho0=wmrho0_compute_vsq(Wp,D,vsq,gamma,gammasq);


#elif(PRIMFROMVSQ==0)

  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

  gammasq=1.0+utsq;
  gamma=sqrt(gammasq);

  wmrho0=wmrho0_compute_utsq(Wp,D,utsq,gamma,gammasq);


#endif

  rho0 = D / gamma;

  //  w = W * (1. - vsq) ;
  // not used, but compute correctly anyways
  // w = \rho_0 + u + p
  w = rho0+wmrho0;

  p = pressure_wmrho0(whicheos,EOSextra,rho0, wmrho0) ;
  //  u = w - (rho0 + p) ;
  u = wmrho0-p;


  *pressure=p; // return pressure since may be non-trivial to obtain from p(rho0,u) or inconsistent with p(rho0,\chi) if computed otherwise

  //  if((nstep==4)&&(ptrgeom->i==0)&&(ptrgeom->j==47)){
  //  dualfprintf(fail_file,"Wp=%g W=%g D=%g vsq=%g gamma=%g wmrho0=%g p=%g u=%g\n",Wp,W,D,vsq,gamma,wmrho0,p,u);
  //    exit(0);
  //  }

  //  dualfprintf(fail_file,"Wp=%g W=%g D=%g vsq=%g gamma=%g wmrho0=%g p=%g u=%g\n",Wp,W,D,vsq,gamma,wmrho0,p,u);

  //  dualfprintf(fail_file,"Wp=%g W=%g D=%g utsq=%g gamma=%g wmrho0=%g p=%g u=%g\n",Wp,W,D,utsq,gamma,wmrho0,p,u);


  // GODMARK
  // fix checks for cold case
  if( (rho0 <= 0.) 
      || ((p < 0.)&&(eomtype!=EOMCOLDGRMHD))
      ) { 
    if(showmessages && debugfail>=2 ) {
      tmpdiff = w - rho0;
      dualfprintf(fail_file,
                  "rho or uu < 0 failure: rho,w,(w-rho),p,u  = %21.15g %21.15g %21.15g %21.15g %21.15g gamma,utsq = %21.15g %21.15g : nstep=%ld steppart=%d :: t=%21.15g :: i=%d j=%d k=%d p=%d\n",rho0,w,tmpdiff,p,u,gamma, utsq,nstep,steppart,t,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);
    }
    // Below changed to pressure since internal energy can be negative due to arbitrary energy per baryon offset
    if((rho0<=0.)&&(p>=0.)) *lpflag=  UTOPRIMFAILRHONEG;
    if((rho0>0.)&&(p<0.)) *lpflag= UTOPRIMFAILUNEG;
    if((rho0<=0.)&&(p<0.)) *lpflag= UTOPRIMFAILRHOUNEG;
    if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNNOTADJUSTED) return(retval) ; // else let assign -- used to check how bad failure is.
  }



  //////////////
  //
  // If didn't exit so far, then assume valid inversion and so assign primitives (i.e. prim[])
  // Assign:
  // 1) rho_0
  // 2) u (either from energy or entropy evolution)
  // 3) \tilde{u}^i[W']
  // 4) [Already inverted field]
  //
  //////////////

  prim[RHO] = rho0 ;



  if(eomtype==EOMCOLDGRMHD){
    prim[UU] = 0;
  }
  else{
    prim[UU] = u ;
  }


#if(DOENTROPY!=DONOENTROPY)
  // put result in both UU and ENTROPY for primitive
  prim[ENTROPY]=prim[UU];
#endif



  //  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn; // Qtcon already defined in good way

  // GODMARK: gamma should be computed from W and conserved quantities directly as \gamma[W] (see notes)
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;


  //  dualfprintf(fail_file,"%d %d : %g %g %g\n",ptrgeom->i,ptrgeom->j,vsq,W,Wp);
  //if(fabs(prim[U1])>1E-30)
  //  for(i=1;i<4;i++) dualfprintf(fail_file,"u[%d]= %g %g\n",i,Qtcon[i],prim[UTCON1+i-1]);
  //  dualfprintf(fail_file,"retval=%d\n",retval);




  ////////////////////
  //
  // make sure nan and inf's don't pass through
  //
  ////////////////////
  foundnan=0;
  PLOOPALLINVERT(pl){
    if(!isfinite(prim[pl])){
      foundnan=1;
    }
  }

  if(foundnan){
    *lpflag=  UTOPRIMFAILNANRESULT;
    
    if(debugfail>=2)  PLOOPALLINVERT(pl2) dualfprintf(fail_file,"prim is nan: prim[%d]=%21.15g -> pr0[%d]=%21.15g\n",pl2,prim[pl2],pl2,pr0[pl2]);

    // remove nan since might contaminate solution or cause super-slowdown
    PLOOPALLINVERT(pl) prim[pl]=pr0[pl];

    return(1);
  }


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file," rho final = %21.15g ,  u final = %21.15g \n", rho0, u);
  }
#endif



  return(0);


}




/// appears to give qualitatively similar, but still qualitatively different results than original phys.ffde.c code.  Gives more oscillatory results for torus DISKFIELD force-free field (i.e. set rho=u=0 but otherwise start with init.fishmon.c type setup)
static int forcefree_inversion(int cleaned, struct of_geom *ptrgeom, FTYPE W, FTYPE *Qtcon, FTYPE Bsq, FTYPE *Bcon, FTYPE *Bcov, FTYPE Qtsq, FTYPE *U, FTYPE *prim, FTYPE *gammaret)
{
  // force-free variables
  FTYPE Qtconclean[NDIM],Qtcovclean[NDIM],Qtsqclean;
  FTYPE alphasq,gammamax,vcon[NDIM],vcov[NDIM],QtdotB,vconff[NDIM],vcovff[NDIM];
  int i;
  FTYPE vsq,gamma;
  FTYPE realBsq;

  // note that doing this means B.B!=Bsq
  Bsq+=SMALL; // in case Bsq=0 identically


  alphasq=GAMMAMAX*GAMMAMAX/(GAMMAMAX*GAMMAMAX-1.0);


  prim[RHO]=prim[UU]=0.0;


  // first 3-velocity
  vcon[TT]=0.0; // Qtcon[TT]=0
  for(i=1;i<4;i++) vcon[i] = Qtcon[i]/Bsq ;

  // first project out parallel velocity from Qtcon
  // note that \tilde{Q}.B = Q.B
  QtdotB=0.0;
  for(i=1;i<4;i++) QtdotB += Qtcon[i]*Bcov[i];

  Qtconclean[TT]=0.0;
  if(cleaned==1){
    for(i=1;i<4;i++) Qtconclean[i] = Qtcon[i] - QtdotB*Bcon[i]/Bsq;

    // define new cleaned 3-velocity
    vconff[TT]=0.0;
    for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/Bsq;
  }
  else{
    for(i=1;i<4;i++) Qtconclean[i] = Qtcon[i];

    // define new non-cleaned 3-velocity
    vconff[TT]=0.0;
    for(i=1;i<4;i++) vconff[i] = (Qtconclean[i] + 0.0*QtdotB*Bcon[i]/fabs(W) )/(fabs(W) + fabs(Bsq));
    
  }


  // v_\alpha cleaned
  lower_vec(vconff,ptrgeom,vcovff) ;

  // v^2 cleaned
  vsq=0.0;
  for(i=1;i<4;i++) vsq+=vconff[i]*vcovff[i]; // vcon[TT]=0

  //  dualfprintf(fail_file,"vsq=%21.15g Bsq=%21.15g alphasq=%21.15g Qtsq=%21.15g\n",vsq,Bsq,alphasq,Qtsq);

  // now limit the 3-velocity to be physical (no change in direction, so QtdotB=0 will be preserved)
  if( (vsq<0.0) || (vsq>1.0/alphasq) ){
    // then need to limit Lorentz factor
    // dissipation is expected to occur without affecting force-free equations

    // define cleaned Qtsq
    lower_vec(Qtconclean,ptrgeom,Qtcovclean) ;
    Qtsqclean=0.0;
    for(i=1;i<4;i++) Qtsqclean+=Qtconclean[i]*Qtcovclean[i];
 
    // define new vconff (overwrites above) that is cleaned AND has dissipation
    // v^\alpha = P^\alpha/(B^2+\sqrt{\alpha^2 P^2})
    vconff[TT]=0.0;
    //    for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/(Bsq+sqrt(alphasq*Qtsqclean)) ;
    //realBsq=sqrt(Bsq+sqrt(alphasq*Qtsqclean));
    // Note that E^2 = P^2/B^2, and then use same limiting procedure as in phys.ffde.c
    realBsq=sqrt(Qtsqclean*alphasq);

    if(cleaned==1){
      for(i=1;i<4;i++) vconff[i] = Qtconclean[i]/realBsq ;
    }
    else{
      for(i=1;i<4;i++) vconff[i] = (Qtconclean[i] + 0.0*QtdotB*Bcon[i]/fabs(W) )/(fabs(W) + fabs(realBsq));
    }

    // v_\alpha cleaned and w/ dissipation
    lower_vec(vconff,ptrgeom,vcovff) ;

    // v^2 cleaned and w/ dissipation
    vsq=0.0;
    for(i=1;i<4;i++) vsq+=vconff[i]*vcovff[i]; // vcon[TT]=0

  }
 
  // notice that vsq<1 enforced now, so should be ok
  if( (vsq<-VSQ_TOO_BIG) || (vsq>1.0/alphasq + 100.0*NUMEPSILON) ){
    dualfprintf(fail_file,"Failed to limit vcon\n");
    if(vsq<-VSQ_TOO_BIG) vsq=-VSQ_TOO_BIG;
    else return(UTOPRIMFAILCONVRET+1);
  }

  // gamma cleaned and w/ dissipation if necessary
  // clearly limits accuracy of gamma -- unavoidable when using B^\alpha instead of b^\alpha as fundamental variable
  gamma=1.0/sqrt(1.0-vsq);

  // set relative 4-velocity
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma*vconff[i];
  
  // set field (not really used)
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  // DEBUG:
  //  PLOOPALLINVERT(i) dualfprintf(fail_file,"vsq=%21.15g prim[%d]=%21.15g\n",vsq,i,prim[i]);

  *gammaret=gamma;

  return(0);

}


/// returns v^2>=1
/// does NOT depend on EOS
/// W^3 (W+2B^2) - S^2 (2W+B^2) \le W^2 (P^2-B^4)
/// GODMARK: vsqgtr1 based upon energy equation information?  What about P^2 equation inversion?
static int vsqgtr1(FTYPE W,FTYPE Bsq,FTYPE QdotBsq, FTYPE Qtsq)
{

  return(W*W*W * ( W + 2.*Bsq ) - QdotBsq*(2.*W + Bsq)  <= W*W*(Qtsq-Bsq*Bsq) );

}






//////////////////////////////////////////////
//
// GENERAL FUNCTIONS USED BY INVERSION
//
//////////////////////////////////////////////


/// returns \tilde{u}^2 = \gamma^2 v^2
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE utsq_calc(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE utsqtop,utsqbottom,utsq;
  FTYPE S,SoW;

  S = QdotB;
  SoW = S/W;

  utsqtop = Qtsq + (Bsq+2.0*W)*SoW*SoW;

  utsqbottom = (Bsq*Bsq-Qtsq) + 2.0*Bsq*W + W*W - SoW*SoW*(Bsq+2.0*W);

  utsq = utsqtop/utsqbottom;

#if(DEBUGINDEX)
  //  if(ifileglobal>90 && CRAZYNEWCHECK){
  //  if(CRAZYNEWCHECK&&0){
  //    dualfprintf(fail_file,"utsqtop=%21.15g utsqbottom=%21.15g Qtsq=%21.15g\n",utsqtop,utsqbottom,Qtsq);
  //  }
#endif

  return(utsq);

}

/// returns \gamma^2
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE gammasq_calc_W(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE gammatop,gammabottom,gammasq;
  FTYPE S,SoW;

  S = QdotB;
  SoW = S/W;

  gammatop = Bsq+W;

  gammabottom = W*W+2.0*Bsq*W + (Bsq*Bsq-Qtsq)-2.0*SoW*SoW*W - Bsq*SoW*SoW;

  gammasq = gammatop*gammatop/gammabottom;

  return(gammasq);


}

/// returns \gamma
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE gamma_calc_W(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE gammasq,gamma;

  gammasq = gammasq_calc_W(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

  // should I check sign of gammasq?  GODMARK

  gamma=sqrt(gammasq);

  return(gamma);

}




/// evaluate v^2 (spatial, normalized velocity) from W = \gamma^2 w
/// does NOT depend on EOS
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE vsq_calc(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Wsq,Xsq,Ssq;
  FTYPE SoW;
 
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);
  // Ssq = QdotBsq / Bsq;

  SoW = QdotB/W;
 
  //return(  Ssq * ( 1./Wsq - 1./Xsq )  +  Qtsq / Xsq  ); 
  // return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );

  return(  (SoW*SoW*(2.0*W+Bsq)  +  Qtsq)/Xsq );

}

///  evaluate dv^2/dW 
/// does NOT depend on EOS
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE dvsq_dW(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE W3,X3,Ssq,Wsq,X;
 
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X3 = X*X*X;

  // if(fabs(Bsq)==0.0) Ssq=0.0;
  // else Ssq = QdotBsq / Bsq;

  //return( -2.*( Ssq * ( 1./W3 - 1./X3 )  +  Qtsq / X3 ) ); 
  // return( -2.*( W3*Qtsq + QdotBsq * ( 3*W*X + Bsq*Bsq ) ) / ( W3 * X3 )   );

  return( -2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3  )  ); // RIGHT (avoids catastrophic cancellation with W^3 term in numerator)

  // return( -2.*( Qtsq/X3  +  QdotBsq/Bsq * (1.0/W3 - 1.0/X3)  )  ); // RIGHT (said was WRONG!)


}


/// re-look at this GODMARK
/// Solution to the cold GRHD equation P^2 = ...[Wp]
static int coldgrhd(PFTYPE *lpflag, FTYPE Qtsq, FTYPE D, FTYPE *Wp)
{
  FTYPE det;

  det = D*D+Qtsq;

  if(det<0.0){
    //    dualfprintf(fail_file,"coldgrhd: no solution: det=%21.15g\n",det);
    *lpflag=UTOPRIMFAILCONVUTSQ;
    return(1);
  }
  else{
    if(D>0.0){
      *Wp = Qtsq/(D+sqrt(det));
    }
    else{
      *Wp = Qtsq/(D-sqrt(det));
    }
    return(0);
  }

  return(0); // shouldn't get here -- only used if coldgrhd() used
}


/// resultP = 0 = -P^2 + P^2 [ Wp ] for cold GRMHD
/// old version with poles
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE Psq_Wp_old(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE S,W,W2,DoW,SoW;
  FTYPE Y;

  W = Wp+D;
  Y = Bsq+2.0*W;
  W2 = W*W;
  S = QdotB;
  DoW = D/W;
  SoW = S/W;
 
#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"Qtsq=%21.15g Bsq=%21.15g D=%21.15g S=%21.15g\n",Qtsq,Bsq,D,S);
  }
#endif

  //  result = -Qtsq + (Bsq+2.0*D)*Bsq + 2.0*(Bsq+D)*Wp + Wp*Wp - (Bsq*DoW*DoW + SoW*SoW)*(Bsq+2.0*W);

  result = -Qtsq + Wp*(D+W)*(1.0 + Bsq*Y/W2) - Y*SoW*SoW;

  return(result);

}


/// old version with poles
/// result = d(resultP)/dWp = + dP^2[Wp]/dWp
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE dPsqdWp_Wp_old(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE S,W,DoW,SoW;

  W = Wp+D;
  S = QdotB;
  DoW = D/W;
  SoW = S/W;

  result = 2.0*(Bsq+W)*(Bsq/W*DoW*DoW + SoW*SoW/W + 1.0);

  return(result);

}

/// new version without poles
/// returns actually W^2 P^2
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static int Psq_Wp(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, FTYPE *resid, FTYPE *norm)
{
  FTYPE S,W,W2,DoW,SoW;
  FTYPE Y;

  W = Wp+D;
  Y = Bsq+2.0*W;
  W2 = W*W;
  S = QdotB;

  *resid = -Qtsq*W2 + Wp*(D+W)*(W2 + Bsq*Y) - Y*S*S;
  *norm = fabs(Qtsq*W2) + fabs(Wp*(D+W)*(W2)) + fabs(Wp*(D+W)*(Bsq*Y)) + fabs(Y*S*S);

  return(0);

}


/// new version without poles
/// result = d(W^2 P^2)/dWp = 
/// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE dPsqdWp_Wp(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE S,W,DoW,SoW;

  W = Wp+D;
  S = QdotB;

  result = 2.0*(D*( (Bsq+D)*(Bsq+D) - Qtsq) - (S*S) + ((Bsq + D)*(Bsq + 5.0*D) - Qtsq)*Wp + 3.0*(Bsq + 2.0*D)*(Wp*Wp) + 2.0*(Wp*Wp*Wp));

  return(result);

}









/// residual for Eprime (old normal residual) when \tilde{u}^2=0, which is minimum value before problems
#define MINRESIDEPRIME(Wp) ( Bsq*0.5+Qdotnp+3.0/5.0*Wp+(Bsq*Qtsq-QdotBsq)/(2.0*(Bsq+D+Wp)*(Bsq+D+Wp)) )

/// 0 = -E' + E'[Wp] for hot GRMHD (can be used for cold GRMHD too)
/// while this has poles near bad roots, the E\propto Wp for ultra relativistic case, so easy to find root
static int Eprime_Wp_unopt(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, FTYPE *resid, FTYPE *norm)
{
  FTYPE X,X2,W;
  FTYPE utsq;
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X2 = X*X;


  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems


  // PERFORMANCEMARK: Tested code with and without these type of if/else (all of them in utoprim_jon.c) and all rest used by func_Eprime_unopt() and found only changed performance from 7% of CPU to 6.90% CPU.  Not interesting change

  if(utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks if nan
    *resid=-VERYBIG; // force residual to be large and negative so indicates a problem
    *norm=-VERYBIG;
  }
  else{
    p_tmp=pressure_Wp_utsq(whicheos,EOSextra,Wp,D,utsq);
    
    // Note that -E'=Qdotnp
    *resid = Qdotnp + Wp - p_tmp + 0.5*Bsq + 0.5*(Bsq*Qtsq-QdotBsq)/X2;
    *norm = fabs(Qdotnp) + fabs(Wp) + fabs(p_tmp) + fabs(0.5*Bsq) + fabs(0.5*(Bsq*Qtsq)/X2) + fabs(0.5*(QdotBsq)/X2);

    if(p_tmp!=p_tmp || *resid!=*resid) *resid=-VERYBIG; // indicates still problem
  }

  return(0);

}




/// dE'/dWp [ Wp ]
/// It is ok that there exists a term (1-dpdW) since even with ideal gas EOS in non-rel limit (1-dpdWp) still order unity
/// as above, while poles near bad roots, linear dependence makes easier to find roots.
/// OPTMARK: Note that I already tried removing if-else and has no effect (05/20/09)
static FTYPE dEprimedWp_unopt(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE X,X3,W;
  FTYPE utsq;
  FTYPE dvsq,dp1,dp2,dpdWp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X3 = X*X*X;

  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems
  if(utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks for nan
    result=0.0; // indicates a problem
  }
  else{
    dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // no precision problem
    dp1 = dpdWp_calc_utsq(whicheos,EOSextra,Wp,D,utsq);
    dp2 = dpdvsq_calc2_Wp(whicheos,EOSextra, Wp, D, utsq );
    dpdWp = dp1  + dp2*dvsq; // dp/dW = dp/dWp
    
    result = (1.0 - dpdWp) - (Bsq*Qtsq - QdotBsq)/X3;

    // check for nan (but not inf)
    if(dvsq!=dvsq || dp1!=dp1 || dp2!=dp2 || dpdWp!=dpdWp || result!=result) result=0.0; // indicates still problem

  }

  return(result);

}


/// this version has no poles
/// assumes X2!=0, in which case previous old version gives inf anyways
/// 0 = (B^2+W)^2 (-E' + E'[Wp]) for hot GRMHD (can be used for cold GRMHD too)
static FTYPE Eprime_Wp_new1(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE X,X2,W;
  FTYPE utsq;
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X2 = X*X;
 
  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems
  if(utsq<UTSQNEGLIMIT || utsq!=utsq){
    result=0.0; // indicates problem
  }
  else{
    p_tmp=pressure_Wp_utsq(whicheos,EOSextra,Wp,D,utsq);
    
    result = X2*(Qdotnp + Wp - p_tmp + 0.5*Bsq) + 0.5*(Bsq*Qtsq-QdotBsq);

    if(p_tmp!=p_tmp || result!=result) result=0.0; // indicates still problem
  }

#if(DEBUGINDEX)
  //  if(ifileglobal>90 && CRAZYNEWCHECK){
  //  if(CRAZYNEWCHECK&&0){
  //    dualfprintf(fail_file,"Eprime: Wp=%21.15g W=%21.15g X2=%21.15g p_tmp=%21.15g utsq=%21.15g D=%21.15g result=%21.15g\n",Wp,W,X2,p_tmp,utsq,D,result);
  //  }
#endif

  return(result);

}


/// this version has no poles
/// d((B^2+W)^2 E')/dWp [ Wp ]
/// It is ok that there exists a term (1-dpdW) since even with ideal gas EOS in non-rel limit (1-dpdWp) still order unity
static FTYPE dEprimedWp_new1(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE X,X3,W;
  FTYPE utsq;
  FTYPE dvsq,dp1,dp2,dpdWp;
  FTYPE p_tmp;


  W = Wp+D;
  X = Bsq + W + SMALL;
  X3 = X*X*X;

  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems
  if(utsq<UTSQNEGLIMIT || utsq!=utsq){
    result=0.0; // indicates problem
  }
  else{

    dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // no precision problem
    dp1 = dpdWp_calc_utsq(whicheos,EOSextra,Wp,D,utsq);
    dp2 = dpdvsq_calc2_Wp(whicheos,EOSextra, Wp, D, utsq );
    dpdWp = dp1  + dp2*dvsq; // dp/dW = dp/dWp
    
    // just copied from Eprime_Wp() to avoid recomputing utsq_calc()
    p_tmp=pressure_Wp_utsq(whicheos,EOSextra,Wp,D,utsq); 
    
    result = X*(2.0*(Qdotnp + Wp - p_tmp + 0.5*Bsq) + (1.0 - dpdWp)*X );

    if(dvsq!=dvsq || dp1!=dp1 || dp2!=dp2 || dpdWp!=dpdWp || p_tmp!=p_tmp || result!=result) result=0.0; // indicates still problem


  }

#if(DEBUGINDEX)
  //  if(ifileglobal>90 && CRAZYNEWCHECK){
  //  if(CRAZYNEWCHECK&&0){
  //    dualfprintf(fail_file,"dEprime: Wp=%21.15g W=%21.15g X3=%21.15g p_tmp=%21.15g utsq=%21.15g dp1=%21.15g dp2=%21.15g dpdWp=%21.15g result=%21.15g\n",Wp,W,X3,p_tmp,utsq,dp1,dp2,dpdWp,result);
  //  }
#endif

  return(result);

}









/// 0 = -Sc + Sc[Wp] for entropy GRMHD
static FTYPE Sc_Wp_old(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE W;
  FTYPE utsq;
  FTYPE Ss_tmp;


  W = Wp+D;
  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems


  if(utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks if nan
    result=-VERYBIG; // force residual to be large and negative so indicates a problem
  }
  else{
    Ss_tmp=Ss_Wp_utsq(whicheos,EOSextra,Wp,D,utsq);
    
    result = -Sc + D*Ss_tmp;

    if(Ss_tmp!=Ss_tmp || result!=result) result=-VERYBIG; // indicates still problem
  }

  // DEBUG:
  //  dualfprintf(fail_file,"Sc_Wp_old: Wp=%21.15g W=%21.15g D=%21.15g utsq=%21.15g :: Ss_tmp=%21.15g result=%21.15g Sc=%21.15g\n",Wp,W,D,utsq,Ss_tmp,result,Sc);

  return(result);

}

/// 0 = Wp(-Sc + Sc[Wp]) for entropy GRMHD
/// Noticed that original residual behaves like log(Wp) near the root, so that the derivative shoots too far.  Multiplying by Wp solves that problem.
static int Sc_Wp_unopt(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, FTYPE *resid, FTYPE *norm)
{
  FTYPE W;
  FTYPE utsq;
  FTYPE Ss_tmp;


  W = Wp+D;
  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems


  if(utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks if nan
    *resid=-VERYBIG; // force residual to be large and negative so indicates a problem
    *norm=-VERYBIG; // force residual to be large and negative so indicates a problem
  }
  else{
    Ss_tmp=Ss_Wp_utsq(whicheos,EOSextra,Wp,D,utsq);
    
    *resid = Wp*(-Sc + D*Ss_tmp);
    *norm = fabs(Wp)*fabs(Sc) + fabs(Wp)*fabs(D*Ss_tmp);

    if(Ss_tmp!=Ss_tmp || *resid!=*resid) *resid=-VERYBIG; // indicates still problem
  }

  // DEBUG:
  //  dualfprintf(fail_file,"Sc_Wp_unopt: Wp=%21.15g W=%21.15g D=%21.15g utsq=%21.15g :: Ss_tmp=%21.15g *resid=%21.15g Sc=%21.15g\n",Wp,W,D,utsq,Ss_tmp,*resid,Sc);

  return(0);

}




/// dSc/dWp [ Wp ] for entropy GRMHD
static FTYPE dScdWp_unopt(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE result;
  FTYPE W;
  FTYPE utsq;
  FTYPE dvsq,dSs1,dSs2,dSsdWp;
  FTYPE Ss;

  W = Wp+D;
  utsq=utsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems
  if(utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks for nan
    result=0.0; // indicates a problem
  }
  else{
    dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // no precision problem
    dSs1 = dSsdWp_calc_utsq(whicheos,EOSextra,Wp,D,utsq);
    dSs2 = dSsdvsq_calc2_Wp(whicheos,EOSextra, Wp, D, utsq );
    dSsdWp = dSs1  + dSs2*dvsq; // dSs/dW = dSs/dWp
    Ss=Ss_Wp_utsq(whicheos,EOSextra,Wp,D,utsq);
    
    // old result for old residual
    //   result = D*dSsdWp;
    
    // new result for new residual
    result = -Sc + D*(dSsdWp*Wp + Ss);

    // check for nan (but not inf)
    if(dvsq!=dvsq || dSs1!=dSs1 || dSs2!=dSs2 || dSsdWp!=dSsdWp || result!=result) result=0.0; // indicates still problem

  }

  // DEBUG:
  //  dualfprintf(fail_file,"dScdWp:  Wp=%21.15g W=%21.15g D=%21.15g utsq=%21.15g :: dvsq=%21.15g dSs1=%21.15g dSs2=%21.15g dSsdWp=%21.15g result=%21.15g\n",Wp,W,D,utsq,dvsq,dSs1,dSs2,dSsdWp,result);


  return(result);

}













/********************************************************************

  x1_of_x0(): Only Used For 2D Inversion
           
    -- calculates x1 from x0 depending on value of vartype;

                          x0   x1
             vartype = 2 : Wp, gamma
             vartype = 3 : Wp, vsq
             vartype = 4 : Wp, utsq
    
    -- asumes x0 is already physical

*********************************************************************/

static int x1_of_x0(FTYPE x0, FTYPE *x1, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra) 
{
  FTYPE vsq;
  //  FTYPE dv = 1.e-15;
  FTYPE W;
  int retval;



  W = x0+D;
  
#if(CHANGEDTOOLDER)
  vsq = fabs(vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra)) ; // guaranteed to be positive (>=0)
#else
  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ; // ok if not positive
#endif

  retval=validate_vsq(vsq,&vsq, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

  //  return( ( vsq > 1. ) ? VSQ_TOO_BIG : vsq   ); 

  *x1=vsq;
  return(retval); // indicates failure of v^2 or not

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon their definitions:

    // added wglobal to properly normalize -- global variable set before find_root_2d()

    // avoided truncating bottom value of x[0] and x[1] since then can't converge properly to slightly unphysical values or temporarily reach out into unphysical space and come back in.

    // avoid truncating x[1] for v^2<0 since all equations are ok for v^2<0 and can help convergence to nearly physical or even physical solutions

    // must avoid x[0] too large and v^2>1, however, since that produces nan's.
    
    *********************************************************************/


/// Wpold is previous Newton iteration's value
/// *Wpnew is updated Newton iterations value
static void validate_Wp(FTYPE Wpold, FTYPE *Wpnew, FTYPE *wglobal, FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Wp;
  FTYPE dv = NUMEPSILON;
  PFTYPE vlpflag;

  /* Always take the absolute value of Wpnew[0] and check to see if it's too big:  */ 
  //  Wpnew[0] = fabs(Wpold);
  // below allows cases when can be negative internal energy density.  Primarily allowed by wglobal[0].
  FTYPE Wpnewsmall=-3.0*MAX(MAX(MAX(MAX(MAX(fabs(D),fabs(Qdotn)),fabs(Bsq)),fabs(Qdotnp)),fabs(wglobal[0])),fabs(wglobal[1]));
  //  dualfprintf(fail_file,"Wpnewsmall=%Lg : %Lg %Lg %Lg %Lg %Lg %Lg\n",Wpnewsmall,D,Qtsq,Qdotn,*wglobal,Bsq,Qdotnp);
  if(Wpnew[0]<Wpnewsmall){
    //  if(Wpnew[0]<-D && D>0.0){
    //    Wpnew[0] = 10.0*Qtsq/(2.0*D+Wp0);
    //    Wpnew[0]=10.0*D+dv;

    // reset above HD solution
    coldgrhd(&vlpflag, Qtsq, D, &Wp); // ignore vlpflag result
    Wpnew[0] = Wp*10.0;

#if(CRAZYDEBUG)
    dualfprintf(fail_file,"Wpnew[0] changed from %21.15g to %21.15g (D=%21.15g Qtsq=%21.15g)\n",Wpold,Wpnew[0],D,Qtsq);
#endif    

  }

  // revert to previous Wp if "new" Wp very bad
  // can't go back to previous Wp otherwise won't go anywhere
  // Wpnew[0]/global = gamma^2
  Wpnew[0] = (fabs(Wpnew[0]/wglobal[1]) > GAMMASQ_TOO_BIG) ?  0.5*(Wpold+Wpnew[0]) : Wpnew[0];


}


/// retvsq=0 : no problem
/// retvsq=1 : >=1.0
/// retvsq=2 : < too negative
/// retvsq=3 : >VSQ_TOO_BIG (forced failure even if converged)
static int validate_vsq(FTYPE vsqold,FTYPE *vsqnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  int retval;


  retval=0; // assume no problem

  if(vsqnew[0]>1.0){
    retval=1;
    vsqnew[0]=VSQ_TOO_BIG;
    return(retval);
  }


  if(vsqnew[0]>VSQ_TOO_BIG){
    retval=3;
    vsqnew[0]=VSQ_TOO_BIG;
    return(retval);
  }

#if(VSQNEGCHECK==1)
  // unnecessary limitation on Newton's method
  if( (vsqnew[0] < 0.) && (fabs(vsqnew[0]) < MAXNEGVSQ) ) { 
    vsqnew[0] = 0.0;
    retval=2;
    return(retval);
  }
#elif(VSQNEGCHECK==2)
  // best can do, but need to indicate that did it since Newton's method will fail to compute residual and derivative correctly
  if(vsqnew[0] <=-VSQ_TOO_BIG) { 
    vsqnew[0] = -VSQ_TOO_BIG;
    retval=2;
    return(retval);
  }
#elif(VSQNEGCHECK==3)
  if(vsqnew[0]<0.0) {
    vsqnew[0]=0.0;
    retval=2;
    return(retval);
  }
#endif


  return(retval);

}


/// validate functions here change where Newton is, so ok to change arbitrarily since next Newton step will be in well-defined place (although may get stuck)
static void validate_x_2d(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  validate_Wp(x0[0], &x[0], wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

  validate_vsq(x0[1], &x[1], wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

}


/// 1D methods using W'
static void validate_x_1d_old(FTYPE x[1], FTYPE x0[1], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE dv = NUMEPSILON;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  //  x[0] = fabs(x[0]);
  FTYPE Wpnewsmall=-3.0*MAX(MAX(MAX(MAX(MAX(fabs(D),fabs(Qdotn)),fabs(Bsq)),fabs(Qdotnp)),fabs(wglobal[0])),fabs(wglobal[1]));
  //  dualfprintf(fail_file,"Wpnewsmall=%Lg : %Lg %Lg %Lg %Lg %Lg %Lg\n",Wpnewsmall,D,Qtsq,Qdotn,*wglobal,Bsq,Qdotnp);
  if(x[0]<Wpnewsmall) x[0]=D+dv;
  // x[0]/global = gamma^2
  x[0] = (x[0]/wglobal[1] > GAMMASQ_TOO_BIG) ?  x0[0] : x[0];

}



/// 1D methods using W'
static void validate_x_1d(FTYPE x[1], FTYPE x0[1], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  validate_Wp(x0[0], &x[0], wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

}


/// set pointer to the correct validate_x function
static void pick_validate_x(int eomtype,void (**ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra))
{

  if(
     (WHICHHOTINVERTER==1 && (eomtype==EOMGRMHD)  ) ||
     (WHICHHOTINVERTER==3 && (eomtype==EOMGRMHD) ) ||
     ( (WHICHCOLDINVERTER==0 || WHICHCOLDINVERTER==1) && eomtype==EOMCOLDGRMHD) ||
     (WHICHCOLDINVERTER==2 && eomtype==EOMCOLDGRMHD)
     || (eomtype==EOMENTROPYGRMHD)
     ){
    *ptr_validate_x = &validate_x_1d;
  }
  else if(WHICHHOTINVERTER==2 && (eomtype==EOMGRMHD)){
    *ptr_validate_x = &validate_x_2d;
  }
  else {
    dualfprintf(fail_file,"pick_validate_x: Shouldn't reach here\n");
    myexit(999);
  }


}

static void validate_x(void (*ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra),FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra) 
{
  (*ptr_validate_x)(x,x0,wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

}



//////////////////////////////////////////////////////
//
//
//////////////  STUFF FOR 1D METHOD
//
//
//////////////////////////////////////////////////////




/**************************************************** 
*****************************************************/

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE vsq,vsqfake,gamma,gamma_sq,rho0,w,W,Wp,drdW,dpdW,Wsq,p_tmp, dvsq,dp1,dp2,a_term,ap_term ;
  FTYPE X,X2,X3;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  int retvalvalidatevsq;

  Wp = x[0];
  W = Wp+D;
  Wsq = W*W;
  X = Bsq+W;
  X2 = X*X;
  X3 = X2*X;

  //dualfprintf(fail_file,"call utsq in residual\n") ;
  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;
  //dualfprintf(fail_file,"done w/ utsq in residual\n") ;

  // GODMARK DANGEROUS 
  //  vsqfake = ( vsq < -1.e-10 ) ?  dv : fabs(vsq) ;
  //  vsqfake = ( vsq > 1. ) ?  (1.-dv) : vsq ;
 
  retvalvalidatevsq=validate_vsq(vsq,&vsqfake, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

#if(0)
  // no longer use below methods 
#if(OPTIMIZED)
  vsqfake = ( vsq < -MAXNEGVSQ ) ?  0.0 : fabs(vsq) ;
  vsqfake = ( vsq > 1. ) ?  VSQ_TOO_BIG : vsq ;
#elif(0)
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<=1.0) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif
#endif

  if(retvalvalidatevsq==0){

    //rho0 = D * sqrt(1. - vsq) ;
    //dualfprintf(fail_file,"rho0,D,gamma: %21.15g %21.15g %21.15g\n",rho0,D,gamma) ;
    //  w = W * (1. - vsq);
    
    p_tmp = pressure_Wp_vsq(whicheos,EOSextra,Wp,D,vsqfake); // needs to be correct
    
    
    
    // Same as previous residual, however, jacobian is calculated using full differentiation w.r.t. W  
    //
    // dp/dW' = dp/dW + dP/dv^2 dv^2/dW
    
    
    dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // no precision problem
    dp1 = dpdWp_calc_vsq(whicheos,EOSextra, Wp, D, vsq ); // vsq can be unphysical
    dp2 = dpdvsq_calc_Wp(whicheos,EOSextra, Wp, D, vsq ); // vsq must be physical
    dpdW = dp1  + dp2*dvsq; // dp/dW = dp/dWp
    
    
#if(0) // Scott
    
    resid[0] = 
      + Wp 
      + 0.5 * Bsq * ( 1. + vsq ) // any vsq
      - 0.5*QdotBsq/Wsq
      + Qdotnp
      - p_tmp;

    norm[0] = 
      + fabs(Wp)
      + fabs(0.5 * Bsq * ( 1. + vsq ))
      + fabs(0.5*QdotBsq/Wsq)
      + fabs(Qdotnp)
      + fabs(p_tmp);
    
    jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;
    
#elif(0)
    //  resid[0] = 
    //    + W*Wsq
    //    + 0.5 * Wsq * Bsq * ( 1. + vsq )
    //    - 0.5*QdotBsq
    //    - Qdotn*Wsq
    //    - p_tmp*Wsq;
    //
    //  jac[0][0] = drdW = W * ( 3.*W  +  Bsq * (1. + vsq + 0.5*W*dvsq)  - 2.*Qdotn - 2.*p_tmp - W*dpdW );
    
#elif(1)  // JON:
    
    resid[0] = Qdotnp + Wp - p_tmp + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;
    norm[0] = fabs(Qdotnp) + fabs(Wp) + fabs(p_tmp) + fabs(0.5*Bsq) + fabs((Bsq*Qtsq)/X2) + fabs((QdotBsq)/X2);
    
    jac[0][0] = drdW = 1. - dpdW + (Bsq*Qtsq - QdotBsq)/X3*(-2.0);
    
#endif
    
    
    dx[0] = -resid[0]/drdW;
    
    *f = 0.5*resid[0]*resid[0];
    *df = -2. * (*f);
    
#if(!OPTIMIZED)
    if( ltrace ) {
      a_term = 0.5*Bsq*(1.+vsq) + Qdotnp - p_tmp;
      ap_term = dvsq * ( 0.5*Bsq - dp2 ) - dp1;
      
      dualfprintf(fail_file,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
      dualfprintf(fail_file,"func_1d_orig(): dvsq = %21.15g,  dp = %21.15g , dp1 = %21.15g , dp2 = %21.15g \n", dvsq, dpdW,dp1,dp2);
      dualfprintf(fail_file,"func_1d_orig(): W = %21.15g , vsq = %21.15g , a = %21.15g,  a_prime = %21.15g  \n", x[0], vsq, a_term, ap_term);
      
    }
#endif
  }
  else{ // indicate failure for this 1-D method since depends on v^2 or \tilde{u}^2 to be well-defined as a function of W'
    resid[0]=-VERYBIG;
    norm[0]=-VERYBIG;
    jac[0][0]=0.0;
    dx[0]=0.0;
    *f=-VERYBIG;
    *df=-VERYBIG;
  }

}





/**************************************************** 
  Routine for line searching for the 1D method . 
*****************************************************/

static FTYPE res_sq_1d_orig(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE vsq,W,Wp,Wsq,p_tmp,resid[1],norm[1];
  FTYPE vsqfake;
  FTYPE X,X2;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  int retvalvalidatevsq;
  

  Wp = x[0];
  W = Wp+D;
  Wsq = W*W;
  X = Bsq+W;
  X2 = X*X;

  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;


  retvalvalidatevsq=validate_vsq(vsq,&vsqfake, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

  if(retvalvalidatevsq==0){

#if(0)
    // GODMARK DANGEROUS 
#if(OPTIMIZED)
    // scn's version
    //  vsqfake = ( vsq < -dv ) ?  dv : fabs(vsq) ;
    //vsqfake = ( vsq > 1. ) ?  (1.-dv) : vsq ;
    // jon's version
    vsqfake = ( vsq < 0.0 ) ?  0.0 : vsq ;
    vsqfake = ( vsq > VSQ_TOO_BIG ) ? VSQ_TOO_BIG  : vsq ; // GODMARK: must let iterate if vsq>1.  Can't just force vsq or else Newton's method won't work.
#else
    if(vsq<0.0){
      if(vsq>-MAXNEGVSQ) vsqfake=0.0;
      else{
        dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
        myexit(1);
      }
    }
    if(vsq>VSQ_TOO_BIG){
      if(vsq<=1.0) vsqfake=VSQ_TOO_BIG;
      else{
        dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
        myexit(1);
      }
    }
#endif
#endif

    p_tmp = pressure_Wp_vsq(whicheos,EOSextra,Wp,D,vsqfake); // vsq<=1 for pressure



#if(0) // Scott
    resid[0] = 
      + Wp 
      + 0.5 * Bsq * ( 1. + vsq ) // residual must use unphysical vsq
      - 0.5*QdotBsq/Wsq
      + Qdotnp
      - p_tmp;

    norm[0] = 
      + fabs(Wp)
      + fabs(0.5 * Bsq * ( 1. + vsq ))
      + fabs(0.5*QdotBsq/Wsq)
      + fabs(Qdotnp)
      + fabs(p_tmp);

#elif(1)  // JON:

    resid[0] = Qdotnp + Wp - p_tmp + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;
    norm[0] = fabs(Qdotnp) + fabs(Wp) + fabs(p_tmp) + fabs(0.5*Bsq) + fabs((Bsq*Qtsq)/X2) + fabs((QdotBsq)/X2);

#endif
  }
  else{
    resid[0]=-VERYBIG;
    norm[0]=-VERYBIG;
  }


  return(  0.5*resid[0]*resid[0] );


}






/////////////////////////////////
//
// SCN versions of 1D method
//
/////////////////////////////////

static void func_1d_orig_scn(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE vsq,vsqfake,gamma,gamma_sq,rho0,w,W,drdW,dpdW,Wsq,p_tmp, dvsq,dp1,dp2,a_term,ap_term ;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  //dualfprintf(fail_file,"call utsq in residual\n") ;
  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;
  //dualfprintf(fail_file,"done w/ utsq in residual\n") ;

  // GODMARK DANGEROUS 
  //  vsq = ( vsq < -1.e-10 ) ?  dv : fabs(vsq) ;
  //  vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  
#if(OPTIMIZED)
  vsqfake = ( vsq < -MAXNEGVSQ ) ?  0.0 : fabs(vsq) ;
  vsqfake = ( vsq > 1. ) ?  VSQ_TOO_BIG : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<=1.0) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  //rho0 = D * sqrt(1. - vsq) ;
  //dualfprintf(fail_file,"rho0,D,gamma: %21.15g %21.15g %21.15g\n",rho0,D,gamma) ;
  //  w = W * (1. - vsq);

  p_tmp = pressure_W_vsq_scn(whicheos,EOSextra, W,D,vsqfake);



  // Same as previous residual, however, jacobian is calculated using full differentiation w.r.t. W  

  dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);
  dp1 = dpdW_calc_vsq_scn(whicheos,EOSextra, W, D, vsq ); // use real vsq
  dp2 = dpdvsq_calc_scn(whicheos,EOSextra, W, D, vsqfake );
  dpdW = dp1  + dp2*dvsq;

  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq ) // use real vsq
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  norm[0] = 
    + fabs(W )
    + fabs(0.5 * Bsq * ( 1. + vsq ))
    + fabs(0.5*QdotBsq/Wsq)
    + fabs(Qdotn)
    + fabs(p_tmp);

  jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;

  //  resid[0] = 
  //    + W*Wsq
  //    + 0.5 * Wsq * Bsq * ( 1. + vsq )
  //    - 0.5*QdotBsq
  //    - Qdotn*Wsq
  //    - p_tmp*Wsq;
  //
  //  jac[0][0] = drdW = W * ( 3.*W  +  Bsq * (1. + vsq + 0.5*W*dvsq)  - 2.*Qdotn - 2.*p_tmp - W*dpdW );

  dx[0] = -resid[0]/drdW;

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

#if(!OPTIMIZED)
  if( ltrace ) {
    a_term = 0.5*Bsq*(1.+vsq) + Qdotn - p_tmp;
    ap_term = dvsq * ( 0.5*Bsq - dp2 ) - dp1;

    dualfprintf(fail_file,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
    dualfprintf(fail_file,"func_1d_orig(): dvsq = %21.15g,  dp = %21.15g , dp1 = %21.15g , dp2 = %21.15g \n", dvsq, dpdW,dp1,dp2);
    dualfprintf(fail_file,"func_1d_orig(): W = %21.15g , vsq = %21.15g , a = %21.15g,  a_prime = %21.15g  \n", x[0], vsq, a_term, ap_term);
    
  }
#endif

}

/**************************************************** 
  Routine for line searching for the 1D method . 
*****************************************************/

static FTYPE res_sq_1d_orig_scn(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE vsq,vsqfake,W,Wsq,p_tmp,resid[1],norm[1];
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  vsq = vsq_calc(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;
  
  // GODMARK DANGEROUS 
#if(OPTIMIZED)
  // scn's version
  //  vsq = ( vsq < -dv ) ?  dv : fabs(vsq) ;
  //vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  // jon's version
  vsqfake = ( vsq < 0.0 ) ?  0.0 : vsq ;
  vsqfake = ( vsq > VSQ_TOO_BIG ) ? VSQ_TOO_BIG  : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsqfake=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<=1.0) vsqfake=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  p_tmp = pressure_W_vsq_scn(whicheos,EOSextra,W,D,vsqfake);


  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq )
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  norm[0] = 
    + fabs(W)
    + fabs(0.5 * Bsq * ( 1. + vsq ))
    + fabs(0.5*QdotBsq/Wsq)
    + fabs(Qdotn)
    + fabs(p_tmp);

  return(  0.5*resid[0]*resid[0] );


}







/// unoptimized vresion of func_Eprime
static void func_Eprime_unopt(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE drdW;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  FTYPE Wp;
  

  Wp = x[0];

  Eprime_Wp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid[0], &norm[0]);
  jac[0][0] = drdW = dEprimedWp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);


  dx[0] = -resid[0]/drdW;
  //  if((fabs(dx[0])/(Wp+D)>0.1){ // don't trust such a big jump
  // dx[0] = sign(dx[0])*0.1;
  // }
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


}


/// res_sq_Eprime() (or any pure residual function) is only called by line searching method
static FTYPE res_sq_Eprime_unopt(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Wp;
  FTYPE resid,norm;

  Wp = x[0];
  Eprime_Wp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid, &norm);

  return(  0.5*resid*resid );


}




/// get first part of kinetic stuff for func() or resid_sq()
/// for optimized version of func()
static void get_kinetics_part1(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra,  FTYPE *W, FTYPE *X, FTYPE *X2, FTYPE *X3, FTYPE *utsq, FTYPE *gamma, FTYPE *gammasq, FTYPE *rho0, FTYPE *wmrho0)
{

  /////////////
  //
  // compute common quantities
  //
  /////////////


  *W = Wp+D;
  *X = Bsq + *W + SMALL;
  *X2 = (*X)  *(*X);
  *X3 = (*X2) *(*X);

  *utsq=utsq_calc(*W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // check for precision problems

  // get \gamma^2, \gamma, \chi, and \rho_0 (no EOS used):
  *gammasq=1.0+ (*utsq);
  *gamma=sqrt(*gammasq);
  *wmrho0=wmrho0_compute_utsq(Wp, D, *utsq, *gamma, *gammasq);
  *rho0=D/(*gamma);



}

/// get kinetic stuff for func() or resid_sq()
/// for optimized version of func()
static void get_kinetics_part2(FTYPE Wp, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra,  FTYPE W, FTYPE X, FTYPE X2, FTYPE X3, FTYPE utsq, FTYPE gamma, FTYPE gammasq, FTYPE rho0, FTYPE wmrho0,FTYPE *dvsq,FTYPE *dwmrho0dW,FTYPE *drho0dW,FTYPE *dwmrho0dvsq,FTYPE *drho0dvsq)
{

  // total derivative dv^2/dW'
  *dvsq = dvsq_dW(W, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra); // no precision problem -- and no EOS calls inside -- purely kinetic

  *dwmrho0dW = 1.0/gammasq; // holding utsq fixed
  *drho0dW = 0.0; // because \rho=D/\gamma and holding utsq fixed

  *dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp); // holding Wp fixed
  *drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma and holding Wp fixed


}





/// optimized version of func_Eprime()
/// optimized in sense that minimizes extra kinetic calculations AND EOS lookups
static void func_Eprime_opt(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  ////////////////
  //
  // Get some kinetic quantities
  //
  ////////////////
  FTYPE Wp;
  FTYPE X,X2,X3,W;
  FTYPE utsq;
  FTYPE gamma,gammasq;
  FTYPE rho0,wmrho0;
  Wp = x[0];
  get_kinetics_part1(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra,  &W, &X, &X2, &X3, &utsq, &gamma, &gammasq, &rho0, &wmrho0);


  ////////////////
  //
  // Now can get thermodynamical quantities and results from them
  //
  ////////////////
  // PERFORMANCEMARK: Tested code with and without these type of if/else (all of them in utoprim_jon.c) and all rest used by func_Eprime() and found only changed performance from 7% of CPU to 6.90% CPU.  Not interesting change
  FTYPE Eprime,dEprimedWp;
  if(gammasq<0.0 || utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks if nan
    Eprime=-VERYBIG; // force residual to be large and negative so indicates a problem
    dEprimedWp=0.0; // indicates a problem

    resid[0]=Eprime;
    norm[0]=Eprime;
    jac[0][0]=dEprimedWp;
    dx[0] = VERYBIG;
    *f = 0.5*(Eprime*Eprime);
    *df = -2. * (*f);

  }
  else{
    ////////////////
    //
    // Now can get pressure and other thermodynamical quantities
    // Optimize for any EOS table lookup by computing P(\rho_0,\chi), dP/d\rho_0 and dP/d\chi now.
    // Otherwise, Eprime_Wp() and dEprimeWp() (as previously coded) repeated these calls and made calls for each EOS term while now can make a pipelined call and obtain all 3 EOS values at once.
    //
    ////////////////
    FTYPE pofchi,dpdrho,dpdchi;
    getall_forinversion(whicheos,EOMGRMHD,CHIDIFF,EOSextra,rho0,wmrho0,&pofchi,&dpdrho,&dpdchi);

    /////////////////
    //
    // get Eprime
    // Note that -E'=Qdotnp
    //
    /////////////////
    Eprime = Qdotnp + Wp - pofchi + 0.5*Bsq + 0.5*(Bsq*Qtsq-QdotBsq)/X2;
    FTYPE normEprime = fabs(Qdotnp) + fabs(Wp) + fabs(pofchi) + fabs(0.5*Bsq) + fabs(0.5*(Bsq*Qtsq)/X2) + fabs(0.5*(QdotBsq)/X2);


    if(pofchi!=pofchi || Eprime!=Eprime) Eprime=-VERYBIG; // indicates still problem

    //////////////////
    //
    // Get derivative-related kinetics
    //
    //////////////////
    FTYPE dvsq,dwmrho0dW,drho0dW;
    FTYPE dwmrho0dvsq,drho0dvsq;
    get_kinetics_part2(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra,  W, X, X2, X3, utsq, gamma, gammasq, rho0, wmrho0, &dvsq,&dwmrho0dW,&drho0dW,&dwmrho0dvsq,&drho0dvsq);

    ////////////////////////
    //
    // get mixed kinetic-thermal derivatives
    //
    ///////////////////////
    FTYPE dpdW,dpdvsq,dpdWp;
    // dvsq = dv^2/dW' total derivative
    dpdW   = drho0dW   *dpdrho  +     dwmrho0dW *dpdchi; // dP/dW' holding utsq fixed
    dpdvsq = drho0dvsq *dpdrho  +   dwmrho0dvsq *dpdchi; // dP/dvsq holding W' fixed where utsq and vsq are same fixed quantity
    dpdWp = dpdW  + dpdvsq*dvsq; // dp/dW = dp/dWp [total derivative]

    ////////////////////////
    //
    // get dEprimedWp
    //
    ///////////////////////
    dEprimedWp = (1.0 - dpdWp) - (Bsq*Qtsq - QdotBsq)/X3;


    // check for nan (but not inf)
    if(dvsq!=dvsq || dpdW!=dpdW || dpdvsq!=dpdvsq || dpdWp!=dpdWp || dEprimedWp!=dEprimedWp) dEprimedWp=0.0; // indicates still problem






    ////////////////////////
    //
    // residual, jacobian, and step
    //
    ///////////////////////
    resid[0]=Eprime;
    norm[0]=normEprime;
    jac[0][0]=dEprimedWp;
    dx[0] = -Eprime/dEprimedWp;
    *f = 0.5*(Eprime*Eprime);
    *df = -2. * (*f);


    // DEBUG:
    //    if(x[0]+dx[0]<SMALL){
    //      dualfprintf(fail_file,"EPRIMEPROBLEM: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g :: %21.15g %21.15g %21.15g :: %21.15g %21.15g %21.15g %21.15g %21.15g :: %21.15g %21.15g %21.15g %21.15g %21.15g :: %21.15g %21.15g\n"
    //    ,Wp, Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, W, X, X2, X3, utsq, gamma, gammasq, rho0, wmrho0
    //    ,pofchi,dpdrho,dpdchi
    //    ,dvsq,dwmrho0dW,drho0dW,dwmrho0dvsq,drho0dvsq
    //    ,dpdW,dpdvsq,dpdWp,Eprime,dEprimedWp
    //            ,x[0],dx[0]);
    //      int extraloop;
    //      for(extraloop=FIRSTEOSGLOBAL;extraloop<=LASTEOSGLOBAL;extraloop++) dualfprintf(fail_file,"EPRIMEPROBELMEXTRA: %21.15g\n",EOSextra[extraloop]);
    //    }


  }





}


/// res_sq_Eprime() (or any pure residual function) is only called by line searching method
static FTYPE res_sq_Eprime_opt(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE jac[NEWT_DIM][NEWT_DIM]; // dummy variable, return value not used
  FTYPE f,df; // dummy, return value not used
  int n=1; // true for this residual
  FTYPE resid[NEWT_DIM];
  FTYPE norm[NEWT_DIM];

  // just call func_Eprime_opt() since more optimized than before -- could optimize further but don't often or at all use this function
  func_Eprime_opt(x, dx, resid, norm, jac, &f, &df, n, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra);

  return(  0.5*resid[0]*resid[0] );


}




/// not that this is well-optimized.  Psq_Wp() and dPsqdWp_Wp() don't use subfunctions like "_unopt" versions for hot and entropy inversions.  Also no EOS functions to call.
static void func_Psq(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE drdW,Wp;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  Wp = x[0];


  Psq_Wp(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid[0], &norm[0]);
  jac[0][0] = drdW = dPsqdWp_Wp(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);


  dx[0] = -resid[0]/drdW;
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


}


static FTYPE res_sq_Psq(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Wp;
  FTYPE resid,norm;

  Wp = x[0];
  Psq_Wp(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid, &norm);

  return(  0.5*resid*resid );


}



static void func_Sc_unopt(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE drdW;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;
  FTYPE Wp;


  Wp = x[0];

  Sc_Wp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid[0], &norm[0]);
  jac[0][0] = drdW = dScdWp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);


  dx[0] = -resid[0]/drdW;
  //  if((fabs(dx[0])/(Wp+D)>0.1){ // don't trust such a big jump
  // dx[0] = sign(dx[0])*0.1;
  // }
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);


}


static FTYPE res_sq_Sc_unopt(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE Wp;
  FTYPE resid,norm;

  Wp = x[0];
  Sc_Wp_unopt(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra, &resid, &norm);

  return(  0.5*resid*resid );


}


/// optimized version of func_Sc()
/// very similar to optimized version of func_Eprime() but different "getall_forinversion()" call and different residual and derivative
static void func_Sc_opt(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{


  ////////////////
  //
  // Get some kinetic quantities
  //
  ////////////////
  FTYPE Wp;
  FTYPE X,X2,X3,W;
  FTYPE utsq;
  FTYPE gamma,gammasq;
  FTYPE rho0,wmrho0;
  Wp = x[0];
  get_kinetics_part1(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra,  &W, &X, &X2, &X3, &utsq, &gamma, &gammasq, &rho0, &wmrho0);


  ////////////////
  //
  // Now can get thermodynamical quantities and results from them
  //
  ////////////////
  // PERFORMANCEMARK: Tested code with and without these type of if/else (all of them in utoprim_jon.c) and all rest used by func_Eprime() and found only changed performance from 7% of CPU to 6.90% CPU.  Not interesting change
  FTYPE Scprime,dScprimedWp;
  if(gammasq<0.0 || utsq<UTSQNEGLIMIT || utsq!=utsq){ // also checks if nan
    Scprime=-VERYBIG; // force residual to be large and negative so indicates a problem
    dScprimedWp=0.0; // indicates a problem

    resid[0]=Scprime;
    norm[0]=Scprime;
    jac[0][0]=dScprimedWp;
    dx[0] = -Scprime/dScprimedWp;
    *f = 0.5*(Scprime*Scprime);
    *df = -2. * (*f);

  }
  else{
    ////////////////
    //
    // Now can get pressure and other thermodynamical quantities
    // Optimize for any EOS table lookup by computing P(\rho_0,\chi), dP/d\rho_0 and dP/d\chi now.
    // Otherwise, Eprime_Wp() and dEprimeWp() (as previously coded) repeated these calls and made calls for each EOS term while now can make a pipelined call and obtain all 3 EOS values at once.
    //
    ////////////////
    FTYPE Ssofchi,dSsdrho,dSsdchi;
    getall_forinversion(whicheos,EOMENTROPYGRMHD,CHIDIFF,EOSextra,rho0,wmrho0,&Ssofchi,&dSsdrho,&dSsdchi);

    /////////////////
    //
    // get Scprime
    //
    /////////////////
    // old residual
    Scprime = -Sc + D*Ssofchi;
    FTYPE Scprimenorm = fabs(Sc) + fabs(D*Ssofchi);

    // new residual (problem near Wp=0 of course that could diverge to Wp=0 and be wrong solution. Was using this "new" residual to help reach solution faster, but in general doesn't help so much and hurts since sometimes with tabulated EOS reaches Wp=0 and results in bad inversion.)
    // Also, for tabulated EOS, new residual more complicated and error doesn't drop.  Even though old residual has many iterations like new residual, at least old residual has very small error.  This is probably because derivative is not exactly consistent with function in the tabulated EOS yet both are used in the derivative expression.
    //////Scprime = Wp*(-Sc + D*Ssofchi);


    if(Ssofchi!=Ssofchi || Scprime!=Scprime) Scprime=-VERYBIG; // indicates still problem

    //////////////////
    //
    // Get derivative-related kinetics
    //
    //////////////////
    FTYPE dvsq,dwmrho0dW,drho0dW;
    FTYPE dwmrho0dvsq,drho0dvsq;
    get_kinetics_part2(Wp, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra,  W, X, X2, X3, utsq, gamma, gammasq, rho0, wmrho0, &dvsq,&dwmrho0dW,&drho0dW,&dwmrho0dvsq,&drho0dvsq);

    ////////////////////////
    //
    // get mixed kinetic-thermal derivatives
    //
    ///////////////////////
    FTYPE dSsdW,dSsdvsq,dSsdWp;
    dSsdW =   drho0dW   *dSsdrho +   dwmrho0dW   *dSsdchi; // dSs/dW' holding utsq fixed
    dSsdvsq = drho0dvsq *dSsdrho +   dwmrho0dvsq *dSsdchi;
    dSsdWp = dSsdW  + dSsdvsq*dvsq; // dSs/dW = dSs/dWp [total derivative]

    ////////////////////////
    //
    // get dScprimeWp
    //
    ///////////////////////
    // old result for old residual
    dScprimedWp = D*dSsdWp;
    
    // new result for new residual
    ////////dScprimedWp = -Sc + D*(dSsdWp*Wp + Ssofchi);

    // check for nan (but not inf)
    if(dvsq!=dvsq || dSsdW!=dSsdW || dSsdvsq!=dSsdvsq || dSsdWp!=dSsdWp || dScprimedWp!=dScprimedWp) dScprimedWp=0.0; // indicates still problem


    ////////////////////////
    //
    // residual, jacobian, and step
    //
    ///////////////////////
    resid[0]=Scprime;
    norm[0]=Scprimenorm;
    jac[0][0]=dScprimedWp;
    dx[0] = -Scprime/dScprimedWp;
    *f = 0.5*(Scprime*Scprime);
    *df = -2. * (*f);

  }



}


/// optimized version of res_sq_Sc()
/// just calls func_Sc_opt() for now since more optimized than original res_sq_Sc() still.  Could optimize further, but don't require only residual unless doing line searching that don't currently do.
static FTYPE res_sq_Sc_opt(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE jac[NEWT_DIM][NEWT_DIM]; // dummy variable, return value not used
  FTYPE f,df; // dummy, return value not used
  int n=1; // true for this residual
  FTYPE resid[NEWT_DIM];
  FTYPE norm[NEWT_DIM];

  // just call func_Eprime_opt() since more optimized than before -- could optimize further but don't often or at all use this function
  func_Sc_opt(x, dx, resid, norm, jac, &f, &df, n, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc, whicheos, EOSextra);

  return(  0.5*resid[0]*resid[0] );


}




/////////////////////////////////////////////////////
//
//////////////// STUFF FOR 2D METHOD
//
/////////////////////////////////////////////////////




/*********************************************************
returns residual  of vsq-1
**********************************************************/

static void func_vsq(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  
  FTYPE  Wp,W, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;

  p_tmp  = pressure_Wp_vsq(whicheos,EOSextra, Wp, D, vsq );
  dPdW   = dpdWp_calc_vsq(whicheos,EOSextra, Wp, D, vsq );
  dPdvsq = dpdvsq_calc_Wp(whicheos,EOSextra, Wp, D, vsq );


#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
  }
#endif


  /* These expressions were calculated using Mathematica */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  dx[0] = (-Bsq/2. + dPdvsq)*(Qtsq - vsq*((Bsq+W)*(Bsq+W)) + 
                              (QdotBsq*(Bsq + 2*W))/Wsq) + 
    ((Bsq+W)*(Bsq+W))*(-Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) -Wp + p_tmp);

  dx[1] = -((-1 + dPdW - QdotBsq/(Wsq*W))*
            (Qtsq - vsq*((Bsq+W)*(Bsq+W)) + (QdotBsq*(Bsq + 2*W))/Wsq)) - 
    2*(vsq + QdotBsq/(Wsq*W))*(Bsq + W)*
    (-Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp);

  detJ = (Bsq + W)*((-1 + dPdW - QdotBsq/(Wsq*W))*(Bsq + W) + 
                    ((Bsq - 2*dPdvsq)*(QdotBsq + vsq*(Wsq*W)))/(Wsq*W));
  
  dx[0] /= -(detJ) ;
  dx[1] /= -(detJ) ;

  jac[0][0] = -2*(vsq + QdotBsq/(Wsq*W))*(Bsq + W);
  jac[0][1] = -(Bsq + W)*(Bsq + W);
  jac[1][0] = -1 - QdotBsq/(Wsq*W) + dPdW;
  jac[1][1] = -Bsq/2. + dPdvsq;
  resid[0]  = Qtsq - vsq*(Bsq + W)*(Bsq + W) + (QdotBsq*(Bsq + 2*W))/Wsq;
  resid[1]  = -Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp;

  norm[0]  = fabs(Qtsq) + fabs(vsq*(Bsq + W)*(Bsq + W)) + fabs((QdotBsq*(Bsq + 2*W))/Wsq);
  norm[1]  = fabs(Qdotnp) + fabs((Bsq*(1 + vsq))/2.) + fabs(QdotBsq/(2.*Wsq)) + fabs(Wp) + fabs(p_tmp);

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"func_vsq(): x[0] = %21.15g , x[1] = %21.15g , dx[0] = %21.15g , dx[1] = %21.15g , f = %21.15g \n", 
                x[0], x[1],dx[0], dx[1], *f);
    
  }
#endif

}


/**********************************************************/

static void func_vsq2(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE norm[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  
  FTYPE  W, Wp,vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE t11;
  FTYPE t16;
  FTYPE t18;
  FTYPE t2;
  FTYPE t21;
  FTYPE t23;
  FTYPE t24;
  FTYPE t25;
  FTYPE t3;
  FTYPE t35;
  FTYPE t36;
  FTYPE t4;
  FTYPE t40;
  FTYPE t9;

  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq(whicheos,EOSextra, Wp, D, vsq );
  dPdW   = dpdWp_calc_vsq(whicheos,EOSextra, Wp, D, vsq );
  dPdvsq = dpdvsq_calc_Wp(whicheos,EOSextra, Wp, D, vsq );

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"W=%21.15g D=%21.15g vsq=%21.15g\n",W,D,vsq);
    dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
  }
#endif

#if(0)
  dualfprintf(fail_file,"p_tmp=%21.15g dPdW=%21.15g dPdvsq=%21.15g\n",p_tmp,dPdW,dPdvsq);
#endif

  /* These expressions were calculated using Mathematica, but made into efficient code using Maple */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  t2 = -0.5*Bsq+dPdvsq;
  t3 = Bsq+W;
  t4 = t3*t3;
  t9 = 1/Wsq;
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9;
  FTYPE t11norm = fabs(Qtsq)+fabs(vsq*t4)+fabs(QdotBsq*(Bsq+2.0*W)*t9);
  t16 = QdotBsq*t9;
  t18 = -Qdotnp-0.5*Bsq*(1.0+vsq)+0.5*t16-Wp+p_tmp;
  FTYPE t18norm = fabs(Qdotnp)+fabs(0.5*Bsq*(1.0+vsq))+fabs(0.5*t16)+fabs(Wp)+fabs(p_tmp);
  //t18 = -Qdotn-0.5*Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t21 = 1/t3;
  t23 = 1/W;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(Bsq-2.0*dPdvsq)*(QdotBsq+vsq*Wsq*W)*t9*t23;
  t36 = 1/t35;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;

#if(0)
  dualfprintf(fail_file,"Qdotn=%21.15g Qdotnp=%21.15g Wp=%21.15g p_tmp=%21.15g\n",Qdotn,Qdotnp,Wp,p_tmp);
  dualfprintf(fail_file,"t2=%21.15g t3=%21.15g t4=%21.15g t9=%21.15g t11=%21.15g t16=%21.15g t18=%21.15g t21=%21.15g t23=%21.15g t24=%21.15g t25=%21.15g t35=%21.15g t36=%21.15g dx[0]=%21.15g\n",t2,t3,t4,t9,t11,t16,t18,t21,t23,t24,t25,t35,t36,dx[0]);
#endif

  t40 = (vsq+t24)*t3;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  detJ = t3*t35;
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;
  norm[0] = t11norm;
  norm[1] = t18norm;

  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"Qdotn=%21.15g QdotBsq=%21.15g Qtsq=%21.15g Bsq=%21.15g vsq=%21.15g W=%21.15g p_tmp=%21.15g -Qdotnp-Wp=%21.15g\n",Qdotn,QdotBsq,Qtsq,Bsq,vsq,W,p_tmp,-Qdotnp-Wp);
    
    dualfprintf(fail_file,"t2=%21.15g\nt3=%21.15g\nt4=%21.15g\nt9=%21.15g\nt11=%21.15g\nt16=%21.15g\nt18=%21.15g\nt21=%21.15g\nt23=%21.15g\nt24=%21.15g\nt25=%21.15g\nt35=%21.15g\nt36=%21.15g\ndx[0]=%21.15g\nt40=%21.15g\ndx[1]=%21.15g\ndetJ=%21.15g\njac[0][0]=%21.15g\njac[0][1]=%21.15g\njac[1][0]=%21.15g\njac[1][1]=%21.15g\nresid[0]=%21.15g\nresid[1]=%21.15g\n",t2,t3,t4,t9,t11,t16,t18,t21,t23,t24,t25,t35,t36,dx[0],t40,dx[1],detJ,jac[0][0],jac[0][1],jac[1][0],jac[1][1],resid[0],resid[1]);
    
    dualfprintf(fail_file,"*df=%21.15g *f=%21.15g Wp=%21.15g W=%21.15g\n",*df,*f,Wp,W);
  }
#endif


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"func_vsq2(): x[0] = %21.15g , x[1] = %21.15g , dx[0] = %21.15g , dx[1] = %21.15g , f = %21.15g \n", 
                x[0], x[1],dx[0], dx[1], *f);
    
  }
#endif

}


/*********************************************************
**********************************************************/
static FTYPE res_sq_vsq(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  
  FTYPE  W, Wp, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE  resid[2],norm[2];
  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq(whicheos,EOSextra, Wp, D, vsq );
  //  dPdW   = dpdW_calc_vsq(whicheos,EOSextra, Wp, D, vsq );
  //  dPdvsq = dpdvsq_calc_Wp(whicheos,EOSextra, Wp, D, vsq );

  /* These expressions were calculated using Mathematica */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  resid[0]  = Qtsq - vsq*(Bsq + W)*(Bsq + W) + (QdotBsq*(Bsq + 2*W))/Wsq;
  resid[1]  = -Qdotnp - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - Wp + p_tmp;

  norm[0]  = fabs(Qtsq) + fabs(vsq*(Bsq + W)*(Bsq + W)) + fabs((QdotBsq*(Bsq + 2*W))/Wsq);
  norm[1]  = fabs(Qdotnp) + fabs((Bsq*(1 + vsq))/2.) + fabs(QdotBsq/(2.*Wsq)) + fabs(Wp) + fabs(p_tmp);


  tmp3 = 0.5 * ( resid[0]*resid[0] + resid[1]*resid[1] ) ;
  
#if(!OPTIMIZED)
  if( ltrace  ) {
    dualfprintf(fail_file,"res_sq_vsq(): W,vsq,resid0,resid1,f =   %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  \n",
                W,vsq,resid[0],resid[1],tmp3);
    dualfprintf(fail_file,"res_sq_vsq(): Qtsq, Bsq, QdotBsq, Qdotn =   %21.15g  %21.15g  %21.15g  %21.15g  \n",
                Qtsq, Bsq, QdotBsq, Qdotn);
  }
#endif
    
  return( tmp3 );
}


/*********************************************************
**********************************************************/
static FTYPE res_sq_vsq2(FTYPE x[], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{

  FTYPE  W, Wp, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE  resid[2],norm[2];
  FTYPE t3, t4 ;
  FTYPE t9 ;
  FTYPE t11;
  FTYPE t16;
  FTYPE t18;

  const int ltrace = 0;


  Wp = x[0];
  W = Wp+D;
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_Wp_vsq(whicheos,EOSextra, Wp, D, vsq );

  /* These expressions were calculated using Mathematica and made into efficient code using Maple*/
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  t3 = Bsq+W;
  t4 = t3*t3;//
  t9 = 1/Wsq;  //
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9; //
  FTYPE t11norm = fabs(Qtsq)+fabs(vsq*t4)+fabs(QdotBsq*(Bsq+2.0*W)*t9); //
  t16 = QdotBsq*t9;//
  t18 = -Qdotnp-0.5*Bsq*(1.0+vsq)+0.5*t16-Wp+p_tmp;//
  FTYPE t18norm = fabs(Qdotnp)+fabs(0.5*Bsq*(1.0+vsq))+fabs(0.5*t16)+fabs(Wp)+fabs(p_tmp);//
  resid[0] = t11;
  resid[1] = t18;
  norm[0] = t11norm;
  norm[1] = t18norm;

  tmp3 = 0.5 * ( resid[0]*resid[0] + resid[1]*resid[1] ) ;
  
#if(!OPTIMIZED)
  if( ltrace  ) {
    dualfprintf(fail_file,"res_sq_vsq(): W,vsq,resid0,resid1,f =   %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  \n",
                W,vsq,resid[0],resid[1],tmp3);
    dualfprintf(fail_file,"res_sq_vsq(): Qtsq, Bsq, QdotBsq, Qdotn =   %21.15g  %21.15g  %21.15g  %21.15g  \n",
                Qtsq, Bsq, QdotBsq, Qdotn);
  }
#endif

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"res: W=%21.15g Wsq=%21.15g p_tmp=%21.15g t3=%21.15g t4=%21.15g t9=%21.15g t11=%21.15g t16=%21.15g t18=%21.15g resid[0]=%21.15g resid[1]=%21.15g tmp3=%21.15g\n",W,Wsq,p_tmp,t3,t4,t9,t11,t16,t18,resid[0],resid[1],tmp3);
  }
#endif

    
  return( tmp3 );
}










/**************************************************** 


GENERAL PROCEDURES BELOW (except validate_x() and must choose function and residual)


*****************************************************/


static int find_root_1D_gen(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_1d[1];
  int retval;


  x_1d[0] = x0;

  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_1d, 1, USE_LINE_SEARCH, func_1d_orig, res_sq_1d_orig, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}



static int find_root_1D_gen_scn(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_1d[1];
  int retval;


  x_1d[0] = x0;

  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_1d, 1, USE_LINE_SEARCH, func_1d_orig_scn, res_sq_1d_orig_scn, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}



static int find_root_1D_gen_Eprime(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_1d[1];
  int retval;


  x_1d[0] = x0;

  // Note u sing func_Eprime_opt() and res_sq_Eprime_opt() : optimized versions.
  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_1d, 1, USE_LINE_SEARCH, func_Eprime_opt, res_sq_Eprime_opt, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}


static int find_root_1D_gen_Psq(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_1d[1];
  int retval;


  x_1d[0] = x0;

  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_1d, 1, USE_LINE_SEARCH, func_Psq, res_sq_Psq, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}




/// below not done yet!
static int find_root_3D_gen_Palpha(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_3d[3];
  int retval;


  x_3d[0] = x0;
  x_3d[1] = x0;
  x_3d[2] = x0;

  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_3d, 3, USE_LINE_SEARCH, func_Psq, res_sq_Psq, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_3d[0] = %21.15g x_3d[1] = %21.15g x_3d[2] = %21.15g \n", x_3d[0], x_3d[1], x_3d[2] );  
#endif
  }

  *xnew = x_3d[0];
  *xnew = x_3d[1];
  *xnew = x_3d[2];

  return(retval);
}

static int find_root_1D_gen_Sc(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE x_1d[1];
  int retval;


  x_1d[0] = x0;

  // Note, using optimized versions of func_Sc() and res_sq_Sc()
  if( (retval=general_newton_raphson(showmessages, lpflag, eomtype, x_1d, 1, USE_LINE_SEARCH, func_Sc_opt, res_sq_Sc_opt, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }

  *xnew = x_1d[0];

  return(retval);
}




/********************************************************************

   find_root_2D_gen(): 
  
         -- performs a 2D Newton-Raphson method for solving the \tilde{Q}^2
             and Q.n equations for two variables in general;
              (usually W and something else like utsq, vsq, or gamma)

         residual vector       = { Qtsq eq. , Q.n eq. }
         indep. var. vector, x = { W, vsq/utsq/gamma/? }

       vartype = 2 : W, gamma
       vartype = 3 : W, vsq
       vartype = 4 : W, utsq
    
     -- Uses  general_newton_raphson();

*********************************************************************/


static int find_root_2D_gen(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x0, FTYPE *xnew, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  
  FTYPE x[2], x_orig[2];
  int ntries = 0;
  int retval, n = 2;
  int it;

  const int ltrace = 0;

  /* Set presets: */
  x_orig[0] = x[0] = fabs(x0) ;
  x1_of_x0( x0, &x[1], wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;
  x_orig[1] = x[1];

#if(CRAZYDEBUG&&DEBUGINDEX)
  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    dualfprintf(fail_file,"x[0]=%21.15g W=%21.15g x[1]=%21.15g\n",x[0],x[0]+D,x[1]);
  }
#endif

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file, "find_root_2D_gen():  x[0] = %21.15g , x[1] = %21.15g \n", x[0], x[1]);
    
  }
#endif

#if(0)
  for(it=0;it<n;it++) dualfprintf(fail_file,"x[%d]=%21.15g xnew[%d]=%21.15g\n",it,x[it],it,xnew[it]);
#endif

  retval = general_newton_raphson(showmessages, lpflag, eomtype, x, n, USE_LINE_SEARCH,  func_vsq2, res_sq_vsq2, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ;
  ntries++;

  while( (retval==1) && (ntries <= MAX_NEWT_RETRIES ) ) {
    //       x[0] = x_orig[0] * (1. + ( (1.*rand())/(1.*RAND_MAX) - 0.5 ) );  
    x[0] = x_orig[0];
    for( it=0; it<ntries; it++)  x[0] *= 10.0;
    x1_of_x0( x[0],&x[1], wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;
    //       retval = general_newton_raphson(showmessages, lpflag, eomtype, x, n, USE_LINE_SEARCH,  func_vsq, res_sq_vsq, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ;
    retval = general_newton_raphson(showmessages, lpflag, eomtype, x, n, USE_LINE_SEARCH,  func_vsq2, res_sq_vsq2, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats) ;

#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file, "find_root_2D_gen():  ntries, x[0,1] = %4i  %21.15g   %21.15g  \n", ntries, x[0], x[1]);
      
    }
#endif
    ntries++;
  }

#if( MAX_NEWT_RETRIES > 0 )
  if( (ntries > MAX_NEWT_RETRIES) &&  (retval==1)  ) {
    if( debugfail >= 2 ) { 
      dualfprintf(fail_file, "find_root_2D_gen():  Bad exit value from general_newton_raphson() !! \n");
      dualfprintf(fail_file, "find_root_2D_gen():  ntries = %d , x[0] = %21.15g , x[1] = %21.15g  \n", ntries, x[0], x[1]); 
    }
  }
#endif

  *xnew = x[0];
  return(retval);

}




/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson(int showmessages, PFTYPE *lpflag, int eomtype, FTYPE x[], int n, int do_line_search,
                                  void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), 
                                  FTYPE (*res_func) (FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra, struct of_newtonstats *newtonstats)
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], x_older[NEWT_DIM], x_olderer[NEWT_DIM], resid[NEWT_DIM], norm[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, errx_old, errx_oldest, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE randtmp, tmp;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE resid_norm, resid_check, grad_check;

  FTYPE res_func_val, res_func_old, res_func_new;
  FTYPE dn[NEWT_DIM], del_f[NEWT_DIM];
  
  FTYPE trueerror[NEWT_DIM],trueerrortotal;

  // Jon's modification
  FTYPE DAMPFACTOR[NEWT_DIM];

  int   keep_iterating, i_increase, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 1;
  int it;
  int foundnan;
  int diddamp=0;
  int didcycle=0;
  void (*ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);
  int newtoniter;


  //////////
  //
  // use inputs from newtonstats that sets up iterations
  // see global.structs.h: of_newtonstats{}
  //
  /////////
  FTYPE NEWT_TOL_VAR;
  if(newtonstats->tryconv>0.0 && newtonstats->tryconv>NEWT_TOL) NEWT_TOL_VAR=newtonstats->tryconv;
  else NEWT_TOL_VAR=NEWT_TOL;

  FTYPE NEWT_TOL_ULTRAREL_VAR;
  if(newtonstats->tryconvultrarel>0.0 && newtonstats->tryconvultrarel>NEWT_TOL_ULTRAREL) NEWT_TOL_ULTRAREL_VAR=newtonstats->tryconvultrarel;
  else NEWT_TOL_ULTRAREL_VAR=NEWT_TOL_ULTRAREL;

  FTYPE MIN_NEWT_TOL_VAR;
  if(newtonstats->mintryconv>0.0 && newtonstats->mintryconv>MIN_NEWT_TOL) MIN_NEWT_TOL_VAR=newtonstats->mintryconv;
  else MIN_NEWT_TOL_VAR=MIN_NEWT_TOL;

  FTYPE MAX_NEWT_ITER_VAR;
  if(newtonstats->maxiter>=0 && newtonstats->maxiter<MAX_NEWT_ITER) MAX_NEWT_ITER_VAR=newtonstats->maxiter;
  else MAX_NEWT_ITER_VAR=MAX_NEWT_ITER;

  int EXTRA_NEWT_ITER_VAR;
  if(newtonstats->extra_newt_iter>=0 && newtonstats->extra_newt_iter<EXTRA_NEWT_ITER) EXTRA_NEWT_ITER_VAR=newtonstats->extra_newt_iter;
  else EXTRA_NEWT_ITER_VAR=EXTRA_NEWT_ITER;

  int EXTRA_NEWT_ITER_ULTRAREL_VAR;
  if(newtonstats->extra_newt_iter_ultrarel>=0 && newtonstats->extra_newt_iter_ultrarel<EXTRA_NEWT_ITER_ULTRAREL) EXTRA_NEWT_ITER_ULTRAREL_VAR=newtonstats->extra_newt_iter_ultrarel;
  else EXTRA_NEWT_ITER_ULTRAREL_VAR=EXTRA_NEWT_ITER_ULTRAREL;


  // setup error list
  FTYPE xlist[NEWT_DIM][MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER],errxlist[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER];
  FTYPE dxlist[NEWT_DIM][MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER],trueerrorlist[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER];
  int listi;
  for(it=0;it<n;it++){
    for(listi=0;listi<MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER;listi++){
      xlist[it][listi]=errxlist[listi]=BIG;
      dxlist[it][listi]=trueerrorlist[listi]=BIG;
    }
  }



  // pick version of validate_x depending upon scheme
  pick_validate_x(eomtype, &ptr_validate_x);


  retval = 0;


  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_olderer[id]=x_older[id]=x_old[id] = x_orig[id] = x[id] ; // save initial guess as orig and old
  vsq_old = vsq = W = W_old = 0.;

  for( id = 0; id < n ; id++) DAMPFACTOR[id]=1.0; // Jon's mod: start out fast
  for( id = 0; id < n ; id++) dx_old[id]=VERYBIG;// Assume bad error at first (notice not normalized to x[id] since x[id] might be 0




  //////////////////////////////////////////////////////////
  //
  /* Start the Newton-Raphson iterations : */
  //
  /////////////////////////////////////////////////////////
  n_iter = 0;
  keep_iterating = 1;
  newtoniter=0;
  while( keep_iterating ) { 
    (newtonstats->nstroke)++;
    (newtonstats->lntries)++;

    if((int)(newtonstats->lntries) > ITERDAMPSTART && diddamp==0 && didcycle==0){
      for( id = 0; id < n ; id++) DAMPFACTOR[id]=0.5; // Jon's mod: try to catch bad closed loop cycles that can occur, e.g., when cubic type behavior
    }



#if(CRAZYDEBUG&&DEBUGINDEX)
    // DEBUG:
    if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
      for(it=0;it<n;it++) dualfprintf(fail_file,"before funcd: x[%d]=%21.15g dx[%d]=%2.15g\n",it,x[it],it,dx[it]);
    }
#endif      


    /* returns with new dx, f, df */
    (*funcd) (x, dx, resid, norm, jac, &f, &df, n, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

    // get true error
    for(it=0;it<n;it++) trueerror[it]=fabs(resid[it]/(SMALL+fabs(norm[it])));
    trueerrortotal=0.0;
    for(it=0;it<n;it++) trueerrortotal = MAX(trueerrortotal,(1.0/(FTYPE)(n))*(fabs(resid[it])/(SMALL+fabs(norm[it]))));

    // store in list
    errxlist[n_iter] = errx;
    trueerrorlist[n_iter] = trueerrortotal;
    for(it=0;it<n;it++){
      dxlist[it][n_iter] = fabs(dx[it]);
      xlist[it][n_iter] = x[it];
    }

    
    // below now used currently
    id=0; if(resid[id] <=-VERYBIG || jac[0][id]==0.0){
      // get error and dx
      newtoniter++;
    }



    // DEBUG:
#if(0)
    //    if(nstep>=26070){
    //      if(1||myid==5 && nstep==1 && steppart==0 && ifileglobal==19 && jfileglobal==15){
    for(it=0;it<n;it++) dualfprintf(fail_file,"lntries=%d after funcd: x[%d]=%26.20g dx[%d]=%26.20g f=%26.20g df=%26.20g errx=%26.20g trueerror=%26.20g (%26.20g %26.20g) diddamp=%d dampfactor=%26.20g didcycle=%d : %g %g %g %g %g %g %g %g %g\n",(int)(newtonstats->lntries),it,x[it],it,dx[it],f,df,errx,resid[it]/norm[it],NEWT_TOL_VAR,NEWT_TOL_ULTRAREL_VAR,diddamp,DAMPFACTOR[it],didcycle,wglobal[2],Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc);
    //      }
    //    }
#endif

#if(CRAZYDEBUG&&DEBUGINDEX)
    // DEBUG:
    if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
      for(it=0;it<n;it++) dualfprintf(fail_file,"after funcd: x[%d]=%21.15g dx[%d]=%2.15g\n",it,x[it],it,dx[it]);
    }
#endif      


#if(!OPTIMIZED)
    /*  Check for bad untrapped divergences : */
    if( (finite(f)==0) ||  (finite(df)==0) ) {
      if( debugfail >= 2 ) { 
        dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
        dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
      }
      return(1);
    }
#endif


#if(!OPTIMIZED)
    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      // for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }
#endif

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_olderer[id]=x_older[id];  //x_olderer contains 2 steps ago
      x_older[id]=x_old[id];  //x_older contains last step's W used by above funcd() to get resid and dx
      x_old[id] = x[id] ; // x_old contains just used W that was used to compute funcd() to get present resid and dx
    }








    ///////////////////////////////////////////////
    //
    ////////////// START LINE SEARCH STUFF

    /* Make the newton step: */
    if( do_line_search == 1 ) { 

      /* Compare the residual to its initial value */ 
      if( n_iter == 0 ) { 
        resid_norm = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_norm += fabs(resid[id]);
        }
        resid_norm /= 1.0*n ;
        if( resid_norm == 0.0 ) resid_norm = 1.0;
      }

      for( id = 0; id < n ; id++) {
        tmp = 0.;
        for( jd = 0; jd < n ; jd++) {
          tmp += jac[jd][id] * resid[jd];
        }
        del_f[id] = tmp;
      }
      for( id = 0; id < n ; id++) {
        dn[id] = dx[id];
      }

      my_lnsrch(eomtype, n, x_old-1, f_old, del_f-1, dn-1, x-1, &f, TOL_LINE_STEP, SCALEMAX, &retval, res_func, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

      /* dx is needed for errx calculation below: */
      for( id = 0; id < n ; id++) {
        dx[id] = x[id] - x_old[id];
      }

#if(!OPTIMIZED)
      if( ltrace ) { 
        res_func_val = res_func(x, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc );
        res_func_old = res_func(x_old, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc );
        dualfprintf(fail_file,"gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
                    f_old, f, res_func_old, res_func_val );
        dualfprintf(fail_file,"gnr(): x_old = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",x_old[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): x     = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",x[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): dn    = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",dn[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): del_f = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",del_f[id]);
        }
        dualfprintf(fail_file,"\n ");
      }
#endif

      /* Check to see if line search problem is because the residual vector is already small enough */
      if( retval == 1 ) {
        resid_check = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_check += fabs(resid[id]);
        }
        resid_check /= 1.0*n;
 
        if( resid_check <= resid_norm * NEWT_FUNC_TOL ) {
          retval = 0;
        }
        if( ltrace && retval ) { 
          dualfprintf(fail_file,"general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
        }   
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH == do_line_search) ) { 
        if(debugfail>=2) dualfprintf(fail_file,"bad first step\n");
#if(!OPTIMIZED)   
        if( ltrace ) { 
          dualfprintf(fail_file,"gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
          dualfprintf(fail_file,"gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
        }
#endif
        retval = general_newton_raphson(showmessages, lpflag, eomtype, x_orig, n, ((do_line_search+1)%2), funcd, res_func, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
        for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
        return( retval );
      }

      /* Check to see if it is trapped in a local minimum, i.e. gradient is too small */ 
      if( retval == 1 ) { 
        grad_check = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_check = (x[id] == 0.) ? 1.0 : fabs(x[id]) ;
          grad_check  +=  del_f[id] * resid_check ;
        }
        resid_check = (f == 0.) ? 1.0 : fabs(f) ;
        grad_check /= resid_check;
 
        /* Then we've most likely found a solution: */
        if( grad_check > GRADMIN ) { 
          retval = -1;
        }
        else if( ltrace ) { 
          dualfprintf(fail_file,"general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval, grad_check);
        }
      }
    }
    ///////////////////////////////////////////////
    //
    ////////////// END LINE SEARCH STUFF


    else {
      /* don't use line search : */
      diddamp=0;

#if(JONGENNEWTSTUFF)
      if(n_iter==0){
        for(id=0;id<n;id++) {
          // start history
          x_olderer[id]=x[id];
          x_older[id]=x[id];
          x_old[id]=x[id];
          dx_old[id]=dx[id];
          //   xold[id]=x[id]; // original W
          //   dxold[id]=dx[id]; // original dW
          //   residold[id]=resid[id]; // original residual
        }
      }
      else{ // only meaningful to consider goodness of W if compute residual from present W, not previous W.
        for(id=0;id<n;id++){
          //  if(sign(residold[id])!=sign(resid[id])){ //  then passed over root -- bad in very relativistic case
          if(resid[id] <=-VERYBIG || jac[0][id]==0.0 ){ // then went too far such that v^2<0 or \tilde{u}^2<0 or Ss=NaN
            //   if(1 || fabs(xold[0])>wglobal[2] || fabs(x_old[0])>wglobal[2]){ // always do now under more careful choice of when returns bad residual indicating problem
            if(1 || fabs(x_old[0])>wglobal[2]){ // always do now under more careful choice of when returns bad residual indicating problem
              // then ultrarelativistic and need to be more careful
              DAMPFACTOR[id]=0.5*DAMPFACTOR[id];
              x[id]=x_older[id]; // since current W is what gave bad residual and already saved it into x_old[0]
              dx[id]=dx_old[id]; // new dx is probably messed up (didn't yet save dx[0], so use dx_old[0]
#if(PRODUCTION==0)
              if(showmessages&&debugfail>=3) dualfprintf(fail_file,"resetW: nstep=%ld steppart=%d :: lntries=%d :: id=%d :: x=%21.15g dx=%21.15g : x/D=%21.15g DAMPFACTOR=%21.15g errx=%21.15g eomtype=%d\n",nstep,steppart,(int)(newtonstats->lntries),id,x[id],dx[id],x[id]/D,DAMPFACTOR[id],newtonstats->lerrx,eomtype);
#endif
              // SUPERGODMARK: Noticed that if Mathematica solution can be found but gives utsq<0 then this damping leads to nearly circular loop till max iterations.
              diddamp=1;
            }
          }
        } // end loop over id's
        if(diddamp){
          // Ensure this newly chosen "old" value gives reasonable result
          (*funcd) (x, dx, resid, norm, jac, &f, &df, n, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);  /* returns with new dx, f, df */
          for(id=0;id<n;id++){
            if(resid[id] <=-VERYBIG || jac[0][id]==0.0 ){ // then went too far such that v^2<0 or \tilde{u}^2<0 or Ss=nan or inf
#if(PRODUCTION==0)
              if(debugfail>=2) dualfprintf(fail_file,"Old value failed too: resid[%d]=%26.20g jac[0][%d]=%26.20g\n",id,resid[id],id,jac[0][id]);
#endif
              return(1); // then old fails to work too
            }
            // reset old values so will use this good value if immediately hit bad solution
            // This restarts history
            x_olderer[id]=x[id];
            x_older[id]=x[id];
            x_old[id]=x[id];
            dx_old[id]=dx[id];
          }
        }
#if(1)
        else{
          // can only use x_olderer if 3 steps, but add 1 extra.
          // Only do if not just doing extra steps after already error is low (so avoid errx=0 case where will cycle).
          if(n_iter>4 && doing_extra==0){
            // see if cycling Newton steps with no change in errx
#if(0)
            // DEBUG:
            if(nstep==19 && steppart==1){// && i==0 && j==0 && k==26){
              dualfprintf(fail_file,"errx=%26.20g errx_old=%26.20g errx_oldest=%26.20g dx=%26.20g dx_old=%26.20g x=%26.20g x_old=%26.20g x_older=%26.20g  x_olderer=%26.20g\n",errx,errx_old,errx_oldest,dx[0],dx_old[0],x[0],x_old[0],x_older[0],x_olderer[0]);
              dualfprintf(fail_file,"differrx=%26.20g is=%d diffx=%26.20g is=%d\n",errx_old-errx_oldest,errx_old==errx_oldest,x[0]-x_olderer[0],x[0]==x_olderer[0]);
            }
#endif

            if(newt_cyclecheck(n, trueerrortotal, errx, errx_old, errx_oldest, dx, dx_old, x, x_old, x_older, x_olderer, n_iter, xlist,dxlist, errxlist, trueerrorlist)){
              //            if( newt_cyclecheck(n,errx,errx_old,errx_oldest,dx,dx_old,x,x_old,x_older,x_olderer) ){
              if(debugfail>=2) dualfprintf(fail_file,"Caught in cycle, trying to damp: n_iter=%d.\n",n_iter);
              // if hits another cycle, damp more
              for(id=0;id<n;id++) DAMPFACTOR[id]=0.1*DAMPFACTOR[id]; // severe damp to try to get good step.
              diddamp=1; // tells start of loop to not reset DAMPFACTOR
              didcycle++; // tells not to reset damp because having trouble and should stay damped
            }
          }// end if avoiding cycle
        }
#endif
      }// end else if not first step


#endif
      

      /////////////////////////
      //
      // Step forward with x+dx
      //
      /////////////////////////
      if(eomtype==EOMENTROPYGRMHD){
        // Function to find root of is like log(Wp) , so near root, Newton can jump to negative Wp.  But for entropy method, this cannot be allowed without using very small DAMPFACTOR.  So instead, somewhat bracket the root between original Wp and 0.0
        for( id = 0; id < n ; id++) {
          x[id] = MAX((x[id]+DAMPFACTOR[id]*dx[id]),(x[id]-0.0)/100.0);

          // DEBUG:
          //   dualfprintf(fail_file,"id=%d newx=%21.15g using dx=%21.15g\n",id,x[id],dx[id]);
        }
      }
      else{
        for( id = 0; id < n ; id++) {
          x[id] +=DAMPFACTOR[id]*dx[id]  ;
        }
      }

      /////
      //
      // undo damping if no longer damping
      // Otherwise too slow eventually, and otherwise error based upon dx thinks error is much smaller than really is.
      //
      /////
      for( id = 0; id < n ; id++) {
        if(diddamp==0) DAMPFACTOR[id]=MIN(1.0,2.0*DAMPFACTOR[id]);
      }


#if(CRAZYDEBUG&&DEBUGINDEX)
      if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
        dualfprintf(fail_file,"no line search: x[0]=%21.15g dx[0]=%21.15g\n",x[0],dx[0]);
      }
#endif

    } /* End of "to do line search or not to do line search..." */







      /////////////////////////////////////////////////////////
      //
      //  Calculate the convergence criterion */
      //
      /////////////////////////////////////////////////////////

      /* For the new criterion, always look at error in "W" : */
      // METHOD specific: 

#if( NEWCONVERGE == 1 )



#if(JONGENNEWTSTUFF==0)
    // SUPERGODMARK: is wglobal always the right normalization for the error?

    //    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);
    //    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/wglobal[0]);
    
    errx  = (wglobal[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/MAX(fabs(wglobal[0]),fabs(x[0])));

#elif(JONGENNEWTSTUFF==1)
    // makes no sense to compute error from present dx for present x.  The error in x[0] is unknown.  errx here will be error in xold
    // actually error of previous x[0]
    // in this way, final funcd() call is actually an error check only, not an update to x[0]!  So final solution returned should be x_old
    errx  = (wglobal[0]==0.) ?  fabs(dx_old[0]) : fabs(dx_old[0]/MAX(fabs(wglobal[0]),fabs(x[0])));
    //    errx  = (wglobal[0]==0.) ?  fabs(dx_old[0]) : fabs(dx_old[0]/fabs(wglobal[0]));
#endif




#if(CRAZYDEBUG&&DEBUGINDEX)
    if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"x[0]=%21.15g dx[0]=%21.15g\n",x[0],dx[0]);
    }
#endif


    /* For the old criterion, look at errors in each indep. variable(s) (except for 5D) : */
#else
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;
#endif


    //dualfprintf(fail_file,"x[0]=%21.15g dx[0]=%21.15g df/f=%21.15g errx=%21.15g out of %21.15g or %21.15g\n",x[0],dx[0],df/f,errx,NEWT_TOL_VAR,newtonstats->tryconv);


    // save errx for information/debug purposes
    newtonstats->lerrx=errx;


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific:

    
    validate_x(ptr_validate_x, x, x_old, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra) ;


    /****************************************/
    /* Check to see if we're in a infinite loop with error function: */
    /****************************************/
#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL_VAR) )  errx = -errx;
#endif 

    /****************************************/
    /* If there's a problem with line search, then stop iterating: */
    /****************************************/
    if( (retval == 1) || (retval == -1) ) errx = -errx;


#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file," general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
                  n_iter,f_old,f,errx_old,errx );
      dualfprintf(fail_file,"gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",x_old[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): x     = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",x[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): dx     = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",dx[id]);
      }
      dualfprintf(fail_file,"\n ");
      
    }
#endif

    /****************************************/
    /* Prepare for the next iteration, set the "old" variables: */
    /****************************************/
    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    /****************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations before stopping */
    /****************************************/
   

    
    // see if reached tolerance
    if(newt_errorcheck(n, NEWT_TOL_VAR, NEWT_TOL_ULTRAREL_VAR, errx, trueerrortotal, x[0],  dx[0], x, dx, wglobal) && (doing_extra == 0) && (newt_extracheck(n, EXTRA_NEWT_ITER_VAR, EXTRA_NEWT_ITER_ULTRAREL_VAR, errx,x[0],dx[0],x, dx, wglobal) > 0) ) {
      doing_extra = 1;
    }

    // see if repeating Newton steps with no change in errx
    if(doing_extra==0 && newt_repeatcheck(n,errx,errx_old,errx_oldest,dx,dx_old,x,x_old,x_older,x_olderer) ){
      if(debugfail>=2) dualfprintf(fail_file,"Repeat found\n");
      doing_extra = 1;
    }


    if( doing_extra == 1 ) i_extra++ ;

    if( (newt_errorcheck(n, NEWT_TOL_VAR, NEWT_TOL_ULTRAREL_VAR, errx, trueerrortotal, x[0],  dx[0],x,dx, wglobal)&&(doing_extra == 0)) || (i_extra > newt_extracheck(n, EXTRA_NEWT_ITER_VAR, EXTRA_NEWT_ITER_ULTRAREL_VAR, errx,x[0],dx[0],x,dx,wglobal)) || (n_iter >= (MAX_NEWT_ITER_VAR-1)) ) {
      keep_iterating = 0;
    }

    // if cycling twice, then jump out and assume repeating cycle with no hope of improvement.  Do so after CYCLESTOP steps
    if(n_iter>CYCLESTOP && didcycle>1){
      keep_iterating = 0;
    }



#if(0)
    // check if error is not dropping sufficiently
#define CHECKDECREASE0 5
#define CHECKDECREASEAVGNUM 3 // should be less than CHECKDECREASE0
#define CHECKDECFACTOR (0.5) // i.e. should drop by this factor compared to older average
    if(0&&newtoniter>CHECKDECREASE0){ // too aggressive, causes no solutions more often when could have gotten solution
      int ci;
      FTYPE avgerror;
      avgerror=0.0;
      for(ci=0;ci<CHECKDECREASEAVGNUM;ci++) if(errorlist[newtoniter-CHECKDECREASE0+ci]!=1.0) avgerror += errorlist[newtoniter-CHECKDECREASE0+ci]; // GODMARK: uses debug info
      avgerror/=(CHECKDECREASEAVGNUM);
      FTYPE currenterror;
      currenterror=errorlist[newtoniter];
      //
      // check both errors to ensure they are decreasing
      int cond=(currenterror>CHECKDECFACTOR*avgerror);
      if(cond){
        keep_iterating=0;
        if(debugfail>=2) dualfprintf(fail_file,"Error not decreasing properly: currenterror=%g avgerror=%g\n",currenterror,avgerror);
      }// end if current error too large compared to older average error
    }// end if large enough iterations so can check how error is trending, to avoid many iterations when error is not dropping enough to matter.
#endif





    n_iter++;

#if(CRAZYDEBUG&&DEBUGINDEX)
    //i=0 j=63 part=2 step=9 :: n_iter=6
    //  if(ifileglobal==0 && jfileglobal==31 && nstep==9 && steppart==2){
    //    if(ifileglobal==0 && jfileglobal==63 && nstep==9 && steppart==2){
    //    if(n_iter>20){
    //if(ifileglobal>90 && CRAZYNEWCHECK){
    if(CRAZYNEWCHECK){
      //    if(ifileglobal==0 && jfileglobal==63 && nstep==12 && steppart==0){
      //    if(ifileglobal==0 && jfileglobal==0 && nstep==0 && steppart==0){
      dualfprintf(fail_file,"i=%d j=%d part=%d step=%ld :: n_iter=%d :: errx=%21.15g minerr=%21.15g :: x[0]=%21.15g dx[0]=%21.15g wglobal0=%21.15g\n",ifileglobal,jfileglobal,steppart,nstep,n_iter,errx,MIN_NEWT_TOL_VAR,x[0],dx[0],wglobal[0]);
    }
    //  }
#endif

  }   // END of while(keep_iterating)










  //  dualfprintf(fail_file,"n_iter=%d\n",n_iter);



  /////////////////////////////////////////////////////////
  //
  //  Check for bad untrapped divergences
  //
  /////////////////////////////////////////////////////////

  foundnan=0; // assume no nan's
  for(id=0;id<n;id++){
    if(finite(x[id])==0)      foundnan=1;
  }

  if( (finite(f)==0) ||  (finite(df)==0) ){  // deprecated: || (real_dimen_newton>=1)&&(finite(x[0])==0) || (real_dimen_newton>=2)&&(finite(x[1])==0)) {
    foundnan=1;
  }


  if(foundnan){
#if(!OPTIMIZED || PRODUCTION==0)
    if(debugfail>=1){
      dualfprintf(fail_file,"naninf: f=%21.15g df=%21.15g\n",f,df);
      for(id=0;id<n;id++){
        dualfprintf(fail_file,"naninf: x[%d]=%21.15g\n",id,x[id]);
      }
    }
    if( debugfail >= 2 ) { 
      dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
    }
#endif
    *lpflag=UTOPRIMFAILNANRESULT;
    return(1);
  }





  /////////////////////////////////////////////////////////
  //
  //  Check error
  //
  /////////////////////////////////////////////////////////

  // assign final result
#if(JONGENNEWTSTUFF==1)
  for(id=0;id<n;id++) {
    x[id]=x_old[id]; // actually evaluated x
  }
#endif





  if(newt_errorcheck(n, NEWT_TOL_VAR, NEWT_TOL_ULTRAREL_VAR, errx, trueerrortotal, x[0],  dx[0], x,dx,wglobal) ){
#if(DOHISTOGRAM)
    bin_newt_data(MAX_NEWT_ITER_VAR, errx, n_iter, 2, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      
    }
    //    dualfprintf(fail_file,"gnr retval5 = %4i \n", 0); 
#endif


    return(0);
  }



  if( (fabs(errx) <= MIN_NEWT_TOL_VAR) && newt_errorcheck(n, NEWT_TOL_VAR, NEWT_TOL_ULTRAREL_VAR, errx, trueerrortotal, x[0],  dx[0], x,dx,wglobal) ){
#if(DOHISTOGRAM)
    bin_newt_data(MAX_NEWT_ITER_VAR, errx, n_iter, 1, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      
    }
    if(ltrace) {
      dualfprintf(fail_file,"general_newton_raphson(): found minimal solution \n");
      
    }
    //    dualfprintf(fail_file,"gnr retval4 = %4i \n", 0); 
#endif


#if(JONGENNEWTSTUFF==1)
    for(id=0;id<n;id++) {
      x[id]=x_old[id]; // actually evaluated x
    }
#endif

    return(0);
  }



  if( fabs(errx) > MIN_NEWT_TOL_VAR){
    if( (do_line_search != USE_LINE_SEARCH) || (USE_LINE_SEARCH < 0) ) { 
#if(DOHISTOGRAM)
      bin_newt_data(MAX_NEWT_ITER_VAR, errx, n_iter, 0, 0 );
#endif

      if(debugfail>=2){
        dualfprintf(fail_file,"fabs(errx)=%21.15g > MIN_NEWT_TOL_VAR=%21.15g n=%d n_iter=%d lntries=%d eomtype=%d\n",fabs(errx),MIN_NEWT_TOL_VAR,n,n_iter,newtonstats->lntries,eomtype);
        //        if(DEBUGINDEX) dualfprintf(fail_file,"i=%d j=%d part=%d step=%ld :: n_iter=%d :: errx=%21.15g minerr=%21.15g :: x[0]=%21.15g dx[0]=%21.15g wglobal0=%21.15g\n",ifileglobal,jfileglobal,steppart,nstep,n_iter,errx,MIN_NEWT_TOL_VAR,x[0],dx[0],wglobal[0]);
      }



#if(!OPTIMIZED)

      if(ltrace2) {
        dualfprintf(fail_file," totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      }
      if(ltrace) {
        dualfprintf(fail_file,"general_newton_raphson():  did not find solution \n");
        if( retval == -1 ) {
          dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x = ");
          for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x[id]);
          dualfprintf(fail_file,"\n");
          dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x_old = ");
          for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x_old[id]);
          dualfprintf(fail_file,"\n");
        }
 
      }
      //      dualfprintf(fail_file,"gnr retval2 = %4i \n", 1); 
#endif
      *lpflag=UTOPRIMFAILCONVRET+4;
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
      //      dualfprintf(fail_file,"gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n", do_line_search, retval, errx );
      //      
      retval2 = general_newton_raphson(showmessages, lpflag, eomtype, x_orig, n, ((do_line_search+1)%2), funcd, res_func, wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra,newtonstats);
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      //      dualfprintf(fail_file,"gnr retval3 = %4i \n", retval2); 

      if(debugfail>=2) dualfprintf(fail_file,"fabs(errx) > MIN_NEWT_TOL_VAR: n_iter=%d lntries=%d\n",n_iter,newtonstats->lntries);

      return( retval2 );
    }
  }




#if(!OPTIMIZED)
  dualfprintf(fail_file,"gnr retval6 = %4i \n", 0);
#endif
  return(0);

}













/// see if meet tolerance
static int newt_errorcheck(int n, FTYPE NEWT_TOL_VAR, FTYPE NEWT_TOL_ULTRAREL_VAR, FTYPE errx, FTYPE trueerrortotal, FTYPE x0, FTYPE dx0, FTYPE *x, FTYPE *dx, FTYPE *wglobal)
{

  int ultrarel=0;

  ultrarel=( (fabs(wglobal[1])>wglobal[2]) || (fabs(x0)>wglobal[2]) );

  // see if dx's are too small to matter, despite condition on errx.
  int dxsmallcount=0;
  int iter;
  for(iter=0;iter<n;iter++) dxsmallcount += (fabs(dx[iter])<NUMEPSILON*x[iter]);

  if(!ultrarel){
    return(fabs(errx) <= NEWT_TOL_VAR || fabs(trueerrortotal) <= NEWT_TOL_VAR || dxsmallcount==n);
  }
  else{
    return(fabs(errx) <= NEWT_TOL_ULTRAREL_VAR || fabs(trueerrortotal) <= NEWT_TOL_ULTRAREL_VAR || dxsmallcount==n);
  }
}



/// see if repeating Newton cycle with no change in errx, so indicates can't do anything (even damping) to reduce error
static int newt_repeatcheck(int n, FTYPE errx, FTYPE errx_old, FTYPE errx_oldest, FTYPE *dx, FTYPE *dx_old, FTYPE *x, FTYPE *x_old, FTYPE *x_older, FTYPE *x_olderer)
{
  int repeat;
  int i;

  repeat=0;
  for(i=0;i<n;i++){
  
    if( (x[i]==x_olderer[i] && fabs(dx[i])==fabs(dx_old[i])) && (errx_oldest==errx_old) ){
      repeat++;
    }
  }

  if(repeat==n) return(1);
  else return(0);


}

//  FTYPE xlist[NEWT_DIM][MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER],errxlist[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER];
//  FTYPE dxlist[NEWT_DIM][MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER],trueerrorlist[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER];


/// see if repeating Newton cycle with no change in errx, so indicates can't do anything (even damping) to reduce error
static int newt_cyclecheck(int n, FTYPE trueerror, FTYPE errx, FTYPE errx_old, FTYPE errx_oldest, FTYPE *dx, FTYPE *dx_old, FTYPE *x, FTYPE *x_old, FTYPE *x_older, FTYPE *x_olderer, int n_iter, FTYPE (*xlist)[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER], FTYPE (*dxlist)[MAX_NEWT_ITER+EXTRA_NEWT_ITER_ULTRAREL+EXTRA_NEWT_ITER], FTYPE *errxlist, FTYPE *trueerrorlist)
{
  int i;

  int cycle=0;
  for(i=0;i<n;i++){
  
    //    dualfprintf(fail_file,"%d %d %g %g %g %g\n",x[i]==x_olderer[i],errx_oldest==errx_old,x[i],x_olderer[i],errx_oldest,errx_old);
    if( (x[i]==x_olderer[i]) && (errx_oldest==errx_old && errx_old!=0.0) ){
      cycle++;
    }
  }


  int cycle2=0;
  int ii;
  for(ii=CYCLESTOP;ii<n_iter;ii++){
    if(errx_old!=0.0 && errx_old==errxlist[ii]) cycle2++;
    if(trueerror!=0.0 && trueerror==trueerrorlist[ii]) cycle2++;
    for(i=0;i<n;i++){
      if(x[i]==xlist[i][ii]) cycle2++;
    }
  }




  if(cycle>=n || cycle2>=NUMCYCLES*n) return(1);
  else return(0);


}

// see if done with extra iterations
static int newt_extracheck(int n, int EXTRA_NEWT_ITER_VAR, int EXTRA_NEWT_ITER_ULTRAREL_VAR, FTYPE errx, FTYPE x0, FTYPE dx0, FTYPE *x, FTYPE *dx, FTYPE *wglobal)
{

  int ultrarel=0;

  ultrarel=( (fabs(wglobal[1])>wglobal[2]) || (fabs(x0)>wglobal[2]) );

  // see if dx's are too small to matter, despite condition on errx.
  int dxsmallcount=0;
  int iter;
  for(iter=0;iter<n;iter++) dxsmallcount += (fabs(dx[iter])<NUMEPSILON*x[iter]);

  if(dxsmallcount) return(0);
  else{
    if(!ultrarel){
      return(EXTRA_NEWT_ITER_VAR);
    }
    else{
      return(EXTRA_NEWT_ITER_ULTRAREL_VAR);
    }
  }
}


#define ALF 1.0e-4

static void my_lnsrch(int eomtype, int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
                      FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE [], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra), FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  int i;
  FTYPE a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  int icount=0;
  int bad_step;
  FTYPE bad_step_factor = 2.0;
  const int ltrace = 0;
  void (*ptr_validate_x)(FTYPE x[NEWT_DIM], FTYPE x0[NEWT_DIM], FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra);



  // pick version of validate_x depending upon scheme
  pick_validate_x(eomtype, &ptr_validate_x);

 
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    //    temp=fabs(p[i])/max(fabs(xold[i]),1.0);
    temp= (xold[i] == 0.) ? fabs(p[i]) :  fabs(p[i]/xold[i]);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"my_lnsrch(): sum, slope, test, alamin =   %21.15g  %21.15g  %21.15g  %21.15g \n",sum,slope,test,alamin); 
  }
#endif

  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

    // METHOD specific:
    validate_x(ptr_validate_x, (x+1), (xold+1), wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);

    *f=(*func)(x+1,wglobal,Bsq,QdotB,QdotBsq,Qtsq,Qdotn,Qdotnp,D,Sc,whicheos,EOSextra);
      
    bad_step = 0;

    if( finite(*f)==0 ) { 
      bad_step = 1;
    }


    //      if( bad_step ) alam /= bad_step_factor;
    //      if (alam < alamin) bad_step = 0;

    if( bad_step ) { 
      *check = 2;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"my_lnsrch(): bad_step = 1,  f = %21.15g , alam = %21.15g \n", *f, alam); 
      for( i = 1; i <= n; i++ ) { 
        dualfprintf(fail_file,"my_lnsrch(): (xold[%d], x[%d], p[%d]) = %21.15g , %21.15g , %21.15g \n", i,i,i, xold[i],x[i],p[i]);
      }
#endif
      return;
    }
      
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
#if(!OPTIMIZED)
      if( ltrace ) { 
        dualfprintf(fail_file,"my_lnsrch(): alam < alamin: alam, alamin = %21.15g  %21.15g \n", alam,alamin); 
      }
#endif
      return;
    } 
    else if (*f <= fold+ALF*alam*slope) {
#if(!OPTIMIZED)
      if( ltrace ) { 
        dualfprintf(fail_file,"my_lnsrch(): good exit:  alam, alamin, f, fold = %21.15g  %21.15g %21.15g  %21.15g \n", alam,alamin, *f, fold); 
      }
#endif
      return;
    }
    else {
      if (alam == 1.0) {
        tmplam = -slope/(2.0*(*f-fold-slope));
#if(!OPTIMIZED)
        if( ltrace ) {
          dualfprintf(fail_file,"my_lnsrch(): setting tmplam!!    tmplam, alam =  %21.15g  %21.15g !!\n", tmplam, alam); 
        }
#endif
      }
      else {
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold2-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0) tmplam = -slope/(2.0*b);
        else {
          disc=b*b-3.0*a*slope;
          if (disc<0.0) {
#if(!OPTIMIZED)
            if( disc < -1.e-10 ) {
              if( ltrace ) dualfprintf(fail_file,"my_lnsrch(): Big Roundoff problem:  disc = %21.15g \n", disc);
            }
#endif       
            disc = 0.;
          }
          else tmplam=(-b+sqrt(disc))/(3.0*a);
        }
        if (tmplam>0.5*alam)
          tmplam=0.5*alam;
#if(!OPTIMIZED)
        if( ltrace ) {
          dualfprintf(fail_file,"my_lnsrch(): rhs1, rhs2, a, b, tmplam, alam =  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g !!\n",
                      rhs1, rhs2, a, b, tmplam, alam );
        }
#endif
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"my_lnsrch(): icount, alam, alam2, tmplam  =   %4i  %21.15g  %21.15g  %21.15g \n",
                  icount, alam, alam2, tmplam);
    }
#endif    
    icount++;
  }
}
#undef ALF


#define N_CONV_TYPES 3
#define N_NITER_BINS (MAX_NEWT_ITER + 3) // want to use MAX_NEWT_ITER_VAR, but not fixed, so constrain to be less.
#define NBINS 200

/**************************************************** 

  Used to gather statistics for the Newton solver during a disk run:

  -- bins are inclusive on their minimum bounds, and exlusive on their max.


      lerrx_[min,max] = the min./max. bounds of the errx historgram
      d_errx          = errx histogram bin width
      n_beyond_range  = number of solutions with errx max. beyond range speciified by lerrx_[min,max]
      n_conv_types    = number of types of ways newton procedure can exit;

      print_now       = 1 if you want to just print the histograms without counting a solution;

*****************************************************/

// OPENMPMARK: Not thread safe due to use of static variables, but not used except for debugging
static void bin_newt_data(FTYPE MAX_NEWT_ITER_VAR, FTYPE errx, int niters, int conv_type, int print_now  ) {

  
  /* General variables */
  int i, j;
  static int first_call   = 1;
  static long int n_bin_calls  = 0L;
  const  long int n_output_freq = 128*128*100;
  FTYPE localerrx;

  /* Variables for the errx histogram */ 
  static const FTYPE localerrx_min = -23.;
  static const FTYPE localerrx_max =   4.;
  static FTYPE d_errx;
  static long int n_errx[N_CONV_TYPES][NBINS];
  static FTYPE xbin[NBINS];
  static long int n_beyond_range = 0L;   /* Number of points that lie out of the bounds of our errx histogram */
  int ibin;

  /* Variables for the histogram of the number of newton iterations : */
  static long int n_niters[N_CONV_TYPES][N_NITER_BINS];
  

#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"bin_newt_data() called in parallel region\n");
    myexit(93586346);
  }
#endif


  /* Clear arrays, set constants : */
  if( first_call ) {
    d_errx    = ((localerrx_max - localerrx_min)/(1.*NBINS));
    
    for( i = 0; i < N_CONV_TYPES; i++ ) { 
      for( j = 0; j < NBINS; j++ ) { 
        n_errx[i][j] = 0;
      }
      for( j = 0; j < N_NITER_BINS; j++ ) { 
        n_niters[i][j] = 0;
      }
    }

    for( j = 0; j < NBINS; j++ ) { 
      xbin[j] = localerrx_min + (1.*j + 0.5)*d_errx;
    }
      
    first_call = 0;
  }


  if( print_now != 1 ) {

    /* Check validity of arguments : */
    errx = fabs(errx) ;
    localerrx = log10(errx + 1.0e-2*pow(10.,localerrx_min));

    if( (niters < 0) || (niters >= N_NITER_BINS) ) {
      dualfprintf(fail_file,"bin_newt_data(): bad value for niters = %d \n", niters );
      fflush(stdout);
      return;
    }
    
    /* Determine ibin */
    if( localerrx < localerrx_min ) {
      ibin = 0;
      n_beyond_range++ ;
    }
    else if( localerrx >= localerrx_max ) {
      ibin = NBINS - 1;
      n_beyond_range++ ;
    }
    else {
      ibin = (int) ( ( localerrx - localerrx_min ) / d_errx );
    }
    
    /* Tally this solution  */
    n_errx[ conv_type][ibin]++;
    n_niters[conv_type][niters]++;

  }


  /* Print out the histograms periodically or when asked to : */
  if( print_now  ||  ( (n_bin_calls % n_output_freq) == 0 ) ) {

    dualfprintf(log_file,"t = %21.15g ,  n_beyond_range = %ld , n_bin_calls = %ld \n", t, n_beyond_range, n_bin_calls);

    /* ERRX */
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");
    dualfprintf(log_file,"                         x");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < NBINS; i++ ) { 
      dualfprintf(log_file,"%21.15g ",xbin[i]);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        dualfprintf(log_file,"%13ld ", n_errx[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");


    /* NITER */

    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");
    dualfprintf(log_file,"        niter");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < N_NITER_BINS; i++ ) { 
      dualfprintf(log_file,"%13d ", i);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        dualfprintf(log_file,"%13ld ", n_niters[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");

  }
    
    
  if( print_now != 1 )  n_bin_calls++;


  return;

}
 
 
#undef N_CONV_TYPES 
#undef N_NITER_BINS 
#undef NBINS 

#undef CRAZYDEBUG

