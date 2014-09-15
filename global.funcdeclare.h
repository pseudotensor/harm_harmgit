

/*! \file global.funcdeclare.h
  \brief Function declarations (used globally) by entire code
 */

/// SPECIAL USE OF PTRDEF: (and in metric.c):
///  int interpX_gcov(FTYPE *X, struct of_compgeom (*compgeom)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*gcovgrid)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][SYMMATRIXNDIM], FTYPE (*gcovpertgrid)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM], FTYPE *gcov, FTYPE *gcovpert);
int interpX_gcov(FTYPE *X, struct of_compgeom PTRDEFMETMACP1A0(compgeom,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3), FTYPE PTRDEFMETMACP1A2(gcovgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM), FTYPE PTRDEFMETMACP1A1(gcovpertgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM), FTYPE *gcov, FTYPE *gcovpert);








extern int main(int argc, char *argv[]);
extern int init(int *argc, char **argv[]);

extern void parainitchecks(void);
extern void myargs(int argc, char *argv[]);


extern void pre_interpolate_and_advance(FTYPE (*pb)[NSTORE2][NSTORE3][NPR]);


extern int set_dt(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE *dt);



#include "copyandinit_functions.funcdeclare.h"


#include "flux.funcdeclare.h"


#include "metric_selfgravity_or_evolvemetric.funcdeclare.h"



extern void get_inversion_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke);
extern void get_stag_startendindices(int *loop, int dir, int *is,int *ie,int *js,int *je,int *ks,int *ke);
extern void get_flux_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke);


extern int avg2cen_interp(int *locpl, int *whichpltoavg,  int *ifnotavgthencopy, int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);



extern void set_defaults_performance_checks_prepreinit(void);
extern void set_defaults_performance_checks_preinit(void);
extern void set_file_versionnumbers(void);

extern int timecheck(int whichlocation, SFTYPE comptstart);
extern int gocheck(int whichlocation);
extern int output_steptimedt_info(SFTYPE comptstart);

#if(PRODUCTION<=1)
extern int error_check(int wherefrom);
#endif

extern int find_horizon(int fromwhere);

/// initialize DUMP stuff
extern int init_dumps(void);
extern void output_nprlist_info(void);
extern void init_dnumcolumns_dnumversion(void);



extern int init_linklists(void);
int setuplinklist(int numcolumns,int which);
extern struct blink * addlink(struct blink * clinkptr);


extern void report_systeminfo(FILE * fileout);
extern int IsLittleEndian(void);
extern void *SwapEndian(void* Addr, const int Nb);

extern void makedirs(void);


/// DUMP file stuff
extern int isenoughfreespace(unsigned long long need);



#include "initbase.funcdeclare.h"


extern int higherorder_set(int whichquantity, int recontype, int*weightsplittype);

extern int get_fluxpldirs(int *Nvec, int dir, int *fluxdir, int* pldir, int *plforflux, FTYPE *signflux);
extern void get_odirs(int dir,int *odir1,int *odir2);
extern int set_location_fluxasemforvpot(int dir, int *numdirs, int *odir1, int *odir2, int *loc);
extern int get_numdirs_fluxasemforvpot(int *numdirs, int *fieldloc);

extern int plstart_set(int whichquantity, int dir, int recontype, int *plstart);



#include "set_grid.funcdeclare.h"





extern FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 );

extern void interpfun(int interptype, int numpoints, int i, FTYPE pos, FTYPE *xfun, FTYPE *fun, FTYPE *answer);





/// log file stuff
#include "mpi_fprintfs.funcdeclare.h"


#include "bounds.funcdeclare.h"

#include "transforms.funcdeclare.h"

#include "coord.funcdeclare.h"

#include "metric.funcdeclare.h"

#include "eos.funcdeclare.h"

#include "phys.funcdeclare.h"
#include "phys.tools.funcdeclare.h"

#include "utoprimgen.funcdeclare.h"




#include "tetrad.funcdeclare.h"


//extern void SHOULDNOTREACHHEREEVERBUGYOUHAVE(void);

#include "wavespeeds.funcdeclare.h"

#include "fixup.funcdeclare.h"

#include "diag.funcdeclare.h"





/// interpolation stuff
extern int get_loop(int pointorlinetype, int interporflux, int dir, struct of_loop *loop);
extern int set_interpalltypes_loop_ranges(int pointorlinetype, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);


/// line types:
extern void set_interp_loop_gen(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
//extern void set_interp_loop(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
//extern void set_interp_loop_expanded(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);

/// point types:
extern int set_interppoint_loop_ranges(int interporflux, int dir, int loc, int continuous, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern int set_interppoint_loop_ranges_3Dextended(int interporflux, int loc, int continuous, int *maxis, int *maxie, int *maxjs, int *maxje, int *maxks, int *maxke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_ranges_2D_EMF_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_ranges_geomcorn_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);

extern void set_interppoint_loop(int interporflux, int dir, int loc, int continuous, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_expanded(int interporflux, int dir, int loc, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_expanded_face2cent(int interporflux, int dir, int loc, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);

#include "fluxvpot.funcdeclare.h"


/// functions for loop stuff
extern void  setup_nprlocalist(int whichprimtype, int *nprlocalstart, int *nprlocalend,int *nprlocallist, int *numprims);

#include "nrutil.funcdeclare.h"


/////////////////////////////////
/////
///// specialty functions
/////
/////////////////////////////////
//extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
//   FTYPE *Edot);
//extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
//     FTYPE rs, FTYPE urs);
//extern void timestep(FTYPE ndtr, FTYPE ndth);
//extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);
//
//extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
//     FTYPE r, FTYPE rs, FTYPE urs);
//extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
//   FTYPE *Urs, FTYPE *Edot);
//extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
//

#include "metric.tools.funcdeclare.h"



extern FTYPE sign_bad(FTYPE a);
extern FTYPE sign_func(FTYPE a);

#ifdef WIN32
// GODMARK: Could refine for a=0
#define sign(a) ((a)>0 ? 1.0 : -1.0)
#else

#if(SUPERLONGDOUBLE)
#define sign(a) (sign_bad(a))
#else
#define sign(a) (copysign(1.0,a))
#endif

#endif



extern FTYPE signavoidzero(FTYPE a);

#ifndef WIN32
extern FTYPE max(FTYPE a, FTYPE b);

extern FTYPE min(FTYPE a, FTYPE b);
#endif

// supplemental trig functions
extern FTYPE mysign(FTYPE x);
extern FTYPE myfabs(FTYPE x);

extern FTYPE mysin(FTYPE th);
extern FTYPE mycos(FTYPE th);

#if(SUPERLONGDOUBLE==0)
extern FTYPE cot(FTYPE arg);
#endif
extern FTYPE csc(FTYPE arg);
extern FTYPE sec(FTYPE arg);

////////////////////////////////////
///
/// SUPERLONGDOUBLE declarations
///
////////////////////////////////////
#if(SUPERLONGDOUBLE)
#include "mconf.h"
extern long double ceill ( long double );
extern long double floorl ( long double );
extern long double atan2l ( long double, long double );
extern int signbitl ( long double );
//
extern long double fabsl ( long double );
extern long double sqrtl ( long double );
extern long double cbrtl ( long double );
extern long double expl ( long double );
extern long double logl ( long double );
extern long double tanl ( long double );
extern long double atanl ( long double );
extern long double sinl ( long double );
extern long double asinl ( long double );
extern long double cosl ( long double );
extern long double acosl ( long double );
extern long double powl ( long double, long double );
extern long double tanhl ( long double );
extern long double atanhl ( long double );
extern long double sinhl ( long double );
extern long double asinhl ( long double );
extern long double coshl ( long double );
extern long double acoshl ( long double );
extern long double exp2l ( long double );
extern long double log2l ( long double );
extern long double exp10l ( long double );
extern long double log10l ( long double );
extern long double gammal ( long double );
extern long double lgaml ( long double );
extern long double jnl ( int, long double );
extern long double ynl ( int, long double );
extern long double ndtrl ( long double );
extern long double ndtril ( long double );
extern long double stdtrl ( int, long double );
extern long double stdtril ( int, long double );
extern long double ellpel ( long double );
extern long double ellpkl ( long double );
long double lgammal(long double);
extern int isfinitel ( long double );
#define finite(arg) isfinitel(arg)
//#define isfinite(arg) isfinitel(arg)
#define copysign( a, b ) ( fabsl(a) * sign(b) )
extern int merror;
#else



#include <math.h>

#ifdef WIN32
#define finite(arg) _finite(arg)
#define isfinite(arg) _finite(arg)
#endif

#ifndef WIN32
//#if USINGICC==0
#if(SUPERLONGDOUBLE==0 && !defined(isfinite))
// needed for Sauron
#define isfinite(arg) finite(arg)  //atch -- on mako, in force-free it would complain about multiply-defined __finite() if not include this line
#endif
#endif // end if not defined WIN32

#endif


#ifdef WIN32
#define copysign( a, b ) ( fabs(a) * sign(b) )
#endif




#if(!DO_ASSERTS)
#define assert assert_func_empty
#else
#define assert assert_func
#endif


extern int assert_func( int is_bad_val, char *s, ... );
extern int assert_func_empty( int is_bad_val, char *s, ... );


#include "global.funcdeclare.rad.h"

#include "global.funcdeclare.user.h"
