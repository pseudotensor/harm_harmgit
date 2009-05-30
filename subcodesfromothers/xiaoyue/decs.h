
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/** here are the few things that we change frequently **/

#define N1	256	/* number of zones */
#define N2	1	/* number of zones */
#define NGhost  3
/** MNEMONICS SECTION **/

/* boundary condition mnemonics */
/*
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3
*/
/* mnemonics for primitive vars; conserved vars */
#define RHO	0	
#define UU	1
#define U1	2
#define U2	3
#define U3	4
#define B1	5
#define B2	6
#define B3	7

/* mnemonics for dimensional indices */
#define TT	0	
#define RR	1
#define TH	2
#define PH	3

/* mnemonics for centering of grid functions */
#define FACE1	0	
#define FACE2	1
#define CORN	2
#define CENT	3

/* mnemonics for slope limiter */
#define MC	0  //regular MC
#define VANL	1  //regular van Leer
#define MINM	2  //regular MINM
#define PMC     3  //limiter used in PPM paper
#define NLIM    4  //no limiter case

/* choices for algorithm */
#define FLUXCT	1
#define FLUXCD	0

/* mnemonics for diagnostic calls */
#define INIT_OUT	0
#define DUMP_OUT	1
#define IMAGE_OUT	2
#define LOG_OUT		3
#define FINAL_OUT	4

/* cooling on or off? */
#define COOLING 0

/** GLOBAL ARRAY SECTION **/

/* size of global arrays */
#define NPR	8	/* number of primitive variables */
#define NDIM	4	/* number of total dimensions.  Never changes */
#define NPG	4	/* number of positions on grid for grid functions */

extern double   a_p[N1+2*NGhost][N2+2*NGhost][NPR] ;	/* space for primitive vars */
extern double  a_dq[N1+2*NGhost][N2+2*NGhost][NPR] ;  /* slopes */
extern double  a_F1[N1+2*NGhost][N2+2*NGhost][NPR] ;	/* fluxes */
extern double  a_F2[N1+2*NGhost][N2+2*NGhost][NPR] ;	/* fluxes */

/*parabolic*/
extern double  a_ph1[N1+2*NGhost][N2+2*NGhost][NPR] ;	/* half-step primitives */
extern double a_ph2[N1+2*NGhost][N2+2*NGhost][NPR];
extern double a_ph3[N1+2*NGhost][N2+2*NGhost][NPR] ;
extern double a_ph4[N1+2*NGhost][N2+2*NGhost][NPR] ;
extern double a_uh[N1+2*NGhost][N2+2*NGhost][NPR] ;
extern double a_Ip_l[N1+2*NGhost][N2+2*NGhost][NPR];
extern double a_Ip_r[N1+2*NGhost][N2+2*NGhost][NPR];

/* for debug */
extern double ivar[N1][N2][NPR] ;
extern short  stat[N1][N2] ;
extern double psave[N1][N2][NPR] ;

/* grid functions */
extern double a_conn[N1+2*NGhost][N2+2*NGhost][NDIM][NDIM][NDIM] ;
extern double a_gcon[N1+2*NGhost][N2+2*NGhost][NPG][NDIM][NDIM] ;
extern double a_gcov[N1+2*NGhost][N2+2*NGhost][NPG][NDIM][NDIM] ;
extern double a_gdet[N1+2*NGhost][N2+2*NGhost][NPG] ;

extern double (*   p)[N2+2*NGhost][NPR] ;
extern double (*  dq)[N2+2*NGhost][NPR] ;
extern double (*  F1)[N2+2*NGhost][NPR] ;
extern double (*  F2)[N2+2*NGhost][NPR] ;
extern double (*  ph1)[N2+2*NGhost][NPR] ;
extern double (*  ph2)[N2+2*NGhost][NPR] ;
extern double (*  ph3)[N2+2*NGhost][NPR] ;
extern double (*  ph4)[N2+2*NGhost][NPR] ;
extern double (*  uh)[N2+2*NGhost][NPR] ;
extern double (*  Ip_l)[N2+2*NGhost][NPR] ;
extern double (*  Ip_r)[N2+2*NGhost][NPR] ;

extern double (* conn)[N2+2*NGhost][NDIM][NDIM][NDIM] ;
extern double (* gcon)[N2+2*NGhost][NPG][NDIM][NDIM] ;
extern double (* gcov)[N2+2*NGhost][NPG][NDIM][NDIM] ;
extern double (* gdet)[N2+2*NGhost][NPG] ;

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern double a ;
extern double gam ;

/* numerical parameters */
extern double L;
extern double Rin,Rout,hslope ;
extern double cour ;
extern double dV,dx[NPR],startx[NPR] ;
extern double dt ;
extern double t,tf ;
extern double x1curr,x2curr ;
extern int nstep ;

/* output parameters */
extern double DTd ;
extern double DTl ;
extern double DTi ;
extern int    DTr ;
extern int    dump_cnt ;
extern int    image_cnt ;
extern int    rdump_cnt ;
extern int    nstroke ;

/* global flags */
extern int failed ;
extern int lim ;
extern double defcon ;

/* diagnostics */
extern double mdot ;
extern double edot ;
extern double ldot ;

/* numerical convenience */
#define SMALL	1.e-30 

/** MACROS **/
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SIGN(a) ( ((a) <0.) ? -1. : 1. )


/* loop over all active zones */
#define ZLOOP for(i=0;i<N1;i++)for(j=0;j<N2;j++)

/* specialty loop */
extern int istart,istop,jstart,jstop ;
#define ZSLOOP(istart,istop,jstart,jstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

/* loop over Primitive variables */
#define PLOOP  for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)

/* set global variables that indicate current local metric, etc. */
extern int icurr,jcurr,pcurr ;
struct of_geom {
	double gcon[NDIM][NDIM] ;
	double gcov[NDIM][NDIM] ;
	double g ;
} ;

struct of_state {
	double ucon[NDIM] ;
	double ucov[NDIM] ;
	double bcon[NDIM] ;
	double bcov[NDIM] ;
} ;

extern double dminarg1,dminarg2; 
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

#define delta(i,j) ( (i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]) 

/* size of step in numerical derivative evaluations */
#define HSTEP	1.e-5

/** FIXUP PARAMETERS **/
#define RHOMIN	1.e-4
#define UUMIN	1.e-6

#define MAXBSQOVERRHO	1.e8
#define MAXBSQOVERUU	1.e8

/* maximum fractional increase in timestep per timestep */
#define SAFE	1.3

/* function declarations */
/** general use **/
double bl_gdet_func(double r, double th) ;
double bsq_calc(double *pr, struct of_geom *geom) ;
double gdet_func(double lgcov[][NDIM]) ;
double mink(int j, int k) ; 
double ranc(int seed) ;
double slope_lim(double y1, double y2, double y3) ;

int restart_init(void) ;

void area_map(int i, int j, double prim[][N2+2*NGhost][NPR]) ;
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) ;
void blgset(int i, int j, struct of_geom *geom);
void bltoks(double *pr, int i, int j) ;
void bl_coord(double *X, double *r, double *th) ;
void bl_gcon_func(double r, double th, double gcov[][NDIM]) ;
void bl_gcov_func(double r, double th, double gcov[][NDIM]) ;
void bound_prim(double pr[][N2+2*NGhost][NPR]) ;
void conn_func(double *X, struct of_geom *geom, double lconn[][NDIM][NDIM]) ;
void coord(int i, int j, int loc, double *X) ;
void diag(int call_code) ;
void diag_flux(double F1[][N2+2*NGhost][NPR], double F2[][N2+2*NGhost][NPR]) ;
void dump(FILE *fp) ;
void dudp_calc(double *pr, struct of_state *q, struct of_geom *geom,
		double **alpha) ;
void fail(int fail_type) ;
void fixup(double (* var)[N2+2*NGhost][NPR]) ;
void flux_ct(double F1[][N2+2*NGhost][NPR],double F2[][N2+2*NGhost][NPR]) ;
void gaussj(double **tmp, int n, double **b, int m) ;
void gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]) ;
void gcov_func(double *X, double lgcov[][NDIM]) ;
void get_geometry(int i, int j, int loc, struct of_geom *geom) ;
void get_state(double *pr, struct of_geom *geom, struct of_state *q) ;
void image(FILE *fp) ;
void init(void) ;
void kstoksp(double *pr, int i, int j) ;
void lower(double *a, struct of_geom *geom, double *b) ;
void ludcmp(double **a, int n, int *indx, double *d) ;
void lubksb(double **a, int n, int *indx, double *d) ;
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)  ;
void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);
void primtoflux(double *pa, struct of_state *q, int dir, struct of_geom *geom,
			double *fl) ;
void primtoU(double *p, struct of_state *q, struct of_geom *geom, double *U) ;
void raise(double *v1, struct of_geom *geom, double *v2) ;
void restart_write(void) ;
void restart_read(FILE *fp) ;
void set_arrays(void) ;
void set_grid(void) ;
void set_points(void) ;
void step_ch(void) ;
void source(double *pa, struct of_geom *geom, int ii, int jj, double *Ua) ;
void timestep(void) ;
void u_to_v(double *pr) ;
void ucon_calc(double *pr, struct of_geom *geom, double *ucon) ;
void usrfun(double *pr,int n,double *beta,double **alpha) ;
void Utoprim(double *Ua, struct of_geom *geom, double *pa) ;
void vchar(double *pr, struct of_state *q, struct of_geom *geom,
		int dir, double *cmax, double *cmin) ;

/* NR routines */
double *dvector(int il, int ih) ;
double **dmatrix(int xl, int xh, int yl, int yh) ;
double ***dtensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh) ;
int *ivector(int il, int ih) ;
int mnewt(int ntrail, double *p, int n, double tolx, double tolf) ;
void free_dvector(double *dvec, int il, int ih) ;
void free_dmatrix(double **dmat, int xl, int xh, int yl, int yh) ;
void free_dtensor(double ***t, int nrl,int nrh,int ncl,int nch,
		int ndl,int ndh) ;
void nrerror(char error_text[]) ;
void free_ivector(int *ivec, int il, int ih) ;

/* failure modes */
#define FAIL_UTOPRIM_NEG	1
#define FAIL_UTOPRIM_TEST	2
#define FAIL_VCHAR_DISCR	3
#define FAIL_COEFF_NEG		4
#define FAIL_COEFF_SUP		5
#define FAIL_UTCALC_DISCR	6

