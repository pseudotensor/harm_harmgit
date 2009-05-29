#define DEBUGINTERP 0 // debug messages

// normal failure to interpolate message
#define SIMPLEDEBUGINTERP 0

// whether to reduce to nearest neighbor at boundaries
// 1: do
// 0: shift points so still interpolate
#define BILINEARREDUCE2NEARESTATBOUNDARY 0

#define LOOPOLDDATA for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)    for(i=0;i<oN1;i++)
#define LOOPINTERP for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)    for(i=0;i<nN1;i++)

// oldgridtype: 0=Cartesian  1=spherical polar 2=cylindrical 3=log(z) vs. log(R)// V in GRMHD code
// newgridtype: -1=nochange (else like oldgridtype) // output coordinate system
#define GRIDTYPENOCHANGE -1
#define GRIDTYPECART 0
#define GRIDTYPESPC 1
#define GRIDTYPECYL 2
#define GRIDTYPELOGLOGCYL 3
#define GRIDTYPELOGSPC 4

#define FLOAT2IMAGE1(x) ( (x<0.0) ? 0.0 : (x>255.0 ? 255.0 : x) )

#define FLOAT2IMAGE(x) ((unsigned char)(round(FLOAT2IMAGE1(ftemp))))

#define sign(a) (copysign(1.0,a))


#define MAXBC 2
// whether to treat \phi direction as periodic (used when N3>1 and oldgridtype==1)
#define PERIODICINPHI 1 // default: 1


#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdarg.h>

#include "global.realdef.h"

#include "coord.h"
#include "global.nondepnmemonics.h"
#include "definit.h"

#define POSDEFMETRIC 0
#define USEMPI 0
#define USEOPENMP 0



#undef MCOORD
//#define MCOORD CYLMINKMETRIC          // coordinates for some metric
#define MCOORD SPCMINKMETRIC

#undef DOSTOREPOSITIONDATA
#define DOSTOREPOSITIONDATA 0

#undef DOGRIDSECTIONING
#define DOGRIDSECTIONING 0

#undef N1
#undef N2
#undef N3
#define N1 2
#define N2 2
#define N3 2

#define ROUND2INT(x) ((int)((x)>0.0 ? (x)+0.5 : (x)-0.5))

// number of dimensions for code (2d or 3d code)
#define COMPDIM 3

#define COORDSINGFIXCYL 0

#define ANALYTICJAC 1 // helps to avoid no solution sometimes -- e.g. when precisionj of numerical jac in question)
// 0: numerical jac
// 1: analytic jac
// applies to newt() and broydn()
#define ROOTMETHOD 2 // 2 has proven best
// 0 : nondamped mnewt() // definitely can fail to find solution even to modest error unless guess is close
// 1 : damped mnewt() // comment is same as 0 above
// 2 : newt() // works with or without good guess, and with or (mostly) without analytic jac --  as long as jacobian is allowed to be unrestricted outside normal domain (see coord.c and singularity/POSDEFMETRIC fix turned off).
// 3 : brodyn() // basically same as comment for 2 above, but can be much worse than 2
#define GOODGUESS 0 // turned on for 2/3 ROOTMETHODS only to speed things up
// 0: simple mid-point guess always
// 1: choose guess with some insight of typical coordinat transformations


#define DOUBLEWORK 1
// 0=use unsigned char during calculations for images (can result in over or under flow)
// 1=use doubles during calculations for images (avoids over/underflow)


//#define dualfprintf fprintf
//#define myfprintf fprintf
#define fail_file stderr
#define logfull_file stderr
#define log_file stderr
#define myexit exit





#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )









#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
//#define SIGN(a) ( ((a) <0.) ? -1. : 1. )

#define SLEPSILON (1.e-6)




#include "global.depnmemonics.h"
#include "global.storage.h"
#include "global.loops.h"

///////////////////
// redefine loops
#undef PLOOP
#undef PDUMPLOOP
#undef PINVERTLOOP
#undef PBOUNDLOOP
#undef PINTERPLOOP

/* loop over all Primitive variables */
#define PLOOP(pliter,pl) for(pl=0;pl<NPR;pl++)
/* loop over all dumped Primitive variables */
#define PDUMPLOOP(pliter,pd) for(pd=0;pd<NPRDUMP;pd++)
/* loop over all inversion Primitive variables */
#define PINVERTLOOP(pliter,pi) for(pi=0;pi<NPRINVERT;pi++)
/* loop over all bounding Primitive variables */
#define PBOUNDLOOP(pliter,pb) for(pb=0;pb<NPRBOUND;pb++)
/* loop over all center to edge variables */
#define PINTERPLOOP(pliter,pl) for(pl=0;pl<NPR2INTERP;pl++)


#include "global.variousmacros.h"
#include "global.fieldmacros.h"
#include "global.structs.h"
#define SIMULBCCALC -1
#include "global.gridsectioning.h"
#include "global.comploops.h"




// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define SCANARG "%f"
#define SCANARGVEC "%f %f %f"
#define SCANARG4VEC "%f %f %f %f"
// 16 args
// 21 args after going to 3D and doing MBH/QBH
#define SCANHEADER "%f %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f"
#elif(REALTYPE==DOUBLETYPE)
#define SCANARG "%lf"
#define SCANARGVEC "%lf %lf %lf"
#define SCANARG4VEC "%lf %lf %lf %lf"
#define SCANHEADER "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf"
#elif(REALTYPE==LONGDOUBLETYPE)
#define SCANARG "%Lf"
#define SCANARGVEC "%Lf %Lf %Lf"
#define SCANARG4VEC "%Lf %Lf %Lf %Lf"
#define SCANHEADER "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf "
#endif

#define PRINTSCANHEADER "%g %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g\n"















///////////////////////////
//
// NR STUFF
//
///////////////////////////
extern FTYPE ranc(int seed);
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);


/* NR routines from nrutil.h used by nrutil.c */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);


extern void free_vector(FTYPE *v, long nl, long nh);
extern FTYPE *vector(long nl, long nh);
extern FTYPE **matrix(long nrl, long nrh, long ncl, long nch);
extern unsigned char **cmatrix(int nrl,int nrh,int ncl,int nch);
extern FTYPE **fmatrix(int nrl,int nrh,int ncl,int nch);
extern FTYPE ***f3matrix(int nzl, int nzh, int nrl, int nrh, int ncl, int nch);
extern void free_cmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch);
extern void free_fmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);

extern void free_f3matrix(FTYPE ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_matrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);

extern void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
extern void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);
extern void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
extern void rotate(FTYPE **r, FTYPE **qt, int n, int i, FTYPE a, FTYPE b);

extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);

extern void lnsrch(int n, FTYPE parms[], FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], FTYPE *f, FTYPE stpmax, int *check, FTYPE (*func)(FTYPE [], FTYPE []));

extern void lubksb(FTYPE **a, int n, int *indx, FTYPE b[]);
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);

extern void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
extern void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
extern void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);


extern FTYPE nrfmin(FTYPE parms[], FTYPE x[]);



extern void newt(int useanalyticjac
		 ,FTYPE parms[]
		 ,FTYPE x[], int n, int *check
		 ,void (*vecfunc)(int n, FTYPE *parms, FTYPE v[], FTYPE f[])
		 ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
		 );
extern void broydn(int useanalyticjac
		   ,FTYPE parms[]
		   ,FTYPE x[], int n, int *check
		   ,void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[])
		   ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
		   );

// functions inside jon_interp_mnewt.c or used by it (i.e. usrfun() )
extern int mnewt(int ntrial, int mintrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf, FTYPE tolxallowed, FTYPE tolfallowed, FTYPE tolxreport, FTYPE tolfreport, FTYPE *parms, int (*usrfun)(int n, FTYPE *, FTYPE *, FTYPE *, FTYPE **, FTYPE*));

extern void fdjac(int n, FTYPE parms[], FTYPE x[], FTYPE fvec[], FTYPE **df,
		  void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]));


// functions inside jon_interp_coord.c (or used by it or for it)
// coordinate stuff
extern void set_coord_parms(int defcoordlocal);
extern void set_coord_parms_nodeps(int defcoordlocal);
extern void set_coord_parms_deps(int defcoordlocal);
extern void write_coord_parms(int defcoordlocal);
extern void read_coord_parms(int defcoordlocal);

#if(COMPDIM==2)
// 2D:
extern void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th);
void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);
extern void coord(int i, int j, int loc, FTYPE *X);
extern void coordf(FTYPE i, FTYPE j, int loc, FTYPE *X);
extern void icoord(FTYPE *X,int loc, int *i, int *j);
#endif

#if(COMPDIM==3)
// 3D:
extern void bl_coord(FTYPE *X, FTYPE *V);
void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void coordf(FTYPE i, FTYPE j, FTYPE k, int loc, FTYPE *X);
extern void icoord(FTYPE *X,int loc, int *i, int *j, int *k);
#endif



extern void coord_ijk(int i, int j, int k, int loc, FTYPE *X);
extern void coord_free(int i, int j, int k, int loc, FTYPE *X);

extern void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V);
extern void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);

extern void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM]);



extern void dxdxp_func(FTYPE *X, FTYPE (*dxdxp)[NDIM]);
extern void set_points();
extern int setihor(void);
extern FTYPE setRin(int ihor);
extern FTYPE rhor_calc(int which);



extern void bcucof(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE d1, FTYPE d2, FTYPE **c);

extern void bcuint(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE x1l,
        FTYPE x1u, FTYPE x2l, FTYPE x2u, FTYPE x1, FTYPE x2, FTYPE *ansy,
	    FTYPE *ansy1, FTYPE *ansy2);


// mpi_fprintfs.c:
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);
extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);


// metric_tools.c:
extern void matrix_inverse(int whichcoord, FTYPE (*generalmatrixlower)[NDIM], FTYPE (*generalmatrixupper)[NDIM]);
extern FTYPE rmso_calc(int which);



// tetrad.c:
extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]);
extern int tetr_func_frommetric(FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);




// various jon_interp*.c files:
extern void copy_old2new(void);
extern int compute_spatial_interpolation(void);

extern void gaussian_filter(int filter,FTYPE sigma,int nx, int ny, int nz, unsigned char***oldimage,FTYPE***olddata);

extern void writeimage(char * name, unsigned char ***image,int nx, int ny, int nz);

extern void refine_data(void);

extern void compute_preprocess(FILE *gdumpin, FTYPE *finaloutput);


extern void setup_newgrid(void);


// memory stuff
extern unsigned char **cmatrix(int a, int b, int c, int d)  ;
extern unsigned char ***c3matrix(int a, int b, int c, int d, int e, int f)  ;
extern FTYPE **fmatrix(int a, int b, int c, int d)  ;
extern FTYPE ***f3matrix(int a, int b, int c, int d, int e, int f)  ;




