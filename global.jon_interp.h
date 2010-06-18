#define DEBUGINTERP 0 // debug messages

// normal failure to interpolate message
#define SIMPLEDEBUGINTERP 0

// whether to reduce to nearest neighbor at boundaries
// 1: do
// 0: shift points so still interpolate
#define BILINEARREDUCE2NEARESTATBOUNDARY 0

// Note: Due to how generally create data as multiple time dumps, easier to read and write time as slowest index even if by N? it's before i
#define LOOPOLDDATA for(h=0;h<oN0;h++) for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)    for(i=0;i<oN1;i++)
#define LOOPINTERP for(h=0;h<nN0;h++) for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)    for(i=0;i<nN1;i++)

// oldgridtype: 0=Cartesian  1=spherical polar 2=cylindrical 3=log(z) vs. log(R)// V in GRMHD code 4 = log for radius (used in Sashas monopole paper) 5=Cartesian, but time is mixed with space to approximate light travel time effects [only makes sense if doing 4D input with oN0>1]
// newgridtype: -1=nochange (else like oldgridtype) // output coordinate system
#define GRIDTYPENOCHANGE -1
#define GRIDTYPECART 0
#define GRIDTYPESPC 1
#define GRIDTYPECYL 2
#define GRIDTYPELOGLOGCYL 3
#define GRIDTYPELOGSPC 4
#define GRIDTYPECARTLIGHT 5

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
#include <time.h>

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

#undef N0
#undef N1
#undef N2
#undef N3
#define N0 2
#define N1 2
#define N2 2
#define N3 2

#define ROUND2INT(x) ((int)((x)>0.0 ? (x)+0.5 : (x)-0.5))

// number of dimensions for code (2d or 3d code)
#define COMPDIM 3 // not important for time (oN0 and nN0>1) stuff

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
// 6 more args
#define SCANHEADER "%f %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %d %d %d %d %d %d %d %d %d"
#elif(REALTYPE==DOUBLETYPE)
#define SCANARG "%lf"
#define SCANARGVEC "%lf %lf %lf"
#define SCANARG4VEC "%lf %lf %lf %lf"
#define SCANHEADER "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %d %d %d %d %d %d %d %d"
#elif(REALTYPE==LONGDOUBLETYPE)
#define SCANARG "%Lf"
#define SCANARGVEC "%Lf %Lf %Lf"
#define SCANARG4VEC "%Lf %Lf %Lf %Lf"
#define SCANHEADER "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf %d %d %d %d %d %d %d %d %d"
#endif

#define PRINTSCANHEADER "%g %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g %d %d %d %d %d %d %d %d %d\n"















#include "nrutil.funcdeclare.h"


// below mnewt() and fdjac() are not like used in HARM code
// functions inside jon_interp_mnewt.c or used by it (i.e. usrfun() )
extern int mnewt(int ntrial, int mintrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf, FTYPE tolxallowed, FTYPE tolfallowed, FTYPE tolxreport, FTYPE tolfreport, FTYPE *parms, int (*usrfun)(int n, FTYPE *, FTYPE *, FTYPE *, FTYPE **, FTYPE*));

extern void fdjac(int n, FTYPE parms[], FTYPE x[], FTYPE fvec[], FTYPE **df,
		  void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]));


#include "coord.funcdeclare.h"




extern void dxdxp_func(FTYPE *X, FTYPE (*dxdxp)[NDIM]);

extern int setihor(void);
extern FTYPE setRin(int ihor);
extern FTYPE rhor_calc(int which);

#include "mpi_fprintfs.funcdeclare.h"

//#include "phys.tools.funcdeclare.h"
#include "metric.tools.funcdeclare.h"

#include "tetrad.funcdeclare.h"

// fake fillers:
extern void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcovinfunc,FTYPE *gcovpertinfunc);
extern void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon);
extern void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *EOMFUNCNAME);
extern void assign_eomfunc(struct of_geom *geom, FTYPE *EOMFUNCNAME);




// various jon_interp*.c files:
extern void copy_old2new(void);
extern int compute_spatial_interpolation(void);

extern void gaussian_filter(int filter,FTYPE sigma, int nt, int nx, int ny, int nz, unsigned char****oldimage,FTYPE****olddata);


extern void writeimage(char * name, unsigned char ****image,int nt, int nx, int ny, int nz);

extern void refine_data(void);

extern void compute_preprocess(FILE *gdumpin, FTYPE *finaloutput);


extern void setup_newgrid(void);





