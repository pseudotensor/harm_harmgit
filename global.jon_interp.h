#define VERSION "1.0"


// 0: Use shifted stencil for higher-order interpolations (This can lead to extrapolating positive quantity into a negative one)
// 1: Fill in boundary zones and don't shift stencil (safest)
#define BOUNDARYEXTRAP 1

// whether to also rotate z-axis towards observer
#define ROTATECARTLIGHT 1


// whether to reduce to nearest neighbor at boundaries
// 1: do
// 0: shift points so still interpolate
#define BILINEARREDUCE2NEARESTATBOUNDARY 0

// whether to:
// 1= allocate memory for newimage or newdata
// or 0 = just output to file the interpolation results directly (can save huge on memory if doing large interpolations)
#define ALLOCATENEWIMAGEDATA 0


// Note: Due to how generally create data as multiple time dumps, easier to read and write time as slowest index even if by N? it's before i
#define LOOPOLDDATA for(h=0;h<oN0;h++) for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)    for(i=0;i<oN1;i++)
#define LOOPOLDDATASPATIAL for(k=0;k<oN3;k++) for(j=0;j<oN2;j++) for(i=0;i<oN1;i++)

///////
// COLMARK: When finally writing the file, we place columns as fastest index, i next, j next, k next, and time next
///////
#define OUTPUTSPACETIMELOOP for(h=0;h<nN0;h++) for(k=0;k<nN3;k++) for(j=0;j<nN2;j++) for(i=0;i<nN1;i++)
#define OUTPUTLOOP OUTPUTSPACETIMELOOP for(coli=0;coli<numoutputcols;coli++)

#if(ALLOCATENEWIMAGEDATA==1)
// choice
#define LOOPINTERP for(h=0;h<nN0;h++) for(k=0;k<nN3;k++) for(j=0;j<nN2;j++)    for(i=0;i<nN1;i++)
#else
// NO choice since have to loop during interpolation just as would output to file so can avoid allocating memory but write to file correctly.
#define LOOPINTERP OUTPUTSPACETIMELOOP
#endif





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

#define FLOAT2IMAGE(x) ((unsigned char)(round(FLOAT2IMAGE1(x))))

#define sign(a) (copysign(1.0,a))

// maximum number of columns for memory allocation of some things
#define MAXCOLS 40

#define MAXINCOLS 30

#define MAXBC 2
// whether to treat \phi direction as periodic (used when N3>1 and oldgridtype==1)
#define PERIODICINPHI 1 // default: 1

#define SLOOP12(j,k) for(j=1;j<NDIM-1;j++)for(k=1;k<NDIM-1;k++)


// field line file quantity list for input file (Sasha has more than 10, but not required for present calculations)
// numcolumns - colini type
#define FLRHO 0
#define FLU 1
#define FLNEGUD0 2
#define FLMUCONST 3
#define FLU0 4
#define FLV1 5
#define FLV2 6
#define FLV3 7
#define FLB1 8
#define FLB2 9
#define FLB3 10
#define FLURAD0 11
#define FLVRAD1 12
#define FLVRAD2 13
#define FLVRAD3 14

// FL's to skip when storing
#define SKIPFL2STORE(x) (x==FLNEGUD0 || x==FLMUCONST)

// only required quantities list for compute_preprocessing() or compute_additionals()
// numcolumnsstorage - colstorei type
#define NUMCOLUMNSSTORE (2+NDIM+(NDIM-1))
#define STORERHO 0
#define STOREU 1
#define STOREU0 2
#define STOREV1 3
#define STOREV2 4
#define STOREV3 5
#define STOREB1 6
#define STOREB2 7
#define STOREB3 8


// output to interpolation and eventually written to disk
// numoutputcols - coli type list
#define OUTRHO 0
#define OUTU 1
#define OUTV1 2
#define OUTV2 3
#define OUTV3 4
#define OUTB1 5
#define OUTB2 6
#define OUTB3 7
#define OUTFEMRAD 8
#define OUTBPHI 9
#define OUTJX0  10
#define OUTJX1  11
#define OUTJX2  12
#define OUTJX3  13



#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <string.h>

#include "global.realdef.h"

#include "coord.h"
#include "global.nondepmnemonics.h"
#include "definit.h"

#define POSDEFMETRIC 0
#define USEMPI 0
#define USEOPENMP 0


// let init.h define this until add MCOORD to header of input files.
// describes original metric.  Only used for coord.c polar axis issue.
//#undef MCOORD
//#define MCOORD CYLMINKMETRIC          // coordinates for some metric
//#define MCOORD CARTMINKMETRIC2          // coordinates for some metric
//#define MCOORD SPCMINKMETRIC

#undef DOSTOREPOSITIONDATA
#define DOSTOREPOSITIONDATA 0

#undef DOGRIDSECTIONING
#define DOGRIDSECTIONING 0

#undef STORETLAB2ORTHO
#define STORETLAB2ORTHO 0

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
#define logdt_file stderr
#define myexit exit





#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )









#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
//#define SIGN(a) ( ((a) <0.) ? -1. : 1. )

#define SLEPSILON (1.e-6)




#include "global.depmnemonics.h"
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


// whether to use new THETAROT header or not
// 2: very old header (21 things: runlocaldipole3dfiducial)
// 1: old header
// 0: new header
//#define OLDERHEADER 2
// NOW this is input variable defaulted to 0


// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define SCANARG "%f"
#define SCANARGVEC "%f %f %f"
#define SCANARG4VEC "%f %f %f %f"
#define SCANFIELDLINE "%f %f %f %f %f %f %f %f %f %f %f" // 11 items
#define SCANHEADER1      " %f %d %d %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %d  %f  %f  %d %d %d %d %d %d %d %d %d" // 30
#define SCANHEADER2      " %f %d %d %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %d  %f  %f" // 21
#define SCANHEADER0      " %f %d %d %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %d  %f  %f  %f  %f %d %d %d %d %d %d %d %d %d" // 32
#elif(REALTYPE==DOUBLETYPE)
#define SCANARG "%lf"
#define SCANARGVEC "%lf %lf %lf"
#define SCANARG4VEC "%lf %lf %lf %lf"
#define SCANFIELDLINE "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" // 11 items
#define SCANHEADER1      "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %d %d %d %d %d %d %d %d" // 30
#define SCANHEADER2      "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf" // 21
#define SCANHEADER0      "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d" // 32
#elif(REALTYPE==LONGDOUBLETYPE)
#define SCANARG "%Lf"
#define SCANARGVEC "%Lf %Lf %Lf"
#define SCANARG4VEC "%Lf %Lf %Lf %Lf"
#define SCANFIELDLINE "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf" // 11 items
#define SCANHEADER1      "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf %d %d %d %d %d %d %d %d %d" // 30
#define SCANHEADER2      "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf" // 21
#define SCANHEADER0      "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf %Lf %Lf %d %d %d %d %d %d %d %d %d" // 32
#endif



#define PRINTSCANHEADER1 "%g  %d %d %d %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %d %g  %g  %d %d %d %d %d %d %d %d %d\n" // 30
#define SCANHEADERARGS1 &tdump,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dX[1],&dX[2],&dX[3],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dtdump,&defcoord,&MBH,&QBH,&is,&ie,&js,&je,&ks,&ke,&whichdump,&whichdumpversion,&numcolumns
#define PRINTHEADERARGS1 tdump,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],readnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numoutputcols

#define PRINTHEADERSTDERR1 "OLD: %22.16g :: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g :: %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n" // 30

#define PRINTHEADERSTDERRARGS1 tdump,oN1,oN2,oN3,startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],realnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns

#define PRINTHEADERSTDOUT1 "%22.16g %d %d %d %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n" // 30

#define PRINTHEADERSTDOUTARGS1 tdump, nN1, nN2, nN3, startxc, startyc, startzc, fakedxc,fakedyc,fakedzc,realnstep,gam,spin,ftemp,endxc,endyc,hslope,dtdump,defcoord,MBH,QBH,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numoutputcols




#define PRINTSCANHEADER2 "%g  %d %d %d %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %d %g  %g\n" // 21
#define SCANHEADERARGS2 &tdump,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dX[1],&dX[2],&dX[3],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dtdump,&defcoord,&MBH,&QBH
#define PRINTHEADERARGS2 tdump,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],readnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH

#define PRINTHEADERSTDERR2 "OLD2: %22.16g :: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g :: %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g\n" // 21

#define PRINTHEADERSTDERRARGS2 tdump,oN1,oN2,oN3,startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],realnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH

#define PRINTHEADERSTDOUT2 "%22.16g %d %d %d %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g\n" // 21

#define PRINTHEADERSTDOUTARGS2 tdump, nN1, nN2, nN3, startxc, startyc, startzc, fakedxc,fakedyc,fakedzc,realnstep,gam,spin,ftemp,endxc,endyc,hslope,dtdump,defcoord,MBH,QBH



#define PRINTSCANHEADER0 "%g  %d %d %d %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %d %g  %g  %g  %g  %d %d %d %d %d %d %d %d %d\n" // 32
#define SCANHEADERARGS0 &tdump,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dX[1],&dX[2],&dX[3],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dtdump,&defcoord,&MBH,&QBH,&EP3,&THETAROT,&is,&ie,&js,&je,&ks,&ke,&whichdump,&whichdumpversion,&numcolumns

#define PRINTHEADERARGS0 tdump,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],readnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numoutputcols
// 32 -- should be same as dump_header_general() in dump.c

#define PRINTHEADERSTDERR0 "OLD: %22.16g :: %d %d %d :: %22.16g %22.16g %22.16g :: %22.16g %22.16g %22.16g :: %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n" // 32

#define PRINTHEADERSTDERRARGS0 tdump,oN1,oN2,oN3,startx[1],startx[2],startx[3],dX[1],dX[2],dX[3],realnstep,gam,spin,R0,Rin,Rout,hslope,dtdump,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns

#define PRINTHEADERSTDOUT0 "%22.16g %d %d %d %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %ld %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %22.16g %d %22.16g %22.16g %22.16g %22.16g %d %d %d %d %d %d %d %d %d\n" // 32

#define PRINTHEADERSTDOUTARGS0 tdump, nN1, nN2, nN3, startxc, startyc, startzc, fakedxc,fakedyc,fakedzc,realnstep,gam,spin,ftemp,endxc,endyc,hslope,dtdump,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numoutputcols














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

extern void gaussian_filter(int filter,FTYPE sigma, int nt, int nx, int ny, int nz, unsigned char*****oldimage,FTYPE*****olddata);


extern void writeimage(char * name, unsigned char *****image,int nt, int nx, int ny, int nz);

extern void refine_data(void);

extern void compute_preprocess(int outputvartypelocal, FILE *gdumpfile, int h, int i, int j, int k, FTYPE*****olddatalocal, FTYPE *finaloutput);



extern void setup_newgrid(void);


extern void interp_bl_coord(FTYPE *X, FTYPE *V);



extern void readelement(int binaryinputlocal, char* inFTYPElocal, FILE *input, FTYPE *datain);
extern void writeelement(int binaryoutputlocal, char* outFTYPElocal, FILE *output, FTYPE dataout);

extern long sizeelement(char* inFTYPElocal);


extern void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom);


extern void compute_gdetFuu(FTYPE gdet, FTYPE *gcov, FTYPE *ucon, FTYPE *Bcon, FTYPE (*Fuu)[NDIM]);

extern void compute_simple_gdetFuu(FTYPE gdet, FTYPE *gcov, FTYPE *ucon, FTYPE *Bcon, FTYPE (*Fuu)[NDIM]);


extern int compute_additionals(void);

extern void apply_boundaryconditions_olddata(int numcols, int oN0local, int numbc0local, int doubleworklocal, unsigned char *****oldimagelocal, FTYPE *****olddatalocal);
extern void apply_boundaryconditions_olddata_cleanpole(int numcols, int oN0local, int numbc0local, int doubleworklocal, unsigned char *****oldimagelocal, FTYPE *****olddatalocal);

extern void gdump_tostartofdata(FILE **gdumpinlocal);
extern void infile_tostartofdata(FILE **infilelocal);

extern void output2file_perpointcoli_postinterpolation(unsigned char newimagelocal, FTYPE newdatalocal);

extern void output2file_perpoint_postinterpolation(int which, int h, int i, int j, int k, unsigned char *newimagelocal, FTYPE *newdatalocal);

extern void raise_vec(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon);

extern void lower_vec(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov);
