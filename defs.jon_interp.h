#include "global.jon_interp.h"

#include "rancdefs.h"

int calledranc;
FTYPE NUMEPSILONPOW23;
FTYPE Xmax[NDIM];

FTYPE Risco,drsing;
FTYPE *rcent,*rcent_tot;
int horizoni,horizoncpupos1;
int numprocs;
int mycpupos[NDIM];
int ncpux1;

// SOME GEOMETRIC VARIABLES (see also global.jon_interp.h for PERIODICINPHI)
int dofull2pi;
int whichdump,whichdumpversion,numcolumns;

int DEBUGINTERP; // detailed (somewhat arbitrary) debug messages

// normal failure to interpolate message
int SIMPLEDEBUGINTERP;

int VERBOSITY;


int oN0,oN1,oN2,oN3,nN0,nN1,nN2,nN3 ;
FTYPE refinefactor;
int roN0,roN1,roN2,roN3;
FTYPE dX[NDIM],Rin,fakeRin,dtc,dxc,dyc,dzc,fakedtc,fakedxc,fakedyc,fakedzc;
FTYPE Zin,Zout,Zeqin;
FTYPE starttc, startxc, startyc, startzc, endtc, endxc, endyc, endzc;
int newgridtype,oldgridtype;
FTYPE gridAAglobal,gridr0global;

FTYPE X[NDIM];
FTYPE Xmetricnew[NDIM],Xmetricold[NDIM]; // used to store time of latest and oldest metric
FTYPE endtdata,starttdata; // for 4D dump inputs


FTYPE tdump,gam,spin,QBH,MBH; // tdump used to be t, like it is in HARM, but now t is used locally for 4D interpolation
int startpos[NDIM];
int totalsize[NDIM];
long realnstep,nstep;
FTYPE readnstep;
int is,ie,js,je,ks,ke;

// number of boundary cells for \phi interpolation
int totalbc,numbc[NDIM];



// NUMREC STUFF
int nn;
FTYPE *fvec;


// 0=image
// 1=data
int DATATYPE;
int INTERPTYPE;
// 0: nearest
// 1: bi-linear (true distances)
// 2: planar (true distances)
// 3: bicubic

int immediateoutput,num4vectors;
int outputvartype; 
// vector component: 0=scalar, 1,2,3 
int vectorcomponent;
int defaultvaluetype,EXTRAPOLATE;
FTYPE totalmin,totalmax;
FTYPE defaultvalue;
int didrefine;
int filter;
FTYPE sigma;
int imagedata;

/////////////////////////////
// file stuff
int READHEADER;
int WRITEHEADER;
int jonheader;
int getgdump;
char gdumpfilename[200];
FILE *gdumpin;


int totalzones;
int defcoord;
FTYPE h_over_r,hslope;
FTYPE jetalpha;

// GODMARK3D -- should these be stored in coordparms.dat?
FTYPE Rin_array[NDIM], Rout_array[NDIM];  //atch -- arrays for a more general way of handling the grid dimensions

FTYPE Rhor,Rout,dx[NDIM],startx[NDIM],endx[NDIM],R0,Diffx[NDIM];
FTYPE dxdxp[NDIM][NDIM];
int myid;
int debugfail;
FTYPE dtdump;
FTYPE dt;

FTYPE spc_target[1+NDIM]; // newt quantity, so starts at index 1, not 0

FTYPE Lunit,dV,dVF;
//int N1,N2,N3;
//int N1BND, N2BND, N3BND;

int BCtype[COMPDIM*2]; // not important for oN0>1 or nN1>1


FTYPE tnrdegrees;


FTYPE a;
int failed;
long nstroke;


int nn;
FTYPE *fvec;

int nrerrorflag;
FTYPE Rchop;

/////////////////////////////
// memory stuff
//
// STILL DO need to use GLOBALPOINT on pointer names since that's why these pointers are here at all, because coord.c uses those names.
// Don't need to worry about shifting them
// Ok to convert N?M->NSTORE?, but not required since never access these
// Ok to also convert SHIFT? -> and N?BND ->
// If time component exists, assume gdump time dimension has no time-dependent values (e.g. no time-dependent metric for anything needing gdump stuff, like vectors -- scalars are always fine without gdump info).
FTYPE (*GLOBALPOINT(pglobal))[NSTORE2][NSTORE3][NPR];
int didstorepositiondata;
FTYPE (*GLOBALPOINT(dxdxpstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM][NDIM];
FTYPE (*GLOBALPOINT(idxdxpstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM][NDIM];

FTYPE (*GLOBALPOINT(Xstore))[NSTORE1+SHIFTSTORE1*3][NSTORE2+SHIFTSTORE2*3][NSTORE3+SHIFTSTORE3*3][NDIM];
FTYPE (*GLOBALPOINT(Vstore))[NSTORE1+SHIFTSTORE1*3][NSTORE2+SHIFTSTORE2*3][NSTORE3+SHIFTSTORE3*3][NDIM];

// Note these below memory things are not affected by global.storage.h since created instead of as global arrays
// So this code presumes [h][i][j][k] format always and is not optimized for arbitrary storage mapping
// That is, matrix() is always (0,1,2,3) and access is always [h][i][j][k] associated with t(h), r(i), theta(j), phi(k)
unsigned char ****oldimage,****oldimage0,****newimage;
FTYPE ****olddata,****olddata0;
FTYPE ****newdata;




// global variable for NUMREC function call
void (*nrfuncv)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]);

