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

int oN1,oN2,oN3,nN1,nN2,nN3 ;
FTYPE refinefactor;
int roN1,roN2,roN3;
FTYPE dX[NDIM],Rin,fakeRin,dxc,dyc,dzc,fakedxc,fakedyc,fakedzc;
FTYPE Zin,Zout,Zeqin;
FTYPE startxc, startyc, startzc, endxc, endyc, endzc;
int newgridtype,oldgridtype;
FTYPE gridAAglobal,gridr0global;

FTYPE X[NDIM];
FTYPE Xmetricnew[NDIM],Xmetricold[NDIM]; // used to store time of latest and oldest metric


FTYPE t,gam,spin,QBH,MBH;
int startpos[NDIM];
int totalsize[NDIM];
long realnstep,nstep;
FTYPE readnstep;

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
FTYPE dt;

FTYPE spc_target[NDIM];

FTYPE Lunit,dV,dVF;
//int N1,N2,N3;
//int N1BND, N2BND, N3BND;

int BCtype[COMPDIM*2];



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
FTYPE (*GLOBALPOINT(pglobal))[NSTORE2][NSTORE3][NPR];
int didstorepositiondata;
FTYPE (*GLOBALPOINT(dxdxpstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM][NDIM];
FTYPE (*GLOBALPOINT(idxdxpstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM][NDIM];
FTYPE (*GLOBALPOINT(Xstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM];
FTYPE (*GLOBALPOINT(Vstore))[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM];

// Note these below memory things are not affected by global.storage.h since created instead of as global arrays
// So this code presumes [i][j][k] format always and is not optimized for arbitrary storage mapping
// That is, matrix() is always (1,2,3) and access is always [i][j][k] associated with r(i), theta(j), phi(k)
unsigned char ***oldimage,***oldimage0,***newimage;
FTYPE ***olddata,***olddata0;
FTYPE ***newdata;




// global variable for NUMREC function call
void (*nrfuncv)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]);

