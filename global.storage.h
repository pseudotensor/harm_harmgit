
/*! \file global.storage.h
    \brief macros and definitions related to storage mapping functions

//////////////
//
// Macrofication of storage for multi-dimensional arrays:
//
// SUPERNOTE: To use old codes (e.g. inits) with this code, either convert them OR set ORDERSTORAGE==0 in definit.h and then the code is just like original in behavior even with macrofication.
//
// This file contains the macrofication of storage accesses and definitions for multi-dimensional arrays
// That is, definition a[N1M][N2M][N3M] with access a[i][j][k] with code association r(i) theta(j) phi(k)
// restricts L1,L2 cache lines to along k-direction.  Code assoation (e.g.) r(i) is well-defined
// through-out code and changes would entail changing user notions and code associations between dir==1 and i.
// Instead, macrofy only how storage arrays are
// defined and accessed so the user always puts in a(i,j,k), but the user can choose which index is fastest.

// Anything that has to do with storing or accessing specific directions of spatial memory in multi-dimensional storage arrays should be here and considered here.
// That is: a[?][?][?] accesses or allocation: FTYPE a[?][?][?] or pointer definitions for those arrays: b=&(a[?][?][?])
// 

// SUPEROPTMARK: Clearly the best way to optimize for multi-core is to do ALL calculations for a single cell at once without ever storing anything.  Then, one just starts with p at (say) +-MAXBND cells and fully determines new p.  There will be MANY repeated calculations, but when memory is a bottleneck that's the way to go.
// Then use OpenMP or pthreads to parallelize over each position.
// Then only memory grabbing is p in multi-D and other very expensive things like the metric that would likely be too expensive to recompute for each index for p.

// See kazfulleos.global.h for some non-spatial array macros for general EOS

*/



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///
/// Define N1M, N2M, N3M and other spatial storage related items
///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///  allocated memory uses this for active zones 0-N1-1 and bc beyond that 
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)


///  NBIGM is bigger of N1M and N2M and N3M 
/// not currently used
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

/// surface areas of sides WITH boundary zones
/// reduces to 0 if that dimension not used
#define MAXSURFA1 (N2M*N3M*N1NOT1)
#define MAXSURFA2 (N1M*N3M*N2NOT1)
#define MAXSURFA3 (N1M*N2M*N3NOT1)

/// maximal surface of boundary exchange
/// notice that this works in any number of dimensions and any N1,N2,N3
/// that is, it reduces correctly when a dimension is degenerate (N=1)
#define NBIGS1M ((MAXSURFA1>MAXSURFA2) ? MAXSURFA1 : MAXSURFA2)
#define NBIGSM ((NBIGS1M>MAXSURFA3) ? NBIGS1M : MAXSURFA3)


//////////////////////////////////
/////////////////////
///
///ORDERSTORAGE
///   : N1M=NSTORE?,  N2M=NSTORE?,  N3M=NSTORE?
/// 0 :           1             2             3 // standard (original) order
/// 1 :           2             3             1
/// 2 :           3             1             2
/// 3 :           1             3             2
/// 4 :           2             1             3
/// 5 :           3             2             1
///
/// ORDERSTORAGE==5 is most efficient order for grid extended in radius, less so in \theta, and even less so in \phi
////////#define ORDERSTORAGE 0 // set in definit.h or init.h
/// here, N1M,N2M,N3M = (r,theta,phi) always for SPC
/// (i,j,k) = (r,theta,\phi) always for SPC
/// But arrays accessed via array[IND1(i,j,k)][IND2(i,j,k)][IND3(i,j,k)]
/// ALL non-array references to i,j,k stay the same.  Only array accesses have changed.
/// Also change LOOP order (e.g. ORDERSTORAGE==5 : LOOPF -> LOOPF3 LOOPF2 LOOPF1)
/// This way (e.g.) coord(i,j,k) and all other such non-array accesses are the same, so little abstract code changes.
/// Only end up changing how arrays are accessed
/// Note that unlike Sasha's mixing of the 1,2,3 directions, this doesn't change (e.g.) r(i) \theta(j) and \phi(k) association.  It only changes how arrays are accessed!  So this doesn't test how different directions are handled in the code.
/// Note that we change the left-hand-side order not the right-hand-side order that should always be 1,2,3 or i,j,k


#define OPENMP3DLOOPBLOCK2IJKSTORAGE0(i,j,k)                    \
  k=nzshift+(int)((blockijk-1)%nzsize);                         \
  j=nyshift+(int)(((blockijk-1)%(nzsize*nysize))/nzsize);       \
  i=nxshift+(int)((blockijk-1)/(nzsize*nysize));

#define OPENMP3DLOOPBLOCK2IJKSTORAGE1(i,j,k)                    \
  i=nxshift+(int)((blockijk-1)%nxsize);                         \
  k=nzshift+(int)(((blockijk-1)%(nxsize*nzsize))/nxsize);       \
  j=nyshift+(int)((blockijk-1)/(nxsize*nzsize));

#define OPENMP3DLOOPBLOCK2IJKSTORAGE2(i,j,k)                    \
  j=nyshift+(int)((blockijk-1)%nysize);                         \
  i=nxshift+(int)(((blockijk-1)%(nysize*nxsize))/nysize);       \
  k=nzshift+(int)((blockijk-1)/(nysize*nxsize));

#define OPENMP3DLOOPBLOCK2IJKSTORAGE3(i,j,k)                    \
  j=nyshift+(int)((blockijk-1)%nysize);                         \
  k=nzshift+(int)(((blockijk-1)%(nysize*nzsize))/nysize);       \
  i=nxshift+(int)((blockijk-1)/(nysize*nzsize));

#define OPENMP3DLOOPBLOCK2IJKSTORAGE4(i,j,k)                    \
  k=nzshift+(int)((blockijk-1)%nzsize);                         \
  i=nxshift+(int)(((blockijk-1)%(nzsize*nxsize))/nzsize);       \
  j=nyshift+(int)((blockijk-1)/(nzsize*nxsize));

#define OPENMP3DLOOPBLOCK2IJKSTORAGE5(i,j,k)                    \
  i=nxshift+(int)((blockijk-1)%nxsize);                         \
  j=nyshift+(int)(((blockijk-1)%(nxsize*nysize))/nxsize);       \
  k=nzshift+(int)((blockijk-1)/(nxsize*nysize));


#if(ORDERSTORAGE==0)
/// 123
#define NSTORE1 N1M
#define NSTORE2 N2M
#define NSTORE3 N3M
#define NSTOREBND1 N1BND
#define NSTOREBND2 N2BND
#define NSTOREBND3 N3BND
#define SHIFTSTORE1 SHIFT1
#define SHIFTSTORE2 SHIFT2
#define SHIFTSTORE3 SHIFT3
#define LOOPORDER1(i,j,k) i
#define LOOPORDER2(i,j,k) j
#define LOOPORDER3(i,j,k) k

#define STO1(i,j,k) (i)
#define STO2(i,j,k) (j)
#define STO3(i,j,k) (k)
#define DEFDIM1(i,j,k) (i)
#define DEFDIM2(i,j,k) (j)
#define DEFDIM3(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET1(i,j,k) (i)
#define STOMET2(i,j,k) (j)
#define STOMET3(i,j,k) (k)
#define DEFDIMMET1(i,j,k) (i)
#define DEFDIMMET2(i,j,k) (j)
#define DEFDIMMET3(i,j,k) (k)
#else
#define STOMET1(i,j,k) (0)
#define STOMET2(i,j,k) (0)
#define STOMET3(i,j,k) (0)
#define DEFDIMMET1(i,j,k) (1)
#define DEFDIMMET2(i,j,k) (1)
#define DEFDIMMET3(i,j,k) (1)
#endif

#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE0(i,j,k)



#elif(ORDERSTORAGE==1)
///231
#define NSTORE2 N1M
#define NSTORE3 N2M
#define NSTORE1 N3M
#define NSTOREBND2 N1BND
#define NSTOREBND3 N2BND
#define NSTOREBND1 N3BND
#define SHIFTSTORE2 SHIFT1
#define SHIFTSTORE3 SHIFT2
#define SHIFTSTORE1 SHIFT3
#define LOOPORDER2(i,j,k) i
#define LOOPORDER3(i,j,k) j
#define LOOPORDER1(i,j,k) k

#define STO2(i,j,k) (i)
#define STO3(i,j,k) (j)
#define STO1(i,j,k) (k)
#define DEFDIM2(i,j,k) (i)
#define DEFDIM3(i,j,k) (j)
#define DEFDIM1(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET2(i,j,k) (i)
#define STOMET3(i,j,k) (j)
#define STOMET1(i,j,k) (k)
#define DEFDIMMET2(i,j,k) (i)
#define DEFDIMMET3(i,j,k) (j)
#define DEFDIMMET1(i,j,k) (k)
#else
#define STOMET2(i,j,k) (0)
#define STOMET3(i,j,k) (0)
#define STOMET1(i,j,k) (0)
#define DEFDIMMET2(i,j,k) (1)
#define DEFDIMMET3(i,j,k) (1)
#define DEFDIMMET1(i,j,k) (1)
#endif

#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE1(i,j,k)


#elif(ORDERSTORAGE==2)
///312
#define NSTORE3 N1M
#define NSTORE1 N2M
#define NSTORE2 N3M
#define NSTOREBND3 N1BND
#define NSTOREBND1 N2BND
#define NSTOREBND2 N3BND
#define SHIFTSTORE3 SHIFT1
#define SHIFTSTORE1 SHIFT2
#define SHIFTSTORE2 SHIFT3
#define LOOPORDER3(i,j,k) i
#define LOOPORDER1(i,j,k) j
#define LOOPORDER2(i,j,k) k

#define STO3(i,j,k) (i)
#define STO1(i,j,k) (j)
#define STO2(i,j,k) (k)
#define DEFDIM3(i,j,k) (i)
#define DEFDIM1(i,j,k) (j)
#define DEFDIM2(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET3(i,j,k) (i)
#define STOMET1(i,j,k) (j)
#define STOMET2(i,j,k) (k)
#define DEFDIMMET3(i,j,k) (i)
#define DEFDIMMET1(i,j,k) (j)
#define DEFDIMMET2(i,j,k) (k)
#else
#define STOMET3(i,j,k) (0)
#define STOMET1(i,j,k) (0)
#define STOMET2(i,j,k) (0)
#define DEFDIMMET3(i,j,k) (1)
#define DEFDIMMET1(i,j,k) (1)
#define DEFDIMMET2(i,j,k) (1)
#endif

#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE2(i,j,k)


#elif(ORDERSTORAGE==3)
///132
#define NSTORE1 N1M
#define NSTORE3 N2M
#define NSTORE2 N3M
#define NSTOREBND1 N1BND
#define NSTOREBND3 N2BND
#define NSTOREBND2 N3BND
#define SHIFTSTORE1 SHIFT1
#define SHIFTSTORE3 SHIFT2
#define SHIFTSTORE2 SHIFT3
#define LOOPORDER1(i,j,k) i
#define LOOPORDER3(i,j,k) j
#define LOOPORDER2(i,j,k) k

#define STO1(i,j,k) (i)
#define STO3(i,j,k) (j)
#define STO2(i,j,k) (k)
#define DEFDIM1(i,j,k) (i)
#define DEFDIM3(i,j,k) (j)
#define DEFDIM2(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET1(i,j,k) (i)
#define STOMET3(i,j,k) (j)
#define STOMET2(i,j,k) (k)
#define DEFDIMMET1(i,j,k) (i)
#define DEFDIMMET3(i,j,k) (j)
#define DEFDIMMET2(i,j,k) (k)
#else
#define STOMET1(i,j,k) (0)
#define STOMET3(i,j,k) (0)
#define STOMET2(i,j,k) (0)
#define DEFDIMMET1(i,j,k) (1)
#define DEFDIMMET3(i,j,k) (1)
#define DEFDIMMET2(i,j,k) (1)
#endif

#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE3(i,j,k)


#elif(ORDERSTORAGE==4)
///213
#define NSTORE2 N1M
#define NSTORE1 N2M
#define NSTORE3 N3M
#define NSTOREBND2 N1BND
#define NSTOREBND1 N2BND
#define NSTOREBND3 N3BND
#define SHIFTSTORE2 SHIFT1
#define SHIFTSTORE1 SHIFT2
#define SHIFTSTORE3 SHIFT3
#define LOOPORDER2(i,j,k) i
#define LOOPORDER1(i,j,k) j
#define LOOPORDER3(i,j,k) k

#define STO2(i,j,k) (i)
#define STO1(i,j,k) (j)
#define STO3(i,j,k) (k)
#define DEFDIM2(i,j,k) (i)
#define DEFDIM1(i,j,k) (j)
#define DEFDIM3(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET2(i,j,k) (i)
#define STOMET1(i,j,k) (j)
#define STOMET3(i,j,k) (k)
#define DEFDIMMET2(i,j,k) (i)
#define DEFDIMMET1(i,j,k) (j)
#define DEFDIMMET3(i,j,k) (k)
#else
#define STOMET2(i,j,k) (0)
#define STOMET1(i,j,k) (0)
#define STOMET3(i,j,k) (0)
#define DEFDIMMET2(i,j,k) (1)
#define DEFDIMMET1(i,j,k) (1)
#define DEFDIMMET3(i,j,k) (1)
#endif

#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE4(i,j,k)


#elif(ORDERSTORAGE==5)
///321
#define NSTORE3 N1M
#define NSTORE2 N2M
#define NSTORE1 N3M
#define NSTOREBND3 N1BND
#define NSTOREBND2 N2BND
#define NSTOREBND1 N3BND
#define SHIFTSTORE3 SHIFT1
#define SHIFTSTORE2 SHIFT2
#define SHIFTSTORE1 SHIFT3
#define LOOPORDER3(i,j,k) i
#define LOOPORDER2(i,j,k) j
#define LOOPORDER1(i,j,k) k

#define STO3(i,j,k) (i)
#define STO2(i,j,k) (j)
#define STO1(i,j,k) (k)
#define DEFDIM3(i,j,k) (i)
#define DEFDIM2(i,j,k) (j)
#define DEFDIM1(i,j,k) (k)

#if(MCOORD!=CARTMINKMETRIC)
#define STOMET3(i,j,k) (i)
#define STOMET2(i,j,k) (j)
#define STOMET1(i,j,k) (k)
#define DEFDIMMET3(i,j,k) (i)
#define DEFDIMMET2(i,j,k) (j)
#define DEFDIMMET1(i,j,k) (k)
#else
#define STOMET3(i,j,k) (0)
#define STOMET2(i,j,k) (0)
#define STOMET1(i,j,k) (0)
#define DEFDIMMET3(i,j,k) (1)
#define DEFDIMMET2(i,j,k) (1)
#define DEFDIMMET1(i,j,k) (1)
#endif


#define OPENMP3DLOOPBLOCK2IJK(i,j,k) OPENMP3DLOOPBLOCK2IJKSTORAGE5(i,j,k)


#endif

/// below doesn't work since #define can't refer to other #'s
///#define OPENMPLOOPGEN(is,ie,js,je,ks,ke,VARPRIVATE,EXTRAFOR) OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP(is,ie,js,je,ks,ke); \
///#pragma omp parallel private(VARPRIVATE) OPENMPGLOBALPRIVATEFULL \
///  {         \
///#pragma omp for schedule(OPENMPSCHEDULE,OPENMPCHUNKSIZE(blocksize)) EXTRAFOR  \
///    OPENMP3DLOOPBLOCK{      \
///      OPENMP3DLOPBLOCK2IJK(i,j,k);


/// Seems I caught everything, but not tested.
///#if(MCOORD==CARTMINKMETRIC)
///#error Not setup yet for CARTMINKMETRIC, need to define macros below for metric accesses so uses MET versions of above so correct access (i.e. i=j=k=0 for metric type stuff).  Then ensure those vars (global or not) are everywhere accessed and defined with such macros, including their pointers
///#endif

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///
/// Define Macros for dealing with global spatial pointers
///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


///See docs regarding # and ## use in C macros as processed by preprocessor: http://gcc.gnu.org/onlinedocs/cpp/Concatenation.html#Concatenation
/// Didn't see how to macrofy a_ and a_s_ so they aren't repeated, but not a big deal since still those defined only here
/// MACP?A? corresponds to P? number of prior array indices and A? to number of after indices

/// generic macros for a_ or a_s_ prefixes corresponding, respectively, to base or shifted pointers
#define GENPOINT(prefix,name) prefix##name
#define GENMAC(prefix,name,i,j,k) prefix##name[STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
#define GENMACP0A0(prefix,name,i,j,k) GENMAC(prefix,name,i,j,k)

#define GENMACP1A0(prefix,name,argp1,i,j,k) prefix##name[argp1][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
#define GENMACP0A1(prefix,name,i,j,k,arga1) prefix##name[STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1]

#define GENMACP2A0(prefix,name,argp1,argp2,i,j,k) prefix##name[argp1][argp2][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
#define GENMACP1A1(prefix,name,argp1,i,j,k,arga1) prefix##name[argp1][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1]
#define GENMACP0A2(prefix,name,i,j,k,arga1,arga2) prefix##name[STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2]

#define GENMACP3A0(prefix,name,argp1,argp2,argp3,i,j,k) prefix##name[argp1][argp2][argp3][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
#define GENMACP2A1(prefix,name,argp1,argp2,i,j,k,arga1) prefix##name[argp1][argp2][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1]
#define GENMACP1A2(prefix,name,argp1,i,j,k,arga1,arga2) prefix##name[argp1][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2]
#define GENMACP0A3(prefix,name,i,j,k,arga1,arga2,arga3) prefix##name[STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2][arga3]

#define GENMACP4A0(prefix,name,argp1,argp2,argp3,argp4,i,j,k) prefix##name[argp1][argp2][argp3][argp4][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
#define GENMACP3A1(prefix,name,argp1,argp2,argp3,i,j,k,arga1) prefix##name[argp1][argp2][argp3][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1]
#define GENMACP2A2(prefix,name,argp1,argp2,i,j,k,arga1,arga2) prefix##name[argp1][argp2][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2]
#define GENMACP1A3(prefix,name,argp1,i,j,k,arga1,arga2,arga3) prefix##name[argp1][STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2][arga3]
#define GENMACP0A4(prefix,name,i,j,k,arga1,arga2,arga3,arga4) prefix##name[STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)][arga1][arga2][arga3][arga4]

/////////////////
/// BASE (a_) wrappers // actual declarations for 3D memory arrays set in superdefs.h
/// 
#define BASEPOINT(name) GENPOINT(a_,name)
#define BASEMAC(name,i,j,k) GENMAC(a_,name,i,j,k)
#define BASEMACP0A0(name,i,j,k) GENMACP0A0(a_,name,i,j,k)

#define BASEMACP1A0(name,argp1,i,j,k) GENMACP1A0(a_,name,argp1,i,j,k)
#define BASEMACP0A1(name,i,j,k,arga1) GENMACP0A1(a_,name,i,j,k,arga1)

#define BASEMACP2A0(name,argp1,argp2,i,j,k) GENMACP2A0(a_,name,argp1,argp2,i,j,k)
#define BASEMACP1A1(name,argp1,i,j,k,arga1) GENMACP1A1(a_,name,argp1,i,j,k,arga1)
#define BASEMACP0A2(name,i,j,k,arga1,arga2) GENMACP0A2(a_,name,i,j,k,arga1,arga2)

#define BASEMACP3A0(name,argp1,argp2,argp3,i,j,k) GENMACP3A0(a_,name,argp1,argp2,argp3,i,j,k)
#define BASEMACP2A1(name,argp1,argp2,i,j,k,arga1) GENMACP2A1(a_,name,argp1,argp2,i,j,k,arga1)
#define BASEMACP1A2(name,argp1,i,j,k,arga1,arga2) GENMACP1A2(a_,name,argp1,i,j,k,arga1,arga2)
#define BASEMACP0A3(name,i,j,k,arga1,arga2,arga3) GENMACP0A3(a_,name,i,j,k,arga1,arga2,arga3)

#define BASEMACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) GENMACP4A0(a_,name,argp1,argp2,argp3,argp4,i,j,k)
#define BASEMACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) GENMACP3A1(a_,name,argp1,argp2,argp3,i,j,k,arga1)
#define BASEMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMACP2A2(a_,name,argp1,argp2,i,j,k,arga1,arga2)
#define BASEMACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) GENMACP1A3(a_,name,argp1,i,j,k,arga1,arga2,arga3)
#define BASEMACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) GENMACP0A4(a_,name,i,j,k,arga1,arga2,arga3,arga4)

/////////////////
/// shifted pointer (a_s_) wrappers
/// Couldn't figure out how to macrofy a_s_ or a_, so when changing have to change all manually
/// If want no header on global array, then just leave as (,name,...) instead of (a_s_,name,...)
///
#define GLOBALPOINT(name) GENPOINT(a_s_,name)
#define GLOBALMAC(name,i,j,k) GENMAC(a_s_,name,i,j,k)
#define GLOBALMACP0A0(name,i,j,k) GENMACP0A0(a_s_,name,i,j,k)

#define GLOBALMACP1A0(name,argp1,i,j,k) GENMACP1A0(a_s_,name,argp1,i,j,k)
#define GLOBALMACP0A1(name,i,j,k,arga1) GENMACP0A1(a_s_,name,i,j,k,arga1)

#define GLOBALMACP2A0(name,argp1,argp2,i,j,k) GENMACP2A0(a_s_,name,argp1,argp2,i,j,k)
#define GLOBALMACP1A1(name,argp1,i,j,k,arga1) GENMACP1A1(a_s_,name,argp1,i,j,k,arga1)
#define GLOBALMACP0A2(name,i,j,k,arga1,arga2) GENMACP0A2(a_s_,name,i,j,k,arga1,arga2)

#define GLOBALMACP3A0(name,argp1,argp2,argp3,i,j,k) GENMACP3A0(a_s_,name,argp1,argp2,argp3,i,j,k)
#define GLOBALMACP2A1(name,argp1,argp2,i,j,k,arga1) GENMACP2A1(a_s_,name,argp1,argp2,i,j,k,arga1)
#define GLOBALMACP1A2(name,argp1,i,j,k,arga1,arga2) GENMACP1A2(a_s_,name,argp1,i,j,k,arga1,arga2)
#define GLOBALMACP0A3(name,i,j,k,arga1,arga2,arga3) GENMACP0A3(a_s_,name,i,j,k,arga1,arga2,arga3)

#define GLOBALMACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) GENMACP4A0(a_s_,name,argp1,argp2,argp3,argp4,i,j,k)
#define GLOBALMACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) GENMACP3A1(a_s_,name,argp1,argp2,argp3,i,j,k,arga1)
#define GLOBALMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMACP2A2(a_s_,name,argp1,argp2,i,j,k,arga1,arga2)
#define GLOBALMACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) GENMACP1A3(a_s_,name,argp1,i,j,k,arga1,arga2,arga3)
#define GLOBALMACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) GENMACP0A4(a_s_,name,i,j,k,arga1,arga2,arga3,arga4)

/////////////////
/// Access to any non-global pointer (not a_s_):
/// note that ok to have nothing as prefix, but macro using prefix##name assumes prefix is defined even inputted as nothing.
#define POINT(name) GENPOINT(,name)
#define MAC(name,i,j,k) GENMAC(,name,i,j,k)
#define MACP0A0(name,i,j,k) GENMACP0A0(,name,i,j,k)

#define MACP1A0(name,argp1,i,j,k) GENMACP1A0(,name,argp1,i,j,k)
#define MACP0A1(name,i,j,k,arga1) GENMACP0A1(,name,i,j,k,arga1)

#define MACP2A0(name,argp1,argp2,i,j,k) GENMACP2A0(,name,argp1,argp2,i,j,k)
#define MACP1A1(name,argp1,i,j,k,arga1) GENMACP1A1(,name,argp1,i,j,k,arga1)
#define MACP0A2(name,i,j,k,arga1,arga2) GENMACP0A2(,name,i,j,k,arga1,arga2)

#define MACP3A0(name,argp1,argp2,argp3,i,j,k) GENMACP3A0(,name,argp1,argp2,argp3,i,j,k)
#define MACP2A1(name,argp1,argp2,i,j,k,arga1) GENMACP2A1(,name,argp1,argp2,i,j,k,arga1)
#define MACP1A2(name,argp1,i,j,k,arga1,arga2) GENMACP1A2(,name,argp1,i,j,k,arga1,arga2)
#define MACP0A3(name,i,j,k,arga1,arga2,arga3) GENMACP0A3(,name,i,j,k,arga1,arga2,arga3)

#define MACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) GENMACP4A0(,name,argp1,argp2,argp3,argp4,i,j,k)
#define MACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) GENMACP3A1(,name,argp1,argp2,argp3,i,j,k,arga1)
#define MACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMACP2A2(,name,argp1,argp2,i,j,k,arga1,arga2)
#define MACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) GENMACP1A3(,name,argp1,i,j,k,arga1,arga2,arga3)
#define MACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) GENMACP0A4(,name,i,j,k,arga1,arga2,arga3,arga4)

/////////////////
/// pointer shift header wrappers
/// Originally just did, e.g.,  GLOBALPOINT(pglobal) = (FTYPE (*)[N2M][N3M][NPR]) (&(BASEPOINT(pglobal)[NSTOREBND1][NSTOREBND2][NSTOREBND3][0]));
/// But this exposes the internal storage, so avoid, and easy to avoid
/// Below is same structurly as above, but with name and the first index removed from resulting macro expansion, since that's how pointers are referenced in (e.g.) set_arrays_multidimen.c
/// Force user to input name and all dimensions to structurly consistent with other related macros
/// PTRMAC used for both shifting pointers and defining pointers (with different indicies)
/// Note that we must use DEFDIM1,2,3 since different than STO in behavior when no CARTMINKMETRIC.  When [DEFDIM1(i,j,k)] does not appear, that could be any i,j,k so that's correct.
/// Note we remove (always) first term of array to create its pointer reference
#define PTRMAC(name,i,j,k) (*)[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRMACP0A0(name,i,j,k) PTRMAC(name,i,j,k)

#define PTRMACP1A0(name,argp1,i,j,k) (*)[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRMACP0A1(name,i,j,k,arga1) (*)[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]

#define PTRMACP2A0(name,argp1,argp2,i,j,k) (*)[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRMACP1A1(name,argp1,i,j,k,arga1) (*)[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRMACP0A2(name,i,j,k,arga1,arga2) (*)[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]

#define PTRMACP3A0(name,argp1,argp2,argp3,i,j,k) (*)[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRMACP2A1(name,argp1,argp2,i,j,k,arga1) (*)[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRMACP1A2(name,argp1,i,j,k,arga1,arga2) (*)[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRMACP0A3(name,i,j,k,arga1,arga2,arga3) (*)[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]

#define PTRMACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) (*)[argp2][argp3][argp4][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRMACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) (*)[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*)[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRMACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) (*)[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]
#define PTRMACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) (*)[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3][arga4]

/////////////
/// Define pointers to multi-dimen arrays for *global* variables (used in superdefs.pointers.h)
///
#define PTRDEFGLOBALMAC(name,i,j,k) (*GLOBALPOINT(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFGLOBALMACP0A0(name,i,j,k) PTRDEFGLOBALMAC(name,i,j,k)

#define PTRDEFGLOBALMACP1A0(name,argp1,i,j,k) (*GLOBALPOINT(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFGLOBALMACP0A1(name,i,j,k,arga1) (*GLOBALPOINT(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]

#define PTRDEFGLOBALMACP2A0(name,argp1,argp2,i,j,k) (*GLOBALPOINT(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFGLOBALMACP1A1(name,argp1,i,j,k,arga1) (*GLOBALPOINT(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFGLOBALMACP0A2(name,i,j,k,arga1,arga2) (*GLOBALPOINT(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]

#define PTRDEFGLOBALMACP3A0(name,argp1,argp2,argp3,i,j,k) (*GLOBALPOINT(name))[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFGLOBALMACP2A1(name,argp1,argp2,i,j,k,arga1) (*GLOBALPOINT(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFGLOBALMACP1A2(name,argp1,i,j,k,arga1,arga2) (*GLOBALPOINT(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRDEFGLOBALMACP0A3(name,i,j,k,arga1,arga2,arga3) (*GLOBALPOINT(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]

#define PTRDEFGLOBALMACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) (*GLOBALPOINT(name))[argp2][argp3][argp4][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFGLOBALMACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) (*GLOBALPOINT(name))[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFGLOBALMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*GLOBALPOINT(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRDEFGLOBALMACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) (*GLOBALPOINT(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]
#define PTRDEFGLOBALMACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) (*GLOBALPOINT(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3][arga4]

/////////////
/// Define pointers to multi-dimen arrays for non-global arrays (used in all other places except superdefs.pointers.h)
///
#define PURENAME(name) GENPOINT(,name)
#define PTRDEFMAC(name,i,j,k) (*PURENAME(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFMACP0A0(name,i,j,k) PTRDEFMAC(name,i,j,k)

#define PTRDEFMACP1A0(name,argp1,i,j,k) (*PURENAME(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFMACP0A1(name,i,j,k,arga1) (*PURENAME(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]

#define PTRDEFMACP2A0(name,argp1,argp2,i,j,k) (*PURENAME(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFMACP1A1(name,argp1,i,j,k,arga1) (*PURENAME(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFMACP0A2(name,i,j,k,arga1,arga2) (*PURENAME(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]

#define PTRDEFMACP3A0(name,argp1,argp2,argp3,i,j,k) (*PURENAME(name))[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFMACP2A1(name,argp1,argp2,i,j,k,arga1) (*PURENAME(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFMACP1A2(name,argp1,i,j,k,arga1,arga2) (*PURENAME(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRDEFMACP0A3(name,i,j,k,arga1,arga2,arga3) (*PURENAME(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]

#define PTRDEFMACP4A0(name,argp1,argp2,argp3,argp4,i,j,k) (*PURENAME(name))[argp2][argp3][argp4][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)]
#define PTRDEFMACP3A1(name,argp1,argp2,argp3,i,j,k,arga1) (*PURENAME(name))[argp2][argp3][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1]
#define PTRDEFMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*PURENAME(name))[argp2][DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2]
#define PTRDEFMACP1A3(name,argp1,i,j,k,arga1,arga2,arga3) (*PURENAME(name))[DEFDIM1(i,j,k)][DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3]
#define PTRDEFMACP0A4(name,i,j,k,arga1,arga2,arga3,arga4) (*PURENAME(name))[DEFDIM2(i,j,k)][DEFDIM3(i,j,k)][arga1][arga2][arga3][arga4]





/// METRIC VERSIONS (only written down if required):
#define GENMETMAC(prefix,name,i,j,k) prefix##name[STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)]
#define BASEMETMAC(name,i,j,k) GENMETMAC(a_,name,i,j,k)
#define GLOBALMETMAC(name,i,j,k) GENMETMAC(a_s_,name,i,j,k)
#define METMAC(name,i,j,k) GENMETMAC(,name,i,j,k)
#define PTRMETMAC(name,i,j,k) (*)[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFGLOBALMETMAC(name,i,j,k) (*GLOBALPOINT(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFMETMAC(name,i,j,k) (*PURENAME(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]

#define GENMETMACP0A1(prefix,name,i,j,k,arga1) prefix##name[STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)][arga1]
#define BASEMETMACP0A1(name,i,j,k,arga1) GENMETMACP0A1(a_,name,i,j,k,arga1)
#define GLOBALMETMACP0A1(name,i,j,k,arga1) GENMETMACP0A1(a_s_,name,i,j,k,arga1)
#define METMACP0A1(name,i,j,k,arga1) GENMETMACP0A1(,name,i,j,k,arga1)
#define PTRMETMACP0A1(name,i,j,k,arga1) (*)[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]
#define PTRDEFGLOBALMETMACP0A1(name,i,j,k,arga1) (*GLOBALPOINT(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]
#define PTRDEFMETMACP0A1(name,i,j,k,arga1) (*PURENAME(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]

#define GENMETMACP1A0(prefix,name,argp1,i,j,k) prefix##name[argp1][STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)]
#define BASEMETMACP1A0(name,argp1,i,j,k) GENMETMACP1A0(a_,name,argp1,i,j,k)
#define GLOBALMETMACP1A0(name,argp1,i,j,k) GENMETMACP1A0(a_s_,name,argp1,i,j,k)
#define METMACP1A0(name,argp1,i,j,k) GENMETMACP1A0(,name,argp1,i,j,k)
#define PTRMETMACP1A0(name,argp1,i,j,k) (*)[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFGLOBALMETMACP1A0(name,argp1,i,j,k) (*GLOBALPOINT(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFMETMACP1A0(name,argp1,i,j,k) (*PURENAME(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]

#define GENMETMACP1A1(prefix,name,argp1,i,j,k,arga1) prefix##name[argp1][STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)][arga1]
#define BASEMETMACP1A1(name,argp1,i,j,k,arga1) GENMETMACP1A1(a_,name,argp1,i,j,k,arga1)
#define GLOBALMETMACP1A1(name,argp1,i,j,k,arga1) GENMETMACP1A1(a_s_,name,argp1,i,j,k,arga1)
#define METMACP1A1(name,argp1,i,j,k,arga1) GENMETMACP1A1(,name,argp1,i,j,k,arga1)
#define PTRMETMACP1A1(name,argp1,i,j,k,arga1) (*)[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]
#define PTRDEFGLOBALMETMACP1A1(name,argp1,i,j,k,arga1) (*GLOBALPOINT(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]
#define PTRDEFMETMACP1A1(name,argp1,i,j,k,arga1) (*PURENAME(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1]

#define GENMETMACP1A2(prefix,name,argp1,i,j,k,arga1,arga2) prefix##name[argp1][STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)][arga1][arga2]
#define BASEMETMACP1A2(name,argp1,i,j,k,arga1,arga2) GENMETMACP1A2(a_,name,argp1,i,j,k,arga1,arga2)
#define GLOBALMETMACP1A2(name,argp1,i,j,k,arga1,arga2) GENMETMACP1A2(a_s_,name,argp1,i,j,k,arga1,arga2)
#define METMACP1A2(name,argp1,i,j,k,arga1,arga2) GENMETMACP1A2(,name,argp1,i,j,k,arga1,arga2)
#define PTRMETMACP1A2(name,argp1,i,j,k,arga1,arga2) (*)[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]
#define PTRDEFGLOBALMETMACP1A2(name,argp1,i,j,k,arga1,arga2) (*GLOBALPOINT(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]
#define PTRDEFMETMACP1A2(name,argp1,i,j,k,arga1,arga2) (*PURENAME(name))[DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]

/// P2A0
#define GENMETMACP2A0(prefix,name,argp1,argp2,i,j,k) prefix##name[argp1][argp2][STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)]
#define BASEMETMACP2A0(name,argp1,argp2,i,j,k) GENMETMACP2A0(a_,name,argp1,argp2,i,j,k)
#define GLOBALMETMACP2A0(name,argp1,argp2,i,j,k) GENMETMACP2A0(a_s_,name,argp1,argp2,i,j,k)
#define METMACP2A0(name,argp1,argp2,i,j,k) GENMETMACP2A0(,name,argp1,argp2,i,j,k)
#define PTRMETMACP2A0(name,argp1,argp2,i,j,k) (*)[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFGLOBALMETMACP2A0(name,argp1,argp2,i,j,k) (*GLOBALPOINT(name))[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]
#define PTRDEFMETMACP2A0(name,argp1,argp2,i,j,k) (*PURENAME(name))[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)]

/// P2A2
#define GENMETMACP2A2(prefix,name,argp1,argp2,i,j,k,arga1,arga2) prefix##name[argp1][argp2][STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)][arga1][arga2]
#define BASEMETMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMETMACP2A2(a_,name,argp1,argp2,i,j,k,arga1,arga2)
#define GLOBALMETMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMETMACP2A2(a_s_,name,argp1,argp2,i,j,k,arga1,arga2)
#define METMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) GENMETMACP2A2(,name,argp1,argp2,i,j,k,arga1,arga2)
#define PTRMETMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*)[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]
#define PTRDEFGLOBALMETMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*GLOBALPOINT(name))[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]
#define PTRDEFMETMACP2A2(name,argp1,argp2,i,j,k,arga1,arga2) (*PURENAME(name))[argp2][DEFDIMMET1(i,j,k)][DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2]

#define GENMETMACP0A3(prefix,name,i,j,k,arga1,arga2,arga3) prefix##name[STOMET1(i,j,k)][STOMET2(i,j,k)][STOMET3(i,j,k)][arga1][arga2][arga3]
#define BASEMETMACP0A3(name,i,j,k,arga1,arga2,arga3) GENMETMACP0A3(a_,name,i,j,k,arga1,arga2,arga3)
#define GLOBALMETMACP0A3(name,i,j,k,arga1,arga2,arga3) GENMETMACP0A3(a_s_,name,i,j,k,arga1,arga2,arga3)
#define METMACP0A3(name,i,j,k,arga1,arga2,arga3) GENMETMACP0A3(,name,i,j,k,arga1,arga2,arga3)
#define PTRMETMACP0A3(name,i,j,k,arga1,arga2,arga3) (*)[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2][arga3]
#define PTRDEFGLOBALMETMACP0A3(name,i,j,k,arga1,arga2,arga3) (*GLOBALPOINT(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2][arga3]
#define PTRDEFMETMACP0A3(name,i,j,k,arga1,arga2,arga3) (*PURENAME(name))[DEFDIMMET2(i,j,k)][DEFDIMMET3(i,j,k)][arga1][arga2][arga3]

/// 1) Check with (for each metric variable from (e.g.) superdefs.pointers.h):
/// make superclean ; grep "MAC" *.c *.h | grep "compgeom" | grep -v "METMAC" 
/// Should result in no results for any file except for commented lines or other things on same line

/// 2) Ensure all METMAC calls are either global or not as should be or not.
/// make superclean ; grep "METMAC" *.c *.h

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///
/// Instructions on how converted code to new form
///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// Now reference any access to an array using: MAC(simplename,i,j,k) and reference a pointer (e.g. for passing via functions) as GLOBALPOINT(simplename)
// Here simplename is the name originally used for the array, that has been replaced by a_s_<originalname>



// emacs regexp: http://jamesthornton.com/emacs/node/emacs_97.html
// In set_arrays_multidimen.c:

// NOTE1) recall that can't put non-functional macro as argument to functional macro since it won't expand.  For example, can't convert just: array[ip1] -> MAC(array,ip1) since ip1 won't expand into (e.g.) i+1 if N1>1.  Must convert to MAC(array,ip1mac(i))
// NOTE2) Note all y[?][?][?] will be multi-dim array of the sort here, and of course some have more indices in front or after spatial indices

// What was done:
//
// 1) I first replaced: a_\([_a-zA-Z0-9]+\)\[ ->  BASEPOINT(\1)[
// 1.5) In set_arrays.c: I replaced (checking each one) all \([_a-zA-Z0-9]+\) *= *( -> POINT(\1) = (
// 2) In set_arrays.c: For pointer assignments, I then replaced all <simplename> -> POINT(simplename)
// 3) In set_arrays.c: For value assignments, I replaced simplename[...][...][i][j][k][...][...] with one of the MAC(name,i,j,k) macros depending upon if additional indices before and/or after spatial indices
//    Emacs regexp replace is (start with longest versions first):

//    a) \([_a-zA-Z0-9]+\)\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\] -> MACP4A0(\1,\2,\3,\4,\5,\6,\7,\8) [pvcorninterp only] MACP0A4(\1,\2,\3,\4,\5,\6,\7,\8) [enodebugarray only] [and other forms]

//    b) \([_a-zA-Z0-9]+\)\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\] -> MACP3A0(\1,\2,\3,\4,\5,\6,\7) [pbcorninterp], MACP1A2(\1,\2,\3,\4,\5,\6,\7) [gcon,gcov,gcovlast,dxdxpstore,idxdxpstore] , MACP0A3(\1,\2,\3,\4,\5,\6,\7) [conn]  [and other forms]

//    c) \([_a-zA-Z0-9]+\)\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\] -> MACP2A0(\1,\2,\3,\4,\5,\6) [fluxstate,wspeed,wspeedcorn] or MACP1A1(\1,\2,\3,\4,\5,\6) [pk,gp_l,gp_r,gcovpert,eomfunc,beta,gcovpertlast,Xstore,Vstore]  or MACP0A2(\1,\2,\3,\4,\5,\6) [failfloorcount,cfaraday]

//    d) \([_a-zA-Z0-9]+\)\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\] -> MACP1A0(\1,\2,\3,\4,\5) [pother,emf,vpotarrayglobal,vpotanalytic,geomcornglobal,weno_prim_lower_order_fraction,EOSglobal,gdet,gdetvol,alphalapse,betasqoalphasq,alphalapselast,compgeom,compgeomlast] or MACP0A1(\1,\2,\3,\4,\5) [pglobal,panalytic,pstaganalytic,ptemparray,utemparray,vconemf,wspeedtemp,uinitglobal,ulastglobal,unewglobal,dUgeomarray,upointglobal,F1,F2,F3,F1EM,F2EM,F3EM,pleft,pright,prc,Bhatglobal,Bhatanalytic,pstagglobal,dq1,dq2,dq3,fluxvectemp,Fa,Fb,stencilvartemp,weno_lower_order_fraction,fluxdump,pflag,dissfunpos,fcon,jcon,<averagethings>,conn2,idxvol]

//    e) \([_a-zA-Z0-9]+\)\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\]\[\([a-zA-Z0-9-+]+\)\] -> MACP0A0(\1,\2,\3,\4) [fluxstatecent]

///////////////////
// For PTRMAC:
//
// 1) For set_arrays_multidimen.c (only need to do once per a,b,c,d,e since assume already have done above that gives macro name:
//    a) (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\), -> (*)PTR\7(\8,FILL,\1,\2,\3,\4,\5,\6)) (&(\7(\8,

//    b) (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\), -> (*)PTR\6(\7,FILL,\1,\2,\3,\4,\5)) (&(\6(\7,

//    c) (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\), -> (*)PTR\5(\6,FILL,\1,\2,\3,\4)) (&(\5(\6,

//    d) (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\), -> (*)PTR\4(\5,FILL,\1,\2,\3)) (&(\4(\5,

//    e) (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\), -> (*)PTR\3(\4,FILL,\1,\2)) (&(\3(\4,

//    f) THEN replace PTRBASE -> PTR

//    g) MUST replace FILL with correct thing (at least true for spatial indices).  This is so reordering knows what to use.
//       Just search/replace: FILL,N2M+SHIFT2 -> N1M+SHIFT1,N2M+SHIFT2  (and then after!)   FILL,N2M -> N1M,N2M
//       That gets all the required FILL's

//    h) Removed setting of pointer to zero since should be zero by default and requires by-hand work I don't want to do.



//////////////////////
// For PTRDEFMAC:

// For superdefs.pointers.h (after partially converted using POINT)

// 1) (*name)[N2M][N3M][NPR] -> PTRDEFMAC(name,i,j,k,pl)
//
//     a) ( *\* *POINT(\([_a-zA-Z0-9+-]+\) *) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *; -> PTRDEFMACP4A0(\1,FILL,\2,\3,\4,\5,\6,\7);   and  PTRDEFMACP0A4(\1,FILL,\2,\3,\4,\5,\6,\7);

//     b) ( *\* *POINT(\([_a-zA-Z0-9+-]+\) *) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *; -> PTRDEFMACP3A0(\1,FILL,\2,\3,\4,\5,\6);   and  PTRDEFMACP1A2(\1,FILL,\2,\3,\4,\5,\6);  and  PTRDEFMACP0A3(\1,FILL,\2,\3,\4,\5,\6);

//     c) ( *\* *POINT(\([_a-zA-Z0-9+-]+\) *) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *; -> PTRDEFMACP2A0(\1,FILL,\2,\3,\4,\5);   and  PTRDEFMACP1A1(\1,FILL,\2,\3,\4,\5);  and  PTRDEFMACP0A2(\1,FILL,\2,\3,\4,\5);

//     d) ( *\* *POINT(\([_a-zA-Z0-9+-]+\) *) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *; -> PTRDEFMACP1A0(\1,FILL,\2,\3,\4);   and  PTRDEFMACP0A1(\1,FILL,\2,\3,\4);

//     e) ( *\* *POINT(\([_a-zA-Z0-9+-]+\) *) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *; -> PTRDEFMACP0A0(\1,FILL,\2,\3);

//     f) Then correct   FILL,N2M+SHIFT2 -> N1M+SHIFT1,N2M+SHIFT2   and then after do   FILL,N2M -> N1M,N2M
//       Remaining FILL's are ok if not associated with spatial indices








//////////////////////
// For rest of files:


// 1) replace all global pointers with GLOBALMACP?A?() or GLOBALPOINT(name) as required.  Ensure not to make global something with same name as global array (e.g. F1,F2,F3).  Code will compile, but will not be correct or general.  Could make global arrays names different so forced to change.

// a) name[loc][i][j][k] -> GLOBALMACP1A0(name,loc,i,j,k):

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     ->    GLOBALMACP0A4(\1,\2,\3,\4,\5,\6,\7,\8) [enodebug]
//           GLOBALMACP4A0(\1,\2,\3,\4,\5,\6,\7,\8) [pvcorn]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> GLOBALMACP3A0(\1,\2,\3,\4,\5,\6,\7)
//        GLOBALMACP2A1(\1,\2,\3,\4,\5,\6,\7)
//        GLOBALMACP1A2(\1,\2,\3,\4,\5,\6,\7) [gcov gcovlast gcon]
//        GLOBALMACP0A3(\1,\2,\3,\4,\5,\6,\7) [full conn access]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> GLOBALMACP2A0(\1,\2,\3,\4,\5,\6)
//        GLOBALMACP1A1(\1,\2,\3,\4,\5,\6)
//        GLOBALMACP0A2(\1,\2,\3,\4,\5,\6) 

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> GLOBALMACP1A0(\1,\2,\3,\4,\5)
//        GLOBALMACP0A1(\1,\2,\3,\4,\5) [very common]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> GLOBALMAC(\1,\2,\3,\4)


// Careful with (e.g.) ip1 since needs to be macro version ip1mac(i) once inside macro.
// Note, I just removed (commented-out) such types of accesses and replaced with macro versions.
// Use macros if want to do (e.g.) ip1 or im1 in future.
// Don't change fun(a) with "a" as argument to function unless global array with global name



// 2) FTYPE (*a)[N2M][N3M] -> FTYPE (*a)[NSTORE2][NSTORE3] , etc.

// Instead of replacing normal pointers (non-global) in function arguments, declarations, etc. with PTRDEFMAC(),
// much simpler and ok for end-user to just have N1M->NSTORE1 N2M->NSTORE2 N3M->NSTORE3, etc.
// Do NOT replace function arguments that have arrays with PTRDEFMACP?A?() similar to above but without the POINT()
//

// N?M macro definitions are here since these are the only safe uses of N?M, except for exceptions listed now.
// All other uses are presumed to actually mean NSTORE?, like in function args and pointer definitions in functions.
// Files excempt from N?M changes to NSTORE?:
//
// ?.h:
// A1) global.storage.h
// A2) global.grmhd.h (just text discussions)
// A3) superdefs.h (because already converted to BASEMACP?A?() form instead of BASEPOINT() form
// A4) superdefs.pointers.h (because already converted to PTRDEFGLOBALMACP?A?() form instead of GLOBALPOINT and using NSTORE
// A5) reconstructeno.superdefs.h (as with superdefs.pointers.h and superdefs.h)
// A6) kazfulleos.superdefs.h ("")
//
// ?.c :
// A7) set_arrays_multidimen.c ("")
// A8) kazfulleos_set_arrays.c ("")
// A9) liaison_set_arrays.c ("")
// A10) reconstructeno_set_arrays.c ("")
//
// Other notes:
// A11) Ensure that reverted rest of files to not use PTRDEF or PTRMAC form, like in global.funcdeclare.h and *.c [true right now]
// A12) Do search-replace across multiple files: http://xahlee.org/emacs/find_replace_inter.html

// http://www.gnu.org/software/emacs/manual/html_node/emacs/Tags.html#Tags
// Recall M-X is ALT-X  and M-? is CTRL-? or ALT-?
//
// to use emacs to search/replace multiple files do:
// B0) create tag list exclusion within file "tagexcludelist.txt" with contents (without //'s):
// jon_interp.c
// jon_interp_computepreprocess.c
// jon_interp_filter.c
// jon_interp_interpolationitself.c
// bin2txt.c
// smcalc.c
// global.h
// defs.h
// defs.general.h
// mpidefs.h
// mpi_set_arrays.c
// nrutil2.c
// tensor.c
// supermpidefs.h
// global.storage.h
// global.grmhd.h
// superdefs.h
// superdecs.h
// superdefs.pointers.h
// superdecs.pointers.h
// kazfulleos.superdefs.h
// set_arrays_multidimen.c
// kazfulleos_set_arrays.c
// liaison_set_arrays.c
// superdefs.liaison.h
// reconstructeno_set_arrays.c
// reconstructeno.superdefs.h
//
// see bottom for final list and comments about each file
// 
// note: make superclean used to remove generated files.
// B1) make superclean ; ctags -e --exclude=@tagexcludelist.txt *.c *.h
// B2) emacs &
// B3) M-x visit-tags-table [hit enter]
// B4) M-x tags-query-replace
// B5) ADD to exclusion list if satisfied arrays present are completely converted as necessary to avoid cycling through non-convertible arrays repeatedly.
//     (Should be easy to process per-point files as indicated by maketail.harm.inc)
//     After some fixes, added other files as at bottom of this file
//     boundmpi.c
//     boundmpiint.c
//     mpi_init.c
//     kazfulleos.c
//
// a) Robust replace: \[N1M\]\[N2M\]\[N3M\] -> [NSTORE1][NSTORE2][NSTORE3]
// b) Robust replace: \[N2M\]\[N3M\] -> [NSTORE2][NSTORE3]  : since for pointers first [N1M] is hidden
// c) Robust replace: \[N1BND\]\[N2BND\]\[N3BND\] -> [NSTOREBND1][NSTOREBND2][NSTOREBND3]
// d) Robust replace: \[N2BND\]\[N3BND\] -> [NSTOREBND2][NSTOREBND3] : since for pointers first [N1M] is hidden [None of these since PTR'ified them]
// e) Robust replace: \[N1M\+SHIFT1\]\[N2M\+SHIFT2\]\[N3M\+SHIFT3\] -> [NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]
// f) Robust replace: \[N2M\+SHIFT2\]\[N3M\+SHIFT3\] -> [NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]  : since for pointers first [N1M] is hidden
// g) [Note: bounds.ns.c's aphifromoutflow shows example where user should use macros for memory, since each dimension not treated similarly]
//    [Note: [-N1BND] is  used for 1-D arrays, and should keep since really refers to radius dimension]
//
// Once done  here, then should no longer have any N1M,N2M,N3M in any files since only used for creating or pointing to storage arrays.
// Still should be use of SHIFT1,SHIFT2,SHIFT3,N1BND,N2BND,N3BND in files if associated with i,j,k and NOT storage arrays.
//
//     Check with things like: grep "\[" *.c *.h | grep N2M   with SHIFT and BND versions too : To ensure didn't miss storage definition (storage accesses inside []'s will be fixed during #3 below [Looks good right now]
//
// B5) C-x s ! [Save all buffers]





// 3) a[i][j][k]: Replace remaining pointers in code that access spatial-related memory with MACP?A?()

//    ANY array[?1][?2][?3] -> macarray(?1,?2,?3) in code.  This is easier and less bulky within normal code than changing all (e.g.) [i][j][k] -> [STO1(i,j,k)][STO2(i,j,k)][STO3(i,j,k)]
//    So every multi-D array should have a macarray() associated with it.  Inside macarray we will use STO?(i,j,k).

// Do same as #1 above for "rest of files," But instead of using GLOBALMACP?A?() use just MACP?A?(), and start from largest and go to smallest access sizes.

// Identify pointer by how it accesses.  Multi-dimen arrays should be accessed via: [i][j][k] in some form.  Don't get confused by 3D accesses that aren't spatially related!  Don't ignore when just (e.g.) [0][j][k] or [j][k] when still spatial.

// Careful with (e.g.) ip1 since needs to be macro version ip1mac(i) once inside macro [already converted by removing definitions and forcing user to always use macro versions]

// Need to track with all pointers that get renamed since while don't have to wrap-up the new name, have to wrap-up the access.

// Watch out for pointer definitions that have p[NSTORE1][...] since these should not be converted to MAC if already have (e.g.) NSTORE1 instead of N1M
// e.g.:
// grep "MAC" *.c *.h | grep STORE
//
// Watch out for p[][...] as well.
// Below replacements miss p[], which is ok, but could catch p[ ].  First do:
//
//     \([_a-zA-Z0-9]+\)\[ *\]\[ -> (*\1)[

// Then do (note below does not catch on [] but does catch on [ ], so why above done so avoids catching on [ ] (even though apparently there were no such instances, there could have been and it would have had MAC(name,NSTORE1,...) and that would be invalid).

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]
//     ->    MACP0A4(\1,\2,\3,\4,\5,\6,\7,\8) [enodebug]
//           MACP4A0(\1,\2,\3,\4,\5,\6,\7,\8) [pvcorn]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]
//     -> MACP3A0(\1,\2,\3,\4,\5,\6,\7) [pbcorn] [ignore EOS tables in kazfulleos.c]
//        MACP2A1(\1,\2,\3,\4,\5,\6,\7)
//        MACP1A2(\1,\2,\3,\4,\5,\6,\7) [gcov gcovlast gcon]
//        MACP0A3(\1,\2,\3,\4,\5,\6,\7) [full conn access]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]
//     -> MACP2A0(\1,\2,\3,\4,\5,\6) [wspeed fluxstate] -> MACP2A0(wspeed,\1,\2,\3,\4,\5)
//        MACP1A1(\1,\2,\3,\4,\5,\6) [fluxvec fluxvecEM fluxveca fluxvecb pl_ct pr_ct dqvec primface_l primface_r]  [ignore dirloopset]
//        MACP0A2(\1,\2,\3,\4,\5,\6) [gcovpert gcovpertlast failfloorcount]
//        [ignore eos interpolation stuff and primfactor]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]
//     -> MACP1A0(\1,\2,\3,\4,\5) [gdet pother (as global) emf A vpot]
//        MACP0A1(\1,\2,\3,\4,\5) [prim ui ucum p2interp panalytic(as global)] [very common]

//    \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\*\/\ ()]+\)\]
//     -> MAC(\1,\2,\3,\4) [fluxstatecent or passing any MACP0A1 type variable to a function]

//   Worry about user storage mapping functions, like in mpi_fileio.c for mapvaluejonio, that presume N3 is fastest, N2 slower, and N1 slowest or the inverse or whatever
//   mapvaluejonio is ok
//   Even ROMIO uses this mapping, so ROMIO good.

//   Finally check that no conversion led to MAC and NSTORE or NSTOREBND on same line:
//   grep MAC *.c *.h | grep NSTORE
//   [appears good right now]

//   Also, check that no remaining types of accesses (grep doesn't have +)
//   grep -e '[_a-zA-Z0-9]\{1,\}\[[_\>a-zA-Z0-9+-\*\/\ ()]*\]\[[_\>a-zA-Z0-9+-\*\/\ ()]*\]\[[_\>a-zA-Z0-9+-\*\/\ ()]*\]' *.c *.h > badness.temp.txt ; grep  --invert-match -e 'dirloopset\[' badness.temp.txt > badness.temp2.txt ; grep  --invert-match -e 'dirgenset\[' badness.temp2.txt > badness.temp3.txt ; grep  --invert-match -e 'workbc\[' badness.temp3.txt > badness.temp4.txt ; grep  --invert-match -e 'inoutlohi\[' badness.temp4.txt > badness.temp5.txt ; grep  --invert-match -e 'localpdotterms\[' badness.temp5.txt > badness.temp6.txt ; grep  --invert-match -e 'ijkcorn\[' badness.temp6.txt > badness.temp7.txt ; grep  --invert-match -e 'primfactor\[' badness.temp7.txt > badness.temp8.txt ; grep  --invert-match -e 'monoindicator\[' badness.temp8.txt > badness.temp9.txt ; grep  --invert-match -e 'pdottermsjet2\[' badness.temp9.txt > badness.temp10.txt ; grep  --invert-match -e 'pdottermsjet2_tot\[' badness.temp10.txt > badness.temp11.txt ; less badness.temp11.txt

// ignore jon_interp stuff, smcalc, or bin2txt stuff.

// Note that sometimes emacs with TAGS wouldn't search/replace in all files for some reason -- or TAGS was missing files.

// Ensure MAC appears to be used in right way:
// grep "MAC(" *.h *.c



// 4) Loops (at least COMP loops in global.comploops.h) need to be reordered depending upon ORDERSTORAGE
//    Change ALL spatial (i,j,k) type loops (2D loops too sometimes) to reorder depending upon ORDERSTORAGE
//    Note, code will work without this change since (e.g.) N1 is always associated with i.
//    But this is the whole point of the macrofication, to allow fastest index to be elongated in space on a CPU so avoid cache-misses.
//    Changes occur in global.loops.h [in some included files] and global.comploops.h only since collected all multi-D loops there including the boundary condition loops.
//    Note some loop parts in within code, and should replace them. (e.g. LOOPF3 LOOPF2 LOOPF1 within code itself not used as macro LOOPF)
//
//    Used LOOPORDER to control how each dimension was ordered in loops.  Nice, but didn't figure out how to macrofy each triple block separated by commas.


// 5) Things not to change:
//
// A) array[i][j] or array[j][k] if related to non-positional arrays like metric
// B) 1darray[i] if related to radius, like dissipatin or luminosity profiles [so only convert mutli-D arrays, not "radial" like 1D arrays

// To help ensure correct full conversion, renamed all multi-D arrays to a_<originalname> and pointers to a_s_<originalname> using macros
// So all original multi-D arrays will be undefined unless converted.


// 6) Compare with previous code for non macrofied parts:
// grep --invert-match -e "MAC" diffglobal.txt | grep --invert-match -e "STORE" | grep -e "^\!" | less



// Final tagexcludelist.txt after working on things:
//
//jon_interp.c
//jon_interp_computepreprocess.c
//jon_interp_filter.c
//jon_interp_interpolationitself.c
//bin2txt.c
//smcalc.c
//global.h
//defs.h
//defs.general.h
//mpidefs.h
//mpi_set_arrays.c
//supermpidefs.h
//nrutil2.c
//tensor.c
//global.storage.h
//global.grmhd.h
//superdefs.h
//superdecs.h
//superdefs.pointers.h
//superdecs.pointers.h
//kazfulleos.superdefs.h
//set_arrays_multidimen.c
//kazfulleos_set_arrays.c
//liaison_set_arrays.c
//superdefs.liaison.h
//reconstructeno_set_arrays.c
//reconstructeno.superdefs.h
//boundmpi.c
//boundmpiint.c
//mpi_init.c
//kazfulleos.c
//bounds.c
//bounds.ff.c
//bounds.fishmon.c
//boundsflux.c
//bounds.grb.c
//boundsint.c
//bounds.ns.c
//bounds.nsold.c
//bounds.ns.backup.c
//bounds.ns.backup2.c
//bounds.rebecca.c
//bounds.sasha.c
//bounds.tools.c
//diag.c
//dump_ener.c
//fixup.c
//fluxctstag.c
//flux.mergedc2ea2cmethod.c
//flux.c
//initbase.boundloop.c
//initbase.enerregions.c
//init.sasha.c
//init.ns.c
//init.ns.backup.c
//init.ns.backup2.c
//init.grb.c
//init.ff.c
//initbase.gridsectioning.c
//initbase.c
//higherorder_pointavg.c
//interpline.c
//interpline.mono.c  [nothing]
//interpline.para.c  [nothing]
//metric.c
//metric_selfgravity_or_evolvemetric.c
//mpi_fileio.c
//phys.ffde.c
//reconstructeno.c
//reconstructeno.weightmin.c
//restart.c
//restart.checks.c
//utoprimgen.c
//global.loops.boundaries.h
//global.variousmacros.h
//reconstructeno_static.h
//reconstructeno.debug.c
//restart.rebeccaoldcode.c
//
