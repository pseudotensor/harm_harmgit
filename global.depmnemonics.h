
/*! \file global.depmnemonics.h
    \brief General code definitions of dependent quantities

    // MNEMONICS or other things that should rarely change, or things that on depend on the above items (e.g. if statements, loops, etc.)
    
*/


/// RK-related macros
/// Doesn't even depend upon N?, N?M, or N?BND, etc.  But does depend upon MAXTIMEORDER
#include "global.stepch.h"


///equals to unity if the interpolation of gamma is performed and if requested to use prim. reduction
#define STORE_GAMMA_PRIM_REDUCTION_FRACTION  (WENO_USE_PRIM_REDUCTION && (VARTOINTERP == PRIMTOINTERP_3VEL_GAMMA || VARTOINTERP == PRIMTOINTERP_RHOV_GAMMA || VARTOINTERP == PRIMTOINTERP_3VELREL_GAMMAREL || VARTOINTERP == PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP) )

#if( WENO_USE_PRIM_REDUCTION ) 
//#error WENO_USE_PRIM_REDUCTION is broken.  A variable is not defined -- need to correct. ** BROKEN **
#endif



#if( MERGEDC2EA2CMETHOD )
#define NUM_CVT_TYPES 8  //4 usual recon types + 4 derivatives
#else
#define NUM_CVT_TYPES 4
#endif


#if(EOMTYPE==EOMGRMHD||EOMTYPE==EOMCOLDGRMHD||EOMTYPE==EOMENTROPYGRMHD)
#define DOEVOLVERHO 1
#else
#define DOEVOLVERHO 0
#endif

#if(EOMTYPE==EOMGRMHD||EOMTYPE==EOMENTROPYGRMHD)
#define DOEVOLVEUU 1
#else
#define DOEVOLVEUU 0
#endif


/// whether to check if rho<=0 or u<=0 when restarting.
#if( DOEVOLVERHO||DOEVOLVEUU)
#define CHECKRHONEGZERORESTART 1
#else
#define CHECKRHONEGZERORESTART 0
#endif


/// how many sets of EOMs (MHD or MHD+RAD)
#if(EOMRADTYPE!=EOMRADNONE)
#define NUMTRUEEOMSETS (2)
//#define NUMEOMSETS (NUMTRUEEOMSETS) // MHD+RAD
#define NUMEOMSETS (NUMTRUEEOMSETS+1) // MHD+RADforflux + RADfordt
#else
#define NUMTRUEEOMSETS (1) // MHD
#define NUMEOMSETS (NUMTRUEEOMSETS)
#endif

#define EOMSETMHD 0
#define EOMSETRAD 1
#define EOMSETRADFORDT 2

#if(DIVERGENCEMETHOD==DIVMETHODPREFLUX)
#define NUMCOND 2
#else
#define NUMCOND 1
#endif

/// for STORESHOCKINDICATOR
#define NUMSHOCKPLS (3*(NUMTRUEEOMSETS*NUMCOND))
#define SHOCKPLDIR1 (0)
#define SHOCKPLDIR2 (SHOCKPLDIR1+1)
#define SHOCKPLDIR3 (SHOCKPLDIR2+1)
#define SHOCKRADPLDIR1 (SHOCKPLDIR3+1)
#define SHOCKRADPLDIR2 (SHOCKRADPLDIR1+1)
#define SHOCKRADPLDIR3 (SHOCKRADPLDIR2+1)
/// below only for DIVERGENCEMETHOD==DIVMETHODPREFLUX
#define DIVPLDIR1 (SHOCKRADPLDIR3+1)
#define DIVPLDIR2 (DIVPLDIR1+1)
#define DIVPLDIR3 (DIVPLDIR2+1)
#define DIVRADPLDIR1 (DIVPLDIR3+1)
#define DIVRADPLDIR2 (DIVRADPLDIR1+1)
#define DIVRADPLDIR3 (DIVRADPLDIR2+1)

/// for STORESHOCKINDICATOR for temp storage
/// assumes no more than NDIM*NUMTRUEEOMSETS in list!
#define SHOCKPLSTOREPTOT (0)
#define SHOCKPLSTOREVEL1 (1)
#define SHOCKPLSTOREVEL2 (SHOCKPLSTOREVEL1+1)
#define SHOCKPLSTOREVEL3 (SHOCKPLSTOREVEL2+1)
#define SHOCKRADPLSTOREPTOT (SHOCKPLSTOREVEL3+1)
#define SHOCKRADPLSTOREVEL1 (SHOCKRADPLSTOREPTOT+1)
#define SHOCKRADPLSTOREVEL2 (SHOCKRADPLSTOREVEL1+1)
#define SHOCKRADPLSTOREVEL3 (SHOCKRADPLSTOREVEL2+1)

#if(DIVERGENCEMETHOD==DIVMETHODPOSTFLUX)
#define NSPECIAL 6
/// choose SPECIALPL
#define SPECIALPL1 RHO
#define SPECIALPL2 UU
#define SPECIALPL3 B1
#define SPECIALPL4 B2
#define SPECIALPL5 B3
#define SPECIALPL6 URAD0
#else
#define NSPECIAL 0
#define SPECIALPL1 RHO // dummy, not used.
#define SPECIALPL2 UU // dummy, not used.
#define SPECIALPL3 B1
#define SPECIALPL4 B2
#define SPECIALPL5 B3
#define SPECIALPL6 URAD0
#endif

#define DOALLPL 0
#define DOBPL 1
#define DONONBPL 2
#define DOSPECIALPL 3

#define BPL(pl) (pl==B1 || pl==B2 || pl==B3)

// YFL for: rho,T10,T13,R10,R13
#define NUMYFL (1 + 2 + 2*(EOMRADTYPE!=EOMRADNONE))
#define YFLPL(pl) (pl==YFL1 || pl==YFL2 || pl==YFL3 || pl==YFL4 || pl==YFL5)
#define POSPL(pl) (pl==RHO || pl==UU || pl==URAD0 || YFLPL(pl) || pl==YL || pl==YNU)

#define SCALARPL(pl) (YFLPL(pl) || pl==YL || pl==YNU)

#define DENSITYPL(pl) (pl==RHO || pl==UU || pl==URAD0 || SCALARPL(pl))

#define RADPL(pl) (pl==PRAD0 || pl==PRAD1 || pl==PRAD2 || pl==PRAD3)
#define RADFULLPL(pl) (pl==PRAD0 || pl==PRAD1 || pl==PRAD2 || pl==PRAD3 || pl==YFL4 || pl==YFL5)

#define NONRADFULLPL(pl) (!RADFULLPL(pl))



#if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2))
#define NUMCURRENTSLOTS 5
#elif(WHICHCURRENTCALC==CURRENTCALC1)
#define NUMCURRENTSLOTS 6
#endif

#if(SPLITPRESSURETERMINFLUXMA)
/// where to put gas pressure part of flux that is T^i_i
#define FLUXSPLITPMA(dir) (B1+dir-1) // put in unused magnetic field part
#else
#define FLUXSPLITPMA(dir) (-100) // so FLUXSPLITPMA(dir)==pl is always 0
#endif

#if(SPLITPRESSURETERMINFLUXEM)
/// where to put magnetic pressure part of flux that is T^i_i
#define FLUXSPLITPEM(dir) (RHO) // put in unused mass advection part
#else
#define FLUXSPLITPEM(dir) (-100) // so FLUXSPLITPEM(dir)==pl is always 0
#endif


#if((USESTOREDSPEEDSFORFLUX==0)||(STOREWAVESPEEDS==0))
#define USEGLOBALWAVE 0
#else
#define USEGLOBALWAVE 1
#endif


/// OPTMARK: Should remove use of b^\mu b_\mu
#define COMPUTE4FIELDforFLUX (MAXWELL==GENMAXWELL || USEGLOBALWAVE==0)

/// OPTMARK: Should remove use of b^\mu b_\mu
#define COMPUTE4FIELDatALL (COMPUTE4FIELDforFLUX  || LIMITDTWITHSOURCETERM || ANALYTICSOURCE || PLINEWITHFIELD || VLINEWITHGDETRHO || TRUEFAST==1)




/// this should be only place eomfunc[] or ieomfuncnosing[] are referred to!
/// special cases: =GLOBALMETMACP1A0(eomfunc,loc,i,j,k);
/// one regexp to use: eomfunc\[\([_\>a-zA-Z0-9+-\ ()]+\)\] -> EOMFUNCMAC(\1)
/// one regexp to use: ieomfuncnosing\[\([_\>a-zA-Z0-9+-\ ()]+\)\] -> IEOMFUNCNOSINGMAC(\1)
#if(WHICHEOM!=WITHGDET)
#define EOMFUNCNAME eomfunc // requires special care within code -- only used when doing eomfunc_func()
#define EOMFUNCASSIGN(pl) eomfunc[pl] // requires special care within code -- only used when doing eomfunc_func()#define EOMFUNCPTR eomfuncptr // requires special care within code -- only used when doing eomfunc_func()
#define EOMFUNCPTR eomfuncptr
#define EOMFUNCMAC(pl) eomfunc[pl]
#define LOCALEOMFUNCMAC(pl) localeomfunc[pl] // used when did GETLOCALMETRIC
#define IEOMFUNCNOSINGMAC(pl) ieomfuncnosing[pl]
#else
#define EOMFUNCNAME gdet // requires special care within code -- only used when doing eomfunc_func()
#define EOMFUNCASSIGN(pl) (*gdet) // requires special care within code -- only used when doing eomfunc_func()
#define EOMFUNCPTR &gdet
#define EOMFUNCMAC(pl) gdet
#define LOCALEOMFUNCMAC(pl) localgdet[0] // used when did GETLOCALMETRIC
#define IEOMFUNCNOSINGMAC(pl) igdetnosing
#endif


///used by mpi_init.c for some bound types to set sign of copy for polar BCs
#define SIGNFLIPGDET (FLIPGDETAXIS==0 ? 1.0 : -1.0)

///used to set sign of U1,B1 across axis
#define SIGNFLIPU1 (FLIPU1AXIS==0 ? 1.0 : -1.0)
#define SIGNFLIPB1 (FLIPB1AXIS==0 ? 1.0 : -1.0)
///used to set sign of U2,B2 across axis
#define SIGNFLIPU2 (FLIPU2AXIS==0 ? 1.0 : -1.0)
#define SIGNFLIPB2 (FLIPB2AXIS==0 ? 1.0 : -1.0)
///used to set sign of U3,B3 across axis
#define SIGNFLIPU3 (FLIPU3AXIS==0 ? 1.0 : -1.0)
#define SIGNFLIPB3 (FLIPB3AXIS==0 ? 1.0 : -1.0)



/// setup for various boundary situations
/// so doesn't produce differences in irrelevant directions, whether boundary zones or not
/// mac(?) macros are for use in definitions of other macros since macro with args needs to be directly a function of what's hardcoded, not some "replacement" since nothing is replaced, not to be used inside code for USING a macro.
/// Commented out non-macro versios of (e.g.) im1 since want to enforce correct macro expansion when used with macrofication of array accesses described in global.storage.h

#if(N1>1)
//#define im1 (i-1)
#define im1mac(i) (i-1)
//#define ip1 (i+1)
#define ip1mac(i) (i+1)

//#define irefshift (2*N1-i-1)
#define irefshiftmac(i) (2*N1-i-1)

#else

//#define im1 (i)
#define im1mac(i) (i)
//#define ip1 (i)
#define ip1mac(i) (i)

//#define irefshift (i)
#define irefshiftmac(i) (i)

#endif

#if(N2>1)
//#define jm1 (j-1)
#define jm1mac(j) (j-1)
//#define jp1 (j+1)
#define jp1mac(j) (j+1)
//#define jrefshift (2*N2-j-1)
#define jrefshiftmac(j) (2*N2-j-1)
#else
//#define jm1 (j)
#define jm1mac(j) (j)
//#define jp1 (j)
#define jp1mac(j) (j)
//#define jrefshift (j)
#define jrefshiftmac(j) (j)
#endif

#if(N3>1)
//#define km1 (k-1)
#define km1mac(k) (k-1)
//#define kp1 (k+1)
#define kp1mac(k) (k+1)
//#define krefshift (2*N3-k-1)
#define krefshiftmac(k) (2*N3-k-1)
#else
//#define km1 (k)
#define km1mac(k) (k)
//#define kp1 (k)
#define kp1mac(k) (k)
//#define krefshift (k)
#define krefshiftmac(k) (k)
#endif


// used to tell if N>1 or not (can just ask directly)
#define N1NOT1 ((N1>1) ? 1 : 0)
#define N2NOT1 ((N2>1) ? 1 : 0)
#define N3NOT1 ((N3>1) ? 1 : 0)

// GODMARK: looks like I can set just above and rest are set for me for any case

// GODMARK: check these new conditions

// restrict loops only over relevant domain in reduced dimension case



// 3 maximum boundary zones needed if doing Parabolic interpolation
// maximum number of boundary zones needed for all calculations

// have to set manually if going to set DOENOFLUX at runtime.



/// x1
#define SHIFT1 N1NOT1
#define N1BND MAXBND*N1NOT1

#define INFULL1 (-N1BND)
#define INFULLP11 (-N1BND+SHIFT1)
#define OUTFULL1 (N1-1+N1BND) 
#define OUTFULLM11 (N1-1+N1BND-SHIFT1) 
#define OUTFULLP11 (N1-1+N1BND+SHIFT1)

#define INHALF1 (-N1BND/2)
#define OUTHALF1 (N1-1+N1BND/2)
#define INP11 (-N1BND+SHIFT1)
#define OUTP11 (N1-1+N1BND-SHIFT1)
#define INM1 -SHIFT1
#define OUTM1 N1-1+SHIFT1

/// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
/// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO1 (-N1BND)
#define INBOUNDHI1 (-1)

#define OUTBOUNDLO1 (N1)
#define OUTBOUNDHI1 (N1+N1BND-1)

//#define INFACEBOUNDLO1 (-N1BND) // (-N1BND+1) // GODMARK: large domain used for easy checking of fluxes after bound_flux().
//#define INFACEBOUNDHI1 (-1+SHIFT1)

/// up to -1 since 0 is actually defined with original primitives
/// generalize (expand) a bit:
#define INFACEBOUNDLO1 (-N1BND)
#define INFACEBOUNDHI1 (-1)


/// from N1+1 since N1 is actually defined with original primitives (only true if not FLUXCTSTAG)
/// generalize (expand) a bit:
//#define OUTFACEBOUNDLO1 (N1+1)
#define OUTFACEBOUNDLO1 (N1)
#define OUTFACEBOUNDHI1 (N1+N1BND-1)

/// x2
#define SHIFT2 N2NOT1
#define N2BND MAXBND*N2NOT1

#define INFULL2 (-N2BND)
#define INFULLP12 (-N2BND+SHIFT2)
#define OUTFULL2 (N2-1+N2BND)
#define OUTFULLM12 (N2-1+N2BND-SHIFT2) 
#define OUTFULLP12 (N2-1+N2BND+SHIFT2)

#define INHALF2 (-N2BND/2)
#define OUTHALF2 (N2-1+N2BND/2)
#define INP12 (-N2BND+SHIFT2)
#define OUTP12 (N2-1+N2BND-SHIFT2)
#define INM2 -SHIFT2
#define OUTM2 N2-1+SHIFT2

/// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
/// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO2 (-N2BND)
#define INBOUNDHI2 (-1)

#define OUTBOUNDLO2 (N2)
#define OUTBOUNDHI2 (N2+N2BND-1)

//#define INFACEBOUNDLO2 (-N2BND)
//#define INFACEBOUNDHI2 (-1+SHIFT2)

/// generalize (expand) a bit:
//#define INFACEBOUNDLO2 (-N2BND+1)
#define INFACEBOUNDLO2 (-N2BND)
#define INFACEBOUNDHI2 (-1)


/// generalize (expand) a bit:
//#define OUTFACEBOUNDLO2 (N2+1)
#define OUTFACEBOUNDLO2 (N2)
#define OUTFACEBOUNDHI2 (N2+N2BND-1)


/// x3
#define SHIFT3 N3NOT1
#define N3BND MAXBND*N3NOT1

#define INFULL3 (-N3BND)
#define INFULLP13 (-N3BND+SHIFT3)
#define OUTFULL3 (N3-1+N3BND)
#define OUTFULLM13 (N3-1+N3BND-SHIFT3)
#define OUTFULLP13 (N3-1+N3BND+SHIFT3)

#define INHALF3 (-N3BND/2)
#define OUTHALF3 (N3-1+N3BND/2)
#define INP13 (-N3BND+SHIFT3)
#define OUTP13 (N3-1+N3BND-SHIFT3)
#define INM3 -SHIFT3
#define OUTM3 N3-1+SHIFT3

/// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
/// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO3 (-N3BND)
#define INBOUNDHI3 (-1)

#define OUTBOUNDLO3 (N3)
#define OUTBOUNDHI3 (N3+N3BND-1)

//#define INFACEBOUNDLO3 (-N3BND)
//#define INFACEBOUNDHI3 (-1+SHIFT3)

/// generalize (expand) a bit:
//#define INFACEBOUNDLO3 (-N3BND+1)
#define INFACEBOUNDLO3 (-N3BND)
#define INFACEBOUNDHI3 (-1)


/// generalize (expand) a bit:
//#define OUTFACEBOUNDLO3 (N3+1)
#define OUTFACEBOUNDLO3 (N3)
#define OUTFACEBOUNDHI3 (N3+N3BND-1)






/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)

#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

/// N?OFF and N?NOT1 are a bit redundant
//#define N1OFF (((N1BND>0)&&(N1>1)) ? 1 : 0)
//#define N2OFF (((N2BND>0)&&(N2>1)) ? 1 : 0)
//#define N3OFF (((N3BND>0)&&(N3>1)) ? 1 : 0)




#if(LIMADJUST!=LIMITERFIXED && FLUXADJUST!=FLUXFIXED)
#define NUMPFLAGSBOUND (NUMPFLAGS) // all
#elif(LIMADJUST!=LIMITERFIXED)
#define NUMPFLAGSBOUND (1+FLAGREALLIM) // only 0,1,2
#elif(EOMRADTYPE!=EOMRADNONE)
#define NUMPFLAGSBOUND (1+FLAGUTOPRIMRADFAIL) // only 0,1
#else
#define NUMPFLAGSBOUND (1+FLAGUTOPRIMFAIL) // only 0
#endif


// number of boundary cells to choose for setting and MPI-copying and using of pflag boundary cells
#define NUMPFLAGBND1 N1BND
#define NUMPFLAGBND2 N2BND
#define NUMPFLAGBND3 N3BND
#define MAXNUMPFLAGBND (MAX(MAX(NUMPFLAGBND1,NUMPFLAGBND2),NUMPFLAGBND3))


/// GODMARK: Could make a volume that is not NBIGBND*NBIGSM but may be smaller?
/// used in init_mpi.c for workbc and workbc_int


/// number of interpolation variables for staggered field method per B and v each
#define NUMCORNINTERP 4

/// for wavespeeds.c and fluxcompute.c
#define NUMCS 2
#define CMIN 0
#define CMAX 1

#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )





#if((DODISS||DODISSVSR)&&(DOENTROPY==DONOENTROPY))
#error Turn on entropy evolution if want dissipation
#endif




/// processed version of npr definition since no immediate evaluation in standard c preprocessor
/// this "global.defnprs.h" is generated by the program generatenprs.c
#include "global.defnprs.h"




#define NMAXBOUND ((NPRBOUND>NFLUXBOUND) ? NPRBOUND : NFLUXBOUND)





// total + pake + en + em + rad
#define NUMPHYSICALFLUXTERMS (1  +  1+1+1 + (EOMRADTYPE!=EOMRADNONE))
#define NUMFLUXESTOSAVE (1 + 2 + (EOMRADTYPE!=EOMRADNONE)*2 + NUMYFL*(DOYFL!=0) + (DOYL!=0) + (DOYNU!=0) )
#define FLUXESTOSAVEPL(pl) (pl==RHO || pl==UU || pl==U3 || pl==URAD0 || pl==URAD3 || YFLPL(pl) || pl==YL || pl==YNU)

#if(FLUXDUMP==0)
#define NUMFLUXDUMP (1)

#elif(FLUXDUMP==1)
/// cent,face1,face2,face3,corn
/// NPR*4 = 1 dUgeom and 3 dUriemanns ; 3 directions for F1,F2,F3 and pl pr and F(pl) and F(pr)
#define NUMFLUXDUMP (NPR*4 + NPR*3*(1+2+2))

#else
/// total,pake,en,em
#define NUMFLUXDUMP (NPR*(NUMPHYSICALFLUXTERMS))

#endif

#if(MODIFYEMFORVPOT==MODIFYVPOT || TRACKVPOT>0 || EVOLVEWITHVPOT>0)
/// 4 space-time directions with only spatial parts used for now
/// vpotarrayglobal holds vpot, vpot0, vpotlast, vpotcum
#define NUMVPOT (NDIM*4)
#define NUMVPOTDUMP (NDIM)
#else
#define NUMVPOT (0)
#define NUMVPOTDUMP (0)
#endif


/// size of certain dumped tavg quantities
/// was 29 =>
#define NUMNORMDUMP (NPR+1+4*4+6) // number of "normal" dump variables
/// for above see diag.c and set_varstavg()
#define NUMFARADAY 6
#define NUMOTHER 1
#define NUMSTRESSTERMS (NUMFLUXTERMS*NDIM*NDIM)

/** GLOBAL ARRAY SECTION **/






/// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1
#define LONGDOUBLETYPE 2
#define LONGLONGINTTYPE 3


#if(REALTYPE>SENSITIVE)
god=deathadflkjasdflkjasdlfkja242424
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON (1.19209290e-07F)
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON (2.2204460492503131e-16L)
#endif

#ifndef LDBL_EPSILON
#define LDBL_EPSILON (1.08420217248550443401e-19L)
#endif


  // need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define MINNUMREPRESENT FLT_MIN
#define NUMEPSILON FLT_EPSILON
#elif(REALTYPE==DOUBLETYPE)
#define MINNUMREPRESENT DBL_MIN
#define NUMEPSILON DBL_EPSILON
#elif(REALTYPE==LONGDOUBLETYPE)
#define MINNUMREPRESENT DBL_MIN
#define NUMEPSILON LDBL_EPSILON
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPE double
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPE long double
#endif


  // used for numerical differencing
#define NUMSQRTEPSILON (sqrt(NUMEPSILON))
#define NUMSQEPSILON (NUMEPSILON*NUMEPSILON)

  // for finite differences.
  // If one uses (say) 1E-5 for DX and has (Vh-Vl)/(Xh-Xl), then even if Xh-Xl is machine representable, Xh and Xl may not be.
#define MY1EM4 (1.0/(2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)) // 1/2^{13} \approx 1.2207E-4\sim 1E-4.

#define MY1EM5 (1.0/(2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)) // 1/2^{17} \approx 7.6E-6\sim 1E-5.

#define MY1EM6 (1.0/(2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)) // 1/2^{20} \approx 7.6E-6\sim 1E-6.


#if(COUNTTYPE==DOUBLETYPE)
#define CTYPE double
#define CTYPEHEADERONEOUT "%26.20g"
#define CTYPEHEADERONEIN "%lf"
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define CTYPE long long int
#define CTYPEHEADERONEOUT "%lld"
#define CTYPEHEADERONEIN "%lld"
#endif


#if(PFLAGTYPE==INTTYPE)
#define PFTYPE signed int
#define PFTYPEHEADERONEOUT "%d"
#define PFTYPEHEADERONEIN "%d"
#elif(PFLAGTYPE==CHARTYPE)
  // on Harvard's BlueGene/L char is by default unsigned!, so force signed as required by our code
#define PFTYPE signed char
#define PFTYPEHEADERONEOUT "%d"
#define PFTYPEHEADERONEIN "%d"
#endif


  // GODMARK: NUMENERVAR outdated?
#define NUMENERVAR (6+NPR+NPR+3)


  /* numerical convenience */
#if(REALTYPE==FLOATTYPE)
#define VERYBIG (1.e37)
#define BIG (1.e+30)
#define SMALL (1.e-35)
#define KINDASMALL (1.e-30)
#else
#define VERYBIG (1.e150)
#define BIG (1.e+100)
#define SMALL (1.e-300)
#define KINDASMALL (1.e-60)
#endif


#define SLEPSILON (1.e-6)


  /* size of step in numerical derivative evaluations */
#define HSTEP (1.e-5)




#if(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPEHEADERONEIN "%Lf"
  // assume sensitive>=realtype in precision
#if(REALTYPE==LONGDOUBLETYPE) // was FLOATTYPE==REALTYPE and SENS=DOUBLETYPE
#define HEADERONEIN "%Lf"
#define HEADER2IN "%Lf %Lf"
#define HEADER3IN "%Lf %Lf %Lf"
#define HEADER4IN "%Lf %Lf %Lf %Lf"
#define HEADER5IN "%Lf %Lf %Lf %Lf %Lf"
#define HEADER6IN "%Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER7IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER8IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER9IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER10IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER14IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER17IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER18IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER19IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER21IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define RESTARTHEADER "%d %d %d "                                       \
    "%Lf %Lf %ld %Lf %Lf %Lf %Lf %Lf "                                  \
    "%Lf %d %d %d %d %d %d %d %d "                                      \
    "%Lf %Lf %Lf %Lf %d "                                               \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#elif(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER10IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER14IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER17IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER18IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER19IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER21IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "                                       \
    "%Lf %Lf %ld %Lf %Lf %Lf %lf %lf "                                  \
    "%Lf %d %d %d %d %d %d %d %d "                                      \
    "%lf %lf %lf %lf %d "                                               \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define HEADER10IN "%f %f %f %f %f %f %f %f %f %f"
#define HEADER14IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER17IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER18IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER19IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER21IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "                                       \
    "%Lf %Lf %ld %Lf %Lf %Lf %f %f "                                    \
    "%Lf %d %d %d %d %d %d %d %d "                                      \
    "%f %f %f %f %d "                                                   \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPEHEADERONEIN "%lf"
  // assume sensitive>=realtype in precision
#if(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER10IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER14IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER17IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER18IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER19IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER21IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "                                       \
    "%lf %lf %ld %lf %lf %lf %lf %lf "                                  \
    "%lf %d %d %d %d %d %d %d %d "                                      \
    "%lf %lf %lf %lf %d "                                               \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define HEADER10IN "%f %f %f %f %f %f %f %f %f %f"
#define HEADER14IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER17IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER18IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER19IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER21IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "                                       \
    "%lf %lf %ld %lf %lf %lf %f %f "                                    \
    "%lf %d %d %d %d %d %d %d %d "                                      \
    "%f %f %f %f %d "                                                   \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==FLOATTYPE)
#define SFTYPEHEADERONEIN "%f"

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define HEADER10IN "%f %f %f %f %f %f %f %f %f %f"
#define HEADER14IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER17IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER18IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER19IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define HEADER21IN "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "                                       \
    "%f %f %ld %f %f %f %f %f "                                         \
    "%f %d %d %d %d %d %d %d %d "                                       \
    "%f %f %f %f %d "                                                   \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#endif


  // 2
  // 6
  // 10
  // 3
  // 5
  // 4
  // SUM=30
  // 23+10=33
  // total=63
#define WRITERESTARTHEADER "%lld %lld %lld "                            \
    "%26.20g %26.20g %ld %26.20g %26.20g %26.20g %26.20g %26.20g "      \
    "%26.20g %d %d %d %d %d %d %d %d "                                  \
    "%26.20g %26.20g %26.20g %26.20g %d "                               \
    "%d %d %d %d %d %d "                                                \
    "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g "

#define HEADERONEOUT "%26.20g"

#if(SENSITIVE==FLOATTYPE)
#define SFTYPEHEADERONEOUT "%31.26g"
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPEHEADERONEOUT "%31.26g"
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPEHEADERONEOUT "%31.26Lg"
#endif



#include "global.depmnemonics.rad.h" // KORAL
