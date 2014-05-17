
/*! \file global.loops.manypoints.h
    \brief Multi-point loop related definitions/macros

    // LOOPS (not computational, but not per-point.  So still over i,j,k but diagnostic or other).
    
*/

//////////////////////////////////////
//////////////////////////////////////
///
/// MOST GENERAL LOOP copy.  While fast, does not allow for control over loop order in memory.
///
//////////////////////////////////////
//////////////////////////////////////


/// Sasha's dad and google say memcopy() and memmove() can be much faster than loop if compiler doesn't recognize this optimization.  So best to force it in simple cases
#define USE_MEMCPY 1

/// applies to continugous memory regions (i.e. feed single pointer for source, size of data to copy in its own dimensions)
#if(USE_MEMCPY)
// iter not used
#define GENFORALL(iter,src,dest,numelem) memcpy(dest, src, numelem*sizeof(src[0]))
#define GENFORALLOVERLAP(iter,src,dest,numelem) memmove(dest, src, numelem*sizeof(src[0]))
#else
/// assumes iteration is single (i.e. all data from start to finish)
#define GENFORALL(iter,src,dest,numelem)        \
  for (iter = 0; iter < numelem; iter++)        \
    {                                           \
      dest[iter] = src[iter];                   \
    }
#define GENFORALLOVERLAP(iter,src,dest,numelem) GENFORALL(iter,src,dest,numelem)
#endif


//////////////////////////////////////
//////////////////////////////////////
///
/// Many-point loops (but not computational)
///
//////////////////////////////////////
//////////////////////////////////////



////////////////////
///
/// below are each per-spatial-dimension loops used later to contruct a multi-dimensionan loop of any order
///
////////////////////


/// these loops used for general purposes
#define LOOPF3 for(k=INFULL3;k<=OUTFULL3;k++)
#define LOOPF2 for(j=INFULL2;j<=OUTFULL2;j++)
#define LOOPF1 for(i=INFULL1;i<=OUTFULL1;i++)

/// only used to initialize emf[]
#define LOOPFPM3 for(k=INFULL3-SHIFT3;k<=OUTFULL3+SHIFT3;k++)
#define LOOPFPM2 for(j=INFULL2-SHIFT2;j<=OUTFULL2+SHIFT2;j++)
#define LOOPFPM1 for(i=INFULL1-SHIFT1;i<=OUTFULL1+SHIFT1;i++)


/// full loop + 1 on outer edge for emf or corner quantities
#define LOOPFP13 for(k=INFULL3;k<=OUTFULLP13;k++)
#define LOOPFP12 for(j=INFULL2;j<=OUTFULLP12;j++)
#define LOOPFP11 for(i=INFULL1;i<=OUTFULLP11;i++)


/// full loop + 1 (shift away from boundary) on inner edge for comptuing emf or corner quantities
//#define LOOPINFP13 for(k=INFULLP13;k<=OUTFULL3;k++)
//#define LOOPINFP12 for(j=INFULLP12;j<=OUTFULL2;j++)
//#define LOOPINFP11 for(i=INFULLP11;i<=OUTFULL1;i++)


//#define LOOPOUTFM13 for(k=INFULL3;k<=OUTFULLM13;k++)
//#define LOOPOUTFM12 for(j=INFULL2;j<=OUTFULLM12;j++)
//#define LOOPOUTFM11 for(i=INFULL1;i<=OUTFULLM11;i++)

#define LOOPH3 for(k=INHALF3;k<=OUTHALF3;k++)
#define LOOPH2 for(j=INHALF2;j<=OUTHALF2;j++)
#define LOOPH1 for(i=INHALF1;i<=OUTHALF1;i++)

#define LOOPP13 for(k=INP13;k<=OUTP13;k++)
#define LOOPP12 for(j=INP12;j<=OUTP12;j++)
#define LOOPP11 for(i=INP11;i<=OUTP11;i++)

#define LOOPN3 for(k=0;k<=N3-1;k++)
#define LOOPN2 for(j=0;j<=N2-1;j++)
#define LOOPN1 for(i=0;i<=N1-1;i++)

#define LOOPFMHP3 for(k=INFULL3;k<=OUTHALF3;k++)
#define LOOPFMHP2 for(j=INFULL2;j<=OUTHALF2;j++)
#define LOOPFMHP1 for(i=INFULL1;i<=OUTHALF1;i++)

#define LOOPHMFP3 for(k=INHALF3;k<=OUTFULL3;k++)
#define LOOPHMFP2 for(j=INHALF2;j<=OUTFULL2;j++)
#define LOOPHMFP1 for(i=INHALF1;i<=OUTFULL1;i++)

#define LOOPHP3 for(k=0;k<=OUTHALF3;k++)
#define LOOPHP2 for(j=0;j<=OUTHALF2;j++)
#define LOOPHP1 for(i=0;i<=OUTHALF1;i++)

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)

/// SUPERGEN and GEN are super general and general, but still to be used with only spatial arrays related to spatial (N1,N2,N3) directions in i,j,k.  Can create non-spatial general loop if required
/// PACKLOOP() in boundmpi.c and PACKLOOP_INT() boundmpiint.c uses these correctly since refer to i,j,k in correct order
/// SUPERGENLOOP() as used in interpline.c is correct since always refer to i,j,k in correct order
/// GENLOOP() in diag.c correctly uses i,j,k in order
#define SUPERGENLOOP1(i,istart,istop,di) for((i)=(istart);(di>0 ? (i)<=(istop) : (i)>=(istop)); (i)+=(di))
#define SUPERGENLOOP2(j,jstart,jstop,dj) for((j)=(jstart);(dj>0 ? (j)<=(jstop) : (j)>=(jstop)); (j)+=(dj))
#define SUPERGENLOOP3(k,kstart,kstop,dk) for((k)=(kstart);(dk>0 ? (k)<=(kstop) : (k)>=(kstop)); (k)+=(dk))

#define GENLOOP1(i,istart,istop) SUPERGENLOOP1(i,istart,istop,1)
#define GENLOOP2(j,jstart,jstop) SUPERGENLOOP2(j,jstart,jstop,1)
#define GENLOOP3(k,kstart,kstop) SUPERGENLOOP3(k,kstart,kstop,1)

#define ZSLOOP1(istart,istop) SUPERGENLOOP1(i,istart,istop,1)
#define ZSLOOP2(jstart,jstop) SUPERGENLOOP2(j,jstart,jstop,1)
#define ZSLOOP3(kstart,kstop) SUPERGENLOOP3(k,kstart,kstop,1)

#define LOOPC3 LOOPN3
#define LOOPC2 LOOPN2
#define LOOPC1 LOOPN1






////////////////////
///
/// below are different ways of combining each spatial direction
///
////////////////////



/// general loop for any indicies
#define SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) LOOPORDER1(SUPERGENLOOP1(i,istart,istop,di),SUPERGENLOOP2(j,jstart,jstop,dj),SUPERGENLOOP3(k,kstart,kstop,dk)) LOOPORDER2(SUPERGENLOOP1(i,istart,istop,di),SUPERGENLOOP2(j,jstart,jstop,dj),SUPERGENLOOP3(k,kstart,kstop,dk)) LOOPORDER3(SUPERGENLOOP1(i,istart,istop,di),SUPERGENLOOP2(j,jstart,jstop,dj),SUPERGENLOOP3(k,kstart,kstop,dk))
/// general loop for any indicies
#define GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) LOOPORDER1(GENLOOP1(i,istart,istop),GENLOOP2(j,jstart,jstop),GENLOOP3(k,kstart,kstop)) LOOPORDER2(GENLOOP1(i,istart,istop),GENLOOP2(j,jstart,jstop),GENLOOP3(k,kstart,kstop)) LOOPORDER3(GENLOOP1(i,istart,istop),GENLOOP2(j,jstart,jstop),GENLOOP3(k,kstart,kstop))
/// general loop, but assumes i,j,k used
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) LOOPORDER1(ZSLOOP1(istart,istop),ZSLOOP2(jstart,jstop),ZSLOOP3(kstart,kstop)) LOOPORDER2(ZSLOOP1(istart,istop),ZSLOOP2(jstart,jstop),ZSLOOP3(kstart,kstop)) LOOPORDER3(ZSLOOP1(istart,istop),ZSLOOP2(jstart,jstop),ZSLOOP3(kstart,kstop))
/// below used for initialization and such, not a computational issue
#define LOOPF LOOPORDER1(LOOPF1,LOOPF2,LOOPF3) LOOPORDER2(LOOPF1,LOOPF2,LOOPF3) LOOPORDER3(LOOPF1,LOOPF2,LOOPF3)

#define LOOPFPM LOOPORDER1(LOOPF1,LOOPFPM2,LOOPFPM3) LOOPORDER2(LOOPFPM1,LOOPFPM2,LOOPFPM3) LOOPORDER3(LOOPFPM1,LOOPFPM2,LOOPFPM3)

#define FULLLOOPPM LOOPFPM

#define LOOPF_12 LOOPORDER1(LOOPF1,LOOPF2,) LOOPORDER2(LOOPF1,LOOPF2,) LOOPORDER3(LOOPF1,LOOPF2,)
#define LOOPF_13 LOOPORDER1(LOOPF1,,LOOPF3) LOOPORDER2(LOOPF1,,LOOPF3) LOOPORDER3(LOOPF1,,LOOPF3)
#define LOOPF_23 LOOPORDER1(,LOOPF2,LOOPF3) LOOPORDER2(,LOOPF2,LOOPF3) LOOPORDER3(,LOOPF2,LOOPF3)
#define LOOPH LOOPORDER1(LOOPH1,LOOPH2,LOOPH3) LOOPORDER2(LOOPH1,LOOPH2,LOOPH3) LOOPORDER3(LOOPH1,LOOPH2,LOOPH3)
#define LOOPP1 LOOPORDER1(LOOPP11,LOOPP12,LOOPP13) LOOPORDER2(LOOPP11,LOOPP12,LOOPP13) LOOPORDER3(LOOPP11,LOOPP12,LOOPP13)
#define LOOP LOOPORDER1(LOOPN1,LOOPN2,LOOPN3) LOOPORDER2(LOOPN1,LOOPN2,LOOPN3) LOOPORDER3(LOOPN1,LOOPN2,LOOPN3)
#define LOOPFMHP LOOPORDER1(LOOPFMHP1,LOOPFMHP2,LOOPFMHP3) LOOPORDER2(LOOPFMHP1,LOOPFMHP2,LOOPFMHP3) LOOPORDER3(LOOPFMHP1,LOOPFMHP2,LOOPFMHP3)
#define LOOPHMFP LOOPORDER1(LOOPHMFP1,LOOPHMFP2,LOOPHMFP3) LOOPORDER2(LOOPHMFP1,LOOPHMFP2,LOOPHMFP3) LOOPORDER3(LOOPHMFP1,LOOPHMFP2,LOOPHMFP3)
#define LOOPHP LOOPORDER1(LOOPHP1,LOOPHP2,LOOPHP3) LOOPORDER2(LOOPHP1,LOOPHP2,LOOPHP3) LOOPORDER3(LOOPHP1,LOOPHP2,LOOPHP3)
#define LOOPC LOOPORDER1(LOOPC1,LOOPC2,LOOPC3) LOOPORDER2(LOOPC1,LOOPC2,LOOPC3) LOOPORDER3(LOOPC1,LOOPC2,LOOPC3)
#define FULLLOOP LOOPF
/// larger loop than full for cornered quantities such as emf defined on corners that need to be initialized for boundary condition reasons
#define FULLLOOPP1 LOOPORDER1(LOOPFP11,LOOPFP12,LOOPFP13) LOOPORDER2(LOOPFP11,LOOPFP12,LOOPFP13) LOOPORDER3(LOOPFP11,LOOPFP12,LOOPFP13)
/// divb loop (for diagnostics only)
//#define LOOPDIVB LOOPORDER1(LOOPP11,LOOPP12,LOOPP13) LOOPORDER2(LOOPP11,LOOPP12,LOOPP13) LOOPORDER3(LOOPP11,LOOPP12,LOOPP13)
/// boundary zones may not require divb=0 since proxy for flux
#define LOOPDIVB LOOPORDER1(LOOPC1,LOOPC2,LOOPC3) LOOPORDER2(LOOPC1,LOOPC2,LOOPC3) LOOPORDER3(LOOPC1,LOOPC2,LOOPC3)









////////////////////
///
/// below are not original multi-D combindations, just renamings that have no control over order of dimensions
///
////////////////////
#define LOOPFC LOOPF
#define LOOPHC LOOPH
#define LOOPFMHPC LOOPFMHP
#define LOOPHMFPC LOOPHMFP
#define LOOPHPC LOOPHP
#define ZLOOP ZSLOOP(0,N1-1,0,N2-1,0,N3-1)











///////////////////
//
// Below are deprecated -- if activate ensure proper general loop behavior as above
//
///////////////////

//#if(BOUNDPLPR&&NOFLUXCTONX1DN)
//#define FULLLOOPflux1 LOOPF1 LOOPF2 LOOPF3 if(!(startpos[1]+i==0))
//#else
//#define FULLLOOPflux1 LOOPF1 LOOPF2 LOOPF3
//#endif

// computing emf for FLUXCT
//#define LOOPINFP1 LOOPINFP11 LOOPINFP12 LOOPINFP13

//#define LOOPINFP1dir1full LOOPF1 LOOPINFP12 LOOPINFP13

//#define LOOPINFP1dir2full LOOPINFP11 LOOPF2 LOOPINFP13

//#define LOOPINFP1dir3full LOOPINFP11 LOOPINFP12 LOOPF3

//#define LOOPINFP1dir23full LOOPINFP11 LOOPF2 LOOPF3

//#define LOOPINFP1dir13full LOOPF1 LOOPINFP12 LOOPF3

//#define LOOPINFP1dir12full LOOPF1 LOOPF2 LOOPINFP13


// computing emf for FLUXCD
//#define LOOPOUTFM1 LOOPOUTFM11 LOOPOUTFM12 LOOPOUTFM13

//#if(BOUNDPLPR&&NOFLUXCTONX1DN)
//#define LOOPOUTFM1dir1fullflux LOOPF1 LOOPOUTFM12 LOOPOUTFM13 if(!(startpos[1]+i==0))
//#else
//#define LOOPOUTFM1dir1fullflux LOOPF1 LOOPOUTFM12 LOOPOUTFM13
//#endif

//#define LOOPOUTFM1dir1full LOOPF1 LOOPOUTFM12 LOOPOUTFM13

//#define LOOPOUTFM1dir2full LOOPOUTFM11 LOOPF2 LOOPOUTFM13

//#define LOOPOUTFM1dir3full LOOPOUTFM11 LOOPOUTFM12 LOOPF3

//#define WSPEEDLOOP ZSLOOP(-SHIFT1,N1-1+2*SHIFT1,-SHIFT2,N2-1+2*SHIFT2,-SHIFT3,N3-1+2*SHIFT3)

//#define PLUSLOOP ZSLOOP(-1, N1, -1, N2, -1, N3)

// below same as FULLLOOP if NBND=2
//#define PLUSPLUSLOOP ZSLOOP(-2, N1+1, -2, N2+1, -2, N3+1)






