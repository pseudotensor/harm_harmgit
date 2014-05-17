
/*! \file global.loops.diagnostics.h
    \brief Diagnostic Loops and definitions/macros
    // DIAGNOSTIC RELATED LOOPS

*/


/// want dump output to be ordered in radius first!!
/// consistent with MPI's use of BUFFERMAP to create buffer for *all* disk writing methods
/// This order should not change for any ORDERSTORAGE, but presumes i,j,k is standard indices for (i.e.) r(i), theta(j), phi(k)
#define DUMPLOOP(istart,istop,jstart,jstop,kstart,kstop)        \
  for(k=kstart;k<=kstop;k++)                                    \
    for(j=jstart;j<=jstop;j++)                                  \
      for(i=istart;i<=istop;i++)

#if(FULLOUTPUT==0)
#define EXTRADUMP1 0
#define EXTRADUMP2 0
#define EXTRADUMP3 0
#else
#define EXTRADUMP1T FULLOUTPUT*N1NOT1
#define EXTRADUMP2T FULLOUTPUT*N2NOT1
#define EXTRADUMP3T FULLOUTPUT*N3NOT1

#define EXTRADUMP1 ((EXTRADUMP1T>N1BND) ? N1BND : EXTRADUMP1T)
#define EXTRADUMP2 ((EXTRADUMP2T>N2BND) ? N2BND : EXTRADUMP2T)
#define EXTRADUMP3 ((EXTRADUMP3T>N3BND) ? N3BND : EXTRADUMP3T)

#endif

#if(FULLOUTPUT==0)
#define DUMPGENLOOP DUMPLOOP(0,N1-1,0,N2-1,0,N3-1)
#else
#define DUMPGENLOOP DUMPLOOP(-EXTRADUMP1,N1-1+EXTRADUMP1,-EXTRADUMP2,N2-1+EXTRADUMP2,-EXTRADUMP3,N3-1+EXTRADUMP3)
#endif


/// defines whether within the enerregion
/// considered loop-related macro
#define WITHINENERREGION(theenerpos,i,j,k) (i>=theenerpos[X1DN])&&(i<=theenerpos[X1UP])&&(j>=theenerpos[X2DN])&&(j<=theenerpos[X2UP])&&(k>=theenerpos[X3DN])&&(k<=theenerpos[X3UP]) 





///#define IMAGELOOP(istart,istop,jstart,jstop,kstart,kstop) for(k=kstart;k<=kstop;k++) for(j=jstart;j<=jstop;j++) for(i=istart;i<=istop;i++)

//#define OLDIMAGELOOP for(j=N2-1;j>=0;j--) for(i=0;i<N1;i++)
// nasty 
// to
// deal 
// with
