
/*! \file global.loops.manypoints1d.h
    \brief Multi-point 1D loop related definitions/macros

    // full CPU 1-D loops
    
*/


#define NUMGRAVPOS (ncpux1*N1+N1BND*2+1)
#define GRAVLOOP(ii) for(ii=-N1BND;ii<ncpux1*N1+N1BND+1;ii++)
#define GRAVLOOPACTIVE(ii) for(ii=0;ii<ncpux1*N1;ii++)
#define DUMPGRAVLOOP(ii) GRAVLOOPACTIVE(ii) // normal
//#define DUMPGRAVLOOP(ii) GRAVLOOP(ii) // abnormal, debugging stuff in SM
