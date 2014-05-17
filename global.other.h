
/*! \file global.other.h
    \brief Other macros and definitions currently reltaed to ENO only
*/

#if( DO_WENO_DEBUG )
///debug atch SASMARK
FILE *debugfp;
#endif //DO_WENO_DEBUG

//defines of the interpolation direction unit vector -- can be used outside of ENODEBUG stuff, so put it here
//#define ITERGLOBALDEF {di = (iterglobal==1); dj = (iterglobal==2); dk = (iterglobal==3);}

#if(DOENODEBUG) //atch 
#define ENODEBUGPARAM_SMONO 0
#define ENODEBUGPARAM_WENO5 1
#define ENODEBUGPARAM_WENO3 2
#define ENODEBUGPARAM_dPP 3
#define ENODEBUGPARAM_LIMCORR 4
#define ENODEBUGPARAM_LIMCORR_PRIM 5
#endif 
