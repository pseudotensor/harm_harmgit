
/*! \file global.structs.h
    \brief macros and definitions related to HARM structures

// structure definitions
*/


/// none of below blink integers should be nearly as large as 32-bit signed integer can handle
struct blink {
  int num; // must stay int since used as argument to MPI functions that assume int -- so this limits number of elements one can pass.
  struct blink * np;
  // only used by cpu=0
  int cpu; // which cpu
  int i,j,k,col; // starting values for cpu=0
  int ri,rj,rk,rcol; // reference values for first cpu in sequence of nodes for a single buffer
  int end;
};



/// if add something, then should set it in (at least) set_grid.c
/// gcon put below gcov,gcovpert,alphalapse since gcon not as often needed
/// store betasqoalphasq to avoid need of gcon in most calculations
#define interiorofgeompart1a                    \
  int i,j,k,p;


#if(WHICHEOM==WITHGDET)
#define interiorofgeompart1b                    \
  FTYPE gdet;
#else
#define interiorofgeompart1b                    \
  FTYPE gdet;                                   \
  FTYPE eomfunc[NPR];
#endif

//  FTYPE gcov[SYMMATRIXNDIM];   \

#define interiorofgeompart1c                    \
  FTYPE gcov[SYMMATRIXNDIM];                    \
  FTYPE gcovpert[NDIM];                         \
  FTYPE alphalapse;                             \
  FTYPE betasqoalphasq;                         \
  FTYPE beta[NDIM];

#if(WHICHEOM==WITHGDET)
/// must always create eomfunc[NPR],ieomfuncnosing[NPR] even if(WHICHEOM!=WITHGDET) since code refers to these arrays[pl] in general
#define interiorofgeompart2                     \
  FTYPE igdetnosing;
#else
#define interiorofgeompart2                     \
  FTYPE igdetnosing;                            \
  FTYPE ieomfuncnosing[NPR];
#endif

//  FTYPE gcon[SYMMATRIXNDIM];

#if(GDETVOLDIFF)
#define interiorofgeompart3                     \
  FTYPE gdetvol;
#else
#define interiorofgeompart3
#endif

#define interiorofgeompart4                     \
  FTYPE gcon[SYMMATRIXNDIM];

// done with parts of structure

#define interiorofgeom                          \
  interiorofgeompart1a                          \
  interiorofgeompart1b                          \
  interiorofgeompart1c                          \
  interiorofgeompart2                           \
  interiorofgeompart3                           \
  interiorofgeompart4

#define interiorofgdetgeom                      \
  interiorofgeompart1b                          \
  interiorofgeompart2



/// stored global geometry used most of the time when previously would call get_geometry()
struct of_compgeom {

  interiorofgeom

};





#if(GETGEOMUSEPOINTER==0)

///typedef struct of_geom struct of_compgeom;
/// force to be the same
#define of_geom of_compgeom

struct of_allgeom {

  interiorofgeom

  // extra in "allgeom"
  FTYPE X[NDIM];
  FTYPE V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
};


#else


struct of_geom {

  // dummy space for gset() version
  FTYPE gengcov[SYMMATRIXNDIM];
  FTYPE gengcovpert[NDIM];
  FTYPE alphalapse;
  FTYPE betasqoalphasq;
  FTYPE beta[NDIM];
  FTYPE gengcon[SYMMATRIXNDIM];

  // bit faster since not all values always used
  FTYPE *gcov;
  FTYPE *gcon;
  FTYPE *gcovpert;

  FTYPE gdet,igdetnosing;
#if(GDETVOLDIFF)
  FTYPE gdetvol;
#endif
#if(WHICHEOM!=WITHGDET)
  FTYPE eomfunc[NPR],ieomfuncnosing[NPR];
#endif
  int i,j,k,p;
};


struct of_allgeom {

#if(GETGEOMUSEPOINTER==0)
  FTYPE gcov[SYMMATRIXNDIM];
  FTYPE gcovpert[NDIM];
  FTYPE alphalapse;
  FTYPE betasqoalphasq;
  FTYPE beta[NDIM];
  FTYPE gcon[SYMMATRIXNDIM];
#else
  // dummy space for gset() version
  FTYPE gengcon[SYMMATRIXNDIM];
  FTYPE gengcov[SYMMATRIXNDIM];
  FTYPE gengcovpert[NDIM];

  // bit faster since not all values always used
  FTYPE *gcov;
  FTYPE *gcovpert;
  FTYPE alphalapse;
  FTYPE betasqoalphasq;
  FTYPE beta[NDIM];
  FTYPE *gcon;
#endif

  FTYPE gdet,igdetnosing;
#if(GDETVOLDIFF)
  FTYPE gdetvol;
#endif
#if(WHICHEOM!=WITHGDET)
  FTYPE eomfunc[NPR],ieomfuncnosing[NPR];
#endif
  int i,j,k,p;

  // extra in "allgeom"
  FTYPE X[NDIM];
  FTYPE V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
};
#endif



#if(NEWMETRICSTORAGE)
/// stored global geometry used most of the time when previously would call get_geometry()
struct of_gdetgeom {

  interiorofgdetgeom

};
#else
#define of_gdetgeom of_compgeom
#endif




///////////////////////
///
/// state structure
///
///////////////////////
struct of_state {
  FTYPE ucon[NDIM];
  FTYPE ucov[NDIM];
#if(EOMRADTYPE!=EOMRADNONE)
  FTYPE uradcon[NDIM];
  FTYPE uradcov[NDIM];
#else
  FTYPE *uradcon;
  FTYPE *uradcov;
#endif
  FTYPE bcon[NDIM];
  FTYPE bcov[NDIM];
  FTYPE pressure; // aux thermodynamical quantity
  FTYPE bsq; // b^2 that is often used
  FTYPE entropy; //aux thermodynamical quantity
  FTYPE ifremoverestplus1ud0elseud0; // 1+u_t
  // OPTMARK: If don't use bcon,bcov, store Bcon/u^t and Bcov/u^t so avoid catastrophic cancellation but still avoid divisions.

  FTYPE others[NUMOTHERSTATERESULTS];

#if(EOMRADTYPE!=EOMRADNONE)
  FTYPE othersrad[NUMOTHERSTATERESULTS];
#else
  FTYPE *othersrad;
#endif


#if(MERGEDC2EA2CMETHOD)
  // for merged method and stored by compute_and_store_???() functions

  FTYPE gdet;
#if(WHICHEOM!=WITHGDET)
  FTYPE eomfunc[NPR]; // eomfunc
#endif
  FTYPE prim[NPR];
  FTYPE Blower[NDIM];
  FTYPE vcon[NDIM];
  FTYPE gdetBcon[NDIM];
  FTYPE overut;

#else
  // avoid allocating since expensive for cache misses
  // don't take more memory than needed for pointer references for code to compile
  //  FTYPE gdet;
  //#if(WHICHEOM!=WITHGDET)
  //  FTYPE *eomfunc;
  //#endif
  //  FTYPE *prim;
  //  FTYPE *Blower;
  //  FTYPE *gdetBcon; // for FLUXB==FLUXCTSTAG
  //  FTYPE *vcon;     // for FLUXB==FLUXCTSTAG
  //  FTYPE overut;

#endif
};




struct of_loop {
  int is, ie;
  int js, je;
  int ks, ke;
  int dir,intdir;
  int ps, pe;
  int bs, be;
  int di,dj,dk;
};

struct of_newtonstats {
  // outputs
  FTYPE lerrx;
  int lntries;
  int nstroke;
  FTYPE invproperty[NUMINVPROPERTY];
  char invpropertytext[NUMINVPROPERTY][10];
  // inputs
  FTYPE tryconv; // default is from u2p_defs.h: NEWT_TOL
  FTYPE tryconvultrarel; // default is from u2p_defs.h: NEWT_TOL_ULTRAREL
  FTYPE mintryconv; // default is from u2p_defs.h: MIN_NEWT_TOL
  int maxiter; // default is from u2p_defs.h: MAX_NEWT_ITER
  int extra_newt_iter; //  EXTRA_NEWT_ITER
  int extra_newt_iter_ultrarel; // EXTRA_NEWT_ITER_ULTRAREL
};

struct of_trueijkp {
  int i,j,k,p,dir,iter,interporflux;
};
