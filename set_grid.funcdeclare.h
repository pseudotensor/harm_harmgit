
/*! \file set_grid.funcdeclare.h
     \brief function declarations for global use of set_grid.c functions
*/

// called in restart.c and initbase.c
extern void set_grid(int whichtime,FTYPE *CUf, FTYPE *Cunew);



extern int assignmetricstorage_new(struct of_compgeom *mygeom, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);
extern int assignmetricstorage_old(int loc, int i, int j, int k, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);
extern int assignmetricstorage_oldlast(int loc, int i, int j, int k, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);

#if(NEWMETRICSTORAGE)
#define GETLOCALMETRIC(loc,i,j,k)     assignmetricstorage_new(&GLOBALMETMACP1A0(compgeom,loc,i,j,k), &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#define GETLASTLOCALMETRIC(loc,i,j,k) assignmetricstorage_new(&GLOBALMETMACP1A0(compgeomlast,loc,i,j,k), &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#else
// uses globals
#define GETLOCALMETRIC(loc,i,j,k)     assignmetricstorage_old(loc, i, j, k, &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#define GETLASTLOCALMETRIC(loc,i,j,k) assignmetricstorage_oldlast(loc, i, j, k, &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#endif

//  FTYPE (*localgcov)[NDIM];   \
//  FTYPE (*localgcon)[NDIM];   \

#define LOCALMETRICTEMPVARS                     \
  FTYPE *localgcov;                             \
  FTYPE *localgcon;                             \
  FTYPE *localgcovpert;                         \
  FTYPE *localgdet,*localgdetvol;               \
  FTYPE *localalphalapse;                       \
  FTYPE *localbetasqoalphasq;                   \
  FTYPE *localbeta;                             \
  FTYPE *localeomfunc;
