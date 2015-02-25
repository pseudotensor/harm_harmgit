


/*! \file interpline.c
     \brief Spatial Interpolation for fluxes based upon providing full 1D line information
     // Instead of acting per point, this acts per line to improve memory efficiency.
*/


#include "decs.h"







#if(MERGEDC2EA2CMETHOD)
#define LINETYPEDEFINESOTHER FTYPE a_youtpolycoef[NPR2INTERP][MAXSPACEORDER][NBIGM];
#error "SUPERGODMARK: Should only need this if doing more than 3 point stencil -- fix"
#else
#define LINETYPEDEFINESOTHER FTYPE a_youtpolycoef[1][1][1]; //VS complained about defining zero-size object
#endif


/// NPR2INTERP always larger than NPR, so can use one memory space for both c2e and others
#define LINETYPEDEFINES1                                                \
  int bs,be,ps,pe;                                                      \
  int di,dj,dk;                                                         \
  int is,ie,js,je,ks,ke;                                                \
  int preforder, whichreduce;                                           \
  int pl,pliter;                                                        \
  int recontype;                                                        \
  extern void pass_1d_line_weno(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp); \
  extern void pass_1d_line_paraline(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp); \
  extern void pass_1d_line_multipl_weno(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp); \
  extern void pass_1d_line_multipl_paraline(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp); \
  int intdir;                                                           \
  int withshifts;                                                       \
  int nprlocalstart,nprlocalend;                                        \
  int nprlocallist[MAXNPR];                                             \
  int pllocal;                                                          \
  int numprims;                                                         \
  int stencilvarisnull;                                                 \
  int doingweno; // last thing has no semicolon



#define GEN3MAC(numypl,pl,numywhich,which,numyi,i) (i + numyi*which + numyi*numywhich*pl)

/// below for drho,dP,etai
#define GEN2MAC(numypl,pl,numyi,i) (i + numyi*pl)
/// below for shockindicator, stiffindicator, Pline, Vline, shift, minorder, maxorder
#define GEN1MAC(numyi,i) (i)

/// below to be used for yin, yout, yprim, ystencilvar that involve [NPR2INTERP][2][NBIGM] type arrays previously
#define YINOUTMAC(numypl,pl,numywhich,which,numyi,i) GEN3MAC(numypl,pl,numywhich,which,numyi,i)



/// For OpenMP must have different memory per thread, so why put inside parallel region
/// However, allocating and deallocating these perijk() is very expensive, so moved back outside perijk() and pass related pointers.
/// OPENMPMARK: Still must allocate within parallel region!  So place this within parallel region.
#define LINETYPEDEFINESMEMORY                                   \
  FTYPE a_yin[NPR2INTERP][2][NBIGM];                            \
  FTYPE a_yout[NPR2INTERP][2][NBIGM];                           \
                                                                \
  FTYPE a_shockindicator[NUMTRUEEOMSETS][NBIGM];                 \
  FTYPE a_stiffindicator[NBIGM];                                \
                                                                \
  FTYPE a_yprim[NPR2INTERP][2][NBIGM];                          \
  FTYPE a_ystencilvar[NPR2INTERP][2][NBIGM];                    \
  FTYPE a_df[NPR2INTERP][NUMDFS][NBIGM];                        \
  FTYPE a_dfformono[NPR2INTERP][NUMDFS][NBIGM];                 \
  FTYPE a_drho[NUMDFS][NBIGM];                                  \
  FTYPE a_dP[NUMDFS][NBIGM];                                    \
  FTYPE a_etai[NPR2INTERP][NUMTRUEEOMSETS][NBIGM];              \
  FTYPE a_Pline[NUMTRUEEOMSETS][NBIGM];                          \
  FTYPE a_Vline[NUMTRUEEOMSETS][NBIGM];                          \
  int a_shift[NBIGM];                                           \
  FTYPE a_monoindicator[NPR2INTERP][NUMMONOINDICATORS][NBIGM];  \
  int a_minorder[NBIGM];                                        \
  int a_maxorder[NBIGM];



/// NPR2INTERP always larger than NPR, so can use one memory space for both c2e and others
/// OPENMPMARK: Must define pointers and shift them within parallel region!
#define LINETYPEDEFINEPOINTERS                          \
  FTYPE (*yin)[2][NBIGM];                               \
  FTYPE (*yout)[2][NBIGM];                              \
                                                        \
  FTYPE (*shockindicator)[NBIGM];               \
  FTYPE *stiffindicator;                                \
                                                        \
  FTYPE (*yprim)[2][NBIGM];                             \
  FTYPE (*ystencilvar)[2][NBIGM];                       \
  FTYPE (*df)[NUMDFS][NBIGM];                           \
  FTYPE (*dfformono)[NUMDFS][NBIGM];                    \
  FTYPE (*drho)[NBIGM];                                 \
  FTYPE (*dP)[NBIGM];                                   \
  FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM];                     \
  FTYPE (*Pline)[NBIGM];                        \
  FTYPE (*Vline)[NBIGM];                        \
  int *shift;                                           \
  FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM];     \
  int *minorder;                                        \
  int *maxorder;                                        \
  FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM];




          
#define LINETYPEDEFINES2                                                \
  int reallim;                                                          \
  int pl,pliter;                                                        \
  int yiter;                                                            \
  int dfiter;                                                           \
  int counter;                                                          \
  int pllocal;                                                          \
  struct of_trueijkp trueijkp;                                          \
  extern int choose_limiter(int dir, int i, int j, int k, int pl);      \
  extern void compute_monotonicity_line(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]); \
  extern void compute_monotonicity_line_multipl(int stencilvarisnull, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]); \
  extern void compute_monotonicity_line_indicatoronly(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]); \
  extern void compute_monotonicity_line_valueonly(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM])

// last thing has no semicolon




/// shift pointers to account for boundary zones (similar to as in set_arrays.c)
#define LINETYPESHIFTS                                                  \
  {yin =(FTYPE (*)[2][NBIGM]) (&(a_yin[0][0][NBIGBND]));                \
    yout=(FTYPE (*)[2][NBIGM]) (&(a_yout[0][0][NBIGBND]));              \
                                                                        \
    shockindicator=(FTYPE (*)[NBIGM]) (&(a_shockindicator[0][NBIGBND]));       \
    stiffindicator=(FTYPE (*)) (&(a_stiffindicator[NBIGBND]));          \
                                                                        \
    yprim =(FTYPE (*)[2][NBIGM]) (&(a_yprim[0][0][NBIGBND]));           \
    ystencilvar =(FTYPE (*)[2][NBIGM]) (&(a_ystencilvar[0][0][NBIGBND])); \
    df=(FTYPE (*)[NUMDFS][NBIGM]) (&(a_df[0][0][NBIGBND]));             \
    dfformono=(FTYPE (*)[NUMDFS][NBIGM]) (&(a_dfformono[0][0][NBIGBND])); \
    drho=(FTYPE (*)[NBIGM]) (&(a_drho[0][NBIGBND]));                    \
    dP=(FTYPE (*)[NBIGM]) (&(a_dP[0][NBIGBND]));                        \
    etai=(FTYPE (*)[NUMTRUEEOMSETS][NBIGM]) (&(a_etai[0][0][NBIGBND]));     \
    Pline=(FTYPE (*)[NBIGM]) (&(a_Pline[0][NBIGBND]));                         \
    Vline=(FTYPE (*)[NBIGM]) (&(a_Vline[0][NBIGBND]));                         \
    shift=(int (*)) (&(a_shift[NBIGBND]));                              \
    monoindicator=(FTYPE (*)[NUMMONOINDICATORS][NBIGM]) (&(a_monoindicator[0][0][NBIGBND])); \
    minorder=(int (*)) (&(a_minorder[NBIGBND]));                        \
    maxorder=(int (*)) (&(a_maxorder[NBIGBND]));                        \
    youtpolycoef=(FTYPE (*)[MAXSPACEORDER][NBIGM]) (&(a_youtpolycoef[0][0][NBIGBND]));}



#define LINETYPEPOINTERSFUNCTIONDECLARE                 \
  FTYPE (*yin)[2][NBIGM],                               \
    FTYPE (*yout)[2][NBIGM],                            \
                                                        \
    FTYPE (*shockindicator)[NBIGM],             \
    FTYPE *stiffindicator,                              \
                                                        \
    FTYPE (*yprim)[2][NBIGM],                           \
    FTYPE (*ystencilvar)[2][NBIGM],                     \
    FTYPE (*df)[NUMDFS][NBIGM],                         \
    FTYPE (*dfformono)[NUMDFS][NBIGM],                  \
    FTYPE (*drho)[NBIGM],                               \
    FTYPE (*dP)[NBIGM],                                 \
    FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM],                               \
    FTYPE (*Pline)[NBIGM],                      \
    FTYPE (*Vline)[NBIGM],                      \
    int *shift,                                         \
    FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM],   \
    int *minorder,                                      \
    int *maxorder,                                      \
    FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]
// no comma on end



#define LINETYPEPOINTERSPASSFUNCTIONARG         \
  yin,                                          \
    yout,                                       \
                                                \
    shockindicator,                             \
    stiffindicator,                             \
                                                \
    yprim ,                                     \
    ystencilvar,                                \
    df,                                         \
    dfformono,                                  \
    drho,                                       \
    dP,                                         \
    etai,                                       \
    Pline,                                      \
    Vline,                                      \
    shift,                                      \
    monoindicator,                              \
    minorder,                                   \
    maxorder,                                   \
    youtpolycoef
// no comma on end







static void slope_lim_linetype_c2e_perijk_precomputeindicators(
                                                               LINETYPEPOINTERSFUNCTIONDECLARE,
                                                               int i, int  j, int k,
                                                               int realisinterp, int whichprimtype, int interporflux,
                                                               int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                                                               int dir, int intdir,
                                                               int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                                                               int recontype, int doingweno,
                                                               void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                                                               void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                                                               int stencilvarisnull, int preforder, int whichreduce,
                                                               int idel, int jdel, int kdel,
                                                               FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);


static void slope_lim_linetype_c2e_perijk(
                                          LINETYPEPOINTERSFUNCTIONDECLARE,
                                          int i, int  j, int k,
                                          int realisinterp, int whichprimtype, int interporflux,
                                          int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                                          int dir, int intdir,
                                          int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                                          int recontype, int doingweno,
                                          void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                                          void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                                          int stencilvarisnull, int preforder, int whichreduce,
                                          int idel, int jdel, int kdel,
                                          FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);


static void slope_lim_linetype_perijk(
                                      LINETYPEPOINTERSFUNCTIONDECLARE,
                                      int i, int  j, int k,
                                      int realisinterp, int whichprimtype, int interporflux,
                                      int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                                      int dir, int intdir,
                                      int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                                      int recontype, int doingweno,
                                      void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                                      void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                                      int stencilvarisnull, int preforder, int whichreduce,
                                      int weightsplittype, int whichquantity,
                                      int idel, int jdel, int kdel,
                                      FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);


static void assign_eno_result(int interporflux, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR], FTYPE (*result1)[NSTORE2][NSTORE3][NPR]);


/// c2e needed functions
static void get_1d_line_c2e_multipl(int whichquantity, int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*yin)[2][NBIGM], struct of_trueijkp *trueijkp);
static void get_1d_line_c2e(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE *yin, struct of_trueijkp *trueijkp);
static void get_1d_line_shockarray(int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*shockarray)[NSTORE1][NSTORE2][NSTORE3], FTYPE (*shockindicator)[NBIGM], struct of_trueijkp *trueijkp);

/// below for shock indicator used on quantities with NPR elements
static void assign_eno_result_c2e_multipl(int whichprimtype, int interporflux, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[2][NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*result1)[NSTORE2][NSTORE3][NPR2INTERP]);
static void assign_eno_result_c2e(int interporflux, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*result1)[NSTORE2][NSTORE3][NPR2INTERP]);
///below for an effective $\gamma^2$ of the flow (only hydro since does not include bsq!) SASMARK
static void get_1d_line_c2e_gammaeffhydro(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *yin, struct of_trueijkp *trueijkp);

/// interior _perijk() functions:
static void get_df_line_gen_new(int realisinterp, int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp);
static void get_df_line_paraline(int realisinterp, int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp);
static int compute_df_line(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM]);
static int compute_df_line_paraline(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM]);
static int compute_df_line_new(int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*yprim)[2][NBIGM], FTYPE *df, FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp);
static int compute_df_line_formono(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM]); 
static void get_1d_line(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*yin)[NBIGM], struct of_trueijkp *trueijkp);
static int get_V_and_P(int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp);
static int get_shock_indicator(int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], FTYPE (*shockindicator)[NBIGM], struct of_trueijkp *trueijkp);
static int get_contact_indicator(int realisinterp, int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM]);
static void causal_shift_order(int whichquantitiy, int interporflux, int dir, int preforder, int bs, int ps, int pe, int be,  int i, int j, int k,  int idel, int jdel, int kdel, int *shift, int *minorder, int *maxorder);


static void set_preforder(int dir, int interporflux, int *preforder, int*whichreduce);
static int get_recon_type(int interporflux);
static void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);
static void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp);


static void set_interp_loop(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
static void set_interp_loop_expanded(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);






/// very similar to slope_lim_linetype, but different number of quantities to handle
void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP])
{

  // VS doesn't compile unless this after all function/variable definitions
  // VS also wants no semicolons back-to-back?
  LINETYPEDEFINES1;



  // RESIDUAL NOTES:
  // DEATHMARK
  // yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
  // shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
  // primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  /////////////////
  //
  // Define loop over starting positions and range of loop for each starting position
  //
  /////////////////


  withshifts=1; // need with shifts since SUPERGENLOOP below has no shifts and shouldn't have shifts since (e.g.) for dir=1 ie=is and shifts would split the loop into 3D type loop instead of over starting positions

  int loc =CENT;
  int continuous=0;
  set_interp_loop_gen(withshifts, interporflux, dir, loc, continuous, &intdir, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk, &bs, &ps, &pe, &be);

  ///////////////////
  //
  // get reconstruction type (c2e, a2c, c2a)
  //
  ////////////////////

  recontype=get_recon_type(interporflux);


  ///////////////////
  //
  // determine which pass_1d_line and pass_1d_line_multipl to use
  //
  // assume recontype is CVT_C2E
  //
  ///////////////////

  doingweno=0;
  if(
     WENOINTERPTYPE(lim[dir])
     ){
    doingweno=1;
    pass_1d_line=&pass_1d_line_weno;
    pass_1d_line_multipl=&pass_1d_line_multipl_weno;
  }
  else if(
          lim[dir]==PARALINE
          ){
    pass_1d_line=&pass_1d_line_paraline;
    pass_1d_line_multipl=&pass_1d_line_multipl_paraline;
  }
  else{
    dualfprintf(fail_file,"No such defined pass_1d_line for lim[dir=%d]=%d\n",dir,lim[dir]);
    myexit(91758726);
  }

  ///////////////////
  //
  // set stencil variable indicator
  //
  ////////////////////

  if(stencilvar==NULL) stencilvarisnull=1;
  else stencilvarisnull=0;


  ///////////////////
  //
  // set preferred order
  //
  ////////////////////

  set_preforder(dir, interporflux, &preforder, &whichreduce);
 

  ////////////////
  //
  // LOOP OVER STARTING POSITIONS: this loop is over starting positions only
  //
  ///////////////
  //  SUPERGENLOOP(i,j,k,is,ie,js,je,ks,ke,di,dj,dk){

#if(STOREFLUXSTATE)
  // don't need to get EOS if already stored since only access already-stored data (i.e. we don't recompute state)
#pragma omp parallel 
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP
#endif
  {
    OPENMP3DLOOPVARSDEFINE;
    int i,j,k;
    LINETYPEDEFINEPOINTERS; // must be within parallel region
    LINETYPEDEFINESOTHER; // must be within parallel region
    LINETYPEDEFINESMEMORY; // must be within parallel region


    /////////////////////////////////
    // perform pointer shifts
    /////////////////////////////////
    LINETYPESHIFTS;


    OPENMP3DLOOPSETUPSUPERGEN(is,ie,js,je,ks,ke,di,dj,dk);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


#if(STORESHOCKINDICATOR)
      // avoids expensive get_V_and_P(), Ficalc, etc. since not necessary to do this for each call to this function since just diffusive-type correction.  Just ensure average (or use max or min with fabs() ) original quantity to right location when finally used.
      slope_lim_linetype_c2e_perijk_precomputeindicators(
                                                         LINETYPEPOINTERSPASSFUNCTIONARG,
                                                         i,j,k,
                                                         realisinterp, whichprimtype, interporflux,
                                                         nprlocalstart, nprlocalend, nprlocallist, numprims,
                                                         dir, intdir,
                                                         is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be,
                                                         recontype, doingweno,
                                                         pass_1d_line, pass_1d_line_multipl,
                                                         stencilvarisnull, preforder, whichreduce,
                                                         idel, jdel, kdel,
                                                         primreal, stencilvar, p2interp, pleft, pright);
#else
      slope_lim_linetype_c2e_perijk(
                                    LINETYPEPOINTERSPASSFUNCTIONARG,
                                    i,j,k,
                                    realisinterp, whichprimtype, interporflux,
                                    nprlocalstart, nprlocalend, nprlocallist, numprims,
                                    dir, intdir,
                                    is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be,
                                    recontype, doingweno,
                                    pass_1d_line, pass_1d_line_multipl,
                                    stencilvarisnull, preforder, whichreduce,
                                    idel, jdel, kdel,
                                    primreal, stencilvar, p2interp, pleft, pright);
#endif

    } // done with main loop over starting points
  }// end parallel region


}




/// interior of main SUPERGENLOOP() inside slope_lim_linetype_c2e()
void slope_lim_linetype_c2e_perijk_precomputeindicators(
                                                        LINETYPEPOINTERSFUNCTIONDECLARE,
                                                        int i, int  j, int k,
                                                        int realisinterp, int whichprimtype, int interporflux,
                                                        int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                                                        int dir, int intdir,
                                                        int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                                                        int recontype, int doingweno,
                                                        void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                                                        void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                                                        int stencilvarisnull, int preforder, int whichreduce,
                                                        int idel, int jdel, int kdel,
                                                        FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP])
{

  // VS doesn't compile unless this after all function/variable definitions
  // VS also wants no semicolons back-to-back?
  LINETYPEDEFINES2;




  

  // RESIDUAL NOTES:
  // DEATHMARK
  // yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
  // shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
  // primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.



  if(doingweno){ // GODMARK: shift,minorder,maxorder only used by weno right now
    ///////////////
    //
    // Figure out shifting of stencil and order of stencil to more closely match with causality.
    //  If no stencil shifting, then set shift and min/max order to default
    //
    ///////////////
    causal_shift_order(whichprimtype, interporflux, dir, preforder, bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, shift, minorder, maxorder);
  }

  // subloop loops from starting position to end of that line


  ///////////////////
  //
  // GET 1D LINES [assume shock indicator (or other things) already computed and exists as part of "primitives"]
  //
  ////////////////////

  // for this we don't need realisinterp and should never use yprim or primreal
  get_1d_line_c2e_multipl(ENOPRIMITIVE, dir, interporflux, bs, ps, pe, be,  i, j, k, p2interp, yin, &trueijkp);

  // get shock indicator from stored array
  if(whichreduce == WENO_REDUCE_TYPE_PPM || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ){
#if(STORESHOCKINDICATOR==0)
    dualfprintf(fail_file,"If using simple c2e method, need to set STORESHOCKINDICATOR=1\n");
    myexit(394834);
#endif

    get_1d_line_shockarray(dir, interporflux, bs, ps, pe, be,  i, j, k, GLOBALPOINT(shockindicatorarray), shockindicator,&trueijkp);
  }



#if(CONTACTINDICATOR)
  dualfprintf(fail_file,"If using simple c2e method with contacts,  need to setup flux.c's compute_and_store_fluxstatecent().\n");
  myexit(394834);
#endif



  if(PARALINEUSESMONO||doingweno){ // let paraline use mono (which needs ystencilvar)
    ///////////////////
    //
    // 1D GET LINE of stencilvar
    //
    ////////////////////
    if(!stencilvarisnull){
      get_1d_line_c2e_multipl(ENOPRIMITIVE, dir, interporflux, bs, ps, pe, be,  i, j, k, stencilvar, ystencilvar,&trueijkp);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      // notice only requesting yin[0], not yin[1], so not good for FLUXSPLIT-- GODMARK
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++) ystencilvar[pl][0][yiter]=yin[pl][0][yiter];
    }
  }


  
  if(doingweno){
    ///////////////////
    //
    // get df's for contact indicator or in general (use this to get stiffindicator too)
    //
    // gets Vline and Pline now too if doingweno==1
    ////////////////////
    get_df_line_gen_new(realisinterp, doingweno,whichprimtype,interporflux,recontype,dir,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, yin, df, &drho, &dP,stiffindicator,Vline,Pline,&trueijkp);
  }
  else{
    get_df_line_paraline(realisinterp, doingweno,whichprimtype,interporflux,recontype,dir,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, yin, df, &drho, &dP,stiffindicator,Vline,Pline,&trueijkp);
  }



  if(PARALINEUSESMONO||doingweno){ // let paraline use mono (which needs dfformono[][][])
    ///////////////////
    //
    // 1D GET LINE of mono df's if needed
    //
    // Point is that MONO (e.g. SMONO) uses df to compute if cusp exists (even if DO_SMONO_???==0).
    // Cusp sets monoindicator that sets effective stencil, so need to use stencilvar for this
    //
    ////////////////////
    if(!stencilvarisnull){
      // notice only requesting yin[0], not yin[1], so not good for FLUXSPLIT-- GODMARK
      NUMPRIMLOOP(pliter,pl) compute_df_line_formono(doingweno,interporflux, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,   yin[pl][0], dfformono[pl]);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      // should be expensive
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++) for(dfiter=0;dfiter<NUMDFS;dfiter++) dfformono[pl][dfiter][yiter]=df[pl][dfiter][yiter];
    }


    ///////////////
    //
    // Compute monotonicity indicator and if monotonic or rough then compute interface value
    // Assume only needed by WENO routines
    //
    ///////////////
      
    if(stencilvar==NULL){
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }
    else{
      // then split indicator from value
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_indicatoronly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], ystencilvar[pl][0], yout[pl], youtpolycoef[pl]);
 
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }
  }




  ///////////////
  //
  // PASS 1D LINE
  //
  ///////////////

  // assume normally want all pl's
  pass_1d_line_multipl( DO_SPLITC2E, ENOPRIMITIVE, dir, ALL_CALC, recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef,&trueijkp);

  //////////////////////////////////////////
  //
  // now have 1D line of point data corresponding to left and right interpolated values
  //
  // Assign result back to pleft/pright (or just pleft if one result)
  //
  /////////////////////////////////////////
  assign_eno_result_c2e_multipl(ENOPRIMITIVE, recontype, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout, pleft, pright);

}




/// interior of main SUPERGENLOOP() inside slope_lim_linetype_c2e()
void slope_lim_linetype_c2e_perijk(
                                   LINETYPEPOINTERSFUNCTIONDECLARE,
                                   int i, int  j, int k,
                                   int realisinterp, int whichprimtype, int interporflux,
                                   int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                                   int dir, int intdir,
                                   int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                                   int recontype, int doingweno,
                                   void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                                   void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                                   int stencilvarisnull, int preforder, int whichreduce,
                                   int idel, int jdel, int kdel,
                                   FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP])
{

  // VS doesn't compile unless this after all function/variable definitions
  // VS also wants no semicolons back-to-back?
  LINETYPEDEFINES2;




  

  // RESIDUAL NOTES:
  // DEATHMARK
  // yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
  // shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
  // primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.



  if(doingweno){ // GODMARK: shift,minorder,maxorder only used by weno right now
    ///////////////
    //
    // Figure out shifting of stencil and order of stencil to more closely match with causality.
    //  If no stencil shifting, then set shift and min/max order to default
    //
    ///////////////
    causal_shift_order(whichprimtype, interporflux, dir, preforder, bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, shift, minorder, maxorder);
  }

  // subloop loops from starting position to end of that line


  ///////////////////
  //
  // GET 1D LINES
  //
  ////////////////////

  if(whichreduce == WENO_REDUCE_TYPE_PPM || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ){
    if(realisinterp){
      // means primreal=p2interp
      // then get all primitives and assign to normal yin
      NUMPRIMLOOP(pliter,pl) get_1d_line_c2e(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, yin[pl][0],&trueijkp);
      // below means can't modify yin without having changed yprim (not general GODMARK)
      yprim=yin;
    }
    else{
      // means primreal!=p2interp
      // then need to get separate p2interp line
      get_1d_line_c2e_multipl(ENOPRIMITIVE, dir, interporflux, bs, ps, pe, be,  i, j, k, p2interp, yin, &trueijkp);

      // get prim line
      if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
        // then need to get primitives separately from interpolated quantities
        PALLREALLOOP(pl) get_1d_line_c2e(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, yprim[pl][0],&trueijkp);
      }
      else if(CONTACTINDICATOR||COMPUTEDRHODP){
        // then only need RHO and UU
        get_1d_line_c2e(dir, interporflux, RHO, bs, ps, pe, be,  i, j, k, primreal, yprim[RHO][0],&trueijkp);
        get_1d_line_c2e(dir, interporflux,  UU, bs, ps, pe, be,  i, j, k, primreal, yprim[UU][0],&trueijkp);
      }

#if(VSQ!=-100)
      //put the effective hydro value of $\gamma^2$ into the VSQ element of the primitive array that will be passed to reconstruction
      get_1d_line_c2e_gammaeffhydro(dir, interporflux, VSQ, bs, ps, pe, be,  i, j, k, idel, jdel, kdel, p2interp, yprim[VSQ][0],&trueijkp);
#endif

    }
  }
  else{
    // then no ppm reduce or contact indicator, so just get normal line
    ///////////////////
    //
    // 1D GET LINE (interpolated quantity is never the same as the prims4shock quantity since not doing c2e here)
    //
    ////////////////////
    get_1d_line_c2e_multipl(ENOPRIMITIVE, dir, interporflux, bs, ps, pe, be,  i, j, k, p2interp, yin, &trueijkp);
  }



  if(PARALINEUSESMONO||doingweno){ // let paraline use mono (which needs ystencilvar)
    ///////////////////
    //
    // 1D GET LINE of stencilvar
    //
    ////////////////////
    if(!stencilvarisnull){
      get_1d_line_c2e_multipl(ENOPRIMITIVE, dir, interporflux, bs, ps, pe, be,  i, j, k, stencilvar, ystencilvar,&trueijkp);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      // notice only requesting yin[0], not yin[1], so not good for FLUXSPLIT-- GODMARK
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++) ystencilvar[pl][0][yiter]=yin[pl][0][yiter];
    }
  }



  if(doingweno){
    ///////////////////
    //
    // get df's for contact indicator or in general (use this to get stiffindicator too)
    //
    // gets Vline and Pline now too if doingweno==1
    ////////////////////
    get_df_line_gen_new(realisinterp, doingweno,whichprimtype,interporflux,recontype,dir,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, yin, df, &drho, &dP,stiffindicator,Vline,Pline,&trueijkp);
  }
  else{
    get_df_line_paraline(realisinterp, doingweno,whichprimtype,interporflux,recontype,dir,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, yin, df, &drho, &dP,stiffindicator,Vline,Pline,&trueijkp);
  }


  if(PARALINEUSESMONO||doingweno){ // let paraline use mono (which needs dfformono[][][])
    ///////////////////
    //
    // 1D GET LINE of mono df's if needed
    //
    // Point is that MONO (e.g. SMONO) uses df to compute if cusp exists (even if DO_SMONO_???==0).
    // Cusp sets monoindicator that sets effective stencil, so need to use stencilvar for this
    //
    ////////////////////
    if(!stencilvarisnull){
      // notice only requesting yin[0], not yin[1], so not good for FLUXSPLIT-- GODMARK
      NUMPRIMLOOP(pliter,pl) compute_df_line_formono(doingweno,interporflux, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,   yin[pl][0], dfformono[pl]);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      // should be expensive
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++) for(dfiter=0;dfiter<NUMDFS;dfiter++) dfformono[pl][dfiter][yiter]=df[pl][dfiter][yiter];
    }
  }


  if(!doingweno){ // need this because get_df_line_gen_new() didn't get Vline,Pline if doingweno==0
    ///////////////////
    //
    // 1D P and V
    //
    ////////////////////
    if(SHOCKINDICATOR || CONTACTINDICATOR){
      get_V_and_P(whichprimtype, interporflux, dir, bs, ps, pe, be,  i, j, k, idel, jdel, kdel, yprim, Vline, Pline,&trueijkp);
    }
  }



  if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
    ///////////////////
    //
    // 1D GET SHOCK INDICATOR
    //
    ////////////////////
    //use primitive values that correspond to the quantities being interpolated
    get_shock_indicator(whichprimtype, interporflux, dir,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yprim, Vline, Pline, shockindicator,&trueijkp);
  }


  if(CONTACTINDICATOR){
    ///////////////////
    //
    // 1D GET CONTACT INDICATOR
    //
    ////////////////////
    //use primitive values that correspond to the quantities being interpolated
    get_contact_indicator(realisinterp, whichprimtype, interporflux, dir,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yin, Vline, Pline, etai);
  }




  if(PARALINEUSESMONO||doingweno){ // let paraline use mono
    ///////////////
    //
    // Compute monotonicity indicator and if monotonic or rough then compute interface value
    // Assume only needed by WENO routines
    //
    ///////////////
      
    if(stencilvar==NULL){
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }
    else{
      // then split indicator from value
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_indicatoronly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], ystencilvar[pl][0], yout[pl], youtpolycoef[pl]);
 
      NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, dfformono[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }
  }




  ///////////////
  //
  // PASS 1D LINE
  //
  ///////////////

  // assume normally want all pl's
  pass_1d_line_multipl( DO_SPLITC2E, ENOPRIMITIVE, dir, ALL_CALC, recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef,&trueijkp);

  //////////////////////////////////////////
  //
  // now have 1D line of point data corresponding to left and right interpolated values
  //
  // Assign result back to pleft/pright (or just pleft if one result)
  //
  /////////////////////////////////////////
  assign_eno_result_c2e_multipl(ENOPRIMITIVE, recontype, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout, pleft, pright);

}










/// linetype assumed to return pleft/pright directly since linetype's are usually higher order
/// this function used for primitive interpolation AND for flux interpolation
/// stencil var can be NULL or some pointer.  If a pointer, then should use to assign all stencil-related things
void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR])
{
  int whichprimtype=whichquantity;
  int realisinterp;
  int weightsplittype;
  // VS doesn't compile unless this after all function/variable definitions
  // VS also wants no semicolons back-to-back?
  LINETYPEDEFINES1;



  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);

  // set realisinterp *after* setting nprlists
  // below not generaly correct if doing c2e, but always doing c2e with slope_lim_linetype_c2e() anyways -- GODMARK
  if(interporflux==ENOINTERPTYPE) set_normal_realisinterp(&realisinterp);
  else realisinterp=0;





  // RESIDUAL NOTES:
  // DEATHMARK
  // yin is of size yin[bf-bs+1] +   yin[bs] is starting point for yin data.  Size:  yout[2][pf-ps+1]  yout[0/1][0] is first data point
  // shift is same size as yout.  It indicates the location of the center of the stencil w.r.t. the point of interest (i).  So a shift if 2 would mean center of stencil is at i+2.
  // primitives of some form are interpolated here and their dimensions are in standard [i][j][k] form.


  /////////////////
  //
  // Define loop over starting positions and range of loop for each starting position
  //
  /////////////////

  withshifts=1; // need with shifts since SUPERGENLOOP below has no shifts and shouldn't have shifts since (e.g.) for dir=1 ie=is and shifts would split the loop into 3D type loop instead of over starting positions

  int loc=CENT;
  int continuous=0;
  set_interp_loop_gen(withshifts,interporflux, dir, continuous, loc, &intdir, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk, &bs, &ps, &pe, &be);

  ///////////////////
  //
  // get reconstruction type (c2e, a2c, c2a)
  //
  ////////////////////

  recontype=get_recon_type(interporflux);

  ///////////////////
  //
  // determine which pass_1d_line and pass_1d_line_multipl to use
  //
  ///////////////////

  doingweno=0;
  if(
     recontype==CVT_C2E && WENOINTERPTYPE(lim[dir]) ||
     recontype==CVT_C2A && WENOINTERPTYPE(avgscheme[dir]) ||
     recontype==CVT_A2C && WENOINTERPTYPE(avgscheme[dir])
     ){
    doingweno=1;
    pass_1d_line=&pass_1d_line_weno;
    pass_1d_line_multipl=&pass_1d_line_multipl_weno;
  }
  else if(
          recontype==CVT_C2E && lim[dir]==PARALINE ||
          recontype==CVT_C2A && avgscheme[dir]==PARALINE ||
          recontype==CVT_A2C && avgscheme[dir]==PARALINE
          ){
    pass_1d_line=&pass_1d_line_paraline;
    pass_1d_line_multipl=&pass_1d_line_multipl_paraline;
  }
  else{
    dualfprintf(fail_file,"No such defined pass_1d_line for lim[dir=%d]=%d avgscheme[1,2,3]=%d %d %d\n",dir,lim[dir],avgscheme[1],avgscheme[2],avgscheme[3]);
    myexit(91758727);
  }



  ///////////////////
  //
  // set stencil variable indicator
  //
  ////////////////////

  if(stencilvar==NULL) stencilvarisnull=1;
  else  stencilvarisnull=0;


  ///////////////////
  //
  // set stencil weight reduction type (some arbitrariness for user definition)
  //
  ////////////////////

  higherorder_set(whichquantity, recontype, &weightsplittype);

  ///////////////////
  //
  // set preferred order
  //
  ////////////////////

  set_preforder(dir, interporflux, &preforder, &whichreduce);


  ////////////////
  //
  // LOOP OVER STARTING POSITIONS: this loop is over starting positions only
  //
  ///////////////
  //  SUPERGENLOOP(i,j,k,is,ie,js,je,ks,ke,di,dj,dk){

#if(STOREFLUXSTATE)
  // don't need to get EOS if already stored since only access already-stored data (i.e. we don't recompute state)
#pragma omp parallel 
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP
#endif
  {
    OPENMP3DLOOPVARSDEFINE;
    int i,j,k;
    LINETYPEDEFINEPOINTERS; // must be within parallel region
    LINETYPEDEFINESOTHER;  // must be within parallel region
    LINETYPEDEFINESMEMORY;  // must be within parallel region


    /////////////////////////////////
    // perform pointer shifts
    /////////////////////////////////
    LINETYPESHIFTS;


    OPENMP3DLOOPSETUPSUPERGEN(is,ie,js,je,ks,ke,di,dj,dk);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      slope_lim_linetype_perijk(
                                LINETYPEPOINTERSPASSFUNCTIONARG,
                                i,j,k,
                                realisinterp, whichprimtype, interporflux,
                                nprlocalstart, nprlocalend, nprlocallist, numprims,
                                dir, intdir,
                                is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be,
                                recontype, doingweno,
                                pass_1d_line, pass_1d_line_multipl,
                                stencilvarisnull, preforder, whichreduce,
                                weightsplittype, whichquantity,
                                idel, jdel, kdel,
                                primreal, stencilvar, p2interpm, p2interpp, pleft, pright);


    } // done with main loop over starting points
  }// end parallel region

}








/// inside of SUPERGENLOOP() in slope_lim_linetype()
void slope_lim_linetype_perijk(
                               LINETYPEPOINTERSFUNCTIONDECLARE,
                               int i, int  j, int k,
                               int realisinterp, int whichprimtype, int interporflux,
                               int nprlocalstart, int nprlocalend, int *nprlocallist, int numprims,
                               int dir, int intdir,
                               int is, int ie, int js, int je, int ks, int ke, int di, int dj, int dk, int bs, int ps, int pe, int be,
                               int recontype, int doingweno,
                               void (*pass_1d_line)(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp),
                               void (*pass_1d_line_multipl)(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystecilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp),
                               int stencilvarisnull, int preforder, int whichreduce,
                               int weightsplittype, int whichquantity,
                               int idel, int jdel, int kdel,
                               FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR])
{
  extern int apply_bc_line(int nprlocalstart, int nprlocalend, int*nprlocallist, int doinverse, int iter, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM],FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]);


  int monoi;
  // VS doesn't compile unless this after all function/variable definitions
  // VS also wants no semicolons back-to-back?
  LINETYPEDEFINES2;






  if(weightsplittype==CONSTANT_ALL_WEIGHTS){

    if(DOMONOINTERP==NOMONOINTERP){
      dualfprintf(fail_file,"Cannot do weightsplittype==CONSTANT_ALL_WEIGHTS if DOMONOINTERP==%d\n",DOMONOINTERP);
      myexit(2469238463);
    }

    // then just do mono with fixed weights

    // get line (don't need yprim,df,shockindicator,etc.)
    NUMPRIMLOOP(pliter,pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interpm, p2interpp, yin[pl],&trueijkp);

    // user-defined modification of the line -- usually adjusting line data so interpolation can be higher-order
    apply_bc_line(nprlocalstart,nprlocalend,nprlocallist,0,trueijkp.iter,recontype,bs,be,yin,NULL,NULL);

    NUMPRIMLOOP(pliter,pl){
      // set all weights equal
      for(monoi=bs;monoi<=be;monoi++){
        monoindicator[pl][MONOINDTYPE][monoi]=1;
        monoindicator[pl][MONOLEFTSET][monoi]=1; // or center
        monoindicator[pl][MONORIGHTSET][monoi]=1;
      }
 
      // compute MONO reconstructions given weights for all "pl"
      // note that shift, shockindicator, and df not used
      compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }





  }
  else{



    ///////////////
    //
    // Figure out shifting of stencil and order of stencil to more closely match with causality.
    //  If no stencil shifting, then set shift and min/max order to default
    //
    ///////////////
    causal_shift_order(whichquantity, interporflux, dir, preforder, bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, shift, minorder, maxorder);
  


    ///////////////////
    //
    // GET 1D LINES
    //
    ////////////////////
    // subloop loops from starting position to end of that line

    if(whichreduce == WENO_REDUCE_TYPE_PPM || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ){
      if((interporflux==ENOINTERPTYPE)&&(realisinterp)){
        // means primreal=p2interpm
        // then get all primitives and assign to normal yin
        NUMPRIMLOOP(pliter,pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yin[pl],&trueijkp);
        yprim=yin;
      }
      else{
        // means primreal!=p2interpm
        // then need to get separate p2interp line
        NUMPRIMLOOP(pliter,pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interpm, p2interpp, yin[pl],&trueijkp);

        // get prim line
        if(whichreduce == WENO_REDUCE_TYPE_PPM || SHOCKINDICATOR ){
          // then need to get primitives separately from interpolated quantities
          PALLREALLOOP(pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[pl],&trueijkp);
        }
        else if(CONTACTINDICATOR||COMPUTEDRHODP){
          // then only need RHO and UU
          get_1d_line(dir, interporflux, RHO, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[RHO],&trueijkp);
          get_1d_line(dir, interporflux,  UU, bs, ps, pe, be,  i, j, k, primreal, NULL, yprim[UU],&trueijkp);
        }
      }
    }
    else{
      // then no ppm reduce or contact indicator
      ///////////////////
      //
      // 1D GET LINE (interpolated quantity is never the same as the prims4shock quantity since not doing c2e here)
      //
      ////////////////////
      NUMPRIMLOOP(pliter,pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, p2interpm, p2interpp, yin[pl],&trueijkp);
    }



    ///////////////////
    //
    // 1D GET LINE of stencilvar
    //
    ////////////////////
    if(!stencilvarisnull){
      NUMPRIMLOOP(pliter,pl) get_1d_line(dir, interporflux, pl, bs, ps, pe, be,  i, j, k, stencilvar, stencilvar, ystencilvar[pl],&trueijkp);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++){
        ystencilvar[pl][0][yiter]=yin[pl][0][yiter];
        ystencilvar[pl][1][yiter]=yin[pl][1][yiter];
      }
    }
    


    // user-defined modification of the line -- usually adjusting line data so interpolation can be higher-order
    apply_bc_line(nprlocalstart,nprlocalend,nprlocallist,0,trueijkp.iter,recontype,bs,be,yin,NULL,NULL);
  

    ///////////////////
    //
    // get df's for contact indicator or in general
    // GODMARK: not done for flux splitting method (i.e. assume only 1 input always)
    // gets Vline and Pline now too if doingweno==1
    //
    ////////////////////
    get_df_line_gen_new(realisinterp, doingweno,whichprimtype,interporflux,recontype,dir,whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, yin, df, &drho, &dP,stiffindicator,Vline,Pline,&trueijkp);



    ///////////////////
    //
    // 1D GET LINE of mono df's if needed
    //
    ////////////////////
    if(!stencilvarisnull){
      // notice only requesting yin[0], not yin[1], so not good for FLUXSPLIT-- GODMARK
      NUMPRIMLOOP(pliter,pl) compute_df_line_formono(doingweno,interporflux, recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,  ystencilvar[pl][0], dfformono[pl]);
    }
    else{
      //  assign rather than pointer assign since yin//ystencilvar may be modified and don't want to have to worry if was
      // should be expensive
      NUMPRIMLOOP(pliter,pl) for(yiter=bs;yiter<=be;yiter++) for(dfiter=0;dfiter<NUMDFS;dfiter++) dfformono[pl][dfiter][yiter]=df[pl][dfiter][yiter];
    }

    

    ///////////////////
    //
    // 1D P and V
    //
    ////////////////////
    if(!doingweno){ // need this because get_df_line_gen_new() didn't get Vline,Pline if doingweno==0
      if(SHOCKINDICATOR || CONTACTINDICATOR){
        get_V_and_P(whichprimtype, interporflux, dir, bs, ps, pe, be,  i, j, k, idel, jdel, kdel, yprim, Vline, Pline,&trueijkp);
      }
    }

    ///////////////////
    //
    // 1D GET SHOCK INDICATOR
    //
    ////////////////////

    if(whichreduce == WENO_REDUCE_TYPE_PPM|| SHOCKINDICATOR ){
      //use primitive values that correspond to the quantities being interpolated
      get_shock_indicator(whichprimtype, interporflux, dir,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yprim, Vline, Pline, shockindicator,&trueijkp);
    }

    ///////////////////
    //
    // 1D GET CONTACT INDICATOR
    //
    ////////////////////
    if(CONTACTINDICATOR){
      //use primitive values that correspond to the quantities being interpolated
      get_contact_indicator(realisinterp, whichprimtype, interporflux, dir,  bs, ps, pe, be,  i, j, k, idel, jdel ,kdel, yin, Vline, Pline, etai);
    }




    ///////////////
    //
    // Compute monotonicity indicator and if monotonic or rough then compute interface value
    // GODMARK: not done for flux splitting method (i.e. assume only 1 input always)
    //
    ///////////////

    if(PARALINEUSESMONO||doingweno){ // let paraline use mono
      // note dfformono used, not df
      compute_monotonicity_line_multipl(stencilvarisnull,weightsplittype,whichquantity,dir,ALL_CALC,recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, dfformono, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef);
    }




#if(0)
    
    if(1||crapdebug){
      NUMPRIMLOOP(pliter,pl){
        if(recontype==CVT_C2E && (monoindicator[pl][0][0]!=monoindicator[pl][0][0+N1*idel+N2*jdel+N3*kdel] || monoindicator[pl][1][0]!=monoindicator[pl][1][0+N1*idel+N2*jdel+N3*kdel] || monoindicator[pl][2][0]!=monoindicator[pl][2][0+N1*idel+N2*jdel+N3*kdel])
           ||
           (monoindicator[pl][0][0]!=monoindicator[pl][0][0+N1*idel+N2*jdel+N3*kdel] || monoindicator[pl][1][0]!=monoindicator[pl][1][0+N1*idel+N2*jdel+N3*kdel])
           ){
          dualfprintf(fail_file,"Detected asymmetry: recontype=%d dir=%d ps=%d pe=%d\n",recontype,dir,ps,pe);
          dualfprintf(fail_file,"i=%d j=%d k=%d di=%d dj=%d dk=%d\n",i,j,k,di,dj,dk);
          dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",monoindicator[pl][0][0],monoindicator[pl][0][0+N1*idel+N2*jdel+N3*kdel], monoindicator[pl][1][0],monoindicator[pl][1][0+N1*idel+N2*jdel+N3*kdel],monoindicator[pl][2][0],monoindicator[pl][2][0+N1*idel+N2*jdel+N3*kdel]);
        }
      }
    }
#endif






    ///////////////
    //
    // PASS 1D LINE
    //
    ///////////////

    // assume default behavior is to send all pl's
    pass_1d_line_multipl(weightsplittype, whichquantity, dir, ALL_CALC, recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef,&trueijkp);


    // user-defined de-modification of the line -- usually adjusting line data so interpolation can be higher-order
    apply_bc_line(nprlocalstart,nprlocalend,nprlocallist,1,trueijkp.iter,recontype,bs,be,yin,yout,youtpolycoef);


  }// end if doing normal variation of weights



  

  //////////////////////////////////////////
  //
  // now have 1D line of point data corresponding to left and right interpolated values
  //
  // Assign result back to pleft/pright (or just pleft if one result)
  //
  /////////////////////////////////////////
  NUMPRIMLOOP(pliter,pl){
#if(0)
    if(crapdebug){
      dualfprintf(fail_file,"pl=%d\n",pl);
    }
    //      if(crapdebug==0 && pl!=3) assign_eno_result(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);
    //assign_eno_result(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);
    if(crapdebug==0 || crapdebug==1 || pl!=3 ) assign_eno_result(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);
    else dualfprintf(fail_file,"crapdebug=%d\n",crapdebug);
#endif

    assign_eno_result(recontype, pl, bs, ps, pe, be, i, j, k, idel, jdel, kdel, yout[pl], pleft, pright);
      
  }



}










/// get type of reconstruction to perform
static int get_recon_type(int interporflux)
{

  // c2a stuff
  if(interporflux==ENOFLUXAVG1TYPE|| interporflux==ENOFLUXAVG2TYPE || interporflux==ENOFLUXAVG3TYPE || interporflux==ENOCENT2AVGTYPE || interporflux==ENOQUASIFIELDFLUXRECONTYPE) return(CVT_C2A);
  // c2e stuff
  else if(interporflux==ENOFLUXSPLITTYPE || interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF) return(CVT_C2E);
  // a2c stuff
  else if(interporflux==ENOFLUXRECONTYPE || interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOAVG2CENTTYPE) return(CVT_A2C);
  else{
    dualfprintf(fail_file,"No such interporflux=%d\n", interporflux);
    myexit(177);
  }

  return(-100);

}



/// sets preferred order based upon scheme and size of stencil
static void set_preforder(int dir, int interporflux, int *preforder, int*whichreduce)
{

  if(interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF){
    //      reallim=choose_limiter(dir,i,j,k,pl); // starting point chooses limiter type
    *preforder=interporder[lim[dir]]; // get order of scheme

    if( lim[dir] == WENO5FLAT ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
      *whichreduce = WENO_REDUCE_TYPE_PPM;
      *preforder = 5;
    }
    else if( WENOBNDPINTERPTYPE(lim[dir]) ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
      *whichreduce = WENO_REDUCE_TYPE_DEFAULT;
      *preforder = 5;  //correct the order of the scheme because the number of points passed to it is larger than its order (the order is 5)
    }
    else if( lim[dir] == PARALINE ) {  // true order
      *whichreduce = WENO_REDUCE_TYPE_DEFAULT;
      *preforder = 5;  //correct the order of the scheme because the number of points passed to it is larger than its order (the order is 5)
    }
    else {
      *whichreduce = WENO_REDUCE_TYPE_DEFAULT;
    }
  }
  else{
    *preforder=interporder[avgscheme[dir]];

    if( avgscheme[dir] == WENO5FLAT ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
      *whichreduce = WENO_REDUCE_TYPE_PPM;
      *preforder -= 2;
    }
    else if( WENOBNDPINTERPTYPE(avgscheme[dir]) ) {  //correct the order for WENO5FLAT: in this case need more points than the order for stencil reduction
      *whichreduce = WENO_REDUCE_TYPE_DEFAULT;
      *preforder = 5;  //correct the order of the scheme because the number of points passed to it is larger than its order (the order is 5)
    }
    else {
      *whichreduce = WENO_REDUCE_TYPE_DEFAULT;
    }
  }


}








/// wrapper for compute_df_line and can be used by both c2e and other methods
/// whichprimtype==ENOPRIMITIVE -> NPR2INTERP type else NPR types
static void get_df_line_gen(int realisinterp, int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM], FTYPE *stiffindicator)
{
  int pl,pliter;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);



  if(CONTACTINDICATOR||COMPUTEDRHODP){
    if((interporflux==ENOINTERPTYPE)&&(realisinterp)){
      // then don't need to get separate drho and dP

      // then get all df's
      NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);

      // then primitives are same as interpolated quantities
      // assign drho and dP.  This overwrites original pointer that pointed to independent memory for drho and dP
      *drhoptr=df[RHO];
      *dPptr=df[UU];
    }
    else{
      // then need to compute separate drho and dP
      ///////////////
      // Compute differentials needed for contactindicator
      ///////////////
      compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, RHO, bs, ps, pe, be, minorder, maxorder, shift, yprim[RHO][0], *drhoptr);
      compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, UU, bs, ps, pe, be, minorder, maxorder, shift, yprim[UU][0], *dPptr);

      // then get all p2interp df's
      NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);
      // if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][1], dfp[pl]);
    }
  }
  else{
    // then don't need drho or dP at all, just get p2interp df's
    // then get all p2interp df's
    NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);
    // if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype, whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][1], dfp[pl]);
  }
}







/// wrapper for compute_df_line and can be used by both c2e and other methods
static void get_df_line_gen_new(int realisinterp, int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp)
{
  int pl,pliter;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;




  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);





  if(SHOCKINDICATOR || CONTACTINDICATOR||COMPUTEDRHODP){
    if((interporflux==ENOINTERPTYPE)&&(realisinterp)){
      // then don't need to get separate drho and dP


      // then get all df's
      NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);

      // weno uses Lorentz factor-based differences (gets stiffindicator as well)
      if(doingweno) compute_df_line_new(doingweno, whichprimtype, interporflux, recontype, dir, whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, (*dPptr)[0], stiffindicator, Vline, Pline,trueijkp);

      // then primitives are same as interpolated quantities
      // assign drho and dP.  This overwrites original pointer that pointed to independent memory for drho and dP
      *drhoptr=df[RHO];
      //      *dPptr=df[UU];
    }
    else{
      // then need to compute separate drho and dP
      ///////////////
      // Compute differentials needed for contactindicator
      ///////////////

      // then get all p2interp df's
      NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);
      // if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][1], dfp[pl]);

      compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, RHO, bs, ps, pe, be, minorder, maxorder, shift, yprim[RHO][0], *drhoptr);

      if(doingweno) compute_df_line_new(doingweno,whichprimtype,interporflux,recontype,dir, whichreduce,preforder, bs, ps, pe, be, minorder, maxorder, shift, yprim, (*dPptr)[0], stiffindicator, Vline, Pline,trueijkp);

    }
  }
  else{
    // then don't need drho or dP at all, just get p2interp df's
    // then get all p2interp df's
    NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);
    // if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][1], dfp[pl]);
  }
}



/// wrapper for compute_df_line and can be used by both c2e and other methods
static void get_df_line_paraline(int realisinterp, int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*yprim)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (**drhoptr)[NBIGM], FTYPE (**dPptr)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp)
{
  int pl,pliter;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;




  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);



  // then don't need drho or dP at all, just get p2interp df's
  // then get all p2interp df's
  NUMPRIMLOOP(pliter,pl) compute_df_line_paraline(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][0], df[pl]);
  // if(interporflux==ENOFLUXSPLITTYPE) NUMPRIMLOOP(pliter,pl) compute_df_line(doingweno,interporflux,recontype,whichreduce,preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, yin[pl][1], dfp[pl]);


}





/// just a wrapper
void set_interp_loop_gen(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{

  if(useghostplusactive){
    // now assume active + active/ghost + ghost layers so only bound primitives during calculation
    set_interp_loop_expanded(withshifts, interporflux, dir, loc, continuous, intdir, is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be);
  }
  else{
    // straight-forward average or de-average along dir and just one ghost layer (original method)
    // GODMARK: later should convert ENOFLUXRECON method to have expanded ghost layer and ghost+active layer
    set_interp_loop(withshifts, interporflux, dir, loc, continuous, intdir, is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be);
  }

}



/// Setup loop over *starting positions* that and define the line of data corresponding to data needed by the given (dir && interporflux) scenario
/// i.e. Define loop over starting positions and range of loop for each starting position
/// Note, function should still return result that is within bounds even for unused directions (e.g. dir==3 should not go out of bounds even if N3==1)
/// Note that interp loop gives output from i=0..N so consistent with requirements for FLUXBSTAG and IF3DSPCTHENMPITRANSFERATPOLE
static void set_interp_loop(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{

  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 1
  //
  //
  ////////////////////
  // determine range of outer loop and range to feed to eno scheme
  if(dir==1){
    *intdir=dir;

    *is=*ie=0; // anything so di=1 iterates, next subloop overwrites i actually useds

    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF)){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL1;
      *be=OUTFULL1;

      *ps=0    -SHIFT1;
      *pe=N1-1 +SHIFT1;
    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL1;
      *be=OUTFULL1;

      *ps=0;
      *pe=N1-1+SHIFT1;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }

      
    }

    *js=fluxloop[dir][FJS]; // 0
    *je=fluxloop[dir][FJE]; // N2-1;

    *ks=fluxloop[dir][FKS]; //0;
    *ke=fluxloop[dir][FKE]; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *bs += SHIFTX1DN;
      *be += SHIFTX1UP;
      *ps += SHIFTX1DN;
      *pe += SHIFTX1UP;
      *js += SHIFTX2DN;
      *je += SHIFTX2UP;
      *ks += SHIFTX3DN;
      *ke += SHIFTX3UP;
    }

  }
  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 2
  //
  //
  ////////////////////
  else if(dir==2){
    *intdir=dir;

    *is=fluxloop[dir][FIS]; //0;
    *ie=fluxloop[dir][FIE]; //N1-1;


    *js=*je=0; // anything so dj=1 iterates, next subloop overwrites j actually useds

    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE) || interporflux==ENOINTERPTYPE4EMF){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      *ps=0    -SHIFT2;
      *pe=N2-1 +SHIFT2;
    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      *ps=0;
      *pe=N2-1+SHIFT2;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }

    }

    *ks=fluxloop[dir][FKS]; // 0;
    *ke=fluxloop[dir][FKE]; // N3-1;

    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *is += SHIFTX1DN;
      *ie += SHIFTX1UP;
      *bs += SHIFTX2DN;
      *be += SHIFTX2UP;
      *ps += SHIFTX2DN;
      *pe += SHIFTX2UP;
      *ks += SHIFTX3DN;
      *ke += SHIFTX3UP;
    }

  }
  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 3
  //
  //
  ////////////////////
  else if(dir==3){
    *intdir=dir;

    *is=fluxloop[dir][FIS]; // 0;
    *ie=fluxloop[dir][FIE]; // N1-1;

    *js=fluxloop[dir][FJS]; // 0;
    *je=fluxloop[dir][FJE]; // N2-1;

    *ks=*ke=0; // anything so dk=1 iterates, next subloop overwrites k actually used


    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE) || interporflux==ENOINTERPTYPE4EMF){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      *ps=0    -SHIFT3;
      *pe=N3-1 +SHIFT3;
    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      *ps=0;
      *pe=N3-1+SHIFT3;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }

    }

    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *is += SHIFTX1DN;
      *ie += SHIFTX1UP;
      *js += SHIFTX2DN;
      *je += SHIFTX2UP;
      *bs += SHIFTX3DN;
      *be += SHIFTX3UP;
      *ps += SHIFTX3DN;
      *pe += SHIFTX3UP;
    }


  }

}





/// Setup loop over *starting positions* that and define the line of data corresponding to data needed by the given (dir && interporflux) scenario
/// i.e. Define loop over starting positions and range of loop for each starting position
/// This function is for any method using the expanded ghost+ghost/active+active layers
/// Note, function should still return result that is within bounds even for unused directions (e.g. dir==3 should not go out of bounds even if N3==1)
static void set_interp_loop_expanded(int withshifts, int interporflux, int dir, int loc, int continuous, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{
  int dir_exception[NDIM];
  int *myUconsloop;


  if(interporflux==ENOINTERPTYPE4EMF){
    myUconsloop=emfUconsloop;
  }
  else{
    myUconsloop=Uconsloop;
  }
  
  //  if((interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE) ){


  // ENOFLUXAVG?TYPE is ?=direction of integration not direction of flux
  if(
     ((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
     ((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
     ((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
     ){
    dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
    myexit(1);
  }

  dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

  // determine range of outer loop and range to feed to eno scheme
  // Effectively below corresponds to 3 cases:
  // For example, for dir==1:
  // 1) not an AVG?TYPE : enter if dir==1
  // 2) AVG<dir>type : enter
  // 3) AVG<other dir>type : do not enter

  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 1
  //
  //
  ////////////////////
  if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){
    *intdir=1; // not generally true that intdir=dir

    *is=*ie=0; // anything so di=1 iterates, next subloop overwrites i actually used

    // input and output at different location

    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE) || interporflux==ENOINTERPTYPE4EMF){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL1;
      *be=OUTFULL1;

      *ps=myUconsloop[FIS]-SHIFT1;
      *pe=myUconsloop[FIE]+SHIFT1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];
    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // input and output at same location

      // then edge quantity and need from 0 .. N
      *bs=INFULL1;
      *be=OUTFULL1;

      //      *ps=myUconsloop[FIS];
      //      *pe=myUconsloop[FIE]+1;
      *ps=0;
      *pe=N1-1+SHIFT1;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }


      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];

    }
    else if( interporflux==ENOAVG2CENTTYPE ){
      // input and output at same location

      // then edge quantity and need from 0 .. N
      *bs=myUconsloop[FIS];  
      *be=myUconsloop[FIE];

      *ps=0;
      *pe=N1-1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];


    }
    else if( interporflux==ENOCENT2AVGTYPE ){
      *bs=INFULL1;
      *be=OUTFULL1;

      *ps=myUconsloop[FIS];
      *pe=myUconsloop[FIE];

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];

    }
    else if(interporflux==ENOFLUXAVG1TYPE){

      // then edge quantity and need from 0 .. N and have ghost zones out to ijkminmax[1][0]+1
      *bs=INFULL1;
      *be=OUTFULL1;

      *ps=myUconsloop[FIS];
      *pe=myUconsloop[FIE];

      if(dir==2){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *js=myUconsloop[FJS];
        *je=myUconsloop[FJE]+SHIFT2;
      }
      else{
        // GODMARK: needed?
        *js=myUconsloop[FJS];
        *je=myUconsloop[FJE];
        // *js=fluxloop[dir][FJS];
        // *je=fluxloop[dir][FJE];
      }

      if(dir==3){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *ks=myUconsloop[FKS];
        *ke=myUconsloop[FKE]+SHIFT3;
      }
      else{
        // GODMARK: needed?
        *ks=myUconsloop[FKS];
        *ke=myUconsloop[FKE];
        //*ks=fluxloop[dir][FKS];
        // *ke=fluxloop[dir][FKE];
      }

    }


    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *bs += SHIFTX1DN;
      *be += SHIFTX1UP;
      *ps += SHIFTX1DN;
      *pe += SHIFTX1UP;
      *js += SHIFTX2DN;
      *je += SHIFTX2UP;
      *ks += SHIFTX3DN;
      *ke += SHIFTX3UP;
    }

  }
  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 2
  //
  //
  ////////////////////
  else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){
    *intdir=2; // not generally true that intdir=dir

    *js=*je=0; // anything so dj=1 iterates, next subloop overwrites j actually used

    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE) || interporflux==ENOINTERPTYPE4EMF){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      *ps=myUconsloop[FJS]-SHIFT2;
      *pe=myUconsloop[FJE]+SHIFT2;

      *is=INFULL1;
      *ie=OUTFULL1;

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];
    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      //      *ps=myUconsloop[FJS];
      //      *pe=myUconsloop[FJE]+1;
      *ps=0;
      *pe=N2-1+SHIFT2;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }


      *is=INFULL1;
      *ie=OUTFULL1;

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];
    }
    else if( interporflux==ENOAVG2CENTTYPE ){
      // input and output at same location

      // then edge quantity and need from 0 .. N
      *bs=myUconsloop[FJS];  //atch correct  SASMARK; ALSO need to set up the same thing for other dimensions-- not sure if this is enough
      *be=myUconsloop[FJE];

      *ps=0;
      *pe=N2 - 1;

      *is=INFULL1;
      *ie=OUTFULL1;

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];
    }
    else if(interporflux==ENOCENT2AVGTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      *ps=myUconsloop[FJS];
      *pe=myUconsloop[FJE];

      *is=INFULL1;
      *ie=OUTFULL1;

      *ks=fluxloop[dir][FKS];
      *ke=fluxloop[dir][FKE];
    }
    else if(interporflux==ENOFLUXAVG2TYPE){
      //is bs and be to be corrected for the shock indicator? SASMARK
      // then edge quantity and need from 0 .. N
      *bs=INFULL2;
      *be=OUTFULL2;

      *ps=myUconsloop[FJS];
      *pe=myUconsloop[FJE];

      if(dir==1){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *is=myUconsloop[FIS];
        *ie=myUconsloop[FIE]+SHIFT1;
      }
      else{
        *is=myUconsloop[FIS];
        *ie=myUconsloop[FIE];
        // GODMARK: needed?
        // *is=INFULL1;
        // *ie=OUTFULL1;
      }

      if(dir==3){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *ks=myUconsloop[FKS];
        *ke=myUconsloop[FKE]+SHIFT3;
      }
      else{
        // GODMARK: needed?
        *ks=myUconsloop[FKS];
        *ke=myUconsloop[FKE];
        // *ks=fluxloop[dir][FKS];
        // *ke=fluxloop[dir][FKE];
      }

    }


    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *is += SHIFTX1DN;
      *ie += SHIFTX1UP;
      *bs += SHIFTX2DN;
      *be += SHIFTX2UP;
      *ps += SHIFTX2DN;
      *pe += SHIFTX2UP;
      *ks += SHIFTX3DN;
      *ke += SHIFTX3UP;
    }

  }
  /////////////////////
  //
  //
  // DIRECTION OF LINE EXTRACTION and INTERPOLATION = 3
  //
  //
  ////////////////////
  else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){
    *intdir=3; // not generally true that intdir=dir


    *ks=*ke=0; // anything so dk=1 iterates, next subloop overwrites k actually used


    if((interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOINTERPTYPE) || interporflux==ENOINTERPTYPE4EMF){
      // then centered quantity, and need from -1 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      *ps=myUconsloop[FKS]-SHIFT3;
      *pe=myUconsloop[FKE]+SHIFT3;

      *is=INFULL1;
      *ie=OUTFULL1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

    }
    else if(interporflux==ENOFLUXRECONTYPE ||  interporflux==ENOFLUXRECONTYPEGHOSTACTIVE || interporflux==ENOQUASIFIELDFLUXRECONTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      //      *ps=myUconsloop[FKS];
      //      *pe=myUconsloop[FKE]+1;
      *ps=0;
      *pe=N3-1+SHIFT3;

      if(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE){
        *ps -= (interporder[avgscheme[dir]]-1)/2;
        *pe += (interporder[avgscheme[dir]]-1)/2;
      }


      *is=INFULL1;
      *ie=OUTFULL1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

    }
    else if( interporflux==ENOAVG2CENTTYPE ){
      // input and output at same location

      // then edge quantity and need from 0 .. N
      *bs=myUconsloop[FKS];  //atch correct  SASMARK; ALSO need to set up the same thing for other dimensions-- not sure if this is enough
      *be=myUconsloop[FKE];

      *ps=0;
      *pe=N3 - 1;

      *is=INFULL1;
      *ie=OUTFULL1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];
    }
    else if(interporflux==ENOCENT2AVGTYPE){
      // then edge quantity and need from 0 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      *ps=myUconsloop[FKS];
      *pe=myUconsloop[FKE];

      *is=INFULL1;
      *ie=OUTFULL1;

      *js=fluxloop[dir][FJS];
      *je=fluxloop[dir][FJE];

    }
    else if(interporflux==ENOFLUXAVG3TYPE){

      // then edge quantity and need from 0 .. N
      *bs=INFULL3;
      *be=OUTFULL3;

      *ps=myUconsloop[FKS];
      *pe=myUconsloop[FKE];

      if(dir==1){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *is=myUconsloop[FIS];
        *ie=myUconsloop[FIE]+SHIFT1;
      }
      else{
        // GODMARK: needed?
        *is=myUconsloop[FIS];
        *ie=myUconsloop[FIE];
        // *is=INFULL1;
        // *ie=OUTFULL1;
      }

      if(dir==2){
        // This method doesn't need "boundary" direction of "fluxes" since they don't actually exist
        *js=myUconsloop[FJS];
        *je=myUconsloop[FJE]+SHIFT2;
      }
      else{
        *js=myUconsloop[FJS];
        *je=myUconsloop[FJE];
        // GODMARK: needed?
        // *js=fluxloop[dir][FJS];
        // *je=fluxloop[dir][FJE];
      }

    }

    *di=1;
    *dj=1;
    *dk=1;

    if(withshifts){
      // shift due to grid sectionion SECTIONMARK
      *is += SHIFTX1DN;
      *ie += SHIFTX1UP;
      *js += SHIFTX2DN;
      *je += SHIFTX2UP;
      *bs += SHIFTX3DN;
      *be += SHIFTX3UP;
      *ps += SHIFTX3DN;
      *pe += SHIFTX3UP;
    }

  }


  //  dualfprintf(fail_file,"pspe: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",dir, interporflux, *is, *ie, *js, *je, *ks, *ke, *di, *dj, *dk, *bs, *ps, *pe, *be);


}





/// Get 1D line of data so can pass it to ENO scheme (used for c2e only)
static void get_1d_line_c2e_multipl(int whichquantity, int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*yin)[2][NBIGM], struct of_trueijkp *trueijkp)
{
  // for NUMPRIMLOOP:
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  // others:
  int yiniter;
  int di2,dj2,dk2;
  int i2,j2,k2;
  int is2,ie2,js2,je2,ks2,ke2;
  int dir_exception[NDIM];



  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);



  trueijkp->dir=dir;
  trueijkp->iter=dir;
  trueijkp->interporflux=interporflux;


  // determine range of outer loop and range to feed to eno scheme
  if(dir==1){

    trueijkp->i=0;
    trueijkp->p=CENT;
    is2=bs;
    ie2=be;

    trueijkp->j=js2=je2=j;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;

  }
  else if(dir==2){

    trueijkp->i=is2=ie2=i;
    
    trueijkp->j=0;
    trueijkp->p=CENT;
    js2=bs;
    je2=be;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else if(dir==3){

    trueijkp->i=is2=ie2=i;

    trueijkp->j=js2=je2=j;

    trueijkp->k=0;
    trueijkp->p=CENT;
    ks2=bs;
    ke2=be;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else{
    dualfprintf(fail_file,"No such dir=%d in get_1d_line_c2e_multipl()\n",dir);
    myexit(246344576);
  }


  if(DOENOFLUX != NOENOFLUX){
    dualfprintf(fail_file,"Reinvestigate why ENO needs this\n");
    myexit(9586);
    // JCM: expensive and not sure why doing it
    // reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
    pl=0; // fake
    for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
      yin[pl][0][yiniter] = 0;
    }
  }



  // 1 input, assumed interporflux==ENOINTERPTYPE
  yiniter=bs;
  SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){
    NUMPRIMLOOP(pliter,pl){ // as with assign_eno_result_c2e_multipl(), more cache friendly to have pl loop inside sinde p2interp has otherwise larger stride to reach other pl's compared to yin
      yin[pl][0][yiniter]=MACP0A1(p2interp,i2,j2,k2,pl);
    }
    yiniter++;
  }


  if(! (interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF) ){
    dualfprintf(fail_file,"get_1d_line_c2e_multipl only handles interporflux==ENOINTERPTYPE or ENOINTERPTYPE4EMF\n");
    myexit(25);
  }



}


/// Get 1D line of data so can pass it to ENO scheme (used for c2e only)
static void get_1d_line_c2e(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *yin, struct of_trueijkp *trueijkp)
{
  int yiniter;
  int di2,dj2,dk2;
  int i2,j2,k2;
  int is2,ie2,js2,je2,ks2,ke2;
  int dir_exception[NDIM];


  trueijkp->dir=dir;
  trueijkp->iter=dir;
  trueijkp->interporflux=interporflux;


  // determine range of outer loop and range to feed to eno scheme
  if(dir==1){

    trueijkp->i=0;
    trueijkp->p=CENT;
    is2=bs;
    ie2=be;

    trueijkp->j=js2=je2=j;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;

  }
  else if(dir==2){

    trueijkp->i=is2=ie2=i;

    trueijkp->j=0;
    trueijkp->p=CENT;
    js2=bs;
    je2=be;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else if(dir==3){

    trueijkp->i=is2=ie2=i;

    trueijkp->j=js2=je2=j;

    trueijkp->k=0;
    trueijkp->p=CENT;
    ks2=bs;
    ke2=be;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else{
    dualfprintf(fail_file,"No such dir=%d in get_1d_line_c2e()\n",dir);
    myexit(246344576);
  }


  if(DOENOFLUX != NOENOFLUX){
    dualfprintf(fail_file,"Reinvestigate why ENO needs this\n");
    myexit(9587);
    // JCM: expensive and not sure why doing it
    // reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
    for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
      yin[yiniter] = 0;
    }
  }


  // 1 input, assumed interporflux==ENOINTERPTYPE
  yiniter=bs;
  SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

    yin[yiniter]=MACP0A1(p2interp,i2,j2,k2,pl);
    yiniter++;
  }

  if(! (interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF) ){
    dualfprintf(fail_file,"get_1d_line_c2e only handles interporflux==ENOINTERPTYPE or ENOINTERPTYPE4EMF\n");
    myexit(25);
  }



}

/// Get 1D line of data so can pass it to ENO scheme (used for c2e only)
void get_1d_line_c2e_gammaeffhydro(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *yin, struct of_trueijkp *trueijkp)
{
  int num,locali,localj,localk,localloc;
  int yiniter;
  int di2,dj2,dk2;
  int i2,j2,k2;
  int is2,ie2,js2,je2,ks2,ke2;
  int dir_exception[NDIM];
  FTYPE rho, u, p, gamma;


  //SASMARK: not sure if this is req'd
  //trueijkp->dir=dir;
  //trueijkp->iter=dir;
  //trueijkp->interporflux=interporflux;


  // determine range of outer loop and range to feed to eno scheme
  if(dir==1){

    trueijkp->i=0;
    trueijkp->p=CENT;
    is2=bs;
    ie2=be;

    trueijkp->j=js2=je2=j;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;

  }
  else if(dir==2){

    trueijkp->i=is2=ie2=i;

    trueijkp->j=0;
    trueijkp->p=CENT;
    js2=bs;
    je2=be;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else if(dir==3){

    trueijkp->i=is2=ie2=i;

    trueijkp->j=js2=je2=j;

    trueijkp->k=0;
    trueijkp->p=CENT;
    ks2=bs;
    ke2=be;

    di2=1;
    dj2=1;
    dk2=1;
  }

  if(DOENOFLUX != NOENOFLUX){
    dualfprintf(fail_file,"Reinvestigate why ENO needs this\n");
    myexit(9588);
    // JCM: expensive and not sure why doing it
    // reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
    for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
      yin[yiniter] = 0;
    }
  }


  // 1 input, assumed interporflux==ENOINTERPTYPE
  yiniter=bs;
  SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

    num=i2*idel+j2*jdel+k2*kdel;

    locali=trueijkp->i+di2*num;
    localj=trueijkp->j+dj2*num;
    localk=trueijkp->k+dk2*num;
    localloc=trueijkp->p;

    gamma = MACP0A1(p2interp,i2,j2,k2,VSQ);  //Lorentz factor
    rho = MACP0A1(p2interp,i2,j2,k2,RHO);    //comoving mass density
    u = MACP0A1(p2interp,i2,j2,k2,UU);       //comoving internal energy
    p = pressure_rho0_u_simple(locali,localj,localk,localloc, rho, u );      //comoving plasma pressure

    //an indicator that is giving an effective $\gamma^2$ for a flow
    //in a shock where the bulk matter comes to rest, $p = \gamma^2 \rho$, hence $p/\rho$ is an estimate of the square of the $\gamma$-factor
    yin[yiniter]= gamma * ( fabs(rho) + fabs(u) + fabs(p) ) / ( fabs(rho) + SQRTMINNUMREPRESENT );
    yiniter++;
  }

  if(! (interporflux==ENOINTERPTYPE || interporflux==ENOINTERPTYPE4EMF) ){
    dualfprintf(fail_file,"get_1d_line_c2e only handles interporflux==ENOINTERPTYPE or ENOINTERPTYPE4EMF\n");
    myexit(25);
  }



}


/// Get 1D line of data so can pass it to ENO scheme
static void get_1d_line(int dir, int interporflux, int pl, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*yin)[NBIGM], struct of_trueijkp *trueijkp)
{
  int yiniter;
  int di2,dj2,dk2;
  int i2,j2,k2;
  int is2,ie2,js2,je2,ks2,ke2;
  int dir_exception[NDIM];


  if(
     ((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
     ((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
     ((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
     ){
    dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
    myexit(1);
  }

  dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

  trueijkp->dir=dir;
  trueijkp->interporflux=interporflux;

  // determine range of outer loop and range to feed to eno scheme
  if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){

    trueijkp->iter=1;

    trueijkp->i=0;
    trueijkp->p=CENT; // GODMARK -- probably off a bit
    is2=bs;
    ie2=be;

    trueijkp->j=js2=je2=j;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;

  }
  else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){

    trueijkp->iter=2;

    trueijkp->i=is2=ie2=i;

    trueijkp->j=0;
    trueijkp->p=CENT;
    js2=bs;
    je2=be;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){

    trueijkp->iter=3;

    trueijkp->i=is2=ie2=i;

    trueijkp->j=js2=je2=j;

    trueijkp->k=0;
    trueijkp->p=CENT;
    ks2=bs;
    ke2=be;

    di2=1;
    dj2=1;
    dk2=1;
  }


  if(DOENOFLUX != NOENOFLUX){
    dualfprintf(fail_file,"Reinvestigate why ENO needs this\n");
    myexit(9589);
    // JCM: expensive and not sure why doing it
    // reset to 0 so eno schemes don't care about values there (assume weights set to also 0 there)
    for(yiniter=-NBIGBND;yiniter<NBIG+NBIGBND;yiniter++){
      yin[0][yiniter] = yin[1][yiniter] = 0;
    }
  }



  if( (interporflux==ENOINTERPTYPE) || (interporflux==ENOINTERPTYPE4EMF)||(interporflux==ENOFLUXRECONTYPE)||(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE)||(interporflux==ENOQUASIFIELDFLUXRECONTYPE)||(interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE)||(interporflux==ENOAVG2CENTTYPE)||(interporflux==ENOCENT2AVGTYPE) ){// these have only 1 input
    yiniter=bs;
    SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

      yin[0][yiniter]=MACP0A1(p2interpm,i2,j2,k2,pl);
      yiniter++;
    }
  }
  else if(interporflux==ENOFLUXSPLITTYPE){ // this method has 2 inputs
    yiniter=bs;
    SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){

      yin[0][yiniter]=MACP0A1(p2interpm,i2,j2,k2,pl);
      if(p2interpp!=NULL)  yin[1][yiniter]=MACP0A1(p2interpp,i2,j2,k2,pl);
      yiniter++;
    }
  }





}








/// Get 1D line of data so can pass it to ENO scheme
static void get_1d_line_shockarray(int dir, int interporflux, int bs, int ps, int pe, int be,  int i, int j, int k, FTYPE (*shockarray)[NSTORE1][NSTORE2][NSTORE3], FTYPE (*shockindicator)[NBIGM], struct of_trueijkp *trueijkp)
{
  int yiniter;
  int di2,dj2,dk2;
  int i2,j2,k2;
  int is2,ie2,js2,je2,ks2,ke2;
  int dir_exception[NDIM];


  if(
     ((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
     ((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
     ((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
     ){
    dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
    myexit(1);
  }

  dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

  trueijkp->dir=dir;
  trueijkp->interporflux=interporflux;

  // determine range of outer loop and range to feed to eno scheme
  if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){

    trueijkp->iter=1;

    trueijkp->i=0;
    trueijkp->p=CENT; // GODMARK -- probably off a bit
    is2=bs;
    ie2=be;

    trueijkp->j=js2=je2=j;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;

  }
  else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){

    trueijkp->iter=2;

    trueijkp->i=is2=ie2=i;

    trueijkp->j=0;
    trueijkp->p=CENT;
    js2=bs;
    je2=be;

    trueijkp->k=ks2=ke2=k;

    di2=1;
    dj2=1;
    dk2=1;
  }
  else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){

    trueijkp->iter=3;

    trueijkp->i=is2=ie2=i;

    trueijkp->j=js2=je2=j;

    trueijkp->k=0;
    trueijkp->p=CENT;
    ks2=bs;
    ke2=be;

    di2=1;
    dj2=1;
    dk2=1;
  }


  yiniter=bs;
  SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){
    /////////////////////
    // ULTRASUPERGODMARK: Really must average input array to correct location.  Use of shockindicator in para currently assumes indicator at the cell face in dir-direction after using the result of ficalc() from surrounding center cells
    /////////////////////
    shockindicator[EOMSETMHD][yiniter]=MACP1A0(shockarray,SHOCKPLDIR1+dir-1,i2,j2,k2);
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)//KORAL
    shockindicator[EOMSETRAD][yiniter]=MACP1A0(shockarray,SHOCKRADPLDIR1+dir-1,i2,j2,k2);
#endif
    yiniter++;
  }




}





/// Figure out shifting of stencil and order of stencil to more closely match with causality
/// (used for both c2e and other routines)
static void causal_shift_order(int whichprimtype, int interporflux, int dir, int preforder, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, int *shift, int *minorder, int *maxorder)
{
  int i3,j3,k3;
  int is3,ie3,js3,je3,ks3,ke3;
  int di3, dj3, dk3;
  int temporder;
  int superdiv;
  FTYPE wspeed0l,wspeed0r,wspeed1l,wspeed1r;
  FTYPE wspeed0ll,wspeed1ll;
  int yiniter;
  FTYPE localspeed[2];
  int shifttemp;

  int dir_exception[NDIM];

  if(
     ((interporflux==ENOFLUXAVG1TYPE)&&(dir==1))||
     ((interporflux==ENOFLUXAVG2TYPE)&&(dir==2))||
     ((interporflux==ENOFLUXAVG3TYPE)&&(dir==3))
     ){
    dualfprintf(fail_file,"No such method with interporflux=%d and dir=%d\n",interporflux,dir);
    myexit(1);
  }

  dir_exception[1] =  (interporflux==ENOFLUXAVG2TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[2] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG3TYPE);
  dir_exception[3] =  (interporflux==ENOFLUXAVG1TYPE) || (interporflux==ENOFLUXAVG2TYPE);

  if( ( (!dir_exception[1]) && (dir==1) ) || (interporflux==ENOFLUXAVG1TYPE) ){
    is3=ps;
    ie3=pe;

    js3=je3=j;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if( ( (!dir_exception[2]) && (dir==2) ) || (interporflux==ENOFLUXAVG2TYPE) ){
    is3=ie3=i;

    js3=ps;
    je3=pe;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if( ( (!dir_exception[3]) && (dir==3) ) || (interporflux==ENOFLUXAVG3TYPE) ){
    is3=ie3=i;

    js3=je3=j;

    ks3=ps;
    ke3=pe;

    di3=1;
    dj3=1;
    dk3=1;
  }


  //#if( (STOREWAVESPEEDS)&& ((VCHARTYPE==GLOBALVCHAR)||(VCHARTYPE==LOCALVCHAR)) ) // this procedure makes no sense with GLOBALVCHAR
#if( (STOREWAVESPEEDS==1)&& ((VCHARTYPE==LOCALVCHAR)||(VCHARTYPE==VERYLOCALVCHAR) ) )
  // GODMARK: This will use precomputed wave speeds in MACP2A0(wspeed,dir,2,i,j,k)
  // wspeed located at cell interface.  Take this into account when forming shifter.
  // therefore wspeed is wave speeds from interface point of view.  For centered quantities should consider average (or max/min) of wspeed.
  //
  // e.g. for  c2e, average wspeed[i] and wspeed[i+1]
  //      for  a2c (as used for ENOFLUXRECONTYPE, not really center, but edge!) then wspeed[i] is correct one
  //      for  a2em/p average wspeed[i] and wspeed[i+1]
  //
  // shift -order/2-1 (e.g. -2 for WENO5) if flow is superRIGHT (when both left/right chars are +)
  // shift order/2-1 (e.g. 2 for WENO5) if flow is superLEFT (when both left/right chars are -)
  //
  // smoothly vary between and feed integer value (which is interpreted correctly for avg2cent vs. avg2edge vs. cent2edge
  yiniter=ps; // note starts at ps not bs since shift only used on points of interest
  SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

    // first determine if should reduce order because of supersonic divergence
    // first assume maximum preferred order
    temporder=preforder;

    // get standard wavespeed (used by any method)
    wspeed0l=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,0,i3,j3,k3);
    wspeed1l=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,1,i3,j3,k3);   


    if((interporflux==ENOINTERPTYPE)||(interporflux==ENOINTERPTYPE4EMF)||(interporflux==ENOFLUXSPLITTYPE)||(interporflux==ENOFLUXAVG1TYPE)||(interporflux==ENOFLUXAVG2TYPE)||(interporflux==ENOFLUXAVG3TYPE)){ // quantities are at CENT-dir (assumes ENOFLUXAVG?TYPE is orthogonal to dir)


      // get grid-on-right wavespeed
      wspeed0r=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,0,i3+idel,j3+jdel,k3+kdel);
      wspeed1r=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,1,i3+idel,j3+jdel,k3+kdel);


      // check for superfast divergence
      superdiv=0;
      if(SUPERFASTDIVREDUCE&&(wspeed1l<0)&&(wspeed0r>0)){
        superdiv=1;
        temporder=1;
      }

      // get correctly positioned wave speed
      localspeed[0]=min(wspeed0l,wspeed0r);
      localspeed[1]=max(wspeed1l,wspeed1r);

    }
    else if((interporflux==ENOFLUXRECONTYPE)||(interporflux==ENOFLUXRECONTYPEGHOSTACTIVE)||(interporflux==ENOQUASIFIELDFLUXRECONTYPE)||(interporflux==ENOAVG2CENTTYPE)||(interporflux==ENOCENT2AVGTYPE)){ // quantities are at FACE-dir

      // get grid-on-right wavespeed
      wspeed0r=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,0,i3+idel,j3+jdel,k3+kdel);
      wspeed1r=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,1,i3+idel,j3+jdel,k3+kdel);

      // get grid-on-left wavespeed
      wspeed0ll=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,0,i3-idel,j3-jdel,k3-kdel);
      wspeed1ll=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,1,i3-idel,j3-jdel,k3-kdel);


      // check for superfast divergence
      superdiv=0;
      if(SUPERFASTDIVREDUCE&&(wspeed1ll<0)&&(wspeed0r>0)){
        superdiv=1;
        temporder=1;
      }


      // get correctly positioned wave speed as average of surroundings
      // wave speeds at interface
      localspeed[0]=wspeed0l;
      localspeed[1]=wspeed1l;
    }

    // only shift if not diverging superfast
    if(superdiv==0){
      // linear interpolation between left and right super"sonic" directed shifts back downstream
      shifttemp=-((temporder+1)/2)-(temporder+1)*(FTYPE)(localspeed[0]/(localspeed[1]-localspeed[0]+SMALL));
      if(shifttemp>(temporder+1)/2) shifttemp=(temporder+1)/2;
      if(shifttemp<-((temporder+1)/2)) shifttemp=-((temporder+1)/2);
    }
    else shifttemp=0; // no shift if super divergence

    // now assign results to arrays of minorder,maxorder, and shifts.
    minorder[yiniter]=MIN(MINPREFORDER,temporder);
    maxorder[yiniter]=temporder;
    shift[yiniter]=shifttemp;
    yiniter++;
  }
#else
  // NO SHIFT OR CHANGE OF ORDER
  yiniter=ps; // GODMARK: was bs
  SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest  //atch correct -- was di3, dj3, dk3
    minorder[yiniter]=MIN(preforder,MINPREFORDER); // minimum preferred order
    maxorder[yiniter]=preforder;
    shift[yiniter]=0;// then no shift since don't have wave speeds (not stored)
    yiniter++;
  }
#endif

}




/// determine shock indicator (used for both c2e and other routines, only operates on quantities with NPR elements)
static int get_V_and_P(int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yrealin)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp)
{
  int num,locali,localj,localk,localloc;
  int di,dj,dk;
  int i3,j3,k3;
  int is3,ie3,js3,je3,ks3,ke3;
  int di3, dj3, dk3;
  FTYPE interplistpl[MAXNPR][MAXSPACEORDER];
  FTYPE *yrealpl[MAXNPR]; // number of pointers
  int plpl;
  int startorderi,endorderi;
  int yiniter;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int iii;
  FTYPE myprim[MAXNPR],bsq;
  struct of_state qdontuse;
  struct of_state *qptr=&qdontuse;
  int pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;



  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);








  if(((interporflux!=ENOFLUXAVG1TYPE)&&(dir==1))||(interporflux==ENOFLUXAVG1TYPE)){
    is3=bs;
    ie3=be;

    js3=je3=j;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG2TYPE)&&(dir==2))||(interporflux==ENOFLUXAVG2TYPE)){
    is3=ie3=i;

    js3=bs;
    je3=be;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG3TYPE)&&(dir==3))||(interporflux==ENOFLUXAVG3TYPE)){
    is3=ie3=i;

    js3=je3=j;

    ks3=bs;
    ke3=be;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else{
    dualfprintf(fail_file,"No such dir=%d in get_V_and_P()\n",dir);
    myexit(246344572);
  }

  // trueijkp->iter should really be passed
  di=(trueijkp->iter==1);
  dj=(trueijkp->iter==2);
  dk=(trueijkp->iter==3);


  yiniter=bs; // per-point, get maximum can get
  SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

    num=iii=i3*idel+j3*jdel+k3*kdel;

    // trueijkp->i,trueijkp->j,trueijkp->k,trueijkp->p should really be passed
    // needed for EOS in general
    locali=trueijkp->i+di*num;
    localj=trueijkp->j+dj*num;
    localk=trueijkp->k+dk*num;
    localloc=trueijkp->p;

#if(PLINEWITHFIELD || VLINEWITHGDETRHO)
    PALLREALLOOP(pl) myprim[pl] = yrealin[pl][0][num];
    get_geometry(locali,localj,localk,localloc,ptrgeom);
    // need u^\mu and Bsq (currently from b^\mu b_\mu)
    // OPTMARK: get_state here is very expensive since only required to have per-centered-point and not per direction, so 3X over computed unless refer to stored data.
    // GODMARK: Using stored data presumes yrealin is at CENT, which is true for any use of get_V_and_P() right now
    // Actually do have data at left-right for each dir in fluxstate[dir].  So can check where yrealin is located using localglobal.
    get_stateforinterpline(myprim,ptrgeom,&qptr); // OPTMARK :Refer to centered state if stored
#endif



#if(VLINEWITHGDETRHO==0)
    Vline[EOMSETMHD][yiniter]=yrealin[UU+dir][0][iii];
#else
    // \detg rho u^dir (approximately)
    //    get_geometry_gdetonly(locali,localj,localk,localloc,ptrgeom);
    //    Vline[yiniter]=(ptrgeom->g)*yrealin[RHO][0][iii]*yrealin[UU+dir][0][iii];
    // \detg rho u^dir
    Vline[EOMSETMHD][yiniter]=(ptrgeom->gdet)*yrealin[RHO][0][num]*(qptr->ucon[dir]);
#endif

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
#if(VLINEWITHGDETRHO==0)
    Vline[EOMSETRAD][yiniter]=yrealin[PRAD0+dir][0][iii];
#else
    Vline[EOMSETRAD][yiniter]=(ptrgeom->gdet)*yrealin[PRAD0][0][num]*(qptr->uradcon[dir]); // approximate KORALTODO
#endif
#endif



    // OPTMARK: Since assuming at CENT with normal primitive, then already have correct pressure
    //    P[yiniter]=pressure_rho0_u_simple(locali,localj,localk,localloc,yrealin[RHO][0][iii],yrealin[UU][0][iii]);
    Pline[EOMSETMHD][yiniter]=qptr->pressure;
#if(PLINEWITHFIELD==0)
    // done
#else
    // need total pressure in general
    // add magnetic pressure
    // OPTMARK: Since assuming at CENT with normal primitive, then already have correct b^2
    //    bsq = dot(q.bcon, q.bcov); // now store bsq
    Pline[EOMSETMHD][yiniter] += 0.5*(qptr->bsq);
#endif

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
    // add radiation pressure to total pressure if optically thick
    FTYPE tautot[NDIM],tautotmax;
    //calcfull_tautot(myprim, ptrgeom, tautot, &tautotmax);
    calc_tautot(myprim, ptrgeom, NULL, tautot, &tautotmax); // very accurate tautot not necessary, so use Tgas=Trad assumption in opacity
    Pline[EOMSETRAD][yiniter]=(4.0/3.0-1.0)*yrealin[PRAD0][0][num]; // radiative pressure
    // add radiation pressure to total pressure if optically thick
    // KORALNOTE: recall pressure just along diagonal and no velocity in R^\mu_\nu
    Pline[EOMSETMHD][yiniter]+=MIN(tautotmax,1.0)*Pline[EOMSETRAD][yiniter];
#endif


    yiniter++;
  }

  return( 0 );

}




/// determine shock indicator (used for both c2e and other routines, only operates on quantities with NPR elements)
static int get_shock_indicator(int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yprim)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], FTYPE (*shockindicator)[NBIGM], struct of_trueijkp *trueijkp)
{
  int i3,j3,k3;
  int is3,ie3,js3,je3,ks3,ke3;
  int di3, dj3, dk3;
#if(0)
  FTYPE yinterplistpl[MAXNPR][MAXSPACEORDER];
  FTYPE *ypl[MAXNPR]; // number of pointers
  int plpl;
  int startorderi,endorderi;
  int l;
#endif
  int yiniter;
  extern FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int num,locali,localj,localk,localloc;
  int di,dj,dk;




  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  // NO!  While Ftilde() doesn't yet use ypl, feed get_shock_indicator() yprim, so need full prim loop
  //  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  // bs+2 .. be-2 since Ficalc needs +-2 points



  if(((interporflux!=ENOFLUXAVG1TYPE)&&(dir==1))||(interporflux==ENOFLUXAVG1TYPE)){
    is3=bs+2;
    ie3=be-2;

    js3=je3=j;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG2TYPE)&&(dir==2))||(interporflux==ENOFLUXAVG2TYPE)){
    is3=ie3=i;

    js3=bs+2;
    je3=be-2;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG3TYPE)&&(dir==3))||(interporflux==ENOFLUXAVG3TYPE)){
    is3=ie3=i;

    js3=je3=j;

    ks3=bs+2;
    ke3=be-2;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else{
    dualfprintf(fail_file,"No such interporflux=%d in get_shock_indicator()\n",interporflux);
    myexit(9894386);
  }



  yiniter=bs+2; // starts at bs+2

  if(MAXBND<4){
    dualfprintf(fail_file,"MAXBND should be 4 for shockindicator???\n");
    myexit(8465684);
  }


#if(0)
  startorderi = - (7)/2; // order=7 fixed for shock detector
  endorderi   = - startorderi;
  // shift pointer
  PALLREALLOOP(plpl){
    ypl[plpl] = yinterplistpl[plpl] - startorderi;
  }
#endif

  // trueijkp->iter should really be passed
  di=(trueijkp->iter==1);
  dj=(trueijkp->iter==2);
  dk=(trueijkp->iter==3);


  SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

#if(0)
    // Ficalc doesn't really use ypl right now
    PALLREALLOOP(plpl){
      // get interpolation points, where y[0] is point of interest for which interpolation is found.
      for(l=startorderi;l<=endorderi;l++){
        ypl[plpl][l]=yprim[plpl][0][i3*idel+j3*jdel+k3*kdel + l];
      }
    }
#endif
    //    dualfprintf(fail_file,"idel=%d jdel=%d kdel=%d :: i3=%d j3=%d k3=%d\n",idel,jdel,kdel,i3,j3,k3);

    //shockindicator[EOMSETMHD][yiniter]=Ficalc(dir,&Vline[EOMSETMHD][yiniter],&Pline[EOMSETMHD][yiniter],ypl);
    shockindicator[EOMSETMHD][yiniter]=Ficalc(dir,&Vline[EOMSETMHD][yiniter],&Pline[EOMSETMHD][yiniter]);

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
    shockindicator[EOMSETRAD][yiniter]=Ficalc(dir,&Vline[EOMSETRAD][yiniter],&Pline[EOMSETRAD][yiniter]); // KORAL
#endif


#if(0 && FLUXDUMP) // DEBUG: (turn off FLUXDUMP in flux.c)
    // trueijkp->i,trueijkp->j,trueijkp->k,trueijkp->p should really be passed
    // ijkcurr needed for EOS in general
    num=i3*idel+j3*jdel+k3*kdel;
    locali=trueijkp->i+di*num;
    localj=trueijkp->j+dj*num;
    localk=trueijkp->k+dk*num;
    localloc=trueijkp->p;

    GLOBALMACP0A1(fluxdump,locali,localj,localk,4*NPR + (dir-1)*NPR*5 + NPR*0 + RHO)=shockindicator[EOMSETMHD][yiniter]; // don't look at EOMSETRAD

#endif


    yiniter++;
  }

  return( 0 );

}



/// currently only sets density etai, rest are 0 GODMARK
static int get_contact_indicator(int realisinterp, int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM])
{
  int i3,j3,k3;
  int is3,ie3,js3,je3,ks3,ke3;
  int di3, dj3, dk3;
  FTYPE yinterplistpl[MAXNPR][MAXSPACEORDER];
  FTYPE *ypl[MAXNPR]; // number of pointers
  int plpl, pliter;
  int startorderi,endorderi;
  int yiniter;
  int l;
  extern FTYPE etaicalc(int pl, FTYPE *y, FTYPE *V, FTYPE *P);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;




  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  // input is expected to be line of data, so correct to use NUMPRIMLOOP, but in the end only apply if realisinterp==1 and doing pl==RHO
  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  // bs+2 .. be-2 since etaicalc needs +-2 points




  if(((interporflux!=ENOFLUXAVG1TYPE)&&(dir==1))||(interporflux==ENOFLUXAVG1TYPE)){
    is3=bs+2;
    ie3=be-2;

    js3=je3=j;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG2TYPE)&&(dir==2))||(interporflux==ENOFLUXAVG2TYPE)){
    is3=ie3=i;

    js3=bs+2;
    je3=be-2;

    ks3=ke3=k;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else if(((interporflux!=ENOFLUXAVG3TYPE)&&(dir==3))||(interporflux==ENOFLUXAVG3TYPE)){
    is3=ie3=i;

    js3=je3=j;

    ks3=bs+2;
    ke3=be-2;

    di3=1;
    dj3=1;
    dk3=1;
  }
  else{
    dualfprintf(fail_file,"No such interporflux=%d in get_contact_indicator()\n",interporflux);
    myexit(9894386);
  }



  yiniter=bs+2; // bs+2

  if(MAXBND<4){
    dualfprintf(fail_file,"MAXBND should be 4 for contactindicator???\n");
    myexit(8465684);
  }

  startorderi = - (7)/2; // order=7 fixed for shock detector (must make sure MAXBND==5)
  endorderi   = - startorderi;

  // shift pointer
  NUMPRIMLOOP(pliter,plpl){
    ypl[plpl] = yinterplistpl[plpl] - startorderi;
  }


  SUPERGENLOOP(i3,j3,k3,is3,ie3,js3,je3,ks3,ke3,di3,dj3,dk3){ // only over points of interest

    NUMPRIMLOOP(pliter,plpl){


      if(plpl==RHO && realisinterp){
        // get interpolation points, where y[0] is point of interest for which interpolation is found.
        for(l=startorderi;l<=endorderi;l++){
          ypl[plpl][l]=yin[plpl][0][i3*idel+j3*jdel+k3*kdel + l];
        }

        etai[plpl][EOMSETMHD][yiniter]=etaicalc(plpl, ypl[plpl], &Vline[EOMSETMHD][yiniter], &Pline[EOMSETMHD][yiniter]);
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
        etai[plpl][EOMSETRAD][yiniter]=etaicalc(plpl, ypl[plpl], &Vline[EOMSETRAD][yiniter], &Pline[EOMSETRAD][yiniter]);
#endif
      }
      else{
        etai[plpl][EOMSETMHD][yiniter]=0.0;
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
        etai[plpl][EOMSETRAD][yiniter]=0.0;
#endif
      }

    }

    yiniter++;
  }

  return( 0 );



}




/// real calculation.  Gets derivatives of input function and passes this to various other function
static int compute_df_line(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM])
{
  int i;
  extern void get_limit_slopes_paraline(FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dq);



  ///////////////////////////
  //
  // check that NUMDFS is sufficient
  //
  ///////////////////////////

  // if doingweno==0 then assume things needed are like what's needed for paraline
  if(NUMDFS<4){
    dualfprintf(fail_file,"paraline requires NUMDFS=4 while set to %d\n",NUMDFS);
    myexit(8923658);
  }
  if(DOMONOINTERP!= NOMONOINTERP && NUMDFS<5){
    dualfprintf(fail_file,"MONO requires NUMDFS=5 while set to %d\n",NUMDFS);
    myexit(8923659);
  }


  // for SWENO: df[0] was one-sided :: df[1] was d^2f :: df[3] was centered diff (but no factor of 0.5!!! GODMARK) df[4] was df centered 2 apart

#if(NUMDFS>DFONESIDED)
  // df
  // df (one-sided difference : used to get Dqp,Dqm)
  // truly centered on face if y on CENT, but effectively used for CENT
  for(i=bs+1;i<=be;i++) df[DFONESIDED][i] = yin[i]-yin[i-1];
#endif

  // df (centered difference : used to get Dqc)
  // On CENT if y CENT
  // df centered, for Sasha's smoothness indicators
#if(NUMDFS>DFCENT)
  for(i=bs+1;i<=be-1;i++) df[DFCENT][i] = 0.5*(yin[i+1]-yin[i-1]);
#endif


#if(NUMDFS>DFCENT2APART)
  // df centered 2 apart, for Sasha's smoothness indicators
  for(i=bs+2;i<=be-2;i++) df[DFCENT2APART][i] = 0.25*(yin[i+2]-yin[i-2]);
#endif

#if(NUMDFS>DF2OFONESIDED)
  // d^2f
  for(i=bs+1;i<=be-1;i++) df[DF2OFONESIDED][i] = df[DFONESIDED][i+1]-df[DFONESIDED][i];
#endif


#if(NUMDFS>DFMONO)
  // PARALIM 3-point MONOTONIZED slope
  // on face if y CENT and so df[DFCENT] CENT
  for(i=bs+1;i<=be-1;i++){
    get_limit_slopes_paraline(&df[DFONESIDED][i], &df[DFONESIDED][i+1], &df[DFCENT][i], &df[DFMONO][i]);
  }
#endif

#if(NUMDFS>DF2OFMONO)
  // d^2f (compute second derivative on *monotonized* slopes used for para)
  // at CENT if y at CENT
  for(i=bs+2;i<=be-1;i++) df[DF2OFMONO][i] = df[DFMONO][i]-df[DFMONO][i-1];
#endif

  // assume all other things computed within PARA-based routines



  return(0);
}


/// real calculation.  Gets derivatives of input function and passes this to various other function
/// Only what paraline needs [removed DFCENT2APART and DF2OFONESIDED calculations as compared to compute_df_line()]
static int compute_df_line_paraline(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM])
{
  int i;
  extern void get_limit_slopes_paraline(FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dq);



  ///////////////////////////
  //
  // check that NUMDFS is sufficient
  //
  ///////////////////////////

  // if doingweno==0 then assume things needed are like what's needed for paraline
#if(NUMDFS<4)
  dualfprintf(fail_file,"paraline requires NUMDFS=4 while set to %d\n",NUMDFS);
  myexit(8923658);
#endif


#if(NUMDFS>DFONESIDED)
  // df
  // df (one-sided difference : used to get Dqp,Dqm)
  // truly centered on face if y on CENT, but effectively used for CENT
  for(i=bs+1;i<=be;i++) df[DFONESIDED][i] = yin[i]-yin[i-1];
#endif

  // df (centered difference : used to get Dqc)
  // On CENT if y CENT
  // df centered, for Sasha's smoothness indicators
#if(NUMDFS>DFCENT)
  for(i=bs+1;i<=be-1;i++) df[DFCENT][i] = 0.5*(yin[i+1]-yin[i-1]);
#endif



#if(NUMDFS>DFMONO)
  // PARALIM 3-point MONOTONIZED slope
  // on face if y CENT and so df[DFCENT] CENT
  for(i=bs+1;i<=be-1;i++){
    get_limit_slopes_paraline(&df[DFONESIDED][i], &df[DFONESIDED][i+1], &df[DFCENT][i], &df[DFMONO][i]);
  }
#endif

#if(NUMDFS>DF2OFMONO)
  // d^2f (compute second derivative on *monotonized* slopes used for para)
  // at CENT if y at CENT
  for(i=bs+2;i<=be-1;i++) df[DF2OFMONO][i] = df[DFMONO][i]-df[DFMONO][i-1];
#endif

  // assume all other things computed within PARA-based routines



  return(0);
}


/// Gets df and ddf for MONO if mono will use different yin than computed in above compute_df_line() function
static int compute_df_line_formono(int doingweno,int interporflux, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE *yin, FTYPE (*df)[NBIGM])
{
  int i;



#if(NUMDFS>DFONESIDED)
  // df
  // df (one-sided difference : used to get Dqp,Dqm)
  // truly centered on face if y on CENT, but effectively used for CENT
  for(i=bs+1;i<=be;i++) df[DFONESIDED][i] = yin[i]-yin[i-1];
#else
  if(DOMONOINTERP==SMONOINTERP){
    dualfprintf(fail_file,"Need NUMDFS=%d bigger for SMONO to be used1\n",NUMDFS);
  }
#endif

#if(NUMDFS>DF2OFONESIDED)
  // d^2f
  for(i=bs+1;i<=be-1;i++) df[DF2OFONESIDED][i] = df[DFONESIDED][i+1]-df[DFONESIDED][i];
#else
  if(DOMONOINTERP==SMONOINTERP){
    dualfprintf(fail_file,"Need NUMDFS=%d bigger for SMONO to be used2\n",NUMDFS);
  }
#endif





  return(0);
}


/// returns the improved pressure jump indicators
int compute_df_line_new(int doingweno, int whichprimtype, int interporflux, int recontype, int dir, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*yprim)[2][NBIGM], FTYPE *df, FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], struct of_trueijkp *trueijkp)
{
  int i;
  int num,locali,localj,localk,localloc;
  int pl,pliter;
  int di,dj,dk;
  FTYPE myprim[MAXNPR];
  FTYPE mypriml[MAXNPR];
  FTYPE myprimr[MAXNPR];
  FTYPE myPgas,myPtot;
  //  FTYPE U[MAXNPR];
  struct of_state qdontuse;
  struct of_state *qptr=&qdontuse;
  FTYPE bsq;
  FTYPE btt;
  FTYPE par,parl,parr;
  FTYPE uparl,uparr,upar;
  FTYPE veffl,veffr,veff;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  FTYPE stifffactor,stifffactor1,stifffactor2;
  FTYPE a_Pgasline[NBIGM];
  FTYPE (*Pgasline);
  FTYPE a_bsqline[NBIGM];
  FTYPE (*bsqline);
  FTYPE a_parline[NBIGM];
  FTYPE (*parline);
  FTYPE a_uparline[NBIGM];
  FTYPE (*uparline);
  FTYPE a_veffline[NBIGM];
  FTYPE (*veffline);
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;





  if(doingweno==0){
    dualfprintf(fail_file,"Shouldn't be in compute_df_line_new if not doing WENO\n");
    myexit(7671515);
  }

  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  // we don't deal with interpolated quantities here, only real primitives
  //  setup_nprlocalist(whichprimtype,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  // shift 1-D arrays
  Pgasline=(FTYPE (*)) (&(a_Pgasline[NBIGBND]));
  bsqline=(FTYPE (*)) (&(a_bsqline[NBIGBND]));
  parline=(FTYPE (*)) (&(a_parline[NBIGBND]));
  uparline=(FTYPE (*)) (&(a_uparline[NBIGBND]));
  veffline=(FTYPE (*)) (&(a_veffline[NBIGBND]));


  // trueijkp->iter should really be passed
  di=(trueijkp->iter==1);
  dj=(trueijkp->iter==2);
  dk=(trueijkp->iter==3);


  for(num=bs;num<=be;num++){

    // trueijkp->i,trueijkp->j,trueijkp->k,trueijkp->p should really be passed
    // needed for EOS in general
    locali=trueijkp->i+di*num;
    localj=trueijkp->j+dj*num;
    localk=trueijkp->k+dk*num;
    localloc=trueijkp->p;


    PALLREALLOOP(pl) myprim[pl] = yprim[pl][0][num];
    get_geometry(locali,localj,localk,localloc,ptrgeom);

    // SUPERGODMARK: OPTMARK: Assume prim at CENT for now, which is consistent with other places but not really true in general when averaging primitive locations for old a2c method
    //    if(interporflux==ENOINTERPTYPE){
    get_stateforinterpline(myprim,ptrgeom,&qptr);
    //    }
    //    else{
    //      // then can't be sure quantity is at CENT
    //      get_state(myprim,ptrgeom,qptr);
    //    }


    Pgasline[num]=pressure_rho0_u_simple(locali,localj,localk,localloc,myprim[RHO],myprim[UU]); // just gas contribution, add magnetic below

#if(VLINEWITHGDETRHO==0)
    //    Vline[EOMSETMHD][num]=myprim[UU+dir];
    //    Vline[EOMSETMHD][num]=(qptr->ucon[dir])/(qptr->ucon[TT]); // more correct generally since wanted w.r.t. grid
    Vline[EOMSETMHD][num]=qptr->ucon[dir]; // more correct generally since wanted w.r.t. grid (need 4-velocity since consistent with equations of motion)
#else
    Vline[EOMSETMHD][num]=(ptrgeom->gdet)*yprim[RHO][0][num]*(qptr->ucon[dir]);
#endif

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
#if(VLINEWITHGDETRHO==0)
    Vline[EOMSETRAD][num]=qptr->uradcon[dir];
#else
    Vline[EOMSETRAD][num]=(ptrgeom->gdet)*yprim[PRAD0][0][num]*(qptr->uradcon[dir]); // approximate KORALTODO
#endif
#endif


    //    bsq = dot(q.bcon, q.bcov);
    bsq = qptr->bsq;
    btt = -qptr->bcon[dir] * qptr->bcov[TT];
    myPgas = Pgasline[num];
    myPtot = myPgas + 0.5*bsq;

    // get total pressure
    Pline[EOMSETMHD][num]=myPtot;
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE) // KORAL
    // add radiation pressure to total pressure if optically thick
    FTYPE tautot[NDIM],tautotmax;
    //    calcfull_tautot(myprim, ptrgeom, tautot, &tautotmax);
    calc_tautot(myprim, ptrgeom, NULL, tautot, &tautotmax); // very accurate tautot not necessary, so use Tgas=Trad assumption in opacity
    Pline[EOMSETRAD][num]=(4.0/3.0-1.0)*yprim[PRAD0][0][num]; // radiative pressure
    // add radiation pressure to total pressure if optically thick
    // KORALNOTE: recall pressure just along diagonal and no velocity in R^\mu_\nu
    Pline[EOMSETMHD][num]+=MIN(tautotmax,1.0)*Pline[EOMSETRAD][num];
#endif

    bsqline[num]=fabs(bsq);

    parline[num]=fabs(qptr->ifremoverestplus1ud0elseud0*(qptr->ucon[dir]));

    // note sign comes from u^{dir}, but otherwise correct velocity to go into -T^t_t
    veffline[num] = sign(Vline[EOMSETMHD][num])*parline[num]/(sqrt(parline[num])+SMALL); // GODMARK: sqrt expensive, not sure how to avoid

    // upar is really just the non-kinetic part of -T^t_t
    // compute it like below so don't repeat expensive plus1ud0 calculation
    uparline[num] = fabs( (myprim[UU] + myPgas  + bsq )*( fabs(qptr->ucon[dir])  *qptr->ucov[TT])  + myPtot  + btt );
  }




  /////////////////////
  //
  // can go all the way to edge of grid for C2E for stiffness indicator
  //
  /////////////////////
  if(recontype==CVT_C2E){
    for(num=bs;num<=be;num++){

      // assignments
      PALLREALLOOP(pl) myprim[pl] = yprim[pl][0][num];
      bsq=bsqline[num];
      par=parline[num];
      upar=uparline[num];
      veff=veffline[num];
      
      stifffactor = max(max(fabs(bsq/myprim[RHO]),fabs(bsq/myprim[UU])),fabs(par));
      // now turn into something that goes from 0 to 1
      // assume moderatively relativistic is ok (stifffactor~3, but higher is sufficiently stiff, so truncate at 3)
      stiffindicator[num] = max(min(stifffactor-2.0,1.0),0.0);
    }
  }



  /////////////////////
  //
  // Get df type shock indicator (and also stiffindicator if not C2E)
  //
  /////////////////////
  for(num=bs+1;num<=be-1;num++){

    // assignments
    PALLREALLOOP(pl) myprim[pl] = yprim[pl][0][num];
    bsq=bsqline[num];
    par=parline[num];
    parl=parline[num-1];
    parr=parline[num+1];
    upar=uparline[num];
    uparl=uparline[num-1];
    uparr=uparline[num+1];
    veff=veffline[num];
    veffl=veffline[num-1];
    veffr=veffline[num+1];


    //The expression is:
    // -(rho v^2/2 + u) = rho * ucon * (1+ucov)  +  gam * u * (ucon*ucov) + (gam-1)*u

    // Original dP/P prescription
    //df[num] = ( fabs(uparr-uparl) ) / ( fabs(upar) + DBL_MIN );

    //normalized difference using the central value -- think may lead to problems
    //df[num] = (fabs(myprim[RHO]*0.5*(parr-parl)) + fabs(0.5*(uparr-uparl)))/(myprim[RHO]*par + upar);

    //normalized difference without using the central value; normalize the difference by the values themselves
    //df[num] = ( fabs(myprim[RHO]*(parr-parl)) + fabs(uparr-uparl) ) / ( fabs(myprim[RHO])*(fabs(parr)+fabs(parl)) + fabs(uparr)+fabs(uparl) + DBL_MIN );

    //divide by one point value, not by the sum of the values themselves
    // below is like: rho delta(v^2)/u
    //    df[num] = ( fabs(myprim[RHO]*(parr-parl)) + fabs(uparr-uparl) ) / ( fabs(myprim[RHO]*par) + fabs(upar) + DBL_MIN );

    // below is now like: rho (delta v)^2/u -- consistent with (\delta Mach)^2 that goes into internal energy error estimate
    df[num] = ( fabs(myprim[RHO]*(veffr-veffl)*(veffr-veffl)) + fabs(uparr-uparl) ) / ( fabs(myprim[RHO]*par) + fabs(upar) + SMALL );



    // equations become stiff when b^2/rho_0, b^2/u, or (\gamma-1) become sufficiently relativistic
    // now also include above-computed df that is for supersonic motion that is stiff
    //    stifffactor = max(max(max(fabs(bsq/myprim[RHO]),fabs(bsq/myprim[UU])),fabs(par)),fabs(df[num]));

    // JCM: for now don't include df[num] since often true and shouldn't really affect C2E that only uses this right now
    // JCM: should directly use df[num] in A2C/C2A to determine if to do those operations or not

    if(recontype==CVT_C2E){
      // then already done above
    }
    else if(recontype==CVT_C2A || recontype==CVT_A2C){
      // CHANGINGMARK: can use above for c2e but need df[num] version for a2c/c2a/etc. to avoid ....?
      //      stifffactor = max(max(fabs(bsq/myprim[RHO]),fabs(bsq/myprim[UU])),fabs(par));
      stifffactor = fabs(par);
      // now turn into something that goes from 0 to 1
      stifffactor1 = max(min(stifffactor-2.0,1.0),0.0);
      //      stifffactor = max(df[num],stifffactor);
      // now turn into something that goes from 0 to 1
      stifffactor2 = max(min(fabs(df[num])*30.0,1.0),0.0);

      // now turn into something that goes from 0 to 1
      stiffindicator[num]=max(min(max(stifffactor1,stifffactor2),1.0),0.0);

      //      stiffindicator[num]=1.0;
    }




#if(0 && FLUXDUMP) // DEBUG: (turn off FLUXDUMP in flux.c)
    // trueijkp->i,trueijkp->j,trueijkp->k,trueijkp->p should really be passed
    // ijkcurr needed for EOS in general
    locali=trueijkp->i+di*num;
    localj=trueijkp->j+dj*num;
    localk=trueijkp->k+dk*num;
    localloc=trueijkp->p;

    GLOBALMACP0A1(fluxdump,locali,localj,localk,4*NPR + (dir-1)*NPR*5 + NPR*0 + RHO)=stifffactor;
    GLOBALMACP0A1(fluxdump,locali,localj,localk,4*NPR + (dir-1)*NPR*5 + NPR*0 + UU)=stiffindicator[num];
    GLOBALMACP0A1(fluxdump,locali,localj,localk,4*NPR + (dir-1)*NPR*5 + NPR*0 + U1)=veff;

    //    dualfprintf(fail_file,"%ld %d :: dir=%d :: ijkcurr=%d %d %d :: %21.15g %21.15g %21.15g\n",nstep,steppart,dir,locali,localj,localk,stifffactor,stiffindicator[num],veff);

#endif





  }

  return(0);
}





/// Assign result of ENO operation to final array
/// for c2e
static void assign_eno_result_c2e_multipl(int whichquantity, int recontype, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[2][NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*result1)[NSTORE2][NSTORE3][NPR2INTERP])
{
  // for NUMPRIMLOOP:
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  // others:
  int l;
  int ii,jj,kk;


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  // where to place in real final array
  // originally had idel*l, jdel*l, kdel*l inside loop.  Was more expensive and was one of most expensive lines in code.

  // note that NUMPRIMLOOP(pliter,pl) is inside for(l) loop since (*result)[][][] has larger memory jumps for l than (*yout)[][]
  // Small difference, since apparently sometimes run and having pl on outside is faster overall, but % usage is better with pl deep inside.
  if(idel==1){
    jj=j;
    kk=k;
    for(l=ps;l<=pe;l++){ // inclusive loop
      ii=i+l;
      NUMPRIMLOOP(pliter,pl) MACP0A1(result0,ii,jj,kk,pl) = yout[pl][0][l];
      NUMPRIMLOOP(pliter,pl) MACP0A1(result1,ii,jj,kk,pl) = yout[pl][1][l];
    }
  }
  else if(jdel==1){
    ii=i;
    kk=k;
    for(l=ps;l<=pe;l++){ // inclusive loop
      jj=j+l;
      NUMPRIMLOOP(pliter,pl) MACP0A1(result0,ii,jj,kk,pl) = yout[pl][0][l];
      NUMPRIMLOOP(pliter,pl) MACP0A1(result1,ii,jj,kk,pl) = yout[pl][1][l];
    }
  }
  else if(kdel==1){
    ii=i;
    jj=j;
    for(l=ps;l<=pe;l++){ // inclusive loop
      kk=k+l;
      NUMPRIMLOOP(pliter,pl) MACP0A1(result0,ii,jj,kk,pl) = yout[pl][0][l];
      NUMPRIMLOOP(pliter,pl) MACP0A1(result1,ii,jj,kk,pl) = yout[pl][1][l];
    }
  }
  else{
    dualfprintf(fail_file,"No such idel=%d jdel=%d kdel=%d in assign_eno_result_c2e()\n",idel,jdel,kdel);
    myexit(3469836);
  }


  if(recontype!=CVT_C2E){
    dualfprintf(fail_file,"assign_eno_result_c2e only handles recontype==CVT_C2E\n");
    myexit(26);
  }

}



/// Assign result of ENO operation to final array
/// for c2e
static void assign_eno_result_c2e(int recontype, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*result1)[NSTORE2][NSTORE3][NPR2INTERP])
{
  int l;
  int ii,jj,kk;

  // where to place in real final array
  // originally had idel*l, jdel*l, kdel*l inside loop.  Was more expensive and was one of most expensive lines in code.

  if(idel==1){
    jj=j;
    kk=k;
    for(l=ps;l<=pe;l++){ // inclusive loop
      ii=i+l;
      MACP0A1(result0,ii,jj,kk,pl) = yout[0][l];
      MACP0A1(result1,ii,jj,kk,pl) = yout[1][l];
    }
  }
  else if(jdel==1){
    ii=i;
    kk=k;
    for(l=ps;l<=pe;l++){ // inclusive loop
      jj=j+l;
      MACP0A1(result0,ii,jj,kk,pl) = yout[0][l];
      MACP0A1(result1,ii,jj,kk,pl) = yout[1][l];
    }
  }
  else if(kdel==1){
    ii=i;
    jj=j;
    for(l=ps;l<=pe;l++){ // inclusive loop
      kk=k+l;
      MACP0A1(result0,ii,jj,kk,pl) = yout[0][l];
      MACP0A1(result1,ii,jj,kk,pl) = yout[1][l];
    }
  }
  else{
    dualfprintf(fail_file,"No such idel=%d jdel=%d kdel=%d in assign_eno_result_c2e()\n",idel,jdel,kdel);
    myexit(3469836);
  }


  if(recontype!=CVT_C2E){
    dualfprintf(fail_file,"assign_eno_result_c2e only handles recontype==CVT_C2E\n");
    myexit(26);
  }

}


/// Assign result of ENO operation to final array
/// only care about how many outputs
static void assign_eno_result(int recontype, int pl, int bs, int ps, int pe, int be, int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yout)[NBIGM], FTYPE (*result0)[NSTORE2][NSTORE3][NPR], FTYPE (*result1)[NSTORE2][NSTORE3][NPR])
{
  int l;
  int ti,tj,tk;


  if(recontype==CVT_C2E){// these methods generated 2 output values
    for(l=ps;l<=pe;l++){ // inclusive loop
      MACP0A1(result0,i+l*idel,j+l*jdel,k+l*kdel,pl) = yout[0][l];
      MACP0A1(result1,i+l*idel,j+l*jdel,k+l*kdel,pl) = yout[1][l];
    }
  }
  else if( (recontype==CVT_C2A)||(recontype==CVT_A2C) ){// these methods generated 1 output value
    for(l=ps;l<=pe;l++){ // inclusive loop
      // result0 can be equal to input p2interpm array since now done with that line and dimensionally split, so never operate in another direction on this quantity

      ti=i+l*idel;
      tj=j+l*jdel;
      tk=k+l*kdel;

      if(crapdebug==0) MACP0A1(result0,i+l*idel,j+l*jdel,k+l*kdel,pl) = yout[0][l];
      else{
        if(pl!=3) MACP0A1(result0,i+l*idel,j+l*jdel,k+l*kdel,pl) = yout[0][l];
        // MACP0A1(result0,i+l*idel,j+l*jdel,k+l*kdel,pl) = yout[0][l];


        if(ti<-N1BND || ti>N1+N1BND-1 || tj<-N2BND || tj>N2+N2BND-1 || tk<-N3BND || tk>N3+N3BND-1 || pl<0 || pl>U3){
          dualfprintf(fail_file,"OUT OF BOUNDS\n");
        }
      }
      
#if(0)
      if(crapdebug && pl==3){
        dualfprintf(fail_file,"ti=%d tj=%d tk=%d :: i=%d j=%d k=%d l=%d idel=%d jdel=%d kdel=%d pl=%d\n",i+l*idel,j+l*jdel,k+l*kdel,i,j,k,l,idel,jdel,kdel,pl);
      }
#endif
    }
  }

}








