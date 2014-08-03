
/*! \file metric.tools.funcdeclare.h
     \brief Function declarations for global use of things in metric.tools.c

*/

extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);

// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;


extern FTYPE rhor_calc(int which);
extern FTYPE rmso_calc(int which) ;
extern FTYPE uphi_isco_calc(int which,FTYPE r);
// metric_tools.c:

extern int gdet_func_metric(int whichcoord, FTYPE *V,FTYPE *gcov, FTYPE *gdet);
extern int gdet_func(int whichcoord, FTYPE *gcov, FTYPE *gdet);
extern int gdet_func_singcheck(int whichcoord, FTYPE *V,FTYPE (*generalmatrixlower)[NDIM], FTYPE *gdet);
extern void metric_sing_check(int whichcoord, FTYPE (*genmatrixlower)[NDIM], int *anglesing, int*centersing, int *truedim);
extern void gdetvol_func(struct of_geom *ptrgeom, FTYPE *gdet, FTYPE *eomfunc, FTYPE *gdetvol);
extern void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *EOMFUNCNAME);

extern void matrix_inverse_metric(int whichcoord, FTYPE *gcov, FTYPE *gcon);
extern void matrix_inverse(int whichcoord, FTYPE (*genmatrixlower)[NDIM], FTYPE (*genmatrixupper)[NDIM]);
extern void matrix_inverse_4d(FTYPE (*genmatrixlower)[NDIM], FTYPE (*genmatrixupper)[NDIM]);
extern void alphalapse_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE *alphalapse);
extern void betasqoalphasq_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE *betasqoalphasq);
extern void beta_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE alphalapse, FTYPE *beta);


extern void gset_genloc(int getprim, int whichcoord, int i, int j, int k, int loc, struct of_geom *geom);
extern void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *geom);
extern void gset_X(int getprim, int whichcoord, int i, int j, int k, int loc, FTYPE *X, struct of_geom *ptrgeom);



extern void get_and_copy_geometry(int ii, int jj, int kk, int pp, struct of_geom *ptrgeom);


extern void set_igdet_old(struct of_geom *geom);

extern void check_rmin(void);

extern void transgcov(FTYPE *gcov, FTYPE (*trans)[NDIM], FTYPE *gcovprim);

extern void transgcovself(FTYPE *gcov, FTYPE (*trans)[NDIM]);

extern void get_gcovpert(FTYPE *gcovprim, FTYPE *gcovpert, FTYPE *gcovpertprim);

extern void transgcovgcovpertself(FTYPE *gcov, FTYPE *gcovpert, FTYPE (*trans)[NDIM]);
