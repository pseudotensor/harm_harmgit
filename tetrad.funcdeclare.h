extern void make_tetr(FTYPE *ucon, FTYPE (*econ)[NDIM]);
extern int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]);
extern int tetr_func_frommetric(FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);
extern void vecX2vecVortho(int concovtype, FTYPE V[],  FTYPE *gcov,  FTYPE (*dxdxp)[NDIM], FTYPE (*tetrcov)[NDIM], FTYPE (*tetrcon)[NDIM], FTYPE *vec, FTYPE *vecortho);
extern int calc_LNRFes(struct of_geom *ptrgeom, FTYPE emuup[][NDIM], FTYPE emulo[][NDIM]);

extern int calc_LNRFes_old(struct of_geom *ptrgeom, FTYPE emuup[][NDIM], FTYPE emulo[][NDIM]);

