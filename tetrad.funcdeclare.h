extern void make_tetr(FTYPE *ucon, FTYPE (*econ)[NDIM]);
extern int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]);
extern int tetr_func_frommetric(FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);
