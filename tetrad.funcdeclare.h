
/*! \file tetrad.funcdeclare.h
    \brief function declarations for global use of tetrad.c functions
*/



extern void make_tetr(FTYPE *ucon, FTYPE (*econ)[NDIM]);
extern int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]);
extern int tetr_func_frommetric(int primcoord, FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);

extern int calc_ORTHOes(int primcoord, struct of_geom *ptrgeom, FTYPE tmuup[][NDIM], FTYPE tmudn[][NDIM]);


extern int calc_generalized_boost_uu(struct of_geom *ptrgeom, FTYPE *wcon, FTYPE *ucon, FTYPE (*lambda)[NDIM]);
extern int calc_ortho_boost_uu(FTYPE *wcon, FTYPE *ucon, FTYPE (*lambda)[NDIM]);


extern int transboost_lab2fluid(int lab2orthofluid, int primcoord, struct of_geom *ptrgeom, FTYPE *uconlab, FTYPE (*transboostup)[NDIM], FTYPE (*transboostlo)[NDIM]);

extern int vector_lab2orthofluidorback(int primcoord, int lab2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, FTYPE v4concovtype, FTYPE *vector4in, FTYPE *vector4out);

extern int tensor_lab2orthofluidorback(int primcoord, int lab2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, int tconcovtypeA, int tconcovtypeB, FTYPE (*tensor4in)[NDIM], FTYPE (*tensor4out)[NDIM]);

extern void vecX2vecVortho(int concovtype, struct of_geom *ptrgeom, FTYPE *vec, FTYPE *vecortho);


extern int vector_harm2orthofluidorback(int whichvector, int harm2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, FTYPE v4concovtype, FTYPE *vector4in, FTYPE *vector4out);





// unused functions:
extern int calc_ZAMOes_old(struct of_geom *ptrgeom, FTYPE tmuup[][NDIM], FTYPE tmudn[][NDIM]);
