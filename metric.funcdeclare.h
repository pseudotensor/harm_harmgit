
/*! \file metric.funcdeclare.h
     \brief Function delcarations for global use of functions in metric.c
*/


// metric stuff
extern int metric_checks(struct of_geom *ptrgeom);

//extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
//extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE *gcov);
//extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE *gcon);
extern void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
                      FTYPE (*lconn)[NDIM][NDIM],FTYPE *conn2);
extern void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol);

extern void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
extern void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon);


extern int rotate_VtoVmetric(int whichcoord, FTYPE ROTANGLE, FTYPE *V, FTYPE *Vmetric);
extern int rotate_VmetrictoV(int whichcoord, FTYPE ROTANGLE, FTYPE *Vmetric, FTYPE *V);

//extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE *gcov);
//extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon);

extern int fix_hp(FTYPE *h, FTYPE *p);

