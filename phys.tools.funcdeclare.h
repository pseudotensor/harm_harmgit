
/*! \file phys.tools.funcdeclare.h
     \brief Function declarations for global use of things in phys.tools.c
     */

extern void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
extern void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);

extern void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
extern void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);

extern int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
                      FTYPE *bcon);
extern void lower_vec(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void lowerf(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise_vec(FTYPE *v1, struct of_geom *geom, FTYPE *v2);



extern int limit_3vel_ffde(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr);


extern int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf);

extern int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf, FTYPE *dualfabs);


