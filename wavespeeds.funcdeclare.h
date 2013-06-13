
extern int get_global_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *pr,FTYPE *wspeedtemp);
extern int get_global_wavespeeds_full(int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*wspeed)[NUMCS][NSTORE1][NSTORE2][NSTORE3]);
extern int global_vchar(FTYPE (*pointspeed)[NSTORE2][NSTORE3][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*wspeed)[NUMCS][NSTORE1][NSTORE2][NSTORE3]);
extern int get_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, struct of_state *state_l, struct of_state * state_r, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE *ctop);
extern int vchar(FTYPE *pr, struct of_state *q, int dir,
                 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE (*ecov)[NDIM],
                             FTYPE (*econ)[NDIM]);
extern void transform(FTYPE *vec, FTYPE (*t)[NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE (*t)[NDIM]);
