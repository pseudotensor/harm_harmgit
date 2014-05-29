
/*! \file global.inits.h
  \brief Function declarations (used globally) from init.tools.c
 */

extern int setRin_withchecks(FTYPE *rin);
extern int user1_prepre_init_specific_init(void);
extern int user1_init_conservatives(int *fieldfrompotential, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*Utemp)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR]);
extern int user1_post_init_specific_init(void);
extern int user1_init_global(void);
extern int user1_init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr);

extern int user1_init_primitives(int inittype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);


extern int user1_init_vpot2field_user(SFTYPE time, int *fieldfrompotential, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
extern int user1_normalize_densities(int eqline, FTYPE *parms, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *rhomax, FTYPE *umax);
extern int user1_getmax_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax);
extern int user1_getmax_densities_full(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax, SFTYPE *uradmax, SFTYPE *utotmax, SFTYPE *pmax, SFTYPE *pradmax, SFTYPE *ptotmax);
extern int user1_get_maxes(int eqslice, FTYPE *parms, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min);
extern int user1_normalize_field(FTYPE beta, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
extern int user1_set_atmosphere(int atmospheretype, int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr);

extern int user1_normalize_field_sigma(FTYPE sigma0, FTYPE bpole, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
extern int user1_normalize_densities_postnormalizefield(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int user1_get_sigmabsq_atpole(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *sigma_pole, FTYPE *bsq_pole);




extern int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);



extern int post_init(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*faraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3]);


extern int coolfunc_user(FTYPE h_over_r, FTYPE *pr, struct of_geom *geom, struct of_state *q,FTYPE (*dUcomp)[NPR]);
