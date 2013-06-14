// fixup stuff

extern int check_pr(FTYPE *pr, FTYPE *prmodel, FTYPE *ucons, struct of_geom *geom, int modelpos,int finalstep);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
                    FTYPE *ucon);


extern int set_atmosphere(int whichcond, int whichvel, struct of_geom *geom, FTYPE *pr);
extern int set_density_floors_default(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int fixup(int stage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int finalstep);
extern int fixup1zone(FTYPE *pr,FTYPE *ucons, struct of_geom *ptrlgeom, int finalstep);

extern int diag_fixup(int docorrectucons, FTYPE *pr0, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep,int whocalled);
int diag_fixup_allzones(int truestep, int finalstep, FTYPE (*pf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);

extern int diag_fixup_Ui_pf(int docorrectucons, FTYPE *Ui, FTYPE *pf, struct of_geom *ptrgeom, int finalstep, int whocalled);
extern int diag_fixup_U(int docorrectucons, FTYPE *Ui, FTYPE *Uf, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep,int whocalled);


extern int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);

extern int limit_gamma(FTYPE gammamax, FTYPE*pr, FTYPE *ucons, struct of_geom *geom, int finalstep);

extern int fixup_checksolution(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],int finalstep);
extern int fixup_utoprim(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int fixup_utoprim_nofixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

extern int post_fixup(int stage, int finalstep, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern int post_fixup_nofixup(int stageit, int finalstep, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);

extern int pre_fixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR]);
extern int get_bsqflags(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR]);

extern int inflow_check_4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_3vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_rel4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);

extern void diag_eosfaillookup(int i, int j, int k);


