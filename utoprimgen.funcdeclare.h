
extern int Utoprimgen(int finalstep, int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprimloop(FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR], struct of_newtonstats *newtonstats);
extern int primtoUloop(FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR]);

extern int Utoprim(int entropyeom, FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_1d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_1d_opt(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_2d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_1d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_2d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);
//extern int Utoprim_2d_final_nonrelcompat_inputnorestmass(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);  //wrong function name, corrected by atch, see below
extern int Utoprim_jon_nonrelcompat_inputnorestmass(int eomtype, FTYPE *EOSextra, FTYPE *U, struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats);
extern int Utoprim_5d2_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats);

extern int Utoprimdiss(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, PFTYPE *otherfail, struct of_newtonstats *newtonstats);

/* // dudp stuff */

/* extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui); */
/* extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon, */
/*                      FTYPE *b, FTYPE (*dbiduj)[NDIM]); */
/* extern void db2dui_calc(FTYPE (*dbiduj)[NDIM], FTYPE *b, */
/*                      FTYPE *db2dui); */
/* extern void duudud_calc(FTYPE *ucon, FTYPE (*duudud)[NDIM]); */

/* extern void dbsqdui_calc(FTYPE (*dbiduj)[NDIM], FTYPE *b, */
/*                       FTYPE *dbsqdui); */
/* extern void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi); */
/* extern void duidvj_calc(FTYPE *dgdv,FTYPE (*duidvj)[NDIM]); */
/* extern void dudduu_calc(FTYPE*dutdui, FTYPE (*dudduu)[NDIM]); */
/* extern void dbdiduj_calc(FTYPE (*dbiduj)[NDIM],FTYPE (*dbdiduj)[NDIM]); */
/* extern void ducon_dv3_calc(struct of_state *q,FTYPE (*ducon_dv)[NDIM]); */
extern int sp_stress_calc(FTYPE *pr, FTYPE (*tens_matt)[NDIM],
                          FTYPE (*tens_em)[NDIM], FTYPE *b,
                          FTYPE *ucon);
