// some useful wrappers for user simplicity:
extern FTYPE pressure_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE cs2_compute_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE compute_entropy_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern void get_EOS_parms_simple(int*numparms, int i, int j, int k, int loc, FTYPE *parlist);
extern FTYPE compute_temp_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE u);
extern int get_extrasprocessed_simple(int doall, int i, int j, int k, int loc, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_u_from_entropy_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE entropy);
extern FTYPE compute_qdot_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdrho0_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdu_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);

// other wrappers
extern int ufromentropy_calc(struct of_geom *ptrgeom, FTYPE entropy, FTYPE *pr);
extern int entropy_calc(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy);
extern int invertentropyflux_calc(struct of_geom *ptrgeom, FTYPE entropyflux,int dir, struct of_state *q, FTYPE*pr);


// for old inversion methods:
extern FTYPE pressure_rho0_w(FTYPE *EOSextra, FTYPE rho0, FTYPE w);

// eos stuff [should be consistent with defs.general.h ptr's and eos.c assignment to pointers and eos.c,idealgaseos.c,etc. versions of those functions.]
extern int pickeos_eomtype(int whicheos, int whicheom);

extern FTYPE pressure_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_u_from_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy);
extern FTYPE u_rho0_p(FTYPE *EOSextra, FTYPE rho0, FTYPE p);
extern FTYPE dpdu_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE dpdrho0_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE cs2_compute(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdrho(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdu(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE pressure_wmrho0(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_qdot(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern int compute_sources_EOS(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR]);
extern void compute_allextras(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u,int *numextrasreturn,FTYPE*extras);
extern int get_extrasprocessed(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_temp(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern void compute_EOS_parms(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void store_EOS_parms(int numparms, FTYPE *EOSextra, FTYPE *parlist);
extern void get_EOS_parms(int*numparms, FTYPE *EOSextra, FTYPE *parlist);
