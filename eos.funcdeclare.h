
/*! \file eos.funcdeclare.h
  \brief Equation of State function declarations
*/


/////
//
// some useful wrappers for user simplicity:
//
////
extern FTYPE pressure_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE pressure_rho0_u_simple_forcheckinversion(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE u_rho0_p_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE p);
extern FTYPE u_rho0_T_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE T);
extern FTYPE cs2_compute_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE compute_entropy_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE compute_entropy_simple_forcheckinversion(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern void get_EOS_parms_simple(int*numparms, int i, int j, int k, int loc, FTYPE *parlist);
extern void fix_primitive_eos_scalars_simple(int i, int j, int k, int loc, FTYPE *pr);


extern void yl2advect_kazfull(FTYPE *EOSextra, FTYPE pryl, FTYPE prynu, FTYPE *prforadvect);
extern void ynu2advect_kazfull(FTYPE *EOSextra, FTYPE pryl, FTYPE prynu, FTYPE *prforadvect);

extern void advect2yl_kazfull(FTYPE *EOSextra, FTYPE ylforadvect, FTYPE ynuforadvect, FTYPE *ye);
extern void advect2ynu_kazfull(FTYPE *EOSextra, FTYPE ylforadvect, FTYPE ynuforadvect,FTYPE *prynu);


extern FTYPE compute_temp_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE u);
extern int get_extrasprocessed_simple(int doall, int i, int j, int k, int loc, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_u_from_entropy_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE entropy);
extern FTYPE compute_qdot_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdrho0_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdu_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);


/////
//
// other wrappers
//
/////
extern int ufromentropy_calc(struct of_geom *ptrgeom, FTYPE entropy, FTYPE *pr);
extern int entropy_calc(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy);
extern int entropy_calc_forcheckinversion(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy);
extern int invertentropyflux_calc(struct of_geom *ptrgeom, FTYPE entropyflux,int dir, struct of_state *q, FTYPE*pr);


// for old inversion methods:
extern FTYPE pressure_rho0_w(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE w);

////
//
// eos stuff [should be consistent with defs.general.h ptr's and eos.c assignment to pointers and eos.c,idealgaseos.c,etc. versions of those functions.]
//
////
extern int pickeos_eomtype(int whicheosinput, int whicheom, int *whicheosoutput);
extern int initeos_eomtype(void);


extern FTYPE pressure_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_u_from_entropy(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE entropy);
extern FTYPE u_rho0_p(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE p);
extern FTYPE u_rho0_T(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE T);
extern FTYPE dpdu_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE dpdrho0_rho0_u(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE cs2_compute(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_entropy(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdrho(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdu(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_specificentropy_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dspecificSdrho_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_dspecificSdwmrho0_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);

extern FTYPE pressure_wmrho0(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_qdot(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern int compute_sources_EOS(int whicheos, FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR]);
extern void compute_allextras(int whicheos, int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u,int *numextrasreturn,FTYPE*extras);
extern int get_extrasprocessed(int whicheos, int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_temp(int whicheos, FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern void compute_EOS_parms(int whicheos, FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void compute_EOS_parms_full(int whicheos, FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void store_EOS_parms(int whicheos, int numparms, FTYPE *EOSextra, FTYPE *parlist);
extern void get_EOS_parms(int whicheos, int*numparms, FTYPE *EOSextra, FTYPE *parlist);
extern void fix_primitive_eos_scalars(int whicheos, FTYPE *EOSextra, FTYPE *pr);
extern void getall_forinversion(int whicheos, int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu);



///
//
// SPECIAL KAZ EOS GLOBAL FUNCTIONS:
//
///
extern void initeos_kazfulleos(void);
extern void read_setup_eostable(void);

