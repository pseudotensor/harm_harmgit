
/*! \file global.funcdeclare.rad.h
  \brief Function declarations (used globally) from RADIATION parts of code
 */

extern void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM]);
extern void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *radstressdir, FTYPE *radstressdirabs);
extern FTYPE my_min(FTYPE a, FTYPE b);
extern FTYPE my_sign(FTYPE x);
extern int inverse_44matrix(FTYPE a[][4], FTYPE ia[][4]);
extern int boost22_fforzamo(int whichdir, FTYPE T1[][4],FTYPE T2[][4],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int boost22_zamo2ff(FTYPE T1[][4],FTYPE T2[][4],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int boost22_ff2zamo(FTYPE T1[][4],FTYPE T2[][4],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int trans22_zamo2lab(FTYPE T1[][4],FTYPE T2[][4],FTYPE elo[][4]);
extern int trans22_lab2zamo(FTYPE T1[][4],FTYPE T2[][4],FTYPE eup[][4]);
extern int trans2_lab2zamo(FTYPE *u1,FTYPE *u2,FTYPE eup[][4]);
extern int trans2_zamo2lab(FTYPE *u1,FTYPE *u2,FTYPE elo[][4]);

extern int boost2_zamo2ff(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int boost2_ff2zamo(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][4]);
extern int boost2_fforzamo(int whichdir, FTYPE A1[4],FTYPE A2[4],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][4]);

extern int calc_Rij(FTYPE *pp, FTYPE Rij[][4]);

extern int indices_2221(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom);
extern int indices_2122(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom);
extern int indices_2212(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom);
extern int indices_21(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom);
extern int indices_12(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom);


extern void koral_source_rad_calc(int computestate, int computeentropy, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Gdpl, FTYPE *Gdabspl, FTYPE *chi, FTYPE *Tgas, FTYPE *Trad, struct of_state *q);
extern int calc_rad_lambda(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q, FTYPE Tgas, FTYPE *lambda, FTYPE *nlambda, FTYPE *kappaemit, FTYPE *kappanemit);


extern void calc_Tandopacityandemission(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE Ruu, FTYPE gammaradgas, FTYPE B, FTYPE *Tgas, FTYPE *Tradff, FTYPE *nradff, FTYPE *kappa, FTYPE *kappan, FTYPE *kappaemit, FTYPE *kappanemit, FTYPE *kappaes, FTYPE *lambda, FTYPE *nlambda);




extern int prad_fforlab(int *whichvel, int *whichcoord, int whichdir, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout);
extern int prad_labtoff(int *whichvel, int *whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout);
extern int prad_fftolab(int *whichvel, int *whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout);

extern int whichfluid_ffrad_to_primeall(int *whichvel, int *whichcoordfluid, int *whichcoordrad, struct of_geom *ptrgeomprimecoords, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout);

extern int primefluid_EVrad_to_primeall(int *whichvel, int *whichcoord, struct of_geom *ptrgeom, FTYPE *pin, FTYPE *pout);


extern int prad_ff2zamo(FTYPE *pp1, FTYPE *pp2, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int f_prad_zamo2ff(FTYPE *ppff, FTYPE *ppzamo, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4],FTYPE *f);
extern int prad_zamo2ff(FTYPE *ppzamo, FTYPE *ppff, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);

extern int u2p_rad(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, int whichcap, FTYPE *uu, FTYPE *pp, struct of_geom *ptrgeom, PFTYPE *lpflag, PFTYPE *lpflagrad);


extern int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);

extern int vchar_rad(FTYPE *pr, struct of_state *q, int dir,
                     struct of_geom *geom, FTYPE *cmax, FTYPE *cmin, FTYPE *cmax2, FTYPE *cmin2,int *ignorecourant);


extern void calc_kappa(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *kappa);
extern void calc_kappaemit(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappaemit);
extern void calc_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappaes);
extern FTYPE calc_kappa_user(FTYPE rho, FTYPE B, FTYPE Tg, FTYPE Tr,FTYPE x,FTYPE y,FTYPE z);
extern FTYPE calc_kappan_user(FTYPE rho, FTYPE B, FTYPE Tg, FTYPE Tr,FTYPE x,FTYPE y,FTYPE z);
extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
extern int calcfull_tautot(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *tautot, FTYPE *tautotmax);
extern int calc_tautot(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q, FTYPE *tautot, FTYPE *tautotmax);


extern FTYPE calc_LTE_NfromT(FTYPE T);
extern FTYPE calc_LTE_NfromE(FTYPE E);


extern FTYPE calc_LTE_EfromT(FTYPE);
extern FTYPE calc_LTE_TfromE(FTYPE);
extern FTYPE calc_LTE_Efromurho(FTYPE E,FTYPE);
extern FTYPE calc_PEQ_ufromTrho(FTYPE,FTYPE);
extern FTYPE calc_PEQ_Tfromurho(FTYPE,FTYPE);
extern FTYPE calc_LTE_EfromT(FTYPE);
extern FTYPE calc_LTE_TfromE(FTYPE);
extern FTYPE calc_LTE_Efromurho(FTYPE E,FTYPE);
extern FTYPE calc_PEQ_ufromTrho(FTYPE,FTYPE);
extern FTYPE calc_PEQ_Tfromurho(FTYPE,FTYPE);







extern int set_ncon_velocity(int whichvel, FTYPE gammamax, FTYPE *ncon, struct of_geom *ptrgeom, FTYPE *uconwhichvel);

