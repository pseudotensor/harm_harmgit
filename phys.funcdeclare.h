/*! \file phys.funcdeclare.h
     \brief Function declarations for global use of things in phys.c
*/

extern int sourcephysics(FTYPE *pi, FTYPE *ph, FTYPE *pf, int *didreturnpf, int *eomtype, struct of_geom *geom, struct of_state *q, FTYPE *Ugeomfreei, FTYPE *Ugeomfreef, FTYPE* CUf, FTYPE *CUimp, FTYPE dissmeasure, FTYPE *dUother, FTYPE (*dUcomp)[NPR]);

extern void postdt(void);
extern int primtoU(int returntype, FTYPE *p, struct of_state *q, struct of_geom *geom,
                   FTYPE *U, FTYPE *Uabs);

extern int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_4vel_bothut(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *ucon2, FTYPE *others);

extern int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);

extern int uconrel(FTYPE *ucon, FTYPE *uconrel, struct of_geom *geom);


#if(RELTYPE==RELEOM)

#if(WHICHVEL==VEL4)
#define ucon_calc ucon_calc_4vel
#define dudp_calc dudp_calc_gen
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#elif(WHICHVEL==VEL3)
#define ucon_calc ucon_calc_3vel
#define dudp_calc dudp_calc_3vel
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#elif(WHICHVEL==VELREL4)
#define ucon_calc ucon_calc_rel4vel
#define dudp_calc dudp_calc_gen
#define compute_1plusud0 compute_1plusud0_rel4vel // uses qsq and gamma from ucon_calc_rel4vel()
#define bsq_calc bsq_calc_rel4vel
#define bsq_calc_fromq bsq_calc_fromq_rel4vel
#elif(RELTYPE==NONRELEOM) // not really right
#define ucon_calc ucon_calc_nonrel
#define dudp_calc dudp_calc_nonrel
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#endif
#endif


extern int ucon_calcother(FTYPE *pr, FTYPE *ucon, FTYPE *others);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,FTYPE *CC, FTYPE *discr);


extern int ucon_calc_whichvel(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);


// physics stuff
extern int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
extern int set_zamo_ucovuconplus1ud0(struct of_geom *ptrgeom, FTYPE *ucov, FTYPE *ucon, FTYPE *plus1ud0);
extern int set_zamo_ucon(struct of_geom *ptrgeom, FTYPE *ucon);

extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern int bsq_calc_fromq(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *b2);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);
extern void bsq_calc_rel4vel_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


extern int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE *gamma, FTYPE *qsq);


extern int dudp_calc_gen(int whicheos, int whichcons, FTYPE *EOSextra, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE **alpha);

extern int dudp_calc_3vel(int whicheos, int whichcons, FTYPE *EOSextra, FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha);


extern int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);

extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
extern void UtoU_evolve2diag(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);


// presume get_geometry() only necessarily feeds back pointer where geometry is located
// if want hard copy to be created using a memory, should use copy_geometry()
// still required to have pointer point to physical allocated memory in general
#if(NEWMETRICSTORAGE)
// overwrites any prior pointer reference
// better than copying since no point in copying if memory already exists and get_geometry() is always for real computational geometry
// This only works if feed pointer directly into get_geometry.  Can't feed &geom.

#if(MCOORD!=CARTMINKMETRIC)

#define get_geometry(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(compgeom,pp,ii,jj,kk);
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP0A1(gdetgeom,ii,jj,kk,pp);
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,ii,jj,kk);
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,ii,jj,kk);
#define set_igdet(arg)
#define set_igdetsimple(arg)

#else// else if CARTMINKMETRIC

#define get_geometry(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(compgeom,pp,0,0,0);
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP0A1(gdetgeom,0,0,0,pp);
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,0,0,0);
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,0,0,0);
#define set_igdet(arg)
#define set_igdetsimple(arg)

#endif // end if CARTMINKMETRIC



#else // else if old metric storage method


#define get_geometry(ii,jj,kk,pp,ptrgeom) get_geometry_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) get_geometry_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) get_geometry_gdetonly_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) get_geometry_geomeonly_old(ii,jj,kk,pp,ptrgeom)
#define set_igdet(arg) set_igdet_old(arg)
#define set_igdetsimple(arg) set_igdetsimple_old(arg)


#endif



#if(NEWMETRICSTORAGE==0)
extern void get_geometry_old(int i, int j, int k, int loc, struct of_geom *geom);
extern void get_geometry_gdetonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void get_geometry_geomeonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void set_igdetsimple_old(struct of_geom *geom);
#endif






extern void get_allgeometry(int i, int j, int k, int loc, struct of_allgeom *allgeom, struct of_geom *geom);



extern int get_state(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_state_norad_part1(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_state_norad_part2(int needentropy, FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_state_radonly(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_state_nofield(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_stateforcheckinversion(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
extern int get_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **q);
extern int get_stateforfluxcalc(int dimen, int isleftright, FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);
extern int get_stateforUdiss(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
extern int get_stateforinterpline(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);
extern int get_stateforglobalwavespeeds(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);


extern void compute_and_store_fluxstatecent(FTYPE (*pr)[NSTORE2][NSTORE3][NPR]);


extern int primtoflux(int returntype, FTYPE *pa, struct of_state *q, int dir,
                      struct of_geom *geom, FTYPE *fl, FTYPE *flabs);
extern int primtoflux_splitmaem(int returntype, FTYPE *pa, struct of_state *q, int fluxdir, int fundir, struct of_geom *geom, FTYPE *flma, FTYPE *flem);


extern int flux_compute_general(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *ctopall);
extern int flux_compute_splitmaem(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *FEM, FTYPE *ctopall);

extern void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom,
                            struct of_state *q,FTYPE *dU);
extern int source(FTYPE *pi, FTYPE *pa,  FTYPE *pf, int *didreturnpf, int *eomtype, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *CUf, FTYPE *CUimp, FTYPE dissmeasure, FTYPE *dUriemann,
                  FTYPE (*Uacomp)[NPR], FTYPE *Ua);

extern FTYPE taper_func(FTYPE R,FTYPE rin) ;

extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
extern void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE (*faraday)[NDIM]);
extern void current_precalc(int which, struct of_geom *geom, struct of_state *q, SFTYPE Dt,FTYPE (*faraday)[3]);
extern void current_calc(FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3]);
extern int current_doprecalc(int which, FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
