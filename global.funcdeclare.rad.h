
extern void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *radstressdir);
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
extern int indices_21(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom);
extern int indices_12(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom);

extern int calc_LNRFes(struct of_geom *ptrgeom, FTYPE emuup[][4], FTYPE emulo[][4]);


extern int p2u_rad(FTYPE *pr, FTYPE *Urad, struct of_geom *ptrgeom, struct of_state *q);


extern int prad_ff2zamo(FTYPE *pp1, FTYPE *pp2, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);
extern int f_prad_zamo2ff(FTYPE *ppff, FTYPE *ppzamo, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4],FTYPE *f);
extern int prad_zamo2ff(FTYPE *ppzamo, FTYPE *ppff, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][4]);

extern int u2p_rad(FTYPE *uu, FTYPE *pp, struct of_geom *ptrgeom);


extern int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);

extern int vchar_rad(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern int vchar_all(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern int vchar_each(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmaxmhd, FTYPE *cminmhd, FTYPE *cmaxrad, FTYPE *cminrad,int *ignorecourant);


#if(0)
//misc.c
int calc_stationary1d_solution()  ;
FTYPE calc_totalmass();
int initialize_arrays();
int free_arrays();
FTYPE my_min(FTYPE a, FTYPE b);
//FTYPE my_max(FTYPE a, FTYPE b);
FTYPE my_sign(FTYPE);
int find_eigenvalues3(FTYPE[3][3],FTYPE*);
int find_eigenvalues(FTYPE*,int);
FTYPE find_max_eigenvalue(FTYPE *data,int N);
int find_max_eigenvalue_lr(FTYPE *data,int N,FTYPE*,FTYPE*);
int my_err(char *);
int inverse_44matrix(FTYPE a[][4], FTYPE ia[][4]);
int convert_out2gif_1d(char *fname,char*,int niter,FTYPE t);
int convert_out2gif_2d(char *fname,char*,int niter,FTYPE t);
int getch(void);
int dosthelse(void);
//finite.c
FTYPE f_calc_fluxes_at_faces(int ix,int iy,int iz);

int f_timeder (FTYPE t, FTYPE dt,FTYPE, FTYPE*, int ifcopy2u0, FTYPE*);
int set_grid(FTYPE*, FTYPE*, FTYPE*,FTYPE*);
int print_grid(FTYPE,FTYPE,FTYPE);
FTYPE fd_flux_limiter(FTYPE r);
FTYPE minmod_fd_flux_limiter(FTYPE ,FTYPE,FTYPE);
FTYPE f_der_kurganovtadmor(int ix,int iy, int yz,FTYPE*);
FTYPE f_der_hlle_obsolete(int ix,int iy, int yz,FTYPE*);
FTYPE f_der_muscl(int ix,int iy, int yz,FTYPE*);
//FTYPE get_x(int,int);
#endif

#if(0)
#define get_x(ic,idim) (idim==0 ? x[ic+NG] : (idim==1 ? x[ic+NG + NX+2*NG] : (idim==2 ? x[ic+NG + NX+2*NG + NY+2*NG ] : 0.)))
//FTYPE get_xb(int,int);
#define get_xb(ic,idim) (idim==0 ? xb[ic+NG] : (idim==1 ? xb[ic+NG + NX+2*NG + 1] : (idim==2 ? xb[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)))
int set_x(int,int,FTYPE);
int set_xb(int,int,FTYPE);
//FTYPE get_u(FTYPE*,int,int,int,int);
#define get_cflag(iflag,ix,iy,iz) (cellflag[iflag + (ix+NG)*NFLAGS + (iy+NG)*(NX+2*NG)*NFLAGS + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NFLAGS])
#define set_cflag(iflag,ix,iy,iz,val) cellflag[iflag + (ix+NG)*NFLAGS + (iy+NG)*(NX+2*NG)*NFLAGS + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NFLAGS]=val
#define get_u(uarr,iv,ix,iy,iz) (uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV])
//int set_u(FTYPE*,int,int,int,int,FTYPE);
#define set_u(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val
FTYPE get_u_scalar(FTYPE*,int,int,int);
//int set_u_scalar(FTYPE*,int,int,int,FTYPE);
#define set_u_scalar(uarr,ix,iy,iz,val) uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)] = val
int copy_u(FTYPE,FTYPE*,FTYPE*);
int add_u(FTYPE f1, FTYPE* u1, FTYPE f2, FTYPE *u2, FTYPE *u3);
#endif

#if(0)
FTYPE f_timeder_source_term(FTYPE t, const FTYPE y[], FTYPE f[],  void *params);
int if_indomain(int ix,int iy,int iz);
int if_outsidegc(int ix,int iy,int iz);
int if_outsidewave(int ix,int iy,int iz);
#endif


#if(0)
FTYPE get_size_x(int ic, int idim);
#endif

#if(0)
int set_ub(FTYPE* uarr,int iv,int ix,int iy,int iz,FTYPE value,int idim);
//#define set_ub(uarr,iv,ix,iy,iz,idim,val) (idim==0 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV]=val : (idim==1 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV]=val : (idim==2 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val : 0.)))
#define set_ubx(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV]=val
#define set_uby(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV]=val
#define set_ubz(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val
//FTYPE get_ub(FTYPE* uarr,int iv,int ix,int iy,int iz,int idim);
#define get_ub(uarr,iv,ix,iy,iz,idim) (idim==0 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV] : (idim==1 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV] : (idim==2 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV] : 0.)))
//FTYPE get_g(FTYPE* uarr, int i,int j, int ix, int iy, int iz);
#define get_g(uarr,i,j,ix,iy,iz) uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE]
int set_g(FTYPE* uarr,int i,int j,int ix,int iy,int iz,FTYPE value);
int set_T(FTYPE* uarr,int i,int j,int ix,int iy,int iz,FTYPE value);

#define get_T(uarr,i,j,ix,iy,iz) uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*16]
#define get_Tb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG+1)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*16] : (idim==1 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*16] : (idim==2 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*16] : 0.)))
int set_Tb(FTYPE* uarr,int i,int j,int ix,int iy,int iz,FTYPE value,int idim);
int set_gb(FTYPE* uarr,int i,int j,int ix,int iy,int iz,FTYPE value,int idim);
//FTYPE get_gb(FTYPE* uarr,int i,int j,int ix,int iy,int iz,int idim);
#define get_gb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG+1)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*gSIZE] : (idim==1 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*gSIZE] : (idim==2 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE] : 0.)))
#define get_gKr(i,j,k,ix,iy,iz) gKr[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*64]
#define get_gKrb(i,j,k,ix,iy,iz,idim) (idim==0 ? gKrbx[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG+1)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*64] : (idim==1 ? gKrby[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG+1)*(NX+2*NG+1)*64] : (idim==2 ? gKrbz[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*64] : 0.)))
#define set_gKr(i,j,k,ix,iy,iz,val) gKr[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*64]=val
int set_Krb(int i,int j,int k,int ix,int iy,int iz,FTYPE value,int idim);
#endif


#if(0)
//fileop.c
int fread_restartfile(FTYPE*);
int fprint_openfiles();
int fprint_closefiles();
int fprint_profiles(FTYPE,FTYPE);
int print_profiles();
#endif

#if(0)
FTYPE max_eigen_Jac(FTYPE *,FTYPE*,int,void*);
int calc_wavespeeds(int,int,int,FTYPE*,FTYPE*,FTYPE*,FTYPE*,FTYPE*,FTYPE*);
int calc_wavespeeds_lr(int,int,int,FTYPE*);
int calc_wavespeeds_lr_faces( int,int,int,int, FTYPE*,FTYPE*);
int max_eigen_lr_Jac(FTYPE *,FTYPE*,int,void*,FTYPE*,FTYPE*);
int calc_Jac_num(FTYPE *xx,FTYPE *ujac, int idim,void *parameters,FTYPE *fd_jac);
int set_initial_profile();
FTYPE f_flux_prime(FTYPE *uu, int,int,int,int,FTYPE *ff);
FTYPE f_diffusion_prime(FTYPE *uu, FTYPE *du, int iv,int,void*);
FTYPE f_grav_potential(FTYPE,FTYPE,FTYPE);
FTYPE f_der_grav_potential(FTYPE,FTYPE,FTYPE,int);
int f_source_term(int,int,int,FTYPE *);
int f_fourforce_source_term(int,int,int,FTYPE *);
int f_implicit_rhou(int ix, int iy, int iz, FTYPE *rho, FTYPE *uint, FTYPE dt);
int f_metric_source_term(int ix, int iy, int iz,FTYPE *ss);
int f_metric_source_term_face(int ix, int iy, int iz,int idim,int ifleft,FTYPE *ss);
int f_source_term_face(int ix, int iy, int iz,int idim,int ifleft,FTYPE *ss);
int initialize_problem();
int calc_maxgravspeed(int ix, int iy, int iz, FTYPE dt);
int calc_maxwavespeed_obsolete(int ix, int iy, int iz,void*);
int calc_maxwavespeed_osbolete_cs(int ix, int iy, int iz,void*);
FTYPE calc_chi(FTYPE rho, FTYPE T);
int rad_fld_factors(int ix,int iy, int iz,FTYPE *,void *rad_param,int);
int rad_fld_factors_arb(int ix,int iy, int iz,FTYPE*,void *rad_param);

int radx_flux_implicit(FTYPE *uu,FTYPE dtt);
int LTE_implicit(FTYPE *uu,FTYPE dtt);
FTYPE max_eigen_radx(FTYPE nx,FTYPE ny,FTYPE nz,int idim);
int max_eigen_lr_radx(FTYPE nx,FTYPE ny,FTYPE nz,int idim,FTYPE *al, FTYPE *ar);
#endif

int gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3);


#if(0)
FTYPE f_der_hlle         (int ix,int iy,int iz, FTYPE *fd_der);
FTYPE calc_kappa(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
FTYPE calc_kappaes(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
FTYPE calc_ufromS(FTYPE S,FTYPE rho);
FTYPE calc_Sfromu(FTYPE S,FTYPE u);
int
avg2point(FTYPE *um2,FTYPE *um1,FTYPE *u0,FTYPE *up1,FTYPE *up2,FTYPE*,FTYPE*,FTYPE dxm2,FTYPE dxm1,FTYPE dx0,FTYPE dxp1,FTYPE dxp2);

//problem.c
FTYPE calc_xb(int i,int idim);
int calc_bc(int,int,int,FTYPE,FTYPE*, FTYPE*);
int pr_tophat_inside(FTYPE x,FTYPE y,FTYPE z);
int my_finger(FTYPE);

//metric.c
int calc_LNRFes(FTYPE g[][5], FTYPE emuup[][4], FTYPE emulo[][4]);
int calc_g(FTYPE*,FTYPE[][5]);
int calc_G(FTYPE*,FTYPE[][5]);
int calc_Krzysie(FTYPE*,FTYPE[][4][4]);
int print_Krzysie(FTYPE g[][4][4]);
int print_g(FTYPE [][5]);
FTYPE calc_gdet(FTYPE *xx);
FTYPE calc_dlgdet(FTYPE *xx, int idim);
#endif


#if(0)
int pick_g(int ix,int iy,int iz,FTYPE gg[][5]);
int pick_G(int ix,int iy,int iz,FTYPE gg[][5]);
int pick_gb(int ix,int iy,int iz,int,FTYPE gg[][5]);
int pick_Gb(int ix,int iy,int iz,int,FTYPE gg[][5]);
int pick_T(FTYPE *arr,int ix,int iy,int iz,FTYPE T[][4]);
int pick_Tb(FTYPE *arr,int ix,int iy,int iz,int,FTYPE T[][4]);

int p2u_Sonly(FTYPE *p, FTYPE *u,FTYPE[][5]);
int print_p(FTYPE *p);
int print_u(FTYPE *p);
int calc_Tmunu( FTYPE *p, FTYPE g[][5], FTYPE T[][4],FTYPE*);
int convert_uold2urel(FTYPE *x,FTYPE *u);
int calc_metric();
int calc_sourceterms(int,int,int);
int calc_primitives(int,int,int);
int calc_conserved(int ix,int iy,int iz);

FTYPE r_horizon_BL(FTYPE a);
FTYPE r_mbound_BL(FTYPE a);
FTYPE r_photon_BL(FTYPE a);
int update_entropy(int ix,int iy,int iz,int u2pflag);

//u2p.c
int u2p(FTYPE *uu, FTYPE *pp, FTYPE gg[][5], FTYPE GG[][5], FTYPE eup[][4], FTYPE elo[][4]);
int u2p_hot(FTYPE*,FTYPE*,FTYPE[][5]);
int u2p_entropy(FTYPE*,FTYPE*,FTYPE[][5]);
int u2p_cold(FTYPE*,FTYPE*,FTYPE[][5]);
int
u2p_rad_num(FTYPE *uu, FTYPE *pp, FTYPE gg[][5], FTYPE eup[][4], FTYPE elo[][4]);
int
u2p_rad(FTYPE *uu, FTYPE *pp, FTYPE gg[][5], FTYPE GG[][5], FTYPE eup[][4], FTYPE elo[][4]);

int
dump_u2p_rad(FTYPE *uu, FTYPE *pp, FTYPE gg[][5], FTYPE eup[][4], FTYPE elo[][4]);


//p2u.c
int p2u(FTYPE *p, FTYPE *u,FTYPE[][5],FTYPE[][4],FTYPE[][4]);
int pff2u(FTYPE *p, FTYPE *u,FTYPE[][5],FTYPE[][4],FTYPE[][4]);
#endif


#if(0)
//frames.c
int boost22_ff2zamo(FTYPE T1[][4],FTYPE T2[][4],FTYPE *pp,FTYPE gg[][5],FTYPE eup[][4]);
int boost22_zamo2ff(FTYPE T1[][4],FTYPE T2[][4],FTYPE *pp,FTYPE gg[][5],FTYPE eup[][4]);
int boost2_zamo2ff(FTYPE A1[4],FTYPE A2[4],FTYPE *pp,FTYPE gg[][5],FTYPE eup[][4]);
int boost2_ff2zamo(FTYPE A1[4],FTYPE A2[4],FTYPE *pp,FTYPE gg[][5],FTYPE eup[][4]);
int trans22_lab2zamo(FTYPE T1[][4],FTYPE T2[][4],FTYPE gg[][5],FTYPE eup[][4]);
int trans22_zamo2lab(FTYPE T1[][4],FTYPE T2[][4],FTYPE gg[][5],FTYPE elo[][4]);
int trans2_lab2zamo(FTYPE *u1,FTYPE *u2,FTYPE e[][4]);
int trans2_zamo2lab(FTYPE *u1,FTYPE *u2,FTYPE e[][4]);
int indices_2221(FTYPE T1[][4],FTYPE T2[][4],FTYPE gg[][5]);
int indices_21(FTYPE A1[4],FTYPE A2[4],FTYPE gg[][5]);
int
indices_12(FTYPE A1[4],FTYPE A2[4],FTYPE GG[][5]);
#endif

#if(0)
int print_tensor(FTYPE T[][4]);
int print_4vector(FTYPE v[4]);
int print_Nvector(FTYPE v[4],int);
int prad_zamo2ff(FTYPE *pp1, FTYPE *pp2, FTYPE gg[][5], FTYPE eup[][4]);
int prad_ff2zamo(FTYPE *pp1, FTYPE *pp, FTYPE gg[][5], FTYPE eup[][4]);
#endif

#if(0)
//rad.c
int calc_Rij(FTYPE *pp, FTYPE  Rij[][4]);
int solve_explicit_ff(int ix,int iy,int iz,FTYPE dt,FTYPE* deltas);
int solve_implicit_ff(int ix,int iy,int iz,FTYPE dt,FTYPE* deltas);
FTYPE calc_LTE_EfromT(FTYPE);
FTYPE calc_LTE_TfromE(FTYPE);
FTYPE calc_LTE_Efromurho(FTYPE E,FTYPE);
FTYPE calc_PEQ_ufromTrho(FTYPE,FTYPE);
FTYPE calc_PEQ_Tfromurho(FTYPE,FTYPE);
int calc_LTE_ff(FTYPE,FTYPE*,FTYPE*,FTYPE,int);
int solve_LTE_ff(int ix,int iy,int iz,FTYPE dt);
int solve_LTE(int ix,int iy,int iz,FTYPE dt);
int solve_radforce_ff(int ix,int iy,int iz,FTYPE dt);
int solve_radforce(int ix,int iy,int iz,FTYPE dt);
int calc_tautot(FTYPE *pp, FTYPE *xx, FTYPE *dl, FTYPE *tautot);
int calc_tauabs(FTYPE *pp, FTYPE *xx, FTYPE *dl, FTYPE *tauabs);
int calc_Gi(FTYPE *pp, FTYPE Gi[4]);
#endif
