
extern void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *radstressdir);
extern ldouble my_min(ldouble a, ldouble b);
extern ldouble my_sign(ldouble x);
extern int inverse_44matrix(ldouble a[][4], ldouble ia[][4]);
extern int boost22_fforzamo(int whichdir, ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);
extern int boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);
extern int boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);
extern int trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble elo[][4]);
extern int trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble eup[][4]);
extern int trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble eup[][4]);
extern int trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble elo[][4]);

extern int boost2_zamo2ff(ldouble A1[],ldouble A2[],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);
extern int boost2_ff2zamo(ldouble A1[],ldouble A2[],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom,ldouble eup[][4]);
extern int boost2_fforzamo(int whichdir, ldouble A1[4],ldouble A2[4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom,ldouble eup[][4]);

extern int calc_Rij(ldouble *pp, ldouble Rij[][4]);

extern int indices_2221(ldouble T1[][NDIM],ldouble T2[][NDIM], struct of_geom *ptrgeom);
extern int indices_21(ldouble A1[NDIM],ldouble A2[NDIM],struct of_geom *ptrgeom);
extern int indices_12(ldouble A1[NDIM],ldouble A2[NDIM],struct of_geom *ptrgeom);

extern int calc_LNRFes(struct of_geom *ptrgeom, ldouble emuup[][4], ldouble emulo[][4]);


extern int p2u_rad(ldouble *pr, ldouble *Urad, struct of_geom *ptrgeom, struct of_state *q);


extern int prad_ff2zamo(ldouble *pp1, ldouble *pp2, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);
extern int f_prad_zamo2ff(ldouble *ppff, ldouble *ppzamo, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4],ldouble *f);
extern int prad_zamo2ff(ldouble *ppzamo, ldouble *ppff, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4]);

extern int u2p_rad(ldouble *uu, ldouble *pp, struct of_geom *ptrgeom);


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
ldouble calc_totalmass();
int initialize_arrays();
int free_arrays();
ldouble my_min(ldouble a, ldouble b);
//ldouble my_max(ldouble a, ldouble b);
ldouble my_sign(ldouble);
int find_eigenvalues3(ldouble[3][3],ldouble*);
int find_eigenvalues(ldouble*,int);
ldouble find_max_eigenvalue(ldouble *data,int N);
int find_max_eigenvalue_lr(ldouble *data,int N,ldouble*,ldouble*);
int my_err(char *);
int inverse_44matrix(ldouble a[][4], ldouble ia[][4]);
int convert_out2gif_1d(char *fname,char*,int niter,ldouble t);
int convert_out2gif_2d(char *fname,char*,int niter,ldouble t);
int getch(void);
int dosthelse(void);
//finite.c
ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz);

int f_timeder (ldouble t, ldouble dt,ldouble, ldouble*, int ifcopy2u0, ldouble*);
int set_grid(ldouble*, ldouble*, ldouble*,ldouble*);
int print_grid(ldouble,ldouble,ldouble);
ldouble fd_flux_limiter(ldouble r);
ldouble minmod_fd_flux_limiter(ldouble ,ldouble,ldouble);
ldouble f_der_kurganovtadmor(int ix,int iy, int yz,ldouble*);
ldouble f_der_hlle_obsolete(int ix,int iy, int yz,ldouble*);
ldouble f_der_muscl(int ix,int iy, int yz,ldouble*);
//ldouble get_x(int,int);
#endif

#if(0)
#define get_x(ic,idim) (idim==0 ? x[ic+NG] : (idim==1 ? x[ic+NG + NX+2*NG] : (idim==2 ? x[ic+NG + NX+2*NG + NY+2*NG ] : 0.)))
//ldouble get_xb(int,int);
#define get_xb(ic,idim) (idim==0 ? xb[ic+NG] : (idim==1 ? xb[ic+NG + NX+2*NG + 1] : (idim==2 ? xb[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)))
int set_x(int,int,ldouble);
int set_xb(int,int,ldouble);
//ldouble get_u(ldouble*,int,int,int,int);
#define get_cflag(iflag,ix,iy,iz) (cellflag[iflag + (ix+NG)*NFLAGS + (iy+NG)*(NX+2*NG)*NFLAGS + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NFLAGS])
#define set_cflag(iflag,ix,iy,iz,val) cellflag[iflag + (ix+NG)*NFLAGS + (iy+NG)*(NX+2*NG)*NFLAGS + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NFLAGS]=val
#define get_u(uarr,iv,ix,iy,iz) (uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV])
//int set_u(ldouble*,int,int,int,int,ldouble);
#define set_u(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val
ldouble get_u_scalar(ldouble*,int,int,int);
//int set_u_scalar(ldouble*,int,int,int,ldouble);
#define set_u_scalar(uarr,ix,iy,iz,val) uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)] = val
int copy_u(ldouble,ldouble*,ldouble*);
int add_u(ldouble f1, ldouble* u1, ldouble f2, ldouble *u2, ldouble *u3);
#endif

#if(0)
ldouble f_timeder_source_term(ldouble t, const ldouble y[], ldouble f[],  void *params);
int if_indomain(int ix,int iy,int iz);
int if_outsidegc(int ix,int iy,int iz);
int if_outsidewave(int ix,int iy,int iz);
#endif


#if(0)
ldouble get_size_x(int ic, int idim);
#endif

#if(0)
int set_ub(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value,int idim);
//#define set_ub(uarr,iv,ix,iy,iz,idim,val) (idim==0 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV]=val : (idim==1 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV]=val : (idim==2 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val : 0.)))
#define set_ubx(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV]=val
#define set_uby(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV]=val
#define set_ubz(uarr,iv,ix,iy,iz,val) uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV]=val
//ldouble get_ub(ldouble* uarr,int iv,int ix,int iy,int iz,int idim);
#define get_ub(uarr,iv,ix,iy,iz,idim) (idim==0 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV] : (idim==1 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV] : (idim==2 ? uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV] : 0.)))
//ldouble get_g(ldouble* uarr, int i,int j, int ix, int iy, int iz);
#define get_g(uarr,i,j,ix,iy,iz) uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE]
int set_g(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);
int set_T(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);

#define get_T(uarr,i,j,ix,iy,iz) uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*16]
#define get_Tb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG+1)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*16] : (idim==1 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*16] : (idim==2 ? uarr[i*4+j + (ix+NG)*16 + (iy+NG)*(NX+2*NG)*16 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*16] : 0.)))
int set_Tb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);
int set_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);
//ldouble get_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,int idim);
#define get_gb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG+1)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*gSIZE] : (idim==1 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*gSIZE] : (idim==2 ? uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE] : 0.)))
#define get_gKr(i,j,k,ix,iy,iz) gKr[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*64]
#define get_gKrb(i,j,k,ix,iy,iz,idim) (idim==0 ? gKrbx[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG+1)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*64] : (idim==1 ? gKrby[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG+1)*(NX+2*NG+1)*64] : (idim==2 ? gKrbz[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*64] : 0.)))
#define set_gKr(i,j,k,ix,iy,iz,val) gKr[i*4*4+j*4+k + (ix+NG)*64 + (iy+NG)*(NX+2*NG)*64 + (iz+NG)*(NY+2*NG)*(NX+2*NG)*64]=val
int set_Krb(int i,int j,int k,int ix,int iy,int iz,ldouble value,int idim);
#endif


#if(0)
//fileop.c
int fread_restartfile(ldouble*);
int fprint_openfiles();
int fprint_closefiles();
int fprint_profiles(ldouble,ldouble);
int print_profiles();
#endif

#if(0)
ldouble max_eigen_Jac(ldouble *,ldouble*,int,void*);
int calc_wavespeeds(int,int,int,ldouble*,ldouble*,ldouble*,ldouble*,ldouble*,ldouble*);
int calc_wavespeeds_lr(int,int,int,ldouble*);
int calc_wavespeeds_lr_faces( int,int,int,int, ldouble*,ldouble*);
int max_eigen_lr_Jac(ldouble *,ldouble*,int,void*,ldouble*,ldouble*);
int calc_Jac_num(ldouble *xx,ldouble *ujac, int idim,void *parameters,ldouble *fd_jac);
int set_initial_profile();
ldouble f_flux_prime(ldouble *uu, int,int,int,int,ldouble *ff);
ldouble f_diffusion_prime(ldouble *uu, ldouble *du, int iv,int,void*);
ldouble f_grav_potential(ldouble,ldouble,ldouble);
ldouble f_der_grav_potential(ldouble,ldouble,ldouble,int);
int f_source_term(int,int,int,ldouble *);
int f_fourforce_source_term(int,int,int,ldouble *);
int f_implicit_rhou(int ix, int iy, int iz, ldouble *rho, ldouble *uint, ldouble dt);
int f_metric_source_term(int ix, int iy, int iz,ldouble *ss);
int f_metric_source_term_face(int ix, int iy, int iz,int idim,int ifleft,ldouble *ss);
int f_source_term_face(int ix, int iy, int iz,int idim,int ifleft,ldouble *ss);
int initialize_problem();
int calc_maxgravspeed(int ix, int iy, int iz, ldouble dt);
int calc_maxwavespeed_obsolete(int ix, int iy, int iz,void*);
int calc_maxwavespeed_osbolete_cs(int ix, int iy, int iz,void*);
ldouble calc_chi(ldouble rho, ldouble T);
int rad_fld_factors(int ix,int iy, int iz,ldouble *,void *rad_param,int);
int rad_fld_factors_arb(int ix,int iy, int iz,ldouble*,void *rad_param);

int radx_flux_implicit(ldouble *uu,ldouble dtt);
int LTE_implicit(ldouble *uu,ldouble dtt);
ldouble max_eigen_radx(ldouble nx,ldouble ny,ldouble nz,int idim);
int max_eigen_lr_radx(ldouble nx,ldouble ny,ldouble nz,int idim,ldouble *al, ldouble *ar);
#endif

int gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3);


#if(0)
ldouble f_der_hlle         (int ix,int iy,int iz, ldouble *fd_der);
ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z);
ldouble calc_kappaes(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z);
ldouble calc_ufromS(ldouble S,ldouble rho);
ldouble calc_Sfromu(ldouble S,ldouble u);
int
avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,ldouble*,ldouble*,ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2);

//problem.c
ldouble calc_xb(int i,int idim);
int calc_bc(int,int,int,ldouble,ldouble*, ldouble*);
int pr_tophat_inside(ldouble x,ldouble y,ldouble z);
int my_finger(ldouble);

//metric.c
int calc_LNRFes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4]);
int calc_g(ldouble*,ldouble[][5]);
int calc_G(ldouble*,ldouble[][5]);
int calc_Krzysie(ldouble*,ldouble[][4][4]);
int print_Krzysie(ldouble g[][4][4]);
int print_g(ldouble [][5]);
ldouble calc_gdet(ldouble *xx);
ldouble calc_dlgdet(ldouble *xx, int idim);
#endif


#if(0)
int pick_g(int ix,int iy,int iz,ldouble gg[][5]);
int pick_G(int ix,int iy,int iz,ldouble gg[][5]);
int pick_gb(int ix,int iy,int iz,int,ldouble gg[][5]);
int pick_Gb(int ix,int iy,int iz,int,ldouble gg[][5]);
int pick_T(ldouble *arr,int ix,int iy,int iz,ldouble T[][4]);
int pick_Tb(ldouble *arr,int ix,int iy,int iz,int,ldouble T[][4]);

int p2u_Sonly(ldouble *p, ldouble *u,ldouble[][5]);
int print_p(ldouble *p);
int print_u(ldouble *p);
int calc_Tmunu( ldouble *p, ldouble g[][5], ldouble T[][4],ldouble*);
int convert_uold2urel(ldouble *x,ldouble *u);
int calc_metric();
int calc_sourceterms(int,int,int);
int calc_primitives(int,int,int);
int calc_conserved(int ix,int iy,int iz);

ldouble r_horizon_BL(ldouble a);
ldouble r_mbound_BL(ldouble a);
ldouble r_photon_BL(ldouble a);
int update_entropy(int ix,int iy,int iz,int u2pflag);

//u2p.c
int u2p(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4], ldouble elo[][4]);
int u2p_hot(ldouble*,ldouble*,ldouble[][5]);
int u2p_entropy(ldouble*,ldouble*,ldouble[][5]);
int u2p_cold(ldouble*,ldouble*,ldouble[][5]);
int
u2p_rad_num(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble eup[][4], ldouble elo[][4]);
int
u2p_rad(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4], ldouble elo[][4]);

int
dump_u2p_rad(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble eup[][4], ldouble elo[][4]);


//p2u.c
int p2u(ldouble *p, ldouble *u,ldouble[][5],ldouble[][4],ldouble[][4]);
int pff2u(ldouble *p, ldouble *u,ldouble[][5],ldouble[][4],ldouble[][4]);
#endif


#if(0)
//frames.c
int boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble eup[][4]);
int boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble eup[][4]);
int boost2_zamo2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble eup[][4]);
int boost2_ff2zamo(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble eup[][4]);
int trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5],ldouble eup[][4]);
int trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5],ldouble elo[][4]);
int trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble e[][4]);
int trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble e[][4]);
int indices_2221(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
int indices_21(ldouble A1[4],ldouble A2[4],ldouble gg[][5]);
int
indices_12(ldouble A1[4],ldouble A2[4],ldouble GG[][5]);
#endif

#if(0)
int print_tensor(ldouble T[][4]);
int print_4vector(ldouble v[4]);
int print_Nvector(ldouble v[4],int);
int prad_zamo2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble eup[][4]);
int prad_ff2zamo(ldouble *pp1, ldouble *pp, ldouble gg[][5], ldouble eup[][4]);
#endif

#if(0)
//rad.c
int calc_Rij(ldouble *pp, ldouble  Rij[][4]);
int solve_explicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas);
int solve_implicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas);
ldouble calc_LTE_EfromT(ldouble);
ldouble calc_LTE_TfromE(ldouble);
ldouble calc_LTE_Efromurho(ldouble E,ldouble);
ldouble calc_PEQ_ufromTrho(ldouble,ldouble);
ldouble calc_PEQ_Tfromurho(ldouble,ldouble);
int calc_LTE_ff(ldouble,ldouble*,ldouble*,ldouble,int);
int solve_LTE_ff(int ix,int iy,int iz,ldouble dt);
int solve_LTE(int ix,int iy,int iz,ldouble dt);
int solve_radforce_ff(int ix,int iy,int iz,ldouble dt);
int solve_radforce(int ix,int iy,int iz,ldouble dt);
int calc_tautot(ldouble *pp, ldouble *xx, ldouble *dl, ldouble *tautot);
int calc_tauabs(ldouble *pp, ldouble *xx, ldouble *dl, ldouble *tauabs);
int calc_Gi(ldouble *pp, ldouble Gi[4]);
#endif
