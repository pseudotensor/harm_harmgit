#define GGG (6.674e-8)
#define CCC 2.998e10
#define MSUNCM 147700.

#define K_BOLTZ (1.3806488e-16L * GGG / CCC / CCC / CCC / CCC  )
#define M_PROTON (1.67262158e-24L * GGG / CCC / CCC)
#define SIGMA_RAD (5.67e-5 * GGG / CCC / CCC / CCC / CCC / CCC * MASSCM * MASSCM  * MASSCM)
//efine A_RAD (4.*5.67e-5/CCC * GGG / CCC / CCC / CCC / CCC)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     

#define SOURCETERMS
#define NUM_SOURCESTEPS 1

#define my_max(x,y) (x>y?x:y)


#define lenCGS2GU(x)    (x/MASSCM )
#define lenGU2CGS(x)    (x*MASSCM )
#define timeCGS2GU(x)   (x/MASSCM*CCC)
#define timeGU2CGS(x)   (x*MASSCM/CCC)
#define velCGS2GU(x)    (x/CCC)
#define velGU2CGS(x)    (x*CCC)
#define rhoCGS2GU(x)    (x*GGG/CCC/CCC*MASSCM*MASSCM*MASSCM)
#define rhoGU2CGS(x)    (x/GGG*CCC*CCC/MASSCM/MASSCM/MASSCM)
#define massCGS2GU(x)    (x*GGG/CCC/CCC)
#define massGU2CGS(x)    (x/GGG*CCC*CCC)
#define kappaGU2CGS(x)  (x*GGG/CCC/CCC*MASSCM*MASSCM)
#define kappaCGS2GU(x)  (x/GGG*CCC*CCC/MASSCM/MASSCM)
#define endenCGS2GU(x) (x*GGG*MASSCM*MASSCM*MASSCM/CCC/CCC/CCC/CCC)
#define endenGU2CGS(x) (x/GGG/MASSCM/MASSCM/MASSCM*CCC*CCC*CCC*CCC)
#define fluxCGS2GU(x) (x*GGG*MASSCM*MASSCM*MASSCM/CCC/CCC/CCC/CCC/CCC)
#define fluxGU2CGS(x) (x/GGG/MASSCM/MASSCM/MASSCM*CCC*CCC*CCC*CCC*CCC)

#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))
#define KAPPA_FF_COEFF (1.7e-25/1.67262158e-24/1.67262158e-24*CCC*CCC*CCC*CCC/GGG/GGG/MASSCM/MASSCM/MASSCM/MASSCM/MASSCM)
#define KAPPA_BF_COEFF (4.8e-24/1.67262158e-24/1.67262158e-24*CCC*CCC*CCC*CCC/GGG/GGG/MASSCM/MASSCM/MASSCM/MASSCM/MASSCM)



#include "problem.h"

#include "mdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>  
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multiroots.h>

#define ldouble long double
#define gSIZE 20 //size of metric arrays = 16 + 1 (gdet) + 3 (dlgdet)

ldouble *u,*x,*xb,*du,*ut1,*ut2,*ut3,*ut4,*ut0,*u_bak,*u_step1,*u_step2,*u_step3,*u_step4,*ahdx,*ahdy,*ahdz,*aradx,*arady,*aradz,
  *ahdxl,*ahdyl,*ahdzl,*aradxl,*aradyl,*aradzl,  *ahdxr,*ahdyr,*ahdzr,*aradxr,*aradyr,*aradzr,*p,*pt0,*px,*py,*pz,*s,*g,*gbx,*gby,*gbz,*Gbx,*Gby,*Gbz,
  *pbLx,*pbRx,*pbLy,*pbRy,*pbLz,*pbRz,*sbLx,*sbRx,*sbLy,*sbRy,*sbLz,*sbRz,*ubLx,*ubRx,*ubLy,*ubRy,*ubLz,*ubRz,
  *flbx,*flby,*flbz,*flLx,*flRx,*flLy,*flRy,*flLz,*flRz,*gKr,*gKrbx,*gKrby,*gKrbz,*G,*emuup,*emulo,*emuupbx,*emulobx,*emuupby,*emuloby,*emuupbz,*emulobz;
int *cellflag,**loop_1,**loop_2,Nloop_1,Nloop_2;
//ldouble ****u;
ldouble Kr_tmp[4][4][4],g_tmp[4][4];

ldouble inputarg[10];
int **gcidx;


ldouble max_ws[3],max_dt;
ldouble min_dx,min_dy,min_dz;
  

FILE *fout1,*fout_totmass;
int nfout1;

//main.c
int solve_all_problems_5(ldouble);
int solve_all_problems_6(ldouble);
gsl_odeiv2_step **odeiv2_step_1;
gsl_odeiv2_step **odeiv2_step_2;
struct evolve_fluxes_1_param
{
  ldouble ix,iy,iz,t,dt;
};
struct evolve_fluxes_2_param
{
  ldouble ix,iy,iz,t,dt;
  struct rad_parameters *rp;
};

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
ldouble f_timeder_source_term(ldouble t, const ldouble y[], ldouble f[],  void *params);
int check_floors(ldouble *uu);
int set_bc(ldouble);
int if_indomain(int ix,int iy,int iz);
int if_outsidegc(int ix,int iy,int iz);
int if_outsidewave(int ix,int iy,int iz);
ldouble get_size_x(int ic, int idim);
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


//fileop.c
int fread_restartfile(ldouble*);
int fprint_openfiles();
int fprint_closefiles();
int fprint_profiles(ldouble,ldouble);
int print_profiles();

//physics.c
struct rad_parameters
{
  ldouble f_Edd[3][3];
  ldouble F0[3];
  ldouble lambda;
  ldouble chi; 
  ldouble kappa;
  ldouble B;
  ldouble tau[3];
  ldouble tautot[3];
  ldouble T;
  ldouble dE0dx;
  ldouble R;
  ldouble f;
  ldouble x,y,z;
};


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
int gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3);
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


//relele.c
#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define kron(i,j) (i == j ? 1. : 0.)
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
int print_tensor(ldouble T[][4]);
int print_4vector(ldouble v[4]);
int print_Nvector(ldouble v[4],int);
int prad_zamo2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble eup[][4]);
int prad_ff2zamo(ldouble *pp1, ldouble *pp, ldouble gg[][5], ldouble eup[][4]);

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

