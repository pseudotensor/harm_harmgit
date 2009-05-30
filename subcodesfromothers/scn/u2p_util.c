
#include "u2p_defs.h"


void primtoU_g( double prim[], double gcov[][4], double gcon[][4],  double U[] );
void ucon_calc_g(double prim[],double gcov[][4],double gcon[][4],double ucon[]);
void raise_g(double vcov[], double gcon[][4], double vcon[]);
void lower_g(double vcon[], double gcov[][4], double vcov[]);
void ncov_calc(double gcon[][4],double ncov[]) ;
void bcon_calc_g(double prim[],double ucon[],double ucov[],double ncov[],double bcon[]); 
double pressure_rho0_u(double rho0, double u);
double pressure_rho0_w(double rho0, double w);
void dualfprintf(FILE* fileptr, char *format, ...);

/* 
 * reference implementation of transformation from
 * primitive to conserved variables.
 *
 * cfg 7-6-04
 *
 * input: 
 * 
 * primitive variables in 8 element array 
 * metric in contravariant and covariant form
 *
 * output:
 * 
 * conserved variables in 8 element array 
 * 
 */

void primtoU_g(
	double prim[8],       /* primitive variables */
	double gcov[4][4],    /* covariant (index dn) form of metric */
	double gcon[4][4],    /* contravariant (index up) form of metric */
	double U[8]           /* matrix of derivatives */
) {
	int i,j ;
	double rho0 ;
	double ucon[4],ucov[4],bcon[4],bcov[4],ncov[4] ;
	double gamma,n_dot_b,bsq,u,p,w ;

	
	/* preliminaries */
	ucon_calc_g(prim,gcov,gcon,ucon) ;
	lower_g(ucon,gcov,ucov) ;
	ncov_calc(gcon,ncov) ;

	gamma = -ncov[0]*ucon[0] ;

	bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
	lower_g(bcon,gcov,bcov) ;

   	n_dot_b = 0. ;
	for(i=0;i<4;i++) n_dot_b += ncov[i]*bcon[i] ;
	bsq = 0. ;
	for(i=0;i<4;i++) bsq += bcov[i]*bcon[i] ;

	rho0 = prim[RHO] ;
	u = prim[UU] ;
	p = pressure_rho0_u(rho0,u) ;
	w = rho0 + u + p ;

	U[RHO] = gamma*rho0 ;

	for(i=0;i<4;i++) 
		U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
			- (p + bsq/2.)*ncov[i] 
			+ n_dot_b*bcov[i] ;

	U[BCON1] = prim[BCON1] ;
	U[BCON2] = prim[BCON2] ;
	U[BCON3] = prim[BCON3] ;

	return ;
}

/* find the contravariant fluid four-velocity from primitive 
   variables plus the metric */
void ucon_calc_g(double prim[8],double gcov[4][4],double gcon[4][4],double ucon[4])
{
	double u_tilde_con[4] ;
	double u_tilde_sq ;
	double gamma,lapse ;
	int i,j ;
	
	u_tilde_con[0] = 0. ;
	u_tilde_con[1] = prim[UTCON1] ;
	u_tilde_con[2] = prim[UTCON2] ;
	u_tilde_con[3] = prim[UTCON3] ;

	u_tilde_sq = 0. ;
	for(i=0;i<4;i++)
	for(j=0;j<4;j++)
		u_tilde_sq += gcov[i][j]*u_tilde_con[i]*u_tilde_con[j] ;
	u_tilde_sq = fabs(u_tilde_sq) ;

	gamma = sqrt(1. + u_tilde_sq) ;

	lapse = sqrt(-1./gcon[0][0]) ;

	for(i=0;i<4;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[0][i] ;

	return ;
}

/* raise covariant vector vcov using gcon, place result in vcon */
void raise_g(double vcov[4], double gcon[4][4], double vcon[4])
{
	int i,j;

	for(i=0;i<4;i++) {
		vcon[i] = 0. ;
		for(j=0;j<4;j++) 
			vcon[i] += gcon[i][j]*vcov[j] ;
	}

	return ;
}
/* lower contravariant vector vcon using gcov, place result in vcov */
void lower_g(double vcon[4], double gcov[4][4], double vcov[4])
{
	int i,j;

	for(i=0;i<4;i++) {
		vcov[i] = 0. ;
		for(j=0;j<4;j++) 
			vcov[i] += gcov[i][j]*vcon[j] ;
	}

	return ;
}

/* set covariant normal observer four-velocity */
void ncov_calc(double gcon[4][4],double ncov[4]) 
{
	double lapse ;

	lapse = sqrt(-1./gcon[0][0]) ;

	ncov[0] = -lapse ;
	ncov[1] = 0. ;
	ncov[2] = 0. ;
	ncov[3] = 0. ;

	return ;
}

/* calculate contravariant magnetic field four-vector b */
void bcon_calc_g(double prim[8],double ucon[4],double ucov[4],double ncov[4],double bcon[4]) 
{
	double Bcon[4] ;
	double u_dot_B ;
	double gamma ;
	int i ;

	Bcon[0] = 0. ;
	for(i=1;i<4;i++) Bcon[i] = prim[BCON1+i-1] ;

	u_dot_B = 0. ;
	for(i=0;i<4;i++) u_dot_B += ucov[i]*Bcon[i] ;

	gamma = -ucon[0]*ncov[0] ;
	for(i=0;i<4;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/* 

pressure as a function of rho0 and u 

this is used by primtoU and Utoprim_?D

*/
double pressure_rho0_u(double rho0, double u)
{
	return((GAMMA - 1.)*u) ;
}


  
/* 

pressure as a function of rho0 and w = rho0 + u + p 

this is used by primtoU and Utoprim_1D

*/
double pressure_rho0_w(double rho0, double w)
{
	return((GAMMA-1.)*(w - rho0)/GAMMA) ;
}


void dualfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;
 
  va_start (arglist, format);
 
  vfprintf (stderr, format, arglist);
  fflush(stderr);

  va_end (arglist);
}
