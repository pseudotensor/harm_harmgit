
/*
 *
 * invert U (conserved variables) to obtain
 * p (primitive variables).  
 * NB: pr must contain initial guess
 *
 *
 * new version using Del Zanna's scheme for reducing
 * complexity of solution.
 *
 * cfg, ldz 16 jan 03
 *
 */

#include "decs.h"

double Dc,Ec,Qsq,Sc,Bsq,normc,Tsq ;
double Wc ;

/* pr *MUST* contain initial guess */
void Utoprim_ldz(double *U, struct of_geom *geom, double *pr)
{

	double alpha, beta[NDIM], Tcov[NDIM], Tcon[NDIM] ;
	double B[NDIM],Q[NDIM],W,vsq,ut,gamma,v[NDIM] ;
	void wvsq_solv_ldz(double *vsq, double *W) ;
	struct of_state q ;
	int j,k ;


	alpha = 1./sqrt(-geom->gcon[0][0]) ;
	beta[0] = 0. ;
	SLOOPA beta[j] = geom->gcon[0][j]*alpha*alpha ;

	/* undo density addition, renormalized conserved variables */
	U[UU] -= U[RHO] ;
	PLOOP U[k] /= geom->g ;

	/* raise index on mhd stress tensor */
	DLOOPA Tcov[j] = U[j+UU] ;
	raise(Tcov,geom,Tcon) ;

	Dc = alpha*U[0] ;
	SLOOPA Q[j] = alpha*(Tcon[j] + beta[j]*Tcon[0]) ;
	Ec = alpha*alpha*Tcon[0] ;
	SLOOPA B[j] = alpha*U[B1+j-1] ;

	/* calculate Q^2 */
	Qsq = 0. ;
	SLOOP Qsq += Q[j]*Q[k]*geom->gcov[j][k] ;

	/* calc. S */
	Sc = 0. ;
	SLOOP Sc += geom->gcov[j][k]*Q[j]*B[k] ;

	/* calc. Bsq */
	Bsq = 0. ;
	SLOOP Bsq += geom->gcov[j][k]*B[j]*B[k] ;

	Tsq = Bsq*Qsq - Sc*Sc ;

	/* calc. guess for W */
	get_state(pr,geom,&q) ;
	vsq = 1. - 1./(alpha*alpha*q.ucon[0]*q.ucon[0]) ;
	normc = U[RHO] ;

	/* now solve */
	wvsq_solv_ldz(&vsq,&W) ;

	/* now invert to find primitives */
	gamma = 1./sqrt(1. - vsq) ;
	pr[RHO] = Dc/gamma ;
	pr[UU] = ((1. - vsq)*W - pr[RHO])/gam ;

	SLOOPA v[j] = (Q[j] + (Sc/W)*B[j])/(W + Bsq) ;

	pr[U1] = alpha*v[1] - beta[1] ;
	pr[U2] = alpha*v[2] - beta[2] ;
	pr[U3] = alpha*v[3] - beta[3] ;

	pr[B1] = U[B1] ;
	pr[B2] = U[B2] ;
	pr[B3] = U[B3] ;

}

void wvsq_solv_ldz(double *vsq, double *W) 
{
	double tol,x1,x2 ;
	double rtsafe(void (*funcd)(), double x1, double x2, double xacc) ;
	extern void func(double x,double *f,double *df) ;
	double nrunsafe(void (*funcd)(), double guess) ;

	x1 = 0. - SMALL ;
	x2 = 1. - 1.e-6 ;

	//x1 = 0.5*(*vsq) ;
	//x2 = 0.5 + 0.5*(*vsq) ;

	*vsq = nrunsafe(func,*vsq) ;

	//tol = 1.e-4 ;
	//*vsq = rtsafe(func,x1,x2,tol) ;
	
	*W = Wc ;

	return ;

}

#define TOL 		1.e-8
#define NITERMAX  	10
#define NITERMIN  	0

double nrunsafe(void (*funcd)(), double guess)
{
	double f,df ;
	int n_iter ;

	(*funcd)(guess,&f,&df);

	n_iter = 0 ;
	while((fabs(f) > TOL && n_iter < NITERMAX)||(n_iter<NITERMIN)) {
		guess -= f/df ;
		(*funcd)(guess,&f,&df);
		n_iter++ ;
		nstroke++;
	}

	if(n_iter == NITERMAX) {
		fprintf(stderr,"max iterations exceeded\n") ;
		exit(1) ;
	}
	else
		return(guess) ;

}

void func(double x,double *f,double *df)
{
	double vsq,glf1,cc,igam1,ee,wsq,dw ;
	double a2,a1,a0,q,r,th,w1,vp2,cosa ;

	vsq = x ;
	glf1 = sqrt(1. - vsq) ;
	igam1 = (gam - 1.)/gam ;
	cc = 1./(1. - igam1*(1. - vsq)) ;
	ee = (Ec - igam1*glf1*Dc - 0.5*Bsq)*cc ;

	/*
	if(ee < 0.) {
		fprintf(stderr,"ee < 0 in func\n") ;
		exit(40) ;
	}
	*/

	if(Tsq > 0.) {
		a2 = (2.*Bsq - ee)/3. ;
		a1 = (Bsq - 2.*ee)*Bsq ;
		a0 = 0.5*Tsq*cc - ee*Bsq*Bsq ;
		q = a1/3. - a2*a2 ;
		r = 0.5*(a1*a2 - a0) - a2*a2*a2 ;
		cosa = r/sqrt(-q*q*q) ;
		if(cosa > 1.) cosa = 1. ;
		th = acos(cosa) ;
		
		//fprintf(stderr,">>> %g %g %g\n",q,r/sqrt(-q*q*q),th) ;

		Wc = 2.*sqrt(-q)*cos(th/3.) - a2 ;
		wsq = Wc*Wc ;
		w1 = 1./(Wc + Bsq) ;
		vp2 = Tsq*w1*w1 ;
		dw = -igam1*cc*(Wc - 0.5*Dc/glf1)/(1. + 2.*w1*(Wc - ee)) ;
		*f =  ( wsq*vsq+vp2*(2.*Wc+Bsq)-Qsq    )/(normc*normc) ;
		*df = ( wsq + 2.*(Wc*dw*(vsq - vp2*w1)))/(normc*normc) ;
	}
	else {
		Wc = ee ;
		wsq = Wc*Wc ;
		dw = -igam1*cc*(Wc - 0.5*Dc/glf1) ;
		*f = wsq*vsq - Qsq ;
		*df = wsq + 2.*Wc*dw*vsq ;
	}
}

#define MAXIT 100

double rtsafe(void (*funcd)(), double x1, double x2, double xacc)
{
	void nrerror();
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	(*funcd)(x1,&fl,&df);
	(*funcd)(x2,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		for(j=0;j<100;j++) {
			(*funcd)(0.01*j,&fl,&df) ;
			fprintf(stderr,"%d %g %g\n",j,fl,df);
		}
		/*
		*/
		nrerror("Root must be bracketed in rtsafe");
	}
	if (fabs(fl) < SMALL) return x1;
	if (fabs(fh) < SMALL) return x2;
	
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,&f,&df);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,&f,&df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT

