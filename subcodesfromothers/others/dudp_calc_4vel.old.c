
/* 
 * calculate du/dp analytically.  
 * p is the full (NPR) element primitive variables
 * alpha is be a 5x5 matrix containing dU(1-5)/dp(1-5).
 * used only by nutoprim
 * 
 *
 * cfg 7-11-01
 *
 * cfg 8-23-01 heavy repairs
 *
 */

#include "decs.h"

struct of_geom *ptrgeom;

int dudp_calc_4vel(FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha)
{
	static FTYPE dutdui[NDIM] ;
	static FTYPE dbtdui[NDIM] ;
	static FTYPE dbsqdui[NDIM] ;
	static FTYPE dbiduj[NDIM][NDIM] ;
	static FTYPE tmp1[NDIM],tmp2[NDIM] ;
	FTYPE w,bsq ;
	int j,k ;

	ptrgeom=geom; // so don't have to rewrite functions yet

	for(j=1;j<=5;j++)
	for(k=1;k<=5;k++) alpha[j][k] = 0. ;


	bsq = dot(q->bcon, q->bcov);

	/*
	fprintf(fail_file,"rcurr,hcurr: %g %g\n",rcurr,hcurr) ;
	fprintf(fail_file,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
			pr[0], pr[1], pr[2], pr[3], pr[4], pr[5], pr[6], pr[7]) ; 
	*/

	alpha[RHO+1][RHO+1] = q->ucon[TT] ;
	alpha[UU+1][RHO+1] = q->ucon[TT]*q->ucon[TT] ;
	alpha[U1+1][RHO+1] = q->ucon[TT]*q->ucon[1] ;
	alpha[U2+1][RHO+1] = q->ucon[TT]*q->ucon[2] ;
	alpha[U3+1][RHO+1] = q->ucon[TT]*q->ucon[3] ;

	alpha[UU+1][UU+1] = gam*q->ucon[TT]*q->ucon[TT] + (gam - 1.)*geom->gcon[TT][TT] ;
	alpha[U1+1][UU+1] = gam*q->ucon[TT]*q->ucon[1] + (gam - 1.)*geom->gcon[TT][1] ;
	alpha[U2+1][UU+1] = gam*q->ucon[TT]*q->ucon[2] + (gam - 1.)*geom->gcon[TT][2] ;
	alpha[U3+1][UU+1] = gam*q->ucon[TT]*q->ucon[3] + (gam - 1.)*geom->gcon[TT][3] ;

	dutdui_calc(q->ucon,dutdui) ;
	dbtdui_calc(dutdui,pr,dbtdui) ;
	dbiduj_calc(dbtdui,dutdui,q->ucon,q->bcon,dbiduj) ;
	dbsqdui_calc(dbiduj,q->bcon,dbsqdui) ;

	w = pr[RHO] + gam*pr[UU] + bsq ;

	alpha[RHO+1][U1+1] = pr[0]*dutdui[1] ;
	alpha[UU+1][U1+1] = dbsqdui[1]*q->ucon[TT]*q->ucon[TT] 
		+ 2.*w*dutdui[1]*q->ucon[TT]
		+ 0.5*dbsqdui[1]*geom->gcon[TT][TT] 
		- 2.*dbtdui[1]*q->bcon[TT]  ;
	alpha[U1+1][U1+1] = dbsqdui[1]*q->ucon[TT]*q->ucon[1] 
		+ w*dutdui[1]*q->ucon[1]
		+ w*q->ucon[TT]
		+ 0.5*dbsqdui[1]*geom->gcon[TT][1] 
		- dbtdui[1]*q->bcon[1] 
		- q->bcon[TT]*dbiduj[1][1] ;
	alpha[U2+1][U1+1] = dbsqdui[1]*q->ucon[TT]*q->ucon[2] 
		+ w*dutdui[1]*q->ucon[2]
		+ 0.5*dbsqdui[1]*geom->gcon[TT][2] 
		- dbtdui[1]*q->bcon[2] 
		- q->bcon[TT]*dbiduj[2][1] ;
	alpha[U3+1][U1+1] = dbsqdui[1]*q->ucon[TT]*q->ucon[3] 
		+ w*dutdui[1]*q->ucon[3]
		+ 0.5*dbsqdui[1]*geom->gcon[TT][3] 
		- dbtdui[1]*q->bcon[3] 
		- q->bcon[TT]*dbiduj[3][1] ;

	alpha[RHO+1][U2+1] = pr[0]*dutdui[2] ;
	alpha[UU+1][U2+1] = dbsqdui[2]*q->ucon[TT]*q->ucon[TT] 
		+ 2.*w*dutdui[2]*q->ucon[TT]
		+ 0.5*dbsqdui[2]*geom->gcon[TT][TT] 
		- dbtdui[2]*q->bcon[TT] 
		- q->bcon[TT]*dbiduj[TT][2] ;
	alpha[U1+1][U2+1] = dbsqdui[2]*q->ucon[TT]*q->ucon[1] 
		+ w*dutdui[2]*q->ucon[1]
		+ 0.5*dbsqdui[2]*geom->gcon[TT][1] 
		- dbtdui[2]*q->bcon[1] 
		- q->bcon[TT]*dbiduj[1][2] ;
	alpha[U2+1][U2+1] = dbsqdui[2]*q->ucon[TT]*q->ucon[2] 
		+ w*dutdui[2]*q->ucon[2]
		+ 0.5*dbsqdui[2]*geom->gcon[TT][2] 
		+ w*q->ucon[TT]
		- dbtdui[2]*q->bcon[2] 
		- q->bcon[TT]*dbiduj[2][2] ;
	alpha[U3+1][U2+1] = dbsqdui[2]*q->ucon[TT]*q->ucon[3] 
		+ w*dutdui[2]*q->ucon[3]
		+ 0.5*dbsqdui[2]*geom->gcon[TT][3] 
		- dbtdui[2]*q->bcon[3] 
		- q->bcon[TT]*dbiduj[3][2] ;

	alpha[RHO+1][U3+1] = pr[0]*dutdui[3] ;
	alpha[UU+1][U3+1] = dbsqdui[3]*q->ucon[TT]*q->ucon[TT] 
		+ 2.*w*dutdui[3]*q->ucon[TT]
		+ 0.5*dbsqdui[3]*geom->gcon[TT][TT] 
		- dbtdui[3]*q->bcon[TT] 
		- q->bcon[TT]*dbiduj[TT][3] ;
	alpha[U1+1][U3+1] = dbsqdui[3]*q->ucon[TT]*q->ucon[1] 
		+ w*dutdui[3]*q->ucon[1]
		+ 0.5*dbsqdui[3]*geom->gcon[TT][1] 
		- dbtdui[3]*q->bcon[1] 
		- q->bcon[TT]*dbiduj[1][3] ;
	alpha[U2+1][U3+1] = dbsqdui[3]*q->ucon[TT]*q->ucon[2] 
		+ w*dutdui[3]*q->ucon[2]
		+ 0.5*dbsqdui[3]*geom->gcon[TT][2] 
		- dbtdui[3]*q->bcon[2] 
		- q->bcon[TT]*dbiduj[2][3] ;
	alpha[U3+1][U3+1] = dbsqdui[3]*q->ucon[TT]*q->ucon[3] 
		+ w*dutdui[3]*q->ucon[3]
		+ 0.5*dbsqdui[3]*geom->gcon[TT][3] 
		+ w*q->ucon[TT]
		- dbtdui[3]*q->bcon[3] 
		- q->bcon[TT]*dbiduj[3][3] ;

	/* mixed indices version of code, i.e. conserved variables
	   2-5 are T^t_j rather than T^{tj}.  Convert from original
	   results for unmixed indices */
	for(k=1;k<=5;k++) {
		for(j=2;j<=5;j++) tmp1[j-2] = alpha[j][k] ;
		lower(tmp1,geom,tmp2) ;
		for(j=2;j<=5;j++) alpha[j][k] = tmp2[j-2] ;
	}

	/* this bit of legacy code can be uncommented if
	   the rest-mass flux is subtracted out from the energy 
	   flux, which may be numerically convenient (tests are unconclusive) */
	j = 2 ; for(k=1;k<=5;k++) alpha[j][k] += alpha[1][k] ;
	/*
	*/

	/* N.B.: all the conserved variables contain a factor of \sqrt{det(g_{\mu\nu})} */
	for(j=1;j<=5;j++)
	for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;

	return(0) ;
}

void dutdui_calc(FTYPE *ucon,FTYPE *dutdui) 
{
	int j ;
	FTYPE ucov[NDIM] ;

	lower(ucon,ptrgeom,ucov) ;
	SLOOPA dutdui[j] = -ucov[j]/ucov[0] ;

	return ;
}

void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui) 
{
	int j ;
	FTYPE B[NDIM],Bcov[NDIM] ;

	B[0] = 0. ;
	SLOOPA B[j] = pr[j+B1-1] ;

	lower(B,ptrgeom,Bcov) ;
	SLOOPA dbtdui[j] = Bcov[j] + Bcov[0]*dutdui[j] ;

	return ;
}

void dbiduj_calc(FTYPE *dbtdui,FTYPE *dutdui,FTYPE *ucon, FTYPE *b, FTYPE dbiduj[][NDIM]) 
{
	int j,k ;

	DLOOP dbiduj[j][k] = 0. ;

	SLOOP dbiduj[j][k] = -b[j]*dutdui[k]/ucon[TT] 
			+ ucon[j]*dbtdui[k]/ucon[TT] ;

	SLOOPA dbiduj[j][j] += b[TT]/ucon[TT] ;

	SLOOPA dbiduj[TT][j] = dbtdui[j] ;

	return ;
}

void dbsqdui_calc(FTYPE dbiduj[][NDIM],FTYPE *b, FTYPE *dbsqdui) 
{
	int j,k ;
	FTYPE bcov[NDIM] ;

	lower(b,ptrgeom,bcov) ;
	DLOOPA dbsqdui[j] = 0. ;
	DLOOP dbsqdui[j] += 2.*bcov[k]*dbiduj[k][j] ;

	return ;
}

FTYPE contract(FTYPE *vcon1, FTYPE *vcon2)
{
	int j,k ;
	FTYPE n ;

	n = 0. ;
	DLOOP n += ptrgeom->gcov[j][k]*vcon1[j]*vcon2[k] ;
	return(n) ;

}

FTYPE covtract(FTYPE *vcov1, FTYPE *vcov2)
{
	int j,k ;
	FTYPE n ;

	n = 0. ;
	DLOOP n += ptrgeom->gcon[j][k]*vcov1[j]*vcov2[k] ;
	return(n) ;

}

void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]) 
{
	int i,j ;
	FTYPE ucov[NDIM] ;

	lower(ucon,ptrgeom,ucov) ;

	for(i=1;i<=3;i++)
	for(j=1;j<=3;j++) 
		duudud[i][j] = ptrgeom->gcon[i][j] + ptrgeom->gcon[i][TT]*(-ucon[j]/ucon[TT]) ;

}
