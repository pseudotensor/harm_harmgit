
#include "decs.h"

void flux_ct(double (* p_h)[N2+2*NG][NP]) 
{
	int i,j ;
	static double emf[N1+1][N2+1] ;
	double v1_av,w1,v2_av,w2,w_12 ;
	double rho_av,sqrtrho,b1_av,b2_av,va1,va2 ;
	double emfmm,emfpm,emfmp,emfpp,alpha ;
	double B1pp,B1pm,B1mp,B1mm,B2pp,B2pm,B2mp,B2mm ;
	double U1pp,U1pm,U1mp,U1mm,U2pp,U2pm,U2mp,U2mm ;
	double cms_func(double *prim_var) ;
	double B1d,B1u,B1l,B1r,B2d,B2u,B2l,B2r ;

#if 0
	/* Toth approach: just average */
	for(i=0;i<N1+1;i++)	
	for(j=0;j<N2+1;j++) {
		emf[i][j] = 0.25*(F1[i][j][B2] + F1[i][j-1][B2]
				- F2[i][j][B1] - F2[i-1][j][B1]) ;
	}
#endif
#if 0
	/* Stone & Gardiner eq. 39 */
	for(i=0;i<N1+1;i++)	
	for(j=0;j<N2+1;j++) {
		emfmp = p_h[i-1][j  ][B2]*p_h[i-1][j  ][U1] -
			p_h[i-1][j  ][B1]*p_h[i-1][j  ][U2] ;
		emfmm = p_h[i-1][j-1][B2]*p_h[i-1][j-1][U1] -
			p_h[i-1][j-1][B1]*p_h[i-1][j-1][U2] ;
		emfpm = p_h[i  ][j-1][B2]*p_h[i  ][j-1][U1] -
			p_h[i  ][j-1][B1]*p_h[i  ][j-1][U2] ;
		emfpp = p_h[i  ][j  ][B2]*p_h[i  ][j  ][U1] -
			p_h[i  ][j  ][B1]*p_h[i  ][j  ][U2] ;

		emf[i][j] = 0.5*(F1[i][j][B2] + F1[i][j-1][B2]
				- F2[i][j][B1] - F2[i-1][j][B1]) 
			    - 0.25*(emfmp + emfmm + emfpm + emfpp) ;
	}
#endif
#if 1
	/* Stone & Gardiner eq. 48 */
	for(i=0;i<N1+1;i++)	
	for(j=0;j<N2+1;j++) {
		emfmp = p_h[i-1][j  ][B2]*p_h[i-1][j  ][U1] -
			p_h[i-1][j  ][B1]*p_h[i-1][j  ][U2] ;
		emfmm = p_h[i-1][j-1][B2]*p_h[i-1][j-1][U1] -
			p_h[i-1][j-1][B1]*p_h[i-1][j-1][U2] ;
		emfpm = p_h[i  ][j-1][B2]*p_h[i  ][j-1][U1] -
			p_h[i  ][j-1][B1]*p_h[i  ][j-1][U2] ;
		emfpp = p_h[i  ][j  ][B2]*p_h[i  ][j  ][U1] -
			p_h[i  ][j  ][B1]*p_h[i  ][j  ][U2] ;



		B1d = 0.5*(
			p_h[i-1][j-1][B1] + 0.5*dq1[i-1][j-1][B1] +
			p_h[i][j-1][B1] - 0.5*dq1[i][j-1][B1]) ;
		B1u = 0.5*(
			p_h[i-1][j][B1] + 0.5*dq1[i-1][j][B1] +
			p_h[i][j][B1] - 0.5*dq1[i][j][B1]) ;

		B2l = 0.5*(
			p_h[i-1][j-1][B2] + 0.5*dq2[i-1][j-1][B2] +
			p_h[i-1][j][B2] - 0.5*dq2[i-1][j][B2]) ;
		B2r = 0.5*(
			p_h[i][j-1][B2] + 0.5*dq2[i][j-1][B2] +
			p_h[i][j][B2] - 0.5*dq2[i][j][B2]) ;

		B1mm = p_h[i-1][j-1][B1] ;
		B1mp = p_h[i-1][j][B1] ;
		B1pm = p_h[i][j-1][B1] ;
		B1pp = p_h[i][j][B1] ;
		B2mm = p_h[i-1][j-1][B2] ;
		B2mp = p_h[i-1][j][B2] ;
		B2pm = p_h[i][j-1][B2] ;
		B2pp = p_h[i][j][B2] ;

		alpha = dx1/dt ;	/* crude approx */

		emf[i][j] = 0.5*(F1[i][j][B2] + F1[i][j-1][B2]
				- F2[i][j][B1] - F2[i-1][j][B1]) 

			    - 0.25*(emfmp + emfmm + emfpm + emfpp) 

			    - 0.125*alpha*(
			       B1d - B1mm - B1u + B1mp
			     + B1d - B1pm - B1u + B1pp
			     + B2r - B2pm - B2l + B2mm
			     + B2r - B2pp - B2l + B2mp) ;


	}
#endif

        for(i=0;i<N1+1;i++)
        for(j=0;j<N2;j++) {
                F1[i][j][B1] = 0. ;
                F1[i][j][B2] = 0.5*(emf[i][j] + emf[i][j+1]) ;
        }
        for(i=0;i<N1;i++)
        for(j=0;j<N2+1;j++) {
                F2[i][j][B1] = -0.5*(emf[i][j] + emf[i+1][j]) ;
                F2[i][j][B2] = 0. ;
        }

}
