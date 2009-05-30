
/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux */
#define HLLF	1.0
#define LAXF	0.0

#define FLUXCTTOTH 1
#define FLUXCTHLL  0

/** end algorithmic choices **/


double advance( double RK_f1, double pi[][N2+2*NGhost][NPR], double RK_f2, double pb[][N2+2*NGhost][NPR], double Dt, double pf[][N2+2*NGhost][NPR],double unew[][N2+2*NGhost][NPR]) ;
double fluxcalc( double pr[][N2+2*NGhost][NPR], double F[][N2+2*NGhost][NPR], 
	int dir ) ;
void   flux_cd(double F1[][N2+2*NGhost][NPR], double F2[][N2+2*NGhost][NPR]) ;

void step_ch()
{
	double ndt, ndt_1,ndt_2,ndt_3 ;
	int i,j,k ;
	struct of_geom geom;
	struct of_state q;
	
	ZLOOP {
	         get_geometry(i,j,CENT,&geom);
		 get_state(p[i][j],&geom,&q);
		 primtoU(p[i][j],&q,&geom,uh[i][j]);
	}

	
	ndt_1 = advance(1., p, 0., p,dt, ph1,uh) ;
	fprintf(stderr,"_1");

	ndt_2 = advance(0.75, p, 0.25, ph1, 0.25*dt, ph2, uh) ;
	fprintf(stderr,"_2") ;

	ZLOOP PLOOP psave[i][j][k] = p[i][j][k] ;
	ndt_3 = advance(1./3., p, 2./3., ph2, 2./3.*dt, p, uh);
	fprintf(stderr,"_3") ;

	ndt= MIN(MIN(ndt_1,ndt_2),ndt_3);
	

        if(dt < 1.e-9) {
                fprintf(stderr,"timestep too small\n") ;
                exit(11) ;
        }

        /* increment time */
        t += dt ;

        /* set next timestep */
        if(ndt > SAFE*dt) ndt = SAFE*dt ;
        dt = ndt ;
        if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */

        /* done! */
}

double advance(
	double RK_f1,
	double pi[][N2+2*NGhost][NPR], 
	double RK_f2,
	double pb[][N2+2*NGhost][NPR], 
	double Dt,
	double pf[][N2+2*NGhost][NPR],
	double unew[][N2+2*NGhost][NPR]
	)
{
	int i,j,k ;
	double ndt,ndt1,ndt2,U[NPR],dU[NPR],Ui[NPR] ;
	struct of_geom geom ;
	struct of_state q ;

	ZLOOP PLOOP pf[i][j][k] = pi[i][j][k] ;        /* needed for Utoprim */

	fprintf(stderr,"0") ;
	ndt1 = fluxcalc(pb, F1, 1) ;
	ndt2 = fluxcalc(pb, F2, 2) ;

	flux_ct(F1,F2) ;

	/* evaluate diagnostics based on fluxes */
	diag_flux(F1,F2) ;

	fprintf(stderr,"1") ;
	/** now update pi to pf **/
	ZLOOP {

		get_geometry(i,j,CENT,&geom) ;

		source(pb[i][j],&geom,i,j,dU) ;

		get_state(pi[i][j],&geom,&q) ;
		primtoU(pi[i][j],&q,&geom,Ui) ;
		
		PLOOP U[k] = Ui[k]*RK_f1 + unew[i][j][k]*RK_f2;
              
		PLOOP {
                        U[k] += Dt*(
                                - (F1[i+1][j][k] - F1[i][j][k])/dx[1]
                                - (F2[i][j+1][k] - F2[i][j][k])/dx[2]
                                + dU[k]
                                ) ;
			unew[i][j][k] = U[k];
                }

                Utoprim(U,&geom,pf[i][j]) ;
	}

	fixup(pf) ;
        bound_prim(pf) ;

	ndt = defcon * 1./(1./ndt1 + 1./ndt2) ;
	fprintf(stderr,"2") ;

	return(ndt) ;
}



double fluxcalc(
	double pr[][N2+2*NGhost][NPR], 
	double F[][N2+2*NGhost][NPR], 
	int dir 
	)
{
	int i,j,k,idel,jdel,face,ipara,jpara ;
	double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
	double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dtij ;
	double ctop ;
	struct of_geom geom ;
	struct of_state state_l,state_r ;

        if     (dir == 1) {idel = 1; jdel = 0; ipara=2; jpara=0;face = FACE1;}
	else if(dir == 2) {idel = 0; jdel = 1; ipara=0; jpara=2;face = FACE2;}
	else { exit(10); }

	/** evaluate slopes of primitive variables_linear reconstruction **/
	/*
	ZSLOOP(1-NGhost,N1+NGhost-2,1-NGhost,N2+NGhost-2) PLOOP {
			dq[i][j][k] = slope_lim(
					pr[i-idel][j-jdel][k],
					pr[i][j][k],
					pr[i+idel][j+jdel][k]
					) ;
	}

	ndt = 1.e9 ;
       
        ZSLOOP(-jdel,N1,-idel,N2) {
                PLOOP {
                        p_l[k] = pr[i-idel][j-jdel][k] 
					+ 0.5*dq[i-idel][j-jdel][k] ;
                        p_r[k] = pr[i][j][k]   
					- 0.5*dq[i][j][k] ;
                }
	*/

	/*parabolic reconstruction*/
	//	ZSLOOP(1-NGhost,N1+NGhost-2,1-NGhost,N2+NGhost-2) 
	ZSLOOP(2-NGhost,N1+NGhost-3,2-NGhost,N2+NGhost-3) 
	PLOOP {
	              para(pr[i-ipara][j-jpara][k],pr[i-idel][j-jdel][k],pr[i][j][k],pr[i+idel][j+jdel][k],pr[i+ipara][j+jpara][k],&Ip_l[i][j][k],&Ip_r[i][j][k]);
	}
	
	ndt = 1.e9;
	
	ZSLOOP(-jdel,N1,-idel,N2) {
	  
	        PLOOP{
		        p_l[k] = Ip_r[i-idel][j-jdel][k];
		        p_r[k] = Ip_l[i][j][k];
		}

		get_geometry(i,j,face,&geom) ;

		get_state(p_l,&geom,&state_l) ;
		get_state(p_r,&geom,&state_r) ;
		
		primtoflux(p_l,&state_l,dir,&geom,F_l) ;
		primtoflux(p_r,&state_r,dir,&geom,F_r) ;

		primtoflux(p_l,&state_l,TT, &geom,U_l) ;
		primtoflux(p_r,&state_r,TT, &geom,U_r) ;

		vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
		vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;

		cmax = fabs(MAX(MAX(0.,cmax_l), cmax_r )) ;
		cmin = fabs(MAX(MAX(0.,-cmin_l), -cmin_r)) ;
		ctop = MAX(cmax,cmin) ;
		//ctop = 0.;

		PLOOP F[i][j][k] = 
			HLLF*(
			(cmax*F_l[k] + cmin*F_r[k] 
				- cmax*cmin*(U_r[k] - U_l[k]))/
					(cmax + cmin + SMALL) 
			) +
			LAXF*(
			0.5*(F_l[k] + F_r[k] 
				- ctop*(U_r[k] - U_l[k])) 
			) ;

                /* evaluate restriction on timestep */
                cmax = MAX(cmax,cmin) ;
                dtij = cour*dx[dir]/cmax ;
		if(dtij < ndt) ndt = dtij ;
	}

	return(ndt) ;

}



void flux_ct(double F1[][N2+2*NGhost][NPR], double F2[][N2+2*NGhost][NPR])
{
	int i,j ;
	static double emf[N1+1][N2+1] ;

	/* calculate EMFs */
	/* Toth approach: just average */
	ZSLOOP(0,N1,0,N2) emf[i][j] = 0.25*(F1[i][j][B2] + F1[i][j-1][B2]
					  - F2[i][j][B1] - F2[i-1][j][B1]) ;

	/* rewrite EMFs as fluxes, after Toth */
        ZSLOOP(0,N1,0,N2-1) {
                F1[i][j][B1] = 0. ;
                F1[i][j][B2] =  0.5*(emf[i][j] + emf[i][j+1]) ;
        }
        ZSLOOP(0,N1-1,0,N2) {
                F2[i][j][B1] = -0.5*(emf[i][j] + emf[i+1][j]) ;
                F2[i][j][B2] = 0. ;
	}

}

