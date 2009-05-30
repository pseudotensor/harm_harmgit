#include "decs.h"

/* select flux scheme
 *
 * current options:
 * HLLF: Harten, Lax, Van Leer
 * LAXF: local lax friedrichs
 * 
 * both time and space are 4-th order accurate
 *
 */

 //direct reconstruction without integration//
#define HLLF    

void step_ch()
{
        double ndt ;
	int i;


        ZLOOP  primtoU(p[i],uh[i]);

	/* quater step */
        ndt = advance(p, p, 0.5*dt, ph1,uh,1./6.*dt) ;
	/* half step */
        ndt = advance(p, ph1, 0.5*dt, ph2,uh,1./3.*dt) ;
	ndt = advance(p, ph2, dt, ph3,uh,1./3.*dt);
	ndt = advance(p, ph3, dt, ph4,uh,1./6.*dt);
  
	ZLOOP Utoprim(uh[i],p[i]);




	bound_prim(p);
        /* increment time */
        t += dt ;
        /* set next timestep */
        if(ndt > SAFE*dt) ndt = SAFE*dt ; /* don't increase too quickly... */
        dt = ndt ;
        if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */
        /* done! */
}

double advance(
        double pi[][NP],
        double pb[][NP],
        double Dt,
        double pf[][NP],
	double unew[][NP],
	double dt2
        )
{
        int i,k ;
        double ndt,U[NP] ;

        ndt = fluxcalc(pb, F) ;

        /** now update pi to pf **/
        ZLOOP {
                primtoU(pi[i],U) ;

                PLOOP {
                        U[k] += Dt*(
                                - (F[i+1][k] - F[i][k])/dx
				/* could insert source terms here */
                                ) ;
			unew[i][k] += dt2*(-(F[i+1][k] - F[i][k])/dx);
                }

                Utoprim(U,pf[i]) ;
        }
        bound_prim(pf) ;

        return(ndt) ;
}


double fluxcalc(
        double pr[][NP],
        double F[][NP]
        )
{
        int i,k ;
        double p_l[NP],p_r[NP],F_l[NP],F_r[NP],U_l[NP],U_r[NP];
        double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dti ;
        double ctop ;
	double fqa,fqb;
	double fqc;
	double wj;

        /** get left and right face value of every zone **/ 
       ZSLOOP(-1,NX) 
                PLOOP {
                        para(pr[i-2][k], pr[i-1][k],pr[i][k],pr[i+1][k],pr[i+2][k],&Ip_l[i][k],&Ip_r[i][k];
			       }

        ndt = 1.e9 ;

        
  
        ZSLOOP(0,NX) {

  	        PLOOP{
		  //$$$$directly reconstruction without integration$$$$$$$//  
		       p_l[k]=Ip_r[i-1][k];
		       p_r[k]=Ip_l[i][k];
		}

                primtoflux(p_l,F_l) ;
                primtoflux(p_r,F_r) ;

                primtoU(p_l,U_l) ;
                primtoU(p_r,U_r) ;

                vchar(p_l,&cmax_l,&cmin_l) ;
                vchar(p_r,&cmax_r,&cmin_r) ;

                cmax = fabs(MAX(MAX(0.,cmax_l), cmax_r )) ;
                cmin = fabs(MAX(MAX(0.,-cmin_l), -cmin_r)) ;
                ctop = MAX(cmax,cmin) ;

                PLOOP F[i][k] =
#ifdef HLLF
                        (cmax*F_l[k] + cmin*F_r[k]
                                - cmax*cmin*(U_r[k] - U_l[k]))/
                                        (cmax + cmin + SMALL)
                         ;
#endif
#ifdef LAXF
                        0.5*(F_l[k] + F_r[k]
                                - ctop*(U_r[k] - U_l[k]))
                        ;
#endif

                /* evaluate restriction on timestep */
                cmax = MAX(cmax,cmin) ;
                dti = cour*dx/cmax ;
                if(dti < ndt) ndt = dti ;
        }

        return(ndt) ;

}


/* calculate characteristic velocities */
void vchar(double *pr, double *cmax, double *cmin)
{
	double cssq, vaxsq, vasq, vfast ;

	cssq = gam*(gam - 1.)*pr[UU]/pr[RHO] ;
	vaxsq = pr[BX]*pr[BX]/pr[RHO] ;
	vasq = (pr[BX]*pr[BX] + pr[BY]*pr[BY] + pr[BZ]*pr[BZ])/pr[RHO] ;

	vfast = sqrt(
		(cssq + vasq 
		 + sqrt((cssq + vasq)*(cssq + vasq) - 4.*cssq*vaxsq)
		 )/2.) ;

	*cmax = pr[UX] + vfast ;
	*cmin = pr[UX] - vfast ;

	return ;
}

