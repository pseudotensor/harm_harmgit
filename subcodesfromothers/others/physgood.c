
#include "decs.h"

/* calculate fluxes in direction dir and conserved variable U; these
   are always needed together, so there is no point in calculated the
   stress tensor twice */

int primtoflux(FTYPE *pr, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *flux)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  int i = 0, j = 0, k = 0;
  FTYPE mhd[NDIM];

  /* particle number flux */
  flux[RHO] = pr[RHO] * q->ucon[dir];

  // GODMARK WTF!
  // if(mhd_calc(pr,dir,q,mhd)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);
  mhd_calc(pr, dir, q, mhd);

  /* MHD stress-energy tensor w/ first index up, second index down. */
  flux[UU] = mhd[0] + flux[RHO];
  flux[U1] = mhd[1];
  flux[U2] = mhd[2];
  flux[U3] = mhd[3];

  /* dual of Maxwell tensor */
  flux[B1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  flux[B2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  flux[B3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];

  PLOOP flux[k] *= geom->g;
#if(WHICHEOM==WITHNOGDET)
  flux[U1]/=geom->g;
  flux[U2]/=geom->g;
#endif
  return (0);
}

/* calculate "conserved" quantities */
int primtoU(FTYPE *pr, struct of_state *q, struct of_geom *geom,
	    FTYPE *U)
{
  int i = 0, j = 0, k = 0;
  MYFUN(primtoflux(pr, q, 0, geom, U) ,"phys.c:primtoU()", "primtoflux_calc() dir=0", 1);

  return (0);
}


/* calculate magnetic field four-vector */
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j <= 3; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}

/* MHD stress tensor, with first index up, second index down */
void mhd_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd)
{
  int j;
  FTYPE r, u, P, w, bsq, eta, ptot;

  r = pr[RHO];
  u = pr[UU];
  P = (gam - 1.) * u;
  w = P + r + u;
  bsq = dot(q->bcon, q->bcov);
  eta = w + bsq;
  ptot = P + bsq*0.5;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  DLOOPA mhd[j] = eta * q->ucon[dir] * q->ucov[j]
      + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];

}

/* add in source terms to equations of motion */
int source(FTYPE *ph, struct of_geom *ptrgeom, int ii, int jj,
	   FTYPE *dU, SFTYPE Dt)
{
  FTYPE mhd[NDIM][NDIM];
  int i = 0, j = 0, k = 0;
  struct of_state q;

  MYFUN(get_state(ph, ptrgeom, &q) ,"phys.c:source()", "get_state() dir=0", 1);
  mhd_calc(ph, 0, &q, mhd[0]);
  mhd_calc(ph, 1, &q, mhd[1]);
  mhd_calc(ph, 2, &q, mhd[2]);
  mhd_calc(ph, 3, &q, mhd[3]);

  /* contract mhd stress tensor with connection */
  PLOOP dU[k] = 0.;
  // mhd^{dir}_{comp} = mhd^j_k
  // dU^{a} = mhd^j_k C^{k}_{a j}
  DLOOP {
    dU[UU] += mhd[j][k] * conn[ii][jj][k][0][j];
    dU[U1] += mhd[j][k] * conn[ii][jj][k][1][j];
    dU[U2] += mhd[j][k] * conn[ii][jj][k][2][j];
    dU[U3] += mhd[j][k] * conn[ii][jj][k][3][j];
  }

#if(WHICHEOM==WITHNOGDET)
  DLOOPA {
    dU[U1] += mhd[j][1] * conn2[ii][jj][j];
    dU[U2] += mhd[j][2] * conn2[ii][jj][j];
  }
#endif

  /* cooling */
  if(cooling){
    dU[UU] += coolfunc(h_over_r, ph, ptrgeom, &q);
  }
  //misc_source(ph, ii, jj, geom, &q, dU, Dt) ;


  PLOOP dU[k] *= ptrgeom->g;
#if(WHICHEOM==WITHNOGDET)
  dU[U1]/=(ptrgeom->g);
  dU[U2]/=(ptrgeom->g);
#endif

  /* done! */
  return (0);
}

/* returns b^2 (i.e., twice magnetic pressure) */
int bsq_calc(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  int i = 0, j = 0, k = 0;
  struct of_state q;

  MYFUN(get_state(pr, ptrgeom, &q) ,"phys.c:bsq_calc()", "get_state() dir=0", 1);
  *bsq = dot(q.bcon, q.bcov);
  return (0);
}


void lower(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov)
{
        ucov[0] = geom->gcov[0][0]*ucon[0]
                + geom->gcov[0][1]*ucon[1]
                + geom->gcov[0][2]*ucon[2]
	  + geom->gcov[0][3]*ucon[3] ;
        ucov[1] = geom->gcov[0][1]*ucon[0]
                + geom->gcov[1][1]*ucon[1]
                + geom->gcov[1][2]*ucon[2]
	  + geom->gcov[1][3]*ucon[3] ;
        ucov[2] = geom->gcov[0][2]*ucon[0]
                + geom->gcov[1][2]*ucon[1]
                + geom->gcov[2][2]*ucon[2]
	  + geom->gcov[2][3]*ucon[3] ;
        ucov[3] = geom->gcov[0][3]*ucon[0]
                + geom->gcov[1][3]*ucon[1]
                + geom->gcov[2][3]*ucon[2]
	  + geom->gcov[3][3]*ucon[3] ;

        return ;
}

void lowerf(FTYPE *fcon, struct of_geom *geom, FTYPE *fcov)
{
  int j,k;
  int jp,kp;
  FTYPE myfcon[NDIM][NDIM],myfcov[NDIM][NDIM];

  myfcon[0][0]=myfcon[1][1]=myfcon[2][2]=myfcon[3][3]=0;
  myfcon[0][1]=fcon[0];
  myfcon[0][2]=fcon[1];
  myfcon[0][3]=fcon[2];
  myfcon[1][2]=fcon[3];
  myfcon[1][3]=fcon[4];
  myfcon[2][3]=fcon[5];
  //
  myfcon[1][0]=-fcon[0];
  myfcon[2][0]=-fcon[1];
  myfcon[3][0]=-fcon[2];
  myfcon[2][1]=-fcon[3];
  myfcon[3][1]=-fcon[4];
  myfcon[3][2]=-fcon[5];
  
  DLOOP{
    myfcov[j][k]=0;
    for(jp=0;jp<NDIM;jp++) for(kp=0;kp<NDIM;kp++){
      myfcov[j][k]+=myfcon[jp][kp]*geom->gcov[j][jp]*geom->gcov[k][kp];
    }
  }
  fcov[0]=myfcov[0][1];
  fcov[1]=myfcov[0][2];
  fcov[2]=myfcov[0][3];
  fcov[3]=myfcov[1][2];
  fcov[4]=myfcov[1][3];
  fcov[5]=myfcov[2][3];

  return ;
}

void raise(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon)
{

        ucon[0] = geom->gcon[0][0]*ucov[0]
                + geom->gcon[0][1]*ucov[1]
                + geom->gcon[0][2]*ucov[2]
	  + geom->gcon[0][3]*ucov[3] ;
        ucon[1] = geom->gcon[0][1]*ucov[0]
                + geom->gcon[1][1]*ucov[1]
                + geom->gcon[1][2]*ucov[2]
	  + geom->gcon[1][3]*ucov[3] ;
        ucon[2] = geom->gcon[0][2]*ucov[0]
                + geom->gcon[1][2]*ucov[1]
                + geom->gcon[2][2]*ucov[2]
	  + geom->gcon[2][3]*ucov[3] ;
        ucon[3] = geom->gcon[0][3]*ucov[0]
                + geom->gcon[1][3]*ucov[1]
                + geom->gcon[2][3]*ucov[2]
	  + geom->gcon[3][3]*ucov[3] ;

        return ;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
int get_state(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int i = 0, j = 0, k = 0;

  /* get ucon */
  
  MYFUN(ucon_calc(pr, ptrgeom, q->ucon) ,"phys.c:get_state()", "ucon_calc()", 1);
  lower(q->ucon, ptrgeom, q->ucov);
  bcon_calc(pr, q->ucon, q->ucov, q->bcon);
  lower(q->bcon, ptrgeom, q->bcov);
  return (0);
}

/* load local geometry into structure geom */
void get_geometry(int ii, int jj, int kk, struct of_geom *geom)
{
  int j, k;

  //  DLOOP geom->gcov[j][k] = gcov[ii][jj][kk][j][k];
  //DLOOP geom->gcon[j][k] = gcon[ii][jj][kk][j][k];
  // let's vectorize it
  for(j=0;j<=NDIM*NDIM-1;j++){
    geom->gcon[0][j] = gcon[ii][jj][kk][0][j];
    geom->gcov[0][j] = gcov[ii][jj][kk][0][j];
  }
  geom->g = gdet[ii][jj][kk];
  geom->i = ii;
  geom->j = jj;
  geom->p = kk;
#if(JONCHECKS)
  icurr = ii;
  jcurr = jj;
  pcurr = kk;
#endif
}


/* find contravariant four-velocity from the relative 4 velocity */
int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
        FTYPE alpha,gamma ;
        FTYPE beta[NDIM] ;
        int j ;

        alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
        SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ;

        MYFUN(gamma_calc(pr,geom,&gamma),"ucon_calc_rel4vel: gamma calc failed\n","phys.c",1);

        ucon[TT] = gamma/alpha ;
        SLOOPA ucon[j] = pr[U1+j-1] - gamma*beta[j]/alpha ;

        return(0) ;
}

/* find gamma-factor wrt normal observer */
int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE*gamma)
{
        FTYPE qsq ;
	int j,k;

        qsq =     geom->gcov[1][1]*pr[U1]*pr[U1]
                + geom->gcov[2][2]*pr[U2]*pr[U2]
                + geom->gcov[3][3]*pr[U3]*pr[U3]
            + 2.*(geom->gcov[1][2]*pr[U1]*pr[U2]
                + geom->gcov[1][3]*pr[U1]*pr[U3]
                + geom->gcov[2][3]*pr[U2]*pr[U3]) ;

#if(JONCHECKS2)
        if(qsq<0.0){
	  if(fabs(qsq)>1E-10){ // then assume not just machine precision
	    dualfprintf(fail_file,"gamma_calc failed: qsq=%g\n",qsq);
	    PLOOP dualfprintf(fail_file,"pr[%d]=%21.15g\n",k,pr[k]);
	    DLOOP dualfprintf(fail_file,"gcov[%d][%d]\n",j,k,gcov[j][k]);
	    if (fail(FAIL_UTCALC_DISCR) >= 1)
	      return (1);
	  }
	  else qsq=1E-10; // set floor
	}
#endif
        *gamma = sqrt(1. + qsq) ;

        return(0) ;
}





/* find contravariant four-velocity */
//int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  FTYPE discr;
  // debug stuff
  int j,k;
  FTYPE r,th,X[NDIM];


  ucon[1] = pr[U1];
  ucon[2] = pr[U2];
  ucon[3] = pr[U3];

  discr = geom->gcov[0][0]
      + geom->gcov[1][1] * ucon[1] * ucon[1]
      + geom->gcov[2][2] * ucon[2] * ucon[2]
      + geom->gcov[3][3] * ucon[3] * ucon[3]
      + 2. * (geom->gcov[0][1]* ucon[1]
	      + geom->gcov[0][2] * ucon[2]
	      + geom->gcov[0][3] * ucon[3]
	      + geom->gcov[1][2] * ucon[1] * ucon[2]
	      + geom->gcov[1][3] * ucon[1] * ucon[3]
	      + geom->gcov[2][3] * ucon[2] * ucon[3]);

  uttdiscr=-discr;

  if (discr > 0.) {
#if(JONCHECKS)
    if(whocalleducon==0){
      // then report on disc
      dualfprintf(fail_file,"disc=%21.15g, should be negative\n",discr);
      PLOOP{
	dualfprintf(fail_file,"uconfailed on pr[%d]=%21.15g\n",k,pr[k]);
      }
      coord(geom->i,geom->j,geom->p,X);
      bl_coord(X,&r,&th);
      dualfprintf(fail_file,"i=%d j=%d pcurr=%d\nx1=%21.15g x2=%21.15g \nr=%21.15g th=%21.15g \ng=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,geom->p,X[1],X[2],r,th,geom->g);
      dualfprintf(fail_file,"\ngcon\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcon[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
      dualfprintf(fail_file,"\ngcov\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcov[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
    }
#endif
    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  ucon[TT] = 1. / sqrt(-discr);
  ucon[1] *= ucon[TT];
  ucon[2] *= ucon[TT];
  ucon[3] *= ucon[TT];

  return (0);
}



/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
	FTYPE AA,BB,CC ;
	FTYPE discr ;
	FTYPE bsq,X[NDIM] ;
	int i=0,j=0,k=0 ;

        ucon[1] = pr[U1] ;
        ucon[2] = pr[U2] ;
        ucon[3] = pr[U3] ;

	AA = geom->gcov[TT][TT] ;
	BB = 2.*(geom->gcov[TT][1]*ucon[1] +
		 geom->gcov[TT][2]*ucon[2] +
		 geom->gcov[TT][3]*ucon[3]) ;
	CC = 1. +
	  geom->gcov[1][1]*ucon[1]*ucon[1] +
	  geom->gcov[2][2]*ucon[2]*ucon[2] +
	  geom->gcov[3][3]*ucon[3]*ucon[3] +
	  2.*(geom->gcov[1][2]*ucon[1]*ucon[2] +
	      geom->gcov[1][3]*ucon[1]*ucon[3] +
	      geom->gcov[2][3]*ucon[2]*ucon[3]) ;
	
	discr = BB*BB - 4.*AA*CC ;
	if(discr < 0.) {
		/*
		fprintf(fail_file,"failure %d %d\n",icurr,jcurr) ;
        	ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
        	ucon[TT] = -BB/(2.*AA) ;
		*/
		fprintf(fail_file,"failure: spacelike four-velocity %g\n",
			discr) ;
		fprintf(fail_file,"%d %d %d\n",geom->i,geom->j,geom->p) ;
		coord(geom->i,geom->j,geom->p,X);
		fprintf(fail_file,"%g %g\n",X[1],X[2]) ;

		/*
		  if(bsq_calc(&bsq,pr[icurr][jcurr])>=1) FAILSTATEMENT("phys.c:ucon_calc()","bsq_calc()",1);

		fprintf(fail_file,"bsq/rho: %g\n",bsq/pr[icurr][jcurr][0]) ;

		*/
		for(k=0;k<NPR;k++) fprintf(fail_file,"%d %15.8g\n",k,pr[k]) ;
		failed=1;
		return(1);
	}

	ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
	return(0) ;
}


FTYPE taper_func(FTYPE R,FTYPE rin)
{

  if(R <= rin)
    return(0.) ;
  else
    return(1. - sqrt(rin/R)) ;

}

// compute the radius of the inner most stable circular orbit
FTYPE rmso_calc(int which)
{
  FTYPE rmso,Z1,Z2,sign ;


  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
                                 pow(1. - a, 1./3.)) ;
  Z2 = sqrt(3.*a*a + Z1*Z1) ;
  rmso=3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;

  return(rmso) ;
}

FTYPE uphi_isco_calc(int which,FTYPE r)
{
  FTYPE uphi;
  FTYPE sign;
  FTYPE Z1,Z2;

  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1=r*r-sign*2.*a*sqrt(r)+a*a;
  Z2=r*(r*r-3.*r+sign*2.*a*sqrt(r));

  uphi=sign*Z1/sqrt(Z2);

  return(uphi);

}

FTYPE rhor_calc(int which)
{
  FTYPE sign;
  if(which==0) sign=1; else sign=-1;

  return(1. +sign*sqrt(1. - a * a));
}




// used Mathematica's MinimumChangePermutations and Signature
FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda)
{
  int i;
  FTYPE lc4sign;
  int l1[24]={1, 2, 3, 1, 2, 3, 4, 2, 1, 4, 2, 1, 1, 3, 4, 1, 3, 4, 4, 3, 2, 4, 3, 2};
  int l2[24]={2, 1, 1, 3, 3, 2, 2, 4, 4, 1, 1, 2, 3, 1, 1, 4, 4, 3, 3, 4, 4, 2, 2, 3};
  int l3[24]={3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 4, 4, 4, 4, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4};
  int l4[24]={4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};

  for(i=0;i<24;i++){
    if((mu==l1[i])&&(nu==l2[i])&&(kappa==l3[i])&&(lambda==l4[i])){
      lc4sign=(i%2) ? -1 : 0;
      if(updown==1) return(-1.0/detg*lc4sign); // up
      else if(updown==0) return(detg*lc4sign); // down
    }
  }
  // if didn't get here, then 0
  return(0.0);
}

void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM])
{
  int nu,mu,kappa,lambda;

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    faraday[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      faraday[mu][nu]+=lc4(which,geom->g,mu,nu,kappa,lambda)*u[kappa]*b[lambda];
    }
  }

}


void current_doprecalc(int which, FTYPE p[][N2+4][NPR])
{
  int i,j;
  int idel, jdel;
  struct of_geom geom;
  struct of_state q;
  FTYPE Dt;
  int face;

  if(which==0){ face=CENT; idel=0; jdel=0; Dt=dt*0.5; }
  else if(which==1){ face=FACE1; idel=1; jdel=0; Dt=dt; }
  else if(which==2){ face=FACE2; idel=0; jdel=1; Dt=dt; }
  else if(which==3){ face=CENT; idel=0; jdel=0; Dt=dt; }

  FZLOOP(-jdel, -idel) {
    get_geometry(i, j, face, &geom);
    MYFUN(get_state(p[i][j], &geom, &q),"phys.c:current_doprecalc()", "get_state()", 1);
    current_precalc(which,&geom,&q,Dt,cfaraday[i][j]);
  }

}

// which: 0: somewhere at half step
// which: 1: doing flux calculation in x1 direction, full step
// which: 2: doing flux calculation in x2 direction, full step
void current_precalc(int which, struct of_geom *geom, struct of_state *q, FTYPE Dt,FTYPE faraday[][3])
{
  // assume outside loop is like flux, from 0..N in r for r, 0..N in h for h.  And must wait till 2nd timestep before computing the current since need time differences
  if(which==0){
    // assume got here when DT==dt/2 and geom and state set at zone center
    // first save old calculation
    faraday[3][0]=faraday[0][0];
    faraday[3][1]=faraday[0][1];
    faraday[3][2]=faraday[0][2];
    // now calculate new version
    faraday[0][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt}
    faraday[0][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[1]-q->bcov[1]*q->ucov[3]); // f^{ht}
    faraday[0][2]=-1.0/geom->g * (-q->bcov[2]*q->ucov[1]+q->bcov[1]*q->ucov[2]); // f^{pt}
  }
  else if(which==1){
    // assume got here with DT=dt and geom and state at radial zone edge
    faraday[1][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}=-f^{rt}
    faraday[1][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[0]+q->bcov[0]*q->ucov[3]); // f^{hr}
    faraday[1][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{pr}
  }
  else if(which==2){
    // assume got here with DT=dt and geom and state at theta zone edge
    faraday[2][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}=-f^{ht}
    faraday[2][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}=-f^{hr}
    faraday[2][2]=-1.0/geom->g * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{ph}
  }
  else if(which==3){
    // DT==dt, but zone center
    fcon[geom->i][geom->j][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}
    fcon[geom->i][geom->j][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}
    fcon[geom->i][geom->j][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[1]-q->bcov[1]*q->ucov[2]); // f^{tp}
    fcon[geom->i][geom->j][3]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}
    fcon[geom->i][geom->j][4]=-1.0/geom->g * (-q->bcov[2]*q->ucov[0]+q->bcov[0]*q->ucov[2]); // f^{rp}
    fcon[geom->i][geom->j][5]=-1.0/geom->g * (q->bcov[1]*q->ucov[0]-q->bcov[0]*q->ucov[1]); // f^{hp}
  }
}

// the current is calculated to end up at the zone and time edge
void current_calc(FTYPE cfaraday[][N2+4][4][3])
{
  int i,j;
  struct of_geom geomt;
  struct of_geom geomr;
  struct of_geom geomh;
  struct of_geom geomrp1;
  struct of_geom geomhp1;
  static FTYPE lastdt;
  static int calls=0;
  FTYPE idtc,idx1,idx2;

  if(calls>0){ // since need 2 times

    idtc=2.0/(lastdt+dt);
    idx1=1.0/dx[1];
    idx2=1.0/dx[2];

    ZLOOP{
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,FACE1,&geomr);
      get_geometry(i,j,FACE2,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,FACE1,&geomrp1);
      get_geometry(i,j+1,FACE2,&geomhp1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomr.g*cfaraday[i][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomh.g*cfaraday[i][j][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      jcon[i][j][1]=
	(cfaraday[i][j][0][0]-cfaraday[i][j][3][0])*idtc+ // F^{rt},t
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomh.g*cfaraday[i][j][2][1])*idx2; // F^{rh},h
      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	(cfaraday[i][j][0][1]-cfaraday[i][j][3][1])*idtc+ // F^{ht},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i][j][1][1])*idx1; // F^{hr},r
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	(cfaraday[i][j][0][2]-cfaraday[i][j][3][2])*idtc+ // F^{pt},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomr.g*cfaraday[i][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomh.g*cfaraday[i][j][2][2])*idx2; // F^{ph},h
    }
  }
  calls++;
  lastdt=dt;

}

