#include "decs.h"
#include "metric.h"		// defines all coordinate metrics and
				// MCOORD

// this file includes metric dependent terms, including for initial
// condition routines for IC coords.

// only needs to be modified if more complex subterms need to be
// defined
// otherwise define metric in metric.h


// /////////////////////////////////////////////////////////////
// 
// Set metric/coord system here!
// 
// /////////////////////////////////////////////////////////////

void gcov_func(FTYPE *X, FTYPE gcov[][NDIM])
{
  int j, k;
  FTYPE sth, cth, s2, rho2;
  FTYPE r, th;
  FTYPE tfac, rfac, hfac, pfac;
  FTYPE dxdxp[NDIM];

  DLOOP gcov[j][k] = 0.;

  bl_coord(X, &r, &th);

  cth = cos(th);
  if(POSDEFMETRIC){
    sth = fabs(sin(th));
  }
  else{
    sth = sin(th);
  }
  if (fabs(sth) < SMALL)
    {
      if(sth>=0) sth=SMALL;
      if(sth<0) sth=-SMALL;
    }
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  // dx/dx' where '=prim coords (i.e. nonuni coords)

  dxdxprim(X, r, th, dxdxp);

  tfac = dxdxp[0];
  rfac = dxdxp[1];
  hfac = dxdxp[2];
  pfac = dxdxp[3];

  // now take term by term:
  // g_{u v} = \vec{e_{\mu}}\cdot\vec{e_{\nu}} 
  //           * (dx/dx')_{mu} * (dx/dx')_{\nu} =
  //          \vec{e'_{\mu}}\cdot\vec{e'_{\nu}} 

  gcov[TT][TT] = gcov00 * tfac * tfac;
  gcov[TT][1] = gcov01 * tfac * rfac;
  gcov[TT][2] = gcov02 * tfac * hfac;
  gcov[TT][3] = gcov03 * tfac * pfac;

  gcov[1][TT] = gcov10 * rfac * tfac;
  gcov[1][1] = gcov11 * rfac * rfac;
  gcov[1][2] = gcov12 * rfac * hfac;
  gcov[1][3] = gcov13 * rfac * pfac;

  gcov[2][TT] = gcov20 * hfac * tfac;
  gcov[2][1] = gcov21 * hfac * rfac;
  gcov[2][2] = gcov22 * hfac * hfac;
  gcov[2][3] = gcov23 * hfac * pfac;

  gcov[3][TT] = gcov30 * pfac * tfac;
  gcov[3][1] = gcov31 * pfac * rfac;
  gcov[3][2] = gcov32 * pfac * hfac;
  gcov[3][3] = gcov33 * pfac * pfac;
}

// ///////////////////////////////////////////////////////
// 
// Boyer-Lindquist ("bl") metric functions */
// The below functions used for bl coords only (r,th,phi)
// 
// bl coords are starting point for most IC and other metric/coords
// 
// ///////////////////////////////////////////////////////////

// find the con/cov forms of the bl metric
void blgset(int i, int j, struct of_geom *geom)
{
  FTYPE r, th, X[NDIM];

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  if (th < 0)
    th *= -1.;
  if (th > M_PI)
    th = 2. * M_PI - th;

  geom->i=i;
  geom->j=j;
  geom->p=CENT;
  icurr=i;
  jcurr=j;
  pcurr=CENT;
  geom->g = bl_gdet_func(r, th);
  bl_gcov_func(r, th, geom->gcov);
  bl_gcon_func(r, th, geom->gcon);

}

// find the determinant of the bl metric
FTYPE bl_gdet_func(FTYPE r, FTYPE th)
{
  FTYPE a2, r2;

  a2 = a * a;
  r2 = r * r;
  return (bl_gdet);
}

// find gcov for bl metric
void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM])
{
  int j, k;
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;

  DLOOP gcov[j][k] = 0.;

  sth = sin(th);
  s2 = sth * sth;
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  gcov[TT][TT] = bl_gcov00;
  gcov[TT][3] = bl_gcov03;
  gcov[1][1] = bl_gcov11;
  gcov[2][2] = bl_gcov22;
  gcov[3][TT] = bl_gcov30;
  gcov[3][3] = bl_gcov33;

}

// find gcon for bl metric
void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM])
{
  int j, k;
  FTYPE sth, cth, a2, r2, r3, DD, mu;

  DLOOP gcon[j][k] = 0.;

  if(POSDEFMETRIC){
    sth = fabs(sin(th));
  }
  else{
    sth = sin(th);
  }
  if (fabs(sth) < SMALL) {
    if(sth>=0) sth=SMALL;
    if(sth<0) sth=-SMALL;
  }
  
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  gcon[TT][TT] = bl_gcon00;
  gcon[TT][3] = bl_gcon03;
  gcon[1][1] = bl_gcon11;
  gcon[2][2] = bl_gcon22;
  gcon[3][TT] = bl_gcon30;
  gcon[3][3] = bl_gcon33;

}

// ///////////////////////////////////////////////////////////
// 
// below are independent of user choice of metric/coords/grid
// 
// ///////////////////////////////////////////////////////////


// find determinant in general of a metric
/* assumes gcov has been set first; returns determinant */
FTYPE gdet_func(FTYPE gcov[][NDIM])
{
  static int firstc = 1;
  static FTYPE **tmp;
  FTYPE d;
  int j, k, indx[NDIM];

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP tmp[j + 1][k + 1] = gcov[j][k];
  ludcmp(tmp, NDIM, indx - 1, &d);
  // below from 1..NDIM due to ludcmp requiring 1..N
  for (j = 1; j <= NDIM; j++)
    d *= tmp[j][j];

  return (sqrt(fabs(d)));
}

/* invert gcov to get gcon */
void gcon_func(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  static int firstc = 1;
  int j, k;
  static FTYPE **tmp;

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP tmp[j + 1][k + 1] = gcov[j][k];
  gaussj(tmp, NDIM, NULL, 0);
  DLOOP gcon[j][k] = tmp[k + 1][j + 1];
}

/* 
   this gives the connection coefficient \Gamma^{i}_{j,k} =
   conn[..][i][j][k] where i,j,k = {0,1,2,3} corresponds to {t,r,theta,phi} 
 */

/*
  we also compute the 2nd connection:
  -d/dj(ln(\detg))
*/

// delta is simply how big the differencing is, should be small, but not so small to lead to errors due to erros in the metric itself (i.e. keep larger than machine precision)
//#define DELTA (NUMEPSILON*1000.0)
//#define DELTA 1.e-5
#define DELTA 1.e-5
// how to generically set this?  Too high, even slightly (10^{-10} for long doubles) and connection is screwed)

/* NOTE: parameter hides global variable */
void conn_func(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE tmp[NDIM][NDIM][NDIM];
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE gh[NDIM][NDIM];
  FTYPE gl[NDIM][NDIM];
  FTYPE lngdeth,lngdetl;

  // gabc_{ijk}=dg_{ij}/dx^k
  for (k = 0; k < NDIM; k++) {
    for (l = 0; l < NDIM; l++)
      Xh[l] = X[l];
    for (l = 0; l < NDIM; l++)
      Xl[l] = X[l];
    Xh[k] += DELTA;
    Xl[k] -= DELTA;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);
    lngdeth=log(gdet_func(gh));
    lngdetl=log(gdet_func(gl));

    conn2[k]= - (lngdeth - lngdetl) / (Xh[k] - Xl[k]);

    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
	conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
  }

  /* now rearrange to find \Gamma_{ijk}=1/2*(gabc_{jik}+gabc_{kij}-gabc_{kji}) */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
	tmp[i][j][k] =
	    0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  /* finally, raise first index */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++) {
	conn[i][j][k] = 0.;
	for (l = 0; l < NDIM; l++)
	  conn[i][j][k] += geom->gcon[i][l] * tmp[l][j][k];
      }

  /* done! */
}

#undef DELTA



/* 
   FTYPE delta(int i, int j) { if(i == j) return(1.) ; else return(0.) 
   ; } */

/* Minkowski metric; signature +2 */
/* 
   FTYPE mink(int i, int j) { if(i == j) { if(i == 0) return(-1.) ;
   else return(1.) ; } else return(0.) ; } */

FTYPE cot(FTYPE arg)
{
  return(1.0/tan(arg));
}

// jon's MKS connection (and conn2)
void mks_conn_func(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE r,th,sigma,dxdxp[NDIM];
  FTYPE cot(FTYPE arg);

  // get bl coordinates
  bl_coord(X,&r,&th);
  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
  }
  else{
    if(th>M_PI) { th=M_PI-th; }
  }
  // avoid singularity at polar axis
  if(fabs(th)<SMALL){
    if(th>=0) th=SMALL;
    if(th<0) th=-SMALL;
  }
  if(fabs(M_PI-th)<SMALL){
    if(th>=M_PI) th=M_PI+SMALL;
    if(th<M_PI) th=M_PI-SMALL;
  }

  // set aux vars
  dxdxprim(X,r,th,dxdxp);
  sigma=r*r+a*a*cos(th)*cos(th);


  conn[0][0][0]=(-2.*r*sigma + 4.*pow(r,3.))*pow(sigma,-3.);
conn[0][0][1]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[0][0][2]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
conn[0][0][3]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[0][1][0]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[0][1][1]=2.*(r + sigma)*pow(dxdxp[1],2.)*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
conn[0][1][2]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
conn[0][1][3]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
conn[0][2][0]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
conn[0][2][1]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
conn[0][2][2]=-2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-1.);
conn[0][2][3]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
conn[0][3][0]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[0][3][1]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
conn[0][3][2]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
conn[0][3][3]=2.*r*pow(sigma,-3.)*pow(sin(th),2.)*
    (-1.*r*pow(sigma,2.) + pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
       pow(sin(th),2.));
conn[1][0][0]=pow(dxdxp[1],-1.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
conn[1][0][1]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
conn[1][0][2]=0.;
conn[1][0][3]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
conn[1][1][0]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
conn[1][1][1]=pow(sigma,-3.)*
    (-1.*dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.)) + 
      pow(sigma,3.) + dxdxp[1]*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
       pow(sin(th),2.));
conn[1][1][2]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
   ;
conn[1][1][3]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
      cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
      2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[1][2][0]=0.;
conn[1][2][1]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
   ;
conn[1][2][2]=-1.*r*pow(dxdxp[1],-1.)*pow(dxdxp[2],2.)*
    pow(sigma,-1.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
conn[1][2][3]=0.;
conn[1][3][0]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
conn[1][3][1]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
      cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
      2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[1][3][2]=0.;
conn[1][3][3]=-1.*pow(dxdxp[1],-1.)*pow(sigma,-3.)*pow(sin(th),2.)*
    (-2.*r + sigma + pow(a,2.)*pow(sin(th),2.))*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[2][0][0]=-1.*r*pow(dxdxp[2],-1.)*pow(a,2.)*pow(sigma,-3.)*
    sin(2.*th);
conn[2][0][1]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][0][2]=0.;
conn[2][0][3]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
conn[2][1][0]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][1][1]=-1.*r*pow(dxdxp[1],2.)*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
conn[2][1][2]=dxdxp[1]*r*pow(sigma,-1.);
conn[2][1][3]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
conn[2][2][0]=0.;
conn[2][2][1]=dxdxp[1]*r*pow(sigma,-1.);
conn[2][2][2]=4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*
     pow(M_PI,2.) - 1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)\
    ;
conn[2][2][3]=0.;
conn[2][3][0]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
conn[2][3][1]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
conn[2][3][2]=0.;
conn[2][3][3]=-1.*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (pow(sigma,3.) + sigma*(4.*r + sigma)*pow(a,2.)*pow(sin(th),2.) + 
      2.*r*pow(a,4.)*pow(sin(th),4.))*sin(th);
conn[3][0][0]=a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
conn[3][0][1]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
   ;
conn[3][0][2]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
conn[3][0][3]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[3][1][0]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
   ;
conn[3][1][1]=a*pow(dxdxp[1],2.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
conn[3][1][2]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
conn[3][1][3]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[3][2][0]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
conn[3][2][1]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
conn[3][2][2]=-1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-1.);
conn[3][2][3]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
conn[3][3][0]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
conn[3][3][1]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
conn[3][3][2]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
conn[3][3][3]=pow(sigma,-3.)*
    (-1.*a*r*pow(sigma,2.)*pow(sin(th),2.) + 
      pow(a,3.)*(-1.*sigma + 2.*pow(r,2.))*pow(sin(th),4.));
conn2[0]=0.;
conn2[1]=-1.*pow(sigma,-1.)*(2.*dxdxp[1]*r + pow(r,2.) + 
      pow(a,2.)*pow(cos(th),2.));
conn2[2]=-1.*dxdxp[2]*cot(th) + 
    4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
    dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th);
conn2[3]=0.;


}
