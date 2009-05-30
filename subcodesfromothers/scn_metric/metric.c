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
  int j, k, jp, kp;
  FTYPE sth, cth, s2, rho2;
  FTYPE r, th;
  FTYPE dx_dxp[NDIM][NDIM];
  static FTYPE gtmp[NDIM][NDIM];

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth=sin(th);

  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  // dx/dx' where '=prim coords (i.e. nonuni coords)

  //-scn  dxdxprim(X, r, th, dxdxp);
  dx_dxp_calc( X, r, th, dx_dxp );

  // Set the metric first that in KS coordinates:
  gtmp[TT][TT] = gcov00  ;  gtmp[TT][1]  = gcov01  ;
  gtmp[TT][2]  = gcov02  ;  gtmp[TT][3]  = gcov03  ;
  gtmp[1][TT]  = gcov10  ;  gtmp[1][1]   = gcov11  ;
  gtmp[1][2]   = gcov12  ;  gtmp[1][3]   = gcov13  ;
  gtmp[2][TT]  = gcov20  ;  gtmp[2][1]   = gcov21  ;
  gtmp[2][2]   = gcov22  ;  gtmp[2][3]   = gcov23  ;
  gtmp[3][TT]  = gcov30  ;  gtmp[3][1]   = gcov31  ;
  gtmp[3][2]   = gcov32  ;  gtmp[3][3]   = gcov33  ;

  // Now transform:
  for( jp = 0 ; jp < NDIM; jp++ ) { 
    for( kp = 0 ; kp < NDIM; kp++ ) { 
      gcov[jp][kp] = 0.;
      DLOOP {
	gcov[jp][kp] +=  gtmp[j][k] * dx_dxp[j][jp] * dx_dxp[k][kp];
      }
    }
  }
      
}

void gcon_func(FTYPE *X, FTYPE gcon[][NDIM])
{
  int j, k, jp, kp;
  FTYPE sth, cth, s2, rho2;
  FTYPE r, th;
  FTYPE dxp_dx[NDIM][NDIM];
  static FTYPE gtmp[NDIM][NDIM];

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth=sin(th);

  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  dxp_dx_calc( X, r, th, dxp_dx );

  // Set the metric first that in KS coordinates:
  gtmp[TT][TT] = gcon00  ;  gtmp[TT][1]  = gcon01  ;
  gtmp[TT][2]  = gcon02  ;  gtmp[TT][3]  = gcon03  ;
  gtmp[1][TT]  = gcon10  ;  gtmp[1][1]   = gcon11  ;
  gtmp[1][2]   = gcon12  ;  gtmp[1][3]   = gcon13  ;
  gtmp[2][TT]  = gcon20  ;  gtmp[2][1]   = gcon21  ;
  gtmp[2][2]   = gcon22  ;  gtmp[2][3]   = gcon23  ;
  gtmp[3][TT]  = gcon30  ;  gtmp[3][1]   = gcon31  ;
  gtmp[3][2]   = gcon32  ;  gtmp[3][3]   = gcon33  ;

  // Now transform:
  for( jp = 0 ; jp < NDIM; jp++ ) { 
    for( kp = 0 ; kp < NDIM; kp++ ) { 
      gcon[jp][kp] = 0.;
      DLOOP {
	gcon[jp][kp] +=  gtmp[j][k] * dxp_dx[jp][j] * dxp_dx[kp][k];
      }
    }
  }
      
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

  if(POSDEFMETRIC){
    if (th < 0)
      th *= -1.;
    if (th > M_PI)
      th = 2. * M_PI - th;
  }

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
#if(COORDSINGFIX)
  if (fabs(sth) < SINGSMALL) {
    if(sth>=0) sth=SINGSMALL;
    if(sth<0) sth=-SINGSMALL;
  }
#endif
  
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
void gcon_func_old(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
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

// Avery mentions that long double trig. functions only return double precision answer.  see ~/research/utils/triglongdouble.c

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

FTYPE csc(FTYPE arg)
{
  return(1.0/sin(arg));
}



/*********************************************************************************************
   Scott's s MKS connection that can be used for any transformation between r,th <-> X1,X2 :
*********************************************************************************************/
#define M (1.)
void mks_conn_func(FTYPE *X, struct of_geom *geom, FTYPE conn[][NDIM][NDIM] )
{
  int i, j, k, l;
  FTYPE r,th,sigma,dx_dxp[NDIM][NDIM],dx_dxp_dxp[NDIM][NDIM][NDIM];

  double t1,   t10,   t102,   t1024,   t1035,   t1037,   t104,   t11,   t114,   t116,   t119;
  double  t12,   t121,   t123,   t126,   t129,   t130,   t132,   t14,   t148,   t149,   t15,   t152;
  double t154,   t156,   t157,   t158,   t159,   t161,   t169,   t17,   t171,   t172,   t175,   t177;
  double t185,   t2,   t203,   t204,   t208,   t209,   t21,   t210,   t212,   t214,   t22,   t221;
  double   t222,   t224,   t227,   t23,   t236,   t24,   t240,   t241,   t242,   t243,   t245,   t246;
  double    t247,   t248,   t25,   t250,   t251,   t258,   t26,   t260,   t261,   t263,   t264,   t271;
  double   t273,   t275,   t276,   t278,   t28,   t280,   t281,   t283,   t284,   t285,   t286,   t288;
  double    t289,   t29,   t297,   t299,   t3,   t30,   t300,   t303,   t305,   t306,   t308,   t309;
  double    t31,   t310,   t313,   t314,   t320,   t325,   t327,   t328,   t329,   t330,   t333,   t336;
  double    t338,   t34,   t340,   t342,   t344,   t346,   t35,   t358,   t361,   t363,   t366,   t367;
  double    t368,   t370,   t372,   t375,   t38,   t380,   t381,   t384,   t385,   t387,   t39,   t399;
  double    t4,   t40,   t400,   t402,   t404,   t405,   t406,   t408,   t409,   t41,   t411,   t412;
  double    t418,   t42,   t421,   t425,   t428,   t431,   t432,   t433,   t434,   t437,   t440,   t442;
  double    t448,   t451,   t453,   t454,   t459,   t462,   t467,   t469,   t480,   t481,   t486,   t487;
  double    t488,   t491,   t492,   t498,   t501,   t504,   t507,   t508,   t510,   t512,   t52,   t521;
  double    t528,   t53,   t530,   t553,   t556,   t56,   t57,   t588,   t60,   t607,   t627,   t628;
  double    t63,   t630,   t631,   t632,   t634,   t636,   t637,   t64,   t651,   t652,   t654,   t656;
  double    t657,   t659,   t661,   t662,   t670,   t673,   t675,   t677,   t686,   t689,   t7,   t712;
  double    t74,   t748,   t75,   t78,   t793,   t794,   t795,   t799,   t8,   t800,   t801,   t803;
  double    t806,  t807,   t813,   t816,   t822,   t83,   t831,   t84,   t845,   t86,   t89,   t891; 
  double    t90,   t91,   t916,   t917,   t920,   t924,   t928,   t940,   t946,   t968,   t97, t970,t991;

  // get bl coordinates
  bl_coord(X,&r,&th);
  // the connection

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dx_dxp_calc(X,r,th,dx_dxp);
  dx_dxp_dxp_calc(X,r,th,dx_dxp_dxp);


  //  t1 = rf(X1,X2);
  t1 = r;
  //  t2 = thf(X1,X2);
  t2 = th;
  t3 = cos(t2);
  t4 = a*t3;
  t7 = (-t1+t4)*(t1+t4);
  t8 = M*M;
  t10 = t1*t1;
  t11 = a*a;
  t12 = t3*t3;
  t14 = t10+t11*t12;
  t15 = t14*t14;
  t17 = 1/t15/t14;
  conn[0][0][0] = -2.0*t7*t8*t1*t17;
  t21 = t10*t10;
  //  t22 = diff(rf(X1,X2),X1);
  t22 = dx_dxp[RR][1];
  t23 = t21*t22;
  t24 = t10*t1;
  t25 = M*t24;
  t26 = t25*t22;
  t28 = t24*t3;
  t29 = sin(t2);
  //  t30 = diff(thf(X1,X2),X1);
  t30 = dx_dxp[TH][1];
  t31 = t29*t30;
  t34 = M*t1;
  t35 = t22*t12;
  t38 = t12*t12;
  t39 = t38*t22;
  t40 = t29*t1;
  t41 = t12*t3;
  t42 = t41*t30;
  conn[0][0][1] = -M*(-t23-2.0*t26+(2.0*t28*t31+2.0*t34*t35+(t39+2.0*t40*t42)*t11)*t11)*t17;
  //  t52 = diff(rf(X1,X2),X2);
  t52 = dx_dxp[RR][2];
  t53 = t21*t52;
  //  t56 = diff(thf(X1,X2),X2);
  t56 = dx_dxp[TH][2];
  t57 = t29*t56;
  t60 = t52*t12;
  t63 = t38*t52;
  t64 = t41*t1;
  conn[0][0][2] = -M*(-t53-2.0*t25*t52+(2.0*t28*t57+2.0*t34*t60+(t63+2.0*t64*t57)*t11)*t11)*t17;
  t74 = -1.0+t3;
  t75 = 1.0+t3;
  t78 = t7*t74*t75;
  conn[0][0][3] = -2.0*t78*a*t1*t8*t17;
  conn[0][1][0] = conn[0][0][1];
  t83 = t30*t30;
  t84 = t21*t10;
  t86 = t22*t22;
  t89 = t24*t22;
  t90 = t3*t29;
  t91 = t90*t30;
  t97 = t86*t12;
  t102 = t22*t29;
  t104 = t64*t102;
  conn[0][1][1] = -2.0*(t83*t84-t21*t86-t25*t86+(2.0*t89*t91+2.0*t83*t21*t12+
						 t34*t97+(t83*t10*t38+t38*t86+2.0*t104*t30)*t11)*t11)*M*t17;
  t114 = t22*t52;
  t116 = t30*t56;
  t119 = t24*t52;
  t121 = t114*t12;
  t123 = t21*t12;
  t126 = t90*t56;
  t129 = t52*t29;
  t130 = t129*t30;
  t132 = t10*t38;
  conn[0][1][2] = -2.0*(-t25*t114+t116*t84-t23*t52+(t119*t91+t34*t121+2.0*
						    t116*t123+t89*t126
						    +(t39*t52+t64*t130+t116*t132+t104*t56)*t11)*t11)*M*t17;
  t148 = 2.0*t28*t30;
  t149 = t102*t12;
  t152 = t41*t24;
  t154 = 2.0*t152*t30;
  t156 = 2.0*t64*t30;
  t157 = t39*t29;
  t158 = t30*t1;
  t159 = t38*t3;
  t161 = 2.0*t158*t159;
  t169 = a*t29*t17;
  conn[0][1][3] = -(2.0*t25*t102+t23*t29+(-t148-2.0*t34*t149+t154+(-t156-t157+t161)*t11)*t11)*M*t169;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = conn[0][1][2];
  t171 = t52*t52;
  t172 = t171*M;
  t175 = t56*t56;
  t177 = t1*t12;
  t185 = t52*t41;
  conn[0][2][2] = -2.0*(-t172*t24-t171*t21+t175*t84+(t172*t177+2.0*t119*t126+
						     2.0*t175*t21*t12
						     +(t171*t38+2.0*t185*t40*t56+t175*t10*t38)*t11)*t11)*M*t17;
  t203 = 2.0*t152*t56;
  t204 = t129*t12;
  t208 = 2.0*t28*t56;
  t209 = t63*t29;
  t210 = t56*t1;
  t212 = 2.0*t210*t159;
  t214 = 2.0*t64*t56;
  conn[0][2][3] = (-2.0*t25*t129-t53*t29+(-t203+2.0*t34*t204+t208+(t209-t212+t214)*t11)*t11)*M*t169;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t221 = t21*t1;
  t222 = t24*t12;
  t224 = t10*t12;
  t227 = t1*t38;
  t236 = (-t221+(-2.0*t222+(-t224+t10)*M+(-t227+(t38-t12)*M)*t11)*t11)*t74*t75;
  conn[0][3][3] = -2.0*t236*t34*t17;
  t240 = t56*M;
  t241 = t240*t24;
  t242 = 2.0*t241;
  t243 = t56*t21;
  t245 = 2.0*t240*t177;
  t246 = t56*t10;
  t247 = t246*t12;
  t248 = t1*t3;
  t250 = 2.0*t248*t129;
  t251 = t56*t12;
  t258 = t30*t52;
  t260 = 1/(-t22*t56+t258);
  t261 = t260*t17;
  conn[1][0][0] = -(-t242+t243+(t245-t247+t250+t246-t251*t11)*t11)*M*t261;
  t263 = M*t22;
  t264 = t56*t38;
  t271 = (-t242+(t245-t247+t250+t246+(-t251+t264)*t11)*t11)*t260*t17;
  conn[1][0][1] = -t263*t271;
  t273 = M*t52;
  conn[1][0][2] = -t273*t271;
  t275 = t24*t29;
  t276 = t240*t275;
  t278 = t119*t3;
  t280 = t243*t29;
  t281 = t40*t12;
  t283 = 2.0*t240*t281;
  t284 = t29*t12;
  t285 = t246*t284;
  t286 = t52*t3;
  t288 = 2.0*t286*t1;
  t289 = t57*t10;
  t297 = a*t29*t260*t17;
  conn[1][0][3] = (-2.0*t276+2.0*t278+t280+(t283-t285+t288+t289-t56*t11*t284)*t11)*M*t297;
  conn[1][1][0] = conn[1][0][1];
  //  t299 = diff(diff(thf(X1,X2),X1),X1);
  t299 = dx_dxp_dxp[TH][1][1];
  t300 = t52*t299;
  t303 = t258*t221*t22;
  t305 = t56*t84;
  //  t306 = diff(diff(rf(X1,X2),X1),X1);
  t306 = dx_dxp_dxp[RR][1][1] ;
  t308 = t83*t56;
  t309 = t21*t24;
  t310 = t308*t309;
  t313 = 2.0*t308*t84;
  t314 = t56*t24;
  t320 = t30*t22;
  t325 = t258*t89*t12;
  t327 = t308*t221;
  t328 = t52*t83;
  t329 = t90*t21;
  t330 = t328*t329;
  t333 = t306*t12;
  t336 = t221*t12;
  t338 = 2.0*t308*t336;
  t340 = 4.0*t308*t123;
  t342 = t248*t29;
  t344 = 2.0*t52*t86*t342;
  t346 = t56*t86;
  t358 = t306*t38;
  t361 = t1*t22;
  t363 = t258*t361*t38;
  t366 = 2.0*t308*t222;
  t367 = t24*t38;
  t368 = t308*t367;
  t370 = t41*t29*t10;
  t372 = 2.0*t328*t370;
  t375 = 2.0*t308*t132;
  t380 = t159*t29;
  t381 = t380*t56;
  t384 = t308*t227;
  t385 = t38*t12;
  t387 = t328*t380;
  conn[1][1][1] = -(-t300*t84-2.0*t303+t305*t306-t310+(-t243*t86+t313-2.0*t314*t86*M)*M
		    +(-2.0*t320*t3*t280-4.0*t325-t327+t330-3.0*t300*t123+3.0*t243*t333
		      -t338+(t340+t344-t246*t97+t346*t10+2.0*t210*t97*M)*M
		      +(-4.0*t320*t41*t289-3.0*t300*t132+3.0*t246*t358-2.0*t363-t366-t368+t372
			+(-t346*t12+t375+2.0*t346*t38)*M
			+(-2.0*t320*t381-t384-t300*t385+t387+t56*t306*t385)*t11)*t11)*t11)*t260*t17;
  //  t399 = diff(diff(thf(X1,X2),X1),X2);
  t399 = dx_dxp_dxp[TH][1][2]; 
  t400 = t52*t399;
  //  t402 = diff(diff(rf(X1,X2),X1),X2);
  t402 = dx_dxp_dxp[RR][1][2]; 
  t404 = t30*t175;
  t405 = t404*t309;
  t406 = t171*t30;
  t408 = t56*t221;
  t409 = t114*t408;
  t411 = 2.0*t404*t84;
  t412 = t52*t56;
  t418 = t402*t12;
  t421 = t404*t221;
  t425 = t114*t314*t12;
  t428 = 2.0*t404*t336;
  t431 = t22*t175;
  t432 = t431*t329;
  t433 = t10*t22;
  t434 = t433*t12;
  t437 = 4.0*t404*t123;
  t440 = 2.0*t171*t22*t342;
  t442 = t210*M;
  t448 = 2.0*t370*t431;
  t451 = t114*t210*t38;
  t453 = 2.0*t404*t222;
  t454 = t402*t38;
  t459 = t404*t367;
  t462 = 2.0*t404*t132;
  t467 = t404*t227;
  t469 = t431*t380;
  conn[1][1][2] = (t400*t84-t305*t402+t405+t406*t221+t409+(-t411+t412*t23+2.0*t114*t241)*M
		   +(-3.0*t243*t418+t421+3.0*t400*t123+2.0*t425+t428+2.0*t406*t222+
		     t432+(t412*t434-t437-t440-t412*t433-2.0*t121*t442)*M
		     +(t448+t406*t227+t451+t453-3.0*t246*t454+3.0*t400*t132+t459
		       +(t114*t251-t462-2.0*t114*t264)*M+(t467+t400*t385+t469
							  -t56*t402*t385)*t11)*t11)*t11)*t260*t17;
  t480 = t286*t21;
  t481 = t57*t221;
  t486 = 2.0*t185*t10;
  t487 = t57*t222;
  t488 = 2.0*t487;
  t491 = t57*t227;
  t492 = t52*t159;
  t498 = (-t491+t492+(-t57*t12+t57*t38)*M)*t11;
  t501 = t480-t481+(2.0*t278-2.0*t276)*M+(t486-t488+(-t285+t289+t288+t283)*M+t498)*t11;
  conn[1][1][3] = t501*t22*t297;
  conn[1][2][0] = conn[1][0][2];
  conn[1][2][1] = conn[1][1][2];
  t504 = t171*t56;
  //  t507 = diff(diff(thf(X1,X2),X2),X2);
  t507 = dx_dxp_dxp[TH][2][2];
  t508 = t52*t507;
  //  t510 = diff(diff(rf(X1,X2),X2),X2);
  t510 = dx_dxp_dxp[RR][2][2];
  t512 = t175*t56;
  t521 = t512*t221;
  t528 = t52*t175;
  t530 = t510*t12;
  t553 = t510*t38;
  t556 = t512*t24;
  conn[1][2][2] = -(-2.0*t504*t221-t508*t84+t305*t510-t512*t309
		    +(2.0*t512*t84-t504*t21-2.0*t504*t25)*M
		    +(-2.0*t521*t12-3.0*t508*t123-4.0*t504*t222-t521-t528*t329+3.0*t243*t530
		      +(-t504*t224+t504*t10+2.0*t342*t171*t52+4.0*t512*t21*t12+2.0*t12*t171*t442)*M
		      +(-3.0*t508*t132-2.0*t504*t227-2.0*t528*t370+3.0*t246*t553-2.0*t556*t12-t556*t38
			+(2.0*t512*t10*t38-t504*t12+2.0*t504*t38)*M
			+(-t528*t380-t512*t1*t38-t508*t385+t56*t510*t385)*t11)*t11)*t11)*t260*t17;
  conn[1][2][3] = t501*t52*t297;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = conn[1][2][3];
  t588 = t84*t29;
  t607 = t29*t38;
  conn[1][3][3] = -(-t56*t309*t29+t286*t84+2.0*t240*t588
		    +(-t481-2.0*t408*t284+t480+2.0*t185*t21+(-4.0*t185*t24+t280+3.0*t243*t284+4.0*t278
							     +(-2.0*t314*t29+2.0*t487)*M)*M
		      +(t492*t10+t486-t314*t607-t488+(-2.0*t492*t1+t289+t288-2.0*t285+3.0*t246*t607
						      +(-2.0*t491+2.0*t210*t284)*M)*M+t498)*t11)*t11)*t29*t261;
  t627 = t30*t21;
  t628 = t30*M;
  t630 = 2.0*t628*t24;
  t631 = t30*t10;
  t632 = t631*t12;
  t634 = 2.0*t628*t177;
  t636 = 2.0*t361*t90;
  t637 = t30*t12;
  conn[2][0][0] = -(-t627+t630+(t632-t634-t631-t636+t637*t11)*t11)*M*t261;
  t651 = (-t630+(t634+t636+t631-t632+(-t637+t30*t38)*t11)*t11)*t260*t17;
  conn[2][0][1] = t263*t651;
  conn[2][0][2] = t273*t651;
  t652 = t628*t275;
  t654 = t89*t3;
  t656 = t627*t29;
  t657 = t631*t284;
  t659 = 2.0*t628*t281;
  t661 = 2.0*t361*t3;
  t662 = t31*t10;
  conn[2][0][3] = (2.0*t652-2.0*t654-t656+(t657-t659-t661-t662+t30*t11*t284)*t11)*M*t297;
  conn[2][1][0] = conn[2][0][1];
  t670 = t30*t86;
  t673 = t30*t84;
  t675 = t22*t299;
  t677 = t83*t30;
  t686 = t677*t221;
  t689 = t83*t22;
  t712 = t677*t24;
  conn[2][1][1] = -(2.0*t670*t221-t673*t306+t675*t84+t677*t309+(-2.0*t677*t84+t627*t86+2.0*t670*t25)*M
		    +(2.0*t686*t12+t689*t329-3.0*t627*t333+t686+4.0*t670*t222+3.0*t675*t123
		      +(-2.0*t86*t22*t1*t90+t631*t97-t670*t10-4.0*t677*t21*t12-2.0*t637*t86*t34)*M
		      +(2.0*t712*t12+t712*t38+2.0*t670*t227+3.0*t675*t132-3.0*t631*t358+2.0*t689*t370
			+(t670*t12-2.0*t677*t10*t38-2.0*t670*t38)*M
			+(t675*t385-t30*t306*t385+t689*t380+t677*t1*t38)*t11)*t11)*t11)*t260*t17;
  t748 = t22*t399;
  conn[2][1][2] = -(t303+t310-t673*t402+t748*t84+t346*t221
		    +(-t313+t258*t23+2.0*t258*t26)*M
		    +(t327+t330+2.0*t325-3.0*t627*t418+2.0*t346*t222+3.0*t748*t123+t338
		      +(-t340-t258*t433-t344+t258*t434-2.0*t60*t30*t361*M)*M
		      +(t346*t227+t372+t363+t366+t368-3.0*t631*t454+3.0*t748*t132
			+(t258*t35-t375-2.0*t258*t39)*M+(t387+t748*t385+t384
							 -t30*t402*t385)*t11)*t11)*t11)*t260*t17;
  t793 = t31*t221;
  t794 = t3*t22;
  t795 = t794*t21;
  t799 = t31*t222;
  t800 = 2.0*t799;
  t801 = t22*t41;
  t803 = 2.0*t801*t10;
  t806 = t31*t227;
  t807 = t22*t159;
  t813 = (-t806+t807+(t31*t38-t31*t12)*M)*t11;
  t816 = -t793+t795+(2.0*t654-2.0*t652)*M+(-t800+t803+(t661-t657+t662+t659)*M+t813)*t11;
  conn[2][1][3] = -t816*t22*t297;
  conn[2][2][0] = conn[2][0][2];
  conn[2][2][1] = conn[2][1][2];
  t822 = t22*t507;
  t831 = t3*t56;
  t845 = t41*t56;
  conn[2][2][2] = -(-t673*t510+2.0*t409+t405+t822*t84+(t406*t21-t411+2.0*t406*t25)*M
		    +(-3.0*t627*t530+t421-t432+2.0*t130*t831*t21+t428+3.0*t822*t123+4.0*t425
		      +(t406*t224-t406*t10-t440-t437-2.0*t406*t177*M)*M
		      +(4.0*t130*t845*t10+2.0*t451+3.0*t822*t132+t453+t459-t448-3.0*t631*t553
			+(t406*t12-t462-2.0*t406*t38)*M+(2.0*t258*t381+t822*t385-t30*t510*t385
							 +t467-t469)*t11)*t11)*t11)*t260*t17;
  conn[2][2][3] = -t816*t52*t297;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = conn[2][2][3];
  t891 = t30*t24;
  conn[2][3][3] = (t794*t84+2.0*t628*t588-t30*t309*t29
		   +(t795-t793-2.0*t30*t221*t284+2.0*t801*t21
		     +(3.0*t627*t284+t656-4.0*t801*t24+4.0*t654+(-2.0*t891*t29+2.0*t799)*M)*M
		     +(t807*t10+t803-t891*t607-t800+(3.0*t631*t607-2.0*t807*t1+t661-2.0*t657+t662
						     +(2.0*t158*t284-2.0*t806)*M)*M+t813)*t11)*t11)*t29*t261;
  t917 = a*M;
  conn[3][0][0] = -t7*t917*t17;
  t920 = t102*t10;
  t924 = 1/t29;
  t916 = t924*t17;
  conn[3][0][1] = -t917*(-t920+t148+(t156+t149)*t11)*t916;
  t928 = t129*t10;
  conn[3][0][2] = -t917*(t208-t928+(t214+t204)*t11)*t916;
  conn[3][0][3] = -t78*M*t11*t17;
  conn[3][1][0] = conn[3][0][1];
  t940 = t83*t29;
  t946 = t86*t29;
  t968 = t916;
  conn[3][1][1] = -(t940*t221+2.0*t794*t627+(4.0*t794*t891-t946*t10)*M
		    +(4.0*t801*t631+2.0*t940*t222+(4.0*t361*t42+t946*t12)*M+(t940*t227+2.0*t807*t30)*t11)
		    *t11)*a*t968;
  t970 = t29*t221;
  t991 = t52*t1;
  conn[3][1][2] = -(t116*t970+t286*t627+t794*t243
		    +(2.0*t286*t891-t102*t52*t10+2.0*t794*t314)*M
		    +(2.0*t185*t631+2.0*t801*t246+2.0*t116*t275*t12+(2.0*t361*t845+2.0*t991*t42+t102*t60)*M
		      +(t116*t40*t38+t492*t30+t807*t56)*t11)*t11)*a*t968;
  t1024 = t38*t41;
  conn[3][1][3] = (t3*t30*t84+t970*t22+(3.0*t42*t21+2.0*t275*t35
					+(t102*t224+t148-t154-t920)*M
					+(t40*t39+3.0*t159*t30*t10+(t156-t161+t149-t157)*M
					  +t1024*t30*t11)*t11)*t11)*t924*t17;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  t1035 = t175*t29;
  t1037 = t171*t29;
  conn[3][2][2] = -(2.0*t286*t243+t1035*t221+(-t1037*t10+4.0*t286*t314)*M
		    +(4.0*t185*t246+2.0*t1035*t222+(t1037*t12+4.0*t991*t845)*M
		      +(t1035*t227+2.0*t492*t56)*t11)*t11)*a*t968;
  conn[3][2][3] = -(-t970*t52-t831*t84+(-2.0*t275*t60-3.0*t845*t21
					+(t203-t208-t129*t224+t928)*M
					+(-t40*t63-3.0*t159*t56*t10+(t209+t212-t204-t214)*M
					  -t1024*t56*t11)*t11)*t11)*t924*t17;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = -t236*a*t17;

  return;

}

#undef M




// jon's MKS connection (and conn2)
void mks_conn_func_old(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE r,th,sigma,dxdxp[NDIM];
  FTYPE cot(FTYPE arg);

  // get bl coordinates
  bl_coord(X,&r,&th);
  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, 
  // but shouldn't/doesn't seem to matter much follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
  }
  else{
    if(th>M_PI) { th=M_PI-th; }
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

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
// jon's MKS source in mid-simplified form
void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom, int ii, int jj,struct of_state *q,
		FTYPE *dU)
{
  int i=0, j=0, k=0, l=0;
  FTYPE r,th,X[NDIM],sigma,dxdxp[NDIM];
  FTYPE cot(FTYPE arg),csc(FTYPE arg);
  FTYPE b[NDIM],u[NDIM],bsq,en,rho;

  bsq = dot(q->bcon, q->bcov);
  u[TT]=q->ucon[TT];
  u[RR]=q->ucon[RR];
  u[TH]=q->ucon[TH];
  u[PH]=q->ucon[PH];

  b[TT]=q->bcon[TT];
  b[RR]=q->bcon[RR];
  b[TH]=q->bcon[TH];
  b[PH]=q->bcon[PH];

  rho=ph[RHO];
  en=ph[UU];

  coord(ptrgeom->i, ptrgeom->j, ptrgeom->p, X);
  // get bl coordinates
  bl_coord(X,&r,&th);

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

dU[U1]+=0.5*dxdxp[1]*pow(sigma,-5.)*
    (r*pow(sigma,4.)*(bsq - 2.*en + 2.*en*gam - 
         2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
         (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
          (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.)) + 
      a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) - 2.*a*
          ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)))*pow(sin(th),2.) + 
      dxdxp[1]*a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
          pow(dxdxp[1],-1.) + 
         sigma*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           - 1.*a*pow(dxdxp[1],-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
            2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       pow(sin(th),2.) - 1.*sigma*
       (4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*pow(sigma,-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         2.*r*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)) - 
         1.*a*(-1.*b[TT]*b[PH] + (bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
      (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
         1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
            1.*dxdxp[1]*b[RR]*b[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
            dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         0.5*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
            2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
         )*pow(sin(th),2.)*(r*pow(sigma,2.) + 
         pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
         ) + 2.*a*sigma*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
            2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
          - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
         2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
            (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.)) + 
      a*sigma*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
         2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
         cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
       pow(sin(th),2.)*(2.*r*(-1.*b[TT]*b[PH] + 
            (bsq + en*gam + rho)*u[TT]*u[PH]) + 
         dxdxp[1]*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*pow(csc(th),2.) - 
            1.*pow(b[PH],2.) + (bsq + en*gam + rho)*pow(u[PH],2.))*
          pow(sin(th),2.)) + sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
         2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
         1.*(bsq + en*gam + rho)*
          (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
         4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
         4.*r*(bsq + en*gam + rho)*u[TT]*
          (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))) + 
      2.*sigma*(2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
            (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
         dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
             ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(u[RR],2.)) - 
         1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*
       (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
         dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
            pow(a,2.)*pow(cos(th),4.) + 
            (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
      4.*dxdxp[1]*sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
         2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
         0.5*(bsq + en*gam + rho)*u[RR]*
          (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
            4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))) - 
      1.*dxdxp[2]*a*cos(th)*pow(sigma,2.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
         1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) + 2.*dxdxp[1]*a*
          (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*sin(th) - 
      2.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,3.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
         dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th) + 
      2.*dxdxp[2]*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
       pow(a,2.)*pow(sigma,3.)*sin(2.*th) + 
      2.*dxdxp[1]*dxdxp[2]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(a,2.)*pow(sigma,3.)*
       sin(2.*th) - 2.*dxdxp[2]*r*pow(a,2.)*pow(sigma,2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
         1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
         2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
         )*sin(2.*th) - 2.*dxdxp[2]*a*
       (b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sigma,3.)*sin(th)*
       (r*(2. + r)*sigma*cos(th) + sigma*pow(a,2.)*pow(cos(th),3.) + 
         r*pow(a,2.)*sin(th)*sin(2.*th)));
dU[U2]+=0.5*(-2.*dxdxp[1]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(dxdxp[2],2.) - 
      1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-2.)*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
         1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) + 2.*dxdxp[1]*a*
          (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*pow(sin(th),2.) - 
      4.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
         1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
         2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
         ) - 2.*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*pow(sigma,-2.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
         dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
      4.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
            2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
          - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
         2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
            (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.))*pow(sin(th),3.)
        - 2.*dxdxp[2]*a*r*cos(th)*pow(sigma,-4.)*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
          (r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.)) - 2.*a*
          ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
         4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
             (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
            (bsq + en*gam + rho)*pow(u[TT],2.)))*sin(th) - 
      1.*dxdxp[1]*dxdxp[2]*a*cos(th)*pow(sigma,-4.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
            dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
          pow(dxdxp[1],-1.) + 
         sigma*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           - 1.*a*pow(dxdxp[1],-1.)*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
            2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       sin(th) - 2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-2.)*
       (2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
            (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
         dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
          (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
             ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
            (bsq + en*gam + rho)*pow(u[RR],2.)) - 
         1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
            0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
           *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th)\
       + (bsq - 2.*en + 2.*en*gam - 2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
         (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
          (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.))*
       (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
         1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) - 
      1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
         2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
         1.*(bsq + en*gam + rho)*
          (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
         4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
         4.*r*(bsq + en*gam + rho)*u[TT]*
          (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))*sin(2.*th) - 
      2.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
         2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
         0.5*(bsq + en*gam + rho)*u[RR]*
          (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
            4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))))*sin(2.*th) + 
      dxdxp[2]*pow(sigma,-2.)*
       (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
         1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
            1.*dxdxp[1]*b[RR]*b[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
            dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
             (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
          (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
         0.5*(r*(2. + 3.*r)*pow(a,2.) + 
            cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
            2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
            2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
         )*pow(sin(th),2.)*(cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)))\
    ;

}
