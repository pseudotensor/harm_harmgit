
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
 * cfg 11-15-02 changed to streamlined form for
 *              velocity primitive variables
 *
 */

#include "decs.h"

int dudp_calc_3vel(int whicheos, int whichcons, FTYPE *EOSextra, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE **alpha)
{
  FTYPE w, bsq;
  FTYPE ducov_dv[NDIM][NDIM], ducon_dv[NDIM][NDIM];
  FTYPE dbcov_dv[NDIM][NDIM], dbcon_dv[NDIM][NDIM];
  FTYPE dbsq_dv[NDIM], dw_dv[NDIM], dptot_dv[NDIM];
  int i, j, k;

#if(REMOVERESTMASSFROMUU==2)
  dualfprintf(fail_file,"This code does not handle REMOVERESTMASSFROMUU==2: dudp_calc_3vel()\n");
  myexit(1);
  // CODE ALSO ONLY FOR IDEAL EOS and no entropy evolution
#endif

#if(WHICHEOS!=IDEALGAS)
  dualfprintf(fail_file,"This code does not handle non-ideal gases: dudp_calc_3vel()\n");
  myexit(1);
#endif

#if(DOENTROPY!=DONOENTROPY)
  dualfprintf(fail_file,"This code does not handle entropy evolution: dudp_calc_3vel()\n");
  myexit(1);
#endif


  bsq = dot(q->bcon, q->bcov);
  w = pr[RHO] + gam * pr[UU] + bsq;

  /* set ducon_dv */
  for (i = 0; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      ducon_dv[i][j] =
        q->ucon[TT] * (q->ucon[i] * q->ucov[j] + delta(i, j));

  /* set ducov_dv */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      ducov_dv[i][j] = 0.;
  for (i = 0; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
        ducov_dv[i][j] += ptrgeom->gcov[GIND(i,k)] * ducon_dv[k][j];

  /* set dbcon_dv */
  for (i = 1; i < NDIM; i++)
    dbcon_dv[TT][i] = pr[B1] * ducov_dv[1][i]
      + pr[B2] * ducov_dv[2][i]
      + pr[B3] * ducov_dv[3][i];
  for (i = 1; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      dbcon_dv[i][j] = -q->bcon[i] * ducon_dv[TT][j] / q->ucon[TT]
        + ducon_dv[i][j] * q->bcon[TT] / q->ucon[TT]
        + q->ucon[i] * dbcon_dv[TT][j] / q->ucon[TT];

  /* set dbcov_dv */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      dbcov_dv[i][j] = 0.;
  for (i = 0; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
        dbcov_dv[i][j] += ptrgeom->gcov[GIND(i,k)] * dbcon_dv[k][j];

  /* set dbsq_dv[i] */
  for (i = 1; i < NDIM; i++)
    dbsq_dv[i] = 0.;
  for (i = 1; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      dbsq_dv[i] += 2. * q->bcov[j] * dbcon_dv[j][i];

  /* set dw_dv */
  for (i = 1; i < NDIM; i++)
    dw_dv[i] = dbsq_dv[i];

  /* set dptot_dv */
  for (i = 1; i < NDIM; i++)
    dptot_dv[i] = 0.5 * dbsq_dv[i];

  /** set alpha matrix **/
  // alpha is dU^i/dp^j = alpha[i][j]
  alpha[RHO + 1][RHO + 1] = q->ucon[TT];
  alpha[RHO + 1][UU + 1] = 0.;
  alpha[RHO + 1][U1 + 1] = pr[RHO] * ducon_dv[TT][1];
  alpha[RHO + 1][U2 + 1] = pr[RHO] * ducon_dv[TT][2];
  alpha[RHO + 1][U3 + 1] = pr[RHO] * ducon_dv[TT][3];

  alpha[UU + 1][RHO + 1] = q->ucon[TT] * q->ucov[TT];
  alpha[U1 + 1][RHO + 1] = q->ucon[TT] * q->ucov[1];
  alpha[U2 + 1][RHO + 1] = q->ucon[TT] * q->ucov[2];
  alpha[U3 + 1][RHO + 1] = q->ucon[TT] * q->ucov[3];

  alpha[UU + 1][UU + 1] = gam * q->ucon[TT] * q->ucov[TT] + (gam - 1.);
  alpha[U1 + 1][UU + 1] = gam * q->ucon[TT] * q->ucov[1];
  alpha[U2 + 1][UU + 1] = gam * q->ucon[TT] * q->ucov[2];
  alpha[U3 + 1][UU + 1] = gam * q->ucon[TT] * q->ucov[3];

  for (i = 0; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      alpha[UU + i + 1][UU + 1 + j] =
        dw_dv[j] * q->ucon[TT] * q->ucov[i]
        + w * ducon_dv[TT][j] * q->ucov[i]
        + w * q->ucon[TT] * ducov_dv[i][j]
        + dptot_dv[j] * delta(TT, i)
        - dbcon_dv[TT][j] * q->bcov[i]
        - q->bcon[TT] * dbcov_dv[i][j];

  /* this bit of legacy code can be uncommented if the rest-mass flux
     is subtracted out from the energy flux, which may be numerically
     convenient */

  // this is handled at the top level in Utoprimgen() now.
  //  if(REMOVERESTMASSFROMUU){
  //    j = 2;
  //    for (k = 1; k <= 5; k++)
  //      alpha[j][k] += alpha[1][k];
  //  }

   
  // alpha is dU^i/dp^j = alpha[i][j]

  // if(WHICHEOM==WITHGDET){
  /* N.B.: all the conserved variables contain a factor of \sqrt{det(g_{\mu\nu})} */
  //   for(j=1;j<=5;j++)
  //     for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  // }
  // else if(WHICHEOM==WITHNOGDET){
  // only the U^t (dU^2/dp^j) and U^\phi (dU^5/dp^j) contain gdet
  //   for(j=1;j<=5;j++){
  //     if(j==1) for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  //     else if((j==2)&&(NOGDET0==0)) for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  //     else if((j==3)&&(NOGDET1==0)) for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  //     else if((j==4)&&(NOGDET2==0)) for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  //     else if((j==5)&&(NOGDET3==0)) for(k=1;k<=5;k++) alpha[j][k] *= ptrgeom->g ;
  //   }
  // }

  return (0);
}
