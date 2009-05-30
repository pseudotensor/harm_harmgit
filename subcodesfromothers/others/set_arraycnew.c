
#include "decs.h"

void set_arrays()
{
  int i, j, k, l;

  p = (FTYPE (*)[N2 + 4][NPR]) (&(a_p[2][2][0]));
  pflag = (int (*)[N2 + 4][NUMFLAGS]) (&(a_pflag[2][2][0]));
  dq1 = (FTYPE (*)[N2 + 4][NPR]) (&(a_dq1[2][2][0]));
  dq2 = (FTYPE (*)[N2 + 4][NPR]) (&(a_dq2[2][2][0]));
  F1 = (FTYPE (*)[N2 + 4][NPR]) (&(a_F1[2][2][0]));
  F2 = (FTYPE (*)[N2 + 4][NPR]) (&(a_F2[2][2][0]));
  F1CT = (FTYPE (*)[N2 + 4][NPR]) (&(a_F1CT[2][2][0]));
  F2CT = (FTYPE (*)[N2 + 4][NPR]) (&(a_F2CT[2][2][0]));
  ph = (FTYPE (*)[N2 + 4][NPR]) (&(a_ph[2][2][0]));
  prc = (FTYPE (*)[N2 + 4][NPR]) (&(a_prc[2][2][0]));

  // this faraday needed for current calculation
  cfaraday =  (FTYPE (*)[N2 + 4][4][3]) (&(a_cfaraday[2][2][0][0]));
  fcon =  (FTYPE (*)[N2 + 4][NUMFARADAY]) (&(a_fcon[2][2][0]));
  jcon = (FTYPE (*)[N2 + 4][NDIM]) (&(a_jcon[2][2][0]));
  

  // assume time average stuff gets zeroed in avg routine
#if(DOAVG)
  normalvarstavg =  (FTYPE (*)[N2 + 4][NUMNORMDUMP]) (&(a_normalvarstavg[2][2][0]));
  anormalvarstavg =  (FTYPE (*)[N2 + 4][NUMNORMDUMP]) (&(a_anormalvarstavg[2][2][0]));

  fcontavg =  (FTYPE (*)[N2 + 4][NUMFARADAY]) (&(a_fcontavg[2][2][0]));
  fcovtavg =  (FTYPE (*)[N2 + 4][NUMFARADAY]) (&(a_fcovtavg[2][2][0]));

  afcontavg =  (FTYPE (*)[N2 + 4][NUMFARADAY]) (&(a_afcontavg[2][2][0]));
  afcovtavg =  (FTYPE (*)[N2 + 4][NUMFARADAY]) (&(a_afcovtavg[2][2][0]));

  massfluxtavg =  (FTYPE (*)[N2 + 4][NDIM]) (&(a_massfluxtavg[2][2][0]));
  amassfluxtavg =  (FTYPE (*)[N2 + 4][NDIM]) (&(a_amassfluxtavg[2][2][0]));

  othertavg =  (FTYPE (*)[N2 + 4][NUMOTHER]) (&(a_othertavg[2][2][0]));
  aothertavg =  (FTYPE (*)[N2 + 4][NUMOTHER]) (&(a_aothertavg[2][2][0]));

  jcontavg = (FTYPE (*)[N2 + 4][NDIM]) (&(a_jcontavg[2][2][0]));
  jcovtavg = (FTYPE (*)[N2 + 4][NDIM]) (&(a_jcovtavg[2][2][0]));

  ajcontavg = (FTYPE (*)[N2 + 4][NDIM]) (&(a_ajcontavg[2][2][0]));
  ajcovtavg = (FTYPE (*)[N2 + 4][NDIM]) (&(a_ajcovtavg[2][2][0]));

  tudtavg = (FTYPE (*)[N2 + 4][NUMSTRESSTERMS]) (&(a_tudtavg[2][2][0]));
  atudtavg = (FTYPE (*)[N2 + 4][NUMSTRESSTERMS]) (&(a_atudtavg[2][2][0]));
#endif  
  

  /* everything must be initialized to zero */
  ZSLOOP(-2, N1 + 1, -2, N2 + 1) {
    for(k=0;k<NUMFLAGS;k++) pflag[i][j][k] = 0;
    PLOOP {
      p[i][j][k] = 0.;
      ph[i][j][k] = 0.;
      prc[i][j][k] = 0.;
      dq1[i][j][k] = 0.;
      dq2[i][j][k] = 0.;
      F1[i][j][k] = 0.;
      F2[i][j][k] = 0.;
      F1CT[i][j][k] = 0.;
      F2CT[i][j][k] = 0.;
    }
    for(k=0;k<4;k++) for(l=0;l<3;l++){
      cfaraday[i][j][k][l]=0.;
    }
    for(k=0;k<NUMFARADAY;k++){
      fcon[i][j][k]=0.;
    }
    for(k=0;k<NDIM;k++){
      jcon[i][j][k]=0.;
    }
  }

  ZLOOP stat[i][j] = 1;

  /* grid functions */
  conn = (FTYPE (*)[N2 + 4][NDIM][NDIM][NDIM])
      (&(a_conn[2][2][0][0][0]));
  gcon = (FTYPE (*)[N2 + 4][NPG][NDIM][NDIM])
      (&(a_gcon[2][2][0][0][0]));
  gcov = (FTYPE (*)[N2 + 4][NPG][NDIM][NDIM])
      (&(a_gcov[2][2][0][0][0]));
  gdet = (FTYPE (*)[N2 + 4][NPG])
      (&(a_gdet[2][2][0]));

}
