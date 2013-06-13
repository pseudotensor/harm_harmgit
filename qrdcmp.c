#include "decs.h"
#define NRANSI
#include "NRUTIL.H"

void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing)
{
  int i,j,k;
  FTYPE scale,sigma,sum,tau;

  *sing=0;
  for (k=1;k<n;k++) {
    scale=0.0;
    for (i=k;i<=n;i++) scale=FMAX(scale,fabs(a[i][k]));
    if (scale == 0.0) {
      *sing=1;
      c[k]=d[k]=0.0;
    } else {
      for (i=k;i<=n;i++) a[i][k] /= scale;
      for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
      sigma=SIGN(sqrt(sum),a[k][k]);
      a[k][k] += sigma;
      c[k]=sigma*a[k][k];
      d[k] = -scale*sigma;
      for (j=k+1;j<=n;j++) {
        for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
        tau=sum/c[k];
        for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
      }
    }
  }
  d[n]=a[n][n];
  if (d[n] == 0.0) *sing=1;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
