#include "decs.h"

void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[])
{
  int i,j;
  FTYPE sum;

  b[n] /= d[n];
  for (i=n-1;i>=1;i--) {
    for (sum=0.0,j=i+1;j<=n;j++) sum += a[i][j]*b[j];
    b[i]=(b[i]-sum)/d[i];
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
