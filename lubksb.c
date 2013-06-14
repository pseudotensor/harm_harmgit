#include "decs.h"

void lubksb(FTYPE **a, int n, int *indx, FTYPE b[])
{
  int i, ii = 0, ip, j;
  FTYPE sum;

  for (i = 1; i <= n; i++) {
#if(PRODUCTION==0)
    // if ludcmp() fails, indx could be 0
    if(indx[i]<1 || indx[i]>n){
      stderrfprintf("Major failure in lubksb: i=%d indx=%d\n",i,indx[i]);
    }
#endif
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n; i >= 1; i--) {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}
