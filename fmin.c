#include "decs.h"
#define NRANSI
#include "NRUTIL.H"


FTYPE nrfmin(FTYPE parms[], FTYPE x[])
{
  int i;
  FTYPE sum;

  (*nrfuncv)(nn,parms,x,fvec);
  for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
  return 0.5*sum;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
