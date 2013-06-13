#include "decs.h"
#define NRANSI
#include "NRUTIL.H"
// smaller doesn't always mean better
//#define EPS (NUMSQRTEPSILON*10)
#define EPS (1E-4)

void fdjac(int n, FTYPE parms[], FTYPE x[], FTYPE fvec[], FTYPE **df,
           void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]))
{
  int i,j;
  FTYPE h,temp,*f;

  f=vector(1,n);
  for (j=1;j<=n;j++) {
    temp=x[j];
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    x[j]=temp+h;
    h=x[j]-temp;
    (*vecfunc)(n,parms,x,f);
    x[j]=temp;
    for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
  }
  free_vector(f,1,n);
}
#undef EPS
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
