#define NRANSI
#include "decs.h"

void bcuint(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE x1l,
            FTYPE x1u, FTYPE x2l, FTYPE x2u, FTYPE x1, FTYPE x2, FTYPE *ansy,
            FTYPE *ansy1, FTYPE *ansy2)
{
  void bcucof(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE d1,
              FTYPE d2, FTYPE **c);
  int i;
  FTYPE t,u,d1,d2,**c;

  c=fmatrix(1,4,1,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  *ansy=(*ansy2)=(*ansy1)=0.0;
  for (i=4;i>=1;i--) {
    *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
  free_fmatrix(c,1,4,1,4);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
