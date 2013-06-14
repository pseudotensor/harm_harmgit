#include "decs.h"

// JCM: removed (unsigned) inside malloc since (unsigned)=(unsigned int) and argument was long int -- makes non sense and will truncate if asking for more than what integer number is

void nrerror(char error_text[])
{
  dualfprintf(fail_file, "Numerical Recipes run-time error...\n");
  dualfprintf(fail_file, "%s\n", error_text);
  dualfprintf(fail_file, "...now exiting to system...\n");
  exit(1);
}



int *ivector(long nl, long nh)
{
  int *v;

  v = (int *) malloc( (nh - nl + 1) * sizeof(int));
  if (!v)
    nrerror("allocation failure in ivector()");
  return v - nl;
}



FTYPE *dvector(long nl, long nh)
{
  FTYPE *v;

  v = (FTYPE *) malloc( (nh - nl + 1) * sizeof(FTYPE));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl;
}
FTYPE *vector(long nl, long nh)
{
  FTYPE *v;

  v = (FTYPE *) malloc( (nh - nl + 1) * sizeof(FTYPE));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl;
}




FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  FTYPE **m;

  m = (FTYPE **) malloc( (nrh - nrl + 1) * sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] =
      (FTYPE *) malloc( (nch - ncl + 1) * sizeof(FTYPE));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

FTYPE **matrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  FTYPE **m;

  m = (FTYPE **) malloc( (nrh - nrl + 1) * sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] =
      (FTYPE *) malloc( (nch - ncl + 1) * sizeof(FTYPE));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  int **m;

  m = (int **) malloc( (nrh - nrl + 1) * sizeof(int *));
  if (!m)
    nrerror("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (int *) malloc( (nch - ncl + 1) * sizeof(int));
    if (!m[i])
      nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}



FTYPE **submatrix(FTYPE **a, long oldrl, long oldrh, long oldcl,
                  long oldch, long newrl, long newcl)
{
  long i, j;
  FTYPE **m;

  m = (FTYPE **) malloc( (oldrh - oldrl + 1) *
                         sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure in submatrix()");
  m -= newrl;

  for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = a[i] + oldcl - newcl;

  return m;
}

FTYPE **dsubmatrix(FTYPE **a, long oldrl, long oldrh, long oldcl,
                   long oldch, long newrl, long newcl)
{
  long i, j;
  FTYPE **m;

  m = (FTYPE **) malloc( (oldrh - oldrl + 1) *
                         sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure in submatrix()");
  m -= newrl;

  for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = a[i] + oldcl - newcl;

  return m;
}



void free_vector(FTYPE *v, long nl, long nh)
{
  free((char *) (v + nl));
  //  free(&v[nl]);
}

void free_ivector(int *v, long nl, long nh)
{
  free((char *) (v + nl));
  //  free(&v[nl]);
}

void free_dvector(FTYPE *v, long nl, long nh)
{
  free((char *) (v + nl));
  //  free(&v[nl]);
}



void free_matrix(FTYPE **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}



void free_submatrix(FTYPE **b, long nrl, long nrh, long ncl, long nch)
{
  free((char *) (b + nrl));
}



FTYPE **convert_matrix(FTYPE *a, long nrl, long nrh, long ncl, long nch)
{
  long i, j, nrow, ncol;
  FTYPE **m;

  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;
  m = (FTYPE **) malloc( (nrow) * sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure in convert_matrix()");
  m -= nrl;
  for (i = 0, j = nrl; i <= nrow - 1; i++, j++)
    m[j] = a + ncol * i - ncl;
  return m;
}



void free_convert_matrix(FTYPE **b, long nrl, long nrh, long ncl,
                         long nch)
{
  free((char *) (b + nrl));
}

#define JMAX 100

FTYPE rtbis(FTYPE (*func)(FTYPE,FTYPE*), FTYPE *parms, FTYPE x1, FTYPE x2, FTYPE xacc)
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be refined until its accuracy is \pm xacc.
//taken from http://gpiserver.dcom.upv.es/Numerical_Recipes/bookcpdf/c9-1.pdf
{
  int j;
  FTYPE dx,f,fmid,xmid,rtb;
  f=(*func)(x1, parms);
  fmid=(*func)(x2, parms);
  if (f*fmid >= 0.0) {
    dualfprintf( fail_file, "f(%g)=%g f(%g)=%g\n", x1, f, x2, fmid );
    nrerror("Root must be bracketed for bisection in rtbis");
  }
  rtb = (f < 0.0) ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0 lies at x+dx.
  for (j=1;j<=JMAX;j++) { 
    fmid=(*func)(xmid=rtb+(dx *= 0.5),parms); //Bisection loop.
    if (fmid <= 0.0) {
      rtb=xmid;
    }
    if (fabs(dx) < xacc || fmid == 0.0) {
      return rtb;
    }
  }
  nrerror("Too many bisections in rtbis");
  return 0.0; //Never get here.
}
