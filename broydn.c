#include "decs.h"

#define NRANSI
#include "NRUTIL.H"
#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);        \
    free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);           \
    free_matrix(r,1,n,1,n);free_matrix(qt,1,n,1,n);free_vector(p,1,n);  \
    free_vector(g,1,n);free_vector(fvcold,1,n);free_vector(d,1,n);      \
    free_vector(c,1,n);return;}


void broydn(int useanalyticjac
            ,FTYPE parms[]
            ,FTYPE x[], int n, int *check
            ,void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[])
            ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
            )
{
  //  FTYPE nrfmin(FTYPE parms[], FTYPE x[]);
  // void lnsrch(int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[],
  //   FTYPE *f, FTYPE stpmax, int *check, FTYPE (*func)(FTYPE []));
  // void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
  // void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
  // void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);
  int i,its,j,k,restrt,sing,skip;
  FTYPE den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
  FTYPE *g,*p,**qt,**r,*s,*t,*w,*xold;
  // DEBUG var
  FTYPE **rtest;
  int qq,pp;
  FTYPE jacdiff,jacavg;


  c=vector(1,n);
  d=vector(1,n);
  fvcold=vector(1,n);
  g=vector(1,n);
  p=vector(1,n);
  qt=matrix(1,n,1,n);
  r=matrix(1,n,1,n);
  // DEBUG line below
  rtest=matrix(1,n,1,n);
  s=vector(1,n);
  t=vector(1,n);
  w=vector(1,n);
  xold=vector(1,n);
  fvec=vector(1,n);
  nn=n;
  nrfuncv=vecfunc;
  f=nrfmin(parms,x);
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
  if (test < 0.01*TOLF) {
    *check=0;
    FREERETURN
      }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(FTYPE)n);
  restrt=1;
  for (its=1;its<=MAXITS;its++) {
    if (restrt) {
      if(useanalyticjac){
        // GODMARK : why doesn't this work?
        //    usrfun2(n,x,fvec,rtest);
        //usrfun_joninterp2(n,x,fvec,r);
        (*usrfun)(n,parms,x,fvec,r);
      }
      else{
        fdjac(n,parms,x,fvec,r,vecfunc);
        /*
          jacdiff=0.0;
          jacavg=0.0;
          for(qq=1;qq<=n;qq++) for(pp=1;pp<=n;pp++){
          jacavg+=fabs(rtest[qq][pp]);
          jacdiff+=fabs(rtest[qq][pp]-r[qq][pp]);
          }
          jacavg/=(n*n);
          stderrfprintf("x=%g %g fvec=%g %g jacdiff=%g %g\n",x[1],x[2],fvec[1],fvec[2],jacdiff,jacdiff/jacavg);
        */
      }
      qrdcmp(r,n,c,d,&sing);
      if (sing) nrerror("singular Jacobian in broydn");
      for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) qt[i][j]=0.0;
        qt[i][i]=1.0;
      }
      for (k=1;k<n;k++) {
        if (c[k]) {
          for (j=1;j<=n;j++) {
            sum=0.0;
            for (i=k;i<=n;i++)
              sum += r[i][k]*qt[i][j];
            sum /= c[k];
            for (i=k;i<=n;i++)
              qt[i][j] -= sum*r[i][k];
          }
        }
      }
      for (i=1;i<=n;i++) {
        r[i][i]=d[i];
        for (j=1;j<i;j++) r[i][j]=0.0;
      }
    } else {
      for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
      for (i=1;i<=n;i++) {
        for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
        t[i]=sum;
      }
      skip=1;
      for (i=1;i<=n;i++) {
        for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
        w[i]=fvec[i]-fvcold[i]-sum;
        if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;
        else w[i]=0.0;
      }
      if (!skip) {
        for (i=1;i<=n;i++) {
          for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
          t[i]=sum;
        }
        for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
        for (i=1;i<=n;i++) s[i] /= den;
        qrupdt(r,qt,n,t,s);
        for (i=1;i<=n;i++) {
          if (r[i][i] == 0.0) nrerror("r singular in broydn");
          d[i]=r[i][i];
        }
      }
    }
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
      g[i]=sum;
    }
    for (i=n;i>=1;i--) {
      for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) {
      xold[i]=x[i];
      fvcold[i]=fvec[i];
    }
    fold=f;
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
      p[i] = -sum;
    }
    rsolv(r,n,d,p);
    lnsrch(n,parms,xold,fold,g,p,x,&f,stpmax,check,nrfmin);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      *check=0;
      FREERETURN
        }
    if (*check) {
      if (restrt) FREERETURN
        else {
          test=0.0;
          den=FMAX(f,0.5*n);
          for (i=1;i<=n;i++) {
            temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
            if (temp > test) test=temp;
          }
          if (test < TOLMIN) FREERETURN
            else restrt=1;
        }
    } else {
      restrt=0;
      test=0.0;
      for (i=1;i<=n;i++) {
        temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
        if (temp > test) test=temp;
      }
      if (test < TOLX) FREERETURN
                         }
  }
  nrerror("MAXITS exceeded in broydn");
  FREERETURN
    }
#undef MAXITS
#undef EPS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
