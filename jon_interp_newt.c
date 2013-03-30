#include "decs.h"
#include "NRUTIL.H"

#define NRANSI

// modifiable
#define MAXITS 200
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0



// for debugging
/*
  #define FREERETURN {fprintf(stderr,"1addr=%d\n",fvec);free_vector(fvec,1,n);fprintf(stderr,"1 n=%d\n",n); \
  stderrfprintf("2addr=%ld\n",xold);free_vector(xold,1,n);fprintf(stderr,"2 n=%d\n",n); \
  stderrfprintf("3addr=%ld\n",p);free_vector(p,1,n);fprintf(stderr,"3 n=%d\n",n);   \
  stderrfprintf("4addr=%ld\n",g);free_vector(g,1,n);fprintf(stderr,"1 n=%d\n",n);   \
  stderrfprintf("5addr=%ld\n",fjac);free_matrix(fjac,1,n,1,n);fprintf(stderr,"1 n=%d\n",n);  \
  stderrfprintf("6addr=%ld\n",indx);free_ivector(indx,1,n);fprintf(stderr,"f n=%d\n",n);\
  return;}
*/
 
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);        \
    free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);    \
    free_ivector(indx,1,n);                                             \
    return;}


/*
  #define TESTDEATH {fprintf(stderr,"1addr=%d\n",fvec);free_vector(fvec,1,n);fprintf(stderr,"1 n=%d\n",n); \
  stderrfprintf("2addr=%ld\n",xold);free_vector(xold,1,n);fprintf(stderr,"2 n=%d\n",n);\
  stderrfprintf("3addr=%ld\n",p);free_vector(p,1,n);fprintf(stderr,"3 n=%d\n",n);\
  stderrfprintf("4addr=%ld\n",g);free_vector(g,1,n);fprintf(stderr,"1 n=%d\n",n);\
  stderrfprintf("5addr=%ld\n",fjac);free_matrix(fjac,1,n,1,n);fprintf(stderr,"1 n=%d\n",n);\
  stderrfprintf("6addr=%ld\n",indx);free_ivector(indx,1,n);fprintf(stderr,"f n=%d\n",n);\
  exit(0);}
*/

void newt(int useanalyticjac
          ,FTYPE parms[]
          ,FTYPE x[], int n, int *check
          ,void (*vecfunc)(int n, FTYPE *parms, FTYPE v[], FTYPE f[])
          ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
          )
{
  //  void fdjac(int n, FTYPE parms[], FTYPE x[], FTYPE fvec[], FTYPE **df,
  //      void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]));
  // void lubksb(FTYPE **a, int n, int *indx, FTYPE b[]);
  //int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
  int i,j,its,*indx;
  FTYPE d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
  // DEBUG VARS
  int qq,pp;




  indx=ivector(1,n);
  fjac=matrix(1,n,1,n);
  g=vector(1,n);
  p=vector(1,n);
  xold=vector(1,n);
  fvec=vector(1,n);
  nn=n;
  nrfuncv=vecfunc;
  f=nrfmin(parms,x);

  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test < 0.01*TOLF) {
    *check=0;
    //stderrfprintf("end tolf0.01\n");
    FREERETURN
      }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(FTYPE)n);


  for (its=1;its<=MAXITS;its++) {
    if(useanalyticjac){
      (*usrfun)(n,parms,x,fvec,fjac);
    }
    else{
      fdjac(n,parms,x,fvec,fjac,vecfunc);
    }


    // DEBUG
    //   for (i=1;i<=n;i++)for (j=1;j<=n;j++)fprintf(stderr,"fjac[%d][%d]=%21.15g\n",i,j,fjac[i][j]);
 

    ///////////////////////////DEBUG
    /*
      for(qq=1;qq<n+1;qq++) for(pp=1;pp<n+1;pp++){
      stderrfprintf("fjac[%d][%d]=%g\n",qq,pp,fjac[qq][pp]);
      }
      fflush(stderr);
      //exit(0);
      */  
    ///////////////////////////DEBUG
    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;
    for (i=1;i<=n;i++) p[i] = -fvec[i];
    ludcmp(fjac,n,indx,&d);
    //  if(its==5) exit(0);

    // DEBUG
    //  for (i=1;i<=n;i++)fprintf(stderr,"indx[%d]=%d\n",i,indx[i]);

    for (i=1;i<=n;i++){
      if(indx[i]<1 || indx[i]>n){
        stderrfprintf("Bad indx[%d]=%d\n",i,indx[i]);
        exit(1);
      }
    }



    // check for singularity
    int sing;
    sing=0;
    int nonzeroelements;

    for (i=1;i<=n;i++){
      nonzeroelements=0;
      for(j=1;j<=n;j++){
        if(fabs(fjac[i][j])>SMALL) nonzeroelements++;
      }
      if(nonzeroelements==0) sing++;
    }
    if(sing>0){
      // found singularity
      *check=2;
      FREERETURN;
    }




    lubksb(fjac,n,indx,p);
    lnsrch(n,parms,xold,fold,g,p,x,&f,stpmax,check,nrfmin);



    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      *check=0;
      //stderrfprintf("end tolf\n");
      FREERETURN
        }
    if (*check) {
      test=0.0;
      den=FMAX(f,0.5*n);
      for (i=1;i<=n;i++) {
        temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
        if (temp > test) test=temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      if(SIMPLEDEBUGINTERP) stderrfprintf("end check\n");
      FREERETURN
        }
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX){
      //    stderrfprintf("end tolx\n");
      FREERETURN
        }
  }
  //nrerror("MAXITS exceeded in newt");
  fprintf( stderr, "MAXITS exceeded in newt");
  FREERETURN;
}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
