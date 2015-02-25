
/*! \file tetlapack_func_doubleonly.c
    \brief External related Tetrad calculations (double only)
*/


#include "decs.h"

/// Charles had wrong function call!
///extern int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
///    double *w, double *work, int *lwork, int *iwork,
///    int *liwork, int *info);
extern int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                  double *w, double *work, int *lwork, int *info);
extern int dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda,
                   double *vl, double *vu, int *il, int *iu, double *abstol, int *M,
                   double *w,
                   double *z, int *ldz, int *isuppz,
                   double *work, int *lwork, int *iwork,
                   int *liwork, int *info);
//extern double dlamch_(char *);


#define LWORKSIZE MAX(NDIM*NDIM*NDIM,26*NDIM)
#define LIWORKSIZE MAX(NDIM*NDIM*NDIM,10*NDIM)
/// native doubles for inside variables
/// find tetrad
int tetlapack_func(double (*metr)[NDIM], double (*tetr)[NDIM], double eigenvalues[])
{
  char jobz,uplo ;
  int n,lda,lwork,info=0 ;
  double a[NDIM][NDIM],w[NDIM]={0},work[LWORKSIZE]={0};
  int chk;
  int j,k ;


  jobz = 'V' ;
  uplo = 'U' ;
  n = NDIM ;
  lda = NDIM ;
  lwork = LWORKSIZE ;


#if(USINGLAPACK)
  DLOOP(j,k) a[j][k] = (double)metr[j][k] ;

  chk = dsyev_(
               &jobz,   /* job: 'V' -> compute eigenvectors too */
               &uplo,  /* which part of a is stored, 'U' -> upper */
               &n,  /* order of matrix a */
               (double *)a, /* matrix (row major order) */
               &lda,  /* leading dimension of a */
               w,  /* eigenvalues, ascending order */
               work,  /* workspace */
               &lwork,  /* size of workspace */
               &info  /* successful? */
               ) ;


  if(info>0 && 0){

    dualfprintf(fail_file,"issue with dsyev: info=%d\n",info);
    
    // doesn't seem to work (gives wrong results)
    // gives right eigenvalues but not right eigenvectors
    int liwork,iwork[LIWORKSIZE] ;
   
    liwork = LIWORKSIZE ;
    // then dsyev failed for some reason, try another algorithm

    // try dsyevr:
    // http://www.gfd-dennou.org/arch/ruby/products/ruby-lapack/doc/dsy.html
    // http://www.netlib.org/lapack/double/dsyevr.f
    DLOOP(j,k) a[j][k] = (double)metr[j][k] ;
   
    char range = 'A';
    double vl=0;
    double vu=1E30; // a big double
    int il=1;
    int iu=NDIM;
    char cmach='s';
    //double abstol=_dlamch(&cmach);
#define NUMEPSILONDBL ((double)(2.2204460492503131e-16))
    double abstol=(double)(100.0*NUMEPSILONDBL); // has to stay as double precision
    int M; // output
    double z[NDIM][NDIM]; // output
    int ldz=NDIM;
    int isuppz[2*NDIM]; // output
    chk = dsyevr_(
                  &jobz,   /* job: 'V' -> compute eigenvectors too */
                  &range,
                  &uplo,  /* which part of a is stored, 'U' -> upper */
                  &n,  /* order of matrix a */
                  (double *)a, /* matrix (row major order) */
                  &lda,  /* leading dimension of a */
                  &vl,&vu,&il,&iu,&abstol,&M,
                  w,  /* eigenvalues, ascending order */
                  (double *)z,&ldz,isuppz,
                  work,  /* workspace */
                  &lwork,  /* size of workspace */
                  iwork,  /* size of iwork */
                  &liwork, /* working array for optimal liwork */
                  &info  /* successful? */
                  ) ;

    if(info!=0) dualfprintf(fail_file,"issue with dsyevr: info=%d\n",info);
   
  }





#else
  info=-1;
  dualfprintf(fail_file,"No LAPACK!\n");
  myexit(873483746);
#endif

  // note a[j][k] corresponds to a[j=which eigenvector][k=which component]
  // however order of j is that of order of eigenvalues that is ascending

  // sign of eigenvalue is necessary to get signature of eigenvector correct so that basis has same handedness as original basis
  // the below fixes the time-component in test, but not r-component where r-r term was negative what it should be
  DLOOP(j,k) tetr[j][k] = (double)(a[j][k]/sqrt(fabs(w[j])+SMALL)*sign(w[j]));

  // still, sign isn't right
  // ensure signature same as minkowski
  DLOOPA(j){
    if(sign(tetr[j][j])!=sign(mink(j,j))){
      // no, in simplest case tetr is Kronecker delta, not Minkowski, so that nothing happens to vector if already in Minkowski
      //    if(sign(tetr[j][j])!=sign(1.0)){
      // well, maybe still so
      DLOOPA(k) tetr[j][k]*=-1.0;
    }
  }

  /*
    DLOOP(j,k) stderrfprintf("chk3 %d %d %g %g\n",j,k,a[j][k],w[j]) ;
  */

  ////////////////////////////
  //
  // copy over eigenvalues
  //
  ////////////////////////////
  DLOOPA(j){
    eigenvalues[j]=(double)w[j];
  }


  return(info);

}
