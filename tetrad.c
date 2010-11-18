
#include "decs.h"

// Charles had wrong function call!
//extern int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
//		  double *w, double *work, int *lwork, int *iwork,
//		  int *liwork, int *info);
extern int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		  double *w, double *work, int *lwork, int *info);
extern int dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda,
		   double *vl, double *vu, int *il, int *iu, double *abstol, int *M,
		   double *w,
		   double *z, int *ldz, int *isuppz,
		   double *work, int *lwork, int *iwork,
		   int *liwork, int *info);
//extern double dlamch_(char *);

// declarations
static int tetlapack_func(FTYPE (*metr)[NDIM], FTYPE (*tetr)[NDIM], FTYPE eigenvalues[]);
static int compute_tetrcon_frommetric_mathematica(FTYPE (*generalmatrix)[NDIM], FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);
static int compute_tetrcon_frommetric(FTYPE (*generalmatrix)[NDIM], FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);


/* find suitable orthonormal tetrad */

// wrapper for tetr_func().  This converts metric using dxdxp so likely simpler form easier to match mathematica version with
int tetr_func_frommetric(FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[])
{
  int jj,kk;
  int ll,pp;
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE newgcov[SYMMATRIXNDIM];
  int info;

  // get MCOORD metric from PRIMECOORD metric
  idxdxprim(dxdxp, idxdxp);

  DLOOP(jj,kk){
    newgcov[GIND(jj,kk)]=0.0;
    DLOOP(ll,pp) {
      newgcov[GIND(jj,kk)] += GINDASSIGNFACTOR(jj,kk)*gcov[GIND(ll,pp)]*idxdxp[ll][jj]*idxdxp[pp][kk];
    }
  }

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d newgcov=%21.15g\n",jj,kk,newgcov[GIND(jj,kk)]);
  //  }

  // get tetrad (feed newgcov rather than gcov since with metric assume best comparison of vectors when have non-twisted metric and likely untwisted when using dxdxp)
  info=tetr_func(METRICTETRAD, newgcov, tetrcov, tetrcon, eigenvalues);

  return(info);

}

// input gcov and get out orthonormal tetrad in covariant and contravariant forms
// if metrictype, assume gcov is inputted as simplest as user could (e.g. using dxdxp's)
int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]) 
{
  FTYPE generalmatrixlower[NDIM][NDIM];
  FTYPE tmpgeneralmatrix[NDIM][NDIM] ;
  FTYPE tmpgeneralvec[NDIM] ;
  int j,k,l ;
  int info;
  int jj,kk;

  // copy to 2D space
  DLOOP(jj,kk) generalmatrixlower[jj][kk] = gcov[GIND(jj,kk)];


  if(inputtype==METRICTETRAD){
    info=compute_tetrcon_frommetric(generalmatrixlower,tetr_con,eigenvalues);
  }
  else{
    info=tetlapack_func(generalmatrixlower,tetr_con,eigenvalues);
  }

  // force tetr_con to have correct signature so when applied it gives back vector as if nothing done (i.e. Kronecker Delta) (e.g. Minkowski to Minkowski)
  DLOOP(j,k) tmpgeneralmatrix[j][k] = tetr_con[j][k];
  DLOOP(j,k) tetr_con[j][k] = 0. ;
  DLOOP(j,k) for(l=0;l<NDIM;l++) tetr_con[j][k] += mink(j,l)*tmpgeneralmatrix[l][k] ;

  // also fix eigenvalues to be consistent
  DLOOPA(j) tmpgeneralvec[j] = eigenvalues[j];
  DLOOPA(j) eigenvalues[j]=0;
  DLOOPA(j) for(l=0;l<NDIM;l++) eigenvalues[j] += mink(j,l)*tmpgeneralvec[l] ;


  // construct tetr_cov
  DLOOP(j,k) tmpgeneralmatrix[j][k] = 0. ;
  DLOOP(j,k) for(l=0;l<NDIM;l++) tmpgeneralmatrix[j][k] += tetr_con[j][l]*gcov[GIND(l,k)] ;
  DLOOP(j,k) tetr_cov[j][k] = 0. ;
  DLOOP(j,k) for(l=0;l<NDIM;l++) tetr_cov[j][k] += mink(j,l)*tmpgeneralmatrix[l][k] ;

#if 0
  /* tetr_ cov & con are inverse transposes of each other */
  DLOOP(j,k) tmpgeneralmatrix[j+1][k+1] = tetr_cov[j][k] ;
  gaussj(tmpgeneralmatrix,NDIM,NULL,0) ;
  DLOOP(j,k) tetr_con[j][k] = tmpgeneralmatrix[k+1][j+1] ;
#endif

  return(info);
}



#define LWORKSIZE MAX(NDIM*NDIM*NDIM,26*NDIM)
#define LIWORKSIZE MAX(NDIM*NDIM*NDIM,10*NDIM)
// native doubles for inside variables
// find tetrad
static int tetlapack_func(FTYPE (*metr)[NDIM], FTYPE (*tetr)[NDIM], FTYPE eigenvalues[])
{
  char jobz,uplo ;
  int n,lda,lwork,info ;
  double a[NDIM][NDIM],w[NDIM],work[LWORKSIZE] ;
  int chk;
  int j,k ;


  jobz = 'V' ;
  uplo = 'U' ;
  n = NDIM ;
  lda = NDIM ;
  lwork = LWORKSIZE ;


#if(USINGLAPACK)
  DLOOP(j,k) a[j][k] = metr[j][k] ;

  chk = dsyev_(
	       &jobz, 		/* job: 'V' -> compute eigenvectors too */
	       &uplo,		/* which part of a is stored, 'U' -> upper */
	       &n,		/* order of matrix a */
	       (double *)a,	/* matrix (row major order) */
	       &lda,		/* leading dimension of a */
	       w,		/* eigenvalues, ascending order */
	       work,		/* workspace */
	       &lwork,		/* size of workspace */
	       &info		/* successful? */
	       ) ;


  if(info>0 && 0){
    // doesn't seem to work (gives wrong results)
    // gives right eigenvalues but not right eigenvectors
    int liwork,iwork[LIWORKSIZE] ;
	  
    liwork = LIWORKSIZE ;
    // then dsyev failed for some reason, try another algorithm

    // try dsyevr:
    // http://www.gfd-dennou.org/arch/ruby/products/ruby-lapack/doc/dsy.html
    // http://www.netlib.org/lapack/double/dsyevr.f
    DLOOP(j,k) a[j][k] = metr[j][k] ;
	  
    char range = 'A';
    double vl=0;
    double vu=1E30;
    int il=1;
    int iu=NDIM;
    char cmach='s';
    //double abstol=_dlamch(&cmach);
    double abstol=100.0*NUMEPSILON;
    int M; // output
    double z[NDIM][NDIM]; // output
    int ldz=NDIM;
    int isuppz[2*NDIM]; // output
    chk = dsyevr_(
		  &jobz, 		/* job: 'V' -> compute eigenvectors too */
		  &range,
		  &uplo,		/* which part of a is stored, 'U' -> upper */
		  &n,		/* order of matrix a */
		  (double *)a,	/* matrix (row major order) */
		  &lda,		/* leading dimension of a */
		  &vl,&vu,&il,&iu,&abstol,&M,
		  w,		/* eigenvalues, ascending order */
		  (double *)z,&ldz,isuppz,
		  work,		/* workspace */
		  &lwork,		/* size of workspace */
		  iwork,		/* size of iwork */
		  &liwork,	/* working array for optimal liwork */
		  &info		/* successful? */
		  ) ;
	  
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
  DLOOP(j,k) tetr[j][k] = a[j][k]/sqrt(fabs(w[j])+SMALL)*sign(w[j]) ;

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
    DLOOP(j,k) fprintf(stderr,"chk3 %d %d %g %g\n",j,k,a[j][k],w[j]) ;
  */

  ////////////////////////////
  //
  // copy over eigenvalues
  //
  ////////////////////////////
  DLOOPA(j){
    eigenvalues[j]=w[j];
  }


  return(info);

}




// GODMARK: only applies for metric as input
// whether to debug lapack by comparing with simpler Mathematica version that doesn't work in all generality but seems to be correct
// GODMARK: for now, always in debug mode because use bootstrap technique to figure out how to reorder lapack result
// 0 : don't use yet since raw lapack isn't ordered correctly
// 1 : fake debug, bootstrap method
// 2 : true debug method
#define DEBUGLAPACK 1

static int compute_tetrcon_frommetric(FTYPE (*generalmatrix)[NDIM], FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[])
{
  FTYPE tetrconother[NDIM][NDIM];
  int jj,kk;
  int ll,pp;
  int info;
  FTYPE eigenvaluesother[NDIM],tempeigenvalues[NDIM];
  FTYPE errorold,errornew,temptetrcon[NDIM][NDIM];
  FTYPE errorlist[NDIM][NDIM];
  FTYPE error1,error2;
  int newlist[NDIM];
  FTYPE signlist[NDIM];





#if(DEBUGLAPACK && (USINGLAPACK))


  // get general tetrad and eigenvalues
  info=tetlapack_func(generalmatrix,tetrcon,eigenvalues);

  // get tetrad and eigenvalues for not-so-general metric
  // (ignore this other info)
  compute_tetrcon_frommetric_mathematica(generalmatrix, tetrconother, eigenvaluesother);

  // try to reorder spatial parts, assuming time part is always negative eigenvalue that will come as first eigenvector always
  // tetrconother[jj][jj] is correct t,r,\theta,\phi order for diagonal part
  // tetrcon[jj][kk] is not correct order in general, with jj (rows) being out of order


  // error for every row compared to every row
  DLOOP(jj,kk){// double over rows
    errorlist[jj][kk]=0.0;
    error1=0.0;
    error2=0.0;
    DLOOPA(ll){ // over columns
      error1+=fabs(tetrcon[jj][ll]-tetrconother[kk][ll]);
    }
    // compensate for potential sign flip for eigenvector
    DLOOPA(ll){ // over columns
      error2+=fabs(tetrcon[jj][ll]+tetrconother[kk][ll]);
    }
    if(error1>error2) errorlist[jj][kk]=-error2; // then use sign to indicate need to flip sign
    else errorlist[jj][kk]=error1; // same sign
  }
  // now have error matrix.


  // start with minimum error pair and increase in error until obtain all unique matches
  DLOOPA(kk) newlist[kk]=-1;
  int notfoundall=1;
  FTYPE signerror=0;
  while(notfoundall){
    int minkk=-1;
    int minjj=-1;
    FTYPE minerror=BIG;
    DLOOP(jj,kk){
      if(fabs(errorlist[jj][kk])<minerror && newlist[kk]==-1){ // only use minimum if minimum AND not already on list (forces result to be unique)
	minerror=fabs(errorlist[jj][kk]);
	minjj=jj;
	minkk=kk;
	signerror=sign(errorlist[jj][kk]);
      }
    }
    if(minjj==-1 || minkk==-1){
      dualfprintf(fail_file,"Problem in finding minimum error between eigenvectors\n");
      myexit(34643634);
    }
    else{
      newlist[minkk]=minjj; // point {t,x,y,z}=minkk to lapack (minjj)
      signlist[minkk]=signerror; // set sign to flip for this given {t,x,y,z}=minkk
    }
    notfoundall=0; // end?
    DLOOPA(kk) if(newlist[kk]==-1) notfoundall=1; // nope, not end
  }

  // store original result
  DLOOP(jj,kk) temptetrcon[jj][kk]=tetrcon[jj][kk];
  DLOOPA(jj) tempeigenvalues[jj]=eigenvalues[jj];
  
  /////////////////
  //
  // reorder
  //
  /////////////////
  DLOOPA(jj){// over rows
    eigenvalues[jj]=tempeigenvalues[newlist[jj]];
  }
  DLOOPA(jj){// over rows
    DLOOPA(kk){
      tetrcon[jj][kk]=temptetrcon[newlist[jj]][kk]*signlist[jj]; // signlist[jj] since jj is t,x,y,z in order
    }
  }

#if(DEBUGLAPACK==2)
  // real debug
  if(1){
      //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
      dualfprintf(fail_file,"\ninfo=%d\n",info);
      DLOOPA(jj){
	dualfprintf(fail_file,"jj=%d eigenorig=%21.15g eigen=%21.15g old=%21.15g\n",jj,tempeigenvalues[jj],eigenvalues[jj],eigenvaluesother[jj]);
      }
      DLOOP(jj,kk){
	dualfprintf(fail_file,"jj=%d kk=%d lapack=%21.15g old=%21.15g\n",jj,kk,tetrcon[jj][kk],tetrconother[jj][kk]);
      }
      dualfprintf(fail_file,"\n");
    
      fflush(fail_file);
      myexit(0);
    }
#endif


#else // else if not debugging





#if(USINGLAPACK)
  // input general metric
  // don't need to convert using idxdxp since result is tetrad used to convert from whatever input metric one started with
  info=tetlapack_func(generalmatrix,tetrcon,eigenvalues);
#else

  // generate orthonormal basis and dxdxp transformation (user-based function that doesn't take general metric)
  info=compute_tetrcon_frommetric_mathematica(generalmatrix, tetrcon, eigenvalues);
#endif




#endif // end over debug/not debug


  return(info);

}

// compute orthonormal tetrad \Lambda^\mu_\nu for converting u^\nu
// Notes:
// For full rotating BH with g_{r\theta}=0 (i.e. use dxdxp first), result can be found analytically:
// See: simple_eigensystem_nonrotbh.nb
// http://en.wikipedia.org/wiki/Cubic_equation
// However, Orthogonalization is quite involved, so probably simpler to do full numerical solution
// Issue always: order of system being out-of-order compared to input order
static int compute_tetrcon_frommetric_mathematica(FTYPE (*generalmatrix)[NDIM], FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[])
{
  FTYPE gtt,grr,grt,ghh,gpp;
  FTYPE blob1,blob1sq;

  // now assume no mixing in r-\theta
  // Right now, only works for a=0 KS or BL

  gtt=generalmatrix[0][0];
  grr=generalmatrix[1][1];
  grt=generalmatrix[1][0];
  ghh=generalmatrix[2][2];
  gpp=generalmatrix[3][3];

  blob1=grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt;
  blob1sq=blob1*blob1;

  //regexp: pow(\([a-zA-Z0-9-+*/ ]+\),2) -> ((\1)*(\1))

  // get orthonormal basis transformation
  tetrcon[0][0]=-(sqrt(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
		  pow((4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*blob1sq,0.25));
  
  tetrcon[0][1] =  (2.0*grt)/(sqrt(4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
			      sqrt(fabs(grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));

  tetrcon[0][2] = 0.0;
  tetrcon[0][3] = 0.0;

  tetrcon[1][0] = (2.0*grt)/sqrt((4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
				 (grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt));

  tetrcon[1][1] = sqrt((grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
		       (sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));

  tetrcon[1][2] = 0.0;
  tetrcon[1][3] = 0.0;

  tetrcon[2][0] = 0.0;
  tetrcon[2][1] = 0.0;
  tetrcon[2][2] = 1.0/fabs(sqrt(ghh));
  tetrcon[2][3] = 0.0;

  tetrcon[3][0] = 0.0;
  tetrcon[3][1] = 0.0;
  tetrcon[3][2] = 0.0;
  tetrcon[3][3] = 1.0/fabs(sqrt(gpp));

  eigenvalues[0]=0.5*(grr - 1.*sqrt(4.*pow(grt,2) + pow(grr - 1.*gtt,2)) + gtt);
  eigenvalues[1]=0.5*(grr + sqrt(4.*pow(grt,2) + pow(grr - 1.*gtt,2)) + gtt);
  eigenvalues[2]=ghh;
  eigenvalues[3]=gpp;

  return(0);
}
