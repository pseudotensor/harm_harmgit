
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
  // Tcon_j^l = Tconorig^{al} \eta_{ja}
  // So first index is orthonormal and second index is old coordinate basis
  DLOOP(j,k) tmpgeneralmatrix[j][k] = tetr_con[j][k];
  DLOOP(j,k) tetr_con[j][k] = 0. ;
  DLOOP(j,k) for(l=0;l<NDIM;l++) tetr_con[j][k] += mink(j,l)*tmpgeneralmatrix[l][k] ;

  // also fix eigenvalues to be consistent
  DLOOPA(j) tmpgeneralvec[j] = eigenvalues[j];
  DLOOPA(j) eigenvalues[j]=0;
  DLOOPA(j) for(l=0;l<NDIM;l++) eigenvalues[j] += mink(j,l)*tmpgeneralvec[l] ;


  // construct tetr_cov
  // Tcov^a_b = Tcon_j^l g_{lb} \eta^{aj}
  // So first index is orthonormal and second index is old coordinate basis.
  // So Tcon and Tcov are inverses but not transposes of each other as would occur if matrix_inverse() operated.
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
    DLOOP(j,k) stderrfprintf("chk3 %d %d %g %g\n",j,k,a[j][k],w[j]) ;
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

// compute orthonormal tetrad \Lambda_\mu[ortho]^\nu[lab] for converting u_\nu[lab] or u^\nu[ortho]
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
  // But, even for a\neq 0, the order is found correctly since spherical polar and KS gtr terms dominate

  if(THETAROT!=0.0){
    dualfprintf(fail_file,"compute_tetrcon_frommetric_mathematica() needs to have rotation added (see jon_interp stuff) so can handle tilts that mix theta and phi.  Then can trust picking up dominate terms.\n");
    myexit(546292153);
  }

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



////////////
//
// Notes:
//
// dxdxp = dV^i/dX^k = \Lambda^i_k
// idxdxp = inverse and tranpose of dxdxp = dX^k/dV^i = (iLambda)^k_i
// Tetrcon_k^j [first index ortho, second index lab]
// Tetrcov^k_j [first index ortho, second index lab (i.e. not transposed!)]


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of ORTHO to transform lab <--> ORTHO
// tmuup : LAB2ORTHO
// tmudn : ORTHO2LAB
int calc_ORTHOes(struct of_geom *ptrgeom, FTYPE tmuup[][NDIM], FTYPE tmudn[][NDIM])
{
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE tetrcovV[NDIM][NDIM];
  FTYPE tetrconV[NDIM][NDIM];
  FTYPE eigenvaluesV[NDIM];
  FTYPE tetrcovX[NDIM][NDIM];
  FTYPE tetrconX[NDIM][NDIM];
  int jj,kk,ll;


  // get dxdxp
  dxdxprim_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,dxdxp);
  idxdxprim(dxdxp, idxdxp);

  // get tetrad (uses dxdxp so that tetrcon and tetrcon and eigenvalues are using V metric not X metric
  tetr_func_frommetric(dxdxp, ptrgeom->gcov, tetrcovV, tetrconV, eigenvaluesV);

  // now convert back to X metric for general internal code use

  // \LambdaXcov^jj[ortho]_kk[labX] = \LambdaVcov^jj[ortho]_ll[labV] idxdxp^ll[labX]_kk[labV]
  DLOOP(jj,kk){
    tetrcovX[jj][kk]=0.0;
    DLOOPA(ll) {
      tetrcovX[jj][kk] += tetrcovV[jj][ll]*idxdxp[kk][ll];
    }
  }

  // \LambdaXcon_jj[ortho]^kk[labX] = \LambdaVcon_jj[ortho]^ll[labV] dxdxp^ll[labV]_kk[labX]
  DLOOP(jj,kk){
    tetrconX[jj][kk]=0.0;
    DLOOPA(ll) {
      tetrconX[jj][kk] += tetrconV[jj][ll]*dxdxp[ll][kk];
    }
  }

  // map to tmuup and tmudn
  DLOOP(jj,kk) tmuup[jj][kk]=tetrcovX[jj][kk]; // \Lambda^jj[ortho]_kk[labX] = tmuup = "LAB2ORTHO" = tetrcovX
  DLOOP(jj,kk) tmudn[jj][kk]=tetrconX[jj][kk]; // \Lambda_jj[ortho]^kk[labX] = tmudn = "ORTHO2LAB" = tetrconX

  return(0);
}


// get stored or compute tetrcov and tetrcon
int get_tetrcovcon(struct of_geom *ptrgeom, FTYPE (**tetrcov)[NDIM],FTYPE (**tetrcon)[NDIM])
{

  ////////////////////////////////////
  // for tetrads
  ////////////////////////////////////
  int ii,jj,kk,pp;
  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;
  pp=ptrgeom->p;

  // get tetrads
  // tmuup ~ tlab2ortho[LAB2ORTHO] ~ tetrcon [which operates on covariant things when in form: ufl_\nu = Tetrcon_\nu^\mu ucovlab_\mu]
  // tmudn ~ tlab2ortho[ORTHO2LAB] ~ tetrcov [which operates on contravariant things when in form: ufl^\nu = Tetrcov^\nu_\mu uconlab^\mu
  // i.e. first and second index each always refer to the same basis without transposition.  So always contract with second index to get thing into form of first index.

  // note that ORTHO2LAB is stored as tetrcon=tmudn
  // note that LAB2ORTHO is stored as tetrcov=tmuup
  if(STORETLAB2ORTHO==1){ // global parameter
    // then just get stored version instead of computing on fly
    *tetrcov=GLOBALMETMACP2A0(tlab2ortho,pp,LAB2ORTHO,ii,jj,kk);
    *tetrcon=GLOBALMETMACP2A0(tlab2ortho,pp,ORTHO2LAB,ii,jj,kk);
  }
  else{// compute on fly here

    // compute (expensive in general!)
    // NOTE ORDER!
    calc_ORTHOes(ptrgeom, *tetrcov, *tetrcon);
  }

  return(0);

}




//**********************************************************************
// Bardeen tensor transforming between ZAMO and LAB frames
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of ZAMO to transform lab <--> ZAMO
// Here, ZAMO is ZAMO frame, and this function only works for BL coords
// Only applies for Boyer-Lindquist coordinates
// KORALTODO: Not used.  Remove? 
int calc_ZAMOes_old(struct of_geom *ptrgeom, FTYPE emuup[][NDIM], FTYPE emudn[][NDIM])
{
  FTYPE e2nu,e2psi,e2mu1,e2mu2,omega;
  FTYPE gtt,gtph,gphph,grr,gthth;
  int i,j;
  // recast as [NDIM][NDIM] matrix
  //  FTYPE emuup[][NDIM]=(FTYPE (*)[NDIM])(&ptremuup[0]);
  //FTYPE emudn[][NDIM]=(FTYPE (*)[NDIM])(&ptremudn[0]);


  gtt=ptrgeom->gcov[GIND(0,0)];
  gtph=ptrgeom->gcov[GIND(0,3)];
  gphph=ptrgeom->gcov[GIND(3,3)];
  grr=ptrgeom->gcov[GIND(1,1)];
  gthth=ptrgeom->gcov[GIND(2,2)];

  //Bardeen's 72 coefficients:
  e2nu=-gtt+gtph*gtph/gphph;
  e2psi=gphph;
  e2mu1=grr;
  e2mu2=gthth;
  omega=-gtph/gphph;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	emuup[i][j]=0.;
	emudn[i][j]=0.;
      }

  emuup[0][0]=sqrt(e2nu);
  emuup[1][1]=sqrt(e2mu1);
  emuup[2][2]=sqrt(e2mu2);
  emuup[0][3]=-omega*sqrt(e2psi);
  emuup[3][3]=sqrt(e2psi);

  emudn[3][0]=omega*1./sqrt(e2nu);
  emudn[0][0]=1./sqrt(e2nu);
  emudn[1][1]=1./sqrt(e2mu1);
  emudn[2][2]=1./sqrt(e2mu2);
  emudn[3][3]=1./sqrt(e2psi);

  // need to return Kerr-Schild prime transformations

  return 0;
}


// Compute general Lorentz boost for an arbitrary metric and arbitrary 4-velocities
// see docs/lorentz.ps.gz by Avery Broderick
//    Lambda^\mu[w]_\nu[u] u^\nu = w^\mu
//    Lambda^\mu[w]_\nu[u] w_\mu = u_\nu
// (iLambda)^\mu[u]_\nu[w] u_\mu = w_\nu
// (iLambda)^\mu[u]_\nu[w] w^\nu = u^\mu
// So if going from u->w (i.e. FF2LAB) frame for  vecff^\mu  , then apply    Lambda^\mu_\nu vecff^nu  = veclab^\mu
// So if going from w->u (i.e. LAB2FF) frame for veclab^\mu  , then apply (iLambda)^\mu_\nu veclab^nu = vecff^\mu
// All this assumes iLambda was formed from matrix_inverse() that gives inverse transpose
int calc_generalized_boost_uu(struct of_geom *ptrgeom, FTYPE *wcon, FTYPE *ucon, FTYPE (*lambda)[NDIM])
{
  int mu,nu;
  FTYPE wcov[NDIM],ucov[NDIM];

  // wcov
  DLOOPA(mu) wcov[mu] = 0.0;
  DLOOP(mu,nu) wcov[nu] += wcon[mu]*(ptrgeom->gcov[GIND(mu,nu)]);

  // ucov
  DLOOPA(mu) ucov[mu] = 0.0;
  DLOOP(mu,nu) ucov[nu] += ucon[mu]*(ptrgeom->gcov[GIND(mu,nu)]);

  // gamma
  FTYPE gamma=0.0;
  DLOOPA(mu) gamma += -wcon[mu]*ucov[mu];

  // lambda
  DLOOP(mu,nu) lambda[mu][nu]= 
    + delta(mu,nu) 
    + (1.0/(1.0+gamma)) * (wcon[mu]*wcov[nu] + ucon[mu]*ucov[nu] - gamma *(ucon[mu]*wcov[nu] + wcon[mu]*ucov[nu]))
    + (ucon[mu]*wcov[nu] - wcon[mu]*ucov[nu]);

  return(0);
}


// calculate boost assuming in orthonormal basis where g_{\mu\nu}=\eta_{\mu\nu} = diag(-1,0,0,0)
// NOTEMARK: Orthonormal does *not* mean in Cartesian coordinates!  If metric was originally in SPC Mink or SPC KS, then final metric is diag(-1,0,0,0) but is still SPC.  Have to do SPC->Cart conversion (see, e.g., jon_interp stuff) to get to Cartesian.
int calc_ortho_boost_uu(FTYPE *wcon, FTYPE *ucon, FTYPE (*lambda)[NDIM])
{
  int mu,nu;
  FTYPE wcov[NDIM],ucov[NDIM];

  // wcov
  DLOOPA(mu) wcov[mu] = wcon[mu];
  wcov[TT]*=-1.0;

  // ucov
  DLOOPA(mu) ucov[mu] = ucon[mu];
  ucov[TT]*=-1.0;

  // gamma
  FTYPE gamma=0.0;
  DLOOPA(mu) gamma += -wcon[mu]*ucov[mu];

  // lambda
  DLOOP(mu,nu) lambda[mu][nu]= 
    + delta(mu,nu) 
    + (1.0/(1.0+gamma)) * (wcon[mu]*wcov[nu] + ucon[mu]*ucov[nu] - gamma *(ucon[mu]*wcov[nu] + wcon[mu]*ucov[nu]))
    + (ucon[mu]*wcov[nu] - wcon[mu]*ucov[nu]);

  return(0);
}



// use lab frame contravariant 4-velocity (uconlab) and get transformation matrix for going to orthonormal basis (same base coordinate system: e.g. SPC, does not convert to Cartesian) or back
// NOTEMARK: If set uconlab=uconZAMO, then no boost and just does Xlab2Vortho
int transboost_lab2fluid(struct of_geom *ptrgeom, FTYPE *uconlab, FTYPE (*transboostup)[NDIM], FTYPE (*transboostlo)[NDIM])
{
  int mu,nu;


  ///////////////////////////////
  // get tetrcov and tetrcon
  ///////////////////////////////
  FTYPE tetrconmem[NDIM][NDIM],tetrcovmem[NDIM][NDIM];
  FTYPE (*tetrcon)[NDIM]=tetrconmem;
  FTYPE (*tetrcov)[NDIM]=tetrcovmem;
  get_tetrcovcon(ptrgeom, &tetrcov,&tetrcon); // pass address of pointer since want to give new pointer address if stored
  

  // get ucovlab
  FTYPE ucovlab[NDIM];
  DLOOPA(mu) ucovlab[mu] = 0.0;
  DLOOP(mu,nu) ucovlab[nu] += uconlab[mu]*(ptrgeom->gcov[GIND(mu,nu)]);


  // set wconlab to LAB frame, which happens to be ZAMO for the equations HARM solves
  // ZAMO frame: \eta_\mu = (-\alpha,0,0,0)
  FTYPE wcovlab[NDIM];
  wcovlab[TT]=-ptrgeom->alphalapse;
  SLOOPA(mu) wcovlab[mu]=0.0;
  // raise to get wconlab = \eta^\mu
  FTYPE wconlab[NDIM];
  DLOOPA(mu) wconlab[mu] = 0.0;
  DLOOP(mu,nu) wconlab[nu] += wcovlab[mu]*(ptrgeom->gcon[GIND(mu,nu)]);
  

  // get orthonormal boost to fluid frame
  // If convert ucon^i[lab] to orthonormal basis, then can just use simple Lorentz boost from special relativity to operate on the tetrad to have something that converts lastly to the fluid frame.
  FTYPE lambda[NDIM][NDIM];
  calc_ortho_boost_uu(wconlab, uconlab, lambda);
  // get inverse lambda
  FTYPE ilambda[NDIM][NDIM];
  // comments for matrix_inverse() say takes lambda^j_k and pops out (ilambda)^k_j such that (lambda)^j_k (ilambda)^k_l = \delta^j_l
  // So need to apply ilambda correctly assuming does transpose
  matrix_inverse(PRIMECOORDS,lambda,ilambda);


  // From earlier in other functions:
  //
  // Tetrcon_k^j [first index ortho, second index lab]
  // Tetrcov^k_j [first index ortho, second index lab (i.e. not transposed!)]
  //
  //    Lambda^\mu[w]_\nu[u] u^\nu = w^\mu
  //    Lambda^\mu[w]_\nu[u] w_\mu = u_\nu
  // (iLambda)^\mu[u]_\nu[w] u_\mu = w_\nu
  // (iLambda)^\mu[u]_\nu[w] w^\nu = u^\mu
  // i.e. iLambda is inverse and transposed
  // So if going from u->w (i.e. FF2LAB) frame for  vecff^\mu  , then apply    Lambda^\mu_\nu vecff^nu  = veclab^\mu
  // So if going from w->u (i.e. LAB2FF) frame for veclab^\mu  , then apply (iLambda)^\mu_\nu veclab^nu = vecff^\mu


  // form transboost
  int aa;
  // apply Lorentz boost on transformation from lab-frame to orthonormal basis

  // TBup_\mu[ff ortho]^\nu[lab coordbasis] = lambda^aa[lab ortho]_\mu[ff ortho] Tetrcon_aa[lab ortho]^\nu[lab coordbasis]
  DLOOP(mu,nu) transboostup[mu][nu]=0.0;
  DLOOP(mu,nu) DLOOPA(aa) transboostup[mu][nu] += lambda[aa][mu]*tetrcon[aa][nu]; // as if operated boost on w_aa

  // TBlo^\mu[ff ortho]_\nu[lab coordbasis] =  ilambda^\mu[ff ortho]_aa[lab ortho] Tetrcov^aa[lab ortho]_\nu[lab coordbasis]
  DLOOP(mu,nu) transboostlo[mu][nu]=0.0;
  DLOOP(mu,nu) DLOOPA(aa) transboostlo[mu][nu] += ilambda[mu][aa]*tetrcov[aa][nu]; // as if operated boost on w^aa

  // for starting in lab frame coordinate basis, use as follows:
  // i.e. ucon^\mu[ff ortho] =  TBlo^\mu[ff ortho]_\nu[lab coordbasis] ucon^\nu[lab coordbasis]
  // i.e. ucov_\mu[ff ortho] =  TBup_\mu[ff ortho]^\nu[lab coordbasis] ucov_\nu[lab coordbasis]

  // if starting with fluid frame orthonormal basis, then use same things but indices are contracted reversely.
  // i.e. ucon^\nu[lab coordbasis] = ucon^\mu[ff ortho]  TBup_\mu[ff ortho]^\nu[lab coordbasis]
  // i.e. ucov_\nu[lab coordbasis] = ucov_\mu[ff ortho]  TBlo^\mu[ff ortho]_\nu[lab coordbasis]

  // that is, as constructed, the TBlo and TBhi have indices always in the same consistent order as far as the meaning of the first and second index with respect to the frame and coordinates.  Unlike those things that use matrix_inverse like idxdxp and ilambda


  return(0);
}




// Wrapper for vector_lab2orthofluidorback()
//
// Correct for HARM tensor quantities that have implicit component set to (e.g.) t, which makes it ambiguous what frame that was measured in.
//
// whichvector:
//
// TYPEUCOV :
//  T^t_\nu type,  so that T^t_\nu  = -[E_\nu]/(-\alpha) = -[-\eta_\mu T^\mu_\nu] /(-\alpha) . So that E_t is negative definite in Minkowski.
//  T^{t\nu} type, so that T^{t\nu} = -[E^\nu]/(-\alpha) = -[-\eta_\mu T^{\mu\nu}]/(-\alpha) . So that E^t is positive definite in Minkowski.

// TYPEUCON:
//  B^i = *F^{it} type, so that B^i = *F^{it}[HARMLAB]  = -[B^\nu]/(-\alpha) = -[\eta_\mu *F^{\mu\nu}]/(-\alpha)
//  B_i = *F_i^t type , so that B_i = *F_i^t[HARMLAB]   = -[B_\nu]/(-\alpha) = -[\eta_\mu *F^\mu_\nu] /(-\alpha)
int vector_harm2orthofluidorback(int whichvector, int harm2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, FTYPE v4concovtype, FTYPE *vector4in, FTYPE *vector4out)
{
  int vector_lab2orthofluidorback(int lab2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, FTYPE v4concovtype, FTYPE *vector4in, FTYPE *vector4out);
  int jj;
  FTYPE vector4incopy[NPR];


  // preserve vector4in
  DLOOPA(jj) vector4incopy[jj]=vector4in[jj];


  // harmlab to ortho fluid
  if(harm2orthofluid==LAB2FF){
    // correct "t" component
    if(whichvector==TYPEUCOV){ // u_\mu type
      DLOOPA(jj) vector4incopy[jj] *= (ptrgeom->alphalapse); // gets E_\nu from T^t_\nu or T^{t\nu} from E^\nu
    }
    else if(whichvector==TYPEUCON){ // u^\nu type
      DLOOPA(jj) vector4incopy[jj] *= (ptrgeom->alphalapse); // gets B^\nu from B^i or B_i from B_\nu
    }

    // transform+boost
    vector_lab2orthofluidorback(harm2orthofluid, ptrgeom, uconcovtype, uconcov, v4concovtype, vector4incopy, vector4out);
    // vector4out is now orthonormalized and boosted into fluid frame
  }


  // ortho fluid to harmlab
  if(harm2orthofluid==FF2LAB){

    // transform+boost
    vector_lab2orthofluidorback(harm2orthofluid, ptrgeom, uconcovtype, uconcov, v4concovtype, vector4incopy, vector4out);
    // vector4out is now orthonormalized and boosted into fluid frame


    // correct "t" component
    if(whichvector==TYPEUCOV){ // u_\mu type
      DLOOPA(jj) vector4out[jj] /= (ptrgeom->alphalapse); // gets T^t_\nu from E_\nu or T^{t\nu} from E^\nu
    }
    else if(whichvector==TYPEUCON){ // u^\nu type
      DLOOPA(jj) vector4out[jj] /= (ptrgeom->alphalapse); // gets B^i from B^\nu or B_i from B_\nu
    }

  }
  
  
  return(0);

}




// convert 4-vector to/from lab-frame coordinate basis to fluid frame orthonormal basis
//
// lab2orthofluid : LAB2FF = then vector4 is lab and vector4out is fluid ortho.  FF2LAB = orthofluid 2 lab
// uconcovtype: TYPEUCON=uconcov is ucon or TYPEUCOV=uconcov is ucov for fluid 4-velocity
// uconcov is in lab-frame as 4-vector of fluid
// v4concovtype: TYPEUCON=vector4 is ucon or TYPEUCOV=vector4 is ucov for 4-vector to transform
// vector4in : inserted 4-vector to transform
// vector4out : returned orthonormal fluid frame 4-vector (returned as same concov as input vector4in)
// NOTEMARK: If insert uconcov as ZAMO, then no boost, so can then use this for just lab2ortho and back
// NOTES:
// 1) "Lab" corresponds to value of uconcov in \eta_\mu = (-\alpha,0,0,0) frame.
// 2) "HARMLAB" corresponds to "frame" used by harm.
// E.g.
//       \rho_harmlab = -\eta_\mu \rho_0 u^\mu /\alpha
//       E_\nu[harm] = -\eta_\mu T^\mu_\nu/\alpha [MA or EM or RAD]
//       B^\nu[harm] = +\eta_\mu *F^{\mu\nu}/\alpha  [i.e. B^i [lab] = *F^{it} is our choice of sign for the magnetic field]
//       [[Note these are without \sqrt{-g}]]
//
int vector_lab2orthofluidorback(int lab2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, FTYPE v4concovtype, FTYPE *vector4in, FTYPE *vector4out)
{
  int mu,nu;
  FTYPE transboostup[NDIM][NDIM],transboostlo[NDIM][NDIM];
  FTYPE ucon[NDIM]; // needed for transboost

  if(uconcovtype==TYPEUCON){ // then uconcov is contravariant u^\mu of fluid frame
    // ucon
    DLOOPA(mu) ucon[mu] = uconcov[mu];
  }
  else if(uconcovtype==TYPEUCOV){ // then uconcov is covariant u_\mu of fluid frame
    // ucon
    DLOOPA(mu) ucon[mu] = 0.0; DLOOP(mu,nu) ucon[nu] += uconcov[mu]*(ptrgeom->gcon[GIND(mu,nu)]);
  }
  else{
    dualfprintf(fail_file,"No such uconcovtype=%d\n",uconcovtype);
    myexit(934627520);
  }


  // get trans boosts (uses ucon always, hence above getting of ucon)
  transboost_lab2fluid(ptrgeom, ucon, transboostup, transboostlo);


  // apply trans boost to 4-vector
  DLOOPA(mu) vector4out[mu]=0.0;
  if(v4concovtype==TYPEUCON){ // then vector4in^\mu is contravariant
    if(lab2orthofluid==LAB2FF){
      // vector4ff^\nu = TBlo^\nu_\mu vector4labcoordbasis^\mu
      DLOOP(mu,nu) vector4out[nu] += transboostlo[nu][mu]*vector4in[mu]; // application on vector4in^\mu
    }
    else{ // i.e. orthoff 2 lab
      // vector4labcoordbasis^\mu = vector4ff^\nu TBup_\nu^\mu 
      DLOOP(mu,nu) vector4out[mu] += vector4in[nu]*transboostup[nu][mu]; // application on vector4in^\nu
    }
  }
  else if(v4concovtype==TYPEUCOV){ // then vector4in_\mu is covariant
    if(lab2orthofluid==LAB2FF){
      //vector4ff_\nu = TBup_\nu^\mu vector4labcoordbasis_\mu
      DLOOP(mu,nu) vector4out[nu] += transboostup[nu][mu]*vector4in[mu]; // application on vector4in_\mu
    }
    else{ // i.e. orthoff 2 lab
      // vector4labcoordbasis_\mu = vector4ff_\nu TBlo^\nu_\mu 
      DLOOP(mu,nu) vector4out[mu] += vector4in[nu]*transboostlo[nu][mu]; // application on vector4in_\nu
    }
  }
  else{
    dualfprintf(fail_file,"No such v4concovtype=%d\n",v4concovtype);
    myexit(934627520);
  }

  return(0);

}



// convert lab frame 4-tensor to fluid frame orthonormal 4-tensor
//
// lab2orthofluid: LAB2FF : lab2orthofluid   FF2LAB: orthofluid 2 lab
// uconcovtype: TYPEUCON=uconcov is ucon or TYPEUCOV=uconcov is ucov for fluid 4-velocity
// uconcov is in lab-frame as 4-vector of fluid
// tconcovtypeA: for 1st index of tensor: TYPEUCON=tconcov is con or TYPEUCOV=tconcov is cov
// tconcovtypeB: for 2nd index of tensor: TYPEUCON=tconcov is con or TYPEUCOV=tconcov is cov
// tensor4in : input lab-frame as 4-tensor
// tensor4out is returned orthonormal fluid frame 4-tensor (same tconcovtypeA and tconcovtypeB as tensor4in)
// NOTEMARK: If insert uconcov as ZAMO, then no boost, so can then use this for just lab2ortho and back
int tensor_lab2orthofluidorback(int lab2orthofluid, struct of_geom *ptrgeom, int uconcovtype, FTYPE *uconcov, int tconcovtypeA, int tconcovtypeB, FTYPE (*tensor4in)[NDIM], FTYPE (*tensor4out)[NDIM])
{
  int mu,nu,aa,bb;
  FTYPE ucon[NDIM];
  FTYPE transboostup[NDIM][NDIM],transboostlo[NDIM][NDIM];


  if(uconcovtype==TYPEUCON){ // then uconcov is contravariant u^\mu of fluid frame
    // ucon
    DLOOPA(mu) ucon[mu] = uconcov[mu];
  }
  else if(uconcovtype==TYPEUCOV){ // then uconcov is covariant u_\mu of fluid frame
    // ucon
	raise_vec(uconcov,ptrgeom,ucon);
  }
  else{
    dualfprintf(fail_file,"No such uconcovtype=%d\n",uconcovtype);
    myexit(934627525);
  }

  // get trans boosts (uses ucon always, hence above getting of ucon)
  transboost_lab2fluid(ptrgeom, ucon, transboostup, transboostlo);


  // apply trans boost to 4-tensor
  DLOOP(mu,nu) tensor4out[mu][nu]=0.0;
  if(tconcovtypeA==TYPEUCON && tconcovtypeB==TYPEUCON){
    if(lab2orthofluid==LAB2FF){
      // tfl^{\nu bb} = TBlo^\nu[ffortho]_\mu[labcoord] TBlo^bb[ffortho]_aa[labcoord] t^{\mu aa}
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[nu][bb] += transboostlo[nu][mu]*transboostlo[bb][aa]*tensor4in[mu][aa]; // application on con con
    }
    else{
      // t^{\mu aa} = tfl^{\nu bb} TBup_\nu^\mu TBup_bb^aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[mu][aa] += tensor4in[nu][bb]*transboostup[nu][mu]*transboostup[bb][aa]; // application on con con
    }
  }
  else if(tconcovtypeA==TYPEUCON && tconcovtypeB==TYPEUCOV){
    if(lab2orthofluid==LAB2FF){
      // tfl^\nu_bb = TBlo^\nu_\mu TBup_bb^aa t^\mu_aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[nu][bb] += transboostlo[nu][mu]*transboostup[bb][aa]*tensor4in[mu][aa]; // application on con cov
    }
    else{
      // t^\mu_aa = t^\nu_bb TBup_\nu^\mu TBlo^bb_aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[mu][aa] += tensor4in[nu][bb]*transboostup[nu][mu]*transboostlo[bb][aa]; // application on con cov
    }
  }
  else if(tconcovtypeA==TYPEUCOV && tconcovtypeB==TYPEUCON){
    if(lab2orthofluid==LAB2FF){
      // tfl_\nu^bb = TBup_\nu^\mu TBlo^bb_aa t_\mu^aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[nu][bb] += transboostup[nu][mu]*transboostlo[bb][aa]*tensor4in[mu][aa]; // application on cov con
    }
    else{
      // t_\mu^aa = t_\nu^bb TBlo^\nu_\mu TBup_bb^aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[mu][aa] += tensor4in[nu][bb]*transboostlo[nu][mu]*transboostup[bb][aa]; // application on cov con
    }
  }
  else if(tconcovtypeA==TYPEUCOV && tconcovtypeB==TYPEUCOV){
    if(lab2orthofluid==LAB2FF){
      // tfl_{\nu bb} = TBup_\nu^\mu TBup_bb^aa t_{\mu aa}
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[nu][bb] += transboostup[nu][mu]*transboostup[bb][aa]*tensor4in[mu][aa]; // application on cov cov
    }
    else{
      // t_{\mu aa} = t_{\nu bb} TBlo^\nu_\mu TBlo^bb_aa
      DLOOP(mu,nu) DLOOP(aa,bb) tensor4out[mu][aa] += tensor4in[nu][bb]*transboostlo[nu][mu]*transboostlo[bb][aa]; // application on cov cov
    }
  }
  else{
    dualfprintf(fail_file,"No such tconcovtypeA=%d tconcovtypeB=%d\n",tconcovtypeA,tconcovtypeB);
    myexit(934627526);
  }

  return(0);

}




// Use as, e.g.:
//    FTYPE vecvortho[NDIM];
//    concovtype=1; // contravariant input in vec[0-3]
//    vecX2vecVortho(concovtype, vecv, vecvortho);

// converts contravariant (concovtype=1 using tetrcov) or covariant (concovtype=2 using tetrcon) X-based vector into orthonormal V-based vector
// no boost here!
void vecX2vecVortho(int concovtype, struct of_geom *ptrgeom, FTYPE *veclab, FTYPE *vecortho)
{

  ///////////////////////////////
  // get tetrcov and tetrcon
  ///////////////////////////////
  FTYPE tetrconmem[NDIM][NDIM],tetrcovmem[NDIM][NDIM];
  FTYPE (*tetrcon)[NDIM]=tetrconmem;
  FTYPE (*tetrcov)[NDIM]=tetrcovmem;
  get_tetrcovcon(ptrgeom, &tetrcov,&tetrcon); // pass address of pointer since want to give new pointer address if stored


  /////////////////////////////////
  // Now setup transformation
  /////////////////////////////////

  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];

  // vector here is in original X coordinates
  int jj;
  DLOOPA(jj) finalvec[jj]=veclab[jj];


  DLOOPA(jj) tempcomp[jj]=0.0;
  int kk;
  // NOTEMARK: here (unlike in jon_interp_computepreprocess.c) tetrcon and tetrcov convert directly from X coords to V-based orthonormal basis
  if(concovtype==TYPEUCON){
    // transform to orthonormal basis for contravariant vector in X coordinates
	// u^kk[ortho] = tetrcov^kk[ortho]_jj[lab] u^jj[lab]
    DLOOP(jj,kk) tempcomp[kk] += tetrcov[kk][jj]*finalvec[jj];
  }
  else if(concovtype==TYPEUCOV){
    // transform to orthonormal basis for covariant vector in X coordinates
	// u_kk[ortho] = tetrcon_kk[ortho]^jj[lab] u_jj[lab]
    DLOOP(jj,kk) tempcomp[kk] += tetrcon[kk][jj]*finalvec[jj];
  }
  else{
    dualfprintf(fail_file,"No such concovtype=%d\n",concovtype);
    myexit(1);
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];


  // Could now apply lambda that transforms from (e.g.) SPC to Cart using normal orthonormal type transformation for both coordinates to get that transformation

  // final answer:
  DLOOPA(jj) vecortho[jj]=finalvec[jj];


}


