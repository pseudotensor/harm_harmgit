
#include "u2p_defs.h"


#define METHODTYPE  10 
#define NEWT_DIM 5
#define USE_LINE_SEARCH -2

#define NORMMETHOD (1)
#define NEWINVERSION (1)

// why so large?
#define TOL_RESID (1.e-4)

static struct of_geom *ptrlgeom;
static FTYPE U_target[NPR];
static int numnormterms;

static PFTYPE *glpflag; // global pflag for local file
static FTYPE *EOSextra; // global with file scope
static int whicheos; // global with file scope


static FTYPE A_orig[NEWT_DIM][NEWT_DIM];



// static declarations
static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
                                   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
                                                  FTYPE *, FTYPE *, int), 
                                   FTYPE (*res_func) (FTYPE []) );

  
static void func_5d2(FTYPE prguess[], FTYPE dx[], FTYPE resid[],
                     FTYPE (*jac)[NEWT_DIM],
                     FTYPE *f, FTYPE *df, int n);

static FTYPE res_sq_5d( FTYPE x[] );
static int LU_decompose( FTYPE (*A)[NEWT_DIM], int permute[] );
static void LU_substitution( FTYPE (*A)[NEWT_DIM], FTYPE B[], int permute[] );
static int dudp_calc_g(FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE (*Am)[NEWT_DIM] );
static void calc_errx_5d( FTYPE *W, FTYPE *vsq, FTYPE x[] );
static void my_lnsrch(int, FTYPE [], FTYPE, FTYPE [], FTYPE [], FTYPE [], FTYPE *, 
                      FTYPE, FTYPE, int *, FTYPE (*res_func) (FTYPE []));

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) ;







/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

int Utoprim_5d2_final(FTYPE U[NPR], struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats)
{

  FTYPE prim_tmp[NPR], gamma_tmp, qsq_tmp;
  int ret, ret_gam; 
  int pl;
  const int ltrace = 0;
  FTYPE rho0,u;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"Utoprim_5d2_final() called in parallel region\n");
    myexit(9366983);
  }
#endif


  // assign global int pointer to lpflag pointer
  glpflag=lpflag;
  EOSextra=GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  whicheos=WHICHEOS;




#if( WHICHVEL != VELREL4 )
  stderrfprintf("Utoprim_5d2() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1);
#endif


  ptrlgeom=ptrgeom;



  PALLLOOP(pl) U_target[pl] = U[pl];

  /* Set the primitive B-fields (same as conserved B-fields): */
  for (pl = BCON1; pl <= BCON3; pl++)
    prim[pl] = U_target[pl];  


  PALLLOOP(pl) prim_tmp[pl] = prim[pl];


  // Solve the 5D linear system with Newton's method:
  ret = general_newton_raphson(prim_tmp, NPRINVERT, USE_LINE_SEARCH, func_5d2, res_sq_5d);


  *pressure=pressure_rho0_u(whicheos,EOSextra,prim_tmp[RHO],prim_tmp[UU]);



  if( ret ){
    *glpflag= UTOPRIMFAILCONV;
    // then old p should be new p, do nothing since initial guess must be final solution
    if(debugfail>=1) dualfprintf(fail_file,"static solution at t=%21.15g step=%ld i=%d j=%d k=%d wtf=%d\n",t,nstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,debugfail);
    return(ret); 

  }// otherwise got good solution
  else{
    // check densities for positivity
    if((prim_tmp[RHO]<=0.0)||(prim_tmp[UU]<0.0)){
      rho0=prim_tmp[RHO];
      u=prim_tmp[UU];

      //      *glpflag= 5;
      // mostly we hit pr[UU]<0, or maybe never pr[RHO]<0 and pr[UU]<0 alot
      if(debugfail>=2) dualfprintf(fail_file,"utoprim found negative density: %21.15g %21.15g\n",prim_tmp[RHO],prim_tmp[UU]);
      if((rho0<=0.)&&(u>=0.)) *glpflag=  UTOPRIMFAILRHONEG;
      if((rho0>0.)&&(u<0.)) *glpflag= UTOPRIMFAILUNEG;
      if((rho0<=0.)&&(u<0.)) *glpflag= UTOPRIMFAILRHOUNEG;
      if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNNOTADJUSTED) return(0) ; // else let assign -- used to check how bad failure is.
    }

    ret_gam = gamma_calc(prim_tmp, ptrgeom, &gamma_tmp, &qsq_tmp);
    
    if( ret_gam == 1 ) { 
      *glpflag= UTOPRIMFAILCONVW;
      if(debugfail>=2) dualfprintf(fail_file,"utoprim found unphysical velocity ! \n");
      return(0);
    }

  }

  // If we got here, then "everything" is ok : 

  *glpflag= UTOPRIMNOFAIL;

  PALLLOOP(pl) prim[pl] = prim_tmp[pl];

  return( 0 );

}

/*******************************************************************/ 
/*   Auxiliary function required by general_newton_raphson.        */
/*******************************************************************/
static void func_5d2(FTYPE prguess[], FTYPE dxm[], FTYPE resid[], FTYPE (*jac)[NEWT_DIM],
                     FTYPE *f, FTYPE *df, int n)
{
  FTYPE prstart[NPR];
  FTYPE pr[NPR];
  static FTYPE U_curr[NPR];
  
  struct of_state q;

  int i, j, k;
  int pl;
  int failreturn=0;
  FTYPE normtmp;
  FTYPE norm;
  FTYPE residnorm;
  FTYPE d;

  int n_subs = 5;
  int i_subs, do_more_subs ;
  FTYPE  max_subs_change, subs_tmp;
  static int permute[NEWT_DIM], permute_orig[NEWT_DIM];
  static FTYPE Am[NEWT_DIM][NEWT_DIM], Am_0[NEWT_DIM][NEWT_DIM];
  static FTYPE Bm[NEWT_DIM], beta[NEWT_DIM];
  

  const int ltrace = 0;


  // Initialized dgesv() quantities:
  for( i = 0 ; i < NEWT_DIM ; i++ ) { 
    permute[i] = 0 ;
  }


  for( i = 0;  i <  BCON1; i++)   pr[i]  =  prguess[i];
  for( i = BCON1; i <= BCON3; i++)   pr[i]  =  U_target[i] ;

  // Doing this constrained iteration may lead to non-convergence, but then we adjust U if that happens until we converge.

  // store old pr
  PALLLOOP(pl) prstart[pl]=pr[pl];


  // Set the contra/co-variant elements of the 4-vel. and magnetic field tensor :
  //  failreturn=get_state(pr, ptrlgeom, &q);
  get_state(pr, ptrlgeom, &q);
#if(!OPTIMIZED)
  if(failreturn>=1)
    FAILSTATEMENTVOID("func_5d2()", "get_state()",0);
#endif

  // Find the conserved variables associated with the current guess for the primitive variables:
  //  failreturn=primtoU(pr, &q, ptrlgeom, U_curr, NULL);
  primtoU(UNOTHING,pr, &q, ptrlgeom, U_curr, NULL); // returns without geometric factors
#if(!OPTIMIZED)
  if(failreturn>=1)
    FAILSTATEMENTVOID("utoprim.c:usrfun2()", "get_state()",1);
#endif

  // Calculate the jacobian of the forward transformation:
  failreturn=dudp_calc_g(pr, &q, ptrlgeom, Am);
#if(!OPTIMIZED)
  if(failreturn>=1)
    FAILSTATEMENTVOID("utoprim.c:usrfun2()", "dudp_calc()",0);
#endif




  /* normalize error = beta to \rho u^t */
  // more general normalization
#if(NORMMETHOD==0)
  norm=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  norm=0.0;
  numnormterms=0;
  for (j = 0; j < NPRINVERT; j++){
    for (k = 0; k < NPRINVERT; k++){
      if(fabs(Am[j][k]) > NUMEPSILON){
        norm+=fabs( Am[j][k] );
        numnormterms++;
      }
    }
  }
  norm=((FTYPE)(numnormterms))/(norm); // (i.e. inverse of average)
#elif(NORMMETHOD==2)
  norm = 1.0;
#endif

  
  for (j = 0; j < NPRINVERT; j++) {
    for (k = 0; k < NPRINVERT; k++) {
      jac[j][k] = Am[j][k];
      Am[j][k] *= (norm);
    }
  }
  
    
  // determine normalized error
  *df = residnorm = 0.;
  for (k = 0; k < NPRINVERT; k++) {
    resid[k] =  U_curr[k] - U_target[k] ;
    beta[k] = -resid[k] *(norm);
    residnorm = 0.5*(fabs(U_curr[k]) + fabs(U_target[k]));
    residnorm = (residnorm == 0.) ? 1. : 1./residnorm;
    *df -= residnorm * resid[k]*resid[k];
  }

  *f = -0.5*(*df);

  //-scn Need to return with normalization factor that normalizes U[] to 1 :
  normtmp = 0.;
  for (k = 0; k < NPRINVERT; k++)
    normtmp += fabs(U_target[k]) ;

  normtmp *= norm;
  norm = normtmp;


  /*******************************************/
  /* Solve the linear system routine:        */
  /*******************************************/

  /******************************************************************************************/
  /**** With simple LU-decomp. and backward/forward subst. */
#if( NEWINVERSION ) 

  for( i = 0 ; i < n ; i++ ) { 
    for( j = 0 ; j < n ; j++ ) { 
      A_orig[i][j] = Am[i][j];
    }
    dxm[i] = beta[i];
  }


  // Get the LU matrix:
  if( LU_decompose( Am,  permute ) != 0  ) { 
    dualfprintf(fail_file, "func_5d2(): singular matrix encountered! \n");
    
  }

  // Solve the linear system: 
  LU_substitution( Am,  dxm, permute );


  /******************************************************************************************/
  /**** With simple LU-decomp. and iterative subst., to minimize erros in subst. procedure */
#else

  for( i = 0 ; i < n ; i++ ) { 
    for( j = 0 ; j < n ; j++ ) { 
      Am_0[i][j] = Am[i][j];
    }
    dxm[i] = beta[i];
  }

  // Get the LU matrix:
  if( LU_decompose( Am,  permute ) != 0  ) { 
    dualfprintf(fail_file, "func_5d2(): singular matrix encountered! \n");
    
  }

  for( i = 0 ; i < n ; i++ ) {
    permute_orig[i] = permute[i];
  }

  // Solve the linear system: 
  LU_substitution( Am,  dxm, permute );

  do_more_subs = 1;
  i_subs = 0;

  while( do_more_subs == 1 ) { 
    i_subs++ ; 

    // Calculate the \delta B: 
    for( i = 0 ; i < n ; i++ ) { 
      Bm[i] = -beta[i];
      for( j = 0 ; j < n ; j++ ) { 
        Bm[i] += Am_0[i][j] * dxm[j] ; 
      }
    }

    // Find the correction to dx[] :
    LU_substitution( Am,  Bm, permute_orig );

    // Add correction to dx[]:
    max_subs_change = 0.;

    dualfprintf(fail_file,"i_subs = %d: ddx/dx[0-4] : ", i_subs);
    for( i = 0 ; i < n ; i++ ) { 
      dxm[i] -= Bm[i];
      subs_tmp = fabs(Bm[i]/dxm[i]);
      if( subs_tmp > max_subs_change ) { 
        max_subs_change = subs_tmp ;
      }
      dualfprintf(fail_file,"%21.15g ", subs_tmp );
    }
    dualfprintf(fail_file,"\n ");
    //    

    if( (i_subs >= n_subs) || (max_subs_change < 1.e-10 ) ) { 
      do_more_subs = 0 ; 
    }
  }

#endif


#if(!OPTIMIZED)
  if( ltrace ) { 
    dualfprintf(fail_file,"usrfun2(): dx[0-4] = ");
    for (i = 1; i <= n; i++) {
      dualfprintf(fail_file," %10.6e ,", dxm[i-1]);
    }
    dualfprintf(fail_file,"\n");
    dualfprintf(fail_file,"usrfun2(): x[0-4]  = ");
    for (i = 1; i <= n; i++) {
      dualfprintf(fail_file," %10.6e ,", prguess[i-1]);
    }
    dualfprintf(fail_file,"\n");
    
  }
#endif


}



static void calc_errx_5d( FTYPE *W, FTYPE *vsq, FTYPE x[] ) 
{
  int id, jd;
  FTYPE gsqtmp;


  *vsq = 0. ;
  for(id=1;id<4;id++)
    for(jd=1;jd<4;jd++) 
      *vsq += ptrlgeom->gcov[GIND(id,jd)]*x[UTCON1+id-1]*x[UTCON1+jd-1] ;      
  
  gsqtmp = 1. + *vsq;
  *vsq /= gsqtmp;

  *W = gsqtmp * ( gam*x[UU] + x[RHO] ) ;
  
}


/***************************************************************/
/* Just calculate the scalar residual for my_lnsrch():  */
/***************************************************************/
static FTYPE res_sq_5d( FTYPE x[] ) 
{
  int i,k;
  FTYPE pr[NPR];
  static FTYPE U_curr[NPR];
  struct of_state q;
  int failreturn=0;
  FTYPE retval, normret; 


  for( i = 0;  i <  BCON1; i++)   pr[i]  =  x[i];
  for( i = BCON1; i <= BCON3; i++)   pr[i]  =  U_target[i] ;


  // Set the contra/co-variant elements of the 4-vel. and magnetic field tensor :
  //  failreturn=get_state(pr, ptrlgeom, &q);
  get_state(pr, ptrlgeom, &q);
#if(!OPTIMIZED)
  if(failreturn>=1)
    FAILSTATEMENT("res_sq_5d()", "get_state()", 1);
#endif

  // Find the conserved variables associated with the current guess for the primitive variables:
  //  failreturn=primtoU(pr, &q, ptrlgeom, U_curr, NULL);
  primtoU(UNOTHING,pr, &q, ptrlgeom, U_curr, NULL);
#if(!OPTIMIZED)
  if(failreturn>=1)
    FAILSTATEMENT("res_sq_5d()", "primtoU()", 1);
#endif


  retval = 0.0;

  for (i = 0; i < NPRINVERT; i++) {
    normret = 0.5*( fabs(U_curr[i]) + fabs(U_target[i]) ) ;
    normret = ( normret == 0. )  ?  1.  : 1./normret;
    retval += normret * (U_curr[i] - U_target[i])*(U_curr[i] - U_target[i]) ;
  }

  return( 0.5*retval );
  
}



// wrapper for 5D2 method to use dudp_calc_gen in dudp_calc.c
static int dudp_calc_g(FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE (*Am)[NEWT_DIM] )
{
  int j,k;
  static int firstc=1;
  static FTYPE **alpha;
 
  // get memory for matrix as setup for dudp_calc()
  if (firstc) {
    firstc = 0;
    alpha = dmatrix(1, NEWT_DIM, 1, NEWT_DIM);
  }

  // get matrix
  dudp_calc(WHICHEOS,EVOLVENOENTROPY,EOSextra,pr,q,geom,alpha);

  // copy offset matrix
  for(j=1;j<=5;j++) for(k=1;k<=5;k++) Am[j-1][k-1]=alpha[j][k];


  return(0) ;
}



/*************************************************************************/
/*************************************************************************/
/* Performs a LU decomposition of the matrix A using Crout's method      */
/* with partial implicit pivoting.  The exact LU decomposition of the    */
/* matrix can be reconstructed from the resultant row-permuted form via  */
/* the integer array permute[]                                           */ 
/*                                                                       */
/* The algorithm closely follows that in ludcmp.c of "Numerical Recipes  */
/* in C" by Press et al. 1992.                                           */
/*                                                                       */
/* This will be used to solve the linear system  A.x = B                 */
/*                                                                       */
/* Returns (1) if a singular matrix is found,  (0) otherwise.            */
/*************************************************************************/


static int LU_decompose( FTYPE (*A)[NEWT_DIM], int permute[] )
{

  const  FTYPE absmin = 1e-30; /* Value to be used instead of 0 for singular matrices */

  static FTYPE row_norm[NEWT_DIM];
  FTYPE  absmax, maxtemp, mintemp;
  

  int i, j, k, max_row;
  int n = NEWT_DIM;



  max_row = 0;


  /* Find the maximum elements per row so that we can pretend later
     we have unit-normalized each equation: */

  for( i = 0; i < n; i++ ) { 
    absmax = 0.;
    
    for( j = 0; j < n ; j++ ) { 
      
      maxtemp = fabs( A[i][j] ); 

      if( maxtemp > absmax ) { 
        absmax = maxtemp; 
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if( absmax == 0. ) { 
      dualfprintf(fail_file, "LU_decompose(): row-wise singular matrix\n");
      //      for(k=0;k<n;k++) { 
      // dualfprintf(fail_file,"A[%d][0-4] = %16.10e %16.10e %16.10e %16.10e %16.10e \n", k, A[k][0],A[k][1],A[k][2],A[k][3],A[k][4]);
      //      }
      return(1);
    }

    row_norm[i] = 1. / absmax ;   /* Set the row's normalization factor. */
  }


  /* The following the calculates the matrix composed of the sum 
     of the lower (L) tridagonal matrix and the upper (U) tridagonal
     matrix that, when multiplied, form the original maxtrix.  
     This is what we call the LU decomposition of the maxtrix. 
     It does this by a recursive procedure, starting from the 
     upper-left, proceding down the column, and then to the next
     column to the right.  The decomposition can be done in place 
     since element {i,j} require only those elements with {<=i,<=j} 
     which have already been computed.
     See pg. 43-46 of "Num. Rec." for a more thorough description. 
  */

  /* For each of the columns, starting from the left ... */
  for( j = 0; j < n; j++ ) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ ) { 
        A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for( i = j; i < n; i++ ) {

      //      if( i==(n-1) && j==(n-1) ) dualfprintf(fail_file,"last aorig = %21.15g \n",A[i][j] );

      for (k = 0; k < j; k++) { 
        A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit 
         unit-normalization (represented by row_norm[i]) of each row: 
      */
      maxtemp = fabs(A[i][j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
        absmax = maxtemp;
        max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
       This is the partial pivoting procedure that ensures we don't divide
       by 0 (or a small number) when we solve the linear system.  
       Also, since the procedure starts from left-right/top-bottom, 
       the pivot values are chosen from a pool involving all the elements 
       of column_j  in rows beneath row_j.  This ensures that 
       a row  is not permuted twice, which would mess things up. 
    */
    if( max_row != j ) {

      /* Don't swap if it will send a 0 to the last diagonal position. 
         Note that the last column cannot pivot with any other row, 
         so this is the last chance to ensure that the last two 
         columns have non-zero diagonal elements.
      */

      if( (j == (n-2)) && (A[j][j+1] == 0.) ) {
        max_row = j;
      }
      else { 
        for( k = 0; k < n; k++ ) { 

          maxtemp       = A[   j   ][k] ; 
          A[   j   ][k] = A[max_row][k] ;
          A[max_row][k] = maxtemp; 

        }

        /* Don't forget to swap the normalization factors, too... 
           but we don't need the jth element any longer since we 
           only look at rows beneath j from here on out. 
        */
        row_norm[max_row] = row_norm[j] ; 
      }
    }

    /* Set the permutation record s.t. the j^th element equals the 
       index of the row swapped with the j^th row.  Note that since 
       this is being done in successive columns, the permutation
       vector records the successive permutations and therefore
       index of permute[] also indexes the chronology of the 
       permutations.  E.g. permute[2] = {2,1} is an identity 
       permutation, which cannot happen here though. 
    */

    permute[j] = max_row;

    if( A[j][j] == 0. ) { 
      A[j][j] = absmin;
      //dualfprintf(fail_file, "LU_decompose(): column-wise singular matrix\n");
      //return(1);
    }


    /* Normalize the columns of the Lower tridiagonal part by their respective 
       diagonal element.  This is not done in the Upper part because the 
       Lower part's diagonal elements were set to 1, which can be done w/o 
       any loss of generality.
    */
    if( j != (n-1) ) { 
      maxtemp = 1. / A[j][j]  ;
      
      for( i = (j+1) ; i < n; i++ ) {
        A[i][j] *= maxtemp;
      }
    }

  }

  return(0);

  /* End of LU_decompose() */

}


/************************************************************************/
/************************************************************************/
/*                                                                      */
/*   Performs the forward (w/ the Lower) and backward (w/ the Upper)    */
/*   substitutions using the LU-decomposed matrix (*A)[] of the original */
/*   matrix A' of the linear equation:  A'.x = B.  Upon entry, (*A)[]    */
/*   is the LU matrix, B[] is the source vector, and permute[] is the  */
/*   array containing order of permutations taken to the rows of the LU */
/*   matrix.  See LU_decompose() for further details.    */
/*         */
/*   Upon exit, B[] contains the solution x[], (*A)[] is left unchanged. */ 
/*         */
/************************************************************************/

static void LU_substitution( FTYPE (*A)[NEWT_DIM], FTYPE B[], int permute[] )
{

  int i, j ;

  int n = NEWT_DIM;

  FTYPE tmpvar,tmpvar2;

  
  //  for( i = 0 ; i < n ; i++ ) { 
  //    dualfprintf(fail_file, "A[%d][0-4] = %21.15g %21.15g %21.15g %21.15g %21.15g \n", i, 
  //  A[i][0], A[i][1], A[i][2], A[i][3], A[i][4]);
  //  }
  //  dualfprintf(fail_file, "B0[0-4] = %21.15g %21.15g %21.15g %21.15g %21.15g \n",  B[0], B[1], B[2], B[3], B[4]);
  //  dualfprintf(fail_file, "permm   = %d %d %d %d %d \n",  permute[0], permute[1], permute[2], permute[3], permute[4]);
  //

  /* Perform the forward substitution using the LU matrix. 
   */
  for(i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the 
       B vector to match the permutation of the LU matrix. 
       Since only the rows above the currrent one matter for 
       this row, we can permute one at a time. 
    */

    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    //    dualfprintf(fail_file, "B0[%d] , B0[%d] = %21.15g %21.15g %21.15g \n", i, permute[i], B[i],B[permute[i]],tmpvar );

    //    for( j = 0; j < i ; j++ ) { 
    for( j = (i-1); j >= 0 ; j-- ) { 
      //      dualfprintf(fail_file, "tmpsum %d %d : before, A , B, sum = %21.15g %21.15g %21.15g ", i, j, tmpvar, A[i][j],B[j] );
      tmpvar -=  A[i][j] * B[j];
      //      dualfprintf(fail_file, " %21.15g \n", tmpvar );
    }

    B[i] = tmpvar; 

    //    dualfprintf(fail_file, "B1[%d] = %21.15g \n", i, B[i]);

  }
    

  /* Perform the backward substitution using the LU matrix. 
   */
  for( i = (n-1); i >= 0; i-- ) { 
    
    for( j = (i+1); j < n ; j++ ) { 
      
      B[i] -=  A[i][j] * B[j];

    }
    
    B[i] /= A[i][i] ; 
    
    //    dualfprintf(fail_file, "B2[%d] = %21.15g \n", i, B[i]);
  }

  /* End of LU_substitution() */

}


/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, int do_line_search,
                                   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), 
                                   FTYPE (*res_func) (FTYPE []) )
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, errx_old, errx_oldest, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE randtmp, tmp;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE resid_norm, resid_check, grad_check;
  FTYPE res_func_val, res_func_old, res_func_new;
  FTYPE dn[NEWT_DIM], del_f[NEWT_DIM];

  int   keep_iterating, i_increase, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 1;


  retval = 0;


  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;


  vsq_old = vsq = W = W_old = 0.;


  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


#if(!OPTIMIZED)
    /*  Check for bad untrapped divergences : */
    if( (finite(f)==0) || (finite(df)==0) ) {
      dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
      return(1);
    }
#endif


#if(!OPTIMIZED)
    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      // for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }
#endif

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    if( do_line_search == 1 ) { 

      /* Compare the residual to its initial value */
      if( n_iter == 0 ) { 
        resid_norm = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_norm += fabs(resid[id]);
        }
        resid_norm /= 1.0*n ;
        if( resid_norm == 0.0 ) resid_norm = 1.0;
      }

      for( id = 0; id < n ; id++) {
        tmp = 0.;
        for( jd = 0; jd < n ; jd++) {
          tmp += jac[jd][id] * resid[jd];
        }
        del_f[id] = tmp;
      }
      for( id = 0; id < n ; id++) {
        dn[id] = dx[id];
      }

      my_lnsrch(n, x_old-1, f_old, del_f-1, dn-1, x-1, &f, TOL_LINE_STEP, SCALEMAX, &retval, res_func);

      /* dx is needed for errx calculation below: */
      for( id = 0; id < n ; id++) {
        dx[id] = x[id] - x_old[id];
      }

#if(!OPTIMIZED)
      if( ltrace ) { 
        res_func_val = res_func(x);
        res_func_old = res_func(x_old);
        dualfprintf(fail_file,"gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
                    f_old, f, res_func_old, res_func_val );
        dualfprintf(fail_file,"gnr(): x_old = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",x_old[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): x     = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",x[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): dn    = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",dn[id]);
        }
        dualfprintf(fail_file,"\n ");
        dualfprintf(fail_file,"gnr(): del_f = ");
        for( id = 0; id < n ; id++) {
          dualfprintf(fail_file," %21.15g ",del_f[id]);
        }
        dualfprintf(fail_file,"\n ");
      }
#endif

      /* Check to see if line search problem is because the residual vector is already small enough */
      if( retval == 1 ) {
        resid_check = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_check += fabs(resid[id]);
        }
        resid_check /= 1.0*n;
 
        if( resid_check <= resid_norm * NEWT_FUNC_TOL ) {
          retval = 0;
        }
        if( ltrace && retval ) { 
          dualfprintf(fail_file,"general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
   
        }   
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH == do_line_search) ) { 
#if(!OPTIMIZED)
        if( ltrace ) { 
          dualfprintf(fail_file,"gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
          dualfprintf(fail_file,"gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
        }
#endif       
        retval = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
        for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
        return( retval );
      }

      /* Check to see if it is trapped in a local minimum, i.e. gradient is too small */ 
      if( retval == 1 ) { 
        grad_check = 0.0e0;
        for( id = 0; id < n ; id++) {
          resid_check = (x[id] == 0.) ? 1.0 : fabs(x[id]) ;
          grad_check  +=  del_f[id] * resid_check ;
        }
        resid_check = (f == 0.) ? 1.0 : fabs(f) ;
        grad_check /= resid_check;
 
        /* Then we've most likely found a solution: */
        if( grad_check > GRADMIN ) { 
          retval = -1;
        }
        else if( ltrace ) { 
          dualfprintf(fail_file,"general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval,grad_check);
        }
      }
    }
    else {
      /* don't use line search : */
      for( id = 0; id < n ; id++) {
        x[id] += dx[id]  ;
      }

    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
#if( NEWCONVERGE == 1 )
    W_old = W;
    vsq_old = vsq;
    calc_errx_5d(&W, &vsq, x);
    errx  = (W==0.) ?  fabs(W-W_old) : fabs((W-W_old)/W);

    /* For the old criterion, look at errors in each indep. variable(s) (except for 5D) : */
#else
    W_old = W;
    vsq_old = vsq;
    calc_errx_5d(&W, &vsq, x);
    errx  += (vsq==0.) ?  fabs(vsq-vsq_old) : fabs((vsq-vsq_old)/vsq);
    errx  += (W==0.) ?  fabs(W-W_old) : fabs((W-W_old)/W);
    errx  /= 2.;
#endif


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific:

    //    x[0] = fabs(x[0]);
    //    x[1] = fabs(x[1]);

    /****************************************/
    /* Check to see if we're in a infinite loop with error function: */
    /****************************************/
#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL) )  errx = -errx;
#endif 

    /****************************************/
    /* If there's a problem with line search, then stop iterating: */
    /****************************************/
    if( (retval == 1) || (retval == -1) ) errx = -errx;


#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file," general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
                  n_iter,f_old,f,errx_old,errx );
      dualfprintf(fail_file,"gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",x_old[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): x     = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",x[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): dx     = ");
      for( id = 0; id < n ; id++) {
        dualfprintf(fail_file," %21.15g ",dx[id]);
      }
      dualfprintf(fail_file,"\n ");
      
    }
#endif

    /****************************************/
    /* Prepare for the next iteration, set the "old" variables: */
    /****************************************/
    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    /****************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations before stopping */
    /****************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (f <= TOL_RESID) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) || (finite(x[0])==0) || (finite(x[1])==0) 
      || (finite(x[2])==0) || (finite(x[3])==0) || (finite(x[4])==0) ) { 
#if(!OPTIMIZED)
    dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
    dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
#endif
    return(1);
  }


  if( (fabs(errx) > MIN_NEWT_TOL) || (f > TOL_RESID)  ){
    if( (do_line_search != USE_LINE_SEARCH) || (USE_LINE_SEARCH < 0) ) { 
      //bin_newt_data( errx, n_iter, 0, 0 );
#if(!OPTIMIZED)
      if(ltrace2) {
        dualfprintf(fail_file," totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
 
      }
      if(ltrace) {
        dualfprintf(fail_file,"general_newton_raphson():  did not find solution \n");
        if( retval == -1 ) {
          dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x = ");
          for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x[id]);
          dualfprintf(fail_file,"\n");
          dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x_old = ");
          for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x_old[id]);
          dualfprintf(fail_file,"\n");
        }
 
      }
      // dualfprintf(fail_file,"gnr retval2 = %4i \n", 1); 
#endif
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
      //      dualfprintf(fail_file,"gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n", do_line_search, retval, errx );
      //      
      retval2 = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      //      dualfprintf(fail_file,"gnr retval3 = %4i \n", retval2); 
      return( retval2 );
    }
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    //bin_newt_data( errx, n_iter, 1, 0 );
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      
    }
    if(ltrace) {
      dualfprintf(fail_file,"general_newton_raphson(): found minimal solution \n");
      
    }
    //    dualfprintf(fail_file,"gnr retval4 = %4i \n", 0); 
#endif
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    //bin_newt_data( errx, n_iter, 2, 0 );
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      
    }
    //    dualfprintf(fail_file,"gnr retval5 = %4i \n", 0); 
#endif
    return(0);
  }

#if(!OPTIMIZED)
  dualfprintf(fail_file,"gnr retval6 = %4i \n", 0);
#endif
  return(0);

}

/**************************************************** 
*****************************************************/

#define ALF 1.0e-4

static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
                       FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) )
{
  int i;
  FTYPE a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  int icount=0;
  int bad_step;
  FTYPE bad_step_factor = 2.0;
  
  const int ltrace = 0;
 
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    //    temp=fabs(p[i])/MYMAX(fabs(xold[i]),1.0);
    temp= (xold[i] == 0.) ? fabs(p[i]) :  fabs(p[i]/xold[i]);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"my_lnsrch(): sum, slope, test, alamin =   %21.15g  %21.15g  %21.15g  %21.15g \n",sum,slope,test,alamin); 
  }
#endif

  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

    *f=(*func)(x+1);
      
    bad_step = 0;

    if( finite(*f)==0 ) { 
      bad_step = 1;
    }

    //      if( bad_step ) alam /= bad_step_factor;
    //      if (alam < alamin) bad_step = 0;

    if( bad_step ) { 
      *check = 2;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"my_lnsrch(): bad_step = 1,  f = %21.15g \n", *f); 
#endif
      return;
    }
      
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
#if(!OPTIMIZED)
      if( ltrace ) { 
        dualfprintf(fail_file,"my_lnsrch(): alam < alamin: alam, alamin = %21.15g  %21.15g \n", alam,alamin); 
      }
#endif
      return;
    } 
    else if (*f <= fold+ALF*alam*slope) {
#if(!OPTIMIZED)
      if( ltrace ) { 
        dualfprintf(fail_file,"my_lnsrch(): good exit:  alam, alamin, f, fold = %21.15g  %21.15g %21.15g  %21.15g \n", alam,alamin, *f, fold); 
      }
#endif
      return;
    }
    else {
      if (alam == 1.0) {
        tmplam = -slope/(2.0*(*f-fold-slope));
#if(!OPTIMIZED)
        if( ltrace ) {
          dualfprintf(fail_file,"my_lnsrch(): setting tmplam!!    tmplam, alam =  %21.15g  %21.15g !!\n", tmplam, alam);
        }
#endif
      }
      else {
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold2-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0) tmplam = -slope/(2.0*b);
        else {
          disc=b*b-3.0*a*slope;
          if (disc<0.0) {
#if(!OPTIMIZED)
            if( disc < -1.e-10 ) { // GODMARK -- needs normalization
              dualfprintf(fail_file,"my_lnsrch(): Big Roundoff problem:  disc = %21.15g \n", disc);
            }
#endif
            disc = 0.;
          }
          else tmplam=(-b+sqrt(disc))/(3.0*a);
        }
        if (tmplam>0.5*alam)
          tmplam=0.5*alam;
#if(!OPTIMIZED)
        if( ltrace ) {
          dualfprintf(fail_file,"my_lnsrch(): rhs1, rhs2, a, b, tmplam, alam =  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g !!\n",
                      rhs1, rhs2, a, b, tmplam, alam );
        }
#endif
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=MYMAX(tmplam,0.1*alam);
#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"my_lnsrch(): icount, alam, alam2, tmplam  =   %4i  %21.15g  %21.15g  %21.15g \n",
                  icount, alam, alam2, tmplam);
    }
#endif
    icount++;
  }
}
#undef ALF


#define N_CONV_TYPES 3
#define N_NITER_BINS (MAX_NEWT_ITER + 3)
#define NBINS 200

/**************************************************** 

  Used to gather statistics for the Newton solver during a disk run:

  -- bins are inclusive on their minimum bounds, and exlusive on their max.


      lerrx_[min,max] = the min./max. bounds of the errx historgram
      d_errx          = errx histogram bin width
      n_beyond_range  = number of solutions with errx max. beyond range speciified by lerrx_[min,max]
      n_conv_types    = number of types of ways newton procedure can exit;

      print_now       = 1 if you want to just print the histograms without counting a solution;

*****************************************************/

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) 
{

  /* General variables */
  int i, j;
  static int first_call   = 1;
  static long int n_bin_calls  = 0L;
  const  long int n_output_freq = 128*128*100;
  FTYPE lerrx;

  /* Variables for the errx histogram */ 
  static const FTYPE lerrx_min = -23.;
  static const FTYPE lerrx_max =   4.;
  static FTYPE d_errx;
  static long int n_errx[N_CONV_TYPES][NBINS];
  static FTYPE xbin[NBINS];
  static long int n_beyond_range = 0L;   /* Number of points that lie out of the bounds of our errx histogram */
  int ibin;

  /* Variables for the histogram of the number of newton iterations : */
  static long int n_niters[N_CONV_TYPES][N_NITER_BINS];
  

  /* Clear arrays, set constants : */
  if( first_call ) {
    d_errx    = ((lerrx_max - lerrx_min)/(1.*NBINS));
    
    for( i = 0; i < N_CONV_TYPES; i++ ) { 
      for( j = 0; j < NBINS; j++ ) { 
        n_errx[i][j] = 0;
      }
      for( j = 0; j < N_NITER_BINS; j++ ) { 
        n_niters[i][j] = 0;
      }
    }

    for( j = 0; j < NBINS; j++ ) { 
      xbin[j] = lerrx_min + (1.*j + 0.5)*d_errx;
    }
      
    first_call = 0;
  }


  if( print_now != 1 ) {

    /* Check validity of arguments : */
    errx = fabs(errx) ;
    lerrx = log10(errx + 1.0e-2*pow(10.,lerrx_min));

    if( (niters < 0) || (niters >= N_NITER_BINS) ) {
      dualfprintf(fail_file,"bin_newt_data(): bad value for niters = %d \n", niters );
      fflush(stdout);
      return;
    }
    
    /* Determine ibin */
    if( lerrx < lerrx_min ) {
      ibin = 0;
      n_beyond_range++ ;
    }
    else if( lerrx >= lerrx_max ) {
      ibin = NBINS - 1;
      n_beyond_range++ ;
    }
    else {
      ibin = (int) ( ( lerrx - lerrx_min ) / d_errx );
    }
    
    /* Tally this solution  */
    n_errx[ conv_type][ibin]++;
    n_niters[conv_type][niters]++;

  }


  /* Print out the histograms periodically or when asked to : */
  if( print_now  ||  ( (n_bin_calls % n_output_freq) == 0 ) ) {

    dualfprintf(log_file,"t = %21.15g ,  n_beyond_range = %ld , n_bin_calls = %ld \n", t, n_beyond_range, n_bin_calls);

    /* ERRX */
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");
    dualfprintf(log_file,"                         x");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < NBINS; i++ ) { 
      dualfprintf(log_file,"%21.15g ",xbin[i]);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        dualfprintf(log_file,"%13ld ", n_errx[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");


    /* NITER */

    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");
    dualfprintf(log_file,"        niter");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < N_NITER_BINS; i++ ) { 
      dualfprintf(log_file,"%13d ", i);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        dualfprintf(log_file,"%13ld ", n_niters[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");

  }
    
    
  if( print_now != 1 )  n_bin_calls++;


  return;

}
  
#undef N_CONV_TYPES 
#undef N_NITER_BINS 
#undef NBINS 

