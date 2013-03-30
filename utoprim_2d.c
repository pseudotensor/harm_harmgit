
#include "decs.h"
#include "utoprim_2d.h" 

/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
static FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;

static int    dodamp          = 0;          // whether to use damped method                                           
static FTYPE dampfactorchange= 0.8;                                                                                  
static FTYPE dampfactor_min  = 1.0e-7;                                                                               
static int    preiterdamp     = 5;                                                                                    
static int    numdampedtot    =5;                                                                                     
static int    numstabletot    =5;                                                                                     
static FTYPE alpha_newt      = 0.;                                                                                   

static FTYPE wglobal;
static PFTYPE *glpflag; // global pflag for local file

/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

int Utoprim_2d(FTYPE U[NPR], struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats)
{

  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, k;
  PFTYPE ret; 
  FTYPE gcov[SYMMATRIXNDIM], gcon[SYMMATRIXNDIM], alpha;
  const int ltrace = 0;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"Utoprim_2d() called in parallel region\n");
    myexit(9366983);
  }
#endif

  newtonstats->lntries=0;
 
  // assign global int pointer to lpflag pointer
  glpflag=lpflag;



#if( WHICHVEL != VELREL4 )
  stderrfprintf("Utoprim_2d() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1);
#endif


  /* Set the geometry variables: */
  for(i = 0; i < NDIM; i++ ) {
    for(j = 0; j < NDIM; j++ ) {
      gcov[GIND(i,j)] = ptrgeom->gcov[GIND(i,j)];
      gcon[GIND(i,j)] = ptrgeom->gcon[GIND(i,j)];
    }
  }
  alpha = 1.0/sqrt(-gcon[GIND(0,0)]);
  

  /* Calculate the transform the CONSERVATIVE variables into the new system */
  // notice -alpha
  U_tmp[RHO] = U[RHO] * alpha;
  for( i = UU; i <= U3; i++ ) {
    U_tmp[i] = -alpha * U[i];
  }
  for( i = B1; i <= B3; i++ ) {
    U_tmp[i] = alpha * U[i] ;
  }


  /* Calculate the transform the PRIMITIVE variables into the new system */
  PALLLOOP(i) prim_tmp[i] = prim[i];

  for( i = B1; i <= B3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, ptrgeom, prim_tmp, pressure);

  /* Check conservative variable transformation: */
  if( ltrace ) {
    primtoU_g(ptrgeom, prim_tmp, gcov, gcon, U_tmp2 ); 
    for( i = 0; i < NPR; i++ ) {
      fprintf( stdout, "Utoprim_new1(): Utmp1[%d] = %21.15g , Utmp2[%d] = %21.15g , dUtmp[%d] = %21.15g \n", 
               i, U_tmp[i], i, U_tmp2[i], i, fabs( (U_tmp[i]-U_tmp2[i]) / ( (U_tmp2[i]!=0.) ? U_tmp2[i] : 1. ) )  ); 
    }
  }

  /* Transform new primitive variables back : */ 
  for( i = 0; i < B1; i++ ) {
    prim[i] = prim_tmp[i];
  }
  for( i = B1; i <= B3; i++ ) {
    prim[i] = prim_tmp[i] / alpha;
  }


  return( ret ) ;

}

/**********************************************************************************

attempt an inversion from U to prim using the initial guess
prim.

  -- assumes that 
             /  rho gamma        \
         U = |  - alpha T^t_\mu  |
             \  alpha B^i        /



             /    rho        \
  P = |    uu         |
             | \tilde{u}^i   |
             \  alpha B^i   /


             *glpflag:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence (occurrence of "nan" or "+/-inf" ;
      
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
      3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: utsq<0 or utsq > UTSQ_TOO_BIG   with new  W;
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR], FTYPE *pressure)
{
  FTYPE Wtest;
  //  extern  void func_1d_orig(FTYPE x[1], FTYPE dx[1], FTYPE *f, FTYPE *df, int n);
  //  extern FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;
  FTYPE x_1d[1];
  FTYPE gcov[SYMMATRIXNDIM], gcon[SYMMATRIXNDIM];
  FTYPE QdotB,Bcon[4],Bcov[4],Qcov[4],Qcon[4],ncov[4],ncon[4],Qsq,Qtcon[4];
  FTYPE rho0,u,p,w,gammasq,gamma,W_last,W,utsq ;
  int i,j, ret, errval ;
  FTYPE Ui[8],primi[8];

  const int ltrace = 0;
  const int ltrace2 = 0;


  if( ltrace ) {
    for(i=0;i<8;i++) fprintf(stdout,"%d %21.15g %21.15g\n",i,prim[i],U[i]) ;
    fflush(stdout);
  }

  // Assume ok initially:
  *glpflag= UTOPRIMNOFAIL;

  for(i = 0; i < NDIM; i++ ) {
    for(j = 0; j < NDIM; j++ ) {
      gcov[GIND(i,j)] = ptrgeom->gcov[GIND(i,j)];
      gcon[GIND(i,j)] = ptrgeom->gcon[GIND(i,j)];
    }
  }

  /* save guess and prims */
  for(i=0;i<8;i++) Ui[i] = U[i] ;
  for(i=0;i<8;i++) primi[i] = prim[i] ;

  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* Set parameters based upoon the conservative variables : */
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;
  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  QdotBsq = QdotB*QdotB ;

  ncov_calc(gcon,ncov) ;
  raise_g(ncov,gcon,ncon);

  Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[GIND(i,j)]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG ) {
    if( debugfail>=2 ) dualfprintf(fail_file,"Utoprim_new1(): utsq < 0 in utoprim_1d attempt, utsq = %21.15g \n", utsq) ;
    *glpflag= UTOPRIMFAILCONVGUESSUTSQ;
    return(0) ;
  }

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u_2d(rho0,u) ;

  *pressure=p;

  w = rho0 + u + p ;
  gammasq = 1. + utsq ;
  W_last = w*gammasq ;

  if( ltrace ) {
    fprintf(stdout,"u = %21.15g,  p = %21.15g,  Bsq = %21.15g, Qsq = %21.15g \n",prim[RHO],prim[UU],Bsq,Qsq);
    fprintf(stdout,"Bcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
    fprintf(stdout,"Bcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcov[0],Bcov[1],Bcov[2],Bcov[3]);
    fprintf(stdout,"Qcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcon[0],Qcon[1],Qcon[2],Qcon[3]);
    fprintf(stdout,"Qcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcov[0],Qcov[1],Qcov[2],Qcov[3]);
    fprintf(stdout,"call find_root\n") ;  fflush(stdout);
    fflush(stdout);
  }

  W = W_last;
  wglobal=w;
  ret = find_root_2D_gen( &W )  ; 

  if( ltrace ) {
    stderrfprintf("(W,W_last,Bsq,Qtsq,QdotB,gammasq,Qdotn) %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
                  W,W_last,
                  Bsq,Qtsq,QdotB,gammasq,Qdotn) ;
  }

  if( ltrace2 ) {
    fprintf(stdout, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    fflush(stdout);
  }

  
  /* Problem with solver, so return denoting error before doing anything further */
  if( (ret != 0) || (W == FAIL_VAL) ) {
    if( debugfail>=2 ) dualfprintf(fail_file, "Failed to find a prim. var. solution!! %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g \n",W_last,Bsq,QdotBsq,Qdotn,D,Qtsq);
    *glpflag= ret+UTOPRIMFAILCONVRET+1; // related to UTOPRIMFAILCONVRET
    return(0);
  }
  else{
    
    Wtest=W/wglobal;
    if(Wtest <= 0. || Wtest > W_TOO_BIG) {
      //if(Wtest <= 0.) {
      if( debugfail>=2 ) dualfprintf(fail_file,"Wtest failure %21.15g \n",Wtest) ;
      *glpflag= UTOPRIMFAILCONVW;
      return(0) ;
    }
    
    // below test for UTSQ_TOO_BIG is good normalized version of the above unnormalized W test.
    utsq = utsq_calc(W) ;
    if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
      utsq = fabs(utsq);
    }
    if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
      if( debugfail>=2 ) dualfprintf(fail_file,"utsq failure:  utsq = %21.15g , W = %21.15g \n",utsq, W) ;
      *glpflag= UTOPRIMFAILCONVUTSQ;
      return(0) ;
    }
  }

  if( ltrace ) {
    fprintf(stdout,"done find_root\n") ; fflush(stdout);
  }


  /* Past all checks, so now use good result to find final values of all primitive variables */

  gamma = sqrt(1. + utsq) ;
  rho0 = D/gamma ;

  w = W/(gamma*gamma) ;
  p = pressure_rho0_w_2d(rho0,w) ;
  u = w - (rho0 + p) ;

  if( (rho0 <= 0.) || (u <= 0.) ) { 
    if( debugfail>=2 ) dualfprintf(fail_file,"utsq failure:  utsq = %21.15g , W = %21.15g \n",utsq, W) ;
    
    if((rho0<=0.)&&(u>=0.)) *glpflag=  UTOPRIMFAILRHONEG;
    if((rho0>=0.)&&(u<=0.)) *glpflag= UTOPRIMFAILUNEG;
    if((rho0<=0.)&&(u<=0.)) *glpflag= UTOPRIMFAILRHOUNEG;
    if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNNOTADJUSTED) return(0) ; // else let assign -- used to check how bad failure is.
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;


  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = 
                     -gamma*Qtcon[i]/(W + Bsq) 
                     - gamma*QdotB*Bcon[i]/(W*(W + Bsq)) ;
 
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  if( ltrace ) {
    fprintf(stdout," rho final = %21.15g ,  u final = %21.15g \n", rho0, u);
  }


  /* done! */
  return(0) ;
}



/*************************************************
  evaluate \tilde{u}^2 from W = \gamma^2 w 
*************************************************/
static FTYPE utsq_calc(FTYPE W)
{
  FTYPE Wsq,W4,utsq ;
  // extern FTYPE Bsq,QdotBsq,Qtsq ;
 
  Wsq = W*W ;
  W4 = Wsq*Wsq ;

  //stderrfprintf("enter utsq %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",utsq,Bsq,QdotBsq,Qtsq,W4,Wsq) ;
  utsq = (Bsq*QdotBsq + W*(2.*QdotBsq + Qtsq*W))/
    (W4 + 2.*Bsq*Wsq*W + Wsq*(Bsq*Bsq - Qtsq) 
     - QdotBsq*(2.*W + Bsq)) ; 

  /*
    if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    //stderrfprintf("utsq failure %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",utsq,Bsq,QdotBsq,Qtsq,W4,Wsq) ;
    return(0.) ;
    }
  */

  return(utsq) ; 
}


/* find the contravariant fluid four-velocity from primitive 
   variables plus the metric */
static void ucon_calc_g(FTYPE prim[8],FTYPE gcov[GIND(4,4)],FTYPE gcon[GIND(4,4)],FTYPE ucon[4])
{
  FTYPE u_tilde_con[4] ;
  FTYPE u_tilde_sq ;
  FTYPE gamma,lapse ;
  int i,j ;
 
  u_tilde_con[0] = 0. ;
  u_tilde_con[1] = prim[UTCON1] ;
  u_tilde_con[2] = prim[UTCON2] ;
  u_tilde_con[3] = prim[UTCON3] ;

  u_tilde_sq = 0. ;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      u_tilde_sq += gcov[GIND(i,j)]*u_tilde_con[i]*u_tilde_con[j] ;
  u_tilde_sq = fabs(u_tilde_sq) ;

  gamma = sqrt(1. + u_tilde_sq) ;

  lapse = sqrt(-1./gcon[GIND(0,0)]) ;

  for(i=0;i<4;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[GIND(0,i)] ;

  return ;
}

/* raise covariant vector vcov using gcon, place result in vcon */
static void raise_g(FTYPE vcov[4], FTYPE gcon[GIND(4,4)], FTYPE vcon[4])
{
  int i,j;

  for(i=0;i<4;i++) {
    vcon[i] = 0. ;
    for(j=0;j<4;j++) 
      vcon[i] += gcon[GIND(i,j)]*vcov[j] ;
  }

  return ;
}
/* lower contravariant vector vcon using gcov, place result in vcov */
static void lower_g(FTYPE vcon[4], FTYPE gcov[GIND(4,4)], FTYPE vcov[4])
{
  int i,j;

  for(i=0;i<4;i++) {
    vcov[i] = 0. ;
    for(j=0;j<4;j++) 
      vcov[i] += gcov[GIND(i,j)]*vcon[j] ;
  }

  return ;
}

/* set covariant normal observer four-velocity */
static void ncov_calc(FTYPE gcon[GIND(4,4)],FTYPE ncov[4]) 
{
  FTYPE lapse ;

  lapse = sqrt(-1./gcon[GIND(0,0)]) ;

  ncov[0] = -lapse ;
  ncov[1] = 0. ;
  ncov[2] = 0. ;
  ncov[3] = 0. ;

  return ;
}

/* calculate contravariant magnetic field four-vector b */
static void bcon_calc_g(FTYPE prim[8],FTYPE ucon[4],FTYPE ucov[4],FTYPE ncov[4],FTYPE bcon[4]) 
{
  FTYPE Bcon[4] ;
  FTYPE u_dot_B ;
  FTYPE gamma ;
  int i ;

  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = prim[BCON1+i-1] ;

  u_dot_B = 0. ;
  for(i=0;i<4;i++) u_dot_B += ucov[i]*Bcon[i] ;

  gamma = -ucon[0]*ncov[0] ;
  for(i=0;i<4;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/*************************************************************** 
   This routine checks to see if the argument is NaN or Inf.
   Tested with an Intel(TM) compiler.
*****************************************************************/

static int is_nan_inf( FTYPE x ) 
{

  char x_char[50];

  sprintf(x_char, "%21.15g", x);
  if( !( strcmp("nan", x_char) && strcmp("inf", x_char) && strcmp("-inf", x_char) ) ) {
    return(1);
  }
  else {
    return(0);
  }

}


/********************************************************************

   find_root_2D_gen(): 
  
         -- performs a 2D Newton-Raphson method for solving the \tilde{Q}^2
             and Q.n equations for two variables in general;
              (usually W and something else like utsq, vsq, or gamma)

         residual vector       = { Qtsq eq. , Q.n eq. }
         indep. var. vector, x = { W, vsq/utsq/gamma/? }

       vartype = 2 : W, gamma
       vartype = 3 : W, vsq
       vartype = 4 : W, utsq
    
     -- Uses  general_newton_raphson();

*********************************************************************/

static int find_root_2D_gen(FTYPE *x0) 
{
  
  FTYPE x[2], x_orig[2];
  int ntries = 0;
  int n = 2;
  int ret;
  const int ltrace = 0;
  //  extern void func_gamma( FTYPE x[2], FTYPE dx[2], FTYPE *f, FTYPE *df, int n);
  //  extern void func_vsq(   FTYPE x[2], FTYPE dx[2], FTYPE *f, FTYPE *df, int n);
  //  extern void func_utsq(  FTYPE x[2], FTYPE dx[2], FTYPE *f, FTYPE *df, int n);

  /* Set presets: */
  x_orig[0] = x[0] = fabs(*x0) ;
  x_orig[1] = x[1] = x1_of_x0( *x0, 3 ) ;

  if( ltrace ) {
    fprintf(stdout, "find_root_2D_gen():  x[0] = %21.15g , x[1] = %21.15g \n", x[0], x[1]);
    fflush(stdout);
  }

  ret = general_newton_raphson( x, n, func_vsq );
  while( (ret != 0)  && (ntries < MAX_NEWT_RETRIES ) ) {
    x[0] = x_orig[0] * (1. + 0.2*(1.*rand())/(1.*RAND_MAX) );  
    x[1] = x1_of_x0( x[0], 3 ) ;
    ntries++;
  }

  if( (ntries >= MAX_NEWT_RETRIES) &&  (MAX_NEWT_RETRIES > 0)  ) {
    stderrfprintf( "find_root_2D_gen():  Bad exit value from general_newton_raphson() !! \n");
    stderrfprintf( "find_root_2D_gen():  ntries = %d , x[0] = %21.15g , x[1] = %21.15g  \n", ntries, x[0], x[1]); fflush(stderr);
  }

  if( ret != 0 ) {
    *x0 = FAIL_VAL;
    return( ret );
  }

  *x0 = x[0];
  return(0);
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates x1 from x0 depending on value of vartype;

                          x0   x1
             vartype = 2 : W, gamma
             vartype = 3 : W, vsq
             vartype = 4 : W, utsq
    
    -- asumes x0 is already physical

*********************************************************************/

static FTYPE x1_of_x0(FTYPE x0, int vartype ) 
{
  FTYPE x1, utsq;
  //  FTYPE dv = 1.e-15;
  

  utsq = fabs(utsq_calc(x0)) ;

  switch( vartype ) {

  case 2 : return( ( utsq > 0. ) ? (sqrt( utsq + 1. ))    : (1. ) );   /* gamma */

  case 3 : return( ( utsq > 0. ) ? (utsq / ( utsq + 1. )) : 0.0        );   /* vsq   */

  case 4 : return( ( utsq > 0. ) ? (utsq)                 : 0.0       );   /* utsq  */

  default : stderrfprintf( "\nx1_of_x0():  only defined for vartype=2,3,4   and vartype = %d !! \n", vartype); exit(1);

  }

  return(utsq);
}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon their definitions:

                          x0   x1
             vartype = 2 : W, gamma
             vartype = 3 : W, vsq
             vartype = 4 : W, utsq
    
*********************************************************************/

static void validate_x(FTYPE x[2], FTYPE x0[2], int vartype ) 
{
  
  //  FTYPE dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */ 
  x[0] = fabs(x[0]);
  x[0] = (x[0]/wglobal > W_TOO_BIG) ?  x0[0] : x[0];
  

  switch( vartype ) {

  case 2:   /* gamma */
    x[1] = (x[1] < 1.)           ?  (1.0) : x[1];    /* if it's too small */
    x[1] = (x[1] > UTSQ_TOO_BIG) ?  x0[1]     : x[1];    /* if it's too big   */
    break;

  case 3:  /* vsq */
    x[1] = (x[1] < 0.) ?   0.0       : x[1];  /* if it's too small */
    x[1] = (x[1] > 1.) ?  VSQ_TOO_BIG : x[1];  /* if it's too big   */
    break;

  case 4:  /* utsq */
    x[1] = (x[1] < 0.)           ?   0.0        : x[1];    /* if it's too small */
    x[1] = (x[1] > UTSQ_TOO_BIG) ?   x0[1]     : x[1];    /* if it's too big   */
    break;

  default : stderrfprintf( "\nvalidate_x0():  only defined for vartype=2,3,4   and vartype = %d !! \n", vartype); exit(1);

  }
}

/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson(
                                  FTYPE x[], int n,
                                  void (*funcd) (FTYPE [], FTYPE [], FTYPE *, FTYPE *, int) )
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], errx, errx_old, errx_oldest;
  int    n_iter, id,numdamped,allowdamp, numstable;
  FTYPE dampfactor, randtmp, xtmp;
  
  const int ltrace = 0;
  const int ltrace2 = 0;
  int    is_nan_inf( FTYPE x );

  dampfactor=1.0;

  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;

  numstable = numdamped = 0;
  allowdamp = 1;

  n_iter = 0;


  while ( ( errx > NEWT_TOL ) && (n_iter < MAX_NEWT_ITER)  ) {

    (*funcd) (x, dx, &f, &df, n);  /* returns with new dx, f, df */

    /*  Check for bad untrapped divergences : */
    xtmp = 0.;
    for( id = 0; id < n ; id++)  xtmp += dx[id] ;
    if( is_nan_inf( f ) || is_nan_inf( df ) || is_nan_inf( xtmp ) ) return(2);


    if( dodamp && ( f >= f_old + alpha_newt*df*dampfactor ) && (allowdamp) && (dampfactor > dampfactor_min) && (n_iter >= preiterdamp) ){
      /* Reduce stepsize and try again */

      dampfactor *= dampfactorchange;

      /* Reset x */
      for( id = 0; id < n ; id++)  x[id] = x_old[id] ;

      numdamped++;
      n_iter--;

      if( numdamped >= numdampedtot ){ allowdamp=0; numdamped=0; numstable=0;}

      if( ltrace ) {
        stderrfprintf("general_newton_raphson():  f_old = %21.15g ,  f = %21.15g , dampfactor = %21.15g \n",f_old, f,dampfactor);
        fflush(stderr);
      }
    }
    else{
      /* Normal Newton-Raphson step with damped stepsize */

      if( dodamp ) {
        if(allowdamp==0) numstable++;
        if(allowdamp) numdamped++; // coming here counts as a general damped run
        if(numstable>=numstabletot){ allowdamp=1; numdamped=0; numstable=0;}
        if(numdamped>=numdampedtot){ allowdamp=0; numdamped=0; numstable=0;}
      }

      if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
        randtmp = ( (1.*rand())/(1.*RAND_MAX) );
        for( id = 0; id < n ; id++) dx[id] *= randtmp;
        // for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
      }

      errx_oldest = errx_old;
      errx_old = errx;
      errx = 0.;
      for( id = 0; id < n ; id++) {
        x_old[id] = x[id] ;
        x[id] += dx[id] * dampfactor ;
        errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
      }
      errx /= 1.*n;


      validate_x( x, x_old, 3 ) ;


#if( CHECK_FOR_STALL )
      if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx < MIN_NEWT_TOL) )  errx = -errx;
#endif 

      if( ltrace ) {
        stderrfprintf(" general_newton_raphson(): niter = %d , f_old = %21.15g , f = %21.15g , errx_old = %21.15g , errx = %21.15g\n",  
                      n_iter,f_old,f,errx_old,errx );
        fflush(stderr);
      }

      for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
      f_old  = f;
      df_old = df;

    }
    n_iter++;
    nstroke++;
    //    lntries++;
  }


  if (n_iter >= MAX_NEWT_ITER) {
    if( errx > MIN_NEWT_TOL){
      if(ltrace2) {
        fprintf(stdout," totalcount = %d   0 \n",n_iter); 
        fflush(stdout);
      }
      if(ltrace) {
        stderrfprintf("general_newton_raphson():  did not find solution \n");
        fflush(stderr);
      }
      return(1);
    }
    else {
      if(ltrace2) {
        fprintf(stdout," totalcount = %d   1 \n",n_iter); 
        fflush(stdout);
      }
      if(ltrace) {
        stderrfprintf("general_newton_raphson(): found minimal solution \n");
        fflush(stderr);
      }
      return(0);
    }
  } 
  else {
    if(ltrace2) {
      fprintf(stdout," totalcount = %d   2 \n",n_iter); 
      fflush(stdout);
    }
    return(0);
  }

}



/*********************************************************

**********************************************************/

static void func_vsq(FTYPE x[2], FTYPE dx[2], FTYPE *f, FTYPE *df, int n)
{

  //  extern FTYPE Qtsq, Bsq, QdotBsq, Qdotn;
  
  FTYPE  W, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  const int ltrace = 0;


  W = x[0];
  vsq = x[1];
  
  Wsq = W*W;
  
  p_tmp  = pressure_W_vsq( W, vsq );
  dPdW   = dpdW_calc_vsq( W, vsq );
  dPdvsq = dpdvsq_calc( W, vsq );


  /* These expressions were calculated using Mathematica */
  /* Since we know the analytic form of the equations, we can explicitly
     calculate the Newton-Raphson step:                  */

  dx[0] = (-Bsq/2. + dPdvsq)*(Qtsq - vsq*((Bsq+W)*(Bsq+W)) + 
                              (QdotBsq*(Bsq + 2*W))/Wsq) + 
    ((Bsq+W)*(Bsq+W))*(Qdotn - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - 
                       W + p_tmp);

  dx[1] = -((-1 + dPdW - QdotBsq/(Wsq*W))*
            (Qtsq - vsq*((Bsq+W)*(Bsq+W)) + (QdotBsq*(Bsq + 2*W))/Wsq)) - 
    2*(vsq + QdotBsq/(Wsq*W))*(Bsq + W)*
    (Qdotn - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - W + p_tmp);

  detJ = (Bsq + W)*((-1 + dPdW - QdotBsq/(Wsq*W))*(Bsq + W) + 
                    ((Bsq - 2*dPdvsq)*(QdotBsq + vsq*(Wsq*W)))/(Wsq*W));
  
  dx[0] /= -(detJ) ;
  dx[1] /= -(detJ) ;

  tmp2 = (Qtsq - vsq*((Bsq+W)*(Bsq+W)) + (QdotBsq*(Bsq + 2*W))/Wsq);
  tmp3 = (Qdotn - (Bsq*(1 + vsq))/2. + QdotBsq/(2.*Wsq) - W + p_tmp);

  *df = -tmp2*tmp2 - tmp3*tmp3; 

  *f = -0.5 * ( *df );

  if( ltrace ) {
    fprintf(stdout,"func_vsq(): x[0] = %21.15g , x[1] = %21.15g , dx[0] = %21.15g , dx[1] = %21.15g , f = %21.15g \n", x[0], x[1],dx[0], dx[1], *f);
    fflush(stdout);
  }
}








/***********************************************************
 Converts primitive variables to conservative variables 
 as they are defined with respect to the new primitive 
 variable solver and not how they are used in the rest of HARM.
************************************************************/

// shouldn't use this function since not optimized to avoid catastrophic cancellations
static void primtoU_g(struct of_geom *ptrgeom, FTYPE *prim,FTYPE gcov[SYMMATRIXNDIM],FTYPE gcon[SYMMATRIXNDIM],FTYPE *U)
{
  int i,j ;
  FTYPE rho0 ;
  FTYPE ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],ncov[NDIM] ;
  FTYPE gamma,n_dot_b,bsq,u,p,w ;

  /* preliminaries */
  ucon_calc_g(prim,gcov,gcon,ucon) ;
  lower_g(ucon,gcov,ucov) ;
  ncov_calc(gcon,ncov) ;

  gamma = -ncov[0]*ucon[0] ;

  bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,gcov,bcov) ;

  n_dot_b = 0. ;
  for(i=0;i<NDIM;i++) n_dot_b += ncov[i]*bcon[i] ;
  bsq = 0. ;
  for(i=0;i<NDIM;i++) bsq += bcov[i]*bcon[i] ;

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho0,u) ;
  w = rho0 + u + p ;

  U[RHO] = gamma*rho0 ;

  for(i=0;i<NDIM;i++) 
    U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
      - (p + bsq/2.)*ncov[i] 
      + n_dot_b*bcov[i] ;

  U[BCON1] = prim[BCON1] ;
  U[BCON2] = prim[BCON2] ;
  U[BCON3] = prim[BCON3] ;

  return ;
}


/***********************************************************
************************************************************/


/*********************************************

 The following routines are based on the 
 particular equation of state being used. 

**********************************************/


/* 

   pressure as a function of rho0 and u 

   this is used by primtoU and Utoprim_?D

*/
static FTYPE pressure_rho0_u_2d(FTYPE rho0, FTYPE u)
{
  return((gam - 1.)*u) ;
}

/* 

   pressure as a function of rho0 and w = rho0 + u + p 

   this is used by primtoU and Utoprim_1D

*/
static FTYPE pressure_rho0_w_2d(FTYPE rho0, FTYPE w)
{
  return((gam-1.)*(w - rho0)/gam) ;
}


/* 

   pressure as a function of W, vsq, and D:

   This is needed by find_root_2D_vsq();

*/

static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq) 
{
  //  extern FTYPE D;
  FTYPE gtmp;
  

  gtmp = 1. - vsq;
  
  return(  (gam - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / gam  );

}


/* 

   partial derivative of pressure with respect to W

   this is needed for find_root_2D_vsq();

*/
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq)
{

  return( (gam - 1.) * (1. - vsq) /  gam ) ;

}

/* 

   partial derivative of pressure with respect to vsq

   this is needed for find_root_2D_vsq();

*/
static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq)
{
  //  extern FTYPE D;

  return( (gam - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / gam  ) ;

}



