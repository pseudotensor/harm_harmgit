
#include "decs.h"
#include "utoprim_1d.h" 


/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */

static FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;
static FTYPE wglobal;
static PFTYPE * glpflag; // global pflag for local file

/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

int Utoprim_1d(FTYPE U[NPR], struct of_geom *ptrgeom,  PFTYPE *lpflag, FTYPE *prim, FTYPE *pressure, struct of_newtonstats *newtonstats)
{

  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, k, ret; 
  FTYPE gcov[SYMMATRIXNDIM], gcon[SYMMATRIXNDIM], alpha;
  const int ltrace = 0;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"Utoprim_1d() called in parallel region\n");
    myexit(9366983);
  }
#endif

#if(WHICHEOS!=IDEALGAS)
  dualfprintf(fail_file,"This code does not handle non-ideal gases: Utoprim_1d()\n");
  myexit(1);
#endif


  // assign global int pointer to lpflag pointer
  glpflag=lpflag;


#if( WHICHVEL != VELREL4 )
  stderrfprintf("Utoprim_1d() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1);
#endif

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(0,0)]);
  

  /* Calculate the transform the CONSERVATIVE variables into the new system */
  // notice the -alpha on divT=0 related terms
  U_tmp[RHO] = alpha * U[RHO];
  for( i = UU; i <= U3; i++ ) {
    U_tmp[i] = -alpha * U[i] ;
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

#if(!OPTIMIZED)
  /* Check conservative variable transformation: */
  if( ltrace ) {
    for(i = 0; i < NDIM; i++ ) {
      for(j = 0; j < NDIM; j++ ) {
        gcov[GIND(i,j)] = ptrgeom->gcov[GIND(i,j)];
        gcon[GIND(i,j)] = ptrgeom->gcon[GIND(i,j)];
        fprintf(stdout,"gcov,gcon %d %d = %21.15g %21.15g \n", i, j, gcov[GIND(i,j)], gcon[GIND(i,j)]);fflush(stdout);
      }
    }
    fprintf(stdout,"gdet = %21.15g \n", ptrgeom->g);fflush(stdout);
    primtoU_g(ptrgeom, prim_tmp, gcov, gcon, U_tmp2 ); 
    for( i = 0; i < NPR; i++ ) {
      fprintf( stdout, "Utoprim_new2(): Utmp1[%d] = %21.15g , Utmp2[%d] = %21.15g , dUtmp[%d] = %21.15g \n", 
               i, U_tmp[i], i, U_tmp2[i], i, fabs( (U_tmp[i]-U_tmp2[i]) / ( (U_tmp2[i]!=0.) ? U_tmp2[i] : 1. ) )  ); 
    }
  }
#endif

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
  // extern FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;
  // extern void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  FTYPE x_1d[NEWT_DIM];
  FTYPE gcov[SYMMATRIXNDIM], gcon[SYMMATRIXNDIM];
  FTYPE QdotB,Bcon[4],Bcov[4],Qcov[4],Qcon[4],ncov[4],ncon[4],Qsq,Qtcon[4];
  FTYPE rho0,u,p,w,gammasq,gamma,W_last,W,utsq,tmpdiff ;
  int i,j, ret ;
  FTYPE Ui[8],primi[8];

  const int ltrace = 0;
  const int ltrace2 = 0;

#if(!OPTIMIZED)
  if( ltrace ) {
    for(i=0;i<8;i++) fprintf(stdout,"%d %21.15g %21.15g\n",i,prim[i],U[i]) ;
    fprintf(stdout,"debugfail = %d \n", debugfail);
    fflush(stdout);
  }
#endif

  // Assume ok initially:
  *glpflag = UTOPRIMNOFAIL;

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

  if( (utsq < 0.) && (fabs(utsq) < MAXNEGUTSQ) ) { 
    utsq = 0.0;
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG ) {
    if( debugfail>=2 ) dualfprintf(fail_file,"Utoprim_new2(): utsq < 0 in utoprim_1d attempt, utsq = %21.15g \n", utsq) ;
    *glpflag = UTOPRIMFAILCONVUTSQ;
    return(0) ;
  }

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u_1d(rho0,u) ;

  *pressure=p;

  w = rho0 + u + p ;
  gammasq = 1. + utsq ;
  W_last = w*gammasq ;

#if(!OPTIMIZED)
  if( ltrace ) {
    fprintf(stdout,"u = %21.15g,  p = %21.15g,  Bsq = %21.15g, Qsq = %21.15g \n",prim[RHO],prim[UU],Bsq,Qsq);
    fprintf(stdout,"Bcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
    fprintf(stdout,"Bcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcov[0],Bcov[1],Bcov[2],Bcov[3]);
    fprintf(stdout,"Qcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcon[0],Qcon[1],Qcon[2],Qcon[3]);
    fprintf(stdout,"Qcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcov[0],Qcov[1],Qcov[2],Qcov[3]);
    fprintf(stdout,"call find_root\n") ;  fflush(stdout);
    fflush(stdout);
  }
#endif

  x_1d[0] = W_last;
  wglobal=w;
  ret = general_newton_raphson( x_1d, 1, USE_LINE_SEARCH_ALWAYS, func_1d_orig, res_sq_1d_orig );
#if(!OPTIMIZED)
  if( ret ) { 
    fprintf( stderr, "general_newton_raphson() failed, x_1d[0] = %21.15g \n", x_1d[0] );  fflush(stderr);
  }
#endif
  W = x_1d[0];

#if(!OPTIMIZED)
  if( ltrace ) {
    stderrfprintf("(W,W_last,Bsq,Qtsq,QdotB,gammasq,Qdotn,gam) %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
                  W,W_last, Bsq,Qtsq,QdotB,gammasq,Qdotn,gam) ;
  }

  if( ltrace2 ) {
    fprintf(stdout, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    fflush(stdout);
  }
#endif
  
  /* Problem with solver, so return denoting error before doing anything further */
  if( (ret != 0) || (W == FAIL_VAL) ) {
    if( debugfail>=2 ) dualfprintf(fail_file, "Failed to find a prim. var. solution!! %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g \n",W_last,Bsq,QdotBsq,Qdotn,D,Qtsq);
    *glpflag = ret+UTOPRIMFAILCONVRET+1;// related to UTOPRIMFAILCONVRET
    return(0);
  }
  else{
    Wtest=W/wglobal;

    if(Wtest <= 0. || Wtest > W_TOO_BIG) {
      //if(Wtest <= 0.) {
      if( debugfail>=2 ) dualfprintf(fail_file,"Wtest failure %21.15g \n",Wtest) ;
      *glpflag = UTOPRIMFAILCONVW;
      return(0) ;
    }

    // below UTSQ_TOO_BIG covers the above TOO_BIG
    utsq = utsq_calc(W) ;
    if( (utsq < 0.) && (fabs(utsq) < MAXNEGUTSQ) ) { 
      utsq = 0.0;
    }
    if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
      if( debugfail>=2 ) dualfprintf(fail_file,"utsq failure:  utsq = %21.15g , W = %21.15g \n",utsq, W) ;
      if( fabs(utsq) > 1.0e-14 ) {
        *glpflag = UTOPRIMFAILCONVGUESSUTSQ;
        return(0) ;
      }
    }
  }

#if(!OPTIMIZED)
  if( ltrace ) {
    fprintf(stdout,"done find_root\n") ; fflush(stdout);
  }
#endif

  /* Past all checks, so now use good result to find final values of all primitive variables */
  gammasq = 1. + utsq ;
  gamma = sqrt(gammasq);
  rho0 = sqrt(D*D/gammasq) ;

  w = W / gammasq;
  p = pressure_rho0_w_1d(rho0,w) ;
  u = w - rho0 - p ;

  if( (rho0 <= 0.) || (u <= 0.) ) { 
    if( debugfail>=2 ) {
      tmpdiff = w - rho0;
      dualfprintf(fail_file,
                  "rho or uu < 0 failure: rho,w,(w-rho),p,u  = %21.15g %21.15g %21.15g %21.15g %21.15g \n",
                  rho0,w,tmpdiff,p,u) ;
      dualfprintf(fail_file,
                  "rho or uu < 0 failure: gamma,utsq = %21.15g %21.15g  \n",  gamma, utsq) ;
    }
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

#if(!OPTIMIZED)
  if( ltrace ) {
    fprintf(stdout," debugfail,  final rho, uu = %d  %21.15g  %21.15g \n", debugfail, rho0, u);
  }
#endif

  /* done! */
  return(0) ;
}



/*************************************************
  evaluate \tilde{u}^2 from W = \gamma^2 w 
*************************************************/
static FTYPE utsq_calc(FTYPE W)
{
  FTYPE Wsq,W4,utsq,tmp ;
  // extern FTYPE Bsq,QdotBsq,Qtsq ;
  const int ltrace = 0;
 
  Wsq = W*W ;
  W4 = Wsq*Wsq ;

#if(!OPTIMIZED)
  if( ltrace ) {
    stderrfprintf("enter utsq %21.15g %21.15g %21.15g %21.15g %21.15g\n",Bsq,QdotBsq,Qtsq,W4,Wsq) ;
    fflush(stderr);
  }
#endif

  // utsq = (Bsq*QdotBsq + W*(2.*QdotBsq + Qtsq*W))/
  //  (W4 + 2.*Bsq*Wsq*W + Wsq*(Bsq*Bsq - Qtsq) 
  //   - QdotBsq*(2.*W + Bsq)) ; 
  //

  tmp = QdotBsq * ( Bsq + 2*W ) + Qtsq * Wsq;


  utsq = tmp / (  Wsq*(W+Bsq)*(W+Bsq) - tmp ) ;

  /*
    if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    stderrfprintf("utsq failure %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",utsq,Bsq,QdotBsq,Qtsq,W4,Wsq) ;
    return(0.) ;
    }
  */

  return(utsq) ; 
}

/*************************************************
  evaluate gamma^2
*************************************************/
static FTYPE gammasq_calc(FTYPE W)
{
  FTYPE Wsq,gammasq,tmp,tmp2 ;
  // extern FTYPE Bsq,QdotBsq,Qtsq ;
  const int ltrace = 0;
 
  Wsq = W*W ;

#if(!OPTIMIZED)
  if( ltrace ) {
    stderrfprintf("enter utsq  %21.15g %21.15g %21.15g %21.15g\n",Bsq,QdotBsq,Qtsq,Wsq) ;
    fflush(stderr);
  }
#endif

  tmp = QdotBsq * ( Bsq + 2*W ) + Qtsq * Wsq;
  tmp2 = Wsq*(W+Bsq)*(W+Bsq);


  gammasq = tmp2 / (  tmp2 - tmp ) ;

  return(gammasq) ; 
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
  // extern void my_lnsrch(int, FTYPE [], FTYPE, FTYPE [], FTYPE [], FTYPE [], FTYPE *,  FTYPE, FTYPE, int *, FTYPE (*res_func) (FTYPE []));
  //void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) ;
  int   keep_iterating, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 0;
  int    is_nan_inf( FTYPE x );

  retval = 0;

  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;


  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 
    nstroke++;

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    if( n_iter == 0 ) { 
      resid_norm = 0.0e0;
      for( id = 0; id < n ; id++) {
        resid_norm += fabs(resid[id]);
      }
      resid_norm /= 1.0*n ;
      if( resid_norm == 0.0 ) resid_norm = 1.0;
    }
       

    /*  Check for bad untrapped divergences : */
    if( is_nan_inf( f ) || is_nan_inf( df ) ) {
      if(debugfail>=2) dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      return(1);
    }

    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      // for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /***************************************************/
    /* Make the newton step:                           */
    /***************************************************/

    /* Perform search along direction of the Newton step: */
    if( do_line_search == 1 ) { 

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
        stderrfprintf("gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
                      f_old, f, res_func_old, res_func_val );
        stderrfprintf("gnr(): x_old = ");
        for( id = 0; id < n ; id++) {
          stderrfprintf(" %21.15g ",x_old[id]);
        }
        stderrfprintf("\n ");
        stderrfprintf("gnr(): x     = ");
        for( id = 0; id < n ; id++) {
          stderrfprintf(" %21.15g ",x[id]);
        }
        stderrfprintf("\n ");
        stderrfprintf("gnr(): dn    = ");
        for( id = 0; id < n ; id++) {
          stderrfprintf(" %21.15g ",dn[id]);
        }
        stderrfprintf("\n ");
        stderrfprintf("gnr(): del_f = ");
        for( id = 0; id < n ; id++) {
          stderrfprintf(" %21.15g ",del_f[id]);
        }
        stderrfprintf("\n ");
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
#if(!OPTIMIZED)
        if( ltrace && retval ) { 
          stderrfprintf("general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
          fflush(stderr);
        }     
#endif
  
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH_ALWAYS == do_line_search) ) { 
        if( ltrace ) { 
#if(!OPTIMIZED)
          stderrfprintf("gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
          stderrfprintf("gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
          fflush(stderr);
#endif
          retval = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
          for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
          return( retval );
        }
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
#if(!OPTIMIZED)
          stderrfprintf("general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval, grad_check);
          fflush(stderr);
#endif
        }
      }
    }
    else {
      /* don't use line search : */
      for( id = 0; id < n ; id++) {
        x[id] += dx[id] ;
      }
    }

    /* Calculate the convergence criterion */
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;

    /* Make sure that the new x[] is physical : */
    x[0] = fabs(x[0]);
    x[0] = (x[0] > W_TOO_BIG) ?  x_old[0] : x[0];
    

#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL) )  errx = -errx;
#endif 

    /* If there's a problem with line search, then stop iterating: */
    if( (retval == 1) || (retval == -1) ) errx = -errx;

#if(!OPTIMIZED)
    if( ltrace ) {
      stderrfprintf(" general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
                    n_iter,f_old,f,errx_old,errx );
      stderrfprintf("gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
        stderrfprintf(" %21.15g ",x_old[id]);
      }
      stderrfprintf("\n ");
      stderrfprintf("gnr(): x     = ");
      for( id = 0; id < n ; id++) {
        stderrfprintf(" %21.15g ",x[id]);
      }
      stderrfprintf("\n ");
      fflush(stderr);
    }
#endif

    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;
  }

  if( (fabs(errx) > MIN_NEWT_TOL) && (i_extra == 0) ){
    if( (do_line_search != USE_LINE_SEARCH_ALWAYS) || (USE_LINE_SEARCH_ALWAYS < 0) ) { 
      if(ltrace2) {
        // stderrfprintf(" totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
        // fflush(stderr);
#if(!OPTIMIZED)
        bin_newt_data( errx, n_iter, 0, 0 );
#endif
      }
#if(!OPTIMIZED)
      if(ltrace) {
        stderrfprintf("general_newton_raphson():  did not find solution \n");
        if( retval == -1 ) {
          stderrfprintf("general_newton_raphson(): lnsrch converged: x = ");
          for( id = 0; id < n ; id++)  stderrfprintf(" %21.15g  ",x[id]);
          stderrfprintf("\n");
          stderrfprintf("general_newton_raphson(): lnsrch converged: x_old = ");
          for( id = 0; id < n ; id++)  stderrfprintf(" %21.15g  ",x_old[id]);
          stderrfprintf("\n");
        }
        fflush(stderr);
      }
#endif
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
#if(!OPTIMIZED)
      stderrfprintf("gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n",
                    do_line_search, retval, errx );
      fflush(stderr);
#endif
      retval2 = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      return( retval2 );
    }
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    if(ltrace2) {
      // stderrfprintf(" totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      // fflush(stderr);
      bin_newt_data( errx, n_iter, 1, 0 );
    }
#if(!OPTIMIZED)
    if(ltrace) {
      stderrfprintf("general_newton_raphson(): found minimal solution \n");
      fflush(stderr);
    }
#endif
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    if(ltrace2) {
      // stderrfprintf(" totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      // fflush(stderr);
      bin_newt_data( errx, n_iter, 2, 0 );
    }
    return(0);
  }

  return(0);

}


/***********************************************************
 Converts primitive variables to conservative variables 
 as they are defined with respect to the new primitive 
 variable solver and not how they are used in the rest of HARM.
************************************************************/


/**************************************************** 
*****************************************************/

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE (*jac)[NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE utsq,gamma,gamma_sq,rho0,w,W,drdW,dpdW,Wsq,p_tmp, dg ;
  // extern FTYPE Bsq,QdotBsq,Qdotn,D ;
  FTYPE dgamma_dW( FTYPE W, FTYPE gamma );
  FTYPE  dpress_dW( FTYPE W, FTYPE gamma, FTYPE dg );
  FTYPE dv = 1.0e-15;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  utsq = fabs(utsq_calc(W)) ;

  gamma = ( utsq < 0. ) ? (1. + dv) : sqrt(1. + utsq) ;
  gamma_sq = gamma*gamma;

  rho0 = D/gamma ;

  w = W/gamma_sq;
  p_tmp = pressure_rho0_w_1d(rho0,w);

  resid[0] = 
    +W 
    +  Bsq * ( 1. - 0.5/gamma_sq )
    - 0.5*QdotBsq/Wsq
    - Qdotn
    - p_tmp;

  dg =   dgamma_dW( W, gamma );

  dpdW = dpress_dW( W, gamma, dg );
  jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + Bsq*dg/(gamma*gamma_sq);

  dx[0] = -resid[0]/drdW;

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

#if(!OPTIMIZED)
  if( ltrace ) {
    fprintf(stdout,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
    fprintf(stdout,"func_1d_orig(): dg = %21.15g,  dp = %21.15g \n", dg, dpdW);
    fflush(stdout);
  }
#endif

}

/**************************************************** 
*****************************************************/

static FTYPE res_sq_1d_orig(FTYPE x[])
{
  FTYPE utsq,gamma,gamma_sq,rho0,w,W,drdW,dpdW,Wsq,p_tmp, dg ;
  /* extern int nresideval ;*/
  // extern FTYPE Bsq,QdotBsq,Qdotn,D ;
  FTYPE dgamma_dW( FTYPE W, FTYPE gamma );
  FTYPE  dpress_dW( FTYPE W, FTYPE gamma, FTYPE dg );
  FTYPE resid[1];
  FTYPE dv = 1.0e-15;

  W = x[0];
  Wsq = W*W;

  utsq = fabs(utsq_calc(W)) ;

  gamma = ( utsq < 0. ) ? (1. + dv) : sqrt(1. + utsq) ;
  gamma_sq = gamma*gamma;

  rho0 = D/gamma ;

  w = W/gamma_sq;
  p_tmp = pressure_rho0_w_1d(rho0,w);

  resid[0] = 
    +W 
    +  Bsq * ( 1. - 0.5/gamma_sq )
    - 0.5*QdotBsq/Wsq
    - Qdotn
    - p_tmp;

  return(  0.5*resid[0]*resid[0] );

}


/**************************************************** 
*****************************************************/

static FTYPE dgamma_dW(FTYPE W, FTYPE gamma)
{
  FTYPE tmp, tmp2, tmp3, Wsq;
  // extern FTYPE Bsq, QdotBsq, Qtsq;
  
  Wsq = W*W;

  tmp =  -QdotBsq * ( Bsq + 2*W ) + 
    Wsq * ( Bsq * (Bsq + 2*W) - Qtsq + Wsq );


  tmp3 = W*(Bsq+W);

  if( tmp > 0. ) {
    tmp2 = tmp3 * sqrt(tmp);

    return( -((W*(Bsq + W)*( Bsq*QdotBsq*( Bsq + 3*W ) + 
                             Wsq*(3*QdotBsq + Qtsq*W))) / (tmp*tmp2))    );
  }
  else {
    return( 0.5*gamma/W ); 
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


/******************************************************************************************
 ******************************************************************************************

 The following routines are based on the 
 particular equation of state being used. 

**********************************************/

/*
  Pressure as a function of only gamma and D 
*/

static FTYPE dpress_dW( FTYPE W, FTYPE gamma, FTYPE dg )
{
  // extern FTYPE D;

  return(  ((gam - 1.)*(-2*dg*W + gamma*(1.0 + D*dg)))/(gamma*gamma*gamma*gam) );

}


/* 
   pressure as a function of rho0 and u 
   this is used by primtoU and Utoprim_?D
*/
static FTYPE pressure_rho0_u_1d(FTYPE rho0, FTYPE u)
{
  return((gam - 1.)*u) ;
}

/* 
   pressure as a function of rho0 and w = rho0 + u + p 
   this is used by primtoU and Utoprim_1D
*/
static FTYPE pressure_rho0_w_1d(FTYPE rho0, FTYPE w)
{
  return((gam-1.)*(w - rho0)/gam) ;
}


/**************************************************** 
*****************************************************/

#define ALF 1.0e-4
/*#define TOLX 1.0e-7*/

static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
                       FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) )
{
  int i;
  FTYPE a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  int icount=0;
  int is_nan_inf( FTYPE );
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
    //    temp=fabs(p[i])/max(fabs(xold[i]),1.0);
    temp= (xold[i] == 0.) ? fabs(p[i]) :  fabs(p[i]/xold[i]);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;

#if(!OPTIMIZED)
  if( ltrace ) {
    stderrfprintf("my_lnsrch(): sum, slope, test, alamin =   %21.15g  %21.15g  %21.15g  %21.15g \n",sum,slope,test,alamin); 
    fflush(stderr);
  }
#endif

  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

    *f=(*func)(x+1);
      
    bad_step = 0;

    if( is_nan_inf(*f) ) { 
      bad_step = 1;
    }

    if( x[1] <= 0. ) { 
      bad_step = 1;
    }

    if( bad_step ) { 
      *check = 2;
      if(debugfail>=2) stderrfprintf("my_lnsrch(): bad_step = 1,  f = %21.15g \n", *f); fflush(stderr);
      return;
    }
      
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
#if(!OPTIMIZED)
      if( ltrace ) { 
        stderrfprintf("my_lnsrch(): alam < alamin: alam, alamin = %21.15g  %21.15g \n", alam,alamin); fflush(stderr);
      }
#endif
      return;
    } 
    else if (*f <= fold+ALF*alam*slope) {
#if(!OPTIMIZED)
      if( ltrace ) { 
        stderrfprintf("my_lnsrch(): good exit:  alam, alamin, f, fold = %21.15g  %21.15g %21.15g  %21.15g \n", alam,alamin, *f, fold); fflush(stderr);
      }
#endif
      return;
    }
    else {
      if (alam == 1.0) {
        tmplam = -slope/(2.0*(*f-fold-slope));
#if(!OPTIMIZED)
        if( ltrace ) {
          stderrfprintf("my_lnsrch(): setting tmplam!!    tmplam, alam =  %21.15g  %21.15g !!\n", tmplam, alam);
          fflush(stderr);
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
            if( disc < -1.e-10 ) {
              if(debugfail>=2) stderrfprintf("my_lnsrch(): Big Roundoff problem:  disc = %21.15g \n", disc);
            }
            disc = 0.;
          }
          else tmplam=(-b+sqrt(disc))/(3.0*a);
        }
        if (tmplam>0.5*alam)
          tmplam=0.5*alam;
#if(!OPTIMIZED)
        if( ltrace ) {
          stderrfprintf("my_lnsrch(): rhs1, rhs2, a, b, tmplam, alam =  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g !!\n",
                        rhs1, rhs2, a, b, tmplam, alam );
          fflush(stderr);
        }
#endif
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);

#if(!OPTIMIZED)
    if( ltrace ) {
      stderrfprintf("my_lnsrch(): icount, alam, alam2, tmplam  =   %4i  %21.15g  %21.15g  %21.15g \n",
                    icount, alam, alam2, tmplam);
      fflush(stderr);
    }
#endif

    icount++;
  }
}
#undef ALF
/* #undef TOLX */



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

  // extern SFTYPE t;
  /* General variables */
  int i, j;
  static int first_call   = 1;
  static long int n_bin_calls  = 0L;
  const  long int n_output_freq = 50000L;
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


  if(OPTIMIZED) return;

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
      fprintf(stdout,"bin_newt_data(): bad value for niters = %d \n", niters );
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

    fflush(stdout);
    fprintf(stdout,"t = %21.15g ,  n_beyond_range = %ld , n_bin_calls = %ld \n", t, n_beyond_range, n_bin_calls);

    /* ERRX */
    fprintf(stdout,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");
    fprintf(stdout,"                         x");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      fprintf(stdout,"            N%d",j);
    }
    fprintf(stdout,"\n");
    
    for( i = 0; i < NBINS; i++ ) { 
      fprintf(stdout,"%21.15g ",xbin[i]);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        fprintf(stdout,"%13ld ", n_errx[j][i]);
      }
      fprintf(stdout,"\n");
    }    
    fprintf(stdout,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");


    /* NITER */

    fprintf(stdout,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");
    fprintf(stdout,"        niter");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      fprintf(stdout,"            N%d",j);
    }
    fprintf(stdout,"\n");
    
    for( i = 0; i < N_NITER_BINS; i++ ) { 
      fprintf(stdout,"%13d ", i);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
        fprintf(stdout,"%13ld ", n_niters[j][i]);
      }
      fprintf(stdout,"\n");
    }    
    fprintf(stdout,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");

    fflush(stdout);
  }
    
    
  if( print_now != 1 )  n_bin_calls++;


  return;

}
  
#undef N_CONV_TYPES 
#undef N_NITER_BINS 
#undef NBINS 

