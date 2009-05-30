
/* 
 *
 Generates tests designed by Komissarov 1999-2004 in 1D force-free

 Todo when using:
1) global.h : resolution to 200x1 probably
2)            EOMTYPE EOMFFDE
3) copy bounds.komtests.c bounds.c
4) use defcoord with uniform grid (modify defcoord==0)
5) probably play with lim, TIMEORDER, FLUXB, dt, etc. in initbase.c
 *
 */

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

// see below
#define TESTNUMBER 8


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  return(0);
}

int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  return(0);
}

int init_grid(void)
{
  SFTYPE rh;
  
  
  // define coordinate type
  defcoord = 0;


  Rin=-.5;
  Rout =1.5 ;

#if(TESTNUMBER==0)
  Rin=-0.5;
  Rout=1.5 ;
#endif  
#if(TESTNUMBER==1)
  Rin=-0.5;
  Rout=1.5 ;
#endif
#if(TESTNUMBER==2)
  Rin=-1.5;
  Rout=0.5 ;
#endif  
#if(TESTNUMBER==3)
  Rin=-0.5;
  Rout=1.5 ;
#endif  
#if(TESTNUMBER==4)
  Rin=-0.5;
  Rout=1.5 ;
#endif  
#if(TESTNUMBER==5)
  Rin=-0.5;
  Rout=1.5 ;
#endif  
#if(TESTNUMBER==6)
  Rin=-0.5;
  Rout=1.5 ;
#endif  
#if(TESTNUMBER==7)
  Rin=-1.0;
  Rout=2.0 ;
#endif  
#if(TESTNUMBER==8)
  Rin=-1.0;
  Rout=2.0 ;
#endif  

  return(0);
}

int init_global(void)
{

  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=OUTFLOW;
  BCtype[X2DN]=OUTFLOW;

  /* output choices */
  tf = 1;

  DTd = .1;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;			/* logfile frequency, in units of M */
  DTi = 10;			/* image file frequ., in units of M */
  DTdebug = DTd; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */

  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
{
  int k;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];

  
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}


// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, FTYPE *pr)
{
  extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr);
  void computeKK(FTYPE *pr, struct of_geom *geom, FTYPE *KK);
  void EBvetatopr(FTYPE *Econ, FTYPE *Bcon, FTYPE *veta, struct of_geom *geom, FTYPE *pr);
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  int k;
  FTYPE E[NDIM],B[NDIM];
  FTYPE x0,dx0;
  FTYPE bcon[NDIM],vcon[NDIM],econ[NDIM];
  FTYPE phi0;
  FTYPE KK;
  FTYPE B0;
  FTYPE muf;

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(i, j, CENT, &geom); // true coordinate system

  pr[RHO]=pr[UU]=0;
  pr[U1]=pr[U2]=pr[U3]=0.0;
  pr[B2]=pr[B3]=0;
  pr[B1]=0;



#if(TESTNUMBER==0) // Fast wave
  tf = 1;
  DTd=tf/10.0;

  //  tf = 1;
  //  DTd=1E-5;
  

  E[1]=0;
  E[2]=0;
  B[3]=0.0;
  B[1]=0.0;
  x0=0.0;
  dx0=0.1;
  if(r-x0<=-dx0) B[2]=1.0;
  else if((r-x0>-dx0)&&(r-x0<dx0)) B[2]=1.0-(0.3/(2.0*dx0))*((r-x0)+dx0);
  else if(r-x0>=dx0) B[2]=0.7;

  muf=1.0;

  E[3]=1.0-muf*B[2];

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  //  pr[U1]=0.9;

  //dualfprintf(fail_file,"pr[U1]=%21.15g pr[U2]=%21.15g\n",pr[U1],pr[U2]);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

#endif
#if(TESTNUMBER==1) // comoving Fast wave (NOT a Komissarov test)
  //tf = 1;
  //  DTd=tf/10.0;

  tf = 1;
  DTd=1E-5;
  

  bcon[3]=0.0;
  bcon[1]=0.0;
  x0=0.0;
  dx0=0.1;
  if(r-x0<=-dx0) bcon[2]=1.0;
  else if((r-x0>-dx0)&&(r-x0<dx0)) bcon[2]=1.0-(0.3/(2.0*dx0))*((r-x0)+dx0);
  else if(r-x0>=dx0) bcon[2]=0.7;

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  vcon[1]=0.9999;
  vcon[2]=vcon[3]=0;

  vbtopr(vcon,bcon,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

#endif
#if(TESTNUMBER==2) // (nondegenerate) Alfven wave
  tf = 2;
  DTd=tf/10.0;

  bcon[1]=bcon[2]=1.0;
  
  x0=0.0;
  if(r-x0<=-0.1) bcon[3]=1.0;
  else if((r-x0>-0.1)&&(r-x0<0.1)) bcon[3]=1.0+3.0/2.0*((r-x0)+0.1);
  else if(r-x0>=0.1) bcon[3]=1.3;

  econ[2]=econ[3]=0.0;

  // can be + or -
  //#define CONSTECON1 1.3
  //  econ[1]=-sqrt(-CONSTECON1+bcon[1]*bcon[1]+bcon[2]*bcon[3]+bcon[3]*bcon[3]);
  //  econ[1]=1-0.5*bcon[3];
#define CONSTECON1 1.0
  econ[1]=sqrt(-CONSTECON1 + bcon[3]*bcon[3]);
  //  econ[1]=0.0;
  //  econ[1]=-bcon[3];

  vcon[1]=-0.5;
  vcon[2]=vcon[3]=0;

  EBvetatopr(econ, bcon, vcon, &geom, pr);
  //  vbtopr(vcon,bcon,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);

#endif

#if(TESTNUMBER==3) // Degenerate Alfven wave
  tf = 2;
  DTd=tf/10.0;
  bcon[1]=0.0;

  
  x0=0.0;
  if(r-x0<=-0.1) phi0=0.0;
  else if((r-x0>-0.1)&&(r-x0<0.1)) phi0=5.0/2.0*M_PI*((r-x0)+0.1);
  else if(r-x0>=0.1) phi0=M_PI*0.5;

  bcon[2]=2.0*cos(phi0);
  bcon[3]=2.0*sin(phi0);


  vcon[1]=0.5;
  vcon[2]=vcon[3]=0;

  vbtopr(vcon,bcon,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif
#if(TESTNUMBER==4) // Three-wave problem
  tf = .75;
  DTd=tf/10.0;

  x0=0.5;
  if(r<x0){
    B[1]=1.0;
    B[2]=1.5;
    B[3]=3.5;
    E[1]=-1.0;
    E[2]=-0.5;
    E[3]=0.5;
  }
  else{
    B[1]=1.0;
    B[2]=2.0;
    B[3]=2.3;
    E[1]=-1.5;
    E[2]=1.3;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif

#if(TESTNUMBER==5) // B^2-E^2<0 problem
  tf = .02;
  DTd=tf/10.0;

  x0=0.5;
  if(r<x0){
    B[0]=0.0;
    B[1]=1.0;
    B[2]=1.0;
    B[3]=1.0;
    E[0]=0.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else{
    B[0]=0.0;
    B[1]=1.0;
    B[2]=-1.0;
    B[3]=-1.0;
    E[0]=0.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif
#if(TESTNUMBER==6) // smoothed B^2-E^2<0 problem
  tf = .02;
  DTd=tf/10.0;

  x0=0.5;
  if(r-x0<-0.1){
    B[1]=1.0;
    B[2]=1.0;
    B[3]=1.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else if(r-x0>0.1){
    B[1]=1.0;
    B[2]=-1.0;
    B[3]=-1.0;
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }
  else if((r-x0>=-0.1)&&(r-x0<=0.1)){
    B[1]=1.0;
    B[2]=1.0+(r-x0+0.1)*(-2.0/0.2);
    B[3]=1.0+(r-x0+0.1)*(-2.0/0.2);
    E[1]=0.0;
    E[2]=0.5;
    E[3]=-0.5;
  }

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention


  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif

#if(TESTNUMBER==7) // Komissarov 2004 C3.1 Alfven wave
  // PARA generates crap on left side, but wave doesn't move
  // MC does very well
  // no obvious difference between HLL and LAXF
  // Athena1/2 ok
  tf = 2.0;
  DTd=tf/10.0;

  //  B[1]=B[2]=E[3]=E[2]=0;
  B[1]=B[2]=E[3]=1.0;
  E[2]=0;


  if(r<=0.5){
    B[3]=1.0;
  }
  else if(r>=0.2+0.5){
    B[3]=1.3;
  }
  else{
    B[3]=1.0+0.15*(1.0+sin(5.0*M_PI*(r-0.1-0.5)));
  }
  E[1]=-B[3];

  //  for(k=1;k<=3;k++) E[k]=-E[k]; // switch for GRFFE formulation sign convention

  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif
#if(TESTNUMBER==8) // Komissarov 2004 C3.2 Current Sheet
  tf = 1.0;
  DTd=tf/10.0;

  E[1]=E[2]=E[3]=0.0;
  B[3]=0.0;
  B[1]=1.0;

  //B0=0.5; // fine
  B0 = 2.0;

  if(r<=0.5){
    B[2]=B0;
  }
  else{
    B[2]=-B0;
  }


  EBtopr(E,B,&geom,pr);
  //EBtopr_2(E,B,&geom,pr);

  computeKK(pr,&geom,&KK);

  dualfprintf(fail_file,"i=%d KK=%21.15g\n",i,KK);


#endif

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}


// assumes normal field in pr
int init_vpot(int i, int j,FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;


  return(0);



}

int init_vpot2field(FTYPE A[][N2M],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(FTYPE A[][N2M],FTYPE p[][N2M][NPR]);


  return(vpot2field(A,pr));

}

// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE X[NDIM],r,th;

  return(0);

}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM];
  FTYPE r,th;

  return(0);

}


#undef SLOWFAC
