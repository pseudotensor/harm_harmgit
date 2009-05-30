
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

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
  
  // metric stuff first
  a = 1.0 ;
  
  // define coordinate type
  defcoord = 0;
  set_coord_parms();

  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.8;
  Rhor=rhor_calc(0);
  Rout = 1E3;
 
  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;
  
  hslope = 0.1;

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
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

  /* output choices */
  tf = 1E3;

  DTd = 50;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;			/* logfile frequency, in units of M */
  DTi = 2;			/* image file frequ., in units of M */
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
int init_dsandvels(int *whichvel, int*whichcoord, int ii, int jj, FTYPE *pr)
{
  FTYPE X[NDIM],r,th;
  struct of_geom geom;

  FTYPE Fcov[NDIM][NDIM];
  FTYPE Mcon[NDIM][NDIM];
  FTYPE Mcov[NDIM][NDIM];
  FTYPE etacov[NDIM],etacon[NDIM];
  FTYPE Ecov[NDIM],Econ[NDIM],Bcov[NDIM],Bcon[NDIM];
  FTYPE  alpha;
  FTYPE Jcon[NDIM];
  int j,k;

  void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM]);
  void Jcon_numerical(FTYPE *X, FTYPE *Jcon);
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  extern void EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern void EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);


  struct of_state q;
  FTYPE faradaytest[NDIM][NDIM];


  if(EOMTYPE!=EOMFFDE){
    dualfprintf(fail_file,"Are you sure?\n");
    myexit(1);
  }
  

  coord(ii, jj, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(ii, jj, CENT, &geom); 

  // 0: split monopole
  // 1: monopole
  // 2: Wald
#define PROBLEMTYPE 2
#define B0 1.0


  pr[U1]=pr[U2]=pr[U3]=0.0;
  pr[B2]=pr[B3]=0;
  if(PROBLEMTYPE==0){
    if(th<M_PI*0.5)  pr[B1]=B0*gdet[horizoni][jj][CENT]/(gdet[ii][jj][CENT]);
    else pr[B1]=-B0*gdet[horizoni][jj][CENT]/(gdet[ii][jj][CENT]);
  }
  else if(PROBLEMTYPE==1){
    // Ruben's talk says they set $\dF^{tr} = C\sin{\theta}/\detg$.
    pr[B1]=B0*gdet[horizoni][N2-1-jj][CENT]/(gdet[ii][N2-1-jj][CENT]);
  }
  else if(PROBLEMTYPE==2){
    //    first get F_{\mu\nu}
    Fcov_numerical(X, Fcov);

  coord(ii, jj, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(ii, jj, CENT, &geom); 


    // check that J^\mu=0
    Jcon_numerical(X, Jcon);


  coord(ii, jj, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(ii, jj, CENT, &geom); 


    //    dualfprintf(fail_file,"Fcov\n");
    //DLOOP dualfprintf(fail_file,"Fcov[%d][%d]=%21.15g\n",j,k,Fcov[j][k]);
    //    dualfprintf(fail_file,"%21.15g %21.15g\n",j,k,Fcov[0][3],Fcov[3][0]);
    
    // lapse
    // define \eta_\alpha
    // assume always has 0 value for space components
    alpha = 1./sqrt(- geom.gcon[0][0]);

    etacov[TT]=-alpha; // any constant will work.
    SLOOPA etacov[j]=0.0; // must be 0
    
    // shift
    // define \eta^\beta
    raise(etacov,&geom,etacon);
    //    dualfprintf(fail_file,"raise\n");

    //    DLOOPA dualfprintf(fail_file,"etacon[%d]=%21.15g etacov[%d]=%21.15g\n",j,etacon[j],j,etacov[j]);
    

    // then get E^\alpha and B^\alpha
    DLOOPA Ecov[j]=0.0;
    DLOOP Ecov[j]+=etacon[k]*Fcov[j][k];
    raise(Ecov,&geom,Econ);
    //    dualfprintf(fail_file,"Ecov[3]=%2.15g\n",Ecov[3]);
    //DLOOPA dualfprintf(fail_file,"Econ[%d]=%2.15g\n",j,Econ[j]);


    MtoF(3,Fcov,&geom,Mcon);
    //    dualfprintf(fail_file,"MtoF\n");
    //    DLOOP dualfprintf(fail_file,"Mcon[%d][%d]=%21.15g\n",j,k,Mcon[j][k]);

    DLOOPA Bcon[j]=0.0;
    DLOOP Bcon[j]+=etacov[k]*Mcon[k][j];

    lower(Bcon,&geom,Bcov);

    //    DLOOPA dualfprintf(fail_file,"Econ[%d]=%21.15g Bcon[%d]=%21.15g\n",j,Econ[j],j,Bcon[j]);

    EBtopr(Econ,Bcon,&geom,pr);
    //    dualfprintf(fail_file,"EBtopr\n");


    // stick J^\mu into dump file
#if(0)
    for(k=U1;k<=U1+3;k++){
      pr[k] = Jcon[k-U1];
    }
#endif
    dualfprintf(fail_file,"ii=%d jj=%d\n",ii,jj);

#if(0)
    // E.B
    k=U1;
    pr[k]=0;
    DLOOPA pr[k]+=Ecov[j]*Bcon[j];

    // B^2-E^2
    k=UU;
    pr[k]=0;
    DLOOPA pr[k]+=Bcov[j]*Bcon[j];
    //    DLOOPA pr[k]+=Bcon[j];
    DLOOPA pr[k]-=Ecov[j]*Econ[j];

    pr[B1]=Econ[1]*geom.g/etacov[TT];
#endif


#if(0)
    // check where faraday changed
    get_state(pr,&geom,&q);

    faraday_calc(0,q.bcon,q.ucon,&geom,faradaytest);
    //    DLOOP dualfprintf(fail_file,"%21.15g  %21.15g\n",faradaytest[j][k],Fcov[j][k]);
    DLOOP{
      if(fabs(faradaytest[j][k]-Fcov[j][k])>1E-10){
	dualfprintf(fail_file,"1 %d %d : %21.15g  %21.15g\n",ii,jj,faradaytest[j][k],Fcov[j][k]);
      }
    }
    if(fabs(faradaytest[0][3])>1E-10) dualfprintf(fail_file,"1 Fcov=%21.15g faraday=%21.15g\n",Fcov[0][3],faradaytest[0][3]);
#endif
 
  }


  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}



// assumes normal field in pr
int init_vpot(int ii, int jj,FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  FTYPE mcov[NDIM],mcon[NDIM],kcov[NDIM],kcon[NDIM];

  if(PROBLEMTYPE!=2) return(0); // otherwise setup poloidal components using vector potential

  //  return(0);

  coord(ii, jj, CORN, X);
  bl_coord(X, &r, &th);
  get_geometry(ii,jj,CORN,&geom);


  mcon[TT]=0;
  mcon[RR]=0;
  mcon[TH]=0;
  mcon[PH]=1.0;

  kcon[TT]=1.0;
  kcon[RR]=0;
  kcon[TH]=0;
  kcon[PH]=0;

  lower(mcon,&geom,mcov);
  lower(kcon,&geom,kcov);


  // A_\phi
  *A = -B0*(mcov[PH]+2.0*a*kcov[PH]);



  return(0);

}

int init_vpot2field(SFTYPE A[][N2+1],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(SFTYPE A[][N2+1],FTYPE p[][N2M][NPR]);

  //  return(0);

#if(PROBLEMTYPE==2)
  return(vpot2field(A,pr));
#endif
  // otherwise assume didn't need vector potential formulation

  return(0);
}


// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][NPR])
{

  return(0);

}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][NPR])
{


  return(0);
}

  /* this is a little test calculation with a radial field, designed to 
     make the scheme fail */

  /*     Br0 = 1.0 ; ZLOOP { GSET(i,j,CENT) p[i][j][B1] = Br0/(rcurr*rcurr) 
     ; p[i][j][B2] = 0. ; p[i][j][B3] = 0. ; } */


#undef SLOWFAC



#define GAMMIEDERIVATIVE 0
#define NUMREC 1

//#define FCOVDERTYPE NUMREC 
#define FCOVDERTYPE GAMMIEDERIVATIVE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define DXDELTA 1E-5
#elif(REALTYPE==LONGDOUBLETYPE)
#define DXDELTA 1E-6
#endif

void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM])
{
  int j,k,l;
  FTYPE Xhk[NDIM], Xlk[NDIM];
  FTYPE Xhj[NDIM], Xlj[NDIM];
  FTYPE mcovhj,mcovlj,kcovhj,kcovlj;
  FTYPE mcovhk,mcovlk,kcovhk,kcovlk;
  FTYPE mcov_func_mcoord(FTYPE* X, int i, int j); // i not used
  FTYPE kcov_func_mcoord(FTYPE* X, int i, int j); // i not used
  extern FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);

  if(FCOVDERTYPE==GAMMIEDERIVATIVE){

    for(k=0;k<NDIM;k++){
      for(j=0;j<NDIM;j++){

	  for(l=0;l<NDIM;l++) Xlk[l]=Xhk[l]=Xlj[l]=Xhj[l]=X[l]; // location of derivative
	  Xhk[k]+=DXDELTA; // shift up
	  Xlk[k]-=DXDELTA; // shift down

	  Xhj[j]+=DXDELTA; // shift up
	  Xlj[j]-=DXDELTA; // shift down

	  //	  dualfprintf(fail_file,"got here1: k=%d j=%d\n",k,j);

	  
	  mcovhj=mcov_func_mcoord(Xhk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.1: k=%d j=%d\n",k,j);
	  mcovlj=mcov_func_mcoord(Xlk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.2: k=%d j=%d\n",k,j);
	  mcovhk=mcov_func_mcoord(Xhj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.3: k=%d j=%d\n",k,j);
	  mcovlk=mcov_func_mcoord(Xlj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.4: k=%d j=%d\n",k,j);

	  kcovhj=kcov_func_mcoord(Xhk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.5: k=%d j=%d\n",k,j);
	  kcovlj=kcov_func_mcoord(Xlk,0,j); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.6: k=%d j=%d\n",k,j);
	  kcovhk=kcov_func_mcoord(Xhj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.7: k=%d j=%d\n",k,j);
	  kcovlk=kcov_func_mcoord(Xlj,0,k); // 0 not used
	  //	  dualfprintf(fail_file,"got here1.8: k=%d j=%d\n",k,j);

	  //	  dualfprintf(fail_file,"got here2\n");

	  Fcov[j][k] = B0*(
	    +(mcovhj - mcovlj) / (Xhk[k] - Xlk[k])
	    -(mcovhk - mcovlk) / (Xhj[j] - Xlj[j])
	    +2.0*a*(
		   +(kcovhj - kcovlj) / (Xhk[k] - Xlk[k])
		   -(kcovhk - kcovlk) / (Xhj[j] - Xlj[j])
		   )
	    );
      }// j
    }// k
  }
  else if(FCOVDERTYPE==NUMREC){

    for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
      // 0 in dfridr not used
      Fcov[j][k]=B0*(
		      +dfridr(mcov_func_mcoord,X,0,j,k)
		      -dfridr(mcov_func_mcoord,X,0,k,j)
		      +2.0*a*(+dfridr(kcov_func_mcoord,X,0,j,k)
			      -dfridr(kcov_func_mcoord,X,0,k,j)
			      )
		      );
    }

  }
}

//#undef GAMMIEDERIVATIVE
//#undef NUMREC
//#undef FCOVDERTYPE
//#undef DXDELTA




// returns MCOORD m_\mu form of m^\mu={0,0,0,1} value for jth element
FTYPE mcov_func_mcoord(FTYPE* X, int ii, int jj) // i not used
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE mcon[NDIM];
  FTYPE mcov[NDIM];
  struct of_geom geom;
  int i,j,k;

  //  dualfprintf(fail_file,"got here3.1: %d %d\n",ii,jj);
  gcov_func(1,MCOORD,X,gcovmcoord);
  
  //  dualfprintf(fail_file,"got here3.2\n");
  DLOOP gengcov[j][k]=gcovmcoord[j][k];
  geom.gcov=gengcov;

  //  dualfprintf(fail_file,"got here3.3\n");
  mcon[TT]=0.0;
  mcon[RR]=0.0;
  mcon[TH]=0.0;
  mcon[PH]=1.0;
  //  dualfprintf(fail_file,"got here3.4\n");

  // lower only needs geom->gcov
  lower(mcon,&geom,mcov);
  //  dualfprintf(fail_file,"got here3.5\n");

  return(mcov[jj]);
}

// returns MCOORD k_\mu form of k^\mu={1,0,0,0} value for jth element
FTYPE kcov_func_mcoord(FTYPE* X, int ii, int jj) // i not used
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE kcon[NDIM];
  FTYPE kcov[NDIM];
  struct of_geom geom;
  int i,j,k;

  gcov_func(1,MCOORD,X,gcovmcoord);
  
  DLOOP gengcov[j][k]=gcovmcoord[j][k];
  geom.gcov=gengcov;

  kcon[TT]=1.0;
  kcon[RR]=0.0;
  kcon[TH]=0.0;
  kcon[PH]=0.0;

  // lower only needs geom->gcov
  lower(kcon,&geom,kcov);

  return(kcov[jj]);
}

void Jcon_numerical(FTYPE *X, FTYPE *Jcon)
{
  int j,k,l;
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE Fconh,Fconl;
  FTYPE Fcon_func_mcoord(FTYPE* X, int i, int j);
  extern FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);


  if(FCOVDERTYPE==GAMMIEDERIVATIVE){

    for(k=0;k<NDIM;k++){
      

      
      
      
      Jcon[k] = 0;
      for(j=0;j<NDIM;j++){

	for(l=0;l<NDIM;l++) Xl[l]=Xh[l]=X[l]; // location of derivative
	Xh[j]+=DXDELTA; // shift up
	Xl[j]-=DXDELTA; // shift down

	// F^{kj}
	Fconh=Fcon_func_mcoord(Xh,k,j);
	Fconl=Fcon_func_mcoord(Xl,k,j);
	
	// J^\mu = {F^{\mu\nu}}_{;\nu}
	// \detg J^k = F^{kj}_{;j} = (\detg F^{kj})_{,j} <--- thing actually computed and returned
	Jcon[k] += (Fconh - Fconl) / (Xh[j] - Xl[j]);
      } // j
    }// k
  }
  else if(FCOVDERTYPE==NUMREC){

    for(k=0;k<NDIM;k++){
      Jcon[k] = 0;
      for(j=0;j<NDIM;j++){
	Jcon[k]+=dfridr(Fcon_func_mcoord,X,0,j,k);
      }
    }

  }
}

#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef FCOVDERTYPE
#undef DXDELTA




// returns MCOORD F^{ii jj}
FTYPE Fcon_func_mcoord(FTYPE* X, int ii, int jj)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE Fcov[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM];
  struct of_geom geom;
  int i,j,k;
  void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM]);
  extern void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM]);

  gcov_func(1,MCOORD,X,gcovmcoord);
  
  DLOOP gengcov[j][k]=gcovmcoord[j][k];

  matrix_inverse(gengcov,gengcon);

  geom.gcov=gengcov;
  geom.gcon=gengcon;


  Fcov_numerical(X, Fcov);
  /* get covariant Maxwell from contravariant */
  raise_A(Fcov, &geom, Fcon) ;

  return(gdet_func(gcovmcoord)*Fcon[ii][jj]);
}
