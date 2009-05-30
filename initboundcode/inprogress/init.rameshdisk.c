/*

 to use this, need to probably change following elsewhere:

1) AVOIDCS in phys.ffde.c
2) EOMTYPE in global.h
3) MCOORD KSCOORDS in metric.h
4) probably defcoord==0's entries in coord.c
5) maybe show boundary zones by modifying dump.c
6) ZEROOUTFLOWFLUX 1 in step_ch.h ?


*/

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

static SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

static FTYPE mm;

#define NTHETAMAX 10000

static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac,mydThetac;
static int NTHETA;


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  Omegastar=0;

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
  
  
#if(MCOORD==CYLMINKMETRIC)

  //  Zin=-1E3;
  //  Zout=1E3;
  //  Rin=1.0;
  //  Rout=1E3;

  Zin=-100;
  Zout=100;
  Rin=22;
  Rout=1E3;


  defcoord=666; // see coord.c for npow parameter

#else
  // metric stuff first
  a = 0.0 ;

  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = -1.0;
  Rhor=rhor_calc(0);
  Rout = 1E3;
 
  //  Rin=setRin(setihor());
  Rin = 0.98 * Rhor;
  
  hslope = 0.1;

  // define coordinate type
  defcoord = 9;
#endif

  return(0);
}

int init_global(void)
{

  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

#if(MCOORD==CYLMINKMETRIC)

  BCtype[X1UP]=FIXED;
  BCtype[X1DN]=FIXED;
  BCtype[X2UP]=FIXED;
  BCtype[X2DN]=FIXED;


#else

#if(1)
  BCtype[X1UP]=FIXED;
  BCtype[X1DN]=FIXED;
#else
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
#endif

  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

#endif

  /* output choices */
  tf = 1.2*Rout;

  DTd = tf/100.0;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;			/* logfile frequency, in units of M */
  DTi = tf/100.0;			/* image file frequ., in units of M */
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
  int j,k;

  void Fcov_numerical(FTYPE *X, FTYPE (*Fcov)[NDIM]);
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  extern int EBtopr(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int EBtopr_2(FTYPE *E,FTYPE *B,struct of_geom *geom, FTYPE *pr);
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);

  struct of_state q;
  FTYPE faradaytest[NDIM][NDIM];


  if(EOMTYPE!=EOMFFDE){
    dualfprintf(fail_file,"Are you sure?\n");
    myexit(1);
  }
  

  coord(ii, jj, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(ii, jj, CENT, &geom); 

  // 0 : split monopole
  // 1 : monopole
  // 2 : Wald
  // 3 : BZ paraboloidal
  // 4 : GRMHD nearly-paraboloidal
  // 5 : BZ para but monopole near horizon
  // 6 : GRMHD nearly-paraboloidal but monopole near horizon
  // 7 : constant velocity Ramesh disk // can do cylindrical coordinatse also
#define PROBLEMTYPE 7
#define B0 1.0

  pr[RHO]=pr[UU]=0.0;
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

    //    DLOOPA dualfprintf(fail_file,"Econ[%d]=%21.15g Bcon[%d]=%21.15g\n",j,Econ[j],j,Bcon[j]);

    EBtopr(Econ,Bcon,&geom,pr);
    //    dualfprintf(fail_file,"EBtopr\n");

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
  else if(PROBLEMTYPE==3){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==4){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==5){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==6){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }
  else if(PROBLEMTYPE==7){
    pr[B1]=0.0; // not used   
    pr[B2]=0.0;
    pr[B3]=0.0;
  }




  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}



//#define NTHETAMAX 10000


// assumes normal field in pr
int init_vpot(int ii, int jj,FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  FTYPE Xc[NDIM],rc,thc;
  struct of_geom geom;
  struct of_geom geomc;
  FTYPE mcov[NDIM],mcon[NDIM],kcov[NDIM],kcon[NDIM];
  FTYPE th2;
  FILE * inTheta;
  //  static FTYPE Thetavstheta[NTHETAMAX],dThetadthetavstheta[NTHETAMAX],theta[NTHETAMAX],dtheta,myfloati,myTheta,myThetac;
  int dumi;
  FTYPE ftemp,ftemp1,ftemp2;
  int i;
  static int firsttime=1;
  FTYPE rprime;
  FTYPE B0para;
  FTYPE rtrans;
  FTYPE nu,ss,ucrit;
  //  int NTHETA;
  FTYPE Ac;
  FTYPE theu;

  if(PROBLEMTYPE<=1) return(0); // otherwise setup poloidal components using vector potential


  coord(ii, jj, CORN, X);
  bl_coord(X, &r, &th);
  get_geometry(ii,jj,CORN,&geom);


  coord(ii, jj, CENT, Xc);
  bl_coord(Xc, &rc, &thc);
  get_geometry(ii,jj,CENT,&geomc);


  if(PROBLEMTYPE==2){

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
  }
  else if((PROBLEMTYPE==3)||(PROBLEMTYPE==5)){
    if(th<=M_PI*0.5){
      th2=th;
    }
    else{
      th2=M_PI-th;
    }

    // BZ paraboloidal
    if(PROBLEMTYPE==3) rprime=r;
    // BZ paraboloidal with monopole near horizon
    else if(PROBLEMTYPE==5) rprime=sqrt(pow(r,2)+pow(2.0*Rhor,2));

    // setup so A=0 on poles
    // B0para setup so B0 is value of aphi on horizon
    B0para=2.0*B0/(Rhor+4.0*log(2.0)-2.0);
    *A = 0.5*B0para*(rprime*(1.0-cos(th2))+2.0*(1.0+cos(th2))*(1.0-log(1.0+cos(th2))) - 4.0*(1.0-log(2.0)) );

  }
  else if((PROBLEMTYPE==4)||(PROBLEMTYPE==6)){
    // GRMHD nearly-paraboloidal

    // for PROBLEMTYPE==4
    NTHETA=10000;
#define INDEXN (1.0/4.0) // requires Theta to be l=-n

    if(firsttime){
      firsttime=0;

      if(myid==0){
	if( (inTheta=fopen("myThetavstheta.txt","rt"))==NULL){
	  dualfprintf(fail_file,"Cannot open myThetavstheta.txt\n");
	  myexit(100);
	}

	// read in data file
	while(fgetc(inTheta)!='\n'); // skip first line, a comment
	for(i=0;i<NTHETA;i++){
	  fscanf(inTheta,"%d %lf %lf %lf\n",&dumi,&theta[i],&Thetavstheta[i],&dThetadthetavstheta[i]);
	}

      }

    // send data to all CPUs
#if(USEMPI)
      MPI_Bcast(&NTHETA,1,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&theta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Thetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&dThetadthetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif
    }

    
    dtheta=theta[1]-theta[0]; // uniform grid


    if(th<0.0){
      // function is even
      th2=-th;
    }
    else if(th<=M_PI*0.5){
      // theta is itself
      th2=th;
    }
    else if(th<=M_PI){
      th2=M_PI-th;
    }
    else{
      // function is even
      th2=-(M_PI-th);
    }

    // linearly interpolate \Theta using \theta
    myfloati=(th2-0.0)/dtheta;

    myTheta=Thetavstheta[(int)myfloati]
      +(Thetavstheta[(int)myfloati+1]-Thetavstheta[(int)myfloati])/(theta[(int)myfloati+1]-theta[(int)myfloati])*(th2-theta[(int)myfloati]);

    rtrans=Rhor;

    // n=1/4 strict
    if(PROBLEMTYPE==4) rprime=r;
    // n=1/4 but relaxes to n=1 near horizon
    else if(PROBLEMTYPE==6) rprime=sqrt(pow(r,2)+pow(rtrans,2));

    *A = 0.5*B0*pow(rprime,-INDEXN+1.0)*myTheta;

  }
  else if(PROBLEMTYPE==7){
    // Ramesh constant velocity disk

    if(firsttime){
      firsttime=0;

      if(myid==0){
	if( (inTheta=fopen("nu.75_m.5.txt","rt"))==NULL){
	  dualfprintf(fail_file,"Cannot open nu.75_m.5.txt\n");
	  myexit(100);
	}
	
	// skip first blank line
	//while(fgetc(inTheta)!='\n'); // skip first line, a comment

	// read in data file

	// read in first row (parameters)
	fscanf(inTheta,"%lf %lf %lf %lf",&nu,&mm,&ss,&ucrit);
	//dualfprintf(fail_file,"got0: %g %g %g %g\n",nu,mm,ss,ucrit);
	trifprintf("Ramesh parameters: %g %g %g %g\n",nu,mm,ss,ucrit);
	
	// for nu.75_m.5.txt
	NTHETA=1314;

	// for nu.75_m.0005.txt
	//	NTHETA=1206;

	for(i=0;i<NTHETA;i++){
	  //	while(!feof(inTheta)){
	  fscanf(inTheta,"%lf %lf %lf",&ftemp,&ftemp1,&ftemp2);
	  // data in reverse order of theta
#if(MCOORD==CYLMINKMETRIC)
	  // Z/R itself
	  theta[i]=ftemp; // true theta since file has ftemp=u= Z/R=tan(\theta)
	  Thetavstheta[i]=ftemp1;
	  dThetadthetavstheta[i]=ftemp2;
#else
	  // theta
	  theta[NTHETA-1-i]=atan(1.0/ftemp); // true theta since file has ftemp=u= Z/R=tan(\theta)
	  Thetavstheta[NTHETA-1-i]=ftemp1;
	  dThetadthetavstheta[NTHETA-1-i]=ftemp2;
#endif
	  //	  dualfprintf(fail_file,"got1: %d %g %g %g\n",i,theta[NTHETA-1-i],Thetavstheta[NTHETA-1-i],dThetadthetavstheta[NTHETA-1-i]);
	  //	  i++;
	}
	//	NTHETA=i;

      }

    // send data to all CPUs
#if(USEMPI)
      MPI_Bcast(&nu,1,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&mm,1,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&ss,1,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&ucrit,1,MPI_FTYPE,0,MPI_COMM_WORLD);


      MPI_Bcast(&NTHETA,1,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&theta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Thetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
      MPI_Bcast(&dThetadthetavstheta,NTHETA,MPI_FTYPE,0,MPI_COMM_WORLD);
#endif
    }



    
#if(MCOORD==CYLMINKMETRIC)
    th2=fabs(th)/r; // Z/R really, but solution symmetric around Z=0
    // R
    rprime=r;
#else
    if(th<0.0){
      // function is even
      th2=-th;
    }
    else if(th<=M_PI*0.5){
      // theta is itself
      th2=th;
    }
    else if(th<=M_PI){
      th2=M_PI-th;
    }
    else{
      // function is even
      th2=-(M_PI-th);
    }
    // R
    rprime=r*sin(th2);
#endif


    for(i=0;i<NTHETA;i++){
      if(theta[i]>th2) break;
    }
    // linearly interpolate \Theta using \theta
    // \theta is not uniform, so loop to find first theta

    
    myTheta=Thetavstheta[i-1] + (Thetavstheta[i]-Thetavstheta[i-1])/(theta[i]-theta[i-1])*(th2-theta[i-1]);    


    
    // A_\phi = P(R,Z) = R^\nu T(u)
    *A = B0*pow(rprime,nu)*myTheta;




    //////////////////////////////
    //
    // set analytic B^\phi.  With this choice and v^i chosen by stationarity, the analytic solution is set so FIXED boundaries are correct
    //
    ///////////////////////////////
#if(MCOORD==CYLMINKMETRIC)
    th2=fabs(thc)/rc; // Z/R really, but solution symmetric around Z=0
    // R
    rprime=rc;
    theu=fabs(th2); // Z/R
#else
    if(thc<0.0){
      // function is even
      th2=-thc;
    }
    else if(thc<=M_PI*0.5){
      // theta is itself
      th2=thc;
    }
    else if(thc<=M_PI){
      th2=M_PI-thc;
    }
    else{
      // function is even
      th2=-(M_PI-thc);
    }
    // R
    // bit different than how setting vector potential, but maybe A_\phi is even and B^\phi not, across polar axes?
    rprime=rc*sin(thc);
    theu=fabs(cos(thc)/sin(thc)); // Z/R
#endif

    // linearly interpolate \Theta using \theta
    // \theta is not uniform, so loop to find first theta

    for(i=0;i<NTHETA;i++){
      if(theta[i]>th2) break;
    }

    // centered myTheta    
    myThetac=Thetavstheta[i-1] + (Thetavstheta[i]-Thetavstheta[i-1])/(theta[i]-theta[i-1])*(th2-theta[i-1]);

    // centered mydTheta    
    mydThetac=dThetadthetavstheta[i-1] + (dThetadthetavstheta[i]-dThetadthetavstheta[i-1])/(theta[i]-theta[i-1])*(th2-theta[i-1]);



    // GODMARK: in this code, A is over same range as p, so ok to do this, else need check on i,j,k
    
    // in orthonormal basis
    p[ii][jj][B3]=B0*(-nu*ss*pow(rprime,nu-2.0)*pow(myThetac,(nu-1.0)/nu));
    // transform to coordinate basis for diagonal metric (at least \phi not connected to any other coordinate -- always true so far in HARM)
    p[ii][jj][B3]/=sqrt(geomc.gcov[PH][PH]);

    // just do directly (not general)
    // original Ramesh formula is super divergent!?
    //    p[ii][jj][B3]=B0*(-nu*ss*pow(rprime,nu-3.0)*pow(myTheta,(nu-1.0)/nu));



    // B3
    // below not general
    //    p[ii][jj][B3]=B0*(-nu*ss*pow(rprime,nu-3.0)*pow(myThetac,(nu-1.0)/nu));


    if(p[ii][jj][B3]>1E30) p[ii][jj][B3]=1E30;
    else if(p[ii][jj][B3]<-1E30) p[ii][jj][B3]=-1E30;

    // B1
    // orthonormal
    p[ii][jj][B1]=B0*(-pow(rprime,nu-2.0)*mydThetac);
    // coordinate
    p[ii][jj][B1]/=sqrt(geomc.gcov[1][1]);
    if(p[ii][jj][B1]>1E30) p[ii][jj][B1]=1E30;
    else if(p[ii][jj][B1]<-1E30) p[ii][jj][B1]=-1E30;

    // B2
    p[ii][jj][B2]=B0*pow(rprime,nu-2.0)*(nu*myThetac  - theu*mydThetac);
    p[ii][jj][B2]/=sqrt(geomc.gcov[2][2]);
    if(p[ii][jj][B2]>1E30) p[ii][jj][B2]=1E30;
    else if(p[ii][jj][B2]<-1E30) p[ii][jj][B2]=-1E30;




#if(MCOORD==CYLMINKMETRIC)
    // B^\phi is antisymmetric across equator
    // thc=Z/R.  Normal range is Z>0
    if(thc<0.0){
      p[ii][jj][B3]*=-1.0;
      p[ii][jj][B1]*=-1.0;
    }
#else
    // B^\phi is antisymmetric across equator
    if(thc>M_PI*0.5){
      p[ii][jj][B3]*=-1.0;
      p[ii][jj][B1]*=-1.0;
    }
#endif
    

    ////////////////
    //
    // need omegaf analytic! (only used at center)
    //
    ////////////////
    Ac=B0*pow(rprime,nu)*myThetac;

#if(MCOORD==CYLMINKMETRIC)
    p[ii][jj][UU]=omegafanalytic[ii][jj][0]=mm*pow(Ac/Thetavstheta[0],-1.0/nu);
    //omegafanalytic[ii][jj][0]=mm*pow(Ac/Thetavstheta[0],-1.0/nu);
    // test below
    //omegafanalytic[ii][jj][0]=mm*pow(B0*pow(rprime,nu)*1E-8,-1.0/nu);
#else
    omegafanalytic[ii][jj][0]=mm*pow(Ac/Thetavstheta[NTHETA-1],-1.0/nu);
#endif

     
    //     dualfprintf(fail_file,"%d %d : th2=%21.15g : theta[%d]=%21.15g myTheta=%21.15g A=%21.15g\n",ii,jj,th2,i,theta[i],myTheta,*A);

  }




  return(0);

}

int init_vpot2field(SFTYPE A[][N2M],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(SFTYPE A[][N2M],FTYPE p[][N2M][NPR]);

#if( PROBLEMTYPE<=1 )
  //#if( (PROBLEMTYPE<=1)||(PROBLEMTYPE==7) ) // GODMARK : test of errors
  return(0);
  // didn't need vector potential formulation
#else
  return(vpot2field(A,pr));
#endif
}



int init_postfield(FTYPE pr[][N2M][NPR])
{
  FTYPE Bccon[NDIM];
  FTYPE prnew[NPR];
  int i,j;
  extern int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  struct of_geom geom;
  FTYPE r,th,X[NDIM];
  FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th);


  ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){

    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
    get_geometry(i, j, CENT, &geom);


    Bccon[0]=0;
    Bccon[1]=pr[i][j][B1];
    Bccon[2]=pr[i][j][B2];
    Bccon[3]=pr[i][j][B3];

    Omegastar=get_omegastar(&geom,r,th);

    if(OBtopr(Omegastar,Bccon,&geom,prnew)>=1){
      dualfprintf(fail_file, "OBtopr(vcon2pr): space-like error in init_postfield()\n");
      dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
      failed=1;
      return(1);
      
    }
    
    pr[i][j][U1]=prnew[U1];
    pr[i][j][U2]=prnew[U2];
    pr[i][j][U3]=prnew[U3];
    
  }

  return(0);

}

FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th)
{
  FTYPE omegak,omegaf,omegah,omegadiskbh;
  FTYPE ftemp1,ftemp2;
  FTYPE rtrans,ntrans;
  FTYPE CON,CON2;
  FTYPE omegatrans;
  FTYPE vtrans;



#if(PROBLEMTYPE<7)

  omegak=1.0/(pow(r,3.0/2.0)+a);
  omegah=(a/(2.0*Rhor));

  //  omegadiskbh=0.6*omegah;
  // does quite well till about dump 11, then has problems and jump to above equator is omegaf2/omegah~0.27.  Leads to crispy omegaf by dump=20

  //  omegadiskbh=0.2652*omegah; // goes uu0~30 before dump 1
  // above chooses sharp changed to omegaf~1 just outside equator by 2nd dump
  // by dump=8 or 11, polar region has settled to standard form with 0.65*omegah just above equator

  //  omegadiskbh=1.0*omegah; // lasts till dump 8 before uu0~30 , but close to disk omegaf~0.27*omegah after settles to non-force-free

  omegadiskbh=0.4*omegah;

#if(0)
  // Just Keplerian for all radii
  omegaf=omegak;
#elif(0)
  if(r>Risco){
    omegaf=omegak;
  }
  else{
    omegaf=1.0/(pow(Risco,3.0/2.0)+a);
  }
#elif(0)
  omegaf = omegak/(1.0 + pow(Rhor/r,2.0));
#elif(0)
  // goes from \Omega_K at large radii to (omegadiskbh) on horizon
  if(r>2.0*Rhor){
    ftemp1=omegak*(1.0-pow(2.0*Rhor/r,3.0));
    ftemp2=(omegadiskbh)*pow(2.0*Rhor/r,3.0);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math1
  rtrans=Rhor;
  ntrans=6.21;
  omegadiskbh=(0.5/0.85)*omegah;

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math2
  rtrans=Rhor;
  ntrans=10.744;
  omegadiskbh=(0.34)*omegah;

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math3
  rtrans=2.0*Rhor;
  //  omegadiskbh=(0.34)*omegah; // works ok
  //  omegadiskbh=(0.315)*omegah;
  omegadiskbh=(0.27)*omegah;
  //  omegadiskbh=(0.4)*omegah;
  //  omegadiskbh=(0.32)*omegah;
  ntrans=3.0/(omegadiskbh*(a+pow(rtrans,3.0/2.0)));

  if(r>rtrans){
    ftemp1=omegak*(1.0-pow(rtrans/r,3.0));
    ftemp2=(omegadiskbh)*pow(rtrans/r,ntrans);
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh;
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see omegaf_vs_omegah_bz77_grffe_grmhd.nb
  // math4
  rtrans=2.0*Rhor;
  omegadiskbh=(0.5/0.85)*omegah;
  CON=-0.0479309;
  CON2=0.0612062;

  if(r>rtrans){
    ftemp1=omegak;
    ftemp2=0;
  }
  else{
    ftemp1=0;
    ftemp2=omegadiskbh+CON*pow(r-Rhor,2.0)+CON2*(r-Rhor);
  }

  omegaf = ftemp1 + ftemp2;
#elif(0)
  // see para_bz_omegaf_compute_plot.nb
  omegadiskbh=0.27*omegah;
  rtrans=2.0*Rhor;
  CON=-3.0*omegadiskbh*pow(Rhor - rtrans,-2); // B
  CON2=2.0*omegadiskbh*pow(-Rhor + rtrans,-3); // C

  if(r>rtrans) omegaf=0;
  else if(r<Rhor) omegaf=omegadiskbh;
  else omegaf=omegadiskbh + CON*(r-Rhor)*(r-Rhor) + CON2*(r-Rhor)*(r-Rhor)*(r-Rhor);

#elif(1)

  // see para_bz_omegaf_compute_plot.nb
  omegadiskbh=0.27*omegah;
  rtrans=2.0*Rhor;
  vtrans=0.5;
  //  omegatrans=(vtrans/rtrans)/(sqrt(geom->gcov[3][3])/rtrans);
  omegatrans=(vtrans/rtrans)/(r/rtrans);

  // fsB
  CON = (-3.0*omegadiskbh*rtrans*rtrans - (Rhor - 4.0*rtrans)*vtrans)/( pow(rtrans-Rhor,2)*rtrans*rtrans );

  // fsC
  CON2 =(+2.0*omegadiskbh*rtrans*rtrans + (Rhor - 3.0*rtrans)*vtrans)/( pow(rtrans-Rhor,3)*rtrans*rtrans );

    
  if(r>=rtrans) omegaf=omegatrans;
  else if(r<=Rhor) omegaf=omegadiskbh;
  else omegaf=omegadiskbh + CON*(r-Rhor)*(r-Rhor) + CON2*(r-Rhor)*(r-Rhor)*(r-Rhor);

#endif


#elif(PROBLEMTYPE==7)

  // mm is global variable with local scope to this file

  omegaf=omegafanalytic[geom->i][geom->j][0];
  // only true at equator
  //#if(MCOORD==CYLMINKMETRIC)
  //  omegaf=mm/r;
  //#else
  //  omegaf=mm/(r*sin(th));
  //#endif

  //  dualfprintf(fail_file,"omegaf=%21.15g mm=%21.15g\n",omegaf,mm);

#endif


  return(omegaf);



  // NON-ROTATING DISK
  //  return(0.0);

}

// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][NPR])
{

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

  if(PROBLEMTYPE==7) return(0); // don't worry about normalization

  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, CENT, &geom);    
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
     
    if((r>Rhor)&&(fabs(th-M_PI*0.5)<0.1)){
      if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
	FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  norm = sqrt((B0*B0)/bsq_max);
  
  bsq_max = 0.;
  ZLOOP {
    p[i][j][B1] *= norm;
    p[i][j][B2] *= norm;
    p[i][j][B3] *= norm;

    get_geometry(i, j, CENT, &geom);
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    if((r>Rhor)&&(fabs(th-M_PI*0.5)<0.1)){
      if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
	FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
    
  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);

  


  return(0);
}

  /* this is a little test calculation with a radial field, designed to 
     make the scheme fail */

  /*     Br0 = 1.0 ; ZLOOP { GSET(i,j,CENT) p[i][j][B1] = Br0/(rcurr*rcurr) 
     ; p[i][j][B2] = 0. ; p[i][j][B3] = 0. ; } */


#undef SLOWFAC



#define GAMMIEDERIVATIVE 0
#define NUMREC 1

//#define DXDERTYPE NUMREC 
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

#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef FCOVDERTYPE
#undef DXDELTA




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
