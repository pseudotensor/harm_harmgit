#define JET6LIKEUSERCOORD 0
#define UNIHALFUSERCOORD 1

#define WHICHUSERCOORD JET6LIKEUSERCOORD


// for defcoord=JET6COORDS like USERCOORDS
static FTYPE npow,r1jet,njet1,njet,r0jet,rsjet,Qjet, ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3, rs, r0,npow2,cpow2,rbr,x1br, h0,cpow3; 


void set_coord_parms_nodeps_user(int defcoordlocal)
{
  if(1){
    // see jet3coords_checknew.nb
    npow=1.0;

    /////////////////////
    // RADIAL GRID SETUP
    /////////////////////
    npow=1.0;  //don't change it, essentially equivalent to changing cpow2

    //radial hyperexponential grid

    //power exponent
    npow2=6.0; // WALD: 6.0->4.0

    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    //    cpow3=0.01;
    cpow3=1.0;
    //radius at which hyperexponentiation kicks in
    //    rbr = 1E3;
    rbr = 5E2; // WALD 5E2->5E7



    // must be same as in dxdxp()
    // GODMARK: Note njet here is overwritten by njet later, but could have been different values if setup variable names differently.
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(1){
      r1jet=2.8;
      njet=0.3;
      r0jet=15.0;
      rsjet=40.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
    }

    // for switches from normal theta to ramesh theta
    rs=40.0; // shift
    r0=20.0; // divisor
 
    // for theta1
    //    hslope=0.3 ; // resolve inner-radial region near equator
    r0jet3=20.0; // divisor
    rsjet3=0.0; // subtractor

    // for theta2
    h0=0.3; // inner-radial "hslope" for theta2
    //h0=0.1; // inner-radial "hslope" for theta2 // for thinner disks, change this.
    // GODMARK: Note that this overwrites above njet!
    // power \theta_j \propto r^{-njet}
    njet=1.0;


    // see fix_3dpoledtissue.nb
#if(0)
    ntheta=21.0;
    htheta=0.15;
    rsjet2=5.0;
    r0jet2=2.0;
#else
    ntheta=5.0;
    htheta=0.15;
    rsjet2=5.0;
    r0jet2=2.0;
#endif

  }

}


void set_coord_parms_deps_user(int defcoordlocal)
{
  if(1){
    x1br = log( rbr - R0 ) / npow;  //the corresponding X[1] value
  }

}

void write_coord_parms_user(int defcoordlocal, FILE *out)
{
  if(1){
    fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,r0jet,rsjet,Qjet,ntheta,htheta,rsjet2,r0jet2,rsjet3,r0jet3,rs,r0,npow2,cpow2,rbr,x1br,cpow3);
  }

}
void read_coord_parms_user(int defcoordlocal, FILE *in)
{

  if(1){
    fscanf(in,HEADER19IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet,&ntheta,&htheta,&rsjet2,&r0jet2,&rsjet3,&r0jet3,&rs,&r0,&npow2,&cpow2,&rbr,&x1br,&cpow3);

  }
}

void read_coord_parms_mpi_user(int defcoordlocal)
{
  if(1){
#if(USEMPI)
    MPI_Bcast(&npow, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&ntheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&htheta, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rsjet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0jet3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rs, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&r0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&npow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow2, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&rbr, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&x1br, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&cpow3, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
#endif
  }

}


void blcoord_user(FTYPE *X, FTYPE *V)
{
  extern FTYPE mysin(FTYPE th);

  if(1){

#if(0) // no change in exponentiation
    // JET3COORDS-like radial grid
    V[1] = R0+exp(pow(X[1],npow)) ;
#elif(WHICHUSERCOORD==UNIHALFUSERCOORD)

    Rout=2000.0;
    theexp = npow*X[1];
    npow=1.0;
    FTYPE gconst1=1.0;
    FTYPE gconst2=gconst1*.000001;
    V[1] = R0 + gconst1*X[1] + gconst2*exp(theexp);


#elif(WHICHUSERCOORD==JET6LIKEUSERCOORD)

    FTYPE theexp = npow*X[1];
    if( X[1] > x1br ) {
      theexp += cpow2 * pow(X[1]-x1br,npow2);
    }
    V[1] = R0+exp(theexp);


    //    FTYPE npowtrue,npowlarger=10.0;
    //    FTYPE npowrs=1E3;
    //    FTYPE npowr0=2E2;
    //    npowtrue = npow + (npowlarger-npow)*(0.5+1.0/M_PI*atan((V[1]-npowrs)/npowr0));
    //    V[1] = R0+exp(pow(X[1],npowtrue)) ;
#elif(0)
    // avoid jump in grid at rbr
    // determine switches
    FTYPE r0rbr=rbr/2.0;
    FTYPE switch0 = 0.5+1.0/M_PI*atan((V[1]-rbr)/r0rbr); // 1 for outer r

    FTYPE V1 = R0+exp(npow*X[1]);
    FTYPE V2 = R0+exp(npow*X[1] + cpow2 * pow(cpow3*(X[1]-x1br*1.0),npow2));

    V[1] = V1*(1.0-switch0) + V2*switch0;

#endif



    FTYPE theta1,theta2,arctan2;


#if(0)
    // JET3COORDS-based:
    FTYPE myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    theta1 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);
#else
    // RAMESH BASED
    // myhslope here is h2 in MCAF paper
    //    // h0 here is h3 in MCAF paper
    //FTYPE njetvsr;
    //if(V[1]<rbr) njetvsr=njet;
    //    else njetvsr=njet/(V[1])*rbr;
    //else njetvsr=
    //njetvsr=njet;

    FTYPE localrbr=Rout; //500.0; // rbr;
    //    FTYPE localrbr=rbr;
    FTYPE localrbrr0=100.0; //MAX(100.0,localrbr/2.0);

    FTYPE switch0 = 0.5+1.0/M_PI*atan((V[1]-localrbr)/localrbrr0); // large r
    FTYPE switch2 = 1.0-switch0; // small r

    //    switch0=0.0; switch2=1.0;

    FTYPE myhslope1=h0 + pow( (V[1]-rsjet3)/r0jet3 , njet);
    FTYPE myhslope2=h0 + pow( (localrbr-rsjet3)/r0jet3 , njet);
    FTYPE myhslope = myhslope1*switch2 + myhslope2*switch0;

    //    myhslope = pow(pow(myhslope1,-2.0) + pow(myhslope2,-2.0),-0.5);

    myhslope=myhslope1;

    // determine theta2
    FTYPE myx2;
    if(X[2]>1.0) myx2=2.0-X[2];
    else if(X[2]<0.0) myx2=-X[2];
    else myx2=X[2];

    FTYPE th2 = 0.5*M_PI*(1.0 + atan(myhslope*(myx2-0.5))/atan(myhslope*0.5));

    if(X[2]>1.0) th2=2.0*M_PI-th2;
    else if(X[2]<0.0) th2=-th2;

    // determine theta0
    // JET3COORDS-based:
    myhslope1=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0jet-rsjet/r0jet)));
    myhslope2=2.0-Qjet*pow(localrbr/r1jet,-njet*(0.5+1.0/M_PI*atan(localrbr/r0jet-rsjet/r0jet)));
    myhslope = myhslope1*switch2 + myhslope2*switch0;
    // myhslope here is h0 in MCAF paper
    FTYPE th0 = M_PI * X[2] + ((1. - myhslope) * 0.5) * mysin(2. * M_PI * X[2]);


    // determine switches (only function of radius and not x2 or theta)
    switch0 = 0.5+1.0/M_PI*atan((V[1]-rs)/r0); // switch in .nb file // switch0->1 as r->infinity
    switch2 = 1.0-switch0; // for inner radial region

    // this works because all functions are monotonic, so final result is monotonic.  Also, th(x2=1)=Pi and th(x2=0)=0 as required
    theta1 = th0*switch2 + th2*switch0; // th0 is activated for small V[1] and th2 is activated at large radii.  Notice that sum of switch2+switch0=1 so normalization correct.
    //    theta1=th0;
    theta1=th2;

#endif

    if(0){
      // fix_3dpoledtissue.nb based:
      theta2 = M_PI*0.5*(htheta*(2.0*X[2]-1.0)+(1.0-htheta)*pow(2.0*X[2]-1.0,ntheta)+1.0);
      
      // generate interpolation factor
      arctan2 = 0.5 + 1.0/M_PI*(atan( (V[1]-rsjet2)/r0jet2) );
      
      // now interpolate between them
      V[2] = theta2 + arctan2*(theta1-theta2);
    }


    //V[2] = theta1;

    if(1){
      FTYPE fraceq=0.3;
      FTYPE fracpole=(1.0-fraceq)/2.0;
      FTYPE x2p1=0.0+fracpole;
      FTYPE x2p2=1.0-fracpole;
      FTYPE swide=0.04; //1E-1;

      // s(x) = 0.5 + 0.5 tanh((x-a)/b)

      //      FTYPE switchh0 = 0.5+1.0/M_PI*atan((X[2]-x2p1)/swide);
      //      FTYPE switchh2 = 0.5+1.0/M_PI*atan((X[2]-x2p2)/swide);

      FTYPE switchh0 = 0.5+0.5*tanh((X[2]-x2p1)/swide);
      FTYPE switchh2 = 0.5+0.5*tanh((X[2]-x2p2)/swide);

      FTYPE eqh=0.1;
      FTYPE theq = M_PI * X[2] + ((1. - eqh) * 0.5) * mysin(2. * M_PI * X[2]);
    
      FTYPE thup=switchh0*theq + (1.0-switchh0)*theta1;

      V[2]=switchh2*theta1 + (1.0-switchh2)*thup;

    }
    else{
      V[2]=theta1;
    }

    if(0){ // Sam Gralla
      FTYPE transwidth=0.06; //1E-1;
      FTYPE xcent=0.5; // fixed
      FTYPE transR=0.5*(1.0-tanh(+(X[2] - xcent)/transwidth));
      FTYPE transL=0.5*(1.0-tanh(-(X[2] - xcent)/transwidth));

      h0=0.0;
      //      FTYPE wpar=h0 + pow( (V[1]-rsjet3)/r0jet3 , -njet);
      FTYPE wpar=pow( (V[1])/1 , -njet);

      FTYPE line1 = wpar*X[2];
      FTYPE line2 = M_PI + wpar*(X[2]-1.0);

      V[2] = line1*transR + line2*transL;

      //      dualfprintf(fail_file,"X[2]=%g V[2]=%g\n",X[2],V[2]);

    }

    if(1){ // Sam Gralla 2

      h0=0.0;

#define cr(x) (exp(-1.0/(x)))
#define tr(x) (cr(x)/(cr(x) + cr(1.0-(x))))
#define trans(x,L,R) ((x)<=(L) ? 0.0 : ((x)>=(R) ? 1.0 : tr(((x)-(L))/((R)-(L)))) )
#define transR(x,center,width) ( 0.5*(1.0-tanh(+((x)-(center))/(width))))
#define transL(x,center,width) ( 0.5*(1.0-tanh(-((x)-(center))/(width))))
#define transM(x,center,width) ( exp(-pow(((x)-(center))/((width)*0.5),2.0) ) )
#define line1(x,w) ((x)*(w))
#define line2(x,w) ((x)*(w)+M_PI-(w))
#define line3(x,w) ((x)*(w))
      //#define wparsam(x,r) (h0 + pow( ((r)-rsjet3)/r0jet3 , -njet))
      //#define wparsam(x,r) (h0 + pow( ((r)-0.0)/4.2 , -njet))
#define wparsam(x,r) (h0 + pow(0.15 + ((r)-0.0)/10.0 , -njet))
#define thetasam(x,r,w,xp1,xp2) (line1(x,w)*(1.0-trans(x,xp1,xp2)) + line2(x,w)*trans(x,xp1,xp2))

        V[2] = thetasam(X[2],V[1],wparsam(X[2],V[1]),0.25,0.75);
      //      V[2] = thetasam(X[2],V[1],1.0/V[1],0.2,0.8);

      //      dualfprintf(fail_file,"tr=%g %g %g %g\n",tr(0.5),line1(0.5,wparsam(0.5,V[1])),trans(X[2],0.2,0.8),line2(0.5,wparsam(0.5,V[1])));
      

    }

    if(0){ // Sam Gralla 3

      h0=0.0;

#define plateau(x,L,R,W) (trans(x,(L)-0.5*(W),(L)+0.5*(W))*(1.0-trans(x,(R)-0.5*(W),(R)+0.5*(W))))
#define lineeq(x,w) ((x)*(w)+(0.5*M_PI)-(0.5*w))
#define linepole(x,w) (line1(x,w))
#define thetaL(x,wp,weq,xp1,xp2) ( linepole(x,wp)*(1.0-trans(x,xp1,xp2)) + lineeq(x,weq)*trans(x,xp1,xp2) )
#define thetasam2(x,wp,weq,xp1,xp2) ( x<0.5 ? thetaL(x,wp,weq,xp1,xp2) : -thetaL(1.0-x,wp,weq,xp1,xp2)+M_PI )

      FTYPE poleslope=wparsam(X[2],V[1]);
      FTYPE xpole=0.25;
      //      FTYPE eqslope=0.1;
      FTYPE eqslope=0.5;
      FTYPE xeq=0.5;
      V[2] = thetasam2(X[2], poleslope, eqslope, xpole, xeq);
      

    }





    // default is uniform \phi grid
    V[3]=2.0*M_PI*X[3];
  }
}


void dxdxp_analytic_user(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM])
{
  dualfprintf(fail_file,"Should not be computing USERCOORDS analytically\n");
  myexit(34698346);
  dxdxp[3][3] = 2.0*M_PI;    

}
void set_points_user(void)
{


  if(WHICHUSERCOORD==UNIHALFUSERCOORD){
    startx[1] = 0.3999985081775278946780799743777598329673;
    startx[2] = 0.;
    startx[3] = 0.;
    
    FTYPE endx1=21.40529883372801383045167738115556702610;
    dx[1] = (endx1-startx[1]) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];
  }

  if(WHICHUSERCOORD==JET6LIKEUSERCOORD){
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    startx[3] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 1.0/totalsize[3];

#if(1)
    startx[1] = log(Rin-R0)/npow;

    trifprintf( "ITERATIVE dx1: Rout=%21.15g R0=%21.15g npow=%21.15g cpow2=%21.15g cpow3=%21.15g npow2=%21.15g x1br=%21.15g rbr=%21.15g\n",Rout,R0,npow,cpow2,cpow3,npow2,x1br,rbr);

    FTYPE x1max0, x1max,dxmax;
    int iter;
    const FTYPE RELACC = NUMEPSILON*100.0;
    const int ITERMAX = 100;


    if( Rout < rbr ) {
      x1max = log(Rout-R0)/npow;
    }
    else {
      x1max0 = 1.1*x1br;
      x1max = 1.2*x1br;

      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {

        // trifprintf( "iter=%d x1max=%21.15g x2max0=%21.15g\n",iter,x1max0,x1max);

        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;

        if(1){
          dxmax= (pow( (log(Rout-R0) - npow*x1max0)/cpow2, 1./npow2 ) + x1br*1.0) - x1max0;
        }
        else{
          // f-f0 = (x-x0)*dfdx -> if f=Rout -> x = (Rout-f0)/dfdx+x0
          
          FTYPE dVdx1=(npow + cpow2*npow2*cpow3*pow(cpow3*(x1max0-x1br*1.0),npow2-1.0)) * exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          FTYPE V0 = R0 + exp(npow*x1max0 + cpow2 * pow(cpow3*(x1max0-x1br*1.0),npow2));
          
          dxmax=(Rout-V0)/dVdx1; // x-x0

          dualfprintf(fail_file,"dVdx1=%g V0=%g dxmax=%g x1max=%g x1max0=%g\n",dVdx1,V0,dxmax,x1max,x1max0);
        }

        // need a slight damping factor
        FTYPE dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;


      }

      if( iter == ITERMAX ) {
        trifprintf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                    x1max, (x1max-x1max0)/x1max, iter );
        exit(1);
      }
      else {
        trifprintf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }

    dx[1] = ( x1max - startx[1] ) /totalsize[1];
#endif


  }

}
FTYPE setRin_user(int ihor, FTYPE ihoradjust)
{
  if(1){
    FTYPE ftemp;

     // see jet3coords_checknew.nb (and fix_3dpolestissue.nb) to have chosen Rin and ihor and compute required R0
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else if(npow2>0){
      return(1.2);
    }
    else{
      dualfprintf(fail_file,"ihoradjust=%21.15g totalsize[1]=%d Rhor=%21.15g R0=%21.15g npow=%21.15g Rout=%21.15g\n",ihoradjust,totalsize[1],Rhor,R0,npow,Rout);
      return(R0+exp( pow((totalsize[1]*pow(log(Rhor-R0),1.0/npow) - ihoradjust*pow(log(Rout-R0),1.0/npow))/(totalsize[1]-ihoradjust),npow)));
    }
  }
}
 
