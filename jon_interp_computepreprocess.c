#include "decs.h"



// local declarations of local functions
static void vec2vecortho(int outputvartypelocal, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho);
static void vB2poyntingdensity(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *compout);
static void vecup2vecdowncomponent(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vec, FTYPE *compout);
static void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom);
static void generate_lambdacoord(int oldgridtype, int newgridtype, FTYPE *V, FTYPE (*lambdacoord)[NDIM]);
static void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon);







// process inputted data
void compute_preprocess(int outputvartypelocal, FILE *gdumpfile, FTYPE *finaloutput)
{
  FTYPE vec[NDIM],vecv[NDIM],vecB[NDIM];
  FTYPE vecortho[NDIM];
  int jj;
  int dir;
  int iter;
  // fastest index is most-right element
  int ti[NDIM];
  FTYPE X[NDIM];
  FTYPE V[NDIM];
  FTYPE conn[NDIM][NDIM][NDIM];
  FTYPE gcon[SYMMATRIXNDIM];
  FTYPE gcov[SYMMATRIXNDIM];
  FTYPE gdet;
  FTYPE ck[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  struct of_geom geom;
  FTYPE s_fvar;
  FTYPE *fvar=&s_fvar;
  extern int Mcon_calc(FTYPE *pr, struct of_state *q, FTYPE (*Mcon)[NDIM]); // in phys.tools.c
  extern void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);



  // if good memory space in finaloutput, then use it
  if(finaloutput!=NULL) fvar=finaloutput;
  

  // no need to recompute X,V if reading in gdump
  // assumes input data at CENT
  //  coord(i,j,k,CENT,X);
  // use bl_coord()
  //bl_coord(X,V);




  // perform local transformation of tensor objects prior to spatial interpolation
  if(outputvartypelocal<=10){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);



    // convert coordinate basis vector compnents to single orthonormal basis component desired
    fscanf(stdin,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
    
    // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use
    vec2vecortho(outputvartypelocal,ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vec, vecortho);
    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(stdout,vecortho[jj]);
    }
    else{ // else store (can only store 1 of them)
      fvar[0]=vecortho[vectorcomponent];
    }

  }
  else if(outputvartypelocal==11){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);

    // input uu0 vu1 vu2 vu3
    fscanf(stdin,SCANARG4VEC,&vecv[0],&vecv[1],&vecv[2],&vecv[3]) ; // vu^i=uu^i/uu0 (i.e. not uu^i as maybe expected)
    SLOOPA(jj) vecv[jj]*=vecv[0]; // now uu[jj]
    // input B^1 B^2 B^3
    vecB[0]=0.0;
    fscanf(stdin,SCANARGVEC,&vecB[1],&vecB[2],&vecB[3]) ;
    // compute
    vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponent, vecv, vecB, &fvar[0]);
    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(stdout,fvar[0]);
    }
    // already stored in fvar[0]

  }
  else if(outputvartypelocal==12){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);

    // input
    fscanf(stdin,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
    // compute
    vecup2vecdowncomponent(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponent, vec, &fvar[0]);
    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(stdout,fvar[0]);
    }
    // already stored in fvar[0]

  }
  else if(outputvartypelocal==13){

    if(immediateoutput!=1){
      dualfprintf(fail_file,"outputvartype==%d should have immediateoutput==1\n",outputvartypelocal);
      exit(1);
    }

    // assumes input is field line file format of:
    // rho0, u, -u_t, -T^r_t/(rho0 u^r), u^t, v^r, v^\theta, v^\phi, B^r, B^\theta, B^\phi

    // compute Full Diag

    // assume always immediateoutput==1 since many quantities and faster/cleaner to first generate all quantities (per grid point) and then do interpolation, integration, or averaging later since then don't have to process gdump or compute many orthonormal conversions.


    FTYPE pr[NPR],vcon[NDIM],negud0,muconst,uu0;
    FTYPE dFuu[NDIM][NDIM];
    FTYPE Tud[NDIM][NDIM];
    FTYPE TudMA[NDIM][NDIM];
    FTYPE TudEM[NDIM][NDIM];
    FTYPE uortho[NDIM], bortho[NDIM], Bortho[NDIM];
    int concovtype;
    struct of_state q;
    //    struct of_geom geom;
    FTYPE ucov[NDIM],ucon[NDIM];
    FTYPE bcon[NDIM],bcov[NDIM],bsq;
    FTYPE Bcon[NDIM],Bcov[NDIM];
    FTYPE Tit[NDIM];
    FTYPE Titsq,dPdw;
    int kk;
    FTYPE flux[NDIM];
    FTYPE fluxdiagpress[NDIM];



    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);
  

    // scan-in FIELDLINE content
    fscanf(stdin,SCANFIELDLINE,&pr[RHO],&pr[UU],&negud0,&muconst,&uu0,&pr[U1],&pr[U2],&pr[U3],&pr[B1],&pr[B2],&pr[B3]);


    ////////////////////
    //
    //  Assignments
    //
    ////////////////////

    // fix u^i to be lab-frame 4-velocity
    vcon[TT]=1.0;
    vcon[RR]=pr[U1];
    vcon[TH]=pr[U2];
    vcon[PH]=pr[U3];

    ucon[TT]=vcon[TT]*uu0;
    ucon[RR]=vcon[RR]*uu0;
    ucon[TH]=vcon[TH]*uu0;
    ucon[PH]=vcon[PH]*uu0;

    // make primitive u^i
    pr[U1]=ucon[RR];
    pr[U2]=ucon[TH];
    pr[U3]=ucon[PH];

    Bcon[TT]=0;
    Bcon[RR]=pr[B1];
    Bcon[TH]=pr[B2];
    Bcon[PH]=pr[B3];

    ////////////////////
    //
    //  PRECOMPUTATIONS
    //
    ////////////////////

    // compute lower components
    DLOOPA(jj) ucov[jj]=0.0;
    DLOOP(jj,kk) ucov[jj] += ucon[kk]*gcov[GIND(jj,kk)];

    // compute b^\mu[pr,ucon,ucov]
    bcon_calc(pr, ucon, ucov, bcon);

    // compute b_\mu
    DLOOPA(jj) bcov[jj]=0.0;
    DLOOP(jj,kk) bcov[jj] += bcon[kk]*gcov[GIND(jj,kk)];

    // compute "fake" B_\mu (i.e. didn't lower upper t in *F^{t\mu} )
    DLOOPA(jj) Bcov[jj]=0.0;
    DLOOP(jj,kk) Bcov[jj] += Bcon[kk]*gcov[GIND(jj,kk)];

    // compute b^2
    bsq = dot(bcon,bcov);


    ////////////////////
    //
    //  STRUCTURE Assignments (must be consistent with global.structs.h or at least how structures used in phys.tools.c as some functions pulled from there are used below)
    //
    ////////////////////

    // assign of_state structure
    DLOOPA(jj){
      q.ucon[jj]=ucon[jj];
      q.ucov[jj]=ucov[jj];
      q.bcon[jj]=bcon[jj];
      q.bcov[jj]=bcov[jj];
    }
    q.pressure=(gam-1)*pr[UU]; // assumes ideal gas
    q.bsq=bsq;
    q.entropy=0; // ignore
    q.ifremoverestplus1ud0elseud0=1.0+ucov[TT];
#if(MERGEDC2EA2CMETHOD)
    // for merged method and stored by compute_and_store_???() functions
    q.gdet=gdet;
#if(WHICHEOM!=WITHGDET)
    for(jj=0;jj<NPR;jj++) eomfunc[jj]=gdet;
#endif
    for(jj=0;jj<NPR;jj++) q.prim[jj]=pr[jj];
    DLOOPA(jj) q.Blower[jj]=Bcov[jj];
    DLOOPA(jj) q.vcon[jj]=vcon[jj];
    DLOOPA(jj) q.gdetBcon[jj]=gdet*Bcon[jj];
    q.overut=1.0/ucon[TT];
#endif



    ////////////////////
    //
    //  COMPUTATIONS
    //
    ////////////////////



    // T^p_t[EM] = b^2 u^p u_t  - b^p b_t
    // EM poloidal flux only
    DLOOPA(jj){
      if(jj==1 || jj==2){
	Tit[jj]=bsq*ucon[jj]*ucov[0] - bcon[jj]*bcov[0] + delta(jj,0)*(bsq*0.5);
      }
      else Tit[jj]=0.0;
    }
    
    // |T^p_t|
    Titsq = 0.0;
    DLOOP(jj,kk) Titsq+=Tit[jj]*Tit[kk]*gcov[GIND(jj,kk)];
    
    // dP/d\omega = \detg T^p_t[EM]/sin(\theta)
    dPdw=sqrt(fabs(Titsq))*fabs(gdet/sin(V[2]));


    // TODO  GODMARK: in phys.tools.c, may need to split it or something since too many dependencies
    // dF^{\mu\nu}
    //    Mcon_calc(pr,&q,dFuu); // used to get *true Bortho* and Omega_F, etc.
    //    mhd_calc_ma(pr, TT, &geom, &q, &flux[UU], &fluxdiagpress[UU]); // fills flux[UU->U3] and fluxdiagonal[UU->U3]



    // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use
    concovtype=1; // means inputting u^\mu
    vec2vecortho(concovtype,ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, ucon, uortho);
    vec2vecortho(concovtype,ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, bcon, bortho);
    vec2vecortho(concovtype,ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, Bcon, Bortho);


    ////////////////////
    //
    //  CONSTRUCT ACTUAL THINGS TO WRITE TO FILE
    //
    ////////////////////

    FTYPE dMdot,dEdot,dLdot,TrtMA,TrtEM,TrphiMAa,TrphiEM;
    FTYPE dPsir, dPsitheta;
    FTYPE Thatrphi,pg,pb,dAp,dAL;
    FTYPE dhor;
    FTYPE omegaf;
    FTYPE dLum;
    FTYPE Deltahattheta,Omega,vahat;
    FTYPE rho0; // and -u_t u^t
    FTYPE gdetB1,gdetB2; // to compute A_\phi // and B^ihat v^ihat
    FTYPE Tphinu[NDIM]; // force-related term that eventually is differenced in x3-direction
    FTYPE term1[NDIM],term2[NDIM],term3[NDIM],term4[NDIM]; // connection-related terms
    
      
    dMdot=gdet*pr[RHO]*pr[U1]*dX[2]*dX[3];




    ////////////////////
    //
    //  OUTPUT
    //
    ////////////////////


    DLOOPA(jj) writeelement(stdout,vecortho[jj]);



  }






  if(immediateoutput==1){
    // output return after entire row is done
    fprintf(stdout,"\n") ;
  }


}



/////////////////////////////////////
//
// coordinate TRANSFORMATION FUNCTIONS
//
/////////////////////////////////////

// ./iinterp 1 3 1 1 128 128 32 1 0 0 1 0 128 128 128 -40 40 -40 40 -40 40 1.1 40 0 0.3 9 0 4 < u0000 > iu0000
// Profile:
// Each sample counts as 0.01 seconds.
//  %   cumulative   self              self     total
// time   seconds   seconds    calls  us/call  us/call  name
// 30.24     59.50    59.50  1099104    54.14    97.52  bicubic_interp_wrap
// 23.99    106.71    47.21   138724   340.33   340.33  dervs_for_bicubic
// 15.94    138.08    31.36                             pow.L
//  6.52    150.90    12.82 13948003     0.92     1.75  dxdxp_numerical
//  5.97    162.65    11.75                             atan.L
//  4.18    170.88     8.23 258855333     0.03     0.04  bl_coord
//  4.14    179.02     8.14                             exp.L
//  2.49    183.92     4.90 482023381     0.01     0.01  mysin
//  2.45    188.75     4.83                             sin.L
//  0.55    189.84     1.09 13948003     0.08     0.08  ludcmp


// coordinate transform of vector and get single component of result
//static void vec2vecortho(FILE *gdumpfile, int oldgridtype, int newgridtype, int i, int j, int k, FTYPE *vec, FTYPE *vecortho)
// concovtype: 1 = inputting u^\mu   2 = inputting u_\mu
static void vec2vecortho(int concovtype, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho)
{
  FTYPE lambdacoord[NDIM][NDIM];
  int jj,kk;
  FTYPE tetrcov[NDIM][NDIM],tetrcon[NDIM][NDIM],eigenvalues[NDIM];
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];


  // vector here is in original X coordinates
  DLOOPA(jj) finalvec[jj]=vec[jj];


  // get SPC -> Cart transformation
  generate_lambdacoord(oldgridtype, newgridtype, V, lambdacoord);

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d lambdacoord=%21.15g\n",jj,kk,lambdacoord[jj][kk]);
  //  }


  // get tetrad (uses dxdxp so that tetrcon and tetrcon and eigenvalues are using V metric not X metric
  tetr_func_frommetric(dxdxp, gcov, tetrcov, tetrcon, eigenvalues);

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d gcov=%21.15g\n",jj,kk,gcov[GIND(jj,kk)]);
  //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d tetrcon=%21.15g\n",jj,kk,tetrcon[jj][kk]);
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d eigenvalues=%21.15g\n",jj,eigenvalues[jj]);
  //  }


  // transform from X to V for contravariant vector
  DLOOPA(jj) tempcomp[jj]=0.0;
  DLOOP(jj,kk){
    tempcomp[jj] += dxdxp[jj][kk]*finalvec[kk];
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];


  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d finalvec(dxdxp)=%21.15g\n",jj,finalvec[jj]);
  //  }


  DLOOPA(jj) tempcomp[jj]=0.0;

  if(concovtype==1){
    // transform to orthonormal basis for contravariant vector in V coordinates
    DLOOP(jj,kk){
      //    tempcomp[kk] += tetrcon[kk][jj]*finalvec[jj];
      tempcomp[kk] += tetrcov[kk][jj]*finalvec[jj];
    }
  }
  else if(concovtype==2){
    // transform to orthonormal basis for covariant vector in V coordinates (GODMARK: unsure about tetrcon[kk][jj] vs. tetrcon[jj][kk])
    DLOOP(jj,kk){
      tempcomp[kk] += tetrcon[kk][jj]*finalvec[jj];
    }
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];



  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d finalvec(tetrcov)=%21.15g\n",jj,finalvec[jj]);
  //  }

  // transform from spherical polar to Cartesian coordinates (assumes vector is in orthonormal basis)
  DLOOPA(jj) tempcomp[jj]=0.0;
  DLOOP(jj,kk){
    tempcomp[kk] += lambdacoord[kk][jj]*finalvec[jj];
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d finalvec(lambdacoord)=%21.15g\n",jj,finalvec[jj]);
  //  }
  

  // final answer:
  DLOOPA(jj) vecortho[jj]=finalvec[jj];


}




static void vB2poyntingdensity(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *compout)
{
  FTYPE newgcov[SYMMATRIXNDIM];
  int jj,kk,ll,pp;
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];
  FTYPE pr[NPR];
  FTYPE ucov[NDIM],ucon[NDIM];
  FTYPE bcon[NDIM],bcov[NDIM],bsq;
  FTYPE Tit[NDIM];
  FTYPE Titsq,dPdw;



  pr[B1]=vecB[1];
  pr[B2]=vecB[2];
  pr[B3]=vecB[3];
  
  // compute lower components
  DLOOPA(jj) ucon[jj]=vecv[jj];
  DLOOPA(jj) ucov[jj]=0.0;
  DLOOP(jj,kk) ucov[jj] += ucon[kk]*gcov[GIND(jj,kk)];

  // compute b^\mu[pr,ucon,ucov]
  bcon_calc(pr, ucon, ucov, bcon);

  // compute b_\mu
  DLOOPA(jj) bcov[jj]=0.0;
  DLOOP(jj,kk) bcov[jj] += bcon[kk]*gcov[GIND(jj,kk)];

  // compute b^2
  bsq = dot(bcon,bcov);

  // T^p_t[EM] = b^2 u^p u_t  - b^p b_t
  // EM poloidal flux only
  DLOOPA(jj){
    if(jj==1 || jj==2){
      Tit[jj]=bsq*ucon[jj]*ucov[0] - bcon[jj]*bcov[0] + delta(jj,0)*(bsq*0.5);
    }
    else Tit[jj]=0.0;
  }

  // |T^p_t|
  Titsq = 0.0;
  DLOOP(jj,kk) Titsq+=Tit[jj]*Tit[kk]*gcov[GIND(jj,kk)];

  // dP/d\omega = \detg T^p_t[EM]/sin(\theta)
  dPdw=sqrt(fabs(Titsq))*fabs(gdet/sin(V[2]));

  // final output
  *compout = dPdw;


}




// input PRIMECOORD coordinate basis v^i and output MCOORD coordinate basis v_{vectorcomponent}
// so no transformation to orthonormal basis or even to Cartesian coordaintes
static void vecup2vecdowncomponent(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vec, FTYPE *compout)
{
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE newgcov[SYMMATRIXNDIM];
  int jj,kk,ll,pp;
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];


  // get inverse coordinate transformation
  idxdxprim(dxdxp, idxdxp); // from coord.c

  // get MCOORD metric from PRIMECOORD metric
  DLOOP(jj,kk){
    newgcov[GIND(jj,kk)]=0.0;
    DLOOP(ll,pp) {
      newgcov[GIND(jj,kk)] += GINDASSIGNFACTOR(jj,kk)*gcov[GIND(ll,pp)]*idxdxp[ll][jj]*idxdxp[pp][kk];
    }
  }


  DLOOPA(jj) finalvec[jj]=vec[jj];

  // compute lower components
  DLOOPA(jj) tempcomp[jj]=0.0;
  DLOOP(jj,kk) tempcomp[jj] += finalvec[kk]*newgcov[GIND(jj,kk)];
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];

  // finally get desired component
  *compout = finalvec[vectorcomponent];


}




// generates transformation matrix from old to new grid type for simple coordinates
static void generate_lambdacoord(int oldgridtype, int newgridtype, FTYPE *V, FTYPE (*lambdacoord)[NDIM])
{
  FTYPE r,th,ph;
  int jj,kk;
  FTYPE lambdatrans[NDIM][NDIM];
  FTYPE lambdarotate[NDIM][NDIM];
  FTYPE lambdatemp[NDIM][NDIM];



  if(oldgridtype==newgridtype || newgridtype==GRIDTYPENOCHANGE){
    // no change
    DLOOP(jj,kk) lambdatrans[jj][kk]=0.0;
    DLOOPA(jj) lambdatrans[jj][jj]=1.0;
  }
  else if(oldgridtype==GRIDTYPESPC && newgridtype==GRIDTYPECART){
    // y=0 is \phi=0
    // x = r sin(\theta) cos(\phi)
    // y = r sin(\theta) sin(\phi)
    // z = r cos(\theta)

    r=V[1];
    th=V[2];
    ph=V[3];

    // Will be assuming only convering contravariant quantities
    // 
    // vout^\mu = \Lambda^\mu_\nu u^\nu
    //
    // where since u^\nu = dxold^\nu/dq then \Lambda^\mu_\nu = dxnew^\mu/dxold^\nu
    //
    // so for going from SPC to Cart with xnew=Cart and xold=SPC, then (e.g.) need terms like \Lambda^x_r = dx/dr

    // assume time doesn't change or mix with space
    lambdatrans[TT][TT]=1.0;
    SLOOPA(jj) lambdatrans[TT][jj] = lambdatrans[jj][TT] = 0.0;

    // rest come from definitions of {x,y,z}(r,\theta,\phi)
    // assumes orthonormal to orhonormal!
    lambdatrans[1][RR] = sin(th)*cos(ph);
    lambdatrans[1][TH] = cos(th)*cos(ph);
    lambdatrans[1][PH] = -sin(ph);

    lambdatrans[2][RR] = sin(th)*sin(ph);
    lambdatrans[2][TH] = cos(th)*sin(ph);
    lambdatrans[2][PH] = cos(ph);

    lambdatrans[3][RR] = cos(th);
    lambdatrans[3][TH] = -sin(th);
    lambdatrans[3][PH] = 0.0;
  }
  else if(oldgridtype==GRIDTYPESPC && newgridtype==GRIDTYPECARTLIGHT){
    // similar to above, but time mixes with space

    r=V[1];
    th=V[2];
    ph=V[3];

    lambdatrans[TT][TT]=1.0;
    SLOOPA(jj) lambdatrans[TT][jj] = 0.0;
    lambdatrans[TT][RR] = -cos(tnrdegrees*M_PI/180.0); // dtobs/dr
    SLOOPA(jj) lambdatrans[jj][TT] = 0.0;

    // rest come from definitions of {x,y,z}(r,\theta,\phi)
    // assumes orthonormal to orhonormal!
    lambdatrans[1][RR] = sin(th)*cos(ph);   // dxhat/drhat
    lambdatrans[1][TH] = cos(th)*cos(ph);   // dxhat/dhhat
    lambdatrans[1][PH] = -sin(ph);   // dxhat/dphhat

    lambdatrans[2][RR] = sin(th)*sin(ph); // dyhat/drhat
    lambdatrans[2][TH] = cos(th)*sin(ph); // dyhat/dhhat
    lambdatrans[2][PH] = cos(ph); // dyhat/dphhat

    lambdatrans[3][RR] = cos(th); // dzhat/drhat
    lambdatrans[3][TH] = -sin(th); // dzhat/dhhat
    lambdatrans[3][PH] = 0.0; // dzhat/dphhat
  }
  else{
    dualfprintf(fail_file,"No transformation setup for oldgridtype=%d newgridtype=%d\n",oldgridtype,newgridtype);
    myexit(246347);
  }


  // initialize rotation matrix
  DLOOP(jj,kk) lambdarotate[jj][kk]=0.0;
  DLOOPA(jj) lambdarotate[jj][jj]=1.0;

  // assume have converted from spc to Cart, now go from Cart to new Cart rotated around y-axis
  // rotate around y-axis, so y-axis is unchanged
  // here 1=x 2=y 3=z (unlike in interpolation code that is x z y)
  // e.g. tnrdegrees=90deg : +znonrot -> +xrot  & +xnonrot -> -zrot  so rotation is from z-axis towards x-axis.  Or with y-axis pointed at you, rotation is counter-clockwise.
  lambdarotate[1][1] = cos(tnrdegrees*M_PI/180.0);
  lambdarotate[1][2] = 0.0;
  lambdarotate[1][3] = sin(tnrdegrees*M_PI/180.0); // e.g. xCartrot = zCartnonrot for tnrdegrees=90deg

  lambdarotate[2][1] = 0.0;
  lambdarotate[2][2] = 1.0;
  lambdarotate[2][3] = 0.0;

  lambdarotate[3][1] = -sin(tnrdegrees*M_PI/180.0); // e.g. zCartrot = -xCartnonrot for tnrdegrees=90deg
  lambdarotate[3][2] = 0.0;
  lambdarotate[3][3] = cos(tnrdegrees*M_PI/180.0);

  // first transform and then rotate
  // Lambda^\mu_\nu = \Lambdarotate^\mu_\kappa \Lambdatrans^\kappa_\nu
  int ll;
  DLOOP(jj,kk){
    lambdacoord[jj][kk]=0.0;
    DLOOPA(ll){
      lambdacoord[jj][kk] += lambdarotate[jj][ll]*lambdatrans[ll][kk];
    }
  }




}





// read a line from the gdump file
// assumed to be used in same order as data and gdump file as looped over in this code
// and of course only can be done once per data point
// fastest index is most-right element
static void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom)
{
  int jj,kk,ll,pp;

  // see HARM's dump.c and SM's gammie.m grid3d
  SLOOPA(jj) fscanf(in,"%d",&ti[jj]);
  SLOOPA(jj) fscanf(in,SCANARG,&X[jj]);
  SLOOPA(jj) fscanf(in,SCANARG,&V[jj]);
  DLOOPA(jj) DLOOPA(kk) DLOOPA(ll) fscanf(in,SCANARG,&conn[jj][kk][ll]);
  DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&gcon[GIND(jj,kk)]);
  DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&gcov[GIND(jj,kk)]);
  fscanf(in,SCANARG,gdet); // gdet already pointer
  DLOOPA(jj) fscanf(in,SCANARG,&ck[jj]);
  DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&dxdxp[jj][kk]);


  // for debug:
  tiglobal[0]=ti[0];
  tiglobal[1]=ti[1];
  tiglobal[2]=ti[2];
  tiglobal[3]=ti[3];



#if(0) // INPROGRESS


  // assign of_geom structure

  // dummy space for gset() version
  FTYPE gengcov[SYMMATRIXNDIM];
  FTYPE gengcovpert[NDIM];
  FTYPE alphalapse;
  FTYPE betasqoalphasq;
  FTYPE beta[NDIM];
  FTYPE gengcon[SYMMATRIXNDIM];

  // bit faster since not all values always used
  FTYPE *gcov;
  FTYPE *gcon;
  FTYPE *gcovpert;

  FTYPE gdet,igdetnosing;
#if(GDETVOLDIFF)
  FTYPE gdetvol;
#endif
#if(WHICHEOM!=WITHGDET)
  FTYPE eomfunc[NPR],ieomfuncnosing[NPR];
#endif
  int i,j,k,p;


#endif



}






/* calculate magnetic field four-vector */
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j <= 3; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}
