#include "decs.h"



// local declarations of local functions

static void vec2vecortho(int concovtype, FTYPE V[],  FTYPE *gcov,  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho);

static void T2Tortho(int concovtype, FTYPE V[],  FTYPE *gcov,  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*T)[NDIM], FTYPE (*Tortho)[NDIM]);


static void vB2poyntingdensity(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *compout);


static void get_beta_gamma_boost(FTYPE *gamma, FTYPE *betavec);

static void getvboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vecvboost);

static void vboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *vecvboosted);

static void Tboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*Tmunuorthocart)[NDIM], FTYPE (*Tmunuboostedcart)[NDIM], FTYPE (*Tmunuboostedspccart)[NDIM], FTYPE (*Tmunuboostedspcspc)[NDIM]);


static void vB2poyntingdensityboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*Tmunuboostedcart)[NDIM], FTYPE (*Tmunuboostedspccart)[NDIM], FTYPE *FEMrad3, FTYPE *FEMrad4);

static void vB2poyntingdensityfull(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vecv, FTYPE *vecB, FTYPE (*Tmunu)[NDIM]);


static int Mcon_calc_uB(FTYPE *ucon, FTYPE *Bcon, FTYPE (*Mcon)[NDIM]);
static int Mcon_calc_ub(FTYPE *ucon, FTYPE *bcon, FTYPE (*Mcon)[NDIM]);
static FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
static void MtoF_simple(int which, FTYPE (*invar)[NDIM], FTYPE gdet, FTYPE (*outvar)[NDIM]);


static void vecup2vecdown(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecout);
static void generate_lambdacoord(int oldgridtype, int newgridtype, FTYPE *V, FTYPE (*lambdacoord)[NDIM]);
static void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon);



static int compute_datatype14(int which, FTYPE *val, FTYPE *fvar, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom);

static int compute_datatype15(int outputvartypelocal, FTYPE *val, FTYPE *fvar, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom);



static FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
static void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE (*faraday)[NDIM]);


// process inputted data
void compute_preprocess(int outputvartypelocal, FILE *gdumpfile, int h, int i, int j, int k, FTYPE*****olddatalocal, FTYPE *finaloutput)
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
  extern void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);



  // if good memory space in finaloutput, then use it
  if(finaloutput!=NULL) fvar=finaloutput;
  

  // no need to recompute X,V if reading in gdump
  // assumes input data at CENT
  // coord(i,j,k,CENT,X);
  // use bl_coord()
  // bl_coord(X,V);




  // perform local transformation of tensor objects prior to spatial interpolation
  if(outputvartypelocal<=10){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);

    // convert coordinate basis vector compnents to single orthonormal basis component desired

    // scan-in appropriate file type
    if(DATATYPE<=10){
      //      fscanf(infile,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
      readelement(binaryinput,inFTYPE,infile,&vec[0]);
      readelement(binaryinput,inFTYPE,infile,&vec[1]);
      readelement(binaryinput,inFTYPE,infile,&vec[2]);
      readelement(binaryinput,inFTYPE,infile,&vec[3]);
    }
    else if(DATATYPE>=1002 && DATATYPE<=1009){
      // fieldline file type input
      FTYPE val[MAXINCOLS];
      int colini;
      for(colini=0;colini<numcolumns;colini++) readelement(binaryinput,inFTYPE,infile,&val[colini]); // fscanf(infile,SCANARG,&val[colini]);
      
      if(DATATYPE>=1002 && DATATYPE<=1005){
        // now assign to vector
        vec[0]=val[4];
        vec[1]=val[4]*val[5];
        vec[2]=val[4]*val[6];
        vec[3]=val[4]*val[7];
      }
      else if(DATATYPE>=1006 && DATATYPE<=1009){
        // now assign to vect
        vec[0]=0.0;
        vec[1]=val[8];
        vec[2]=val[9];
        vec[3]=val[10];
      }
      else dualfprintf(fail_file,"BAD1 DATATYPE=%d\n",DATATYPE);
    }
    else dualfprintf(fail_file,"BAD2 DATATYPE=%d\n",DATATYPE);
      
    
    // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use
    vec2vecortho(outputvartypelocal,V,gcov,dxdxp,oldgridtype, newgridtype, vec, vecortho);

    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(binaryoutput,outFTYPE,outfile,vecortho[jj]);
    }
    else{ // else store (can only store 1 of them)
      fvar[0]=vecortho[vectorcomponent];
    }

  }
  else if(outputvartypelocal==11){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);

    if(DATATYPE==11){
      // input uu0 vu1 vu2 vu3
      //fscanf(infile,SCANARG4VEC,&vecv[0],&vecv[1],&vecv[2],&vecv[3]) ; // vu^i=uu^i/uu0 (i.e. not uu^i as maybe expected)
      readelement(binaryinput,inFTYPE,infile,&vecv[0]);
      readelement(binaryinput,inFTYPE,infile,&vecv[1]);
      readelement(binaryinput,inFTYPE,infile,&vecv[2]);
      readelement(binaryinput,inFTYPE,infile,&vecv[3]);
      SLOOPA(jj) vecv[jj]*=vecv[0]; // now uu[jj]
      // input B^1 B^2 B^3
      vecB[0]=0.0;
      //      fscanf(infile,SCANARGVEC,&vecB[1],&vecB[2],&vecB[3]) ;
      readelement(binaryinput,inFTYPE,infile,&vecB[1]);
      readelement(binaryinput,inFTYPE,infile,&vecB[2]);
      readelement(binaryinput,inFTYPE,infile,&vecB[3]);
    }
    else if(DATATYPE==1011){
      // fieldline file type input
      FTYPE val[MAXINCOLS];
      int colini;
      for(colini=0;colini<numcolumns;colini++) readelement(binaryinput,inFTYPE,infile,&val[colini]); // fscanf(infile,SCANARG,&val[colini]);
      // now assign to vector
      vecv[0]=val[4];
      vecv[1]=val[4]*val[5];
      vecv[2]=val[4]*val[6];
      vecv[3]=val[4]*val[7];
      vecB[0]=0.0;
      vecB[1]=val[8];
      vecB[2]=val[9];
      vecB[3]=val[10];
    }
    else dualfprintf(fail_file,"BAD3 DATATYPE=%d\n",DATATYPE);

    // compute
    vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponent, vecv, vecB, &fvar[0]);
    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(binaryoutput,outFTYPE,outfile,fvar[0]);
    }
    // already stored in fvar[0]

  }
  else if(outputvartypelocal==12){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);

    if(DATATYPE==12){
      // input
      //      fscanf(infile,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
      readelement(binaryinput,inFTYPE,infile,&vec[0]);
      readelement(binaryinput,inFTYPE,infile,&vec[1]);
      readelement(binaryinput,inFTYPE,infile,&vec[2]);
      readelement(binaryinput,inFTYPE,infile,&vec[3]);
    }
    else if(DATATYPE==1012){
      // fieldline file type input
      FTYPE val[MAXINCOLS];
      int colini;
      //      for(colini=0;colini<numcolumns;colini++) fscanf(infile,SCANARG,&val[colini]);
      for(colini=0;colini<numcolumns;colini++) readelement(binaryinput,inFTYPE,infile,&val[colini]); // fscanf(infile,SCANARG,&val[colini]);
      // now assign to vector
      vec[0]=0.0;
      vec[1]=val[8];
      vec[2]=val[9];
      vec[3]=val[10];
    }
    else dualfprintf(fail_file,"BAD4 DATATYPE=%d\n",DATATYPE);
    
    // compute
    FTYPE vecdowntemp[NDIM];
    vecup2vecdown(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vec, vecdowntemp);
    fvar[0]=vecdowntemp[vectorcomponent];
    
    if(immediateoutput==1){ // then immediately write to output
      DLOOPA(jj) writeelement(binaryoutput,outFTYPE,outfile,fvar[0]);
    }
    // already stored in fvar[0]

  }
  else if(outputvartypelocal==13){ // full diag output (GODMARK: NOT DONE)

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
    //    fscanf(infile,SCANFIELDLINE,&pr[RHO],&pr[UU],&negud0,&muconst,&uu0,&pr[U1],&pr[U2],&pr[U3],&pr[B1],&pr[B2],&pr[B3]);
    readelement(binaryinput,inFTYPE,infile,&pr[RHO]);
    readelement(binaryinput,inFTYPE,infile,&pr[UU]);
    readelement(binaryinput,inFTYPE,infile,&negud0);
    readelement(binaryinput,inFTYPE,infile,&muconst);
    readelement(binaryinput,inFTYPE,infile,&uu0);
    readelement(binaryinput,inFTYPE,infile,&pr[U1]);
    readelement(binaryinput,inFTYPE,infile,&pr[U2]);
    readelement(binaryinput,inFTYPE,infile,&pr[U3]);
    readelement(binaryinput,inFTYPE,infile,&pr[B1]);
    readelement(binaryinput,inFTYPE,infile,&pr[B2]);
    readelement(binaryinput,inFTYPE,infile,&pr[B3]);


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
    //    mhd_calc_ma(pr, TT, &geom, &q, &flux[UU], NULL, &fluxdiagpress[UU], NULL); // fills flux[UU->U3] and fluxdiagonal[UU->U3]



    // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use
    concovtype=1; // means inputting u^\mu
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, ucon, uortho);
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, bcon, bortho);
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, Bcon, Bortho);


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


    DLOOPA(jj) writeelement(binaryoutput,outFTYPE,outfile,vecortho[jj]);



  }
  else if(outputvartypelocal==14){ // reading in field line data and outputting "numoutputcols" columns of results for interpolation


    if(docurrent==0){ // then not computing current, so do normal procedure

      if(numoutputcols!=10){
        fprintf(stderr,"numoutputcols=%d != %d as expected for datatype==14 with docurrent=%d\n",numoutputcols,10,docurrent);
        myexit(1);
      }


      // first get gdump data (only once per call to compute_preprocess() !!)
      read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);


      // normal read-in
      int colini;
      FTYPE val[MAXINCOLS];
      for(colini=0;colini<numcolumns;colini++) readelement(binaryinput,inFTYPE,infile,&val[colini]); // numcolumns is input number of columns

      // 1 means compute all things except current
      compute_datatype14(1, val, fvar, ti,  X,  V,  conn,  gcon,  gcov,  gdet,  ck,  dxdxp, &geom);

    }
    else{

      if(numoutputcols!=14){
        fprintf(stderr,"numoutputcols=%d != %d as expected for datatype==14 with docurrent=%d\n",numoutputcols,14,docurrent);
        myexit(1);
      }


      // compute_additionls() already computed per-point stuff and put it in olddata0[]
      // just do assignment from olddata0 to fvar (this just flips back to oldata0 after out of this function)
      int coli; for(coli=0;coli<numoutputcols;coli++) fvar[coli]=olddata0[coli][h][i][j][k];

    }


    //DEBUG
    //for(coli=0;coli<numoutputcols;coli++) dualfprintf(fail_file,"coli=%d fvar=%g\n",coli,fvar[coli]);


    if(immediateoutput==1){ // then immediately write to output
      int coli; for(coli=0;coli<numoutputcols;coli++) writeelement(binaryoutput,outFTYPE,outfile,fvar[coli]);
    }
    else{
      // already stored in fvar
    }
  }
  else if(outputvartypelocal==15 || outputvartypelocal==16  || outputvartypelocal==17  || outputvartypelocal==18 || outputvartypelocal==19){ // reading in field line data and outputting "numoutputcols" columns of results for interpolation

    if(docurrent==1){
      dualfprintf(fail_file,"Can't do outputvartypelocal==15 with docurrent=1\n");
      myexit(5235);
    }

    //    if(numoutputcols!=4 && numoutputcols!=6){
    //      fprintf(stderr,"numoutputcols=%d != %d (%d) as expected for datatype==15 or 16\n",numoutputcols,4,6);
    //      myexit(1);
    //    }


    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);
    
    
    // normal read-in
    int colini;
    FTYPE val[MAXINCOLS];
    for(colini=0;colini<numcolumns;colini++) readelement(binaryinput,inFTYPE,infile,&val[colini]); // numcolumns is input number of columns
    
    compute_datatype15(outputvartypelocal, val, fvar, ti,  X,  V,  conn,  gcon,  gcov,  gdet,  ck,  dxdxp, &geom);
    

    //DEBUG
    //for(coli=0;coli<numoutputcols;coli++) dualfprintf(fail_file,"coli=%d fvar=%g\n",coli,fvar[coli]);


    if(immediateoutput==1){ // then immediately write to output
      int coli; for(coli=0;coli<numoutputcols;coli++) writeelement(binaryoutput,outFTYPE,outfile,fvar[coli]);
    }
    else{
      // already stored in fvar
    }



  }
  else{
    fprintf(stderr,"No such outputvartypelocal=%d\n",outputvartypelocal);
    exit(1);
  }





  if(immediateoutput==1){
    // output return after entire row is done
    if(binaryoutput==0) fprintf(outfile,"\n") ;
  }


}














// compute datatype==14 stuff
// val is input stuff from fieldline file (colini - numcolumns)
// fvar is final stuff (coli - numoutputcols)
static int compute_datatype14(int which, FTYPE *val, FTYPE *fvar, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom)
{
  int concovtype;
  int jj;




  /////////////////////
  // also assign to coli type with numoutputcols data values
  int coli;




  if(which==0 || which==2){
    // recover this position's 4-current (J^\mu) from 3D storage in "olddatacurrent" (same as olddatalocal passed to this function, but just access global)
    FTYPE Jcon[NDIM];
    FTYPE Jortho[NDIM];
    if(olddatacurrent!=NULL){
      int href=0;
      DLOOPA(jj) Jcon[jj]=olddatacurrent[jj][href][ti[RR]][ti[TH]][ti[PH]];
      
      // get orthonormal J^\mu
      concovtype=1; // contravariant
      vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, Jcon, Jortho);
    }


    if(olddatacurrent!=NULL){
      fvar[OUTJX0]=Jortho[TT];
      fvar[OUTJX1]=Jortho[RR];
      fvar[OUTJX2]=Jortho[TH];
      fvar[OUTJX3]=Jortho[PH];
    }
    else if(numoutputcols==14){ // assume this is a testing/debug mode or disabled current calculation
      fvar[OUTJX0]=0.0;
      fvar[OUTJX1]=0.0;
      fvar[OUTJX2]=0.0;
      fvar[OUTJX3]=0.0;
    }


  }



  if(which==0 || which==1){
    // now do normal assign to vector (for middle time if doing3time==1)
    FTYPE vecv[NDIM],vecB[NDIM];
    vecv[0]=val[FLU0]; vecv[1]=val[FLV1]; vecv[2]=val[FLV2]; vecv[3]=val[FLV3];
    SLOOPA(jj) vecv[jj]*=vecv[TT]; // now uu[jj]
    vecB[TT]=0.0; SLOOPA(jj) vecB[jj]=val[FLB1+jj-1];

  


    // DEBUG:
    //    SLOOPA(jj) dualfprintf(fail_file,"jj=%d vecv=%g vecB=%g\n",jj,vecv[jj],vecB[jj]);


    // convert coordinate basis vector compnents to single orthonormal basis component desired

    // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use

    // do vecv
    FTYPE vecvortho[NDIM];
    concovtype=1; // contravariant
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vecv, vecvortho);

    // do vecB
    FTYPE vecBortho[NDIM];
    concovtype=1; // contravariant
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vecB, vecBortho);

    int vectorcomponentlocal;

    // compute FEMrad
    FTYPE FEMrad;
    vectorcomponentlocal=1; // radial component
    vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponentlocal, vecv, vecB, &FEMrad);

    // compute B_{\phi}
    vectorcomponentlocal=3; // \phi component
    FTYPE vecBdown[NDIM],Bphi;
    vecup2vecdown(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vecB, vecBdown);
    Bphi=vecBdown[vectorcomponentlocal];


    fvar[OUTRHO]=val[FLRHO]; // rho0
    fvar[OUTU]=val[FLU]; // ug
    fvar[OUTV1]=vecvortho[RR]; // vx
    fvar[OUTV2]=vecvortho[TH]; // vy
    fvar[OUTV3]=vecvortho[PH]; // vz
    fvar[OUTB1]=vecBortho[RR]; // Bx
    fvar[OUTB2]=vecBortho[TH]; // By
    fvar[OUTB3]=vecBortho[PH]; // Bz
    fvar[OUTFEMRAD]=FEMrad; // FEMrad = radial energy flux per unit sin(\theta)
    fvar[OUTBPHI]=Bphi; // B_{\phi} == poloidal current


  }


  return(0);
}

// compute datatype==15 stuff
// val is input stuff from fieldline file (colini - numcolumns)
// fvar is final stuff (coli - numoutputcols)
static int compute_datatype15(int outputvartypelocal, FTYPE *val, FTYPE *fvar, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom)
{
  int concovtype;
  int jj;




  /////////////////////
  // also assign to coli type with numoutputcols data values
  int coli;




  // now do normal assign to vector (for middle time if doing3time==1)
  FTYPE vecv[NDIM],vecB[NDIM];
  vecv[0]=val[FLU0]; vecv[1]=val[FLV1]; vecv[2]=val[FLV2]; vecv[3]=val[FLV3];
  SLOOPA(jj) vecv[jj]*=vecv[TT]; // now uu[jj]
  vecB[TT]=0.0; SLOOPA(jj) vecB[jj]=val[FLB1+jj-1];

  FTYPE vecvrad[NDIM];
  if(outputvartypelocal==19){
    vecvrad[0]=val[FLURAD0]; vecvrad[1]=val[FLVRAD1]; vecvrad[2]=val[FLVRAD2]; vecvrad[3]=val[FLVRAD3];
    SLOOPA(jj) vecvrad[jj]*=vecvrad[TT]; // now uu[jj]
  }


  // DEBUG:
  //    SLOOPA(jj) dualfprintf(fail_file,"jj=%d vecv=%g vecB=%g\n",jj,vecv[jj],vecB[jj]);


  // convert coordinate basis vector compnents to single orthonormal basis component desired

  // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use

  // do vecv
  FTYPE vecvortho[NDIM];
  concovtype=1; // contravariant
  vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vecv, vecvortho);

  FTYPE vecvradortho[NDIM];
  if(outputvartypelocal==19){
    // do vecvrad
    concovtype=1; // contravariant
    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vecvrad, vecvradortho);
  }

  // do vecB
  FTYPE vecBortho[NDIM];
  concovtype=1; // contravariant
  vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vecB, vecBortho);

  int vectorcomponentlocal;

  // compute FEMrad
  FTYPE FEMrad;
  vectorcomponentlocal=1; // radial component
  vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponentlocal, vecv, vecB, &FEMrad);

  // compute FEMrad2
  FTYPE FEMrad2;
  vectorcomponentlocal=1; // radial component
  FTYPE vecvboosted[NDIM];
  vboost(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponentlocal, vecv, vecB, vecvboosted);
  vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponentlocal, vecvboosted, vecB, &FEMrad2);


  // get dF^{\mu\nu}
  FTYPE Mcon[NDIM][NDIM],Mconorthocart[NDIM][NDIM];
  Mcon_calc_uB(vecv, vecB, Mcon);
  concovtype=1; // M^{\mu\nu} form
  T2Tortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, Mcon, Mconorthocart);
  //  FTYPE vtemp[NDIM],vtemportho[NDIM];
  int kk;
//  for(kk=0;kk<NDIM;kk++){  // bit of overkill, since gets tetrad 4X as many times as required
    //    DLOOPA(jj) vtemp[jj]=Mcon[jj][kk];
  //    vec2vecortho(concovtype,V,gcov,dxdxp,oldgridtype, newgridtype, vtemp, vtemportho);
  //    DLOOPA(jj) Mconorthocart[jj][kk]=vtemportho[jj];
  //  }

  FTYPE r=V[1];
  FTYPE h=V[2];
  FTYPE ph=V[3];
  FTYPE R=fabs(r*sin(h));


  // *F^{it} = B^i
  FTYPE Bx=Mconorthocart[1][0];
  FTYPE By=Mconorthocart[2][0];
  FTYPE Bz=Mconorthocart[3][0];

  FTYPE Br=sin(h)*cos(ph)*Bx + sin(h)*sin(ph)*By + cos(h)*Bz;
  FTYPE Bh=cos(h)*cos(ph)*Bx + cos(h)*sin(ph)*By - sin(h)*Bz;
  FTYPE Bp=-sin(ph)*Bx + cos(ph)*By;

  FTYPE Bdx=-Mconorthocart[1][0];
  FTYPE Bdy=-Mconorthocart[2][0];
  FTYPE Bdz=-Mconorthocart[3][0];

  FTYPE Bdr=sin(h)*cos(ph)*Bdx + sin(h)*sin(ph)*Bdy + cos(h)*Bdz;
  FTYPE Bdh=cos(h)*cos(ph)*Bdx + cos(h)*sin(ph)*Bdy - sin(h)*Bdz;
  FTYPE Bdp=-sin(ph)*Bdx + cos(ph)*Bdy;


  FTYPE vx=vecvortho[1];
  FTYPE vy=vecvortho[2];
  FTYPE vz=vecvortho[3];

  FTYPE vr=sin(h)*cos(ph)*vx + sin(h)*sin(ph)*vy + cos(h)*vz;
  FTYPE vh=cos(h)*cos(ph)*vx + cos(h)*sin(ph)*vy - sin(h)*vz;
  FTYPE vp=-sin(ph)*vx + cos(ph)*vy;


  // Get F_{munu}
  FTYPE Fcovorthocart[NDIM][NDIM];
  MtoF_simple(0,Mconorthocart,1.0,Fcovorthocart); // gdet=1

  // Get \Omega_F in coordinate (e.g. BH) frame
  FTYPE omegaf1b = (vp/R) - (Bp/R)*(vr*Br + vh*Bh)/(Br*Br + Bh*Bh);


  // compute B_{\phi}
  FTYPE Bphi;
  FTYPE vecBdown[NDIM];
  if(0){
    // not true B_\phi since time component still up
    vectorcomponentlocal=3; // \phi component
    vecup2vecdown(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vecB, vecBdown);
    Bphi=vecBdown[vectorcomponentlocal];
  }
  else{
    // \dF^{it} = B^i
    // dF_{it}= B_i  B_\phi = Bphi
    // orthonormal so flat space-time and only time comonent going down gives -1
    //    Bphi=-sin(ph)*Bdx + cos(ph)*Bdy;

    // poloidal current
    Bphi=R*Bdp;

  }
  FTYPE truevecBortho[NDIM];
  DLOOPA(jj) truevecBortho[jj] = Mconorthocart[jj][0];

  // T^\mu_nu in coordinate basis
  FTYPE Tmunu[NDIM][NDIM];
  vB2poyntingdensityfull(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vecv, vecB, Tmunu);

  // setup flat stuff
  FTYPE gcovflat[SYMMATRIXNDIM]={-1,0,0,0,1,0,0,1,0,1},gconflat[SYMMATRIXNDIM]={-1,0,0,0,1,0,0,1,0,1};
  FTYPE gdetflat=1.0;
  FTYPE dxdxpflat[NDIM][NDIM]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

  // T^\mu_nu in Cartesian orthonormal basis
  FTYPE Tmunuorthocart[NDIM][NDIM];
  vB2poyntingdensityfull(ti,X,V,conn,gconflat,gcovflat,gdetflat,ck,dxdxpflat,oldgridtype, newgridtype, vecvortho, truevecBortho, Tmunuorthocart);

  // get boosted T^\mu_nu in cart (T^{txyz}_{txyz}) and first index SPC and second cart (T^{trhp}_{txyz})
  FTYPE Tmunuboostedcart[NDIM][NDIM];
  FTYPE Tmunuboostedspccart[NDIM][NDIM];
  FTYPE Tmunuboostedspcspc[NDIM][NDIM];
  Tboost(ti,  X,  V,  conn,  gconflat,  gcovflat,  gdetflat,  ck,  dxdxpflat, oldgridtype, newgridtype, Tmunuorthocart, Tmunuboostedcart, Tmunuboostedspccart, Tmunuboostedspcspc);

  // T^r_x
  FTYPE Fmomx=Tmunuboostedspccart[1][1];



  // get boosted M^\mu_nu in cart (M^{txyz}_{txyz}) and first index SPC and second cart (M^{trhp}_{txyz})
  FTYPE Mconboostedcart[NDIM][NDIM];
  FTYPE Mconboostedspccart[NDIM][NDIM];
  FTYPE Mconboostedspcspc[NDIM][NDIM];
  Tboost(ti,  X,  V,  conn,  gconflat,  gcovflat,  gdetflat,  ck,  dxdxpflat, oldgridtype, newgridtype, Mconorthocart, Mconboostedcart, Mconboostedspccart, Mconboostedspcspc);

  // B^i = *F^{it}
  FTYPE Bbr=Mconboostedspcspc[1][0];
  FTYPE Bbh=Mconboostedspcspc[2][0];
  FTYPE Bbp=Mconboostedspcspc[3][0];

  FTYPE Bdbr=-Mconboostedspcspc[1][0];
  FTYPE Bdbh=-Mconboostedspcspc[2][0];
  FTYPE Bdbp=-Mconboostedspcspc[3][0];

  // boosted poloidal current
  FTYPE Bphiboosted = R*Bdbp;


  // get boosted M^\mu_nu in cart (M^{txyz}_{txyz}) and first index SPC and second cart (M^{trhp}_{txyz})
  FTYPE Fcovboostedcart[NDIM][NDIM];
  FTYPE Fcovboostedspccart[NDIM][NDIM];
  FTYPE Fcovboostedspcspc[NDIM][NDIM];
  Tboost(ti,  X,  V,  conn,  gconflat,  gcovflat,  gdetflat,  ck,  dxdxpflat, oldgridtype, newgridtype, Fcovorthocart, Fcovboostedcart, Fcovboostedspccart, Fcovboostedspcspc);

  // E_i = F_{it}
  FTYPE Ebr=Fcovboostedspcspc[1][0];
  FTYPE Ebh=Fcovboostedspcspc[2][0];
  FTYPE Ebp=Fcovboostedspcspc[3][0];

  // E^i
  FTYPE Eubr=-Fcovboostedspcspc[1][0];
  FTYPE Eubh=-Fcovboostedspcspc[2][0];
  FTYPE Eubp=-Fcovboostedspcspc[3][0];

  // Nielson et al. version.  Sam says it's  E_theta^2 + B_theta^2 or r^2 |phi_2-phi_0|^2 = F_{t theta}^2 + (F_{r phi} / sin theta)^2
  FTYPE FEMrad5 = r*r*(Ebh*Ebh + Bbh*Bbh);


  // get velocity in SPC that is the boosted frame
  FTYPE vecvboost[NDIM];
  getvboost(ti,  X,  V,  conn,  gcon,  gcov,  gdet,  ck,  dxdxp, oldgridtype, newgridtype, vecvboost);

  // get boosted \Omega_F using 3-vel and 3-field
  FTYPE vbr=vecvboost[1]/vecvboost[0];
  FTYPE vbh=vecvboost[2]/vecvboost[0];
  FTYPE vbp=vecvboost[3]/vecvboost[0];
  FTYPE omegaf1bboosted= (vbp/R) - (Bbp/R)*(vbr*Bbr + vbh*Bbh)/(Bbr*Bbr + Bbh*Bbh);
  


  // compute FEMrad3 and FEMrad4
  FTYPE FEMrad3,FEMrad4;
  // feed in cart-ortho vecv and true cart-ortho vecB
  vB2poyntingdensityboost(ti,X,V,conn,gconflat,gcovflat,gdetflat,ck,dxdxpflat,oldgridtype, newgridtype, Tmunuboostedcart, Tmunuboostedspccart, &FEMrad3, &FEMrad4);
  if(r<3.2){
    FEMrad3=FEMrad; // avoid boosting near BH
    FEMrad4=FEMrad; // avoid boosting near BH
  }

  
  



  // 3-velocity magnitude: v^2=vx^2+vy^2+vz^2 (i.e. just spatials are squared)
  // Lorentz factor: ut = 1/sqrt(1-v^2)
  // 4-velocity: u={ut,ut*vx,ut*vy,ut*vz)

  // v^2 = u^i u^i / (u^t u^t)
  FTYPE vsq=(vecvortho[RR]*vecvortho[RR] + vecvortho[TH]*vecvortho[TH] + vecvortho[PH]*vecvortho[PH])/vecvortho[TT]*vecvortho[TT];
  if(vsq<0.0) vsq=0.0;
  if(vsq>1.0-1E-10) vsq=1.0-1E-10;


  // B^2
  FTYPE Bsq=0.0;
  DLOOPA(jj) Bsq += vecBortho[jj];

  // field along flow: u.B = ux*Bx + uy*By + uz*Bz (i.e. Bt=0)
  FTYPE udotB=0.0;
  DLOOPA(jj) udotB += vecBortho[jj]*vecvortho[jj]; // B.u (B 3-frame magfield and u 4-vel)
  // comoving magnetic field: b^\mu = (B^\mu + (u.B)u^\mu)/ut
  FTYPE vecBorthoco[NDIM];
  DLOOPA(jj) vecBorthoco[jj] = (vecBortho[jj] + udotB*vecvortho[jj])/vecvortho[TT];
  // comoving mag energy: b^2/2 = 0.5*(B^2 + (u.B)^2)/ut^2
  //  FTYPE bsq = (Bsq + udotB*udotB)/(vecvortho[TT]*vecvortho[TT]); // using B,u, but same as using b directly
  FTYPE bsq=0.0;
  DLOOPA(jj) bsq += vecBorthoco[jj]*vecBorthoco[jj]; // directly


  FTYPE rho=val[FLRHO];
  FTYPE ug=val[FLU];
  FTYPE uu0ortho=vecvortho[TT];

  FTYPE prad0=val[FLURAD0];

#if(0)
  // not necessary for Sasha models
  // fix-up jet region density
  if(bsq/rho>30.0){
    rho=ug=1E-10;
  }
#endif


  FTYPE Rcyl=V[1]*sin(V[2]);


  ///////////////////
  // set outputs
  if(outputvartypelocal==15){
    fvar[0]=rho;
    fvar[1]=ug;
    fvar[2]=uu0ortho;
    fvar[3]=bsq;
  }
  else if(outputvartypelocal==16){
    fvar[0]=rho;
    fvar[1]=ug;
    fvar[2]=uu0ortho;
    fvar[3]=bsq;
    fvar[4]=log10(rho);
    fvar[5]=-log10(rho);
    fvar[6]=log10(bsq);
    fvar[7]=Rcyl;
  }
  else if(outputvartypelocal==17){
    fvar[0]=rho;
    fvar[1]=ug;
    fvar[2]=uu0ortho;
    //    fvar[3]=bsq;
    fvar[3]=MIN(bsq/rho,1E2);
    fvar[4]=log10(rho);
    fvar[5]=-log10(rho);
    fvar[6]=log10(bsq);
    fvar[7]=Rcyl;
    fvar[8]=vecvortho[1]/vecvortho[TT];
    fvar[9]=vecvortho[2]/vecvortho[TT];
    fvar[10]=vecvortho[3]/vecvortho[TT];
    fvar[11]=vecBortho[1];
    fvar[12]=vecBortho[2];
    fvar[13]=vecBortho[3];
  }
  else if(outputvartypelocal==18){
    fvar[0]=rho; // rest-mass density
    fvar[1]=ug;  // internal energy density
    fvar[2]=uu0ortho; // Lorentz factor
    //    fvar[3]=bsq;
    fvar[3]=MIN(bsq/rho,1E2); // limited bsq/rho
    fvar[4]=log10(rho);
    fvar[5]=-log10(rho);
    fvar[6]=log10(bsq);
    fvar[7]=Rcyl; // Cylindrical radius
    fvar[8]=vecvortho[1]/vecvortho[TT]; // vx
    fvar[9]=vecvortho[2]/vecvortho[TT]; // vy
    fvar[10]=vecvortho[3]/vecvortho[TT]; // vz
    fvar[11]=vecBortho[1]; // Bx
    fvar[12]=vecBortho[2]; // By
    fvar[13]=vecBortho[3]; // Bz
    fvar[14]=V[1]; // spherical polar radius r
    fvar[15]=V[2]; // spherical polar angle \theta = 0..\pi
    fvar[16]=V[3]; // spherical polar angle \phi = 0..2\pi
    fvar[17]=V[1]*sin(V[2])*cos(V[3]); // x
    fvar[18]=V[1]*sin(V[2])*sin(V[3]); // y
    fvar[19]=V[1]*cos(V[2]); // z

  }
  else if(outputvartypelocal==19){
    fvar[0]=rho; // rest-mass density
    fvar[1]=ug;  // internal energy density
    fvar[2]=uu0ortho; // Lorentz factor
    //    fvar[3]=bsq;
    fvar[3]=MIN(bsq/rho,1E2); // limited bsq/rho
    fvar[4]=log10(rho);
    fvar[5]=-log10(rho);
    fvar[6]=log10(bsq);
    fvar[7]=Rcyl; // Cylindrical radius
    fvar[8]=vecvortho[1]/vecvortho[TT]; // vx
    fvar[9]=vecvortho[2]/vecvortho[TT]; // vy
    fvar[10]=vecvortho[3]/vecvortho[TT]; // vz
    fvar[11]=truevecBortho[1];//vecBortho[1]; // Bx
    fvar[12]=truevecBortho[2];//vecBortho[2]; // By
    fvar[13]=truevecBortho[3];//vecBortho[3]; // Bz
    fvar[14]=V[1]; // spherical polar radius r
    fvar[15]=V[2]; // spherical polar angle \theta = 0..\pi
    fvar[16]=V[3]; // spherical polar angle \phi = 0..2\pi
    fvar[17]=prad0;
    fvar[18]=vecvradortho[1]/vecvradortho[TT]; // vradx
    fvar[19]=vecvradortho[2]/vecvradortho[TT]; // vrady
    fvar[20]=vecvradortho[3]/vecvradortho[TT]; // vradz
    fvar[21]=FEMrad; // FEMrad = lab-frame radial energy flux: T^r_t
    fvar[22]=FEMrad2; // FEMrad using boosted velocity
    fvar[23]=FEMrad3; // boosted FEMrad
    fvar[24]=FEMrad4; // Sam's FEMrad
    fvar[25]=FEMrad5; // Rezzolla FEMrad
    fvar[26]=Fmomx; // T^r_x
    fvar[27]=Bphi; // *F_{\phi t} = B_{\phi} == poloidal current
    fvar[28]=Bphiboosted; // *F_{\phi t} = B_{\phi} == poloidal current boosted
    fvar[29]=omegaf1b; // \Omega_F
    fvar[30]=omegaf1bboosted; // \Omega_F boosted
  }
  else{
    dualfprintf(fail_file,"No such outputvartypelocal=%d\n",outputvartypelocal);
    myexit(1);
  }




  return(0);
}









// compute additional (or all) things in cases when need full temporal-spatial information to compute somthing.  In that case, often just compute everything here.
int compute_additionals(void)
{
  int colini,h,i,j,k;
  int jj; // for tensor indices
  int kprior;
  int colstorei;





  // default to indicate nothing done
  olddata3time=olddatagdump=olddatacurrent=NULL;


  /////////////
  // SETUP 3-TIME DATA READ
  int doing3time=0;
  if(infilem1!=NULL && infilep1!=NULL && DATATYPE==14 && outputvartype!=0 && docurrent==1){
    doing3time=1;
  }



  if(doing3time){
    fprintf(stderr,"BEGIN Computing additionals (currently only Jcon, which requires spatio-temporal derivatives)\n"); fflush(stderr);

    int doubleworkfake;




    //////////////
    //
    // SETUP MEMORY
    //
    ///////////////

    // then store full data for 3 times at once
    int oN0fake=3; // 3-time
    int numbc0fake=0; // no boundary conditions
    // for storing (only required) fieldline quantities
    // only require rho0,ug,u^\mu,B^i
    // olddata3time NUMCOLUMNSSTORE related to colstorei loops
    //
    olddata3time = f5matrix(0,NUMCOLUMNSSTORE-1,-numbc0fake,oN0fake-1+numbc0fake,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;   // olddata3time[colini][h][i][j][k]
    if(olddata3time==NULL){
      fprintf(stderr,"Couldn't allocate olddata3time\n"); fflush(stderr);
      myexit(1);
    }

    // only for storing gcov(SYMMATRIXNDIM) + gdet(1) + V(NDIM-1) + dxdxp((NDIM-2)*(NDIM-2))
    // NOTEMARK: only store V[spatial] since V[0] fixed for now.  Also only store dxdxp[1,2] since dxdx[0,3] are fixed.
    int numgdumps=SYMMATRIXNDIM  + 1 + (NDIM-1) + (NDIM-2)*(NDIM-2);
    olddatagdump = f5matrix(0,numgdumps-1,0,0,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ; // no temporal gdump
    if(olddatagdump==NULL){
      fprintf(stderr,"Couldn't allocate olddatagdump\n"); fflush(stderr);
      myexit(1);
    }

    olddatacurrent = f5matrix(0,NDIM-1,0,0,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ; // no temporal info required
    if(olddatacurrent==NULL){
      fprintf(stderr,"Couldn't allocate olddatacurrent\n"); fflush(stderr);
      myexit(1);
    }



    // old column space (i.e numcolumns with colini loops)
    FTYPE *dataorigin=NULL;
    dataorigin=(FTYPE*)malloc((unsigned)(numcolumns)*sizeof(FTYPE));
    if(dataorigin==NULL){
      fprintf(stderr,"Couldn't allocate dataorigin\n");
      exit(1);
    }
    

    // temp vars for gdump
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
    FTYPE val[MAXINCOLS]; // to hold numcolumns - coli type data (even if have to set a value to zero if not used, keep in fieldline file format)


    if(MAXINCOLS<numcolumns || MAXINCOLS<NUMCOLUMNSSTORE){
      dualfprintf(fail_file,"Not enough MAXINCOLS=%d for given numcolumns=%d or NUMCOLUMNSSTORE=%d\n",MAXINCOLS,numcolumns,NUMCOLUMNSSTORE);
      myexit(1);
    }


    // fvar holds coli -- numoutputcols data (i.e. results of preprocessing)
    FTYPE *fvar=(FTYPE*)malloc((unsigned)(numoutputcols)*sizeof(FTYPE));
    if(fvar==NULL){
      fprintf(stderr,"Couldn't allocate fvar\n");
      exit(1);
    }


    if(olddata0==NULL){
      fprintf(stderr,"olddata0 should be allocated before calling compute_additionals() for DATATYPE==14\n");
      myexit(1);
    }







    //////////////
    //
    // LOOP OVER fieldline file data for each of the 3 times
    //
    ///////////////

    //for(i=1;i<=4;i++) while(fgetc(infilem1)!='\n'); // asume at this point that all 3 files have had header read if header exists
    kprior=-1000;
    fprintf(stderr,"BEGIN Reading(numcolumns=%d)/Storing (and Computing local datatype==14 stuff) 3-time data and gdump data with oN3=%d dots produced\n",numcolumns,oN3);fflush(stderr);
    //
    LOOPOLDDATASPATIAL{
      if(k!=kprior){ fprintf(stderr,"."); fflush(stderr); kprior=k; }

      // READ-IN such that columns are fastest index, then i, then j, then k as in data files
      colstorei=0;
      for(colini=0;colini<numcolumns;colini++){

        // get m1
        h=0; readelement(binaryinput,inFTYPE,infilem1,&dataorigin[colini]);
        if(!SKIPFL2STORE(colini)) olddata3time[colstorei][h][i][j][k]=dataorigin[colini];

        // get normal middle time value
        h=1; readelement(binaryinput,inFTYPE,infile,&dataorigin[colini]);
        if(!SKIPFL2STORE(colini)) olddata3time[colstorei][h][i][j][k]=dataorigin[colini];
        val[colini]=dataorigin[colini]; // full store for compute_datatype14(1) to use (store middle time)

        // get p1
        h=2; readelement(binaryinput,inFTYPE,infilep1,&dataorigin[colini]);
        if(!SKIPFL2STORE(colini)) olddata3time[colstorei][h][i][j][k]=dataorigin[colini];


        // only store required information, so only iterate colstorei if storing this quantity
        if(!SKIPFL2STORE(colini)) colstorei++;
      }

      // check that stored correct number of things
      if(colstorei!=NUMCOLUMNSSTORE){
        fprintf(stderr,"Didn't properly store olddata3time: %d %d\n",colstorei,NUMCOLUMNSSTORE);
        myexit(1);
      }

      //DEBUG
      //      fprintf(stderr,"PRECOMPUTE0: %d %d %d\n",i,j,k);
      //      for(colstorei=0;colstorei<NUMCOLUMNSSTORE;colstorei++) fprintf(stderr,"PRECOMPUTE1: %d : %g %g %g\n",colstorei,olddata3time[colstorei][0][i][j][k],olddata3time[colstorei][1][i][j][k],olddata3time[colstorei][2][i][j][k]);

      //////////////
      //
      // GET gdump data for "normal" time (same for all times currently)
      //
      //////////////
      read_gdumpline(gdumpin, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp, &geom);




      // SUPERGODMARK: TODOMARK OPTMARK: Can speed things up alot if know interpolation box won't use portion of original data, which can then be skipped.  Can use r,\theta,\phi to detect (e.g. for Cart output) if need certain data.  Just compute x,y,z and see if inside new box.  If not, skip.  Will only generally exclude radial parts, and so maybe save factor of 2-4 or so.  Should do it!



      ///////////////////
      //
      // COMPUTE PER-POINT DATATYPE==14
      //
      ///////////////////
      // so that don't have to reload data or gdump file for better speed
      /////////////// 
      // 1 means compute everything but Jconortho
      compute_datatype14(1, val, fvar, ti,  X,  V,  conn,  gcon,  gcov,  gdet,  ck,  dxdxp, &geom);
      // store ONLY non-current results
      h=0; int coli; for(coli=OUTRHO;coli<=OUTBPHI;coli++) olddata0[coli][h][i][j][k]=fvar[coli];
      

      // put required elements into storage for J^\mu computation
      h=0; int gdumpjj=0;
      for(jj=0;jj<SYMMATRIXNDIM;jj++){ olddatagdump[gdumpjj][h][i][j][k]=gcov[jj]; gdumpjj++;}
      olddatagdump[gdumpjj][h][i][j][k]=gdet; gdumpjj++;
      SLOOPA(jj){ olddatagdump[gdumpjj][h][i][j][k]=V[jj]; gdumpjj++; }
      int lll,mmm; SLOOP12(lll,mmm){ olddatagdump[gdumpjj][h][i][j][k]=dxdxp[lll][mmm];  gdumpjj++;}
      if(gdumpjj!=numgdumps){
        fprintf(stderr,"Memory issue with olddatagdump assignment\n");
        myexit(1);
      }
      // now have all that's required to compute F^{\mu\nu} locally for this grid point
      // and also have enough to compute orthonormal versions of things (e.g. for J^\mu later)
    }
    fprintf(stderr,"\nEND Reading/Storing 3-time data and gdump data with oN3=%d dots produced\n",oN3);fflush(stderr);



    //////////////
    //
    // Apply boundary conditions
    //
    ///////////////

    // apply boundary conditions on this 3-time data
    doubleworkfake=1; // just forces use of olddata3time
    apply_boundaryconditions_olddata(NUMCOLUMNSSTORE,oN0fake,numbc0fake,doubleworkfake,oldimage0,olddata3time);

    // apply boundary conditions to gdump data
    doubleworkfake=1; // just forces use of olddatagdump
    apply_boundaryconditions_olddata(numgdumps,1,0,doubleworkfake,oldimage0,olddatagdump);



    //////////////
    //
    // Store cubical data and then compute space-time derivatives to get current
    //
    ///////////////
    

    // compute F^{\mu\nu} with only local cube of memory
    FTYPE uconcube[3][3][3][3][NDIM];
    FTYPE Bconcube[3][3][3][3][NDIM];
    FTYPE gcovcube[3][3][3][3][SYMMATRIXNDIM];
    FTYPE gdetcube[3][3][3][3];
    FTYPE gdetFuucube[3][3][3][3][NDIM][NDIM];
    FTYPE Jcon[NDIM];

    // can order sub-loop however one wants, as long as memory access to arrays is correct.  Order for faster memory access (k fastest index)
    //#define SUBLOOPOLDDATA for(hhh=0;hhh<3;hhh++)  for(iii=0;iii<3;iii++) for(jjj=0;jjj<3;jjj++) for(kkk=0;kkk<3;kkk++)
    // limit to only required positions
#define SUBLOOPOLDDATA for(hhh=0;hhh<3;hhh++)  for(iii=0;iii<3;iii++) for(jjj=0;jjj<3;jjj++) for(kkk=0;kkk<3;kkk++)

    kprior=-1000;
    fprintf(stderr,"BEGIN Computing Jcon with oN3=%d dots produced\n",oN3);fflush(stderr);
    //
    int hhh,iii,jjj,kkk;
    h=1; // middle time from 3-time data is reference position where current is eventually located effectively
    int hgdump=0; // gdump has no temporal information (yet)
    LOOPOLDDATASPATIAL{ // 3D loop. In the end, don't need temporal information and just use +-1 values offset from normal dataset


      // SUPERGODMARK: TODOMARK OPTMARK: Can speed things up alot if know interpolation box won't use portion of original data, which can then be skipped.  Can use r,\theta,\phi to detect (e.g. for Cart output) if need certain data.  Just compute x,y,z and see if inside new box.  If not, skip.  Will only generally exclude radial parts, and so maybe save factor of 2-4 or so.  Should do it!


      if(k!=kprior){ fprintf(stderr,"."); fflush(stderr); kprior=k; }

      SUBLOOPOLDDATA{ // full 4D loop
 
        // restrict loop to only required points for computing derivatives
        if(abs(hhh-1)+abs(iii-1)+abs(jjj-1)+abs(kkk-1)>1) continue; // i.e. only 1 direction offset at a time.

 
        // u^\mu
        DLOOPA(jj) uconcube[hhh][iii][jjj][kkk][jj]=olddata3time[STOREU0+jj][h+hhh-1][i+iii-1][j+jjj-1][k+kkk-1];
        SLOOPA(jj) uconcube[hhh][iii][jjj][kkk][jj]*=uconcube[hhh][iii][jjj][kkk][TT];
        // B^\mu
        Bconcube[hhh][iii][jjj][kkk][TT]=0.0;
        SLOOPA(jj) Bconcube[hhh][iii][jjj][kkk][jj]=olddata3time[STOREB1+jj-1][h+hhh-1][i+iii-1][j+jjj-1][k+kkk-1];
        // g_{\mu\nu} for SYMMATRIXNDIM unique elements
        int gdumpjj=0;
        for(jj=0;jj<SYMMATRIXNDIM;jj++){ gcovcube[hhh][iii][jjj][kkk][jj]=olddatagdump[gdumpjj][hgdump][i+iii-1][j+jjj-1][k+kkk-1]; gdumpjj++;}
        // \detg = \sqrt{-g}
        gdetcube[hhh][iii][jjj][kkk]=olddatagdump[gdumpjj][hgdump][i+iii-1][j+jjj-1][k+kkk-1]; gdumpjj++;

        // DEBUG:
        //fprintf(stderr,"Before compute_gdetFuu: %d %d %d %d : %d %d %d %d\n",h,i,j,k,hhh,iii,jjj,kkk);fflush(stderr);

        // \detg F^{\mu\nu}
        compute_simple_gdetFuu(gdetcube[hhh][iii][jjj][kkk], gcovcube[hhh][iii][jjj][kkk], uconcube[hhh][iii][jjj][kkk], Bconcube[hhh][iii][jjj][kkk], gdetFuucube[hhh][iii][jjj][kkk]);


        // DEBUG:
#if(0)
        fprintf(stderr,"COMPUTE0: %d %d %d %d : %d %d %d %d \n",h,i,j,k,hhh,iii,jjj,kkk);
        for(jj=0;jj<SYMMATRIXNDIM;jj++) fprintf(stderr,"COMPUTE1: %d %g\n",jj,gcovcube[hhh][iii][jjj][kkk][jj]);
        fprintf(stderr,"COMPUTE2:  %g\n",gdetcube[hhh][iii][jjj][kkk]);
        DLOOPA(jj) fprintf(stderr,"COMPUTE3: %d  %g\n",jj,uconcube[hhh][iii][jjj][kkk][jj]);
        DLOOPA(jj) fprintf(stderr,"COMPUTE4: %d  %g\n",jj,Bconcube[hhh][iii][jjj][kkk][jj]);
        int kk; DLOOP(jj,kk) fprintf(stderr,"COMPUTE5: %d %d %g\n",jj,kk,gdetFuucube[hhh][iii][jjj][kkk][jj][kk]);
#endif


        // DEBUG:
        //fprintf(stderr,"After compute_gdetFuu\n");fflush(stderr);

      }// over hhh,iii,jjj,kkk

      // now have \detg F^{\mu\nu} over spatio-temporal cube
      // So can compute J^\mu = (1/\detg) (\detg F^{\mu\nu})_{,\nu}  (sign won't matter in the end when get J^2)
      // note that only using grid-aligned derivatives using centered difference, which is as accurate as parabolae because that term cancels.
      FTYPE fakeDT = endtdata0 - starttdata0; // because dX[TT]=dt=1 when oN0=1 and avoid using endtdata and starttdata that get modified since used for another purpose (true 4D data)
      
      Jcon[TT] =
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][2][1][1][TT][RR] - gdetFuucube[1][0][1][1][TT][RR])/dX[RR]
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][2][1][TT][TH] - gdetFuucube[1][1][0][1][TT][TH])/dX[TH]
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][1][2][TT][PH] - gdetFuucube[1][1][1][0][TT][PH])/dX[PH]
        ;

      Jcon[RR] =
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[2][1][1][1][RR][TT] - gdetFuucube[0][1][1][1][RR][TT])/fakeDT
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][2][1][RR][TH] - gdetFuucube[1][1][0][1][RR][TH])/dX[TH]
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][1][2][RR][PH] - gdetFuucube[1][1][1][0][RR][PH])/dX[PH]
        ;

      Jcon[TH] =
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[2][1][1][1][TH][TT] - gdetFuucube[0][1][1][1][TH][TT])/fakeDT
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][2][1][1][TH][RR] - gdetFuucube[1][0][1][1][TH][RR])/dX[RR]
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][1][2][TH][PH] - gdetFuucube[1][1][1][0][TH][PH])/dX[PH]
        ;

      Jcon[PH] =
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[2][1][1][1][PH][TT] - gdetFuucube[0][1][1][1][PH][TT])/fakeDT
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][2][1][1][PH][RR] - gdetFuucube[1][0][1][1][PH][RR])/dX[RR]
        +(1.0/gdetcube[1][1][1][1])*(gdetFuucube[1][1][2][1][PH][TH] - gdetFuucube[1][1][0][1][PH][TH])/dX[TH]
        ;


      // DEBUG:
      //      DLOOPA(jj) fprintf(stderr,"Jcon[%d]=%g : Dx=%g fakeDT=%g gdet=%g\n",jj,Jcon[jj],dX[jj],fakeDT,gdetcube[1][1][1][1]); fflush(stderr);
          


      int href=0;
      DLOOPA(jj) olddatacurrent[jj][href][i][j][k]=Jcon[jj];
      // now have current stored in global "olddatacurrent" location.  This data can be used just like other single-point versions of data.

      // V^i (assume V^t=t didn't change from older assignment during gdump read-in since constant)
      int gdumpjj=SYMMATRIXNDIM+1;
      SLOOPA(jj){ V[jj]=olddatagdump[gdumpjj][href][i][j][k]; gdumpjj++; }
      // dx^\mu/dxp^\nu (only modify 1-2 components, since rest didn't change from assignment of constant from gdump read-in in prior loop)
      int lll,mmm; SLOOP12(lll,mmm){ dxdxp[lll][mmm]=olddatagdump[gdumpjj][href][i][j][k]; gdumpjj++; }
      if(gdumpjj!=numgdumps){
        fprintf(stderr,"Memory issue with olddatagdump assignment\n");
        myexit(1);
      }


      // setup ti (no need to store, so just regenerate from loop)
      ti[TT]=0;ti[RR]=i;ti[TH]=j;ti[PH]=k;
      // 2 means compute Jconortho
      // "val" below is not filled but not used.
      compute_datatype14(2, val, fvar, ti,  NULL,  V,  NULL,  NULL,  gcovcube[1][1][1][1],  gdetcube[1][1][1][1],  NULL,  dxdxp, NULL); // don't really need gdet
      // store ONLY current results (so don't overwrite non-current results with non-computed fvar)
      int coli; for(coli=OUTJX0;coli<=OUTJX3;coli++) olddata0[coli][href][i][j][k]=fvar[coli];


      
    } // over h,k,j,i
    fprintf(stderr,"\nEND Computing Jcon with oN3=%d dots produced\n",oN3);fflush(stderr);


    // for thickdisk7 using DATATYPE==14, without current calculation, takes only 750MB
    // with currents, after reducing to only needed things, now 4.4GB resident.
    // with DATATYPE==14 and docurrent==1, code takes 5 minutes for Phase1 (reading/storing), then taking 16 minutes (5:00 - 21:00) for Phease 2 (computing Jcon), then quickly does Phase 3 (re-storing data) in <1 minute (21:00 - 21:30), and then interpolation at 100^3 takes 2 more minutes (21:30 - 22:45):  Total time: 23 minutes
    // So since computing Jcon is so expensive, might as well do 512^3 that takes about 20 minutes so total time is about 40 minutes.
    //
    // density only 100^3: 2 minutes
    // no current 100^3:  4minutes preprocess + 2 minutes interp = 6  minutes total
    // no current 512^3:  4minutes preprocess +   minutes interp =    minutes total
    //    current 100^3: 22minutes preprocess + 1 minutes interp = 23 minutes total
    //    current 512^3: 22minutes preprocess +   minutes interp =    minutes total


    /////////////////////////
    // free temp space
    //
    fprintf(stderr,"BEGIN Freeing space used by computing additionals\n"); fflush(stderr);
    fprintf(stderr,"olddata3time: %d %d %d %d %d %d %d %d %d %d",0,NUMCOLUMNSSTORE-1,-numbc0fake,oN0fake-1+numbc0fake,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]);
    fprintf(stderr,"olddatagdump: %d %d %d %d %d %d %d %d %d %d",0,numgdumps-1,0,0,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]);

    //    olddata3time = f5matrix(0,NUMCOLUMNSSTORE-1,-numbc0fake,oN0fake-1+numbc0fake,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ;    
    //    olddatagdump = f5matrix(0,numgdumps-1,0,0,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]) ; // no temporal gdump

    free_f5matrix(olddata3time,0,NUMCOLUMNSSTORE-1,-numbc0fake,oN0fake-1+numbc0fake,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]);
    free_f5matrix(olddatagdump,0,numgdumps-1,0,0,-numbc[1]+0,oN1-1+numbc[1],-numbc[2]+0,oN2-1+numbc[2],-numbc[3]+0,oN3-1+numbc[3]);
  
    fprintf(stderr,"DONE Freeing space used by computing additionals\n"); fflush(stderr);


    fprintf(stderr,"BEGIN rewinds\n"); fflush(stderr);

    // now rewind fieldline files and gdump files and pass header if exists, so starting at normal starting point as if this routine were never called
    infile_tostartofdata(&infile);
    // no need to rewind infilem1 and infilep1
    gdump_tostartofdata(&gdumpin);

    fprintf(stderr,"END rewinds\n"); fflush(stderr);


    fprintf(stderr,"END Computing additionals\n"); fflush(stderr);
  }// end if reading-in special 3-time data set



 

  return(0);
}








// compute \detg F^{\mu\nu}
void compute_gdetFuu(FTYPE gdet, FTYPE *gcov, FTYPE *ucon, FTYPE *Bcon, FTYPE (*Fuu)[NDIM])
{
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM],bcov[NDIM],bsq;
  FTYPE Bcov[NDIM];
  struct of_geom geom;
  int jj,kk;


  /// assign gdet
  geom.gdet=gdet;

  // compute lower components
  DLOOPA(jj) ucov[jj]=0.0;
  DLOOP(jj,kk) ucov[jj] += ucon[kk]*gcov[GIND(jj,kk)];
  
  // make primitive u^i
  FTYPE pr[NPR];
  pr[U1]=ucon[RR];
  pr[U2]=ucon[TH];
  pr[U3]=ucon[PH];
  pr[B1]=Bcon[RR];
  pr[B2]=Bcon[TH];
  pr[B3]=Bcon[PH];

  // compute b^\mu[pr,ucon,ucov]
  bcon_calc(pr, ucon, ucov, bcon);
  
  // compute b_\mu
  DLOOPA(jj) bcov[jj]=0.0;
  DLOOP(jj,kk) bcov[jj] += bcon[kk]*gcov[GIND(jj,kk)];
  
  // compute "fake" B_\mu (i.e. didn't lower upper t in *F^{t\mu} )
  DLOOPA(jj) Bcov[jj]=0.0;
  DLOOP(jj,kk) Bcov[jj] += Bcon[kk]*gcov[GIND(jj,kk)];
  
  // compute b^2 (not used yet)
  bsq = dot(bcon,bcov);

  // now get F^{\mu\nu}
  int updown=1;
  faraday_calc(updown,bcov,ucov,&geom,Fuu);


  DLOOP(jj,kk) Fuu[jj][kk]*=gdet;



}



// compute \detg F^{\mu\nu}
void compute_simple_gdetFuu(FTYPE gdet, FTYPE *gcov, FTYPE *ucon, FTYPE *Bcon, FTYPE (*Fuu)[NDIM])
{
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM],bcov[NDIM],bsq;
  FTYPE Bcov[NDIM];
  struct of_geom geom;
  int jj,kk;


  /// assign gdet
  geom.gdet=gdet;

  // compute lower components
  DLOOPA(jj) ucov[jj]=0.0;
  DLOOP(jj,kk) ucov[jj] += ucon[kk]*gcov[GIND(jj,kk)];
  
  // make primitive u^i
  FTYPE pr[NPR];
  pr[U1]=ucon[RR];
  pr[U2]=ucon[TH];
  pr[U3]=ucon[PH];
  pr[B1]=Bcon[RR];
  pr[B2]=Bcon[TH];
  pr[B3]=Bcon[PH];

  // compute b^\mu[pr,ucon,ucov]
  bcon_calc(pr, ucon, ucov, bcon);

  // compute b_\mu
  DLOOPA(jj) bcov[jj]=0.0;
  DLOOP(jj,kk) bcov[jj] += bcon[kk]*gcov[GIND(jj,kk)];
  
  // now get F^{\mu\nu}
  int updown=1; // this uses bcov and ucov
  faraday_calc(updown,bcov,ucov,&geom,Fuu);


  DLOOP(jj,kk) Fuu[jj][kk]*=gdet;



}




/////////////////////////////////////////////////////
// lc4() and faraday_calc() COPIED FROM phys.tools.c
/////////////////////////////////////////////////////

// used Mathematica's MinimumChangePermutations and Signature
// updown refers to whether \epsilon^{\mu\nu\kappa\lambda} was fully up or fully down (it does not refer to the faraday)
// updown = 0 : down (i.e. \epsilon_{\alpha\beta\delta\gamma} = \sqrt{-g} [\alpha\beta\gamma\delta]
// updown = 1 : up (i.e. \epsilon^{\alpha\beta\delta\gamma} = (-1/\sqrt{-g}) [\alpha\beta\gamma\delta]
static FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda)
{
  int i;
  FTYPE lc4sign; // 1,-1,1,-1... for all 24 entires
  int l1[24]={1, 2, 3, 1, 2, 3, 4, 2, 1, 4, 2, 1, 1, 3, 4, 1, 3, 4, 4, 3, 2, 4, 3, 2};
  int l2[24]={2, 1, 1, 3, 3, 2, 2, 4, 4, 1, 1, 2, 3, 1, 1, 4, 4, 3, 3, 4, 4, 2, 2, 3};
  int l3[24]={3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 4, 4, 4, 4, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4};
  int l4[24]={4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};

  for(i=0;i<24;i++){
    if((1+mu==l1[i])&&(1+nu==l2[i])&&(1+kappa==l3[i])&&(1+lambda==l4[i])){
      lc4sign=(i%2) ? -1 : 1;
      if(updown==1) return(-1.0/detg*lc4sign); // upper epsilon
      else if(updown==0) return(detg*lc4sign); // lower epsilon
    }
  }
  // if didn't get here, then 0
  return(0.0);
}

// assume anti-symmetric
// assumes b and u are inputted as bcov&ucov for F^{\mu\nu} and bcon&ucon for F_{\mu\nu}
static void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE (*faraday)[NDIM])
{
  int nu,mu,kappa,lambda;

  // initialize
  for(nu=0;nu<NDIM;nu++){
    for(mu=0;mu<NDIM;mu++){
      faraday[mu][nu]=0.0;
    }
  }

  // get unique part (nu=0 mu=1,2,3  then nu=1 mu=2,3  then nu=2 mu=3)
  for(nu=0;nu<NDIM;nu++){
    for(mu=nu+1;mu<NDIM;mu++){
   
      faraday[mu][nu]=0.0;
      for(kappa=0;kappa<NDIM;kappa++){
        for(lambda=0;lambda<NDIM;lambda++){
 
          // faraday_calc(which) refers to whether faraday is fully up or down, and lc4(updown) refers to updown being \epsilon fully up or down.  And these are the same.
          // F^{\alpha\beta} = -b_\gamma u_\delta \epsilon^{\alpha\beta\gamma\delta}
          // F^{\mu\nu} = \epsilon^{\mu\nu\kappa\lambda} u_\kappa b_\lambda
          // So sign below is correct.
          faraday[mu][nu] += lc4(which,geom->gdet,mu,nu,kappa,lambda)*u[kappa]*b[lambda];
   
          // DEBUG:
          //fprintf(stderr,"faraday[%d][%d]=%g : lc4=%g : %d %g %d %d %g %g\n",mu,nu,faraday[mu][nu],lc4(which,geom->gdet,mu,nu,kappa,lambda),which,geom->gdet,kappa,lambda,u[kappa],b[lambda]); fflush(stderr);
        }
      }
    }
  }

  // assign anti-symmetric part that's non-zero
  // nu=1 mu=0  then nu=2 mu=0,1  then nu=3 mu=0,1,2)
  for(nu=0;nu<NDIM;nu++){
    for(mu=0;mu<nu;mu++){
      faraday[mu][nu]=-faraday[nu][mu];
    }
  }

  // DEBUG:
  //  DLOOP(mu,nu) fprintf(stderr,"faraday[%d][%d]=%g\n",mu,nu,faraday[mu][nu]); fflush(stderr);


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
// concovtype: 1 = inputting u^\mu   2 = inputting u_\mu
static void vec2vecortho(int concovtype, FTYPE V[],  FTYPE *gcov,  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho)
{
  FTYPE lambdacoord[NDIM][NDIM];
  int jj,kk;
  FTYPE tetrcov[NDIM][NDIM],tetrcon[NDIM][NDIM],eigenvalues[NDIM];
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];


  // vector here is in original X coordinates
  DLOOPA(jj) finalvec[jj]=vec[jj];


  // get tetrad (uses dxdxp so that tetrcov and tetrcon and eigenvalues are using V metric not X metric)
  int primcoord=1; // tells to use dxdxp to make simpler
  tetr_func_frommetric(primcoord, dxdxp, gcov, tetrcov, tetrcon, eigenvalues);

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
      // u^kk[ortho] = tetrcov^kk[ortho]_jj[lab] u^jj[lab]
      tempcomp[kk] += tetrcov[kk][jj]*finalvec[jj];
    }
  }
  else if(concovtype==2){
    // transform to orthonormal basis for covariant vector in V coordinates
    DLOOP(jj,kk){
      // u_kk[ortho] = tetrcon_kk[ortho]^jj[lab] u_jj[lab]
      tempcomp[kk] += tetrcon[kk][jj]*finalvec[jj];
    }
  }
  else{
    dualfprintf(fail_file,"No such concovtype=%d\n",concovtype);
    myexit(1);
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];



  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d finalvec(tetrcov)=%21.15g\n",jj,finalvec[jj]);
  //  }

  // get SPC -> Cart transformation
  generate_lambdacoord(oldgridtype, newgridtype, V, lambdacoord);

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d lambdacoord=%21.15g\n",jj,kk,lambdacoord[jj][kk]);
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

// coordinate transform of vector and get single component of result
// concovtype: 1 = inputting T^{\mu\nu}   2 = inputting T_{\mu\nu} 3 = inputting T^\mu_\nu 4 = inputting T_\mu^\nu
static void T2Tortho(int concovtype, FTYPE V[],  FTYPE *gcov,  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*T)[NDIM], FTYPE (*Tortho)[NDIM])
{
  int jj,kk,ll,mm;
  FTYPE tetrcov[NDIM][NDIM],tetrcon[NDIM][NDIM],eigenvalues[NDIM];


  // get tetrad (uses dxdxp so that tetrcov and tetrcon and eigenvalues are using V metric not X metric)
  int primcoord=1; // tells to use dxdxp to make simpler
  tetr_func_frommetric(primcoord, dxdxp, gcov, tetrcov, tetrcon, eigenvalues);

  // DEBUG
  //  if(tiglobal[1]==200 && tiglobal[2]==10 && tiglobal[3]==0){
  //    //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d gcov=%21.15g\n",jj,kk,gcov[GIND(jj,kk)]);
  //    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d tetrcon=%21.15g\n",jj,kk,tetrcon[jj][kk]);
  //    DLOOPA(jj) dualfprintf(fail_file,"jj=%d eigenvalues=%21.15g\n",jj,eigenvalues[jj]);
  //  }

  FTYPE idxdxp[NDIM][NDIM];
  // get inverse coordinate transformation
  idxdxprim(dxdxp, idxdxp); // from coord.c


  // transform from X to V
  FTYPE TinV[NDIM][NDIM];
  DLOOP(jj,kk) TinV[jj][kk]=0.0;
  if(concovtype==1){
    // transform from X to V for contravariant tensor
    int ll,mm;
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        TinV[jj][kk] += dxdxp[jj][ll]*dxdxp[kk][mm]*T[ll][mm];
      }
    }
  }
  else if(concovtype==2){
    // transform from X to V for covariant tensor
    int ll,mm;
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        TinV[jj][kk] += idxdxp[ll][jj]*idxdxp[mm][kk]*T[ll][mm];
      }
    }
  }
  else if(concovtype==3){
    // transform from X to V for T^\mu_\nu
    int ll,mm;
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        TinV[jj][kk] += dxdxp[jj][ll]*idxdxp[mm][kk]*T[ll][mm];
      }
    }
  }
  else if(concovtype==4){
    // transform from X to V for T_\mu^\nu
    int ll,mm;
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        TinV[jj][kk] += idxdxp[ll][jj]*dxdxp[kk][mm]*T[ll][mm];
      }
    }
  }
  else{
    dualfprintf(fail_file,"not yet: concovtype=%d\n",concovtype);
    myexit(9834353);
  }


  // transform to orthonormal basis in SPC
  FTYPE Torthospc[NDIM][NDIM];
  DLOOP(jj,kk) Torthospc[jj][kk]=0.0;
  if(concovtype==1){
    // transform to orthonormal basis for contravariant tensor in V coordinates
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        // u^kk[ortho] = tetrcov^kk[ortho]_jj[lab] u^jj[lab]
        Torthospc[jj][kk] += tetrcov[jj][ll]*tetrcov[kk][mm]*TinV[ll][mm];
      }
    }
  }
  else if(concovtype==2){
    // transform to orthonormal basis for covariant tensor in V coordinates
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        // u_kk[ortho] = tetrcon_kk[ortho]^jj[lab] u_jj[lab]
        Torthospc[jj][kk] += tetrcon[jj][ll]*tetrcon[kk][mm]*TinV[ll][mm];
      }
    }
  }
  else if(concovtype==3){
    // transform to orthonormal basis for T^\mu_\nu in V coordinates
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        Torthospc[jj][kk] += tetrcov[jj][ll]*tetrcon[kk][mm]*TinV[ll][mm];
      }
    }
  }
  else if(concovtype==4){
    // transform to orthonormal basis for T_\mu^\nu in V coordinates
    DLOOP(jj,kk){
      DLOOP(ll,mm){
        Torthospc[jj][kk] += tetrcon[jj][ll]*tetrcov[kk][mm]*TinV[ll][mm];
      }
    }
  }
  else{
    dualfprintf(fail_file,"No such concovtype=%d\n",concovtype);
    myexit(1);
  }



  // get SPC -> Cart transformation
  FTYPE lambdacoord[NDIM][NDIM];
  generate_lambdacoord(oldgridtype, newgridtype, V, lambdacoord);


  // transform from spherical polar to Cartesian coordinates (assumes tensor is in orthonormal basis)
  DLOOP(jj,kk) Tortho[jj][kk]=0.0;
  DLOOP(jj,kk){
    DLOOP(ll,mm){
      Tortho[jj][kk] += lambdacoord[jj][ll]*lambdacoord[kk][mm]*Torthospc[ll][mm];
    }
  }


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
    //    if(jj==1 || jj==2){
    if(jj==1){
      Tit[jj]=bsq*ucon[jj]*ucov[0] - bcon[jj]*bcov[0] + delta(jj,0)*(bsq*0.5);
    }
    else Tit[jj]=0.0;
  }

  // |T^p_t|
  Titsq = 0.0;
  DLOOP(jj,kk) Titsq+=Tit[jj]*Tit[kk]*gcov[GIND(jj,kk)];

  // dP/d\omega = \detg T^p_t[EM]/sin(\theta)
  //dPdw=sqrt(fabs(Titsq))*fabs(gdet/sin(V[2]));
  //dPdw=sqrt(fabs(Titsq))*fabs(gdet);
  dPdw=sqrt(fabs(Titsq));

  // final output
  *compout = dPdw;


}

static void vB2poyntingdensityfull(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vecv, FTYPE *vecB, FTYPE (*Tmunu)[NDIM])
{

  FTYPE newgcov[SYMMATRIXNDIM];
  int jj,kk,ll,mm,pp;
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];
  FTYPE pr[NPR];
  FTYPE ucov[NDIM],ucon[NDIM];
  FTYPE bcon[NDIM],bcov[NDIM],bsq;


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

  // T^\mu_\nu[EM] = b^2 u^\mu u_\nu  - b^\mu b_\nu
  DLOOP(jj,kk){
    Tmunu[jj][kk]=bsq*ucon[jj]*ucov[kk] - bcon[jj]*bcov[kk] + delta(jj,kk)*(bsq*0.5);
  }

}


// assume input Cartesian orthonormal quantities
static void vB2poyntingdensityboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*Tmunuboostedcart)[NDIM], FTYPE (*Tmunuboostedspccart)[NDIM], FTYPE *FEMrad3, FTYPE *FEMrad4)
{
  FTYPE Titsq,dPdw;

  FTYPE Tr=Tmunuboostedspccart[1][0]; // T^r_t
  FTYPE Txr=Tmunuboostedspccart[1][1]; // T^r_x

  // now extract specific direction.
  Titsq = Tr*Tr;

  // dP/d\omega = \detg T^p_t[EM]/sin(\theta)
  //dPdw=sqrt(fabs(Titsq))*fabs(gdet/sin(V[2]));
  //dPdw=sqrt(fabs(Titsq))*fabs(gdet);
  dPdw=sqrt(fabs(Titsq));

  // final output
  *FEMrad3 = dPdw;

  FTYPE gamma,betavec[NDIM];
  get_beta_gamma_boost(&gamma, betavec);

  // Sam's version
  *FEMrad4 = Tr - betavec[1]*Txr;
  

}


// assume input Cartesian orthonormal quantities
static void Tboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE (*Tmunu)[NDIM], FTYPE (*Tmunuboostedcart)[NDIM], FTYPE (*Tmunuboostedspccart)[NDIM], FTYPE (*Tmunuboostedspcspc)[NDIM])
{
  int jj,kk,ll,mm,pp;


  FTYPE r=V[1];
  FTYPE th=V[2];
  FTYPE ph=V[3];


  FTYPE lambdatrans[NDIM][NDIM];
  if(0){
    lambdatrans[TT][TT]=1.0;
    SLOOPA(jj) lambdatrans[TT][jj] = lambdatrans[jj][TT] = 0.0;


    // rest come from definitions of {x,y,z}(r,\theta,\phi)
    // assumes orthonormal to orhonormal!
    lambdatrans[0][TT] = 1.0;
    lambdatrans[0][RR] = 0.0;
    lambdatrans[0][TH] = 0.0;
    lambdatrans[0][PH] = 0.0;

    lambdatrans[1][TT] = 0.0;
    lambdatrans[1][RR] = sin(th)*cos(ph);
    lambdatrans[1][TH] = cos(th)*cos(ph);
    lambdatrans[1][PH] = -sin(ph);

    lambdatrans[2][TT] = 0.0;
    lambdatrans[2][RR] = sin(th)*sin(ph);
    lambdatrans[2][TH] = cos(th)*sin(ph);
    lambdatrans[2][PH] = cos(ph);

    lambdatrans[3][TT] = 0.0;
    lambdatrans[3][RR] = cos(th);
    lambdatrans[3][TH] = -sin(th);
    lambdatrans[3][PH] = 0.0;
  }
  else{
    // if inputed ortho, then already Cart as well
  }

  FTYPE gamma,betavec[NDIM];
  get_beta_gamma_boost(&gamma, betavec);
  FTYPE beta=betavec[1]; // for our special case

  FTYPE vlambdatrans[NDIM][NDIM];
  // assume time doesn't change or mix with space

  // Cartesian Lorentz boost in x-direction only
  vlambdatrans[0][0] = gamma;
  vlambdatrans[0][1] = -beta*gamma;
  vlambdatrans[0][2] = 0.0;
  vlambdatrans[0][3] = 0.0;

  vlambdatrans[1][0] = -beta*gamma;
  vlambdatrans[1][1] = gamma;
  vlambdatrans[1][2] = 0.0;
  vlambdatrans[1][3] = 0.0;

  vlambdatrans[2][0] = 0.0;
  vlambdatrans[2][1] = 0.0;
  vlambdatrans[2][2] = 1.0;
  vlambdatrans[2][3] = 0.0;

  vlambdatrans[3][0] = 0.0;
  vlambdatrans[3][1] = 0.0;
  vlambdatrans[3][2] = 0.0;
  vlambdatrans[3][3] = 1.0;


  FTYPE boostA[NDIM][NDIM];
  if(0){
    FTYPE ilambdatrans[NDIM][NDIM];
    matrix_inverse_4d(lambdatrans,ilambdatrans); // inverse and transpose.
    
    // transform Cart boost to SPC boost
    DLOOP(jj,kk) boostA[jj][kk]=0.0;
    DLOOP(kk,mm){
      DLOOP(jj,ll){
        // B^kk_mm = L^kk_jj (iL)^ll_mm V^jj_ll
        boostA[kk][mm] += lambdatrans[kk][jj]*ilambdatrans[ll][mm]*vlambdatrans[jj][ll];
      }
    }
  }
  else{
    DLOOP(jj,kk) boostA[jj][kk]=vlambdatrans[jj][kk];
  }


  FTYPE iboostA[NDIM][NDIM];
  matrix_inverse_4d(boostA,iboostA); // inverse and transpose.

  // transform Tmunu in orthonormal 
  //  FTYPE Tmunuboostedcart[NDIM][NDIM];
  DLOOP(jj,kk) Tmunuboostedcart[jj][kk]=0.0;
  DLOOP(kk,mm){
    DLOOP(jj,ll){
      // T^kk_mm = B^kk_jj (iB)^ll_mm   T^jj_ll
      Tmunuboostedcart[kk][mm] += boostA[kk][jj]*iboostA[ll][mm]*Tmunu[jj][ll];
    }
  }

  // get ortho SPC version: T^i_t
  FTYPE Ttt=Tmunuboostedcart[0][0]; // T^t_t
  FTYPE Txt=Tmunuboostedcart[1][0]; // T^x_t
  FTYPE Tyt=Tmunuboostedcart[2][0]; // T^y_t
  FTYPE Tzt=Tmunuboostedcart[3][0]; // T^x_t
  // get SPC
  FTYPE Trt=sin(th)*cos(ph)*Txt + sin(th)*sin(ph)*Tyt + cos(th)*Tzt; // T^r_t
  FTYPE Tht=cos(th)*cos(ph)*Txt + cos(th)*sin(ph)*Tyt - sin(th)*Tzt; // T^h_t
  FTYPE Tpt=-sin(ph)*Txt + cos(ph)*Tyt; // T^p_t

  // get ortho SPC version: T^i_x
  FTYPE Ttx=Tmunuboostedcart[0][1]; // T^t_x
  FTYPE Txx=Tmunuboostedcart[1][1]; // T^x_x
  FTYPE Tyx=Tmunuboostedcart[2][1]; // T^y_x
  FTYPE Tzx=Tmunuboostedcart[3][1]; // T^z_x
  // get SPC
  FTYPE Trx=sin(th)*cos(ph)*Txx + sin(th)*sin(ph)*Tyx + cos(th)*Tzx; // T^r_x
  FTYPE Thx=cos(th)*cos(ph)*Txx + cos(th)*sin(ph)*Tyx - sin(th)*Tzx; // T^h_x
  FTYPE Tpx=-sin(ph)*Txx + cos(ph)*Tyx; // T^p_x
  
  // get ortho SPC version: T^i_y
  FTYPE Tty=Tmunuboostedcart[0][2]; // T^t_y
  FTYPE Txy=Tmunuboostedcart[1][2]; // T^x_y
  FTYPE Tyy=Tmunuboostedcart[2][2]; // T^y_y
  FTYPE Tzy=Tmunuboostedcart[3][2]; // T^z_y
  // get SPC
  FTYPE Try=sin(th)*cos(ph)*Txy + sin(th)*sin(ph)*Tyy + cos(th)*Tzy; // T^r_y
  FTYPE Thy=cos(th)*cos(ph)*Txy + cos(th)*sin(ph)*Tyy - sin(th)*Tzy; // T^h_y
  FTYPE Tpy=-sin(ph)*Txy + cos(ph)*Tyy; // T^p_y

  // get ortho SPC version: T^i_z
  FTYPE Ttz=Tmunuboostedcart[0][3]; // T^t_z
  FTYPE Txz=Tmunuboostedcart[1][3]; // T^x_z
  FTYPE Tyz=Tmunuboostedcart[2][3]; // T^y_z
  FTYPE Tzz=Tmunuboostedcart[3][3]; // T^z_z
  // get SPC
  FTYPE Trz=sin(th)*cos(ph)*Txz + sin(th)*sin(ph)*Tyz + cos(th)*Tzz; // T^r_z
  FTYPE Thz=cos(th)*cos(ph)*Txz + cos(th)*sin(ph)*Tyz - sin(th)*Tzz; // T^h_z
  FTYPE Tpz=-sin(ph)*Txz + cos(ph)*Tyz; // T^p_z


  // GET full SPC
  FTYPE Ttr=sin(th)*cos(ph)*Ttx + sin(th)*sin(ph)*Tty + cos(th)*Ttz; // T^t_r
  FTYPE Tth=cos(th)*cos(ph)*Ttx + cos(th)*sin(ph)*Tty - sin(th)*Ttz; // T^t_h
  FTYPE Ttp=-sin(ph)*Ttx + cos(ph)*Tty; // T^t_p

  FTYPE Trr=sin(th)*cos(ph)*Trx + sin(th)*sin(ph)*Try + cos(th)*Trz; // T^r_r
  FTYPE Trh=cos(th)*cos(ph)*Trx + cos(th)*sin(ph)*Try - sin(th)*Trz; // T^r_h
  FTYPE Trp=-sin(ph)*Trx + cos(ph)*Try; // T^r_p

  FTYPE Thr=sin(th)*cos(ph)*Thx + sin(th)*sin(ph)*Thy + cos(th)*Thz; // T^h_r
  FTYPE Thh=cos(th)*cos(ph)*Thx + cos(th)*sin(ph)*Thy - sin(th)*Thz; // T^h_h
  FTYPE Thp=-sin(ph)*Thx + cos(ph)*Thy; // T^h_p

  FTYPE Tpr=sin(th)*cos(ph)*Tpx + sin(th)*sin(ph)*Tpy + cos(th)*Tpz; // T^p_r
  FTYPE Tph=cos(th)*cos(ph)*Tpx + cos(th)*sin(ph)*Tpy - sin(th)*Tpz; // T^p_h
  FTYPE Tpp=-sin(ph)*Tpx + cos(ph)*Tpy; // T^p_p


  // T^{trhp}_t
  Tmunuboostedspccart[0][0]=Ttt;
  Tmunuboostedspccart[1][0]=Trt;
  Tmunuboostedspccart[2][0]=Tht;
  Tmunuboostedspccart[3][0]=Tpt;

  // T^{trhp}_x
  Tmunuboostedspccart[0][1]=Ttx;
  Tmunuboostedspccart[1][1]=Trx;
  Tmunuboostedspccart[2][1]=Thx;
  Tmunuboostedspccart[3][1]=Tpx;

  // T^{trhp}_y
  Tmunuboostedspccart[0][2]=Tty;
  Tmunuboostedspccart[1][2]=Try;
  Tmunuboostedspccart[2][2]=Thy;
  Tmunuboostedspccart[3][2]=Tpy;

  // T^{trhp}_z
  Tmunuboostedspccart[0][3]=Ttz;
  Tmunuboostedspccart[1][3]=Trz;
  Tmunuboostedspccart[2][3]=Thz;
  Tmunuboostedspccart[3][3]=Tpz;
  

  // T^{trhp}_t
  Tmunuboostedspcspc[0][0]=Ttt;
  Tmunuboostedspcspc[1][0]=Trt;
  Tmunuboostedspcspc[2][0]=Tht;
  Tmunuboostedspcspc[3][0]=Tpt;

  // T^{trhp}_r
  Tmunuboostedspcspc[0][1]=Ttr;
  Tmunuboostedspcspc[1][1]=Trr;
  Tmunuboostedspcspc[2][1]=Thr;
  Tmunuboostedspcspc[3][1]=Tpr;

  // T^{trhp}_h
  Tmunuboostedspcspc[0][2]=Tth;
  Tmunuboostedspcspc[1][2]=Trh;
  Tmunuboostedspcspc[2][2]=Thh;
  Tmunuboostedspcspc[3][2]=Tph;

  // T^{trhp}_p
  Tmunuboostedspcspc[0][3]=Ttp;
  Tmunuboostedspcspc[1][3]=Trp;
  Tmunuboostedspcspc[2][3]=Thp;
  Tmunuboostedspcspc[3][3]=Tpp;
  

}







// input PRIMECOORD coordinate basis v^i and output MCOORD coordinate basis v_{vectorcomponent}
// so no transformation to orthonormal basis or even to Cartesian coordaintes
static void vecup2vecdown(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecdown)
{
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE newgcov[SYMMATRIXNDIM];
  int jj,kk,ll,pp;
  FTYPE tempcomp[NDIM];


  // get inverse coordinate transformation
  idxdxprim(dxdxp, idxdxp); // from coord.c

  // get MCOORD metric from PRIMECOORD metric
  DLOOP(jj,kk){
    newgcov[GIND(jj,kk)]=0.0;
    DLOOP(ll,pp) {
      newgcov[GIND(jj,kk)] += GINDASSIGNFACTOR(jj,kk)*gcov[GIND(ll,pp)]*idxdxp[ll][jj]*idxdxp[pp][kk];
    }
  }


  // compute lower components
  DLOOPA(jj) vecdown[jj]=0.0;
  DLOOP(jj,kk) vecdown[jj] += vec[kk]*newgcov[GIND(jj,kk)];

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
void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], struct of_geom *ptrgeom)
{
  int jj,kk,ll,pp;

  // see HARM's dump.c and SM's gammie.m grid3d
  if(0){
    SLOOPA(jj) fscanf(in,"%d",&ti[jj]);
    SLOOPA(jj) fscanf(in,SCANARG,&X[jj]);
    SLOOPA(jj) fscanf(in,SCANARG,&V[jj]);
    DLOOPA(jj) DLOOPA(kk) DLOOPA(ll) fscanf(in,SCANARG,&conn[jj][kk][ll]);
    DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&gcon[GIND(jj,kk)]); // notice that read-in 16 elements but only store 10 unique ones
    DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&gcov[GIND(jj,kk)]); // notice that read-in 16 elements but only store 10 unique ones
    fscanf(in,SCANARG,gdet); // gdet already pointer
    DLOOPA(jj) fscanf(in,SCANARG,&ck[jj]);
    DLOOPA(jj) DLOOPA(kk) fscanf(in,SCANARG,&dxdxp[jj][kk]);
  }
  else{
    FTYPE tiFTYPE[NDIM];
    // 3
    SLOOPA(jj){ readelement(binaryinputgdump,inFTYPEgdump,in,&tiFTYPE[jj]); ti[jj]=(int)tiFTYPE[jj];} // file stores gdump columns as FTYPE when binary format, and ok to read as float if text format.
    
    // 6
    SLOOPA(jj) readelement(binaryinputgdump,inFTYPEgdump,in,&X[jj]);
    // 9
    SLOOPA(jj) readelement(binaryinputgdump,inFTYPEgdump,in,&V[jj]);
    // 9+64=73
    DLOOPA(jj) DLOOPA(kk) DLOOPA(ll) readelement(binaryinputgdump,inFTYPEgdump,in,&conn[jj][kk][ll]);
    // 73+16=89
    DLOOPA(jj) DLOOPA(kk) readelement(binaryinputgdump,inFTYPEgdump,in,&gcon[GIND(jj,kk)]); // notice that read-in 16 elements but only store 10 unique ones
    // 89+16=105
    DLOOPA(jj) DLOOPA(kk) readelement(binaryinputgdump,inFTYPEgdump,in,&gcov[GIND(jj,kk)]); // notice that read-in 16 elements but only store 10 unique ones
    // 106
    readelement(binaryinputgdump,inFTYPEgdump,in,gdet); // gdet already pointer
    // 110
    DLOOPA(jj) readelement(binaryinputgdump,inFTYPEgdump,in,&ck[jj]);
    // 126
    DLOOPA(jj) DLOOPA(kk) readelement(binaryinputgdump,inFTYPEgdump,in,&dxdxp[jj][kk]);


    // DEBUG:
    //    dualfprintf(fail_file,"readgdump: %d %d %d : dxdxp00=%g\n",ti[1],ti[2],ti[3],dxdxp[0][0]);


  }



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








void raise_vec(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon)
{

  ucon[0] = geom->gcon[GIND(0,0)]*ucov[0]
    + geom->gcon[GIND(0,1)]*ucov[1]
    + geom->gcon[GIND(0,2)]*ucov[2]
    + geom->gcon[GIND(0,3)]*ucov[3] ;
  ucon[1] = geom->gcon[GIND(0,1)]*ucov[0]
    + geom->gcon[GIND(1,1)]*ucov[1]
    + geom->gcon[GIND(1,2)]*ucov[2]
    + geom->gcon[GIND(1,3)]*ucov[3] ;
  ucon[2] = geom->gcon[GIND(0,2)]*ucov[0]
    + geom->gcon[GIND(1,2)]*ucov[1]
    + geom->gcon[GIND(2,2)]*ucov[2]
    + geom->gcon[GIND(2,3)]*ucov[3] ;
  ucon[3] = geom->gcon[GIND(0,3)]*ucov[0]
    + geom->gcon[GIND(1,3)]*ucov[1]
    + geom->gcon[GIND(2,3)]*ucov[2]
    + geom->gcon[GIND(3,3)]*ucov[3] ;

  return ;
}


void lower_vec(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov)
{
  ucov[0] = geom->gcov[GIND(0,0)]*ucon[0]
    + geom->gcov[GIND(0,1)]*ucon[1]
    + geom->gcov[GIND(0,2)]*ucon[2]
    + geom->gcov[GIND(0,3)]*ucon[3] ;
  ucov[1] = geom->gcov[GIND(0,1)]*ucon[0]
    + geom->gcov[GIND(1,1)]*ucon[1]
    + geom->gcov[GIND(1,2)]*ucon[2]
    + geom->gcov[GIND(1,3)]*ucon[3]
    ;
  ucov[2] = geom->gcov[GIND(0,2)]*ucon[0]
    + geom->gcov[GIND(1,2)]*ucon[1]
    + geom->gcov[GIND(2,2)]*ucon[2]
#if(DOMIXTHETAPHI)
    + geom->gcov[GIND(2,3)]*ucon[3]
#endif
    ;
  ucov[3] = geom->gcov[GIND(0,3)]*ucon[0]
    + geom->gcov[GIND(1,3)]*ucon[1]
#if(DOMIXTHETAPHI)
    + geom->gcov[GIND(2,3)]*ucon[2]
#endif
    + geom->gcov[GIND(3,3)]*ucon[3] ;

  return ;
}


// get boosted velocity (not Lorentz correct, just correct in non-rel limit)
static void vboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *vecvboosted)
{
  int jj,kk;

  FTYPE r=V[1];
  FTYPE th=V[2];
  FTYPE ph=V[3];

  int BOOSTFIELD=1; // for moving BH problem


  if(BOOSTFIELD){
 
    
    FTYPE vecvboost[NDIM]={0.0};
    getvboost(ti,  X,  V,  conn,  gcon,  gcov,  gdet,  ck,  dxdxp, oldgridtype, newgridtype, vecvboost);


    // add non-ortho coordinate basis velocity to ortho
    FTYPE finalvec[NDIM];
    finalvec[TT]=0.0;
    if(r<3.2){
      //if(0){
      finalvec[RR]=vecv[1];
      finalvec[TH]=vecv[2];
      finalvec[PH]=vecv[3];
    }
    else{
      // no TT change since non-rel
      finalvec[RR]=vecv[1] + vecvboost[1]/sqrt(gcov[GIND(1,1)]);
      finalvec[TH]=vecv[2] + vecvboost[2]/sqrt(gcov[GIND(2,2)]);
      finalvec[PH]=vecv[3] + vecvboost[3]/sqrt(gcov[GIND(3,3)]);
    }

    vecvboosted[0] = vecv[0]; // non-rel
    vecvboosted[1] = finalvec[RR];
    vecvboosted[2] = finalvec[TH];
    vecvboosted[3] = finalvec[PH];
  }    
}




// gives vecvboosted: Cartesian frame specific value as written in SPC (no boosting, just SPC value)
static void getvboost(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vecvboost)
{
  int jj,kk;

  FTYPE r=V[1];
  FTYPE th=V[2];
  FTYPE ph=V[3];

  int BOOSTFIELD=1; // for moving BH problem


  if(BOOSTFIELD){
    // BOOST of field
    FTYPE xx=r*sin(th)*cos(ph);
    FTYPE yy=r*sin(th)*sin(ph);
    FTYPE zz=r*cos(th);
    FTYPE lambdatrans[NDIM][NDIM];
    FTYPE ilambdatrans[NDIM][NDIM];
    // assume time doesn't change or mix with space
    lambdatrans[TT][TT]=1.0;
    SLOOPA(jj) lambdatrans[TT][jj] = lambdatrans[jj][TT] = 0.0;

    ilambdatrans[TT][TT]=1.0;
    SLOOPA(jj) ilambdatrans[TT][jj] = ilambdatrans[jj][TT] = 0.0;

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

    // Cart 2 SPC
    ilambdatrans[1][RR] = sin(th)*cos(ph);
    ilambdatrans[1][TH] = cos(th)*cos(ph);
    ilambdatrans[1][PH] = -sin(ph);

    ilambdatrans[2][RR] = sin(th)*sin(ph);
    ilambdatrans[2][TH] = cos(th)*sin(ph);
    ilambdatrans[2][PH] = cos(ph);

    ilambdatrans[3][RR] = cos(th);
    ilambdatrans[3][TH] = -sin(th);
    ilambdatrans[3][PH] = 0.0;

    // quasi-orthonormal
    FTYPE gamma,betavec[NDIM];
    get_beta_gamma_boost(&gamma, betavec);

    FTYPE finalvec[NDIM];
    finalvec[TT]=gamma;
    finalvec[RR]=-gamma*betavec[1]; // x // OPPOSITE of what's in init.c
    finalvec[TH]=-gamma*betavec[2]; // y
    finalvec[PH]=-gamma*betavec[3]; // z


    // transform from ortho Cart to ortho SPC
    FTYPE tempcomp[NDIM];
    DLOOPA(jj) tempcomp[jj]=0.0;
    DLOOP(jj,kk){
      tempcomp[kk] += ilambdatrans[jj][kk]*finalvec[jj];
    }
    DLOOPA(jj) finalvec[jj]=tempcomp[jj]; // spc


    vecvboost[0] = finalvec[TT];
    vecvboost[1] = finalvec[RR];
    vecvboost[2] = finalvec[TH];
    vecvboost[3] = finalvec[PH];
  }    
}




/// dual of Maxwell tensor
/// returns \dF^{\mu \nu}
static int Mcon_calc_uB(FTYPE *ucon, FTYPE *Bcon, FTYPE (*Mcon)[NDIM])
{
  int j,k;
  FTYPE vcon[NDIM];


  // diagonal is 0
  DLOOPA(j) Mcon[j][j]=0.0;

  // space-time terms
  SLOOPA(k) {
    // \dF^{it} = B^i = Bcon[i]
    Mcon[k][0] = Bcon[k] ; 
    Mcon[0][k] = - Mcon[k][0] ;
  }

  // get v^i
  SLOOPA(k) vcon[k]= ucon[k]/ucon[TT];

  // space-space terms
  //  SLOOP(j,k) Mcon[j][k] = (Bcon[1+j-1] * vcon[k] - Bcon[1+k-1] * vcon[j]);
  // optimize
  Mcon[1][2] = (Bcon[1] * vcon[2] - Bcon[2] * vcon[1]);
  Mcon[1][3] = (Bcon[1] * vcon[3] - Bcon[3] * vcon[1]);
  Mcon[2][3] = (Bcon[2] * vcon[3] - Bcon[3] * vcon[2]);
  Mcon[2][1] = -Mcon[1][2];
  Mcon[3][1] = -Mcon[1][3];
  Mcon[3][2] = -Mcon[2][3];
  Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0.0;


  return(0);

}



/// dual of Maxwell tensor
/// returns \dF^{\mu \nu}
static int Mcon_calc_ub(FTYPE *ucon, FTYPE *bcon, FTYPE (*Mcon)[NDIM])
{
  int j,k;
  FTYPE vcon[NDIM];

  DLOOP(j,k) Mcon[j][k] = bcon[j] * ucon[k] - bcon[k] * ucon[j];


  return(0);

}




/// Maxwell to Faraday
/// which=0 : Mcon -> Fcov (for clean Mcon, Fcov has \detg)
/// which=1 : Mcov -> Fcon (for clean Mcov)
/// which=2 : Fcon -> Mcov
/// which=3 : Fcov -> Mcon
/// copies faraday_calc() in phys.c
static void MtoF_simple(int which, FTYPE (*invar)[NDIM], FTYPE gdet, FTYPE (*outvar)[NDIM])
{
  int nu,mu,kappa,lambda;
  FTYPE prefactor;
  int whichlc;

  if((which==0)||(which==1)){
    prefactor=-0.5;
    if(which==0) whichlc=0;
    if(which==1) whichlc=1;
  }
  if((which==2)||(which==3)){
    prefactor=0.5;
    if(which==2) whichlc=0;
    if(which==3) whichlc=1;
  }

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
      outvar[mu][nu]=0.0;
      for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
          outvar[mu][nu]+=prefactor*lc4(whichlc,gdet,mu,nu,kappa,lambda)*invar[kappa][lambda];
        }
    }


}


static void get_beta_gamma_boost(FTYPE *gamma, FTYPE *betavec)
{

  FTYPE beta=+0.3; // what's in init.c
  *gamma = 1.0/sqrt(1.0-beta*beta);

  betavec[1]=beta; // in x direction (note Loretnz boost in code is purely x-directional)
  betavec[2]=0.0;
  betavec[3]=0.0;
}

