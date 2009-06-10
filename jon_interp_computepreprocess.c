#include "decs.h"



// local declarations of local functions
static void vec2vecortho(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho);
static void vB2poyntingdensity(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vecv, FTYPE *vecB, FTYPE *compout);
static void vecup2vecdowncomponent(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, int vectorcomponent, FTYPE *vec, FTYPE *compout);
static void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM]);
static void generate_lambdacoord(int oldgridtype, int newgridtype, FTYPE *V, FTYPE (*lambdacoord)[NDIM]);
static void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon);







// process inputted data
void compute_preprocess(FILE *gdumpfile, FTYPE *finaloutput)
{
  FTYPE vec[NDIM],vecv[NDIM],vecB[NDIM];
  FTYPE vecortho[NDIM];
  int jj;
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
  //

  // no need to recompute X,V if reading in gdump
  // assumes input data at CENT
  //  coord(i,j,k,CENT,X);
  // use bl_coord()
  //bl_coord(X,V);

  // perform local transformation of tensor objects prior to spatial interpolation
  if(outputvartype==1){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp);

    for(iter=0;iter<num4vectors;iter++){
      // convert coordinate basis vector compnents to single orthonormal basis component desired
      fscanf(stdin,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
      // instantly transform vector from original to new coordinate system while reading in to avoid excessive memory use
      vec2vecortho(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vec, vecortho);
      if(immediateoutput){
	// immediately output result
	if(sizeof(FTYPE)==sizeof(double)){
	  DLOOPA(jj) fprintf(stdout,"%22.16g ",vecortho[jj]);
	}
	else if(sizeof(FTYPE)==sizeof(float)){
	  DLOOPA(jj) fprintf(stdout,"%15.7g ",vecortho[jj]);
	}
      }
      else{
	finaloutput[iter]=vecortho[vectorcomponent];
      }
    }
    if(immediateoutput){
      // output return after entire row is done
      fprintf(stdout,"\n") ;
    }
  }
  else if(outputvartype==2){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp);

    for(iter=0;iter<num4vectors;iter++){
      // input uu0 vu1 vu2 vu3
      fscanf(stdin,SCANARG4VEC,&vecv[0],&vecv[1],&vecv[2],&vecv[3]) ; // vu^i=uu^i/uu0 (i.e. not uu^i as maybe expected)
      SLOOPA(jj) vecv[jj]*=vecv[0]; // now uu[jj]
      // input B^1 B^2 B^3
      vecB[0]=0.0;
      fscanf(stdin,SCANARGVEC,&vecB[1],&vecB[2],&vecB[3]) ;
      // compute
      vB2poyntingdensity(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponent, vecv, vecB, &finaloutput[iter]);
    }
  }
  else if(outputvartype==3){

    // first get gdump data (only once per call to compute_preprocess() !!)
    read_gdumpline(gdumpfile, ti,  X,  V,  conn,  gcon,  gcov,  &gdet,  ck,  dxdxp);

    for(iter=0;iter<num4vectors;iter++){
      // input
      fscanf(stdin,SCANARG4VEC,&vec[0],&vec[1],&vec[2],&vec[3]) ;
      // compute
      vecup2vecdowncomponent(ti,X,V,conn,gcon,gcov,gdet,ck,dxdxp,oldgridtype, newgridtype, vectorcomponent, vec, &finaloutput[iter]);
    }
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
static void vec2vecortho(int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM], int oldgridtype, int newgridtype, FTYPE *vec, FTYPE *vecortho)
{
  FTYPE lambdacoord[NDIM][NDIM];
  int jj,kk;
  FTYPE tetrcov[NDIM][NDIM],tetrcon[NDIM][NDIM],eigenvalues[NDIM];
  FTYPE tempcomp[NDIM];
  FTYPE finalvec[NDIM];


  DLOOPA(jj) finalvec[jj]=vec[jj];


  // get SPC -> Cart transformation
  generate_lambdacoord(oldgridtype, newgridtype, V, lambdacoord);

  // get tetrad
  tetr_func_frommetric(dxdxp, gcov, tetrcov, tetrcon, eigenvalues);


  // transform to orthonormal basis
  DLOOPA(jj) tempcomp[jj]=0.0;
  DLOOP(jj,kk){
    tempcomp[kk] += tetrcon[kk][jj]*finalvec[jj];
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];


  // transform from spherical polar to Cartesian coordinates (assumes vector is in orthonormal basis)
  DLOOPA(jj) tempcomp[jj]=0.0;
  DLOOP(jj,kk){
    tempcomp[kk] += lambdacoord[kk][jj]*finalvec[jj];
  }
  DLOOPA(jj) finalvec[jj]=tempcomp[jj];
  

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




// input coordinate basis v^i and output coordinate basis v_{vectorcomponent}
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





static void generate_lambdacoord(int oldgridtype, int newgridtype, FTYPE *V, FTYPE (*lambdacoord)[NDIM])
{
  FTYPE r,th,ph;
  int jj,kk;


  if(oldgridtype==newgridtype || newgridtype==GRIDTYPENOCHANGE){
    // no change
    DLOOP(jj,kk) lambdacoord[jj][kk]=0.0;
    DLOOPA(jj) lambdacoord[jj][jj]=1.0;
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
    lambdacoord[TT][TT]=1.0;
    SLOOPA(jj) lambdacoord[TT][jj] = lambdacoord[jj][TT] = 0.0;

    // rest come from definitions of {x,y,z}(r,\theta,\phi)
    lambdacoord[1][RR] = sin(th)*cos(ph);
    lambdacoord[1][TH] = r*cos(th)*cos(ph);
    lambdacoord[1][PH] = -r*sin(th)*sin(ph);

    lambdacoord[2][RR] = sin(th)*sin(ph);
    lambdacoord[2][TH] = r*cos(th)*sin(ph);
    lambdacoord[2][PH] = r*sin(th)*cos(ph);

    lambdacoord[3][RR] = cos(th);
    lambdacoord[3][TH] = -r*sin(th);
    lambdacoord[3][PH] = 0.0;
  }
  else{
    dualfprintf(fail_file,"No transformation setup for oldgridtype=%d newgridtype=%d\n",oldgridtype,newgridtype);
    myexit(246347);
  }
}





// read a line from the gdump file
// assumed to be used in same order as data and gdump file as looped over in this code
// and of course only can be done once per data point
// fastest index is most-right element
static void read_gdumpline(FILE *in, int ti[],  FTYPE X[],  FTYPE V[],  FTYPE (*conn)[NDIM][NDIM],  FTYPE *gcon,  FTYPE *gcov,  FTYPE *gdet,  FTYPE ck[],  FTYPE (*dxdxp)[NDIM])
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
