
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define EXTRAP 1
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  SFTYPE Dt=-1.0; // no boundary zone should ever be in accounting
  int i, j, k;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j
  FTYPE prescale[NPR];
  int bound_disk(FTYPE prim[][N2M][NPR]);


  if((boundstage==STAGE0)||(boundstage==STAGEM1)){





    if (mycpupos[1] == 0) {
      if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
	/* inner r boundary condition: u, just copy */
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--)	  PBOUNDLOOP{
	    //prim[i][j][k] = exp((log(prim[ri+1][rj][k])-log(prim[ri][rj][k]))*(i-ri)+log(prim[ri][rj][k]));
	    //prim[i][j][k] = (prim[ri+1][rj][k]-prim[ri][rj][k])*(i-ri)+prim[ri][rj][k];
	    prim[i][j][k] = prim[ri][rj][k];
	  }
#elif(EXTRAP==1)

	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    k=U1;	    // treat U1 as special
	    prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    for(k=U2;k<=U3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. + (i-ri)*dx[1]) ;
	    }
	    k=B1; // treat B1 special
	    prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    for(k=B2;k<=B3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. + (i-ri)*dx[1]) ;
	    }
	  }
#elif(EXTRAP==2)
	  ri=0;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess
	    PBOUNDLOOP prim[i][j][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);	    
	    rescale(-1,1,prim[i][j],&geom,prescale);
	  }
#endif
	  for(i=-1;i>=-N1BND;i--){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom, Dt) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom, Dt) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#endif
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom,Dt) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom,Dt)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	  }
	}
      }
      else if(BCtype[X1DN]==FIXED){
	// inner radial fixed boundary conditions
	for (j = -N2BND; j < N2+N2BND; j++)
	  for(i=-1;i>=-N1BND;i--)
	    PBOUNDLOOP prim[i][j][k] = panalytic[i][j][k];
	//PBOUNDLOOP prim[i][j][k] = p[i][j][k];
      }

    }// end if mycpupos==0







    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
	/* outer r BC: outflow */
      
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++)	  PBOUNDLOOP{
	    prim[i][j][k] = (prim[ri][rj][k]-prim[ri-1][rj][k])*(i-ri)+prim[ri][rj][k];
	  }
#elif(EXTRAP==1)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    k=U1; // treat U1 as special
	    prim[i][j][k] = prim[ri][rj][k] * (1. - 2*(i-ri)*dx[1]) ;
	    for(k=U2;k<=U3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    }
	    k=B1; // treat B1 special
	    prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    for(k=B2;k<=B3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    }
	  }
#elif(EXTRAP==2)
	  ri=N1-1;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  for(i=N1;i<=N1-1+N1BND;i++){
	    // set guess
	    PBOUNDLOOP prim[i][j][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    rescale(-1,1,prim[i][j],&geom,prescale);
	  }
#endif

	  for(i=N1;i<=N1-1+N1BND;i++){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom,Dt) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom,Dt) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#endif
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom,Dt) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom, Dt)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
	  }
	}
      }
      else if(BCtype[X1UP]==FIXED){
	// outer radial fixed boundary conditions
	for (j = -N2BND; j < N2+N2BND; j++)
	  for(i=N1;i<=N1-1+N1BND;i++)
	    PBOUNDLOOP prim[i][j][k] = panalytic[i][j][k];
	//PBOUNDLOOP prim[i][j][k] = p[i][j][k];
      }
    }


    // to help protect the pole from death blows to the computational grid
    // a sort of crushing regularization
#define POLEDEATH 1 // choice
    // causes problems with stability at just beyond pole
    // for field line plots, can just set B^\theta=0 along pole



    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      if(BCtype[X2DN]==POLARAXIS){
	for (i = -N1BND; i <=N1-1+N1BND; i++){
	  ri=i;
	  rj=0;
	  for(j=-N2BND;j<=-1;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j-1)][k];
	}

	/* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	/* inner pole */
	for (i = -N1BND; i <= N1-1+N1BND; i++) {
	  for (j = -N2BND; j < 0; j++) {
	    if(POSDEFMETRIC==0){
	      // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	      // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	      prim[i][j][U2] *= -1.;
	      prim[i][j][U3] *= 1.;
	      prim[i][j][B2] *= -1.;
	      prim[i][j][B3] *= 1.;
	    }
	    else{
	      prim[i][j][U2] *= -1.;
	      prim[i][j][U3] *= 1.;
	      prim[i][j][B2] *= -1.;
	      prim[i][j][B3] *= 1.;
	    }
	  }
	}


	for (i = -N1BND; i <= N1-1+N1BND; i++) {
	  for (j = -N2BND; j < 0+POLEDEATH; j++) {
	    if(POSDEFMETRIC==0){
	      // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	      // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	      prim[i][j][U2] *= 0;
	      prim[i][j][B2] *= 0.;
	    }
	    else{
	      prim[i][j][U2] *= 0.;
	      prim[i][j][B2] *= 0.;
	    }
	  }
	}

      }
      else if(BCtype[X2DN]==FIXED){
	for (i = -N1BND; i <=N1-1+N1BND; i++){
	  for(j=-N2BND;j<=-1;j++)  PBOUNDLOOP prim[i][j][k] = panalytic[i][j][k];
	}
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
      if(BCtype[X2UP]==POLARAXIS){
	for (i = -N1BND; i <=N1-1+N1BND; i++){
	  ri=i;
	  rj=N2-1;
	  for(j=N2;j<=N2-1+N2BND;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j+1)][k];
	}

	for (i = -N1BND; i <= N1-1+N1BND; i++) {
	  for (j = N2; j <= N2-1+N2BND; j++) {
	    // normal axis
#if(1)
	    prim[i][j][U1] *= 1.;
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B1] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
#endif
	    // monopole
#if(0)
	    prim[i][j][U1] *= 1.;
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B1] *= -1.;
	    prim[i][j][B2] *= 1.;
	    prim[i][j][B3] *= -1.;
#endif
#if(0)
	    // split monopole
	    prim[i][j][U1] *= 1.;
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B1] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
#endif
	  }
	}

	for (i = -N1BND; i <= N1-1+N1BND; i++) {
	  for (j = N2-POLEDEATH; j <= N2-1+N2BND; j++) {
	    if(POSDEFMETRIC==0){
	      prim[i][j][U2] *= 0.;
	      prim[i][j][B2] *= 0.;
	    }
	    else{
	      prim[i][j][U2] *= 0.;
	      prim[i][j][B2] *= 0.;
	    }
	  }
	}


      }
      else if(BCtype[X2UP]==FIXED){
	for (i = -N1BND; i <=N1-1+N1BND; i++){
	  for(j=N2;j<=N2-1+N2BND;j++) PBOUNDLOOP prim[i][j][k] = panalytic[i][j][k];
	}
      }
    }


  

  }// end if stage0 or stagem1



  if (USEMPI) bound_mpi(boundstage, prim);



  // real boundary zones, but analytic.  Put after MPI so can use data.  Assumes analytic sets ghost zones as well analytically
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    // DISK BOUNDARY CONDITION, modifies evolution region
    if(1){
      bound_disk(prim);
    }
  }


  return (0);
}






int bound_disk(FTYPE prim[][N2M][NPR])
{
  int i, j, k;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j
  FTYPE prescale[NPR];
  FTYPE primtest[NPR];

  FTYPE Bcon[NDIM],Bcov[NDIM],ucon[NDIM],Bsq;
  int sanitycheck;
  FTYPE X[NDIM];
  FTYPE r,th;
  extern FTYPE get_omegastar(struct of_geom *geom, FTYPE r, FTYPE th);
  FTYPE Mcon[NDIM][NDIM],Fcov[NDIM][NDIM],romegaf2;
  extern void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  extern void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);
  FTYPE realgammamax;
  int tscale;
  FTYPE fun;
  static int firsttime=1;
  int limit_3vel(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr);


#define DELTAJ 2

  //#define RDISKINNER (2.0*Rhor)
#define RDISKINNER (0.0*Rhor)
#define RDISKOUTER (1E30)

#define RNEARBLACKHOLE (0.0) // outflow \Omega_F for r<RNEARBLACKHOLE (doesn't work apparently)

  // time while static boundary conditions
#define TSTATIC (15.0)

  // whether to use new method to limit v^i to be time-like
#define LIMIT3VELBOUND 1

  // check boundary situation (can't copy from nothing)
  if( ((startpos[2]==totalsize[2]/2)||(startpos[2]==totalsize[2]/2-1)) && ( (DELTAJ>=N1BND)||(DELTAJ>=N2BND)) ){
    dualfprintf(fail_file,"DELTAJ out of bounds\n");
    myexit(1);
  }

#if(1)
  if(!firsttime){
    // make the underlying evolved state exactly symmetric
    LOOPF2 LOOPF1{ // loop over entire domain looking for equator
      if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<0)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ) ){ // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
	// one side of equator

	coord(i, j, CENT, X);
	bl_coord(X, &r, &th);

	if((r>RDISKINNER)&&(r<RDISKOUTER)){

	  ri=i;
	  // symmetric partner across equator (for even or odd sized grids)
	  // global j
	  rj=totalsize[2]/2-(startpos[2]+j) + (totalsize[2]-1)/2;
	  // local j
	  rj=rj - startpos[2];

	  // average rather than choose one
	  // GODMARK: really doesn't apply at t=0

	  PLOOP{
	    if((k==B1)||(k==B3)||(k==U2)){
	      // antisymmetric
	      prim[i][j][k]=0.5*(prim[i][j][k]-prim[ri][rj][k]);
	      prim[ri][rj][k]=0.5*(prim[ri][rj][k]-prim[i][j][k]);
	    }
	    else{
	      // symmetric
	      prim[i][j][k]=0.5*(prim[i][j][k]+prim[ri][rj][k]);
	      prim[ri][rj][k]=0.5*(prim[i][j][k]+prim[ri][rj][k]);
	    }
	  }
	  //	dualfprintf(fail_file,"i=%d j=%d rj=%d\n",i,j,rj);
	}
      }
    }
  }
#endif





  LOOPF2 LOOPF1{ // loop over entire domain looking for equator

    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
    get_geometry(i, j, CENT, &geom);



    // one side of equator and outside horizon
    if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<0)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ) ){ // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
      // one side of equator

      if((r>RDISKINNER)&&(r<RDISKOUTER)){
	// OUTFLOW some quantities into conductor that is the current sheet
	ri=i;
	rj=-(DELTAJ+1)-startpos[2]+totalsize[2]/2;

	  // outflow B^r
	  //	  prim[i][j][B1]=prim[ri][rj][B1]+(prim[ri][rj][B1]-prim[ri][rj-1][B1])*(j-rj);
	//	  prim[i][j][B1]=prim[ri][rj][B1];
	  // outflow B^\phi
	  //prim[i][j][B3]=prim[ri][rj][B3]+(prim[ri][rj][B3]-prim[ri][rj-1][B3])*(j-rj);
	//	  prim[i][j][B3]=prim[ri][rj][B3];
	
 	// define function that becomes 1.0 at TSTATIC
	//	if(t<=TSTATIC) fun = 1.0*sin(M_PI*0.5*t/TSTATIC);
	if(t<=TSTATIC) fun = 0.0;
	else fun=1.0;
	if(fun<0.0) fun=0.0;
	if(fun>1.0) fun=1.0;
	
	// smoothly turn on "outflow" boundary condition to avoid transients
	// fix B^r
	prim[i][j][B1]=(1.0-fun)*panalytic[i][j][B1] + fun*prim[ri][rj][B1];
	// fix B^\phi
	prim[i][j][B3]=(1.0-fun)*panalytic[i][j][B3] + fun*prim[ri][rj][B3];

	//	if(!((r>RDISKINNER)&&(r<RDISKOUTER))){
	//	  prim[i][j][B2]=prim[ri][rj][B2];
	//	}

#if(0)
	if(r<2.0*Rhor){
	  // outflow B^\theta
	  //	  prim[i][j][B2]=prim[ri][rj][B2];
	  prim[i][j][B2]=0; // monopole
	}
	else{
	  // fix B^\theta
	  prim[i][j][B2]=panalytic[i][j][B2];
	}
#endif


	// fix B^\theta
       prim[i][j][B2]=panalytic[i][j][B2];



#if(0)// outflow instead of evolved state
	// instead of letting evolution be default, let outflow'ed 4-velocity be default
	prim[i][j][U1]=prim[ri][rj][U1];
	//	prim[i][j][U2]=prim[ri][rj][U2];
	prim[i][j][U2]=0.0;
	prim[i][j][U3]=prim[ri][rj][U3];
#endif

	if(r<=RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)

	  get_geometry(ri, rj, CENT, &rgeom);
	   
	  prim[i][j][B2]=panalytic[i][j][B2];

 
	  Max_con(prim[ri][rj],&rgeom,Mcon);
	  MtoF(3,Fcov,&rgeom,Mcon);
	  romegaf2=Fcov[TT][TH]/Fcov[TH][PH];
	}

      }// end if within radius of disk
    }// end if one side of equator


    // other side of equator
    if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<=DELTAJ-1)&&(startpos[2]+j-totalsize[2]/2>=0) ) ){ // then at equator for an even-sized grid with 2 zones around equator for 4 zones total
	  
      if((r>RDISKINNER)&&(r<RDISKOUTER)){
	// OUTFLOW some quantities into conductor that is the current sheet
	ri=i;
	rj=(DELTAJ)-startpos[2]+totalsize[2]/2;

	// outflow B^r
	//	  prim[i][j][B1]=prim[ri][rj][B1]+(prim[ri][rj][B1]-prim[ri][rj-1][B1])*(j-rj);
	//	  prim[i][j][B1]=prim[ri][rj][B1];
	// outflow B^\phi
	//prim[i][j][B3]=prim[ri][rj][B3]+(prim[ri][rj][B3]-prim[ri][rj-1][B3])*(j-rj);
	//	  prim[i][j][B3]=prim[ri][rj][B3];
	
	// define function that becomes 1.0 at TSTATIC
	//	if(t<=TSTATIC) fun = 1.0*sin(M_PI*0.5*t/TSTATIC);
	if(t<=TSTATIC) fun = 0.0;
	else fun=1.0;
	if(fun<0.0) fun=0.0;
	if(fun>1.0) fun=1.0;
	
	// smoothly turn on "outflow" boundary condition to avoid transients
	// fix B^r
	prim[i][j][B1]=(1.0-fun)*panalytic[i][j][B1] + fun*prim[ri][rj][B1];
	// fix B^\phi
	prim[i][j][B3]=(1.0-fun)*panalytic[i][j][B3] + fun*prim[ri][rj][B3];


	//	if(!((r>RDISKINNER)&&(r<RDISKOUTER))){
	//	  prim[i][j][B2]=prim[ri][rj][B2];
	//	}

#if(0)
	if(r<2.0*Rhor){
	  // outflow B^\theta
	  //	  prim[i][j][B2]=prim[ri][rj][B2];
	  prim[i][j][B2]=0; // monopole
	}
	else{
	  // fix B^\theta
	  prim[i][j][B2]=panalytic[i][j][B2];
	}
#endif

	// fix B^\theta
	prim[i][j][B2]=panalytic[i][j][B2];


#if(0)// outflow instead of evolved state
	// instead of letting evolution be default, let outflow'ed 4-velocity be default
	prim[i][j][U1]=prim[ri][rj][U1];
	//	prim[i][j][U2]=prim[ri][rj][U2];
	prim[i][j][U2]=0.0;
	prim[i][j][U3]=prim[ri][rj][U3];
#endif


	if(r<=RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)

	  get_geometry(ri, rj, CENT, &rgeom);

	  prim[i][j][B2]=panalytic[i][j][B2];

	  Max_con(prim[ri][rj],&rgeom,Mcon);
	  MtoF(3,Fcov,&rgeom,Mcon);
	  romegaf2=Fcov[TT][TH]/Fcov[TH][PH];
	}


      }

#if(0)
      if((r>RDISKINNER)&&(r<RDISKOUTER)&&(i>=-N1BND)&&(i<N1+N1BND)&&( (abs(startpos[2]+j-(totalsize[2]-1)/2)<=DELTAJ) ) ){

	// OUTFLOW some quantities into conductor that is the current sheet

	// just average B^r
	prim[i][j][B1]=0.5*(prim[i][j+1][B1]+prim[i][j-1][B1]);
	prim[i][j][B3]=0.5*(prim[i][j+1][B3]+prim[i][j-1][B3]);
	  
	if(!((r>RDISKINNER)&&(r<RDISKOUTER))){
	  prim[i][j][B2]=0.5*(prim[i][j+1][B2]+prim[i][j-1][B2]);
	}

      }
#endif
    }


  }// end LOOPF2/1







  LOOPF2 LOOPF1{ // loop over entire domain looking for equator
    //	if(abs(startpos[2]+j-totalsize/2)<=2){ // then at equator
    if((i>=-N1BND)&&(i<N1+N1BND)&&( (startpos[2]+j-totalsize[2]/2<=DELTAJ-1)&&(startpos[2]+j-totalsize[2]/2>=-DELTAJ) ) ){ // then at equator for an even-sized grid with 2 zones around equator for 4 zones total

      if((r>RDISKINNER)&&(r<RDISKOUTER)){
      // outflow inner radial region till t=30, then try to fix
      /*
      if(
	 ((r>2.0*Rhor) || (t>30.0))&&
	 //	 ((r>3.0) || (t>30.0))&&
	 ((r>RDISKINNER)&&(r<RDISKOUTER))
	 ){
      */

	coord(i, j, CENT, X);
	bl_coord(X, &r, &th);
	get_geometry(i, j, CENT, &geom);

	if((r>RDISKINNER)&&(r<RDISKOUTER)){ // assume done above in split loops
	  // fix B^r, assume other field can do whatever they want.
	  //	  prim[i][j][B1]=panalytic[i][j][B1];
	  // fix B^\theta, assume other field can do whatever they want.
	  //	  prim[i][j][B2]=panalytic[i][j][B2];
	}

	// fix E_\phi=0, which completely determines all components of v^i
	// try it, but don't have to force it except at late times.
	// since near horizon there is the light surface and the grid goes through this surface, we can't guarantee that our v^i is time-like.
	Bcon[0]=0;
	Bcon[1]=prim[i][j][B1];
	Bcon[2]=prim[i][j][B2];
	Bcon[3]=prim[i][j][B3];
	  
	lower(Bcon,&geom,Bcov);
	Bsq=0.0;
	for(k=0;k<=3;k++) Bsq+=Bcon[k]*Bcov[k];

	////////////////
	// damp the calculation of v^i by damping Bsq
	//	if(r<2.0*Rhor) Bsq*=1.0/(1.0-exp(-t/10)+SMALL);
 
	if(r>RNEARBLACKHOLE){ // outflow omegaf2 (use it below to set v^i)
	  // must be same as in init.c
	  Omegastar=get_omegastar(&geom,r,th);
	}
	else{
	  Omegastar=romegaf2;
	}
 
	vcon[1]=-Bcon[1]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	if((r>RDISKINNER)&&(r<RDISKOUTER)) vcon[2]=-Bcon[2]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	else vcon[2]=0.0; // like AVOIDCS in phys.ffde.c (kinda stationary since v^\theta=0 on equator)

	vcon[3]=Omegastar-Bcon[3]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;

#if(LIMIT3VELBOUND)

	// guarenteed to return prim[i][j] that's time-like with Lorentz factor limited by GAMMAMAX
	limit_3vel(Bcon, &geom, vcon, prim[i][j]);

#else
	// check to make sure that chosen v^i is time-like.
	// allow modification of v^i since evolution toward stationary solution may not be adiabatic

	// get electric field from v and B ? Not sure how to get E^2(v^i,B^i)


	//	    for(k=1;k<=3;k++) dualfprintf(fail_file,"Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",k,Bcon[k],k,Bcov[k]);
	  
	//	    for(k=1;k<=3;k++) dualfprintf(fail_file,"vcon[%d]=%21.15g Bsq=%21.15g\n",k,vcon[k],Bsq);

	//	    for(k=1;k<=3;k++) if(isnan(vcon[k])){
	//	      dualfprintf(fail_file,"nan encountered: k=%d vcon=%21.15g\n",k,vcon[k]);
	//	    }
	// sticks primitive velocity into prim[i][j][U1->U3]
	whocalleducon=1;
	sanitycheck=1;
	if(vcon2pr(WHICHVEL,vcon,&geom,primtest)>=1){
	  failed=0; // don't fail
	  //	    dualfprintf(fail_file,"v^\phi[fixed]=%g\n",vcon[PH]);
	  //return(1);
	  // keep old velocity
	  sanitycheck=0;
	    // mark when this type of problem
	  if(DODEBUG){
	    TSCALELOOP failfloorcount[i][j][tscale][COUNTUTOPRIMFAILRHONEG]++;
	  }

	  ///////////////////////
	  //
	  // keeping old velocity leads to non-stationary solution and noise from current sheet dissipating leading to large u^t in evolved solution
	  //

	}
	else{ // then primtest is ok, but do sanity check first
	  for(k=U1;k<=U3;k++) if(isnan(primtest[k])){
	    dualfprintf(fail_file,"nan encountered: i=%d j=%d k=%d prim=%21.15g\n",i,j,k,primtest[k]);
	    sanitycheck=0;
	    // mark when this type of problem
	    if(DODEBUG){
	      TSCALELOOP failfloorcount[i][j][tscale][COUNTUTOPRIMFAILUNEG]++;
	    }
	  }
	  if(sanitycheck){ // then primtest ok, so assign
	    prim[i][j][U1]=primtest[U1];
	    prim[i][j][U2]=primtest[U2];
	    prim[i][j][U3]=primtest[U3];


	  }
	  // otherwise still avoid primtest

	}
	whocalleducon=0;

	// if no good stationary v^i, then reduce to AVOIDCS like behavior
	if(!sanitycheck){
	  prim[i][j][U2]=0.0;
	}

	// if set boundary condition, then ignore failure from inversion
	// code in FFDE mode doesn't actually use pflag.
	if(sanitycheck){
	  pflag[i][j][FLAGUTOPRIMFAIL]=UTOPRIMNOFAIL;
	}

#endif // end if old way of dealing with when v^i is not good (LIMIT3VELBOUND)


#if(0)
	/////////////////////////////
	//
	// additional limiting?
	//
	////////////////////////////

	if(sanitycheck){ // make sure still ok gamma
	  // force flow to not move too fast inside ergosphere
	  if(r<2) realgammamax=3;
	  else realgammamax=GAMMAMAX;
	      
	  // limit gamma
	  limit_gamma(realgammamax,prim[i][j],&geom,-1.0);
	}
#endif



      }

    }// end if at equator
  } // end loop over domain


  firsttime=0;
  return(0);
}



