
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define EXTRAP 0
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
  FTYPE primtest[NPR];

  FTYPE Bcon[NDIM],Bcov[NDIM],ucon[NDIM],Bsq;
  int sanitycheck;




  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    // DISK BOUNDARY CONDITION, modifies evolution region
    if(1){

      LOOPF2 LOOPF1{ // loop over entire domain looking for equator
	//	if(abs(startpos[2]+j-totalsize/2)<=2){ // then at equator
	if((i>=0)&&(i<N1)&&( (startpos[2]+j-totalsize[2]/2<=1)&&(startpos[2]+j-totalsize[2]/2>=-2) ) ){ // then at equator for an even-sized grid with 2 zones around equator for 4 zones total

	  get_geometry(i, j, CENT, &geom);


	  // fix B^\theta, assume other field can do whatever they want.
	  prim[i][j][B1]=panalytic[i][j][B1];
	  prim[i][j][B2]=panalytic[i][j][B2];


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

	  vcon[1]=-Bcon[1]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	  vcon[2]=-Bcon[2]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	  vcon[3]=Omegastar-Bcon[3]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;

	  //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",k,Bcon[k],k,Bcov[k]);
	  
	  //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"vcon[%d]=%21.15g Bsq=%21.15g\n",k,vcon[k],Bsq);

	  //	    for(k=1;k<=3;k++) if(isnan(vcon[k])){
	  //	      dualfprintf(fail_file,"nan encountered: k=%d vcon=%21.15g\n",k,vcon[k]);
	  //	    }
	  // sticks primitive velocity into prim[i][j][U1->U3]
	  if(vcon2pr(WHICHVEL,vcon,&geom,primtest)>=1){
	    failed=0; // don't fail
//	    dualfprintf(fail_file,"v^\phi[fixed]=%g\n",vcon[PH]);
	    //return(1);
	    // keep old primitives
	  }
	  else{ // then primtest is ok, but do sanity check first
	    sanitycheck=1;
	    for(k=U1;k<=B3;k++) if(isnan(primtest[k])){
	      dualfprintf(fail_file,"nan encountered: i=%d j=%d k=%d prim=%21.15g\n",i,j,k,primtest[k]);
	      sanitycheck=0;
	    }
	    if(sanitycheck){ // then primtest ok, so assign
	      prim[i][j][U1]=primtest[U1];
	      prim[i][j][U2]=primtest[U2];
	      prim[i][j][U3]=primtest[U3];
	    }
	    // otherwise still avoid primtest
	    
	  }
	  
	    
	  //	    dualfprintf(fail_file,"i=%d j=%d\n",geom.i,geom.j);
	  

	}// end if at equator
      } // end loop over domain
    }// end if DISK BOUNDARY CONDITION











    if (mycpupos[1] == 0) {
      if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
	/* inner r boundary condition: u, just copy */
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--)	  PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
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
    }// end if mycpupos==0

    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
	/* outer r BC: outflow */
      
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++)	  PBOUNDLOOP prim[i][j][k] = prim[ri][rj][k];
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
      else if((BCtype[X1UP]==FIXED)||(BCtype[X1UP]==FIXED)){
	// outer radial fixed boundary conditions
	for (j = -N2BND; j < N2+N2BND; j++)
	  for(i=N1;i<=N1-1+N1BND;i++)
	    PBOUNDLOOP prim[i][j][k] = p[i][j][k];
      }
    }

    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=0;
	for(j=-N2BND;j<=-1;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j-1)][k];
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
       for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=N2-1;
	for(j=N2;j<=N2-1+N2BND;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j+1)][k];
      }
    }
    /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
    /* inner pole */
    if (mycpupos[2] == 0) {
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
    }
    /* outer pole */
    if (mycpupos[2] == ncpux2 - 1) {
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
    }


  
  

  }// end if stage0 or stagem1

  if (USEMPI) bound_mpi(boundstage, prim);


  return (0);
}
