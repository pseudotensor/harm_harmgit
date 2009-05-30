
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define EXTRAP 2
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *ptrgeom,FTYPE *p2interp);
  SFTYPE Dt=-1.0; // no boundary zone should ever be in accounting
  int i, j, k;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
  FTYPE ucon[NDIM];
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j
  FTYPE prescale[NPR];
  FTYPE rescaled[NBIGBND][NPR];
  FTYPE Bcon[NDIM],Bcov[NDIM];
  FTYPE copied[NBIGBND][NPR];


  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

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
      else if((BCtype[X1DN]==NSSURFACE)){
	// same as above but don't modify magnetid field
	// fix field
	// outflow v and densities 


	/* inner r boundary condition: u, just copy */
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++)  prim[i][j][k] = prim[ri][rj][k];
	    vcon[RR]=0; // surface that completely dissipates normal direction momentum
	    vcon[TH]=0; // "" for this component
	    vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	    get_geometry(i, j, CENT, &geom);
	    vcon2pr(WHICHVEL, vcon,&geom,prim[i][j]); // get in terms of primitive velocity

	    /*
	    // outflow field, NS field completely reconnects through surface
	    for(k=B3;k<=B3;k++){
	      //	      prim[i][j][k] = prim[ri][rj][k];
	      if(i==-1) prim[i][j][k] = prim[ri][rj][k] -(prim[ri+1][rj][k]- prim[ri][rj][k]) ;
	      else if(i==-2) prim[i][j][k] = prim[ri][rj][k]-2.0(prim[ri+2][rj][k]- prim[ri+1][rj][k]) ;
	      else if(i==-3) prim[i][j][k] = prim[ri][rj][k]-3.0(prim[ri+3][rj][k]- prim[ri+2][rj][k]) ;

	    }
	    */
	    // fixed field
	    //	    for(k=B1;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];

	    for(k=B3;k<=B3;k++)  prim[i][j][k] = prim[ri][rj][k];
	  }
#elif(EXTRAP==1)

	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    // set  	    ucon[TT,etc.]
	    vcon[RR]=0; // surface that completely dissipates normal direction momentum
	    vcon[TH]=0; // "" for this component
	    vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	    get_geometry(i, j, CENT, &geom);
	    vcon2pr(WHICHVEL,vcon,&geom,prim[i][j]); // get in terms of primitive velocity

	    /*
	    // outflow field, NS field completely reconnects through surface
	    for(k=B3;k<=B3;k++){
	      //	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	      if(i==-1) prim[i][j][k] = prim[ri][rj][k] -(prim[ri+1][rj][k] * gdet[ri+1][rj][CENT]- prim[ri][rj][k] * gdet[ri][rj][CENT])/gdet[i][j][CENT] ;
	      else if(i==-2) prim[i][j][k] = prim[ri][rj][k] -2.0*(prim[ri+2][rj][k] * gdet[ri+2][rj][CENT]- prim[ri+1][rj][k] * gdet[ri+1][rj][CENT])/gdet[i][j][CENT] ;
	      else if(i==-3) prim[i][j][k] = prim[ri][rj][k] -3.0*(prim[ri+3][rj][k] * gdet[ri+3][rj][CENT]- prim[ri+2][rj][k] * gdet[ri+2][rj][CENT])/gdet[i][j][CENT] ;
	    }
	    */
	    // fixed field
	    for(k=B1;k<=B3;k++)  prim[i][j][k] = p[i][j][k];

	  }
#elif(EXTRAP==2)
	  ri=0;
	  rj=j;
	  //	  dualfprintf(fail_file,"got here1\n");
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  //	  dualfprintf(fail_file,"got here1.5\n");
	  //	  PBOUNDLOOP dualfprintf(fail_file,"list: k=%d prim=%g prescale=%g\n",k,prim[ri][rj][k],prescale[k]);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess and get copied version
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    if(rescale(-1,1,rescaled[i+N1BND],&geom,prescale)>=1){
	      // then interpolate of conserved quantities didn't work out, so revert to copied
	      PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k];
	    }
	  }
	  //	  dualfprintf(fail_file,"got here2\n");
	  // now rescaled[][NPR] contains interpolated quantities as primitive-type
	  // these are guarenteed to preserve v.B=0 and E.B=0, can't screw that up

	  // now work with them

	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){
	    prim[i][j][B3]=rescaled[i+N1BND][B3];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
	  }

	  // forced
	  vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	  
	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);
	    
	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    Bcon[3]=prim[i][j][B3];
	    lower(Bcon,&geom,Bcov);

	    // get u^\mu from interpolation
	    if(ucon_calc(rescaled[i+N1BND],&geom,ucon)>=1){
	      // then get from copying primitives, which cannot fail
	      ucon_calc(copied[i+N1BND],&geom,ucon);
	    }

	    //	    dualfprintf(fail_file,"got here3\n");

	    if(Bcov[RR]!=0.0){
	      // then copy v^\theta and fix v^r
	      vcon[TH]=ucon[TH]/ucon[TT];
	      vcon[RR]=-(vcon[TH]*Bcov[TH]+vcon[PH]*Bcov[PH])/(Bcov[RR]);
	    }
	    else if(Bcov[TH]!=0.0){
	      // then copy v^r and fix v^\theta
	      vcon[RR]=ucon[RR]/ucon[TT];
	      vcon[TH]=-(vcon[RR]*Bcov[RR]+vcon[PH]*Bcov[PH])/(Bcov[TH]);
	    }
	    else{
	      dualfprintf(fail_file,"So Bcov[TH,RR]=0, but v^\phi!=0, so B_\phi=0?\n");
	      myexit(1);
	    }

	    // sticks primitive velocity into prim[i][j][U1->U3]
	    vcon2pr(WHICHVEL,vcon,&geom,prim[i][j]); 
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

    }

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
      /* if fixed BC: do nothing */
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
	    prim[i][j][U2] *= 1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= 1.;
	    prim[i][j][B3] *= 1.;
	  }
	  else{
	    prim[i][j][U2] *= 1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= 1.;
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
	  prim[i][j][U2] *= 1.;
	  prim[i][j][U3] *= 1.;
	  prim[i][j][B1] *= 1.;
	  prim[i][j][B2] *= 1.;
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
