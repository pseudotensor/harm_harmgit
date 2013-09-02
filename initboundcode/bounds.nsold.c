
#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: something seriously wrong with EXTRAP=1 (EOMFFDE)

// FILE SIMILAR TO for fishmon but with NSSURFACE for inner radial boundary

#define EXTRAP 1
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
#define POLEDEATH 0
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole


// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPN2 LOOPN3 LOOPBOUNDIN1 & LOOPBOUNDOUT1
// X2 in&out: LOOPF1 LOOPN3 LOOPBOUNDIN2 & LOOPBOUNDOUT2  // LOOPF1 ok if X1 dimension not there, then LOOPF1->LOOPN1
// X3 in&out: LOOPF1 LOOPF2 LOOPBOUNDIN3 & LOOPBOUNDOUT3  // as above



int bound_prim_user(int boundstage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE V[NDIM],X[NDIM];
  FTYPE myomega;
  FTYPE Dt;


  myomega=Omegastar/(2.0*M_PI);
  k=0;
  Dt=1;


  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    if (mycpupos[1] == 0) {
      if((BCtype[X1DN]==NSSURFACE)){
        // same as above but don't modify magnetid field
        // fix field
        // outflow v and densities 


        /* inner r boundary condition: u, just copy */
        for (j = 0; j < N2; j++) {

          ri=0;
          rj=j;
          for(i=-1;i>=-N1BND;i--){
            for(pl=RHO;pl<=UU;pl++){
              MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,k)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
            }
            // set       ucon[TT,etc.]
            vcon[RR]=0; // surface that completely dissipates normal direction momentum
            vcon[TH]=0; // "" for this component
            vcon[PH]=myomega; // surface rotates with angular frequency myomega to observer at infinity
            get_geometry(i, j, k, CENT, &geom);
            vcon2pr(WHICHVEL,vcon,&geom,MAC(prim,i,j,k)); // get in terms of primitive velocity

            // outflow field, NS field completely reconnects through surface
            for(pl=B1;pl<=B3;pl++){
              MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,k)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
            }

          }


          for(i=-1;i>=-N1BND;i--){
#if(WHICHVEL==VEL4)
            get_geometry(i, j, k, CENT, &geom);
            inflow_check_4vel(1,MAC(prim,i,j,k),&geom, Dt) ;
#elif(WHICHVEL==VEL3)
            get_geometry(i, j, k, CENT, &geom);
            inflow_check_3vel(1,MAC(prim,i,j,k),&geom, Dt) ;
            // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
            if(jonchecks){
              //fixup1zone(MAC(prim,i,j,k),&geom,0);
              failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),&geom,-3);
              if(failreturn){
                dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
                if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
              }
            }
#endif
#elif(WHICHVEL==VELREL4)
            get_geometry(i,j,k,CENT,&geom) ;
            inflow_check_rel4vel(1,MAC(prim,i,j,k),&geom,Dt) ;
            if(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,MAC(prim,i,j,k),&geom,Dt)>=1)
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
          ri=N1-1;
          rj=j;
          for(i=N1;i<=N1-1+N1BND;i++){
            for(pl=RHO;pl<=UU;pl++){
              MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,k)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
            }
            pl=U1; // treat U1 as special
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * (1. - 2*(i-ri)*dx[1]) ;
            for(pl=U2;pl<=U3;pl++){
              MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * (1. - (i-ri)*dx[1]) ;
            }
            pl=B1; // treat B1 special
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,k)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
            for(pl=B2;pl<=B3;pl++){
              MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl) * (1. - (i-ri)*dx[1]) ;
            }
          }




          for(i=N1;i<=N1-1+N1BND;i++){
#if(WHICHVEL==VEL4)
            get_geometry(i, j, k, CENT, &geom);
            inflow_check_4vel(1,MAC(prim,i,j,k),&geom,Dt) ;
#elif(WHICHVEL==VEL3)
            get_geometry(i, j, k, CENT, &geom);
            inflow_check_3vel(1,MAC(prim,i,j,k),&geom,Dt) ;
            // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
            if(jonchecks){
              //fixup1zone(MAC(prim,i,j,k),&geom,0);
              failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),&geom,-3);
              if(failreturn){
                dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
                if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
              }
            }
#endif
#elif(WHICHVEL==VELREL4)
            get_geometry(i,j,k,CENT,&geom) ;
            inflow_check_rel4vel(1,MAC(prim,i,j,k),&geom,Dt) ;
            if(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,MAC(prim,i,j,k),&geom, Dt)>=1)
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
        for(j=-N2BND;j<=-1;j++) PLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j-1),k,pl);
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -N1BND; i <=N1-1+N1BND; i++){
        ri=i;
        rj=N2-1;
        for(j=N2;j<=N2-1+N2BND;j++) PLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j+1),k,pl);
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
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
        }
      }
    }
    /* outer pole */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -N1BND; i <= N1-1+N1BND; i++) {
        for (j = N2; j <= N2-1+N2BND; j++) {
          if(POSDEFMETRIC==0){
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= -1.;
            MACP0A1(prim,i,j,k,U3) *= 1.;
            MACP0A1(prim,i,j,k,B2) *= -1.;
            MACP0A1(prim,i,j,k,B3) *= 1.;
          }
        }
      }
    }

    // to help protect the pole from death blows to the computational grid
    // a sort of crushing regularization
#define POLEDEATH 0
    // causes problems with stability at just beyond pole
    // for field line plots, can just set B^\theta=0 along pole


#if(POLEDEATH)

    /* inner pole */
    if (mycpupos[2] == 0) {
      for (i = -N1BND; i <= N1-1+N1BND; i++) {
        for (j = -N2BND; j < 0+POLEDEATH; j++) {
          if(POSDEFMETRIC==0){
            // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
            // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
            MACP0A1(prim,i,j,k,U2) *= 0;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= 0.;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
        }
      }
    }
    /* outer pole */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -N1BND; i <= N1-1+N1BND; i++) {
        for (j = N2-POLEDEATH; j <= N2-1+N2BND; j++) {
          if(POSDEFMETRIC==0){
            MACP0A1(prim,i,j,k,U2) *= 0.;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
          else{
            MACP0A1(prim,i,j,k,U2) *= 0.;
            MACP0A1(prim,i,j,k,B2) *= 0.;
          }
        }
      }
    }
#endif




  }// end if stage0 or stagem1





  return (0);
}




// see interpline.c
int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM])
{
  int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM]);

  if(doinverse==0){
    flip_y(iterglobal, recontype, bs, be, yin);
  }
  else{
    flip_y(iterglobal, recontype, bs, be, yin);
    flip_y(iterglobal, recontype, bs, be, yout);
  }

  return(0);

}


#include "reconstructeno.h"

int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM])
{
  int pl,myi;


#if( WENO_DIR_FLIP_CONS_SIGN_DN )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_DN && (recontype == CVT_C2A || recontype == CVT_A2C) && mycpupos[iterglobal] == 0 ) { 
    PLOOP(pliter,pl) 
      for( myi = bs; myi < 0; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif
 
#if( WENO_DIR_FLIP_CONS_SIGN_UP )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_UP && (recontype == CVT_C2A || recontype == CVT_A2C)  && mycpupos[iterglobal] == numbercpu[iterglobal] - 1 ) { 
    PLOOP(pliter,pl) 
      for( myi = N1*(iterglobal==1) + N2*(iterglobal==2) + N3*(iterglobal==3); myi <= be; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif


  return(0);

}
