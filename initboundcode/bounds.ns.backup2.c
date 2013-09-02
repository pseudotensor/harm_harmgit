
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
//#define POLEDEATH N2BND
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
  int i, j, k, pl, pl2;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE rdxdxp[NDIM][NDIM];
  FTYPE V[NDIM],X[NDIM];
  FTYPE rV[NDIM],rX[NDIM];
  extern int OBtopr_general(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr);
  FTYPE Bcon[NDIM],prnew[NPR];
  FTYPE Bconnew[NDIM];
  extern void getnewucon(FTYPE *uconmetpin, FTYPE *rV, struct of_geom *rptrgeom, FTYPE (*rdxdxp)[NDIM], FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetpout);

  //  dualfprintf(fail_file,"primmem: %ld %d\n",prim,prim);

  ///////////////////////////
  //
  // X1 inner OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////

  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
      /* inner r boundary condition: u, just copy */
      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#elif(EXTRAP==1)

        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN{
          for(pl=RHO;pl<=UU;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          }
          pl=U1;     // treat U1 as special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
          for(pl=U2;pl<=U3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. + (i-ri)*dx[1]) ;
          }
          pl=B1; // treat B1 special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          for(pl=B2;pl<=B3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. + (i-ri)*dx[1]) ;
          }
        }
#elif(EXTRAP==2)
        ri=0;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,MAC(prim,ri,rj,rk),&rgeom,prescale);
        LOOPBOUND1IN{
          // set guess
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,k,pl);
          get_geometry(i, j, k, CENT, &geom);     
          rescale(-1,1,MAC(prim,i,j,k),&geom,prescale);
        }
#endif



        LOOPBOUND1IN{
#if(WHICHVEL==VEL4)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_4vel(1,MAC(prim,i,j,k),&geom, 0) ;
#elif(WHICHVEL==VEL3)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_3vel(1,MAC(prim,i,j,k),&geom, 0) ;
          // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
          if(jonchecks){
            //fixup1zone(MAC(prim,i,j,k),&geom,0);
            failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),&geom,-3);
            if(failreturn){
              dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
              if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
            }
          }
#endif
#elif(WHICHVEL==VELREL4)
          get_geometry(i,j,k,CENT,&geom) ;
          inflow_check_rel4vel(1,MAC(prim,i,j,k),&geom,0) ;
          if(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,MAC(prim,i,j,k),&geom,0)>=1)
            FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif 
        }
      }// end 2 3
    }// end if correct bound type
    else if((BCtype[X1DN]==NSSURFACE)){
      // same as above but don't modify magnetid field
      // fix field
      // outflow v and densities 


      /* inner r boundary condition: u, just copy */
      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN{


          get_geometry(i, j, k, CENT, &geom);
          // get in terms of primitive velocity
          coord(i, j, k, CENT, X);
          bl_coord( X, V );
          dxdxprim(X, V, dxdxp);

          //   for(pl=RHO;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          for(pl=UU;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          // outflow field, NS field completely reconnects through surface
          //   for(pl=B1;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl);

          // outflow B^\theta B^\phi
          for(pl=B2;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          // fix B^r
          MACP0A1(prim,i,j,k,B1) = MACP0A1(panalytic,i,j,k,B1);
   
          // fix \rho_0
          MACP0A1(prim,i,j,k,RHO) = MACP0A1(panalytic,i,j,k,RHO);


          Bcon[0]=0;
          Bcon[1]=MACP0A1(prim,i,j,k,B1);
          Bcon[2]=MACP0A1(prim,i,j,k,B2);
          Bcon[3]=MACP0A1(prim,i,j,k,B3);


          ///////////////////////////////
          //
          // new way to get velocity
          // surface rotates with angular frequency Omegastar to observer at infinity
          if(OBtopr_general(Omegastar/dxdxp[3][3],Bcon,&geom,prnew)>=1){
            dualfprintf(fail_file, "OBtopr(bounds_disk): space-like error in init_postfield()\n");
            dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
            failed=1;
            return(1);
          }
          // assign answer
          MACP0A1(prim,i,j,k,U1)=prnew[U1];
          MACP0A1(prim,i,j,k,U2)=prnew[U2];
          MACP0A1(prim,i,j,k,U3)=prnew[U3];

#if(0)

          vcon[RR]=0; // surface that completely dissipates normal direction momentum
          vcon[TH]=0; // "" for this component

          // below assumes no phi mixing with other directions in grid
          vcon[PH]=Omegastar/dxdxp[3][3]; // surface rotates with angular frequency Omegastar to observer at infinity

          //   trifprintf("i=%d j=%d k=%d ri=%d rj=%d rk=%d  :: MACP0A1(prim,-2,8,0,RHO)=%21.15g\n",i,j,k,ri,rj,rk,MACP0A1(prim,-2,8,0,RHO));

          MYFUN(vcon2pr(WHICHVEL,vcon,&geom,MAC(prim,i,j,k)),"bounds.ns.c:bound_prim_user()", "vcon2pr()", 0);
#endif


        }
#elif(EXTRAP==1)
        ri=0;
        rj=j;
        rk=k;
        LOOPBOUND1IN{

          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord( X, V );
          dxdxprim(X, V, dxdxp);

          get_geometry(ri, rj, rk, CENT, &rgeom);
          coord(ri, rj, rk, CENT, rX);
          bl_coord( rX, rV );
          dxdxprim(rX, rV, rdxdxp);


          //   for(pl=RHO;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          for(pl=UU;pl<=UU;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          // outflow field, NS field completely reconnects through surface
          //   for(pl=B1;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,k,pl);

          // outflow B^\theta B^\phi
          // SUPERGODMARK : removed B2 -- worried about divb ... probably just normalization is crazy
          //   MACP0A1(prim,i,j,k,B2) = MACP0A1(panalytic,i,j,k,B2);

          //   for(pl=B2;pl<=B3;pl++){
          //  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          // }
          //pl = B2;
          // for dipole, constant is B^\theta \propto 1/(\detg r^2)
          // below seems best
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)*rV[1]*rV[1]/(GLOBALMACP1A0(gdet,CENT,i,j,k)*V[1]*V[1]) ;
          // seems to do ok
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/(GLOBALMACP1A0(gdet,CENT,i,j,k)) ;
          // does worst as for polar artifacts
          //   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) ;

#if(0)
          // reference field
          Bcon[0]=0;
          Bcon[1]=MACP0A1(prim,ri,rj,rk,B1);
          Bcon[2]=MACP0A1(prim,ri,rj,rk,B2);
          Bcon[3]=MACP0A1(prim,ri,rj,rk,B3);

          getnewucon(Bcon, rV, &rgeom, rdxdxp, V, &geom, dxdxp, Bconnew);

          MACP0A1(prim,i,j,k,B2)=Bconnew[TH];
          MACP0A1(prim,i,j,k,B3)=Bconnew[PH];
#endif
#if(0)
          pl = B3;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/(GLOBALMACP1A0(gdet,CENT,i,j,k)) ;

          pl=B2;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;

          //pl=B3;
          //MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;


          //   for(pl=B2;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#endif
#if(1)
          // MOST ROBUST
          pl = B3;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          pl=B2;
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
#endif

          // fix B^r
          MACP0A1(prim,i,j,k,B1) = MACP0A1(panalytic,i,j,k,B1);
   
          // fix \rho_0
          MACP0A1(prim,i,j,k,RHO) = MACP0A1(panalytic,i,j,k,RHO);


          Bcon[0]=0;
          Bcon[1]=MACP0A1(prim,i,j,k,B1);
          Bcon[2]=MACP0A1(prim,i,j,k,B2);
          Bcon[3]=MACP0A1(prim,i,j,k,B3);

#if(0)

          //   dualfprintf(fail_file,"Omegastar=%21.15g dxdxp[3][3]=%21.15g B1=%21.15g B2=%21.15g B3=%21.15g\n",Omegastar,dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);

          ///////////////////////////////
          //
          // new way to get velocity
          // surface rotates with angular frequency Omegastar to observer at infinity
          if(OBtopr_general(Omegastar/dxdxp[3][3],Bcon,&geom,prnew)>=1){
            dualfprintf(fail_file, "OBtopr(bounds): space-like error in init_postfield()\n");
            dualfprintf(fail_file,"Cannot continue without 4-velocity!\n");
            failed=1;
            return(1);
          }
          // assign answer
          MACP0A1(prim,i,j,k,U1)=prnew[U1];
          MACP0A1(prim,i,j,k,U2)=prnew[U2];
          MACP0A1(prim,i,j,k,U3)=prnew[U3];

          //   if(t>1.9 && t<2.1){
          //dualfprintf(fail_file,"t=%21.15g i=%d j=%d\n",t,i,j);
          //  dualfprintf(fail_file,"newus: %21.15g %21.15g %21.15g\n",prnew[U1],prnew[U2],prnew[U3]);
          //  dualfprintf(fail_file,"Omegastar'=%21.15g Bcon1=%21.15g Bcon2=%21.15g Bcon3=%21.15g\n",Omegastar/dxdxp[3][3],Bcon[1],Bcon[2],Bcon[3]);
          // }

#endif

#if(1)
          // set       ucon[TT,etc.]
          vcon[RR]=0; // surface that completely dissipates normal direction momentum
          vcon[TH]=0; // "" for this component
          // below assumes no phi mixing with other directions in grid
          vcon[PH]=Omegastar/dxdxp[3][3]; // surface rotates with angular frequency Omegastar to observer at infinity

          // get in terms of primitive velocity
          MYFUN(vcon2pr(WHICHVEL,vcon,&geom,MAC(prim,i,j,k)),"bounds.ns.c:bound_prim_user()", "vcon2pr()", 1);
#endif

          PLOOP(pliter,pl){
            if(!isfinite(MACP0A1(prim,i,j,k,pl))){
              dualfprintf(fail_file,"i=%d j=%d steppart=%d nstep=%ld\n",i,j,steppart,nstep);
              PLOOP(pliter,pl2) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl2,MACP0A1(prim,i,j,k,pl2));
              myexit(2525);
            }
          }


        }
#elif(EXTRAP==2)
        ri=0;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,MAC(prim,ri,rj,rk),&rgeom,MAC(prim,ri,rj,rk));
        LOOPBOUND1IN{
          for(pl=RHO;pl<=UU;pl++)   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
          for(pl=U1;pl<=U3;pl++)   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl); // not really used

          // outflow all (crushing field)
          //   for(pl=B1;pl<=B3;pl++)   MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          // outflow B^\theta B^\phi
          for(pl=B2;pl<=B3;pl++)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          get_geometry(i, j, k, CENT, &geom);
          rescale(-1,1,MAC(prim,i,j,k),&geom,MAC(prim,i,j,k)); // assumes B^r not used to rescale if B^r fixed
        }
        rescale(-1,1,MAC(prim,ri,rj,rk),&rgeom,MAC(prim,ri,rj,rk)); 

        // fix B^r
        MACP0A1(prim,i,j,k,B1) = MACP0A1(panalytic,i,j,k,B1);

        LOOPBOUND1IN{
          get_geometry(i, j, k, CENT, &geom);
          coord(i, j, k, CENT, X);
          bl_coord( X, V );
          dxdxprim(X, V, dxdxp);

          // now fix velocity
          vcon[RR]=0; // surface that completely dissipates normal direction momentum
          vcon[TH]=0; // "" for this component
          // below assumes no phi mixing with other directions in grid
          vcon[PH]=Omegastar/dxdxp[3][3]; // surface rotates with angular frequency Omegastar to observer at infinity
   

          // get in terms of primitive velocity
          MYFUN(vcon2pr(WHICHVEL,vcon,&geom,MAC(prim,i,j,k)),"bounds.ns.c:bound_prim_user()", "vcon2pr()", 2);
        }

#endif


      }
    }




  }// end if mycpupos[1]==0


  ///////////////////////////
  //
  // X1 outer OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////


  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
      /* outer r BC: outflow */

      LOOPN2 LOOPN3{
#if(EXTRAP==0)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
#elif(EXTRAP==1)
        ri=N1-1;
        rj=j;
        rk=k;
        LOOPBOUND1OUT{
          for(pl=RHO;pl<=UU;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          }
          pl=U1; // treat U1 as special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - 2*(i-ri)*dx[1]) ;
          for(pl=U2;pl<=U3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
          }
          pl=B1; // treat B1 special
          MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * GLOBALMACP1A0(gdet,CENT,ri,rj,rk)/GLOBALMACP1A0(gdet,CENT,i,j,k) ;
          for(pl=B2;pl<=B3;pl++){
            MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl) * (1. - (i-ri)*dx[1]) ;
          }
        }
#elif(EXTRAP==2)
        ri=N1-1;
        rj=j;
        rk=k;
        get_geometry(ri, rj, rk, CENT, &rgeom);
        rescale(1,1,MAC(prim,ri,rj,rk),&rgeom,prescale);
        LOOPBOUND1OUT{
          // set guess
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(prim,ri,rj,rk,pl);
          get_geometry(i, j, k, CENT, &geom);
          rescale(-1,1,MAC(prim,i,j,k),&geom,prescale);
        }
#endif

        LOOPBOUND1OUT{
#if(WHICHVEL==VEL4)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_4vel(1,MAC(prim,i,j,k),&geom,0) ;
#elif(WHICHVEL==VEL3)
          get_geometry(i, j, k, CENT, &geom);
          inflow_check_3vel(1,MAC(prim,i,j,k),&geom,0) ;
          // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
          if(jonchecks){
            //fixup1zone(MAC(prim,i,j,k),&geom,0);
            failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),&geom,-3);
            if(failreturn){
              dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
              if (fail(i,j,k,FAIL_BCFIX) >= 1) return (1);
            }
          }
#endif
#elif(WHICHVEL==VELREL4)
          get_geometry(i,j,k,CENT,&geom) ;
          inflow_check_rel4vel(1,MAC(prim,i,j,k),&geom,0) ;
          if(limit_gamma(0,GAMMAMAX,GAMMAMAXRAD,MAC(prim,i,j,k),&geom, 0)>=1)
            FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif 
        }
      }// end 2 3
    }// end if correct bound type
  }// end if mycpu is correct


  ///////////////////////////
  //
  // X2 inner POLARAXIS
  //
  ///////////////////////////


  /* inner polar BC (preserves u^t rho and u) */
  if (mycpupos[2] == 0) {
    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
      LOOPF1 LOOPN3{
        ri=i;
        rj=0;
        rk=k;
        LOOPBOUND2IN{
          PBOUNDLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j-1),rk,pl);
          // dualfprintf(fail_file,"i=%d j=%d ri=%d rj=%d  :: rj+(rj-j-1)=%d\n",i,j,ri,rj,rj+(rj-j-1));
        }
      }
    }

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPN3{
        LOOPBOUND2IN {
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
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPN3 {
        for (j = 0-POLEDEATH; j < 0+POLEDEATH; j++) {
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
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==0


  ///////////////////////////
  //
  // X2 outer POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == ncpux2-1) {
    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
      LOOPF1 LOOPN3{
        ri=i;
        rj=N2-1;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pliter,pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j+1),rk,pl);
      }
    }

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPF1 LOOPF3{
        LOOPBOUND2OUT {
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
      }// end loop 13

#if(POLEDEATH)
      // fixup
      LOOPF1 LOOPN3 {
        for (j = N2-1+1-POLEDEATH; j <= N2-1+POLEDEATH; j++) {
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
      }// end loop 13
#endif

    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==ncpux2-1



  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another

      LOOPF1 LOOPF2{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=j;
        rk=N3;
        LOOPBOUND3IN PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+k,pl);

        // copy from lower side to upper boundary zones
        ri=i;
        rj=j;
        rk=0;
        LOOPBOUND3OUT PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(k-N3),pl);
      }
    }
  }





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


///Called after the MPI boundary routines
int bound_prim_user_after_mpi(int boundstage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  return(0);
}

