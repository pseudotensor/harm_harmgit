
/*! \file fluxcompute.c
    \brief Per-point calculations for getting flux.

    // ideas:
    // GODMARK: Consider ZIP-average for non-diffusive part of LAXF or HLL flux
    // Chacon (2004) Computer Physics Communications

*/


#include "decs.h"



static void diag_fluxdump_1(int dir, int i, int j, int k, FTYPE *p_l, FTYPE *p_r, FTYPE *F_l, FTYPE *F_r);


///////////////////////////////////
///
/// actually compute the flux and ctop for dt calculation
///
//////////////////////////////////
/// actually compute flux from given data
int flux_compute_general(int i, int j, int k, int dir, struct of_geom *ptrgeom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *ctopallptr)
{
  int flux_compute(int i, int j, int k, int dir, struct of_geom *geom, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE ctopmhd, FTYPE *cminmaxrad_l, FTYPE *cminmaxrad_r, FTYPE *cminmaxrad, FTYPE ctoprad, FTYPE CUf, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, FTYPE *F);
  int p2SFUevolve(int dir, int isleftright, FTYPE *p, struct of_geom *geom, struct of_state **ptrstate, FTYPE *F, FTYPE *U);
  FTYPE cminmax_l[NUMCS], cminmax_r[NUMCS], cminmax[NUMCS];
  FTYPE cminmaxrad_l[NUMCS], cminmaxrad_r[NUMCS], cminmaxrad[NUMCS];
  FTYPE cminmaxrad2_l[NUMCS], cminmaxrad2_r[NUMCS], cminmaxrad2[NUMCS];
  FTYPE ctopmhd,ctoprad,ctoprad2;
  struct of_state state_c, state_l, state_r;
  struct of_state *ptrstate_c, *ptrstate_l, *ptrstate_r;
  FTYPE F_c[NPR], F_l[NPR], F_r[NPR];
  FTYPE U_c[NPR], U_l[NPR], U_r[NPR];
  int pl,pliter;
  int jj,kk;

  
  // default
  ptrstate_c = &state_c;
  ptrstate_l = &state_l;
  ptrstate_r = &state_r;



  //////////////////////
  //
  // setup flux calculation based upon interpolated primitives
  //
  //////////////////////
  MYFUN(p2SFUevolve(dir, ISLEFT, p_l, ptrgeom, &ptrstate_l, F_l, U_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
  MYFUN(p2SFUevolve(dir, ISRIGHT, p_r, ptrgeom, &ptrstate_r, F_r, U_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);


#if(FLUXDUMP==2)
  FTYPE FPAKE[NPR], FPAKE_l[NPR], FPAKE_r[NPR], UPAKE_l[NPR], UPAKE_r[NPR];
  FTYPE FEN[NPR], FEN_l[NPR], FEN_r[NPR], UEN_l[NPR], UEN_r[NPR];
  FTYPE FEM[NPR], FEM_l[NPR], FEM_r[NPR], UEM_l[NPR], UEM_r[NPR];
  FTYPE FRAD[NPR], FRAD_l[NPR], FRAD_r[NPR], URAD_l[NPR], URAD_r[NPR];
  int p2SFUevolveall(int dir, int isleftright, FTYPE *p, struct of_geom *ptrgeom, struct of_state **ptrstate, FTYPE *FPAKE, FTYPE *FEN, FTYPE *FEM, FTYPE *FRAD, FTYPE *UPAKE, FTYPE *UEN, FTYPE *UEM, FTYPE *URAD);
  MYFUN(p2SFUevolveall(dir, ISLEFT, p_l, ptrgeom, &ptrstate_l, FPAKE_l, FEN_l, FEM_l, FRAD_l, UPAKE_l, UEN_l, UEM_l, URAD_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
  MYFUN(p2SFUevolveall(dir, ISRIGHT, p_r, ptrgeom, &ptrstate_r, FPAKE_r, FEN_r, FEM_r, FRAD_r, UPAKE_r, UEN_r, UEM_r, URAD_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);
#endif



  // usually "always" need cminmax_l cminmax_r and always need ctop
  get_wavespeeds(dir, ptrgeom, p_l, p_r, U_l, U_r, F_l, F_r, ptrstate_l, ptrstate_r, cminmax_l, cminmax_r, cminmax, &ctopmhd, cminmaxrad_l, cminmaxrad_r, cminmaxrad, &ctoprad, cminmaxrad2_l, cminmaxrad2_r, cminmaxrad2, &ctoprad2);

  ///////////////////
  // choose ctopall that will go into choosing timestep
  if(EOMRADTYPE!=EOMRADNONE){
    // get maximum over both mhd and radiation 
    *ctopallptr=MAX(ctopmhd,ctoprad2); // choose ctoprad2 for setting timestep

    if(FORCESOLVEL){
      *ctopallptr=1.0/sqrt(ptrgeom->gcov[GIND(dir,dir)]);
    }

  }
  else *ctopallptr=ctopmhd;



  //  dualfprintf(fail_file,"ctopmhdprim=%g ctopradprim=%g\n",ctopmhd*sqrt(ptrgeom->gcov[GIND(dir,dir)]),ctoprad*sqrt(ptrgeom->gcov[GIND(dir,dir)]));



  // cminmaxrad and ctoprad are used for fluxes (not "2" version that's only for setting timestep)
  MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, U_l, U_r, F_l, F_r, F),"step_ch.c:fluxcalc()", "flux_compute", 1);


#if(FLUXDUMP==2)
  if(dir==1){ // only for radial fluxes
    // no accounting for RK, so only correct for RK2 midpoint method right now
    PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,0*NPR + pl)=F[pl];


    MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, UPAKE_l, UPAKE_r, FPAKE_l, FPAKE_r, FPAKE),"step_ch.c:fluxcalc()", "flux_compute", 1);
    PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,1*NPR + pl)=FPAKE[pl];

    MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, UEN_l, UEN_r, FEN_l, FEN_r, FEN),"step_ch.c:fluxcalc()", "flux_compute", 1);
    PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,2*NPR + pl)=FEN[pl];

    MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, UEM_l, UEM_r, FEM_l, FEM_r, FEM),"step_ch.c:fluxcalc()", "flux_compute", 1);
      PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,3*NPR + pl)=FEM[pl];

    if(EOMRADTYPE!=EOMRADNONE){
      MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, URAD_l, URAD_r, FRAD_l, FRAD_r, FRAD),"step_ch.c:fluxcalc()", "flux_compute", 1);
      PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + pl)=FRAD[pl];
    }
  }
#endif

#if(FLUXDUMP==1)
  diag_fluxdump_1(dir,i,j,k,p_l, p_r, F_l, F_r);
#endif  

  return(0);
}



/// actually compute flux from given data
/// here F is MA only and FEM is EM only
/// only should be needed for old a2c method
int flux_compute_splitmaem(int i, int j, int k, int dir, struct of_geom *ptrgeom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *FEM, FTYPE *ctopallptr)
{
  int flux_compute(int i, int j, int k, int dir, struct of_geom *geom, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE ctopmhd, FTYPE *cminmaxrad_l, FTYPE *cminmaxrad_r, FTYPE *cminmaxrad, FTYPE ctoprad, FTYPE CUf, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, FTYPE *F);
  int p2SFUevolve_splitmaem(int dir, int isleftright, FTYPE *p, struct of_geom *geom, struct of_state **ptrstate, FTYPE *F, FTYPE *FEM, FTYPE *U, FTYPE *UEM);
  FTYPE cminmax_l[NUMCS], cminmax_r[NUMCS], cminmax[NUMCS];
  FTYPE cminmaxrad_l[NUMCS], cminmaxrad_r[NUMCS], cminmaxrad[NUMCS];
  FTYPE cminmaxrad2_l[NUMCS], cminmaxrad2_r[NUMCS], cminmaxrad2[NUMCS];
  FTYPE ctopmhd,ctoprad,ctoprad2;
  struct of_state state_c, state_l, state_r;
  struct of_state *ptrstate_c, *ptrstate_l, *ptrstate_r;
  FTYPE F_c[NPR], F_l[NPR], F_r[NPR];
  FTYPE U_c[NPR], U_l[NPR], U_r[NPR];
  FTYPE FEM_c[NPR], FEM_l[NPR], FEM_r[NPR];
  FTYPE UEM_c[NPR], UEM_l[NPR], UEM_r[NPR];
  int pl,pliter;
  int jj,kk;


  // default
  ptrstate_c = &state_c;
  ptrstate_l = &state_l;
  ptrstate_r = &state_r;


  //////////////////////
  //
  // setup flux calculation based upon interpolated primitives
  //
  //////////////////////
  MYFUN(p2SFUevolve_splitmaem(dir, ISLEFT, p_l, ptrgeom, &ptrstate_l, F_l, FEM_l, U_l, UEM_l),"step_ch.c:fluxcalc()", "p2SFUevolve()", 1);
  MYFUN(p2SFUevolve_splitmaem(dir, ISRIGHT, p_r, ptrgeom, &ptrstate_r, F_r, FEM_r, U_r, UEM_r),"step_ch.c:fluxcalc()", "p2SFUevolve()", 2);


  // usually "always" need cminmax_l cminmax_r and always need ctop
  get_wavespeeds(dir, ptrgeom, p_l, p_r, U_l, U_r, F_l, F_r, ptrstate_l, ptrstate_r, cminmax_l, cminmax_r, cminmax, &ctopmhd, cminmaxrad_l, cminmaxrad_r, cminmaxrad, &ctoprad, cminmaxrad2_l, cminmaxrad2_r, cminmaxrad2, &ctoprad2);


  ///////////////////
  // choose ctopall that will go into choosing timestep
  if(EOMRADTYPE!=EOMRADNONE){
    // get maximum over both mhd and radiation for setting dt later
    *ctopallptr=MAX(ctopmhd,ctoprad2); // choose ctoprad2 for setting timestep
  }
  else *ctopallptr=ctopmhd;


  // cminmaxrad and ctoprad are used for fluxes (not "2" version that's only for setting timestep)

  // GODMARK:
  // below assumes flux_compute is linear in FMA and FEM, which it is right now
  // if used Riemann solver, would need to have it output FMA and FEM parts
  MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, U_l, U_r, F_l, F_r, F),"step_ch.c:fluxcalc()", "flux_compute", 1);

  MYFUN(flux_compute(i, j, k, dir, ptrgeom, cminmax_l,cminmax_r, cminmax, ctopmhd, cminmaxrad_l,cminmaxrad_r, cminmaxrad, ctoprad, CUf, p_l, p_r, UEM_l, UEM_r, FEM_l, FEM_r, FEM),"step_ch.c:fluxcalc()", "flux_compute", 1);


  // Note that if splitmaem==1 then pressure term was moved as separate quasi-flux/conserved term.  flux remains, and conserved quantity not used


#if(FLUXDUMP==1)
  diag_fluxdump_1(dir,i,j,k,p_l, p_r, F_l, F_r);
#endif  


  return(0);
}



/// P->sate,flux,U
int p2SFUevolve(int dir, int isleftright, FTYPE *p, struct of_geom *ptrgeom, struct of_state **ptrstate, FTYPE *F, FTYPE *U)
{
  MYFUN(get_stateforfluxcalc(dir,isleftright, p, ptrgeom, ptrstate),"flux.c:p2SFUevolve", "get_state()", 1);


  // DEBUG:
  //  if(ptrgeom->i==26 && ptrgeom->j==40 && dir==1){
  //    dualfprintf(fail_file,"NORMAL: INp2SFUevolve: gdet=%21.15g : B1=%21.15g B2=%21.15g B3=%21.15g uu0=%21.15g uu1=%21.15g uu2=%21.15g uu3=%21.15g\n",ptrgeom->gdet,p[B1],p[B2],p[B3],(*ptrstate)->ucon[TT],(*ptrstate)->ucon[1],(*ptrstate)->ucon[2],(*ptrstate)->ucon[3]);
  //  }


  MYFUN(primtoflux(UEVOLVE,p, *ptrstate, dir, ptrgeom, F, NULL),"flux.c:p2SFUevolve()","primtoflux_calc() dir=1/2 l", 1);
  MYFUN(primtoflux(UEVOLVE,p, *ptrstate, TT, ptrgeom, U, NULL),"flux.c:p2SFUevolve()", "primtoflux_calc() dir=l0", 1);

  return(0);
}


// p2SFUevolve() for each physical term
int p2SFUevolveall(int dir, int isleftright, FTYPE *p, struct of_geom *ptrgeom, struct of_state **ptrstate, FTYPE *FPAKE, FTYPE *FEN, FTYPE *FEM, FTYPE *FRAD, FTYPE *UPAKE, FTYPE *UEN, FTYPE *UEM, FTYPE *URAD)
{


  // DEBUG:
  //  if(ptrgeom->i==26 && ptrgeom->j==40 && dir==1){
  //    dualfprintf(fail_file,"NORMAL: INp2SFUevolve: gdet=%21.15g : B1=%21.15g B2=%21.15g B3=%21.15g uu0=%21.15g uu1=%21.15g uu2=%21.15g uu3=%21.15g\n",ptrgeom->gdet,p[B1],p[B2],p[B3],(*ptrstate)->ucon[TT],(*ptrstate)->ucon[1],(*ptrstate)->ucon[2],(*ptrstate)->ucon[3]);
  //  }

  MYFUN(get_stateforfluxcalc(dir,isleftright, p, ptrgeom, ptrstate),"flux.c:p2SFUevolve", "get_state()", 1);

  int pliter,pl,jj;
  struct of_state qtemp;
  FTYPE ptemp[NPR];

  PLOOP(pliter,pl) ptemp[pl] = p[pl]; qtemp=**ptrstate;
  ptemp[UU]=qtemp.pressure=qtemp.entropy=SMALL; //ug
  ptemp[B1]=ptemp[B2]=ptemp[B3]=0.0; DLOOPA(jj) qtemp.bcon[jj]=qtemp.bcov[jj]=0.0; qtemp.bsq=0.0; // b
  if(URAD0>=0) ptemp[URAD0]=SMALL; // urad
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, dir, ptrgeom, FPAKE,NULL),"flux.c:p2SFUevolve()","primtoflux_calc() dir=1/2 l", 1);
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, TT, ptrgeom, UPAKE,NULL),"flux.c:p2SFUevolve()", "primtoflux_calc() dir=l0", 1);

  ptemp[RHO]=SMALL; // rho
  ptemp[B1]=ptemp[B2]=ptemp[B3]=0.0; DLOOPA(jj) qtemp.bcon[jj]=qtemp.bcov[jj]=0.0; qtemp.bsq=0.0; // b
  if(URAD0>=0) ptemp[URAD0]=SMALL; // urad
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, dir, ptrgeom, FEN,NULL),"flux.c:p2SFUevolve()","primtoflux_calc() dir=1/2 l", 1);
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, TT, ptrgeom, UEN,NULL),"flux.c:p2SFUevolve()", "primtoflux_calc() dir=l0", 1);

  ptemp[RHO]=SMALL; // rho
  ptemp[UU]=qtemp.pressure=qtemp.entropy=SMALL; //ug
  if(URAD0>=0) ptemp[URAD0]=SMALL; // urad
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, dir, ptrgeom, FEM,NULL),"flux.c:p2SFUevolve()","primtoflux_calc() dir=1/2 l", 1);
  MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, TT, ptrgeom, UEM,NULL),"flux.c:p2SFUevolve()", "primtoflux_calc() dir=l0", 1);

  if(EOMRADTYPE!=EOMRADNONE){
    ptemp[RHO]=SMALL; // rho
    ptemp[UU]=qtemp.pressure=qtemp.entropy=SMALL; //ug
    ptemp[B1]=ptemp[B2]=ptemp[B3]=0.0; DLOOPA(jj) qtemp.bcon[jj]=qtemp.bcov[jj]=0.0; qtemp.bsq=0.0; // b
    MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, dir, ptrgeom, FRAD,NULL),"flux.c:p2SFUevolve()","primtoflux_calc() dir=1/2 l", 1);
    MYFUN(primtoflux(UEVOLVE,ptemp, &qtemp, TT, ptrgeom, URAD,NULL),"flux.c:p2SFUevolve()", "primtoflux_calc() dir=l0", 1);
  }


  return(0);
}


/// P->sate,flux,U for splitmaem method
int p2SFUevolve_splitmaem(int dir, int isleftright, FTYPE *p, struct of_geom *ptrgeom, struct of_state **ptrstate, FTYPE *F, FTYPE *FEM, FTYPE *U, FTYPE *UEM)
{

#if(ROEAVERAGEDWAVESPEED && ((USESTOREDSPEEDSFORFLUX==0)||(STOREWAVESPEEDS==0)) )
  dualfprintf(fail_file,"Only supposed to use this splitting of F and U pulling out particular linear terms if not assuming U is not true conserved quantity\n");
  myexit(72679262);
#endif

  MYFUN(get_stateforfluxcalc(dir,isleftright, p, ptrgeom, ptrstate),"flux.c:p2SFUevolve", "get_state()", 1);

  MYFUN(primtoflux_splitmaem(UEVOLVE, p, *ptrstate, dir, dir, ptrgeom, F, FEM),"flux.c:p2SFUevolve_splitmaem()","primtoflux_ma() dir=1/2 l", 1);
  MYFUN(primtoflux_splitmaem(UEVOLVE, p, *ptrstate, dir, TT, ptrgeom, U, UEM),"flux.c:p2SFUevolve_splitmaem()", "primtoflux_ma() dir=l0", 1);




  return(0);
}




/// DEBUG flux store
void diag_fluxdump_1(int dir, int i, int j, int k, FTYPE *p_l, FTYPE *p_r, FTYPE *F_l, FTYPE *F_r)
{
  int pl,pliter;

#if(FLUXDUMP==1)
  PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*1 + pl)=F_l[pl];
  PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*2 + pl)=F_r[pl];
  PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*3 + pl)=p_l[pl];
  PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*4 + pl)=p_r[pl];

  if(0 && dir==2 && j==0 && i==N1/2){ // just pick one point on polar axis
    PLOOP(pliter,pl){
      //      dualfprintf(fail_file,"n=%ld sp=%d pl=%d :: F_l=%21.15g F_r=%21.15g U_l=%21.15g U_r=%21.15g p_l=%21.15g p_r=%21.15g\n",nstep,steppart,pl,F_l[pl],F_r[pl],U_l[pl],U_r[pl],p_l[pl],p_r[pl]);
    }
  }
#endif
}






/// Harten-Lax-van Leer
#define HLLCOMPUTE(cmin,cmax,U_l,U_r,F_l,F_r) ((  (cmax * F_l + cmin * F_r) - (FLUXDISSIPATION)*(cmax * cmin) * (U_r - U_l) ) / (cmax + cmin) )

/// Lax-Friedrichs with Rusanov wave speed
#define LAXFCOMPUTE(ctop,U_l,U_r,F_l,F_r) (0.5 * ( (F_l + F_r) - (FLUXDISSIPATION)*ctop * (U_r - U_l) ) )
//#define LAXFCOMPUTE(ctop,U_l,U_r,F_l,F_r) (0.5 * ( (F_l + F_r) ) )
//#define LAXFCOMPUTE(ctop,U_l,U_r,F_l,F_r) (0.5 * ( (F_l + F_r) - (FLUXDISSIPATION)*fabs(ctop) * fabs(U_r - U_l) ) )
//#define LAXFCOMPUTE(ctop,U_l,U_r,F_l,F_r) (0.5 * ( (F_l + F_r) - (FLUXDISSIPATION)*1E-5 * (U_r - U_l) ) )

/// Lax-Friedrichs flux
#define LFCOMPUTE(crus,U_l,U_r,F_l,F_r) (0.5*( (F_l+F_r) - (FLUXDISSIPATION)*crus*(U_r - U_l)))

/// mid-step for 2-step LW-flux
#define UHALF(crus,U_l, U_r, F_l, F_r) (0.5*((U_l+U_r)-(FLUXDISSIPATION)*(F_r-F_l)/crus))

/// 2-step Lax-Wendroff flux
#define LWCOMPUTE(crus,U_l,U_r,F_l,Fhalf,F_r) (Fhalf)

/// FORCE
#define FORCECOMPUTE(crus,U_l,U_r,F_l,Fhalf,F_r) (0.25*( (F_l+F_r) + 2.0*Fhalf - (FLUXDISSIPATION)*crus*(U_r - U_l)))

/// GFORCE // Toro & Titarev JCP 2006
#define GFORCECOMPUTE(cg,crus,U_l,U_r,F_l,Fhalf,F_r) ( (LWCOMPUTE(crus,U_l,U_r,F_l,Fhalf,F_r) + cg*LFCOMPUTE(crus,U_l,U_r,F_l,F_r))/(1.0+cg) )


/// use more HLL as f_s from 0 to 1, use more LAXF as f_s from 1 to 0
#define HLLLAXF1(cmin,cmax,ctop,f_s, U_l, U_r, F_l, F_r)  ( HLLCOMPUTE(cmin,cmax,U_l,U_r,F_l,F_r)*f_s + LAXFCOMPUTE(ctop,U_l,U_r,F_l,F_r)*(1.0-f_s) )





/// actually compute flux from given data
/// GODMARK: right now if splitmaem==1 assumes F is linear in FMA and FEM. In general Riemann solver should return FMA and FEM in parts
int flux_compute(int i, int j, int k, int dir, struct of_geom *geom, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE ctopmhd, FTYPE *cminmaxrad_l, FTYPE *cminmaxrad_r, FTYPE *cminmaxrad, FTYPE ctoprad, FTYPE CUf, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, FTYPE *F)
{
  int pl,pliter;
  void choose_flux(int fluxmethodlocal, int i, int j, int k, int pl, FTYPE *laxffrac,FTYPE *hllfrac);
  FTYPE laxffrac[NPR],hllfrac[NPR];
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  int forceflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int mustaflux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  FTYPE crus[NPR];
  FTYPE cforce[NPR];
  FTYPE cmin_l[NPR], cmax_l[NPR], cmin_r[NPR], cmax_r[NPR], cmax[NPR], cmin[NPR];
  FTYPE ctop[NPR];
  FTYPE dPoP,f_s;

  // defaults
  int fluxmethodlocal[NPR];
  PALLLOOP(pl) fluxmethodlocal[pl]=fluxmethod[pl];
  FTYPE cmaxfactor=1.0;
  FTYPE cminfactor=1.0;
#if(OUTERRADIALSUPERFAST)
  // force effective superfast condition on outer radial boundaries
  if(ISBLACKHOLEMCOORD(MCOORD)){
    if(dir==1 && startpos[1]+i==totalsize[1]){ PALLLOOP(pl) fluxmethodlocal[pl]=HLLFLUX; cminfactor=0.0; }
    if(dir==1 && startpos[1]+i==0){ PALLLOOP(pl) fluxmethodlocal[pl]=HLLFLUX; cmaxfactor=0.0; }
  }
#endif

  // assign cmin/cmax/ctop for each conserved quantity
  PLOOP(pliter,pl){
    if(RADFULLPL(pl)==0){
      // MHD terms
      cmin_l[pl] = cminmax_l[CMIN];
      cmax_l[pl] = cminmax_l[CMAX];
      cmin_r[pl] = cminmax_r[CMIN];
      cmax_r[pl] = cminmax_r[CMAX];
      cmin[pl] = cminmax[CMIN]*cminfactor;
      cmax[pl] = cminmax[CMAX]*cmaxfactor;
      ctop[pl] = ctopmhd;

    }
    else{
      // radiation terms
      cmin_l[pl] = cminmaxrad_l[CMIN];
      cmax_l[pl] = cminmaxrad_l[CMAX];
      cmin_r[pl] = cminmaxrad_r[CMIN];
      cmax_r[pl] = cminmaxrad_r[CMAX];
      cmin[pl] = cminmaxrad[CMIN]*cminfactor;
      cmax[pl] = cminmaxrad[CMAX]*cmaxfactor;
      ctop[pl] = ctoprad;
    }
    // if(steppart==0) crus[pl]=dx[dir]/(dt*0.5);
    // else  crus[pl]=dx[dir]/(dt);
    crus[pl]=dx[dir]/(dt*CUf);
    // crus[pl]=dx[dir]/(dt);
    // dualfprintf(fail_file,"CUf=%21.15g dt=%21.15g dx[%d]=%21.15g\n",CUf,dt,dir,dx[dir]);
    // crus[pl]=dx[dir]/(dt);
  }


  //////////////////
  //
  // store non-dissipative contribution
  //
  ////////////////////
  int plsp;
  if(NSPECIAL>=1){
    plsp=NPR;
    F[plsp] = 0.5 * (F_l[SPECIALPL1] + F_r[SPECIALPL1]);
  }
  if(NSPECIAL>=2){
    plsp=NPR+1;
    F[plsp] = 0.5 * (F_l[SPECIALPL2] + F_r[SPECIALPL2]);
  }
  if(NSPECIAL>=3){
    plsp=NPR+2;
    F[plsp] = 0.5 * (F_l[SPECIALPL3] + F_r[SPECIALPL3]);
  }
  if(NSPECIAL>=4){
    plsp=NPR+3;
    F[plsp] = 0.5 * (F_l[SPECIALPL4] + F_r[SPECIALPL4]);
  }
  if(NSPECIAL>=5){
    plsp=NPR+4;
    F[plsp] = 0.5 * (F_l[SPECIALPL5] + F_r[SPECIALPL5]);
  }
  if(NSPECIAL>=6){
    plsp=NPR+5;
    F[plsp] = 0.5 * (F_l[SPECIALPL6] + F_r[SPECIALPL6]);
  }


  //////////////////
  //
  // store full non-dissipative + dissipative contributions
  //
  ////////////////////
  if(fluxmethodlocal[RHO]==LAXFFLUX || fluxmethodlocal[RHO]==HLLFLUX){
    PLOOP(pliter,pl){
      if(fluxmethodlocal[pl]==LAXFFLUX) F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
      else hllflux_compute(dir,geom,cmin,cmax,ctop,p_l,p_r,U_l,U_r,F_l,F_r,F);
    }     
  }
  else if(fluxmethodlocal[RHO]==HLLFLUX){
    //////////////////////////////////
    //
    // decide which flux formula to use
    //
    /////////////////////////////////
    //    PLOOP(pliter,pl) choose_flux(fluxmethodlocal, i,j,k,pl,laxffrac,hllfrac);


    //////////////////////////////////
    //
    // Decide if using different flux calculation on boundary zones
    //
    /////////////////////////////////
    //#if(HLLBOUNDARY)
    //    if(
    //       ((dir == 3) && ( ((startpos[3]+k == 0)&&(BCtype[X3DN]==OUTFLOW)) || ((startpos[3]+k == totalsize[3])&&(BCtype[X3UP]==OUTFLOW))  ))||
    //       ((dir == 2) && ( ((startpos[2]+j == 0)&&(BCtype[X2DN]==POLARAXIS)) || ((startpos[2]+j == totalsize[2])&&(BCtype[X2UP]==POLARAXIS))  ))||
    //       ((dir == 1) && ( ((startpos[1]+i == 0)&&(BCtype[X1DN]==OUTFLOW)) || ((startpos[1]+i == totalsize[1])&&(BCtype[X1UP]==OUTFLOW))  ))
    //       )
    //      {
    // PLOOP(pliter,pl){
    //   hllfrac[pl]=1.0;
    //   laxffrac[pl]=0.0;
    // }
    //      }
    //#endif

    hllflux_compute(dir,geom,cmin,cmax,ctop,p_l,p_r,U_l,U_r,F_l,F_r,F);
    
  }
  else if(fluxmethodlocal[RHO]==FORCEFLUX){

    // normal Rusanov speed
    //cforce=crus;
    // LAXF speed
    PLOOP(pliter,pl) cforce[pl]=ctop[pl];


    forceflux_compute(dir,geom,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);

  }
  else if(fluxmethodlocal[RHO]==MUSTAFLUX){

    PLOOP(pliter,pl) {
      //cforce[pl]=crus[pl];
      cforce[pl]=ctop[pl];
    }
    mustaflux_compute(dir,geom,cmin_l,cmin_r,cmax_l,cmax_r,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);

  }
  else if(fluxmethodlocal[RHO]==HLLLAXF1FLUX){

#define MINDPOP (0.0)
#define MAXDPOP (0.3)

    dPoP = fabs( (p_r[UU]-p_l[UU])/(0.5*(p_r[UU]+p_l[UU])) );
    f_s = (1.0-0.0)/(MAXDPOP-MINDPOP) * (dPoP - MINDPOP) + 0.0;
    if(f_s>1.0) f_s=1.0;
    if(f_s<0.0) f_s=0.0;

    PLOOP(pliter,pl){
      if(cmax[pl]+cmin[pl]!=0.0){
        F[pl] = HLLLAXF1(cmin[pl],cmax[pl],ctop[pl],f_s, U_l[pl], U_r[pl], F_l[pl], F_r[pl]);
      }
      else{
        F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
      }
    }
  
    
    
    //PLOOP(pliter,pl) dualfprintf(fail_file,"p_l[%d]=%g p_r[%d]=%g U_l[%d]=%g U_r[%d]=%g F_l[%d]=%g F_r[%d]=%g :: F[%d]=%g\n",pl,p_l[pl],pl,p_r[pl],pl,U_l[pl],pl,U_r[pl],pl,F_l[pl],pl,F_r[pl],pl,F[pl]);
    
    //    dualfprintf(fail_file,"cmin=%g cmax=%g ctop=%g\n",cmin,cmax,ctop);

    //    dualfprintf(fail_file,"i=%d dPoP=%g f_s=%g :: %g %g\n",geom->i,dPoP,f_s,p_l[UU],p_r[UU]);
    

  }









  return(0);
}



/// e.g. double rarefaction Einfelt test (test2 in Liska & Wendroff 2003) leads to F at interface corresponding to negative pressure.  This is supposed to correct for that.  Leads to slightly better solution except for near center
#define USE_CORRECTED_STATES 0

/// compute HLL flux
int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F)
{
  int pl,pliter;
  FTYPE vmin[NPR],vmax[NPR];
  FTYPE cminreal[NPR],cmaxreal[NPR];

  PLOOP(pliter,pl) {


#if(USE_CORRECTED_STATES)
    // get vmin/vmax
    vmin[pl]=min(p_l[UU+dir],p_r[UU+dir]);
    vmax[pl]=max(p_l[UU+dir],p_r[UU+dir]);
    vmin[pl] = fabs(max(0., -vmin[pl]));
    vmax[pl] = fabs(max(0., vmax[pl]));

    // get HLL+correction  
    if(cmax[pl]+cmin[pl]!=0.0){
      if( pl == UU + dir && (vmin[pl]+vmax[pl]!=0.0) ) {
        //      if(vmin[pl]+vmax[pl]!=0.0){
        cminreal[pl] = vmin[pl];
        cmaxreal[pl] = vmax[pl];
      }
      else {
        cminreal[pl] = cmin[pl];
        cmaxreal[pl] = cmax[pl];
      }
      F[pl] =  HLLCOMPUTE(cminreal[pl],cmaxreal[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
    }
    else F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);

#else

    if(cmax[pl]+cmin[pl]!=0.0){
      F[pl] =  HLLCOMPUTE(cmin[pl],cmax[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
    }
    else F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);

#endif
  
  } // over pl's
  
  return(0);

}


#define MUSTACOEF (1.0)
/// exactly 1.0 works good with MUSTAHLL to give stationary contact
/// works good with MUSTAHLL

//#define MUSTACOEF (1.0)
// HLL+0.9 gives nice and smooth Noh, but blurs stationary contact

//#define MUSTACOEF (0.49)
// HLL+0.49 seems no better than without MUSTA

//#define MUSTACOEF (0.4) // ok test2 velocity with MUSTAHLL

//#define MUSTACOEF (0.9) // ok test2 velocity with MUSTAHLL

//#define MUSTACOEF (0.96) // ok test2 velocity with MUSTAHLL

//#define MUSTACOEF (0.99) // death with MUSTAHLL

//#define MUSTACOEF (0.999) // death with MUSTAHLL

//#define MUSTACOEF (0.9990)
// 0.9 gives nice and smooth Noh, but blurs stationary contact
// 0.99 gives oscillatory Noh and blurs stationary a bit
// 0.9990 gives good double ratefaction as long as a2c is turned off, and for a2c on or off a nearly stationary contact and decent Noh


//#define MUSTACOEF (2.0)
// MUSTAFLUX==MUSTAFORCE w/ MUSTACOEF==2.0 has stationary contact
// does ok with rarefaction, but fails a bit like MUSTAHLL w/ MUSTACOEF=1.0
// bit more oscillatory in flat region for Noh

/// force flux
int forceflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F)
{
  int pl,pliter;
  FTYPE umid[NPR],pmid[NPR],fmid[NPR];
  struct of_state state;
  struct of_state *ptrstate;
  int doforceflux;
  FTYPE vmin[NPR],vmax[NPR];
  FTYPE vminorig[NPR], vmaxorig[NPR];
  FTYPE cminreal[NPR], cmaxreal[NPR];
  int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F); 
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; 

  ////////////////////////////
  //
  // get middle umid
  //
  ////////////////////////////
  PLOOP(pliter,pl) umid[pl] = UHALF(cforce[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"%d : cforce=%21.15g : U_l[%d]=%21.15g U_r[%d]=%21.15g umid[%d]=%21.15g\n",geom->i,cforce,pl,U_l[pl],pl,U_r[pl],pl,umid[pl]);
  
  // set guess for inversion
  PLOOP(pliter,pl) pmid[pl]=0.5*(p_l[pl]+p_r[pl]);
  // get primitive pmid(umid)
  int eomtype=EOMDEFAULT;
  FTYPE dissmeasure=-1.0; // assume energy try ok
  int whichcap=CAPTYPEBASIC;
  int whichmethod=MODEDEFAULT;
  int modprim=0;
  int checkoninversiongas=CHECKONINVERSION;
  int checkoninversionrad=CHECKONINVERSIONRAD;
  MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, umid, NULL, geom, dissmeasure, pmid, pmid,&newtonstats),"flux.c:flux_compute()", "Utoprimgen", 1);
  doforceflux=1;
  if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
    if(debugfail>=1) dualfprintf(fail_file,"Failed to find inversion for FORCEFLUX, trying p_l : nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
    PLOOP(pliter,pl) pmid[pl]=p_l[pl];
    // get primitive pmid(umid)
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, umid, NULL, geom, dissmeasure, pmid, pmid,&newtonstats),"flux.c:flux_compute()", "Utoprimgen", 1);
    if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
      if(debugfail>=1) dualfprintf(fail_file,"Failed to find inversion for FORCEFLUX, trying p_r : nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
      PLOOP(pliter,pl) pmid[pl]=p_r[pl];
      // get primitive pmid(umid)
      MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, umid, NULL, geom, dissmeasure, pmid, pmid,&newtonstats),"flux.c:flux_compute()", "Utoprimgen", 1);
      if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
        if(debugfail>=1) dualfprintf(fail_file,"No initial guess worked, rejecting FORCEFLUX method : nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
        doforceflux=0;
      }
    }
  }

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"%d : pmid[%d]=%21.15g p_l[%d]=%21.15g p_r[%d]=%21.15g\n",geom->i,pl,pmid[pl],pl,p_l[pl],pl,p_r[pl]);

  /////////////////////////////////////
  //
  // get flux for pmid Fmid(pmid)
  //
  ////////////////////////////////////
  if(doforceflux){
    ptrstate=&state; // default
    MYFUN(get_stateforfluxcalc(dir,ISMIDDLE, pmid, geom, &ptrstate),"flux.c:flux_compute()", "get_state()", 1);
    MYFUN(primtoflux(UEVOLVE, pmid, ptrstate, dir, geom, fmid, NULL),"flux.c:flux_compute()","primtoflux_calc() dir=1/2 l", 1);
    
    // compute FORCE flux
    //    PLOOP(pliter,pl) F[pl] = FORCECOMPUTE(cforce,U_l[pl],U_r[pl],F_l[pl],fmid[pl],F_r[pl]);
    PLOOP(pliter,pl) F[pl] = GFORCECOMPUTE(MUSTACOEF,cforce[pl],U_l[pl],U_r[pl],F_l[pl],fmid[pl],F_r[pl]);

    //    PLOOP(pliter,pl) dualfprintf(fail_file,"%d : F[%d]=%21.15g Fhll[%d]=%21.15g\n",geom->i,pl,F[pl],pl,HLLCOMPUTE(cmin,cmax,U_l[pl],U_r[pl],F_l[pl],F_r[pl]));
  }
  else{
    // reduce to HLL
    hllflux_compute(dir,geom,cmin,cmax,ctop,p_l,p_r,U_l,U_r,F_l,F_r,F);
  }

  return(0);
}



/// number of single-cell musta iterations
#define NUMMUSTAITERS 1

/// whether to do mutli-cell version of MUSTA
#define DOMULTICELL 0
/// number of multicell musta iterations
#define NUMMULTIMUSTAITERS 2 // 1=expensive way to do no musta, 2 minimum for interesting calculation
/// NUMMULTIMUSTAITERS=2 and NUMLOCALCELLS=1 is simplest scheme
/// multicell much much slower than single cell version
/// number of musta cells
#define NUMLOCALCELLS 2
#define MULTIMUSTACOEF (0.9)


#define MUSTAFORCE 0
#define MUSTALAXF 1
#define MUSTAHLL 2

#define WHICHFLUX MUSTAHLL
//#define WHICHFLUX MUSTALAXF
//#define WHICHFLUX MUSTAFORCE
// MUSTAFORCE worse than MUSTAHLL. Contacts are blurred (at one point this wasn't true!)
// MUSTAFORCE doesn't fail on rarefaction however, while MUSTAHLL does a bit
// MUSTAFORCE behaves better with a2c in FV method


// whether to limit MUSTA by some fraction of HLL
#define HLLBOUNDMUSTA 0

/// musta flux
int mustaflux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F)
{
  int musta1flux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int musta2flux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int domustaflux;
  FTYPE cmaxorig[NPR],cminorig[NPR],ctoporig[NPR];
  FTYPE Uorig_l[NPR],Uorig_r[NPR],Forig_l[NPR],Forig_r[NPR],porig_l[NPR],porig_r[NPR];
  FTYPE Forig[NPR];
  FTYPE Fother[NPR];
  int pl,pliter;
  FTYPE shockstrength;
  FTYPE rarestrength;
  FTYPE contactstrength;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; 



  if(DOEVOLVEUU){
    // check for strong shock that MUSTA can't handle
    shockstrength=fabs((p_r[UU]-p_l[UU])/max(min(fabs(p_l[UU]),fabs(p_r[UU])),SMALL));
    
    // contact strength
    contactstrength=fabs((p_r[RHO]-p_l[RHO])/max(min(fabs(p_l[RHO]),fabs(p_r[RHO])),SMALL));
  }
  else{
    shockstrength=contactstrength=0;
  }

  // rarefaction strength
  rarestrength=fabs((p_r[UU+dir]-p_l[UU+dir])/max(min(fabs(p_l[UU+dir]),fabs(p_r[UU+dir])),SMALL));


  // keep original left/right states in case of musta failure
  PLOOP(pliter,pl){

    cmaxorig[pl]=cmax[pl];
    cminorig[pl]=cmin[pl];
    ctoporig[pl]=ctop[pl];

    Uorig_l[pl]=U_l[pl];
    Uorig_r[pl]=U_r[pl];
    Forig_l[pl]=F_l[pl];
    Forig_r[pl]=F_r[pl];
    porig_l[pl]=p_l[pl];
    porig_r[pl]=p_r[pl];
  }


  
  if(DOMULTICELL==0){
    domustaflux=musta1flux_compute(dir,geom,cmin_l,cmin_r,cmax_l,cmax_r,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);
  }
  else if(DOMULTICELL==1){
    domustaflux=musta2flux_compute(dir,geom,cmin_l,cmin_r,cmax_l,cmax_r,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);
  }

  /////////////////
  //
  //  DONE WITH MUSTA.  F is final musta solution or musta failed
  //
  ////////////////

  // check if musta failed during inversion
  if(!domustaflux){
    // reduce to HLL
    hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,F);

    //    if(debugfail>=1) dualfprintf(fail_file,"DIDNOTDO_MUSTAFLUX: nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
  }
  else{
    //    if(debugfail>=1) dualfprintf(fail_file,"DIDDO_MUSTAFLUX: nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
  }
    
#if(0)
  //  PLOOP(pliter,pl) F[pl] = HLLCOMPUTE(cmin,cmax,U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
  hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,Fother);
  // average HLL and MUSTAHLL
  PLOOP(pliter,pl) F[pl] = 0.5*(Fother[pl]+F[pl]);
#endif



#if(HLLBOUNDMUSTA)
  // don't allow for a change in the flux by more than 200%
  hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,Fother);
  //  PLOOP(pliter,pl) if( (Fother[pl]-F[pl])/(fabs(Fother[pl])+SMALL)>0.9) F[pl]=Fother[pl];
  PLOOP(pliter,pl) if( fabs(Fother[pl]-F[pl])/(fabs(Fother[pl])+SMALL)>2.0) F[pl]=Fother[pl];

#endif



  // check sign of flux and don't allow to switch
#if(0)
  hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,Fother);
  PLOOP(pliter,pl) if( ((Fother[pl]<0.0)&&(F[pl]>0.0))||((Fother[pl]>0.0)&&(F[pl]<0.0)) ) F[pl]=Fother[pl];
#endif



#if(0)
  hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,Fother);
  //F[RHO]=Fother[RHO];
  if( (shockstrength>0.5)||(rarestrength>0.5)){
    F[UU]=Fother[UU];
    F[U1]=Fother[U1];// leads to inconsistent rho*v with F[RHO]
    F[U2]=Fother[U2];
    F[U3]=Fother[U3];
  }
  //  }
  //  F[B1]=Fother[B1];
  //  F[B2]=Fother[B2];
  //  F[B3]=Fother[B3];


#endif

#if(1) // not too bad, but generates problem at moving contact.
  hllflux_compute(dir,geom,cminorig,cmaxorig,ctoporig,porig_l,porig_r,Uorig_l,Uorig_r,Forig_l,Forig_r,Fother);

  F[UU]=Fother[UU];
  F[U1]=Fother[U1];// leads to inconsistent rho*v with F[RHO]
  F[U2]=Fother[U2];
  F[U3]=Fother[U3];
#endif

 
  return(0);
}



/// first version of MUSTA flux
int musta1flux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F)
{
  int forceflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int pl,pliter;
  FTYPE umid[NPR],plnew[NPR],prnew[NPR],fmid[NPR];
  struct of_state state, state_l, state_r;
  struct of_state *ptrstate, *ptrstate_l, *ptrstate_r;
  int domustaflux;
  int mustaloop;
  FTYPE cmusta[NPR];
  FTYPE cmineach_l[NUMEOMSETS],cmaxeach_l[NUMEOMSETS];
  FTYPE cmineach_r[NUMEOMSETS],cmaxeach_r[NUMEOMSETS];
  int ignorecourant;
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  FTYPE Fother[NPR];
  FTYPE correctionl[NPR],correctionr[NPR],mymustacoef[NPR],mymustacoeffinal;
  //  FTYPE correctionl,correctionr,mymustacoef;
  FTYPE Ucl,Ucr;
  int pllargest;
  FTYPE fracl[NPR],fracr[NPR];
  FTYPE mustaf;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  FTYPE cminmhd,cmaxmhd,ctopmhd;
  FTYPE cminrad,cmaxrad,ctoprad;
  FTYPE cminrad2,cmaxrad2,ctoprad2;
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; 


#if(EOMRADTYPE!=EOMRADNONE)
  dualfprintf(fail_file,"musta1flux_compute not setup for radiation.");
  myexit(1);
#endif

  // default
  ptrstate = &state;
  ptrstate_l = &state_l;
  ptrstate_r = &state_r;


  // default is to do musta unless musta fails
  domustaflux=1;


  // get predictor flux
#if(WHICHFLUX==MUSTAFORCE)
  PLOOP(pliter,pl) cmusta[pl]=ctop[pl];
  forceflux_compute(dir,geom,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);
#elif(WHICHFLUX==MUSTAHLL)
  PLOOP(pliter,pl) cmusta[pl]=ctop[pl];
  hllflux_compute(dir,geom,cmin,cmax,ctop,p_l,p_r,U_l,U_r,F_l,F_r,F);
#elif(WHICHFLUX==MUSTALAXF)
  PLOOP(pliter,pl) cmusta[pl]=ctop[pl];
  PLOOP(pliter,pl) F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
#endif



  // do musta stage loop
  for(mustaloop=0;mustaloop<NUMMUSTAITERS;mustaloop++){

    PLOOP(pliter,pl){
      // get wave speeds for opening Riemann fan
#if(WHICHFLUX==MUSTAFORCE)
      //      cmusta[pl]=cforce[pl];
      cmusta[pl]=ctop[pl];
#elif(WHICHFLUX==MUSTAHLL)
      cmusta[pl]=ctop[pl];
#elif(WHICHFLUX==MUSTALAXF)
      cmusta[pl]=ctop[pl];
#endif
    }


    ////////////////////////////////////////////
    //
    //
#define MUSTAVERSION 0

    // update U_l U_r (Open Riemann fan)
#if(MUSTAVERSION==0)
    // does bad with Noh problem
    PLOOP(pliter,pl){
      U_l[pl] -= MUSTACOEF*(F[pl]-F_l[pl])/cmusta[pl];
      U_r[pl] -= MUSTACOEF*(F_r[pl]-F[pl])/cmusta[pl];
    }
#elif(MUSTAVERSION==1)
    // does good with Noh problem

    PLOOP(pliter,pl){
      correctionl[pl]=MUSTACOEF*(F[pl]-F_l[pl])/cmusta[pl];
      correctionr[pl]=MUSTACOEF*(F_r[pl]-F[pl])/cmusta[pl];
      // check for monotonicity
      //      if(fabs(U_l[pl]-U_r[pl])<fabs((U_l[pl]-correctionl[pl])-(U_r[pl]-correctionr[pl]))) return(0);
      fracl[pl]=fabs(U_l[pl]-correctionl[pl])/(fabs(U_l[pl])+SMALL);
      fracr[pl]=fabs(U_r[pl]-correctionr[pl])/(fabs(U_r[pl])+SMALL);

      mymustacoef[pl]=(U_l[pl]-U_r[pl])/(correctionl[pl]-correctionr[pl]);
      // check that corrections are possibly tunable to be monotonic and so are at least in the right direction
    }
    
    // find maximum of coefficients
    mymustacoeffinal=0.0;
    PLOOP(pliter,pl){
      if(pl==RHO){
        if(mymustacoeffinal<=mymustacoef[pl]){
          mymustacoeffinal=mymustacoef[pl];
        }
      }
    }
    if(mymustacoeffinal>1.0) mymustacoeffinal=1.0;
    if(mymustacoeffinal<0.9) mymustacoeffinal=0.9;
    
    // find maximum fractional change
    //    mustaf=0.0;
    //    pllargest=-1;
    //    PLOOP(pliter,pl){
    //      if(pl<=RHO){
    // if(mustaf<=max(fracl[pl],fracr[pl])){
    //   mustaf=max(fracl[pl],fracr[pl]);
    //   pllargest=pl;
    // }
    //      }
    //    }
    //    if((mustaf>1.0)||(mustaf<0.0)) return(0);
    // MUSTA is now correctable


    // now use single musta coefficient
    PLOOP(pliter,pl){
      //      mymustacoef=(mymustacoef>1.0) ? 1.0 : ((mymustacoef<0.0) ? 0.0 : mymustacoef);
      //      mymustacoef=MUSTACOEF;
      U_l[pl] -= mymustacoeffinal*correctionl[pl];
      U_r[pl] -= mymustacoeffinal*correctionr[pl];
    }
#elif(MUSTAVERSION==2)
    PLOOP(pliter,pl){
      correctionl[pl]=MUSTACOEF*(F[pl]-F_l[pl])/cmusta[pl];
      correctionr[pl]=MUSTACOEF*(F_r[pl]-F[pl])/cmusta[pl];
      // check for monotonicity
      //      if(fabs(U_l[pl]-U_r[pl])<fabs((U_l[pl]-correctionl[pl])-(U_r[pl]-correctionr[pl]))) return(0);

      if(fabs(correctionl[pl]-correctionr[pl])!=0.0){
        mymustacoef[pl]=(U_l[pl]-U_r[pl])/(correctionl[pl]-correctionr[pl]);
      }
      else mymustacoef[pl]=BIG;
      // check that corrections are possibly tunable to be monotonic and so are at least in the right direction
      //      if((mymustacoef[pl]>1.0)||(mymustacoef[pl]<0.0)) return(0);

      // if coef<0, then when using coef>0, this one will not be monotone
      //      if(mymustacoef[pl]<0.0) return(0);
      //      if(mymustacoef[pl]<-1.0) return(0);
      if(mymustacoef[pl]<0.0) return(0);
    }
    // MUSTA is now correctable

#define MUSTACOEFTYPE 1

#if(MUSTACOEFTYPE==0)
    if(1){
      PLOOP(pliter,pl){
        dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d coef[%d]=%21.15g\n",nstep,steppart,geom->i,pl,mymustacoef[pl]);
      }
    }
    if(1){
      // find minimum of coefficients
      mymustacoeffinal=BIG;
      PLOOP(pliter,pl){
        if(mymustacoeffinal>mymustacoef[pl]) mymustacoeffinal=mymustacoef[pl];
      }
      if(mymustacoeffinal>1.0) return(0);
      //      if(fabs(mymustacoef[RHO]-1.0)>2.0) return(0);
      //      if(mymustacoeffinal>1.0) mymustacoeffinal=1.0;
      PLOOP(pliter,pl){
        mymustacoef[pl]=mymustacoeffinal;
      }
    }
#elif(MUSTACOEFTYPE==1)
    // find maximum below 1
    mymustacoeffinal=0;
    PLOOP(pliter,pl){
      if( (mymustacoeffinal<mymustacoef[pl])&&(mymustacoef[pl]<=1.0)&&(mymustacoef[pl]>0.9) ) mymustacoeffinal=mymustacoef[pl];
    }
    //    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d coef=%21.15g\n",nstep,steppart,geom->i,mymustacoeffinal);
    //    dualfprintf(fail_file,"p_l[rho]=%21.15g p_r[rho]=%21.15g\n",p_l[RHO],p_r[RHO]);
    //    if(mymustacoeffinal>1.0) return(0);
    //    if(mymustacoeffinal>1.0) mymustacoeffinal=1.0;
    PLOOP(pliter,pl){
      mymustacoef[pl]=mymustacoeffinal;
    }
#elif(MUSTACOEFTYPE==2)
    // find good coefficients (this doesn't make much sense)
    PLOOP(pliter,pl){
      if(mymustacoef[pl]<0.0) mymustacoef[pl]=0.0;
      if(mymustacoef[pl]>1.0) mymustacoef[pl]=1.0;
    }
#elif(MUSTACOEFTYPE==3)
    // find good coefficients
    PLOOP(pliter,pl){
      mymustacoef[pl]=MUSTACOEF;
    }
#endif

    
    //    if(mymustacoeffinal<0.9) mymustacoeffinal=0.9;
    //    mymustacoeffinal=1.0;

    // now use single musta coefficient
    PLOOP(pliter,pl){
      //      mymustacoef=(mymustacoef>1.0) ? 1.0 : ((mymustacoef<0.0) ? 0.0 : mymustacoef);
      //      mymustacoef=MUSTACOEF;

      // get corrected value
      Ucl=U_l[pl]-mymustacoef[pl]*correctionl[pl];
      Ucr=U_r[pl]-mymustacoef[pl]*correctionr[pl];

      // check for monotone result, can use final or pl coefficient
      // checking if final point is between original values
      if((U_r[pl] - U_l[pl])*(Ucl - U_l[pl])<0.0) return(0);
      if((U_r[pl] - U_l[pl])*(Ucr - U_r[pl])>0.0) return(0);

      // then good
      U_l[pl] = Ucl;
      U_r[pl] = Ucr;
    }
#elif(MUSTAVERSION==3)
    PLOOP(pliter,pl){
      correctionl=(F[pl]-F_l[pl])/cmusta[pl];
      correctionr=(F_r[pl]-F[pl])/cmusta[pl];
      mymustacoef=MUSTACOEF;

      // fail      
      U_l[pl] -= mymustacoef*correctionl;
      U_r[pl] -= mymustacoef*correctionr;

      //      U_l[pl] -= MUSTACOEF*correctionl;
      //      U_r[pl] -= MUSTACOEF*correctionr;

      // succeed
      //U_l[pl] -= MUSTACOEF*(F[pl]-F_l[pl])/cmusta;
      //U_r[pl] -= MUSTACOEF*(F_r[pl]-F[pl])/cmusta;

      //      U_l[pl] -= mymustacoef*(F[pl]-F_l[pl])/cmusta[pl];
      //      U_r[pl] -= mymustacoef*(F_r[pl]-F[pl])/cmusta[pl];

    }
#endif
    //
    //
    ////////////////////////////////////////////




    // invert to get p_l p_r so can get F_l F_r and U_l U_r
    // get new primitive p_l
    int eomtype=EOMDEFAULT;
    FTYPE dissmeasure=-1.0; // assume energy try ok
    int whichcap=CAPTYPEBASIC;
    int whichmethod=MODEDEFAULT; 
    int modprim=0;
    int checkoninversiongas=CHECKONINVERSION;
    int checkoninversionrad=CHECKONINVERSIONRAD;
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, U_l, ptrstate_l, geom, dissmeasure, p_l, p_l,&newtonstats),"flux.c:mustaflux_compute()", "Utoprimgen", 1);
    if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
      if(debugfail>=1) dualfprintf(fail_file,"Failed to find inversion for MUSTAFORCEFLUX(left): nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
      domustaflux=0;
      break;
    }
    // get new primitive p_r
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, U_r, ptrstate_r, geom, dissmeasure, p_r, p_r,&newtonstats),"flux.c:mustaflux_compute()", "Utoprimgen", 1);
    if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
      if(debugfail>=1) dualfprintf(fail_file,"Failed to find inversion for MUSTAFORCEFLUX(right): nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
      domustaflux=0;
      break;
    }

    // get fluxes (state used for vchar() below so need full state)
    MYFUN(get_state(p_l, geom, ptrstate_l),"flux.c:flux_compute()", "get_state()", 1);
    MYFUN(primtoflux(UEVOLVE, p_l, ptrstate_l, dir, geom, F_l, NULL),"flux.c:flux_compute()","primtoflux_calc() dir=1/2 l", 1);

    MYFUN(get_state(p_r, geom, ptrstate_r),"flux.c:flux_compute()", "get_state()", 1);
    MYFUN(primtoflux(UEVOLVE, p_r, ptrstate_r, dir, geom, F_r, NULL),"flux.c:flux_compute()","primtoflux_calc() dir=1/2 l", 1);

    // now have p_l, p_r, U_l, U_r, F_l, F_r


#if(1||((WHICHFLUX==MUSTAHLL)||(WHICHFLUX==MUSTALAXF)))
    // wave speeds for new left-right states
    MYFUN(vchar_each(p_l, ptrstate_l, dir, geom, &cmaxeach_l[EOMSETMHD], &cmineach_l[EOMSETMHD], &cmaxeach_l[EOMSETRAD], &cmineach_l[EOMSETRAD], &cmaxeach_l[EOMSETRADFORDT], &cmineach_l[EOMSETRADFORDT],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
    MYFUN(vchar_each(p_r, ptrstate_r, dir, geom, &cmaxeach_r[EOMSETMHD], &cmineach_r[EOMSETMHD], &cmaxeach_r[EOMSETRAD], &cmineach_r[EOMSETRAD], &cmaxeach_r[EOMSETRADFORDT], &cmineach_r[EOMSETRADFORDT],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
    cminmax_calc(cmineach_l[EOMSETMHD],cmineach_r[EOMSETMHD],cmaxeach_l[EOMSETMHD],cmaxeach_r[EOMSETMHD],&cminmhd,&cmaxmhd,&ctopmhd);
    cminmax_calc(cmineach_l[EOMSETRAD],cmineach_r[EOMSETRAD],cmaxeach_l[EOMSETRAD],cmaxeach_r[EOMSETRAD],&cminrad,&cmaxrad,&ctoprad);
    cminmax_calc(cmineach_l[EOMSETRADFORDT],cmineach_r[EOMSETRADFORDT],cmaxeach_l[EOMSETRADFORDT],cmaxeach_r[EOMSETRADFORDT],&cminrad2,&cmaxrad2,&ctoprad2);
    // cmin, cmax, ctop for fluxes
    PLOOP(pliter,pl){
      if(pl<URAD0 || pl>URAD3){
        cmin[pl]=cminmhd;
        cmax[pl]=cmaxmhd;
        ctop[pl]=ctopmhd;
      }
      else{
        cmin[pl]=cminrad;
        cmax[pl]=cmaxrad;
        ctop[pl]=ctoprad;
      }
    }
#endif


    // get F
    // get flux using corrected left/right states
#if(WHICHFLUX==MUSTAFORCE)
    PLOOP(pliter,pl) cforce[pl]=ctop[pl];
    forceflux_compute(dir,geom,cmin,cmax,ctop,cforce,p_l,p_r,U_l,U_r,F_l,F_r,F);
#elif(WHICHFLUX==MUSTAHLL)
    hllflux_compute(dir,geom,cmin,cmax,ctop,p_l,p_r,U_l,U_r,F_l,F_r,F);
#elif(WHICHFLUX==MUSTALAXF)
    PLOOP(pliter,pl) F[pl] = LAXFCOMPUTE(ctop[pl],U_l[pl],U_r[pl],F_l[pl],F_r[pl]);
#endif

  }// end musta loop


  return(domustaflux);
}





/// second (multi-cell) version of MUSTA flux
int musta2flux_compute(int dir,struct of_geom *geom, FTYPE *cmin_l, FTYPE *cmin_r, FTYPE *cmax_l, FTYPE *cmax_r, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F)
{
  int forceflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *cforce, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int hllflux_compute(int dir,struct of_geom *geom, FTYPE *cmin, FTYPE *cmax, FTYPE *ctop, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r,FTYPE *F_l,FTYPE *F_r, FTYPE *F);
  int pl,pliter;
  FTYPE umid[NPR],plnew[NPR],prnew[NPR],fmid[NPR];
  struct of_state state;
  struct of_state *ptrstate;
  int domustaflux;
  int mustaloop;
  FTYPE cmusta[NPR];
  FTYPE *Unow,*Fnow,*pnow;
  FTYPE *Unow_l, *Unow_r, *Fnow_l, *Fnow_r, *pnow_l, *pnow_r;
  FTYPE *cmaxnow,*cminnow,*ctopnow;
  FTYPE Umusta[2*(NUMLOCALCELLS+1)][NPR],Fmusta[2*(NUMLOCALCELLS+1)][NPR],Fmusta_edge[2*(NUMLOCALCELLS+1)][NPR],pmusta[2*(NUMLOCALCELLS+1)][NPR];
  FTYPE cmaxmusta[2*(NUMLOCALCELLS+1)][NPR],cminmusta[2*(NUMLOCALCELLS+1)][NPR];
  FTYPE cmaxnow_l[NPR], cmaxnow_r[NPR], cminnow_l[NPR], cminnow_r[NPR];
  int ignorecourant;
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  FTYPE otherF[NPR];
  int mustacellloop;
  int ms,me;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; 


  // setup multi-cell MUSTA
  // left state
  for(mustacellloop=0;mustacellloop<=NUMLOCALCELLS;mustacellloop++){
    PLOOP(pliter,pl){
      Umusta[mustacellloop][pl]=U_l[pl];
      Fmusta[mustacellloop][pl]=F_l[pl];
      pmusta[mustacellloop][pl]=p_l[pl];
      cmaxmusta[mustacellloop][pl]=cmax_l[pl];
      cminmusta[mustacellloop][pl]=cmin_l[pl];
    }
  }
  // right state
  for(mustacellloop=NUMLOCALCELLS+1;mustacellloop<=2*NUMLOCALCELLS+1;mustacellloop++){
    PLOOP(pliter,pl){
      Umusta[mustacellloop][pl]=U_r[pl];
      Fmusta[mustacellloop][pl]=F_r[pl];
      pmusta[mustacellloop][pl]=p_r[pl];
      cmaxmusta[mustacellloop][pl]=cmax_r[pl];
      cminmusta[mustacellloop][pl]=cmin_r[pl];
    }
  }

  
  // default is to do musta unless musta fails
  domustaflux=1;



  ///////////////////
  //
  // loop over musta stages
  //
  // final flux is located on edge, while left/right states are "cell centered"
  //
  ///////////////////
  for(mustaloop=0;mustaloop<NUMMULTIMUSTAITERS;mustaloop++){

    /*
    // range of fluxes required
    if(mustaloop<(NUMMULTIMUSTAITERS-1)/2){
    ms=NUMLOCALCELLS+1-mustaloop;
    me=NUMLOCALCELLS+2+mustaloop;
    }
    else if(mustaloop>=(NUMMULTIMUSTAITERS-1)/2){
    ms=mustaloop+2;
    me=NUMLOCALCELLS+1-mustaloop;
    }
    */

    ////////////////////////
    //
    // get predictor edge flux
    //
    ////////////////////////
    for(mustacellloop=1;mustacellloop<=2*NUMLOCALCELLS+1;mustacellloop++){
      //    for(mustacellloop=ms;mustacellloop<=me;mustacellloop++){

      // pointers to memory locations GETS RETRIEVED
      Unow_l=Umusta[mustacellloop-1];
      Unow_r=Umusta[mustacellloop];
      Fnow_l=Fmusta[mustacellloop-1];
      Fnow_r=Fmusta[mustacellloop];
      pnow_l=pmusta[mustacellloop-1];
      pnow_r=pmusta[mustacellloop];
      // assignments GETS RETRIEVED
      PLOOP(pliter,pl){
        cmaxnow_l[pl]=cmaxmusta[mustacellloop-1][pl];
        cmaxnow_r[pl]=cmaxmusta[mustacellloop][pl];
        cminnow_l[pl]=cminmusta[mustacellloop-1][pl];
        cminnow_r[pl]=cminmusta[mustacellloop][pl];
      }
      // Fnow on edge GETS SET
      Fnow=Fmusta_edge[mustacellloop];


      // get wave speeds and min/max/top versions of left/right states GETS RETRIEVED
      //      cminmax_calc(cminnow_l,cminnow_r,cmaxnow_l,cmaxnow_r,&cmin,&cmax,&ctop);
      // SUPERGODMARK: Not setup for radiation
      if(EOMRADTYPE!=EOMRADNONE){
        dualfprintf(fail_file,"musta2 not setup for radiation\n");
        myexit(1);
      }
#if(WHICHFLUX==MUSTAFORCE)
      PLOOP(pliter,pl) cmusta[pl]=cforce[pl];
      forceflux_compute(dir,geom,cmin,cmax,ctop,cforce,pnow_l,pnow_r,Unow_l,Unow_r,Fnow_l,Fnow_r,Fnow);
#elif(WHICHFLUX==MUSTAHLL)
      PLOOP(pliter,pl) cmusta[pl]=ctop[pl];
      hllflux_compute(dir,geom,cmin,cmax,ctop,pnow_l,pnow_r,Unow_l,Unow_r,Fnow_l,Fnow_r,Fnow);
#elif(WHICHFLUX==MUSTALAXF)
      PLOOP(pliter,pl) cmusta[pl]=ctop[pl];
      PLOOP(pliter,pl) Fnow[pl] = LAXFCOMPUTE(ctop[pl],Unow_l[pl],Unow_r[pl],Fnow_l[pl],Fnow_r[pl]);
#endif






    }// end loop over musta cells

    // Fnow in middle is final flux (NUMMULTIMUSTAITERS==1 gives same as no MUSTA)
    //  (NUMMULTIMUSTAITERS==2 gives MUSTA-1)
    if(mustaloop==NUMMULTIMUSTAITERS-1) break;

    ///////////////////////////
    //
    // update U_l U_r (Open Riemann fan)
    //
    /////////////////////////
    for(mustacellloop=1;mustacellloop<=2*NUMLOCALCELLS;mustacellloop++){


      // Unow, pnow, and Fnow at center GETS SET
      Unow=Umusta[mustacellloop];
      pnow=pmusta[mustacellloop];
      Fnow=Fmusta[mustacellloop];
      PLOOP(pliter,pl){
        cmaxnow[pl]=cmaxmusta[mustacellloop][pl];
        cminnow[pl]=cminmusta[mustacellloop][pl];
      }
      // Fnow_l and Fnow_r on edge GETS RETRIEVED
      Fnow_l=Fmusta_edge[mustacellloop];
      Fnow_r=Fmusta_edge[mustacellloop+1];



      // set wave speed
      PLOOP(pliter,pl){
#if(WHICHFLUX==MUSTAFORCE)
        cmusta[pl]=cforce[pl];
#elif(WHICHFLUX==MUSTAHLL)
        cmusta[pl]=ctop[pl];
#elif(WHICHFLUX==MUSTALAXF)
        cmusta[pl]=ctop[pl];
#endif
      }

      // Open Riemann fan
      PLOOP(pliter,pl) Unow[pl] -= MULTIMUSTACOEF*(Fnow_r[pl]-Fnow_l[pl])/cmusta[pl];


      // get new primitive pnow from Unow
      int eomtype=EOMDEFAULT;
      FTYPE dissmeasure=-1.0; // assume energy try ok
      int whichcap=CAPTYPEBASIC;
      int whichmethod=MODEDEFAULT;
      int modprim=0;
      int checkoninversiongas=CHECKONINVERSION;
      int checkoninversionrad=CHECKONINVERSIONRAD;
      MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UEVOLVE, Unow, NULL, geom, dissmeasure, pnow, pnow,&newtonstats),"flux.c:mustaflux_compute()", "Utoprimgen", 1);
      if(GLOBALMACP0A1(pflag,geom->i,geom->j,geom->k,FLAGUTOPRIMFAIL)){
        if(debugfail>=1) dualfprintf(fail_file,"Failed to find inversion for MUSTAFORCEFLUX(right): nstep=%ld t=%21.15g i=%d j=%d k=%d\n",nstep,t,geom->i,geom->j,geom->k);
        domustaflux=0;
        break;
      }
      

      // get corrected centered fluxes (state needed for vchar() so need full state)
      ptrstate=&state; // default
      MYFUN(get_state(pnow, geom, ptrstate),"flux.c:flux_compute()", "get_state()", 1);
      MYFUN(primtoflux(UEVOLVE, pnow, ptrstate, dir, geom, Fnow, NULL),"flux.c:flux_compute()","primtoflux_calc() dir=1/2 l", 1);
      // now have p, U, F at center again

      // get min and max wave speeds
      MYFUN(vchar(pnow, ptrstate, dir, geom, cmaxnow, cminnow, &ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);
      // SUPERGODMARK: Not setup for radiation

    }// end loop over musta cells


    if(domustaflux==0) break;

    ///////////////////////////
    //
    // Set boundary conditions (zeroth order copy -- outflow)
    //
    /////////////////////////

    PLOOP(pliter,pl){
      Umusta[0][pl]=Umusta[1][pl];
      Fmusta[0][pl]=Fmusta[1][pl];
      pmusta[0][pl]=pmusta[1][pl];
      
      Umusta[2*NUMLOCALCELLS+1][pl]=Umusta[2*NUMLOCALCELLS][pl];
      Fmusta[2*NUMLOCALCELLS+1][pl]=Fmusta[2*NUMLOCALCELLS][pl];
      pmusta[2*NUMLOCALCELLS+1][pl]=pmusta[2*NUMLOCALCELLS][pl];
    }

  }// end loop over musta stages


  ////////////////
  //
  // assign final edge flux as MUSTA solution
  //
  //////////////

  PLOOP(pliter,pl) F[pl]=Fmusta_edge[NUMLOCALCELLS+1][pl];



  return(domustaflux);
}



/// choose flux form
void choose_flux(int fluxmethodlocal, int i, int j, int k, int pl, FTYPE *laxffrac,FTYPE *hllfrac)
{
#if(FLUXADJUST)

#if(HYDROFLUXADJUSTONLY)
  if(pl<B1){
    if(GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)==HLLFLUX){ hllfrac[pl]=1.0; laxffrac[pl]=0.0; }
    else { hllfrac[pl]=0.0; laxffrac[pl]=1.0; }
  }
  else{
    if(fluxmethodlocal==HLLFLUX){
      hllfrac[pl]=1.0;
      laxffrac[pl]=0.0;
    }
    else if(fluxmethodlocal==LAXFFLUX){
      hllfrac[pl]=0.0;
      laxffrac[pl]=1.0;
    }
  }
#else
  if(GLOBALMACP0A1(pflag,i,j,k,FLAGREALFLUX)==HLLFLUX){ hllfrac[pl]=1.0; laxffrac[pl]=0.0; }
  else { hllfrac[pl]=0.0; laxffrac[pl]=1.0; }
#endif

#else
  if(fluxmethodlocal==HLLFLUX){
    hllfrac[pl]=1.0;
    laxffrac[pl]=0.0;
  }
  else if(fluxmethodlocal==LAXFFLUX){
    hllfrac[pl]=0.0;
    laxffrac[pl]=1.0;
  }
#endif

}

