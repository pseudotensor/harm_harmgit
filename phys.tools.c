
/*! \file phys.tools.c
     \brief All locally-related physics calculations
     // THINGS IN HERE ARE PER-POINT and use only LOCAL variables (no globals)
     */
#include "decs.h"





// making variables static inside functions made things slower if anything
// would have thought consruction/destruction operations would be important
// why static makes slower?
//#define VARSTATIC static
#define VARSTATIC 











/// Calculate fluxes in direction dir and conserved variable U
/// returntype==0 : flux with geometric factor geom->e (used by evolution code)
/// returntype==1 : flux with physical geometry factor geom->gdet (used by diagnostics)
/// see UtoU and source_conn()
/// Note that if MAXWELL==PRIMMAXWELL then primtoflux doesn't use b^\mu or b_\mu (bcon and bcov)
int primtoflux(int returntype, FTYPE *pr, struct of_state *q, int dir,
               struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs)
{
  int primtoflux_ma(int needentropy,int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs, FTYPE *fluxdiag, FTYPE *fluxdiagabs);
  int primtoflux_rad(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  void UtoU_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC FTYPE fluxinput[NPR],fluxinputma[NPR],fluxinputrad[NPR],fluxinputem[NPR];
  VARSTATIC FTYPE fluxinputabs[NPR],fluxinputmaabs[NPR],fluxinputradabs[NPR],fluxinputemabs[NPR];
  VARSTATIC FTYPE fluxdiag;
  VARSTATIC FTYPE fluxdiagabs;
  VARSTATIC int pl,pliter;



  // initialize fluxinputma and fluxinputem so individual functions only have to define non-zero terms
  PLOOP(pliter,pl) fluxinput[pl]=fluxinputma[pl]=fluxinputrad[pl]=fluxinputem[pl]=0.0;
  PLOOP(pliter,pl) fluxinputabs[pl]=fluxinputmaabs[pl]=fluxinputradabs[pl]=fluxinputemabs[pl]=0.0;
  fluxdiag=0.0;


  // define MA terms
  primtoflux_ma(1,&returntype, pr, q, dir, geom, fluxinputma, fluxinputmaabs, &fluxdiag, &fluxdiagabs);
  fluxinputma[UU+dir]+=fluxdiag; // add back to normal term
  fluxinputmaabs[UU+dir]+=fluxdiagabs; // add back to normal term
  // add up MA
  PLOOP(pliter,pl){
    fluxinput[pl] += fluxinputma[pl];
    fluxinputabs[pl] += fluxinputmaabs[pl];
  }

  if(EOMRADTYPE!=EOMRADNONE){
    // define RAD terms
    primtoflux_rad(&returntype, pr, q, dir, geom, fluxinputrad, fluxinputradabs);
    // add up RAD
    PLOOP(pliter,pl){
      fluxinput[pl] += fluxinputrad[pl];
      fluxinputabs[pl] += fluxinputradabs[pl];
    }
  }

  // define EM terms
  primtoflux_em(&returntype, pr, q, dir, geom, fluxinputem, fluxinputemabs);
  // add up EM
  PLOOP(pliter,pl){
    fluxinput[pl] += fluxinputem[pl];
    fluxinputabs[pl] += fluxinputemabs[pl];
  }


  // DEBUG:
  //  if((fabs(pr[U1])>1.0 || fabs(fluxinput[U1])>1.0) && (geom->i==10 || geom->i==9 || geom->i==11)&&(geom->k==0 || geom->k==-1 || geom->k==1)) PALLLOOP(pl) dualfprintf(fail_file,"ALLBEFORE: ijk=%d %d %d : pl=%d pr=%g flux=%21.15g\n",geom->i,geom->j,geom->k,pl,pr[pl],fluxinput[pl]);


  // convert from UNOTHING->returntype
  // notice that geometry comes after subtractions/additions of EOMs
  UtoU_fromunothing(returntype,geom,fluxinput,flux);
  if(fluxabs!=NULL) UtoU_fromunothing(returntype,geom,fluxinputabs,fluxabs);



  // DEBUG:
  //  PALLLOOP(pl) dualfprintf(fail_file,"ALLAFTER: pl=%d flux=%21.15g\n",pl,flux[pl]);


  return(0);

}




int primtoflux_nonradonly(int needentropy, FTYPE *pr, struct of_state *q, int dir,
               struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs)
{
  int primtoflux_ma(int needentropy,int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs, FTYPE *fluxdiag, FTYPE *fluxdiagabs);
  int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  VARSTATIC FTYPE fluxma[NPR],fluxem[NPR];
  VARSTATIC FTYPE fluxmaabs[NPR],fluxemabs[NPR];
  //  VARSTATIC FTYPE fluxdiag;
  //  VARSTATIC FTYPE fluxdiagabs;
  VARSTATIC int pl,pliter;
  int returntype=UNOTHING;


  // initialize fluxma and fluxem so individual functions only have to define non-zero terms
  PLOOP(pliter,pl) flux[pl]=fluxma[pl]=fluxem[pl]=0.0;
  PLOOP(pliter,pl) fluxmaabs[pl]=fluxemabs[pl]=0.0;
  if(fluxabs!=NULL) PLOOP(pliter,pl) fluxabs[pl]=0.0;
  //  fluxdiag=0.0;
  //  fluxdiagabs=0.0;


  // define MA terms
  primtoflux_ma(needentropy,&returntype, pr, q, dir, geom, fluxma, fluxmaabs, NULL, NULL);
  //  fluxma[UU+dir]+=fluxdiag; // add back to normal term
  //  fluxmaabs[UU+dir]+=fluxdiagabs; // add back to normal term
  // add up MA
  PLOOP(pliter,pl){
    flux[pl] += fluxma[pl];
    if(fluxabs!=NULL) fluxabs[pl] += fluxmaabs[pl];
  }


  // define EM terms
  primtoflux_em(&returntype, pr, q, dir, geom, fluxem, fluxemabs);
  // add up EM
  PLOOP(pliter,pl){
    flux[pl] += fluxem[pl];
    if(fluxabs!=NULL) fluxabs[pl] += fluxemabs[pl];
  }


  // returns UNOTHING form

  return(0);
}

int primtoflux_radonly(FTYPE *pr, struct of_state *q, int dir,
               struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs)
{
  int primtoflux_rad(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  VARSTATIC FTYPE fluxrad[NPR];
  VARSTATIC FTYPE fluxradabs[NPR];
  //  VARSTATIC FTYPE fluxdiag;
  VARSTATIC int pl,pliter;
  int returntype=UNOTHING;



  PLOOP(pliter,pl) flux[pl]=fluxrad[pl]=0.0;
  if(fluxabs!=NULL) PLOOP(pliter,pl) fluxabs[pl]=0.0;


  if(EOMRADTYPE!=EOMRADNONE){
    // define RAD terms
    primtoflux_rad(&returntype, pr, q, dir, geom, fluxrad, fluxradabs);
    // add up RAD
    PLOOP(pliter,pl){
      flux[pl] += fluxrad[pl];
      if(fluxabs!=NULL) fluxabs[pl] += fluxradabs[pl];
    }
  }

  // returns UNOTHING form

  return(0);
}




/* calculate fluxes in direction dir and conserved variable U; these
   are always needed together, so there is no point in calculated the
   stress tensor twice */

/// returntype==0 : flux with geometric factor geom->e (used by evolution code)
/// returntype==1 : flux with physical geometry factor geom->gdet (used by diagnostics)
/// see UtoU and source_conn()
/// fluxdir = true flux direction even if passing back conserved quantity (fundir==TT)
int primtoflux_splitmaem(int returntype, FTYPE *pr, struct of_state *q, int fluxdir, int fundir, struct of_geom *geom, FTYPE *fluxma, FTYPE *fluxem)
{
  int primtoflux_ma(int needentropy,int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs, FTYPE *fluxdiag, FTYPE *fluxdiagabs);
  int primtoflux_rad(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
  void UtoU_ma_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_rad_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_em_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC FTYPE fluxinput[NPR],fluxinputma[NPR],fluxinputrad[NPR],fluxinputem[NPR];
  VARSTATIC FTYPE fluxdiag;
  VARSTATIC int pl,pliter;
  FTYPE fluxrad[NPR];

  // initialize fluxinputma and fluxinputrad+fluxinputem so individual functions only have to define non-zero terms
  PLOOP(pliter,pl) fluxinputma[pl]=fluxinputrad[pl]=fluxinputem[pl]=0.0;

  // define MA terms
  primtoflux_ma(1,&returntype, pr, q, fundir, geom, fluxinputma, NULL, &fluxdiag, NULL);

  // SUPERGODMARK CHANGINGMARK
  // Note that pressure part of flux has 0 conserved quantity associated with it from the point of view of this output and correct HLL/LAXF calculation
#if(1)
  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag;
  else fluxinputma[UU]+=fluxdiag;

  //  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fundir)]=fluxdiag; // only need to return single diagonal value and store it inside field part for now
  //  else fluxinputma[UU]+=fluxdiag; // then really part of conserved quantity and then don't separate it

#elif(0)

  // not as robust to put conserved quantity here such that dissipative term is not together with momentum-flux_dir_dir term
  fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag; // only need to return single diagonal value and store it inside field part for now
#else
  // not as robust to put conserved quantity here such that dissipative term is not together with momentum-flux_dir_dir term
  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag;
  else{
    fluxinputma[UU]+=0.5*fluxdiag;
    fluxinputma[FLUXSPLITPMA(fluxdir)]=0.5*fluxdiag;
  }
#endif
  // convert from UNOTHING->returntype
  // notice that geometry comes after subtractions/additions of EOMs
  //  UtoU_ma(UNOTHING,returntype,geom,fluxinputma,fluxma); // properly converts separate diagonal flux in FLUXSPLITPMA(fluxdir)
  UtoU_ma_fromunothing(returntype,geom,fluxinputma,fluxma);


  if(EOMRADTYPE!=EOMRADNONE){

    // define RAD terms (as separate fluid from MHD fluid)
    primtoflux_rad(&returntype, pr, q, fundir, geom, fluxinputrad, NULL);

    //  UtoU_rad(UNOTHING,returntype,geom,fluxinputrad,fluxrad);
    UtoU_rad_fromunothing(returntype,geom,fluxinputrad,fluxrad);

    PLOOP(pliter,pl) fluxma[pl]+=fluxrad[pl]; // KORAL: Just add RAD to MA for this split approach
  }

  // define EM terms
  primtoflux_em(&returntype, pr, q, fundir, geom, fluxinputem, NULL);
  //  UtoU_em(UNOTHING,returntype,geom,fluxinputem,fluxem);
  UtoU_em_fromunothing(returntype,geom,fluxinputem,fluxem);



  return(0);

}



/// matter only terms (as if B=0)
int primtoflux_ma(int needentropy,int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs, FTYPE *fluxdiag, FTYPE *fluxdiagabs)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  int ynuflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum);
  int ylflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum);
  int yflflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum);
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs);
  int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux, FTYPE *entropyfluxabs);
  VARSTATIC FTYPE fluxdiagpress[NPR]; // temp var
  VARSTATIC FTYPE fluxdiagpressabs[NPR]; // temp var
  VARSTATIC int pl,pliter;

  // USE of SPLITNPR here is simplified for speed
  // assume if quantity is bracketed that including it


#if(SPLITNPR)
  if(nprlist[nprstart]<=RHO && nprlist[nprend]>=RHO)
#endif
    massflux_calc(pr, dir, q, &flux[RHO], &fluxabs[RHO]); // fills RHO only


  // GODMARK WTF!  Problems with code (compiling?) with this
  // if(mhd_calc(pr,dir,geom,q,mhd,mhdabs)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);

  // MHD stress-energy tensor w/ first index up, second index down.
#if(SPLITNPR)
  if(nprlist[nprstart]<=UU && nprlist[nprend]>=U3)
#endif
    {
      mhd_calc_ma(pr, dir, geom, q, &flux[UU], &fluxabs[UU], &fluxdiagpress[UU], &fluxdiagpressabs[UU]); // fills flux[UU->U3] and fluxdiagonal[UU->U3]
      if(fluxdiag!=NULL) *fluxdiag = fluxdiagpress[UU+dir];
      if(fluxdiagabs!=NULL) *fluxdiagabs = fluxdiagpressabs[UU+dir];
    }

  //  dualfprintf(fail_file,"fluxdiagabs=%g\n",*fluxdiagabs);




#if(YFL1>=0)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YFL1 && nprlist[nprend]>=YFL1)
#endif
    yflflux_calc(geom,pr, dir, q, &flux[YFL1], &fluxabs[YFL1],YFL1); // fills YFL1 only
#endif

#if(YFL2>=0)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YFL2 && nprlist[nprend]>=YFL2)
#endif
    yflflux_calc(geom,pr, dir, q, &flux[YFL2], &fluxabs[YFL2],YFL2); // fills YFL2 only
#endif

#if(YFL3>=0)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YFL3 && nprlist[nprend]>=YFL3)
#endif
    yflflux_calc(geom,pr, dir, q, &flux[YFL3], &fluxabs[YFL3],YFL3); // fills YFL3 only
#endif





#if(DOYL!=DONOYL)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YL && nprlist[nprend]>=YL)
#endif
    ylflux_calc(geom,pr, dir, q, &flux[YL], &fluxabs[YL],YL); // fills YL only
#endif
#if(DOYNU!=DONOYNU)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YNU && nprlist[nprend]>=YNU)
#endif
    ynuflux_calc(geom, pr, dir, q, &flux[YNU], &fluxabs[YNU],YNU); // fills YNU only
#endif


#if(DOENTROPY!=DONOENTROPY)
  if(needentropy){
#if(SPLITNPR)
    if(nprlist[nprstart]<=ENTROPY && nprlist[nprend]>=ENTROPY)
#endif
      entropyflux_calc(pr, dir, q, &flux[ENTROPY], &fluxabs[ENTROPY]); // fills ENTROPY only
    
    // below is special for utoprim() 5D version for full entropy evolution and inversion
    if(*returntype==UENTROPY){
      flux[UU]=flux[ENTROPY]; // overwrite for utoprim()
      fluxabs[UU]=fluxabs[ENTROPY]; // overwrite for utoprim()
      if(fluxdiag!=NULL) *fluxdiag = 0.0; // overwrite for utoprim()
      if(fluxdiagabs!=NULL) *fluxdiagabs = 0.0; // overwrite for utoprim()
      *returntype=UNOTHING; // reset returntype for UtoU
    }
  }

#endif


  // DEBUG:
  //  PALLLOOP(pl) dualfprintf(fail_file,"ALL: pl=%d flux=%21.15g\n",pl,flux[pl]);


  return (0);
}



/// radiation terms (as if rho=u=p=0)
int primtoflux_rad(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs)
{
  int yflflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum);

  // Radiation stress-energy tensor w/ first index up, second index down.
#if(SPLITNPR)
#error "primtoflux_rad not setup for SPLITNPR"
#endif

  if(EOMRADTYPE!=EOMRADNONE){
    mhd_calc_rad(pr, dir, geom, q, &flux[URAD0], &fluxabs[URAD0]); // fills URAD0->URAD3
  }
  // else don't fill flux[RAD0->RAD3] since assume entries don't exist

#if(YFL4>=0)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YFL4 && nprlist[nprend]>=YFL4)
#endif
    yflflux_calc(geom,pr, dir, q, &flux[YFL4], &fluxabs[YFL4],YFL4); // fills YFL4 only
#endif

#if(YFL5>=0)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YFL5 && nprlist[nprend]>=YFL5)
#endif
    yflflux_calc(geom,pr, dir, q, &flux[YFL5], &fluxabs[YFL5],YFL5); // fills YFL5 only
#endif



  return (0);
}


/// electromagnetic terms (as if rho=u=p=0)
int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  //  FTYPE dualf[NDIM];
  //  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs);
  //  int advectedscalarflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum);


  // USE of SPLITNPR here is simplified for speed
  // assume if quantity is bracketed that including it


  // GODMARK WTF!  Problems with code (compiling?) with this
  // if(mhd_calc(pr,dir,geom,q,mhd,mhdabs)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);

  // MHD stress-energy tensor w/ first index up, second index down.
#if(SPLITNPR)
  if(nprlist[nprstart]<=UU && nprlist[nprend]>=U3)
#endif
    mhd_calc_em(pr, dir, geom, q, &flux[UU], &fluxabs[UU]); // fills UU->U3



#if(SPLITNPR)
  if(nprlist[nprstart]<=B1 && nprlist[nprend]>=B3)
#endif
    dualfaradayspatial_calc(pr,dir,q,&flux[B1],&fluxabs[B1]); // fills B1->B3


#if(DEBUGNSBH)
  // DEBUG:
  if(geom->i==26 && geom->j==40 && dir==1){
    dualfprintf(fail_file,"INprimetoflux_em: %21.15g %21.15g %21.15g\n",flux[B1],flux[B2],flux[B3]);
  }
#endif


  return (0);
}


int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs)
{
  //  particle number flux 
  *massflux = pr[RHO] * q->ucon[dir];
  *massfluxabs=fabs(*massflux);

  // DEBUG:
  //  dualfprintf(fail_file,"massflux: %d %21.15g %21.15g\n",dir,pr[RHO],q->ucon[dir]);

  return(0);
}


/// flux associated with Y_fl variable
int yflflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum)
{
  FTYPE prforadvect;

  FTYPE udir=0.0;
  if(pnum==YFL1 || pnum==YFL2 || pnum==YFL3) udir=q->ucon[dir];
  if(pnum==YFL4 || pnum==YFL5) udir=q->uradcon[dir];

  prforadvect = pr[pnum];

#if(DOYFL==1)
  // equation is d_t(\rho_0 u^t y) = d_i (\rho_0 u^i y)
  // get flux associated with Y_L
  *advectedscalarflux = prforadvect * pr[RHO]*udir; // y form
#elif(DOYFL==2)
  // equation is d_t(rhofl u^t) = d_i (rhofl u^i)
  // equation is d_t(T^t_t/u^t u^t) = d_i (T^t_t/u^t u^i)
  // equation is d_t(T^t_\phi/u^t u^t) = d_i (T^t_\phi/u^t u^i)
  // etc.
  *advectedscalarflux = prforadvect * udir; // rho form
#endif


  *advectedscalarfluxabs = fabs(*advectedscalarflux);

  return(0);
}

/// flux associated with Y_L variable
int ylflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum)
{
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs);
  VARSTATIC FTYPE massflux;
  VARSTATIC FTYPE massfluxabs;
  FTYPE prforadvect;

  // get mass flux
  massflux_calc(pr, dir, q, &massflux, &massfluxabs);

#if(WHICHEOS==KAZFULL)
  yl2advect_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),pr[YL],pr[YNU],&prforadvect);
#else
  prforadvect = pr[pnum];
#endif

  // get flux associated with Y_L
  *advectedscalarflux = prforadvect * massflux;

  *advectedscalarfluxabs = fabs(*advectedscalarflux);

  return(0);
}


/// flux asociated with Ynu variable
int ynuflux_calc(struct of_geom *ptrgeom, FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum)
{
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs);
  VARSTATIC FTYPE massflux, massfluxabs;
  FTYPE prforadvect;


  // get mass flux
  massflux_calc(pr, dir, q, &massflux, &massfluxabs);

#if(WHICHEOS==KAZFULL)
  ynu2advect_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),pr[YL],pr[YNU],&prforadvect);
#else
  prforadvect = pr[pnum];
#endif

  // get flux associated with Y_\nu
  *advectedscalarflux = prforadvect * massflux;

  *advectedscalarfluxabs = fabs(*advectedscalarflux);

  return(0);
}



/// flux of scalar
int advectedscalarflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, FTYPE *advectedscalarfluxabs, int pnum)
{
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux, FTYPE *massfluxabs);
  VARSTATIC FTYPE massflux;
  VARSTATIC FTYPE massfluxabs;

  //  yflx/yl/ynu/etc. per unit rest-mass flux 

  // entropy=entropy per unit volume, where conserved quantity is specific entropy:
  // d/d\tau(entropy/rho)=0
  // -> \nabla_\mu(entropy u^\mu)=0

  massflux_calc(pr, dir, q, &massflux, &massfluxabs);
  *advectedscalarflux = pr[pnum] * massflux;
  *advectedscalarfluxabs = fabs(*advectedscalarflux);

  return(0);
}



/// flux of specific entropy
int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux, FTYPE *entropyfluxabs)
{

  // get entropy
  //  entropy_calc(ptrgeom,pr,&entropy); // now done in get_state_thermodynamics()
  //  entropy per unit rest-mass flux 
  // entropy=entropy per unit volume, where conserved quantity is specific entropy:
  // d/d\tau(entropy/rho)=0
  // -> \nabla_\mu(entropy u^\mu)=0

  // SUPERKORALTODO: With an non-zero \kappa, thermal flux exists between fluid and radiation.  Below does not account for this.  See http://adsabs.harvard.edu/abs/2012MNRAS.426.1613R  APPENDIX A.  They don't give general expression, but should be a G_\mu term that enters.
  //
  // d/d\tau(entropy/rho)=d/\tau SOURCE
  // -> \nabla_\mu(entropy u^\mu)= ...


  // DEBUG:
  //dualfprintf(fail_file,"entropy=%21.15g dir=%d ucondir=%21.15g\n",entropy,dir,q->ucon[dir]);

  *entropyflux = (q->entropy) * (q->ucon[dir]);

  *entropyfluxabs = fabs(*entropyflux);

  return(0);
}




/// spatial part of dualfaraday
int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf, FTYPE *dualfabs)
{
  VARSTATIC FTYPE dualffull[NDIM];


  dualfullfaraday_calc(pr,dir,q,dualffull);
  dualf[0]=dualffull[1];
  dualf[1]=dualffull[2];
  dualf[2]=dualffull[3];

  dualfabs[0]=fabs(dualf[0]);
  dualfabs[1]=fabs(dualf[1]);
  dualfabs[2]=fabs(dualf[2]);


  return(0);

}



/// Notation for HARM paper and JCM's GRFFE paper is that:
/// B^\mu = \eta_\nu *F^{\nu\mu} and for lab-frame we chose \eta_\nu={-1,0,0,0}
///
/// One can then show that b^i u^j - b^j u^i = B^i v^j - B^j v^i where
/// b^\mu = u_\nu *F^{\nu\mu} and v^i = u^i/u^t
///
/// The form using B^i and v^i avoids catastrophic cancellation because otherwise
/// starting with B^i (as primitive) and converting to b^i leads to u^i u^j (u.B)/u^t term that cancels exactly
/// Since this is quite a high order term, then for highly relativistic flows this causes catastrophic
/// cancellation issues, so this is why we use the primitives directly even if more complicated looking and more expensive to compute stress tensor or maxwell tensor
///
/// returns \dF^{\mu dir}
/// well, actually returns dualffull[dir], so gives columns instead of rows
int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualffull)
{

#if(MAXWELL==GENMAXWELL)
  //  dual of Maxwell tensor 
  dualffull[0] = q->bcon[0] * q->ucon[dir] - q->bcon[dir] * q->ucon[0];
  dualffull[1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  dualffull[2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  dualffull[3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];
#elif(MAXWELL==PRIMMAXWELL)
  if(dir>0){
    //  dual of Maxwell tensor 
    // dir refers to the direction of the derivative of the dualffull
    // B1,B2,B3 refers to LHS of equation dB^i/dt
    // due to antisymmetry, dir==i is 0
    dualffull[0] = - pr[B1+dir-1] ; // dualffull[i]=\dF^{i dir} where \dF^{0 dir} =-B^{dir}
    dualffull[1] = (pr[B1] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[1])/q->ucon[0];
    dualffull[2] = (pr[B2] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[2])/q->ucon[0];
    dualffull[3] = (pr[B3] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[3])/q->ucon[0];
  }
  else{
    dualffull[0] = 0;
    dualffull[1] = pr[B1];
    dualffull[2] = pr[B2];
    dualffull[3] = pr[B3];
  }
#endif

  return(0);

}


/// dual of Maxwell tensor
/// returns \dF^{\mu \nu}
int Mcon_calc(FTYPE *pr, struct of_state *q, FTYPE (*Mcon)[NDIM])
{
  VARSTATIC int j,k;
  VARSTATIC FTYPE vcon[NDIM];


#if(MAXWELL==GENMAXWELL)

  DLOOP(j,k) Mcon[j][k] = q->bcon[j] * q->ucon[k] - q->bcon[k] * q->ucon[j];

#elif(MAXWELL==PRIMMAXWELL)

  // diagonal is 0
  DLOOPA(j) Mcon[j][j]=0.0;

  // space-time terms
  SLOOPA(k) {
    // \dF^{it} = B^i = pr[B1+i-1]
    Mcon[k][0] = pr[B1+k-1] ; 
    Mcon[0][k] = - Mcon[k][0] ;
  }

  // get v^i
  SLOOPA(k) vcon[k]= q->ucon[k]/q->ucon[TT];

  // space-space terms
  //  SLOOP(j,k) Mcon[j][k] = (pr[B1+j-1] * vcon[k] - pr[B1+k-1] * vcon[j]);
  // optimize
  Mcon[1][2] = (pr[B1] * vcon[2] - pr[B2] * vcon[1]);
  Mcon[1][3] = (pr[B1] * vcon[3] - pr[B3] * vcon[1]);
  Mcon[2][3] = (pr[B2] * vcon[3] - pr[B3] * vcon[2]);
  Mcon[2][1] = -Mcon[1][2];
  Mcon[3][1] = -Mcon[1][3];
  Mcon[3][2] = -Mcon[2][3];
  Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0.0;

#endif


  return(0);

}



/// returns entire space-time(NDIM in size) / EOM(NPR in size) matrix
int primtofullflux(int returntype, FTYPE *pr, struct of_state *q,
                   struct of_geom *ptrgeom, FTYPE (*flux)[NPR], FTYPE (*fluxabs)[NPR])
{
  VARSTATIC int j;
  
  // j=0,1,2,3 corresponding to U^j_\nu , where \nu corresponds to all EOMs and j to space-time for each
  // 1 stands for obey nprlist
  DLOOPA(j) primtoflux(returntype,pr,q,j,ptrgeom,flux[j],fluxabs[j]);

  return(0);
}


///  calculate "conserved" quantities 
int primtoU(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *geom,
            FTYPE *U,
            FTYPE *Uabs)
{
  MYFUN(primtoflux(returntype,pr, q, 0, geom, U, Uabs) ,"phys.c:primtoU()", "primtoflux_calc() dir=0", 1);

  return (0);
}




void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_gen_allgdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC int removemass;

#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
  removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#else
  // otherwise removemass = 0 effectively
  removemass = 0;
#endif

#if(WHICHEOM!=WITHGDET)
  UtoU_gen_gengdet(removemass, whichmaem, inputtype, returntype,ptrgeom,Uin, Uout);
#else
  UtoU_gen_allgdet(removemass, whichmaem, inputtype, returntype,ptrgeom,Uin, Uout);
#endif

}


void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC int removemass;


#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
  removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#else
  // otherwise removemass = 0 effectively
  removemass = 0;
#endif


#if(WHICHEOM!=WITHGDET)
  UtoU_gen_gengdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Uin, Uout);
#else
  UtoU_gen_allgdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Uin, Uout);
#endif

}


/// standardized U form is geometry free and 
/// \rho u^t , T^t_\nu , *F^{it}
/// convert one form of U(or component of Flux) to another form
/// UtoU controls meaning of todo and REMOVERESTMASSFROMUU.
/// present order means start with geometry-free EOMs, add/subtract them, THEN geometry is assigned to that list of new EOMs.
/// can choose to change order so that add geometry terms, THEN add/subtract them.  Rest of code shouldn't care (except source_conn()'s first connection)
/// if SPLITNPR, operates on other primitives but only changes Uout, not Uin, so doesn't matter if change all conserved quantities
/// whichmaem == ISEMONLY ISMAONLY ISMAANDEM
void UtoU_gen_gengdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout);
  VARSTATIC FTYPE Ugeomfree[NPR];
  VARSTATIC int pl,pliter;




  if(inputtype==returntype){
    // then just copy
    PLOOP(pliter,pl) Uout[pl]=Uin[pl];
  }
  else{

    /////////////////////
    //
    // input
    //
    
    // set basic transformation
    PLOOP(pliter,pl) Ugeomfree[pl]=Uin[pl];
    
    // now fine-tune transformation
    switch(inputtype){
      
    case UEVOLVE:
      
      // get igdet
      set_igdet(ptrgeom);
      
      PLOOP(pliter,pl) Ugeomfree[pl] *= ptrgeom->IEOMFUNCNOSINGMAC(pl);
      
      //    dualfprintf(fail_file,"diff: %21.15g %21.15g\n",ptrgeom->EOMFUNCMAC(UU),MACP0A1(gdetvol,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p)); 
#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
      if(removemass){
        // go back to standard stress-energy tensor form
        Ugeomfree[UU]  +=  - Ugeomfree[RHO] ; // - means adding back rest-mass
      }
#endif
      break;
    case UDIAG:
      // get igdet
      set_igdet(ptrgeom);
      
      PLOOP(pliter,pl) Ugeomfree[pl] *= ptrgeom->igdetnosing;
      break;
    case UNOTHING:
      break;
    default:
      dualfprintf(fail_file,"UtoU: No such inputtype=%d\n",inputtype);
      myexit(246364);
      break;
    }
    

    // at this point, Ugeomfree is geometry-free standard form of conserved quantities
    // output
    UtoU_gen_gengdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Ugeomfree, Uout);
    
  }

}





void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout)
{    
  VARSTATIC int pl,pliter;


  /////////////////////////
  //
  // output
  //
  switch(returntype){
  case UEVOLVE:
#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
    if(removemass){ // diagnostics want normal stress-energy tensor
      // "subtract" rest-mass
      // should be done on geometry-free version
      Ugeomfree[UU] += Ugeomfree[RHO];
    }
#endif
    PLOOP(pliter,pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->EOMFUNCMAC(pl);
    break;
  case UDIAG:
    PLOOP(pliter,pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->gdet;
    break;
  case UNOTHING:
    PLOOP(pliter,pl) Uout[pl]=Ugeomfree[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such returntype=%d\n",returntype);
    myexit(724365);
    break;
  }

}




/// simplified version that assumes WHICHEOM==WITHGDET
void UtoU_gen_allgdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  VARSTATIC FTYPE Ugeomfree[NPR];
  VARSTATIC int pl,pliter;
  void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout);



#if(WHICHEOM!=WITHGDET)
  dualfprintf(fail_file,"Using UtoU_gen_allgdet() but WHICHEOM!=WITHGDET\n");
  myexit(342968346);
#endif
  
  
  
#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
  removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#endif// otherwise removemass = 0 effectively
  
  
  /////////////////////
  //
  // input
  //
  
  // now fine-tune transformation
  switch(inputtype){
    
  case UEVOLVE:
  case UDIAG:
    
    // get igdet
    set_igdetsimple(ptrgeom);
    
    PLOOP(pliter,pl) Ugeomfree[pl] = Uin[pl]*ptrgeom->igdetnosing;
    
#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
    if(removemass){
      // go back to standard stress-energy tensor form
      Ugeomfree[UU]  +=  - Ugeomfree[RHO] ; // - means adding back rest-mass
    }
#endif
    break;
  case UNOTHING:
    PLOOP(pliter,pl) Ugeomfree[pl] = Uin[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such inputtype=%d\n",inputtype);
    myexit(6746364);
    break;
  }

  
  // at this point, Ugeomfree is geometry-free standard form of conserved quantities
  // output
  UtoU_gen_allgdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Ugeomfree, Uout);
  

}



void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout)
{    
  VARSTATIC int pl,pliter;

  /////////////////////////
  //
  // output
  //
  switch(returntype){
  case UEVOLVE:
  case UDIAG:
#if((REMOVERESTMASSFROMUU==1)&&(DOEVOLVERHO))
    if(removemass){ // diagnostics want normal stress-energy tensor
      // "subtract" rest-mass
      // should be done on geometry-free version
      Ugeomfree[UU] += Ugeomfree[RHO];
    }
#endif
    PLOOP(pliter,pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->gdet;
    break;
  case UNOTHING:
    PLOOP(pliter,pl) Uout[pl]=Ugeomfree[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such returntype=%d\n",returntype);
    myexit(924636);
    break;
  }

}


void UtoU_evolve2diag(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int pl,pliter;

#if(WHICHEOM==WITHGDET)
  // then UEVOLVE and UDIAG have same geometry (geom.g)
  PALLLOOP(pl) Uout[pl]=Uin[pl];
#else
  UtoU_gen(ISMAANDEM, inputtype, returntype, ptrgeom, Uin, Uout);
#endif

}


void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISMAANDEM, inputtype, returntype, ptrgeom, Uin, Uout);

}

void UtoU_ma(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISMAONLY, inputtype, returntype, ptrgeom, Uin, Uout);

}

void UtoU_rad(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISRADONLY, inputtype, returntype, ptrgeom, Uin, Uout);

}

void UtoU_em(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISEMONLY, inputtype, returntype, ptrgeom, Uin, Uout);

}


void UtoU_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISMAANDEM, returntype, ptrgeom, Uin, Uout);

}

void UtoU_ma_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISMAONLY, returntype, ptrgeom, Uin, Uout);

}

void UtoU_rad_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISRADONLY, returntype, ptrgeom, Uin, Uout);

}

void UtoU_em_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISEMONLY, returntype, ptrgeom, Uin, Uout);

}



/// standardized primitive form is assumed to be
/// \rho , u, \tilde{u}^\mu , *F^{it}=B^i
/// where \tilde{u} is relative 4-velocity, as relative to $n_\mu = (-\alpha,0,0,0)$ and $\alpha^2=-1/g^{tt}$.
/// For any space-time with no time-like curves this 4-velocity is always single-valued (i.e. unique).  It can also take on any value, so a reasonable primitive quantity.
/// convert from one primitive form to another
void PtoP(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *pin, FTYPE *pout)
{
  VARSTATIC FTYPE pstandard[NPR];
  VARSTATIC int pl,pliter;



}


/// calculate magnetic field four-vector
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  VARSTATIC int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j <= 3; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}

/// inverse of bcon_calc()
void Bcon_calc(struct of_state *q, FTYPE*B)
{
  VARSTATIC FTYPE uu0,uu1,uu2,uu3;
  VARSTATIC FTYPE ud0,ud1,ud2,ud3;
  VARSTATIC FTYPE bu1,bu2,bu3;
  VARSTATIC FTYPE denom;


  uu0=q->ucon[TT];
  uu1=q->ucon[RR];
  uu2=q->ucon[TH];
  uu3=q->ucon[PH];

  ud0=q->ucov[TT];
  ud1=q->ucov[RR];
  ud2=q->ucov[TH];
  ud3=q->ucov[PH];

  bu1=q->bcon[RR];
  bu2=q->bcon[TH];
  bu3=q->bcon[PH];
 
  denom=1.0/(1.0+ud1*uu1+ud2*uu2+ud3*uu3);
  
  B[1]=uu0*(-(bu2*ud2+bu3*ud3)*uu1+bu1*(1.0+ud2*uu2+ud3*uu3))*denom;
  B[2]=uu0*(-(bu1*ud1+bu3*ud3)*uu2+bu2*(1.0+ud1*uu1+ud3*uu3))*denom;
  B[3]=uu0*(-(bu2*ud2+bu2*ud2)*uu3+bu3*(1.0+ud2*uu2+ud1*uu1))*denom;


}


/// convert (e^\mu=0 case) b^\mu and (3-velocity in coordinate lab frame) v^\mu to pr
void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr)
{
  void Bcon_calc(struct of_state *q, FTYPE*B);
  VARSTATIC int pl,pliter;
  VARSTATIC struct of_state q;
  VARSTATIC FTYPE prim[NPR];
  VARSTATIC FTYPE ucon[NDIM];


  // go ahead and get pr velocity
  PLOOP(pliter,pl) prim[pl]=0.0;
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  //  vcon2pr(WHICHVEL,vcon,geom,pr); // need u^\mu, so do below instead
  ucon_calc_3vel(prim,geom,q.ucon,q.others);
  ucon2pr(WHICHVEL,q.ucon,geom,pr); // fills pr[U1->U3]
  
  //  q.ucon[TT]=ucon[TT];
  //  q.ucon[RR]=ucon[RR];
  //  q.ucon[TH]=ucon[TH];
  //  q.ucon[PH]=ucon[PH];
  
  lower_vec(q.ucon,geom,q.ucov);

  //  q.bcon[TT]=bcon[TT]; // not used below
  q.bcon[RR]=bcon[RR];
  q.bcon[TH]=bcon[TH];
  q.bcon[PH]=bcon[PH];
  
  Bcon_calc(&q,&pr[B1-1]); // &pr[B1-1] since Bcon_calc() fills 1-3



}





/// MHD stress tensor, with first index up, second index down
/// mhd^dir_j
void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
  void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);

#if(REMOVERESTMASSFROMUU==2)
  mhd_calc_norestmass(pr, dir, geom, q, mhd, mhdabs);
#else
  mhd_calc_0(pr, dir, geom, q, mhd, mhdabs);
#endif

}

/// MHD stress tensor, with first index up, second index down
/// mhd^dir_j
/// understood that mhddiagpress only contains non-zero element on mhddiagpress[dir] and all others should be 0.0
void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs)
{
  void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);
  void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);

#if(REMOVERESTMASSFROMUU==2)
  mhd_calc_norestmass_ma(pr, dir, geom, q, mhd, mhdabs, mhddiagpress, mhddiagpressabs);
#else
  mhd_calc_0_ma(pr, dir, q, mhd, mhdabs, mhddiagpress, mhddiagpressabs);
#endif

}

/// MHD stress tensor, with first index up, second index down
/// mhd^dir_j
void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  void mhd_calc_0_em(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
  void mhd_calc_primfield_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);

#if(MAXWELL==GENMAXWELL)
  mhd_calc_0_em(pr, dir, q, mhd, mhdabs);
#elif(MAXWELL==PRIMMAXWELL)
  mhd_calc_primfield_em(pr, dir, geom, q, mhd, mhdabs);
#else
#error No such MAXWELL
#endif


#if(DEBUGNSBH)
  // DEBUG:
  if(geom->i==26 && geom->j==40 && dir==1){
    dualfprintf(fail_file,"INEM1: ucondir=%21.15g %21.15g\n",q->ucon[dir],mhd[3]);
  }
#endif



}


///  MHD stress tensor, with first index up, second index down 
void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);
  void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
  VARSTATIC int j;
  VARSTATIC FTYPE mhdma[NDIM],mhdem[NDIM];
  VARSTATIC FTYPE mhdmaabs[NDIM],mhdemabs[NDIM];
  VARSTATIC FTYPE mhddiagpress[NDIM];
  VARSTATIC FTYPE mhddiagpressabs[NDIM];


  mhd_calc_0_ma(pr, dir, q, mhdma, mhdmaabs, mhddiagpress, mhddiagpressabs);
  mhd_calc_em(pr, dir, geom, q, mhdem, mhdemabs);

  
  // add up MA+EM (no RAD here in T because R kept as separate fluid)
  DLOOPA(j) mhd[j] = (mhdma[j] + mhddiagpress[j]) + mhdem[j];
  if(mhdabs!=NULL) DLOOPA(j) mhdabs[j] = (mhdmaabs[j] + mhddiagpressabs[j]) + mhdemabs[j];

}



///  MHD stress tensor, with first index up, second index down 
void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs)
{
  VARSTATIC int j;
  VARSTATIC FTYPE rho, u, P, w, eta, ptot;

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(DOEVOLVERHO)
  rho = pr[RHO];
#else
  rho = 0.0;
#endif

#if(DOEVOLVEUU)
  u = pr[UU];
  P = q->pressure;
#else
  u = P = 0.0;
#endif

  w = P + rho + u;
  eta = w;
  ptot = P;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j];
  if(mhdabs!=NULL) DLOOPA(j) mhdabs[j] = fabs(mhd[j]);

  if(mhddiagpress!=NULL) DLOOPA(j) mhddiagpress[j] = 0.0;
#if(SPLITPRESSURETERMINFLUXMA==0)
  mhd[dir] += ptot;
  if(mhdabs!=NULL) mhdabs[dir] += fabs(ptot);
#else
  // below equivalent to ptot * delta(dir,j)
  if(mhddiagpress!=NULL) mhddiagpress[dir] = ptot; // NOTEMARK: assumes if here, then mdhdiagpress better not be NULL
  if(mhddiagpressabs!=NULL) mhddiagpressabs[dir] = fabs(ptot);
#endif

}


/// EM part of stress-energy tensor
void mhd_calc_0_em(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  VARSTATIC int j;
  VARSTATIC FTYPE r, u, P, w, bsq, eta, ptot;

  bsq = dot(q->bcon, q->bcov);
  eta = bsq;
  ptot = bsq*0.5;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  //  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j] + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];
  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j] - q->bcon[dir] * q->bcov[j];
  if(mhdabs!=NULL) DLOOPA(j) mhdabs[j] = fabs(mhd[j]);
  mhd[dir] += ptot;
  if(mhdabs!=NULL) mhdabs[dir] += fabs(ptot);

}


/// MHD stress tensor, with first index up, second index down
/// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
/// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
/// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
/// T^dir_\mu
void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhdma, FTYPE *mhdmaabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs);
  void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhdem, FTYPE *mhdemabs);
  VARSTATIC FTYPE mhdma[NDIM];
  VARSTATIC FTYPE mhdem[NDIM];
  VARSTATIC FTYPE mhdmaabs[NDIM];
  VARSTATIC FTYPE mhdemabs[NDIM];
  VARSTATIC int j;
  VARSTATIC FTYPE mhddiagpress[NDIM];
  VARSTATIC FTYPE mhddiagpressabs[NDIM];

  mhd_calc_norestmass_ma(pr, dir, geom, q, mhdma, mhdmaabs, mhddiagpress, mhddiagpressabs);
  mhd_calc_em(pr, dir, geom, q, mhdem, mhdemabs);

  // add up MA and RAD and EM parts
  DLOOPA(j) mhd[j] = (mhdma[j] + mhddiagpress[j]) + mhdem[j];

  if(mhdabs!=NULL) DLOOPA(j) mhdabs[j] = (mhdmaabs[j] + mhddiagpressabs[j]) + mhdemabs[j];

}



/// MHD stress tensor, with first index up, second index down
/// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
/// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
/// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
/// T^dir_\mu
void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs, FTYPE *mhddiagpress, FTYPE *mhddiagpressabs)
{
  VARSTATIC int j;
  VARSTATIC FTYPE rho, u, P, w, bsq, eta, ptot;


  //////////////////////
  //
  // Hydro part
  //
  ///////////////////////

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(DOEVOLVERHO)
  rho = pr[RHO];
#else
  rho = 0.0;
#endif

#if(DOEVOLVEUU)
  u = pr[UU];
  P = q->pressure;
#else
  u = P = 0.0;
#endif

  w = P + rho + u;
  eta = w;
  ptot = P;



  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  // eta u^dir u_j + rho u^dir = (p+u+b^2) u^dir u_j + rho u^dir u_j + rho u^dir
  // = (p+u+b^2) u^dir u_j + rho u^dir (u_j + 1)

  if(mhddiagpress!=NULL) DLOOPA(j) mhddiagpress[j] = 0.0;
  if(mhddiagpressabs!=NULL) DLOOPA(j) mhddiagpressabs[j] = 0.0;

  // T^dir_0
  j=0;
  FTYPE term1,term2;
  term1 = (P+u) * q->ucon[dir] *q->ucov[j];
  term2 = rho * q->ucon[dir] * q->ifremoverestplus1ud0elseud0;
  mhd[j] = term1 + term2 ;
  if(mhdabs!=NULL) mhdabs[j] = fabs(term1) + fabs(term2) ;

  // T^dir_j
  SLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j];
  if(mhdabs!=NULL) SLOOPA(j) mhdabs[j] = fabs(mhd[j]);

#if(SPLITPRESSURETERMINFLUXMA==0)
  // below equivalent to ptot * delta(dir,j)
  mhd[dir] += ptot;
  if(mhdabs!=NULL) mhdabs[dir] += fabs(ptot);
#else
  // below equivalent to ptot * delta(dir,j)
  if(mhddiagpress!=NULL) mhddiagpress[dir] = ptot; // NOTEMARK: assume mhddiagpress!=NULL if here
  if(mhddiagpressabs!=NULL) mhddiagpressabs[dir] = fabs(ptot);
#endif


#if(DEBUGNSBH)
  // DEBUG:
  if(geom->i==26 && geom->j==40 && dir==1){
    dualfprintf(fail_file,"INMA: ucondir=%21.15g %21.15g\n",q->ucon[dir],mhd[3]);
  }
#endif


}



/// MHD stress tensor, with first index up, second index down
/// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
/// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
/// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
/// T^dir_\mu
/// SHOULD NOT USE b^\mu b_\mu here
void mhd_calc_primfield_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs)
{
  int Mcon_calc(FTYPE *pr, struct of_state *q, FTYPE (*Mcon)[NDIM]);
  void ffdestresstensor_dir(int dir, FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE *TEMdir);
  VARSTATIC FTYPE Mcon[NDIM][NDIM];
  VARSTATIC FTYPE TEMdir[NDIM];
  int mu;
  FTYPE Bsq,udotB,Bcon[NDIM],Bcov[NDIM],oneovergammasq;
  int jj,kk;

  //////////////////////
  //
  // Electrommagnetic part
  //
  ///////////////////////
  // gives \dF^{\mu \nu}

#if(0) // old way
  // GODMARK: need full Mcon to get partial TEM -- can one optimize this?
  Mcon_calc(pr, q, Mcon);

  // gives T^dir_\mu
  ffdestresstensor_dir(dir, Mcon, geom, mhd);
#else
  // optimization that still avoids cancellation issue

  // T^\mu_\nu[EM] = b^2 u^\mu u_\nu - b^\mu b_\nu + (b^2/2)\delta^\mu_\nu
  // b^\mu = P^\mu_\nu B^\nu\gamma ; \gamma = -\eta.u where here \eta={-1,0,0,0} so \gamma=u^t
  // b^2 = (B^2 + (udotB)^2)/\gamma^2
  // --> 2 terms cancel
  // T^\mu_\nu = { B^2 u^\mu u_\nu - B^\mu B_\nu - udotB [ B^\mu u_\nu + B_\nu u^\mu] + (1/2)\delta^\mu_\nu[B^2 + udotB^2] }/\gamma^2


  // set B^\mu
  Bcon[TT]=0.0; SLOOPA(jj) Bcon[jj] = pr[B1+jj-1];

  // Get B_\mu = B^i g_{\mu i}
  DLOOPA(jj){
    Bcov[jj] = 0.0;
    SLOOPA(kk) Bcov[jj] += Bcon[kk]*(geom->gcov[GIND(jj,kk)]);
  }
  Bsq=0.0;   SLOOPA(kk) Bsq   += Bcon[kk]*Bcov[kk]; // B^i B_i
  udotB=0.0; SLOOPA(kk) udotB += Bcon[kk]*(q->ucov[kk]); // B^i u_i
  

  // T^dir_\mu, where dir is fixed and mu varies
  DLOOPA(mu) mhd[mu] = (  Bsq*(q->ucon[dir])*(q->ucov[mu]) - Bcon[dir]*Bcov[mu] - udotB*(Bcon[dir]*(q->ucov[mu]) + Bcov[mu]*(q->ucon[dir])) );
  // \delta^dir_\mu term:
  mhd[dir] += 0.5*(Bsq + udotB*udotB);

  oneovergammasq=1.0/((q->ucon[TT])*(q->ucon[TT]));

  DLOOPA(mu) mhd[mu] *= oneovergammasq ;

#endif


  if(mhdabs!=NULL) DLOOPA(mu) mhdabs[mu] = fabs(mhd[mu]);



#if(DEBUGNSBH)
  // DEBUG:
  if(geom->i==26 && geom->j==40 && dir==1){
    dualfprintf(fail_file,"INEM2: ucondir=%21.15g Bcondir=%21.15g (%21.15g %21.15g %21.15g) mhd3=%21.15g\n",q->ucon[dir],Bcon[dir],Bsq,udotB,oneovergammasq,mhd[3]);
    DLOOPA(mu) dualfprintf(fail_file,"mu=%d ucov[mu]=%21.15g Bcov[mu]=%21.15g\n",mu,q->ucov[mu],Bcov[mu]);
  }
#endif


}




/// plus1ud0=(1+q->ucov[TT])
/// avoids non-relativistic velocitiy issue with machine precision, but introduces relativistic limit problem
void compute_1plusud0_old(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  int j,k;
  FTYPE plus1gv00;
  FTYPE AA,BB,alpha;
  FTYPE vcon[NDIM];

  // 3-velocity in coordinate basis
  SLOOPA(j) vcon[j]=q->ucon[j]/q->ucon[TT];

  //  plus1gv00=1.0+geom->gcov[GIND(TT,TT)];
  plus1gv00=geom->gcovpert[TT];

  AA=0.0;
  SLOOPA(j) AA+=2.0*geom->gcov[GIND(TT,j)]*vcon[j];
  SLOOP(j,k) AA+=geom->gcov[GIND(j,k)]*vcon[j]*vcon[k];
  //AA/=geom->gcov[GIND(TT,TT)];
  BB=geom->gcov[GIND(TT,TT)];

  alpha=0.0;
  SLOOPA(j) alpha+=geom->gcov[GIND(j,TT)]*q->ucon[j];

  //  *plus1ud0=(plus1gv00+(2.0*alpha+alpha*alpha)*(1.0+AA)+AA)/((1.0+alpha)*(1.0+AA)+sqrt(-geom->gcov[GIND(TT,TT)]*(1.0+AA)));

  *plus1ud0=(plus1gv00*BB+(2.0*alpha+alpha*alpha)*(BB+AA)+AA)/((1.0+alpha)*(BB+AA)+BB*sqrt(-(BB+AA)));

}

/// plus1ud0=(1+q->ucov[TT])
/// avoids both non-rel and rel limit issues with machine precision
/// GODMARK: Slow, can one optmize this somehow?
/// This computes 1+u_t , which for nonrelativistic cases is ~0 .  If computed as 1+u_t, then residual will be large error if small residual.
/// old2 is newer than old
void compute_1plusud0_general(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  VARSTATIC int j,k;
  VARSTATIC FTYPE plus1gv00;
  VARSTATIC FTYPE vsq,gvtt,alpha;
  VARSTATIC FTYPE vcon[NDIM];
  VARSTATIC FTYPE uu0;


  // 3-velocity in coordinate basis
  SLOOPA(j) vcon[j]=q->ucon[j]/q->ucon[TT];

  //  plus1gv00=1.0+geom->gcov[GIND(TT,TT)];
  plus1gv00=geom->gcovpert[TT];

  vsq=geom->gcovpert[TT];
  SLOOPA(j) vsq+=2.0*geom->gcov[GIND(TT,j)]*vcon[j];
  SLOOP(j,k) vsq+=geom->gcov[GIND(j,k)]*vcon[j]*vcon[k];

  gvtt=geom->gcov[GIND(TT,TT)];

  alpha=0.0;
  SLOOPA(j) alpha+=geom->gcov[GIND(j,TT)]*q->ucon[j];

  uu0 = q->ucon[TT];

  *plus1ud0=alpha + ((1.0-gvtt)*plus1gv00 - uu0*uu0*vsq*gvtt*gvtt)/(1.0-gvtt*uu0);

  //  dualfprintf(fail_file,"%g %g %g : wrong: %g\n",alpha,plus1gv00,gvtt,*plus1ud0);

  //  *plus1ud0=1.0+q->ucov[TT];

  //  dualfprintf(fail_file,"right: %g\n",*plus1ud0);

}


/// latest method to compute 1+u_t without catastrophic cancellation
void compute_1plusud0_rel4vel(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  VARSTATIC int j,k;
  VARSTATIC FTYPE ud0tilde,alpha,alphasq;
  VARSTATIC FTYPE gamma,qsq;

  // from q->others[] only generated for WHICHVEL=REL4VEL
  gamma=q->others[OTHERGAMMA];
  qsq=q->others[OTHERQSQ];
  

  // \tilde{u}_t = \tilde{u}^i g_{ti} since \tilde{u}^t=0
  ud0tilde = 0.0;
  SLOOPA(j) ud0tilde += pr[U1+j-1]*(geom->gcov[GIND(TT,j)]);

  alpha = geom->alphalapse;
  alphasq = alpha*alpha;

  *plus1ud0 = ud0tilde + (geom->gcovpert[TT] - alphasq*( (geom->betasqoalphasq) + qsq))/(1.0+gamma*alpha);

}


/// add in source terms to equations of motion
/// ui and dUriemann in UEVOLVE form
/// assume q(pr) so consistent, but p or ui don't yet account for dUriemann!
int source(FTYPE *pi, FTYPE *pr, FTYPE *pf, int *didreturnpf, int *eomtype, struct of_geom *ptrgeom, struct of_state *q, FTYPE *ui, FTYPE *uf, FTYPE *CUf, FTYPE *CUimp, FTYPE dissmeasure, FTYPE *dUriemann, FTYPE (*dUcomp)[NPR], FTYPE *dU)
{
  //  double (*)[8]
  VARSTATIC int i,j,sc;
  VARSTATIC int pl,pliter;
  int source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU);
  VARSTATIC FTYPE dUother[NPR];
  VARSTATIC FTYPE Ugeomfreei[NPR],Ugeomfreef[NPR];


  ////
  // initialize source terms to be zero
  ////
  PLOOP(pliter,pl){
    SCLOOP(sc) dUcomp[sc][pl] = 0.;
    dU[pl] = 0.;
  }

#if(SPLITNPR)
#if((WHICHEOM==WITHNOGDET)&&(NOGDETB1>0)||(NOGDETB2>0)||(NOGDETB3>0))
#error Not correct SPLITNPR use#1
#endif
  // if only doing B1-B3, then assume no source term for magnetic field
  if(nprlist[nprstart]==B1 && nprlist[nprend]==B3) return(0);
#endif



  ////////////////
  //
  // First get geometry source (already does contain geometry prefactor term)
  // KORALTODO: Should I use dUriemann to already update p(U) and q(p) for more accurate geometry?
  source_conn(pr,ptrgeom, q, dUcomp[GEOMSOURCE]);


  ////////////////
  //
  // Second get physics source terms (using dUother for constraint)
  // get update pre- additional physics
  // Assumes all quantities are at cell centers
  PLOOP(pliter,pl){
    dUother[pl] = (dUriemann[pl] + dUcomp[GEOMSOURCE][pl])*(ptrgeom->IEOMFUNCNOSINGMAC(pl)); // remove geometry factor for dUother
    // ui is current timestep's ui.
    Ugeomfreei[pl] = ui[pl]*(ptrgeom->IEOMFUNCNOSINGMAC(pl)); // expect ui to be UEVOLVE form
    // uf is "previous" timestep's uf, which is not generally ui for arbitrary RK methods, and should not include dUother, etc.
    Ugeomfreef[pl] = uf[pl]*(ptrgeom->IEOMFUNCNOSINGMAC(pl)); // expect uf to be UEVOLVE form
  }

  if(FLUXB==FLUXCTSTAG){
    // if staggered field, then ui and uf are staggered field locations, but for *physics* processes, we want to assume all at same location, so send in physical version
    // this overwrites above assignments
    PLOOPBONLY(pl){
      // assume already got field update in advance_standard() [as opposed to advance_standard_orig()] and no geometry for field as required for that method.
      // but want to be able to have Ui by itself mean no changes, so that's pr. Uf is only used in some RK methods, that can be pf.  But then want dUother to be so that when using IFSET() with full dt that get pf
      //      Ugeomfreei[pl]=pi[pl];  // unnecessary
      //      Ugeomfreef[pl]=pf[pl]; // No, should stay as from uf as *previous* timestep's Uf, so that RK3/4 are consistent with below.
      // KORALTODO SUPERGODMARK: This means need advance.c source() to have uf as previous field uf, not updated uf from dUriemann *and* not tempucum for finalstep=1
      dUother[pl]=dUfromUFSET(CUf,dt,Ugeomfreei[pl],Ugeomfreef[pl],pf[pl]);
      // Also, update "guess" pr with new field, since pr is used as default guess and want field to be correct.  So pr is no longer what computed flux.
      //      pr[pl]=pf[pl]; // only for implicit schemes, and handle this inside koral_source_rad_implicit() now.
    }
    // now sourcephysics() call will have all CENT quantities
  }

  sourcephysics(pi, pr, pf, didreturnpf, eomtype, ptrgeom, q, Ugeomfreei, Ugeomfreef, CUf, CUimp, dissmeasure, dUother, dUcomp);

  //////////////////
  //
  // Third, deal with equation of motion factors that don't depend upon any additional physics
  // don't add geometry prefactor onto geometry source since already has it -- only add to additional physics terms
  SCPHYSICSLOOP(sc) PLOOP(pliter,pl) dUcomp[sc][pl] *= ptrgeom->EOMFUNCMAC(pl);


  //////////////////
  //
  // Fourth, compute total since that's all the evolution cares about (comp left just for diagnostics).
  SCLOOP(sc) PLOOP(pliter,pl) dU[pl]+=dUcomp[sc][pl];
  


  //  done! 
  return (0);
}



///  returns b^2 (i.e., twice magnetic pressure) 
int bsq_calc_general(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  VARSTATIC struct of_state q;

  MYFUN(get_state(pr, ptrgeom, &q) ,"phys.c:bsq_calc()", "get_state() dir=0", 1);
  *bsq = dot(q.bcon, q.bcov);
  return (0);
}

///  returns b^2 (i.e., twice magnetic pressure) 
int bsq_calc_fromq_general(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq)
{

  *bsq = dot(q->bcon, q->bcov);

  return (0);
}

///  returns b^2 (i.e., twice magnetic pressure) 
int bsq_calc_rel4vel(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  VARSTATIC struct of_state q;
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int bsq_calc_fromq_rel4vel(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);
  

  MYFUN(get_state_uconucovonly(pr, ptrgeom, &q) ,"phys.c:bsq_calc()", "get_state_uconucovonly() dir=0", 1);

  bsq_calc_fromq_rel4vel(pr, ptrgeom, &q, bsq);

  return (0);
}


int bsq_calc_fromq_rel4vel(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq)
{
  FTYPE Bsq,udotB,Bcon[NDIM],Bcov[NDIM],oneovergammasq;
  int jj,kk;

  // set B^\mu
  Bcon[TT]=0.0; SLOOPA(jj) Bcon[jj] = pr[B1+jj-1];

  // Get B_\mu = B^i g_{\mu i}
  DLOOPA(jj){
    Bcov[jj] = 0.0;
    SLOOPA(kk) Bcov[jj] += Bcon[kk]*(ptrgeom->gcov[GIND(jj,kk)]);
  }
  Bsq=0.0;   SLOOPA(kk) Bsq   += Bcon[kk]*Bcov[kk]; // B^i B_i
  udotB=0.0; SLOOPA(kk) udotB += Bcon[kk]*(q->ucov[kk]); // B^i u_i
  
  *bsq = (Bsq + udotB*udotB);

  oneovergammasq=1.0/((q->ucon[TT])*(q->ucon[TT]));

  *bsq *= oneovergammasq;

  return(0);
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

void lowerf(FTYPE *fcon, struct of_geom *geom, FTYPE *fcov)
{
  VARSTATIC int j,k;
  VARSTATIC int jp,kp;
  VARSTATIC FTYPE myfcon[NDIM][NDIM],myfcov[NDIM][NDIM];


  myfcon[0][0]=myfcon[1][1]=myfcon[2][2]=myfcon[3][3]=0;
  myfcon[0][1]=fcon[0];
  myfcon[0][2]=fcon[1];
  myfcon[0][3]=fcon[2];
  myfcon[1][2]=fcon[3];
  myfcon[1][3]=fcon[4];
  myfcon[2][3]=fcon[5];
  //
  myfcon[1][0]=-fcon[0];
  myfcon[2][0]=-fcon[1];
  myfcon[3][0]=-fcon[2];
  myfcon[2][1]=-fcon[3];
  myfcon[3][1]=-fcon[4];
  myfcon[3][2]=-fcon[5];
  
  DLOOP(j,k){
    myfcov[j][k]=0;
    for(jp=0;jp<NDIM;jp++) for(kp=0;kp<NDIM;kp++){
        myfcov[j][k]+=myfcon[jp][kp]*geom->gcov[GIND(j,jp)]*geom->gcov[GIND(k,kp)];
      }
  }
  fcov[0]=myfcov[0][1];
  fcov[1]=myfcov[0][2];
  fcov[2]=myfcov[0][3];
  fcov[3]=myfcov[1][2];
  fcov[4]=myfcov[1][3];
  fcov[5]=myfcov[2][3];

  return ;
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


/// find ucon, ucov, bcon, bcov from primitive variables
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int get_state(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

  get_state_bconbcovonly(pr, ptrgeom, q);

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);

  return (0);
}


/// find ucon, ucov, bcon, bcov from primitive variables
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int get_state_radonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif
  return (0);
}

/// find ucon, ucov, bcon, bcov from primitive variables
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int get_state_norad_part1(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

  get_state_bconbcovonly(pr, ptrgeom, q);

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  return (0);
}

/// find ucon, ucov, bcon, bcov from primitive variables
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int get_state_norad_part2(int needentropy, FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);


  get_state_thermodynamics(needentropy, ptrgeom, pr, q);

  return (0);
}



/// all get_state() things except the field quantities
int get_state_nofield(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

  get_state_thermodynamics(1,ptrgeom, pr, q);

  return (0);
}


/// used to check inversion, which has consistent pressure used to get things as functions of \chi instead of u in case of using jon's inversion
int get_stateforcheckinversion(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics_forcheckinversion(struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

  get_state_bconbcovonly(pr, ptrgeom, q);

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics_forcheckinversion(ptrgeom, pr, q);

  return (0);
}




/// find ucon, ucov, bcon, bcov from primitive variables
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure,g,e[pl],prim,Blower}
int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);
  int get_state_prims(FTYPE *pr, struct of_state *q);
  int get_state_geom(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_Blower(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_vcon_gdetBcon_overut(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


#if(MERGEDC2EA2CMETHOD)
  // only needed if doing merged method
  get_state_geom(pr, ptrgeom, q);

  get_state_prims(pr, q);

  get_state_Blower(pr, ptrgeom, q);
#endif

  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

#if(MERGEDC2EA2CMETHOD)
  get_state_vcon_gdetBcon_overut(pr, ptrgeom, q);
#endif

#if(COMPUTE4FIELDforFLUX)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);


  return (0);
}


///  find ucon, ucov, bcon, bcov from primitive variables 
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int pureget_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif


  // below needed for dUtodt() and compute_dt_fromsource()
#if(COMPUTE4FIELDatALL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);

  return (0);
}

///  find ucon, ucov, bcon, bcov from primitive variables 
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int pureget_stateforinterpline(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

#if(COMPUTE4FIELDatALL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);

  return (0);
}

///  find ucon, ucov, bcon, bcov from primitive variables 
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int pureget_stateforglobalwavespeeds(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

#if(COMPUTE4FIELDatALL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);

  return (0);
}

///  find ucon, ucov, bcon, bcov from primitive variables 
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure,g,e[pl],prim,Blower}
int pureget_stateforfluxcalcorsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);
  int get_state_prims(FTYPE *pr, struct of_state *q);
  int get_state_geom(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_Blower(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_vcon_gdetBcon_overut(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


#if(MERGEDC2EA2CMETHOD)
  // only needed if doing merged method
  get_state_geom(pr, ptrgeom, q);

  get_state_prims(pr, q);

  get_state_Blower(pr, ptrgeom, q);
#endif

  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

#if(MERGEDC2EA2CMETHOD)
  get_state_vcon_gdetBcon_overut(pr, ptrgeom, q);
#endif

#if(COMPUTE4FIELDatALL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);


  return (0);
}



///  find ucon, ucov, bcon, bcov from primitive variables 
/// when calling get_state, users of this function expect to get q->{ucon,ucov,bcon,bcov,pressure}
int get_stateforUdiss(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_thermodynamics(int needentropy,struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
  int bsq_calc_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(EOMRADTYPE!=EOMRADNONE)
  get_state_uradconuradcovonly(pr, ptrgeom, q);
#endif

#if(COMPUTE4FIELDatALL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif

  bsq_calc_fromq(pr,ptrgeom,q,&(q->bsq));

  get_state_thermodynamics(1,ptrgeom, pr, q);


  return (0);
}



///  find ucon, ucov, bcon, bcov from primitive variables 
int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{

  // b^\mu
  bcon_calc(pr, q->ucon, q->ucov, q->bcon);
  // b^\mu
  lower_vec(q->bcon, ptrgeom, q->bcov);

  return (0);
}


/// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  void compute_1plusud0(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])


  // u^\mu
  MYFUN(ucon_calc(pr, ptrgeom, q->ucon,q->others) ,"phys.c:get_state()", "ucon_calc()", 1);


  // u_\mu
  lower_vec(q->ucon, ptrgeom, q->ucov);

#if(REMOVERESTMASSFROMUU==1 || REMOVERESTMASSFROMUU==2)

#if(UD0PLUS1FIX==0)
  q->ifremoverestplus1ud0elseud0=1.0+(q->ucov[TT]);
#else
  compute_1plusud0(pr,ptrgeom,q,&(q->ifremoverestplus1ud0elseud0)); // plus1ud0=(1+q->ucov[TT])
#endif

#else
  q->ifremoverestplus1ud0elseud0=q->ucov[TT];
#endif

  return (0);
}




/// get v^i and \detg B^i used for EMF computations for FLUXB==FLUXCTSTAG
int get_state_vcon_gdetBcon_overut(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int jj;
  FTYPE overut;
  
  overut=1.0/q->ucon[TT];

#if(MERGEDC2EA2CMETHOD)
  q->overut=overut;

  q->vcon[TT]=1.0; SLOOPA(jj) q->vcon[jj]=q->ucon[jj]*overut;

  q->gdetBcon[TT]=0.0; SLOOPA(jj) q->gdetBcon[jj]=(q->prim[B1-1+jj])*(q->EOMFUNCMAC(B1-1+jj));
#endif

  return (0);
}



/// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_Blower(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int jj;
  FTYPE Bcon[NDIM];
  
  Bcon[TT]=0.0; SLOOPA(jj) Bcon[jj]=pr[B1-1+jj];

#if(MERGEDC2EA2CMETHOD)
  // get B_\mu
  lower_vec(Bcon, ptrgeom, q->Blower);
#endif

  return (0);
}



/// separate function for getting thermodynamical quantities
int get_state_thermodynamics(int needentropy, struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q)
{

  // does assume ptrgeom or something related set indices for advanced EOSs
  q->pressure = pressure_rho0_u_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pr[RHO],pr[UU]);

#if(DOENTROPY!=DONOENTROPY)
  if(needentropy) entropy_calc(ptrgeom,pr,&(q->entropy));
#endif

  return (0);

}


/// separate function for getting thermodynamical quantities
int get_state_thermodynamics_forcheckinversion(struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q)
{

  // does assume ptrgeom or something related set indices for advanced EOSs
  q->pressure = pressure_rho0_u_simple_forcheckinversion(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pr[RHO],pr[UU]);

#if(DOENTROPY!=DONOENTROPY)
  entropy_calc_forcheckinversion(ptrgeom,pr,&(q->entropy));
#endif

  return (0);

}

/// separate function for getting primitives
int get_state_prims(FTYPE *pr, struct of_state *q)
{
  int pl,pliter;

#if(MERGEDC2EA2CMETHOD)
  // just copies  
  PALLLOOP(pl) q->prim[pl] = pr[pl];
#endif

  return (0);

}

/// separate function for getting geometry
int get_state_geom(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int pl,pliter;

#if(MERGEDC2EA2CMETHOD)
  q->gdet=ptrgeom->gdet;
  PALLLOOP(pl) q->EOMFUNCMAC(pl)=ptrgeom->EOMFUNCMAC(pl); //  note that this used even if WITHEOM==WITHGDET
#endif

  return (0);
}


///ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p
int invertentropyflux_calc(struct of_geom *ptrgeom, FTYPE entropyflux,int dir, struct of_state *q, FTYPE *pr)
{
  VARSTATIC FTYPE entropy;

  // get entropy [entropy can be any value]
  entropy=entropyflux/q->ucon[dir];
  ufromentropy_calc(ptrgeom, entropy, pr);

  return(0);
}

/// u from entropy (uses pr[RHO])
/// wrapper
int ufromentropy_calc(struct of_geom *ptrgeom, FTYPE entropy, FTYPE *pr)
{

  // entropy version of ie
  pr[ENTROPY]=compute_u_from_entropy_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pr[RHO],entropy);

  return(0);
}


/// find u^\mu from \tilde{u}^\mu
/// fill full 4-vector
/// OPTMARK: Before storing \beta, this was the most expensive function over the entire code [determined by avoiding inlining]
int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{
  VARSTATIC FTYPE gamma;
  VARSTATIC FTYPE qsq;
  VARSTATIC int j ;
  int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);


  MYFUN(gamma_calc_fromuconrel(uconrel,geom,&gamma,&qsq),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.c",1);

  // aux results
  others[OTHERGAMMA]=gamma;
  others[OTHERQSQ]=qsq;

  //  alpha = 1./sqrt(-geom->gcon[GIND(TT,TT)]) ;
  ucon[TT] = gamma/(geom->alphalapse) ;
  // stored beta as well since otherwise have to use gcon and expensive to look that up from memory if cache-miss
  SLOOPA(j) ucon[j] = uconrel[j] - ucon[TT]*(geom->beta[j]);

  // hence v^j = uconrel^j/u^t - beta^j

  return(0) ;
}


/// find \tilde{u}^\mu from u^\mu
/// only fill spatial parts so can feed in 3-vector
int uconrel(FTYPE *ucon, FTYPE *uconrel, struct of_geom *geom)
{
  VARSTATIC int j ;

  // stored beta as well since otherwise have to use gcon and expensive to look that up from memory if cache-miss
  SLOOPA(j) uconrel[j] = ucon[j] + ucon[TT]*(geom->beta[j]);

  // hence v^j = uconrel^j/u^t - beta^j

  return(0) ;
}


///  find contravariant four-velocity from the relative 4 velocity 
int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{
  VARSTATIC FTYPE uconrel[NDIM];
  VARSTATIC int j;
  int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon, FTYPE *others);

  uconrel[0]=0;
  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];



  //  SLOOPA(j) dualfprintf(fail_file,"ucon_calc_rel4vel: uconrel[%d]=%21.15g\n",j,uconrel[j]); // CHANGINGMARK

  MYFUN(ucon_calc_rel4vel_fromuconrel(uconrel, geom, ucon, others),"ucon_calc_rel4vel: ucon_calc_rel4vel_fromuconrel failed\n","phys.c",1);
  

  return(0) ;
}

int ucon_calc_whichvel(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{
  if(whichvel==VEL4){
    return(ucon_calc_4vel(pr,geom,ucon,others));
  }
  else if(whichvel==VEL3){
    return(ucon_calc_3vel(pr,geom,ucon,others));
  }
  else if(whichvel==VELREL4){
    return(ucon_calc_rel4vel(pr,geom,ucon,others));
  }
  else{
    dualfprintf(fail_file,"No such whichvel=%d for ucon_calc_whichvel()\n",whichvel);
    myexit(384852);
  }

  return(1); // shouldn't get here

}


///  find gamma-factor wrt normal observer 
int gamma_calc(FTYPE *pr, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq)
{
  VARSTATIC FTYPE uconrel[NDIM];
  VARSTATIC int j;
  int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);

#if(WHICHVEL!=VELREL4)
  dualfprintf(fail_file,"gamma_calc() designed for WHICHVEL=VELREL4\n");
  myexit(77);
#endif

  // assumes input pr is WHICHVEL=VELREL4
  uconrel[0]=0;
  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];

  //  SLOOPA(j) dualfprintf(fail_file,"gamma_calc: uconrel[%d]=%21.15g\n",j,uconrel[j]); // CHANGINGMARK


  // get gamma
  MYFUN(gamma_calc_fromuconrel(uconrel, geom, gamma, qsq),"gamma_calc: gamma_calc_fromuconrel failed\n","phys.c",1);

  return(0);
}


///  find gamma-factor wrt normal observer 
/// This function and qsq_calc() have about the same cache miss amount now
int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq)
{
  int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq);


  // get q^2 = \tilde{u}^i \tilde{u}^j g_{ij}
  qsq_calc(uconrel,geom,qsq);


#if(JONCHECKS && PRODUCTION==0)
  VARSTATIC int j,k;
  if(*qsq<0.0){
    if(*qsq>-NUMEPSILON*100){ // then assume not just machine precision
      *qsq=0.0;
    }
    else{
      if(failed!=-1){
        dualfprintf(fail_file,"gamma_calc failed: i=%d j=%d k=%d qsq=%21.15g\n",geom->i,geom->j,geom->k,*qsq);
        SLOOPA(j) dualfprintf(fail_file,"uconrel[%d]=%21.15g\n",j,uconrel[j]);
        DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,geom->gcov[GIND(j,k)]);
      }
      if (fail(geom->i,geom->j,geom->k,geom->p,FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
  }
#else
  // minimal check -- forces qsq=0 if qsq<0.0
  *qsq = (*qsq)*((*qsq)>=0.0);
#endif

  *gamma = sqrt(1. + (*qsq)) ;

  //  SLOOPA(j) dualfprintf(fail_file,"uconrel[%d]=%21.15g\n",j,uconrel[j]);
  //  dualfprintf(fail_file,"qsq=%21.15g, gamma=%21.15g\n",*qsq,*gamma); // CHANGINGMARK

  return(0) ;
}


/// get \tilde{u}^i \tilde{u}^j g_{ij}
/// OPTMARK: This is currently the most expensive function [disabled inlining]
int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq)
{
  //  VARSTATIC int j,k;
  
  *qsq =
    + geom->gcov[GIND(1,1)]*uconrel[1]*uconrel[1]
    + geom->gcov[GIND(2,2)]*uconrel[2]*uconrel[2]
    + geom->gcov[GIND(3,3)]*uconrel[3]*uconrel[3]
#if(DOMIXTHETAPHI)
    + 2.*(geom->gcov[GIND(1,2)]*uconrel[1]*uconrel[2]
          + geom->gcov[GIND(1,3)]*uconrel[1]*uconrel[3]
          + geom->gcov[GIND(2,3)]*uconrel[2]*uconrel[3])
#else
    + 2.*(geom->gcov[GIND(1,2)]*uconrel[1]*uconrel[2]
          + geom->gcov[GIND(1,3)]*uconrel[1]*uconrel[3])
#endif
    ;

  //  DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g uconrel[%d]=%21.15g uconrel[%d]=%21.15g\n",j,k,geom->gcov[GIND(j,k)],j,uconrel[j],k,uconrel[k]); // CHANGINGMARK

  return(0) ;
}






///  find contravariant four-velocity 
///int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{
  VARSTATIC FTYPE negdiscr;
  VARSTATIC FTYPE velterm;
  VARSTATIC FTYPE vcon[NDIM];
  // debug stuff
  VARSTATIC int j,k;
  VARSTATIC FTYPE V[NDIM],X[NDIM];
  //
  int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm);



  SLOOPA(j) vcon[j] = pr[U1+j-1];

  get_3velterm(vcon, geom, &velterm);

  negdiscr = geom->gcov[GIND(0,0)] + velterm ;


#if(JONCHECKS && WHICHVEL==VEL3)  
  uttdiscr=-negdiscr;
#endif

  if (negdiscr > 0.) {
#if(JONCHECKS)
    if(failed!=-1){
      if(whocalleducon==0){
        // then report on disc
        dualfprintf(fail_file,"negdisc=%21.15g, should be negative\n",negdiscr);
        for(k=U1;k<=U3;k++){
          dualfprintf(fail_file,"uconfailed on pr[%d]=%21.15g\n",k,pr[k]);
        }
        bl_coord_ijk_2(geom->i,geom->j,geom->k,geom->p,X, V);
        dualfprintf(fail_file,"i=%d j=%d k=%d pcurr=%d\nx1=%21.15g x2=%21.15g x3=%21.15g \nV1=%21.15g V2=%21.15g V3=%21.15g \ng=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p,X[1],X[2],X[3],V[1],V[2],V[3],geom->gdet);
        dualfprintf(fail_file,"\ngcon\n");
        dualfprintf(fail_file,"{");
        for(j=0;j<NDIM;j++){
          dualfprintf(fail_file,"{");
          for(k=0;k<NDIM;k++){
            dualfprintf(fail_file,"%21.15g",geom->gcon[GIND(j,k)]);
            if(k!=NDIM-1) dualfprintf(fail_file," , ");
          }
          dualfprintf(fail_file,"}"); 
          if(j!=NDIM-1) dualfprintf(fail_file," , ");
        }
        dualfprintf(fail_file,"}");
        dualfprintf(fail_file,"\ngcov\n");
        dualfprintf(fail_file,"{");
        for(j=0;j<NDIM;j++){
          dualfprintf(fail_file,"{");
          for(k=0;k<NDIM;k++){
            dualfprintf(fail_file,"%21.15g",geom->gcov[GIND(j,k)]);
            if(k!=NDIM-1) dualfprintf(fail_file," , ");
          }
          dualfprintf(fail_file,"}"); 
          if(j!=NDIM-1) dualfprintf(fail_file," , ");
        }
        dualfprintf(fail_file,"}");
      }
    }
#endif
    if (fail(geom->i,geom->j,geom->k,geom->p,FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  ucon[TT] = 1. / sqrt(-negdiscr);
  SLOOPA(j) ucon[j]=vcon[j]*ucon[TT];

  return (0);
}


/// 2v^i g_{it} + v^i v^j g_{ij}
int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm)
{
  VARSTATIC int j;


  *velterm =
    + geom->gcov[GIND(1,1)] * vcon[1] * vcon[1]
    + geom->gcov[GIND(2,2)] * vcon[2] * vcon[2]
    + geom->gcov[GIND(3,3)] * vcon[3] * vcon[3]
    + 2. * (geom->gcov[GIND(0,1)]* vcon[1]
            + geom->gcov[GIND(0,2)] * vcon[2]
            + geom->gcov[GIND(0,3)] * vcon[3]
            + geom->gcov[GIND(1,2)] * vcon[1] * vcon[2]
            + geom->gcov[GIND(1,3)] * vcon[1] * vcon[3]
            + geom->gcov[GIND(2,3)] * vcon[2] * vcon[3]
            );
  return(0);
}






/// for FFDE, check to make sure 3-velocity is good, otherwise will have to limit it
/// Bcon and vcon are code versions
/// returns code primitives such as WHICHVEL for velocity (although coordinates will still be whatever geom was)
int limit_3vel_ffde(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr)
{
  int dualf_calc(FTYPE *Bcon, FTYPE *vcon, FTYPE (*dualffull)[NDIM]);
  void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM]);
  FTYPE dualf[NDIM][NDIM], T[NDIM][NDIM];
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  FTYPE U[NPR];

#if(EOMTYPE!=EOMFFDE && EOMTYPE!=EOMFFDE2)
  dualfprintf(fail_file,"Only call limit_3vel_ffde() if doing EOMTYPE==EOMFFDE || EOMTYPE==EOMFFDE2\n");
  myexit(346983463);
#endif


  ////////////////////////
  //
  // get \dF^{\mu\nu}
  //
  ////////////////////////
  dualf_calc(Bcon, vcon, dualf);

  //////////////////////
  //
  // get T^\mu_\nu
  //
  //////////////////////
  ffdestresstensor(dualf, geom, T);


  ///////////////////
  //
  // setup conserved quantities
  //
  ///////////////////
  U[RHO] = 0;
  U[UU] = T[0][0]; // T^t_t
  U[U1] = T[0][1]; // T^t_x
  U[U2] = T[0][2]; // T^t_y
  U[U3] = T[0][3]; // T^t_z
  U[B1] = Bcon[1];
  U[B2] = Bcon[2];
  U[B3] = Bcon[3];



  ///////////////////
  //
  // filter through inversion so v^i is limited to have Lorentz factor <=GAMMAMAX
  //
  ///////////////////


  //  Utoprim_ffde(U, geom, pr);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int finalstep=1; // doesn't matter for ffde
  int eomtype=EOMDEFAULT;
  FTYPE dissmeasure=-1.0; // assume trying energy ok
  int whichcap=CAPTYPEBASIC;
  int whichmethod=MODEDEFAULT;
  int modprim=0;
  int checkoninversiongas=CHECKONINVERSION;
  int checkoninversionrad=CHECKONINVERSIONRAD;
  MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep,&eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM,UNOTHING,U, NULL, geom, dissmeasure, pr, pr,&newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
  //  nstroke+=newtonstats.nstroke; newtonstats.nstroke=newtonstats.lntries=0;



  return(0);

}




///  find contravariant time component of four-velocity from the 4velocity (3 terms)
int ucon_calc_4vel_bothut(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *ucon2, FTYPE *others)
{
  VARSTATIC int i,j,k;
  VARSTATIC FTYPE AA,BB,CCM1 ;
  VARSTATIC FTYPE discr ;
  VARSTATIC FTYPE bsq,X[NDIM] ;
  VARSTATIC FTYPE sqrtdiscr;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  ucon2[1] = pr[U1] ;
  ucon2[2] = pr[U2] ;
  ucon2[3] = pr[U3] ;

  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr < 0.) {
    /*
      dualfprintf(fail_file,"failure") ;
      ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
      ucon[TT] = -BB/(2.*AA) ;
    */
    dualfprintf(fail_file,"failure: spacelike four-velocity %21.15g\n",
                discr) ;
    dualfprintf(fail_file,"i=%d j=%d k=%d p=%d\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p) ;
    coord(geom->i,geom->j,geom->k,geom->p,X);
    dualfprintf(fail_file,"%21.15g %21.15g %21.15g\n",X[1],X[2],X[3]) ;

    for(k=0;k<NPR;k++) dualfprintf(fail_file,"%d %21.15g\n",k,pr[k]) ;
    // GODMARK -- why did we have failed=1?
    //  failed=1;
    return(1);
  }

  sqrtdiscr=sqrt(discr);

  ucon[TT] = (-BB - sqrtdiscr)/(2.*AA) ; // standard solution always valid outside ergosphere
  ucon2[TT] = (-BB + sqrtdiscr)/(2.*AA) ; // can be a valid solution within ergosphere

  return(0) ;
}


///  find contravariant time component of four-velocity from the 4velocity (3 terms)
int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{
  VARSTATIC int i,j,k;
  VARSTATIC FTYPE AA,BB,CCM1 ;
  VARSTATIC FTYPE discr ;
  VARSTATIC FTYPE bsq,X[NDIM] ;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr < 0.) {
    /*
      dualfprintf(fail_file,"failure\n") ;
      ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
      ucon[TT] = -BB/(2.*AA) ;
    */
    dualfprintf(fail_file,"failure: spacelike four-velocity %21.15g\n",
                discr) ;
    dualfprintf(fail_file,"i=%d j=%d k=%d p=%d\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p) ;
    coord(geom->i,geom->j,geom->k,geom->p,X);
    dualfprintf(fail_file,"%21.15g %21.15g %21.15g\n",X[1],X[2],X[3]) ;

    for(k=0;k<NPR;k++) dualfprintf(fail_file,"%d %21.15g\n",k,pr[k]) ;
    // GODMARK -- why did we have failed=1?
    //  failed=1;
    return(1);
  }

  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;

  return(0) ;
}


int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr)
{
  VARSTATIC FTYPE CC;

  // g_{tt}
  *AA = geom->gcov[GIND(TT,TT)] ;

  // 2 u^i g_{it}
  *BB = 2.*(geom->gcov[GIND(TT,1)]*ucon[1] +
            geom->gcov[GIND(TT,2)]*ucon[2] +
            geom->gcov[GIND(TT,3)]*ucon[3]) ;

  // u^i u^j g_{ij}
  *CCM1 = geom->gcov[GIND(1,1)]*ucon[1]*ucon[1] +
    geom->gcov[GIND(2,2)]*ucon[2]*ucon[2] +
    geom->gcov[GIND(3,3)]*ucon[3]*ucon[3] +
    2.*(geom->gcov[GIND(1,2)]*ucon[1]*ucon[2] +
        geom->gcov[GIND(1,3)]*ucon[1]*ucon[3] +
        geom->gcov[GIND(2,3)]*ucon[2]*ucon[3]) ;

  // 1 + u^i u^j g_{ij}
  CC = *CCM1 + 1.0;
  
  // B^2 - 4AC
  *discr = (*BB)*(*BB) - 4.*(*AA)*CC ;

  return(0);
}


///  find contravariant time component of four-velocity from the 4velocity (3 terms)
int ucon_calc_nonrel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others)
{

  // GODMARK: this isn't really right: neglects kinetic energy term.  Need to really re-write T^\mu_\nu

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;
  ucon[TT] = 1.0;

  return(0) ;
}


FTYPE taper_func(FTYPE R,FTYPE rin)
{

  if(R <= rin)
    return(0.) ;
  else
    return(1. - sqrt(rin/R)) ;

}





/// used Mathematica's MinimumChangePermutations and Signature
/// updown = 0 : down
/// updown = 1 : up
FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda)
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
      if(updown==1) return(-1.0/detg*lc4sign); // up
      else if(updown==0) return(detg*lc4sign); // down
    }
  }
  // if didn't get here, then 0
  return(0.0);
}

/// Compute faraday
/// assumes b and u are inputted as bcov&ucov for F^{\mu\nu} and bcon&ucon for F_{\mu\nu}
void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE (*faraday)[NDIM])
{
  int nu,mu,kappa,lambda;

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
      faraday[mu][nu]=0.0;
      for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
          faraday[mu][nu]+=lc4(which,geom->gdet,mu,nu,kappa,lambda)*u[kappa]*b[lambda];
        }
    }

}




/// assumes below get inlined
/// much faster than macro using ? :
/// super slow for get_geometry()'s sign() call!
FTYPE sign_bad(FTYPE a)
{
  if(a>0.0) return 1.0;
  else if(a<0.0) return -1.0;
  else return 0.0;
}


/// not any faster than above (except when used alot)
FTYPE sign_func(FTYPE a)
{
#if(SUPERLONGDOUBLE)
  return(sign(a));
  // no such function
#else
  return(copysign(1.0,a));
#endif
}

/// much faster than macro using ? :
/// Generally avoid using, instead do:
/// igdet = sign(geom.gdet)/(fabs(geom.gdet)+SMALL)
FTYPE signavoidzero(FTYPE a)
{
  if(a>0.0) return 1.0;
  else if(a<0.0) return -1.0;
  else return 1.0; // arbitrarily choose 1.0
}



#ifndef WIN32

#ifndef max
FTYPE max(FTYPE a, FTYPE b)
{
  if(a>b) return a;
  else return b;
  // if equal, then above is fine
}
#endif

#ifndef min
FTYPE min(FTYPE a, FTYPE b)
{
  if(a>b) return b;
  else return a;
  // if equal, then above is fine
}
#endif

#endif

/// compute speed of light 3-velocity in particular direction assuming other direction velocities fixed
/// dx^dir/dt assuming other direction's velocities are fixed
int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin)
{
  int i,j,k;
  FTYPE vu[NDIM],BB,CC,vsol1p,vsol1m;
  FTYPE ftemp1,ftemp2,ftemp3;
  FTYPE disc;
  int diro1,diro2;

  /* 
     m%3+1 gives next 1->2,2->3,3->1
     3-(4-m)%3 gives previous 1->3,2->1,3->2
  */
  diro1=dir%3+1;
  diro2=3-(4-dir)%3;

  // 3-velocity in coordinate frame
  SLOOPA(j) vu[j]=q->ucon[j]/q->ucon[TT];

  ftemp1=0;
  ftemp2=0;
  ftemp3=0;
  SLOOPA(j){
    if(j!=dir){
      ftemp1+=vu[j]*geom->gcov[GIND(dir,j)];
      ftemp2+=vu[j]*geom->gcov[GIND(TT,j)];
      ftemp3+=vu[j]*vu[j]*geom->gcov[GIND(j,j)];
    }
  }

  BB=2.0*(ftemp1+geom->gcov[GIND(0,dir)])/geom->gcov[GIND(dir,dir)];
  CC=(geom->gcov[GIND(TT,TT)] + 2.0*ftemp2 + ftemp3 + 2.0*vu[diro1]*vu[diro2]*geom->gcov[GIND(diro1,diro2)])/geom->gcov[GIND(dir,dir)];
  
  disc=BB*BB-4.0*CC;
  if(disc>0){
    *vmax=0.5*(-BB+sqrt(disc));
    *vmin=0.5*(-BB-sqrt(disc));
  }
  else{
    dualfprintf(fail_file,"disc=%21.15g < 0\n",disc);
    return(1);
  }

  return(0);

}


/// limit the 3-velocity to a physically valid velocity (i.e. less than c ), assuming all other velocity directions are the same.
int limitv3(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *v)
{
  FTYPE vmax,vmin;
  int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);
  FTYPE ratv;
  
  // get speed of light 3-velocity
  MYFUN(sol(pr,q,dir,geom,&vmax,&vmin),"phys.c:limitv3()", "sol()", 1);

  // get ratio of given 3-velocity to speed of light 3-velocity for appropriate direction of coordinate-based velocity
  ratv=(*v>0) ? *v/vmax : *v/vmin;

  // limit 3-velocity to speed of light
  if(ratv>1.0){
    if(*v>0.0) *v=vmax;
    else *v=vmin;
  }

  return(0);
}

/// take projection of v onto u, both are 4-vectors
/// vcon=0 or 1 : whether or not ucon (1=true, 0=ucov)
/// same for vresultcon
void projectionvec(int vcon,int vresultcon, struct of_state *q, struct of_geom *geom,FTYPE *v,FTYPE*vresult)
{
  FTYPE proj[NDIM][NDIM];
  int j,k;


  if((vcon)&&(vresultcon)){ // vresult^\mu = P^\mu_\nu v^\nu
    DLOOP(j,k) proj[j][k]=delta(j,k) + q->ucon[j]*q->ucov[k];
  }
  if((!vcon)&&(!vresultcon)){ // vresult_\mu = P_\mu^\nu v_\nu
    DLOOP(j,k) proj[j][k]=delta(j,k) + q->ucov[j]*q->ucon[k];
  }
  else if((!vcon)&&(vresultcon)){ // vresult^\mu = P^{\mu\nu} v_\nu
    DLOOP(j,k) proj[j][k]=geom->gcon[GIND(j,k)] + q->ucon[j]*q->ucon[k];
  }
  else if((vcon)&&(!vresultcon)){ // vresult_\mu = P_{\mu\nu} v^\nu
    DLOOP(j,k) proj[j][k]=geom->gcov[GIND(j,k)] + q->ucov[j]*q->ucov[k];
  }
  DLOOPA(j) vresult[j]=0.0;
  DLOOP(j,k) vresult[j]+=proj[j][k]*v[k];
  

}


/// g^{tt}+1 accurate for non-rel gravity to order v^2 without machine precision problems
void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1)
{
  int j,k;
  FTYPE gconttsq;

  // accurate for non-rel gravity since gconttsq is order v^4
  // could store gconttplus1 instead
  gconttsq = (geom->gcon[GIND(TT,TT)]*geom->gcon[GIND(TT,TT)]);
  *gconttplus1= (geom->gcovpert[TT]) * gconttsq + (1.0-gconttsq);
  SLOOPA(j) *gconttplus1 += 2.0*geom->gcov[GIND(TT,j)]*(geom->gcon[GIND(j,TT)]*geom->gcon[GIND(j,TT)]);
  SLOOP(j,k) *gconttplus1 += geom->gcov[GIND(j,k)]*geom->gcon[GIND(j,TT)]*geom->gcon[GIND(k,TT)];

}

/// compute (u^t)^2 - 1 which is v^2 in non-rel regime.
/// Has no non-rel problem and can be used to rescale both non-rel and rel velocities
int quasivsq_3vel(FTYPE *vcon, struct of_geom *geom, FTYPE *quasivsq)
{
  int i,j,k;
  FTYPE gcovttplus1, gconttplus1;
  FTYPE velterm;
  int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm);
  void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1);


  // since (u^t)^2 = -1/ (g_{tt} + 2v^i g_{it} + v^2 )
  // where v^2 = v^i v^j g_{ij}
  //
  // then (u^t)^2 - 1 = [(1 + g_{tt}) + (2v^i g_{it} + v^2)]/[ - g_{tt} - (2v^i g_{it} + v^2) ]

  // later this can be initialized as it's own variable when wanting to do non-rel gravity
  //  gcovttplus1 = geom->gcov[GIND(TT,TT)] + 1.0;
  gcovttplus1 = geom->gcovpert[TT];

  //  gconttplus1= geom->gcon[GIND(TT,TT)] + 1.0;
  compute_gconttplus1(geom,&gconttplus1);


  get_3velterm(vcon, geom, &velterm);

  // this has a good machine representable value in the limit of non-rel velocities
  // *utsqm1 = (gcovttplus1 + velterm) / (-geom->gcov[GIND(TT,TT)] - velterm);

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1)
  *quasivsq =  (gcovttplus1 + velterm) / (-geom->gcov[GIND(TT,TT)] - velterm) + gconttplus1 ;
  
  return(0);


}


int quasivsq_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq)
{
  FTYPE AA,BB,CCM1,CCPLUSAA ;
  FTYPE discr ;
  FTYPE bsq,X[NDIM] ;
  int i=0,j=0,k=0 ;
  FTYPE gcovttplus1,gconttplus1;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);
  FTYPE ucon[NDIM];
  void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  // coefficients for quadratic solution  
  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr<0){
    dualfprintf(fail_file,"Problem with utsqm1_4vel\n");
    return(1);
  }

  // 1.0 + g_{tt}
  gcovttplus1 = geom->gcovpert[TT];

  //  gconttplus1= geom->gcon[GIND(TT,TT)] + 1.0;
  compute_gconttplus1(geom,&gconttplus1);

  // v^2 in non-rel case
  CCPLUSAA = gcovttplus1 + CCM1;

  // first term is genuinely relativistic, so no correction needed for non-rel case
  // second term is combination of |discr| and -4A^2
  // in the non-rel limit this reduces (for BB=0 and AA=-1) to CCPLUSAA, which is v^2
  //  *utsqm1 = ( BB*(BB + sqrt(discr)) + fabs(BB*BB - 4.0*AA*CCPLUSAA) ) /(4.0*AA*AA) ;

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1)
  *quasivsq = ( BB*(BB + sqrt(discr)) + fabs(BB*BB - 4.0*AA*CCPLUSAA) ) /(4.0*AA*AA) + gconttplus1;

  return(0) ;


}


/// interestingly, this can be negative in some reasonable way if we allow qsq<0, which is a similar idea behind the v^2<0 for 3-velocities
/// but I don't see how to write stress-energy tensor as I did for 3-velocity
int quasivsq_rel4vel(FTYPE *uconrel, struct of_geom *geom, FTYPE *quasivsq)
{
  FTYPE alphasq;
  FTYPE qsq;
  FTYPE timeterm;
  int j ;
  int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq);
  
  //        alpha = 1./sqrt(-geom->gcon[GIND(TT,TT)]) ;

  alphasq=1./(-geom->gcon[GIND(TT,TT)]); // positive definite and 1 in non-rel case
  
  //timeterm = -(geom->gcon[GIND(TT,TT)] + 1.0 ); // ~0 in non-rel case, and equal to 1+g_{tt} in non-rel case, which is always positive
  
  qsq_calc(uconrel,geom,&qsq);
  
  // *utsqm1 = qsq/alphasq  + timeterm ;

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1) = utsqm1 - timeterm
  // already accurate to order v^2 when non-rel velocity/gravity used
  *quasivsq = qsq/alphasq;
  
  return(0) ;
}


/// (u^t)^2 - 1 + (g^{tt} + 1)
int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq)
{
  int quasivsq_3vel(FTYPE *vcon, struct of_geom *geom, FTYPE *quasivsq);
  int quasivsq_4vel(FTYPE *ucon, struct of_geom *geom, FTYPE *quasivsq);
  int quasivsq_rel4vel(FTYPE *uconrel, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE ucon[NDIM],vcon[NDIM],uconrel[NDIM];

  
#if(WHICHVEL==VEL4)

  ucon[1]=pr[U1];
  ucon[2]=pr[U2];
  ucon[3]=pr[U3];
  return(quasivsq_4vel(ucon, geom, quasivsq));

#elif(WHICHVEL==VEL3)

  vcon[1]=pr[U1];
  vcon[2]=pr[U2];
  vcon[3]=pr[U3];
  return(quasivsq_3vel(vcon, geom, quasivsq));

#elif(WHICHVEL==VELREL4)

  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];
  return(quasivsq_rel4vel(uconrel, geom, quasivsq));

#endif

  //  return(0);

}

int limit_quasivsq(FTYPE quasivsqnew, struct of_geom *geom, FTYPE *pr)
{
  int limit_quasivsq_3vel(FTYPE quasivsqold,    FTYPE quasivsqnew, struct of_geom *geom, FTYPE *vcon);
  int limit_quasivsq_4vel(FTYPE quasivsqold,    FTYPE quasivsqnew, struct of_geom *geom, FTYPE *ucon);
  int limit_quasivsq_rel4vel(FTYPE quasivsqold, FTYPE quasivsqnew, struct of_geom *geom, FTYPE *uconrel);
  FTYPE ucon[NDIM],vcon[NDIM],uconrel[NDIM];
  int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE quasivsqold;


  // get pr's quasivsq
  if(quasivsq_compute(pr, geom, &quasivsqold)>=1){
    dualfprintf(fail_file,"Problem with limit_quasivsq using quasivsq_compute\n");
    myexit(3);
  }

  //  dualfprintf(fail_file,"IN: v^x = %21.15g v^y = %21.15g v^z = %21.15g :: old = %21.15g new=%21.15g\n",pr[U1],pr[U2],pr[U3],quasivsqold,quasivsqnew);

  // now rescale the velocities to agree with quasivsqnew  
#if(WHICHVEL==VEL4)

  ucon[1]=pr[U1];
  ucon[2]=pr[U2];
  ucon[3]=pr[U3];
  limit_quasivsq_4vel(quasivsqold,quasivsqnew, geom, ucon);
  pr[U1]=ucon[1];
  pr[U2]=ucon[2];
  pr[U3]=ucon[3];

#elif(WHICHVEL==VEL3)

  vcon[1]=pr[U1];
  vcon[2]=pr[U2];
  vcon[3]=pr[U3];
  limit_quasivsq_3vel(quasivsqold,quasivsqnew, geom, vcon);
  pr[U1]=vcon[1];
  pr[U2]=vcon[2];
  pr[U3]=vcon[3];

#elif(WHICHVEL==VELREL4)

  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];
  limit_quasivsq_rel4vel(quasivsqold,quasivsqnew, geom, uconrel);
  pr[U1]=uconrel[1];
  pr[U2]=uconrel[2];
  pr[U3]=uconrel[3];

#endif

  //  dualfprintf(fail_file,"OUT: v^x = %21.15g v^y = %21.15g v^z = %21.15g :: old = %21.15g new=%21.15g\n",pr[U1],pr[U2],pr[U3],quasivsqold,quasivsqnew);

  return(0);

}

int limit_quasivsq_3vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *vcon)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // only true in non-rel case
  vcon[1]*=pref;
  vcon[2]*=pref;
  vcon[3]*=pref;

  return(0);
}

int limit_quasivsq_4vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *ucon)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // only true in non-rel case
  ucon[1]*=pref;
  ucon[2]*=pref;
  ucon[3]*=pref;

  return(0);
}

int limit_quasivsq_rel4vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *uconrel)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // relativistically correct
  uconrel[1]*=pref;
  uconrel[2]*=pref;
  uconrel[3]*=pref;

  return(0);
}


/// input \Omega_F and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];


  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  ftemp=(Bccov[TT]+omegaf*Bccov[PH])/Bsq;
  // ftemp2 is set so that no round off error in limit where toroidal field dominates
  // Does this cause problem for frame-dragged fields?
  //ftemp2=omegaf*(1.0 - Bccov[PH]*Bccon[3]/Bsq);
  // below more accurate than above?
  ftemp2=omegaf*(Bccov[RR]*Bccon[RR]+Bccov[TH]*Bccon[TH])/Bsq;

  vcon[1] = -Bccon[1]*ftemp;
  vcon[2] = -Bccon[2]*ftemp;
  //  vcon[3] =  omegaf  - Bccon[3]*ftemp;
  // designed so to avoid catastrophic cancellation like above has problems with
  vcon[3] =  ftemp2 - (Bccon[3]*Bccov[TT]/Bsq);


  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  //  if(t>1.9 && t<2.1){
  //    dualfprintf(fail_file,"ftemp=%21.15g Bsq=%21.15g ftemp2=%21.15g :: v1=%21.15g v2=%21.15g v3=%21.15g pru1=%21.15g pru2=%21.15g pru3=%21.15g\n",ftemp,Bsq,ftemp2,vcon[1],vcon[2],vcon[3],pr[U1],pr[U2],pr[U3]);
  // }

  return(0);

}



/// input \Omega_F, 3-vel along field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general2(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;

  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  v0oB=v0/sqrt(Bsq);

  ftemp=(Bccov[TT]+omegaf*Bccov[PH])/Bsq;
  // ftemp2 is set so that no round off error in limit where toroidal field dominates
  // Does this cause problem for frame-dragged fields?
  //ftemp2=omegaf*(1.0 - Bccov[PH]*Bccon[3]/Bsq);
  // below more accurate than above?
  ftemp2=omegaf*(Bccov[RR]*Bccon[RR]+Bccov[TH]*Bccon[TH])/Bsq;

  vcon[1] = (v0oB-ftemp)*Bccon[1];
  vcon[2] = (v0oB-ftemp)*Bccon[2];
  //  vcon[3] =  omegaf  - Bccon[3]*ftemp;
  // designed so to avoid catastrophic cancellation like above has problems with
  vcon[3] =  ftemp2 + Bccon[3]*(v0oB-Bccov[TT]/Bsq);


  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  //  if(t>1.9 && t<2.1){
  //    dualfprintf(fail_file,"ftemp=%21.15g Bsq=%21.15g ftemp2=%21.15g :: v1=%21.15g v2=%21.15g v3=%21.15g pru1=%21.15g pru2=%21.15g pru3=%21.15g\n",ftemp,Bsq,ftemp2,vcon[1],vcon[2],vcon[3],pr[U1],pr[U2],pr[U3]);
  // }

  return(0);

}

/// input \Omega_F, extra 3-vel along field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general3(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE absB;
  FTYPE v0other;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;


  //ensure outflow
  if (Bccon[1]<0) {
    v0 *= -1.;
  }


  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  absB=sqrt(Bsq);

  vcon[1] =        v0*Bccon[1]/absB;
  vcon[2] =        v0*Bccon[2]/absB;
  vcon[3] = omegaf+v0*Bccon[3]/absB;

  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  return(0);

}

/// input \Omega_F, extra 3-vel along the poloidal field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general3p(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq_poloidal;
  FTYPE absB_poloidal;
  FTYPE v0other;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;

  lower_vec(Bccon,geom,Bccov);

  Bsq_poloidal=0.0+SMALL;
  SLOOPA(j) if( j <= 2 ) Bsq_poloidal+=Bccon[j]*Bccov[j];  //only poloidal components

  absB_poloidal=sqrt(Bsq_poloidal);

  vcon[1] =        v0*Bccon[1]/absB_poloidal;
  vcon[2] =        v0*Bccon[2]/absB_poloidal;
  vcon[3] = omegaf+v0*Bccon[3]/absB_poloidal;

  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general3p()", "vcon2pr() dir=0", 1);

  return(0);

}




/// input \Omega_F, extra 3-vel along the normal field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general3n(FTYPE omegaf, FTYPE v0, FTYPE *Bccon, FTYPE *normalvec,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq_normal;
  FTYPE absB_normal;
  FTYPE v0other;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;

  FTYPE normalveccov[NDIM];


  lower_vec(Bccon,geom,Bccov);
  lower_vec(normalvec,geom,normalveccov);

  Bsq_normal=0.0;
  SLOOPA(j) Bsq_normal+=Bccon[j]*normalveccov[j];

  absB_normal=sqrt(Bsq_normal);

  vcon[1] =        v0*Bccon[1]/absB_normal;
  vcon[2] =        v0*Bccon[2]/absB_normal;
  vcon[3] = omegaf+v0*Bccon[3]/absB_normal;

  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general3p()", "vcon2pr() dir=0", 1);

  return(0);

}


void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM])
{
  int i,j,k ;
  void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);


  FTYPE Mcov[NDIM][NDIM] ; //  covariant Maxwell 
  FTYPE Msq ;

  //  get covariant Maxwell from contravariant 
  lower_A(Mcon, geom, Mcov) ;

  //  out with the old 
  DLOOP(j,k)  T[j][k] = 0. ;

  /* in with the new:
   *
   * a, b, c, d run 0 to 4:
   *
   * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
   *
   */

  Msq = 0. ;
  DLOOP(j,k) Msq += Mcov[j][k]*Mcon[j][k] ;

  DLOOP(j,k) {
    for(i = 0; i < NDIM; i++) {
      T[j][k] += Mcon[j][i]*Mcov[k][i] ;
    }
  }
  DLOOPA(j) T[j][j] -= 0.25 * Msq ;
}


/// we use this in order to avoid catastrophic cancellation in the non-rel regime and ultrarel regime
/// give back T^dir_\mu
void ffdestresstensor_dir(int dir, FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE *TEMdir)
{
  int i,j,k ;
  void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);


  FTYPE Mcov[NDIM][NDIM] ; //  covariant Maxwell 
  FTYPE Msq ;

  //  get covariant Maxwell from contravariant 
  lower_A(Mcon, geom, Mcov) ;

  //  out with the old 
  DLOOPA(k)  TEMdir[k] = 0. ;

  /* in with the new:
   *
   * a, b, c, d run 0 to 4:
   *
   * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
   *
   */

  Msq = 0. ;
  DLOOP(j,k) Msq += Mcov[j][k]*Mcon[j][k] ;

  DLOOPA(k) {
    for(i = 0; i < NDIM; i++) {
      TEMdir[k] += Mcon[dir][i]*Mcov[k][i] ;
    }
  }
  TEMdir[dir] -= 0.25 * Msq ;
}


void raise_A(FTYPE (*Acov)[NDIM], struct of_geom *geom, FTYPE (*Acon)[NDIM])
{
  int j,k ;
  FTYPE localgcon[SYMMATRIXNDIM];

  //  out with the old 
  DLOOP(j,k){
    Acon[j][k] = 0. ;
    // store since use all terms many times -- compiler doesn't do this for some reason
    localgcon[GIND(j,k)]=geom->gcon[GIND(j,k)];
  }

  //  in with the new 
  DLOOP(j,k) {
    Acon[0][1] += (localgcon[GIND(0,j)])*(localgcon[GIND(1,k)])*Acov[j][k] ;
    Acon[0][2] += (localgcon[GIND(0,j)])*(localgcon[GIND(2,k)])*Acov[j][k] ;
    Acon[0][3] += (localgcon[GIND(0,j)])*(localgcon[GIND(3,k)])*Acov[j][k] ;
    Acon[1][2] += (localgcon[GIND(1,j)])*(localgcon[GIND(2,k)])*Acov[j][k] ;
    Acon[2][3] += (localgcon[GIND(2,j)])*(localgcon[GIND(3,k)])*Acov[j][k] ;
    Acon[3][1] += (localgcon[GIND(3,j)])*(localgcon[GIND(1,k)])*Acov[j][k] ;
  }

  /*
   * diagonal is already set to zero,
   * just copy the rest
   */

  Acon[1][0] = - Acon[0][1] ;
  Acon[2][0] = - Acon[0][2] ;
  Acon[3][0] = - Acon[0][3] ;
  Acon[2][1] = - Acon[1][2] ;
  Acon[3][2] = - Acon[2][3] ;
  Acon[1][3] = - Acon[3][1] ;
}

void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM])
{
  int j,k ;
  FTYPE localgcov[SYMMATRIXNDIM];

  //  out with the old 
  DLOOP(j,k){
    Acov[j][k] = 0. ;
    // store since use all terms many times -- compiler doesn't do this for some reason and lower_A() appears to be very costly
    localgcov[GIND(j,k)]=geom->gcov[GIND(j,k)];
  }

  //  in with the new 
  DLOOP(j,k) {
    Acov[0][1] += (localgcov[GIND(0,j)])*(localgcov[GIND(1,k)])*Acon[j][k] ;
    Acov[0][2] += (localgcov[GIND(0,j)])*(localgcov[GIND(2,k)])*Acon[j][k] ;
    Acov[0][3] += (localgcov[GIND(0,j)])*(localgcov[GIND(3,k)])*Acon[j][k] ;
    Acov[1][2] += (localgcov[GIND(1,j)])*(localgcov[GIND(2,k)])*Acon[j][k] ;
    Acov[2][3] += (localgcov[GIND(2,j)])*(localgcov[GIND(3,k)])*Acon[j][k] ;
    Acov[3][1] += (localgcov[GIND(3,j)])*(localgcov[GIND(1,k)])*Acon[j][k] ;
  }

  /*
   * diagonal is already set to zero,
   * just copy the rest
   */

  Acov[1][0] = - Acov[0][1] ;
  Acov[2][0] = - Acov[0][2] ;
  Acov[3][0] = - Acov[0][3] ;
  Acov[2][1] = - Acov[1][2] ;
  Acov[3][2] = - Acov[2][3] ;
  Acov[1][3] = - Acov[3][1] ;
}




/// Maxwell to Faraday
/// which=0 : Mcon -> Fcov (for clean Mcon, Fcov has \detg)
/// which=1 : Mcov -> Fcon (for clean Mcov)
/// which=2 : Fcon -> Mcov
/// which=3 : Fcov -> Mcon
/// copies faraday_calc() in phys.c
void MtoF(int which, FTYPE (*invar)[NDIM],struct of_geom *geom, FTYPE (*outvar)[NDIM])
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
          outvar[mu][nu]+=prefactor*lc4(whichlc,geom->gdet,mu,nu,kappa,lambda)*invar[kappa][lambda];
        }
    }


}


/// Get dual of faraday for given B^i and v^i
int dualf_calc(FTYPE *Bcon, FTYPE *vcon, FTYPE (*dualffull)[NDIM])
{
  int j,k;

  SLOOPA(j){
    // \dF^{it} = B^i
    dualffull[j][0] = Bcon[j];
    dualffull[0][j] = -Bcon[j];
  }

  DLOOPA(j) dualffull[j][j] = 0;

  // \dF^{ij} = B^i v^j - B^j v^i
  SLOOP(j,k) dualffull[j][k] = Bcon[j] * vcon[k] - Bcon[k] * vcon[j] ;

  return(0);

}


/// sets velocity U1-U3 part of primitive to be zamo velocity
int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  int jj;

  // bl-normal observer (4-vel components)

  FTYPE etacov[NDIM],etacon[NDIM];
  etacov[TT]=-ptrgeom->alphalapse;
  etacov[RR]=0;
  etacov[TH]=0;
  etacov[PH]=0;

  raise_vec(etacov,ptrgeom,etacon);
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    SLOOPA(jj) pr[U1+jj-1] = etacon[jj];
  }
  else if(whichvel==VEL3){
    SLOOPA(jj) pr[U1+jj-1] = etacon[jj]/etacon[TT];
  }
  else if(whichvel==VELREL4){
    SLOOPA(jj) pr[U1+jj-1] = 0.0;
  }

  return(0);

}


int set_zamo_ucovuconplus1ud0(struct of_geom *ptrgeom, FTYPE *ucov, FTYPE *ucon, FTYPE *plus1ud0)
{
  FTYPE ialpha;
  FTYPE alpha;
  int j;
  FTYPE gconpert;

  // set alpha -- fabs just for roundoff error
  //  ialpha=sqrt(fabs(-ptrgeom->gcon[GIND(TT,TT)]));
  alpha = ptrgeom->alphalapse;
  ialpha = 1.0/alpha;

  ucov[TT]=-alpha;
  SLOOPA(j) ucov[j]=0.0;

  raise_vec(ucov,ptrgeom,ucon);

  // 1+u_t
  //  gconpert = ptrgeom->gconpert[TT];  // should be GODMARK
  gconpert = 0.0; // for now since this function now used GODMARK

  *plus1ud0 = fabs(gconpert)*alpha/(ialpha+1.0);

  return(0);

}

int set_zamo_ucon(struct of_geom *ptrgeom, FTYPE *ucon)
{
  FTYPE ucov[NDIM];
  FTYPE plus1ud0;

  set_zamo_ucovuconplus1ud0(ptrgeom, ucov,ucon,&plus1ud0);

  return(0);

}






/// entropy wrapper
/// this function should NOT be called by utoprim_jon.c inversion
int entropy_calc(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy)
{
  
  *entropy = compute_entropy_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pr[RHO],pr[UU]);

  return(0);
}

/// entropy wrapper
/// this function should NOT be called by utoprim_jon.c inversion
int entropy_calc_forcheckinversion(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy)
{
  
  *entropy = compute_entropy_simple_forcheckinversion(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pr[RHO],pr[UU]);

  return(0);
}



/// wrapper [assumed not called by utoprim_jon.c that could change EOS type]
FTYPE pressure_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(pressure_rho0_u(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u));

}

/// wrapper [assumed not called by utoprim_jon.c that could change EOS type]
/// used for inversion check
FTYPE pressure_rho0_u_simple_forcheckinversion(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  if(WHICHEOS==KAZFULL && loc==CENT){
    // then avoid inconsistency by using stored pressure that was perhaps p(rho0,chi) instead of p(rho0,u)
    return(GLOBALMACP0A1(EOSextraglobal,i,j,k,PGASGLOBAL));
  }
  else{
    // then didn't store anything
    return(pressure_rho0_u(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u));
  }

}

/// wrapper [assumed not called by utoprim_jon.c that could change EOS type]
FTYPE u_rho0_p_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE p)
{

  return(u_rho0_p(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,p));

}

/// wrapper [assumed not called by utoprim_jon.c that could change EOS type]
FTYPE u_rho0_T_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE T)
{

  return(u_rho0_T(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,T));

}

/// wrapper
FTYPE cs2_compute_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(cs2_compute(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u));

}

/// wrapper
/// this function should NOT be called by utoprim_jon.c inversion
FTYPE compute_entropy_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{
  FTYPE entropy;


  // unlike inversion, require non-NaN entropy, so force rho,u to be positive
  if(rho<SMALL) rho=SMALL;

#if(WHICHEOS!=KAZFULL)
  // for Kaz EOS, u<0 is ok
  if(u<SMALL) u=SMALL;
#endif
  
  entropy=compute_entropy(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u);

  if(!isfinite(entropy)) entropy=-BIG*MAX(rho,SMALL); // very negative corresponding to very little entropy.

  return(entropy);

}

/// wrapper
/// this function should NOT be called by utoprim_jon.c inversion
FTYPE compute_entropy_simple_forcheckinversion(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  if(WHICHEOS==KAZFULL && loc==CENT){
    // then avoid inconsistency by using stored pressure
    FTYPE chi = u + GLOBALMACP0A1(EOSextraglobal,i,j,k,PGASGLOBAL);
    return(rho*compute_specificentropy_wmrho0(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,chi));
  }
  else{
    // unlike inversion, require non-NaN entropy, so force rho,u to be positive
    if(rho<SMALL) rho=SMALL;
#if(WHICHEOS!=KAZFULL)
    // for Kaz EOS, u<0 is ok
    if(u<SMALL) u=SMALL;
#endif
    
    return(compute_entropy(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u));
  }

}

/// wrapper
void get_EOS_parms_simple(int*numparms, int i, int j, int k, int loc, FTYPE *parlist)
{

  get_EOS_parms(WHICHEOS,numparms, GLOBALMAC(EOSextraglobal,i,j,k),parlist);


}

/// wrapper
void fix_primitive_eos_scalars_simple(int i, int j, int k, int loc, FTYPE *pr)
{

  fix_primitive_eos_scalars(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),pr);


}

FTYPE compute_temp_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(compute_temp(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k),rho,u));

}


int get_extrasprocessed_simple(int doall, int i, int j, int k, int loc, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{

  return(get_extrasprocessed(WHICHEOS,1, GLOBALMAC(EOSextraglobal,i,j,k), GLOBALMAC(pdump,i,j,k), extras, processed));

}


/// this function should NOT be called by utoprim_jon.c inversion
FTYPE compute_u_from_entropy_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE entropy)
{

  // unlike inversion, require non-NaN entropy, so force rho,u to be positive
  if(rho<SMALL) rho=SMALL;
  // entropy can be any value

  return(compute_u_from_entropy(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k), rho, entropy));

}

FTYPE compute_qdot_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(compute_qdot(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k), rho, u));

}

FTYPE dpdrho0_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(dpdrho0_rho0_u(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k), rho, u));

}

FTYPE dpdu_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u)
{

  return(dpdu_rho0_u(WHICHEOS,GLOBALMAC(EOSextraglobal,i,j,k), rho, u));

}
