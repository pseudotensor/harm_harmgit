
/*! \file wavespeeds.c
     \brief Wrapper for using vchar and storing wavespeed results
*/


#include "decs.h"



///////////////////////////////////
///
/// get wave speeds for flux calculation and for dt calculation
///
//////////////////////////////////
int vchar_all(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxall, FTYPE *vminall,int *ignorecourant)
{
  FTYPE vminmhd,vmaxmhd;
  FTYPE vminrad,vmaxrad;
  FTYPE vminrad2,vmaxrad2;
  
  vchar_each(pr, q, dir, geom, &vmaxmhd, &vminmhd, &vmaxrad, &vminrad, &vmaxrad2, &vminrad2, ignorecourant);
  // below correct even if EOMRADTYPE==EOMRADNONE because vchar_each() sets both mhd and rad to be mhd and so below always chooses the mhd values.
  *vminall=MIN(vminmhd,vminrad2); // vminrad2 for dt
  *vmaxall=MAX(vmaxmhd,vmaxrad2); // vmaxrad2 for dt
  

  return(0);
}

int vchar_each(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxmhd, FTYPE *vminmhd, FTYPE *vmaxrad, FTYPE *vminrad, FTYPE *vmaxrad2, FTYPE *vminrad2,int *ignorecourant)
{
  
  vchar(pr, q, dir, geom, vmaxmhd, vminmhd,ignorecourant);
  if(EOMRADTYPE!=EOMRADNONE){
    vchar_rad(pr, q, dir, geom, vmaxrad, vminrad, vmaxrad2, vminrad2, ignorecourant);
  }
  else{// default as if no other values for wave speeds
    *vmaxrad2=*vmaxrad=*vmaxmhd;
    *vminrad2=*vminrad=*vminmhd;
  }

  //  dualfprintf(fail_file,"SPEEDS: %g %g  : %g %g\n",*vmaxmhd,*vminmhd,*vmaxrad2,*vminrad2);

  return(0);
}


/// get wave speeds for flux calculation
int get_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, struct of_state *state_l, struct of_state * state_r, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE *ctopptr, FTYPE *cminmaxrad_l, FTYPE *cminmaxrad_r, FTYPE *cminmaxrad, FTYPE *ctopradptr, FTYPE *cminmaxrad2_l, FTYPE *cminmaxrad2_r, FTYPE *cminmaxrad2, FTYPE *ctoprad2ptr)
{
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  int cminmax_calc_1(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax);
  int cminmax_calc_2(FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  int failreturn;
  FTYPE ftemp;
  FTYPE p_roe[NPR],cminmax_roe[NUMCS];
  struct of_state state_roe;
  FTYPE cminmax_o[NUMCS];
  FTYPE Hroe,cminmaxnonrel_roe[NUMCS];
  void get_roe_averaged_state(int dir, FTYPE *p_l, struct of_state *state_l, FTYPE *Ul, FTYPE * p_r, struct of_state *state_r, FTYPE *Ur, struct of_geom *geom, FTYPE * p_roe, FTYPE *Hroe, FTYPE *cminnonrel_roe, FTYPE*cmaxnonrel_roe);
  int q;
  int ignorecourant;
  int i,j,k,loc;

  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;



#if(USEGLOBALWAVE==0 || STOREWAVESPEEDS==2)

  // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
  MYFUN(vchar_each(p_l, state_l, dir, ptrgeom, &cminmax_l[CMAX], &cminmax_l[CMIN], &cminmaxrad_l[CMAX], &cminmaxrad_l[CMIN], &cminmaxrad2_l[CMAX], &cminmaxrad2_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
  MYFUN(vchar_each(p_r, state_r, dir, ptrgeom, &cminmax_r[CMAX], &cminmax_r[CMIN], &cminmaxrad_r[CMAX], &cminmaxrad_r[CMIN], &cminmaxrad2_r[CMAX], &cminmaxrad2_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);

  cminmax_calc(cminmax_l[CMIN],cminmax_r[CMIN],cminmax_l[CMAX],cminmax_r[CMAX],&cminmax[CMIN],&cminmax[CMAX],ctopptr);
  //  cminmax_calc_1(cminmax_l[CMIN],cminmax_r[CMIN],cminmax_l[CMAX],cminmax_r[CMAX],&cminmax[CMIN],&cminmax[CMAX]);
  if(EOMRADTYPE!=EOMRADNONE){
    cminmax_calc(cminmaxrad_l[CMIN],cminmaxrad_r[CMIN],cminmaxrad_l[CMAX],cminmaxrad_r[CMAX],&cminmaxrad[CMIN],&cminmaxrad[CMAX],ctopradptr);

    cminmax_calc(cminmaxrad2_l[CMIN],cminmaxrad2_r[CMIN],cminmaxrad2_l[CMAX],cminmaxrad2_r[CMAX],&cminmaxrad2[CMIN],&cminmaxrad2[CMAX],ctoprad2ptr);
    //  cminmax_calc_1(cminmaxrad_l[CMIN],cminmaxrad_r[CMIN],cminmaxrad_l[CMAX],cminmaxrad_r[CMAX],&cminmaxrad[CMIN],&cminmaxrad[CMAX]);
  }
  // have cmin,cmax,ctop here

#if(STOREWAVESPEEDS==2)
  // Use fact that always compute flux before need wspeed [GODMARK: Not sure for old interpline use of wspeed, but that's deprecated code]
  // remove sign information for wspeed by calling cminmax_calc() above already
  GLOBALMACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k)=cminmax[CMIN];
  GLOBALMACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k)=cminmax[CMAX];
  if(EOMRADTYPE!=EOMRADNONE){
    GLOBALMACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k)=cminmaxrad[CMIN];
    GLOBALMACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k)=cminmaxrad[CMAX];

    GLOBALMACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k)=cminmaxrad2[CMIN];
    GLOBALMACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k)=cminmaxrad2[CMAX];
  }
#endif

  //  cminmax_calc_2(&cminmax[CMIN],&cminmax[CMAX],ctopptr);

#else
  // other non-local estimate
  cminmax_l[CMIN]=cminmax_r[CMIN]=cminmax[CMIN]=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k);
  cminmax_l[CMAX]=cminmax_r[CMAX]=cminmax[CMAX]=GLOBALMACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k);
  if(EOMRADTYPE!=EOMRADNONE){
    cminmaxrad_l[CMIN]=cminmaxrad_r[CMIN]=cminmaxrad[CMIN]=GLOBALMACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k);
    cminmaxrad_l[CMAX]=cminmaxrad_r[CMAX]=cminmaxrad[CMAX]=GLOBALMACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k);

    cminmaxrad2_l[CMIN]=cminmaxrad2_r[CMIN]=cminmaxrad2[CMIN]=GLOBALMACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k);
    cminmaxrad2_l[CMAX]=cminmaxrad2_r[CMAX]=cminmaxrad2[CMAX]=GLOBALMACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k);
  }
  // change so feeds into Riemann solver with "standard" sign
  cminmax[CMIN] = fabs(max(0., -cminmax[CMIN]));
  cminmax[CMAX] = fabs(max(0., cminmax[CMAX]));
  *ctopptr=max(cminmax[CMIN],cminmax[CMAX]);

  if(EOMRADTYPE!=EOMRADNONE){
    cminmaxrad[CMIN] = fabs(max(0., -cminmaxrad[CMIN]));
    cminmaxrad[CMAX] = fabs(max(0., cminmaxrad[CMAX]));
    *ctopradptr=max(cminmaxrad[CMIN],cminmaxrad[CMAX]);

    cminmaxrad2[CMIN] = fabs(max(0., -cminmaxrad2[CMIN]));
    cminmaxrad2[CMAX] = fabs(max(0., cminmaxrad2[CMAX]));
    *ctoprad2ptr=max(cminmaxrad2[CMIN],cminmaxrad2[CMAX]);
  }
  // have cmin,cmax,ctop here
#endif



  //////////////////////////////////
  //
  // get Roe-averaged versions of wave speed for flux calculation and for dt calculation
  //
  /////////////////////////////////

#if(ROEAVERAGEDWAVESPEED)

  if(EOMRADTYPE!=EOMRADNONE){
    dualfprintf(fail_file,"ROEAVERAGE not setup for radiation\n");
    myexit(1);
  }

  // get Roe-averaged primitive state
  get_roe_averaged_state(dir, p_l, state_l, U_l, p_r, state_r, U_r, ptrgeom, p_roe, &Hroe, &cminmaxnonrel_roe[CMIN], &cminmaxnonrel_roe[CMAX]);


#if(ATHENAROE==0) // doing Jon method

  // get HARM state
  MYFUN(get_state(p_roe, ptrgeom, &state_roe),"flux.c:p2SFUevolve", "get_state()", 1);
  // get Roe-averaged state's wavespeeds
  MYFUN(vchar(p_roe, &state_roe, dir, ptrgeom, &cminmax_roe[CMAX], &cminmax_roe[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 3);


#if(1)  // Jon's version of MIN/MAXing

  for(q=0;q<NUMCS;q++){
    // fastest left-going wave (q=0) and then fastest right-going wave (q=1)
    cminmax[q]=MINMAX(q,cminmax_roe[q],MINMAX(q,cminmax_l[q],cminmax_r[q]));
  }


#else // doing Athena-like method of min/maxing (fails to work for some tests)

  // New version of cminmax
  /*
    for(q=0;q<NUMCS;q++){
    cminmax_l[q] = MINMAX(q,cminmax_l[q],0.0);
    cminmax_r[q] = MINMAX(q,cminmax_r[q],0.0);
    cminmax[q]=2.0*cminmax_l[q]*cminmax_r[q]/(cminmax_l[q]+cminmax_r[q]+SMALL);
    }
  */

  cminmax[CMIN] = cminmax_l[CMIN];
  cminmax[CMAX] = cminmax_r[CMAX];
  for(q=0;q<NUMCS;q++){
    cminmax[q]=MINMAX(q,cminmax_roe[q],cminmax[q]);
  }
#endif




#else // else if ATHENAROE==1  (fails to work for some tests)

  // non-rel HD version
  for(q=0;q<NUMCS;q++){
    cminmax_roe[q]=cminmaxnonrel_roe[q];
  }


#if(1) // Athena way
  // This is Athena's way of min/maxing
  cminmax[CMAX]=MAX(cminmax_roe[CMAX],cminmax_r[CMAX]);
  cminmax[CMIN]=MIN(cminmax_roe[CMIN],cminmax_l[CMIN]);

#else // testing way (same as Athena right now)

  // New version of cmin/cmax
  /*
    for(q=0;q<NUMCS;q++){
    cminmax_l[q] = MINMAX(q,cminmax_l[q],0.0);
    cminmax_r[q] = MINMAX(q,cminmax_r[q],0.0);
    cminmax[q] = 2.0*cminmax_l[q]*cminmax_r[q]/(cminmax_l[q] + cminmax_r[q] + SMALL);
    }
  */

  // Jon's way even with ATHENAROE==1
  for(q=0;q<NUMCS;q++){
    // fastest left-going wave (q=0) and then fastest right-going wave (q=1)
    cminmax[q]=MINMAX(q,cminmax_roe[q],MINMAX(q,cminmax_l[q],cminmax_r[q]));
  }
#endif



#endif // end over ATHENAROE=0/1



  // fastest eigenvalue
  // Athena sets this, but never uses it apparently
  //ctop_roe=MAX(fabs(cminmax_roe[CMAX]),fabs(cminmax_roe[CMIN]));


  // change so feeds into Riemann solver with "standard" sign
  cminmax[CMIN] = fabs(max(0., -cminmax[CMIN]));
  cminmax[CMAX] = fabs(max(0., cminmax[CMAX]));
  *ctopptr=max(cminmax[CMIN],cminmax[CMAX]);

  // have Roe versions of cmin,cmax,ctop here
#endif // end over ROEAVERAGEDWAVESPEED=1

  return(0);

}




/// get Roe-averaged primitive state
/// based upon Athena2's flux_hlle.c
void get_roe_averaged_state(int dir, FTYPE *p_l, struct of_state *state_l, FTYPE *Ul, FTYPE * p_r, struct of_state *state_r, FTYPE *Ur, struct of_geom *geom, FTYPE * p_roe, FTYPE *Hroe, FTYPE *cminnonrel_roe, FTYPE*cmaxnonrel_roe)
{
  FTYPE sqrtrhol,sqrtrhor,isqrtrholr;
  int j,k;
  FTYPE Pl, Pr, Hl, Hr; // specific enthalpy
  FTYPE bsql,bsqr;
  FTYPE vsqroe;
  FTYPE croe;
  int ii=geom->i;
  int jj=geom->j;
  int kk=geom->k;
  int loc=geom->p;

  Pl=pressure_rho0_u_simple(ii,jj,kk,loc,p_l[RHO],p_l[UU]);
  Pr=pressure_rho0_u_simple(ii,jj,kk,loc,p_r[RHO],p_r[UU]);


  // get weighted density
  sqrtrhol=sqrt(p_l[RHO]);
  sqrtrhor=sqrt(p_r[RHO]);
  isqrtrholr=1.0/(sqrtrhol+sqrtrhor);

  // Roe-averaged density
  p_roe[RHO] = sqrtrhol*sqrtrhor;

  // Roe-averaged internal energy (no consideration of geometry)
  p_roe[UU] = p_roe[RHO] * ( p_l[UU] / sqrtrhol + p_r[UU] / sqrtrhor ) * isqrtrholr;

  // Roe-averaged velocity (no consideration of geometry)
  SLOOPA(j) p_roe[UU+j] = (sqrtrhol*p_l[UU+j] + sqrtrhor*p_r[UU+j])*isqrtrholr;

  // Roe-averaged magnetic field (no consideration of geometry)
  SLOOPA(j) p_roe[B1-1+j] = (sqrtrhor*p_l[B1-1+j] + sqrtrhol*p_r[B1-1+j])*isqrtrholr;


  // b^2
  bsql = dot(state_l->bcon, state_l->bcov);
  bsqr = dot(state_r->bcon, state_r->bcov);


#if(!ATHENAROE)
  // Jon version (relativistically correct)
  Hl = (p_l[UU] + Pl + 0.5*bsql)/p_l[RHO];
  Hr = (p_r[UU] + Pr + 0.5*bsqr)/p_r[RHO];
  *Hroe=(sqrtrhol*Hl + sqrtrhor*Hr)*isqrtrholr;
#else
  // Athena version (relativistically correct)
  *Hroe = ((Ul[UU] + Pl + 0.5*bsql)/sqrtrhol + (Ur[UU] + Pr + 0.5*bsqr)/sqrtrhor)*isqrtrholr;
#endif

  /////////////////////////////////
  // GET NON-REL CASE wave speeds
  // try Einfeldt (1988) 5.1a - 5.3 for non-rel HD
  // non-rel HD Roe-averaged version of wavespeeds

  // v^2 in non-rel case
  vsqroe=0.0;
  SLOOP(j,k) vsqroe += p_roe[UU+j]*p_roe[UU+k] * geom->gcov[GIND(j,k)];

  // ideal gas only GODMARK (gam)
  croe=sqrt( (gam-1.0)*(*Hroe - 0.5*vsqroe));
  *cminnonrel_roe = p_roe[UU+dir] - croe;
  *cmaxnonrel_roe = p_roe[UU+dir] + croe;

  
  


}



/// complete storage of wave speeds per dimension
/// only called for STOREWAVESPEEDS==1
int get_global_wavespeeds_full(int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*finalwspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3])
{


  // OPTMARK: bit excessive
  //  COMPFULLLOOP{
#if(STOREFLUXSTATE)
  // if stored state, then don't need EOS stuff
#pragma omp parallel 
#else
#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM
#endif
  {

    int i,j,k;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;

    // generally ptr's are different inside parallel block
    ptrgeom=&geomdontuse;

#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      // get geometry for center pre-interpolated values
      get_geometry(i, j, k, CENT, ptrgeom); 

      MYFUN(get_global_wavespeeds(dir,ptrgeom, MAC(prim,i,j,k),GLOBALMACP1A0(wspeedtemp,EOMSETMHD,i,j,k),GLOBALMACP1A0(wspeedtemp,EOMSETRAD,i,j,k),GLOBALMACP1A0(wspeedtemp,EOMSETRADFORDT,i,j,k)),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
    }
  }// end parallel region



  // get char. velocity estimate as some average over local or global zones
  global_vchar(GLOBALPOINT(wspeedtemp), dir, is, ie, js, je, ks, ke, idel, jdel, kdel, POINT(finalwspeed));
  // now wspeed contains left/right fastest wave speed (with sign preserved)

  return(0);

}




/// store wavespeeds somewhere
int get_global_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *pr,FTYPE *output,FTYPE *outputrad,FTYPE *outputrad2)
{
  struct of_state qdontuse;
  struct of_state *qptr=&qdontuse;
  int ignorecourant;

  // wave speed for cell centered value
  // uses b^\mu b_\mu so need full state
  // OPTMARK: Should avoid use of b^\mu and b_\mu in vchar and if STOREFLUXSTATE, then use that data instead of recomputing.
  MYFUN(get_stateforglobalwavespeeds(pr, ptrgeom, &qptr),"step_ch.c:fluxcalc()", "get_state()", 0);
  MYFUN(vchar_each(pr, qptr, dir, ptrgeom, &output[CMAX], &output[CMIN], &outputrad[CMAX], &outputrad[CMIN], &outputrad2[CMAX], &outputrad2[CMIN],&ignorecourant),"wavespeeds.c:get_global_wavespeeds()", "vchar() dir=1or2", 0);
  
  // uses output as temporary space since not yet needed and not used before global_vchar() below
  
  return(0);
}



/// GODMARK: something wrong with comparing multiple velocities since grid/metric changes in space (i.e. v=dx/dt means something different at each grid point)
/// GODMARK: assumes boundary zones exist (flux method of bounding won't work) -- have to apply extra limits on values (i,j,k) used here
///
/// defines an effective maximum wave speed centered on the cell interface (FACE)
///
/// might choose wavespeeds that correspond to interpolation stencil, which to first approximation is a symmetric stencil of size interporder[reallim]
int global_vchar(FTYPE (*pointspeed)[NSTORE1][NSTORE2][NSTORE3][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*wspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3])
{


#pragma omp parallel
  {
    int i,j,k;
    int reallim;
    int m;
    FTYPE ftemp[NUMCS];
    FTYPE ctop,ctoprad,ctoprad2;
    extern int choose_limiter(int dir, int i, int j, int k, int pl);
    int cminmax_calc_2(FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ); // due to averaging of this face quantity to get centered i/j/k=0 back


    // initialize ftemp if needed
    ftemp[CMIN]=+BIG;
    ftemp[CMAX]=-BIG;


    // COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) { // due to averaging of this face quantity to get centered i/j/k=0 back
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


#if(VCHARTYPE==VERYLOCALVCHAR) // then get maximum around interface

      // get minimum/maximum wspeed over this stencil (but centered on edge, not zone center)
      MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k)=+BIG; // assume all zones are smaller than this
      MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k)=-BIG; // assume all zones are larger than this
      if(EOMRADTYPE!=EOMRADNONE){
        MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k)=+BIG; // assume all zones are smaller than this
        MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k)=-BIG; // assume all zones are larger than this
        MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k)=+BIG; // assume all zones are smaller than this
        MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k)=-BIG; // assume all zones are larger than this
      }
      // GODMARK: 1D, could do multi-D stencil even if interpolation is 1D
      for(m=-1;m<=0;m++){
        MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETMHD,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
        MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETMHD,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
      }
      if(EOMRADTYPE!=EOMRADNONE){
        for(m=-1;m<=0;m++){
          MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETRAD,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
          MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETRAD,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
          MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETRADFORDT,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
          MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETRADFORDT,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
        }
      }

#elif(VCHARTYPE==LOCALVCHAR)

      reallim=choose_limiter(dir,i,j,k,RHO); // base stencil size on density limiter GODMARK

      // get minimum/maximum wspeed over this stencil (but centered on edge, not zone center)
      MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k)=+BIG; // assume all zones are smaller than this
      MACP2A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k)=-BIG; // assume all zones are larger than this
      // GODMARK: 1D, could do multi-D stencil even if interpolation is 1D
      // e.g. if dir=1, then expandi=0 and so COMPFZLOOP from i=0..N inclusive.  So can grab from (relative to i) -2 .. 1 inclusive for average centered on i edge
      for(m=-reallim/2;m<=reallim/2-1;m++){
        MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETMHD,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
        MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETMHD,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
      }
      if(EOMRADTYPE!=EOMRADNONE){
        for(m=-reallim/2;m<=reallim/2-1;m++){
          MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETRAD,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
          MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETRAD,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
          MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k)=MIN(MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k),MACP1A1(pointspeed,EOMSETRADFORDT,i+idel*m,j+jdel*m,k+kdel*m,CMIN));
          MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k)=MAX(MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k),MACP1A1(pointspeed,EOMSETRADFORDT,i+idel*m,j+jdel*m,k+kdel*m,CMAX));
        }
      }
      
#elif(VCHARTYPE==GLOBALVCHAR)

      dualfprintf(fail_file,"VCHARTYPE==GLOBALVCHAR deprecated\n");
      myexit(346983);

      // COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) { // due to averaging of this face quantity to get centered i/j/k=0 back
      //    // get minimum/maximum wspeed over entire grid containing the fluxes
      //    ftemp[CMIN]=MIN(ftemp[CMIN],MACP1A1(pointspeed,i,j,k,CMIN));
      //    ftemp[CMAX]=MAX(ftemp[CMAX],MACP1A1(pointspeed,i,j,k,CMAX));
      //  }
      // COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) {
      //    MACP2A0(wspeed,dir,CMIN,i,j,k)=ftemp[CMIN];
      //    MACP2A0(wspeed,dir,CMAX,i,j,k)=ftemp[CMAX];
      //  }
#endif

      // go ahead and also restrict cmin and cmax
      cminmax_calc_2(&MACP3A0(wspeed,EOMSETMHD,dir,CMIN,i,j,k),&MACP3A0(wspeed,EOMSETMHD,dir,CMAX,i,j,k),&ctop); // ctop not used
      if(EOMRADTYPE!=EOMRADNONE){
        cminmax_calc_2(&MACP3A0(wspeed,EOMSETRAD,dir,CMIN,i,j,k),&MACP3A0(wspeed,EOMSETRAD,dir,CMAX,i,j,k),&ctoprad); // ctoprad not used
        cminmax_calc_2(&MACP3A0(wspeed,EOMSETRADFORDT,dir,CMIN,i,j,k),&MACP3A0(wspeed,EOMSETRADFORDT,dir,CMAX,i,j,k),&ctoprad2); // ctoprad2 not used
      }

    }// end 3D loop
  }// end parallel region


  return(0);

}




/// really HARM is currently using VERY local lax Friedrich.
/// maybe try local lax Friedrich, using max wave speed from zones used to reconstruct the zone (most common?)
/// also can try more global wave speed, or even speed of light.
int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop)
{
  int cminmax_calc_1(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax);
  int cminmax_calc_2(FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);


  cminmax_calc_1(cmin_l,cmin_r,cmax_l,cmax_r,cmin,cmax);
  cminmax_calc_2(cmin,cmax,ctop);

  return(0);

}

/// determine cmin,cmax,ctop
int cminmax_calc_1(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax)
{
  FTYPE lmin,lmax,ltop;

  // need to make use of ROE averages
      
      
  // As per the HLLE method, one must use the wave speed at the interface.  The below is arbitrarily choosing the largest of the interpolated states.

  // Apparently one should use Roe's linearisation to compare
  // against, not the left/right interpolations.  Then compare the
  // left most going and right most going Rho linearisation
  // eigenvalues with the left and right numerical approximations
  // and choose the minimum and maximum values as cmin and cmax,
  // all respectively.
  // Then one defines the new cmin as the minumum of cmin and 0, and new cmax as maximum of cmax and 0.
  // Thus this is slightly different and doesn't avoid negative mass/ie densities like HLLE.


  // LAXF here is really the Rusanov flux.  Lax-Friedrichs would say ctop->dx/dt


#if(1)
  lmin=min(cmin_l,cmin_r);
  lmax=max(cmax_l,cmax_r);
#elif(0)
  lmin=cmin_l;
  lmax=cmax_r;
#else
  lmin=2.0*cmin_l*cmin_r/(cmin_l+cmin_r);
  lmax=2.0*cmax_l*cmax_r/(cmax_l+cmax_r);
#endif

  // returns signed speed still
  *cmin=lmin;
  *cmax=lmax;

  return(0);
      


}




/// determine cmin,cmax,ctop
int cminmax_calc_2(FTYPE *cmin,FTYPE *cmax,FTYPE *ctop)
{
  FTYPE lmin=*cmin;
  FTYPE lmax=*cmax;

#if(HYPERHLL==0)
  // sharp change at c=0
  lmin = fabs(max(0., -lmin));
  lmax = fabs(max(0., lmax));
#else

#define HYPERFUN(c,x0) (0.5*( (c)+(M_SQRT2l-1.0)*sqrt(2.0*(M_SQRT2l+1.0)*(x0)*(x0)+(3.0+2.0*M_SQRT2l)*(c)*(c)) ) )
  // GODMARK
  // need a good way to set x0 based upon fraction of local light speed
  // or something...

  // function is always positive for given lmin/lmax as original sign
  lmin=HYPERFUN(-lmin,1E-4);
  lmax=HYPERFUN(lmax,1E-4);
#endif

  // returns positive cmin,cmax,ctop
  *cmin=lmin;
  *cmax=lmax;
  *ctop=max(lmax,lmin);


  return(0);
      


}

