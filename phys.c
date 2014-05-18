
/*! \file phys.c
     \brief All globally-related physics calculations
     // THINGS IN HERE ARE NOT PER POINT OR USE GLOBAL VARIABLES or call get_geometry()
*/


#include "decs.h"



/// making variables static inside functions made things slower if anything
/// would have thought consruction/destruction operations would be important
/// why static makes slower?
//#define VARSTATIC static
#define VARSTATIC 


/// macros to avoid messy function calls or restoring them as i,j,k
#define MYII ptrgeom->i
#define MYJJ ptrgeom->j
#define MYKK ptrgeom->k


/// wrapper for choosing to get state by computing it or from global array
/// notice that ptrgeom->p inputted not used
int get_stateforfluxcalc(int dimen, int isleftright, FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr)
{
  int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  struct of_state *qptrlocal;


#if(STOREFLUXSTATE)
  // retrieve state from global array


  // JCM: I see no point in copying the data, just should change pointer address
  // pointer assign (overrides the existing default assignment of the pointer address in the calling function, so assumes calling function only uses pointer form subsequently)
  if(isleftright==ISMIDDLE){
    *qptr=&(GLOBALMAC(fluxstatecent,MYII,MYJJ,MYKK));
  }
  else{
    *qptr=&(GLOBALMACP1A1(fluxstate,dimen,MYII,MYJJ,MYKK,isleftright));
  }

#else // else if not storing fluxstate
  // compute state right now
  pureget_stateforfluxcalc(pr, ptrgeom, *qptr);
#endif

  return (0);
}


/// wrapper for choosing to get state by computing it or from global array
/// pointer assign (overrides the existing default assignment of the pointer address in the calling function, so assumes calling function only uses pointer form subsequently)
int get_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr)
{
  int pureget_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  //  int i,j,k;



#if(STOREFLUXSTATE)
  // retrieve state from global array
  //  i=ptrgeom->i;
  //  j=ptrgeom->j;
  //  k=ptrgeom->k;

  *qptr=&(GLOBALMAC(fluxstatecent,MYII,MYJJ,MYKK)); // note that fluxstatecent for source has more than needed for source since was computed for flux+source and flux>source for all terms
#else
  // compute state right now
  pureget_stateforsource(pr, ptrgeom, *qptr);
#endif


  return (0);
}

/// wrapper for choosing to get state by computing it or from global array
/// pointer assign (overrides the existing default assignment of the pointer address in the calling function, so assumes calling function only uses pointer form subsequently)
int get_stateforinterpline(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr)
{
  int pureget_stateforinterpline(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  //  int i,j,k,jj;



#if(STOREFLUXSTATE)
  //if(ptrgeom->p==CENT) (GODMARK: assume interpline yrealin always at CENT, as it is now)

  // retrieve state from global array
  //  i=ptrgeom->i;
  //  j=ptrgeom->j;
  //  k=ptrgeom->k;

  *qptr=&(GLOBALMAC(fluxstatecent,MYII,MYJJ,MYKK)); // note that fluxstatecent for interpline has more than needed for interpline since was computed for flux+interpline and flux>interpline for all terms
#else
  // compute state right now
  pureget_stateforinterpline(pr, ptrgeom, *qptr);
#endif


  return (0);
}


/// wrapper for choosing to get state by computing it or from global array
/// pointer assign (overrides the existing default assignment of the pointer address in the calling function, so assumes calling function only uses pointer form subsequently)
int get_stateforglobalwavespeeds(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr)
{
  int pureget_stateforglobalwavespeeds(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  //  int i,j,k,jj;



#if(STOREFLUXSTATE)
  //if(ptrgeom->p==CENT) (GODMARK: assume interpline yrealin always at CENT, as it is now)

  // retrieve state from global array
  //  i=ptrgeom->i;
  //  j=ptrgeom->j;
  //  k=ptrgeom->k;

  *qptr=&(GLOBALMAC(fluxstatecent,MYII,MYJJ,MYKK)); // note that fluxstatecent for globalwavespeds has more than needed for interpline since was computed for flux+(globalwavespeeds) and flux>interpline for all terms
#else
  // compute state right now
  pureget_stateforglobalwavespeeds(pr, ptrgeom, *qptr);
#endif


  return (0);
}



/* add in source terms related to geometry to equations of motion */
int source_conn(FTYPE *pr, struct of_geom *ptrgeom,
                struct of_state *q,FTYPE *dU)
{
  VARSTATIC int i, j, k, l;
  VARSTATIC int pl,pliter;
  VARSTATIC FTYPE dUconn[NPR];
  VARSTATIC FTYPE todo[NPR];
  VARSTATIC FTYPE ftemp;
  VARSTATIC FTYPE mhd[NDIM][NDIM];
  VARSTATIC FTYPE mhdrad[NDIM][NDIM];
  VARSTATIC FTYPE flux[NDIM][NPR];
  int primtofullflux(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE (*flux)[NPR], FTYPE (*fluxabs)[NPR]);
  void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
  void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhdabs);
  VARSTATIC int myii,myjj,mykk;



  if((ANALYTICSOURCE)&&(defcoord==LOGRSINTH)&&(MCOORD==KSCOORDS)&&(EOMRADTYPE==EOMRADNONE)){
    // have both WHICHEOM==WITHGDET and WHICHEOM==WITHNOGDET
    // avoided if doing radiation since can't use mks_source_conn() for radiation without explicitly creating new function (can't just substitute q->uradcon/cov for q->ucon/cov)
    mks_source_conn(pr,ptrgeom, q, dU); // returns without geometry prefactor
    PLOOP(pliter,pl) dU[pl]*=ptrgeom->EOMFUNCMAC(pl);
  }
  else { // general non-analytic form, which uses an analytic or numerical origin for conn/conn2 (see set_grid.c)


    /////////////////////
    //
    // define first connection
    //
    // notice mhd_calc(), rather than primtoflux(), is used because this is a special connection only operated on the (at most) 4 energy-momentum EOMs.
    //
    // notice that mhd_calc_0 is used, which includes rest-mass since connection coefficient source term has rest-mass.  The rest-mass was only subtracted out of conserved quantity and flux term.
    // SUPERGODMARK: problem in non-rel gravity for below term
    //DLOOPA(j)  mhd_calc_0(pr, j, ptrgeom, q, mhd[j], NULL);
    // 


    // Get MHD stress-energy tensor
    DLOOPA(j) mhd_calc(pr, j, ptrgeom, q, mhd[j], NULL);

    // Get radiation stress-energy tensor (separate fluid from MHD fluid)
    if(EOMRADTYPE!=EOMRADNONE){
      DLOOPA(j) mhd_calc_rad(pr, j, ptrgeom, q, mhdrad[j], NULL);
    }


    /* contract mhd stress tensor with connection */
    // mhd^{dir}_{comp} = mhd^j_k
    // dU^{l} = T^j_k C^{k}_{l j}
    PLOOP(pliter,pl) dUconn[pl]=0;


#if(MCOORD!=CARTMINKMETRIC)
    myii=ptrgeom->i;
    myjj=ptrgeom->j;
    mykk=ptrgeom->k;
#else
    myii=0;
    myjj=0;
    mykk=0;
#endif




    // For MHD and radiation stress-energy tensor
    // Group these since conn[] repeats so should be easier on memory accesses
    if(EOMRADTYPE!=EOMRADNONE){
      FTYPE connterm;
      for(l=0;l<NDIM;l++)  DLOOP(j,k){
          connterm=GLOBALMETMACP0A3(conn,myii,myjj,mykk,k,l,j);
          dUconn[UU+l]    += mhd[j][k]    * connterm;
          dUconn[URAD0+l] += mhdrad[j][k] * connterm;
        }
    }
    else{
      // For MHD stress-energy tensor
      for(l=0;l<NDIM;l++)  DLOOP(j,k) dUconn[UU+l] += mhd[j][k] * GLOBALMETMACP0A3(conn,myii,myjj,mykk,k,l,j);
    }


#if(REMOVERESTMASSFROMUU==2)
    // then need to add-in density term to source (used to avoid non-rel problems)
    // GODMARK: for no non-rel problems one must make sure conn^t_{lj} is computed accurately.
    for(l=0;l<NDIM;l++)  DLOOPA(j){
        dUconn[UU+l] += - pr[RHO] * q->ucon[j] * GLOBALMETMACP0A3(conn,myii,myjj,mykk,TT,l,j);
      }
#endif


    ////////////////////////////
    //
    // use first connection
    //
    // true as long as UU->U3 are not added/subtracted from other EOMs
    // Note that UtoU is ignorant of this connection term and additions/subtractions of various EOMs to/from another must account for this connection here.
    PLOOP(pliter,pl) dU[pl]+=ptrgeom->EOMFUNCMAC(pl)*dUconn[pl]; 


 

    ///////////////////
    //
    // second connection
    //

#if(WHICHEOM!=WITHGDET)

    // d_t(f T^t_\nu) = -d_x^j (f F^j_\nu) + (f F^j_\nu)(d_x^j(\ln(f/\detg))) + f T^\lambda_\kappa \Gamma^\kappa_{\nu\lambda}

    // conn2_j = (d_x^j(\ln(f/\detg)))

    // deal with second connection.  Now can use general flux as defined by primtofullflux since all EOMs operate similarly with respect to second connection
    primtofullflux(UEVOLVE,pr,q,ptrgeom,flux,NULL); // returns with geometry prefactor (f F^j_\nu)


    // todo = whether that EOM has the NOGDET form.  If so, then need 2nd connection.  Do this instead of different connection2 for each EOM since each spatial component is actually the same.
    PLOOP(pliter,pl) todo[pl] = (nogdetlist[pl]>0 ? 1 : 0);

    // conn2 is assumed to take care of sign
    // conn2 has geom->e and normal d(ln(gdet)) factors combined

    if(REMOVERESTMASSFROMUU){
      if(todo[RHO]!=todo[UU]){
        dualfprintf(fail_file,"Mixed form of REMOVERESTMASSFROMUU and NOGDET's not allowed.\n");
        myexit(1);
      }
    }

    //////////////////////////////////
    //
    // notice that we assume equations are differenced first, then one manipulates the connections.
    // Thus, the manipulation of connections applies to final form of EOMs afer differencing or adding.
    // This agrees with how code evolves conserved quantities, so simplest choice to make.
    //
    // Alternatively stated, UtoU() controls this ordering issue.

    PLOOP(pliter,pl) DLOOPA(j) dU[pl] += todo[pl]*(flux[j][pl] * GLOBALMETMACP0A1(conn2,myii,myjj,mykk,j)); // (f F^j_\nu) conn2_j



#endif // end if (WHICHEOM!=WITHGDET)

  }



  /* done! */
  return (0);
}



/// load local geometry into structure geom
/// fills struct of_geom with PRIMECOORDS metric data
#if(NEWMETRICSTORAGE==0)
void get_geometry_old(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int j, k;
  void assign_eomfunc(struct of_geom *geom, FTYPE *EOMFUNCNAME);
  VARSTATIC int myii,myjj,mykk;
  VARSTATIC int pl,pliter;
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=ii;
  myjj=jj;
  mykk=kk;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  

#if(GETGEOMUSEPOINTER==0)
  //  DLOOP(j,k) geom->gcov[GIND(j,k)] = GLOBALMETMACP1A2(gcov,pp,myii,myjj,mykk,j,k);
  //DLOOP(j,k) geom->gcon[GIND(j,k)] = GLOBALMETMACP1A2(gcon,pp,myii,myjj,mykk,j,k);
  // let's vectorize it
  for(j=0;j<=NDIM*NDIM-1;j++){
    geom->gcon[GIND(0,j)] = GLOBALMETMACP1A2(gcon,pp,myii,myjj,mykk,0,j);
    geom->gcov[GIND(0,j)] = GLOBALMETMACP1A2(gcov,pp,myii,myjj,mykk,0,j);
  }
  for(j=0;j<NDIM;j++){
    geom->gcovpert[j]=GLOBALMETMACP1A1(gcovpert,pp,myii,myjj,mykk,j);
  }
#else
  // let's pointer it, faster by a bit
  geom->gcov=(FTYPE (*)[NDIM])(&(GLOBALMETMACP1A2(gcov,pp,myii,myjj,mykk,0,0))); // pointer

  geom->gcovpert=(FTYPE (*))(&(GLOBALMETMACP1A0(gcovpert,pp,myii,myjj,mykk))); // pointer

  geom->gcon=(FTYPE (*)[NDIM])(&(GLOBALMETMACP1A2(gcon,pp,myii,myjj,mykk,0,0))); // pointer
#endif

  geom->gdet = GLOBALMETMACP1A0(gdet,pp,myii,myjj,mykk);

  // get eomfunc (see eomfunc_func() and lngdet_func_mcoord() in metric.c)
  assign_eomfunc(geom,GLOBALMETMACP1A0(eomfunc,pp,myii,myjj,mykk));

  geom->alphalapse = GLOBALMETMACP1A0(alphalapse,pp,myii,myjj,mykk);

  ///////////
  //
  // these should be real grid i,j,k positions
  geom->i = ii;
  geom->j = jj;
  geom->k = kk;
  geom->p = pp;
}

// load local geometry into structure geom
/// fills struct of_geom with PRIMECOORDS metric data
void get_geometry_gdetonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int j, k;
  void assign_eomfunc(struct of_geom *geom, FTYPE *EOMFUNCNAME);
  VARSTATIC int myii,myjj,mykk;
  VARSTATIC int pl,pliter;
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=ii;
  myjj=jj;
  mykk=kk;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  

  geom->gdet = GLOBALMETMACP1A0(gdet,pp,myii,myjj,mykk);


  ///////////
  //
  // these should be real grid i,j,k positions
  geom->i = ii;
  geom->j = jj;
  geom->k = kk;
  geom->p = pp;


}


/// load local geometry into structure geom
/// fills struct of_geom with PRIMECOORDS metric data
void get_geometry_geomeonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int pl,pliter;
  //  void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom);
  //  void get_geometry_gdetonly(int ii, int jj, int kk, int pp, struct of_geom *geom);
  


#if(WHICHEOM!=WITHGDET)
  // more complicated ptrgeom->e
  get_geometry(ii, jj, kk, pp, geom);

#else
  // then ptrgeom->e is just ptrgeom->g
  get_geometry_gdetonly(ii, jj, kk, pp, geom);
  PALLLOOP(pl) geom->EOMFUNCMAC(pl)=geom->gdet;
#endif


}


#endif // end if old metric method







/* load local geometry into structure geom */
void get_allgeometry(int ii, int jj, int kk, int pp, struct of_allgeom *allgeom, struct of_geom *ptrgeom)
{
  VARSTATIC int j, k;
  VARSTATIC int pl,pliter;
  

  get_geometry(ii, jj, kk, pp, ptrgeom);

  allgeom->i=ptrgeom->i;
  allgeom->j=ptrgeom->j;
  allgeom->k=ptrgeom->k;
  allgeom->p=ptrgeom->p;

#if(GETGEOMUSEPOINTER==0 || NEWMETRICSTORAGE==1)
  for(j=0;j<=NDIM*NDIM-1;j++){
    allgeom->gcov[GIND(0,j)]=ptrgeom->gcov[GIND(0,j)];
    allgeom->gcon[GIND(0,j)]=ptrgeom->gcon[GIND(0,j)];
  }
  for(j=0;j<NDIM;j++){
    allgeom->gcovpert[j]=ptrgeom->gcovpert[j];
  }
#else
  allgeom->gcov=ptrgeom->gcov;
  allgeom->gcon=ptrgeom->gcon;
  allgeom->gcovpert=ptrgeom->gcovpert;
#endif

  allgeom->gdet=ptrgeom->gdet;
  PLOOP(pliter,pl) allgeom->EOMFUNCMAC(pl)=ptrgeom->EOMFUNCMAC(pl);
  allgeom->alphalapse = ptrgeom->alphalapse;

  PALLLOOP(pl) allgeom->IEOMFUNCNOSINGMAC(pl) = ptrgeom->IEOMFUNCNOSINGMAC(pl);
  allgeom->igdetnosing = ptrgeom->igdetnosing;


  coord_ijk(ii, jj, kk, pp, allgeom->X);
  bl_coord_ijk( ii, jj, kk, pp , allgeom->V );
  dxdxprim_ijk( ii, jj, kk, pp , allgeom->dxdxp);


}


///GLOBALMETMACP1A0(eomfuncarray,pp,ii,jj,kk);
/// stationary factor in equation of motion that multiplies T^\mu_\nu
void assign_eomfunc(struct of_geom *geom, FTYPE *EOMFUNCNAME)
{
  int pliter,pl;

#if(WHICHEOM==WITHGDET)
  return; // now nothing to do since using new EOMFUNCMAC stuff
#endif

  // now set each EOM
  PLOOP(pliter,pl){
    if(nogdetlist[pl]==0) geom->EOMFUNCMAC(pl)=geom->gdet;
    else geom->EOMFUNCMAC(pl)=EOMFUNCASSIGN(pl);
  }

  if(NEWMETRICSTORAGE){
    // assign inverses (should already be set if PRIMECOORDS)
    // this only changes (overrides) to use igdetnosing since otherwise already should be set
    
    // now set each EOM
    PLOOP(pliter,pl){
      if(nogdetlist[pl]==0) geom->IEOMFUNCNOSINGMAC(pl)=geom->igdetnosing;
      // else will use igdet
    }
  } // end if NEWMETRICSTORAGE==1

}



/*
Kraken stopped when generated core dump:

Core was generated by `/lustre/scratch/jmckinne/rad1b//grmhd 32 8 4 1 0'.
Program terminated with signal 11, Segmentation fault.
#0  0x00000000005359fb in current_doprecalc (which=-131040, p=0xfffffffffffffffc) at phys.c:549
549           current_precalc(which,ptrgeom,&q,Dt,GLOBALMAC(cfaraday,i,j,k));
(gdb) bt
#0  0x00000000005359fb in current_doprecalc (which=-131040, p=0xfffffffffffffffc) at phys.c:549
#1  0x000000000054d57d in step_ch_full (truestep=-131040, prim=0xfffffffffffffffc, pstag=0xfffffffffffffffc, ucons=0xfffffffffffffffc, vpot=0xffffffffffffff40, Bhat=0xffffffffffff3400, pl_ct=0x302ca00, pr_ct=0x1415d20, F1=0x1f2f968, F2=0x2e4cd68, F3=0xb84cc8, Atemp=0x196dca0, 
    uconstemp=0xa60c60) at step_ch.c:51
#2  0x00000000005059a6 in main (argc=32700776, argv=0x302ca00) at main.c:115

*/


/// setup faraday for either instrinsic use (which==CURTYPEFARADAY) or for computing the currents
int current_doprecalc(int which, FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{


#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif




#pragma omp parallel OPENMPGLOBALPRIVATEFORSTATEANDGEOM // requires full state (OPTMARK: But shouldn't have to since just need field part of state)
  {
    int i,j,k;
    int idel, jdel,kdel;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    struct of_state q;
    FTYPE Dt;
    int face;
    OPENMP3DLOOPVARSDEFINE;


    if(WHICHCURRENTCALC==CURRENTCALC0){ // then need to save 2 times
      // which==1,2 should be using face primitives, but probably not
      if(which==CURTYPET){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt*0.5; }
      else if(which==CURTYPEX){ face=FACE1; idel=1; jdel=0; kdel=0; Dt=dt; }
      else if(which==CURTYPEY){ face=FACE2; idel=0; jdel=1; kdel=0; Dt=dt; }
      else if(which==CURTYPEZ){ face=FACE3; idel=0; jdel=0; kdel=1; Dt=dt; }
      else if(which==CURTYPEFARADAY){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt; }
    }
    else if(WHICHCURRENTCALC==CURRENTCALC1){ // time and space centered on present time (best method)
      face=CENT;
      idel=0;
      jdel=0;
      kdel=0;
      Dt=dt;
    }
    else if(WHICHCURRENTCALC==CURRENTCALC2){ // then need to save 2 times
      if(which==CURTYPET){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt*0.5; }
      else if(which==CURTYPEX){ face=CENT; idel=1; jdel=0; kdel=0; Dt=dt; }
      else if(which==CURTYPEY){ face=CENT; idel=0; jdel=1; kdel=0; Dt=dt; }
      else if(which==CURTYPEZ){ face=CENT; idel=0; jdel=0; kdel=1; Dt=dt; }
      else if(which==CURTYPEFARADAY){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt; }
    }
    else{
      dualfprintf(fail_file,"No such WHICHCURRENTCALC=%d\n",WHICHCURRENTCALC);
    }

    // assume no other conditions GODMARK

    //  COMPFULLLOOP{
    OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      get_geometry(i, j, k, face, ptrgeom);
      MYFUN(get_state(MAC(p,i,j,k), ptrgeom, &q),"phys.c:current_doprecalc()", "get_state()", 1);
      current_precalc(which,ptrgeom,&q,Dt,GLOBALMAC(cfaraday,i,j,k));
    }
  }// end parallel region



  return(0);
}


/// which: 0: somewhere at half step
/// which: 1: doing flux calculation in x1 direction, full step
/// which: 2: doing flux calculation in x2 direction, full step
/// faraday[0-3][0-3]=[0=mid-time,1=radial edge final time, 2= theta edge final time, 3=old mid-time][whatever 3 things needed to compute relevant current]
void current_precalc(int which, struct of_geom *geom, struct of_state *q, FTYPE Dt,FTYPE (*faraday)[3])
{
  // assume outside loop is like flux, from 0..N in r for r, 0..N in h for h.  And must wait till 2nd timestep before computing the current since need time differences

#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif


  if(which==CURTYPET){ // for d/dt

    // assume got here when DT==dt/2 and geom and state set at zone center
    // first save old calculation
    if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2)){ // then need to save 2 times
      faraday[4][0]=faraday[0][0];
      faraday[4][1]=faraday[0][1];
      faraday[4][2]=faraday[0][2];
    }
    else if(WHICHCURRENTCALC==CURRENTCALC1){ // then need to save 3 times and have [4] as 2 times ago
      // 2 times ago
      faraday[4][0]=faraday[5][0];
      faraday[4][1]=faraday[5][1];
      faraday[4][2]=faraday[5][2];

      // 1 time ago
      faraday[5][0]=faraday[0][0];
      faraday[5][1]=faraday[0][1];
      faraday[5][2]=faraday[0][2];
    }
    // now calculate new version

    // Need F^{xt} F^{yt} F^{zt}
    faraday[0][0]=-1.0/geom->gdet * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt}
    faraday[0][1]=-1.0/geom->gdet * (q->bcov[3]*q->ucov[1]-q->bcov[1]*q->ucov[3]); // f^{ht}
    faraday[0][2]=-1.0/geom->gdet * (-q->bcov[2]*q->ucov[1]+q->bcov[1]*q->ucov[2]); // f^{pt}
  }
  else if(which==CURTYPEX){// for d/dx

    // Need F^{tx} F^{yx} F^{zx}

    // assume got here with DT=dt and geom and state at radial zone edge
    faraday[1][0]=-1.0/geom->gdet * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}=-f^{rt}
    faraday[1][1]=-1.0/geom->gdet * (-q->bcov[3]*q->ucov[0]+q->bcov[0]*q->ucov[3]); // f^{hr}
    faraday[1][2]=-1.0/geom->gdet * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{pr}
  }
  else if(which==CURTYPEY){ // for d/dy

    // Need F^{ty} F^{xy} F^{zy}

    // assume got here with DT=dt and geom and state at theta zone edge
    faraday[2][0]=-1.0/geom->gdet * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}=-f^{ht}
    faraday[2][1]=-1.0/geom->gdet * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}=-f^{hr}
    faraday[2][2]=-1.0/geom->gdet * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{ph}
  }
  else if(which==CURTYPEZ){ // for d/dz

    // Need F^{tz} F^{xz} F^{yz}

    faraday[3][0]=1.0/geom->gdet * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt} // f^{tr}=-f^{rt}
    faraday[3][1]=1.0/geom->gdet * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{rp}=-f^{pr}
    faraday[3][2]=1.0/geom->gdet * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{hp}= -f^{ph}

  }
  else if(which==CURTYPEFARADAY){
    // DT==dt, but zone center
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,0)=-1.0/geom->gdet * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,1)=-1.0/geom->gdet * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,2)=-1.0/geom->gdet * (q->bcov[2]*q->ucov[1]-q->bcov[1]*q->ucov[2]); // f^{tp}
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,3)=-1.0/geom->gdet * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,4)=-1.0/geom->gdet * (-q->bcov[2]*q->ucov[0]+q->bcov[0]*q->ucov[2]); // f^{rp}
    GLOBALMACP0A1(fcon,geom->i,geom->j,geom->k,5)=-1.0/geom->gdet * (q->bcov[1]*q->ucov[0]-q->bcov[0]*q->ucov[1]); // f^{hp}
  }
  else{
    dualfprintf(fail_file,"No such which in current_doprecalc()\n");
  }
}


/// choose type of current calculation
void current_calc(FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3])
{
  void current_calc_0(FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3]);
  void current_calc_1(int which, FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3]);


#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif

  if(WHICHCURRENTCALC==CURRENTCALC0){
    current_calc_0(cfaraday);
  }
  else if((WHICHCURRENTCALC==CURRENTCALC1)||(WHICHCURRENTCALC==CURRENTCALC2)){
    current_calc_1(WHICHCURRENTCALC,cfaraday);
  }

}



/// the current is calculated to end up at the zone and time edge
/// point is to have J^t at same time as rest of J's, although different spacial points
/// OPENMPMARK: Assume this function not called by multiple threads
void current_calc_0(FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3])
{
  static FTYPE lastdt;
  static int calls=0;



#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif


#pragma omp parallel 
  {
    int i,j,k;
    struct of_geom geomtdontuse;
    struct of_geom *ptrgeomt=&geomtdontuse;

    struct of_geom geomrdontuse;
    struct of_geom *ptrgeomr=&geomrdontuse;

    struct of_geom geomhdontuse;
    struct of_geom *ptrgeomh=&geomhdontuse;

    struct of_geom geompdontuse;
    struct of_geom *ptrgeomp=&geompdontuse;

    struct of_geom geomrp1dontuse;
    struct of_geom *ptrgeomrp1=&geomrp1dontuse;

    struct of_geom geomhp1dontuse;
    struct of_geom *ptrgeomhp1=&geomhp1dontuse;

    struct of_geom geompp1dontuse;
    struct of_geom *ptrgeompp1=&geompp1dontuse;

    FTYPE idtc,idx1,idx2,idx3;
    FTYPE timeterm[NDIM];
    OPENMP3DLOOPVARSDEFINE;    


    idtc=2.0/(lastdt+dt);
    idx1=1.0/dx[1];
    idx2=1.0/dx[2];
    idx3=1.0/dx[3];

    /////////COMPLOOPP1{ // largest possible loop for this differencing (could isolate directions)
    OPENMP3DLOOPSETUP(INP11,OUTP11,INP12,OUTP12,INP13,OUTP13);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      get_geometry(i,j,k,CENT,ptrgeomt);
      get_geometry(i,j,k,FACE1,ptrgeomr);
      get_geometry(i,j,k,FACE2,ptrgeomh);
      get_geometry(i,j,k,FACE3,ptrgeomp);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(ip1mac(i),j,k,FACE1,ptrgeomrp1);
      get_geometry(i,jp1mac(j),k,FACE2,ptrgeomhp1);
      get_geometry(i,j,kp1mac(k),FACE3,ptrgeompp1);

      if(calls>0){ // since need 2 times
        timeterm[0]=0;
        timeterm[1]=+(GLOBALMACP0A2(cfaraday,i,j,k,0,0)-GLOBALMACP0A2(cfaraday,i,j,k,4,0))*idtc; // F^{rt},t
        timeterm[2]=+(GLOBALMACP0A2(cfaraday,i,j,k,0,1)-GLOBALMACP0A2(cfaraday,i,j,k,4,1))*idtc; // F^{ht},t
        timeterm[3]=+(GLOBALMACP0A2(cfaraday,i,j,k,0,2)-GLOBALMACP0A2(cfaraday,i,j,k,4,2))*idtc; // F^{pt},t
      }
      else{
        timeterm[0]=0;
        timeterm[1]=0;
        timeterm[2]=0;
        timeterm[3]=0;
      }
      
      // J^t = F^{tr},r + F^{th},h + F^{tp},p
      GLOBALMACP0A1(jcon,i,j,k,0)=
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,0)-ptrgeomr->gdet*GLOBALMACP0A2(cfaraday,i,j,k,1,0))*idx1 // F^{tr},r
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,0)-ptrgeomh->gdet*GLOBALMACP0A2(cfaraday,i,j,k,2,0))*idx2 // F^{th},h
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,0)-ptrgeomp->gdet*GLOBALMACP0A2(cfaraday,i,j,k,3,0))*idx3 // F^{tp},p
        ;
      
      // J^r = F^{rt},t + F^{rh},h + F^{rp},p
      GLOBALMACP0A1(jcon,i,j,k,1)=
        +timeterm[1]
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,1)-ptrgeomh->gdet*GLOBALMACP0A2(cfaraday,i,j,k,2,1))*idx2 // F^{rh},h
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,1)-ptrgeomp->gdet*GLOBALMACP0A2(cfaraday,i,j,k,3,1))*idx3 // F^{rp},p
        ;
      
      // J^h = F^{ht},t + F^{hr},r + F^{hp},p
      GLOBALMACP0A1(jcon,i,j,k,2)=
        +timeterm[2]
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,1)-ptrgeomr->gdet*GLOBALMACP0A2(cfaraday,i,j,k,1,1))*idx1 // F^{hr},r
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,2)-ptrgeomp->gdet*GLOBALMACP0A2(cfaraday,i,j,k,3,2))*idx3 // F^{hp},p
        ;
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      GLOBALMACP0A1(jcon,i,j,k,3)=
        +timeterm[3]
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,2)-ptrgeomr->gdet*GLOBALMACP0A2(cfaraday,i,j,k,1,2))*idx1 // F^{pr},r
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,2)-ptrgeomh->gdet*GLOBALMACP0A2(cfaraday,i,j,k,2,2))*idx2 // F^{ph},h
        ;

    }// end loops


  }// end over parallel region


  // static (should be outside parallel region)
  calls++;
  lastdt=dt;

}



/// the current is calculated to end up at the zone and time edge
/// point is to have J^t at same time as rest of J's and at same spatial location
/// the temporal value of things is obtained at same time of everything else by fitting to parabola in time
/// OPENMPMARK: Assume this function not called by multiple threads
void current_calc_1(int which, FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3])
{
  static FTYPE lastdt;
  static long long calls=0;


#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif



#pragma omp parallel 
  {
    int i,j,k;
    struct of_geom geomtdontuse;
    struct of_geom *ptrgeomt=&geomtdontuse;

    struct of_geom geomcdontuse;
    struct of_geom *ptrgeomc=&geomcdontuse;

    struct of_geom geomrp1dontuse;
    struct of_geom *ptrgeomrp1=&geomrp1dontuse;

    struct of_geom geomhp1dontuse;
    struct of_geom *ptrgeomhp1=&geomhp1dontuse;

    struct of_geom geompp1dontuse;
    struct of_geom *ptrgeompp1=&geompp1dontuse;

    struct of_geom geomrm1dontuse;
    struct of_geom *ptrgeomrm1=&geomrm1dontuse;

    struct of_geom geomhm1dontuse;
    struct of_geom *ptrgeomhm1=&geomhm1dontuse;

    struct of_geom geompm1dontuse;
    struct of_geom *ptrgeompm1=&geompm1dontuse;

    FTYPE idtc,idx1,idx2,idx3;
    FTYPE AA,BB;
    FTYPE dtl,dtr;
    FTYPE fl,fr,f0,derf;
    FTYPE timeterm[NDIM];
    OPENMP3DLOOPVARSDEFINE;    


    if(which==CURRENTCALC1){
      dtl=lastdt;
      dtr=dt;
    }
    else if(which==CURRENTCALC2){
      idtc=2.0/(lastdt+dt);
    }
  
    idx1=1.0/(2.0*dx[1]);
    idx2=1.0/(2.0*dx[2]);
    idx3=1.0/(2.0*dx[3]);


  
    /////////COMPLOOPP1{ // largest possible loop for this differencing (could isolate directions)
    OPENMP3DLOOPSETUP(INP11,OUTP11,INP12,OUTP12,INP13,OUTP13);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);


      get_geometry(i,j,k,CENT,ptrgeomt);

      get_geometry(i,j,k,CENT,ptrgeomc);

      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(ip1mac(i),j,k,CENT,ptrgeomrp1);
      get_geometry(i,jp1mac(j),k,CENT,ptrgeomhp1);
      get_geometry(i,j,kp1mac(k),CENT,ptrgeompp1);

      get_geometry(im1mac(i),j,k,CENT,ptrgeomrm1);
      get_geometry(i,jm1mac(j),k,CENT,ptrgeomhm1);
      get_geometry(i,j,km1mac(k),CENT,ptrgeompm1);

      if(which==CURRENTCALC1){
      
        if(calls>=2){ // since need 3 times to properly time center without having to worry about what RK is doing

          timeterm[0]=0;
 
          fl=GLOBALMACP0A2(cfaraday,i,j,k,4,0);
          f0=GLOBALMACP0A2(cfaraday,i,j,k,5,0);
          fr=GLOBALMACP0A2(cfaraday,i,j,k,0,0);
          AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          derf=2.0*AA*dtr+BB;
 
          timeterm[1]=derf;
 
          fl=GLOBALMACP0A2(cfaraday,i,j,k,4,1);
          f0=GLOBALMACP0A2(cfaraday,i,j,k,5,1);
          fr=GLOBALMACP0A2(cfaraday,i,j,k,0,1);
          AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          derf=2.0*AA*dtr+BB;
 
          timeterm[2]=derf;

          fl=GLOBALMACP0A2(cfaraday,i,j,k,4,2);
          f0=GLOBALMACP0A2(cfaraday,i,j,k,5,2);
          fr=GLOBALMACP0A2(cfaraday,i,j,k,0,2);
          AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
          derf=2.0*AA*dtr+BB;
 
          timeterm[3]=derf;
        }
        else{
          timeterm[0]=0;
          timeterm[1]=0;
          timeterm[2]=0;
          timeterm[3]=0;
        }
      }
      else if(which==CURRENTCALC2){
        // the current is calculated to end up at the zone and time edge
        // point is to have J^t at same time as rest of J's, and same spatial points.  time for J is one timestep back for all components.
        // like current_calc_0 in time and current_calc_1 in space
      
        if(calls>0){ // since need 2 times
          timeterm[0]=0;
          timeterm[1]=(GLOBALMACP0A2(cfaraday,i,j,k,0,0)-GLOBALMACP0A2(cfaraday,i,j,k,4,0))*idtc; // F^{rt},t
          timeterm[2]=(GLOBALMACP0A2(cfaraday,i,j,k,0,1)-GLOBALMACP0A2(cfaraday,i,j,k,4,1))*idtc; // F^{ht},t
          timeterm[3]=(GLOBALMACP0A2(cfaraday,i,j,k,0,2)-GLOBALMACP0A2(cfaraday,i,j,k,4,2))*idtc; // F^{pt},t
        }
        else{
          timeterm[0]=0;
          timeterm[1]=0;
          timeterm[2]=0;
          timeterm[3]=0;
        }

      }

      // similar to other current_calc's except position of current and spacing of derivatives

      // J^t = F^{tr},r + F^{th},h + F^{tp},p
      GLOBALMACP0A1(jcon,i,j,k,0)=
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,0)-ptrgeomrm1->gdet*GLOBALMACP0A2(cfaraday,im1mac(i),j,k,1,0))*idx1 // F^{tr},r
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,0)-ptrgeomhm1->gdet*GLOBALMACP0A2(cfaraday,i,jm1mac(j),k,2,0))*idx2 // F^{th},h
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,0)-ptrgeompm1->gdet*GLOBALMACP0A2(cfaraday,i,j,km1mac(k),3,0))*idx3 // F^{tp},p
        ;

      // J^r = F^{rt},t + F^{rh},h + F^{rp},p
      GLOBALMACP0A1(jcon,i,j,k,1)=
        +timeterm[1]
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,1)-ptrgeomhm1->gdet*GLOBALMACP0A2(cfaraday,i,jm1mac(j),k,2,1))*idx2 // F^{rh},h
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,1)-ptrgeompm1->gdet*GLOBALMACP0A2(cfaraday,i,j,km1mac(k),3,1))*idx3 // F^{rp},p
        ;
      
      // J^h = F^{ht},t + F^{hr},r + F^{hp},p
      GLOBALMACP0A1(jcon,i,j,k,2)=
        +timeterm[2]
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,1)-ptrgeomrm1->gdet*GLOBALMACP0A2(cfaraday,im1mac(i),j,k,1,1))*idx1 // F^{hr},r
        +1./ptrgeomt->gdet*(ptrgeompp1->gdet*GLOBALMACP0A2(cfaraday,i,j,kp1mac(k),3,2)-ptrgeompm1->gdet*GLOBALMACP0A2(cfaraday,i,j,km1mac(k),3,2))*idx3 // F^{hp},p
        ;
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      GLOBALMACP0A1(jcon,i,j,k,3)=
        +timeterm[3]
        +1./ptrgeomt->gdet*(ptrgeomrp1->gdet*GLOBALMACP0A2(cfaraday,ip1mac(i),j,k,1,2)-ptrgeomrm1->gdet*GLOBALMACP0A2(cfaraday,im1mac(i),j,k,1,2))*idx1 // F^{pr},r
        +1./ptrgeomt->gdet*(ptrgeomhp1->gdet*GLOBALMACP0A2(cfaraday,i,jp1mac(j),k,2,2)-ptrgeomhm1->gdet*GLOBALMACP0A2(cfaraday,i,jm1mac(j),k,2,2))*idx2 // F^{ph},h
        ;

    
    }// end loops

  }// end parallel region



  // static (should be outside parallel region)
  calls++;
  lastdt=dt;

}




/// compute vorticity
///#define COMPLOOPVORT COMPLOOPP11 COMPLOOPP12 COMPLOOPP13
/// use (*p)[][][U1,U2,U3] to obtain (*pvort)[][][whichpl]
int compute_vorticity(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE (*pvort)[NSTORE2][NSTORE3][NPR],int whichpl)
{


#pragma omp parallel 
  {
    int i,j,k;
    FTYPE X[NDIM],V[NDIM];
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    FTYPE dxdxp[NDIM][NDIM];
    FTYPE vxm,vxp,vym,vyp;
    OPENMP3DLOOPVARSDEFINE;    
  
  
    //  COMPLOOPVORT{
    OPENMP3DLOOPSETUP(INP11,OUTP11,INP12,OUTP12,INP13,OUTP13);
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);

      bl_coord_ijk_2(i, jm1mac(j), k, CENT, X, V);
      get_geometry(i,jm1mac(j),k,CENT,ptrgeom) ;
      // dx/dx' where '=prim coords (i.e. nonuni coords)
      vxm = MACP0A1(p,i,jm1mac(j),k,U1)*sqrt(ptrgeom->gcov[GIND(1,1)]);

      bl_coord_ijk_2(i, jp1mac(j), k, CENT, X, V);
      get_geometry(i,jp1mac(j),k,CENT,ptrgeom) ;
      // dx/dx' where '=prim coords (i.e. nonuni coords)
      vxp = MACP0A1(p,i,jp1mac(j),k,U1)*sqrt(ptrgeom->gcov[GIND(1,1)]);



      bl_coord_ijk_2(im1mac(i), j, k, CENT, X, V);
      get_geometry(im1mac(i),j,k,CENT,ptrgeom) ;
      // dx/dx' where '=prim coords (i.e. nonuni coords)
      vym = MACP0A1(p,im1mac(i),j,k,U2)*sqrt(ptrgeom->gcov[GIND(2,2)]);

      bl_coord_ijk_2(ip1mac(i), j, k, CENT, X, V);
      get_geometry(ip1mac(i),j,k,CENT,ptrgeom) ;
      // dx/dx' where '=prim coords (i.e. nonuni coords)
      vyp = MACP0A1(p,ip1mac(i),j,k,U2)*sqrt(ptrgeom->gcov[GIND(2,2)]);

      // center
      bl_coord_ijk_2(i, j, k, CENT, X, V);
      get_geometry(i,j,k,CENT,ptrgeom) ;
      dxdxprim_ijk(i,j,k, CENT, dxdxp);

      MACP0A1(pvort,i,j,k,whichpl) = (vyp-vym)/(2.0*dx[1])/dxdxp[1][1];

      MACP0A1(pvort,i,j,k,whichpl) += -(vxp-vxm)/(2.0*dx[2])/dxdxp[2][2];

    }
  }// end parallel region


  return(0);

}
