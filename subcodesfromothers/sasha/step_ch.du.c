

// pi: initial values at t=t0 to compute Ui
// pb: values used to compute flux/source
// pf: solution using flux(pb) from pi's Ui -> Uf

// pi, pb, and pf can all be the same since
// 1) pb used first on a stencil, not modified, to compute fluxes
// 2) pf=pi is assigned by value at each zone
// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
//
// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)

// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance_eno_du(int stage,
	    FTYPE pi[][N2M][NPR],
	    FTYPE ulast[][N2M][NPR], 
	    FTYPE pb[][N2M][NPR],
	    FTYPE *CUf, FTYPE pf[][N2M][NPR],
	    FTYPE *Cunew, FTYPE unew[][N2M][NPR],
	    int stagenow, int numstages,
	    FTYPE *ndt)
{
  int i, j, k, sc;
  FTYPE ndt1, ndt2;
  FTYPE Uf[NPR], Ui[NPR], Ub[NPR], dU[NPR],dUcomp[NUMSOURCES][NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  FTYPE (*dUriemannavg)[N2M][NPR];
  FTYPE dUriemann[NPR], dUgeom[NPR];
  //  FTYPE *dUriemannavg[N2M][NPR];
  SFTYPE dt4diag;
  void flux2dUavg(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dUavg);
  void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *ulast, FTYPE *Uf, FTYPE *unew);


  dUriemannavg=ua;

  //   dUriemannavg=(FTYPE (*[N2M][NPR]))(&ua[0][0][0][0]);
  //   dUriemannavg=(FTYPE *[N2M][NPR])(&ua[2][2][0]);
  


  // tells diagnostics functions if should be accounting or not
  if(stagenow==numstages-1) dt4diag=dt; else dt4diag=-1.0;


  trifprintf( "#0f");
  // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
  MYFUN(fluxcalc(stage, pb, F1, 1, &ndt1),"step_ch.c:advance()", "fluxcalc", 1);
  MYFUN(fluxcalc(stage, pb, F2, 2, &ndt2),"step_ch.c:advance()", "fluxcalc", 2);
  fix_flux(F1, F2);
  flux_ct(stage, F1, F2);



  trifprintf( "1f");
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil


  /** now update pi to pf **/

  // get dU^{n+1}
  CZLOOP flux2dUavg(i,j,F1,F2,dUriemannavg[i][j]);


  // get U^{n+1}  
  CZLOOP {

    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    ////////////
    PLOOP pf[i][j][k] = pi[i][j][k];

    // initialize ulast and unew if very first time here
    if(stagenow==0) PLOOP ulast[i][j][k]=unew[i][j][k]=0.0;


    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom);

    // find Ui(pi)
    MYFUN(get_state(pi[i][j], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
    MYFUN(primtoU(pi[i][j], &q, &geom, Ui),"step_ch.c:advance()", "primtoU()", 1);


    // find dU(pb)
    MYFUN(source(pb[i][j], &geom, i, j, dUcomp, dUgeom, CUf[2]),"step_ch.c:advance()", "source", 1);
    diag_source(i,j,dUcomp,dt4diag);

    //convert to cell centered ones
    a2cij(1,i,j,dUriemannavg, dUriemann );


    dUtoU(dUgeom, dUriemann, CUf, Cunew, Ui, ulast[i][j], Uf, unew[i][j]);


    // invert U->p (but do conversion from averages to cell-centered quantities first!)
    if(stagenow==numstages-1){ // last call, so unew is cooked and ready to eat!
      // instead of unew  need to supply cell centered quantities
      MYFUN(Utoprimgen(unew[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
    else{ // otherwise still iterating on primitives
      // instead of Uf  need to supply cell centered quantities
      MYFUN(Utoprimgen(Uf, &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }

    // immediate local (i.e. 1-zone) fix
#if(FIXUP1ZONE)
      MYFUN(fixup1zone(pf[i][j],&geom,dt4diag),"fixup.c:fixup()", "fixup1zone()", 1);
#endif

  }
  

  //atch: fixup_zones: removed to avoid any fixup; not the proper way to solve: need a switch rather than 
  //      comment out
  //#if(!FIXUP1ZONE)
  //fixup(stage,pf,dt4diag);
  //#enedif  


  *ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2);

  trifprintf( "2f");

  return (0);
}



// convert point versions of U_i^{n} and dU -> U_i^{n+1} and other versions
void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *ulast, FTYPE *Uf, FTYPE *unew)
{
    int k;
    // finally assign new Uf and unew
    // store ulast to avoid recomputing U(pf) used later as pb for advance()
    PLOOP ulast[k]=Uf[k]   = CUf[0]*Ui[k] + CUf[1]*ulast[k] + CUf[2]*dt*(dUriemann[k]+dUgeom[k]);

    // how much of Ui, dU, and Uf to keep for final solution
    // ultimately unew is actual solution used to find final pf
    PLOOP unew[k] += Cunew[0]*Ui[k] + Cunew[1]*dt*(dUriemann[k]+dUgeom[k]) + Cunew[2]*Uf[k];



}

// get dUavg
void flux2dUavg(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dUavg)
{
  FTYPE idx1,idx2;
  int k;

#if(VOLUMEDIFF==0)
    idx1=1.0/dx[RR];
    idx2=1.0/dx[TH];
#else
    idx1=idxvol[i][j][RR];
    idx2=idxvol[i][j][TH];
#endif


#if(FLUXB==FLUXCD) // don't use volume reg. since differencing is large
    for(k=0;k<=U3;k++){
      dUavg[k]=(
	      -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
	      - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
	      );
    }
    k=B1;
    dUavg[k]=(
	   - (F2[i][j + 1][k] - F2[i][j - 1][k]) *idx2
	   );
    k=B2;
    dUavg[k]=(
	   - (F1[i+1][j][k] - F1[i-1][j][k]) *idx1
	   );
    for(k=B3;k<NPR;k++){
      dUavg[k]=(
	      -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
	      - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
	      );
    }
#else
    // other (normal) FLUXB methods
    PLOOP {
      dUavg[k]=(
	      -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
	      - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
	      );
    }
#endif




}






