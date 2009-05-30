
// pi: initial values at t=t0 to compute Ui
// pb: values used to compute flux/source
// pf: solution using flux(pb) from pi's Ui -> Uf

// pi, pb, and pf can all be the same since
// 1) pb used first on a stencil, not modified, to compute fluxes
// 2) pf=pi is assigned by value at each zone
// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
//
// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance_eno_u(int stage,
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
  SFTYPE dt4diag;

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
  CZLOOP {

    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    ////////////
    PLOOP pf[i][j][k] = pi[i][j][k];

    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom);

    // find dU(pb)
    //MYFUN(source(pb[i][j], &geom, i, j, dUcomp, dU,CUf[2]),"step_ch.c:advance()", "source", 1);
    //diag_source(i,j,dUcomp,dt4diag);

    // find Ui(pi)
    MYFUN(get_state(pi[i][j], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
    MYFUN(primtoU(pi[i][j], &q, &geom, Ui),"step_ch.c:advance()", "primtoU()", 1);
    
    // convert zone centered conserved quantities in Ui to zone averaged ones (atchekho)
    PLOOP uc[i][j][k] = Ui[k];
  }
  
  c2a( ua, uc );
  //now ua contains cell averaged values
  
  CZLOOP {
    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom); //added by atch

    // find dU(pb) //added by atch
    MYFUN(source(pb[i][j], &geom, i, j, dUcomp, dU,CUf[2]),"step_ch.c:advance()", "source", 1);
    diag_source(i,j,dUcomp,dt4diag);

    PLOOP Ui[k] = ua[i][j][k];  //replace with cell averaged values
  
    // initialize ulast and unew if very first time here
    if(stagenow==0) PLOOP ulast[i][j][k]=unew[i][j][k]=0.0;

    // find ulast==Uf and additional terms to unew
    flux2U(i,j,F1,F2,dU,CUf,Cunew,Ui, ulast[i][j],Uf,unew[i][j]);

    // invert U->p (but do conversion from averages to cell-centered quantities first!)
    if(stagenow==numstages-1){ // last call, so unew is cooked and ready to eat!
      // instead of unew  need to supply cell centered quantities
      PLOOP ua[i][j][k] = unew[i][j][k];
      //MYFUN(Utoprimgen(unew[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
    else{ // otherwise still iterating on primitives
      // instead of Uf  need to supply cell centered quantities
      PLOOP ua[i][j][k] = Uf[k];
      //MYFUN(Utoprimgen(Uf, &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
  }
  
  //now ua contains cell averaged conserved quantities, need to convert to cell-centered ones
  
  //convert to cell centered ones
  a2c( uc, ua );

  CZLOOP {
    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom); //added by atch

    // find dU(pb) //added by atch
    MYFUN(source(pb[i][j], &geom, i, j, dUcomp, dU,CUf[2]),"step_ch.c:advance()", "source", 1);
    diag_source(i,j,dUcomp,dt4diag);

    // invert U->p
    MYFUN(Utoprimgen(uc[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    
    // immediate local (i.e. 1-zone) fix
#if(FIXUP1ZONE)
      MYFUN(fixup1zone(pf[i][j],&geom,dt4diag),"fixup.c:fixup()", "fixup1zone()", 1);
#endif
  }// end CZLOOP

  //atch: fixup_zones: removed to avoid any fixup; not the proper way to solve: need a switch rather than 
  //      comment out
  //#if(!FIXUP1ZONE)
  //fixup(stage,pf,dt4diag);
  //#enedif  


  *ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2);

  trifprintf( "2f");

  return (0);
}







