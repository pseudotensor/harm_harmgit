
/*! \file step_ch.c
     \brief Code to take full RK step
*/


#include "decs.h"


static void setup_rktimestep(int truestep, int *numtimeorders,
                             FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR],
                             FTYPE (*pk)[NSTORE1][NSTORE2][NSTORE3][NPR],
                             FTYPE (*pii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*pbb[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*pff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],
                             FTYPE (*uii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*uff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*ucum[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],
                             FTYPE (*CUf)[NUMDTCUFS],FTYPE (*CUnew)[NUMDTCUFS]);




static int pre_stepch(int *dumpingnext, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
static int post_stepch(int *dumpingnext, FTYPE fullndt, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
static int step_ch(int truestep, int *dumpingnext, FTYPE *fullndt,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);
static int post_advance(int truestep, int *dumpingnext, int timeorder, int numtimeorders, int finalstep, SFTYPE boundtime, SFTYPE fluxtime, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);









/// take full step (called from main.c)
int step_ch_full(int truestep, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  FTYPE fullndt;
  // dumpingnext[0] = dumping just after this step
  // dumpingnext[1] = dumping just after the step after this step
  int dumpingnext[2];



  // things to do before taking full step
  if(truestep) pre_stepch(dumpingnext,prim);

  // take full step
  step_ch(truestep,dumpingnext, &fullndt,prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct, F1, F2, F3,Atemp,uconstemp);

  // things to do after taking full step
  if(truestep) post_stepch(dumpingnext, fullndt, prim, ucons);


  if(truestep){ // don't do if just passing through -- otherwise would end up looping endlessly!
    // add up contributions to flux through horizon (really inner boundary)
    if(DOEVOLVEMETRIC && (EVOLVEMETRICSUBSTEP==0 || EVOLVEMETRICSUBSTEP==2) ){
      compute_new_metric_longsteps(prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct, F1, F2, F3,Atemp,uconstemp);
    }
    
    // below must come after longstep update so that don't change before metric step taken!
    control_metric_update();
    
    postdt(); // here one can alter variables and try to restart, or implement any post step operations


  }
  
  
  // must check before MPI operation (since asymmetries would desynchronize cpus)
  MYFUN(error_check(2),"main.c","error_check",1);
  //^^ otherwise ok
  // eventually all cpus come here, either in failure mode or not, and cleanly tell others if failed/exit/dump/etc.


  
  
  
  return(0);
}



/// take step itself
int step_ch(int truestep, int *dumpingnext, FTYPE *fullndt,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  int step_ch_simplempi(int truestep, int *dumpingnext, FTYPE *fullndt,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);
  int step_ch_supermpi(int truestep, int *dumpingnext, FTYPE *fullndt,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);



#if(SIMULBCCALC==-1)
  MYFUN(step_ch_simplempi(truestep, dumpingnext, fullndt,prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct,F1,F2,F3,Atemp,uconstemp),"step_ch.c:step_ch()", "step_ch_simplempi()", 1);
#else
  MYFUN(step_ch_supermpi(truestep, dumpingnext, fullndt,prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct,F1,F2,F3,Atemp,uconstemp),"step_ch.c:step_ch()", "step_ch_supermpi()", 1);
#endif


  /* done! */
  return (0);
}



/// things to do before taking a full timestep
/// assume not called when dt==0.0
int pre_stepch(int *dumpingnext, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  


#if(PRODUCTION==0)
  trifprintf( "\n#step=%ld",realnstep);
#endif



  // if not doing diagnostics, then dumpingnext will be 0
  // first check if dumping next step
  // note dt is already set to correct value here
  dumpingnext[0]=(diag(FUTURE_OUT,t+dt,nstep+1,realnstep+1)==DOINGFUTUREOUT);
  dumpingnext[1]=(diag(FUTURE_OUT,t+dt+dt/2,nstep+2,realnstep+2)==DOINGFUTUREOUT); // estimate of t+dt+dt/2 good with 1<SAFEFACTOR<2



    
  if(dumpingnext[0] || dumpingnext[1]){
#if(CALCFARADAYANDCURRENTS)
    if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2)){
      // for this method, all J components are located back in time
      current_doprecalc(CURTYPEX,prim);
      current_doprecalc(CURTYPEY,prim);
      current_doprecalc(CURTYPEZ,prim);
      // get faraday in case used for time averaging or dumping it
      current_doprecalc(CURTYPEFARADAY,prim);
    }
#endif
  }


  return(0);


}




/// things to do after taking a full timestep
/// assume not called when dt==0.0
int post_stepch(int *dumpingnext, FTYPE fullndt, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{




  if(dumpingnext[0]){
#if(CALCFARADAYANDCURRENTS)
    // compute final faradays
    if(WHICHCURRENTCALC==CURRENTCALC1){
      // compute faraday components needed for time centering of J
      current_doprecalc(CURTYPET,prim);
      // J is located at this time
      current_doprecalc(CURTYPEX,prim);
      current_doprecalc(CURTYPEY,prim);
      current_doprecalc(CURTYPEZ,prim);
      current_doprecalc(CURTYPEFARADAY,prim); // for diagnostics
    }
    // compute current
    current_calc(GLOBALPOINT(cfaraday));
#endif
  }


#if(ACCURATEDIAG==0)
  // compute flux diagnostics
  // this doesn't exactly make conservation work -- should have in middle step point using full step.  But if PARA, no middle point that's exact.
  // think about where to put this
  // GODMARK : use of globals
  diag_flux(prim,F1, F2, F3, dt); // should use REAL dt, not within a timeorderd RK integration step
  //    horizon_flux(F1,dt); // subsumed
#endif

  // general flux only done on full steps since no requirement for accuracy and code can be expensive computationally
  diag_flux_general(prim,dt);// Should be full dt, not substep dt.




  if(DOONESTEPDUACCOUNTING){
    // do one-step accounting for fixup related things (failures, floors, fixups, etc.)
    // ucons is in UEVOLVE form (originates from unewglobal in step_ch_full() called in main.c)
    diag_fixup_allzones(prim, ucons);
  }



#if(PRODUCTION==0)
  trifprintf( "[d]");
#endif



  //////////////////////////
  //
  // increment time
  //
  //////////////////////////
  t += dt;
  tstepparti = tsteppartf = t;
    
  realnstep++;
    
  ////////////////////////////
  //
  // set next timestep
  //
  //
  // Notes:
  // 
  // dt is computed in a few different places for different steps of the calculation:
  // 1) flux.c (per dimension over all cells):
  //    A) dtij = cour * dx[dir] / ctop;
  //    B) if (dtij < *ndt){ *ndt = dtij; }
  //
  // 2) advance.c: prepare_globaldt() ("min" over each dimension and operator):
  //    A) wavedt = 1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3);
  //    B) *ndt = defcon * MIN(wavedt , accdt);
  //    C) etc...
  // 
  // 3) step_ch.c:  ("min" over each substep of RK steps):
  //    A) step_ch_simplempi(): if(ndt<lastndt) lastndt=ndt;
  //    B) HERE: to find global minimum over all CPUs
  //    C) HERE: to use safety factor
  //    D) HERE: constrain new dt if near end of calculation
  //
  // But, should do #1's MIN over all cells per CPU and #3's global MIN over all CPUs together for perfect MPI consistency.
  //  Further, should compute "wavedt" per cell since otherwise non-local dt_i can impact minimum dt.  Can lead to dt smaller by as much as 3X (with no correct numerical method effect) !
  // 
  ////////////////////////////
  // find global minimum value of ndt over all cpus
  mpifmin(&fullndt);
  // note that this assignment to fullndt doesn't get returned outside function
  if (fullndt > SAFE * dt)    fullndt = SAFE * dt;
  dt = fullndt;

  ///////////////////////////////////
  //
  // don't step beyond end of run
  //
  ////////////////////////////////////
  if(onemorestep){
    dualfprintf(fail_file,"Got onemorestep dt=%g\n",dt);
    // check if previous step was onemorestep==1
    reallaststep=1;
    dt=SMALL;
  }
  else{
    if (t + dt > tf){
      dualfprintf(fail_file,"Got t+dt>tf : %g %g %g\n",t,dt,tf);
      dt = tf - t;
      onemorestep=1;
    }
    else if (t + dt == tf){
      dualfprintf(fail_file,"Got t+dt==tf : %g %g %g\n",t,dt,tf);
      reallaststep=1;
    }
    // make sure don't get pathological case of dt=0 on last step
    if(dt<SMALL){
      dualfprintf(fail_file,"Got dt<SMALL : %g\n",dt);
      reallaststep=1;
      dt=SMALL;
    }
  }




  if(dt==0.0 && t>=tf){
    // somehow got here, then finish
    dualfprintf(fail_file,"SOMEHOW GOT HERE\n");
    reallaststep=1;
  }



  return(0);



}





/// things to do during substep before advance()
int pre_advance(int timeorder, int numtimeorders, int finalstep, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR])
{


  //////////////////
  //
  // prefixup
  //
  //////////////////
#if(EOMTYPE!=EOMFFDE && EOMTYPE!=EOMFFDE2)
  // force-free and cold GRMHD don't use pre_fixup, currently, even if could
  MYFUN(pre_fixup(STAGEM1, pi),"step_ch.c:step_ch_simple()", "pre_fixup()", 1);
#endif


  return(0);

}






/// things to do after advance()
int post_advance(int truestep, int *dumpingnext, int timeorder, int numtimeorders, int finalstep, SFTYPE boundtime, SFTYPE fluxtime, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{





  ////////////////
  //
  // Use tracked A_i to update magnetic field
  // Required to keep A_i in synch with field and only is different at the machine error level (which grows over long times, so why this is required and hardless)
  //
  ////////////////
  if(EVOLVEWITHVPOT && TRACKVPOT && (finalstep)){
    // if evolving with vpot, then assume good enough to use final full timestep A_i to obtain new point fields
    // less expensive than doing every substep and only wrong at machine error with factors extra for number of substeps
    // SUPERGODMARK: Check convergence rate and check errors!!  SUPERCHANGINGMARK

    //  only do this on the final step where A_i has been set as the cumulative Acum  (like ucum) so that not just an arbitrary intermediate step that redefines B's.

    // NOTEMARK: No, actually want to do on every substep that defines vpot (vpot holds substep (Uf-like) value on substeps except final step where it holds cumulative final value (Ucum-like)) since need staggered field (say inside NS for init.nsbh.c) to be consistent with A_i so that the smooth extrapolation of A_i translates into smooth Bstag^i so that interpolation of Bstag->Bcent gives smooth results.
    // NOTEMARK: If A_i is set inconsistent with EMF_i, then this step will reveal A_i overridding EMF_i results and probably producing undesirable results.
    // NOTEMARK: No again, if vpot update is consistent with EMFs (as it should be as code is setup in fluxvpot.c), then Bstag only off by machine error from A_i, so only need to do this infrequently.
    evolve_withvpot(boundtime, pf, pstag, ucons, vpot, Bhat, F1, F2, F3, Atemp,uconstemp); // boundtime since this time is used for a boundary call inside this function
  }




  //////////////////////
  //
  // need to compute various things for EOS
  // use primitive located where EOS will be used (pb)
  //
  //////////////////////
  compute_EOS_parms(WHICHEOS,GLOBALPOINT(EOSextraglobal),pb);


  /////////////////////////////////////
  //
  // post-advance bound and fixup
  //
  /////////////////////////////////////

#if( (1 == CHECKSOLUTION || UTOPRIMADJUST == UTOPRIMAVG) )
  // if CHECKSOLUTION==1, then need values to be bounded right now, since use them to check whether even good solutions are really good.
  // post_fixup() will use previous time step pff boundary values to fixup_utoprim() if this is not called.
  // bound advanced values before post_fixup() so fixup_utoprim() has updated boundary values to base fixup on.


  if(failed>0) dualfprintf(fail_file,"1failed=%d\n",failed);

  ////////////
  //
  // bounding 1 layer so fixups have that layer available for fixing up inversion failures
  // not done if don't care about perfect MPI consistency with failures -- then fixups don't use values across MPI boundaries.  Will be less robust in general, but still consistent behavior across MPI boundaries.
  //
  ////////////

#if(PRODUCTION==0)
  trifprintf("[b1");
#endif

  MYFUN(bound_beforeevolveprim(STAGEM1,finalstep,boundtime,pf,pstag,ucons,UTOPRIMFIXMPICONSISTENT>0),"step_ch.c:post_advance()", "bound_beforeevolveprim()", 1);

#if(PRODUCTION==0)
  trifprintf("]");
#endif




  if(failed>0) dualfprintf(fail_file,"2failed=%d\n",failed);


#if(PRODUCTION==0)
  trifprintf("[x");
#endif


  /////////////
  // done when all timeorders are completed, so stencil used doesn't matter
  //
  // postfixup
  // post_fixup: If this modifies values and then uses the modified values for other modified values, then must bound_prim() after this
  // if one doesn't care about MPI being same as non-MPI, then can move bound_prim() above to below error_check() and then remove prc<-pv in fixup_utoprim()
  //
  /////////////
  MYFUN(post_fixup(STAGEM1,finalstep,boundtime, pf,pb,ucons),"step_ch.c:advance()", "post_fixup()", 1);

#if(PRODUCTION==0)
  trifprintf( "]");
#endif




#else

  ////////////
  //
  // just report problems, don't fix them
  // Then no need for boundary calls
  //
  ////////////
  MYFUN(post_fixup_nofixup(STAGEM1,finalstep,boundtime, pf,pb,ucons),"step_ch.c:advance()", "post_fixup_nofixup()", 1);

#endif


  /////////////////////////////////////
  //
  // Synch CPUs
  //
  /////////////////////////////////////

  // must check before MPI operation (since asymmetries would
  // desynchronize cpus
  MYFUN(error_check(1),"step_ch.c", "error_check", 1);

  // bound final values (comes after post_fixup() since changes made by post_fixup)
  //#if(MPIEQUALNONMPI)



#if(ASYMDIAGCHECK)
  dualfprintf(fail_file,"1before bound\n");
  asym_compute_2(pf);
  dualfprintf(fail_file,"2before bound\n");
#endif




#if(PRODUCTION==0)
  trifprintf("[rf");
#endif



  if( DOGRIDSECTIONING){
    // this must come before last bound() call so boundary conditions set consistently with new section so next step has all values needed in ghost cells
    // redo all such enerregions since may be time-dependent
    recompute_fluxpositions(0,timeorder, numtimeorders, nstep,t);
  }

#if(PRODUCTION==0)
  trifprintf("]");
#endif



  if(timeorder==numtimeorders-1){
    // ensure total e-m conservation for final solution after all sub-steps
    int finaluu=1;
    consfixup_allzones(finaluu,pf,ucons);
  }




#if(PRODUCTION==0)
  trifprintf( "[b2");
#endif


  /////////////////////////////////////
  //
  // normal bondary call (only bounds p (centered) and pstag (staggered field) and not A_i or F or other things)
  // required in general
  //
  /////////////////////////////////////
  MYFUN(bound_evolveprim(STAGEM1,finalstep,boundtime,pf,pstag,ucons,USEMPI),"step_ch.c:post_advance()", "bound_evolveprim()", 2);

  //#endif

#if(PRODUCTION==0)
  trifprintf( "]");
#endif






  /////////////////////////////////////
  //
  // Compute divb for point-field method
  //
  /////////////////////////////////////
  if(truestep){
    // update Bhat so later can compute divb
    if(FLUXB==FLUXCTSTAG && extrazones4emf==0 && dofluxreconevolvepointfield==1 && DOENOFLUX == ENOFLUXRECON){
      //bound_uavg(STAGEM1,unew); // DEBUG
      // here utoinvert is point conserved field at staggered FACE position
      field_Bhat_fluxrecon(pf,ucons,Bhat);
    }
  }


  /////////////////////////////////////
  //
  // Compute current update
  //
  /////////////////////////////////////
  if(truestep){ // don't do if just passing through

    if(dumpingnext[0] || dumpingnext[1]){
#if(CALCFARADAYANDCURRENTS)

#if(PRODUCTION==0)
      trifprintf( "[cu");
#endif

      if((WHICHCURRENTCALC==CURRENTCALC0)||(WHICHCURRENTCALC==CURRENTCALC2)){
        // puts J at the time center, but hard to know if RK is at mid point in time except for midpoint method
        // compute current_doprecalc if near half-step in time
        if(
           ((numtimeorders>=3)&&(timeorder==1))
           ||((numtimeorders<=2)&&(timeorder==0))
           )
          current_doprecalc(CURTYPET,pf); // should be called using half-time step data
      }

#if(PRODUCTION==0)
      trifprintf( "]");
#endif

#endif
    }// end if dumping next




  }// end if truestep






  /////////////////////////////////////
  //
  // Diagnostics
  //
  /////////////////////////////////////
  if(truestep){ // don't do if just passing through
    /* perform diagnostics */
    // no error check since assume if step_ch passed, diag(1) will pass
    if (DODIAGS && DODIAGEVERYSUBSTEP ){ //SASMARK -- moved the diags calls here
      GLOBALPOINT(pdump) = pf;
      diag(DUMP_OUT,boundtime,nstep,realnstep); // boundtime corresponds to pf for this timeorder
#if(PRODUCTION==0)
      trifprintf( "D");
#endif
    }
  }

  return(0);
}






/// take full step using normal MPI method
/// This method transfers ALL boundary cells AT ONCE after updating ALL active cells
int step_ch_simplempi(int truestep, int *dumpingnext, FTYPE *fullndt, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  int advance(int truestep, int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
              FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
              FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
              FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
              FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
              FTYPE (*ui)[NSTORE2][NSTORE3][NPR],FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
              FTYPE *CUf,FTYPE *CUnew,SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime, int timeorder, int numtimeorders, FTYPE *ndt);
  
  int pre_advance(int timeorder, int numtimeorders, int finalstep, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR]);
  int asym_compute_2(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
  int pl,pliter;
  //
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pb)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pf)[NSTORE2][NSTORE3][NPR];
  FTYPE (*prevpf)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*pbb[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*pff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*uii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*uff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*ucum[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE CUf[MAXTIMEORDER][NUMDTCUFS],CUnew[MAXTIMEORDER][NUMDTCUFS];


  int i, j, k;
  int numtimeorders;
  int timeorder;
  int finalstep,initialstep;
  SFTYPE fluxdt[MAXTIMEORDER], boundtime[MAXTIMEORDER], fluxtime[MAXTIMEORDER];





  ////////////////////////
  //
  // setup time-stepping
  //
  ////////////////////////
  setup_rktimestep(truestep, &numtimeorders,prim,pstag,ucons,vpot,Bhat,GLOBALPOINT(pk),pii,pbb,pff,uii,uff,ucum,CUf,CUnew);


  /////////////////////////////////////
  //
  // Obtain initial time of substep, final time of substep, and true dt used for flux conservation that is used to iterate ucum in advance.c
  //
  /////////////////////////////////////
  get_truetime_fluxdt(numtimeorders, dt, CUf, CUnew, fluxdt, boundtime, fluxtime, NULL,NULL);


  // global debug tracking var
  numstepparts=numtimeorders;
  // Initialize old and new dt to be very large so minimization procedure always obtains new dt
  ndt=lastndt=BIG;




  /////////////////////////////
  //
  // general temporal loop
  //
  /////////////////////////////
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
    

    /////////////////////////////////////
    //
    //the starting and the ending times of the current substep
    //
    /////////////////////////////////////
    if(timeorder!=numtimeorders-1){
      tstepparti = t + CUf[timeorder][3] * dt;
      tsteppartf = t + CUf[timeorder][2] * dt +  CUf[timeorder][3] * dt;
    }
    else{
      tstepparti=t;
      tsteppartf=t+dt;
    }


    /////////////////////////////////////
    //
    // global debug tracking var
    //
    /////////////////////////////////////
    steppart=timeorder;
    if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;
    if(timeorder==0) initialstep=1; else initialstep=0;

#if(PRODUCTION==0)
    trifprintf("|%ds",timeorder);
#endif



    /////////////////////////////////////
    //
    // BEFORE ADVANCE
    //
    /////////////////////////////////////
    pre_advance(timeorder, numtimeorders, finalstep, pii[timeorder],pbb[timeorder],pff[timeorder]);
 


    /////////////////////////////////////
    //
    // ADVANCE
    //
    /////////////////////////////////////
    if(SPLITNPR){


      /////////////////////////////////////
      //
      // SPLITNPR ADVANCE
      //
      /////////////////////////////////////
    

#if(SPLITNPR&&(USESTOREDSPEEDSFORFLUX==0 || STOREWAVESPEEDS==0))
#error Not correct use of SPLITNPR and wavespeeds
#endif
      // Note that even though PLOOP is split-up, some things are still repeated for these 2 advance() calls:
      // 1) vchar for Riemann solver still needed -- could store vchar on first pass and use on second pass, but then using old B for vchar() on second pass and want to use new B everywhere for second pass
      // 2) dt is minimized twice -- seems second pass will overwrite first pass
      // 3) getting of geometry and state everywhere is duplicated
      //    NOTE OF IMPORTANCE: get_state needs all primitives, so must feed non-PLOOP quantities
      // 4) get_wavespeeds requires ALL quantities during interpolation -- must set USESTOREDSPEEDSFORFLUX==1 and STOREWAVESPEEDS>0 so computes wave speeds across 2 advance() calls and uses centered speed for interface to avoid needing during interpolation for p_l,p_r->wavespeed@interface
      //    Basically flux_compute_general(p_l,p_r needs velocities too) for first pass
      // 5) p2SFUevolve's get_state needs velocity -- that is, flux_B(B,v) so have to still interpolate v with B -- save v-interpolation (all interpolations) and avoid interpolations on second pass.
      //    Made pleft and pright [3]X in size to store across advance(), dq is already this way
      // 6) checked NPR and NPR2INTERP uses.  All such explicit uses are ok
      //    e.g. grep -e "[^\[]NPR" *.c

      // Things NOT done on first pass:
      // 1) source term if NOGDETBi=0
      // 2) metric update if evolving metric
      // 3) all non-B quantities if doing PLOOP or PLOOPINTERP
      // 4) diag_flux, diag_source, diag_source_all -- all use PDIAGLOOP
      // 5) fixup1zone() or fixup()

      // Things to fix:
      // 3) compute_df_line_new during interpolation needs get_state for u and b -- only for pressure indicator for WENO
      // 4) calculation of bsq requires b^\mu requires B and v/u/urel -- called PLOOP or what?
      // 5) Check all PLOOP's everywhere to ensure doing ALL NPR if required
      // 6) could bound only field when bounding between



      // choice for range of PLOOP
      // choose to only update magnetic field
      nprstart=0;
      nprend=2;
      nprlist[0]=B1;
      nprlist[1]=B2;
      nprlist[2]=B3;
      // choice for range of PLOOPINTERP
      // but interpolate everything on first pass (GODMARK: in reality only need to do field on timeorder=0 since last field update is next field pbb??? -- probably only true if doing simple RK methods (1,2))
#pragma omp parallel private(pl)
      {
        npr2interpstart=0;
        npr2interpend=NPR2INTERP-1;
        for(pl=npr2interpstart;pl<=npr2interpend;pl++)  npr2interplist[pl]=pl;
        // choice for range of PLOOPNOTINTERP
        npr2notinterpstart=0;
        npr2notinterpend=-1;
        npr2notinterplist[0]=0;
      }
      advancepassnumber=0;


      // advance (field only)
      // Only field parts of pff, uff, and ucum updated
      // on final timeorder, ucum used to get final B that will be used as B for final timeorder for non-field quantities
      MYFUN(advance(truestep, STAGEM1, pii[timeorder], pbb[timeorder], pff[timeorder], pstag, pl_ct, pr_ct, F1, F2, F3, vpot, uii[timeorder], uff[timeorder], ucum[timeorder], CUf[timeorder], CUnew[timeorder], fluxdt[timeorder], boundtime[timeorder], fluxtime[timeorder], timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_simplempi()", "advance()", 1);

      // only need to bound field, so control PLOOPMPI
      nprboundstart=0;
      nprboundend=2;
      nprboundlist[0]=B1;
      nprboundlist[1]=B2;
      nprboundlist[2]=B3;
      // with magnetic field, only  need to bound
      MYFUN(bound_evolveprim(STAGEM1,finalstep,boundtime[timeorder],pff[timeorder],pstag,uff[timeorder],USEMPI),"step_ch.c:step_ch_simplempi()", "bound_evolveprim()", 1);



      // now copy over pff to pbb so Bnew is used in flux calculation of non-field quantities
      COMPFULLLOOP{
        PLOOPBONLY(pl) MACP1A1(pbb,timeorder,i,j,k,pl)=MACP1A1(pff,timeorder,i,j,k,pl);
      }


      // choice for range of PLOOP
      // don't update field on second pass (result overwritten in utoprimgen.c even if done)
      nprstart=0;
      nprend=NPR-1-3; // no field
      for(pl=nprstart;pl<=nprend;pl++){
        if(pl>=B1) nprlist[pl]=pl+3; // skip field
        else nprlist[pl]=pl;
      }
      // choice for range of PLOOPINTERP
      // Interpolation of only fields on second pass since want to use updated field to compute flux and so at end of both steps first updated field and non-field quantites are fed to inversion for consistent P(Bnew) and Bnew is used
#pragma omp parallel private(pl)
      {
        npr2interpstart=0;
        npr2interpend=2;
        npr2interplist[0]=B1;
        npr2interplist[1]=B2;
        npr2interplist[2]=B3;
        // choice for range of PLOOPNOTINTERP
        // what not interpolating (all non-field):
        npr2notinterpstart=0;
        npr2notinterpend=NPR-1-3; // no field
        for(pl=npr2notinterpstart;pl<=npr2notinterpend;pl++){
          if(pl>=B1) npr2notinterplist[pl]=pl+3; // skip field
          else npr2notinterplist[pl]=pl;
        }
      }
      advancepassnumber=1;


      // advance (non-field quantities)
      // only non-field parts of pff, uff, ucum updated
      MYFUN(advance(truestep, STAGEM1, pii[timeorder], pbb[timeorder], pff[timeorder], pstag, pl_ct, pr_ct, F1, F2, F3, vpot, uii[timeorder], uff[timeorder], ucum[timeorder], CUf[timeorder], CUnew[timeorder], fluxdt[timeorder], boundtime[timeorder], fluxtime[timeorder], timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_simplempi()", "advance()", 1);



      //////////////////////////
      // go back to defaults

      // choice for range of PLOOP
      nprstart=0;
      nprend=NPR-1;
      for(pl=nprstart;pl<=nprend;pl++) nprlist[pl]=pl;
#pragma omp parallel private(pl)
      {
        // choice for range of PLOOPINTERP
        npr2interpstart=0;
        npr2interpend=NPR2INTERP-1;
        for(pl=npr2interpstart;pl<=npr2interpend;pl++) npr2interplist[pl]=pl;
        // choice for range of PLOOPNOTINTERP
        npr2notinterpstart=0;
        npr2notinterpend=-1;
        npr2notinterplist[0]=0;
      }
      advancepassnumber=-1;

      // default choice for range of PLOOPMPI
      nprboundstart=0;
      nprboundend=NPRBOUND-1;
      for(pl=nprboundstart;pl<=nprboundend;pl++) nprboundlist[pl]=pl;




    }
    else{


      /////////////////////////////////////
      //
      // NORMAL (not SPLITNPR) ADVANCE
      //
      /////////////////////////////////////

      advancepassnumber=-1; // implies do all things, no split
      // advance (all)
      MYFUN(advance(truestep, STAGEM1, pii[timeorder], pbb[timeorder], pff[timeorder], pstag, pl_ct, pr_ct, F1, F2, F3, vpot, uii[timeorder], uff[timeorder], ucum[timeorder], CUf[timeorder], CUnew[timeorder], fluxdt[timeorder], boundtime[timeorder], fluxtime[timeorder], timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_simplempi()", "advance()", 1);



    } // end else if normal advance()







    ///////////////
    //    
    // need to minimize dt over all substeps rather than just last step
    // lastndt holds cumulative-all-substep ndt
    //
    ///////////////
    if(ndt<lastndt) lastndt=ndt;



    /////////////////////////////////////
    //
    // POST ADVANCE
    //
    /////////////////////////////////////
    post_advance(truestep, dumpingnext, timeorder, numtimeorders, finalstep, boundtime[timeorder], fluxtime[timeorder], pii[timeorder],pbb[timeorder],pff[timeorder],pstag,ucons,vpot,Bhat, F1, F2, F3, Atemp, uconstemp);



  }// end timeorder


  
  /////////////////////////
  //
  // pass back the final minimal dt over all substeps
  //
  /////////////////////////
  *fullndt = lastndt;



  /* done! */
  return (0);
}



/// Obtain initial time of substep, final time of substep, and true dt used for flux conservation that is used to iterate ucum in advance.c
void get_truetime_fluxdt(int numtimeorders, SFTYPE localdt, FTYPE (*CUf)[NUMDTCUFS], FTYPE (*CUnew)[NUMDTCUFS], SFTYPE *fluxdt, SFTYPE *boundtime, SFTYPE *fluxtime, SFTYPE *tstepparti, SFTYPE *tsteppartf)
{
  int timeorder;
  SFTYPE ufdt[MAXTIMEORDER],ucumdt[MAXTIMEORDER];
  SFTYPE oldufdt,olducumdt;
  FTYPE Ui, dUriemann, dUgeom ;


  //  NOT YET:
  /////////////////////////////////////
  //
  //the starting and the ending times of the current substep
  //
  /////////////////////////////////////
  //  if(timeorder!=numtimeorders-1){
  //    *tstepparti = t + CUf[timeorder][3] * localdt;
  //    *tsteppartf = t + CUf[timeorder][2] * localdt +  CUf[timeorder][3] * localdt;
  //  }
  //  else{
  //    *tstepparti=t;
  //    *tsteppartf=t+localdt;
  //  }





  // initialize
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
    fluxdt[timeorder] = 0.0;
    boundtime[timeorder] = 0.0;
    fluxtime[timeorder] = 0.0;
  }
  //////////////////////
  //
  // Get fluxdt (for up to explicit RK5)
  // See setup_rktimestep() for how Uf^i and unew^i are defined, such that fluxdt[] below just adds all dUe terms appropriately.
  //
  //////////////////////
  for(timeorder=4;timeorder<numtimeorders;timeorder++){ // for >=RK5
    fluxdt[timeorder-4] += CUnew[timeorder][2]*CUf[timeorder][1]*CUf[timeorder-1][1]*CUf[timeorder-2][1]*CUf[timeorder-3][1]*CUf[timeorder-4][2]*localdt;
  }

  for(timeorder=3;timeorder<numtimeorders;timeorder++){ // for >=RK4
    fluxdt[timeorder-3] += CUnew[timeorder][2]*CUf[timeorder][1]*CUf[timeorder-1][1]*CUf[timeorder-2][1]*CUf[timeorder-3][2]*localdt;
  }

  for(timeorder=2;timeorder<numtimeorders;timeorder++){ // for >=RK3
    fluxdt[timeorder-2] += CUnew[timeorder][2]*CUf[timeorder][1]*CUf[timeorder-1][1]*CUf[timeorder-2][2]*localdt;
  }

  for(timeorder=1;timeorder<numtimeorders;timeorder++){ // for >=RK2
    fluxdt[timeorder-1] += CUnew[timeorder][2]*CUf[timeorder][1]*CUf[timeorder-1][2]*localdt;
  }

  for(timeorder=0;timeorder<numtimeorders;timeorder++){ // for >=RK1
    fluxdt[timeorder] += (CUnew[timeorder][1] + CUnew[timeorder][2]*CUf[timeorder][2])*localdt;
  }

  if(nstep==0){
    for(timeorder=0;timeorder<numtimeorders;timeorder++){
      dualfprintf(fail_file,"timeorder=%d fluxdt/dt=%g\n",timeorder,fluxdt[timeorder]/localdt);
    }
  }


  // assuming fluxdt used before any calls in advance.c
  // then represents *amount* of flux'ed stuff in time dt
  // This is NOT the time of the flux evaluation
  // This is to be used with diag_flux and other things multiplied by dt for diagnostics only
  // Other things (e.g. sources()) use same dt as dUriemann automatically, so fluxdt is not for anything but diagnostics
  // Only exception is updating vpot that is out of place and uses ucum-type updates


  // loop up to current substep to get current fluxdt
  Ui=dUgeom=0.0; // don't care about update from non-flux terms
  FTYPE uf;
  FTYPE dUnongeomall[MAXTIMEORDER]={0.0};
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
    if(timeorder==0){
      oldufdt=0.0;
    }
    else{
      oldufdt=ufdt[timeorder-1];
    }

    // follows dUtoU() in advance.c
    ufdt[timeorder] = UFSET(CUf[timeorder],localdt,Ui,oldufdt,dUriemann,dUgeom,dUnongeomall);
    // below is NOT += since just want current change, not all prior changes
    // if did +=, then get 1 for timeorder==numtimeorders-1 as required!
    ucumdt[timeorder] = UCUMUPDATE(CUnew[timeorder],localdt,Ui,ufdt[timeorder],dUriemann,dUgeom,dUnongeomall);

    //  time of pf at end of substep
    boundtime[timeorder] = t + ufdt[timeorder];

    // time of pb when used to compute flux
    // and time of pb corresponds to time of previous pf (except for first substep where it's just t)
    //     tstepparti = t + CUf[timeorder][3] * dt;
    //      tsteppartf = t + CUf[timeorder][2] * dt +  CUf[timeorder][3] * dt;
    if(timeorder>0) fluxtime[timeorder] =  t + CUf[timeorder-1][2] * dt +  CUf[timeorder-1][3] * dt;
    else fluxtime[timeorder] =  t;

  }


  

#if(0)
  //////////////////////
  //
  // Get boundtime
  //
  // Get current time for pf
  // The below works because pi=p always and so always just use 
  // assuming bound called after advance(), then need to get time of next flux evaluation
  // This corresponds to time where pb *will be* located
  // 
  //
  //////////////////////
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
    if(timeorder<numtimeorders-1){
      boundtime[timeorder] = t + localdt*CUf[timeorder][3]+ localdt*CUf[timeorder][2];
    }
    else{
      // if timeorder==numtimeorders-1, then final step and bound should be set for time of flux of next step, which is full t+dt
      boundtime[timeorder] = t + localdt;
    }
  }
#endif


#if(0)
  // DEBUG:
  dualfprintf(fail_file,"FLUXDT/BOUNDTIME: nstep=%ld ",nstep);
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
    dualfprintf(fail_file,"to=%d fluxdt=%21.15g fluxdtperdt=%21.15g boundtime=%21.15g fluxtime=%21.15g\n",timeorder,fluxdt[timeorder],fluxdt[timeorder]/dt,boundtime[timeorder],fluxtime[timeorder]);
  }
  dualfprintf(fail_file,"\n");
#endif


}








/// take full step using super MPI method
///
/// Method 1) 
/// This method transfers SOME boundary cells BETWEEN EACH update of active cells
/// Simplest version:
/// 1) Iterate a layer of size NBND? (Dependent Active = DA) cells that are just inside of ghost cells but are still active cells.  These are needed by other CPUs to do next iteration
/// 2) Start transfer of DA cells to other CPUs
/// 3) Iterate region never used by other CPUs identified as (Super Active = SA) cells
/// 4) See if DA cells finished transfer
/// 5) Once finished, go back to #1 above and repeat.
///
/// Method 2) 
/// Better version (noted by Sasha and noted as synched with MPI by Jon):  Idea: Keep MPI updates in each dimension and direction in synch with computations
///  1) See if received THIS core's X1DN data from other X1DN CPU
///  2) Iterate X1BND (NBND1 X full x2 X full x3) of DA cells on X1DN side
///  3) Send #2 to X1DN CPU using single X1DN MPI call (already split off like this but just looped over)
///  4) See if received THIS core's X1UP data from other X1UP CPU
///  5) Iterate X1BND (NBND1 X full x2 X full x3) of DA cells on X1UP side
///  6) Send #5 to X1UP CPU using single X1UP MPI call ("")
///  7) See if received THIS core's X2DN data from other X2DN CPU
///  8) Iterate X2BND (NBND2 X active x1 X full x3) of DA cells on X2DN side
///  9) Send #8 to X2DN CPU using single X2DN MPI call ("") -- which only transfers NBND2 X active x1 X full x3 as exactly required
/// 10) See if received THIS core's X2UP data from other X2UP CPU
/// 11) Iterate X2BND (NBND2 X active x1 X full x3) of DA cells on X2UP side
/// 12) Send #11 to X2UP CPU using single X2UP MPI call ("") -- ""
/// 13) See if received THIS core's X3DN data from other X3DN CPU
/// 14) Iterate X3BND (NBND3 X active x1 X active x2) of DA cells on X3DN side
/// 15) Send #14 to X3DN CPU using single X3DN MPI call ("") -- which only transfers NBND3 X active x1 X active x2 as exactly required
/// 16) See if received THIS core's X3UP data from other X3UP CPU
/// 17) Iterate X3BND (NBND3 X active x1 X active x2) of DA cells on X3UP side
/// 18) Send #17 to X3UP CPU using single X3UP MPI call ("") -- ""
/// 19) Iterate SA cells -- no need to check since SA cells are independent of MPI computations.
/// 20) Repeat, but reverse UP and DN labels.  This way correctly first waits for first things should have recieved.  May be minor difference.
int step_ch_supermpi(int truestep, int *dumpingnext, FTYPE *fullndt, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  int advance(int truestep, int stage, FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*pf)[NSTORE2][NSTORE3][NPR],
              FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
              FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
              FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],
              FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
              FTYPE (*ui)[NSTORE2][NSTORE3][NPR],FTYPE (*uf)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR],
              FTYPE *CUf,FTYPE *CUnew,SFTYPE fluxdt, SFTYPE boundtime, SFTYPE fluxtime, int timeorder, int numtimeorders, FTYPE *ndt);
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  int timeorder;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pb)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pf)[NSTORE2][NSTORE3][NPR];
  FTYPE (*prevpf)[NSTORE2][NSTORE3][NPR];
  FTYPE (*pii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*pbb[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*pff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*uii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*uff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  FTYPE (*ucum[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE CUf[MAXTIMEORDER][NUMDTCUFS],CUnew[MAXTIMEORDER][NUMDTCUFS];
  int i, j, k, pl, pliter;
  int numtimeorders;
  int finalstep;
  //  extern int horizon_flux(FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], SFTYPE Dt);
  SFTYPE fluxdt[MAXTIMEORDER], boundtime[MAXTIMEORDER], fluxtime[MAXTIMEORDER];





  ////////////////////////
  //
  // setup time-stepping
  //
  ////////////////////////
  setup_rktimestep(truestep, &numtimeorders,prim,pstag,ucons,vpot,Bhat,GLOBALPOINT(pk),pii,pbb,pff,uii,uff,ucum,CUf,CUnew);


  /////////////////////////////////////
  //
  // Obtain initial time of substep, final time of substep, and true dt used for flux conservation that is used to iterate ucum in advance.c
  //
  /////////////////////////////////////
  get_truetime_fluxdt(numtimeorders, dt, CUf, CUnew, fluxdt, boundtime, fluxtime, NULL, NULL);



  // SPECIAL BOUNDARY/COMPUTATION MPI METHOD (out of date, and doesn't yet work right even if essentially complete code)
  /* check timestep */
  //  if (dt < MINDT) {
  //  trifprintf( "timestep too small\n");
  //  myexit(11);
  // }
  
  lastndt=1E30; // initialize lastndt
  for(timeorder=1;timeorder<=numtimeorders;timeorder++){


    /////////////////////////////////////
    //
    //the starting and the ending times of the current substep
    //
    /////////////////////////////////////
    if(timeorder!=numtimeorders-1){
      tstepparti = t + CUf[timeorder][3] * dt;
      tsteppartf = t + CUf[timeorder][2] * dt +  CUf[timeorder][3] * dt;
    }
    else{
      tstepparti=t;
      tsteppartf=t+dt;
    }

    if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;



#if(PRODUCTION==0)
    trifprintf("-to%d/%d-",timeorder,numtimeorders);
#endif
    if(numtimeorders==2){
      // note that pb (used for flux calc which has a stencil
      // calculation) must be different from pf so new stencils in
      // different stages won't affect stencil calculations -- must
      // use old values, not new from most previous temporary stage
      //
      // pi however can be the same as pf since each pi is replaced 1
      // zone at a time with a 0 stencil.
      if(timeorder==1){
        pi=prim;
        pb=prim;
        pf=GLOBALPOINT(pk)[0]; // different already, so good for simulbccalc
        prevpf=prim; // previous final true array
        mydt=0.5*dt;
      }
      else if(timeorder==2){
        pi=prim;
        pb=GLOBALPOINT(pk)[0];
        pf=prim;
        prevpf=GLOBALPOINT(pk)[0];
        mydt=dt;
      }
    }
    else if(numtimeorders==1){
      pi=prim;
      pb=prim;
      if(SIMULBCCALC<=0) pf=prim; else pf=GLOBALPOINT(pk)[0]; // need to be different if doing simulbccalc
      prevpf=prim;
      mydt=dt;
    }
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}

    
    // initialize bound stage
    if(SIMULBCCALC>=1) boundstage=STAGE0;
    else boundstage=STAGEM1;
    for(stage=stagei;stage<=stagef;stage++){
#if(PRODUCTION==0)
      if(SIMULBCCALC>=1) trifprintf("!s%d!",stage);
#endif
      // setup stage loop
#if(SIMULBCCALC==2)
#if(TYPE2==1)
      // GODMARK: isf1, etc. are NOT defined?!
      STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
      STAGECONDITION(0,N1,-1,N2,isf1,ief1,jsf1,jef1);
      STAGECONDITION(-1,N1,0,N2,isf2,ief2,jsf2,jef2);
      STAGECONDITION(0,N1,0,N2,ise,iee,jse,jee);
      STAGECONDITION(0,N1,0,N2-1,isf1ct,ief1ct,jsf1ct,jef1ct);
      STAGECONDITION(0,N1-1,0,N2,isf2ct,ief2ct,jsf2ct,jef2ct);
      STAGECONDITION(-1,N1,-1,N2,isdq,iedq,jsdq,jedq);
      STAGECONDITION(-2,N1+1,-2,N2+1,ispdq,iepdq,jspdq,jepdq);
      // GODMARK : probably not right for general boundary condition size
#endif
#endif

      // only bounding if safe zones, unsafe needs bz complete
      if(stage<STAGE2){
        bound_beforeevolveprim(boundstage, finalstep, boundtime[timeorder], prevpf,pstag,ucons,USEMPI); // pstag,ucons not right for supermpi
        if(stage!=STAGEM1) boundstage++;
      }

      // done here instead of local since pseudo-complicated calculation that might slow the dq calculation if done locally per zone
      MYFUN(pre_fixup(stage, prevpf),"step_ch.c:advance()", "pre_fixup()", 1);

      // go from previous solution to new solution
      partialstep=timeorder;      
      // not right for numtimeorders==4 // GODMARK
      // advance
      MYFUN(advance(truestep,-1, pii[timeorder], pbb[timeorder], pff[timeorder], pstag, pl_ct, pr_ct, F1, F2, F3, vpot, uii[timeorder], uff[timeorder], ucum[timeorder],CUf[timeorder], CUnew[timeorder], fluxdt[timeorder], boundtime[timeorder], fluxtime[timeorder], timeorder,numtimeorders,&ndt),"step_ch.c:step_ch_supermpi()", "advance()", 1);
      // must check before MPI operation (since asymmetries would desynchronize cpus)
      if(stage<STAGE2){
        MYFUN(error_check(1),"step_ch.c", "error_check", 1);
      }
      if(stage!=STAGEM1){
        if(stage<STAGE2){
          bound_evolveprim(boundstage,finalstep, boundtime[timeorder], prevpf,pstag,ucons,USEMPI);
          boundstage++;
        }
      }
      if(timeorder==numtimeorders){
        if(ndt>lastndt) ndt=lastndt; // don't change if last was lower
        else lastndt=ndt; // new is lower, keep it
      }
    }
    if(timeorder==numtimeorders){// only do on full step
      // find global minimum value of ndt over all cpus
      mpifmin(&ndt);
    }
    // done when all stages are completed, so stencil used doesn't matter
    
    MYFUN(post_fixup(STAGEM1,finalstep, boundtime[timeorder], pf,pb,ucons),"step_ch.c:advance()", "post_fixup()", 1);
  }


  // pass back the final minimal dt over all substeps
  *fullndt = lastndt;


  // copy the contents to the final working array
  if((numtimeorders==1)&&(SIMULBCCALC>=1)) COMPFULLLOOP PLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(pf,i,j,k,pl);
  

  /* done! */
  return (0);
}








/// Setup RK steps
/// for the ith stage:
///
/// Uf^i = ulast^i = CUf^{i0} Ui^i + CUf^{i1} ulast^i + CUf^{i2} dUexplicit^i + CUf^{i4} dUimplicit^i
///
/// unew^i = CUnew^{i0} Ui^i + CUnew^{i1} dUexplicit^i + CUnew^{i2} Uf^i  + CUnew^{i4} dUimplicit^i
///
/// [SUPERNOTE: currently only TIMEORDER==4 has final unew different from final uf.  This factis usd in utoprimgen.c to avoid inversion if requested, unless cannot because unew must itself be inverted.]
///
/// (how also used) CUf[timeorder][2] : t + dt*CUf[timeorder][3]+ dt*CUf[timeorder][2] = time of pf
///
/// CUf[timeorder][3] : t + dt*CUf[timeorder][3] = time of pi
///
/// CUnew[timeorder][3] : t + dt*CUnew[timeorder][3] = approximate time of pb
///
/// GODMARK: Currently RK methods all feed pf into pb on next step, so only need one pstag.  In general should have multiple pstag's.
void setup_rktimestep(int truestep, int *numtimeorders,
                      FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR],
                      FTYPE (*pk)[NSTORE1][NSTORE2][NSTORE3][NPR],
                      FTYPE (*pii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*pbb[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*pff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],
                      FTYPE (*uii[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*uff[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],FTYPE (*ucum[MAXTIMEORDER])[NSTORE2][NSTORE3][NPR],
                      FTYPE (*CUf)[NUMDTCUFS],FTYPE (*CUnew)[NUMDTCUFS])
{


  // initialize CUf and CUnew to be zero
  int ii,jj;
  for(ii=0;ii<MAXTIMEORDER;ii++){
    for(jj=0;jj<NUMDTCUFS;jj++){
      CUf[ii][jj]=CUnew[ii][jj]=0.0;
    }
  }


  
  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // EXPLICIT TIME STEPPING
  //
  //////////////////////////////////////////////////////////////////////////////////////////


  if(TIMETYPE==TIMEEXPLICIT){



    // to avoid special copying of final pff->p, always use p as final pff
    if(TIMEORDER==4 && truestep){
      // RK4 stepping
      *numtimeorders=4;

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0;  CUf[0][1]=0.0;      CUf[0][2]=0.5;  CUf[0][3] = 0.0;  CUf[0][4+0] = CUf[0][2];
      CUf[1][0]=1.0;  CUf[1][1]=0.0;      CUf[1][2]=0.5;  CUf[1][3] = 0.0;  CUf[1][4+1] = CUf[1][2];
      CUf[2][0]=1.0;  CUf[2][1]=0.0;      CUf[2][2]=1.0;  CUf[2][3] = 0.0;  CUf[2][4+2] = CUf[2][2];
      CUf[3][0]=1.0;  CUf[3][1]=0.0;      CUf[3][2]=1.0;  CUf[3][3] = 0.0;  CUf[3][4+3] = CUf[3][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      CUnew[0][0]=1.0;  CUnew[0][1]=1.0/6.0;      CUnew[0][2]=0.0;  CUnew[0][3] = 0.0;  CUnew[0][4+0] = CUnew[0][1];
      CUnew[1][0]=0.0;  CUnew[1][1]=1.0/3.0;      CUnew[1][2]=0.0;  CUnew[1][3] = 0.5;  CUnew[1][4+1] = CUnew[1][1];
      CUnew[2][0]=0.0;  CUnew[2][1]=1.0/3.0;      CUnew[2][2]=0.0;  CUnew[2][3] = 0.5;  CUnew[2][4+2] = CUnew[2][1];
      CUnew[3][0]=0.0;  CUnew[3][1]=1.0/6.0;      CUnew[3][2]=0.0;  CUnew[3][3] = 1.0;  CUnew[3][4+3] = CUnew[3][1];

      //primitive values used for initial state, fluxes, final state (where you output)
      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0]; // produces U1
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
      pii[2]=p;    pbb[2]=pk[1];   pff[2]=pk[0]; // produces U3
      pii[3]=p;    pbb[3]=pk[0];   pff[3]=p; // produces U4 (only dUe part used)
    
      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
      uii[2]=GLOBALPOINT(uinitialglobal);  uff[2]=GLOBALPOINT(ulastglobal); ucum[2]=ucons;
      uii[3]=GLOBALPOINT(uinitialglobal);  uff[3]=GLOBALPOINT(ulastglobal); ucum[3]=ucons;

      // GODMARK: note that pbstag (staggered field from conserved de-averaging and inversion to primitive no geometry version) is always same memory space and comes from operating on final inverted quantity (ulastglobal or ucons), so just use same quantity for now and avoid adding extra code for pbstag[]
      //    pbstag[0]=pstagglobal;
      //    pbstag[1]=pstagglobal;
      //    pbstag[2]=pstagglobal;
      //    pbstag[3]=pstagglobal;

    }
    else if(TIMEORDER==3  && truestep){
      // TVD optimal RK3 method as in Shu's report
      *numtimeorders=3;
    
      // Ui ulastglobal dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0;      CUf[0][1]=0.0;      CUf[0][2]=1.0;      CUf[0][3] = 0.0;  CUf[0][4+0] = CUf[0][2];
      CUf[1][0]=3.0/4.0;  CUf[1][1]=1.0/4.0;  CUf[1][2]=1.0/4.0;  CUf[1][3] = 0.0;  CUf[1][4+1] = CUf[1][2];
      CUf[2][0]=1.0/3.0;  CUf[2][1]=2.0/3.0;  CUf[2][2]=2.0/3.0;  CUf[2][3] = 0.0;  CUf[2][4+2] = CUf[2][2];
    
      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U3
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=0.0;  CUnew[0][3] = 0.0;     CUnew[0][4+0] = CUnew[0][1];
      CUnew[1][0]=0.0;   CUnew[1][1]=0.0;      CUnew[1][2]=0.0;  CUnew[1][3] = 1.0;     CUnew[1][4+1] = CUnew[1][1];
      CUnew[2][0]=0.0;   CUnew[2][1]=0.0;      CUnew[2][2]=1.0;  CUnew[2][3] = 1.0/4.0; CUnew[2][4+2] = CUnew[2][1];
    
      //always starting the substeps from the initial time
      pii[0]=p;      pbb[0]=p;       pff[0]=pk[0]; // produces U1
      pii[1]=p;      pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
      pii[2]=p;      pbb[2]=pk[1];   pff[2]=p; // produces U3

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
      uii[2]=GLOBALPOINT(uinitialglobal);  uff[2]=GLOBALPOINT(ulastglobal); ucum[2]=ucons;

    }
    else if(TIMEORDER==2  && truestep){
#if(1)
      // midpoint method

      *numtimeorders=2;

      // old ucons not used for this method (i.e. [?][1]=0)
      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=0.5; CUf[0][3] = 0.0;  CUf[0][4+0] = CUf[0][2];
      CUf[1][0]=1.0; CUf[1][1]=0.0; CUf[1][2]=1.0; CUf[1][3] = 0.0;  CUf[1][4+1] = CUf[1][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U2
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=0.0;  CUnew[0][3] = 0.0;  CUnew[0][4+0] = CUnew[0][1];
      CUnew[1][0]=0.0;   CUnew[1][1]=0.0;      CUnew[1][2]=1.0;  CUnew[1][3] = 0.5;  CUnew[1][4+1] = CUnew[1][1];

      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;

#else
      *numtimeorders=2;
      // TVD RK2 (Chi-Wang Shu 1997 - eq 4.10)
      // actually less robust than generic midpoint method above

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0; CUf[0][3] = 0.0;  CUf[0][4+0] = CUf[0][2];
      CUf[1][0]=0.5; CUf[1][1]=0.5; CUf[1][2]=0.5; CUf[1][3] = 0.0;  CUf[1][4+1] = CUf[1][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U2
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=0.0;   CUnew[0][3] = 0.0;  CUnew[0][4+0] = CUnew[0][1];
      CUnew[1][0]=0.0;   CUnew[1][1]=0.0;      CUnew[1][2]=1.0;   CUnew[1][3] = 1.0;  CUnew[1][4+1] = CUnew[1][1];

      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;

#endif
    }
    else if(TIMEORDER==1 ||  dt==0.0){ // dt==0.0 case is case when just passing through
      // Euler method
      *numtimeorders=1;

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0; CUf[0][3] = 0.0;   CUf[0][4] = CUf[0][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U1
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=1.0;   CUnew[1][3] = 0.0;     CUnew[1][4] = CUnew[0][1];

      pii[0]=p;    pbb[0]=p;    pff[0]=p;

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;

    }
  }



  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // IMPLICIT TIME STEPPING
  //
  // First stage is just implicit solution (no flux needed)
  // Last stage is just flux (no implicit source term needed (i.e. no M[timeorder-1] needed), but geometry or any other explicits still needed)
  //
  // to avoid special copying of final pff->p, always use p as final pff
  //
  //////////////////////////////////////////////////////////////////////////////////////////
  if(TIMETYPE==TIMEIMPLICIT){

    if(TIMEORDER==5 && truestep){
      // IMEX3
      *numtimeorders=5;

      static FTYPE imp_alpha=0.24169426078821;
      static FTYPE imp_beta=0.06042356519705;
      static FTYPE imp_eta=0.12915286960590;

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0;     CUf[0][1]=0.0;      CUf[0][2]=0.0;     CUf[0][3] = 0.0;  CUf[0][4+0] = imp_alpha; // M0 no flux step : U0 = [Un] + \alpha dt M0
      CUf[1][0]=2.0;     CUf[1][1]=-1.0;     CUf[1][2]=0.0;     CUf[1][3] = 0.0;  CUf[1][4+1] = imp_alpha; // M1 no flux step : U1 = [2Un - U0] + \alpha dt M1
      CUf[2][0]=1.0;     CUf[2][1]=0.0;      CUf[2][2]=1.0;     CUf[2][3] = 0.0;  CUf[2][4+1] = (1.0-imp_alpha);  CUf[2][4+2] = imp_alpha; // M2 and F1 : U2 = [Un + dt F1 + (1-\alpha) dt M1] + \alpha dt M2
      CUf[3][0]=3.0/4.0; CUf[3][1]=1.0/4.0;  CUf[3][2]=1.0/4.0; CUf[3][3] = 0.0;  CUf[3][4+0] = imp_beta; CUf[3][4+1] = (-1.0+imp_alpha+4.0*imp_eta)/4.0;  CUf[3][4+2] = (2.0-5.0*imp_alpha-4.0*(imp_beta+imp_eta))/4.0; CUf[3][4+3] = imp_alpha; // M3 and F2 : U3 = [(3/4)Un + (1/4)U2 + (dt/4)F2 + \beta dt M0 + (1/4)(-1+\alpha+4\eta) dt M1 + (1/4)(2-5\alpha-4(\beta+\eta)) dt M2] + \alpha dt M3
      CUf[4][0]=1.0/3.0; CUf[4][1]=2.0/3.0;  CUf[4][2]=2.0/3.0; CUf[4][3] = 0.0;  CUf[4][4+0] = (-2.0*imp_beta/3.0); CUf[4][4+1] = (1.0-4.0*imp_eta)/6.0; CUf[4][4+2] = (-1.0+4.0*imp_alpha+4.0*(imp_beta+imp_eta))/6.0; CUf[4][4+3] = 4.0*(1.0-imp_alpha)/6.0; // F3 : U4 = [(1/3)Un + (2/3)U3 + (2dt/3) F3 + (-2beta dt/3) M0 + (1/6)(1-4\eta) dt M1 + (1/6)(-1 + 4\alpha + 4(\beta+\eta)) dt M2 + (1/6)(4(1-\alpha))dt M3]

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // NOTE: Could instead only use final U^{n+1} = U4 as possible for this case and setup above for final step, but current way is more local and really would require less memory (but not implemented to do that)
      CUnew[0][0]=1.0/3.0; CUnew[0][1]=0.0;          CUnew[0][2]=0.0;     CUnew[0][3] = 0.0;  CUnew[0][4+0] = (-2.0*imp_beta/3.0); // + (1/3)Un + (-2\beta/3) dt M0
      CUnew[1][0]=0.0;     CUnew[1][1]=0.0;          CUnew[1][2]=0.0;     CUnew[1][3] = 0.0;  CUnew[1][4+1] = (1.0-4.0*imp_eta)/6.0; // + ((1-4\eta)/6) dt M1
      CUnew[2][0]=0.0;     CUnew[2][1]=0.0;          CUnew[2][2]=0.0;     CUnew[2][3] = 0.0;  CUnew[2][4+2] = (-1.0+4.0*imp_alpha+4.0*(imp_beta+imp_eta))/6.0; // + ( (-1 + 4\alpha + 4(\beta+\eta))/6) dt M2
      CUnew[3][0]=0.0;     CUnew[3][1]=0.0;          CUnew[3][2]=2.0/3.0; CUnew[3][3] = 0.0;  CUnew[3][4+3] = 4.0*(1.0-imp_alpha)/6.0; // + (2/3)U3 + (4(1-\alpha)/6) dt M3
      CUnew[4][0]=0.0;     CUnew[4][1]=2.0/3.0;      CUnew[4][2]=0.0;     CUnew[4][3] = 0.0;  CUnew[4][4+4] = 0.0; // + (2/3)dt F3

      //primitive values used for initial state, fluxes, final state (where you output)
      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0]; // produces U1
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
      pii[2]=p;    pbb[2]=pk[1];   pff[2]=pk[0]; // produces U3
      pii[3]=p;    pbb[3]=pk[0];   pff[3]=pk[1]; // produces U4 (only dUe part used)
      pii[4]=p;    pbb[4]=pk[1];   pff[4]=p; // produces U5 (only dUe part used)
    
      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
      uii[2]=GLOBALPOINT(uinitialglobal);  uff[2]=GLOBALPOINT(ulastglobal); ucum[2]=ucons;
      uii[3]=GLOBALPOINT(uinitialglobal);  uff[3]=GLOBALPOINT(ulastglobal); ucum[3]=ucons;
      uii[4]=GLOBALPOINT(uinitialglobal);  uff[4]=GLOBALPOINT(ulastglobal); ucum[4]=ucons;

    }
    else if(TIMEORDER==4 && truestep){
      // IMEX2B
      *numtimeorders=4;

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0;     CUf[0][1]=0.0;      CUf[0][2]=0.0;     CUf[0][3] = 0.0;  CUf[0][4+0] = 1.0/4.0;   // M0 no flux step : U0 = [Un] + (dt/4)M0
      CUf[1][0]=1.0;     CUf[1][1]=0.0;      CUf[1][2]=0.5;     CUf[1][3] = 0.0;  CUf[1][4+1] = 1.0/4.0;   // F0 and M1 : U1 = [Un + (dt/2)F0] + (dt/4)M1
      CUf[2][0]=0.0;     CUf[2][1]=1.0;      CUf[2][2]=1.0/2.0; CUf[2][3] = 0.0;  CUf[2][4+0] = 1.0/3.0; CUf[2][4+1] = 1.0/12.0; CUf[2][4+2] = 1.0/3.0; // F1 and M2 : U2 = [U1 + (dt/2)F1 + (dt/3)M0 + (dt/12)M1] + (dt/3)M2
      CUf[3][0]=1.0/3.0; CUf[3][1]=2.0/3.0;  CUf[3][2]=1.0/3.0; CUf[3][3] = 0.0;  CUf[3][4+0] = 1.0/9.0; CUf[3][4+1] = 1.0/9.0;  CUf[3][4+2] = 1.0/9.0; // F2 no new implicit step // U3 = [(1/3)Un + (2/3)U2 + (dt/3)F2 + (dt/9)(M0 + M1 + M2)]

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // NOTE: Could instead only use final U^{n+1} = U3 as possible for this case and setup above for final step, but current way is more local and really would require less memory (but not implemented to do that)
      CUnew[0][0]=1.0/3.0; CUnew[0][1]=0.0;          CUnew[0][2]=0.0;     CUnew[0][3] = 0.0;  CUnew[0][4+0] = (1.0/9.0);   // + (1/3)Un + (1/9)M0
      CUnew[1][0]=0.0;     CUnew[1][1]=0.0;          CUnew[1][2]=0.0;     CUnew[1][3] = 0.0;  CUnew[1][4+1] = (1.0/9.0);   // + (1/9)dt M1
      CUnew[2][0]=0.0;     CUnew[2][1]=0.0;          CUnew[2][2]=2.0/3.0; CUnew[2][3] = 0.0;  CUnew[2][4+2] = (1.0/9.0);   // + (2/3)U2 + (1/9)dt M2
      CUnew[3][0]=0.0;     CUnew[3][1]=1.0/3.0;      CUnew[3][2]=0.0;     CUnew[3][3] = 0.0;  CUnew[3][4+3] = 0.0;         // + (1/3)dt F2

      //primitive values used for initial state, fluxes, final state (where you output)
      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0]; // produces U1
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
      pii[2]=p;    pbb[2]=pk[1];   pff[2]=pk[0]; // produces U3
      pii[3]=p;    pbb[3]=pk[0];   pff[3]=p; // produces U4 (only dUe part used)
    
      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
      uii[2]=GLOBALPOINT(uinitialglobal);  uff[2]=GLOBALPOINT(ulastglobal); ucum[2]=ucons;
      uii[3]=GLOBALPOINT(uinitialglobal);  uff[3]=GLOBALPOINT(ulastglobal); ucum[3]=ucons;

    }
    else if(TIMEORDER==3  && truestep){
      // IMEX2
      *numtimeorders=3;
      //      static FTYPE sqrt2=1.414213562373095048802;
      static FTYPE imp_gamma = 0.2928932188134524755992; // 1.0 - 1.0/sqrt2;
 
      // Ui ulastglobal dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0;                            CUf[0][1]=0.0;                            CUf[0][2]=0.0;      CUf[0][3] = 0.0;  CUf[0][4+0] = imp_gamma; // M0 no flux step : U0 = [Un] + \gamma dt M0
      CUf[1][0]=(3.0*imp_gamma-1.0)/imp_gamma;  CUf[1][1]=(1.0-2.0*imp_gamma)/imp_gamma;  CUf[1][2]=1.0;      CUf[1][3] = 0.0;  CUf[1][4+1] = imp_gamma; // computes F0 and M1 : U1 = [ (1/gamma)(3\gamma-1) Un + (1/\gamma)(1-2\gamma) U0 + dt F0] + \gamma dt M1
      CUf[2][0]=1.0/2.0;                        CUf[2][1]=1.0/2.0;                        CUf[2][2]=1.0/2.0;  CUf[2][3] = 0.0;  CUf[2][4+0] = imp_gamma; CUf[2][4+1]=(1.0-imp_gamma)/2.0; // F1 no new implicit step : U2 = [(1/2)Un + (1/2)U1 + (dt/2)F1 + \gamma dt M0 + (1/2)(1-\gamma) dt M1]
    
      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U3
      // NOTE: Could instead only use final U^{n+1} = U2 as possible for this case and setup above for final step, but current way is more local and really would require less memory (but not implemented to do that)
      CUnew[0][0]=0.5;   CUnew[0][1]=0.0;      CUnew[0][2]=0.0;  CUnew[0][3] = 0.0;  CUnew[0][4+0]=imp_gamma;             // + (1/2)Un + dt\gamma M0
      CUnew[1][0]=0.0;   CUnew[1][1]=0.0;      CUnew[1][2]=0.5;  CUnew[1][3] = 0.0;  CUnew[1][4+1]=(1.0-imp_gamma)/2.0;   // + (1/2)U1 + (1/2)(1-\gamma)M1
      CUnew[2][0]=0.0;   CUnew[2][1]=1.0/2.0;  CUnew[2][2]=0.0;  CUnew[2][3] = 0.0;                                       // + (1/2)dt F1
    
      //always starting the substeps from the initial time
      pii[0]=p;      pbb[0]=p;       pff[0]=pk[0]; // produces U1
      pii[1]=p;      pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
      pii[2]=p;      pbb[2]=pk[1];   pff[2]=p; // produces U3

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
      uii[2]=GLOBALPOINT(uinitialglobal);  uff[2]=GLOBALPOINT(ulastglobal); ucum[2]=ucons;

    }
    else if(TIMEORDER==2  && truestep){
      // IMEX1B = Backward Euler (1st order implicit, 2nd order explicit)
      *numtimeorders=2;

      // old ucons not used for this method (i.e. [?][1]=0)
      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=0.5; CUf[0][3] = 0.0;  CUf[0][4+0] = CUf[0][2];
      CUf[1][0]=1.0; CUf[1][1]=0.0; CUf[1][2]=1.0; CUf[1][3] = 0.0;  CUf[1][4+1] = CUf[1][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U2
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=0.0;  CUnew[0][3] = 0.0;  CUnew[0][4+0] = CUnew[0][1];
      CUnew[1][0]=0.0;   CUnew[1][1]=0.0;      CUnew[1][2]=1.0;  CUnew[1][3] = 0.0;  CUnew[1][4+1] = CUnew[1][1];

      pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
      pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;
      uii[1]=GLOBALPOINT(uinitialglobal);  uff[1]=GLOBALPOINT(ulastglobal); ucum[1]=ucons;
    }
    else if(TIMEORDER==1 ||  dt==0.0){ // dt==0.0 case is case when just passing through
      // Euler method = IMEX1 (1st order implicit, 1st order explicit)
      *numtimeorders=1;

      // Ui ulast dUe(pb) <timestuff> dUi(Uf)
      CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0; CUf[0][3] = 0.0;   CUf[0][4] = CUf[0][2];

      // Ui dUe(Ub) Uf <timestuff> dUi(Uf)
      // ucons=U1
      CUnew[0][0]=0.0;   CUnew[0][1]=0.0;      CUnew[0][2]=1.0;   CUnew[1][3] = 0.0;     CUnew[1][4] = CUnew[0][1];

      pii[0]=p;    pbb[0]=p;    pff[0]=p;

      // GODMARK: use of globals: just scratch anyways
      uii[0]=GLOBALPOINT(uinitialglobal);  uff[0]=GLOBALPOINT(ulastglobal); ucum[0]=ucons;

    }


  }
  



}









/////////////////////////////
//
// System boundary routines that calls other system routines and user routines
//
/////////////////////////////



// Who calls bound_anypstag():
// v1)
// step_ch.c: bound_evolveprim calls bound_allprim calls bound_anyallprim calls bound_anypstag
// step_ch.c: bound_beforeevolveprim calls bound_anyallprim calls bound_anypstag
// fluxctstag.c : bound_pstag calls bound_anypstag
//
// v2)
// But if bound pstag in fluxctstag.c so Bcent can be interpolated from Bstag, then don't need to do it again!
// fluxctstag.c : bound_pstag calls bound_anypstag




/// normal full bounding during evolution
/// NOTEMARK: Was previously also bounding pstag here, but not necessary since done by calling bound_pstag in fluxctstag.c
int bound_evolveprim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int boundvartype=BOUNDPRIMTYPE;



  bound_anyprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);
  if(unewisavg && BOUNDUNEW){
    // SUPERGODMARK:
    // if not outflow boundary, then can bound u as p
    // desirable since machine errors can be different and lead to secular changes
    // esp. in MPI
    bound_uavg(boundstage, finalstep, boundtime, boundvartype, ucons,pstag, ucons,doboundmpi);
  }


  //  bound_allprim(boundstage,finalstep, boundtime, prim,pstag,ucons,doboundmpi);

  return(0);

}


/// simple bounding during evolution
/// NOTEMARK: Was previously also bounding pstag here, but not necessary since done by calling bound_pstag in fluxctstag.c
int bound_beforeevolveprim(int boundstage,int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int boundvartype=BOUNDPRIMSIMPLETYPE; // tells boundary routines that only bounding 1 layer deep


  bound_anyprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);
  if(unewisavg && BOUNDUNEW){
    // SUPERGODMARK:
    // if not outflow boundary, then can bound u as p
    // desirable since machine errors can be different and lead to secular changes
    // esp. in MPI
    bound_uavg(boundstage, finalstep, boundtime, boundvartype, ucons,pstag, ucons,doboundmpi);
  }



  return(0);

}

/// normal bounding of non-staggered primitives
int bound_prim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int boundvartype=BOUNDPRIMTYPE;


  bound_anyprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);

  return(0);

}

/// normal bounding of staggered primitives
int bound_pstag(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int boundvartype=BOUNDPSTAGTYPE;

  bound_anypstag(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);

  return(0);

}


/// normal bounding of all primitives
int bound_allprim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int boundvartype=BOUNDPRIMTYPE;


  bound_anyallprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag,ucons,doboundmpi);

  return(0);

}




/// bound all prims
int bound_anyallprim(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int doboundmpi)
{
  int mystagboundvar;



  if(FLUXB==FLUXCTSTAG){
    bound_anyprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);
    if(boundvartype==BOUNDPRIMTYPE) mystagboundvar=BOUNDPSTAGTYPE;
    else if(boundvartype==BOUNDPRIMSIMPLETYPE) mystagboundvar=BOUNDPSTAGSIMPLETYPE;
    else  mystagboundvar=boundvartype;
    bound_anypstag(boundstage, finalstep, boundtime, mystagboundvar, prim, pstag, ucons,doboundmpi);
  }
  else{
    bound_anyprim(boundstage, finalstep, boundtime, boundvartype, prim, pstag, ucons,doboundmpi);
  }


  if(unewisavg && BOUNDUNEW){
    // SUPERGODMARK:
    // if not outflow boundary, then can bound u as p
    // desirable since machine errors can be different and lead to secular changes
    // esp. in MPI
    bound_uavg(boundstage, finalstep, boundtime, boundvartype, ucons,pstag, ucons,doboundmpi);
  }

  return(0);

}



/// bound uavg
int bound_uavg(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int mystagboundvar;
  FTYPE (*uavg)[NSTORE2][NSTORE3][NPR];




  // assign uavg
  uavg=ucons;


  // for now only worry about bounding staggered fields since that's what can do given "general" boundary condition code
  if(1||DOENOFLUX != NOENOFLUX){

    // CHANGINGMARK: DEBUG:
    // can modify for diag call if choose to avoid if outflow condition
    // or do something simple for outflow just for diagnostics
#if(FULLOUTPUT!=0 && PRODUCTION==0)
    bound_anyprim(boundstage, finalstep, boundtime, boundvartype, uavg,pstag, uavg,doboundmpi);
#endif


    if(DOENOFLUX == ENOFINITEVOLUME){
      // then need to bound all conservative quantities
      // SUPERGODMARK -- CHANGINGMARK -- need to generalize so bound knows if p or U quantity
      // other methods have only field "averaged" so only need to modify it
      bound_anyprim(boundstage, finalstep, boundtime, boundvartype, uavg, pstag, uavg,doboundmpi);
    }

    if(FLUXB==FLUXCTSTAG){
      // bound unew for staggered fields
      // SUPERGODMARK -- CHANGINGMARK: need to tell bound if p or U
      if(boundvartype==BOUNDPRIMTYPE) mystagboundvar=BOUNDPSTAGTYPE;
      else if(boundvartype==BOUNDPRIMSIMPLETYPE) mystagboundvar=BOUNDPSTAGSIMPLETYPE;
      else  mystagboundvar=boundvartype;

      bound_anypstag(boundstage, finalstep, boundtime, mystagboundvar, uavg, pstag, uavg,doboundmpi);
    }
  }

  return(0);

}





/// bound any prim type
/// CALLS directly real and MPI boundary functions
/// pstag = pstagglobal
/// ucons = unewglobal
int bound_anyprim(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pl,pliter;
  int dir;




  if(DOGRIDSECTIONING){
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      MYFUN(bound_gridsectioning(CENTEREDPRIM,prim,pstag,ucons,finalstep),"step_ch.c:bound_pstag()", "bound_pstag_user()", 1);
    }
  }



  DIMENLOOP(dir){

    if(0 && FLUXB==FLUXCTSTAG){ // apparently need to bound both fields SUPERGODMARK

      PALLLOOP(pl) whichpltoavg[pl]=1;
      PLOOPBONLY(pl) whichpltoavg[pl]=0;

      PALLLOOP(pl) ifnotavgthencopy[pl]=0;
      PLOOPBONLY(pl) ifnotavgthencopy[pl]=0;

      addremovefromanynpr(REMOVEFROMNPR, whichpltoavg, ifnotavgthencopy, &nprboundstart, &nprboundend, nprboundlist, &nprlocalstart, &nprlocalend, nprlocallist, NULL, NULL);
    }
    else{
      // no change then, doing all variables
    }



    // pre-post MPI recv's
    if(doboundmpi){
      MYFUN(bound_mpi_dir(boundstage, finalstep, -dir, boundvartype, prim, NULL, NULL, NULL,NULL),"step_ch.c:bound_prim()", "bound_mpi()", 1);
    }

    // real boundary zones
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      MYFUN(bound_prim_user_dir(boundstage, finalstep, boundtime, dir, boundvartype, prim),"step_ch.c:bound_prim()", "bound_prim_user()", 1);
    }// end if stage0 or stagem1
  
    if(doboundmpi){
      MYFUN(bound_mpi_dir(boundstage, finalstep, +dir, boundvartype, prim, NULL, NULL, NULL,NULL),"step_ch.c:bound_prim()", "bound_mpi()", 1);
    }

    // real boundary zones
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      int ispstag=BOUNDPRIMLOC;
      MYFUN(bound_prim_user_after_mpi_dir(boundstage, finalstep, boundtime, dir, boundvartype, ispstag, prim),"step_ch.c:bound_prim()", "bound_prim_user_after_mpi()", 1);
    }// end if stage0 or stagem1


    ////////////////////////////////////////////
    //
    // restore choice for bounding

    if(0 && FLUXB==FLUXCTSTAG){ // apparently need to bound both fields SUPERGODMARK
      addremovefromanynpr(RESTORENPR, whichpltoavg, ifnotavgthencopy, &nprboundstart, &nprboundend, nprboundlist, &nprlocalstart, &nprlocalend, nprlocallist, NULL, NULL);
    }
  }

  
  return(0);
}


/// bound any pstag type
/// CALLS directly real and MPI boundary functions
/// pstag is at FACE1,2,3 for fields, so user bound is different
/// MPI bounding is the same as CENT quantities
/// used when restarting in initbase.c to bound unewglobal for FLUXB==FLUXCTSTAG
int bound_anypstag(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi)
{
  int pl,pliter;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int dir;
  int mystagboundvar;


  if(boundvartype==BOUNDPRIMTYPE) mystagboundvar=BOUNDPSTAGTYPE;
  else if(boundvartype==BOUNDPRIMSIMPLETYPE) mystagboundvar=BOUNDPSTAGSIMPLETYPE;
  else  mystagboundvar=boundvartype;





  if(DOGRIDSECTIONING){
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      MYFUN(bound_gridsectioning(STAGGEREDPRIM,prim,pstag,ucons,finalstep),"step_ch.c:bound_pstag()", "bound_pstag_user()", 1);
    }
  }




  DIMENLOOP(dir){
  
    ////////////////////////////////////////////
    //
    // save choice for bounding
    nprlocalstart=nprboundstart;
    nprlocalend=nprboundend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=nprboundlist[pl];

    // choose range of PBOUNDLOOP and PLOOPMPI
    nprboundstart=0;
    nprboundend=2;
    nprboundlist[0]=B1;
    nprboundlist[1]=B2;
    nprboundlist[2]=B3;


    // pre-post MPI recvs
    if(doboundmpi){
      MYFUN(bound_mpi_dir(boundstage, finalstep, -dir, mystagboundvar, pstag, NULL, NULL, NULL,NULL),"step_ch.c:bound_pstag()", "bound_mpi()", 1);
    }

    // real boundary zones
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      MYFUN(bound_pstag_user_dir(boundstage, finalstep, boundtime, dir,mystagboundvar,pstag),"step_ch.c:bound_pstag()", "bound_pstag_user()", 1);
    }// end if stage0 or stagem1


    if(doboundmpi){
      MYFUN(bound_mpi_dir(boundstage, finalstep, +dir, mystagboundvar, pstag, NULL, NULL, NULL,NULL),"step_ch.c:bound_pstag()", "bound_mpi()", 1);
    }

    // real boundary zones
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      int ispstag=BOUNDPSTAGLOC;
      MYFUN(bound_prim_user_after_mpi_dir(boundstage, finalstep, boundtime, dir, mystagboundvar, ispstag, pstag),"step_ch.c:bound_pstag()", "bound_prim_user_after_mpi()", 1);
    }// end if stage0 or stagem1


    ////////////////////////////////////////////
    //
    // restore choice for bounding
    nprboundstart=nprlocalstart;
    nprboundend=nprlocalend;
    PMAXNPRLOOP(pl) nprboundlist[pl]=nprlocallist[pl];
  }


  return(0);
}


/// bound any flux type
/// CALLS directly real and MPI boundary functions
/// special bound for flux that only bounds in direction of flux itself (so not full cross-flux information)
/// used for finite difference version of ENO
int bound_flux(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], int doboundmpi)
{
  int dir;



  if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE){
    // then can stay
  }
  else{
    dualfprintf(fail_file,"Shouldn't be in bound_flux() with boundvartype=%d\n",boundvartype);
    myexit(2467367);
  }


  if(DOENOFLUX==ENOFLUXRECON && BOUNDFLUXRECON==0){
    dualfprintf(fail_file,"DEBUG: Assuming bound_flux() called only for debugging purposes\n");
  }


  // pre-post recv's
  if(doboundmpi){
    MYFUN(bound_mpi(boundstage, finalstep, -1, boundvartype, NULL, F1, F2, F3, NULL),"step_ch.c:bound_flux()", "bound_mpi()", 1);
  }

  
  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    MYFUN(bound_flux_user(boundstage, finalstep, boundtime, boundvartype,F1,F2,F3),"step_ch.c:bound_flux()", "bound_flux_user()", 1);
  }// end if stage0 or stagem1


  if(doboundmpi){
    MYFUN(bound_mpi(boundstage, finalstep, +1, boundvartype, NULL, F1, F2, F3, NULL),"step_ch.c:bound_flux()", "bound_mpi()", 1);
  }


  return(0);
}


/// used when restarting in initbase.c
int bound_vpot(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int doboundmpi, int doboundnonmpi)
{
  int dir;



  if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE){
    // then can stay
  }
  else{
    dualfprintf(fail_file,"Shouldn't be in bound_vpot() with boundvartype=%d\n",boundvartype);
    myexit(3483466);
  }


  // pre-post MPI recv's
  if(doboundmpi){
    MYFUN(bound_mpi(boundstage, finalstep, -1, boundvartype, NULL, NULL, NULL, NULL, vpot),"step_ch.c:bound_flux()", "bound_mpi()", 1);
  }


  // real boundary zones
  if(doboundnonmpi){
    if((boundstage==STAGE0)||(boundstage==STAGEM1)){
      MYFUN(bound_vpot_user(boundstage, finalstep, boundtime, boundvartype,vpot),"step_ch.c:bound_vpot()", "bound_vpot_user()", 1);
    }// end if stage0 or stagem1
  }

  if(doboundmpi){
    MYFUN(bound_mpi(boundstage, finalstep, +1, boundvartype, NULL, NULL, NULL, NULL, vpot),"step_ch.c:bound_flux()", "bound_mpi()", 1);
  }


  return(0);
}







/// bound pflag type
/// CALLS directly real and MPI boundary functions
int bound_pflag(int boundstage, int finalstep, SFTYPE boundtime, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS], int doboundmpi)
{
  int boundvartype=BOUNDINTTYPE;
  



  if(doboundmpi){

    if(UTOPRIMFIXMPICONSISTENT==1){
      MYFUN(bound_mpi_int(boundstage, finalstep, -1, boundvartype, prim),"step_ch.c:bound_pflag()", "bound_mpi_int()", 1);
    }
    else{
      // need to fill boundary cells with failure
      // otherwise, would have to go deeper into fixups and dependency chain for the UTOPRIMFIXMPICONSISTENT chocie would become too deep
      MYFUN(bound_mpi_int_fakeutoprimmpiinconsisent(boundstage, finalstep, -1, boundvartype, prim,UTOPRIMFAILFAKEVALUE),"step_ch.c:bound_pflag()", "bound_mpi_int()", 1);
    }

  }




  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    MYFUN(bound_pflag_user(boundstage, finalstep, boundtime, boundvartype, prim),"step_ch.c:bound_pflag()", "bound_pflag_user()", 1);

  }// end bound stage


  if(doboundmpi){

    if(UTOPRIMFIXMPICONSISTENT==1){
      MYFUN(bound_mpi_int(boundstage, finalstep, +1, boundvartype, prim),"step_ch.c:bound_pflag()", "bound_mpi_int()", 1);
    }
    else{
      // need to fill boundary cells with failure
      // otherwise, would have to go deeper into fixups and dependency chain for the UTOPRIMFIXMPICONSISTENT chocie would become too deep
      MYFUN(bound_mpi_int_fakeutoprimmpiinconsisent(boundstage, finalstep, +1, boundvartype, prim,UTOPRIMFAILFAKEVALUE),"step_ch.c:bound_pflag()", "bound_mpi_int()", 1);
    }

  }

  return(0);

}









