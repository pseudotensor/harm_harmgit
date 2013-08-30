#include "decs.h"

// TODO:
// 1) redo horizon enerregion to use Sasha method for moving horizon


// OPENMPMARK: None of these functions should be called by multiple threads.

// static declarations
static int get_new_metric_parms(SFTYPE *MBHpar, SFTYPE *apar, SFTYPE *QBHpar, SFTYPE *EP3par, SFTYPE *THETAROTpar);
static int compute_phi_self_gravity_simple(FTYPE (*pb)[NSTORE2][NSTORE3][NPR]);
static int get_new_metric_parms(SFTYPE *MBHpar, SFTYPE *apar, SFTYPE *QBHpar, SFTYPE *EP3par, SFTYPE *THETAROTpar);
static int interpX_phi(FTYPE *X, SFTYPE *phivsr_tot, SFTYPE *phi);



//#define MAXSKIP 2
#define MAXSKIP 100 // DEBUG: GODMARK: testing -- should work if get gravitydt to be correct

// how much not to trust smallness of gravitydt
#define ENHANCEFACTOR (1.0)

// used to control metric update method between longsteps and substeps
// only changed after each full time step
void control_metric_update(void)
{
  int oldgravityskipstep;
  int newgravityskipstep;

  if(DOEVOLVEMETRIC && (EVOLVEMETRICSUBSTEP==2 || EVOLVEMETRICSUBSTEP==0) ){
    // decide whether to use longsteps or substeps method

    // first MINIMIZE over all CPUs
    //    sourcedtglobal and  wavedtglobal not used since already that information is in dt that is already known to all CPUs the same
    // find global minimum value of gravitydtglobal over all cpus
    mpifmin(&gravitydtglobal);
    
    oldgravityskipstep=gravityskipstep;
    // no more than MAXSKIP skips since after that not computationally a burden
    // trying to account for metric changes when not doing self-gravity (changes in metric parameters)
    // alternative is to modify the skipstep factor also in the longstep part that checks on metric updates
    newgravityskipstep=(int)(MIN(MAX((gravitydtglobal/(ENHANCEFACTOR*dt)),0),MAXSKIP));

    // don't increase skipsteps by more than a factor of 2 per timestep
    if(newgravityskipstep>2*oldgravityskipstep+1) gravityskipstep=2*oldgravityskipstep+1;
    else gravityskipstep=newgravityskipstep;

#if(1)
    if(gravityskipstep>=5){ // not much point if skipping every other step, so limit to 5 steps to skip
      // if here, then gravityskipstep is used for longsteps to help determine when to update metric
      doevolvemetricsubsteps=0;
    }
    else{
      // then gravity is changing as much as anything else
      doevolvemetricsubsteps=1;
    }
#else
    // DEBUG
    // GODMARK DEBUG
    // fix for now as test
    gravityskipstep=2;
    doevolvemetricsubsteps=0;    
#endif
  }

  //  trifprintf("gravitydt %21.15g dt=%21.15g ratio=%21.15g\n",gravitydtglobal,dt,gravitydtglobal/dt);
  //  trifprintf("gravity: old=%d new=%d set=%d\n",oldgravityskipstep,newgravityskipstep,gravityskipstep);
  
}





// set_grid uses this
// notice that store_old_metric called directly there initially to store first metric
// This function sets metric time and chooses whether to store metric
// If this function is called, then about to compute new metric
void control_time_store_metric(int whichtime, FTYPE *Cunew)
{
  extern int store_old_metric(void);
  static FTYPE lastCunew, lastdt;
  static int firsttime=1;
  FTYPE dtcheck;
  FTYPE currenttime;
  static long nsteplaststore;


  if(whichtime==0 || DOEVOLVEMETRIC==0){
    // definition of Xmetric's should come before any calls to coord()
    // then initial time and (e.g.) don't have old metric if DOEVOLVEMETRIC==1
    Xmetricnew[TT] = Xmetricold[TT] = t;
    lastCunew = 0.0;
  }// otherwise Xnew/Xold should have been set correctly



  if(DOEVOLVEMETRIC && whichtime==1){

#if(0)
    // then last time computed metric is at Xmetricnew, so set as old since to be overwritten with new
    if(doevolvemetricsubsteps || EVOLVEMETRICSUBSTEP==1){

      // old metric is only stored once on first substep
      // always want reasonable old metric to take time derivatives
      // needs to never be duplicate of new time and can't be 1.0 since next step will give wrong dg/dt (i.e. dg/0)
      if(Cunew[3]!=1.0 && ( (steppart>0 && fabs(Cunew[3]-lastCunew)>1E-3) || (steppart==0) ) ){
        //      if(steppart==0){// try only storing on first substep
        Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
        // store old grid if doing DOEVOLVEMETRIC==1 
        // do only when changing metric
        //dualfprintf(fail_file,"Storing metric if changing metric: doevolvemetricsubsteps=%d steppart=%d Cunew[3]=%21.15g lastCunew=%21.15g\n",doevolvemetricsubsteps,steppart,Cunew[3],lastCunew);
        store_old_metric();
      }
      // new metric is always at new time
      Xmetricnew[TT] = t + Cunew[3]*dt; // present g is at this time
      //      Xmetricnew[TT] = t;
      // dualfprintf(fail_file,"Xmetricold=%21.15g Xmetricnew=%21.15g t=%21.15g Cunew[3]=%21.15g dt=%21.15g\n",Xmetricold[TT],Xmetricnew[TT],t,Cunew[3],dt);
      // update lastCunew
      lastCunew = Cunew[3];
    }
    else{
      // if evolving metric, then store old metric before computing new one
      // store metric, needed to have dg/dt terms in connection coefficients
      // store old grid if doing DOEVOLVEMETRIC==1 
      // do only when changing metric
      //dualfprintf(fail_file,"Storing metric if changing metric (not on substep)\n");
      store_old_metric();


      Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
      Xmetricnew[TT] = t; // present g
    }
#elif(0)
    // old-way where metric only affects things (wrongly) on first substep
    // GODMARK: works for now, but should use above that doesn't yet work
    //dualfprintf(fail_file,"Storing metric if changing metric (not on substep -- testing)\n");
    store_old_metric();
    
    Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
    Xmetricnew[TT] = t; // present g
#elif(0)
    // old-way where metric only affects things (wrongly) on first substep
    // GODMARK: works for now, but should use above that doesn't yet work
    //dualfprintf(fail_file,"Storing metric if changing metric (not on substep -- testing)\n");
    store_old_metric();
    
    Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
    Xmetricnew[TT] = t + Cunew[3]*dt; // present g
#elif(0)

    if(steppart==0 && Cunew[3]!=1.0){
      // old-way where metric only affects things (wrongly) on first substep
      // GODMARK: works for now, but should use above that doesn't yet work
      //dualfprintf(fail_file,"Storing metric if changing metric (not on substep -- testing)\n");
      store_old_metric();
      
      Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
    }
    Xmetricnew[TT] = t + Cunew[3]*dt; // present g
#elif(0)

    if(firsttime){
      lastdt = dt;
      firsttime=0;
    }

    currenttime=t + Cunew[3]*dt; // present g
    // update metric if gone sufficiently far in time (use lastdt in case dt increases alot)
    dtcheck=1.0*MIN(dt,lastdt);
    if(currenttime - Xmetricold[TT] > dtcheck ){
      // old-way where metric only affects things (wrongly) on first substep
      // GODMARK: works for now, but should use above that doesn't yet work
      //dualfprintf(fail_file,"Storing metric if changing metric (not on substep -- testing)\n");
      store_old_metric();
      
      Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
      lastdt=dt;
    }
    Xmetricnew[TT] = currenttime;

    if(Xmetricnew[TT]==Xmetricold[TT]){
      dualfprintf(fail_file,"WTF\n");
      myexit(2166);
    }
#elif(1)

    // NO NEED TO DIFFERENTIATE BETWEEN STEPPING METHODS, THIS GENERALLY STORES IF NECESSARY
    if(firsttime){
      lastdt = dt;
      nsteplaststore=nstep;
    }

    currenttime=t + Cunew[3]*dt; // present g

    // update metric if gone sufficiently far in time (use lastdt in case dt increases alot)
    dtcheck=0.01*MIN(dt,lastdt);
    if(nsteplaststore<nstep && currenttime - Xmetricold[TT] > dtcheck && currenttime - Xmetricnew[TT] > dtcheck ){
      // old-way where metric only affects things (wrongly) on first substep
      // GODMARK: works for now, but should use above that doesn't yet work
      //dualfprintf(fail_file,"Storing metric if changing metric (not on substep -- testing)\n");
      store_old_metric();
      //dualfprintf(fail_file,"Storing old metric during evolution\n");
      
      Xmetricold[TT] = Xmetricnew[TT]; // old g is at that time (time computed previously)
      nsteplaststore=nstep;
      lastdt=dt;
    }
    Xmetricnew[TT] = currenttime;

    //    if(firsttime==0 && Xmetricnew[TT]==Xmetricold[TT]){
    //  dualfprintf(fail_file,"WTF: steppart=%d nstep=%ld :: %21.15g %21.15g %d :: %21.15g\n",steppart,nstep,Xmetricold[TT],Xmetricnew[TT],firsttime,dtcheck);
    //  myexit(2166);
    // }

    firsttime=0;
#endif


  }




}





// whether to use true energy accreted (SHOULD! -- especially when magnetic torques exist)
//
// use 0 gives perfect conservation of mass through conservation of baryon number
// but while baryon number must be conserved it is NOT the conserved (ADM) mass
// 
//
// have to synch metric Mvsr, phivsr, and this so get good conservation of desired quantity
//
// if 0, then probably using USEUNEW==2
// if 1, then probably using (for now) (USENEW==1) or (USENEW==0 with MHDTTTT==1,2,3) until do more correct thing
#define USETRUEENERGYACCRETED 1

// functions that define type of energy accreted or written into metric
//
// 1) static int get_new_metric_parms(SFTYPE *MBHpar, SFTYPE *apar, SFTYPE *QBHpar)
//
// below just uses Mvsr as rest-mass and energy terms
// 2)  int action_ifchange_horizoni(int horizontiold,int horizontinew)
//
// 3) static int compute_phi_self_gravity_simple(FTYPE (*pb)[NSTORE2][NSTORE3][NPR])
//
// 4) static int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself)
//
// 5) static int set_gcov_bl_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself)


static int get_new_metric_parms(SFTYPE *MBHpar, SFTYPE *apar, SFTYPE *QBHpar, SFTYPE *EP3par, SFTYPE *THETAROTpar)
{
  int enerregion;  
  FTYPE j,dj;
  SFTYPE (*localpcum)[NPR];
  SFTYPE (*localpcum_tot)[NPR];
  SFTYPE (*localpdot)[NPR];
  SFTYPE (*localpdot_tot)[NPR];
  int dir;




  /////////////////////////////////////
  //
  // now integrate over all CPUs
  // also done in dump_ener.c so diagonstics updated (doesn't change the evolution to update this more frequently)
  //
  enerregion=OUTSIDEHORIZONENERREGION;
  localpdot=pdotreg[enerregion];
  localpdot_tot=pdotreg_tot[enerregion];
  localpcum=pcumreg[enerregion];
  localpcum_tot=pcumreg_tot[enerregion];
  //
  // really only need dir=X1DN here, but doesn't hurt to do all in case needed later
  //DIRLOOP(dir)
  if(integrate((COMPDIM*2)*NPR,&localpdot[0][0],&localpdot_tot[0][0],SURFACETYPE,enerregion)>=1) return(1);
  //DIRLOOP(dir)
  if(integrate((COMPDIM*2)*NPR,&localpcum[0][0],&localpcum_tot[0][0],SURFACETYPE,enerregion)>=1) return(1);
  // below have been subsumed into own full enerregion
  //  if(integrate(NPR,&horizonflux[0],&horizonflux_tot[0],SURFACETYPE,enerregion)>=1) return(1);
  //if(integrate(NPR,&horizoncum[0],&horizoncum_tot[0],SURFACETYPE,enerregion)>=1) return(1);

  // now can evolve MBH, a, and QBH using localpcum
  // evolve MBH = MBH0 + localpcum_tot[UU]
  // evolve a = a0 + localpcum_tot[U3]

  // if present actual value of mass accreted is greater than 1% more than present metric value used, then update metric value

  // need to get units right
  // \int_t \dot{E} dt = localpcum_tot[UU];
  // \int_t \dot{L} dt = localpcum_tot[U3];


  if(myid==0){
    //////////////////////////////
    //
    // get metric change
    // GODMARK : Unsure how MBH and a should change with accretion
    // issue is mass appearing in metric vs. gravitational mass

    // dE=dM since C=1
#if(USETRUEENERGYACCRETED)
    // T^r_t \sim \rho u^r u_t \sim -\rho u^r, so if u^r<0, then T^r_t\sim \rho.
    // So positive T^r_t means adding energy to black hole
    // dE is in energyunit, while MBH is in Lunit with G=c=1, so need to convert dE to length unit
    dE = localpcum_tot[X1DN][UU]*Mfactor;
  
    if(REMOVERESTMASSFROMUU==2){
      // then need to add back in the mass
      dE += -localpcum_tot[X1DN][RHO]*Mfactor;
    }
    else{
      // then already correct
    }
#else

    // total energy changes total gravitational mass of BH
    // proxy for ADM mass
    dE = -localpcum_tot[X1DN][RHO]*Mfactor;

#endif

    // T^r_\phi \sim \rho u^r u_\phi \sim \rho u^r J, so if u^r<0 and J>0, then adding angular momentum to black hole
    // So negative T^r_\phi means adding positive angular momentum to black hole
    // dJ is angular momentum unit, convert to Length^2 type unit
    dJ = -localpcum_tot[X1DN][U3]*Jfactor;
  
    // dj = dJ/M^2 - 2j dE/M // from GSM04, where dJ is minus mine
    // dj = dJ M^2 + 2j dE/M // what I get for dJ and dE added to BH (seems to have overall sign difference from above)
    // da = -dJ + 2a dE // from gammie.m (probably a here is dimensionless and M=1)
  
    // da = dJ/M - a/M dE // what I get
    // all terms should be in metric-Length units
    // dJ in Ltilde^2 units, MBH in Ltilde units, a in Ltilde units, and dE in Ltilde units
    //  dabh = dJ/MBH - a*(dE/MBH) ; // first term is adding spin, while second term accounts for adding raw energy onto present spinning object (which slows the spin down)
  
    // really need discrete (finite) difference form, as below
    dabh = (a0*MBH0+dJ)/(MBH0+dE+SMALL) - a0;


    j = a/(MBH+SMALL); // present dimensionless Kerr parameter
    // below in discrete form
    dj = (a0*MBH0+dJ)/pow(MBH0+dE+SMALL,2) - a0/(MBH0+SMALL) ;
    // below in infinitesimal form
    //  dj = dJ/(MBH*MBH) - 2.0*j*(dE/MBH);



    // unewglobal holds final averaged conserved quantities from last full step, so just reprocess this after changing grid
    // unewglobal is = $\detg U$ = conserved quantity within cell
    

    // DEBUG
    //    trifprintf("MetricParms:: MBH=%21.15g to %21.15g and a=%21.15g to %21.15g and j=%21.15g to %21.15g\n",MBH,MBH0+dE,a,a0+dabh,j,a0/(MBH0+SMALL)+dj);


    // update MBH, a, and QBH
    *MBHpar=MBH0 + dE;
    *apar=a0 + dabh; // a = J/MBH.  Since used MBH to normalize dJ and dE above, then a is in correct per Lunit with G=c=1
    *QBHpar=QBH0;
    *EP3par=EP30;
    *THETAROTpar=THETAROT0;
  }

#if(USEMPI)
  MPI_Bcast(MBHpar,1,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(apar,1,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(QBHpar,1,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(EP3par,1,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(THETAROTpar,1,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
#endif



  


  return(0);

}





#define FRACCHANGE (0.01)
#define ABSJCHANGE (0.1)

// called every full timestep inside step_ch.c:step_ch_full()
// used to determine value of mass, spin, and charge
int compute_new_metric_longsteps(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  int dochange, dochangesend;
  SFTYPE MBHpar, apar, QBHpar, EP3par, THETAROTpar;
  static long lastupdatenstep;
  static int firsttime=1;


  if(firsttime){
    lastupdatenstep=nstep;
  }



#if(DOEVOLVEMETRIC==0 || EVOLVEMETRICSUBSTEP==1)
  return(0);
#endif
#if(DOEVOLVEMETRIC==1 && EVOLVEMETRICSUBSTEP==2)
  // then check if wanting to do longstep method
  if(doevolvemetricsubsteps==1) return(0);
#endif



  ///////////
  //
  // get new Kerr metric parameters as determined by flux through inner-radial boundary
  // if self-gravity is on, then this is just a way to get parameters to detect whether need to update metric rather than being update itself
  //
  ////////////
  get_new_metric_parms(&MBHpar, &apar, &QBHpar, &EP3par, &THETAROTpar);





  ///////////////////////////////
  //
  // see if metric change is sufficient to actually change the metric
  // dochange will become 1 of either metric on-grid changed due to d/dt part of gravity or just due to metric parameters changing
  // assume no change by default
  //
  ///////////////////////////////
  dochange=0;

  if(DOEVOLVEMETRIC){
    if(
       (fabs(MBH0+dE)>MBH*(1+FRACCHANGE))||
       //       ((fabs(a0+dabh)>a*(1+FRACCHANGE))&&(fabs(dj)>ABSJCHANGE)) // change if fractionally larger by 1% or a changes by dj=0.1, since j has an absolute scale
       (fabs(a0+dabh)>a*(1+FRACCHANGE)) // change if fractionally larger by 1%
       ){
      dochange=1;
    }
    if(DOEVOLVEMETRIC && (EVOLVEMETRICSUBSTEP==2 || EVOLVEMETRICSUBSTEP==0) ){
      // check if nstep is such that time to perform metric update
      // Notice that if gravityskipstep radically decreases, then condition will get triggered sooner.  That is, gravityskipstep as chosen when last did metric update does NOT control next update
      // this allows more flexibility so that gravityskipstep can change alot and we can still keep up with metric changes
      if(nstep>lastupdatenstep+gravityskipstep){
        // e.g. if lastupdatenstep=0 and nstep=3 and gravityskipstep=2, then time to change!
        dochange=1;
      }
    }

    if(!dochange && nstep%100==0){
      trifprintf("NOT Changing metric from MBH=%21.15g to %21.15g and a=%21.15g to %21.15g\n",MBH,MBH0+dE,a,a0+dabh);
    }

  }


#if(USEMPI)
  // condition above should be reasonable, but problematic near machine errors, so enforce all CPUs to have same condition
  dochangesend=dochange;
  MPI_Allreduce(&dochangesend, &dochange,1, MPI_INT, MPI_MAX,MPI_COMM_GRMHD);
#endif






  /////////////
  //
  // If BH should change mass or spin, then change it
  //
  /////////////

  if(dochange){
    lastupdatenstep=nstep;

    // unewglobal holds final averaged conserved quantities from last full step, so just reprocess this after changing grid
    // unewglobal is = $\detg U$ = conserved quantity within cell
    
    trifprintf("AM Changing metric from MBH=%21.15g to %21.15g and a=%21.15g to %21.15g\n",MBH,MBH0+dE,a,a0+dabh);

    // 1 below tells us not the initial time
    compute_new_metric_and_prims(1,MBHpar,apar,QBHpar,EP3par,THETAROTpar,prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct, F1, F2, F3, Atemp, uconstemp);

  }
  else{
    // DEBUG
    // GODMARK
    // for now, compute self-gravitational potential so can see in SM how evolves
    // ??should set output file of phivsr, etc. to as often as need to update metric??

    //    if(DOSELFGRAVVSR){
    // must be called after MHD solution set
    // determine phi(r) , which will be used by set_grid to determine metric
    //      compute_phi_self_gravity_simple(p);// use standard primitive
    //}
  }


  firsttime=0;
  return(0);
}














///////////////////                                                                                                            
//                                                                                                                             
// NOTES ABOUT compute_new_metric_and_prims() in initbase.c                                                                    
//                                                                                                                             
//  if(DOSELFGRAVVSR && EVOLVEMETRICSUBSTEP==0){ // since if EVOLVEMETRICSUBSTEP==1 then getting proper metric at beginning of substep                                                                                                                         
// apparently time-change of metric too large if avoiding this step when EVOLVEMETRICSUBSTEP==1 since inversion failures occur when evolving                                                                                                                   

// now that have conserved quantities and primitive MHD state, can determine full metric with self-gravity if wanted           
// if no self-gravity then this will produce same result as original set_grid() above, so not needed                           
// this needs to be done after init and after restart so that have metric with self-gravity                                    

// 0 above tells function that we are still at initial time                                                                    
// initial time does 2 things:                                                                                                 
//    a) dg/dt not computed between set_grid(0) and this function that will call set_grid(0)                                   
//    b) initial conditions of primitive quantities (rho,u,v^i) are assumed to NOT change since assumed user knew what potential was                                                                                                                           
//    c) GODMARK: B needs to have divB=0 preserved, so should recompute B from A_\phi using new grid                           
//                                                                                                                             
// if used 1, then dg/dt would be computed from difference between metric with and without self-gravity terms                  

// if RESTARTMODE==1, then assume primitives and conserved quantities are correctly read-in and just need to compute metric without changing conserved quantities or primitive quantities                                                                      

// since if EVOLVEMETRICSUBSTEP==1 then getting proper metric at beginning of substep                                          







// used to evolve metric using new Kerr parameters and updating self-gravity if wanted
// after metric is updated, then use conserved quantities to obtain corrected primitives that are consistent with that metric
// called in function above when metric should change
// called in initbase.c at t=0 when doing self-gravity.  That is, this function is called after primitives and conserved quantities are assigned so self-gravity potential can be determined.
// this can be called for only self-gravity or only evolving metric
//
// GODMARK: currently not using input parameters
// used with EVOLVEMETRICSUBSTEP==0 || ==2
int compute_new_metric_and_prims(int whichtime, SFTYPE MBHpar, SFTYPE apar, SFTYPE QBHpar, SFTYPE EP3par, SFTYPE THETAROTpar, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR])
{
  SFTYPE dtold;
  FTYPE CUf=0.0, Cunew=0.0;


  // tells code, if necessary, that doing special non-temporal update
  specialstep=1; // signifies to rest of code below that evolving on long steps
  compute_new_metric_anystep(whichtime, &CUf, &Cunew, prim, ucons);



  if(whichtime!=0){
    ////////////////////////////////////////////////
    //
    // NOW: take fake step so grid effect is processed without using non-updated primitives for fluxes so we are conservative
    // just way to correct primitives with new conserved quantities
    // can't just use old primitives since that's the effect we are trying to account for
    // note that standard advance is forced to use prior ucum that becomes ucum and is processed
    // note that finitevolume advance already uses prior ucum for Ui which becomes new ucum and is processed
    dtold=dt;
    dt=0.0;
    int truestep=0;
    step_ch_full(truestep,prim,pstag,ucons,vpot,Bhat,pl_ct, pr_ct,F1,F2,F3,Atemp,ucons);
    //    dump(9999);
    dt=dtold;
    //    realnstep--;  // not necessary since don't iterate if dt==0.0
    //nstep--;


    // GODMARK: if using WENO method need to make conserved quantities consistent with primitives when avoiding converting to new primitives at t=0
    // basically need to recompute "conserved quantities" post self-gravity rather than only before
  }

  //////////////////
  // now set some things post-new primitives

  //////
  // set dt since primitive changed and self-gravity may have changed
  set_dt(prim,&dt);



  return(0);

}


// returns whether doing metric substep method or not
int doingmetricsubstep(void)
{

#if(DOEVOLVEMETRIC==0 || EVOLVEMETRICSUBSTEP==0)
  return(0);
#endif

#if(DOEVOLVEMETRIC==1 && EVOLVEMETRICSUBSTEP==2)
  // then check if wanting to do substep method
  // can only switch after full step
  if(doevolvemetricsubsteps==0) return(0);
#endif

  return(1);
}





// compute new metric if doing it on a substep of the hydrodynamical step
int compute_new_metric_substep(FTYPE *CUf, FTYPE *Cunew, FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{
  int retvalue;

  retvalue=doingmetricsubstep();
  if(retvalue==0) return(0);


  specialstep=2; // signifies to rest of code below that evolving on substeps
  compute_new_metric_anystep(1,CUf, Cunew,pb,ucons); // 1 tells us not initial time
  specialstep=0; // back to normal


  return(0);

}





// like compute_new_metric_and_prims(), but instead is called during substep after conserved quantities are updated and before inversion
// this then avoids need to recompute primitives as above function that's used only infrequently
// used with EVOLVEMETRICSUBSTEP==1 || ==2

// CUf and Cunew are for if ever actually evolving metric
// metric update doesn't directly give dt limit.  dt limit computed in advance.c from how source terms change stress-energy tensor
int compute_new_metric_anystep(int whichtime,FTYPE *CUf, FTYPE *Cunew,FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{
  extern int post_par_set(void);
  int ii;
  static int mbhcount=500;
  int action_ifchange_horizoni(int horizontiold,int horizontinew);
  int horizontiold, horizontinew;

  




  ////////////////////////////////////
  //
  // compute self-gravity contribution to space-time
  //
  ////////////////////////////////////
  if(DOSELFGRAVVSR){
    // must be called after MHD solution set
    // determine phi(r) , which will be used by set_grid to determine metric
    while(compute_phi_self_gravity_simple(pb)==-1){
    }
  }





  // DEBUG
  //  if(nstep>1051){
  //  GRAVLOOP(ii){
  //    dualfprintf(fail_file,"ii=%d phivsr=%21.15g Mvsr=%21.15g Jvsr=%21.15g vrsqvsr=%21.15g\n",ii,phivsr_tot[ii],Mvsr_tot[ii],Jvsr_tot[ii],vrsqvsr_tot[ii]);
  //  }
  // }






  ///////////////////////////////////////////////
  //
  // Find horizon's grid position
  // This loop iterates until grid position of horizon doesn't change
  //
  ///////////////////////////////////////////////
  do{

    ////////////////////////////////////
    //
    // get new Kerr parameters from fluxes through inner-radial boundary (includes contribution from self-gravity)
    get_new_metric_parms(&MBH, &a, &QBH, &EP3, &THETAROT);


    // get old and new horizon
    horizontiold=horizoni+horizoncpupos1*N1;
    find_horizon(2); // with new metric parameters, find new wanted horizoni
    horizontinew=horizoni+horizoncpupos1*N1; // wanted value of horizonti, but need to correct for mass enclosed between locations


    // update horizon position in an exactly conservative way
    if(horizontiold!=horizontinew){
      
      // below too may change black hole parameters since horizonti may have moved with new mass
      // horizoni can change due to accreting into black hole or forming black hole due to self-gravity
      // change here is through "effective accretion" when horizon position on grid-scale moves
      // GODMARK: Notice that need Mvsrtot_face1 defined for both self-gravity and if not (i.e. just evolving due to accretion)
      action_ifchange_horizoni(horizontiold,horizontinew);
      
      ////////////////////////////////////
      //
      // must recompute self-gravity to be consistent since may have lost mass to black hole in above operation
      if(DOSELFGRAVVSR){
        // must be called after MHD solution set
        // determine phi(r) , which will be used by set_grid to determine metric
        while(compute_phi_self_gravity_simple(pb)==-1){
        }
      }
    }

  } while(horizontiold!=horizontinew);







  ///////////////////////////////////////////
  //
  // NOW: recompute grid (set_grid()?)
  // 1 tells set_grid that beyond first set_grid call
  // 0 would tell set_grid to set Xmetricold[TT]=t so that no old metric used in connection calculation (i.e. stationarity assumed until iterations are underway)
  // Also do other things that normally need to be done after set_grid()
  //
  ///////////////////////////////////////////
  set_grid(whichtime,CUf, Cunew);

  //  specialstep=2; 

  // user post_set_grid function (Rhor and such things that depend upon metric values)
  init_grid_post_set_grid(pb,GLOBALPOINT(pstagglobal),ucons,GLOBALPOINT(vpotarrayglobal),GLOBALPOINT(Bhatglobal),GLOBALPOINT(panalytic),GLOBALPOINT(pstaganalytic),GLOBALPOINT(vpotanalytic),GLOBALPOINT(Bhatanalytic),GLOBALPOINT(F1),GLOBALPOINT(F2),GLOBALPOINT(F3),GLOBALPOINT(emf));
  // once all parameters are set, now can set dependent items that may be used to set primitives or conservatives
  // assume grid-related parameters same -- just values of metric-related quantities different
  //    post_par_set();


  // unlike compute_new_metric_and_prims(), this function is done once the metric is updated
  // e.g. primitives are computed self-consistently using new metric
  // e.g. dt is set correctly
  // 

  // after all parameters and primitives are set, then can set these items
  // post_init calls store_old_metric()
  //  post_init(); // GODMARK: don't need to do all the things in here, but presently doing them doesn't hurt anything
  // user post_init function
  post_init_specific_init();  // GODMARK: nothing done so far is a problem, ok and good to be done in case user depends on how grid set






  // DEBUG
  //  if(horizoni>0){
  //  gdump(mbhcount++);
  // }



  return(0);

}




// Things to do when horizon position changes due to accretion.  Then absorb those grid cell's mass, etc.
// assume Mvsrface1_tot well-defined globally
// insert desired new horizon due to self-gravity
int action_ifchange_horizoni(int horizontiold,int horizontinew)
{
  int horizonmatches;
  int current_horizonti, test_horizonti;
  int dir,enerregion;
  SFTYPE (*localpcum)[NPR];
  int faketimeorder,fakenumtimeorders;



  if(horizontiold==horizontinew){
    // then no changes necessary
    return(0);
  }
  // otherwise must make changes to black hole


  if(myid==0){
    // prepare to add mass to black hole
    enerregion=OUTSIDEHORIZONENERREGION;
    localpcum=pcumreg[enerregion];
  }


  current_horizonti = horizontinew; // try desired horizon position
  horizonmatches=0;
  do{

   

    if(myid==0){

      // USE FACE1 since treating like flux @ FACE1 and then conservation is explicit
      
      // gradual accumulation of accreted mass leads to Rhor increasing and horizontinew increasing
      // Once horizontinew changes, now flux surface for horizon changes and effectively suddenly we added mass of a number of cells to black hole.


#if(DOSELFGRAVVSR)
      // Since self-gravitating part has non-BH mass, this is what's added up to horizontinew
      
      // using new horizon position to accumulate all the mass within a cell
      //      horizontinew = horizoni + horizoncpupos1*N1;
      
      // pretend added mass are in flux form, so get signs right
      // doesn't matter which cpu we add this to since just cumulates into total
      // Mfactor divisor because localpcum is in normal units not metric units (later converted back if computing MBH)
      localpcum[X1DN][RHO] += -Mvsrface1_tot[current_horizonti]/Mfactor; // not used, but keep track
      if(REMOVERESTMASSFROMUU!=2) localpcum[X1DN][UU] += Mvsrface1_tot[current_horizonti]/Mfactor;
      localpcum[X1DN][U3] += -Jvsrface1_tot[current_horizonti]/Jfactor;
#endif
    }
    
    // set metric parameters for all CPUs
    get_new_metric_parms(&MBH, &a, &QBH, &EP3, &THETAROT);
    
    // now MBH has changed and need to compute new current_horizonti.
    // possible that current_horizonti will change again because of this extra mass, so need to continue growing until done
    
    // when called by black hole formation:
    // have to recompute since redefined where black hole and self-gravitating regions are and have to recompute self-gravity with this fact
    //    recompute_horizonflux_quantities(2); // 2 tells it to recompute everything and assume horizoni on upside rather than downside
    find_horizon(2);
    test_horizonti = horizoni + horizoncpupos1*N1;


#if(DOSELFGRAVVSR)
    horizonmatches=0;
    // need to make sure horizon matches so account for amount of stuff into black hole
    if(test_horizonti!=current_horizonti){
      trifprintf(" (metric units) MBH=%21.15g\n",MBH);
      if(REMOVERESTMASSFROMUU==2) trifprintf("Tried but failed to use: (metric units) Mnew=%21.15g Jnew=%21.15g\n",(-localpcum[X1DN][RHO]+localpcum[X1DN][UU])*Mfactor,localpcum[X1DN][U3]*Jfactor);
      // then remove mass and try again with new larger position
      if(myid==0){
        localpcum[X1DN][RHO] -= -Mvsrface1_tot[current_horizonti]/Mfactor;
        if(REMOVERESTMASSFROMUU!=2) localpcum[X1DN][UU] -= Mvsrface1_tot[current_horizonti]/Mfactor;
        localpcum[X1DN][U3] -= -Jvsrface1_tot[current_horizonti]/Jfactor;
      }
      trifprintf("Trying to move horizon from %d to %d\n",horizontiold, horizontinew);
      trifprintf("Tried horizon: %d but led to real horizon %d, so trying %d\n",current_horizonti,test_horizonti,current_horizonti+1);
      if(REMOVERESTMASSFROMUU==2) trifprintf("Mold=%21.15g Jold=%21.15g\n",(-localpcum[X1DN][RHO]+localpcum[X1DN][UU])*Mfactor,localpcum[X1DN][U3]*Jfactor);
      current_horizonti++;
      if(current_horizonti>=ncpux1*N1){
        dualfprintf(fail_file,"Never found consistent horizon\n");
        myexit(7236);
      }
      horizonmatches=0;
    }
    else{
      horizonmatches=1;
    }
#else
    horizonmatches=1;
    // GODMARK: for now:
    // when not doing self-gravity, don't worry about mass-energy of shells of matter subsumbed by moving of horizoni
    // when horizon moves due to accretion, just move horizon -- there is no inconsistency since no extra mass added from grid
#endif


  } while(horizonmatches==0);



  // put here anything when changing horizonti (apart from recomputing metric, which is done outside)
  // if got here then horizonti changed, so need new positions of fluxes
  faketimeorder=-1;
  fakenumtimeorders=-1;
  recompute_fluxpositions(0,faketimeorder,fakenumtimeorders,nstep,t);




  // recompute_horizonflux_quantities() will get new horizoni and horizoncpupos1 for rest of code
  // localpcum will be updated correctly
  // get_new_metric_parms updates metric parameters for black hole
  return(-1); // -1 signifies that made some changes so may have to recompute metric or something


}





// stores metric into memory
// use by DOEVOLVEMETRIC==1 in order to get old metric so can compute time derivative of metric
// notice we store metric in PRIMECOORDS so assumes in some weak sense that the metric changes while the coordinate labels due to dxdxp remain the same (although probably error is comparable so ok if dxdxp changes)
// usually comes AFTER set_grid() since set_grid needs old metric and just newly computed metric in order to compute time-derivatives
// called through compute_new_metric_and_prims() via post_init after new metric is computed so next time metric is computed the old metric can be used with the new metric to determine the connection coefficients that have dgdt in them.
// the post_init function is in initbase.c
// function can be accessed outside this file
int store_old_metric(void)
{
  int i,j,k;
  int jj,kk;
  int loc;


#if(DOEVOLVEMETRIC==0)
  return(0); // shouldn't be here if DOEVOLVEMETRIC==1
#endif

  ///////////////////////////////
  //
  // store old metric so can compute time derivatives of metric


  // over grid locations needing these quantities
  for (loc = NPG - 1; loc >= 0; loc--) {
    // as in set_grid()
    COMPFULLLOOPP1{

#if(NEWMETRICSTORAGE)
      // copy entire structure at once
      GLOBALMETMACP1A0(compgeomlast,loc,i,j,k)=GLOBALMETMACP1A0(compgeom,loc,i,j,k);
#else
      DLOOP(jj,kk) GLOBALMETMACP1A2(gcovlast,loc,i,j,k,jj,kk)=GLOBALMETMACP1A2(gcov,loc,i,j,k,jj,kk);
      DLOOPA(jj) GLOBALMETMACP1A1(gcovpertlast,loc,i,j,k,jj)=GLOBALMETMACP1A1(gcovpert,loc,i,j,k,jj);
      GLOBALMETMACP1A0(alphalapselast,loc,i,j,k)=GLOBALMETMACP1A0(alphalapse,loc,i,j,k);
#endif
      
    }
  }


  return(0);


}




// initialize self-gravity terms
int init_selfgrav(void)
{
  int ii;

  GRAVLOOP(ii){

    dMvsr_tot[ii]=0;
    dVvsr_tot[ii]=0;
    vrsqvsr_tot[ii]=0;
    dTrrvsr_tot[ii]=0;
    Mvsr_tot[ii]=0;
    Mvsrface1_tot[ii]=0;
    MOrvsr_tot[ii]=0;
    phivsr_tot[ii]=0;
    dJvsr_tot[ii]=0;
    Jvsr_tot[ii]=0;
    Jvsrface1_tot[ii]=0;

    // position data
    rcent_tot[ii]=0;
  }

  return(0);
}




// compute jacobian used for self-gravity
void compute_jacvol(int i, int j, int k, FTYPE *jacvol)
{
  FTYPE dr,dh,dp;
  FTYPE jacvolnew;
  void compute_dr(int i, int j, int k, FTYPE *dr);
  void compute_dh(int i, int j, int k, FTYPE *dh);
  void compute_dp(int i, int j, int k, FTYPE *dp);



  compute_dr(i,j,k,&dr);
  compute_dh(i,j,k,&dh);
  compute_dp(i,j,k,&dp);

  jacvolnew = dr*dh*dp;
  *jacvol=jacvolnew;

  //  dualfprintf(fail_file,"dr=%21.15g dh=%21.15g dp=%21.15g\n",dr,dh,dp);


}




// compute volume version of d\phi 
void compute_dp(int i, int j, int k, FTYPE *dp)
{
  FTYPE Xkm[NDIM],Vkm[NDIM],Xkp1[NDIM],Vkp1[NDIM];
  

  //  coord_ijk(i, j, k, FACE3, Xkm);
  //  bl_coord_ijk(i, j, k, FACE3, Vkm);
  bl_coord_coord(i, j, k, FACE3, Xkm, Vkm);

  bl_coord_coord(i, j, k+1, FACE3, Xkp1, Vkp1); // must use bl_coord_coord() not _ijk versions
  


#if(N3>1)
  *dp = fabs(Vkp1[3]-Vkm[3]); // just d\phi 
#else
  *dp = 2.0*M_PI;
#endif

}



// compute volume version of d\theta
void compute_dh(int i, int j, int k, FTYPE *dh)
{
  FTYPE Xjm[NDIM],Vjm[NDIM],Xjp1[NDIM],Vjp1[NDIM];

  //  coord_ijk(i, j, k, FACE2, Xjm);
  //  bl_coord_ijk(i, j, k, FACE2, Vjm);
  bl_coord_coord(i, j, k, FACE2, Xjm, Vjm);

  bl_coord_coord(i, j+1, k, FACE2, Xjp1, Vjp1); // must use bl_coord_coord() not _ijk versions


#if(N2>1)
  // dh doesn't reduce to Pi
  *dh = fabs(-(cos(Vjp1[2])-cos(Vjm[2]))); // really \delta(-cos(h)) = sinh dh
  // below is inconsistent with rest of code
  //    dh = (N2==1) ? M_PI : -(cos(Vjp1[2])-cos(Vjm[2])); // really \delta(-cos(h)) = sinh dh
#else
  *dh = 2.0; // full integral reduces to this
#endif


}



// compute volume version of dr
void compute_dr(int i, int j, int k, FTYPE *dr)
{
  FTYPE Xim[NDIM],Vim[NDIM],Xip1[NDIM],Vip1[NDIM];


  // FINITE VOLUME -- dh is too different in limit of totalsize[2]==1
  //  coord_ijk(i, j, k, FACE1, Xim);
  //  bl_coord_ijk(i, j, k, FACE1, Vim);
  bl_coord_coord(i, j, k, FACE1, Xim, Vim);

  bl_coord_coord(i+1, j, k, FACE1, Xip1, Vip1); // must use bl_coord_coord() not _ijk versions

  
#if(N1>1)
  *dr = fabs(THIRD*(pow(Vip1[1],3)-pow(Vim[1],3))); // really \delta(r^3/3) = r^2 dr
#else
  *dr = 1.0; // shouldn't matter
#endif



}















// whether to use unewglobal[] in computing M(r)
// ==0 means compute stress-energy tensor components necessary
// ==1 means use conserved energy
// ==2 means use conserved particle number
// ==3 use conserved particle number without gravitational term in \detg=\alpha\sqrt{\gamma}
// can use unewglobal or primitive for substep or skipping step method since unewglobal is kept up-to-date in advance.c relative to where metric is computed
//
// check consistency with USETRUEENERGYACCRETED
#define USEUNEW 0

// whether to use gdetvol_func and gdetvol[] in integrating self-gravity potential
#define SELFGRAVVOLDIFF 1

// whether to analytically integrate constant density part of density
// can only be used with non-relativistic gravity where density comes in linearly and not non-linearly
// once doing trapezium rule, USERHOREF is not as necessary but still helps near r=0
#define USERHOREF 0
// above NOT setup


// if USEUNEW==0,then do one of the below
//
// whether to use T^t_t as explicitly defined or just density as static case gives
// extra factor of \gamma^2 means mass enclosed diverges when "black hole" forms
// but how then should it match onto mass of black hole?
// seems more reasonable to use density instead of -T^t_t since then space-times more closely match
// energy of motion goes where?
// ==1 means use -T^t_t
// ==2 means use rho
// ==3 means use rho u^t (like USEUNEW==2/3 but for substeps) ( *** somehow conserved baryon number becomes not conserved even with new advance() method *** )
// == 4 means using -T^t_t/(-u_t) and T^t_\phi/(-u_t) with jacvol
// == 5 use u.T.u as energy and u.T.P_\phi as angular momentum
// == 6 use Komar mass
// == 7 like 5 but with static frame at infinity
// == 8 like 5 but with stationary ZAMO frame
// == 9 Shibata 1997 -LIKE form for M*
// == 10 extention of #9 for all energies (seems to be most correct for static and most general for non-static and with all energies)
#define USEMHDTTTT 10

// When \Gamma<0 or \Gamma>Gammagrav_max, then assume black hole has formed
#define Gammagrav_max (100.0)

// whether to use vrsq!=0
// comes self-consistently into metric and rest is correct given that metric
// doesn't seem to work
#define COMPUTE_VRSQ 0

// activated by DOSELFGRAVVSR
// Compute self gravitational potential assuming only radial dependence
// only correct when r(x_1) and not r(x_1,x_2), otherwise have to make gravity 2-D
// called by compute_new_metric_and_prims() at start of function to use conserved quantities/primitives to define the self-gravity potential
static int compute_phi_self_gravity_simple(FTYPE (*pb)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int ii,jj,kk;
  int pl,pliter;
  FTYPE dr,dh,dp,jacvol,jacvolnew,jacvolold;
  FTYPE rface1,dphi,dphiim1;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q;
  FTYPE mhdnative[NDIM][NDIM];
  FTYPE mhdprime[NDIM][NDIM];
  FTYPE V[NDIM], X[NDIM];
  FTYPE uconnative[NDIM],ucovnative[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE Gamma;
  FTYPE rhoref,rhoref_send;
  int get_rhoref(int i, int j, int k, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE *rhoref);
  FTYPE vrsq;
  static int firsttime=1;
  static FTYPE truephiouter;
  FTYPE origphiouter;
  static FTYPE trueMouter;
  FTYPE origMouter;
  int ri,startphiiter;
  FTYPE Mrealouter_tot,phirealouter_tot;
  static int formedblackhole=0,whereformedblackhole=-100,didformblackhole=0;
  FTYPE alpha;
  FTYPE Jrealouter_tot,trueJouter,origJouter;
  extern int recompute_horizonflux_quantities(int fromwhere);
  FTYPE Gammatot,MBHOr;
  int horizonmatches;
  FTYPE phibh;
  FTYPE phibh_compute(FTYPE M, FTYPE a, FTYPE r, FTYPE th);
  int horizontiold,horizontinew;
  int action_ifchange_horizoni(int horizontiold, int horizontinew);
  int dir,enerregion;
  FTYPE udotudotT, udotPdotT;
  FTYPE Ttr, etad[NDIM],xiu[NDIM], Tetaxi, Tother;
  void compute_jacvol(int i, int j, int k, FTYPE *jacvol);
  FTYPE dvsq;
  FTYPE uu0ud0;
  FTYPE rho,u;
  FTYPE dM,dJ,dTrr;
  FTYPE bsq,P,Ptot;
  FTYPE rhoco,rhouu0eff,Jco,Jeff;
  FTYPE unewuueff;
  FTYPE absrcent;
  extern int find_RinRout(FTYPE *localRin, FTYPE *localRout);
  static FTYPE localRin, localRout;
  // LOOPBOUND stuff
  int inboundloop[NDIM];
  int outboundloop[NDIM];
  int innormalloop[NDIM];
  int outnormalloop[NDIM];
  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  int riin,riout,rjin,rjout,rkin,rkout;
  int dosetbc[COMPDIM*2];
  int boundvartype=BOUNDPRIMTYPE;
  int loc;




  ////////////////////////
  //
  // If not doing self-gravity, then shouldn't be here, so leave.
  //
  ///////////////////////

#if(!DOSELFGRAVVSR)
  return(0);
#endif





  ////////////////////////
  //
  // set bound loop (used to set boundary condition values for Mvsr_tot, etc.)
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);








  ////////////////////////
  //
  // Initialize self-gravity quantities
  //
  ///////////////////////

  // non-relativistic
  // M(x_1) = \int_0^{x_1} \rho \detg dx1 dx2 dx3

  // more relativistic and consistent with mass-energy accreted into black hole
  // lab-frame volume energy
  // M(x_1) = \int_0^{x_1} T^t_t \detg dx1 dx2 dx3


  // force = -GM(x_1)/r^2

  // Potential = \phi = -\int_0^{x_1} f dr =  -\int_0^{x_1} f (dx1/dr) dr
  // above assumes dr/dx2=0, which is also assumption that gravity can be represented by 1D function

  // default is no black hole formed
  if(firsttime){
    // presume form black hole only once and value stays 1 once set to 1
    // then growth of black hole only proceeds through accretion
    didformblackhole=0;
    formedblackhole=0;
    whereformedblackhole=-100;

    // get Rout without dependence on coord.c parameters
    find_RinRout(&localRin, &localRout);
  }
  else{
    formedblackhole=0; // assume no newly formed black hole
  }

  rhoref=0.0; // not using this method right now
  //////////////////////////////////////////////////////////////////
  //
  // per CPU -- determine the per shell mass as a function of radius
  // initialize dMvsr
  // notice that this sets dMvsr outside CZLOOP, so 0 contribution there
  GRAVLOOP(ii){
    dMvsr[ii]=0; // used to determine Mvsr
    dJvsr[ii]=0; // used to determine a (black hole spin) when forming black hole
    dTrrvsr[ii]=0; // used to determine \phi (grav. pot.) only
    dVvsr[ii]=0; // used to determine vrsq

    // only thing directly used for metric
    vrsqvsr[ii]=0; // used to determine vrsq used in \phi and metric
  }








  /////////////////////////////////////////////////////////
  //
  // Compute vrsqvsr
  //
  /////////////////////////////////////////////////////////



  ///////////////////////
  //
  // check values
  //
  ///////////////////////
#if(PRODUCTION==0)
  if(debugfail>=2){
    enerregion=OUTSIDEHORIZONENERREGION;
    enerpos=enerposreg[enerregion];
    
    ZSLOOP(enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]){
      PALLLOOP(pl){
        if(!isfinite(MACP0A1(pb,i,j,k,pl)) || !isfinite(GLOBALMACP0A1(unewglobal,i,j,k,pl))){
          dualfprintf(fail_file,"pl=%d pb=%21.15g ucumformetric=%21.15g\n",pl,MACP0A1(pb,i,j,k,pl),GLOBALMACP0A1(unewglobal,i,j,k,pl));
        }
      }
    }
  }
#endif



  
  
#if(COMPUTE_VRSQ)
  // for self-gravity, only loop over region outside horizon
  enerregion=OUTSIDEHORIZONENERREGION;
  enerpos=enerposreg[enerregion];
  
  loc=CENT;
  // loop over active region outside horizon
  // // can't use COMPZLOOP with enerpos: (COMPZSLOOP(enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]){)
  // SECTIONMARK: This goes beyond ACTIVEREGION
  ZSLOOP(enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]){


    // get jacvol
    compute_jacvol(i, j, k, &jacvol);

    // determine mass-weighted vrsq
    get_geometry(i,j,k,loc,ptrgeom);
    MYFUN(get_state(MAC(pb,i,j,k), ptrgeom, &q) ,"phys.c:source()", "get_state() dir=0", 1);

    // I think the way v^2 is introduced in prior papers is wrong -- should really be a dr dt metric term as in Lemos, J. (1992, Phys. Rev. Letters)
    vrsq = fabs((q.ucon[1])*(q.ucov[1])/((q.ucov[0])*(q.ucon[0])));
    // density-weighted square of radial 3-velocity

    rho = GLOBALMACP0A1(pb,i,j,k,RHO);

    vrsqvsr[startpos[1] + i] += rho*vrsq*jacvol;

    dVvsr[startpos[1]+i] += rho*jacvol;

  }


  // now cumulate over all CPUs and give result to CPU=0 to do simple integrals
  enerregion=OUTSIDEHORIZONENERREGION;
  // enerregion not really used for this integration, so assume per-CPU loops used enerregion information
  // diagnostics do not do this integrate() so they may be a substep out of "synch" with other diagnostics that are fully up to date
  // old way:
  //  GRAVLOOP(ii){
  //    if(integrate(&dVvsr[ii],&dVvsr_tot[ii],CUMULATIVETYPE,enerregion)>=1) return(1);
  //    if(integrate(&vrsqvsr[ii],&vrsqvsr_tot[ii],CUMULATIVETYPE,enerregion)>=1) return(1);
  //  }

  // new way:
  if(integrate(NUMGRAVPOS,&dVvsr[-N1BND],&dVvsr_tot[-N1BND],CUMULATIVETYPE,enerregion)>=1) return(1);
  if(integrate(NUMGRAVPOS,&vrsqvsr[-N1BND],&vrsqvsr_tot[-N1BND],CUMULATIVETYPE,enerregion)>=1) return(1);


  //////////////////////////////////////////////////////////////////
  //
  // Now compute vrsq(r)
  if(myid==0){

    /////////////////////////////
    //
    // Compute actual local vrsq and set boundary conditions
    //
    ////////////////////////////
    GRAVLOOP(ii){
      // obtain volume-averaged vrsq
      // SMALL assumes vrsqvsr_tot is 0 when dVvsr_tot is 0 simply because out of integration range
      vrsqvsr_tot[ii]/=(dVvsr_tot[ii]+SMALL); // \int vrsq jac / \int jac

    }

    //////////////////
    // now set boundary conditions on vrsqvsr_tot
    if(BCtype[X1DN]==R0SING){
      startphiiter=0;
      ri=riin;
      LOOPBOUND1IN{ 
        vrsqvsr_tot[i]=vrsqvsr_tot[ri+(ri-i-1)];
      }
    }
    else{
      startphiiter=-N1BND;
      // then assume outflow, which for vrsqvsr means copy from active region
      LOOPBOUND1IN{
        vrsqvsr_tot[i]=vrsqvsr_tot[0];
      }
    }
  }

  // now broadcast total vrsqvsr_tot to all CPUs
#if(USEMPI)
  MPI_Bcast(&(vrsqvsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
#else
  // then already have full phi in right memory space
#endif


  // now have vrsqvsr_tot on all CPUs
#else
  GRAVLOOP(ii) vrsqvsr_tot[ii] = 0.0;
#endif













  /////////////////////////////////////////////////////////
  //
  // Compute dMvsr, dTrrvsr
  //
  /////////////////////////////////////////////////////////


  // for self-gravity, only loop over region outside horizon
  enerregion=OUTSIDEHORIZONENERREGION;
  enerpos=enerposreg[enerregion];

  // find dMvsr in each shell
  // this is the only loop in this function that is over more than just radius (i)
  //  COMPZLOOP{    // consistent with GRAVLOOPACTIVE(ii) used below for RHOREF
  
  //  COMPZSLOOP(-N1BND,N1+N1BND-1,0,N2-1,0,N3-1){ // over entire domain (was needed when no boundary conditions on Mvsr and phivsr)
  

  // loop over active region outside horizon
  //COMPZSLOOP(enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]){
  //   // SECTIONMARK: Note that this integrates outside ACTIVEREGION
  loc=CENT;



  ZSLOOP(enerpos[X1DN],enerpos[X1UP],enerpos[X2DN],enerpos[X2UP],enerpos[X3DN],enerpos[X3UP]){

    // avoid measurements outside valid CPU regions so don't end up cumulating multiple zones -- want only unique zones when doing MPI
    //    if(i<0 && startpos[1]+i>=0 || i>N1-1 && startpos[1]+i<=totalsize[1]-1) break; // was used when looping over entire domain
 
    // correct only when purely radial dependence
    // assumes unewglobal updated for both DOENOFLUX==ENOFINITEVOLUME and (DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFLUXSPLIT)
    // unewglobal includes \detg factor
    // flipped sign of both UU and RHO terms since energy comes in with negative sign and we want M(r) to be positive


#if(USEUNEW==1)
    get_geometry(i,j,k,loc,ptrgeom);
    //    alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(TT,TT)]);
    alpha = ptrgeom->alphalapse;

    // when EVOLVEMETRICSUBSTEP==0, then unewglobal is final conserved quantity from RK'ing
    // when EVOLVEMETRICSUBSTEP==1, then unewglobal is either final or so-far accumulated conserved quantity.  Not using "Uf" from evolution since Uf represents conserved state at some intermediate state that isn't necessarily going to accumulate to get final RK value.
    // However, non-linearity means that using unewglobal each time is more correct with RK method so should converge to fourth order in correct way
    dM = -GLOBALMACP0A1(unewglobal,i,j,k,UU)*dVF/(alpha*(q.ucov[TT]));

    if(REMOVERESTMASSFROMUU==2){
      // then need to add back in the rest-mass
      dM += GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF/(alpha);
    }
    else{
      // then already correct
    }

    dJ = GLOBALMACP0A1(unewglobal,i,j,k,U3)*dVF/(alpha*(q.ucov[TT]));

#elif(USEUNEW==2)
    // when EVOLVEMETRICSUBSTEP==0, then unewglobal is final conserved quantity from RK'ing
    // when EVOLVEMETRICSUBSTEP==1, then unewglobal is either final or so-far accumulated conserved quantity.  Not using "Uf" from evolution since Uf represents conserved state at some intermediate state that isn't necessarily going to accumulate to get final RK value.
    // However, non-linearity means that using unewglobal each time is more correct with RK method so should converge to fourth order in correct way

    dM = GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF;
    // includes gravitational mass, which overestimates correct mass

    dJ = GLOBALMACP0A1(unewglobal,i,j,k,U3)*dVF;


#elif(USEUNEW==3)
    // try to remove gravitational mass such that baryon number is conserved
    get_geometry(i,j,k,loc,ptrgeom);
    //    alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(TT,TT)]);
    alpha = ptrgeom->alphalapse;
    // since \detg = \alpha \sqrt(\gamma) and only want \sqrt(\gamma) part
    dM = GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF/(alpha);

    dJ = GLOBALMACP0A1(unewglobal,i,j,k,U3)*dVF/(alpha*(q.ucov[TT]));




#elif(USEUNEW==0)
    

    get_geometry(i,j,k,loc,ptrgeom);
    coord_ijk(i, j, k, loc, X);
    bl_coord_ijk(i, j, k, loc, V);
    dxdxprim_ijk(i, j, k, loc, dxdxp);
    idxdxprim_ijk(i, j, k, loc, idxdxp);

    // get MHD stress-tensor in PRIMECOORDS
    MYFUN(get_state(MAC(pb,i,j,k), ptrgeom, &q) ,"phys.c:source()", "get_state() dir=0", 1);
    DLOOPA(jj) mhd_calc_0(MAC(pb,i,j,k), jj, ptrgeom, &q, mhdprime[jj], NULL); // includes rest-mass

    // now need to convert from PRIMECOORDS to MCOORD (for use natively with spherical polar coorindates to determine metric)
    DLOOP(jj,kk) mhdnative[jj][kk] = mhdprime[jj][kk];
    metptomet_simple_Tud(dxdxp,idxdxp,mhdnative);

    // get u^\mu in native form
    DLOOPA(jj) uconnative[jj] = q.ucon[jj];
    metptomet_simple(dxdxp,uconnative);

    // get u_\mu in native form
    DLOOPA(jj) ucovnative[jj] = q.ucov[jj];
    metptomet_ucov_simple(idxdxp,ucovnative);

    // thermo quantities
    rho = MACP0A1(pb,i,j,k,RHO);
    u = MACP0A1(pb,i,j,k,UU);

    // now do integral for M(r)
    //
    // M(r) = -\int_0^r [4\pi s^2 T^t_t(s) - vrsq/2] ds
    //
    // M(r) = -\int_0^r [s^2 sin(\theta) T^t_t ds d\theta d\phi - vrsq/2 ds]
    //
    // Standard r(x_1) and \phi independent of other coordinates
    // M(r) = -\int_0^r s^2 sin(\theta) T^t_t dx1 (ds/dx1) (dx2 d\theta/dx2 + dx1 d\theta/dx1) dx3 (d\phi/dx3)
    //
    // Alternative is to form gdet from spherical polar Minkowski metric in PRIMECOORDS (i.e. mapped)

#if(SELFGRAVVOLDIFF)
    compute_jacvol(i,j,k,&jacvol);
#else
    dr = fabs(dx[1]*dxdxp[1][1]);
    dh = (totalsize[2]==1) ? 2.0 : fabs(dx[2]*dxdxp[2][2]+dx[1]*dxdxp[2][1]);
    //dh = (dx[2]*dxdxp[2][2]);
    dp = fabs(dx[3]*dxdxp[3][3]);
    jacvolold = fabs(V[1]*V[1]*sin(V[2])*dr*dh*dp);
    jacvol=jacvolold; // disabled as above since not used
#endif


    //    dualfprintf(fail_file,"i=%d jacvolold=%21.15g jacvolnew=%21.15g gdet=%21.15g\n",i,jacvolold,jacvolnew,GLOBALMACP1A0(gdet,loc,i,j,k)*dx[1]*dx[2]*dx[3]);

    //    if(startpos[i]+i==0 && j==totalsize[2]/2 && k==0) rhoref=-mhdnative[TT][TT]; // assume get here first for now
    //    if(startpos[i]+i==0 && j==totalsize[2]/2 && k==0) rhoref=0.0; // assume get here first for now

    vrsq = vrsqvsr_tot[startpos[1] + i];
    // problem with dvsq!=0
    dvsq = vrsqvsr_tot[startpos[1] + i] - vrsqvsr_tot[startpos[1]+i-1];
    dvsq = 0.0;

    //    dualfprintf(fail_file,"%d mhdnative[tt][tt]=%21.15g\n",i,mhdnative[TT][TT]);
    // below is generalized version
    // vrsq already dimensionless metric form if c=1
#if(USEMHDTTTT==1)
    // this is formally most naive extension of the static TOV solution
    dM = -(mhdnative[TT][TT]+rhoref)*jacvol+vrsq*0.5*dr;

    // dJ(r)/dr in metric units (most naive)
    dJ = mhdnative[TT][PH]*jacvol;

#elif(USEMHDTTTT==2)
    // below is standard static version so minimizes errors in treating motion in some exaggerated way
    dM = -(-rho+rhoref)*jacvol + vrsq*0.5*dr + 0.5*V[1]*dvsq;


    // dJ(r)/dr in metric units
    // consistent with baryon conservation
    dJ = mhdnative[TT][PH]/(uconnative[TT])*jacvol;

#elif(USEMHDTTTT==3)

    // doesn't include other energies when should

    // conservation of particle number introduced so that final black hole has same mass as Mvsr within some radius
    //    dM = (rho*(uconnative[TT])+rhoref)*jacvol;
    dM = (rho*(uconnative[TT])+rhoref)*jacvol + vrsq*0.5*dr + 0.5*V[1]*dvsq;
    // dJ(r)/dr in metric units
    // consistent with how treating baryon conservation
    dJ = (rho*(uconnative[TT]*(ucovnative[PH])/(-ucovnative[TT])))*jacvol;

#elif(USEMHDTTTT==4)
    // consistent with baryon conservation
    dM = -(mhdnative[TT][TT]/(-ucovnative[TT])+rhoref)*jacvol+vrsq*0.5*dr;

    // dJ(r)/dr in metric units
    // consistent with baryon conservation
    dJ = mhdnative[TT][PH]/(-ucovnative[TT])*jacvol;

#elif(USEMHDTTTT==5)


    // u.T.u
    udotudotT = 0.0;
    DLOOP(jj,kk) udotudotT += mhdprime[jj][kk]*(q.ucov[jj])*(q.ucon[kk]);

    // this is formally 2nd most naive extension of the static TOV solution
    dM = jacvol*udotudotT;

    // u.T.P
    udotPdotT = 0.0;
    kk=PH;
    DLOOPA(jj) udotPdotT += udotudotT*(q.ucov[kk]) + (q.ucov[jj])*mhdprime[jj][kk];


    // dJ(r)/dr in metric units (most naive)
    dJ = jacvol*udotPdotT;

#elif(USEMHDTTTT==6)



    // u.T.u
    Ttr = 0.0;
    DLOOPA(jj) Ttr += mhdprime[jj][jj];

    // coordinate observer eta_\mu
    //    etad[TT] = -1.0;
    // ZAMO frame
    //    alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(TT,TT)]);
    //    etad[TT] = -alpha;
    //etad[RR] = etad[TH] = etad[PH] = 0.0;

    // approximate time-like Killing vector
    xiu[TT] = 1.0;
    xiu[RR] = xiu[TH] = xiu[PH] = 0.0;

    // unit normal
    lower_vec(xiu,ptrgeom,etad);



    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    Tother=0.0;
    DLOOP(jj,kk) Tother += Ttr*delta(jj,kk)*etad[jj]*xiu[kk];

    // Komar version
    dM = jacvol*(2.0*(Tetaxi-0.5*Tother));


    // approximate \phi-like Killing vector
    xiu[TT] = 0.0;
    xiu[RR] = xiu[TH] = 0.0 ;
    xiu[PH] = 1.0;

    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    // dJ(r)/dr in metric units (most naive)
    // Komar version
    dJ = jacvol*(-Tetaxi);

#elif(USEMHDTTTT==7)

    // approximate time-like Killing vector
    xiu[TT] = 1.0;
    xiu[RR] = xiu[TH] = xiu[PH] = 0.0;

    // unit normal
    lower_vec(xiu,ptrgeom,etad);

    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    dM = jacvol*Tetaxi;


    // approximate \phi-like Killing vector
    xiu[TT] = 0.0;
    xiu[RR] = xiu[TH] = 0.0 ;
    xiu[PH] = 1.0;

    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    // dJ(r)/dr in metric units (most naive)
    // Komar version
    dJ = jacvol*(-Tetaxi);

#elif(USEMHDTTTT==8)

    // approximate time-like Killing vector
    xiu[TT] = 1.0;
    xiu[RR] = xiu[TH] = xiu[PH] = 0.0;

    // coordinate observer eta_\mu
    //    etad[TT] = -1.0;
    // ZAMO frame
    //    alpha = 1.0/sqrt(-ptrgeom->gcon[GIND(TT,TT)]);
    alpha = ptrgeom->alphalapse;
    etad[TT] = -alpha;
    etad[RR] = etad[TH] = etad[PH] = 0.0;
    // unit normal

    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    dM = jacvol*Tetaxi;


    // approximate \phi-like Killing vector
    xiu[TT] = 0.0;
    xiu[RR] = xiu[TH] = 0.0 ;
    xiu[PH] = 1.0;

    Tetaxi=0.0;
    DLOOP(jj,kk) Tetaxi += mhdprime[jj][kk]*etad[jj]*xiu[kk];

    // dJ(r)/dr in metric units (most naive)
    // Komar version
    dJ = jacvol*(-Tetaxi);

#elif(USEMHDTTTT==9)

    // -1 if static
    uu0ud0=(uconnative[TT])*(ucovnative[TT]);
    
    
    dM = (
          jacvol*rho*(-1.0/uu0ud0)+
          GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF*(1.0+1.0/uu0ud0)
          );
    
    
    dJ = (
          jacvol*rho*(ucovnative[U3])*(-1.0/uu0ud0)+
          GLOBALMACP0A1(unewglobal,i,j,k,U3)*dVF*(1.0+1.0/uu0ud0)
          );
    

    

    /*
    // equivalent to USEUNEW==2
    dM = (
    GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF
    );


    dJ = (
    GLOBALMACP0A1(unewglobal,i,j,k,U3)*dVF
    );
    */

#elif(USEMHDTTTT==10)


    // -1 if static
    uu0ud0=(uconnative[TT])*(ucovnative[TT]);


    // u.T.u
    udotudotT = 0.0;
    DLOOP(jj,kk) udotudotT += mhdprime[jj][kk]*(q.ucov[jj])*(q.ucon[kk]);
   
    // comoving density (reduced to rho for dust)
    rhoco = udotudotT;



#if(0) // 1st submethod



#if(0)
    // can't use unewglobal[UU] since it had added to it the gravitational changes
    unewuueff = GLOBALMACP0A1(unewglobal,i,j,k,UU)*dVF;
    if(REMOVERESTMASSFROMUU==2){
      // then need to add back in the rest-mass
      unewuueff -= GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF;
    }
    // effective rho u^t from MHD stress-tensor with gdet in it
    rhouu0eff = unewuueff/(q.ucov[TT]);
#else
    unewuueff = mhdprime[TT][TT]*(ptrgeom->gdet)*dVF;
    rhouu0eff = unewuueff/(q.ucov[TT]);
#endif

    // see /mnt/data1/jon/codebackups/crazyformationofbh.zip (and turn on/off DEBUG below -- different results)
    // sometimes (bug!) huge black hole masses (thought bug since inclusiion or removal of DEBUG for jacvol and geomgdV below changed behavior!!!!) GODMARK!!!
    dM = (
          jacvol*rhoco*(-1.0/uu0ud0)+
          rhouu0eff*(1.0+1.0/uu0ud0) // causes mass to be too large (why?)
          // pressure term too large (u+p) u^t + p when trying to form "effective" rho u^t
          );



#elif(0) // 2nd submethod


    // even after mhdnative fixed, still leads to unstable consentrated region
    // mhdnative still has pressure term in it, so probably the problem since include pressure separately
    dM = (
          jacvol*rhoco*(-1.0/uu0ud0)+
          mhdnative[TT][TT]/(ucovnative[TT])*jacvol*(1.0+1.0/uu0ud0)
          );


#elif(0) // 3rd submethod
    // seems to keep things stable.  Recall that pressure included in separate term
    dM = (
          jacvol*rhoco*(-1.0/uu0ud0)+
          rho*(uconnative[TT])*jacvol*(1.0+1.0/uu0ud0)
          );


#elif(0) // 4th submethod
    // see if mhd tensor messed up with dxdxp's
    // stable
    // but when using UNI2LOG with large Afactor the mass is very not conserved and big (~1) BH forms 
    dM = (
          jacvol*rho*(-1.0/uu0ud0)+
          rho*(uconnative[TT])*jacvol*(1.0+1.0/uu0ud0)
          );


#elif(1)  // 5th submethod
    // try using actual limited baryon number since \detg can shrink allowed mass and may limit errors
    dM = (
          jacvol*rho*(-1.0/uu0ud0)+
          // below term is dM~ \rho u^t \detg dx1dx2dx3
          (GLOBALMACP0A1(unewglobal,i,j,k,RHO)*dVF)*(1.0+1.0/uu0ud0)
          );


#elif(0) // 6th submethod
    // doesn't cause mass to be too large, but too small if uu0 large
    dM = (
          jacvol*rhoco
          );
#endif




#if(0) // DEBUG
    if((jacvol-(ptrgeom->gdet)*dVF)/(jacvol+SMALL)>1.2){
      dualfprintf(fail_file,"jacvol=%21.15g geomgdV=%21.15g\n",jacvol,(ptrgeom->gdet)*dVF);
    }
#endif






    // u.T.P (uses u.T.u from above)
    udotPdotT = 0.0;
    kk=PH;
    DLOOPA(jj) udotPdotT += udotudotT*(q.ucov[kk]) + (q.ucov[jj])*mhdprime[jj][kk];


    // comoving angular momentum
    Jco = udotPdotT;

    // effective angular momentum from MHD stress-tensor with gdet in it
    //    Jeff = (mhdprime[TT][PH]*(ptrgeom->gdet))/(ucovnative[TT]);
    //
    // unlike energy, angular momentum is conserved
    // still use mhd instead of unewglobal since unewglobal might eventually have sources
    Jeff = (mhdprime[TT][PH]*(ptrgeom->gdet))*dVF;

    
    dJ = (
          jacvol*Jco*(-1.0/uu0ud0)+
          Jeff*(1.0+1.0/uu0ud0)
          );
    


#endif












    //    dualfprintf(fail_file,"diff=%21.15g mhd=%21.15g rhoref=%21.15g\n",mhdnative[TT][TT]+rhoref,mhdnative[TT][TT],rhoref);

    //    dualfprintf(fail_file,"i=%d jacvol=%21.15g dr=%21.15g dh=%21.15g dp=%21.15g\n",startpos[1]+i,jacvol,dr,dh,dp);

    //    PLOOP(pliter,pl) dualfprintf(fail_file,"pb[%d]=%21.15g\n",pl,MACP0A1(pb,i,j,k,pl));
    //    dualfprintf(fail_file,"dMvsr[%d]=%21.15g jacvol=%21.15g r=%21.15g mhdtt=%21.15g\n",startpos[1]+i,dM,jacvol,V[1],mhdnative[TT][TT]);



    // used itself and in:
    //
    // \Phi(r) = \int_0^r [M(s)/s^2\Gamma ds + 4\pi s T^r_r \Gamma ds]
    //
    // where \Gamma = 1/(1-2M(s)/s)
    //
    
    // differential term for above equation
    // most naive extension of static TOV solution
    // exaggerates energy content
    //    dTrr = mhdnative[RR][RR]*jacvol;


    // slightly less naive perhaps
    P = pressure_rho0_u_simple(i,j,k,loc,rho,u);
    bsq = dot(q.bcon, q.bcov);
    Ptot = P + bsq*0.5 ; 
    dTrr = P*jacvol;



#if(0)
    // DEBUG:
    dualfprintf(fail_file,"i=%d :: rho0=%21.15g u=%21.15g :: dM=%21.15g dTrr=%21.15g :: P=%21.15g jacvol=%21.15g\n",i,rho,u,dM,dTrr,P,jacvol);
#endif


    
    
#endif // end if USENEW==0






    // now apply factores
    dMvsr[startpos[1] + i] += dM*Mfactor; // rho0unittype controls if M is mass unit or energy unit and Mfactor changes in concordance
    dJvsr[startpos[1] + i] += dJ*Jfactor;
    dTrrvsr[startpos[1] + i] += dTrr*Mfactor;



  }// end loop over domain












  //////////////////////////////////
  //  
  // now cumulate over all CPUs and give result to CPU=0 to do simple integrals
  //
  //////////////////////////////////
  enerregion=OUTSIDEHORIZONENERREGION;
  // enerregion not really used for this integration, so assume per-CPU loops used enerregion information
  // diagnostics do not do this integrate() so they may be a substep out of "synch" with other diagnostics that are fully up to date
  // old way:
  //  GRAVLOOP(ii){
  //    if(integrate(&dMvsr[ii],&dMvsr_tot[ii],CUMULATIVETYPE,enerregion)>=1) return(1);
  //    if(integrate(&dJvsr[ii],&dJvsr_tot[ii],CUMULATIVETYPE,enerregion)>=1) return(1);
  //    if(integrate(&dTrrvsr[ii],&dTrrvsr_tot[ii],CUMULATIVETYPE,enerregion)>=1) return(1);
  //  }

  // new way:
  if(integrate(NUMGRAVPOS,&dMvsr[-N1BND],&dMvsr_tot[-N1BND],CUMULATIVETYPE,enerregion)>=1) return(1);
  if(integrate(NUMGRAVPOS,&dJvsr[-N1BND],&dJvsr_tot[-N1BND],CUMULATIVETYPE,enerregion)>=1) return(1);
  if(integrate(NUMGRAVPOS,&dTrrvsr[-N1BND],&dTrrvsr_tot[-N1BND],CUMULATIVETYPE,enerregion)>=1) return(1);
  



#if(0)
  // DEBUG:
  GRAVLOOP(i){
    dualfprintf(fail_file,"i=%d :: dMvsr_tot=%21.15g dJvsr_tot=%21.15g dTrrvsr_tot=%21.15g\n",i,dMvsr_tot[i],dJvsr_tot[i],dTrrvsr_tot[i]);
  }
#endif










  //////////////////////////////////////////////////////////////////
  //
  // Now compute Mvsr/r and Phivsr on CPU myid=0
  //
  //////////////////////////////////////////////////////////////////
  if(myid==0){



    //////////////////////////////
    //
    // Bound dTrrvsr_tot since used directly as computed
    //
    //////////////////////////////

    //////////////////
    // now set boundary conditions on dTrrvsr_tot
    if(BCtype[X1DN]==R0SING){
      ri=0;
      LOOPBOUND1IN{ 
        dTrrvsr_tot[i]=dTrrvsr_tot[ri+(ri-i-1)]; // at CENT
      }
    }
    else{
      // then assume outflow, which for dTrrvsr means no more mass enclosed (i.e. ddTrrvsr=0)
      LOOPBOUND1IN{ 
        dTrrvsr_tot[i]=0.0;
      }
    }





    //////////////////////////////
    //
    // Compute Mvsr_tot from dMvsr_tot
    //
    //////////////////////////////

    // account for fact that dMvsr and Mvsr are both centered and at same location
    // assume true Mvsr_tot starts cumulating at startpos[1]+i==0
    // assumes outer boundary condition is outflow of some kind so nothing special to do
    Mvsr_tot[0] = 0.5*dMvsr_tot[0];
    Mvsrface1_tot[0] =0.0;


    for(ii=1;ii<ncpux1*N1+N1BND+1;ii++){ // GRAVLOOP +1 away from first iteration
      // Trapezium rule
      Mvsr_tot[ii]  = Mvsr_tot[ii-1] + 0.5*dMvsr_tot[ii-1]+0.5*dMvsr_tot[ii];

      // perfect rule since dM is considered to be an integral over cell
      Mvsrface1_tot[ii]  = Mvsrface1_tot[ii-1] + dMvsr_tot[ii-1];
    }
    // determine ON-GRID Mvsr(r) @ outer r
    // this includes all active mass (explicitly includes extra half grid cell at outer edge)
    ii=ncpux1*N1-1;
    Mvsrface1_tot[ii] = Mrealouter_tot=Mvsr_tot[ii]+0.5*dMvsr_tot[ii]; // same as Mvsrface1_tot[ii]+dMvsr_tot[ii-1]


    //////////////////
    // now set boundary conditions on Mvsr_tot
    if(BCtype[X1DN]==R0SING){
      startphiiter=0;
      ri=0;
      LOOPBOUND1IN{ 
        Mvsr_tot[i]=Mvsr_tot[ri+(ri-i-1)]; // at CENT
        Mvsrface1_tot[i]=Mvsrface1_tot[ri+(ri-i)]; // at FACE1
      }
    }
    else{
      startphiiter=-N1BND;
      // then assume outflow, which for Mvsr means no more mass enclosed (i.e. dMvsr=0)
      LOOPBOUND1IN{ 
        Mvsr_tot[i]=0.0;
        Mvsrface1_tot[i]=0.0;
      }
    }





    //////////////////////
    //
    // determine Mvsr/r since want to interpolate that quantity rather than Mvsr directly
    //
    //////////////////////
    GRAVLOOP(ii){
      // the below causes the local mass at ii to be accelerated by itself
      // that's the correct behavior!  Pressure would support it
      MOrvsr_tot[ii] = Mvsr_tot[ii]/rcent_tot[ii];
      // 
      // mass within ii will only be affected by mass inside instead of itself
      // GODMARK
      //      MOrvsr_tot[ii] = Mvsrface1_tot[ii]/rcent_tot[ii];
    }



    //////////////////
    //
    // now set boundary conditions on MOrvsr_tot
    // this is antisymmetric due to division by $r$
    //
    //////////////////
    if(BCtype[X1DN]==R0SING){
      startphiiter=0;
      ri=0;
      LOOPBOUND1IN{ 
        MOrvsr_tot[i]=-MOrvsr_tot[ri+(ri-i-1)]; // at CENT
      }
    }
    else{
      startphiiter=-N1BND;
      // then assume outflow, which for Mvsr means no more mass enclosed (i.e. dMvsr=0)
      LOOPBOUND1IN{ 
        MOrvsr_tot[i]=0.0;
      }
    }


    /////////////////////
    //
    // get outer Mvsr
    //
    /////////////////////
    if(firsttime){
      trueMouter = Mrealouter_tot;
    }
    origMouter = Mrealouter_tot;



    
    //////////////////////////////////////////////
    //
    // now relevel potential so consistent with phiouter
    //
    //////////////////////////////////////////////
    GRAVLOOP(ii){
      // match so that Mvsr at vacuum is constant (sprinkles mass everywhere that's missing from integration)
      // problem: small errors in integration cause mass in inner-radial region to go negative (form BH)
      // new Mvsr set before phivsr set causes huge spike at r=0 and run takes much longer to run 2X longer due to big spike at center
      //      Mvsr_tot[ii] += (trueMouter-origMouter);

      // TRY?
      // M can only grow, not shrink, then shouldn't create negative mass in inner-radial region
      //      Mvsr_tot[ii] += MAX(trueMouter-origMouter,0.0);

#if(0)
      // DEBUG:
      dualfprintf(fail_file,"Mvsr_tot[%d]=%21.15g dMvsr_tot=%21.15g rho=%21.15g\n",ii,Mvsr_tot[ii],dMvsr_tot[ii],MACP0A1(pb,ii,0,0,RHO));
#endif


    }










    //////////////////////////////
    //
    // Compute Jvsr
    //
    //////////////////////////////

    // account for fact that dJvsr and Jvsr are both centered and at same location
    // assume true Jvsr_tot starts cumulating at startpos[1]+i==0
    // assumes outer boundary condition is outflow of some kind so nothing special to do
    Jvsr_tot[0] = 0.5*dJvsr_tot[0];
    Jvsrface1_tot[0] = 0.0;
    for(ii=1;ii<ncpux1*N1+N1BND+1;ii++){ // GRAVLOOP +1 away from first iteration
      // Trapezium rule
      Jvsr_tot[ii]  = Jvsr_tot[ii-1] + 0.5*dJvsr_tot[ii-1]+0.5*dJvsr_tot[ii];

      // perfect rule
      Jvsrface1_tot[ii]  = Jvsrface1_tot[ii-1] + dJvsr_tot[ii-1];
    }
    // determine ON-GRID Jvsr(r) @ outer r
    // this includes all active mass (explicitly includes extra half grid cell at outer edge)
    ii=ncpux1*N1-1;
    Jvsrface1_tot[ii] = Jrealouter_tot=Jvsr_tot[ii]+0.5*dJvsr_tot[ii];


    //////////////////
    // now set boundary conditions on Jvsr_tot
    if(BCtype[X1DN]==R0SING){
      startphiiter=0;
      ri=0;
      LOOPBOUND1IN{ 
        Jvsr_tot[i]=Jvsr_tot[ri+(ri-i-1)];
        Jvsrface1_tot[i]=Jvsr_tot[ri+(ri-i)];
      }
    }
    else{
      startphiiter=-N1BND;
      // then assume outflow, which for Jvsr means no more mass enclosed (i.e. dJvsr=0)
      LOOPBOUND1IN{ 
        Jvsr_tot[i]=0.0;
        Jvsrface1_tot[i]=0.0;
      }
    }


    if(firsttime){
      trueJouter = Jrealouter_tot;
    }
    origJouter = Jrealouter_tot;

    // now relevel potential so consistent with phiouter
    GRAVLOOP(ii){
      // match so that Jvsr at vacuum is constant (sprinkles mass everywhere that's missing from integration)
      // problem: small errors in integration cause mass in inner-radial region to go negative (form BH)
      // new Jvsr set before phivsr set causes huge spike at r=0 and run takes much longer to run 2X longer due to big spike at center
      //      Jvsr_tot[ii] += (trueJouter-origJouter);
      //      dualfprintf(fail_file,"Jvsr_tot[%d]=%21.15g dJvsr=%21.15g rho=%21.15g\n",ii,Jvsr_tot[ii],dJvsr_tot[ii],MACP0A1(pb,ii,0,0,RHO));
    }











    //////////////////////////////
    // now compute phivsr (assumes Mfactor above included correct G dependence)
    // phi is defined everywhere needed to be defined for metric and connection coefficients
    //
    // GODMARK: Apparently might be able to add Mvsr(r) directly to Kerr metric
    //          instead of computing phi first
    //
    // The below is what I was already doing, basically
    // GODMARK: Also, apparently able to compute phi with 2 single integrals over rho
    //
    // i.e. for spherical potential:
    //
    // \phi = G(1/r \int_0^r \rho(s) 4\pi s^2 ds + \int_r^\infty \rho 4\pi s^2 ds)
    //
    //
    //
    // center potential
    // first value

#define GAMMAGRAV(MOr) (1.0/(1.0-2.0*(MOr)+vrsqvsr_tot[ii])) // relativistic TOV factor

    ii=startphiiter;
    absrcent=myfabs(rcent_tot[ii]);
    // if Gamma~1 then can use RHOREF, otherwise cannot since potential is non-linear in density
    Gamma = GAMMAGRAV(MOrvsr_tot[ii]);
    Gammatot = GAMMAGRAV(MBH/absrcent + MOrvsr_tot[ii]);
    //Gamma = 1.0/(1.0-2.0*MvsrOr_tot[ii]); // relativistic TOV factor
    if( didformblackhole==0 && (Gammatot<0.0 || Gammatot>Gammagrav_max) && ii>=1+horizoni+horizoncpupos1*N1){
      formedblackhole=1;
      whereformedblackhole=ii;
    }
    dr=0.5*myfabs(rcent_tot[ii+1]-rcent_tot[ii]);
    dphi=MOrvsr_tot[ii]/(absrcent)*dr + (dTrrvsr_tot[ii] - vrsqvsr_tot[ii]*0.5*dr)/absrcent;
    dphi*=Gamma;
    phivsr_tot[ii] = dphi;
    
    // rest of values
    for(ii=startphiiter+1;ii<ncpux1*N1+N1BND+1;ii++){ // GRAVLOOP +1 away from first iteration

      // see tov_solution_timeindepsol.nb
      Gamma = GAMMAGRAV(MOrvsr_tot[ii]);
      absrcent=myfabs(rcent_tot[ii]);
      Gammatot = GAMMAGRAV(MBH/absrcent + MOrvsr_tot[ii]);
      if( didformblackhole==0 && (Gammatot<0.0 || Gammatot>Gammagrav_max) && ii>=1+horizoni+horizoncpupos1*N1){
        formedblackhole=1;
        whereformedblackhole=ii;
      }
      dr=myfabs(rcent_tot[ii]-rcent_tot[ii-1]);
      dphi=MOrvsr_tot[ii]/(absrcent)*dr + (dTrrvsr_tot[ii] - vrsqvsr_tot[ii]*0.5*dr)/absrcent;
      dphi*=Gamma;
      


      jj=ii-1;
      Gamma = GAMMAGRAV(MOrvsr_tot[jj]);
      absrcent=myfabs(rcent_tot[jj]);
      Gammatot = GAMMAGRAV(MBH/absrcent + MOrvsr_tot[jj]);
      if( didformblackhole==0 && (Gammatot<0.0 || Gammatot>Gammagrav_max)  && jj>=1+horizoni+horizoncpupos1*N1){
        formedblackhole=1;
        whereformedblackhole=jj;
      }
      dr=myfabs(rcent_tot[jj]-rcent_tot[jj-1]);
      dphiim1=MOrvsr_tot[jj]/(absrcent)*dr + (dTrrvsr_tot[jj] - vrsqvsr_tot[jj]*0.5*dr)/absrcent;
      dphiim1*=Gamma;


      //      phivsr_tot[ii] = phivsr_tot[ii-1] + dphi;
      // Trapezium rule
      phivsr_tot[ii] = phivsr_tot[ii-1] +0.5*dphi+0.5*dphiim1;

      //dualfprintf(fail_file,"ii=%d dphi=%g phi=%g ratio=%g\n",ii,dphi,phivsr_tot[ii],phivsr_tot[ii]/(2.0*M_PI/3.0*V[1]*V[1]*rhoref*Mfactor));

      // DEBUG
      //phivsr_tot[ii] = 0.0;
      //phivsr_tot[ii] = 1.0/pow(absrcent,1);


    } // end for look over 0 through ncpux1*N1

    // add last bit for real phi outer that's located at Rout (outer FACE1)
    ii=ncpux1*N1-1;
    Gamma = GAMMAGRAV(MOrvsr_tot[ii]);
    absrcent=myfabs(rcent_tot[ii]);
    Gammatot = GAMMAGRAV(MBH/absrcent + MOrvsr_tot[ii]);
    if( didformblackhole==0 && (Gammatot<0.0 || Gammatot>Gammagrav_max)  && ii>=1+horizoni+horizoncpupos1*N1){
      // new black hole forms only if Gammatot large and position of BH horizon is a grid zone beyond present value
      formedblackhole=1;
      whereformedblackhole=ii;
    }
    dr=0.5*myfabs(rcent_tot[ii+1]-rcent_tot[ii]);
    dphi=MOrvsr_tot[ii]/(absrcent)*dr + (dTrrvsr_tot[ii] - vrsqvsr_tot[ii]*0.5*dr)/absrcent;
    dphi*=Gamma;

    // set total outer phi located at FACE1 at outer radial edge
    phirealouter_tot = phivsr_tot[ii] + dphi;


    //////////////////////////
    //
    // now set boundary conditions on phivsr_tot
    if(BCtype[X1DN]==R0SING){
      ri=0;
      LOOPBOUND1IN{ 
        phivsr_tot[i]=phivsr_tot[ri+(ri-i-1)];
      }
    }
    else{
      // then assume outflow, which for phivsr means no more mass enclosed (i.e. dMvsr=0 so the-so-far-computed phivsr is 0)
      LOOPBOUND1IN{ 
        phivsr_tot[i]=0.0;
      }
    }
  
  




    // now enforce bondary condition that pseudo-relativistic Scharzschild outside (i.e. no extra mass)
    //ii = ncpux1*N1+N1BND; // causes mass at outer ghost cells to matter for Mvsr(r) everywhere
    //
    //    ii = ncpux1*N1-1; // mass enclosed doesn't change on active grid
    //
    //presume that mass enclosed on grid doesn't change
    // otherwise dphi/dt will be non-negligible at arbitrary locations where there is no physical changes
    // suggests one should use a mass-conserving way of integrating
    // but Mvsr-type mass IS NOT conserved!  Only gravitational mass conserved?
    // well, but baryon # is conserved and in static case just integrating to get that so self-consistent
    //
    // so force effective potential of some kind to be fixed at the outer radius, which is consistent with Birkoff's theorem that no monopolar changes can be noticed and exterior space-time remains constant and Schwartzschild
    if(firsttime){
      //trueMouter = Mrealouter_tot;
      //truephiouter = -trueMouter/(r-2.0*Mfactor);

      // Note that black hole + TOV solution can be written as
      // g_{tt} = -\exp{2\phi}, so linear in \phi, so at large radius MBH will introduce it's own potential of MBH/Rout, so remove that from required offset of self-gravitating potential
      // black hole \phi is exactly the below term


      phibh = phibh_compute(MBH, a, localRout,M_PI*0.5);
      // trying to be smart and know answer (ok, matches to expected boundary value)
      truephiouter = phibh_compute(trueMouter,0.0,localRout,M_PI*0.5) + phibh;
      //truephiouter = phirealouter_tot; // should be correct and consistent (no, because no offset here!)
      trifprintf("trueMouter=%21.15g truephiouter=%21.15g router=%21.15g\n",trueMouter,truephiouter,localRout);

      // check whether ever possible to resolve collapse
      // assumes spherical polar coords and myid==0 is on the inner-radial boundary
      if(BCtype[X1DN]==R0SING){
#define GAMMAGRAVTEST(MOr) (1.0/(1.0-2.0*(MOr))) // relativistic TOV factor
        coord_ijk(SHIFT1*1, 0, 0, FACE1, X);
        bl_coord_ijk(SHIFT1*1, 0, 0, FACE1, V);
        Gamma = GAMMAGRAVTEST(trueMouter/V[1]);
        if(Gamma<Gammagrav_max && Gamma>0.0){
          trifprintf("WARNING: Never possible to form black hole at this resolution\n");
        }
        trifprintf("Gamma(max possible) = %21.15g M/r = %21.15g Gammagrav_max=%21.15g\n",Gamma,trueMouter/V[1],Gammagrav_max);
        //
        coord_ijk(SHIFT1*3, 0, 0, FACE1, X);
        bl_coord_ijk(SHIFT1*3, 0, 0, FACE1, V);
        Gamma = GAMMAGRAVTEST(trueMouter/V[1]);
        if(Gamma<Gammagrav_max && Gamma>0.0){
          trifprintf("WARNING: Difficult to form black hole at this resolution\n");
        }
        trifprintf("Gamma(nearly max possible) = %21.15g M/r = %21.15g Gammagrav_max=%21.15g\n",Gamma,trueMouter/V[1],Gammagrav_max);
      }

    }// end if firsttime==1






    origphiouter = phirealouter_tot;
    //    origMouter = Mrealouter_tot;

    //    trifprintf("origMouter=%21.15g origphiouter=%21.15g router=%21.15g\n",origMouter,origphiouter,rcent_tot[ii]);


    // now relevel potential so consistent with phiouter
    phibh = phibh_compute(MBH, a, localRout,M_PI*0.5);
    GRAVLOOP(ii){
      // match onto exterior vacuum (phivsr here is only self-gravitating part, so remove black hole part)
      phivsr_tot[ii] += (truephiouter - phibh - origphiouter);

      // match so that Mvsr at vacuum is constant (sprinkles mass everywhere that's missing from integration)
      // problem: small errors in integration cause mass in inner-radial region to go negative (form BH)
      //      Mvsr_tot[ii] += (trueMouter-origMouter);
      //      dualfprintf(fail_file,"Mvsr_tot[%d]=%21.15g dMvsr=%21.15g rho=%21.15g\n",ii,Mvsr_tot[ii],dMvsr_tot[ii],MACP0A1(pb,ii,0,0,RHO));
    }// end GRAVLOOP
 


  }// end if cpu==0











  ///////////////////////////////////////
  //
  // now broadcast totals to all CPUs
  //
  ///////////////////////////////////////
#if(USEMPI)
  MPI_Bcast(&(Mvsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(Mvsrface1_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(Jvsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(Jvsrface1_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(MOrvsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(phivsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);
  //  MPI_Bcast(&(vrsqvsr_tot[-N1BND]),NUMGRAVPOS,MPI_SFTYPE,MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&formedblackhole,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&whereformedblackhole,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);

#else
  // then already have full phi in right memory space
#endif






#if(0)
  // DEBUG:
  dualfprintf(fail_file,"Jfactor=%21.15g Mfactor=%21.15g\n",Jfactor,Mfactor);

  dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",truephiouter,phibh,origphiouter,trueJouter,origJouter,trueMouter,origMouter);


  GRAVLOOP(ii){
    dualfprintf(fail_file,"ii=%d Mvsr_tot=%21.15g  Mvsrface1_tot=%21.15g  MOrvsr_tot=%21.15g  phivsr_tot=%21.15g dTrrvsr_tot=%21.15g vrsqvsr_tot=%21.15g Jvsr_tot=%21.15g Jvsrface1_tot=%21.15g\n",ii,Mvsr_tot[ii],Mvsrface1_tot[ii],MOrvsr_tot[ii],phivsr_tot[ii],dTrrvsr_tot[ii],vrsqvsr_tot[ii],Jvsr_tot,Jvsrface1_tot);
  }
  myexit(0);
#endif










    
  if(formedblackhole){
    // presently this is setup that this only happens once.

    // DEBUG GODMARK
    //    DODIAGEVERYSUBSTEP=1;
    //    DTd=1E-5; // get every dump


#if(0)
    if(didformblackhole==0){
      GRAVLOOP(ii) dualfprintf(fail_file,"r=%21.15g Mvsr=%21.15g Mvsrface1=%21.15g\n",rcent_tot[ii],Mvsr_tot[ii],Mvsrface1_tot[ii]);
    }
#endif


    Gamma = 1.0/(1.0-2.0*MOrvsr_tot[whereformedblackhole]+vrsqvsr_tot[whereformedblackhole]); // relativistic TOV factor
    Gammatot = 1.0/(1.0-2.0*(MBH/myfabs(rcent_tot[whereformedblackhole]) + MOrvsr_tot[whereformedblackhole])+vrsqvsr_tot[whereformedblackhole]); // relativistic TOV factor

    trifprintf("Formed black hole ( %d / %d ) at whereformedblackhole=%d r=%21.15g Mcent(r)=%21.15g Mface1(r)=%21.15g M/r(r)=%21.15g Gamma=%21.15g Gammatot=%21.15g\n",formedblackhole,didformblackhole,whereformedblackhole,rcent_tot[whereformedblackhole],Mvsr_tot[whereformedblackhole],Mvsrface1_tot[whereformedblackhole],MOrvsr_tot[whereformedblackhole],Gamma,Gammatot);


    if(MCOORD==KS_BH_TOV_COORDS){

      // whereformedblackhole+1 is used because including all mass in that cell that generated horizon 
      horizontinew=whereformedblackhole++ ; // expand horizon a bit to be at FACE1 -- same location as horizoni
      // if horizon moves, make changes
      horizontiold=horizoni+horizoncpupos1*N1;
      action_ifchange_horizoni(horizontiold,horizontinew);


      // below should never happen since above forces counter inequality
      //      if(whereformedblackhole>horizoni+horizoncpupos1*N1){
      // dualfprintf(fail_file,"black hole shifted to a lower radius compared to estimate - fluxes will be wrong and can change evolution and conservation\n");
      //      }


      trifprintf("whereformedblackhole=%d horizoncpupos1=%d horizoni=%d\n",whereformedblackhole,horizoncpupos1,horizoni);

      trifprintf("Using KS BH+TOV mixed metric and turning on BH with mass MBH=%21.15g, spin=%21.15g, and charge=%21.15g EP3=%21.15g THETAROT=%21.15g in metric units\n",MBH,a,QBH,EP3,THETAROT);

      // DEBUG:
      //      myexit(9826);


    }
    else if(MCOORD==KS_TOV_COORDS){
      // then must stop
      trifprintf("Using TOV self-gravity and once formed BH can't continue\n");
      myexit(987);
    }
    else if(MCOORD==BLCOORDS ||
            MCOORD==KSCOORDS ||
            MCOORD==KS_JP1_COORDS ||
            MCOORD==HTMETRIC ||
            MCOORD==HTMETRICACCURATE ||
            MCOORD==SPCMINKMETRIC){
      // then using perturbed self-gravity and shouldn't form black hole!
      trifprintf("Using perturbed self-gravity and shouldn't form BH\n");
      myexit(988);      
    }
     
  }






  ///////////////////////////
  //
  // things to do before exiting function with good status
  // firsttime is used by ALL CPUs
  firsttime=0;







  // DEBUG
  //  if(didformblackhole){
  //  if(nstep<1100){
  //    GRAVLOOP(ii) dualfprintf(fail_file,"ii=%d rcent=%21.15g phivsr=%21.15g Mvsr=%21.15g Jvsr=%21.15g vrsqvsr=%21.15g\n",ii,rcent_tot[ii],phivsr_tot[ii],Mvsr_tot[ii],Jvsr_tot[ii],vrsqvsr_tot[ii]);
  //  }
  // }






  /////////////////////////
  //
  // determine how to exit
  //
  /////////////////////////
  if(formedblackhole==0 || didformblackhole)  return(0); // then normal self-gravity
  else if(formedblackhole==1 && didformblackhole==0){
    didformblackhole=1;
    // then self-gravitating fluid formed black hole and need to recompute self-gravity with location of horizon known
    // presumes that space-time is ok-modelled by self-gravitating fluid as proxy for black hole while computing new self-gravity that cuts at horizon
    // if don't recompute self-gravity, then can be severely off in metric and next metric update will lead to large temporal changes in the metric causing the connection coefficients to be large
    return(-1); 
  }
  else{
    dualfprintf(fail_file,"Unknown formedblackhole==%d\n",formedblackhole);
    myexit(9827);
    return(1);
  }




  return(1); // shouldn't get here

}









FTYPE mysign(FTYPE x)
{
  if(x<0) return(-1.0);
  else return(1.0); // unlike sign, if positive or 0 then return positive 1.0
}

FTYPE myfabs(FTYPE x)
{
  // return(fabs(x));
  return(x);

}


// g_{tt} = -exp(2\phi) for Kerr black hole in BL or KS coords
FTYPE phibh_compute(FTYPE Ms, FTYPE as, FTYPE r, FTYPE th)
{
  FTYPE phibh;
  FTYPE r2,a2,ch,ch2;

  r2=r*r;
  a2=as*as;
  ch=cos(th);
  ch2=ch*ch;

  phibh = 0.5*log(1.0-2.0*Ms*r/(r2+a2*ch2));
  
  return(phibh);

}













//////////////////////////////////////////////////////////
//
// not sure how to incorporate below into above code
//
/////////////////////////////////////////////////////////

#if(USERHOREF)
// first determine reference density
rhoref_send=-BIG;
COMPZLOOP{
  if(startpos[1]+i==0 && startpos[2]+j==totalsize[2]/2 && startpos[3]+k==totalsize[3]/2){
    get_rhoref(i, j, k, pb, &rhoref_send);
  }
}
// now maximize over all CPUs since assume only 1 CPU got above 
#if(USEMPI)
MPI_Allreduce(&rhoref_send,&rhoref,1,MPI_FTYPE,MPI_MAX,MPI_COMM_GRMHD);
#else
rhoref=rhoref_send;
#endif

#else
//  rhoref=0.0;
#endif





//DEBUG
//    GRAVLOOP(ii) phivsr_tot[ii]=0.0;
//  GRAVLOOP(ii) phivsr_tot[ii]=-0.05/pow(rcent_tot[ii],1);

// The above boundary conditions imply no force beyond the computational domain in the interior and exterior

////////////////////////////////////
//
// add in constant density term if using RHOREF (only for non-rel case)
#if(USERHOREF)
GRAVLOOP(ii){ // now dMvsr done in ghost cells
  //   GRAVLOOPACTIVE(ii){// because don't determine dMvsr in ghost cells, so don't offset potential or Mvsr
  // for rho=constant this is true for spherical polar coordinates
  //    phivsr_tot[ii] += 2.0*M_PI/3.0*(rcent_tot[ii]*rcent_tot[ii]-rcent_tot[0]*rcent_tot[0])*rhoref*Mfactor;

  // since differenced potential above doesn't include inside Rin, we must subtract off interior part of constant density potential
  phivsr_tot[ii] += 2.0*M_PI/3.0*(rcent_tot[ii]*rcent_tot[ii]-localRin*localRin)*rhoref*Mfactor;

  //    dualfprintf(fail_file,"ii=%d rho=%g phi=%g\n",ii,rho,phivsr_tot[ii]);

      
  // since differenced potential above doesn't include inside Rin, we must subtract off interior part of constant density potential
  Mvsr_tot[ii] += 4.0*M_PI/3.0*(myfabs(rcent_tot[ii]*rcent_tot[ii]*rcent_tot[ii])-localRin*localRin*localRin)*rhoref*Mfactor;

  MOrvsr_tot[ii] = Mvsr_tot[ii]/myfabs(rcent_tot[ii]);

}

// recompute MOrvsr

// correct FACE-FACE totals
phirealouter_tot += 2.0*M_PI/3.0*(localRout*localRout-localRin*localRin)*rhoref*Mfactor;
Mrealouter_tot += 4.0*M_PI/3.0*(localRout*localRout*localRout-localRin*localRin*localRin)*rhoref*Mfactor;
#endif













// i,j,k of reference density
// very similar code to within loop just after this is called
int get_rhoref(int i, int j, int k, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE *rhoref)
{
  int jj;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q;
  FTYPE mhd[NDIM][NDIM];
  FTYPE V[NDIM], X[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE idxdxp[NDIM][NDIM];
  int loc=CENT;


  get_geometry(i,j,k,loc,ptrgeom);
  coord_ijk(i, j, k, loc, X);
  bl_coord_ijk(i, j, k, loc, V);
  dxdxprim_ijk(i, j, k, loc, dxdxp);
  idxdxprim_ijk(i, j, k, loc, idxdxp);
    
  //    dualfprintf(fail_file,"mi in selfgravity\n");
  //  matrix_inverse(PRIMECOORDS, dxdxp,idxdxp);
  
  MYFUN(get_state(MAC(pb,i,j,k), ptrgeom, &q) ,"phys.c:source()", "get_state() dir=0", 1);
  DLOOPA(jj) mhd_calc_0(MAC(pb,i,j,k), jj, ptrgeom, &q, mhd[jj], NULL); // includes rest-mass
  
  // now need to convert from PRIMECOORDS to MCOORD
  metptomet_simple_Tud(dxdxp,idxdxp,mhd);

  *rhoref=-mhd[TT][TT];

  return(0);

}







// similar to bounds.c R0SING inside-horizon boundary on primitives
// This sets primitives AND connection inside non-active region in case loops go over non-active region
// With ACTIVEREGION, shouldn't be as important
int bound_spacetime_inside_horizon(void)
{
  int locali,globali;
  int enerregion;
  int i,j,k;
  int ri,rj,rk;
  int jj,kk,ll;
  int loc;
  int pl,pliter;


  if(BCtype[X1DN]==R0SING){ // global condition that all CPUs know about
    // for self-gravity, only loop over region outside horizon
    enerregion=OUTSIDEHORIZONENERREGION;
    enerpos=enerposreg[enerregion];


    
    // copy from very outer boundary or from location on grid defined by where last ghost cell is inside horizon
    // if no black hole, then globali is beyond any real grid points
    // globali is location of fake boundary inside horizon, within which there is no followed dynamics
    // globali indicates position of first followed ghost cell
    globali=N1*horizoncpupos1+horizoni-N1BND;
    if(globali<startpos[1]-N1BND) locali=FLUXNOTONGRID; // indicates not acting on this CPU
    else if(globali>endpos[1]+N1BND-1) locali=N1+N1BND-1; // indicates all of this CPU is acted on
    else locali=horizoni-N1BND; // indicates portion of this CPU acted on
    
    ri=locali;

    //    dualfprintf(fail_file,"locali=%d globali=%d startpos[1]=%d\n",locali,globali,startpos[1]);


    LOOPF{// SECTIONMARK: Should this be FULL or COMPFULL?

      rj=j;
      rk=k;
      
      if(WITHINENERREGION(enerpos,i,j,k)){
        // then do nothing
      }
      else if(locali!=FLUXNOTONGRID){
        // assume horizon on negative side of "i" so don't modify right side of "i" that happens to be in the outer radial boundary

        if(startpos[1]+i < globali){
          // then inside horizon and inside N1BND more grid cells for interpolation
   
          // just outflow all metric and connection quantities in PRIMECOORDS so that outflowed primitives are consistent and have little influence
          // effectively all other points reduce to a duplicate of copied point

          // problem is that connection can change rapidly in time, so for connection set to 0 (i.e. no more acceleration)

          DLOOP(jj,kk) DLOOPA(ll){
            //     GLOBALMETMACP0A3(conn,i,j,k,jj,kk,ll) = GLOBALMETMACP0A3(conn,ri,rj,rk,jj,kk,ll);
            // decided to make space-time effectively flat within fake boundary region
            GLOBALMETMACP0A3(conn,i,j,k,jj,kk,ll) = 0.0;
          }

          DLOOPA(jj){
            //GLOBALMETMACP0A1(conn2,i,j,k,jj) = GLOBALMETMACP0A1(conn2,ri,rj,rk,jj);
            // decided to make space-time effectively flat within fake boundary region
            GLOBALMETMACP0A1(conn2,i,j,k,jj) = 0.0;
          }

   
#if(NEWMETRICSTORAGE)

          // copy over entire structure at once
          for (loc = NPG - 1; loc >= 0; loc--) GLOBALMETMACP1A0(compgeom,loc,i,j,k) = GLOBALMETMACP1A0(compgeom,loc,ri,rj,rk);
          for (loc = NPG - 1; loc >= 0; loc--) GLOBALMETMACP1A0(compgeomlast,loc,i,j,k) = GLOBALMETMACP1A0(compgeomlast,loc,ri,rj,rk);

#else

          for (loc = NPG - 1; loc >= 0; loc--){
            DLOOP(jj,kk){
              GLOBALMETMACP1A2(gcon,loc,i,j,k,jj,kk) = GLOBALMETMACP1A2(gcon,loc,ri,rj,rk,jj,kk);
              // setting last to be equal to current do dg/dt = 0 in this region
              GLOBALMETMACP1A2(gcov,loc,i,j,k,jj,kk) = GLOBALMETMACP1A2(gcov,loc,ri,rj,rk,jj,kk);
#if(DOEVOLVEMETRIC)
              GLOBALMETMACP1A2(gcovlast,loc,i,j,k,jj,kk) = GLOBALMETMACP1A2(gcov,loc,i,j,k,jj,kk) ;
#endif
            }

            DLOOPA(jj){
              // setting last to be equal to current do dg/dt = 0 in this region
              GLOBALMETMACP1A1(gcovpert,loc,i,j,k,jj) = GLOBALMETMACP1A1(gcovpert,loc,ri,rj,rk,jj);
#if(DOEVOLVEMETRIC)
              GLOBALMETMACP1A1(gcovpertlast,loc,i,j,k,jj) = GLOBALMETMACP1A1(gcovpert,loc,i,j,k,jj) ;
#endif
            }

            GLOBALMETMACP1A0(alphalapse,loc,i,j,k) = GLOBALMETMACP1A0(alphalapse,loc,i,j,k) ;
#if(DOEVOLVEMETRIC)
            GLOBALMETMACP1A0(alphalapselast,loc,i,j,k) = GLOBALMETMACP1A0(alphalapse,loc,i,j,k) ;
#endif

            GLOBALMETMACP1A0(gdet,loc,i,j,k) = GLOBALMETMACP1A0(gdet,loc,ri,rj,rk);
            PLOOP(pliter,pl) GLOBALMETMACP1A1(eomfunc,loc,i,j,k,pl) = GLOBALMETMACP1A1(eomfunc,loc,ri,rj,rk,pl);
#if(GDETVOLDIFF)
            GLOBALMETMACP1A0(gdetvol,loc,i,j,k) = GLOBALMETMACP1A0(gdetvol,loc,ri,rj,rk);
#endif
          }
     
     
#endif // end if old metric storage



#if(VOLUMEDIFF)
          DLOOPA(jj) GLOBALMETMACP0A1(idxvol,i,j,k,jj) = GLOBALMETMACP0A1(idxvol,ri,rj,rk,jj);
#endif

   
        }// end if inside globali
      }// end if inside horizoni
    }// end loop over all cells


  } // end if R0SING


  return(0);

}







// determine modification to diagonal components of metric from self-gravity
// assumes original metric is diagonalized and spherical polar coordinates
// assumes r(x_1) only and not r(x_1,x_2)
//
// called by gcov_func() in order to add self-gravity to metric
//
// GODMARK: Apparently might be able to add Mvsr(r) directly to Kerr metric instead of here
//
//
int set_gcov_selfspcmetric(FTYPE *X, FTYPE *V, FTYPE *gcovselfpert)
{
  SFTYPE phi;
  int jj;

  // get interpolated value of phi from grid value of phi
  interpX_phi(X,phivsr_tot,&phi);

  // assign gravitational perturbation in spherical polar coordinates
  gcovselfpert[TT] = -2.0*phi;
  gcovselfpert[RR] = -2.0*phi;
  gcovselfpert[TH] = 0.0;
  gcovselfpert[PH] = 0.0;


  // DEBUG
  //  DLOOPA(jj){
  //    dualfprintf(fail_file,"inside_gcovselfpert[%d]=%21.15g\n",jj,gcovselfpert[jj]);
  //  }


  return(0);
}





// TOV in KS form with velocity
// See grmhd-ksksp.nb, grmhd-connectiononly.nb, tov_solution_timedepsol.nb, tov_solution_timeindepsol.nb, ks_from_bl.nb
int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself)
{
  int jj,kk;
  FTYPE r,r2,th,sh,sh2;
  FTYPE genpot;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;

  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif



  th=V[2];

#if(DOSELFGRAVVSR)
  // get interpolated value of phi from grid value of phi
  interpX_phi(X,phivsr_tot,phiself);

  // get interpolated value of Mvsr/r from grid value of phi
  interpX_phi(X,MOrvsr_tot,MOrself);

  // get interpolated value of vrsqvsr from grid value of phi
  interpX_phi(X,vrsqvsr_tot,vrsqself);
#else
  *phiself=0.0;
  *MOrself=0.0;
  *vrsqself=0.0;
#endif

  r2=rsharp*rsharp;
  sh=sin(th);
  sh2=sh*sh;

  // 0-out metric
  DLOOP(jj,kk) gcov[GIND(jj,kk)]=0.0;

  gcov[GIND(TT,TT)] = -exp(2.0*phiself[0]);
  genpot=vrsqself[0] - 2.0*MOrself[0];
  gcov[GIND(RR,RR)] = 1.0+fabs(-genpot); // [RR][RR] term is positive and be 1 in limit that genpot=0
  // [RR][TT] term must carry sign of r
  gcov[GIND(RR,TT)] = gcov[GIND(TT,RR)] = sqrt(-gcov[GIND(TT,TT)]/(1.0-fabs(genpot)))*(-genpot);
  gcov[GIND(TH,TH)] = r2;
  gcov[GIND(PH,PH)] = r2*sh2;


  gcovpert[TT] = gcov[GIND(TT,TT)]+1.0; // no obvious way to decouple non-linear piece
  gcovpert[RR] = fabs(-genpot);
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;


#if(0)
  // DEBUG:
  DLOOP(jj,kk) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",jj,kk,gcov[GIND(jj,kk)]);
  dualfprintf(fail_file,"r=%21.15g th=%21.15g\n",r,th);
  dualfprintf(fail_file,"phiself=%21.15g MOrself=%21.15g vrsqself=%21.15g :: genpot=%21.15g\n",phiself[0],MOrself[0],vrsqself[0],genpot);
#endif


  return(0);
}







// TOV in BL form with velocity
// See grmhd-ksksp.nb, grmhd-connectiononly.nb, tov_solution_timedepsol.nb, tov_solution_timeindepsol.nb, ks_from_bl.nb
int set_gcov_bl_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself)
{
  int jj,kk;
  FTYPE r,r2,th,sh,sh2;
  FTYPE genpot;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;

  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];

#if(DOSELFGRAVVSR)
  // get interpolated value of phi from grid value of phi
  interpX_phi(X,phivsr_tot,phiself);

  // get interpolated value of Mvsr/r from grid value of phi
  interpX_phi(X,MOrvsr_tot,MOrself);

  // get interpolated value of vrsqvsr from grid value of phi
  interpX_phi(X,vrsqvsr_tot,vrsqself);
#else
  *phiself=0.0;
  *MOrself=0.0;
  *vrsqself=0.0;
#endif


  r2=rsharp*rsharp;
  sh=sin(th);
  sh2=sh*sh;

  // 0-out metric
  DLOOP(jj,kk) gcov[GIND(jj,kk)]=0.0;

  gcov[GIND(TT,TT)] = -exp(2.0*phiself[0]);
  genpot=vrsqself[0] - 2.0*MOrself[0];
  gcov[GIND(RR,RR)] = 1.0/(1.0+fabs(genpot));
  gcov[GIND(RR,TT)] = gcov[GIND(TT,RR)] = 0.0;
  gcov[GIND(TH,TH)] = r2;
  gcov[GIND(PH,PH)] = r2*sh2;


  gcovpert[TT] = gcov[GIND(TT,TT)]+1.0; // no obvious way to decouple non-linear piece
  gcovpert[RR] = gcov[GIND(RR,RR)]-1.0; // no obvious way ...
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;


  //  DLOOP(jj,kk) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",jj,kk,gcov[GIND(jj,kk)]);

  //  dualfprintf(fail_file,"r=%21.15g th=%21.15g\n",r,th);
  //  dualfprintf(fail_file,"phiself=%21.15g MOrself=%21.15g vrsqself=%21.15g :: genpot=%21.15g\n",phiself[0],MOrself[0],vrsqself[0],genpot);

  return(0);
}







// used by set_gcov_selfspcmetric() to determine potential at any location given potential at CENT
// must be a continuous function (i.e. so nearest neighbord won't work -- must use at least linear interpolation)
//
//
// 0 = nearest neighbord: NOT ALLOWED
// 1 = linear : minimum allowed so phi continuous and so Connection well-defined crossing interpolated regions
// 2 = parabolic : Allows Connection to be (itself) continuous (might be good to achieve higher order even if higher order will be obtained via connecting connection (source) at each point by a higher-order function)
// INTERPPHIORDER can't be larger than interpolating using for flux interpolation
#define INTERPPHIORDER 2
//
static int interpX_phi(FTYPE *X, SFTYPE *phigrid, SFTYPE *phi)
{
  int i,j,k;
  extern void icoord(FTYPE *X,int loc, int *i, int *j, int *k);
  extern void icoord_round(FTYPE *X,int loc, int *i, int *j, int *k);
  int ip,ipp;
  FTYPE Xi[NDIM],Xip[NDIM],Xipp[NDIM];
  SFTYPE slope;
  int loc;
  SFTYPE f0,f1,f2,AA,BB;
  int shift1,loweri,upperi;



  loc=CENT; // location where original gridded quantity is defined

  // use icoord() to get i,j,k from X so can get potential at any location from interpolation of existing potential at certain locations

  // get centered i,j,k
#if(0)
  // never got this version to work (for i<0 gets wrong values)
  icoord(X,loc,&i,&j,&k); // truncates to lower i
  if(BCtype[X1DN]==R0SING){
    if(i<0){
      // need to preserve symmetry, so truncate up
      i=i+1-INTERPPHIORDER;
      // will use same 3 points as if on other i-side
    }
  }
#else
  // this version works and is tested
  // with rounding, no need to shift
  icoord_round(X,loc,&i,&j,&k); // rounds so that connection calculation will be accurate since connection around CENT
  if(BCtype[X1DN]==R0SING){
    if(i<0){
      // need to preserve symmetry, so truncate up
      i=i-INTERPPHIORDER;
      // will use same values if quantity is symmetric
    }
  }
#endif




  // defaults
  shift1=SHIFT1;
  loweri=i;
  upperi=loweri+shift1*INTERPPHIORDER;


  // limit lower value of i
  if(startpos[1]+loweri<-N1BND){
    loweri=i=-N1BND;
    upperi=loweri+shift1*INTERPPHIORDER;
  }
  // limit upper value of i
  // startpos[1]+i==ncpux1*N1+N1BND+1 is location of last centered phi
  // +shift1*INTERPPHIORDER is because need phi[i+INTERPPHIORDER] to stay in range of phi array
  if(startpos[1] + upperi > ncpux1*N1+N1BND+1){
    // last possible position of phi given interpolation order
    loweri=i=ncpux1*N1+N1BND+1-shift1*INTERPPHIORDER-startpos[1];
    upperi=loweri+shift1*INTERPPHIORDER;
  }

  ///////////////////////
  // get ith position
  coord_ijk(i,j,k, loc, Xi);


  // get interpolated value

#if(INTERPPHIORDER==0)

  *phi = phigrid[startpos[1]+i];// nearest neighbor (can't do this since leaves jump in potential and dphi/dx is large near jump and will be picked up by connection calculation -- need phi to be smooth)
  
  dualfprintf(fail_file,"Cannot do nearest neighbord with potential since connection large across regions of interpolation\n");
  myexit(2366);

#elif(INTERPPHIORDER==1)

  // get X-position of coordinates i and i+1
  ip=i+shift1;
  coord_ijk(ip,j,k, loc, Xip);

  // use this position to interpolate in X (which happens to be the same as interpolating in i)
  // notice that each cpu accesses a total potential using startpos[1] as the offset
  slope=(phigrid[startpos[1]+ip]-phigrid[startpos[1]+i])/(Xip[1]-Xi[1]);
  *phi = phigrid[startpos[1]+i] + slope*(X[1]-Xi[1]);

  //  dualfprintf(fail_file,"slope=%21.15g X=%21.15g Xi=%21.15g Xip=%21.15g\n",slope,X[1],Xi[1],Xip[1]);

#elif(INTERPPHIORDER==2)

  // get X-position of coordinates i and i+1
  ip=i+shift1;
  coord_ijk(ip,j,k, loc, Xip);

  // get X-position of coordinates i and i+1
  ipp=i+shift1+shift1;
  coord_ijk(ipp,j,k, loc, Xipp);

  // use this position to interpolate in X (which happens to be the same as interpolating in i)
  // notice that each cpu accesses a total potential using startpos[1] as the offset
 
  f0=phigrid[startpos[1]+i];
  f1=phigrid[startpos[1]+ip];
  f2=phigrid[startpos[1]+ipp];

  AA=(f0-f1)/(Xi[1]-Xip[1]) + (f0-f2)/(Xi[1]-Xipp[1]) + (f2-f1)/(Xip[1]-Xipp[1]) ;

  BB=( (f0-f2)/(Xi[1]-Xipp[1]) + (f2-f1)/(Xip[1]-Xipp[1]) )/(Xi[1]-Xip[1]);
  
  *phi = f0 + AA*(X[1]-Xi[1]) + BB*(X[1]-Xi[1])*(X[1]-Xi[1]);

  // DEBUG:
  //  dualfprintf(fail_file,"ip=%d ipp=%d f0=%21.15g f1=%21.15g f2=%21.15g AA=%21.15g BB=%21.15g *phi=%21.15g\n",ip,ipp,f0,f1,f2,AA,BB,*phi);

#endif

  return(0);

}








