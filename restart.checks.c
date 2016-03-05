
/*! \file restart.checks.c
     \brief Functions related to checking if restarting was correct/reasonable
*/


#include "decs.h"


static int restart_init_point_check_pglobal(int which, int i, int j, int k);
static int restart_init_point_check_unewglobal(int which, int i, int j, int k);
static int restart_init_point_check_pstagglobal(int which, int i, int j, int k);


/// only check basic important inputs from restart dump file
int restart_init_simple_checks(int which)
{
  int gotnan;
  int i,j,k;

  //////////////
  //
  // make sure all zones are not nan just as read-in from file
  //
  //////////////
  // make sure all zones are not nan
  gotnan=0;
  // OPENMPOPTMARK: Don't optimize since critical region

  if(which<=2){ // see initbase.c
    LOOP{  gotnan+=restart_init_point_check_pglobal(which,i,j,k); }
    LOOP{  gotnan+=restart_init_point_check_unewglobal(which,i,j,k); }
  }
  else if(which<=3){
    FULLLOOP{  gotnan+=restart_init_point_check_pglobal(which,i,j,k); }
    LOOP{  gotnan+=restart_init_point_check_unewglobal(which,i,j,k); }
  }
  else{
    FULLLOOP{  gotnan+=restart_init_point_check_pglobal(which,i,j,k); }
    LOOP{  gotnan+=restart_init_point_check_unewglobal(which,i,j,k); }
    FULLLOOP{  gotnan+=restart_init_point_check_pstagglobal(which,i,j,k); }
  }



  if(gotnan) myexit(39476346);

  return(0);

}


/// check pglobal
static int restart_init_point_check_pglobal(int which, int i, int j, int k)
{
  int pliter,pl;
  int gotnan;

  gotnan=0;

  PDUMPLOOP(pliter,pl){
    if(!finite(GLOBALMACP0A1(pglobal,i,j,k,pl)) ){
      dualfprintf(fail_file,"restart_init(%d): restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d : pglobal=%21.15g\n",which,i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl,GLOBALMACP0A1(pglobal,i,j,k,pl));
      if(SCALARPL(pl)){
        dualfprintf(fail_file,"scalar went nan, reset to floor: pl=%d\n",pl);
        GLOBALMACP0A1(pglobal,i,j,k,pl)=NUMEPSILON;
      } 
     else {
     // myexit(24968346);
      gotnan++;
   }
   }
   }

  return(gotnan);

}


/// also check unewglobal for that portion that's used
static int restart_init_point_check_unewglobal(int which, int i, int j, int k)
{
  int pliter,pl;
  int gotnan;

  gotnan=0;

  PDUMPLOOP(pliter,pl){
    if(DOENOFLUX != NOENOFLUX || (FLUXB==FLUXCTSTAG &&(pl==B1 && i>=-N1BND+SHIFT1 || pl==B2 && j>=-N2BND+SHIFT2 || pl==B3 && k>=-N3BND+SHIFT3)) ){
      if(!finite(GLOBALMACP0A1(unewglobal,i,j,k,pl)) ){
        dualfprintf(fail_file,"restart_init(%d): restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d : unewglobal=%21.15g\n",which,i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl,GLOBALMACP0A1(unewglobal,i,j,k,pl));
        // myexit(24968346);
        gotnan++;
      }
    }
  }

  return(gotnan);

}


/// also check pstagglobal once computed
static int restart_init_point_check_pstagglobal(int which, int i, int j, int k)
{
  int pliter,pl;
  int gotnan;

  gotnan=0;

  
  if(FLUXB==FLUXCTSTAG){
    PDUMPLOOP(pliter,pl){
      if(pl==B1 && i>=-N1BND+SHIFT1 || pl==B2 && j>=-N2BND+SHIFT2 || pl==B3 && k>=-N3BND+SHIFT3){
        if(!finite(GLOBALMACP0A1(pstagglobal,i,j,k,pl)) ){
          dualfprintf(fail_file,"restart_init(%d): restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d : pstagglobal=%21.15g\n",which,i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl,GLOBALMACP0A1(pstagglobal,i,j,k,pl));
          // myexit(24968346);
          gotnan++;
        }
      }
    }
  }


  return(gotnan);

}




//////////////////////
///
/// At this point all grid type parameters should be set as if done with init()
///
/// Perform some extra checks to ensure restart file read-in is reasonable
///
/// OPENMPOPTMARK: Don't optimize since many critical regions
///
//////////////////////
int restart_init_checks(int which, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR])
{
  char ans[100];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int failreturn;
  FTYPE ucon[NDIM];
  FTYPE utmax=0.0;
  int i,j,k,pl,pliter;
  int failflag=0;
  extern int set_dt(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE *dt);
  int gotnan;



  ////////////////////////////////////////////////////////////////
  //
  // NOW have all parameters and data.
  //
  /////////



  //////////////
  //
  // make sure all zones are not nan before fixup and bound
  //
  //////////////
  // make sure all zones are not nan
  gotnan=0;
  LOOP{// OPENMPOPTMARK: Don't optimize since critical region
    PDUMPLOOP(pliter,pl){
      if(!finite(MACP0A1(prim,i,j,k,pl))){
        dualfprintf(fail_file,"before fixup & bound: restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d\n",i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
        // myexit(24968346);
        gotnan=1;
      }
    }
  }
  if(gotnan) myexit(24968341);



  //////////
  // NOW CHECK THE DATA and apply "fixups" if necessary
  //
  ////////////////////////////////////////////////////////////////

#if(CHECKRHONEGZERORESTART)

  // report if data has rho<0 or u<0
  ZLOOP{
    if(STEPOVERNEGRHO!=NEGDENSITY_NEVERFIXUP){
      if(MACP0A1(prim,i,j,k,RHO)<=0.0){
        dualfprintf(fail_file,"restart data has negative mass density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
        failflag++;
      }
    }
    if(STEPOVERNEGU!=NEGDENSITY_NEVERFIXUP){
      if(MACP0A1(prim,i,j,k,UU)<=0.0){
        dualfprintf(fail_file,"restart data has negative ie density at i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
        failflag++;
      }
    }

    // check 3-velocity
#if(WHICHVEL==VEL3)
    if(jonchecks){
      get_geometry(i,j,k,CENT,ptrgeom);
      failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),MAC(ucons,i,j,k), ptrgeom,-2,-1);
      if(failreturn){
        dualfprintf(fail_file,"restart data has large or imaginary u^t=%21.15g at i=%d j=%d k=%d.  I will attempt to correct.\n",1.0/sqrt(uttdiscr),startpos[1]+i,startpos[2]+j,startpos[3]+k);
      }
      if(1.0/sqrt(uttdiscr)>utmax) utmax=1.0/sqrt(uttdiscr);
      // need to settle over limit u^t's
      failreturn=check_pr(MAC(prim,i,j,k),MAC(prim,i,j,k),MAC(ucons,i,j,k), ptrgeom,-1,-1);
      if(failreturn){
        dualfprintf(fail_file,"restart data has imaginary u^t at i=%d j=%d k=%d.  Unable to correct.\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
        return(1);
      }
    }
#endif    
  }
  if(WHICHVEL==VEL3){
    if(jonchecks){
      dualfprintf(fail_file,"max u^t of restart data=%21.15g\n",utmax);
    }
  }
#endif

  if(failflag>0){
    dualfprintf(fail_file,"Restart data has at least %d failures, please correct or accept that fixup will fix them.\n",failflag);
    //myexit(1); // aborts MPI
    //    return(1);  // currently doesn't abort properly in MPI
  }




  /////////////////////
  //
  // fixup() during restart
  //
  /////////////////////
#if(FIXUPAFTERRESTART)
  if(fixup(STAGEM1,prim,ucons,-1)>=1)
    FAILSTATEMENT("restart.c:restart_init()", "fixup()", 1);

  trifprintf( "proc: %d fixup restart completed: failed=%d\n", myid,failed);
#endif



  //////////////
  //
  // make sure all zones are not nan after fixup
  //
  //////////////
  // make sure all zones are not nan
  gotnan=0;
  LOOP{// OPENMPOPTMARK: Don't optimize since critical region
    PDUMPLOOP(pliter,pl){
      if(!finite(MACP0A1(prim,i,j,k,pl))){
        dualfprintf(fail_file,"after fixup & before bound: restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d\n",i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
        // myexit(24968346);
        gotnan=1;
      }
    }
  }
  if(gotnan) myexit(24968341);



  /////////////////////
  //
  // BOUND during restart
  //
  /////////////////////
  int finalstep=1; // user would want to know about changes to conserved quants during restart
  if (bound_allprim(STAGEM1,finalstep,t,prim,pstag,ucons, USEMPI) >= 1) {
    dualfprintf(fail_file, "restart_init:bound_allprim: failure\n");
    fflush(fail_file);
    return (1);
  }

  trifprintf( "proc: %d bound restart completed: failed=%d\n", myid,failed);


  //////////////
  //
  // make sure all zones are not nan after bound
  //
  //////////////
  // make sure all zones are not nan
  gotnan=0;
  FULLLOOP{// OPENMPOPTMARK: Don't optimize since critical region
    PDUMPLOOP(pliter,pl){
      if(!finite(MACP0A1(prim,i,j,k,pl))){
        dualfprintf(fail_file,"after fixup & bound: restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d\n",i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
        // myexit(24968346);
        gotnan=1;
      }
    }
  }
  if(gotnan) myexit(24968346);

  /////////////////////
  //
  // pre_fixup() during restart
  //
  /////////////////////
  
  if(pre_fixup(STAGEM1,prim)>=1)
    FAILSTATEMENT("init.c:init()", "postbc_fixup()", 1);



  /////////////////////
  //
  // make sure all zones are good now
  //
  /////////////////////
#if(CHECKRHONEGZERORESTART)
  failflag=0;

  if(STEPOVERNEGRHO!=NEGDENSITY_NEVERFIXUP){
    FULLLOOP{
      if(MACP0A1(prim,i,j,k,RHO)<=0.0){
        dualfprintf(fail_file,"restart data has negative mass density at i=%d j=%d k=%d : %21.15g\n",startpos[1]+i,startpos[2]+j,startpos[3]+k,MACP0A1(prim,i,j,k,RHO));
        // return(1);
        failflag++;
      }
    }// end fullloop
  }// end negrho check

  if(STEPOVERNEGU!=NEGDENSITY_NEVERFIXUP){
    FULLLOOP{
      if(MACP0A1(prim,i,j,k,UU)<=0.0){
        dualfprintf(fail_file,"restart data has negative ie density at i=%d j=%d k=%d : %21.15g\n",startpos[1]+i,startpos[2]+j,startpos[3]+k,MACP0A1(prim,i,j,k,UU));
        // return(1);
        failflag++;
      }
    }// end fullloop
  }

  if(failflag>0){
    dualfprintf(fail_file,"Restart data has at least %d failures -- even after fixup should have been applied!\n",failflag);
    myexit(1); // aborts MPI (abort, not return, so no stall on other cores)
  }

#endif



  ////////////////////////
  //
  // test read by looking at images
  //
  ////////////////////////
  if(image_dump(-3)>=1) return(1);



  ////////////////////////
  //
  // Now adjust timestep since may be dt=0 or whatever if last run ended on tf
  // NOW done whether restarting or not at end of init()
  //
  ///////////////////////
  // set_dt(prim,&dt);


  trifprintf("end restart_init_checks\n");

  /* done! */
  return (0);

}
