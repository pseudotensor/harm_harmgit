


#include "decs.h"

// Note that if using primtoU() or primtoflux(), these respect nprlist for PLOOP, hence why these tests use PLOOP while others use PALLLOOP



// fractional error above which inversion is reported on as having made a signficant error
#if(PRODUCTION)
#define CHECKONINVFRAC (1E-1)
#else
#define CHECKONINVFRAC (1E-2)
#endif

// whether to fail if check on inversion fails
#define FAILIFBADCHECK 1

// fraction upon if greater will treat as failure if FAILIFBADCHECK is triggered
#define CHECKONINVFRACFAIL (1E-1)

static int negdensitycheck(int finalstep, FTYPE *prim, PFTYPE *pflag);
static int check_on_inversion(int usedhotinversion,int usedentropyinversion,int usedcoldinversion,int usedffdeinversion, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
static int debug_utoprimgen(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew);
static int compare_ffde_inversions(int showmessages, int allowlocalfailurefixandnoreport, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_geom *ptrgeom, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);

static int tryentropyinversion(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
static int trycoldinversion(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);


int Utoprimgen(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  int whichentropy;
  //  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  int otherfail;
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int Utoprimgen_compare(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int Utoprimgen_tryagain(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int Utoprimgen_tryagain2(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  void convert_U_removerestmassfromuu(int utoprimverison, int removerestmassfromuu, FTYPE *U);
  int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree,FTYPE *pr);
  int pl,pliter;
  PFTYPE lpflag,lpflagrad;
  int usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion;
  FTYPE pressuremem;
  FTYPE *pressure=&pressuremem;
  int eomtypelocal;







  // DEBUG:
  //  dualfprintf(fail_file,"Doing inversion for ijk=%d %d %d nstep=%ld steppart=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);


 
  ////////////////////////
  //
  // INVERT Magnetic field or restore magnetic field for SPLITNPR case
  //
  ////////////////////////

#if(SPLITNPR)
  if(nprlist[nprstart]==B1 && nprlist[nprend]==B3){ // then assume doing B1,B2,B3
    // invert magnetic field only
    PLOOPBONLY(pl) pr[pl]=U[pl]/ptrgeom->e[pl]; // direct inversion from geometry-free conserved quantity
    return(0);
  }
  else{
    // otherwise assume want to invert ONLY non-B components and assume using same memory space
    // can achieve this by setting Ugeomfree to *existing* pr (assume memory space for SPLITNPR step 0 and 1 didn't change
    // overwrites update (that is no field update done on second pass) to magnetic field
    // Note that field U not updated in this pass, so need to fill with something.
    // Assume pr here is pbb that contains field-update's Bnew
    PLOOPBONLY(pl) U[pl]=pr[pl]*ptrgeom->e[pl];
  }
#else
  // invert magnetic field
  // direct inversion from geometry-free conserved quantity
  // assume normal inversion takes care of this
  //  PLOOPBONLY(pl) pr[pl]=U[pl]/ptrgeom->e[pl];
#endif









  ///////////
  //
  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  //
  ///////////
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);
  // e.g. if REMOVERESTMASSFROMUU==1 and inputtype==UEVOLVE, then Ugeomfree has energy T^t_t that includes rest-mass
#if(SPLITNPR)
  if(!(nprlist[nprstart]==B1 && nprlist[nprend]==B3)){
    // then need to update field quantities since UtoU uses PLOOP not PALLLOOP -- GODMARK: could use PALLLOOP in UtoU, but then adds excessive code to phys.c related to SPLITNPR
    PLOOPBONLY(pl) Ugeomfree[pl]=U[pl]/ptrgeom->e[pl];
  }
#endif



  /////////////////////////////////////////////////////////
  //
  // Copy over Ugeomfree/pr to Ugeomfree0/pr0
  //
  ////////////////////////////////////////////////////////
  PALLLOOP(pl){
    Uold[pl]=Ugeomfree0[pl]=Ugeomfree[pl];
    pr0[pl]=pr[pl];
  }




  /////////////////////////////////////////////////////////////
  //
  // convert U into form appropriate for inversion routine
  //
  /////////////////////////////////////////////////////////////
  convert_U_removerestmassfromuu(UTOPRIMVERSION,REMOVERESTMASSFROMUU,Ugeomfree);
  convert_U_removerestmassfromuu(UTOPRIMVERSION,REMOVERESTMASSFROMUU,Ugeomfree0);




  ///////////////
  //
  //  if(nstep==4 && steppart==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"PRIMUTOPRIMGEN0(%d): pl=%d pp=%21.15g uu=%21.15g\n",*eomtype,pl,pr[pl],Ugeomfree[pl]);
  //  }

  static long long int didnothing=0,didsomething=0;
  ///////////////
  if(EOMDONOTHING(*eomtype) && (finalstep==0 || TIMEORDER<=3)){


    // if finalstep==0, then don't do inversion.  If finalstep==1, only need to do inversion if ucum is different than final uf.
    // then do nothing since assume already pr=pr[U].
    // also don't change any failure flags, assuming prior flags are correct.

    // still can check inversion
    if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)<=UTOPRIMNOFAIL){
      usedhotinversion=usedentropyinversion=usedcoldinversion=usedffdeinversion=0;
      if(*eomtype==EOMDIDGRMHD) usedhotinversion=1;
      if(*eomtype==EOMDIDENTROPYGRMHD) usedentropyinversion=1;
      if(*eomtype==EOMDIDCOLDGRMHD) usedcoldinversion=1;
      if(*eomtype==EOMDIDFFDE) usedffdeinversion=1;

      // check cold inversion if inversion thinks it was successful
      check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }
    if(DOEVOLVERHO){
      negdensitycheck(finalstep, pr, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL));
    }
    didnothing++;
    if(debugfail>=2) if(didnothing%totalzones==0) dualfprintf(fail_file,"DIDNOTHING: %lld : %ld %d : %d\n",didnothing,nstep,steppart,*eomtype);
    return(0);
  }
  else{
    // setup default eomtype
    if(*eomtype==EOMDEFAULT) eomtypelocal=EOMTYPE; // use default
    // if did certain inversion on last sub-step, then assume want same inversion type (e.g. for finalstep==1 for TIMEORDER>2 cases)
    else if(*eomtype==EOMDIDGRMHD) eomtypelocal=EOMGRMHD;
    else if(*eomtype==EOMDIDENTROPYGRMHD) eomtypelocal=EOMENTROPYGRMHD;
    else if(*eomtype==EOMDIDCOLDGRMHD) eomtypelocal=EOMCOLDGRMHD;
    else if(*eomtype==EOMDIDFFDE) eomtypelocal=EOMFFDE;
    else eomtypelocal=*eomtype; // force eomtype

    didsomething++;
    if(debugfail>=2) if(didsomething%totalzones==0) dualfprintf(fail_file,"DIDSOMETHING: %lld : %ld %d : eomtype=%d -> eomtypelocal=%d\n",didsomething,nstep,steppart,*eomtype,eomtypelocal);
  }













  ////////////////////////////////////////////////////////
  //
  // START INVERSION PROCESS
  //
  //
  ////////////////////////////////////////////////////////





  ////////////////////////
  //
  // INVERT pseudo-passive scalars
  //
  // all pseudo-passive scalars are trivially inverted
  //
  // So don't include passive scalars in check_on_inversion() [Since modify passive scalars]
  //
  ////////////////////////
  invert_scalars(ptrgeom, Ugeomfree,pr);





  /////////////////////////////////////////////////////////
  //
  // backup (in inversion routine format with already-inverted scalars)
  //
  ////////////////////////////////////////////////////////
  PALLLOOP(pl){
    Uold[pl]=Ugeomfree0[pl]=Ugeomfree[pl];
    pr0[pl]=pr[pl];
  }
  










  ////////////////////////
  //
  // DO non-scalar INVERSION (from here after, only modify RHO,UU,U1,U2,U3,B1,B2,B3.  Any initial guess should be using pr0 that contains other interpolated quantities.)
  //
  ////////////////////////

  // defaults
  usedhotinversion=0;
  usedentropyinversion=0;
  usedcoldinversion=0;
  usedffdeinversion=0;


  // if some process already set fail flag (e.g. inversions somewhere else that wasn't reset) (e.g. like in phys.tools.rad.c for implicit or explicit solvers with radiation), then assume want to treat as failure for purposes of fixing up.  That is, get inversion, but that inversion will be for undefined 4-force and only used in worst-case scenario for fixups.
  PFTYPE preexistingfailgas=0;
  if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)>UTOPRIMNOFAIL) preexistingfailgas=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);

  // Also see if radiation inversion got locally corrected.  Then, might still do better to do multi-point averaging in fixups like with MHD variables.  So capture this flag for fixups.
  PFTYPE preexistingfailrad=0;
  if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)>UTOPRIMRADNOFAIL) preexistingfailrad=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);





  if(eomtypelocal==EOMGRMHD){
    ///////////////////////////////////////////////////
    //
    ///////////// HOT GRMHD
    //
    ///////////////////////////////////////////////////

    // report back that used this
    *eomtype=eomtypelocal;


    if(UTOPRIMVERSION!=UTOPRIMCOMPARE) Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMVERSION, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    else Utoprimgen_compare(showmessages, allowlocalfailurefixandnoreport, eomtypelocal, EVOLVENOENTROPY,Ugeomfree,ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    usedhotinversion=1;
    

    //    if(pr[UU]>1.0) dualfprintf(fail_file,"0DEATH\n");
    //    if(pr[UU]<1E-13) dualfprintf(fail_file,"0DEATH2\n");
    //    if(Ugeomfree0[ENTROPY]<-10) dualfprintf(fail_file,"0DEATH3\n");




    // try other methods (assumes all methods can handle WHICHVEL, etc. used for primary model)
    // right now, all can handle WHICHVEL==VELREL4 and energy equation evolution and REMOVERESTMASSFROMUU=0,1
#if((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU<=1)&&(UTOPRIMTRYAGAIN))
    Utoprimgen_tryagain(showmessages, allowlocalfailurefixandnoreport, eomtypelocal, EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
#elif((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU==2)&&(UTOPRIMTRYAGAIN))
    // Can only try again using same type of U since tryagain code doesn't convert U 
    Utoprimgen_tryagain2(showmessages, allowlocalfailurefixandnoreport, eomtypelocal, EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
#endif


    PFTYPE hotpflag;
    // get failure flag
    hotpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
    int hardfailure=(IFUTOPRIMFAIL(hotpflag) && IFUTOPRIMFAILSOFT(hotpflag)==0);

    if(!hardfailure){
      // check on hot inversion if it thinks it was successful
      check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }


    // normally don't ever do this unless really debugging inversion.
    if(0&&debugfail>=2&&hardfailure){ // DEBUG: only report if hard failure.  Seems when gets negative density or internal energy, jon's inversion is fine.
      // first report info so can check on inversion
      extern int mathematica_report_check(int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, FTYPE errorabs, int iters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);
      static long long int failnum=0;
      FTYPE fakedt=0.0; // since no 4-force
      FTYPE fakeCUf[4]={0}; // fake
      FTYPE dUother[NPR]={0};// fake
      struct of_state *qptr=NULL; // fake
      failnum++;
      // fake ppfirst as pr
      mathematica_report_check(2, failnum, (int)hotpflag, eomtypelocal, 0.0, -1, fakedt, ptrgeom, pr, pr, pr, pr0, pr0, pr0, Ugeomfree0, Ugeomfree, Ugeomfree0, Ugeomfree0, fakeCUf, qptr, dUother);
    }


    ////////////////////
    // If hot GRMHD failed or gets suspicious solution, revert to entropy GRMHD if solution
    // If radiation, then this redoes radiation inversion since entropy would give new velocity and local corrections in u2p_rad() might use velocity.
    ///////////////////
    if(HOT2ENTROPY){
      // get failure flag
      hotpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);


      if(hotpflag){
        tryentropyinversion(showmessages, allowlocalfailurefixandnoreport,finalstep, hotpflag, pr0, pr, pressure, Ugeomfree, Ugeomfree0, ptrgeom,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)<=UTOPRIMNOFAIL){
          // report back that used this
          *eomtype=EOMENTROPYGRMHD;

          usedhotinversion=0;
          usedentropyinversion=1;
          usedcoldinversion=0;
          usedffdeinversion=0;

          // check entropy inversion if inversion thinks it was successful
          check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        }
          
      }
    }




    ////////////////////
    //  If hot GRMHD failed or gets suspicious solution, revert to cold GRMHD if solution is cold
    ///////////////////
    if(HOT2COLD){
      // get failure flag
      hotpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);

      if(hotpflag){
        trycoldinversion(showmessages, allowlocalfailurefixandnoreport,finalstep, hotpflag, pr0, pr, pressure, Ugeomfree, Ugeomfree0, ptrgeom,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)<=UTOPRIMNOFAIL){
          usedhotinversion=0;
          usedentropyinversion=0;
          usedcoldinversion=1;
          usedffdeinversion=0;

          // report back that used this
          *eomtype=EOMCOLDGRMHD;
        
          // check cold inversion if inversion thinks it was successful
          check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        }

      }// end if hotpflag

    }




  }
  else if(eomtypelocal==EOMENTROPYGRMHD){
    ///////////////////////////////////////////////////
    //
    ///////////// ENTROPY GRMHD
    //
    // direct entropy evolution (can use old Utoprim() or new code, but not all codes have entropy inversion)
    //
    ///////////////////////////////////////////////////

    // report back that used this
    *eomtype=eomtypelocal;

    // setup as if full entropy evolution
    whichentropy=EVOLVEFULLENTROPY;


#if(0)
    // GODMARK: I noticed that old 5D method is very poor in finding the solution in large gradients so fails alot leading to averaging and run-away field growths in semi-static regions leading to untrustable solution.
    // original entropy inversion method (works fine)

    // do inversion for entropy version of EOMs
    // only one inversion is setup to handle this
    PALLLOOP(pl) prother[pl]=pr0[pl];

    
    ////////////////////////
    // get entropy evolution (don't use failure -- otherfail)
    Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtypelocal, whichentropy, Uold, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));

#else


    // new faster code within utoprim_jon.c
    // This method seems very robust, hardly every failing to find solution, indicating that 5D method's inability to find solution was intrinsic to the 5D method NOT to the new conserved quantities.


    // get entropy evolution inversion
    Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtypelocal, whichentropy, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));


    usedentropyinversion=1;

    // check entropy inversion
    check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
 

#if(0)
    // DEBUG
    // If utoprim_jon.c fails to find solution, then see if old 5D method finds solution.  If so, then complain so JCM can improve new inversion
    lpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
    if(IFUTOPRIMFAIL(lpflag)){

      // copy over utoprim_jon result for check_on_inversion below
      FTYPE prorig[NPR],pr0orig[NPR],Uoldorig[NPR],Uneworig[NPR];
      PLOOP(pliter,pl){
        pr0orig[pl]=pr0[pl];
        prorig[pl]=pr[pl];
        Uoldorig[pl]=Uold[pl];
        Uneworig[pl]=Unew[pl];
      }

      // Get original inversion for entropy
      Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtypelocal, whichentropy, Uold, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));

      if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)==0){
        check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&lpflag, pr0orig, prorig, pressure, ptrgeom, Uoldorig, Uneworig,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)); // checks/outputs utoprim_jon.c, not original.  But only wanted outputted if original method succeeds where new fails
        dualfprintf(fail_file,"Old inversion method worked while new failed\n");
        myexit(0);
      }
    }
#endif






    ////////////////////
    //  If entropy GRMHD failed or gets suspicious solution, revert to cold GRMHD if solution is cold
    //  This is the same trycoldinversion() as for HOT2COLD
    ///////////////////
    if(ENTROPY2COLD){
      PFTYPE entropypflag;
      // get failure flag
      entropypflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
 
      if(entropypflag){
        trycoldinversion(showmessages, allowlocalfailurefixandnoreport,finalstep, entropypflag, pr0, pr, pressure, Ugeomfree, Ugeomfree0, ptrgeom,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)<=UTOPRIMNOFAIL){
          usedhotinversion=0;
          usedentropyinversion=0;
          usedcoldinversion=1;
          usedffdeinversion=0;

          // report back that used this
          *eomtype=EOMCOLDGRMHD;

          // check cold inversion
          check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
        }

      }// end if entropypflag
 
    }




#endif

  }
  else if(eomtypelocal==EOMCOLDGRMHD){
    ///////////////////////////////////////////////////
    //
    ///////////// COLDGRMHD
    //
    ///////////////////////////////////////////////////

    // report back that used this
    *eomtype=eomtypelocal;


    if(1){
      // Jon's inversion
      Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }
    else if(0){ // not working yet
      Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMCOLDGRMHD, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }

    usedcoldinversion=1;

    // check cold inversion
    check_on_inversion(usedhotinversion,usedentropyinversion,usedcoldinversion,usedffdeinversion,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));



  }
  else if(eomtypelocal==EOMFFDE){
    ///////////////////////////////////////////////////
    //
    ///////////// FORCE FREE
    //
    ///////////////////////////////////////////////////

    // report back that used this
    *eomtype=eomtypelocal;



    // GODMARK: inversions lead to different behavior (start with torus with rho=u=0 but loop of field)!
    
    if(0){ // Jon's old inversion
      Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMFFDE, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }
    else if(1){
      // Jon's inversion
      Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr,pressure,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }
    else if(0){
      compare_ffde_inversions(showmessages, allowlocalfailurefixandnoreport,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr0, pr, pressure, ptrgeom, Ugeomfree0, Ugeomfree, Uold, Unew,newtonstats,&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL));
    }

    usedffdeinversion=1;


  }//end if EOMFFDE






  ////////////
  //
  // modify failure flag if necessary
  //
  ////////////

#if(DOEVOLVERHO)
  negdensitycheck(finalstep, pr, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL));
#endif






  ///////////////////////////////////////////////////
  //
  ///////////// Report failure
  //
  // complain if pflag set
  // only output if failed
  //
  ///////////////////////////////////////////////////

  // for now only report if not just negative density failure
  lpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);

  lpflagrad=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);



  if(IFUTOPRIMFAILSOFT(lpflag) || IFUTOPRIMRADFAIL(lpflagrad)){
    // then don't report info for now SUPERGODMARK
    if(showmessages&&debugfail>=1) dualfprintf(fail_file, "SOFT Failed to find a prim. var. solution on finalstep=%d!! t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d : errx=%21.15g lpflagrad=%d\n",finalstep,t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag,newtonstats->lerrx,lpflagrad);
  }
  else if((IFUTOPRIMFAIL(lpflag) || IFUTOPRIMRADFAIL(lpflagrad)) &&(debugfail>=1)){
    if(showmessages) dualfprintf(fail_file, "Failed to find a prim. var. solution on finalstep=%d!! t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d : errx=%21.15g lpflagrad=%d\n",finalstep,t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag,newtonstats->lerrx,lpflagrad);
  }
  else{
    //    if(showmessages) dualfprintf(fail_file, "NO FAIL t=%21.15g steppart=%d nstep=%ld i=%d j=%d k=%d : fail=%d : errx=%21.15g\n",t,steppart,realnstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag,newtonstats->lerrx,lpflagrad);
  }






#if(PRODUCTION==0)
  ///////////////////////////////////////////////////
  //
  ///////////// LOTS OF DEBUG STUFF
  //
  ///////////////////////////////////////////////////
  debug_utoprimgen(&lpflag, pr0, pr, ptrgeom, Uold, Unew);
#endif




  ///////////////////////////////////////////////////
  //
  ///////////// override fail flags with preexisting conditions, with assumption that wanted to get inversion just for fixups to have something to revert.  In case of radiation-mhd simulations, this will be without the correct 4-force, which isn't necessarily a good backup and can itself be a bad backup.
  //
  ///////////////////////////////////////////////////
  if(preexistingfailgas) GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=preexistingfailgas;
  if(preexistingfailrad) GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=preexistingfailrad;
  
  // KORALTODO: Could further decide to use original pf from some point in source explicit/implicit updates instead of this no-force solution.



  ///////////////
  //
  //  if(nstep==4 && steppart==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"PRIMUTOPRIMGEN1 (%d): pl=%d pp=%21.15g uu=%21.15g\n",*eomtype,pl,pr[pl],Ugeomfree[pl]);
  //  }

     
  return(0);

}











// use entropy grmhd if hot grmhd gives u<zerouuperbaryon*rho since then entropy is valid
#define USEENTROPYIFHOTUNEG 1
// use entropy grmhd if hot grmhd gives rho<0 since rho<0 is nearly no solution
#define USEENTROPYIFHOTRHONEG 1
// use entropy grmhd if hot grmhd leads to convergence or other "no solution" type failure
#define USEENTROPYIFHOTFAILCONV 1

// try entropy inversion of hot one fails
int tryentropyinversion(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int pl;
  FTYPE prhot[NPR],prentropy[NPR];
  PFTYPE entropypflag;
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int eomtypelocal=EOMENTROPYGRMHD;

  //  dualfprintf(fail_file,"Got here in tryentropyinversion\n");

  if( (USEENTROPYIFHOTRHONEG==0 && IFUTOPRIMFAILSOFTRHORELATED(hotpflag)) || (USEENTROPYIFHOTUNEG==0 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) ){
    // then maybe not so bad failure
    // e.g. if get here, do nothing when either rho<=0 or u<=0
  }
  else if( (USEENTROPYIFHOTRHONEG==1 && IFUTOPRIMFAILSOFTRHORELATED(hotpflag)) || (USEENTROPYIFHOTUNEG==1 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) || (USEENTROPYIFHOTFAILCONV==1 && IFUTOPRIMFAIL(hotpflag)) ){
    // get here if want to fix up rho<=0, u<=0, or convergence failure
    // First if() is needed since last conditional or else if() would trigger on any failure type
    

    // then bad failure, so try to use entropy grmhd
    // restore backup in case previous inversion changed things
    PALLLOOP(pl){
      prhot[pl]=pr[pl];
      Ugeomfree[pl]=Ugeomfree0[pl];

      // setup input guess and other already-inverted solutions
      prentropy[pl]=pr0[pl];
    }
    
    
    // get entropy evolution inversion
    Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &entropypflag, prentropy,pressure,newtonstats, lpflagrad);

    
    ///////////////////////////////////
    //
    // check if entropy solution is good
    if(IFUTOPRIMNOFAILORFIXED(entropypflag)){


      // set pflag so diag_fixup knows used entropy inversion
      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMFAILFIXEDENTROPY;


      // then keep entropy solution
      PALLLOOP(pl) pr[pl]=prentropy[pl];


      // accounting (since fixup.c accounting doesn't know what original pr0 is and actual prhot is undefined since Uhot->prhot wasn't possible)
      // GODMARK: For DOENOFLUX>0, should modify conserved quantity in some way.  For DOENOFLUX==0, primitives form basis of conserved quantities, so once primitives are modified all is done.  So should probably pass in U[] from  Utoprimgen() or whatever is the full conserved quantity later used.
      int modcons=0;
      FTYPE Uievolve[NPR];
      // ucons not modified (i.e. modcons=0), but ucons may be used by diag_fixup()
      UtoU(UNOTHING,UEVOLVE,ptrgeom,Ugeomfree0,Uievolve);
      // account for change to hot MHD conserved quantities
      diag_fixup_Ui_pf(modcons,Uievolve,pr,ptrgeom,finalstep,COUNTENTROPY);
      
      // KORALTODO: allowlocalfailurefixandnoreport could be used below to avoid reset, but for now let implicit solver be ok with entropy inversion
      // reset pflag, unless still need to do some kind of averaging
      if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)==UTOPRIMFAILFIXEDENTROPY){
        // default then is that no failure
        // already accounted for and just a soft failure that doesn't require fixup anymore (and can't leave as UTOPRIMFAILFIXEDUTOPRIM since fixup would complain about not resetting the flag)
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
      }



      // but set internal energy to previous value (should really evolve with entropy equation, but if negligible and no strong shocks, then ok )
      // GODMARK: another ok way is to set special failure that only averages internal energy!  Then can evolve at least -- in some diffusive way
#if(PRODUCTION==0)      
      if(debugfail>=2 && showmessages) dualfprintf(fail_file,"Tried entropy and good on finalstep=%d! steppart=%d nstep=%ld i=%d j=%d k=%d :: hotpflag=%d entropypflag=%d\n",finalstep,steppart,nstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,hotpflag,entropypflag);

      if(debugfail>=2){
        //        int pliter;
        //        if(pr[UU]>1.0) dualfprintf(fail_file,"PDEATH\n");
        //        if(pr[UU]<1E-13) dualfprintf(fail_file,"PDEATH2\n");
        //        if(Ugeomfree0[ENTROPY]<-10) dualfprintf(fail_file,"PDEATH3\n");
        //        PLOOP(pliter,pl) dualfprintf(fail_file,"POST GO TO ENTROPY: ijk=%d %d %d nstep=%ld steppart=%d pl=%d pr0=%g pr=%g U0=%g U=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,pl,pr0[pl],pr[pl],Ugeomfree0[pl],Ugeomfree[pl]);

      }
#endif

      

    }// else if entropypflag is bad
    else{
      // then both hot and entropy are bad, so keep hot
      // GODMARK: Could try going to force-free and "failing" the parallel velocity so it gets averaged like internal energy in entropy case!
      // only can go to force-free if b^2/rho>>1 as well
      // keep hotpflag and keep hot solution
      PALLLOOP(pl) pr[pl]=prhot[pl];

#if(PRODUCTION==0)      
      if(debugfail>=2 && showmessages){
        dualfprintf(fail_file,"Tried entropy and bad on finalstep=%d! t=%g steppart=%d nstep=%ld i=%d j=%d k=%d :: hotpflag=%d entropypflag=%d\n",finalstep,t,steppart,nstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,hotpflag,entropypflag);
        //        int pliter;
        //        if(pr[UU]>1.0) dualfprintf(fail_file,"ADEATH\n");
        //        if(pr[UU]<1E-13) dualfprintf(fail_file,"ADEATH2\n");
        //        if(Ugeomfree0[ENTROPY]<-10) dualfprintf(fail_file,"ADEATH3\n");
        //        PLOOP(pliter,pl) dualfprintf(fail_file,"POST DID NOT GO TO ENTROPY: ijk=%d %d %d nstep=%ld steppart=%d pl=%d pr0=%g pr=%g U0=%g U=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,pl,pr0[pl],pr[pl],Ugeomfree0[pl],Ugeomfree[pl]);
      }

#endif

    }
        
  }// if bad failure of some kind





  return(0);
}








// use cold grmhd if hot grmhd gives u<zerouuperbaryon*rho since then cold is valid
#define USECOLDIFHOTUNEG 1
// use cold grmhd if hot grmhd gives rho<0 since rho<0 is nearly no solution
#define USECOLDIFHOTRHONEG 1
// use cold grmhd if hot grmhd leads to convergence or other "no solution" type failure
#define USECOLDIFHOTFAILCONV 1

// try cold inversion of hot one fails
int trycoldinversion(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, PFTYPE hotpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, FTYPE *Ugeomfree, FTYPE *Ugeomfree0, struct of_geom *ptrgeom, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int pl;
  FTYPE prhot[NPR],prcold[NPR];
  PFTYPE coldpflag;
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);


  //  dualfprintf(fail_file,"Got here in trycoldinversion\n");

  if( (USECOLDIFHOTRHONEG==0 && IFUTOPRIMFAILSOFTRHORELATED(hotpflag)) || (USECOLDIFHOTUNEG==0 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) ){
    // then maybe not so bad failure
    // e.g. if get here, do nothing when either rho<=0 or u<=0
  }
  else if( (USECOLDIFHOTRHONEG==1 && IFUTOPRIMFAILSOFTRHORELATED(hotpflag)) || (USECOLDIFHOTUNEG==1 && IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)) || (USECOLDIFHOTFAILCONV==1 && IFUTOPRIMFAIL(hotpflag)) ){
    // get here if want to fix up rho<=0, u<=0, or convergence failure
    // First if() is needed since last conditional or else if() would trigger on any failure type
    

    // then bad failure, so try to use cold grmhd
    // restore backup in case previous inversion changed things
    PALLLOOP(pl){
      prhot[pl]=pr[pl];
      Ugeomfree[pl]=Ugeomfree0[pl];

      // set guess and sets any pre-inverted quantities
      prcold[pl]=pr0[pl];
    }
    
    
    // get cold inversion
    Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, EOMCOLDGRMHD, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &coldpflag, prcold,pressure,newtonstats, lpflagrad);

    
    ///////////////////////////////////
    //
    // check if cold solution is good
    if(IFUTOPRIMNOFAILORFIXED(coldpflag)){

      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMFAILFIXEDCOLD; // default then is that no failure (or fixed failure)

      // then keep cold solution
      PALLLOOP(pl) pr[pl]=prcold[pl];


      // but set internal energy to previous value (should really evolve with entropy equation, but if negligible and no strong shocks, then ok )
      // GODMARK: another ok way is to set special failure that only averages internal energy!  Then can evolve at least -- in some diffusive way
      
#if(PRODUCTION==0)      
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and good on finalstep=%d! i=%d j=%d k=%d :: hotpflag=%d coldpflag=%d\n",finalstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,hotpflag,coldpflag);
#endif

      ///////////////////////////////
      //
      // decide how to use cold inversion solution    
      if(0&&IFUTOPRIMFAILSOFTNOTRHORELATED(hotpflag)){ // avoided now since reverting to 0 can introduce extra structure.  Want to keep to using averaging if possible that will generate better solution.
        // since cold approximation is very good, then use cold solution and just set internal energy to 0
        // if internal energy is actually small, then just set it to 0
        // works for Hubble flow!
        pr[UU]=zerouuperbaryon*pr[RHO];
      }
      else{
        //////////////
        //  if internal energy is not negligible or unknown, then should average or evolve!
        pr[UU]=pr0[UU];
        // then very bad failure, so try cold inversion and average internal energy for now
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMFAILU2AVG1FROMCOLD; // assume only internal energy needs correcting by averaging
 
      }// end if (else) trying cold inversion





      // accounting (since fixup.c accounting doesn't know what original pr0 is and actual prhot is undefined since Uhot->prhot wasn't possible)
      // GODMARK: For DOENOFLUX>0, should modify conserved quantity in some way.  For DOENOFLUX==0, primitives form basis of conserved quantities, so once primitives are modified all is done.  So should probably pass in U[] from  Utoprimgen() or whatever is the full conserved quantity later used.
      int modcons=0;
      FTYPE Ui[NPR];
      // ucons not modified (i.e. modcons=0), but ucons may be used by diag_fixup()
      UtoU(UNOTHING,UDIAG,ptrgeom,Ugeomfree0,Ui);
      // account for change to hot MHD conserved quantities
      int counttype;
      if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)==UTOPRIMFAILFIXEDCOLD) counttype=COUNTCOLD;
      else counttype=COUNTNOTHING; // do correction, but don't do counting until later
      diag_fixup_Ui_pf(modcons,Ui,pr,ptrgeom,finalstep,counttype);

      // KORALTODO: allowlocalfailurefixandnoreport could be used below to avoid reset, but for now let implicit solver be ok with cold inversion
      // reset pflag since above does full accounting, unless need to average-out internal energy still
      if(GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)==UTOPRIMFAILFIXEDCOLD){
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
      }




    }// else if coldpflag is bad
    else{
      // then both hot and cold are bad, so keep hot
      // GODMARK: Could try going to force-free and "failing" the parallel velocity so it gets averaged like internal energy in cold case!
      // only can go to force-free if b^2/rho>>1 as well
      // keep hotpflag and keep hot solution
      PALLLOOP(pl) pr[pl]=prhot[pl];

#if(PRODUCTION==0)      
      if(debugfail>=2) dualfprintf(fail_file,"Tried cold and bad on finalstep=%d! t=%g steppart=%d nstep=%ld i=%d j=%d k=%d :: hotpflag=%d coldpflag=%d\n",finalstep,t,steppart,nstep,ptrgeom->i,ptrgeom->j,ptrgeom->k,hotpflag,coldpflag);
#endif

    }
        
  }// if bad failure of some kind

  return(0);
}










// Function to check on success of inversion by comparing original conserved quantities (Uold) with new conserved quantities (Unew(p(Uold))) as computed from primitives (p(Uold)) obtained from the inversion.
// Changes failure flag if bad inversion is bad enough
// Note that inversion might imply small error, but still very wrong.  This occurs when no solution but approach solution from one side.

// Like force-free that doesn't limit |vrad|<c, radiation primitives might not correspond to conserved because of such limiters in u2p_rad(), so should generally have this turned off and only use for debugging

static int check_on_inversion(int usedhotinversion,int usedentropyinversion,int usedcoldinversion,int usedffdeinversion, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int badinversion,badinversionfail;
  int badinversionrad,badinversionfailrad;
  FTYPE Unormalnew[NPR],Unormalold[NPR];
  FTYPE fdiff[NPR];
  struct of_state q;
  int pl,pliter;
  int j,k,jj;
  FTYPE errornorm;


  

  // see if should check at all
  if(CHECKONINVERSION==0 && CHECKONINVERSIONRAD==0) return(0);



  ////////////
  //
  // Store extra quantities to enforce consistency when using tabulated EOS
  //
  ////////////

  if(WHICHEOS==KAZFULL && ptrgeom->p==CENT){
    // then store pressure to be used by get_stateforcheckinversion()
    // assume standard inversion at loc=CENT
    GLOBALMACP0A1(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k,PGASGLOBAL)=(*pressure);
  }

  

  // assume not bad
  badinversion=badinversionfail=0;
  badinversionrad=badinversionfailrad=0;



  //  if(1 || !(*lpflag)){ // only report if not failure
  if(!(*lpflag)){ // only report if not failure

    /////////////
    //
    // get new U (only applicable if hot inversion used)
    //
    /////////////
    MYFUN(get_stateforcheckinversion(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included


    /////////////
    //
    // Create a normalized (geometrically) conserved quantity
    //
    /////////////
    // assume default value assignment for normalized version of U
    PLOOP(pliter,pl){
      Unormalold[pl] = Uold[pl];
      Unormalnew[pl] = Unew[pl];
    }      

#if(REMOVERESTMASSFROMUU==2)
    // now make U_0 and U_i comparable
    // like \rho\gamma^2 v^2/(1+\gamma) + (u+p)\gamma^2 \sim \rho \gamma v^2 + (u+p)\gamma^2
    // In the limit as v->0, this can be compared against the normalized version of the U_i terms
    Unormalold[UU] = Uold[UU];
    Unormalnew[UU] = Unew[UU];

    // No longer doing below, just use gcon or gcov to get correct dimensional comparison
    // now deal with momentum terms and correct for geometry so all U_i's are comparable
    // (\rho+u+p)\gamma u_i -> (\rho+u+p)\gamma^2 v^2
    //    SLOOPA(jj) Unormalold[UU+jj] = Uold[UU+jj]*(q.ucon[jj]/q.ucon[TT]);
    //    SLOOPA(jj) Unormalnew[UU+jj] = Unew[UU+jj]*(q.ucon[jj]/q.ucon[TT]);
#endif



    //////////////////////////////
    //
    // now check for errors
    //
    //////////////////////////////
    PLOOP(pliter,pl){

      if(CHECKONINVERSION==0 && (pl!=PRAD0 && pl!=PRAD1 && pl!=PRAD2 && pl!=PRAD3)) continue; // don't check MHD (non-radiation) inversion if didn't want to
      if(CHECKONINVERSIONRAD==0 && (pl==PRAD0 && pl==PRAD1 && pl==PRAD2 && pl==PRAD3)) continue; // don't check radiation (non-MHD) inversion if didn't want to
  

      // default is to assume nothing wrong
      fdiff[pl]=0.0;


      // pl here is conserved quantity, so checking if conserved quantity not same as used to get primitive
      // in force-free projection on velocity always occurs that makes momenta terms not the same
      // no point checking if inversion doesn't handle or is inconsistent with conservation of that quantity
      if(usedffdeinversion && (pl==RHO || pl==UU || pl==U1 || pl==U2 || pl==U3 || pl==ENTROPY || pl==YNU || pl==YL) ) continue; // if ffde, no mass or thermal/hot components
      if(usedcoldinversion  && (pl==UU || pl==ENTROPY || pl==YNU || pl==YL) ) continue; // if cold, can't evolve any hot or thermal component
      // inversion either uses energy or entropy and can't use both at once inside inversion routine
      if(usedentropyinversion  && (pl==UU) ) continue; // entropy doesn't use energy equation, but does use conserved entropy
      if(usedhotinversion && (pl==ENTROPY) ) continue; // hot doesn't use entropy, but does use conserved energy

      //(EOMRADTYPE==EOMRADNONE && (pl==URAD0 || pl==URAD1 || pl==URAD2 || pl==URAD3)) || // no need to check pl if no such pl's
      // lpflagrad: Checks that if u2p placed limiter on p (e.g. velocity), then should skip this check since won't be accurate inversion
      if(EOMRADTYPE!=EOMRADNONE && (*lpflagrad==1)) continue;
      // If doing Eddington approximation, actually ignore conserved flux evolution.
      if(EOMRADTYPE==EOMRADEDD && (pl==URAD1 || pl==URAD2 || pl==URAD3)) continue;


      if(pl==YNU || pl==YL){
        // only check that new U is finite
        if(!isfinite(Unormalnew[pl])){
          fdiff[pl]=BIG; // indicates was not finite
        }
        else continue; // then avoid checking passive scalars in case manipulated
      }

      // leave geometry out of it
      //      Unormalnew[pl]*=ptrgeom->gdet;
      //      Unormalold[pl]*=ptrgeom->gdet;
      if(pl==RHO || pl==UU || pl==URAD0&&(*lpflagrad==0)){
        errornorm=(fabs(Unormalnew[pl])+fabs(Unormalold[pl])+SMALL);
        // KORALTODO: when URAD0<<UU, can't expect radiation error to be small relative to only itself when interaction between radiation and fluid. 
        if(pl==URAD0) errornorm+=(fabs(Unormalnew[UU])+fabs(Unormalold[UU])+SMALL);
        fdiff[pl] = fabs(Unormalnew[pl]-Unormalold[pl])/errornorm;
      }
      else if(pl==ENTROPY){// can be + or -, so use positive definite exp(S/rho)=exp(s) version
        errornorm=(fabs(exp(Unormalnew[pl]/(SMALL+fabs(Unormalnew[RHO]))))+fabs(exp(Unormalold[pl]/(SMALL+fabs(Unormalold[RHO]))))+SMALL);
        fdiff[pl] = fabs(exp(Unormalnew[pl]/(SMALL+fabs(Unormalnew[RHO])))-exp(Unormalold[pl]/(SMALL+fabs(Unormalold[RHO]))))/errornorm;
      }
      else if(pl==U1 || pl==U2 || pl==U3){

        errornorm  = THIRD*(fabs(Unormalnew[U1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)]))) + fabs(Unormalold[U1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)]))) + fabs(Unormalnew[U2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalold[U2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalnew[U3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)])))+fabs(Unormalold[U3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)]))));
#if(REMOVERESTMASSFROMUU==2)
        // only valid comparison if rest-mass taken out of energy term and modify U_i term to be comparable with U_t term
        errornorm = MAX(errornorm,0.5*(fabs(Unormalold[UU])+fabs(Unormalnew[UU])));
#endif
     
        fdiff[pl] = sqrt(fabs(ptrgeom->gcon[GIND(pl-UU,pl-UU)]))*fabs(Unormalnew[pl]-Unormalold[pl]) / (errornorm+SMALL);

      }
      else if( (pl==URAD1 || pl==URAD2 || pl==URAD3)&&(*lpflagrad==0) ){
        //        errornorm  = THIRD*(fabs(Unormalnew[URAD1])+fabs(Unormalold[URAD1])+fabs(Unormalnew[URAD2])+fabs(Unormalold[URAD2])+fabs(Unormalnew[URAD3])+fabs(Unormalold[URAD3]));
        errornorm  = THIRD*(fabs(Unormalnew[URAD1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)])))+fabs(Unormalold[URAD1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)])))+fabs(Unormalnew[URAD2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalold[URAD2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalnew[URAD3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)])))+fabs(Unormalold[URAD3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)]))));
        // KORALTODO: when URAD0<<UU, can't expect radiation error to be small relative to only itself when interaction between radiation and fluid. 
        errornorm  += THIRD*(fabs(Unormalnew[U1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)])))+fabs(Unormalold[U1]*sqrt(fabs(ptrgeom->gcon[GIND(1,1)])))+fabs(Unormalnew[U2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalold[U2]*sqrt(fabs(ptrgeom->gcon[GIND(2,2)])))+fabs(Unormalnew[U3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)])))+fabs(Unormalold[U3]*sqrt(fabs(ptrgeom->gcon[GIND(3,3)]))));
        errornorm = MAX(errornorm,0.5*(fabs(Unormalold[URAD0])+fabs(Unormalnew[URAD0])));

        fdiff[pl] = sqrt(fabs(ptrgeom->gcon[GIND(pl-URAD0,pl-URAD0)]))*fabs(Unormalnew[pl]-Unormalold[pl]) / (errornorm+SMALL);
      }
      else if(pl==B1 || pl==B2 || pl==B3){

        errornorm  = THIRD*(fabs(Unormalnew[B1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)])))+fabs(Unormalold[B1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)])))+fabs(Unormalnew[B2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)])))+fabs(Unormalold[B2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)])))+fabs(Unormalnew[B3]*sqrt(fabs(ptrgeom->gcov[GIND(3,3)])))+fabs(Unormalold[B3]*sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))));
        fdiff[pl] = sqrt(fabs(ptrgeom->gcov[GIND(pl-B1+1,pl-B1+1)]))*fabs(Unormalnew[pl]-Unormalold[pl]) / (errornorm+KINDASMALL); // for field order 1E-100 or less, geometry division and multiplication causes this to shift in value by order unity.
      }
    }


    
    // broke loop to check multiple directions


    int plcheck;
    PLOOP(pliter,pl){

      if(CHECKONINVERSION==0 && (pl!=PRAD0 && pl!=PRAD1 && pl!=PRAD2 && pl!=PRAD3)) continue; // don't check MHD (non-radiation) inversion if didn't want to
      if(CHECKONINVERSIONRAD==0 && (pl==PRAD0 && pl==PRAD1 && pl==PRAD2 && pl==PRAD3)) continue; // don't check radiation (non-MHD) inversion if didn't want to


      plcheck=(pl>=RHO)&&(pl<=B3 || pl<=ENTROPY && usedentropyinversion || (*lpflagrad==0)&&(EOMRADTYPE!=EOMRADNONE && (pl==URAD0 || pl==URAD1&&(EOMRADTYPE!=EOMRADEDD) || pl==URAD2&&(EOMRADTYPE!=EOMRADEDD) || pl==URAD3&&(EOMRADTYPE!=EOMRADEDD) )));

      if(IFUTOPRIMFAIL(*lpflag) || fdiff[pl]>CHECKONINVFRAC){
        if(
           ( (plcheck)&&((fabs(Unormalold[pl])>SMALL)&&(fabs(Unormalnew[pl])>SMALL)) )
           ){
          if(pl<URAD0 && pl>URAD3) badinversion++;
          else  badinversionrad++;

          if(pl==ENTROPY) dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g : %21.15g %21.15g : %21.15g %21.15g\n",pl,fdiff[pl],Unormalold[pl],Unormalnew[pl],exp(Unormalold[pl]/Unormalold[RHO]),exp(Unormalnew[pl]/Unormalnew[RHO]),Unormalold[RHO],Unormalnew[RHO]);
          else dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",pl,fdiff[pl],Unormalold[pl],Unormalnew[pl]);

        }
      }
      if(fdiff[pl]>CHECKONINVFRACFAIL){
        if( (plcheck)&&((fabs(Unormalold[pl])>SMALL)&&(fabs(Unormalnew[pl])>SMALL)) ){
          if(pl<URAD0 && pl>URAD3) badinversionfail++;
          else  badinversionfailrad++;
        }
      }


      //      dualfprintf(fail_file,"pl=%d Uold=%26.20g Unew=%26.20g\n",pl,Unormalold[pl],Unormalnew[pl]);


    }

    ////////////
    //
    // change failure flag if really bad check
    //
    /////////////
    if(CHECKONINVERSION==1 && FAILIFBADCHECK && badinversionfail){
      if(debugfail>=2) dualfprintf(fail_file,"Changing flag since sufficiently bad inversion.\n");
      (*lpflag)=UTOPRIMFAILCONVBADINVERTCOMPARE;
    }
    // CHECKONINVERSIONRAD==1 case: No, would not redo inversion since no reduction to another inversion.



    ////////////
    //
    // Account for any conserved quantity change (could be just difference created by primtoU, but still tells order of issue)
    //
    // only makes sense to account for changes in U if inversion is treated as success
    //
    // Also, only makes sense if using primitives as basis.  Since, if use U as basis, then each iteration relies upon old correct U without any error introduced due to inversion.
    //
    /////////////



    ////////////
    //
    // Report bad inversion
    //
    /////////////

    //    if(myid==5 && nstep==1 && steppart==0 && ptrgeom->i==19 && ptrgeom->j==15){
    //      badinversion=1;
    //    }
    
    if(badinversion || badinversionrad){
      dualfprintf(fail_file,"Bad inversion (or possibly Bad U(p) calculation):\n");
      dualfprintf(fail_file,"Inversion types: %d %d %d %d\n",usedffdeinversion,usedcoldinversion,usedentropyinversion,usedhotinversion);

      dualfprintf(fail_file,"t=%21.15g nstep=%ld stepart=%d :: i=%d j=%d k=%d :: lntries=%d lerrx=%21.15g\n",t,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,newtonstats->lntries,newtonstats->lerrx);
      PLOOP(pliter,pl) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Unewgeomfree[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g) fdiff=%21.15g\n",pl,Uold[pl],pl,Unew[pl],pl,pr[pl],pl,pr0[pl],fdiff[pl]);
      dualfprintf(fail_file,"special inversion pressure=%21.15g statepressure=%21.15g stateentropy=%21.15g\n",*pressure,q.pressure,q.entropy);
      
      dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->gdet);
      
      DLOOP(j,k) dualfprintf(fail_file,"gcov=%21.15g gcon=%21.15g\n",ptrgeom->gcov[GIND(j,k)],ptrgeom->gcon[GIND(j,k)]);

      DLOOPA(k) dualfprintf(fail_file,"k=%d : q.ucon=%21.15g q.ucov=%21.15g :  q.uradcon=%21.15g q.uradcov=%21.15g  : q.bcon=%21.15g q.bcov=%21.15g\n",k,q.ucon[k],q.ucov[k],q.uradcon[k],q.uradcov[k],q.bcon[k],q.bcov[k]);

      // only really need the below quantities to check on inversion in mathematica
      // Use Eprime_inversion.nb to check on utoprim_jon.c inversion
      for(j=0;j<NUMINVPROPERTY;j++) dualfprintf(fail_file,"%sstr=\"%21.15g\";\n",newtonstats->invpropertytext[j],newtonstats->invproperty[j]);
    }

  }

  //  if(ptrgeom->i==39){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d : pl=%d Uold=%21.15g Unew=%21.15g pr0=%21.15g pr=%21.15g\n",ptrgeom->i,pl,Uold[pl],Unew[pl],pr0[pl],pr[pl]);
  //  }


  return(0);


}










      






// compare ffde inversions
static int compare_ffde_inversions(int showmessages, int allowlocalfailurefixandnoreport, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_geom *ptrgeom, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, FTYPE *Uold, FTYPE *Unew, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int j,k;
  int pl;
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  struct of_state q;
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int eomtypelocal=EOMTYPE; // as chosen before introduced eomtypelocal


  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMFFDE, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), pr, pressure, newtonstats,lpflagrad);

  PALLLOOP(pl) Ugeomfree[pl]=Ugeomfree0[pl]; // make sure good conserved quantity
      
  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtypelocal, EVOLVENOENTROPY, Ugeomfree, ptrgeom, &GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL), prother, pressure, newtonstats, lpflagrad);


  // now get conserved quantity from pr to check
  // find U(p)
  MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Upr),"step_ch.c:advance()", "primtoU()", 1);

  MYFUN(get_state(prother, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  MYFUN(primtoU(UNOTHING,prother, &q, ptrgeom, Uprother),"step_ch.c:advance()", "primtoU()", 1);


  // now compare

  //PALLLOOP(pl)
  dualfprintf(fail_file,"nstep=%ld stepart=%d i=%d j=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j);
  for(pl=U1;pl<B3;pl++){
    /*
      if(
      ((fabs(Ugeomfree[pl]-Upr[pl])/(fabs(Ugeomfree[pl])+fabs(Upr[pl])+SMALL)) > 1E-12)||
      ((fabs(Ugeomfree[pl]-Uprother[pl])/(fabs(Ugeomfree[pl])+fabs(Uprother[pl])+SMALL)) > 1E-12)
      ){
      dualfprintf(fail_file,"DIFF: %21.15g %21.15g ::pl=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(Ugeomfree[pl]-Upr[pl])/(fabs(Ugeomfree[pl])+fabs(Upr[pl])+SMALL)),(fabs(Ugeomfree[pl]-Uprother[pl])/(fabs(Ugeomfree[pl])+fabs(Uprother[pl])+SMALL)),pl,Ugeomfree[pl],pr[pl],prother[pl],Upr[pl],Uprother[pl]);
      }
    */
    if(
       ((fabs(pr[pl]-prother[pl])/(fabs(pr[pl])+fabs(prother[pl])+SMALL)) > 1E-12)&&( (fabs(pr[pl])>1E-15)||(fabs(prother[pl])>1E-15) )
       ){
      dualfprintf(fail_file,"DIFF: %21.15g :: pl=%d Ugeomfree=%21.15g prold=%21.15g prnew=%21.15g Upr=%21.15g Uprother=%21.15g\n",(fabs(pr[pl]-prother[pl])/(fabs(pr[pl])+fabs(prother[pl])+SMALL)),pl,Ugeomfree[pl],pr[pl],prother[pl],Upr[pl],Uprother[pl]);
    }
  }


  return(0);

}




// arbitrary debugging code created at any time
static int debug_utoprimgen(PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Uold, FTYPE *Unew)
{
  int j,k;
  int pl;
  FTYPE fdiff[NPR];






#define CRAZYDEBUG 0
#if(CRAZYDEBUG) // begin crazy debug stuff
  if(1|| newtonstats->lntries>3){
    if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){
      //    if(nstep==0 && steppart==0 && ptrgeom->i==4 && ptrgeom->j==36){
      //  if(nstep==0 && steppart==0 && ptrgeom->i == 5 && ptrgeom->j == 31 ){
      dualfprintf(fail_file,"nstep=%ld stepart=%d :: i=%d j=%d :: lntries=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j,newtonstats->lntries);
    
      //    PALLLOOP(pl) dualfprintf(fail_file,"Ugeomfree[%d]=%21.15g pr[%d]=%21.15g\n",pl,Ugeomfree[pl],pl,pr[pl]);
      PALLLOOP(pl) dualfprintf(fail_file,"Uoldgeomfree[%d]=%21.15g Uold[%d]=%21.15g pr[%d]=%21.15g (pr0[%d]=%21.15g)\n",pl,Uold[pl],pl,Uold[pl]*ptrgeom->gdet,pl,pr[pl],pl,pr0[pl]);

      dualfprintf(fail_file,"g=%21.15g\n",ptrgeom->gdet);

      DLOOP(j,k) dualfprintf(fail_file,"gcon=%21.15g gcov=%21.15g\n",ptrgeom->gcov[GIND(j,k)],ptrgeom->gcon[GIND(j,k)]);
      //    myexit(777);
      //  }
    }
  }


  // test inversion for accuracy in conservative space
  //  if(nstep==0 && steppart==0 && ptrgeom->i == 4 && ptrgeom->j == 36 ){
  if(nstep==9 && steppart==2 && ptrgeom->i==0 && ptrgeom->j==31){

#if(0)
    pr[0]= 7.22714301361038e-06 ;
    pr[1]=-2.55753797775927e-07 ;
    pr[2]= 0.000892168434512681 ;
    pr[3]=  -0.0349621334882251 ;
    pr[4]=                   -0 ;
    pr[5]= -0.00101880685475446 ;
    pr[6]=   0.0399035382308458 ;
    pr[7]=                    0 ;
#elif(0)
    pr[0]= 7.46157819677347e-06;
    pr[1]=  7.3407120498644e-06;
    pr[2]= 0.000367464781319951;
    pr[3]=  -0.0144105143008605;
    pr[4]=                    0;
    pr[5]= -0.00101880685475446;
    pr[6]=   0.0399035382308458;
    pr[7]=                    0;
#elif(0)
    pr[0]=  7.5124289176258e-06 ;
    pr[1]= 1.33752037209996e-08 ;
    pr[2]= 1.33529579432262e-07 ;
    pr[3]=-5.23276639757274e-06 ;
    pr[4]=                    0 ;
    pr[5]= -0.00101880685475444 ;
    pr[6]=   0.0399035382308459 ;
    pr[7]=                    0 ;

#endif

    MYFUN(get_state(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included

    for(k=0;k<4;k++){
      dualfprintf(fail_file,"q.ucon[%d]=%21.15g q.ucov[%d]=%21.15g q.bcon[%d]=%21.15g q.bcov[%d]=%21.15g\n",k,q.ucon[k],k,q.ucov[k],k,q.bcon[k],k,q.bcov[k]);
    }

    PLOOP(pliter,pl){
      Unew[pl]*=ptrgeom->gdet;
      Uold[pl]*=ptrgeom->gdet;
      fdiff[pl] = fabs(Unew[pl]-Uold[pl])/(fabs(Unew[pl]+Uold[pl])+SMALL);
      //    if(fdiff[pl]>1E-10){
      if((pl>=RHO)&&(pl<=B3)&&((fabs(Uold[pl])>1E-20)||(fabs(Unew[pl])>1E-20))){
        dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %21.15g %21.15g\n",pl,fdiff[pl],Uold[pl],Unew[pl]);
      }
      //    }
    }

    dualfprintf(fail_file,"dt=%21.15g\n",dt);


    myexit(124);
  }
#endif// end crazy debug stuff
#undef CRAZYDEBUG


  return(0);


}




// deal with negative density in special way to tell fixup routine to perform special "average" of density
static int negdensitycheck(int finalstep, FTYPE *prim, PFTYPE *pflag)
{


  //Inversion from the average value succeeded or has a negative density or internal energy
  if(IFUTOPRIMFAILSOFT(*pflag)) {

    if( prim[UU] < zerouuperbaryon*prim[RHO] && (finalstep&&STEPOVERNEGU==NEGDENSITY_FIXONFULLSTEP) || STEPOVERNEGU==NEGDENSITY_ALWAYSFIXUP) {
      *pflag = UTOPRIMFAILU2AVG2;
    }
  }

  return(0);
}





// perform U->p inversion on advected scalars
int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr)
{
  FTYPE myrhouu0,oneOmyrhouu0;
  int i,j,k,loc;
  FTYPE prforadvect;
  FTYPE ylforadvect,ynuforadvect;



  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;
  


  /////////////
  //
  // avoid division by 0 but allow sign and 1/0 ->0 assumed 0 geometry  means value can be taken to be 0
  //
  /////////////
  myrhouu0=Ugeomfree[RHO];
  oneOmyrhouu0=sign(Ugeomfree[RHO])/(fabs(Ugeomfree[RHO])+SMALL);




  ///////////////
  //
  // Invert U->direct Primitive for scalars
  //
  ///////////////

#if(DOYL!=DONOYL)
  ylforadvect = Ugeomfree[YL]*oneOmyrhouu0;
#else
  ylforadvect=0.0;
#endif

#if(DOYNU!=DONOYNU)
  ynuforadvect = Ugeomfree[YNU]*oneOmyrhouu0;
#else
  ynuforadvect=0.0;
#endif

  ///////////////
  //
  // Invert U->P for scalars
  //
  ///////////////


#if(DOYL!=DONOYL)
#if(WHICHEOS==KAZFULL)
  advect2yl_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),ylforadvect,ynuforadvect,&pr[YE]); // pr[YE] is pr[YL] memory space
#else
  pr[YL] = ylforadvect;
#endif
#endif

#if(DOYNU!=DONOYNU)
#if(WHICHEOS==KAZFULL)
  advect2ynu_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),ylforadvect,ynuforadvect,&pr[YNU]);
#else
  pr[YNU] = ynuforadvect;
#endif
#endif



  //////////////
  //
  // Change the primitives to be constrained or fixed-up.  Also perform any extra operations required by the EOS used.
  //
  /////////////
  fix_primitive_eos_scalars_simple(i,j,k,loc,pr);



  //////////////////
  //
  // Adjust conserved quantities based upon fixed-up primitives
  //
  //////////////////

#if(DOYL!=DONOYL)
#if(WHICHEOS==KAZFULL)
  yl2advect_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),pr[YL],pr[YNU],&prforadvect);
#else
  prforadvect = pr[YL];
#endif
  Ugeomfree[YL] = prforadvect*myrhouu0;
#endif

#if(DOYNU!=DONOYNU)
#if(WHICHEOS==KAZFULL)
  ynu2advect_kazfull(GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),pr[YL],pr[YNU],&prforadvect);
#else
  prforadvect = pr[YNU];
#endif
  Ugeomfree[YNU] = prforadvect*myrhouu0;
#endif


  // SUPER TODO: Consider if need to recompute KAZ EOSextra stuff for fluxes.  Maybe now with primitives as lookups, interpolations correct.  But need to ensure assignments to EOSextra for Y_e and Ynu0 is correctly done.  Can leave neutrino stuff as "DONOR" cell

  return(0);
}





// convert U (conserved quantities) to form acceptable to utoprimversion for a given removerestmassfromuu
void convert_U_removerestmassfromuu(int utoprimversion, int removerestmassfromuu, FTYPE *U)
{

  ///////////////////////////////////////////////////
  //
  ///////////// Setup U[UU] with or without rest-mass for different inversions and setting of REMOVERESTMASSFROMUU
  //
  // check if inversion is chosen consistent with REMOVERESTMASSFROMUU
  // At this point, Ugeomfree[UU] in "standard form with rest-mass" for REMOVERESTMASSFROMUU=0,1.  If ==2, then no rest-mass
  //
  ///////////////////////////////////////////////////

  if(removerestmassfromuu==2){
    // 5D1 handles removerestmassfromuu internally
    // NONRELCOMPAT only accepts UU without rest-mass, so good to go
    if( utoprimversion!=UTOPRIM5D1 && (utoprimversion!=UTOPRIMJONNONRELCOMPAT) ){
      //    dualfprintf(fail_file,"This code does not handle removerestmassfromuu==2: other Utoprim's besides 5D1 and 2DFINAL\n");
      //    myexit(1);
      // then change conserved quantity so can work
      U[UU]-=U[RHO]; // "add in rest-mass" so in correct form
    }
  }
  else if(removerestmassfromuu==0 || removerestmassfromuu==1){
    if(utoprimversion==UTOPRIMJONNONRELCOMPAT){
      // then change conserved quantity so can work for UTOPRIMJONNONRELCOMPAT
      U[UU]+=U[RHO]; // "subtract out rest-mass" so in correct form for NONRELCOMPAT
    }
  }


}










// Used for dissipation calculation only
int Utoprimdiss(int showmessages, int allowlocalfailurefixandnoreport, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, PFTYPE *otherpflag, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  FTYPE Upr[NPR],Uprother[NPR];
  int whichentropy;
  struct of_state q;
  FTYPE Uold[NPR],Unew[NPR];
  FTYPE fdiff[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int pl,pliter;
  FTYPE pressuremem;
  FTYPE *pressure=&pressuremem;
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int eomtypelocal=EOMTYPE; // as chosen before introduced eomtypelocal





  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);


  /////////////////////////////////////////////////////////
  //
  // backup
  //
  ////////////////////////////////////////////////////////
  PALLLOOP(pl){
    Uold[pl]=Ugeomfree0[pl]=Ugeomfree[pl];
    pr0[pl]=pr[pl];
  }
  

  /////////////////////////////////////////////////////////////
  //
  // convert U into form appropriate for inversion routine (probably does nothing since assume UTOPRIM5D1 used) and feed Utprim Uold anyways
  //
  /////////////////////////////////////////////////////////////
  convert_U_removerestmassfromuu(UTOPRIM5D1,REMOVERESTMASSFROMUU,Ugeomfree);
  convert_U_removerestmassfromuu(UTOPRIM5D1,REMOVERESTMASSFROMUU,Ugeomfree0);


  ////////////////////////
  //
  // DO INVERSION

  // assume solution is good unless flagged as bad
  *otherpflag=UTOPRIMNOFAIL;

  // do separate inversion for entropy version of EOMs
  // only one inversion is setup to handle this
  whichentropy=WHICHENTROPYEVOLVE;

  ////////////////////////
  // get entropy evolution (don't use failure -- otherpflag)
  // this inversion knows about global settings for DOENTROPY and increases tolerance for doing comparison
  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtypelocal, whichentropy, Uold, ptrgeom, otherpflag, pr, pressure, newtonstats, lpflagrad);

  // DEBUG:
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"otherpflag=%d Uold[%d]=%21.15g pr0[%d]=%21.15g pr[%d]=%21.15g\n",*otherpflag,pl,Uold[pl],pl,pr0[pl], pl,pr[pl]);

  
  // now get conserved quantity from pr to check
  // find U(p)
  //  MYFUN(get_state(pr, ptrgeom, &q),"step_ch.c:advance()", "get_state()", 1);
  //  MYFUN(primtoU(UNOTHING,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1);
  
  
  //  if(0&& *otherpflag){
  //    PLOOP(pliter,pl){
  //  fdiff[pl] = fabs(Unew[pl]-Uold[pl])/(fabs(Unew[pl]+Uold[pl])+1E-30);
  //      if(fdiff[pl]>1E-10){
  // if((pl>=U1)&&(pl<=B3)&&((fabs(Uold[pl])>1E-20)||(fabs(Unew[pl])>1E-20))){
  //   dualfprintf(fail_file,"fdiff[%d]=%21.15g :: %g %g\n",pl,fdiff[pl],Uold[pl],Unew[pl]);
  // }
  //      }
  // }
  //      }

     
  return(0);
}






//////////////////////////////
//
// COMPARISON of 2 (currently only 2) methods
int Utoprimgen_compare(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  int pl,pliter;
  FTYPE Uo1[NPR], Uo2[NPR];
  FTYPE ptest1[NPR],ptest2[NPR];
  FTYPE Uo[NPR], po[NPR];
  int test;
  PFTYPE lpflag1,lpflag2;
  int lntries1,lntries2;
  FTYPE lerrx1,lerrx2;
  FTYPE pressure1mem,pressure2mem;
  FTYPE *pressure1=&pressure1mem,*pressure2=&pressure2mem;

  // backup
  PALLLOOP(pl){
    ptest1[pl]=ptest2[pl]=pr[pl];
    Uo1[pl]=Uo2[pl]=Ugeomfree[pl];
  }

  // currently comparing 5D1 and LDZ
  //  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtype, EVOLVENOENTROPY, Uo1, ptrgeom, lpflag, ptest1,newtonstats, lpflagrad);
  //  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMLDZ, eomtype, EVOLVENOENTROPY, Uo2, ptrgeom, lpflag, ptest2,newtonstats, lpflagrad);

  // remove rest-mass
  Uo1[UU]+=Uo1[RHO];
  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtype, EVOLVENOENTROPY, Uo1, ptrgeom, &lpflag1, ptest1,pressure1,newtonstats,lpflagrad);
  lntries1=newtonstats->lntries; newtonstats->lntries=0; lerrx1=newtonstats->lerrx;
  Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, UTOPRIM2DFINAL, eomtype, EVOLVENOENTROPY, Uo2, ptrgeom, &lpflag2, ptest2,pressure2,newtonstats,lpflagrad);
  lntries2=newtonstats->lntries; lerrx2=newtonstats->lerrx;

#define ERRORCONST (1E-10)

  PALLLOOP(pl){
    if(ptest1[pl]!=0.0) test=fabs(ptest1[pl]-ptest2[pl])/ptest1[pl]>ERRORCONST;
    else test=(ptest1[pl]-ptest2[pl])>ERRORCONST;
    if(test){
      dualfprintf(fail_file,"utoprimdiff: %d %ld :: %d %d %d :: %d ::  %21.15g   %21.15g   %21.15g %21.15g :: %d %d :: %21.15g %21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,pl,ptest1[pl]-ptest2[pl],(ptest1[pl]!=0.0) ? fabs(ptest1[pl]-ptest2[pl])/ptest1[pl] : ptest1[pl]-ptest2[pl],ptest1[pl],ptest2[pl],lntries1,lntries2,lerrx1,lerrx2);
    }
  }
  if(IFUTOPRIMFAIL(lpflag1) || IFUTOPRIMFAIL(lpflag2)){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2);
  }
  if(lntries1>10 || lntries2>10){
    dualfprintf(fail_file,"%d %ld :: %d %d %d :: lpflag1=%d lpflag2=%d errx1=%21.15g errx2=%21.15g lntries1=%d lntries2=%d\n",steppart,nstep,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,startpos[3]+ptrgeom->k,lpflag1,lpflag2,lerrx1,lerrx2,lntries1,lntries2);
  }

  
  // always do (use old utoprim)
  if(IFUTOPRIMNOFAILORFIXED(lpflag1)) PALLLOOP(pl){
      pr[pl]=ptest1[pl];
      *lpflag=lpflag1;
    }
  else{
    PALLLOOP(pl) pr[pl]=ptest2[pl];
    *lpflag=lpflag2;
  }

  

  return(0);
}




// just picks the algorithm to invert
int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats,PFTYPE *lpflagrad)
{
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  extern int Utoprim_coldgrmhd(FTYPE *U, struct of_geom *geom, FTYPE *pr, int *positivityproblem);
  int positivityproblem;
  int pliter,pl;



  //  FTYPE Ugeomfree0[NPR];
  //  if(ptrgeom->i==10 && ptrgeom->k==0){
  //    PLOOP(pliter,pl) Ugeomfree0[pl]=Ugeomfree[pl];
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"before inversion: pl=%d pr=%g\n",pl,pr[pl]);
  //    if(pr[UU]>0.1) dualfprintf(fail_file,"beforeADEATHUU: ijk=%d %d %d u=%g %g (%g)\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[UU],Ugeomfree[UU],Ugeomfree0[UU]);
  //    if(fabs(pr[U1])>1.0) dualfprintf(fail_file,"beforeADEATHU1: ijk=%d %d %d u=%g %g (%g)\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[U1],Ugeomfree[U1],Ugeomfree0[U1]);
  //  }


  /////////////
  //
  // do fluid inversion
  //
  /////////////

  // below is only one that uses parameter to choose whether entropy inversion or not
  if(which==UTOPRIM5D1){      MYFUN(Utoprim(parameter,Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim", 1);}
  else if(which==UTOPRIMLDZ){ MYFUN(Utoprim_ldz(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_ldz", 1);}
  else if(which==UTOPRIM2D){ MYFUN(Utoprim_2d(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d", 1);}
  else if(which==UTOPRIM1D){ MYFUN(Utoprim_1d(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d", 1);}
  else if(which==UTOPRIM1DOPT){ MYFUN(Utoprim_1d_opt(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d_opt", 1);}
  else if(which==UTOPRIM1DFINAL){ MYFUN(Utoprim_1d_final(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_1d_final", 1);}
  else if(which==UTOPRIM2DFINAL){ MYFUN(Utoprim_2d_final(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d_final", 1);}
  // Below Jon's method handles full  hotMHD, entropy, and coldMHD via setting eomtype (doesn't use parameter)
  else if(which==UTOPRIMJONNONRELCOMPAT){ MYFUN(Utoprim_jon_nonrelcompat_inputnorestmass(showmessages,eomtype,GLOBALMAC(EOSextraglobal,ptrgeom->i,ptrgeom->j,ptrgeom->k),Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_2d_final_nonrelcompat_inputnorestmass", 1);}
  else if(which==UTOPRIM5D2){ MYFUN(Utoprim_5d2_final(Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats),"step_ch.c:Utoprimgen()", "Utoprim_5d2_final", 1);}
  // alt (not really working) inversion for coldmhd (alt compared to UTOPRIMJONNONRELCOMPAT)
  else if(which==UTOPRIMCOLDGRMHD){ MYFUN(Utoprim_coldgrmhd(Ugeomfree, ptrgeom, pr,&positivityproblem),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);}
  // alt force-free inversion (alt compared to UTOPRIMJONNONRELCOMPAT)
  else if(which==UTOPRIMFFDE){ MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);}
  else{
    dualfprintf(fail_file,"No such which=%d in Utoprimgen_pick()\n");
    myexit(1769384215);
  }


  //  if(ptrgeom->i==10 && ptrgeom->k==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"after inversion: pl=%d pr=%g\n",pl,pr[pl]);
  //    if(pr[UU]>0.1) dualfprintf(fail_file,"ADEATHUU: ijk=%d %d %d u=%g %g (%g)\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[UU],Ugeomfree[UU],Ugeomfree0[UU]);
  //    if(fabs(pr[U1])>1.0) dualfprintf(fail_file,"ADEATHU1: ijk=%d %d %d u=%g %g (%g)\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[U1],Ugeomfree[U1],Ugeomfree0[U1]);
  //  }

  /////////////
  //
  // now do other non-fluid inversions
  //
  /////////////

  // KORAL
  // NOTEMARK: u2p_rad() uses pr, which will have updated velocities in case radiation inversion wants to use fluid frame reduction.  But need to know if got good solution, so pass that flag to u2p_rad()
  if(EOMRADTYPE!=EOMRADNONE) u2p_rad(showmessages, allowlocalfailurefixandnoreport,Ugeomfree,pr,ptrgeom,lpflag,lpflagrad);
  //*lpflagrad=0; // test that check_on_inversion triggered where velocity limiter applies


  //  if(ptrgeom->i==10 && ptrgeom->k==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"after rad inversion: pl=%d pr=%g\n",pl,pr[pl]);
  //    if(pr[UU]>0.1) dualfprintf(fail_file,"BDEATHUU: ijk=%d %d %d u=%g %g (%g)\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[UU],Ugeomfree[UU],Ugeomfree0[UU]);
  //    if(fabs(pr[U1])>1.0) dualfprintf(fail_file,"BDEATHU1: ijk=%d %d %d u=%g %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,pr[U1],Ugeomfree[U1],Ugeomfree0[U1]);
  //  }
  
  return(0);
}




// Utoprimgen try again.  Can be whatever sequence or number of inversions
int Utoprimgen_tryagain(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int Utoprimgen_tryagain_substep(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);

  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM2DFINAL, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  // can't trust below really (causes some kind of superslow-down when doing torus problem with  MCSTEEP)
  //  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  // LDZ as normally used generates nan's
  //  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIMLDZ, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM1DFINAL, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM1DOPT, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D2, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM2D, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM1D, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);

  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}


// Utoprimgen try again.  Can be whatever sequence or number of inversions
// only those algorithms designed directly for REMOVERESTMASSFROMUU==2
int Utoprimgen_tryagain2(int showmessages, int allowlocalfailurefixandnoreport, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int Utoprimgen_tryagain_substep(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);
  void convert_U_removerestmassfromuu(int utoprimverison, int removerestmassfromuu, FTYPE *U);
  FTYPE Uorig0[NPR],Uorig[NPR];
  int pl,pliter;


  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIMJONNONRELCOMPAT, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM5D1, eomtype, parameter, Ugeomfree0, Ugeomfree, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);

  PALLLOOP(pl){
    Uorig0[pl]=Ugeomfree0[pl];
    Uorig[pl]=Ugeomfree[pl];
  }
  // assume at this point that U is such that setup for REMOVERESTMASSFROMUU==2 so that UU has rest-mass subtracted.
  // so need to add it back in for 2D FINAL method
  // convert U into form appropriate for inversion routine 
  convert_U_removerestmassfromuu(UTOPRIM2DFINAL,REMOVERESTMASSFROMUU,Uorig);
  convert_U_removerestmassfromuu(UTOPRIM2DFINAL,REMOVERESTMASSFROMUU,Uorig0);
  // now in form usable by UTOPRIM2DFINAL
  Utoprimgen_tryagain_substep(showmessages, allowlocalfailurefixandnoreport, UTOPRIM2DFINAL, eomtype, parameter, Uorig0, Uorig, ptrgeom, lpflag, pr0, pr,pressure,newtonstats,lpflagrad);


  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}





// need to keep this function of to date if going to add other "non-critical" failures.  Point is that some inversions can't handle non-rel case and might get positive results when shouldn't have.
int Utoprimgen_tryagain_substep(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr0, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad)
{
  int pl;
  int Utoprimgen_pick(int showmessages, int allowlocalfailurefixandnoreport, int which, int eomtype, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, PFTYPE *lpflag, FTYPE *pr, FTYPE *pressure, struct of_newtonstats *newtonstats, PFTYPE *lpflagrad);

  // if really a bad failure that don't want / can't handle, then try again
  if(
     (IFUTOPRIMFAIL(*lpflag))
     &&(!((IFUTOPRIMFAILSOFTNOTRHORELATED(*lpflag))&&(STEPOVERNEGU)))
     &&(!((*lpflag==UTOPRIMFAILRHONEG)&&(STEPOVERNEGRHO)))
     &&(!((*lpflag==UTOPRIMFAILRHOUNEG)&&(STEPOVERNEGRHOU)))
     ){
    // restore
    PALLLOOP(pl){
      Ugeomfree[pl]=Ugeomfree0[pl];
      pr[pl]=pr0[pl];
    }
    // try again
    MYFUN(Utoprimgen_pick(showmessages, allowlocalfailurefixandnoreport, which,eomtype,parameter,Ugeomfree, ptrgeom, lpflag, pr,pressure,newtonstats,lpflagrad),"step_ch.c:Utoprimgen_tryagain_substep()", "Utoprimgen_pick", 1);
  }

  return(0);
}




// there may be something wrong with this function -- didn't work in TIMEORDER==4, had to do standard method
// could have just been that I wasn't bounding after using this
int Utoprimloop(FTYPE (*U)[NSTORE2][NSTORE3][NPR],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], struct of_newtonstats *newtonstats)
{
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int i,j,k;
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1;

  ZLOOP{
    get_geometry(i, j, k, CENT, ptrgeom);
    // invert True U->p
    int eomtype=EOMDEFAULT;
    MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport,0,&eomtype,EVOLVEUTOPRIM,UEVOLVE,MAC(U,i,j,k), ptrgeom, MAC(prim,i,j,k),newtonstats),"step_ch.c:advance()", "Utoprimgen", 1);
  }
  return(0);
}



// loop for P->U
int primtoUloop(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q;
  int i,j,k;

  ZLOOP{
    MYFUN(get_state(MAC(prim,i,j,k), ptrgeom, &q),"step_ch.c:primtoUloop()", "get_state()", 1);
    get_geometry(i, j, k, CENT, ptrgeom);
    // forward calculate U(p)
    MYFUN(primtoU(UEVOLVE,MAC(prim,i,j,k), &q, ptrgeom, MAC(U,i,j,k)),"step_ch.c:primtoUloop()", "primtoU()", 1);
  }
  return(0);
}



// filter out velocity along field line
void filterffde(int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  //  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  struct of_newtonstats newtonstats;
  int finalstep=0;
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1;

  // now filter velocity
  get_geometry(i,j,k,CENT,ptrgeom);
  get_state(pr,ptrgeom,&q);
  primtoU(UNOTHING,pr,&q,ptrgeom,U);

  //  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion
  int eomtype=EOMDEFAULT;
  Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtype,EVOLVEUTOPRIM,UNOTHING,U,ptrgeom,prout,&newtonstats);

  PALLLOOP(pl) pr[pl]=prout[pl];
  // kill densities
  pr[RHO]=pr[UU]=0.0;



}
