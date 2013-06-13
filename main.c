// #include "decs.h"
#include "defs.h"


// OPENMPMARK: static ok here
int main(int argc, char *argv[])
{
  //  extern void testffdeinversion(void);
  static SFTYPE comptstart;
  static SFTYPE diagtstart;

  
  // uncomment if get SIGFPE 8 error at run with floats and no result problem
  // signal(8,SIG_IGN);


  crapdebug=0; // CHANGINGMARK: to debug stuff


  // start timing for init -> computation start
  timecheck(INITSTARTTIME,0);



  //////////////////////////
  //
  // perform initializations
  //
  //////////////////////////
  if (init(&argc,&argv) >= 1) {
    dualfprintf(fail_file, "main:init: failure\n");
    myexit(1);
  }


  /////////////////////
  //
  // write version header info for perf and step logs
  //
  /////////////////////
  if((myid<=0)&&(appendold==0)){
    if(DOLOGSTEP){
      myfprintf(logstep_file,"#%10s\n%10d %10d\n#","STEPVER",STEPVER,STEPTYPE);
    }
    if(DOLOGPERF){
      myfprintf(logperf_file,"#%10s\n%10d %10d\n#","PERFVER",PERFVER,PERFTYPE);
    }
  }





  // end timing for init -> computation start
  timecheck(INITSTOPTIME,0);

  /////////////////////////
  //
  // set initial computational time [start before diags since report how diags affect computational time]
  //
  /////////////////////////

  // GODMARK: Should store timestart and initial "t" time so running-average can properly be computed when restart simulation.
  //  comptstart=t; // not right
  comptstart=0; // assumes normal simulation started t=0 at nstep=0
  // start timing
  timecheck(STARTTIME,comptstart);
          

  //////////////////////
  //
  // Do initial diagnostics
  //
  //////////////////////
  if (DODIAGS) {
    diagtstart=0; timecheck(DIAGSTARTTIME,diagtstart);

    trifprintf("proc: %04d : Start initial diagnostics\n", myid);
    // no error_check since if init passed, diag(0) should pass
    diag(INIT_OUT,t,nstep,realnstep);
    trifprintf("proc: %04d : End initial diagnostics\n", myid);

    timecheck(DIAGSTOPTIME,diagtstart);
  }

  //  testffdeinversion();
  //  myexit(0);




  //////////////////////
  //
  // START MAJOR WHILE LOOP OVER ITERATIONS
  //
  //////////////////////

  trifprintf("proc: %04d : Start computation\n", myid);

  onemorestep=reallaststep=0; // global variable set also in step_ch.c
  while (reallaststep==0) {

    /* step variables forward in time */
    nstroke = 0;

    // take full step with all non-physics evolutionary code as well
    // SUPERGODMARK: Should insert *all* global variables here, including non-array one's, so that call is modular
    // Then can vary ACTIVEREGION or ACTIVESUBREGION or something on this scale to do AMR and ATR
    // currently ok to have emf as Atemp since both are used as temp variables to get flux or vpot independently
    // ulastglobal used for evolving vpot, which occurs after all steps, so ok to use ulastglobal

    int truestep=1; // indicates true time step and not fake pass-through (as used by, e.g., metric gravity update)

    step_ch_full(truestep,GLOBALPOINT(pglobal),GLOBALPOINT(pstagglobal),GLOBALPOINT(unewglobal),GLOBALPOINT(vpotarrayglobal),GLOBALPOINT(Bhatglobal),GLOBALPOINT(gp_l),GLOBALPOINT(gp_r),GLOBALPOINT(F1),GLOBALPOINT(F2),GLOBALPOINT(F3),GLOBALPOINT(emf),GLOBALPOINT(ulastglobal));

    // get total number of inversion steps
    mpiisum0(&nstroke,0);

    
    /* restart dump */
    // if(nstep == 130) restart_write(1) ;
    //    if(nstep==50) break;

    nstep++;
    // restartsteps[whichrestart]=realnstep;


    /* perform diagnostics */  
    // no error check since assume if step_ch passed, diag(1) will pass


    // start diagnostic timing
    if (DODIAGS && !DODIAGEVERYSUBSTEP){
      // start diagnostic timing
      diagtstart=0; timecheck(DIAGSTARTTIME,diagtstart);

      GLOBALPOINT(pdump) = GLOBALPOINT(pglobal);
      diag(DUMP_OUT,t,nstep,realnstep);
#if(PRODUCTION==0)
      trifprintf( "D");
#endif
      // stop diagnostic timing
      timecheck(DIAGSTOPTIME,diagtstart);
    }


    // time check (should come after first nstep++ to setup things first)
    timecheck(CHECKTIME,comptstart);
    
    // speed check
    timecheck(SPEEDTIME,comptstart);
    
    
    // output timestep info
    output_steptimedt_info(comptstart);
    

  }
  trifprintf("proc: %04d : End computation\n", myid);

  //////////////////////
  //
  // END MAJOR WHILE LOOP OVER ITERATIONS
  //
  //////////////////////


  // stop timing
  timecheck(STOPTIME,comptstart);


  /////////////////////
  //
  // do final diagnostics
  //
  /////////////////////
  if (DODIAGS){
    diagtstart=0; timecheck(DIAGSTARTTIME,diagtstart);

    GLOBALPOINT(pdump) = GLOBALPOINT(pglobal);
    diag(FINAL_OUT,t,nstep,realnstep);

    timecheck(DIAGSTOPTIME,diagtstart);
  }

  ////////////
  //
  // report if go initiated stopping of simulation
  //
  ////////////
  gocheck(STOPTIME);

  // report final timing and other simulation performance characteristics
  timecheck(REPORTTIME,comptstart);


  myexit(0);
  return (0);
}






int gocheck(int whichlocation)
{
  static int goend=0; // init to 0
  // non-static
  char stemp[MAXFILENAME];
  char goch;
  FILE *gogo_file,*gocont_file;



  if(whichlocation==STARTTIME){

    if(myid<=0){
      if(CHECKCONT){
      
        sprintf(stemp,"%sgo.cont",DATADIR);
      
        if((gocont_file=fopen(stemp,"rt"))==NULL){
          dualfprintf(fail_file,"WARNING: Could not open go.cont file: %s , assume user doesn't want to use it\n",stemp);
          //    myexit(1); // can't exit yet if want clean MPI exit
          goch='z';
        }
        else goch='a';
      }
    }

#if(USEMPI)
    MPI_Bcast(&goch,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
#endif

    if(goch=='z'){
      // myexit(1);
      // for now just assume if file doesn't exist that user didn't want to restart
    }
    else{
      if(CHECKCONT){
        if(myid<=0){
          goch=fgetc(gocont_file);
          if( (goch=='y')||(goch=='Y')){
            gocont=1;
            trifprintf("#go.cont called\n");

            fscanf(gocont_file,"%d",&runtype); // can be used to specify which restart file to use among other things
          }
          fclose(gocont_file);
        }
      
    
        if(numprocs>1){
#if(USEMPI)
          MPI_Bcast(&gocont,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
          MPI_Bcast(&runtype,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
#endif
        }
      }// end if checking cont file

    }// end if cont file exits
  }
  else if(whichlocation==CHECKTIME){
    // check if user wants to stop or not(go.go)
    if(myid==0){

      if(!(nstep%NGOCHECK)){
        sprintf(stemp,"%sgo.go",DATADIR);
        
        if((gogo_file=fopen(stemp,"rt"))==NULL){
          //    dualfprintf(fail_file,"Could not open go file: %s\n",stemp);
          //    myexit(1);
          // just assume user didn't want to use this file
        }
        else{
          goch=fgetc(gogo_file);
          if( (goch=='n')||(goch=='N')){
            goend=1;
            trifprintf("#go.go called\n");
          }
          fclose(gogo_file);      
          // if myid!=0 and numprocs>1 then could deal with this, but messy due to timers
        
        }// end if go file exists
      }// end if time to check go file

    }// end myid==0

#if(USEMPI)
    MPI_Bcast(&goend,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);
#endif
    trifprintf("#proc: %s go.go called\n",myidtxt);
    if(goend==1) reallaststep=1;

  }
  else if(whichlocation==STOPTIME){
    if(goend==1) trifprintf("proc: %s Go end called(go.go)\n",myidtxt) ;
  }




  return(0);

}












