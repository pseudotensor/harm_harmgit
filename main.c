/*! \mainpage HARM Documentation Page
 *
 *


 * \section intro_sec Introduction
 *
 * HARM (and HARMARD) solve the GRMHD (radiative) equations of motion.  This page itself gives some documentation/tutorials about how to use harm.

 * The "Files" Link and "search" functions and "Data Structures" Link allow one to see the code layout, direct code documentation, and all other doxygen related context like the list of all global variables.

 * \section code_sec Code
 * Code is present at: <a href="https://harm.unfuddle.com">HARM Unfuddle</a>

  <a href="svngit_8txt.html">SVN GIT Notes</a>


 *
 * \section install_sec Installation and Quick Start Guides

 * Note, if already running on Ubuntu with pre-installed packages by Jon or others like on many supercomputers, then can skip all apt-get or similar package install commands.

<a href="quick__start__guide_21_8txt.html">note: git code and compile and run</a>


<a href="_o_s_xinstallation_8txt.html">note: OSX installation issues#1</a>

<a href="_o_s_xinstallation__changestocode_8txt.html">note: OSX installation issues#2</a>

See also harmgit/makefiles.other

<a href="quick__start__guide_25_8txt.html">note: Fieldline files to Viz5D file</a>

<a href="quick__start__guide_22_8txt.html">note: install r8 stuff to view images from harm</a>

See harmgit/r8toras directory.

<a href="quick__start__guide_26_8txt.html">note: How to setup new problem</a>

<a href="quick__start__guide_23_8txt.html">note: About diagnostics outputted by HARM</a>

<a href="quick__start__guide_27_8txt.html">note: How to use SM</a>

<a href="quick__start__guide_24_8txt.html">note: Compile and run Viz5D</a>

  See harmgit/docs/ for other docs that aren't in txt format (pdfs, latex, png, etc.)


 * \section data_sec HARM data description

  <a href="datanewdesc_8txt.html">new data description</a>

  <a href="tousedata_8txt.html">To use data, follow this</a>

  <a href="datadesc_8txt.html">data description</a>

  <a href="datadesc__mb09_8txt.html">MB09 data description</a>


 * \section viz_sec Analysis and Viz Stuff

<a href="general__plotting__guide_8txt.html">note: Jon's general plotting guide using Python with full detailed Tutorial</a>


<a href="guide_8txt.html">note: Sasha Python help</a>

<a href="viz_21_8txt.html">note: Viz routines/scripts for vis5d</a>


 * \section debug_sec Debugging

 <a href="debug_8txt.html">Debug</a>


 * \section emacs_sec Efficient use of emacs

 <a href="emacsefficient_8txt.html">Emacs</a>



 * \section filetransfer_sec File transfer

 Use globusconnect as part of globusonline:  <a href="http://globusonline.org/">GlobusOnline</a>

 Some example globusonline commands: <a href="globusonline_8txt.html">Example GlobusOnline commands</a>

 I used to use bbcp: <a href="bbcp_8txt.html">bbcp</a> and  <a href="tocopy_8txt.html">bbcp more notes</a>

 Nothing is reliable and nothing is as reliable as globusconnect.

 * \section callgraph_sec Optimizations and Call Graphs

<a href="optimizations_8txt.html">Optimization Notes</a>

<a href="callgraphs_8txt.html">Callgraphs</a>

<a href="installperfstuff_8txt.html">Performance profile software</a>

See harmgit/performancedata directory.

 * \section parallel_sec MPI and OpenMP Notes and SuperComputers

<a href="parallel_8txt.html">Parallel Notes</a>

<a href="supercomputertips.html">SuperCompute Tips</a>

  See harmgit/environmentfiles/ for environment files for other computers, including supercomputers.  Note that there are "hidden" . files in these subdirectories.


  See harmgit/batches for batch system files for various supercomputers running harm.

<a href="batch1_8txt.html">note: Batch queue dependency lists</a>


<a href="_xvfb_8txt.html">Using X remotely</a>

 * \section eos_sec Equation of State Notes

<a href="eos_21_8txt.html">note: Install and Compile EOS stuff and generate stellar model</a>

<a href="eos_22_8txt.html">note: General EOS table from EOS Fortran code</a>

<a href="eos_23_8txt.html">note: Running harm with EOS</a>

<a href="eos_24_8txt.html">note: About Ynu variable</a>

<a href="eos_8c.html">note: About </a>

See also harmgit/eosstuff


 * \section othercode_sec Other code for harm

See harmgit/initboundcode and harmgit/initbounddata for code/data for other initial conditions/boundary conditions.

See harmgit/subcodesfromothers for codes by other people.


 * \section pnmhd_sec Docs on related PNMHD code

<a href="pnmhd_21_8txt.html">note: Avery wind</a>

<a href="pnmhd_24_8txt.html">note: Initial and boundary conditions</a>

<a href="pnmhd_27_8txt.html">note: Details on Code</a>

<a href="pnmhd_22_8txt.html">note: 2D simulations</a>

<a href="pnmhd_25_8txt.html">note: Mathematica files</a>

<a href="pnmhd_23_8txt.html">note: Gravity Potential</a>

<a href="pnmhd_26_8txt.html">note: Modifying Viscosity</a>

* \section scripts_sec Scripts locations

See harmgit/scripts for many scripts that do many things



* \section utils_sec Utilities locations

See harmgit/utils for many utilities that do many things.

See also harmgit/homescripts for other scripts.


* \section license_sec License info

<a href="license_8txt.html">License</a>

<a href="harm__authorship__policy_8txt.html">Harm Authorship Policy</a>




 */




/*! \file main.c
    \brief Main file
    
    Includes main() and gocheck()
*/


////////////////////////////////////////////////////////////////////////
//
////// makefile descriptions
//
//////////////////////////////////////

/*! \file makehead.inc
    \brief Header for makefile.  Sets primary code compile conditionals.
*/

/*! \file maketail.harm.inc
    \brief Tail for makefile for harm related files
*/

/*! \file maketail.inc
    \brief Tail for makefile with various final links
*/

/*! \file maketail.ldouble.deps.inc
    \brief Tail including dependencies for long double makefile.
    See scripts/longdouble2double.sh
*/

/*! \file maketail.ldouble.inc
    \brief Tail for long double makefile.
    See scripts/longdouble2double.sh
*/

/*! \file maketailsuperlong.inc
    \brief Tail for super long double makefile.
    See scripts/longdouble2double.sh
*/

////////////////////////////////////////////////////////////////////////




// #include "decs.h"
#include "defs.h"

///
///
///    Main function that primarily calls init(), diag(), and step_ch_full() in loop.
///
int main(int argc, char *argv[])
{
  //  extern void testffdeinversion(void);
  //OPENMPMARK: static ok here
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

  int fakenstep=0;
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

    if(nstep%10 && PRODUCTION==0){
      dualfprintf(fail_file,"nstroke: %ld\n",nstroke);
    }

    // get total number of inversion steps
    mpiisum0(&nstroke,0);

    
    /* restart dump */
    // if(nstep == 130) restart_write(1) ;
    //if(fakenstep==50) break;


    nstep++;
    fakenstep++;
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





/*! \fn gocheck
    \brief gocheck sees if should stop code

    go.cont or go.go exists and controls stopping behavior of code while it runs.

*/
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
          // myexit(1); // can't exit yet if want clean MPI exit
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
          // dualfprintf(fail_file,"Could not open go file: %s\n",stemp);
          // myexit(1);
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












