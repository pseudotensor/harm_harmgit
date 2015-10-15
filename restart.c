
/*! \file restart.c
     \brief Functions related to restarting code and reading/writing rdump's
     restart functions; restart_init and restart_dump
*/


#include "decs.h"

#include "restart.h"


static int restart_process_extra_variables(void);

static int restartupperpole_set(void);


int extrarestartfunction_new(void)
{
  // nothing new to do
  return(0);
}


/// This function reads header and data (primitives and conserved quantities)
/// GODMARK: fills in *global* quantities
/// It then outputs restart file for user to check consistency
/// CHANGINGMARK: Note that for evolving metric should really store old metric in restart file
int restart_init(int which)
{



  trifprintf("begin restart init\n");

  ////////////////
  //
  // set random seed
  //
  ////////////////
  ranc(1,0);

  ////////////////
  //
  // get restart header and grid data
  //
  ////////////////
  restart_read(which);


  ////////////////
  //
  // set any extra things (only used with restart.rebeccaoldcode.c)
  //
  ////////////////
  trifprintf("before extrarestartfunction()\n");
  extrarestartfunction();


  ////////////////
  //
  // get restart old metric
  //
  ////////////////
  if(DOEVOLVEMETRIC){
    trifprintf("before restartmetric_read(which)\n");
    restartmetric_read(which);
  }

  ////////////////
  //
  // process variables not written to file
  //
  ////////////////
  trifprintf("before restart_process_extra_variables()\n");
  restart_process_extra_variables();
    

  ////////////////
  //
  // set all CPUs to have same restart header stuff (needed in mpicombine==1 mode)
  //
  ////////////////
  trifprintf("before restart_read_defs()\n");
  restart_read_defs();

  ////////////////
  //
  // set coordinate parameters by reading coordparms.dat
  // assumes doesn't overwrite what just set (i.e. what's written to restart file)
  //
  ////////////////
  trifprintf("before read_coord_parms()\n");
  read_coord_parms(defcoord);

  ////////////////
  //
  // write header to screen to all CPUs
  //
  ////////////////
  trifprintf("before write_restart_header(TEXTOUTPUT,log_file)\n");
  logfprintf("header contents below\n"); 
  if(log_file) write_restart_header(RESTARTDUMPTYPE,dnumversion[RESTARTDUMPTYPE],dnumcolumns[RESTARTDUMPTYPE],TEXTOUTPUT,log_file);

  ////////////////
  //
  // check if failed inside rdump file
  //
  ////////////////
  if(failed!=0){
    dualfprintf(log_file,"WARNING: failed=%d in rdump file.  Must clense the way between us.\n",failed);
    dualfprintf(log_file,"Setting failed=0\n");
    failed=0;
  }

  // certain global.h parameters may override restart parameters.  This happens, say, when the user decides not to DOAVG after restart, but DOAVGDIAG is still set.  Without this the code will segfault in diag.c
  // corresponds to "if" statements in initbase.c as well, where user has either choice or no choice.
  if(DOAVG==0) DOAVGDIAG=0;



  ////////////////////////
  //
  // test read by looking at written file
  //
  ////////////////////////
  trifprintf("before restart_write(3)\n");
  restart_write(3);

  // can't check (e.g.) images here if wanting to look at conserved quantities

  trifprintf("proc: %d t=%21.15g failed=%d\n", myid, t, failed);


  // very basic checks on primary grid-based quantities read from restart dump
  restart_init_simple_checks(1);


  trifprintf("end restart_init\n");


  return(0);

}







/// assign any MAC(quantity,i,j,k) that isn't set by restart_read() that want to set so later checks or later computational uses can occur without errors
/// GODMARK: must keep this up to date with what one is evolving not written to dump
/// OPENMPOPTMARK: Assume not important to optimize
int restart_process_extra_variables(void)
{
  int i,j,k;

  // if evolving entropy need to reset primitive to avoid error checks detecting problem
  if(DOENTROPY!=DONOENTROPY){
    FULLLOOP GLOBALMACP0A1(pglobal,i,j,k,ENTROPY) = GLOBALMACP0A1(pglobal,i,j,k,UU);
  }


  return(0);

}



int restart_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  long truedump_cnt;


  // get special upperpole restart header and grid data (do inside restart_read() since always want this file with standard restart file)
  if(FLUXB==FLUXCTSTAG && special3dspc==1 && N3>1) restartupperpole_write(dump_cnt);



  trifprintf("begin dumping rdump# %ld ... ",dump_cnt);

  whichdump=RESTARTDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rdump");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restart_header,rdump_content)>=1) return(1);

  trifprintf("end dumping rdump# %ld ... ",dump_cnt);


  return(0);

}



/// number of columns for restart
void set_rdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  // always NPR
  //*numcolumns=NPR*2; // primitives and conservatives
  //*numcolumns=NPR; // primitives only


  // spatial debug counters not crucial
  //  *numcolumns=NPR*2 + dnumcolumns[VPOTDUMPTYPE] + dnumcolumns[FAILFLOORDUDUMPTYPE] + dnumcolumns[DEBUGDUMPTYPE] ;

  *numcolumns=NPR*2 + dnumcolumns[DISSDUMPTYPE] + dnumcolumns[FAILFLOORDUDUMPTYPE] ;

  
  if(EVOLVEWITHVPOT||TRACKVPOT){
    // even with TRACKVPOT, with vpot as diagnostic, don't regenerate vpot from B, so need to store in restart file so can continue updating it.s
    *numcolumns += dnumcolumns[VPOTDUMPTYPE];
  }

  *numversion=1;
}


/// Note that in initbase.c that many required things are computed during restart setup
/// For example, post_init() calls compute_EOS_parms_full() that fills EOSextraglobal[] if doing WHICHEOS==KAZFULL, so don't have to store that in restart file.
int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // ensure what is written (even at t=0) for fake entropy primitive is just internal energy, as assumed through-out evolution.
  if(ENTROPY>=0){
    GLOBALMACP0A1(pglobal,i,j,k,ENTROPY)=GLOBALMACP0A1(pglobal,i,j,k,UU);
  }

  // always NPR
  myset(datatype,GLOBALMAC(pglobal,i,j,k),0,NPR,writebuf);
  myset(datatype,GLOBALMAC(unewglobal,i,j,k),0,NPR,writebuf);

  // NOTEMARK: see also dump.c

  if(EVOLVEWITHVPOT||TRACKVPOT){
    if(dnumcolumns[VPOTDUMPTYPE]>0){
      int jj;
      for(jj=0;jj<dnumcolumns[VPOTDUMPTYPE];jj++){
        myset(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,j,k),0,1,writebuf); // 1 each
      }
    }
  }
  
  if(dnumcolumns[DISSDUMPTYPE]>0){
    myset(datatype,&GLOBALMAC(dissfunpos,i,j,k),0,dnumcolumns[DISSDUMPTYPE],writebuf);
  }


  if(dnumcolumns[FAILFLOORDUDUMPTYPE]>0){
    myset(datatype,GLOBALMAC(failfloordu,i,j,k),0,dnumcolumns[FAILFLOORDUDUMPTYPE],writebuf);
  }

  // too many of these and not crucial since just counters
  //  if(dnumcolumns[DEBUGDUMPTYPE]>0){
  //    myset(datatype,GLOBALMAC(failfloorcount,i,j,k),0,dnumcolumns[DEBUGDUMPTYPE],writebuf);
  //  }


  return(0);
}








/// read restart file
int restart_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  int bintxt;


  // get special upperpole restart header and grid data (do inside restart_read() since always want this file with standard restart file)
  if(FLUXB==FLUXCTSTAG && special3dspc==1 && N3>1) restartupperpole_read(dump_cnt);
  else restartupperpole_set();


  trifprintf("begin reading rdump# %ld ... ",dump_cnt);

  whichdump=RESTARTDUMPTYPE;
  datatype=MPI_FTYPE;
  bintxt=binaryoutput;
  strcpy(fileprefix,"dumps/rdump");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
 
  int failreturn;
  failreturn=dump_gen(READFILE,dump_cnt,bintxt,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restart_header,rdump_read_content);
  
  if(failreturn==FILENOTFOUND){
    dualfprintf(fail_file,"restart file not found\n");
    myexit(87343363);
  }
  else if(failreturn>0){
    dualfprintf(fail_file,"restart file: other read error\n");
    myexit(7165766);
  }



  // NOTE: for FLUXB==FLUXCTSTAG, unewglobal and vpot are bounded in initbase.c for boundary edges and inter-MPI edges

  trifprintf("end reading rdump# %ld ... ",dump_cnt);


  return(0);

}


/// read-restart content
int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // always NPR
  myget(datatype,GLOBALMAC(pglobal,i,j,k),0,NPR,writebuf);
  myget(datatype,GLOBALMAC(unewglobal,i,j,k),0,NPR,writebuf);

  // NOTEMARK: see also dump.c

  if(EVOLVEWITHVPOT||TRACKVPOT){
    if(dnumcolumns[VPOTDUMPTYPE]>0){
      int jj;
      for(jj=0;jj<dnumcolumns[VPOTDUMPTYPE];jj++){
        myget(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,j,k),0,1,writebuf); // 1 each
      }
    }
  }

  if(dnumcolumns[DISSDUMPTYPE]>0){
    myget(datatype,&GLOBALMAC(dissfunpos,i,j,k),0,dnumcolumns[DISSDUMPTYPE],writebuf);
  }


  if(dnumcolumns[FAILFLOORDUDUMPTYPE]>0){
    myget(datatype,GLOBALMAC(failfloordu,i,j,k),0,dnumcolumns[FAILFLOORDUDUMPTYPE],writebuf);
  }

  // counters not crucial
  //  if(dnumcolumns[DEBUGDUMPTYPE]>0){
  //    myget(datatype,GLOBALMAC(failfloorcount,i,j,k),0,dnumcolumns[DEBUGDUMPTYPE],writebuf);
  //  }


  

  return(0);
}











/// read-restart file with upper pole information for FLUXB==FLUXCTSTAG and special3dspc==1
/// Since both write *and* read full 3D file with j shifted up, read-in puts quantities in correct location, so no remapping needed.
/// All j=N2 values will be overwritten by normal MPI calls, except the last one for mycpupos[2]==ncpux2-1 that will use the read-in values
/// Note that order of normal restart read/write and upperpole read/write doesn't matter since all values are correctly positioned during read-write.
/// This is all done because normal MPI fileio routines assume standard j=0..N-1 or j=1..N block, but staggered stuff has j=0..N block for spherical polar coordinates.  No other directions (i.e. i or k) require this.
/// For simplicity, we just output two full 3D files with j shifted.
int restartupperpole_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  int bintxt;
  int restartupperpole_set(void);


  trifprintf("begin reading rdumpupperpole# %ld ... ",dump_cnt);

  whichdump=RESTARTUPPERPOLEDUMPTYPE;
  datatype=MPI_FTYPE;
  bintxt=binaryoutput;
  strcpy(fileprefix,"dumps/rdumpupperpole");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
 
  int failreturn;
  failreturn=dump_gen(READFILE,dump_cnt,bintxt,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restartupperpole_header,rupperpoledump_read_content);

  if(failreturn==FILENOTFOUND){
    dualfprintf(fail_file,"SUPERWARNING: resetting upperpole to zero since no restart file for upperpole is available\n");
    restartupperpole_set();
  }
  else if(failreturn>0) return(1);
    


#if(0)
  // NOT USED ANYMORE, but, if did remapping, it would look something like the below.
  {
    // SPECIAL CASE:
    // now that have read-in shifted values, must assign them to correct locations
    
    int i,j,k,pliter,pl;
    jshifted=j+SHIFT2;
    LOOPF1 LOOPF3{
      PLOOP(pliter,pl){
        // if(pl==B1 || pl==B2 || pl==B3){
        if(pl==B2){
          GLOBALMACP0A1(unewglobal,i,jshifted-SHIFT2,k,pl)=GLOBALMACP0A1(unewglobal,i,jshifted,k,pl);
        }
      }
    }
  }
#endif


  trifprintf("end reading rdumpupperpole# %ld ... ",dump_cnt);


  return(0);

}




/// restart needs to set B2=0 along outer pole if not being used.
static int restartupperpole_set(void)
{

  if(mycpupos[2]==ncpux2-1){// only need to operate on true upper pole
    // then assume user knows what they are doing and just set array to zero
    int i,j,k,jshifted;
    DUMPGENLOOP{ // same loop as dump_gen() uses
      if(j!=N2-1) continue; // force only assignment right at j==N2 so still doesn't matter what order the normal and upperpole restart calls are made
      jshifted=j+SHIFT2;
      GLOBALMACP0A1(unewglobal,i,jshifted,k,B2)=0.0;

      if(EVOLVEWITHVPOT||TRACKVPOT){
        if(dnumcolumns[VPOTDUMPTYPE]>0){
          int jj;
          for(jj=0;jj<dnumcolumns[VPOTDUMPTYPE];jj++){
            if(jj==2) continue; // skip A_2 that is not on pole, so not needed
            GLOBALMACP1A0(vpotarraydump,jj,i,jshifted,k)=0.0;
          }
        }
      }
    }// end FULLLOOP
  }// end if true upper pole

  return(0);
}


/// read-restart upperpole header
int read_restartupperpole_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{

  if(bintxt==BINARYOUTPUT){
  }
  else{
    // nothing so far, but must have new line if not NULL function
    fprintf(headerptr,"\n");
    fflush(headerptr);
  }

  return(0);
}

/// write-restart upperpole header
int write_restartupperpole_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{

  if(bintxt==BINARYOUTPUT){
  }
  else{
    // nothing so far, but must have new line if not NULL function
    fprintf(headerptr,"\n");
    fflush(headerptr);
  }

  return(0);
}


/// read-restart upperpole content
int rupperpoledump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  int jshifted;

  // j+SHIFT2 so that fake 3D read-in fills j=N2 for mycpupos[2]==ncpux2-1 (tj=totalsize[2]) when MPI routines think they are accessing j=N2-SHIFT2
  jshifted=j+SHIFT2;

  // Just B2 is on pole
  myget(datatype,GLOBALMAC(unewglobal,i,jshifted,k),B2,1,writebuf);

  // NOTEMARK: see also dump.c
  if(EVOLVEWITHVPOT||TRACKVPOT){
    if(dnumcolumns[VPOTDUMPTYPE]>0){
      int jj;
      for(jj=0;jj<dnumcolumns[VPOTDUMPTYPE];jj++){
        if(jj==2) continue; // skip A_2 that is not on pole, so not needed
        myget(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,jshifted,k),0,1,writebuf); // 1 each
      }
    }
  }

 

  return(0);
}





/// write-restart file with upper pole information for FLUXB==FLUXCTSTAG and special3dspc==1
int restartupperpole_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  long truedump_cnt;

  trifprintf("begin dumping rdumpupperpole# %ld ... ",dump_cnt);

  whichdump=RESTARTUPPERPOLEDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rdumpupperpole");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restartupperpole_header,rupperpoledump_content)>=1) return(1);

  trifprintf("end dumping rdumpupperpole# %ld ... ",dump_cnt);


  return(0);

}




/// read-restart upperpole numcolumns
/// number of columns for restart
void set_rupperpoledump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  // Just B2 is on pole
  *numcolumns=1;

  
  if(EVOLVEWITHVPOT||TRACKVPOT){
    // even with TRACKVPOT, with vpot as diagnostic, don't regenerate vpot from B, so need to store in restart file so can continue updating it.s
    *numcolumns += NDIM-1; // only A_0, A_1, and A_3 are on pole
  }

  *numversion=0;
}

/// write-restart upperpole numcolumns
void set_rupperpoledump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{
  extern void set_rupperpoledump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);

  // same as read:
  set_rupperpoledump_read_content_dnumcolumns_dnumversion(numcolumns,numversion);

}



/// write-restart upperpole content
int rupperpoledump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  int jshifted;

  // j+SHIFT2 so that fake 3D read-in fills j=N2 for mycpupos[2]==ncpux2-1 (tj=totalsize[2]) when MPI routines think they are accessing j=N2-SHIFT2
  jshifted=j+SHIFT2;

  // Only B2 is on pole
  myset(datatype,GLOBALMAC(unewglobal,i,jshifted,k),B2,1,writebuf);

  // NOTEMARK: see also dump.c
  if(EVOLVEWITHVPOT||TRACKVPOT){
    if(dnumcolumns[VPOTDUMPTYPE]>0){
      int jj;
      for(jj=0;jj<dnumcolumns[VPOTDUMPTYPE];jj++){
        if(jj==2) continue; // skip A_2 that is not on pole, so not needed
        myset(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,jshifted,k),0,1,writebuf); // 1 each
      }
    }
  }


  return(0);
}



















/// write metric restart file
/// This is done to keep track of older time's metric so can take temporal difference to get connection.  Otherwise, restart would be missing the temporal change in metric contribution to the connection.
int restartmetric_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  long truedump_cnt;

  trifprintf("begin dumping rmetricdump# %ld ... ",dump_cnt);

  whichdump=RESTARTMETRICDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rmetricdump");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");


  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restartmetric_header,rmetricdump_content)>=1) return(1);

  trifprintf("end dumping rmetricdump# %ld ... ",dump_cnt);


  return(0);

}


void set_rmetricdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{
  extern void set_rmetricdump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);

  // same as read:
  set_rmetricdump_read_content_dnumcolumns_dnumversion(numcolumns,numversion);

}


/// must be consistent with dnumcolumns[RESTARTMETRICDUMPTYPE]
int rmetricdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int gridpos;
  struct of_compgeom *tempcompgeom;
  FTYPE generalmatrixlower[NDIM][NDIM];
  FTYPE generalmatrixupper[NDIM][NDIM];
  int jj,kk;


  GRIDLOOP(gridpos){ // note that in older code gridpos was quickest index!  Now slowest compared to matrix indices
    
#if(NEWMETRICSTORAGE)

    tempcompgeom = &GLOBALMETMACP1A0(compgeomlast,gridpos,i,j,k);


    DLOOP(jj,kk){
      generalmatrixlower[jj][kk]=tempcompgeom->gcov[GIND(jj,kk)];
      generalmatrixupper[jj][kk]=tempcompgeom->gcon[GIND(jj,kk)];
    }

    myset(datatype,generalmatrixlower,0,NDIM*NDIM,writebuf);
    myset(datatype,generalmatrixupper,0,NDIM*NDIM,writebuf);
    myset(datatype,tempcompgeom->gcovpert,0,NDIM,writebuf);
    myset(datatype,&(tempcompgeom->gdet),0,1,writebuf);
#if(WHICHEOM!=WITHGDET)
    myset(datatype,&(tempcompgeom->EOMFUNCMAC(0)),0,NPR,writebuf);
    myset(datatype,&(tempcompgeom->IEOMFUNCNOSINGMAC(0)),0,NPR,writebuf);
#endif
#if(GDETVOLDIFF)
    myset(datatype,&(tempcompgeom->gdetvol),0,1,writebuf);
#endif
    myset(datatype,&(tempcompgeom->igdetnosing),0,1,writebuf);
    myset(datatype,&(tempcompgeom->alphalapse),0,1,writebuf);
    myset(datatype,&(tempcompgeom->betasqoalphasq),0,1,writebuf);
    myset(datatype,&(tempcompgeom->beta),0,NDIM,writebuf);

#else

    DLOOP(jj,kk){
      generalmatrixlower[jj][kk]=GLOBALMETMACP1A1(gcovlast,gridpos,i,j,k,GIND(jj,kk));
    }
      
    myset(datatype,generalmatrixlower,0,NDIM*NDIM,writebuf);
    myset(datatype,GLOBALMETMACP1A0(gcovpertlast,gridpos,i,j,k),0,NDIM,writebuf);
    myset(datatype,&(GLOBALMETMACP1A0(alphalapselast,gridpos,i,j,k)),0,1,writebuf);

#endif

  }// end over gridpos  


  return(0);
}







int write_restartmetric_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{

  if(bintxt==BINARYOUTPUT){
  }
  else{
    // nothing so far, but must have new line if not NULL function
    fprintf(headerptr,"\n");
    fflush(headerptr);
  }

  return(0);
}

int read_restartmetric_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{

  if(bintxt==BINARYOUTPUT){
  }
  else{
    // nothing so far, but must have new line if not NULL function
    fprintf(headerptr,"\n");
    fflush(headerptr);
  }

  return(0);
}


int restartmetric_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
  int bintxt;


  trifprintf("begin reading rmetricdump# %ld ... ",dump_cnt);

  dualfprintf(fail_file,"Reading in old metric, but so far need to bound somehow since no BCs but needed!\n");
  myexit(2496736);

  whichdump=RESTARTMETRICDUMPTYPE;
  datatype=MPI_FTYPE;
  bintxt=binaryoutput;
  strcpy(fileprefix,"dumps/rmetricdump");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
 
  int failreturn;
  failreturn=dump_gen(READFILE,dump_cnt,bintxt,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restartmetric_header,rmetricdump_read_content);
  
  if(failreturn==FILENOTFOUND){
    dualfprintf(fail_file,"restartmetric file not found\n");
    myexit(26525667);
  }
  else if(failreturn>0){
    dualfprintf(fail_file,"restartmetric file: other read error\n");
    myexit(71857346);
  }





  trifprintf("end reading rmetricdump# %ld ... ",dump_cnt);


  return(0);

}



/// number of columns
void set_rmetricdump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{


  if(DOEVOLVEMETRIC){
#if(NEWMETRICSTORAGE)
    *numcolumns=NPG*(NDIM*NDIM + NDIM*NDIM + NDIM + 1 + 1 + 1 + 1 + NDIM); // See restart.c's rmetricdump_content()
#if(WHICHEOM!=WITHGDET)
    *numcolumns+=NPG*(2*NPR);
#endif
#if(GDETVOLDIFF)
    *numcolumns+=NPG*(1);
#endif

#else
    *numcolumns=NPG*(NDIM*NDIM+NDIM+1); // See restart.c's rmetricdump_content()
#endif
  }
  else{
    *numcolumns=0;
  }

  *numversion=0;


}

/// almost the same as the write function except here we use myget() instead of myset()
int rmetricdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int gridpos;
  struct of_compgeom *tempcompgeom;
  FTYPE generalmatrixlower[NDIM][NDIM];
  FTYPE generalmatrixupper[NDIM][NDIM];
  int jj,kk;


  GRIDLOOP(gridpos){ // note that in older code gridpos was quickest index!  Now slowest compared to matrix indices
    
#if(NEWMETRICSTORAGE)

    tempcompgeom = &GLOBALMETMACP1A0(compgeomlast,gridpos,i,j,k);

    myget(datatype,generalmatrixlower,0,NDIM*NDIM,writebuf);
    myget(datatype,generalmatrixupper,0,NDIM*NDIM,writebuf);
    myget(datatype,tempcompgeom->gcovpert,0,NDIM,writebuf);
    myget(datatype,&(tempcompgeom->gdet),0,1,writebuf);
#if(WHICHEOM!=WITHGDET)
    myget(datatype,&(tempcompgeom->EOMFUNCMAC(0)),0,NPR,writebuf);
    myget(datatype,&(tempcompgeom->IEOMFUNCNOSINGMAC(0)),0,NPR,writebuf);
#endif
#if(GDETVOLDIFF)
    myget(datatype,&(tempcompgeom->gdetvol),0,1,writebuf);
#endif
    myget(datatype,&(tempcompgeom->igdetnosing),0,1,writebuf);
    myget(datatype,&(tempcompgeom->alphalapse),0,1,writebuf);
    myget(datatype,&(tempcompgeom->betasqoalphasq),0,1,writebuf);
    myget(datatype,&(tempcompgeom->beta),0,NDIM,writebuf);

    // assign
    DLOOP(jj,kk){
      tempcompgeom->gcov[GIND(jj,kk)]=generalmatrixlower[jj][kk];
      tempcompgeom->gcon[GIND(jj,kk)]=generalmatrixupper[jj][kk];
    }

#else

    myget(datatype,generalmatrixlower,0,NDIM*NDIM,writebuf);
    myget(datatype,GLOBALMETMACP1A0(gcovpertlast,gridpos,i,j,k),0,NDIM,writebuf);
    myget(datatype,&(GLOBALMETMACP1A0(alphalapselast,gridpos,i,j,k)),0,1,writebuf);

    // assign
    DLOOP(jj,kk){
      GLOBALMETMACP1A1(gcovlast,gridpos,i,j,k,GIND(jj,kk))=generalmatrixlower[jj][kk];
    }
#endif      

  }// end over gridpos  


  return(0);



}














/// headerptr created and only used here OR passed a given pointer
int bcast_restart_header(void)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=1;
  
  // 0 and NULL not used
  readwrite_restart_header(NOTHINGHEAD, 0, bcasthead, NULL);

  return(0);
}


/// headerptr created and only used here OR passed a given pointer
int read_restart_header_new(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=0;
  
  readwrite_restart_header(READHEAD, bintxt, bcasthead, headerptr);

  return(0);
}



int write_restart_header_new(int whichdump, int whichdumpversion, int numcolumns, int bintxt,FILE*headerptr)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=0;
  
  readwrite_restart_header(WRITEHEAD, bintxt, bcasthead, headerptr);

  return(0);
}



/// here bintxt is only 
/// GODMARK: For super-consistent constant-format restart read/write, need to have maximum sizes for things like NPR, etc. anything that changes in size when changing paramters.  This way can always read other restart files with different parameters
/// GODMARK: But at same time, if parameters change, restart may not be valid if data doesn't exist, so checks are needed then.
int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr)
{
  int ii;
  int k,dir,pl;
  int enerregion, floor, tscale;
  int idum1,idum2,idum3;
  int dissloop;
  int dtloop;
  char headerone[MAXFILENAME];
  char sheaderone[MAXFILENAME];
  char ctypeheaderone[MAXFILENAME];


  if(readwrite==READHEAD) trifprintf("begin reading header of restart file\n");
  else if(readwrite==WRITEHEAD) trifprintf("begin writing header of restart file\n");


  // setup idum1,2,3 so read and write are processed same way
  if(readwrite==WRITEHEAD){
    idum1=totalsize[1];
    idum2=totalsize[2];
    idum3=totalsize[3];
  }


  if(readwrite==READHEAD){
    sprintf(headerone,"%s",HEADERONEIN);
    sprintf(sheaderone,"%s",SFTYPEHEADERONEIN);
    sprintf(ctypeheaderone,"%s",CTYPEHEADERONEIN);
  }
  else if(readwrite==WRITEHEAD){
    sprintf(headerone,"%s",HEADERONEOUT);
    sprintf(sheaderone,"%s",SFTYPEHEADERONEOUT);
    sprintf(ctypeheaderone,"%s",CTYPEHEADERONEOUT);
  }


  /////////////////
  //
  // START HEADER read/write list
  //
  /////////////////
  int headercount=0;

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&idum1,sizeof(int),"%d",1,MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&idum2,sizeof(int),"%d",1,MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&idum3,sizeof(int),"%d",1,MPI_INT,headerptr); // 3D thing

  // all cpus read the rest of header the same
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&t, sizeof(SFTYPE), headerone, 1, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&tf, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&nstep, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&a, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&MBH, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&QBH, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  if(readwrite!=READHEAD||1){
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&EP3, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&THETAROT, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  }
  else{
    EP3=0.0;
    
    // radians
    THETAROT=1.5708; // as if restarted with this set
  }

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&gam, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&gamideal, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&cour, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);


  // for metric evolution
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Xmetricold, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Xmetricnew, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);

      
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&dt, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&lim[1], sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&lim[2], sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&lim[3], sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&TIMEORDER, sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,fluxmethod, sizeof(int), "%d", NPR, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&FLUXB, sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&UTOPRIMVERSION, sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&failed, sizeof(int), "%d", 1, MPI_INT, headerptr);
  
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&defcoord, sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&R0, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Rin, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Rout, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&hslope, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Zin, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Zout, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Rin_array, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Rout_array, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&BCtype[0],sizeof(int), "%d", COMPDIM*2, MPI_INT, headerptr);

  // new May 6, 2003
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&realnstep, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&debugfail,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&whichrestart,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&cooling,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&restartsteps[0],sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&restartsteps[1],sizeof(long), "%ld", 1, MPI_LONG, headerptr);

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMIEDUMP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMIEIMAGE,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMIEENER,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DODIAGS,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOENERDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOGDUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DORDUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DODUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOAVGDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOIMAGEDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOAREAMAPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DODIAGEVERYSUBSTEP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOENODEBUGEVERYSUBSTEP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  if(readwrite==READHEAD&&0){
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOCOLSPLIT,sizeof(int), "%d", NUMDUMPTYPES-1, MPI_INT, headerptr);
    int l5,inputl5=NUMDUMPTYPES-2;
    for(l5=NUMDUMPTYPES-1;l5>=0;l5--){
      if(l5!=RESTARTUPPERPOLEDUMPTYPE){
        DOCOLSPLIT[l5]=DOCOLSPLIT[inputl5];
        inputl5--;
      }
    }
  }
  else headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOCOLSPLIT,sizeof(int), "%d", NUMDUMPTYPES, MPI_INT, headerptr);


  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&POSDEFMETRIC,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&periodicx1,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&periodicx2,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&periodicx3,sizeof(int), "%d", 1, MPI_INT, headerptr); // 3D thing
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&dofull2pi,sizeof(int), "%d", 1, MPI_INT, headerptr); // 3D thing
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&binaryoutput,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&sortedoutput,sizeof(int), "%d", 1, MPI_INT, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&defcon,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&SAFE,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&RHOMIN,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&UUMIN,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&RHOMINLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&UUMINLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&UULIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&BSQORHOLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&BSQOULIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&UORHOLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMAMAX,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMADAMP,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&GAMMAFAIL,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);      
  // end new May 6, 2003
  
  // reorganized order for DT related stuff
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DTr, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  if(readwrite==READHEAD&&0){
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DTdumpgen[0], sizeof(SFTYPE), sheaderone, NUMDUMPTYPES-1, MPI_SFTYPE, headerptr);
    int l5,inputl5=NUMDUMPTYPES-2;
    for(l5=NUMDUMPTYPES-1;l5>=0;l5--){
      if(l5!=RESTARTUPPERPOLEDUMPTYPE){
        DTdumpgen[l5]=DTdumpgen[inputl5];
        inputl5--;
      }
    }
  }
  else headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DTdumpgen[0], sizeof(SFTYPE), sheaderone, NUMDUMPTYPES, MPI_SFTYPE, headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&rdump_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  if(readwrite==READHEAD&&0){
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&dumpcntgen[0], sizeof(long), "%ld", NUMDUMPTYPES-1, MPI_LONG, headerptr);
    int l5,inputl5=NUMDUMPTYPES-2;
    for(l5=NUMDUMPTYPES-1;l5>=0;l5--){
      if(l5!=RESTARTUPPERPOLEDUMPTYPE){
        dumpcntgen[l5]=dumpcntgen[inputl5];
        inputl5--;
      }
    }
  }
  else headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&dumpcntgen[0], sizeof(long), "%ld", NUMDUMPTYPES, MPI_LONG, headerptr);

  
  
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&prMAX[0],sizeof(FTYPE), headerone, NPR, MPI_FTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&prfloorcoef[0],sizeof(FTYPE), headerone, NPR, MPI_FTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&rescaletype,sizeof(int), "%d", 1, MPI_INT,headerptr);
  
  // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&avgscheme[1],sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&avgscheme[2],sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&avgscheme[3],sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&dofluxreconevolvepointfield,sizeof(int), "%d", 1, MPI_INT,headerptr);

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&do_transverse_flux_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&do_source_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&do_conserved_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&INVERTFROMAVERAGEIFFAILED,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&LIMIT_AC_FRAC_CHANGE,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&MAX_AC_FRAC_CHANGE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOENOFLUX,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&PARAMODWENO,sizeof(int), "%d", 1, MPI_INT,headerptr);


  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&CHECKCONT,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOTSTEPDIAG,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOLOGSTEP,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&DOLOGPERF,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&NDTCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&NZCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&NDTDOTCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&NGOCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&NTIMECHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&PERFWALLTIME,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&ZCPSESTIMATE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);


  // new June 6, 2003 (cumulatives)
  // always NPR
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&pcumreg_tot[0][0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*(COMPDIM*2)*NPR, MPI_SFTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&fladdreg_tot[0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NPR, MPI_SFTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&fladdtermsreg_tot[0][0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NUMFAILFLOORFLAGS*NPR, MPI_SFTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&sourceaddreg_tot[0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NPR, MPI_SFTYPE,headerptr);
  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&sourceaddtermsreg_tot[0][0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NUMSOURCES*NPR, MPI_SFTYPE,headerptr);

  headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&Ureg_init_tot[0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NPR, MPI_SFTYPE,headerptr);
  //FAILFLOORLOOP(indexfinalstep,tscale,floor)
  if(readwrite==READHEAD&&0){
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&failfloorcountlocal_tot[0][0][0],sizeof(CTYPE),ctypeheaderone,2*NUMTSCALES*(NUMFAILFLOORFLAGS-3),MPI_CTYPE,headerptr);
  }
  else{ //new format
    headercount+=header1_gen(!DONOTACCESSMEMORY,readwrite,bintxt,bcasthead,&failfloorcountlocal_tot[0][0][0],sizeof(CTYPE),ctypeheaderone,2*NUMTSCALES*NUMFAILFLOORFLAGS,MPI_CTYPE,headerptr);
  }

  // end new June 6,2003

  headercount+=header1_gen(DODISS,readwrite,bintxt,bcasthead,&dissreg_tot[0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NUMDISSVERSIONS, MPI_SFTYPE,headerptr);
    
  headercount+=header1_gen(DOLUMVSR,readwrite,bintxt,bcasthead,&lumvsr_tot[0],sizeof(SFTYPE),sheaderone, ncpux1*N1, MPI_SFTYPE, headerptr);

  // for consistent restart file, assume NUMDISSVERSIONS doesn't change
  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){// this loop is over pointers, not a continuous memory space!
    headercount+=header1_gen(DODISSVSR,readwrite,bintxt,bcasthead,&dissvsr_tot[dissloop][0],sizeof(SFTYPE),sheaderone, ncpux1*N1, MPI_SFTYPE, headerptr);
  }

 
  // Aug 16, 2007 -- active region parameters (updated by JCM 07/24/08)
  headercount+=header1_gen(DOGRIDSECTIONING,readwrite,bintxt,bcasthead,&global_enerregiondef,sizeof(int), "%d", NUMENERREGIONS*NUMUPDOWN*NDIM, MPI_INT,headerptr);
  headercount+=header1_gen(DOGRIDSECTIONING,readwrite,bintxt,bcasthead,&t_transition_in,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  headercount+=header1_gen(DOGRIDSECTIONING,readwrite,bintxt,bcasthead,&t_transition_out,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  //end Aug 16, 2007 (updated by JCM 07/24/08)
 

  trifprintf("\nheadercount=%d\n",headercount);


  // BELOW moved to dump_gen
  // now read of tail is controlled by dump_gen()
  //  if(bintxt==TEXTOUTPUT){
  // flush to just after the header line in case binary read of data
  //    if(readwrite==READHEAD) while(fgetc(headerptr)!='\n');
  // if(readwrite==READHEAD){
  // do nothing
  // }
  // else if(readwrite==WRITEHEAD) fprintf(headerptr,"\n");
  // }

  
  if(readwrite==WRITEHEAD) {
    if(bintxt==TEXTOUTPUT){
      // flush to just after the header line in case binary read of data
      fprintf(headerptr,"\n");
    }

    // flux header so written to disk fully
    // Only flush on writing as flushing on reading is undefined on VS
    fflush(headerptr);
  }




  if(readwrite==READHEAD){
    /////////////////
    //
    // some checks
    if (idum1 != totalsize[1]) {
      dualfprintf(fail_file, "error reading restart file; N1 differs\n");
      dualfprintf(fail_file, "got totalsize[1]=%d needed totalsize[1]=%d\n",idum1,totalsize[1]);
      myexit(3);
    }
    if (idum2 != totalsize[2]) {
      dualfprintf(fail_file, "error reading restart file; N2 differs\n");
      dualfprintf(fail_file, "got totalsize[2]=%d needed totalsize[2]=%d\n",idum2,totalsize[2]);
      myexit(4);
    }
    if (idum3 != totalsize[3]) {
      dualfprintf(fail_file, "error reading restart file; N3 differs\n");
      dualfprintf(fail_file, "got totalsize[3]=%d needed totalsize[3]=%d\n",idum3,totalsize[3]);
      myexit(4);
    }
  }



  // final things need to set but lock to rdump header content since don't want to add new header entry.  GODMARK: Eventually should add to restart header.
  if(readwrite==READHEAD){
    fakesteps[0]=restartsteps[0];
    fakesteps[1]=restartsteps[1];
    whichfake=whichrestart;
    DTfake=MAX(1,DTr/10); // only thing that matters currently.
  }
  if(bcasthead){
    MPI_Bcast(&fakesteps, 2, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&whichfake, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
    MPI_Bcast(&DTfake, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  }



  if(readwrite==READHEAD) trifprintf("end reading header of restart file\n");
  else if(readwrite==WRITEHEAD) trifprintf("end writing header of restart file\n");












  return(0);



}


/// setup global accounting for CPU=0
int restart_read_defs_new(void)
{
  int enerregion;
  int indexfinalstep,floor,tscale,sc;
  int dissloop;
  int dir,pl,pliter;
  int ii;
  int dtloop;
  int bcast_restart_header(void);
  

  if(myid==0){
    ////////////
    //
    // Define for cpu=0 only, which will continue to keep track of the total after restart
    //
    // Recall that PDUMPLOOP is only for primitives, while PLOOP should always be used for conserved quantities or fluxes
    // pdot and pdotterms computed each time and not cumulative, so no need to store.
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PLOOP(pliter,pl) pcumreg[enerregion][dir][pl]=pcumreg_tot[enerregion][dir][pl];
    ENERREGIONLOOP(enerregion) PLOOP(pliter,pl) fladdreg[enerregion][pl]=fladdreg_tot[enerregion][pl];
    ENERREGIONLOOP(enerregion) PLOOP(pliter,pl) FLOORLOOP(floor) fladdtermsreg[enerregion][floor][pl]=fladdtermsreg_tot[enerregion][floor][pl];
    ENERREGIONLOOP(enerregion) PLOOP(pliter,pl) sourceaddreg[enerregion][pl]=sourceaddreg_tot[enerregion][pl];
    ENERREGIONLOOP(enerregion) PLOOP(pliter,pl) SCLOOP(sc) sourceaddtermsreg[enerregion][sc][pl]=sourceaddtermsreg_tot[enerregion][sc][pl];
    ENERREGIONLOOP(enerregion) PLOOP(pliter,pl) Ureg_init[enerregion][pl]=Ureg_init_tot[enerregion][pl];

    if(DODEBUG){
      FAILFLOORLOOP(indexfinalstep,tscale,floor){
        failfloorcountlocal[indexfinalstep][tscale][floor]=failfloorcountlocal_tot[indexfinalstep][tscale][floor];
        // failfloorcountlocal overwritten by counttotal by integratel in dump_ener.c, so also put in spatial spot.  Need this if not tracking counters spatially.  Or shouldn't reset failfloorcountlocal to zero in counttotal in dump_ener.c and reset counters to zero elsewhere at start of simulation.
        // So just stick it somewhere we can easily track down later on this myid==0 core.
        int i=-N1NOT1;
        int j=-N2NOT1;
        int k=-N3NOT1;
        GLOBALMACP0A3(failfloorcount,i,j,k,indexfinalstep,tscale,floor) = failfloorcountlocal[indexfinalstep][tscale][floor];
      }
    }

    if(DODISS) ENERREGIONLOOP(enerregion) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) dissreg[enerregion][dissloop]=dissreg_tot[enerregion][dissloop];
    // assume (*dissfunpos)[][] not restored since zeroed out each dump

    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=lumvsr_tot[ii];
    if(DODISSVSR) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) for(ii=0;ii<ncpux1*N1;ii++) dissvsr[dissloop][ii]=dissvsr_tot[dissloop][ii];

  }


  /////////////
  //
  // broadcast things read from header to all CPUs
  //
  /////////////
  bcast_restart_header();



  return(0);
}



///////////
///
/// old header functions:
///
///////////
#include "restart.rebeccaoldcode.c"



