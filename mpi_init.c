
/*! \file mpi_init.c
     \brief Initialize MPI functions
*/


#include "decs.h"




/////////////////
///
/// General Initialize MPI
/// Only include non-blocking MPI commands
///
/////////////////
int init_MPI_general(int *argc, char **argv[])
{



#if(USEMPI)
  int ierr;
#if(USEOPENMP)
  int provided,required=MPI_THREAD_MULTIPLE;
  // lonestar4 locked-up here for some reason.   Had to set USEOPENMP->0 in makehead.inc.  Ranger was fine with openmp, but wasn't using openmp with lonestar3 with production runs, so unsure what situation is.
  ierr=MPI_Init_thread(argc, argv,required,&provided);
  stderrfprintf("Using MPI_Init_thread with required=%d and provided=%d\n",required,provided);
#else
  stderrfprintf( "Begin: MPI_Init\n"); fflush(stderr);
  ierr=MPI_Init(argc, argv);
  stderrfprintf( "End: MPI_Init\n"); fflush(stderr);
#endif



  if(ierr!=0){
    stderrfprintf("MPI Error during MPI_Init\n");
    exit(1);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &truenumprocs); // WORLD total number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_world); // WORLD proc id
  MPI_Get_processor_name(processor_name, &procnamelen); // to ensure really on certain nodes

  stderrfprintf( "WORLD proc: %d of %d on %s\n", myid_world,truenumprocs,processor_name);

  stderrfprintf( "end: init_MPI\n");
  fflush(stderr);
#else
  truenumprocs=1;
  myid=0;
  myid_world=0;
#endif


  // allocate things that are truenumprocs in size
  MPIid=(int*)malloc(sizeof(int)*truenumprocs);
  if(MPIid==NULL){
    stderrfprintf("Problem allocating memory for MPIid with truenumprocs=%d\n",truenumprocs); fflush(stderr);
    myexit(679438212);
  }


  return(0);

}



/// general initialization for OpenMP
/// Nothing required before user arguments read-in setting number of threads
int init_OPENMP_general(FILE *out)
{

  // nothing to do for OpenMP so far

  return(0);
}




/// general initialization for OpenMP
/// Nothing required before user arguments read-in setting number of threads
int init_OPENMP_sets_fromargs(void)
{

#if(USEOPENMP)
  // Set number of threads either from initialization somehow or from user arguments that set numopenmpthreads variable
  if(numopenmpthreads!=numopenmpthreadsorig) omp_set_num_threads(numopenmpthreads);
#endif

  return(0);

}




void init_MPI_setupfilesandgrid(int argc, char *argv[])
{
  int size;



  ////////////////////
  //
  // choose to combine files or not
  //
  //
  ///////////////////
  if(USEMPI){
    mpicombine = 1;    // choice
    //mpicombine=0;

    // 
    if(mpicombine){
      if(USEROMIO==0){
        // choice
        if(sortedoutput==SORTED) mpicombinetype=MPICOMBINEMINMEM;
        else if(sortedoutput==UNSORTED) mpicombinetype=MPICOMBINESIMPLE; //forced to happen since no unsorted method for the advanced combine technique
        //mpicombinetype=MPICOMBINESIMPLE; // forced for testing
      }
      else truempicombinetype=mpicombinetype=MPICOMBINEROMIO;
    }
  }
  else{
    // no choice
    mpicombine = 0;
  }




  ///////////////////
  //
  // Setup files and grid information
  // always done
  //
  ///////////////////

  // after below can use dualfprintf, and trifprintf, etc.
  init_genfiles(0);

  // after below, can access intra-CPU MPI related setup (all but boundary (i.e. inter-CPU MPI related) stuff)
  init_placeongrid_gridlocation();



  ///////////////////
  //
  // Some MPI type size checks
  //
  ///////////////////

#if(USEMPI)
  
  if(sizeof(CTYPE)==sizeof(long long int)){
    MPI_Type_size(MPI_LONG_LONG_INT,&size);
    if(size!=8){
      dualfprintf(fail_file,"size of the long long int in MPI=%d, should be 8\n",size);
      myexit(1000);
    }
  }
  
  if((sizeof(REALTYPE)==sizeof(long double))||(sizeof(SENSITIVE)==sizeof(long double))){
    MPI_Type_size(MPI_LONG_DOUBLE,&size);
    if(size!=16){
      dualfprintf(fail_file,"size of the long double in MPI=%d, should be 16\n",size);
      myexit(1000);
    }
  }
#endif


  trifprintf("done with init_MPI_setupfilesandgrid()\n");

}





/// myargs() reads in arguments for all of GRMHD codes
/// globally reads and sets ncpux1,ncpux2,ncpux3,RESTARTMODE, WHICHFILE for all CPUs
/// also globally sets numprocs and makes sure numprocs makes sense
void myargs(int argc, char *argv[])
{
  int argi,numargs,numextraargs;


  // default:
  ncpux1 = 1;
  ncpux2 = 1;
  ncpux3 = 1;

  // default unless user adds below
  numextraargs=2;
  RESTARTMODE=0;
  WHICHFILE=0;




  if(USEMPI && USEOPENMP){
    numargs=1+3+1;

    if(! (argc==numargs || argc==numargs+numextraargs) ){
      if(myid==0){
        stderrfprintf("proc: %04d : Incorrect command line: argc: %d needed at least=%d, please specify:\n",myid,argc,numargs);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> numopenmpithreads ncpux1 ncpux2 ncpux3\n",myid);
        stderrfprintf("proc: %04d : OR\n",myid);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> numopenmpithreads ncpux1 ncpux2 ncpux3 RESTARTMODE WHICHFILE\n",myid);
      }
      exit(1);
    }
    argi=1;
    numopenmpthreads=atoi(argv[argi++]);
    ncpux1=atoi(argv[argi++]);
    ncpux2=atoi(argv[argi++]);
    ncpux3=atoi(argv[argi++]);
  }
  else if(USEMPI){
    numargs=1+3;

    if(! (argc==numargs || argc==numargs+numextraargs) ){
      if(myid==0){
        stderrfprintf("proc: %04d : Incorrect command line: argc: %d needed at least=%d, please specify:\n",myid,argc,numargs);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2 ncpux3\n",myid);
        stderrfprintf("proc: %04d : OR\n",myid);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2 ncpux3 RESTARTMODE WHICHFILE\n",myid);
      }
      exit(1);
    }
    argi=1;
    ncpux1=atoi(argv[argi++]);
    ncpux2=atoi(argv[argi++]);
    ncpux3=atoi(argv[argi++]);
  }
  else if(USEOPENMP){
    numargs=1+1;

    if(! (argc==numargs || argc==numargs+numextraargs) ){
      if(myid==0){
        stderrfprintf("proc: %04d : Incorrect command line: argc: %d needed at least=%d, please specify:\n",myid,argc,numargs);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> numopenmpthreads\n",myid);
        stderrfprintf("proc: %04d : OR\n",myid);
        stderrfprintf("proc: %04d : mpirun <mpirunoptions> <progname> numopenmpthreads RESTARTMODE WHICHFILE\n",myid);
      }
      exit(1);
    }
    argi=1;
    numopenmpthreads=atoi(argv[argi++]);
  }
  else{
    numargs=1;

    if(! (argc==numargs || argc==numargs+numextraargs) ){
      if(myid==0){
        stderrfprintf("<progname>\n");
        stderrfprintf("OR\n");
        stderrfprintf("<progname> RESTARTMODE WHICHFILE\n");
      }
      exit(1);
    }
  }// end if single CPU mode

  ///////////////////////
  //
  // get user-defined restart options
  //
  ///////////////////////
  if(argc==numargs+numextraargs){
    RESTARTMODE=atoi(argv[argi++]);
    WHICHFILE=atoi(argv[argi++]);
  }


  // no failure variables inputted at run-time, since have to recompile anyways for failure tracking



  //////////////////
  //
  // General checks (that work if MPI or OpenMP not used)
  //
  /////////////////

  if(ncpux1>1 && N1==1){
    stderrfprintf("Cannot have ncpux1=%d>N1=%d\n",ncpux1,N1);
    exit(1);
  }
  if(ncpux2>1 && N2==1){
    stderrfprintf("Cannot have ncpux2=%d>N2=%d\n",ncpux2,N2);
    exit(1);
  }
  if(ncpux3>1 && N3==1){
    stderrfprintf("Cannot have ncpux3=%d>N3=%d\n",ncpux3,N3);
    exit(1);
  }


  if(numopenmpthreads>N1*N2*N3) stderrfprintf("OpenMP threads (%d) larger than total number of points (%d), so parallelization will be poor\n",numopenmpthreads,N1*N2*N3);



  
  //////////////////////
  //
  // get total number of processors
  //
  //////////////////////
  numprocs = ncpux1*ncpux2*ncpux3;

  // test
  if(sizeproclist_grmhd!=numprocs){
    stderrfprintf( "Got (sizeproclist_grmhd=%d) != (numprocs=ncpux1*ncpux2*ncpux3=%d). ncpux1=%d ncpux2=%d ncpux3=%d.  Did you run without mpirun?\n",sizeproclist_grmhd,numprocs,ncpux1,ncpux2,ncpux3);
    exit(1);
  }


  myfprintf(stderr,
            "numprocs=%d ncpux1=%d ncpux2=%d ncpux3=%d numopenmpthreads=%d :: percpusize: N1=%d N2=%d N3=%d\n",
            numprocs, ncpux1, ncpux2, ncpux3, numopenmpthreads, N1, N2,N3);


}




void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];
  char binarytype[MAXFILENAME];
  void set_binarytype(char *binarytype);


  stderrfprintf( "begin: init_genfiles ... ");
  fflush(stderr);

  if (gopp == 1) {
    strcpy(extension, PPEXT);
  } else if (gopp == 0) {
    strcpy(extension, OUTEXT);
  }
  // always have fail and general log open

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    sprintf(temps, "%s0_fail%s%s", DATADIR, extension, myidtxt);
    if ((fail_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "fail: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);
  }


  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    sprintf(temps, "%s0_log%s%s", DATADIR, extension, myidtxt);
    if ((log_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "log: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);
    fprintf(log_file, "fail_file: %d log_file: %d\n", (int)fail_file,
            (int)log_file);
    fflush(log_file);
  }


  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    sprintf(temps, "%s0_logdt%s%s", DATADIR, extension, myidtxt);
    if ((logdt_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "logdt: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);
    fflush(logdt_file);
  }





  if (myid == 0) {

    set_binarytype(binarytype);


    sprintf(temps, "%s0_logfull%s%s", DATADIR, binarytype, extension);
    if ((logfull_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "logfull: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);
    fprintf(logfull_file, "logfull_file: %d \n", (int)logfull_file);
    fflush(logfull_file);


    sprintf(temps, "%s0_logdtfull%s%s", DATADIR, binarytype, extension);
    if ((logdtfull_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "logdtfull: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);
    fprintf(logdtfull_file, "logdtfull_file: %d \n", (int)logdtfull_file);
    fflush(logdtfull_file);

    sprintf(temps, "%s0_logstep%s%s", DATADIR, binarytype, extension);
    if ((logstep_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "logstep: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);

    sprintf(temps, "%s0_logperf%s%s", DATADIR, binarytype, extension);
    if ((logperf_file = fopen(temps, "at")) == NULL) {
      stderrfprintf( "logperf: Cannot open: %s\n", temps);
      exit(1);
    }
    stderrfprintf( "opened: %s\n", temps);



  }


  // ok now
  trifprintf("end: init_genfiles\n");
}



void set_binarytype(char *binarytype)
{
  if(DOINGGRMHDTYPECODE){
    sprintf(binarytype,".grmhd");
  }
  else if(DOINGGRRAYTYPECODE){
    sprintf(binarytype,".grray");
  }
  else if(DOINGLIAISONTYPECODE){
    sprintf(binarytype,".liaison");
  }
}



/// determine where quantities are located in cent or stag sense w.r.t. direction
/// in reality only should depend upon per dimension not per direction, but eventually this info is accessed by "dir" not "dimen"
void set_primgridpos(void)
{
  int dir,pl,pliter;

  // all centered
  DIRLOOP(dir) PALLLOOP(pl) primgridpos[BOUNDPRIMTYPE][dir][pl]=primgridpos[BOUNDPRIMSIMPLETYPE][dir][pl]=CENTGRID;

  // pstag (centered except along dir)
  DIRLOOP(dir) PALLLOOP(pl){
    if(
       ((dir==X1DN || dir==X1UP) && pl==B1)
       ||((dir==X2DN || dir==X2UP) && pl==B2)
       ||((dir==X3DN || dir==X3UP) && pl==B3)
       ){
      primgridpos[BOUNDPSTAGTYPE][dir][pl]=primgridpos[BOUNDPSTAGSIMPLETYPE][dir][pl]=STAGGRID;
    }
    else{
      primgridpos[BOUNDPSTAGTYPE][dir][pl]=primgridpos[BOUNDPSTAGSIMPLETYPE][dir][pl]=CENTGRID;
    }
  }


  // all staggered (since F1 along dir=1, etc.)
  DIRLOOP(dir) PALLLOOP(pl) primgridpos[BOUNDFLUXTYPE][dir][pl]=primgridpos[BOUNDFLUXSIMPLETYPE][dir][pl]=STAGGRID;


  // vpot (centered except along perp to dir for A_i)
  DIRLOOP(dir) DLOOPA(pl){
    if(
       ((dir==X2DN || dir==X2UP || dir==X3DN || dir==X3UP) && pl==1)
       ||((dir==X1DN || dir==X1UP || dir==X3DN || dir==X3UP) && pl==2)
       ||((dir==X1DN || dir==X1UP || dir==X2DN || dir==X2UP) && pl==3)
       ){
      primgridpos[BOUNDVPOTTYPE][dir][pl]=primgridpos[BOUNDVPOTSIMPLETYPE][dir][pl]=STAGGRID;
    }
    else{
      primgridpos[BOUNDVPOTTYPE][dir][pl]=primgridpos[BOUNDVPOTSIMPLETYPE][dir][pl]=CENTGRID;
    }
  }


  // all centered
  DIRLOOP(dir) PALLLOOP(pl) primgridpos[BOUNDINTTYPE][dir][pl]=CENTGRID;


}




/// setup per-cpu parameters
void init_placeongrid_gridlocation(void)
{
  // 3's were COMPDIM, but below code is fixed to require all 3 now
  int i, j, m, l;
  int N[3 + 1];


  trifprintf("begin: init_placeongrid_gridlocation ... ");



  /////////////////
  //
  // Set CPU and grid size information
  //
  /////////////////
  N[1] = N1;
  N[2] = N2;
  N[3] = N3;

  numbercpu[1] = ncpux1;
  numbercpu[2] = ncpux2;
  numbercpu[3] = ncpux3;

  mycpupos[1]=myid%ncpux1;
  mycpupos[2]=(int)((myid%(ncpux1*ncpux2))/ncpux1);
  mycpupos[3]=(int)(myid/(ncpux1*ncpux2));



  for (m = 1; m <= COMPDIM; m++) {
    startpos[m] = mycpupos[m] * N[m];
    endpos[m] = (mycpupos[m] + 1) * N[m] - 1;

    // add up sizes for total size of grid
    totalsize[m] = 0;
    itotalsize[m] = 0;
    for (i = 0; i < numbercpu[m]; i++) {
      totalsize[m] += N[m];
      itotalsize[m] += N[m];
    }
  }

  for(m=1;m<=COMPDIM;m++){
    if((mycpupos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate mycpupos0[%d]\n",m);
      myexit(2845725);
    }
    if((startpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate startpos0[%d]\n",m);
      myexit(2845726);
    }
    if((endpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate endpos0[%d]\n",m);
      myexit(2845727);
    }
  }
  // for cpu=0 as master to rest, needs this info
  for(i=0;i<numprocs;i++){
    mycpupos0[1][i]=i%ncpux1;
    mycpupos0[2][i]=(int)((i%(ncpux1*ncpux2))/ncpux1);
    mycpupos0[3][i]=(int)(i/(ncpux1*ncpux2));

    
    for (m = 1; m <= COMPDIM; m++) {
      startpos0[m][i] = mycpupos0[m][i] * N[m];
      endpos0[m][i] = (mycpupos0[m][i] + 1) * N[m] - 1;
    }
  }



  realtotalzones = totalzones = totalsize[1] * totalsize[2] * totalsize[3];
  realtotalcompzones=realtotalzones;
  itotalzones = itotalsize[1] * itotalsize[2] * itotalsize[3];





  /////////////////
  //
  // output those things that were defined
  //
  /////////////////

#if(USEMPI)
  logfprintf("myid=%d node_name=%s procnamelen=%d\n",myid,processor_name,procnamelen);
#endif
  trifprintf("\nnumprocs(MPI)=%d ncpux1=%d ncpux2=%d ncpux3=%d numopenmpthreads=%d\n",numprocs,ncpux1,ncpux2,ncpux3,numopenmpthreads);
  trifprintf("\n Per MPI Task: N1=%d N2=%d N3=%d\n",N1,N2,N3);
  logfprintf("per: %d %d %d\n", periodicx1, periodicx2, periodicx3);

  for (m = 1; m <= COMPDIM; m++) {
    logfprintf("mycpupos[%d]: %d\n", m, mycpupos[m]);
    logfprintf( "startpos[%d]: %d\n", m, startpos[m]);
    logfprintf( "endpos[%d]: %d\n", m, endpos[m]);
    logfprintf( "totalsize[%d]: %lld\n", m, totalsize[m]);
  }

  trifprintf("totalzones: %d\n", totalzones);


  




  trifprintf("end: init_placeongrid_gridlocation\n");


}



/// setup inter-CPU grid information for message passing for domain decompsition
/// These thigns only need to be setup by first call to boundary conditions code
void init_placeongrid_griddecomposition(void)
{
  // 3's were COMPDIM, but below code is fixed to require all 3 now
  int stage,stagei,stagef;
  int i, j, m, l;

  int dir, bti;
  int opp[3*2];
  int pl,pliter;

  int numbnd[NDIM],surfa[NDIM];
  int numnpr;
  int gridpos;


  trifprintf("begin: init_placeongrid_griddecomposition ... ");



  /////////////////////
  //
  // standard interior MPI data transfer setup
  //
  /////////////////////
  for(bti=0;bti<NUMBOUNDTYPES;bti++) for(dir=0;dir<COMPDIM*2;dir++) for(j=0;j<DIRGENNUMVARS;j++){
        dirgenset[bti][dir][j]=0;
      }

  for(bti=0;bti<NUMBOUNDTYPES;bti++) for(dir=0;dir<COMPDIM*2;dir++) for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++) for(j=0;j<DIRLOOPNUMVARS;j++){
          dirloopset[bti][dir][gridpos][j]=0;
        }
  // see where this cpu needs to send/recv



  ////////////////
  //
  // set variable grid positions as CENT or STAG
  //
  ////////////////
  set_primgridpos();


  ////////////////
  //
  // set whether doing particular boundary transfer
  //
  ////////////////
  for(bti=0;bti<NUMBOUNDTYPES;bti++) { // same for any bounding type

    if(N1>1){
      // figure out left/right send/recv
      if (mycpupos[1] > 0) {
        dirgenset[bti][X1DN][DIRIF] = 1;  // do -x1 dir
      }
      if (mycpupos[1] < ncpux1 - 1) {
        dirgenset[bti][X1UP][DIRIF] = 1;  // do +x1 dir
      }
      
      // only do periodic mpi if 
      if(periodicx1&&(ncpux1>1)){
        if(mycpupos[1]==0) dirgenset[bti][X1DN][DIRIF]=1;
        else if(mycpupos[1]==ncpux1-1) dirgenset[bti][X1UP][DIRIF]=1;
      }
      

    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirgenset[bti][X1UP][DIRIF] = 0;
      dirgenset[bti][X1DN][DIRIF] = 0;
    }
    
    // figure out up/down send/recv
    if(N2>1){
      if (mycpupos[2] > 0) {
        dirgenset[bti][X2DN][DIRIF] = 1;  // -x2 dir
      }
      if (mycpupos[2] < ncpux2 - 1) {
        dirgenset[bti][X2UP][DIRIF] = 1;  // towards and from +x2 dir
      }

      // can't have both periodicx2 and SPC periodicx3      
      if(periodicx2&&(ncpux2>1)){
        if(mycpupos[2]==0) dirgenset[bti][X2DN][DIRIF]=1;
        else if(mycpupos[2]==ncpux2-1) dirgenset[bti][X2UP][DIRIF]=1;
      }
      else if(special3dspc&&(ncpux3>1)){
        // non-reflecting, transmissive boundary conditions
        if(mycpupos[2]==0) dirgenset[bti][X2DN][DIRIF]=1;
        else if(mycpupos[2]==ncpux2-1) dirgenset[bti][X2UP][DIRIF]=1;
      }
      
      
    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirgenset[bti][X2UP][DIRIF] = 0;
      dirgenset[bti][X2DN][DIRIF] = 0;
    }
    
    // figure out out/in send/recv
    if(N3>1){
      if (mycpupos[3] > 0) {
        dirgenset[bti][X3DN][DIRIF] = 1;  // -x3 dir
      }
      if (mycpupos[3] < ncpux3 - 1) {
        dirgenset[bti][X3UP][DIRIF] = 1;  // towards and from +x3 dir
      }
      
      if(periodicx3&&(ncpux3>1)){
        if(mycpupos[3]==0) dirgenset[bti][X3DN][DIRIF]=1;
        else if(mycpupos[3]==ncpux3-1) dirgenset[bti][X3UP][DIRIF]=1;
      }
      
      
    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirgenset[bti][X3UP][DIRIF] = 0;
      dirgenset[bti][X3DN][DIRIF] = 0;
    }
  }
  


  //////////////
  //
  // define which direction the other guy communicates.  Just a simple way to store this obvious fact
  //
  //////////////
  opp[X1DN]=X1UP;
  opp[X1UP]=X1DN;
  opp[X2DN]=X2UP;
  opp[X2UP]=X2DN;
  opp[X3DN]=X3UP;
  opp[X3UP]=X3DN;
  

  ////////////////////////////////////
  //
  // set which CPUs communicate to eachother 
  // same for any method (bti)
  // In setting up communication, assumes the fastest index is x1 then x2 then x3 as slowest index
  //
  ////////////////////////////////////
  for(bti=0;bti<NUMBOUNDTYPES;bti++){
    for(dir=0;dir<COMPDIM*2;dir++){
      // sets opposite direction
      dirgenset[bti][dir][DIROPP]=opp[dir];
    

      // matching CPU to transfer to/from
    
      // x1
      if((dir==X1UP)||(dir==X1DN)){
        if(ncpux1>1 &&
           (
            ((mycpupos[1]>0)&&(mycpupos[1]<ncpux1-1)) // interior CPUs
            || (mycpupos[1]==0 && dir==X1UP) // inner-CPUs pointing up
            || (mycpupos[1]==ncpux1-1 && dir==X1DN) // outer-CPUs pointing down
            )
           ){
          if(dir==X1UP) dirgenset[bti][dir][DIROTHER]=myid+1;
          if(dir==X1DN) dirgenset[bti][dir][DIROTHER]=myid-1;
        }
        else if(periodicx1 && ncpux1>1){ // here X1DN/X1UP are implicitly associated with 0/ncpux1-1
          if(mycpupos[1]==0 && dir==X1DN) dirgenset[bti][dir][DIROTHER]=myid+(ncpux1-1);  // for X1DN
          else if(mycpupos[1]==ncpux1-1 && dir==X1UP) dirgenset[bti][dir][DIROTHER]=myid-(ncpux1-1); // for X1UP
        }
      }
    
      // x2
      if((dir==X2UP)||(dir==X2DN)){
        if(ncpux2>1 &&
           (
            (mycpupos[2]>0 && mycpupos[2]<ncpux2-1) // interior CPU
            || (mycpupos[2]==0 && dir==X2UP) // exterior CPU connected to interior
            || (mycpupos[2]==ncpux2-1 && dir==X2DN) // exterior CPU connected to interior
            )
           ){
          if(dir==X2UP) dirgenset[bti][dir][DIROTHER]=myid+ncpux1;
          if(dir==X2DN) dirgenset[bti][dir][DIROTHER]=myid-ncpux1;
        }
        else if(periodicx2 && ncpux2>1){
          if(mycpupos[2]==0 && dir==X2DN) dirgenset[bti][dir][DIROTHER]=myid+(ncpux2-1)*ncpux1;
          else if(mycpupos[2]==ncpux2-1 && dir==X2UP) dirgenset[bti][dir][DIROTHER]=myid-(ncpux2-1)*ncpux1;
        }
        else if(special3dspc&&(ncpux3>1)){

          if(ncpux3%2){
            // must have ncpux3 as even
            dualfprintf(fail_file,"ncpux3=%d must be even for polar transmissive BCs\n",ncpux3);
            myexit(33676958);
          }
          // see placeongrid.nb
          // spherical polar wrapping
          int othercpupos1 = mycpupos[1];
          int othercpupos2 = mycpupos[2];
          int othercpupos3 = (mycpupos[3] + (int)ncpux3/2)%ncpux3;
          int othermyid = othercpupos1 + othercpupos2*ncpux1 + othercpupos3*ncpux1*ncpux2;
          if(mycpupos[2]==0 && dir==X2DN){
            dirgenset[bti][dir][DIROTHER] = othermyid;
            dirgenset[bti][dir][DIROPP]=X2DN; // X2DN communicates with X2DN on other CPU
          }
          else if(mycpupos[2]==ncpux2-1 && dir==X2UP){
            dirgenset[bti][dir][DIROTHER] = othermyid;
            dirgenset[bti][dir][DIROPP]=X2UP; // X2UP communicates with X2UP on other CPU
          }
        }
      }
    
      // x3
      if((dir==X3UP)||(dir==X3DN)){
        if(ncpux3>1 &&
           (
            ((mycpupos[3]>0)&&(mycpupos[3]<ncpux3-1))
            || (mycpupos[3]==0 && dir==X3UP)
            || (mycpupos[3]==ncpux3-1 && dir==X3DN)
            )
           ){
          if(dir==X3UP) dirgenset[bti][dir][DIROTHER]=myid+ncpux1*ncpux2;
          if(dir==X3DN) dirgenset[bti][dir][DIROTHER]=myid-ncpux1*ncpux2;
        }
        else if(periodicx3 && ncpux3>1){
          if(mycpupos[3]==0 && dir==X3DN) dirgenset[bti][dir][DIROTHER]=myid+(ncpux3-1)*ncpux1*ncpux2;
          else if(mycpupos[3]==ncpux3-1 && dir==X3UP) dirgenset[bti][dir][DIROTHER]=myid-(ncpux3-1)*ncpux1*ncpux2;
        }
      }

      // MPI tags that label transfer, must be unique while doing multiple transfers
      // Send and Receive tags must match expected communications (Note that first level choice of CPU transfers is choosing which CPUs, but this tag identifies which transfer among all transfers as relevant between those CPUs)
      // Normally all communications are X1DN (cpu#1) -> X1UP (cpu#2), and likewise.
      // For SPC polar transfer, we have instead X1DN (cpu#1) -> X1DN (cpu#2).  This is accounted for by overriding default DIROPP
      // send tag (set send tag as base tag identified by current CPU id (myid) and direction pointing out of grid (dir)
      dirgenset[bti][dir][DIRTAGS]= myid     * COMPDIM * 2 + dir;
      // receive tag (must match how send tag set so that receive using send tag according to other CPU)
      dirgenset[bti][dir][DIRTAGR]= dirgenset[bti][dir][DIROTHER] * COMPDIM * 2 + dirgenset[bti][dir][DIROPP];

      //////////////
      //
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are
      // like: (fromid,otherid*COMPDIM*2+*) where ? and * are
      // opposites(i.e. 0 and 2, 1 and 3, 4 and 5 or non-opposites if doing polar transmissive BCs)
      //
      //////////////

    }// end over dir
  }// end over bti



  /////////////////////
  //
  // set transfer sizes
  //
  /////////////////////
  for(bti=0;bti<NUMBOUNDTYPES;bti++){

    ///////////////////
    //
    // get number of boundary cells and number of quantities to bound
    //
    ///////////////////
    set_numbnd(bti, numbnd, &numnpr);


    // loop over directions
    for(dir=0;dir<COMPDIM*2;dir++){

      /////////////////
      //
      // enter below if doing that particular direction for this CPU
      //
      /////////////////
      if(dirgenset[bti][dir][DIRIF]){

  
        //////////////////////
        //
        // set number of variable types to transfer
        //
        //////////////////////
        if(bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE || bti==BOUNDPSTAGTYPE || bti==BOUNDPSTAGSIMPLETYPE || bti==BOUNDINTTYPE ){
          dirgenset[bti][dir][DIRNUMPR]=NPRBOUND; // not used if SPLITNPR==1 or doing general range for quantities // not used for BOUNDINTTYPE
        }
        else if(bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE){
          dirgenset[bti][dir][DIRNUMPR]=NFLUXBOUND;// not used if SPLITNPR==1 or doing general range for quantities
        }
        else if(bti==BOUNDVPOTTYPE || bti==BOUNDVPOTSIMPLETYPE){
          dirgenset[bti][dir][DIRNUMPR]=NDIM;// used
        }
        else{
          dualfprintf(fail_file,"No such bti=%d setup in set number of variable types in mpi_init.c\n",bti);
          myexit(246346769);
        }

        //////////////////////
        //
        // set transfer size
        //
        //////////////////////

        // surface area must be consistent with loops
        surfa[1]=(N2+numbnd[2]*2)*(N3+numbnd[3]*2)*N1NOT1;
        surfa[2]=(N1+numbnd[1]*2)*(N3+numbnd[3]*2)*N2NOT1;
        surfa[3]=(N1+numbnd[1]*2)*(N2+numbnd[2]*2)*N3NOT1;

        if(
           bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE
           || bti==BOUNDPSTAGTYPE || bti==BOUNDPSTAGSIMPLETYPE
           || bti==BOUNDINTTYPE
           || bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE
           || bti==BOUNDVPOTTYPE || bti==BOUNDVPOTSIMPLETYPE
           ){
          // sets size of transfer for primitive
          if((dir==X1UP)||(dir==X1DN)) dirgenset[bti][dir][DIRSIZE]=numbnd[1]*surfa[1]*numnpr;
          else if((dir==X2UP)||(dir==X2DN)) dirgenset[bti][dir][DIRSIZE]=numbnd[2]*surfa[2]*numnpr;
          else if((dir==X3UP)||(dir==X3DN)) dirgenset[bti][dir][DIRSIZE]=numbnd[3]*surfa[3]*numnpr;
        }
        // GODMARK: Why was I doing the below?
        //      else if(
        //       ){
        // // sets size of transfer for fluxes
        // // (different for "left" and "right") for flux-types
        // if(dir==X1UP) dirgenset[bti][dir][DIRSIZE]=numbnd[1]*surfa[1]*numnpr;
        // else if(dir==X1DN) dirgenset[bti][dir][DIRSIZE]=(numbnd[1]-1)*surfa[1]*numnpr;
        // else if(dir==X2UP) dirgenset[bti][dir][DIRSIZE]=numbnd[2]*surfa[2]*numnpr;
        // else if(dir==X2DN) dirgenset[bti][dir][DIRSIZE]=(numbnd[2]-1)*surfa[2]*numnpr;
        // else if(dir==X3UP) dirgenset[bti][dir][DIRSIZE]=numbnd[3]*surfa[3]*numnpr;
        // else if(dir==X3DN) dirgenset[bti][dir][DIRSIZE]=(numbnd[3]-1)*surfa[3]*numnpr;
        //      }
        else{
          dualfprintf(fail_file,"No such bti=%d setup in set transfer size in mpi_init.c\n",bti);
          myexit(246346770);
        }

      }// end if DIRIF
    }// over dir's
  }// end over bti




  /////////////////
  //
  // Below code implements the following mappings:
  //
  // CENT = centered quantities in the direction of information copying
  // STAG = staggered on faces (pstag for B^{dir} or fluxes F1 for dir=1, etc.) in direction of information copying
  //
  // For internal CPU boundary (any i=0..NBND)
  // CENT: 0+i -> N+i  and   N-1-i -> -1-i
  // STAG: 0+i -> N+i  and   N-1-i -> -1-i (assumes i=N set by first operation)
  //
  // For periodic BCs  (any i=0..NBND)
  // Note that same as interior CPU copy
  // CENT: 0+i -> N+i  and   N-1-i -> -1-i (so lower CPU controls i=0 boundary)
  // STAG: 0+i -> N+i  and   N-1-i -> -1-i (with i=N set by first operation)
  //
  // At pole in 3D (for j=0..NBND)
  // Note copying is in opposite j directions
  // CENT: 0+j -> -1-j  and   N-1-j -> N+j
  // STAG: 0+j -> -0-j  and   N-0-j -> N+j
  //
  // dirloopset[bti] sets limits on all inclusive for loops used in boundmpi.c
  //
  // In reduced dimensions (i.e. N=1), below reduces to 0 positions,
  // but above DIRIF's control whether that direction is bounded or
  // not and should take care not MPI bounding the reduced dimensions.
  // for reduced dimensions, normal boundary zones are handled by LOOPS defined in global.h
  //
  /////////////////

  ///////////////////////////////////////////////////////////////////
  //
  // first set BOUNDPRIMTYPE : CENT or FACE1,2,3 quantities that need full bounding
  //
  // BOUNDPRIMSIMPLETYPE as well
  //
  // then set BOUNDINTTYPE : CENT or FACE1,2,3 quantities that need full bounding
  //
  // set BOUNDFLUXTYPE : FACE(1/2/3) quantities -- used to only bound along flux direction (all that's needed for ENO to de-average flux using finite difference method)
  //
  // set BOUNDVPOTTYPE : perp to dir is staggered except A_0
  //
  // and BOUNDFLUXSIMPLETYPE too
  ///////////////////////////////////////////////////////////////////

  for(bti=0;bti<NUMBOUNDTYPES;bti++){

    if(bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE || bti==BOUNDPSTAGTYPE || bti==BOUNDPSTAGSIMPLETYPE || bti==BOUNDINTTYPE){
      // then can stay
    }
    else if(bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE){
      // then can stay
    }
    else if(bti==BOUNDVPOTTYPE || bti==BOUNDVPOTSIMPLETYPE){
      // then can stay
    }
    else{
      dualfprintf(fail_file,"No such bti=%d defined for dirloopset[]\n",bti);
      myexit(3468346);
    }

    ///////////////////
    //
    // get number of boundary cells and number of quantities to bound for this bti
    //
    ///////////////////
    set_numbnd(bti, numbnd, &numnpr);


    // loop over directions
    for(dir=0;dir<COMPDIM*2;dir++){


      /////////////////
      //
      // enter below if doing that particular direction for this CPU
      //
      /////////////////
      if(dirgenset[bti][dir][DIRIF]){

        // default factor by which to multiply data
        // allows simple transformations on MPI copies
        PALLLOOP(pl){
          for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){
            primfactor[bti][dir][gridpos][PACK][pl]  =1.0;
            primfactor[bti][dir][gridpos][UNPACK][pl]=1.0;
          }
        }

        // NOTEMARK: Must ensure that # of elements copied is same for PACK and UNPACK (generally numbnd[] BCs in all cases for a given direction along that direction)


        ///////////////////
        //
        // PACKING quantities (inclusive range for loops)
        //
        //////////////////
        // zones to copy from (packing -- where to copy FROM)
        if(dir==X1UP){ // right
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART1]=(N1-1)-(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRPSTOP1] =(N1-1);
          dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART1]=(N1-1)-(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRPSTOP1] =(N1-1);
          dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;
        }
        else if(dir==X1DN){ // left
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART1]=+0;
          dirloopset[bti][dir][gridpos][DIRPSTOP1] =+0+(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART1]=+0;
          dirloopset[bti][dir][gridpos][DIRPSTOP1] =+0+(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;
        }

        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X1UP)||(dir==X1DN)){
            dirloopset[bti][dir][gridpos][DIRPSTART2]=0   -numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =N2-1+numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;
            dirloopset[bti][dir][gridpos][DIRPSTART3]=0   -numbnd[3];
            dirloopset[bti][dir][gridpos][DIRPSTOP3] =N3-1+numbnd[3];
            dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;
          }
        }



        // below also correct for "if(periodicx3&&(ncpux3>1)&&ISSPCMCOORDNATIVE(MCOORD))"
        // for ISSPC, assumes "j=0" at \theta=0 and "j=N2" at \theta=\pi is copied to other CPUs.
        // mycpupos[3]<ncpux3/2 CPUs dominate others for polar value of B2, but all consistent in the end!
        // That should only affect things to machine accuracy since both poles should evolve polar B2 the same.
        //
        // Note that this special polar copy is unlike was setup, where copied j=0 and j=N2-1 effectively.
        if(special3dspc&&(ncpux3>1)&&(mycpupos[2]==0 && dir==X2DN || mycpupos[2]==ncpux2-1 && dir==X2UP) ){
          if(dir==X2UP){ // up
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=(N2-1); // inverted order
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =(N2-1)-(numbnd[2]-SHIFT2); //N2-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPDIR2]=-1;

            gridpos=STAGGRID;
            // mycpupos[3]<ncpux3/2 packs j=N2

            if(mycpupos[3]<ncpux3/2) dirloopset[bti][dir][gridpos][DIRPSTART2]=(N2-1+SHIFT2);  // inverted order // diff compared to non-pole // includes j=N2 right at pole
            else  dirloopset[bti][dir][gridpos][DIRPSTART2]=(N2-1+SHIFT2)-(SHIFT2);  // inverted order // do not pack j=N2 right at pole since will come from matching CPU

            dirloopset[bti][dir][gridpos][DIRPSTOP2] =(N2-1+SHIFT2)-(numbnd[2]-SHIFT2); //N2-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPDIR2]=-1;


            for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){
              if(bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE || bti==BOUNDPSTAGTYPE || bti==BOUNDPSTAGSIMPLETYPE ){
                // at both poles we flip signature of B2 and U2 only
                // Here we flip upon packing
                primfactor[bti][dir][gridpos][PACK][URAD1]=primfactor[bti][dir][gridpos][PACK][U1]=SIGNFLIPU1;
                primfactor[bti][dir][gridpos][PACK][B1]=SIGNFLIPB1;
                primfactor[bti][dir][gridpos][PACK][URAD2]=primfactor[bti][dir][gridpos][PACK][U2]=SIGNFLIPU2;
                primfactor[bti][dir][gridpos][PACK][B2]=SIGNFLIPB2;
                // NOTEMARK: if only interpolate U3 and B3 across pole and not \detg U3 and \detg B3 (with FLIPGDETAXIS==1), then have to flip sign across pole to avoid jump in U3 and B3 at pole.  Then, U3 and B3 will be an extremum and reduce to lower order but not have a dissipation term in the EMF-type flux calculation.
                primfactor[bti][dir][gridpos][PACK][URAD3]=primfactor[bti][dir][gridpos][PACK][U3]=SIGNFLIPU3;
                primfactor[bti][dir][gridpos][PACK][B3]=SIGNFLIPB3;
              }
              else if(bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE){
                // process sign while packing
                PALLLOOP(pl) primfactor[bti][dir][gridpos][PACK][pl]=SIGNFLIPGDET; // (e.g. gdet T^2_1. Assuming primitives are correct on active domain, then T^2_1 would be opposite signs for continuous flow through pole, but gdet has kink, so product has no kink)
                primfactor[bti][dir][gridpos][PACK][URAD2]=primfactor[bti][dir][gridpos][PACK][U2]=-SIGNFLIPGDET; // \detg T^2_2
                primfactor[bti][dir][gridpos][PACK][B2]=-SIGNFLIPGDET; // Note that F^2_{B2) = 0, so doesn't matter, but maintain consistency
                primfactor[bti][dir][gridpos][PACK][URAD3]=primfactor[bti][dir][gridpos][PACK][U3]=-SIGNFLIPGDET; // \detg T^2_3 is like \detg T^2_2 as far as sign if don't interpolate \detg U3 and \detg B3 across pole.
                primfactor[bti][dir][gridpos][PACK][B3]=-SIGNFLIPGDET; // F^2_{B3) like T^2_3 like T^2_2
                // No need to handle T^2_3 
              }
              else if(bti==BOUNDVPOTTYPE || bti==BOUNDVPOTSIMPLETYPE){
                // flip while packing
                DLOOPA(pl) primfactor[bti][dir][gridpos][PACK][pl]=-SIGNFLIPGDET; // A_1 A_3 : These point in 1 and 3 directions like scalars, but gdet-compressed at pole with kink, so need to unkink
                primfactor[bti][dir][gridpos][PACK][2]=SIGNFLIPGDET; // A_2 (points into axis but with gdet, so as if gdet*B2)
              }
            }// end over gridpos
          }
          else if(dir==X2DN){ // down
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=0;
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =0+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;

            gridpos=STAGGRID;
            // mycpupos[3]<ncpux3/2 packs j=0
            if(mycpupos[3]<ncpux3/2) dirloopset[bti][dir][gridpos][DIRPSTART2]=0; // includes j=0 right at pole
            else dirloopset[bti][dir][gridpos][DIRPSTART2]=SHIFT2; // doesn't includes j=0 right at pole

            dirloopset[bti][dir][gridpos][DIRPSTOP2] =0+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;
          }
        }
        else{
          // Old treatment of pole
          if(dir==X2UP){ // up
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=(N2-1)-(numbnd[2]-SHIFT2); //N2-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =(N2-1);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;

            gridpos=STAGGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=(N2-1)-(numbnd[2]-SHIFT2); //N2-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =(N2-1);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;
          }
          else if(dir==X2DN){ // down
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=0;
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =0+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;

            gridpos=STAGGRID;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=0;
            dirloopset[bti][dir][gridpos][DIRPSTOP2] =0+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;
          }
        }

        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X2UP)||(dir==X2DN)){
            dirloopset[bti][dir][gridpos][DIRPSTART1]=-numbnd[1];
            dirloopset[bti][dir][gridpos][DIRPSTOP1]=N1-1+numbnd[1];
            dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;
            dirloopset[bti][dir][gridpos][DIRPSTART3]=-numbnd[3];
            dirloopset[bti][dir][gridpos][DIRPSTOP3]=N3-1+numbnd[3];
            dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;
          }
        }



        if(dir==X3UP){ // up
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART3]=(N3-1)-(numbnd[3]-SHIFT3); //N3-numbnd[3];
          dirloopset[bti][dir][gridpos][DIRPSTOP3] =(N3-1);
          dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART3]=(N3-1)-(numbnd[3]-SHIFT3); //N3-numbnd[3];
          dirloopset[bti][dir][gridpos][DIRPSTOP3] =(N3-1);
          dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;
        }
        else if(dir==X3DN){ // down
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART3]=+0;
          dirloopset[bti][dir][gridpos][DIRPSTOP3] =+0+(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRPSTART3]=+0;
          dirloopset[bti][dir][gridpos][DIRPSTOP3] =+0+(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRPDIR3]=+1;
        }

        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X3UP)||(dir==X3DN)){
            dirloopset[bti][dir][gridpos][DIRPSTART1]=-numbnd[1];
            dirloopset[bti][dir][gridpos][DIRPSTOP1]=N1-1+numbnd[1];
            dirloopset[bti][dir][gridpos][DIRPDIR1]=+1;
            dirloopset[bti][dir][gridpos][DIRPSTART2]=-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPSTOP2]=N2-1+numbnd[2];
            dirloopset[bti][dir][gridpos][DIRPDIR2]=+1;
          }
        }






        ///////////////////
        //
        // UNPACKING quantities
        //
        //////////////////
        // zones to copy into (unpacking -- where to copy INTO)

        // x1
        if(dir==X1UP){ // right
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART1]=(N1-1+SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUSTOP1] =(N1-1+SHIFT1)+(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART1]=(N1-1+SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUSTOP1] =(N1-1+SHIFT1)+(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;
        }
        else if(dir==X1DN){ // left
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART1]=-SHIFT1-(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUSTOP1] =-SHIFT1;
          dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART1]=-SHIFT1-(numbnd[1]-SHIFT1);
          dirloopset[bti][dir][gridpos][DIRUSTOP1] =-SHIFT1;
          dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;
        }
        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X1UP)||(dir==X1DN)){
            dirloopset[bti][dir][gridpos][DIRUSTART2]=-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRUSTOP2]=N2-1+numbnd[2];
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;
            dirloopset[bti][dir][gridpos][DIRUSTART3]=-numbnd[3];
            dirloopset[bti][dir][gridpos][DIRUSTOP3]=N3-1+numbnd[3];
            dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;
          }
        }

        // normal bc has 0,1,2 -> N,N+1,N+2 for both pcent and pstag
        // for "if(periodicx3&&(ncpux3>1)&&ISSPCMCOORDNATIVE(MCOORD))" below gives: 0,1,2 -> 0,-1,-2 for pstag and 0,1,2 -> -1,-2,-3 for pcent
        // When inverting v^\theta -> -v^\theta, this effectively inverts order in x2 direction too
        // x2
        if(special3dspc&&(ncpux3>1)&&(mycpupos[2]==0 && dir==X2DN || mycpupos[2]==ncpux2-1 && dir==X2UP)){
          if(dir==X2UP){ // up
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=(N2-1+SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =(N2-1+SHIFT2)+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;

            gridpos=STAGGRID;
            // mycpupos[3]<ncpux3/2 packs j=N2,0, so mycpupos[3]>=ncpux3/2 unpacks j=N2,0
            if(mycpupos[3]<ncpux3/2) dirloopset[bti][dir][gridpos][DIRUSTART2]=(N2-1+SHIFT2)+(SHIFT2);
            else dirloopset[bti][dir][gridpos][DIRUSTART2]=(N2-1+SHIFT2); // includes from j=N2

            dirloopset[bti][dir][gridpos][DIRUSTOP2] =(N2-1+SHIFT2)+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;
          }
          else if(dir==X2DN){ // down
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=-SHIFT2; // inverted order compared to packing
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =-SHIFT2-(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=-1;

            gridpos=STAGGRID;
            // mycpupos[3]<ncpux3/2 packs j=N2,0, so mycpupos[3]>=ncpux3/2 unpacks j=N2,0
            if(mycpupos[3]<ncpux3/2) dirloopset[bti][dir][gridpos][DIRUSTART2]=-SHIFT2; // inverted order compared to packing
            else dirloopset[bti][dir][gridpos][DIRUSTART2]=-0; // diff due to pole copy // inverted order compared to packing

            dirloopset[bti][dir][gridpos][DIRUSTOP2] =-0-(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=-1;

            for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){
              if(bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE || bti==BOUNDPSTAGTYPE || bti==BOUNDPSTAGSIMPLETYPE ){
                // at both poles we flip signature of B2 and U2 only
                // Here we flip upon unpacking
                primfactor[bti][dir][gridpos][UNPACK][URAD1]=primfactor[bti][dir][gridpos][UNPACK][U1]=SIGNFLIPU1;
                primfactor[bti][dir][gridpos][UNPACK][B1]=SIGNFLIPB1;
                primfactor[bti][dir][gridpos][UNPACK][URAD2]=primfactor[bti][dir][gridpos][UNPACK][U2]=SIGNFLIPU2;
                primfactor[bti][dir][gridpos][UNPACK][B2]=SIGNFLIPB2;
                // again assumes U3 and B3 interpolated across pole and not \detg U3 and \detg B3
                primfactor[bti][dir][gridpos][UNPACK][URAD3]=primfactor[bti][dir][gridpos][UNPACK][U3]=SIGNFLIPU3;
                primfactor[bti][dir][gridpos][UNPACK][B3]=SIGNFLIPB3;
              }
              else if(bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE){
                // Here we flip upon unpacking
                PALLLOOP(pl) primfactor[bti][dir][gridpos][UNPACK][pl]=SIGNFLIPGDET; // gdet T^2_1, so like gdet*B2, and kink avoided if don't flip sign since B2 standard in active domains with sign change itself in active domains.
                // override for symmetric quantities
                primfactor[bti][dir][gridpos][UNPACK][URAD2]=primfactor[bti][dir][gridpos][UNPACK][U2]=-SIGNFLIPGDET; // \detg T^2_2 , avoid kink must flip sign
                primfactor[bti][dir][gridpos][UNPACK][B2]=-SIGNFLIPGDET; // Note that F^2_{B2) = 0, so doesn't matter, but maintain consistency
                primfactor[bti][dir][gridpos][UNPACK][URAD3]=primfactor[bti][dir][gridpos][UNPACK][U3]=-SIGNFLIPGDET; // \detg T^2_3 like \detg T^2_3 if U3&B3 (not \detg U3&B3) interpolated
                primfactor[bti][dir][gridpos][UNPACK][B3]=-SIGNFLIPGDET; // F^2_{B3) like T^2_3 like T^2_2
              }
              else if(bti==BOUNDVPOTTYPE || bti==BOUNDVPOTSIMPLETYPE){
                // Here we flip upon unpacking
                DLOOPA(pl) primfactor[bti][dir][gridpos][UNPACK][pl]=-SIGNFLIPGDET; // A_0 A_1 A_3 like scalars, but compressed by gdet.  So flip sign so no kink
                primfactor[bti][dir][gridpos][UNPACK][2]=SIGNFLIPGDET; // A_2 like gdet B2.  A_2 will have opposite sign across pole in active domains, but gdet + in both, so avoid flipping so that A_2 has no kink at pole.
              }
            }// end over gridpos

          }// end dir==X2DN
        }
        else{
          // Old treatment of pole
          if(dir==X2UP){ // up
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=(N2-1+SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =(N2-1+SHIFT2)+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;

            gridpos=STAGGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=(N2-1+SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =(N2-1+SHIFT2)+(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;
          }
          else if(dir==X2DN){ // down
            gridpos=CENTGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=-SHIFT2-(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =-SHIFT2;
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;

            gridpos=STAGGRID;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=-SHIFT2-(numbnd[2]-SHIFT2);
            dirloopset[bti][dir][gridpos][DIRUSTOP2] =-SHIFT2;
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;
          }
        }

        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X2UP)||(dir==X2DN)){
            dirloopset[bti][dir][gridpos][DIRUSTART1]=-numbnd[1];
            dirloopset[bti][dir][gridpos][DIRUSTOP1]=N1-1+numbnd[1];
            dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;
            dirloopset[bti][dir][gridpos][DIRUSTART3]=-numbnd[3];
            dirloopset[bti][dir][gridpos][DIRUSTOP3]=N3-1+numbnd[3];
            dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;
          }
        }


        // x3
        if(dir==X3UP){ // up
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART3]=(N3-1+SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUSTOP3] =(N3-1+SHIFT3)+(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART3]=(N3-1+SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUSTOP3] =(N3-1+SHIFT3)+(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;
        }
        else if(dir==X3DN){ // down
          gridpos=CENTGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART3]=-SHIFT3-(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUSTOP3] =-SHIFT3;
          dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;

          gridpos=STAGGRID;
          dirloopset[bti][dir][gridpos][DIRUSTART3]=-SHIFT3-(numbnd[3]-SHIFT3);
          dirloopset[bti][dir][gridpos][DIRUSTOP3] =-SHIFT3;
          dirloopset[bti][dir][gridpos][DIRUDIR3]=+1;
        }

        for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++){ // stag and cent same for off-dir directions.  Both are equivalent to CENTGRID
          if((dir==X3UP)||(dir==X3DN)){
            dirloopset[bti][dir][gridpos][DIRUSTART1]=-numbnd[1];
            dirloopset[bti][dir][gridpos][DIRUSTOP1]=N1-1+numbnd[1];
            dirloopset[bti][dir][gridpos][DIRUDIR1]=+1;
            dirloopset[bti][dir][gridpos][DIRUSTART2]=-numbnd[2];
            dirloopset[bti][dir][gridpos][DIRUSTOP2]=N2-1+numbnd[2];
            dirloopset[bti][dir][gridpos][DIRUDIR2]=+1;
          }
        }


      }// end if DIRIF
    }// end DIRLOOP
  }// end bti loop


















  /////////////////
  //
  // output those things that were defined
  //
  /////////////////

  for(bti=0;bti<NUMBOUNDTYPES;bti++) {
    for (m = 0; m < COMPDIM*2; m++) {
      for(l = 0 ; l < DIRGENNUMVARS ; l++) {
        logfprintf( "dirgenset[%d][%d][%d]: %d\n", bti, m, l, dirgenset[bti][m][l]);
      }
    }
  }

  for(bti=0;bti<NUMBOUNDTYPES;bti++) {
    for (m = 0; m < COMPDIM*2; m++) {
      for(gridpos=0;gridpos<NUMPRIMGRIDPOS;gridpos++) {
        for(l = 0 ; l < DIRLOOPNUMVARS ; l++) {
          logfprintf( "dirloopset[%d][%d][%d][%d]: %d\n", bti, m, gridpos, l, dirloopset[bti][m][gridpos][l]);
        }
      }
    }
  }





  /////////////////// 
  //
  // Setup supermpi method (not working right now)
  //
  ///////////////////


#if(SIMULBCCALC!=-1)
  // this is definitely not setup for 3D, and never fully worked...still interesting idea.

  if(SIMULBCCALC==2){
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}
    
    if(SIMULBCCALC>=1){
      for(stage=stagei;stage<=stagef;stage++){
        STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
        logfprintf("CZLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(0,N1,-1,N2,isc,iec,jsc,jec);
        logfprintf("F1LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(-1,N1,0,N2,isc,iec,jsc,jec);
        logfprintf("F2LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(0,N1,0,N2,isc,iec,jsc,jec);
        logfprintf("EMFLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(0,N1,0,N2-1,isc,iec,jsc,jec);
        logfprintf("F1CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(0,N1-1,0,N2,isc,iec,jsc,jec);
        logfprintf("F2CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(-1,N1,-1,N2,isc,iec,jsc,jec);
        logfprintf("DQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        STAGECONDITION(-2,N1+1,-2,N2+1,isc,iec,jsc,jec);
        logfprintf("PREDQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
        logfprintf("\n");
      }
    }
  }
#endif





  trifprintf("end: init_placeongrid_griddecomposition\n");
}




int myexit(int call_code)
{
  int i, j, k, l;
  int cleanfinish,dofaildump;
  FILE *faildump;
  char mysys[MAXFILENAMELONG];
  char binarytype[MAXFILENAME];
  void set_binarytype(char *binarytype);
  int inparallel,tid;



  trifprintf("proc: %s : Exiting cc: %d nstep: %ld\n", myidtxt, call_code, nstep);


#if(USEOPENMP)
  if(omp_in_parallel()){
    inparallel=1;
    tid=omp_get_thread_num();
  }
  else{
    tid=0;
    inparallel=0;    
  }
#else
  tid=0;
  inparallel=0;
#endif


  if(tid==0){ // only 1 thread does below


#if(MAILWHENDONE && !MPIAVOIDFORK)
    if(myid==0){


      set_binarytype(binarytype);

      sprintf(mysys,"echo \"%s : done with `pwd`, inparallel=%d\" > done%s.txt",EMAILMESSAGE,inparallel,binarytype);
      system(mysys);
      if(MAILFROMREMOTE){
        sprintf(mysys,"scp done.txt %s ; ssh %s \"mail %s < done%s.txt\"",REMOTEHOST,REMOTEHOST,EMAILADDRESS,binarytype);
        system(mysys);
      }
      else{
        sprintf(mysys,"mail %s < done%s.txt",EMAILADDRESS,binarytype);
        system(mysys);
      }
    }
#endif

    dofaildump=0;
    if (call_code > 0) {
      stderrfprintf(
                    "proc: %s : Failure.  Please check failure file: cc: %d\n",
                    myidtxt, call_code);

      if(call_code<ERRORCODEBELOWCLEANFINISH) cleanfinish = 1;
      else cleanfinish=0; // assume this means dump procedure failed, so don't get into infinite failure loop
      // should never have non-clean finish, but sometimes do have them in code, but not marked right now
      if(cleanfinish) dofaildump=1;
      if(!cleanfinish){
#if(USEMPI)
        // must abort since no clear to communicate to other cpus now
        MPI_Abort(MPI_COMM_GRMHD, 1);
#endif
      }
    }
    else{
      dofaildump=0;
      cleanfinish=1;
    }



    if (dofaildump) {
      stderrfprintf( "proc: %s : dumping failure dump with callcode=2\n",
                     myidtxt);

      // assume want previous timestep data, not bad just-computed
      // data\n");
      // now diag should not fail if last timestep was non-fail type
      if (DODIAGS) diag(FINAL_OUT,t,nstep,realnstep);
    }


    ////////////////////////////////
    //
    // first cleanup any prior MPI non-blocking calls by making fake write call
    //
    ////////////////////////////////
#if(USEMPI)
    if(cleanfinish && USEROMIO==1 && MPIVERSION==2 ){
      fakedump(0);
    }
#endif



    ////////////////////////////////
    //
    // must close AFTER diag()
    //
    ////////////////////////////////
    if (call_code >= 0) {
      if (fail_file) fclose(fail_file);
      if (log_file) fclose(log_file);
      myfclose(&logfull_file,"Can't close logfull_file\n");
    }



    if(cleanfinish){
      stderrfprintf( "Ending Computation on proc: %s, holding for other cpus\n", myidtxt);
    }


    myfprintf(stderr, "Ended Computation on all processors\n");
    //final_myexit(); // Don't want to Abort if don't have to
    stderrfprintf( "END\n");
    fflush(stderr);
    exit(0);



  }// end if master thread



  return (0);
}



#if(PRODUCTION<=1)

/// note, this may be called in different locations of the code by
/// different CPUs
int error_check(int wherefrom)
{
  int i, j, k;
  int errorsend = 0;
  // check if error exists and exit if so

  if (failed > 0) {
    dualfprintf(fail_file,
                "Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld steppart=%d :: t: %21.15g wherefrom = %d\n",
                myid, failed, nstep, realnstep, steppart, t,wherefrom);
  }

  if (numprocs > 1) {
    errorsend = failed;
#if(USEMPI)
    // dualfprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
    MPI_Allreduce(&errorsend, &failed, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_GRMHD);
    // dualfprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
#endif
  }
  if (failed > 0) {
    dualfprintf(fail_file,
                "Result: Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld steppart=%d :: t: %21.15g\n",
                myid, failed, nstep, realnstep, steppart, t);
    // control behavior of failure here (i.e. could return(1) and
    // continue or something)
    // if(failed==1) myexit(1);
    // if(failed==2) myexit(1);
    // if(failed==3) myexit(1);
    myexit(wherefrom);
    return (1);
  }
  return (0);
}

#endif





/// just copied from pnmhd code
#if(0)
void init_MPIgroup(void)
{
  int *ranks;
  int i,j,k,numranks;


  // allocate things that are truenumprocs in size
  ranks=(int*)malloc(sizeof(int)*truenumprocs);
  if(ranks==NULL){
    stderrfprintf("Problem allocating memory for ranks with truenumprocs=%d\n",truenumprocs); fflush(stderr);
    myexit(915213756);
  }


  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

  stderrfprintf("begin: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);

  // x1-inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==0){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    stderrfprintf("problem with inner x1-group: numranks: %d ncpux2: %d ncpux3: %d ncpux2*ncpux3: %d\n",numranks,ncpux2,ncpux3,ncpux2*ncpux3);
    myexit(97834683);
  }
  // now ranks holds inner x1 boundary of cpus, and numranks holds number of such ranks

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[2]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[2], &combound[2]); 

  // x1 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==ncpux1-1){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    stderrfprintf("problem with outer x1-group: numranks: %d ncpux2*ncpux3: %d\n",numranks,ncpux2*ncpux3);
    myexit(92787621);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[0]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[0], &combound[0]); 

  // x2 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==0){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    stderrfprintf("problem with inner x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(83649545);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[1]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[1], &combound[1]); 

  // x2 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==ncpux2-1){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    stderrfprintf("problem with outer x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(28364888);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[3]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[3], &combound[3]); 



  // x3 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==0){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    stderrfprintf("problem with inner x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(18758365);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[5]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[5], &combound[5]); 

  // x3 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==ncpux3-1){
      ranks[j]=MPIid[i];
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    stderrfprintf("problem with outer x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(29776546);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[4]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[4], &combound[4]); 


  free(ranks);

  // 0: right 1: up 2: left 3: down 4: out 5: in(as in bound.c)

  // when using these communicators, must make sure the call to a communication using it isn't done by a non-member cpu!
  // (above: stupid, I know, should just skip if non-member cpu tries a function)
  stderrfprintf("end: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);
}
#endif







// Environment variables (in bash, set like: export <variable>=<value> )
//
// OMP_SCHEDULE : Which schedule method is set (e.g. export OMP_SCHEDULE="guided, 4")
//
// OMP_NUM_THREADS : Sets the maximum number of threads to use during execution (e.g. export OMP_NUM_THREADS=8)
//
// OMP_DYNAMIC :  Enables or disables dynamic adjustment of the number of threads available for execution of parallel regions. Valid values are TRUE or FALSE. (e.g. export OMP_DYNAMIC=TRUE )
//
// OMP_NESTED : Enables or disables nested parallelism. Valid values are TRUE or FALSE. (e.g. export OMP_NESTED=TRUE )
//   Implementation notes:        * Your implementation may or may not support nested parallelism and/or dynamic threads. If nested parallelism is supported, it is often only nominal, in that a nested parallel region may only have one thread.  * Consult your implementation's documentation for details - or experiment and find out for yourself if you can't find it in the documentation. 
//
// OMP_STACKSIZE    New with OpenMP 3.0. Controls the size of the stack for created (non-Master) threads. Examples:
//    setenv OMP_STACKSIZE 2000500B
//    setenv OMP_STACKSIZE "3000 k "
//    setenv OMP_STACKSIZE 10M
//    setenv OMP_STACKSIZE " 10 M "
//    setenv OMP_STACKSIZE "20 m "
//    setenv OMP_STACKSIZE " 1G"
//    setenv OMP_STACKSIZE 20000
//
// Default is about 4-8MB on modern systems.
//
// OMP_WAIT_POLICY     New with OpenMP 3.0. Provides a hint to an OpenMP implementation about the desired behavior of waiting threads. A compliant OpenMP implementation may or may not abide by the setting of the environment variable. Valid values are ACTIVE and PASSIVE. ACTIVE specifies that waiting threads should mostly be active, i.e., consume processor cycles, while waiting. PASSIVE specifies that waiting threads should mostly be passive, i.e., not consume processor cycles, while waiting. The details of the ACTIVE and PASSIVE behaviors are implementation defined. Examples:
//    setenv OMP_WAIT_POLICY ACTIVE
//    setenv OMP_WAIT_POLICY active
//    setenv OMP_WAIT_POLICY PASSIVE
//    setenv OMP_WAIT_POLICY passive
//
//OMP_MAX_ACTIVE_LEVELS    New with OpenMP 3.0. Controls the maximum number of nested active parallel regions. The value of this environment variable must be a non-negative integer. The behavior of the program is implementation defined if the requested value of OMP_MAX_ACTIVE_LEVELS is greater than the maximum number of nested active parallel levels an implementation can support, or if the value is not a non-negative integer. Example:
//   setenv OMP_MAX_ACTIVE_LEVELS 2
//
// OMP_THREAD_LIMIT   New with OpenMP 3.0. Sets the number of OpenMP threads to use for the whole OpenMP program. The value of this environment variable must be a positive integer. The behavior of the program is implementation defined if the requested value of OMP_THREAD_LIMIT is greater than the number of threads an implementation can support, or if the value is not a positive integer. Example:
//   setenv OMP_THREAD_LIMIT 8 


// To enable OpenMP:
// Compiler  Flag
///------------------------
// IBM          -qsmp=omp
// Intel  -openmp
// PathScale  -mp
// PGI          -mp
// GNU          -fopenmp




/// get number of OpenMP threads expected to operate in real pragma calls
void get_report_openmp_thread_info(FILE * out)
{
  int tid;

#if(USEOPENMP) // need inside cpp conditional since using OpenMP functions, not just directives.

#pragma omp parallel private(tid)
  {
    // Obtain and print thread id
    fprintf(out,"proc: %d : Thread = %d activated\n", myid,omp_get_thread_num());

    tid = omp_get_thread_num();
    if (tid == 0){
      // Only master thread does this
      numopenmpthreadsorig = omp_get_num_threads();

#if(OPENMPVERSION==3)
      fprintf(out,"OpenMP 3.0 activated: Maximum number of threads available=%d\n",omp_get_thread_limit());
#endif

      fprintf(out,"Master MPI proc=%d reports: Number of threads originaly=%d out of maximum=%d out of procs=%d\n", myid, numopenmpthreadsorig,omp_get_max_threads(),omp_get_num_procs());
      if(omp_get_dynamic()){
        fprintf(out,"Dynamic thread adjustment is enabled\n");
      }
      else{
        fprintf(out,"Dynamic thread adjustment is disabled\n");
      }

      // Use omp_set_nested() to enable if desired

      if(omp_get_nested()){
        fprintf(out,"Nested parallelism is enabled, so allowed\n");
      }
      else{
        fprintf(out,"Nested parallelism is disabled, so is NOT allowed\n");
      }

      // Note that when numopenmpthreads is requested in the next parallel region, all threads will have the value since the numopenmpthreads value is shared.
    }


  }// end parallel region

#endif

}





