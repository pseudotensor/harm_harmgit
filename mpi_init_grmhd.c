#include "decs.h"



// Initialize MPI for GRMHD code
int init_MPI_GRMHD(int *argc, char **argv[])
{


#if(USEMPI)
  fprintf(stderr, "begin: init_MPI_GRMHD\n");
  fflush(stderr);
  // init MPI (assumes nothing in set_arrays.c used here) : always done
  // non-blocking:
  init_MPI_general(argc, argv);
#else
  fprintf(stderr, "Did NOT init_MPI_GRMHD\n");
  fflush(stderr);
#endif


#if(USEOPENMP)
  init_OPENMP_general(0);
#endif


  // always do below since just sets defaults if not doing liaisonmode
  grmhd_init_mpi_liaisonmode_globalset();


#if(USEMPI)
  // this is non-blocking local operation
  MPI_Comm_rank(MPI_COMM_GRMHD, &myid); // proc id within GRMHD context only
#endif


  // currently INIT provides args to rest of processes
  myargs(*argc,*argv);


  // set default MPIid (must come after myargs())
  init_default_MPI_GRMHD_myid();


#if(USEOPENMP)
  // Setup OpenMP (just number of threads currently that was set on user arguments)
  init_OPENMP_sets_fromargs();
#endif


  // Create "myid" for HARM code not to be used in MPI functions
  // Allows translation and overloading of nodes in desirable way
#if(USEMPI)
  init_MPI_GRMHD_myid();
#endif



  // for file names
  sprintf(myidtxt, ".grmhd.%04d", MPIid[myid]);

  // rest of initialization
  init_MPI_setupfilesandgrid(*argc, *argv);

#if(USEOPENMP)
  // output to logfull_file
  init_OPENMP_general(1);
#endif


#if(DOINGLIAISON)
  // liaison-related test code:
  test_nonliaison();
#endif

  return (0);


}  


// for testing LIAISON+GRMHD code communication
void test_nonliaison(void)
{
  int myint;

  myint=myid;

#if(USEMPI&&DOINGLIAISON)
  MPI_Bcast(&myint,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD_LIAISON);
  dualfprintf(fail_file,"myid=%d myint=%d\n",myid,myint);
#endif

  myexit(0);

  
}



// Set default MPIid[] mapping
// must come after assigning rank to myid and after getting numprocs
int init_default_MPI_GRMHD_myid(void)
{
  int proc;
  
  for(proc=0;proc<numprocs;proc++){
    MPIid[proc]=proc; // True MPIid
  }


  if(MPIid[myid]!=myid){
    fprintf(stderr,"Failure to setup default MPIid[myid]=%d: myid=%d numprocs=%d\n",MPIid[myid],myid,numprocs);
    for(proc=0;proc<numprocs;proc++){
      fprintf(stderr,"MPIid[proc=%d]=%d\n",proc,MPIid[proc]);
    }
    myexit(1486754);
  }

  return(0);

}

int init_MPI_GRMHD_myid(void)
{
  int truempiid;

  //////////////////////////////////////////////////////////////
  //
  // Note, rank's only used in code in following MPI functions so far:
  //
  // MPI_Reduce() : 2nd to last is rank
  // MPI_Bcast() : 2nd to last is rank of root
  // MPI_Irecv() : 4th is rank of source
  // MPI_Isend() : 4th is rank of destination
  // MPI_Issend(): 4th is rank of destination
  // MPI_Group_incl(): 3rd is array of ranks, but only used BEFORE GRMHD ranks are created, so no mapping.  Used in mpi_init.c, so "ranks" array required to be set of mapped ranks [as it is currently].
  //  MPI_Comm_rank(): 2nd is rank, but get rank, not set or use rank
  //
  // And note that we are only mapping *ranks* NOT tags!
  //
  // Ensure all converted: (e.g.):
  // grep "MPI_Irecv" *.c *.h | grep -v MPIid | less
  //
  // DO NOT CHANGE RANK FOR:
  //
  // MPI_Type_create_darray(): 2nd is rank [rank here is rank within array geometry, which is fixed to be GRMHD "myid" rank]
  //
  //////////////////////////////////////////////////////////////

  // depending upon gridsectioning or whatever, allow user function call here that chooses some way to distribute MPI processes
  // For example, to keep L2 cache coherency require (say) 32x4x4 in 3D per *node*.  But if doing totalsize[1]=1024 and grid sectioning, then those outer procs are not used!

  // However, can create look-up table so that when MPI command is called, any given MPI process is fed myid but returns MPIid desired.  So one can redistribute (even in real-time if transfered on-grid data around) how MPI procs are arranged to overload certain nodes with (say) half of the MPI procs not being used on each node.

  // Note that we are below creating a certain direction of the mapping.  The mapping returns the MPI id for a given HARM id.  We do this direction because HARM id is setup to have specific meaning w.r.t. spatial grid.  Also, HARM id appears in more places than MPI id, that only appears for MPI commands when realling needing to know true MPI rank.

  // Note default was already set as MPIid[myid]=myid , where MPIid is true MPI id used for MPI commands.

  // true id for this proc:
  truempiid=MPIid[myid];

  // One way to distribute is to take the 2nd the grid (in real space) and doulbe it onto the 1st half.
  // So we initiate (say) 512 MPI procs using "mpirun -np 512 ./grmhd 32 4 4" and we request from the system 512 procs: (Ranger: -pe 1way 512)
  // This enables us to only have 256 processors available, but 512 processes going.  We know that half will not be used, so this is ok.
  //
  // Then we check all 256
  

  return(0);
  

}

