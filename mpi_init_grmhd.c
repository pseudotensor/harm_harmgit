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
  init_OPENMP_general(stderr);
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
  // report MPIid[myid] ordering
  report_myid(stderr);


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

  // report MPIid[myid] ordering
  if(myid==0) report_myid(logfull_file);
  report_myid(log_file);


#if(USEOPENMP)
  // output to logfull_file
  if(myid==0) init_OPENMP_general(logfull_file);
  init_OPENMP_general(log_file);
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

  // One way to distribute is to take the 2nd the grid (in real space) and doulbe it onto the 1st half.


  /////////////////////////////////////////
  //
  // Take Ranger for example:
  // 
  // http://www.tacc.utexas.edu/services/userguides/ranger/
  // 
  // When reordering tasks, one should ensure to overload each core in the correct way with the correct socket and core affinity for a given MPI task.
  // 
  // Ranger has 16 cores/node. One sets:
  // 
  // -pe <tpn>way <nonx16>
  // 
  // where <tpn> is MPI tasks per node and NoN=Number of nodes. Then we run the HARM and choose:
  // 
  // ./grmhd OMPPERTASK ncpux1 ncpux2 ncpux3
  // 
  // Also note that on Ranger, each socket has its own memory! I didn't know this. The documentation pages also explain how to assign socket/core affinity so that memory accesses are more efficient. This is important, for example, in order to avoid a thread or task on one socket trying to use or allocate memory on another socket since that would use the PCI bus that is very slow. Also, when overloading a socket, we want to have MPI tasks assigned with the correct affinity to a socket -- otherwise we'll have NO gain!
  // 
  // So now, let us set:
  // 
  // -pe 8way 256
  // 
  // and one sets (although code can override this and the affinity file mentioned below doesn't use it):
  // 
  // export OMP_NUM_THREADS=4
  // 
  // and run using:
  // 
  // ./grmhd 4 8 2 1
  // 
  // Ranger will allocate 256/16=16 nodes with 8 MPI tasks per node. If we run using:
  // 
  // ibrun tacc_affinity ./grmhd 4 8 2 1
  // 
  // then the affinity will be ensured to use round robin to assign MPI tasks to sockets. That is, it will assign (for 8way) for the first node 8 MPI tasks with:
  // 
  // rank0 -> Socket0 rank1 -> Socket0 rank2 -> Socket1 rank3 -> Socket1 rank4 -> Socket2 rank5 -> Socket2 rank6 -> Socket3 rank7 -> Socket4
  // 
  // That's round robin with 8way set. Clearly, then, for each node with 4 sockets we must KNOW that this is infact the ordering. If this is the ordering, then we can take those 8 MPI tasks and ensure that rank1,3,5,7 are associated with the outer radial regions beyond the initially active grid sectioning. Only then will each socket be able to launch 4 OpenMP threads that efficiently stay on a single socket for each MPI task.
  // 
  // BTW, one can write your own ordering of affinity very easily. The script is just /share/sge/default/pe_scripts/tacc_affinity and appears to be very easy to control. However, let's assume the default "tacc_affinity" affinity method and modify our code as below.
  // 
  // Within the code, we reorder the MPIid[]'s (see new code and mpi_init_grmhd.c) so that (e.g. for 1D for demonstration purposes) MPI tasks (ranks) 0-15 are instead ordered as 0,2,4,6,8,10,12,14 in physical model space. This assumes that the MPI has done a good job of affinity itself, which on TACC is true.
  // 
  // Then assuming grid sectioning operates on roughly half the grid at any one time, then MPI task 0 and 1 won't be going at the same time. All MPI tasks will use 4 cores per task so they only access their own memory channel and don't go across the PCI bus on Ranger.
  // 
  // So this reordering of MPI task numbers depends critically (at least for multi-socket systems where each socket accesses different memory) on knowing how the default affinity is defined. And it depended upon that one chose "8way" but tricked the system into a setup that really means "4way" since ultimately we used 4 OpenMP threads per MPI task.
  //
  //////////////////////////////



  // true id for this proc:
  //  truempiid=MPIid[myid];


  // Call user function
  theproblem_set_myid();

  

  return(0);
  

}


// report listings of MPIid[myid]
int report_myid(FILE *out)
{
  int ranki,rankj,rankk,origid;

  fprintf(out,"BEGIN Rank orders in physical model space\n");

  fprintf(out,"\n");
  for(rankk=0;rankk<ncpux3;rankk++){
    if(rankk==0) fprintf(out,"rankk=%d::\n",rankk);
    for(rankj=0;rankj<ncpux2;rankj++){
      for(ranki=0;ranki<ncpux1;ranki++){
	origid=ranki + rankj*ncpux1 + rankk*ncpux1*ncpux2;
	fprintf(out,"%04d",MPIid[origid]);
	if(ranki!=ncpux1-1) fprintf(out," ");
	else fprintf(out,"\n");
      }
      if(rankj==ncpux2-1) fprintf(out,"\n");
    }
    if(rankk==ncpux3-1) fprintf(out,"\n");
  }
  fprintf(out,"\n");

  fprintf(out,"END Rank orders in physical model space\n");

  return(0);
}


