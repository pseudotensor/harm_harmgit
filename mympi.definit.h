
/*! \file mympi.definit.h
     \brief Default MPI definitions/macros
*/



// 1 = MPI-1
// 2 = MPI-2 (for non-blocking MPI-IO)
#define MPIVERSION 2 // choice

// OpenMP spec version
#define OPENMPVERSION 2 // 2 or 3 (3.0 not often implemented) 

// whether to use ROMIO and whether to avoid fork()
#if(USEMPI==1)
#include <mpi.h>
// can't use ROMIO unless file system is shared across all CPUs (i.e. not BH)
// ROMIO for tungsten since problems with minmem method.
// Some systems problem with ROMIO as happens to some people: File locking failed in ADIOI_Set_lock.
// TACC's Lonestar can NOT do ROMIO properly (no file locking)
#define USEROMIO  1   // choice, whether to use ROMIO parallel I/O package
// below comes from compiler so tied to machine's MPI setup type
#define MPIAVOIDFORK (USINGMPIAVOIDFORK)   // choice (avoids system/fork/etc calls)
#else
#define USEROMIO  0   // no choice
#define MPIAVOIDFORK 0   // always 0, can't have GM without MPI
#endif


// whether to put barrier before ROMIO write/read start so that avoids large unexpected messages on (e.g.) Cray Kraken NICS.
#define BARRIERROMIOPRE 1

// see boundmpi.c sendrecv() for details
// 0 : no strong flow control
// 1 : simple handshake to ensure recv posts
// 2 : post recv before computations so very likely all recv's posted before send. [NOT YET -- requires careful handling of which boundary calls repeat and also need workbc space for each type want to have pre-post recv's.]
#define MPIFLOWCONTROL 0




// whether to simultaneously compute and transfer bc (i.e. actually use non-blocking with a purpose).
// -1: use old loop even, super NONONONO!
// 0: no
// 1: yes with if type loops over fewer loops
// 2: yes with more loops without if, more blocks
#define SIMULBCCALC -1


// first 1 is choice, but no choice when USINGOPENMP==1
#define MPIEQUALNONMPI (1 || USINGOPENMP)
// 0= relax condition for MPI to equal non-MPI
// 1= guarentee not only that MPI boundaries are transfered but that an MPI run is identical to non-MPI run.


// this is the total stencil half width, where sts is full stencil size as: (sts-1)/2 which represents the one-sided safetey size in number of zones
// 2 even for parabolic, due to magnetic field stencil!
#define SAFESIZE (2)

