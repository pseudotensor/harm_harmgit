
/*! \file mympi.h
     \brief General MPI wrapper for macros/definitions

   To start mpi run (without openmp), first do:

   mpirun -np 4 ./grmhd

   e.g. 4 cpus using mpi:

   rm nohup.out ; nohup sh -c 'mpirun -np 4 ./grmhd' &

   Note: cannot have any cpu write to same file or pipe to same file
   without parallel I/O or blocking through each CPU

// TO GET MPI GOING
// 0) Choose system in makefile (e.g. USEORANGE=1)
// 1) set USEMPI 1 in makehead.inc
// 2) set USEMPI 1 below
// 3) ensure FULLOUTPUT 0
// 4) Choose file output method (e.g. ROMIO or not -- see also mpi_init.c)
// 4) make sure MAXBND set so not excessive (fail file will report warning)
// 5) recompile
// 6) submit job (see batch.*) or run mpirun directly
*/





/// macros that depend on no other macros
#include "mympi.global.nondepmnemonics.h"

/// default settings
#include "mympi.definit.h"

/// always include these settings(must come after mympi.definit.h and mympi.global.nondepmnemonics)
#include "global.mpi_grmhd_grray_liaison.h"

/// macros, etc. that depend on other macros
#include "mympi.global.depmnemonics.h"

/// some MPI-related per-point loops
#include "mympi.global.loops.h"

/// SIMULBCCALC-related loops and related functions that doesn't yet work even in 2D, and not setup for 3D.
#include "mympi.simulbccalcstuff.h"

/// MPI-related function declarations
#include "mympi.global.funcdeclare.h"



// END
//
//////////////////////////

