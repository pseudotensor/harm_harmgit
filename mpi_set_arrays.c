
/*! \file mpi_set_arrays.c
     \brief Initialize MPI arrays
*/

#include "decs.h"

// set mpi arrays (both perpoint and multidimen)
void mpi_set_arrays(void)
{
  void mpi_set_arrays_multidimen(void);
  void mpi_set_arrays_perpoint_perline(void);


  mpi_set_arrays_perpoint_perline();
  mpi_set_arrays_multidimen();


}


void mpi_set_arrays_multidimen(void)
{

  // no true multi-dimensional arrays so far

}


void mpi_set_arrays_perpoint_perline(void)
{

  // these arrays don't  have dependence on directions (i,j,k) in storage mapping function
#if(USEMPI)
  workbc = (FTYPE(*)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM]) (&(workbca[-1][0][0]));
  
  workbc_int =(PFTYPE(*)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM]) (&(workbc_inta[-1][0][0]));
#endif


}

