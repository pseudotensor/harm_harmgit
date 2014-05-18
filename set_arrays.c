
/*! \file set_arrays.c
     \brief Wrapper for setting up (allocation/pointer shifting, dummy assignment) of all per-point, 1D, and multi-D arrays
*/


#include "decs.h"




void set_arrays()
{
  extern void set_arrays_perpoint_perline(void);
  extern void set_arrays_multidimen(void);
  extern void reconstructeno_set_arrays(void);


  set_arrays_perpoint_perline();
  set_arrays_multidimen();

  // below is here since no initialilzation routine for reconstructeno.c
  reconstructeno_set_arrays();

  // below is here because no initialization routine for reconstruct.c and wouldn't be allowed since need parallel entrance allowed.
  //  reconstructeno_set_vars();

}
