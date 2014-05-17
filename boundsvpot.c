
#include "decs.h"


/*! \file boundsvpot.c
  \brief User Boundary conditions for vector potential

// For fluxes, e.g. F1, assume fluxes exist everywhere -- including j/k boundary zones.  Only i-boundary zones need to be bounded.
// This assumesCOMPZSLOOP(is,ie,js,je,ks,ke) is over boundary zones in flux.c, which in general to be compatible with any flux method (including finite volume) this is how it should be.

// With fluxes, only need to bound each dir-flux along that direction (as presently used by ENO-type schemes)

// Assume flux at 0 through N are computed correctly. So only need fluxes in other boundary zones.
// Self-assigns for 0 or N for simplicity of coding

// OUTFLOW leaves true edge of boundary unchanged
// Therefore, if FIXEDOUTFLOW, then extrapolation is always ok.
// if OUTFLOW, then extrapolation is ok as long as flux is from active zones out of boundary
*/

// order of outflow extrap
// 0: none/ copy
// 1: first order
#define EXTRAP 0 //atch



int inboundloop[NDIM];
int outboundloop[NDIM];
int innormalloop[NDIM];
int outnormalloop[NDIM];
int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
int riin,riout,rjin,rjout,rkin,rkout;
int dosetbc[COMPDIM*2];




int bound_vpot_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{

  // nothing for now.  Assume don't need vpot in ghost cells.  Or assume set true ghost cells and just getting MPI cells set.

  return (0);
}
