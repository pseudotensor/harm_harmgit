#include "decs.h"


///////////////////////////
// General method outline:
// 1) Downsample full MPI grid to all CPUs so every CPU has same downsampled dataset
// 2) Each CPU ray traces an (e.g.) grid-based fraction of rays so have a ray-trace from every point.
// 3) Compute H at every downsampled point on the grid-based fraction of points
// 4) Compute extras and use H+extras->processed at every grid-based fraction of points
// 5) Exhange full or just needed processed quantities to all CPUs.  If used grid-based fractions to start ray-tracings, then already each CPU has what it needs except for boundary issues
// 6) Interpolate P's to get P at required true grid point for the given CPU.

// This setup will be useful in connecting with Avery's code for the ray-tracing aspect.


////////////////////////////////
// Geodesic integration method:
// Charles says he uses velocity Verlet (see grmonty paper) and reduces to standard RK4 if no convergence.
// They also tried GSL's RK4:  http://www.gnu.org/software/gsl/
// Apparently standard RK4 is fast and accurate.  Want error control, so use NR's error control version of RK4.
// Use cubic interpolation of quantities on grid (e.g. \Gamma^\mu_{\nu\kappa}, primitives, etc.).

#define NUMINTERACTIONS 2 // for now, just Rosseland-Mean for total (elastic+inelastic) and inelastic processes
//#define NRAYTRACE1 50 
//#define NRAYTRACE2 50
#define NRAYTRACE1 1
#define NRAYTRACE2 1
#define NRAYTRACE3 1

FTYPE opacity[NUMINTERACTIONS][NRAYTRACE1][NRAYTRACE2][NRAYTRACE3];

// inputs:
// x0[NDIM]: emission position in internal coordinates
// u0[NDIM]: emission 4-velocity in internal coordinates
// opacity[NUMINTERACTIONS][][][]: opacities for each types of interaction

// use generic RK4 since apparently most robust according to grmonty paper.
// for now, use linear interpolation on grid quantities -- probably require cubic in general to match to 2nd order so RK4 stepping is resolved.  For now, hope that RK4 can handle jumpy steps.


// keep track of:
// 1) \tau
// 2) energy and number deposition for dU/d\tau and dYe/dt
// 3) 
int trace_geodesic(FTYPE *x0, FTYPE *u0)
{

  // \nu = - C k^\mu u_\mu (GRMONTY20)
  // d\tau_a = D (\nu \alpha) d\lambda (GRMONTY18)
  // I_\nu/\nu^3 \propto \omega with d\omega/d\tau_a = -\omega
  //
  // So propagate \omega_{n+1,prop} = \omega_{n,prop} \exp(-\tau_a)
  // Desosit: \omega_{n,dep} = \omega_{n,prop} (1-\exp(-\tau_a))
  // where \tau_a = D 1/2( (\nu\alpha)_n + (\nu\alpha)_{n+1} ) \delta\lambda

  // also track cumulative optical depth for pressure, etc.
  // \tau_{cum,a} = \int (d\tau/dL)_{co} \nu d\lambda

  // 

  return(0);

}






int downsample_harmquantities(void)
{

  return(0);




}
