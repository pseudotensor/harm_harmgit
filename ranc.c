//#include "global.h"
#include "decs.h"

/* ranc -- return a random deviate between 0 and 1 FTYPE ranc(iseed) * 
   iseed = integer seed for random number generator - on first call,
   will seed with iseed (if > 0) or with time(0) - will reseed anytime
   iseed > 0 */


#define RND ( 0x7fff & rand() )
#define MINN ( 2 << 8 )

// OPENMPMARK: constant static, so ok
static int P[NRANC] = {
  46337, 46327, 46309, 46307, 46301, 46279, 46273, 46271,
  46261, 46237, 46229, 46219, 46199, 46187, 46183, 46181,
  46171, 46153, 46147, 46141, 46133, 46103, 46099, 46093,
  46091, 46073, 46061, 46051, 46049, 46027, 46021, 45989,
  45979, 45971, 45959, 45953, 45949, 45943, 45893, 45887,
  45869, 45863, 45853, 45841, 45833, 45827, 45823, 45821,
  45817, 45779, 45767, 45763, 45757, 45751, 45737, 45707,
  45697, 45691, 45677, 45673, 45667, 45659, 45641, 45631
};




// ranc made thread safe by using global rancaa, rancS, and rancvaln that are shared
FTYPE ranc(int initialize, int iseed)
{
  int i;
  
  
  // thread shared random table
#pragma omp critical
  {
    
    // seed the random number generator if first time called or if iseed > 0
    if (initialize || iseed != 0) {
      if (iseed == 0) iseed = time(0);
      srand(iseed);
      rancvaln = 0;
      for (i = 0; i < NRANC; i++) {
 
        rancaa[i] = 0.0;
        rancS[i] = 0.0;
        while ((rancaa[i] = RND) < MINN || rancaa[i] > P[i] - MINN);
        while ((rancS[i] = RND) < 1 || rancaa[i] > P[i] - 1);
      }
    }
  
    rancvaln = rancS[rancvaln] & (NRANC - 1);
    rancS[rancvaln] = (rancS[rancvaln] * rancaa[rancvaln]) % P[rancvaln];
  }

  // return final random number
  return (FTYPE) rancS[rancvaln] / (FTYPE) P[rancvaln];

}

#undef NRANC
#undef RND
#undef MINN
