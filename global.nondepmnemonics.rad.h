
/*! \file global.nondepmnemonics.rad.h
    \brief General code definitions of independent quantities for RADIATION parts of code

    // all things here don't depend on anything else, just names for numbers or purely functional macros
    // Various physics and model setup parameters that are macros either for performance reasons or since no need to change them at runtime.
*/

/// radiation equations of motion
#define EOMRADNONE 0 // no source term
#define EOMRADEDD 1 // true Edd
#define EOMRADEDDWITHFLUX 2 // fake Edd with extra flux (KORALTODO, using Fragile 2012 paper inversion -- simple to do)
#define EOMRADM1CLOSURE 3 // M1 Closure


/// KORAL
#define BOOSTGRIDPOS (NDIM) // CENT,FACE1,FACE2,FACE3 == assumes they are 0,1,2,3
#define BOOSTDIRS (2)
#define ORTHO2LAB (0) // elo
#define LAB2ORTHO (1) // eup

#define ZAMO2FF (0)
#define FF2ZAMO (1)

#define LAB2FF (0)
#define FF2LAB (1)
#define HARM2FF (2)
#define FF2HARM (3)

#define TYPEUCON (0)
#define TYPEUCOV (1)

/// array indices for prioritermethod
#define NUMPRIORITERMETHODINDEX 8
#define PRIORITERMETHODNOTSET -1
#define BASEITERMETHODINDEX 0
#define ITERMODEINDEX 1
#define IMPTRYCONVINDEX 2
#define IMPMAXITERINDEX 3
#define NUMDAMPINDEX 4
#define MODPRIMINDEX 5
#define CHECKRADINVINDEX 6
#define EOMTYPEINDEX 7



//relele.c
// same as dot(a,b) in global.variousmacros.h
//#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define kron(i,j) (i == j ? 1. : 0.)
 

#include "koral.mdefs.h"

/*********************/
//passive definitions
/*********************/



////////////////////
///
/// Mathematical or Methods constants (not physical constants)
///
////////////////////
#define Pi (M_PI)     

#define RADSHOCKFLAT 1 // 0 or 1.  Whether to include radiation in shock flatener

/// whether to fixup inversion failures using harm fixups
/// can lead to issues because diffuses, so across sharp boundary radiation can be given quite "wrong" values that don't match what solution "wants" 
#define DORADFIXUPS 0 // for RADDONUT ok, since no sharp edges.
/// KORALTODO SUPERGODMARK: Turn this on and rest all tests and see if makes worse or better.  Makes things worse at failure boundary.  Leads to very bad results for (e.g.) RADDONUT.


/// what \tau to imply mhd fluid ang rad gas are coupled enough to say single fluid as far as shock detectors or dissipation measures.
#define TAUTOTMAXSWITCH (0.5)



#define GGG0 (6.674e-8)
#define CCCTRUE0 (2.99792458e10) // cgs in cm/s
#define MSUN (1.989E33) //cgs in grams
#define ARAD0 (7.56593E-15) // cgs in erg/(K^4 cm^3)
#define ARAD (ARAD0) // only in koral would the code version depend upon G,c
#define K_BOLTZ (1.3806488e-16) // cgs in erg/K
#define M_PROTON (1.67262158e-24) // proton mass in cgs in grams
#define MB (1.66054E-24) // = 1/N_A = 1/(Avogadro's number) = baryon mass in cgs in grams (as often used in general EOSs)
#define MPOME (1836.15)
#define MELE (M_PROTON/MPOME) // electron mass in cgs in grams
//#define SIGMA_RAD (5.67e-5) // cgs in erg/(cm^2 s K^4)
#define HPLANCK (6.62607E-27) // cgs
#define QCHARGE (4.8029E-10) // cgs




////////////////
///
/// scaling for c and G constants (can be overwridden by user in their init.h)
///
////////////////
#define cTILDA (1.0) // like koral
//#define cTILDA (1E-5)
//#define gTILDA (1E-10) // like koral
#define gTILDA (1.0)


