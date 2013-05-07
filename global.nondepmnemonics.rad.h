
// radiation equations of motion
#define EOMRADNONE 0 // no source term
#define EOMRADEDD 1 // true Edd
#define EOMRADEDDWITHFLUX 2 // fake Edd with extra flux (KORALTODO, using Fragile 2012 paper inversion -- simple to do)
#define EOMRADM1CLOSURE 3 // M1 Closure

// KORAL
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

//relele.c
// same as dot(a,b) in global.variousmacros.h
//#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define kron(i,j) (i == j ? 1. : 0.)
 

#include "koral.mdefs.h"

/*********************/
//passive definitions
/*********************/


#define COURRADEXPLICIT (0.1) // Effective Courant-like factor for stiff explicit radiation source term.  Required to not only avoid failure of explicit scheme, but also that explicit scheme is really accurate compared to implicit.  E.g., near \tau\sim 1, explicit won't fail with RADPULSEPLANAR but will not give same results as implicit.  So only use explicit if really in optically thin regime.


///////////////////
//
// Mathematical or Methods constants (not physical constants)
//
///////////////////
#define Pi (M_PI)     

// KORALTODO: The below need to be chosen intelligently

// below this \tau, no source term applied.
// KORALTODO: Need to fix implicit solver so avoids dU-del in fluid if no radiatoin-fluid interaction, else overestimates effect and inversion failures occur.
#define MINTAUSOURCE (NUMEPSILON)


// IMPLICIT SOLVER TOLERANCES or DERIVATIVE SIZES
 // for used implicit solver (needs to be chosen more generally.  KORALTODO: 1E-8 too small in general).  Could start out with higher, and allow current checks to avoid inversion failure.
//#define IMPEPS (1.e-8)
// use large, and it'll go smaller if no inversion, but can't start out with too small since then Jac will have diag() terms =0
// KORALTODO: For difficult iterations, there can be solution but Jacobian is too rough and jump around alot in primitive space for small changes in U.  Should really modify IMPEPS in such cases when pr changes alot for such changes in U.
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define IMPEPS (MY1EM5)
#elif(REALTYPE==LONGDOUBLETYPE)
#define IMPEPS (MY1EM6)
#endif

// maximum EPS for getting Jacobian
#define MAXIMPEPS (0.3)

// maximum number of times to (typically) increase EPS in getting Jacobian for implicit scheme.  This might generally override MAXIMPEPS.
#define MAXJACITER (10)

#if(0)
// RADPULSEPLANAR: ~5 f1iters and ~8 iters on average
// RADTUBE NTUBE=31: ~0 f1iters and ~5 iters
// and each f1iter does 1 inversion, while each iter does 16 inversions!
// below too hard to get for more realistic problems like RADFLATDISK
#define IMPTRYCONV (1.e-12)  // for used implicit solver
#define IMPALLOWCONV (1.e-3)  // for used implicit solver
#elif(1)
// 1E-9 is common ok first iteration for RADFLATDISK.  More is too hard.
// So Choose 1E-8 as good enough solution.
#define IMPTRYCONV (1.e-8)  // for used implicit solver
#define IMPALLOWCONV (1.e-5)  // for used implicit solver KORALTODO: Have to be more careful since f/fnorm~1E-3 might mean large changes in primitives.
//#define IMPALLOWCONV (1.e-1) // KORALTODO SUPERGODMARK
#else
// RADPULSEPLANAR: below leads to ~5 f1iters and ~7 iters on average
// RADTUBE NTUBE=31: ~0 f1iters and ~1.5-2 iters
#define IMPTRYCONV (1.e-6)  // for used implicit solver
#define IMPALLOWCONV (1.e-3)  // for used implicit solver
#endif

#define IMPMAXITER (100) // for used implicit solver

// 1 : normalize radiation error by only radiation thermal energy
// 2 : normalize radiation error by max(radiation,gas) thermal energy
// 3 : normalize using radiation URAD0 but also fnorm from actual f
// 4 : URAD0, fnorm, and UU
// normalize error.  Can't expect radiation to be relatively accurate to itself if UU>>URAD0 due to G between them
// 3 is safest, but more expensive than 4.  4 should be fine for real systems.
#define IMPLICITERRORNORM 4


#define MAXSUBCYCLES (2000) // for explicit sub-cycles when doing reversion

// if tries more than this number of sub-cycles, just fail and assume no 4-force since probably due to no actual solution even for implicit scheme due to sitting at radiative failure (e.g. gamma->gammamax or Erf->ERADLIMIT)
#define MAXSUBCYCLESFAIL (MAXSUBCYCLES*100)


#define MAXF1TRIES 50 // 50 might sound like alot, but Jacobian does 4 to 5 inversions each iteration, and that amount is only typically needed for very first iteration.
  // goes to f1iter=10 for RADPULSE KAPPAES=1E3 case.  Might want to scale maximum iterations with \tau, although doubling of damping means exponential w.r.t. f1iter, so probably 50 is certainly enough since 2^(-50) is beyond machine precision for doubles.

#define RADDAMPDELTA (0.5) // choose, but best to choose 1/Integer so no machine precision issues.
#define RADDAMPUNDELTA (1.0/RADDAMPDELTA)



#define GAMMASMALLLIMIT (1.0-1E-10) // at what point above which assume gamma^2=1.0

#define RADSHOCKFLAT 1 // 0 or 1.  Whether to include radiation in shock flatener
// RADSHOCKFLAT 1 causes excessive oscillations in RADBEAMFLAT at injection point

// whether to choose Jon or Olek way of handling u2p_rad inversion failures
#define JONCHOICE 0
#define OLEKCHOICE 1

#define CASECHOICE JONCHOICE // choose
//#define CASECHOICE OLEKCHOICE // choose

#define TOZAMOFRAME 0 // reduce to ZAMO gammarel=1 frame (e.g. in non-GR that would be grid frame or v=0 frame or gammarel=1).
#define TOFLUIDFRAME 1 // reduce to using fluid frame (probably more reasonable in general).
#define TOOPACITYDEPENDENTFRAME 2

#define M1REDUCE TOOPACITYDEPENDENTFRAME // choose


// whether to fixup inversion failures using harm fixups
// can lead to issues because diffuses, so across sharp boundary radiation can be given quite "wrong" values that don't match what solution "wants" 
#define DORADFIXUPS 0 // KORALTODO SUPERGODMARK: Turn this on and rest all tests and see if makes worse or better.  Makes things worse at failure boundary.  Leads to very bad results for (e.g.) RADDONUT.

#define TAUFAILLIMIT (2.0/3.0) // at what \tau below which to assume "failure1" in u2p_rad() means should be moving at gammamax rather than not moving.

// whether to revert to sub-cycle explicit if implicit fails.  Only alternative is die.
#define IMPLICITREVERTEXPLICIT 0 // FUCK -- problem not a good idea.

// like SAFE for normal dt step, don't allow explicit substepping to change dt too fast to avoid instabilities.
#define MAXEXPLICITSUBSTEPCHANGE 1.e-2

// 0 : tau suppression
// 1 : space-time merged
// 2 : all space merged but separate from time
// 3 : full split
// 4 : split even mhd and rad
#define TAUSUPPRESS 0 // makes physical sense, but might be wrong in some limits (i.e. Gd can still be large relative), but try to account for lambda as well.
#define SPACETIMESUBSPLITNONE 1 // DON'T USE! (gets inversion failures because overly aggressive)
#define SPACETIMESUBSPLITTIME 2 // probably not ok -- need to split off mhd and rad
#define SPACETIMESUBSPLITALL 3 // probably not ok -- need to split off mhd and rad
#define SPACETIMESUBSPLITSUPERALL 4 // OK TO USE sometimes, but not always
#define SPACETIMESUBSPLITMHDRAD 5 // KINDA OK TO USE (generates noise in velocity where E is very small, but maybe ok since sub-dominant and don't care about velocity where E is small.  E evolves fine, but Rtx eventually shows issues.)
#define SPACETIMESUBSPLITTIMEMHDRAD 6 // OK TO USE sometimes (works fine and no noise in velocity because split it off.  Might have trouble in multiple dimensions if sub-dominant momentum dimension requires implicit step -- but general diffusive would suggest unlikely.  But, not efficient in optically thick regime once radiation velocity is non-trivial in magnitude)

#define WHICHSPACETIMESUBSPLIT TAUSUPPRESS // only tausuppress works in general.

//SPACETIMESUBSPLITTIMEMHDRAD // TAUSUPPRESS
// SPACETIMESUBSPLITTIMEMHDRAD  //  SPACETIMESUBSPLITMHDRAD // SPACETIMESUBSPLITSUPERALL // TAUSUPPRESS



// whether to use dUriemann and dUgeom or other dU's sitting in dUother for radiation update
#define USEDUINRADUPDATE 1






#define GGG0 (6.674e-8)
#define CCCTRUE0 (2.99792458e10) // cgs in cm/s
#define MSUN (1.989E33) //cgs in grams
#define ARAD0 (7.56593E-15) // cgs in erg/(K^4 cm^3)
#define ARAD (ARAD0) // only in koral would the code version depend upon G,c
#define K_BOLTZ (1.3806488e-16) // cgs in erg/K
#define M_PROTON (1.67262158e-24) // proton mass in cgs in grams
#define MB (1.66054E-24) // = 1/N_A = 1/(Avogadro's number) = baryon mass in cgs in grams (as often used in general EOSs)
//#define SIGMA_RAD (5.67e-5) // cgs in erg/(cm^2 s K^4)



///////////////
//
// scaling for c and G constants (can be overwridden by user in their init.h)
//
///////////////
#define cTILDA (1.0) // like koral
//#define cTILDA (1E-5)
//#define gTILDA (1E-10) // like koral
#define gTILDA (1.0)
