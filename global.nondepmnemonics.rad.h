
// radiation equations of motion
#define EOMRADNONE 0
#define EOMRADEDD 1
#define EOMRADM1CLOSURE 2

// KORAL
#define BOOSTGRIDPOS (NDIM) // CENT,FACE1,FACE2,FACE3 == assumes they are 0,1,2,3
#define BOOSTDIRS (2)
#define ORTHO2LAB (0) // elo
#define LAB2ORTHO (1) // eup

#define ZAMO2FF (0)
#define FF2ZAMO (1)

#define LAB2FF (0)
#define FF2LAB (1)

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

#define RADSOURCEMETHODEXPLICIT 1
#define RADSOURCEMETHODIMPLICIT 2
#define RADSOURCEMETHODNONE 3


///////////////////
//
// Mathematical or Methods constants (not physical constants)
//
///////////////////
#define Pi (M_PI)     

// KORALTODO: The below need to be chosen intelligently
#define RADEPS (1.e-6)
#define RADCONV (1.e-7)
#define PRADEPS (1.e-6)
#define PRADCONV (1.e-8)
#define IMPEPS (1.e-6)
#define IMPCONV (1.e-6)
#define IMPMAXITER (50)
#define GAMMASMALLLIMIT (1.0-1E-10) // at what point above which assume gamma^2=1.0




///////////////
//
// Some physical constants
//
///////////////
#define cTILDA (1.0) // like koral
//#define cTILDA (1E-5)
//#define gTILDA (1E-10) // like koral
#define gTILDA (1.0)
#define GGG0 (6.674e-8)
#define GGG (GGG0/gTILDA) // cgs in cm^3/(kg s^2)
#define CCCTRUE0 (2.99792458e10) // cgs in cm/s
#define CCCTRUE (CCCTRUE0/cTILDA) // cgs in cm/s
#define MSUN (1.989E33) //cgs in grams
#define ARAD (7.56593E-15) // cgs in erg/(K^4 cm^3)
#define K_BOLTZ (1.3806488e-16) // cgs in erg/K
#define M_PROTON (1.67262158e-24) // proton mass in cgs in grams
#define MB (1.66054E-24) // = 1/N_A = 1/(Avogadro's number) = baryon mass in cgs in grams (as often used in general EOSs)
#define SIGMA_RAD (5.67e-5) // cgs in erg/(cm^2 s K^4)

/////////////////////
//
// derived constants
//
/////////////////////
#define MSUNCM (GGG*MSUN/(CCCTRUE*CCCTRUE)) // Msun in cm
