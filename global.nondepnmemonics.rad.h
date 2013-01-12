
// KORAL
#define BOOSTGRIDPOS (NDIM) // CENT,FACE1,FACE2,FACE3 == assumes they are 0,1,2,3
#define BOOSTDIRS (2)
#define ZAMO2LAB (0) // elo
#define LAB2ZAMO (1) // eup

#define ZAMO2FF (0)
#define FF2ZAMO (1)



#define my_max(x,y) (x>y?x:y)
#define ldouble FTYPE
#define gSIZE 20 //size of metric arrays = 16 + 1 (gdet) + 3 (dlgdet)
//relele.c
// same as dot(a,b) in global.variousmacros.h
//#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define kron(i,j) (i == j ? 1. : 0.)
#define SOURCETERMS
#define NUM_SOURCESTEPS 1
  

#include "koral.mdefs.h"

/*********************/
//passive definitions
/*********************/


#ifndef IMPLABPREC
#define IMPLABPREC 1.e-7 //precision for the numerical solver in solve_implicit_lab()
#endif

#ifndef EXPLICIT_RAD_SOURCE
#ifndef IMPLICIT_FF_RAD_SOURCE
#ifndef IMPLICIT_LAB_RAD_SOURCE
#define IMPLICIT_LAB_RAD_SOURCE
#endif
#endif
#endif

#define RADSOURCEMETHODEXPLICIT 1
#define RADSOURCEMETHODIMPLICIT 2
