
/*! \file utoprim_jon.h
    \brief U->P inversion definitions for utoprim_jon.c
    
*/



/// whether to turn on jon's modificatoins to Scott's general Newton's method
#define JONGENNEWTSTUFF 1

#define CRAZYNEWCHECK 0

/// 1: 1-D E'(W')  (really has 3 roots, but seemingly always positive-most root found AND correct root)
/// 2: 2-D E'(W',v^2) and P^2(W',v^2) (unsure how many roots)
/// 3: 1-D E'(W') Jon's method without intermediate v^2 and using \tilde{u}^2 instead of v^2 -- avoids catastrophic cancellation (still has 3 roots)
#define WHICHHOTINVERTER 3

/// whether to do final calculation of primitives from v^2 or \tilde{u}^2
/// 0: from \tilde{u}^2 -- avoids catastrophic cancellation at high Lorentz factors
/// 1: from v^2 -- requires computing \gamma from 1/sqrt(1-v^2) which has problems when v^2\sim 1
#define PRIMFROMVSQ 0


/// 0: original 1-D but with E' (W and v^2, but can give nan's when Wp<0 slightly -- maybe fixed with new normalization?) (apparently also not giving unique solution)
/// 1: Jon's simple 1-D E' (no use of v^2, but obtains an inconsistent v^i[P^2,E'])  (apparently also not giving unique solution)
///    Uses same as WHICHHOTINVERTER==3
/// 2: simple 1-D P^2 (no use of v^2, obtains a consistent v^i[P^2]) (but 4 roots and never sure which one is right -- not always positive root.  If D<0, then can be negative root)
/// 3: 1-D P^\alpha error on all 3 components of P^i so there is a unique answer (NOT DONE YET)
//#define WHICHCOLDINVERTER 2
#define WHICHCOLDINVERTER 2


/// 0: 1-D method using equation Sc = \rho_0 \gamma Ss = D Ss with Ss specific entropy to find W'
#define WHICHENTROPYINVERTER 0



#define METHODTYPE  3
#define NEWT_DIM 2 // actually maximum # of dimensions
#define USE_LINE_SEARCH -2

#define MAX_NEWT_RETRIES 0    /*Max. # of retries of N-R procedure, while increasing guess for W by *10 after each time*/


