

////////////////////////////////////
//
// Setup EOSextra[]
//
// Some global position variables used to determine EOS
//
/////////////////////////////////////

// should be 4
#define NUMHDIRECTIONS 4

#if(NUMHDIRECTIONS!=4)
#error "NUMHDIRECTIONS should be 4"
#endif



////////////////////////////////////
//
// Setup macros to variables in EOSextra[]
//
/////////////////////////////////////

// these should be ordered and numbered such that correspond to EOS table independent variables
// do NOT correspond to expanded independent variables list from EOS as read-in (i.e. not RHOEOS, UEOS, PEOS, CHIEOS, SEOS,  YEEOS, YNUEOS,  HEOS)
#define RHOGLOBAL (-2) // dummy reference
#define UGLOBAL (-1) // dummy reference
#define TDYNORYEGLOBAL (0)        // Tdyn or Y_e depending upon whichrnpmethod
#define YNU0GLOBAL     (TDYNORYEGLOBAL+1) // Tdyn or Y_\nu depending upon whichynumethod

#define YNU0OLDGLOBAL  (YNU0GLOBAL+1)     // Tdyn or Y_\nu depending upon whichynumethod
#define YNUOLDGLOBAL   (YNU0OLDGLOBAL+1)  // Tdyn or Y_\nu depending upon whichynumethod (is actually Ynu latest for WHICHEVOLVEYNU==EVOLVEYNUNOTRAD)

// scale-height (for method that uses this for EOS, some averaged version of H)
// 2,3,4 are other directions for axisymmetric emission
#define HGLOBAL  (YNUOLDGLOBAL+1)
#define H2GLOBAL (HGLOBAL+1)
#define H3GLOBAL (H2GLOBAL+1) 
#define H4GLOBAL (H3GLOBAL+1)

// SUPERTODO: Need 4 more H's for inelastic scattering terms.

#define UNUGLOBAL (H4GLOBAL+1)       // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define PNUGLOBAL (UNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define SNUGLOBAL (PNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define PGASGLOBAL (SNUGLOBAL+1) // Pgas -- used because want to use consistent lookup of pressure instead of using p(rho0,u) mixed with p(rho0,chi).  Since inversion forces p(rho0,chi), generally stay consistent with chi version by storing it and using it later for phys.tools.c
// below 3 used to indicate position if EOS coming from grid
#define IGLOBAL (PGASGLOBAL+1)
#define JGLOBAL (IGLOBAL+1)
#define KGLOBAL (JGLOBAL+1)

// Ye should always be first
// Note that RHO,U,etc. GLOBAL should resolve to 1,2,3,4,5
// below is whatever comes after RHO and U
#define FIRSTEOSGLOBAL (TDYNORYEGLOBAL)
#define LASTEOSGLOBAL (KGLOBAL)

// NOTE: must be in same order and number as EOS independent vars
// number of per CPU position-based data for EOS
#define NUMEOSGLOBALS (LASTEOSGLOBAL-FIRSTEOSGLOBAL+1)
#define MAXPARLIST (NUMEOSGLOBALS) // for get_EOS_parms()
