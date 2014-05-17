
/*! \file global.openmpthreadprivates.h
    \brief All macros and definitions related to OpenMP private quantities

// Those variables in defs.general.h (or any global variable) that per-point or per-time may change, so required to make private for OpenMP
//
// OPENMPNOTE: If global variable that should be private makes no sense to be private, then must redo code!
//
// Notes:
// 1) *All* variables in OPENMPGLOBALPRIVATEFULL are threadprivate, so each thread has own private global copy, which is not initialized unless use copyin().
//     Hence, if make something threadprivate then MUST copyin() if that is used *at all* regardless of if it was changed within parallel region.
// only need to copyin() things that *change*, rest are *shared* and set outside parallel regions
//
// 2) All threadprivate variables should NOT be *set* to something *outside* a parallel region.  Should always set them inside parallel region otherwise tid!=0 threads will not be set since outside parallel region only sets tid==0 thread.

*/


//#define OPENMPGLOBALPRIVATEPLOOPFULL nprstart,nprend,nprlist
#define OPENMPGLOBALPRIVATEPLOOP2INTERPFULL npr2interpstart,npr2interpend,npr2interplist
#define OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL npr2notinterpstart,npr2notinterpend,npr2notinterplist
//#define OPENMPGLOBALPRIVATEPLOOPBOUNDFULL nprboundstart,nprboundend,nprboundlist
//#define OPENMPGLOBALPRIVATEPLOOPFLUXBOUNDFULL nprfluxboundstart,nprfluxboundend,nprfluxboundlist
//#define OPENMPGLOBALPRIVATEPLOOPDUMPFULL nprdumpstart,nprdumpend,nprdumplist
//#define OPENMPGLOBALPRIVATEPLOOPINVERTFULL nprinvertstart,nprinvertend,nprinvertlist


/// whocalleducon only used for diagnostics, depracated check_pr() code, or diag_flux_general() that is not called in parallel region
/// Removed global use of icurr,jcurr,kcurr,pcurr
//#define OPENMPGLOBALPRIVATEOTHER2 icurr,jcurr,kcurr,pcurr,whocalleducon
//#define OPENMPGLOBALPRIVATEOTHER2 icurr,jcurr,kcurr,pcurr
//#define OPENMPGLOBALPRIVATEOTHER3 ifail,jfail,kfail // these don't change in parallel regions

#if(WHICHVEL==VEL3 && USEOPENMP==1)
#define OPENMPGLOBALPRIVATEOTHER uttdiscr // ignored for now
#error Setup openmpthreadprivates for the above quantity, neglected for now in the below copyins
#endif




/// OPENMPNOTE: Should figure out way to make the below local to eos functions rahter than global, but this is ok for now
/// OPENMPNOTE: Redid how set EOS functions (oustide parallel regions) so no longer need the below as threadprivate
/// Below were included in "state" and "full" threadprivate before moved outside parallel regions
#define INDEXPARAMETERSNAMES kaziiwhichd,kazjjwhichd,kazkkwhichd,kazllwhichd,kazmmwhichd, \
    kaziiowhichd,kazjjowhichd,kazkkowhichd,kazllowhichd,kazmmowhichd,   \
    kazstartiiiwhichd,kazstartjjjwhichd,kazstartkkkwhichd,kazstartlllwhichd,kazstartmmmwhichd, \
    kazendiiiwhichd,kazendjjjwhichd,kazendkkkwhichd,kazendlllwhichd,kazendmmmwhichd, \
    kazdiwhichd,kazdjwhichd,kazdkwhichd,kazdlwhichd,kazdmwhichd

/// below defined in kazfulleos.c [enabled or disabled below using ALLOWKAZEOS]
//#define OPENMPKAZEOSPRIVATE INDEXPARAMETERSNAMES,indexarray,qoldarray,whichtable,resultold,repeatedfun,qoldarrayextras,extrasold,processedold,doallextrasold
#define OPENMPKAZEOSPRIVATE indexarray,qoldarray,whichtable,resultold,repeatedfun,qoldarrayextras,extrasold,processedold,doallextrasold



////////////
///
/// Below are some macros for various parts of the code [no new actual variables listed on the below lines]
///
////////////
//#define OPENMPGLOBALPRIVATEFORGEOM
//#define OPENMPGLOBALPRIVATEFORUCONANDGEOM

#if(ALLOWKAZEOS)
#define OPENMPGLOBALPRIVATELIST OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL,OPENMPKAZEOSPRIVATE
#define OPENMPGLOBALPRIVATEFULL copyin(OPENMPGLOBALPRIVATELIST)
/// geom and state (state is most expensive due to EOS stuff):
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOM copyin(OPENMPKAZEOSPRIVATE)
/// geom and ucon (but not state (i.e not EOS)):
/// geom and state but with ploop, interp, and no interp varaibles only
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL,OPENMPKAZEOSPRIVATE)
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERPFULLNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL,OPENMPKAZEOSPRIVATE)
/// need "state" stuff (i.e. EOS stuff) since inversion processes EOS
/// Need EOS for inversion itself!
#define OPENMPGLOBALPRIVATEFORINVERSION copyin(OPENMPKAZEOSPRIVATE)
/// geom but not state:
#define OPENMPGLOBALPRIVATEFORGEOMNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)
#define OPENMPGLOBALPRIVATEFORUCONANDGEOMNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)
#define OPENMPGLOBALPRIVATEPLOOPINTERPONLY copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)

#else
#define OPENMPGLOBALPRIVATELIST OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL
#define OPENMPGLOBALPRIVATEFULL copyin(OPENMPGLOBALPRIVATELIST)
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOM 
/// need "state" stuff (i.e. EOS stuff) since inversion processes EOS
/// Need EOS for inversion itself!
#define OPENMPGLOBALPRIVATEFORINVERSION
/// geom and state but with ploop, interp, and no interp varaibles only
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL)
#define OPENMPGLOBALPRIVATEFORSTATEANDGEOMINTERPFULLNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL,OPENMPGLOBALPRIVATEPLOOP2NOTINTERPFULL)
/// geom but not state:
#define OPENMPGLOBALPRIVATEFORGEOMNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)
#define OPENMPGLOBALPRIVATEFORUCONANDGEOMNPR2INTERP copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)
#define OPENMPGLOBALPRIVATEPLOOPINTERPONLY copyin(OPENMPGLOBALPRIVATEPLOOP2INTERPFULL)
#endif



