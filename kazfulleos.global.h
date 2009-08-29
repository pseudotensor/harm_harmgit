

// TODO:
//
// 1) choose case where checking if within EOS is not done (assume always there and truncate instead of extend)
//    TWO if's removed then: if(iswithin_eostable) and in each type of calculation the call: if(get_eos_fromtable()) are removed : might speed things up
// 2) Is my interpolation best/fastest/correct?  Thompson et al. (2003) use bivariate interpolation that looks alot more complicated


// TODO NEW:
// 1) New H calculation (MAKE MPI)
// 3) check additional code that's presently in debug mode
// 4) check MPI stuff here and for gravity

// whether to allow Kaz EOS table
#define ALLOWKAZEOS 0 // expensive for OpenMP due to many large globals, so normally disable unless required

//////////////////////////////
//
// Some often changable variables
//
///////////////////////////////

// whether to allow use of full table (if 0, then others must be turned on)
#define ALLOWFULLTABLE 1

// whether to only use full table (0) for allow use of simple tables if can (1)
#define ALLOWSIMPLETABLE 0

// whether to use simplezoom table if can
// zoom not needed anymore with new degen offset method
// indeed, set to 0 unless make new simplezoom table
#define ALLOWSIMPLEZOOMTABLE 0

// how many dimensions to consider.  Other dimensions' values will be consider as the dimension's lowest value
#define WHICHEOSDIMEN 4


// whether to use degen offset (otherwise assume degen offset from file is 0 even if read-in differently)
#define ALLOWDEGENOFFSET 1


// whether to check if table returns a valid EOS value by using existence of stored inversion to temperature
// so far only setup for F(rho0,u)
// if some invalid, then don't use those data points.  If all surrounding points are invalid, then use them as if they were valid
#define CHECKIFVALIDEOSDATA 1


// using log interpolation results in much smoother results, and consistent with eos_extract.m for interpolation
// That is, using integer position is log-interp since all independents are log on the grid
// And those functions in eos_extract.m interpolated as log are here interpolated as log
// 0 or 1
#define DOLOGINTERP 1


// which EOS to reduce to if beyond table
// mignone doesn't make sense
#define REDUCE2WHICHEOS IDEALGAS
// ensure that gamideal is chosen
// GODMARK: Could choose nearest tabulated value of dp/du|rho0 and dp/dchi|rho0 for gamideal when indeps are rho0,u and rho0,chi



#include "kazfulleos.global.tablecolumnsizes.h"

#include "kazfulleos.global.tablesizes.h"

#include "kazfulleos.global.eosextra.h"


// GODMARK: could have a table for Ynu=thermalized and have an array that stores when source term forces Ynu to be perfectly thermal, and use that table in that case.
// generating it now




///////////////
//
// Table limits indices

#define UPDOWN 2 // 0=down 1=up

#define TBLITEMS (UPDOWN+2+2)
#define TBLLINEARITEMS (UPDOWN+2)


// for calling usereduced_eos()
#define REDUCENOOFFSET 0
#define REDUCEUSEOFFSET 1






///////////////////////////
//
// Some constants, tolerances, etc.
//
///////////////////////////

// tolerance to check whether repeated case for i,j,k,rho0,u
#define OLDTOLERANCE (1E-14)

// tolerance for checks on input values of table
#define TABLETOL (1E-14)

// tolerance for truncation error-level checks
// beyond 30% accuracy start reporting issues
#define TABLETOLTRUNCATION (0.3)

// value of read-in temperature such that below this is treated as indicating an invalid (rho0,u) EOS pair
// actual read-in value is 1E-20, but using 5E-20 guarantees no machine-error choices and works with floats too
// also, generally is more accurate as temperature since problems with inversion are near T~0
#define INVALIDTEMP (5E-20)



// initialize kaziio, etc. with this so first call has no old index used
#define INITKAZINDEX (-100)



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// define EOS array macros to avoid some confusion with spatial array macros
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define GENEOSPOINT(prefix,name) prefix##name
#define GENEOSTABLEMAC(prefix,name,a1,a2,a3,a4,a5,a6,a7) prefix##name[a1][a2][a3][a4][a5][a6][a7]
//
#define EOSBASEPOINT(name) GENEOSPOINT(a_,name) // not used
#define BASEEOSMAC(name,a1,a2,a3,a4,a5,a6,a7) GENEOSTABLEMAC(a_,name,a1,a2,a3,a4,a5,a6,a7)
//
#define EOSPOINT(name) GENEOSPOINT(,name)
#define EOSMAC(name,a1,a2,a3,a4,a5,a6,a7) GENEOSTABLEMAC(,name,a1,a2,a3,a4,a5,a6,a7)
#define PTRDEFEOSMAC(name,a1,a2,a3,a4,a5,a6,a7) (*EOSPOINT(name))[a2][a3][a4][a5][a6][a7]
#define PTREOSMAC(name,a1,a2,a3,a4,a5,a6,a7) (*)[a2][a3][a4][a5][a6][a7]




// e.g. superdefs.h like: double BASEEOSMAC(name,....)
//      superdefs.pointers.h like: double PTRDEFEOSMAC(name,....)
// set_arrays_multidimen.c like: EOSPOINT(name) = (double PTREOSMAC(name,..)) (&(BASEEOSMAC(name,...)));
//
// \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> BASEEOSMAC(\1,\2,\3,\4,\5,\6,\7) [kazfulleos.c at top and kazfulleos_set_arrays.c for most-RHS of pointer shifting code]
//  OR -> EOSMAC(\1,\2,\3,\4,\5,\6,\7) [kazfulleos.c in code]
//
// (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\),
//     -> (*)PTR\6(\7,FILL,\1,\2,\3,\4,\5)) (&(\6(\7,   [in kazfulleos_set_arrays.c]
// Then :
//     -> (\*)PTRBASEEOSMAC -> PTREOSMAC    [in kazfulleos_set_arrays.c]
// Then:
// Replace (e.g.) eostable with EOSPOINT
//
// 
// ( *\* *\([_a-zA-Z0-9+-]+\) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *;
//     -> PTRDEFEOSMAC(\1,FILL,\2,\3,\4,\5,\6); [defining pointer in kazfulleos.c just after BASEEOSMAC defines global array]
//
