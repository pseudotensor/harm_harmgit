

// TODO:
// 1) Is my interpolation best/fastest/correct?  Thompson et al. (2003) use bivariate interpolation that looks alot more complicated
// 2) New H calculation (MAKE MPI)
// 3) check additional code that's presently in debug mode
// 4) check MPI stuff here and for gravity

// whether to allow Kaz EOS table
#define ALLOWKAZEOS 0 // expensive for OpenMP due to many large globals, so normally disable unless required

//////////////////////////////
//
// Some often changable variables
//
///////////////////////////////



// SUPERNOTE: If turn on off ALLOW????TABLE, then probably need to change kazfulleos.c:which_eostable() since particular to types of tables and not general!

// SUPERNOTE: Also should ensure EOS???N? are correct in kazfulleos.global.tablesizes.h


// whether to allow use of full table (if 0, then others must be turned on)
#define ALLOWFULLTABLE 1

// whether to only use full table (0) for allow use of simple tables if can (1)
#define ALLOWSIMPLETABLE 1

// whether to use simplezoom table if can
// zoom not needed anymore with new degen offset method
// indeed, set to 0 unless make new simplezoom table
#define ALLOWSIMPLEZOOMTABLE 0

// how many dimensions to consider.  Other dimensions' values will be consider as the dimension's lowest value
#define WHICHEOSDIMEN 4

// define which Ynu to evolve and keep as primitive (i.e. pr[YNU])
// OLD COMMENT: pr[YNU] = Ynu[orig] = Ynu (i.e. not Ynu0 so that conservation law equation for Y_\nu is correctly using radiative transfer version of Y_\nu, while EOSextra[YNU] is Ynu0 used for table lookup
// note that if ynu changes meaning from Ynu to Ynu0 depending upon WHICHEVOLVEYNU==EVOLVEYNUNOTRAD or WHICHEVOLVEYNU==EVOLVEYNURAD
#define EVOLVEYNUNOTRAD 0
#define EVOLVEYNURAD 1
#define WHICHEVOLVEYNU EVOLVEYNUNOTRAD




// whether to use degen offset (otherwise assume degen offset from file is 0 even if read-in differently)
#define ALLOWDEGENOFFSET 1


// whether to check if table returns a valid EOS value by using existence of stored inversion to temperature
// so far only setup for F(rho0,u)
// if some invalid, then don't use those data points.  If all surrounding points are invalid, then use them as if they were valid
#define CHECKIFVALIDEOSDATA 1


// whether to truncate \rho_0 and u,p,\chi independent variables in table lookup when they are beyond the table.
// Note that we always truncate in Ye and Ynu0
#define TRUNCATEHIGHRHOU 1


// SUPERTODO: need extrapolation for pipelined lookup!
// whether to extrapolate at rho or u,p,chi,s beyond table (overrides truncation)
#define EXTRAPOLATEHIGHRHOU 1

// using log interpolation results in much smoother results, and consistent with eos_extract.m for interpolation
// That is, using integer position is log-interp since all independents are log on the grid
// And those functions in eos_extract.m interpolated as log are here interpolated as log
// 0 or 1
#define DOLOGINTERP 1
// whether to prelogify table values that should be log interpolated
#define DOPRELOGIFY DOLOGINTERP
#define OUTOFBOUNDSPRELOGIFY BIG // 10^(BIG) is used to trigger that value is out of bounds and was unable to be logified (i.e. <=0.0 before log).  Impossible to have had value logify to BIG, so works.  Note that this is in code units.


// which EOS to reduce to if beyond table
// mignone doesn't make sense
#define REDUCE2WHICHEOS IDEALGAS
// ensure that gamideal is chosen
// GODMARK: Could choose nearest tabulated value of dp/du|rho0 and dp/dchi|rho0 for gamideal when indeps are rho0,u and rho0,chi

// pick size of data type for EOS stuff
#define REALTYPEEOS DOUBLETYPE

#if(REALTYPEEOS==DOUBLETYPE)
#define FTYPEEOS double
#define MPI_FTYPEEOS MPI_DOUBLE
#define EOSHEADERONEIN "%lf"
#elif(REALTYPEEOS==FLOATTYPE)
#define FTYPEEOS float
#define MPI_FTYPEEOS MPI_FLOAT
#define EOSHEADERONEIN "%f"
#elif(REALTYPEEOS==LONGDOUBLETYPE)
#define FTYPEEOS long double
#define MPI_FTYPEEOS MPI_LONG_DOUBLE
#define EOSHEADERONEIN "%Lf"
#endif


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

#define YELOOPTYPELOG 0
#define YELOOPTYPELINEAR 1
#define YELOOPTYPESPECIAL 2

// number of Ynu0 iterations per Ynu0 computation
// need to do at least more than 1 -- 2 minimum at t=0 to avoid derivative guess error
#define NUMBEROFYNU0ITERATIONS 1




///////////////////////////
//
// Some constants, tolerances, etc.
//
///////////////////////////

// tolerance to check whether repeated case for i,j,k,rho0,u
#define OLDTOLERANCE (1E-14)

// tolerance for checks on input values of table
#define TABLETOL (1E-14)

// tolerance for ye index check
#define TABLETOLYEINDEX (1E-13)

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




// e.g. superdefs.h like: FTYPEEOS BASEEOSMAC(name,....)
//      superdefs.pointers.h like: FTYPEEOS PTRDEFEOSMAC(name,....)
// set_arrays_multidimen.c like: EOSPOINT(name) = (FTYPEEOS PTREOSMAC(name,..)) (&(BASEEOSMAC(name,...)));
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
