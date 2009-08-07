#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

#include "SwapEndian.h"



#ifndef WIN32
#include<sys/stat.h>
#endif

///////////////////
//
// pre-global stuff:
//
///////////////////

// whether doing performance testing (also see pnmhd code)
#define DOINGLIAISON 0 // choice of 0 or 1 (should always be 1 for liaison mode and liaison code comilation)
#define NCSA 0
#define PERFTEST 0
#include "mytime.h"

#define USEOPENMP (USINGOPENMP) // choice (set through makehead.inc)


#include "kazfulleos.global.h" // for array sizes

//////////////////
//
// GLOBAL STUFF:
//
//////////////////


#ifndef GLOBAL_H
#define GLOBAL_H

// define FTYPE, SFTYPE, etc.
#include "global.realdef.h"


// all pure nmenomics that don't depend upon anything else but what's inside file
// Doesn't even depend upon N?, N?M, or N?BND, etc.
#include "global.nondepnmemonics.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Default global and user choices for various code options that can change for each run without significant modifications
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// default values
// Must not make storage items
#include "definit.h"

// user specific values
// Must not make storage items
#include "init.h"

// OpenMP "threadprivate" related macros
#include "global.openmpthreadprivates.h"

// define nmenomics that depend on other mnemonics
// Must not make storage items
#include "global.depnmemonics.h"

// Determines how we store arrays and how we loop
// Depends upon N?,N?BND,SHIFT? created in global.depnmemonics.h
#include "global.storage.h"

// UP TO THIS POINT SHOULD NOT CREATE MEMORY ITEMS that depend upon grid sizes, since depend upon global.storage.h setting the storage sizes

// loops use ORDERSTORAGE set in definit.h or init.h by user
#include "global.loops.h"

#include "global.variousmacros.h"

#include "global.fieldmacros.h"

#include "global.structs.h"


// now that all hashes have been defined, get mpi header
#include "mympi.h"
#include "global.gridsectioning.h"
#include "global.comploops.h"
#include "global.openmploops.h"

// all global external function declarations
#include "global.funcdeclare.h"

#include "global.dump.h" // func declarations

#include "global.bounds.h" // func declarations

// some inits function declarations
#include "global.inits.h"

#include "global.other.h"

// put in the below file things not to be converted by double2longdouble.sh
#include "unconverted_by_double2longdouble.inc"

#endif// endif for #ifndef GLOBAL_H



/////////////////////////////////////////
//
// test declarations
// see checkexterns.sh
// comment below 3 lines when not wanting to test
// uncomment when wanting to test, and then look at make.log
//
//#if(!defined(doublereal))
//typedef double doublereal;
//#endif
//#include "utoprim_jon.h" // extra defs
//#include "temptempfinalallc.txt"
//
//
/////////////////////////////////////////





