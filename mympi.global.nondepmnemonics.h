
/*! \file mympi.global.nondepmnemonics.h
     \brief MPI independent macros/definitions
*/

// whether to use MPI
#define USEMPI (USINGMPI) // choice (set through makehead.inc)
#define USEMPIGRMHD USEMPI // always this way
#define USEMPIGRRAY USEMPI // always this way

///////////////////
//
// MNEMOICS
//
//////////////////////

#define INITROMIO 0
#define WRITECLOSEROMIO 1
#define READROMIO 2
#define READFREEROMIO 3
#define WRITEENDROMIO 4




#define NUMPACKUNPACK 2 // starting from 1
#define PACK 1
#define UNPACK 2

#define REQRECV 0
#define REQSEND 1

#define DIRGENNUMVARS 7
#define DIRIF      0
#define DIRSIZE    1
#define DIROTHER   2
#define DIRTAGS    3
#define DIRTAGR    4
#define DIROPP     5
#define DIRNUMPR   6

#define DIRLOOPNUMVARS 18
#define DIRPSTART1 0
#define DIRPSTOP1  1
#define DIRPDIR1   2
#define DIRUSTART1 3
#define DIRUSTOP1  4
#define DIRUDIR1   5
#define DIRPSTART2 6
#define DIRPSTOP2  7
#define DIRPDIR2   8
#define DIRUSTART2 9
#define DIRUSTOP2  10
#define DIRUDIR2   11
#define DIRPSTART3 12
#define DIRPSTOP3  13
#define DIRPDIR3   14
#define DIRUSTART3 15
#define DIRUSTOP3  16
#define DIRUDIR3   17



#define TEXTOUTPUT 0
#define BINARYOUTPUT 1
#define MIXEDOUTPUT 2 // means header is text and dump is binary (handled by dump_gen()

#define UNSORTED 0
#define SORTED 1


/// simple algorithm, but eats alot of memory on cpu=0 (unbounded) if doing sorted output
#define MPICOMBINESIMPLE 0
/// more parallel:
#define MPICOMBINEMINMEM 1 // homebrew, but buggy on tungsten/mako, no problem on BH cluster -- ever.
#define MPICOMBINEROMIO 2 // requires romio package

/// for various uses (dumping and special boundary/comp interchange routine
#define STAGEM1 (-1)
#define STAGE0 0
#define STAGE1 1
#define STAGE2 2
#define STAGE3 3
#define STAGE4 4
#define STAGE5 5
#define STAGE6 6
#define STAGE7 7

/// #define DATADIR "./"
#define DATADIR ""

/// extention for data files
#define DATEXT ".dat"
#define PAREXT ".par"
#define INEXT ".in"
#define OUTEXT ".out"
#define PPEXT ".pp"

#define CPUTXT ".%04d"


#define MYOUT stderr  // normally stderr
