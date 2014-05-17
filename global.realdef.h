
/*! \file global.realdef.h
    \brief macros and definitions related to variable types

// autoparallelization
// use icc and -parallel to try autoparallelization (must  have -O2 or -O3)
// doesn't work except for a few trivial loops

// for multi-core, use icc and -openmp as compiler option for directives
// export OMP_NUM_THREADS=2
//#pragma omp parallel for private( privIndx, privDbl ) reduction( + : globalCount )
*/



/// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1
#define LONGDOUBLETYPE 2
#define LONGLONGINTTYPE 3
#define INTTYPE 4
#define CHARTYPE 5

//////////////////////////////////////////////////////
///
/// fundamental to set real type's for variables
/// define your user type here
/// (normal non-sensitive or performance critical datatypes)
#define REALTYPE DOUBLETYPE
/// (non-perf critical or sensitive data types) 
#define SENSITIVE DOUBLETYPE
/// WE ASSUME SENSITIVE>=REALTYPE !

/// counter (integer) stuff where counts can exceed integer (2 billion)
#define COUNTTYPE DOUBLETYPE // can't make long long int work, so use double
//#define COUNTTYPE LONGLONGINTTYPE // can't make long long int work, so use double

/// type for pflags
#define PFLAGTYPE CHARTYPE


/// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define FTYPE float
#elif(REALTYPE==DOUBLETYPE)
#define FTYPE double
#elif(REALTYPE==LONGDOUBLETYPE)
#define FTYPE long double
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPE double
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPE long double
#endif

