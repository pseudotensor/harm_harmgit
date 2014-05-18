
/*! \file mympi.global.loops.h
     \brief MPI loop definitions/macros
*/


// also includes 3 flags: 2 b^2 and 1 utoprim fail flag

#define PLOOPMPIORIG(pr,num) for(pr=0;pr<num;pr++)

#if(1)
// now always control range
// ignore num and use specified range
#define PLOOPMPI(pl,num) for(plmpiglobal=nprboundstart,pl=nprboundlist[plmpiglobal];plmpiglobal<=nprboundend;plmpiglobal++,pl=nprboundlist[plmpiglobal])
#else
#define PLOOPMPI(pr,num) PLOOPMPIORIG(pr,num)
#endif
