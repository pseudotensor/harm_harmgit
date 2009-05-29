// also includes 3 flags: 2 b^2 and 1 utoprim fail flag
#if(1)
// now always control range
// ignore num and use specified range
#define PLOOPMPI(pl,num) for(plmpiglobal=nprboundstart,pl=nprboundlist[plmpiglobal];plmpiglobal<=nprboundend;plmpiglobal++,pl=nprboundlist[plmpiglobal])
#else
#define PLOOPMPI(pr,num) for(pr=0;pr<num;pr++)
#endif
