
/*! \file global.openmploops.h
    \brief All macros and definitions related to OpenMP loop wrappers

    // Macros and definitions related to OpenMP loops

    // Note that OPENMP3DLOOPBLOCK2IJK defined in global.storage.h because depends upon ORDERSTORAGE

    // SUPERNOTE: Assume these loops are only used for non-geometry items (i.e. not Xstore, Vstore, etc., gcov, gcon, etc.) since those only exist at one point if unless (MCOORD!=CARTMINKMETRIC)

// For general instructions about OpenMP, see:

// https://computing.llnl.gov/tutorials/openMP/
// http://perfdynamics.blogspot.com/2009/02/poor-scalability-on-multicore.html
// http://community.edc.intel.com/t5/Multicore-Virtualization-Blog/Multi-Core-Performance-and-Plumbing/bc-p/511;jsessionid=CA7E9DBCCE03E6734F5111DE91504905
// http://www.intel.com/design/pentiumii/manuals/245127.htm
// http://www.sandia.gov/news/resources/releases/2009/multicore.html
// http://www.google.com/url?q=http://www.acumem.com/images/stories/articles/HP-CAST.pdf&sa=U&start=6&ei=bM4RSoviA47msgODneX4CQ&sig2=iMevZJjH6it4ZUeqfT1jog&usg=AFQjCNENB7jrMj9bPLnAHeaezqjLPFA-qw
// http://www.google.com/url?q=http://www.ecmwf.int/newsevents/meetings/workshops/2008/high_performance_computing_13th/presentations/Yelick.pdf&sa=U&start=18&ei=QtIRStm1E47msgODneX4CQ&sig2=MdBB5BWnAwNpzRY0ix1bnQ&usg=AFQjCNGeQ7-gqqCR49giZJpnZ33UYGwe9Q
// http://www.multicoreinfo.com/category/hpc/
// http://crd.lbl.gov/%7Eoliker/papers/ipdps07.pdf
// http://www.cs.berkeley.edu/~samw/research/papers/ipdps08.pdf
// http://www.google.com/url?q=http://crd.lbl.gov/~oliker/papers/ipdps07.pdf&sa=U&start=1&ei=19MRSt2xDo7msgODneX4CQ&sig2=Iu6MRVS6j2owklJxk-rI8A&usg=AFQjCNFgeq7xM_gtde0LJYs0HGDNHpDrSQ
// http://www.hpcwire.com/offthewire/Multicore_Code_Booster.html
// http://software.intel.com/en-us/articles/performance-tools-for-software-developers-auto-parallelization-and-qpar-threshold/
// http://software.intel.com/en-us/articles/performance-tools-for-software-developers-auto-parallelization-and-qpar-threshold/
// http://www.cs.ucsb.edu/~tyang/class/pthreads/index_sgi.html
// https://computing.llnl.gov/tutorials/openMP/exercise.html
// http://msdn.microsoft.com/en-us/magazine/cc163717.aspx
// https://computing.llnl.gov/tutorials/openMP/#ThreadprivateExamples
// http://sam.zoy.org/writings/programming/gprof.html
// http://books.google.com/books?id=18CmnqIhbhUC&pg=PA59&lpg=PA59&dq=openmp+reduction+operator+on+element+of+structure&source=bl&ots=sUk-5CRZhV&sig=WFGKz6IdgUDI6dNKG8HdGW8ACls&hl=en
// http://www.cita.utoronto.ca/~merz/pi/

*/



/// below should be larger than # of threads or cores used!
/// SUPERGODMARK: OPENMPOPTMARK:
#define OPENMPNUMCHUNKS 100 // when large blocksize, if too small chunks, then overhead is large
#define MINCHUNKSIZE 10 // if too few chunks, then overhead is large (GODMARK: should worry about if <10 iterations? -- not for now)
/// OPENMPNOTE: If really blocksize<OPENMPMINCHUNKNUMBER, then other threads just stall as sufficient
/// Note that by using not too many chunks, each thread is well-spaced in memory, so avoids false sharing problem.


/// below seems best on average
#if(N1*N2*N3<200) // overhead becomes a problem for too few iterations in loops
#define OPENMPSCHEDULE(arg) static
#define OPENMPCHUNKSIZE(blocksize) (MAX(blocksize/numopenmpthreads,MINCHUNKSIZE))
#else
#define OPENMPSCHEDULE(arg) guided
#define OPENMPCHUNKSIZE(blocksize) (MAX(blocksize/OPENMPNUMCHUNKS,MINCHUNKSIZE))
#endif

/// below allows compiler and/or run-time system to decide, which may be best.
///#define OPENMPSCHEDULE(arg) auto // option doesn't seem to exist in icc.

/// below is for loops ensured not to vary in how long each iteration takes
/// Used to avoid overhead from guided that is not needed for simple loops (e.g. simple = just setting to 0 or just taking simple difference of variables)
#define OPENMPNOVARYSCHEDULE(arg) static

/// below used when don't want to provide CHUNK argument 
#define OPENMPFULLNOVARYSCHEDULE(arg) static

/// below is for loops with very different times for each iteration (e.g. inversion loop)
#define OPENMPVARYENDTIMESCHEDULE(arg) guided

/// Oddly, performance changes by 10% when changing how start and stop blockijk (and fixing how used too of course)
//#define OPENMP3DLOOPBLOCK     for(blockijk=0;blockijk<blocksize;blockijk++)
#define OPENMP3DLOOPBLOCK     for(blockijk=1;blockijk<=blocksize;blockijk++)



#define OPENMP3DLOOPVARSDEFINE int nxsize, nxshift, nysize, nyshift, nzsize, nzshift, blocksize, blockijk


/// COMP version of SUPERGENLOOP for OpenMP : Not used anywhere currently
#define OPENMP3DLOOPSETUPSUPERGENCOMP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) \
  {nxsize=(  ((istop+SHIFTX1UP) - (istart+SHIFTX1DN))*di +1 );          \
    nxshift=istart+SHIFTX1DN;                                           \
    nysize=( ((jstop+SHIFTX2UP) - (jstart+SHIFTX2DN))*dj +1 );          \
    nyshift=jstart+SHIFTX2DN;                                           \
    nzsize=( ((kstop+SHIFTX3UP) - (kstart+SHIFTX3DN))*dk +1 );          \
    nzshift=kstart+SHIFTX3DN;                                           \
    blocksize=nxsize*nysize*nzsize;}

/// non-comp loop (or used in interlpline.c where start and stops already have grid section shifts
#define OPENMP3DLOOPSETUPSUPERGEN(istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) \
  {nxsize=(  (istop - istart)*di +1 );                                  \
    nxshift=istart;                                                     \
    nysize=( (jstop - jstart)*dj +1 );                                  \
    nyshift=jstart;                                                     \
    nzsize=( (kstop - kstart)*dk +1 );                                  \
    nzshift=kstart;                                                     \
    blocksize=nxsize*nysize*nzsize;}

/// forced to be a computational loop, so start/stop cannot already have shifts in them (which is normal for everywhere except interpline.c)
#define OPENMP3DLOOPSETUP(istart,istop,jstart,jstop,kstart,kstop)       \
  {nxsize=((istop+SHIFTX1UP) - (istart+SHIFTX1DN) +1 );                 \
    nxshift=istart+SHIFTX1DN;                                           \
    nysize=((jstop+SHIFTX2UP) - (jstart+SHIFTX2DN) +1 );                \
    nyshift=jstart+SHIFTX2DN;                                           \
    nzsize=((kstop+SHIFTX3UP) - (kstart+SHIFTX3DN) +1 );                \
    nzshift=kstart+SHIFTX3DN;                                           \
    blocksize=nxsize*nysize*nzsize;}


/// COMPFULLLOOP equivalent for OpenMP
#define OPENMP3DLOOPSETUPFULL OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND,-N3BND,N3-1+N3BND)

#define OPENMP3DLOOPSETUPFULLINOUT2 OPENMP3DLOOPSETUP(-N1BND+2*SHIFT1,N1-1+N1BND-2*SHIFT1,-N2BND+2*SHIFT2,N2-1+N2BND-2*SHIFT2,-N3BND+2*SHIFT3,N3-1+N3BND-2*SHIFT3)

#define OPENMP3DLOOPSETUPFULLINOUT2DIR1 OPENMP3DLOOPSETUP(-N1BND+2*SHIFT1,N1-1+N1BND-2*SHIFT1,-N2BND,N2-1+N2BND,-N3BND,N3-1+N3BND)

#define OPENMP3DLOOPSETUPFULLINOUT2DIR2 OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,-N2BND+2*SHIFT2,N2-1+N2BND-2*SHIFT2,-N3BND,N3-1+N3BND)

#define OPENMP3DLOOPSETUPFULLINOUT2DIR3 OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND,-N3BND+2*SHIFT3,N3-1+N3BND-2*SHIFT3)


/// COMPFULLLOOPP1 equivalent for OpenMP
#define OPENMP3DLOOPSETUPFULLP1 OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND+SHIFT1,-N2BND,N2-1+N2BND+SHIFT2,-N3BND,N3-1+N3BND+SHIFT3)

#define OPENMP3DLOOPSETUPFULLP2 OPENMP3DLOOPSETUP(-N1BND-SHIFT1,N1-1+N1BND+SHIFT1*2,-N2BND-SHIFT2,N2-1+N2BND+SHIFT2*2,-N3BND-SHIFT3,N3-1+N3BND+SHIFT3*2)

#define OPENMP3DLOOPSETUPFULLP1EXCEPTX2 OPENMP3DLOOPSETUP(-N1BND,N1-1+N1BND+SHIFT1,-N2BND,N2-1+N2BND,-N3BND,N3-1+N3BND+SHIFT3)

/// COMPZLOOP equivalent for OpenMP
#define OPENMP3DLOOPSETUPZLOOP OPENMP3DLOOPSETUP(0,N1-1,0,N2-1,0,N3-1)



/// This is used for BC's that have SHIFT's included or not automatically, so shouldn't appear here.
#define OPENMPBCLOOPSETUP(istart,istop,jstart,jstop,kstart,kstop)       \
  {nxsize=((istop) - (istart) +1 );                                     \
    nxshift=istart;                                                     \
    nysize=((jstop) - (jstart) +1 );                                    \
    nyshift=jstart;                                                     \
    nzsize=((kstop) - (kstart) +1 );                                    \
    nzshift=kstart;                                                     \
    blocksize=nxsize*nysize*nzsize;}


#define OPENMPBCLOOPVARSDEFINELOOPX1DIR int nxsize, nxshift, nysize, nyshift, nzsize, nzshift, blocksize, blockijk, fooi

#define OPENMPBCLOOPVARSDEFINELOOPX2DIR int nxsize, nxshift, nysize, nyshift, nzsize, nzshift, blocksize, blockijk, fooj

#define OPENMPBCLOOPVARSDEFINELOOPX3DIR int nxsize, nxshift, nysize, nyshift, nzsize, nzshift, blocksize, blockijk, fook

/// LOOPX1dir equivalent (assumes use foo variable for i iterator to not conflict with true i iterator)
/// Note like in global.loops.boundaries.h, assume X1-direction is done first since uses innormal and outnormal first
#define OPENMPBCLOOPSETUPLOOPX1DIR OPENMPBCLOOPSETUP(0,0,innormalloop[2],outnormalloop[2],innormalloop[3],outnormalloop[3])

/// LOOPX2dir equivalent (assumes use foo variable for j iterator to not conflict with true j iterator)
/// Note like in global.loops.boundaries.h, assume X2-direction is done second
#define OPENMPBCLOOPSETUPLOOPX2DIR OPENMPBCLOOPSETUP(inboundloop[1],outboundloop[1],0,0,innormalloop[3],outnormalloop[3])

/// LOOPX3dir equivalent (assumes use foo variable for k iterator to not conflict with true k iterator)
/// Note like in global.loops.boundaries.h, assume X3-direction is done last
#define OPENMPBCLOOPSETUPLOOPX3DIR OPENMPBCLOOPSETUP(inboundloop[1],outboundloop[1],inboundloop[2],outboundloop[2],0,0)

#define OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k) OPENMP3DLOOPBLOCK2IJK(fooi,j,k);

#define OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k) OPENMP3DLOOPBLOCK2IJK(i,fooj,k);

#define OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j) OPENMP3DLOOPBLOCK2IJK(i,j,fook);

/// basic loop block is the same
#define OPENMPBCLOOPBLOCK OPENMP3DLOOPBLOCK








