
/*! \file freespace.c
  \brief Get free space on disk
 */

#include "decs.h"

#ifndef WIN32

#if( 0 == OSX )
// below for free HD space check in dump.c
#include <sys/vfs.h>
//#include <sys/statfs.h>
#else
#include <sys/param.h>
#include <sys/mount.h>
#endif

/// Check if enough free space on disk
/// this doesn't compile on some icc versions (8.0 seems fine), works fine in gcc
int isenoughfreespace(unsigned long long need)
{
  // assume all cpus check this and get reasonbly same answer

  struct statfs buf;
  unsigned long  long x;
  int result;
  long waititer;
  int errorsend;
  char mysys[MAXFILENAME];
  unsigned long long minspace=1024*1024*100; // min 100MB

  x=0;
  waititer=0;
  //result= 0: good 1: bad (not enough space)
  while((result=((2*x)<need+minspace))){  // factor of 2 for safety

    if(waititer>0){
      
      if (numprocs > 1) {
        errorsend = result;
#if(USEMPI)
        MPI_Allreduce(&errorsend, &result, 1, MPI_INT, MPI_MAX,MPI_COMM_GRMHD);
#endif
      }
      if(result>0){
        dualfprintf(fail_file,"Ran out of hard drive space -- waiting 60 seconds , waititer=%ld\n",waititer);
        if(myid==0 && !MPIAVOIDFORK){
          if(MAILFROMREMOTE){
            sprintf(mysys,"resulthdspace=`echo \"\\`pwd\\` ran out of HD space\"` ; ssh %s@%s \"echo $resulthdspace | mail %s\"",REMOTEUSER,REMOTEHOST,EMAILADDRESS);
            stderrfprintf("system call: %s\n",mysys);
            system(mysys);
          }
          else{
            sprintf(mysys,"echo \"`pwd` ran out of HD space | mail %s\"",EMAILADDRESS);
            stderrfprintf("system call: %s\n",mysys);
            system(mysys);
          }
        }
        // 2nd time here, so no space, so sleep before checking space again
        sleep(60); // wait a minute
        // man 3 sleep
      }
    }
    // find free space on local disk
    statfs(".",&buf);
    x=(unsigned long long) buf.f_bfree*(unsigned long long)buf.f_bsize;
    //  printf("%lld\n",x);
    
    waititer++; // and then try again
  }
  return(0);

}


#else
int isenoughfreespace(unsigned long long need)
{
  return(0);
}

#endif

