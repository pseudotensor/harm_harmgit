
/*! \file mpi_fprintfs.c
     \brief Various fprintf() calls that account for which core we are on
     Generally, myid=0 writes general full-core info, while all cores write specific core info
*/


#include "decs.h"

/// only to myid==0
void myfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;

  if  (myid==0) {
    va_start (arglist, format);

    if(fileptr==NULL){
      stderrfprintf("tried to print to null file pointer: %s\n",format);
      fflush(stderr);
    }
    else{
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);
    }
    va_end(arglist);
  }
}

#ifdef WIN32
#define va_copy(a,b) (a=b)
#endif

/// prints to stderr(only cpu=0) AND file pointer of choice (all cpus)
void dualfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist,arglistcopy;


  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(fileptr==NULL){
      stderrfprintf("tried to print to null file pointer: %s\n",format);
      fflush(stderr);
    }
    else{
      va_copy(arglistcopy,arglist);
      vfprintf (fileptr, format, arglistcopy);
      fflush(fileptr);
      va_end(arglistcopy);
    }
  }
  if(myid==0){
    va_copy(arglistcopy,arglist);
    vfprintf (stderr, format, arglistcopy);
    fflush(stderr);
    va_end(arglistcopy);
  }
  va_end(arglist);
}

/// prints to both logfull_file(cpu=0) and log_file(all cpus)
void logsfprintf(char *format, ...)
{
  va_list arglist,arglistcopy;


  

  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(log_file){
      va_copy(arglistcopy,arglist);
      vfprintf (log_file, format, arglistcopy);
      fflush(log_file);
      va_end(arglistcopy);
    }
  }
  if  ((myid==0)&&(logfull_file)){
    va_copy(arglistcopy,arglist);
    vfprintf (logfull_file, format, arglistcopy);
    fflush(logfull_file);
    va_end(arglistcopy);
  }
  va_end(arglist);
}



/// prints to log_file(all cpus)
void logfprintf(char *format, ...)
{
  va_list arglist,arglistcopy;


  

  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(log_file){
      va_copy(arglistcopy,arglist);
      vfprintf (log_file, format, arglistcopy);
      fflush(log_file);
      va_end(arglistcopy);
    }
  }
  va_end(arglist);
}

/// prints to logdt_file(all cpus)
void logdtfprintf(char *format, ...)
{
  va_list arglist,arglistcopy;


  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(logdt_file){
      va_copy(arglistcopy,arglist);
      vfprintf (logdt_file, format, arglistcopy);
      fflush(logdt_file);
      va_end(arglistcopy);
    }
  }
  va_end(arglist);
}

/// prints to stderr if PRODUCTION<=1 for all CPUs but only myid==0 if PRODUCTION>=2
void stderrfprintf(char *format, ...)
{
  va_list arglist,arglistcopy;



  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(stderr){
      va_copy(arglistcopy,arglist);
      vfprintf (stderr, format, arglistcopy);
      fflush(stderr);
      va_end(arglistcopy);
    }
  }
  va_end(arglist);
}

/// prints to logfull_file, log_file, and stderr (but only using cpu=0)
void trifprintf(char *format, ...)
{
  va_list arglist, arglistcopy;


  va_start (arglist, format);

  if(PRODUCTION<=2 && myid==0 || PRODUCTION<=1){
    if(log_file){
      va_copy(arglistcopy,arglist);
      vfprintf (log_file, format, arglistcopy);
      fflush(log_file);
      va_end(arglistcopy);
    }
  }
  if(myid==0){
    if(logfull_file){
      va_copy(arglistcopy,arglist);
      vfprintf (logfull_file, format, arglistcopy);
      fflush(logfull_file);
      va_end(arglistcopy);
    }
    va_copy(arglistcopy,arglist);
    vfprintf (stderr, format, arglistcopy);
    fflush(stderr);
    va_end(arglistcopy);
  }
  va_end(arglist);
}


/// cpu==0 opens a file
void myfopen(char*fname, char*fmt,char*message,FILE**fileptrptr)
{
  if(myid==0){
    *fileptrptr = fopen(fname, fmt);
    if (*fileptrptr == NULL) {
      dualfprintf(fail_file, message);
      myexit(10100);
    }
  }
}

/// cpu==0 closes a file
void myfclose(FILE ** fileptrptr,char*message)
{
  int reterror;

  if(myid==0){
    if(*fileptrptr!=NULL){
      reterror = fclose(*fileptrptr);
      if (reterror == EOF) {
        dualfprintf(fail_file, message);
        myexit(10101);
      }
    }
    else{
      dualfprintf(fail_file,"file already closed: %s\n",message);
      myexit(10102);
    }
  }
}
