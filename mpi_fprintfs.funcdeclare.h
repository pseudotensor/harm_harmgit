
/*! \file mpi_fprintfs.funcdeclare.h
     \brief Function declarations for global use of functions in mpi_fprintfs.c
*/




// mpi_fprintfs.c:
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void logfprintf(char *format, ...);
extern void logdtfprintf(char *format, ...);
extern void stderrfprintf(char *format, ...);
extern void trifprintf(char *format, ...);
extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);
