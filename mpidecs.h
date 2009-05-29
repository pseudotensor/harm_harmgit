


#ifndef USINGLIAISON
// for only grmhd code
#include "supermpidecs.h"
#endif


// for both grmhd and liaison codes whether or not DOINGGRMHDLIAISON==1
#include "mpidecs.mpi_grmhd_grray_liaison.h"




extern int romiocoliter;
extern int periodicx1, periodicx2, periodicx3;
extern int mpiperiodicx1, mpiperiodicx2, mpiperiodicx3;
extern int skipix1, reflectix1, reflectox1;
extern int skipix2, reflectix2, reflectox2;
extern int skipix3, reflectix3, reflectox3;
extern int intix1, intox1, intix2, intox2, intix3, intox3;
extern int skipintix1, skipintix2, skipintix3;
extern int ncpux1, ncpux2, ncpux3;
extern int truenumprocs;
extern int myid_world, numprocs;
extern int myid; // GRMHD CODE rank for non-MPI commands
extern int MPIid[MAXCPUS]; // GRMHD MPI internal rank for MPI commands
extern char myidtxt[MAXFILENAME];
extern int totalzones, realtotalzones,realtotalcompzones;
extern int rtotalzones;
extern int itotalzones;
//int sizes[COMPDIM + 1][MAXCPUS];
//int isizes[COMPDIM + 1][MAXCPUS];
extern int totalsize[COMPDIM + 1];
extern int itotalsize[COMPDIM + 1];
extern int mycpupos[COMPDIM + 1];		// my position amongst the cpus
extern int primgridpos[NUMBOUNDTYPES][COMPDIM*2][MAXNPR];
// NUMPACKUNPACK+1 because start at 1
extern FTYPE primfactor[NUMBOUNDTYPES][COMPDIM*2][NUMPRIMGRIDPOS][1+NUMPACKUNPACK][MAXNPR];
extern int dirgenset[NUMBOUNDTYPES][COMPDIM*2][DIRGENNUMVARS];
extern int dirloopset[NUMBOUNDTYPES][COMPDIM*2][NUMPRIMGRIDPOS][DIRLOOPNUMVARS];
extern int srdir[3*2];			// which direction this cpu
				// sends/receives normal interior data
extern int startpos[COMPDIM + 1];
extern int endpos[COMPDIM + 1];		// startj and endj are where this CPU
				// located on full grid 
extern int *startpos0[COMPDIM+1];
extern int *endpos0[COMPDIM+1];
extern int *mycpupos0[COMPDIM+1];

extern int plmpiglobal;


extern int procnamelen;
#if(USEMPI)
extern MPI_Group MPI_GROUP_WORLD;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern MPI_Status mpichstatus;
#endif

// MPI transmit vars, so minimum local code changes
extern FTYPE ndtsend, bsq_maxsend;

// for data output
extern int nextbuf,numcolumns;
extern int bufferoffset;
extern int joniosize,writebufsize;


/////////////////////////////
//
// OpenMP stuff
//
/////////////////////////////

extern int numopenmpthreads;
extern int numopenmpthreadsorig;
