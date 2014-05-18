
/*! \file mympi.global.depmnemonics.h
     \brief MPI dependent macros/definitions
*/

#if(USEMPI==0)

// still use MPI data type for communicating types to functions
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
#define MPI_BYTE           ((MPI_Datatype)3)
#define MPI_SHORT          ((MPI_Datatype)4)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_UNSIGNED       ((MPI_Datatype)7)
#define MPI_LONG           ((MPI_Datatype)8)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)13)

// fake non-MPI types to avoid USEMPI==0 conditionals on function calls
#define MPI_Request int
#define MPI_Status int

#define MPI_COMM_WORLD 0
#define MPI_COMM_GRMHD 1

#define MPI_STATUS_IGNORE 0

#define MPI_Irecv(buf,size,datatype,id,tag,comm,request)
#define MPI_Isrecv(buf,size,datatype,id,tag,comm,request)
#define MPI_Isend(buf,size,datatype,id,tag,comm,request)
#define MPI_Issend(buf,size,datatype,id,tag,comm,request)
#define MPI_Sendrecv(addrs,sizes,datatypes,ids,tags, addrr,sizer,datatyper,idr,tagr,comm,status)
#define MPI_Wait(req,status)
#define MPI_Bcast(add,size,datatype,id,comm)

#endif

// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define MPI_FTYPE MPI_FLOAT
#elif(REALTYPE==DOUBLETYPE)
#define MPI_FTYPE MPI_DOUBLE
#elif(REALTYPE==LONGDOUBLETYPE)
#define MPI_FTYPE MPI_LONG_DOUBLE
#endif


#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define MPI_SFTYPE MPI_FLOAT
#elif(SENSITIVE==DOUBLETYPE)
#define MPI_SFTYPE MPI_DOUBLE
#elif(SENSITIVE==LONGDOUBLETYPE)
#define MPI_SFTYPE MPI_LONG_DOUBLE
#endif

#if(COUNTTYPE==DOUBLETYPE)
#define MPI_CTYPE MPI_DOUBLE
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define MPI_CTYPE MPI_LONG_LONG_INT
#endif


#if(PFLAGTYPE==INTTYPE)
#define MPI_PFTYPE MPI_INT
#elif(PFLAGTYPE==CHARTYPE)
#define MPI_PFTYPE MPI_CHAR
#endif



// This forces file dumping and reading to be always in same format w.r.t. user-basis (N1,N2,N3) rather than how memory stored that could be arbitrary with ORDERSTORAGE
// This forces storage to be i fastest, j slower, and k slowest
//
#define BUFFERMAP ((long long int)bufferoffset+(long long int)(k*N1*N2+j*N1+i)*(long long int)numcolumns+(long long int)nextbuf++)
#define BUFFERMAP2 ((long long int)k*N1*N2+(long long int)j*N1+(long long int)i)
#define BUFFERINIT0 bufferoffset=0
// mpi uses BUFFERINIT in global.h as well





// MPIFLOWCONTROL setup
// for MPIFLOWCONTROL==2, setup global tag space *and* non-overlapping buffer spaces for recv's.
// So can pre-post recv's long before sends, so to avoid unexpected buffers filling up and/or direct write to application buffer.
#if(MPIFLOWCONTROL==2 || 1)
// bound_flux requires global tag space for even simple separate pre-post recv's.

// not setup because requires workbc separate for each bound call in normal computational loop

#define TAGSTARTBOUNDMPI (0) // numprocs*COMPDIM*2 in size
#define TAGSTARTBOUNDMPIINT (TAGSTARTBOUNDMPI + numprocs*COMPDIM*2) // numprocs*COMPDIM*2 in size
#define TAGSTARTBOUNDMPIPOLESMOOTH (TAGSTARTBOUNDMPIINT + numprocs*COMPDIM*2) // 2*numprocs*ncpux3 in size
#define TAGSTARTFRDOT (TAGSTARTBOUNDMPIPOLESMOOTH + 2*numprocs*ncpux3) // numprocs in size

#else

#define TAGSTARTBOUNDMPI (0)
#define TAGSTARTBOUNDMPIINT (0)
#define TAGSTARTBOUNDMPIPOLESMOOTH (0)
#define TAGSTARTFRDOT (0) // numprocs in size

#endif
