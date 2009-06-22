// still use MPI data type for communicating types to functions
#if(USEMPI==0)
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
