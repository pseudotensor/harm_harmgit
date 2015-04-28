
/*! \file global.dump.h
  \brief Function declarations (used globally) for dump.c and dumpgen.c
 */



extern int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump,MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr),int (*content) (int i, int j, int k, MPI_Datatype datatype, void*setbuf));

extern int header1_gen(int accessmemory, int readwrite, int bintxt, int bcasthead, void *ptr, size_t size, char *format, size_t nmemb, MPI_Datatype datatype, FILE *stream);


extern int dump(long dump_cnt);
extern int dump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int dump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);

extern int avgdump(long avg_cnt);
extern int avg_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int avg2dump(long avg_cnt);
extern int avg2_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int debugdump(long debug_cnt);
extern int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int enodebugdump(long dump_cnt);
extern int enodebug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int eno_dump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);


extern int gdump(long gdump_cnt);
extern int gdump_content(int i, int j, int k, MPI_Datatype datatype, void *writebuf);

extern int fieldlinedump(long fieldline_cnt);
extern int fieldline_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int dissdump(long dump_cnt);
extern int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int dumpother(long dump_cnt);
extern int dumpother_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);


extern int fluxdumpdump(long dump_cnt);
extern int fluxdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int eosdump(long dump_cnt);
extern int eosdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int raddump(long dump_cnt);
extern int raddump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int vpotdump(long dump_cnt);
extern int vpotdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int failfloordudump(long dump_cnt);
extern int failfloordudump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int fluxsimpledump(long dump_cnt);
extern int fluxsimpledump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);


extern int dissmeasuredump(long dump_cnt);
extern int dissmeasuredump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);



extern int fakedump(long dump_cnt);
extern int fakedump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int fakedump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);


extern int image_dump(long image_cnt);
extern int imagedefs(int whichk, int scale, int limits, int vartype);
extern int image(long dump_cnt, int whichk, int scale, int limits, int vartype);
extern int image_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);
extern int image_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern void prminmaxsum(FTYPE (*p)[NSTORE2][NSTORE3][NPR], int start,int nmemb, FTYPE *max, FTYPE*min,FTYPE*sum);

extern int restart_init(int which);
extern int restart_init_simple_checks(int which);
extern int restart_init_checks(int which, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);



// restart dump
extern int restart_read(long which);
extern int check_fileformat(int readwrite, int bintxt, int whichdump, int numcolumns, int docolsplit, int mpicombine, int sizeofdatatype, FILE *stream);
extern int read_restart_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE* headerptr);
extern int restart_read_defs(void);
extern int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int restart_write(long dump_cnt);
extern int write_restart_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE* headerptr);
extern int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

// restart upperpole dump
extern int restartupperpole_read(long dump_cnt);
extern int read_restartupperpole_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);
extern int rupperpoledump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int restartupperpole_write(long dump_cnt);
extern int write_restartupperpole_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr);
extern int rupperpoledump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);




// old metric restart dump
extern int restartmetric_read(long which);
extern int read_restartmetric_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE* headerptr);
extern int rmetricdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int restartmetric_write(long dump_cnt);
extern int write_restartmetric_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE* headerptr);
extern int rmetricdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);



extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);

extern void myset(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);
extern void myget(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);

extern void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);

extern void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);
