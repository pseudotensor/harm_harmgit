
/*! \file restart.h
     \brief Function declarations for restart.c
*/


// whether to assume restart file to READ is in old format or not -- assume written format is new
// note that Rebecca's code never changes number of output primitives (NPRDUMP) unlike Jon's code
// This means only have to change header information and choose to change read/write as below:
#define OLD3DMODEREAD 0
// assume want to write in new mode so have dissipation stuff in case need to restart yet again
#define OLD3DMODEWRITE 0


#if(OLD3DMODEREAD)
// old assumes DISS and DISSVSR off
#define read_restart_header read_restart_header_old
#define restart_read_defs restart_read_defs_old
#define extrarestartfunction extrarestartfunction_old
#else
#define read_restart_header read_restart_header_new
#define restart_read_defs restart_read_defs_new
#define extrarestartfunction extrarestartfunction_new
#endif

#if(OLD3DMODEWRITE)
// old assumes DISS and DISSVSR off
#define write_restart_header write_restart_header_old
#else
#define write_restart_header write_restart_header_new
#endif


int extrarestartfunction_old(void);
int read_restart_header_old(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr);
int restart_read_defs_old(void);
int write_restart_header_old(int whichdump, int whichdumpversion, int numcolumns, int bintxt,FILE*headerptr);

int extrarestartfunction_new(void);
int read_restart_header_new(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr);
int restart_read_defs_new(void);
int write_restart_header_new(int whichdump, int whichdumpversion, int numcolumns, int bintxt,FILE*headerptr);

