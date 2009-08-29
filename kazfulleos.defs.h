
////////////////
//
// globals (OPENMPMARK: that are fixed in value for all threads)
//
// SUPERNOTE: Must create corresponding MPI_Bcast() for these in kazfulleos_read_setup_eostable.c
//
////////////////
static int numeosquantitiesinstandardlist[NUMTBLS]; // standard number of of canonical list of columns corresponding to "in" format.  Should be same for standard and extra tables, for example.
static int numeosdegenquantitiesinstandardlist[NUMTBLS]; // (for degen table) standard number of of canonical list of columns corresponding to "in" format.  Should be same for standard and extra tables, for example.
static int numeosquantitiesinfile[NUMTBLS]; // stored number of quantities per independent variable in read-in table
static int numeosdegenquantitiesinfile[NUMTBLS]; // stored number of quantities per independent variable in read-in degen table
static int numallquantitiesinfile[NUMTBLS]; // total number of quantities per independent variable in read-in table
static int numalldegenquantitiesinfile[NUMTBLS]; // total number of quantities per independent variable in read-in degen table
static int extralimits[NUMTBLS][2]; // limits for extras

// first [4] : above 4 types of independent variables
// second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = log10 value of offset for log_base stepping so can control how resolved
static FTYPE inputtablelimits[NUMTBLS][NUMEOSINDEPS][TBLITEMS];
// below is same as input, but converted for easy use
static FTYPE tablelimits[NUMTBLS][NUMEOSINDEPS][TBLITEMS];
// first [4] as above, second [2] : 0=lower linear limit, 1=upper linear limit, 2 = log(-base) 3 = linear offset
static FTYPE lineartablelimits[NUMTBLS][NUMEOSINDEPS][TBLLINEARITEMS]; // makes no sense for temperature-like quantities
// first [4] : as first [4] above
static int tablesize[NUMTBLS][NUMEOSINDEPS];

static int vartypeeosextraarray[NUMINDEPDIMENSMEM];
static int vartypearray[NUMINDEPDIMENSMEM]; // unchanging

static int vartypeheightarray[NUMHDIRECTIONS+1];

static int vardegentypearray[NUMEOSDEGENINDEPS+1];
static int varnormalcompare2degentypearray[NUMEOSDEGENINDEPS+1]; // to be used to compare against degen table


static int numcolintablesubtype[NUMTABLESUBTYPES];
static int whichdintablesubtype[NUMTABLESUBTYPES];
static int firsteosintablesubtype[NUMTABLESUBTYPES];

static int whichtablesubtypeinquantity[NUMEOSQUANTITIESMEM];
static int whichcolinquantity[NUMEOSQUANTITIESMEM];

static int dologinterp_sub_coli[NUMTABLESUBTYPES][MAXEOSPIPELINE];

// code value of invalid temperature and log10 version
static FTYPE invalidtempcode,invalidlogtempcode;

static int whichrnpmethod[NUMTBLS], whichynumethod[NUMTBLS], whichhcmmethod[NUMTBLS];
static int whichdatatype[NUMTBLS],utotdegencut[NUMTBLS],numc[NUMTBLS],numextras[NUMTBLS];

static int primarytable=NOTABLE; // NOTABLE indicates no EOS setup/read-in yet



static FTYPE FAKE2IDEALNUCLEAROFFSET;
static FTYPE TRUENUCLEAROFFSET, DEGENNUCLEAROFFSET;
static FTYPE lsoffset, fakelsoffset;

static int didsetupkazeos; // OPENMPMARK: changes once, but master thread only calls setup







