
////////////////
//
// globals made private with OpenMP so thread safe reentrance in global.nondepmnemonics.h in definition OPENMPKAZEOSPRIVATE
//
// variables need to keep so can quickly get EOS values if repeated input q1-q5
//
////////////////

// [NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1] used because each "whichd" can have independent repeated lookup
//int kaziiwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazjjwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazkkwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazllwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazmmwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1];
//int kaziiowhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazjjowhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazkkowhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazllowhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazmmowhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1];
//int kazstartiiiwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazstartjjjwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazstartkkkwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazstartlllwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazstartmmmwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1];
//int kazendiiiwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazendjjjwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazendkkkwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazendlllwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1],kazendmmmwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1];
//FTYPEEOS kazdiwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2],kazdjwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2],kazdkwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2],kazdlwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2],kazdmwhichd[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2];

// each "whichd" has its own indexarray and qoldarray so independently repeatable
FTYPEEOS indexarray[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][NUMINDEPDIMENSMEM]; // [UTOTDIFF, etc.][q1-q5]
FTYPE qoldarray[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][NUMINDEPDIMENSMEM];
// whichtable: which table qarray led to (is ==NOTABLE for not within any table)
int whichtable[NUMEXTRATABLETYPES][NUMEOSDEGENQUANTITIESMEM1][2];

//////////
//
// Old function values : must access pipeline using firsteosintablesubtype[whichtablesubtype]+coli and not just coli
//
//////////
FTYPE resultold[MAXEOSPIPELINE];
// repeatedfun: whether this function has a previously stored value.  Since subtables are linked to whichd, then indexarray linked to different sub tables.
int repeatedfun[MAXEOSPIPELINE];

///////////////
// old EOS table position and results for extras and processed
// only whichd==UTOTDIFF required
FTYPE qoldarrayextras[NUMINDEPDIMENSMEM];
FTYPE extrasold[MAXNUMEXTRAS];
FTYPE processedold[MAXPROCESSEDEXTRAS];
int doallextrasold;
