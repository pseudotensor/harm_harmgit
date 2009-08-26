
////////////////
//
// globals made private with OpenMP so thread safe reentrance in global.nondepnmemonics.h in definition OPENMPKAZEOSPRIVATE
//
// variables need to keep so can quickly get EOS values if repeated input q1-q5
//
////////////////

// [NUMEOSDEGENQUANTITIESMEM1] used because each "whichd" can have independent repeated lookup
extern int kaziiwhichd[NUMEOSDEGENQUANTITIESMEM1],kazjjwhichd[NUMEOSDEGENQUANTITIESMEM1],kazkkwhichd[NUMEOSDEGENQUANTITIESMEM1],kazllwhichd[NUMEOSDEGENQUANTITIESMEM1],kazmmwhichd[NUMEOSDEGENQUANTITIESMEM1];
extern int kaziiowhichd[NUMEOSDEGENQUANTITIESMEM1],kazjjowhichd[NUMEOSDEGENQUANTITIESMEM1],kazkkowhichd[NUMEOSDEGENQUANTITIESMEM1],kazllowhichd[NUMEOSDEGENQUANTITIESMEM1],kazmmowhichd[NUMEOSDEGENQUANTITIESMEM1];
extern int kazstartiiiwhichd[NUMEOSDEGENQUANTITIESMEM1],kazstartjjjwhichd[NUMEOSDEGENQUANTITIESMEM1],kazstartkkkwhichd[NUMEOSDEGENQUANTITIESMEM1],kazstartlllwhichd[NUMEOSDEGENQUANTITIESMEM1],kazstartmmmwhichd[NUMEOSDEGENQUANTITIESMEM1];
extern int kazendiiiwhichd[NUMEOSDEGENQUANTITIESMEM1],kazendjjjwhichd[NUMEOSDEGENQUANTITIESMEM1],kazendkkkwhichd[NUMEOSDEGENQUANTITIESMEM1],kazendlllwhichd[NUMEOSDEGENQUANTITIESMEM1],kazendmmmwhichd[NUMEOSDEGENQUANTITIESMEM1];
extern FTYPE kazdiwhichd[NUMEOSDEGENQUANTITIESMEM1][2],kazdjwhichd[NUMEOSDEGENQUANTITIESMEM1][2],kazdkwhichd[NUMEOSDEGENQUANTITIESMEM1][2],kazdlwhichd[NUMEOSDEGENQUANTITIESMEM1][2],kazdmwhichd[NUMEOSDEGENQUANTITIESMEM1][2];

// each "whichd" has its own indexarray and qoldarray so independently repeatable
extern FTYPE indexarray[NUMEOSDEGENQUANTITIESMEM1][NUMINDEPDIMENS+1]; // [UTOTDIFF, etc.][q1-q5]
extern FTYPE qoldarray[NUMEOSDEGENQUANTITIESMEM1][NUMINDEPDIMENS+1];
// whichtable: which table qarray led to (is ==NOTABLE for not within any table)
extern int whichtable[NUMEOSDEGENQUANTITIESMEM1][2];

//////////
// Old function values
extern FTYPE resultold[NUMTABLESUBTYPES][MAXEOSPIPELINE];
// repeatedfun: whether this function has a previously stored value.  Since subtables are linked to whichd, then indexarray linked to different sub tables.
extern int repeatedfun[NUMTABLESUBTYPES][MAXEOSPIPELINE];

///////////////
// old EOS table position and results for extras and processed
extern FTYPE qoldarrayextras[NUMINDEPDIMENS+1];
extern FTYPE extrasold[MAXNUMEXTRAS];
extern FTYPE processedold[MAXPROCESSEDEXTRAS];
extern int doallextrasold;
