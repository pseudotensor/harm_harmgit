
////////////////
//
// globals made private with OpenMP so thread safe reentrance in global.nondepnmemonics.h in definition OPENMPKAZEOSPRIVATE
//
////////////////
extern int kazii,kazjj,kazkk,kazll,kazmm;
extern int kaziio,kazjjo,kazkko,kazllo,kazmmo;
extern int kazstartiii,kazstartjjj,kazstartkkk,kazstartlll,kazstartmmm;
extern int kazendiii,kazendjjj,kazendkkk,kazendlll,kazendmmm;
extern FTYPE kazdi[2],kazdj[2],kazdk[2],kazdl[2],kazdm[2];
extern int gottable[2]; //={0,0}; // for degen and non-degen case : should end up same!
extern int whichtable[2];
extern FTYPE indexarray[NUMINDEPDIMENS+1];
extern FTYPE qoldarray[NUMINDEPDIMENS+1];
extern FTYPE resultold[NUMEOSQUANTITIESMEM];
extern int repeatedfun[NUMEOSQUANTITIESMEM];
extern FTYPE extrasold[MAXNUMEXTRAS];
extern FTYPE processedold[MAXPROCESSEDEXTRAS];
extern FTYPE qoldarray[NUMINDEPDIMENS+1];
