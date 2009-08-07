
////////////////
//
// globals made private with OpenMP so thread safe reentrance in global.nondepnmemonics.h in definition OPENMPKAZEOSPRIVATE
//
////////////////
int kazii,kazjj,kazkk,kazll,kazmm;
int kaziio,kazjjo,kazkko,kazllo,kazmmo;
int kazstartiii,kazstartjjj,kazstartkkk,kazstartlll,kazstartmmm;
int kazendiii,kazendjjj,kazendkkk,kazendlll,kazendmmm;
FTYPE kazdi[2],kazdj[2],kazdk[2],kazdl[2],kazdm[2];
int gottable[2]; //={0,0}; // for degen and non-degen case : should end up same!
int whichtable[2];
FTYPE indexarray[NUMINDEPDIMENS+1];
FTYPE resultold[NUMEOSQUANTITIESMEM];
int repeatedfun[NUMEOSQUANTITIESMEM];
FTYPE extrasold[MAXNUMEXTRAS];
FTYPE processedold[MAXPROCESSEDEXTRAS];
FTYPE qoldarray[NUMINDEPDIMENS+1];
FTYPE qoldarrayextras[NUMINDEPDIMENS+1];
int doallextrasold;
