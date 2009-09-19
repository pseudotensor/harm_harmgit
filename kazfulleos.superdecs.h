// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h

// EOSextraglobal is per CPU data for EOS
#if(WHICHEOS==KAZFULL)
extern FTYPE BASEMACP0A1(EOSextraglobal,N1M,N2M,N3M,NUMEOSGLOBALS);
#endif
extern FTYPE PTRDEFGLOBALMACP0A1(EOSextraglobal,N1M,N2M,N3M,NUMEOSGLOBALS);


///////////
// global! :
// was in kazfulleos.decs.h, but there is no decs version of that created
FTYPEEOS TRUENUCLEAROFFSETPRIMARY;
