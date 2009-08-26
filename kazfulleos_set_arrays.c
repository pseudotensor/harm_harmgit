// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h


void kazfulleos_set_arrays(void)
{
  void kazfulleos_set_arrays_multidimen(void);
  void kazfulleos_set_arrays_perpoint_perline(void);


  kazfulleos_set_arrays_perpoint_perline();
  kazfulleos_set_arrays_multidimen();


}


void kazfulleos_set_arrays_perpoint_perline(void)
{

  ///////////////////////////////////////////////
  //
  // assign global pointer to eostable (static assignment)
  //
  // Take stuff in within WHICHEOS==KAZFULL in kazfulleos.eostablesdefs.h and regexp:   double BASEEOSMAC(\([a-zA-Z]+\),\([a-zA-Z0-9,]+\));    -> EOSPOINT(\1) = (double PTREOSMAC(\1,\2)) (&(BASEEOSMAC(\1,0,0,0,0,0,0,0)));
  // BUT, then must replace "0,0,0,0,0,0,0" with correct offsets.  For example, standard table should have offset at end of -FIRSTEOSSTANDARD so when accessed as standard[FIRSTEOSSTANDARD] it's really base_standard[0]
  //
  //////////////////////////////////////////////
#if(WHICHEOS==KAZFULL)
  // full
  EOSPOINT(eosfulltabledegen) = (double PTREOSMAC(eosfulltabledegen,NUMEOSDEGENQUANTITIESMEM1,EOSFULLDEGENN5,EOSFULLDEGENN4,EOSFULLDEGENN3,EOSFULLDEGENN2,EOSFULLDEGENN1,NUMEOSDEGENQUANTITIESMEM2)) (&(BASEEOSMAC(eosfulltabledegen,0,0,0,0,0,0,-FIRSTEOSDEGEN)));
  EOSPOINT(eosfulltablestandard) = (double PTREOSMAC(eosfulltablestandard,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSSTANDARDQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltablestandard,0,0,0,0,0,0,-FIRSTEOSSTANDARD)));
  EOSPOINT(eosfulltableguess) = (double PTREOSMAC(eosfulltableguess,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSGUESSQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltableguess,0,0,0,0,0,0,-FIRSTEOSGUESS)));
  EOSPOINT(eosfulltablediss) = (double PTREOSMAC(eosfulltablediss,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSDISSQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltablediss,0,0,0,0,0,0,-FIRSTEOSDISS)));
  EOSPOINT(eosfulltabledp) = (double PTREOSMAC(eosfulltabledp,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSDPQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltabledp,0,0,0,0,0,0,-FIRSTEOSDP)));
  EOSPOINT(eosfulltablesden) = (double PTREOSMAC(eosfulltablesden,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSSDENQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltablesden,0,0,0,0,0,0,-FIRSTEOSSDEN)));
  EOSPOINT(eosfulltablesspec) = (double PTREOSMAC(eosfulltablesspec,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSSSPECQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltablesspec,0,0,0,0,0,0,-FIRSTEOSSSPEC)));
  EOSPOINT(eosfulltablepofchi) = (double PTREOSMAC(eosfulltablepofchi,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSPOFCHIQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltablepofchi,0,0,0,0,0,0,-FIRSTEOSPOFCHI)));
  EOSPOINT(eosfulltabletemp) = (double PTREOSMAC(eosfulltabletemp,NUMEOSTEMPQUANTITIESMEM1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSTEMPQUANTITIESMEM2)) (&(BASEEOSMAC(eosfulltabletemp,0,0,0,0,0,0,-FIRSTEOSTEMP)));
  EOSPOINT(eosfulltableextra) = (double PTREOSMAC(eosfulltableextra,1,EOSFULL5,EOSFULL4,EOSFULL3,EOSFULL2,EOSFULL1,NUMEOSEXTRAQUANTITIESMEM)) (&(BASEEOSMAC(eosfulltableextra,0,0,0,0,0,0,-FIRSTEOSEXTRA)));

  // simple
  EOSPOINT(eossimpletabledegen) = (double PTREOSMAC(eossimpletabledegen,NUMEOSDEGENQUANTITIESMEM1,EOSSIMPLEDEGENN5,EOSSIMPLEDEGENN4,EOSSIMPLEDEGENN3,EOSSIMPLEDEGENN2,EOSSIMPLEDEGENN1,NUMEOSDEGENQUANTITIESMEM2)) (&(BASEEOSMAC(eossimpletabledegen,0,0,0,0,0,0,-FIRSTEOSDEGEN)));
  EOSPOINT(eossimpletablestandard) = (double PTREOSMAC(eossimpletablestandard,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSSTANDARDQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletablestandard,0,0,0,0,0,0,-FIRSTEOSSTANDARD)));
  EOSPOINT(eossimpletableguess) = (double PTREOSMAC(eossimpletableguess,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSGUESSQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletableguess,0,0,0,0,0,0,-FIRSTEOSGUESS)));
  EOSPOINT(eossimpletablediss) = (double PTREOSMAC(eossimpletablediss,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSDISSQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletablediss,0,0,0,0,0,0,-FIRSTEOSDISS)));
  EOSPOINT(eossimpletabledp) = (double PTREOSMAC(eossimpletabledp,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSDPQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletabledp,0,0,0,0,0,0,-FIRSTEOSDP)));
  EOSPOINT(eossimpletablesden) = (double PTREOSMAC(eossimpletablesden,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSSDENQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletablesden,0,0,0,0,0,0,-FIRSTEOSSDEN)));
  EOSPOINT(eossimpletablesspec) = (double PTREOSMAC(eossimpletablesspec,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSSSPECQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletablesspec,0,0,0,0,0,0,-FIRSTEOSSSPEC)));
  EOSPOINT(eossimpletablepofchi) = (double PTREOSMAC(eossimpletablepofchi,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSPOFCHIQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletablepofchi,0,0,0,0,0,0,-FIRSTEOSPOFCHI)));
  EOSPOINT(eossimpletabletemp) = (double PTREOSMAC(eossimpletabletemp,NUMEOSTEMPQUANTITIESMEM1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSTEMPQUANTITIESMEM2)) (&(BASEEOSMAC(eossimpletabletemp,0,0,0,0,0,0,-FIRSTEOSTEMP)));
  EOSPOINT(eossimpletableextra) = (double PTREOSMAC(eossimpletableextra,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSEXTRAQUANTITIESMEM)) (&(BASEEOSMAC(eossimpletableextra,0,0,0,0,0,0,-FIRSTEOSEXTRA)));

  // simple zoom
  EOSPOINT(eossimplezoomtabledegen) = (double PTREOSMAC(eossimplezoomtabledegen,NUMEOSDEGENQUANTITIESMEM1,EOSSIMPLEZOOMDEGENN5,EOSSIMPLEZOOMDEGENN4,EOSSIMPLEZOOMDEGENN3,EOSSIMPLEZOOMDEGENN2,EOSSIMPLEZOOMDEGENN1,NUMEOSDEGENQUANTITIESMEM2)) (&(BASEEOSMAC(eossimplezoomtabledegen,0,0,0,0,0,0,-FIRSTEOSDEGEN)));
  EOSPOINT(eossimplezoomtablestandard) = (double PTREOSMAC(eossimplezoomtablestandard,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSSTANDARDQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtablestandard,0,0,0,0,0,0,-FIRSTEOSSTANDARD)));
  EOSPOINT(eossimplezoomtableguess) = (double PTREOSMAC(eossimplezoomtableguess,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSGUESSQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtableguess,0,0,0,0,0,0,-FIRSTEOSGUESS)));
  EOSPOINT(eossimplezoomtablediss) = (double PTREOSMAC(eossimplezoomtablediss,1,EOSSIMPLE5,EOSSIMPLE4,EOSSIMPLE3,EOSSIMPLE2,EOSSIMPLE1,NUMEOSDISSQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtablediss,0,0,0,0,0,0,-FIRSTEOSDISS)));
  EOSPOINT(eossimplezoomtabledp) = (double PTREOSMAC(eossimplezoomtabledp,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSDPQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtabledp,0,0,0,0,0,0,-FIRSTEOSDP)));
  EOSPOINT(eossimplezoomtablesden) = (double PTREOSMAC(eossimplezoomtablesden,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSSDENQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtablesden,0,0,0,0,0,0,-FIRSTEOSSDEN)));
  EOSPOINT(eossimplezoomtablesspec) = (double PTREOSMAC(eossimplezoomtablesspec,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSSSPECQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtablesspec,0,0,0,0,0,0,-FIRSTEOSSSPEC)));
  EOSPOINT(eossimplezoomtablepofchi) = (double PTREOSMAC(eossimplezoomtablepofchi,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSPOFCHIQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtablepofchi,0,0,0,0,0,0,-FIRSTEOSPOFCHI)));
  EOSPOINT(eossimplezoomtabletemp) = (double PTREOSMAC(eossimplezoomtabletemp,NUMEOSTEMPQUANTITIESMEM1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSTEMPQUANTITIESMEM2)) (&(BASEEOSMAC(eossimplezoomtabletemp,0,0,0,0,0,0,-FIRSTEOSTEMP)));
  EOSPOINT(eossimplezoomtableextra) = (double PTREOSMAC(eossimplezoomtableextra,1,EOSSIMPLEZOOM5,EOSSIMPLEZOOM4,EOSSIMPLEZOOM3,EOSSIMPLEZOOM2,EOSSIMPLEZOOM1,NUMEOSEXTRAQUANTITIESMEM)) (&(BASEEOSMAC(eossimplezoomtableextra,0,0,0,0,0,0,-FIRSTEOSEXTRA)));
#endif


}




void kazfulleos_set_arrays_multidimen(void)
{

#if(WHICHEOS==KAZFULL)

  GLOBALPOINT(EOSextraglobal) = (FTYPE PTRMACP0A1(EOSextraglobal,N1M,N2M,N3M,NUMEOSGLOBALS)) (&(BASEMACP0A1(EOSextraglobal,N1BND,N2BND,N3BND,-FIRSTEOSGLOBAL))); // -FIRSTEOSGLOBAL so EOSextraglobal[FIRSTEOSGLOBAL] is a_EOSextraglobal[0]

  // indicate that IGLOBAL and other things have not been set yet.
  // Used to indicate if starting iterating on Ynu0 yet
  int i,j,k;
  COMPFULLLOOP{
    GLOBALMACP0A1(EOSextraglobal,i,j,k,IGLOBAL) = -100;
    GLOBALMACP0A1(EOSextraglobal,i,j,k,JGLOBAL) = -100;
    GLOBALMACP0A1(EOSextraglobal,i,j,k,KGLOBAL) = -100;
  }


#endif





}
