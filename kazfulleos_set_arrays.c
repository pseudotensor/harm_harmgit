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
  //////////////////////////////////////////////
#if(WHICHEOS==KAZFULL)
  EOSPOINT(eostable) = (double PTREOSMAC(a_eostable,FILL,EOSN5,EOSN4,EOSN3,EOSN2,EOSN1)) (&(BASEEOSMAC(a_eostable,0,0,0,0,0,0)));
  EOSPOINT(eosdegentable) = (double PTREOSMAC(a_eosdegentable,FILL,EOSN5,EOSN4,EOSN3,1,EOSN1)) (&(BASEEOSMAC(a_eosdegentable,0,0,0,0,0,0)));
  
  EOSPOINT(eossimpletable) = (double PTREOSMAC(a_eossimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,EOSSIMPLEN2,EOSSIMPLEN1)) (&(BASEEOSMAC(a_eossimpletable,0,0,0,0,0,0)));
  EOSPOINT(eosdegensimpletable) = (double PTREOSMAC(a_eosdegensimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,1,EOSSIMPLEN1)) (&(BASEEOSMAC(a_eosdegensimpletable,0,0,0,0,0,0)));

  EOSPOINT(eossimplezoomtable) = (double PTREOSMAC(a_eossimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,EOSSIMPLEZOOMN2,EOSSIMPLEZOOMN1)) (&(BASEEOSMAC(a_eossimplezoomtable,0,0,0,0,0,0)));
  EOSPOINT(eosdegensimplezoomtable) = (double PTREOSMAC(a_eosdegensimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,1,EOSSIMPLEZOOMN1)) (&(BASEEOSMAC(a_eosdegensimplezoomtable,0,0,0,0,0,0)));
#endif


}




void kazfulleos_set_arrays_multidimen(void)
{

#if(WHICHEOS==KAZFULL)

  GLOBALPOINT(EOSextraglobal) = (FTYPE PTRMACP0A1(EOSextraglobal,N1M,N2M,N3M,NUMEOSGLOBALS)) (&(BASEMACP0A1(EOSextraglobal,N1BND,N2BND,N3BND,-FIRSTEOSGLOBAL))); // -FIRSTEOSGLOBAL so EOSextraglobal[FIRSTEOSGLOBAL] is a_EOSextraglobal[0]

#endif





}
