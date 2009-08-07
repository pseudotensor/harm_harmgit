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
  EOSPOINT(eostable) = (double PTREOSMAC(eostable,FILL,EOSN5,EOSN4,EOSN3,EOSN2,EOSN1)) (&(BASEEOSMAC(eostable,0,0,0,0,0,0)));
  EOSPOINT(eosdegentable) = (double PTREOSMAC(eosdegentable,FILL,EOSN5,EOSN4,EOSN3,1,EOSN1)) (&(BASEEOSMAC(eosdegentable,0,0,0,0,0,0)));
  
  EOSPOINT(eossimpletable) = (double PTREOSMAC(eossimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,EOSSIMPLEN2,EOSSIMPLEN1)) (&(BASEEOSMAC(eossimpletable,0,0,0,0,0,0)));
  EOSPOINT(eosdegensimpletable) = (double PTREOSMAC(eosdegensimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,1,EOSSIMPLEN1)) (&(BASEEOSMAC(eosdegensimpletable,0,0,0,0,0,0)));

  EOSPOINT(eossimplezoomtable) = (double PTREOSMAC(eossimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,EOSSIMPLEZOOMN2,EOSSIMPLEZOOMN1)) (&(BASEEOSMAC(eossimplezoomtable,0,0,0,0,0,0)));
  EOSPOINT(eosdegensimplezoomtable) = (double PTREOSMAC(eosdegensimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,1,EOSSIMPLEZOOMN1)) (&(BASEEOSMAC(eosdegensimplezoomtable,0,0,0,0,0,0)));
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
