////////////////////////////////////
//
// Allocate space for table
//
///////////////////////////////////////////////
//
// This entire file assumes all things are FTYPE except eostable itself that is explicitly double
//
///////////////////////////////////////////////

// presumed to be several functions as functions of 4 other quantities that are the independent variables (UEOS,PEOS,CHIEOS share same index)
// Notice that eostable has density as EOSN1 so fastest indicies are related to density and pressure-related quantities since this will result in fastest lookup for fixed H and Tdynorye
// due to float limit of 1E38, need doubles
// the 1 on N2 is where the DEGEN table is stored for this EOS table




// OPENMPMARK: These can be static since pointers and globally shared values of EOS table
#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eostable,NUMEOSQUANTITIESMEM,EOSN5,EOSN4,EOSN3,EOSN2,EOSN1);
#endif
static double PTRDEFEOSMAC(eostable,FILL,EOSN5,EOSN4,EOSN3,EOSN2,EOSN1);

#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eosdegentable,NUMEOSDEGENQUANTITIESMEM,EOSN5,EOSN4,EOSN3,1,EOSN1);
#endif
static double PTRDEFEOSMAC(eosdegentable,FILL,EOSN5,EOSN4,EOSN3,1,EOSN1);


// simple density-internal energy dependent EOS table
#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eossimpletable,NUMEOSQUANTITIESMEM,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,EOSSIMPLEN2,EOSSIMPLEN1);
#endif
static double PTRDEFEOSMAC(eossimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,EOSSIMPLEN2,EOSSIMPLEN1);

#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eosdegensimpletable,NUMEOSDEGENQUANTITIESMEM,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,1,EOSSIMPLEN1);
#endif
static double PTRDEFEOSMAC(eosdegensimpletable,FILL,EOSSIMPLEN5,EOSSIMPLEN4,EOSSIMPLEN3,1,EOSSIMPLEN1);


// simplezoom density-internal energy dependent EOS table
#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eossimplezoomtable,NUMEOSQUANTITIESMEM,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,EOSSIMPLEZOOMN2,EOSSIMPLEZOOMN1);
#endif
static double PTRDEFEOSMAC(eossimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,EOSSIMPLEZOOMN2,EOSSIMPLEZOOMN1);

#if(WHICHEOS==KAZFULL)
double BASEEOSMAC(eosdegensimplezoomtable,NUMEOSDEGENQUANTITIESMEM,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,1,EOSSIMPLEZOOMN1);
#endif
static double PTRDEFEOSMAC(eosdegensimplezoomtable,FILL,EOSSIMPLEZOOMN5,EOSSIMPLEZOOMN4,EOSSIMPLEZOOMN3,1,EOSSIMPLEZOOMN1);


