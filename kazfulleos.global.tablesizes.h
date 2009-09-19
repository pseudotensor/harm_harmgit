/////////////////////////////
//
// Setup INPUT and HARM EOS table sizes for N1,N2,N3,N4,N5.
// Both INPUT and HARM EOS table sizes are same for these dimensions associated with independent variables
//
// EOS?N1=rho0
// EOS?N2=u or p or chi
// EOS?N3=tdynorye in seconds or Y_e
// EOS?N4=tdynorynu in seconds or Y_\nu
// EOS?N5=height in cm
//
////////////////////////////

////////////
//
// full EOS table w/ neutrino part to be computed during run-time
//
// normal and degen tables are separate files.  However, temp table is not separate file, but memory is different size because only read-in unique temperature data.  For example, for whichdatatype==4, degen table and temp table have no Ynu dependence.  So don't store redundant data.
// So have to ignore Ynu or H dependence when reading-in file for degen tables and temp-part of standard tables
//
////////////

// without extras:
#define EOSFULLN1 200
#define EOSFULLN2 200
#define EOSFULLN3 100
#define EOSFULLN4 1
#define EOSFULLN5 1         // H not tabulated

// temperature table slightly different in general:
#define EOSFULLTEMPN1 EOSFULLN1
#define EOSFULLTEMPN2 EOSFULLN2
#define EOSFULLTEMPN3 EOSFULLN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSFULLTEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSFULLTEMPN4 EOSFULLN4
#endif
#define EOSFULLTEMPN5 EOSFULLN5

#define EOSFULLDEGENN1 EOSFULLN1
#define EOSFULLDEGENN2 (1) // always 1
#define EOSFULLDEGENN3 EOSFULLN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSFULLDEGENN4 (1) // Then degen table is not function of neutrino term Y_\nu
#else
#define EOSFULLDEGENN4 EOSFULLN4
#endif
#define EOSFULLDEGENN5 EOSFULLN5 // should be 1 anyways

// extras table:
#define EOSFULLEXTRAN1 100
#define EOSFULLEXTRAN2 50
#define EOSFULLEXTRAN3 100
#define EOSFULLEXTRAN4 10
#define EOSFULLEXTRAN5 1         // H not tabulated

// temperature table slightly different in general:
#define EOSFULLEXTRATEMPN1 EOSFULLEXTRAN1
#define EOSFULLEXTRATEMPN2 EOSFULLEXTRAN2
#define EOSFULLEXTRATEMPN3 EOSFULLEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSFULLEXTRATEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSFULLEXTRATEMPN4 EOSFULLEXTRAN4
#endif
#define EOSFULLEXTRATEMPN5 EOSFULLEXTRAN5


#define EOSFULLEXTRADEGENN1 EOSFULLEXTRAN1
#define EOSFULLEXTRADEGENN2 (1) // always 1
#define EOSFULLEXTRADEGENN3 EOSFULLEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSFULLEXTRADEGENN4 (1) // Then degen table is not function of neutrino term Y_\nu
#else
#define EOSFULLEXTRADEGENN4 EOSFULLEXTRAN4
#endif
#define EOSFULLEXTRADEGENN5 EOSFULLEXTRAN5 // should be 1 anyways



/////////////////
//
// EOS table with assumptions
// If EOSSIMPLEN5==1, then assumed height such that optically thin problem
// If EOSSIMPLEN3==1, then assumed Y_e(\rhob,T) (not usual except perhaps for testing)
// degenerate table assumes to be same size except N2=1
//
/////////////////

// no extras
#define EOSSIMPLEN1 200
#define EOSSIMPLEN2 50
#define EOSSIMPLEN3 50
#define EOSSIMPLEN4 1   // here Y_\nu~0 (optically thin)
#define EOSSIMPLEN5 1   // H not tabulated

// temperature table slightly different in general:
#define EOSSIMPLETEMPN1 EOSSIMPLEN1
#define EOSSIMPLETEMPN2 EOSSIMPLEN2
#define EOSSIMPLETEMPN3 EOSSIMPLEN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLETEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSSIMPLETEMPN4 EOSSIMPLEN4
#endif
#define EOSSIMPLETEMPN5 EOSSIMPLEN5

#define EOSSIMPLEDEGENN1 EOSSIMPLEN1
#define EOSSIMPLEDEGENN2 (1) // always 1
#define EOSSIMPLEDEGENN3 EOSSIMPLEN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEDEGENN4 (1)
#else
#define EOSSIMPLEDEGENN4 EOSSIMPLEN4
#endif
#define EOSSIMPLEDEGENN5 EOSSIMPLEN5 // should be 1


// extras
#define EOSSIMPLEEXTRAN1 200
#define EOSSIMPLEEXTRAN2 50
#define EOSSIMPLEEXTRAN3 50
#define EOSSIMPLEEXTRAN4 1   // here Y_\nu~0 (optically thin)
#define EOSSIMPLEEXTRAN5 1   // H not tabulated

// temperature table slightly different in general:
#define EOSSIMPLEEXTRATEMPN1 EOSSIMPLEEXTRAN1
#define EOSSIMPLEEXTRATEMPN2 EOSSIMPLEEXTRAN2
#define EOSSIMPLEEXTRATEMPN3 EOSSIMPLEEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEEXTRATEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSSIMPLEEXTRATEMPN4 EOSSIMPLEEXTRAN4
#endif
#define EOSSIMPLEEXTRATEMPN5 EOSSIMPLEEXTRAN5

#define EOSSIMPLEEXTRADEGENN1 EOSSIMPLEEXTRAN1
#define EOSSIMPLEEXTRADEGENN2 (1) // always 1
#define EOSSIMPLEEXTRADEGENN3 EOSSIMPLEEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEEXTRADEGENN4 (1)
#else
#define EOSSIMPLEEXTRADEGENN4 EOSSIMPLEEXTRAN4
#endif
#define EOSSIMPLEEXTRADEGENN5 EOSSIMPLEEXTRAN5 // should be 1


//////////////////////
//
// EOS table with assumed Height (small) and assumed tdynorye (large)
// degenerate table assumes to be same size except N2=1
// NOT USED RIGHT NOW
//
//////////////////////
// no extras:
#define EOSSIMPLEZOOMN1 1
#define EOSSIMPLEZOOMN2 1
#define EOSSIMPLEZOOMN3 1
#define EOSSIMPLEZOOMN4 1
#define EOSSIMPLEZOOMN5 1

// temperature table slightly different in general:
#define EOSSIMPLEZOOMTEMPN1 EOSSIMPLEZOOMN1
#define EOSSIMPLEZOOMTEMPN2 EOSSIMPLEZOOMN2
#define EOSSIMPLEZOOMTEMPN3 EOSSIMPLEZOOMN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEZOOMTEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSSIMPLEZOOMTEMPN4 EOSSIMPLEZOOMN4
#endif
#define EOSSIMPLEZOOMTEMPN5 EOSSIMPLEZOOMN5

#define EOSSIMPLEZOOMDEGENN1 EOSSIMPLEZOOMN1
#define EOSSIMPLEZOOMDEGENN2 (1)
#define EOSSIMPLEZOOMDEGENN3 EOSSIMPLEZOOMN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEZOOMDEGENN4 (1)
#else
#define EOSSIMPLEZOOMDEGENN4 EOSSIMPLEZOOMN4
#endif
#define EOSSIMPLEZOOMDEGENN5 EOSSIMPLEZOOMN5

// extras:
#define EOSSIMPLEZOOMEXTRAN1 1
#define EOSSIMPLEZOOMEXTRAN2 1
#define EOSSIMPLEZOOMEXTRAN3 1
#define EOSSIMPLEZOOMEXTRAN4 1
#define EOSSIMPLEZOOMEXTRAN5 1

// temperature table slightly different in general:
#define EOSSIMPLEZOOMEXTRATEMPN1 EOSSIMPLEZOOMEXTRAN1
#define EOSSIMPLEZOOMEXTRATEMPN2 EOSSIMPLEZOOMEXTRAN2
#define EOSSIMPLEZOOMEXTRATEMPN3 EOSSIMPLEZOOMEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEZOOMEXTRATEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSSIMPLEZOOMEXTRATEMPN4 EOSSIMPLEZOOMEXTRAN4
#endif
#define EOSSIMPLEZOOMEXTRATEMPN5 EOSSIMPLEZOOMEXTRAN5

#define EOSSIMPLEZOOMEXTRADEGENN1 EOSSIMPLEZOOMEXTRAN1
#define EOSSIMPLEZOOMEXTRADEGENN2 (1)
#define EOSSIMPLEZOOMEXTRADEGENN3 EOSSIMPLEZOOMEXTRAN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEZOOMEXTRADEGENN4 (1)
#else
#define EOSSIMPLEZOOMEXTRADEGENN4 EOSSIMPLEZOOMEXTRAN4
#endif
#define EOSSIMPLEZOOMEXTRADEGENN5 EOSSIMPLEZOOMEXTRAN5








/////////////////////////////
//
// Setup input file table names
//
/////////////////////////////

// FULLTABLE:
#define EOSFULLHEADNAME "eosnew.head"
#define EOSFULLTABLENAME "eosnew.dat"
#define EOSFULLTABLEDEGENNAME "eosdegennew.dat"
// FULLTABLEEXTRA:
#define EOSFULLEXTRAHEADNAME "eosextranew.head"
#define EOSFULLTABLEEXTRANAME "eosextranew.dat"
#define EOSFULLTABLEEXTRADEGENNAME "eosextradegennew.dat"
// SIMPLETABLE:
#define EOSSIMPLEHEADNAME "eossimplenew.head"
#define EOSSIMPLETABLENAME "eossimplenew.dat"
#define EOSSIMPLETABLEDEGENNAME "eossimpledegennew.dat"
// SIMPLETABLEEXTRA:
#define EOSSIMPLEEXTRAHEADNAME "eossimpleextranew.head"
#define EOSSIMPLETABLEEXTRANAME "eossimpleextranew.dat"
#define EOSSIMPLETABLEEXTRADEGENNAME "eossimpleextradegennew.dat"
// SIMPLEZOOMTABLE:
#define EOSSIMPLEZOOMHEADNAME "eossimplezoomnew.head"
#define EOSSIMPLEZOOMTABLENAME "eossimplezoomnew.dat"
#define EOSSIMPLEZOOMTABLEDEGENNAME "eossimplezoomdegennew.dat"
// EXTRASIMPLEZOOMTABLE:
#define EOSSIMPLEZOOMEXTRAHEADNAME "eossimplezoomextranew.head"
#define EOSSIMPLEZOOMEXTRATABLENAME "eossimplezoomextranew.dat"
#define EOSSIMPLEZOOMEXTRATABLEDEGENNAME "eossimplezoomextradegennew.dat"

