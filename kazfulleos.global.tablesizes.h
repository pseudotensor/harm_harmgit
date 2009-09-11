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
#define EOSEXTRAFULLN1 100
#define EOSEXTRAFULLN2 50
#define EOSEXTRAFULLN3 100
#define EOSEXTRAFULLN4 10
#define EOSEXTRAFULLN5 1         // H not tabulated

// temperature table slightly different in general:
#define EOSEXTRAFULLTEMPN1 EOSEXTRAFULLN1
#define EOSEXTRAFULLTEMPN2 EOSEXTRAFULLN2
#define EOSEXTRAFULLTEMPN3 EOSEXTRAFULLN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRAFULLTEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSEXTRAFULLTEMPN4 EOSEXTRAFULLN4
#endif
#define EOSEXTRAFULLTEMPN5 EOSEXTRAFULLN5


#define EOSEXTRAFULLDEGENN1 EOSEXTRAFULLN1
#define EOSEXTRAFULLDEGENN2 (1) // always 1
#define EOSEXTRAFULLDEGENN3 EOSEXTRAFULLN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRAFULLDEGENN4 (1) // Then degen table is not function of neutrino term Y_\nu
#else
#define EOSEXTRAFULLDEGENN4 EOSEXTRAFULLN4
#endif
#define EOSEXTRAFULLDEGENN5 EOSEXTRAFULLN5 // should be 1 anyways



/////////////////
//
// EOS table with assumptions
// If EOSSIMPLEN5==1, then assumed height such that optically thin problem
// If EOSSIMPLEN3==1, then assumed Y_e(\rhob,T) (not usual except perhaps for testing)
// degenerate table assumes to be same size except N2=1
//
/////////////////

// no extras
#define EOSSIMPLEN1 1
#define EOSSIMPLEN2 1
#define EOSSIMPLEN3 1
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
#define EOSEXTRASIMPLEN1 1
#define EOSEXTRASIMPLEN2 1
#define EOSEXTRASIMPLEN3 1
#define EOSEXTRASIMPLEN4 1   // here Y_\nu~0 (optically thin)
#define EOSEXTRASIMPLEN5 1   // H not tabulated

// temperature table slightly different in general:
#define EOSEXTRASIMPLETEMPN1 EOSEXTRASIMPLEN1
#define EOSEXTRASIMPLETEMPN2 EOSEXTRASIMPLEN2
#define EOSEXTRASIMPLETEMPN3 EOSEXTRASIMPLEN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRASIMPLETEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSEXTRASIMPLETEMPN4 EOSEXTRASIMPLEN4
#endif
#define EOSEXTRASIMPLETEMPN5 EOSEXTRASIMPLEN5

#define EOSEXTRASIMPLEDEGENN1 EOSEXTRASIMPLEN1
#define EOSEXTRASIMPLEDEGENN2 (1) // always 1
#define EOSEXTRASIMPLEDEGENN3 EOSEXTRASIMPLEN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRASIMPLEDEGENN4 (1)
#else
#define EOSEXTRASIMPLEDEGENN4 EOSEXTRASIMPLEN4
#endif
#define EOSEXTRASIMPLEDEGENN5 EOSEXTRASIMPLEN5 // should be 1


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
#define EOSEXTRASIMPLEZOOMN1 1
#define EOSEXTRASIMPLEZOOMN2 1
#define EOSEXTRASIMPLEZOOMN3 1
#define EOSEXTRASIMPLEZOOMN4 1
#define EOSEXTRASIMPLEZOOMN5 1

// temperature table slightly different in general:
#define EOSEXTRASIMPLEZOOMTEMPN1 EOSEXTRASIMPLEZOOMN1
#define EOSEXTRASIMPLEZOOMTEMPN2 EOSEXTRASIMPLEZOOMN2
#define EOSEXTRASIMPLEZOOMTEMPN3 EOSEXTRASIMPLEZOOMN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRASIMPLEZOOMTEMPN4 (1) // Then temperature not dependent upon neutrino Y_\nu since table function of du,dp,dchi,ds instead of u,p,chi,s
#else
#define EOSEXTRASIMPLEZOOMTEMPN4 EOSEXTRASIMPLEZOOMN4
#endif
#define EOSEXTRASIMPLEZOOMTEMPN5 EOSEXTRASIMPLEZOOMN5

#define EOSEXTRASIMPLEZOOMDEGENN1 EOSEXTRASIMPLEZOOMN1
#define EOSEXTRASIMPLEZOOMDEGENN2 (1)
#define EOSEXTRASIMPLEZOOMDEGENN3 EOSEXTRASIMPLEZOOMN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSEXTRASIMPLEZOOMDEGENN4 (1)
#else
#define EOSEXTRASIMPLEZOOMDEGENN4 EOSEXTRASIMPLEZOOMN4
#endif
#define EOSEXTRASIMPLEZOOMDEGENN5 EOSEXTRASIMPLEZOOMN5








/////////////////////////////
//
// Setup input file table names
//
/////////////////////////////

// FULLTABLE:
#define EOSFULLHEADNAME "eosnew.head"
#define EOSFULLTABLENAME "eosnew.dat"
#define EOSFULLTABLEDEGENNAME "eosdegennew.dat"
// EXTRAFULLTABLE:
#define EOSEXTRAFULLHEADNAME "eosextranew.head"
#define EOSEXTRAFULLTABLENAME "eosextranew.dat"
#define EOSEXTRAFULLTABLEDEGENNAME "eosextradegennew.dat"
// SIMPLETABLE:
#define EOSSIMPLEHEADNAME "eossimplenew.head"
#define EOSSIMPLETABLENAME "eossimplenew.dat"
#define EOSSIMPLETABLEDEGENNAME "eossimpledegennew.dat"
// EXTRASIMPLETABLE:
#define EOSEXTRASIMPLEHEADNAME "eosextrasimplenew.head"
#define EOSEXTRASIMPLETABLENAME "eosextrasimplenew.dat"
#define EOSEXTRASIMPLETABLEDEGENNAME "eosextrasimpledegennew.dat"
// SIMPLEZOOMTABLE:
#define EOSSIMPLEZOOMHEADNAME "eossimplezoomnew.head"
#define EOSSIMPLEZOOMTABLENAME "eossimplezoomnew.dat"
#define EOSSIMPLEZOOMTABLEDEGENNAME "eossimplezoomdegennew.dat"
// EXTRASIMPLEZOOMTABLE:
#define EOSEXTRASIMPLEZOOMHEADNAME "eosextrasimplezoomnew.head"
#define EOSEXTRASIMPLEZOOMTABLENAME "eosextrasimplezoomnew.dat"
#define EOSEXTRASIMPLEZOOMTABLEDEGENNAME "eosextrasimplezoomdegennew.dat"

