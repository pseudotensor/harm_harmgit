

// TODO:
//
// 1) choose case where checking if within EOS is not done (assume always there and truncate instead of extend)
//    TWO if's removed then: if(iswithin_eostable) and in each type of calculation the call: if(get_eos_fromtable()) are removed : might speed things up
// 2) Is my interpolation best/fastest/correct?  Thompson et al. (2003) use bivariate interpolation that looks alot more complicated


// TODO NEW:
// 1) New H calculation (MAKE MPI)
// 3) check additional code that's presently in debug mode
// 4) check MPI stuff here and for gravity

// whether to allow Kaz EOS table
#define ALLOWKAZEOS 0 // expensive for OpenMP due to many large globals, so normally disable unless required

//////////////////////////////
//
// Some often changable variables
//
///////////////////////////////

// whether to allow use of full table (if 0, then others must be turned on)
#define ALLOWFULLTABLE 1

// whether to only use full table (0) for allow use of simple tables if can (1)
#define ALLOWSIMPLETABLE 0

// whether to use simplezoom table if can
// zoom not needed anymore with new degen offset method
// indeed, set to 0 unless make new simplezoom table
#define ALLOWSIMPLEZOOMTABLE 0

// how many dimensions to consider.  Other dimensions' values will be consider as the dimension's lowest value
#define WHICHEOSDIMEN 4


// whether to use degen offset (otherwise assume degen offset from file is 0 even if read-in differently)
#define ALLOWDEGENOFFSET 1


// whether to check if table returns a valid EOS value by using existence of stored inversion to temperature
// so far only setup for F(rho0,u)
// if some invalid, then don't use those data points.  If all surrounding points are invalid, then use them as if they were valid
#define CHECKIFVALIDEOSDATA 1


// using log interpolation results in much smoother results, and consistent with eos_extract.m for interpolation
// That is, using integer position is log-interp since all independents are log on the grid
// And those functions in eos_extract.m interpolated as log are here interpolated as log
// 0 or 1
#define DOLOGINTERP 1


// which EOS to reduce to if beyond table
// mignone doesn't make sense
#define REDUCE2WHICHEOS IDEALGAS
// ensure that gamideal is chosen
// GODMARK: Could choose nearest tabulated value of dp/du|rho0 and dp/dchi|rho0 for gamideal when indeps are rho0,u and rho0,chi


//////////////////////////////
//
// More permenant variables
//
///////////////////////////////



//////////////
//
// here:
// [independent variables]
// rhob = rho0 -- rest-mass density in g/cc
// tk = Temperature in Kelvin
// hcm = height of medium in cm
// tdynorye = dynamical time in seconds (assumed NSE time) OR Y_e
//
// [dependent variables]
// p_tot = total pressure in cgs units
// u_tot = u = internal relativistic comoving energy (no rest-mass) in cgs units
// s_tot = total entropy density in comoving frame in cgs units



// 5 true dimensions (rhob,u/p/chi/sden/(tk),tdynorye,tdynorynu,hcm)
// Recall that u/p/chi/sden/tk are all in same independent variable slot for different functions
#define NUMINDEPDIMENS 5
// what is contained within file (rhob,utotdiff/ptotdiff/chidiff/stotdiff,tdynorye,tdynorynu,hcm)


#define NUMEOSINDEPS 8 // independent-like quantities that have read-in table limits [i.e. rho,u,p,chi,s,T,ynu,H]

// these are fixed in order and number as consistent with what's read-in from file
#define RHOEOS 0    // rest-mass density
#define UEOS   1    // internal energy density: used for u
#define PEOS   2    // pressure energy density: used for p
#define CHIEOS 3    // enthalpy energy density: used for \chi
#define SEOS   4    // entropy density: used for sden
#define YEEOS  5    // tdynorye: dynamical time or Y_e
#define YNUEOS 6    // tdynorynu: dynamical time or Y_\nu
#define HEOS   7    // height (not used normally anymore)

// begin and end types of temperature-originated quantities
// which temperature-like independent is chosen using vartypearray[2]
#define FIRSTTKLIKE UEOS
#define LASTTKLIKE  SEOS







/////////////////////////////
//
// Setup input file table names
//
/////////////////////////////


#define EOSHEADNAME "eosnew.head"
#define EOSTABLENAME "eosnew.dat"
#define EOSDEGENTABLENAME "eosdegennew.dat"

#define EOSSIMPLEHEADNAME "eossimplenew.head"
#define EOSSIMPLETABLENAME "eossimplenew.dat"
#define EOSDEGENSIMPLETABLENAME "eosdegensimplenew.dat"

#define EOSSIMPLEZOOMHEADNAME "eossimplezoomnew.head"
#define EOSSIMPLEZOOMTABLENAME "eossimplezoomnew.dat"
#define EOSDEGENSIMPLEZOOMTABLENAME "eosdegensimplezoomnew.dat"









////////////////////
//
// Setup HARM EOS version of degeneracy table
//
///////////////////


// rho, tdynorye, tdynorynu, H
// independents for degen table
// For accessing degen table N1,N3,N4,N5 (N2 skipped via use of vardegentablearray[]
#define NUMEOSDEGENINDEPS (NUMINDEPDIMENS-1)
#define RHOEOSDEGEN 0
#define YEEOSDEGEN 1
#define YNUEOSDEGEN 2
#define HEOSDEGEN 3

// degen offset quantities for independent variables utot,ptot, chi
// utotoffset,in,out, ptotoffset,in,out, chioffset,in,out, stotoffset,in,out
#define NUMEOSDEGENQUANTITIESMEM1 (4)
// these are fixed in order and number from read-in file
// so can access functional degentable of independent variables by whichindep-1
// below is also mapped to "whichd" variable
#define UTOTDIFF (UEOS-1)   // should always resolve to: 0
#define PTOTDIFF (PEOS-1)   // :1
#define CHIDIFF  (CHIEOS-1) // :2
#define STOTDIFF (SEOS-1)   // :3

// fastest indexed quantities in degen table
#define NUMEOSDEGENQUANTITIESMEM2 (3)
// table is accessed via: (like) r= R0 + exp(x1[i]), with r ranging from Rin to Rout, then:

#define EOSOFFSET 0 // like R0
#define EOSIN 1 // like Rin
#define EOSOUT 2 // like Rout
#define FIRSTEOSDEGEN EOSOFFSET
#define LASTEOSDEGEN EOSOUT
#define FIRSTEOSQUANTITY FIRSTEOSDEGEN // first in list of quantities can grab per whichd



// for reading-in degen table
#define NUMEOSDEGENQUANTITIESNOTSTOREDin (NUMEOSDEGENINDEPS*2) // not stored degen quantities

#define NUMEOSDEGENQUANTITIESin (NUMEOSDEGENQUANTITIESMEM1*NUMEOSDEGENQUANTITIESMEM2) // stored degen quantities
// labels for inputting degen table (should be NUMEOSDEGENQUANTITIESin of them)
// must have offset's as first 4 in order in order for set_arrays_eostable(), etc. to be correct
#define UTOTOFFSETin 0
#define PTOTOFFSETin UTOTOFFSETin+1
#define CHIOFFSETin  UTOTOFFSETin+2
#define STOTOFFSETin UTOTOFFSETin+3
#define UTOTINin NUMEOSDEGENQUANTITIESMEM1
#define PTOTINin UTOTINin+1
#define CHIINin  UTOTINin+2
#define STOTINin UTOTINin+3
#define UTOTOUTin NUMEOSDEGENQUANTITIESMEM1*2
#define PTOTOUTin UTOTOUTin+1
#define CHIOUTin  UTOTOUTin+2
#define STOTOUTin UTOTOUTin+3

// for seteostable() and kazfulleos.c
#define ISNOTDEGENTABLE 0
#define ISDEGENTABLE 1

// below used for table lookup in case want to access degen quantities directly for some reason
#define NUMEOSDEGENQUANTITIESMEM NUMEOSDEGENQUANTITIESMEM2

///////////////////////
//
// Setup INPUT for normal (non-degen) table macro names
// these are FIXED in order and number from *read-in* file
//
////////////////////////

#define NUMEOSQUANTITIESNOTSTOREDin (NUMINDEPDIMENS+NUMEOSINDEPS+NUMEOSDEGENQUANTITIESMEM1)

#define PofRHOUin 0      // p(rho0,u)
#define UofRHOPin 1      // u(rho0,p)
// below used for simple dissipation entropy-inversion tracking same fluid element as energy-based inversion
#define UofRHOSin 2      // U(rho0,Sden)

// below used for 5D inversion and sources.c
#define DPDRHOofRHOUin 3 // dpdrho0 |u (rho0,u)
#define DPDUofRHOUin 4   // dp/du |rho0 (rho0,u)

// below used for wave speeds for Riemann solution's dissipation term
#define CS2ofRHOUin 5    // cs^2(rho0,u)

// below used for utoprim.orig.c (dudp_calc.c) entropy inversion
#define SofRHOUin 6      // S(rho0,u)
#define DSDRHOofRHOUin 7 // dS/drho0 |u (rho0,u)
#define DSDUofRHOUin 8   // dS/du |\rho0 (rho0,u)

// below used for utoprim_jon.c entropy inversion
#define SSofRHOCHIin 9      // specificS(rho0,\chi)
#define DSSDRHOofRHOCHIin 10 // dspecificS/drho0 |\chi (rho0,\chi)
#define DSSDCHIofRHOCHIin 11   // dspecificS/d\chi |\rho0 (rho0,\chi)

// below used for utoprim_jon.c hot inversion
#define PofRHOCHIin 12    // P(rho0,\chi)  \chi = u+p
#define IDRHO0DPin 13     // 1/(d\rho0/dp) |\chi
#define IDCHIDPin 14      // 1/(d\chi/dp) |\rho0

// below used to define temperature and if <invalidtemperature this indicates not in valid part of table (i.e. interpolation from T->u leads to general u range that is mapped onto fixed u range, so final table bounds original table)
#define TEMPUin 15 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPPin 16 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPCHIin 17 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPSin 18 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)

// last index of quantities above
#define LASTEOSQUANTITIESBASEin TEMPSin

// number of base quantities to *store* from table made by eos_extract.m starting from PofRHOU
#define NUMEOSQUANTITIESBASEin (1+LASTEOSQUANTITIESBASEin) // for memory, so 1+lastindex

// so-called "extras" in eos_extract.m: Those things that didn't require extra processing as independent variables or derivatives -- just interpolated from T -> U only
// extras:

#define EXTRA1in  LASTEOSQUANTITIESBASEin+1
#define EXTRA2in  EXTRA1in+1
#define EXTRA3in  EXTRA1in+2
#define EXTRA4in  EXTRA1in+3
#define EXTRA5in  EXTRA1in+4
#define EXTRA6in  EXTRA1in+5
#define EXTRA7in  EXTRA1in+6
#define EXTRA8in  EXTRA1in+7
#define EXTRA9in  EXTRA1in+8
#define EXTRA10in EXTRA1in+9
#define EXTRA11in EXTRA1in+10
#define EXTRA12in EXTRA1in+11
#define EXTRA13in EXTRA1in+12
#define EXTRA14in EXTRA1in+13
#define EXTRA15in EXTRA1in+14
#define EXTRA16in EXTRA1in+15
#define EXTRA17in EXTRA1in+16
#define EXTRA18in EXTRA1in+17
#define EXTRA19in EXTRA1in+18
#define EXTRA20in EXTRA1in+19
#define EXTRA21in EXTRA1in+20
#define EXTRA22in EXTRA1in+21
#define EXTRA23in EXTRA1in+22
#define EXTRA24in EXTRA1in+23

#define FIRSTEXTRAin EXTRA1in
#define LASTEXTRAin EXTRA24in


// below ending # corresponds to whichrnpmethod
// monotonized is 22 + extras
// for input, total is 7 + extras
// totals for input are 8 23 16 is present original data for version1,2,3, then eos_extract adds:
// for normal table, eos_extract.m adds 3 (degens) + 3 (tk's) = 6 total normal added
// totals for input after eos_extract are: 7+extra+6
// eos_extract always makes 24 + extras = 10 iterators + 14 eosquantities + extras
// for degen table, eox_extract.m has 9 total


// for memory optimization, specifiy which datatype
#define WHICHDATATYPEGENERAL 4


#if(WHICHEOSDIMEN==4 && WHICHDATATYPEGENERAL!=4)
#error WHICHEOSDIMEN and WHICHDATATYPEGENERAL inconsistent
#endif

// for different datatypes have different extra things
#define MAXNUMDATATYPES 4
// below used for both "in"put tables and HARM EOS tables
#define NUMEXTRAEOSQUANTITIESTYPE1 (1)  // (Full EOS, rho,T,H)
#define NUMEXTRAEOSQUANTITIESTYPE2 (16) // (non-neutrino EOS, rho,T,Y_e)
#define NUMEXTRAEOSQUANTITIESTYPE3 (11) // (Full EOS, rho,T,Y_e,Y_nu but H-dependent)
#define NUMEXTRAEOSQUANTITIESTYPE4 (24) // (non-neutrino EOS, rho,T,Y_e,Y_\nu)

// number of columns to read-in for a table
#define NUMEOSQUANTITIESTYPE1in (NUMEOSQUANTITIESBASEin+NUMEXTRAEOSQUANTITIESTYPE1)
#define NUMEOSQUANTITIESTYPE2in (NUMEOSQUANTITIESBASEin+NUMEXTRAEOSQUANTITIESTYPE2)
#define NUMEOSQUANTITIESTYPE3in (NUMEOSQUANTITIESBASEin+NUMEXTRAEOSQUANTITIESTYPE3)
#define NUMEOSQUANTITIESTYPE4in (NUMEOSQUANTITIESBASEin+NUMEXTRAEOSQUANTITIESTYPE4)

// for simple arrays handling input tables
#define MAXNUMEOSQUANTITIESin (MAX(MAX(MAX(NUMEOSQUANTITIESTYPE1in,NUMEOSQUANTITIESTYPE2in),NUMEOSQUANTITIESTYPE3in),NUMEOSQUANTITIESTYPE4in))



// maximum size of pipeline is maximum number of things stored during read-in
#define MAXEOSPIPELINE (MAXNUMEOSQUANTITIESin+NUMEOSDEGENQUANTITIESin)


///////////////////
//
// Setup HARM EOS tables
//
///////////////////

// assume degen table is always stored along with corresponding normal table
#define NUMTBLS 3
#define NOTABLE -1 // just used to indicate no table setup, not to be iterated over, so don't include in NUMTBLS
#define FULLTABLE 0
#define SIMPLETABLE 1
#define SIMPLEZOOMTABLE 2

// number of table subtypes per table
// Degeneracy is unique table compared to others because size is 1 for temperature-like dimension N2, but still include in list
#define NUMTABLESUBTYPES 10
#define SUBTYPEDEGEN 0
#define SUBTYPESTANDARD 1
#define SUBTYPEGUESS 2
#define SUBTYPEDISS 3
#define SUBTYPEDP 4
#define SUBTYPESDEN 5
#define SUBTYPESSPEC 6
#define SUBTYPEPOFCHI 7
#define SUBTYPETEMP 8
#define SUBTYPEEXTRA 9

///////////////////////
//
// Setup HARM EOS table names
// These can be in any order per table subtype except TEMP types
// Note that each of these tables should be function of *same* independent variables (except TEMP table where only one of them is to be accessed), while one should have tables (even with same independent variables) broken-up enough if some functions are basically not used together ever
//
////////////////////////
#define NUMEOSSTANDARDQUANTITIESMEM 2
#define PofRHOU (LASTEOSDEGEN+1)
#define CS2ofRHOU (PofRHOU+1)
#define FIRSTEOSSTANDARD PofRHOU
#define LASTEOSSTANDARD CS2ofRHOU

#define NUMEOSGUESSQUANTITIESMEM 1
#define UofRHOP (LASTEOSSTANDARD+1)
#define FIRSTEOSGUESS UofRHOP
#define LASTEOSGUESS UofRHOP

#define NUMEOSDISSQUANTITIESMEM 1
#define UofRHOS (LASTEOSGUESS+1)
#define FIRSTEOSDISS UofRHOS
#define LASTEOSDISS UofRHOS

#define NUMEOSDPQUANTITIESMEM 2
#define DPDRHOofRHOU (LASTEOSDISS+1)
#define DPDUofRHOU (LASTEOSDISS+2)
#define FIRSTEOSDP DPDRHOofRHOU
#define LASTEOSDP DPDUofRHOU

#define NUMEOSSDENQUANTITIESMEM 3
#define SofRHOU (LASTEOSDP+1)
#define DSDRHOofRHOU (LASTEOSDP+2)
#define DSDUofRHOU (LASTEOSDP+3)
#define FIRSTEOSSDEN DSDRHOofRHOU
#define LASTEOSSDEN DSDUofRHOU

#define NUMEOSSSPECQUANTITIESMEM 3
#define SSofRHOCHI (LASTEOSSDEN+1)
#define DSSDRHOofRHOCHI (LASTEOSSDEN+2)
#define DSSDCHIofRHOCHI (LASTEOSSDEN+3)
#define FIRSTEOSSSPEC SSofRHOCHI
#define LASTEOSSSPEC DSSDCHIofRHOCHI

#define NUMEOSPOFCHIQUANTITIESMEM 3
#define PofRHOCHI (LASTEOSSSPEC+1)
#define IDRHO0DP (LASTEOSSSPEC+2)
#define IDCHIDP (LASTEOSSSPEC+3)
#define FIRSTEOSPOFCHI PofRHOCHI
#define LASTEOSPOFCHI IDCHIDP

// TEMP? must be in same order as  UEOS, PEOS, CHIEOS, SEOS and DIFF versions for easy array access
#define NUMEOSTEMPQUANTITIESMEM2 1 // number per whichd
#define NUMEOSTEMPQUANTITIESMEM NUMEOSTEMPQUANTITIESMEM2
#define TEMPGEN (LASTEOSPOFCHI+1)
#define FIRSTEOSTEMP TEMPGEN
#define LASTEOSTEMP TEMPGEN

#define NUMEOSTEMPQUANTITIESMEM1 NUMEOSDEGENQUANTITIESMEM1 // total number over all whichd
#define TEMPU UTOTDIFF
#define TEMPP PTOTDIFF
#define TEMPCHI CHIDIFF
#define TEMPS STOTDIFF

#if(WHICHDATATYPEGENERAL==1)
#define NUMEOSEXTRAQUANTITIES NUMEXTRAEOSQUANTITIESTYPE1
#elif(WHICHDATATYPEGENERAL==2)
#define NUMEOSEXTRAQUANTITIES NUMEXTRAEOSQUANTITIESTYPE2
#elif(WHICHDATATYPEGENERAL==3)
#define NUMEOSEXTRAQUANTITIES NUMEXTRAEOSQUANTITIESTYPE3
#elif(WHICHDATATYPEGENERAL==4)
#define NUMEOSEXTRAQUANTITIES NUMEXTRAEOSQUANTITIESTYPE4
#endif


#define NUMEOSEXTRAQUANTITIESMEM NUMEOSEXTRAQUANTITIES

// maximum number of extra variables in kazfulleos.c
// some constant so dump files don't change
#define MAXNUMEXTRAS 24

#if(MAXNUMEXTRAS<NUMEOSEXTRAQUANTITIES)
#error "Need to make MAXNUMEXTRAS larger."
#endif


// up to 24 extras so far with WHICHDATATYPE==4
#define EXTRA1  (LASTEOSTEMP+1)
#define EXTRA2  EXTRA1+1
#define EXTRA3  EXTRA1+2
#define EXTRA4  EXTRA1+3
#define EXTRA5  EXTRA1+4
#define EXTRA6  EXTRA1+5
#define EXTRA7  EXTRA1+6
#define EXTRA8  EXTRA1+7
#define EXTRA9  EXTRA1+8
#define EXTRA10 EXTRA1+9
#define EXTRA11 EXTRA1+10
#define EXTRA12 EXTRA1+11
#define EXTRA13 EXTRA1+12
#define EXTRA14 EXTRA1+13
#define EXTRA15 EXTRA1+14
#define EXTRA16 EXTRA1+15
#define EXTRA17 EXTRA1+16
#define EXTRA18 EXTRA1+17
#define EXTRA19 EXTRA1+18
#define EXTRA20 EXTRA1+19
#define EXTRA21 EXTRA1+20
#define EXTRA22 EXTRA1+21
#define EXTRA23 EXTRA1+22
#define EXTRA24 EXTRA1+23

#define FIRSTEOSEXTRA EXTRA1
#define LASTEOSEXTRA EXTRA24


#define LASTEOSQUANTITY LASTEOSEXTRA // last in list of quantities can grab per whichd

// number of quantities can look up per whichd
#define NUMEOSQUANTITIES (LASTEOSQUANTITY-FIRSTEOSQUANTITY+1)
#define NUMEOSQUANTITIESMEM (LASTEOSQUANTITY+1) // don't offset variables that use this, but want to allow first quantity to be non-zero, so use this amount (more than required)


// used to map request to correct EXTRA for given table that uses certain whichdatatype
#define LAMBDATOT -100
#define QDOTNU -101


// BELOW FOR HARM EOS TABLES, not inputs
// whichdatatype==1
// EXTRA1: Neutrino cooling rate (erg/s/cm^2)
#define DATATYPE1_EXTRAFINAL EXTRA1

// whichdatatype==2
// EXTRA1:  qtautelohcm
// EXTRA2:  qtauaelohcm
// EXTRA3:  qtautmuohcm
// EXTRA4:  qtauamuohcm
// EXTRA5:  qtauttauohcm
// EXTRA6:  qtauatauohcm
// EXTRA7:  ntautelohcm
// EXTRA8:  ntauaelohcm
// EXTRA9:  ntautmuohcm
// EXTRA10:  ntauamuohcm
// EXTRA11:  ntauttauohcm
// EXTRA12:  ntauatauohcm
// EXTRA13:  gammapeglobal+gammaAeglobal
// EXTRA14:  gammapnuglobal+gammapenuglobal
// EXTRA15:  gammanglobal + gammaneglobal
// EXTRA16: gammannuglobal
#define DATATYPE2_EXTRAFINAL EXTRA16

// whichdatatype==3
// for now this is opimal choice for simplicity, although big table unresolved in H and Ynu
// EXTRA1:  Qphoton
// EXTRA2:  Qm
// EXTRA3:  graddotrhouyl
// EXTRA4:  Tthermaltot
// EXTRA5:  Tdifftot
// EXTRA6:  lambdatot
// EXTRA7:  lambdaintot
// EXTRA8:  Enuglobal
// EXTRA9:  Enueglobal
// EXTRA10:  Enuebarglobal
// EXTRA11: Ynuthermal
#define DATATYPE3_EXTRAFINAL EXTRA11

// whichdatatype==4
// GODMARK: if this is going to work, need also the energy density parts as functions of \chi, but for now doesn't seem this method will be useful due to need to iterate within table to find du from u
// EXTRA1:  qtautnueohcm 5E7
// EXTRA2:  qtauanueohcm 1E-50
// EXTRA3:  qtautnuebarohcm 1E36
// EXTRA4:  qtauanuebarohcm 1E36
// EXTRA5:  qtautmuohcm 1E-9
// EXTRA6:  qtauamuohcm 1E-14
// EXTRA7:  ntautnueohcm 5E7
// EXTRA8:  ntauanueohcm 3E-48
// EXTRA9:  ntautnuebarohcm 2E35
// EXTRA10:  ntauanuebarohcm 2E35
// EXTRA11:  ntautmuohcm 1E-9
// EXTRA12:  ntauamuohcm 1E-14
// EXTRA13:  unue0 8E59
// EXTRA14:  unuebar0 3E-26
// EXTRA15:  unumu0 4E25
// EXTRA16:  nnue0 1E57
// EXTRA17:  nnuebar0 3E-20
// EXTRA18:  nnumu0 1E31
// EXTRA19:  lambdatot 3E-34
// EXTRA20:  lambdaintot 3E-34
// EXTRA21:  tauphotonohcm 2E9
// EXTRA22:  tauphotonabsohcm 2E9
// EXTRA23:  nnueth0
// EXTRA24:  nnuebarth0
#define DATATYPE4_EXTRAFINAL EXTRA24




// names for processed quantities
#define NUMPROCESSEDEXTRAS 13 // for using get_extrasprocessed().
#define QPHOTON 0
#define QNEUTRINO 1
#define GRADDOTRHOUYL 2
#define TTHERMAL 3
#define TDIFF 4
#define RHONU 5
#define PNU 6
#define SNU 7
#define YNULOCAL 8
#define YNUTHERMAL 9
#define ENUAVG 10
#define ENUE 11
#define ENUEBAR 12

// Below for dumping processed quantities.  Needs to be constant (not changing alot) so consistent dump files
#define MAXPROCESSEDEXTRAS 13

#if(MAXPROCESSEDEXTRAS<NUMPROCESSEDEXTRAS)
#error "Need to make MAXPROCESSEDEXTRAS larger."
#endif

















/////////////////////////////
//
// Setup INPUT and HARM EOS table sizes for N1,N2,N3,N4,N5.
// Both INPUT and HARM EOS table sizes are same for these dimensions associated with independent variables
//
// EOSN1=rho0
// EOSN2=u or p or chi
// EOSN3=tdynorye in seconds or Y_e
// EOSN4=tdynorynu in seconds or Y_\nu
// EOSN5=height in cm
//
////////////////////////////

// full EOS table w/ neutrino part to be computed during run-time
#define EOSFULLN1 200
#define EOSFULLN2 200
#define EOSFULLN3 50
#define EOSFULLN4 10
#define EOSFULLN5 1         // H not tabulated

#define EOSFULLDEGENN1 EOSFULLN1
#define EOSFULLDEGENN2 (1) // always 1
#define EOSFULLDEGENN3 EOSFULLN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSFULLDEGENN4 (1) // Then degen table is not function of neutrino term Y_\nu
#else
#define EOSFULLDEGENN4 EOSFULLN4
#endif
#define EOSFULLDEGENN5 EOSFULLN5 // should be 1 anyways

// EOS table with assumptions
// If EOSSIMPLEN5==1, then assumed height such that optically thin problem
// If EOSSIMPLEN3==1, then assumed Y_e(\rhob,T) (not usual except perhaps for testing)
// degenerate table assumes to be same size except N2=1
#define EOSSIMPLEN1 1
#define EOSSIMPLEN2 1
#define EOSSIMPLEN3 1
#define EOSSIMPLEN4 1   // here Y_\nu~0 (optically thin)
#define EOSSIMPLEN5 1   // H not tabulated

#define EOSSIMPLEDEGENN1 EOSSIMPLEN1
#define EOSSIMPLEDEGENN2 (1) // always 1
#define EOSSIMPLEDEGENN3 EOSSIMPLEN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEDEGENN4 (1)
#else
#define EOSSIMPLEDEGENN4 EOSSIMPLEN4
#endif
#define EOSSIMPLEDEGENN5 EOSSIMPLEN5 // should be 1

// EOS table with assumed Height (small) and assumed tdynorye (large)
// degenerate table assumes to be same size except N2=1
// NOT USED RIGHT NOW
#define EOSSIMPLEZOOMN1 1
#define EOSSIMPLEZOOMN2 1
#define EOSSIMPLEZOOMN3 1
#define EOSSIMPLEZOOMN4 1
#define EOSSIMPLEZOOMN5 1

#define EOSSIMPLEZOOMDEGENN1 EOSSIMPLEZOOMN1
#define EOSSIMPLEZOOMDEGENN2 (1)
#define EOSSIMPLEZOOMDEGENN3 EOSSIMPLEZOOMN3
#if(WHICHDATATYPEGENERAL==4)
#define EOSSIMPLEZOOMDEGENN4 (1)
#else
#define EOSSIMPLEZOOMDEGENN4 EOSSIMPLEZOOMN4
#endif
#define EOSSIMPLEZOOMDEGENN5 EOSSIMPLEZOOMN5

// GODMARK: could have a table for Ynu=thermalized and have an array that stores when source term forces Ynu to be perfectly thermal, and use that table in that case.
// generating it now




///////////////
//
// Table limits indices

#define UPDOWN 2 // 0=down 1=up

#define TBLITEMS (UPDOWN+2+2)
#define TBLLINEARITEMS (UPDOWN+2)






////////////////////////////////////
//
// Setup EOSextra[]
//
// Some global position variables used to determine EOS
//
/////////////////////////////////////

// should be 4
#define NUMHDIRECTIONS 4

#if(NUMHDIRECTIONS!=4)
#error "NUMHDIRECTIONS should be 4"
#endif



////////////////////////////////////
//
// Setup macros to variables in EOSextra[]
//
/////////////////////////////////////

// these should be ordered and numbered such that correspond to EOS table independent variables
// do NOT correspond to expanded independent variables list from EOS as read-in (i.e. not RHOEOS, UEOS, PEOS, CHIEOS, SEOS,  YEEOS, YNUEOS,  HEOS)
#define RHOGLOBAL (-2) // dummy reference
#define UGLOBAL (-1) // dummy reference
#define TDYNORYEGLOBAL (0)        // Tdyn or Y_e depending upon whichrnpmethod
#define YNU0GLOBAL     (TDYNORYEGLOBAL+1) // Tdyn or Y_\nu depending upon whichynumethod
#define YNU0OLDGLOBAL  (YNU0GLOBAL+1)     // Tdyn or Y_\nu depending upon whichynumethod
#define YNUOLDGLOBAL   (YNU0OLDGLOBAL+1)  // Tdyn or Y_\nu depending upon whichynumethod
#define HGLOBAL        (YNUOLDGLOBAL+1)   // scale-height (for method that uses this for EOS, some averaged version of H
#define H2GLOBAL (HGLOBAL+1)         // 2,3,4 are other directions for axisymmetric emission
#define H3GLOBAL (H2GLOBAL+1) 
#define H4GLOBAL (H3GLOBAL+1) 
#define UNUGLOBAL (H4GLOBAL+1)       // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define PNUGLOBAL (UNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define SNUGLOBAL (PNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
// below 3 used to indicate position if EOS coming from grid
#define IGLOBAL (SNUGLOBAL+1)
#define JGLOBAL (IGLOBAL+1)
#define KGLOBAL (JGLOBAL+1)

// Ye should always be first
// Note that RHO,U,etc. GLOBAL should resolve to 1,2,3,4,5
// below is whatever comes after RHO and U
#define FIRSTEOSGLOBAL (TDYNORYEGLOBAL)
#define LASTEOSGLOBAL (KGLOBAL)

// NOTE: must be in same order and number as EOS independent vars
// number of per CPU position-based data for EOS
//#define NUMEOSGLOBALS (1+3+4+3+3)
#define NUMEOSGLOBALS (LASTEOSGLOBAL-FIRSTEOSGLOBAL+1)
#define MAXPARLIST (NUMEOSGLOBALS) // for get_EOS_parms()





///////////////////////////
//
// Some constants, tolerances, etc.
//
///////////////////////////

// tolerance to check whether repeated case for i,j,k,rho0,u
#define OLDTOLERANCE (1E-14)

// tolerance for checks on input values of table
#define TABLETOL (1E-14)

// value of read-in temperature such that below this is treated as indicating an invalid (rho0,u) EOS pair
// actual read-in value is 1E-20, but using 5E-20 guarantees no machine-error choices and works with floats too
// also, generally is more accurate as temperature since problems with inversion are near T~0
#define INVALIDTEMP (5E-20)



// initialize kaziio, etc. with this so first call has no old index used
#define INITKAZINDEX -100



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// define EOS array macros to avoid some confusion with spatial array macros
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define GENEOSPOINT(prefix,name) prefix##name
#define GENEOSTABLEMAC(prefix,name,a1,a2,a3,a4,a5,a6,a7) prefix##name[a1][a2][a3][a4][a5][a6][a7]
//
#define EOSBASEPOINT(name) GENEOSPOINT(a_,name) // not used
#define BASEEOSMAC(name,a1,a2,a3,a4,a5,a6,a7) GENEOSTABLEMAC(a_,name,a1,a2,a3,a4,a5,a6,a7)
//
#define EOSPOINT(name) GENEOSPOINT(,name)
#define EOSMAC(name,a1,a2,a3,a4,a5,a6,a7) GENEOSTABLEMAC(,name,a1,a2,a3,a4,a5,a6,a7)
#define PTRDEFEOSMAC(name,a1,a2,a3,a4,a5,a6,a7) (*EOSPOINT(name))[a2][a3][a4][a5][a6][a7]
#define PTREOSMAC(name,a1,a2,a3,a4,a5,a6,a7) (*)[a2][a3][a4][a5][a6][a7]




// e.g. superdefs.h like: double BASEEOSMAC(name,....)
//      superdefs.pointers.h like: double PTRDEFEOSMAC(name,....)
// set_arrays_multidimen.c like: EOSPOINT(name) = (double PTREOSMAC(name,..)) (&(BASEEOSMAC(name,...)));
//
// \([_a-zA-Z0-9]+\)\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\]
//     -> BASEEOSMAC(\1,\2,\3,\4,\5,\6,\7) [kazfulleos.c at top and kazfulleos_set_arrays.c for most-RHS of pointer shifting code]
//  OR -> EOSMAC(\1,\2,\3,\4,\5,\6,\7) [kazfulleos.c in code]
//
// (\*) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *) *( *& *(\([_a-zA-Z0-9]+\)(\([_a-zA-Z0-9]+\),
//     -> (*)PTR\6(\7,FILL,\1,\2,\3,\4,\5)) (&(\6(\7,   [in kazfulleos_set_arrays.c]
// Then :
//     -> (\*)PTRBASEEOSMAC -> PTREOSMAC    [in kazfulleos_set_arrays.c]
// Then:
// Replace (e.g.) eostable with EOSPOINT
//
// 
// ( *\* *\([_a-zA-Z0-9+-]+\) *) *\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\]\[\([_a-zA-Z0-9+-]+\)\] *;
//     -> PTRDEFEOSMAC(\1,FILL,\2,\3,\4,\5,\6); [defining pointer in kazfulleos.c just after BASEEOSMAC defines global array]
//
