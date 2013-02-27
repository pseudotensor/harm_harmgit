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
// note these start at 1 not 0, so memory should always be allocated with NUMINDEPDIMENS+1
#define NUMINDEPDIMENSMEM (NUMINDEPDIMENS+1)
//
// below fed into procedures as function arguments
#define NONEINDEP (-1000)
#define RHOINDEP 1
#define TEMPLIKEINDEP 2
// below stored in EOSextra[]
#define YEINDEP 3
#define YNUINDEP 4
#define HINDEP 5
#define FIRSTINDEPDIMEN RHOINDEP
#define LASTINDEPDIMENUSED (MIN(WHICHEOSDIMEN,HINDEP)) // only works if rho,u,ye,ynu,h start at 1
#define LASTINDEPDIMEN HINDEP

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
// which temperature-like independent is chosen using vartypearraylocal[2]
#define FIRSTTKLIKE UEOS
#define LASTTKLIKE  SEOS





////////////////////
//
// Setup HARM EOS version of degeneracy table
//
///////////////////

// for seteostable() and kazfulleos.c
#define ISNOTDEGENTABLE 0
#define ISDEGENTABLE 1


// utotdegencut==DEGENCUTLASTOLDVERSION is last of old versions of cutting that only used utotoffset
#define DEGENCUTLASTOLDVERSION 1

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
#define NOSUCHDIFF (-1)
#define UTOTDIFF (UEOS-1)   // should always resolve to: 0
#define PTOTDIFF (PEOS-1)   // :1
#define CHIDIFF  (CHIEOS-1) // :2
#define STOTDIFF (SEOS-1)   // :3

#define TTOTDIFF (-1) // SUPERGODMARK: Not yet.

#define SSPECTOTDIFF (STOTDIFF+1) // :4 : just label, not used except as label and not used for memory items since nothing is function of SSPECTOTDIFF

// fastest indexed quantities in degen table
#define NUMEOSDEGENQUANTITIESMEM2 (3)
// table is accessed via: (like) r= R0 + exp(x1[i]), with r ranging from Rin to Rout, then:

// below used for table lookup in case want to access degen quantities directly for some reason
#define NUMEOSDEGENQUANTITIESMEMNEW NUMEOSDEGENQUANTITIESMEM2
#define NUMEOSDEGENQUANTITIESMEMOLD 1


#define EOSOFFSET 0 // like R0
#define EOSIN 1 // like Rin
#define EOSOUT 2 // like Rout
#define FIRSTEOSDEGEN EOSOFFSET
#define LASTEOSDEGEN EOSOUT
#define FIRSTEOSQUANTITY FIRSTEOSDEGEN // FIRST in list of quantities can grab per whichd



// for reading-in degen table
#define NUMEOSDEGENQUANTITIESNOTSTOREDin (NUMEOSDEGENINDEPS*2) // not stored degen quantities

#define NUMEOSDEGENQUANTITIESin (NUMEOSDEGENQUANTITIESMEM1*NUMEOSDEGENQUANTITIESMEM2) // stored degen quantities
// labels for inputting degen table (should be NUMEOSDEGENQUANTITIESin of them)
// must have offset's as first 4 in order in order for set_arrays_eostable(), etc. to be correct
#define UTOTOFFSETin 0
#define PTOTOFFSETin (UTOTOFFSETin+1)
#define CHIOFFSETin  (UTOTOFFSETin+2)
#define STOTOFFSETin (UTOTOFFSETin+3)
#define UTOTINin NUMEOSDEGENQUANTITIESMEM1
#define PTOTINin (UTOTINin+1)
#define CHIINin  (UTOTINin+2)
#define STOTINin (UTOTINin+3)
#define UTOTOUTin NUMEOSDEGENQUANTITIESMEM1*2
#define PTOTOUTin (UTOTOUTin+1)
#define CHIOUTin  (UTOTOUTin+2)
#define STOTOUTin (UTOTOUTin+3)


///////////////////////
//
// Setup INPUT for normal (non-degen) table macro names
// these are FIXED in order and number from *read-in* file
//
////////////////////////

// SUPERNOTE: all tables should have the non-stored data used for checking table consistency
#define NUMEOSQUANTITIESNOTSTOREDin (NUMINDEPDIMENS+NUMEOSINDEPS+NUMEOSDEGENQUANTITIESMEM1)

// Stored-data labels:
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
#define NUMTEMPin NUMEOSDEGENQUANTITIESMEM1
#define TEMPUin 15 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPPin (TEMPUin+1) // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPCHIin (TEMPUin+2) // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPSin (TEMPUin+3) // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)


// last index of quantities above
#define FIRSTEOSQUANTITIESBASEin PofRHOUin
#define LASTEOSQUANTITIESBASEin TEMPSin

// number of base quantities to *store* from table made by eos_extract.m starting from PofRHOU
#define NUMEOSQUANTITIESBASEin (1+LASTEOSQUANTITIESBASEin) // for memory, so 1+lastindex

// so-called "extras" in eos_extract.m: Those things that didn't require extra processing as independent variables or derivatives -- just interpolated from T -> U only
// extras:

#define EXTRA1in  (LASTEOSQUANTITIESBASEin+1)
#define EXTRA2in  (EXTRA1in+1)
#define EXTRA3in  (EXTRA1in+2)
#define EXTRA4in  (EXTRA1in+3)
#define EXTRA5in  (EXTRA1in+4)
#define EXTRA6in  (EXTRA1in+5)
#define EXTRA7in  (EXTRA1in+6)
#define EXTRA8in  (EXTRA1in+7)
#define EXTRA9in  (EXTRA1in+8)
#define EXTRA10in (EXTRA1in+9)
#define EXTRA11in (EXTRA1in+10)
#define EXTRA12in (EXTRA1in+11)
#define EXTRA13in (EXTRA1in+12)
#define EXTRA14in (EXTRA1in+13)
#define EXTRA15in (EXTRA1in+14)
#define EXTRA16in (EXTRA1in+15)
#define EXTRA17in (EXTRA1in+16)
#define EXTRA18in (EXTRA1in+17)
#define EXTRA19in (EXTRA1in+18)
#define EXTRA20in (EXTRA1in+19)
#define EXTRA21in (EXTRA1in+20)
#define EXTRA22in (EXTRA1in+21)
#define EXTRA23in (EXTRA1in+22)
#define EXTRA24in (EXTRA1in+23)

#define FIRSTEXTRAin EXTRA1in
#define LASTEXTRAin EXTRA24in



//////////////////////////////////
// for "pure" extra table:
// when reading-in "pure" extra table
// only used to convert to canonical format for these variables
#define TEMPUinextra 0
#define TEMPPinextra 1
#define TEMPCHIinextra 2
#define TEMPSinextra 3
#define FIRSTEOSQUANTITIESBASEinextra TEMPUinextra
#define LASTEOSQUANTITIESBASEinextra TEMPSinextra
#define NUMEOSQUANTITIESBASEinextra (LASTEOSQUANTITIESBASEinextra-FIRSTEOSQUANTITIESBASEinextra+1)

#define EXTRA1inextra (LASTEOSQUANTITIESBASEinextra+1)
#define EXTRA24inextra (EXTRA1inextra+23)
#define FIRSTEXTRAinextra EXTRA1inextra
#define LASTEXTRAinextra EXTRA24inextra





// for memory optimization, must specifiy which datatype
// so input whichdatatype no longer can be different than this!  Too complicated to allow different whichdatatypes at once.
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
// number of tables to track for repeated lookup
#define NUMEXTRATABLETYPES 2
#define ISNOTEXTRATABLETYPE 0
#define ISEXTRATABLETYPE 1

// assume degen table is always stored along with corresponding normal table
#if(ALLOWFULLTABLE&&ALLOWSIMPLETABLE&&ALLOWSIMPLEZOOMTABLE)
#define NUMTBLS 6
#elif(ALLOWFULLTABLE&&ALLOWSIMPLETABLE==0&&ALLOWSIMPLEZOOMTABLE==0)
#define NUMTBLS 2 // more memory friendly
#elif(ALLOWFULLTABLE&&ALLOWSIMPLETABLE&&ALLOWSIMPLEZOOMTABLE==0)
#define NUMTBLS 4 // more memory friendly
#else
#error "No setup for NUMTBLS"
#endif

#define NOTABLE (-1) // just used to indicate no table setup, not to be iterated over, so don't include in NUMTBLS
#define FULLTABLE 0
#define FULLTABLEEXTRA 1
#define SIMPLETABLE 2
#define SIMPLETABLEEXTRA 3
#define SIMPLEZOOMTABLE 4
#define EXTRASIMPLEZOOMTABLE 5

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
#define SUBTYPEEXTRA 9 // extra is separate table when WHICHDATATYPEGENERAL==4



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

#define UofRHOT -1 // SUPERGODMARK: Not setup yet.

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
#define FIRSTEOSSDEN SofRHOU
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
#define EXTRA2  (EXTRA1+1)
#define EXTRA3  (EXTRA1+2)
#define EXTRA4  (EXTRA1+3)
#define EXTRA5  (EXTRA1+4)
#define EXTRA6  (EXTRA1+5)
#define EXTRA7  (EXTRA1+6)
#define EXTRA8  (EXTRA1+7)
#define EXTRA9  (EXTRA1+8)
#define EXTRA10 (EXTRA1+9)
#define EXTRA11 (EXTRA1+10)
#define EXTRA12 (EXTRA1+11)
#define EXTRA13 (EXTRA1+12)
#define EXTRA14 (EXTRA1+13)
#define EXTRA15 (EXTRA1+14)
#define EXTRA16 (EXTRA1+15)
#define EXTRA17 (EXTRA1+16)
#define EXTRA18 (EXTRA1+17)
#define EXTRA19 (EXTRA1+18)
#define EXTRA20 (EXTRA1+19)
#define EXTRA21 (EXTRA1+20)
#define EXTRA22 (EXTRA1+21)
#define EXTRA23 (EXTRA1+22)
#define EXTRA24 (EXTRA1+23)

#define FIRSTEOSEXTRA EXTRA1
#define LASTEOSEXTRA EXTRA24


#define LASTEOSQUANTITY LASTEOSEXTRA // last in list of quantities can grab per whichd

// number of quantities can look up per whichd
#define NUMEOSQUANTITIES (LASTEOSQUANTITY-FIRSTEOSQUANTITY+1)
#define NUMEOSQUANTITIESMEM (LASTEOSQUANTITY+1) // don't offset variables that use this, but want to allow first quantity to be non-zero, so use this amount (more than required)


// used to map request to correct EXTRA for given table that uses certain whichdatatype
#define LAMBDATOT (-100)
#define QDOTNU (-101)


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
//// EXTRA12: Ynuthermal0 // GODMARK: Should be added if doing this!
#define DATATYPE3_EXTRAFINAL EXTRA11 // should be EXTRA12! GODMARK

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
// EXTRA19:  lambdatot 3E-34 // coli=18
// EXTRA20:  lambdaintot 3E-34
// EXTRA21:  tauphotonohcm 2E9
// EXTRA22:  tauphotonabsohcm 2E9
// EXTRA23:  nnueth0
// EXTRA24:  nnuebarth0
#define DATATYPE4_EXTRAFINAL EXTRA24




// names for processed quantities
// most are based upon optical depth correction, YNUTHERMAL0 is not.
#define NUMPROCESSEDEXTRAS 14 // for using get_extrasprocessed().
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
#define YNUTHERMAL0 10
#define ENUAVG 11
#define ENUE 12
#define ENUEBAR 13

// Below for dumping processed quantities.  Needs to be constant (not changing alot) so consistent dump files
#define MAXPROCESSEDEXTRAS (NUMPROCESSEDEXTRAS)

#if(MAXPROCESSEDEXTRAS<NUMPROCESSEDEXTRAS)
#error "Need to make MAXPROCESSEDEXTRAS larger."
#endif
