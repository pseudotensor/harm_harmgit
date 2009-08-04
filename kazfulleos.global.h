

// TODO:
//
// 1) choose case where checking if within EOS is not done (assume always there and truncate instead of extend)
//    TWO if's removed then: if(iswithin_eostable) and in each type of calculation the call: if(get_eos_fromtable()) are removed : might speed things up
// 2) Is my interpolation best/fastest/correct?  Thompson et al. (2003) use bivariate interpolation that looks alot more complicated


// TODO NEW:
// 1) New H calculation (MAKE MPI)
// 3) check additional code that's presently in debug mode
// 4) check MPI stuff here and for gravity
// 5) change bounds so loops over bounding vars rather than only up to B3
// 6) rescale() depends often on final quantity being B3

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

// assume degen table is always stored along with corresponding normal table
#define NUMTBLS 3
#define FULLTABLE 0
#define SIMPLETABLE 1
#define SIMPLEZOOMTABLE 2


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



///////////////////////
//
// these are fixed in order and number from read-in file

// below used generally in code
#define PofRHOU 0      // p(rho0,u)
#define UofRHOP 1      // u(rho0,p)
// below used for simple dissipation entropy-inversion tracking same fluid element as energy-based inversion
#define UofRHOS 2      // U(rho0,Sden)

// below used for 5D inversion and sources.c
#define DPDRHOofRHOU 3 // dpdrho0 |u (rho0,u)
#define DPDUofRHOU 4   // dp/du |rho0 (rho0,u)

// below used for wave speeds for Riemann solution's dissipation term
#define CS2ofRHOU 5    // cs^2(rho0,u)

// below used for utoprim.orig.c (dudp_calc.c) entropy inversion
#define SofRHOU 6      // S(rho0,u)
#define DSDRHOofRHOU 7 // dS/drho0 |u (rho0,u)
#define DSDUofRHOU 8   // dS/du |\rho0 (rho0,u)

// below used for utoprim_jon.c entropy inversion
#define SSofRHOCHI 9      // specificS(rho0,\chi)
#define DSSDRHOofRHOCHI 10 // dspecificS/drho0 |\chi (rho0,\chi)
#define DSSDCHIofRHOCHI 11   // dspecificS/d\chi |\rho0 (rho0,\chi)

// below used for utoprim_jon.c hot inversion
#define PofRHOCHI 12    // P(rho0,\chi)  \chi = u+p
#define IDRHO0DP 13     // 1/(d\rho0/dp) |\chi
#define IDCHIDP 14      // 1/(d\chi/dp) |\rho0

// below used to define temperature and if <invalidtemperature this indicates not in valid part of table (i.e. interpolation from T->u leads to general u range that is mapped onto fixed u range, so final table bounds original table)
#define TEMPU 15 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPP 16 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPCHI 17 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPS 18 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)

// last index of quantities above
#define LASTEOSQUANTITIESBASE TEMPS



// below ending # corresponds to whichrnpmethod
// monotonized is 22 + extras
// for input, total is 7 + extras
// totals for input are 8 23 16 is present original data for version1,2,3, then eos_extract adds:
// for normal table, eos_extract.m adds 3 (degens) + 3 (tk's) = 6 total normal added
// totals for input after eos_extract are: 7+extra+6
// eos_extract always makes 24 + extras = 10 iterators + 14 eosquantities + extras
// for degen table, eox_extract.m has 9 total

// number of base quantities to *store* from table made by eos_extract.m starting from PofRHOU
#define NUMEOSQUANTITIESBASE (1+LASTEOSQUANTITIESBASE) // for memory, so 1+lastindex

// for memory optimization, specifiy which datatype
#define WHICHDATATYPEGENERAL 4


#if(WHICHEOSDIMEN==4 && WHICHDATATYPEGENERAL!=4)
#error WHICHEOSDIMEN and WHICHDATATYPEGENERAL inconsistent
#endif

// for different datatypes have different extra things
#define MAXNUMDATATYPES 4
#define NUMEOSQUANTITIESTYPE1 (NUMEOSQUANTITIESBASE + 1)  // (Full EOS, rho,T,H)
#define NUMEOSQUANTITIESTYPE2 (NUMEOSQUANTITIESBASE + 16) // (non-neutrino EOS, rho,T,Y_e)
#define NUMEOSQUANTITIESTYPE3 (NUMEOSQUANTITIESBASE + 11) // (Full EOS, rho,T,Y_e,Y_nu but H-dependent)
#define NUMEOSQUANTITIESTYPE4 (NUMEOSQUANTITIESBASE + 24) // (non-neutrino EOS, rho,T,Y_e,Y_\nu)

#if(WHICHDATATYPEGENERAL==3)
#define NUMEOSQUANTITIESMEM NUMEOSQUANTITIESTYPE3
#else
#define NUMEOSQUANTITIESMEM (MAX(MAX(MAX(NUMEOSQUANTITIESTYPE1,NUMEOSQUANTITIESTYPE2),NUMEOSQUANTITIESTYPE3),NUMEOSQUANTITIESTYPE4))
#endif


// so-called "extras" in eos_extract.m: Those things that didn't require extra processing as independent variables or derivatives -- just interpolated from T -> U only
// extras:

#define EXTRA1  LASTEOSQUANTITIESBASE+1
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

// used to map request to correct EXTRA for given table that uses certain whichdatatype
#define LAMBDATOT -100
#define QDOTNU -101


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




// processed quantities
// MAXPROCESSEDEXTRAS is set in nondepnmenmics.h
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



// degen offset quantities for independent variables utot,ptot, chi
#define NUMEOSDEGENQUANTITIESMEM (4) // utotdiff, ptotdiff, chidiff, stotdiff

#define NUMEOSDEGENINDEPS (NUMINDEPDIMENS-1) // rho, tdynorye, tdynorynu, H

#define NUMEOSDEGENQUANTITIESBASE (NUMEOSDEGENINDEPS*2)
#define NUMEOSDEGENQUANTITIES (NUMEOSDEGENQUANTITIESBASE+NUMEOSDEGENQUANTITIESMEM)

// independents for degen table
#define RHOEOSDEGEN 0
#define YEEOSDEGEN 1
#define YNUEOSDEGEN 2
#define HEOSDEGEN 3


// these are fixed in order and number from read-in file
// so can access functional degentable of independent variables by whichindep-1
#define UTOTDIFF (UEOS-1)   // should always resolve to: 0
#define PTOTDIFF (PEOS-1)   // :1
#define CHIDIFF  (CHIEOS-1) // :2
#define STOTDIFF (SEOS-1)   // :3







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






///////////////
//
// EOSN1=rho0
// EOSN2=u or p or chi
// EOSN3=tdynorye in seconds or Y_e
// EOSN4=tdynorynu in seconds or Y_\nu
// EOSN5=height in cm

// full EOS table w/ neutrino part to be computed during run-time
#define EOSN1 100
#define EOSN2 50
#define EOSN3 50
#define EOSN4 2
#define EOSN5 1         // H not tabulated

// EOS table with assumptions
// If EOSSIMPLEN5==1, then assumed height such that optically thin problem
// If EOSSIMPLEN3==1, then assumed Y_e(\rhob,T) (not usual except perhaps for testing)
// degenerate table assumes to be same size except N2=1
#define EOSSIMPLEN1 1
#define EOSSIMPLEN2 1
#define EOSSIMPLEN3 1
#define EOSSIMPLEN4 1   // here Y_\nu~0 (optically thin)
#define EOSSIMPLEN5 1   // H not tabulated

// EOS table with assumed Height (small) and assumed tdynorye (large)
// degenerate table assumes to be same size except N2=1
// NOT USED RIGHT NOW
#define EOSSIMPLEZOOMN1 1
#define EOSSIMPLEZOOMN2 1
#define EOSSIMPLEZOOMN3 1
#define EOSSIMPLEZOOMN4 1
#define EOSSIMPLEZOOMN5 1

// GODMARK: could have a table for Ynu=thermalized and have an array that stores when source term forces Ynu to be perfectly thermal, and use that table in that case.
// generating it now




///////////////
//
// Table limits

#define UPDOWN 2 // 0=down 1=up

#define TBLITEMS (UPDOWN+2+2)
#define TBLLINEARITEMS (UPDOWN+2)




////////////////////////////////////
//
// Some global position variables used to determine EOS
//

// should be 4
#define NUMHDIRECTIONS 4

#define NUMEXTRAINDEP (1+2) // other besides rho and u-versions [i.e. 1 is just offset for memory to start at index=1.  2 is for rho,u]

// NOTE: must be in same order and number as EOS independent vars
// GODMARK: must also change MAXPARLIST in nondepnmemonics.h
#define NUMNONSTANDARD (6+3+3)
#define NUMEOSGLOBALS (NUMEXTRAINDEP+NUMNONSTANDARD)   // number of per CPU position-based data for EOS

// these should be ordered and numbered such that correspond to EOS table independent variables
// rho, u/p/chi/sden, Y_e, Y_\nu, H
// 1    2              3     4    5 == 5D
// do NOT correspond to expanded independent variables list from EOS as read-in (i.e. not RHOEOS, UEOS, PEOS, CHIEOS, SEOS,  YEEOS, YNUEOS,  HEOS)
#define FIRSTEOSGLOBAL (NUMEXTRAINDEP) // 1=rho, 2=u/p/chi/sden, 3=Y_e, 4=Y_\nu 5=H
#define TDYNORYEGLOBAL (FIRSTEOSGLOBAL)             // Tdyn or Y_e depending upon whichrnpmethod
#define YNUGLOBAL (TDYNORYEGLOBAL+1) // Tdyn or Y_\nu depending upon whichynumethod
#define HGLOBAL (YNUGLOBAL+1)        // scale-height (for method that uses this for EOS, some averaged version of H
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


// tolerance to check whether repeated case for i,j,k,rho0,u
#define OLDTOLERANCE (1E-14)

// tolerance for checks on input values of table
#define TABLETOL (1E-14)

// value of read-in temperature such that below this is treated as indicating an invalid (rho0,u) EOS pair
// actual read-in value is 1E-20, but using 5E-20 guarantees no machine-error choices and works with floats too
// also, generally is more accurate as temperature since problems with inversion are near T~0
#define INVALIDTEMP (5E-20)


#define EOSHEADNAME "eosnew.head"
#define EOSTABLENAME "eosnew.dat"
#define EOSDEGENTABLENAME "eosdegennew.dat"

#define EOSSIMPLEHEADNAME "eossimplenew.head"
#define EOSSIMPLETABLENAME "eossimplenew.dat"
#define EOSDEGENSIMPLETABLENAME "eosdegensimplenew.dat"

#define EOSSIMPLEZOOMHEADNAME "eossimplezoomnew.head"
#define EOSSIMPLEZOOMTABLENAME "eossimplezoomnew.dat"
#define EOSDEGENSIMPLEZOOMTABLENAME "eosdegensimplezoomnew.dat"
