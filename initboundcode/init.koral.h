
/*! \file init.koral.h
    \brief USER KORAL/RAD choices/switches
    
*/

//undefs
#undef DOYFL
#undef MAXWELL
#undef TRACKVPOT
#undef EVOLVEWITHVPOT
#undef DOGRIDSECTIONING
#undef MERGEDC2EA2CMETHODEM
#undef MERGEDC2EA2CMETHODMA
#undef MERGEDC2EA2CMETHOD
#undef ACCURATESINCOS
#undef ACCURATESOURCEDIAG
#undef ACCURATEDIAG
#undef REMOTEHOST
#undef WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS
#undef WENO_USE_PRIM_REDUCTION
#undef LIMIT_FLUXC2A_PRIM_CHANGE
#undef COMPDIM
#undef SPLITNPR
#undef FIELDSTAGMEM
#undef HIGHERORDERMEM
#undef MAXBND
#undef MCOORD
#undef ALLOWMETRICROT
#undef PERCELLDT
#undef PRODUCTION
#undef FULLOUTPUT
#undef MAILWHENDONE
#undef EMAILMESSAGE
#undef EMAILADDRESS 
#undef PERFTEST
#undef DOAVG
#undef DOJETDIAG
#undef DOAVG2
#undef DODEBUG
#undef DO_WENO_DEBUG
#undef DOENODEBUG
#undef DODISS
#undef DOEVOLVEMETRIC
#undef EVOLVEMETRICSUBSTEP
#undef LIMITDTWITHSOURCETERM
#undef LIMITSOURCES
#undef USEGRAVITYDTINDTLIMIT
#undef RESTRICTDTSETTINGINSIDEHORIZON
#undef DOLUMVSR
#undef DODISSVSR
#undef DOSELFGRAVVSR
#undef DOFIELDLINE
#undef ROEAVERAGEDWAVESPEED
#undef ATHENAROE
#undef STOREWAVESPEEDS
#undef USESTOREDSPEEDSFORFLUX
#undef VCHARTYPE
#undef PRECISEINVERSION
#undef WHICHVEL
#undef WHICHEOM
#undef REMOVERESTMASSFROMUU
#undef RELTYPE
#undef EOMTYPE
#undef UTOPRIMTRYAGAIN
#undef WHICHEOS
#undef CHECKONINVERSION
#undef CHECKONINVERSIONRAD
#undef DOENTROPY
#undef WHICHENTROPYEVOLVE
#undef FIXUPAFTERINIT
#undef FIXUPAFTERRESTART
#undef CHECKSOLUTION
#undef GAMMAPERCDIFFMAX
#undef UPERCDIFFMAX
#undef DOEXTRAINTERP
#undef LIMADJUST
#undef HYDROLIMADJUSTONLY
#undef FLUXADJUST
#undef HYDROFLUXADJUSTONLY
#undef STEPOVERNEGU
#undef STEPOVERNEGRHO
#undef STEPOVERNEGRHOU
#undef UTOPRIMADJUST
#undef UTOPRIMFAILRETURNTYPE
#undef SMOOTHSING
#undef COORDSINGFIX
#undef SINGSMALL
#undef DOSTOREPOSITIONDATA
#undef CONNDERTYPE
#undef VOLUMEDIFF
#undef GDETVOLDIFF
#undef FIXGDETSPC_WHEN_1DRADIAL
#undef MINDT
#undef JONCHECKS
#undef JONCHECKS2
#undef FLOORDIAGS
#undef ANALYTICGCON
#undef ANALYTICCONNECTION
#undef ANALYTICSOURCE
#undef OUTFLOWAVOIDBC
#undef FLUXDIMENSPLIT
#undef A2CDIMENSPLIT
#undef DODQMEMORY
#undef BOUNDFLUXRECON
#undef DOENOFLUXMEMORY
#undef BOUNDARYINTERPADJUST
#undef COMPUTEFRDOT
#undef CALCFARADAYANDCURRENTS
#undef WHICHCURRENTCALC
#undef FARADAYT0
#undef CURRENTST0
#undef EVOLVECHECKS
#undef FIXUPZONES
#undef HLLBOUNDARY
#undef FIXUPFLUX
#undef ZEROOUTFLOWFLUX
#undef ZEROPOLEFLUX
#undef RESCALEINTERP
#undef BDIRCONT
#undef HYPERHLL
#undef HORIZONSUPERFAST
#undef VARTOINTERP
#undef USEAVGPRIMITIVEFORWENOFLAT
#undef USEPRIMITIVEFROMAVGCONSERVED
#undef CONTACTINDICATOR
#undef COMPUTEDRHODP
#undef SUPERFASTDIVREDUCE
#undef MINPREFORDER
#undef SHOCKINDICATOR
#undef WHICHPARA
#undef HOT2ENTROPY
#undef HOT2COLD
#undef ENTROPY2COLD

//undefine the grid size parameters if they have already been defined
#ifdef N1 
#undef N1
#endif

#ifdef N2 
#undef N2
#endif

#ifdef N3 
#undef N3
#endif







//****************************************//
//****************************************//
//****************************************//
//****************************************//
//****************************************//

#define DOYFL 2 // track floors as scalars

#define ACCURATESOURCEDIAG 2 // 2 means full component decomposition is accurate as well as sum.
#define ACCURATEDIAG 1 // 1 means fluxes are accurate


#define MAXWELL PRIMMAXWELL
#define TRACKVPOT 1 // now on by default
#define EVOLVEWITHVPOT 0 // not on by default
#define DOGRIDSECTIONING 0 // not on by default
#define MERGEDC2EA2CMETHODEM 0
#define MERGEDC2EA2CMETHODMA 0
#define MERGEDC2EA2CMETHOD 0
#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 1
//#define WENO_USE_LIM_PRIM_CORRECTION_FOR_FLUX_LIMITING 1
#define WENO_USE_PRIM_REDUCTION 1
#define LIMIT_FLUXC2A_PRIM_CHANGE 0




#define ALLOWMETRICROT 0 // WALD->1
#if(ALLOWMETRICROT==1)
#undef CONNAXISYMM
#define CONNAXISYMM 0 //required to be 0 if really rotating metric
#undef DOMIXTHETAPHI
#define DOMIXTHETAPHI 1
#endif

#undef DOPOLEDEATH
#undef DOPOLESMOOTH
#undef DOPOLEGAMMADEATH
// needed to avoid random death at pole at large distances when grid focuses on axis and so makes-up information a bit.
// no, causes injection of lots of radiation and violates energy-momentum conservation alot, even changes solution alot.
#define DOPOLEDEATH 1
//#define DOPOLEDEATH 0 // WALD
#define DOPOLESMOOTH 0 // GODMARK: Need to reject outliers
#define DOPOLEGAMMADEATH 1
//#define DOPOLEGAMMADEATH 0 // WALD
// Note that if DOPOLESMOOTH>=DOPOLEGAMMADEATH or DOPOLESMOOTH>=DOPOLEDEATH, then DOPOLEGAMMADEATH or DOPOLEDEATH do nothing -- they are overwritten by DOPOLESMOOTH.

#undef IF3DSPCTHENMPITRANSFERATPOLE
#define IF3DSPCTHENMPITRANSFERATPOLE 0 // need to reject outliers and use full 3D info before this is used again.  Otherwise (more) problems at r>rbr when hyperexpoential grid is used.

#define PERCELLDT 0
#define COMPDIM 3
#define SPLITNPR 0 // TESTING
#define FIELDSTAGMEM 1 // testing
#define HIGHERORDERMEM 0
#define MAXBND 4 // 4 for PARAFLAT, 6 for WENO5BND wo/a2c stuff : 11 for full point-field FLUXRECON method
#define PRODUCTION 1 // WALDPROD
//#define FULLOUTPUT MAXBND // TESTING BCs
#define FULLOUTPUT 0

#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "relativity.cfa.harvard.edu"
#define EMAILADDRESS "jmckinney@cfa.harvard.edu"
#define EMAILMESSAGE "Done with GRMHD run DEFAULT"
#define PERFTEST 0
#define DOAVG 0
#define DOJETDIAG 0
#define DOAVG2 0
#define DODEBUG 1
#define DO_WENO_DEBUG 0
#define DOENODEBUG 0
#define DOEVOLVEMETRIC 0
#define EVOLVEMETRICSUBSTEP 1 // evolve metric every substep
#define LIMITSOURCES 0
#define LIMITDTWITHSOURCETERM 0 // causes problems, drops dt too low
#define USEGRAVITYDTINDTLIMIT 0
#define RESTRICTDTSETTINGINSIDEHORIZON 2
#define DODISS 0
#define DOLUMVSR 0
#define DODISSVSR 0
#define DOSELFGRAVVSR 0
#define DOFIELDLINE 1
#define ROEAVERAGEDWAVESPEED 0
#define ATHENAROE 0

//set this and the following one to unity to use the DONOR interpolated states for computing wavespeeds
#if(1 || SPLITNPR==1 || FIELDSTAGMEM==1) // should also be on if FLUXB==FIELDSTAG
#define STOREWAVESPEEDS 2 // no choice
#else
#define STOREWAVESPEEDS 0 // choice
#endif

#define USESTOREDSPEEDSFORFLUX (STOREWAVESPEEDS>0) // choice really
#define VCHARTYPE VERYLOCALVCHAR
#define PRECISEINVERSION 1
#define WHICHVEL VELREL4
#define WHICHEOM WITHGDET
//#define WHICHEOM (ISSPCMCOORD(MCOORD) ? WITHNOGDET : WITHGDET) // now default is WITHNOGDET for normal problems -- assumes half or full \theta hemispheres since main benefit is near poles. // still seems wrong -- need to test.
#define REMOVERESTMASSFROMUU 2
#define RELTYPE RELEOM
#define EOMTYPE EOMGRMHD
#undef EOMRADTYPE
#define EOMRADTYPE EOMRADM1CLOSURE
//#define EOMTYPE EOMFFDE
//#define EOMTYPE EOMCOLDGRMHD
#define UTOPRIMTRYAGAIN 0
#define WHICHEOS IDEALGAS

#define CHECKONINVERSION 1 // can slow things down
#define CHECKONINVERSIONRAD 1 // can independently check radiation inversion when no corrections

#if(DODISS || DOLUMVSR || DODISSVSR)
// for diss: testing CHANGINGMARK
#define DOENTROPY DOEVOLVEENTROPY
#define WHICHENTROPYEVOLVE EVOLVEFULLENTROPY
#else
// no diss/entropy
#define DOENTROPY DONOENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY
#endif

// force entropy variable enabled so can use HOT2ENTROPY
#undef DOENTROPY
#define DOENTROPY DOEVOLVEENTROPY

#undef EVOLVENRAD
#define EVOLVENRAD 0

#define FIXUPAFTERINIT 1
#define FIXUPAFTERRESTART 1

#define CHECKSOLUTION 0 // can cause erratic behavior near BH -- when gamma jumps are relatively large this averages causing large heating -- could just use internal energy check

#define GAMMAPERCDIFFMAX (2.0)
#define UPERCDIFFMAX (1E3) // 10.0 too restrictive
#define LIMADJUST LIMITERFIXED
#define HYDROLIMADJUSTONLY 0
#define FLUXADJUST FLUXFIXED
#define HYDROFLUXADJUSTONLY 0
//#define STEPOVERNEGU NEGDENSITY_NEVERFIXUP
//#define STEPOVERNEGRHO NEGDENSITY_NEVERFIXUP
//#define STEPOVERNEGRHOU NEGDENSITY_NEVERFIXUP
//#define STEPOVERNEGU NEGDENSITY_FIXONFULLSTEP
//#define STEPOVERNEGRHO NEGDENSITY_FIXONFULLSTEP
//#define STEPOVERNEGRHOU NEGDENSITY_FIXONFULLSTEP

// to be thermodynamically consistent on sub-steps as required to make sense of TVD type RK2 or RK3 methods or RK4 method.
#define STEPOVERNEGU NEGDENSITY_ALWAYSFIXUP
#define STEPOVERNEGRHO NEGDENSITY_ALWAYSFIXUP
#define STEPOVERNEGRHOU NEGDENSITY_ALWAYSFIXUP

#define UTOPRIMADJUST UTOPRIMAVG
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED
#define SMOOTHSING 0 // near BH
#define COORDSINGFIX (1) // for FLUXB==FLUXCTSTAG
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
//#define SINGSMALL (1E-3)
#define SINGSMALL (10000*NUMEPSILON) // must be larger than machine precision to work for outer M_PI boundary!
// Hawley uses 0.06283 (0.02Pi)

#define DOSTOREPOSITIONDATA 1 // DEBUG
#define CONNDERTYPE DIFFGAMMIE // DEBUG
//#define CONNDERTYPE DIFFNUMREC
#define VOLUMEDIFF 0
#define GDETVOLDIFF 0 // doesn't help much
#define FIXGDETSPC_WHEN_1DRADIAL 1

#define MINDT 1.e-20 
#define JONCHECKS 1    //SASMARK - do I need this?
#define JONCHECKS2 1   //SASMARK - do I need this?
#define FLOORDIAGS 1
#define ANALYTICGCON 0
#define ANALYTICCONNECTION 0  //SASMARK - Don't I need this?
#define ANALYTICSOURCE 0
#define OUTFLOWAVOIDBC 0
#define FLUXDIMENSPLIT PERFECTUNSPLIT
#define A2CDIMENSPLIT PERFECTUNSPLIT
#define DODQMEMORY 1
#define BOUNDFLUXRECON 0 // can set this to 1 if want to bound fluxes instead for FLUXRECON method (may be useful near poles)
#define DOENOFLUXMEMORY 0
#define BOUNDARYINTERPADJUST 0  //should be set to zero always
#define COMPUTEFRDOT 0
#define CALCFARADAYANDCURRENTS 0 // WALD->1
#define WHICHCURRENTCALC CURRENTCALC1
#define FARADAYT0 1
#define CURRENTST0 1


#define EVOLVECHECKS 1
#define FIXUPZONES FIXUP1ZONE
#define HLLBOUNDARY 0
#define FIXUPFLUX 0
#define ZEROOUTFLOWFLUX 0
#define ZEROPOLEFLUX 0
#define BDIRCONT 1
#define HYPERHLL 0
#define HORIZONSUPERFAST 0

//#define VARTOINTERP PRIMTOINTERP
#define VARTOINTERP PRIMTOINTERP_GDETFULLVERSION
//#define VARTOINTERP PRIMTOINTERP_GDETFULLVERSION_WALD // WALD
//#define VARTOINTERP PRIMTOINTERP_RHOU
//#define VARTOINTERP PRIMTOINTERP_VSQ
// #define VARTOINTERP PRIMTOINTERP_3VELREL_GAMMAREL (used in Sasha tests)
#undef VARTOINTERPFIELD
//#define VARTOINTERPFIELD GDETVERSION
#define VARTOINTERPFIELD NOFIELDRESCALE
#define RESCALEINTERP 1
#undef RESCALEINTERPFLUXCTSTAG
#define RESCALEINTERPFLUXCTSTAG 1 // WALD: 0->1
#define DOEXTRAINTERP 0

#define USEAVGPRIMITIVEFORWENOFLAT 1
#define USEPRIMITIVEFROMAVGCONSERVED 0
#define CONTACTINDICATOR 0
#define COMPUTEDRHODP 1
#define SUPERFASTDIVREDUCE 0
#define MINPREFORDER 3
#define SHOCKINDICATOR 1
#define WHICHPARA PARA4

#undef DO_VORTICITY_IMAGE
#define DO_VORTICITY_IMAGE 0

#define HOT2ENTROPY 1
#define HOT2COLD 1
#define ENTROPY2COLD 1

#define ACCURATESINCOS 1

#undef FLIPGDETAXIS
#define FLIPGDETAXIS 1
//#define FLIPGDETAXIS 0

#undef BOUNDPLPR
#define BOUNDPLPR 0

#undef NOFLUXCTONX1DN
#define NOFLUXCTONX1DN 0

#undef NUMPANALYTICOTHER
#undef DODUMPOTHER

#define NUMPOTHER 0
#define DODUMPOTHER 0

#undef FLUXDUMP
#define FLUXDUMP 2

#undef OUTERRADIALSUPERFAST
#define OUTERRADIALSUPERFAST 0 // can be better, but can also be much worse.


struct Ccoordparams {
  double timescalefactor;
}  coordparams;


// problem-dependent code activation
#undef USERRESETREGION
#define USERRESETREGION 0

///////////////////////////////////////
//
// disable things that are not really needed because they are debugging type things
//
///////////////////////////////////////

#if(PRODUCTION>=2)
#undef DOEOSDUMP
#define DOEOSDUMP 0
#undef DODISSMEASUREDUMP
#define DODISSMEASUREDUMP 0
#undef DOVPOTDUMP
#define DOVPOTDUMP 0
#undef DODEBUGDUMP // very large, only for speed debug
#define DODEBUGDUMP 0
//#undef DOJETDIAG
//#define DOJETDIAG 0
#endif

#if(PRODUCTION>=3)
#undef DOFLOORDUMP // for accounting for energy from floor and controls if radial fluxes are dumped
#define DOFLOORDUMP 0
#endif

#if(PRODUCTION>=4)

// only needed files for python, not SM
#undef DOMAINDUMP
#undef DORADDUMP
#undef DOIMAGEDUMP

#define DOMAINDUMP 0
#define DORADDUMP 0
#define DOIMAGEDUMP 0

#endif








//////////////////////////////////
// Set or Override with RADIATION settings
//////////////////////////////////









//problem names
#define RADBEAM2D (1) // beam of light in SPC
#define RADTUBE (6) // radiative shock tubes as in Farris et al 09 - assumes Edd.approximation which is currently not handled
#define RADBONDI (7) // like in Fragile's paper (called BONDI in koral)
#define RADPULSE (10) //  radiative blob spreading around
#define RADSHADOW (11)  // radiative shadow
#define RADATM (12) // atmosphere enlighted
#define RADPULSEPLANAR (1000) // like RADPULSE but with scattering
#define RADWAVE (15) // 1d linear rad wave with periodic BC
#define RADPULSE3D (16) // radiative blob spreading around
#define RADDBLSHADOW (17) // radiative shadow with two beams inclined
#define ATMSTATIC (18) // simple hydrostatic atmosphere in SPC
#define RADBEAM2DKS (19) // like RADBEAM2D, just chooses MCOORD KSCOORDS
#define RADBEAMFLAT (24) //  beam of light in Cartesian
#define RADDONUT (25) // 2d radiative Polish donut in KS (called RDONUT in koral.  Similar setup to RADNT.)
#define RADBEAM2DKSVERT (26) // 2d radiative beam in r,theta plane (does more than KS)
#define FLATNESS (27) // flat  (koral: but with non-zero four-force)
#define RADWALL (29) // flat with wall
#define RADNT (30) // emission from midplane
#define RADFLATDISK (31) // emission from flat disk (called FLATDISK in koral.  Very similar to RADNT.)
#define RADCYLBEAM (32) // beam towards the axis in cylindrical (called CYLBEAM in koral.  Somewhat similar to RADFLATDISK but in CYL coords.)
#define RADDOT (33) // radiating dot (Olek changes this while I was testing)
#define RADCYLBEAMCART (40) //  similar to RADCYLBEAM but in cartesian




// TOTRY : optically thin dot in SPC or pulse in SPC.

// TOTRY: Maybe need to avoid bounding if not in PBOUNDLOOP?  Generally true.


//TODO:
#define RADDOTFLAT (41) // similar to RADDOT but in cartesian (well, RADDOT was similar, but still different after Olek changes)
#define RVDONUT (42) // radiative and viscous 
#define RVDONUTIN (43) //  radiative and viscous dougnut inflowing
#define RADNTCYL (44) // emission from midplane in cylindrical

// non-implemented NON-radiative problems in KORAL that are semi-duplicated by some other radiative tests
#define RADINFALL (2) // RADBEAM2D with FLATBACKGROUND=0 is like this
#define DONUT (3) // initboundcode/*fishmon* similar
#define GEODESICINFALL (4) //like RADINFALL but with blobs
#define HDTUBE (8) // as in HARM paper
#define HDTUBE2 (9) // 2D of HDTUBE
#define DONUTOSC (13) // 2d Polish donut oscillating
#define ATMKS (20) //  radial atmosphere infalling in KS
#define DONUTKS (21) // 2d Polish donut in KS (like DONUT)
#define DONUTMKS1 (22) // 2d Polish donut in MKS1 (like DONUT)
#define ATMMKS1 (23) //  radial atmosphere infalling in MKS1 (like ATMKS)
#define BOWSHOCK (28) // bow shock hydro test

// non-implemented radiative problems in KORAL.  Olek says not interesting pre-test versions of other actual tests.
#define RADWAVEBC (14) // 1d linear rad wave imposed on boundary (not setup in koral yet -- looks like time-dep BC for density on left boundary)
#define EDDINFALL (5) // infall with flux from inside


// other additional tests
#define KOMIPROBLEM 50

   // other applications
#define RADCYLJET (51)




// RADDONUT types
#define NODONUT 0
#define DONUTOLEK 1
#define DONUTOHSUGA 2
#define DONUTTHINDISK 3
#define DONUTTHINDISK2 4


////////////////
// other BCTypes beyond those in definit.h (can't overlap numbers from there)
//////////////
#define RADBEAMFLATINFLOW 201
#define RADSHADOWINFLOW 202
#define RADSHADOWINFLOWX2UP 203
#define RADSHADOWINFLOWX2DN 204
#define RADBEAM2DBEAMINFLOW 205
#define RADBEAM2DFLOWINFLOW 206
#define RADATMBEAMINFLOW 207
#define RADWALLINFLOW 208
#define RADBONDIINFLOW 209
#define RADNTBC 210
#define RADCYLBEAMBC 211
#define RADBEAM2DKSVERTBEAMINFLOW 212
#define RADCYLBEAMCARTBC 213
#define HORIZONOUTFLOWSTATIC 214
#define OUTFLOWSTATIC 215
#define RADCYLJETBC 216

#define WALDMONOBC 300


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///// CHOOSE PROBLEM
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

//#define WHICHPROBLEM FLATNESS
//#define WHICHPROBLEM RADBEAMFLAT
//#define WHICHPROBLEM RADPULSE
//#define WHICHPROBLEM RADPULSEPLANAR
//#define WHICHPROBLEM RADPULSE3D
//#define WHICHPROBLEM RADTUBE
//#define WHICHPROBLEM RADSHADOW
//#define WHICHPROBLEM RADDBLSHADOW
//#define WHICHPROBLEM ATMSTATIC
//#define WHICHPROBLEM RADATM
//#define WHICHPROBLEM RADBEAM2D
//#define WHICHPROBLEM RADWALL
//#define WHICHPROBLEM RADWAVE
//#define WHICHPROBLEM RADBONDI
//#define WHICHPROBLEM RADDOT
//#define WHICHPROBLEM RADNT
//#define WHICHPROBLEM RADFLATDISK
#define WHICHPROBLEM RADDONUT
//#define WHICHPROBLEM RADCYLBEAM
//#define WHICHPROBLEM RADBEAM2DKSVERT
//#define WHICHPROBLEM RADCYLBEAMCART
//#define WHICHPROBLEM RADCYLJET








//****************************************//
//****************************************//
//****************************************//
// If set any dimensional constants, should not convert to code units here since conversion not yet defined.  Just set MPERSUN here and convert in init.c
//****************************************//
//****************************************//

// number of solar masses to define units

//#define MPERSUN (1.0)

// odd default choice by koral that gives Mass in cm as 1cm.
#define MPERSUN (6.77255E-1)
//#define MPERSUN (67725.2) // GM/c^2=1cm if gTILDA=1E-10 and cTILDA=1
//#define MPERSUN (6.77255E-11) // with koral's gTILDA=1E-10







//****************************************//
//****************************************//

#if(WHICHPROBLEM==FLATNESS)

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT

#define N1 20
#define N2 20
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif

//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR)

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef FORCESOLVEL
#define FORCESOLVEL 0 // for testing against koral

#define N1 100
#define N2 1 
#define N3 1

#endif

#if(WHICHPROBLEM==RADPULSE3D)

// due to memory per core limited by many variables, can't do 50^3.  Roughly can't do more than 32^3 per core.
#define N1 32
#define N2 32
#define N3 32

#endif

#if(WHICHPROBLEM==RADPULSE || WHICHPROBLEM==RADPULSEPLANAR || WHICHPROBLEM==RADPULSE3D)

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICITSUBCYCLE

#define MCOORD CARTMINKMETRIC2


#endif


//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADBEAMFLAT)

#undef FORCESOLVEL
#define FORCESOLVEL 0 // to compare against koral

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK // works!

//#define N1 20
//#define N2 20
//#define N1 30
//#define N2 30
#define N1 31 // making like problem24 in koral code
#define N2 31 // making like problem24 in koral code
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif



//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADTUBE)

#undef EOMRADTYPE
//#define EOMRADTYPE EOMRADEDD // used by calc_Rij_ff() to set IC so IC use Eddington approximation with Prad=(1/3)Irad (intensity)
#define EOMRADTYPE EOMRADM1CLOSURE

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICITSUBCYCLE
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 800
#define N2 1
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif



//****************************************//
//****************************************//



#if(WHICHPROBLEM==RADSHADOW)

#undef EOMRADTYPE
//#define EOMRADTYPE EOMRADEDD // used by calc_Rij_ff() to set IC so IC use Eddington approximation with Prad=(1/3)Irad (intensity)
#define EOMRADTYPE EOMRADM1CLOSURE

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICITSUBCYCLE
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 100
#define N2 50
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif



//****************************************//
//****************************************//



#if(WHICHPROBLEM==RADDBLSHADOW)

#undef EOMRADTYPE
//#define EOMRADTYPE EOMRADEDD // used by calc_Rij_ff() to set IC so IC use Eddington approximation with Prad=(1/3)Irad (intensity)
#define EOMRADTYPE EOMRADM1CLOSURE

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICITSUBCYCLE
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

//#define N1 120 // code
#define N1 100 // paper
//#define N2 20 // code
#define N2 50 // paper
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif




//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADBEAM2D || WHICHPROBLEM==RADBEAM2DKS)

#undef FORCESOLVEL
#define FORCESOLVEL 0 // doesn't seem to help avoid failures for this test.


#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

// KORALNOTE: Paper says 30x60 for rin-rout and phi=0..pi/2, which is same as 30x30 for rin-rout and phi=0..pi/4 as setup in koral
#define N1 30
#define N2 1
//#define N3 30
#define N3 60 // so like koral paper.

// can choose any spherical polar coordinate system
#if(WHICHPROBLEM==RADBEAM2D)
//#define MCOORD SPCMINKMETRIC
#define MCOORD BLCOORDS // default koral is a=0 BLCOORDS
#elif(WHICHPROBLEM==RADBEAM2DKS)
#define MCOORD KSCOORDS
#endif

#endif


//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADBEAM2DKSVERT)

// below required to avoid runaway energy gains when putting singulary on the grid.  With singularity on grid, can either fail or not in very sensitive way due to machine precision issues right around the singularity.
#undef FORCEGDETPOSITIVE
#define FORCEGDETPOSITIVE 1

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 30
#define N2 30
#define N3 1

//#define MCOORD SPCMINKMETRIC
#define MCOORD KSCOORDS

#endif



//****************************************//
//****************************************//

#if(WHICHPROBLEM==ATMSTATIC)

#undef ANALYTICMEMORY
#define ANALYTICMEMORY 1

//#define DOSTOREPOSITIONDATA 0

//#undef VARTOINTERP
//#define VARTOINTERP PRIMTOINTERP
//#define VARTOINTERP PRIMTOINTERP_GDETFULLVERSION


#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 400
#define N2 1
#define N3 1

// can choose any spherical polar coordinate system
//#define MCOORD SPCMINKMETRIC
//#define MCOORD KSCOORDS
#define MCOORD BLCOORDS

#endif


//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADATM)

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 40
#define N2 1
#define N3 1

// can choose any spherical polar coordinate system with gravity
//#define MCOORD KSCOORDS
#define MCOORD BLCOORDS

#undef MPERSUN
#define MPERSUN (1.0) // So mass=1 as in koral for gTILDE=1.0


#endif

//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADWALL)

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 60
#define N2 20
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif

//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADWAVE)

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#define N1 100
#define N2 1
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif


//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADBONDI)

#undef MPERSUN
#define MPERSUN (3.0)

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

//#define N1 112 // KORALTODO: 512 in paper
#define N1 512 // KORALTODO: 512 in paper
#define N2 1
#define N3 1

// can choose any spherical polar coordinate system with gravity
//#define MCOORD BLCOORDS
#define MCOORD KSCOORDS

#endif



//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADDOT)

//#undef MPERSUN
//#define MPERSUN (1.0/MSUN)

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

// choose odd so DOT is located at center of single cell symmetrically around grid rather than at edge of grid or offset.
#define N1 41
#define N2 41
//#define N3 41 // koral original is 3D, but ok to test in 2D
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif


//****************************************//
//****************************************//

#if(WHICHPROBLEM==RADNT || WHICHPROBLEM==RADFLATDISK || WHICHPROBLEM==RADDONUT || WHICHPROBLEM==RADCYLBEAM || WHICHPROBLEM==RADCYLBEAMCART)

#undef OUTERDEATH
#define OUTERDEATH 1 // do it
//#define OUTERDEATH 0 // don't do it // WALD
#undef OUTERDEATHRADIUS
#define OUTERDEATHRADIUS (3E3)
#undef OUTERDEATHGAMMAMAX
#define OUTERDEATHGAMMAMAX (6.0)
#undef OUTERDEATHGAMMAMAXRAD
#define OUTERDEATHGAMMAMAXRAD (GAMMAMAXRAD)

#undef MPERSUN
#define MPERSUN (10.0)

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK

#undef ANALYTICMEMORY
#define ANALYTICMEMORY 1 // set disk BC using analytical result (at least partially so don't duplicate code.)


#if(WHICHPROBLEM==RADNT)

#define N1 30
#define N2 30
#define N3 1

// can choose any spherical polar coordinate system
//#define MCOORD SPCMINKMETRIC
//#define MCOORD BLCOORDS
#define MCOORD KSCOORDS

#elif(WHICHPROBLEM==RADDONUT)

#undef MUMEAN
#define MUMEAN (MUMEANIONIZED) // ASSUMPTION: fully ionized // CHOICE

#undef WHICHFIT
#define WHICHFIT ISFITNEW // use new fits for radiation donut (fits depend upon other things in global.depmnemonics.rad.h)

#undef DOCOMPTON
#define DOCOMPTON 1 // enable thermal Comptonization

#undef ENSURECONS
#define ENSURECONS 1

#undef DOPERF
#define DOPERF 1 // enable performance enhancements

#undef BORROWENTROPY
#define BORROWENTROPY 1

#undef ENFORCEMHDCONS2RADCONS
#define ENFORCEMHDCONS2RADCONS 1 // assumes not true that pmhd>>prad beyond machine precision

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE // WALD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICITSUBCYCLECHECKSFROMIMPLICIT // least stable result since doesn't use time-advanced prnew to get force since assumes will be doing accurate substeps.
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK // SUPERKORALTODO: Actually doesn't work -- donut heats-up improperly and grows and unsteady.  Even with using time-advanced prnew, leads to problems eventually and noisy overall.

// N1=30 if using log coords from r=1.7 to r=50
// N1=60 if using 1.5*hor - 40 (or 27.8)
// N1=70 if using 1.5*hor - 30 (or 27.8)
#define N1 32
#define N2 16
#define N3 1

   //#define N1 128
//#define N2 64
//#define N3 1

// can choose any spherical polar coordinate system
//#define MCOORD SPCMINKMETRIC
//#define MCOORD BLCOORDS
#define MCOORD KSCOORDS

#undef cTILDA
#define cTILDA (1.0) // like koral
#undef gTILDA
//#define gTILDA (1E-10) // like koral (no longer)
#define gTILDA (1.0)

#undef MPERSUN
#define MPERSUN (10.0*gTILDA) // due to koral fixing MSUNCM, have to do this.

#undef RADIUSMOREDEATH
#define RADIUSMOREDEATH (300.0)



#elif(WHICHPROBLEM==RADFLATDISK)

//#define N1 120 // older koral
#define N1 40 // new koral
#define N2 40
#define N3 1

//#undef WHICHRADSOURCEMETHOD// DEBUG
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT


#define MCOORD SPCMINKMETRIC // i.e. RADFLATDISK

#undef cTILDA
#define cTILDA (1.0) // like koral
#undef gTILDA
#define gTILDA (1E-10) // like koral (no longer)
//#define gTILDA (1.0)

#undef MPERSUN
#define MPERSUN (10.0*gTILDA) // due to koral fixing MSUNCM, have to do this.

//#undef FORCESOLVEL
//#define FORCESOLVEL 1 


#undef ARAD
//#define ARAD (ARAD0*gTILDA*gTILDA) // stupid koral issue with units
#define ARAD (ARAD0) // stupid koral issue with units

#elif(WHICHPROBLEM==RADCYLBEAM)

#define N1 50 // R // 120 for defcoord=UNIFORMCOORDS
#define N2 1 // z
#define N3 30 // \phi

#define MCOORD CYLMINKMETRIC

#elif(WHICHPROBLEM==RADCYLBEAMCART)

#define N1 80
#define N2 80
#define N3 1

#define MCOORD CARTMINKMETRIC2

#endif


#endif





//choose which Komissarov's test number to use (1-9 and 101-106)
#define WHICHKOMI 1


#if(WHICHPROBLEM==KOMIPROBLEM)

#define MCOORD CARTMINKMETRIC2
//#undef EOMTYPE
//#define EOMTYPE EOMGRMHD
////#define EOMTYPE EOMCOLDGRMHD
//#undef EOMRADTYPE
//#define EOMRADTYPE EOMRADNONE // EOMRADM1CLOSURE
//#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK


#define N2 1
#define N3 1

#if(WHICHKOMI==1)
#define N1 40

#elif(WHICHKOMI==2)
#define N1 200

#elif(WHICHKOMI==3)
#define N1 150

#elif(WHICHKOMI==4)
#define N1 150

#elif(WHICHKOMI==5)
#define N1 200

#elif(WHICHKOMI==6)
#define N1 200

#elif(WHICHKOMI==7)
#define N1 400

#elif(WHICHKOMI==8)
#define N1 500

#elif(WHICHKOMI==9)
#define N1 200

#elif(WHICHKOMI>=101 && WHICHKOMI<=109)
#define N1 200 // all use 200, see McKinney (2006) FFDE code paper

#undef EOMTYPE
//#define EOMTYPE EOMGRMHD // can use this as long as bsq/rho large and rho,P constant
#define EOMTYPE EOMFFDE // what these problems really are for

#endif
#endif





#if(WHICHPROBLEM==RADCYLJET)

#undef CONNMACHINEBODY
#define CONNMACHINEBODY 1

//#undef ANALYTICMEMORY
//#define ANALYTICMEMORY 1

#undef ADJUSTFLUX
#define ADJUSTFLUX 1

#undef DOYFL
#define DOYFL 0


#undef OUTERDEATH
#define OUTERDEATH 0 // don't do it
#undef OUTERDEATHRADIUS
#define OUTERDEATHRADIUS (1E3)
#undef OUTERDEATHGAMMAMAX
#define OUTERDEATHGAMMAMAX (6.0)
#undef OUTERDEATHGAMMAMAXRAD
#define OUTERDEATHGAMMAMAXRAD (GAMMAMAXRAD)

#undef MPERSUN
#define MPERSUN (10.0)

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 1

#undef WHICHRADSOURCEMETHOD
//#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
//#define WHICHRADSOURCEMETHOD SOURCEMETHODEXPLICIT
#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICIT
//#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK


//#undef DOCOMPTON
//#define DOCOMPTON 1 // enable thermal Comptonization
//#undef ENSURECONS
//#define ENSURECONS 1
//#undef DOPERF
//#define DOPERF 1 // enable performance enhancements
//#undef BORROWENTROPY
//#define BORROWENTROPY 1
//#undef ENFORCEMHDCONS2RADCONS
//#define ENFORCEMHDCONS2RADCONS 1 // assumes not true that pmhd>>prad beyond machine precision

// RADCYLJET_TYPE==6 needs 2D, rest 1D

#define N1 64 // R // 120 for defcoord=UNIFORMCOORDS
#define N2 32 // z
#define N3 1 // \phi

#define MCOORD CYLMINKMETRIC

#endif






// DEBUG sometimes with the below to check code and geometry issues
#if(0)

#undef WHICHEOM
#define WHICHEOM WITHNOGDET

#undef CONNMACHINEBODY
#define CONNMACHINEBODY 0

#undef DOSTOREPOSITIONDATA
#define DOSTOREPOSITIONDATA 0

#undef STOREFLUXSTATE
#define STOREFLUXSTATE 0

#undef STORESHOCKINDICATOR
#define STORESHOCKINDICATOR 0

#undef VARTOINTERP
#define VARTOINTERP PRIMTOINTERP
//#define VARTOINTERP PRIMTOINTERP_GDETFULLVERSION

#undef RADSHOCKFLAT
#define RADSHOCKFLAT 0

#endif










