//init.h 
//12-feb-06 by alexandre tchekhovskoi
//1-D tests initialization headers
//Define the test number, the grid size for every test number
#ifndef INIT_H

#define INIT_H

//TESTNUMBER:
// -3 -- 1d cartesian expansion on [0, 1]
// -2 -- 2d, vphi == 0 (2d-analogs of 0 and -3)
// -1 -- 2d, vphi != 0 (rotating solution)
//  0 -- 1d cartesian expansion on [-1, 1]
//the rest of the tests, #1 - #9, are as in Liska and Wendroff; 
//test 10 is test 3a from L & W; 
//tests 11 & 12 are the Lax problem at different resolutions.
//test 13 is Jon's smoothed version of the double rarefaction problem for testing the symmetry of code's evolution
//tests 14 (quarter of the space with reflecting walls) & 15 (full space) are the 2D cylindrical Noh problem
//#define TESTNUMBER 666
#define TESTNUMBER 1102

// whether a 1D problem is slab-like or full 2D (will then be tilted at angle such that not moving along diagonal)
#define FULL2D 0


// 1101 full FLUXRECON WENO5BND FLUXCTSTAG gets 3.6kzcps for 128x64

//WENO DEBUG
#define DO_WENO_DEBUG 0  //enables the output of WENO debug information (what it interpolates and which weights it chooses) to wenodebugoutput.out

//define permutation (X, Y, Z) -> (dir1, dir2, dir3)
//can be used to change orientation of test cases
#define DIRX 1
#define DIRY 2
#define DIRZ 3

#define RESOLUTION_FACTOR (1)  //how many times the resolution is additonally to be increased
#define NPROCX 1
#define NPROCY 1
#define NPROCZ 1


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

#define FFDETEST( no ) (no >= 200 && no < 250)

// atch adjusts
//flip the sign of the conservative quantities in the ghost zones on the other side of the singularity to get analytic behaviour of conserved quantities
#undef MAXWELL
#undef TRACKVPOT
#undef EVOLVEWITHVPOT
#undef DOGRIDSECTIONING
#undef MERGEDC2EA2CMETHODEM
#undef MERGEDC2EA2CMETHODMA
#undef MERGEDC2EA2CMETHOD
#undef ACCURATESINCOS

#undef LIMITSOURCES
#undef LIMITDTWITHSOURCETERM
#undef USEGRAVITYDTINDTLIMIT
#undef RESTRICTDTSETTINGINSIDEHORIZON
#undef CHECKONINVERSION
#undef DOYL

#undef DOSTOREPOSITIONDATA
#undef CONNDERTYPE

#undef BONDI_BOUNDARY_SET_PL_PR
#undef WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS
#undef WENO_USE_PRIM_REDUCTION
#undef LIMIT_FLUXC2A_PRIM_CHANGE
#undef COMPDIM
#undef FIELDSTAGMEM
#undef HIGHERORDERMEM
#undef MAXBND
#undef MCOORD
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
#undef DOENODEBUG
#undef DODISS
#undef DOLUMVSR
#undef DODISSVSR
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
#undef COORDSINGFIX
#undef COORDSINGFIXCYL
#undef SINGSMALL
#undef VOLUMEDIFF
#undef MINDT
#undef JONCHECKS
#undef JONCHECKS2
#undef FLOORDIAGS
#undef ANALYTICCONNECTION
#undef ANALYTICSOURCE
#undef OUTFLOWAVOIDBC
#undef FLUXDIMENSPLIT
#undef A2CDIMENSPLIT
#undef DODQMEMORY
#undef BOUNDFLUXRECON
#undef BOUNDARYINTERPADJUST
#undef COMPUTEFRDOT
#undef CALCFARADAYANDCURRENTS
#undef WHICHCURRENTCALC
#undef FARADAYT0
#undef CURRENTST0

#undef OUTFLOW
#undef SYMM
#undef ASYMM
#undef FIXED
#undef POLARAXIS
#undef FIXEDOUTFLOW
#undef NSSURFACE
#undef PERIODIC
#undef BCEXTRAP
#undef CYLAXIS
#undef BCEXTRAP_VEL3

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
#undef HOT2COLD

#define MAXWELL PRIMMAXWELL
#define TRACKVPOT 0 // not on by default
#define EVOLVEWITHVPOT 0 // not on by default
#define DOGRIDSECTIONING 0 // not on by default
#define MERGEDC2EA2CMETHODEM 0
#define MERGEDC2EA2CMETHODMA 0
#define MERGEDC2EA2CMETHOD 0

#define USEGRAVITYDTINDTLIMIT 0
#define RESTRICTDTSETTINGINSIDEHORIZON 0
#define CHECKONINVERSION 1
#define DOYL 0

#define DOSTOREPOSITIONDATA 1 // DEBUG
#define CONNDERTYPE DIFFFINITE // DEBUG

#define ACCURATESINCOS 1
#define LIMITSOURCES 0
#define LIMITDTWITHSOURCETERM 0



// atch adjusts

#define BONDI_BOUNDARY_SET_PL_PR 0  //do analytically set p_l & p_r at the outer boundary for the Bondi problem

#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 0
#define WENO_USE_PRIM_REDUCTION 1


#define LIMIT_FLUXC2A_PRIM_CHANGE 0
#define COMPDIM 3
#define FIELDSTAGMEM 1 // debug
#define HIGHERORDERMEM 1
#define MAXBND 11 // CHANGINGMARK // 13 for WENOBNDPLUSMIN
//#define MAXBND 4
#define PRODUCTION 0
//#define FULLOUTPUT MAXBND
#define FULLOUTPUT 0
#define MAILWHENDONE 0
#define EMAILMESSAGE "Done with Sasha Run #1"
#define EMAILADDRESS "jmckinne@stanford.edu"
#define PERFTEST 0
#define DOAVG 0
#define DOJETDIAG 0
#define DOAVG2 0
#define DODEBUG 1
#define DOENODEBUG 0
#define DODISS 0
#define DOLUMVSR 0
#define DODISSVSR 0
#define DOFIELDLINE 1
#define ROEAVERAGEDWAVESPEED 0
#define ATHENAROE 0
#define STOREWAVESPEEDS 1  //set this and the following one to unity to use the DONOR interpolated states for computing wavespeeds
#define USESTOREDSPEEDSFORFLUX STOREWAVESPEEDS
#define VCHARTYPE VERYLOCALVCHAR
#define PRECISEINVERSION 1
#define WHICHVEL VELREL4
#define WHICHEOM WITHGDET
#define REMOVERESTMASSFROMUU 2
#define RELTYPE RELEOM

#if( FFDETEST(TESTNUMBER) )
#define EOMTYPE EOMFFDE
#else
#define EOMTYPE EOMGRMHD
#endif

#define UTOPRIMTRYAGAIN 0

#if(DODISS || DOLUMVSR || DODISSVSR)
// for diss: testing CHANGINGMARK
#define DOENTROPY DOEVOLVECOMPAREENTROPY
#define WHICHENTROPYEVOLVE EVOLVEFULLENTROPY
#else
// no diss/entropy
#define DOENTROPY DONOENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY
#endif

#if( FFDETEST(TESTNUMBER) || TESTNUMBER==154) // Torus needs initial atmosphere
#define FIXUPAFTERINIT 1
#define FIXUPAFTERRESTART 1
#define CHECKSOLUTION 1
#define GAMMAPERCDIFFMAX (2.0)
#define UPERCDIFFMAX (10.0)
#else
#define FIXUPAFTERINIT 0
#define FIXUPAFTERRESTART 0
#define CHECKSOLUTION 0
#define GAMMAPERCDIFFMAX 1E30
#define UPERCDIFFMAX 1E30
#endif

#define LIMADJUST LIMITERFIXED
#define HYDROLIMADJUSTONLY 0
#define FLUXADJUST FLUXFIXED
#define HYDROFLUXADJUSTONLY 0


#define STEPOVERNEGU NEGDENSITY_NEVERFIXUP
#define STEPOVERNEGRHO NEGDENSITY_NEVERFIXUP
#define STEPOVERNEGRHOU NEGDENSITY_NEVERFIXUP
#define UTOPRIMADJUST UTOPRIMAVG  //controls the behaviour of fixups:  UTOPRIMAVG means fix it up, UTOPRIMSTATIC means do not do it
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED  //controls the behaviour of inversion: does allow the return of solutions with negative densities, etc.
//#define COORDSINGFIX (TESTNUMBER == 153 || TESTNUMBER==154)
#define COORDSINGFIX (0)
#define COORDSINGFIXCYL (MCOORD == CYLMINKMETRIC)
#define SINGSMALL (1E-15)  //make it larger than machine precision so that it is numerically representable when added to pi
#define VOLUMEDIFF 0
#define MINDT 1.e-20 
#define JONCHECKS 1    //SASMARK - do I need this?
#define JONCHECKS2 1   //SASMARK - do I need this?
#define FLOORDIAGS 1
#define ANALYTICCONNECTION ((TESTNUMBER==154)?(0):(1))  //Don't use analytic connection/source terms for spinning BH's -- torus problem works incorrectly with anal. conn./source terms
#define ANALYTICSOURCE ((TESTNUMBER==154)?(0):(1))
#define OUTFLOWAVOIDBC 0
//#define FLUXDIMENSPLIT QUASISTRANG
#define FLUXDIMENSPLIT PERFECTUNSPLIT
//#define FLUXDIMENSPLIT UNSPLIT
//#define A2CDIMENSPLIT QUASISTRANG
#define A2CDIMENSPLIT PERFECTUNSPLIT
//#define A2CDIMENSPLIT UNSPLIT
#define DODQMEMORY 1
#define BOUNDFLUXRECON 0 // choose default to bound only primitives not flux
#define BOUNDARYINTERPADJUST 0  //should be set to zero always
#define COMPUTEFRDOT 0
#define CALCFARADAYANDCURRENTS 1
#define WHICHCURRENTCALC CURRENTCALC1
#define FARADAYT0 1
#define CURRENTST0 1

#define OUTFLOW 0
#define SYMM 1
#define ASYMM 2
#define FIXED 3
#define POLARAXIS 4
#define FIXEDOUTFLOW 5
#define NSSURFACE 6
#define PERIODIC 7
#define BCEXTRAP 100
#define CYLAXIS 101
#define BCEXTRAP_VEL3 102

#define EVOLVECHECKS 0
#define FIXUPZONES ((TESTNUMBER==154||(TESTNUMBER>=200&&TESTNUMBER<250))?(FIXUP1ZONE):(FIXUPNOZONES)) //enable density (rho & u) floors for the torus problem, #154, but disable floors for the rest of the problems
#define HLLBOUNDARY 0
#define FIXUPFLUX 0
#define ZEROOUTFLOWFLUX 0
#define ZEROPOLEFLUX 0
#define BDIRCONT 1
#define HYPERHLL 0
#define HORIZONSUPERFAST 0

#if( FFDETEST(TESTNUMBER) )
#define VARTOINTERP PRIMTOINTERP
#define RESCALEINTERP 0
#define DOEXTRAINTERP 0
#else

#if(1)
#define VARTOINTERP PRIMTOINTERP
//#define VARTOINTERP PRIMTOINTERP_RHOU
//#define VARTOINTERP PRIMTOINTERP_ABSVELREL4
//#define VARTOINTERP PRIMTOINTERP_3VELREL_GAMMAREL
#define RESCALEINTERP 0
#define DOEXTRAINTERP 0
#else // testing
#define VARTOINTERP PRIMTOINTERP_3VELREL_GAMMAREL
#define RESCALEINTERP 1
#define DOEXTRAINTERP 1
#endif

#endif




#define USEAVGPRIMITIVEFORWENOFLAT 1
#define USEPRIMITIVEFROMAVGCONSERVED 0
#define CONTACTINDICATOR 1
#define COMPUTEDRHODP 1
#define SUPERFASTDIVREDUCE 0
#define MINPREFORDER 3
#define SHOCKINDICATOR 1
#define WHICHPARA PARA4

#undef DO_VORTICITY_IMAGE
#define DO_VORTICITY_IMAGE ((FFDETEST(TESTNUMBER))?(0):(1))

#define HOT2COLD 0


//index (in BCtype array) of the outer boundary condition for direction dir
#define BC_TYPE_UP_INDEX( dir ) ((X1UP + (dir - 1) * 2))

//index (in BCtype array) of the inner boundary condition for direction dir
#define BC_TYPE_DN_INDEX( dir ) ((X1DN + (dir - 1) * 2))

#if( defined(DIRX) && defined(DIRY) && defined(DIRZ) )  
//set mnemonics for X1UP and X1DN along the direction of tests; assumes a particular indexing of X1UP, X1DN, X2UP, etc.
#define X_UP   BC_TYPE_UP_INDEX( DIRX )
#define X_DN   BC_TYPE_DN_INDEX( DIRX )
#define Y_UP   BC_TYPE_UP_INDEX( DIRY )
#define Y_DN   BC_TYPE_DN_INDEX( DIRY )
#define Z_UP   BC_TYPE_UP_INDEX( DIRZ )
#define Z_DN   BC_TYPE_DN_INDEX( DIRZ )
#endif

/* size of global arrays */
#if( TESTNUMBER == -3 )
#define NXTEST 32  
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == -2 )
#define NXTEST 32  
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == -1 )
#define NXTEST 32  
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 0 ) // Hubble
#define NXTEST 64  
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 3 ) // 1-D Noh
#define NXTEST 100  
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 666 )// 1d analog of 2d test case 4
//#define NXTEST 100  /// play around with test
#define NXTEST 400  /// play around with test
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 4)
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 10 ) // test3a
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 7 ) // peak
#define NXTEST 800
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 99 ) //asymmetry problem
#define NXTEST 560
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 8 || TESTNUMBER == 9 ) //blast wave or entropy shock wave 
#define NXTEST 400
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 11 ) // Lax's problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 12 ) // high-res Lax
#define NXTEST 400
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 13 )
#define NXTEST 512
#define NYTEST 1
#define NZTEST 1
#elif( TESTNUMBER == 14 ) // 2-D Noh
#define NXTEST 400
#define NYTEST 400
//#define NXTEST 20
//#define NYTEST 80
#define NZTEST 1

#elif( TESTNUMBER == 15 ) // full grid version of 2D Noh
#define NXTEST 64
#define NYTEST 64
#define NZTEST 1

#elif( TESTNUMBER == 16 )//2d test case 3 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 17 )//2d test case 4 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 18 )//2d test case 6 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 19 )//2d test case 12 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 20 )//2d test case 15 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 21 )//2d test case 17 from L&W
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 22 )  //Implosion problem
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 23 )  //Explosion problem
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1

#elif( TESTNUMBER == 24 )  //Smooth problem
#define NXTEST 25
#define NYTEST 25
#define NZTEST 1

#elif( TESTNUMBER == 25 )  //Odd-Even decoupling, blast wave in 2D
#define NXTEST 800
#define NYTEST 10
#define NZTEST 1

#elif( TESTNUMBER == 26 )  //Stationary vortex, Gresho problem
#define NXTEST 40
#define NYTEST 40
#define NZTEST 1

#elif( TESTNUMBER == 27 )  //Moving vortex, Gresho problem
#define NXTEST 160
#define NYTEST 40
#define NZTEST 1

#elif( TESTNUMBER == 28 )  // RT instability problem
#define NXTEST 100  //from Jim Stone's website; should correspond to L&W
#define NYTEST 600
//#define NXTEST 20
//#define NYTEST 80
#define NZTEST 1

#elif( TESTNUMBER == 29 )  // Smooth problem 4.1 from L&W
#define NXTEST 800
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  

#elif( TESTNUMBER == 30 )  //Smooth sound wave problem
#define NXTEST 25
#define NYTEST 25
#define NZTEST 1

#elif( TESTNUMBER == 31 )  //Smooth 1d sound wave problem
#define NXTEST 25
#define NYTEST 1
#define NZTEST 1

#elif( TESTNUMBER == 32 )  //Smooth sound wave problem
#define NXTEST 64
#define NYTEST (NXTEST/2)
#define NZTEST 1

#elif( TESTNUMBER == 33 )  //Smooth density wave problem
#define NXTEST 64
#define NYTEST (NXTEST/2)
#define NZTEST 1

#elif( TESTNUMBER == 49 ) // 1D Caustic problem
#define NXTEST 128
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic) // jon

#elif( TESTNUMBER == 51 )  // Sod's 1D Riemann problem
#define NXTEST 150
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic) // jon



#elif( TESTNUMBER == 52 )  // Sod's 1D Riemann problem
#define NXTEST 1350
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic) // jon


#elif( TESTNUMBER == 101 )  // 1D Riemann problem #1 from RAM paper
#define NXTEST 400
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic) // jon

#elif( TESTNUMBER == 102 )  // 1D Riemann problem #2 from RAM paper
#define NXTEST 400
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 103 )  // 1D Riemann problem #3 from RAM paper
#define NXTEST 400 
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 104 )  // 1D Riemann problem #4 from RAM paper
#define NXTEST 400 
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 105 )  // 1D Riemann problem #5 from RAM paper
#define NXTEST 100  //note lower resolution
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1055 )  // 1D Riemann problem (RPSR (table 3 and figure 5) -- 1-v = 1E-1 to 1e-11  (largest gives \gamma=2.24E5) from Aloy 1999
#define NXTEST 200  //note slightly lower resolution
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1056 )  // harder than Aloy test
#define NXTEST 200  //note slightly lower resolution
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1057 )  // very much harder than Aloy test
#define NXTEST 200  //note slightly lower resolution
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 106 )  // 1D Riemann problem #6.1 from RAM paper (Hard test)
#define NXTEST 400  //note lower resolution
//#define NXTEST 3000
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 107 )  // 1D Riemann problem #3 (reduced) from RAM paper
#define NXTEST 400 
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 151 )  // 2D Riemann shock tube problem from RAM paper, fig. 8
#define NXTEST 400
#define NYTEST 400
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 152 )  // 2D Riemann shock tube problem from RAM paper, fig. 8
#define NXTEST 384
#define NYTEST 1152
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 153 )  // 1D Bondi problem
#define NXTEST 16
#define NYTEST 16
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 154 )  // torus problem
#define NXTEST 16
#define NYTEST 16
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 200 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 201 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 202 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 203 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 204 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 205 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 206 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 207 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 208 )  // torus problem
#define NXTEST 200
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#elif( TESTNUMBER == 1001 )  // 1D Mignone mild blast wave
#define NXTEST 800 
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1002 )  // 1D Mignone strong blast wave
#define NXTEST 800
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1003 )  // RJ95A
#define NXTEST 400
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1100 )  // 2D Current sheet
// Original Athena is 256x256, Athena web version is 100x100 as here
// Athena seems to generate more islands
#define NXTEST 100
#define NYTEST 100
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1101 )  // 2D Field loop
#define NXTEST 128
#define NYTEST 64
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1102 )  // Circ Alfven
#define NXTEST 256 // 256 is Athena standard resolution
#define NYTEST (NXTEST/2) // gives square cells (ala Athena website) if Lx=2Ly as designed in init.c for this problem
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 1103 )  // 1D 2D MHD wave tests
// if make 2D, then presently just does 1D slab in X-direction
#define NXTEST 256
#define NYTEST (NXTEST/2) // not necessary to have this NXTEST/2 -- can be anything
//#define NYTEST 1 // not necessary to have this NXTEST/2 -- can be anything
#define NZTEST 1
#undef FULL2D
#define FULL2D 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU  2

#elif( TESTNUMBER == 667 )  // uniform density distribution in cyl.coords -- not stationary?
#define NXTEST 384
#define NYTEST 1
#define NZTEST 1
#undef REMOVERESTMASSFROMUU
#define REMOVERESTMASSFROMUU 2  //redefine to use rest mass so that can use the UTOPRIM2D inversion (can use it because the problem is relativistic)

#else
#define NXTEST 100   //all other tests: 1D with this default number of grid cells
#define NYTEST 1
#define NZTEST 1
#endif


#if( TESTNUMBER == 28 )  // RT instability problem
#define MCOORD UNIGRAVITY
#elif( TESTNUMBER == 152  || TESTNUMBER == 667) //slab jet or constant density distrib. in cyl. coord.
#define MCOORD CYLMINKMETRIC
#elif( TESTNUMBER == 153 || TESTNUMBER == 154) 
#define MCOORD KSCOORDS
#else
#define MCOORD CARTMINKMETRIC
#endif




#if( defined( NXTEST ) && defined( NYTEST ) && defined( NZTEST ) && !defined( N1 ) && !defined( N2 ) && !defined( N3 ) )
//computational grid size has not been defined, so define it


//decrease the resolution per processor to get the same combined (MPI) resolution;
//and enlarge the resolution (used for convergence testing) <-- done only for non-degenerate dimensions
#define RESOLUTION_FACTOR_X (( NXTEST == 1 ) ? ( 1 ) : ( RESOLUTION_FACTOR ))
#define RESOLUTION_FACTOR_Y (( NYTEST == 1 ) ? ( 1 ) : ( RESOLUTION_FACTOR ))
#define RESOLUTION_FACTOR_Z (( NZTEST == 1 ) ? ( 1 ) : ( RESOLUTION_FACTOR ))


#define NXTEST_PER_PROC (NXTEST * RESOLUTION_FACTOR_X /NPROCX)
#define NYTEST_PER_PROC (NYTEST * RESOLUTION_FACTOR_Y /NPROCY)
#define NZTEST_PER_PROC (NZTEST * RESOLUTION_FACTOR_Z /NPROCZ)

#if( ((RESOLUTION_FACTOR_X * NXTEST)%NPROCX!=0)||((RESOLUTION_FACTOR_Y * NYTEST)%NPROCY!=0)||((RESOLUTION_FACTOR_Z * NZTEST)%NPROCZ!=0) )
#error "Need to have integral Ntotal/Nproc"
#endif


//expands into N?TEST where '?' is X, Y, or Z is such that DIR?_OF_TESTS == dir.
#define NDIR_TEST_FUNC( dir ) (NXTEST_PER_PROC * ( DIRX == dir ) + NYTEST_PER_PROC * ( DIRY == dir ) + NZTEST_PER_PROC * ( DIRZ == dir ))

//NOTE: additionally increase resolution by RESOLUTION_FACTOR
#define N1  NDIR_TEST_FUNC( 1 )
#define N2  NDIR_TEST_FUNC( 2 )
#define N3  NDIR_TEST_FUNC( 3 )

#else
#error "init.h: unsupported test case"
#endif


#define WHICH_INITVEL VEL3

struct Ccoordparams {
  FTYPE timescalefactor, rho0, u0, Omega0, vz0, tstart, L0;
};


#endif

#undef FLUXDUMP
#define FLUXDUMP 0

#undef ASYMDIAGCHECK
#define ASYMDIAGCHECK 0











