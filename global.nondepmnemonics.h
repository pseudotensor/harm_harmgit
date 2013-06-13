// all things here don't depend on anything else, just names for numbers or purely functional macros


#include "metric.h"
#include "coord.h"

// RK-related macros
// Doesn't even depend upon N?, N?M, or N?BND, etc.
#include "global.stepch.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various physics and model setup parameters that are macros either for performance reasons or since no need to change them at runtime.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






// define how to access symmetric matrices with size 4x4 without redundant elements

// 1) regexp: gcov\[\([_\>a-zA-Z0-9+-\ ()]+\)\]\[\([_\>a-zA-Z0-9+-\ ()]+\)\] -> gcov[GIND(\1,\2)]
// Then revert: gcov\[GIND( *NDIM *, *NDIM *)\] -> gcov[SYMMATRIXNDIM] since otherwise will be 1 larger than required
// also need to catch: localgcon gcon gcovinfunc gcovtovks gcovbhks gcovmcoord gcovmid tmpgcov gcovprim tmpgcon glgen ghgen
// 2) Then need to replace any multi-D pointer arg type with simple arg type:
// FTYPE (**localgcov)[NDIM] -> FTYPE **localgcov
// 3) Then need to get all global variables with [NDIM][NDIM]
// 4) Also get *gcov type things:  (\*gcov)\[NDIM\] -> *gcov
//    also need to catch: gcon gcovinfunc gcovinfuncprim gcovprim gconprim gcovmcoord gcovptr gconptr
// 5) If making assignment TO something using GIND(), then must control loops.
//    E.g. tetrad.c:tetr_func_frommetric(): newgcov[GIND(jj,kk)] +=
//   Find some maybe by doing: grep -e "+ \{0,\}=" *.c *.h | grep GIND
// Compared new and old codes by doing:
// for fil in `ls -d *.c *.h`; do echo $fil; diff -bBdpy -W 600 --suppress-common-lines $fil $1/$fil >>newcodediff.txt ; done ; less newcodediff.txt
// Check for name mangling: grep "gcov" newcodediff.txt | grep "gcon" | less

// below are non-conditional ways of getting same result as:
//#define GCOVI(i,j) (i>=j) ? i : j
//#define GCOVJ(i,j) (i>=j) ? j : i
#define GCOVI(i,j) ((i>=j)*(i-j) + j)
#define GCOVJ(i,j) ((i>=j)*(j-i) + i)
#define GIND(i,j) GCOVJ(i,j)*4 + GCOVI(i,j) - MAX(GCOVJ(i,j),0) - MAX(GCOVJ(i,j)-1,0) - MAX(GCOVJ(i,j)-2,0)
// must multiply assignments by the below so don't duplicate sums
#define GINDASSIGNFACTOR(i,j) (1.0*(i>=j))
//#if(PRODUCTION==0)
//#define GINDASSIGNMAC(name,i,j) (i>=j ? name[GIND(i,j)] : SHOULDNOTREACHHEREEVERBUGYOUHAVE())
//#else
//#define GINDASSIGNMAC(name,i,j) name[GIND(i,j)]
//#endif



// define how many values in table
#define NRANC 64



//////////////////////
//
// nmenomics associated with definit.h
//
//
//
////////////////////////
#define NUMJETS 2
#define INNERJET 0
#define OUTERJET 1


#define QUASISTRANG 0
#define UNSPLIT 1
#define PERFECTUNSPLIT 2


// set enerregion types
#define NUMENERREGIONS 8

#define NULLENERREGIONS -2 // indicates avoid processes that operate on this region
#define ALLENERREGIONS -1 // to indicate not just one region
#define GLOBALENERREGION 0 // standard computational grid
// "outside" horizon means r>=r_+
#define OUTSIDEHORIZONENERREGION 1 // "outside" horizon: i>ihorizon
// TRUEGLOBAL???? seeks to include all cells that could be active
// That is, horizon not likely to be at standard boundary edge yet need roughly nbnd cells inside horizon always to evolve there properly, so at least have those active but no less than i=0 itself.  If max(0) part used, then probably issue with number of cells inside boundary, however.
#define TRUEGLOBALENERREGION 2 // max(0,ihorizon-nbnd)
#define TRUEGLOBALWITHBNDENERREGION 3 // outside horizon with boundary cells: max(-nbnd,ihorizon-2*nbnd)
// ACTIVE??? is directly used by GRIDSECTIONING
#define ACTIVEREGION 4  // only active computational cells
#define ACTIVEWITHBNDREGION 5  // ACTIVEREGION+required boundary cells
// Jet regions:
#define INNERJETREGION 6 // lower-\theta jet region
#define OUTERJETREGION 7 // upper-\theta jet region


// for VCHARTYPE
#define VERYLOCALVCHAR 0
#define LOCALVCHAR 1
#define GLOBALVCHAR 2


// number of extra things in state that ucon_calc computes
// right now used for rel4vel storing gamma and qsq to avoid computing multiple times
#define NUMOTHERSTATERESULTS 2
#define OTHERGAMMA 0
#define OTHERQSQ 1


#define ERRORCODEBELOWCLEANFINISH 1000

// maximum char size of filenames or directories or commands
#define MAXFILENAME 200
// for longer file names
#define MAXFILENAMELONG 2000

// return code used to indicate failure of dump_gen() was file not found instead of any other error
#define FILENOTFOUND 2490834

// for WHICHVEL
// 0: 4-velocity (leads to ambiguous u^t +- discr part)
// 1: 3-velocity (unambiguous u^t but interpolation is not constrained to be a good 3-velocity)
// 2: relative 4-velocity (unambiguous u^t and any interpolation gives good value)
#define VEL4 0
#define VEL3 1
#define VELREL4 2

// for WHICHEOM
#define WITHGDET 0
#define WITHNOGDET 1
#define WITHSINSQ 2

// This stays naturally, simply consistent with how code evolves conserved quantities.
// for WHICHEOM
#define NOGDETRHO 0
#define NOGDETU0 0
#define NOGDETU1 1 // U1 and U2 are only reasonable choices in SPC
#define NOGDETU2 1
#define NOGDETU3 0
#define NOGDETB1 0
#define NOGDETB2 0
#define NOGDETB3 0
#define NOGDETENTROPY 0
#define NOGDETYL 0
#define NOGDETYNU 0




// defines how one forms the EM stress-energy tensor
#define GENMAXWELL 0
#define PRIMMAXWELL 1


// for RELTYPE
#define RELEOM 0
#define NONRELEOM 1 // NOT FINISHED // NOT RIGHT
// whether relativistic or nonrelativistic EOMs (speed of light limitation)

// for EOMTYPE
// 0 = FF(D)E force-free electrodynamics
// 1 = cold GRMHD
// 2 = entropy conservation version of GRMHD
// 3 = GRMHD
// for force-free, must turn off:
// ok now, but effectively setup already the below 2 lines implicitly
// global.h : FIXUPAFTERINIT, FIXUPAFTERRESTART,CHECKSOLUTION,LIMADJUST,FLUXADJUST
// global.h FIXUPZONES->FIXUPNOZONES
#define EOMFFDE 0
#define EOMCOLDGRMHD 1
#define EOMENTROPYGRMHD 2
#define EOMGRMHD 3





// macros for defining which fluxcalc method to use in flux.c
#define ORIGINALFLUXCALC 0
#define NEWFLUXCALC 1



////////////////////////////
//
// Choices for cooling function
//
////////////////////////////

#define NOCOOLING 0
#define COOLGAMMIETHINDISK 1
#define COOLEOSGENERAL 2
#define COOLREBECCATHINDISK 3
#define COOLUSER 10 // user-defined cooling function


////////////////////////////
//
// Choices for EOS
//
////////////////////////////


#define NUMEOSS 5 // number of EOS types

#define COLDEOS 0
#define IDEALGAS 1
#define MIGNONE 2
#define GRBPWF99 3
#define KAZFULL 4



// mnenomics
// for DOENTROPY
#define DONOENTROPY  0
#define DOEVOLVEENTROPY 1 // generic activation of entropy variable in conservation laws, etc.

#define DONOYL 0
#define DOEVOLVEYL 1

#define DONOYNU 0
#define DOEVOLVEYNU 1


// for WHICHENTROPYEVOLVE
#define EVOLVENOENTROPY 0
#define EVOLVESIMPLEENTROPY 1 // should be used with DOENTROPY==DOEVOLVECOMPAREENTOPY
#define EVOLVEFULLENTROPY 2 // should only be used with DOENTROPY==DOEVOLVEDIRECTENTROPY or DOENTROPY==DOEVOLVECOMPAREENTOPY



// defines how Utoprimgen is used
#define EVOLVEUTOPRIM 0
#define OTHERUTOPRIM 1


// defines data return types for primtoU() and primtoflux()
#define UEVOLVE 0
#define UDIAG 1
#define UNOTHING 2
#define UENTROPY 3 // implicit UNOTHING but cons UU is overwritten by cons entropy (see primtoflux() in phys.c as used by utoprim() in utoprim.c)




// for LIMADJUST
// 0: use fixed limiter
// 1: use limiter based upon b^2/rho
// 2: use limiter based upon b^2/u
// 3: use limiter based upon both b^2/rho or b^2/u
#define LIMITERFIXED 0
#define LIMITERBSQORHO 1
#define LIMITERBSQOU 2
#define LIMITERBSQORHOANDU 3

// for FLUXADJUST
#define FLUXFIXED 0 // (see get_bsqflags() in fixup.c)
#define FLUXBSQORHO 1
#define FLUXBSQOU 2
#define FLUXBSQORHOANDU 3


// for UTOPRIMFAILRETURNTYPE  --  controls the behaviour of inversion: does allow the return of solutions with negative densities, etc.
#define UTOPRIMRETURNNOTADJUSTED 0
#define UTOPRIMRETURNADJUSTED 1



// for UTOPRIMADJUST  -- controls the behaviour of fixups:  UTOPRIMAVG means fix it up, UTOPRIMSTATIC means do not do it
// 0=just use static solution
// 1=use average surrounding solution, and if no good surrounding solution use the normal observer velocity with static densities
#define UTOPRIMSTATIC 0
#define UTOPRIMAVG 1


// interpolation function used by init.readdata.c and other code
#define LINEARTYPE 0
#define LOGTYPE 1
#define QUADRATICTYPE 2 // includes limiters to ovoid overshoots and new extremums within 3 point domain


// for MODIFYEMFORVPOT
#define MODIFYEMF 0
#define MODIFYVPOT 1

////////////////////////////////////////////////
//
// Some things related to higher-order interpolations
//
////////////////////////////////////////////////
// used to choose which method interpline.c uses
#define NUMENOINTERPTYPES 12

#define NONENOINTERPTYPE 0
#define ENOINTERPTYPE 1 // ce2
#define ENOINTERPTYPE4EMF 2 // f2corn
#define ENOFLUXRECONTYPE 3
#define ENOFLUXRECONTYPEGHOSTACTIVE 4 // used at t=0 where *construct* Uavg in ghost+active region and for dofluxreconevolvepointfield==0 for EMFs
#define ENOFLUXSPLITTYPE 5
#define ENOAVG2CENTTYPE 6
#define ENOCENT2AVGTYPE 7
#define ENOFLUXAVG1TYPE 8
#define ENOFLUXAVG2TYPE 9
#define ENOFLUXAVG3TYPE 10
#define ENOQUASIFIELDFLUXRECONTYPE 11


//quantities to interp
#define ENOSOURCETERM 0
#define ENOCONSERVED  1
#define ENOPRIMITIVE  2
#define ENOFLUX       3
#define ENOMAFLUX     4
#define ENOSMOOTHFLUX 5
#define ENOSMOOTHCONSERVED  6


// for DOENOFLUX:
#define NOENOFLUX 0
#define ENOFLUXRECON 1
#define ENOFLUXSPLIT 2
#define ENOFINITEVOLUME 3
// 0: no ENO flux reconstruction
// 1: reconstruct F for finite difference rep. of U
// 2 : flux splitting (not done yet)
// 3: reconstruct dU for finite volume rep. of U

// df and monoindicator sizes for interpline.c and reconstruct.c
// 0,1,2,3 for paraline and 0,4 for SMONO and 0 (dP) for WENO
#define NUMDFS 5
// currently don't need DFCENT2APART, so this is why NUMDFS is 5 and not 6
// within calculation NUMDFS is checked so no access to that array element

#define DFONESIDED 0
#define DFCENT 1
#define DFMONO 2
#define DF2OFMONO 3
#define DF2OFONESIDED 4
#define DFCENT2APART 5

#define NUMMONOINDICATORS 3

#define MONOINDTYPE 0 // -1,0,1 for rough, ambiguous, monotonic
#define MONOLEFTSET 1 // whether left interface (or central for a2c/c2a) was set to MONO value
#define MONORIGHTSET 2 // whether right interface was set to MONO value



////////
//
//  possible settings for pass_1d_line() in interpline.c
//
#define WEIGHT_CALC 1
#define RECON_CALC 2
#define ALL_CALC (WEIGHT_CALC | RECON_CALC)  //compute everything


// defines types of high order interpolations
#define CVT_A2C 0
#define CVT_C2A 1
#define CVT_C2L 2
#define CVT_C2R 3
#define CVT_C2DER1 4
#define CVT_C2DER2 5
#define CVT_C2DER3 6
#define CVT_C2DER4 7



#define CVT_C2E CVT_C2L //use the same number as CVT_C2L because does not add a new reconstruction type

//0 -- don't do weighs minimization
// -- do 1st version of weights minimization
#define NOSPLITA2C 0
#define MINIMIZE_ALL_WEIGHTS 1
#define ENERGY_CONTROLS_ALL_WEIGHTS 2 // as below
#define ENERGY_IS_ALL_WEIGHTS 3 // means T^i_i for flux and T^t_t for energy
#define MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS 4 // means lock \rho u^i , T^i_i and T^i_t all together to help with relativistic flows to maintain consistency
#define MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD 5 // means lock \rho u^i , T^i_i and T^i_t all together to help with relativistic flows to maintain consistency
#define GENFUN_IS_ALL_WEIGHTS 5
#define CONSTANT_ALL_WEIGHTS 6 // forces MONO and equal weights for all


////////
//
//  MONOINTERP 
//
#define NOMONOINTERP 0
#define JMONOINTERP 1
#define SMONOINTERP 2

// defines types of input to flux_point2avg()
#define ISMAONLY 0
#define ISEMONLY 1
#define ISMAANDEM 2



#define DISSSIMPLEINVCO 0
#define DISSFULLINVCO 1
#define DISSENTROPYCO 2
#define DISSSIMPLEINVCONOMAX 3 
#define DISSFULLINVCONOMAX 4
#define DISSENTROPYCONOMAX 5
#define DISSSIMPLEINVLAB1 6
#define DISSFULLINVLAB1 7
#define DISSENTROPYLAB1 8
#define DISSSIMPLEINVLAB1NOMAX 9
#define DISSFULLINVLAB1NOMAX 10
#define DISSENTROPYLAB1NOMAX 11
#define DISSSIMPLEINVLAB2 12
#define DISSFULLINVLAB2 13
#define DISSENTROPYLAB2 14
#define DISSSIMPLEINVLAB2NOMAX 15
#define DISSFULLINVLAB2NOMAX 16
#define DISSENTROPYLAB2NOMAX 17
// failure indicator:
#define DISSFAILUREINV 18

// totals:
#define NUMDISSVERSIONS 18 // for consistent restart output, this number can't change
#define NUMDISSFUNPOS (NUMDISSVERSIONS+1)  // includes failure indicator


// EMF loop related things
#define NUMPOS4EMF 3
#define LEFT4EMF 0 // left or -
#define RIGHT4EMF 1 // right or +
#define CENT4EMF 2  // centered position
#define START4EMF LEFT4EMF // which is 0
#define END4EMF CENT4EMF // which is last value


// see interp_loop_set() in initbase.c
#define NUMFLUXLOOPNUMBERS 10
#define FIDEL 0
#define FJDEL 1
#define FKDEL 2
#define FFACE 3
#define FIS 4
#define FIE 5
#define FJS 6
#define FJE 7
#define FKS 8
#define FKE 9


// number of inversion quantities to report when inversion fails if CHECKONINVERSION = 1
// See utoprim_jon.c:check_on_inversion()
#define NUMINVPROPERTY 13

// for WHICHCURRENTCALC
// 0: original time is on edge and spatial on edge, but spatials are different locations.  old time.
// 1: all centered in space and all time, present time (best)
// 2: like 0, but spatially centered (i.e. old time)
#define CURRENTCALC0 0
#define CURRENTCALC1 1
#define CURRENTCALC2 2



#define CURRENTPRECALCTYPES 5

#define CURTYPET 0
#define CURTYPEX 1
#define CURTYPEY 2
#define CURTYPEZ 3
#define CURTYPEFARADAY 4


// whether and which type of fixups to be used
#define FIXUP1ZONE 0
#define FIXUPALLZONES 1
#define FIXUPNOZONES 2



/* mnemonics for flux method (Riemann solver) */
// ordered from most diffusive to least diffusive, so can back track
// 0 should be reasonable most diffusive
#define LAXFFLUX 0
#define HLLFLUX 1
#define FORCEFLUX 2
#define MUSTAFLUX 3 // not yet working
#define HLLLAXF1FLUX 4


// DIVB constraint method
#define FLUXCTHLL 0
#define FLUXCTTOTH 1
#define FLUXCD 2
#define ATHENA1 3
#define ATHENA2 4
#define FLUXCTSTAG 5
/* these are different ways of calculating the EMFs */
//#define FLUXB FLUXCTTOTH
// 0: HLL
// 1: FLUXCT TOTH version (toth 2000 eq. 25)
// 2: FLUXCD TOTH version (toth 2000 eq. 31)
// 3: Athena type eq 39
// 4: Athena type eq 48
// 5: Jon's staggered grid #1

//#define UTOPRIMVERSION 6
// 0: original gammie 5D method
#define UTOPRIM5D1 0
// 1: ldz method
#define UTOPRIMLDZ 1
// 2: SCN 2D method
#define UTOPRIM2D 2
// 3: SCN 1D method
#define UTOPRIM1D 3
// 4: SCN 1D OPTIMIZED method -- not sure if identical to 3 otherwise
#define UTOPRIM1DOPT 4
// 5: SCN 1D final and optimized
#define UTOPRIM1DFINAL 5
// 6: SCN 2D final and optimized and recommended by Scott
#define UTOPRIM2DFINAL 6
// 7: SCN 5D final -- bit less accurate compared to 1D and 2D
#define UTOPRIM5D2 7
// 8: Jon 1D/2D final version -- can handle non-rel problems
#define UTOPRIMJONNONRELCOMPAT 8
// 100: use 5D, but compare with ldz in runtime
#define UTOPRIMCOMPARE 100


/* mnemonics for slope limiter */
// ordered from most diffusive to least diffusive, so can start high and go down if needed
// 0 should be reasonble most diffusive, highest should be least diffusive

#define DONOR   0

// POSINTERPS:
#define VANL    1
#define MINM    2
#define MC      3
#define PARA    4
#define PARAFLAT 5
#define MCSTEEP 6 // uses 3-point limiter, but other features of PARAFLAT
#define CSSLOPE      7 // not tested/compared against others

// assume here and beyond all higher numbers are using WENO or ENO
#define WENO3 8
#define WENO4 9
#define WENO5  10
#define WENO6  11
#define WENO7  12
#define WENO8  13
#define WENO9  14
#define ENO3 15
#define ENO5 16
#define WENO5FLAT 17
#define WENO5BND 18
#define WENO5BNDPLUSMIN 19

#define PARALINE 20

#define FIRSTWENO WENO3
#define LASTWENO WENO5BNDPLUSMIN

#define FIRSTINTERPLINE WENO3
#define LASTINTERPLINE PARALINE

#define NUMPOSINTERPS LASTINTERPLINE
#define NUMNEGINTERPS 4
// 1+ for DONOR
#define NUMINTERPS (1 + NUMPOSINTERPS + NUMNEGINTERPS)


// defines limiters that are WENO/ENO
#define WENOINTERPTYPE(lim) (lim>=FIRSTWENO && lim<=LASTWENO)

// defines limiters that are WENO/ENO
#define WENOBNDPINTERPTYPE(lim) (lim==WENO5BND || lim==WENO5BNDPLUSMIN)

// defines which limiters are for interpline.c (rest are for interppoint.c)
#define LINEINTERPTYPE(lim) (lim>=FIRSTINTERPLINE && lim<=LASTINTERPLINE)

// negative versions for testing only
#define NLIM    -1 // no limiter
#define NLIMCENT    -2 // no limiter
#define NLIMUP    -3 // no limiter
#define NLIMDOWN    -4 // no limiter


// see orders_set() in initbase.c
#define MAXSPACEORDER 15 // maximum number of points in stencil
//#define MAXSPACESHIFT ((MAXSPACEORDER-1)/2) // center point for symmetric stencil

// for timing
#define STARTTIME 0
#define CHECKTIME 1
#define SPEEDTIME 2
#define STOPTIME 3
#define REPORTTIME 4
#define DIAGSTARTTIME 5
#define DIAGSTOPTIME 6
#define INITSTARTTIME 7
#define INITSTOPTIME 8


#define MAXTIMEORDER 4

//#define TIMEORDER 3
// order of algorithm in time from 1 to 4.
// TIMEORDER: 1 : single step (Euler method -- error term is 2nd order for smooth flows)
// TIMEORDER: 2 : 2 steps in halfs (midpoint method -- error term is 3rd order for smooth flows)
// TIMEORDER: 3 : 4 steps (classic RK3 method -- error term is 4th order for smooth flows)
// TIMEORDER: 4 : 4 steps (classic RK4 method -- error term is 5th order for smooth flows)

// tetrad.c stuff:
#define METRICTETRAD 0
#define NONMETRICTETRIC 1


//////////////////////////////////
//
// which variable to interpolate
//
/////////////////////////////////

#define PRIMTOINTERP -1
#define PRIMTOINTERP_JONRESCALED1 0
#define CONSTOINTERP 1
#define PRIMTOINTERPLGDEN 2
#define PRIMTOINTERP_LGDEN_RHOU 3
#define PRIMTOINTERP_RHOU 4
#define PRIMTOINTERP_VSQ 5
#define PRIMTOINTERP_3VEL_GAMMA 6
#define PRIMTOINTERP_RHOV_GAMMA 7
#define PRIMTOINTERP_VELREL4SQ 8
#define PRIMTOINTERP_3VELREL_GAMMAREL 9
#define PRIMTOINTERP_RAMESH1 10
#define PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP 11
#define PRIMTOINTERP_GDETFULLVERSION 12


#define NOFIELDRESCALE -1
#define NOSPECIALFIELD 0
#define PULSARFIELD 1
#define PULSARFIELD2 2
#define PULSARFIELD3 3
#define GDETVERSION 4
#define GDETFULLVERSION 5



#define WENO_REDUCE_TYPE_DEFAULT 0 
#define WENO_REDUCE_TYPE_PPM 1

// definition of minmod operator
#define MINMODB(a,b) ( (fabs(a)<fabs(b)) ? (a) : (b) )
#define MINMOD(a,b) ( ((a)*(b)<=0) ? 0.0 : MINMODB(a,b) )
//#define MINMOD3( x, y, z )   ( 0.25 * (sign(x) + sign(y)) * (sign(x) + sign(z)) * MIN( MIN(fabs(x), fabs(y)), fabs(z)) )    

#define MINMODB(a,b) ( (fabs(a)<fabs(b)) ? (a) : (b) )
#define MINMODGEN(extremeallow,a,b) ( (!extremeallow && (a)*(b)<=0) ? 0.0 : MINMODB(a,b) )


#define REMOVEFROMNPR 0
#define RESTORENPR 1


///////////////////////////////
//
// parabolic interpolation stuff
//
////////////////////////////////

#define PARA1 0 // old
#define PARA2 1 // works
#define PARA3 2 // broken
#define PARA4 3 // latest
#define PARAJON 4 // Created to do well with high \sigma monopole


// GODMARK: wth NUMREC had problems with large run (jetnewnoenv,jetnew on sauron)
// GODMARK: NUMREC not working right now after trying to make accurate
//#define CONNDERTYPE GAMMIEDERIVATIVE

#define DIFFGAMMIE 0 // use infinitesimal differences of analytical metric
#define DIFFNUMREC 1 // use advanced (but presently broken) uniformly accurate Numerical Recipies numerical derivatives
#define DIFFFINITE 2 // use previously defined gridded values of metric to compute finite differences (fastest)


#define INTERPPOINTTYPE 0
#define INTERPLINETYPE 1



//////////////////////////////////////
//
// PURE mnemonics
//
///////////////////////////////////
#define NUMBOUNDTYPES 9
//
#define BOUNDPRIMTYPE 0
#define BOUNDPRIMSIMPLETYPE 1
#define BOUNDPSTAGTYPE 2
#define BOUNDPSTAGSIMPLETYPE 3
#define BOUNDINTTYPE 4 // always simple
#define BOUNDFLUXTYPE 5
#define BOUNDFLUXSIMPLETYPE 6
#define BOUNDVPOTTYPE 7
#define BOUNDVPOTSIMPLETYPE 8


// ispstag:
#define BOUNDPRIMLOC 0
#define BOUNDPSTAGLOC 1


// -------------> r
// |         3    
// |        1-0   
// |         2    
// v            
// theta      
// and likewise for 4,5 (4=out,5=in)
// directions:
#define X1UP    0
#define X1DN    1
#define X2UP    2
#define X2DN    3
#define X3UP    4
#define X3DN    5

#define NUMUPDOWN 2
#define POINTUP 0
#define POINTDOWN 1

#define NUMLEFTRIGHT 2
#define ISLEFT 0
#define ISRIGHT 1
#define ISMIDDLE 2 // just macro, not used to access memory space and so why NUMLEFTRIGHT is still 2

// used by doflux[] to see if flux surface on grid or not (each CPU)
#define FLUXNOTONGRID -100


// direction (-1,+1) for a given direction as defined above
#define DIRSIGN(dir) (1-2*((dir)%2 ))
// dimension=1,2,3 for given direction defined above
#define DIMEN(dir) (1+(dir)/2)

#define DIRFROMDIMEN(dimen,dirsign) (( (dirsign==-1) + (dimen - 1) * 2))

// direction (0,1) for a given dir=X1DN,etc.
#define POINTFROMDIR(dir) (DIRSIGN(dir)==-1 ? POINTDOWN : POINTUP)

// long double constants
# define M_El           2.7182818284590452353602874713526625L  /* e */
# define M_LOG2El       1.4426950408889634073599246810018922L  /* log_2 e */
# define M_LOG10El      0.4342944819032518276511289189166051L  /* log_10 e */
# define M_LN2l         0.6931471805599453094172321214581766L  /* log_e 2 */
# define M_LN10l        2.3025850929940456840179914546843642L  /* log_e 10 */

# define M_PIl          3.1415926535897932384626433832795029L  /* pi */
# define M_PI_2l        1.5707963267948966192313216916397514L  /* pi/2 */
# define M_PI_4l        0.7853981633974483096156608458198757L  /* pi/4 */
# define M_1_PIl        0.3183098861837906715377675267450287L  /* 1/pi */
# define M_2_PIl        0.6366197723675813430755350534900574L  /* 2/pi */
# define M_2_SQRTPIl    1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
# define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
# define M_SQRT1_2l     0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
# define SIXTH          0.1666666666666666666666666666666666L  /* 1/6 */
# define FOURTHIRD      1.3333333333333333333333333333333333L  /* 4/3 */
# define THIRD          0.3333333333333333333333333333333333L  /* 1/3 */
# define ONE            1.0000000000000000000000000000000000L
# define PTFIVE         0.5L
# define TWO            2.0L
# define ONEPT25        1.25L
# define THREE          3.0L
# define SIX            6.0L
# define EIGHT          8.0L

#ifdef WIN32
# define M_PI           3.1415926535897932384626433832795029L  /* pi */
#endif



#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SIGNSINGLE(a) ( ((a) <0.) ? -1. : 1. )
// rounds to definite integer (round() returns double and so isn't useful as an integer value)
#define ROUND2INT(x) ((int)((x)>0.0 ? (x)+0.5 : (x)-0.5))
#define ROUND2LONGLONGINT(x) ((long long int)((x)>0.0 ? (x)+0.5 : (x)-0.5))


// restart macro stuff
// or use DODISS, etc. that are 0 or non-zero
// assume any non-zero will work in code
#define DONOTACCESSMEMORY 0


#define PROGRADERISCO 0
#define RETROGRADERISCO 1

#define NUMTSCALES 4
// number of times scales to watch failure rates at
#define ALLTS 0 // full cumulative
#define ENERTS 1 // cumulative each dump_ener (over all grid)
#define IMAGETS 2 // cumulative each image dump (full grid)
#define DEBUGTS 3 // debug dump time scale (full grid)



//STEPOVERNEGXXX: possible modes of stepping over occurences of negative densities: controls when to revert to inversion from an average conserved quantity
//For how these occurences are reported, see fixup.c: DOCOUNTUNEG, etc.
#define NEGDENSITY_NEVERFIXUP -1
#define NEGDENSITY_ALWAYSFIXUP 0
#define NEGDENSITY_FIXONFULLSTEP 1




// see failfloorcount counter
#define COUNTNOTHING -1
#define NUMFAILFLOORFLAGS 15
//  mnemonics
#define COUNTUTOPRIMFAILCONV 0 // if failed to converge
#define COUNTFLOORACT 1 // if floor activated
#define COUNTLIMITGAMMAACT 2 // if Gamma limiter activated
#define COUNTINFLOWACT 3 // if inflow check activated
#define COUNTUTOPRIMFAILRHONEG 4
#define COUNTUTOPRIMFAILUNEG 5
#define COUNTUTOPRIMFAILRHOUNEG 6
#define COUNTGAMMAPERC 7 // see fixup_checksolution()
#define COUNTUPERC 8 // see fixup_checksolution()
#define COUNTENTROPY 9
#define COUNTCOLD 10
#define COUNTEOSLOOKUPFAIL 11
#define COUNTBOUND1 12 // see bounds.tools.c (used when boundary code actually affects active zone values)
#define COUNTBOUND2 13
#define COUNTONESTEP 14


// below 3 used to indicate when eos lookup failure shouldn't report failure since (e.g.) was not at a particular grid location
#define AVOIDI -100
#define AVOIDJ -100
#define AVOIDK -100

// failure codes for utoprim failures
// NOTE: PFLAGTYPE is probably "char" so can't use value of pflag beyond -127..127
#define NANPFLAG -100 // bad pflag
#define UTOPRIMFAILFIXEDONESTEP -6
#define UTOPRIMFAILFIXEDBOUND2 -5
#define UTOPRIMFAILFIXEDBOUND1 -4
#define UTOPRIMFAILFIXEDCOLD -3
#define UTOPRIMFAILFIXEDENTROPY -2
#define UTOPRIMFAILFIXEDUTOPRIM -1
#define UTOPRIMNOFAIL 0
#define UTOPRIMFAILCONV 1
#define UTOPRIMFAILCONVW     2
#define UTOPRIMFAILCONVUTSQ  3
#define UTOPRIMFAILCONVGUESSUTSQ 4
#define UTOPRIMFAILCONVUTSQVERYBAD 5
#define UTOPRIMFAILCONVBADINVERTCOMPARE 6
#define UTOPRIMFAILNANGUESS 7
#define UTOPRIMFAILNANRESULT 8
#define UTOPRIMFAILRHONEG 9
#define UTOPRIMFAILUNEG 10
#define UTOPRIMFAILRHOUNEG 11
#define UTOPRIMFAILGAMMAPERC 12
#define UTOPRIMFAILUPERC 13
#define UTOPRIMFAILU2AVG1 14
#define UTOPRIMFAILU2AVG2 15
#define UTOPRIMFAILU2AVG1FROMCOLD 16
#define UTOPRIMFAILU2AVG2FROMCOLD 17
#define UTOPRIMFAILFAKEVALUE 18
#define UTOPRIMFAILCONVRET   50

#define IFUTOPRIMFAILSOFTRHORELATED(pflag) (pflag==UTOPRIMFAILRHONEG || pflag==UTOPRIMFAILRHOUNEG)

#define IFUTOPRIMFAILSOFTNOTRHORELATED(pflag) (pflag==UTOPRIMFAILUNEG || pflag==UTOPRIMFAILGAMMAPERC  || pflag==UTOPRIMFAILUPERC  || pflag==UTOPRIMFAILU2AVG1  || pflag==UTOPRIMFAILU2AVG2 || pflag==UTOPRIMFAILU2AVG1FROMCOLD  || pflag==UTOPRIMFAILU2AVG2FROMCOLD)

#define IFUTOPRIMFAILSOFT(pflag) (IFUTOPRIMFAILSOFTRHORELATED(pflag)||IFUTOPRIMFAILSOFTNOTRHORELATED(pflag))

#define IFUTOPRIMFAILFIXED(pflag) (pflag==UTOPRIMFAILFIXEDCOLD || pflag==UTOPRIMFAILFIXEDENTROPY || pflag==UTOPRIMFAILFIXEDUTOPRIM)

#define IFUTOPRIMNOFAILORFIXED(pflag) (pflag==UTOPRIMFAILFIXEDCOLD || pflag==UTOPRIMFAILFIXEDENTROPY || pflag==UTOPRIMFAILFIXEDUTOPRIM || pflag==UTOPRIMNOFAIL)

#define IFUTOPRIMNOFAIL(pflag) (pflag==UTOPRIMNOFAIL)

#define IFUTOPRIMFAIL(pflag) (pflag>UTOPRIMNOFAIL)


/* failure modes */
#define FAIL_UTOPRIM_NEG        1
#define FAILSTR01 "UTOPRIM_NEG"
#define FAIL_UTOPRIM_TEST       2
#define FAILSTR02 "UTOPRIM_TEST"
#define FAIL_VCHAR_DISCR        3
#define FAILSTR03 "VCHAR_DISCR"
#define FAIL_COEFF_NEG          4
#define FAILSTR04 "COEFF_NEG"
#define FAIL_COEFF_SUP          5
#define FAILSTR05 "COEFF_SUP"
#define FAIL_UTCALC_DISCR       6
#define FAILSTR06 "UTCALC_DISCR"
#define FAIL_LDZ                7
#define FAILSTR07 "FAIL_LDZ"
#define FAIL_BCFIX              8
#define FAILSTR08 "FAIL_BCFIX"
#define FAIL_VSQ_NEG            9
#define FAILSTR09 "FAIL_VSQ_NEG"

/* mnemonics for primitive vars; conserved vars */
#define VARNOTDEFINED -100
#define RHO     0
#define UU      1
#define U1      2
#define U2      3
#define U3      4
#define B1      5
#define B2      6
#define B3      7

// for STORESHOCKINDICATOR
#define NUMSHOCKPLS 3
#define SHOCKPLDIR1 0
#define SHOCKPLDIR2 SHOCKPLDIR1+1
#define SHOCKPLDIR3 SHOCKPLDIR2+1

// for STORESHOCKINDICATOR for temp storage
// assumes no more than NDIM in list!
#define SHOCKPLSTOREPTOT 0
#define SHOCKPLSTOREVEL1 1
#define SHOCKPLSTOREVEL2 SHOCKPLSTOREVEL1+1
#define SHOCKPLSTOREVEL3 SHOCKPLSTOREVEL2+1


// primitive type
#define CENTEREDPRIM 0
#define STAGGEREDPRIM 1


#define NDIM    4               /* number of total dimensions.  Never changes */
#define SYMMATRIXNDIM 10 // total number of independent elements of a symmetric matrix

// flag failures/problems for correction/check in fixup
#define NUMPFLAGS (5)
// the below needs to be bounded since one CPU doesn't know if the other failed, and neighbor failure determines nature of how failure is treated
// also, bounded values at real boundaries need to identify if bad copy
#define FLAGUTOPRIMFAIL 0 // changes behavior of fixup()
// the below flags are done after bound_prim, and can be determined at any time, so just come after bound.
#define FLAGREALLIM 1 // value of limiter to be used
#define FLAGBSQORHO 2 // set when B^2/RHO > BSQORHOLIMIT ; currently changes  behavior of slope_lim
#define FLAGBSQOU 3 // set when B^2/u > BSQOULIMIT
#define FLAGREALFLUX 4 // type of flux to use




#define NUMSOURCES 3
// these get ADDED UP, not independently treated
// number of source terms.  Currently includes: 0) geometry, 1) radiative cooling, 2) radiative heating
#define GEOMSOURCE 0 // SHOULD ALWAYS BE 0 !
#define RADSOURCE 1
#define RADSOURCE2 2





// max number of terms in stress tensor (for term-level flux diagnostic)
#define NUMFLUXTERMS (7)


#define ORDERDEBUG 3
//#define NUMENODEBUGS (NPR*(3+ORDERDEBUG*3 + 1 + 2 + ORDERDEBUG*2 + 1))
// p p_l p_r  order*3 per point reduce
// NPR*3 order*3*NPR   NPR
// Uavg Upoint  order*2 per point reduce
// NPR*2 order*2*NPR NPR

//#define NUMENODEBUGS 21
// see email

// short switches:
//1)       SMONO (0,1)
//2)       WENO5 (0,1)
//3)       WENO3 (0,1)
//4)       -> dP/P (0,1)
//5)       limit c2e/c2a/a2c correction (0,1) through checking the change of the quantity being interpolated
//6)       limit c2e/c2a/a2c correction (0,1) through checking the change of primitives

#define NUMENODEBUGS 6


// allow for pk[0] and pk[1]
#define MAXDTSTAGES 2
// maximum number of allowed temporal integration stages


// time period between dumps for various types of dumps
// NO, use "DUMPTYPE" names now
//#define NUMDTDS 11
//#define DTDUMP 0
//#define DTAVG 1
//#define DTENER 2
//#define DTIMAGE 3
//#define DTDEBUG 4
//#define DTDISS 5
//#define DTFIELDLINE 6
//#define DTFLUX 7
//#define DTOTHER 8
//#define DTEOS 9
//#define DTVPOT 10





/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3


/* mnemonics for centering of grid functions */
// GODMARK: is there a way to pick and choose the dimension and number of grid positions?
/* number of positions on grid for grid 
   functions */
#define NPG 8
#define NOWHERE -1 // tells a function that not necessarily requesting value at a standard grid location (unusual) -- this stores no memroy, so doesn't increase NPG
#define CENT    0
#define FACE1   1
#define FACE2   2
#define FACE3   3
#define CORN1   4 // corner in 2-3 plane
#define CORN2   5 // corner in 1-3 plane
#define CORN3   6 // corner in 1-2 plane
#define CORNT   7 // true corner: full 3D corner (only required for 3D)

// used for primgridpos[]
#define NUMPRIMGRIDPOS 2
#define CENTGRID 0
#define STAGGRID 1

/* mnemonics for diagnostic calls */
#define INIT_OUT        0
#define DUMP_OUT        1
#define IMAGE_OUT       1
#define LOG_OUT         1
#define FINAL_OUT       2
#define FUTURE_OUT      3

// below should not be >=0
#define DOINGFUTUREOUT -100


#define SURFACETOTAL 0
#define VOLUMETOTAL 1

#define CONSTYPE 0
#define SURFACETYPE 1
#define CUMULATIVETYPE 2
#define CONSJETINNERTYPE 3
#define CONSJETOUTERTYPE 4
#define CUMULATIVETYPE2 5
#define CONSTYPE2 6
#define CUMULATIVETYPE3 7

#define NOTHINGHEAD -1
#define WRITEHEAD 0
#define READHEAD 1

#define TIMESERIESAREAMAP 0
#define FINALTDUMPAREAMAP 1

#define WRITEFILE 0
#define READFILE 1

#define ENERFNAME "ener.out"
#define GENERFNAME "gener.out"



// these dump types also control period of output
// Period can be controlled for non-spatial dumps such as ENER outputs, in which case dump.c doesn't have to be setup for that type of "DUMPTYPE"
#define NUMDUMPTYPES 19 // number of dump types listed below

#define IMAGEDUMPTYPE 0
#define RESTARTDUMPTYPE 1
#define RESTARTUPPERPOLEDUMPTYPE 2
#define RESTARTMETRICDUMPTYPE 3
#define MAINDUMPTYPE 4
#define GRIDDUMPTYPE 5
#define AVG1DUMPTYPE 6
#define AVG2DUMPTYPE 7 // used when needing AVG2DUMPTYPE to avoid too large a file size for avgdump
#define DEBUGDUMPTYPE 8
#define FIELDLINEDUMPTYPE 9
#define ENODEBUGDUMPTYPE 10
#define DISSDUMPTYPE 11
#define OTHERDUMPTYPE 12
#define FLUXDUMPTYPE 13
#define EOSDUMPTYPE 14
#define VPOTDUMPTYPE 15
#define FAILFLOORDUDUMPTYPE 16
#define ENERDUMPTYPE 17
#define FAKEDUMPTYPE 18

#define MYDUMPNAMELIST {"IMAGEDUMPTYPE","RESTARTDUMPTYPE","RESTARTUPPERPOLEDUMPTYPE","RESTARTMETRICDUMPTYPE","MAINDUMPTYPE","GRIDDUMPTYPE","AVG1DUMPTYPE","AVG2DUMPTYPE","DEBUGDUMPTYPE","FIELDLINEDUMPTYPE","ENODEBUGDUMPTYPE","DISSDUMPTYPE","OTHERDUMPTYPE","FLUXDUMPTYPE","EOSDUMPTYPE","VPOTDUMPTYPE","FAILFLOORDUDUMPTYPE","ENERDUMPTYPE","FAKEDUMPTYPE"}


#define NUMIMAGEPARMS 3

#define ORIGIN 0
#define LMIN 1
#define LMAX 2











