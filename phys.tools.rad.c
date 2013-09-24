#include "decs.h"

//////////////////////////////
//
// BEGIN: f2c stuff for Ramesh's solver code in fortran
//
//////////////////////////////

// f2c prototype
#include "f2c.h"


#ifdef KR_headers
//double floor();
integer i_dnnt(x) doublereal *x;
#else
                  //#undef abs
                  //#include "math.h"
integer i_dnnt(doublereal *x)
#endif
{
return( (*x)>=0 ?
	floor(*x + .5) : -floor(.5 - *x) );
}

#include "testfpp.P"
// not linking with libf2c since don't want that dependence and conversion doesn't need it since the original code was simple

int get_rameshsolution(int whichcall, int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, FTYPE errorabs, FTYPE errorabsbestexternal, int iters, int totaliters, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent);

int get_rameshsolution_wrapper(int whichcall, int eomtype, FTYPE errorabs, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *dUother, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE (*dUcompeng)[NPR], FTYPE (*dUcompent)[NPR], struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent);

//////////////////////////////////////////////
//
// END: f2c stuff for Ramesh's solver code in fortran
//
//////////////////////////////////////////////




////////////////////////////////
//
// Local options (used to be in global.nondepmnemonics.rad.h)
//
///////////////////////////////

#define COURRADEXPLICIT (0.1) // Effective Courant-like factor for stiff explicit radiation source term.  Required to not only avoid failure of explicit scheme, but also that explicit scheme is really accurate compared to implicit.  E.g., near \tau\sim 1, explicit won't fail with RADPULSEPLANAR but will not give same results as implicit.  So only use explicit if really in optically thin regime.


////////////////////////////////////
#define IMPTRYCONVHIGHTAU (NUMEPSILON*5.0)  // for used implicit solver

// Funny, even 1E-5 does ok with torus, no worse at Erf~SMALL instances.  Also, does ~3 iterations, but not any faster than using 1E-12 with ~6 iterations.
#define IMPTRYCONV (1.e-8) // works generally to avoid high iterations
#define IMPTRYCONV2 (1.e-8) // use for high gamma or high tau cases
// error for comparing to sum over all absolute errors
#define IMPTRYCONVABS ((FTYPE)(NDIM+2)*IMPTRYCONV)

// too allowing to allow 1E-4 error since often solution is nuts at even errors>1E-8
#define IMPALLOWCONV (MAX(IMPTRYCONV,1.e-3))
#define IMPALLOWCONV2 IMPALLOWCONV
#define IMPALLOWCONVABS ((FTYPE)(NDIM+2)*IMPALLOWCONV)

// what tolerance to use for saying can switch to entropy when u_g is suggested to be bad for energy
#define IMPOKCONV (MAX(IMPTRYCONV,1E-10))
#define IMPOKCONVABS ((FTYPE)(NDIM+2)*IMPOKCONV)

// tolerance above which say energy solution is probably bad even if not very large error.  These have tended (or nearly 100%) to be cases where actual solution has u_g<0 but harm gets error u_g>0 and error not too large.
#define IMPBADENERGY (MIN(IMPALLOWCONV,1E-7))


// whether to abort even the backup if error is not reducing.
#define ABORTBACKUPIFNOERRORREDUCE 1
#define IMPTRYCONVALT (MAX(1E-8,IMPTRYCONVABS)) // say that if error isn't reducing, ok to abort with this error.   Only time saver, but realistic about likelihood of getting smaller error.

// tolerance above which to continue to try damping
#define IMPTRYDAMPCONV (5.0*IMPTRYCONVABS)

////////////////////////////////



// IMPLICIT SOLVER TOLERANCES or DERIVATIVE SIZES
 // for used implicit solver (needs to be chosen more generally.  KORALTODO: 1E-8 too small in general).  Could start out with higher, and allow current checks to avoid inversion failure.
//#define IMPEPS (1.e-8)
// use large, and it'll go smaller if no inversion, but can't start out with too small since then Jac will have diag() terms =0
// KORALTODO: For difficult iterations, there can be solution but Jacobian is too rough and jump around alot in primitive space for small changes in U.  Should really modify IMPEPS in such cases when pr changes alot for such changes in U.
// roughly (NUMEPSILON)**(1/3) as in NR5.7 on numerical derivatives
#if(REALTYPE==FLOATTYPE)
#define IMPEPSLARGE (1E-4) // on small side
#define IMPEPSSMALL (1E-4) // on small side
#define ERRORFORIMPEPSSMALL (1E-5)
#elif(REALTYPE==DOUBLETYPE)
#define IMPEPSLARGE (1E-6)
#define IMPEPSSMALL (1E-6)
//#define IMPEPSSMALL (1E-8)
#define ERRORFORIMPEPSSMALL (1E-9)
#elif(REALTYPE==LONGDOUBLETYPE)
#define IMPEPSLARGE (1E-8)
#define IMPEPSSMALL (1E-10)
#define ERRORFORIMPEPSSMALL (1E-9)
#endif

#define IMPEPS IMPEPSSMALL

// maximum EPS for getting Jacobian
#define MAXIMPEPS (0.3)

// maximum number of times to (typically) increase EPS in getting Jacobian for implicit scheme.  This might generally override MAXIMPEPS.
#define MAXJACITER (10)



//#define IMPMAXITER (15) // for used implicit solver // For others
#define IMPMAXITER (100) // for used implicit solver
#define IMPMAXITER2 (100) // for used implicit solver

#define IMPMAXITERMEDIUM (40)

#define IMPMAXITERQUICK (20)

#define IMPMINABSERROR (1E-100) // minimum absolute error (or value) below which don't treat as bad error and just avoid 4-force.  Otherwise will "fail" implicit solver even if impossible to reach smaller relative error due to absolute machine precision issues.

// 1 : normalize radiation error by only radiation thermal energy
// 2 : normalize radiation error by max(radiation,gas) thermal energy
// 3 : normalize using radiation URAD0 but also fnorm from actual f
// 4 : URAD0, fnorm, and UU
// normalize error.  Can't expect radiation to be relatively accurate to itself if UU>>URAD0 due to G between them
// 3 is safest, but more expensive than 4.  4 should be fine for real systems.
#define IMPLICITERRORNORM 4


#define MAXSUBCYCLES (2000) // for explicit sub-cycles when doing reversion

// if tries more than this number of sub-cycles, just fail and assume no 4-force since probably due to no actual solution even for implicit scheme due to sitting at radiative failure (e.g. gamma->gammamax or Erf->ERADLIMIT)
#define MAXSUBCYCLESFAIL (MAXSUBCYCLES*100)


#define MAXF1TRIES 50 // 50 might sound like alot, but Jacobian does 4 to 5 inversions each iteration, and that amount is only typically needed for very first iteration.
  // goes to f1iter=10 for RADPULSE KAPPAES=1E3 case.  Might want to scale maximum iterations with \tau, although doubling of damping means exponential w.r.t. f1iter, so probably 50 is certainly enough since 2^(-50) is beyond machine precision for doubles.

#define RADDAMPDELTA (0.5) // choose, but best to choose 1/Integer so no machine precision issues.
#define RADDAMPUNDELTA (1.0/RADDAMPDELTA)


#define TAUFAILLIMIT (2.0/3.0) // at what \tau below which to assume "failure1" in u2p_rad() means should be moving at gammamax rather than not moving.

// whether to revert to sub-cycle explicit if implicit fails.  Only alternative is die.
#define IMPLICITREVERTEXPLICIT 0 // problem.  Not a good idea. -- should try implicit again, starting out with more damping.

// like SAFE for normal dt step, don't allow explicit substepping to change dt too fast to avoid instabilities.
#define MAXEXPLICITSUBSTEPCHANGE 1.e-2

// 0 : tau suppression
// 1 : space-time merged
// 2 : all space merged but separate from time
// 3 : full split
// 4 : split even mhd and rad
#define TAUSUPPRESS 0 // makes physical sense, but might be wrong in some limits (i.e. Gd can still be large relative), but try to account for lambda as well.
#define SPACETIMESUBSPLITNONE 1 // DON'T USE! (gets inversion failures because overly aggressive)
#define SPACETIMESUBSPLITTIME 2 // probably not ok -- need to split off mhd and rad
#define SPACETIMESUBSPLITALL 3 // probably not ok -- need to split off mhd and rad
#define SPACETIMESUBSPLITSUPERALL 4 // OK TO USE sometimes, but not always
#define SPACETIMESUBSPLITMHDRAD 5 // KINDA OK TO USE (generates noise in velocity where E is very small, but maybe ok since sub-dominant and don't care about velocity where E is small.  E evolves fine, but Rtx eventually shows issues.)
#define SPACETIMESUBSPLITTIMEMHDRAD 6 // OK TO USE sometimes (works fine and no noise in velocity because split it off.  Might have trouble in multiple dimensions if sub-dominant momentum dimension requires implicit step -- but general diffusive would suggest unlikely.  But, not efficient in optically thick regime once radiation velocity is non-trivial in magnitude)

#define WHICHSPACETIMESUBSPLIT TAUSUPPRESS // only tausuppress works in general.

//SPACETIMESUBSPLITTIMEMHDRAD // TAUSUPPRESS
// SPACETIMESUBSPLITTIMEMHDRAD  //  SPACETIMESUBSPLITMHDRAD // SPACETIMESUBSPLITSUPERALL // TAUSUPPRESS



// whether to use dUriemann and dUgeom or other dU's sitting in dUother for radiation update
#define USEDUINRADUPDATE 1



#define GAMMASMALLLIMIT (1.0-1E-10) // at what point above which assume gamma^2=1.0


// whether to choose Jon or Olek way of handling u2p_rad inversion failures
#define JONCHOICE 0
#define OLEKCHOICE 1

#define CASECHOICE JONCHOICE // choose
//#define CASECHOICE OLEKCHOICE // choose

#define TOZAMOFRAME 0 // reduce to ZAMO gammarel=1 frame (e.g. in non-GR that would be grid frame or v=0 frame or gammarel=1).
#define TOFLUIDFRAME 1 // reduce to using fluid frame (probably more reasonable in general).
#define TOOPACITYDEPENDENTFRAME 2

#define M1REDUCE TOOPACITYDEPENDENTFRAME // choose

// KORALTODO: The below need to be chosen intelligently

// below this \tau, no source term applied.
// KORALTODO: Need to fix implicit solver so avoids dU-del in fluid if no radiatoin-fluid interaction, else overestimates effect and inversion failures occur.
#define MINTAUSOURCE (NUMEPSILON)



///////////////////////////////
//
// START SOME LOCAL OPTIONS
//
///////////////////////////////

// whether to do subjac iter-dependent solver.
// 0 : normal full 4D method
// 1 : Invert sub Jacobian method
#define DOSUBJAC 0
#define ENERGYSIMPLEFACTOR (NDIM) // NDIM=4 times simpler than full step


#define ITERMODENORMAL 0
#define ITERMODESTAGES 1
#define ITERMODECOLD 2


////////////
// for ITERSTAGES
#define BEGINMOMSTEPS0 1
#define ENDMOMSTEPS0 2

#define BEGINENERGYSTEPS0 3
#define ENDENERGYSTEPS0 13

#define BEGINFULLSTEPS0 14
#define ENDFULLSTEPS0 (IMPMAXITER*2)

#define BEGINNORMALSTEPS0 BEGINFULLSTEPS0
////////////



// number of prior iterations to see if error has dropped enough to seem relevant
#define NUMPRIORERRORSITER0 7
#define NUMPRIORERRORS 5
#define PRIORERRORCHANGEREQUIRED (0.5) // damping is included below when used

// whether to store steps for primitive and so debug max iteration cases
#define DEBUGMAXITER 1

// 0: show primitive
// 1: show u^\mu and urad^\mu
#define DEBUGMAXITERVELOCITY 1

#define DEBUGLEVELIMPSOLVER 3 // which debugfail>=# to use for some common debug stuff
//#define DEBUGLEVELIMPSOLVER 2 // which debugfail>=# to use for some common debug stuff

#define DEBUGLEVELIMPSOLVERMORE 3 // which debugfail>=# to use for some common debug stuff
//#define DEBUGLEVELIMPSOLVERMORE 2 // which debugfail>=# to use for some common debug stuff

// how many holds on u_g to apply while stepping velocity.
#define RAMESHFIXEARLYSTEPS (DOSUBJAC==1 ? 0 : 3) // 3 if used is default

// whether to apply Jon's hold on u_g or rho from going negative
#define JONHOLDPOS 0

// number of times allowed to hold u_g as positive
#define NUMHOLDTIMES 6


// stop iterating energy if pbenergy[UU]<0.5*pbentropy[UU] consistently starting aafter below number of iterations and lasting for 2nd below number of iterations
#define RAMESHSTOPENERGYIFTOOOFTENBELOWENTROPY 3

//#define SWITCHTOENTROPYIFCHANGESTOENTROPY (*implicitferr==QTYUMHD ? 0 : 1)
#define SWITCHTOENTROPYIFCHANGESTOENTROPY (0)


// whether to try again with 10* higher u_g for entropy case since important as backup.
// Sometimes initial guess for u_g is too low, and when too low very hard for NR to recover.
#define TRYENTROPYHARDER 1

#define TRYENERGYHARDER 1

// whether to use ramesh solver as backup
#define USERAMESH 0 // FUCK: Can only use if QTYPMHD method used.  URAD will mess up if Erf<<ug.

// error below which to feed best guess into Ramesh solver
#define TRYHARDERFEEDGUESSTOL (1E-6)

// error below which will use entropy as guess for energy if entropy didn't hard fail.
#define ERRORTOUSEENTROPYFORENERGYGUESS (1E-4)


// whether to get lowest error solution instead of final one.
#define GETBEST 0


// need to do final check since get f1 and then do step.
#define DOFINALCHECK 1


// below 1 if reporting cases when MAXITER reached, but allowd error so not otherwise normally reported.
#define REPORTMAXITERALLOWED (PRODUCTION==0)

// whether to ensure rho and u_g in Jacobian calculation difference do not cross over 0 and stay on same side as origin point.
#define FORCEJDIFFNOCROSS 0

// whether to check pp-ppp
// 1: directly check post pp-ppp relative error and see if changes by LOCALPREIMPCONVX
// 2: directly check if any changes to pp during Newton step are bigger than DIFFXLIMIT.
#define POSTNEWTONCONVCHECK 0

// below which sum of all primitives is taken as no interesting change.
#define DIFFXLIMIT (10.0*NUMEPSILON)

#define LOCALPREIMPCONVX (10.0*NUMEPSILON)
#define LOCALPREIMPCONVXABS ((FTYPE)(NDIM+2)*LOCALPREIMPCONVX)


// number of iterations by which to check (1st) whether after some number of times (2nd) error rose instead of reduced.
#define NUMNOERRORREDUCE0 (5+BEGINNORMALSTEPS)
#define NUMNOERRORREDUCE 5


// whether to use EOMDONOTHING if error is good enough.
// 1: always avoid external inversion (so no longer can do cold MHD, but cold MHD in \tau\gtrsim 1 places is very bad).  Or avoid energy switching to entropy, which also is bad.
#define SWITCHTODONOTHING 0

  // whether to change damp factor during this instance.
#define CHANGEDAMPFACTOR 1
#define NUMDAMPATTEMPTS 3

#define NUMDAMPATTEMPTSQUICK 3


// factor by which error jumps as indication that u_g stepped to was very bad choice.
#define FACTORBADJUMPERROR (1.0E2)


// 0 : old Jon  method
// 2 : Jon's paper draft method
#define WHICHU2PRAD 0

// during implicit solver, don't limit gamma so much as normally.  Otherwise, solution may not be found and solver struggles and leads to high errors and iterations.  If limit gammarad but not gammafluid, then gammafluid can be too high.  If limit both, no solutions can be found.   So just limit afterwards for now.
//#define GAMMAMAXRADIMPLICITSOLVER (1E5)
#define GAMMAMAXRADIMPLICITSOLVER (GAMMAMAXRAD) // for radiation, seek actual solution with this limit.  Solver can find solutions, while harder when limiting gamma_{gas} for some reason.


// whether to avoid computing entropy during iterations if not needed
#define ENTROPYOPT 0

///////////////////////////////
//
// END SOME LOCAL OPTIONS
//
///////////////////////////////






//////// implicit stuff
static int koral_source_rad_implicit(int *eomtype, FTYPE *pb, FTYPE *pf, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother ,FTYPE (*dUcomp)[NPR]);

static int koral_source_rad_implicit_mode(int modprim, int havebackup, int didentropyalready, int *eomtype, int whichcap, int itermode, int trueimpmaxiter, int truenumdampattempts, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pb, FTYPE *uub, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR], FTYPE *errorabs, FTYPE errorabsbestexternal, int *iters, int *f1iters);

static int f_implicit(int iter, int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pp, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *f, FTYPE *fnorm, int *goexplicit);

static int get_implicit_iJ_old(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, FTYPE impepsjac, int iter, FTYPE errorabs, FTYPE *dimfactU, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *fpl1, FTYPE *fpl1norm, FTYPE (*iJ)[NPR]);

static int get_implicit_iJ(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, FTYPE impepsjac, int iter, FTYPE errorabs, FTYPE *dimfactU, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NPR]);
static int inverse_33matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM]);
static int inverse_11matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM]);


static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE convabs, FTYPE realdt, int dimtypef, int eomtype, int itermode, FTYPE fracenergy, FTYPE dissmeasure, FTYPE *dimfactU, FTYPE *pp, FTYPE *piin, FTYPE *f1, FTYPE *f1norm, FTYPE *f1report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom, FTYPE *errorabs);

static int f_error_check_wrap(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, FTYPE *f1, FTYPE *f1norm, FTYPE *f1report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom);


static int Utoprimgen_failwrapper(int doradonly, int *radinvmod, int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int whichcap, int evolvetype, int inputtype,FTYPE *U,  struct of_state *qptr, struct of_geom *ptrgeom, FTYPE dissmeasure, FTYPE *pr, struct of_newtonstats *newtonstats);

static void define_method(int iter, int *eomtype, int itermode, FTYPE fracenergy, FTYPE dissmeasure, int *implicititer, int *implicitferr, int *BEGINMOMSTEPS, int *ENDMOMSTEPS, int *BEGINENERGYSTEPS, int *ENDENERGYSTEPS, int *BEGINFULLSTEPS, int *ENDFULLSTEPS, int *BEGINNORMALSTEPS);
static void get_refUs(int *numdims, int *startjac, int *endjac, int *implicititer, int *implicitferr, int *irefU, int *iotherU, int *erefU, int *eotherU, int *signgd2, int *signgd4, int *signgd6, int *signgd7);


// debug stuff
static void showdebuglist(int debugiter, FTYPE (*pppreholdlist)[NPR],FTYPE (*ppposholdlist)[NPR],FTYPE (*f1reportlist)[NDIM],FTYPE (*f1list)[NDIM],FTYPE *errorabsf1list, int *realiterlist, FTYPE *jac00list);
int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, FTYPE errorabs, FTYPE errorabsbestexternal, int iters, int iterstotal, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);

// explicit stuff
static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gdpl, FTYPE chi, FTYPE *Gdplabs, struct of_geom *ptrgeom, FTYPE *dtsub);
static void koral_source_dtsub_rad_calc(int method, FTYPE *pr, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gdpl, struct of_geom *ptrgeom, FTYPE *dtsub);
static int source_explicit(int whichsc, int whichradsourcemethod, int methoddtsub,int *eomtype,
                           void (*sourcefunc)(int method, FTYPE *pr, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gpl, struct of_geom *ptrgeom, FTYPE *dtsub),
                           FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother, FTYPE (*dUcomp)[NPR]);

// RAD inversion stuff
static int get_m1closure_gammarel2_old(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn);
static int get_m1closure_gammarel2(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn);

static int get_m1closure_gammarel2_cold_old(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel);
static int get_m1closure_gammarel2_cold(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel);

static int get_m1closure_Erf(struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE gammarel2, FTYPE *Erfreturn);

static int get_m1closure_urfconrel_old(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad);
static int get_m1closure_urfconrel(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad);
static int get_m1closure_urfconrel_olek(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad);

static int opacity_interpolated_urfconrel(FTYPE tautotmax, FTYPE *pp,struct of_geom *ptrgeom,FTYPE *Avcon, FTYPE Erf,FTYPE gammarel2,FTYPE *Erfnew, FTYPE *urfconrel);

// general stuff
static FTYPE compute_dt(FTYPE *CUf, FTYPE dtin);

static void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G, FTYPE *Tgas, FTYPE *chieffreturn, FTYPE *Gabs);
static void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu, FTYPE *Tgas, FTYPE *chieffreturn, FTYPE *Gabs);
void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM]);
static int simplefast_rad(int dir, struct of_geom *geom,struct of_state *q, FTYPE vrad2,FTYPE *vmin, FTYPE *vmax);

static void calc_kappa_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa, FTYPE *kappaes, FTYPE *Tgas);









// KORALTODO:  If solve for full RHO+MHD+RAD solution and iterate primitives instead, then can nominally better avoid out of bounds p(U) inversion.  While involves 1+4+4=9 dimensional Newton's method, use of p(U) is avoided completely so saves lots of time.  Only ever need to call U(p).  But then doesn't take advantage of accurate(and reductions) for p(U).

#define DIMTYPEFCONS 0
#define DIMTYPEFPRIM 1

// mnemonics for return modes so schemes know how failed and what to do.
// worse failure should be larger number
#define UTOPRIMGENWRAPPERRETURNNOFAIL  (UTOPRIMNOFAIL)
#define UTOPRIMGENWRAPPERRETURNFAILRAD (1)
#define UTOPRIMGENWRAPPERRETURNFAILMHD (2)

// wrapper for Utoprimgen() that returns non-zero if failed in some way so know can't continue with that method
// doradonly: ==1: Do only radiative inversion.  ==0: do full inversion.
// showmessages : 0 or 1 : whether to show messages for issues
// allowlocalfailurefixandnoreport : 0 or 1 : whether to have inversion avoid report and just use local fix
// finalstep : whether this is the final step of RK substeps
// evolvetype :
// inputtype :
// U : conserved quantity
// ptrgeom : geometry pointer
// pr : primitive (acts as guess for inversion and holds output for U->P)
// newtonstats: for inversion method report
static int Utoprimgen_failwrapper(int doradonly, int *radinvmod, int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int whichcap, int evolvetype, int inputtype,FTYPE *U,  struct of_state *qptr,  struct of_geom *ptrgeom, FTYPE dissmeasure, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int failreturn;

  // defaults
  failreturn=0;



  // KORALTODO: 
  // flag needs to be reset to preexistingfail(gas/rad) is not a failure.  Only use preexisting catches in utoprimgen.c if done with 4-force and report error in pflag and eventually go to the final advance.c inversion.
  GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
  GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
  

  PFTYPE *lpflag,*lpflagrad;
  lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
  lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);


  if(doradonly==1){
    u2p_rad(showmessages, allowlocalfailurefixandnoreport,GAMMAMAXRADIMPLICITSOLVER,whichcap,U,pr,ptrgeom,lpflag,lpflagrad);
    *radinvmod=(int)(*lpflagrad);
  }
  else{
    //calculating primitives  
    // OPTMARK: Should optimize this to  not try to get down to machine precision
    int whichmethod=MODEDEFAULT; // means don't change method from eomtype.
    MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, whichmethod, evolvetype, inputtype, U, qptr, ptrgeom, dissmeasure, pr, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);
    //    MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, eomtype, evolvetype, inputtype, U, NULL, ptrgeom, dissmeasure, pr, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);
    *radinvmod=0;  //KORALTODO: Not using method that needs this call, so for now don't pass radinvmod through to Utoprimgen().
    nstroke+=(newtonstats->nstroke);
    // this can change eomtype
  }

  // check how inversion did.  If didn't succeed, then check if soft failure and pass.  Else if hard failure have to return didn't work.
  if(IFUTOPRIMFAILSOFT(*lpflag)){
    // assume soft failure ok, but reset
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
    if(PRODUCTION==0 && showmessages && debugfail>=2) dualfprintf(fail_file,"Got soft MHD failure inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(IFUTOPRIMRADFAIL(*lpflagrad)){
    // can reduce Newton step if getting failure.
    // reset pflag for radiation to no failure, but treat here locally
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
    if(PRODUCTION==0&&showmessages && debugfail>=2) dualfprintf(fail_file,"Got some radiation inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILRAD;
  }
  else if( IFUTOPRIMFAIL(*lpflag) || IFUTOPRIMRADFAIL(*lpflagrad) ){
    // these need to get fixed-up, but can't, so return failure
    if(PRODUCTION==0&&showmessages && debugfail>=2) dualfprintf(fail_file,"Got hard failure of inversion (MHD part only considered as hard) in f_implicit(): ijk=%d %d %d : %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,*lpflag,*lpflagrad);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILMHD;
  }
  else if(PRODUCTION==0){
    // no failure
    // dualfprintf(fail_file,"No failure in Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }


  //DEBUG:
  if(PRODUCTION==0&&debugfail>=2 && showmessages){
    struct of_state q;
    MYFUN(get_stateforcheckinversion(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    int outputtype=inputtype;
    FTYPE Unew[NPR];
    MYFUN(primtoU(outputtype,pr, &q, ptrgeom, Unew, NULL),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included
    int pliter,pl;
    PLOOP(pliter,pl) dualfprintf(fail_file,"COMPARE: pl=%d pr=%g U=%g Unew=%g\n",pl,pr[pl],U[pl],Unew[pl]);
    int jj;
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d uradcon=%g uradcov=%g\n",jj,q.uradcon[jj],q.uradcov[jj]);
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d ucon=%g ucov=%g\n",jj,q.ucon[jj],q.ucov[jj]);
  }
  //DEBUG:
  if(PRODUCTION==0 && (showmessages || debugfail>=2)){
    static int maxlntries=0,maxnstroke=0;
    int diff;
    diff=0;
    // For RADSHADOW, gets up to 5
    if(newtonstats->lntries>maxlntries){ maxlntries=newtonstats->lntries; diff=1;}
    if(newtonstats->nstroke>maxnstroke){ maxnstroke=newtonstats->nstroke; diff=1;}
    // only report if grew beyond prior maximum
    if(diff) dualfprintf(fail_file,"newtonsteps: lntries=%d (max=%d) nstroke=%d (max=%d) logerror=%g\n",newtonstats->lntries,maxlntries,newtonstats->nstroke,maxnstroke,newtonstats->lerrx);
  }

  // return failure mode of inversion U->P
  return(failreturn);
}




















//////////////////////////////////////////////////////
//
// whether iterate or have as error function: MHD T^t_\nu or RAD R^t_\nu, etc.
//
//////////////////////////////////////////////////////
#define QTYUMHD 0 // iter or ferr
#define QTYUMHDENERGYONLY 1 // ferr only for now
#define QTYUMHDMOMONLY 2 // ferr only for now
#define QTYURAD 3 // iter or ferr
#define QTYURADENERGYONLY 4 // ferr only for now
#define QTYURADMOMONLY 5 // ferr only for now
#define QTYPMHD 6 // only iter
#define QTYPMHDENERGYONLY 7 // iter
#define QTYPMHDMOMONLY 8 // iter
#define QTYPRAD 9 // only iter
#define QTYPRADENERGYONLY 10 // only iter
#define QTYPRADMOMONLY 11 // only iter
#define QTYENTROPYUMHD 12 // iter or ferr
#define QTYENTROPYUMHDENERGYONLY 13 // ferr only for now
#define QTYENTROPYUMHDMOMONLY 14 // ferr only for now
#define QTYENTROPYPMHD 15 // iter (not used)
#define QTYENTROPYPMHDENERGYONLY 16 // iter (not used)
#define QTYENTROPYPMHDMOMONLY 17 // iter (not used)


// P or U types
#define IMPPTYPE(implicititer) (implicititer==QTYPMHD ||implicititer==QTYPMHDENERGYONLY ||implicititer==QTYPMHDMOMONLY || implicititer==QTYPRAD ||implicititer==QTYPRADENERGYONLY ||implicititer==QTYPRADMOMONLY   || implicititer==QTYENTROPYPMHD ||implicititer==QTYENTROPYPMHDENERGYONLY ||implicititer==QTYENTROPYPMHDMOMONLY)

#define IMPUTYPE(implicititer) (implicititer==QTYUMHD ||implicititer==QTYUMHDENERGYONLY ||implicititer==QTYUMHDMOMONLY  || implicititer==QTYURAD ||implicititer==QTYURADENERGYONLY ||implicititer==QTYURADMOMONLY || implicititer==QTYENTROPYUMHD ||implicititer==QTYENTROPYUMHDENERGYONLY ||implicititer==QTYENTROPYUMHDMOMONLY)

// MHD or RAD type
#define IMPMHDTYPE(implicititer) (implicititer==QTYPMHD ||implicititer==QTYPMHDENERGYONLY ||implicititer==QTYPMHDMOMONLY || implicititer==QTYUMHD ||implicititer==QTYUMHDENERGYONLY ||implicititer==QTYUMHDMOMONLY || implicititer==QTYENTROPYUMHD ||implicititer==QTYENTROPYUMHDENERGYONLY ||implicititer==QTYENTROPYUMHDMOMONLY)

#define IMPRADTYPE(implicititer) (implicititer==QTYPRAD ||implicititer==QTYPRADENERGYONLY ||implicititer==QTYPRADMOMONLY || implicititer==QTYURAD ||implicititer==QTYURADENERGYONLY ||implicititer==QTYURADMOMONLY)

// PMHD types
#define IMPPMHDTYPE(implicititer) (implicititer==QTYPMHD ||implicititer==QTYPMHDENERGYONLY ||implicititer==QTYPMHDMOMONLY)



// set default method
// can control method per iteration
static void define_method(int iter, int *eomtype, int itermode, FTYPE fracenergy, FTYPE dissmeasure, int *implicititer, int *implicitferr, int *BEGINMOMSTEPS, int *ENDMOMSTEPS, int *BEGINENERGYSTEPS, int *ENDENERGYSTEPS, int *BEGINFULLSTEPS, int *ENDFULLSTEPS, int *BEGINNORMALSTEPS)
{
  int eomtypelocal=*eomtype; // default, but not changing it so far.


  if(eomtypelocal==EOMDEFAULT){
    eomtypelocal=EOMTYPE; // override
  }


  if(PRODUCTION==0&&EOMDONOTHING(eomtypelocal)){
    dualfprintf(fail_file,"Can't have EOMDONOTHING in radiation code.\n");
    myexit(938463651);
  }



  *BEGINMOMSTEPS=-1;
  *ENDMOMSTEPS=-1;

  *BEGINENERGYSTEPS=-1;
  *ENDENERGYSTEPS=-1;

  *BEGINFULLSTEPS=1;
  *ENDFULLSTEPS=(IMPMAXITER*2);

  *BEGINNORMALSTEPS=1; // has to be consistent with actual first iteration, which is currently 1.

  *implicititer=(QTYURAD); // choice
  *implicitferr=(QTYURAD); // choice


}



#define JACLOOP(jj,startjj,endjj) for(jj=startjj;jj<=endjj;jj++)
#define JACLOOPALT(jj,startjj,endjj) DLOOPA(jj) //for(jj=startjj;jj<=endjj;jj++) // for those things might or might not want to do all terms
#define JACLOOPFULLERROR(itermode,jj,startjj,endjj) for(jj=(itermode==ITERMODECOLD ? startjj : 0);jj<=(itermode==ITERMODECOLD ? endjj : NDIM-1);jj++)
#define JACLOOPSUBERROR(jj,startjj,endjj) JACLOOP(jj,startjj,endjj)
#define JAC2DLOOP(ii,jj,startjj,endjj) JACLOOP(ii,startjj,endjj) JACLOOP(jj,startjj,endjj)

static void get_refUs(int *numdims, int *startjac, int *endjac, int *implicititer, int *implicitferr, int *irefU, int *iotherU, int *erefU, int *eotherU, int *signgd2, int *signgd4, int *signgd6, int *signgd7)
{
  int jj;

  // default
  *numdims=NDIM;
  *signgd7= (+1.0); // not used for PMHD
  DLOOPA(jj) irefU[jj]=UU+jj;
  DLOOPA(jj) iotherU[jj]=URAD0+jj;
  *startjac=0; *endjac=NDIM-1;


  *numdims=NDIM;
  *signgd7= (-1.0);
  DLOOPA(jj) irefU[jj]=URAD0+jj;
  DLOOPA(jj) iotherU[jj]=UU+jj;
  *startjac=0; *endjac=NDIM-1;
  
  //default
  *numdims=NDIM;
  DLOOPA(jj) erefU[jj]=UU+jj;
  DLOOPA(jj) eotherU[jj]=URAD0+jj;

  *numdims=NDIM;
  DLOOPA(jj) erefU[jj]=URAD0+jj;
  DLOOPA(jj) eotherU[jj]=UU+jj;

  // sign that goes into implicit differencer that's consistent with sign for *signgd of -1 when using the radiative uu to measure f.
  *signgd2=(+1.0);
  *signgd4=(+1.0); // for entropy alone for Gdpl in error function // Appears for QTYUMHD,QTYENTROPYUMHD this sign is the right one.  But both cases have lots of cold MHD inversions.
  *signgd6=(-1.0); // for entropy as goes into GS from dUrad or dUmhd  //  // KORALTODO SUPERGODMARK:  -- unsure about sign!

}

//        failreturnferr=f_implicit(iter, failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, radinvmod, pp, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, &goexplicit); // modifies uu and pp

static int f_implicit(int iter, int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom, struct of_state *q,  FTYPE *fpl, FTYPE *fplnorm, int *goexplicit)
//static int f_implicit_lab(int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f, FTYPE *fnorm)
{
  //  FTYPE pp[NPR];
  int pliter, pl;
  int iv;
  FTYPE f[NDIM],fnorm[NDIM];
  struct of_newtonstats newtonstats;
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.


  FTYPE pp[NPR];
  PLOOP(pliter,pl) pp[pl]=pp0[pl];

  // get primitive (don't change uu0).  This pp is used for inversion guess and to hold final inversion answer.
  //  PLOOP(pliter,pl) pp[pl] = pp0[pl];

  // get change in conserved quantity between fluid and radiation (equal and opposite 4-force)
  // required for inversion to get P(U) for MHD and RAD variables
  // this preserves machine accurate conservation instead of applying 4-force on each fluid and radiation separately that can accumulate errors
  DLOOPA(iv) uu[UU+iv] = uu0[UU+iv] - (uu[URAD0+iv]-uu0[URAD0+iv]);

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"f_implicit_lab: wc=%d i=%d j=%d pl=%d uu=%g\n",whichcall, ptrgeom->i,ptrgeom->j,pl,uu[pl]);
  
  // Get P(U)
  int failreturn;
  //  dualfprintf(fail_file,"BEFORE 0\n");
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d uu=%21.15g pp=%21.15g\n",pl,uu[pl],pp[pl]);
  failreturn=Utoprimgen_failwrapper(0,radinvmod,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype,whichcap,EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure,pp, &newtonstats);
  if(failreturn && failreturn>failreturnallowable){
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of f_implicit_lab(): %d vs. %d\n",failreturn,failreturnallowable);
    return(failreturn);
  }


  // get 4-force for all pl due to radiation
  FTYPE Gdpl[NPR];
  koral_source_rad_calc(1,0,pp, ptrgeom, Gdpl, NULL, NULL, NULL, NULL);


  // compute difference vector between original and new 4-force's effect on conserved radiative quantities
  // NR1992 Eq. 16.6.16: y_{n+1} = y_n + h f(y_{n+1}) , so error function is f = (y_{n+1} - y_n) - h f(y_{n+1})
  // i.e. f->0 as change in conserved quantity approaches the updated value of 4-force
#define SIGNGD2 (1.0) // sign that goes into implicit differencer that's consistent with sign for SIGNGD of -1 when using the radiative uu to measure f.
  DLOOPA(iv) f[iv] = (uu[URAD0+iv] - uu0[URAD0+iv]) + (SIGNGD2 * localdt * Gdpl[URAD0+iv]);

  // get error normalization that involves actual things being differenced
  DLOOPA(iv) fnorm[iv] = 0.5*(fabs(uu[URAD0+iv]) + fabs(uu0[URAD0+iv]) + fabs(SIGNGD2 * localdt * Gdpl[URAD0+iv]));

  struct of_state qalt;
  get_state(pp, ptrgeom, &qalt);
  primtoU(UNOTHING,pp,&qalt,ptrgeom, uu, NULL);


  // save better guess for later inversion (including this inversion above) from this inversion
  //  PLOOP(pliter,pl) pp0[pl]=pp[pl];

  PLOOP(pliter,pl) pp0[pl]=pp[pl];


  FTYPE sign[NPR],extrafactor[NPR];
  PLOOP(pliter,pl){
    sign[pl]=+1.0;
    extrafactor[pl]=1.0;
  }
  //  FTYPE fpl[NPR],fplnorm[NPR];
  PLOOP(pliter,pl){
    iv=pl-URAD0;
    if(iv>=0 && iv<NDIM){
      fpl[pl]=f[iv];
      fplnorm[pl]=fnorm[iv];
    }
  }
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"FUCK: pl=%d uu0=%21.15g uu=%21.15g Gdpl=%21.15g f=%21.15g fnorm=%21.15g\n",pl,uu0[pl],uu[pl],sign[pl]*localdt*Gdpl[pl],fpl[pl],fplnorm[pl]);


  //  PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d pl=%d uu=%g pp=%g\n",ptrgeom->i,pl,uu[pl],pp[pl]);


  //  dualfprintf(fail_file,"i=%d Gd=%Lg %Lg %Lg %Lg : uuG=%Lg %Lg %Lg %Lg\n",ptrgeom->i,Gdpl[URAD0+0],Gdpl[URAD0+1],Gdpl[URAD0+2],Gdpl[URAD0+3],SIGNGD2 * localdt * Gdpl[URAD0],SIGNGD2 * localdt * Gdpl[URAD1],SIGNGD2 * localdt * Gdpl[URAD2],SIGNGD2 * localdt * Gdpl[URAD3]);


  *goexplicit=0;

  return 0;
} 




// compute dt for this sub-step
static FTYPE compute_dt(FTYPE *CUf, FTYPE dtin)
{
  // what's applied to source and flux terms to get update (see global.stepch.h and step_ch.c:get_truetime_fluxdt() and step_ch.c:setup_rktimestep()) to get Uf
  // We don't use the ucum update version of dt.  As part of the RK method, the ucum update is separate from the substeps used to get information on updates(?). GODMARK.
  return(CUf[2]*dtin);
}








// compute changes to U (both T and R) using implicit method
// KORALTODO: If doing implicit, should also add geometry source term that can sometimes be stiff.  Would require inverting sparse 8x8 matrix (or maybe 6x6 since only r-\theta for SPC).  Could be important for very dynamic radiative flows.
//static int koral_source_rad_implicit(FTYPE *pin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
static int koral_source_rad_implicit(int *eomtype, FTYPE *pin, FTYPE *pf, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int i1,i2,i3,iv,ii,jj,kk,pliter,sc;
  FTYPE iJ[NDIM][NDIM];
  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uu[NPR]; 
  FTYPE uuporig[NPR],uu0orig[NPR];
  FTYPE f1[NDIM],f1norm[NDIM],f3report[NDIM];
  FTYPE x[NDIM];
  FTYPE realdt;
  FTYPE radsource[NPR], deltas[NDIM]; 
  int pl;
  FTYPE bestuu[NPR],lowestfreport[NDIM];
  int gotbest;


  // static counter for diagnosing issues
  static long long int numimplicits=0;
  static long long int numoff1iter=0,numofiter=0;
  static long long int failnum=0;
#define NUMNUMHIST (20)
  static long long int numhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
  static long long int numhistiter[IMPMAXITER+1]={0}; // histogram of error for implicit solver to be reported infrequently
  numimplicits++;
  int showmessages=0; // by default 0, don't show any messages for inversion stuff during implicit solver, unless debugging.  Assume any moment of inversion failure is corrected for now unless failure of final inversion done outside implicit solver.
  int showmessagesheavy=0;  // very detailed for common debugging
  int allowlocalfailurefixandnoreport=0; // must be 0 so implicit method knows when really failure



  
  int doingitsomecpu=0;
  int doingit=0;
#if(0)
  if(nstep==15 && steppart==0 && ptrgeom->i==6 && ptrgeom->j==8 && ptrgeom->k==0){
    doingitsomecpu=1;
    if(myid==0){ // so similar situation and grid at least
      dualfprintf(fail_file,"DOINGIT\n");
      doingit=1;

      dt=  0.0039533717184373549738;CUf[2]=                       0.5;pin[0]= 5.8671573019816740585e-10;Uiin[0]= 1.8599071962339443935e-09;Ufin[0]=                         0;dUother[0]=-4.8369293236147600752e-09;pin[1]= 0.00014782262868803244442;Uiin[1]= -0.0016713894984997005364;Ufin[1]=                         0;dUother[1]=   0.021472966823648780533;pin[2]=   -0.97252338339604188115;Uiin[2]= -0.0024519334986580419117;Ufin[2]=                         0;dUother[2]=   0.032932795156698776105;pin[3]=  -0.018379994918252814202;Uiin[3]=-0.00024591049864253674949;Ufin[3]=                         0;dUother[3]=   0.002978006528940405943;pin[4]=   -0.11277198666607900795;Uiin[4]= 0.00035374677585359898144;Ufin[4]=                         0;dUother[4]= -0.0062075118422483781515;pin[5]=                         0;Uiin[5]=                         0;Ufin[5]=                         0;dUother[5]=                         0;pin[6]=                         0;Uiin[6]=                         0;Ufin[6]=                         0;dUother[6]=                         0;pin[7]=                         0;Uiin[7]=                         0;Ufin[7]=                         0;dUother[7]=                         0;pin[8]=     8.8758605434764514036;Uiin[8]=    -8.9393883704115929432;Ufin[8]=                         0;dUother[8]=    -5.0519156634008387076;pin[9]= -0.0011231891505075573715;Uiin[9]=   -0.20132916332265980042;Ufin[9]=                         0;dUother[9]=    -1.8634463818029622854;pin[10]= 0.00031504056189326296922;Uiin[10]=    0.30458156264522790929;Ufin[10]=                         0;dUother[10]=     3.6127343678807874675;pin[11]= 0.00034068784347714866871;Uiin[11]=     1.0039678561251826477;Ufin[11]=                         0;dUother[11]=    -1.4134718906839051758;pin[12]= 0.00014782262868803244442;Uiin[12]= 1.0280003452952356973e-07;Ufin[12]=                         0;dUother[12]=  -4.04686258718244875e-07;ptrgeom->gcov[GIND(0,0)]=     0.1112097304662429864;ptrgeom->gcon[GIND(0,0)]=    -2.1112097304662429864;ptrgeom->gcov[GIND(0,1)]=     1.8347557720715785915;ptrgeom->gcon[GIND(0,1)]=    0.67299805449784306957;ptrgeom->gcov[GIND(0,2)]=                         0;ptrgeom->gcon[GIND(0,2)]=  0.0076778296331190674917;ptrgeom->gcov[GIND(0,3)]=    -4.7169723403409002565;ptrgeom->gcon[GIND(0,3)]=                         0;ptrgeom->gcov[GIND(1,0)]=     1.8347557720715785915;ptrgeom->gcon[GIND(1,0)]=    0.67299805449784306957;ptrgeom->gcov[GIND(1,1)]=      5.762700808386148114;ptrgeom->gcon[GIND(1,1)]=    0.03738451761499675204;ptrgeom->gcov[GIND(1,2)]=   -0.61634522990604687586;ptrgeom->gcon[GIND(1,2)]=  0.0004264974545557230198;ptrgeom->gcov[GIND(1,3)]=    -14.797246198131100303;ptrgeom->gcon[GIND(1,3)]=   0.030408359722302656999;ptrgeom->gcov[GIND(2,0)]=                         0;ptrgeom->gcon[GIND(2,0)]=  0.0076778296331190674917;ptrgeom->gcov[GIND(2,1)]=   -0.61634522990604687586;ptrgeom->gcon[GIND(2,1)]=  0.0004264974545557230198;ptrgeom->gcov[GIND(2,2)]=     54.025572387868411831;ptrgeom->gcon[GIND(2,2)]=   0.018514618641898817149;ptrgeom->gcov[GIND(2,3)]=                         0;ptrgeom->gcon[GIND(2,3)]=  0.0003469106690726516806;ptrgeom->gcov[GIND(3,0)]=    -4.7169723403409002565;ptrgeom->gcon[GIND(3,0)]=                         0;ptrgeom->gcov[GIND(3,1)]=    -14.797246198131100303;ptrgeom->gcon[GIND(3,1)]=   0.030408359722302656999;ptrgeom->gcov[GIND(3,2)]=                         0;ptrgeom->gcon[GIND(3,2)]=  0.0003469106690726516806;ptrgeom->gcov[GIND(3,3)]=     122.58803675566230759;ptrgeom->gcon[GIND(3,3)]=   0.011827907711600354189;ptrgeom->gcovpert[0]=     1.1112097304662429864;ptrgeom->beta[0]=                         0;ptrgeom->gcovpert[1]=      4.762700808386148114;ptrgeom->beta[1]=     0.3187736608002545943;ptrgeom->gcovpert[2]=     53.025572387868411831;ptrgeom->beta[2]=  0.0036366967820972874747;ptrgeom->gcovpert[3]=     121.58803675566230759;ptrgeom->beta[3]=                         0;ptrgeom->alphalapse=    0.68823112663654761439;ptrgeom->betasqoalphasq=     1.2347870650828603863;ptrgeom->gdet=     111.59068711266855677;ptrgeom->igdetnosing=  0.0089613212883109220732;ptrgeom->i=6;ptrgeom->j=8;ptrgeom->k=0;ptrgeom->p=0;q->ucon[0]=      3.170031244271817163;q->ucov[0]=    -2.7539241537565477709;q->ucon[1]=    -1.9830458479837551403;q->ucov[1]=    -3.9243180145269111608;q->ucon[2]=  -0.029908437343443991947;q->ucov[2]=   -0.39357959761648096219;q->ucon[3]=   -0.11277198666607900795;q->ucov[3]=    0.56617149112032208284;q->uradcon[0]=     1.4530283094901191383;q->uradcov[0]=   -0.69191220961803263223;q->uradcon[1]=   -0.46431034261307814926;q->uradcov[1]=  -0.011708018669792006875;q->uradcon[2]= -0.0049691828155257148169;q->uradcov[2]=   0.017712518956876775202;q->uradcon[3]= 0.00034068784347714866871;q->uradcov[3]=   0.058384360265514081798;q->pressure= 4.9274209562677481479e-05;q->entropy=   3.24287133495233427e-08;q->ifremoverestplus1ud0elseud0=     -1.753924153756547771;
      showmessages=showmessagesheavy=1;
    }// end on doing it core
  }

#endif
  //  showmessages=showmessagesheavy=1;





  if(showmessagesheavy){
    dualfprintf(fail_file,"DOING: nstep=%ld steppart=%d ijk=%d %d %d\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }


  //////////////
  // setup reversion to best solution for uu in case iterations lead to worse error and reach maximum iterations.
  gotbest=0;
  DLOOPA(jj) lowestfreport[jj]=BIG;
  // setup locally-used pinuse that can pass back as pin if good solution
  FTYPE pinuse[NPR];
  int gotfirstnofail=0;
  FTYPE pinusefirst[NPR];
  PLOOP(pliter,pl) pinuse[pl]=pinusefirst[pl]=pin[pl];

  ///////////////////
  // setup implicit iteration procedure and loops
  realdt = compute_dt(CUf,dt);
  FTYPE DAMPFACTOR=1.0; // factor by which step Newton's method.
  FTYPE fracdtuu0=1.0,fracdtG=1.0,fracuup=1.0; // initially try full realstep step
  FTYPE fracdtuu0p=fracdtuu0;
  FTYPE fracdtuu0pp=fracdtuu0p;

  FTYPE fracdtGp=fracdtG;
  FTYPE fracdtGpp=fracdtGp;

  // default is to allow no failure
  int failreturnallowable=UTOPRIMNOFAIL;
  int failreturnallowableuse=failreturnallowable;
  int failreturnallowablefirst=-1;

  // see if uu0->p possible and what type of failure one gets as reference for what failure to allow.
  // KORALTODO: Need to check if UFSET with no dUother fails.  How it fails, must allow since can do nothing better.  This avoids excessive attempts to get good solution without that failure!  Should speed-up things.  But what about error recovery?  If goes from CASE to no case!
  if(USEDUINRADUPDATE){
    // uu0 will hold original vector of conserved
    // here original means U[before fluxes, geometry, etc.] + dU[due to fluxes, geometry, etc. already applied and included in dUother]
    // This is required for stiff source term so immediately have balance between fluxes+geometry+radiation.  Otherwise, radiation diffuses.
    // I'm guessing that even though one uses RK2 or RK3, the first step generates large radiative velocities without any balanced source term since U isn't updated yet.  One would hope RK2 would recover on the final substep, but it doesn't!  In RK2, upon the final substep, that velocity is present for the radiation source term.  But it's also present for the fluxes!  That is, if there were no flux update on the final substep, then the source would have balanced the previous flux, but yet another flux is done, so there can be no balance.  This leads to a run-away velocity that would be similar to the \tau\sim 1 case.
    // NOTE: If this gives radiation or mhd failure, then less likely that will be actual solution since not even original uu0 has inversion.
    PLOOP(pliter,pl) uu[pl]=uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
    // Note that "q" isn't used in this function or used in function call, so don't have to update it here.

    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)
    struct of_newtonstats newtonstats;
    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    int failreturn;
    int finalstep = 1;
    FTYPE prtest[NPR];
    PLOOP(pliter,pl) prtest[pl]=pinuse[pl]; // initial guess
    //    dualfprintf(fail_file,"BEFORE 1\n");
    //    failreturnallowable=Utoprimgen_failwrapper(showmessages,allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, Uiin, ptrgeom, prtest, &newtonstats);
    int doradonly=0;
    int radinvmod=0;
    //    eomtypelocal=*eomtype; // stick with default choice so far
    int whichcap=CAPTYPEFIX1;
    //    FTYPE dissmeasure=-1.0;
    failreturnallowable=Utoprimgen_failwrapper(doradonly,&radinvmod,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu0, q, ptrgeom, dissmeasure, prtest, &newtonstats);
    if(failreturnallowable!=UTOPRIMGENWRAPPERRETURNNOFAIL){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is already a problem with %d\n",failreturnallowable);
    }
    else{
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is NOT a problem with %d\n",failreturnallowable);
    }

    // set first pinuse that is p(uu0), assuming that can do inversion
    PLOOP(pliter,pl) pinusefirst[pl]=prtest[pl];
    if(failreturnallowable==UTOPRIMGENWRAPPERRETURNNOFAIL) gotfirstnofail=1;
    else gotfirstnofail=0;

    failreturnallowablefirst=failreturnallowable;

  }
  else{
    // then (not recommended) just using Uiin as uu0
    PLOOP(pliter,iv) uu[iv] = uu0[iv] = Uiin[iv];
  }  

  

  


  ////////////////////////////////
  // START IMPLICIT ITERATIONS
  ////////////////////////////////
  int iter=0;
  int failreturn;
  int f1iter;
  int checkconv,changeotherdt;
  FTYPE impepsjac=IMPEPS;
  FTYPE errorabs=BIG;
  FTYPE errorabsp=errorabs;
  int convreturn;
  FTYPE f3[NDIM],f3norm[NDIM];

  // initialize previous 'good inversion' based uu's
  PLOOP(pliter,pl)  uupp[pl]=uuporig[pl]=uup[pl]=uu0orig[pl]=uu[pl];
  
  do{
    iter++;


    
    if(iter>10){ // KORALTODO: improve upon this later
      // assume trying hard and failing to work, then allow CASE radiation errors
      failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
      failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
    }

    // cautious first step to get reasonable error measurement. Allows often only 1 iteration to get sufficiently small error, while DAMPFACTOR=1 would already go beyond point where error actually increases.
    if(iter==1){
      DAMPFACTOR=0.37;
    }


#if(0)
    // KORALTODO: improve on this later.
    // Problem is can start jumping too far in steps for near tough spots.
    // While in normal inversion routine damping is to avoid breaching into unphysical solution space with nan/inf, below is to avoid oscillating around solution

    //    if(iter>10){
    //      impepsjac=IMPEPS*1E-10;
    //    }
    if(iter>10){
      DAMPFACTOR=0.37;
    }
    if(iter>30){
      DAMPFACTOR=0.23;
    }
    if(iter>40){
      DAMPFACTOR=0.17;
    }
    if(iter>50){
      DAMPFACTOR=0.07;
    }
#endif


    //vector of conserved at the previous two iterations
    PLOOP(pliter,pl)  uupp[pl]=uup[pl]; // uupp will have solution for inversion: P(uupp)
    PLOOP(pliter,pl)  uup[pl]=uu[pl]; // uup will not necessarily have P(uup) because uu used Newton step.
    PLOOP(pliter,pl)  uuporig[pl]=uu[pl];

    fracdtuu0pp=fracdtuu0p; // fracdtuu0 used when computing previous f1 and f2's
    fracdtuu0p=fracdtuu0; // fracdtuu0 used when computing previous f1 and f2's

    fracdtGpp=fracdtGp; // fracdtG used when computing previous f1 and f2's
    fracdtGp=fracdtG; // fracdtG used when computing previous f1 and f2's

    errorabsp=errorabs;



    // KORALTODO: Once use certain uu0 value and succeed in getting f1, not enough.  But if take a step and *then* good f1, then know that original U was good choice and can restart from that point instead of backing up uu0 again.
    // KORALTODO: Look at whether bounding error by sign as in bisection (only true in 1D!!), and if approaching 0 slowly enough and still hit CASE, then must be real CASE not a stepping issue.
    // KORALTODO: Getting f1 is just about f1(U) as far as radiation is concerned.  So as far as CASE issues, getting f1(U) means we are good with that U and we can certainly stick with the used uu0.


    /////////////////
    //
    // get error function (f1) and inversion (uu->pinuse) using uu
    //
    /////////////////
    for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
      int whichcall=1;
      //      PLOOP(pliter,pl) dualfprintf(fail_file,"BEFORE f1 iter=%d f1iter=%d pl=%d uu=%21.15g pinuse=%21.15g pin=%21.15g\n",iter,f1iter,pl,uu[pl],pinuse[pl],pin[pl]);
      int radinvmod=0;
      FTYPE fpl1[NPR],fpl1norm[NPR];
      //static int f_implicit(int iter, int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom, struct of_state *q,  FTYPE *fpl, FTYPE *fplnorm, int *goexplicit)
      int goexplicit=0;
      failreturn=f_implicit(iter,failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, eomtype, CAPTYPEFIX1, ITERMODENORMAL, 1.0,dissmeasure, &radinvmod, pinuse, uu0, uu, fracdtG*realdt, ptrgeom, q, fpl1, fpl1norm,&goexplicit); // modifies uu and pinuse
      DLOOPA(ii){
        f1[ii] = fpl1[URAD0+ii];
        f1norm[ii] = fpl1norm[URAD0+ii];
      }

      if(failreturn){

#define BACKUPRELEASEFAIL (1E-5)
        // if backing up alot, then allow same failure as original Uiin in hopes that can recover that way (otherwise would have hoped would recover via flux update to Uiin)
        if(fracdtuu0<BACKUPRELEASEFAIL){
          failreturnallowableuse=failreturnallowable;
        }


        if(iter==1){
          // if initial uu failed, then should take smaller jump from Uiin->uu until settled between fluid and radiation.
          // If here, know original Uiin is good, so take baby steps in case uu needs heavy raditive changes.
          // if f1 fails, try going back to Uiin a bit
          fracdtuu0*=RADDAMPDELTA; // DAMP Uiin->uu0 step that may be too large and generated too large G

          // modifies uup and uu as if starting over, and then another call to f_implicit_lab(f1) will change uu by generally smaller amount.
          PLOOP(pliter,pl) uup[pl]=uu[pl]=uu0[pl];
        }
        else{
          // if here, then assume prior uup was good in sense that no failures like P(uup) is good inversion.  And uu=uup before f_implicit_lab(f1) is called.
          // No, not necessarily true, because uu updated with some-sized Newton step without checking if inversion is good for that uu.
          // So need to damp between original uup and uupp .  This essentially damps Newton step now that have knowledge the Newton step was too large as based upon P(U) failure.

          // Avoid G-damping because not needed so far.  If added, competes in non-trivial way with fracuup damping that's required separately for large forces in some problems beyond iter=1 (e.g. NTUBE=31).
          //          fracdtG*=RADDAMPDELTA; // DAMP give only fraction of 4-force to let uu to catch-up
          //           fracdtG=0.5; // DAMP give only fraction of 4-force to let uu to catch-up

          fracuup*=RADDAMPDELTA; // DAMP in case Newton step is too large after iter>1 and stuck with certain uu from Newton step that gets stored in uup above.

          PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uuporig[pl];
          //          PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uup[pl];
          PLOOP(pliter,pl) uup[pl]=uu[pl]; // store new version of prior Newton step
          
          // get interpolated fracdtuu0 so using fracdtuu0 that was used with the corresponding uu
          fracdtuu0=(1.0-fracuup)*fracdtuu0pp + fracuup*fracdtuu0p;
          // same for fracdtG
          fracdtG=(1.0-fracuup)*fracdtGpp + fracuup*fracdtGp;

        }


        // get uu0 (which may be changing)
        PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);



        // keep below so can count inversion failures against retry successes in the failure file.
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f1 failed: iter=%d  Backing-up both uu0 and G.: f1iter=%d fracdtuu0=%g fracdtG=%g fracuup=%g\n",iter,f1iter,fracdtuu0,fracdtG,fracuup);
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
      }// end if failed inversion in f_implicit_lab()
      else{
        // then success, so was able to do inversion P(U) with change in radiation on fluid: P(uu[fluid] = uu0[fluid] - (uu[rad]-uu0[rad]))
        // This doesn't necessarily mean could do P(uu0) except for iter=1.
        // Corresponds to success for a certain P(uu0,uu) pair.
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"SUCCESS: pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g: pinuse(pnew)=%g : fracdtuu0=%g fracdtG=%g fracuup=%g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl],pinuse[pl],fracdtuu0,fracdtG,fracuup);
        break;
      }
    }// end loop over f1iter
    if(f1iter==MAXF1TRIES){
      if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES: nstep=%ld steppart=%d ijk=%d %d %d : iter=%d\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter);
      failnum++;
      //      mathematica_report_check(0, failnum, gotfirstnofail, realdt, ptrgeom, pinuse,pin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
      if(doingit==1) myexit(10000000); // DEBUG
      // Note that if inversion reduces to entropy or cold, don't fail, so passes until reached this point.  But convergence can be hard if flipping around which EOMs for the inversion are actually used.
      return(1);
    }
    else{
      // restore fracuup back to 1 since this is only meant to adjust how much go back to previous uu to be able to get f1 computed.
      // fracuup doesn't stay <1.0 because each attempt to get f1 is independent.
      // KORALNOTE: Perhaps reasonable and safer to keep fracuup as not reverted back to 1, since apparently unable to take full steps.  This effectively damps stepping.
      fracuup=1.0;
    }

    // diagnose
    numoff1iter += f1iter;
    if(showmessagesheavy) dualfprintf(fail_file,"i=%d f1: %g %g %g %g\n",ptrgeom->i,f1[0],f1[1],f1[2],f1[3]);



    // see if pre-convergence (happens if force is small or no force at all.  Can't necessarily continue since Jacobian can require arbitrarily large dU on fluid and fail to invert even if no fluid-radiation interaction!
    //test pre-convergence using initial |dU/U|
    // KORALTODO: This isn't a completely general error check since force might be large for fluid that needs itself to have more accuracy, but if using ~NUMEPSILON, won't resolve 4-force of radiation on fluid to better than that.
    FTYPE LOCALPREIMPCONV=(10.0*NUMEPSILON); // more strict than later tolerance
    if(f_error_check_wrap(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,realdt,f1,f1norm,f3report,Uiin, uu0,uu,ptrgeom)){
      dualfprintf(fail_file,"Initial error near machine precision.\n");
      break;
    }


    int notfinite= !isfinite(uu[0])|| !isfinite(uu[1])|| !isfinite(uu[2])|| !isfinite(uu[3]) || !isfinite(uup[0])|| !isfinite(uup[1])|| !isfinite(uup[2])|| !isfinite(uup[3]);
    if(!notfinite){
    
      /////////
      //
      // get Jacobian and inversion Jacobian 
      //
      /////////
      //      int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, impepsjac, uu, uup, uu0, pinuse, fracdtG, realdt, ptrgeom, f1, f1norm, iJ);
      FTYPE fpl1[NPR],fpl1norm[NPR],iJpl[NPR][NPR];
      DLOOPA(ii){
        fpl1[URAD0+ii]=f1[ii];
        fpl1norm[URAD0+ii]=f1norm[ii];
      }
      int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, eomtype, CAPTYPEFIX1, ITERMODENORMAL, 1.0, -1.0, IMPEPS, iter, 1E-6, NULL, uu, uup, uu0, pinuse, pinuse, fracdtG, realdt, ptrgeom, q, fpl1, fpl1norm, iJpl);
      DLOOP(ii,jj){
        iJ[ii][jj]=iJpl[URAD0+ii][URAD0+jj];
      }

      if(failreturniJ!=0) return(failreturniJ);

      if(showmessagesheavy){
        int iii,jjj;
        DLOOP(iii,jjj) dualfprintf(fail_file,"iJ[%d][%d]=%g\n",iii,jjj,iJ[iii][jjj]);
      }


      /////////
      //
      //updating x, start with previous uu = uup
      //
      /////////
      DLOOPA(ii) x[ii]=uup[ii+URAD0];

      /////////
      //
      // step forward uu=x
      // DAMPFACTOR unused so far because don't know a priori whether to damp.  fracuup does post-inversion effective damping of this Newton step.
      // Newton step: U = U0 - (df/dU)^{-1}|_{U=U0} f(U0)
      //
      /////////
      DLOOPA(ii){
        DLOOPA(jj){
          x[ii] -= DAMPFACTOR*iJ[ii][jj]*f1[jj];
        }
      }

      /////////
      //
      // assign new uu
      //
      /////////
      PLOOP(pliter,pl) uu[pl]=uu0[pl];
      DLOOPA(ii) uu[ii+URAD0]=x[ii];
    
    
      if(showmessagesheavy){
        dualfprintf(fail_file,"POSTDX: uu: %g %g %g %g : uup=%g %g %g %g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
      }

      /////////
      //
      // try to increase DAMPs if was damping or check convergence if no damping active.
      //
      /////////
      checkconv=1;
      changeotherdt=1;
      //    if(fracuup!=1.0){
      if(fabs(fracuup-1.0)>10.0*NUMEPSILON){
        // try increasing amount of uu used
        fracuup*=RADDAMPUNDELTA;
        checkconv=0;
        changeotherdt=0; // ensure fracuup back to 1.0 first before reverting others.
      }
      //    if(fracdtuu0!=1.0){
      if(fabs(fracdtuu0-1.0)>10.0*NUMEPSILON && changeotherdt){
        // try increasing uu0 away from Uiin to account for full dUother
        fracdtuu0*=RADDAMPUNDELTA;
        PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0
        checkconv=0;
      }
      //    if(fracdtG!=1.0){
      if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
        // try increasing amount of G applied
        fracdtG*=RADDAMPUNDELTA;
        checkconv=0;
      }      
    }// end if finite



    /////////
    //
    // test convergence after Newton step
    // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
    // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
    /////////
    //test convergence using |dU/U|
    DLOOPA(ii){
      f3[ii]=(uu[ii+URAD0]-uup[ii+URAD0]);
      f3norm[ii]=fabs(uu[ii+URAD0])+fabs(uup[ii+URAD0]);
    }
  
    convreturn=f_error_check_wrap(showmessages, showmessagesheavy, iter, IMPTRYCONV,realdt,f3,f3norm,f3report,Uiin,uup,uu,ptrgeom);
    // store error and solution in case eventually lead to max iterations and actually get worse error
    errorabs=0.0;     DLOOPA(jj) errorabs     += fabs(f3report[jj]);


    // DAMP CONTROL
    if(iter==1){
      // un-damp first step, but still don't go back to 1.0.
      DAMPFACTOR=0.7;
    }
#define LOWESTDAMP (0.05)
    // see if need to damp, but don't damp below some point.
    if(errorabs>errorabsp && DAMPFACTOR>LOWESTDAMP) DAMPFACTOR*=0.5;
  


    // check convergence
    if(checkconv){
      if(convreturn) break;
      else{
        // store error and solution in case eventually lead to max iterations and actually get worse error
        FTYPE errorabsbest=0.0; DLOOPA(jj) errorabsbest += fabs(lowestfreport[jj]);
        if(errorabsbest>errorabs && isfinite(errorabs)){
          PLOOP(pliter,pl) bestuu[pl]=uu[pl];
          DLOOPA(jj) lowestfreport[jj]=f3report[jj];
          gotbest=1;
        }

      }
    } // end if checking convergence


    /////////
    // see if took too many Newton steps or not finite results
    /////////
    int itermaxed=iter>IMPMAXITER;
    if(itermaxed || notfinite ){
      convreturn=f_error_check_wrap(showmessages, showmessagesheavy, iter, IMPALLOWCONV,realdt,f3,f3norm,f3report,Uiin,uup,uu,ptrgeom);

      FTYPE errorabsbest=0.0;
      if(gotbest){
        // see if should revert to prior best
        DLOOPA(jj) errorabsbest += fabs(lowestfreport[jj]);
        if(errorabsbest<errorabs || !isfinite(errorabs) ){
          PLOOP(pliter,pl) uu[pl]=bestuu[pl];
          errorabs=errorabsbest;
        }
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"Using best: %g %g\n",errorabs,errorabsbest);
        // get new convreturn
        convreturn=(lowestfreport[0]<IMPALLOWCONV && lowestfreport[1]<IMPALLOWCONV && lowestfreport[2]<IMPALLOWCONV && lowestfreport[3]<IMPALLOWCONV);
      }

      // KORALTODO: If convreturn doesn't work, but still (say) 10% error, might want to hold onto result in case explicit backup fails as well (which is likely), in which case *much* better to use 10% error because otherwise 4-force not accounted for, which can lead to very big changes in fluid behavior due to large flux from previous step.
      // KORALTODO: Or, perhaps should really take it as a failure and use fixups.  Probably should allow for result to be written if error<10%, but only use as super-backup in fixups.  So should set pflag still.

      if(convreturn){
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab().  But f3 was allowed error. checkconv=%d (if checkconv=0, could be issue!) : %g %g %g %g : %g %g %g %g : errorabs=%g : %g %g %g\n",IMPMAXITER,checkconv,f3report[0],f3report[1],f3report[2],f3report[3],lowestfreport[0],lowestfreport[1],lowestfreport[2],lowestfreport[3],errorabs,fracdtuu0,fracuup,fracdtG);
        // NOTE: If checkconv=0, then wasn't ready to check convergence and smallness of f3 might only mean smallness of fracuup.  So look for "checkconv=0" cases in fail output.
        break;
      }
      else{
        // KORALTODO: Need backup that won't fail.
        if(debugfail>=2){
          if(itermaxed) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab(). nstep=%ld steppart=%d ijk=%d %d %d :  Bad error.\n",IMPMAXITER,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          if(notfinite) dualfprintf(fail_file,"IMPGOTNAN at iter=%d : in solve_implicit_lab(). ijk=%d %d %d :  Bad error.\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          dualfprintf(fail_file,"checkconv=%d failreturnallowable=%d: %g %g %g %g : %g %g %g %g\n",checkconv,failreturnallowable,f3report[0],f3report[1],f3report[2],f3report[3],lowestfreport[0],lowestfreport[1],lowestfreport[2],lowestfreport[3]);
          if(1||showmessages){
            failnum++;
            //            mathematica_report_check(1, failnum, gotfirstnofail, realdt, ptrgeom, pinuse,pin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
            if(doingit==1) myexit(100000000); // DEBUG
          }
        }        
        return(1);
      }
    }//end if maximum iterations

  }// end do
  while(1);


  

  // diagnose
  numofiter+=iter;
  // i=j=k=0 just to show infrequently
  if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0) || showmessagesheavy) dualfprintf(fail_file,"numimplicits=%lld averagef1iter=%g averageiter=%g\n",numimplicits,(FTYPE)numoff1iter/(FTYPE)numimplicits,(FTYPE)numofiter/(FTYPE)numimplicits);




  ///////////////////
  //
  // get source update as "dU" = dU/dt using real dt that used during implicit iterations, and will eventually use to update U in advance.c.
  DLOOPA(jj) deltas[jj]=(uu[URAD0+jj]-uu0[URAD0+jj])/realdt;

  // apply source update as force
  PLOOP(pliter,pl) radsource[pl] = 0;
#define SIGNGD3 (1.0) // sign that goes into implicit solver (just that U += U_{new} - U_{old} so that U = U_{new})
  DLOOPA(jj) radsource[UU+jj]    = -SIGNGD3*deltas[jj];
  DLOOPA(jj) radsource[URAD0+jj] = +SIGNGD3*deltas[jj];


  // DEBUG:
  //  DLOOPA(jj) dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d implicitGd[%d]=%g %g\n",nstep,steppart,ptrgeom->i,jj,radsource[UU+jj],radsource[URAD0+jj]);



  // store source update in dUcomp for return.
  sc = RADSOURCE;
  PLOOP(pliter,pl) dUcomp[sc][pl] += radsource[pl];

  // save better guess for later inversion from this inversion
  // pinuse was modified by f_implicit_lab(f1,pinuse) with output from inversion returned through pp0
  // only use pinuse if successful with implicit method, since if not successful can be various bad reasons with no good pin
  PLOOP(pliter,pl) pin[pl]=pinuse[pl];



  if(debugfail>=2){
    // then do some diagnostics and reporting
    errorabs=0.0;
    if(gotbest) DLOOPA(jj) errorabs += fabs(lowestfreport[jj]);
    else DLOOPA(jj) errorabs += fabs(f3report[jj]);
    //    dualfprintf(fail_file,"errorabs=%g\n",errorabs);
    numhisterr[MAX(MIN((int)(-log10(SMALL+errorabs)),NUMNUMHIST-1),0)]++;
    numhistiter[MAX(MIN(iter,IMPMAXITER),0)]++;
#define HISTREPORTSTEP (20)
    if(nstep%HISTREPORTSTEP==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0){
      int histi;
      for(histi=0;histi<NUMNUMHIST;histi++){
        dualfprintf(fail_file,"numhisterr%d=%lld\n",histi,numhisterr[histi]);
      }
      for(histi=0;histi<=IMPMAXITER;histi++){
        dualfprintf(fail_file,"numhistiter%d=%lld\n",histi,numhistiter[histi]);
      }
    }
  }


  return(0);
  
}





#define FAILRETURNGOEXPLICIT -1
#define FAILRETURNNOFAIL 0
#define FAILRETURNGENERAL 1
#define FAILRETURNJACISSUE 2
#define FAILRETURNMODESWITCH 3
#define FAILRETURNNOTTOLERROR 4

#define ACCEPTASNOFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNNOTTOLERROR || failreturn==FAILRETURNGOEXPLICIT)
//#define GOODNOFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL)
#define NOTACTUALFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNMODESWITCH)
#define NOTBADFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNMODESWITCH  || failreturn==FAILRETURNNOTTOLERROR)

#define ACTUALHARDFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE)
#define ACTUALHARDORSOFTFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE || failreturn==FAILRETURNNOTTOLERROR)
#define SWITCHGOODIDEAFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE || failreturn==FAILRETURNNOTTOLERROR || failreturn==FAILRETURNMODESWITCH)



#define MODEMETHOD MODEENERGY


// wrapper for mode method
static int koral_source_rad_implicit_new(int *eomtype, FTYPE *pb, FTYPE *pf, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int pliter,pl;
  int sc;

  int failreturn,noprims;
  int havebackup,didentropyalready;
  int usedenergy=0,usedentropy=0,usedboth=0,usedcold=0,usedimplicit=0,usedexplicitgood=0,usedexplicitkindabad=0,usedexplicitbad=0;

  // set backups that might change and contaminate a fresh start
  // piin, Uiin, Ufin, CUf, ptrgeom, dUother don't change, rest can.
  int failfinalreturn;
  int eomtypelocal;

  // for implicit scheme, set pf->pb in case pf contains final time information like the magnetic field
  PLOOP(pliter,pl){
    pb[pl] = pf[pl];
  }
  // but we end up modifying pb, not pf.

  FTYPE uub[NPR]; // holds returned uu from implicit solver
  FTYPE pbbackup[NPR];
  FTYPE dUcompbackup[NUMSOURCES][NPR];
  struct of_state qbackup;
  PLOOP(pliter,pl){
    pbbackup[pl]=pb[pl];
    SCLOOP(sc) dUcompbackup[sc][pl]=dUcomp[sc][pl];
    qbackup=*q;
  }

  // It's up to inversion method to set failure flags, not utoprimgen() that just checks them mostly (it might modify them based upon doing reductions).
  // setup pflags
  PFTYPE *lpflag,*lpflagrad;
  lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
  lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);
  // set default
  *lpflag=UTOPRIMNOFAIL;
  *lpflagrad=UTOPRIMRADNOFAIL;

  // default is didn't get good primitives.  Similar, but slightly different from, failreturn
  failreturn=1;
  noprims=1;
  usedenergy=0;
  usedentropy=0;
  usedboth=0;
  usedcold=0;
  usedimplicit=0;

  // whether doing energy at all
  int eomtypecond=(*eomtype==EOMGRMHD || *eomtype==EOMDEFAULT && EOMTYPE==EOMGRMHD);


  // diags
  FTYPE errorabs;
  int iters=0;
  int f1iters=0;
  FTYPE fracenergy;
  int itermode;
  int radinvmod=0;
  int whichcap;
  int trueimpmaxiter=IMPMAXITER;
  int truenumdampattempts=NUMDAMPATTEMPTS;
  int goexplicit;

  //////////////////////////////
  //
  // MODEENERGY
  //
  //////////////////////////////

  if(MODEMETHOD==MODEENERGY && eomtypecond || MODEMETHOD==MODEDEFAULT && *eomtype==EOMGRMHD){
    havebackup=0;
    didentropyalready=0;
    eomtypelocal=*eomtype;
    errorabs=1.0;
    fracenergy=1.0;
    itermode=ITERMODENORMAL;
    whichcap=CAPTYPEFIX1;
    failreturn=koral_source_rad_implicit_mode(0,havebackup, didentropyalready, &eomtypelocal, whichcap, itermode, trueimpmaxiter,  truenumdampattempts, fracenergy, dissmeasure, &radinvmod, pb, uub, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother ,dUcomp, &errorabs, errorabs, &iters, &f1iters);
    if(ACTUALHARDFAILURE(failreturn)){
      failfinalreturn=1;
      *lpflag=UTOPRIMFAILCONV;
      *lpflagrad=UTOPRIMRADFAILCASE1A;
      // restore backups in case got contaminated
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        SCLOOP(sc) dUcomp[sc][pl]=dUcompbackup[sc][pl];
      }
      *q=qbackup;
      goexplicit=0;
      usedimplicit=1;
    }
    else if(failreturn>=0){
      failfinalreturn=0;
      noprims=0;
      *eomtype=eomtypelocal; // can be EOMDONOTHING if successful and small enough error
      usedenergy=1;
      goexplicit=0;
      usedimplicit=1;
    }
    else{
      failfinalreturn=-1; // indicates to koral_source_rad() that implicit says just do trivial explicit
      noprims=1;
      *eomtype=EOMGRMHD;
      goexplicit=1;
      usedexplicitgood=1;
      dualfprintf(fail_file,"WTF=%d\n",failreturn);
    }
  }





  // whether set some primitives (implies also failfinalreturn=0)
  if(noprims==0){
    if((pb[RHO]<=0.)&&(pb[UU]>=0.)) *lpflag= UTOPRIMFAILRHONEG;
    if((pb[RHO]>0.)&&(pb[UU]<0.))   *lpflag= UTOPRIMFAILUNEG;
    if((pb[RHO]<=0.)&&(pb[UU]<0.))  *lpflag= UTOPRIMFAILRHOUNEG;
    if(pb[PRAD0]<=0.) *lpflagrad = UTOPRIMRADFAILERFNEG;
  }
  


  // DEBUG:
  //  failfinalreturn=1;
  //  *eomtype=EOMGRMHD;



  





  if(PRODUCTION==0 && debugfail>=2){


    ////////////////
    //
    // static counter for diagnosing issues
    //
    ////////////////
#define NUMNUMHIST (20)

    static long long int numenergy=0;
    static long long int numentropy=0;
    static long long int numboth=0;
    static long long int numcold=0;
    static long long int numbad=0;
    static long long int numramesh=0;
    static long long int numrameshenergy=0;
    static long long int numrameshentropy=0;

    static long long int numimplicits=0;
    static long long int numexplicitsgood=0;
    static long long int numexplicitskindabad=0;
    static long long int numexplicitsbad=0;
    static long long int numoff1iter=0,numofiter=0;
    static long long int numhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
    static long long int numhistiter[IMPMAXITER+1]={0}; // histogram of error for implicit solver to be reported infrequently

    // static counter for diagnosing issues
    static long long int totalnumenergy=0;
    static long long int totalnumentropy=0;
    static long long int totalnumcold=0;
    static long long int totalnumboth=0;
    static long long int totalnumbad=0;
    static long long int totalnumramesh=0;
    static long long int totalnumrameshenergy=0;
    static long long int totalnumrameshentropy=0;

    static long long int totalnumimplicits=0;
    static long long int totalnumexplicitsgood=0;
    static long long int totalnumexplicitskindabad=0;
    static long long int totalnumexplicitsbad=0;
    static long long int totalnumoff1iter=0,totalnumofiter=0;
    static long long int totalnumhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
    static long long int totalnumhistiter[IMPMAXITER+1]={0}; // histogram of error for implicit solver to be reported infrequently


    ////////////////////
    //
    // Do some diagnostics and reporting.  Done if ACCEPTASNOFAILURE(failreturn)==0 or not.
    //
    // KORALNOTE: If set IMPALLOWCONV to be smaller, energy solution can be then invalidated (while previously would have been used because it has no u_g issue) and go to entropy can succeed, and then appears that fewer low-energy events.
    // But then no FAILINFO will be reported to check why energy got high error!
    //
    //////////////////////

    // simple counters
    numimplicits+=usedimplicit;
    numexplicitsgood+=usedexplicitgood;
    numexplicitskindabad+=usedexplicitkindabad;
    numexplicitsbad+=usedexplicitbad;
    numofiter+=iters;
    numoff1iter+=f1iters;
    numenergy+=usedenergy;
    numentropy+=usedentropy;
    numboth+=usedboth;
    numcold+=usedcold;
    numbad+=(usedenergy==0 && usedentropy==0 && usedboth==0 && usedcold==0 && usedimplicit==1);
    numramesh+=0;
    numrameshenergy+=0;
    numrameshentropy+=0;


    // i=j=k=0 just to show infrequently
    if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0)) dualfprintf(fail_file,"numimplicits=%lld numexplicitsgood=%lld numexplicitskindabad=%lld numexplicitsbad=%lld numenergy=%lld numentropy=%lld numboth=%lld numcold=%lld numbad=%lld numramesh=%lld numrameshenergy=%lld numrameshentropy=%lld averagef1iter=%g averageiter=%g\n",numimplicits,numexplicitsgood,numexplicitskindabad,numexplicitsbad,numenergy,numentropy,numboth,numcold,numbad,numramesh,numrameshenergy,numrameshentropy,(FTYPE)numoff1iter/(SMALL+(FTYPE)numimplicits),(FTYPE)numofiter/(SMALL+(FTYPE)numimplicits));
    
    numhisterr[MAX(MIN((int)(-log10l(SMALL+errorabs)),NUMNUMHIST-1),0)]++;
    numhistiter[MAX(MIN(iters,IMPMAXITER),0)]++;
#define HISTREPORTSTEP (20)
    if(nstep%HISTREPORTSTEP==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0){
      int histi;
      for(histi=0;histi<NUMNUMHIST;histi++){
        dualfprintf(fail_file,"numhisterr%d=%lld\n",histi,numhisterr[histi]);
      }
      for(histi=0;histi<=IMPMAXITER;histi++){
        dualfprintf(fail_file,"numhistiter%d=%lld\n",histi,numhistiter[histi]);
      }
    }

#define SHOWCOLLECTIVETOTALS 1

    if(SHOWCOLLECTIVETOTALS){
      // over all cores
      if(USEMPI){
        MPI_Reduce(&numimplicits, &totalnumimplicits, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numexplicitsgood, &totalnumexplicitsgood, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numexplicitskindabad, &totalnumexplicitskindabad, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numexplicitsbad, &totalnumexplicitsbad, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numoff1iter, &totalnumoff1iter, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numofiter, &totalnumofiter, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numenergy, &totalnumenergy, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numentropy, &totalnumentropy, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numboth, &totalnumboth, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numcold, &totalnumcold, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numbad, &totalnumbad, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numramesh, &totalnumramesh, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numrameshenergy, &totalnumrameshenergy, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numrameshentropy, &totalnumrameshentropy, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);

        MPI_Reduce(numhisterr, totalnumhisterr, NUMNUMHIST, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(numhistiter, totalnumhistiter, IMPMAXITER+1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
      }
      if(myid==MPIid[0]){
        // i=j=k=0 just to show infrequently
        if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0)&&steppart==0) trifprintf("totalnumimplicits=%lld totalnumexplicitsgood=%lld totalnumexplicitskindabad=%lld totalnumexplicitsbad=%lld totalnumenergy=%lld totalnumentropy=%lld totalnumboth=%lld totalnumcold=%lld totalnumbad=%lld totalnumramesh=%lld totalnumrameshenergy=%lld totalnumrameshentropy=%lld totalaveragef1iter=%g totalaverageiter=%g\n",totalnumimplicits,totalnumexplicitsgood,totalnumexplicitskindabad,totalnumexplicitsbad,totalnumenergy,totalnumentropy,totalnumboth,totalnumcold,totalnumbad,totalnumramesh,totalnumrameshenergy,totalnumrameshentropy,(FTYPE)totalnumoff1iter/(SMALL+(FTYPE)totalnumimplicits),(FTYPE)totalnumofiter/(SMALL+(FTYPE)totalnumimplicits));

        if(nstep%HISTREPORTSTEP==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0&&steppart==0){
          int histi;
          for(histi=0;histi<NUMNUMHIST;histi++){
            trifprintf("totalnumhisterr%d=%lld\n",histi,totalnumhisterr[histi]);
          }
          for(histi=0;histi<=IMPMAXITER;histi++){
            trifprintf("totalnumhistiter%d=%lld\n",histi,totalnumhistiter[histi]);
          }
        }
      }// end if myid==0
    }
  }


  if(failfinalreturn!=0) dualfprintf(fail_file,"BOG\n");

  return(failfinalreturn);
}







// compute changes to U (both T and R) using implicit method
// KORALTODO: If doing implicit, should also add geometry source term that can sometimes be stiff.  Would require inverting sparse 8x8 matrix (or maybe 6x6 since only r-\theta for SPC).  Could be important for very dynamic radiative flows.
static int koral_source_rad_implicit_mode(int modprim, int havebackup, int didentropyalready, int *eomtype, int whichcap, int itermode, int trueimpmaxiter, int truenumdampattempts, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pb, FTYPE *uub, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR], FTYPE *errorabsreturn, FTYPE errorabsbestexternal, int *itersreturn, int *f1itersreturn)
{
  // some geometry stuff to store pre-step instead of for each step.
  int pliter,pl;
  int jjdim;
  FTYPE dimfactU[NPR];
  PLOOP(pliter,pl) dimfactU[pl]=1.0; // default
  DLOOPA(jjdim) dimfactU[UU+jjdim]=dimfactU[URAD0+jjdim]=sqrt(fabs(ptrgeom->gcon[GIND(jjdim,jjdim)]));
  SLOOPA(jjdim) dimfactU[B1+jjdim-1] = 1.0/dimfactU[U1+jjdim-1];


  int i1,i2,i3,iv,ii,jj,kk,sc;
  FTYPE realdt;
  int gotbest,bestfailreturnf;
  FTYPE iJ[NPR][NPR];
  FTYPE iJ2[NPR][NPR];

  // store pb as might have didentropyalready=1 and then can use pborig[UU] as entropy's solution for u_g, and that can be used to create condition to avoid over-iterating with energy solver.
  FTYPE pborig[NPR];
  PLOOP(pliter,pl) pborig[pl]=pb[pl];

#if(DEBUGMAXITER)
  FTYPE pppreholdlist[IMPMAXITER+2][NPR]={{0}}; // for debug
  FTYPE ppposholdlist[IMPMAXITER+2][NPR]={{0}}; // for debug
  FTYPE f1reportlist[IMPMAXITER+2][NDIM]={{0}}; // for debug
  FTYPE f1list[IMPMAXITER+2][NDIM]={{0}}; // for debug
  FTYPE errorabsf1list[IMPMAXITER+2]={0}; // for debug
  int realiterlist[IMPMAXITER+2]={{-1}}; // for debug
  FTYPE jac00list[IMPMAXITER+2]={{BIG}}; // for debug
#else
  FTYPE (*pppreholdlist)[NPR];
  FTYPE (*ppposholdlist)[NPR];
  FTYPE (*f1reportlist)[NDIM];
  FTYPE (*f1list)[NDIM];
  FTYPE *errorabsf1list;
  int *realiterlist;
  FTYPE *jac00list;
#endif

  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uuppp[NPR],uu[NPR],uuporig[NPR],uu0orig[NPR],bestuu[NPR];
  FTYPE pp0[NPR],ppp[NPR],pppp[NPR],ppppp[NPR],pp[NPR],ppporig[NPR],pp0orig[NPR],bestpp[NPR];
  FTYPE f1[NPR],f1norm[NPR],f1report[NPR],f3report[NPR],lowestfreportf1[NPR],lowestfreportf3[NPR];
  FTYPE f1p[NPR];
  FTYPE pppriorsteptype[NPR],uupriorsteptype[NPR];

  FTYPE radsource[NPR], deltas[NPR]; 
  extern int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, FTYPE errorabs, FTYPE errorabsbestexternal, int iters, int iterstotal, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);
  int mathfailtype;


  int eomtypelocal=*eomtype; // default choice
  int allowlocalfailurefixandnoreport=0; // must be 0 so implicit method knows when really failure
  int failreturn=0;



  // DEBUG VARS
  int showmessages=0; // by default 0, don't show any messages for inversion stuff during implicit solver, unless debugging.  Assume any moment of inversion failure is corrected for now unless failure of final inversion done outside implicit solver.
  int showmessagesheavy=0;  // very detailed for common debugging
  static long long int failnum=0;
  
  int doingitsomecpu=0;
  int doingit=0;
#if(0)
  if(nstep==15 && steppart==0 && ptrgeom->i==6 && ptrgeom->j==8 && ptrgeom->k==0){
    doingitsomecpu=1;
    if(myid==0){ // so similar situation and grid at least
      dualfprintf(fail_file,"DOINGIT\n");
      doingit=1;

      ///////////
      // insert FAILRETURN data here.
      ///////////

      showmessages=showmessagesheavy=1;
    }// end on doing it core
  }

#endif





  ////////////////////
  //
  // initialize errors
  //
  ////////////////////
  PLOOP(pliter,pl){
    f1[pl]=BIG;
    f1p[pl]=BIG;
  }


  //////////////
  //
  // setup reversion to best solution for uu in case iterations lead to worse error and reach maximum iterations.
  //
  /////////////
  gotbest=0;
  bestfailreturnf=0;
  PLOOP(pliter,pl){
    lowestfreportf1[pl]=BIG;
    lowestfreportf3[pl]=BIG;
  }
  // setup locally-used ppfirst that can pass back as pb if good solution
  int gotfirstnofail=0;
  FTYPE ppfirst[NPR];



  ///////////////////
  //
  // setup implicit iteration procedure and loops
  //
  ///////////////////
  realdt = compute_dt(CUf,dt);
  FTYPE fracdtuu0=1.0,fracdtG=1.0,fracuup=1.0; // initially try full realstep step
  FTYPE fracdtuu0p=fracdtuu0;
  FTYPE fracdtuu0pp=fracdtuu0p;

  FTYPE fracdtGp=fracdtG;
  FTYPE fracdtGpp=fracdtGp;

  int numdims,startjac,endjac,implicititer,implicitferr,BEGINMOMSTEPS,ENDMOMSTEPS,BEGINENERGYSTEPS,ENDENERGYSTEPS,BEGINFULLSTEPS,ENDFULLSTEPS,BEGINNORMALSTEPS,irefU[NDIM],iotherU[NDIM],erefU[NDIM],eotherU[NDIM],signgd2,signgd4,signgd6,signgd7;
  int fakeitermethod=IMPMAXITER;// used by normal stepping, so just make maximum
  define_method(fakeitermethod, &eomtypelocal, itermode, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
  // no need to define numdims,jacs,refs, or signs yet.



  /////////////////
  //
  // set uu0 = "initial+flux" contribution to uu
  //
  ////////////////
  PLOOP(pliter,pl) uu[pl]=uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);




  /////////////////////
  //
  // get default failure state from Uiin
  //
  //////////////////////
  FTYPE prtestUiin[NPR];
  PLOOP(pliter,pl) prtestUiin[pl]=piin[pl]; // initial guess (should be easy with piin=ppin(Uiin))

  //#define GETDEFAULTFAILURESTATE 1 // 1 was normal before.
#define GETDEFAULTFAILURESTATE (IMPPMHDTYPE(implicititer)==0 && 0) // with new rad inversion scheme, seems to push through without issue to just let hit gamma ceiling temporarily -- even with RADPULSEPLANAR
  // Otherwise, takes *many* iterations, not just f1iters, to get solution.  Even though that was old code state that was working, not required now.


  // default is to allow no failure unless iterating MHD primitives in which case a radiation failure is ok.
  int failreturnallowable;

  // no need to worry about RAD failure if not iterating rad quantities
  // As stated above, with new rad inversion, seems ok to temporarily hit ceiling.
  // With tests like RADTUBE, can't allow radiative inversion cieling else dies.  While with tests like RADPULSE, fastest to converge to good solution with allowing hitting the ceiling.  So mixed issue.
  if(IMPMHDTYPE(implicititer)) failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
  else failreturnallowable=UTOPRIMNOFAIL;
  //  if(IMPMHDTYPE(implicititer)) failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
  //  else failreturnallowable=UTOPRIMNOFAIL;

  int failreturnallowableuse=failreturnallowable;
  int failreturnallowablefirst=-1;

  if(GETDEFAULTFAILURESTATE && IMPPMHDTYPE(implicititer)==0){    // No need for failure state if doing MHD iteration using primitives
    // KORALTODO: Instead of doing this, just keep inversion error from previously as stored in pflag.  So don't reset pflag elsewhere in fixup_utoprim() for elsewhere
    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)
    struct of_newtonstats newtonstats;
    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    int finalstep = 1;
    int doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    failreturnallowable=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Uiin, q, ptrgeom, dissmeasure, prtestUiin, &newtonstats);
    if(failreturnallowable!=UTOPRIMGENWRAPPERRETURNNOFAIL){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is already a problem with %d\n",failreturnallowable);
    }
    else{
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is NOT a problem with %d\n",failreturnallowable);
    }

    if(failreturnallowable==UTOPRIMGENWRAPPERRETURNNOFAIL) gotfirstnofail=1;
    else gotfirstnofail=0;
    failreturnallowablefirst=failreturnallowable;

    //    if(myid==5 && nstep==1 && steppart==0 && ptrgeom->i==19 && ptrgeom->j==15){
    //      PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d prtestUiin=%21.15g Uiin=%21.15g\n",pl,prtestUiin[pl],Uiin[pl]);
    //    }

  }
  else{
    failreturnallowablefirst=failreturnallowable;
    gotfirstnofail=1;//default is no failure
  }






  /////////////////////
  //
  // see if uu0->p possible and desired
  //
  //////////////////////
  FTYPE prtestUU0[NPR];
  PLOOP(pliter,pl) prtestUU0[pl]=pb[pl]; // initial guess

  // KORALTODO: Need to check if UFSET with no dUother fails.  How it fails, must allow since can do nothing better.  This avoids excessive attempts to get good solution without that failure!  Should speed-up things.  But what about error recovery?  If goes from CASE to no case!
  if(USEDUINRADUPDATE==2){
    // uu0 will hold original vector of conserved
    // here original means U[before fluxes, geometry, etc.] + dU[due to fluxes, geometry, etc. already applied and included in dUother]
    // This is required for stiff source term so immediately have balance between fluxes+geometry+radiation.  Otherwise, radiation diffuses.
    // I'm guessing that even though one uses RK2 or RK3, the first step generates large radiative velocities without any balanced source term since U isn't updated yet.  One would hope RK2 would recover on the final substep, but it doesn't!  In RK2, upon the final substep, that velocity is present for the radiation source term.  But it's also present for the fluxes!  That is, if there were no flux update on the final substep, then the source would have balanced the previous flux, but yet another flux is done, so there can be no balance.  This leads to a run-away velocity that would be similar to the \tau\sim 1 case.
    // NOTE: If this gives radiation or mhd failure, then less likely that will be actual solution since not even original uu0 has inversion.
    // Note that "q" isn't used in this function or used in function call, so don't have to update it here.

    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)
    struct of_newtonstats newtonstats;
    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    int failreturninversion;
    int finalstep = 1;
    int doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    failreturninversion=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, uu0, q, ptrgeom, dissmeasure, prtestUU0, &newtonstats);

    // get pp0(uu0) so uu,uu0 properly associated with pp,pp0
    // set first pp that is p(uu0), assuming that can do inversion
    PLOOP(pliter,pl) ppfirst[pl]=pp[pl]=pp0[pl]=prtestUU0[pl];

  }
  else if(USEDUINRADUPDATE==1){
    // normal behavior, where uu0 has no known primitive inversion yet.
    // if entropy got solution and then energy is being done here, then pb is good entropy guess.
    PLOOP(pliter,pl){
      uu[pl] = uu0[pl];
      ppfirst[pl] = pp[pl] = pp0[pl] = pb[pl];
    }
  }
  else{
    // then (not recommended) just using Uiin as uu0 .  KORALNOTE: But uu0 reset later anyways.
    PLOOP(pliter,pl){
      uu[pl] = uu0[pl] = Uiin[pl];
      ppfirst[pl] = pp[pl] = pp0[pl] = piin[pl];
    }
  }  



  /////////////////////////////
  //
  // Fix u_g so entropy not smaller than guess
  // Only do this if u_g is just guess and is being solved for.
  // When eomtype==EOMCOLDGRMHD, then u_g is static so shouldn't modify as if was guess.
  //
  /////////////////////////////

  if(*eomtype!=EOMCOLDGRMHD && IMPMHDTYPE(implicititer)==1){
    if(ENTROPYFIXGUESS){
      entropyfixguess(q, ptrgeom, uu0, pp);
      // piin is sometimes even a bit higher, and want to start high, and helps to avoid lack of convergence issue.
      pp[UU]=MAX(pp[UU],piin[UU]);
    }
    
    if(modprim==1){
      pp[UU]=10.0*MAX(pp[UU],piin[UU]); // raise u_g a bit.
    }
  }
  





  /////////////////////////////
  //
  // SETUP Damping loop
  //
  /////////////////////////////

  FTYPE ppdampbackup[NPR],uudampbackup[NPR],uu0dampbackup[NPR];
  // backup those things that below can change, yet we want to start fresh each damp attempt (except best and backup stuff values and errors)
  PLOOP(pliter,pl){
    ppdampbackup[pl]=pp[pl];
    uudampbackup[pl]=uu[pl];
    uu0dampbackup[pl]=uu0[pl];
  }
  // backups and bests
  int gotbackup=0;
  int totaliters=0;
  int iter=0;
  int debugiter=0;
  int debugiteratteempts[NUMDAMPATTEMPTS];
  int momiters=0,energyiters=0,fulliters=0;
  FTYPE errorabsf1=BIG;
  FTYPE suberrorabsf1=BIG;
  FTYPE errorabsf3=BIG;
  FTYPE suberrorabsf3=BIG;
  
  FTYPE trueimptryconv=IMPTRYCONV;
  FTYPE trueimptryconvabs=IMPTRYCONVABS;
  FTYPE trueimptryconvalt=IMPTRYCONVALT;

  // best is over all damps as well
  FTYPE errorabsbest=BIG;
  int failreturnbest=FAILRETURNGENERAL;

  ////////////////
  //
  // THE DAMP LOOP ITSELF
  //
  ////////////////
  int dampattempt;
  for(dampattempt=0;dampattempt<truenumdampattempts;dampattempt++){
    FTYPE DAMPFACTOR0;
    // dampattempt>0 refers to any attempt beyond the very first.  Uses iter from end of region inside this loop.
    if(dampattempt>0 && (errorabsf1<IMPTRYDAMPCONV && ACCEPTASNOFAILURE(failreturn) || failreturn==FAILRETURNMODESWITCH) ){ // try damping if any bad failure or not desired tolerance when damping
      if(dampattempt>=2 && failreturn!=FAILRETURNMODESWITCH){ // dampattempt>=2 refers to attempts with at least 1 damp attempt
        if(debugfail>=2) dualfprintf(fail_file,"Damping worked to reach desired tolerance: errorabsf1=%g (IMPTRYDAMPCONV=%g), so should have lower error: dampattempt=%d iter=%d\n",errorabsf1,IMPTRYDAMPCONV,dampattempt,iter);
      }
      break; // if didn't hit problem, so no need to damp since got tolerance requested or returned because will just switch to another scheme.
    }
    else{
      // control factor by which step Newton's method.
      DAMPFACTOR0=1.0/pow(2.0,(FTYPE)(dampattempt));
      if(dampattempt>0) if(debugfail>=2) dualfprintf(fail_file,"Trying dampattempt=%d DAMPFACTOR0=%g failreturn=%d errorabsf1=%g iter=%d\n",dampattempt,DAMPFACTOR0,failreturn,errorabsf1,iter);

      // start fresh
      iter=debugiter=0;
      failreturn=FAILRETURNNOFAIL; // default is no failure
      mathfailtype=-1; // indicates not set
      PLOOP(pliter,pl){
        uu0[pl]=uu0dampbackup[pl];
      }
      if(errorabsf1<TRYHARDERFEEDGUESSTOL){
        // then keep pp and uu as better starting point
      }
      else{
        PLOOP(pliter,pl){
          pp[pl]=ppdampbackup[pl];
          uu[pl]=uudampbackup[pl];
        }
      }

    }



    ////////////////////////////////
    //
    // SETUP IMPLICIT ITERATIONS
    //
    ////////////////////////////////
    int f1iter;
    int checkconv,changeotherdt;
    FTYPE impepsjac=IMPEPSLARGE; // default
    FTYPE errorabspf1[NUMPRIORERRORS]={BIG};
    FTYPE errorabspf3=BIG;

    // initialize previous 'good inversion' based uu's
    PLOOP(pliter,pl){
      uupriorsteptype[pl]=uuppp[pl]=uupp[pl]=uuporig[pl]=uup[pl]=uu0orig[pl]=uu[pl];
      pppriorsteptype[pl]=ppppp[pl]=pppp[pl]=ppporig[pl]=ppp[pl]=pp0orig[pl]=pp[pl];
    }


    // setup debug so can see starting guess
    if(DEBUGMAXITER&& dampattempt==0){
      int iterlist=0;
      PLOOP(pliter,pl) pppreholdlist[iterlist][pl]=pp[pl];
      PLOOP(pliter,pl) ppposholdlist[iterlist][pl]=pp[pl];
      if(DEBUGMAXITERVELOCITY==1){
        SLOOPA(jj){
          pppreholdlist[iterlist][U1+jj-1]=q->ucon[jj];
          pppreholdlist[iterlist][URAD1+jj-1]=q->uradcon[jj];
          ppposholdlist[iterlist][U1+jj-1]=q->ucon[jj];
          ppposholdlist[iterlist][URAD1+jj-1]=q->uradcon[jj];
        }
      }
      DLOOPA(jj) f1reportlist[iterlist][jj]=BIG;
      DLOOPA(jj) f1list[iterlist][jj]=BIG;
      errorabsf1list[iterlist]=BIG;
      realiterlist[iterlist]=-1;
      jac00list[iterlist]=BIG;
    }
    // DEBUG:
    *errorabsreturn=BIG;

    // whether holding as positive and outher counts
    int holdingaspositive=0,iterhold=0;
    int countholdpositive=0;
    int countbadenergy=0;
    int counterrorrose=0;
    int priorerrorscount=0;
    int canbreak=0;
    int notfinite=0;
    int convreturnf3limit=0;
    int notholding=1;
    FTYPE DAMPFACTOR;
    int earlylowerror;

    ////////////////////////////////
    //
    // IMPLICIT LOOP ITSELF
    //
    ////////////////////////////////

    do{
      debugiter++; // never skips, goes every step
      debugiteratteempts[dampattempt]=debugiter;
      iter++; // below might skip some iter, used to control which equations used
      if(iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
        momiters++;
      }
      if(iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
        energyiters++;
      }
      if(iter>=BEGINFULLSTEPS && iter<=ENDFULLSTEPS){
        fulliters++;
      }


      // setup method and signs
      define_method(iter, &eomtypelocal, itermode, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
      get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);


      ////////////////////
      //
      // things to reset each iteration start for any steps
      //
      ///////////////////
      canbreak=0; // reset each start of iteration
      convreturnf3limit=0;
      notholding=1; // default is no hold
      earlylowerror=0;

    
      if(iter>10){ // KORALTODO: improve upon this later.  Only matters if not doing PMHD method
        // assume trying hard and failing to work, then allow CASE radiation errors
        failreturnallowable=failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
      }

      if(CHANGEDAMPFACTOR&&trueimpmaxiter==IMPMAXITERQUICK && iter>IMPMAXITERQUICK/2){
        DAMPFACTOR=0.5*DAMPFACTOR0;
      }
      else if(CHANGEDAMPFACTOR&&trueimpmaxiter==IMPMAXITER && iter>MIN(IMPMAXITER/2,20)){
        DAMPFACTOR=0.5*DAMPFACTOR0;
      }
      else DAMPFACTOR=DAMPFACTOR0;


      ///////////
      //
      // vector of conserved, primitives, and fractions of steps used at the previous few iterations
      //
      //////////
      PLOOP(pliter,pl){
        uuppp[pl]=uupp[pl]; // uuppp...
        uupp[pl]=uup[pl]; // uupp will have solution for inversion: P(uupp)
        uup[pl]=uu[pl]; // uup will not necessarily have P(uup) because uu used Newton step.
        uuporig[pl]=uu[pl];
        
        ppppp[pl]=pppp[pl]; // ppppp will have knowledge of 2 prior errorabsf1's
        pppp[pl]=ppp[pl]; // pppp will have solution for inversion
        ppp[pl]=pp[pl]; // ppp will not necessarily have solution because pp used Newton step.
        ppporig[pl]=pp[pl];
      }

      fracdtuu0pp=fracdtuu0p; // fracdtuu0 used when computing previous f1 and f2's
      fracdtuu0p=fracdtuu0; // fracdtuu0 used when computing previous f1 and f2's

      fracdtGpp=fracdtGp; // fracdtG used when computing previous f1 and f2's
      fracdtGp=fracdtG; // fracdtG used when computing previous f1 and f2's

      errorabspf3=errorabsf3;


      /////////////
      //
      // reset at start of new use of equations
      //
      ////////////
      int itererror;
      if(iter==BEGINMOMSTEPS || iter==BEGINENERGYSTEPS || iter==BEGINFULLSTEPS){
        iterhold=0;
        countbadenergy=0;
        counterrorrose=0;
        priorerrorscount=0;
        holdingaspositive=0;
        countholdpositive=0;
        for(itererror=0;itererror<NUMPRIORERRORS;itererror++) errorabspf1[itererror]=BIG;
        // store uu,pp before modified using new step
        PLOOP(pliter,pl){
          uupriorsteptype[pl]=uu[pl];
          pppriorsteptype[pl]=pp[pl];
        }
      }
      else if(iter>1){
        for(itererror=MIN(priorerrorscount,NUMPRIORERRORS)-1;itererror>=1;itererror--) errorabspf1[itererror]=errorabspf1[itererror-1]; // shift to higher
        errorabspf1[0]=suberrorabsf1; // only for equations that are currently iterating.
      }
      // store previous f1 values
      PLOOP(pliter,pl) f1p[pl]=f1[pl];




      // KORALTODO: Once use certain uu0 value and succeed in getting f1, not enough.  But if take a step and *then* good f1, then know that original U was good choice and can restart from that point instead of backing up uu0 again.
      // KORALTODO: Look at whether bounding error by sign as in bisection (only true in 1D!!), and if approaching 0 slowly enough and still hit CASE, then must be real CASE not a stepping issue.
      // KORALTODO: Getting f1 is just about f1(U) as far as radiation is concerned.  So as far as CASE issues, getting f1(U) means we are good with that U and we can certainly stick with the used uu0.




      /////////////////
      //
      // get error function (f1) and inversion (uu->pp) using uu
      //
      /////////////////
      int failreturnferr;
      for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
        int whichcall=1;
      
        eomtypelocal=*eomtype; // re-chose default each time.  If this reduces to a new eomtype, then Jacobian will stick with that for consistency!
        int goexplicit=0;
        failreturnferr=f_implicit(iter, failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, radinvmod, pp, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, &goexplicit); // modifies uu and pp

        // see if 4-force negligible
        if(goexplicit){
          // immediate return.
          dualfprintf(fail_file,"SHIT: failreturnferr=%d goexplicit=%d\n",failreturnferr,goexplicit);
          return(-1);
        }

        // f1 error calculation updated non-iterated pp or uu, so store.  Ok to re-store iterated uu and pp as well
        // KORALTODO: This might mess up when f1 backs-up uu and pp
        PLOOP(pliter,pl){
          uup[pl]=uu[pl];
          ppp[pl]=pp[pl];
        }


#if(0)
        if(pp[PRAD0]<10.0*ERADLIMIT){
          // try smaller tolerance
          trueimptryconv=10.0*NUMEPSILON;
          trueimptryconvabs=((FTYPE)(NDIM+2)*trueimptryconv);
          trueimptryconvalt=(MAX(1E-8,trueimptryconvabs));
        }
#endif


        if(failreturnferr){

          //          if(debugfail>=2) dualfprintf(fail_file,"Failed f1: %d\n",failreturnferr);

#define BACKUPRELEASEFAIL0 (1E-4)
#define BACKUPRELEASEFAIL1 (1E-6)
#define BACKUPRELEASEFAIL2 (1E-7)
#define MAXF1ITERRADTRY (10) // number of iterations after which to release and allow radiation to "fail"

          // if backing up alot, then allow same failure as original Uiin in hopes that can recover that way (otherwise would have hoped would recover via flux update to Uiin)
          if(fracdtuu0<BACKUPRELEASEFAIL0){
            failreturnallowableuse=failreturnallowable;
          }
          if(fracdtuu0<BACKUPRELEASEFAIL1 || f1iter>=MAXF1ITERRADTRY){
            failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
          }


          if(iter==1){
            // if initial uu failed, then should take smaller jump from Uiin->uu until settled between fluid and radiation.
            // If here, know original Uiin is good, so take baby steps in case uu needs heavy raditive changes.
            // if f1 fails, try going back to Uiin a bit
            fracdtuu0*=RADDAMPDELTA; // DAMP Uiin->uu0 step that may be too large and generated too large G
            fracdtuu0=MIN(1.0,fracdtuu0);

            //          if(fracdtuu0<BACKUPRELEASEFAIL2) fracdtuu0=0.0; // just stop trying to start with some uu0 and revert to Uiin

            // modifies uup (ppp) and uu (pp) as if starting over, and then another call to f_implicit(f1) will change uu by generally smaller amount.
            PLOOP(pliter,pl){
              uup[pl]=uu[pl]=uu0[pl];
              ppp[pl]=pp[pl]=pp0[pl];
            }
          }
          else{
            // if here, then assume prior uup was good in sense that no failures like P(uup) is good inversion.  And uu=uup before f_implicit(f1) is called.
            // No, not necessarily true, because uu updated with some-sized Newton step without checking if inversion is good for that uu.
            // So need to damp between original uup and uupp .  This essentially damps Newton step now that have knowledge the Newton step was too large as based upon P(U) failure.

            // Avoid G-damping because not needed so far.  If added, competes in non-trivial way with fracuup damping that's required separately for large forces in some problems beyond iter=1 (e.g. NTUBE=31).
            //          fracdtG*=RADDAMPDELTA; // DAMP give only fraction of 4-force to let uu to catch-up
            //           fracdtG=0.5; // DAMP give only fraction of 4-force to let uu to catch-up

            fracuup*=RADDAMPDELTA; // DAMP in case Newton step is too large after iter>1 and stuck with certain uu from Newton step that gets stored in uup above.

            PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uuporig[pl];
            //          PLOOP(pliter,pl) uu[pl]=(1.0-fracuup)*uupp[pl] + fracuup*uup[pl];
            PLOOP(pliter,pl) uup[pl]=uu[pl]; // store new version of prior Newton step

            PLOOP(pliter,pl) pp[pl]=(1.0-fracuup)*pppp[pl] + fracuup*ppporig[pl];
            //          PLOOP(pliter,pl) pp[pl]=(1.0-fracuup)*pppp[pl] + fracuup*ppp[pl];
            PLOOP(pliter,pl) ppp[pl]=pp[pl]; // store new version of prior Newton step
          
            // get interpolated fracdtuu0 so using fracdtuu0 that was used with the corresponding uu
            fracdtuu0=(1.0-fracuup)*fracdtuu0pp + fracuup*fracdtuu0p;
            // same for fracdtG
            fracdtG=(1.0-fracuup)*fracdtGpp + fracuup*fracdtGp;

          }


          // get uu0 (which may be changing)
          PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
          // KORALNOTE: pp0 not used, setting uu0 above fixes error function where only uu0 is needed.



          // keep below so can count inversion failures against retry successes in the failure file.
          if(showmessages && debugfail>=2) dualfprintf(fail_file,"f_implicit for f1 failed: iter=%d  Backing-up both uu0 and G.: f1iter=%d fracdtuu0=%g fracdtG=%g fracuup=%g\n",iter,f1iter,fracdtuu0,fracdtG,fracuup);
          if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Ui=%26.20g uu0=%26.20g uu0orig=%26.20g uu=%26.20g uup=%26.20g dUother=%26.20g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
        }// end if failed inversion in f_implicit()
        else{
          // then success, so was able to do inversion P(U) with change in radiation on fluid: P(uu[fluid] = uu0[fluid] - (uu[rad]-uu0[rad]))
          // This doesn't necessarily mean could do P(uu0) except for iter=1.
          // Corresponds to success for a certain P(uu0,uu) pair.
          if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"SUCCESS: pl=%d Ui=%26.20g uu0=%26.20g uu0orig=%26.20g uu=%26.20g uup=%26.20g dUother=%26.20g: pp(pnew)=%g : fracdtuu0=%g fracdtG=%g fracuup=%g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl],pp[pl],fracdtuu0,fracdtG,fracuup);
          break;
        }

        // if during f1iter-ations u_g or rho is negative, switch to entropy
#define F1ITERSWITCHONNEG 20
        if(havebackup){
          if(f1iter>F1ITERSWITCHONNEG && (pp[RHO]<0 || pp[UU]<0)){
            failreturn=FAILRETURNMODESWITCH; mathfailtype=70;
            if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Switched modes during f1iter=%d : rho=%21.15g ug=%21.15g\n",f1iter,pp[RHO],pp[UU]);
            break;
          }
        }
      

      }// end loop over f1iter
      *f1itersreturn += f1iter;


      //////////////
      // break again out of total loop if broke in f1iter loop
      if(failreturn){
        if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Breaking out of loop as think f1iter wanted us to.\n");
        break;
      }
      else{// else if good f1

        /////////////////////////////////
        //
        // Check if reached max f1 iterations
        //
        ////////////////////////////////
        if(f1iter==MAXF1TRIES){

          if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES: fracdtuu0=%g nstep=%ld steppart=%d ijk=%d %d %d : iter=%d eomtype=%d failreturn=%d\n",fracdtuu0,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter,eomtypelocal,failreturn);
          if(havebackup){
            failreturn=FAILRETURNMODESWITCH; mathfailtype=20;
            if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"SWITCHING MODE: Detected MAXF1TRIES\n");
            break;
          }
          else{
            failreturn=FAILRETURNGENERAL; mathfailtype=2;
            if(doingit==1) myexit(10000000); // DEBUG
            // Note that if inversion reduces to entropy or cold, don't fail, so passes until reached this point.  But convergence can be hard if flipping around which EOMs for the inversion are actually used.
            break;
          }
        }
        else{
          // restore fracuup back to 1 since this is only meant to adjust how much go back to previous uu to be able to get f1 computed.
          // fracuup doesn't stay <1.0 because each attempt to get f1 is independent.
          // KORALNOTE: Perhaps reasonable and safer to keep fracuup as not reverted back to 1, since apparently unable to take full steps.  This effectively damps stepping.
          fracuup=1.0;
        }
  
        // diagnose
        if(showmessagesheavy) JACLOOPALT(ii,startjac,endjac) dualfprintf(fail_file,"i=%d ii=%d erefU[ii]=%d f1=%g\n",ptrgeom->i,ii,erefU[ii],f1[erefU[ii]]);
      }// else if f1 calculation didn't fail.



    



      /////////
      //
      // see if should check convergence or check how solution is behaving.
      //
      /////////
      checkconv=1;
      changeotherdt=1;
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracuup!=1.0){
      //      if(fabs(fracuup-1.0)>10.0*NUMEPSILON){
      if(fracuup<1.0){
        // try increasing amount of uu or pp used
        fracuup*=RADDAMPUNDELTA;
        fracuup=MIN(1.0,fracuup);
        checkconv=0;
        changeotherdt=0; // ensure fracuup back to 1.0 first before reverting others.
      }
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracdtuu0!=1.0){
      //      if(fabs(fracdtuu0-1.0)>10.0*NUMEPSILON && changeotherdt){
      if(fracdtuu0<1.0 && changeotherdt){
        // try increasing uu0 away from Uiin to account for full dUother
        fracdtuu0*=RADDAMPUNDELTA;
        fracdtuu0=MIN(1.0,fracdtuu0);
        PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0
        // KORALNOTE: No need to get pp0, since never used.  uu0 only used in error function.
        checkconv=0;
      }
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracdtG!=1.0){
      //      if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
      if(fracdtG<1.0 && changeotherdt){
        // try increasing amount of G applied
        fracdtG*=RADDAMPUNDELTA;
        fracdtG=MIN(1.0,fracdtG);
        checkconv=0;
      }      
      ////////////////////////////////////////////////////////////////////////////
      if(iter<BEGINNORMALSTEPS) checkconv=0; // don't check actual convergence till doing full steps
      ////////////////////////////////////////////////////////////////////////////





      //////////////
      //
      // get error using f1 and f1norm
      //
      //////////////
      int convreturnf1=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs, realdt, DIMTYPEFCONS,eomtypelocal ,itermode,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1);
      // but don't break, since need to iterate a bit first and check |dU/U| and need to see if checkconv==1
      suberrorabsf1=0.0;  JACLOOPSUBERROR(jj,startjac,endjac) suberrorabsf1     += fabs(f1report[erefU[jj]]); // e.g. may only be energy error or only momentum error.



      // DEBUG STUFF
      if(DEBUGMAXITER&& dampattempt==0){
        PLOOP(pliter,pl) pppreholdlist[debugiter][pl]=pp[pl]; // just default dummy value in case break
        PLOOP(pliter,pl) ppposholdlist[debugiter][pl]=pp[pl]; // just default dummy value in case break
        if(DEBUGMAXITERVELOCITY==1){
          SLOOPA(jj){
            pppreholdlist[debugiter][U1+jj-1]=q->ucon[jj];
            pppreholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
            ppposholdlist[debugiter][U1+jj-1]=q->ucon[jj];
            ppposholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
          }
        }
        JACLOOPALT(jj,startjac,endjac) f1reportlist[debugiter][jj]=f1report[erefU[jj]];
        JACLOOPALT(jj,startjac,endjac) f1list[debugiter][jj]=f1[erefU[jj]];
        errorabsf1list[debugiter]=errorabsf1;
        realiterlist[debugiter]=iter;
        jac00list[debugiter]=BIG; // just default dummy value in case break
      }



      ////////////////
      //
      // If error dropped below tolerance for this sub-matrix iteration mode, then continue to next level.
      //
      /////////////////
      // check if energy only iteration has error that has dropped below tolerance, then can move on to 
      if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
        if(fabs(f1report[irefU[1]])<trueimptryconv && fabs(f1report[irefU[2]])<trueimptryconv && fabs(f1report[irefU[3]])<trueimptryconv){
          if(iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1; continue;} // force as if already doing energy steps.  If already next iteration is to be this energy step, then no skipping needed.
        }
      }
      // check if energy only iteration has error that has dropped below tolerance, then can move on to 
      if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
        //        if(f1report[irefU[0]]==BIG){ dualfprintf(fail_file,"FUDGE\n"); }
        // SUPERGODMARK: valgrind says belw is undefined, but don't see it.
        if(fabs(f1report[irefU[0]])<trueimptryconv){
          if(iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1; continue;} // force as if already doing normal steps.  If already next iteration is to be normal step, no need to skip.
        }
      }







      ///////////////////////////
      //
      //  PRE NEWTON ADJUSTMENTS
      //
      ///////////////////////////

      // check if doing energy stepping and error jumped up too much
      if(iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS || iter>=BEGINFULLSTEPS && iter<=ENDFULLSTEPS){// now all steps beyond energy
        if(fabs(f1[erefU[0]]/f1p[erefU[0]])>FACTORBADJUMPERROR && fabs(f1report[erefU[0]])>trueimptryconv){
          // then pseudo-bisect (between zero and previous ok error case
          if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"pseudo-bisect: iter=%d f1=%g f1p=%g pp=%g ppp=%g  pppp=%g  ppppp=%g\n",iter,f1[erefU[0]],f1p[erefU[0]],pp[erefU[0]],ppp[erefU[0]],pppp[erefU[0]],ppppp[erefU[0]]);
          //          pp[erefU[0]] = ppp[erefU[0]] = pppp[erefU[0]] = 0.5*(fabs(ppppp[erefU[0]]));
          pp[erefU[0]] = ppp[erefU[0]] = 0.5*fabs(pppp[erefU[0]]);
          // update debug with modifications
          if(DEBUGMAXITER&& dampattempt==0){
            PLOOP(pliter,pl) pppreholdlist[debugiter][pl]=pp[pl];
            PLOOP(pliter,pl) ppposholdlist[debugiter][pl]=pp[pl];
            if(DEBUGMAXITERVELOCITY==1){
              SLOOPA(jj){
                pppreholdlist[debugiter][U1+jj-1]=q->ucon[jj];
                pppreholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
                ppposholdlist[debugiter][U1+jj-1]=q->ucon[jj];
                ppposholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
              }
            }
            jac00list[debugiter]=iJ[irefU[0]][erefU[0]];
          }
          // need to get new error function so can take step based upon this as reference!
          continue; // head to start of loop to iter++ and get new error function.
        }
      }





      // only check convergence or check the properties of the solution if checkconv==1
      if(checkconv){
        /////////////////
        //
        // see if pre-convergence (happens if force is small or no force at all.  Can't necessarily continue since Jacobian can require arbitrarily large dU on fluid and fail to invert even if no fluid-radiation interaction!
        //test pre-convergence using initial |dU/U|
        // KORALTODO: This isn't a completely general error check since force might be large for fluid that needs itself to have more accuracy, but if using ~NUMEPSILON, won't resolve 4-force of radiation on fluid to better than that.
        //
        /////////////////
        FTYPE LOCALPREIMPCONV=MIN(10.0*NUMEPSILON,trueimptryconv); // more strict than later tolerance
        FTYPE LOCALPREIMPCONVABS=(FTYPE)(NDIM+2)*LOCALPREIMPCONV; // more strict than later tolerance
        if(f_error_check(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,LOCALPREIMPCONVABS,realdt,DIMTYPEFCONS,eomtypelocal,itermode,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin, uu0,uu,ptrgeom,&errorabsf1)){
          earlylowerror=1;
        }

        //////////////
        //
        // try to get best solution (should be based upon immediate f_error_check(uu,pp,errorabsf1)
        //
        // store error and solution in case eventually lead to max iterations and actually get worse error
        // f1 based
        //
        //////////////
        errorabsbest=0.0; JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsbest += fabs(lowestfreportf1[erefU[jj]]);
        if(errorabsbest>errorabsf1 && isfinite(errorabsf1) && (itermode==ITERMODECOLD || pp[RHO]>0.0 && pp[UU]>0.0 && pp[PRAD0]>0.0)){
          PLOOP(pliter,pl) bestuu[pl]=uu[pl];
          PLOOP(pliter,pl) bestpp[pl]=pp[pl];
          errorabsbest=errorabsf1;
          JACLOOPALT(jj,startjac,endjac) lowestfreportf1[erefU[jj]]=f1report[erefU[jj]];
          gotbest=1;
        }


        if(earlylowerror){
          if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Early low error=%g iter=%d\n",errorabsf1,iter);
          //  not failure.
          break;
        }


        // check if error repeatedly rises
        // do this even if damping, since damping should only help get continuous reduction in error.
        if(NUMNOERRORREDUCE && iter>=BEGINNORMALSTEPS){
          if(iter>NUMNOERRORREDUCE0 && suberrorabsf1>trueimptryconvabs){ // no need to do this if actually error is below desired tolerance,  hence second argument
            if(suberrorabsf1>=errorabspf1[0]) counterrorrose++;
            int allowedtoabort=(havebackup || ABORTBACKUPIFNOERRORREDUCE==1 && suberrorabsf1<trueimptryconvalt);
            if(counterrorrose>=NUMNOERRORREDUCE && allowedtoabort){ // would be risky to do abort if don't have backup, even if enter into limit cycle.

              if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
                if(iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1; continue;} // force as if already doing energy steps.  If already next iteration is to be this energy step, then no skipping needed.
              }// end if doing momentum steps and error not decreasing fast enough, then skip to energy steps
              else if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
                if(iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1; continue;} // force as if already doing normal steps.  If already next iteration is to be normal step, no need to skip.
              }// end if doing energy steps and error not decreasing fast enough, then skip to normal steps
              else{// else if normal step

                if(havebackup){
                  failreturn=FAILRETURNMODESWITCH; mathfailtype=80;
                  if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"SWITCHING MODE: Detected did not decrease error %d times at iter=%d : suberrorabsf1=%g errorabspf1=%g\n",counterrorrose,iter,suberrorabsf1,errorabspf1[0]);
                  // if want to ensure should have gotten solution, should still report
                }
                else{
                  // then aborting due to error alone even without backup
                  canbreak=2;
                  mathfailtype=81;
                  if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"Aborting even without backup because error oscillated (iter=%d) and suberrorabsf1=%g errorabsf1=%g\n",iter,suberrorabsf1,errorabsf1);
                  // no failure or switch.
                }
                break;
              }
            }// end if normal step
          }
        }      


        // check if error isn't decreasing enough for be interesting
        // but don't check on error if holding on u_g>0 since error can be bad until settle to near root.
        // and don't do this check if damping, since if damping really want to try harder.
        if(NUMPRIORERRORS>0 && holdingaspositive==0 && dampattempt==0){
          if(priorerrorscount>=NUMPRIORERRORSITER0){
            FTYPE erroraverage=0.0; int numerroraverage=0;
            for(itererror=1;itererror<MIN(priorerrorscount,NUMPRIORERRORS);itererror++){ erroraverage += errorabspf1[itererror]; numerroraverage++;} erroraverage/=((FTYPE)numerroraverage);
            FTYPE changerequired;
            if(dampattempt==0) changerequired=PRIORERRORCHANGEREQUIRED;
            else if(dampattempt==1) changerequired=0.7;
            else if(dampattempt==2) changerequired=0.8;
            else if(dampattempt==3) changerequired=0.9;
            else changerequired=0.95;
            if(suberrorabsf1/erroraverage>changerequired && suberrorabsf1>trueimptryconvabs){

              if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
                if(iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1; continue;} // force as if already doing energy steps.  If already next iteration is to be this energy step, then no skipping needed.
              }// end if doing momentum steps and error not decreasing fast enough, then skip to energy steps
              else if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
                if(iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1; continue;} // force as if already doing normal steps.  If already next iteration is to be normal step, no need to skip.
              }// end if doing energy steps and error not decreasing fast enough, then skip to normal steps
              else{
                canbreak=3;
                if(debugfail>=DEBUGLEVELIMPSOLVER){
                  dualfprintf(fail_file,"Error is not decreasing sufficiently fast: iter=%d priorerrorscount=%d suberrorabsf1=%g\n",iter,priorerrorscount,suberrorabsf1);
                  for(itererror=1;itererror<MIN(priorerrorscount,NUMPRIORERRORS);itererror++) dualfprintf(fail_file,"errorabspf1[%d]=%g\n",itererror,errorabspf1[itererror]);
                }
                break;
              }// end if doing normal steps and error not decreasing fast enough
            }// end if error not decreasing enough and also higher than desired tolerance in error
          }// end if enough prior errors to make average measurement of past error
          priorerrorscount++;
        }





      }// end if checkconv==1




      /////////////////////////
      //
      // See if solution has nan'ed or inf'ed out.
      //
      ////////////////////////
      if(IMPUTYPE(implicititer)){
        notfinite = !isfinite(uu[irefU[0]])|| !isfinite(uu[irefU[1]])|| !isfinite(uu[irefU[2]])|| !isfinite(uu[irefU[3]]) || !isfinite(uup[irefU[0]])|| !isfinite(uup[irefU[1]])|| !isfinite(uup[irefU[2]])|| !isfinite(uup[irefU[3]]);
      }
      else if(IMPPTYPE(implicititer)){
        notfinite = !isfinite(pp[irefU[0]])|| !isfinite(pp[irefU[1]])|| !isfinite(pp[irefU[2]])|| !isfinite(pp[irefU[3]]) || !isfinite(ppp[irefU[0]])|| !isfinite(ppp[irefU[1]])|| !isfinite(ppp[irefU[2]])|| !isfinite(ppp[irefU[3]]);
      }


      /////////////////
      //
      // continue with computing Jacobian and Newton step if origin point for error function didn't nan'out or inf'out.
      //
      /////////////////
      if(!notfinite){
    
        /////////
        //
        // get Jacobian and inversion Jacobian 
        //
        /////////
        //      eomtypelocal=*eomtype; // re-chose default each time.  No, stick with what f1 reduced to for consistency.
        //        dualfprintf(fail_file,"iJ call: iter=%d\n",iter);

        // assume as error gets small, function becomes linear and can use smaller delta for Jacobian
        if(errorabsf1<ERRORFORIMPEPSSMALL) impepsjac=IMPEPSSMALL;
        else impepsjac=IMPEPSLARGE;
        //        int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, impepsjac, iter, errorabsf1, dimfactU, uu, uup, uu0, pp, ppp, fracdtG, realdt, ptrgeom, q, f1, f1norm, iJ);
        int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, impepsjac, iter, errorabsf1, dimfactU, uu, uup, uu0, pp, ppp, fracdtG, realdt, ptrgeom, q, f1, f1norm, iJ);
        //        get_implicit_iJ_old(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, impepsjac, iter, errorabsf1, dimfactU, uu, uup, uu0, pp, ppp, fracdtG, realdt, ptrgeom, q, f1, f1norm, iJ2);

        //        DLOOP(ii,jj) if(iJ[URAD0+ii][URAD0+jj]-iJ2[URAD0+ii][URAD0+jj]!=0.0) dualfprintf(fail_file,"ii= %d jj= %d iJ= %21.15g iJ2= %21.15g diff= %21.15g\n",ii,jj,iJ[URAD0+ii][URAD0+jj],iJ2[URAD0+ii][URAD0+jj],iJ[URAD0+ii][URAD0+jj]-iJ2[URAD0+ii][URAD0+jj]);


        if(failreturniJ!=0){
          if(havebackup){
            failreturn=FAILRETURNMODESWITCH; mathfailtype=30;
            if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"SWITCHING MODE: Detected bad Jacobian\n");
            break;
          }
          else{
            failreturn=FAILRETURNJACISSUE; mathfailtype=12;
            break;
          }
        }

        if(showmessagesheavy){
          int iii,jjj;
          JAC2DLOOP(iii,jjj,startjac,endjac)  dualfprintf(fail_file,"iJ[i %d][e %d]=%g\n",iii,jjj,iJ[irefU[iii]][erefU[jjj]]);
        }




        /////////
        //
        // Newton step (uup or ppp)
        //
        // DAMPFACTOR unused so far because don't know a priori whether to damp.  fracuup does post-inversion effective damping of this Newton step.
        // Newton step: x = x0 - (df/dx)^{-1}|_{x=x0} f(x0)
        //
        // Only updates 4D part of NPR data
        //
        /////////


        //////////////
        //
        // ITERATING U
        //
        ///////////////
        if(IMPUTYPE(implicititer)){
          PLOOP(pliter,pl) uu[pl] = uup[pl];
          JAC2DLOOP(ii,jj,startjac,endjac) uu[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];

          if(POSTNEWTONCONVCHECK==2){
            // check if any actual changes in primitives.  If none, then have to stop.
            FTYPE diffuu=0.0,sumuu=0.0;
            PLOOP(pliter,pl){
              diffuu += fabs(uu[pl]-uup[pl]);
              sumuu += fabs(uu[pl])+fabs(uup[pl]);
            }
            if(diffuu<DIFFXLIMIT*sumuu){
              convreturnf3limit=1;
            }
          }

          if(showmessagesheavy) dualfprintf(fail_file,"POSTDX: uu: %g %g %g %g : uup=%g %g %g %g\n",uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],uup[irefU[0]],uup[irefU[1]],uup[irefU[2]],uup[irefU[3]]);
        }// end iterating U


        //////////////
        //
        // ITERATING P
        //
        ///////////////
        else if(IMPPTYPE(implicititer)){

          PLOOP(pliter,pl) pp[pl]=ppp[pl];
          JAC2DLOOP(ii,jj,startjac,endjac){
            pp[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
            if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"added to ppp=%21.15g ii=%d jj=%d irefU=%d an amount of negative %21.15g\n",ppp[irefU[ii]],ii,jj,irefU[ii],DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]]);
          }



          
          // DEBUG: store steps in case hit max iter and want to debug
          if(DEBUGMAXITER&& dampattempt==0){
            PLOOP(pliter,pl) pppreholdlist[debugiter][pl]=pp[pl];
            if(DEBUGMAXITERVELOCITY==1){
              SLOOPA(jj){
                pppreholdlist[debugiter][U1+jj-1]=q->ucon[jj];
                pppreholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
              }
            }
            jac00list[debugiter]=iJ[irefU[0]][erefU[0]];
          }




          ///////////////////////////
          //
          //  POST NEWTON ADJUSTMENTS
          //
          ///////////////////////////


          // HOLD TT-iterated quantity or momentum-iterated quantity in initial steps
          if(RAMESHFIXEARLYSTEPS){
            // RAMESH  HOLD
            if(iter<RAMESHFIXEARLYSTEPS) pp[irefU[0]]=ppp[irefU[0]]; // don't trust first Newton step in u_g, Erf, or S
            else if(iter==RAMESHFIXEARLYSTEPS) SLOOPA(jj) pp[irefU[jj]]=ppp[irefU[jj]]; // don't trust second Newton step in velocity-momentum.

            if(pp[RHO]<=0.0||pp[UU]<=0.0){
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Detected negative rho=%21.15g ug=%21.15g\n",iter,pp[RHO],pp[UU]);
            }
          }// end Ramesh hold




          // check if u_g<0.  Do even if RAMESHFIXEARLYSTEPS going or even if checkconv==0
          int eomcond=(eomtypelocal==EOMGRMHD || eomtypelocal==EOMDEFAULT && EOMDEFAULT==EOMGRMHD);
          //          FTYPE umin=10.0*calc_PEQ_ufromTrho(TEMPMIN,fabs(pp[RHO]));
          if(pp[irefU[0]]<=0.0 && (IMPPTYPE(implicititer)&&implicititer!=QTYENTROPYPMHD)){ // don't consider implicititer==QTYENTROPYPMHD since S can be positive or negative.  Would only be unphysical or absolute-limited if the related u_g<0 or rho<0.
            if(JONHOLDPOS){
#if(1)
              holdingaspositive=0; // default
              if(countholdpositive<NUMHOLDTIMES){
                if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"HOLDING: Detected unphysical iter=%d countholdpositive=%d : pp[irefU[0]]=%g : ijknstepsteppart=%d %d %d %ld %d\n",iter,countholdpositive,pp[irefU[0]],ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
                //pp[irefU[0]]=MAX(100.0*NUMEPSILON*fabs(pp[RHO]),fabs(ppp[irefU[0]])); // hold as positive -- ppp might be too large so hold might be too aggressive to be useful.
                //                if(gotbest) umin=MAX(umin,0.5*fabs(pp[UU])); // override with last best version of u_g, because using very minimum leads to jump in behavior. // causes problems, leads to many high error events.
                //pp[irefU[0]]=umin; // hold as positive.
                pp[irefU[0]]=0.5*fabs(ppp[irefU[0]]); // just drop by half of positive value of *previous* value, not of negative value.
                countholdpositive++;
                holdingaspositive=1;
              }
#else
#define NUMITERHOLD (eomcond ? 2 : 4)
              if(holdingaspositive==0 || holdingaspositive==1 && iter<iterhold+NUMITERHOLD){
                if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"HOLDING: Detected unphysical iter=%d iterhold+NUMITERHOLD=%d holdingaspositive=%d : pp[irefU[0]]=%g : ijknstepsteppart=%d %d %d %ld %d\n",iter,iterhold+NUMITERHOLD,holdingaspositive,pp[irefU[0]],ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
                if(holdingaspositive==0) iterhold=iter;
                //              else holdingaspositive=1;
                holdingaspositive=1;
                //                pp[irefU[0]]=100.0*NUMEPSILON*fabs(pp[RHO]); // hold as positive just one iteration
                //                if(gotbest) umin=MAX(umin,0.5*fabs(pp[UU])); // override with last best version of u_g, because using very minimum leads to jump in behavior. // causes problems, leads to many high error events.
                //pp[irefU[0]]=umin; // hold as positive
                pp[irefU[0]]=0.5*fabs(ppp[irefU[0]]); // just drop by half of positive value of *previous* value, not of negative value.
              }
#endif

              else{// then exceeding hold attempts

                // if pre-normal step, skip to next type of step because non-normal step seems to want wrong/bad u_g anyways.
                if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
                  // revert to previous stepping's u_g, since using modified u_g (after u_g<0) should be bad.
                  PLOOP(pliter,pl){
                    pp[pl] = pppriorsteptype[pl];
                    uu[pl] = uupriorsteptype[pl];
                  }                    
                  iter=BEGINENERGYSTEPS-1;
                  continue;
                }
                else if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
                  // revert previous stepping's u_g, since using modified u_g (after u_g<0) should be bad.
                  PLOOP(pliter,pl){
                    pp[pl] = pppriorsteptype[pl];
                    uu[pl] = uupriorsteptype[pl];
                  }                    
                  iter=BEGINFULLSTEPS-1;
                  continue;
                }
                else{// else doing normal steps
                  pp[irefU[0]]=0.5*fabs(ppp[irefU[0]]); // just drop by half of positive value of *previous* value, not of negative value.
                  if(DOFINALCHECK){
                    if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"Unable to hold off u_g<0, setting canbreak=1 and letting finalchecks confirm error is good or bad.  iter=%d\n",iter);
                    canbreak=1; // just break since might be good (or at least allowable) error still.  Let final error check handle this.
                    // not fail.
                    mathfailtype=89;
                    break;
                  }
                  else if(havebackup){
                    failreturn=FAILRETURNMODESWITCH; mathfailtype=90;
                    if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"SWITCHING MODE: Detected bad u_g\n");
                    break;
                  }
                  else{
                    // then full failure
                    failreturn=FAILRETURNGENERAL; mathfailtype=10;
                    break;
                  }
                }// end else normal step
              } // if beyond hold counts allowed
            } // end if JONHOLDPOS
            else{//  not using JONHOLDPOS
              // if pre-normal step, skip to next type of step because non-normal step may have wrong u_g anyways.
              if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1;    continue;    }
              else if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1;   continue; }
              else{// else doing normal steps
                if(havebackup){
                  failreturn=FAILRETURNMODESWITCH; mathfailtype=90;
                  if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"SWITCHING MODE: Detected unphysical pp[irefU[0]]: iter=%d\n",iter);
                  break;
                }
              }
            }// else if not using JONHOLDPOS
          }// end if u_g<0 and iterating primitives
          else{
            // if u_g>0, not holding.
            holdingaspositive=0;
          }


          // DEBUG: store steps in case hit max iter and want to debug
          if(DEBUGMAXITER&& dampattempt==0){
            PLOOP(pliter,pl) ppposholdlist[debugiter][pl]=pp[pl];
            if(DEBUGMAXITERVELOCITY==1){
              SLOOPA(jj){
                ppposholdlist[debugiter][U1+jj-1]=q->ucon[jj];
                ppposholdlist[debugiter][URAD1+jj-1]=q->uradcon[jj];
              }
            }
          }


          // determine if holding so post-newton checks know.
          notholding=(RAMESHFIXEARLYSTEPS && iter>=RAMESHFIXEARLYSTEPS || RAMESHFIXEARLYSTEPS==0) && (JONHOLDPOS && holdingaspositive==0 || JONHOLDPOS==0);



          ///////////////////////////
          //
          //  POST NEWTON CHECKS
          //
          ///////////////////////////


          // only do post-newton checks if checking convergence allowed
          if(checkconv==1){

            // check if energy u_g too often bad compared to entropy u_g
            // assume this is only done after 
            if(RAMESHSTOPENERGYIFTOOOFTENBELOWENTROPY){
              if((implicititer==QTYPMHD || implicititer==QTYPMHDENERGYONLY)  && didentropyalready){
                if(notholding==1 || implicititer==QTYPMHD){ // if holding on energy equation, don't  use this u_g check.  But, if normal steps, then ignore holding and assume holding means u_g bad if on normal steps.
                  if(BADENERGY(pp[irefU[0]],pborig[irefU[0]])) countbadenergy++;
                  if(countbadenergy>=RAMESHSTOPENERGYIFTOOOFTENBELOWENTROPY){

                    // if pre-normal step, skip to next type of step because non-normal step may have wrong u_g anyways.
                    if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1;    continue;    }
                    else if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1;   continue; }
                    else{// else doing normal steps
                      // "switch" to entropy by just stopping trying to get energy solution
                      failreturn=FAILRETURNMODESWITCH; mathfailtype=100;
                      if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"SWITCHING MODE: Detected entropy u_g preferred consistently: iter=%d: %g %g\n",iter,pp[irefU[0]],pborig[irefU[0]]);
                      break;
                    }
                  }// end if normal step
                }//end if at point where iterations can be considered
              }// whether have information necessary to check
            }// whether to check u_g energy vs. entropy



            if(POSTNEWTONCONVCHECK==2 && notholding==1){
              // check if any actual changes in primitives.  If none, then have to stop.
              FTYPE diffpp=0.0,sumpp=0.0;
              PLOOP(pliter,pl){
                diffpp += fabs(pp[pl]-ppp[pl]);
                sumpp += fabs(pp[pl])+fabs(ppp[pl]);
              }
              if(diffpp<DIFFXLIMIT*sumpp){
                convreturnf3limit=1;
              }
            }
          }// end if checkconv==1




          if(showmessagesheavy) dualfprintf(fail_file,"POSTDX: pp: %g %g %g %g : ppp=%g %g %g %g\n",pp[irefU[0]],pp[irefU[1]],pp[irefU[2]],pp[irefU[3]],ppp[irefU[0]],ppp[irefU[1]],ppp[irefU[2]],ppp[irefU[3]]);



        }// end if iterating primitves


   



        ///////////////////////////
        //
        //  POST NEWTON ADJUSTMENTS part 2
        //
        ///////////////////////////


        // only do post-newton checks if checking convergence allowed
        if(checkconv==1){
          FTYPE f3[NPR]={0},f3norm[NPR]={0};
          if(POSTNEWTONCONVCHECK==1 && notholding==1){
            int dimtypef;
            /////////
            //
            // test convergence after Newton step
            // test convergence using |dU/U|
            // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
            // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
            /////////
            // 0 = conserved R^t_\nu type, 1 = primitive R^{ti} type
            if(IMPUTYPE(implicititer)){ // still considering iterate error
              dimtypef=DIMTYPEFCONS;
              JACLOOPALT(ii,startjac,endjac){
                f3[erefU[ii]]=(uu[irefU[ii]]-uup[irefU[ii]]);
                f3norm[erefU[ii]]=fabs(uu[irefU[ii]])+fabs(uup[irefU[ii]]);
              }
            }
            else if(IMPPTYPE(implicititer)){ // still considering iterate error
              dimtypef=DIMTYPEFPRIM;
              JACLOOPALT(ii,startjac,endjac){
                f3[erefU[ii]]=(pp[irefU[ii]]-ppp[irefU[ii]]);
                f3norm[erefU[ii]]=fabs(pp[irefU[ii]])+fabs(ppp[irefU[ii]]);
              }
            }
  
            // store error and solution in case eventually lead to max iterations and actually get worse error
            // f_error_check(uu0,uu) is ok to use since it just normalizes error
            int convreturnf3;
            convreturnf3=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs,realdt, dimtypef,eomtypelocal ,itermode,fracenergy,dissmeasure,dimfactU,pp,piin,f3,f3norm,f3report,Uiin,uu0,uu,ptrgeom,&errorabsf3);
            // while using f1 for true error, can't do better if f3 error is below near machine precision.
            convreturnf3limit=(errorabsf3<LOCALPREIMPCONVXABS);
          }

          if(POSTNEWTONCONVCHECK==0 || notholding==0){
            convreturnf3limit=0;
          }


          /////////////////
          //
          // check convergence
          // then can check convergence: using f1 and f3limit
          //
          /////////////////
          if(convreturnf1 || convreturnf3limit || canbreak){
            if(debugfail>=DEBUGLEVELIMPSOLVERMORE){
              if(convreturnf3limit && debugfail>=3){
                dualfprintf(fail_file,"f3limit good\n");
                if(POSTNEWTONCONVCHECK==1) JACLOOPALT(ii,startjac,endjac) dualfprintf(fail_file,"ii=%d erefU[ii]=%d f3=%21.15g f3norm=%21.15g f3report=%21.15g\n",ii,erefU[ii],f3[erefU[ii]],f3norm[erefU[ii]],f3report[erefU[ii]]);          
              }
              if(convreturnf1) dualfprintf(fail_file,"f1 good: ijknstepsteppart=%d %d %d %ld %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
              if(convreturnf3limit) dualfprintf(fail_file,"f3 good: ijknstepsteppart=%d %d %d %ld %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
              if(canbreak) dualfprintf(fail_file,"canbreak=%d good: ijknstepsteppart=%d %d %d %ld %d\n",canbreak,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
            }
            // so done.
            break;
          }
          else{
            // then not done
          }

        } // end if checkconv==1

  
      }// end if finite


      /////////
      // see if took too many Newton steps or not finite results
      /////////
      if(iter>trueimpmaxiter || notfinite ){
        if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"iter=%d>%d or notfinite=%d\n",iter,trueimpmaxiter,notfinite);
        //      failreturn=FAILRETURNGENERAL; // no don't fail, might be allowable error.
        break;
      }

    }// end do
    while(1);













  


    //////////////
    //
    // if no failure, then see if really still failed or not with some final checks.
    //
    ///////////////
    int failreturnf=0;
    if(failreturn==0 && earlylowerror==0){


      int fakeiter;
      for(fakeiter=0;fakeiter<=0;fakeiter++){

        int convreturn=1,convreturnallow=1; // default is solution is acceptable.


        if(DOFINALCHECK){
          //////////////////////////
          //
          // check and get error for last iteration or any mods from above that are post-iteration
          //
          // The call to f_implicit() also ensures uu is consistent with new pp
          //
          ////////////////////////
          int whichcall=2;
          //  eomtypelocal=*eomtype; // re-chose default each time. No, stick with what f1 (last call to f1) chose
          int goexplicitfake; // not used here
          failreturnf=f_implicit(iter,failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, fracenergy, dissmeasure, radinvmod, pp, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, &goexplicitfake); // modifies uu and pp and q
          // radinvmod contains whether radiative inversion modified process.

          int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive (u,v^i) type, i.e. v^i has no energy density term
          convreturn=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv,trueimptryconvabs,realdt,dimtypef,eomtypelocal,itermode,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1);
          convreturnallow=(errorabsf1<IMPALLOWCONVABS);
          if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"DOFINALCHECK: convreturn=%d convreturnallow=%d (IMPALLOWCONV=%g) f1report: %g %g %g %g : %g\n",convreturn,convreturnallow,IMPALLOWCONV,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],errorabsf1);
        }// end if doing final check
        else{
          // kinda risky to rely upon last step but not checking its error

          if(IMPPTYPE(implicititer)){
            // since iterated pp or uu but didn't call f_implicit(), uu or pp (respectively) is no longer consistent.
            // So get uu(pp)
            struct of_state qcons;
            get_state(pp, ptrgeom, &qcons);
            primtoU(UNOTHING,pp,&qcons,ptrgeom, uu, NULL);
          }

          // now get error
          int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive (u,v^i) type, i.e. v^i has no energy density term
          convreturn=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs,realdt,dimtypef,eomtypelocal,itermode,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1);
          convreturnallow=(errorabsf1<IMPALLOWCONVABS);
        }



        if(GETBEST){
      
          if(gotbest){
            // f1-based
            // using old uu,uup, but probably ok since just helps normalize error
            errorabsf1=0.0;     JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsf1     += fabs(f1report[erefU[jj]]);

            // f1 based
            // see if should revert to prior best
            errorabsbest=0.0; JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsbest += fabs(lowestfreportf1[erefU[jj]]);
            if(errorabsbest<errorabsf1 || !isfinite(errorabsf1) ){
              PLOOP(pliter,pl) uu[pl]=bestuu[pl];
              PLOOP(pliter,pl) pp[pl]=bestpp[pl];

              JACLOOPALT(jj,startjac,endjac) f1report[erefU[jj]] = lowestfreportf1[erefU[jj]];
              convreturn=1; JACLOOPFULLERROR(itermode,jj,startjac,endjac) convreturn*=fabs(f1report[erefU[jj]])<trueimptryconv; // like doing &&
              convreturnallow=1; JACLOOPFULLERROR(itermode,jj,startjac,endjac) convreturnallow*=fabs(f1report[erefU[jj]])<IMPALLOWCONV; // like doing &&
              if(showmessages && debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Using best: %g %g\n",errorabsf1,errorabsbest);
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"GETBEST: convreturn=%d convreturnallow=%d (IMPALLOWCONV=%g) f1report: %g %g %g %g : %g\n",convreturn,convreturnallow,IMPALLOWCONV,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],errorabsf1);
            }
            else{
              PLOOP(pliter,pl) bestuu[pl]=uu[pl];
              PLOOP(pliter,pl) bestpp[pl]=pp[pl];
              errorabsbest=errorabsf1;
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"gotbest=%d but errorabsbest=%g while errorabsf1=%g\n",gotbest,errorabsbest,errorabsf1);
            }
          }
        }// end GETBEST












        // KORALTODO: If convreturnallow doesn't work, but still (say) 10% error, might want to hold onto result in case explicit backup fails as well (which is likely), in which case *much* better to use 10% error because otherwise 4-force not accounted for, which can lead to very big changes in fluid behavior due to large flux from previous step.
        // KORALTODO: Or, perhaps should really take it as a failure and use fixups.  Probably should allow for result to be written if error<10%, but only use as super-backup in fixups.  So should set pflag still.

        ///////////////
        //
        // See if reached desired tolerance
        //
        ///////////////
        if(convreturn){
          // then really done
          failreturn=FAILRETURNNOFAIL; mathfailtype=0;
          break;
        }
        else{ // convreturn==0
          
          // see if at least allowed error
          if(convreturnallow){

            // set as soft allowable failure
            failreturn=FAILRETURNNOTTOLERROR; mathfailtype=202;

            if(iter>trueimpmaxiter){// then reached maximum iterations
              dualfprintf(fail_file,"trueimpmaxiter=%d eomtype=%d MAXcheckconv=%d havebackup=%d failreturnallowable=%d: f1report=%g %g %g %g : f1=%g %g %g %g\n",trueimpmaxiter,*eomtype,checkconv,havebackup,failreturnallowable,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]]);

              if(showmessages && debugfail>=2) dualfprintf(fail_file,"iter>trueimpmaxiter=%d : iter exceeded in solve_implicit_lab().  But f1 was allowed error. checkconv=%d (if checkconv=0, could be issue!) : %g %g %g %g : %g %g %g %g : errorabs=%g : %g %g %g\n",trueimpmaxiter,checkconv,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]],errorabsf1,fracdtuu0,fracuup,fracdtG);
              if(REPORTMAXITERALLOWED){
                if(havebackup){
                  if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"SWITCHING MODE: Detected MAXITER\n");
                  // don't break, just reporting or not
                  mathfailtype=50;
                }
                else{
                  mathfailtype=(eomtypelocal==EOMGRMHD ? 6 : 600);
                }
              }
            }

            // then nothing else to do
            break;
          }
          else{ // didn't reach allowable error
            // not allowable failure
            failreturn=FAILRETURNGENERAL; mathfailtype=203;

            // KORALTODO: Need backup that won't fail.
            if(debugfail>=2){
              if(canbreak==1 && havebackup==0) dualfprintf(fail_file,"Held u_g, couldn't hold anymore and broke, but error still larger than allowed : iter=%d ijknstepsteppart=%d %d %d %ld %d\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
              if(canbreak==2 && havebackup==0) dualfprintf(fail_file,"Aborted due to oscillatory error despite not having backup.\n");
              if(canbreak==3 && havebackup==0) dualfprintf(fail_file,"Aborted due to error not decreasing fast enough: errorabsf1=%g\n",errorabsf1);
              if(iter>trueimpmaxiter && havebackup==0) dualfprintf(fail_file,"iter>trueimpmaxiter=%d : iter exceeded in solve_implicit_lab(). nstep=%ld steppart=%d ijk=%d %d %d :  Bad error.\n",trueimpmaxiter,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
              if(notfinite && havebackup==0) dualfprintf(fail_file,"IMPGOTNAN at iter=%d : in solve_implicit_lab(). ijk=%d %d %d :  Bad error.\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);
              if(havebackup==0) dualfprintf(fail_file,"checkconv=%d havebackup=%d failreturnallowable=%d: f1report=%g %g %g %g : f1=%g %g %g %g\n",checkconv,havebackup,failreturnallowable,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]]);
              if(1||showmessages){
                if(havebackup){
                  // don't break, just don't report.

                  mathfailtype=60;
                }
                else{
                  mathfailtype=1;
                }
              }

            }// debug

            // nothing else to do except leave
            break;
          }// end if convreturnallow=0 (didn't reach allowed error)
        }// end if didn't reach desired error



      

      }// end fake loop that can break out of
    }// end if failreturn==0 originally


    // get best failreturn
    if(GETBEST){
      if(gotbest){
        failreturnbest=failreturn;
      }
    }


    // estimate effective work based on iterations done for each equation type
    FTYPE fndim=(FTYPE)NDIM;
    totaliters += (int)( (3.0/fndim)*(FTYPE)momiters + (1.0/fndim)*(FTYPE)energyiters + (fndim/fndim)*(FTYPE)fulliters);

  }// end loop over damping
  if(dampattempt==NUMDAMPATTEMPTS){
    if(debugfail>=2) dualfprintf(fail_file,"Damping failed to avoid max iterations (but error might have dropped): failreturn=%d dampattempt=%d eomtypelocal=%d *eomtype=%d\n",failreturn,dampattempt,eomtypelocal,*eomtype);
  }

  ///////////
  //
  // Once damping done, ensure choose best over all damp tries
  //
  ///////////
  if(GETBEST){
    if(gotbest){
      PLOOP(pliter,pl) uu[pl]=bestuu[pl];
      PLOOP(pliter,pl) pp[pl]=bestpp[pl];
      errorabsf1=errorabsbest;
      failreturn=failreturnbest;
    }
  }

  ///////////
  //
  // have to compute final uu[ENTROPY] if doing entropy optimization where avoid it during iterations if not needed.
  //
  ///////////
  if(ENTROPYOPT){
    int needentropy=1;
    get_state_norad_part2(needentropy, pp, ptrgeom, q); // where entropy would be computed
    //    get_state(pp, ptrgeom, q);
    extern int primtoflux_nonradonly(int needentropy, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
    FTYPE uuentropy[NPR];
    primtoflux_nonradonly(needentropy,pp,q,TT,ptrgeom, uuentropy, NULL);
    uu[ENTROPY]=uuentropy[ENTROPY];
  }



  /////////////////////
  //
  // if didn't fail to get some reasonable solution, then now can use it.
  //
  /////////////////////
  if(ACCEPTASNOFAILURE(failreturn)==1){

    ///////////////////
    //
    // get source update as "dU" = dU/dt using real dt that used during implicit iterations, and will eventually use to update U in advance.c.
    // apply source update as force
    // KORALNOTE: As long as f_implicit() updates both full primitive and full U, using uu below is fine and good.
    PLOOP(pliter,pl) radsource[pl] = +(uu[pl]-uu0[pl])/realdt;
    // KORALNOTE: Could re-enforce energy conservation here, but would be inconsistenet with how applied error function.
    PLOOPBONLY(pl) radsource[pl]  = 0.0; // force to machine accuracy

    DLOOPA(ii) radsource[UU+ii] = -radsource[URAD0+ii]; // force energy conservation


    // OLD, but misses rho changes due to u^t changes:
    //  DLOOPA(jj) radsource[iotherU[jj]]  = -(uu[irefU[jj]]-uu0[irefU[jj]])/realdt;
    //  DLOOPA(jj) radsource[irefU[jj]]    = +(uu[irefU[jj]]-uu0[irefU[jj]])/realdt;


    // DEBUG:
    //  DLOOPA(jj) dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d implicitGd[%d]=%g %g\n",nstep,steppart,ptrgeom->i,jj,radsource[iotherU[jj]],radsource[irefU[jj]]);


    /////////////
    //
    // store source update in dUcomp for return.
    //
    //////////////
    sc = RADSOURCE;
    PLOOP(pliter,pl) dUcomp[sc][pl] += radsource[pl];

    ////////////////
    //
    // save better guess for later inversion from this inversion
    // pp was modified by f_implicit(f1,pp) with output from inversion returned through pp0
    // only use pp if successful with implicit method, since if not successful can be various bad reasons with no good pb
    //
    ////////////////
    PLOOP(pliter,pl){
      pb[pl]=pp[pl]; // actual solution that's used
      uub[pl]=uu[pl]; // used for getting whether solution really worked and switching methods, etc.
    }

    // DEBUG:
    //  PLOOP(pliter,pl) dualfprintf(fail_file,"POOP2: pl=%d uu=%21.15g uu0=%21.15g piin=%21.15g pb=%21.15g\n",pl,uu[pl],uu0[pl],piin[pl],pb[pl]);

    ///////////
    //
    // choose new eomtype for external inverison or other checks
    //
    // This ensures whatever implicit solver settled on using, the external inversion doesn't switch.
    //
    //////////
    *eomtype=eomtypelocal;



    if(SWITCHTODONOTHING){
      // always do nothing, so entropy won't try to revert to cold.
      if(EOMENTROPYGRMHD) *eomtype=EOMDIDENTROPYGRMHD;
      if(EOMCOLDGRMHD) *eomtype=EOMDIDCOLDGRMHD;
      if(EOMGRMHD) *eomtype=EOMDIDGRMHD;
    }

  }// end if didn't fail, so can set final solution.

 


  //////////////
  //
  // report any bad failure (using previously set mathfailtype value)
  //
  //////////////

  // for checking cases where tau>=1 but still Erf<0
  //  FTYPE tautot[NDIM],tautotmax;
  //  calc_tautot(pp, ptrgeom, tautot, &tautotmax);
  //  //  if(tautotmax>1 && pp[PRAD0]<10.0*ERADLIMIT){
  //  if(tautotmax>2 && pp[PRAD0]<10.0*ERADLIMIT){


  if(PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1>=IMPTRYCONVALT || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){ // as in previous code

  // for seeing Erf<0 and small errors not tol errors.
  //  if(failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1>=IMPTRYCONVALT || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){

  // for catching oscillators at small error but still >tol.
  //  if(PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0){


  //    if(REPORTERFNEG && failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT || *radinvmod ) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1>=trueimptryconvalt || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){
  //  if(REPORTERFNEG && failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1>=trueimptryconvalt || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){
    //    if(NOTBADFAILURE(failreturn)==0){
    struct of_state qcheck; get_state(pp, ptrgeom, &qcheck);  primtoU(UNOTHING,pp,&qcheck,ptrgeom, uu, NULL);
    failnum++; mathematica_report_check(*radinvmod, mathfailtype, failnum, gotfirstnofail, eomtypelocal, itermode, errorabsf1, errorabsbestexternal, iter, totaliters, realdt, ptrgeom, ppfirst,pp,pb,piin,prtestUiin,prtestUU0,uu0,uu,Uiin,Ufin, CUf, q, dUother);
    int usedebugiter=debugiteratteempts[0];
    showdebuglist(usedebugiter,pppreholdlist,ppposholdlist,f1reportlist,f1list,errorabsf1list,realiterlist,jac00list);
  }




  ////////////////
  //
  // Return error and iterations and fail mode
  //
  ////////////////

  // report error no matter whether got solution or not.
  *errorabsreturn=errorabsf1;

  // report iters not matter what the error.
  *itersreturn=totaliters;

  return(failreturn);
  
}












// 0 : full
// 1: optimal
#define DEBUGMAXMODE 1

// DEBUGMAXITER stuff
static void showdebuglist(int debugiter, FTYPE (*pppreholdlist)[NPR],FTYPE (*ppposholdlist)[NPR],FTYPE (*f1reportlist)[NDIM],FTYPE (*f1list)[NDIM],FTYPE *errorabsf1list, int *realiterlist, FTYPE *jac00list)
{

  if(DEBUGMAXITER==0) return;

  int listiter;
  if(DEBUGMAXMODE==0) dualfprintf(fail_file,"%3s : %3s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s %21s %21s : %21s : %21s\n","li","ri","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","f1rep0","f1rep1","f1rep2","f1rep3","errorabs","umin");
  else if(DEBUGMAXMODE==1) dualfprintf(fail_file,"%3s : %3s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s : %21s %21s %21s %21s : %21s %21s %21s %21s : %21s : %21s : %21s\n","li","ri","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","rho","ug","f10","f11","f12","f13","f1rep0","f1rep1","f1rep2","f1rep3","errorabs","umin","jac00");

  for(listiter=0;listiter<=debugiter;listiter++){
    FTYPE umin=calc_PEQ_ufromTrho(TEMPMIN,pppreholdlist[listiter][RHO]);
    if(DEBUGMAXMODE==0){
      // full, but excessive
      dualfprintf(fail_file
                  ,"%3d : %3d : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g : %21.15g : %21.15g\n"
                  ,listiter,realiterlist[listiter]
                  ,pppreholdlist[listiter][RHO],pppreholdlist[listiter][UU],pppreholdlist[listiter][U1],pppreholdlist[listiter][U2],pppreholdlist[listiter][U3],pppreholdlist[listiter][PRAD0],pppreholdlist[listiter][PRAD1],pppreholdlist[listiter][PRAD2],pppreholdlist[listiter][PRAD3]
                  ,ppposholdlist[listiter][RHO],ppposholdlist[listiter][UU],ppposholdlist[listiter][U1],ppposholdlist[listiter][U2],ppposholdlist[listiter][U3],ppposholdlist[listiter][PRAD0],ppposholdlist[listiter][PRAD1],ppposholdlist[listiter][PRAD2],ppposholdlist[listiter][PRAD3]
                  ,f1reportlist[listiter][0],f1reportlist[listiter][1],f1reportlist[listiter][2],f1reportlist[listiter][3]
                  ,errorabsf1list[listiter]
                  ,umin
                  );
    }
    else if(DEBUGMAXMODE==1){
      // optimal
      dualfprintf(fail_file
                  ,"%3d : %3d : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g : %21.15g : %21.15g : %21.15g\n"
                  ,listiter,realiterlist[listiter]
                  ,pppreholdlist[listiter][RHO],pppreholdlist[listiter][UU],pppreholdlist[listiter][U1],pppreholdlist[listiter][U2],pppreholdlist[listiter][U3],pppreholdlist[listiter][PRAD0],pppreholdlist[listiter][PRAD1],pppreholdlist[listiter][PRAD2],pppreholdlist[listiter][PRAD3]
                  ,ppposholdlist[listiter][RHO],ppposholdlist[listiter][UU]
                  ,f1list[listiter][0],f1list[listiter][1],f1list[listiter][2],f1list[listiter][3]
                  ,f1reportlist[listiter][0],f1reportlist[listiter][1],f1reportlist[listiter][2],f1reportlist[listiter][3]
                  ,errorabsf1list[listiter]
                  ,umin,jac00list[listiter]
                );
    }
  }// end listiter loop

}







                             
int get_rameshsolution_wrapper(int whichcall, int eomtype, FTYPE errorabs, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *dUother, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE (*dUcompeng)[NPR], FTYPE (*dUcompent)[NPR], struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent)
{
  // BEGIN get ramesh solution
  int radinvmod=0; // fake
  int failreturnentropy=0; // fake
  //  int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITER,itersent=IMPMAXITER;
  int failnum=0,gotfirstnofail=0;
  int itermode=0,iters=0;
  FTYPE errorabsbestexternal=BIG;
  // itersentropy is last, but might feed in best
  FTYPE realdt = compute_dt(CUf,dt);
  // get uu0
  FTYPE uu0[NPR];
  FTYPE fracdtuu0=1.0;
  int pliter,pl;
  PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
  FTYPE uu[NPR]; // not used except when doing diagnostics in ramesh code
  FTYPE uueng[NPR],uuent[NPR]; // filled with answer if successful
  //
  get_rameshsolution(whichcall, radinvmod, failreturnentropy, failnum,  gotfirstnofail,  eomtype,  itermode, errorabs, errorabs, iters, iters, realdt, ptrgeom, pp, pp, piin, uu0, uu, Uiin, Ufin, CUf, q, ppeng, ppent, uueng, uuent, qeng, qent, failtypeeng, errorabseng, iterseng, radinvmodeng, failtypeent, errorabsent, itersent, radinvmodent);
  //
  // pp and q are assigned, but external call might just use only *eng and *ent versions of these and other things.
  if(eomtype==EOMGRMHD){
    PLOOP(pliter,pl) pp[pl]=ppeng[pl];
    *q=*qeng;
  }
  if(eomtype==EOMENTROPYGRMHD){
    PLOOP(pliter,pl) pp[pl]=ppent[pl];
    *q=*qent;
  }
  //
  int sc;
  sc = RADSOURCE;
  PLOOP(pliter,pl){
    dUcompeng[sc][pl] = + (uueng[pl]-uu0[pl])/realdt;
    dUcompent[sc][pl] = + (uuent[pl]-uu0[pl])/realdt;
  }
  // END WITH RAMESH SOLUTION

#if(0)
  // DEBUG
  sc = RADSOURCE;
  PLOOP(pliter,pl){
    if(1||!isfinite(dUcompeng[sc][pl])) dualfprintf(fail_file,"POOPENG %g %g %g\n",uueng[pl],uu0[pl],realdt);
    if(1||!isfinite(dUcompent[sc][pl])) dualfprintf(fail_file,"POOPENT %g %g %g\n",uuent[pl],uu0[pl],realdt);
  }
#endif


  return(0);
}

//#define WHICHVELRAMESH VEL4 // should be same as WHICHVEL in test.f
#define WHICHVELRAMESH VELREL4 // should be same as WHICHVEL in test.f

int get_rameshsolution(int whichcall, int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, FTYPE errorabs, FTYPE errorabsbestexternal, int iters, int totaliters, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent)
{

  return(0);
}



int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, FTYPE errorabs, FTYPE errorabsbestexternal, int iters, int totaliters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother)
{
  int jj,kk;
  int pliter,pl;


  if(0){ // old mathematica style
    dualfprintf(fail_file,"FAILINFO: %d %d %lld %d\ndt=%21.15g\n",failtype, myid, failnum, gotfirstnofail,realdt);
    DLOOP(jj,kk) dualfprintf(fail_file,"gn%d%d=%21.15g\n",jj+1,kk+1,ptrgeom->gcon[GIND(jj,kk)]);
    DLOOP(jj,kk) dualfprintf(fail_file,"gv%d%d=%21.15g\n",jj+1,kk+1,ptrgeom->gcov[GIND(jj,kk)]);
    // shows first pp(uu0)
    PLOOP(pliter,pl) dualfprintf(fail_file,"pp%d=%21.15g\npb%d=%21.15g\nuu0%d=%21.15g\nuu%d=%21.15g\nuui%d=%21.15g\n",pl,pp[pl],pl,pb[pl],pl,uu0[pl],pl,uu[pl],pl,Uiin[pl]);
    struct of_state qreport;
    get_state(pp,ptrgeom,&qreport);
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"uradcon%d=%21.15g\nuradcov%d=%21.15g\n",jj,qreport.uradcon[jj],jj,qreport.uradcov[jj]);
    else DLOOPA(jj) dualfprintf(fail_file,"uradcon%d=%21.15g\nuradcov%d=%21.15g\n",jj,0.0,jj,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"ucon%d=%21.15g\nucov%d=%21.15g\n",jj,qreport.ucon[jj],jj,qreport.ucov[jj]);
    // then do:
    // 1) grep -A 134 --text FAILINFO 0_fail.out.grmhd* > fails.txt
    // 2) emacs regexp:  \([0-9]\)e\([-+]*[0-9]+\)   ->   \1*10^(\2)
  }
  else{

    /////////////////////////
    //
    // Fix-up some terms to avoid nan or inf in output.
    //
    /////////////////////////

    if(!isfinite(pp[UU])) pp[UU]=BIG;
    if(!isfinite(pp[ENTROPY])) pp[ENTROPY]=BIG;
    if(!isfinite(errorabs)) errorabs=BIG;
    if(!isfinite(errorabsbestexternal)) errorabsbestexternal=BIG;

    /////////////////////////
    //
    // Get state for pp,pb,piin
    //
    /////////////////////////


    struct of_state qpp;
    get_state(pp,ptrgeom,&qpp);
    struct of_state qpb;
    get_state(pb,ptrgeom,&qpb);
    struct of_state qpiin;
    get_state(piin,ptrgeom,&qpiin);

    /////////////////////////
    //
    // Output stuff for mathematica
    //
    /////////////////////////

    // 211+11 things
    dualfprintf(fail_file,"\nFAILINFO: ");
    dualfprintf(fail_file,"%d %d %d ",1,1,1); // 3
    dualfprintf(fail_file,"%d %d %lld %d %d %d %21.15g %21.15g %d %d %21.15g %lld %d %21.15g ",failtype,myid,failnum,gotfirstnofail,eomtypelocal,itermode,errorabs,errorabsbestexternal,iters,totaliters,realdt,nstep,steppart,gam); // 14
    dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g ",GAMMAMAXRAD,ERADLIMIT,IMPTRYCONVABS,IMPALLOWCONVABS); // 4
    //    dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g ",ARAD_CODE,KAPPA_ES_CODE(1.0,1.0),KAPPA_FF_CODE(1.0,1.0),KAPPA_BF_CODE(1.0,1.0)); // 4
    dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g ",ARAD_CODE,calc_kappaes_user(1.0,1.0,0,0,0),calc_kappa_user(1.0,1.0,0,0,0),0.0); // 4
    DLOOP(jj,kk) dualfprintf(fail_file,"%21.15g ",ptrgeom->gcon[GIND(jj,kk)]); // 16
    DLOOP(jj,kk) dualfprintf(fail_file,"%21.15g ",ptrgeom->gcov[GIND(jj,kk)]); // 16
    PLOOP(pliter,pl) dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ",pp[pl],ppfirst[pl],pb[pl],piin[pl],prtestUiin[pl],prtestUU0[pl],uu0[pl],uu[pl],Uiin[pl]);  // 9*13
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpp.uradcon[jj],qpp.uradcov[jj]); // 4*2=8
    else DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpp.ucon[jj],qpp.ucov[jj]); // 4*2=8
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpb.uradcon[jj],qpb.uradcov[jj]); // 4*2=8
    else dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpb.ucon[jj],qpb.ucov[jj]); // 4*2=8
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpiin.uradcon[jj],qpiin.uradcov[jj]); // 4*2=8
    else dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpiin.ucon[jj],qpiin.ucov[jj]); // 4*2=8
    dualfprintf(fail_file,"\n");



    // get ramesh solution
    int whichcall=0; // debug call
    FTYPE ppeng[NPR]={0},ppent[NPR]={0},errorabseng=BIG,errorabsent=BIG;
    FTYPE uueng[NPR]={0},uuent[NPR]={0};
    int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITER,itersent=IMPMAXITER, radinvmodeng=-1, radinvmodent=-1;
    struct of_state qeng=*q,qent=*q;
    get_rameshsolution(whichcall, radinvmod, failtype, failnum,  gotfirstnofail,  eomtypelocal,  itermode, errorabs, errorabsbestexternal,  iters,  totaliters,realdt, ptrgeom, pp, pb, piin, uu0, uu, Uiin, Ufin, CUf, q, ppeng, ppent, uueng, uuent, &qeng, &qent, &failtypeeng, &errorabseng, &iterseng, &radinvmodeng, &failtypeent, &errorabsent, &itersent, &radinvmodent);



    /////////////////////////
    //
    // FAILRETURNABLE
    //
    /////////////////////////
    dualfprintf(fail_file,"\nFAILREPEATABLE: %d %d %lld : ",failtype,myid,failnum);
    dualfprintf(fail_file,"dt=%21.15g;CUf[2]=%21.15g;gam=%21.15g;",realdt,CUf[2],gam);
    PLOOP(pliter,pl) dualfprintf(fail_file,"pp[%d]=%21.15g;ppfirst[%d]=%21.15g;pb[%d]=%21.15g;piin[%d]=%21.15g;Uiin[%d]=%21.15g;Ufin[%d]=%21.15g;dUother[%d]=%21.15g;",pl,pp[pl],pl,ppfirst[pl],pl,pb[pl],pl,piin[pl],pl,Uiin[pl],pl,Ufin[pl],pl,dUother[pl]);
     // ptrgeom stuff
    DLOOP(jj,kk) dualfprintf(fail_file,"ptrgeom->gcov[GIND(%d,%d)]=%21.15g;ptrgeom->gcon[GIND(%d,%d)]=%21.15g;",jj,kk,ptrgeom->gcov[GIND(jj,kk)],jj,kk,ptrgeom->gcon[GIND(jj,kk)]);
    DLOOPA(jj) dualfprintf(fail_file,"ptrgeom->gcovpert[%d]=%21.15g;ptrgeom->beta[%d]=%21.15g;",jj,ptrgeom->gcovpert[jj],jj,ptrgeom->beta[jj]);
    dualfprintf(fail_file,"ptrgeom->alphalapse=%21.15g;ptrgeom->betasqoalphasq=%21.15g;ptrgeom->gdet=%21.15g;ptrgeom->igdetnosing=%21.15g;ptrgeom->i=%d;ptrgeom->j=%d;ptrgeom->k=%d;ptrgeom->p=%d;",ptrgeom->alphalapse,ptrgeom->betasqoalphasq,ptrgeom->gdet,ptrgeom->igdetnosing,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);
    if(q!=NULL){
      DLOOPA(jj) dualfprintf(fail_file,"q->ucon[%d]=%21.15g;q->ucov[%d]=%21.15g;",jj,q->ucon[jj],jj,q->ucov[jj]);
      if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"q->uradcon[%d]=%21.15g;q->uradcov[%d]=%21.15g;",jj,q->uradcon[jj],jj,q->uradcov[jj]);
      else DLOOPA(jj) dualfprintf(fail_file,"q->uradcon[%d]=%21.15g;q->uradcov[%d]=%21.15g;",jj,0.0,jj,0.0);
      dualfprintf(fail_file,"q->pressure=%21.15g;q->entropy=%21.15g;q->ifremoverestplus1ud0elseud0=%21.15g;",q->pressure,q->entropy,q->ifremoverestplus1ud0elseud0);
    }
    dualfprintf(fail_file,"\n");
    

    // then do:
    // 1) grep -h --text FAILINFO 0_fail.out.grmhd* | sed 's/FAILINFO: //g'| sort -T ./ -r -g -k 10 > fails.txt
    //
    // or:
    // grep -h --text FAILINFO 0_fail.out.grmhd* | grep -v "FAILINFO: 100" |grep -v "FAILINFO: 80"| sed 's/FAILINFO: //g' | sort -T ./ -r -g -k 7 > failshigherror.txt ; head -100 failshigherror.txt > failshigherror100.txt
    //

    // grep --text BAD 0_fail.out.grmhd.00*|wc -l ;  grep FAILINFO 0_fail.out.grmhd.00*|wc -l ; grep MAXF1 0_fail.out.grmhd.00*| wc -l ; grep MAXITER 0_fail.out.grmhd.00*|wc -l ; grep "also failed" 0_fail.out.grmhd.00*|wc -l 

    // TO CHECK ON DAMPING:
    //     grep -i "Damping worked" 0_fail.out.grmhd.00*|wc -l ; grep -i "Damping failed" 0_fail.out.grmhd.00*|wc -l

    //   see if any MAXF1ITER: less -S fails.txt| awk '{print $1}'|sort|less

    // 2) Choose numfails in below mathematica script

    // 3) ~/bin/math < /data/jon/harm_math/solveimplicit_superwrapper.m
    // or:  nohup ~/bin/math < /data/jon/harm_math/solveimplicit_superwrapper.m &> math.out &
    // or: preferred:  ~/bin/math < /data/jon/harm_math/solveimplicit_superwrapper.m &> math.out &

    // 4) grep 0Good math.out
    // if any Good's appear, mathematica found solution when harm did not.  Fix it!

    // If any appear in  grep 0WGood math.out   , then precision issue with long doubles in harm even.

    // If any appear in  grep 0MGood math.out   , then gammamax case should work!

    // 5) More advanced version of #4.  Use check.sh script:
    // bash /data/jon/harmgit/scripts/check.sh math.out .  Gives overall report.
  }



  return(0);
}


static int f_error_check_wrap(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, FTYPE *f1, FTYPE *f1norm, FTYPE *f1report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom)
{

  //  int convreturn=f_error_check(showmessages, showmessagesheavy, iter, IMPTRYCONV,realdt,f3,f3norm,f3report,Uiin,uup,uu,ptrgeom);
  FTYPE fpl1[NDIM],fpl1norm[NDIM],fpl1report[NDIM];

  int ii;
  DLOOPA(ii){
    fpl1[URAD0+ii]=f1[ii];
    fpl1norm[URAD0+ii]=f1norm[ii];
    fpl1report[URAD0+ii]=f1report[ii];
  }
  FTYPE errorabsf1=BIG;
  int convreturn=f_error_check(showmessages, showmessagesheavy, iter, IMPTRYCONV, IMPTRYCONVABS,realdt,0,EOMGRMHD,ITERMODENORMAL,1.0,-1.0,NULL,NULL,NULL,fpl1,fpl1norm,fpl1report,Uiin,uu0,uu,ptrgeom,&errorabsf1);

  DLOOPA(ii){
    f1report[ii]=fpl1report[URAD0+ii];
  }

  return(convreturn);

}

// use f and check the error
//static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, FTYPE *f1, FTYPE *f1norm, FTYPE *f3report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom)
static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE convabs, FTYPE realdt, int dimtypef, int eomtype, int itermode, FTYPE fracenergy, FTYPE dissmeasure, FTYPE *dimfactU, FTYPE *pp, FTYPE *piin, FTYPE *fpl1, FTYPE *fpl1norm, FTYPE *fpl3report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom, FTYPE *errorabs)
{
  int ii,jj;
  FTYPE f3[NDIM];
  FTYPE f3a[NDIM],f3b[NDIM],f3c[NDIM],f3d[NDIM];
  int passedconv;

  FTYPE f1[NDIM],f1norm[NDIM],f3report[NDIM];

  DLOOPA(ii){
    f1[ii]=fpl1[URAD0+ii];
    f1norm[ii]=fpl1norm[URAD0+ii];
    f3report[ii]=fpl3report[URAD0+ii];
  }

  // default
  passedconv=0;

  // get error
  // NOTE: use of gcov[ii,ii] so comparable dimensionally to f1[ii] and f1norm[ii] that are like R^t_\nu and so need sqrt(gcon[nu,nu]) multiplied on them.  This ensures error is non-dimensional (or, really only ^t dimensional)
  FTYPE dimfact[NDIM];
  DLOOPA(ii){
    dimfact[ii]=sqrt(fabs(ptrgeom->gcon[GIND(ii,ii)]));
  }    
  DLOOPA(ii){
    f3a[ii]=fabs(f1[ii]*dimfact[ii]/(SMALL+fabs(uu0[URAD0]*dimfact[0])));
    f3b[ii]=fabs(f1[ii]*dimfact[ii]/(SMALL+MAX(fabs(uu0[UU]*dimfact[0]),fabs(uu0[URAD0]*dimfact[0]))));
    f3c[ii]=fabs(f1[ii]*dimfact[ii]/(SMALL+MAX(fabs(f1norm[ii]*dimfact[ii]),fabs(uu0[URAD0]*dimfact[0]))));
    f3d[ii]=fabs(f1[ii]*dimfact[ii]/(SMALL+MAX(fabs((uu0[UU]-Uiin[UU])*dimfact[0]),MAX(fabs(f1norm[ii]*dimfact[ii]),fabs(uu0[URAD0]*dimfact[0])))));
  }

  // evaluate error
  if(IMPLICITERRORNORM==1){
    if(f3a[0]<conv && f3a[1]<conv && f3a[2]<conv && f3a[3]<conv) passedconv=1;
    DLOOPA(ii) f3report[ii]=f3a[ii];
  }
  else if(IMPLICITERRORNORM==2){
    if(f3b[0]<conv && f3b[1]<conv && f3b[2]<conv && f3b[3]<conv) passedconv=1;
    DLOOPA(ii) f3report[ii]=f3b[ii];
  }
  else if(IMPLICITERRORNORM==3){
    if(f3c[0]<conv && f3c[1]<conv && f3c[2]<conv && f3c[3]<conv) passedconv=1;
    DLOOPA(ii) f3report[ii]=f3c[ii];
  }
  else if(IMPLICITERRORNORM==4){
    if(f3d[0]<conv && f3d[1]<conv && f3d[2]<conv && f3d[3]<conv) passedconv=1;
    DLOOPA(ii) f3report[ii]=f3d[ii];
  }

  // get absolute error over all terms
  *errorabs=0.0;
  DLOOPA(jj) *errorabs     += fabs(f3report[jj]); // always full error.
  


  DLOOPA(ii){
    fpl3report[URAD0+ii]=f3report[ii];
  }


  // see if passed convergence test
  if(passedconv){
    if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g realdt=%g i=%d iter=%d DONE1 for conv=%g : f3report=%g %g %g %g\n",nstep,steppart,dt,realdt,ptrgeom->i,iter,conv,f3report[0],f3report[1],f3report[2],f3report[3]);
    return(1);
  }
  else{
    // report if didn't pass
    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTF1 (conv=%26.20g): uu: %26.20g %26.20g %26.20g %26.20g : uu0=%26.20g %26.20g %26.20g %26.20g\n",conv,uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uu0[URAD0],uu0[URAD1],uu0[URAD2],uu0[URAD3]);
      int iii;
      DLOOPA(iii) dualfprintf(fail_file,"iii=%d f1=%26.20g f1norm=%26.20g\n",iii,f1[iii],f1norm[iii]);
      dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3report[0],f3report[1],f3report[2],f3report[3]);
    }
    return(0);
  }

  return(0);
}


// calculating approximate Jacobian: dUresid(dUrad,G(Urad))/dUrad = dy(x)/dx
// then compute inverse Jacobian
static int get_implicit_iJ(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtypelocal, int whichcap, int itermode, FTYPE fracenergy, FTYPE dissmeasure, FTYPE impepsjac, int iter, FTYPE errorabs, FTYPE *dimfactU, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *fpl1, FTYPE *fpl1norm, FTYPE (*iJ)[NPR])
//static int get_implicit_iJ(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, FTYPE impepsjac, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pin, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NDIM])
{
  int ii,jj;
  FTYPE J[NDIM][NDIM],f2[NDIM],f2norm[NDIM];
  FTYPE f1[NDIM],f1norm[NDIM];

  DLOOPA(ii) f1[ii]=fpl1[URAD0+ii];
  DLOOPA(ii) f1norm[ii]=fpl1norm[URAD0+ii];

  // for scaling del's norm
  FTYPE sqrtgcov[NDIM];
  sqrtgcov[0]=1.0;
  SLOOPA(jj) sqrtgcov[jj] = sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]));


  // ensure pin doesn't get modified by del-shifts to get Jacobian, which can change primitives to order unity at high radiation gamma
  FTYPE pinjac[NPR];
  int pliter,pl;
  PLOOP(pliter,pl) pinjac[pl]=pp[pl];


  int failreturn;
  int fulljaciter=0;
  FTYPE FRACIMPEPSCHANGE=0.1;
  FTYPE del;
  FTYPE IMPEPSSTART=impepsjac;
  while(1){ // ensuring that Jacobian is non-singular if only because del too small (and then if singular, increase del)

    FTYPE localIMPEPS=IMPEPSSTART; // start with fresh del

    DLOOPA(jj){

      while(1){

        // when |URAD0|>>|URAD1|, then can't get better than machine error on URAD0, not URAD1, so using small del just for URAD1 makes no sense, so avoid above
        del = localIMPEPS*MAX(fabs(uup[jj+URAD0]), fabs(uup[URAD0]*sqrtgcov[jj]))  ;
          
        // offset uu (KORALTODO: How to ensure this doesn't have machine precision problems or is good enough difference?)
        uu[jj+URAD0]=uup[jj+URAD0]-del;
 
        // get dUresid for this offset uu
        int whichcall=2;
        //        PLOOP(pliter,pl) dualfprintf(fail_file,"BEFORE Jf1 jj=%d pl=%d uu=%21.15g pinjac=%21.15g\n",jj,pl,uu[pl],pinjac[pl]);
        int radinvmod=0;
        int goexplicit=0;
        int eomtypelocallocal=*eomtypelocal;
        FTYPE fpl2[NPR],fpl2norm[NPR];
        failreturn=f_implicit(iter,failreturnallowableuse, whichcall,showmessages,allowlocalfailurefixandnoreport, &eomtypelocallocal, whichcap, itermode, fracenergy, dissmeasure, &radinvmod, pinjac,uu0,uu,fracdtG*realdt,ptrgeom,q, fpl2,fpl2norm,&goexplicit);
        DLOOPA(ii){
          f2[ii]=fpl2[URAD0+ii];
          f2norm[ii]=fpl2norm[URAD0+ii];
        }
        
        if(failreturn){
          if(showmessages&& debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed: jj=%d.  Trying smaller localIMPEPS=%g (giving del=%g) to %g\n",jj,localIMPEPS,del,localIMPEPS*FRACIMPEPSCHANGE);
          localIMPEPS*=FRACIMPEPSCHANGE;
          // try making smaller until no error, unless doesn't work out
          // see if will be able to resolve differences
          int baddiff = fabs(uu[jj+URAD0]-uup[jj+URAD0])/(fabs(uu[jj+URAD0])+fabs(uup[jj+URAD0])) < 10.0*NUMEPSILON;
          if(localIMPEPS<10.0*NUMEPSILON || baddiff){
            // then probably can't resolve difference due to too small 
            if(failreturnallowableuse>=UTOPRIMGENWRAPPERRETURNFAILRAD){
              if(debugfail>=2) dualfprintf(fail_file,"Bad error: f_implicit_lab for f2 failed: jj=%d with localIMPEPS=%g (giving del=%g)\n",jj,localIMPEPS,del);
              return(1); // can't go below machine precision for difference else will be 0-0 and no reversion to do.
            }
            else{
              // instead of failing, allow radiation error, and restart process
              failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD; // just changes value in this function only.  Assume that once do this, applies to all further terms in Jacobian.
              localIMPEPS=IMPEPSSTART; // start with fresh del
            }
          }
        }// end if failreturn!=0
        else{
          // didn't fail
          break;
        }
      }


      // get Jacobian (uncentered, ok?  Probably actually best.  Don't want to go back along unknown trajectory in U that might lead to bad P(U))
      DLOOPA(ii) J[ii][jj]=(f2[ii] - f1[ii])/(uu[jj+URAD0]-uup[jj+URAD0]);

      //      DLOOPA(ii) dualfprintf(fail_file,"OLD: ii=%d jj=%d J=%26.21g : %26.21g %26.21g : %26.21g %26.21g\n",ii,jj,J[ii][jj], f2[ii],f1[ii],uu[jj+URAD0],uup[jj+URAD0]);


      if(debugfail>=2){
        DLOOPA(ii) if(showmessagesheavy || !isfinite(J[ii][jj])){
          dualfprintf(fail_file,"JAC: uu: %26.20g %26.20g %26.20g %26.20g : uup=%26.20g %26.20g %26.20g %26.20g (del=%26.20g localIMPEPS=%26.20g)\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3],del,localIMPEPS);
          dualfprintf(fail_file,"i=%d jj=%d f2: %26.20g %26.20g %26.20g %26.20g\n",ptrgeom->i,jj,f2[0],f2[1],f2[2],f2[3]);
          dualfprintf(fail_file,"JISNAN: %d %d : %26.20g : %26.20g %26.20g : %26.20g %26.20g\n",ii,jj,J[ii][jj],f2[ii],f1[ii],uu[jj+URAD0],uup[jj+URAD0]);
        }
      }

      // restore uu after getting changed by f_implicit_lab(f2)
      uu[jj+URAD0]=uup[jj+URAD0];
    }
    


    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTJAC: uu: %26.20g %26.20g %26.20g %26.20g : uup=%26.20g %26.20g %26.20g %26.20g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
      int iii,jjj;
      DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%26.20g\n",iii,jjj,J[iii][jjj]);
    }

    

    /////////////////////
    //
    //invert Jacobian
    //
    /////////////////////
    FTYPE iJsub[NDIM][NDIM];
    failreturn=inverse_44matrix(J,iJsub);

    // copy back inverse matrix from sub-version
    DLOOP(ii,jj){
      iJ[URAD0+ii][URAD0+jj]=iJsub[ii][jj];
    }
    {
      //      int iii,jjj; DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%26.20g iJsub=%26.20g\n",iii,jjj,J[iii][jjj],iJsub[iii][jjj]);
    }



    /////////////////////////////////////
    //
    // check if inversion was successful
    //
    /////////////////////////////////////
    if(failreturn){
      // try increasing localIMPEPS
      IMPEPSSTART/=FRACIMPEPSCHANGE;
      int condnotdiff;
      condnotdiff=IMPEPSSTART > MAXIMPEPS;

      if(condnotdiff){ // KORALTODO: But error relative to uu needs to be accounted for!
        if(debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed to be different enough from f1 and gave singular Jacobian: IMPEPSSTART=%g (giving del=%g)\n",IMPEPSSTART,del);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC1: uu: %26.20g %26.20g %26.20g %26.20g : uup=%26.20g %26.20g %26.20g %26.20g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%26.20g\n",iii,jjj,J[iii][jjj]);
        }
        return(1); // can't expect good derivative above ~0.3, so just return as failure of implicit method.
      }
      else{
        if(debugfail>=2) dualfprintf(fail_file,"inverse_44matrix(J,iJ) failed, trying IMPEPSSTART=%g :: ijk=%d %d %d\n",IMPEPSSTART,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC2: uu: %26.20g %26.20g %26.20g %26.20g : uup=%26.20g %26.20g %26.20g %26.20g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%26.20g\n",iii,jjj,J[iii][jjj]);
        }
      }
    }// end if failred to invert J
    else break; // good Jacobian

    // check if trying too many times to get Jacobian
    fulljaciter++;
    if(fulljaciter>MAXJACITER){
      // this is a catch in case bouncing back and forth between singular Jac and no inversion for P(U) to get f2
      if(debugfail>=2) dualfprintf(fail_file,"Failed to get inverse Jacobian with fulljaciter=%d with IMPEPSSTART=%g (giving del=%g)\n",fulljaciter,IMPEPSSTART,del);
      if(debugfail>=2 || showmessagesheavy){
        dualfprintf(fail_file,"POSTJAC3: uu: %g %g %g %g : uup=%g %g %g %g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
        int iii,jjj;
        DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%g\n",iii,jjj,J[iii][jjj]);
      }
      return(1);
    }
  }// end over ensuring Jacobian is non-singular for the given del


  return(0);

}




#define JDIFFONESIDED 0
#define JDIFFCENTERED 1




















// get dt for explicit sub-cyclings
static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother,  FTYPE *CUf, FTYPE *Gdpl, FTYPE chi, FTYPE *Gdplabs, struct of_geom *ptrgeom, FTYPE *dtsub)
{
  int jj;
  int pliter,pl;
  FTYPE idtsub0,idtsub;
  //
  FTYPE Umhd,Urad,Gmhd,Grad,iUmhd,iUrad;
  //
  FTYPE idtsubs,idtsubt;
  FTYPE idtsubmhd,idtsubrad;
  FTYPE Usmhd,Usrad,Gsmhd,Gsrad,iUsmhd,iUsrad;
  FTYPE Utmhd,Utrad,Gtmhd,Gtrad,iUtmhd,iUtrad;
  FTYPE Gddtpl[NPR];


  if(REMOVERESTMASSFROMUU!=2){
    dualfprintf(fail_file,"get_dtsub() assumes removed rest mass from UU so can compare G and U[UU]\n");
    myexit(9285345);
  }

  // KORALTODO: If Umhd is very small and dynamically unimportant, should probably ignore trying to subcycle, but at least limit number of sub-cycles.  May be more of a problem when deal with MHD in force-free magnetosphere.

  // NOTE: The timestep does not fllow NR1992 S19.2L on diffusion equation step size.  That applies if diffusion was part of flux calculation, not source term.  And it applies to the radiative velocity limiter for the advection's effective wave speed.
  // The relevant NR1992 is S16.6L on stiff source terms.

  // get G*dt that can be compared to Umhd or Urad
  FTYPE realdt=compute_dt(CUf,dt);

  // choose if using actual 4-force or absolutified 4-force
  if(0){
    PLOOP(pliter,pl) Gddtpl[pl] = Gdpl[pl]*realdt;
  }
  else{
    PLOOP(pliter,pl) Gddtpl[pl] = Gdplabs[pl]*realdt;
  }


  //////////
  //
  // get updated uu from dUother in case leads to different result
  FTYPE U0[NPR];
  PLOOP(pliter,pl) U0[pl]=UFSET(CUf,dt,Ui[pl],Uf[pl],dUother[pl],0.0);

  //  PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d U0=%g realdt=%g dt=%g Ui=%g Uf=%g dUother=%g\n",pl,U0[pl],realdt,dt,Ui[pl],Uf[pl],dUother[pl]);



  // get original U
  FTYPE U[NPR];
  PLOOP(pliter,pl) U[pl]=Ui[pl];




  // get smallest timestep for stiff source terms of 8 equations with a single source term vector.
  // Based upon NR 16.6.6L with removal of factor of two
  if(method==TAUSUPPRESS){
    // KORALTODO: \tau suppression not general enough because G\propto R\tau + \gamma R \tau + \gamma B, and further, this assumes T\sim R.  If T<<R, then suppression by \tau won't be enough to treat stiffness effect on fluid's T.
    // dR^t_t/R^t_t \sim c \gamma_{fluid} \chi dt = \sim c \gamma_{fluid} \tau dt/dx , but might as well just use first version
    

    if(0){
      // Older Olek method
      // use approximate dt along each spatial direction.  chi is based in orthonormal basis
      // get maximum \tau for all relevant directions
      FTYPE dxortho[NDIM],tautotdir[NDIM];
      SLOOPA(jj) dxortho[jj] = (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])));
      // KORALTODO: Should this only depend upon kappa and not kappaes?  Stiffness in source term comes from full chi, so seems to be full chi.
      SLOOPA(jj) tautotdir[jj] = chi * dxortho[jj];
      // only include relevant directions
      FTYPE taumax=SMALL+MAX(MAX(tautotdir[1]*N1NOT1,tautotdir[2]*N2NOT1),tautotdir[3]*N3NOT1);
      
      idtsub0=taumax/realdt;
    }
    else{
      // New Jon method (problem is this only makes sense in perfectly LTE.  If T gets high quickly, then G gets high before the density reacts.)
      FTYPE ucon0=q->ucon[TT]; // what enters G for dR^t_t/R^t_t from time part
      //      FTYPE ratchangeRtt=chi * ucon0 * realdt * 1.0; // 1.0 = c (1st term)
      FTYPE ratchangeRtt=SMALL+fabs(chi * ucon0 * ucon0 * realdt * 1.0); // 1.0 = c (2nd term with chi instead of kappaes to be sever and account for B-based term)
      idtsub0 = ratchangeRtt/realdt; // if ratchange=1, then in principle right at edge of big change.
      // this is like having a "source speed" of v_s\sim \tau \gamma^2 c and limiting the timestep so the source wave only reaches across a cell dxortho in time dt.

      //      dualfprintf(fail_file,"ucon0=%g chi=%g ratchangeRtt=%g idtsub=%g\n",ucon0,chi,ratchangeRtt,idtsub0);


    }


    // pre-ratio idtsub
    idtsub=idtsub0;



    //    dualfprintf(fail_file,"i=%d dtsub0=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);

    // account for case where effect on fluid is more than on radiation (where above would only account for effect on radiation)

    // first compare to original U
    jj=TT; Umhd=SMALL+fabs(U[UU+jj]);
    jj=TT; Urad=fabs(U[URAD0+jj]);
    idtsub=MAX(idtsub,idtsub0*Urad/Umhd);

    //    dualfprintf(fail_file,"i=%d dtsub1=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);
       
    // also compare against changed U=U0
    jj=TT; Umhd=SMALL+fabs(U0[UU+jj]);
    jj=TT; Urad=fabs(U0[URAD0+jj]);
    idtsub=MAX(idtsub,idtsub0*Urad/Umhd);

    //    dualfprintf(fail_file,"i=%d dtsub2=%g (realdt=%g)  Urad=%g Umhd=%g\n",ptrgeom->i,1/idtsub,realdt,Urad,Umhd);
       
  }
  // below is if Diffusion is part of source term as in Koral
  // source term should lead to small (<1/2) change in conserved quantities
  else if(method==SPACETIMESUBSPLITNONE){
    // merged space-time to avoid negligible total momentum with large update needing to be resolved.
    Umhd=Urad=Gmhd=Grad=0.0;
    DLOOPA(jj) Umhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Urad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Gmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Grad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    iUmhd=1.0/(fabs(Umhd)+SMALL);
    iUrad=1.0/(fabs(Urad)+SMALL);
    idtsub=SMALL+fabs(Gmhd*iUmhd);
    idtsub=MAX(idtsub,SMALL+fabs(Grad*iUrad));
    idtsub=sqrt(idtsub);

    //    if(1||realdt/(COURRADEXPLICIT/idtsub)>1.0) dualfprintf(fail_file,"UMHD: Umhdrad=%g %g : G=%g %g %g %g : Gmhdrad= %g %g :: iUmhdrad=%g %g ::: dtsub=%g realdt/dtsub=%g\n",Umhd,Urad,Gddtpl[UU],Gddtpl[U1],Gddtpl[U2],Gddtpl[U3],Gmhd,Grad,iUmhd,iUrad,COURRADEXPLICIT/idtsub,realdt/(COURRADEXPLICIT/idtsub));

  }
  else if(method==SPACETIMESUBSPLITTIME){
    // won't be efficient if v~0
    // if v<<1 and G is mid-range but still negligible, then dt will be incredibly small and code will halt if sub-cycling.
    Usmhd=Usrad=Gsmhd=Gsrad=0.0;
    Utmhd=Utrad=Gtmhd=Gtrad=0.0;
    SLOOPA(jj) Usmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Utmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Usrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Utrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Gsmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Gtmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Gsrad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Gtrad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    iUsmhd=1.0/(fabs(Usmhd)+SMALL);
    iUtmhd=1.0/(fabs(Utmhd)+SMALL);
    iUsrad=1.0/(fabs(Usrad)+SMALL);
    iUtrad=1.0/(fabs(Utrad)+SMALL);
    idtsubs=SMALL+fabs(Gsmhd*iUsmhd);
    idtsubs=MAX(idtsubs,SMALL+fabs(Gsrad*iUsrad));
    idtsubt=SMALL+fabs(Gtmhd*iUtmhd);
    idtsubt=MAX(idtsubt,SMALL+fabs(Gtrad*iUtrad));
    idtsub=MAX(idtsubs,idtsubt);
    idtsub=sqrt(idtsub);
  }
  else if(method==SPACETIMESUBSPLITALL){
    // won't be efficient if flow becomes grid-aligned or if v~0
    Usmhd=Usrad=Gsmhd=Gsrad=0.0;
    idtsub=0.0;
    DLOOPA(jj){
      Umhd = fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Urad = fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Gmhd = fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Grad = fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      iUmhd=1.0/(fabs(Umhd)+SMALL);
      iUrad=1.0/(fabs(Urad)+SMALL);
      idtsub=MAX(idtsub,SMALL+fabs(Gmhd*iUmhd));
      idtsub=MAX(idtsub,SMALL+fabs(Grad*iUrad));
    }
  }
  else if(method==SPACETIMESUBSPLITSUPERALL){
    // won't be efficient if flow becomes grid-aligned or if v~0 or if radiation neglibile contribution to fluid dynamics
    Usmhd=Usrad=Gsmhd=Gsrad=0.0;
    idtsub=0.0;
    DLOOPA(jj){
      Umhd = fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Urad = fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Gmhd = fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      Grad = fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
      iUmhd=1.0/(fabs(Umhd)+SMALL);
      iUrad=1.0/(fabs(Urad)+SMALL);
      idtsub=MAX(idtsub,SMALL+fabs(Gmhd*iUmhd));
      idtsub=MAX(idtsub,SMALL+fabs(Grad*iUrad));
    }
  }
  else if(method==SPACETIMESUBSPLITMHDRAD){
    // merged space-time to avoid negligible total momentum with large update needing to be resolved.
    Umhd=Urad=Gmhd=Grad=0.0;
    idtsub=0.0;
    DLOOPA(jj) Umhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Urad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Gmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    DLOOPA(jj) Grad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    iUmhd=1.0/(fabs(Umhd)+SMALL);
    iUrad=1.0/(fabs(Urad)+SMALL);
    idtsub=MAX(idtsub,SMALL+fabs(Gmhd*iUmhd));
    idtsub=MAX(idtsub,SMALL+fabs(Grad*iUrad));

    //    dualfprintf(fail_file,"UMHD: %g %g %g %g %g %g\n",Umhd,Urad,Gtot,iUmhd,iUrad);

  }
  else if(method==SPACETIMESUBSPLITTIMEMHDRAD){
    // won't be efficient if v~0
    // if v<<1 and G is mid-range but still negligible, then dt will be incredibly small and code will halt.
    Usmhd=Usrad=Gsmhd=Gsrad=0.0;
    Utmhd=Utrad=Gtmhd=Gtrad=0.0;
    SLOOPA(jj) Usmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Utmhd += fabs(U[UU+jj]*U[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Usrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Utrad += fabs(U[URAD0+jj]*U[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Gsmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Gtmhd += fabs(Gddtpl[UU+jj]*Gddtpl[UU+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    SLOOPA(jj) Gsrad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    jj=TT;     Gtrad += fabs(Gddtpl[URAD0+jj]*Gddtpl[URAD0+jj]*ptrgeom->gcon[GIND(jj,jj)]);
    iUsmhd=1.0/(fabs(Usmhd)+SMALL);
    iUtmhd=1.0/(fabs(Utmhd)+SMALL);
    iUsrad=1.0/(fabs(Usrad)+SMALL);
    iUtrad=1.0/(fabs(Utrad)+SMALL);
    idtsub=SMALL;
    idtsub=MAX(idtsub,SMALL+fabs(Gsmhd*iUsmhd));
    idtsub=MAX(idtsub,SMALL+fabs(Gsrad*iUsrad));
    idtsub=MAX(idtsub,SMALL+fabs(Gtmhd*iUtmhd));
    idtsub=MAX(idtsub,SMALL+fabs(Gtrad*iUtrad));
  }


  
  // what to return
  *dtsub=COURRADEXPLICIT/idtsub;



  //  dualfprintf(fail_file,"*dtsub=%g idtsub=%g method=%d\n",*dtsub,idtsub,method);
 
}

#define EXPLICITFAILEDBUTWENTTHROUGH -2
#define EXPLICITNOTNECESSARY -1
#define EXPLICITNOTFAILED 0 // should stay zero
#define EXPLICITFAILED 1 // should stay one


#define GETADVANCEDUNEW0FOREXPLICIT 1 // Use this to check if single explicit step was really allowable, but get_dtsub() already uses advanced U.  But chi will be not updated for fluid dUriemann update, so still might want to do this (with proper code changes) in order to get chi good.

// Based upon size of Gd, sub-cycle this force.
// 1) calc_Gd()
// 2) locally set dtsub~dt/\tau or whatever it should be
// 3) update T^t_\nu and R^t_\nu
// 4) U->P locally
// 5) repeat.
// Only change dUcomp, and can overwrite prnew, Unew, and qnew since "prepare" function isolated original values already
static int source_explicit(int whichsc, int whichradsourcemethod, int methoddtsub,int *eomtype,
                           void (*sourcefunc)(int methoddtsub, FTYPE *pr, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gpl, struct of_geom *ptrgeom, FTYPE *dtsub),
                           FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother, FTYPE (*dUcomp)[NPR])
{
  int pliter, pl;

  int eomtypelocal=*eomtype; // default

  ////////////////
  //
  // SETUP LOOPS
  //
  ///////////////
  int showmessages=0;
  int showmessagesheavy=0;
  int allowlocalfailurefixandnoreport=0; // need to see if any failures.
  struct of_newtonstats newtonstats;
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.


  //  if(1||nstep>=800){
  //    showmessages=showmessagesheavy=1;
  //  }

  ////////////////
  //
  // SETUP U and P and q
  //
  ///////////////

  FTYPE pb0[NPR],Uiin0[NPR];
  FTYPE prnew[NPR],Unew[NPR],Unew0[NPR];
  FTYPE prforG[NPR];
  struct of_state q0,qnew;
  FTYPE Gpl[NPR];
  FTYPE chi;

  // backup pb, Uiin, and q and setup "new" versions to be iterated
  PLOOP(pliter,pl) prnew[pl]=pb0[pl]=pb[pl];
  PLOOP(pliter,pl) Unew[pl]=Uiin0[pl]=Uiin[pl];
  qnew=q0=*q;


  // get updated U (try getting full update)
  FTYPE fracdtuu0=1.0; // try full uu0 at first
  PLOOP(pliter,pl) Unew[pl]=Unew0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);

  // if reversion from implicit, then no choice but to push through CASE radiation errors and hope the reductions there are ok.  Would be worse to have no reversion solution!
  int pushthroughraderror=0;
  if(whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
    pushthroughraderror=1;
  }



  // Get prnew(Unew) using largest fracdtuu0 possible in order to get realistic estimate of dtsub
  // used to use this as starting point for U, but that wasn't consistent with explicit stepping
  if(GETADVANCEDUNEW0FOREXPLICIT){
    //////////////
    //
    // Get good Unew0 that has P(Unew0) solution
    //
    //////////////
    while(1){ // loop bounded by fracdtuu0 becoming too small
      // Get pnew from Unew
      // OPTMARK: Should optimize this to  not try to get down to machine precision
      // initialize counters
      newtonstats.nstroke=newtonstats.lntries=0;
      int doradonly=0;
      int whichcap=CAPTYPEFIX1;
      eomtypelocal=*eomtype; // re-default
      int radinvmod=0;
      FTYPE dissmeasure=-1.0; // assume ok to try energy
      int failutoprim=Utoprimgen_failwrapper(doradonly,&radinvmod,showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Unew, q, ptrgeom, dissmeasure, prnew, &newtonstats);

      if(failutoprim){

        if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
          // then ok to be here
        }
        else{
          // if here, then must be doing implicit checks, so return that should just do implicit instead of any more expensive calculations here.
          return(EXPLICITFAILED);
        }

        // backing off dU
        fracdtuu0*=0.5;

        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"Backing off fracdtuu0=%g\n",fracdtuu0);

        if(fracdtuu0<NUMEPSILON){
          if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"In explicit, backed-off to very small level of fracdtuu0=%g, so must abort\n",fracdtuu0);
          if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
            // just use initial Unew0 then
            fracdtuu0=0.0;
            break;
          }
          else return(EXPLICITFAILED);
        }
        else{
          // recompute use of full dU so Unew0 is updated
          PLOOP(pliter,pl) Unew[pl]=Unew0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
          // reset prnew
          PLOOP(pliter,pl) prnew[pl]=pb0[pl]=pb[pl];
        }
            
      }
      else{
        // then found good Unew0
        break;
      }
    }// end while trying to get good Unew0

  }

  //////////
  //
  // Get future force and dtsub so don't overestimate dtsub for first sub-cycle steps and end-up possibly jumping too far
  //
  //////////
  // get prforG to compute G and chi for calc_dtsub() to ensure future update doesn't have radically different opacity and so 4-force and underpredict that implicit or sub-cycles are needed.
  PLOOP(pliter,pl) prforG[pl]=prnew[pl];

 
  // get dtsubforG (don't use Gpl from this)
  FTYPE dtsubforG;
  sourcefunc(methoddtsub, prforG, Uiin, Ufin, dUother, CUf, Gpl, ptrgeom, &dtsubforG);
  //  dualfprintf(fail_file,"dtsubforG=%g\n",dtsubforG);

  if(!(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT)){
    // then if sub-cycling, really want to start with beginning pb so consistently do sub-steps for effective flux force and full-pl fluid force in time.
    // But if end-up doing just one explcit step, then probably wanted to use time-advanced prnew as estimate.  That gives more stable result.
    // then prforG is only used to ensure not getting bad guess for whether *should* sub-cycle.
    PLOOP(pliter,pl) prnew[pl]=pb[pl];
  }


  ////////////////
  //
  // SETUP explicit sub-cycle LOOP
  //
  ///////////////
  FTYPE dttrue=0.0,dtcum=0.0;  // cumulative sub-cycle time
  FTYPE dtdiff;
  FTYPE dtsub,dtsubold,dtsubuse;
  FTYPE realdt=compute_dt(CUf,dt);
  FTYPE fracdtG;

  int jj;
  FTYPE Gplprevious[NPR]={0}, sourcepl[NPR];


  // initialize source update
  PLOOP(pliter,pl) sourcepl[pl] = 0;

  //////////////
  //
  // explicit source LOOP
  //
  //////////////
  int itersub=0;
  int done=0;
  while(1){

  
    // get 4-force for full pl set
    PLOOP(pliter,pl) Gplprevious[pl]=Gpl[pl];
    // get Gpl and dtsub for sub-cycling
    sourcefunc(methoddtsub, prnew, Uiin, Ufin, dUother, CUf, Gpl, ptrgeom, &dtsub);
    if(itersub==0)  PLOOP(pliter,pl) Gplprevious[pl]=Gpl[pl]; // constant interpolation rather than just using zero for midpoint method
    if(itersub==0 && dtsub>dtsubforG) dtsub=dtsubforG; // ensure initial sub-stepping is not too large due to large state changes not accounted for yet.  Assuming slow safety factor growth in dtsub can occur from then on and that's ok since will catch updated state change to some fraction.
    //    dualfprintf(fail_file,"itersub=%d dtsub=%g\n",itersub,dtsub);


    // if no solution for implicit and come to explicit, then failure can manifest as T large and then Gpl->nan or inf.  Must fail this scenario.
    // This can happen even when have gotten quite close to end of step, but just no actually solution for the accurate value of U0 and G
    PLOOP(pliter,pl) if(!isfinite(Gpl[pl])) return(EXPLICITFAILED);


    if(showmessagesheavy&&debugfail>=2){
      PLOOP(pliter,pl) dualfprintf(fail_file,"SOURCE: pl=%d Gpl=%g dtsub=%g realdt=%g prnew=%g Uiin=%g Ufin=%g dUother=%g\n",pl,Gpl[pl],dtsub,realdt,prnew[pl],Uiin[pl],Ufin[pl],dUother[pl]);
    }

    if(whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLECHECKSFROMIMPLICIT){
      // then just check if need sub-cycles, and exit if so
      if(realdt>dtsub) return(EXPLICITFAILED);
      // else can do explicit step (or sub-cycles if happen to switch to them) and can avoid implicit step
    }
    // if still here, then even with implicit checks, doing either explicit or sub-cycle explicit stepping

    //////////////////
    // get fracdtG
    //////////////////
      
    if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT){
      fracdtG=1.0;
    }
    else{

      if(realdt/dtsub>MAXSUBCYCLESFAIL){
        // then some major issue
        if(debugfail>=2) dualfprintf(fail_file,"MAJOR explicit issue: realdt=%g dtsub=%g and aborted assuming no real solution and that best approximation is no force.\n");
        return(EXPLICITFAILED);
      }
      else if(realdt/dtsub>MAXSUBCYCLES && whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
        // Use MAXSUBCYCLES if otherwise was trying implicit.
        // Impractical to assume if really revert to explicit (or really trying to use it) then rarely occurs or want to solve for actual solution, so do all needed sub-cycles!
        // Semi-required to limit number of cycles for non-simple "methoddtsub" procedure that can produce arbitrarily small dtsub due to (e.g.) momentum term or something like that.
        // NOTE: For high \tau, rad velocity entering chars goes like 1/\tau, so that timestep is higher.  But dtsub remains what it should be for explicit stepping, and so in high-tau case, explicit steps required per actual step goes like \tau^2.  So "kinda" ok that takes long time for explicit sub-stepping since ultimately reaching longer time.
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"itersub=%d dtsub very small: %g with realdt=%g and only allowing MAXSUBCYCLES=%d subcycles, so limit dtsub: ijk=%d %d %d\n",itersub,dtsub,realdt,MAXSUBCYCLES,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        dtsub=realdt/(FTYPE)MAXSUBCYCLES;
      }
      else if(NUMEPSILON*dtsub>=realdt && itersub==0){
        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit not necessary\n");
        // then no need for source term at all
        return(EXPLICITNOTNECESSARY);
      }


      if(itersub==0) dtsubold=dtsubuse=dtsub;
      else{
        // override if dtsub is larger than realdt, indicating really done with iterations and reached some equilibrium, so no longer necessary to check vs. dtsubold
        // No, too speculative.
        //        if(dtsub>realdt) dtsub=realdt;
        
        // ensure don't change step too fast.  Sometimes first guess for dtsub can be small, and very next iteration suggests very large.  Not trustable, so stay slow.
        if(dtsub>dtsubold*(1.0+MAXEXPLICITSUBSTEPCHANGE)) dtsubuse=dtsubold*(1.0+MAXEXPLICITSUBSTEPCHANGE);
        else dtsubuse=dtsub;

        // need to compare with previous actual dt
        dtsubold=dtsubuse;
      }
        
      // time left to go in sub-cycling
      dtdiff=MAX(realdt-dtcum,0.0);
      dttrue=MIN(dtsubuse,dtdiff);
      fracdtG=MIN(1.0,dttrue/realdt); // expect fraction of G that can be handled is similar to what sub-cycle requires for stable explicit stepping

      
      if(showmessagesheavy&&debugfail>=2){
        dualfprintf(fail_file,"DoingSUBCYCLE: itersub=%d : dtsub=%g dtsubuse=%g dtdiff=%g dttrue=%g dtcum=%g realdt=%g fracdtG=%g ijk=%d %d %d\n",itersub,dtsub,dtsubuse,dtdiff,dttrue,dtcum,realdt,fracdtG,ptrgeom->i,ptrgeom->j,ptrgeom->k);
      }
      if(debugfail>=2 && (1||showmessages)&&(1.0/fracdtG>MAXSUBCYCLES && itersub==0)){ // have to use itersub==0 since already might be done with itersub==1 and then fracdtG=inf (but itersub=0 might be using bad force)
        dualfprintf(fail_file,"DoingLOTSofsub-cycles: ijk=%d %d %d  1/fracdtG=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,1.0/fracdtG);
      }
    }// end else if not explicit
      


    // add to final result if starting point is full U0
    // forward step integration
    //    PLOOP(pliter,pl) sourcepl[pl] += Gpl[pl]*fracdtG;
    // Trapezoidal Rule (midpoint method):
    PLOOP(pliter,pl) sourcepl[pl] += 0.5*(Gplprevious[pl]+Gpl[pl])*fracdtG;



    ///////////////
    //
    // see if done
    //
    ///////////////
    if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT){
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit done\n");
      break;
    }
    else if(whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLECHECKSFROMIMPLICIT){
      if(dtcum>=realdt){
        // then done
        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit sub-cycle done\n");
        break;
      }
      else{
        // then keep sub-cycling _or_ getting balance of U and G
      }
    }
    else{
      if(fracdtG!=1.0){
        // if here, then must be doing implicit checks, so return that should just do implicit instead of sub-cycling.
        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit sub-cycle can't be done.\n");
        return(EXPLICITFAILED);
      }
      else{
        // then implicit testing, and had to do step, but only 1 step, so consider success.
        // still need to fill dUcomp[] below
        if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit success, just one step.\n");
        break;
      }
    }

    ///////
    // If not done, get prnew from Unew for next step
    ///////

    // get new Unew0
    FTYPE xint=(dtcum/realdt);
    // NOTE: Below interpolates from Unew0 could use up to final, but not consistent with explicit stepping and leads to erroneous 0-force results.
    //    FTYPE fakefracdtuu0=fracdtuu0*(1.0-xint) + 1.0*(xint);
    // NOTE: Below interpolates from step's true starting Unew0 to step's final Unew0 assuming linear interpolation between -- *most* consistent with explicit stepping!!
    FTYPE fakefracdtuu0=1.0*(xint);
    // NOTE: Below sticks to the Unew0 that could use, but not consistent with explicit stepping and leads to erroneous 0-force results.
    //    FTYPE fakefracdtuu0=fracdtuu0;
    FTYPE tempdt= fakefracdtuu0*dt; // uses dt here, because below UFSET() computes "realdt" using CUf internally
    PLOOP(pliter,pl) Unew0[pl]=UFSET(CUf,tempdt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
    

    // get new Unew using 1) current Unew0 (so Unew updates a bit towards final Unew0 as if fracdtuu0=1) and 2) cumulative 4-force so far
    PLOOP(pliter,pl) Unew[pl] = Unew0[pl] + sourcepl[pl] * realdt;

    // get prnew(Unew)
    newtonstats.nstroke=newtonstats.lntries=0;
    int doradonly=0;
    eomtypelocal=*eomtype; // re-default
    int whichcap=CAPTYPEFIX1;
    int radinvmod=0;
    FTYPE dissmeasure=-1.0; // assume ok to try energy
    int failutoprim=Utoprimgen_failwrapper(doradonly,&radinvmod,showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Unew, q, ptrgeom, dissmeasure, prnew, &newtonstats);
    // push through inversion failure if just radiation inversion failure since have local fixups that can be ok or even recovered from.  Bad to just stop if doing reversion from implicit.
    if(pushthroughraderror==0 && failutoprim==UTOPRIMGENWRAPPERRETURNFAILRAD || pushthroughraderror==1 && failutoprim==UTOPRIMGENWRAPPERRETURNFAILMHD){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"BAD: Utoprimgen_wrapper() failed during explicit sub-stepping.  So sub-cycling failed.\n");
      return(EXPLICITFAILED);
    }


    // DEBUG:
    if(showmessagesheavy &&debugfail>=2) PLOOP(pliter,pl) dualfprintf(fail_file,"POSTEXSTEP: pl=%2d Unew0=%21.15g Unew=%21.15g sourcepl*realdt=%21.15g fakefracdtuu0=%g\n",pl,Unew0[pl],Unew[pl],sourcepl[pl]*realdt,fakefracdtuu0);


    // step      
    dtcum += realdt*fracdtG; // cumulative true time
    itersub++;


  }// done looping


  ////////////
  //
  // apply 4-force as update in dUcomp[][]
  // only changed this quantity, none of other among function arguments
  //
  ////////////
  PLOOP(pliter,pl) dUcomp[whichsc][pl] += sourcepl[pl];

  // save better guess for later inversion from this inversion
  PLOOP(pliter,pl) pb[pl]=prnew[pl];

  // save eomtype settled on
  *eomtype=eomtypelocal;


  return(EXPLICITNOTFAILED);
  
}








// General radiation source term calculation
// NOTE: source_explicit() takes as first argument a form of function like general koral_source_rad_calc() .  It doesn't have to be just used for radiation.
// NOTE: koral_source_rad_implicit() currently only works for radiation where only 4 equations involved since 4-force of rad affects exactly mhd.  So only invert 4x4 matrix.
// For recursion of other consistencies, should keep koral_source_rad() same function arguments as explicit and implicit functions.  Once make koral_source_rad() general, can use this function as general source function instead of it getting called just for radiation.
int koral_source_rad(int whichradsourcemethod, FTYPE *piin, FTYPE *pb, FTYPE *pf, int *didreturnpf, int *eomtype, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother, FTYPE (*dUcomp)[NPR])
{
  int pliter,pl;
  int showmessages=0; // 0 ok if not debugging and think everything works.
  int showmessagesheavy=0;

  //  if(1||nstep>=800){
  //    showmessages=showmessagesheavy=1;
  //  }


  // make code use "orig" values that below can modify, so return preserves no changes to these no matter what internal functions do.
  // save pb, Uiin, Ufin, and q to avoid being modified upon return
  FTYPE pborig[NPR],pforig[NPR],piinorig[NPR],Uiinorig[NPR],Ufinorig[NPR];
  struct of_state qorigmem;
  struct of_state *qorig=&qorigmem;
  PLOOP(pliter,pl){
    pborig[pl]=pb[pl];
    pforig[pl]=pf[pl];
    piinorig[pl]=piin[pl];
    Uiinorig[pl]=Uiin[pl];
    Ufinorig[pl]=Ufin[pl];
  }
  *qorig=*q;


  /////////////////
  //
  // Check energy density to see if even any radiation or will be any radiation
  // Note, can't just check size of Erf, because in non-LTE, B can be large even if E is not.
  //
  /////////////////



  /////////////////
  //
  // EXPLICIT TYPES
  //
  /////////////////
  if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLECHECKSFROMIMPLICIT){

    int methoddtsub;
    // SPACETIMESUBSPLITMHDRAD doesn't work -- generates tons of noise in prad1 with COURRADEXPLICIT=0.2, and was asymmetric in x.
    methoddtsub=TAUSUPPRESS; // forced -- only method that is efficient and effective and noise free at moderate optical depths.
    //    methoddtsub=SPACETIMESUBSPLITNONE;
    //    methoddtsub=SPACETIMESUBSPLITTIME;


    int whichsc = RADSOURCE;
    // try explicit (with or without sub-cycling)
    //    dualfprintf(fail_file,"Trying explicit: whichradsourcemethod=%d\n",whichradsourcemethod);
    int failexplicit=source_explicit(whichsc, whichradsourcemethod,methoddtsub,eomtype,koral_source_dtsub_rad_calc,pborig, piinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);
    if(failexplicit==EXPLICITFAILED){
      if(whichradsourcemethod==SOURCEMETHODEXPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITREVERSIONFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
        // still do explicit anyways, since best can do with the choice of method -- will fail possibly to work if stiff regime, but ok in non-stiff.
        // assume nothing else to do, BUT DEFINITELY report this.
        if(debugfail>=2) dualfprintf(fail_file,"BAD: explicit failed: ijk=%d %d %d : whichradsourcemethod=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,whichradsourcemethod);
        *didreturnpf=0;
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE3A; // must set as failure in case can fixup.
        return(EXPLICITFAILEDBUTWENTTHROUGH);
      }
      else if(whichradsourcemethod==SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLECHECKSFROMIMPLICIT){
        // tells that explicit didn't work for implicit checks
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"explicit failed for implicit check. ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
        *didreturnpf=0;
        return(EXPLICITFAILED);
      }
      else{
        // tells that explicit didn't work
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"explicit failed. ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE3B; // must set as failure in case can fixup.
        *didreturnpf=0;
        return(EXPLICITFAILED);
      }
    }
    else if(failexplicit==EXPLICITNOTNECESSARY){
      // then don't need any source term
      if(debugfail>=2) dualfprintf(fail_file,"ODD: explicit found not necessary while implicit failed: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
      *didreturnpf=0;
      return(0);
    }

    //else explicit succeeded, so just return
    if(whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLE || whichradsourcemethod==SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT){
      // if sub-cycled, then have better pf than pb assumed saved in pborig[].
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"explicit didn't fail. ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
      PLOOP(pliter,pl) pf[pl]=pborig[pl];
      *didreturnpf=1;
    }
    return(EXPLICITNOTFAILED);
  }
  /////////////////
  //
  // IMPLICIT TYPES (implicit function should handle pflag (rad or gas) error issues), but still leave explicit as handling error mode here
  //
  /////////////////
  else if(whichradsourcemethod==SOURCEMETHODIMPLICIT){
 
    int failimplicit=koral_source_rad_implicit(eomtype, pborig, pforig, piinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dissmeasure, dUother, dUcomp);
 
    if(failimplicit>0){
      if(IMPLICITREVERTEXPLICIT){ // single level recusive call (to avoid duplicate confusing code)
        // assume if revert from implicit, then need to do sub-cycles
        int failexplicit=koral_source_rad(SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dissmeasure, dUother, dUcomp);
        if(failexplicit==EXPLICITFAILED){
          // nothing else to revert to, but just continue and report
          *didreturnpf=0;
          if(debugfail>=2) dualfprintf(fail_file,"BAD: explicit failed while implicit failed: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
          GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE1A; // must set as failure in case can fixup.
          return(0);
        }
        else if(failexplicit==EXPLICITNOTNECESSARY){
          // then don't need any source term
          if(debugfail>=2) dualfprintf(fail_file,"ODD: explicit found not necessary while implicit failed: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
          GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE1B; // must set as failure in case can fixup.
          *didreturnpf=0;
          return(0);
        }
        else if(failexplicit==EXPLICITFAILEDBUTWENTTHROUGH){
          // then had issues, but nothing else can do.
          if(debugfail>=2) dualfprintf(fail_file,"HMM: explicit found necessary and had problems while implicit failed: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
          GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE2A; // must set as failure in case can fixup.
          *didreturnpf=0;
          return(0);
        }
        else{
          // if sub-cycled, then have better pf than pb assumed saved in pborig[].
          if(debugfail>=2) dualfprintf(fail_file,"GOOD: explicit worked while implicit failed (%d): ijk=%d %d %d\n",failexplicit,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          PLOOP(pliter,pl) pf[pl]=pborig[pl];
          *didreturnpf=1;
          return(0);
        }
      }// end if reverting to explicit
      else{
        *didreturnpf=0;
        if(debugfail>=2) dualfprintf(fail_file,"BAD: implicit failed and didn't choose to revert: ijknstepsteppart=%d %d %d %ld %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart);
        //        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE2B; // must set as failure in case can fixup. // NO, leave implicit call as handling error mode.
        return(0);        
      }
    }// end if failed to do implicit
    else if(failimplicit==0){
      // no failure in implicit, then just return
      // and if did implicit, then better pf guess
      //PLOOP(pliter,pl) pf[pl]=pborig[pl];
      *didreturnpf=0;
    }
    else{
      *didreturnpf=0;
      // e.g. if failimplicit==FAILRETURNGOEXPLICIT, then aborted implicit and letting trivial explicit operate.
    }
    //    if(debugfail>=2) dualfprintf(fail_file,"Good: Imlicit good.: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    return(0);

  }
  /////////////////
  //
  // IMPLICIT WITH EXPLICIT CHECK TYPES
  //
  /////////////////
  else if(whichradsourcemethod==SOURCEMETHODIMPLICITEXPLICITCHECK){

    // try explicit (or see if no source at all required)
    // Just check using explicit method, since if sub-cycles required then should just do implicit
    int failreturn=koral_source_rad(SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dissmeasure, dUother, dUcomp);

    // determine if still need to do implicit
    // don't set didreturnpf since already was set
    int doimplicit;
    if(failreturn==EXPLICITFAILED || failreturn==EXPLICITFAILEDBUTWENTTHROUGH || GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)>UTOPRIMRADNOFAIL){
      doimplicit=1;
      GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL; // reset and let implicit set this
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Tried explicit step, but wasn't good choice or failed: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
      // don't return until do implicit
    }
    else if(failreturn==EXPLICITNOTNECESSARY){
      // then no source at all required
      doimplicit=0;
      return(0);
    }
    else{
      doimplicit=0;
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Was able to take explicit step: ijk=%d %d %d : failreturn=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,failreturn);
      return(0);
    }
  

    if(doimplicit){
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"NOTE: Had to take implicit step: %d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);

      // one-deep recursive call to implicit scheme
      return(koral_source_rad(SOURCEMETHODIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dissmeasure, dUother, dUcomp));
    }// end if doimplicit==1

  }
  /////////////////
  //
  // NO SOURCE TYPE
  //
  /////////////////
  else if(whichradsourcemethod==SOURCEMETHODNONE){
    // then no source applied even if required
    *didreturnpf=0;
    return(0);
  }
  /////////////////
  //
  // UNKNOWN TYPE
  //
  /////////////////
  else{

    dualfprintf(fail_file,"3 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
    myexit(18754363);

  }

  // KORALTODO: SUPERGODMARK: Need to add NR 2007 page940 17.5.2L StepperSie method here as higher-order alternative if 1st order Newton breaks


  return(0);

}



//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************
//absorption
void calc_kappa(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa)
{

  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappa()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  *kappa = calc_kappa_user(rho,T,xx,yy,zz);
  //  dualfprintf(fail_file,"kappaabs=%g\n",*kappa);
}

//scattering
void calc_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappaes)
{  
  extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappaes()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  *kappaes = calc_kappaes_user(rho,T,xx,yy,zz);
  //  dualfprintf(fail_file,"kappaes=%g\n",*kappa);
}

// get \chi = \kappa_{abs} + \kappa_{es}
void calc_chi(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *chi)
{
  FTYPE kappa,kappaes;
  calc_kappa(pr,ptrgeom,&kappa);
  calc_kappaes(pr,ptrgeom,&kappaes);
  
  *chi=kappa+kappaes;
}

// get \kappa_{abs} and \kappa_{es}
static void calc_kappa_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa, FTYPE *kappaes, FTYPE *Tgas)
{
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappa()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u); // KORALNOTE: Currently, primary location where Tgas is computed for speed purposes.
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  *kappa = calc_kappa_user(rho,T,xx,yy,zz);
  *kappaes = calc_kappaes_user(rho,T,xx,yy,zz);
  *Tgas = fabs(T) + TEMPMIN;
  //  dualfprintf(fail_file,"kappaabs=%g kappaes=%g\n",*kappa,*kappaes);
}

// get G_\mu
static void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G, FTYPE *Tgas, FTYPE* chieffreturn, FTYPE *Gabs)
{
  calc_Gu(pp, ptrgeom, q, G, Tgas, chieffreturn,Gabs);
  indices_21(G, G, ptrgeom);
}



// get 4-force for all pl's
void koral_source_rad_calc(int computestate, int computeentropy, FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Gdpl, FTYPE *Gdplabs, FTYPE *chi, FTYPE *Tgas, struct of_state *q)
{
  int jj;
  int pliter,pl;
  FTYPE Gd[NDIM],Gdabs[NDIM];
  struct of_state qlocal;
  FTYPE chilocal,Tgaslocal;

  if(q==NULL){ q=&qlocal; computestate=1; }
  if(chi==NULL) chi=&chilocal;
  if(Tgas==NULL) Tgas=&Tgaslocal;


  // no, thermodynamics stuff can change since MHD fluid U changes, so must do get_state() as above
  //  get_state_uconucovonly(pr, ptrgeom, q);
  //  get_state_uradconuradcovonly(pr, ptrgeom, q);
  if(computestate) get_state(pr,ptrgeom,q);

  calc_Gd(pr, ptrgeom, q, Gd, Tgas, chi, Gdabs);

  PLOOP(pliter,pl) Gdpl[pl] = 0.0;
  // equal and opposite forces on fluid and radiation due to radiation 4-force
  // sign of G that goes between Koral determination of G and HARM source term (e.g. positive \lambda is a cooling of the fluid and heating of the photons, and gives G_t>0 so -G_t<0 and adds to R^t_t such that R^t_t - G_t becomes more negative and so more photon energy density)
  // That is, equation is d_t R^t_t + Gdpl = 0
#define SIGNGD (1.0)
  // keep SIGNGD as 1.0.  Just apply signgd2 in front of Gdpl in other places.
  // sign fixed-linked as + for URAD0 case.
  DLOOPA(jj) Gdpl[UU+jj]     = -SIGNGD*Gd[jj];
  DLOOPA(jj) Gdpl[URAD0+jj]  = +SIGNGD*Gd[jj];

  if(Gdplabs!=NULL){
    PLOOP(pliter,pl) Gdplabs[pl] = 0.0;
    DLOOPA(jj) Gdplabs[UU+jj]     = Gdabs[jj];
    DLOOPA(jj) Gdplabs[URAD0+jj]  = Gdabs[jj];
  }

#if(DOENTROPY!=DONOENTROPY && ENTROPY!=-100)
  if(computeentropy){
    pl=ENTROPY;
    // The equation is (1/\sqrt{-g})*d_\mu(\sqrt{-g} s\rho_0 u^\mu) + Gdpl[mu].ucon[mu] = 0
    FTYPE Gdplentropycontribs[NDIM];
    // -Gdpl[UU+jj] is so heating (so lowering of T^t_t to be more negative) implies increases entropy.
    // assumes Gpl includes kappa already with rho so that Gpl is energy per unit volume per unit time.  Dividing by T (energy) gives a dimensionless thing (entropy) per unit volume.
    DLOOPA(jj) Gdplentropycontribs[jj] = (-Gdpl[UU+jj])*(q->ucon[jj])/(*Tgas);

    Gdpl[pl] = 0.0;
    DLOOPA(jj) Gdpl[pl] += Gdplentropycontribs[jj];

    if(Gdplabs!=NULL){
      Gdplabs[pl] = 0.0;
      DLOOPA(jj) Gdplabs[pl] += fabs(Gdplentropycontribs[jj]);
    }
  }
#endif



}


// get 4-force and dtsub for all pl's
static void koral_source_dtsub_rad_calc(int method, FTYPE *pr, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gdpl, struct of_geom *ptrgeom, FTYPE *dtsub)
{
  FTYPE Gdplabs[NPR];
  FTYPE chi,Tgas;
  struct of_state q;

  int computestate=1;
  int computeentropy=1;
  koral_source_rad_calc(computestate,computeentropy,pr,ptrgeom,Gdpl,Gdplabs,&chi,&Tgas,&q);

  if(dtsub!=NULL){
    // then assume expect calculation of dtsub
    get_dtsub(method, pr, &q, Ui, Uf, dUother, CUf, Gdpl, chi, Gdplabs, ptrgeom, dtsub);
  }
  // else "method" can be anything and it doesn't matter


}

// compute G^\mu 4-force
static void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu, FTYPE *Tgas, FTYPE* chieffreturn, FTYPE *Gabs) 
{
  int i,j,k;
  
  //radiative stress tensor in the lab frame
  FTYPE Rij[NDIM][NDIM];

  //this call returns R^i_j, i.e., the first index is contra-variant and the last index is co-variant
  mhdfull_calc_rad(pp, ptrgeom, q, Rij);

  //the four-velocity of fluid in lab frame
  FTYPE *ucon,*ucov;
  ucon = q->ucon;
  ucov = q->ucov;
  
 
  // get opacities
  FTYPE kappa,kappaes;
  calc_kappa_kappaes(pp,ptrgeom,&kappa,&kappaes,Tgas);

  // get cooling rate
  FTYPE lambda;
  calc_rad_lambda(pp, ptrgeom, kappa, kappaes,*Tgas, &lambda);

  // get chi
  FTYPE chi=kappa+kappaes;
  
  // compute contravariant four-force in the lab frame
  
  //R^a_b u_a u^b
  FTYPE Ruu=0.; DLOOP(i,j) Ruu+=Rij[i][j]*ucov[i]*ucon[j];



  FTYPE Ru,term1,term2,term3;
  DLOOPA(i){
    Ru=0.; DLOOPA(j) Ru+=Rij[i][j]*ucon[j];

    // group by independent terms
    term1 = -(kappa*Ru + lambda*ucon[i]);
    term2 = -kappaes*(Ru + Ruu*ucon[i]);

    // actual source term
    //    Gu[i]=-chi*Ru - (kappaes*Ruu + lambda)*ucon[i];
    Gu[i] = term1 + term2;
    
    // absolute magnitude of source term that can be used for estimating importance of 4-force relative to existing conserved quantities to get dtsub.  But don't split kappa terms because if those cancel then physically no contribution.
    Gabs[i] = fabs(term1) + fabs(term2);

#if(0)
    // DEBUG:
    if(ptrgeom->i==3 && ptrgeom->j==26){
      dualfprintf(fail_file,"i=%d term1=%g term2=%g kappa=%g lambda=%g kappaes=%g ucon=%g Gu=%g Gabs=%g\n",i,term1,term2,kappa,lambda,kappaes,ucon[i],Gu[i],Gabs[i]);
    }
#endif

  }

  // really a chi-effective that also includes lambda term in case cooling unrelated to absorption
  *chieffreturn=chi + lambda/(SMALL+fabs(pp[PRAD0])); // if needed


}


// energy density loss rate integrated over frequency and solid angle
int calc_rad_lambda(FTYPE *pp, struct of_geom *ptrgeom, FTYPE kappa, FTYPE kappaes, FTYPE Tgas, FTYPE *lambda)
{

  // get gas properties
  FTYPE rho=pp[RHO];
  FTYPE u=pp[UU];
  //  FTYPE T=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho,u);


  // This is aT^4/(4\pi) that is the specific black body emission rate in B_\nu d\nu d\Omega corresponding to energy density rate per unit frequency per unit solid angle, which has been integrated over frequency.
  // More generally, kappa*4*Pi*B can be replaced by some \Lambda that is some energy density rate
  // But, have to be careful that "kappa rho" is constructed from \Lambda/(u*c) or else balance won't occur.
  // This is issue because "kappa" is often frequency integrated directly, giving different answer than frequency integrating j_v -> \Lambda/(4\pi) and B_\nu -> (aT^4)/(4\pi) each and then taking the ratio.
  // Note if T is near maximum for FTYPE, then aradT^4 likely too large.
  FTYPE B=0.25*ARAD_CODE*pow(Tgas,4.)/Pi;


  // energy density loss rate integrated over frequency and solid angle
  *lambda = kappa*4.*Pi*B;

  return(0);
}


// compute radiative characteristics as limited by opacity
int vchar_rad(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin, FTYPE *vmax2, FTYPE *vmin2,int *ignorecourant)
{

  
  // compute chi
  // Assume computed as d\tau/dorthonormallength as defined by user.
  // Assume \chi defined in fluid frame (i.e. not radiation frame).
  FTYPE kappa,chi;
  calc_chi(pr,geom,&chi);
  // KORALTODO: in paper, suggests only kappaes should matter?
  //  calc_kappa(pr,geom,&kappa);
  //  chi=kappa;


  //characterisitic wavespeed in the radiation rest frame
  FTYPE vrad2=THIRD;
  FTYPE vrad2limited;

  if(chi>0.0){// && WHICHRADSOURCEMETHOD==SOURCEMETHODIMPLICIT){
    // NOT DOING THIS:
    // compute tautot assuming chi is optical depth per unit grid dx[1-3].  I.e. calc_chi() computes grid-based opacity
    // tautot is the total optical depth of the cell in dim dimension
    //  tautot = chi * dx[dir];

    // DOING THIS:
    // KORALTODO: Approximation to any true path, but approximation is sufficient for approximate wave speeds.
    // \tau_{\rm tot}^2 \approx \chi^2 [dx^{dir} \sqrt{g_{dirdir}}]^2 
    FTYPE tautotsq,vrad2tau;
    // Note that tautot is frame independent once multiple \chi by the cell length.  I.e. it's a Lorentz invariant.
    tautotsq = chi*chi * dx[dir]*dx[dir]*fabs(geom->gcov[GIND(dir,dir)]);

    //    dualfprintf(fail_file,"chi=%g dx=%g dir=%d tautot=%g\n",chi,dx[dir],dir,sqrt(tautotsq));
  
    vrad2tau=(4.0/3.0)*(4.0/3.0)/tautotsq; // KORALTODO: Why 4.0/3.0 ?  Seems like it should be 2.0/3.0 according to NR1992 S19.2.6L or NR2007 S20.2L with D=1/(3\chi), but twice higher speed is only more robust.
    vrad2limited=MIN(vrad2,vrad2tau);

    // NOTEMARK: For explicit method, this will lead to very large dt relative to step desired by explicit method, leading to ever more sub-cycles for WHICHRADSOURCEMETHOD==SOURCEMETHODEXPLICITSUBCYCLE method.

    // TODOMARK: I wonder if another possibility is to use a speed limiter in the advection equation.  With my pseudo-Newtonian code is has a limiter on the sound and Alfven speeds following the idea of limiting the Alfven speed by Miller & Stone (2000, http://adsabs.harvard.edu/abs/2000ApJ...534..398M).  That is, there must be a way to insert a term into the radiation advection equations to limit the velocity to ~c/\tau that only becomes effective at and beyond that speed.  Then the Jacobian would be modified (Or thinking of how the Jacobian could get modified, one gets a different equation of motion).

  }
  else{
    vrad2limited=vrad2;
  }


  // for setting flux so diffusive term is not exaggerated in high \tau regions
  simplefast_rad(dir,geom,q,vrad2limited,vmin,vmax);

  
  if(FORCESOLVELFLUX){
    FTYPE ftemp=1.0/sqrt(fabs(geom->gcov[GIND(dir,dir)]));
    *vmin=-ftemp;
    *vmax=+ftemp;
  }
  //    cminmaxrad_l[CMIN]=-ftemp;
  //    cminmaxrad_l[CMAX]=+ftemp;
  //      cminmax_calc(cminmaxrad_l[CMIN],cminmaxrad_r[CMIN],cminmaxrad_l[CMAX],cminmaxrad_r[CMAX],&cminmaxrad[CMIN],&cminmaxrad[CMAX],ctopradptr);
  //      ctoprad=ftemp;

#if(1)
  *vmin2=*vmin;
  *vmax2=*vmax;
#elif(0)
  // for setting timestep since advective part has no knowledge of \tau-limited velocity
  simplefast_rad(dir,geom,q,vrad2,vmin2,vmax2);
#else
  // works even if not using damping of implicit solver
  simplefast_rad(dir,geom,q,2.0/3.0,vmin2,vmax2);
#endif

  
  return(0);
}

// get lab-frame 3-velocity for radiative emission in radiative frame
static int simplefast_rad(int dir, struct of_geom *geom,struct of_state *q, FTYPE vrad2,FTYPE *vmin, FTYPE *vmax)
{
  extern int simplefast(int dir, struct of_geom *geom,struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax);

  //need to substitute ucon,ucov with uradcon,uradcov to fool simplefast
  FTYPE ucon[NDIM],ucov[NDIM];
  int ii;
  DLOOPA(ii){
    ucon[ii]=q->ucon[ii];
    ucov[ii]=q->ucov[ii];
    q->ucon[ii]=q->uradcon[ii];
    q->ucov[ii]=q->uradcov[ii];
  }

  //calculating vmin, vmax
  simplefast(dir,geom,q,vrad2,vmin,vmax);

  //restoring gas 4-velocities
  DLOOPA(ii){
    q->ucon[ii]=ucon[ii];
    q->ucov[ii]=ucov[ii];
  }


#if(0)
  // Cartesian-Minkowski speed-of-light limit of radiation velocity
  FTYPE dxdxp[NDIM][NDIM];
  dxdxprim_ijk(geom->i, geom->j, geom->k, geom->p, dxdxp);
  // characeristic wavespeeds are 3-velocity in lab-frame
  *vmin=-1.0/dxdxp[dir][dir]; // e.g. dxdxp = dr/dx1
  *vmax=+1.0/dxdxp[dir][dir];
#endif


  return(0);
}


// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  void compute_1plusud0(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])

  // urad^\mu
  // ucon_calc() assumes primitive velocities are in U1 through U3, but otherwise calculation is identical for radiation velocity, so just shift effective list of primitives so ucon_calc() operates on U1RAD through U3RAD
  MYFUN(ucon_calc(&pr[URAD1-U1], ptrgeom, q->uradcon,q->othersrad) ,"phys.c:get_state()", "ucon_calc()", 1);
  // urad_\mu
  lower_vec(q->uradcon, ptrgeom, q->uradcov);


  return (0);
}


void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM])
{
  int jj,kk;
  
  if(EOMRADTYPE!=EOMRADNONE){
    DLOOPA(jj){
      mhd_calc_rad( pr, jj, ptrgeom, q, &(radstressdir[jj][0]) , NULL );
    }
  }
  else DLOOP(jj,kk) radstressdir[jj][kk]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.
}

// compute radiation stres-energy tensor assuming M1 closure
void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *ptrgeom, struct of_state *q, FTYPE *radstressdir, FTYPE *radstressdirabs)
{
  int jj;
  FTYPE term1[NDIM],term2[NDIM];

  // R^{dir}_{jj} radiation stress-energy tensor
  if(EOMRADTYPE==EOMRADEDD){
    // force radiation frame to be fluid frame
    DLOOPA(jj){
      term1[jj]=THIRD*pr[PRAD0]*(4.0*q->ucon[dir]*q->ucov[jj]);
      term2[jj]=THIRD*pr[PRAD0]*(delta(dir,jj));
      radstressdir[jj]=term1[jj]+term2[jj];
      if(radstressdirabs!=NULL) radstressdirabs[jj]=fabs(term1[jj])+fabs(term2[jj]);
    }
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){
    DLOOPA(jj) radstressdir[jj]=THIRD*pr[PRAD0]*(4.0*q->uradcon[dir]*q->uradcov[jj] + delta(dir,jj));

  }
  else{
    // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.
    DLOOPA(jj){
      radstressdir[jj]=0.0;
      if(radstressdirabs!=NULL) radstressdirabs[jj]=0.0;
    }
  }


}

// compute fluid frame orthonormal basis radiation stress-energy tensor assuming M1 closure
int calc_Rij_ff(FTYPE *pp, FTYPE Rij[][NDIM])
{
  FTYPE E=pp[PRAD0];
  FTYPE F[NDIM-1]={pp[PRAD1],pp[PRAD2],pp[PRAD3]};

  FTYPE nx,ny,nz,nlen,f;

  nx=F[0]/E;
  ny=F[1]/E;
  nz=F[2]/E;
  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
  if(EOMRADTYPE==EOMRADEDD){
    f=1./3.; // f and Rij are both as if nx=ny=nz=0
    //  f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){

    if(nlen>=1.) f=1.; // KORALTODO: limiter, but only used so far for IC
    else  f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  //M1
  }
  else if(EOMRADTYPE==EOMRADNONE){

  }
  else{
    dualfprintf(fail_file,"1 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
    myexit(837453242);
  }
  
  ////////// Get R^{ij} in orthonormal fluid frame 
  Rij[0][0]=E;

  if(EOMRADTYPE==EOMRADEDD){
    // KORALTODO: Below 3 should be zero for Eddington approximation, but only if F=0 exactly.
    Rij[0][1]=Rij[1][0]=0.0;
    Rij[0][2]=Rij[2][0]=0.0;
    Rij[0][3]=Rij[3][0]=0.0;
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){
    Rij[0][1]=Rij[1][0]=F[0];
    Rij[0][2]=Rij[2][0]=F[1];
    Rij[0][3]=Rij[3][0]=F[2];
  }
  else if(EOMRADTYPE==EOMRADNONE){

  }
  else{
    dualfprintf(fail_file,"2 No Such EOMRADTYPE=%d\n",EOMRADTYPE);
    myexit(837453243);
  }


  // normalize n^i for Rij calculation
  if(nlen>0){
    nx/=nlen;
    ny/=nlen;
    nz/=nlen;
  }
  else{
    ;
  }

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}



FTYPE my_min(FTYPE a, FTYPE b)
{
  if(a<b) return a;
  else return b;
}

FTYPE my_sign(FTYPE x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 0.;
  return 0;
}




//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 4by4 matrix
// gives inverse transpose matrix
int inverse_44matrix(FTYPE a[][NDIM], FTYPE ia[][NDIM])
{
  FTYPE mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  FTYPE tmp[12]; FTYPE src[16]; FTYPE det, idet;
  /* transpose matrix */
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

 
  /* calculate matrix inverse */
  idet = 1.0/det;

  //  if(isnan(idet)){
  if(!isfinite(idet) || !isfinite(det)){
    dualfprintf(fail_file,"idet (det=%26.15g idet=%26.15g) in inverse 4x4 zero or nan\n",det,idet);
    return(1); // indicates failure
    //    myexit(13235);
  }

  for (j = 0; j < 16; j++)
    dst[j] *= idet;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      ia[i][j]= dst[i*4+j];

  return 0;
}



//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 3by3 matrix
// gives inverse transpose matrix
static int inverse_33matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM])
{

  FTYPE det = +a[sj+0][sj+0]*(a[sj+1][sj+1]*a[sj+2][sj+2]-a[sj+2][sj+1]*a[sj+1][sj+2])
    -a[sj+0][sj+1]*(a[sj+1][sj+0]*a[sj+2][sj+2]-a[sj+1][sj+2]*a[sj+2][sj+0])
    +a[sj+0][sj+2]*(a[sj+1][sj+0]*a[sj+2][sj+1]-a[sj+1][sj+1]*a[sj+2][sj+0]);
  FTYPE idet = 1.0/det;
  ia[sj+0][sj+0] =  (a[sj+1][sj+1]*a[sj+2][sj+2]-a[sj+2][sj+1]*a[sj+1][sj+2])*idet;
  ia[sj+1][sj+0] = -(a[sj+0][sj+1]*a[sj+2][sj+2]-a[sj+0][sj+2]*a[sj+2][sj+1])*idet;
  ia[sj+2][sj+0] =  (a[sj+0][sj+1]*a[sj+1][sj+2]-a[sj+0][sj+2]*a[sj+1][sj+1])*idet;
  ia[sj+0][sj+1] = -(a[sj+1][sj+0]*a[sj+2][sj+2]-a[sj+1][sj+2]*a[sj+2][sj+0])*idet;
  ia[sj+1][sj+1] =  (a[sj+0][sj+0]*a[sj+2][sj+2]-a[sj+0][sj+2]*a[sj+2][sj+0])*idet;
  ia[sj+2][sj+1] = -(a[sj+0][sj+0]*a[sj+1][sj+2]-a[sj+1][sj+0]*a[sj+0][sj+2])*idet;
  ia[sj+0][sj+2] =  (a[sj+1][sj+0]*a[sj+2][sj+1]-a[sj+2][sj+0]*a[sj+1][sj+1])*idet;
  ia[sj+1][sj+2] = -(a[sj+0][sj+0]*a[sj+2][sj+1]-a[sj+2][sj+0]*a[sj+0][sj+1])*idet;
  ia[sj+2][sj+2] =  (a[sj+0][sj+0]*a[sj+1][sj+1]-a[sj+1][sj+0]*a[sj+0][sj+1])*idet;

  if(!isfinite(det) || !isfinite(idet)){
    dualfprintf(fail_file,"inverse_33matrix got singular det=%g idet=%g\n",det,idet);
    return(1); // indicates failure
    //    myexit(13235);
  }

  return(0);
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 1by1 matrix
// gives inverse transpose matrix (for 1by1, transpose does nothing)
static int inverse_11matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM])
{
  // trivial inversion, and can't fail unless divide by zero
  // sj==endjac

  ia[sj][sj]=1.0/a[sj][sj];

  if(!isfinite(ia[sj][sj])){
    dualfprintf(fail_file,"inverse 1x1 zero or nan\n",a[sj][sj]);
    return(1); // indicates failure
    //    myexit(13235);
  }
  return(0);

}





/*********************************************************************************/
/****** radiative ortonormal ff primitives (E,F^i) <-> primitives in lab frame  *******/
// Used only for initial conditions
/*********************************************************************************/
// whichvel: input vel type for U1-U3
// whichcoord: input coord type for both U1-U3 and URAD1-URAD3
// whichdir: LAB2FF or FF2LAB  . In addition, here lab means HARM-lab different by alpha factor from true lab.
// i,j,k,loc = standard grid location
// ptrgeom: any input geometry for the lab frame (ptrgeom could be from MCOORD, PRIMECOORDS, etc.) (same for pin's velocity as well as orthonormal basis)
//          If ptrgeom==NULL, then use i,j,k,loc to get geometry in whichcoord coordinates
// pradffortho: radiation primitives (PRAD0-3) should be fluid-frame orthonormal basis values (i.e. E,F in fluid frame orthonormal basis)
// pin: inputs for primitives (i.e. whichvel for U1-U3 and whichcoord for U1-U3,URAD1-URAD3)
// pout: outputs for primitives ("")
int prad_fforlab(int *whichvel, int *whichcoord, int whichdir, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{

  if(whichdir==FF2LAB) prad_fftolab(whichvel, whichcoord, i, j, k, loc, ptrgeom, pradffortho, pin, pout);
  else if(whichdir==LAB2FF) prad_labtoff(whichvel, whichcoord, i, j, k, loc, ptrgeom, pradffortho, pin, pout);
  else{
    dualfprintf(fail_file,"prad_fforlab() not yet setup for whichdir=%d.",whichdir);
    myexit(652526624);
  }

  return(0);

}

// like prad_fforlab() but for only whichdir=LAB2FF
// used for dumps or diags
int prad_labtoff(int *whichvel, int *whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{
  int jj;

  //  DLOOPA(jj) dualfprintf(fail_file,"ijk=%d %d %d : jj=%d pin=%g\n",i,j,k,jj,pin[PRAD0+jj]);

  // assume ptrgeom is PRIMECOORDS lab-frame geometry
  struct of_geom geomdontuse;
  if(ptrgeom==NULL){
    ptrgeom=&geomdontuse;
    get_geometry(i, j, k, loc, ptrgeom);
  }


  
  // get state
  struct of_state q;
  get_state(pin,ptrgeom,&q);

  //  DLOOPA(jj) dualfprintf(fail_file,"ijk=%d %d %d : jj=%d uradcon=%g uradcov=%g\n",i,j,k,jj,q.uradcon[jj],q.uradcov[jj]);


  // get lab-frame R^\mu_\nu
  FTYPE Rijlab[NDIM][NDIM];
  mhdfull_calc_rad(pin, ptrgeom, &q, Rijlab);

#if(0) // STAY AS ZERO

  // get U=R^t_\mu [harm type]
  FTYPE U[NDIM];
  DLOOPA(jj) U[jj]=Rijlab[TT][jj];

  //  DLOOPA(jj) dualfprintf(fail_file,"ijk=%d %d %d : jj=%d U=%g ucon=%g ucov=%g\n",i,j,k,jj,U[jj],q.ucon[jj],q.ucov[jj]);

  // transform lab-frame R^t_\nu [harm type] to fluid-frame version
  FTYPE Uff[NDIM];
  vector_harm2orthofluidorback(TYPEUCOV, LAB2FF, ptrgeom, TYPEUCON, q.ucon, TYPEUCOV, U, Uff);

  //  DLOOPA(jj) dualfprintf(fail_file,"ijk=%d %d %d : jj=%d Uff=%g alpha=%g\n",i,j,k,jj,Uff[jj],ptrgeom->alphalapse);

  //                            00 01 02 03 11 12 13 22 23 33
  //  FTYPE etamink[SYMMATRIXNDIM]={-1 ,0 ,0 ,0, 1 ,0 ,0 ,1 ,0 ,1};
  Uff[TT]*=1.0; // if original was R_\mu get R^\mu in fluid frame orthonormal basis

  DLOOPA(jj) pradffortho[PRAD0+jj] = Uff[jj];

#else

  //  int kk;
  //  DLOOP(jj,kk) dualfprintf(fail_file,"gn%d%d=%21.15g\n",jj+1,kk+1,ptrgeom->gcon[GIND(jj,kk)]);
  //  DLOOP(jj,kk) dualfprintf(fail_file,"gv%d%d=%21.15g\n",jj+1,kk+1,ptrgeom->gcov[GIND(jj,kk)]);

  indices_2122(Rijlab,Rijlab,ptrgeom);

  // Need to use full Rijlab since can be mixing between components in general
  // transform and boost (ultimately converts pin -> Rijlab-> pradffortho -> Rijff effectively)
  int tconcovtypeA=TYPEUCON;
  int tconcovtypeB=TYPEUCON;
  int primcoord=1;// 1 so that will use optimal way to get tetrads
  FTYPE Rijff[NDIM][NDIM];
  tensor_lab2orthofluidorback(primcoord, LAB2FF, ptrgeom, TYPEUCON, q.ucon, tconcovtypeA, tconcovtypeB, Rijlab, Rijff);

  //  DLOOPA(jj) dualfprintf(fail_file,"ijk=%d %d %d : jj=%d Rijff[TT]=%g alpha=%g\n",i,j,k,jj,Rijff[TT][jj],ptrgeom->alphalapse);

  // get in pradffortho form
  DLOOPA(jj) pradffortho[PRAD0+jj] = Rijff[TT][jj];
  
#endif

  // just copy pout
  int pliter,pl;
  PLOOP(pliter,pl) pout[pl]=pin[pl];

  return(0);
}

// like prad_fforlab() but for only whichdir=FF2LAB
// used for IC
int prad_fftolab(int *whichvel, int *whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{
  FTYPE Rijff[NDIM][NDIM],Rijlab[NDIM][NDIM],U[NPR]={0};
  int pliter,pl;
  int primcoord;
  int jj,kk;
  struct of_geom geomtousedontuse;
  struct of_geom *ptrgeomtouse=&geomtousedontuse;


  if(ptrgeom==NULL){
    if(*whichcoord!=PRIMECOORDS){
      // get metric grid geometry for these ICs
      int getprim=0;
      gset_genloc(getprim,*whichcoord,i,j,k,loc,ptrgeomtouse);
    }
    else{
      get_geometry(i, j, k, loc, ptrgeomtouse);
    }
  }
  else{
    // then assumes ptrgeom is in *whichcoord coordinates
    ptrgeomtouse=ptrgeom;
  }


  // set primitive that can use as pre-existing fluid velocity if need to use for reduction
  // also use pout instead of pin so preserves pin no matter what (unless user set pin=pout)
  PLOOP(pliter,pl) pout[pl]=pin[pl];


  // radiative stress tensor in the fluid frame orthonormal basis
  // assuming input pradffortho for radiation is in fluid frame orthonormal basis, but in "primitive" format so using pradffortho[PRAD0-PRAD3]
  // gets R^{ij} in fluid frame orthonormal basis from primitive quantities in fluid frame orthonormal basis
  calc_Rij_ff(pradffortho,Rijff);
  
  //  PLOOPRADONLY(pl) dualfprintf(fail_file,"pl=%d pout=%g\n",pl,pout[pl]);
  //  DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijff=%g\n",jj,kk,Rijff[jj][kk]);
  //  DLOOP(jj,kk) dualfprintf(fail_file,"gn%d%d=%21.15g\n",jj+1,kk+1,ptrgeomtouse->gcon[GIND(jj,kk)]);
  //  DLOOP(jj,kk) dualfprintf(fail_file,"gv%d%d=%21.15g\n",jj+1,kk+1,ptrgeomtouse->gcov[GIND(jj,kk)]);


  // get ucon (assumed primitive velocity in ptrgeomtouse coordinates)
  FTYPE ucon[NDIM],others[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(*whichvel,pout,ptrgeomtouse,ucon,others);

  //  DLOOPA(jj) dualfprintf(fail_file,"jj=%d ucon=%g\n",jj,ucon[jj]);


  // also convert whichvel ucon to VELREL4 primitive velocity for use by u2p_rad() and as needed for consistent final output from this function and as possible backup value
  if(*whichvel!=VELREL4) ucon2pr(VELREL4,ucon,ptrgeomtouse,pout);

  //  SLOOPA(jj) dualfprintf(fail_file,"jj=%d u4rel=%g\n",jj,pout[UU+jj]);
  
  // transform and boost (ultimately converts pradffortho -> Rijff -> Rijlab -> U)
  int tconcovtypeA=TYPEUCON;
  int tconcovtypeB=TYPEUCON;
  if(*whichcoord==PRIMECOORDS) primcoord=1;
  else primcoord=0;
  tensor_lab2orthofluidorback(primcoord, FF2LAB, ptrgeomtouse, TYPEUCON, ucon, tconcovtypeA, tconcovtypeB, Rijff, Rijlab);

  //  DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijlab=%g\n",jj,kk,Rijlab[jj][kk]);

  //R^munu -> R^mu_nu so in standard form to extract conserved quantity R^t_\nu
  indices_2221(Rijlab,Rijlab,ptrgeomtouse);

  //  DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Ridownjlab=%g\n",jj,kk,Rijlab[jj][kk]);

  // Store radiation conserved quantity from R^t_\nu .  u2p_rad() below only uses radiation U's.
#if(0) // STAY ZERO NOW
  // for true lab to fake-harm lab, end up dividing by alpha (see vector_harm2orthofluidorback() in tetrad.c)
  FTYPE alpha=ptrgeomtouse->alphalapse;
  U[URAD0]=Rijlab[TT][TT]/alpha;
  U[URAD1]=Rijlab[TT][RR]/alpha;
  U[URAD2]=Rijlab[TT][TH]/alpha;
  U[URAD3]=Rijlab[TT][PH]/alpha;
#else
  U[URAD0]=Rijlab[TT][TT];
  U[URAD1]=Rijlab[TT][RR];
  U[URAD2]=Rijlab[TT][TH];
  U[URAD3]=Rijlab[TT][PH];
#endif

  //  DLOOPA(jj) dualfprintf(fail_file,"jj=%d URAD=%g\n",jj,U[URAD0+jj]);




  PFTYPE lpflag=UTOPRIMNOFAIL,lpflagrad=UTOPRIMRADNOFAIL;
  int showmessages=1; // LEAVE on (not normal debugging)
  int allowlocalfailurefixandnoreport=1;
  // NOTEMARK: lpflag=UTOPRIMNOFAIL means accept input pout for velocity to maybe be used in local reductions to fluid frame.
  // u2p_rad() only uses U[URAD0-URAD3]
  // generally u2p_rad() could use all of pout[] except only assigns pout[PRAD0-PRAD3] and doesn't use that for anything except as "static" solution (i.e. uses pin effectively)
  u2p_rad(showmessages, allowlocalfailurefixandnoreport, GAMMAMAXRAD, CAPTYPEFIX1, U, pout, ptrgeomtouse, &lpflag, &lpflagrad);

  //  DLOOPA(jj) dualfprintf(fail_file,"u2p_rad: jj=%d pout=%g\n",jj,pout[PRAD0+jj]);



  // get back to whichvel
  FTYPE uconback[NDIM],othersback[NUMOTHERSTATERESULTS];
  // for fluid
  ucon_calc_whichvel(VELREL4,pout,ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,pout);
  // KORALTODO: for radiation (always returned as VELREL4 so far.
  ucon_calc_whichvel(VELREL4,&pout[URAD1-U1],ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,&pout[URAD1-U1]);



  // DEBUG:
  if(lpflag!=UTOPRIMNOFAIL || lpflagrad!=UTOPRIMRADNOFAIL){ // DEBUG with 1||
    dualfprintf(fail_file,"Failed to invert during prad_fftolab().  Assuming fixups won't be applied: %d %d\n",lpflag,lpflagrad);
    dualfprintf(fail_file,"ijk=%d %d %d : %d\n",ptrgeomtouse->i,ptrgeomtouse->j,ptrgeomtouse->k,ptrgeomtouse->p);
    PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pin=%g U=%g\n",pl,pin[pl],U[pl]);
    DLOOPA(jj) dualfprintf(fail_file,"jj=%d ucon=%g\n",jj,ucon[jj]);
    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijff=%g Rijlab=%g\n",jj,kk,Rijff[jj][kk],Rijlab[jj][kk]);
    DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d gcov=%g gcon=%g\n",jj,kk,ptrgeomtouse->gcov[GIND(jj,kk)],ptrgeomtouse->gcon[GIND(jj,kk)]);
    PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pout=%g\n",pl,pout[pl]);
    myexit(189235);
    // KORALTODO: Check whether really succeeded?  Need to call fixups?  Probably, but need per-cell fixup.  Hard to do if other cells not even set yet as in ICs.  Should probably include fixup process during initbase.c stuff.
  }



  return 0;
} 












// for BCs, to take E[radiation frame] and u^i as radiation primitives in whichvel/whichcoord
// obtains WHICHVEL/PRIMECOORD primitives
int primefluid_EVrad_to_primeall(int *whichvel, int *whichcoord, struct of_geom *ptrgeom, FTYPE *pin, FTYPE *pout)
{
  int pliter,pl;
  int i=ptrgeom->i;
  int j=ptrgeom->j;
  int k=ptrgeom->k;
  int loc=ptrgeom->p;

  // copy over
  PLOOP(pliter,pl) pout[pl]=pin[pl];

  // get metric grid geometry for these ICs
  int getprim=0;
  struct of_geom geomrealdontuse;
  struct of_geom *ptrgeomreal=&geomrealdontuse;
  gset_genloc(getprim,*whichcoord,i,j,k,loc,ptrgeomreal);

  FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(*whichvel,&pout[URAD1-U1],ptrgeomreal,uradcon,othersrad);

  // now convert velocity so in PRIMECOORDS assuming whichcoord=MCOORD
  mettometp_genloc(i,j,k,loc,uradcon);

  if(*whichcoord!=MCOORD){
    dualfprintf(fail_file,"primefluid_EVrad_to_primeall() needs whichcoord (%d) to be MCOORD (%d)\n",whichcoord,MCOORD);
    myexit(87345246);
  }

  // assumed already inputted PRIMECOORDS WHICHVEL for fluid velocity, so no conversion for the fluid velocity

  // now go from ucon[PRIMECOORDS] -> primitive[PRIMECOORDS] for radiation velocity and get WHICHVEL version
  ucon2pr(WHICHVEL,uradcon,ptrgeom,&pout[URAD1-U1]);

  // now all PRIMECOORDS WHICHVEL type assuming ptrgeom inputted PRIMECOORDS version as expected
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  
  return(0);
}


// Input: start with pin [with fluid in whichvel velocity and whichcoordfluid coordinates (PRIMECOORDS or MCOORD) and radiation as E,F in fluid frame orthonormal basis in whichcoordrad coordinates]
// Output: pout [with all WHICHVEL PRIMECOORDS and radiation using velocity primitive]
//
// Useful for BCs when have (say) VEL3,MCOORD for fluid velocity as well as E,F in ff for radiation and need normal WHICHVEL PRIMECOORDS fluid velocity as well as normal velocity primitive for radiation
int whichfluid_ffrad_to_primeall(int *whichvel, int *whichcoordfluid, int *whichcoordrad, struct of_geom *ptrgeomprimecoords, FTYPE *pradffortho, FTYPE *pin, FTYPE *pout)
{
  int pliter,pl;
  int i=ptrgeomprimecoords->i;
  int j=ptrgeomprimecoords->j;
  int k=ptrgeomprimecoords->k;
  int loc=ptrgeomprimecoords->p;


  //  PLOOP(pliter,pl) dualfprintf(fail_file,"ijk=%d %d %d pl=%d pin0=%g pout0=%g\n",i,j,k,pl,pin[pl],pout[pl]);

  // prad_fforlab() should only use radiation primitives, but copy all primitives so can form ucon for transformation
  PLOOP(pliter,pl) pout[pl]=pin[pl];

  // 4 cases:
  // rad   fluid
  // PRIME PRIME
  // PRIME other
  // other PRIME
  // other other

  // get real geometry if needed
  struct of_geom geomrealraddontuse;
  struct of_geom *ptrgeomrealrad=&geomrealraddontuse;
  struct of_geom geomrealfluiddontuse;
  struct of_geom *ptrgeomrealfluid=&geomrealfluiddontuse;
  if(*whichcoordrad!=PRIMECOORDS){
    int getprim=0;
    gset_genloc(getprim,*whichcoordrad,i,j,k,loc,ptrgeomrealrad);
  }
  else ptrgeomrealrad=ptrgeomprimecoords;
  if(*whichcoordfluid!=PRIMECOORDS){
    int getprim=0;
    gset_genloc(getprim,*whichcoordfluid,i,j,k,loc,ptrgeomrealfluid);
  }
  else ptrgeomrealfluid=ptrgeomprimecoords;



  // make whichcoord for fluid same as for rad before continuing (make fluid same as radiation by only changing fluid whichcoord)
  if(*whichcoordfluid!=*whichcoordrad){
    FTYPE ucon[NDIM];
    FTYPE others[NUMOTHERSTATERESULTS];

    // pr->ucon for fluid
    if (pr2ucon(*whichvel,pout, ptrgeomrealfluid, ucon) >= 1) FAILSTATEMENT("bounds.koral.c:bl2met2metp2v_genloc() for radiation", "pr2ucon()", 2);

    if(*whichcoordfluid==PRIMECOORDS){ // then radiation is not PRIMECOORDS, so go to radiation coords
      metptomet_genloc(i,j,k,loc,ucon); // now ucon is in MCOORD
      // convert MCOORD->whichcoordrad
      coordtrans(MCOORD,*whichcoordrad,i,j,k,loc,ucon);
    }
    else if(*whichcoordrad==PRIMECOORDS){ // then go to PRIMECOORDS for fluid
      // convert whichcoordfluid->MCOORD
      coordtrans(*whichcoordfluid,MCOORD,i,j,k,loc,ucon);
      mettometp_genloc(i,j,k,loc,ucon); // now ucon is in PRIMECOORDS
    }
    else{ // then neither is PRIMECOORDS, so just transform and skip mettometp() or metptomet()
      // convert whichcoordfluid->whichcoordrad
      coordtrans(*whichcoordfluid,*whichcoordrad,i,j,k,loc,ucon);
    }

    // get whichvel primitive
    ucon2pr(*whichvel,ucon,ptrgeomrealrad,pout);

    // changed fluid to have same whichcoord as radiation (set by radiation), so set that
    *whichcoordfluid=*whichcoordrad;

  }
  else{
    // otherwise, whichcoord same for fluid and radiation, so can continue
  }


  //  PLOOP(pliter,pl) dualfprintf(fail_file,"PRE:  ijk=%d %d %d pl=%d pout=%g\n",i,j,k,pl,pout[pl]);



  // get WHICHVEL primitives (still will be in whichcoord coordinates)
  // whichvel here is for fluid velocity (prad_fforlab() converts velocity to WHICHVEL for consistency with only currently allowed output of radiation velocity)
  // pradffortho assumed as in orthonormal fluid frame, but coordinates of whichcoordrad
  int whichframedir=FF2LAB; // fluid frame orthonormal to lab-frame
  prad_fforlab(whichvel, whichcoordrad, whichframedir, i, j, k, loc, ptrgeomrealrad, pradffortho, pout, pout);

  // output from prad_fforlab() is always WHICHVEL for both fluid and radiation primitives
  // changed whichvel's, so report that back if needed
  //  *whichvel=WHICHVEL;
  // above change of whichvel no longer true (and anyways, whichvel was changed in prad_fforlab() directly)
 
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"POST: ijk=%d %d %d pl=%d pout=%g\n",i,j,k,pl,pout[pl]);

 
  // output from prad_fforlab() not yet necessarily PRIMECOORDS.
  if(*whichcoordrad==MCOORD){
    // Get all primitives in WHICHVEL/PRIMECOORDS (no conversion for WHICHVEL since prad_fforlab() already put quantities in WHICHVEL due to u2p_rad() only setup for WHICHVEL)
    if (bl2met2metp2v_genloc(*whichvel, *whichcoordrad, pout, i,j,k,loc) >= 1){
      FAILSTATEMENT("bounds.koral.c:bound_radatmbeaminflow()", "bl2ks2ksp2v()", 1);
    }
  }

  // changed coordinates to PRIMECOORDS, so set that as the case
  *whichcoordfluid=*whichcoordrad=PRIMECOORDS;

  return(0);

}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// T^ij -> T^i_j
int indices_2221(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
        {
          Tt[i][j]=0.;
          for(k=0;k<NDIM;k++)
            {
              Tt[i][j]+=T1[i][k]*ptrgeom->gcov[GIND(k,j)];
            }   
        }
    }

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
        {
          T2[i][j]=Tt[i][j];
        }
    }

  return 0;
}

// T^i_j -> T^{ij}
int indices_2122(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
        {
          Tt[i][j]=0.;
          for(k=0;k<NDIM;k++)
            {
              Tt[i][j]+=T1[i][k]*ptrgeom->gcon[GIND(k,j)];
            }   
        }
    }

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
        {
          T2[i][j]=Tt[i][j];
        }
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A^i -> A^_j
int indices_21(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
        {
          At[i]+=A1[k]*ptrgeom->gcov[GIND(i,k)];
        }   
    }

  for(i=0;i<NDIM;i++)
    {
      A2[i]=At[i];
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A_i -> A^_j
int indices_12(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
        {
          At[i]+=A1[k]*ptrgeom->gcon[GIND(i,k)];
        }   
    }

  for(i=0;i<NDIM;i++)
    {
      A2[i]=At[i];
    }

  return 0;
}



int u2p_rad(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, int whichcap, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  int u2p_rad_orig(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad);
  int toreturn;

  toreturn=u2p_rad_orig(showmessages, allowlocalfailurefixandnoreport, gammamaxrad, uu, pin, ptrgeom,lpflag, lpflagrad);

  return(toreturn);
}




//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//**********************************************************************
//**********************************************************************
//
///////////////
//
// Invert U->direct Primitive for radiation
// OLD (i.e. no longer true): (must come after HD or MHD or whatever sets velocity of fluid, because radiation needs to have updated velocity so that can define fluid frame)
// old code inside utoprimgen.c was:
//    struct of_state qrad;
// this uses new pr to get only ucon and ucov
//get_state_uconucovonly(pr, ptrgeom, &qrad); // OLD
// get new radiation primitives
//
// NEW (currently true): fluid frame no longer needed because go directly from lab-frame conserved quantities to lab-frame primitive quantities.
//
//
// uu: Conserved quantities with URAD0,1,2,3 as radiation conserved quantities
// pp: primitives with PRAD0,1,2,3 as radiation primitive quantities
// ptrgeom: Standard pointer to geometry
// lpflag: see gobal.nondepmnemonics.h .  Tells u2p_rad() if can use/trust fluid velocity.
// lpflagrad: Should be set to indicate success of u2p_rad() inversion
//
// NOTES:
//
// Using *lpflag<=UTOPRIMNOFAIL to check for fluid inversion success rather than a SOFTer condition (e.g. no fail or IFUTOPRIMFAILSOFT==1) because only want to trust fluid as reduction of M1 in case where velocity is accurate with non-negative densities.

// 0 or 1
// generally, should have TRYCOLD=1 as most general way to deal with failure
#define TRYCOLD 1

// for debugging
FTYPE globaluu[NPR];
FTYPE globalpin[NPR];

//
///////////////
int u2p_rad_orig(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  int jj,kk;
  FTYPE pp[NPR];
  int pliter,pl;

  PLOOP(pliter,pl) globaluu[pl]=uu[pl];
  PLOOP(pliter,pl) globalpin[pl]=pin[pl];

  if(WHICHVEL!=VELREL4){
    dualfprintf(fail_file,"u2p_rad() only setup for relative 4-velocity, currently.\n");
    myexit(137432636);
  }


  // copy over pin so pin isn't modified until end
  PLOOP(pliter,pl) pp[pl]=pin[pl];

  //////////////////////
  //
  // Prepare inversion from U->p for radiation assuming M1 closure
  //
  //////////////////////

  *lpflagrad=UTOPRIMRADNOFAIL;


  //conserved - R^t_mu
  FTYPE Avcov[NDIM]={uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3]};
  //indices up - R^tmu
  FTYPE Avcon[NDIM];
  indices_12(Avcov,Avcon,ptrgeom);


  FTYPE gammarel2,delta,numerator,divisor;
  FTYPE Erf;
  FTYPE urfconrel[NDIM];


  if(EOMRADTYPE==EOMRADEDD){
    // NOTEMARK: Can't use normal inversion that assumes R^t_i are independently evolved because they will generally lead to different velocity than fluid.

    // radiation is same as fluid gamma (assume fluid has already been inverted)
    urfconrel[1]=pp[PRAD1]=pp[U1];
    urfconrel[2]=pp[PRAD2]=pp[U2];
    urfconrel[3]=pp[PRAD3]=pp[U3];
 
    // get gammarel2
    FTYPE gammarel,qsq;
    gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammarel,&qsq);
    gammarel2=gammarel*gammarel;
 
    FTYPE alpha=ptrgeom->alphalapse; //sqrt(-1./ptrgeom->gcon[GIND(0,0)]);
    // get energy density in fluid frame from lab-frame
    Erf=3.*Avcon[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){

    // get \gamma^2 for relative 4-velocity
    get_m1closure_gammarel2(showmessages,ptrgeom,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor);

    if(0){
      // testing
      FTYPE Avconnew[NDIM]={Avcon[0],Avcon[1],Avcon[2],Avcon[3]};
      FTYPE urfconrelnew[NDIM];
      FTYPE gammarel2new,deltanew,numeratornew,divisornew,Erfnew;
      //      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avconnew,&gammarel2new,&deltanew,&numeratornew,&divisornew,&Erfnew,urfconrelnew);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avconnew,Avcov,NULL,&deltanew,&numeratornew,&divisornew,&Erfnew,urfconrelnew);
    }



    // get E in radiation frame
    get_m1closure_Erf(ptrgeom,Avcon,gammarel2,&Erf);
    FTYPE Erforig=Erf;

    // get relative 4-velocity
    if(CASECHOICE==JONCHOICE) get_m1closure_urfconrel(showmessages,allowlocalfailurefixandnoreport,ptrgeom,pp,Avcon,Avcov,gammarel2,delta,numerator,divisor,&Erf,urfconrel,lpflag,lpflagrad);
    else if(CASECHOICE==OLEKCHOICE) get_m1closure_urfconrel_olek(showmessages,allowlocalfailurefixandnoreport,ptrgeom,pp,Avcon,Avcov,gammarel2,delta,&Erf,urfconrel,lpflag,lpflagrad);

#if(0)
    // TESTING:
    FTYPE Erf2=Erforig,urfconrel2[NDIM];
    get_m1closure_urfconrel_olek(showmessages,allowlocalfailurefixandnoreport,ptrgeom,pp,Avcon,Avcov,gammarel2,delta,&Erf2,urfconrel2,lpflag,lpflagrad);
    FTYPE ERRORCHECK;
    ERRORCHECK=1E-1;
    if( fabs(Erf2-Erf)/(fabs(Erf2)+fabs(Erf))>ERRORCHECK || fabs(urfconrel2[1]-urfconrel[1])/(fabs(urfconrel2[1])+fabs(urfconrel[1]))>ERRORCHECK || fabs(urfconrel2[2]-urfconrel[2])/(fabs(urfconrel2[2])+fabs(urfconrel[2]))>ERRORCHECK || fabs(urfconrel2[3]-urfconrel[3])/(fabs(urfconrel2[3])+fabs(urfconrel[3]))>ERRORCHECK){
      dualfprintf(fail_file,"JONVSOLEK: ijk=%d %d %d : nstep=%ld steppart=%d : %g %g %g %g : %g %g %g %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,Erf,urfconrel[1],urfconrel[2],urfconrel[3],Erf2,urfconrel2[1],urfconrel2[2],urfconrel2[3]);
    }
    ERRORCHECK=0.4;
    if( fabs(Erf2-Erf)/(fabs(Erf2)+fabs(Erf))>ERRORCHECK || fabs(urfconrel2[1]-urfconrel[1])/(fabs(urfconrel2[1])+fabs(urfconrel[1]))>ERRORCHECK || fabs(urfconrel2[2]-urfconrel[2])/(fabs(urfconrel2[2])+fabs(urfconrel[2]))>ERRORCHECK || fabs(urfconrel2[3]-urfconrel[3])/(fabs(urfconrel2[3])+fabs(urfconrel[3]))>ERRORCHECK){
      dualfprintf(fail_file,"JONVSOLEK: ijk=%d %d %d : nstep=%ld steppart=%d : %g %g %g %g : %g %g %g %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,Erf,urfconrel[1],urfconrel[2],urfconrel[3],Erf2,urfconrel2[1],urfconrel2[2],urfconrel2[3]);
    }
#endif

  }// end if M1
  else{
    dualfprintf(fail_file,"No such EOMRADTYPE=%d in u2p_rad()\n",EOMRADTYPE);
    myexit(368322162);
  }


  //new primitives (only uses urfcon[1-3])
  pin[PRAD0]=Erf;
  pin[PRAD1]=urfconrel[1];
  pin[PRAD2]=urfconrel[2];
  pin[PRAD3]=urfconrel[3];

  //  DLOOPA(jj){
  //    if(!isfinite(pin[PRAD0+jj])){
  //      dualfprintf(fail_file,"caughtnan: jj=%d : ijk=%d %d %d\n",jj,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  //    }
  //  }

  if(DORADFIXUPS==1 || allowlocalfailurefixandnoreport==0){
    // KORALTODO: Problem is fixups can average across shock or place where (e.g.) velocity changes alot, and averaging diffuses shock and can leak-out more failures.
  }
  else{
    // CASE reductions (so set as no failure so fixups don't operate -- but might also want to turn off CHECKINVERSIONRAD else that routine won't know when to ignore bad U->P->U cases.)
    *lpflagrad=UTOPRIMRADNOFAIL;
  }

  return 0;
}





// interpolate between optically thick and thin limits when no u2p_rad() inversion solution
static int opacity_interpolated_urfconrel(FTYPE tautotmax, FTYPE *pp,struct of_geom *ptrgeom,FTYPE *Avcon, FTYPE Erf,FTYPE gammarel2,  FTYPE *Erfnew, FTYPE *urfconrel)
{
  int jj;
  FTYPE alpha=ptrgeom->alphalapse; //sqrt(-1./ptrgeom->gcon[GIND(0,0)]);

  //  dualfprintf(fail_file,"Erf=%g gammarel2=%g\n",Erf,gammarel2);

  FTYPE gammafluid,gammarel2fluid,qsqfluid,Erffluid;
  gamma_calc_fromuconrel(&pp[U1-1],ptrgeom,&gammafluid,&qsqfluid);
  gammarel2fluid=gammafluid*gammafluid;
  get_m1closure_Erf(ptrgeom, Avcon, gammarel2fluid, &Erffluid);
  if(Erffluid<ERADLIMIT) Erffluid=ERADLIMIT;

  FTYPE gammarad,gammarel2rad,qsqrad,Erfrad;
  gamma_calc_fromuconrel(&pp[URAD1-1],ptrgeom,&gammarad,&qsqrad);
  gammarel2rad=gammarad*gammarad;
  get_m1closure_Erf(ptrgeom, Avcon, gammarel2rad, &Erfrad);
  if(Erfrad<ERADLIMIT) Erfrad=ERADLIMIT;

  // now set urfconrel.  Choose fluid if tautotmax>=2/3 (updated fluid value), while choose previous radiation value (i.e. static!)
  // limit for interpolation below
  // below makes no sense because even for tau<<1 Erf and uradcon^i can be very different, so can end-up using too much of fluid version
  //  FTYPE tautotmaxlim=MIN(fabs(tautotmax),1.0);
  FTYPE tautotmaxlim=(1.0 - 1.0/(1.0+fabs(tautotmax)));
  // done with Erf, so get Erfnew (Erf and *Erfnew might be same variable, but Erf passed by value so changing *Erfnew won't change Erf anyways)
  *Erfnew = (1.0-tautotmaxlim)*Erfrad + tautotmaxlim*Erffluid;
  SLOOPA(jj) urfconrel[jj] = (1.0-tautotmaxlim)*pp[URAD1+jj-1] + tautotmaxlim*pp[U1+jj-1];

  dualfprintf(fail_file,"i=%d tautotmax=%g tautotmaxlim=%g\n",ptrgeom->i,tautotmax,tautotmaxlim);
  SLOOPA(jj) dualfprintf(fail_file,"jj=%d Erfrad=%g Erffluid=%g gammarad=%g gammafluid=%g Erfnew=%g urfconrel=%g\n",jj,Erfrad,Erffluid,gammarad,gammafluid,*Erfnew,urfconrel[jj]);

  return(0);
}



// get's gamma^2 for lab-frame gamma
static int get_m1closure_gammarel2_old(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn)
{
  FTYPE gamma2,gammarel2,delta,numerator,divisor;

  if(0){
    // has some catastrophic cancellation issue for non-moving velocity at very low E\sim 1E-92 (as in RADPULSE test if no temperature conversion)

    //g_munu R^tmu R^tnu
    int jj,kk;
    FTYPE gRR=0.0;
    DLOOP(jj,kk) gRR += ptrgeom->gcov[GIND(jj,kk)]*Avcon[jj]*Avcon[kk];

    //the quadratic equation for u^t of the radiation rest frame (urf[0])
    // Formed as solution for solving two equations (R^{t\nu} R^t_\nu(E,ut) and R^{tt}(E,ut)) for ut
    //supposed to provide two roots for (u^t)^2 of opposite signs
    FTYPE a,b,c;
    a=16.*gRR;
    b=8.*(gRR*ptrgeom->gcon[GIND(0,0)]+Avcon[0]*Avcon[0]);
    c=ptrgeom->gcon[GIND(0,0)]*(gRR*ptrgeom->gcon[GIND(0,0)]-Avcon[0]*Avcon[0]);
    delta=b*b-4.*a*c;

    numerator=0.5*(-b-sqrt(delta));
    divisor=a;

    gamma2=numerator/divisor; // lab-frame gamma^2
    //if unphysical try the other root
    if(gamma2<=0.){
      numerator=0.5*(-b+sqrt(delta));
      divisor=a;
      gamma2=  numerator/divisor; 
    }
    
    *numeratorreturn=numerator;
    *divisorreturn=divisor;
  }
  //    dualfprintf(fail_file,"GAMMA2CHECK: ijk=%d %d %d : %g %g : a=%g b=%g c=%g : delta=%g gRR=%g Avcon0123=%g %g %g %g : gamma2=%g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,0.5*(-b-sqrt(delta))/a,0.5*(-b+sqrt(delta))/a,a,b,c,delta,gRR,Avcon[0],Avcon[1],Avcon[2],Avcon[3],gamma2);


  else{
    // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
    // well, then had problems for R~1E-14 for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
    // see u2p_inversion.nb
    static FTYPE gctt, gv11, gv12,  gv13,  gv14,  gv22,  gv23,  gv24,  gv33,  gv34,  gv44,  Rtt,  Rtx,  Rty,  Rtz;
    gv11=ptrgeom->gcov[GIND(0,0)];
    gv12=ptrgeom->gcov[GIND(0,1)];
    gv13=ptrgeom->gcov[GIND(0,2)];
    gv14=ptrgeom->gcov[GIND(0,3)];
    gv22=ptrgeom->gcov[GIND(1,1)];
    gv23=ptrgeom->gcov[GIND(1,2)];
    gv24=ptrgeom->gcov[GIND(1,3)];
    gv33=ptrgeom->gcov[GIND(2,2)];
    gv34=ptrgeom->gcov[GIND(2,3)];
    gv44=ptrgeom->gcov[GIND(3,3)];
    Rtt=Avcon[0];
    Rtx=Avcon[1];
    Rty=Avcon[2];
    Rtz=Avcon[3];
    gctt=ptrgeom->gcon[GIND(0,0)];

    delta = (1. + 3.*gctt*gv11)*((Rtt)*(Rtt)) + 
      6.*gctt*Rtt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 
      3.*gctt*(gv22*((Rtx)*(Rtx)) + 2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 
               2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*((Rtz)*(Rtz)));

    divisor=(gv11*((Rtt)*(Rtt)) + 2.*gv12*Rtt*Rtx + gv22*((Rtx)*(Rtx)) + 2.*gv13*Rtt*Rty + 
             2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 2.*(gv14*Rtt + gv24*Rtx + gv34*Rty)*Rtz + 
             gv44*((Rtz)*(Rtz)));

    numerator=(-0.25*((1. + gctt*gv11)*((Rtt)*(Rtt)) + 
                      gctt*(gv22*((Rtx)*(Rtx)) + 2.*gv23*Rtx*Rty + gv33*((Rty)*(Rty)) + 
                            2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*((Rtz)*(Rtz))) + 
                      Rtt*(2.*gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 
                           Sqrt(delta))));

    gamma2 = numerator/divisor;
  }


  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  FTYPE alpha=ptrgeom->alphalapse;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  // check for machine error away from 1.0 that happens sometimes
  if(gammarel2>GAMMASMALLLIMIT && gammarel2<1.0){
    // if(debugfail>=2) dualfprintf(fail_file,"Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);
    gammarel2=1.0;
  }

  //  dualfprintf(fail_file,"gammarel2=%g gamma2=%g delta=%21.15g\n",gammarel2,gamma2,delta);

  *gammarel2return=gammarel2;
  *deltareturn=delta;
  *numeratorreturn=numerator;
  *divisorreturn=divisor;
  return(0);

}












// get's gamma^2 for lab-frame gamma  using Rd and gcon
static int get_m1closure_gammarel2(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn)
{
  FTYPE gamma2,gammarel2,delta,numerator,divisor;
  FTYPE gamma2a,gamma2b;

  // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
  // well, then had problems for R~1E-14L for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
  // see u2p_inversion.nb
  static FTYPE gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=ptrgeom->gcon[GIND(0,0)];
  gn12=ptrgeom->gcon[GIND(0,1)];
  gn13=ptrgeom->gcon[GIND(0,2)];
  gn14=ptrgeom->gcon[GIND(0,3)];
  gn22=ptrgeom->gcon[GIND(1,1)];
  gn23=ptrgeom->gcon[GIND(1,2)];
  gn24=ptrgeom->gcon[GIND(1,3)];
  gn33=ptrgeom->gcon[GIND(2,2)];
  gn34=ptrgeom->gcon[GIND(2,3)];
  gn44=ptrgeom->gcon[GIND(3,3)];

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  gamma2a=(-0.25*(2.*Power(gn11,2)*Power(Rdtt,2) + (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*
        (gn12*Rdtx + gn13*Rdty + gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
            gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
               8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))) + 
       gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 
          2.*gn34*Rdty*Rdtz + gn44*Power(Rdtz,2) + Rdtt*
           (4.*gn13*Rdty + 4.*gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
               gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                  8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))))))/
   (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
    2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));


  if( gamma2a<GAMMASMALLLIMIT || !isfinite(gamma2a) ){
    gamma2b=(0.25*(-2.*Power(gn11,2)*Power(Rdtt,2) - 1.*gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 
                                                              Rdty*(4.*gn13*Rdtt + 2.*gn23*Rdtx + gn33*Rdty) + 2.*(2.*gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2)) + 
                   gn11*Rdtt*Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                  gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                        8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))) + 
                   (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*(-1.*gn12*Rdtx - 1.*gn13*Rdty - 1.*gn14*Rdtz + 
                                                        Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                                             gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                                                   8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))))))/
      (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
       2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));
    gamma2=gamma2b;
  }
  else{
    // choose
    gamma2=gamma2a;
  }

  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  FTYPE alpha=ptrgeom->alphalapse;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  // check for machine error away from 1.0 that happens sometimes
  if(gammarel2>GAMMASMALLLIMIT && gammarel2<1.0){
    // if(debugfail>=2) dualfprintf(fail_file,"Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);
    gammarel2=1.0;
  }

  //  dualfprintf(fail_file,"gammarel2=%g gamma2=%g delta=%21.15g\n",gammarel2,gamma2,delta);

  *gammarel2return=gammarel2;
  *deltareturn=delta=0;
  *numeratorreturn=numerator=0;
  *divisorreturn=divisor=0;
  return(0);

}







// get Erf
static int get_m1closure_Erf(struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE gammarel2, FTYPE *Erfreturn)
{
  FTYPE alpha=ptrgeom->alphalapse;

  ////////////
  //
  // get initial attempt for Erf
  // If delta<0, then gammarel2=nan and Erf<RADLIMIT check below will fail as good.
  //
  ////////////
  *Erfreturn = 3.*Avcon[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  return(0);
}



// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel_old(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  FTYPE Erf=*Erfreturn; // get initial Erf
  FTYPE gammamax=GAMMAMAXRAD;
  FTYPE gammamaxfail=GAMMAMAXRADFAIL;
  int jj,kk;


  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  //////////////////////
  //
  // First case is if gammarel>gammamax, then set gammarel=gammamax unless Erf<ERADLIMIT (~0) in which case set Erf=ERADLIMIT and gammarel=1.
  // Note, can't set urfcon[0]=gammamax in case gammamax still remains space-like, e.g. inside horizon if gammamax isn't big enough.
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  int nonfailure=gammarel2>=1.0 && Erf>ERADLIMIT && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADLIMIT (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Avcon[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25L && delta>=0.0 && divisor!=0.0) || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADLIMIT;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=gammarel2>gammamax*gammamax && Erf>=ERADLIMIT || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;



  if(nonfailure){
    // get good relative velocity
    FTYPE gammarel=sqrt(gammarel2);
    FTYPE alpha=ptrgeom->alphalapse;

    SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);

    *Erfreturn=Erf; // pass back new Erf to pointer
    return(0);

    //        dualfprintf(fail_file,"NO failure: %g %g ijk=%d %d %d\n",Erf,gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(failure1){
    if(TRYCOLD){
      gammarel2=pow(1.0+10.0*NUMEPSILON,2.0);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
    }
    else{
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;
     
      SLOOPA(jj) urfconrel[jj] = 0.0; // consistent with gammarel2=1
    }
    if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE3A;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE3A: normal gamma, but Erf<ERADLIMIT. ijk=%d %d %d : %ld %d %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);

  }    
  else if(failure2){
    if(TRYCOLD){
      gammarel2=pow(1.0+10.0*NUMEPSILON,2.0);
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
    }
    else{
      FTYPE gammarel2orig=gammarel2;
      // override
      gammarel2=1.0;
      FTYPE gammarel=1.0;  // use this below

      // get new Erf(gammarel)
      get_m1closure_Erf(ptrgeom, Avcon, gammarel2, &Erf);
      if(Erf<ERADLIMIT)  Erf=ERADLIMIT;
    
      SLOOPA(jj) urfconrel[jj] = 0.0;

    }
    if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2A;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2A: normal gamma, but Erf<ERADLIMIT. ijk=%d %d %d : %ld %d %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
  }
  else{
    if(TRYCOLD){
      gammarel2=gammamax*gammamax;
      get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor,&Erf,urfconrel);
      if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax (cold): gammarel2=%g Erf=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,Erf,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
    }
    else{
      FTYPE gammarel2orig=gammarel2;
      FTYPE gammarel=gammamax;
      gammarel2=gammamax*gammamax;

      // get new Erf(gammarel)
      get_m1closure_Erf(ptrgeom, Avcon, gammarel2, &Erf);


      // Check if Erf is too small with gamma->gammamax
      if(Erf<ERADLIMIT){
        if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1A;
        // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
        Erf=ERADLIMIT;

        // can't use normal velocity with small Erf -- fails with inf or nan
        SLOOPA(jj) urfconrel[jj] = 0.0;

        if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
      }
      else{
        // if Erf normal, assume ok to have gammamax for radiation.  This avoids fixups, which can generate more oscillations.
        // KORALTODO: But note that then check_on_inversion() won't know that failure and will check and report issue.
        if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammarel2=%g : i=%d j=%d k=%d : %ld %d %g\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);

        // regardless of Erf value, now that have some Erf, ensure gamma=gammamax
        // lab-frame radiation relative 4-velocity
        FTYPE alpha=ptrgeom->alphalapse;
        SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
        
        // compute \gammarel using this (gammatemp can be inf if Erf=ERADLIMIT, and then rescaling below will give urfconrel=0 and gammarel=1
        FTYPE gammatemp,qsqtemp;
        int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);
        MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);

        if(!isfinite(gammatemp)){
          SLOOPA(jj) urfconrel[jj] =0.0;
        }
        else if(0&&gammatemp<=gammamax){
          // do nothing, don't make gamma larger just to get consistency
        }
        else{
          // now rescale urfconrel[i] so will give desired \gammamax
          SLOOPA(jj) urfconrel[jj] *= (gammamax/gammatemp);
        }
 
#if(PRODUCTION==0)
        // check that gamma really correctly gammamax
        FTYPE gammatemp2,qsqtemp2;
        MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp2,&qsqtemp2),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
        if(showmessages) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammarel2orig=%g gammamax=%g gammatemp=%g gammatemp2=%g ijk=%d %d %d : %ld %d %g\n",gammarel2orig,gammamax,gammatemp,gammatemp2,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
#endif
      }
      //    if(showmessages && debugfail>=2) DLOOPA(jj) dualfprintf(fail_file,"CASE1B: urfconrel[%d]=%g uu[%d]=%g\n",jj,urfconrel[jj],jj,uu[URAD0+jj]);

      //    SLOOPA(jj) urfconrel[jj] = 0.0; // consistent with gammarel2=1
    }
  }


  // if here, then one failure mode.  See if optically thick or thin and reduce to (e.g.) fluid frame if thick

  // can't use normal velocity with small Erf -- fails with inf or nan
  // setup "old" pp in case used
  pp[PRAD0] = Erf;
  SLOOPA(jj) pp[PRAD1+jj-1] = urfconrel[jj];


#if(0)
  // KORALTODO: Problems when tau<<1 and gamma->gammamax
  // KORALTODO: DUH, also should only do if failure, not if no failure.
  FTYPE tautot[NDIM],tautotmax;
  if(M1REDUCE==TOOPACITYDEPENDENTFRAME){
    // then will possibly need tautotmax
    // get tautot based upon previous pp in order to determine what to do in case of failure
    calc_tautot(pp, ptrgeom, tautot, &tautotmax);
  }

  
  if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
  else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;
  else if(M1REDUCE==TOOPACITYDEPENDENTFRAME) opacity_interpolated_urfconrel(tautotmax,pp,ptrgeom,Avcon,Erf,gammarel2,&Erf,urfconrel);
#endif     

  *Erfreturn=Erf; // pass back new Erf to pointer
  return(0);
}




// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE numerator, FTYPE divisor, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  FTYPE Erf=*Erfreturn; // get initial Erf
  FTYPE gammamax=GAMMAMAXRAD;
  FTYPE gammamaxfail=GAMMAMAXRADFAIL;
  int jj,kk;


  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  //  int nonfailure=gammarel2>=1.0 && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  int nonfailure=gammarel2>=1.0L && Erf>ERADLIMIT && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADLIMIT (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Avcon[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25L && delta>=0.0 && divisor!=0.0) || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADLIMIT;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=gammarel2>gammamax*gammamax && Erf>=ERADLIMIT || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;

  // any failure
  int failure=!nonfailure || !isfinite(gammarel2) || !isfinite(Erf);

  if(failure && (failure1==0 && failure2==0 && failure3==0)){
    if(debugfail>=2) dualfprintf(fail_file,"Undetected failure, now considered\n");
  }


  FTYPE Erf0=MAX(NUMEPSILON*fabs(Avcon[TT])/gammamax/gammamax,ERADLIMIT);
  if(nonfailure){
    // get good relative velocity
    FTYPE gammarel=sqrt(gammarel2);
    FTYPE alpha=ptrgeom->alphalapse;

    SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);


    *Erfreturn=Erf; // pass back new Erf to pointer
    return(0);

    //        dualfprintf(fail_file,"NO failure: %g %g ijk=%d %d %d\n",Erf,gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else{
    FTYPE Avconorig[NDIM],Avcovorig[NDIM];
    DLOOPA(jj){
      Avconorig[jj]=Avcon[jj];
      Avcovorig[jj]=Avcov[jj];
    }
    FTYPE gammarel2orig;
    gammarel2orig=gammarel2;
    FTYPE Erforig;
    Erforig=Erf;

    // get \gammarel=1 case
    FTYPE gammarel2slow=pow(1.0+10.0*NUMEPSILON,2.0);
    FTYPE Avconslow[NDIM],Avcovslow[NDIM],Erfslow,urfconrelslow[NDIM];
    DLOOPA(jj){
      Avconslow[jj]=Avcon[jj];
      Avcovslow[jj]=Avcov[jj];
    }
    Erfslow=Erf;
    get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avconslow,Avcovslow,&gammarel2slow,&delta,&numerator,&divisor,&Erfslow,urfconrelslow);

    // get \gammarel=gammamax case
    FTYPE gammarel2fast=gammamax*gammamax;
    FTYPE Avconfast[NDIM],Avcovfast[NDIM],Erffast,urfconrelfast[NDIM];
    DLOOPA(jj){
      Avconfast[jj]=Avcon[jj];
      Avcovfast[jj]=Avcov[jj];
    }
    Erffast=Erf;
    get_m1closure_gammarel2_cold(showmessages,ptrgeom,Avconfast,Avcovfast,&gammarel2fast,&delta,&numerator,&divisor,&Erffast,urfconrelfast);

    //    dualfprintf(fail_file,"JONVSOLEK: Avconorig: %g %g %g %g\n",Avcon[0],Avcon[1],Avcon[2],Avcon[3]);

    int usingfast=1;
    // choose by which Avcov[0] is closest to original
    if( fabs(Avconslow[0]-Avcon[0])>fabs(Avconfast[0]-Avcon[0]) ){ // compare Avcon that has positive sign always
      usingfast=1;
      Erf=Erffast;
      gammarel2=gammarel2fast;
      DLOOPA(jj){
        Avcon[jj]=Avconfast[jj];
        Avcov[jj]=Avcovfast[jj];
        urfconrel[jj]=urfconrelfast[jj];
      }
    }
    else{
      usingfast=0;
      Erf=Erfslow;
      gammarel2=gammarel2slow;
      DLOOPA(jj){
        Avcon[jj]=Avconslow[jj];
        Avcov[jj]=Avcovslow[jj];
        urfconrel[jj]=urfconrelslow[jj];
      }
    }

    // catch bad issue for when using fast or slow will be bad because probably momentum is bad if inverted energy
    if(Avcovorig[TT]>0.0){
      SLOOPA(jj) urfconrel[jj]=0.0;
    }

    // report
    if(1||allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1A;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASEGEN: gammarel>gammamax (cold, usingfast=%d): gammarel2=%g Erf=%g : i=%d j=%d k=%d : %ld %d %g\n",usingfast,gammarel2,Erf,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,t);
  }


  // if here, then one failure mode.  See if optically thick or thin and reduce to (e.g.) fluid frame if thick

  // can't use normal velocity with small Erf -- fails with inf or nan
  // setup "old" pp in case used
  pp[PRAD0] = Erf;
  SLOOPA(jj) pp[PRAD1+jj-1] = urfconrel[jj];



#if(0)
  // KORALTODO: Problems when tau<<1 and gamma->gammamax
  FTYPE tautot[NDIM],tautotmax;
  if(M1REDUCE==TOOPACITYDEPENDENTFRAME){
    // then will possibly need tautotmax
    // get tautot based upon previous pp in order to determine what to do in case of failure
    calc_tautot(pp, ptrgeom, tautot, &tautotmax);
  }

  
  if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
  else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;
  else if(M1REDUCE==TOOPACITYDEPENDENTFRAME) opacity_interpolated_urfconrel(tautotmax,pp,ptrgeom,Avcon,Erf,gammarel2,&Erf,urfconrel);
#endif     

  *Erfreturn=Erf; // pass back new Erf to pointer


  // catch any nan/inf's:
  int notfinite=(!isfinite(Erf) || !isfinite(urfconrel[1])|| !isfinite(urfconrel[2])|| !isfinite(urfconrel[3]));
  if(notfinite){
    // nothing else to do unless want to use nan/inf as indicator that should abort something
    // using such a small Erf can lead itself to problems due to precision issues, so assume will fixup this
    Erf=ERADLIMIT;
    SLOOPA(jj) urfconrel[jj]=0.0; // ZAMO
    if(1||allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
  }



  // DEBUG:
  if(debugfail>=2){
    if(notfinite){
      dualfprintf(fail_file,"JONNAN: ijk=%d %d %d :  %g %g : %g %g %g : %d %d %d %d : %g %g %g %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,Erf,gammarel2,urfconrel[1],urfconrel[2],urfconrel[3],failure1,failure2,failure3,failure,Avcon[0],Avcon[1],Avcon[2],Avcon[3]);
    }
  }


  return(0);
}




// get contravariant relative 4-velocity in lab frame using Olek's koral choices
static int get_m1closure_urfconrel_olek(int showmessages, int allowlocalfailurefixandnoreport, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *Avcon, FTYPE *Avcov, FTYPE gammarel2, FTYPE delta, FTYPE *Erfreturn, FTYPE *urfconrel, PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  FTYPE Erf=*Erfreturn; // get initial Erf
  FTYPE gammamax=GAMMAMAXRAD;
  FTYPE gammamaxfail=GAMMAMAXRADFAIL;
  int jj,kk;


  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  int failure1=gammarel2>1.01*gammamax*gammamax || gammarel2<0. || delta<0.;
  int failure2=gammarel2<1. || delta<0. || !isfinite(gammarel2); // NOTE: first failure1 already catches delta<0.



  if(failure1){
    //    if(failure1 && tautotmax<TAUFAILLIMIT){ // works for DBLSHADOW

    FTYPE gammarel=gammamax;
    gammarel2=gammamax*gammamax;

    // get new Erf(gammarel)
    get_m1closure_Erf(ptrgeom, Avcon, gammarel2, &Erf);


    // Check if Erf is too small with gamma->gammamax
    if(Erf<ERADLIMIT || !isfinite(Erf)){
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1A;
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;

      // can't use normal velocity with small Erf -- fails with inf or nan
      SLOOPA(jj) urfconrel[jj] = 0.0;

      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }
    else{
      // if Erf normal, assume ok to have gammamax for radiation.  This avoids fixups, which can generate more oscillations.
      if(allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE1B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);

      // regardless of Erf value, now that have some Erf, ensure gamma=gammamax
      // lab-frame radiation relative 4-velocity
      FTYPE alpha=ptrgeom->alphalapse;
      SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
        
      // compute \gammarel using this (gammatemp can be inf if Erf=ERADLIMIT, and then rescaling below will give urfconrel=0 and gammarel=1
      FTYPE gammatemp,qsqtemp;
      int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);
      MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
        
      // now rescale urfconrel[i] so will give desired \gammamax
      SLOOPA(jj) urfconrel[jj] *= (gammamax/gammatemp);
 
#if(PRODUCTION==0)
      // check that gamma really correctly gammamax
      FTYPE gammatemp2,qsqtemp2;
      MYFUN(gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammatemp2,&qsqtemp2),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
      if(showmessages) dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammamax=%g gammatemp=%g gammatemp2=%g ijk=%d %d %d\n",gammamax,gammatemp,gammatemp2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif
    }

  }
  //////////////////////
  //
  // Second case is if gammarel<1 or delta<0, then set gammarel=1.  If Erf<ERADLIMIT (~0), then set Erf=ERADLIMIT and gammarel=1.
  // Can't assume this condition is equivalent to large gamma, because if not, then leads to crazy boost of energy.
  //
  //////////////////////
  else if(failure2){


    FTYPE gammarel2orig=gammarel2;
    // override
    gammarel2=1.0;
    FTYPE gammarel=1.0;  // use this below

    // get new Erf(gammarel)
    get_m1closure_Erf(ptrgeom, Avcon, gammarel2, &Erf);
    SLOOPA(jj) urfconrel[jj] = 0.0;


    if(Erf<ERADLIMIT || !isfinite(Erf)){ // JCM
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2A;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2A: gamma<1 or delta<0 and Erf<ERADLIMIT : gammarel2=%g : i=%d j=%d k=%d\n",gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }
    else{
      // normal Erf
      if(1 || allowlocalfailurefixandnoreport==0) *lpflagrad=UTOPRIMRADFAILCASE2B;
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"CASE2B: gamma<1 or delta<0 and Erf normal : gammamax=%g gammarel2orig=%21.15g gammarel2=%21.15g delta=%g : i=%d j=%d k=%d\n",gammamax,gammarel2orig,gammarel2,delta,ptrgeom->i,ptrgeom->j,ptrgeom->k);
    }


      
  }
  //////////////////////
  //
  // Third case is if no bad conditions, then try regular calculation.  If Erf<ERADLIMIT, then already caught with first condition
  //
  //////////////////////
  else{

    if(Erf<ERADLIMIT || !isfinite(Erf)){
      Erf=ERADLIMIT;
      SLOOPA(jj) urfconrel[jj] = 0.0;
      // must use above because if use ERADLIMIT in normal urfconrel, then urfconrel will be HUGE and probably give inf or nan due to Avcon/(Erf*gammarel) term.
    }
    else{
      // get good relative velocity
      FTYPE gammarel=sqrt(gammarel2);
      FTYPE alpha=ptrgeom->alphalapse;
      
      // NOTEMARK: This overwrites choice above for urfconrel when Erf<ERADLIMIT.
      SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
    }


    //        dualfprintf(fail_file,"NO failure: %g %g ijk=%d %d %d\n",Erf,gammarel2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }


  if(debugfail>=2){
    if(!isfinite(Erf) || !isfinite(gammarel2) || !isfinite(urfconrel[0])|| !isfinite(urfconrel[1])|| !isfinite(urfconrel[2])|| !isfinite(urfconrel[3]) ){
      dualfprintf(fail_file,"OLEKNAN: ijk=%d %d %d :  %g %g : %g %g %g : %d %d : %g %g %g %g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,Erf,gammarel2,urfconrel[1],urfconrel[2],urfconrel[3],failure1,failure2,Avcon[0],Avcon[1],Avcon[2],Avcon[3]);
    }
  }
      
  *Erfreturn=Erf; // pass back new Erf to pointer
  return(0);
}








// get's gamma assuming fixed E rather than using original R^{tt} that we assume is flawed near floor regions.  We want to preserve R^{ti} (i.e momentum)
static int get_m1closure_gammarel2_cold_old(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel)
{
  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=ptrgeom->alphalapse;
  int jj;

  static FTYPE gctt, gv11, gv12,  gv13,  gv14,  gv22,  gv23,  gv24,  gv33,  gv34,  gv44,  Rtt,  Rtx,  Rty,  Rtz;
  gv11=ptrgeom->gcov[GIND(0,0)];
  gv12=ptrgeom->gcov[GIND(0,1)];
  gv13=ptrgeom->gcov[GIND(0,2)];
  gv14=ptrgeom->gcov[GIND(0,3)];
  gv22=ptrgeom->gcov[GIND(1,1)];
  gv23=ptrgeom->gcov[GIND(1,2)];
  gv24=ptrgeom->gcov[GIND(1,3)];
  gv33=ptrgeom->gcov[GIND(2,2)];
  gv34=ptrgeom->gcov[GIND(2,3)];
  gv44=ptrgeom->gcov[GIND(3,3)];
  FTYPE Rttold=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];
  gctt=ptrgeom->gcon[GIND(0,0)];


  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamaxfail=GAMMAMAXRADFAIL;
    FTYPE gammamax=GAMMAMAXRAD;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);



  // but check if utsq is too small
  FTYPE utsqmina=(0.5*(-8.*gctt*Power(gv12,2)*Power(Rtx,2) + 8.*gv22*Power(Rtx,2) + 
                       8.*gctt*gv11*gv22*Power(Rtx,2) - 16.*gctt*gv12*gv13*Rtx*Rty + 
                       16.*gv23*Rtx*Rty + 16.*gctt*gv11*gv23*Rtx*Rty - 
                       8.*gctt*Power(gv13,2)*Power(Rty,2) + 8.*gv33*Power(Rty,2) + 
                       8.*gctt*gv11*gv33*Power(Rty,2) - 16.*gctt*gv12*gv14*Rtx*Rtz + 
                       16.*gv24*Rtx*Rtz + 16.*gctt*gv11*gv24*Rtx*Rtz - 
                       16.*gctt*gv13*gv14*Rty*Rtz + 16.*gv34*Rty*Rtz + 
                       16.*gctt*gv11*gv34*Rty*Rtz - 8.*gctt*Power(gv14,2)*Power(Rtz,2) + 
                       8.*gv44*Power(Rtz,2) + 8.*gctt*gv11*gv44*Power(Rtz,2) - 
                       1.*Sqrt(Power(8.*gctt*Power(gv12,2)*Power(Rtx,2) - 8.*gv22*Power(Rtx,2) - 
                                     8.*gctt*gv11*gv22*Power(Rtx,2) + 16.*gctt*gv12*gv13*Rtx*Rty - 
                                     16.*gv23*Rtx*Rty - 16.*gctt*gv11*gv23*Rtx*Rty + 
                                     8.*gctt*Power(gv13,2)*Power(Rty,2) - 8.*gv33*Power(Rty,2) - 
                                     8.*gctt*gv11*gv33*Power(Rty,2) + 16.*gctt*gv12*gv14*Rtx*Rtz - 
                                     16.*gv24*Rtx*Rtz - 16.*gctt*gv11*gv24*Rtx*Rtz + 
                                     16.*gctt*gv13*gv14*Rty*Rtz - 16.*gv34*Rty*Rtz - 
                                     16.*gctt*gv11*gv34*Rty*Rtz + 8.*gctt*Power(gv14,2)*Power(Rtz,2) - 
                                     8.*gv44*Power(Rtz,2) - 8.*gctt*gv11*gv44*Power(Rtz,2),2) - 
                               4.*(16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
                                   32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
                                   16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
                                   32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 
                                   32.*gv13*gv14*Rty*Rtz - 32.*gv11*gv34*Rty*Rtz + 
                                   16.*Power(gv14,2)*Power(Rtz,2) - 16.*gv11*gv44*Power(Rtz,2))*
                               (Power(gctt,2)*Power(gv12,2)*Power(Rtx,2) + gctt*gv22*Power(Rtx,2) - 
                                1.*Power(gctt,2)*gv11*gv22*Power(Rtx,2) + 
                                2.*Power(gctt,2)*gv12*gv13*Rtx*Rty + 2.*gctt*gv23*Rtx*Rty - 
                                2.*Power(gctt,2)*gv11*gv23*Rtx*Rty + 
                                Power(gctt,2)*Power(gv13,2)*Power(Rty,2) + gctt*gv33*Power(Rty,2) - 
                                1.*Power(gctt,2)*gv11*gv33*Power(Rty,2) + 
                                2.*Power(gctt,2)*gv12*gv14*Rtx*Rtz + 2.*gctt*gv24*Rtx*Rtz - 
                                2.*Power(gctt,2)*gv11*gv24*Rtx*Rtz + 
                                2.*Power(gctt,2)*gv13*gv14*Rty*Rtz + 2.*gctt*gv34*Rty*Rtz - 
                                2.*Power(gctt,2)*gv11*gv34*Rty*Rtz + 
                                Power(gctt,2)*Power(gv14,2)*Power(Rtz,2) + gctt*gv44*Power(Rtz,2) - 
                                1.*Power(gctt,2)*gv11*gv44*Power(Rtz,2)))))/
    (16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
     32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
     16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
     32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 32.*gv13*gv14*Rty*Rtz - 
     32.*gv11*gv34*Rty*Rtz + 16.*Power(gv14,2)*Power(Rtz,2) - 
     16.*gv11*gv44*Power(Rtz,2));

  FTYPE utsqminb=(0.5*(-8.*gctt*Power(gv12,2)*Power(Rtx,2) + 8.*gv22*Power(Rtx,2) + 
                       8.*gctt*gv11*gv22*Power(Rtx,2) - 16.*gctt*gv12*gv13*Rtx*Rty + 
                       16.*gv23*Rtx*Rty + 16.*gctt*gv11*gv23*Rtx*Rty - 
                       8.*gctt*Power(gv13,2)*Power(Rty,2) + 8.*gv33*Power(Rty,2) + 
                       8.*gctt*gv11*gv33*Power(Rty,2) - 16.*gctt*gv12*gv14*Rtx*Rtz + 
                       16.*gv24*Rtx*Rtz + 16.*gctt*gv11*gv24*Rtx*Rtz - 
                       16.*gctt*gv13*gv14*Rty*Rtz + 16.*gv34*Rty*Rtz + 
                       16.*gctt*gv11*gv34*Rty*Rtz - 8.*gctt*Power(gv14,2)*Power(Rtz,2) + 
                       8.*gv44*Power(Rtz,2) + 8.*gctt*gv11*gv44*Power(Rtz,2) + 
                       Sqrt(Power(8.*gctt*Power(gv12,2)*Power(Rtx,2) - 8.*gv22*Power(Rtx,2) - 
                                  8.*gctt*gv11*gv22*Power(Rtx,2) + 16.*gctt*gv12*gv13*Rtx*Rty - 
                                  16.*gv23*Rtx*Rty - 16.*gctt*gv11*gv23*Rtx*Rty + 
                                  8.*gctt*Power(gv13,2)*Power(Rty,2) - 8.*gv33*Power(Rty,2) - 
                                  8.*gctt*gv11*gv33*Power(Rty,2) + 16.*gctt*gv12*gv14*Rtx*Rtz - 
                                  16.*gv24*Rtx*Rtz - 16.*gctt*gv11*gv24*Rtx*Rtz + 
                                  16.*gctt*gv13*gv14*Rty*Rtz - 16.*gv34*Rty*Rtz - 
                                  16.*gctt*gv11*gv34*Rty*Rtz + 8.*gctt*Power(gv14,2)*Power(Rtz,2) - 
                                  8.*gv44*Power(Rtz,2) - 8.*gctt*gv11*gv44*Power(Rtz,2),2) - 
                            4.*(16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
                                32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
                                16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
                                32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 
                                32.*gv13*gv14*Rty*Rtz - 32.*gv11*gv34*Rty*Rtz + 
                                16.*Power(gv14,2)*Power(Rtz,2) - 16.*gv11*gv44*Power(Rtz,2))*
                            (Power(gctt,2)*Power(gv12,2)*Power(Rtx,2) + gctt*gv22*Power(Rtx,2) - 
                             1.*Power(gctt,2)*gv11*gv22*Power(Rtx,2) + 
                             2.*Power(gctt,2)*gv12*gv13*Rtx*Rty + 2.*gctt*gv23*Rtx*Rty - 
                             2.*Power(gctt,2)*gv11*gv23*Rtx*Rty + 
                             Power(gctt,2)*Power(gv13,2)*Power(Rty,2) + gctt*gv33*Power(Rty,2) - 
                             1.*Power(gctt,2)*gv11*gv33*Power(Rty,2) + 
                             2.*Power(gctt,2)*gv12*gv14*Rtx*Rtz + 2.*gctt*gv24*Rtx*Rtz - 
                             2.*Power(gctt,2)*gv11*gv24*Rtx*Rtz + 
                             2.*Power(gctt,2)*gv13*gv14*Rty*Rtz + 2.*gctt*gv34*Rty*Rtz - 
                             2.*Power(gctt,2)*gv11*gv34*Rty*Rtz + 
                             Power(gctt,2)*Power(gv14,2)*Power(Rtz,2) + gctt*gv44*Power(Rtz,2) - 
                             1.*Power(gctt,2)*gv11*gv44*Power(Rtz,2)))))/
    (16.*Power(gv12,2)*Power(Rtx,2) - 16.*gv11*gv22*Power(Rtx,2) + 
     32.*gv12*gv13*Rtx*Rty - 32.*gv11*gv23*Rtx*Rty + 
     16.*Power(gv13,2)*Power(Rty,2) - 16.*gv11*gv33*Power(Rty,2) + 
     32.*gv12*gv14*Rtx*Rtz - 32.*gv11*gv24*Rtx*Rtz + 32.*gv13*gv14*Rty*Rtz - 
     32.*gv11*gv34*Rty*Rtz + 16.*Power(gv14,2)*Power(Rtz,2) - 
     16.*gv11*gv44*Power(Rtz,2));

  dualfprintf(fail_file,"utsq=%g utsqmina=%g utsqminb=%g\n",utsq,utsqmina,utsqminb);
  // KORALTODO: override (only applicable for first root)  Unsure if 2nd root used for GR in ergosphere.  e.g. gv11 switches sign!
  if(utsq<utsqmina && utsqmina>utsqminb) utsq=utsqmina;
  if(utsq<utsqminb && utsqminb>utsqmina) utsq=utsqminb;


  FTYPE Avcovorig[NDIM];
  DLOOPA(jj) Avcovorig[jj]=Avcov[jj];
  

  // get new Avcon[0]=R^{tt}
  
  Avcon[0]=(-1.*(gctt + 4.*utsq)*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
                               0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
                                                        36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
                                                        (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
    (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));

  Erf=(-3.*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
            0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
                                     36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
                                     (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
    (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));


  dualfprintf(fail_file,"NOR SOL: Avcon0new=%g Avcon0old=%g Erf=%g :: %g %g %g\n",Avcon[0],Rttold,Erf,Rtx,Rty,Rtz);
 
  FTYPE Avcovnew[NDIM];
  indices_21(Avcon,Avcovnew,ptrgeom);
  DLOOPA(jj) dualfprintf(fail_file,"jj=%d Avcovorig=%g Avcovnew=%g\n",jj,Avcovorig[jj],Avcovnew[jj]);

 
  delta=0; // not yet

  if(1){
    // alt solution
    FTYPE Avalt=(-1.*(gctt + 4.*utsq)*(gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) + 4.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq - 
                                       0.16666666666666666*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
                                                                36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
                                                                (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq)))))/
      (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));


    FTYPE Erfalt=(-3.*gctt*(gv12*Rtx + gv13*Rty + gv14*Rtz) - 12.*(gv12*Rtx + gv13*Rty + gv14*Rtz)*utsq + 
                  0.5*Sqrt(36.*Power(gv12*Rtx + gv13*Rty + gv14*Rtz,2)*Power(gctt + 4.*utsq,2) - 
                           36.*(gv22*Power(Rtx,2) + 2.*gv23*Rtx*Rty + gv33*Power(Rty,2) + 2.*gv24*Rtx*Rtz + 2.*gv34*Rty*Rtz + gv44*Power(Rtz,2))*
                           (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq))))/
      (Power(gctt,2)*gv11 + 8.*utsq*(1. + 2.*gv11*utsq) + gctt*(-1. + 8.*gv11*utsq));

    dualfprintf(fail_file,"ALT SOL: Avalt=%g Av0old=%g Erfalt=%g : %g %g %g\n",Avalt,Rttold,Erfalt,Rtx,Rty,Rtz);
  }


  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  dualfprintf(fail_file,"NORM ROOT 4-vel: %g %g %g : %g\n",urfconrel[1],urfconrel[2],urfconrel[3],ptrgeom->gdet);

  
  *Erfreturn=Erf; // pass back new Erf to pointer


  return(0);
}











// get's gamma assuming fixed E rather than using original R^t_t that we assume is flawed near floor regions.  We want to preserve R^t_i (i.e conserved momentum)
static int get_m1closure_gammarel2_cold(int showmessages, struct of_geom *ptrgeom, FTYPE *Avcon, FTYPE *Avcov, FTYPE *gammarel2return, FTYPE *deltareturn, FTYPE *numeratorreturn, FTYPE *divisorreturn, FTYPE *Erfreturn, FTYPE *urfconrel)
{
  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=ptrgeom->alphalapse;
  int jj;

  static FTYPE gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=ptrgeom->gcon[GIND(0,0)];
  gn12=ptrgeom->gcon[GIND(0,1)];
  gn13=ptrgeom->gcon[GIND(0,2)];
  gn14=ptrgeom->gcon[GIND(0,3)];
  gn22=ptrgeom->gcon[GIND(1,1)];
  gn23=ptrgeom->gcon[GIND(1,2)];
  gn24=ptrgeom->gcon[GIND(1,3)];
  gn33=ptrgeom->gcon[GIND(2,2)];
  gn34=ptrgeom->gcon[GIND(2,3)];
  gn44=ptrgeom->gcon[GIND(3,3)];

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];


  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamax=GAMMAMAXRAD;
    FTYPE gammamaxfail=GAMMAMAXRADFAIL;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);


  FTYPE Avcovorig[NDIM],Avconorig[NDIM];
  DLOOPA(jj) Avcovorig[jj]=Avcov[jj];
  DLOOPA(jj) Avconorig[jj]=Avcon[jj];

  // get new Avcov[0]=R^t_t

  // NOTEMARK: Note that Sqrt() is only ever negative when gammarel2<0, so never has to be concern.
  Avcov[0]=(0.25*(-4.*(gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*utsq*(gn11 + utsq) + 
                        Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                              1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                       gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2))))/(gn11*utsq*(gn11 + utsq));

  Erf=(0.75*Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                           1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                    gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2)))/(utsq*(gn11 + utsq)*(gn11 + 4.*utsq));

  if(0&&showmessages && debugfail>=2) dualfprintf(fail_file,"NOR SOL: Avcov0new=%g Avcov0old=%g Erf=%g :: %g %g %g\n",Avcov[0],Avcovorig[0],Erf,Rtx,Rty,Rtz);

  //modify Avcon
  indices_12(Avcov,Avcon,ptrgeom);
  if(0&&showmessages && debugfail>=2) DLOOPA(jj) dualfprintf(fail_file,"jj=%d Avconorig=%g Avconnew=%g\n",jj,Avconorig[jj],Avcon[jj]);

 
  delta=0; // not yet

  if(0&&showmessages && debugfail>=2){
    // alt solution

    FTYPE Avcovalt = (-0.25*(4.*(gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*utsq*(gn11 + utsq) + 
       Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
           1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                    gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2))))/(gn11*utsq*(gn11 + utsq));
    
    FTYPE Erfalt=(-0.75*Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                        1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                  gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2)))/(utsq*(gn11 + utsq)*(gn11 + 4.*utsq));

    if(showmessages && debugfail>=2) dualfprintf(fail_file,"ALT SOL: Avcovalt=%g Avcov0old=%g Erfalt=%g : %g %g %g\n",Avcovalt,Avcovorig[0],Erfalt,Rdtx,Rdty,Rdtz);
  }


  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*ptrgeom->gcon[GIND(0,jj)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  if(0&&showmessages && debugfail>=2) dualfprintf(fail_file,"NORM ROOT 4-vel: %g %g %g : %g\n",urfconrel[1],urfconrel[2],urfconrel[3],ptrgeom->gdet);

  
  *Erfreturn=Erf; // pass back new Erf to pointer


  return(0);
}




































//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
int calc_tautot(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *tautot, FTYPE *tautotmax)
{
  //xx[0] holds time
  FTYPE kappa,kappaes,chi;
  calc_kappa(pp,ptrgeom,&kappa);
  calc_kappaes(pp,ptrgeom,&kappaes);
  chi=kappa+kappaes;
  int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1}; // want to ignore non-used dimensions

  int jj;
  *tautotmax=0.0;
  SLOOPA(jj){
    tautot[jj]=chi * (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])))*NxNOT1[jj];
    *tautotmax=MAX(*tautotmax,tautot[jj]);
  }

  return 0;
}

//**********************************************************************
//******* calculates abs opacity over dx[] ***************************
//**********************************************************************
int calc_tauabs(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *tauabs, FTYPE *tauabsmax)
{
  FTYPE kappa;
  calc_kappa(pp,ptrgeom,&kappa);

  int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1}; // want to ignore non-used dimensions

  int jj;
  *tauabsmax=0.0;
  SLOOPA(jj){
    tauabs[jj]=kappa * (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])))*NxNOT1[jj];
    *tauabsmax=MAX(*tauabsmax,tauabs[jj]);
  }

  return 0;
}






//**********************************************************************
//suplementary routines for conversions
//**********************************************************************
FTYPE calc_PEQ_ufromTrho(FTYPE T,FTYPE rho)
{
  // if use local function function instead of below directly,
  // then assume user doesn't care about position for EOS.
  FTYPE u=u_rho0_T_simple(0, 0, 0, CENT, rho, T);
  return u;
}

FTYPE calc_PEQ_Tfromurho(FTYPE u,FTYPE rho)
{
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
  return T;
}

// E=urad=arad T^4 (this is LTE only if put in T was gas T)
FTYPE calc_LTE_EfromT(FTYPE T)
{
  //  return 4.*SIGMA_RAD*T*T*T*T;
  return (ARAD_CODE*T*T*T*T);
}

// E=urad=arad T^4 and just solve for T  (this is LTE only if assume resulting T is gas T).  If put in fluid-frame E, then correct T for radiation in fluid frame.
FTYPE calc_LTE_TfromE(FTYPE E )
{
  //  return sqrt(sqrt((E/4./SIGMA_RAD)));
  return (sqrt(sqrt((E/(SMALL+ARAD_CODE)))));
}

// This will really give back only LTE E
FTYPE calc_LTE_Efromurho(FTYPE u,FTYPE rho)
{
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
  return (calc_LTE_EfromT(T));
}



// set velocity based upon ncon and gammamax and return in whichvel format for the ptrgeom geometry/coords
int set_ncon_velocity(int whichvel, FTYPE gammamax, FTYPE *ncon, struct of_geom *ptrgeom, FTYPE *uconwhichvel)
{
  int ii,jj;
  FTYPE ncondefault[NDIM]={0.0,-1.0, 0.0, 0.0}; //for radially flowing photons
  FTYPE *nconuse;

  // default is radial motion in zamo frame
  if(ncon==NULL) nconuse=ncondefault;
  else nconuse=ncon;

  // compute \gammarel using nconuse
  FTYPE gammatemp,qsq;
  gamma_calc_fromuconrel(nconuse,ptrgeom,&gammatemp,&qsq);

  // now rescale nconuse[i] so will give desired \gammamax
  FTYPE prtemp[NPR];
  SLOOPA(ii) prtemp[U1+ii-1] = nconuse[ii]*(gammamax/gammatemp);

  // now get u^\mu[lab]
  FTYPE uconlab[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  ucon_calc_whichvel(WHICHVEL,prtemp,ptrgeom,uconlab,others);

  // now get other whichvel type
  FTYPE prtemp2[NPR];
  ucon2pr(whichvel,uconlab,ptrgeom,prtemp2);

  SLOOPA(jj) uconwhichvel[jj] = prtemp2[U1+jj-1];

  return(0);

}



//#include "phys.tools.rad.notused.c"

