// KORALTODO: B Not setup for iter=0 2 1 0 : 1 2 : 3 13 : 14 200 : 14

// include globals
#include "decs.h"

//////////////////////////////
//
// BEGIN: f2c stuff for Ramesh's solver code in fortran
//
//////////////////////////////

// f2c prototype
#include "f2c.h"

#ifdef KR_headers
double d_sign(a,b) doublereal *a, *b;
#else
double d_sign(doublereal *a, doublereal *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}

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


#ifdef KR_headers
extern void f_exit();
int s_stop(s, n) char *s; ftnlen n;
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
extern "C" {
#endif
void f_exit(void);

int s_stop(char *s, ftnlen n)
#endif
{
int i;

if(n > 0)
	{
	fprintf(stderr, "STOP ");
	for(i = 0; i<n ; ++i)
		putc(*s++, stderr);
	fprintf(stderr, " statement executed\n");
	}
#ifdef NO_ONEXIT
f_exit();
#endif
exit(0);

/* We cannot avoid (useless) compiler diagnostics here:		*/
/* some compilers complain if there is no return statement,	*/
/* and others complain that this one cannot be reached.		*/

return 0; /* NOT REACHED */
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
}
#endif


#include "testfpp.P"
// not linking with libf2c since don't want that dependence and conversion doesn't need it since the original code was simple

int get_rameshsolution(int whichcall, int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, int baseitermethod, FTYPE *errorabs, FTYPE *errorabsbestexternal, int iters, int totaliters, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent);

int get_rameshsolution_wrapper(int whichcall, int eomtype, FTYPE *errorabs, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *dUother, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, FTYPE (*dUcompeng)[NPR], FTYPE (*dUcompent)[NPR], struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent);

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

// Funny, even 1E-5 does ok with torus, no worse at Erf~ERADLIMIT instances.  Also, does ~3 iterations, but not any faster than using 1E-12 with ~6 iterations.
#define IMPTRYCONV (1.e-12) // works generally to avoid high iterations
#define IMPTRYCONVQUICK (1.e-9) // less greedy so doesn't slow things down so much.
// error for comparing to sum over all absolute errors
#define IMPTRYCONVABS ((FTYPE)(NDIM+2)*trueimptryconv)

// what tolerance to use for saying can switch to entropy when u_g is suggested to be bad for energy
#define IMPOKCONVCONST (1E-9)
#define IMPOKCONVCONSTABS ((FTYPE)(NDIM+2)*IMPOKCONVCONST)
#define IMPOKCONV (MAX(trueimptryconv,IMPOKCONVCONST))
#define IMPOKCONVABS ((FTYPE)(NDIM+2)*IMPOKCONV)

// too allowing to allow 1E-4 error since often solution is nuts at even errors>1E-8
#define IMPALLOWCONVCONST (1.e-7)
#define IMPALLOWCONVCONSTABS ((FTYPE)(NDIM+2)*IMPALLOWCONVCONST)
//#define IMPALLOWCONV (MAX(trueimptryconv,IMPALLOWCONVCONST))
#define IMPALLOWCONV (trueimpallowconv)
#define IMPALLOWCONV2 (IMPALLOWCONV)
#define IMPALLOWCONVABS ((FTYPE)(NDIM+2)*IMPALLOWCONV)

// tolerance above which say energy solution is probably bad even if not very large error.  These have tended (or nearly 100%) to be cases where actual solution has u_g<0 but harm gets error u_g>0 and error not too large.
#define IMPBADENERGY (MIN(IMPALLOWCONV,1E-7))

// how many iterations before we try harder to get better 1D MHD inversion solution
#define ITERMHDINVTRYHARDER 5
#define MINTRYCONVFORMHDINVERSION (1E-4) // assume not failure if got down to this much. -- don't have to be related to implicit allowance.


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
#define IMPEPSLARGE (1E-8) // Was 1E-6 (see next line)
#define IMPEPSSMALL (1E-10) // Was 1E-7 (but with PMHD method and disk-jet problem, led to very often not getting solution.  Not sure why)
//#define IMPEPSSMALL (1E-8)
#define ERRORFORIMPEPSSMALL (1E-9)
#elif(REALTYPE==LONGDOUBLETYPE)
#define IMPEPSLARGE (1E-8)
#define IMPEPSSMALL (1E-10)
#define ERRORFORIMPEPSSMALL (1E-9)
#endif

// maximum EPS for getting Jacobian
#define MAXIMPEPS (0.3)

// maximum number of times to (typically) increase EPS in getting Jacobian for implicit scheme.  This might generally override MAXIMPEPS.
#define MAXJACITER (10)



//#define IMPMAXITER (15) // for used implicit solver // For others
#define IMPMAXITERLONG (100) // for used implicit solver
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


#define MAXF1TRIES 20 // 20 might sound like alot, but Jacobian does 4 to 5 inversions each iteration, and that amount is only typically needed for very first iteration.
  // goes to f1iter=10 for RADPULSE KAPPAES=1E3 case.  Might want to scale maximum iterations with \tau, although doubling of damping means exponential w.r.t. f1iter, so probably 50 is certainly enough since 2^(-50) is beyond machine precision for doubles.

#define MAXMHDSTEPS (trueimpmaxiter*6) // limit total number of MHD inversion steps


#define RADDAMPDELTA (0.5) // choose, but best to choose 1/Integer so no machine precision issues.
#define RADDAMPUNDELTA (1.0/RADDAMPDELTA)


#define TAUFAILLIMIT (2.0/3.0) // at what \tau below which to assume "failure1" in u2p_rad() means should be moving at gammamax rather than not moving.

#define TAUSWITCHPBVSPIIN (2.0/3.0) // at what \tau to switch using pb(uu0) vs. piin(Uiin).


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

// determine which error will use when deciding if converged or not.
// If only use iterated, then rest can be large error, and that's not desired.  So generally should use WHICHERROR 1
#define NUMERRORTYPES 2 // 0: over iterated   1: over all relevant terms
#define WHICHERROR 1 // choose, but generally should be 1.


//SPACETIMESUBSPLITTIMEMHDRAD // TAUSUPPRESS
// SPACETIMESUBSPLITTIMEMHDRAD  //  SPACETIMESUBSPLITMHDRAD // SPACETIMESUBSPLITSUPERALL // TAUSUPPRESS

// whether to get pb(uu0) instead of using pb that was used to compute F(pb).  Better guess to use pb(uu0) in optically thin regions.
// : 0 don't
// : 1 do for RAD and MHD methods under all cases
// : 2 do only for RAD methods
// : 3 do for both methods but only for STAGES for MHD methods
//#define GETRADINVFROMUU0FORPB 3 // makes use more ITERMODESTAGES
#define GETRADINVFROMUU0FORPB 1  // might be expensive of uu0 has no solution for the MHD inversion.

// whether to use dUriemann and dUgeom or other dU's sitting in dUother for radiation update
// -1 : use Uiin,piin
// 0: use Uiin,piin [initial U and initial p] // bad choice to use, e.g., uu0[RHO] should always be used
// 1: use uu0 for non-iterated U's and pb for non-iterated.  And Uiin's for iterated U's and piin for iterated p's. // balanced choice
// 2: use uu0,pb [initial+flux U and pb used to get flux] // maybe best for dynamical flows
// 3: use uu0,puu0 [initial+flux U and p] // costly since needs inversion, and probably no better than 1
#define USEDUINRADUPDATE 1

// whether to use inputted uub, and pb as guess if errorabsreturn inputted is small enough
// makes sense in general only if WHICHERROR=1
#define USEINPUTASGUESSIFERRORSMALL (WHICHERROR==1)


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
#define DOSUBJAC 1
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
#define ENDFULLSTEPS0 (IMPMAXITERLONG*2)

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
#define JONHOLDPOS 1

#define NEWJONHOLDPOS 0 // FUCK: WORKING ON IT

// number of times allowed to hold u_g as positive
#define NUMHOLDTIMES 6


// stop iterating energy if pbenergy[UU]<0.5*pbentropy[UU] consistently starting aafter below number of iterations and lasting for 2nd below number of iterations
#define RAMESHSTOPENERGYIFTOOOFTENBELOWENTROPY 3

//#define SWITCHTOENTROPYIFCHANGESTOENTROPY (*implicitferr==QTYUMHD ? 0 : 1)
#define SWITCHTOENTROPYIFCHANGESTOENTROPY (0)



// whether to use ramesh solver as backup
#define USERAMESH 0 // too slow if used too often, and rarely result really used.

// error below which to feed best guess into next attempt
#define TRYHARDERFEEDGUESSTOL (1E-4)

// error below which will use entropy as guess for energy if entropy didn't hard fail.
#define ERRORTOUSEENTROPYFORENERGYGUESS (1E-4)


// whether to get lowest error solution instead of final one.
#define GETBEST 1


// need to do final check since get f1 and then do step.
#define DOFINALCHECK 1


// below 1 if reporting cases when MAXITER reached, but allowd error so not otherwise normally reported.
#define REPORTMAXITERALLOWED (PRODUCTION==0)

// whether to ensure rho and u_g in Jacobian calculation difference do not cross over 0 and stay on same side as origin point.
#define FORCEJDIFFNOCROSS 1

// whether to check pp-ppp
// 1: directly check post pp-ppp relative error and see if changes by LOCALPREIMPCONVX
// 2: directly check if any changes to pp during Newton step are bigger than DIFFXLIMIT.
#define POSTNEWTONCONVCHECK 1

// below which sum of all primitives is taken as no interesting change.
#define DIFFXLIMIT (10.0*NUMEPSILON)

#define LOCALPREIMPCONVX (10.0*NUMEPSILON)
#define LOCALPREIMPCONVXABS ((FTYPE)(NDIM+2)*LOCALPREIMPCONVX)


// number of iterations by which to check (1st) whether after some number of times (2nd) error rose instead of reduced.
#define NUMNOERRORREDUCE0 (5+BEGINNORMALSTEPS)
#define NUMNOERRORREDUCE 5


// whether to use EOMDONOTHING if error is good enough.
// 1: always avoid external inversion (so no longer can do cold MHD, but cold MHD in \tau\gtrsim 1 places is very bad).  Or avoid energy switching to entropy, which also is bad.
#define SWITCHTODONOTHING 1

  // whether to change damp factor during this instance.
#define CHANGEDAMPFACTOR 2 // bit risky to set to 1 since changing DAMPFACTOR for no good reason limits ability to converge at normal rate.
#define NUMDAMPATTEMPTS 3

#define NUMDAMPATTEMPTSQUICK 1


// factor by which error jumps as indication that u_g stepped to was very bad choice.
#define FACTORBADJUMPERROR (1.0E2)


// 0 : old Jon  method
// 1 : Jon's paper draft method
#define WHICHU2PRAD 1

// during implicit solver, don't limit gamma so much as normally.  Otherwise, solution may not be found and solver struggles and leads to high errors and iterations.  If limit gammarad but not gammafluid, then gammafluid can be too high.  If limit both, no solutions can be found.   So just limit afterwards for now.
//#define GAMMAMAXRADIMPLICITSOLVER (1E5)
#define GAMMAMAXRADIMPLICITSOLVER (GAMMAMAXRAD) // for radiation, seek actual solution with this limit.  Solver can find solutions, while harder when limiting gamma_{gas} for some reason.


// whether to avoid computing entropy during iterations if not needed
#define ENTROPYOPT 1

// whether to try using uualt (see f_implicit())
#define ALLOWUSEUUALT 0

// whether to use CAPTYPEFIX2 for f1 (central Jac and error estimate).  Still will use CAPTYPEBASIC for final check, but this allows faster convergence or non-convergence as required.
#define USECAPTYPEFIX2FORF1 1
// "" for final check.
#define USECAPTYPEFIX2FORFINALCHECK 1
// whether to avoid including URAD0 in total error when hitting radinvmod!=0.  If using CAPTYPEFIX2FORF1=1, then can't get error better than what CAPTYPEFIX2 provides for non-URAD0 terms.
#define AVOIDURAD0IFRADINVMODANDPMHDMETHOD (USECAPTYPEFIX2FORF1!=0)
// whether to avoid back-tracing f1 calculation if rad inv hits cap, just push through if ==1.
#define TREATRADINVCAPASNONFAILUREFORPMHDMETHOD (USECAPTYPEFIX2FORF1!=0)

// whether to just let PMHD fail and try it first no matter whether primary considerations say otherwise.
#define LETPMHDFAIL 1

// whether to use history of methods to determine current method -- KORALTODO: Needs work.
#define USEPRIORITERMETHOD 0

// stop if WHICHERROR is low error
#define STOPIFVERYLOWERROR (1)
// stop if iter error is low error
#define STOPIFITERLOWERROR (0)

///////////////////////////////
//
// END SOME LOCAL OPTIONS
//
///////////////////////////////






//////// implicit stuff
static int koral_source_rad_implicit(int *eomtype, FTYPE *pb, FTYPE *pf, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother ,FTYPE (*dUcomp)[NPR]);

static int koral_source_rad_implicit_mode(int allowbaseitermethodswitch, int modprim, int havebackup, int didentropyalready, int *eomtype, int whichcap, int itermode, int *baseitermethod, FTYPE trueimptryconv, FTYPE trueimpokconv, FTYPE trueimpallowconv, int trueimpmaxiter, int truenumdampattempts, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pb, FTYPE *uub, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR], FTYPE *errorabs, FTYPE *errorabsbestexternal, int *iters, int *f1iters, int *nummhdinvs, int *nummhdsteps);

static int f_implicit(int allowbaseitermethodswitch, int iter, int f1iter, int failreturnallowable, int whichcall, FTYPE impeps, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, int *baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE conv, FTYPE convabs, FTYPE allowconvabs, int maxiter, FTYPE realdt, int dimtypef, FTYPE *dimfactU, FTYPE *ppprev, FTYPE *pp, FTYPE *piin, FTYPE *uuprev, FTYPE *Uiin, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom, struct of_state *q,  FTYPE *f, FTYPE *fnorm, FTYPE *freport, int *goexplicit, FTYPE *errorabs, FTYPE *errorallabs, int whicherror, int *convreturn, int *nummhdinvsreturn);


static int calc_tautot_chieff(FTYPE *pp, FTYPE chieff, struct of_geom *ptrgeom, FTYPE *tautot, FTYPE *tautotmax);


static int get_implicit_iJ(int allowbaseitermethodswitch, int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtypelocal, int whichcap, int itermode, int *baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, FTYPE impepsjac, FTYPE trueimptryconv, FTYPE trueimptryconvabs, FTYPE trueimpallowconvabs, int trueimpmaxiter, int iter, FTYPE errorabs, FTYPE errorallabs, int whicherror, int dimtypef, FTYPE *dimfactU, FTYPE *Uiin, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *piin, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NPR], int *nummhdinvsreturn);

static int inverse_33matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM]);
static int inverse_11matrix(int sj, int ej, FTYPE a[][NDIM], FTYPE ia[][NDIM]);


static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE convabs, FTYPE realdt, int dimtypef, int eomtype, int radinvmod, int itermode, int baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, FTYPE *dimfactU, FTYPE *pp, FTYPE *piin, FTYPE *f1, FTYPE *f1norm, FTYPE *f1report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom, FTYPE *errorabs, FTYPE *errorallabs, int whicherror);

static int compute_ZAMORAD(FTYPE *uu, struct of_geom *ptrgeom, FTYPE *Er, FTYPE *Utildesq, FTYPE *Utildecon);


static int Utoprimgen_failwrapper(int doradonly, int *radinvmod, int showmessages, int checkoninversiongas, int checkoninversionrad, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int whichcap, int evolvetype, int inputtype,FTYPE *U,  struct of_state *qptr, struct of_geom *ptrgeom, FTYPE dissmeasure, FTYPE *pr, struct of_newtonstats *newtonstats);

static void define_method(int iter, int *eomtype, int itermode, int baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, int *implicititer, int *implicitferr, int *BEGINMOMSTEPS, int *ENDMOMSTEPS, int *BEGINENERGYSTEPS, int *ENDENERGYSTEPS, int *BEGINFULLSTEPS, int *ENDFULLSTEPS, int *BEGINNORMALSTEPS);
static void get_refUs(int *numdims, int *startjac, int *endjac, int *implicititer, int *implicitferr, int *irefU, int *iotherU, int *erefU, int *eotherU, int *signgd2, int *signgd4, int *signgd6, int *signgd7);


// debug stuff
static void showdebuglist(int debugiter, FTYPE (*pppreholdlist)[NPR],FTYPE (*ppposholdlist)[NPR],FTYPE (*f1reportlist)[NPR],FTYPE (*f1list)[NPR],FTYPE *errorabsf1list,FTYPE *errorallabsf1list, int *realiterlist, FTYPE (*jaclist)[NDIM][NDIM], FTYPE *fracdamplist, int *implicititerlist, int *implicitferrlist);
int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, int baseitermethod, FTYPE *errorabs, FTYPE *errorabsbestexternal, int iters, int iterstotal, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);

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


// f_implicit() call types
#define FIMPLICITCALLTYPEF1 1
#define FIMPLICITCALLTYPEFINALCHECK 2
#define FIMPLICITCALLTYPEJAC 3
#define FIMPLICITCALLTYPEFINALCHECK2 4





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
static int Utoprimgen_failwrapper(int doradonly, int *radinvmod, int showmessages, int checkoninversiongas, int checkoninversionrad, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int whichcap, int evolvetype, int inputtype,FTYPE *U,  struct of_state *qptr,  struct of_geom *ptrgeom, FTYPE dissmeasure, FTYPE *pr, struct of_newtonstats *newtonstats)
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
    int modprim=0;
    //    int checkoninversiongas=CHECKONINVERSION;
    //    int checkoninversionrad=CHECKONINVERSIONRAD;
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad, allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, whichmethod, modprim, evolvetype, inputtype, U, qptr, ptrgeom, dissmeasure, pr, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);
    *radinvmod=(int)(*lpflagrad);
    //    *radinvmod=0;  //KORALTODO: Not using method that needs this call, so for now don't pass radinvmod through to Utoprimgen().
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

// MHD or RAD baseitermethod types
#define IMPMHDTYPEBASE(baseitermethod) (baseitermethod==QTYPMHD || baseitermethod==QTYUMHD || baseitermethod==QTYENTROPYUMHD)

#define IMPRADTYPEBASE(baseitermethod) (baseitermethod==QTYPRAD || baseitermethod==QTYURAD)

// PMHD types
#define IMPPMHDTYPE(implicititer) (implicititer==QTYPMHD ||implicititer==QTYPMHDENERGYONLY ||implicititer==QTYPMHDMOMONLY)



// set default method
// can control method per iteration
static void define_method(int iter, int *eomtype, int itermode, int baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, int *implicititer, int *implicitferr, int *BEGINMOMSTEPS, int *ENDMOMSTEPS, int *BEGINENERGYSTEPS, int *ENDENERGYSTEPS, int *BEGINFULLSTEPS, int *ENDFULLSTEPS, int *BEGINNORMALSTEPS)
{
  int eomtypelocal=*eomtype; // default, but not changing it so far.


  if(eomtypelocal==EOMDEFAULT){
    eomtypelocal=EOMTYPE; // override
  }


  if(PRODUCTION==0&&EOMDONOTHING(eomtypelocal)){
    dualfprintf(fail_file,"Can't have EOMDONOTHING in radiation code.\n");
    myexit(938463651);
  }

  if(0){
    if(fracenergy>0.0 && fracenergy<1.0){
      // override as "energy" method, so erefU[TT]=UU to keep things simple and avoid confusions as to where error function is put and avoid duplicates.
      if(EOMTYPE==EOMGRMHD) eomtypelocal=EOMTYPE;
      else if(PRODUCTION==0){
        dualfprintf(fail_file,"Shouldn't be in define_method with fracenergy=%g and EOMTYPE=%d\n",fracenergy,EOMTYPE);
        myexit(2875346346);
      }
    }
  }


  if(itermode==ITERMODESTAGES){
    *BEGINMOMSTEPS=BEGINMOMSTEPS0;
    *ENDMOMSTEPS=ENDMOMSTEPS0;

    *BEGINENERGYSTEPS=BEGINENERGYSTEPS0;
    *ENDENERGYSTEPS=ENDENERGYSTEPS0;

    *BEGINFULLSTEPS=BEGINFULLSTEPS0;
    *ENDFULLSTEPS=ENDFULLSTEPS0;

    *BEGINNORMALSTEPS=BEGINNORMALSTEPS0;
  }
  else if(itermode==ITERMODENORMAL){
    *BEGINMOMSTEPS=-1;
    *ENDMOMSTEPS=-1;

    *BEGINENERGYSTEPS=-1;
    *ENDENERGYSTEPS=-1;

    *BEGINFULLSTEPS=1;
    *ENDFULLSTEPS=(IMPMAXITERLONG*2);

    *BEGINNORMALSTEPS=1; // has to be consistent with actual first iteration, which is currently 1.
  }
  else if(itermode==ITERMODECOLD){
    *BEGINMOMSTEPS=1;
    *ENDMOMSTEPS=IMPMAXITERLONG*2;

    *BEGINENERGYSTEPS=-1;
    *ENDENERGYSTEPS=-1;

    *BEGINFULLSTEPS=-1;
    *ENDFULLSTEPS=-1;

    *BEGINNORMALSTEPS=1;
  }
  else if(PRODUCTION==0){
    dualfprintf(fail_file,"No such itermode=%d\n",itermode);
    myexit(30486346);
  }

  /// PMHD:
  if(eomtypelocal==EOMGRMHD||eomtypelocal==EOMCOLDGRMHD){

    /// OLD SETUP, works even with entropy source
    // V1: Works fine, but slow.  Gives hot fluid near jumps like torus edge.  Computing UU[ENTROPY] with or without source doesn't change much.
    //*implicititer=(QTYURAD); // choice
    //*implicitferr=(QTYURAD); // choice

    /// UMHD: no setup works unless entropy source turned off, and then QTYENTROPYUMHD as FERR can't be chosen.
    // V1: Ok unless entropy source added.
    //*implicititer=(QTYUMHD); // choice
    //*implicitferr=(QTYUMHD); // choice

    /// PRAD:
    // V1: unstable with entropy source
    //*implicititer=(QTYPRAD);
    //*implicitferr=(QTYURAD);


    if(DOSUBJAC==1){
      // assume iter=1 is first iteration
      if(iter>=*BEGINMOMSTEPS && iter<=*ENDMOMSTEPS){
        if(baseitermethod==QTYPMHD){
          *implicititer=(QTYPMHDMOMONLY); // choice
          *implicitferr=(QTYUMHDMOMONLY);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURADMOMONLY); // choice
          *implicitferr=(QTYURADMOMONLY); // choice
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRADMOMONLY); // choice
          *implicitferr=(QTYURADMOMONLY); // choice
        }
      }
      else if(iter>=*BEGINENERGYSTEPS && iter<=*ENDENERGYSTEPS){
        if(baseitermethod==QTYPMHD){
          *implicititer=(QTYPMHDENERGYONLY); // choice
          *implicitferr=(QTYUMHDENERGYONLY);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURADENERGYONLY); // choice
          *implicitferr=(QTYURADENERGYONLY);
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRADENERGYONLY); // choice
          *implicitferr=(QTYURADENERGYONLY);
        }
      }
      else if(iter>=*BEGINFULLSTEPS && iter<=*ENDFULLSTEPS){
        if(baseitermethod==QTYPMHD){
          // V1: noisy in atmosphere and torus edge, but otherwise ok.
          *implicititer=(QTYPMHD); // choice
          *implicitferr=(QTYUMHD);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURAD); // choice
          *implicitferr=(QTYURAD); // choice
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRAD); // choice
          *implicitferr=(QTYURAD); // choice
        }
      }
      else if(PRODUCTION==0){
        dualfprintf(fail_file,"A Not setup for iter=%d %d %d %g : %d %d : %d %d : %d %d : %d : %d\n",iter, *eomtype, itermode, fracenergy, *BEGINMOMSTEPS, *ENDMOMSTEPS, *BEGINENERGYSTEPS, *ENDENERGYSTEPS, *BEGINFULLSTEPS, *ENDFULLSTEPS, *BEGINNORMALSTEPS,baseitermethod);
        myexit(3498346457);
      }
    }
    else{
      if(baseitermethod==QTYPMHD){
        *implicititer=(QTYPMHD); // choice
        *implicitferr=(QTYUMHD);
      }
      else if(baseitermethod==QTYURAD){
        *implicititer=(QTYURAD); // choice
        *implicitferr=(QTYURAD); // choice
      }
      else if(baseitermethod==QTYPRAD){
        *implicititer=(QTYPRAD); // choice
        *implicitferr=(QTYURAD); // choice
      }
    }
  }
  else if(eomtypelocal==EOMENTROPYGRMHD){

    /// UMHD: no setup works unless entropy source turned off, and then QTYENTROPYUMHD as FERR can't be chosen.
    // V2: Like V1 but goes unstable.
    //*implicititer=(QTYUMHD); // choice
    //*implicitferr=(QTYENTROPYUMHD); // choice


    if(DOSUBJAC==1){
      // assume iter=1 is first iteration
      if(iter>=*BEGINMOMSTEPS && iter<=*ENDMOMSTEPS){
        if(baseitermethod==QTYPMHD){
          *implicititer=(QTYPMHDMOMONLY); // choice
          *implicitferr=(QTYENTROPYUMHDMOMONLY);
        }
        else if(baseitermethod==QTYENTROPYUMHD){
          *implicititer=(QTYENTROPYUMHDMOMONLY);
          *implicitferr=(QTYENTROPYUMHDMOMONLY);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURADMOMONLY); // choice
          *implicitferr=(QTYURADMOMONLY); // choice
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRADMOMONLY); // choice
          *implicitferr=(QTYURADMOMONLY); // choice
        }
      }
      else if(iter>=*BEGINENERGYSTEPS && iter<=*ENDENERGYSTEPS){
        if(baseitermethod==QTYPMHD){
          *implicititer=(QTYPMHDENERGYONLY); // choice
          *implicitferr=(QTYENTROPYUMHDENERGYONLY);
        }
        else if(baseitermethod==QTYENTROPYUMHD){
          *implicititer=(QTYENTROPYUMHDENERGYONLY);
          *implicitferr=(QTYENTROPYUMHDENERGYONLY);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURADENERGYONLY); // choice
          *implicitferr=(QTYURADENERGYONLY); // choice
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRADENERGYONLY); // choice
          *implicitferr=(QTYURADENERGYONLY);
        }
      }
      else if(iter>=*BEGINFULLSTEPS && iter<=*ENDFULLSTEPS){
        if(baseitermethod==QTYPMHD){
          // V2: perfectly fine as entropy method
          *implicititer=(QTYPMHD); // choice
          *implicitferr=(QTYENTROPYUMHD); // choice
        }
        else if(baseitermethod==QTYENTROPYUMHD){
          // V1: works perfectly fine and acts like PMHD,ENTROPYUMHD method just slower:
          *implicititer=(QTYENTROPYUMHD);
          *implicitferr=(QTYENTROPYUMHD);
        }
        else if(baseitermethod==QTYURAD){
          *implicititer=(QTYURAD); // choice
          *implicitferr=(QTYURAD); // choice
        }
        else if(baseitermethod==QTYPRAD){
          *implicititer=(QTYPRAD); // choice
          *implicitferr=(QTYURAD);
        }
      }
      else if(PRODUCTION==0){
        dualfprintf(fail_file,"B Not setup for iter=%d %d %d %g : %d %d : %d %d : %d %d : %d : %d\n",iter, *eomtype, itermode, fracenergy, *BEGINMOMSTEPS, *ENDMOMSTEPS, *BEGINENERGYSTEPS, *ENDENERGYSTEPS, *BEGINFULLSTEPS, *ENDFULLSTEPS, *BEGINNORMALSTEPS,baseitermethod);
        myexit(3498346458);
      }
    }
    else{
      if(baseitermethod==QTYPMHD){
        *implicititer=(QTYPMHD); // choice
        *implicitferr=(QTYENTROPYUMHD); // choice
      }
      else if(baseitermethod==QTYENTROPYUMHD){
        *implicititer=(QTYENTROPYUMHD);
        *implicitferr=(QTYENTROPYUMHD);
      }
      else if(baseitermethod==QTYURAD){
        *implicititer=(QTYURAD); // choice
        *implicitferr=(QTYURAD); // choice
      }
      else if(baseitermethod==QTYPRAD){
        *implicititer=(QTYPRAD); // choice
        *implicitferr=(QTYURAD);
      }
    }
  }
  else if(PRODUCTION==0){
    dualfprintf(fail_file,"No such eomtypelocal=%d in define_method\n",eomtypelocal);
    myexit(938463653);
  }

  // TESTING:
  //*implicititer=(QTYURAD); // choice
  //*implicitferr=(QTYURAD); // choice


}



#define JACLOOP(jj,startjj,endjj) for(jj=startjj;jj<=endjj;jj++)
#define JACLOOPALT(jj,startjj,endjj) DLOOPA(jj) //for(jj=startjj;jj<=endjj;jj++) // for those things might or might not want to do all terms
#define JACLOOPSUPERFULL(pliter,pl,eomtype,baseitermethod,radinvmod) PLOOP(pliter,pl) if(pl!=ENTROPY && pl!=UU && pl!=URAD0 || (eomtype==EOMDEFAULT && EOMTYPE==EOMENTROPYGRMHD || eomtype==EOMENTROPYGRMHD || eomtype==EOMDIDENTROPYGRMHD) && (pl==ENTROPY || pl==URAD0 && IMPMHDTYPEBASE(baseitermethod)==0 || IMPMHDTYPEBASE(baseitermethod)==1 && (pl==URAD0 && AVOIDURAD0IFRADINVMODANDPMHDMETHOD==0 || pl==URAD0 && AVOIDURAD0IFRADINVMODANDPMHDMETHOD==1 && radinvmod==0)) || (eomtype==EOMDEFAULT && EOMTYPE==EOMGRMHD || eomtype==EOMGRMHD || eomtype==EOMDIDGRMHD) && (pl==UU || pl==URAD0 && IMPMHDTYPEBASE(baseitermethod)==0 || IMPMHDTYPEBASE(baseitermethod)==1 && (pl==URAD0 && AVOIDURAD0IFRADINVMODANDPMHDMETHOD==0 || pl==URAD0 && AVOIDURAD0IFRADINVMODANDPMHDMETHOD==1 && radinvmod==0) )) // over f's, not primitives.
#define JACLOOPFULLERROR(itermode,jj,startjj,endjj) for(jj=(itermode==ITERMODECOLD ? startjj : 0);jj<=(itermode==ITERMODECOLD ? endjj : NDIM-1);jj++)
#define JACLOOPSUBERROR(jj,startjj,endjj) JACLOOP(jj,startjj,endjj)
#define JACLOOP2D(ii,jj,startjj,endjj) JACLOOP(ii,startjj,endjj) JACLOOP(jj,startjj,endjj)

static void get_refUs(int *numdims, int *startjac, int *endjac, int *implicititer, int *implicitferr, int *irefU, int *iotherU, int *erefU, int *eotherU, int *signgd2, int *signgd4, int *signgd6, int *signgd7)
{
  int jj;

  // default
  *numdims=NDIM;
  *signgd7= (+1.0); // not used for PMHD
  DLOOPA(jj) irefU[jj]=UU+jj;
  DLOOPA(jj) iotherU[jj]=URAD0+jj;
  *startjac=0; *endjac=NDIM-1;

  // same list and numbers in array for both primitives and conserved
  if(*implicititer==QTYUMHD || *implicititer==QTYPMHD){
  }
  else if(*implicititer==QTYUMHDENERGYONLY || *implicititer==QTYPMHDENERGYONLY){
    *startjac=0; *endjac=0;
  }
  else if(*implicititer==QTYUMHDMOMONLY || *implicititer==QTYPMHDMOMONLY){
    *startjac=1; *endjac=3;
  }
  else if(*implicititer==QTYENTROPYUMHD || *implicititer==QTYENTROPYPMHD){
    irefU[TT]=ENTROPY;
  }
  else if(*implicititer==QTYENTROPYUMHDENERGYONLY){
    irefU[TT]=ENTROPY;
    *startjac=0; *endjac=0;
  }
  else if(*implicititer==QTYENTROPYUMHDMOMONLY){
    irefU[TT]=ENTROPY;
    *startjac=1; *endjac=3;
  }
  else if(*implicititer==QTYURAD || *implicititer==QTYPRAD){
    *numdims=NDIM;
    *signgd7= (+1.0); // required to make URAD method work for (e.g.) RADSHADOW if using Gddt-based GS
    DLOOPA(jj) irefU[jj]=URAD0+jj;
    DLOOPA(jj) iotherU[jj]=UU+jj;
    *startjac=0; *endjac=NDIM-1;
  }
  else if(*implicititer==QTYURADENERGYONLY || *implicititer==QTYPRADENERGYONLY){
    *numdims=NDIM;
    *signgd7= (+1.0);
    DLOOPA(jj) irefU[jj]=URAD0+jj;
    DLOOPA(jj) iotherU[jj]=UU+jj;
    *startjac=0; *endjac=0;
  }
  else if(*implicititer==QTYURADMOMONLY || *implicititer==QTYPRADMOMONLY){
    *numdims=NDIM;
    *signgd7= (+1.0);
    DLOOPA(jj) irefU[jj]=URAD0+jj;
    DLOOPA(jj) iotherU[jj]=UU+jj;
    *startjac=1; *endjac=3;
  }
  else if(PRODUCTION==0){
    dualfprintf(fail_file,"No such implicititer=%d\n",*implicititer);
    myexit(468346321);
  }

  // same list and numbers in array for both primitives and conserved
  
  //default
  *numdims=NDIM;
  DLOOPA(jj) erefU[jj]=UU+jj;
  DLOOPA(jj) eotherU[jj]=URAD0+jj;

  if(*implicitferr==QTYUMHD){
  }
  else if(*implicitferr==QTYUMHDENERGYONLY){
  }
  else if(*implicitferr==QTYUMHDMOMONLY){
  }
  else if(*implicitferr==QTYENTROPYUMHD){
    erefU[TT]=ENTROPY;
  }
  else if(*implicitferr==QTYENTROPYUMHDENERGYONLY){
    erefU[TT]=ENTROPY;
  }
  else if(*implicitferr==QTYENTROPYUMHDMOMONLY){
    erefU[TT]=ENTROPY;
  }
  else if(*implicitferr==QTYURAD){
    *numdims=NDIM;
    DLOOPA(jj) erefU[jj]=URAD0+jj;
    DLOOPA(jj) eotherU[jj]=UU+jj;
  }
  else if(*implicitferr==QTYURADENERGYONLY){
    *numdims=NDIM;
    DLOOPA(jj) erefU[jj]=URAD0+jj;
    DLOOPA(jj) eotherU[jj]=UU+jj;
  }
  else if(*implicitferr==QTYURADMOMONLY){
    *numdims=NDIM;
    DLOOPA(jj) erefU[jj]=URAD0+jj;
    DLOOPA(jj) eotherU[jj]=UU+jj;
  }
  else if(PRODUCTION==0){
    dualfprintf(fail_file,"No such implicitferr=%d\n",*implicitferr);
    myexit(468346322);
  }

  // sign that goes into implicit differencer that's consistent with sign for *signgd of -1 when using the radiative uu to measure f.
  *signgd2=(+1.0);
  *signgd4=(+1.0); // for entropy alone for Gdpl in error function // Appears for QTYUMHD,QTYENTROPYUMHD this sign is the right one.  But both cases have lots of cold MHD inversions.
  *signgd6=(-1.0); // for entropy as goes into GS from dUrad or dUmhd  //  // KORALTODO SUPERGODMARK:  -- unsure about sign!

}



//uu0 - original cons. qty
//uu -- current iteration
//f - (returned) errors

// returns error function for which we seek to be zero using Newton's method, which solves the implicit problem.
// failreturnallowable : what failure level to allow so that push through and continue despite the failure
// whichcall : which call to this function
// showmessages:
// allowlocalfailurefixandnoreport:
// pp : primitive (associated with uu) used as guess for inversion as well as for returning inversion result for later quicker inversion
// uu0 : reference initial conserved quantity representing U_n for implicit problem
// uu : current conserved quantity representing U_{n+1} that solves implicit problem
// localdt : timestep for 4-force
// ptrgeom:
// f : error function returrned
// fnorm : norm of error function for esimating relative error in f.
// 
// returns: eomtype, radinvmod,pp,uu,q,f,fnorm,freport,goexplicit,errorabs,converturn

static int f_implicit(int allowbaseitermethodswitch, int iter, int f1iter, int failreturnallowable, int whichcall, FTYPE impeps, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtype, int whichcap, int itermode, int *baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE conv, FTYPE convabs, FTYPE allowconvabs, int maxiter, FTYPE realdt, int dimtypef, FTYPE *dimfactU, FTYPE *ppprev, FTYPE *pp, FTYPE *piin, FTYPE *uuprev, FTYPE *Uiin, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom, struct of_state *q,  FTYPE *f, FTYPE *fnorm, FTYPE *freport, int *goexplicit, FTYPE *errorabs, FTYPE *errorallabs, int whicherror, int *convreturn, int *nummhdinvsreturn)
{
  int pliter, pl;
  int iv;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  // set inputs for errors, maxiters, etc.
  if(iter>=ITERMHDINVTRYHARDER || whichcall==FIMPLICITCALLTYPEFINALCHECK || whichcall==FIMPLICITCALLTYPEFINALCHECK2){
    // try lowest error allowed, may be raised a bit in MHD inversion code.
    // min between normal desired error and iterated-quantity error.  This seeks low error if iterated got low error, while avoids excessive attempt if not.
    newtonstats.tryconv=1E-2*MIN(convabs,*errorabs); // NUMEPSILON
    newtonstats.tryconvultrarel=1E-2*MIN(convabs,*errorabs); // NUMEPSILON
    newtonstats.extra_newt_iter=1; // KORALNOTE: apparently should keep this as >=1 to ensure error really drops
    newtonstats.extra_newt_iter_ultrarel=2; // KORALNOTE: apparently should keep this as >=1 to ensure error really drops
  }
  else{
    // KORALNOTE: If make newtonstats.tryconv~convabs, then if convabs~1E-12, then MHD inversion may return error~1E-10 in terms of how measured with f_error_check(), so must try harder than expected.
    newtonstats.tryconv=convabs*1E-2;
    newtonstats.tryconvultrarel=convabs*1E-2; // just bit smaller, not as extreme as default
    newtonstats.extra_newt_iter=1; // KORALNOTE: apparently should keep this as >=1 to ensure error really drops
    newtonstats.extra_newt_iter_ultrarel=1; // KORALNOTE: apparently should keep this as >=1 to ensure error really drops
  }
  //  newtonstats.mintryconv=allowconvabs;
  newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION;
  newtonstats.maxiter=maxiter;
  // override with less strict error for Jacobian calculation
  if(whichcall==FIMPLICITCALLTYPEJAC){
    newtonstats.tryconv=MAX(impeps*1E-2,newtonstats.tryconv);
    newtonstats.tryconvultrarel=MAX(impeps*1E-3,newtonstats.tryconvultrarel);
  }

  // set whether should check inversion result inside Utoprimgen()
  int checkoninversiongas;
  int checkoninversionrad;
  if(whichcall==FIMPLICITCALLTYPEFINALCHECK || whichcall==FIMPLICITCALLTYPEFINALCHECK2){
    checkoninversiongas=checkoninversionrad=1;
  }
  else{
    // don't check since slows down code and could be good enough solution if original error says ok.
    checkoninversiongas=checkoninversionrad=0;
  }



 



  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.
  FTYPE Gdpl[NPR]={0.0},Gdplabs[NPR]={0.0}, Tgas={0.0};
  int failreturn;
  FTYPE uuabs[NPR]={0.0};

  // default is no failure
  failreturn=0;

  ////////
  //
  // 0) backup pp and uu and prepare alternative versions
  //
  /////////
  struct of_state qbackup;
  FTYPE ppbackup[NPR],uubackup[NPR]={0.0};
  FTYPE ppalt[NPR]={0.0},uualt[NPR]={0.0};
  int radinvmodbackup,radinvmodalt,failreturnalt;
  PLOOP(pliter,pl){
    ppalt[pl]=ppbackup[pl]=pp[pl];
    uualt[pl]=uubackup[pl]=uu[pl];
  }
  qbackup=*q;
  radinvmodbackup=*radinvmod;




  /////////
  //
  // setup method and signs
  //
  ////////
  int numdims,startjac,endjac,implicititer,implicitferr,BEGINMOMSTEPS,ENDMOMSTEPS,BEGINENERGYSTEPS,ENDENERGYSTEPS,BEGINFULLSTEPS,ENDFULLSTEPS,BEGINNORMALSTEPS,irefU[NDIM],iotherU[NDIM],erefU[NDIM],eotherU[NDIM],signgd2,signgd4,signgd6,signgd7;
  define_method(iter, eomtype, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
  get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);



  // optimize whether need to really compute entropy with log/pow so slow
  int needentropy=1; // default get uu[entropy] and q->entropy
  if(ENTROPYOPT){
    // whichcall==FIMPLICITCALLTYPEFINALCHECK means final check where if wasn't computing entropy during iterations, need at end so next RK substeps have it ready
    if(whichcall==FIMPLICITCALLTYPEFINALCHECK || whichcall==FIMPLICITCALLTYPEFINALCHECK2 || *eomtype==EOMENTROPYGRMHD || (implicititer==QTYENTROPYUMHDMOMONLY)||(implicititer==QTYENTROPYUMHDENERGYONLY)||(implicititer==QTYENTROPYUMHD || implicititer==QTYENTROPYPMHD) || (implicitferr==QTYENTROPYUMHD || implicitferr==QTYENTROPYUMHDENERGYONLY || implicitferr==QTYENTROPYUMHDMOMONLY) || (fracenergy>0.0 && fracenergy<1.0)){
      needentropy=1;
    }
    else needentropy=0;
  }

  //1) Reorder URAD method so:
  // a) URAD
  // b) URAD->PRAD
  // c) Re-get URAD(PRAD) in case floors/ceilings
  // d) Apply relative 4-force condition using this new URAD -> UGAS
  // e) UGAS->PGAS
  // f) Re-get UGAS(PGAS) in case floors/celings/failures/errors/etc.
  //
  // As long as MHD doesn't fail, then even if RAD hits ceilings, the solution will be in relative force balance. -- i.e. total energy-momentum conservation will hold despite change in URAD.
  //
  // Right now, I use pre-corrected URAD to get dUrad -> dUgas, so if rad hits ceilings while gas does not, then relative force balance is lost when could have been maintained.
  //
  // Need to have Utoprimgen() call without doing radiation inversion to save time .. doradonly=-1 ?
  //
  // BUT, if hit radinvmod with change in energy, then change in U would be dumped into GAS even if gas<<rad or tau\sim 0.

  if(
     implicititer==QTYUMHD || implicititer==QTYUMHDENERGYONLY || implicititer==QTYUMHDMOMONLY
     || implicititer==QTYURAD || implicititer==QTYURADENERGYONLY || implicititer==QTYURADMOMONLY
     ){ // then don't yet have updated primitives, so invert to get them
    // 1) get change in conserved quantity between fluid and radiation (equal and opposite 4-force)
    // required for inversion to get P(U) for MHD and RAD variables
    // this preserves machine accurate conservation instead of applying 4-force on each fluid and radiation separately that can accumulate errors
    FTYPE Gddt[NDIM]; DLOOPA(iv) Gddt[iv]=(uu[irefU[iv]]-uu0[irefU[iv]]);
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - Gddt[iv];
    // reject what radiation says gas should be if forced into T^t_t<0 regime.
    int badchange=0;
    if(implicititer==QTYURAD || implicititer==QTYURADENERGYONLY || implicititer==QTYURADMOMONLY){
      if(
         REMOVERESTMASSFROMUU==2 && (uu[RHO]-uu[iotherU[0]]<=0.0)
         ||REMOVERESTMASSFROMUU!=2 && (-uu[iotherU[0]]<=0.0)
         ){
        badchange=1;
      }
      // can also go bad when momentum exceeds energy, but more complicated check for GRMHD -- so just do inversion
    }
    if(implicititer==QTYUMHD || implicititer==QTYUMHDENERGYONLY || implicititer==QTYUMHDMOMONLY){
      if(-uu[iotherU[0]]<=0.0){
        badchange=1;
      }
      // inversion will catch when momentum exceeds energy.
    }
    // "revert" uu and ignore 4-force if badchange and just fail directly without expense of inversion.
    if(badchange && iter==1){
      // if first iteration, assume guess as (e.g.) Uiin, but can now change to uu0.
      PLOOP(pliter,pl) uu[pl] = uu0[pl];
    }
    else if(badchange){
      // only change other to uubackup (i.e. before 4-force applied)
      DLOOPA(iv) uu[iotherU[iv]] = uubackup[iotherU[iv]];
    }
    // 2) Get estimated U[entropy]
    // mathematica has Sc = Sc0 + dt *GS with GS = -u.G/T
    FTYPE GS=0.0;
    FTYPE Tgaslocal=0.0;
    if(badchange==0){
      if(0){
        Tgaslocal=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pp[RHO],pp[UU]);
        get_state(pp, ptrgeom, q);
        DLOOPA(iv) GS += (-q->ucon[iv]*signgd2*(signgd7*Gddt[iv]))/(Tgaslocal+TEMPMIN); // maybe more accurate than just using entropy from pp and ucon[TT] from state from pp. FUCK: Can't be right that this is the same (signgd2=signgd7=1) for both URAD and UMHD methods.
      }
      else{
        // Get GS completely consisent with primitives, in case using entropy error function, because then shouldn't use Gddt[TT] related to energy equation.
        // Below rad inv may not be completely necessary, but not too expensive, so ok.
        int doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
        int computestate=1;
        int computeentropy=1;
        koral_source_rad_calc(computestate,computeentropy,pp, ptrgeom, Gdpl, Gdplabs, NULL, &Tgaslocal, q);
        GS = -signgd4 * localdt * Gdpl[ENTROPY]/(signgd6); // Consistent with uu = uu0 + signgd6*GS and how used when getting error function later
      }
    }
    uu[ENTROPY] = uu0[ENTROPY] + signgd6*GS; // KORALTODO SUPERGODMARK: Problem with UMHD,UMHD no matter signgd7.  Ok with URAD,URAD.   Ok with UMHD,ENTROPYUMHD if signgd7 +1 and signgd4 +1.
    // 3) Do MHD+RAD Inversion
    //    PLOOP(pliter,pl) dualfprintf(fail_file,"BEFORE: pl=%d pr=%g\n",pl,pp[pl]);
    int doradonly=0; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
    *nummhdinvsreturn++;

    if(failreturn==UTOPRIMGENWRAPPERRETURNFAILRAD && whichcall==FIMPLICITCALLTYPEJAC){
      // then try other uu, even if total energy is not conserved, better than raditive failure leading to problems or innacuracy.
      PLOOP(pliter,pl) uu[pl] = uuprev[pl];
      doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
      if(failreturn==UTOPRIMGENWRAPPERRETURNNOFAIL){
        int jjj;
        //        DLOOPA(jjj) dualfprintf(fail_file,"DIDIT: jjj=%d uu=%21.15g\n",jjj,uu[URAD0+jjj]);
      }
      else{
        //        dualfprintf(fail_file,"DIDITNOT\n");      
        int jjj;
        SLOOPA(jjj){
          uu[URAD0+jjj]*=0.1;
        }
        //        DLOOPA(jjj) dualfprintf(fail_file,"trying: jjj=%d uu=%21.15g\n",jjj,uu[URAD0+jjj]);
        doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
        //doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu0, q, ptrgeom, dissmeasure, pp, &newtonstats);
        //        if(failreturn==UTOPRIMGENWRAPPERRETURNNOFAIL) dualfprintf(fail_file,"DIDIT02\n");
        //        else dualfprintf(fail_file,"DIDITNOT02\n");      
      }
    }



    // When failreturn==UTOPRIMGENWRAPPERRETURNFAILMHD, Utoprimgen() passes back pp as initial guess, and when u(pp) is computed below, effectively didn't move uu[gas] away from uu0[gas]
    if(badchange) failreturn==UTOPRIMGENWRAPPERRETURNFAILMHD; // indicate that messed-up MHD inversion, without having to do the inversion.
    //    PLOOP(pliter,pl) dualfprintf(fail_file,"AFTER: pl=%d pr=%g\n",pl,pp[pl]);
    if(debugfail>=3) dualfprintf(fail_file,"NEWTON: %d : ijk=%d %d %d : %d %g\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k,newtonstats.lntries,newtonstats.lerrx);
    radinvmodalt=*radinvmod; // default
    failreturnalt=failreturn; // default

    if(*eomtype==EOMCOLDGRMHD){// restore uu because this inversion assumes u=0 and sets u=0, but if cold valid then ok to keep u~0 and constant as long as fixups applied.
      pp[UU]=ppbackup[UU];
      if(ENTROPY>=0) pp[ENTROPY]=pp[UU];
    }
    // KORALTODO: now can check if actually did eomtype==EOMGRMHD or EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function if using QTYUMHD.
    // 4) now get consistent uu[] based upon actual final primitives.
    // Needed in case inversion reduced to entropy or cold.
    // Also needed to get correct UU[ENTROPY] if did full MHD inversion or get correct U[UU] if did entropy inversion, etc.
    get_state(pp, ptrgeom, q);
    primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
    // now have full primitives and full U including entropy and these are consistent with each other.
    // 5)
    // set alternative that doesn't keep any changes to iterated quantities
    PLOOP(pliter,pl) ppalt[pl] = pp[pl];
    PLOOP(pliter,pl) uualt[pl] = uu[pl];
    DLOOPA(iv) uualt[irefU[iv]] = uubackup[irefU[iv]];
    if(implicititer==QTYURAD || implicititer==QTYURADENERGYONLY || implicititer==QTYURADMOMONLY){
      radinvmodalt=0; // if iterating URAD, then ignore radinvmod==1 for alternative error but need to still use new primitives.
      if(failreturn==UTOPRIMGENWRAPPERRETURNFAILRAD) failreturnalt=UTOPRIMGENWRAPPERRETURNNOFAIL; // only change failure if only rad failure
    }
  }
  else if(implicititer==QTYENTROPYUMHD || implicititer==QTYENTROPYUMHDENERGYONLY || implicititer==QTYENTROPYUMHDMOMONLY){
    // so iterating U[ENTROPY,U1,U2,U3]
    //    FTYPE uuorig[NPR]; PLOOP(pliter,pl) uuorig[pl]=uu[pl];
    // 1) Do pure ENTROPYMHD inversion to get pmhd
    int doradonly=0; int eomtypetemp=EOMENTROPYGRMHD; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, &eomtypetemp, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
    *nummhdinvsreturn++;
    radinvmodalt=*radinvmod; // default
    failreturnalt=failreturn; // default
    if(*eomtype==EOMCOLDGRMHD){// restore uu because this inversion assumes u=0 and sets u=0, but if cold valid then ok to keep u~0 and constant as long as fixups applied.
      pp[UU]=ppbackup[UU];
      if(ENTROPY>=0) pp[ENTROPY]=pp[UU];
    }
    // set eomtype!  So for entire process for this cell, the entropy equations are used.
    *eomtype=EOMENTROPYGRMHD;
    // KORALTODO: now can check if actually did eomtypetemp==EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function in case reduced to cold MHD (maybe G=0?)
    // 2) get state (mhd and entropy are up-to-date, while rad is out-of-date, but fixed below)
    get_state(pp, ptrgeom, q);
    // 3) Compute U[UU] so have it for below (computes U[Ui,ENTROPY,RAD], but those not used or already known.)
    FTYPE uuentropy[NPR],uuentropyabs[NPR];
    primtoU(UNOTHING,pp,q,ptrgeom, uuentropy, uuentropyabs);
    uu[ENTROPY]=uuentropy[ENTROPY];
    uuabs[ENTROPY]=uuentropyabs[ENTROPY];
    //    PLOOP(pliter,pl) if(pl==U1||pl==U2||pl==U3||pl==ENTROPY) pp[pl]=pporig[pl]; // get back U[ENTROPY,Ui] so no machine error introduced from u(p(u)) for actually known quantities that won't change due to fixups or anything. -- no, primitives from Utoprimgen might change things and act as Newton step fix.
    // 4) Get actual Urad [note uses UU not irefU that is ENTROPY for irefU[0]]
    DLOOPA(iv) uu[URAD0+iv] = uu0[URAD0+iv] - (uu[UU+iv]-uu0[UU+iv]);
    // 5) Do RAD-ONLY inversion
    int failreturn2;  doradonly=1; failreturn2=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
    //  no need to concern with eomtype in RAD only case.  i.e. eomtype won't change.
    if(failreturn2>failreturn) failreturn=failreturn2; // use worst case.
    // 6) now get consistent uu[] based upon actual final primitives.
    // Needed in case inversion reduced to entropy or cold or radiative inversion used fixups.
    get_state(pp, ptrgeom, q);
    primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
    // now have full primitives and full U including entropy and these are consistent with each other.
    // 7)
    // set alternative that doesn't keep any changes to iterated quantities
    PLOOP(pliter,pl) ppalt[pl] = pp[pl];
    PLOOP(pliter,pl) uualt[pl] = uu[pl];
    DLOOPA(iv) uualt[irefU[iv]] = uubackup[irefU[iv]];
  }
  else if(implicititer==QTYENTROPYPMHD){
    // not setup
    myexit(92846534);
  }
  else if(implicititer==QTYPMHD || implicititer==QTYPMHDENERGYONLY || implicititer==QTYPMHDMOMONLY){
    // 0) Have pmhd={ug,uvel1,uvel2,uvel3}
    //
    // 0.5) trivially invert field
    PLOOPBONLY(pl) pp[pl] = uu0[pl];
    // set p[ENTROPY] in case evaluated as a diagnostic.
    if(ENTROPY>=0) pp[ENTROPY] = pp[UU];
    //
    // 1) get state of non-density related things (mhd state needed for q->ucon[TT] and primtoU) and Umhd, Uentropy [also computes Urad, but overwritten next and not expensive]
    get_state_norad_part1(pp, ptrgeom, q);
    //
    // 2) trivially invert to get rho
    pp[RHO]= uu0[RHO]/q->ucon[TT];
    //
    if(1){
      // Don't restrict rho or u_g except as by iteration catches.  Don't restrict gamma or gammarad during iterations so can smoothly go out of bounds if required temporarily.
      // KORALTODO: SUPERGODMARK: causes many more high error - high iteration events.  fixup1zone() has no radiation, so not consistent.  But eventually applied!  Is it due to limit_gamma() or density floors?
      // 3) Apply any fixups, like floors or limit_gamma's.  Won't help reduce error, but even if solution with high fluid gamma, later apply limit_gamma, but then balance between fluid and radiation lost.  So better to see as high error event than accurate high gamma event.
      int finalstepfixup=0; // treat as if not needing to diagnose, just act.  KORALTODO: Although, Utoprimgen(finalstep) used to get change in fluid primitives.
      // Note, uu is old, but only used for diagnostics, and don't care about diagnostics in this stepping.
      FTYPE ppfixup[NPR],ppfloor[NPR],uufixup[NPR];
      PLOOP(pliter,pl){
        ppfixup[pl]=ppfloor[pl]=pp[pl];
        uufixup[pl]=uu[pl];
      }
      // bsq is accurate using below
      FTYPE bsq; bsq_calc_fromq(ppfixup, ptrgeom, q, &bsq);
      // uu isn't exactly like pfixup here, but close enough
      set_density_floors_alt(ptrgeom, q, ppfixup, uu, bsq, ppfloor);
      //      fixup1zone(ppfloor,uufixup, ptrgeom,finalstepfixup); // too complicated for implicit stepping given how rare shoul be used.
      if(pp[RHO]<0.0) pp[RHO]=ppfloor[RHO]; // only fix RHO if really went negative.  Not smooth, but avoids problems in difficult regimes.

      //      limit_gamma(0,GAMMAMAX,GAMMAMAXRADIMPLICITSOLVER,pp,NULL,ptrgeom,0);
      // fix uu[RHO] to be consistent, since uu[RHO] used to get inverted scalars
      uu[RHO] = pp[RHO]*q->ucon[TT];
    }
    //
    // 4) Invert other scalars (only uses uu[RHO], not pp[RHO])
    extern int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr);
    invert_scalars(ptrgeom, uu,pp);
    //
    // save final fixed-up pp's
    FTYPE pporig[NPR]; PLOOP(pliter,pl) pporig[pl]=pp[pl];
    // NOW all MHD pp's are set and backed-up
    //
    // 5) Do rest of get_state that used rho, and other scalars, in case used.
    // This computes pressure (as required for T^t_\mu) and entropy (as required for primtoflux_nonradonly below for entropy flux)
    // KORALTODO: Although, don't need entropy if doing implicitferr==UMHD
    get_state_norad_part2(needentropy, pp, ptrgeom, q); // where entropy would be computed

    //
    // 6) Compute Umhd and Uentropy (keeps Urad as zero, but Urad set next)
    //primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
    extern int primtoflux_nonradonly(int needentropy, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
    FTYPE uumhd[NPR],uumhdabs[NPR];
    primtoflux_nonradonly(needentropy,pp,q,TT,ptrgeom, uumhd, uumhdabs); // anything not set is set as zero, which is rad.
    PLOOP(pliter,pl) if(!RADUPL(pl)) uu[pl]=uumhd[pl];
    PLOOP(pliter,pl) if(!RADUPL(pl)) uuabs[pl]=uumhdabs[pl];
    //    primtoflux_nonradonly(1,pp,q,TT,ptrgeom, uu, uuabs); // doesn't actually compute entropy again, just multiplies existing things.
    //    if(needentropy==0) uu[ENTROPY]=sqrt(-1.0);
    
    // 7) Get actual Urad(G) via energy conservation (correct even if using entropy as error function, because just computed correct U[ENTROPY] consistent with U[UU].
    // KORALNOTE: uu set by p->uu, which has an error of order machice precision.  So even if primitives didn't change, uu-uu0 can be order machine precision.  So when fluid uu>>rad uu, this can lead to the below giving huge radiation changes even if actually fluid background is not changing.
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - (uu[irefU[iv]]-uu0[irefU[iv]]);
    //
    // 8) Do RAD-ONLY Inversion (eomtype not used)
    int whichcapnew;
    if(USECAPTYPEFIX2FORF1 && whichcall==FIMPLICITCALLTYPEF1) whichcapnew=CAPTYPEFIX2;
    else if(USECAPTYPEFIX2FORFINALCHECK && whichcall==FIMPLICITCALLTYPEFINALCHECK) whichcapnew=CAPTYPEFIX2;
    else whichcapnew=whichcap; // for jacobian, likely to be CAPTYPEFIX2 anyways.
    int doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcapnew, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
    radinvmodalt=*radinvmod; // default
    failreturnalt=failreturn; // default
    //  no need to concern with eomtype in RAD only case.  i.e. eomtype won't change.
    //    if(*radinvmod>0) dualfprintf(fail_file,"radinvmodinside=%d  whichcall=%d : u=%g Erf=%g\n",*radinvmod,whichcall,pp[UU],pp[PRAD0]);


    if(0){
      // deal with Erf<0 and possibly gammarad caps or E_r<0 issues
      if(iter<BEGINNORMALSTEPS+4 && pp[PRAD0]<10.0*ERADLIMIT && ppbackup[PRAD0]>10.0*ERADLIMIT){
        // then can play with Erf in case negative E_r to avoid bad NR due to floor on Erf.
        PLOOP(pliter,pl) if(RADPL(pl)) pp[pl]=ppbackup[pl]; // only concern is error might coincidentally be small
        dualfprintf(fail_file,"Caught\n");
      }
      if(myid==8){
        dualfprintf(fail_file,"bob: iter=%d (%d) %g %g (%g)\n",iter,BEGINNORMALSTEPS+4,pp[PRAD0],ppbackup[PRAD0],10.0*ERADLIMIT);
      }
    }

    //
    // 9) Get consistent RAD state
    // only changes radiation state, not fluid state, but keeps old q's for fluid state for next step.
    get_state_radonly(pp, ptrgeom, q);
    //
    // fix-up primitives to avoid violent steps in temperature
    //
    // 10) Get new uu[RAD] since original uu[RAD]->pp[RAD] might have had fixups applied and then uu[RAD] no longer consistent.
    // This violates total energy conservation in favor of consistency of radiative quantities between pp and uu
    //    primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
    extern int primtoflux_radonly(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxabs);
    FTYPE uurad[NPR],uuradabs[NPR];
    primtoflux_radonly(pp,q,TT,ptrgeom, uurad,uuradabs); // all non-rad stuff is set to zero.
    // write new uurad's to uu
    PLOOP(pliter,pl) if(RADUPL(pl)) uu[pl]=uurad[pl];
    PLOOP(pliter,pl) if(RADUPL(pl)) uuabs[pl]=uuradabs[pl];
    //
    // 11) Recover actual iterated pmhd to avoid machine related differences between original pp and pp(U(pp)) for pmhd quantities
    // This assumes that iterated pmhd is optimal and not modified except by iteration by Newton step, which is currently true.
    // This gives machine error priority to mhd primitives rather than mhd U's.
    // below not necessary since don't overwrite pp[non-rad]
    //PLOOP(pliter,pl) if(!RADPL(pl)) pp[pl]=pporig[pl];
    //
    // 12) overwrite any setting of uu[B1,B2,B3], so machine accurate field
    PLOOPBONLY(pl) uu[pl]=pp[pl];
    //
    // now have full primitives (pp) and full U (uu) including entropy and these (pp and uu) are fully consistent with each other.
    // 13)
    // set alternative that doesn't keep any changes to iterated quantities
    PLOOP(pliter,pl) ppalt[pl] = pp[pl];
    DLOOPA(iv) ppalt[irefU[iv]] = ppbackup[irefU[iv]];
    PLOOP(pliter,pl) uualt[pl] = uu[pl];

  }
  else if(implicititer==QTYPRAD || implicititer==QTYPRADENERGYONLY || implicititer==QTYPRADMOMONLY){
    // 0.5) trivially invert field
    PLOOPBONLY(pl) pp[pl] = uu0[pl];
    // Have prad={Erf,uradvel1,uradvel2,uradvel3}
    // 1) get state (rad state needed for mhd_calc_rad, while old mhd state needed to estimate uu[ENTROPY] below)
    get_state(pp, ptrgeom, q);
    // 2) Compute Urad[prad0,uradcon,uradcov] [also computes old Umhd and old Uentropy, but not used]
    FTYPE uurad[NPR],uuradabs[NPR];
    primtoU(UNOTHING,pp,q,ptrgeom, uurad, uuradabs);
    PLOOP(pliter,pl) if(RADPL(pl)){
      uu[pl]=uurad[pl];
      uuabs[pl]=uuradabs[pl];
    }
    // 3) Get actual Umhd(G) via energy conservation:
    FTYPE Gddt[NDIM]; DLOOPA(iv) Gddt[iv]=(uu[irefU[iv]]-uu0[irefU[iv]]);
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - Gddt[iv];
    // uu0[RHO] doesn't change
    // 4) Estimate U[ENTROPY](G) using old rho,u
    // UNSURE what to do, since blows up either way:
    //    FTYPE GS=0.0; DLOOPA(iv) GS += (-q->ucon[iv]*signgd2*(signgd7*Gddt[iv]))/(Tgaslocal+TEMPMIN); // more accurate than just using entropy from pp and ucon[TT] from state from pp.
    FTYPE GS=0.0;
    FTYPE Tgaslocal;
    if(0){
      Tgaslocal=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pp[RHO],pp[UU]);
      DLOOPA(iv) GS += (-q->ucon[iv]*signgd2*(signgd7*Gddt[iv]))/(Tgaslocal+TEMPMIN); // more accurate than just using entropy from pp and ucon[TT] from state from pp.
    }
    else{
      // Get GS completely consisent with primitives, in case using entropy error function, because then shouldn't use Gddt[TT] related to energy equation.
      // Below rad inv may not be completely necessary, but not too expensive, so ok.
      int computestate=0; // already computed above
      int computeentropy=1;
      koral_source_rad_calc(computestate,computeentropy,pp, ptrgeom, Gdpl, Gdplabs, NULL, &Tgaslocal, q);
      GS = -signgd4 * localdt * Gdpl[ENTROPY]/(signgd6); // so uu = uu0 + signgd6*GS is consistent with how Gdpl included in error function later.
    }
    uu[ENTROPY] = uu0[ENTROPY] + signgd6*GS;
    // 5) Invert to get pmhd (also does rad inversion, but not expensive so ok)
    int doradonly=0; failreturn=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, eomtype, whichcap, EVOLVEUTOPRIM, UNOTHING, uu, q, ptrgeom, dissmeasure, pp, &newtonstats);
    *nummhdinvsreturn++;
    radinvmodalt=*radinvmod; // default
    failreturnalt=failreturn; // default
    if(*eomtype==EOMCOLDGRMHD){// restore uu because this inversion assumes u=0 and sets u=0, but if cold valid then ok to keep u~0 and constant as long as fixups applied.
      pp[UU]=ppbackup[UU];
      if(ENTROPY>=0) pp[ENTROPY]=pp[UU];
    }
    // KORALTODO: now can check if actually did eomtype==EOMGRMHD or EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function
    // 6) Recover actual iterated prad to avoid machine related differences between original pp and pp(U(pp)) for prad quantities
    // This assumes that iterated prad is optimal and not modified except by iteration by Newton step
    PLOOP(pliter,pl) if(RADPL(pl)) pp[pl]=ppbackup[pl];
    // 7) Get consistent Urad [also computes Umhd and Uentropy, which is ok]
    get_state(pp, ptrgeom, q);
    primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
    // 8) overwrite any setting of uu[B1,B2,B3], so machine accurate field
    PLOOPBONLY(pl) uu[pl]=pp[pl];
    // now have full primitives and full U including entropy and these are consistent with each other.
    // 9) 
    // set alternative that doesn't keep any changes to iterated quantities
    PLOOP(pliter,pl) ppalt[pl] = pp[pl];
    DLOOPA(iv) ppalt[irefU[iv]] = ppbackup[irefU[iv]];
    PLOOP(pliter,pl) uualt[pl] = uu[pl];
  }




  ///////////
  //
  // check which baseitermethod we should really be using based upon any inversions done above
  //
  /////////
  // only check if not in Jacobian or final check
  if(allowbaseitermethodswitch && whichcall==FIMPLICITCALLTYPEF1){//FIMPLICITCALLTYPEJAC

    int doswitchbaseitermethod=0;
    FTYPE rdU[NPR];
    //  PLOOP(pliter,pl) rdU[pl]=Gallabs[pl]/uuallabs[pl];
    // Note that we use uu-uu0 since that's exactly what's pushed into other non-iterated set of equations.
    FTYPE uuallabs1;
    PLOOP(pliter,pl){
      uuallabs1 = THIRD*(fabs(uuabs[pl]) + fabs(uu[pl]) + fabs(uu0[pl]));
      rdU[pl]=(uu[pl]-uu0[pl])/uuallabs1;
    }
    // if dU[UU] changes relatively near machine precision, then could be fake change induced by machine precision errors.  Then actual real evolution of URAD? would be destroyed by those machine errors.
    // For iter=1, using initial guess's uu.  But by second iteration, have used Jacobian to move iterates.  If dU[UU]/U[UU] near machine precision but dU[URAD0]/U[URAD0] far from it, then should use radiation iterate.
    if(iter>1){
      if(rdU[UU]<rdU[URAD0] && IMPMHDTYPEBASE(*baseitermethod)==1){
        if(debugfail>=3) dualfprintf(fail_file,"Should switch base to QTYURAD or QTYPRAD\n");
      }
      if(rdU[UU]>rdU[URAD0] && (*baseitermethod==QTYURAD||*baseitermethod==QTYPRAD)){
        if(debugfail>=3) dualfprintf(fail_file,"Should switch base to QTYPMHD\n");
      }
    }
    if(pp[URAD0]<10.0*ERADLIMIT && IMPMHDTYPEBASE(*baseitermethod)==1){
      if(debugfail>=3) dualfprintf(fail_file,"Should switch base to QTYURAD or QTYPRAD: Erf=%21.15g<10*%21.15g baseitermethod=%d\n",pp[URAD0],ERADLIMIT,*baseitermethod);
      if(1){
        // then Erf dropped-out.  Probably because using QTYPMHD and too big change
        *baseitermethod=QTYPRAD; // prad safest and fastest
        //*baseitermethod=QTYURAD; // URAD somehow more robust (e.g. on RADTUBE PRAD methods fails a few times if switched from QTYPMHD, but fails more if used directly!)
        doswitchbaseitermethod=1;
        PLOOP(pliter,pl){
          pp[pl] = ppprev[pl]; // revert primitives from old primitive iteration to new primitive iteration
          uu[pl] = uuprev[pl]; // start with previous uu's before any pp->uu calculations
        }
        // get state since not synched with q necessarily
        get_state(pp, ptrgeom, q);
        primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
        //*q=qbackup; // revert to q before any pp->q calculations
        *radinvmod=radinvmodbackup; // revert to previous radinvmod
      }
      else{
        // problem with switching for RADTUBE, so just report if want.  If had backup method that was radiation iterate, then might want to break, but not always as QTYPMHD method can recover.  So just stick with normal backup procedure...?
      }
    }
    if(pp[UU]<10.0*UUMINLIMIT &&  (*baseitermethod==QTYURAD||*baseitermethod==QTYPRAD)){
      // Then u_g dropped-out.  Probably because iteration radiation quantity
      if(debugfail>=3) dualfprintf(fail_file,"Should switch base to QTYMHD: ug=%21.15g\n",pp[UU]);
    }

    if(doswitchbaseitermethod){
      // in case required, change settings
      define_method(iter, eomtype, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
      get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);
    }
  }


  /////////
  //
  // At this point, must have pp, uu, uuabs, and q all consistent and defined
  //
  /////////



  ////////////////
  // get 4-force for all pl due to radiation as due to pp[uu]
  ////////////////
  int computestate=0;// already computed above
  int computeentropy=needentropy;
  FTYPE chieff;
  koral_source_rad_calc(computestate,computeentropy,pp, ptrgeom, Gdpl, Gdplabs, &chieff, &Tgas, q);
  FTYPE tautot,tautotmax;
  calc_tautot_chieff(pp, chieff, ptrgeom, &tautot, &tautotmax);







  ///////////////////////////////
  // GET ERROR FUNCTION
  ///////////////////////////////


  // get abs versions for all pl
  FTYPE uuallabs[NPR]={0.0},Gallabs[NPR]={0.0};

  FTYPE sign[NPR],extrafactor[NPR];
  PLOOP(pliter,pl){
    sign[pl]=signgd2;
    extrafactor[pl]=1.0;
  }
  if(ENTROPY>=0){
    pl=ENTROPY;
    sign[pl]=signgd4;
    // replace original equation with dS*T equation
    // error function is T*dS so no actual division by T.  Found in mathematica that this works best in difficult precision cases.
    extrafactor[pl]=fabs(Tgas)+TEMPMIN;
  }

  // get f, uuallabs, Gallabs, and fnorm
  // also get falt that uses (in case iterated) originally iterated value of uu -- instead of post-inversion modified version that would be consistent with primitives.
  FTYPE falt[NPR];
  //
  PLOOP(pliter,pl){
    f[pl] = ((uu[pl] - uu0[pl]) + (sign[pl] * localdt * Gdpl[pl]))*extrafactor[pl];

    falt[pl] = ((uualt[pl] - uu0[pl]) + (sign[pl] * localdt * Gdpl[pl]))*extrafactor[pl];

    // get error normalization that involves actual things being differenced
    // KORALNOTE: Use Gdplabs to ensure absolute value over more terms in source in case coincidental cancellation without physical significance.
    uuallabs[pl] = THIRD*(fabs(uuabs[pl]) + fabs(uu[pl]) + fabs(uu0[pl]))*extrafactor[pl];
    Gallabs[pl] = fabs(sign[pl] * localdt * Gdplabs[pl])*extrafactor[pl];
    fnorm[pl] = uuallabs[pl] + Gallabs[pl];
  }




  // compute difference vector between original and new 4-force's effect on conserved radiative quantities
  // NR1992 Eq. 16.6.16: y_{n+1} = y_n + h f(y_{n+1}) , so error function is f = (y_{n+1} - y_n) - h f(y_{n+1})
  // i.e. f->0 as change in conserved quantity approaches the updated value of 4-force
  if(0){
    if(fracenergy>0.0 && fracenergy<1.0){ // NOT USED!
      // then erefU=UU is set as default error term, so add entropy to this.

      FTYPE fentropy[NPR];
      FTYPE fnormentropy[NPR];

      //    fracenergy=1.0;
      if(startjac==TT){
        // Get final interpolated energy-entropy-term error function
        //      dualfprintf(fail_file,"fracenergy=%g fUU=%g fnormUU=%g fE=%g fnormE=%g\n",fracenergy,f[UU],fnorm[UU],fentropy[ENTROPY],fnormentropy[ENTROPY]);
        f[UU] = fracenergy*fabs(f[UU]) + (1.0-fracenergy)*fabs(f[ENTROPY]);
        fnorm[UU] = fracenergy*fnorm[UU] + (1.0-fracenergy)*fnorm[ENTROPY];
      }
    }
  }


  ///////
  //
  // get single number that measures error for iterated quantities
  //
  //////
  *convreturn=f_error_check(showmessages, showmessagesheavy, iter, conv, convabs, realdt, dimtypef,*eomtype , *radinvmod, itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,f,fnorm,freport,Uiin,uu0,uu,ptrgeom,errorabs,errorallabs,whicherror);


  /////////
  //
  // get alternative error
  //
  // Used in case optically thin and no interaction between RAD and GAS and radiation hits ceiling/floor and then must treat as ok since no force balance will help since no G and dUgas from dUrad will be wrong (and if u_g<<Erf, then would wrongly affect gas despite actually G<<U -- and if u_g>>Erf would no effect on gas be caused by that happening).
  //
  ////////
  if(ALLOWUSEUUALT){ // purpose of this is equivalent to whether one can go explicit, but overly complicated and if reduces to using alternative when not desired, then not energy-momentum conserving
    if(whichcall==FIMPLICITCALLTYPEF1 && (iter==1 && tautotmax<NUMEPSILON) && failreturn<=UTOPRIMGENWRAPPERRETURNFAILRAD){//FIMPLICITCALLTYPEJAC
      // Have to check if really small error for alternative error, since otherwise f1iter loop won't work.
      // Only check on iter==1 since if moved beyond iter=1, then must be relevant 4-force.  Only check if f1iter>1 because for f1iter=1 4-force is not set yet.

      // have to check if unmodified iterated quantities actually have smaller error, in which case probably U->P->Unew led to Unew different from U for iterated quantities due to no normal solution without fixups.  When gas~rad in energy density, the other compnent can abosorb that error, but in limit that (say) RAD<<GAS, then dRAD won't register to GAS values and RAD will never converge if iterated on.
      FTYPE errorabsalt,errorallabsalt;
      FTYPE faltreport[NPR];
      int convreturnalt=f_error_check(showmessages, showmessagesheavy, iter, conv, convabs, realdt, dimtypef,*eomtype , *radinvmod, itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,falt,fnorm,faltreport,Uiin,uu0,uualt,ptrgeom,&errorabsalt,&errorallabsalt,whicherror);

      //  if(nstep>=32 && nstep<=46){
      //    PLOOP(pliter,pl) dualfprintf(fail_file,"DEBUGIT: wc=%d iter=%d pl=%d uu0=%21.15g uu=%21.15g uualt=%21.15g Gdpl=%21.15g f=%21.15g falt=%21.15g errorabs=%21.15g errorabsalt=%21.15g radinvmod=%d radinvmodalt=%d failreturn=%d\n",whichcall,iter,pl,uu0[pl],uu[pl],uualt[pl],sign[pl]*localdt*Gdpl[pl],f[pl],falt[pl],*errorabs,errorabsalt,*radinvmod,radinvmodalt,failreturn);
      //  }



      // see if alternative error is smaller *and* in converged limit.  Catches cases when optically thin and ok that hits radiation ceiling.
      if(errorabsalt<*errorabs && convreturnalt==1){
        //dualfprintf(fail_file,"UsingALT\n");

        // then switch to use falt.
        // This typically controls what uu is used for the error
        // Eventually want uu to be consistent with pp
        // But allow for possible iteration upon uualt, so want to keep uu=uualt even though not consistent with pp
        // So have to fix uu(pp) once iterations are done.
        PLOOP(pliter,pl){
          uu[pl]=uualt[pl];
          pp[pl]=ppalt[pl];
          f[pl]=falt[pl];
          freport[pl]=faltreport[pl];
        }
        if(IMPUTYPE(implicititer)){
          // q based upon primitives, so if uu-method, then stuck with that q
        }
        else{
          // get new q(ppalt) now q(pp)
          get_state(pp, ptrgeom, q);
          // get new uu(ppalt) now uu(pp)
          primtoU(UNOTHING,pp,q,ptrgeom, uu, uuabs);
        }
        *convreturn=convreturnalt;
        *errorabs=errorabsalt;
        *errorallabs=errorallabsalt;
        *radinvmod=radinvmodalt;
        failreturn=failreturnalt;
      }
    }
  }



  /////////
  //
  // At this point, even if first iteration, know whether source term is what contributes to changes in uu.
  // If no absolute force to machine precision for each absolute uu, then implicit stepping can be avoided.
  // Even if inversions led to no consistent inversion (e.g. raditive inversion uses ceilings and so uu!=uu0 even for G=0), the below is correct.
  // This even accounts for case where entropy or energy suggest need implicit
  // This even accounts for when RAD quantities or MHD quantities differ on whether need explicit, since go over all pl always.
  //
  // NO: Apparently guess can lead to small G, but next iteration may not.  But generally iterations can slowly grow G, so can't use any iter condition either.
  // Ok, if iter=1, then check both tau and G (not either, but both) and also ensure only rad inversion changes, no gas inversion changes/issues/failures.
  ////////
  *goexplicit=0; // default
#define ITERCHECKEXPLICITSAFE 1 // iteration by which assume G has settled and can test if can go explicit.
  if(1){
    if(whichcall==FIMPLICITCALLTYPEF1){//FIMPLICITCALLTYPEJAC)

      FTYPE ucon0=q->ucon[TT]; // what enters G for dR^t_t/R^t_t from time part
      FTYPE ratchangeRtt=SMALL+fabs(chieff * ucon0 * ucon0 * realdt * 1.0); // 1.0 = c (2nd term with chi instead of kappaes to be sever and account for B-based term)

      //      if( (iter>ITERCHECKEXPLICITSAFE || iter==1 && tautotmax<NUMEPSILON ) && failreturn<=UTOPRIMGENWRAPPERRETURNFAILRAD){
      if( (iter>ITERCHECKEXPLICITSAFE || iter==1 && tautotmax<NUMEPSILON ) ){
        // iter>1 so at least have estimate of G even if not great.
        // At iter=1, U->p->G can give G=0, while dUrad=-dUgas can still lead to changes that upon next iteration lead to G!=0.
        *goexplicit=1;
        PLOOP(pliter,pl) if(fabs(Gallabs[pl])>NUMEPSILON*fabs(uuallabs[pl])) *goexplicit=0;
        pl=URAD0;
        if(fabs(ratchangeRtt*uu[pl])>NUMEPSILON*fabs(uuallabs[pl])) *goexplicit=0;
        if(fabs(ratchangeRtt*uu0[pl])>NUMEPSILON*fabs(uuallabs[pl])) *goexplicit=0;
        //  if(*goexplicit) dualfprintf(fail_file,"Went explicit\n");
        //  else dualfprintf(fail_file,"Stayed implicit\n");
      }
    }
  }







  if(debugfail>=3){
    ///////////
    //
    // get stats on any Utoprimgen() newton calls
    //
    /////////
    static long long int newtoncounttotal=0;
    newtoncounttotal+=newtonstats.lntries;
    static long long int newtoncounthere=0;
    newtoncounthere++;
    dualfprintf(fail_file,"Newtonstat: local=%d total=%d average=%21.15g\n",newtonstats.lntries,newtoncounttotal,(FTYPE)newtoncounttotal/(FTYPE)newtoncounthere);
  }  
  


  if(failreturn && failreturn>failreturnallowable){
    if(PRODUCTION==0&&showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of f_implicit(): %d vs. %d\n",failreturn,failreturnallowable);
    return(failreturn);
  }
  //  else{
  //    // save better guess for later inversion (including this inversion above) from this inversion
  //    PLOOP(pliter,pl) pp0[pl]=pp[pl];
  //  }





  return(0);
} 

// compute dt for this sub-step
static FTYPE compute_dt(FTYPE *CUf, FTYPE dtin)
{
  // what's applied to source and flux terms to get update (see global.stepch.h and step_ch.c:get_truetime_fluxdt() and step_ch.c:setup_rktimestep()) to get Uf
  // We don't use the ucum update version of dt.  As part of the RK method, the ucum update is separate from the substeps used to get information on updates(?). GODMARK.
  return(CUf[2]*dtin);
}




#define FAILRETURNGOEXPLICIT -1
#define FAILRETURNNOFAIL 0
#define FAILRETURNGENERAL 1
#define FAILRETURNJACISSUE 2
#define FAILRETURNMODESWITCH 3
#define FAILRETURNNOTTOLERROR 4

#define ACCEPTASNOFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNNOTTOLERROR || failreturn==FAILRETURNGOEXPLICIT)
//#define GOODNOFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNGOEXPLICIT)
#define NOTACTUALFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNMODESWITCH || failreturn==FAILRETURNGOEXPLICIT)
#define NOTBADFAILURE(failreturn) (failreturn==FAILRETURNNOFAIL || failreturn==FAILRETURNMODESWITCH  || failreturn==FAILRETURNNOTTOLERROR || failreturn==FAILRETURNGOEXPLICIT)

#define ACTUALHARDFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE)
#define ACTUALHARDORSOFTFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE || failreturn==FAILRETURNNOTTOLERROR)
#define SWITCHGOODIDEAFAILURE(failreturn) (failreturn==FAILRETURNGENERAL || failreturn==FAILRETURNJACISSUE || failreturn==FAILRETURNNOTTOLERROR || failreturn==FAILRETURNMODESWITCH)

#define RADINVBAD(radinvmod) (radinvmod==UTOPRIMRADFAILERFNEG || radinvmod==UTOPRIMRADFAILBAD1)
#define RADINVOK(radinvmod) (RADINVBAD(radinvmod)==0)


// choose to switch to entropy only if energy fails or gives u_g<0.  Or choose to always do both and use best solution.
//#define MODEMETHOD MODEPICKBEST // general switching method
//#define MODEMETHOD MODEPICKBESTSIMPLE // switches with only PMHD method
//#define MODEMETHOD MODEPICKBESTSIMPLE2 // switches with all methods but no ITERMODESTAGES attempted
//#define MODEMETHOD MODESWITCH
//#define MODEMETHOD MODEENERGY
#define MODEMETHOD MODEENTROPY
//#define MODEMETHOD MODEENERGYRAMESH


// KORALTODOMAYBE: average good neighbor if can.  only use alternative backup if no good neighbor.


// wrapper for mode method
static int koral_source_rad_implicit(int *eomtype, FTYPE *pb, FTYPE *pf, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE dissmeasure, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int i=ptrgeom->i;
  int j=ptrgeom->j;
  int k=ptrgeom->k;
  int p=ptrgeom->p;

  int pliter,pl;
  int sc;

  int failreturn,noprims;
  int havebackup,didentropyalready;
  int usedenergy=0,usedentropy=0,usedboth=0,usedcold=0,usedimplicit=0,usedexplicitgood=0,usedexplicitkindabad=0,usedexplicitbad=0;

  // set backups that might change and contaminate a fresh start
  // piin, Uiin, Ufin, CUf, ptrgeom, dUother don't change, rest can.
  int failfinalreturn;
  int eomtypelocal;



  /////////
  //
  // Setup pb and uub that hold final primitive and conserved solutions
  //
  /////////

  //////
  // set full uu0
  FTYPE fracdtuu0=1.0;
  FTYPE uu0[NPR],dUtot[NPR],rdUtot[NPR];
  PLOOP(pliter,pl){
    uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // initial+flux value of U
    dUtot[pl] = uu0[pl]-Uiin[pl]; // absolute change due to flux-advection step
    rdUtot[pl] = fabs(uu0[pl]-Uiin[pl])/(fabs(uu0[pl])+fabs(Uiin[pl])); // relative change to energy due to flux-advection step
  }

  // but we end up modifying pb, not pf.
  FTYPE uub[NPR]; // holds returned uu from implicit solver
  // for implicit scheme, set pf->pb in case pf contains final time information like the magnetic field
  // at this point, pb=pf was used to compute Flux(pb=pf), while uu0 is result of that flux.  Best for optically thin regime, but better guess limits used in actual solver.
  PLOOP(pliter,pl){
    pb[pl] = pf[pl]; // piin[pl];
    uub[pl] = uu0[pl]; //Uiin[pl];  // at this point, best approximation can make, but do better in solver itself
  }



  ////////////
  //
  // setup backup's for reversions of starting guess or solution
  //
  ///////////
  FTYPE pbbackup[NPR];
  FTYPE uubbackup[NPR];
  FTYPE dUcompbackup[NUMSOURCES][NPR];
  struct of_state qbackup;
  PLOOP(pliter,pl){
    pbbackup[pl]=pb[pl];  // at this point, best approximation can make, but do better in solver itself
    uubbackup[pl]=uub[pl]; // at this point, best approximation can make, but do better in solver itself
    SCLOOP(sc) dUcompbackup[sc][pl]=dUcomp[sc][pl]; // stores dUcomp so far computed from other physical terms
  }
  qbackup=*q;  // at this point, best approximation can make, but do better in solver itself


  ////////
  //
  // It's up to inversion method to set failure flags, not utoprimgen() that just checks them mostly (it might modify them based upon doing reductions).
  // setup pflags
  //
  ////////
  PFTYPE *lpflag,*lpflagrad;
  lpflag=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
  lpflagrad=&GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);
  // set default
  *lpflag=UTOPRIMNOFAIL;
  *lpflagrad=UTOPRIMRADNOFAIL;


  //////
  //
  // default is didn't get good primitives.  Similar, but slightly different from, failreturn
  //
  //////
  failreturn=1;
  noprims=1;
  usedenergy=0;
  usedentropy=0;
  usedboth=0;
  usedcold=0;
  usedimplicit=0;


  // whether doing energy at all
  int eomtypecond[NUMEOMTYPES];
  eomtypecond[EOMFFDE]=(*eomtype==EOMFFDE || *eomtype==EOMDEFAULT && EOMTYPE==EOMFFDE);
  eomtypecond[EOMCOLDGRMHD]=(*eomtype==EOMCOLDGRMHD || *eomtype==EOMDEFAULT && EOMTYPE==EOMCOLDGRMHD);
  eomtypecond[EOMENTROPYGRMHD]=(*eomtype==EOMENTROPYGRMHD || *eomtype==EOMDEFAULT && EOMTYPE==EOMENTROPYGRMHD);
  eomtypecond[EOMGRMHD]=(*eomtype==EOMGRMHD || *eomtype==EOMDEFAULT && EOMTYPE==EOMGRMHD);
  int testeom,counteom=0;
  for(testeom=0;testeom<NUMEOMTYPES;testeom++) if(eomtypecond[testeom]) counteom++;
  if(counteom>1){
    dualfprintf(fail_file,"Too many default eomtypeconds\n");
    for(testeom=0;testeom<NUMEOMTYPES;testeom++) dualfprintf(fail_file,"%d %d\n",testeom,eomtypecond[testeom]);
    myexit(843968364);
  }


  /////////
  //
  // diags
  //
  ////////
  FTYPE errorabs[NUMERRORTYPES]; // 0: over iterated pl and 1: over full relavant pl
  int iters=0;
  int f1iters=0;
  int nummhdinvs=0;
  int nummhdsteps=0;
  FTYPE fracenergy;
  int itermode;
  int baseitermethod;
  int radinvmod=0;
  int whichcap;
  FTYPE trueimptryconv=IMPTRYCONV;
  FTYPE trueimpokconv=IMPOKCONV;
  FTYPE trueimpallowconv=IMPALLOWCONVCONST;
  int trueimpmaxiter=IMPMAXITERLONG;
  int truenumdampattempts=NUMDAMPATTEMPTS;
  int goexplicit;







  //////////////////////////////
  //
  // MODEENERGY
  //
  //////////////////////////////

  if(MODEMETHOD==MODEENERGY && eomtypecond[EOMGRMHD]==1 || MODEMETHOD==MODEDEFAULT && *eomtype==EOMGRMHD){
    havebackup=0;
    didentropyalready=0;
    eomtypelocal=*eomtype;
    errorabs[0]=errorabs[1]=1.0;
    fracenergy=1.0;
    //itermode=ITERMODESTAGES;
    //baseitermethod=QTYPMHD;
    itermode=ITERMODENORMAL;
    baseitermethod=QTYURAD;

    // NOTES: ITERMODESTAGES with QTYPMHD always goes to damp and sometimes fails with RADTUBE.  Lots of Jsub issues, etc.  Probably not right.  Compared to ramesh code, settles on different u -- far too large even though error in f1[0] small.
    // Actually, if don't switch to PRAD, goes kinda ok with only 2 early failures, but needs to damp every time.  Can't be right.  So 2 issues.  Very broad iteration tail.

    // but, ITERMODESTAGES with QTYURAD does fine, even if takes more iterations.
    // but, ITERMODESTAGES with QTYPRAD does kinda ok, 6 early bads and ~20 damps, but then settles.

    whichcap=CAPTYPEBASIC;
    failreturn=koral_source_rad_implicit_mode(1,0,havebackup, didentropyalready, &eomtypelocal, whichcap, itermode, &baseitermethod, trueimptryconv, trueimpokconv, trueimpallowconv, trueimpmaxiter,  truenumdampattempts, fracenergy, dissmeasure, &radinvmod, pb, uub, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother ,dUcomp, errorabs, errorabs, &iters, &f1iters, &nummhdinvs, &nummhdsteps);
    if(ACTUALHARDFAILURE(failreturn)){
      failfinalreturn=1;
      *lpflag=UTOPRIMFAILCONV;
      *lpflagrad=UTOPRIMRADFAILCASE1A;
      // restore backups in case got contaminated
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        uub[pl]=uubbackup[pl];
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
    }
  }

  //////////////////////////////
  //
  // MODEENTROPY
  //
  //////////////////////////////

  if(MODEMETHOD==MODEENTROPY || MODEMETHOD==MODEDEFAULT && *eomtype==EOMENTROPYGRMHD){
    havebackup=0;
    didentropyalready=0;
    eomtypelocal=EOMENTROPYGRMHD;
    errorabs[0]=errorabs[1]=1.0;
    fracenergy=0.0;
    itermode=ITERMODESTAGES;
    whichcap=CAPTYPEBASIC;
    baseitermethod=QTYPMHD;
    failreturn=koral_source_rad_implicit_mode(1,0,havebackup, didentropyalready, &eomtypelocal, whichcap, itermode, &baseitermethod, trueimptryconv, trueimpokconv, trueimpallowconv, trueimpmaxiter,  truenumdampattempts, fracenergy, dissmeasure, &radinvmod, pb, uub, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother ,dUcomp, errorabs, errorabs, &iters, &f1iters, &nummhdinvs, &nummhdsteps);
    if(ACTUALHARDFAILURE(failreturn)){
      failfinalreturn=1;
      *lpflag=UTOPRIMFAILCONV;
      *lpflagrad=UTOPRIMRADFAILCASE1A;
      // restore backups in case got contaminated
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        uub[pl]=uubbackup[pl];
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
      usedentropy=1;
      goexplicit=0;
      usedimplicit=1;
    }
    else{
      failfinalreturn=-1; // indicates to koral_source_rad() that implicit says just do trivial explicit
      noprims=1;
      *eomtype=EOMENTROPYGRMHD;
      goexplicit=1;
      usedexplicitgood=1;
    }
  }





  if(MODEMETHOD==MODEENERGYRAMESH&&USERAMESH){
    // these ramesh solutions are here so both entropy and energy can have chance to fill them.
    int gotrameshsolution=0;
    FTYPE ppeng[NPR]={0},ppent[NPR]={0},errorabseng[NUMERRORTYPES],errorabsent[NUMERRORTYPES];
    set_array(errorabseng,NUMERRORTYPES,MPI_FTYPE,1.0);
    set_array(errorabsent,NUMERRORTYPES,MPI_FTYPE,1.0);
    int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITERLONG,itersent=IMPMAXITERLONG,radinvmodeng=-1,radinvmodent=-1;
    struct of_state qeng, qent;
    FTYPE uueng[NPR]={0.0},uuent[NPR]={0.0};
    qeng=qbackup;
    qent=qbackup;
    FTYPE dUcompeng[NUMSOURCES][NPR],dUcompent[NUMSOURCES][NPR];


    goexplicit=0; // force since no explicit check
    FTYPE errorabsforramesh[NUMERRORTYPES];
    set_array(errorabsforramesh,NUMERRORTYPES,MPI_FTYPE,1.0);

    PLOOP(pliter,pl){
      pb[pl]=pbbackup[pl];
      uub[pl]=uubbackup[pl];
      SCLOOP(sc) dUcompeng[sc][pl]=dUcomp[sc][pl]=dUcompbackup[sc][pl];
    }
    *q=qbackup;
    //
    // BEGIN GET RAMESH SOLUTION
    failreturn=FAILRETURNGENERAL;// default to fail
    *eomtype=EOMGRMHD;
    int whichcall=*eomtype;
    get_rameshsolution_wrapper(whichcall, *eomtype, errorabsforramesh, ptrgeom, pb, piin, Uiin, Ufin, dUother, CUf, q, ppeng, ppent, uueng, uuent, dUcompeng, dUcompent, &qeng, &qent, &failtypeeng, errorabseng, &iterseng, &radinvmodeng, &failtypeent, errorabsent, &itersent, &radinvmodent);
    gotrameshsolution=1; // indicates did at least attempt ramesh soltion call
    // translate
    PLOOP(pliter,pl){
      pb[pl]=ppeng[pl];
      uub[pl]=uueng[pl];
      SCLOOP(sc) dUcomp[sc][pl]=dUcompeng[sc][pl];
    }
    *q=qeng;
    //
    *lpflag=*lpflagrad=(PFTYPE)failtypeeng; // need better translation
    radinvmod=radinvmodeng;
    //    radErfneg=0; // not allowed, considered BADNEG type
    iters=iterseng;
    errorabs[0]=errorabseng[0];
    errorabs[1]=errorabseng[1];
    if(failtypeeng || errorabs[WHICHERROR]>IMPALLOWCONVCONSTABS){
      failfinalreturn=1;
      failreturn=FAILRETURNGENERAL; *eomtype=EOMGRMHD;
      // restore backups in case got contaminated
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        uub[pl]=uubbackup[pl];
        SCLOOP(sc) dUcomp[sc][pl]=dUcompbackup[sc][pl];
      }
      *q=qbackup;
      goexplicit=0;
      usedimplicit=1;
    }
    else if(errorabs[WHICHERROR]<IMPTRYCONVABS){
      failfinalreturn=0;
      failreturn=FAILRETURNNOFAIL;
      noprims=0;
      *eomtype=EOMDIDGRMHD;
      usedenergy=1;
      goexplicit=0;
      usedimplicit=1;
    }
    else{
      failreturn=FAILRETURNNOTTOLERROR; *eomtype=EOMDIDGRMHD;
      failfinalreturn=1;
      noprims=1;
      *eomtype=EOMGRMHD;
      usedenergy=1;
      goexplicit=0;
      usedimplicit=1;
    }
    // END GET RAMESH SOLUTION
  }




  //////////////////////////////
  //
  // MODESWITCH
  //
  //////////////////////////////

  if(MODEMETHOD==MODESWITCH && eomtypecond[EOMGRMHD]==1){
    FTYPE errorabsenergy[NUMERRORTYPES],errorabsentropy[NUMERRORTYPES];
    set_array(errorabsenergy,NUMERRORTYPES,MPI_FTYPE,1.0);
    set_array(errorabsentropy,NUMERRORTYPES,MPI_FTYPE,1.0);

    int itersenergy=0,itersentropy=0;
    int f1itersenergy=0,f1itersentropy=0;
    int nummhdinvsenergy=0,nummhdinvsentropy=0;
    int nummhdstepsenergy=0,nummhdstepsentropy=0;
    int itermodeenergy,itermodeentropy;
    int baseitermethodenergy,baseitermethodentropy;
    int whichcapenergy,whichcapentropy;
    FTYPE trueimptryconvenergy,trueimptryconventropy;
    FTYPE trueimpokconvenergy,trueimpokconventropy;
    FTYPE trueimpallowconvenergy,trueimpallowconventropy;
    int trueimpmaxiterenergy,trueimpmaxiterentropy;
    int truenumdampattemptsenergy,truenumdampattemptsentropy;
    int failreturnenergy=FAILRETURNGENERAL; // default
    int failreturnentropy=FAILRETURNGENERAL; // default

    eomtypelocal=*eomtype;
    havebackup=1; // only time this is used is here where we tell energy that we have backup method, so can give up quickly.
    didentropyalready=0;
    fracenergy=1.0;
    errorabsenergy[0]=errorabsenergy[1]=1.0;
    itermodeenergy=ITERMODESTAGES;
    whichcapenergy=CAPTYPEBASIC;
    trueimptryconvenergy=IMPTRYCONV;
    trueimpokconvenergy=IMPOKCONVCONST;
    trueimpallowconvenergy=IMPALLOWCONVCONST;
    trueimpmaxiterenergy=IMPMAXITERLONG;
    truenumdampattemptsenergy=NUMDAMPATTEMPTS;
    baseitermethodenergy=QTYPMHD;
    failreturnenergy=koral_source_rad_implicit_mode(1,0,havebackup, didentropyalready, &eomtypelocal, whichcapenergy, itermodeenergy, &baseitermethodenergy, trueimptryconvenergy, trueimpokconvenergy, trueimpallowconvenergy, trueimpmaxiterenergy,  truenumdampattemptsenergy, fracenergy, dissmeasure, &radinvmod, pb, uub, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother ,dUcomp, errorabsenergy, errorabsenergy, &itersenergy, &f1itersenergy, &nummhdinvsenergy, &nummhdstepsenergy);

    if(SWITCHGOODIDEAFAILURE(failreturn) && eomtypecond[EOMGRMHD]){
      // if failed with GRMHD or return reported switching to entropy is preferred, then do entropy method
      eomtypelocal=EOMENTROPYGRMHD;
      // restore backups for fresh start
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        uub[pl]=uubbackup[pl];
        SCLOOP(sc) dUcomp[sc][pl]=dUcompbackup[sc][pl];
      }
      *q=qbackup;
      // get fresh start entropy solution
      havebackup=0;
      didentropyalready=0;
      fracenergy=0.0;
      errorabsentropy[0]=errorabsentropy[1]=1.0;
      itermodeentropy=ITERMODESTAGES;
      whichcapentropy=CAPTYPEBASIC;
      trueimptryconventropy=IMPTRYCONV;
      trueimpokconventropy=IMPOKCONVCONST;
      trueimpallowconventropy=IMPALLOWCONVCONST;
      trueimpmaxiterentropy=IMPMAXITERLONG;
      truenumdampattemptsentropy=NUMDAMPATTEMPTS;
      baseitermethodentropy=QTYPMHD;
      failreturnentropy=koral_source_rad_implicit_mode(1,0,havebackup, didentropyalready, &eomtypelocal, whichcapentropy, itermodeentropy, &baseitermethodentropy, trueimptryconventropy, trueimpokconventropy, trueimpallowconventropy, trueimpmaxiterentropy,  truenumdampattemptsentropy, fracenergy, dissmeasure, &radinvmod, pb, uub, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother ,dUcomp, errorabsentropy, errorabsentropy, &itersentropy, &f1itersentropy, &nummhdinvsentropy, &nummhdstepsentropy);
      if(ACTUALHARDFAILURE(failreturnentropy)){
        failfinalreturn=1;
        *lpflag=UTOPRIMFAILCONV;
        *lpflagrad=UTOPRIMRADFAILCASE1A;
        goexplicit=0;
        usedimplicit=1;
        if(debugfail>=2) dualfprintf(fail_file,"Entropy also failed: energy=%d entropy=%d\n",failreturnenergy,failreturnentropy);
      }
      else if(failreturnentropy>=0){
        // use entropy
        failfinalreturn=0;
        usedentropy=1;
        noprims=0;
        errorabs[0]=errorabsentropy[0];
        errorabs[1]=errorabsentropy[1];
        iters=itersentropy;
        f1iters=f1itersentropy;
        // tell an externals to switch to entropy
        //*eomtype=EOMENTROPYGRMHD;
        *eomtype=eomtypelocal; // EOMDONOTHING if successful call to koral_source_rad_implicit_mode()
        goexplicit=0;
        usedimplicit=1;
      }
      else{
        failfinalreturn=-1; // indicates to koral_source_rad() that implicit says just do trivial explicit
        noprims=1;
        *eomtype=EOMDEFAULT;
        iters=itersentropy;
        f1iters=f1itersentropy;
        goexplicit=1;
        usedexplicitgood=1;
      }
    }
    else if(failreturnenergy>=0){// failreturnenergy==0 or eomtypecond[EOMGRMHD]==0
      // can stick with energy
      failfinalreturn=0;
      usedenergy=1;
      noprims=0;
      errorabs[0]=errorabsenergy[0];
      errorabs[1]=errorabsenergy[1];
      iters=itersenergy;
      f1iters=f1itersenergy;
      // switch to whatever solver suggested if didn't meet go-entropy condition
      *eomtype=eomtypelocal; // can also be EOMDONOTHING if successful and good enough error.
      goexplicit=0;
      usedimplicit=1;
      if(ACTUALHARDFAILURE(failreturnenergy) && debugfail>=2) dualfprintf(fail_file,"Decided didn't meet go-entropy condition but failed: failreturn=%d eomtypelocal=%d\n",failreturn,eomtypelocal);
    }
    else{
      failfinalreturn=-1; // indicates to koral_source_rad() that implicit says just do trivial explicit
      noprims=1;
      *eomtype=EOMDEFAULT;
      iters=itersentropy;
      f1iters=f1itersentropy;
      goexplicit=1;
      usedexplicitgood=1;
    }
  }





// KORALTODO: SUPERGODMARK: run normal koral tests.  Figure out all entropy signs.
// KORALTODO: SUPERGODMARK: Need backup to entropy since really dies if no backup.  E.g. cold backup.   But maybe using high accurate cold bad compared to lower accuracy entropy or lower accuracy energy.  Not sure should always prefer entropy if didn't reach desired tolerance.  But, currently if allowed tolerance, treated as ok solution and not failure to reject.  So this issue is ok relative to chosen IMPALLOWCONV.







  //////////////////////////////
  //
  // MODEPICKBEST / MODEPICKBESTSIMPLE / MODEPICKBESTSIMPLE2
  //
  // KORALTODO: SUPERGODMARK: Need to still think of energy-entropy implicit solver for intermediate regime.  Also, need to add conditions to avoid entropy or energy if Erf<0 or gammarad reaches limit or entropy-based effective absorption opacity is too large."
  //
  //////////////////////////////
  // 1) start with entropy
  // 2) Then if entropy didn't fail, use as guess for energy.
  //    During energy iteration, stop if repeatedly have u_g[entropy]>2*u_g[energy] using entropy guess stored for reference.
  // 3) If energy failed, stick with non-failed entropy.  If both failed, revert to failure modes or G=0
  //    If both succeeded, use entropy if u_g[entropy]>2*u_g[energy]

  // NEW ORDER (because entropy dies too much in shocks and takes too many trials)
  // 1) Try energy with QTYPMHD, QTYURAD, QTYPRAD.  Try quickly with no damping, ITERMODENORMAL.  Only try URAD/PRAD if PMHD fails or gives radinv=1
  // 2) Following #1, but then use ITERMODESTAGES, damping, higher uu if applicable.

#define NUMPHASES (6)
#define NUMPHASESENT (8)
#define NUMPHASESCOLD (1)

  // counters for which method was *attempted* even if not used
  static long long int tryphaselistenergy[NUMPHASES]={0};
  static long long int tryphaselistentropy[NUMPHASESENT]={0};
  static long long int tryphaselistcold[NUMPHASESCOLD]={0};


  int gotrameshsolution=0,usedrameshenergy=0,usedrameshentropy=0;
  if(MODEMETHOD==MODEPICKBEST || MODEMETHOD==MODEPICKBESTSIMPLE  || MODEMETHOD==MODEPICKBESTSIMPLE2){



    /////////
    //
    // ENTROPY VARIABLES
    //
    ////////

    // choices per attempt
    int itermodeentropy=ITERMODENORMAL; // start with normal
    int baseitermethodentropy=QTYPMHD; // default
    FTYPE trueimptryconventropy;
    FTYPE trueimpokconventropy;
    FTYPE trueimpallowconventropy;
    int trueimpmaxiterentropy;
    int truenumdampattemptsentropy;
    int whichcapentropy;

    // hold error and iterations from last attempt
    FTYPE errorabsentropyold[NUMERRORTYPES];
    set_array(errorabsentropyold,NUMERRORTYPES,MPI_FTYPE,1.0);

    int radinvmodentropyold=UTOPRIMRADFAILBAD1;
    int itersentropyold;
    // hold iterations for total entropy attempts
    int itersentropy=0;
    int f1itersentropy=0;
    int nummhdinvsentropy=0;
    int nummhdstepsentropy=0;


    // set that no entropy backup yet.
    havebackup=0;
    didentropyalready=0;
    fracenergy=0.0;

    // default in case no best solution with no error below 1.0 that is set as default best

    // holds best entropy solution
    FTYPE pbentropybest[NPR]; PLOOP(pliter,pl) pbentropybest[pl]=pbbackup[pl];
    FTYPE uubentropybest[NPR]; PLOOP(pliter,pl) uubentropybest[pl]=uubbackup[pl];
    FTYPE dUcompentropybest[NUMSOURCES][NPR]; PLOOP(pliter,pl) SCLOOP(sc) dUcompentropybest[sc][pl]=dUcompbackup[sc][pl];
    struct of_state qentropybest=qbackup;
    PFTYPE lpflagentropybest=1;
    PFTYPE lpflagradentropybest=1;
    int radinvmodentropybest=UTOPRIMRADFAILBAD1;
    int radErfnegentropybest=1;
    int failreturnentropybest=FAILRETURNGENERAL;
    int eomtypeentropybest=EOMENTROPYGRMHD;
    FTYPE errorabsentropybest[NUMERRORTYPES];
    set_array(errorabsentropybest,NUMERRORTYPES,MPI_FTYPE,1.0);
    int goexplicitentropybest=0;

    // holds latest entropy solution
    FTYPE pbentropy[NPR];
    FTYPE uubentropy[NPR]; // holds returned uu from implicit solver
    FTYPE dUcompentropy[NUMSOURCES][NPR]={{BIG}};
    struct of_state qentropy=qbackup;
    PFTYPE lpflagentropy=1;
    PFTYPE lpflagradentropy=1;
    int radinvmodentropy=UTOPRIMRADFAILBAD1;
    int radErfnegentropy=1;
    int failreturnentropy=FAILRETURNGENERAL;  // default to fail in case energy not to be done at all
    int eomtypeentropy=EOMENTROPYGRMHD;
    FTYPE errorabsentropy[NUMERRORTYPES]; // default is high error
    set_array(errorabsentropy,NUMERRORTYPES,MPI_FTYPE,1.0);
    int goexplicitentropy=0;





   
    ////////////////
    //
    // set fracenergy
    //
    ////////////////
    set_fracenergy(ptrgeom->i,ptrgeom->j,ptrgeom->k,dissmeasure, &fracenergy);


    /////////
    //
    // ENERGY VARIABLES
    //
    ////////

    // hold old error and iterations
    FTYPE errorabsenergyold[NUMERRORTYPES];
    set_array(errorabsenergyold,NUMERRORTYPES,MPI_FTYPE,1.0);
    int radinvmodenergyold=UTOPRIMRADFAILBAD1;
    int itersenergyold;

    // latest energy iterations
    int itersenergy=0;
    int f1itersenergy=0;
    int nummhdinvsenergy=0;
    int nummhdstepsenergy=0;

    // setting for each attempt
    int itermodeenergy=ITERMODENORMAL; // start with normal
    int baseitermethodenergy=QTYPMHD; // default
    FTYPE trueimptryconvenergy;
    FTYPE trueimpokconvenergy;
    FTYPE trueimpallowconvenergy;
    int trueimpmaxiterenergy;
    int truenumdampattemptsenergy;
    int whichcapenergy;


    // holds best result from implicit solver
    FTYPE pbenergybest[NPR]; PLOOP(pliter,pl) pbenergybest[pl]=pbbackup[pl];
    FTYPE uubenergybest[NPR]; PLOOP(pliter,pl) uubenergybest[pl]=uubbackup[pl];
    FTYPE dUcompenergybest[NUMSOURCES][NPR];  PLOOP(pliter,pl) SCLOOP(sc) dUcompenergybest[sc][pl]=dUcompbackup[sc][pl];
    struct of_state qenergybest=qbackup;
    // default is best is fail situation
    PFTYPE lpflagenergybest=1;
    PFTYPE lpflagradenergybest=1;
    int radinvmodenergybest=UTOPRIMRADFAILBAD1;
    int radErfnegenergybest=1;
    int failreturnenergybest=FAILRETURNGENERAL;
    int eomtypeenergybest=EOMGRMHD;
    FTYPE errorabsenergybest[NUMERRORTYPES];
    set_array(errorabsenergybest,NUMERRORTYPES,MPI_FTYPE,1.0);
    int goexplicitenergybest=0;

    // holds result from implicit solver
    FTYPE pbenergy[NPR];
    FTYPE uubenergy[NPR];
    FTYPE dUcompenergy[NUMSOURCES][NPR];
    set_array(dUcompenergy,NUMSOURCES*NPR,MPI_FTYPE,BIG);
    struct of_state qenergy=qbackup;
    PFTYPE lpflagenergy=1;
    PFTYPE lpflagradenergy=1;
    int radinvmodenergy=UTOPRIMRADFAILBAD1;
    int radErfnegenergy=1;
    int failreturnenergy=FAILRETURNGENERAL; // default
    int eomtypeenergy=EOMGRMHD;
    FTYPE errorabsenergy[NUMERRORTYPES];
    set_array(errorabsenergy,NUMERRORTYPES,MPI_FTYPE,1.0);
    int goexplicitenergy=0;



    /////////
    //
    // RAMESH VARIABLES AND DEFAULTS
    //
    // these ramesh solutions are here so both entropy and energy can have chance to fill them.
    //
    /////////
    FTYPE ppeng[NPR],ppent[NPR];
    FTYPE uueng[NPR],uuent[NPR];
    FTYPE dUcompeng[NUMSOURCES][NPR],dUcompent[NUMSOURCES][NPR];
    PLOOP(pliter,pl){
      ppeng[pl]=ppent[pl]=pbbackup[pl];
      uueng[pl]=uuent[pl]=uubbackup[pl];
      SCLOOP(sc) dUcompeng[sc][pl]=dUcompent[sc][pl]=dUcompbackup[sc][pl];
    }
    struct of_state qeng=qbackup, qent=qbackup;
    FTYPE errorabseng[NUMERRORTYPES],errorabsent[NUMERRORTYPES];
    set_array(errorabseng,NUMERRORTYPES,MPI_FTYPE,1.0);
    set_array(errorabsent,NUMERRORTYPES,MPI_FTYPE,1.0);
    int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITERLONG,itersent=IMPMAXITERLONG,radinvmodeng=UTOPRIMRADFAILBAD1,radinvmodent=UTOPRIMRADFAILBAD1;





    //////
    //
    // check if should re-order method attempt list in case where primary evolving quantity is radiation.
    //
    // assume at this stage that radiation primarily evolves if Erf sufficiently large compared to u_g or changes in URAD0 are sufficiently large
    //
    // check
    static FTYPE sqrtnumepsilon1;
    static FTYPE sqrtnumepsilon2;
    static int firsttimeset=1;
    if(firsttimeset){
      //      sqrtnumepsilon=10.0*pow(NUMEPSILON,1.0/3.0);
      //      sqrtnumepsilon=1E-1; // playing -- required for RADBONDI to be fast by using QTYURAD method first as QTYPMHD method fails more.
      sqrtnumepsilon1=MIN(1.0,10.0*NUMEPSILON/IMPTRYCONV);
      sqrtnumepsilon2=pow(NUMEPSILON,1.0/3.0);
      firsttimeset=0;
    }
    int radprimaryevolves=0;
    if(fabs(rdUtot[UU])<sqrtnumepsilon1*fabs(rdUtot[URAD0]) || fabs(uu0[URAD0])<sqrtnumepsilon1*fabs(uu0[UU]) || fabs(pb[URAD0])<sqrtnumepsilon1*fabs(pb[UU])){
      radprimaryevolves=1;
    }
    else{
      radprimaryevolves=0;
    }
    int radextremeprimaryevolves=0;
    if(fabs(rdUtot[UU])<10.0*NUMEPSILON*fabs(dUtot[URAD0]) || fabs(uu0[URAD0])<10.0*NUMEPSILON*fabs(uu0[UU]) || fabs(pb[URAD0])<10.0*NUMEPSILON*fabs(pb[UU])){
      radextremeprimaryevolves=1;
    }
    else{
      radextremeprimaryevolves=0;
    }
    int gasprimaryevolves=0;
    //    if(radprimaryevolves==0 && (pb[URAD0]<sqrtnumepsilon2*pb[UU] || (-uu0[URAD0])<sqrtnumepsilon2*(-uu0[UU]))){
    //    if(pb[URAD0]<sqrtnumepsilon2*pb[UU] || fabs(dUtot[URAD0])<sqrtnumepsilon2*fabs(dUtot[UU])){
    if(radprimaryevolves==0 && (fabs(rdUtot[URAD0])<sqrtnumepsilon2*fabs(rdUtot[UU]) || fabs(uu0[UU])<sqrtnumepsilon2*fabs(uu0[URAD0]) || fabs(pb[UU])<sqrtnumepsilon2*fabs(pb[URAD0]) )){
      gasprimaryevolves=1;
    }
    else{
      gasprimaryevolves=0;
    }
    int gasextremeprimaryevolves=0;
    if(radprimaryevolves==0 && radextremeprimaryevolves==0 && (fabs(rdUtot[URAD0])<10.0*NUMEPSILON*fabs(dUtot[UU]) || fabs(uu0[UU])<10.0*NUMEPSILON*fabs(uu0[URAD0]) || fabs(pb[UU])<10.0*NUMEPSILON*fabs(pb[URAD0]) )){
      gasextremeprimaryevolves=1;
    }
    else{
      gasextremeprimaryevolves=0;
    }

    if(LETPMHDFAIL){
      // just let pmhd fail or use total error in way that is setup now...
      // still allow extreme catches
      radprimaryevolves=0;
      gasprimaryevolves=0;
    }

    if(MODEMETHOD==MODEPICKBESTSIMPLE){
      // forcing PMHD method, so must use gas
      radprimaryevolves=radextremeprimaryevolves=0;
      gasprimaryevolves=gasextremeprimaryevolves=1;
    }






    // DEBUG:
    //    dualfprintf(fail_file,"PRIMARYEVOLVES: %d %d %d %d : pb=%g %g uu0=%g %g dUtot=%g %g : sqrtnumepsilon1=%g sqrtnumepsilon2=%g\n",radprimaryevolves,radextremeprimaryevolves,gasprimaryevolves,gasextremeprimaryevolves,pb[UU],pb[URAD0],-uu0[UU],-uu0[URAD0],dUtot[UU],dUtot[URAD0],sqrtnumepsilon1,sqrtnumepsilon2);



    /////////
    //
    // store prior itermethod, before we change it.
    //
    /////////
    int prioritermethodlist[NUMPRIORITERMETHODINDEX];
    int oo;
    for(oo=0;oo<NUMPRIORITERMETHODINDEX;oo++){
      prioritermethodlist[oo]=GLOBALMACP0A1(prioritermethod,i,j,k,oo);
    }



    /////////////
    //
    // GET ENERGY
    //
    // Now do energy if would use only energy or if doing interpolation
    // also check energy if entropy thinks we should do explicit
    //
    /////////////
    //    if(prioritermethodlist[EOMINDEX] == EOMGRMHD || prioritermethodlist[EOMINDEX] == PRIORITERMETHODNOTSET && eomtypecond[EOMGRMHD] && (fracenergy!=0.0 || RADINVBAD(radinvmodentropy) ||  ACTUALHARDORSOFTFAILURE(failreturnentropy)==1 || goexplicitentropy==1)){
    if(eomtypecond[EOMGRMHD] && (fracenergy!=0.0 || RADINVBAD(radinvmodentropy) ||  ACTUALHARDORSOFTFAILURE(failreturnentropy)==1 || goexplicitentropy==1)){


      // quickly try QTYPMHD then QTYURAD
      int tryphase1;

      int baseitermethodlist[NUMPHASES]={QTYPMHD,QTYURAD,QTYPRAD,QTYPMHD,QTYURAD,QTYPRAD}; int whichfirstpmhd=0,whichfirsturad=1,whichfirstprad=2;
      int itermodelist[NUMPHASES]={ITERMODENORMAL,ITERMODENORMAL,ITERMODENORMAL,ITERMODESTAGES,ITERMODESTAGES,ITERMODESTAGES};
      FTYPE trueimptryconvlist[NUMPHASES]={IMPTRYCONV,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK};
      FTYPE trueimpokconvlist[NUMPHASES]={IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK};
      FTYPE trueimpallowconvconstlist[NUMPHASES]={IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST};
      int trueimpmaxiterlist[NUMPHASES]={IMPMAXITERQUICK,IMPMAXITERQUICK,IMPMAXITERQUICK,IMPMAXITERMEDIUM,IMPMAXITERMEDIUM,IMPMAXITERMEDIUM};
      int truenumdampattemptslist[NUMPHASES]={NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK};
      int modprimlist[NUMPHASES]={0,0,0,1,0,0};
      int checkradinvlist[NUMPHASES]={0,1,1,0,1,1};
      int eomtypelist[NUMPHASES]={EOMGRMHD,EOMGRMHD,EOMGRMHD,EOMGRMHD,EOMGRMHD,EOMGRMHD};

      // results in list
      FTYPE errorabslist[NUMPHASES][NUMERRORTYPES];
      set_array(errorabslist,NUMPHASES*NUMERRORTYPES,MPI_FTYPE,1.0);
      int radinvmodlist[NUMPHASES];
      set_array(radinvmodlist,NUMPHASES,MPI_INT,UTOPRIMRADFAILBAD1);

      // reorder method if desired
      // KORALNOTE: radinv check would nominally catch if PMHD method failed due to machine errors in GAS leading to huge changes in RAD leading to E_r<0 or gamma>gammaradmax, but might as well use desired method first.
      // normal checkradinvlist will catch if prad method doesn't lead to radinv but urad method does.
      // normal checkradinvlist will catch if somehow rad method leads to radinv but gas method doesn't.
      if(gasprimaryevolves==0){ // if gas primary evolves, then keep original order because PMHD is fastest method.
        if(radprimaryevolves){
          if(
             (USEPRIORITERMETHOD && (prioritermethodlist[BASEITERMETHODINDEX] == PRIORITERMETHODNOTSET ||  prioritermethodlist[BASEITERMETHODINDEX] == QTYURAD))
             || USEPRIORITERMETHOD==0
             ){
            tryphase1=-1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD; modprimlist[tryphase1]=0; whichfirsturad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD; modprimlist[tryphase1]=0; whichfirstprad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD; modprimlist[tryphase1]=0; whichfirstpmhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD; modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD; modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD; modprimlist[tryphase1]=1;
          }
          if(
             USEPRIORITERMETHOD && prioritermethodlist[BASEITERMETHODINDEX] == QTYPRAD
             ){
            tryphase1=-1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD; modprimlist[tryphase1]=0; whichfirstprad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD; modprimlist[tryphase1]=0; whichfirsturad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD; modprimlist[tryphase1]=0; whichfirstpmhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD; modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD; modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD; modprimlist[tryphase1]=1;
          }
        }// if rad is expected to need to be evolved
      }

      




      // ENERGY PHASE LOOP
      int firsttryphase1used=-1;
      for(tryphase1=0;tryphase1<NUMPHASES;tryphase1++){

        //        dualfprintf(fail_file,"TRYING: tryphase1=%d : %d %d %d\n",tryphase1,gasextremeprimaryevolves,radextremeprimaryevolves,baseitermethodlist[tryphase1]);

        // pick best simple method avoids all solvers except PMHD
        if(MODEMETHOD==MODEPICKBESTSIMPLE && baseitermethodlist[tryphase1]!=QTYPMHD) continue;

        // pick best simple 2 method avoids all itermodestages methods
        if(MODEMETHOD==MODEPICKBESTSIMPLE2 && itermodelist[tryphase1]==ITERMODESTAGES) continue;

        // avoid method in case very non-dominant since then would give errorneous (critically bad even) results.  If methods that can use fail, have to revert to entropy or fixups.
        if(radextremeprimaryevolves && baseitermethodlist[tryphase1]==QTYPMHD) continue;
        if(gasextremeprimaryevolves && (baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD)) continue;
        
        if(radextremeprimaryevolves==0 && (baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD) ){
          // if not in extreme radiative regime, then don't do damping for raditive schemes that are slow.
          truenumdampattemptslist[tryphase1]=NUMDAMPATTEMPTSQUICK;
          // and don't try to get too good of error.
          trueimptryconvlist[tryphase1]=IMPTRYCONVQUICK;
          // and avoid stages
          itermodelist[tryphase1]=ITERMODENORMAL;
        }
        if(gasextremeprimaryevolves==0 && (baseitermethodlist[tryphase1]==QTYPMHD)){
          // then still try to get good error since fast method, so no changes.
        }


        if(MODEMETHOD==MODEPICKBESTSIMPLE){
          // do nothing, don't skip stages.  Even when iteration-error is small, often stages can get both iteration and total errors to be small.
        }
        else{
          // If already tried QTYPMHD and that got error(0)<tol but error(1)>tol, then skip QTYPMHD with stages since should switch to QTYURAD or QTYPRAD because STAGES won't improve on error(1) if error(0)<tol.
          // independent of radinvmod.  If ==0 or ==1, still should skip if error in that condition.
          if(baseitermethodlist[tryphase1]==QTYPMHD && errorabslist[whichfirstpmhd][0]<IMPTRYCONV && errorabslist[whichfirstpmhd][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon radiative solvers
            continue;
          }

          if(baseitermethodlist[tryphase1]==QTYURAD && errorabslist[whichfirsturad][0]<IMPTRYCONV && errorabslist[whichfirsturad][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon pmhd or prad
            continue;
          }

          if(baseitermethodlist[tryphase1]==QTYPRAD && errorabslist[whichfirstprad][0]<IMPTRYCONV && errorabslist[whichfirstprad][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon pmhd or urad
            continue;
          }

          //        if(prioritermethodlist[EOMTYPEINDEX]==EOMENTROPYGRMHD && itermodelist[tryphase1]==ITERMODESTAGES && errorabsenergybest[0]<IMPTRYCONV && errorabsenergybest[1]>IMPTRYCONV && WHICHERROR==1){
          //          // then assume went through all ITERMODENORMAL cases and best error had good iteration tolerance, but best error has bad total error, and previous solution was entropy, so skip itermodestages.
          //          // this assumes ITERMODESTAGES is expensive and should only do if ........
          //          continue;
          //        }
        }


        //        dualfprintf(fail_file,"MAYBEREALLYTRYING: tryphase1=%d : %d %d %d\n",tryphase1,radinvmodenergybest,checkradinvlist[tryphase1],failreturnenergybest);



        // FUCK
        //        if(baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD) continue; // skip this method for now.



        // consider radinvmod only if error bad for original approach.  Avoids excessive attempts when should hit radiative ceiling and error is small.
        // KORALTODO: KORALNOTE: If explicit was triggered (failreturnenergy) then could move on, but go ahead and test using other methods in case radinvmod!=0 can be avoided.
        //radinvmodenergybest!=0
        if(RADINVBAD(radinvmodenergybest) && checkradinvlist[tryphase1] || ACTUALHARDORSOFTFAILURE(failreturnenergybest) && failreturnenergybest!=FAILRETURNMODESWITCH){
          if(firsttryphase1used==-1) firsttryphase1used=tryphase1;
          tryphaselistenergy[tryphase1]++;

          //          dualfprintf(fail_file,"REALLYTRYING: tryphase1=%d\n",tryphase1);
     

          errorabsenergyold[0]=errorabsenergy[0];
          errorabsenergyold[1]=errorabsenergy[1];
          radinvmodenergyold=radinvmodenergy;
          itersenergyold=itersenergy;
          //
          whichcapenergy=CAPTYPEBASIC;
          //whichcapenergy=CAPTYPEBASIC;
          baseitermethodenergy=baseitermethodlist[tryphase1];
          itermodeenergy=itermodelist[tryphase1];
          trueimptryconvenergy=trueimptryconvlist[tryphase1];
          trueimpokconvenergy=(MAX(trueimptryconvlist[tryphase1],trueimpokconvlist[tryphase1]));
          trueimpallowconvenergy=(MAX(trueimptryconvlist[tryphase1],trueimpallowconvconstlist[tryphase1]));
          trueimpmaxiterenergy=trueimpmaxiterlist[tryphase1];
          truenumdampattemptsenergy=truenumdampattemptslist[tryphase1];
          int modprim=0;
          if(havebackup==0) modprim=modprimlist[tryphase1];
          //
          // start fresh
          *lpflag=UTOPRIMNOFAIL;
          *lpflagrad=UTOPRIMRADNOFAIL;
          radinvmodenergy=0;
          radErfnegenergy=0;
          failreturnenergy=FAILRETURNGENERAL; // default to fail in case energy not to be done at all
          eomtypeenergy=eomtypelist[tryphase1];
          errorabsenergy[0]=errorabsenergy[1]=1.0;
        
          // setup solver conditions for existence of backup solutions or entropy solution
          if(ACTUALHARDFAILURE(failreturnentropy)==1){
            havebackup=0; // no entropy solver solution, so no backup.
            didentropyalready=0;
          }
          else{
            havebackup=1; // here havebackup=1 means can break out of energy solver early if issues, since we do have entropy solver solution as backup.
            didentropyalready=1; // so can use entropy solution as reference for whether to stop energy iteration
          }

          // setup guess
          if(ACTUALHARDFAILURE(failreturnentropy)==1 || errorabsentropy[WHICHERROR]>ERRORTOUSEENTROPYFORENERGYGUESS){
            if(errorabsenergybest[WHICHERROR]>TRYHARDERFEEDGUESSTOL){
              PLOOP(pliter,pl){
                pbenergy[pl]=pbbackup[pl];
                uubenergy[pl]=uubbackup[pl];
                SCLOOP(sc) dUcompenergy[sc][pl]=dUcompbackup[sc][pl];
              }
              qenergy=qbackup;
              errorabsenergy[0]=errorabsenergy[1]=1.0;
              radinvmodenergy=UTOPRIMRADFAILBAD1;
            }
            else{
              PLOOP(pliter,pl){
                pbenergy[pl]=pbenergybest[pl];
                uubenergy[pl]=uubenergybest[pl];
                SCLOOP(sc) dUcompenergy[sc][pl]=dUcompbackup[sc][pl];
              }
              qenergy=qenergybest;
              errorabsenergy[0]=errorabsenergybest[0];
              errorabsenergy[1]=errorabsenergybest[1];
              radinvmodenergy=radinvmodenergybest;
            }
          }
          else{ // start with entropy
            PLOOP(pliter,pl){
              pbenergy[pl]=pbentropy[pl];
              uubenergy[pl]=uubentropy[pl];
              SCLOOP(sc) dUcompenergy[sc][pl]=dUcompbackup[sc][pl];
            }
            qenergy=qentropy;
            errorabsenergy[0]=errorabsentropy[0]; // assume error estimate only used to decide if guess is used, not if converged with energy
            errorabsenergy[1]=errorabsentropy[1]; // assume error estimate only used to decide if guess is used, not if converged with energy
            radinvmodenergy=radinvmodentropy;
          }
          //
          failreturnenergy=koral_source_rad_implicit_mode(0,modprim,havebackup, didentropyalready, &eomtypeenergy, whichcapenergy, itermodeenergy, &baseitermethodenergy, trueimptryconvenergy, trueimpokconvenergy, trueimpallowconvenergy, trueimpmaxiterenergy,  truenumdampattemptsenergy, fracenergy, dissmeasure, &radinvmodenergy, pbenergy, uubenergy, piin, Uiin, Ufin, CUf, ptrgeom, &qenergy, dUother ,dUcompenergy, errorabsenergy, errorabsenergybest, &itersenergy, &f1itersenergy, &nummhdinvsenergy, &nummhdstepsenergy);
          
          // store error, no matter if using solution or not or even if explicit
          errorabslist[tryphase1][0]=errorabsenergy[0];
          errorabslist[tryphase1][1]=errorabsenergy[1];

          if(failreturnenergy==FAILRETURNGOEXPLICIT){
            lpflagenergybest=*lpflag;
            lpflagradenergybest=*lpflagrad;
            radinvmodenergybest=radinvmodenergy;
            radErfnegenergybest=radErfnegenergy;
            failreturnenergybest=failreturnenergy;
            eomtypeenergybest=eomtypelist[tryphase1];
            goexplicitenergybest=1;

            // save which method ended up using
            int prioriter; for(prioriter=0;prioriter<NUMPRIORITERMETHODINDEX;prioriter++){
              GLOBALMACP0A1(prioritermethod,i,j,k,prioriter) = PRIORITERMETHODNOTSET; // back to unknown
            }
          }
          else{
            // see if want to keep
            //        if(errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && ACTUALHARDFAILURE(failreturnenergy)==0 || failreturnenergy==FAILRETURNMODESWITCH){
            //        if((errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && (RADINVBAD(radinvmodenergy)==0 || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)) || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)==0 && (errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR]||errorabsenergy[WHICHERROR]<IMPOKCONVABS)) && ACTUALHARDFAILURE(failreturnenergy)==0 || failreturnenergy==FAILRETURNMODESWITCH){
            if(
               (
                // best and no radinv problem
                errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && RADINVBAD(radinvmodenergy)==0
                // best because prior best also had radinv problem and still smaller error
                || errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && (RADINVBAD(radinvmodenergy) && RADINVBAD(radinvmodenergybest))
                // best because even though has radinv problem (while best doesn't) much smaller error
                || errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && errorabsenergy[WHICHERROR]<IMPOKCONVABS && errorabsenergybest[WHICHERROR]>100.0*IMPOKCONVABS && (RADINVBAD(radinvmodenergy) && RADINVBAD(radinvmodenergybest)==0)
                )
               && ACTUALHARDFAILURE(failreturnenergy)==0){
              // store result in case better than latter results
              lpflagenergybest=*lpflag;
              lpflagradenergybest=*lpflagrad;
              radinvmodenergybest=radinvmodenergy;
              radErfnegenergybest=radErfnegenergy;
              failreturnenergybest=failreturnenergy;
              eomtypeenergybest=eomtypeenergy;
              goexplicitenergybest=0;
              errorabsenergybest[0]=errorabsenergy[0];
              errorabsenergybest[1]=errorabsenergy[1];
              PLOOP(pliter,pl){
                pbenergybest[pl]=pbenergy[pl];
                uubenergybest[pl]=uubenergy[pl];
                SCLOOP(sc) dUcompenergybest[sc][pl]=dUcompenergy[sc][pl];
              }
              qenergybest=qenergy;

              // save which method ended up using
              GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) = baseitermethodlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) = itermodelist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,IMPTRYCONVINDEX) = (trueimptryconvlist[tryphase1]==IMPTRYCONV);
              GLOBALMACP0A1(prioritermethod,i,j,k,IMPMAXITERINDEX) = trueimpmaxiterlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,NUMDAMPINDEX) = truenumdampattemptslist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,MODPRIMINDEX) = modprimlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,CHECKRADINVINDEX) = checkradinvlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = eomtypelist[tryphase1];
            }
          }

          // still cost more iters
          itersenergy+=itersenergyold;
          
          int noproblem=(ACTUALHARDORSOFTFAILURE(failreturnenergy)==0 || failreturnenergy==FAILRETURNMODESWITCH);

          
          if(tryphase1!=firsttryphase1used && noproblem){
            if(debugfail>=2) dualfprintf(fail_file,"Recovered using tryphase1=%d (energy: failreturnenergy=%d radinvmod=%d -> %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g (best=%21.15g) iters: %d->%d baseitermethod=%d eomtype=%d\n",tryphase1,failreturnenergy,radinvmodenergyold,radinvmodenergy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsenergyold[WHICHERROR],errorabsenergy[WHICHERROR],errorabsenergybest[WHICHERROR],itersenergyold,itersenergy,baseitermethodenergy,eomtypeenergy);
          }
          if(noproblem==0){
            if(debugfail>=2) dualfprintf(fail_file,"Failed to: <Recovered> using tryphase1=%d (energy: failreturnenergy=%d radinvmod=%d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g (best=%21.15g) iters: %d->%d baseitermethod=%d eomtype=%d\n",tryphase1,failreturnenergy,radinvmodenergy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsenergyold[WHICHERROR],errorabsenergy[WHICHERROR],errorabsenergybest[WHICHERROR],itersenergyold,itersenergy,baseitermethodenergy,eomtypeenergy);
          }

         
        }// end if doing _mode() call


      }// end over try loop








      // see if want to try harder
      if(
         !(MODEMETHOD==MODEPICKBESTSIMPLE2) // Ramesh method is like stages, so avoid if simple2 method
         && (ACTUALHARDORSOFTFAILURE(failreturnenergy) && failreturn!=FAILRETURNMODESWITCH && USERAMESH && radextremeprimaryevolves==0)
         && !(errorabslist[whichfirstpmhd][0]<IMPTRYCONV && errorabslist[whichfirstpmhd][1]>IMPTRYCONV) // avoid ramesh (pmhd-type) method if already did it and got iterate-small error.
         ){
        errorabsenergyold[0]=errorabsenergy[0];
        errorabsenergyold[1]=errorabsenergy[1];
        radinvmodenergyold=radinvmodenergy;
        itersenergyold=itersenergy;
        goexplicitenergy=0; // force since no explicit check
        //
        if(gotrameshsolution){// then just recall solution
          PLOOP(pliter,pl){
            pbenergy[pl]=ppeng[pl];
            uubenergy[pl]=uueng[pl];
            SCLOOP(sc) dUcompenergy[sc][pl]=dUcompeng[sc][pl];
          }
          qenergy=qeng;
          //
          *lpflag=*lpflagrad=(PFTYPE)failtypeeng; // need better translation
          radinvmodenergy=radinvmodeng;
          radErfnegenergy=0; // not allowed, considered BADNEG type
          itersenergy=iterseng;
          errorabsenergy[0]=errorabseng[0];
          errorabsenergy[1]=errorabseng[1];
          if(failtypeeng || errorabsenergy[WHICHERROR]>IMPALLOWCONVCONSTABS){ failreturnenergy=FAILRETURNGENERAL; eomtypeenergy=EOMGRMHD; }
          else if(errorabsenergy[WHICHERROR]<IMPTRYCONVABS){ failreturnenergy=FAILRETURNNOFAIL; eomtypeenergy=EOMDIDGRMHD;}
          else{ failreturnenergy=FAILRETURNNOTTOLERROR; eomtypeenergy=EOMDIDGRMHD;}
        }
        else{ // else get new ramesh solution
          FTYPE errorabsforramesh[NUMERRORTYPES];
          set_array(errorabsforramesh,NUMERRORTYPES,MPI_FTYPE,1.0);
          // start fresh with ramesh
          if(errorabsenergy[WHICHERROR]>TRYHARDERFEEDGUESSTOL){
            if(ACTUALHARDFAILURE(failreturnentropy)==1 || errorabsentropy[WHICHERROR]>ERRORTOUSEENTROPYFORENERGYGUESS){
              PLOOP(pliter,pl){
                pbenergy[pl]=pbbackup[pl];
                uubenergy[pl]=uubbackup[pl];
                SCLOOP(sc) dUcompeng[sc][pl]=dUcompenergy[sc][pl]=dUcompbackup[sc][pl];
              }
              qenergy=qbackup;
              errorabseng[0]=errorabseng[1]=1.0;
              errorabsforramesh[0]=errorabsforramesh[1]=1.0;
              radinvmodeng=UTOPRIMRADFAILBAD1;
            }
            else{
              PLOOP(pliter,pl){
                pbenergy[pl]=pbentropy[pl];
                uubenergy[pl]=uubentropy[pl];
                SCLOOP(sc) dUcompeng[sc][pl]=dUcompenergy[sc][pl]=dUcompbackup[sc][pl];
              }
              qenergy=qentropy;
              errorabseng[0]=errorabsentropybest[0];
              errorabseng[1]=errorabsentropybest[1];
              errorabsforramesh[0]=errorabsforramesh[1]=1.0;
              radinvmodeng=radinvmodentropybest;
            }
          }
          else{
            PLOOP(pliter,pl){
              pbenergy[pl]=pbenergybest[pl];
              uubenergy[pl]=uubenergybest[pl];
              SCLOOP(sc) dUcompeng[sc][pl]=dUcompenergy[sc][pl]=dUcompbackup[sc][pl]; // always backup
            }
            qenergy=qenergybest;
            errorabseng[0]=errorabsenergybest[0];
            errorabseng[1]=errorabsenergybest[1];
            errorabsforramesh[0]=errorabsenergybest[0];
            errorabsforramesh[1]=errorabsenergybest[1];
            radinvmodeng=radinvmodenergybest;
          }
          //
          // BEGIN GET RAMESH SOLUTION
          failreturnenergy=FAILRETURNGENERAL;// default to fail
          eomtypeenergy=EOMGRMHD;
          int whichcall=eomtypeenergy;
          get_rameshsolution_wrapper(whichcall, eomtypeenergy, errorabsforramesh, ptrgeom, pbenergy, piin, Uiin, Ufin, dUother, CUf, q, ppeng, ppent, uueng, uuent, dUcompeng, dUcompent, &qeng, &qent, &failtypeeng, errorabseng, &iterseng, &radinvmodeng, &failtypeent, errorabsent, &itersent, &radinvmodent);
          gotrameshsolution=1; // indicates did at least attempt ramesh soltion call
          // translate
          PLOOP(pliter,pl){
            pbenergy[pl]=ppeng[pl];
            uubenergy[pl]=uueng[pl];
            SCLOOP(sc) dUcompenergy[sc][pl]=dUcompeng[sc][pl];
          }
          qenergy=qeng;
          //
          *lpflag=*lpflagrad=(PFTYPE)failtypeeng; // need better translation
          radinvmodenergy=radinvmodeng;
          radErfnegenergy=0; // not allowed, considered BADNEG type
          itersenergy=iterseng;
          errorabsenergy[0]=errorabseng[0];
          errorabsenergy[1]=errorabseng[1];
          if(failtypeeng || errorabsenergy[WHICHERROR]>IMPALLOWCONVCONSTABS){ failreturnenergy=FAILRETURNGENERAL; eomtypeenergy=EOMGRMHD; }
          else if(errorabsenergy[WHICHERROR]<IMPTRYCONVABS){ failreturnenergy=FAILRETURNNOFAIL; eomtypeenergy=EOMDIDGRMHD;}
          else{ failreturnenergy=FAILRETURNNOTTOLERROR; eomtypeenergy=EOMDIDGRMHD;}
          // END GET RAMESH SOLUTION
        }// end else if need to get ramesh solution

        // see if want to keep
        //        if(errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && ACTUALHARDFAILURE(failreturnenergy)==0){// || failreturnenergy==FAILRETURNMODESWITCH){ // no switch mode in ramesh solver yet.
        //        if((errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && (RADINVBAD(radinvmodenergy)==0 || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)) || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)==0 && (errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR]||errorabsenergy[WHICHERROR]<IMPOKCONVABS)) && ACTUALHARDFAILURE(failreturnenergy)==0){

        //        if((errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && (RADINVBAD(radinvmodenergy)==0 || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)) || RADINVBAD(radinvmodenergybest) && RADINVBAD(radinvmodenergy)==0 && (errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR]||errorabsenergy[WHICHERROR]<IMPOKCONVABS) || RADINVBAD(radinvmodenergybest)==0 && RADINVBAD(radinvmodenergy) && (errorabsenergybest[WHICHERROR]>IMPOKCONVABS && errorabsenergy[WHICHERROR]<IMPOKCONVABS)) && ACTUALHARDFAILURE(failreturnenergy)==0){

        if(
           (
            // best and no radinv problem
            errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && RADINVBAD(radinvmodenergy)==0
            // best because prior best also had radinv problem and still smaller error
            || errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && (RADINVBAD(radinvmodenergy) && RADINVBAD(radinvmodenergybest))
            // best because even though has radinv problem (while best doesn't) much smaller error
            || errorabsenergy[WHICHERROR]<errorabsenergybest[WHICHERROR] && errorabsenergy[WHICHERROR]<IMPOKCONVABS && errorabsenergybest[WHICHERROR]>100.0*IMPOKCONVABS && (RADINVBAD(radinvmodenergy) && RADINVBAD(radinvmodenergybest)==0)
            )
           && ACTUALHARDFAILURE(failreturnenergy)==0){

          if(ACCEPTASNOFAILURE(failreturnenergy)) usedrameshenergy=1; // means will use this actually, not just best yet no good enough
          // store result in case better than latter results
          lpflagenergybest=*lpflag;
          lpflagradenergybest=*lpflagrad;
          radinvmodenergybest=radinvmodenergy;
          radErfnegenergybest=radErfnegenergy;
          failreturnenergybest=failreturnenergy;
          eomtypeenergybest=eomtypeenergy;
          errorabsenergybest[0]=errorabsenergy[0];
          errorabsenergybest[1]=errorabsenergy[1];
          PLOOP(pliter,pl){
            pbenergybest[pl]=pbenergy[pl];
            uubenergybest[pl]=uubenergy[pl];
            SCLOOP(sc) dUcompenergybest[sc][pl]=dUcompenergy[sc][pl];
          }
          qenergybest=qenergy;

          // save which method ended up using
          GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) = QTYPMHD;
          GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) = ITERMODESTAGES;
          GLOBALMACP0A1(prioritermethod,i,j,k,IMPMAXITERINDEX) = IMPMAXITERLONG;
          GLOBALMACP0A1(prioritermethod,i,j,k,NUMDAMPINDEX) = NUMDAMPATTEMPTSQUICK;
          GLOBALMACP0A1(prioritermethod,i,j,k,MODPRIMINDEX) = 0;
          GLOBALMACP0A1(prioritermethod,i,j,k,CHECKRADINVINDEX) = 1;
          GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMGRMHD;

        }

        // still cost more iters
        itersenergy+=itersenergyold;

        if(ACTUALHARDORSOFTFAILURE(failreturnenergy)==0){// || failreturnenergy==FAILRETURNMODESWITCH){
          if(debugfail>=2) dualfprintf(fail_file,"(energy: failreturnenergy=%d): Recovered using ramesh (%d %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g iters: %d->%d\n",failreturnenergy,gotrameshsolution,usedrameshenergy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsenergyold[WHICHERROR],errorabsenergy[WHICHERROR],itersenergyold,itersenergy);
        }
        else{
          if(debugfail>=2) dualfprintf(fail_file,"Failed to  (energy: failreturnenergy=%d): Recovered using ramesh (%d %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g iters: %d->%d\n",failreturnenergy,gotrameshsolution,usedrameshenergy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsenergyold[WHICHERROR],errorabsenergy[WHICHERROR],itersenergyold,itersenergy);
        }
      }// end if using ramesh

      




      ///////////
      //
      // use best result from trying harder over all energy attempts
      //
      ///////////
      *lpflag=lpflagenergybest;
      *lpflagrad=lpflagradenergybest;
      radinvmodenergy=radinvmodenergybest;
      radErfnegenergy=radErfnegenergybest;
      failreturnenergy=failreturnenergybest;
      eomtypeenergy=eomtypeenergybest;
      errorabsenergy[0]=errorabsenergybest[0];
      errorabsenergy[1]=errorabsenergybest[1];
      PLOOP(pliter,pl){
        pbenergy[pl]=pbenergybest[pl];
        uubenergy[pl]=uubenergybest[pl];
        SCLOOP(sc) dUcompenergy[sc][pl]=dUcompenergybest[sc][pl];
      }
      qenergy=qenergybest;
      goexplicitenergy=goexplicitenergybest;

      // store these in case energy ultimately used
      lpflagenergy=*lpflag;
      lpflagradenergy=*lpflagrad;

      



     
    }// end if doing GRMHD inversion
    else{
      // if didn't do energy inversion, treat as failure of said inversion
      failreturnenergy=FAILRETURNGENERAL;
    }







    if(
       eomtypecond[EOMENTROPYGRMHD]==1
       || (pbenergy[UU]<=0.0 && ACCEPTASNOFAILURE(failreturnenergy)==1 || ACTUALHARDFAILURE(failreturnenergy)==1)
       ){
      //////////////////////////////////
      //
      // GET ENTROPY (currently, only if failure for energy solver)
      //
      //////////////////////////////////


      // quickly try QTYPMHD then QTYURAD
      int tryphase1;

      int baseitermethodlist[NUMPHASESENT]={QTYPMHD,QTYENTROPYUMHD,QTYURAD,QTYPRAD,QTYPMHD,QTYENTROPYUMHD,QTYURAD,QTYPRAD}; int whichfirstpmhd=0,whichfirstentropyumhd=1,whichfirsturad=2,whichfirstprad=3;
      int itermodelist[NUMPHASESENT]={ITERMODENORMAL,ITERMODENORMAL,ITERMODENORMAL,ITERMODENORMAL,ITERMODESTAGES,ITERMODESTAGES,ITERMODESTAGES,ITERMODESTAGES};
      FTYPE trueimptryconvlist[NUMPHASESENT]={IMPTRYCONV,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK};
      FTYPE trueimpokconvlist[NUMPHASESENT]={IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK,IMPTRYCONVQUICK};
      FTYPE trueimpallowconvconstlist[NUMPHASESENT]={IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST,IMPALLOWCONVCONST};
      int trueimpmaxiterlist[NUMPHASESENT]={IMPMAXITERQUICK,IMPMAXITERQUICK,IMPMAXITERQUICK,IMPMAXITERQUICK,IMPMAXITERMEDIUM,IMPMAXITERMEDIUM,IMPMAXITERMEDIUM,IMPMAXITERMEDIUM};
      int truenumdampattemptslist[NUMPHASESENT]={NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK,NUMDAMPATTEMPTSQUICK};
      int modprimlist[NUMPHASESENT]={0,0,0,0,1,0,0,0};
      int checkradinvlist[NUMPHASESENT]={0,1,1,1,0,1,1,1};
      int eomtypelist[NUMPHASESENT]={EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD,EOMENTROPYGRMHD};

      // results in list
      FTYPE errorabslist[NUMPHASES][NUMERRORTYPES];
      set_array(errorabslist,NUMPHASES*NUMERRORTYPES,MPI_FTYPE,1.0);
      int radinvmodlist[NUMPHASES];
      set_array(radinvmodlist,NUMPHASES,MPI_INT,UTOPRIMRADFAILBAD1);

      if(gasprimaryevolves==0){ // if gas primary evolves, then keep original order because PMHD is fastest method.
        if(radprimaryevolves){
          if(
             (USEPRIORITERMETHOD && (prioritermethodlist[BASEITERMETHODINDEX] == PRIORITERMETHODNOTSET ||  prioritermethodlist[BASEITERMETHODINDEX] == QTYURAD))
             || USEPRIORITERMETHOD==0
             ){
            tryphase1=-1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD;        modprimlist[tryphase1]=0; whichfirsturad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD;        modprimlist[tryphase1]=0; whichfirstprad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD;        modprimlist[tryphase1]=0; whichfirstpmhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYENTROPYUMHD; modprimlist[tryphase1]=0; whichfirstentropyumhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD;        modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD;        modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD;        modprimlist[tryphase1]=1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYENTROPYUMHD; modprimlist[tryphase1]=0;
          }
          if(
             (USEPRIORITERMETHOD && prioritermethodlist[BASEITERMETHODINDEX] == QTYPRAD)
             ){
            tryphase1=-1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD;        modprimlist[tryphase1]=0; whichfirstprad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD;        modprimlist[tryphase1]=0; whichfirsturad=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD;        modprimlist[tryphase1]=0; whichfirstpmhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYENTROPYUMHD; modprimlist[tryphase1]=0; whichfirstentropyumhd=tryphase1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPRAD;        modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYURAD;        modprimlist[tryphase1]=0;
            tryphase1++; baseitermethodlist[tryphase1]=QTYPMHD;        modprimlist[tryphase1]=1;
            tryphase1++; baseitermethodlist[tryphase1]=QTYENTROPYUMHD; modprimlist[tryphase1]=0;
          }
        }
      }


      int firsttryphase1used;
      for(tryphase1=0;tryphase1<NUMPHASESENT;tryphase1++){

        // pick best simple method avoids all solvers except PMHD
        if(MODEMETHOD==MODEPICKBESTSIMPLE && baseitermethodlist[tryphase1]!=QTYPMHD) continue;

        // pick best simple 2 method avoids all itermodestages methods
        if(MODEMETHOD==MODEPICKBESTSIMPLE2 && itermodelist[tryphase1]==ITERMODESTAGES) continue;

        // avoid method in case very non-dominant since then would give errorneous (critically bad even) results.  If methods that can use fail, have to revert to entropy or fixups.
        if(radextremeprimaryevolves && (baseitermethodlist[tryphase1]==QTYPMHD || baseitermethodlist[tryphase1]==QTYENTROPYUMHD) ) continue;
        if(gasextremeprimaryevolves && (baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD)) continue;

        if(radextremeprimaryevolves==0 && (baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD) ){
          // if not in extreme radiative regime, then don't do damping for raditive schemes that are slow.
          truenumdampattemptslist[tryphase1]=NUMDAMPATTEMPTSQUICK;
          // and don't try to get too good of error.
          trueimptryconvlist[tryphase1]=IMPTRYCONVQUICK;
          // and avoid stages
          itermodelist[tryphase1]=ITERMODENORMAL;
        }
        if(gasextremeprimaryevolves==0 && (baseitermethodlist[tryphase1]==QTYPMHD)){
          // then still try to get good error since fast method, so no changes.
        }


        if(MODEMETHOD==MODEPICKBESTSIMPLE){
          // do nothing, don't skip stages.
        }
        else{
          // If already tried QTYPMHD and that got error(0)<tol but error(1)>tol, then skip QTYPMHD with stages since should switch to QTYURAD or QTYPRAD because STAGES won't improve on error(1) if error(0)<tol.
          if(baseitermethodlist[tryphase1]==QTYPMHD && errorabslist[whichfirstpmhd][0]<IMPTRYCONV && errorabslist[whichfirstpmhd][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon radiative solvers
            continue;
          }

          if(baseitermethodlist[tryphase1]==QTYENTROPYUMHD && errorabslist[whichfirstentropyumhd][0]<IMPTRYCONV && errorabslist[whichfirstentropyumhd][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon others
            continue;
          }

          if(baseitermethodlist[tryphase1]==QTYURAD && errorabslist[whichfirsturad][0]<IMPTRYCONV && errorabslist[whichfirsturad][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon others
            continue;
          }

          if(baseitermethodlist[tryphase1]==QTYPRAD && errorabslist[whichfirstprad][0]<IMPTRYCONV && errorabslist[whichfirstprad][1]>IMPTRYCONV && WHICHERROR==1){
            // then should skip this case and rely upon others
            continue;
          }
        }


        // FUCK
        //if(baseitermethodlist[tryphase1]==QTYENTROPYUMHD || baseitermethodlist[tryphase1]==QTYURAD || baseitermethodlist[tryphase1]==QTYPRAD) continue; // skip this method for now.
        if(baseitermethodlist[tryphase1]==QTYENTROPYUMHD) continue; // skip this method for now.



        // consider radinvmod only if error bad for original approach.  Avoids excessive attempts when should hit radiative ceiling and error is small.
        if(RADINVBAD(radinvmodentropybest) && checkradinvlist[tryphase1] || ACTUALHARDORSOFTFAILURE(failreturnentropybest) && failreturnentropybest!=FAILRETURNMODESWITCH){
          if(firsttryphase1used==-1) firsttryphase1used=tryphase1;
          tryphaselistentropy[tryphase1]++;


          itersentropyold=itersentropy;
          errorabsentropyold[0]=errorabsentropy[0];
          errorabsentropyold[1]=errorabsentropy[1];
          radinvmodentropyold=radinvmodentropy;
          //
          whichcapentropy=CAPTYPEBASIC;
          baseitermethodentropy=baseitermethodlist[tryphase1];
          itermodeentropy=itermodelist[tryphase1];
          trueimptryconventropy=trueimptryconvlist[tryphase1];
          trueimpokconventropy=(MAX(trueimptryconvlist[tryphase1],trueimpokconvlist[tryphase1]));
          trueimpallowconventropy=(MAX(trueimptryconvlist[tryphase1],trueimpallowconvconstlist[tryphase1]));
          trueimpmaxiterentropy=trueimpmaxiterlist[tryphase1];
          truenumdampattemptsentropy=truenumdampattemptslist[tryphase1];
          int modprim=0;
          if(havebackup==0) modprim=modprimlist[tryphase1];
          //
          // start fresh
          *lpflag=UTOPRIMNOFAIL;
          *lpflagrad=UTOPRIMRADNOFAIL;
          radinvmodentropy=0;
          radErfnegentropy=0;
          failreturnentropy=FAILRETURNGENERAL;// default to fail
          eomtypeentropy=eomtypelist[tryphase1];
          errorabsentropy[0]=errorabsentropy[1]=1.0;
          //
          // setup guess
          if(ACTUALHARDFAILURE(failreturnenergy)==1 || errorabsenergy[WHICHERROR]>ERRORTOUSEENTROPYFORENERGYGUESS || pbenergy[UU]<=0.0){
            if(errorabsentropybest[WHICHERROR]>TRYHARDERFEEDGUESSTOL){
              PLOOP(pliter,pl){
                pbentropy[pl]=pbbackup[pl];
                uubentropy[pl]=uubbackup[pl];
                SCLOOP(sc) dUcompentropy[sc][pl]=dUcompbackup[sc][pl];
              }
              qentropy=qbackup;
              errorabsentropy[0]=errorabsentropy[1]=1.0;
              radinvmodentropy=UTOPRIMRADFAILBAD1;
            }
            else{
              PLOOP(pliter,pl){
                pbentropy[pl]=pbentropybest[pl];
                uubentropy[pl]=uubentropybest[pl];
                SCLOOP(sc) dUcompentropy[sc][pl]=dUcompbackup[sc][pl];
              }
              qentropy=qentropybest;
              errorabsentropy[0]=errorabsentropybest[0];
              errorabsentropy[1]=errorabsentropybest[1];
              radinvmodentropy=radinvmodentropybest;
            }
          }
          else{ // start with energy
            PLOOP(pliter,pl){
              pbentropy[pl]=pbenergy[pl];
              SCLOOP(sc) dUcompentropy[sc][pl]=dUcompbackup[sc][pl];
            }
            qentropy=qenergy;
            errorabsentropy[0]=errorabsenergy[0]; // just for setting guess, not accepted error for entropy
            errorabsentropy[1]=errorabsenergy[1]; // just for setting guess, not accepted error for entropy
            radinvmodentropy=radinvmodenergy;
          }
          //
          //
          FTYPE fracenergyentropy=0;
          failreturnentropy=koral_source_rad_implicit_mode(0,modprim,havebackup, didentropyalready, &eomtypeentropy, whichcapentropy, itermodeentropy, &baseitermethodentropy, trueimptryconventropy, trueimpokconventropy, trueimpallowconventropy, trueimpmaxiterentropy,  truenumdampattemptsentropy, fracenergyentropy, dissmeasure, &radinvmodentropy, pbentropy, uubentropy, piin, Uiin, Ufin, CUf, ptrgeom, &qentropy, dUother ,dUcompentropy, errorabsentropy, errorabsentropybest, &itersentropy, &f1itersentropy, &nummhdinvsentropy, &nummhdstepsentropy);

          // store error, no matter if using solution or not or even if explicit
          errorabslist[tryphase1][0]=errorabsentropy[0];
          errorabslist[tryphase1][1]=errorabsentropy[1];

          if(failreturnentropy==FAILRETURNGOEXPLICIT){
            lpflagentropybest=*lpflag;
            lpflagradentropybest=*lpflagrad;
            radinvmodentropybest=radinvmodentropy;
            radErfnegentropybest=radErfnegentropy;
            failreturnentropybest=failreturnentropy;
            eomtypeentropybest=eomtypelist[tryphase1];
            goexplicitentropybest=1;

            // save which method ended up using
            int prioriter; for(prioriter=0;prioriter<NUMPRIORITERMETHODINDEX;prioriter++){
              GLOBALMACP0A1(prioritermethod,i,j,k,prioriter) = PRIORITERMETHODNOTSET; // back to unknown
            }
          }
          else{
            // see if want to keep
            //            if((errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && (RADINVBAD(radinvmodentropy)==0 || RADINVBAD(radinvmodentropybest) && RADINVBAD(radinvmodentropy)) || RADINVBAD(radinvmodentropybest) && RADINVBAD(radinvmodentropy)==0 && (errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR]||errorabsentropy[WHICHERROR]<IMPOKCONVABS) || RADINVBAD(radinvmodentropybest)==0 && RADINVBAD(radinvmodentropy) && (errorabsentropybest[WHICHERROR]>IMPOKCONVABS && errorabsentropy[WHICHERROR]<IMPOKCONVABS)) && ACTUALHARDFAILURE(failreturnentropy)==0){
            if(
               (
                // best and no radinv problem
                errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && RADINVBAD(radinvmodentropy)==0
                // best because prior best also had radinv problem and still smaller error
                || errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && (RADINVBAD(radinvmodentropy) && RADINVBAD(radinvmodentropybest))
                // best because even though has radinv problem (while best doesn't) much smaller error
                || errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && errorabsentropy[WHICHERROR]<IMPOKCONVABS && errorabsentropybest[WHICHERROR]>100.0*IMPOKCONVABS && (RADINVBAD(radinvmodentropy) && RADINVBAD(radinvmodentropybest)==0)
                )
               && ACTUALHARDFAILURE(failreturnentropy)==0){

              // store result in case better than latter results
              lpflagentropybest=*lpflag;
              lpflagradentropybest=*lpflagrad;
              radinvmodentropybest=radinvmodentropy;
              radErfnegentropybest=radErfnegentropy;
              failreturnentropybest=failreturnentropy;
              eomtypeentropybest=eomtypeentropy;
              goexplicitentropybest=0;
              errorabsentropybest[0]=errorabsentropy[0];
              errorabsentropybest[1]=errorabsentropy[1];
              PLOOP(pliter,pl){
                pbentropybest[pl]=pbentropy[pl];
                uubentropybest[pl]=uubentropy[pl];
                SCLOOP(sc) dUcompentropybest[sc][pl]=dUcompentropy[sc][pl];
              }
              qentropybest=qentropy;

              // save which method ended up using
              GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) = baseitermethodlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) = itermodelist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,IMPTRYCONVINDEX) = (trueimptryconvlist[tryphase1]==IMPTRYCONV);
              GLOBALMACP0A1(prioritermethod,i,j,k,IMPMAXITERINDEX) = trueimpmaxiterlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,NUMDAMPINDEX) = truenumdampattemptslist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,MODPRIMINDEX) = modprimlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,CHECKRADINVINDEX) = checkradinvlist[tryphase1];
              GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = eomtypelist[tryphase1];

            }
          }

          // regardless, still cost more iters
          itersentropy+=itersentropyold;

          int noproblem=ACTUALHARDORSOFTFAILURE(failreturnentropy)==0;

          if(tryphase1!=firsttryphase1used && noproblem){
            if(debugfail>=2) dualfprintf(fail_file,"Recovered using tryphase1=%d (entropy: %d radinvmod=%d -> %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g (best=%21.15g) iters: %d->%d\n",tryphase1,failreturnentropy,radinvmodentropyold,radinvmodentropy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsentropyold[WHICHERROR],errorabsentropy[WHICHERROR],errorabsentropybest[WHICHERROR],itersentropyold,itersentropy);
          }
          if(noproblem==0){
            if(debugfail>=2) dualfprintf(fail_file,"Failed to: <Recovered> using tryphase1=%d (entropy: %d radinvmod=%d -> %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g (best=%21.15g) iters: %d->%d\n",tryphase1,failreturnentropy,radinvmodentropyold,radinvmodentropy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsentropyold[WHICHERROR],errorabsentropy[WHICHERROR],errorabsentropybest[WHICHERROR],itersentropyold,itersentropy);
          }

        }// done trying harder
      }// end over tryphase1 loop





      if(
         !(MODEMETHOD==MODEPICKBESTSIMPLE2) // Ramesh method is like stages, so avoid if simple2 method
         && ACTUALHARDORSOFTFAILURE(failreturnentropy) && USERAMESH && radextremeprimaryevolves==0
         && !(errorabslist[whichfirstpmhd][0]<IMPTRYCONV && errorabslist[whichfirstpmhd][1]>IMPTRYCONV)
         ){ // try ramesh
        errorabsentropyold[0]=errorabsentropy[0];
        errorabsentropyold[1]=errorabsentropy[1];
        radinvmodentropyold=radinvmodentropy;
        itersentropyold=itersentropy;
        goexplicitentropy=0; // doesn't check, so force implicit
        //
        //
        if(gotrameshsolution){// then just recall solution
          PLOOP(pliter,pl){
            pbenergy[pl]=ppeng[pl];
            uubenergy[pl]=uueng[pl];
            SCLOOP(sc) dUcompenergy[sc][pl]=dUcompeng[sc][pl];
          }
          qenergy=qeng;
          //
          *lpflag=*lpflagrad=(PFTYPE)failtypeeng; // need better translation
          radinvmodenergy=radinvmodeng;
          radErfnegenergy=0; // not allowed, considered BADNEG type
          itersenergy=iterseng;
          errorabsenergy[0]=errorabseng[0];
          errorabsenergy[1]=errorabseng[1];
          if(failtypeeng || errorabsenergy[WHICHERROR]>IMPALLOWCONVCONSTABS){ failreturnenergy=FAILRETURNGENERAL; eomtypeenergy=EOMGRMHD; }
          else if(errorabsenergy[WHICHERROR]<IMPTRYCONVABS){ failreturnenergy=FAILRETURNNOFAIL; eomtypeenergy=EOMDIDGRMHD;}
          else{ failreturnenergy=FAILRETURNNOTTOLERROR; eomtypeenergy=EOMDIDGRMHD;}
        }
        else{ // else get new ramesh solution
          FTYPE errorabsforramesh[NUMERRORTYPES];
          set_array(errorabsforramesh,NUMERRORTYPES,MPI_FTYPE,1.0);
          // start fresh with ramesh
          if(errorabsentropy[WHICHERROR]>TRYHARDERFEEDGUESSTOL){
            PLOOP(pliter,pl){
              pbentropy[pl]=pbbackup[pl];
              uubentropy[pl]=uubbackup[pl];
              SCLOOP(sc) dUcompent[sc][pl]=dUcompentropy[sc][pl]=dUcompbackup[sc][pl];
            }
            qentropy=qbackup;
            errorabsent[0]=errorabsent[1]=1.0;
            errorabsforramesh[0]=errorabsforramesh[1]=1.0;
            radinvmodent=UTOPRIMRADFAILBAD1;
          }
          else{
            PLOOP(pliter,pl){
              pbentropy[pl]=pbentropybest[pl];
              SCLOOP(sc) dUcompent[sc][pl]=dUcompentropy[sc][pl]=dUcompbackup[sc][pl]; // always backup
            }
            qentropy=qentropybest;
            errorabsforramesh[0]=errorabsentropybest[0];
            errorabsforramesh[1]=errorabsentropybest[1];
            errorabsent[0]=errorabsentropybest[0];
            errorabsent[1]=errorabsentropybest[1];
            radinvmodent=radinvmodentropybest;
          }
          //
          // BEGIN GET RAMESH SOLUTION
          failreturnentropy=FAILRETURNGENERAL;// default to fail
          eomtypeentropy=EOMENTROPYGRMHD;
          int whichcall=eomtypeentropy; // real entropy call
          get_rameshsolution_wrapper(whichcall, eomtypeentropy, errorabsforramesh, ptrgeom, pbentropy, piin, Uiin, Ufin, dUother, CUf, q, ppeng, ppent, uueng, uuent, dUcompeng, dUcompent, &qeng, &qent, &failtypeeng, errorabseng, &iterseng, &radinvmodeng, &failtypeent, errorabsent, &itersent, &radinvmodent);
          gotrameshsolution=1; // indicates did at least attempt ramesh solution call
          // translate
          PLOOP(pliter,pl){
            pbentropy[pl]=ppent[pl];
            uubentropy[pl]=uuent[pl];
            SCLOOP(sc) dUcompentropy[sc][pl]=dUcompent[sc][pl];
          }
          qentropy=qent;
          //
          *lpflag=*lpflagrad=(PFTYPE)failtypeent; // need better translation
          radinvmodentropy=radinvmodent;
          radErfnegentropy=0; // not allowed, considered BADNEG type
          itersentropy=itersent;
          errorabsentropy[0]=errorabsent[0];
          errorabsentropy[1]=errorabsent[1];
          if(failtypeent || errorabsentropy[WHICHERROR]>IMPALLOWCONVCONSTABS){ failreturnentropy=FAILRETURNGENERAL; eomtypeentropy=EOMENTROPYGRMHD; }
          else if(errorabsentropy[WHICHERROR]<IMPTRYCONVABS){ failreturnentropy=FAILRETURNNOFAIL; eomtypeentropy=EOMDIDENTROPYGRMHD;}
          else{ failreturnentropy=FAILRETURNNOTTOLERROR; eomtypeentropy=EOMDIDENTROPYGRMHD;}
          // END GET RAMESH SOLUTION
          //
        }// end else if need to get ramesh solution
        // see if want to keep
        //        if((errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && (RADINVBAD(radinvmodentropy)==0 || RADINVBAD(radinvmodentropybest) && RADINVBAD(radinvmodentropy)) || RADINVBAD(radinvmodentropybest) && RADINVBAD(radinvmodentropy)==0 && (errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR]||errorabsentropy[WHICHERROR]<IMPOKCONVABS) || RADINVBAD(radinvmodentropybest)==0 && RADINVBAD(radinvmodentropy) && (errorabsentropybest[WHICHERROR]>IMPOKCONVABS && errorabsentropy[WHICHERROR]<IMPOKCONVABS)) && ACTUALHARDFAILURE(failreturnentropy)==0){
        if(
           (
            // best and no radinv problem
            errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && RADINVBAD(radinvmodentropy)==0
            // best because prior best also had radinv problem and still smaller error
            || errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && (RADINVBAD(radinvmodentropy) && RADINVBAD(radinvmodentropybest))
            // best because even though has radinv problem (while best doesn't) much smaller error
            || errorabsentropy[WHICHERROR]<errorabsentropybest[WHICHERROR] && errorabsentropy[WHICHERROR]<IMPOKCONVABS && errorabsentropybest[WHICHERROR]>100.0*IMPOKCONVABS && (RADINVBAD(radinvmodentropy) && RADINVBAD(radinvmodentropybest)==0)
            )
           && ACTUALHARDFAILURE(failreturnentropy)==0){

          if(ACCEPTASNOFAILURE(failreturnentropy)) usedrameshentropy=1; // means will use this actually, not just best yet no good enough
          // store result in case better than latter results
          lpflagentropybest=*lpflag;
          lpflagradentropybest=*lpflagrad;
          radinvmodentropybest=radinvmodentropy;
          radErfnegentropybest=radErfnegentropy;
          failreturnentropybest=failreturnentropy;
          eomtypeentropybest=eomtypeentropy;
          errorabsentropybest[0]=errorabsentropy[0];
          errorabsentropybest[1]=errorabsentropy[1];
          PLOOP(pliter,pl){
            pbentropybest[pl]=pbentropy[pl];
            uubentropybest[pl]=uubentropy[pl];
            SCLOOP(sc) dUcompentropybest[sc][pl]=dUcompentropy[sc][pl];
          }
          qentropybest=qentropy;

          // save which method ended up using
          GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) = QTYPMHD;
          GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) = ITERMODESTAGES;
          GLOBALMACP0A1(prioritermethod,i,j,k,IMPMAXITERINDEX) = IMPMAXITERLONG;
          GLOBALMACP0A1(prioritermethod,i,j,k,NUMDAMPINDEX) = NUMDAMPATTEMPTSQUICK;
          GLOBALMACP0A1(prioritermethod,i,j,k,MODPRIMINDEX) = 0;
          GLOBALMACP0A1(prioritermethod,i,j,k,CHECKRADINVINDEX) = 1;
          GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMENTROPYGRMHD;

        }
        // still cost more iters
        itersentropy+=itersentropyold;
        if(ACTUALHARDORSOFTFAILURE(failreturnentropy)==0){
          if(debugfail>=2) dualfprintf(fail_file,"Recovered using ramesh(%d) (entropy: radinvmod=%d -> %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g iters: %d->%d\n",failtypeent,radinvmodentropyold,radinvmodentropy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsentropyold[WHICHERROR],errorabsentropy[WHICHERROR],itersentropyold,itersentropy);
        }
        else{
          if(debugfail>=2) dualfprintf(fail_file,"Failed to: Recovered using ramesh(%d) (entropy: radinvmod=%d -> %d): ijknstepsteppart=%d %d %d %ld %d : error: %21.15g->%21.15g iters: %d->%d\n",failtypeent,radinvmodentropyold,radinvmodentropy,ptrgeom->i,ptrgeom->j,ptrgeom->k,nstep,steppart,errorabsentropyold[WHICHERROR],errorabsentropy[WHICHERROR],itersentropyold,itersentropy);
        }
      }





      ////////////
      //
      // use best result
      //
      /////////////
      *lpflag=lpflagentropybest;
      *lpflagrad=lpflagradentropybest;
      radinvmodentropy=radinvmodentropybest;
      radErfnegentropy=radErfnegentropybest;
      failreturnentropy=failreturnentropybest;
      eomtypeentropy=eomtypeentropybest;
      errorabsentropy[0]=errorabsentropybest[0];
      errorabsentropy[1]=errorabsentropybest[1];
      PLOOP(pliter,pl){
        pbentropy[pl]=pbentropybest[pl];
        uubentropy[pl]=uubentropybest[pl];
        SCLOOP(sc) dUcompentropy[sc][pl]=dUcompentropybest[sc][pl];
      }
      qentropy=qentropybest;
      goexplicitentropy=goexplicitentropybest;


      // store these in case entropy ultimately used (if goexplicit, then these are set by last call to _mode() function)
      lpflagentropy=*lpflag;
      lpflagradentropy=*lpflagrad;


      // eomtypeentropy can become EOMDONOTHING if this call was successful
    }
    else{
      failreturnentropy=FAILRETURNGENERAL;// fail if wasn't done
    }



  







    ///////////////////
    //
    // See if want to try getting "cold" solution where only momentum exchange is considered
    //
    //////////////////////    


    /////////
    //
    // COLD VARIABLES
    //
    ////////

    // holds iterations
    int iterscold=0;
    int f1iterscold=0;
    int nummhdinvscold=0;
    int nummhdstepscold=0;

    // implicit cold solution held in these variables
    FTYPE pbcold[NPR]; PLOOP(pliter,pl) pbcold[pl]=pbbackup[pl];
    FTYPE uubcold[NPR];  PLOOP(pliter,pl) uubcold[pl]=uubbackup[pl]; // holds returned uu from implicit solver
    FTYPE dUcompcold[NUMSOURCES][NPR]; PLOOP(pliter,pl) SCLOOP(sc) dUcompcold[sc][pl]=dUcompbackup[sc][pl];
    struct of_state qcold=qbackup;
    PFTYPE lpflagcold=1;
    PFTYPE lpflagradcold=1;
    int radinvmodcold=UTOPRIMRADFAILBAD1;
    int radErfnegcold=1;
    int failreturncold=FAILRETURNGENERAL; // default to fail in case cold not to be done at all
    int eomtypecold=EOMCOLDGRMHD;
    FTYPE errorabscold[NUMERRORTYPES];
    set_array(errorabscold,NUMERRORTYPES,MPI_FTYPE,1.0);
    int goexplicitcold=0;


    // try cold if energy and entropy both failed or the successfully converged versions gave u_g<=0.0
    if( 
       (pbenergy[UU]<=0.0 && ACCEPTASNOFAILURE(failreturnenergy)==1 || ACCEPTASNOFAILURE(failreturnenergy)==0)
       &&
       (pbentropy[UU]<=0.0 && ACCEPTASNOFAILURE(failreturnentropy)==1 || ACCEPTASNOFAILURE(failreturnentropy)==0)
        ){

      //////
      //
      // measure whether flow really cold enough to even use cold inversion
      //
      /////
      int iscoldflow;

      int jj;
      FTYPE usq=0.0;  SLOOPA(jj) usq+=q->ucon[jj]*q->ucov[jj];
      FTYPE ucovgas[NDIM],ucovrad[NDIM];
      DLOOPA(jj) ucovgas[jj]=uu0[UU+jj];
      DLOOPA(jj) ucovrad[jj]=uu0[URAD0+jj];
      FTYPE ucongas[NDIM],uconrad[NDIM];
      raise_vec(ucovgas,ptrgeom,ucongas);
      raise_vec(ucovrad,ptrgeom,uconrad);
      FTYPE ugas=0.0;  SLOOPA(jj) ugas+=ucovgas[jj]*ucongas[jj];
      FTYPE urad=0.0;  SLOOPA(jj) urad+=ucovrad[jj]*uconrad[jj];
#define COLDFACTOR (0.1)
      iscoldflow = (pb[UU]<COLDFACTOR*pb[RHO]*fabs(usq) && fabs(ucongas[TT]*ucovgas[TT])<COLDFACTOR*fabs(ugas) && fabs(uconrad[TT]*ucovrad[TT])<COLDFACTOR*fabs(urad));

      //      iscoldflow=1;

      //      dualfprintf(fail_file,"iscoldflow: %d : ug=%21.15g Erf=%21.15g rhovsq=%21.15g : Esqgas=%21.15g usqgas=%21.15g Esqrad=%21.15g usqrad=%21.15g\n",iscoldflow,pb[UU],pb[URAD0],pb[RHO]*fabs(usq),fabs(ucongas[TT]*ucovgas[TT]),fabs(ugas),fabs(uconrad[TT]*ucovrad[TT]),fabs(urad));


      // if flow is cold, do inversion
      if(iscoldflow==1){

        // count attempt to use cold
        tryphaselistcold[0]++;

        failreturncold=FAILRETURNGENERAL; // if doing cold, default is fail.
        int whichcapcold=CAPTYPEBASIC;
        int itermodecold=ITERMODECOLD;
        int baseitermethodcold=QTYPMHD;
        FTYPE trueimptryconvcold=IMPTRYCONVQUICK;
        FTYPE trueimpokconvcold=IMPTRYCONVQUICK;
        FTYPE trueimpallowconvcold=IMPTRYCONVQUICK;
        int trueimpmaxitercold=IMPMAXITERQUICK;
        int truenumdampattemptscold=NUMDAMPATTEMPTSQUICK;
      

        // start fresh or use entropy as starting point
        *lpflag=UTOPRIMNOFAIL;
        *lpflagrad=UTOPRIMRADNOFAIL;
        radinvmodcold=0;
        radErfnegcold=0;
        failreturncold=FAILRETURNGENERAL; // default to fail in case cold not to be done at all
        eomtypecold=EOMCOLDGRMHD;
        goexplicitcold=0;
        errorabscold[0]=errorabscold[1]=1.0;
        FTYPE fracenergycold=0.0;
        FTYPE errorabscoldbest[NUMERRORTYPES];
        set_array(errorabscoldbest,NUMERRORTYPES,MPI_FTYPE,1.0);

        havebackup=0; // no backup since entropy and energy failed
        didentropyalready=0;
        PLOOP(pliter,pl){
          pbcold[pl]=pbbackup[pl];
          uubcold[pl]=uubbackup[pl];
          if(0){
            // but force cold setup
            pbcold[UU]=UUMINLIMIT;
            pbcold[URAD0]=ERADLIMIT;
          }
          else{
            // try just using static u_g and solving momentum equation as if not changing.
          }
          SCLOOP(sc) dUcompcold[sc][pl]=dUcompbackup[sc][pl];
        }
        qcold=qbackup;
        errorabscold[0]=errorabscold[1]=1.0;


        failreturncold=koral_source_rad_implicit_mode(0,0,havebackup, didentropyalready, &eomtypecold, whichcapcold, itermodecold, &baseitermethodcold, trueimptryconvcold, trueimpokconvcold, trueimpallowconvcold, trueimpmaxitercold, truenumdampattemptscold, fracenergycold, dissmeasure, &radinvmodcold, pbcold, uubcold, piin, Uiin, Ufin, CUf, ptrgeom, &qcold, dUother ,dUcompcold, errorabscold, errorabscoldbest, &iterscold, &f1iterscold, &nummhdinvscold, &nummhdstepscold);

        if(failreturncold==FAILRETURNGOEXPLICIT) goexplicitcold=1;
        // store these in case cold ultimately used
        lpflagcold=*lpflag;
        lpflagradcold=*lpflagrad;

        // save which method ended up using
        GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) = QTYPMHD;
        GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) = ITERMODECOLD;
        GLOBALMACP0A1(prioritermethod,i,j,k,IMPMAXITERINDEX) = IMPMAXITERQUICK;
        GLOBALMACP0A1(prioritermethod,i,j,k,NUMDAMPINDEX) = NUMDAMPATTEMPTS;
        GLOBALMACP0A1(prioritermethod,i,j,k,MODPRIMINDEX) = 0;
        GLOBALMACP0A1(prioritermethod,i,j,k,CHECKRADINVINDEX) = 1;
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMCOLDGRMHD;


      }// end if cold flow
    }// end if not failure












    /////////////
    //
    // see if should use the entropy solution
    //
    // Conditions to use entropy include:
    // 1) dissmeasure says should avoid energy
    // 2) energy failed
    // 3) energy predicts smaller u_g than entropy
    // 4) error in energy is large
    // 5) etc.
    // but keep changing conditions as learn more..
    //
    /////////////

    int doentropy=whetherdoentropy(ptrgeom, fracenergy, ACCEPTASNOFAILURE(failreturnentropy),ACCEPTASNOFAILURE(failreturnenergy), radinvmodentropy, radinvmodenergy, IMPTRYCONVABS, IMPOKCONVABS, IMPBADENERGY, errorabsentropy[WHICHERROR], errorabsenergy[WHICHERROR], pbentropy, pbenergy);
    if(goexplicitentropy==1) doentropy=0; // override


    // check if doing implicit (don't care about cold check)
    if(goexplicitentropy==0 && goexplicitenergy==0){
      usedimplicit=1;

      /////////////
      //
      // cases where must use entropy
      //
      /////////////
      if(doentropy){
        //      dualfprintf(fail_file,"USING ENTROPY\n");
        // tell an externals to switch to entropy
        // store these in case energy ultimately used
        *lpflag=lpflagentropy;
        *lpflagrad=lpflagradentropy;
        *eomtype=eomtypeentropy; // can be EOMDONOTHING if successful and small enough error
        // set result as entropy result
        PLOOP(pliter,pl){
          pb[pl]=pbentropy[pl];
          uub[pl]=uubentropy[pl];
          SCLOOP(sc) dUcomp[sc][pl]=dUcompentropy[sc][pl];
          *q=qentropy;
        }
        usedentropy=1;
        noprims=0;
        errorabs[0]=errorabsentropy[0];
        errorabs[1]=errorabsentropy[1];
        iters=itersentropy+itersenergy; // count both since did both
        f1iters=f1itersentropy+f1itersenergy; // count both since did both
        failreturn=failreturnentropy;
        failfinalreturn=0;

        // save
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMENTROPYGRMHD;

      }
      else if(ACCEPTASNOFAILURE(failreturnenergy)==1){ // automatically also done when fracenergy==1.0 (now, or when also fracenergy>0.0)
        //      dualfprintf(fail_file,"USING ENERGY: errorabsenergy=%g errorabsentropy=%g\n",errorabsenergy[WHICHERROR],errorabsentropy[WHICHERROR]);
        *lpflag=lpflagenergy;
        *lpflagrad=lpflagradenergy;
        *eomtype=eomtypeenergy; // can be EOMDONOTHING if successful
        // set result as energy result
        PLOOP(pliter,pl){
          pb[pl]=pbenergy[pl];
          uub[pl]=uubenergy[pl];
          SCLOOP(sc) dUcomp[sc][pl]=dUcompenergy[sc][pl];
          *q=qenergy;
        }
        usedenergy=1;
        noprims=0;
        errorabs[0]=errorabsenergy[0];
        errorabs[1]=errorabsenergy[1];
        iters=itersentropy+itersenergy; // count both since did both
        f1iters=f1itersentropy+f1itersenergy; // count both since did both
        failreturn=failreturnenergy;
        failfinalreturn=0;

        // save
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMGRMHD;

      }
      else if(ACCEPTASNOFAILURE(failreturncold)==1){
        //      dualfprintf(fail_file,"USING COLD: errorabscold=%g errorabsentropy=%g\n",errorabscold[WHICHERROR],errorabsentropy[WHICHERROR]);
        *lpflag=lpflagcold;
        *lpflagrad=lpflagradcold;
        *eomtype=eomtypecold; // can be EOMDONOTHING if successful
        // set result as cold result
        PLOOP(pliter,pl){
          pb[pl]=pbcold[pl];
          uub[pl]=uubcold[pl];
          SCLOOP(sc) dUcomp[sc][pl]=dUcompcold[sc][pl];
          *q=qcold;
        }
        usedcold=1;
        noprims=0;
        errorabs[0]=errorabscold[0];
        errorabs[1]=errorabscold[1];
        iters=itersentropy+itersenergy+iterscold; // count all since did all
        f1iters=f1itersentropy+f1itersenergy+f1iterscold; // ""
        failreturn=failreturncold;
        failfinalreturn=0;

        // need to fail on u_g and Erf so averages
        *lpflag=UTOPRIMFAILU2AVG1FROMCOLD;
        *lpflagrad=UTOPRIMFAILU2AVG1FROMCOLD;
        // only keep u_g and Erf static as backup
        pbcold[UU]=pbbackup[UU];
        pbcold[URAD0]=pbbackup[URAD0]; 

        // save
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMCOLDGRMHD;

      }
      else{ // force even if was default
        usedimplicit=0;

        // save
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) = EOMDEFAULT;

      }
    }



    if(usedimplicit==0){// No source either because goexplicitenergy==1 && goexplicitentropy==1 or failed
      if(goexplicitenergy==1 || goexplicitentropy==1 || *lpflagrad>UTOPRIMRADNOFAIL && *lpflag==UTOPRIMNOFAIL){ // then assume just radinvmod (i.e. any radiation failure is always fixable), so revert to explicit if couldn't find solution.
        *lpflag=UTOPRIMNOFAIL;
        *lpflagrad=UTOPRIMRADNOFAIL;
        noprims=1;
        failfinalreturn=1;
        *eomtype=EOMDEFAULT;
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) == *eomtype;

        if(debugfail>=3) dualfprintf(fail_file,"Went explicit: eenergy=%g eentropy=%g ienergy=%d ientropy=%d\n",errorabsenergy[WHICHERROR],errorabsentropy[WHICHERROR],itersenergy,itersentropy);
        if(goexplicitenergy==1 || goexplicitentropy==1){ usedexplicitgood=1; failfinalreturn=-1;}
        else{ usedexplicitkindabad=1; failfinalreturn=1;} // FUCK: might want to treat as actual failure if QTYPMHD mode since lpflag never set.
      }
      else{ // actual failure
        // if no source, then will do normal inversion (no change to *eomtype) as if G=0.
        GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) == *eomtype;
        *lpflag=UTOPRIMFAILCONV;
        *lpflagrad=UTOPRIMRADFAILCASE1A;
        if(errorabsenergy[WHICHERROR]<errorabsentropy[WHICHERROR]){
          // indicates error that could have had if chose to raise IMPALLOWCONV
          errorabs[0]=errorabsenergy[0];
          errorabs[1]=errorabsenergy[1];
        }
        else{
          errorabs[0]=errorabsentropy[0];
          errorabs[1]=errorabsentropy[1];
        }
        if(debugfail>=2) dualfprintf(fail_file,"No source: eenergy=%g eentropy=%g ienergy=%d ientropy=%d\n",errorabsenergy[WHICHERROR],errorabsentropy[WHICHERROR],itersenergy,itersentropy);
        failfinalreturn=1;
        usedexplicitbad=1;
      }
      // set prims and dU, but shouldn't be used
      PLOOP(pliter,pl){
        pb[pl]=pbbackup[pl];
        uub[pl]=uubbackup[pl];
        SCLOOP(sc) dUcomp[sc][pl]=dUcompbackup[sc][pl];
        *q=qbackup;
      }
      noprims=1;
      iters=itersentropy+itersenergy; // count both since did both
      f1iters=f1itersentropy+f1itersenergy; // count both since did both
      //      failreturn=0;
      // KORALTODO: But might want to fail more aggressively and report total failure.  Need to have estimate of whether G was important.
      failreturn=failreturnenergy;
    } // end if using implicit solution


    if(*eomtype>=0 && failreturnenergy>=0 && failreturnentropy>=0 && failreturncold>=0) dualfprintf(fail_file,"WTF: %d %d %d : %d : %d : %d\n",failreturnentropy,failreturnenergy,failreturncold,failfinalreturn,noprims,*eomtype);

  }// end MODEPICKBEST || MODEPICKBESTSIMPLE || MODEPICKBESTSIMPLE2



  // DEBUG:
  //  int prioriter; for(prioriter=0;prioriter<NUMPRIORITERMETHODINDEX;prioriter++){
  //    dualfprintf(fail_file,"INDEX: prioriter=%d : %d\n",prioriter,GLOBALMACP0A1(prioritermethod,i,j,k,prioriter));
  //  }






  // whether set some primitives (implies also failfinalreturn=0)
  if(noprims==0){
    if((pb[RHO]<=0.)&&(pb[UU]>=0.)) *lpflag= UTOPRIMFAILRHONEG;
    if((pb[RHO]>0.)&&(pb[UU]<0.))   *lpflag= UTOPRIMFAILUNEG;
    if((pb[RHO]<=0.)&&(pb[UU]<0.))  *lpflag= UTOPRIMFAILRHOUNEG;
    if(pb[PRAD0]<=0.) *lpflagrad = UTOPRIMRADFAILERFNEG;
  }


  // force field to evolve as directly.
  PLOOPBONLY(pl){
    pf[pl]=pb[pl]=uu0[pl];
    dUcomp[RADSOURCE][pl]=0.0; // always has to be.
  }
  


  // DEBUG:
  //  failfinalreturn=1;
  //  *eomtype=EOMGRMHD;




  // check if uncaught nan/inf
  int caughtnan=0;
  PLOOP(pliter,pl){
    if(!isfinite(pb[pl])) caughtnan++;
    if(!isfinite(pf[pl])) caughtnan++;
    if(!isfinite(uub[pl])) caughtnan++;
    if(!isfinite(dUcomp[RADSOURCE][pl])) caughtnan++;
  }
  if(caughtnan){
    // Doesn't seem to happen, even on Kraken
    if(debugfail>=2){
      dualfprintf(fail_file,"implicit solver generated nan result and it wasn't caught\n");
      PLOOP(pliter,pl) dualfprintf(fail_file,"1implicit solver: pl=%d pb=%21.15g pf=%21.15g dU=%21.15g\n",pl,pb[pl],pf[pl],dUcomp[RADSOURCE][pl]);
      int jj;
      DLOOPA(jj) dualfprintf(fail_file,"2implicit solver: jj=%d ucon=%21.15g ucov=%21.15g uradcon=%21.15g uradcov=%21.15g\n",jj,q->ucon[jj],q->ucov[jj],q->uradcon[jj],q->uradcov[jj]);
    }

    // choose as bad solution

    // if no source, then will do normal inversion (no change to *eomtype) as if G=0.
    GLOBALMACP0A1(prioritermethod,i,j,k,EOMTYPEINDEX) == *eomtype;
    *lpflag=UTOPRIMFAILCONV;
    *lpflagrad=UTOPRIMRADFAILCASE1A;
    if(debugfail>=2) dualfprintf(fail_file,"No sourceNAN\n");
    failfinalreturn=1;
    usedexplicitbad=1;
  
    // set prims and dU, but shouldn't be used
    PLOOP(pliter,pl){
      pb[pl]=pbbackup[pl];
      uub[pl]=uubbackup[pl];
      SCLOOP(sc) dUcomp[sc][pl]=dUcompbackup[sc][pl];
      *q=qbackup;
    }
    noprims=1;
    //      failreturn=0;
  }
  



  
  if(PRODUCTION==0&&errorabs[WHICHERROR]>IMPALLOWCONVCONSTABS && (usedenergy||usedentropy||usedcold)){
    dualfprintf(fail_file,"WTF2: %g : %d %d %d : %d %d\n",errorabs[WHICHERROR],usedenergy,usedentropy,usedcold,usedrameshenergy,usedrameshentropy);
    myexit(666);
  }





  if(PRODUCTION==0 && debugfail>=2 || PRODUCTION<=1){
    // this debug MPI stuff is very expensive


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

    static long long int numqtypmhd=0;
    static long long int numqtyumhd=0;
    static long long int numqtyprad=0;
    static long long int numqtyurad=0;
    static long long int numqtyentropyumhd=0;
    static long long int numqtyentropypmhd=0;
    static long long int numitermodenormal=0;
    static long long int numitermodestages=0;
    static long long int numitermodecold=0;

    static long long int numimplicits=0;
    static long long int numexplicitsgood=0;
    static long long int numexplicitskindabad=0;
    static long long int numexplicitsbad=0;
    static long long int numoff1iter=0,numofiter=0;
    static long long int numhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
    static long long int numhistiter[IMPMAXITERLONG+1]={0}; // histogram of error for implicit solver to be reported infrequently
    //    static long long int numindex[NUMPRIORITERMETHODINDEX]={0};

    // static counter for diagnosing issues
    static long long int totalnumenergy=0;
    static long long int totalnumentropy=0;
    static long long int totalnumcold=0;
    static long long int totalnumboth=0;
    static long long int totalnumbad=0;
    static long long int totalnumramesh=0;
    static long long int totalnumrameshenergy=0;
    static long long int totalnumrameshentropy=0;

    static long long int totalnumqtypmhd=0;
    static long long int totalnumqtyumhd=0;
    static long long int totalnumqtyprad=0;
    static long long int totalnumqtyurad=0;
    static long long int totalnumqtyentropyumhd=0;
    static long long int totalnumqtyentropypmhd=0;
    static long long int totalnumitermodenormal=0;
    static long long int totalnumitermodestages=0;
    static long long int totalnumitermodecold=0;


    static long long int totalnumimplicits=0;
    static long long int totalnumexplicitsgood=0;
    static long long int totalnumexplicitskindabad=0;
    static long long int totalnumexplicitsbad=0;
    static long long int totalnumoff1iter=0,totalnumofiter=0;
    static long long int totalnumhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
    static long long int totalnumhistiter[IMPMAXITERLONG+1]={0}; // histogram of error for implicit solver to be reported infrequently
    //    static long long int totalnumindex[NUMPRIORITERMETHODINDEX]={0};


    static long long int totaltryphaselistenergy[NUMPHASES];
    static long long int totaltryphaselistentropy[NUMPHASESENT];
    static long long int totaltryphaselistcold[NUMPHASESCOLD];


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
    numramesh+=gotrameshsolution;
    numrameshenergy+=usedrameshenergy;
    numrameshentropy+=usedrameshentropy;

    numqtypmhd += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYPMHD);
    numqtyumhd += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYUMHD);
    numqtyprad += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYPRAD);
    numqtyurad += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYURAD);
    numqtyentropyumhd += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYENTROPYUMHD);
    numqtyentropypmhd += (GLOBALMACP0A1(prioritermethod,i,j,k,BASEITERMETHODINDEX) == QTYENTROPYPMHD);
    numitermodenormal += (GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) == ITERMODENORMAL);
    numitermodestages += (GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) == ITERMODESTAGES);
    numitermodecold += (GLOBALMACP0A1(prioritermethod,i,j,k,ITERMODEINDEX) == ITERMODECOLD);




    // i=j=k=0 just to show infrequently
    if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0)){
      dualfprintf(fail_file,"nstep=%ld numimplicits=%lld numexplicitsgood=%lld numexplicitskindabad=%lld numexplicitsbad=%lld : numenergy=%lld numentropy=%lld numboth=%lld numcold=%lld : numbad=%lld : numramesh=%lld numrameshenergy=%lld numrameshentropy=%lld : averagef1iter=%g averageiter=%g  : numqtypmhd=%lld numqtyumhd=%lld numqtyprad=%lld numqtyurad=%lld numqtyentropyumhd=%lld numqtyentropypmhd=%lld numitermodenormal=%lld numitermodestages=%lld numitermodecold=%lld\n",nstep,numimplicits,numexplicitsgood,numexplicitskindabad,numexplicitsbad,numenergy,numentropy,numboth,numcold,numbad,numramesh,numrameshenergy,numrameshentropy,(FTYPE)numoff1iter/(SMALL+(FTYPE)numimplicits),(FTYPE)numofiter/(SMALL+(FTYPE)numimplicits),numqtypmhd,numqtyumhd,numqtyprad,numqtyurad,numqtyentropyumhd,numqtyentropypmhd,numitermodenormal,numitermodestages,numitermodecold);
      // counters for which method was *attempted* even if not used
      int oo;
      dualfprintf(fail_file,"tryenergy: ");
      for(oo=0;oo<NUMPHASES;oo++) dualfprintf(fail_file,"%lld ",tryphaselistenergy[oo]);
      dualfprintf(fail_file,"\n");
      dualfprintf(fail_file,"tryentropy: ");
      for(oo=0;oo<NUMPHASESENT;oo++) dualfprintf(fail_file,"%lld ",tryphaselistentropy[oo]);
      dualfprintf(fail_file,"\n");
      dualfprintf(fail_file,"trycold: ");
      for(oo=0;oo<NUMPHASESCOLD;oo++) dualfprintf(fail_file,"%lld ",tryphaselistcold[oo]);
      dualfprintf(fail_file,"\n");
    }
    
    numhisterr[MAX(MIN((int)(-log10l(SMALL+errorabs[WHICHERROR])),NUMNUMHIST-1),0)]++;
    numhistiter[MAX(MIN(iters,IMPMAXITERLONG),0)]++;
#define HISTREPORTSTEP (20)
    if(debugfail>=2 && nstep%HISTREPORTSTEP==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0){
      int histi;
      for(histi=0;histi<NUMNUMHIST;histi++){
        dualfprintf(fail_file,"numhisterr%d=%lld\n",histi,numhisterr[histi]);
      }
      for(histi=0;histi<=IMPMAXITERLONG;histi++){
        dualfprintf(fail_file,"numhistiter%d=%lld\n",histi,numhistiter[histi]);
      }
    }

    // only need to get totals over cores when need to show result
    // i=j=k=0 just to show infrequently
#define GETCOLLECTIVETOTALS ( (debugfail>=2 || PRODUCTION<=1 && nstep%HISTREPORTSTEP==0 && steppart==0 ) && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0) )

    if(GETCOLLECTIVETOTALS){
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

        MPI_Reduce(&numqtypmhd, &totalnumqtypmhd, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numqtyumhd, &totalnumqtyumhd, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numqtyprad, &totalnumqtyprad, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numqtyurad, &totalnumqtyurad, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numqtyentropyumhd, &totalnumqtyentropyumhd, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numqtyentropypmhd, &totalnumqtyentropypmhd, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numitermodenormal, &totalnumitermodenormal, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(&numitermodestages, &totalnumitermodestages, 1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);

        MPI_Reduce(numhisterr, totalnumhisterr, NUMNUMHIST, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(numhistiter, totalnumhistiter, IMPMAXITERLONG+1, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);

        // attempt counters
        MPI_Reduce(tryphaselistenergy, totaltryphaselistenergy, NUMPHASES, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(tryphaselistentropy, totaltryphaselistentropy, NUMPHASESENT, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
        MPI_Reduce(tryphaselistcold, totaltryphaselistcold, NUMPHASESCOLD, MPI_LONG_LONG_INT, MPI_SUM, MPIid[0], MPI_COMM_GRMHD);
      }
      if(myid==MPIid[0]){// only show result on one core that got the final result     
        dualfprintf(fail_file,"nstep=%ld totalnumimplicits=%lld totalnumexplicitsgood=%lld totalnumexplicitskindabad=%lld totalnumexplicitsbad=%lld : totalnumenergy=%lld totalnumentropy=%lld totalnumboth=%lld totalnumcold=%lld : totalnumbad=%lld : totalnumramesh=%lld totalnumrameshenergy=%lld totalnumrameshentropy=%lld : totalaveragef1iter=%g totalaverageiter=%g  : totalnumqtypmhd=%lld totalnumqtyumhd=%lld totalnumqtyprad=%lld totalnumqtyurad=%lld totalnumqtyentropyumhd=%lld totalnumqtyentropypmhd=%lld totalnumitermodenormal=%lld totalnumitermodestages=%lld totalnumitermodecold=%lld\n",nstep,totalnumimplicits,totalnumexplicitsgood,totalnumexplicitskindabad,totalnumexplicitsbad,totalnumenergy,totalnumentropy,totalnumboth,totalnumcold,totalnumbad,totalnumramesh,totalnumrameshenergy,totalnumrameshentropy,(FTYPE)totalnumoff1iter/(SMALL+(FTYPE)totalnumimplicits),(FTYPE)totalnumofiter/(SMALL+(FTYPE)totalnumimplicits),totalnumqtypmhd,totalnumqtyumhd,totalnumqtyprad,totalnumqtyurad,totalnumqtyentropyumhd,totalnumqtyentropypmhd,totalnumitermodenormal,totalnumitermodestages,totalnumitermodecold);
        // counters for which method was *attempted* even if not used
        int oo;
        dualfprintf(fail_file,"totaltryenergy: ");
        for(oo=0;oo<NUMPHASES;oo++) dualfprintf(fail_file,"%lld ",totaltryphaselistenergy[oo]);
        dualfprintf(fail_file,"\n");
        dualfprintf(fail_file,"totaltryentropy: ");
        for(oo=0;oo<NUMPHASESENT;oo++) dualfprintf(fail_file,"%lld ",totaltryphaselistentropy[oo]);
        dualfprintf(fail_file,"\n");
        dualfprintf(fail_file,"totaltrycold: ");
        for(oo=0;oo<NUMPHASESCOLD;oo++) dualfprintf(fail_file,"%lld ",totaltryphaselistcold[oo]);
        dualfprintf(fail_file,"\n");
      
        if(nstep%HISTREPORTSTEP==0){
          int histi;
          for(histi=0;histi<NUMNUMHIST;histi++){
            trifprintf("totalnumhisterr%d=%lld\n",histi,totalnumhisterr[histi]);
          }
          for(histi=0;histi<=IMPMAXITERLONG;histi++){
            trifprintf("totalnumhistiter%d=%lld\n",histi,totalnumhistiter[histi]);
          }
        }
      }// end if myid==0
    }// end if getting totals
  }// end if ok production level to get these diags



  return(failfinalreturn);
}







// compute changes to U (both T and R) using implicit method
// KORALTODO: If doing implicit, should also add geometry source term that can sometimes be stiff.  Would require inverting sparse 8x8 matrix (or maybe 6x6 since only r-\theta for SPC).  Could be important for very dynamic radiative flows.
static int koral_source_rad_implicit_mode(int allowbaseitermethodswitch, int modprim, int havebackup, int didentropyalready, int *eomtype, int whichcap, int itermode, int *baseitermethod, FTYPE trueimptryconv, FTYPE trueimpokconv, FTYPE trueimpallowconv, int trueimpmaxiter, int truenumdampattempts, FTYPE fracenergy, FTYPE dissmeasure, int *radinvmod, FTYPE *pb, FTYPE *uub, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR], FTYPE *errorabsreturn, FTYPE *errorabsbestexternal, int *itersreturn, int *f1itersreturn, int *nummhdinvsreturn, int *nummhdstepsreturn)
{
  int nstrokeorig=nstroke;
  *nummhdinvsreturn=0;
  *nummhdstepsreturn=0; // global number of strokes for this core at this point, add negative so can add position later and get difference additional steps for this _mode() call

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

  // store pb as might have didentropyalready=1 and then can use pborig[UU] as entropy's solution for u_g, and that can be used to create condition to avoid over-iterating with energy solver.
  FTYPE pborig[NPR];
  PLOOP(pliter,pl) pborig[pl]=pb[pl];
  int radinvmodorig=*radinvmod;

#if(DEBUGMAXITER)
  FTYPE pppreholdlist[IMPMAXITERLONG+2][NPR]={{0}}; // for debug
  FTYPE ppposholdlist[IMPMAXITERLONG+2][NPR]={{0}}; // for debug
  FTYPE f1reportlist[IMPMAXITERLONG+2][NPR]={{0}}; // for debug
  FTYPE f1list[IMPMAXITERLONG+2][NPR]={{0}}; // for debug
  FTYPE errorabsf1list[IMPMAXITERLONG+2]={0}; // for debug
  FTYPE errorallabsf1list[IMPMAXITERLONG+2]={0}; // for debug
  int realiterlist[IMPMAXITERLONG+2]; // for debug
  set_array(realiterlist,IMPMAXITERLONG+2,MPI_INT,-1);
  FTYPE jaclist[IMPMAXITERLONG+2][NDIM][NDIM]; // for debug
  set_array(jaclist,(IMPMAXITERLONG+2)*NDIM*NDIM,MPI_FTYPE,BIG);
  int implicititerlist[IMPMAXITERLONG+2]={0}; // for debug
  int implicitferrlist[IMPMAXITERLONG+2]={0}; // for debug
#define NUMFRACDAMP 10
  FTYPE fracdamplist[NUMFRACDAMP]={{0}};

#else
  FTYPE (*pppreholdlist)[NPR];
  FTYPE (*ppposholdlist)[NPR];
  FTYPE (*f1reportlist)[NPR];
  FTYPE (*f1list)[NPR];
  FTYPE *errorabsf1list;
  FTYPE *errorallabsf1list;
  int *realiterlist;
  FTYPE (*jaclist)[NDIM][NDIM];
  FTYPE *fracdamplist;
  int *implicititerlist;
  int *implicitferrlist;
#endif

  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uuppp[NPR],uu[NPR],uuporig[NPR],uu0orig[NPR],bestuu[NPR];
  FTYPE pp0[NPR],ppp[NPR],pppp[NPR],ppppp[NPR],pp[NPR],ppporig[NPR],pp0orig[NPR],bestpp[NPR];
  FTYPE f1[NPR],f1norm[NPR],f1report[NPR],f3report[NPR],lowestf1[NPR],lowestf1norm[NPR],lowestf1report[NPR],lowestf3report[NPR];
  FTYPE f1p[NPR];
  FTYPE pppriorsteptype[NPR],uupriorsteptype[NPR];

  FTYPE radsource[NPR], deltas[NPR]; 
  extern int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, int baseitermethod, FTYPE *errorabs, FTYPE *errorabsbestexternal, int iters, int iterstotal, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);
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
    lowestf1[pl]=BIG;
    lowestf1norm[pl]=BIG;
    lowestf1report[pl]=BIG;
    lowestf3report[pl]=BIG;
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
  int fakeitermethod=IMPMAXITERLONG;// used by normal stepping, so just make maximum
  define_method(fakeitermethod, &eomtypelocal, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
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
#define GETDEFAULTFAILURESTATE (IMPPMHDTYPE(implicititer)==0) // with new rad inversion scheme, seems to push through without issue to just let hit gamma ceiling temporarily -- even with RADPULSEPLANAR
  // Otherwise, takes *many* iterations, not just f1iters, to get solution.  Even though that was old code state that was working, not required now.


  // default is to allow no failure unless iterating MHD primitives in which case a radiation failure is ok.
  int failreturnallowable;

  // no need to worry about RAD failure if not iterating rad quantities
  // As stated above, with new rad inversion, seems ok to temporarily hit ceiling.
  // With tests like RADTUBE, can't allow radiative inversion cieling else dies.  While with tests like RADPULSE, fastest to converge to good solution with allowing hitting the ceiling.  So mixed issue.
  // With RADTUBE and QTYPMHD, must also avoid rad failure if possible to avoid lack of solution, hence &&0 below.
  if(TREATRADINVCAPASNONFAILUREFORPMHDMETHOD && IMPMHDTYPEBASE(*baseitermethod)==1){
    if(IMPMHDTYPE(implicititer)) failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
    else failreturnallowable=UTOPRIMNOFAIL;
  }
  else{
    if(IMPMHDTYPE(implicititer)&&0) failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
    else failreturnallowable=UTOPRIMNOFAIL;
  }

  int failreturnallowableuse=failreturnallowable;
  int failreturnallowablefirst=-1;

  if(GETDEFAULTFAILURESTATE && IMPPMHDTYPE(implicititer)==0){    // No need for failure state if doing MHD iteration using primitives
    // KORALTODO: Instead of doing this, just keep inversion error from previously as stored in pflag.  So don't reset pflag elsewhere in fixup_utoprim() for elsewhere
    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)

    ////////
    // be quick about checking this
    struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
    newtonstats.nstroke=newtonstats.lntries=0;
    // set inputs for errors, maxiters, etc.
#define IMPOKCONVCONSTFORDEFAULT (1E-6)
#define IMPALLOWCONVCONSTFORDEFAULT (1E-4)
#define IMPMAXITERFORDEFAULT (10)
    newtonstats.tryconv=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
    newtonstats.tryconvultrarel=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
    newtonstats.extra_newt_iter=1;
    newtonstats.extra_newt_iter_ultrarel=2;
    newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION;//IMPALLOWCONVCONSTFORDEFAULT;
    newtonstats.maxiter=MIN(trueimpmaxiter,IMPMAXITERFORDEFAULT);
    //
    int finalstep = 1;
    int doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    int checkoninversiongas;
    int checkoninversionrad;
    // don't check since slows down code and could be good enough solution if original error says ok.
    checkoninversiongas=checkoninversionrad=0;
    //
    failreturnallowable=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Uiin, q, ptrgeom, dissmeasure, prtestUiin, &newtonstats);
    *nummhdinvsreturn++;
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
  // Get p(uu0) for radiative guess for pb that is used in optically thin regime when USEDUINRADUPDATE==1
  //
  // get better version of pb that agrees with uu0.  This is only good for optically thin region, so best used with USEDUINRADUPDATE==1 for general optical depths.
  //
  // Getting prad(uu0_{rad parts}) is only relevant for PRAD method.  All other methods overwrite the radiative pb and pp.
  //
  //////////////////////
  if(GETRADINVFROMUU0FORPB==1 || GETRADINVFROMUU0FORPB==2 && IMPRADTYPEBASE(*baseitermethod) || GETRADINVFROMUU0FORPB==3 && (IMPRADTYPEBASE(*baseitermethod) || IMPMHDTYPEBASE(*baseitermethod) && itermode==ITERMODESTAGES)){

      struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
    if(0){
      // initialize counters
      newtonstats.nstroke=newtonstats.lntries=0;
      // set inputs for errors, maxiters, etc.
      newtonstats.tryconv=trueimptryconv;
      newtonstats.tryconvultrarel=trueimptryconv*1E-1; // just bit smaller, not as extreme as default
      newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION; //IMPALLOWCONV;
      newtonstats.maxiter=trueimpmaxiter;
      newtonstats.extra_newt_iter=0;
      newtonstats.extra_newt_iter_ultrarel=1;
    }
    else{
      ////////
      // be quick about checking this
      newtonstats.nstroke=newtonstats.lntries=0;
      // set inputs for errors, maxiters, etc.
      newtonstats.tryconv=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
      newtonstats.tryconvultrarel=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
      newtonstats.extra_newt_iter=1;
      newtonstats.extra_newt_iter_ultrarel=2;
      newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION; //IMPALLOWCONVCONSTFORDEFAULT;
      newtonstats.maxiter=MIN(trueimpmaxiter,IMPMAXITERFORDEFAULT);
      //
    }
    //
    int finalstep = 1;
    int doradonly=0;
    if(*baseitermethod==QTYPRAD){
      doradonly=1;
    }
    else doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    FTYPE prtest[NPR];
    int whichcapnew=CAPTYPEFIX2; // so Erf doesn't drop out for guess.
    int checkoninversiongas;
    int checkoninversionrad;
    // don't check since slows down code and could be good enough solution if original error says ok.
    checkoninversiongas=checkoninversionrad=0;
    //
    PLOOP(pliter,pl) prtest[pl]=pb[pl];
    failreturnallowable=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcapnew, EVOLVEUTOPRIM, UNOTHING, uu0, q, ptrgeom, dissmeasure, prtest, &newtonstats);
    if(doradonly==0) *nummhdinvsreturn++;
    if(*baseitermethod==QTYPRAD && (1||failreturnallowable==UTOPRIMGENWRAPPERRETURNNOFAIL)){ // 1|| so assigns if failed or not, because with CAPTYEPFIX2 won't give dropped-out answer.
      PLOOP(pliter,pl) if(RADPL(pl)) pb[pl]=prtest[pl];
    }
    else if(*baseitermethod==QTYPMHD && (failreturnallowable==UTOPRIMGENWRAPPERRETURNNOFAIL)){
      PLOOP(pliter,pl) if(RADPL(pl)) pb[pl]=prtest[pl];
    }
    else{
      // then just leave pb as pb.
    }
  }





  /////////////////////
  //
  // see if uu0->p possible and desired.  Or just assign pp,pb,uu in terms of piin,Uiin,uu0 (based possibly upon optical depth).
  //
  //////////////////////
  FTYPE prtestUU0[NPR];
  PLOOP(pliter,pl) prtestUU0[pl]=pb[pl]; // initial guess

  // KORALTODO: Need to check if UFSET with no dUother fails.  How it fails, must allow since can do nothing better.  This avoids excessive attempts to get good solution without that failure!  Should speed-up things.  But what about error recovery?  If goes from CASE to no case!
  if(USEDUINRADUPDATE==3){
    // uu0 will hold original vector of conserved
    // here original means U[before fluxes, geometry, etc.] + dU[due to fluxes, geometry, etc. already applied and included in dUother]
    // This is required for stiff source term so immediately have balance between fluxes+geometry+radiation.  Otherwise, radiation diffuses.
    // I'm guessing that even though one uses RK2 or RK3, the first step generates large radiative velocities without any balanced source term since U isn't updated yet.  One would hope RK2 would recover on the final substep, but it doesn't!  In RK2, upon the final substep, that velocity is present for the radiation source term.  But it's also present for the fluxes!  That is, if there were no flux update on the final substep, then the source would have balanced the previous flux, but yet another flux is done, so there can be no balance.  This leads to a run-away velocity that would be similar to the \tau\sim 1 case.
    // NOTE: If this gives radiation or mhd failure, then less likely that will be actual solution since not even original uu0 has inversion.
    // Note that "q" isn't used in this function or used in function call, so don't have to update it here.

    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)
    struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
    if(0){
      // initialize counters
      newtonstats.nstroke=newtonstats.lntries=0;
    }
    else{
      ////////
      // be quick about checking this
      newtonstats.nstroke=newtonstats.lntries=0;
      // set inputs for errors, maxiters, etc.
      newtonstats.tryconv=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
      newtonstats.tryconvultrarel=1E-2*MAX(trueimptryconv,IMPOKCONVCONSTFORDEFAULT);
      newtonstats.extra_newt_iter=1;
      newtonstats.extra_newt_iter_ultrarel=2;
      newtonstats.mintryconv=MINTRYCONVFORMHDINVERSION; //IMPALLOWCONVCONSTFORDEFAULT;
      newtonstats.maxiter=MIN(trueimpmaxiter,IMPMAXITERFORDEFAULT);
    }
    //
    int failreturninversion;
    int finalstep = 1;
    int doradonly=0;
    int checkoninversiongas;
    int checkoninversionrad;
    // don't check since slows down code and could be good enough solution if original error says ok.
    checkoninversiongas=checkoninversionrad=0;
    //
    eomtypelocal=*eomtype; // stick with default choice so far
    failreturninversion=Utoprimgen_failwrapper(doradonly,radinvmod,showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, uu0, q, ptrgeom, dissmeasure, prtestUU0, &newtonstats);
    *nummhdinvsreturn++;


    // get pp0(uu0) so uu,uu0 properly associated with pp,pp0
    // set first pp that is p(uu0), assuming that can do inversion
    PLOOP(pliter,pl) ppfirst[pl]=pp[pl]=pp0[pl]=prtestUU0[pl];

  }
  else if(USEDUINRADUPDATE==2){
    // normal behavior, where uu0 has no known primitive inversion yet.
    // if entropy got solution and then energy is being done here, then pb is good entropy guess.
    PLOOP(pliter,pl){
      uu[pl] = uu0[pl];
      ppfirst[pl] = pp[pl] = pp0[pl] = pb[pl];
    }
  }
  else if(USEDUINRADUPDATE==1){


    // use tau as guess for which guess is best
    FTYPE tautot[NDIM]={0.0},tautotmax=0.0;
    calc_tautot(pb, ptrgeom, tautot, &tautotmax);

    // iterated, so keep as initial (i.e. previous full solution, not just initial+flux)
    PLOOP(pliter,pl){
      if(tautotmax>=TAUSWITCHPBVSPIIN){
        uu[pl] = Uiin[pl]; // assume somewhat static and Uiin is best guess
        ppfirst[pl] = pp[pl] = pp0[pl] = piin[pl];
      }
      else{
        uu[pl] = uu0[pl]; // assume dynamic and should use full uu0 to avoid pushing any errors into other fluid component.
        ppfirst[pl] = pp[pl] = pp0[pl] = pb[pl];
      }
    }
    // non-iterated, so keep as initial+flux (since all there is)
    PLOOP(pliter,pl){
      if(!(pl>=UU && pl<=U3 || RADPL(pl) || pl==ENTROPY)){
        uu[pl] = uu0[pl];
        ppfirst[pl] = pp[pl] = pp0[pl] = pb[pl];
      }
    }

  }
  else if(USEDUINRADUPDATE==0){
    // bad choice for non-iterated quantities, like uu[RHO] should be uu0[RHO]
    PLOOP(pliter,pl){
      uu[pl] = Uiin[pl];
      ppfirst[pl] = pp[pl] = pp0[pl] = piin[pl];
    }
  }
  else{
    // then (not recommended) just using Uiin as uu0 .  KORALNOTE: uu0 reset later, but intermediates would use this uu0.
    PLOOP(pliter,pl){
      uu[pl] = uu0[pl] = Uiin[pl];
      ppfirst[pl] = pp[pl] = pp0[pl] = piin[pl];
    }
  }  


  if(USEINPUTASGUESSIFERRORSMALL){
    // if error small enough and no fixups used with the radative inversion, override guess to use previous solution as the new guess
    if(errorabsreturn[WHICHERROR]<TRYHARDERFEEDGUESSTOL && RADINVBAD(radinvmodorig)==0){
      PLOOP(pliter,pl){
        pp[pl] = pp0[pl] = ppfirst[pl] = pb[pl]; // not just F(pb) anymore, holds better guess
        uu[pl] = uub[pl]; // not just Uiin or uu0 anymore, holds better guess
      }
    }
  }


  /////////////////////////////
  //
  // Fix u_g so entropy not smaller than guess
  // Only do this if u_g is just guess and is being solved for.
  // When eomtype==EOMCOLDGRMHD, then u_g is static so shouldn't modify as if was guess.
  //
  /////////////////////////////

  if(*eomtype!=EOMCOLDGRMHD && IMPMHDTYPE(implicititer)==1 && modprim==1){
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
  FTYPE errorabsf1[NUMERRORTYPES];
  set_array(errorabsf1,NUMERRORTYPES,MPI_FTYPE,BIG);
  FTYPE suberrorabsf1=BIG;
  FTYPE errorabsf3=BIG;
  FTYPE suberrorabsf3=BIG;
  // things that can happen that then no matter what other conditions, avoid damping or other attempts.
  int lowitererror=0;
  int earlylowerror=0;

  //  FTYPE trueimptryconv=IMPTRYCONV;
  FTYPE trueimptryconvabs=IMPTRYCONVABS;
  FTYPE trueimptryconvalt=IMPTRYCONVALT;
  FTYPE trueimpokconvabs=IMPOKCONVABS;
  FTYPE trueimpallowconvabs=IMPALLOWCONVABS;

  // best is over all damps as well
  FTYPE errorabsbest[NUMERRORTYPES];
  set_array(errorabsbest,NUMERRORTYPES,MPI_FTYPE,BIG);
  int failreturnbest=FAILRETURNGENERAL;
  int radinvmodbest=UTOPRIMRADFAILBAD1;


  ////////////////
  //
  // THE DAMP LOOP ITSELF
  //
  ////////////////
  int dampattempt;
  for(dampattempt=0;dampattempt<truenumdampattempts;dampattempt++){
    FTYPE DAMPFACTOR0;
    // dampattempt>0 refers to any attempt beyond the very first.  Uses iter from end of region inside this loop.
    if(dampattempt>0 && (lowitererror || earlylowerror || errorabsf1[WHICHERROR]<IMPTRYDAMPCONV && ACCEPTASNOFAILURE(failreturn) || failreturn==FAILRETURNMODESWITCH) ){ // try damping if any bad failure or not desired tolerance when damping
      if(dampattempt>=2 && failreturn!=FAILRETURNMODESWITCH){ // dampattempt>=2 refers to attempts with at least 1 damp attempt
        if(debugfail>=2) dualfprintf(fail_file,"Damping worked to reach desired tolerance: errorabsf1=%g errorallabsf1=%g (IMPTRYDAMPCONV=%g), so should have lower error: dampattempt=%d iter=%d ijk=%d %d %d\n",errorabsf1[0],errorabsf1[1],IMPTRYDAMPCONV,dampattempt,iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);
      }
      break; // if didn't hit problem, so no need to damp since got tolerance requested or returned because will just switch to another scheme.
    }
    else{
      // control factor by which step Newton's method.
      DAMPFACTOR0=1.0/pow(2.0,(FTYPE)(dampattempt));
      if(dampattempt>0) if(debugfail>=2) dualfprintf(fail_file,"Trying dampattempt=%d DAMPFACTOR0=%g failreturn=%d errorabsf1=%g errorallabsf1=%g iter=%d ijk=%d %d %d\n",dampattempt,DAMPFACTOR0,failreturn,errorabsf1[0],errorabsf1[1],iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);

      // start fresh
      iter=debugiter=0;
      failreturn=FAILRETURNNOFAIL; // default is no failure
      mathfailtype=-1; // indicates not set
      {
        //        PLOOP(pliter,pl) uu0[pl]=uu0dampbackup[pl];
      }
      if(errorabsf1[WHICHERROR]<TRYHARDERFEEDGUESSTOL){
        // then keep pp and uu as better starting point
      }
      else{
        PLOOP(pliter,pl){
          pp[pl]=ppdampbackup[pl];
          uu[pl]=uudampbackup[pl];
          uu0[pl]=uu0dampbackup[pl];
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
    FTYPE errorabspf1[NUMPRIORERRORS];
    set_array(errorabspf1,NUMPRIORERRORS,MPI_FTYPE,BIG);
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
      PLOOP(pliter,pl) f1reportlist[iterlist][pl]=BIG;
      PLOOP(pliter,pl) f1list[iterlist][pl]=BIG;
      errorabsf1list[iterlist]=BIG;
      errorallabsf1list[iterlist]=BIG;
      realiterlist[iterlist]=-1;
      DLOOP(jj,kk) jaclist[iterlist][jj][kk]=BIG;
      implicititerlist[iterlist]=-1;
      implicitferrlist[iterlist]=-1;
      fracdamplist[0]=fracdtuu0;
      fracdamplist[1]=fracdtG;
      fracdamplist[2]=0;
    }
    // DEBUG:
    errorabsreturn[0]=errorabsreturn[1]=BIG;

    // whether holding as positive and outher counts
    // or things that happen that mean need to break out of this attempt, but not break out of damping loop.
    int holdingaspositive=0,iterhold=0;
    int countholdpositive=0;
    int countugnegative=0;
    int countbadenergy=0;
    int counterrorrose=0;
    int priorerrorscount=0;
    int canbreak=0;
    int notfinite=0;
    int convreturnf3limit=0;
    int notholding=1;
    FTYPE DAMPFACTOR;
    int numjumpchecks=0;

    ////////////////////////////////
    //
    // IMPLICIT LOOP ITSELF
    //
    ////////////////////////////////

    do{
      debugiter++; // never skips, goes every step
      debugiteratteempts[dampattempt]=debugiter;
      iter++; // below might skip some iter, used to control which equations used

      //      dualfprintf(fail_file,"iter=%d debugiter=%d\n",iter,debugiter);

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
      define_method(iter, &eomtypelocal, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
      get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);


      ////////////////////
      //
      // things to reset each iteration start for any steps
      //
      ///////////////////
      canbreak=0; // reset each start of iteration
      convreturnf3limit=0;
      notholding=1; // default is no hold
      earlylowerror=0; // reset each iteration
      lowitererror=0;
    
      if(iter>10){ // KORALTODO: improve upon this later.  Only matters if not doing PMHD method
        // assume trying hard and failing to work, then allow CASE radiation errors
        failreturnallowable=failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
      }

      if(CHANGEDAMPFACTOR==1){
        if(trueimpmaxiter==IMPMAXITERQUICK && iter>IMPMAXITERQUICK/2){
          DAMPFACTOR=0.5*DAMPFACTOR0;
        }
        else if(trueimpmaxiter==IMPMAXITERLONG && iter>MIN(IMPMAXITERLONG/2,20)){
          DAMPFACTOR=0.5*DAMPFACTOR0;
        }
        else DAMPFACTOR=DAMPFACTOR0;
      }
      else if(CHANGEDAMPFACTOR==2){
        if(itermode==ITERMODESTAGES && iter==BEGINMOMSTEPS) DAMPFACTOR=0.5*DAMPFACTOR0; else DAMPFACTOR=DAMPFACTOR0;
        if(itermode==ITERMODESTAGES && iter==BEGINENERGYSTEPS) DAMPFACTOR=0.5*DAMPFACTOR0; else DAMPFACTOR=DAMPFACTOR0;
        if(itermode==ITERMODESTAGES && iter==BEGINFULLSTEPS) DAMPFACTOR=0.25*DAMPFACTOR0; else DAMPFACTOR=DAMPFACTOR0; // 0.25 to be very careful since can change energy too much when suddenly turning back on momentum
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
        countugnegative=0;
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
      int convreturnf1;
      for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
        
        //        dualfprintf(fail_file,"iter=%d debugiter=%d f1iter=%d\n",iter,debugiter,f1iter);


        int whichcall=FIMPLICITCALLTYPEF1;
        eomtypelocal=*eomtype; // re-chose default each time.  If this reduces to a new eomtype, then Jacobian will stick with that for consistency!
        int goexplicit;
        int dimtypef=DIMTYPEFCONS;

        // get original baseitermethod
        int baseitermethodorig=*baseitermethod;
        // use ppppp and uuppp as backup since one previous step is often off alot even if not already hitting point at which one should change baseitermethod
        failreturnferr=f_implicit(allowbaseitermethodswitch, iter, f1iter, failreturnallowableuse, whichcall, impepsjac, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, baseitermethod, fracenergy, dissmeasure, radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, ppppp, pp, piin, uuppp, Uiin, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, f1report, &goexplicit, &errorabsf1[0], &errorabsf1[1], WHICHERROR, &convreturnf1, nummhdinvsreturn); // modifies uu and pp, f1poret, goexplicit, errorabsf1[0,1], convreturnf1
        // if baseitermethod changed, then need to deal with fracdtuu0 that may have different properties.  So while we don't necessarily go back too far on pp,uu, we need to pretend starting over with how we deal with uu0.
        if(*baseitermethod!=baseitermethodorig) iter=1;


        // setup method and signs in case changed baseitermethod
        define_method(iter, &eomtypelocal, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
        get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);


        // see if 4-force negligible
        if(goexplicit){
          // immediate return.
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
          {
            FTYPE uuiterback[NPR];
            PLOOP(pliter,pl) uuiterback[pl] = uu[pl]; // hold actual iterated quantities
            PLOOP(pliter,pl) uu[pl] = uu0[pl]; // "iterate" non-iterated quantities
            JACLOOP(ii,startjac,endjac) uu[irefU[ii]] = uuiterback[irefU[ii]]; // overwrite with actual previous step for iterated quantities
          }
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


        *nummhdstepsreturn = (int)(nstroke-nstrokeorig); // get number of mhd inversion steps
        if(*nummhdstepsreturn>MAXMHDSTEPS){
          failreturn=0;
          break;
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
        *nummhdstepsreturn = (int)(nstroke-nstrokeorig); // get number of mhd inversion steps
        if(f1iter==MAXF1TRIES || *nummhdstepsreturn>MAXMHDSTEPS){

          if(f1iter==MAXF1TRIES) if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES=%d (nummhdstepsreturn=%d: %d %d trueimpmaxiter=%d): fracdtuu0=%g nstep=%ld steppart=%d ijk=%d %d %d : iter=%d eomtype=%d baseitermethod=%d failreturn=%d\n",MAXF1TRIES,*nummhdstepsreturn,nstrokeorig,nstroke,trueimpmaxiter,fracdtuu0,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter,eomtypelocal,*baseitermethod,failreturn);
          if(*nummhdstepsreturn>MAXMHDSTEPS) if(debugfail>=2) dualfprintf(fail_file,"Reached MAXMHDSTEPS=%d: fracdtuu0=%g nstep=%ld steppart=%d ijk=%d %d %d : iter=%d eomtype=%d baseitermethod=%d failreturn=%d\n",*nummhdstepsreturn,fracdtuu0,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter,eomtypelocal,*baseitermethod,failreturn);
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
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d ii=%d pl=%d f1=%g\n",ptrgeom->i,ii,pl,f1[pl]);
      }// else if f1 calculation didn't fail.


    
    
      ////////
      //
      // get type of EOMTYPE
      //
      ////////
      int eomcond=(eomtypelocal==EOMGRMHD || eomtypelocal==EOMDEFAULT && EOMDEFAULT==EOMGRMHD);


      /////////
      //
      // see if should check convergence or check how solution is behaving.
      //
      /////////
      checkconv=1;
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracuup!=1.0){
      //      if(fabs(fracuup-1.0)>10.0*NUMEPSILON){
      if(fracuup<1.0){
        // try increasing amount of uu or pp used
        checkconv=0;
      }
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracdtuu0!=1.0){
      //      if(fabs(fracdtuu0-1.0)>10.0*NUMEPSILON && changeotherdt){
      if(fracdtuu0<1.0){
        checkconv=0;
      }
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracdtG!=1.0){
      //      if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
      if(fracdtG<1.0 && changeotherdt){
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
      //      int convreturnf1=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs, realdt, DIMTYPEFCONS,eomtypelocal , *radinvmod, itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1[0],&errorabsf1[1],WHICHERROR);
      // but don't break, since need to iterate a bit first and check |dU/U| and need to see if checkconv==1
      FTYPE numsub=0;
      suberrorabsf1=0.0;  JACLOOPSUBERROR(jj,startjac,endjac){
        suberrorabsf1     += fabs(f1report[erefU[jj]]); // e.g. may only be energy error or only momentum error.
        numsub += 1.0;
      }



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
        PLOOP(pliter,pl) f1reportlist[debugiter][pl]=f1report[pl];
        PLOOP(pliter,pl) f1list[debugiter][pl]=f1[pl];
        errorabsf1list[debugiter]=errorabsf1[0];
        errorallabsf1list[debugiter]=errorabsf1[1];
        realiterlist[debugiter]=iter;
        DLOOP(jj,kk) jaclist[debugiter][jj][kk]=BIG; // just default dummy value in case break
        implicititerlist[debugiter]=implicititer;
        implicitferrlist[debugiter]=implicitferr;
        fracdamplist[0]=fracdtuu0;
        fracdamplist[1]=fracdtG;
        fracdamplist[2]=DAMPFACTOR;
      }



      ////////////////
      //
      // If error dropped below tolerance for this sub-matrix iteration mode, then continue to next level.
      //
      /////////////////
      // check if energy only iteration has error that has dropped below tolerance, then can move on to 
      if(itermode==ITERMODESTAGES && iter>=BEGINMOMSTEPS && iter<=ENDMOMSTEPS){
        FTYPE errorabs3=0.0; int testi;
        for(testi=1;testi<=3;testi++) errorabs3 += fabs(f1report[irefU[testi]]);
        FTYPE errorneed=((FTYPE)(3+2)*trueimptryconv);
        if(errorabs3<errorneed){
          if(iter<=ENDMOMSTEPS){ iter=BEGINENERGYSTEPS-1; continue;} // force as if already doing energy steps.  If already next iteration is to be this energy step, then no skipping needed.
          // continue assumes not triggered when iter>trueimpmaxiter
        }
      }
      // check if energy only iteration has error that has dropped below tolerance, then can move on to 
      if(itermode==ITERMODESTAGES && iter>=BEGINENERGYSTEPS && iter<=ENDENERGYSTEPS){
        //        if(f1report[irefU[0]]==BIG){ dualfprintf(fail_file,"FUDGE\n"); }
        // SUPERGODMARK: valgrind says belw is undefined, but don't see it.
        FTYPE errorabs1=0.0; int testi;
        for(testi=1;testi<=1;testi++) errorabs1 += fabs(f1report[irefU[testi]]);
        FTYPE errorneed=((FTYPE)(1+2)*trueimptryconv);
        if(errorabs1<errorneed){
          if(iter<=ENDENERGYSTEPS){ iter=BEGINFULLSTEPS-1; continue;} // force as if already doing normal steps.  If already next iteration is to be normal step, no need to skip.
          // continue assumes not triggered when iter>trueimpmaxiter
        }
      }







      ///////////////////////////
      //
      //  PRE NEWTON ADJUSTMENTS
      //
      ///////////////////////////

#define ERRORJUMPCHECK 0 // no longer needed with error trend check, and expensive if u_g or Erf try to be <0 since hits continue.
// whether to check if error jumps up, between last and current step, in irefU[0] -- u_g for QTYPMHD method.  If so, backs-up step a bit for all quantities iterated and try to get error again
#define BUFFERITER 0 // how many iterations to wait until start to check how error is doing.  When switching iteration methods, error will often rise initially in f1[0], but that's ok.  But can temper jump by bridging used pp,uu, so ok to keep as 0 perhaps.
      // check if doing energy stepping and error jumped up too much
#define NUMJUMPCHECKSMAX 5 // must limit, else if really drops-out and can't help, need to just accept.
      if(ERRORJUMPCHECK && numjumpchecks<NUMJUMPCHECKSMAX){
        if(iter>=BEGINENERGYSTEPS+BUFFERITER && iter<=ENDENERGYSTEPS || iter>=BEGINFULLSTEPS+BUFFERITER && iter<=ENDFULLSTEPS){// now all steps beyond energy
          if(
             (fabs(f1[erefU[0]]/f1p[erefU[0]])>FACTORBADJUMPERROR && fabs(f1report[erefU[0]])>trueimptryconv)
             || (pp[URAD0]<10.0*ERADLIMIT) // expensive
             || (pp[UU]<10.0*UUMINLIMIT) // expensive
             ){
            // then pseudo-bisect (between zero and previous ok error case)
            if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"pseudo-bisect: iter=%d f1=%g f1p=%g pp=%g ppp=%g  pppp=%g  ppppp=%g\n",iter,f1[erefU[0]],f1p[erefU[0]],pp[irefU[0]],ppp[irefU[0]],pppp[irefU[0]],ppppp[irefU[0]]);
            if(0){
              // doens't make sense in general
              //          pp[irefU[0]] = ppp[irefU[0]] = pppp[irefU[0]] = 0.5*(fabs(ppppp[irefU[0]]));
              pp[irefU[0]] = ppp[irefU[0]] = 0.5*fabs(pppp[irefU[0]]);
              uu[irefU[0]] = uup[irefU[0]] = 0.5*fabs(uupp[irefU[0]]);
            }
            else{
              // half-way between current (pp and ppp are same current primitve) and last primitive
              // uu and pp won't be consistent, but when get to f_implicit(), as continue forces, this will be done.
              PLOOP(pliter,pl) pp[pl] = 0.75*pppp[pl] + 0.25*ppp[pl];
              PLOOP(pliter,pl) uu[pl] = 0.75*uupp[pl] + 0.25*uup[pl];
            }
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
              DLOOP(jj,kk) jaclist[debugiter][jj][kk]=iJ[irefU[jj]][erefU[kk]];
              implicititerlist[debugiter]=implicititer;
              implicitferrlist[debugiter]=implicitferr;
              fracdamplist[0]=fracdtuu0;
              fracdamplist[1]=fracdtG;
              fracdamplist[2]=DAMPFACTOR;
            }
            // need to get new error function so can take step based upon this as reference!
            if(iter>trueimpmaxiter){
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"iter=%d>%d\n",iter,trueimpmaxiter);
            }
            else{
              // continue assumes not triggered when iter>trueimpmaxiter
              numjumpchecks++;
              if(numjumpchecks<NUMJUMPCHECKSMAX) continue; // head to start of loop to iter++ and get new error function.
            }
          }
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
        if(STOPIFVERYLOWERROR && errorabsf1[WHICHERROR]<=LOCALPREIMPCONVABS){
          earlylowerror=1;
        }
        // see if error on iterated quantities has already gone near/below machine precision.  If so, can't do any better with total error, so stop.
        if(STOPIFITERLOWERROR && WHICHERROR==1 && errorabsf1[0]<=LOCALPREIMPCONVABS){
          lowitererror=1; // FUCK: makes use urad alot and why?
        }

        //////////////
        //
        // try to get best solution (should be based upon immediate f_error_check(uu,pp,errorabsf1)
        //
        // store error and solution in case eventually lead to max iterations and actually get worse error
        // f1 based
        //
        //////////////
        errorabsbest[0]=0.0; JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsbest[0] += fabs(lowestf1report[erefU[jj]]);
        errorabsbest[1]=0.0; JACLOOPSUPERFULL(pliter,pl,*eomtype,*baseitermethod,*radinvmod) errorabsbest[1] += fabs(lowestf1report[pl]);
        if(errorabsbest[WHICHERROR]>errorabsf1[WHICHERROR] && isfinite(errorabsf1[WHICHERROR]) && (itermode==ITERMODECOLD || pp[RHO]>0.0 && pp[UU]>0.0 && pp[PRAD0]>0.0)){
          PLOOP(pliter,pl) bestuu[pl]=uu[pl];
          PLOOP(pliter,pl) bestpp[pl]=pp[pl];
          PLOOP(pliter,pl) lowestf1[pl]=f1[pl];
          PLOOP(pliter,pl) lowestf1norm[pl]=f1norm[pl];
          PLOOP(pliter,pl) lowestf1report[pl]=f1report[pl];
          errorabsbest[0]=errorabsf1[0];
          errorabsbest[1]=errorabsf1[1];
          radinvmodbest=*radinvmod;
          gotbest=1;
        }


        if(earlylowerror){
          if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Early low error=%g %g : iter=%d\n",errorabsf1[0],errorabsf1[1],iter);
          //  not failure.
          break;
        }
        if(lowitererror){
          if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Low iter error=%g %g : iter=%d\n",errorabsf1[0],errorabsf1[1],iter);
          //  not failure.
          break;
        }


#define CHECKDECREASE0 5
#define CHECKDECREASEAVGNUM 3 // should be less than CHECKDECREASE0
#define CHECKDECFACTOR (0.5) // i.e. should drop by this factor compared to older average
        if(debugiter>=BEGINNORMALSTEPS+CHECKDECREASE0){
          int ci;
          FTYPE avgerror[NUMERRORTYPES];
          avgerror[0]=avgerror[1]=0.0;
          for(ci=0;ci<CHECKDECREASEAVGNUM;ci++) avgerror[0] += errorabsf1list[debugiter-CHECKDECREASE0+ci]; // GODMARK: uses debug info
          for(ci=0;ci<CHECKDECREASEAVGNUM;ci++) avgerror[1] += errorallabsf1list[debugiter-CHECKDECREASE0+ci]; // GODMARK: uses debug info
          avgerror[0]/=(CHECKDECREASEAVGNUM);
          avgerror[1]/=(CHECKDECREASEAVGNUM);
          FTYPE currenterror[NUMERRORTYPES];
          currenterror[0]=errorabsf1list[debugiter];
          currenterror[1]=errorallabsf1list[debugiter];
          //
          // check both errors to ensure they are decreasing
          int cond1=(currenterror[0]>CHECKDECFACTOR*avgerror[0] || currenterror[1]>CHECKDECFACTOR*avgerror[1])&&(WHICHERROR==1);
          int cond2=(currenterror[0]>CHECKDECFACTOR*avgerror[0])&&(WHICHERROR==0);
          if(cond1 || cond2){

            if(havebackup){
              failreturn=FAILRETURNMODESWITCH; mathfailtype=80;
              if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"SWITCHING MODE: Detected did not decrease error\n");
              // if want to ensure should have gotten solution, should still report
                }
            else{
              // then aborting due to error alone even without backup
              canbreak=4; // use same canbreak as below
              mathfailtype=81;
              if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"Aborting even without backup because error not decreasing\n");
              // no failure or switch.
            }
            break;

          }// end if current error too large compared to older average error
        }// end if large enough iterations so can check how error is trending, to avoid many iterations when error is not dropping enough to matter.


        // check if error repeatedly rises
        // do this even if damping, since damping should only help get continuous reduction in error.
        if(NUMNOERRORREDUCE && iter>=BEGINNORMALSTEPS){
          FTYPE errorneed=((FTYPE)(numsub+2)*trueimptryconv);
          FTYPE errorneedalt=(MAX(1E-8,errorneed));
          if(iter>NUMNOERRORREDUCE0 && suberrorabsf1>errorneed){ // no need to do this if actually error is below desired tolerance,  hence second argument
            if(suberrorabsf1>=errorabspf1[0]) counterrorrose++;
            int allowedtoabort=(havebackup || havebackup==0 && ABORTBACKUPIFNOERRORREDUCE==1 && suberrorabsf1<errorneedalt);
            //+int allowedtoabort=(havebackup || havebackup==0 && ABORTBACKUPIFNOERRORREDUCE==1 && suberrorabsf1<errorneedalt || MODEMETHOD==MODEPICKBEST || MODEMETHOD==MODEPICKBESTSIMPLE  || MODEMETHOD==MODEPICKBESTSIMPLE2);
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
                  if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"Aborting even without backup because error oscillated (iter=%d) and suberrorabsf1=%g errorabsf1[0,1]=%g %g\n",iter,suberrorabsf1,errorabsf1[0],errorabsf1[1]);
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
            FTYPE errorneed=((FTYPE)(numsub+2)*trueimptryconv);
            if(suberrorabsf1/erroraverage>changerequired && suberrorabsf1>errorneed){

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
        if(errorabsf1[WHICHERROR]<ERRORFORIMPEPSSMALL) impepsjac=IMPEPSSMALL;
        else impepsjac=IMPEPSLARGE;
        int dimtypef=DIMTYPEFCONS;
        int failreturniJ=get_implicit_iJ(allowbaseitermethodswitch, failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, baseitermethod, fracenergy, dissmeasure, impepsjac, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, iter, errorabsf1[0], errorabsf1[1], WHICHERROR, dimtypef, dimfactU, Uiin, uu, uup, uu0, piin, pp, ppp, fracdtG, realdt, ptrgeom, q, f1, f1norm, iJ, nummhdinvsreturn);

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
          JACLOOP2D(iii,jjj,startjac,endjac)  dualfprintf(fail_file,"iJ[i %d][e %d]=%g\n",iii,jjj,iJ[irefU[iii]][erefU[jjj]]);
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
          JACLOOP2D(ii,jj,startjac,endjac) uu[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];

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

          if(NEWJONHOLDPOS==0){
            PLOOP(pliter,pl) pp[pl]=ppp[pl];
            JACLOOP2D(ii,jj,startjac,endjac){
              pp[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"added to ppp=%21.15g ii=%d jj=%d irefU=%d an amount of negative %21.15g\n",ppp[irefU[ii]],ii,jj,irefU[ii],DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]]);
            }
          }
          else{
            // if u_g (or Erf)<0, then shift entire jacobian to make smaller changes

            FTYPE dpp[NPR]={0.0};
            PLOOP(pliter,pl) pp[pl]=ppp[pl];
            JACLOOP2D(ii,jj,startjac,endjac){
              dpp[irefU[ii]] += -DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"added to ppp=%21.15g ii=%d jj=%d irefU=%d an amount of negative %21.15g\n",ppp[irefU[ii]],ii,jj,irefU[ii],DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]]);
            }
            JACLOOP(ii,startjac,endjac){
              pp[irefU[ii]] += dpp[irefU[ii]];
            }

            if(startjac==0){
              ii=0;
              if(pp[irefU[ii]]<=0.0){
#if(0)
                FTYPE ppbad[NPR];
                PLOOP(pliter,pl) ppbad[pl]=pp[pl];
                // rescale iJ so u_g -> u_g0/2 instead (or for Erf)
                PLOOP(pliter,pl) pp[pl]=ppp[pl];
                FTYPE REDAMP=fabs(+0.5*ppp[UU]/dpp[UU])*DAMPFACTOR; // KORALTOD: inverted damp.  Sometimes works, sometimes bad.
                JACLOOP2D(ii,jj,startjac,endjac){
                  pp[irefU[ii]] += REDAMP*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
                }
#else
                pp[irefU[ii]]=ppp[irefU[ii]]; // super hold, works just as well as inverted damp.
#endif
              }
            }

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
            DLOOP(jj,kk) jaclist[debugiter][jj][kk]=iJ[irefU[jj]][erefU[kk]];
            implicititerlist[debugiter]=implicititer;
            implicitferrlist[debugiter]=implicitferr;
            fracdamplist[0]=fracdtuu0;
            fracdamplist[1]=fracdtG;
            fracdamplist[2]=DAMPFACTOR;
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






          ////////////
          //
          // check if u_g,Erf<0.  Do even if RAMESHFIXEARLYSTEPS going or even if checkconv==0
          //
          ////////////
          //          FTYPE umin=10.0*calc_PEQ_ufromTrho(TEMPMIN,fabs(pp[RHO]));
          if(pp[irefU[0]]<=0.0 && (IMPPTYPE(implicititer)&& *baseitermethod!=QTYENTROPYPMHD)){ // don't consider implicititer==QTYENTROPYPMHD since S can be positive or negative.  Would only be unphysical or absolute-limited if the related u_g<0 or rho<0.
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
                  if(badenergy(ptrgeom,pp,pborig)) countbadenergy++;
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

          ////////////
          //
          // check if u_g<0 too often and know have entropy backup
          //
          ////////////

          if(IMPRADTYPEBASE(*baseitermethod)){//FUCK, might be too aggressive for PMHD method.
            // abort if too often hit negative u_g or Erf and know have entropy backup method.
#define NUMUGNEGENERGY (3)

            if(pp[UU]<=10.0*UUMINLIMIT && *baseitermethod!=QTYENTROPYPMHD && iter>=BEGINNORMALSTEPS){
              countugnegative++;
            }
            else{
              //            dualfprintf(fail_file,"NOTNEG: countugnegative=%d\n",countugnegative);
            }

            if(countugnegative>=NUMUGNEGENERGY && eomcond && (MODEMETHOD==MODEPICKBEST || MODEMETHOD==MODEPICKBESTSIMPLE  || MODEMETHOD==MODEPICKBESTSIMPLE2) && iter>=BEGINNORMALSTEPS){
              if(DOFINALCHECK){
                if(debugfail>=DEBUGLEVELIMPSOLVER) dualfprintf(fail_file,"2Unable to hold off u_g<0, setting canbreak=1 and letting finalchecks confirm error is good or bad.  iter=%d\n",iter);
                canbreak=5; // just break since might be good (or at least allowable) error still.  Let final error check handle this.
                // not fail.
                mathfailtype=89;
                break;
              }
              else if(havebackup){
                failreturn=FAILRETURNMODESWITCH; mathfailtype=90;
                if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"2SWITCHING MODE: Detected bad u_g\n");
                break;
              }
              else{
                // then full failure
                failreturn=FAILRETURNGENERAL; mathfailtype=10;
                break;
              }
            }
          }

          FTYPE f3[NPR]={0},f3norm[NPR]={0};
          if(POSTNEWTONCONVCHECK==1 && notholding==1){
            int dimtypef3;
            /////////
            //
            // test convergence after Newton step
            // test convergence using |dU/U|
            // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
            // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
            /////////
            // 0 = conserved R^t_\nu type, 1 = primitive R^{ti} type
            if(IMPUTYPE(implicititer)){ // still considering iterate error
              dimtypef3=DIMTYPEFCONS;
              JACLOOPALT(ii,startjac,endjac){
                f3[erefU[ii]]=(uu[irefU[ii]]-uup[irefU[ii]]);
                f3norm[erefU[ii]]=fabs(uu[irefU[ii]])+fabs(uup[irefU[ii]]);
              }
            }
            else if(IMPPTYPE(implicititer)){ // still considering iterate error
              dimtypef3=DIMTYPEFPRIM;
              JACLOOPALT(ii,startjac,endjac){
                f3[erefU[ii]]=(pp[irefU[ii]]-ppp[irefU[ii]]);
                f3norm[erefU[ii]]=fabs(pp[irefU[ii]])+fabs(ppp[irefU[ii]]);
              }
            }
  
            // store error and solution in case eventually lead to max iterations and actually get worse error
            // f_error_check(uu0,uu) is ok to use since it just normalizes error
            int convreturnf3;
            FTYPE errorallabsf3; // not used
            convreturnf3=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs, realdt, dimtypef3,eomtypelocal ,*radinvmod, itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,f3,f3norm,f3report,Uiin,uu0,uu,ptrgeom,&errorabsf3,&errorallabsf3,WHICHERROR);
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
      //
      // revert any back-ups if ok to do so
      //
      // NOTEMARK: Can only modify uu0 or frac's after they are used consistently to get f1, iJ, take step, and check convergence with error function
      //
      /////////
      changeotherdt=1;
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracuup!=1.0){
      //      if(fabs(fracuup-1.0)>10.0*NUMEPSILON){
      if(fracuup<1.0){
        // try increasing amount of uu or pp used
        fracuup*=RADDAMPUNDELTA;
        fracuup=MIN(1.0,fracuup);
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
        {
          FTYPE uuiterback[NPR];
          PLOOP(pliter,pl) uuiterback[pl] = uu[pl]; // hold actual iterated quantities
          PLOOP(pliter,pl) uu[pl] = uu0[pl]; // "iterate" non-iterated quantities
          JACLOOP(ii,startjac,endjac) uu[irefU[ii]] = uuiterback[irefU[ii]]; // overwrite with actual previous step for iterated quantities
        }
        // KORALNOTE: No need to get pp0, since never used.  uu0 only used in error function.
      }
      ////////////////////////////////////////////////////////////////////////////
      //    if(fracdtG!=1.0){
      //      if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
      if(fracdtG<1.0 && changeotherdt){
        // try increasing amount of G applied
        fracdtG*=RADDAMPUNDELTA;
        fracdtG=MIN(1.0,fracdtG);
      }      


      /////////
      // see if took too many Newton steps or not finite results
      /////////
      *nummhdstepsreturn = (int)(nstroke-nstrokeorig); // get number of mhd inversion steps
      if(iter>trueimpmaxiter || *nummhdstepsreturn>MAXMHDSTEPS || notfinite ){
        if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"iter=%d>%d mhdsteps=%d>%d or notfinite=%d\n",iter,trueimpmaxiter,*nummhdstepsreturn,MAXMHDSTEPS,notfinite);
        //      failreturn=FAILRETURNGENERAL; // no don't fail, might be allowable error.
        break;
      }

    }// end do
    while(1);








    ////////////
    //
    // once done iterating, regardless of result or failure, need to ensure non-iterated quantities are normal
    //
    /////////////
    fracdtuu0=1.0;
    //
    PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0
    PLOOP(pliter,pl){
      if(!(pl>=UU && pl<=U3 || RADPL(pl) || pl==ENTROPY)){
        uu[pl] = uu0[pl];
      }
    }
    // regardless, invert non-iterated quantities
    // 1) trivially invert field
    PLOOPBONLY(pl) bestpp[pl]=bestuu[pl]=pp[pl] = uu0[pl];
    // 2) trivially invert to get rho assuming q up-to-date
    //    pp[RHO]= uu0[RHO]/q->ucon[TT]; // q might not be up-to-date
    // 3) Invert other scalars (only uses uu[RHO], not pp[RHO])
    extern int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr);
    invert_scalars(ptrgeom, uu,pp);
    PLOOP(pliter,pl){
      if(pl!=RHO && pl<UU && pl>U3 && pl<B1 && pl>B3 & RADPL(pl)==0 && pl!=ENTROPY){
        bestpp[pl]=pp[pl];
        bestuu[pl]=uu[pl];
      }
    }
    //
    // end trivial inversion
    //
    //////////////

  


    //////////////
    //
    // if no failure, then see if really still failed or not with some final checks.
    //
    ///////////////
    int failreturnf=0;
    //    if(failreturn==0 && (earlylowerror==0 && lowitererror==0) ){
    if(failreturn==0 && earlylowerror==0){ // still check if lowitererror==1 because might have best total error at another iter.


      // if didn't fail, shouldn't need to set fracdtuu0=1, but do so just in caes.
      fracdtuu0=1.0;


      int fakeiter;
      for(fakeiter=0;fakeiter<=0;fakeiter++){

        int convreturn=1,convreturnok=1,convreturnallow=1; // default is solution is acceptable.


        if(DOFINALCHECK){
          //////////////////////////
          //
          // check and get error for last iteration or any mods from above that are post-iteration
          //
          // uses final uu0 with fracdtuu0=1 to ensure really updated non-iterated quantities correctly
          //
          // The call to f_implicit() also ensures uu is consistent with new pp
          //
          ////////////////////////
          int whichcall=FIMPLICITCALLTYPEFINALCHECK;
          //  eomtypelocal=*eomtype; // re-chose default each time. No, stick with what f1 (last call to f1) chose
          int goexplicitfake; // not used here
          int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive (u,v^i) type, i.e. v^i has no energy density term
          int fakef1iter=-1;
          failreturnf=f_implicit(allowbaseitermethodswitch, iter,fakef1iter,failreturnallowableuse, whichcall,impepsjac,showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, baseitermethod, fracenergy, dissmeasure, radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, pp, pp, piin, uu, Uiin, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, f1report, &goexplicitfake, &errorabsf1[0], &errorabsf1[1], WHICHERROR, &convreturn, nummhdinvsreturn); // modifies uu and pp and q and f1report and goexplicitfake and errorabsf1
          // radinvmod contains whether radiative inversion modified process.

          //          convreturn=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv,trueimptryconvabs,realdt,dimtypef,eomtypelocal,*radinvmod, itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1,WHICHERROR);
          convreturnallow=(errorabsf1[WHICHERROR]<trueimpallowconvabs);
          convreturnok=(errorabsf1[WHICHERROR]<trueimpokconvabs);
          if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"DOFINALCHECK: convreturn=%d convreturnok=%d convreturnallow=%d (IMPOKCONV=%g IMPALLOWCONV=%g) f1report: %g %g %g %g : %g %g\n",convreturn,convreturnok,convreturnallow,IMPOKCONV,IMPALLOWCONV,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],errorabsf1[0],errorabsf1[1]);
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
          convreturn=f_error_check(showmessages, showmessagesheavy, iter, trueimptryconv, trueimptryconvabs,realdt,dimtypef,eomtypelocal,*radinvmod,itermode,*baseitermethod,fracenergy,dissmeasure,dimfactU,pp,piin,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom,&errorabsf1[0],&errorabsf1[1],WHICHERROR);
          convreturnallow=(errorabsf1[WHICHERROR]<trueimpallowconvabs);
          convreturnok=(errorabsf1[WHICHERROR]<trueimpokconvabs);
        }



        if(GETBEST){
      
          if(gotbest){
            // f1-based
            // using old uu,uup, but probably ok since just helps normalize error
            errorabsf1[0]=0.0;   JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsf1[0]   += fabs(f1report[erefU[jj]]);
            errorabsf1[1]=0.0;   JACLOOPSUPERFULL(pliter,pl,*eomtype,*baseitermethod,*radinvmod)  errorabsf1[1]   += fabs(f1report[pl]);
            errorabsbest[0]=0.0; JACLOOPFULLERROR(itermode,jj,startjac,endjac) errorabsbest[0] += fabs(lowestf1report[erefU[jj]]);
            errorabsbest[1]=0.0; JACLOOPSUPERFULL(pliter,pl,*eomtype,*baseitermethod,*radinvmod)  errorabsbest[1] += fabs(lowestf1report[pl]);

            // see if should revert to prior best
            if(errorabsbest[WHICHERROR]<errorabsf1[WHICHERROR] || !isfinite(errorabsf1[WHICHERROR]) ){
              if(showmessages && debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"Using best: %g %g : %g %g\n",errorabsf1[0],errorabsf1[1],errorabsbest[0],errorabsbest[1]);

              PLOOP(pliter,pl) uu[pl]=bestuu[pl];
              PLOOP(pliter,pl) pp[pl]=bestpp[pl];
              errorabsf1[0]=errorabsbest[0];
              errorabsf1[1]=errorabsbest[1];
              PLOOP(pliter,pl) f1[pl] = lowestf1[pl];
              PLOOP(pliter,pl) f1norm[pl] = lowestf1norm[pl];
              PLOOP(pliter,pl) f1report[pl] = lowestf1report[pl];
              *radinvmod = radinvmodbest;

              // get whether converged
              convreturn=(errorabsf1[WHICHERROR]<trueimptryconvabs);
              convreturnallow=(errorabsf1[WHICHERROR]<trueimpallowconvabs);
              convreturnok=(errorabsf1[WHICHERROR]<trueimpokconvabs);
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"GETBEST: convreturn=%d convreturnok=%d convreturnallow=%d (IMPOKCONV=%g IMPALLOWCONV=%g) f1report: %g %g %g %g : %g %g\n",convreturn,convreturnok,convreturnallow,IMPOKCONV,IMPALLOWCONV,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],errorabsf1[0],errorabsf1[1]);
            }
            else{
              if(debugfail>=DEBUGLEVELIMPSOLVERMORE) dualfprintf(fail_file,"gotbest=%d but errorabsbest=%g %g while errorabsf1=%g %g\n",gotbest,errorabsbest[0],errorabsbest[1],errorabsf1[0],errorabsf1[1]);
              PLOOP(pliter,pl) bestuu[pl]=uu[pl];
              PLOOP(pliter,pl) bestpp[pl]=pp[pl];
              errorabsbest[0]=errorabsf1[0];
              errorabsbest[1]=errorabsf1[1];
              radinvmodbest = *radinvmod;
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
        if(convreturn || convreturnok){
          // then really done
          failreturn=FAILRETURNNOFAIL; mathfailtype=0;
          break;
        }
        else{ // convreturn==0 && convreturnok==0
          
          // see if at least allowed error
          if(convreturnallow){

            // set as soft allowable failure
            failreturn=FAILRETURNNOTTOLERROR; mathfailtype=202;

            if(iter>trueimpmaxiter){// then reached maximum iterations
              if(debugfail>=2) dualfprintf(fail_file,"trueimpmaxiter=%d eomtype=%d MAXcheckconv=%d havebackup=%d failreturnallowable=%d: f1report=%g %g %g %g : f1=%g %g %g %g\n",trueimpmaxiter,*eomtype,checkconv,havebackup,failreturnallowable,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]]);

              if(showmessages && debugfail>=2) dualfprintf(fail_file,"iter>trueimpmaxiter=%d : iter exceeded in solve_implicit_lab().  But f1 was allowed error. checkconv=%d (if checkconv=0, could be issue!) : %g %g %g %g : %g %g %g %g : errorabs=%g %g : %g %g %g\n",trueimpmaxiter,checkconv,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]],errorabsf1[0],errorabsf1[1],fracdtuu0,fracuup,fracdtG);
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
              if(canbreak==2 && havebackup==0) dualfprintf(fail_file,"Aborted due to oscillatory error despite not having backup: iter=%d\n",iter);
              if(canbreak==3 && havebackup==0) dualfprintf(fail_file,"Aborted due to error not decreasing fast enough: iter=%d errorabsf1=%g %g\n",iter,errorabsf1[0],errorabsf1[1]);
              if(canbreak==4 && havebackup==0) dualfprintf(fail_file,"Aborted due to error not decreasing fast enough: iter=%d errorabsf1=%g %g\n",iter,errorabsf1[0],errorabsf1[1]);
              if(canbreak==5 && havebackup==0) dualfprintf(fail_file,"Aborted due to error not decreasing fast enough: iter=%d errorabsf1=%g %g\n",iter,errorabsf1[0],errorabsf1[1]);
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
        radinvmodbest = *radinvmod;
      }
    }


    // estimate effective work based on iterations done for each equation type
    FTYPE fndim=(FTYPE)NDIM;
    totaliters += (int)( (3.0/fndim)*(FTYPE)momiters + (1.0/fndim)*(FTYPE)energyiters + (fndim/fndim)*(FTYPE)fulliters);

  }// end loop over damping
  if(dampattempt==truenumdampattempts && truenumdampattempts>1){
    if(debugfail>=2) dualfprintf(fail_file,"Damping failed to avoid max iterations (but error might have dropped: %21.15g %21.15g): failreturn=%d dampattempt=%d eomtypelocal=%d *eomtype=%d ijk=%d %d %d\n",errorabsf1[0],errorabsf1[1],failreturn,dampattempt,eomtypelocal,*eomtype,ptrgeom->i,ptrgeom->j,ptrgeom->k);
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
      errorabsf1[0]=errorabsbest[0];
      errorabsf1[1]=errorabsbest[1];
      failreturn=failreturnbest;
      *radinvmod = radinvmodbest;
    }
  }




  if(ACTUALHARDFAILURE(failreturn)==0){ // deal with radinv issue but only didn't fail

    ///////
    //
    // But, in general, for QTYPMHD methods, might have use CAPTYPEFIX2 or CAPTYPEFIX1, but should use whichcap (probably CAPTYPEBASIC) to be conservative on value of Erf in general.  So while converged", could have avoided URAD0 error in total error for QTYPMHD methods.  So recover whichcap result before moving onto setting dUcomp.
    //
    ////////
    if(RADINVBAD(*radinvmod) && IMPMHDTYPEBASE(*baseitermethod)==1){
      int whichcall=FIMPLICITCALLTYPEFINALCHECK2; // KEY CHOICE IS THIS, which will use whichcap
      int goexplicitfake;
      int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive (u,v^i) type, i.e. v^i has no energy density term
      int fakef1iter=-1;
      FTYPE fakefracdtG=1.0;
      FTYPE f1fake[NPR],f1normfake[NPR],f1reportfake[NPR];
      FTYPE errorabsf1fake[NUMERRORTYPES];
      int convreturnfake=1;
      FTYPE fakeimpepsjac=1E-6;
      int fakefailreturnf=f_implicit(allowbaseitermethodswitch, iter,fakef1iter,failreturnallowableuse, whichcall,fakeimpepsjac,showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, baseitermethod, fracenergy, dissmeasure, radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, pp, pp, piin, uu, Uiin, uu0, uu, fakefracdtG*realdt, ptrgeom, q, f1fake, f1normfake, f1reportfake, &goexplicitfake, &errorabsf1fake[0], &errorabsf1fake[1], WHICHERROR, &convreturnfake, nummhdinvsreturn);
      // ONLY modifies uu and pp and q using true whichcap
      // *and* radinvmod, which contains whether radiative inversion modified solution.
    }

    ///////
    //
    // In general, uu might still not be uu=uu(pp) because might have used uualt in f_implicit()
    // But now that done with determining which solution we want to use with what error, now make consistent.
    // This needs to come *before* radsource or dUcomp is computed using uu, because uualt was only for error and not actual value that we want to fix-up to no longer require uualt in sense of making U change enough to be like primitive p.
    //
    ///////
    if(ALLOWUSEUUALT){
      get_state(pp, ptrgeom, q);
      primtoU(UNOTHING,pp,q,ptrgeom, uu, NULL);
    }


    if(0){ // no longer do this since do above full primtoU
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
    }

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

    
    //    DLOOPA(ii) radsource[iotherU[ii]] = -radsource[erefU[ii]]; // force energy conservation (NOTEMARK: Unsure if optimal since loop iterates and might break this for good reason -- e.g. no exactly valid corresponding MHD solution if iterating URAD)


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







  ////////////////
  //
  // Return error and iterations and fail mode
  //
  ////////////////

  // report error no matter whether got solution or not.
  errorabsreturn[0]=errorabsf1[0];
  errorabsreturn[1]=errorabsf1[1];

  // report iters not matter what the error.
  *itersreturn=totaliters;

  *nummhdstepsreturn = (int)(nstroke-nstrokeorig); // get number of mhd inversion steps




  // check if uncaught nan/inf
  int caughtnan=0;
  PLOOP(pliter,pl){
    if(!isfinite(pb[pl])) caughtnan++;
    if(!isfinite(uub[pl])) caughtnan++;
    if(!isfinite(dUcomp[RADSOURCE][pl])) caughtnan++;
  }
  if(!isfinite(errorabsreturn[0])) caughtnan++;  
  if(!isfinite(errorabsreturn[1])) caughtnan++;  

  if(caughtnan){
    // this doesn't seem to be hit even on Kraken
    if(debugfail>=2){
      dualfprintf(fail_file,"per mode implicit solver generated nan result and it wasn't caught\n");
      dualfprintf(fail_file,"per mode implicit solver: %d %d %d %d %d %d %d %d : %g %g %g : %d %d : %g %g : %d\n",allowbaseitermethodswitch, modprim, havebackup, didentropyalready, *eomtype, whichcap, itermode, *baseitermethod, trueimptryconv, trueimpokconv, trueimpallowconv, trueimpmaxiter, truenumdampattempts, fracenergy, dissmeasure, *radinvmod);
      PLOOP(pliter,pl) dualfprintf(fail_file,"0implicit solver: pl=%d Uiin=%21.15g dUother=%21.15g dU=%21.15g\n",pl,Uiin[pl],dUother[pl],dUcomp[RADSOURCE][pl]);
      PLOOP(pliter,pl) dualfprintf(fail_file,"1implicit solver: pl=%d pb=%21.15g piin=%21.15g Ufin=%21.15g dU=%21.15g uub=%21.15g\n",pl,pb[pl],piin[pl],Ufin[pl],dUcomp[RADSOURCE][pl],uub[pl]);
      int jjj;
      DLOOPA(jjj) dualfprintf(fail_file,"2implicit solver: jj=%d ucon=%21.15g ucov=%21.15g uradcon=%21.15g uradcov=%21.15g\n",jj,q->ucon[jj],q->ucov[jj],q->uradcon[jj],q->uradcov[jj]);
    }
    // reset solution as bad error, so conditions don't get confused by nan always giving false
    errorabsreturn[0]=1.0;
    errorabsreturn[1]=1.0;

    failreturn=FAILRETURNGENERAL;
  }
 


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


  if(PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1[WHICHERROR]>=trueimptryconvalt || PRODUCTION==0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){ // as in previous code

    PLOOP(pliter,pl) dualfprintf(fail_file,"ERRORCHECK: pl=%d f1=%21.15g f1norm=%21.15g f1report=%21.15g errorabsf1=%21.15g errorallabsf1=%21.15g\n",pl,f1[pl],f1norm[pl],f1report[pl],errorabsf1[0],errorabsf1[1]);
 


  // for seeing Erf<0 and small errors not tol errors.
  //  if(failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1[WHICHERROR]>=trueimptryconvalt || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){

  // for catching oscillators at small error but still >tol.
  //  if(PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0){


  //    if(REPORTERFNEG && failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT || RADINVBAD(*radinvmod) ) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1[WHICHERROR]>=trueimptryconvalt || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){
  //  if(REPORTERFNEG && failreturn!=FAILRETURNMODESWITCH && (pp[PRAD0]<10.0*ERADLIMIT) || PRODUCTION==0 && NOTACTUALFAILURE(failreturn)==0 && errorabsf1[WHICHERROR]>=trueimptryconvalt || PRODUCTION>0 && NOTBADFAILURE(failreturn)==0 && havebackup==0){
    //    if(NOTBADFAILURE(failreturn)==0){
    struct of_state qcheck; get_state(pp, ptrgeom, &qcheck);  primtoU(UNOTHING,pp,&qcheck,ptrgeom, uu, NULL);
    failnum++; mathematica_report_check(*radinvmod, mathfailtype, failnum, gotfirstnofail, eomtypelocal, itermode, *baseitermethod, errorabsf1, errorabsbestexternal, iter, totaliters, realdt, ptrgeom, ppfirst,pp,pb,piin,prtestUiin,prtestUU0,uu0,uu,Uiin,Ufin, CUf, q, dUother);
    int usedebugiter=debugiteratteempts[0];
    showdebuglist(usedebugiter,pppreholdlist,ppposholdlist,f1reportlist,f1list,errorabsf1list,errorallabsf1list,realiterlist,jaclist,fracdamplist,implicititerlist, implicitferrlist);


    // report how much MHD inversion used when failed
    dualfprintf(fail_file,"nummhdinvsreturn=%d nummhdstepsreturn=%d nstroke=%d\n",*nummhdinvsreturn,*nummhdstepsreturn,nstroke);
  }





  /// return
  return(failreturn);
  
}












// 0 : full
// 1: optimal
#define DEBUGMAXMODE 1

// DEBUGMAXITER stuff
static void showdebuglist(int debugiter, FTYPE (*pppreholdlist)[NPR],FTYPE (*ppposholdlist)[NPR],FTYPE (*f1reportlist)[NPR],FTYPE (*f1list)[NPR],FTYPE *errorabsf1list,FTYPE *errorallabsf1list, int *realiterlist, FTYPE (*jaclist)[NDIM][NDIM], FTYPE *fracdamplist, int *implicititerlist, int *implicitferrlist)
{

  if(DEBUGMAXITER==0) return;

  int listiter;
  if(DEBUGMAXMODE==0) dualfprintf(fail_file,"%3s : %3s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s %21s %21s : %21s %21s : %21s\n","li","ri","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","f1rep0","f1rep1","f1rep2","f1rep3","errorabs","errorallabs","umin");
  else if(DEBUGMAXMODE==1) dualfprintf(fail_file,"%3s : %3s : %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s : %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s : %21s %21s : %21s : %21s %21s %21s %21s : %21s %21s %21s : %2d %2d\n","li","ri","rho","ug","v1","v2","v3","Erf","vr1","vr2","vr3","rho","ug","f0","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","frep0","frep1","frep2","frep3","frep4","frep5","frep6","frep7","frep8","frep9","frep10","frep11","errorabs","errorallabs","umin","jac00","jac11","jac22","jac33","fracdtuu0","fracdtG","DAMPFACTOR","it","fe");

  for(listiter=0;listiter<=debugiter;listiter++){
    FTYPE umin=calc_PEQ_ufromTrho(TEMPMIN,pppreholdlist[listiter][RHO]);
    if(DEBUGMAXMODE==0){
      // full, but excessive
      dualfprintf(fail_file
                  ,"%3d : %3d : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g : %21.15g\n"
                  ,listiter,realiterlist[listiter]
                  ,pppreholdlist[listiter][RHO],pppreholdlist[listiter][UU],pppreholdlist[listiter][U1],pppreholdlist[listiter][U2],pppreholdlist[listiter][U3],pppreholdlist[listiter][PRAD0],pppreholdlist[listiter][PRAD1],pppreholdlist[listiter][PRAD2],pppreholdlist[listiter][PRAD3]
                  ,ppposholdlist[listiter][RHO],ppposholdlist[listiter][UU],ppposholdlist[listiter][U1],ppposholdlist[listiter][U2],ppposholdlist[listiter][U3],ppposholdlist[listiter][PRAD0],ppposholdlist[listiter][PRAD1],ppposholdlist[listiter][PRAD2],ppposholdlist[listiter][PRAD3]
                  ,f1reportlist[listiter][0],f1reportlist[listiter][1],f1reportlist[listiter][2],f1reportlist[listiter][3]
                  ,errorabsf1list[listiter]
                  ,errorallabsf1list[listiter]
                  ,umin
                  );
    }
    else if(DEBUGMAXMODE==1){
      // optimal
      dualfprintf(fail_file
                  ,"%3d : %3d : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g : %21.15g  %21.15g : %21.15g : %21.15g %21.15g %21.15g %21.15g : %21.15g %21.15g %21.15g : %d %d\n"
                  ,listiter,realiterlist[listiter]
                  ,pppreholdlist[listiter][RHO],pppreholdlist[listiter][UU],pppreholdlist[listiter][U1],pppreholdlist[listiter][U2],pppreholdlist[listiter][U3],pppreholdlist[listiter][PRAD0],pppreholdlist[listiter][PRAD1],pppreholdlist[listiter][PRAD2],pppreholdlist[listiter][PRAD3]
                  ,ppposholdlist[listiter][RHO],ppposholdlist[listiter][UU]
                  ,f1list[listiter][0],f1list[listiter][1],f1list[listiter][2],f1list[listiter][3],f1list[listiter][4],f1list[listiter][5],f1list[listiter][6],f1list[listiter][7],f1list[listiter][8],f1list[listiter][9],f1list[listiter][10],f1list[listiter][11]
                  ,f1reportlist[listiter][0],f1reportlist[listiter][1],f1reportlist[listiter][2],f1reportlist[listiter][3],f1reportlist[listiter][4],f1reportlist[listiter][5],f1reportlist[listiter][6],f1reportlist[listiter][7],f1reportlist[listiter][8],f1reportlist[listiter][9],f1reportlist[listiter][10],f1reportlist[listiter][11]
                  ,errorabsf1list[listiter]
                  ,errorallabsf1list[listiter]
                  ,umin
                  ,jaclist[listiter][0][0]
                  ,jaclist[listiter][1][1]
                  ,jaclist[listiter][2][2]
                  ,jaclist[listiter][3][3]
                  ,fracdamplist[0] // fracdtuu0
                  ,fracdamplist[1] // fracdtG
                  ,fracdamplist[2] // DAMPFACTOR
                  ,implicititerlist[listiter], implicitferrlist[listiter]
                );
    }
  }// end listiter loop

}






                             
int get_rameshsolution_wrapper(int whichcall, int eomtype, FTYPE *errorabs, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *dUother, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, FTYPE (*dUcompeng)[NPR], FTYPE (*dUcompent)[NPR], struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent)
{
  // BEGIN get ramesh solution
  int radinvmod=0; // fake
  int failreturnentropy=0; // fake
  //  int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITERLONG,itersent=IMPMAXITERLONG;
  int failnum=0,gotfirstnofail=0;
  int itermode=ITERMODENORMAL,iters=IMPMAXITERLONG;
  int baseitermethod=QTYPMHD;
  FTYPE errorabsbestexternal[2];
  set_array(errorabsbestexternal,2,MPI_FTYPE,BIG);
  // itersentropy is last, but might feed in best
  FTYPE realdt = compute_dt(CUf,dt);
  // get uu0
  FTYPE uu0[NPR];
  FTYPE fracdtuu0=1.0;
  int pliter,pl;
  PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);
  FTYPE uu[NPR]; // not used except when doing diagnostics in ramesh code
  //  FTYPE uueng[NPR],uuent[NPR]; // filled with answer if successful
  //
  get_rameshsolution(whichcall, radinvmod, failreturnentropy, failnum,  gotfirstnofail,  eomtype,  itermode, baseitermethod, errorabs, errorabs, iters, iters, realdt, ptrgeom, pp, pp, piin, uu0, uu, Uiin, Ufin, CUf, q, ppeng, ppent, uueng, uuent, qeng, qent, failtypeeng, errorabseng, iterseng, radinvmodeng, failtypeent, errorabsent, itersent, radinvmodent);
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

int get_rameshsolution(int whichcallramesh, int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, int baseitermethod, FTYPE *errorabs, FTYPE *errorabsbestexternal, int iters, int totaliters, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *ppeng, FTYPE *ppent, FTYPE *uueng, FTYPE *uuent, struct of_state *qeng, struct of_state *qent, int *failtypeeng, FTYPE *errorabseng, int *iterseng, int *radinvmodeng, int *failtypeent, FTYPE *errorabsent, int *itersent, int *radinvmodent)
{
  //  int failtype=1;
  //  long long int failnum=0;
  //  int gotfirstnofail=0;
  //  int eomtypelocal=EOMGRMHD;
  //  int itermode=ITERMODENORMAL;
  //  FTYPE errorabs=BIG;
  //  FTYPE iters=IMPMAXITERLONG;
  //  FTYPE errorabsbestexternal=BIG;
  //  int totaliters=IMPMAXITERLONG;
  FTYPE *ppfirst=pp,*prtestUiin=piin,*prtestUU0=piin;
  struct of_state *qpp=q;
  struct of_state *qpb=q;
  struct of_state *qpiin=q;
  int jj,kk;
  int pliter,pl;



  // 11 vars, failcode, error, iterations
  // MUST BE SAME AS IN test.f code
#define NUMRESULTS 15
  doublereal resultseng[NUMRESULTS]={0},resultsent[NUMRESULTS]={0};


  if(WHICHVELRAMESH==VEL4 && ptrgeom->gcov[GIND(TT,TT)]>=0.0 || WHICHEOS!=IDEALGAS){
    // KORALTODO SUPERGODMARK: Also only can use ramesh if lambda cooling only includes kappaff inverse and only normal ff and es opacities.
    if(WHICHVELRAMESH==VEL4 && ptrgeom->gcov[GIND(TT,TT)]>=0.0) dualfprintf(fail_file,"Wanted to call ramesh, but inside ergosphere\n");
    *failtypeeng = 1;
    *radinvmodent = 0;
    errorabseng[0] = errorabseng[1] = BIG;
    *iterseng = 0;
    *failtypeent = 1;
    *radinvmodent = 0;
    errorabsent[0] = errorabsent[1] = BIG;
    *itersent = 0;
    for(pl=0;pl<NUMRESULTS;pl++) resultseng[pl]=resultsent[pl]=-1.0;
    PLOOP(pliter,pl) ppeng[pl]=ppent[pl]=-1.0;
  }
  else{

    /////////////////////////
    //
    // Call Ramesh's solver
    //
    /////////////////////////


    // below things must be same order as in test.f
    // MUST BE SAME AS IN test.f code
#define NUMARGS (211+11)

    // call fortran code
    doublereal args[NUMARGS]={0};
    int na;

    na=-1;
    na++; args[na]=(doublereal)NUMARGS;
    na++; args[na]=(doublereal)NUMRESULTS;
    na++; args[na]=(doublereal)WHICHVELRAMESH;
    na++; args[na]=(doublereal)failtype;
    na++; args[na]=(doublereal)myid;
    na++; args[na]=(doublereal)failnum;
    na++; args[na]=(doublereal)gotfirstnofail;
    na++; args[na]=(doublereal)eomtypelocal;
    na++; args[na]=(doublereal)itermode;
    na++; args[na]=errorabs[WHICHERROR];
    na++; args[na]=errorabsbestexternal[WHICHERROR];
    na++; args[na]=(doublereal)iters;
    na++; args[na]=(doublereal)totaliters;
    na++; args[na]=realdt;
    na++; args[na]=(doublereal)nstep;
    na++; args[na]=(doublereal)steppart;
    na++; args[na]=gam;
    na++; args[na]=GAMMAMAXRAD;
    na++; args[na]=ERADLIMIT;
    FTYPE trueimptryconv=IMPTRYCONV;
    na++; args[na]=IMPTRYCONVABS;
    na++; args[na]=IMPALLOWCONVCONSTABS;
    na++; args[na]=ARAD_CODE;
    // so as really code uses:
    na++; args[na]=calc_kappaes_user(1.0,1.0,0,0,0);
    na++; args[na]=calc_kappa_user(1.0,1.0,0,0,0);
    na++; args[na]=0.0;
    DLOOP(jj,kk){ na++; args[na]=ptrgeom->gcon[GIND(jj,kk)];}
    DLOOP(jj,kk){ na++; args[na]=ptrgeom->gcov[GIND(jj,kk)];}
    PLOOP(pliter,pl){ na++; args[na]=pp[pl]; na++; args[na]=ppfirst[pl]; na++; args[na]=pb[pl]; na++; args[na]=piin[pl]; na++; args[na]=prtestUiin[pl]; na++; args[na]=prtestUU0[pl]; na++; args[na]=uu0[pl]; na++; args[na]=uu[pl]; na++; args[na]=Uiin[pl]; }
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpp->uradcon[jj]; na++; args[na]=qpp->uradcov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpp->ucon[jj]; na++; args[na]=qpp->ucov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }

    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpb->uradcon[jj]; na++; args[na]=qpb->uradcov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpb->ucon[jj]; na++; args[na]=qpb->ucov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }

    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpiin->uradcon[jj]; na++; args[na]=qpiin->uradcov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj){ na++; args[na]=qpiin->ucon[jj]; na++; args[na]=qpiin->ucov[jj]; }
    else DLOOPA(jj){ na++; args[na]=0.0; na++; args[na]=0.0; }

    if(na!=NUMARGS-1){
      dualfprintf(fail_file,"Wrong number of args=%d\n",na);
      myexit(304583453);
    }

    // calling f2c generated code, which assumes single array as input of 211 items.
    rameshsolver_(args,resultseng,resultsent);

    // process results from rameshsolver()

    *failtypeeng = (int)resultseng[11];
    if(*failtypeeng==0){
      PLOOPBONLY(pl) ppeng[pl]=pp[pl]; // ramesh solver doesn't pass this back
      FTYPE uconeng[NDIM],uradconeng[NDIM];
      // return solution in harm format
      ppeng[RHO] = resultseng[0];
      ppeng[UU] = ppeng[ENTROPY] = resultseng[1];
      if(WHICHVELRAMESH==VEL4){
        uconeng[0] = resultseng[2];
        uconeng[1] = resultseng[3];
        uconeng[2] = resultseng[4];
        uconeng[3] = resultseng[5];
        uconrel(uconeng,&ppeng[UU],ptrgeom); // get \tilde{u}^i
      }
      else{
        SLOOPA(jj) ppeng[UU+jj] = resultseng[2+jj];
      }
      ppeng[URAD0] = resultseng[6];
      if(WHICHVELRAMESH==VEL4){
        uradconeng[0] = resultseng[7];
        uradconeng[1] = resultseng[8];
        uradconeng[2] = resultseng[9];
        uradconeng[3] = resultseng[10];
        uconrel(uradconeng,&ppeng[URAD0],ptrgeom); // get \tilde{u}^i_{\rm rad}
      }
      else{
        SLOOPA(jj) ppeng[URAD0+jj] = resultseng[7+jj];
      }
      // get full state (thermo and other stuff needed)
      get_state(ppeng,ptrgeom,qeng);
      primtoU(UNOTHING,ppeng, qeng, ptrgeom, uueng, NULL);
      extern int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr);
      invert_scalars(ptrgeom, uueng,ppeng);
    }
    errorabseng[0] = resultseng[12];
    *iterseng = (int)resultseng[13];
    *radinvmodeng = (int)resultseng[14];

    // override ramesh error in case meaningless, if there is a failure.
    if(*failtypeeng) errorabseng[0] = BIG;

    // return solution in harm format
    *failtypeent = (int)resultsent[11];
    if(*failtypeent==0){
      PLOOPBONLY(pl) ppent[pl]=pp[pl]; // ramesh solver doesn't pass this back
      FTYPE uconent[NDIM],uradconent[NDIM];
      ppent[RHO] = resultsent[0];
      ppent[UU] = ppent[ENTROPY] = resultsent[1];
      if(WHICHVELRAMESH==VEL4){
        uconent[0] = resultsent[2];
        uconent[1] = resultsent[3];
        uconent[2] = resultsent[4];
        uconent[3] = resultsent[5];
        uconrel(uconent,&ppent[UU],ptrgeom); // get \tilde{u}^i
      }
      else{
        SLOOPA(jj) ppent[UU+jj] = resultsent[2+jj];
      }
      ppent[URAD0] = resultsent[6];
      if(WHICHVELRAMESH==VEL4){
        uradconent[0] = resultsent[7];
        uradconent[1] = resultsent[8];
        uradconent[2] = resultsent[9];
        uradconent[3] = resultsent[10];
        uconrel(uradconent,&ppent[URAD0],ptrgeom); // get \tilde{u}^i_{\rm rad}
      }
      else{
        SLOOPA(jj) ppent[URAD0+jj] = resultsent[7+jj];
      }
      // get full state (thermo and other stuff needed)
      get_state(ppent,ptrgeom,qent);
      primtoU(UNOTHING,ppent, qent, ptrgeom, uuent, NULL);
      extern int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr);
      invert_scalars(ptrgeom, uuent,ppent);
    }
    errorabsent[0] = resultsent[12];
    *itersent = (int)resultsent[13];
    *radinvmodent = (int)resultsent[14];

    // override ramesh error in case meaningless, if there is a failure.
    if(*failtypeent) errorabsent[0] = BIG;
  }



  //////////
  //
  // Double check ramesh error with f_error_check() or f_implicit().  Also uses error as consistent with harm -- i.e. total error if WHICHERROR==1
  //
  //////////
  int failreturnferr;
  int allowbaseitermethodswitch=0;
  int fakeiter=1;
  int fakef1iter=-1;
  int failreturnallowableuse=0;
  int whichcall=FIMPLICITCALLTYPEFINALCHECK2;
  FTYPE impepsjac=0;
  int showmessages=0;
  int showmessagesheavy=0;
  int allowlocalfailurefixandnoreport=0;
  //int &eomtypelocal;
  int whichcap=CAPTYPEBASIC;
  //  int itermode;
  //  int baseitermethod;
  FTYPE fracenergy=1.0;
  FTYPE dissmeasure=-1.0;
  radinvmod=*radinvmodeng; // default
  FTYPE trueimptryconv=IMPTRYCONV;
  FTYPE trueimptryconvabs=IMPTRYCONVABS;
  FTYPE trueimpallowconvabs=IMPALLOWCONVCONSTABS;
  int trueimpmaxiter=IMPMAXITERLONG;
  //  FTYPE realdt;
  int dimtypef=DIMTYPEFCONS;
  // some geometry stuff to store pre-step instead of for each step.
  int jjdim;
  FTYPE dimfactU[NPR];
  PLOOP(pliter,pl) dimfactU[pl]=1.0; // default
  DLOOPA(jjdim) dimfactU[UU+jjdim]=dimfactU[URAD0+jjdim]=sqrt(fabs(ptrgeom->gcon[GIND(jjdim,jjdim)]));
  SLOOPA(jjdim) dimfactU[B1+jjdim-1] = 1.0/dimfactU[U1+jjdim-1];
  FTYPE ppppp[NPR];
  //  FTYPE pp[NPR];
  //  FTYPE piin[NPR];
  FTYPE uuppp[NPR];
  //FTYPE Uiin[NPR];
  //  FTYPE uu0[NPR];
  //  FTYPE uu[NPR];
  FTYPE fracdtG=1.0;
  //FTYPE ptrgeom;
  //FTYPE q;
  FTYPE f1[NPR];
  FTYPE f1norm[NPR];
  FTYPE f1report[NPR];
  int goexplicit=0;
  FTYPE errorabsf1[NUMERRORTYPES];
  // WHICHERROR;
  int convreturnf1;


  int nummhdinvsreturn=0;
  if(*failtypeeng==0){
    ////////
    //
    // test ppeng
    //
    ////////
    PLOOP(pliter,pl){
      ppppp[pl] = pp[pl] = ppeng[pl];
      uuppp[pl] = uu[pl] = uueng[pl];
    }
 
    failreturnferr=f_implicit(allowbaseitermethodswitch, fakeiter, fakef1iter, failreturnallowableuse, whichcall, impepsjac, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, &baseitermethod, fracenergy, dissmeasure, &radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, ppppp, pp, piin, uuppp, Uiin, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, f1report, &goexplicit, &errorabsf1[0], &errorabsf1[1], WHICHERROR, &convreturnf1, &nummhdinvsreturn); // modifies uu and pp, f1poret, goexplicit, errorabsf1[0,1], convreturnf1
    // translate result
    //  *qeng=*q; // already done
    *failtypeeng=failreturnferr; // trust my error report
    errorabseng[0]=errorabsf1[0]; // trust my error measure
    errorabseng[1]=errorabsf1[1]; // trust my error measure
    *radinvmodeng=radinvmod; // trust my rad inv check
  }

  if(*failtypeent==0){
    ////////
    //
    // test ppent
    //
    ////////
    PLOOP(pliter,pl){
      ppppp[pl] = pp[pl] = ppent[pl];
      uuppp[pl] = uu[pl] = uuent[pl];
    }
 
    failreturnferr=f_implicit(allowbaseitermethodswitch, fakeiter, fakef1iter, failreturnallowableuse, whichcall, impepsjac, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, whichcap, itermode, &baseitermethod, fracenergy, dissmeasure, &radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, ppppp, pp, piin, uuppp, Uiin, uu0, uu, fracdtG*realdt, ptrgeom, q, f1, f1norm, f1report, &goexplicit, &errorabsf1[0], &errorabsf1[1], WHICHERROR, &convreturnf1, &nummhdinvsreturn); // modifies uu and pp, f1poret, goexplicit, errorabsf1[0,1], convreturnf1
    // translate result
    //  *qent=*q; // already done
    *failtypeent=failreturnferr; // trust my error report
    errorabsent[0]=errorabsf1[0]; // trust my error measure
    errorabsent[1]=errorabsf1[1]; // trust my error measure
    *radinvmodent=radinvmod; // trust my rad inv check
  }

  


  // DEBUG and only when real call (for now) and only show if one is not a bad solution
  if(debugfail>=2 && whichcallramesh>0 && (*failtypeent==0 || *failtypeeng==0)){ // NORMAL
  //  if(debugfail>=2 && whichcallramesh==0){ // DEBUGGING Erf~0 solutions
  //  if(PRODUCTION==0 && debugfail>=2 && (*failtypeeng==0 || *failtypeent==0)){ // DEBUGGING Erf~0 solutions

    FTYPE resultsjon[NUMRESULTS];
    resultsjon[0] = pp[RHO];
    resultsjon[1] = pp[UU];
    if(WHICHVELRAMESH==VEL4){
      resultsjon[2] = q->ucon[0];
      resultsjon[3] = q->ucon[1];
      resultsjon[4] = q->ucon[2];
      resultsjon[5] = q->ucon[3];
    }
    else{
      resultsjon[2] = 0.0;
      SLOOPA(jj) resultsjon[2+jj] = pp[UU+jj];
    }
    resultsjon[6] = pp[URAD0];
    if(WHICHVELRAMESH==VEL4){
      resultsjon[7] = q->uradcon[0];
      resultsjon[8] = q->uradcon[1];
      resultsjon[9] = q->uradcon[2];
      resultsjon[10] = q->uradcon[3];
    }
    else{
      resultsjon[7]=0.0;
      SLOOPA(jj) resultsjon[7+jj] = pp[URAD0+jj];
    }
    resultsjon[11] = (FTYPE)failtype;
    resultsjon[12] = errorabs[WHICHERROR];
    resultsjon[13] = (FTYPE)iters;
    resultsjon[14] = radinvmod;

    for(pl=0;pl<NUMRESULTS;pl++) dualfprintf(fail_file,"RAMESH1: pl=%d | resultsjon=%21.15g | resultseng=%21.15g resultsent=%21.15g\n",pl,resultsjon[pl],resultseng[pl],resultsent[pl]);
    SLOOPA(jj) dualfprintf(fail_file,"RAMESH2: jj=%d utildegasjon=%21.15g utildegaseng=%21.15g utildegasent=%21.15g\n",jj,pp[U1+jj-1],ppeng[U1+jj-1],ppent[U1+jj-1]);
    SLOOPA(jj) dualfprintf(fail_file,"RAMESH3: jj=%d utilderadjon=%21.15g utilderadeng=%21.15g utilderadent=%21.15g\n",jj,pp[URAD1+jj-1],ppeng[URAD1+jj-1],ppent[URAD1+jj-1]);
    PLOOP(pliter,pl) dualfprintf(fail_file,"RAMESH3.5: pl=%d ppeng=%21.15g uueng=%21.15g ppent=%21.15g uuent=%21.15g\n",pl,ppeng[pl],uueng[pl],ppent[pl],ppeng[pl]);

    // get radiative inversion stuff for Ramesh
    FTYPE Er,Utildesq,Utildecon[NDIM];

    compute_ZAMORAD(uu0, ptrgeom, &Er, &Utildesq, Utildecon);
    dualfprintf(fail_file,"RAMESH4: E_r=%21.15g Utildesq=%21.15g Utildecon0=%21.15g Utildecon1=%21.15g Utildecon2=%21.15g Utildecon3=%21.15g\n",Er,Utildesq,Utildecon[0],Utildecon[1],Utildecon[2],Utildecon[3]);

    if(*failtypeeng==0){
      compute_ZAMORAD(uueng, ptrgeom, &Er, &Utildesq, Utildecon);
      dualfprintf(fail_file,"RAMESH5: E_r=%21.15g Utildesq=%21.15g Utildecon0=%21.15g Utildecon1=%21.15g Utildecon2=%21.15g Utildecon3=%21.15g\n",Er,Utildesq,Utildecon[0],Utildecon[1],Utildecon[2],Utildecon[3]);
    }

    if(*failtypeent==0){
      compute_ZAMORAD(uuent, ptrgeom, &Er, &Utildesq, Utildecon);
      dualfprintf(fail_file,"RAMESH6: E_r=%21.15g Utildesq=%21.15g Utildecon0=%21.15g Utildecon1=%21.15g Utildecon2=%21.15g Utildecon3=%21.15g\n",Er,Utildesq,Utildecon[0],Utildecon[1],Utildecon[2],Utildecon[3]);
    }

  }// end if DEBUG

  return(0);
}



int mathematica_report_check(int radinvmod, int failtype, long long int failnum, int gotfirstnofail, int eomtypelocal, int itermode, int baseitermethod, FTYPE *errorabs, FTYPE *errorabsbestexternal, int iters, int totaliters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *prtestUiin, FTYPE *prtestUU0, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother)
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
    if(!isfinite(errorabs[0])) errorabs[0]=BIG;
    if(!isfinite(errorabs[1])) errorabs[1]=BIG;
    if(!isfinite(errorabsbestexternal[0])) errorabsbestexternal[0]=BIG;
    if(!isfinite(errorabsbestexternal[1])) errorabsbestexternal[1]=BIG;

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
    dualfprintf(fail_file,"%d %d %d ",NUMARGS,NUMRESULTS,WHICHVELRAMESH); // 3
    dualfprintf(fail_file,"%d %d %lld %d %d %d %21.15g %21.15g %d %d %21.15g %lld %d %21.15g ",failtype,myid,failnum,gotfirstnofail,eomtypelocal,itermode,errorabs[WHICHERROR],errorabsbestexternal[WHICHERROR],iters,totaliters,realdt,nstep,steppart,gam); // 14
    FTYPE trueimptryconv=IMPTRYCONV;
    dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g ",GAMMAMAXRAD,ERADLIMIT,IMPTRYCONVABS,IMPALLOWCONVCONSTABS); // 4
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
    FTYPE ppeng[NPR]={0},ppent[NPR]={0},errorabseng[NUMERRORTYPES],errorabsent[NUMERRORTYPES];
    set_array(errorabseng,NUMERRORTYPES,MPI_FTYPE,BIG);
    set_array(errorabsent,NUMERRORTYPES,MPI_FTYPE,BIG);

    FTYPE uueng[NPR]={0},uuent[NPR]={0};
    int failtypeeng=1,failtypeent=1,iterseng=IMPMAXITERLONG,itersent=IMPMAXITERLONG, radinvmodeng=UTOPRIMRADFAILBAD1, radinvmodent=UTOPRIMRADFAILBAD1;
    struct of_state qeng=*q,qent=*q;
    get_rameshsolution(whichcall, radinvmod, failtype, failnum,  gotfirstnofail,  eomtypelocal,  itermode, baseitermethod, errorabs, errorabsbestexternal,  iters,  totaliters,realdt, ptrgeom, pp, pb, piin, uu0, uu, Uiin, Ufin, CUf, q, ppeng, ppent, uueng, uuent, &qeng, &qent, &failtypeeng, errorabseng, &iterseng, &radinvmodeng, &failtypeent, errorabsent, &itersent, &radinvmodent);



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


// use f and check the error
// note that eomtype is just integer, not pointer as often case when might want to change eomtype.
static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE convabs, FTYPE realdt, int dimtypef, int eomtype, int radinvmod, int itermode, int baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, FTYPE *dimfactU, FTYPE *pp, FTYPE *piin, FTYPE *fin, FTYPE *finnorm, FTYPE *finreport, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom, FTYPE *errorabs, FTYPE *errorallabs, int whicherror)
{
  int ii,jj;
  int pliter,pl;

  // setup method and signs
  int numdims,startjac,endjac,implicititer,implicitferr,BEGINMOMSTEPS,ENDMOMSTEPS,BEGINENERGYSTEPS,ENDENERGYSTEPS,BEGINFULLSTEPS,ENDFULLSTEPS,BEGINNORMALSTEPS,irefU[NDIM],iotherU[NDIM],erefU[NDIM],eotherU[NDIM],signgd2,signgd4,signgd6,signgd7;
  define_method(iter, &eomtype, itermode, baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
  get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);


  int passedconv;

  // default
  passedconv=0;

  // get error
  // NOTE: use of gcov[ii,ii] so comparable dimensionally to fin[erefU[ii]] and finnorm[erefU[ii]] that are like R^t_\nu and so need sqrt(gcon[nu,nu]) multiplied on them.  This ensures error is non-dimensional (or, really only ^t dimensional)
  FTYPE dimfactferr[NPR];
  PLOOP(pliter,pl) dimfactferr[pl]=dimfactU[pl]; // default
  if(dimtypef==DIMTYPEFCONS){
    // assume fin and finnorm are conservative or source (R^t_\nu form)
    DLOOPA(ii) dimfactferr[erefU[ii]]=dimfactU[erefU[ii]];
    DLOOPA(ii) dimfactferr[eotherU[ii]]=dimfactU[eotherU[ii]];
  }
  else{
    // assume fin and finnorm are primitve (u,vel^i form)
    DLOOPA(ii) dimfactferr[erefU[ii]]=1.0/dimfactU[erefU[ii]];
    DLOOPA(ii) dimfactferr[eotherU[ii]]=1.0/dimfactU[eotherU[ii]];
  }

  // replace finnorm -> finnormnew that's already non-dimensionalized
  FTYPE finnormnew[NPR];
  // Tds\rho_0 u^t \propto \rho_0 (v/c)^2 or higher  \propto \gamma in ultrarel limit
  // T^t_t + \rho_0 u^t \propto \rho_0 (v/c)^2 or higher  \propto \rho_0\gamma^2 in ultrarel limit
  ii=TT;
  FTYPE fnormtime = fabs(finnorm[erefU[ii]]*dimfactferr[erefU[ii]]);
  FTYPE fnormtimeother = fabs(finnorm[eotherU[ii]]*dimfactferr[eotherU[ii]]);
  // T^t_i \propto \rho_0 (v/c)^1 or higher  \propto \gamma^2 \rho_0 (v/c) in ultrarel limit
  // get spatial contributions as total term so not dominated by small errors in some dimensions
  FTYPE fnormspace=0.0; SLOOPA(ii) fnormspace += fabs(finnorm[erefU[ii]]*dimfactferr[erefU[ii]]);
  FTYPE fnormspaceother=0.0; SLOOPA(ii) fnormspaceother += fabs(finnorm[eotherU[ii]]*dimfactferr[eotherU[ii]]);

  // need to account for errors when momentum~0
#if(0)
#if(0)
  // only kinda applicable for QTYPMHD
  FTYPE rhoref=MAX(pp[RHO],piin[RHO]);
  FTYPE fnormspace2 = rhoref*pow(fabs(fnormtime)/rhoref,0.5); // energy term in correct scale with \rho_0(v/c)^1 to get reference in case v\sim 0 NOTE: This assumes rho v term dominates
  FTYPE fnormtime2=rhoref*pow(fabs(fnormspace/rhoref),2.0); // momentum term in correct scale with \rho_0(v/c)^2 in case E\sim 0.  NOTE: assumes rho v term dominates.
#else
  FTYPE fakevel = MIN(1.0,fnormspace/(SMALL+fnormtime));
  FTYPE fnormspace2 = fabs(fnormtime)*fakevel;
  //  FTYPE fnormtime2 = fabs(fnormspace)/fakevel;
  FTYPE fnormtime2 = fnormtime; // no problem with this since now include absolute sum version
#endif
  //  dualfprintf(fail_file,"fnormspacetime: %g %g %g %g: uu0=%g fakevel=%g\n",fnormspace,fnormspace2,fnormtime,fnormtime2,fabs(uu0[erefU[0]]*dimfactU[erefU[0]]),fakevel);
#else
  // These suggest that velocity can be no better than NUMEPSILON*c=NUMEPSILON, but non-rel velocities can be much smaller.  But with rad inversion, v^i comes from mixed-up R^t_\mu -- catastrophic cancellation issue.  What about MHD?
  FTYPE fnormspace2 = fnormtime; // maybe ok SUPERGODMARK -- problem for very small non-rel velocities. FUCK
  FTYPE fnormtime2 = fnormtime; // no problem with this since now include absolute sum version
  FTYPE fnormspaceother2 = fnormtimeother; // maybe ok SUPERGODMARK -- problem for very small non-rel velocities. FUCK
  FTYPE fnormtimeother2 = fnormtimeother; // no problem with this since now include absolute sum version
#endif

  // now assign
  PLOOP(pliter,pl) finnormnew[pl] = finnorm[pl]; // default
  ii=TT; finnormnew[erefU[ii]] = fnormtime + fnormtime2;
  SLOOPA(ii) finnormnew[erefU[ii]] = fnormspace + fnormspace2;
  ii=TT; finnormnew[eotherU[ii]] = fnormtimeother + fnormtimeother2;
  SLOOPA(ii) finnormnew[eotherU[ii]] = fnormspaceother + fnormspaceother2;


  // get non-dimensionalized fin
  FTYPE finnew[NPR];
  PLOOP(pliter,pl) finnew[pl]=fin[pl]*dimfactferr[pl];

  // get relative errors (keep sign)
  //  JACLOOPALT(ii,startjac,endjac)
  PLOOP(pliter,pl) finreport[pl]=finnew[pl]/fabs(IMPMINABSERROR+fabs(finnormnew[pl]));

  // get absolute error over all (baseitermethod)-iterated terms
  // NOTE: the SUBJACJ methods directly use freport[] as needed, not errorabs or passedconv
  *errorabs=0.0;
  JACLOOPFULLERROR(itermode,jj,startjac,endjac)      *errorabs     += fabs(finreport[erefU[jj]]); // always full error.

  // completely full relevant error
  *errorallabs=0.0; JACLOOPSUPERFULL(pliter,pl,eomtype,baseitermethod,radinvmod) *errorallabs += fabs(finreport[pl]);

  //////////////
  //
  // see if passed convergence test criteria
  //
  ////////////////
  if(0){ // old way of measuring if passed error test
    if(whicherror==0){
      passedconv=1; JACLOOPFULLERROR(itermode,ii,startjac,endjac) passedconv *= fabs(finreport[erefU[ii]])<conv;
    }
    else if(whicherror==1){
      passedconv=1; JACLOOPSUPERFULL(pliter,pl,eomtype,baseitermethod,radinvmod) passedconv *= fabs(finreport[pl])<conv;
    }
  }
  else{ // new way, so easier to test outside with single number
    if(whicherror==0){
      passedconv=(*errorabs<convabs);
    }
    else if(whicherror==1){
      passedconv=(*errorallabs<convabs);
    }
  }
  

  // report if passed convergence test
  if(passedconv){
    if(PRODUCTION==0&&showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g realdt=%g i=%d iter=%d DONE1 for conv=%g : finreport=%g %g %g %g\n",nstep,steppart,dt,realdt,ptrgeom->i,iter,conv,finreport[erefU[0]],finreport[erefU[1]],finreport[erefU[2]],finreport[erefU[3]]);
    return(1);
  }
  else{
    // report if didn't pass
    if(PRODUCTION==0&&showmessagesheavy){
      dualfprintf(fail_file,"POSTFIN (conv=%21.15g): uu: %21.15g %21.15g %21.15g %21.15g : uu0=%21.15g %21.15g %21.15g %21.15g\n",conv,uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],uu0[irefU[0]],uu0[irefU[1]],uu0[irefU[2]],uu0[irefU[3]]);
      PLOOP(pliter,pl) dualfprintf(fail_file,"iii=%d fin=%21.15g finnorm=%21.15g\n",ii,fin[pl],finnorm[pl]);
      dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,finreport[erefU[0]],finreport[erefU[1]],finreport[erefU[2]],finreport[erefU[3]]);
    }
    return(0);
  }

  return(0);
}


#define JDIFFONESIDED 0
#define JDIFFCENTERED 1



// calculating approximate Jacobian: dUresid(dUrad,G(Urad))/dUrad = dy(x)/dx
// then compute inverse Jacobian
static int get_implicit_iJ(int allowbaseitermethodswitch, int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtypelocal, int whichcap, int itermode, int *baseitermethod, FTYPE fracenergy, FTYPE dissmeasure, FTYPE impepsjac, FTYPE trueimptryconv, FTYPE trueimptryconvabs, FTYPE trueimpallowconvabs, int trueimpmaxiter, int iter, FTYPE errorabs, FTYPE errorallabs, int whicherror, int dimtypef, FTYPE *dimfactU, FTYPE *Uiin, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *piin, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, struct of_state *q, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NPR], int *nummhdinvsreturn)
{
  int ii,jj;
  struct of_state qjac=*q; // not required as input, but set as output.

  int eomtypelocallocal=*eomtypelocal; // default

  // setup method and signs
  // ok to do this heree, because eomtypelocallocal is re-defaulted to *eomtypelocal for every call to the error function.
  int numdims,startjac,endjac,implicititer,implicitferr,BEGINMOMSTEPS,ENDMOMSTEPS,BEGINENERGYSTEPS,ENDENERGYSTEPS,BEGINFULLSTEPS,ENDFULLSTEPS,BEGINNORMALSTEPS,irefU[NDIM],iotherU[NDIM],erefU[NDIM],eotherU[NDIM],signgd2,signgd4,signgd6,signgd7;
  define_method(iter,&eomtypelocallocal, itermode, *baseitermethod, fracenergy, dissmeasure, &implicititer, &implicitferr, &BEGINMOMSTEPS, &ENDMOMSTEPS, &BEGINENERGYSTEPS, &ENDENERGYSTEPS, &BEGINFULLSTEPS, &ENDFULLSTEPS, &BEGINNORMALSTEPS);
  get_refUs(&numdims, &startjac, &endjac, &implicititer, &implicitferr, irefU, iotherU, erefU, eotherU, &signgd2, &signgd4, &signgd6, &signgd7);

  int JDIFFTYPE;
  if(IMPPTYPE(implicititer)){
    // with implicititer==QTYPMHD, no longer expensive so can do JDIFFCENTERED
    // choose:
    //    if(itermode==ITERMODENORMAL) JDIFFTYPE=JDIFFCENTERED; // and only helps rarely and makes 2X slower.
    //    else if(itermode==ITERMODESTAGES) JDIFFTYPE=JDIFFONESIDED; // seems accurate enough
    if(itermode==ITERMODESTAGES) JDIFFTYPE=JDIFFCENTERED; // and only helps rarely and makes 2X slower.
    else if(itermode==ITERMODENORMAL) JDIFFTYPE=JDIFFONESIDED; // seems accurate enough
    //    JDIFFTYPE=JDIFFCENTERED; // and only helps rarely and makes 2X slower.
    //    JDIFFTYPE=JDIFFONESIDED; // and only helps rarely and makes 2X slower.
  }
  else{
    JDIFFTYPE=JDIFFONESIDED;
  }


  // ensure uu and pp don't get modified by del-shifts to get Jacobian, which can change primitives to order unity at high radiation gamma
  FTYPE uujac[NPR],ppjac[NPR];
  FTYPE uujacalt[NPR],ppjacalt[NPR];
  int pliter,pl;


  ///////////////
  //
  // Setup which quantitiy iterating
  // Scale-out dimensional stuff before forming predel
  //
  ///////////////
  // for scaling del's norm.  Applies to time component to make as if like space component.
  FTYPE upitoup0[NPR], upitoup0U[NPR], upitoup0P[NPR];
  FTYPE x[NPR],xp[NPR],xjac[2][NPR],xjacalt[NPR];
  FTYPE velmomscale;

  //  FTYPE dimfactU[NPR];
  //  PLOOP(pliter,pl) dimfactU[pl]=1.0; // default
  //  DLOOPA(jjdim) dimfactU[UU+jjdim]=dimfactU[URAD0+jjdim]=sqrt(fabs(ptrgeom->gcon[GIND(jjdim,jjdim)]));
  //  SLOOPA(jjdim) dimfactU[B1+jjdim-1] = 1.0/dimfactU[U1+jjdim-1];

  // U
  // \rho_0 u^t * sqrt(-1/g^{tt})  and most things will be scalars like S u^t as well
  // R^t_nu * sqrt(g^{ii}/g^{tt}) = R^t_orthonu
  PLOOP(pliter,pl) upitoup0U[pl] = dimfactU[pl]*ptrgeom->alphalapse ;

  // P
  PLOOP(pliter,pl) upitoup0P[pl] = 1.0; // comoving quantities
  // v^i / sqrt(g^{ii}) = vortho^i
  SLOOPA(jj) upitoup0P[UU+jj] =upitoup0P[URAD0+jj] = 1.0/dimfactU[UU+jj];


  /////////
  //
  // use uucopy as uu instead of modifying uu in case uu needs adjusting
  //
  /////////
  FTYPE uucopy[NPR],ppcopy[NPR];
  PLOOP(pliter,pl){
    uucopy[pl]=uu[pl];
    ppcopy[pl]=pp[pl];
  }

#define MAXSIGN(x,y) (fabs(x)>fabs(y) ? (x) : (y))
#define MINSIGN(x,y) (fabs(x)<fabs(y) ? (x) : (y))

  ////////
  //
  // check whether origin point is too small relative to uu0
  //
  ////////
  //  JACLOOP(jj,startjac,endjac){
  DLOOPA(jj){ // so all related U's are changed.
    //    if(fabs(uucopy[irefU[jj]])<100.0*NUMEPSILON*fabs(uu0[irefU[jj]])){
    // then set "floor" on uucopy used for Jacobian because otherwise uu is unresolved by f_implicit() error function
    uucopy[irefU[jj]]=MAXSIGN(MAXSIGN(uucopy[irefU[jj]],100.0*NUMEPSILON*uu0[irefU[jj]]),100.0*NUMEPSILON*uup[irefU[jj]]);
    // modify associated pp in reasonable way in case P method since otherwise error will be unresolved using that p associated with that uu
    // only modify primitives associated with those conserved quantities (i.e. gas uu and pp or rad uu and pp)
    ppcopy[irefU[jj]]=MAXSIGN(MAXSIGN(ppcopy[irefU[jj]],100.0*NUMEPSILON*piin[irefU[jj]]),100.0*NUMEPSILON*ppp[irefU[jj]]);
  }



  ///////////
  //
  // setup origin point for difference
  //
  ///////////

  // U
  if(IMPUTYPE(implicititer)){
    PLOOP(pliter,pl) x[pl]=uucopy[pl];
    PLOOP(pliter,pl) xp[pl]=uup[pl];
    PLOOP(pliter,pl) xjac[0][pl]=xjac[1][pl]=xjacalt[pl]=uucopy[pl];
    PLOOP(pliter,pl) upitoup0[pl] = upitoup0U[pl];
    // velmomscale is set such that when (e.g.) T^t_i is near zero, we use T^t_t as reference since we consider T^t_i/T^t_t to be velocity scale that is up to order unity.
    //    velmomscale=pow(fabs(x[irefU[TT]]*upitoup0[irefU[TT]]),1.5);
    velmomscale=pow(fabs(x[irefU[TT]]*upitoup0[irefU[TT]]),1.0); // FUCK
  }
  // P
  else if(IMPPTYPE(implicititer)){
    PLOOP(pliter,pl) x[pl]=ppcopy[pl];
    PLOOP(pliter,pl) xp[pl]=ppp[pl];
    PLOOP(pliter,pl) xjac[0][pl]=xjac[1][pl]=xjacalt[pl]=ppcopy[pl];
    PLOOP(pliter,pl) upitoup0[pl] = upitoup0P[pl];
    // for velocity, assume ortho-scale is order unity (i.e. v=1.0*c)
    //    velmomscale=1.0; // reference scale considered to be order unity  KORALTODO: Not sure if should use something like T^t_i/T^t_t with denominator something more non-zero-ish.
    if(IMPMHDTYPE(implicititer)) velmomscale=MAX(SMALL,sqrt(fabs(x[UU])/MAX(ppp[RHO],ppcopy[RHO]))); // u_g/\rho_0\propto (v/c)^2, so this gives \propto (v/c) .  Leads to more problems for RADTUBE
    else velmomscale=1.0; // FUCK

    // limit
    velmomscale=MIN(1.0,velmomscale);

    //    velmomscale=1.0; // FUCK
  }

  //////////////////////////
  //
  // form pre-del
  //
  //////////////////////////
  FTYPE delspace,deltime,predel[NDIM];

  // get everything in terms of quasi-orthonormal quantities
  FTYPE vsqnorm=0.0;
  JACLOOPALT(jj,startjac,endjac){
    predel[jj] = fabs(x[irefU[jj]]*upitoup0[irefU[jj]]);
    if(irefU[jj]==UU || irefU[jj]==URAD0 || IMPUTYPE(implicititer)){
      // below makes sense because U and ferr are linear in x and same dimensional units
      // TT term (e.g. u_g) can be too small primitive or otherwise, so use maximum of ferr, uu, and x.
      // This avoids issue with Jacobian giving J44[0][0]=0 so that can't iterate u_g because u_g itself is very small.
      predel[jj] = MAX(predel[jj],fabs(uup[irefU[jj]]*upitoup0U[irefU[jj]]));
      predel[jj] = MAX(predel[jj],fabs(uu0[irefU[jj]]*upitoup0U[irefU[jj]]));
      // using error function to normalize is too risky since error is not linear in irefU.  e.g., error is not linear in u_g.
      //      predel[jj] = MAX(predel[jj],fabs(f1[erefU[jj]]*upitoup0U[erefU[jj]]));
      //      predel[jj] = MAX(predel[jj],fabs(f1norm[erefU[jj]]*upitoup0U[erefU[jj]]));
    }
    else{
      // if primitive velocity, then different units than conserved or error function that have energy density scale
    }

    if((jj==1 || jj==2 || jj==3) && IMPPTYPE(implicititer)){
      // v^2 quasi-orthonormal
      vsqnorm += fabs(x[irefU[jj]]*upitoup0[irefU[jj]])*fabs(x[irefU[jj]]*upitoup0[irefU[jj]]);
    }
  }
  // limit
  vsqnorm=MIN(1.0,vsqnorm);

  ///////////
  // form delspace that absorbs all spatial values into a single dimensionless scale to avoid one dimensions smallness causing issues.
  // KORALTODO: Maybe causes problems if Jacobian becomes singular because no change in small velocity-momentum values when should be change.  But can't resolve it really, so should be ok.
  // NOTE: If JACLOOPALT doesn't include momentum, the below doesn't hurt.
  //////////

  // Add all spatial terms in quasi-orthonormal way
  // Also add (u_g)^{1/2} \propto \rho_0(v/c)
  delspace=0.0; SLOOPA(jj) delspace = MAX(delspace,MAX(fabs(predel[jj]) , velmomscale )); // dimensionless-ortho

  //FUCK: Needs to improve and be more general
  if(IMPPMHDTYPE(implicititer)){
    // u_g goes like \rho_0 v^2
    jj=TT; deltime = MAX(fabs(predel[jj]),MAX(ppp[RHO],ppcopy[RHO])*vsqnorm); // FUCK : Leads to more problems for RADTUBE
    //jj=TT; deltime = fabs(predel[jj]);
  }
  else if(IMPPTYPE(implicititer)){
    jj=TT; deltime = fabs(predel[jj]);
  }
  else{
    jj=TT; deltime = fabs(predel[jj]);
  }

  /////////////////
  // back to actual spatial dimension-space scale for application to x and xjac
  SLOOPA(jj) predel[jj] = delspace/upitoup0[irefU[jj]];



  // back to actual time dimension scale
  // NOTE: If JACLOOPALT doesn't include energy, the below doesn't hurt.
  jj=TT;
  predel[jj] = predel[jj]/upitoup0[irefU[jj]];



  FTYPE J[NPR][NPR];
  FTYPE f2[2][NPR],f2norm[2][NPR],f2report[2][NPR],errorabsf2[NUMERRORTYPES];
  set_array(errorabsf2,NUMERRORTYPES,MPI_FTYPE,BIG);

  // set defaults for one-sided case.
  if(JDIFFTYPE==JDIFFONESIDED){
    PLOOP(pliter,pl) f2[0][pl]=f1[pl];
    //  Note that xjac[0] already set correctly.
  }

  int failreturn;
  int fulljaciter=0;
  FTYPE FRACIMPEPSCHANGE=0.1;
  FTYPE del;
  FTYPE IMPEPSSTART=impepsjac;
  while(1){ // ensuring that Jacobian is non-singular if only because del too small (and then if singular, increase del)

    FTYPE localIMPEPS=IMPEPSSTART; // start with fresh del
    
    JACLOOP(jj,startjac,endjac){

      int sided,signside;
      int numsides=2; // fixed at 2
      for(sided=0;sided<numsides;sided++){
        if(JDIFFTYPE==JDIFFONESIDED && sided==0) continue;

        // want chosen energy to have controlled drop, so other energy rises and doesn't drop out of bounds
        if(IMPPTYPE(implicititer)){
          if(sided==1) signside=-1.0;
          else signside=+1.0;
        }
        else{
          if(sided==1) signside=+1.0;
          else signside=-1.0;
        }

        while(1){

          // when |irefU[0]|>>|irefU[1]|, then can't get better than machine error on irefU[0], not irefU[1], so using small del just for irefU[1] makes no sense, so avoid above
          del = localIMPEPS*predel[jj];

         
          // origin point
          PLOOP(pliter,pl) xjac[sided][pl]=xjacalt[pl]=x[pl];

         
          if(JDIFFTYPE==JDIFFONESIDED && sided==1){
            if(*baseitermethod==QTYPMHD || *baseitermethod==QTYPRAD){
              // then choose direction so that decrease u_g and decreases magnitude of \gamma
              if(x[irefU[jj]]>0.0){
                xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del;
                xjacalt[irefU[jj]]=x[irefU[jj]] - signside*del;
              }
              else{
                xjac[sided][irefU[jj]]=x[irefU[jj]] - signside*del;
                xjacalt[irefU[jj]]=x[irefU[jj]] + signside*del;
              }
            }
            else if(*baseitermethod==QTYURAD){
              // ensure R^t_t decreases so Erf doesn't drop out
              if(x[irefU[jj]]>0.0){
                xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del*(jj==TT ? -1.0 : 1.0);
                xjacalt[irefU[jj]]=x[irefU[jj]] - signside*del*(jj==TT ? -1.0 : 1.0);
              }
              else{
                xjac[sided][irefU[jj]]=x[irefU[jj]] - signside*del*(jj==TT ? -1.0 : 1.0);
                xjacalt[irefU[jj]]=x[irefU[jj]] + signside*del*(jj==TT ? -1.0 : 1.0);
              }
            }
            else if(*baseitermethod==QTYUMHD){
              // STILL for this method, ensure R^t_t decreases so Erf doesn't drop out
              if(x[irefU[jj]]>0.0){
                xjac[sided][irefU[jj]]=x[irefU[jj]] - signside*del*(jj==TT ? -1.0 : 1.0);
                xjacalt[irefU[jj]]=x[irefU[jj]] + signside*del*(jj==TT ? -1.0 : 1.0);
              }
              else{
                xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del*(jj==TT ? -1.0 : 1.0);
                xjacalt[irefU[jj]]=x[irefU[jj]] - signside*del*(jj==TT ? -1.0 : 1.0);
              }
            }
            else{
              xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del;
              xjacalt[irefU[jj]]=x[irefU[jj]] - signside*del;
            }
          }
          else{
            // offset xjac (KORALTODO: How to ensure this doesn't have machine precision problems or is good enough difference?)
            xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del; // KORALNOTE: Not sure why koral was using uup or xp here.  Should use x (or uu or pp) because as updated from f_implicit() and uup or ppp for xp hasn't been set yet, so not consistent with desired jacobian or ferr for Newton step.
            xjacalt[irefU[jj]]=x[irefU[jj]] - signside*del;
            //          dualfprintf(fail_file,"NEW: jj=%d del=%g xjac=%g x=%g\n",jj,del,xjac[sided][irefU[jj]],x[irefU[jj]]);
          }

          // set uujac and ppjac using xjac
          if(IMPUTYPE(implicititer)){
            PLOOP(pliter,pl){
              uujac[pl]=xjac[sided][pl];
              uujacalt[pl]=xjacalt[pl];
              ppjac[pl]=ppjacalt[pl]=ppcopy[pl];
            }
          }
          else if(IMPPTYPE(implicititer)){
            PLOOP(pliter,pl){
              uujac[pl]=uujacalt[pl]=uucopy[pl];
              ppjac[pl]=xjac[sided][pl];
              ppjacalt[pl]=xjacalt[pl];
            }
            if(FORCEJDIFFNOCROSS){
              // constrain that which enters jacobian
              int subjj;
              JACLOOP(subjj,startjac,endjac){
                pl=irefU[subjj];
                FTYPE ppmin=100.0*NUMEPSILON*fabs(predel[subjj]);
                //            FTYPE umin=calc_PEQ_ufromTrho(TEMPMIN,ppjac[RHO]);
                //
                if(POSPL(pl)){
                  if(xjac[sided][pl]<ppmin && x[pl]>0.0){
                    xjac[sided][pl]=ppmin;
                    if(debugfail>=3) dualfprintf(fail_file,"Got 1: sided=%d\n",sided);
                  }
                  else if(xjac[sided][pl]>-ppmin && x[pl]<0.0){
                    xjac[sided][pl]=-ppmin;
                    if(debugfail>=3) dualfprintf(fail_file,"Got 2: sided=%d\n",sided);
                  }
                  // now fix ppjac
                  ppjac[pl]=xjac[sided][pl];
                }// end if supposed to normally be positive value quantity
              }
            }
          } 


          // get dUresid for this offset xjac
          int whichcall=FIMPLICITCALLTYPEJAC;
          eomtypelocallocal=*eomtypelocal; // re-default
          int fakeiter=iter;
          int fakef1iter=-1;
          int radinvmod=0; // ignore, assume normal error check will be where this information is used.
          //          dualfprintf(fail_file,"iJ calls f_implicit: iter=%d\n",iter);
          int goexplicitfake; // in Jacobian, so don't try to abort.
          int convreturnf2;
          int whichcapnew=CAPTYPEFIX2; // so Erf stays well-defined even if innaccurate a bit.  So J stays well-defined so iJ doesn't nan-out
          //          PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d ppjac=%21.15g uu0=%21.15g uujac=%21.15g\n",pl,ppjac[pl],uu0[pl],uujac[pl]);
          failreturn=f_implicit(allowbaseitermethodswitch, fakeiter,fakef1iter,failreturnallowableuse, whichcall,localIMPEPS,showmessages,showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocallocal, whichcapnew,itermode, baseitermethod, fracenergy, dissmeasure, &radinvmod, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, realdt, dimtypef, dimfactU, ppjacalt, ppjac,piin,uujacalt, Uiin,uu0,uujac,fracdtG*realdt,ptrgeom,&qjac,f2[sided],f2norm[sided],f2report[sided], &goexplicitfake, &errorabsf2[0], &errorabsf2[1], whicherror, &convreturnf2, nummhdinvsreturn);
          if(failreturn){ // FUCK: Noticed if ==1 (radinv failure), then hits this and tries again.  Maybe costly, and maybe not required with CAPTYPEFIX2 used at least for PMHD,PRAD methods.
            if(showmessages&& debugfail>=2) dualfprintf(fail_file,"f_implicit for f2 failed: jj=%d.  Trying smaller localIMPEPS=%g (giving del=%g) to %g\n",jj,localIMPEPS,del,localIMPEPS*FRACIMPEPSCHANGE);
            localIMPEPS*=FRACIMPEPSCHANGE;
            // try making smaller until no error, unless doesn't work out
            // see if will be able to resolve differences
            int baddiff = fabs(xjac[sided][irefU[jj]]-x[irefU[jj]])/(fabs(xjac[sided][irefU[jj]])+fabs(x[irefU[jj]])) < 10.0*NUMEPSILON;
            if(localIMPEPS<10.0*NUMEPSILON || baddiff){
              // then probably can't resolve difference due to too small 
              if(failreturnallowableuse>=UTOPRIMGENWRAPPERRETURNFAILRAD){
                if(debugfail>=2) dualfprintf(fail_file,"Bad error: f_implicit for f2 failed: jj=%d with localIMPEPS=%g (giving del=%g)\n",jj,localIMPEPS,del);
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
        }// end while(1), checking if ferr failed due to (e.g.) inversion issue.
      }// end sided=0,1
      

      // get Jacobian
      JACLOOP(ii,startjac,endjac) J[erefU[ii]][irefU[jj]]=(f2[1][erefU[ii]] - f2[0][erefU[ii]])/(xjac[1][irefU[jj]]-xjac[0][irefU[jj]]);

      //JACLOOP(ii,startjac,endjac) dualfprintf(fail_file,"NEW: ii=%d jj=%d J=%g : %g %g : %g %g\n",ii,jj,J[erefU[ii]][irefU[jj]], f2[1][erefU[ii]],f2[0][erefU[ii]],xjac[1][irefU[jj]],xjac[0][irefU[jj]]);


      // debug info
      if(debugfail>=2){
        int badcond=0;
        JACLOOP(ii,startjac,endjac){
          if(showmessagesheavy || !isfinite(J[erefU[ii]][irefU[jj]])|| fabs(J[erefU[jj]][irefU[jj]])<SMALL  ) {
            badcond++;
          }
        }
        if(badcond){
          JACLOOP(ii,startjac,endjac){
            dualfprintf(fail_file,"JISNAN: iter=%d startjac=%d endjac=%d ii=%d jj=%d irefU[jj]=%d erefU[ii]=%d : xjac[0]: %21.15g :  xjac[1]: %21.15g : x=%21.15g (del=%21.15g localIMPEPS=%21.15g) : f2[0]=%21.15g f2[1]=%21.15g J=%21.15g : f2norm[0]=%21.15g f2norm[1]=%21.15g\n",
                        iter,startjac,endjac,
                        ii,jj,irefU[jj],erefU[ii],
                        xjac[0][irefU[jj]],
                        xjac[1][irefU[jj]],
                        x[irefU[jj]],
                        del,localIMPEPS,
                        f2[0][erefU[ii]],f2[1][erefU[ii]],J[erefU[ii]][irefU[jj]],
                        f2norm[0][erefU[ii]],f2norm[1][erefU[ii]]
                        );
          }
          JACLOOP(ii,startjac,endjac){
            dualfprintf(fail_file,"NEW: ii=%d jj=%d J=%21.15g : %21.15g %21.15g : %21.15g %21.15g\n",ii,jj,J[erefU[ii]][irefU[jj]], f2[1][erefU[ii]],f2[0][erefU[ii]],xjac[1][irefU[jj]],xjac[0][irefU[jj]]);
          }
          PLOOP(pliter,pl) dualfprintf(fail_file,"JISNAN2: pl=%d ppjac=%21.15g uu0=%21.15g uujac=%21.15g\n",pl,ppjac[pl],uu0[pl],uujac[pl]);
          JACLOOP(ii,startjac,endjac) dualfprintf(fail_file,"ii=%d uucopy=%21.15g uu=%21.15g\n",ii,uucopy[irefU[ii]],uu[irefU[ii]]);
        }//end if badcond>0
      }//end debug


    }
    


    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTJAC: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
      int iii,jjj;
      JACLOOP2D(iii,jjj,startjac,endjac) dualfprintf(fail_file,"J[%d][%d]=%21.15g\n",iii,jjj,J[erefU[iii]][irefU[jjj]]);
    }



    /////////////////////
    //
    //invert Jacobian
    //
    /////////////////////


    // copy over matrix to sub (up to) 4x4 version
    FTYPE Jsub[NDIM][NDIM];
    FTYPE iJsub[NDIM][NDIM];
    JACLOOP2D(ii,jj,startjac,endjac) Jsub[ii][jj]=J[erefU[ii]][irefU[jj]];

    if(endjac-startjac+1==NDIM){
      // inverse *and* transpose index order, so J[erefU][irefU] ->  iJ[irefU][erefU]
      failreturn=inverse_44matrix(Jsub,iJsub);
    }    
    else if(endjac-startjac+1==NDIM-1){ // probably momentum only
      // inverse *and* transpose index order, so J[erefU][irefU] ->  iJ[irefU][erefU]
      failreturn=inverse_33matrix(startjac,endjac,Jsub,iJsub);
    }    
    else if(endjac-startjac+1==1){ // probably energy only
      // inverse *and* transpose index order, so J[erefU][irefU] ->  iJ[irefU][erefU]
      failreturn=inverse_11matrix(startjac,endjac,Jsub,iJsub);
    }    
    
    // copy back inverse matrix from sub-version
    JACLOOP2D(ii,jj,startjac,endjac){
      iJ[irefU[ii]][erefU[jj]]=iJsub[ii][jj]; // iJsub has been transposed from input Jsub, so iJsub[irefU][erefU]
      // dualfprintf(fail_file,"NEW: ii=%d jj=%d iJ=%g\n",ii,jj,iJ[irefU[ii]][erefU[jj]]);
    }


    if(failreturn){

    // debug:
      if(debugfail>=2){
        dualfprintf(fail_file,"Tried to invert Jacobian with %d %d %d %d : %g %g %g : %g %g %g : %d %d : %g %g : %d : %g %g\n",*eomtypelocal, whichcap, itermode, *baseitermethod, fracenergy, dissmeasure, impepsjac, trueimptryconv, trueimptryconvabs, trueimpallowconvabs, trueimpmaxiter, iter, errorabs, errorallabs, whicherror,fracdtG,realdt);
        PLOOP(pliter,pl) dualfprintf(fail_file,"1Tried: pl=%d Uiin=%g uu=%g uup=%g uu0=%g piin=%g pp=%g ppp=%g f1=%g f1norm=%g\n", pl, Uiin[pl], uu[pl], uup[pl], uu0[pl], piin[pl], pp[pl], ppp[pl], f1[pl], f1norm[pl]);
        JACLOOP2D(ii,jj,startjac,endjac){
          dualfprintf(fail_file,"2Tried: ii=%d jj=%d (%d %d) : j=%g iJ=%g\n",ii,jj,startjac,endjac,Jsub[ii][jj],iJsub[ii][jj]);
        }
      }

      // act on failure?

    }// end if failreturn!=0



    /////////////////////////////////////
    //
    // check if Jacobian inversion was successful
    //
    /////////////////////////////////////
    if(failreturn){
      // try increasing localIMPEPS
      IMPEPSSTART/=FRACIMPEPSCHANGE;
      int condnotdiff;
      condnotdiff=IMPEPSSTART > MAXIMPEPS;

      if(condnotdiff){ // KORALTODO: But error relative to x needs to be accounted for!
        if(debugfail>=2) dualfprintf(fail_file,"f_implicit for f2 failed to be different enough from f1 and gave singular Jacobian: IMPEPSSTART=%g (giving del=%g)\n",IMPEPSSTART,del);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC1: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"Jsub[%d][%d]=%21.15g\n",iii,jjj,Jsub[iii][jjj]);
        }
        return(1); // can't expect good derivative above ~0.3, so just return as failure of implicit method.
      }
      else{
        if(debugfail>=2) dualfprintf(fail_file,"inverse_44matrix(J,iJ) failed with eomtypelocallocal=%d, trying IMPEPSSTART=%g :: ijk=%d %d %d\n",eomtypelocallocal,IMPEPSSTART,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC2: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
          int iii,jjj;
          DLOOPA(iii) dualfprintf(fail_file,"predel[%d]=%21.15g\n",iii,predel[iii]);
          DLOOP(iii,jjj) dualfprintf(fail_file,"Jsub[%d][%d]=%21.15g\n",iii,jjj,Jsub[iii][jjj]);
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
        dualfprintf(fail_file,"POSTJAC3: x: %g %g %g %g : x=%g %g %g %g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
        int iii,jjj;
        DLOOP(iii,jjj) dualfprintf(fail_file,"Jsub[%d][%d]=%g\n",iii,jjj,Jsub[iii][jjj]);
      }
      return(1);
    }
  }// end over ensuring Jacobian is non-singular for the given del



  return(0);

}





















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
    jj=TT; Umhd=UUMINLIMIT+fabs(U[UU+jj]);
    jj=TT; Urad=fabs(U[URAD0+jj]);
    idtsub=MAX(idtsub,idtsub0*Urad/Umhd);

    //    dualfprintf(fail_file,"i=%d dtsub1=%g (realdt=%g)\n",ptrgeom->i,1/idtsub,realdt);
       
    // also compare against changed U=U0
    jj=TT; Umhd=UUMINLIMIT+fabs(U0[UU+jj]);
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
    iUmhd=1.0/(fabs(Umhd)+UUMINLIMIT);
    iUrad=1.0/(fabs(Urad)+UUMINLIMIT);
    idtsub=UUMINLIMIT+fabs(Gmhd*iUmhd);
    idtsub=MAX(idtsub,UUMINLIMIT+fabs(Grad*iUrad));
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
    iUsmhd=1.0/(fabs(Usmhd)+UUMINLIMIT);
    iUtmhd=1.0/(fabs(Utmhd)+UUMINLIMIT);
    iUsrad=1.0/(fabs(Usrad)+UUMINLIMIT);
    iUtrad=1.0/(fabs(Utrad)+UUMINLIMIT);
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
      iUmhd=1.0/(fabs(Umhd)+UUMINLIMIT);
      iUrad=1.0/(fabs(Urad)+UUMINLIMIT);
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
      iUmhd=1.0/(fabs(Umhd)+UUMINLIMIT);
      iUrad=1.0/(fabs(Urad)+UUMINLIMIT);
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
    iUmhd=1.0/(fabs(Umhd)+UUMINLIMIT);
    iUrad=1.0/(fabs(Urad)+UUMINLIMIT);
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
    iUsmhd=1.0/(fabs(Usmhd)+UUMINLIMIT);
    iUtmhd=1.0/(fabs(Utmhd)+UUMINLIMIT);
    iUsrad=1.0/(fabs(Usrad)+UUMINLIMIT);
    iUtrad=1.0/(fabs(Utrad)+UUMINLIMIT);
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
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
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
      int whichcap=CAPTYPEBASIC;
      eomtypelocal=*eomtype; // re-default
      int radinvmod=0;
      FTYPE dissmeasure=-1.0; // assume ok to try energy
      int checkoninversiongas;
      int checkoninversionrad;
      // don't check since slows down code and could be good enough solution if original error says ok.
      checkoninversiongas=checkoninversionrad=0;
      //
      int failutoprim=Utoprimgen_failwrapper(doradonly,&radinvmod,showmessages,checkoninversiongas,checkoninversionrad, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Unew, q, ptrgeom, dissmeasure, prnew, &newtonstats);

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
    int whichcap=CAPTYPEBASIC;
    int radinvmod=0;
    FTYPE dissmeasure=-1.0; // assume ok to try energy
    int checkoninversiongas=CHECKONINVERSION;
    int checkoninversionrad=CHECKONINVERSIONRAD;
    int failutoprim=Utoprimgen_failwrapper(doradonly,&radinvmod,showmessages,checkoninversiongas,checkoninversionrad, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, whichcap, EVOLVEUTOPRIM, UNOTHING, Unew, q, ptrgeom, dissmeasure, prnew, &newtonstats);
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
      PLOOP(pliter,pl) pf[pl]=pborig[pl];
      *didreturnpf=1;
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


#if(0)

  // Ru^\mu = R^\mu_\nu u^\nu
  // Ruu = R^a_b u_a u^b

  // Ruuu^\mu = Ru^\mu + Ruu u^\mu
  //          = R^\mu_c u^c + R^a_b u_a u^b u^\mu

  // Ruuu^t = R^t_t u^t + R^t_t u_t u^t u^t   + R^t_i u^i   + R^i_t u_i u^t u^t + R^t_j u_t u^j u^t

  // Ruuu^t = R^t_t u^t (1 + u_t u^t)         + Rus + Ruuss

  // (1 + u_t u^t) = 1 + u^t (u^t g_{tt} + u^i g_{ti}) = 1+ u^t^2 g_{tt} + u^t u^i g_{ti} = 1 - (\gamma/\alpha)^2 (-g_{tt}) + u^t u^i g_{ti}

  // = 1 - (\gamma/\alpha)^2 (-(1 + g_{tt} -1)) = 1 - (\gamma/\alpha)^2 + (\gamma/\alpha)^2 (1+g_{tt})

  // 

  // get R^t_t u^t + (R^t_t u^t u_t)u^t and avoid catastrophic cancellation
  FTYPE Ruuss=0.; DLOOP(i,j) if(i!=TT && j!=TT) Ruuss+=Rij[i][j]*ucov[i]*ucon[j];
  // FUCK: Check again.
  FTYPE fact=(-ptrgeom->gcov[GIND(TT,TT)])*(-ptrgeom->gcon[GIND(TT,TT)]);
  FTYPE fact2=1.0-fact;
  FTYPE utildecon[NDIM]={0.0,pp[URAD1],pp[URAD2],pp[URAD3]};
  FTYPE utsq = 0.0,utildecov[NDIM]; lower_vec(utildecon,ptrgeom,utildecov); SLOOPA(j) utsq+=utildecon[j]*utildecov[j];
  FTYPE ucontucovt = ( fact2 - utsq*fact + ucon[TT]*(ucon[1]*ptrgeom->gcov[GIND(TT,1)]+ucon[2]*ptrgeom->gcov[GIND(TT,2)]+ucon[3]*ptrgeom->gcov[GIND(TT,3)]));
  FTYPE Rut=Rij[TT][TT]*ucon[TT] * ucontucovt;
  FTYPE Rus;
#endif

  FTYPE Ru,Ruuu,term1,term2,term3;
  DLOOPA(i){
    Ru=0.; DLOOPA(j) Ru+=Rij[i][j]*ucon[j];
#if(1)
    Ruuu=(Ru + Ruu*ucon[i]);
#else
    if(i!=TT) Ruuu=(Ru + Ruu*ucon[i]);
    else{
      Rus=0.; DLOOPA(j) if(j!=TT) Rus+=Rij[i][j]*ucon[j];
      Ruuu=Ruuss + Rus + Rut;
    }
#endif

    // group by independent terms
    term1 = -(kappa*Ru + lambda*ucon[i]);
    term2 = -kappaes*Ruuu;

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
  *chieffreturn=chi + lambda/(ERADLIMIT+fabs(pp[PRAD0])); // if needed


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
  extern int simplefast(int whichcall, int dir, struct of_geom *geom,struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax);

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
  simplefast(0, dir,geom,q,vrad2,vmin,vmax); // simplefast(0) means first call rather than an resursive attempt.

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
    DLOOPA(jj){
      term1[jj]=THIRD*pr[PRAD0]*(4.0*q->uradcon[dir]*q->uradcov[jj]);
      term2[jj]=THIRD*pr[PRAD0]*(delta(dir,jj));
      radstressdir[jj]=term1[jj]+term2[jj];
      if(radstressdirabs!=NULL) radstressdirabs[jj]=fabs(term1[jj])+fabs(term2[jj]);
    }
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
    if(debugfail>=2) dualfprintf(fail_file,"idet (det=%26.15g idet=%26.15g) in inverse 4x4 zero or nan\n",det,idet);
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
    if(debugfail>=2) dualfprintf(fail_file,"inverse_33matrix got singular det=%g idet=%g\n",det,idet);
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
    if(debugfail>=2) dualfprintf(fail_file,"inverse 1x1 zero or nan\n",a[sj][sj]);
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
  u2p_rad(showmessages, allowlocalfailurefixandnoreport, GAMMAMAXRAD, CAPTYPEBASIC, U, pout, ptrgeomtouse, &lpflag, &lpflagrad);

  //  DLOOPA(jj) dualfprintf(fail_file,"u2p_rad: jj=%d pout=%g\n",jj,pout[PRAD0+jj]);



  // get back to whichvel
  FTYPE uconback[NDIM],othersback[NUMOTHERSTATERESULTS];
  // for fluid
  ucon_calc_whichvel(VELREL4,pout,ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,pout);
  // KORALTODO: for radiation (always returned as VELREL4 so far).
  ucon_calc_whichvel(VELREL4,&pout[URAD1-U1],ptrgeomtouse,uconback,othersback);
  ucon2pr(*whichvel,uconback,ptrgeomtouse,&pout[URAD1-U1]);



  // DEBUG:
  if(lpflag!=UTOPRIMNOFAIL || lpflagrad!=UTOPRIMRADNOFAIL){ // DEBUG with 1||
    // allows fixups to be applied, such as gamma radiation limiting
    if(lpflag>UTOPRIMNOFAIL || lpflagrad>UTOPRIMRADNOFAIL){ // DEBUG with 1||
      if(debugfail>=2){
        dualfprintf(fail_file,"Failed to invert during prad_fftolab().  Assuming fixups won't be applied: %d %d\n",lpflag,lpflagrad);
        dualfprintf(fail_file,"ijk=%d %d %d : %d\n",ptrgeomtouse->i,ptrgeomtouse->j,ptrgeomtouse->k,ptrgeomtouse->p);
        PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pin=%g U=%g\n",pl,pin[pl],U[pl]);
        DLOOPA(jj) dualfprintf(fail_file,"jj=%d ucon=%g\n",jj,ucon[jj]);
        DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d Rijff=%g Rijlab=%g\n",jj,kk,Rijff[jj][kk],Rijlab[jj][kk]);
        DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d gcov=%g gcon=%g\n",jj,kk,ptrgeomtouse->gcov[GIND(jj,kk)],ptrgeomtouse->gcon[GIND(jj,kk)]);
        PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pout=%g\n",pl,pout[pl]);
      }
      myexit(189235);
      // KORALTODO: Check whether really succeeded?  Need to call fixups?  Probably, but need per-cell fixup.  Hard to do if other cells not even set yet as in ICs.  Should probably include fixup process during initbase.c stuff.
    }
    else{
      if(debugfail>=2){
        dualfprintf(fail_file,"Fixups applied during invert during prad_fftolab(): %d %d\n",lpflag,lpflagrad);
        dualfprintf(fail_file,"ijk=%d %d %d : %d\n",ptrgeomtouse->i,ptrgeomtouse->j,ptrgeomtouse->k,ptrgeomtouse->p);
      }
    }
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
  int u2p_rad_new(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, int whichcap, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad);
  int u2p_rad_orig(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad);
  int u2p_rad_new_pre(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad);
  int toreturn;
  int pliter,pl;
  FTYPE prorig[NPR];


  ///////////////
  //
  // CHECK if should abort inversion attempt if already dropped-out value in uu
  //
  //////////////
  if(0){
    // put in a catch for when inputted uu[URAD0] has dropped-out already and don't try to invert.
    if(fabs(uu[URAD0]<=2.0*10.0*ERADLIMIT)){ // often 10*ERADLIMIT is used to set as above ERADLIMIT, so here a higher catch is 2*10*ERADLIMIT
      // force to be reasonable
      // currently always return WHICHVEL=VELREL4, so just set to floor values
      pin[URAD0]=ERADLIMIT;
      int jj;
      SLOOPA(jj) pin[URAD1+jj-1] = 0.0;
      return(0);
    }
  }


  ///////////////
  //
  // NORMAL ATTEMPT where we compute inversion
  //
  //////////////


  // store orig
  PLOOP(pliter,pl) prorig[pl] = pin[pl];

#if(WHICHU2PRAD==0)
  toreturn=u2p_rad_orig(showmessages, allowlocalfailurefixandnoreport, gammamaxrad, uu, pin, ptrgeom,lpflag, lpflagrad);
#else
  //toreturn=u2p_rad_new_pre(showmessages, allowlocalfailurefixandnoreport, gammamaxrad, uu, pin, ptrgeom,lpflag, lpflagrad);
  toreturn=u2p_rad_new(showmessages, allowlocalfailurefixandnoreport, gammamaxrad, whichcap, uu, pin, ptrgeom,lpflag, lpflagrad);
#endif

  int caughtnan=0;
  if(!finite(pin[URAD0]) || !finite(pin[URAD1]) || !finite(pin[URAD2]) || !finite(pin[URAD3])){
    // FUCK: Shouldn't happen, but does on Kraken
    caughtnan=1;
  }

  if(caughtnan){

    if(debugfail>=2){
      dualfprintf(fail_file,"u2p_rad() generated nan result: %d %d %g %d\n",showmessages, allowlocalfailurefixandnoreport, gammamaxrad, whichcap);
      PLOOP(pliter,pl) dualfprintf(fail_file,"u2p_rad: pl=%d prorig=%21.15g uu=%21.15g pin=%21.15g\n",pl,prorig[pl],uu[pl],pin[pl]);
    }


    if(0){ // FUCK: if doing iterations, need to let fail with nan so aborts and stops trying right away.  Otherwise huge waste.
      // force to be reasonable
      // currently always return WHICHVEL=VELREL4, so just set to floor values
      pin[URAD0]=ERADLIMIT;
      int jj;
      SLOOPA(jj) pin[URAD1+jj-1] = 0.0;
    }

  }


  return(toreturn);
}



///////////////
//
// Like u2p_rad_orig(), but uses Jon's paper draft ZAMO RAD version
//
//////////////
int u2p_rad_new_pre(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  static long long int numyvarneg,numyvarbig,numErneg,nummod;

  if(WHICHVEL!=VELREL4){
    dualfprintf(fail_file,"u2p_rad() only setup for relative 4-velocity, currently.\n");
    myexit(137432636);
  }


  // copy over pin so pin isn't modified until end
  int pliter,pl;
  FTYPE pp[NPR];
  PLOOP(pliter,pl) pp[pl]=pin[pl];

  //////////////////////
  //
  // Prepare inversion from U->p for radiation assuming M1 closure
  //
  //////////////////////

  *lpflagrad=UTOPRIMRADNOFAIL;

  FTYPE Er,Utildesq,Utildecon[NDIM];
  compute_ZAMORAD(uu, ptrgeom, &Er, &Utildesq, Utildecon);
  


  FTYPE Ersq,yvar;
  int didmod=0;

  //  if(1||gammamaxrad>0.9*GAMMAMAXRADIMPLICITSOLVER || Er>=ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
  if(Er>=ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
    // Er^2
    Ersq=Er*Er;
    // y
    yvar = Utildesq / (ERADLIMIT*ERADLIMIT + Ersq); // ERADLIMIT*ERADLIMIT better be machine representable in case Ersq is not.
  }
  else{// then bad solution
    //    dualfprintf(fail_file,"Er=%26.20g<ERADLIMIT=%26.20g yvar=%26.20g Utildesq=%26.20g Ersq=%26.20g\n",Er,ERADLIMIT,yvar,Utildesq,Ersq);
    Ersq=ERADLIMIT*ERADLIMIT;
    yvar = 0.0;
    didmod=1;
    numErneg++;
  }

  //  dualfprintf(fail_file,"Er=%g Utildesq=%g\n",Er,Utildesq);

  // \gamma_{\rm rad}^2 :  only 1 root
  FTYPE gammasq,gamma;
  FTYPE gammamax=gammamaxrad;
  FTYPE gammamaxsq=gammamax*gammamax;
  FTYPE ylimit = 16.0*gammamaxsq*(gammamaxsq-1.0)/((4.0*gammamaxsq-1.0)*(4.0*gammamaxsq-1.0));
  if(yvar<0.0){
    //    dualfprintf(fail_file,"Er=%26.20g yvar=%26.20g<0.0 Utildesq=%26.20g Ersq=%26.20g\n",Er,yvar,Utildesq,Ersq);
    yvar=0.0;
    gammasq = 1.0;
    gamma = 1.0;
    didmod=1;
    numyvarneg++;
  }
  else if(yvar>ylimit){ // beyond gamma limit, then rescale gamma
    yvar=ylimit;
    gammasq = gammamaxsq;
    gamma = gammamax;
    didmod=1; *lpflagrad=UTOPRIMRADFAILCASE2A; // used to detec if modified primitives to not be consistent with inputted uu
    numyvarbig++;
    //    dualfprintf(fail_file,"yvar=%g>%g Ersq=%g gamma=%g\n",yvar,ylimit,Ersq,gamma);
  }
  else{ // normal solution
    gammasq = (2.0 - yvar + sqrt(4.0-3.0*yvar))/ (4.0*(1.0-yvar));
    gamma=sqrt(gammasq);
    //    dualfprintf(fail_file,"yvar=%g Ersq=%g gamma=%g\n",yvar,Ersq,gamma);
  }

  FTYPE Erf;
  FTYPE urfconrel[NDIM]={0.0};
  int jj;
  //  if(1||gammamaxrad>0.9*GAMMAMAXRADIMPLICITSOLVER || Er>=ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
  if(Er>=ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
    // now obtain primitives
    FTYPE pr = Er/(4.0*gammasq-1.0);
    // radiation frame energy density
    //    Erf = pr/(4.0/3.0-1.0);
    Erf = 3.0*pr;
    
    // radiation frame relativity 4-velocity
    SLOOPA(jj) urfconrel[jj] = (Utildecon[jj]/(4.0*pr*gamma));
  }
  else{
    Erf = ERADLIMIT;
    gamma=1.0;
    // radiation frame relativity 4-velocity
    SLOOPA(jj) urfconrel[jj] = 0.0;
    didmod=1; *lpflagrad=UTOPRIMRADFAILERFNEG; // used to detect if modified primitives to not be consistent with inputted uu
  }


  /////////////////
  //
  //new primitives (only uses urfcon[1-3])
  //
  /////////////////
  pin[PRAD0]=Erf;
  pin[PRAD1]=urfconrel[1];
  pin[PRAD2]=urfconrel[2];
  pin[PRAD3]=urfconrel[3];

  //  dualfprintf(fail_file,"didmod=%d\n",didmod);

  // make sure E_r no larger than starting value
  if(didmod==1){
    nummod++;
    

    // First, ensure \gamma correct
    FTYPE gammanew,qsqnew;
    gamma_calc_fromuconrel(&pin[URAD1-1],ptrgeom,&gammanew,&qsqnew);
    //    dualfprintf(fail_file,"gamma=%g gammanew=%g\n",gamma,gammanew);

    // rescale, assuming want to be gamma
    FTYPE fvar=sqrt((gamma*gamma-1.0)/(gammanew*gammanew-1.0));
    if(gammanew>1.0){
      SLOOPA(jj) urfconrel[jj] *= fvar;
    }
    else{
      SLOOPA(jj) urfconrel[jj] *= 0.0;
    }
    //    dualfprintf(fail_file,"urfconrel=%g %g %g\n",urfconrel[1],urfconrel[2],urfconrel[3]);

    
    // new prims
    pin[PRAD1]=urfconrel[1];
    pin[PRAD2]=urfconrel[2];
    pin[PRAD3]=urfconrel[3];

    if(0){ // causes more problems for implicit solver.

      // Second, ensure not creating energy in ZAMO frame
      struct of_state q;
      get_state_uradconuradcovonly(pin, ptrgeom, &q);
      FTYPE Rtnu[NDIM];
      mhd_calc_rad( pin, TT, ptrgeom, &q, Rtnu, NULL );
      //    dualfprintf(fail_file,"Rtnu=%g %g %g %g\n",Rtnu[0],Rtnu[1],Rtnu[2],Rtnu[3]);
      FTYPE Ernew,Utildesqnew,Utildeconnew;
      compute_ZAMORAD(&Rtnu[0-URAD0], ptrgeom, &Ernew, &Utildesqnew, &Utildeconnew); // out of range warning ok.
      Erf = Erf*MIN(1.0,Er/Ernew);
    }

    *lpflagrad=UTOPRIMRADFAILCASE2A;

    // could continue iterating, or should find closed form expressions for all this.
    //    dualfprintf(fail_file,"Ernew=%g Utildesqnew=%g\n",Ernew,Utildesqnew);
  }


  /////////////////
  //
  // really new primitives (only uses urfcon[1-3])
  //
  /////////////////
  pin[PRAD0]=Erf;
  pin[PRAD1]=urfconrel[1];
  pin[PRAD2]=urfconrel[2];
  pin[PRAD3]=urfconrel[3];


  if(debugfail>=2){
    static long int nstepold=-1;
    if(nstep!=nstepold && nstep%100==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0 && steppart==0){
      nstepold=nstep;
      dualfprintf(fail_file,"numyvarneg=%lld numyvarbig=%lld numErneg=%lld nummod=%lld : nstep=%ld\n",numyvarneg,numyvarbig,numErneg,nummod,nstep);
    }
  }


  if(DORADFIXUPS==1 || allowlocalfailurefixandnoreport==0){
    // KORALTODO: Problem is fixups can average across shock or place where (e.g.) velocity changes alot, and averaging diffuses shock and can leak-out more failures.
  }
  else{
    // CASE reductions (so set as no failure so fixups don't operate -- but might also want to turn off CHECKINVERSIONRAD else that routine won't know when to ignore bad U->P->U cases.)
    if(*lpflagrad==0) *lpflagrad=UTOPRIMRADNOFAIL;
    else *lpflagrad=UTOPRIMRADFAILFIXEDUTOPRIMRAD; //UTOPRIMRADNOFAIL;
  }

  return 0;



}


///////////////
//
// Like u2p_rad_orig(), but uses Jon's paper draft ZAMO RAD version
//
//////////////
int u2p_rad_new(int showmessages, int allowlocalfailurefixandnoreport, FTYPE gammamaxrad, int whichcap, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  static long long int numyvarneg,numyvarbig,numErneg,nummod;

  if(WHICHVEL!=VELREL4){
    dualfprintf(fail_file,"u2p_rad() only setup for relative 4-velocity, currently.\n");
    myexit(137432636);
  }


  // copy over pin so pin isn't modified until end
  int pliter,pl;
  FTYPE pp[NPR];
  PLOOP(pliter,pl) pp[pl]=pin[pl];

  //////////////////////
  //
  // Prepare inversion from U->p for radiation assuming M1 closure
  //
  //////////////////////

  *lpflagrad=UTOPRIMRADNOFAIL;

  FTYPE Er,Utildesq,Utildecon[NDIM];
  compute_ZAMORAD(uu, ptrgeom, &Er, &Utildesq, Utildecon);
  
  // \gamma_{\rm rad}^2 :  only 1 root
  FTYPE gammasq,gamma,qsq;
  FTYPE gammamax=gammamaxrad;
  FTYPE gammamaxsq=gammamax*gammamax;
  FTYPE ylimit = 16.0*gammamaxsq*(gammamaxsq-1.0)/((4.0*gammamaxsq-1.0)*(4.0*gammamaxsq-1.0));

  FTYPE yvar;
  int didmod=0;
  int didmodEr=0,didmody=0;

  ///////////////////
  //
  // Get y
  //
  ///////////////////
  if(Er>ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
    // E_r^2
    FTYPE Ersq=Er*Er;
    // y
    yvar = Utildesq / (ERADLIMIT*ERADLIMIT+Ersq);
  }
  else{// then bad solution
    //    dualfprintf(fail_file,"Er=%26.20g<ERADLIMIT=%26.20g yvar=%26.20g Utildesq=%26.20g Ersq=%26.20g\n",Er,ERADLIMIT,yvar,Utildesq,Ersq);
    Er=ERADLIMIT;
    yvar = ylimit; // used
    didmod=1;
    didmodEr=1;
    numErneg++;
  }


  ///////////////////
  //
  // Get \gamma and \gamma^2 from y
  //
  ///////////////////
  FTYPE Erf,pr;
  if(yvar<0.0){
    //    dualfprintf(fail_file,"Er=%26.20g yvar=%26.20g<0.0 Utildesq=%26.20g Ersq=%26.20g\n",Er,yvar,Utildesq,Ersq);
    gammasq = 1.0;
    gamma = 1.0;
    if(didmodEr) didmod=didmody=1;
    //    didmod=1; // assume when y<0 that don't need to modify how Erf computed (i.e. Er is ok)
    numyvarneg++;
    pr = Er/(4.0*gammasq-1.0);
    // radiation frame energy density
    Erf = pr/(4.0/3.0-1.0);
  }
  else if(yvar>ylimit){ // beyond gamma limit, then rescale gamma
    gammasq = gammamaxsq;
    gamma = gammamax;
    didmod=1;
    didmody=1;
    numyvarbig++;

    pr = Er/(4.0*gammasq-1.0);
    // radiation frame energy density
    Erf = pr/(4.0/3.0-1.0);
    //    dualfprintf(fail_file,"yvar=%g>%g Ersq=%g gamma=%g\n",yvar,ylimit,Ersq,gamma);
  }
  else{ // normal solution
    gammasq = (2.0 - yvar + sqrt(4.0-3.0*yvar))/ (4.0*(1.0-yvar));
    gamma=sqrt(gammasq);

    pr = Er/(4.0*gammasq-1.0);
    // radiation frame energy density
    Erf = pr/(4.0/3.0-1.0);
  }



  ///////////////////
  //
  // Get uconrel and Erf from gamma,gammasq, and Utildecon
  //
  ///////////////////
  FTYPE urfconrel[NDIM]={0.0};
  int jj;
  //  if(1||gammamaxrad>0.9*GAMMAMAXRADIMPLICITSOLVER || Er>=ERADLIMIT){ // then good solution.  Avoid caps during implicit solver to allow progress on solution in smooth way.
  if(didmody==0 && didmodEr==0){ // then good solution
   
    // radiation frame relativity 4-velocity
    SLOOPA(jj) urfconrel[jj] = gamma*(Utildecon[jj]/(4.0*pr*gammasq));

    //if(startpos[1]+ptrgeom->i==131 && startpos[2]+ptrgeom->j==19) dualfprintf(fail_file,"0didmod=%d : urfconrel=%g %g %g : %g : pr=%g Er=%g Ersq=%g yvar=%g Utildesq=%g\n",didmod,urfconrel[1],urfconrel[2],urfconrel[3],gamma,pr,Er,Er*Er,yvar,Utildesq);
  }
  else{ // fixes in case when Er<ERADLIMIT (whether or not y is modified)
    if(whichcap==CAPTYPEFIX1 || whichcap==CAPTYPEFIX2){

      // override if Er<0 originally since otherwise out of control rise in Erf
      //      if(didmodEr || yvar<-NUMEPSILON || yvar>=2.0){
      if(didmodEr && whichcap==CAPTYPEFIX1){
        Er = ERADLIMIT;
        Erf = ERADLIMIT;
        gamma=1.0;
        // radiation frame relativity 4-velocity
        SLOOPA(jj) urfconrel[jj] = 0.0;
        // If E_r<0, then can't trust E_r and Erf is arbitrary.  Choose instead first urfconrel to be maximum gamma as pointed in same 4-velocity direction as from Utildecon over previous Erf or u_g just to set scale
        // But better to not trust to avoid run-away energy creation
      }
      else{
        //        if(yvar<1.0&&0||1){
        if(yvar<1.0&&0){
          // below sticks to what would end up like normal velocity scale.
          SLOOPA(jj) urfconrel[jj] = gamma*(Utildecon[jj]/(4.0*pr*gammasq));
          if(yvar>=1.0) SLOOPA(jj) urfconrel[jj]*=gamma; // fake to pretend very high gamma that we are coming from.          
        }
        else{
          // if yvar>=1.0, then can't trust that gamma will rescale urfconrel into proper gamma>>1 range, so force.
          // Get urfconrel in same direction as Utildecon
          FTYPE Utildeabs=0.5*(sqrt(fabs(Utildesq))+fabs(Er)+ERADLIMIT);
          SLOOPA(jj) urfconrel[jj] = gamma*Utildecon[jj]/Utildeabs; // gives something that gives back ~gamma when using gamma_calc_fromuconrel()
        }

        // now get gamma for this fake urfconrel that is so-far only very roughly expected to be correct.
        FTYPE gammanew,qsqnew;
        gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammanew,&qsqnew);
        //        if(gammanew<gammamax || gammanew<10.0){
        //        if(gammanew<10.0 || gammanew>1000.0){
        //          dualfprintf(fail_file,"gamma=%g gammanew=%g %g %g %g\n",gamma,gammanew,Utildecon[1],Utildecon[2],Utildecon[3]);
        //        }

        // rescale, assuming want to be gammamax
        if(gammanew>1.0){
          FTYPE fvar=sqrt((gammamax*gammamax-1.0)/(gammanew*gammanew-1.0));
          SLOOPA(jj) urfconrel[jj] *= fvar;
          // verify
          FTYPE gammaneworig=gammanew;
          gamma_calc_fromuconrel(urfconrel,ptrgeom,&gammanew,&qsqnew);
          //          dualfprintf(fail_file,"VERIFY: gamma=%g gammanew=%g->%g fvar=%g\n",gamma,gammaneworig,gammanew,fvar);
          gamma=gammanew;
          gammasq=gammanew*gammanew;
          qsq=qsqnew;
        }
        else{
          SLOOPA(jj) urfconrel[jj] *= 0.0;
          gamma=1.0;
          gammasq=1.0;
          qsq=0.0;
        }


        // Get Er independent of R^t_t when hit limits when wanting Jacobian, so Jacobian doesn't NAN-out when taking differences and drop-outs in Erf lead to err function dominated by uu0 (and less often, G) and unaffected by the uu we are changing.
        if(whichcap==CAPTYPEFIX2 || yvar>=1.0-100.0*NUMEPSILON){ // don't use, just trust Er>0.  Kinda works to set to if(1), but leads to some artifacts near where suddenly R^t_t would give different result for Er.
          //FTYPE utildesq=1.0+qsq;
          //pr=gamma*Utildesq/(4.0*gammasq) / (utildesq);
          
          // This determination of Er and pr connects continuously with y<=ylimit case no matter what original Er was.
          Er = ERADLIMIT + sqrt(fabs(Utildesq)/ylimit);
          pr=Er/(4.0*gammasq-1.0);
          // Get Erf
          FTYPE Erforig=Erf;
          Erf = pr/(4.0/3.0-1.0);

          // only modify Erf if really used as solution, not just Jacobian
          if(Erforig<Erf && whichcap==CAPTYPEFIX1) Erf=Erforig;
        }
      }

    }// endif capfixtype1
    else if(whichcap==CAPTYPEBASIC){


      if(didmodEr){
        // This just rejects entire radiative solution and makes it up.  Bit extreme, but works.  But leaves Erf having lowest values in spots where radiation just slightly went beyond speed of light.
        // However, this is most reasonable since if E_r<0, that means radiative energy is in another cell.  If fill this cell with gammamax version of Erf using Utildesq, then adding energy, and situation can blow-up fast.
        Erf = ERADLIMIT;
        gamma=1.0;
        // radiation frame relativity 4-velocity
        SLOOPA(jj) urfconrel[jj] = 0.0;
      }
      else{
        // radiation frame relativity 4-velocity (using pr>0 and chosen gamma,gammasq)
        SLOOPA(jj) urfconrel[jj] = gamma*(Utildecon[jj]/(4.0*pr*gammasq));

        // Get resulting gamma and fix \tilde{u}^i, but don't modify Erf.
        FTYPE gammanew,qsqnew;
        gamma_calc_fromuconrel(&pin[URAD1-1],ptrgeom,&gammanew,&qsqnew);
        //    dualfprintf(fail_file,"didmod: gamma=%g gammanew=%g\n",gamma,gammanew);
        
        // rescale, assuming want to be gamma that chose in previous section
        if(gammanew>1.0){
          FTYPE fvar=sqrt((gammasq-1.0)/(gammanew*gammanew-1.0));
          SLOOPA(jj) urfconrel[jj] *= fvar;
        }
        else{
          SLOOPA(jj) urfconrel[jj] *= 0.0;
        }

        // don't modify Erf if E_r was positive so that only had to rescale y (even if rescaled down from y>1 or gamma^2<0)
        //        pr = Er/(4.0*gammanew*gammanew-1.0);
        //        // radiation frame energy density
        //        Erf = pr/(4.0/3.0-1.0);


        //    if(startpos[1]+ptrgeom->i==131 && startpos[2]+ptrgeom->j==19) dualfprintf(fail_file,"didmod=%d : urfconrel=%g %g %g : %g %g\n",didmod,urfconrel[1],urfconrel[2],urfconrel[3],gamma,gammanew);
      }
    }
    else{
      dualfprintf(fail_file,"No such whichcap=%d\n",whichcap);
      myexit(234534634);
    }
  }

  //  if(startpos[1]+ptrgeom->i==131 && startpos[2]+ptrgeom->j==19){
  //    dualfprintf(fail_file,"AFTER urfconrel=%g %g %g : %g\n",urfconrel[1],urfconrel[2],urfconrel[3],gamma);
  //  dualfprintf(fail_file,"AFTER2 urfconrel2=%g %g %g : %g\n",urfconrel[1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)])),urfconrel[2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)])),urfconrel[3]*sqrt(fabs(ptrgeom->gcov[GIND(3,3)])),gamma);
  // }


  /////////////////
  //
  //new primitives (only uses urfcon[1-3])
  //
  /////////////////
  pin[PRAD0]=Erf;
  pin[PRAD1]=urfconrel[1];
  pin[PRAD2]=urfconrel[2];
  pin[PRAD3]=urfconrel[3];

  //  dualfprintf(fail_file,"didmod=%d\n",didmod);

  // make sure E_r no larger than starting value
  if(didmod==1){
    nummod++;
    if(didmodEr){
      *lpflagrad=UTOPRIMRADFAILERFNEG; // used to detec if modified primitives to not be consistent with inputted uu
    }
    else{
      *lpflagrad=UTOPRIMRADFAILCASE2A; // used to detec if modified primitives to not be consistent with inputted uu
    }
  }


  if(debugfail>=2){
    static long int nstepold=-1;
    if(nstep!=nstepold && nstep%100==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0 && steppart==0){
      nstepold=nstep;
      dualfprintf(fail_file,"numyvarneg=%lld numyvarbig=%lld numErneg=%lld nummod=%lld : nstep=%ld\n",numyvarneg,numyvarbig,numErneg,nummod,nstep);
    }
  }


  if(DORADFIXUPS==1 || allowlocalfailurefixandnoreport==0){
    // KORALTODO: Problem is fixups can average across shock or place where (e.g.) velocity changes alot, and averaging diffuses shock and can leak-out more failures.
  }
  else{
    // CASE reductions (so set as no failure so fixups don't operate -- but might also want to turn off CHECKINVERSIONRAD else that routine won't know when to ignore bad U->P->U cases.)
    if(*lpflagrad==0) *lpflagrad=UTOPRIMRADNOFAIL;
    else *lpflagrad=UTOPRIMRADFAILFIXEDUTOPRIMRAD; //UTOPRIMRADNOFAIL;
  }

  return 0;



}

// compute ZAMO version of radiation quantities as in paper draft
static int compute_ZAMORAD(FTYPE *uu, struct of_geom *ptrgeom, FTYPE *Er, FTYPE *Utildesq, FTYPE *Utildecon)
{
  int jj,kk;
  FTYPE etacov[NDIM],etacon[NDIM];
  FTYPE Ucon[NDIM],Ucov[NDIM],Utildecov[NDIM];

  //  dualfprintf(fail_file,"uu=%g %g %g %g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3]);

  // \eta_\mu
  etacov[TT] = -ptrgeom->alphalapse;
  SLOOPA(jj) etacov[jj]=0.0;
  // \eta^\mu
  raise_vec(etacov,ptrgeom,etacon); // could use ptrgeom->beta

  // U_\mu = -R^\nu_\mu \eta_\nu = \alpha R^t_\mu
  DLOOPA(jj) Ucov[jj] = ptrgeom->alphalapse*uu[URAD0+jj];
  // U^\mu
  raise_vec(Ucov,ptrgeom,Ucon);

  // \tilde{U}^\mu = j^\mu_\nu U^\nu = (\delta^\mu_\nu + \eta^\mu \eta_\nu) U^\nu : ZAMO frame momentum
  DLOOPA(jj) Utildecon[jj]=0.0;
  DLOOP(jj,kk) Utildecon[jj] += (delta(jj,kk) + etacon[jj]*etacov[kk])*Ucon[kk];
  // \tilde{U}_\mu
  lower_vec(Utildecon,ptrgeom,Utildecov);
  
  // \tilde{U}^2 = \tilde{U}^\mu \tilde{U}_\mu
  *Utildesq=0.0;
  DLOOPA(jj) *Utildesq += Utildecon[jj]*Utildecov[jj];

  // -Er = -R^\nu_\mu \eta_\nu \eta^\mu = U_\mu \eta^\mu = alpha R^t_\mu \eta^\mu : ZAMO frame energy
  *Er=0.0;
  DLOOPA(jj) *Er += -Ucov[jj]*etacon[jj];


  return(0);
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
    if(*lpflagrad==0) *lpflagrad=UTOPRIMRADNOFAIL;
    else *lpflagrad=UTOPRIMRADFAILFIXEDUTOPRIMRAD; //UTOPRIMRADNOFAIL;
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

#if(0)
    //    if(Erf<ERADLIMIT) Erf=ERADLIMIT; // case when velocity fine and probably just Erf slightly negative
    //    if(Erf<ERADLIMIT) Erf=MAX(NUMEPSILON*fabs(Erf),ERADLIMIT); // case when velocity fine and probably just Erf slightly negative
    if(Erf<ERADLIMIT) Erf=Erf0; // case when velocity fine and probably just Erf slightly negative
#endif

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
      //dualfprintf(fail_file,"THIS ONE1\n");
    }
    else if(Avcov[TT]>0.0){
      //      Erf=ERADLIMIT;
      Erf=Erf0;
      SLOOPA(jj) urfconrel[jj]=0.0;
      //dualfprintf(fail_file,"THIS ONE2\n");
    }
    else{
      //    Erf=ERADLIMIT;
      //    Erf=MAX(MIN(Erf,Erforig),ERADLIMIT);
      if(gammarel2orig>=1.0 && isfinite(gammarel2orig) && isfinite(Erforig)){
        Erf=MAX(MIN(Erf,Erforig),Erf0);
        //        Erf=ERADLIMIT;
        //        dualfprintf(fail_file,"THIS ONE3\n");
      }
#define AVCOVRELDIFFALLOWED 1E-2 // KORALTODO: only use new "cold" solution if relatively close Avcov.
      else if(fabs(Avcovorig[TT]-Avcov[TT])/fabs(fabs(Avcovorig[TT])+fabs(Avcov[TT]))<AVCOVRELDIFFALLOWED ){
        //dualfprintf(fail_file,"THIS ONE4\n");
        Erf=MAX(MIN(Erf,Erf*(-Avcovorig[TT])/(ERADLIMIT+fabs(-Avcov[TT]))),Erf0);
        Erf=MAX(MIN(Erf,Erf*(-Avcov[TT])/(ERADLIMIT+fabs(-Avcovorig[TT]))),Erf0);
        //dualfprintf(fail_file,"nstep=%ld steppart=%d ijk=%d %d %d : Erforig=%g Erf=%g urfconrel=%g %g %g : Avcovorig=%g Avcov=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,Erforig,Erf,urfconrel[1],urfconrel[2],urfconrel[3],Avcovorig[0],Avcov[0]);
        //        Erf=6E-15;
      }
      else{
        //dualfprintf(fail_file,"THIS ONE5\n");
        Erf=Erf0;
      }
#if(0)
      // for RADBEAM2DKSVERT, very tricky and very sensitive (i.e. whether really fails) at coordinate singularity.
      //      if(gammarel2orig<1.0){
      //      if(gammarel2orig<0.5){
      if(gammarel2orig<=0.0){
        dualfprintf(fail_file,"THIS ONE6\n");
        Erf=Erf0;
        //        Erf=ERADLIMIT;
      }
#endif
    }
    
    //    dualfprintf(fail_file,"JONVSOLEK: usingfast=%d Avconfast: %g %g %g %g : Avconslow: %g %g %g %g : Erffast=%g Erfslow=%g urfconfast=%g %g %g urfconslow=%g %g %g\n",usingfast,Avconfast[0],Avconfast[1],Avconfast[2],Avconfast[3],Avconslow[0],Avconslow[1],Avconslow[2],Avconslow[3],Erffast,Erfslow,urfconrelfast[1],urfconrelfast[2],urfconrelfast[3],urfconrelslow[1],urfconrelslow[2],urfconrelslow[3]);
    

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
  // normally don't ever do this unless really debugging inversion.
  //  if((ptrgeom->j==14 || ptrgeom->j==13 || ptrgeom->j==15) &&nstep>=195){// || *lpflagrad!=0 && debugfail>=2){
  if(nstep>=223){
  //  if(0&&nstep>=223){// || *lpflagrad!=0 && debugfail>=2){
    // first report info so can check on inversion
    static long long int failnum=0;
    FTYPE fakedt=0.0; // since no 4-force
    FTYPE fakeCUf[4]={0}; // fake
    FTYPE dUother[NPR]={0};// fake
    struct of_state *qptr=NULL; // fake
    failnum++;
    globalpin[ENTROPY]=0.0;
    globaluu[ENTROPY]=0.0;
    pp[ENTROPY]=0.0;

    globalpin[PRAD0] = Erf;
    SLOOPA(jj) globalpin[PRAD1+jj-1] = urfconrel[jj];

    // ppfirst is faked as pp
    mathematica_report_check(0, 3, failnum, *lpflagrad, BIG, -1, -1, fakedt, ptrgeom, pp, pp, pp, globalpin, globalpin, globalpin, globaluu, globaluu, globaluu, globaluu, fakeCUf, qptr, dUother);
  }
  //  if(nstep==224) exit(0);
#endif

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
int calc_tautot_chieff(FTYPE *pp, FTYPE chieff, struct of_geom *ptrgeom, FTYPE *tautot, FTYPE *tautotmax)
{
  //xx[0] holds time
  FTYPE chi;
  chi=chieff;
  int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1}; // want to ignore non-used dimensions

  int jj;
  *tautotmax=0.0;
  SLOOPA(jj){
    tautot[jj]=chi * (dx[jj]*sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)])))*NxNOT1[jj];
    *tautotmax=MAX(*tautotmax,tautot[jj]);
  }

  return 0;
}
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

