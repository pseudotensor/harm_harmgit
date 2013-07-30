#include "decs.h"

// whether to use older method before new iter-error approaches
#define DOOLD 0

static int f_implicit_lab_old(int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f, FTYPE *fnorm);
static int get_implicit_iJ_old(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, FTYPE impepsjac, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pb, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NDIM]);
static int f_error_check_old(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, FTYPE *f1, FTYPE *f1norm, FTYPE *f3report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom);
static int Utoprimgen_failwrapper_old(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats);


//////// implicit stuff
static int koral_source_rad_implicit(int *eomtype, FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR]);
static int koral_source_rad_implicit_perdampstrategy(int dampstrategy, FTYPE imptryconv, FTYPE impallowconv, int impmaxiter, int *eomtype, FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE *radsource, FTYPE *errorabs, int *iterreturn, int *returntype);

static int f_implicit_lab(int iter, int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, int *eomtype, FTYPE *pp, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f, FTYPE *fnorm);
static int get_implicit_iJ(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtype, FTYPE impepsjac, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NPR]);
static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, int dimtypef, FTYPE *f1, FTYPE *f1norm, FTYPE *f1report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom);

static int Utoprimgen_failwrapper(int doradonly, int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats);

// debug stuff
int mathematica_report_check(int failtype, long long int failnum, int gotfirstnofail, FTYPE errorabs, int iters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);

// explicit stuff
static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gdpl, FTYPE chi, FTYPE *Gdabspl, struct of_geom *ptrgeom, FTYPE *dtsub);
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
static int Utoprimgen_failwrapper(int doradonly, int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int *eomtype, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  if(DOOLD){
    int failreturnallowable=Utoprimgen_failwrapper_old(showmessages,allowlocalfailurefixandnoreport, finalstep, evolvetype, inputtype, U, ptrgeom, pr, newtonstats);
    return(failreturnallowable);
  }


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
    u2p_rad(showmessages, allowlocalfailurefixandnoreport,U,pr,ptrgeom,lpflag,lpflagrad);
  }
  else{
    //calculating primitives  
    // OPTMARK: Should optimize this to  not try to get down to machine precision
    MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, eomtype, evolvetype, inputtype, U, ptrgeom, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);
    // this can change eomtype
  }

  // check how inversion did.  If didn't succeed, then check if soft failure and pass.  Else if hard failure have to return didn't work.
  if(IFUTOPRIMFAILSOFT(*lpflag)){
    // assume soft failure ok, but reset
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got soft MHD failure inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(IFUTOPRIMRADFAIL(*lpflagrad)){
    // can reduce Newton step if getting failure.
    // reset pflag for radiation to no failure, but treat here locally
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got some radiation inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILRAD;
  }
  else if( IFUTOPRIMFAIL(*lpflag) || IFUTOPRIMRADFAIL(*lpflagrad) ){
    // these need to get fixed-up, but can't, so return failure
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got hard failure of inversion (MHD part only considered as hard) in f_implicit_lab(): ijk=%d %d %d : %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,*lpflag,*lpflagrad);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILMHD;
  }
  else{
    // no failure
    // dualfprintf(fail_file,"No failure in Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  


  //DEBUG:
  if(debugfail>=2 && showmessages){
    struct of_state q;
    MYFUN(get_stateforcheckinversion(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    int outputtype=inputtype;
    FTYPE Unew[NPR];
    MYFUN(primtoU(outputtype,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included
    int pliter,pl;
    PLOOP(pliter,pl) dualfprintf(fail_file,"COMPARE: pl=%d pr=%g U=%g Unew=%g\n",pl,pr[pl],U[pl],Unew[pl]);
    int jj;
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d uradcon=%g uradcov=%g\n",jj,q.uradcon[jj],q.uradcov[jj]);
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d ucon=%g ucov=%g\n",jj,q.ucon[jj],q.ucov[jj]);
  }

  //DEBUG:
  if(showmessages || debugfail>=2){
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
#define QTYURAD 1 // iter or ferr
#define QTYPMHD 2 // only iter
#define QTYPRAD 3 // only iter
#define QTYENTROPYUMHD 4 // iter or ferr
#define QTYENTROPYPMHD 5 // iter (not used)

/// OLD SETUP, works even with entropy source
// V1: Works fine, but slow.  Gives hot fluid near jumps like torus edge.  Computing UU[ENTROPY] with or without source doesn't change much.
//#define IMPLICITITER (QTYURAD) // choice
//#define IMPLICITFERR (QTYURAD) // choice

/// UMHD: no setup works unless entropy source turned off, and then QTYENTROPYUMHD as FERR can't be chosen.
// V1: Ok unless entropy source added.
//#define IMPLICITITER (QTYUMHD) // choice
//#define IMPLICITFERR (QTYUMHD) // choice
// V2: Like V1 but goes unstable.
//#define IMPLICITITER (QTYUMHD) // choice
//#define IMPLICITFERR (QTYENTROPYUMHD) // choice

/// PMHD:
// V1: noisy in atmosphere and torus edge, but otherwise ok.
#define IMPLICITITER (QTYPMHD) // choice
#define IMPLICITFERR (QTYUMHD)
// V2: perfectly fine as entropy method
//#define IMPLICITITER (QTYPMHD) // choice
//#define IMPLICITFERR (QTYENTROPYUMHD) // choice

/// PRAD:
// V1: unstable with entropy source
//#define IMPLICITITER (QTYPRAD)
//#define IMPLICITFERR (QTYURAD)

/// ENTROPY:
// V1: works perfectly fine and acts like PMHD,ENTROPYUMHD method just slower:
//#define IMPLICITITER (QTYENTROPYUMHD)
//#define IMPLICITFERR (QTYENTROPYUMHD)



// same list and numbers in array for both primitives and conserved
#if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYPMHD)
#define SIGNGD7 (+1.0)
int irefU[NDIM]={UU,U1,U2,U3};
int iotherU[NDIM]={URAD0,URAD1,URAD2,URAD3};
#elif(IMPLICITITER==QTYENTROPYUMHD || IMPLICITITER==QTYENTROPYPMHD)
#define SIGNGD7 (+1.0) // not used
int irefU[NDIM]={ENTROPY,U1,U2,U3};
int iotherU[NDIM]={URAD0,URAD1,URAD2,URAD3};
#elif(IMPLICITITER==QTYURAD || IMPLICITITER==QTYPRAD)
#define SIGNGD7 (-1.0)
int irefU[NDIM]={URAD0,URAD1,URAD2,URAD3};
int iotherU[NDIM]={UU,U1,U2,U3};
#endif

// same list and numbers in array for both primitives and conserved
#if(IMPLICITFERR==QTYUMHD)
int erefU[NDIM]={UU,U1,U2,U3};
int eotherU[NDIM]={URAD0,URAD1,URAD2,URAD3};
#elif(IMPLICITFERR==QTYENTROPYUMHD)
int erefU[NDIM]={ENTROPY,U1,U2,U3};
int eotherU[NDIM]={URAD0,URAD1,URAD2,URAD3};
#elif(IMPLICITFERR==QTYURAD)
int erefU[NDIM]={URAD0,URAD1,URAD2,URAD3};
int eotherU[NDIM]={UU,U1,U2,U3};
#endif

// sign that goes into implicit differencer that's consistent with sign for SIGNGD of -1 when using the radiative uu to measure f.
#define SIGNGD2 (+1.0)
#define SIGNGD4 (+1.0) // for entropy alone for Gdpl in error function // Appears for QTYUMHD,QTYENTROPYUMHD this sign is the right one.  But both cases have lots of cold MHD inversions.
#define SIGNGD6 (-1.0) // for entropy as goes into GS from dUrad or dUmhd  //  // FUCK  -- unsure about sign!

// whether to use Ramesh's fix
// 0, 1, 2
#define RAMESHTRADTGASFIX 0

//#define SWITCHTOENTROPYIFCHANGESTOENTROPY (IMPLICITFERR==QTYUMHD ? 0 : 1)
#define SWITCHTOENTROPYIFCHANGESTOENTROPY (1)

// whether to avoid solutions where u_g<0 or rho<0 and use backup that didn't have that problem.
// KORALTODO: Bad to use backup since can be high error, but can recheck error or just use this to indicate to try entropy.
#define AVOIDSUCK 1


// need to do final check since get f1 and then do step.
#define DOFINALCHECK 1

// need better normalization for sanity check before use it again
// The current Erf,u normalization kinda checks whether Erf or u_g are going to have high *relative self* error.
#define DOSANITYCHECK 0



// below 1 if reporting cases when MAXITER reached, but allowd error so not otherwise normally reported.
#define REPORTMAXITERALLOWED 0

//FUCK: Do basic koral damping if error high or no solution.
// FUCK: entropy switch condition.
// FUCK: run normal koral tests.  Figure out all entropy signs.
// FUCK: Need backup to entropy since really dies if no backup.  Maybe backup has to be energy still.

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
static int f_implicit_lab(int iter, int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, int *eomtype, FTYPE *pp, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f, FTYPE *fnorm)
{
  if(DOOLD){
    FTYPE f4[NDIM],fnorm4[NDIM];
    int failreturn=f_implicit_lab_old(failreturnallowable, whichcall,showmessages, allowlocalfailurefixandnoreport, pp, uu0, uu, localdt, ptrgeom, f4, fnorm4); // modifies uu and pp
    int ii;
    DLOOPA(ii){
      f[erefU[ii]]=f4[ii];
      fnorm[erefU[ii]]=fnorm4[ii];
    }
    return(failreturn);
  }

  static FTYPE Tgasold=0,Tgasiter1=0;
  static FTYPE Tradold=0,Traditer1=0;



  int pliter, pl;
  int iv;
  struct of_newtonstats newtonstats;
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.
  FTYPE Gdpl[NPR],Gdabspl[NPR], Tgas;
  int failreturn;


  // default is no failure
  failreturn=0;
  // get primitive (don't change uu0).  This pp is used for inversion guess and to hold final inversion answer.
  //  PLOOP(pliter,pl) pp[pl] = pp0[pl];
  // initialize Gdpl
  PLOOP(pliter,pl) Gdpl[pl] = 0.0;




  //  PLOOP(pliter,pl) dualfprintf(fail_file,"f_implicit_lab: wc=%d i=%d j=%d pl=%d uu=%g\n",whichcall, ptrgeom->i,ptrgeom->j,pl,uu[pl]);



  if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD){ // then don't yet have updated primitives, so invert to get them
    // 1) get change in conserved quantity between fluid and radiation (equal and opposite 4-force)
    // required for inversion to get P(U) for MHD and RAD variables
    // this preserves machine accurate conservation instead of applying 4-force on each fluid and radiation separately that can accumulate errors
    FTYPE Gddt[NDIM]; DLOOPA(iv) Gddt[iv]=(uu[irefU[iv]]-uu0[irefU[iv]]);
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - Gddt[iv];
    // 2) Get estimated U[entropy]
    // mathematica has Sc = Sc0 + dt *GS with GS = -u.G/T
    FTYPE Tgaslocal=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pp[RHO],pp[UU]);
    struct of_state q; get_state(pp, ptrgeom, &q);
    FTYPE GS=0.0; DLOOPA(iv) GS += (-q.ucon[iv]*SIGNGD2*(SIGNGD7*Gddt[iv]))/(Tgaslocal+TEMPMIN); // more accurate than just using entropy from pp and ucon[TT] from state from pp.
    uu[ENTROPY] = uu0[ENTROPY] + SIGNGD6*GS; // FUCK: Problem with UMHD,UMHD no matter SIGNGD7.  Ok with URAD,URAD.   Ok with UMHD,ENTROPYUMHD if SIGNGD7 +1 and SIGNGD4 +1.
    // 3) Do MHD+RAD Inversion
    int doradonly=0; failreturn=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
    // KORALTODO: now can check if actually did eomtype==EOMGRMHD or EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function if using QTYUMHD.
    // 4) now get consistent uu[] based upon actual final primitives.
    // Needed in case inversion reduced to entropy or cold.
    // Also needed to get correct UU[ENTROPY] if did full MHD inversion or get correct U[UU] if did entropy inversion, etc.
    get_state(pp, ptrgeom, &q);
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // now have full primitives and full U including entropy and these are consistent with each other.
  }
  else if(IMPLICITITER==QTYENTROPYUMHD){
    // so iterating U[ENTROPY,U1,U2,U3]
    //    FTYPE uuorig[NPR]; PLOOP(pliter,pl) uuorig[pl]=uu[pl];
    // 1) Do pure ENTROPYMHD inversion to get pmhd
    int doradonly=0; int eomtypetemp=EOMENTROPYGRMHD; failreturn=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypetemp, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
    // set eomtype!  So for entire process for this cell, the entropy equations are used.
    *eomtype=EOMENTROPYGRMHD;
    // KORALTODO: now can check if actually did eomtypetemp==EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function in case reduced to cold MHD (maybe G=0?)
    // 2) get state (mhd and entropy are up-to-date, while rad is out-of-date, but fixed below)
    struct of_state q; get_state(pp, ptrgeom, &q);
    // 3) Compute U[UU] so have it for below (computes U[Ui,ENTROPY,RAD], but those not used or already known.)
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    //    PLOOP(pliter,pl) if(pl==U1||pl==U2||pl==U3||pl==ENTROPY) pp[pl]=pporig[pl]; // get back U[ENTROPY,Ui] so no machine error introduced from u(p(u)) for actually known quantities that won't change due to fixups or anything. -- no, primitives from Utoprimgen might change things and act as Newton step fix.
    // 4) Get actual Urad [note uses UU not irefU that is ENTROPY for irefU[0]]
    DLOOPA(iv) uu[URAD0+iv] = uu0[URAD0+iv] - (uu[UU+iv]-uu0[UU+iv]);
    // 5) Do RAD-ONLY inversion
    int failreturn2;  doradonly=1; failreturn2=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
    //  no need to concern with eomtype in RAD only case.  i.e. eomtype won't change.
    if(failreturn2>failreturn) failreturn=failreturn2; // use worst case.
    // 6) now get consistent uu[] based upon actual final primitives.
    // Needed in case inversion reduced to entropy or cold or radiative inversion used fixups.
    get_state(pp, ptrgeom, &q);
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // now have full primitives and full U including entropy and these are consistent with each other.
  }
  else if(IMPLICITITER==QTYENTROPYPMHD){
    // not setup
    myexit(92846534);
  }
  else if(IMPLICITITER==QTYPMHD){
    extern int invert_scalars(struct of_geom *ptrgeom, FTYPE *Ugeomfree, FTYPE *pr);
    // Have pmhd={ug,uvel1,uvel2,uvel3}
    // 1) get state (mhd state needed for primtoU) and Umhd, Uentropy [also computes Urad, but overwritten next and not expensive]
    struct of_state q; get_state(pp, ptrgeom, &q);
    // 1.5) trivially invert to get rho and other scalars
    pp[RHO]= uu0[RHO]/q.ucon[TT];
    invert_scalars(ptrgeom, uu,pp);
    // set p[ENTROPY] in case evaluated as a diagnostic.
    if(ENTROPY>=0) pp[ENTROPY] = pp[UU];
    // 1.6) save these primitives as original pmhd's
    FTYPE pporig[NPR]; PLOOP(pliter,pl) pporig[pl]=pp[pl];
    // 1.7) Do rest of get_state that used rho and other scalars
    extern int get_state_thermodynamics(struct of_geom *ptrgeom, FTYPE *pr, struct of_state *q);
    get_state_thermodynamics(ptrgeom, pp, &q);
    // 2) Compute Umhd and Uentropy [happens to compute old Urad from old pprad, but overwritten next step]
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // 3) Get actual Urad(G) via energy conservation (correct even if using entropy as error function, because just computed correct U[ENTROPY] consistent with U[UU].
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - (uu[irefU[iv]]-uu0[irefU[iv]]);
    // 4) Do RAD-ONLY Inversion
    int doradonly=1; failreturn=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
    //  no need to concern with eomtype in RAD only case.  i.e. eomtype won't change.
    // 5) Get consistent Urad [also computes Umhd and Uentropy, which is ok]
    get_state_uradconuradcovonly(pp, ptrgeom, &q); // only changes radiation state, not fluid state, but keeps old q's for fluid state for next step.
    // fix-up primitives to avoid violent steps in temperature
#if(RAMESHTRADTGASFIX)
    if(IMPLICITFERR==QTYUMHD){
      FTYPE Tgasnew=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pp[RHO],pp[UU]);
      FTYPE Tradnew=calc_LTE_TfromE(pp[PRAD0]);
      if(iter>1){
        if(
           (RAMESHTRADTGASFIX==1 && sign(Tgasnew-Tradnew)!=sign(Tgasold-Tradold))
           || (RAMESHTRADTGASFIX==2 && sign(Tgasnew-Tradnew)!=sign(Tgasiter1-Traditer1))
           ){
          // restrict u_g to have Tgas=Trad
          pp[UU]=calc_PEQ_ufromTrho(Tradnew,pp[RHO]);
          Tgasnew=Tradnew;
        }
      }
      if(iter==1){
        Tgasiter1=Tgasnew;
        Traditer1=Tradnew;
      }
      Tgasold=Tgasnew;
      Tradold=Tradnew;
    }
#endif
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // 6) Recover actual iterated pmhd to avoid machine related differences between original pp and pp(U(pp)) for pmhd quantities
    // This assumes that iterated pmhd is optimal and not modified except by iteration by Newton step, which is currently true.
    // This gives machine error priority to mhd primitives rather than mhd U's.
    // below not necessary since don't overwrite pp[non-rad]
    //    PLOOP(pliter,pl) if(!RADPL(pl)) pp[pl]=pporig[pl];
    // now have full primitives and full U including entropy and these are consistent with each other.
  }
  else if(IMPLICITITER==QTYPRAD){
    FTYPE pporig[NPR]; PLOOP(pliter,pl) pporig[pl]=pp[pl];
    // Have prad={Erf,uradvel1,uradvel2,uradvel3}
    // 1) get state (rad state needed for mhd_calc_rad, while old mhd state needed to estimate uu[ENTROPY] below)
    struct of_state q; get_state(pp, ptrgeom, &q);
    // 2) Compute Urad[prad0,uradcon,uradcov] [also computes old Umhd and old Uentropy, but overwritten next]
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // 3) Get actual Umhd(G) via energy conservation:
    FTYPE Gddt[NDIM]; DLOOPA(iv) Gddt[iv]=(uu[irefU[iv]]-uu0[irefU[iv]]);
    DLOOPA(iv) uu[iotherU[iv]] = uu0[iotherU[iv]] - Gddt[iv];
    // uu0[RHO] doesn't change
    // 4) Estimate U[ENTROPY](G) using old rho,u
    FTYPE Tgaslocal=compute_temp_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,pp[RHO],pp[UU]);
    // UNSURE what to do, since blows up either way:
    //    FTYPE GS=0.0; DLOOPA(iv) GS += (-q.ucon[iv]*SIGNGD2*(SIGNGD7*Gddt[iv]))/(Tgaslocal+TEMPMIN); // more accurate than just using entropy from pp and ucon[TT] from state from pp.
    FTYPE GS=0.0; DLOOPA(iv) GS += (-q.ucon[iv]*SIGNGD2*(1.0*Gddt[iv]))/(Tgaslocal+TEMPMIN); // more accurate than just using entropy from pp and ucon[TT] from state from pp.
    uu[ENTROPY] = uu0[ENTROPY] + SIGNGD6*GS;
    // 5) Invert to get pmhd (also does rad inversion, but not expensive so ok)
    int doradonly=0; failreturn=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, eomtype, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);
    // KORALTODO: now can check if actually did eomtype==EOMGRMHD or EOMENTROPYGRMHD or EOMCOLDGRMHD and apply correct error function
    // 6) Recover actual iterated prad to avoid machine related differences between original pp and pp(U(pp)) for prad quantities
    // This assumes that iterated prad is optimal and not modified except by iteration by Newton step
    PLOOP(pliter,pl) if(RADPL(pl)) pp[pl]=pporig[pl];
    // 7) Get consistent Urad [also computes Umhd and Uentropy, which is ok]
    get_state(pp, ptrgeom, &q);
    primtoU(UNOTHING,pp,&q,ptrgeom, uu);
    // now have full primitives and full U including entropy and these are consistent with each other.
  }




  ////////////////
  // get 4-force for all pl due to radiation as due to pp[uu]
  ////////////////
  koral_source_rad_calc(pp, ptrgeom, Gdpl, Gdabspl, NULL, &Tgas, NULL);




  ///////////////////////////////
  // GET ERROR FUNCTION
  ///////////////////////////////
  // compute difference vector between original and new 4-force's effect on conserved radiative quantities
  // NR1992 Eq. 16.6.16: y_{n+1} = y_n + h f(y_{n+1}) , so error function is f = (y_{n+1} - y_n) - h f(y_{n+1})
  // i.e. f->0 as change in conserved quantity approaches the updated value of 4-force
  DLOOPA(iv){
    pl=erefU[iv];
    f[pl] = (uu[pl] - uu0[pl]) + (SIGNGD2 * localdt * Gdpl[pl]);

    // get error normalization that involves actual things being differenced
    // KORALNOTE: Use Gdabspl to ensure absolute value over more terms in source in case coincidental cancellation without physical significance.
    fnorm[pl] = 0.5*(fabs(uu[pl]) + fabs(uu0[pl]) + fabs(SIGNGD2 * localdt * Gdabspl[pl]));
  }

  ///////////////////////////////
  // GET ENTROPY ERROR FUNCTION
  ///////////////////////////////
  if(IMPLICITFERR==QTYENTROPYUMHD || IMPLICITFERR==QTYUMHD && *eomtype==EOMENTROPYGRMHD && SWITCHTOENTROPYIFCHANGESTOENTROPY){ // risky perhaps.
  //  if(IMPLICITFERR==QTYENTROPYUMHD){
    // replace original equation with dS*T equation
    pl=erefU[0];
    // error function is T*dS so no actual division by T.  Found in mathematica that this works best in difficult precision cases.
    f[pl] = ((uu[pl] - uu0[pl]) + (SIGNGD4 * localdt * Gdpl[pl]))*Tgas;
    fnorm[pl] = (0.5*(fabs(uu[pl]) + fabs(uu0[pl]) + fabs(SIGNGD4 * localdt * Gdabspl[pl])))*Tgas;

    // force external inversion or any further processing to use entropy equation instead of energy equation since that's only consistent with our choice here.
    *eomtype=EOMENTROPYGRMHD;
    // assume if got here with UMHD, then each implicit step for computing f1 tries to revert back to UMHD.
  }
  else if(IMPLICITFERR==QTYUMHD){
    // force external inversion or any further processing to first try to use energy equation since that's only consistent with our choice here.
    // No, UMHD inversion call still might have reduced to entropy during 1D MHD inversion, and then should really have used entropy error function for consistency.  But that's not what one has chosen.  But, we can switch if that was the case.  This is done above.  So if here, then really sticking with EOMGRMHD.
    *eomtype=EOMGRMHD;
   
  }
  else if(IMPLICITFERR==QTYURAD){
    // URAD method forces no specific constraint on eomtype other than what came out of MHD inversion stage.
    // i.e. *eomtype=*eomtype and is can be either energy or entropy equation and that's fine.
  }


  //  dualfprintf(fail_file,"i=%d fnorm=%g %g %g %g : Gd=%g %g %g %g : uuG=%g %g %g %g\n",ptrgeom->i,fnorm[erefU[0]],fnorm[erefU[1]],fnorm[erefU[2]],fnorm[erefU[3]],Gdpl[erefU[0],Gdpl[erefU[1],Gdpl[erefU[2],Gdpl[erefU[3],SIGNGD2 * localdt * Gdpl[erefU[0],SIGNGD2 * localdt * Gdpl[erefU[1],SIGNGD2 * localdt * Gdpl[erefU[2],SIGNGD2 * localdt * Gdpl[erefU[3]);


  if(failreturn && failreturn>failreturnallowable){
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of f_implicit_lab(): %d vs. %d\n",failreturn,failreturnallowable);
    return(failreturn);
  }
  //  else{
  //    // save better guess for later inversion (including this inversion above) from this inversion
  //    PLOOP(pliter,pl) pp0[pl]=pp[pl];
  //  }



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
static int koral_source_rad_implicit(int *eomtype, FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int i1,i2,i3,iv,ii,jj,kk,pliter,sc;
  int pl;
  FTYPE realdt;
  int gotbest;
  FTYPE iJ[NPR][NPR];
  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uu[NPR],uuporig[NPR],uu0orig[NPR],bestuu[NPR];
  FTYPE pp0[NPR],ppp[NPR],pppp[NPR],pp[NPR],ppporig[NPR],pp0orig[NPR],bestpp[NPR];
  FTYPE f1[NPR],f1norm[NPR],f1report[NPR],f3report[NPR],lowestfreportf1[NPR],lowestfreportf3[NPR];
  FTYPE uubackup[NPR]={0},ppbackup[NPR]={0},fbackup[NPR];
  FTYPE radsource[NPR], deltas[NPR]; 
  extern int mathematica_report_check(int failtype, long long int failnum, int gotfirstnofail, FTYPE errorabs, int iters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);

  int eomtypelocal=*eomtype; // default choice
  // override default and selection of eomtype itself right now if error function will be entropy-based always.
  if(IMPLICITFERR==QTYENTROPYUMHD){
    // force external inversion or any further processing to use entropy equation instead of energy equation since that's only consistent with our choice here.
    eomtypelocal=*eomtype=EOMENTROPYGRMHD;
  }

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
  // FUCK
  //  showmessages=showmessagesheavy=1;



  // whether to change damp factor
#define CHANGEDAMPFACTOR 0

  //////////////
  // setup reversion to best solution for uu in case iterations lead to worse error and reach maximum iterations.
  gotbest=0;
  DLOOPA(jj) lowestfreportf1[erefU[jj]]=BIG;
  DLOOPA(jj) lowestfreportf3[erefU[jj]]=BIG;
  // setup locally-used ppfirst that can pass back as pb if good solution
  int gotfirstnofail=0;
  FTYPE ppfirst[NPR];

  ///////////////////
  // setup implicit iteration procedure and loops
  realdt = compute_dt(CUf,dt);
  FTYPE DAMPFACTOR=1.0; // factor by which step Newton's method.
  FTYPE fracdtuu0=1.0,fracdtG=1.0,fracuup=1.0; // initially try full realstep step
  FTYPE fracdtuu0p=fracdtuu0;
  FTYPE fracdtuu0pp=fracdtuu0p;

  FTYPE fracdtGp=fracdtG;
  FTYPE fracdtGpp=fracdtGp;



  /////////////////
  //
  // set uu0
  //
  ////////////////
  PLOOP(pliter,pl) uu[pl]=uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0);




#define GETDEFAULTFAILURESTATE 1 // 1 was normal before.
  // FUCK: Instead of doing this, just keep inversion error from previously as stored in pflag.  So don't reset pflag elsewhere in fixup_utoprim() for elsewhere

  // default is to allow no failure
  int failreturnallowable=UTOPRIMNOFAIL;
  int failreturnallowableuse=failreturnallowable;
  int failreturnallowablefirst=-1;

  // get default failure state
  if(GETDEFAULTFAILURESTATE){
    // Need to get default failure state.  Can allow such an error if having trouble with convergence (e.g. backing up too much)
    struct of_newtonstats newtonstats;
    // initialize counters
    newtonstats.nstroke=newtonstats.lntries=0;
    int failreturn;
    int finalstep = 1;
    FTYPE prtest[NPR];
    PLOOP(pliter,pl) prtest[pl]=piin[pl]; // initial guess (should be easy with piin=ppin(Uiin))
    int doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    failreturnallowable=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, EVOLVEUTOPRIM, UNOTHING, Uiin, ptrgeom, prtest, &newtonstats);
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
    //      PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d prtest=%21.15g Uiin=%21.15g\n",pl,prtest[pl],Uiin[pl]);
    //    }

  }
  else{
    failreturnallowablefirst=failreturnallowable;
  }


  // see if uu0->p possible and what type of failure one gets as reference for what failure to allow.
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
    int failreturn;
    int finalstep = 1;
    FTYPE prtest[NPR];
    PLOOP(pliter,pl) prtest[pl]=pb[pl]; // initial guess
    int doradonly=0;
    eomtypelocal=*eomtype; // stick with default choice so far
    failreturn=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, EVOLVEUTOPRIM, UNOTHING, uu0, ptrgeom, prtest, &newtonstats);

    // get pp0(uu0) so uu,uu0 properly associated with pp,pp0
    // set first pp that is p(uu0), assuming that can do inversion
    PLOOP(pliter,pl) ppfirst[pl]=pp[pl]=pp0[pl]=prtest[pl];

  }
  else if(USEDUINRADUPDATE==1){
    // normal behavior, where uu0 has no known primitive inversion yet.
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

  

  


  ////////////////////////////////
  // START IMPLICIT ITERATIONS
  ////////////////////////////////
  int gotbackup=0;
  int iter=0;
  int failreturn;
  int f1iter;
  int checkconv,changeotherdt;
  FTYPE impepsjac=IMPEPS;
  FTYPE errorabsf1=0.0;
  FTYPE errorabspf1=BIG;
  FTYPE errorabsf3=0.0;
  FTYPE errorabspf3=BIG;

  // initialize previous 'good inversion' based uu's
  PLOOP(pliter,pl){
    uupp[pl]=uuporig[pl]=uup[pl]=uu0orig[pl]=uu[pl];
    pppp[pl]=ppporig[pl]=ppp[pl]=pp0orig[pl]=pp[pl];
  }
  
  do{
    iter++;


    
    if(iter>10){ // KORALTODO: improve upon this later
      // assume trying hard and failing to work, then allow CASE radiation errors
      failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
      failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
    }


    //    impepsjac=1E-8;

    if(CHANGEDAMPFACTOR){
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
    }


    //vector of conserved at the previous two iterations
    PLOOP(pliter,pl)  uupp[pl]=uup[pl]; // uupp will have solution for inversion: P(uupp)
    PLOOP(pliter,pl)  uup[pl]=uu[pl]; // uup will not necessarily have P(uup) because uu used Newton step.
    PLOOP(pliter,pl)  uuporig[pl]=uu[pl];

    PLOOP(pliter,pl)  pppp[pl]=ppp[pl]; // pppp will have solution for inversion
    PLOOP(pliter,pl)  ppp[pl]=pp[pl]; // ppp will not necessarily have solution because pp used Newton step.
    PLOOP(pliter,pl)  ppporig[pl]=pp[pl];

    fracdtuu0pp=fracdtuu0p; // fracdtuu0 used when computing previous f1 and f2's
    fracdtuu0p=fracdtuu0; // fracdtuu0 used when computing previous f1 and f2's

    fracdtGpp=fracdtGp; // fracdtG used when computing previous f1 and f2's
    fracdtGp=fracdtG; // fracdtG used when computing previous f1 and f2's

    errorabspf3=errorabsf3;
    errorabspf1=errorabsf1;



    // KORALTODO: Once use certain uu0 value and succeed in getting f1, not enough.  But if take a step and *then* good f1, then know that original U was good choice and can restart from that point instead of backing up uu0 again.
    // KORALTODO: Look at whether bounding error by sign as in bisection (only true in 1D!!), and if approaching 0 slowly enough and still hit CASE, then must be real CASE not a stepping issue.
    // KORALTODO: Getting f1 is just about f1(U) as far as radiation is concerned.  So as far as CASE issues, getting f1(U) means we are good with that U and we can certainly stick with the used uu0.


    /////////////////
    //
    // get error function (f1) and inversion (uu->pp) using uu
    //
    /////////////////
    for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
      int whichcall=1;
      
      eomtypelocal=*eomtype; // re-chose default each time.  If this reduces to a new eomtype, then Jacobian will stick with that for consistency!
      failreturn=f_implicit_lab(iter, failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, pp, uu0, uu, fracdtG*realdt, ptrgeom, f1, f1norm); // modifies uu and pp
      if(failreturn){

#define BACKUPRELEASEFAIL0 (1E-4)
#define BACKUPRELEASEFAIL1 (1E-6)
#define BACKUPRELEASEFAIL2 (1E-7)
        // if backing up alot, then allow same failure as original Uiin in hopes that can recover that way (otherwise would have hoped would recover via flux update to Uiin)
        if(fracdtuu0<BACKUPRELEASEFAIL0){
          failreturnallowableuse=failreturnallowable;
        }
        if(fracdtuu0<BACKUPRELEASEFAIL1){
          failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
        }


        if(iter==1){
          // if initial uu failed, then should take smaller jump from Uiin->uu until settled between fluid and radiation.
          // If here, know original Uiin is good, so take baby steps in case uu needs heavy raditive changes.
          // if f1 fails, try going back to Uiin a bit
          fracdtuu0*=RADDAMPDELTA; // DAMP Uiin->uu0 step that may be too large and generated too large G

          //          if(fracdtuu0<BACKUPRELEASEFAIL2) fracdtuu0=0.0; // just stop trying to start with some uu0 and revert to Uiin

          // modifies uup (ppp) and uu (pp) as if starting over, and then another call to f_implicit_lab(f1) will change uu by generally smaller amount.
          PLOOP(pliter,pl){
            uup[pl]=uu[pl]=uu0[pl];
            ppp[pl]=pp[pl]=pp0[pl];
          }
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
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f1 failed: iter=%d  Backing-up both uu0 and G.: f1iter=%d fracdtuu0=%g fracdtG=%g fracuup=%g\n",iter,f1iter,fracdtuu0,fracdtG,fracuup);
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Ui=%26.20g uu0=%26.20g uu0orig=%26.20g uu=%26.20g uup=%26.20g dUother=%26.20g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
      }// end if failed inversion in f_implicit_lab()
      else{
        // then success, so was able to do inversion P(U) with change in radiation on fluid: P(uu[fluid] = uu0[fluid] - (uu[rad]-uu0[rad]))
        // This doesn't necessarily mean could do P(uu0) except for iter=1.
        // Corresponds to success for a certain P(uu0,uu) pair.
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"SUCCESS: pl=%d Ui=%26.20g uu0=%26.20g uu0orig=%26.20g uu=%26.20g uup=%26.20g dUother=%26.20g: pp(pnew)=%g : fracdtuu0=%g fracdtG=%g fracuup=%g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl],pp[pl],fracdtuu0,fracdtG,fracuup);
        break;
      }
    }// end loop over f1iter
    if(f1iter==MAXF1TRIES){
      if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES: fracdtuu0=%g nstep=%ld steppart=%d ijk=%d %d %d : iter=%d\n",fracdtuu0,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter);
      failnum++;
      mathematica_report_check(0, failnum, gotfirstnofail, BIG, iter, realdt, ptrgeom, ppfirst,pp,pb,piin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
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
    if(showmessagesheavy) dualfprintf(fail_file,"i=%d f1: %g %g %g %g\n",ptrgeom->i,f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]]);



    ///////////////////////
    //
    // set backup solution
    //
    // avoid arbitrary sucking on thermal energy densities
    // but only get backup if already below error
    // at this point, uu and pp are fully up-to-date and consistent
    //
    ///////////////////////
    if(pp[RHO]>=0.0 && pp[UU]>=0.0 && pp[PRAD0]>=0.0){
      gotbackup=1;
      PLOOP(pliter,pl){
        fbackup[pl]=f1[pl];
        uubackup[pl]=uu[pl];
        ppbackup[pl]=pp[pl];
      }
    }

    /////////////////
    //
    // see if pre-convergence (happens if force is small or no force at all.  Can't necessarily continue since Jacobian can require arbitrarily large dU on fluid and fail to invert even if no fluid-radiation interaction!
    //test pre-convergence using initial |dU/U|
    // KORALTODO: This isn't a completely general error check since force might be large for fluid that needs itself to have more accuracy, but if using ~NUMEPSILON, won't resolve 4-force of radiation on fluid to better than that.
    //
    /////////////////
    FTYPE LOCALPREIMPCONV=(10.0*NUMEPSILON); // more strict than later tolerance
    int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive R^{ti} type
    if(f_error_check(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,realdt,dimtypef,f1,f1norm,f1report,Uiin, uu0,uu,ptrgeom)){
      break;
    }

    //////////////
    //
    // get error using f1 and f1norm
    //
    //////////////
    int convreturnf1;
    convreturnf1=f_error_check(showmessages, showmessagesheavy, iter, IMPTRYCONV,realdt, dimtypef ,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom);
    errorabsf1=0.0;     DLOOPA(jj) errorabsf1     += fabs(f1report[erefU[jj]]);
    if(convreturnf1){
      // but don't break, since need to iterate a bit first and check |dU/U|
    }

    // store error and solution in case eventually lead to max iterations and actually get worse error
    // f1 based
    FTYPE errorabsbest=0.0;
    DLOOPA(jj) errorabsbest += fabs(lowestfreportf1[erefU[jj]]);
    if(errorabsbest>errorabsf1 && isfinitel(errorabsf1)){
      PLOOP(pliter,pl) bestuu[pl]=uu[pl];
      PLOOP(pliter,pl) bestpp[pl]=pp[pl];
      DLOOPA(jj) lowestfreportf1[erefU[jj]]=f1report[erefU[jj]];
      gotbest=1;
    }



    int notfinite;
    if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
      notfinite = !isfinitel(uu[irefU[0]])|| !isfinitel(uu[irefU[1]])|| !isfinitel(uu[irefU[2]])|| !isfinitel(uu[irefU[3]]) || !isfinitel(uup[irefU[0]])|| !isfinitel(uup[irefU[1]])|| !isfinitel(uup[irefU[2]])|| !isfinitel(uup[irefU[3]]);
    }
    else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
      notfinite = !isfinitel(pp[irefU[0]])|| !isfinitel(pp[irefU[1]])|| !isfinitel(pp[irefU[2]])|| !isfinitel(pp[irefU[3]]) || !isfinitel(ppp[irefU[0]])|| !isfinitel(ppp[irefU[1]])|| !isfinitel(ppp[irefU[2]])|| !isfinitel(ppp[irefU[3]]);
    }

    if(!notfinite){
    
      /////////
      //
      // get Jacobian and inversion Jacobian 
      //
      /////////
      //      eomtypelocal=*eomtype; // re-chose default each time.  No, stick with what f1 reduced to for consistency.
      int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, impepsjac, uu, uup, uu0, pp, ppp, fracdtG, realdt, ptrgeom, f1, f1norm, iJ);


      if(failreturniJ!=0) return(failreturniJ);

      if(showmessagesheavy){
        int iii,jjj;
        DLOOP(iii,jjj) dualfprintf(fail_file,"iJ[i %d][e %d]=%g\n",iii,jjj,iJ[irefU[iii]][erefU[jjj]]);
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
      if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
        PLOOP(pliter,pl) uu[pl] = uup[pl];
        DLOOP(ii,jj) uu[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];

        if(showmessagesheavy) dualfprintf(fail_file,"POSTDX: uu: %g %g %g %g : uup=%g %g %g %g\n",uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],uup[irefU[0]],uup[irefU[1]],uup[irefU[2]],uup[irefU[3]]);
      }
      else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
        PLOOP(pliter,pl) pp[pl]=ppp[pl];
        DLOOP(ii,jj) pp[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];

        if(showmessagesheavy) dualfprintf(fail_file,"POSTDX: pp: %g %g %g %g : ppp=%g %g %g %g\n",pp[irefU[0]],pp[irefU[1]],pp[irefU[2]],pp[irefU[3]],ppp[irefU[0]],ppp[irefU[1]],ppp[irefU[2]],ppp[irefU[3]]);

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
        // try increasing amount of uu or pp used
        fracuup*=RADDAMPUNDELTA;
        checkconv=0;
        changeotherdt=0; // ensure fracuup back to 1.0 first before reverting others.
      }
      //    if(fracdtuu0!=1.0){
      if(fabs(fracdtuu0-1.0)>10.0*NUMEPSILON && changeotherdt){
        // try increasing uu0 away from Uiin to account for full dUother
        fracdtuu0*=RADDAMPUNDELTA;
        PLOOP(pliter,pl) uu0[pl]=UFSET(CUf,fracdtuu0*dt,Uiin[pl],Ufin[pl],dUother[pl],0.0); // modifies uu0
        // KORALNOTE: No need to get pp0, since never used.  uu0 only used in error function.
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
    // test convergence using |dU/U|
    // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
    // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
    /////////
    FTYPE f3[NPR],f3norm[NPR];
    // 0 = conserved R^t_\nu type, 1 = primitive R^{ti} type
    if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){ // still considering iterate error
      dimtypef=DIMTYPEFCONS;
      DLOOPA(ii){
        f3[erefU[ii]]=(uu[irefU[ii]]-uup[irefU[ii]]);
        f3norm[erefU[ii]]=fabs(uu[irefU[ii]])+fabs(uup[irefU[ii]]);
      }
    }
    else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){ // still considering iterate error
      dimtypef=DIMTYPEFPRIM;
      DLOOPA(ii){
        f3[erefU[ii]]=(pp[irefU[ii]]-ppp[irefU[ii]]);
        f3norm[erefU[ii]]=fabs(pp[irefU[ii]])+fabs(ppp[irefU[ii]]);
      }
    }
  
    // store error and solution in case eventually lead to max iterations and actually get worse error
    // f_error_check(uu0,uu) is ok to use since it just normalizes error
    int convreturnf3;
    convreturnf3=f_error_check(showmessages, showmessagesheavy, iter, IMPTRYCONV,realdt, dimtypef ,f3,f3norm,f3report,Uiin,uu0,uu,ptrgeom);
    errorabsf3=0.0;     DLOOPA(jj) errorabsf3     += fabs(f3report[erefU[jj]]);
    // while using f1 for true error, can't do better if f3 error is below near machine precision.
    int convreturnf3limit;
    convreturnf3limit=f_error_check(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,realdt,DIMTYPEFCONS,f1,f1norm,f1report,Uiin, uu0,uu,ptrgeom);





    if(CHANGEDAMPFACTOR){
      /////////
      //
      // DAMP if error changes too much
      //
      /////////
      if(iter==1){
        DAMPFACTOR=0.7;
      }
      else{
        if(errorabsf3>errorabspf3 && DAMPFACTOR>0.1) DAMPFACTOR*=0.5;
      }
    }

    /////////////////
    //
    // check convergence
    //
    /////////////////
    if(checkconv){
      // then can check convergence: using f1 and f3limit (not f3 so far anymore)
      if(convreturnf1 || convreturnf3limit){
         // then done
        break;
      }
      else{
        // then not done
      }

    } // end if checking convergence




    /////////
    // see if took too many Newton steps or not finite results
    /////////
    int itermaxed=iter>IMPMAXITER;
    int convreturnf1allow;
    if(itermaxed || notfinite ){
      // f1-based
      convreturnf1allow=f_error_check(showmessages, showmessagesheavy, iter, IMPALLOWCONV,realdt,dimtypef,f1,f1norm,f1report,Uiin,uup,uu,ptrgeom);
      
      if(gotbest){
        // f1 based
        // see if should revert to prior best
        errorabsbest=0.0;
        DLOOPA(jj) errorabsbest += fabs(lowestfreportf1[erefU[jj]]);
        if(errorabsbest<errorabsf1 || !isfinitel(errorabsf1) ){
          PLOOP(pliter,pl) uu[pl]=bestuu[pl];
          PLOOP(pliter,pl) pp[pl]=bestpp[pl];
          errorabsf1=errorabsbest;
        }
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"Using best: %g %g\n",errorabsf1,errorabsbest);
        // get new convreturnf1allowlow
        convreturnf1allow=(lowestfreportf1[erefU[0]]<IMPALLOWCONV && lowestfreportf1[erefU[1]]<IMPALLOWCONV && lowestfreportf1[erefU[2]]<IMPALLOWCONV && lowestfreportf1[erefU[3]]<IMPALLOWCONV);
      }

      // KORALTODO: If convreturnf1allow doesn't work, but still (say) 10% error, might want to hold onto result in case explicit backup fails as well (which is likely), in which case *much* better to use 10% error because otherwise 4-force not accounted for, which can lead to very big changes in fluid behavior due to large flux from previous step.
      // KORALTODO: Or, perhaps should really take it as a failure and use fixups.  Probably should allow for result to be written if error<10%, but only use as super-backup in fixups.  So should set pflag still.

      if(convreturnf1allow){
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab().  But f1 was allowed error. checkconv=%d (if checkconv=0, could be issue!) : %g %g %g %g : %g %g %g %g : errorabs=%g : %g %g %g\n",IMPMAXITER,checkconv,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],lowestfreportf1[erefU[0]],lowestfreportf1[erefU[1]],lowestfreportf1[erefU[2]],lowestfreportf1[erefU[3]],errorabsf1,fracdtuu0,fracuup,fracdtG);
        if(REPORTMAXITERALLOWED){
          failnum++;
          mathematica_report_check(6, failnum, gotfirstnofail, errorabsf1, iter, realdt, ptrgeom, ppfirst,pp,pb,piin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
        }

        break;
      }
      else{
        // KORALTODO: Need backup that won't fail.
        if(debugfail>=2){
          if(itermaxed) dualfprintf(fail_file,"iter>IMPMAXITER=%d : iter exceeded in solve_implicit_lab(). nstep=%ld steppart=%d ijk=%d %d %d :  Bad error.\n",IMPMAXITER,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          if(notfinite) dualfprintf(fail_file,"IMPGOTNAN at iter=%d : in solve_implicit_lab(). ijk=%d %d %d :  Bad error.\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          dualfprintf(fail_file,"checkconv=%d failreturnallowable=%d: %g %g %g %g : %g %g %g %g\n",checkconv,failreturnallowable,f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],lowestfreportf1[erefU[0]],lowestfreportf1[erefU[1]],lowestfreportf1[erefU[2]],lowestfreportf1[erefU[3]]);
          if(1||showmessages){
            failnum++;
            mathematica_report_check(1, failnum, gotfirstnofail, errorabsf1, iter, realdt, ptrgeom, ppfirst,pp,pb,piin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
            if(doingit==1) myexit(100000000); // DEBUG
          }
        }        
        return(1);
      }
    }//end if maximum iterations

  }// end do
  while(1);




  if(AVOIDSUCK){
    ////////////////////
    //
    // See if final solution is sucking on energy.  If so, backup to solution where that wasn't happening.
    //
    ////////////////////
    if(gotbackup && !(pp[RHO]>=0.0 && pp[UU]>=0.0 && pp[PRAD0]>=0.0)){
      // then might be sucking into radiation or into thermal energy density nothingness, so see if backup is good enough error
      if(debugfail>=2) dualfprintf(fail_file,"Thought was good error, but gotbackup=%d and pp=%g %g %g : f1report: %g %g %g %g : nstep=%ld steppart=%d ijk=%d %d %d\n",gotbackup,pp[RHO],pp[UU],pp[PRAD0],f1report[erefU[0]],f1report[erefU[1]],f1report[erefU[2]],f1report[erefU[3]],nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
      PLOOP(pliter,pl){
        uu[pl]=uubackup[pl];
        pp[pl]=ppbackup[pl];
      }
      if(debugfail>=2) dualfprintf(fail_file,"Used backup: pp=%g %g %g : fbackup=%g %g %g %g\n",pp[RHO],pp[UU],pp[PRAD0],fbackup[erefU[0]],fbackup[erefU[1]],fbackup[erefU[2]],fbackup[erefU[3]]);
    }
  }


  if(DOSANITYCHECK||DOFINALCHECK){
    //////////////////////////
    //
    // sanity check on error, and final inversion
    //
    ////////////////////////
    int whichcall=1;
    //  eomtypelocal=*eomtype; // re-chose default each time. No, stick with what f1 (last call to f1) chose
    failreturn=f_implicit_lab(iter,failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, pp, uu0, uu, fracdtG*realdt, ptrgeom, f1, f1norm); // modifies uu and pp
    if(DOSANITYCHECK){
      // can't just check (uup-uu)/(|uu|-|uup|) because uup and uu  might diverge, leading to huge G, while uup-uu might be similar enough.
      // In that case, also using f1 and f1norm not good enough.
      // Have to compare with static norm
      // overwrite f1norm with static choice to avoid run-away G
      jj=TT; f1norm[erefU[jj]]=MAX(pb[PRAD0],pb[UU]);
      SLOOPA(jj) f1norm[erefU[jj]]=MAX(pb[PRAD0],pb[UU])/sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)]));  // make into conserved type: R^t_\nu ~ {prad0,uu} / sqrt(gcon[nu,nu])
    }
    int dimtypef=DIMTYPEFCONS; // 0 = conserved R^t_\nu type, 1 = primitive (u,v^i) type, i.e. v^i has no energy density term
    int convreturn=f_error_check(showmessages, showmessagesheavy, iter, IMPALLOWCONV,realdt,dimtypef,f1,f1norm,f1report,Uiin,uu0,uu,ptrgeom);
    errorabsf1=0.0;     DLOOPA(jj) errorabsf1     += fabs(f1report[erefU[jj]]);
    if(convreturn==0){
      if(debugfail>=2){
        dualfprintf(fail_file,"Thought was good error (gotbackup=%d : %g %g %g), but final-check or re-check showed otherwise: nstep=%ld steppart=%d ijk=%d %d %d\n",gotbackup,pp[RHO],pp[UU],pp[PRAD0],nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        DLOOPA(jj) dualfprintf(fail_file,"Issue with jj=%d %g %g %g : %g %g\n",jj,f1[erefU[jj]],f1norm[erefU[jj]],f1report[erefU[jj]],uup[irefU[jj]],uu[irefU[jj]]);
        failnum++;
        mathematica_report_check(4, failnum, gotfirstnofail, errorabsf1, iter, realdt, ptrgeom, ppfirst,pp,pb,piin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
      }
      return(1);
    }
  }

  /////////////
  //
  // diagnose
  //
  ////////////
  numofiter+=iter;
  // i=j=k=0 just to show infrequently
  if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0) || showmessagesheavy) dualfprintf(fail_file,"numimplicits=%lld averagef1iter=%g averageiter=%g\n",numimplicits,(FTYPE)numoff1iter/(FTYPE)numimplicits,(FTYPE)numofiter/(FTYPE)numimplicits);




  ///////////////////
  //
  // get source update as "dU" = dU/dt using real dt that used during implicit iterations, and will eventually use to update U in advance.c.
  // apply source update as force
  // KORALNOTE: As long as f_implicit_lab() updates both full primitive and full U, using uu below is fine and good.
  PLOOP(pliter,pl) radsource[pl] = +(uu[pl]-uu0[pl])/realdt;
  // KORALNOTE: Could re-enforce energy conservation here, but would be inconsistenet with how applied error function.

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
  // pp was modified by f_implicit_lab(f1,pp) with output from inversion returned through pp0
  // only use pp if successful with implicit method, since if not successful can be various bad reasons with no good pb
  //
  ////////////////
  PLOOP(pliter,pl) pb[pl]=pp[pl];

  ///////////
  //
  // choose new eomtype for external inverison or other checks
  //
  // This ensures whatever implicit solver settled on using, the external inversion doesn't switch.
  //
  //////////
  *eomtype=eomtypelocal;

  // can further just choose to avoid inversion entirely since this step fully inverted (even if somewhat inaccurately) all parts used in advance.c:
  // especially for entropy inversion, no need to be very accurate as normal inversion does, since no need to have exact entropy conservation, unlike desirable to have exact energy conservation.
  //  if(IMPLICITITER==QTYPMHD && IMPLICITFERR==QTYENTROPYUMHD){
  // more general case when can just avoid external entropy inversion
  if(*eomtype==EOMENTROPYGRMHD && SWITCHTOENTROPYIFCHANGESTOENTROPY || IMPLICITFERR==QTYENTROPYUMHD){ // risky perhaps
    *eomtype=EOMDONOTHING;
  }


  ////////////////////
  //
  // Do some diagnostics and reporting
  //
  //////////////////////
  if(debugfail>=2){
    errorabsf1=0.0;
    if(gotbest) DLOOPA(jj) errorabsf1 += fabs(lowestfreportf1[erefU[jj]]);
    else DLOOPA(jj) errorabsf1 += fabs(f1report[erefU[jj]]);
    //    dualfprintf(fail_file,"errorabs=%g\n",errorabsf1);
    numhisterr[MAX(MIN((int)(-log10l(errorabsf1)),NUMNUMHIST-1),0)]++;
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

  //  if(nstep==4 && steppart==0){
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"PRIMIMPLICIT(%d,%g): pl=%d pp=%21.15g uu=%21.15g : Uiin=%21.15g Ufin=%21.15g uu0=%21.15g : dUother=%21.15g : dU=%21.15g\n",*eomtype,realdt,pl,pp[pl],uu[pl],Uiin[pl],Ufin[pl],uu0[pl],dUother[pl]*realdt,dUcomp[sc][pl]*realdt);
  //  }


  return(0);
  
}

















#define RETURNIMPLICITFAIL 0
#define RETURNIMPLICITTRY 1
#define RETURNIMPLICITALLOW 2
#define RETURNIMPLICITFAILALLOW 3

// wrapper for koral_source_rad_implicit_perdampstrategy()
// CHANGED: pb (if successful), dUcomp[RADSOURCE] (if successful)
static int koral_source_rad_implicit_new(int *eomtype, FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE (*dUcomp)[NPR])
{
  int failreturn;
#define NUMDAMPSTRATEGY (12)
  // 1: works for RADPULSEPLANAR  2: works for RADDONUT // 1?: for reduced DAMPFACTOR each ?.
  // start with simpler #1, try #2, then just reduce DAMPFACTOR sequentially like in koral code.
  int dampstrategyarray[NUMDAMPSTRATEGY]={2,0, 1,13,14,15,16,17,18,19,20,21}; // if add/remove, change number
  int dampstrategy;
  FTYPE imptryconvarray[NUMDAMPSTRATEGY];//={IMPTRYCONV,IMPTRYCONV2, IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2,IMPTRYCONV2}; // ""
  FTYPE imptryconv;
  FTYPE impallowconvarray[NUMDAMPSTRATEGY];//={IMPALLOWCONV,IMPALLOWCONV2, IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2,IMPALLOWCONV2}; // ""
  FTYPE impallowconv;
  int impmaxiterarray[NUMDAMPSTRATEGY];//={IMPMAXITER,IMPMAXITER2, IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2,IMPMAXITER2}; // ""
  int impmaxiter;
  FTYPE errorabsbest=BIG;
  FTYPE radsource[NPR]={0},radsourcebest[NPR]={0};
  FTYPE pbsource[NPR]={0},pbsourcebest[NPR]={0};
  int sc;
  int pliter,pl;
  int returntype;
  int iterreturn,itertotal;
  FTYPE errorabs;
  int strati;
  int gotsomesolution=0;



  int showmessages=0;
  int showmessagesheavy=0;


  // get defaults
  itertotal=0;
  for(strati=0;strati<NUMDAMPSTRATEGY;strati++){
    imptryconvarray[strati]=IMPTRYCONV2;
    impallowconvarray[strati]=IMPALLOWCONV2;
    impmaxiterarray[strati]=IMPMAXITER2;
  }

  // get tau to know how hard should try for implicit solver.  The higher the tau, the harder it is, but the more likely a problem if innacurate solution (even accurate to 1E-10 is not good enough!)
  FTYPE tautot[NDIM],tautotmax;
  calc_tautot(pb, ptrgeom, tautot, &tautotmax);
#define TAUIMPLIMIT (0.05)
  if(tautotmax>=TAUIMPLIMIT){
    imptryconvarray[0]=IMPTRYCONVHIGHTAU;
    impallowconvarray[0]=IMPALLOWCONV;
    impmaxiterarray[0]=IMPMAXITER;
  }
  else{
    imptryconvarray[0]=IMPTRYCONV;
    impallowconvarray[0]=IMPALLOWCONV;
    impmaxiterarray[0]=IMPMAXITER;
  }
  


    
#if(1)
  impmaxiter=IMPMAXITER2;
  //    failreturn=koral_source_rad_implicit_perdampstrategy(1, IMPTRYCONV, IMPALLOWCONV, eomtype, pb, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother, radsource, &errorabs, &iterreturn, &returntype);
  failreturn=koral_source_rad_implicit_perdampstrategy(2, IMPTRYCONV, IMPALLOWCONV, impmaxiter, eomtype, pb, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother, radsource, &errorabs, &iterreturn, &returntype);
  if(failreturn==0){
    sc = RADSOURCE;
    PLOOP(pliter,pl){
      // store source update in dUcomp for return.
      dUcomp[sc][pl] += radsource[pl];
    }
  }
  itertotal+=iterreturn;

#elif(0)
  for(strati=0;strati<NUMDAMPSTRATEGY;strati++){
    dampstrategy=dampstrategyarray[strati];
    imptryconv=imptryconvarray[strati];
    impallowconv=impallowconvarray[strati];
    impmaxiter=impmaxiterarray[strati];

    returntype=RETURNIMPLICITFAIL;//default
    PLOOP(pliter,pl) pbsource[pl]=pb[pl]; // reset to defalut input
    failreturn=koral_source_rad_implicit_perdampstrategy(3, imptryconv, impallowconv, impmaxiter, eomtype,pbsource, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother, radsource, &errorabs, &iterreturn, &returntype);
    itertotal+=iterreturn;

    returntype=RETURNIMPLICITFAIL;//default
    PLOOP(pliter,pl) pbsource[pl]=pb[pl]; // reset to defalut input
    failreturn=koral_source_rad_implicit_perdampstrategy(1, imptryconv, impallowconv, impmaxiter, eomtype, pbsource, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother, radsource, &errorabs, &iterreturn, &returntype);
    itertotal+=iterreturn;

    if(failreturn==0){
      gotsomesolution=1;
      // get best
      if(errorabs<errorabsbest){
        errorabsbest=errorabs;
        PLOOP(pliter,pl){
          radsourcebest[pl]=radsource[pl];
          pbsourcebest[pl]=pbsource[pl];
        }
      }
     
      dualfprintf(fail_file,"nstep=%ld steppart=%d ijk=%d %d %d strati=%d errorabs=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,strati,errorabs);
 
      // if no failure and successful at getting desired error, then break.  Otherwise continue.
      //      if(strati==0 && returntype==RETURNIMPLICITTRY || strati>0 && (returntype==RETURNIMPLICITALLOW || returntype==RETURNIMPLICITTRY)) break;
      if(strati==0) break;
    }// end if no failure in implicit solver
  }


  if(gotsomesolution){
    // use best got
    errorabs=errorabsbest;
    sc = RADSOURCE;
    PLOOP(pliter,pl){
      // store source update in dUcomp for return.
      dUcomp[sc][pl] += radsourcebest[pl];
      pb[pl]=pbsourcebest[pl];
    }
  }
  else{
    if(debugfail>=2) dualfprintf(fail_file,"Got no solution.\n");
  }


#else
  for(strati=0;strati<NUMDAMPSTRATEGY;strati++){
    dampstrategy=dampstrategyarray[strati];
    imptryconv=imptryconvarray[strati];
    impallowconv=impallowconvarray[strati];
    impmaxiter=impmaxiterarray[strati];

    returntype=RETURNIMPLICITFAIL;//default
    PLOOP(pliter,pl) pbsource[pl]=pb[pl]; // reset to defalut input
    failreturn=koral_source_rad_implicit_perdampstrategy(dampstrategy, imptryconv, impallowconv, impmaxiter, eomtype, pbsource, piin, Uiin, Ufin, CUf, ptrgeom, q, dUother, radsource, &errorabs, &iterreturn, &returntype);
    itertotal+=iterreturn;

    if(failreturn==0){
      // get best
      if(errorabs<errorabsbest){
        errorabsbest=errorabs;
        PLOOP(pliter,pl){
          radsourcebest[pl]=radsource[pl];
          pbsourcebest[pl]=pbsource[pl];
        }
        if(strati>0) dualfprintf(fail_file,"new best: nstep=%ld steppart=%d ijk=%d %d %d strati=%d errorabs=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,strati,errorabs);
      }
     
 
      // if no failure and successful at getting desired error, then break.  Otherwise continue.
      if(strati==0 && returntype==RETURNIMPLICITTRY || strati>0 && (returntype==RETURNIMPLICITALLOW || returntype==RETURNIMPLICITTRY)) break;
      //      if(strati==1) break;
    }// end if no failure in implicit solver
  }

  // use best got
  errorabs=errorabsbest;
  sc = RADSOURCE;
  PLOOP(pliter,pl){
    // store source update in dUcomp for return.
    dUcomp[sc][pl] += radsourcebest[pl];
    pb[pl]=pbsourcebest[pl];
  }

#endif


  // Do some diagnostics and reporting
  if(debugfail>=2){

    // report result of best answer
    if(showmessagesheavy&&errorabs>IMPALLOWCONV) dualfprintf(fail_file,"Never got IMPALLOWCONV: %g %g\n",IMPALLOWCONV,errorabs);
    if(showmessagesheavy&&errorabs>IMPALLOWCONV2) dualfprintf(fail_file,"Never got IMPALLOWCONV2: %g %g\n",IMPALLOWCONV,errorabs);

    if(showmessagesheavy&&errorabs>IMPTRYCONV) dualfprintf(fail_file,"Never got IMPTRYCONV: %g %g\n",IMPTRYCONV,errorabs);

#define NUMNUMHIST (20)
    static long long int numhisterr[NUMNUMHIST]={0}; // histogram of error for implicit solver to be reported infrequently
    static long long int numhistiter[MAX(IMPMAXITER,IMPMAXITER2)+1]={0}; // histogram of error for implicit solver to be reported infrequently
    //    dualfprintf(fail_file,"errorabs=%g\n",errorabs);
    numhisterr[MAX(MIN((int)(-log10l(errorabs+SMALL)),NUMNUMHIST-1),0)]++;
    numhistiter[MAX(MIN(itertotal,impmaxiter),0)]++;
#define HISTREPORTSTEP (20)
    if(nstep%HISTREPORTSTEP==0 && ptrgeom->i==0 && ptrgeom->j==0 && ptrgeom->k==0){
      int histi;
      for(histi=0;histi<NUMNUMHIST;histi++){
        dualfprintf(fail_file,"numhisterr%d=%lld\n",histi,numhisterr[histi]);
      }
      for(histi=0;histi<=MAX(IMPMAXITER,IMPMAXITER2);histi++){
        dualfprintf(fail_file,"numhistiter%d=%lld\n",histi,numhistiter[histi]);
      }
    }
  }



  return(failreturn);

}




// compute changes to U (both T and R) using implicit method
// KORALTODO: If doing implicit, should also add geometry source term that can sometimes be stiff.  Would require inverting sparse 8x8 matrix (or maybe 6x6 since only r-\theta for SPC).  Could be important for very dynamic radiative flows.
static int koral_source_rad_implicit_perdampstrategy(int dampstrategy, FTYPE imptryconv, FTYPE impallowconv, int impmaxiter, int *eomtype, FTYPE *pb, FTYPE *piin, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q, FTYPE *dUother ,FTYPE *radsource, FTYPE *errorabs, int *iterreturn, int *returntype)
{
  int i1,i2,i3,iv,ii,jj,kk,pliter,sc;
  int pl;
  FTYPE realdt;
  FTYPE iJ[NPR][NPR];
  FTYPE uu0[NPR],uup[NPR],uupp[NPR],uu[NPR],uuporig[NPR],uu0orig[NPR],bestuu[NPR];
  FTYPE pp0[NPR],ppp[NPR],pppp[NPR],pp[NPR],ppporig[NPR],pp0orig[NPR],bestpp[NPR];
  FTYPE f1[NPR],f1norm[NPR],f3report[NPR],lowestfreport[NPR];
  FTYPE uubackup[NPR]={0},ppbackup[NPR]={0};
  FTYPE deltas[NPR]; 
  extern int mathematica_report_check(int failtype, long long int failnum, int gotfirstnofail, FTYPE errorabs, int iters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother);

  FTYPE f3reportuu[NPR];
  int gotbest=0,iterbest=0;

  int eomtypelocal=*eomtype; // default
  // override default and selection of eomtype itself right now if error function will be entropy-based always.
  if(IMPLICITFERR==QTYENTROPYUMHD){
    // force external inversion or any further processing to use entropy equation instead of energy equation since that's only consistent with our choice here.
    eomtypelocal=*eomtype=EOMENTROPYGRMHD;
  }

  // static counter for diagnosing issues
  static long long int numimplicits=0;
  static long long int numoff1iter=0,numofiter=0;
  static long long int failnum=0;
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
  DLOOPA(jj) lowestfreport[erefU[jj]]=BIG;
  // setup locally-used pp that can pass back as pb if good solution
  int gotfirstnofail=0;
  FTYPE ppfirst[NPR];

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
    PLOOP(pliter,pl) prtest[pl]=pb[pl]; // initial guess
    int doradonly=0;
    eomtypelocal=*eomtype; // default
    failreturnallowable=Utoprimgen_failwrapper(doradonly,showmessages,allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, EVOLVEUTOPRIM, UNOTHING, Uiin, ptrgeom, prtest, &newtonstats);
    if(failreturnallowable!=UTOPRIMGENWRAPPERRETURNNOFAIL){
      if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is already a problem with %d\n",failreturnallowable);
    }
    else{
      if(showmessagesheavy && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() says that Uiin is NOT a problem with %d\n",failreturnallowable);
    }

    // set first pp that is p(uu0), assuming that can do inversion
    PLOOP(pliter,pl) pp0[pl]=prtest[pl]; // gets correct pp0
    PLOOP(pliter,pl) ppfirst[pl]=prtest[pl];
    if(failreturnallowable==UTOPRIMGENWRAPPERRETURNNOFAIL) gotfirstnofail=1;
    else gotfirstnofail=0;

    failreturnallowablefirst=failreturnallowable;

  }
  else{
    // then (not recommended) just using Uiin as uu0
    PLOOP(pliter,pl){
      uu[pl] = uu0[pl] = Uiin[pl];
      pp[pl] = pp0[pl] = piin[pl];
      ppfirst[pl] = piin[pl];
    }
  }  

  

  


  ////////////////////////////////
  // START IMPLICIT ITERATIONS
  ////////////////////////////////
  int iter=0;
  int failreturn;
  int f1iter;
  int checkconv,changeotherdt;
  FTYPE impepsjac=IMPEPS;
  FTYPE errorabslocal=BIG;
  FTYPE errorabsp=errorabslocal;
  FTYPE errorabsbest=BIG;
  int convreturn;
  FTYPE f3[NPR],f3norm[NPR];
  int f1break;
  FTYPE LOCALPREIMPCONV=(NUMEPSILON); // more strict than later tolerance

  // initialize previous 'good inversion' based uu's
  PLOOP(pliter,pl)  uupp[pl]=uuporig[pl]=uup[pl]=uu0orig[pl]=uu[pl];
  PLOOP(pliter,pl)  pppp[pl]=ppporig[pl]=ppp[pl]=pp0orig[pl]=pp[pl];
  
  do{
    iter++;

    
    if(iter>10){ // KORALTODO: improve upon this later
      // assume trying hard and failing to work, then allow CASE radiation errors
      failreturnallowableuse=UTOPRIMGENWRAPPERRETURNFAILRAD;
      failreturnallowable=UTOPRIMGENWRAPPERRETURNFAILRAD;
    }

    if(dampstrategy==2){
      // cautious first step to get reasonable error measurement. Allows often only 1 iteration to get sufficiently small error, while DAMPFACTOR=1 would already go beyond point where error actually increases.
      if(iter==1){
        DAMPFACTOR=0.37;
      }
    }

    if(dampstrategy>=10){
      DAMPFACTOR=pow(2.0,9.0-(FTYPE)dampstrategy);
    }


    if(dampstrategy==1||dampstrategy==2){
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
    }

    if(dampstrategy==3){
      DAMPFACTOR=1.0*(pow((FTYPE)iter,-1.3));
    }

    //vector of conserved at the previous two iterations
    PLOOP(pliter,pl)  uupp[pl]=uup[pl]; // uupp will have solution for inversion: P(uupp)
    PLOOP(pliter,pl)  uup[pl]=uu[pl]; // uup will not necessarily have P(uup) because uu used Newton step.
    PLOOP(pliter,pl)  uuporig[pl]=uu[pl];

    PLOOP(pliter,pl)  pppp[pl]=ppp[pl]; // pppp will have solution for inversion
    PLOOP(pliter,pl)  ppp[pl]=pp[pl]; // ppp will not necessarily have inversion because pp used Newton step.
    PLOOP(pliter,pl)  ppporig[pl]=pp[pl];

    fracdtuu0pp=fracdtuu0p; // fracdtuu0 used when computing previous f1 and f2's
    fracdtuu0p=fracdtuu0; // fracdtuu0 used when computing previous f1 and f2's

    fracdtGpp=fracdtGp; // fracdtG used when computing previous f1 and f2's
    fracdtGp=fracdtG; // fracdtG used when computing previous f1 and f2's

    errorabsp=errorabslocal;



    // KORALTODO: Once use certain uu0 value and succeed in getting f1, not enough.  But if take a step and *then* good f1, then know that original U was good choice and can restart from that point instead of backing up uu0 again.
    // KORALTODO: Look at whether bounding error by sign as in bisection (only true in 1D!!), and if approaching 0 slowly enough and still hit CASE, then must be real CASE not a stepping issue.
    // KORALTODO: Getting f1 is just about f1(U) as far as radiation is concerned.  So as far as CASE issues, getting f1(U) means we are good with that U and we can certainly stick with the used uu0.

    
    /////////////////
    //
    // get error function (f1) and inversion (uu->pp) using uu
    //
    /////////////////
    f1break=0;

    for(f1iter=0;f1iter<MAXF1TRIES;f1iter++){
      int whichcall=1;
      eomtypelocal=*eomtype; // default for f1 at first.
      failreturn=f_implicit_lab(iter,failreturnallowableuse, whichcall,showmessages, allowlocalfailurefixandnoreport, &eomtypelocal, pp, uu0, uu, fracdtG*realdt, ptrgeom, f1, f1norm); // modifies uu and pp

      // regardless of failure, check if below machine precision for 4-force and abort if so.
      // below ==1.0 assumes powers of 2 for changing those things so definitely comes back to 1.0 exactly.
      int dimtypef=DIMTYPEFCONS;
      if(fracdtuu0==1.0 && fracdtG==1.0 && f_error_check(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,realdt,dimtypef,f1,f1norm,f3report,Uiin, uu0,uu,ptrgeom)){
        if(showmessagesheavy&&debugfail>=2) dualfprintf(fail_file,"Very Initial error near machine precision.\n");
        // this is also required to catch cases when G is below machine error and won't affect result.
        // break, and next coming f_error_check() post-f1iter loop will catch this and break again.
        break;
      }
  
      if(failreturn){

#define BACKUPRELEASEFAIL (100.0*NUMEPSILON)  // (1E-30L)
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
          PLOOP(pliter,pl) ppp[pl]=pp[pl]=pp0[pl];
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
        // KORALNOTE: No need to get pp0.  uu0 just used in error function and in setting final dU from uu and uu0.


        // keep below so can count inversion failures against retry successes in the failure file.
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f1 failed: iter=%d  Backing-up both uu0 and G.: f1iter=%d fracdtuu0=%g fracdtG=%g fracuup=%g\n",iter,f1iter,fracdtuu0,fracdtG,fracuup);
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl]);
      }// end if failed inversion in f_implicit_lab()
      else{
        // then success, so was able to do inversion P(U) with change in radiation on fluid: P(uu[fluid] = uu0[fluid] - (uu[rad]-uu0[rad]))
        // This doesn't necessarily mean could do P(uu0) except for iter=1.
        // Corresponds to success for a certain P(uu0,uu) pair.
        if(showmessagesheavy) PLOOP(pliter,pl) dualfprintf(fail_file,"SUCCESS: pl=%d Ui=%21.15g uu0=%21.15g uu0orig=%21.15g uu=%21.15g uup=%21.15g dUother=%21.15g: pp(pnew)=%g : fracdtuu0=%g fracdtG=%g fracuup=%g\n",pl,Uiin[pl],uu0[pl],uu0orig[pl],uu[pl],uup[pl],dUother[pl],pp[pl],fracdtuu0,fracdtG,fracuup);
        break;
      }
    }// end loop over f1iter
    if(f1iter==MAXF1TRIES){
      if(debugfail>=2) dualfprintf(fail_file,"Reached MAXF1TRIES: nstep=%ld steppart=%d ijk=%d %d %d : iter=%d : fracdtuu0=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k,iter,fracdtuu0);
      failnum++;
      *returntype=RETURNIMPLICITFAIL;
      break;
    }
    else{
      // restore fracuup back to 1 since this is only meant to adjust how much go back to previous uu to be able to get f1 computed.
      // fracuup doesn't stay <1.0 because each attempt to get f1 is independent.
      // KORALNOTE: Perhaps reasonable and safer to keep fracuup as not reverted back to 1, since apparently unable to take full steps.  This effectively damps stepping.
      //      fracuup=1.0; // KORALTODO.
    }

    // diagnose
    numoff1iter += f1iter;
    if(showmessagesheavy) dualfprintf(fail_file,"i=%d f1: %g %g %g %g\n",ptrgeom->i,f1[erefU[0]],f1[erefU[1]],f1[erefU[2]],f1[erefU[3]]);



    // see if pre-convergence (happens if force is small or no force at all.  Can't necessarily continue since Jacobian can require arbitrarily large dU on fluid and fail to invert even if no fluid-radiation interaction!
    //test pre-convergence using initial |dU/U|
    // KORALTODO: This isn't a completely general error check since force might be large for fluid that needs itself to have more accuracy, but if using ~NUMEPSILON, won't resolve 4-force of radiation on fluid to better than that.
    int dimtypef=DIMTYPEFCONS;
    if(f_error_check(showmessages, showmessagesheavy, iter, LOCALPREIMPCONV,realdt,dimtypef,f1,f1norm,f3report,Uiin, uu0,uu,ptrgeom)){
      if(showmessagesheavy&&debugfail>=2) dualfprintf(fail_file,"Initial error near machine precision.\n");
      // this is also required to catch cases when G is below machine error and won't affect result.
      f1break=1;
      *returntype=RETURNIMPLICITTRY;
      break;
    }


    int notfinite;
    if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
      notfinite = !isfinitel(uu[irefU[0]])|| !isfinitel(uu[irefU[1]])|| !isfinitel(uu[irefU[2]])|| !isfinitel(uu[irefU[3]]) || !isfinitel(uup[irefU[0]])|| !isfinitel(uup[irefU[1]])|| !isfinitel(uup[irefU[2]])|| !isfinitel(uup[irefU[3]]);
    }
    else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
      notfinite = !isfinitel(pp[irefU[0]])|| !isfinitel(pp[irefU[1]])|| !isfinitel(pp[irefU[2]])|| !isfinitel(pp[irefU[3]]) || !isfinitel(ppp[irefU[0]])|| !isfinitel(ppp[irefU[1]])|| !isfinitel(ppp[irefU[2]])|| !isfinitel(ppp[irefU[3]]);
    }


    if(!notfinite){
    
      /////////
      //
      // get Jacobian and inversion Jacobian 
      //
      /////////
      //      eomtypelocal=*eomtype; // default .  No, use same as f1 determined.
      int failreturniJ=get_implicit_iJ(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, &eomtypelocal, impepsjac, uu, uup, uu0, pp, ppp, fracdtG, realdt, ptrgeom, f1, f1norm, iJ);
      if(failreturniJ!=0){
        *returntype=RETURNIMPLICITFAIL;
        break;
      }

      if(showmessagesheavy){
        int iii,jjj;
        DLOOP(iii,jjj) dualfprintf(fail_file,"iJ[i %d][e %d]=%g\n",iii,jjj,iJ[irefU[iii]][erefU[jjj]]);
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
      if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
        PLOOP(pliter,pl) uu[pl] = uup[pl];
        DLOOP(ii,jj) uu[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
      }
      else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
        PLOOP(pliter,pl) pp[pl]=ppp[pl];
        DLOOP(ii,jj) pp[irefU[ii]] -= DAMPFACTOR*iJ[irefU[ii]][erefU[jj]]*f1[erefU[jj]];
      }
    
    
      if(showmessagesheavy){
        dualfprintf(fail_file,"POSTDX: uu: %g %g %g %g : uup=%g %g %g %g\n",uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],uup[irefU[0]],uup[irefU[1]],uup[irefU[2]],uup[irefU[3]]);
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
        // KORALNOTE: No need to get pp0.  uu0 just used in error function and in setting final dU from uu and uu0.
        checkconv=0;
      }
      //    if(fracdtG!=1.0){
      if(fabs(fracdtG-1.0)>10.0*NUMEPSILON && changeotherdt){
        // try increasing amount of G applied
        fracdtG*=RADDAMPUNDELTA;
        checkconv=0;
      }      



      /////////
      //
      // test convergence after Newton step
      // KORALTODO: This isn't a completely general error check since force might be large for fluid.  So using (e.g.) 1E-6 might still imply a ~1 or larger error for the fluid.  Only down to ~NUMEPSILON will radiation 4-force be unresolved as fluid source term.
      // NOTE: Have to be careful with decreasing DAMPFACTOR or fracdtuu0 because can become small enough that apparently fake convergence with below condition, so only check for convergence if all DAMPs are 1.0.
      /////////
      // 0 = conserved R^t_\nu type, 1 = primitive R^{ti} type
      if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){ // still considering iterate error
        dimtypef=DIMTYPEFCONS;
        DLOOPA(ii){
          f3[erefU[ii]]=(uu[irefU[ii]]-uup[irefU[ii]]);
          f3norm[erefU[ii]]=fabs(uu[irefU[ii]])+fabs(uup[irefU[ii]]);
        }
      }
      else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){ // still considering iterate error
        dimtypef=DIMTYPEFPRIM;
        DLOOPA(ii){
          f3[erefU[ii]]=(pp[irefU[ii]]-ppp[irefU[ii]]);
          f3norm[erefU[ii]]=fabs(pp[irefU[ii]])+fabs(ppp[irefU[ii]]);
        }
      }
  
      convreturn=f_error_check(showmessages, showmessagesheavy, iter, imptryconv,realdt,dimtypef,f3,f3norm,f3reportuu,Uiin,uup,uu,ptrgeom);
      // store error and solution in case eventually lead to max iterations and actually get worse error
      // errorabslocal=0.0;      DLOOPA(jj) errorabslocal     += fabs(f3reportuu[erefU[jj]]);// use for machine precision test of uu-uup.  FUCK
      errorabslocal=0.0;     DLOOPA(jj) errorabslocal     += fabs(f3report[erefU[jj]]); // try using real error
      //DLOOPA(jj) f3report[erefU[jj]]=f3reportuu[erefU[jj]]; // FUCK
      DLOOPA(jj) f3reportuu[erefU[jj]]=f3report[erefU[jj]]; // FUCK

      if(checkconv){ // only store best if settled into full solution
        // store error and solution in case eventually lead to max iterations and actually get worse error
        errorabsbest=0.0; DLOOPA(jj) errorabsbest += fabs(lowestfreport[erefU[jj]]);
        if(errorabsbest>errorabslocal && isfinitel(errorabslocal)){
          //        dualfprintf(fail_file,"Got best: %ld %d : %d %d %d : jj=%d uu=%g bestuu=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->j,jj,uu[irefU[jj]],bestuu[irefU[jj]]);
          //          dualfprintf(fail_file,"Choosing best: %ld %d : %d %d %d  iter=%d : errl=%g errb=%g : uu=%g %g %g %g : uubest=%g %g %g %g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->j,iter,errorabslocal,errorabsbest,uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],bestuu[irefU[0]],bestuu[irefU[1]],bestuu[irefU[2]],bestuu[irefU[3]]);
          PLOOP(pliter,pl) bestuu[pl]=uu[pl];
          PLOOP(pliter,pl) bestpp[pl]=pp[pl];
          DLOOPA(jj) lowestfreport[erefU[jj]]=f3report[erefU[jj]];
          gotbest=1;
          iterbest=iter;
        }
      }


      // DAMP CONTROL
      if(dampstrategy==2){
        if(iter==1){
          // un-damp first step, but still don't go back to 1.0.
          DAMPFACTOR=0.7;
        }
#define LOWESTDAMP (0.05L)
        // see if need to damp, but don't damp below some point.
        if(errorabslocal>errorabsp && DAMPFACTOR>LOWESTDAMP) DAMPFACTOR*=0.5;
      }
  


      // check convergence
      if(checkconv){
        if(convreturn){
          *returntype=RETURNIMPLICITTRY;
          break;
        }
      } // end if checking convergence

    }// end if finite


    /////////
    // see if took too many Newton steps or not finite results
    /////////
    int itermaxed=iter>impmaxiter;
    if(itermaxed || notfinite ){
      if(gotbest){
        if(showmessages && debugfail>=2) dualfprintf(fail_file,"Using best: %g %g\n",errorabslocal,errorabsbest);
        // get new convreturn
        convreturn=(lowestfreport[erefU[0]]<impallowconv && lowestfreport[erefU[1]]<impallowconv && lowestfreport[erefU[2]]<impallowconv && lowestfreport[erefU[3]]<impallowconv);
      }
      else{
        convreturn=(f3report[erefU[0]]<impallowconv && f3report[erefU[1]]<impallowconv && f3report[erefU[2]]<impallowconv && f3report[erefU[3]]<impallowconv);
      }

      // KORALTODO: If convreturn doesn't work, but still (say) 10% error, might want to hold onto result in case explicit backup fails as well (which is likely), in which case *much* better to use 10% error because otherwise 4-force not accounted for, which can lead to very big changes in fluid behavior due to large flux from previous step.
      // KORALTODO: Or, perhaps should really take it as a failure and use fixups.  Probably should allow for result to be written if error<10%, but only use as super-backup in fixups.  So should set pflag still.

      if(convreturn){
        if(1||showmessages && debugfail>=2) dualfprintf(fail_file,"iter>impmaxiter=%d : iter exceeded in solve_implicit_lab().  But f3 was allowed error. checkconv=%d (if checkconv=0, could be issue!) : %g %g %g %g : %g %g %g %g : errorabslocal=%g : %g %g %g\n",impmaxiter,checkconv,f3report[erefU[0]],f3report[erefU[1]],f3report[erefU[2]],f3report[erefU[3]],lowestfreport[erefU[0]],lowestfreport[erefU[1]],lowestfreport[erefU[2]],lowestfreport[erefU[3]],errorabslocal,fracdtuu0,fracuup,fracdtG);
        // NOTE: If checkconv=0, then wasn't ready to check convergence and smallness of f3 might only mean smallness of fracuup.  So look for "checkconv=0" cases in fail output.
        *returntype=RETURNIMPLICITALLOW;
        break;
      }
      else{
        // KORALTODO: Need backup that won't fail.
        if(debugfail>=2){
          if(itermaxed) dualfprintf(fail_file,"iter>impmaxiter=%d : iter exceeded in solve_implicit_lab(). nstep=%ld steppart=%d ijk=%d %d %d :  Bad error.\n",impmaxiter,nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          if(notfinite) dualfprintf(fail_file,"IMPGOTNAN at iter=%d : in solve_implicit_lab(). ijk=%d %d %d :  Bad error.\n",iter,ptrgeom->i,ptrgeom->j,ptrgeom->k);
          dualfprintf(fail_file,"checkconv=%d failreturnallowable=%d: %g %g %g %g : %g %g %g %g\n",checkconv,failreturnallowable,f3report[erefU[0]],f3report[erefU[1]],f3report[erefU[2]],f3report[erefU[3]],lowestfreport[erefU[0]],lowestfreport[erefU[1]],lowestfreport[erefU[2]],lowestfreport[erefU[3]]);
        }
        *returntype=RETURNIMPLICITFAILALLOW;
        break;
      }
    }//end if maximum iterations
      
  }// end do
  while(1);



  // diagnose
  numofiter+=iter;
  // i=j=k=0 just to show infrequently
  if(debugfail>=2 && (ptrgeom->i==0 && ptrgeom->j==0  && ptrgeom->k==0) || showmessagesheavy) dualfprintf(fail_file,"numimplicits=%lld averagef1iter=%g averageiter=%g\n",numimplicits,(FTYPE)numoff1iter/(FTYPE)numimplicits,(FTYPE)numofiter/(FTYPE)numimplicits);
  


  if(*returntype>RETURNIMPLICITFAIL){

    // final actual errorabs
    errorabslocal=0.0;
    if(f1break==0 && gotbest){
      DLOOPA(jj){
        //        dualfprintf(fail_file,"Got best: %ld %d : %d %d %d : jj=%d uu=%g bestuu=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->j,jj,uu[irefU[jj]],bestuu[irefU[jj]]);
        errorabslocal += fabs(lowestfreport[erefU[jj]]);
      }
      PLOOP(pliter,pl){
        uu[pl]=bestuu[pl];
        pp[pl]=bestpp[pl];
      }
      FTYPE errorabslocaltest=0.0;
      DLOOPA(jj) errorabslocaltest += fabs(f3report[erefU[jj]]);
      //      dualfprintf(fail_file,"Got best final: %ld %d : %d %d %d : error=%g lasterror=%g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->j,errorabslocal,errorabslocaltest);
      //      dualfprintf(fail_file,"Best: %ld %d : %d %d %d  iter=%d : uu=%g %g %g %g : uubest=%g %g %g %g\n",nstep,steppart,ptrgeom->i,ptrgeom->j,ptrgeom->j,iter,uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],bestuu[irefU[0]],bestuu[irefU[1]],bestuu[irefU[2]],bestuu[irefU[3]]);
      //      dualfprintf(fail_file,"iter=%d iterbest=%d\n",iter,iterbest);
      //      if(iter!=iterbest) dualfprintf(fail_file,"HEREIT\n");
    }
    else DLOOPA(jj) errorabslocal += fabs(f3report[erefU[jj]]);



    ///////////////////
    //
    // get source update as "dU" = dU/dt using real dt that used during implicit iterations, and will eventually use to update U in advance.c.
    // apply source update as force
    PLOOP(pliter,pl) radsource[pl] = +(uu[pl]-uu0[pl])/realdt;

    // save better guess for later inversion from this inversion
    // pp was modified by f_implicit_lab(f1,pp) with output from inversion returned through pp0
    // only use pp if successful with implicit method, since if not successful can be various bad reasons with no good pb
    PLOOP(pliter,pl) pb[pl]=pp[pl];

    ///////////
    //
    // choose new eomtype for external inverison or other checks
    //
    // This ensures whatever implicit solver settled on using, the external inversion doesn't switch.
    //
    //////////
    *eomtype=eomtypelocal;



    // for returning of error
    *errorabs=errorabslocal;
    *iterreturn=iter;
  }


  if(*returntype==RETURNIMPLICITFAIL || *returntype==RETURNIMPLICITFAILALLOW){
    mathematica_report_check(0, failnum, gotfirstnofail, *errorabs, iter, realdt, ptrgeom, ppfirst,pp,pb,piin,uu0,uu,Uiin,Ufin, CUf, q, dUother);
    //      if(doingit==1) myexit(10000000); // DEBUG
    // Note that if inversion reduces to entropy or cold, don't fail, so passes until reached this point.  But convergence can be hard if flipping around which EOMs for the inversion are actually used.
  }

  if(*returntype==RETURNIMPLICITFAIL){
    return(1);
  }
  

  return(0);
  
}

int mathematica_report_check(int failtype, long long int failnum, int gotfirstnofail, FTYPE errorabs, int iters, FTYPE realdt,struct of_geom *ptrgeom, FTYPE *ppfirst, FTYPE *pp, FTYPE *pb, FTYPE *piin, FTYPE *uu0, FTYPE *uu, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_state *q, FTYPE *dUother)
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
    // 181 things
    dualfprintf(fail_file,"\nFAILINFO: %d %d %lld %d %21.15g %d %21.15g %lld %d %21.15g ",failtype,myid,failnum,gotfirstnofail,errorabs,iters,realdt,nstep,steppart,gam); // 8
    DLOOP(jj,kk) dualfprintf(fail_file,"%21.15g ",ptrgeom->gcon[GIND(jj,kk)]); // 16
    DLOOP(jj,kk) dualfprintf(fail_file,"%21.15g ",ptrgeom->gcov[GIND(jj,kk)]); // 16
    PLOOP(pliter,pl) dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ",pp[pl],ppfirst[pl],pb[pl],piin[pl],uu0[pl],uu[pl],Uiin[pl]);  // 7*13
    struct of_state qpp;
    get_state(pp,ptrgeom,&qpp);
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpp.uradcon[jj],qpp.uradcov[jj]); // 4*2=8
    else DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpp.ucon[jj],qpp.ucov[jj]); // 4*2=8
    struct of_state qpb;
    get_state(pb,ptrgeom,&qpb);
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpb.uradcon[jj],qpb.uradcov[jj]); // 4*2=8
    else dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpb.ucon[jj],qpb.ucov[jj]); // 4*2=8
    struct of_state qpiin;
    get_state(piin,ptrgeom,&qpiin);
    if(EOMRADTYPE!=EOMRADNONE) DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpiin.uradcon[jj],qpiin.uradcov[jj]); // 4*2=8
    else dualfprintf(fail_file,"%21.15g %21.15g ",0.0,0.0);
    DLOOPA(jj) dualfprintf(fail_file,"%21.15g %21.15g ",qpiin.ucon[jj],qpiin.ucov[jj]); // 4*2=8
    dualfprintf(fail_file,"\n");

    /////////////
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
    // 1) grep -h --text FAILINFO 0_fail.out.grmhd* | sed 's/FAILINFO: //g' > fails.txt
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
static int f_error_check(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, int dimtypef, FTYPE *fin, FTYPE *finnorm, FTYPE *finreport, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom)
{
  int ii,jj;

  if(DOOLD){
    FTYPE fin4[NDIM],finnorm4[NDIM],finreport4[NDIM];
    DLOOPA(ii){
      fin4[ii]=fin[erefU[ii]];
      finnorm4[ii]=finnorm[erefU[ii]];
    }
    int convreturn=f_error_check_old(showmessages, showmessagesheavy, iter, conv,realdt,fin4,finnorm4,finreport4,Uiin,uu0,uu,ptrgeom);
    DLOOPA(ii) finreport[erefU[ii]]=finreport4[ii];
    return(convreturn);
  }



  FTYPE fina[NPR],finb[NPR],finc[NPR],find[NPR];
  int passedconv;

  // default
  passedconv=0;

  // get error
  // NOTE: use of gcov[ii,ii] so comparable dimensionally to fin[erefU[ii]] and finnorm[erefU[ii]] that are like R^t_\nu and so need sqrt(gcon[nu,nu]) multiplied on them.  This ensures error is non-dimensional (or, really only ^t dimensional)
  FTYPE dimfactferr[NPR],dimfactU[NPR];
  DLOOPA(ii) dimfactU[erefU[ii]]=sqrt(fabs(ptrgeom->gcon[GIND(ii,ii)]));
  if(dimtypef==DIMTYPEFCONS) DLOOPA(ii) dimfactferr[erefU[ii]]=dimfactU[erefU[ii]]; // assume fin and finnorm are conservative or source (R^t_\nu form)
  else DLOOPA(ii) dimfactferr[erefU[ii]]=1.0/dimfactU[erefU[ii]]; // assume fin and finnorm are primitve (u,vel^i form)

  // replace finnorm -> finnormnew that's already non-dimensionalized
  FTYPE finnormnew[NPR];
  ii=TT; finnormnew[erefU[ii]] = fabs(finnorm[erefU[ii]]*dimfactferr[erefU[ii]]);
  // get spatial contributions as total term so not dominated by small errors in some dimensions
  FTYPE fnormspace=0.0; SLOOPA(ii) fnormspace += fabs(finnorm[erefU[ii]]*dimfactferr[erefU[ii]]);
  SLOOPA(ii) finnormnew[erefU[ii]] = fnormspace;

  // get non-dimensionalized fin
  FTYPE finnew[NPR];
  DLOOPA(ii) finnew[erefU[ii]]=fin[erefU[ii]]*dimfactferr[erefU[ii]];

  // get relative errors
  DLOOPA(ii){
    fina[erefU[ii]]=fabs(finnew[erefU[ii]]/(IMPMINABSERROR+fabs(uu0[erefU[0]]*dimfactU[erefU[0]])));
    finb[erefU[ii]]=fabs(finnew[erefU[ii]]/(IMPMINABSERROR+MAX(fabs(uu0[iotherU[0]]*dimfactU[erefU[0]]),fabs(uu0[erefU[0]]*dimfactU[erefU[0]]))));
    finc[erefU[ii]]=fabs(finnew[erefU[ii]]/(IMPMINABSERROR+MAX(finnormnew[erefU[ii]],fabs(uu0[erefU[0]]*dimfactU[erefU[0]]))));
    // SMALL was (uu0[iotherU[0]]-Uiin[iotherU[0]])
    find[erefU[ii]]=fabs(finnew[erefU[ii]]/(IMPMINABSERROR+MAX(fabs((SMALL)*dimfactU[erefU[0]]),MAX(finnormnew[erefU[ii]],fabs(uu0[erefU[0]]*dimfactU[erefU[0]])))));
  }

  // evaluate error
  if(IMPLICITERRORNORM==1){
    if(fina[erefU[0]]<conv && fina[erefU[1]]<conv && fina[erefU[2]]<conv && fina[erefU[3]]<conv) passedconv=1;
    DLOOPA(ii) finreport[erefU[ii]]=fina[erefU[ii]];
  }
  else if(IMPLICITERRORNORM==2){
    if(finb[erefU[0]]<conv && finb[erefU[1]]<conv && finb[erefU[2]]<conv && finb[erefU[3]]<conv) passedconv=1;
    DLOOPA(ii) finreport[erefU[ii]]=finb[erefU[ii]];
  }
  else if(IMPLICITERRORNORM==3){
    if(finc[erefU[0]]<conv && finc[erefU[1]]<conv && finc[erefU[2]]<conv && finc[erefU[3]]<conv) passedconv=1;
    DLOOPA(ii) finreport[erefU[ii]]=finc[erefU[ii]];
  }
  else if(IMPLICITERRORNORM==4){
    if(find[erefU[0]]<conv && find[erefU[1]]<conv && find[erefU[2]]<conv && find[erefU[3]]<conv) passedconv=1;
    DLOOPA(ii) finreport[erefU[ii]]=find[erefU[ii]];
  }
  

  // see if passed convergence test
  if(passedconv){
    if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g realdt=%g i=%d iter=%d DONE1 for conv=%g : finreport=%g %g %g %g\n",nstep,steppart,dt,realdt,ptrgeom->i,iter,conv,finreport[erefU[0]],finreport[erefU[1]],finreport[erefU[2]],finreport[erefU[3]]);
    return(1);
  }
  else{
    // report if didn't pass
    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTFIN (conv=%21.15g): uu: %21.15g %21.15g %21.15g %21.15g : uu0=%21.15g %21.15g %21.15g %21.15g\n",conv,uu[irefU[0]],uu[irefU[1]],uu[irefU[2]],uu[irefU[3]],uu0[irefU[0]],uu0[irefU[1]],uu0[irefU[2]],uu0[irefU[3]]);
      DLOOPA(ii) dualfprintf(fail_file,"iii=%d fin=%21.15g finnorm=%21.15g\n",ii,fin[erefU[ii]],finnorm[erefU[ii]]);
      dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,finreport[erefU[0]],finreport[erefU[1]],finreport[erefU[2]],finreport[erefU[3]]);
    }
    return(0);
  }

  return(0);
}


#define JDIFFONESIDED 0
#define JDIFFCENTERED 1

// choose:
#if(IMPLICITITER==QTYPMHD)
// with IMPLICITITER==QTYPMHD, no longer expensive so can do JDIFFCENTERED
#define JDIFFTYPE JDIFFCENTERED
//#define JDIFFTYPE JDIFFONESIDED
#else
#define JDIFFTYPE JDIFFONESIDED
#endif

// calculating approximate Jacobian: dUresid(dUrad,G(Urad))/dUrad = dy(x)/dx
// then compute inverse Jacobian
static int get_implicit_iJ(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, int *eomtypelocal, FTYPE impepsjac, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pp, FTYPE *ppp, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NPR])
{
  int ii,jj;


  if(DOOLD){
    FTYPE iJ44[NDIM][NDIM],f14[NDIM],f1norm4[NDIM];
    DLOOPA(ii){
      f14[ii]=f1[erefU[ii]];
      f1norm4[ii]=f1norm[erefU[ii]];
    }
    int failreturniJ=get_implicit_iJ_old(failreturnallowableuse, showmessages, showmessagesheavy, allowlocalfailurefixandnoreport, impepsjac, uu, uup, uu0, pp, fracdtG, realdt, ptrgeom, f14, f1norm4, iJ44);
    DLOOP(ii,jj) iJ[irefU[ii]][erefU[jj]]=iJ44[ii][jj];
    return(failreturniJ);
  }

  int eomtypelocallocal=*eomtypelocal; // default



  // ensure uu and pp don't get modified by del-shifts to get Jacobian, which can change primitives to order unity at high radiation gamma
  FTYPE uujac[NPR],ppjac[NPR];
  int pliter,pl;

  // for scaling del's norm.  Applies to time component to make as if like space component.
  FTYPE upitoup0[NPR];
  FTYPE x[NPR],xp[NPR],xjac[2][NPR];
  FTYPE velmomscale;
  if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
    PLOOP(pliter,pl) x[pl]=uu[pl];
    PLOOP(pliter,pl) xp[pl]=uup[pl];
    PLOOP(pliter,pl) xjac[0][pl]=xjac[1][pl]=uu[pl];
    DLOOPA(jj) upitoup0[jj] = sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)]));  // R^t_nu * sqrt(g^{ii}) = R^t_orthonu
    // velmomscale is set such that when (e.g.) T^t_i is near zero, we use T^t_t as reference since we consider T^t_i/T^t_t to be velocity scale that is up to order unity.
    velmomscale=fabs(x[irefU[0]]*upitoup0[0]);
  }
  else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
    PLOOP(pliter,pl) x[pl]=pp[pl];
    PLOOP(pliter,pl) xp[pl]=ppp[pl];
    PLOOP(pliter,pl) xjac[0][pl]=xjac[1][pl]=pp[pl];
    jj=TT; upitoup0[jj] = 1.0; // comoving quantity, not v^t
    //DLOOPA(jj) upitoup0[jj] = 1.0/sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])); // v^i / sqrt(g^{ii}) = vortho^i
    SLOOPA(jj) upitoup0[jj] = 1.0/sqrt(fabs(ptrgeom->gcon[GIND(jj,jj)])); // v^i / sqrt(g^{ii}) = vortho^i
    // for velocity, assume ortho-scale is order unity (i.e. v=1.0*c)
    velmomscale=1.0; // reference scale considered to be order unity
    //velmomscale=fabs(x[irefU[0]]*upitoup0[0]);
  }

  // form pre-del
  FTYPE delspace,predel[NDIM];
  jj=TT; predel[jj] = fabs(x[irefU[jj]]);
  // form delspace that absorbs all spatial values into a single dimensionless scale to avoid one dimensions smallness causing issues.
  // KORALTODO: Maybe causes problems if Jacobian becomes singular because no change in small velocity-momentum values when should be change.  But can't resolve it really, so should be ok.
  delspace=0.0; SLOOPA(jj) delspace = MAX(delspace,MAX(fabs(x[irefU[jj]]*upitoup0[jj]) , velmomscale )); // dimensionless-ortho
  SLOOPA(jj) predel[jj] = delspace/upitoup0[jj]; // back to actual spatial dimension-space scale for application to x and xjac


  FTYPE J[NPR][NPR];
  FTYPE f2[2][NPR],f2norm[2][NPR];
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

    DLOOPA(jj){

      int sided,signside;
      int numsides=2;
      for(sided=0;sided<numsides;sided++){
        if(JDIFFTYPE==JDIFFONESIDED && sided==0) continue;
        if(sided==1) signside=-1.0;
        else signside=+1.0;

        while(1){

          // when |irefU[0]|>>|irefU[1]|, then can't get better than machine error on irefU[0], not irefU[1], so using small del just for irefU[1] makes no sense, so avoid above
          del = localIMPEPS*predel[jj];
         
          // origin point
          PLOOP(pliter,pl) xjac[sided][pl]=x[pl];
          // offset xjac (KORALTODO: How to ensure this doesn't have machine precision problems or is good enough difference?)
          xjac[sided][irefU[jj]]=x[irefU[jj]] + signside*del; // KORALNOTE: Not sure why koral was using uup or xp here.  Should use x (or uu or pp) because as updated from f_implicit_lab() and uup or ppp for xp hasn't been set yet, so not consistent with desired jacobian or ferr for Newton step.

          //          dualfprintf(fail_file,"NEW: jj=%d del=%g xjac=%g x=%g\n",jj,del,xjac[sided][irefU[jj]],x[irefU[jj]]);

          // set uujac and ppjac using xjac
          if(IMPLICITITER==QTYUMHD || IMPLICITITER==QTYURAD || IMPLICITITER==QTYENTROPYUMHD){
            PLOOP(pliter,pl){
              uujac[pl]=xjac[sided][pl];
              ppjac[pl]=pp[pl];
            }
          }
          else if(IMPLICITITER==QTYPMHD || IMPLICITITER==QTYPRAD || IMPLICITITER==QTYENTROPYPMHD){
            PLOOP(pliter,pl){
              uujac[pl]=uu[pl];
              ppjac[pl]=xjac[sided][pl];
            }
          }
 
          // get dUresid for this offset xjac
          int whichcall=2;
          eomtypelocallocal=*eomtypelocal; // re-default
          int fakeiter=1; // ok to let Trad and Tgas switch a bit due to small delta to get Jacobian.
          failreturn=f_implicit_lab(fakeiter,failreturnallowableuse, whichcall,showmessages,allowlocalfailurefixandnoreport, &eomtypelocallocal, ppjac,uu0,uujac,fracdtG*realdt,ptrgeom,f2[sided],f2norm[sided]);
          if(failreturn){
            if(showmessages&& debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed: jj=%d.  Trying smaller localIMPEPS=%g (giving del=%g) to %g\n",jj,localIMPEPS,del,localIMPEPS*FRACIMPEPSCHANGE);
            localIMPEPS*=FRACIMPEPSCHANGE;
            // try making smaller until no error, unless doesn't work out
            // see if will be able to resolve differences
            int baddiff = fabs(xjac[sided][irefU[jj]]-x[irefU[jj]])/(fabs(xjac[sided][irefU[jj]])+fabs(x[irefU[jj]])) < 10.0*NUMEPSILON;
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
        }// end while(1)
      }// end sided=0,1
      

      // get Jacobian
      DLOOPA(ii) J[erefU[ii]][irefU[jj]]=(f2[1][erefU[ii]] - f2[0][erefU[ii]])/(xjac[1][irefU[jj]]-xjac[0][irefU[jj]]);

      //      DLOOPA(ii) dualfprintf(fail_file,"NEW: ii=%d jj=%d J=%g : %g %g : %g %g\n",ii,jj,J[erefU[ii]][irefU[jj]], f2[1][erefU[ii]],f2[0][erefU[ii]],xjac[1][irefU[jj]],xjac[0][irefU[jj]]);


      // debug info
      if(debugfail>=2){
        DLOOPA(ii) if(showmessagesheavy || !isfinitel(J[erefU[ii]][irefU[jj]])){
          dualfprintf(fail_file,"JAC: xjac[0]: %21.15g %21.15g %21.15g %21.15g :  xjac[1]: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g (del=%21.15g localIMPEPS=%21.15g)\n",
                      xjac[0][irefU[0]],xjac[0][irefU[1]],xjac[0][irefU[2]],xjac[0][irefU[3]],
                      xjac[1][irefU[0]],xjac[1][irefU[1]],xjac[1][irefU[2]],xjac[1][irefU[3]],
                      x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],del,localIMPEPS);
          dualfprintf(fail_file,"i=%d jj=%d f2[0]: %21.15g %21.15g %21.15g %21.15g\n",ptrgeom->i,jj,f2[0][erefU[0]],f2[0][erefU[1]],f2[0][erefU[2]],f2[0][erefU[3]]);
          dualfprintf(fail_file,"i=%d jj=%d f2[1]: %21.15g %21.15g %21.15g %21.15g\n",ptrgeom->i,jj,f2[1][erefU[0]],f2[1][erefU[1]],f2[1][erefU[2]],f2[1][erefU[3]]);
          dualfprintf(fail_file,"JISNAN: %d %d : J=%21.15g : f2[1]=%21.15g f2[0]=%21.15g : xjac[1]=%21.15g xjac[0]=%21.15g dx=%g\n",ii,jj,J[erefU[ii]][irefU[jj]], f2[1][erefU[ii]],f2[0][erefU[ii]],xjac[1][irefU[jj]],xjac[0][irefU[jj]],xjac[1][irefU[jj]]-xjac[0][irefU[jj]]);
        }
      }

    }
    


    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTJAC: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
      int iii,jjj;
      DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%21.15g\n",iii,jjj,J[erefU[iii]][irefU[jjj]]);
    }


    

    /////////////////////
    //
    //invert Jacobian
    //
    /////////////////////
    // copy over matrix to pure 4x4 version
    FTYPE J44[NDIM][NDIM];
    FTYPE iJ44[NDIM][NDIM];
    DLOOP(ii,jj) J44[ii][jj]=J[erefU[ii]][irefU[jj]];

    failreturn=inverse_44matrix(J44,iJ44); // inverse *and* transpose index order, so J[irefU][erefU] ->  iJ[erefU][irefU]

    // copy back inverse matrix from 4x4 version
    DLOOP(ii,jj) iJ[irefU[ii]][erefU[jj]]=iJ44[ii][jj];


    //    DLOOP(ii,jj) dualfprintf(fail_file,"NEW: ii=%d jj=%d iJ=%g\n",ii,jj,iJ[irefU[ii]][erefU[jj]]);
    
      


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
        if(debugfail>=2) dualfprintf(fail_file,"f_implicit_lab for f2 failed to be different enough from f1 and gave singular Jacobian: IMPEPSSTART=%g (giving del=%g)\n",IMPEPSSTART,del);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC1: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J44[%d][%d]=%21.15g\n",iii,jjj,J44[iii][jjj]);
        }
        return(1); // can't expect good derivative above ~0.3, so just return as failure of implicit method.
      }
      else{
        if(debugfail>=2) dualfprintf(fail_file,"inverse_44matrix(J,iJ) failed, trying IMPEPSSTART=%g :: ijk=%d %d %d\n",IMPEPSSTART,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC2: x: %21.15g %21.15g %21.15g %21.15g : x=%21.15g %21.15g %21.15g %21.15g\n",x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]],x[irefU[0]],x[irefU[1]],x[irefU[2]],x[irefU[3]]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J44[%d][%d]=%21.15g\n",iii,jjj,J44[iii][jjj]);
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
        DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%g\n",iii,jjj,J44[iii][jjj]);
      }
      return(1);
    }
  }// end over ensuring Jacobian is non-singular for the given del


  return(0);

}













//////////////////////////////////////////////////////////////////





static int f_implicit_lab_old(int failreturnallowable, int whichcall, int showmessages, int allowlocalfailurefixandnoreport, FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE localdt, struct of_geom *ptrgeom,  FTYPE *f, FTYPE *fnorm)
{
  FTYPE pp[NPR];
  int pliter, pl;
  int iv;
  struct of_newtonstats newtonstats;
  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.
  FTYPE Gdpl[NPR],Gdabspl[NPR], Tgas;


  // get primitive (don't change uu0).  This pp is used for inversion guess and to hold final inversion answer.
  PLOOP(pliter,pl) pp[pl] = pp0[pl];
  // initialize Gdpl
  PLOOP(pliter,pl) Gdpl[pl] = 0.0;

#define SIGNGD2OLD (1.0) // sign that goes into implicit differencer that's consistent with sign for SIGNGD of -1 when using the radiative uu to measure f.


  // get change in conserved quantity between fluid and radiation (equal and opposite 4-force)
  // required for inversion to get P(U) for MHD and RAD variables
  // this preserves machine accurate conservation instead of applying 4-force on each fluid and radiation separately that can accumulate errors
  DLOOPA(iv) uu[UU+iv] = uu0[UU+iv] - (uu[URAD0+iv]-uu0[URAD0+iv]);


#if(0&&DOENTROPY!=DONOENTROPY && ENTROPY!=-100)
  // always get uu[ENTROPY] because even MHD inversion might reduce to use entropy.
  // KORALTODO: Problem is updated U[MHD] from U[RAD] but don't know corresponding primitives yet, so below uu[ENTROPY] will be off.
  koral_source_rad_calc(pp, ptrgeom, Gdpl, NULL, NULL, &Tgas, NULL);
  pl=ENTROPY;
  // sign below consistent with mathematica and definition of Gdpl[ENTROPY]=-u.G/T
  uu[pl] = uu0[pl]; // + (SIGNGD2OLD*localdt*Gdpl[pl]); // I think this is + based upon  mathematica.

  //  dualfprintf(fail_file,"ENTROPY: %g %g %g :  %g %g\n",uu[pl],uu0[pl],(SIGNGD2OLD*localdt*Gdpl[pl]),Gdpl[pl],localdt);
#endif



  //  PLOOP(pliter,pl) dualfprintf(fail_file,"f_implicit_lab: wc=%d i=%d j=%d pl=%d uu=%g\n",whichcall, ptrgeom->i,ptrgeom->j,pl,uu[pl]);
  
  // Get P(U)
  int failreturn;
  failreturn=Utoprimgen_failwrapper_old(showmessages,allowlocalfailurefixandnoreport, finalstep, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats);


  // get 4-force for all pl due to radiation as due to pp[uu]
  koral_source_rad_calc(pp, ptrgeom, Gdpl, Gdabspl, NULL, &Tgas, NULL);


  // compute difference vector between original and new 4-force's effect on conserved radiative quantities
  // NR1992 Eq. 16.6.16: y_{n+1} = y_n + h f(y_{n+1}) , so error function is f = (y_{n+1} - y_n) - h f(y_{n+1})
  // i.e. f->0 as change in conserved quantity approaches the updated value of 4-force
  DLOOPA(iv){
    pl=URAD0+iv;
    f[iv] = (uu[pl] - uu0[pl]) + (SIGNGD2OLD * localdt * Gdpl[pl]);

    // get error normalization that involves actual things being differenced
    // KORALNOTE: Use Gdabspl to ensure absolute value over more terms in source in case coincidental cancellation without physical significance.
    fnorm[iv] = 0.5*(fabs(uu[pl]) + fabs(uu0[pl]) + fabs(SIGNGD2OLD * localdt * Gdabspl[pl]));
  }

#if(0)
  // no, still using radiative quantities.  Only need to set uu[ENTROPY] and get MHD inversion as above.
  if(doingentropy==1){
    // replace energy equation with dS*T equation
    iv=0; pl=ENTROPY;
    // error function is T*dS so no actual division by T.  Found in mathematica that this works best in difficult precision cases.
    f[iv] = ((uu[pl] - uu0[pl]) + (SIGNGD2OLD * localdt * Gdpl[pl]))*Tgas;
    fnorm[iv] = (0.5*(fabs(uu[pl]) + fabs(uu0[pl]) + fabs(SIGNGD2OLD * localdt * Gdabspl[pl])))*Tgas;
  }
#endif



  //  dualfprintf(fail_file,"i=%d fnorm=%g %g %g %g : Gd=%g %g %g %g : uuG=%g %g %g %g\n",ptrgeom->i,fnorm[0],fnorm[1],fnorm[2],fnorm[3],Gdpl[URAD0+0],Gdpl[URAD0+1],Gdpl[URAD0+2],Gdpl[URAD0+3],SIGNGD2OLD * localdt * Gdpl[URAD0],SIGNGD2OLD * localdt * Gdpl[URAD1],SIGNGD2OLD * localdt * Gdpl[URAD2],SIGNGD2OLD * localdt * Gdpl[URAD3]);


  if(failreturn && failreturn>failreturnallowable){
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Utoprimgen_wrapper() failed, must return out of f_implicit_lab(): %d vs. %d\n",failreturn,failreturnallowable);
    return(failreturn);
  }
  else{
    // save better guess for later inversion (including this inversion above) from this inversion
    PLOOP(pliter,pl) pp0[pl]=pp[pl];
  }

  return 0;
} 




// calculating approximate Jacobian: dUresid(dUrad,G(Urad))/dUrad = dy(x)/dx
// then compute inverse Jacobian
static int get_implicit_iJ_old(int failreturnallowableuse, int showmessages, int showmessagesheavy, int allowlocalfailurefixandnoreport, FTYPE impepsjac, FTYPE *uu, FTYPE *uup, FTYPE *uu0, FTYPE *pb, FTYPE fracdtG, FTYPE realdt, struct of_geom *ptrgeom, FTYPE *f1, FTYPE *f1norm, FTYPE (*iJ)[NDIM])
{
  int ii,jj;
  FTYPE J[NDIM][NDIM],f2[NDIM],f2norm[NDIM];

  // for scaling del's norm
  FTYPE sqrtgcov[NDIM];
  sqrtgcov[0]=1.0;
  SLOOPA(jj) sqrtgcov[jj] = sqrt(fabs(ptrgeom->gcov[GIND(jj,jj)]));


  // ensure pb doesn't get modified by del-shifts to get Jacobian, which can change primitives to order unity at high radiation gamma
  FTYPE pbjac[NPR];
  int pliter,pl;
  PLOOP(pliter,pl) pbjac[pl]=pb[pl];


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
 
        //        dualfprintf(fail_file,"OLD: jj=%d del=%g xjac=%g x=%g\n",jj,del,uu[jj+URAD0],uup[jj+URAD0]);
        
        // get dUresid for this offset uu
        int whichcall=2;
        if(0){
          failreturn=f_implicit_lab_old(failreturnallowableuse, whichcall,showmessages,allowlocalfailurefixandnoreport, pbjac,uu0,uu,fracdtG*realdt,ptrgeom,f2,f2norm); // old calling old is ok and necessary
        }
        else{
          FTYPE f2pl[NPR],f2normpl[NPR];
          int eomtype=EOMDEFAULT;
          int fakeiter=1;
          failreturn=f_implicit_lab(fakeiter,failreturnallowableuse, whichcall,showmessages,allowlocalfailurefixandnoreport, &eomtype, pbjac,uu0,uu,fracdtG*realdt,ptrgeom,f2pl,f2normpl);
          int kk;
          DLOOPA(kk){
            f2[kk]=f2pl[erefU[kk]];
            f2norm[kk]=f2normpl[erefU[kk]];
          }
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

      //      DLOOPA(ii) dualfprintf(fail_file,"OLD: ii=%d jj=%d J=%g : %g %g : %g %g\n",ii,jj,J[ii][jj],f2[ii],f1[ii],uu[jj+URAD0],uup[jj+URAD0]);

      if(debugfail>=2){
        DLOOPA(ii) if(showmessagesheavy || !isfinitel(J[ii][jj])){
          dualfprintf(fail_file,"JAC: uu: %21.15g %21.15g %21.15g %21.15g : uup=%21.15g %21.15g %21.15g %21.15g (del=%21.15g localIMPEPS=%21.15g)\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3],del,localIMPEPS);
          dualfprintf(fail_file,"i=%d jj=%d f2: %21.15g %21.15g %21.15g %21.15g\n",ptrgeom->i,jj,f2[0],f2[1],f2[2],f2[3]);
          dualfprintf(fail_file,"JISNAN: %d %d : %21.15g : %21.15g %21.15g : %21.15g %21.15g\n",ii,jj,J[ii][jj],f2[ii],f1[ii],uu[jj+URAD0],uup[jj+URAD0]);
        }
      }

      // restore uu after getting changed by f_implicit_lab(f2)
      uu[jj+URAD0]=uup[jj+URAD0];
    }
    


    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTJAC: uu: %21.15g %21.15g %21.15g %21.15g : uup=%21.15g %21.15g %21.15g %21.15g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
      int iii,jjj;
      DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%21.15g\n",iii,jjj,J[iii][jjj]);
    }

    

    /////////////////////
    //
    //invert Jacobian
    //
    /////////////////////
    failreturn=inverse_44matrix(J,iJ);

    
    //    DLOOP(ii,jj) dualfprintf(fail_file,"OLD: ii=%d jj=%d iJ=%g\n",ii,jj,iJ[ii][jj]);


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
          dualfprintf(fail_file,"POSTJAC1: uu: %21.15g %21.15g %21.15g %21.15g : uup=%21.15g %21.15g %21.15g %21.15g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%21.15g\n",iii,jjj,J[iii][jjj]);
        }
        return(1); // can't expect good derivative above ~0.3, so just return as failure of implicit method.
      }
      else{
        if(debugfail>=2) dualfprintf(fail_file,"inverse_44matrix(J,iJ) failed, trying IMPEPSSTART=%g :: ijk=%d %d %d\n",IMPEPSSTART,ptrgeom->i,ptrgeom->j,ptrgeom->k);
        if(debugfail>=2 || showmessagesheavy){
          dualfprintf(fail_file,"POSTJAC2: uu: %21.15g %21.15g %21.15g %21.15g : uup=%21.15g %21.15g %21.15g %21.15g\n",uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uup[URAD0],uup[URAD1],uup[URAD2],uup[URAD3]);
          int iii,jjj;
          DLOOP(iii,jjj) dualfprintf(fail_file,"J[%d][%d]=%21.15g\n",iii,jjj,J[iii][jjj]);
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





// use f and check the error
static int f_error_check_old(int showmessages, int showmessagesheavy, int iter, FTYPE conv, FTYPE realdt, FTYPE *f1, FTYPE *f1norm, FTYPE *f3report, FTYPE *Uiin, FTYPE *uu0, FTYPE *uu, struct of_geom *ptrgeom)
{
  int ii,jj;
  FTYPE f3[NDIM];
  FTYPE f3a[NDIM],f3b[NDIM],f3c[NDIM],f3d[NDIM];
  int passedconv;

  // default
  passedconv=0;

  // get error
  // NOTE: use of gcov[ii,ii] so comparable dimensionally to f1[ii] and f1norm[ii] that are like R^t_\nu and so need sqrt(gcon[nu,nu]) multiplied on them.  This ensures error is non-dimensional (or, really only ^t dimensional)
  FTYPE dimfact[NDIM];
  DLOOPA(ii){
    dimfact[ii]=sqrt(fabs(ptrgeom->gcon[GIND(ii,ii)]));
  }    
  DLOOPA(ii){
    f3a[ii]=fabs(f1[ii]*dimfact[ii]/(IMPMINABSERROR+fabs(uu0[URAD0]*dimfact[0])));
    f3b[ii]=fabs(f1[ii]*dimfact[ii]/(IMPMINABSERROR+MAX(fabs(uu0[UU]*dimfact[0]),fabs(uu0[URAD0]*dimfact[0]))));
    f3c[ii]=fabs(f1[ii]*dimfact[ii]/(IMPMINABSERROR+MAX(fabs(f1norm[ii]*dimfact[ii]),fabs(uu0[URAD0]*dimfact[0]))));
    //uu0[UU]-Uiin[UU]
    f3d[ii]=fabs(f1[ii]*dimfact[ii]/(IMPMINABSERROR+MAX(fabs((SMALL)*dimfact[0]),MAX(fabs(f1norm[ii]*dimfact[ii]),fabs(uu0[URAD0]*dimfact[0])))));
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
  

  // see if passed convergence test
  if(passedconv){
    if(showmessagesheavy) dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g realdt=%g i=%d iter=%d DONE1 for conv=%g : f3report=%g %g %g %g\n",nstep,steppart,dt,realdt,ptrgeom->i,iter,conv,f3report[0],f3report[1],f3report[2],f3report[3]);
    return(1);
  }
  else{
    // report if didn't pass
    if(showmessagesheavy){
      dualfprintf(fail_file,"POSTF1 (conv=%21.15g): uu: %21.15g %21.15g %21.15g %21.15g : uu0=%21.15g %21.15g %21.15g %21.15g\n",conv,uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3],uu0[URAD0],uu0[URAD1],uu0[URAD2],uu0[URAD3]);
      int iii;
      DLOOPA(iii) dualfprintf(fail_file,"iii=%d f1=%21.15g f1norm=%21.15g\n",iii,f1[iii],f1norm[iii]);
      dualfprintf(fail_file,"nstep=%ld steppart=%d dt=%g i=%d iter=%d : %g %g %g %g\n",nstep,steppart,dt,ptrgeom->i,iter,f3report[0],f3report[1],f3report[2],f3report[3]);
    }
    return(0);
  }

  return(0);
}



static int Utoprimgen_failwrapper_old(int showmessages, int allowlocalfailurefixandnoreport, int finalstep, int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats)
{
  int failreturn;

  // defaults
  failreturn=0;

  // KORALTODO: 
  // flag needs to be reset to preexistingfail(gas/rad) is not a failure.  Only use preexisting catches in utoprimgen.c if done with 4-force and report error in pflag and eventually go to the final advance.c inversion.
  GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
  GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
  


  //calculating primitives  
  // OPTMARK: Should optimize this to  not try to get down to machine precision
  int eomtype=EOMDEFAULT;
  MYFUN(Utoprimgen(showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtype, evolvetype, inputtype, U, ptrgeom, pr, newtonstats),"phys.tools.rad.c:Utoprimgen_failwrapper()", "Utoprimgen", 1);

  // check how inversion did.  If didn't succeed, then check if soft failure and pass.  Else if hard failure have to return didn't work.
  int lpflag,lpflagrad;
  lpflag=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL);
  lpflagrad=GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL);
  if(IFUTOPRIMFAILSOFT(lpflag)){
    // assume soft failure ok, but reset
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMFAIL)=UTOPRIMNOFAIL;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got soft MHD failure inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  else if(IFUTOPRIMRADFAIL(lpflagrad)){
    // can reduce Newton step if getting failure.
    // reset pflag for radiation to no failure, but treat here locally
    GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADNOFAIL;
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got some radiation inversion failure during Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILRAD;
  }
  else if( IFUTOPRIMFAIL(lpflag) || IFUTOPRIMRADFAIL(lpflagrad) ){
    // these need to get fixed-up, but can't, so return failure
    if(showmessages && debugfail>=2) dualfprintf(fail_file,"Got hard failure of inversion (MHD part only considered as hard) in f_implicit_lab(): ijk=%d %d %d : %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,lpflag,lpflagrad);
    failreturn=UTOPRIMGENWRAPPERRETURNFAILMHD;
  }
  else{
    // no failure
    // dualfprintf(fail_file,"No failure in Utoprimgen_failwrapper: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
  }
  


  //DEBUG:
  if(debugfail>=2 && showmessages){
    struct of_state q;
    MYFUN(get_stateforcheckinversion(pr, ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
    int outputtype=inputtype;
    FTYPE Unew[NPR];
    MYFUN(primtoU(outputtype,pr, &q, ptrgeom, Unew),"step_ch.c:advance()", "primtoU()", 1); // UtoU inside doesn't do anything...therefore for REMOVERESTMASSFROMUU==1, Unew[UU] will have rest-mass included
    int pliter,pl;
    PLOOP(pliter,pl) dualfprintf(fail_file,"COMPARE: pl=%d pr=%g U=%g Unew=%g\n",pl,pr[pl],U[pl],Unew[pl]);
    int jj;
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d uradcon=%g uradcov=%g\n",jj,q.uradcon[jj],q.uradcov[jj]);
    DLOOPA(jj) dualfprintf(fail_file,"COMPARE: jj=%d ucon=%g ucov=%g\n",jj,q.ucon[jj],q.ucov[jj]);
  }

  //DEBUG:
  if(showmessages || debugfail>=2){
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










// get dt for explicit sub-cyclings
static void get_dtsub(int method, FTYPE *pr, struct of_state *q, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother,  FTYPE *CUf, FTYPE *Gdpl, FTYPE chi, FTYPE *Gdabspl, struct of_geom *ptrgeom, FTYPE *dtsub)
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
    PLOOP(pliter,pl) Gddtpl[pl] = Gdabspl[pl]*realdt;
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
      eomtypelocal=*eomtype; // re-default
      int failutoprim=Utoprimgen_failwrapper(doradonly,showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, EVOLVEUTOPRIM, UNOTHING, Unew, ptrgeom, prnew, &newtonstats);

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
    PLOOP(pliter,pl) if(!isfinitel(Gpl[pl])) return(EXPLICITFAILED);


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
    int failutoprim=Utoprimgen_failwrapper(doradonly,showmessages, allowlocalfailurefixandnoreport, finalstep, &eomtypelocal, EVOLVEUTOPRIM, UNOTHING, Unew, ptrgeom, prnew, &newtonstats);
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
int koral_source_rad(int whichradsourcemethod, FTYPE *piin, FTYPE *pb, FTYPE *pf, int *didreturnpf, int *eomtype, FTYPE *Uiin, FTYPE *Ufin, FTYPE *CUf, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *dUother, FTYPE (*dUcomp)[NPR])
{
  int pliter,pl;
  int showmessages=0; // 0 ok if not debugging and think everything works.
  int showmessagesheavy=0;

  //  if(1||nstep>=800){
  //    showmessages=showmessagesheavy=1;
  //  }


  // make code use "orig" values that below can modify, so return preserves no changes to these no matter what internal functions do.
  // save pb, Uiin, Ufin, and q to avoid being modified upon return
  FTYPE pborig[NPR],piinorig[NPR],Uiinorig[NPR],Ufinorig[NPR];
  struct of_state qorigmem;
  struct of_state *qorig=&qorigmem;
  PLOOP(pliter,pl){
    pborig[pl]=pb[pl];
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
  // IMPLICIT TYPES
  //
  /////////////////
  else if(whichradsourcemethod==SOURCEMETHODIMPLICIT){
 
    int failimplicit=koral_source_rad_implicit(eomtype, pborig, piinorig, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);
 
    if(failimplicit){
      if(IMPLICITREVERTEXPLICIT){ // single level recusive call (to avoid duplicate confusing code)
        // assume if revert from implicit, then need to do sub-cycles
        int failexplicit=koral_source_rad(SOURCEMETHODEXPLICITSUBCYCLEREVERSIONFROMIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);
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
        if(debugfail>=2) dualfprintf(fail_file,"BAD: implicit failed and didn't choose to revert: ijk=%d %d %d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k);
        GLOBALMACP0A1(pflag,ptrgeom->i,ptrgeom->j,ptrgeom->k,FLAGUTOPRIMRADFAIL)=UTOPRIMRADFAILCASE2B; // must set as failure in case can fixup.
        return(0);        
      }
    }// end if failed to do implicit

    // no failure in implicit, then just return
    // and if did implicit, then better pf guess
    PLOOP(pliter,pl) pf[pl]=pborig[pl];
    *didreturnpf=1;
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
    int failreturn=koral_source_rad(SOURCEMETHODEXPLICITCHECKSFROMIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp);

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
      return(koral_source_rad(SOURCEMETHODIMPLICIT, piinorig, pborig, pf, didreturnpf, eomtype, Uiinorig, Ufinorig, CUf, ptrgeom, qorig, dUother, dUcomp));
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
void koral_source_rad_calc(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *Gdpl, FTYPE *Gdabspl, FTYPE *chi, FTYPE *Tgas, struct of_state *q)
{
  int jj;
  int pliter,pl;
  FTYPE Gd[NDIM],Gdabs[NDIM];
  struct of_state qlocal;
  FTYPE chilocal,Tgaslocal;

  if(q==NULL) q=&qlocal;
  if(chi==NULL) chi=&chilocal;
  if(Tgas==NULL) Tgas=&Tgaslocal;


  // no, thermodynamics stuff can change since MHD fluid U changes, so must do get_state() as above
  //  get_state_uconucovonly(pr, ptrgeom, q);
  //  get_state_uradconuradcovonly(pr, ptrgeom, q);
  get_state(pr,ptrgeom,q);

  calc_Gd(pr, ptrgeom, q, Gd, Tgas, chi, Gdabs);

  PLOOP(pliter,pl) Gdpl[pl] = 0.0;
  // equal and opposite forces on fluid and radiation due to radiation 4-force
  // sign of G that goes between Koral determination of G and HARM source term (e.g. positive \lambda is a cooling of the fluid and heating of the photons, and gives G_t>0 so -G_t<0 and adds to R^t_t such that R^t_t - G_t becomes more negative and so more photon energy density)
  // That is, equation is d_t R^t_t + Gdpl = 0
#define SIGNGD (1.0)
  // keep SIGNGD as 1.0.  Just apply SIGNGD2 in front of Gdpl in other places.
  // sign fixed-linked as + for URAD0 case.
  DLOOPA(jj) Gdpl[UU+jj]     = -SIGNGD*Gd[jj];
  DLOOPA(jj) Gdpl[URAD0+jj]  = +SIGNGD*Gd[jj];

  if(Gdabspl!=NULL){
    PLOOP(pliter,pl) Gdabspl[pl] = 0.0;
    DLOOPA(jj) Gdabspl[UU+jj]     = Gdabs[jj];
    DLOOPA(jj) Gdabspl[URAD0+jj]  = Gdabs[jj];
  }

#if(DOENTROPY!=DONOENTROPY && ENTROPY!=-100)
  pl=ENTROPY;
  // The equation is (1/\sqrt{-g})*d_\mu(\sqrt{-g} s\rho_0 u^\mu) + Gdpl[mu].ucon[mu] = 0
  FTYPE Gdplentropycontribs[NDIM];
  // -Gdpl[UU+jj] is so heating (so lowering of T^t_t to be more negative) implies increases entropy.
  // assumes Gpl includes kappa already with rho so that Gpl is energy per unit volume per unit time.  Dividing by T (energy) gives a dimensionless thing (entropy) per unit volume.
  DLOOPA(jj) Gdplentropycontribs[jj] = (-Gdpl[UU+jj])*(q->ucon[jj])/(*Tgas);

  Gdpl[pl] = 0.0;
  DLOOPA(jj) Gdpl[pl] += Gdplentropycontribs[jj];

  if(Gdabspl!=NULL){
    Gdabspl[pl] = 0.0;
    DLOOPA(jj) Gdabspl[pl] += fabs(Gdplentropycontribs[jj]);
  }
#endif



}


// get 4-force and dtsub for all pl's
static void koral_source_dtsub_rad_calc(int method, FTYPE *pr, FTYPE *Ui, FTYPE *Uf, FTYPE *dUother, FTYPE *CUf, FTYPE *Gdpl, struct of_geom *ptrgeom, FTYPE *dtsub)
{
  FTYPE Gdabspl[NPR];
  FTYPE chi,Tgas;
  struct of_state q;

  koral_source_rad_calc(pr,ptrgeom,Gdpl,Gdabspl,&chi,&Tgas,&q);

  if(dtsub!=NULL){
    // then assume expect calculation of dtsub
    get_dtsub(method, pr, &q, Ui, Uf, dUother, CUf, Gdpl, chi, Gdabspl, ptrgeom, dtsub);
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
      mhd_calc_rad( pr, jj, ptrgeom, q, &(radstressdir[jj][0]) );
    }
  }
  else DLOOP(jj,kk) radstressdir[jj][kk]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.
}

// compute radiation stres-energy tensor assuming M1 closure
void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *ptrgeom, struct of_state *q, FTYPE *radstressdir)
{
  int jj;

  // R^{dir}_{jj} radiation stress-energy tensor
  if(EOMRADTYPE==EOMRADEDD){
    // force radiation frame to be fluid frame
    DLOOPA(jj) radstressdir[jj]=THIRD*pr[PRAD0]*(4.0*q->ucon[dir]*q->ucov[jj] + delta(dir,jj));
  }
  else if(EOMRADTYPE==EOMRADM1CLOSURE){
    DLOOPA(jj) radstressdir[jj]=THIRD*pr[PRAD0]*(4.0*q->uradcon[dir]*q->uradcov[jj] + delta(dir,jj));
  }
  else DLOOPA(jj) radstressdir[jj]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.


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
  if(!isfinitel(idet)){
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
  u2p_rad(showmessages, allowlocalfailurefixandnoreport, U, pout, ptrgeomtouse, &lpflag, &lpflagrad);

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
int u2p_rad(int showmessages, int allowlocalfailurefixandnoreport, FTYPE *uu, FTYPE *pin, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
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
  //    if(!isfinitel(pin[PRAD0+jj])){
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


  if( gamma2a<GAMMASMALLLIMIT || !isfinitel(gamma2a) ){
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

        if(!isfinitel(gammatemp)){
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
  int failure=!nonfailure || !isfinitel(gammarel2) || !isfinitel(Erf);

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
        Erf=MAX(MIN(Erf,Erf*(-Avcovorig[TT])/(SMALL+fabs(-Avcov[TT]))),Erf0);
        Erf=MAX(MIN(Erf,Erf*(-Avcov[TT])/(SMALL+fabs(-Avcovorig[TT]))),Erf0);
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
    mathematica_report_check(3, failnum, *lpflagrad, BIG, -1, fakedt, ptrgeom, pp, pp, pp, globalpin, globaluu, globaluu, globaluu, globaluu, fakeCUf, qptr, dUother);
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
  int notfinite=(!isfinitel(Erf) || !isfinitel(urfconrel[1])|| !isfinitel(urfconrel[2])|| !isfinitel(urfconrel[3]));
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
  int failure2=gammarel2<1. || delta<0. || !isfinitel(gammarel2); // NOTE: first failure1 already catches delta<0.



  if(failure1){
    //    if(failure1 && tautotmax<TAUFAILLIMIT){ // works for DBLSHADOW

    FTYPE gammarel=gammamax;
    gammarel2=gammamax*gammamax;

    // get new Erf(gammarel)
    get_m1closure_Erf(ptrgeom, Avcon, gammarel2, &Erf);


    // Check if Erf is too small with gamma->gammamax
    if(Erf<ERADLIMIT || !isfinitel(Erf)){
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


    if(Erf<ERADLIMIT || !isfinitel(Erf)){ // JCM
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

    if(Erf<ERADLIMIT || !isfinitel(Erf)){
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
    if(!isfinitel(Erf) || !isfinitel(gammarel2) || !isfinitel(urfconrel[0])|| !isfinitel(urfconrel[1])|| !isfinitel(urfconrel[2])|| !isfinitel(urfconrel[3]) ){
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

