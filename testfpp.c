/* testfpp.temp.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal gn[16]	/* was [4][4] */, gv[16]	/* was [4][4] */;
} metric_;

#define metric_1 metric_

struct {
    doublereal eps, epsbis, dvmin, tol, uminfac, dlogmax, gammaradceiling;
    integer itermax;
} accuracy_;

#define accuracy_1 accuracy_

union {
    struct {
	doublereal gam, gam1, en, en1, rhou0c, sc, ttcovc[4], rtcovc[4], bbc[
		4], dt, fnumargsharm, fnumresultsharm, whichvelrameshharm, 
		gammamaxrad, eradlimit, toltry, tolallow, arad_code__, 
		okappa_es_code11__, okappa_ff_code11__, okappa_bf_code11__;
    } _1;
    struct {
	doublereal gam, gam1, en, en1, rhou0, s, ttcov[4], rtcov[4], bb[4], 
		dt, fnumargsharm, fnumresultsharm, whichvelrameshharm, 
		gammamaxrad, eradlimit, toltry, tolallow, arad_code__, 
		okappa_es_code11__, okappa_ff_code11__, okappa_bf_code11__;
    } _2;
} conserved_;

#define conserved_1 (conserved_._1)
#define conserved_2 (conserved_._2)

struct {
    doublereal ffkap, eskap, arad, ucon[4], ucov[4], rho, bcon[4], bcov[4], 
	    bsq, tmunu[16]	/* was [4][4] */, e, ehat, urcon[4], urcov[4],
	     rmunu[16]	/* was [4][4] */, tgas, trad, b4pi, gamma, dtau, gcon[
	    4], gcov[4];
} funcmhdd_;

#define funcmhdd_1 funcmhdd_

union {
    struct {
	doublereal ffkap, eskap, arad, ucon[4], ucov[4], rho, bcon[4], bcov[4]
		, bsq, tmunu[16]	/* was [4][4] */, e, ehat, urcon[4], 
		urcov[4], rmunu[16]	/* was [4][4] */, tgas, trad, b4pi, 
		gamma, dtau, gcon[4], gcov[4], u0;
	integer iuerr;
    } _1;
    struct {
	doublereal ffkap, eskap, arad, ucon[4], ucov[4], rho, bcon[4], bcov[4]
		, bsq, tmunu[16]	/* was [4][4] */, e, ehat0, urcon[4], 
		urcov[4], rmunu[16]	/* was [4][4] */, tgas, trad, b4pi, 
		gamma, dtau, gcon[4], gcov[4], u0;
	integer iuerr;
    } _2;
} funcradd_;

#define funcradd_1 (funcradd_._1)
#define funcradd_2 (funcradd_._2)

/* Table of constant values */

 doublereal c_b19 = 1.;
 integer c__4 = 4;
 integer c__3 = 3;
 doublereal c_b57 = 3.5;
 doublereal c_b63 = -3.5;
 doublereal c_b64 = -4.5;
 doublereal c_b66 = -.5;
 doublereal c_b73 = .25;

/* ccccccccccccccccccccccccc */

/* As stand-alone program: */

/*     1) Set PRODUCTION 0, 1, or 2 */
/*     2) To compile only: */
/*            gfortran -Wall -cpp -g -O2 test.f -o test.e */
/*        To compile/run with debug info, do: */
/*            cp ../test.f . ; gfortran -Wall -cpp -g -O2 test.f -o test.e ; ./test.e > test.e.out */
/*        To compile/run with debug info and gdb easy reading, do: */
/*            cp ../test.f . ; gfortran -Wall -cpp -g -O0 test.f -o test.e ; ./test.e > test.e.out */
/*        If crashes, then do: */
/*            gdb ./test.e core */
/*     To test for bad warnings: */
/*     cp ../test.f . ; gfortran -Wall -cpp -g -O2 test.f -o test.e &> test.e.makelog ; grep -v "Unused dummy"
 test.e.makelog | grep -v "Unused variable" | grep -i warning */
/* ccccccccccccccccccccccccc */

/* As harm subroutine */

/* 1) Set PRODUCTION 3 */
/* 2) fpp -P test.f > testfpp.temp.f ; f2c -Wall -f -P testfpp.temp.f ; sed 's///g' testfpp.temp.c > tes
tfpp.c ;sed 's///g' testfpp.temp.P > testfpp.P */
/*    i.e. testfpp.c from f2c: MUST remove  in front of variables! */
/* about pre-processor directives: */
/* http://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html */

/* 0 means normal full output */
/* 1 means very simple Jon's version of err and sol output */
/* 2 means only Jon's err output */
/* 3 means no output or input (harm mode) */
/* #define PRODUCTION 0 */
/* #define PRODUCTION 1 */
/* #define PRODUCTION 2 */
/* 0 : iterate 4-velocities u^i=prim(2,3,4) (radiation is inverted separately and values for u^i_{rad} compute
d separately) */
/* 1 : iterate 3-velocities (not used) */
/* 2 : iterate relative 4-velocities \tilde{u}^i = prim(2,3,4) ("") */
/*     WHICHVEL should be set same as in harm */
/* #define 2 0 */
/* error below which will consider BAD solution and not even compute final solution */
/* #define FAILLIMIT (1D-6) */
/* Choose actual harm choice even if pretty strict. */
/* #define FAILLIMIT (1D-8) */
/* 11 vars, failcode, error, iterations */
/* Subroutine */ int rameshsolver_(doublereal *args, doublereal *resultseng, 
	doublereal *resultsent)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer i_dnnt(doublereal *);

    /* Local variables */
    extern /* Subroutine */ int getfinal_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
     doublereal uradconf[4], uradconi[4], ugasconf[4], ugasconi[4], 
	    uradcovf[4];
     integer iproblem;
     doublereal uradconp[4], uradcovi[4], ugascovf[4], ugasconp[4], 
	    ugascovi[4], errorabs, uradcovp[4], ugascovp[4];
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer showboth, ientropy, i__, j;
    extern /* Subroutine */ int readtype1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , readtype2_(doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *), readtype3_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *);
     doublereal s[5], bconfinal[4], bcovfinal[4], b1[5], b2[5];
     integer idatatype;
    extern /* Subroutine */ int radsource_(doublereal *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), calctmunu_(doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    mhdinvert_(doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *), solveucon_(doublereal *, doublereal *
	    );
     doublereal guesstype, ef, hf[16]	/* was [4][4] */, ep, hp[16]	
	    /* was [4][4] */, sf, uf, si, ui, sp, up, rtcovfinal[4], 
	    ttcovfinal[4], rmunufinal[16]	/* was [4][4] */, tmunufinal[
	    16]	/* was [4][4] */, origintype, bbf[4], bbp[4], r00f, r01f, 
	    t00f, t01f, t02f, t00i;
     integer isc[4];
     doublereal t01i, t02i, t03i, t03f, t00p, t01p, t02p, t03p, r00i, 
	    r00p, r01i, r01p, r02i, r02f, r02p, r03i, r03f, r03p, guessstype;
     integer imhd;
     doublereal bsqf, rhof, bsqi, rhoi;
     integer iter;
     doublereal bsqp, prim[4], rhop;
    extern /* Subroutine */ int calcbconbcov_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
     doublereal uradconfinal[4], ugasconfinal[4], uradcovfinal[4], 
	    ugascovfinal[4];
     integer iflag, jflag;
     doublereal bconf[4], bconi[4], bcovf[4], bconp[4], vradf[4], bcovi[
	    4];
     integer which;
     doublereal vgasf[4], bcovp[4], vradp[4], vgasp[4], error[4];
    extern doublereal mymax_(doublereal *, doublereal *);
     doublereal ttotc[4], rhou0f, rhou0i, rhou0p, rtcovi[4], ttcovi[4], 
	    rmunuf[16]	/* was [4][4] */, tmunuf[16]	/* was [4][4] */, 
	    rmunui[16]	/* was [4][4] */, tmunui[16]	/* was [4][4] */, 
	    rmunup[16]	/* was [4][4] */, tmunup[16]	/* was [4][4] */;
     integer ifinish;
     doublereal primeng[4], problem, priment[4], itertot;

/*     Reads in Jon's error file fails.dat and tests our ideas for */
/*     applying the radiation source term */
/*     variables with 'p' at the end refer to the previous time step, */
/*     those with 'i' at the end refer to the current time step */
/*     pre-radiation, and those with 'f' at the end are for the current */
/*     time step post-radiation. */
/*     variables with 'c' at the end are conserved quantities */
/*     corresponding to the problem currently being solved */
/*     variables with 'final' at the end correspond to the final solution */
/*     inputs */
/*     results */
/*     internals */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     origintype=0 means read from file */
/*     origintype=1 means read from function args as if called by C code using latest number of columns, in wh
ich case don't write anything, just pass back result. */
    /* Parameter adjustments */
    --resultsent;
    --resultseng;
    --args;

    /* Function Body */
    if (TRUE_) {
	origintype = 1.;
    } else {
	origintype = 0.;
    }
/* cccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     idatatype=1 means Jon's old data format with 134 numbers */
/*     idatatype=2 means Jon's new data format with 181 numbers */
/*     idatatype=3 means Jon's new data format with 211+11 numbers */
/*      write (*,*) ' which type data file? old(1) new(2) ' */
/*      read (*,*) idatatype */
    idatatype = 3;
/* cccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     iMHD=0 means don't do initial MHD inversion before radiation */
/*     iMHD=1 means do initial MHD inversion before radiation */
    imhd = 0;
    itertot = 0.;
    iproblem = 0;
/*     Read in data from Jon's datafile */
    for (i__ = 1; i__ <= 1; ++i__) {
/*     If ientropy=0, we will try the energy equation. If it is 1, we */
/*     proceed directly to the entropy equation */
	ientropy = 0;
	if (idatatype == 1) {
	    readtype1_(metric_1.gn, metric_1.gv, &rhof, &rhop, &rhou0i, &
		    rhou0f, &rhou0p, &uf, &up, &t00i, &t00f, &t00p, vgasf, 
		    vgasp, &t01i, &t01f, &t01p, &t02i, &t02f, &t02p, &t03i, &
		    t03f, &t03p, bbf, bbp, &ef, &ep, &r00i, &r00f, &r00p, 
		    vradf, vradp, &r01i, &r01f, &r01p, &r02i, &r02f, &r02p, &
		    r03i, &r03f, &r03p, s, &si, &sf, &sp, uradconf, uradcovf, 
		    ugasconf, ugascovf, uradconp, uradcovp, ugasconp, 
		    ugascovp, &ifinish, &errorabs);
	    conserved_1.fnumargsharm = 222.;
	    conserved_1.fnumresultsharm = 15.;
	    conserved_1.whichvelrameshharm = 2.;
	    conserved_1.gam = 1.3333333333333333;
	    conserved_1.gammamaxrad = 50.;
	    conserved_1.eradlimit = 9.9999999999999999e-300;
	    conserved_1.toltry = 1e-12;
	    conserved_1.tolallow = 1e-8;
	    conserved_1.arad_code__ = 1.18316e17;
	    conserved_1.okappa_es_code11__ = 590799.;
	    conserved_1.okappa_ff_code11__ = 3.46764e-17;
	    conserved_1.okappa_bf_code11__ = 6.93528e-15;
	} else if (idatatype == 2) {
	    readtype2_(metric_1.gn, metric_1.gv, &rhof, &rhop, &rhou0i, &
		    rhou0f, &rhou0p, &uf, &up, &t00i, &t00f, &t00p, vgasf, 
		    vgasp, &t01i, &t01f, &t01p, &t02i, &t02f, &t02p, &t03i, &
		    t03f, &t03p, bbf, bbp, &ef, &ep, &r00i, &r00f, &r00p, 
		    vradf, vradp, &r01i, &r01f, &r01p, &r02i, &r02f, &r02p, &
		    r03i, &r03f, &r03p, s, &si, &sf, &sp, uradconf, uradcovf, 
		    ugasconf, ugascovf, uradconp, uradcovp, ugasconp, 
		    ugascovp, &ifinish, &errorabs);
	    conserved_1.fnumargsharm = 222.;
	    conserved_1.fnumresultsharm = 15.;
	    conserved_1.whichvelrameshharm = 2.;
/*         Gam=4.d0/3.d0 */
	    conserved_1.gammamaxrad = 50.;
	    conserved_1.eradlimit = 9.9999999999999999e-300;
	    conserved_1.toltry = 1e-12;
	    conserved_1.tolallow = 1e-8;
	    conserved_1.arad_code__ = 1.18316e17;
	    conserved_1.okappa_es_code11__ = 590799.;
	    conserved_1.okappa_ff_code11__ = 3.46764e-17;
	    conserved_1.okappa_bf_code11__ = 6.93528e-15;
	} else if (idatatype == 3) {
	    readtype3_(&origintype, &args[1], metric_1.gn, metric_1.gv, &rhof,
		     &rhop, &rhou0i, &rhou0f, &rhou0p, &uf, &up, &t00i, &t00f,
		     &t00p, vgasf, vgasp, &t01i, &t01f, &t01p, &t02i, &t02f, &
		    t02p, &t03i, &t03f, &t03p, bbf, bbp, &ef, &ep, &r00i, &
		    r00f, &r00p, vradf, vradp, &r01i, &r01f, &r01p, &r02i, &
		    r02f, &r02p, &r03i, &r03f, &r03p, s, &si, &sf, &sp, 
		    uradconf, uradcovf, ugasconf, ugascovf, uradconp, 
		    uradcovp, ugasconp, ugascovp, &ifinish, &errorabs);
	}
/*     Initialize constants and read in data from Jon's datafile */
/*     eps is fractional shift in primitives for numerical derivatives */
/*     epsbib is fractional accuracy desired in bisection */
/*     dvmin: if abs(uconmu)<1d-4 then shift used is dvmin*eps */
/*     tol: all fractional errors must be below tol for convergence */
/*     uminfac: minimum u_g allowed is uminfrac*rho */
/*     dlogmax: minimum log() in entropy is log() conserved - dlogmin */
/*     itermax is maximum no. of iterations with u_g, entropy below minimum */
	accuracy_1.eps = 1e-6;
	accuracy_1.epsbis = .001;
	accuracy_1.dvmin = 1e-4;
/*     tol=1.d-10 */
/*     harm-like error */
	accuracy_1.tol = conserved_1.toltry;
	accuracy_1.uminfac = 1e-10;
	accuracy_1.dlogmax = log(2.);
	accuracy_1.itermax = 3;
	accuracy_1.gammaradceiling = conserved_1.gammamaxrad;
/*     Gam=adiabatic index Gamma, Gam1=Gamma-1 */
/*     en=polytropic index = 1/Gamma, en1 = n+1 */
	conserved_1.gam1 = conserved_1.gam - 1.;
	conserved_1.en = 1. / conserved_1.gam1;
	conserved_1.en1 = conserved_1.en + 1.;
	if (ifinish == 1) {
	    goto L10;
	}
	++iproblem;
/*     Calculate b^mu and b_mu for 'p' and 'f' */
	calcbconbcov_(bbp, ugasconp, ugascovp, bconp, bcovp, &bsqp);
	calcbconbcov_(bbf, ugasconf, ugascovf, bconf, bcovf, &bsqf);
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Set up initial guess primitives ('p' quantities from prevous time */
/*     step). Make sure up is reasonable. If not, reset up to */
/*     umin=uminfac*rhop. */
/*      CHOOSE to use previous primitve or harm solution as starting point. */
/*     0 = use original prims (normal mode) */
/*     1 = use harm prims and still do normal Ramesh stages */
/*     2 = use harm prims and go straight to 4D iterations to avoid mis-steps (works to get some entropy solut
ions not otherwise found by Ramesh code, but misses *many* energy solutions even though I'm providing 
very close guess.) */
	guesstype = 0.;
/*         guesstype=2 */
/*         guesstype=1 */
	if (guesstype == 0. || errorabs > 1e-9f) {
	    d__1 = accuracy_1.uminfac * rhop;
	    prim[0] = mymax_(&up, &d__1);
	    if (FALSE_) {
		prim[1] = ugasconp[1];
		prim[2] = ugasconp[2];
		prim[3] = ugasconp[3];
	    } else if (TRUE_) {
		prim[1] = vgasp[1];
		prim[2] = vgasp[2];
		prim[3] = vgasp[3];
	    } else {
		s_stop("", (ftnlen)0);
	    }
	    guessstype = 0.;
	} else {
	    d__1 = accuracy_1.uminfac * rhof;
	    prim[0] = mymax_(&uf, &d__1);
	    if (FALSE_) {
		prim[1] = ugasconf[1];
		prim[2] = ugasconf[2];
		prim[3] = ugasconf[3];
	    } else if (TRUE_) {
		prim[1] = vgasf[1];
		prim[2] = vgasf[2];
		prim[3] = vgasf[3];
	    } else {
		s_stop("", (ftnlen)0);
	    }
	}
/*     put any overrides for guesses here */
/*         prim(1)=1.1081585946596330690836004556281327032809D-8 */
/*         prim(2)=-0.075584957687290353262066821240292555500 */
/*         prim(3)=-0.003678584407270331563545276169512108269 */
/*         prim(4)=0.0923803797665372578214676962826934511966 */
/*     Set up conserved quantities 'c' of the current time step ('i' */
/*     state). The values go into common block /conserved/ */
	conserved_1.rhou0c = rhou0i;
	conserved_1.sc = si;
	conserved_1.ttcovc[0] = t00i;
	conserved_1.ttcovc[1] = t01i;
	conserved_1.ttcovc[2] = t02i;
	conserved_1.ttcovc[3] = t03i;
	conserved_1.rtcovc[0] = r00i;
	conserved_1.rtcovc[1] = r01i;
	conserved_1.rtcovc[2] = r02i;
	conserved_1.rtcovc[3] = r03i;
	conserved_1.bbc[0] = 0.;
	conserved_1.bbc[1] = bbf[1];
	conserved_1.bbc[2] = bbf[2];
	conserved_1.bbc[3] = bbf[3];
	if (imhd == 1) {
/*     Carry out initial MHD inversion of the 'i' state */
	    mhdinvert_(prim, &iflag, &jflag, &ientropy, &itertot, &guesstype);
/*     MHD inversion done. Update all the relevant quantities */
/*     corresponding to the 'i' solution */
/*     Calculate u, u^0 and rho */
	    ui = prim[0];
	    solveucon_(prim, ugasconi);
	    contocov_(ugasconi, ugascovi);
	    rhoi = rhou0i / ugasconi[0];
/*     Calculate b^mu and b_mu, and update T^t_mu and R^t_mu so as to be */
/*     consistent with primitives */
	    for (j = 1; j <= 4; ++j) {
		ttotc[j - 1] = conserved_1.ttcovc[j - 1] + conserved_1.rtcovc[
			j - 1];
	    }
	    calcbconbcov_(conserved_1.bbc, ugasconi, ugascovi, bconi, bcovi, &
		    bsqi);
	    calctmunu_(&rhoi, &ui, &bsqi, &conserved_1.gam, ugasconi, 
		    ugascovi, bconi, bcovi, tmunui);
	    for (j = 1; j <= 4; ++j) {
		ttcovi[j - 1] = tmunui[(j << 2) - 4];
		conserved_1.ttcovc[j - 1] = ttcovi[j - 1];
		rtcovi[j - 1] = ttotc[j - 1] - ttcovi[j - 1];
		conserved_1.rtcovc[j - 1] = rtcovi[j - 1];
	    }
/*     Initial MHD inversion done */
	} else {
/*     No MHD inversion done. Proceed directly to solve the radiation */
/*     problem. */
	}
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Now it is time to apply the implicit radiation source term and */
/*     solve for the post-radiation primitives */
/*     common block /conserved/ already contains the relevant conserved */
/*     quantities for the pre-radiation state. For the primitives, use */
/*     the 'i' state as the initial guess. Then solve the lab frame */
/*     energy equation using Newton-Raphson */
	radsource_(prim, &iter, &iflag, &jflag, &ientropy, &itertot, &
		guesstype, primeng, priment, &resultseng[1], &resultsent[1]);
/*     Radiation inversion done. Calculate relevant quantities */
/*     corresponding to the final solutions */
	showboth = i_dnnt(&resultseng[12]) + i_dnnt(&resultsent[12]);
	if (i_dnnt(&resultseng[12]) == 0) {
	    which = 3;
	    getfinal_(&which, &showboth, &resultseng[1], isc, metric_1.gn, 
		    metric_1.gv, hp, hf, vgasp, vgasf, b1, b2, bbp, bbf, 
		    conserved_1.bbc, vradp, vradf, s, ugasconf, ugascovf, 
		    ugasconp, ugascovp, uradconf, uradcovf, uradconp, 
		    uradcovp, ugasconi, ugascovi, uradconi, uradcovi, bconp, 
		    bcovp, bconf, bcovf, bconi, bcovi, tmunup, tmunuf, rmunup,
		     rmunuf, tmunui, rmunui, ttcovi, rtcovi, 
		    conserved_1.ttcovc, conserved_1.rtcovc, ttotc, primeng, 
		    error, primeng, priment, ugasconfinal, ugascovfinal, 
		    uradconfinal, uradcovfinal, bconfinal, bcovfinal, 
		    ttcovfinal, rtcovfinal, tmunufinal, rmunufinal, &
		    conserved_1.rhou0c, &conserved_1.gam);
/*     prim->primeng */
	}
	if (i_dnnt(&resultsent[12]) == 0) {
	    which = 2;
	    getfinal_(&which, &showboth, &resultsent[1], isc, metric_1.gn, 
		    metric_1.gv, hp, hf, vgasp, vgasf, b1, b2, bbp, bbf, 
		    conserved_1.bbc, vradp, vradf, s, ugasconf, ugascovf, 
		    ugasconp, ugascovp, uradconf, uradcovf, uradconp, 
		    uradcovp, ugasconi, ugascovi, uradconi, uradcovi, bconp, 
		    bcovp, bconf, bcovf, bconi, bcovi, tmunup, tmunuf, rmunup,
		     rmunuf, tmunui, rmunui, ttcovi, rtcovi, 
		    conserved_1.ttcovc, conserved_1.rtcovc, ttotc, priment, 
		    error, primeng, priment, ugasconfinal, ugascovfinal, 
		    uradconfinal, uradcovfinal, bconfinal, bcovfinal, 
		    ttcovfinal, rtcovfinal, tmunufinal, rmunufinal, &
		    conserved_1.rhou0c, &conserved_1.gam);
/*     prim->priment */
	}
/*     enddo below is loop over cases */
    }
L10:
    problem = (doublereal) iproblem;
    return 0;
} /* rameshsolver_ */

/*     Once radiation inversion done, calculate relevant quantities */
/*     corresponding to the final solutions */
/* Subroutine */ int getfinal_(integer *which, integer *showboth, doublereal *
	results, integer *isc, doublereal *gn, doublereal *gv, doublereal *hp,
	 doublereal *hf, doublereal *vgasp, doublereal *vgasf, doublereal *b1,
	 doublereal *b2, doublereal *bbp, doublereal *bbf, doublereal *bbc, 
	doublereal *vradp, doublereal *vradf, doublereal *s, doublereal *
	ugasconf, doublereal *ugascovf, doublereal *ugasconp, doublereal *
	ugascovp, doublereal *uradconf, doublereal *uradcovf, doublereal *
	uradconp, doublereal *uradcovp, doublereal *ugasconi, doublereal *
	ugascovi, doublereal *uradconi, doublereal *uradcovi, doublereal *
	bconp, doublereal *bcovp, doublereal *bconf, doublereal *bcovf, 
	doublereal *bconi, doublereal *bcovi, doublereal *tmunup, doublereal *
	tmunuf, doublereal *rmunup, doublereal *rmunuf, doublereal *tmunui, 
	doublereal *rmunui, doublereal *ttcovi, doublereal *rtcovi, 
	doublereal *ttcovc, doublereal *rtcovc, doublereal *ttotc, doublereal 
	*prim, doublereal *error, doublereal *primeng, doublereal *priment, 
	doublereal *ugasconfinal, doublereal *ugascovfinal, doublereal *
	uradconfinal, doublereal *uradcovfinal, doublereal *bconfinal, 
	doublereal *bcovfinal, doublereal *ttcovfinal, doublereal *rtcovfinal,
	 doublereal *tmunufinal, doublereal *rmunufinal, doublereal *rhou0c, 
	doublereal *gam)
{
    /* Builtin functions */
    double sqrt(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
     doublereal gammarad, gammagas, bsqfinal, rhofinal;
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer j, radinvmod;
    extern /* Subroutine */ int calctmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), solveucon_(doublereal *
	    , doublereal *);
     doublereal rhou0final;
    extern /* Subroutine */ int rmunuinvert_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , calcbconbcov_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
     doublereal alpha, turadconfinal[4], tugasconfinal[4], efinal, 
	    ufinal;
    extern /* Subroutine */ int uconrel_(doublereal *, doublereal *);

/*     prim that's used */
    /* Parameter adjustments */
    rmunufinal -= 5;
    tmunufinal -= 5;
    --rtcovfinal;
    --ttcovfinal;
    --bcovfinal;
    --bconfinal;
    --uradcovfinal;
    --uradconfinal;
    --ugascovfinal;
    --ugasconfinal;
    --priment;
    --primeng;
    --error;
    --prim;
    --ttotc;
    --rtcovc;
    --ttcovc;
    --rtcovi;
    --ttcovi;
    rmunui -= 5;
    tmunui -= 5;
    rmunuf -= 5;
    rmunup -= 5;
    tmunuf -= 5;
    tmunup -= 5;
    --bcovi;
    --bconi;
    --bcovf;
    --bconf;
    --bcovp;
    --bconp;
    --uradcovi;
    --uradconi;
    --ugascovi;
    --ugasconi;
    --uradcovp;
    --uradconp;
    --uradcovf;
    --uradconf;
    --ugascovp;
    --ugasconp;
    --ugascovf;
    --ugasconf;
    --s;
    --vradf;
    --vradp;
    --bbc;
    --bbf;
    --bbp;
    --b2;
    --b1;
    --vgasf;
    --vgasp;
    hf -= 5;
    hp -= 5;
    gv -= 5;
    gn -= 5;
    --isc;
    --results;

    /* Function Body */
    ufinal = prim[1];
    solveucon_(&prim[1], &ugasconfinal[1]);
    contocov_(&ugasconfinal[1], &ugascovfinal[1]);
    rhou0final = *rhou0c;
    rhofinal = *rhou0c / ugasconfinal[1];
/*     Update T^mu_nu so as to be consistent with the final */
/*     primitives. Adjust R^t_mu so that the total energy and momentum */
/*     density are unchaned, then calculate the full R^mu_nu */
    for (j = 1; j <= 4; ++j) {
	ttotc[j] = ttcovc[j] + rtcovc[j];
    }
    calcbconbcov_(&bbc[1], &ugasconfinal[1], &ugascovfinal[1], &bconfinal[1], 
	    &bcovfinal[1], &bsqfinal);
    calctmunu_(&rhofinal, &ufinal, &bsqfinal, gam, &ugasconfinal[1], &
	    ugascovfinal[1], &bconfinal[1], &bcovfinal[1], &tmunufinal[5]);
/*     whether to fix cons or not */
    for (j = 1; j <= 4; ++j) {
	ttcovfinal[j] = tmunufinal[(j << 2) + 1];
	rtcovfinal[j] = ttotc[j] - ttcovfinal[j];
    }
    rmunuinvert_(&rtcovfinal[1], &efinal, &uradconfinal[1], &uradcovfinal[1], 
	    &ugasconfinal[1], &rmunufinal[5], &radinvmod);
    alpha = 1. / sqrt(-gn[5]);
    gammagas = ugasconfinal[1] * alpha;
    gammarad = uradconfinal[1] * alpha;
    if (TRUE_) {
/*     get \tilde{u}^\mu_{\rm rad} */
	uconrel_(&uradconfinal[1], turadconfinal);
/*     get \tilde{u}^\mu_{\rm gas} */
	uconrel_(&ugasconfinal[1], tugasconfinal);
    }
/*     HARM order */
    results[1] = rhofinal;
    results[2] = ufinal;
    if (FALSE_) {
	results[3] = ugasconfinal[1];
	results[4] = ugasconfinal[2];
	results[5] = ugasconfinal[3];
	results[6] = ugasconfinal[4];
    } else if (TRUE_) {
	results[3] = tugasconfinal[0];
	results[4] = tugasconfinal[1];
	results[5] = tugasconfinal[2];
	results[6] = tugasconfinal[3];
    } else {
	s_stop("", (ftnlen)0);
    }
    results[7] = efinal;
    if (FALSE_) {
	results[8] = uradconfinal[1];
	results[9] = uradconfinal[2];
	results[10] = uradconfinal[3];
	results[11] = uradconfinal[4];
    } else if (TRUE_) {
	results[8] = turadconfinal[0];
	results[9] = turadconfinal[1];
	results[10] = turadconfinal[2];
	results[11] = turadconfinal[3];
    } else {
	s_stop("", (ftnlen)0);
    }
    results[15] = (doublereal) radinvmod;
/*     DEBUG TEST */
/*      write (13,*) 'RAMESH2: ' */
/*     &     ,results(1),results(2),results(3),results(4),results(5) */
/*     &     ,results(6),results(7),results(8),results(9),results(10) */
/*     &     ,results(11) */
    return 0;
} /* getfinal_ */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int readtype1_(doublereal *gn, doublereal *gv, doublereal *
	rhof, doublereal *rhop, doublereal *rhou0i, doublereal *rhou0f, 
	doublereal *rhou0p, doublereal *uf, doublereal *up, doublereal *t00i, 
	doublereal *t00f, doublereal *t00p, doublereal *vgasf, doublereal *
	vgasp, doublereal *t01i, doublereal *t01f, doublereal *t01p, 
	doublereal *t02i, doublereal *t02f, doublereal *t02p, doublereal *
	t03i, doublereal *t03f, doublereal *t03p, doublereal *bbf, doublereal 
	*bbp, doublereal *ef, doublereal *ep, doublereal *r00i, doublereal *
	r00f, doublereal *r00p, doublereal *vradf, doublereal *vradp, 
	doublereal *r01i, doublereal *r01f, doublereal *r01p, doublereal *
	r02i, doublereal *r02f, doublereal *r02p, doublereal *r03i, 
	doublereal *r03f, doublereal *r03p, doublereal *s, doublereal *si, 
	doublereal *sf, doublereal *sp, doublereal *uradconf, doublereal *
	uradcovf, doublereal *ugasconf, doublereal *ugascovf, doublereal *
	uradconp, doublereal *uradcovp, doublereal *ugasconp, doublereal *
	ugascovp, integer *ifinish, doublereal *errorabs)
{
     doublereal tolallow, errorabsbestexternal;

/*     Read in data in Jon's old format (134 numbers) */
    /* Parameter adjustments */
    --ugascovp;
    --ugasconp;
    --uradcovp;
    --uradconp;
    --ugascovf;
    --ugasconf;
    --uradcovf;
    --uradconf;
    --s;
    --vradp;
    --vradf;
    --bbp;
    --bbf;
    --vgasp;
    --vgasf;
    gv -= 5;
    gn -= 5;

    /* Function Body */
    *ifinish = 0;
    *errorabs = 1.f;
    errorabsbestexternal = 1.f;
    tolallow = 1e-8f;
/*      set \tilde{u}^t=0 */
    vgasf[1] = 0.;
    vgasp[1] = 0.;
    vradf[1] = 0.;
    vradp[1] = 0.;
    return 0;
/* L10: */
    *ifinish = 1;
    return 0;
} /* readtype1_ */

/* Subroutine */ int readtype2_(doublereal *gn, doublereal *gv, doublereal *
	rhof, doublereal *rhop, doublereal *rhou0i, doublereal *rhou0f, 
	doublereal *rhou0p, doublereal *uf, doublereal *up, doublereal *t00i, 
	doublereal *t00f, doublereal *t00p, doublereal *vgasf, doublereal *
	vgasp, doublereal *t01i, doublereal *t01f, doublereal *t01p, 
	doublereal *t02i, doublereal *t02f, doublereal *t02p, doublereal *
	t03i, doublereal *t03f, doublereal *t03p, doublereal *bbf, doublereal 
	*bbp, doublereal *ef, doublereal *ep, doublereal *r00i, doublereal *
	r00f, doublereal *r00p, doublereal *vradf, doublereal *vradp, 
	doublereal *r01i, doublereal *r01f, doublereal *r01p, doublereal *
	r02i, doublereal *r02f, doublereal *r02p, doublereal *r03i, 
	doublereal *r03f, doublereal *r03p, doublereal *s, doublereal *si, 
	doublereal *sf, doublereal *sp, doublereal *uradconf, doublereal *
	uradcovf, doublereal *ugasconf, doublereal *ugascovf, doublereal *
	uradconp, doublereal *uradcovp, doublereal *ugasconp, doublereal *
	ugascovp, integer *ifinish, doublereal *errorabs)
{
     integer eomtype;

/*     Read in data in Jon's new format (181 numbers) */
    /* Parameter adjustments */
    --ugascovp;
    --ugasconp;
    --uradcovp;
    --uradconp;
    --ugascovf;
    --ugasconf;
    --uradcovf;
    --uradconf;
    --s;
    --vradp;
    --vradf;
    --bbp;
    --bbf;
    --vgasp;
    --vgasf;
    gv -= 5;
    gn -= 5;

    /* Function Body */
    *ifinish = 0;
    eomtype = 0;
    conserved_1.tolallow = 1e-8f;
/*      set \tilde{u}^t=0 */
    vgasf[1] = 0.;
    vgasp[1] = 0.;
    vradf[1] = 0.;
    vradp[1] = 0.;
    if (*errorabs < conserved_1.tolallow) {
    } else {
    }
    return 0;
/* L10: */
    *ifinish = 1;
    return 0;
} /* readtype2_ */

/* Subroutine */ int readtype3_(doublereal *origintype, doublereal *args, 
	doublereal *gn, doublereal *gv, doublereal *rhof, doublereal *rhop, 
	doublereal *rhou0i, doublereal *rhou0f, doublereal *rhou0p, 
	doublereal *uf, doublereal *up, doublereal *t00i, doublereal *t00f, 
	doublereal *t00p, doublereal *vgasf, doublereal *vgasp, doublereal *
	t01i, doublereal *t01f, doublereal *t01p, doublereal *t02i, 
	doublereal *t02f, doublereal *t02p, doublereal *t03i, doublereal *
	t03f, doublereal *t03p, doublereal *bbf, doublereal *bbp, doublereal *
	ef, doublereal *ep, doublereal *r00i, doublereal *r00f, doublereal *
	r00p, doublereal *vradf, doublereal *vradp, doublereal *r01i, 
	doublereal *r01f, doublereal *r01p, doublereal *r02i, doublereal *
	r02f, doublereal *r02p, doublereal *r03i, doublereal *r03f, 
	doublereal *r03p, doublereal *s, doublereal *si, doublereal *sf, 
	doublereal *sp, doublereal *uradconf, doublereal *uradcovf, 
	doublereal *ugasconf, doublereal *ugascovf, doublereal *uradconp, 
	doublereal *uradcovp, doublereal *ugasconp, doublereal *ugascovp, 
	integer *ifinish, doublereal *errorabs)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
     doublereal uradconb[4], ugasconb[4], uradcovb[4];
     integer itermode;
     doublereal ugascovb[4], failtype, steppart, eb;
     integer na;
     doublereal ub;
     integer totaliters;
     doublereal scr, src, rhob;
     integer myid, iters, nstep;
     doublereal errorabsbestexternal, gotfirstnofail, failnum;
     integer eomtype;

/*     Read in data in Jon's new format (211+11 numbers) */
    /* Parameter adjustments */
    --ugascovp;
    --ugasconp;
    --uradcovp;
    --uradconp;
    --ugascovf;
    --ugasconf;
    --uradcovf;
    --uradconf;
    --s;
    --vradp;
    --vradf;
    --bbp;
    --bbf;
    --vgasp;
    --vgasf;
    gv -= 5;
    gn -= 5;
    --args;

    /* Function Body */
    *ifinish = 0;
/*     then pull from array of doubles from args */
    na = 0;

    ++na;
    conserved_1.fnumargsharm = args[na];
    ++na;
    conserved_1.fnumresultsharm = args[na];
    ++na;
    conserved_1.whichvelrameshharm = args[na];
    ++na;

    failtype = args[na];
    ++na;
    myid = (integer) args[na];
    ++na;
    failnum = args[na];
    ++na;
    gotfirstnofail = args[na];
    ++na;
    eomtype = (integer) args[na];
    ++na;
    itermode = (integer) args[na];
    ++na;
    *errorabs = args[na];
    ++na;
    errorabsbestexternal = args[na];
    ++na;
    iters = (integer) args[na];
    ++na;
    totaliters = (integer) args[na];
    ++na;
    conserved_1.dt = args[na];
    ++na;
    nstep = (integer) args[na];
    ++na;
    steppart = args[na];
    ++na;
    conserved_1.gam = args[na];
    ++na;

    conserved_1.gammamaxrad = args[na];
    ++na;
    conserved_1.eradlimit = args[na];
    ++na;
    conserved_1.toltry = args[na];
    ++na;
    conserved_1.tolallow = args[na];
    ++na;
    conserved_1.arad_code__ = args[na];
    ++na;
    conserved_1.okappa_es_code11__ = args[na];
    ++na;
    conserved_1.okappa_ff_code11__ = args[na];
    ++na;
    conserved_1.okappa_bf_code11__ = args[na];

    ++na;
    gn[5] = args[na];
    ++na;
    gn[9] = args[na];
    ++na;
    gn[13] = args[na];
    ++na;
    gn[17] = args[na];
    ++na;
    gn[6] = args[na];
    ++na;
    gn[10] = args[na];
    ++na;
    gn[14] = args[na];
    ++na;
    gn[18] = args[na];
    ++na;
    gn[7] = args[na];
    ++na;
    gn[11] = args[na];
    ++na;
    gn[15] = args[na];
    ++na;
    gn[19] = args[na];
    ++na;
    gn[8] = args[na];
    ++na;
    gn[12] = args[na];
    ++na;
    gn[16] = args[na];
    ++na;
    gn[20] = args[na];
    ++na;
    gv[5] = args[na];
    ++na;
    gv[9] = args[na];
    ++na;
    gv[13] = args[na];
    ++na;
    gv[17] = args[na];
    ++na;
    gv[6] = args[na];
    ++na;
    gv[10] = args[na];
    ++na;
    gv[14] = args[na];
    ++na;
    gv[18] = args[na];
    ++na;
    gv[7] = args[na];
    ++na;
    gv[11] = args[na];
    ++na;
    gv[15] = args[na];
    ++na;
    gv[19] = args[na];
    ++na;
    gv[8] = args[na];
    ++na;
    gv[12] = args[na];
    ++na;
    gv[16] = args[na];
    ++na;
    gv[20] = args[na];
    ++na;
    *rhof = args[na];
    ++na;
    scr = args[na];
    ++na;
    rhob = args[na];
    ++na;
    *rhop = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *rhou0i = args[na];
    ++na;
    *rhou0f = args[na];
    ++na;
    *rhou0p = args[na];
    ++na;
    *uf = args[na];
    ++na;
    scr = args[na];
    ++na;
    ub = args[na];
    ++na;
    *up = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *t00i = args[na];
    ++na;
    *t00f = args[na];
    ++na;
    *t00p = args[na];
    ++na;
    vgasf[2] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vgasp[2] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *t01i = args[na];
    ++na;
    *t01f = args[na];
    ++na;
    *t01p = args[na];
    ++na;
    vgasf[3] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vgasp[3] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *t02i = args[na];
    ++na;
    *t02f = args[na];
    ++na;
    *t02p = args[na];
    ++na;
    vgasf[4] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vgasp[4] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *t03i = args[na];
    ++na;
    *t03f = args[na];
    ++na;
    *t03p = args[na];
    ++na;
    bbf[2] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    bbp[2] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    bbf[3] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    bbp[3] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    bbf[4] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    bbp[4] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    *ef = args[na];
    ++na;
    scr = args[na];
    ++na;
    eb = args[na];
    ++na;
    *ep = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *r00i = args[na];
    ++na;
    *r00f = args[na];
    ++na;
    *r00p = args[na];
    ++na;
    vradf[2] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vradp[2] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *r01i = args[na];
    ++na;
    *r01f = args[na];
    ++na;
    *r01p = args[na];
    ++na;
    vradf[3] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vradp[3] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *r02i = args[na];
    ++na;
    *r02f = args[na];
    ++na;
    *r02p = args[na];
    ++na;
    vradf[4] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    vradp[4] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *r03i = args[na];
    ++na;
    *r03f = args[na];
    ++na;
    *r03p = args[na];
    ++na;
    s[1] = args[na];
    ++na;
    scr = args[na];
    ++na;
    scr = args[na];
    ++na;
    s[2] = args[na];
    ++na;
    src = args[na];
    ++na;
    src = args[na];
    ++na;
    *si = args[na];
    ++na;
    *sf = args[na];
    ++na;
    *sp = args[na];
/* con/cov in 1,2,3,4 */
    ++na;
    uradconf[1] = args[na];
    ++na;
    uradcovf[1] = args[na];
    ++na;
    uradconf[2] = args[na];
    ++na;
    uradcovf[2] = args[na];
    ++na;
    uradconf[3] = args[na];
    ++na;
    uradcovf[3] = args[na];
    ++na;
    uradconf[4] = args[na];
    ++na;
    uradcovf[4] = args[na];

    ++na;
    ugasconf[1] = args[na];
    ++na;
    ugascovf[1] = args[na];
    ++na;
    ugasconf[2] = args[na];
    ++na;
    ugascovf[2] = args[na];
    ++na;
    ugasconf[3] = args[na];
    ++na;
    ugascovf[3] = args[na];
    ++na;
    ugasconf[4] = args[na];
    ++na;
    ugascovf[4] = args[na];

    ++na;
    uradconb[0] = args[na];
    ++na;
    uradcovb[0] = args[na];
    ++na;
    uradconb[1] = args[na];
    ++na;
    uradcovb[1] = args[na];
    ++na;
    uradconb[2] = args[na];
    ++na;
    uradcovb[2] = args[na];
    ++na;
    uradconb[3] = args[na];
    ++na;
    uradcovb[3] = args[na];

    ++na;
    ugasconb[0] = args[na];
    ++na;
    ugascovb[0] = args[na];
    ++na;
    ugasconb[1] = args[na];
    ++na;
    ugascovb[1] = args[na];
    ++na;
    ugasconb[2] = args[na];
    ++na;
    ugascovb[2] = args[na];
    ++na;
    ugasconb[3] = args[na];
    ++na;
    ugascovb[3] = args[na];

    ++na;
    uradconp[1] = args[na];
    ++na;
    uradcovp[1] = args[na];
    ++na;
    uradconp[2] = args[na];
    ++na;
    uradcovp[2] = args[na];
    ++na;
    uradconp[3] = args[na];
    ++na;
    uradcovp[3] = args[na];
    ++na;
    uradconp[4] = args[na];
    ++na;
    uradcovp[4] = args[na];

    ++na;
    ugasconp[1] = args[na];
    ++na;
    ugascovp[1] = args[na];
    ++na;
    ugasconp[2] = args[na];
    ++na;
    ugascovp[2] = args[na];
    ++na;
    ugasconp[3] = args[na];
    ++na;
    ugascovp[3] = args[na];
    ++na;
    ugasconp[4] = args[na];
    ++na;
    ugascovp[4] = args[na];

    if (na != 222) {
	s_stop("", (ftnlen)0);
    }
/*      set \tilde{u}^t=0 */
    vgasf[1] = 0.;
    vgasp[1] = 0.;
    vradf[1] = 0.;
    vradp[1] = 0.;
/*     HARMEOMTYPE=2 is entropy, 3 is energy.  If 3 fails, could have reverted to entropy. */
/*     If itermode=0, then in default harm mode, this was not the last attempt for the solver, so just looking
 at why this strategy failed -- not failure of harm ultimately. */
/*     If itermode=1 shows up, then check whether PREVBESTHARMERR was ok/good enough even if not <tol.  If ite
rmode=1 and both current and best error is bad, actually BAD case for harm. */
    if (*errorabs < conserved_1.tolallow) {
    } else {
    }
    return 0;
/* L10: */
    *ifinish = 1;
    return 0;
} /* readtype3_ */

doublereal mymax_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern integer myisnan_(doublereal *);

    ret_val = max(*x,*y);
    if (*x != *x) {
	ret_val = *x;
    }
    if (*y != *y) {
	ret_val = *y;
    }
    if (myisnan_(x) == 1) {
	ret_val = *x;
    }
    if (myisnan_(y) == 1) {
	ret_val = *y;
    }
    return ret_val;
} /* mymax_ */

doublereal myabs_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern integer myisnan_(doublereal *);

    ret_val = abs(*x);
    if (myisnan_(x) == 1) {
	ret_val = *x;
    }
    return ret_val;
} /* myabs_ */

integer myisnan_(doublereal *x)
{
    /* System generated locals */
    integer ret_val;

/* for C / harm */
    if (*x != *x) {
	ret_val = 1;
    } else {
	ret_val = 0;
    }
    return ret_val;
} /* myisnan_ */

doublereal mydiv_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal myabs_(doublereal *);
    extern integer myisnan_(doublereal *);

/*     Avoid division by zero */
    ret_val = *x * d_sign(&c_b19, y) / (myabs_(y) + 1e-300);
    if (*x != *x) {
	ret_val = *x;
    }
    if (*y != *y) {
	ret_val = *y;
    }
    if (myisnan_(x) == 1) {
	ret_val = *x;
    }
    if (myisnan_(y) == 1) {
	ret_val = *y;
    }
    return ret_val;
} /* mydiv_ */

/* Subroutine */ int solveucon_(doublereal *prim, doublereal *con)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
     integer j;
    extern /* Subroutine */ int solveu0old_(doublereal *), solveu0new_(
	    doublereal *, doublereal *);
     doublereal tcon[4];

/*     Calculates full u^\mu given the primitives */
    /* Parameter adjustments */
    --con;
    --prim;

    /* Function Body */
    if (FALSE_) {
	for (j = 2; j <= 4; ++j) {
	    con[j] = prim[j];
	}
	solveu0old_(&con[1]);
    } else if (TRUE_) {
	tcon[0] = 0.;
	for (j = 2; j <= 4; ++j) {
	    tcon[j - 1] = prim[j];
	}
	solveu0new_(tcon, &con[1]);
    } else {
	s_stop("", (ftnlen)0);
    }
    return 0;
} /* solveucon_ */

/* Subroutine */ int uconrel_(doublereal *con, doublereal *tcon)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
     integer j;
     doublereal alphalapse;

/*     Calculates \tilde{u}^i from full u^\mu */
/*     \alpha = 1/sqrt(-g^{tt}) */
    /* Parameter adjustments */
    --tcon;
    --con;

    /* Function Body */
    alphalapse = 1.f / sqrt((d__1 = -metric_1.gn[0], abs(d__1)));
    tcon[1] = 0.;
    for (j = 2; j <= 4; ++j) {
	tcon[j] = con[j] + con[1] * (metric_1.gn[(j << 2) - 4] * alphalapse * 
		alphalapse);
    }
    return 0;
} /* uconrel_ */

/* Subroutine */ int solveu0old_(doublereal *con)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
     integer j, k;
     doublereal au0, bu0, cu0;

/*     Calculates u^0 given the other three components of the 4-velocity, */
/*     u^1, u^2, u^3 */
    /* Parameter adjustments */
    --con;

    /* Function Body */
    au0 = metric_1.gv[0];
    bu0 = 0.;
    cu0 = 1.;
    for (j = 2; j <= 4; ++j) {
	bu0 += metric_1.gv[(j << 2) - 4] * 2. * con[j];
	for (k = 2; k <= 4; ++k) {
	    cu0 += metric_1.gv[j + (k << 2) - 5] * con[j] * con[k];
	}
    }
    con[1] = (-bu0 - sqrt(bu0 * bu0 - au0 * 4. * cu0)) / (au0 * 2.);
    return 0;
} /* solveu0old_ */

/* Subroutine */ int solveu0new_(doublereal *tcon, doublereal *con)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
     integer j, k;
     doublereal alphalapse, qsq, gamma;

/*     Calculates u^\mu given \tilde{u}^i */
/*      q^2 = \tilde{u}^\mu \tilde{u}^\nu g_{\mu\nu} and \tilde{u}^t=0 */
    /* Parameter adjustments */
    --con;
    --tcon;

    /* Function Body */
    tcon[1] = 0.f;
    qsq = 0.f;
    for (j = 2; j <= 4; ++j) {
	for (k = 2; k <= 4; ++k) {
	    qsq += metric_1.gv[j + (k << 2) - 5] * tcon[j] * tcon[k];
	}
    }
/*     \gamma = \sqrt{1+q^2} */
    gamma = sqrt(qsq + 1.f);
/*     \alpha = 1/sqrt(-g^{tt}) */
    alphalapse = 1.f / sqrt((d__1 = -metric_1.gn[0], abs(d__1)));
/*     u^t = \gamma/\alpha */
    con[1] = gamma / alphalapse;
/*     u^j = \tilde{u}^j - u^t \beta^j */
    for (j = 2; j <= 4; ++j) {
	con[j] = tcon[j] - con[1] * (metric_1.gn[(j << 2) - 4] * alphalapse * 
		alphalapse);
    }
    return 0;
} /* solveu0new_ */

/* Subroutine */ int contocov_(doublereal *vecin, doublereal *vecout)
{
     integer i__, j;

/*     Converts a contravariant 4-vector to a covariant 4-vector */
    /* Parameter adjustments */
    --vecout;
    --vecin;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	vecout[i__] = 0.;
	for (j = 1; j <= 4; ++j) {
	    vecout[i__] += metric_1.gv[i__ + (j << 2) - 5] * vecin[j];
	}
    }
    return 0;
} /* contocov_ */

/* Subroutine */ int covtocon_(doublereal *vecin, doublereal *vecout)
{
     integer i__, j;

/*     Converts a covariant 4-vector to a contravariant 4-vector */
    /* Parameter adjustments */
    --vecout;
    --vecin;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	vecout[i__] = 0.;
	for (j = 1; j <= 4; ++j) {
	    vecout[i__] += metric_1.gn[i__ + (j << 2) - 5] * vecin[j];
	}
    }
    return 0;
} /* covtocon_ */

doublereal condotcov_(doublereal *con, doublereal *cov)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
     integer i__;

/*     Calculates the dot product con^mu cov_mu */
    /* Parameter adjustments */
    --cov;
    --con;

    /* Function Body */
    ret_val = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	ret_val += con[i__] * cov[i__];
    }
    return ret_val;
} /* condotcov_ */

/* Subroutine */ int calcbconbcov_(doublereal *bb, doublereal *ucon, 
	doublereal *ucov, doublereal *bcon, doublereal *bcov, doublereal *bsq)
{
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer i__;
    extern doublereal condotcov_(doublereal *, doublereal *);
     doublereal udotb;

/*     Given B(1), B(2), B(3), and gas 4-velocity, calculates the */
/*     four-vectors b^mu and b_mu, and bsq = b^mub_mu */
    /* Parameter adjustments */
    --bcov;
    --bcon;
    --ucov;
    --ucon;
    --bb;

    /* Function Body */
    bb[1] = 0.;
    udotb = condotcov_(&bb[1], &ucov[1]);
    for (i__ = 1; i__ <= 4; ++i__) {
	bcon[i__] = (bb[i__] + udotb * ucon[i__]) / ucon[1];
    }
    contocov_(&bcon[1], &bcov[1]);
    *bsq = condotcov_(&bcon[1], &bcov[1]);
    return 0;
} /* calcbconbcov_ */

/* Subroutine */ int calctmunu_(doublereal *rho, doublereal *u, doublereal *
	bsq, doublereal *gam, doublereal *ucon, doublereal *ucov, doublereal *
	bcon, doublereal *bcov, doublereal *tmunu)
{
     integer i__, j;
     doublereal p;
    extern /* Subroutine */ int makedelta_(doublereal *);
     doublereal bsq2, delta[16]	/* was [4][4] */;

/*     Given rho, u, 4-velocity, magnetic 4-vector, calculates the gas */
/*     stress energy tensor T^mu_nu */
    /* Parameter adjustments */
    tmunu -= 5;
    --bcov;
    --bcon;
    --ucov;
    --ucon;

    /* Function Body */
    makedelta_(delta);
    p = (*gam - 1.) * *u;
    bsq2 = *bsq * .5;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    if (i__ == 1 && j == 1) {
/*     Replace T^0_0 -> T^0_0 + rho*u^0 */
		tmunu[i__ + (j << 2)] = *rho * ucon[i__] * (ucov[j] + 1.) + (*
			u + p + *bsq) * ucon[i__] * ucov[j] + (p + bsq2) * 
			delta[i__ + (j << 2) - 5] - bcon[i__] * bcov[j];
	    } else {
		tmunu[i__ + (j << 2)] = (*rho + *u + p + *bsq) * ucon[i__] * 
			ucov[j] + (p + bsq2) * delta[i__ + (j << 2) - 5] - 
			bcon[i__] * bcov[j];
	    }
	}
    }
    return 0;
} /* calctmunu_ */

/* Subroutine */ int calcrmunu_(doublereal *e, doublereal *ucon, doublereal *
	ucov, doublereal *rmunu)
{
     integer i__, j;
    extern /* Subroutine */ int makedelta_(doublereal *);
     doublereal one3, four3, delta[16]	/* was [4][4] */;

/*     Given the radiation frame energy density E and 4-velocity, */
/*     calculates the radiation stress-energy tensor R^mu_nu */
    /* Parameter adjustments */
    rmunu -= 5;
    --ucov;
    --ucon;

    /* Function Body */
    makedelta_(delta);
    four3 = 1.3333333333333333;
    one3 = .33333333333333331;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    rmunu[i__ + (j << 2)] = four3 * *e * ucon[i__] * ucov[j] + one3 * 
		    *e * delta[i__ + (j << 2) - 5];
	}
    }
    return 0;
} /* calcrmunu_ */

/* Subroutine */ int makedelta_(doublereal *delta)
{
     integer i__, j;

/*     Calculate Kronecker delta */
    /* Parameter adjustments */
    delta -= 5;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    if (i__ == j) {
		delta[i__ + (j << 2)] = 1.;
	    } else {
		delta[i__ + (j << 2)] = 0.;
	    }
	}
    }
    return 0;
} /* makedelta_ */

/* Subroutine */ int rmunuinvert_(doublereal *rtcov, doublereal *e, 
	doublereal *ucon, doublereal *ucov, doublereal *ugascon, doublereal *
	rmunu, integer *radinvmod)
{
    extern /* Subroutine */ int rmunuinvert1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;

/*      call Rmunuinvert0(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod) */
    /* Parameter adjustments */
    rmunu -= 5;
    --ugascon;
    --ucov;
    --ucon;
    --rtcov;

    /* Function Body */
    rmunuinvert1_(&rtcov[1], e, &ucon[1], &ucov[1], &ugascon[1], &rmunu[5], 
	    radinvmod);
    return 0;
} /* rmunuinvert_ */

/*     Old routine */
/* Subroutine */ int rmunuinvert0_(doublereal *rtcov, doublereal *e, 
	doublereal *ucon, doublereal *ucov, doublereal *ugascon, doublereal *
	rmunu, integer *radinvmod)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int covtocon_(doublereal *, doublereal *), 
	    contocov_(doublereal *, doublereal *);
     integer i__, j;
    extern /* Subroutine */ int calcrmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
     doublereal sum, disc, aquad, bquad, cquad, rtcon[4];

/*     Given the row R^t_mu of the radiation tensor, solves for the */
/*     radiation frame energy density E and 4-velocity and calculates the */
/*     full tensor R^mu_nu */
/*     Convert R^0_mu to R^0^mu */
    /* Parameter adjustments */
    rmunu -= 5;
    --ugascon;
    --ucov;
    --ucon;
    --rtcov;

    /* Function Body */
    covtocon_(&rtcov[1], rtcon);
/*     Set up and solve quadratic equation for E */
/*     default is no ceiling hit on gamma */
    *radinvmod = 0;

    aquad = metric_1.gn[0];
    bquad = rtcon[0] * -2.;
    cquad = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    cquad -= metric_1.gv[i__ + (j << 2) - 5] * rtcon[i__ - 1] * rtcon[
		    j - 1];
	}
    }
    cquad *= 3.;
    disc = bquad * bquad - aquad * 4. * cquad;
/*     Make sure signs are okay for a physical solution. If not, go to 10 */
/*     for alternative calculation */
    if (aquad * cquad >= 0. || disc < 0.) {
	goto L10;
    }
/*     Use negative sign of discriminant and solve for E and ucon */
    *e = cquad * 2. / (-bquad + sqrt(disc));
/*     Make sure E is positive. If not, go to 10 */
    if (*e < 0.) {
	goto L10;
    }
    ucon[1] = sqrt(rtcon[0] * 3. / *e - metric_1.gn[0]) * .5;
    for (i__ = 2; i__ <= 4; ++i__) {
	ucon[i__] = (rtcon[i__ - 1] * 3. - *e * metric_1.gn[(i__ << 2) - 4]) /
		 (*e * 4. * ucon[1]);
    }
/*     Make sure the Lorentz factor is below the ceiling. If not go to 10 */
    if (ucon[1] <= accuracy_1.gammaradceiling * sqrt(-metric_1.gn[0])) {
	goto L20;
    } else {
	goto L10;
    }
/*     This segment is for problem cases. We set gamma_radiation equal to */
/*     its ceiling value and solve for E and u^i without using R^00. */
L10:
    ucon[1] = accuracy_1.gammaradceiling * sqrt(-metric_1.gn[0]);
    *radinvmod = 1;
    aquad = 0.;
    bquad = 0.;
    cquad = metric_1.gv[0] * ucon[1] * ucon[1] + 1.;
    for (i__ = 2; i__ <= 4; ++i__) {
	bquad += metric_1.gv[(i__ << 2) - 4] * 1.5 * rtcon[i__ - 1];
	cquad -= metric_1.gv[(i__ << 2) - 4] * .5 * metric_1.gn[(i__ << 2) - 
		4];
	for (j = 2; j <= 4; ++j) {
	    aquad += metric_1.gv[i__ + (j << 2) - 5] * 9. * rtcon[i__ - 1] * 
		    rtcon[j - 1] / (ucon[1] * 16. * ucon[1]);
	    bquad -= metric_1.gv[i__ + (j << 2) - 5] * 3. * (rtcon[i__ - 1] * 
		    metric_1.gn[(j << 2) - 4] + rtcon[j - 1] * metric_1.gn[(
		    i__ << 2) - 4]) / (ucon[1] * 16. * ucon[1]);
	    cquad += metric_1.gv[i__ + (j << 2) - 5] * metric_1.gn[(i__ << 2) 
		    - 4] * metric_1.gn[(j << 2) - 4] / (ucon[1] * 16. * ucon[
		    1]);
	}
    }
    disc = bquad * bquad - aquad * 4. * cquad;
/*     Check if disc > 0. If not, we have trouble again, and we need to */
/*     use yet another scheme! */
    if (disc < 0.) {
	goto L30;
    }
    *e = (-bquad - sqrt(disc)) / (cquad * 2.);
    for (i__ = 2; i__ <= 4; ++i__) {
	ucon[i__] = (rtcon[i__ - 1] * 3. - *e * metric_1.gn[(i__ << 2) - 4]) /
		 (*e * 4. * ucon[1]);
    }
    sum = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    sum += metric_1.gv[i__ + (j << 2) - 5] * ucon[i__] * ucon[j];
	}
    }
    goto L20;
/*     Third try! What should we do here? */
L30:
/*     Last-ditch effort. Set radiation velocity equal to gas velocity */
    for (i__ = 1; i__ <= 4; ++i__) {
	ucon[i__] = ugascon[i__];
    }
/* Computing 2nd power */
    d__1 = ucon[1];
    *e = rtcon[0] * 3. / (d__1 * d__1 * 4. - metric_1.gn[0]);
/*     Calculate urad_mu */
L20:
    contocov_(&ucon[1], &ucov[1]);
/*     Calculate the full tensor R^mu_nu */
    calcrmunu_(e, &ucon[1], &ucov[1], &rmunu[5]);
    return 0;
} /* rmunuinvert0_ */

/*     New routine */
/* Subroutine */ int rmunuinvert1_(doublereal *rtcov, doublereal *e, 
	doublereal *ucon, doublereal *ucov, doublereal *ugascon, doublereal *
	rmunu, integer *radinvmod)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int covtocon_(doublereal *, doublereal *), 
	    contocov_(doublereal *, doublereal *);
     integer i__, j;
    extern /* Subroutine */ int calcrmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
     doublereal gc, gammaradsqr, disc, rsqr, alpha, aquad, bquad, cquad,
	     rtcon[4], xisqr, etacon[4], etacov[4], xisqrc, alphasq;

/*     Given the row R^t_mu of the radiation tensor, solves for the */
/*     radiation frame 4-velocity u^mu and energy density E and then */
/*     calculates the full tensor R^mu_nu */
/*     default is no ceiling hit on gamma */
    /* Parameter adjustments */
    rmunu -= 5;
    --ugascon;
    --ucov;
    --ucon;
    --rtcov;

    /* Function Body */
    *radinvmod = 0;
/*      arad=1.18316d17 */
/*     Convert R^0_mu to R^0^mu */
    covtocon_(&rtcov[1], rtcon);
/*      write (*,*) ' Rtcov: ',(Rtcov(i),i=1,4) */
/*      write (*,*) ' Rtcon: ',(Rtcon(i),i=1,4) */
/*     Calculate lapse alpha and ZAMO four velocity eta_mu and eta^mu */
    alphasq = 1. / (-metric_1.gn[0]);
    alpha = sqrt(alphasq);
    etacov[0] = -alpha;
    etacov[1] = 0.;
    etacov[2] = 0.;
    etacov[3] = 0.;
    covtocon_(etacov, etacon);
    if (rtcon[0] < 0.) {
/*     Check if R^tt is positive. If not, set energy density to a very */
/*     small value and set 4-velocity equal to ZAMO velocity. */
	ucon[1] = 1. / alpha;
	ucon[2] = etacon[1];
	ucon[3] = etacon[2];
	ucon[4] = etacon[3];
	*e = conserved_1.eradlimit;
/*         E=arad*1.d-36 */
    } else {
/*     Compute |R|^2, xisqr=|R|^2*g^tt/(R^tt)^2, and xisqrceiling.  Check */
/*     which regime we are in and choose appropriate solution. */
	rsqr = 0.;
	for (i__ = 1; i__ <= 4; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		rsqr += metric_1.gv[i__ + (j << 2) - 5] * rtcon[i__ - 1] * 
			rtcon[j - 1];
	    }
	}
/* Computing 2nd power */
	d__1 = rtcon[0];
	xisqr = rsqr * metric_1.gn[0] / (d__1 * d__1);
	gc = accuracy_1.gammaradceiling;
/* Computing 2nd power */
	d__1 = gc;
/* Computing 4th power */
	d__2 = gc, d__2 *= d__2;
/* Computing 2nd power */
	d__3 = gc;
	xisqrc = (d__1 * d__1 * 8. + 1.) / (d__2 * d__2 * 16. - d__3 * d__3 * 
		8. + 1.);
	if (xisqr >= 1.) {
	    *radinvmod = 1;
/*     If xisqr>=1, set gammarad=1, solve for E, and obtain remaining */
/*     urad^i. */
	    ucon[1] = 1. / alpha;
	    ucon[2] = etacon[1];
	    ucon[3] = etacon[2];
	    ucon[4] = etacon[3];
/* Computing 2nd power */
	    d__1 = ucon[1];
	    *e = rtcon[0] * 3. / (d__1 * d__1 * 4. + metric_1.gn[0]);
	} else if (xisqr <= xisqrc) {
/*     If xisqr<=xisqrceiling, set gammarad=gammaradceiling and compute */
/*     the rest of the solution using our old method */
	    *radinvmod = 1;
	    ucon[1] = accuracy_1.gammaradceiling * sqrt(-metric_1.gn[0]);
	    aquad = 0.;
	    bquad = 0.;
	    cquad = metric_1.gv[0] * ucon[1] * ucon[1] + 1.;
	    for (i__ = 2; i__ <= 4; ++i__) {
		bquad += metric_1.gv[(i__ << 2) - 4] * 1.5 * rtcon[i__ - 1];
		cquad -= metric_1.gv[(i__ << 2) - 4] * .5 * metric_1.gn[(i__ 
			<< 2) - 4];
		for (j = 2; j <= 4; ++j) {
		    aquad += metric_1.gv[i__ + (j << 2) - 5] * 9. * rtcon[i__ 
			    - 1] * rtcon[j - 1] / (ucon[1] * 16. * ucon[1]);
		    bquad -= metric_1.gv[i__ + (j << 2) - 5] * 3. * (rtcon[
			    i__ - 1] * metric_1.gn[(j << 2) - 4] + rtcon[j - 
			    1] * metric_1.gn[(i__ << 2) - 4]) / (ucon[1] * 
			    16. * ucon[1]);
		    cquad += metric_1.gv[i__ + (j << 2) - 5] * metric_1.gn[(
			    i__ << 2) - 4] * metric_1.gn[(j << 2) - 4] / (
			    ucon[1] * 16. * ucon[1]);
		}
	    }
	    disc = bquad * bquad - aquad * 4. * cquad;
	    *e = (-bquad - sqrt(disc)) / (cquad * 2.);
	    for (i__ = 2; i__ <= 4; ++i__) {
		ucon[i__] = (rtcon[i__ - 1] * 3. - *e * metric_1.gn[(i__ << 2)
			 - 4]) / (*e * 4. * ucon[1]);
	    }
	} else {
/*     We have a physically valid situation with xisqr within the */
/*     acceptable range. Calculate the solution. */
	    gammaradsqr = (xisqr + 1. + sqrt(xisqr * 3. + 1.)) / (xisqr * 4.);
	    ucon[1] = sqrt(gammaradsqr) / alpha;
/* Computing 2nd power */
	    d__1 = ucon[1];
	    *e = rtcon[0] * 3. / (d__1 * d__1 * 4. + metric_1.gn[0]);
	    for (i__ = 2; i__ <= 4; ++i__) {
		ucon[i__] = (rtcon[i__ - 1] * 3. - *e * metric_1.gn[(i__ << 2)
			 - 4]) / (*e * 4. * ucon[1]);
	    }
	}
    }
/*     Calculate the full tensor R^mu_nu */
    contocov_(&ucon[1], &ucov[1]);
    calcrmunu_(e, &ucon[1], &ucov[1], &rmunu[5]);
    return 0;
} /* rmunuinvert1_ */

/* Subroutine */ int mhdinvert_(doublereal *prim, integer *iflag, integer *
	jflag, integer *ientropy, doublereal *itertot, doublereal *guesstype)
{
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal primsave[4], pressure;
     integer i__, j;
    extern /* Subroutine */ int usolvemhd_(doublereal *, integer *, integer *,
	     integer *, U_fp, U_fp), solveucon_(doublereal *, doublereal *);
     doublereal err4, rhoi;
     integer iter;
     doublereal ucon[4], entr;
    extern /* Subroutine */ int umhd1_();
    extern /* Subroutine */ int newton3_(doublereal *, integer *, integer *, 
	    integer *, U_fp), newton4_(doublereal *, integer *, integer *, 
	    integer *, U_fp, doublereal *);
    extern /* Subroutine */ int funcmhd1_(), funcmhd2_();

/*     Given initial primitives prim(4), solves for primitives that */
/*     satisfy the initial pre-radiation conserved quantities: rhou0c, s, */
/*     Ttcovc(4), Rtcovc(4) */
    /* Parameter adjustments */
    --prim;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	primsave[i__ - 1] = prim[i__];
    }
    if (*guesstype <= 1.) {
/*     First do one round of Newton-Raphson just on the velocities */
	newton3_(&prim[1], &iter, iflag, jflag, (U_fp)funcmhd1_);
	*itertot += iter;
/*     Next work only on u_g using the energy equation */
	usolvemhd_(&prim[1], &iter, iflag, jflag, (U_fp)funcmhd1_, (U_fp)
		umhd1_);
	for (i__ = 1; i__ <= 4; ++i__) {
	    primsave[i__ - 1] = prim[i__];
	}
    }
/*     Carry out the full Newton-Raphson MHD inversion via the */
/*     energy equation */
    err4 = 1e5;
    newton4_(&prim[1], &iter, iflag, jflag, (U_fp)funcmhd1_, &err4);
    *itertot += iter;
/*     Check if the iterations converged. If not, re-solve using the */
/*     entropy equation. */
    if (*iflag == 0 && prim[1] < 0.) {
/*     Check if the internal energy is okay. If not, re-solve using the */
/*     entropy equation. */
	*iflag = 2;
    } else {
/*     Check if entropy looks okay. If not, re-solve using the entropy */
/*     equation. */
	pressure = conserved_1.gam1 * prim[1];
	solveucon_(&prim[1], ucon);
	rhoi = conserved_1.rhou0c / ucon[0];
	entr = log(pow_dd(&pressure, &conserved_1.en) / pow_dd(&rhoi, &
		conserved_1.en1));
	if (entr - conserved_1.sc / conserved_1.rhou0c < -accuracy_1.dlogmax) 
		{
	    *iflag = 3;
	}
    }
    if (*iflag >= 1) {
/*     Restore saved primitives to post-3+1 solution before proceeding to */
/*     solving the entropy equation */
/*     Reset jflag */
	*jflag = 0;
	for (i__ = 1; i__ <= 4; ++i__) {
	    prim[i__] = primsave[i__ - 1];
	}
	*ientropy = 1;
/*     Next work only on u_g using the entropy equation */
/*      call usolveMHD(prim,iter,iflag,jflag,funcMHD2,uMHD2) */
/*     Carry out full Newton-Raphson inversion with the entropy equation */
	err4 = 1e6;
	newton4_(&prim[1], &iter, iflag, jflag, (U_fp)funcmhd2_, &err4);
	*itertot += iter;
	if (*iflag == 9) {
	    for (j = 1; j <= 4; ++j) {
		prim[j] = primsave[j - 1];
	    }
	} else {
	}
	if (*iflag == 0) {
	} else {
	    for (j = 1; j <= 4; ++j) {
		prim[j] = primsave[j - 1];
	    }
	    *iflag = 9;
	}
    }
    return 0;
} /* mhdinvert_ */

/* Subroutine */ int radsource_(doublereal *prim, integer *iter, integer *
	iflag, integer *jflag, integer *ientropy, doublereal *itertot, 
	doublereal *guesstype, doublereal *primeng, doublereal *priment, 
	doublereal *resultseng, doublereal *resultsent)
{
     doublereal primsave[4];
     integer i__;
    extern /* Subroutine */ int usolverad_(doublereal *, integer *, integer *,
	     integer *, U_fp, U_fp), solveucon_(doublereal *, doublereal *);
     doublereal err4;
    extern /* Subroutine */ int uerr1_();
     doublereal erreng, errent;
    extern /* Subroutine */ int newton3_(doublereal *, integer *, integer *, 
	    integer *, U_fp), newton4_(doublereal *, integer *, integer *, 
	    integer *, U_fp, doublereal *);
     integer itereng;
     doublereal ugascon[4], engconv;
     integer iterent;
     doublereal envconv;
    extern /* Subroutine */ int funcrad1_(), funcrad2_();

/*     Given initial primitives prim(4), solves for primitives that */
/*     satisfy the post-radiation equations */
    /* Parameter adjustments */
    --resultsent;
    --resultseng;
    --priment;
    --primeng;
    --prim;

    /* Function Body */
    itereng = 0;
    erreng = 1e3;
    iterent = 0;
    errent = 1e3;
    for (i__ = 1; i__ <= 4; ++i__) {
	primsave[i__ - 1] = prim[i__];
    }
    if (*guesstype <= 1.) {
/*     First do one round of Newton-Raphson just on the velocities */
	newton3_(&prim[1], iter, iflag, jflag, (U_fp)funcrad1_);
	*itertot += *iter;
	for (i__ = 1; i__ <= 4; ++i__) {
	    primsave[i__ - 1] = prim[i__];
	}
	if (*ientropy == 0) {
/*     ientropy=0, so we will first try the energy equation. Initially */
/*     work only on u_g using the energy equation. If ientropy=1, we skip */
/*     all this and go directly to working with the entropy equation. */
	    usolverad_(&prim[1], iter, iflag, jflag, (U_fp)funcrad1_, (U_fp)
		    uerr1_);
	    if (prim[1] > 0.) {
		primsave[0] = prim[1];
	    }
	}
/*   endif guesstype.le.1 */
    }
/*     Carry out full Newton-Raphson on all four primitives using */
/*     Newton-Raphson */
    envconv = 0.;
    err4 = 1e8;
    newton4_(&prim[1], iter, iflag, jflag, (U_fp)funcrad1_, &err4);
    *itertot += *iter;
    itereng = *iter;
    erreng = err4;
/*     default is failed */
    resultseng[12] = 1.;
    if (*iflag == 0 && erreng < conserved_1.tolallow) {
	solveucon_(&prim[1], ugascon);
	if (erreng < conserved_1.tolallow && prim[1] > 0.f && *jflag == 0 && 
		prim[1] == prim[1] && prim[2] == prim[2] && prim[3] == prim[3]
		 && prim[4] == prim[4] && ugascon[0] == ugascon[0]) {
	    resultseng[12] = 0.;
	    engconv = 1.;
	} else if (erreng < conserved_1.tolallow && prim[1] == prim[1] && 
		prim[2] == prim[2] && prim[3] == prim[3] && prim[4] == prim[4]
		 && ugascon[0] == ugascon[0]) {
/*     only bad because negative something */
	    resultseng[12] = 2.;
	    engconv = 1.;
	} else {
/*            treat as bad because error not small or nan'ed out */
	    resultseng[12] = 1.;
	    engconv = 0.;
	}
/*     commenting below to see entropy even if energy converged. */
/*         return */
    } else {
	resultseng[12] = 1.;
	engconv = 0.;
/*     If energy equation does not converge, calculate using the entropy */
/*     equation. */
    }
    resultseng[13] = erreng;
    resultseng[14] = (doublereal) itereng;
    for (i__ = 1; i__ <= 4; ++i__) {
	primeng[i__] = prim[i__];
    }
/*     Restore saved primitives and do Newton-Raphson with the entropy */
/*     equation */
    for (i__ = 1; i__ <= 4; ++i__) {
	prim[i__] = primsave[i__ - 1];
    }
/*     Next work only on u_g using the entropy equation. Since we are */
/*     caurrently using the entropy equation even for the energy equation */
/*     step, we do not need to repeat the 1D search. If and when we */
/*     modify uerr1 to do the proper energy equation, we need this second */
/*     round of 1D search using uerr2. */
/*      call usolverad(prim,iter,iflag,jflag,funcrad2,uerr2) */
/*      default is ent failed */
    resultsent[12] = 1.;
    err4 = 1e8;
    newton4_(&prim[1], iter, iflag, jflag, (U_fp)funcrad2_, &err4);
    *itertot += *iter;
    iterent = *iter;
    errent = err4;
    if (*iflag == 0 && errent < conserved_1.tolallow) {
	solveucon_(&prim[1], ugascon);
	if ((erreng < conserved_1.tolallow || errent < conserved_1.tolallow) 
		&& prim[1] > 0.f && *jflag == 0 && prim[1] == prim[1] && prim[
		2] == prim[2] && prim[3] == prim[3] && prim[4] == prim[4] && 
		ugascon[0] == ugascon[0]) {
	    resultsent[12] = 0.;
	} else if ((erreng < conserved_1.tolallow || errent < 
		conserved_1.tolallow) && prim[1] == prim[1] && prim[2] == 
		prim[2] && prim[3] == prim[3] && prim[4] == prim[4] && 
		ugascon[0] == ugascon[0]) {
	    resultsent[12] = 2.;
	} else {
	    resultsent[12] = 1.;
	}
    } else {
	resultsent[12] = 1.;
    }
    resultsent[13] = errent;
    resultsent[14] = (doublereal) iterent;
    for (i__ = 1; i__ <= 4; ++i__) {
	priment[i__] = prim[i__];
    }
/*     We should never reach this point. Unclear what to do in this case! */
/*     We could keep the solution or restore the saved primitives. */
/*      do i=1,4 */
/*         prim(i)=primsave(i) */
/*      enddo */
    return 0;
} /* radsource_ */

/* Subroutine */ int newton4_(doublereal *prim0, integer *iter, integer *
	iflag, integer *jflag, S_fp func, doublereal *err4)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
     doublereal d__;
     integer i__, j, k;
     doublereal errornorm[4], sum, err3, ajac[16]	/* was [4][4] */;
     integer indx[4];
     doublereal prim[4];
    extern doublereal myabs_(doublereal *);
     doublereal dprim;
     integer niter;
     doublereal error[4];
    extern doublereal mymax_(doublereal *, doublereal *);
     doublereal error0[4];
    extern /* Subroutine */ int lubksb_(doublereal *, integer *, integer *, 
	    integer *, doublereal *), ludcmp_(doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
     doublereal retval, primold[4];

/*     Calculates iteratively by the Newton-Raphson technique the */
/*     solution to a set of four non-linear equations. The initial */
/*     primitives are in the array prim0, and the final solution is */
/*     returned in the same array. func is the function that computes the */
/*     four error terms for a given set of primitives. */
/*     niter is the maximum number of Newton-Raphson iterations */
/*     iflag=0 means that a good solution was found */
    /* Parameter adjustments */
    --prim0;

    /* Function Body */
    niter = 100;
/*      niter=20 */
    *iflag = 0;
    *jflag = 0;
/*     Do Newton-Raphson until err is smaller than tolerance tol */
    *iter = 0;
    i__1 = niter;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     Make sure u == prim0(1) is in safe territory */
	d__1 = accuracy_1.uminfac * conserved_2.rhou0;
	prim0[1] = mymax_(&prim0[1], &d__1);
	(*func)(&prim0[1], error0, errornorm, err4, &err3, iflag, jflag);
/*     iflag=9 means a serious error in the value of u. This can be */
/*     tolerated for a few steps (up to iter=itermax). After that, return */
/*     with iflag=9. */
	if (*iflag == 9 && *iter >= accuracy_1.itermax) {
	    return 0;
	}
/*     If all the four equations give fractional errors less than tol, */
/*     return */
	if (*err4 < accuracy_1.tol) {
	    return 0;
	}
	if (i__ > 1) {
	    sum = 0.;
	    for (j = 1; j <= 4; ++j) {
		d__1 = prim0[j] - primold[j - 1];
		sum += myabs_(&d__1);
	    }
	    if (sum == 0.) {
		return 0;
	    }
	}
	for (j = 1; j <= 4; ++j) {
	    primold[j - 1] = prim0[j];
	}
/*     Calculate the Jacobian numerically by shifting the primitives one */
/*     by one by a fraction eps of their current values and calculating */
/*     the errors. */
	for (j = 1; j <= 4; ++j) {
	    for (k = 1; k <= 4; ++k) {
		prim[k - 1] = prim0[k];
	    }
	    if (j == 1 || myabs_(&prim0[j]) > accuracy_1.dvmin) {
		dprim = prim0[j] * accuracy_1.eps;
	    } else {
		dprim = accuracy_1.dvmin * accuracy_1.eps;
	    }
	    prim[j - 1] = prim0[j] + dprim;
/*         prim(j)=prim0(j)-myabs(dprim) */
	    (*func)(prim, error, errornorm, err4, &err3, iflag, jflag);
	    for (k = 1; k <= 4; ++k) {
		ajac[k + (j << 2) - 5] = (error[k - 1] - error0[k - 1]) / 
			dprim;
	    }
	}
/*     Invert the Jacobian using subroutines from Numerical Recipes and */
/*     compute the shifts to the primitives */
	ludcmp_(ajac, &c__4, &c__4, indx, &d__, &retval);
	if (retval != 0.) {
	    *err4 = 128.f;
	    return 0;
	}
	lubksb_(ajac, &c__4, &c__4, indx, error0);
/*     Apply the Newton-Raphson shifts */
	for (j = 1; j <= 4; ++j) {
	    prim0[j] -= error0[j - 1];
	}
	++(*iter);
    }
    *iflag = 1;
    return 0;
} /* newton4_ */

/* Subroutine */ int newton3_(doublereal *prim0, integer *iter, integer *
	iflag, integer *jflag, U_fp func)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
     doublereal primsave, d__;
     integer i__, j, k;
     doublereal errornorm[3], err3, err4, ajac[9]	/* was [3][3] */;
     integer indx[3];
    extern /* Subroutine */ int func3_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     U_fp);
     doublereal prim3[3], prim30[3];
    extern doublereal myabs_(doublereal *);
     doublereal dprim;
     integer niter;
    extern doublereal mymax_(doublereal *, doublereal *);
     doublereal error3[3];
    extern /* Subroutine */ int lubksb_(doublereal *, integer *, integer *, 
	    integer *, doublereal *), ludcmp_(doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
     doublereal error30[3], retval;

/*     Calculates iteratively by the Newton-Raphson technique the */
/*     solution to a set of three non-linear equations. The initial */
/*     primitives are in the 4-array prim0, of which the first element is */
/*     not varied and only the other three are solved for. The final */
/*     solution is returned in the same 4-array with the first element */
/*     unchanged. func is the function that computes the four error terms */
/*     for a given set of primitives. Only the last three errors are */
/*     used. */
/*     niter is the maximum number of Newton-Raphson */
/*     iterations. Currently we do only one iteration of Newton3. */
/*     iflag=0 measn that a good solution was found */
    /* Parameter adjustments */
    --prim0;

    /* Function Body */
    niter = 1;
/*      niter=5 */
    *iflag = 0;
    *jflag = 0;
/*     Make sure u == primsave(1) is in safe territory */
    d__1 = accuracy_1.uminfac * conserved_2.rhou0;
    primsave = mymax_(&prim0[1], &d__1);
    for (i__ = 1; i__ <= 3; ++i__) {
	prim30[i__ - 1] = prim0[i__ + 1];
    }
/*     Do Newton-Raphson until err is less than tolerance tol */
    *iter = 0;
    for (i__ = 1; i__ <= 100; ++i__) {
	func3_(&primsave, prim30, error30, errornorm, &err4, &err3, iflag, 
		jflag, (U_fp)func);
	if (*iter >= niter) {
	    prim0[1] = primsave;
	    for (j = 1; j <= 3; ++j) {
		prim0[j + 1] = prim30[j - 1];
	    }
	    return 0;
	}
/*     If all the three equations give fractional errors less than tol, */
/*     return */
	if (err3 < accuracy_1.tol) {
	    prim0[1] = primsave;
	    for (j = 1; j <= 3; ++j) {
		prim0[j + 1] = prim30[j - 1];
	    }
	    return 0;
	}
/*     Calculate the Jacobian numerically by shifting the primitives one */
/*     by one by a fraction eps of their current values and calculating */
/*     the errors. */
	for (j = 1; j <= 3; ++j) {
	    for (k = 1; k <= 3; ++k) {
		prim3[k - 1] = prim30[k - 1];
	    }
	    if (myabs_(&prim30[j - 1]) > accuracy_1.dvmin) {
		dprim = prim30[j - 1] * accuracy_1.eps;
	    } else {
		dprim = accuracy_1.dvmin * accuracy_1.eps;
	    }
	    prim3[j - 1] = prim30[j - 1] + dprim;
/*         prim3(j)=prim30(j)-myabs(dprim) */
	    func3_(&primsave, prim3, error3, errornorm, &err4, &err3, iflag, 
		    jflag, (U_fp)func);
	    for (k = 1; k <= 3; ++k) {
		ajac[k + j * 3 - 4] = (error3[k - 1] - error30[k - 1]) / 
			dprim;
	    }
	}
/*     Invert the Jacobian using subroutines from Numerical Recipes and */
/*     compute the shifts to the primitives */
	ludcmp_(ajac, &c__3, &c__3, indx, &d__, &retval);
	if (retval != 0.) {
	    err3 = 1e4;
	    return 0;
	}
	lubksb_(ajac, &c__3, &c__3, indx, error30);
/*     Apply the Newton-Raphson shifts */
	for (j = 1; j <= 3; ++j) {
	    prim30[j - 1] -= error30[j - 1];
	}
	++(*iter);
    }
    *iflag = 1;
    prim0[1] = primsave;
    for (i__ = 1; i__ <= 3; ++i__) {
	prim0[i__ + 1] = prim30[i__ - 1];
    }
    return 0;
} /* newton3_ */

/* Subroutine */ int usolvemhd_(doublereal *prim, integer *iter, integer *
	iflag, integer *jflag, S_fp funcmhd, S_fp umhd)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
     integer i__;
     doublereal u0, errornorm[4], el, er, ul, ur, rho0, err0, err3, 
	    err4, emid, umid, ehat0;
    extern doublereal myabs_(doublereal *);
     doublereal error[4];

/*     Solves a 1D equation for u, using either the energy or entropy */
/*     equation without radiation source term */
/*     Call funcMHD and obtain basic parameters needed for solving the 1D */
/*     energy equation: rho, Ehat */
    /* Parameter adjustments */
    --prim;

    /* Function Body */
    u0 = prim[1];
    (*funcmhd)(&prim[1], error, errornorm, &err4, &err3, iflag, jflag);
    rho0 = funcmhdd_1.rho;
    ehat0 = funcmhdd_1.ehat;
/*     Calculate initial error */
    (*umhd)(&u0, &u0, &rho0, &ehat0, &funcmhdd_1.ffkap, &funcmhdd_1.arad, &
	    funcmhdd_1.dtau, &err0);
/*     Decide whether u is too low or too high and search accordingly to */
/*     bracket the solution for u */
    if (err0 < 0.) {
	ur = u0;
	er = err0;
	for (i__ = 1; i__ <= 50; ++i__) {
	    ul = ur;
	    el = er;
	    ur *= 2.;
	    (*umhd)(&ur, &u0, &rho0, &ehat0, &funcmhdd_1.ffkap, &
		    funcmhdd_1.arad, &funcmhdd_1.dtau, &er);
	    if (er >= 0.) {
		goto L10;
	    }
	}
	return 0;
    } else {
	ul = u0;
	el = err0;
	for (i__ = 1; i__ <= 50; ++i__) {
	    ur = ul;
	    er = el;
	    ul *= .5;
	    (*umhd)(&ul, &u0, &rho0, &ehat0, &funcmhdd_1.ffkap, &
		    funcmhdd_1.arad, &funcmhdd_1.dtau, &el);
	    if (er < 0.) {
		goto L10;
	    }
	}
	return 0;
    }
/*     Bracketing is done. Now solve for u by bisection */
L10:
    for (i__ = 1; i__ <= 50; ++i__) {
	umid = (ul + ur) * .5;
	(*umhd)(&umid, &u0, &rho0, &ehat0, &funcmhdd_1.ffkap, &
		funcmhdd_1.arad, &funcmhdd_1.dtau, &emid);
	if (emid >= 0.) {
	    ur = umid;
	    er = emid;
	} else {
	    ul = umid;
	    el = emid;
	}
	d__1 = ur - ul;
	if (myabs_(&d__1) < accuracy_1.epsbis * ul) {
	    prim[1] = (ul + ur) * .5;
	    return 0;
	}
    }
    prim[1] = (ul + ur) * .5;
    return 0;
} /* usolvemhd_ */

/* Subroutine */ int umhd1_(doublereal *u, doublereal *u0, doublereal *rho, 
	doublereal *ehat, doublereal *ffkap, doublereal *arad, doublereal *
	dtau, doublereal *err)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal tgas, entropy;

/*     Error function for 1D search in u for the energy equation without */
/*     radiation source term. Currently, this has the entropy equation. */
/*     Compute entropy and compute deviation from conserved entropy */
    tgas = conserved_2.gam1 * *u / *rho;
/*      B4pi=arad*Tgas**4 */
/*      ff=ffkap*rho*rho/Tgas**(3.5d0) */
/*      Gdt=ff*(Ehat-B4pi)*dt */
    d__1 = conserved_2.gam1 * *u;
    entropy = conserved_2.rhou0 * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(
	    rho, &conserved_2.en1));
    *err = entropy - conserved_2.s;
/*      err=entropy-s-Gdtau */
    return 0;
} /* umhd1_ */

/* Subroutine */ int umhd2_(doublereal *u, doublereal *u0, doublereal *rho, 
	doublereal *ehat, doublereal *ffkap, doublereal *arad, doublereal *
	dtau, doublereal *err)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal tgas, entropy;

/*     Error function for 1D search in u for the entropy equation without */
/*     radiation source term. */
/*     Compute entropy and compute deviation from conserved entropy */
    tgas = conserved_2.gam1 * *u / *rho;
/*      B4pi=arad*Tgas**4 */
/*      ff=ffkap*rho*rho/Tgas**(3.5d0) */
/*      Gdt=ff*(Ehat-B4pi)*dt */
    d__1 = conserved_2.gam1 * *u;
    entropy = conserved_2.rhou0 * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(
	    rho, &conserved_2.en1));
    *err = entropy - conserved_2.s;
/*      err=entropy-s-Gdtau */
    return 0;
} /* umhd2_ */

/* Subroutine */ int usolverad_(doublereal *prim, integer *iter, integer *
	iflag, integer *jflag, S_fp funcrad, U_fp uerr)
{
     integer i__;
     doublereal errornorm[4], el, er, ul, ur, rho0, err0, err3, err4, 
	    uacc, derr, ugas, ehat0, error[4];
    extern /* Subroutine */ int uuerr_(doublereal *, doublereal *, doublereal 
	    *);
    extern doublereal rtsafe_(S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*     Solves a 1D equation for u, using either the energy or entropy */
/*     equation including the radiation source term */
/*     Call funcrad and obtain basic parameters needed for solving the 1D */
/*     energy equation: rho, Ehat */
    /* Parameter adjustments */
    --prim;

    /* Function Body */
    funcradd_1.u0 = prim[1];
    (*funcrad)(&prim[1], error, errornorm, &err4, &err3, iflag, jflag);
    rho0 = funcradd_1.rho;
    ehat0 = funcradd_1.ehat;
/*     Calculate initial error */
    funcradd_1.iuerr = 0;
    uuerr_(&funcradd_1.u0, &err0, &derr);
/*      call uerr(u0,u0,rho0,Ehat0,ffkap,arad,dtau,err0) */
/*     Decide whether u is too low or too high and search accordingly to */
/*     bracket the solution for u */
    if (err0 < 0.) {
	ur = funcradd_1.u0;
	er = err0;
	for (i__ = 1; i__ <= 50; ++i__) {
	    ul = ur;
	    el = er;
	    ur *= 3.;
	    uuerr_(&ur, &er, &derr);
	    if (er >= 0.) {
		goto L10;
	    }
	}
	return 0;
    } else {
	ul = funcradd_1.u0;
	el = err0;
	for (i__ = 1; i__ <= 50; ++i__) {
	    ur = ul;
	    er = el;
	    ul *= .3;
	    uuerr_(&ul, &el, &derr);
	    if (el < 0.) {
		goto L10;
	    }
	}
	return 0;
    }
/*     Bracketing is done. Now solve for u using rtsafe */
L10:
    funcradd_1.iuerr = 1;
    uacc = accuracy_1.eps * min(ul,ur);
    ugas = rtsafe_((S_fp)uuerr_, &ul, &ur, &el, &er, &uacc);
    prim[1] = ugas;
/*      do i=1,50 */
/*         umid=0.5d0*(ul+ur) */
/*         call uerr(umid,u0,rho0,Ehat0,ffkap,arad,dtau,emid) */
/*         if (emid.ge.0.d0) then */
/*            ur=umid */
/*            er=emid */
/*         else */
/*            ul=umid */
/*            el=emid */
/*         endif */
/*         if (myabs(ur-ul).lt.epsbis*ul) then */
/*            prim(1)=0.5d0*(ul+ur) */
/*            return */
/*         endif */
/*      enddo */
/*      prim(1)=0.5d0*(ul+ur) */
    return 0;
} /* usolverad_ */

/* Subroutine */ int uerr1_(doublereal *u, doublereal *u0, doublereal *rho, 
	doublereal *ehat0, doublereal *ffkap, doublereal *arad, doublereal *
	dtau, doublereal *err)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal ff, b4pi, ehat, tgas, gdtau, entropy;

/*     Error function for 1D search in u for the energy equation */
/*     including the radiation source term. Currently, this has the */
/*     entropy equation. */
/*     Compute entropy and compute deviation from conserved entropy */
    tgas = conserved_2.gam1 * *u / *rho;
/* Computing 4th power */
    d__1 = tgas, d__1 *= d__1;
    b4pi = *arad * (d__1 * d__1);
    ff = *ffkap * *rho * *rho / pow_dd(&tgas, &c_b57);
    ehat = *u0 + *ehat0 - *u;
    gdtau = ff * (ehat - b4pi) * *dtau;
    d__1 = conserved_2.gam1 * *u;
    entropy = conserved_2.rhou0 * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(
	    rho, &conserved_2.en1));
    *err = entropy - conserved_2.s - gdtau;
    return 0;
} /* uerr1_ */

/* Subroutine */ int uerr2_(doublereal *u, doublereal *u0, doublereal *rho, 
	doublereal *ehat0, doublereal *ffkap, doublereal *arad, doublereal *
	dtau, doublereal *err)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal ff, b4pi, ehat, tgas, gdtau, entropy;

/*     Error function for 1D search in u for the entropy equation */
/*     including the radiation source term */
/*     Compute entropy and compute deviation from conserved entropy */
    tgas = conserved_2.gam1 * *u / *rho;
/* Computing 4th power */
    d__1 = tgas, d__1 *= d__1;
    b4pi = *arad * (d__1 * d__1);
    ff = *ffkap * *rho * *rho / pow_dd(&tgas, &c_b57);
    ehat = *u0 + *ehat0 - *u;
    gdtau = ff * (ehat - b4pi) * *dtau;
    d__1 = conserved_2.gam1 * *u;
    entropy = conserved_2.rhou0 * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(
	    rho, &conserved_2.en1));
    *err = entropy - conserved_2.s - gdtau;
    return 0;
} /* uerr2_ */

/* Subroutine */ int uuerr_(doublereal *u, doublereal *err, doublereal *derr)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
     doublereal c1, c2, c3, c4, ff, ehat, gdtau, entropy;

/*     The fluid frame error is: */
/*        rho*ln(p^n/rho^(n+1)) - s - (ffkap*rho^2/T^3.5)*dtau*(Ehat-B4pi) */
/*     We compute coefficients corresponding to these terms and then */
/*     compute the error and its derivative wrt u */
    c1 = conserved_2.gam1 / funcradd_2.rho;
    funcradd_2.tgas = c1 * *u;
/* Computing 4th power */
    d__1 = c1, d__1 *= d__1;
    c2 = funcradd_2.arad * (d__1 * d__1);
/* Computing 4th power */
    d__1 = *u, d__1 *= d__1;
    funcradd_2.b4pi = c2 * (d__1 * d__1);
    c3 = funcradd_2.ffkap * funcradd_2.rho * funcradd_2.rho / pow_dd(&c1, &
	    c_b57);
    ff = c3 / pow_dd(u, &c_b57);
    c4 = funcradd_2.u0 + funcradd_2.ehat0;
    ehat = c4 - *u;
/* Computing 4th power */
    d__1 = *u, d__1 *= d__1;
    gdtau = c3 * funcradd_2.dtau * pow_dd(u, &c_b63) * (c4 - *u - c2 * (d__1 *
	     d__1));
    d__1 = conserved_2.gam1 * *u;
    entropy = funcradd_2.rho * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(&
	    funcradd_2.rho, &conserved_2.en1));
    *err = entropy - conserved_2.s / funcradd_2.ucon[0] - gdtau;
    if (funcradd_2.iuerr == 1) {
	*derr = funcradd_2.rho * conserved_2.en / *u + c3 * 3.5 * 
		funcradd_2.dtau * c4 * pow_dd(u, &c_b64) - c3 * 2.5 * 
		funcradd_2.dtau * pow_dd(u, &c_b63) + c3 * .5 * 
		funcradd_2.dtau * c2 * pow_dd(u, &c_b66);
/* L10: */
    } else {
    }
    return 0;
} /* uuerr_ */

/* Subroutine */ int funcmhd1_(doublereal *prim, doublereal *error, 
	doublereal *errornorm, doublereal *err4, doublereal *err3, integer *
	iflag, integer *jflag)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer i__;
     doublereal u;
    extern /* Subroutine */ int calctmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), solveucon_(doublereal *
	    , doublereal *), calcbconbcov_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal myabs_(doublereal *), mydiv_(doublereal *, doublereal *)
	    , mymax_(doublereal *, doublereal *);
    extern integer myisnan_(doublereal *);
     doublereal entropy;

/*     This subroutine calculates errors for the MHD inversion problem */
/*     without radiation source term using the energy equation */
/*     Save u, compute rho, ucon, ucov, bcon, bcov, Tmunu */
    /* Parameter adjustments */
    --errornorm;
    --error;
    --prim;

    /* Function Body */
    u = prim[1];
    solveucon_(&prim[1], funcmhdd_1.ucon);
    contocov_(funcmhdd_1.ucon, funcmhdd_1.ucov);
    funcmhdd_1.rho = conserved_2.rhou0 / funcmhdd_1.ucon[0];
    if (funcmhdd_1.rho < 0.f || myisnan_(&funcmhdd_1.rho) == 1) {
	*jflag = 1;
    } else {
	*jflag = 0;
    }
    calcbconbcov_(conserved_2.bb, funcmhdd_1.ucon, funcmhdd_1.ucov, 
	    funcmhdd_1.bcon, funcmhdd_1.bcov, &funcmhdd_1.bsq);
    calctmunu_(&funcmhdd_1.rho, &u, &funcmhdd_1.bsq, &conserved_2.gam, 
	    funcmhdd_1.ucon, funcmhdd_1.ucov, funcmhdd_1.bcon, 
	    funcmhdd_1.bcov, funcmhdd_1.tmunu);
/*     Compute normalized err3 from the three momentum equations and the */
/*     normalized err4 from all four equations */
    *err3 = 0.;
    for (i__ = 2; i__ <= 4; ++i__) {
	error[i__] = funcmhdd_1.tmunu[(i__ << 2) - 4] - conserved_2.ttcov[i__ 
		- 1];
	d__2 = myabs_(&error[i__]) + myabs_(&funcmhdd_1.tmunu[(i__ << 2) - 4])
		 + myabs_(&conserved_2.ttcov[i__ - 1]);
	d__1 = mydiv_(&error[i__], &d__2);
	errornorm[i__] = myabs_(&d__1) + 1e-300;
	*err3 = mymax_(err3, &errornorm[i__]);
    }
/*     Use lab frame energy equation */
    error[1] = funcmhdd_1.tmunu[0] - conserved_2.ttcov[0];
    d__2 = myabs_(&error[1]) + myabs_(funcmhdd_1.tmunu) + myabs_(
	    conserved_2.ttcov);
    d__1 = mydiv_(&error[1], &d__2);
    errornorm[1] = myabs_(&d__1) + 1e-300;
    *err4 = mymax_(err3, &errornorm[1]);
/*     If u has an unreasonable value, set iflag=9 and reset u */
    if (u < 0.) {
	*iflag = 9;
    } else {
	d__1 = conserved_2.gam1 * u;
	entropy = log(pow_dd(&d__1, &conserved_2.en) / pow_dd(&funcmhdd_1.rho,
		 &conserved_2.en1));
	if (entropy < conserved_2.s / conserved_2.rhou0 - accuracy_1.dlogmax) 
		{
	    *iflag = 9;
	}
    }
    if (*iflag == 9) {
	d__1 = accuracy_1.uminfac * funcmhdd_1.rho;
	prim[1] = mymax_(&prim[1], &d__1);
    }
    return 0;
} /* funcmhd1_ */

/* Subroutine */ int funcmhd2_(doublereal *prim, doublereal *error, 
	doublereal *errornorm, doublereal *err4, doublereal *err3, integer *
	iflag, integer *jflag)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer i__;
     doublereal u;
    extern /* Subroutine */ int calctmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), solveucon_(doublereal *
	    , doublereal *);
     doublereal bsq, rho, bcon[4], bcov[4], ucon[4], ucov[4];
    extern /* Subroutine */ int calcbconbcov_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal myabs_(doublereal *), mydiv_(doublereal *, doublereal *)
	    , mymax_(doublereal *, doublereal *);
     doublereal tmunu[16]	/* was [4][4] */;
    extern integer myisnan_(doublereal *);
     doublereal entropy;

/*     This subroutine calculates errors for the MHD inversion problem */
/*     without radiation source term using the entropy equation */
/*     Save u, compute rho, ucon, ucov, bcon, bcov, Tmunu */
    /* Parameter adjustments */
    --errornorm;
    --error;
    --prim;

    /* Function Body */
    u = prim[1];
    solveucon_(&prim[1], ucon);
    contocov_(ucon, ucov);
    rho = conserved_2.rhou0 / ucon[0];
    if (rho < 0.f || myisnan_(&rho) == 1) {
	*jflag = 1;
    } else {
	*jflag = 0;
    }
    calcbconbcov_(conserved_2.bb, ucon, ucov, bcon, bcov, &bsq);
    calctmunu_(&rho, &u, &bsq, &conserved_2.gam, ucon, ucov, bcon, bcov, 
	    tmunu);
/*     Compute normalized err3 from the three momentum equations and the */
/*     normalized err4 from all four equations */
    *err3 = 0.;
    for (i__ = 2; i__ <= 4; ++i__) {
	error[i__] = tmunu[(i__ << 2) - 4] - conserved_2.ttcov[i__ - 1];
	d__2 = myabs_(&error[i__]) + myabs_(&tmunu[(i__ << 2) - 4]) + myabs_(&
		conserved_2.ttcov[i__ - 1]);
	d__1 = mydiv_(&error[i__], &d__2);
	errornorm[i__] = myabs_(&d__1) + 1e-300;
	*err3 = mymax_(err3, &errornorm[i__]);
    }
/*     Use entropy equation */
    d__1 = conserved_2.gam1 * u;
    entropy = ucon[0] * rho * log(pow_dd(&d__1, &conserved_2.en) / pow_dd(&
	    rho, &conserved_2.en1));
    error[1] = entropy - conserved_2.s;
    d__2 = myabs_(&error[1]) + myabs_(&entropy) + myabs_(&conserved_2.s);
    d__1 = mydiv_(&error[1], &d__2);
    errornorm[1] = myabs_(&d__1) + 1e-300;
    *err4 = mymax_(err3, &errornorm[1]);
/*     If u has an unreasonable value, set iflag=9 and reset u */
    if (prim[1] < 0.) {
	*iflag = 9;
	d__1 = accuracy_1.uminfac * rho;
	prim[1] = mymax_(&prim[1], &d__1);
    }
    return 0;
} /* funcmhd2_ */

/* Subroutine */ int func3_(doublereal *primsave, doublereal *prim3, 
	doublereal *error3, doublereal *errornorm3, doublereal *err4, 
	doublereal *err3, integer *iflag, integer *jflag, S_fp func)
{
     integer i__;
     doublereal errornorm4[4], prim4[4], error4[4];

/*     This function is called by Newton3. It takes a 3-array with */
/*     primitives prim3(3), transfers the value primsave to prim4(1) and */
/*     the three velocity components to prim4(2-4) and computes the */
/*     4-array of errors error4(4). Then transfers appropriate elements */
/*     to error3(3) and returns this along with overall error err3. */
/*     Transfer primsave and prim3(3) to prim4(4) */
    /* Parameter adjustments */
    --errornorm3;
    --error3;
    --prim3;

    /* Function Body */
    prim4[0] = *primsave;
    for (i__ = 1; i__ <= 3; ++i__) {
	prim4[i__] = prim3[i__];
    }
/*     Call appropriate function to calculate error4(4) */
    (*func)(prim4, error4, errornorm4, err4, err3, iflag, jflag);
/*     Transfer the momentum equation errors to error3(3) and return */
    for (i__ = 1; i__ <= 3; ++i__) {
	error3[i__] = error4[i__];
	errornorm3[i__] = errornorm4[i__];
    }
    return 0;
} /* func3_ */

/* Subroutine */ int funcrad1_(doublereal *prim, doublereal *error, 
	doublereal *errornorm, doublereal *err4, doublereal *err3, integer *
	iflag, integer *jflag)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);

    /* Local variables */
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer i__, j;
     doublereal u;
     integer radinvmod;
    extern /* Subroutine */ int calctmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), solveucon_(doublereal *
	    , doublereal *);
     doublereal ff, es, rt[4], tt[4], chi1, chi2;
    extern /* Subroutine */ int rmunuinvert_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , calcbconbcov_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
     doublereal alpha;
    extern doublereal myabs_(doublereal *), mydiv_(doublereal *, doublereal *)
	    , mymax_(doublereal *, doublereal *);
    extern integer myisnan_(doublereal *);
     doublereal entropy;

/*     This subroutine calculates errors for the radiation inversion */
/*     problem including the radiation source term using the energy */
/*     equation */
/*      ffkap=3.46764d-17 */
/*      eskap=5.90799d5 */
/*      arad=1.18316d17 */
    /* Parameter adjustments */
    --errornorm;
    --error;
    --prim;

    /* Function Body */
    funcradd_1.ffkap = conserved_2.okappa_ff_code11__;
    funcradd_1.eskap = conserved_2.okappa_es_code11__;
    funcradd_1.arad = conserved_2.arad_code__;
/*     Save u and solve for lab frame conserved quantities corresponding */
/*     to the given primitives: rho, ucon, ucov */
    u = prim[1];
/*     Calculate u^0 and rho */
    solveucon_(&prim[1], funcradd_1.ucon);
    contocov_(funcradd_1.ucon, funcradd_1.ucov);
    funcradd_1.rho = conserved_2.rhou0 / funcradd_1.ucon[0];
/*     Calculate b^mu, b_mu, bsq */
    calcbconbcov_(conserved_2.bb, funcradd_1.ucon, funcradd_1.ucov, 
	    funcradd_1.bcon, funcradd_1.bcov, &funcradd_1.bsq);
/*     Calculate the full gas stress-energy tensor T^mu_nu */
    calctmunu_(&funcradd_1.rho, &u, &funcradd_1.bsq, &conserved_2.gam, 
	    funcradd_1.ucon, funcradd_1.ucov, funcradd_1.bcon, 
	    funcradd_1.bcov, funcradd_1.tmunu);
/*     Evaluate the first row of the radiation stress-energy tensor: */
/*     R^0_mu */
    for (i__ = 1; i__ <= 4; ++i__) {
	tt[i__ - 1] = funcradd_1.tmunu[(i__ << 2) - 4];
	rt[i__ - 1] = conserved_2.ttcov[i__ - 1] + conserved_2.rtcov[i__ - 1] 
		- tt[i__ - 1];
    }
/*     Calculate the full R^mu_nu tensor */
    rmunuinvert_(rt, &funcradd_1.e, funcradd_1.urcon, funcradd_1.urcov, 
	    funcradd_1.ucon, funcradd_1.rmunu, &radinvmod);
/*     Calculate radiation energy density in the gas frame \hat{E}, and */
/*     the gas and radiation temperatures */
    funcradd_1.ehat = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    funcradd_1.ehat += funcradd_1.rmunu[i__ + (j << 2) - 5] * 
		    funcradd_1.ucov[i__ - 1] * funcradd_1.ucon[j - 1];
	}
    }
    d__1 = funcradd_1.ehat / funcradd_1.arad;
    funcradd_1.trad = pow_dd(&d__1, &c_b73);
    funcradd_1.tgas = conserved_2.gam1 * u / funcradd_1.rho;
    if (funcradd_1.trad < 0.f || myisnan_(&funcradd_1.trad) == 1) {
	*jflag = 1;
    } else {
	*jflag = 0;
    }
/*     Calculate quantities needed to compute the radiation source term */
    ff = funcradd_1.ffkap * funcradd_1.rho * funcradd_1.rho / pow_dd(&
	    funcradd_1.tgas, &c_b57);
    es = funcradd_1.eskap * funcradd_1.rho;
/* Computing 4th power */
    d__1 = funcradd_1.tgas, d__1 *= d__1;
    funcradd_1.b4pi = funcradd_1.arad * (d__1 * d__1);
/*     The following side calculation is to estimate quantities that are */
/*     relevant for deciding whether the radiation source term can be */
/*     handled explicitly. Currently, all calculations are done */
/*     implicitly. */
    alpha = 1. / sqrt(-metric_1.gn[0]);
    funcradd_1.gamma = funcradd_1.ucon[0] * alpha;
/*     Decide how to set dtau: either dt/gamma or dt/ucon(1) */
/*      dtau=dt/gamma */
    funcradd_1.dtau = conserved_2.dt / funcradd_1.ucon[0];
    chi1 = ff * funcradd_1.dtau * (funcradd_1.b4pi * 4. / u + 1.);
    chi2 = (ff + es) * funcradd_1.dtau * (funcradd_1.ehat / (funcradd_1.rho + 
	    conserved_2.gam * u) + 1.);
/*     Compute the radiation source term G^mu and G_mu directly in the */
/*     lab frame (formula taken from Jon's draft of the paper, with some */
/*     corrections) */
    for (i__ = 1; i__ <= 4; ++i__) {
	funcradd_1.gcon[i__ - 1] = -(ff * funcradd_1.b4pi + es * 
		funcradd_1.ehat) * funcradd_1.ucon[i__ - 1];
	for (j = 1; j <= 4; ++j) {
	    funcradd_1.gcon[i__ - 1] -= (ff + es) * funcradd_1.rmunu[i__ + (j 
		    << 2) - 5] * funcradd_1.ucon[j - 1];
	}
    }
    contocov_(funcradd_1.gcon, funcradd_1.gcov);
/*     Compute normalized err3 from the three momentum equations and the */
/*     normalized err4 from all four equations */
    *err3 = 0.;
    for (i__ = 2; i__ <= 4; ++i__) {
	error[i__] = tt[i__ - 1] - conserved_2.ttcov[i__ - 1] - 
		funcradd_1.gcov[i__ - 1] * conserved_2.dt;
	d__3 = funcradd_1.gcov[i__ - 1] * conserved_2.dt;
	d__2 = myabs_(&d__3) + myabs_(&tt[i__ - 1]) + myabs_(&
		conserved_2.ttcov[i__ - 1]);
	d__1 = mydiv_(&error[i__], &d__2);
	errornorm[i__] = myabs_(&d__1) + 1e-300;
	*err3 = mymax_(err3, &errornorm[i__]);
    }
/*     Calculate error(1) from lab frame energy equation */
    error[1] = tt[0] - conserved_2.ttcov[0] - funcradd_1.gcov[0] * 
	    conserved_2.dt;
    d__3 = funcradd_1.gcov[0] * conserved_2.dt;
    d__4 = (conserved_2.gam * u + funcradd_1.bsq) * funcradd_1.ucon[0] * 
	    funcradd_1.ucov[0];
    d__2 = myabs_(&d__3) + myabs_(tt) + myabs_(conserved_2.ttcov) + myabs_(&
	    d__4);
    d__1 = mydiv_(&error[1], &d__2);
    errornorm[1] = myabs_(&d__1);
    *err4 = mymax_(err3, &errornorm[1]);
/*     If u has an unreasonable value, set iflag=9 and reset u */
    if (u < 0.) {
	*iflag = 9;
    } else {
	d__1 = conserved_2.gam1 * u;
	entropy = log(pow_dd(&d__1, &conserved_2.en) / pow_dd(&funcradd_1.rho,
		 &conserved_2.en1));
	if (entropy < conserved_2.s / conserved_2.rhou0 - accuracy_1.dlogmax) 
		{
	    *iflag = 9;
	}
    }
    if (*iflag == 9) {
	d__1 = accuracy_1.uminfac * funcradd_1.rho;
	prim[1] = mymax_(&prim[1], &d__1);
    }
    funcradd_1.u0 = prim[1];
    return 0;
} /* funcrad1_ */

/* Subroutine */ int funcrad2_(doublereal *prim, doublereal *error, 
	doublereal *errornorm, doublereal *err4, doublereal *err3, integer *
	iflag, integer *jflag)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);

    /* Local variables */
     doublereal ghatdtau;
    extern /* Subroutine */ int contocov_(doublereal *, doublereal *);
     integer i__, j;
     doublereal u;
     integer radinvmod;
    extern /* Subroutine */ int calctmunu_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), solveucon_(doublereal *
	    , doublereal *);
     doublereal ff, es, gt, rt[4], tt[4], chi1, chi2;
    extern /* Subroutine */ int rmunuinvert_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , calcbconbcov_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
     doublereal alpha;
    extern doublereal myabs_(doublereal *), mydiv_(doublereal *, doublereal *)
	    , mymax_(doublereal *, doublereal *);
    extern integer myisnan_(doublereal *);
     doublereal entropy;

/*     This subroutine calculates errors for the radiation inversion */
/*     problem including the radiation source term using the entropy */
/*     equation */
/*      ffkap=3.46764d-17 */
/*      eskap=5.90799d5 */
/*      arad=1.18316d17 */
    /* Parameter adjustments */
    --errornorm;
    --error;
    --prim;

    /* Function Body */
    funcradd_1.ffkap = conserved_2.okappa_ff_code11__;
    funcradd_1.eskap = conserved_2.okappa_es_code11__;
    funcradd_1.arad = conserved_2.arad_code__;
/*     Save u and solve for lab frame conserved quantities corresponding */
/*     to the given primitives: rho, ucon, ucov */
    u = prim[1];
    if (prim[1] < 0.) {
	*iflag = 9;
    }
/*     Calculate u^0 and rho */
    solveucon_(&prim[1], funcradd_1.ucon);
    contocov_(funcradd_1.ucon, funcradd_1.ucov);
    funcradd_1.rho = conserved_2.rhou0 / funcradd_1.ucon[0];
/*     Calculate b^mu, b_mu, bsq */
    calcbconbcov_(conserved_2.bb, funcradd_1.ucon, funcradd_1.ucov, 
	    funcradd_1.bcon, funcradd_1.bcov, &funcradd_1.bsq);
/*     Calculate the full gas stress-energy tensor T^mu_nu */
    calctmunu_(&funcradd_1.rho, &u, &funcradd_1.bsq, &conserved_2.gam, 
	    funcradd_1.ucon, funcradd_1.ucov, funcradd_1.bcon, 
	    funcradd_1.bcov, funcradd_1.tmunu);
/*     Evaluate the first row of the radiation stress-energy tensor: */
/*     R^0_mu */
    for (i__ = 1; i__ <= 4; ++i__) {
	tt[i__ - 1] = funcradd_1.tmunu[(i__ << 2) - 4];
	rt[i__ - 1] = conserved_2.ttcov[i__ - 1] + conserved_2.rtcov[i__ - 1] 
		- tt[i__ - 1];
    }
/*     Calculate the full R^mu_nu tensor */
    rmunuinvert_(rt, &funcradd_1.e, funcradd_1.urcon, funcradd_1.urcov, 
	    funcradd_1.ucon, funcradd_1.rmunu, &radinvmod);
/*     Calculate radiation energy density in the gas frame \hat{E}, and */
/*     the gas and radiation temperatures */
    funcradd_1.ehat = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    funcradd_1.ehat += funcradd_1.rmunu[i__ + (j << 2) - 5] * 
		    funcradd_1.ucov[i__ - 1] * funcradd_1.ucon[j - 1];
	}
    }
    d__1 = funcradd_1.ehat / funcradd_1.arad;
    funcradd_1.trad = pow_dd(&d__1, &c_b73);
    funcradd_1.tgas = conserved_2.gam1 * u / funcradd_1.rho;
/*     Calculate quantities needed to compute the radiation source term */
    ff = funcradd_1.ffkap * funcradd_1.rho * funcradd_1.rho / pow_dd(&
	    funcradd_1.tgas, &c_b57);
    es = funcradd_1.eskap * funcradd_1.rho;
/* Computing 4th power */
    d__1 = funcradd_1.tgas, d__1 *= d__1;
    funcradd_1.b4pi = funcradd_1.arad * (d__1 * d__1);
    if (funcradd_1.trad < 0.f || myisnan_(&funcradd_1.trad) == 1) {
	*jflag = 1;
    } else {
	*jflag = 0;
    }
/*     The following side calculation is to estimate quantities that are */
/*     relevant for deciding whether the radiation source term can be */
/*     explicitly. Currently, all calculations are done implicitly. */
    alpha = 1. / sqrt(-metric_1.gn[0]);
    funcradd_1.gamma = funcradd_1.ucon[0] * alpha;
/*     Decide how to set dtau: either dt/gamma or dt/ucon(1) */
/*      dtau=dt/gamma */
    funcradd_1.dtau = conserved_2.dt / funcradd_1.ucon[0];
    ghatdtau = ff * (funcradd_1.ehat - funcradd_1.b4pi) * funcradd_1.dtau;
    chi1 = ff * funcradd_1.dtau * (funcradd_1.b4pi * 4. / u + 1.);
    chi2 = (ff + es) * funcradd_1.dtau * (funcradd_1.ehat / (funcradd_1.rho + 
	    conserved_2.gam * u) + 1.);
/*     Compute the radiation source term G^mu and G_mu directly in the */
/*     lab frame (formula taken from Jon's draft of the paper, with some */
/*     corrections) */
    for (i__ = 1; i__ <= 4; ++i__) {
	funcradd_1.gcon[i__ - 1] = -(ff * funcradd_1.b4pi + es * 
		funcradd_1.ehat) * funcradd_1.ucon[i__ - 1];
	for (j = 1; j <= 4; ++j) {
	    funcradd_1.gcon[i__ - 1] -= (ff + es) * funcradd_1.rmunu[i__ + (j 
		    << 2) - 5] * funcradd_1.ucon[j - 1];
	}
    }
    contocov_(funcradd_1.gcon, funcradd_1.gcov);
/*     Compute normalized err3 from the three momentum equations and the */
/*     normalized err4 from all four equations */
    *err3 = 0.;
    for (i__ = 2; i__ <= 4; ++i__) {
	error[i__] = tt[i__ - 1] - conserved_2.ttcov[i__ - 1] - 
		funcradd_1.gcov[i__ - 1] * conserved_2.dt;
	d__3 = funcradd_1.gcov[i__ - 1] * conserved_2.dt;
	d__2 = myabs_(&d__3) + myabs_(&tt[i__ - 1]) + myabs_(&
		conserved_2.ttcov[i__ - 1]);
	d__1 = mydiv_(&error[i__], &d__2);
	errornorm[i__] = myabs_(&d__1) + 1e-300;
	*err3 = mymax_(err3, &errornorm[i__]);
    }
/*     Calculate error(1) corresponding to the lab frame entropy equation */
    d__1 = conserved_2.gam1 * prim[1];
    entropy = funcradd_1.ucon[0] * funcradd_1.rho * log(pow_dd(&d__1, &
	    conserved_2.en) / pow_dd(&funcradd_1.rho, &conserved_2.en1));
    gt = funcradd_1.gcov[0] * conserved_2.dt;
    error[1] = entropy - conserved_2.s - gt;
    d__2 = myabs_(&entropy) + myabs_(&conserved_2.s) + myabs_(&gt);
    d__1 = mydiv_(&error[1], &d__2);
    errornorm[1] = myabs_(&d__1) + 1e-300;
/*      entropy=rho*log((Gam1*prim(1))**en/rho**en1) */
/*      Gt=Ghatdtau */
/*      error(1)=entropy-(s/ucon(1))-Gt */
/*      errornorm(1)=(1D-300)+myabs(mydiv(error(1),(myabs(entropy)+ */
/*     &     myabs(s/ucon(1))+myabs(Gt)))) */
/*     Alternatively, use the fluid frame entropy equation */
/*      entropy=ucon(1)*rho*log((Gam1*prim(1))**en/rho**en1) */
/*      Gdtau=ff*(Ehat-B4pi)*dtau */
/*      Gt=Gdtau */
/*      error(1)=entropy-s-Gdtau */
    *err4 = mymax_(err3, &errornorm[1]);
    funcradd_1.u0 = prim[1];
    return 0;
} /* funcrad2_ */

/* Subroutine */ int lubksb_(doublereal *a, integer *n, integer *np, integer *
	indx, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
     integer i__, j, ii, ll;
     doublereal sum;

/*     Matrix inversion subroutine 1 from Numerical Recipes */
    /* Parameter adjustments */
    --b;
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ll = indx[i__];
	sum = b[ll];
	b[ll] = b[i__];
	if (ii != 0) {
	    i__2 = i__ - 1;
	    for (j = ii; j <= i__2; ++j) {
		sum -= a[i__ + j * a_dim1] * b[j];
/* L11: */
	    }
	} else if (sum != 0.) {
	    ii = i__;
	}
	b[i__] = sum;
/* L12: */
    }
    for (i__ = *n; i__ >= 1; --i__) {
	sum = b[i__];
	i__1 = *n;
	for (j = i__ + 1; j <= i__1; ++j) {
	    sum -= a[i__ + j * a_dim1] * b[j];
/* L13: */
	}
	b[i__] = sum / a[i__ + i__ * a_dim1];
/* L14: */
    }
    return 0;
} /* lubksb_ */

/* Subroutine */ int ludcmp_(doublereal *a, integer *n, integer *np, integer *
	indx, doublereal *d__, doublereal *retval)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
     integer i__, j, k;
     doublereal vv[500], dum, sum;
     integer imax;
     doublereal aamax;
    extern doublereal myabs_(doublereal *);

/*     Matrix inversion subroutine 2 from Numerical Recipes */
    /* Parameter adjustments */
    --indx;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *d__ = 1.;
    imax = -100;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aamax = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (myabs_(&a[i__ + j * a_dim1]) > aamax) {
		aamax = myabs_(&a[i__ + j * a_dim1]);
	    }
/* L11: */
	}
	if (aamax == 0.) {
	    *retval = 1.;
	    return 0;
	}
	vv[i__ - 1] = 1. / aamax;
/* L12: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = a[i__ + j * a_dim1];
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
/* L13: */
	    }
	    a[i__ + j * a_dim1] = sum;
/* L14: */
	}
	aamax = 0.;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum = a[i__ + j * a_dim1];
	    i__3 = j - 1;
	    for (k = 1; k <= i__3; ++k) {
		sum -= a[i__ + k * a_dim1] * a[k + j * a_dim1];
/* L15: */
	    }
	    a[i__ + j * a_dim1] = sum;
	    dum = vv[i__ - 1] * myabs_(&sum);
	    if (dum >= aamax) {
		imax = i__;
		aamax = dum;
	    }
/* L16: */
	}
	if (imax > *n || imax < 1 || imax == -100) {
	    *retval = 2.;
	    return 0;
	}
	if (j != imax) {
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		dum = a[imax + k * a_dim1];
		a[imax + k * a_dim1] = a[j + k * a_dim1];
		a[j + k * a_dim1] = dum;
/* L17: */
	    }
	    *d__ = -(*d__);
	    vv[imax - 1] = vv[j - 1];
	}
	indx[j] = imax;
	if (a[j + j * a_dim1] == 0.) {
	    a[j + j * a_dim1] = 1e-80;
	}
	if (j != *n) {
	    dum = 1. / a[j + j * a_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= dum;
/* L18: */
	    }
	}
/* L19: */
    }
    *retval = 0.;
    return 0;
} /* ludcmp_ */

doublereal rtsafe_(S_fp funcd, doublereal *x1, doublereal *x2, doublereal *fl,
	 doublereal *fr, doublereal *xacc)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
     doublereal f;
     integer j;
     doublereal df, dx, xh, xl, temp, dxold;
    extern doublereal myabs_(doublereal *);

/*      call funcd(x1,fl,df) */
/*      call funcd(x2,fh,df) */
/*      if(fl.eq.0.d0)then */
/*        rtsafe=x1 */
/*        return */
/*      else if(fh.eq.0.d0)then */
/*        rtsafe=x2 */
/*        return */
    if (*fl < 0.) {
	xl = *x1;
	xh = *x2;
    } else {
	xh = *x1;
	xl = *x2;
    }
    ret_val = (*x1 + *x2) * .5;
    (*funcd)(&ret_val, &f, &df);
    if (f < 0.) {
	xl = ret_val;
    } else {
	xh = ret_val;
    }
    d__1 = xh - xl;
    dxold = myabs_(&d__1);
    dx = dxold;
    for (j = 1; j <= 100; ++j) {
	d__1 = f * 2.;
	d__2 = dxold * df;
	if (((ret_val - xh) * df - f) * ((ret_val - xl) * df - f) >= 0. || 
		myabs_(&d__1) > myabs_(&d__2)) {
	    dxold = dx;
	    dx = (xh - xl) * .5;
	    ret_val = xl + dx;
	    if (xl == ret_val) {
		return ret_val;
	    }
	} else {
	    dxold = dx;
	    dx = f / df;
	    temp = ret_val;
	    ret_val -= dx;
	    if (temp == ret_val) {
		return ret_val;
	    }
	}
	if (myabs_(&dx) < *xacc) {
	    return ret_val;
	}
	(*funcd)(&ret_val, &f, &df);
	if (f < 0.) {
	    xl = ret_val;
	} else {
	    xh = ret_val;
	}
/* L11: */
    }
    return ret_val;
} /* rtsafe_ */

