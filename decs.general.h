// BELOW USED BY makeopenmpsharedlist.sh to generate OPENMPSHAREDLIST, so should come last among all things in this file.
extern int BEGINOPENMPSHAREDLIST;



#include "mpidecs.h"


#include "kazfulleos.decsglobalprivate.h" // put here so OpenMP private globals are defined before global.nondepnmemonics.h sets up thread private pragma's



/////////////////////////////////////////////////////////////////
//
//
// GLOBAL VARIABLES that are not for every point in space
//
// OPENMPNOTE: Note that any global variable *written to* inside parallel region must be made private and dealt with or neglected in value.
//             Otherwise race condition occurs and if results depend upon that variable, then result will be undefined and wrong in general.
//
//
/////////////////////////////////////////////////////////////////


// OPENMPMARK: below aa,S,n were static and calledranc was global, but ranc() may be called by multiple threads, so required to be global and use omp critical for that region.  Otherwise would have to make threadprivate and copyin() everytime think used, which is nasty.
extern int rancaa[NRANC];
extern int rancS[NRANC];
extern int rancvaln;


extern FTYPE Xmetricnew[NDIM],Xmetricold[NDIM]; // used to store time of latest and oldest metric

extern SFTYPE *lumvsr,*lumvsr_tot;

extern SFTYPE *dissvsr[NUMDISSVERSIONS],*dissvsr_tot[NUMDISSVERSIONS];

extern FTYPE *rcent,*rcent_tot;

extern SFTYPE *dVvsr,*dVvsr_tot,*vrsqvsr,*vrsqvsr_tot,*dMvsr,*dMvsr_tot,*dTrrvsr, *dTrrvsr_tot,*Mvsr_tot,*Mvsrface1_tot,*MOrvsr_tot,*phivsr_tot,*dJvsr,*dJvsr_tot,*Jvsr_tot,*Jvsrface1_tot;


/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern FTYPE gam,gamideal;

/* numerical parameters */
extern int defcoord;
extern FTYPE Rin, R0, Rout, hslope, Zin, Zout;
extern FTYPE Rin_array[NDIM], Rout_array[NDIM];  //atch -- arrays for a more general way of handling the grid dimensions
extern FTYPE Risco,Rhor;
extern FTYPE cour;
extern FTYPE dV, dVF, dx[NDIM], startx[NDIM], endx[NDIM], Diffx[NDIM];
extern SFTYPE dt,t,tf,tstepparti,tsteppartf;
extern FTYPE TDYNORYEglobal,Hglobal;
extern FTYPE rcurr, hcurr;
extern FTYPE drsing;

//int istart, istop, jstart, jstop;
#if(SIMULBCCALC!=-1) 
extern int isc,iec,jsc,jec;
extern int isf1,ief1,jsf1,jef1,ksf1,kef1;
extern int isf2,ief2,jsf2,jef2,ksf2,kef2;
extern int isf3,ief3,jsf3,jef3,ksf3,kef3;
extern int ise,iee,jse,jee;
extern int isf1ct,ief1ct,jsf1ct,jef1ct;// GODMARK: other stage type requires more
extern int isf2ct,ief2ct,jsf2ct,jef2ct;
extern int isf3ct,ief3ct,jsf3ct,jef3ct;
extern int isdq,iedq,jsdq,jedq;
extern int ispdq,iepdq,jspdq,jepdq;
#endif

extern FTYPE mydminarg1, mydminarg2;
extern long nstep;
extern int specialstep;

extern int steppart,numstepparts;

extern int gocont; // used to continue running(runtype, directinput, timereenter)
extern int runtype;

/* output parameters */
extern FILE* log_file;
extern FILE* fail_file;
extern FILE* logfull_file;
extern FILE* logdt_file;
extern FILE* logdtfull_file;
extern FILE* logstep_file;
extern FILE* logperf_file;
extern SFTYPE DTstep,DTstepdot,DTperf,DTgocheck,DTtimecheck,DTperfdump;
extern int reallaststep,onemorestep;

// file version stuff:
extern int PVER,GRIDVER,DVER,FLVER,NPVER,AVG1DVER,AVG2DVER,ENERVER,MODEVER,LOSSVER,SPVER,TSVER,LOGDTVER,STEPVER,PERFVER,ADVER,PDVER,CALCVER,FLINEVER;

extern int PTYPE,GRIDTYPE,DTYPE,FLTYPE,NPTYPE,AVG2DTYPE,AVG1DTYPE,ENERTYPE,LOSSTYPE,SPTYPE,TSTYPE,LOGDTTYPE,STEPTYPE,PERFTYPE,ADTYPE,PDTYPE,CALCTYPE,FLINETYPE,MODETYPE,EXPANDTYPE,NPCOMPUTETYPE;


extern SFTYPE DTdumpgen[NUMDTDS];
extern long dumpcntgen[NUMDTDS];
//SFTYPE DTd;
//SFTYPE DTavg;
//SFTYPE DTener;
//SFTYPE DTi;
//SFTYPE DTdebug;
extern long DTr;
//long dump_cnt;
//long avg_cnt;
//long debug_cnt;
//long image_cnt;
extern long rdump_cnt;
//long fieldline_cnt; // assumed to keep track with images (as in diag.c), so no need to include in restart()

extern int nstroke; // OPENMPMARK: Bad in some inversion codes


//for holding absolute values of indices of regions -- for restarting
extern FTYPE t_transition_in,t_transition_out;
extern int global_enerregiondef[NUMENERREGIONS][NUMUPDOWN][NDIM];




/* global flags */
extern int failed;
extern int lim[NDIM],fluxmethod,FLUXB,UTOPRIMVERSION,TIMEORDER,DOENOFLUX,avgscheme[NDIM];
extern int dofluxreconevolvepointfield,emffixedstencil,extrazones4emf,splitmaem,unewisavg;
extern int do_transverse_flux_integration[NPR],do_conserved_integration[NPR],do_source_integration[NPR];
extern int useghostplusactive;
extern FTYPE defcon;

/* diagnostics */
// don't track this separately in other regions except global region
extern SFTYPE frdot[N1][NPR];
extern SFTYPE pdottermsjet2[COMPDIM*2][NUMFLUXTERMS][NPR];
extern CTYPE failfloorcountlocal[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet
extern CTYPE failfloorcountlocal_tot[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet

// general stuff for ener.out file for regions to completely track, including terms within flux
extern int dothisenerreg[NUMENERREGIONS];
extern int dofluxreg[NUMENERREGIONS][COMPDIM*2];
extern int enerposreg[NUMENERREGIONS][COMPDIM*2];
// these quantities contain diagnostics
// all these require writing to restart file
// other _tot quantities appear in dump_ener.c that don't need to be written to restart file since easily computed from existing data.
extern SFTYPE fladdreg[NUMENERREGIONS][NPR];
extern SFTYPE fladdreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE fladdtermsreg[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE fladdtermsreg_tot[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE Ureg_init[NUMENERREGIONS][NPR];
extern SFTYPE Ureg_init_tot[NUMENERREGIONS][NPR];
extern SFTYPE pcumreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pcumreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdotreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdottermsreg[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
extern SFTYPE sourceaddreg[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddtermsreg[NUMENERREGIONS][NUMSOURCES][NPR];
extern SFTYPE sourceaddtermsreg_tot[NUMENERREGIONS][NUMSOURCES][NPR];
extern SFTYPE dissreg[NUMENERREGIONS][NUMDISSVERSIONS];
extern SFTYPE dissreg_tot[NUMENERREGIONS][NUMDISSVERSIONS];
//SFTYPE horizonflux[NPR];
//SFTYPE horizoncum[NPR];
//SFTYPE horizonflux_tot[NPR];
//SFTYPE horizoncum_tot[NPR];

// below quantities are not kept in restart file since easily recomputed
// kept global so can always access current value throughout code
extern SFTYPE pdotreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdottermsreg_tot[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];

// used for each region, related to global quantities
// _tot quantities here are global since used in restart.
// dangerous as global quantity since could overlap if nesting use of enerregion stuff! GODMARK
extern int *doflux;
extern int *enerpos;
extern SFTYPE *fladd;
extern SFTYPE *fladd_tot;
extern SFTYPE (*fladdterms)[NPR];
extern SFTYPE (*fladdterms_tot)[NPR];
extern SFTYPE *U_init;
extern SFTYPE *U_init_tot;
extern SFTYPE (*pcum)[NPR];
extern SFTYPE (*pcum_tot)[NPR];
extern SFTYPE (*pdot)[NPR];
extern SFTYPE (*pdotterms)[NUMFLUXTERMS][NPR];
extern SFTYPE *sourceadd;
extern SFTYPE *sourceadd_tot;
extern SFTYPE (*sourceaddterms)[NPR];
extern SFTYPE (*sourceaddterms_tot)[NPR];
extern SFTYPE *diss;
extern SFTYPE *diss_tot;

// kept global
extern SFTYPE (*pdot_tot)[NPR];
extern SFTYPE (*pdotterms_tot)[NUMFLUXTERMS][NPR];


// end changes after ...


/* Jon's addition */
extern int horizoni,horizoncpupos1;
extern long realnstep;
extern int partialstep;
extern int mpicombine;
extern int mpicombinetype;
extern int truempicombinetype;
extern int halftimep;
extern int whichrestart;
extern int appendold;
extern int whocalleducon; // OPENMPNOTE: Ensure those are set as threadprivate [noted only called outside parallel regions]
// global flags
extern long restartsteps[2];
extern int binaryoutput,sortedoutput;
extern int CHECKCONT,DOTSTEPDIAG,DOLOGSTEP,DOLOGPERF;
extern int NDTCCHECK,NZCCHECK,NDTDOTCCHECK,NGOCHECK,NTIMECHECK,NDTPERFDUMPCHECK;
extern SFTYPE PERFWALLTIME,ZCPSESTIMATE;

extern long steptofaildump,steptofailmap;
extern int ifail,jfail,kfail; // OPENMPNOTE: Ensure those are set as private [noted only for diagnostics outside parallel regions]
extern int dofailmap,dofaildump,restartonfail;
// IC
extern FTYPE h_over_r;
// BC
extern FTYPE h_over_r_jet;
extern int BCtype[COMPDIM*2];
extern int rescaletype;
extern int cooling;
extern int DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG, DOIMAGEDIAG,DOAREAMAPDIAG;
extern int GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC,DOENODEBUGEVERYSUBSTEP,DODIAGEVERYSUBSTEP;
extern int INVERTFROMAVERAGEIFFAILED,LIMIT_AC_PRIM_FRAC_CHANGE,LIMIT_AC_FRAC_CHANGE; //atch
extern int PARAMODWENO;
extern FTYPE MAX_AC_FRAC_CHANGE,MAX_AC_PRIM_FRAC_CHANGE;
extern FTYPE RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT;
extern FTYPE prMAX[NPR];
extern FTYPE prfloorcoef[NPR];
extern FTYPE BSQORHOLIMIT,BSQOULIMIT,UORHOLIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL;
extern FTYPE SAFE;
extern int debugfail;
extern FTYPE uttdiscr;  // OPENMPNOTE: Ensure those are set as threadprivate [noted only for WHICHVEL==VEL3]
extern int jonchecks;
extern int dnumcolumns[NUMDUMPTYPES];
extern struct blink * blinkptr0[NUMDUMPTYPES];
extern struct blink * cpulinkptr0[NUMDUMPTYPES];
extern int DOCOLSPLIT[NUMDUMPTYPES];
extern int docolsplit; // global var for now CHANGINGMARK
extern int nextcol;
extern int doevolvemetricsubsteps, gravityskipstep;
extern FTYPE gravitydtglobal, sourcedtglobal, wavedtglobal;
extern int waveglobaldti[NDIM],waveglobaldtj[NDIM],waveglobaldtk[NDIM];
extern int didstorepositiondata,didstoremetricdata;

/* physical consts */
extern FTYPE msun,lsun,rsun,G,H,C,qe,Na,malpha,mn,me,kb,arad,sigmasb,sigmamat,mevocsq,ergPmev,mp,Q,R,Re,hpl,hbar,K,K2;
extern SFTYPE a,MBH,QBH;
extern FTYPE Mfactor,Jfactor,rhofactor;
extern SFTYPE dabh,dE,dJ,dEold,dJold;
extern FTYPE mb,mbcsq,mbwithrhounit,amu,a0,MBH0,QBH0,Mdot,Mdotc,Mcgs,Ccode;
extern FTYPE Lunit,Tunit,Vunit,rhounit,rhomassunit,Munit,mdotunit,energyunit,edotunit,Pressureunit,Tempunit,Bunit,massunitPmsun;
extern int rho0unittype;
extern FTYPE ledd,leddcode;

extern int NUMBUFFERS;

extern int advancepassnumber;

// OPENMPNOTE: Ensure all pl's are set as private as required [only npr2interp and npr2notinterp changed in parallel regions] and removed use of "plglobal" type globals
// for SPLITNPR (NOW generally used):
// for choosing range of PLOOP type arrays
extern int nprstart,nprend; // normally 0 and NPR-1
extern int nprlist[MAXNPR]; // maximum is NPR elements

// for choosing range of PLOOPINTERP type arrays
extern int npr2interpstart,npr2interpend; // normally 0 and NPR2INTERP-1
extern int npr2interplist[MAXNPR]; // maximum is NPR2INTERP elements

// for choosing range of PLOOPNOTINTERP type arrays
extern int npr2notinterpstart,npr2notinterpend; // normally 0 and -1
extern int npr2notinterplist[MAXNPR]; // maximum is NPR2INTERP elements

// for choosing range of PBOUNDLOOP and PLOOPMPI type arrays
extern int nprboundstart,nprboundend; // normally 0 and NPRBOUND-1
extern int nprboundlist[MAXNPR]; // maximum is NPRBOUND elements

// for choosing range of PFLUXBOUNDLOOP and PLOOPMPI type arrays
extern int nprfluxboundstart,nprfluxboundend; // normally 0 and NPRFLUXBOUND-1
extern int nprfluxboundlist[MAXNPR]; // maximum is NPRFLUXBOUND elements

// for choosing range of PDUMPLOOP
extern int nprdumpstart,nprdumpend; // normally 0 and NPRDUMP-1
extern int nprdumplist[MAXNPR]; // maximum is NPRDUMP elements

// for choosing range of PINVERTLOOP
extern int nprinvertstart,nprinvertend; // normally 0 and NPRINVERT-1
extern int nprinvertlist[MAXNPR]; // maximum is NPRINVERT elements


extern int fluxloop[NDIM][NUMFLUXLOOPNUMBERS];
extern int emffluxloop[NDIM][NUMFLUXLOOPNUMBERS];
extern int Uconsloop[NUMFLUXLOOPNUMBERS];
extern int emfUconsloop[NUMFLUXLOOPNUMBERS];
extern int Uconsevolveloop[NUMFLUXLOOPNUMBERS];
extern int a_interporder[NUMINTERPS];
extern int *interporder;


// ENO DEBUG GLOBAL VARS
//int dirglobal,locglobal,iglobal,jglobal,kglobal,iterglobal,interporfluxglobal;

// SOME GEOMETRIC VARIABLES
extern int special3dspc;

extern int numbercpu[ 3+1 ];

// Ramesh stuff
extern FTYPE nu,ss,ucrit,Ttpow,jetalpha;


// EOS related functions
extern FTYPE (*ptr_pressure_rho0_u)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_u_from_entropy)(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy);
extern FTYPE (*ptr_u_rho0_p)(FTYPE *EOSextra, FTYPE rho0, FTYPE p);
extern FTYPE (*ptr_dpdu_rho0_u)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_dpdrho0_rho0_u)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_cs2_compute)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_entropy)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_dSdrho)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_compute_dSdu)(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE (*ptr_pressure_wmrho0) (FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE (*ptr_compute_idwmrho0dp) (FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE (*ptr_compute_idrho0dp) (FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE (*ptr_compute_qdot) (FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern int (*ptr_compute_sources_EOS) (FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR]);
extern void (*ptr_compute_allextras) (int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras);
extern int (*ptr_get_extrasprocessed) (int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE (*ptr_compute_temp) (FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern void (*ptr_compute_EOS_parms) (FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS],  FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void (*ptr_store_EOS_parms) (int numparms, FTYPE *EOSextra, FTYPE *parlist);
extern void (*ptr_get_EOS_parms) (int*numparms, FTYPE *EOSextra, FTYPE *parlist);

extern FTYPE SQRTMINNUMREPRESENT;

extern FTYPE NUMEPSILONPOW23;

// some grid section or loop things
extern int AVOIDADVANCESHIFTX1DN,AVOIDADVANCESHIFTX1UP,AVOIDADVANCESHIFTX2DN,AVOIDADVANCESHIFTX2UP,AVOIDADVANCESHIFTX3DN,AVOIDADVANCESHIFTX3UP,GLOBALBCMOVEDWITHACTIVESECTION;


#include "decs.user.h"



extern int crapdebug;



// BELOW USED BY makeopenmpsharedlist.sh to generate OPENMPSHAREDLIST, so should come last among all things in this file.
extern int ENDOPENMPSHAREDLIST;


#pragma omp threadprivate(OPENMPGLOBALPRIVATELIST)
