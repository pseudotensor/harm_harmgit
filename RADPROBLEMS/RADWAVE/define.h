#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MINKOWSKI

#define RK3STEPPING
#define INT_ORDER 1
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 2.

#define ALLSTEPSOUTPUT 0

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

//#define EDDINGTON_APR
#define KAPPAES 0.
#define KAPPA 0.
#define NX 100

//#define RADIATION

#undef SIGMA_RAD

#define NWAVE 5
//#define RADOUTPUTINZAMO


#if (NWAVE==5) //sound wave with radiation set up according to Jiang+12
//#define FLUXDISSIPATIONOFF

#define NUMERO 41

#if (NUMERO==41)
#define PP 100.
#define CC 1.e2
#undef KAPPA
#define KAPPA 10.
#define RHOFAC 0.01
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define DVRE (4.06372e-6*RHOFAC)
#define DVIM (6.90937e-6*RHOFAC)
#define DURE (9.88671e-4*RHOFAC)
#define DUIM (6.97077e-6*RHOFAC)
#define DERE (-4.52724e-5*RHOFAC)
#define DEIM (2.78566e-5*RHOFAC)
#define DFRE (-5.83678e-6*RHOFAC)
#define DFIM (-9.48194e-6*RHOFAC)
#define OMRE 0.0255331
#define OMIM 0.0434128
#define DTOUT1 1.e-0
#endif

#if (NUMERO==11)
#define PP 0.01
#define CC 1.e2
#undef KAPPA
#define KAPPA 0.01
#define RHOFAC 0.01
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define DVRE (9.99998e-6*RHOFAC)
#define DVIM (8.48878e-9*RHOFAC)
#define DURE (1.66666e-3*RHOFAC)
#define DUIM (2.82938e-6*RHOFAC)
#define DERE (1.95853e-8*RHOFAC)
#define DEIM (1.91123e-7*RHOFAC)
#define DFRE (-1.33508e-5*RHOFAC)
#define DFIM (4.23463e-6*RHOFAC)
#define OMRE 6.28317e-2
#define OMIM 5.33366e-5
#define DTOUT1 1.e-0
#endif


#if (NUMERO==1)
#define PP 0.01
#define CC 1.e4
#undef KAPPA
#define KAPPA 0.01
#define DRRE 1.e-3
#define DRIM 0.
#define DVRE 9.7548e-8
#define DVIM 7.92788e-9
//#define DPRE 1.61075e-3
//#define DPIM 2.07402e-4
#define DURE 1.57546e-3
#define DUIM 2.57783e-4
#define DERE 1.6874e-8//1.79137e-8
#define DEIM 9.48966e-9//8.56498e-9
#define DFRE -1.77115e-6//-1.32035e-6
#define DFIM 3.65291e-6//3.88814e-6
#define OMRE 6.12912e-4//7.99077
#define OMIM 4.98123e-5//0.512336
#define DTOUT1 1.e-2
#endif

#define RHO (ldouble)1.
#define KK 2.*Pi
#define UINT (ldouble)((1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)) //to get proper sound speed
#define TEMP (ldouble)(calc_PEQ_Tfromurho(UINT,RHO)) //temperature from rho and uint
#define SIGMA_RAD (ldouble)((1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)) //to get the proper radiation to gas pressure ratio, PP=4 sig T^4 / P
#define ERAD (ldouble)(calc_LTE_EfromT(TEMP)) //to get thermal equilibrium, E=4 sig T^4

#define RADIATION
#endif


#if (NWAVE==1) //density wave advected with the gas
//#define FLUXDISSIPATIONOFF //switches of the dissipative term in LAXF 
#define PP 0.1
#define CC 1.e6
#define VX 1.e-3
#define DTOUT1 (.05/VX)
#define RHO 1.
#define AAA 1.e-5
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#endif

#if (NWAVE==2) //hydro sound wave
//#define FLUXDISSIPATIONOFF
#define PP 0.01
#define CC 1.e6
#define DTOUT1 (.05*CC)
#define VX 0.
#define RHO 1.
#define AAA 1.e-5
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#endif

#if (NWAVE==3) //radiative density wave advected with the gas
#define FLUXDISSIPATIONOFF
#define PP 10.
#define CC 1.e6
#define VX 1.e-2
#define DTOUT1 (.0005/VX)
#define RHO 1.
#define AAA 1.e-5
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#define ERAD calc_LTE_EfromT(TEMP)
#define RADIATION
#undef KAPPAES
#define KAPPAES 10.
#endif


#if (NWAVE==4) //sound wave with radiation, set up without the phase shifts etc.
#define FLUXDISSIPATIONOFF
#define PP 1.
#define CC 1.e2
#define DTOUT1 (.005*CC)
#define VX 0.
#define RHO 1.
#define AAA 1.e-1
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#define ERAD calc_LTE_EfromT(TEMP)
#define RADIATION
#undef KAPPA
#define KAPPA 100.
#define ERADFACTOR .5
#define GASFACTOR .5
#endif




#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40
