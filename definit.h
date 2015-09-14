
/*! \file definit.h
    \brief Default values for code choices/switches appearing in user's init.h
    
*/


/// Set array storage order as defined in global.storage.h
/// 0 corresponds to original HARM
/// 5 optimal for extended radial grid, less so in \theta, and even less so in \phi
#define ORDERSTORAGE 5

/// whether to use new metric storage method that directly creates geom-> structure that is directly used rather than doing indirect memory copies and calculations
#define NEWMETRICSTORAGE 1

/// whether to check symmetry along diagonal (for test==151 in Sasha's tests)
#define ASYMDIAGCHECK 0

/// whether to use accurate (bit more expensive) diagnostics.  Only necessary for timeorder>2
#define ACCURATEDIAG 1
#define FLUXDIAGPARTS 0 // whether to compute parts of fluxes

/// whether to accurately update diagnostics related to source() terms.
/// Not really necessary since those source terms don't lead to perfect conservation anyways
/// Calling diag_source_comp() can be somewhat expensive, so just do on final steps for now (i.e. below set to 0)
// 0 : not accurate
// 1 : cumulative source accurate
// 2 : each component accurate
#define ACCURATESOURCEDIAG 1
#define DIAGSOURCECOMPSTEP 4 // how many fullsteps to wait to compute diag_source_comp() that is somewhat expensive but just diagnostics and not incredibly needed compared to diag_source_all()


/// whether to use sin/cos exactly symmetrized around pi/2
/// if one doesn't do this, then \theta=\pi pole can lead to \detg=small number instead of 0.0
#define ACCURATESINCOS 1

/// whether to try entropy inversion if hot fails
#define HOT2ENTROPY 1
/// whether to try cold inversion if hot fails (only used if cold makes sense)
#define HOT2COLD 1
/// whether to try FFDE inversion if hot fails (only used if FFDE makes sense)
#define HOT2FFDE 1

/// whether to try cold inversion if entropy fails
#define ENTROPY2COLD 1
/// whether to try HOT inversion if entropy fails
#define ENTROPY2HOT 1
/// whether to try FFDE inversion if entropy fails (only used if FFDE makes sense)
#define ENTROPY2FFDE 1




#define BADENERGYMAC(ugenergy,ugentropy) ((ugenergy) < 0.5*(ugentropy))
#define BADENERGY2MAC(ugenergy,ugentropy) ((ugenergy) < (ugentropy))

/// which method to store divergence condition
#define DIVMETHODPREFLUX 0
#define DIVMETHODPOSTFLUX 1

#define DIVERGENCEMETHOD 1 // choose, but normally 1 for best results.


/// *NUMBER* OF DIMENSIONS FOR COMPUTATION
/// Can choose 3 and still do 1D optimized problems
/// Choosing 2 or even 1 reduces some memory load for those things that are accessed by using 1,2,3 instead of arbitrarily accessed
/// not sure if 1 works.
#define COMPDIM 3

/// whether to have get_geometry() point to global gcov,gcon,gcovpert -- faster and no limitations
/// apparently faster now to choose 0
///#define GETGEOMUSEPOINTER 1
/// not used when NEWMETRICSTORAGE==1
#define GETGEOMUSEPOINTER 0

/// whether to split NPR PLOOPs
#define SPLITNPR 0

/// whether to store analytical values -- required for some boundary condition methods
#define ANALYTICMEMORY 1

/// whether want to use weno memory (required for c2e or a2c or c2a for weno/eno)
#define WENOMEMORY 0

/// whether need higher order memory stuff
#define HIGHERORDERMEM 1

/// whether to split off MA and EM parts of flux for higher order
/// for now always do if doing higher order
/// GODMARK: If SPLITMAEM==0 and pointfield method, then only EMF terms have fixed stencil
///#define SPLITMAEMMEM (HIGHERORDERMEM)
#define SPLITMAEMMEM (0) // CHANGINGMARK

///#define DOMONOINTERP NOMONOINTERP
///#define DOMONOINTERP JMONOINTERP
#define DOMONOINTERP SMONOINTERP

/// whether to split isotropic pressure off of flux for separate interpolation
/// this splitting only used ultimately if splitmaem==1 (even if computed separately)
#define SPLITPRESSURETERMINFLUXMA 0 // CHANGINGMARK
/// below not setup yet:
#define SPLITPRESSURETERMINFLUXEM 0 // CHANGINGMARK


/// whether to do a2c during c2e
/// see interpline.c
#define MERGEDC2EA2CMETHODMA 0
#define MERGEDC2EA2CMETHODEM 0
#define MERGEDC2EA2CMETHOD (MERGEDC2EA2CMETHODMA || MERGEDC2EA2CMETHODEM)


#define HIGHERORDERROUGH 0
#define HIGHERORDERSMOOTH 1
#define HIGHERORDERVERYSMOOTH 2

///////////////////////////////////////////////////////////////////
///
/// CHOOSE SMOOTHNESS EXPECTED OF SOLUTION for HIGHERORDER TERMS:
///
///////////////////////////////////////////////////////////////////
#define HIGHERORDERTYPE HIGHERORDERROUGH // original method
///#define HIGHERORDERTYPE HIGHERORDERSMOOTH



#if(HIGHERORDERTYPE==HIGHERORDERROUGH)// non-smooth version:
/// for pointfield method doesn't keep divb=0 properly
/// but original WHAM method for a2c



#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C ENERGY_CONTROLS_ALL_WEIGHTS // c2a and a2c for conserved quants
#define DO_SPLITA2CSMOOTH DO_SPLITA2C
#define DO_SPLITA2C4FLUX DO_SPLITA2C // more generally
#define DO_SPLITA2C4MAFLUX MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS // couples weights together
#define DO_SPLITA2C4SMOOTHFLUX DO_SPLITA2C
#define DO_SPLITSOURCE NOSPLITA2C // c2a for ENOSOURCETERM

#define EMFLUXFIELDTERMTYPE ENOFLUX
#define FIELDTERMTYPE ENOFLUX
#define FIELDFVTYPE ENOCONSERVED



#elif(HIGHERORDERTYPE==HIGHERORDERSMOOTH)  // smooth version


/// should be used with pointfield method for general problems
/// SMOOTH-named variables used for magnetic field operations

#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C ENERGY_CONTROLS_ALL_WEIGHTS // c2a and a2c for conserved quants
#define DO_SPLITA2CSMOOTH CONSTANT_ALL_WEIGHTS
///#define DO_SPLITA2C4FLUX CONSTANT_ALL_WEIGHTS // needed for TESTNUMBER==49 if doing FLUXRECON -- GODMARK: perhaps good for smooth problems
///#define DO_SPLITA2C4FLUX NOSPLITA2C // so pressure changes don't control other terms (for example pressure can be constant but field or velocity can vary) -- causes problems for test=106 since E-M not coupled enough when correction occurs
///#define DO_SPLITA2C4FLUX ENERGY_IS_ALL_WEIGHTS
///#define DO_SPLITA2C4FLUX MINIMIZE_ALL_WEIGHTS
///#define DO_SPLITA2C4FLUX ENERGY_CONTROLS_ALL_WEIGHTS
/// CHANGINGMARK
#define DO_SPLITA2C4FLUX MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS
#define DO_SPLITA2C4MAFLUX DO_SPLITA2C4FLUX
///#define DO_SPLITA2C4FLUX MINIMIZE_ALL_WEIGHTS
#define DO_SPLITA2C4SMOOTHFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITSOURCE NOSPLITA2C // c2a for ENOSOURCETERM

#define EMFLUXFIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDFVTYPE ENOSMOOTHCONSERVED



#elif(HIGHERORDERTYPE==HIGHERORDERVERYSMOOTH) // supersmooth version
/// for problems with no very strong discontinuities (often good even for strong discontinuities)



#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2CSMOOTH CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4FLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4MAFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4SMOOTHFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITSOURCE NOSPLITA2C

#define EMFLUXFIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDFVTYPE ENOSMOOTHCONSERVED

#endif





////////////////
/// for MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS
/// used in reconstruct.c and interpline.mono.c
#define VELTERMSMINIMIZE(pl) (pl==RHO || pl==UU || pl==UU+dir)
#define ORTHOVEL1TERMSMINIMIZE(pl) (pl==UU+odir1)
#define ORTHOVEL2TERMSMINIMIZE(pl) (pl==UU+odir2)
#define PRESSUREMINIMIZE(pl) (pl==FLUXSPLITPMA(dir)) // pressure term when splitmaem==1 and doing MA term
#define EMFTERMS(pl) (pl>=B1 && pl<=B3) // only used if splitmaem==0

/// ALLOTHERS corresponds to all those pl's not specially handled by the minimization routine.  EMFTERMS only handled specially if emffixedstencil==1
/// below is always true unless hit certain terms that are true
/// 
#define ALLOTHERSMINIMIZE(pl) (! (VELTERMSMINIMIZE(pl) || ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl) || PRESSUREMINIMIZE(pl) || ((emffixedstencil==1)&&EMFTERMS(pl))  ) )
///#define OTHERADVECTMINIMIZE(pl) ( YFLPL(pl) && DOYFL!=DONOYFL || pl==YL && DOYL!=DONOYL || pl==YNU && DOYNU!=DONOYNU || pl==ENTROPY && DOENTROPY!=DONOENTROPY)



/// whether to use staggered fields (which requires memory for it)
#define FIELDSTAGMEM 1

/// whether to store get_state() data into array for later use since expensive to compute get_geometry() and ucon_calc() type things
#define STOREFLUXSTATE 1

/// whether to compute and store Shock (and other indicators) just once rather than for every call to slope_lim().  Since there can be many (e.g. 33 in 3D for staggered field method) that would be overly expensive for a simple diffusive check or correction that it itself is not higher-order!
#define STORESHOCKINDICATOR 1

/// whether to define Toth memory stuff
///#define FIELDTOTHMEM (!FIELDSTAGMEM) // CHANGINGMARK
#define FIELDTOTHMEM 1 // CHANGINGMARK


/// whether to track and dump A_i the vector potential
#define TRACKVPOT 0

/// whether to redefine fields from A_i during evolution
/// currently only done after each full timestep to avoid being expensive
/// only machine error different from evolution of field itself
#define EVOLVEWITHVPOT TRACKVPOT // choice

/// more often makes sense to modify A_i and not EMF_i, since want A_i to be smooth so Bstag^i is computed properly (and Bcent^i based upon Bstag^i)
#define MODIFYEMFORVPOT MODIFYVPOT

/// whether to specify gdet at end when setting EMF or to have internal to variables before averaging.
/// point is that gdet at end is probably better, esp. at coordinate singularities.
#define CORNGDETVERSION 1
///#define CORNGDETVERSION 0


/// size of global arrays
/// ALSO CONTROLS WHETHER MEMORY IS ALLOCATED FOR SPECIAL 1D/2D/3D directional arrays
/// e.g. if N3==1, then second dimension is neglected in all calculations

#define N1 64  // number of zones in 1-direction
#define N2 64        // number of zones in 2-direction
#define N3      1               // number of zones in 3-direction

///#if(DOENOFLUX==ENOFINITEVOLUME)
#define MAXBND 5
///#else
///#define MAXBND 3
///#endif


/// choose metric
#define MCOORD KSCOORDS


/// defines how one forms the EM stress-energy tensor
/// GENMAXWELL is perfectly fine for relativistic problems that are not near machine error for highly relativistic flows that usually is only for tests
///#define MAXWELL GENMAXWELL
/// prim version avoids catastrophic cancellation in non-rel limit
/// OPTMARK: MAXWELL GENMAXWELL faster if don't care about cancellation issues (i.e. not in ultra rel or non-rel regimes)
#define MAXWELL PRIMMAXWELL


/// whether to avoid catastrophic cancellations for non-rel and ultra-rel speed and gravity regimes with u_t+1
#define UD0PLUS1FIX 1


#define DO_VORTICITY_IMAGE 0


#define PRODUCTION 0
/// 0: full images, dumps, etc., few more key failure stderr/fail_file messages
/// 1: only log density image since too many images (takes alot of time), no utoprim failure messages -- assume debug.out and debug???? will have info if needed
/// 2: #1 but also avoid error_check() and avoid per-MPI-proc log and fail files (good for running on >1000 cores)
/// 3: #2 but also avoid dumps not needed except by python (remove all SM-related dumps)


/// 0: normal computational zones outputted on diagnostics
/// else, FULLOUTPUT # of extra boundary zones on each side (if allowed by dimensionality)
/// If # of requested boundary zones is larger than real #, then real # used

/// ONLY CAN BE USED WITH USEMPI==0
#define FULLOUTPUT 0


#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "ki-rh42.slac.stanford.edu"
#define REMOTEUSER "jon"
#define EMAILADDRESS "jmckinne@stanford.edu"
#define EMAILMESSAGE "Done with GRMHD run DEFAULT"



/// whether doing super long doubles
#define SUPERLONGDOUBLE 0


#define PERFTEST 0
/// 0: don't perform performance test
/// 1: do

#define DOAVG 0
#define DOAVGDUMP (DOAVG)
/// 0: don't do time average dumps, so don't allocate memory
/// 1: do

/// GODMARK: only setup for full Pi theta grid, not Pi/2 or any other theta grid.
/// e.g. will fail for defcoord=3 Pi/2 grid.  Causes code to crash
#define DOJETDIAG 1
/// 0: don't do jet diagnostics (ener, not flener that is always done)
/// 1: do

#define DOAVG2 0 // only make 1 if DOAVG 1 above
/// 0: don't split AVG file
/// 1: do


#define DOGDUMP 1
/// 0: don't output grid dump
/// 1: do

#define DORDUMP 1
/// 0: don't output restart dumps
/// 1: do

#define DOMAINDUMP 1
#define DOMAINDUMPDIAG (DOMAINDUMP)
/// 0: don't output main dump files
/// 1: do

#define DODEBUGDUMP 1
/// 0: don't output debug dump file or ener file(ener is based on dump counts)
/// 1: do

#define DODEBUG 1
/// 0: don't output debug dump file or ener file(ener is based on dump counts)
/// 1: do

#define DOFLOORDUMP 1
#define DOFLOORDIAG (DOFLOORDUMP)
/// 0: don't output file with diag_fixup() activated changes in conserved quantities
/// 1: do

#define DODISSMEASUREDUMP 0
#define DODISSMEASUREDIAG (DODISSMEASUREDUMP)
/// 0: don't do dissmeasure dump files
/// 1: do

#define DOEOSDUMP 1
#define DOEOSDIAG (DOEOSDUMP)
/// 0: don't do eosdump files
/// 1: do

#define DORADDUMP 1
#define DORADDIAG (DORADDUMP)
/// 0: don't do raddump files
/// 1: do


/// whether to dump vector potential
#define DOVPOTDUMP (1 && TRACKVPOT)

#define DOIMAGEDUMP (1)

#define DOAREAMAPDUMP (1)


#define DOENODEBUG 0
/// whether to do ENO debug




#define DOGRIDSECTIONING 0
/// 0: don't do grid sectioning
/// 1: do (only part of the grid in the active section is evolved)



#define WENO_USE_PRIM_REDUCTION 1
/// 0: do not perform the c2a limiting on the fluxes
/// 1: perform the limiting of c2a flux reconstruction based on the limiting of the a2c primitive correction:
///    if the a2c reconstruction of the conserved quantities leads to very different primitives, then no a2c reconstruction is done; in this case
///    no reconstruction is done on the fluxes either.

#define WENO_EXTRA_A2C_MINIMIZATION 0 // CHANGINGMARK

#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 1
/// 0: do not perform additional reduction
/// 1: perform additional reduction:
///    for NDIM > 1, subsequent 1d a2c reconstructions reduce to lower order (=no reconstruction) if any of the previous reconstructions reduced

#define LIMIT_FLUXC2A_PRIM_CHANGE 0
/// 0: do not limit primtive correction due to c2a done on fluxes
/// 1: limit c2a correction
 
#define DO_WENO_DEBUG 0
/// whether to debug WENO

#define DOSUPERDEBUG 0
/// 0: don't output super debug output
/// 1: do

/// whether to allow metric rotations
/// 0: don't
/// 1: do
#define ALLOWMETRICROT 0
/// TODO:
/// * coordinate transformations still valid?
/// * SET RESET THETAROT for IC so IC uses normal metric?


/// whether metric mixes r-\phi or \theta-\phi.  Alows for optimizations since not often mixing with \phi
#if(ALLOWMETRICROT)
#define DOMIXTHETAPHI 1 // for g_{\theta\phi} // no choice
#else
#define DOMIXTHETAPHI 0 // choice
#endif
/// careful not to use DOMIXTHETAPHI in .h files since init.h might change it.




/// whether to evolve metric value of M and a and Q
#define DOEVOLVEMETRIC 0

/// 0 = evolve only infrequently after many steps
/// 1 = evolve on every substep in full RK form
/// 2 = use accdt and ndt1/2/3 to determine whether one should evolve metric on substeps or longsteps
#define EVOLVEMETRICSUBSTEP 0

/// whether to limit timestep using source term at t=0 and full flux term for future times
#define LIMITDTWITHSOURCETERM 0

/// whether to limit source update by some constraint (gives dUriemann to source())
#define LIMITSOURCES 0


/// whether to use gravitydt to limit dt
#define USEGRAVITYDTINDTLIMIT 0


/// 0: no restriction beyond normal ghost+active bounds
/// 1: don't limit via accdt or gravitydt inside horizon
/// 2: in addition to #1, don't limit due to velocity inside horizon (can't be used if boundary condition drives solution and would limit dt more than active region)
#define RESTRICTDTSETTINGINSIDEHORIZON 1


#define DODISS 0

#define DODISSMEASURE 0
/// 0: don't do calculation
/// 1: do

/// see diag_source()
#define DOLUMVSR 0

/// see diag_source()
#define DODISSVSR 0

/// below only applies to dissipation inversion
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY
///#define WHICHENTROPYEVOLVE EVOLVEFULLENTROPY


/// see metric.c
#define DOSELFGRAVVSR 0

#define DOFIELDLINE 1
/// 0: don't output energy@infinity and field line stuff
/// 1: do

/// whether to use Roe-averaged state in determining diffusive term in flux
#define ROEAVERAGEDWAVESPEED 0
/// whether to use Athena's version of Roe averaging and estimating wave speeds
#define ATHENAROE 0

/// whether to store wave speeds over whole grid before computing flux
/// useful to avoid extra calculations if doing "local" LAXF or "global" LAXF.
/// default HARM was using VERY local LAXF (only wavespeeds from primitives interpolated to the edge).
#define STOREWAVESPEEDS 0

/// whether to compute per cell $dt$ by storing $dt$ per dimension and then computing minimum per cell rather than minimum per dimension.  Can give up to a factor of 3X improvement in speed (just changes effective Courant factor).
#define PERCELLDT 1

/// whether to use stored wave speeds for flux calculation (allows one to store wave speeds for interp.c but use true VERYLOCALVCHAR that is vchar's estimated from boundary as in standard HARM -- rather than maximum from center zones as done by STORED version of VERYLOCALVCHAR)
/// silly choice would be 0 if VCHARTYPE=GLOBALVCHAR since interp.c doesn't use the stored wave speeds if this is the case.  So shouldn't store in this case since no point.
#define USESTOREDSPEEDSFORFLUX (STOREWAVESPEEDS) // choice really independent of STOREWAVESPEEDS, but generall normally want to couple them

#define VCHARTYPE VERYLOCALVCHAR

/// TRUEFAST=0 : use overall larger phase velocity
/// TRUEFAST=1 : use true phase velocity (has problems with catastrophic cancellation)
#define TRUEFAST 0

/// whether to force dt to use velocity based upon speed of light
/// Doesn't change how fluxes are computed!
#define FORCESOLVEL 0
/// below changes how fluxes are computed
#define FORCESOLVELFLUX 0


/// whether to check on inversion and report problem
/// Also, use diag_fixup() to account for changes in U when inversion success (tracks changes due to Newton iteration error and machine error)
#define CHECKONINVERSION 1
#define CHECKONINVERSIONRAD 0 // check on radiation inversion.  Default 0 because generally limiters on |v| mean not the same, so only use for debugging.

#define PRECISEINVERSION 1
/// whether we use ultimately precise inversion or "workable" inversion precision (see utoprim's for how defined)

#define WHICHVEL VELREL4
/// which velocity to compute in (init can be anything (see init.c's for transformations)


/// whether to subtract rest-mass from energy equation for divT=0 equation of motion
/// 0: use MHD stress tensor with explicit rest-mass included
/// 1: as 0, but subtract out rest-mass from conservation and flux terms for evolution
/// 2: use MHD stress tensor withOUT rest-mass
/// this changes mhd_calc() in phys.c and assumes rest of code uses mhd stress-energy tensor without restmass also!!
/// DIAGNOSTICS also without rest-mass for UU terms.
#define REMOVERESTMASSFROMUU 1

#define RELTYPE RELEOM

#define EOMTYPE EOMGRMHD
///#define EOMTYPE EOMCOLDGRMHD

#define EOMRADTYPE EOMRADNONE // EOMRADM1CLOSURE
#define EVOLVENRAD 0 // whether to evolve number density of radiation

#define WHICHRADSOURCEMETHOD SOURCEMETHODNONE
///#define WHICHRADSOURCEMETHOD SOURCEMETHODIMPLICITEXPLICITCHECK



/// whether to try other methods for the inversion if primary choices fails
/// created because utoprim_2d_final fails for large b^2/rho when other methods (even utoprim_2d) do not fail.
/// seems that in some cases causes code to stall while 1D method takes forever on many iterations
#define UTOPRIMTRYAGAIN 0


/// which EOS (EoS) (equation of state) to use
#define WHICHEOS IDEALGAS
///#define WHICHEOS MIGNONE
///#define WHICHEOS GRBPWF99
///#define WHICHEOS KAZFULL



/// whether to activate entropy variable
#define DOENTROPY DONOENTROPY // normal total energy equation

/// whether to call fixup() after initialization
#define FIXUPAFTERINIT 1

/// whether to call fixup() after restart
#define FIXUPAFTERRESTART 1

/// checks if solution is reasonble and fails a point if not
/// only checks if b^2/\rho>>BSQORHOMAX, etc.
/// can cause erratic issues near BH, so default is turned off
/// KORALNOTEMARK: Not setup for radiation, but ok since not generally used.
#define CHECKSOLUTION 0



/// factor by which zone-zone values can be different for \gamma and internal energy, used of CHECKSOLUTION==1
#define GAMMAPERCDIFFMAX 2.0
#define UPERCDIFFMAX 10.0


/// whether to interpolate extra quantity (e.g. v^2) and use to renormalize velocities after they are interpolated
#define DOEXTRAINTERP 0
/// must also set RESCALEINTERP=1

/// whether to evolve Y_fl (see fixup.c)
#define DOYFL 0
/// whether to evolve Y_l (see kazeosfull.c)
#define DOYL 0
/// whether to evolve Y_\nu (see kazeosfull.c)
#define DOYNU DOYL


/// whether to control the limiter
#define LIMADJUST LIMITERFIXED
///#define LIMADJUST LIMITERBSQORHOANDU

/// whether to limit all variables or just hydro variables
#define HYDROLIMADJUSTONLY 1

/// determine how flux is computed per zone
#define FLUXADJUST FLUXFIXED
///#define FLUXADJUST FLUXBSQORHOANDU

/// whether to change flux calculation for all variables or just hydro variables
#define HYDROFLUXADJUSTONLY 1

/// whether to allow negative internal energy in substep
/// UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED should be set in global.h
#define STEPOVERNEGU NEGDENSITY_FIXONFULLSTEP
/// seems to make things worse (more failures, but also worse failures)

#define STEPOVERNEGRHO NEGDENSITY_FIXONFULLSTEP

#define STEPOVERNEGURHO NEGDENSITY_FIXONFULLSTEP




///#define UTOPRIMADJUST UTOPRIMSTATIC
#define UTOPRIMADJUST UTOPRIMAVG

#define UTOPRIMFIXMPICONSISTENT 1

/// 0: return non-updated quantities
/// 1: if u<0 or rho<0, then return updated quantities, else return non-updated quantities
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED

/// whether to smooth singularity of black hole over grid scale
#define SMOOTHSING 1

#define COORDSINGFIXCYL 0 //whether perform the same fix for CYLMINKMETRIC's z axis

#define COORDSINGFIX 1
/// whether to move polar axis to a bit larger theta
/// theta value where singularity is displaced to


/// SINGSMALL can't be smaller than DXDELTA in dxdxp (i.e. currently ~1E-5 for DOUBLE)
///#define SINGSMALL (1E-3)
///#define SINGSMALL (1E-20)
#define SINGSMALL (1000*NUMEPSILON) // must be larger than machine precision to work for outer M_PI boundary!.  1E-14 works, but need some insurance so use 1E-13

/// Hawley uses 0.06283 (0.02Pi)

#define VOLUMEDIFF 0
/// whether to use volume regularization or not

#define DOSTOREPOSITIONDATA 1
/// whether to store X,V,dxdxp in case are expensive to compute

/// whether to store tetrads (uses tons of memory, and so far not needed -- shouldn't be needed in general if do things covariantly)
#define STORETLAB2ORTHO 0

/// which derivative type to use when computing connection coefficients
///#define CONNDERTYPE DIFFGAMMIE
#define CONNDERTYPE DIFFNUMREC // improved now and much more accurate then DIFFGAMMIE in general. However this is too slow to be used when time-dependent metric is cycling near substeps.

/// whether to correct connection so body forces are 0 when p=constant
/// only makes sense for non-higher order scheme
/// does help near the pole to avoid failures
#define CONNMACHINEBODY 1

/// whether connection coefficients are computed as being axisymmetric
/// Speeds-up connection calculation when starting computations
/// currently all setups are axisymmetric metrics, but don't assume list won't diverge.
#define CONNAXISYMM 1


/// if WHICHEOM==WITHNOGDET, then below determines which EOMs get what geometric prefactor.  Notice (as described in phys.c's source_conn() ) that geometry issue applies AFTER additions/subtractions of EOMs (as done by REMOVERESTMASSFROMUU).
#define WHICHEOM WITHGDET
///#define WHICHEOM WITHNOGDET


/// whether to include memory and fill with gdetvol_func()
#define GDETVOLDIFF 0

/// whether to make gdet finite volume-like when reverting to 1D spherical polar coordinates
#define FIXGDETSPC_WHEN_1DRADIAL 1

/// whether to force gdet to be positive everywhere.
#define FORCEGDETPOSITIVE 0

#define MINDT 1.e-20 // minimum dt

#define JONCHECKS 1 // for vel=3 extra checks (standard things)
#define JONCHECKS2 1 // check_pr() crazy thing for vel3 (crazy u^t check and fix -- 2D only)


#define FLOORDIAGS 1 // whether to compute floor diagnostics

/// whether to use analytic g^{\mu\nu} if available
#define ANALYTICGCON 0

#define ANALYTICCONNECTION 0 // whether to use analytic connection
/// only applies to certain metric's and coordinates, see set_grid.c

#define ANALYTICSOURCE 0 // whether to use analytic source function
/// very slow with gcc, some extend pgcc, not a bit problem with icc
/// only applies to certain metric's and coordaintes, see phys.c

/// whether to (with ENO method) avoid boundary zones if doing outflow on that boundary
#define OUTFLOWAVOIDBC 0

#define FLUXDIMENSPLIT PERFECTUNSPLIT
#define A2CDIMENSPLIT PERFECTUNSPLIT

/// whether to bound flux if doing fluxrecon
#define BOUNDFLUXRECON 0

/// whether to bound unew (GODMARK: no separate algorithm yet, so only works for non-outflow boundary conditions) -- only needed if unewisavg==1
#define BOUNDUNEW 0

/// whether to have dq's : only used by second order methods
#define DODQMEMORY 1

/// whether to adjust interpolation based upon boundary conditions
/// defunct.  Not used anymore by ENO scheme and shouldn't be activated for other schemes.
/// LEAVE 0!
#define BOUNDARYINTERPADJUST 0

/// whether to compute \dot{F_r} as in diag.c
///disable the analysis if N1 <= 1 because leads to segfault if F1 == NULL (as is in the case of N1 == 1)
#define COMPUTEFRDOT 0

/// whether to compute faraday and currents and output them to dumps
#define CALCFARADAYANDCURRENTS 0 // FUCK: Some memory leak maybe because Kraken crashes if this is 1.  Didn't crash before, so unsure what's wrong.  core dump gives some info, but not much.

#define WHICHCURRENTCALC CURRENTCALC1


/// whether want faraday at t=0 for dump
#define FARADAYT0 1

/// whether want partial currents for t=0 dump
#define CURRENTST0 1


#define BOUNDPLPR 0
#define NOFLUXCTONX1DN 0


/// boundary condition mnemonics
/// can be reset/added to by user init.h
#define OUTFLOW 0
#define SYMM 1
#define ASYMM 2
#define FIXED 3
#define POLARAXIS 4
#define FIXEDOUTFLOW 5 // means fixed inflow but allows outflow -- basically outflow if no inflow, but if inflow then set values to fixed quantities
#define NSSURFACE 6 // whatever in bounds.c for NS surface
#define PERIODIC 7 // periodic boundary conditions
#define OUTFLOWNOINFLOW 8  //copy if velocity directed outward; force velocity to zero and copy if directed inward (into the grid from the boundary)
#define RAMESHOUTFLOW 9 //OUTFLOW quantities according the Ramesh's power-law solution
#define HORIZONOUTFLOW 10

#define BCEXTRAP 100  //atch: extrapolation into the boundary with high order (can be specified)
#define CYLAXIS 101   //atch: cylindrical axis BC
#define BCEXTRAP_VEL3 102  //atch: the same as BCEXTRAP but extrapolates 3-velocity rather than 4-velocity; important for Hubble flow
#define JETINJECTION 103 //Fixed boundary condition for jet injection surrounded by an outflow condition
#define BCU1EXTRAPOTHERFIXED 104
#define BCEXTRAPCONSTRAINED 105
#define RESCALEOUTFLOW 106
#define RESCALEFIXEDALLOUTFLOWU1 107
#define FIXED_RESCALEOUTFLOWU1 108
#define BONDIMDOTOUTFLOW 109
#define BONDIINTOUTFLOW 110
#define DISKSURFACE 111
#define FREEOUTFLOW 112
#define R0SING 113
#define FIXEDUSEPANALYTIC 200


#define EVOLVECHECKS 0
/// whether to check boundary conditions and limit gamma during advance()

/// whether and which type of fixups to be used
#define FIXUPZONES FIXUP1ZONE

/** algorithmic choices **/

#define HLLBOUNDARY 0
/// use HLL on boundary fluxes


#define FIXUPFLUX 0
/// fix up the flux using fix_flux() in fixup.c
/// caution...
/// GODMARK: Not sure how metric and how emf works around axes.  Not unlike \Omega and such things that are "symmetric" while v^\phi might be antisymmetric?


#define ZEROOUTFLOWFLUX 0
/// 0: don't do anything special, let boundary conditions control fux
/// 1: zero the x1-velocity at the outflow boundary if inflow in flux routine so interpolation is overridden
/// seems to cause problems at outer edge in internal energy

/// seems to cause problems with EOMFFDE at boudaries.
/// GODMARK





/// assumes inner and outer theta boundaries are at 0 and Pi.
#define ZEROPOLEFLUX 0
/// 0: don't do anything special
/// 1: zero polar theta flux since should be 0


/// seems to cause problems at 0 and pi boundaries with FFDE, or even just one boundary???
/// GODMARK

/// REASON: Cannot just set F=0.  The flux of the kinetic energy still has the pressure term even if the velocity into the wall is 0!


/// whether to flip gdet sign over coordinate singularities
/// completely generally, this should be 1 so that \detg is smooth across the axis.  So then standard boundary conditions on primitives give correct non-kinked behavior through polar axis (including for ram pressure flux term).
#define FLIPGDETAXIS 1

/// whether to flip sign of U3,B3 across the pole
///Another thing that might have helped is how I treat the BCs.  Recall I now flip U2,B2,U3,B3.  U2,B2 makes sense so interpolation sees continuous function for (e.g.) a Cartesian flow through the pole.  In addition, as I mentioned before, I flip U3,B3 because in axisymmetry that just gives the same result of a DONOR-like interpolation.  So there's no change.  In addition, there is no EMF on the pole in that case because U2,B2 on the pole itself is zero.  In non-axisymmetry, flow through the axis would lead to a sign flip and singularity right on the pole.  This would lead to a highly dissipative EMF right at the pole.  By flipping U3,B3 I'm choosing to make the region a numerical "core" instead of a sign-changed singularity.  This core will use DONOR, but have no dissipation term across the pole.
///As we discussed, one could modulate U3,B3 by \sin\theta and achieve a higher-order result.  But the flip or modulation is required to avoid a dissipation-dominated result at the pole.  Most generally, some scheme should be capable of arbitrary high order even in SPC, and I'm guessing modulation by \sin\theta is probably the right thing to do given U3,B3\propto \pm 1\theta near the pole when U3,B3 near the pole matters.
#if(0)
#define FLIPU3AXIS 1
#define FLIPB3AXIS 1
/// should always be 1
#define FLIPU2AXIS 1
#define FLIPB2AXIS 1
/// should always be 0
#define FLIPU1AXIS 0
#define FLIPB1AXIS 0

#else
/// See coord.c and "void dxdxprim() with regard to dxdxp[2][1] and gv311 and uu1 and other symmetry issues.
#define FLIPU3AXIS 0
#define FLIPB3AXIS 0
#define FLIPU2AXIS 0
#define FLIPB2AXIS 0
#define FLIPU1AXIS 0
#define FLIPB1AXIS 0
#endif

/// control bounds.tools.c for SPC coordinates polar axis fixups
#define DOPOLEDEATH 0
#define DOPOLESMOOTH 0 // can choose 1, but probably not necessary and generally won't treat flow correctly near pole even if possibly more robust.
#define DOPOLEGAMMADEATH 0



/// if(periodicx3&&(ncpux3>1)&&ISSPCMCOORDNATIVE(MCOORD)) and below is 1, then do polar MPI boundary transfer
#define IF3DSPCTHENMPITRANSFERATPOLE 1 // working fine now that wavespeed bug in fluxctstag.c was fixed, extrapfunc B1,B2 bug fixed, and extrap gdet B3 instead of Bd3 that exaggerates extrapolation near poles and inconsistent with interpolation.  Also using VARTOINTERPFIELD GDETVERSION.




/// whether to rescale interpolation
#define RESCALEINTERP 0
/// 0: don't rescale
/// 1: do rescale

#define RESCALEINTERPFLUXCTSTAG 0
/// 0: don't rescale
/// 1: do rescale

#define BDIRCONT 1
/// 0: don't
/// 1: make field along flux direction continuous

/// whether to use hyperbolic function for cmin and cmax instead of sharp max(0,l,r)
#define HYPERHLL 0


/// whether to shift stencil inside horizon
#define HORIZONSUPERFAST 0



///////////////////////////////////
///
/// which variable to interpolate
///
//////////////////////////////////


///#define VARTOINTERP PRIMTOINTERP_VSQ
#define VARTOINTERP PRIMTOINTERP

///#define VARTOINTERPFIELD NOSPECIALFIELD
#define VARTOINTERPFIELD GDETVERSION // most consistent with fluxctstag.c and standard extrapfunc in bounds.tools.c



//////////////////////////////////////////////////////////
///
/// STUFF FOR INTERPLINE.C
///
//////////////////////////////////////////////////////////

/// \detg rho u^i used or not (then just u^i or primvel^i)
#define VLINEWITHGDETRHO 1

/// whether Pline in interpline.c for PARALINE without weno has b^2/2 term in pressure (should!)
#define PLINEWITHFIELD 1

/// whether PARALINE uses MONO or not
#define PARALINEUSESMONO 0


/// whether to use appropriately centered primitive
/// GODMARK: averaging is diffusive, so diffuses shock -- may make indicator weak
#define USEAVGPRIMITIVEFORWENOFLAT 1


/// whether to invert average conserved quantity to get a primitive to use as shock indicator
/// then correctly positioned in time, but more diffused value
/// may be a bad idea since inversion is expensive, may find no solution, or may lead to negative internal energy in shock.
#define USEPRIMITIVEFROMAVGCONSERVED 0


/// whether to use contact discontinuity indicator to steepen contacts in interpline.c
#define CONTACTINDICATOR 0

/// send Sasha dP
#define COMPUTEDRHODP 1


/// whether to reduce to 1-point stencil in case of superfast divergence at that point
#define SUPERFASTDIVREDUCE 0 //atch

/// minimum preferred order
#define MINPREFORDER 3

/// whether to compute shock indicator
#define SHOCKINDICATOR 1

/// interp general stuff
#define WHICHPARA PARA4


#define BONDI_BOUNDARY_SET_PL_PR 0  //do not analytically set p_l & p_r at the outer boundary for the Bondi problem




#define NUMPANALYTICOTHER 0
#define DODUMPOTHER 0 // whether to dump other stuff

/// whether to do flux dumps and have memory for it
/// 0: no flux dump files
/// 1 : flux dump files with decomposition by part of LAXF flux
/// 2 : flux dump files with decomposition by physical terms for full LAXF flux
/// decomposition of only needed total fluxes dumped along with fieldline file
#define FLUXDUMP 0

/// default number of extra things needed
#define NUMPOTHER 0





/// delta is simply how big the differencing is, should be small, but not so small to lead to errors due to erros in the metric itself (i.e. keep larger than machine precision)
///#define DELTA (NUMEPSILON*1000.0)
///#define DELTA 1.e-5
///#define DELTA 1.e-5
///#define DELTA (pow(NUMEPSILON,1.0/3.0))
///#define DELTA NUMSQRTEPSILON // as in NR's fdjac() -- smaller isn't always better
/// how to generically set this?  Too high, even slightly (10^{-10} for long doubles) and connection is screwed)

/// Avery mentions that long double trig. functions only return double precision answer.  see ~/research/utils/triglongdouble.c

#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
///#define CONNDELTA (dx[1])
#define CONNDELTA (MY1EM5) // default -- seems to work pretty good generally to reduce max error
///#define CONNDELTA 5.4E-5 // default -- seems to work pretty good
///#define CONNDELTA 4.6E-5 // min of error for a specific case, but apparently not generally good
#elif(REALTYPE==LONGDOUBLETYPE)
///#define CONNDELTA 7.17E-6 // based on min of error for specific case
#define CONNDELTA (MY1EM5) // based on min of error for specific case
/// polar region likes 6.5E-8 (min of error for specific case)
#endif
/// GODMARK: Should really choose CONNDELTA as (totalsize[jj]*dx[jj]*MY1EM5)


/// problem-dependent code activation
#define USERERSETREGION 0


/// whether poledeath in bounds.tools.c keeps sigma constant in some cases
#define BCSIGMACONSTATPOLE 0

/// whether to do one-step diag_fixup accounting
/// for dU (failfloordudump) best to have as 1.  Counting in each function still done.
#define DOONESTEPDUACCOUNTING 1

#define FIELDLINEGDETB 0

#define DO_ASSERTS 0










///////////////////////////////////////////////
///
/// this for bounds.tools.c that user might want to change
///
///////////////////////////////////////////////
#define ADJUSTFLUXCT 0 // whether to adjust fluxCT

#define ADJUSTFLUX 0 // whether to adjust normal flux

/// GODMARK: something seriously wrong with OUTEREXTRAP=1 (EOMFFDE)

#define DEBUGINOUTLOOPS 0

#define OUTEREXTRAP 3
/// 0: just copy
/// 1: gdet or other extrapolation
/// 2: copy (with rescale())
/// 3: Treat same as HORIZONEXTRAP==3

/// how to bound near horizon
#define HORIZONEXTRAP 3
/// as above and also:
/// 3: jon version

/// number of iterations to get v^\phi\sim constant
#define NUMITERVPHI 5


///////////////////////////////
/// number of zones to use pole crushing regularizations
/// to help protect the pole from death blows to the computational grid
/// a sort of crushing regularization
/// causes problems with stability at just beyond pole
/// for field line plots, can just set B^\theta=0 along pole
#define POLEDEATH (MIN(DOPOLEDEATH,N2BND)) // with expansion by 1 point if detects jumps in densities or Lorentz factor (see poldeath())
///#define MAXPOLEDEATH N2BND // can't be larger than N2BND
#define MAXPOLEDEATH (MIN(DOPOLEDEATH+1,N2BND)) // can't be larger than N2BND
#define DEATHEXPANDAMOUNT 0

#define POLEINTERPTYPE 3 // 0=set uu2=bu2=0, 1=linearly interpolate uu2,bu2  2=interpolate B_\phi into pole  3 =linearly for uu2 unless sucking on pole





///////////////////////////////////////////
/// number of zones to enforce Lorentz factor to be small
/// notice that at pole uu1 and uu2 are artificially large and these regions can lead to runaway low densities and so even higher uu1,uu3
/// problem with POLEGAMMADEATH is that at large radius as fluid converges toward pole the fluid stagnates and can fall back at larger angles for no reason -- even for simple torus problem this happens when GAMMAPOLE=1.001
#define POLEGAMMADEATH (MIN(DOPOLEGAMMADEATH,N2BND))
/// maximum allowed Lorentz factor near the pole (set to something large that should be allowed by solution -- problem and grid dependent)
///#define GAMMAPOLE (2.0)

#define GAMMAPOLEOUTGOING 1.1 // keep low
#define GAMMAPOLEINGOINGOUT 1.3 // keep low
#define GAMMAPOLEOUTGOINGPOWER 0.1 // keep low unless know will grow alot
#define GAMMAPOLEOUTGOINGRADIUS 6.0 // very model dependent
#define GAMMAPOLEOUTGOINGRADIUSIN 4.0 // very model dependent
#define GAMMAPOLEOUTGOINGMAX (5.0)
#define GAMMAPOLEINGOING GAMMAMAX

#define GAMMARADPOLEOUTGOING 1.3 // keep low
#define GAMMARADPOLEINGOINGOUT 1.7 // keep low
#define GAMMARADPOLEOUTGOINGPOWER 0.5
#define GAMMARADPOLEOUTGOINGRADIUS 6.0 // very model dependent
#define GAMMARADPOLEOUTGOINGRADIUSIN 4.0 // very model dependent
#define GAMMARADPOLEOUTGOINGMAX (3.0)
#define GAMMARADPOLEINGOING GAMMAMAXRAD

/// factor by which to allow quantities to jump near pole
#define POLEDENSITYDROPFACTOR 5.0
#define POLEGAMMAJUMPFACTOR 2.0

#define RADIUSMOREDEATH (300.0)

////////////////////////////////////////////
/// whether to average in radius for poledeath
#define AVERAGEINRADIUS 0 // not correct  across MPI boundaries since have to shift near boundary yet need that last cell to be consistent with as if no MPI boundary // OPENNPMARK: Also not correct for OpenMP
#define RADIUSTOSTARTAVERAGING 7 // should be beyond horizon so doesn't diffuse across horizon
#define RADIUSTOAVOIDRADIALSUCK (2.0*Rhor)


/// whether if doing full special 3d (i.e. special3dspc==1) that should only do poledeath for inflow
/// 0 : no limit
/// 1 : limit to poledeath acting if radial inflow
/// 2 : limit to poledeath acting on flow within r=RADIUSLIMITPOLEDEATHIN
/// 3 : limit if inflow OR out to r=RADIUSLIMITPOLEDEATHIN (in case inflow only starts near horizon, still poledeath out to that radius
#define IFLIMITPOLEDEATH 0

// whether to limit poledeath to inside horizon and inside death region
#define LIMITEDPOLEDEATHINRADIUS 1


/// radius within which to use poledeath if have IFLIMITPOLEDEATH==3
#define RADIUSLIMITPOLEDEATHIN (3.0) // choose r=3M since always close to BH but always slightly outside horizon to help control stability.

/// how many zones to use poledeath at outer *physical* edge
#define IFLIMITPOLEDEATHIOUT (-100)



////////////////////////////////////////////
/// number of zones to smooth pole
#define POLESMOOTH (MIN(DOPOLESMOOTH,N2BND))

// defaults for bounds.tools.c
#define OUTERDEATH 0
#define OUTERDEATHRADIUS (1E3)
#define OUTERDEATHGAMMAMAX (6.0)
#define OUTERDEATHGAMMAMAXRAD (GAMMAMAXRAD)

// whether to use old non-continuous method for Bi along i interpolation in staggered field method
#define OLDNONCONT 0

/// whether for BH problems to have inner and outer radial boundary as HLL superfast
#define OUTERRADIALSUPERFAST 0

/// whether to do comptonization (upscatter or downscatter)
#define DOCOMPTON 0
// whether to (during utoprimgen()) to do mhd inversion, recompute umhd, use energy conservation to get urad, and then invert that urad->prad.  Doesn't work if mhd>>rad when non-QTYPMHD methods used.
#define ENFORCEMHDCONS2RADCONS 0
#define DOPERF 0 // whether to enable performance opts, including enforcing only QTYPMHD method.
#define ENSURECONS 0 // when done with full step, whether to ensure energy conservation by modifying rad result if possible
#define BORROWENTROPY 0 // whether to try borrowing energy/entropy from radiation to get total energy conservation with entropy method
