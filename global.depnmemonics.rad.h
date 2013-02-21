

///////////////
// below 6 should only be used for initial condition setup for HARMUNITS
#define GGG_IC (6.674e-8)
#define CCCTRUE_IC 2.998e10
#define MSUNCM_IC 147700.
#define K_BOLTZ_IC (1.3806488e-16L * GGG_IC / CCCTRUE_IC / CCCTRUE_IC / CCCTRUE_IC / CCCTRUE_IC  )
#define M_PROTON_IC (1.67262158e-24L * GGG_IC / CCCTRUE_IC / CCCTRUE_IC)
#define SIGMA_RAD_IC (5.67e-5 * GGG_IC / CCCTRUE_IC / CCCTRUE_IC / CCCTRUE_IC / CCCTRUE_IC / CCCTRUE_IC * MASSCM_IC * MASSCM_IC  * MASSCM_IC)
//////////////////////////////

/////////////
// Dimensionless things or (e.g. sigma) sets scale for density that would be otherwise invariant without radiation
////////////
#define SIGMA_RAD (1.0) // KORALTODO: HARMUNITS?  SUPERGODMARK.
//efine A_RAD (4.*5.67e-5/CCCTRUE * GGG / CCCTRUE / CCCTRUE / CCCTRUE / CCCTRUE)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     

#define EFLOOR SMALL
#define RADEPS (1.e-6)
#define RADCONV (1.e-7)

#define PRADEPS (1.e-6)
#define PRADCONV (1.e-8)




#define MSUNCMFAKE (1.0)
#define CCCFAKE (1.0)
#define GGGFAKE (1.0)

#ifndef MASS
#define MASS 1./MSUNCMFAKE //default mass of the BH used to calibrate radiation constant, Solar mass units
#endif

#define MASSCM (MASS*MSUNCMFAKE) //mass in cm

#define LCM (MASSCMFAKE) //unit of length in cm

#define TSEC (MASSCMFAKE/CCCFAKE) //unit of time in seconds

#define GMC2CM (MASSCMFAKE) //gravitational radius in cm


#define KAPPA_ES_COEFF_IC (kappaCGS2GU(0.4))
#define KAPPA_FF_COEFF_IC (1.7e-25/1.67262158e-24/1.67262158e-24*CCCTRUE_IC*CCCTRUE_IC*CCCTRUE_IC*CCCTRUE_IC/GGG_IC/GGG_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC)
#define KAPPA_BF_COEFF_IC (4.8e-24/1.67262158e-24/1.67262158e-24*CCCTRUE_IC*CCCTRUE_IC*CCCTRUE_IC*CCCTRUE_IC/GGG_IC/GGG_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC/MASSCM_IC)



// KORALTODO HARMUNITS : Needs to be figured out SUPERGODMARK
#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))
#define KAPPA_FF_COEFF (1.7e-25/1.67262158e-24/1.67262158e-24*CCCFAKE*CCCFAKE*CCCFAKE*CCCFAKE/GGGFAKE/GGGFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE)
#define KAPPA_BF_COEFF (4.8e-24/1.67262158e-24/1.67262158e-24*CCCFAKE*CCCFAKE*CCCFAKE*CCCFAKE/GGGFAKE/GGGFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE/MASSCMFAKE)
