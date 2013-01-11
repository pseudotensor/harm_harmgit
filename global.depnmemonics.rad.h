

#define GGG (6.674e-8)
#define CCCTRUE 2.998e10
#define MSUNCM 147700.

#define K_BOLTZ (1.3806488e-16L * GGG / CCCTRUE / CCCTRUE / CCCTRUE / CCCTRUE  )
#define M_PROTON (1.67262158e-24L * GGG / CCCTRUE / CCCTRUE)
#define SIGMA_RAD (5.67e-5 * GGG / CCCTRUE / CCCTRUE / CCCTRUE / CCCTRUE / CCCTRUE * MASSCM * MASSCM  * MASSCM)
//efine A_RAD (4.*5.67e-5/CCCTRUE * GGG / CCCTRUE / CCCTRUE / CCCTRUE / CCCTRUE)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     

#define EFLOOR SMALL
#define RADEPS (1.e-6)
#define RADCONV (1.e-7)

#define PRADEPS (1.e-6)
#define PRADCONV (1.e-8)





#ifndef MASS
#define MASS 1./MSUNCM //default mass of the BH used to calibrate radiation constant, Solar mass units
#endif

#define MASSCM (MASS*MSUNCM) //mass in cm

#define LCM (MASSCM) //unit of length in cm

#define TSEC (MASSCM/CCC) //unit of time in seconds

#define GMC2CM (MASSCM) //gravitational radius in cm


#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))
#define KAPPA_FF_COEFF (1.7e-25/1.67262158e-24/1.67262158e-24*CCCTRUE*CCCTRUE*CCCTRUE*CCCTRUE/GGG/GGG/MASSCM/MASSCM/MASSCM/MASSCM/MASSCM)
#define KAPPA_BF_COEFF (4.8e-24/1.67262158e-24/1.67262158e-24*CCCTRUE*CCCTRUE*CCCTRUE*CCCTRUE/GGG/GGG/MASSCM/MASSCM/MASSCM/MASSCM/MASSCM)
