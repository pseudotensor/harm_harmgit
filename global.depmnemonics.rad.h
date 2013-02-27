


///////////////////
//
// Mathematical or Methods constants (not physical constants)
//
///////////////////
#define EFLOOR SMALL

// KORALTODO: The below need to be chosen intelligently
#define RADEPS (1.e-6)
#define RADCONV (1.e-7)
#define PRADEPS (1.e-6)
#define PRADCONV (1.e-8)
#define Pi (M_PI)     



///////////////
//
// Some physical constants
//
///////////////
#define GGG (6.674e-8) // cgs in cm^3/(kg s^2)
#define CCCTRUE (2.99792458e10) // cgs in cm/s
#define MSUN (1.989E33) //cgs
#define ARAD (7.56593E-15) // cgs in erg/(K cm^3)
#define K_BOLTZ (1.3806488e-16) // cgs in erg/K
#define M_PROTON (1.67262158e-24) // proton mass in cgs in grams
#define MB (1.66054E-24) // = 1/N_A = 1/(Avogadro's number) = baryon mass in cgs in grams (as often used in general EOSs)
#define SIGMA_RAD (5.67e-5) // cgs in erg/(cm^2 s K^4)

/////////////////////
//
// derived constants
//
/////////////////////
#define MSUNCM (GGG*MSUN/(CCCTRUE*CCCTRUE)) // Msun in cm


//////////////////////////
//
// Define Length, Time, and Mass (and Temperature) units
//
////////////////////////////
#define LBAR (GGG*MPERSUN*MSUN/(CCCTRUE*CCCTRUE)) // cgs
#define TBAR (LBAR/CCCTRUE) // cgs
#define VBAR (LBAR/TBAR) // cgs // HARM requires this be CCCTRUE!
#define RHOBAR (1.0) //cgs so density unit is 1g/cm^3
#define MBAR (RHOBAR*LBAR*LBAR*LBAR) // cgs
#define UBAR (RHOBAR) // energy density is in g/cm^3 so no $c$ floating in code
#define TEMPBAR (M_PROTON*CCCTRUE*CCCTRUE/K_BOLTZ) // cgs unit of temperature in Kelvin (used to make Kelvin dimensionless)
// so for (e.g.) ideal gas, ucode = rhocode * Tcode using u=\rho_0 k_b T / (m_b c^2)  gives both ucode and rhocode in g/cm^3


////////////////
// Note: Interaction requires specific mass unit, which determines code value of ARAD, KAPPAs, and how initial conditions (IC) are chosen.
/////////////////


////////////////////////
//
// Code versions of constants required for radiation interaction
//
////////////////
#define ARAD_CODE (ARAD/K_BOLTZ*LBAR*LBAR*LBAR*UBAR) // so ucode=aradcode*Tcode^4 gives ucode in UBAR units


//////////////
//
// Physical choices for plasma composition
//
//////////////

#define XFACT (1.0) // hydrogen mass fraction
#define ZATOM (1.0) // Atomic number
#define AATOM (1.0) // Atomic weight
#define MUE (1.0) // n_b = rho/m_b = n_e/Y_e = n_e \mu_e
#define MUI (1.0) // n_I = rho/m_b = n_i/Y_i = n_i \mu_i
// MUE,MUI=1,1 for pure hydrogen


////////////////////
//
// Opacities in code units
//
/////////////////////

#define OPACITYBAR (LBAR*LBAR/MBAR) // cgs in cm^2/g
// non-relativistic ES:
#define KAPPA_ES_CODE(rhocode,Tcode) (0.2*(1.0+XFACT)/OPACITYBAR)
#define KAPPA_FF_CODE(rhocode,Tcode) (1.0E23*ZATOM*ZATOM/(MUE*MUI)*(rhocode*RHOBAR)*pow(Tcode*TEMPBAR,-7.0/2.0)/OPACITYBAR)
#define KAPPA_BF_CODE(rhocode,Tcode) (1.0E25*ZATOM*(1.0+XFACT)*(rhocode*RHOBAR)*pow(Tcode*TEMPBAR,-7.0/2.0)/OPACITYBAR)



