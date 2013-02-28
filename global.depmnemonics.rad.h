
//////////////////////////
//
// Define Length, Time, and Mass (and Temperature) units
//
//  (depends upon MPERSUN set by user)
//
////////////////////////////
#define LBAR (GGG*MPERSUN*MSUN/(CCCTRUE*CCCTRUE)) // cgs in cm
#define TBAR (LBAR/CCCTRUE) // cgs in s
#define VBAR (LBAR/TBAR) // cgs // HARM requires this be CCCTRUE!
#define RHOBAR (1.0) //cgs so density unit is 1g/cm^3
#define MBAR (RHOBAR*LBAR*LBAR*LBAR) // cgs in grams
#define ENBAR (MBAR*VBAR*VBAR) // cgs energy in ergs
#define UBAR (RHOBAR*VBAR*VBAR) // cgs energy density in ergs/cm^3
#define TEMPBAR (M_PROTON*CCCTRUE*CCCTRUE/K_BOLTZ) // cgs unit of temperature in Kelvin (used to make Kelvin dimensionless)
// so for (e.g.) ideal gas, ucode = rhocode * Tcode using u=\rho_0 k_b T / (m_b c^2)  gives both ucode and rhocode in g/cm^3

// NOTEMARK: If IC has RHO~UU in terms of units, then need to use RHOBAR to make each dimensionless.  Else add explicit CCCTRUE^2 prefactor in IC value of UU (i.e. energy density) before normalizing with UBAR.

////////////////
// Note: Interaction requires specific mass unit, which determines code value of ARAD, KAPPAs, and how initial conditions (IC) are chosen.
/////////////////


////////////////////////
//
// Code versions of constants required for radiation interaction
//
////////////////
// Below to be used in formula for Urad_{code} = arad_{code} * Tcode^4
#define ARAD_CODE (ARAD*(TEMPBAR*TEMPBAR*TEMPBAR*TEMPBAR)/UBAR)


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




//////////////////////////////////
//
// Notes on differences with Koral setup
//
// 1) Koral sets mass scale as M (i.e. black hole mass), where I set it via the density scale of 1g/cm^3
// 2) Koral hard codes MASSCM, which makes it hard to match for similar constants when gTILDA is changed.
// 3) Koral currently uses v=c as the wavespeed limit when including radiation
// 4) Koral is defaulted to use OpenMP , and default is 4 cores, so any output from parallel loops will be randomly ordered.
// 5) Koral is not fully unsplit for forces at each step.  Koral updates primitives with geometry+MHD sources before applying radiation forces.  This means never will be able to reach exact equilibrium. Note that this is different from operator split.  While the code can be operator split, this doesn't mean one has to use new primitives in each operation as in koral or zeus.  One can remove the "calc_prmitives() before the source update to ensure unsplit (this is what's done for implicit_lab but for some reason not explicit_lab by default)
// 6) Koral typically uses MINM with theta=1
// 7) Koral has vectors in an orthonormal basis
// 8) Koral's F and U in LAXF have no \detg
//
//////////////////////////////////
