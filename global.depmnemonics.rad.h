
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
#define ARAD_CODE_DEF (ARAD*(TEMPBAR*TEMPBAR*TEMPBAR*TEMPBAR)/UBAR)


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






//
//At this point actually everything seems to be working very well qualitatively.  Some effort went into adding tests, fixed the boundary conditions, adding all the missing pieces of code (opacity limiter), getting u2p_rad() to be more generally robust via HARM-like fixups (both inside u2p_rad() and using fixup_utoprim()), coding-up tetrad stuff that might be used, getting units to make sense, debugging the opacity and source terms, folding koral into HARM's general EOS framework, etc.  Check the git log for some other details.
//
//Many things left to do:
//
//0) Setup all the rest of the problems.  Now that code actually works (Hurray!) and is no longer incomplete (whew), we can add the other test problems.  We also need to think of MHD test problems for the new paper.
//
//1) See KORALTODO or SUPERGODMARK anywhere, and you'll see an issue.
//
//2) Compute_dt() needs to be using information passed from the place where dt for each step is set
//
//3) Implicit use of Uptorprimgen() needs to allow only a few steps to be taken -- no point in getting machine accurate U->p there if taking multiple implicit steps.  But often only take 1-2 implicit steps, so need to be careful.
//
//4) If Utoprimgen fails in implicit, won't converge, so need a backup method.  Koral backup won't work for rel flows, so need to sub-cycle with explicit scheme or use a different semi-implicit scheme like one I mentioned from numerical recipes.
//
//6) Need to check factors of Pi and 1/4 for B in calc_Gu as well as 4\pi in Gu.  Are those really dimensionless and should be there?
//
//7) I'm unsure about Olek's velocity limiter for tau>1.  I'm not sure the 4/3 is right.
//
//8) Not sure if reversions in u2p_rad() are all agreeable.  Better than prior choices, but e.g. CASE2B is a concern.  It's what also caused the implicit method to fail when CASE2B is hit.  Need some backup approach for CASE2B situation.
//
//9) Iteration and other constants as in global.depmnemonics.rad.h need to be chosen intelligently.  Can't always just be 1E-6.  For the MHD inversion, machine precision accuracy is always sought.  Maybe required for radiation sometimes.
//
//10) #9 is only possible if the equations are written to avoid catastrophic cancellation.  Maybe G and other terms have catastrophic cancellation issues.  E.g., kappaes cancels in Gu for static flow, but maybe other actual serious cancellation issues somewhere.  

