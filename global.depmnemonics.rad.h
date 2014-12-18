
/*! \file global.depmnemonics.rad.h
    \brief RADIATION code definitions of dependent quantities

    // MNEMONICS or other things that should rarely change, or things that on depend on the above items (e.g. if statements, loops, etc.)
    
*/


////////////////
///
/// Some physical constants (here in dep since might want to change cTILDA and gTILDA per problem)
///
////////////////
#define GGG (GGG0/gTILDA) // cgs in cm^3/(kg s^2)
#define CCCTRUE (CCCTRUE0/cTILDA) // cgs in cm/s

//////////////////////
///
/// derived constants
///
//////////////////////
#define MSUNCM (GGG*MSUN/(CCCTRUE*CCCTRUE)) // Msun in cm

///////////////////////////
///
/// Define Length, Time, and Mass (and Temperature) units
///
///  (depends upon MPERSUN set by user)
///
/////////////////////////////
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
/// Below to be used in formula for Urad_{code} = arad_{code} * Tcode^4
#define ARAD_CODE_DEF (ARAD*(TEMPBAR*TEMPBAR*TEMPBAR*TEMPBAR)/UBAR)


//////////////
//
// Physical choices for plasma composition
//
//////////////

// Sun-like stars: X=0.75, Y=0.23, Z=0.02, so may be applicable to X-ray binaries
// Solar: X=0.70 Y=0.28 Z=0.02 by mass, neutral MUMEAN=1.3, ionized MUMEAN=0.62
#define XSOLAR (0.7)
#define YSOLAR (0.28)
#define ZSOLAR (1.0-(XSOLAR+YSOLAR))
// Big-Bang: X=0.77 Y=0.23 Z=0.00
// post-Big-Bang: X=0.6 Y=0.4 Z=0, so applicable to GRBs from Pop3 stars.

#define XFACT (0.7) // Hydrogen mass fraction
#define YFACT (0.28) // Helium mass fraction
#define AVG10AJ (1.0/15.5) // for solar composition.  Not used unless MUMEAN is chosen to be for neutrals.


#define ZFACT (1.0-(XFACT+YFACT)) // Rest.  Called "metals".
// http://www.astro.wisc.edu/~townsend/resource/teaching/astro-310-F08/21-eos.pdf
#define MUELE (2.0/(1.0+XFACT)) //(1.0/(XFACT + 0.5*(YFACT+ZFACT))) // inverse of Ye = electron fraction by baryon mass
#define YELE (1.0/MUELE)
#define MUMEANNEUTRAL (1.0/(XFACT + 0.25*YFACT + AVG1OAJ*ZFACT))
#define MUMEANIONIZED (1.0/(2.0*XFACT + 0.75*YFACT + 0.5*ZFACT))
#define MUMEAN (MUMEANIONIZED) // ASSUMPTION: fully ionized // CHOICE

/// while avoids singular behavior, can make inversion unable to reach solution and get locked in cycles due to bad Jacobian, etc.
//#define TEMPMINKELVIN (1.0E+2) // Kelvin // Problem with consistency in error function and entropy estimate for URAD method.
#define TEMPMINKELVIN (1.0E-10) // Kelvin
#define TEMPMIN (TEMPMINKELVIN/TEMPBAR) // Code units


////////////////////
//
// Opacities in code units
//
/////////////////////

/// whether to do comptonization (upscatter or downscatter)
#define DOCOMPTON 1


#define OPACITYBAR (LBAR*LBAR/MBAR) // cgs in cm^2/g


/// PURE ELASTIC SCATTERING
//#define KAPPA_ES_KNCORRF(f) (0.75*((-1.*(1. + 3.*(f)))/Power(1. + 2.*(f),2) +  (0.5*Log(1. + 2.*(f)))/(f) + ((1. + (f))*((2. + 2.*(f))/(1. + 2.*(f)) - (1.*Log(1. + 2.*(f)))/(f)))/Power((f),2)))
//#define KAPPA_ES_KNCORR(rhocode,Tcode) (KAPPA_ES_KNCORREP(K_BOLTZ*(Tcode)*TEMPBAR/(MELE*CCCTRUE*CCCTRUE)))
#define KAPPA_ES_FERMICORR(rhocode,Tcode) (1.0/(1.0+2.7E11*((rhocode)*RHOBAR)/pow((Tcode)*TEMPBAR,2.0))) // Buchler and Yueh 1976 (just Fermi part). Fewer electrons when near Fermi fluid limit.
#define KAPPA_ES_KNCORR(rhocode,Tcode) (1.0/(1.0+pow((Tcode)*TEMPBAR/4.5E8,0.86)))  // Buchler and Yueh 1976 .  Klein-Nishina for thermal electrons.
/// kappaes = sigma_T n_e = sigma_T n_b (n_e/n_b) = sigma_T rho/mb (ne/nb)
#define KAPPA_ES_CODE(rhocode,Tcode) (0.2*(1.0+XFACT)*KAPPA_ES_FERMICORR(rhocode,Tcode)*KAPPA_ES_KNCORR(rhocode,Tcode)/OPACITYBAR)

// INELASTIC COMPTON TERM
// see DOCOMPTON in physics.tools.rad.c:
// in term2, doesn't change photon energy, so that in scattering-dominated atmospheres, photons move through unchanged by temperature of gas.
// Eq A7 in http://adsabs.harvard.edu/abs/2012ApJ...752...18K used first in http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0904.4123 as based upon a calculation in http://adsabs.harvard.edu/abs/2000thas.book.....P .
// Also interesting for next steps:
//1) Conservative form of Kompaneet's equation and dealing with the diffusion term implicitly: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.13.798 and related: http://www.osti.gov/scitech/biblio/891567 .  We should ensure the 4-force is consistent (numerically and analytically) with what's used in this equation for n.
//2) Relativistic corrections:  http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1201.5606 and http://adsabs.harvard.edu/abs/2002nmgm.meet.2329S .  It looks like nothing more difficult as far as actually using the expressions in place of non-relativistic version.
//One semi-relevant application: http://www.aanda.org/articles/aa/full_html/2009/45/aa12061-09/aa12061-09.html
#define KAPPA_FORCOMPT_RELCORREP(ep) ((1.0 + 3.683*(ep)+4.0*(ep)*(ep))/(1.0 + (ep))) // Sadowski et al. (2014) Eq 26 and 27.
#define KAPPA_FORCOMPT_RELCORR(rhocode,Tcode) (KAPPA_FORCOMPT_RELCORREP(K_BOLTZ*(Tcode)*TEMPBAR/(MELE*CCCTRUE*CCCTRUE)))
#define KAPPA_FORCOMPT_CODE(rhocode,Tcode) (0.2*(1.0+XFACT)*KAPPA_ES_FERMICORR(rhocode,Tcode)*KAPPA_FORCOMPT_RELCORR(rhocode,Tcode)/OPACITYBAR)

#define TEMPELE (MELE*CCCTRUE*CCCTRUE/K_BOLTZ)

/// EMISSION (Tr=Tg) or ABSORBPTION (Tr different from Tg)
#define KAPPA_ZETA(Tgcode,Trcode) ((TEMPMIN+Trcode)/(TEMPMIN+Tgcode))
#define KAPPA_FF_CODE(rhocode,Tgcode,Trcode) (4.0E22*(1.0+XFACT)*(1.0-ZFACT)*((rhocode)*RHOBAR)*pow((Tgcode)*TEMPBAR,-0.5)*pow((Trcode)*TEMPBAR,-3.0)*log(1.0+1.6*KAPPA_ZETA(Tgcode,Trcode))*(1.0+4.4E-10*(Tgcode*TEMPBAR))/OPACITYBAR)  // ASSUMPTION: Thermal ele and no pairs.  See Rybicki & Lightman Eq~5.25 and McKinney & Uzdensky (2012) .  For Tr,Tg split, see Ramesh notes.
#define KAPPA_BF_CODE(rhocode,Tgcode,Trcode) (3.0E25*(ZFACT)*(1.0+XFACT+0.75*YFACT)*((rhocode)*RHOBAR)*pow((Tgcode)*TEMPBAR,-0.5)*pow((Trcode)*TEMPBAR,-3.0)*log(1.0+1.6*KAPPA_ZETA(Tgcode,Trcode))/OPACITYBAR) // ASSUMPTION: Number of electrons similar to for solar abundances for 1+X+(3/4)Y term.  For Tr,Tg split, see Ramesh notes.
#define KAPPA_CHIANTIBF_CODE(rhocode,Tgcode,Trcode) (4.0E34*((rhocode*RHOBAR))*(ZFACT/ZSOLAR)*YELE*pow((Tgcode)*TEMPBAR,-1.7)*pow((Trcode)*TEMPBAR,-3.0)/OPACITYBAR) // *XFACT literally from Fig 34.1 in Draine book, but for solar n_H\sim n_b\sim 1/cm^3 only
#define KAPPA_HN_CODE(rhocode,Tgcode,Trcode) (1.1E-25*pow(ZFACT,0.5)*pow((rhocode)*RHOBAR,0.5)*pow((Tgcode)*TEMPBAR,7.7+3.0)*pow((Trcode)*TEMPBAR,-3.0)/OPACITYBAR)
#define KAPPA_MOL_CODE(rhocode,Tgcode,Trcode) (0.1*ZFACT/OPACITYBAR)
// see opacities.nb
#define KAPPA_GENFF_CODE(rhocode,Tgcode,Trcode) (1.0/(1.0/(KAPPA_MOL_CODE(rhocode,Tgcode,Trcode)+KAPPA_HN_CODE(rhocode,Tgcode,Trcode)) + 1.0/(KAPPA_CHIANTIBF_CODE(rhocode,Tgcode,Trcode)+KAPPA_BF_CODE(rhocode,Tgcode,Trcode)+KAPPA_FF_CODE(rhocode,Tgcode,Trcode)))) // for 1.3E3K \le T \le 1E9K or higher.  Numerically better to have kappa bottom out at low T so no diverent opacity as T->0

// whether to allow kappa to depend explicitly upon position, which would require getting position and can be expensive.
#define ALLOWKAPPAEXPLICITPOSDEPENDENCE 0

//////////////////////////////////
//
// Notes on differences with Koral setup
//
// 1) Units:
// 1a) Koral sets mass scale as M (i.e. black hole mass), where I set it via the density scale of 1g/cm^3
// 1b) Koral's M is mass per unit Msun, just like my MPERSUN, but as 1a, my mass scale is independent of MPERSUN
// 1c) RADOUTPUTINZAMO sets output of koral in orthonormal basis in zamo frame.
// 1d) Koral hard codes MASSCM, which makes it hard to match for similar constants when gTILDA is changed.
// 2) Algorithmic differences:
// 2a) Koral currently uses v=c as the wavespeed limit when including radiation
// 2b) Koral (and now koralinsert) updates radiation source term with dUgeom+dUriemann.  Required for high opacity cases (e.g. RADPULSEPLANAR) to avoid diffusion.  See notes elsewhere.
// 2c) Koral typically uses MINM with theta=1, but uses MP5 when needing more accuracy
// 3) Details:
// 3a) Koral's F and U in LAXF have no \detg
// 4) Runtime differences:
// 4a) Koral is defaulted to use OpenMP , and default is full use of all cores, so any output from parallel loops will be randomly ordered.  remove -fopenmp to get linear debug output.
//
//////////////////////////////////


////////////////////////
//
// Koral units conversions
//
////////////////
// Olek keeps changing gTILDA and cTILDA in different code versions, and so actually changing the problems despite PROBLEM files wanted to describe specific problem as in (e.g.) the paper.
#define MASSCM (GGG*MPERSUN*MSUN/pow(CCCTRUE,2.0))
#define KORAL2HARMRHO(rho) ((rho/GGG*CCCTRUE*CCCTRUE/MASSCM/MASSCM)/RHOBAR)






//
//At this point actually everything seems to be working very well qualitatively.  Some effort went into adding tests, fixed the boundary conditions, adding all the missing pieces of code (opacity limiter), getting u2p_rad() to be more generally robust via HARM-like fixups (both inside u2p_rad() and using fixup_utoprim()), coding-up tetrad stuff that might be used, getting units to make sense, debugging the opacity and source terms, folding koral into HARM's general EOS framework, etc.  Check the git log for some other details.
//
//Many things left to do:
//
//0) Setup all the rest of the problems.  Now that code actually works (Hurray!) and is no longer incomplete (whew), we can add the other test problems.  We also need to think of MHD test problems for the new paper.
//
//1) See KORALTODO or SUPERGODMARK anywhere, and you'll see an issue.
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
//10) #9 is only possible if the equations are written to avoid catastrophic cancellation.  Maybe G and other terms have catastrophic cancellation issues.  E.g., kappaes cancels in Gu for static flow, but maybe other actual serious cancellation issues somewhere.  Yes, catastrophic issue in u2p_rad() for non-rel limit is fixed.  Maybe issue in ultrarel limit!






















/* Hi guys, */

/* Regarding M1's imperfection, there are many studies on using radiative transfer in planetary and stellar atmospheres, and Shane is probably fully aware of the issues.  It's well known what fundamental things one misses with various approximations.  Even Wikipedia goes through some of those things (http://en.wikipedia.org/wiki/Two-stream_approximation_(radiative_transfer)!  But apart from those papers I mentioned that are astro-centric, there are *many* papers (just google) that go through why the 2-stream approximation is so great and what we would miss. */

/* We'd want to get that effective 2-stream approximation per grid cell, and we are *far* from it.  So we'll miss many basic physics effects.  This includes the ability to treat multiple scatterings when we involve Compton scattering. */

/* But forgetting about scattering, I think the main problem is directional independence in the optically thin limit.  Per-grid cell independence for each face would be a basic freedom.  Here's how the method would work: */

/* 1) In each cell, in the averaged energy-weighted radiation rest frame (average of, e.g., R^t_\mu over bins), we split the angles up into bins/bands.  Easiest is 2N to 4N for N-dimensions giving 8N to 16N new variables.  Only increases total number of variables from current 12 to 24 (or 48) -- so factor of only two (or up to four). */

/* 2) Once the problem has been setup (one now can control (per-cell) multiple radiation directions and energy densities), those are evolved with the flux conservative approach like normal to get the new U.  There's no interaction between bins at the flux level, but radiation might want to change its bin -- but that's not done yet (maybe it's a good idea to rebin here already, but maybe not). */

/* 3) The source term is applied.  The cumulative radiation source still applies to the fluid.  Each radiation bin has its own independent source value determination/step as if the other bins don't exist.  After each internal implicit/explicit source step, the new value of R^t_\mu is compared to the average R^t_\mu and we rebin each bin at each implicit step if the radiation wants to be in a new bin.  This allows for interaction between bins at the source term level (maybe this isn't required, but it probably is). */

/* 4) We do the final U->P and determine if rebinning (e.g. just even due to fluxes) is required by using each bin's R^t_\mu as compared to the average.  If already rebinnned at #2 and #3, then no rebinning required here except as related to how harm or koral updates variables (i.e. harm and koral don't update U during implicit and instead keep it as a dU that later gets added). */


/* This takes the flux step (during which no rebinning yet occurs) and source step (during which rebinning can occur) and doing the inversion (during which rebinning can occur) to get new values for the radiation for each bin.  This allows for direct interaction between beams via the sources during stage 3, or the fluxes evolve the beam to get rebinned during stage 4. */

/* Notes: */

/* 1) The 4-force and u2p_rad() assume isotropic emission in the radiation frame.  We don't violate this because each beam bin is consistent with this.  The only "violation" occurs during rebinning, which just merges or splits beams and conserves energy-momentum -- there's no violation of isotropy per beam when isotropy is assumed. */

/* 2) The only non-linearity for R^t_\mu enters in the source step.  But as far as a solution to the equations of motion that determines whether U->P can be done, by the linearity of the equations, if R^t_\mu is a solution, so is any other sum of other solutions, so we can't run into the sum being a non-solution (i.e. no primitive).  Even analytically, the source step could have a crazy cooling or heating that ends up giving bad U->P, but not the advection+geometry part. */

/* Let me know what you think -- if there are gaps or errors. */
