/* tau_neededbyharm.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
   -lf2c -lm   (in that order)
*/

#include "f2c.h"

/*     f2c -f -P tau_neededbyharm.f ; cp tau_neededbyharm.c tau_neededbyharm.P ~/latestcode */

/*     Currently, ki-rh42 does not have f2c, so use ki-rh39: */
/*     f2c -f -P tau_neededbyharm.f */
/*     scp tau_neededbyharm.c tau_neededbyharm.P jon@ki-rh42:/data/jon/latestcode/harm_harm/ */
/*     tau_neededbyharm.c from f2c: MUST remove static in front of variables! */
/*     Duplicates those operations from tau_calc() that involve things dependent upon hcm so that hcm isn't in
       dependent variable in table */
/* ====================================================================== */
/*     ENSURE this list is the same for how used in tau_calc.f */
/* Subroutine */ int computefinal_fromhcm__(doublereal *clight, doublereal *
                                            mb, doublereal *rhob, doublereal *kbtk, doublereal *hcm, doublereal *
                                            unue0, doublereal *unuebar0, doublereal *unumu0, doublereal *
                                            qtautnueohcm, doublereal *qtautnuebarohcm, doublereal *qtautmuohcm, 
                                            doublereal *qtauanueohcm, doublereal *qtauanuebarohcm, doublereal *
                                            qtauamuohcm, doublereal *nnue0, doublereal *nnuebar0, doublereal *
                                            nnumu0, doublereal *ntautnueohcm, doublereal *ntautnuebarohcm, 
                                            doublereal *ntautmuohcm, doublereal *ntauanueohcm, doublereal *
                                            ntauanuebarohcm, doublereal *ntauamuohcm, doublereal *lambdatot, 
                                            doublereal *lambdaintot, doublereal *tauphotonohcm, doublereal *
                                            tauphotonabsohcm, doublereal *nnueth0, doublereal *nnuebarth0, 
                                            doublereal *qphoton, doublereal *qm, doublereal *graddotrhouye, 
                                            doublereal *tthermaltot, doublereal *tdifftot, doublereal *rho_nu__, 
                                            doublereal *p_nu__, doublereal *s_nu__, doublereal *ynulocal, 
                                            doublereal *ynuthermal, doublereal *ynuthermal0, doublereal *enu, 
                                            doublereal *enue, doublereal *enuebar)
{
  /* System generated locals */
  doublereal d__1;

  /* Local variables */
  doublereal n_nutau0__, u_nutau0__, n_nuebar__, u_nuebar__, 
    n_nuelth__, ntauatau, qtauatau, ntauttau, qtauttau, n_nuebar0__, 
    u_nuebar0__, h__, u_photon0__, nm_nuebar__, qm_nuebar__, 
    ntaua_nue__, qtaua_nue__, ntaut_nue__, qtaut_nue__, tauphoton, nb,
    nm, n_nuebarth__;
  extern doublereal tthermalfromlambda_(doublereal *, doublereal *, 
                                        doublereal *);
  doublereal n_nuebarth0__;
  extern doublereal rate_2stream__(doublereal *, doublereal *, doublereal *,
                                   doublereal *, doublereal *);
  doublereal nmel, qmel, nmmu, qmmu, ntaua_nuebar__, qtaua_nuebar__, 
    ntauatauohcm, qtaut_nuebar__, qtauatauohcm, ntaut_nuebar__, 
    tauphotonabs, ntauttauohcm, qtauttauohcm, n_nue__, small, u_nue__,
    nmtau, qmtau, n_nue0__, u_nue0__, n_nuel__, nm_nue__, qm_nue__, 
    u_nuel__, n_numu__, u_numu__, n_numu0__, u_numu0__, n_nueth__, 
    n_nutau__;
  extern doublereal density_2stream__(doublereal *, doublereal *, 
                                      doublereal *);
  doublereal ntauamu, qtauamu;
  extern doublereal tdifffromlambda_(doublereal *, doublereal *, doublereal 
                                     *);
  doublereal u_nutau__, ntautmu, qtautmu, n_nueth0__;

  /* ====================================================================== */
  /* outputs */
  /* more outputs */
  /*     Passed to here */
  /* functions */
  /*     Locals */
  /* outputs */
  /* more outputs */
  /*     Set some things */
  nb = *rhob / *mb;
  small = 1e-150;
  /*     Set some equivalences */
  h__ = *hcm;
  u_photon0__ = *unumu0;
  u_nue0__ = *unue0;
  u_nuebar0__ = *unuebar0;
  u_numu0__ = *unumu0;
  u_nutau0__ = *unumu0;
  n_nue0__ = *nnue0;
  n_nuebar0__ = *nnuebar0;
  n_numu0__ = *nnumu0;
  n_nutau0__ = *nnumu0;
  n_nueth0__ = *nnueth0;
  n_nuebarth0__ = *nnuebarth0;
  qtauttauohcm = *qtautmuohcm;
  qtauatauohcm = *qtauamuohcm;
  /*     Non-optical depth version of Ynuthermal: */
  *ynuthermal0 = (n_nueth0__ - n_nuebarth0__) / nb;
  /*     Set 2-stream approximation solutions */
  /* ccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*     Photon emission energy-volume-rate */
  tauphoton = *tauphotonohcm * *hcm;
  tauphotonabs = *tauphotonabsohcm * *hcm;
  *qphoton = rate_2stream__(clight, hcm, &u_photon0__, &tauphoton, &
                            tauphotonabs);
  /* ccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*     Neutrino emission energy-volume-rate */
  qtaut_nue__ = *qtautnueohcm * *hcm;
  qtaua_nue__ = *qtauanueohcm * *hcm;
  qtaut_nuebar__ = *qtautnuebarohcm * *hcm;
  qtaua_nuebar__ = *qtauanuebarohcm * *hcm;
  qtautmu = *qtautmuohcm * *hcm;
  qtauamu = *qtauamuohcm * *hcm;
  qtauttau = qtauttauohcm * *hcm;
  qtauatau = qtauatauohcm * *hcm;
  /*     energy rates */
  qm_nue__ = rate_2stream__(clight, hcm, &u_nue0__, &qtaut_nue__, &
                            qtaua_nue__);
  qm_nuebar__ = rate_2stream__(clight, hcm, &u_nuebar0__, &qtaut_nuebar__, &
                               qtaua_nuebar__);
  qmel = qm_nue__ + qm_nuebar__;
  /*      write(*,*) 'qminus',qminus_nue,Qm_nue,qtaut_nue,qtaua_nue */
  /*      write(*,*) 'qminus2',u_nue0,etanu,hcm */
  qmmu = rate_2stream__(clight, hcm, &u_numu0__, &qtautmu, &qtauamu);
  qmtau = rate_2stream__(clight, hcm, &u_nutau0__, &qtauttau, &qtauatau);
  /*     does NO LONGER include photons!  But never use in any other way except as total cooling */
  /*     e.g. might want to compute Qm/Nm for energy in HARM */
  *qm = qmel + qmmu + qmtau;
  /* ccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*     Neutrino emission number-volume-rate */
  /*     Define similar quantities */
  ntauttauohcm = *ntautmuohcm;
  ntauatauohcm = *ntauamuohcm;
  ntaut_nue__ = *ntautnueohcm * *hcm;
  ntaua_nue__ = *ntauanueohcm * *hcm;
  ntaut_nuebar__ = *ntautnuebarohcm * *hcm;
  ntaua_nuebar__ = *ntauanuebarohcm * *hcm;
  ntautmu = *ntautmuohcm * *hcm;
  ntauamu = *ntauamuohcm * *hcm;
  ntauttau = ntauttauohcm * *hcm;
  ntauatau = ntauatauohcm * *hcm;
  /*     number rates (in 1 direction) */
  nm_nue__ = rate_2stream__(clight, hcm, &n_nue0__, &ntaut_nue__, &
                            ntaua_nue__);
  nm_nuebar__ = rate_2stream__(clight, hcm, &n_nuebar0__, &ntaut_nuebar__, &
                               ntaua_nuebar__);
  nmel = nm_nue__ + nm_nuebar__;
  nmmu = rate_2stream__(clight, hcm, &n_numu0__, &ntautmu, &ntauamu);
  nmtau = rate_2stream__(clight, hcm, &n_nutau0__, &ntauttau, &ntauatau);
  nm = nmel + nmmu + nmtau;
  /*     Per volume (of size H^3) change: */
  /*     Per area is obtained as mb*(Nm_nuebar-Nm_nue)*hcm */
  /*     and then this is flux per face */
  /*     Can be positive or negative */
  *graddotrhouye = *mb * (nm_nuebar__ - nm_nue__);
  /* ccccccccccccccccccccccccccccc */

  /*     Compute final neutrino energy and number densities */

  /* ccccccccccccccccccccccccccccc */
  /*     energy densities */
  u_nue__ = density_2stream__(&u_nue0__, &qtaut_nue__, &qtaua_nue__);
  u_nuebar__ = density_2stream__(&u_nuebar0__, &qtaut_nuebar__, &
                                 qtaua_nuebar__);
  u_nuel__ = u_nue__ + u_nuebar__;
  /* diag only */
  u_numu__ = density_2stream__(&u_numu0__, &qtautmu, &qtauamu);
  u_nutau__ = density_2stream__(&u_nutau0__, &qtauttau, &qtauatau);
  /*     number densities */
  n_nue__ = density_2stream__(&n_nue0__, &ntaut_nue__, &ntaua_nue__);
  n_nuebar__ = density_2stream__(&n_nuebar0__, &ntaut_nuebar__, &
                                 ntaua_nuebar__);
  n_nuel__ = n_nue__ + n_nuebar__;
  /* diag only */
  n_numu__ = density_2stream__(&n_numu0__, &ntautmu, &ntauamu);
  n_nutau__ = density_2stream__(&n_nutau0__, &ntauttau, &ntauatau);
  /*     thermal number densities (using standard optical depth rather than thermalizing optical depth) GODMARK 
         -- approximation to avoid tabulating thermalizing optical depths separately -- somewhat inconsistent w
         ith using lambdaintot for Tthermaltot */
  n_nueth__ = density_2stream__(&n_nueth0__, &ntaut_nue__, &ntaua_nue__);
  n_nuebarth__ = density_2stream__(&n_nuebarth0__, &ntaut_nuebar__, &
                                   ntaua_nuebar__);
  n_nuelth__ = n_nueth__ + n_nuebarth__;
  /*     total internal energy density */
  /* diag only */
  *rho_nu__ = u_nuel__ + u_numu__ + u_nutau__;
  /*     total pressure */
  *p_nu__ = *rho_nu__ / 3.;
  /*     entropy density in 1/cc */
  *s_nu__ = *rho_nu__ * 1.3333333333333333 / (*kbtk + small);
  /*     [k_b] form is s_nu/(rho/m_b) is "per baryon" entropy for rho and S given in cgs. */
  /*     Electron neutrino fraction (other species cancel since their chemical potential is 0) */
  /*     Should always be larger than 0.  For convergence process, let smallest be SMALL */
  /*     can be negative due to optical depth supression */
  /* now in 1/cc */
  *ynulocal = (n_nue__ - n_nuebar__) / nb;
  /*     Could now check of Ynulocal is same as expected Ynu */
  /*     Thermalized faction Y_\nu */
  /*     can be negative due to optical depth supression */
  /* here this is just a diagnostic check */
  *ynuthermal = (n_nueth__ - n_nuebarth__) / nb;
  *tdifftot = tdifffromlambda_(clight, &h__, lambdatot);
  *tthermaltot = tthermalfromlambda_(clight, &h__, lambdaintot);
  /*     These are energies of *escaping* neutrinos */
  /*     These energies are used for neutrino annihilation heating rates */
  /* Computing MAX */
  d__1 = *qm / (nm + small);
  *enu = max(d__1,0.);
  /* Computing MAX */
  d__1 = qm_nue__ / (nm_nue__ + small);
  *enue = max(d__1,0.);
  /* Computing MAX */
  d__1 = qm_nuebar__ / (nm_nuebar__ + small);
  *enuebar = max(d__1,0.);
  return 0;
} /* computefinal_fromhcm__ */

/*     Very similar to above but only computes energy densities (i.e. rho_nu, p_nu, s_nu) */
/* ====================================================================== */
/* Subroutine */ int computefinal_justdensities_fromhcm__(doublereal *clight, 
                                                          doublereal *mb, doublereal *rhob, doublereal *kbtk, doublereal *hcm, 
                                                          doublereal *unue0, doublereal *unuebar0, doublereal *unumu0, 
                                                          doublereal *qtautnueohcm, doublereal *qtautnuebarohcm, doublereal *
                                                          qtautmuohcm, doublereal *qtauanueohcm, doublereal *qtauanuebarohcm, 
                                                          doublereal *qtauamuohcm, doublereal *rho_nu__, doublereal *p_nu__, 
                                                          doublereal *s_nu__)
{
  doublereal u_nutau0__, u_nuebar__, qtauatau, qtauttau, u_nuebar0__,
    h__, qtaua_nue__, qtaut_nue__, nb, qtaua_nuebar__, 
    qtaut_nuebar__, qtauatauohcm, qtauttauohcm, small, u_nue__, 
    u_nue0__, u_nuel__, u_numu__, u_numu0__;
  extern doublereal density_2stream__(doublereal *, doublereal *, 
                                      doublereal *);
  doublereal qtauamu, u_nutau__, qtautmu;

  /* ====================================================================== */
  /* outputs */
  /*     Passed to here */
  /* functions */
  /*     Locals */
  /*     Set some things */
  nb = *rhob / *mb;
  small = 1e-150;
  /*     Set some equivalences */
  h__ = *hcm;
  u_nue0__ = *unue0;
  u_nuebar0__ = *unuebar0;
  u_numu0__ = *unumu0;
  u_nutau0__ = *unumu0;
  qtauttauohcm = *qtautmuohcm;
  qtauatauohcm = *qtauamuohcm;
  /*     Set 2-stream approximation solutions */
  /* ccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*     Neutrino emission energy-volume-rate */
  qtaut_nue__ = *qtautnueohcm * *hcm;
  qtaua_nue__ = *qtauanueohcm * *hcm;
  qtaut_nuebar__ = *qtautnuebarohcm * *hcm;
  qtaua_nuebar__ = *qtauanuebarohcm * *hcm;
  qtautmu = *qtautmuohcm * *hcm;
  qtauamu = *qtauamuohcm * *hcm;
  qtauttau = qtauttauohcm * *hcm;
  qtauatau = qtauatauohcm * *hcm;
  /* ccccccccccccccccccccccccccccc */

  /*     Compute final neutrino energy and number densities */

  /* ccccccccccccccccccccccccccccc */
  /*     energy densities */
  u_nue__ = density_2stream__(&u_nue0__, &qtaut_nue__, &qtaua_nue__);
  u_nuebar__ = density_2stream__(&u_nuebar0__, &qtaut_nuebar__, &
                                 qtaua_nuebar__);
  u_nuel__ = u_nue__ + u_nuebar__;
  /* diag only */
  u_numu__ = density_2stream__(&u_numu0__, &qtautmu, &qtauamu);
  u_nutau__ = density_2stream__(&u_nutau0__, &qtauttau, &qtauatau);
  /*     total internal energy density */
  *rho_nu__ = u_nuel__ + u_numu__ + u_nutau__;
  /*     total pressure */
  *p_nu__ = *rho_nu__ / 3.;
  /*     entropy density in 1/cc */
  *s_nu__ = *rho_nu__ * 1.3333333333333333 / (*kbtk + small);
  /*     [k_b] form is s_nu/(rho/m_b) is "per baryon" entropy for rho and S given in cgs. */
  /* now in 1/cc */
  return 0;
} /* computefinal_justdensities_fromhcm__ */

/*     Compute Ynuthermal0.  No optical depth corrections. */
/* ====================================================================== */
/*     ENSURE this list is the same for how used in tau_calc.f */
/* Subroutine */ int compute_ynuthermal0_fromhcm__(doublereal *clight, 
                                                   doublereal *mb, doublereal *rhob, doublereal *nnueth0, doublereal *
                                                   nnuebarth0, doublereal *ynuthermal0)
{
  doublereal nb, n_nuebarth0__, small, n_nueth0__;

  /* ====================================================================== */
  /* outputs */
  /*     Passed to here */
  /*     Locals */
  /* outputs */
  /*     Set some things */
  nb = *rhob / *mb;
  small = 1e-150;
  n_nueth0__ = *nnueth0;
  n_nuebarth0__ = *nnuebarth0;
  /*     Thermalized faction Y_\nu */
  *ynuthermal0 = (n_nueth0__ - n_nuebarth0__) / nb;
  return 0;
} /* compute_ynuthermal0_fromhcm__ */

/* ===================================================================== */
doublereal rate_2stream__(doublereal *clight, doublereal *h__, doublereal *
                          density, doublereal *tautot, doublereal *tauabs)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  doublereal prefactor, small, osqrt3;

  /*     2-stream approximation for volume rate */
  /* ===================================================================== */
  /* ===================================================================== */
  /* ===================================================================== */
  small = 1e-150;
  /*      SMALL = 0.0d0 */
  /*     /H means rates are per unit volume */
  prefactor = *clight / 4.f / .75f / (*h__ + small);
  osqrt3 = .5773502691896f;
  ret_val = prefactor * *density / (*tautot * .5 + osqrt3 + 1. / (*tauabs * 
                                                                  3. + small) + small);
  return ret_val;
} /* rate_2stream__ */

/* ===================================================================== */
doublereal density_2stream__(doublereal *density, doublereal *tautot, 
                             doublereal *tauabs)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  doublereal twosqrt3, small;

  /*     2-stream approximation for density */
  /*     As tau->0, final density->0 */
  /*     As tau->\infty, final density -> density */
  /* ===================================================================== */
  /* ===================================================================== */
  /* ===================================================================== */
  small = 1e-150;
  /*      SMALL = 0.0d0 */
  /*      twosqrt3 = 2.0*sqrt(3.0) */
  twosqrt3 = 3.4641016151377545f;
  ret_val = *density * (*tautot * 3.f + twosqrt3) / (*tautot * 3.f + 
                                                     twosqrt3 + 2.f / (*tauabs + small) + small);
  return ret_val;
} /* density_2stream__ */

/* ===================================================================== */
doublereal tdifffromlambda_(doublereal *clight, doublereal *h__, doublereal *
                            lambda)
{
  /* System generated locals */
  doublereal ret_val, d__1, d__2;

  /* Local variables */
  doublereal small;

  /*     diffusion time limited by speed of light */
  /* ===================================================================== */
  /* ===================================================================== */
  /* ===================================================================== */
  small = 1e-150;
  /*      SMALL = 0.0d0 */
  /* Computing MAX */
  d__1 = *h__ * 3.f * *h__ / (small + *clight * 1.f * *lambda), d__2 = *h__ 
    / *clight;
  ret_val = max(d__1,d__2);
  return ret_val;
} /* tdifffromlambda_ */

/* ===================================================================== */
doublereal tthermalfromlambda_(doublereal *clight, doublereal *h__, 
                               doublereal *lambda)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  doublereal small;

  /*     thermalization time limited by speed of light */
  /* ===================================================================== */
  /* ===================================================================== */
  /* ===================================================================== */
  small = 1e-150;
  /*      SMALL = 0.0d0 */
  /*     Timescale to undergo a collision that thermalizes the particle.  Assume input lambda only thermalizing 
         (i.e. inelastic) opacities. */
  /*     H doesn't enter */
  ret_val = *lambda / *clight;
  /*     Assume it takes at least 3 inelastic events -- handled by exponential solution so that complete thermal
         ization takes about 3X the above time! */
  /*      tthermalfromlambda=3.0*lambda/clight */
  return ret_val;
} /* tthermalfromlambda_ */

