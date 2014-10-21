// actual EOS physics stuff
// stuff that should never directly access global EOS arrays
// all uses of lineartablelimits[] should be on non-temperature-like quantities (i.e. only on rho,ye,ynu,h)


static int needspostoffset(int whichd, int whichtablesubtype, int coli, int *whichdfake);
static int offsetquant2(FTYPE signature, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out);
static int offsetquant2_general(int whichd, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out);
static int offsetquant2_general_inverse(int whichd, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out);

static void constrain_H(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
static void constrain_TDYNORYE_YNU(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

static void get_lambdatot(FTYPE *EOSextra, FTYPE rho0, FTYPE u, FTYPE *lambdatotptr);
static FTYPE pressure_dp_rho0_u_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u, FTYPE *dp);
static FTYPE pressure_dp_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0, FTYPE *dp);


static int get_dquant2(int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *dquant2, FTYPE *unu, FTYPE *pnu, FTYPE *chinu, FTYPE *snu, FTYPE *quant2nu);
static int usereduced_eos(int domod, int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *final);


static FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
static int dfun2fun_inputdquant2_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE dquant2, FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, int badlookup, FTYPE *dfinal, FTYPE *final);


static int prepare_fudgefrac_kazfull(
                                     int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                     FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, // inputs
                                     FTYPE *utot, FTYPE *ugas, // outputs
                                     FTYPE *ptot, FTYPE *pgas, // outputs
                                     FTYPE *chitot, FTYPE *chigas // outputs
                                     );
static int prepare_fudgefrac_input_dpofwhichd_kazfull(
                                                      int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                                      FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, // inputs
                                                      FTYPE dpfun, // inputs
                                                      FTYPE *utot, FTYPE *ugas, // outputs
                                                      FTYPE *ptot, FTYPE *pgas, // outputs
                                                      FTYPE *chitot, FTYPE *chigas // outputs
                                                      );

static int preparepart2_fudgefrac_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, // inputs
                                          FTYPE quant2nu, FTYPE quant2, // inputs
                                          FTYPE *fakeneutrino, FTYPE *faketotal // outputs
                                          );
static int finish_fudgefrac_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                    FTYPE utot, FTYPE ugas, FTYPE unu, // inputs
                                    FTYPE ptot, FTYPE pgas, FTYPE pnu, // inputs
                                    FTYPE chitot, FTYPE chigas, FTYPE chinu, // inputs
                                    FTYPE fakeneutrino, FTYPE faketotal, // inputs
                                    int badlookup, FTYPE nonneutrino, // inputs
                                    FTYPE *final // outputs
                                    );
static FTYPE fudgefracsingle_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2);
static FTYPE compute_tabulated_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE rho0, FTYPE u);



static void reduce_extras(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra,FTYPE rho0, FTYPE du, FTYPE *extras, int *badlookups);

static void compute_allextras_du_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras);
static void compute_notallextras1_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras);
static void compute_notallextras1_du_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras);



static void compute_IJKglobal(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS]);
static void compute_Hglobal(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
static void compute_upsnu_global(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);



static int iterateupsnu(FTYPE *processed,FTYPE *pr,FTYPE *EOSextra);
static int iterateynu0(FTYPE *processed,FTYPE *pr,FTYPE *EOSextra);
static int compute_and_iterate_upsnu(int whichd, FTYPE *EOSextra, FTYPE *pr);
static int compute_and_iterate_ynu0_upsnu(int whichd, FTYPE *EOSextra, FTYPE *pr);
static int get_extrasprocessed_manyiterations(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);



// report whether need to process with offset and return the effective associated "whichd"
static int needspostoffset(int whichd, int whichtablesubtype, int coli, int *whichdfake)
{

  if(whichtablesubtype==SUBTYPEDEGEN){
    // if here, then called for degen value directly so need to process it
    // this is unlike internally used degen values that are all table-based with offset and so consistently used
    *whichdfake=whichd;
    // call coli are same type of quantity
    return(1);
  }
  else if(whichtablesubtype==SUBTYPESTANDARD){
    return(0);
  }
  else if(whichtablesubtype==SUBTYPEGUESS){
    if(coli==0){ *whichdfake=UTOTDIFF; return(1); }
    return(0);
  }
  else if(whichtablesubtype==SUBTYPEDISS){
    if(coli==0){ *whichdfake=UTOTDIFF; return(1); }
    return(0);
  }
  else if(whichtablesubtype==SUBTYPEDP){
    // GODMARK: derivatives!
  }
  else if(whichtablesubtype==SUBTYPESDEN){
    if(coli==0){ *whichdfake=STOTDIFF; return(1); }
    // GODMARK: derivatives!
    return(0);
  }
  else if(whichtablesubtype==SUBTYPESSPEC){
    // GODMARK: derivatives!
    if(coli==0){ *whichdfake=SSPECTOTDIFF; return(1); }
    return(0);
  }
  else if(whichtablesubtype==SUBTYPEPOFCHI){
    // GODMARK: derivatives!
    if(coli==0){ return(0);}
    return(0);
  }
  else if(whichtablesubtype==SUBTYPETEMP){
    return(0);
  }
  else if(whichtablesubtype==SUBTYPEEXTRA){
    // no neutrino quantities need correcting since correction is only on nuclear term
    return(0);
  }
  else{
    dualfprintf(fail_file,"No such whichtablesubtype=%d setup in needspostoffset()\n",whichtablesubtype);
    myexit(34968346);
  }
  
  return(0);

}



// general offset for energy per baryon to be applied always before accessing table with nuclear EOS
// quant2out is the \delta u that is tabulated in the EOS table itself.
// inputted quant2in is some true internal energy (maybe no neutrinos)
// So we add the offset before looking up the results from the table that has the offset
// signature : +1 for mapping HARM values to TABLE values
// signature : -1 for mapping TABLE values to HARM values

// only valid if rest-mass density and other energy densities are in same units w.r.t. speed of light
// so requires rho0unittype==1 in init.grb.c as much of code requires!
// Also, note that quant1 is in g/cc if rho0unittype==1
static int offsetquant2(FTYPE signature, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out)
{


  if(whichd==UTOTDIFF || whichd==CHIDIFF){
    // nuclear offset
    // assumes lsoffset in helm/jon_lsbox.f is offsetting only energy/baryon = u/\rho_0
    *quant2out = quant2in + signature*TRUENUCLEAROFFSET[primarytable]*quant1;
  }
  else if(whichd==STOTDIFF){// Sden
    *quant2out = quant2in + signature*TRUEENTROPYNUCLEAROFFSET[primarytable]*quant1;

    // DEBUG:
    //    dualfprintf(fail_file,"quant2in=%21.15g sig=%21.15g true=%21.15g quant1=%21.15g quant2out=%21.15g\n",quant2in,signature,TRUEENTROPYNUCLEAROFFSET[primarytable],quant1,*quant2out);



  }
  else if(whichd==SSPECTOTDIFF){// Ss=SSPEC
    *quant2out = quant2in + signature*TRUEENTROPYNUCLEAROFFSET[primarytable];
  }


  return(0);

}



// currently unused extra degen offset.  Was used when degens were negative due to small table size
static int offsetquant2_general(int whichd, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out)
{

  *quant2out = quant2in + DEGENNUCLEAROFFSET[primarytable]*quant1;

  return(0);

}

// currently unused extra degen offset.  Was used when degens were negative due to small table size
static int offsetquant2_general_inverse(int whichd, FTYPE quant1, FTYPE quant2in, FTYPE *quant2out)
{

  *quant2out = quant2in - DEGENNUCLEAROFFSET[primarytable]*quant1;

  return(0);

}




// get dquant2 = quant2 - quant2nu, the difference between full temperature-like quantity and neutrino version
// dquant2 used directly for lookup
static int get_dquant2(int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *dquant2, FTYPE *unu, FTYPE *pnu, FTYPE *chinu, FTYPE *snu, FTYPE *quant2nu)
{

  if(whichdatatype[primarytable]==4){

    *unu = EOSextra[UNUGLOBAL];
    *pnu = EOSextra[PNUGLOBAL];
    *chinu = (*unu + *pnu);
    *snu = EOSextra[SNUGLOBAL];

    if(whichd==UTOTDIFF){
      *quant2nu = *unu;
    }
    else if(whichd==PTOTDIFF){
      *quant2nu = *pnu;
    }
    else if(whichd==CHIDIFF){
      *quant2nu = *chinu;
    }
    else if(whichd==STOTDIFF){
      *quant2nu = *snu;
    }

    // get dquant2
    *dquant2  = quant2 - *quant2nu;

  }
  else{
    // then neutrino inside table already
    *unu= *pnu= *chinu= *snu=0.0;
    *dquant2 = quant2;
    *quant2nu = *chinu = *pnu = *unu = 0.0;
  }



  return(0);

}





// tabulated dfun[du] and need to get u=u+unu and then fun= dfun+fun_nu
static FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn)
{
  FTYPE final,dfinal;
  FTYPE unu,pnu,chinu,snu;
  FTYPE dquant2,quant2nu;
  int badlookup;

  ///////////////////
  //
  // get dquant2 for lookup and get associated \nu quantities for the final answer
  //
  ///////////////////
  get_dquant2(whichd, EOSextra, quant1, quant2, &dquant2, &unu, &pnu, &chinu, &snu, &quant2nu);


  badlookup=getsingle_eos_fromtable(whichfun,whichd,EOSextra,quant1,dquant2,&dfinal); // input quant1,dquant2 and get dfinal

  // get final and possibly different dfinal results
  dfun2fun_inputdquant2_kazfull(whichfun, whichd, EOSextra, quant1, quant2, dquant2, unu, pnu, chinu, snu, badlookup, &dfinal, &final);

  //////////////
  //
  // sometimes need pure "gas" (non-neutrino) part of answer
  //
  //////////////
  *dfinalreturn=dfinal;
 

  // DEBUG:
  //  if(whichfun==PofRHOU){
  //    dualfprintf(fail_file,"whichd=%d EOSe[0]=%21.15g ig=%d q1=%21.15g q2=%21.15g dquant2=%21.15g unu=%21.15g pnu=%21.15g chinu=%21.15g snu=%21.15g quant2unu=%21.15g\n",whichd,EOSextra[TDYNORYEGLOBAL],(int)EOSextra[IGLOBAL],quant1,quant2,dquant2,unu,pnu,chinu,snu,quant2nu);
  //    dualfprintf(fail_file,"dfinal=%21.15g final=%21.15g : %d\n",dfinal,final,badlookup);
  //  }

 
  //////////////
  //
  // return final full answer
  //
  //////////////
  return(final);

}







// tabulated dfun[du] and need to get u=u+unu and then fun= dfun+fun_nu
// no lookup here!  Just process to get "final" full EOS value
static int dfun2fun_inputdquant2_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE dquant2, FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, int badlookup, FTYPE *dfinal, FTYPE *final)
{


  ///////////////////
  //
  // do lookup
  //
  ///////////////////
  if(badlookup){
    usereduced_eos(REDUCEUSEOFFSET,whichfun,whichd,EOSextra,quant1,quant2,final); // final is pointer for output
    *dfinal=*final; // correct dfinal
  }
  else{
    // sum up to get final answer for function value
    if(whichfun==PofRHOU)          *final = *dfinal + pnu;
    else if(whichfun==UofRHOP)     *final = *dfinal + unu;
    else if(whichfun==UofRHOS)     *final = *dfinal + unu;
    else if(whichfun==SofRHOU)     *final = *dfinal + snu;
    else if(whichfun==SSofRHOCHI)  *final = *dfinal + snu/quant1; // specific entropy is snu/rho0 in code units
    else if(whichfun==PofRHOCHI)   *final = *dfinal + pnu;
  }


  return(0);
}



// use reduced ideal-like EOS since (probably) bad table lookup
static int usereduced_eos(int domod, int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *final)
{
  FTYPE quant2mod;


  if(domod==REDUCEUSEOFFSET){
    if(whichd==UTOTDIFF || whichd==CHIDIFF){
      quant2mod = quant2 + FAKE2IDEALNUCLEAROFFSET[primarytable]*quant1; // set nuclear per baryon offset so can smoothly connect to ideal gas EOS
    }
    else if(whichd==STOTDIFF){
      // SUPERGODMARK: Need to offset by u/(kb*t)
      //      quant2mod = quant2 + FAKE2IDEALNUCLEAROFFSET[primarytable]*quant1; // set nuclear per baryon offset so can smoothly connect to ideal gas EOS     
    }
  }
  else{
    // no mod
    quant2mod = quant2;
  }


  // otherwise use TM EOS
#if(REDUCE2WHICHEOS==MIGNONE)
  if(whichfun==PofRHOU)        *final = pressure_rho0_u_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==UofRHOP)   *final = u_rho0_p_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==UofRHOT)   *final = u_rho0_T_mignone(EOSextra,quant1, quant2mod);// SUPERGODMARK: NOT setup yet
  else if(whichfun==UofRHOS)   *final = compute_u_from_entropy_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==SofRHOU)   *final = compute_entropy_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==SSofRHOCHI)   *final = compute_specificentropy_wmrho0_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==PofRHOCHI) *final = pressure_wmrho0_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==CS2ofRHOU)           *final = cs2_compute_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DPDUofRHOU)     *final = dpdu_rho0_u_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DPDRHOofRHOU)   *final = dpdrho0_rho0_u_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==IDRHO0DP)       *final = compute_idrho0dp_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==IDCHIDP)        *final = compute_idwmrho0dp_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DSDRHOofRHOU)   *final = compute_dSdrho_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DSDUofRHOU)     *final = compute_dSdu_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DSSDRHOofRHOCHI) *final = compute_dspecificSdrho_wmrho0_mignone(EOSextra,quant1, quant2mod);
  else if(whichfun==DSSDCHIofRHOCHI) *final = compute_dspecificSdwmrho0_wmrho0_mignone(EOSextra,quant1, quant2mod);
#elif(REDUCE2WHICHEOS==IDEALGAS)
  // use ideal EOS
  if(whichfun==PofRHOU)        *final = pressure_rho0_u_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==UofRHOP)   *final = u_rho0_p_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==UofRHOT)   *final = u_rho0_T_idealgas(EOSextra,quant1, quant2mod);// SUPERGODMARK: NOT setup yet
  else if(whichfun==UofRHOS)   *final = compute_u_from_entropy_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==SofRHOU)   *final = compute_entropy_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==SSofRHOCHI)   *final = compute_specificentropy_wmrho0_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==PofRHOCHI) *final = pressure_wmrho0_idealgas(EOSextra,quant1, quant2mod);    
  else if(whichfun==CS2ofRHOU)           *final = cs2_compute_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DPDUofRHOU)     *final = dpdu_rho0_u_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DPDRHOofRHOU)   *final = dpdrho0_rho0_u_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==IDRHO0DP)       *final = compute_idrho0dp_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==IDCHIDP)        *final = compute_idwmrho0dp_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DSDRHOofRHOU)   *final = compute_dSdrho_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DSDUofRHOU)     *final = compute_dSdu_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DSSDRHOofRHOCHI) *final = compute_dspecificSdrho_wmrho0_idealgas(EOSextra,quant1, quant2mod);
  else if(whichfun==DSSDCHIofRHOCHI) *final = compute_dspecificSdwmrho0_wmrho0_idealgas(EOSextra,quant1, quant2mod);
#endif

  return(0);
}





// setup for fudgefrac()
// Introduces another lookup, so kinda expensive for something like the derivative that probably does not often have dominant neutrino contribution, so use prepare_fudgefrac_input_dpofwhichd_kazfull() if can.
static int prepare_fudgefrac_kazfull(
                                     int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                     FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, // inputs
                                     FTYPE *utot, FTYPE *ugas, // outputs
                                     FTYPE *ptot, FTYPE *pgas, // outputs
                                     FTYPE *chitot, FTYPE *chigas // outputs
                                     )
{
  FTYPE dpfun;


  // first get total pressure in order to frac-fudge the answer
  if(whichd==UTOTDIFF){
    // pressure_dp_rho0_u() takes in P(rho0,u) not P(rho0,du), so feed in quant1,quant2
    *ptot = pressure_dp_rho0_u_kazfull(EOSextra, quant1, quant2, pgas); // need pgas to form chi since otherwise don't have it
    dpfun=*pgas;
  }
  else if(0&&whichd==PTOTDIFF){
    // Not setup.  Not needed yet.  Would require (utot,ugas)[rho0,p]
  }
  else if(whichd==CHIDIFF){
    // pressure_dp_rho0_u() takes in P(rho0,u) not P(rho0,du), so feed in quant1,quant2
    *ptot =pressure_dp_wmrho0_kazfull(EOSextra, quant1, quant2, pgas); // don't need separate pgas
    dpfun=*pgas;
  }
  else if(0&&whichd==STOTDIFF){
    // don't need fudgefrac(STOTDIFF) since only used for complicated functions like derivatives that have mixed derivatives between gas and neutrino terms
    // pressure_dp_rho0_s() takes in P(rho0,Sden) not P(rho0,dSden), so feed in quant1,quant2
    //    utot =u_du_rho0_s_kazfull(EOSextra, quant1, quant2, &ugas); // need ugas to form chi since otherwise don't have it
    //    ptot =pressure_dp_rho0_s_kazfull(EOSextra, quant1, quant2, &pgas); // need pgas to form chi since otherwise don't have it
    //    chigas = ugas+pgas;
  }
  else{
    dualfprintf(fail_file,"prepare_fudgefrac_kazfull not setup for whichd=%d\n",whichd);
    myexit(19672606);
  }


  // get other dependent results
  prepare_fudgefrac_input_dpofwhichd_kazfull(whichd, EOSextra, quant1, quant2,
                                             unu, pnu, chinu, snu,
                                             dpfun,
                                             utot, ugas,
                                             ptot, pgas,
                                             chitot, chigas);



  return(0);

}


// setup for fudgefrac()
// input dpfun instead of extra lookup
static int prepare_fudgefrac_input_dpofwhichd_kazfull(
                                                      int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                                      FTYPE unu, FTYPE pnu, FTYPE chinu, FTYPE snu, // inputs
                                                      FTYPE dpfun, // inputs
                                                      FTYPE *utot, FTYPE *ugas, // outputs
                                                      FTYPE *ptot, FTYPE *pgas, // outputs
                                                      FTYPE *chitot, FTYPE *chigas // outputs
                                                      )
{

  // first get total pressure in order to frac-fudge the answer
  if(whichd==UTOTDIFF){
    // pressure_dp_rho0_u() takes in P(rho0,u) not P(rho0,du), so feed in quant1,quant2
    *pgas = dpfun;
    *ptot = *pgas + pnu;

    *utot = quant2;
    *ugas = *utot - unu;
    *chigas = *ugas + *pgas;
    *chitot = *chigas + chinu;
  }
  else if(0&&whichd==PTOTDIFF){
    // Not setup.  Not needed yet.  Would require (utot,ugas)[rho0,p]
  }
  else if(whichd==CHIDIFF){
    // pressure_dp_rho0_u() takes in P(rho0,u) not P(rho0,du), so feed in quant1,quant2
    *pgas = dpfun;
    *ptot = *pgas + pnu;

    *chitot = quant2;
    *chigas = *chitot - chinu;
    *utot = *chitot - *ptot;
    *ugas = *utot - unu;
  }
  else if(0&&whichd==STOTDIFF){
    // don't need fudgefrac(STOTDIFF) since only used for complicated functions like derivatives that have mixed derivatives between gas and neutrino terms
    // pressure_dp_rho0_s() takes in P(rho0,Sden) not P(rho0,dSden), so feed in quant1,quant2
    //    utot =u_du_rho0_s_kazfull(EOSextra, quant1, quant2, &ugas); // need ugas to form chi since otherwise don't have it
    //    ptot =pressure_dp_rho0_s_kazfull(EOSextra, quant1, quant2, &pgas); // need pgas to form chi since otherwise don't have it
    //    chigas = ugas+pgas;
  }
  else{
    dualfprintf(fail_file,"prepare_fudgefrac_input_dpofwhichd_kazfull not setup for whichd=%d\n",whichd);
    myexit(19672606);
  }





  return(0);

}





static int preparepart2_fudgefrac_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, // inputs
                                          FTYPE quant2nu, FTYPE quant2, // inputs
                                          FTYPE *fakeneutrino, FTYPE *faketotal // outputs
                                          )
{


  ////////////////////////////////////////
  //
  // determine fake-neutrino part
  // GODMARK: Assume at least that if neutrino pressure is dominant then treat like mixture of gas+neutrinos
  //
  ////////////////////////////////////////
  usereduced_eos(REDUCENOOFFSET,whichfun,whichd, EOSextra, quant1, quant2nu, fakeneutrino);

  // reduce to using totals if not within table
  usereduced_eos(REDUCEUSEOFFSET,whichfun,whichd, EOSextra, quant1, quant2, faketotal);

  return(0);
}





// finish-up for fudgefrac()
static int finish_fudgefrac_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, // inputs
                                    FTYPE utot, FTYPE ugas, FTYPE unu, // inputs
                                    FTYPE ptot, FTYPE pgas, FTYPE pnu, // inputs
                                    FTYPE chitot, FTYPE chigas, FTYPE chinu, // inputs
                                    FTYPE fakeneutrino, FTYPE faketotal, // inputs
                                    int badlookup, FTYPE nonneutrino, // inputs
                                    FTYPE *final // outputs
                                    )
{
  FTYPE frac;



  if(badlookup){
    *final=faketotal;
  }
  else{

    // then estimate importance of neutrinos and its control over sound speed
    if(whichfun==CS2ofRHOU){
      // for sound speed more accurate to do:
      // c_s^2[total] = c_s^2[gas]*(hgas/htotal) + c_s^2[neutrino]*(hneutrino/htotal)
      // hgas = (rho_0 + ugas + pgas)/rho_0
      // hneutrino = (rho_0 + unu + pnu)/rho_0
      *final = (nonneutrino*(quant1+chigas) + fakeneutrino*(quant1+chinu))/(quant1+chitot);

      //      if(!isfinite(final)){
      // dualfprintf(fail_file,"IF: final=%21.15g chitot=%21.15g\n",final,chitot);
      // dualfprintf(fail_file,"%21.15g%21.15g%21.15g%21.15g%21.15g%21.15g%21.15g%21.15g\n",nonneutrino,quant1,chigas,fakeneutrino,quant1,chinu,quant1,chitot);
      //      }
      
      
    }
    else{
      //    frac = (pnu/(ptot-pnu));
      frac = fabs(pnu)/(fabs(ptot)+SMALL); // -> 1 if ptot=pnu.  Suppose photons+electrons+neutrinos equally dominate.  Then cs2 still correct.
      // choose cs2total if neutrino-dominated, otherwise choose non-neutrino cs2
      *final = faketotal*frac + nonneutrino*(1.0-frac);
    }
  }// end else if good lookup

  ///////////////
  //
  // return that was success
  //
  ///////////////
  return(0);


}





// general function to interpolate between non-neutrino and neutrino values
// used for those quantities not yet setup for exactly correct answer (more specifically, fudgefrac() used for derivatives that would otherwise involved mixed derivatives between gas and neutrino terms that are not easily computed/stored without many extra things to store
// never used on pressure itself, but pressure looked-up in prepare_fudgefrac_kazfull() so expensive in that sense
static FTYPE fudgefracsingle_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2)
{
  FTYPE ptot, pgas, pnu;
  FTYPE utot, ugas, unu;
  FTYPE chitot, chigas, chinu;
  FTYPE snu;
  FTYPE dquant2,quant2nu,quant2mod;
  FTYPE fakeneutrino,faketotal;

  FTYPE nonneutrino;
  int badlookup;

  FTYPE final;



  ///////////////////
  //
  // get dquant2 for lookup and get associated \nu quantities for the final answer
  //
  ///////////////////
  get_dquant2(whichd, EOSextra, quant1, quant2, &dquant2, &unu, &pnu, &chinu, &snu, &quant2nu);


  /////////////////////////////////
  //
  // now get non-neutrino value
  // nonneutrino will then be only "gas" part without neutrinos
  //
  /////////////////////////////////
  badlookup=getsingle_eos_fromtable(whichfun,whichd,EOSextra,quant1,dquant2,
                                    &nonneutrino);


  ////////////////////////
  //
  // setup fudgefrac()
  //
  ////////////////////////
  prepare_fudgefrac_kazfull(whichd, EOSextra, quant1, quant2,
                            unu, pnu, chinu, snu,
                            &utot, &ugas,
                            &ptot, &pgas,
                            &chitot, &chigas
                            );

  ////////////////////////
  //
  // get estimated full neutrino and fake total values
  //
  ////////////////////////
  preparepart2_fudgefrac_kazfull(whichfun, whichd, EOSextra, quant1,
                                 quant2nu, quant2,
                                 &fakeneutrino, &faketotal
                                 );

  ////////////////////////
  //
  // finish fudgefrac()
  //
  ////////////////////////
  finish_fudgefrac_kazfull(whichfun, whichd, EOSextra, quant1, quant2,
                           utot, ugas, unu,
                           ptot, pgas, pnu,
                           chitot, chitot, chinu,
                           fakeneutrino, faketotal,
                           badlookup, nonneutrino,
                           &final);

  ////////////////////////
  //
  // return single value
  //
  ////////////////////////
  return(final);

}




// returns all 3 quantities for inversion based upon p(\rho_0,\chi) or sspec(\rho_0,\chi)
// note that no optimization for 5D method based upon p(\rho_0,u), but could set that up.  Would have to add P to dP table for simplicity.
// Note that quant2 must coincide with quantity obtaining.  So far quant2=wmrho0=\chi for whichd=CHIDIFF quantities.
void getall_forinversion_kazfull(int eomtype, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2, FTYPE *fun, FTYPE *dfunofrho, FTYPE *dfunofu)
{
  FTYPE utot, ugas, unu;
  FTYPE ptot, pgas, pnu;
  FTYPE chitot, chigas, chinu;
  FTYPE snu;
  FTYPE dquant2,quant2nu;
  FTYPE fakeneutrino,faketotal;

  FTYPE deltafun; // pgas really

  FTYPE nonneutrinos[MAXEOSPIPELINE];
  FTYPE finals[MAXEOSPIPELINE];
  int badlookups[MAXEOSPIPELINE];
  int iffun[MAXEOSPIPELINE];
  int whichtablesubtype;
  int numcols,coli;
  int firsteos;

  int returnfun,returndfunofrho,returndfunofu;
  int failreturn;
  int whichfun;
  int simplefudge;



  ///////////////
  //
  // Quickly return if doing cold or force-free
  //
  ///////////////
  if(eomtype==EOMCOLDGRMHD || eomtype==EOMFFDE || eomtype==EOMFFDE2){
    *fun=*dfunofrho=*dfunofu=0.0;
    return; // done!
  }


  ///////////////
  //
  // get number of columns possible for this whichtablesubtype
  //
  // returnfun, returndfunofrho, and returndfunofu are in "whichfun" format, not "coli" format
  //
  ///////////////
  if(eomtype==EOMGRMHD){
    // whichd==CHIDIFF
    whichtablesubtype=SUBTYPEPOFCHI;
    returnfun = PofRHOCHI;
    returndfunofrho = IDRHO0DP;
    returndfunofu = IDCHIDP;
  }
  else if(eomtype==EOMENTROPYGRMHD){
    // whichd==CHIDIFF
    whichtablesubtype=SUBTYPESSPEC;
    returnfun = SSofRHOCHI;
    returndfunofrho = DSSDRHOofRHOCHI;
    returndfunofu = DSSDCHIofRHOCHI;
  }
  firsteos=firsteosintablesubtype[whichtablesubtype];
  numcols = numcolintablesubtype[whichtablesubtype]; // ok use of numcolintablesubtype

  // setup iffun
  for(coli=0;coli<numcols;coli++) iffun[coli]=1; // want all!


  ///////////////////
  //
  // get dquant2 for lookup and get associated \nu quantities for the final answer
  //
  ///////////////////
  get_dquant2(whichd, EOSextra, quant1, quant2, &dquant2, &unu, &pnu, &chinu, &snu, &quant2nu);






  /////////////////////////////////
  //
  // now get non-neutrino value
  // nonneutrino will then be only "gas" part without neutrinos
  //
  /////////////////////////////////
  failreturn=get_eos_fromtable(whichtablesubtype, iffun, whichd,EOSextra,quant1,dquant2,nonneutrinos,badlookups);

  

  if(whichtablesubtype==SUBTYPEPOFCHI){
    // get pgas = p - pnu
    deltafun = nonneutrinos[returnfun-firsteos];
    simplefudge=1;
  }
  else{
    // don't have pressure to use to interpolate result, so have to do full fudge
    // only have to do for coli==0 and have for rest of coli>0
    simplefudge=0;
  }

  ////////////////////////
  //
  // get final results for function and its derivatives
  //
  ////////////////////////

  // used to indicate if already hit derivative "coli" so must already have computed deltafun
  int hitdercoli=0;

  for(coli=0;coli<numcols;coli++){ // don't need to deal with iffun since want all

    whichfun = coli+firsteos;

    if(whichfun==returndfunofrho || whichfun==returndfunofu){

      ////////////////////////
      //
      // setup fudgefrac() for dpofchidrho0 and dpofchidchi
      //
      ////////////////////////

      if(simplefudge==0 && hitdercoli==0){
        // then need to look up pressure still
        prepare_fudgefrac_kazfull(whichd, EOSextra, quant1, quant2,
                                  unu, pnu, chinu, snu,
                                  &utot, &ugas,
                                  &ptot, &pgas,
                                  &chitot, &chigas
                                  );

        deltafun=pgas; // after coli==0,  have this for all other coli
        hitdercoli=1; // so never come back into this if-else segment
      }
      else{
        // special prepare_fudgefrac that inputs already-looked-up p-pnu
        prepare_fudgefrac_input_dpofwhichd_kazfull(whichd, EOSextra, quant1, quant2,
                                                   unu, pnu, chinu, snu,
                                                   deltafun,
                                                   &utot, &ugas,
                                                   &ptot, &pgas,
                                                   &chitot, &chigas
                                                   );
      }


      // get preparepart2_fudgefrac(), which is done per final derivative quantity
      preparepart2_fudgefrac_kazfull(whichfun, whichd, EOSextra, quant1,
                                     quant2nu, quant2,
                                     &fakeneutrino, &faketotal
                                     );

      
      if(failreturn==0){
        // finish fudgefrac(), which is done per final derivative quantity
        // this call takes care of if EOS lookup returned badlookups==1
        finish_fudgefrac_kazfull(whichfun, whichd, EOSextra, quant1, quant2,
                                 utot, ugas, unu,
                                 ptot, pgas, pnu,
                                 chitot, chitot, chinu,
                                 fakeneutrino, faketotal,
                                 badlookups[coli], nonneutrinos[coli],
                                 &finals[coli]);
      }
      else{
        // reduce to fake total if true failure in lookup
        finals[coli]=faketotal;
      }

      ////////
      //
      // assign to final values
      //
      ///////
      if(whichfun==returndfunofrho){
        *dfunofrho=finals[coli];
      }
      else if(whichfun==returndfunofu){
        *dfunofu=finals[coli];
      }
    }// end if a derivative
    else if(whichfun==returnfun){
      // then process non-derivative using finish_dfun2fun()

      // finish dffun2fun(), which is done for the non-derivative quantity
      FTYPE dfun=nonneutrinos[coli];
      // below can change dfun, but change not used
      // this call takes care of if EOS lookup returned badlookups==1
      dfun2fun_inputdquant2_kazfull(whichfun, whichd, EOSextra, quant1, quant2, dquant2, unu, pnu, chinu, snu, badlookups[coli], &dfun, &finals[coli]);

      // assign final avlues
      *fun=finals[coli];

    }
  }// end over coli


  ////////////////////////
  //
  // return whether successful
  //
  ////////////////////////
  if(failreturn){
    *fun=sqrt(-1); // force failure for now
  }


  // BEGIN DEBUG
  //  dualfprintf(fail_file,"whichd=%d q1=%21.15g q2=%21.15g fun=%21.15g dfunofrho=%21.15g dfunofu=%21.15g\n",whichd,quant1,quant2,*fun,*dfunofrho,*dfunofu);
  //  dualfprintf(fail_file,"whichtablesubtype=%d :: nonneutrinos0=%21.15g nonneutrinos0=%21.15g nonneutrinos0=%21.15g\n",whichtablesubtype,nonneutrinos[0],nonneutrinos[1],nonneutrinos[2]);
  //  dualfprintf(fail_file,"finals0=%21.15g finals0=%21.15g finals0=%21.15g\n",finals[0],finals[1],finals[2]);
  //  dualfprintf(fail_file,"%d %21.15g %21.15g %21.15g : %21.15g : %21.15g %21.15g %21.15g %21.15g %21.15g\n",whichd, EOSextra[0], quant1, quant2, dquant2, unu, pnu, chinu, snu, quant2nu);
  //  dualfprintf(fail_file,"fakes: %21.15g %21.15g\n",fakeneutrino, faketotal);
  // END DEBUG

  return;

}





// used for fully tabulated quantities that are functions of du/dp/dchi with whichdatatype==4
FTYPE compute_tabulated_kazfull(int whichfun, int whichd, FTYPE *EOSextra, FTYPE quant1, FTYPE quant2)
{
  FTYPE final;
  FTYPE dquant2;
  FTYPE unu, pnu, chinu, snu, quant2nu;

  ///////////////////
  //
  // get dquant2 for lookup and get associated \nu quantities for the final answer
  //
  ///////////////////
  get_dquant2(whichd, EOSextra, quant1, quant2, &dquant2, &unu, &pnu, &chinu, &snu, &quant2nu);

  
  // do lookup
  if(getsingle_eos_fromtable(whichfun,whichd,EOSextra,quant1,dquant2,&final)){ // uses dquant2
    if(whichfun==TEMPGEN) final=0.0;
    else if(whichfun==QDOTNU) final=0.0;
    else final=0.0; // use mignone?
  }  

  return(final);

}










// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE dp;
  return(dfun2fun_kazfull(PofRHOU, UTOTDIFF, EOSextra, rho0, u, &dp));
}

// p(rho0, u) (used internally)
static FTYPE pressure_dp_rho0_u_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u, FTYPE *dp)
{
  return(dfun2fun_kazfull(PofRHOU, UTOTDIFF, EOSextra, rho0, u, dp));
}



// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE p)
{
  FTYPE du;
  return(dfun2fun_kazfull(UofRHOP, PTOTDIFF, EOSextra, rho0, p, &du));
}

// u(rho0, T) (used for initial conditions)
// SUPERGODMARK: NOT setup yet
FTYPE u_rho0_T_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE T)
{
  FTYPE dT; // maybe not right.
  return(dfun2fun_kazfull(UofRHOT, TTOTDIFF, EOSextra, rho0, T, &dT));
}


// p(rho0, w-rho0 = u+p)
// Notice that using EOSextraglobal bypasses need to have other quantities as functions of wmrho0 unless want more direct (by iteration) result compared to what EOSextraglobal gives
FTYPE pressure_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dp;
  return(dfun2fun_kazfull(PofRHOCHI, CHIDIFF, EOSextra, rho0, wmrho0, &dp));
}

// used internally
static FTYPE pressure_dp_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0, FTYPE *dp)
{
  return(dfun2fun_kazfull(PofRHOCHI, CHIDIFF, EOSextra, rho0, wmrho0, dp));
}






// frac-fudged because dpdu is complicated with p=dp+pnu and u=du+unu
FTYPE dpdu_rho0_u_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(fudgefracsingle_kazfull(DPDRHOofRHOU,UTOTDIFF,EOSextra,rho0, u));
}


// for this could store dp_\nu/drho0, which for whichdatatype==4 is 0 given other independent variables are rho0,T,Y_e,Y_\nu,H
// for whichdatatype==4 check that dpnu/drho0 = 0 for whichdatatype==4 if ever want to use that method (not now)
// dp(rho0, u)/drho0
// frac-fudged since dp/drho0 would otherwise involve dpnu/drho0 that would require storing many neutrino-related derivatives
FTYPE dpdrho0_rho0_u_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(fudgefracsingle_kazfull(DPDUofRHOU,UTOTDIFF,EOSextra,rho0, u));
}


// sound speed squared (for vchar.c) -- important for treatment of shocks
// frac-fudged
FTYPE cs2_compute_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(fudgefracsingle_kazfull(CS2ofRHOU,UTOTDIFF,EOSextra,rho0, u));
}


// entropy as function of rho0 and internal energy (u)
// Sden(rho0,u)
// tabulated ds(du), so first compute du and then ds and then s
FTYPE compute_entropy_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  FTYPE ds;
  return(dfun2fun_kazfull(SofRHOU, UTOTDIFF, EOSextra, rho0, u, &ds));
}

// u(rho0,Sden)
// here S is entropy is entropy per unit volume
FTYPE compute_u_from_entropy_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy)
{
  FTYPE du;
  return(dfun2fun_kazfull(UofRHOS, STOTDIFF, EOSextra, rho0, entropy, &du));
}


// used for dudp_calc
// frac-fudged
FTYPE compute_dSdrho_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(fudgefracsingle_kazfull(DSDRHOofRHOU,UTOTDIFF,EOSextra,rho0, u));
}


// used for dudp_calc
// frac-fudged
FTYPE compute_dSdu_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(fudgefracsingle_kazfull(DSDUofRHOU,UTOTDIFF, EOSextra, rho0, u));
}

// exactly correct answer (not frac-fudged)
// entropy as function of rho0 and internal energy (u)
// S(rho0,\chi)
// tabulated ds(d\chi), so first compute d\chi and then ds and then s
FTYPE compute_specificentropy_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  FTYPE ds;
  return(dfun2fun_kazfull(SSofRHOCHI, CHIDIFF, EOSextra, rho0, wmrho0, &ds));
}

// used for utoprim_jon entropy inversion
// frac-fudged
FTYPE compute_dspecificSdrho_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(fudgefracsingle_kazfull(DSSDRHOofRHOCHI,CHIDIFF,EOSextra,rho0, wmrho0));
}


// used for utoprim_jon entropy inversion
// frac-fudged
FTYPE compute_dspecificSdwmrho0_wmrho0_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(fudgefracsingle_kazfull(DSSDCHIofRHOCHI,CHIDIFF, EOSextra, rho0, wmrho0));
}




// 1 / (drho0/dp) holding wmrho0 fixed
// frac-fudged
FTYPE compute_idrho0dp_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(fudgefracsingle_kazfull(IDRHO0DP,CHIDIFF, EOSextra, rho0, wmrho0));
}



// 1 / (d(u+p)/dp) holding rho0 fixed
// frac-fudged
FTYPE compute_idwmrho0dp_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0)
{
  return(fudgefracsingle_kazfull(IDRHO0DP,CHIDIFF, EOSextra, rho0, wmrho0));
}




// volume heating rate(rho0,u)
FTYPE compute_qdot_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  int whichfun;

  if(whichdatatype[primarytable]==1){
    whichfun=EXTRA1;
  }
  else if(whichdatatype[primarytable]==3){
    whichfun=EXTRA1;
  }
  else{
    dualfprintf(fail_file,"Shouldn't request whichfun0=%d if primarytable=%d\n",QDOTNU,primarytable);
  }

  return(compute_tabulated_kazfull(whichfun, UTOTDIFF, EOSextra, rho0, u));
}


// temperature(rho0,u)
FTYPE compute_temp_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u)
{
  return(compute_tabulated_kazfull(TEMPGEN, UTOTDIFF, EOSextra, rho0, u));
}


// temperature (direct lookup from differential quant2: dquant2)
FTYPE compute_temp_whichd_kazfull(int whichd, FTYPE *EOSextra, FTYPE rho0, FTYPE dquant2)
{
  FTYPE temp;


  if(getsingle_eos_fromtable(TEMPGEN,whichd,EOSextra,rho0,dquant2,&temp)){
    temp=0.0;
  }

  return(temp);



}



// not pipelined version
// uses global variable "numextras" that should be set before this function is called
// f(rho0,du)
// direct lookup
static void compute_allextras_du_kazfull_old(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras)
{
  int whichtablesubtype=SUBTYPEEXTRA;
  int whichd=UTOTDIFF;
  int whichfun;
  int coli;
  int numcols;


  *numextrasreturn = numcols = numcolintablesubtype[whichtablesubtype];// ok use of numcolintablesubtype
  
  
  if(justnum==0){
    // assume all tables have same number of extras or else this doesn't make sense in general
    for(coli=0;coli<numcols;coli++){
      whichfun=EXTRA1+coli;
      if(getsingle_eos_fromtable(whichfun,whichd,EOSextra,rho0,du,&(extras[coli]))){
        extras[coli]=0.0;
        *numextrasreturn=0; // tells calling function that no tabulated extras
        if(WHICHDATATYPEGENERAL==4){
          extras[whichcolinquantity[EXTRA19]]=extras[whichcolinquantity[EXTRA20]]=BIG;
        }
        else{
          dualfprintf(fail_file,"Not setup asdfasdfasd\n");
          myexit(2093851);
        }
      }
    }

    // set rest to 0 for diagnostics to be uniform
    // These won't be used during calculation
    for(coli=numcols;coli<MAXNUMEXTRAS;coli++){
      extras[coli] = 0.0;
    }
  }// end if justnum==0

}



// pipelined version
static void compute_allextras_du_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras)
{
  int badlookups[MAXEOSPIPELINE];
  int iffun[MAXEOSPIPELINE];
  int whichtablesubtype=SUBTYPEEXTRA;
  int whichd=UTOTDIFF;
  int numcols,coli;
  int failreturn;


  *numextrasreturn = numcols = numcolintablesubtype[whichtablesubtype];// ok use of numcolintablesubtype

  
  if(justnum==0){

    // setup iffun
    // want all extras
    for(coli=0;coli<numcols;coli++) iffun[coli]=1;

    // do pipelined lookup
    failreturn=get_eos_fromtable(whichtablesubtype,iffun,whichd,EOSextra,rho0,du,extras,badlookups);

    // deal with bad lookups for extras
    reduce_extras(whichtablesubtype,iffun,whichd,EOSextra,rho0,du,extras,badlookups);

    int numbad=0;
    for(coli=0;coli<numcols;coli++) if(iffun[coli] && badlookups[coli]) numbad++;
    if(numbad==numcols) *numextrasreturn=0; // indicate that no point in doing more complex computations based upon these results -- can reduce to some simpler estimates


    // set rest to 0 for diagnostics to be uniform
    // These won't be used during calculation
    for(coli=numcols;coli<MAXNUMEXTRAS;coli++){
      extras[coli] = 0.0;
    }

    // DEBUG:
    //    dualfprintf(fail_file,"compute_allextras_du_kazfull: rho0=%21.15g du=%21.15g\n",rho0,du);
    //    for(coli=0;coli<numcols;coli++) dualfprintf(fail_file,"extras[%d]=%21.15g bad=%d\n",coli,extras[coli],badlookups[coli]);


  }// end if justnum==0

}


// deal with bad lookups (i.e. reduce to ideal gas effectively for extras, but idealgaseos.c is not aware of meaning of extras and how they should reduce)
static void reduce_extras(int whichtablesubtype, int *iffun, int whichd, FTYPE *EOSextra,FTYPE rho0, FTYPE du, FTYPE *extras, int *badlookups)
{
  int numcols,coli;

  numcols = numcolintablesubtype[whichtablesubtype];// ok use of numcolintablesubtype


  for(coli=0;coli<numcols;coli++){
    if(iffun[coli] && badlookups[coli]){ // only check if iffun==1

      // make zero
      extras[coli]=0.0;

      // except things that should be made large:
      if(WHICHDATATYPEGENERAL==4){
        // then not in table, so revert to estimated non-tabulated values (i.e. optically thin)
        // extra quantities were assumed set to 0 if not in table, so only need to change non-zero things
        // optical depths and densities are complicated and if out of table then they aren't determined easily
        extras[whichcolinquantity[EXTRA19]]=extras[whichcolinquantity[EXTRA20]]=BIG;
      }
      else{
        dualfprintf(fail_file,"Not setup for other data types in reduce_extras()\n");
        myexit(2093852);
      }

    }// end if badlookup
  }// end over coli

}



// uses global variable "numextras" that should be set before this function is called
// f(rho0,u)
// global function
void compute_allextras_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  int i;
  FTYPE du;

  if(whichdatatype[primarytable]==4){
    du = u - EOSextra[UNUGLOBAL];
  }
  else{
    du = u;
  }

  compute_allextras_du_kazfull(justnum, EOSextra, rho0, du, numextrasreturn, extras);


  // DEBUG:
  //  dualfprintf(fail_file,"compute_allextras_kazfull: rho0=%21.15g u=%21.15g du=%21.15g\n",rho0,u,du);
  //  int ei;
  //  for(ei=0;ei<MAXNUMEXTRAS;ei++) dualfprintf(fail_file,"extras[%d]=%21.15g\n",ei,extras[ei]);


}



// similar to compute_allextras_kazfull()
static void compute_notallextras1_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  FTYPE du;


  if(whichdatatype[primarytable]==4){
    du = u - EOSextra[UNUGLOBAL];
  }
  else{
    du = u;
  }

  compute_notallextras1_du_kazfull(justnum, EOSextra, rho0, du, numextrasreturn, extras);

  // DEBUG:
  //  dualfprintf(fail_file,"compute_notallextras1_kazfull: rho0=%21.15g u=%21.15g du=%21.15g\n",rho0,u,du);
  //  int ei;
  //  for(ei=0;ei<MAXNUMEXTRAS;ei++) dualfprintf(fail_file,"notallextras[%d]=%21.15g\n",ei,extras[ei]);

}



// similar to compute_allextras_du_kazfull()
static void compute_notallextras1_du_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras)
{
  int badlookups[MAXEOSPIPELINE];
  int iffun[MAXEOSPIPELINE];
  int whichtablesubtype=SUBTYPEEXTRA;
  int whichd=UTOTDIFF;
  int numcols,coli;
  int failreturn;


  *numextrasreturn = numcols = numcolintablesubtype[whichtablesubtype];// ok use of numcolintablesubtype

  
  if(justnum==0){

    // setup iffun
    // want not all extras
    for(coli=0;coli<numcols;coli++) iffun[coli]=0;
    // EXTRA1:  qtautnueohcm
    // EXTRA2:  qtauanueohcm
    // EXTRA3:  qtautnuebarohcm
    // EXTRA4:  qtauanuebarohcm
    // EXTRA5:  qtautmuohcm
    // EXTRA6:  qtauamuohcm
    // EXTRA13:  unue0
    // EXTRA14:  unuebar0
    // EXTRA15:  unumu0
    iffun[whichcolinquantity[EXTRA1]]=1;
    iffun[whichcolinquantity[EXTRA2]]=1;
    iffun[whichcolinquantity[EXTRA3]]=1;
    iffun[whichcolinquantity[EXTRA4]]=1;
    iffun[whichcolinquantity[EXTRA5]]=1;
    iffun[whichcolinquantity[EXTRA6]]=1;
    iffun[whichcolinquantity[EXTRA13]]=1;
    iffun[whichcolinquantity[EXTRA14]]=1;
    iffun[whichcolinquantity[EXTRA15]]=1;

    // do pipelined lookup
    failreturn=get_eos_fromtable(whichtablesubtype,iffun,whichd,EOSextra,rho0,du,extras,badlookups);

    // deal with bad lookups for extras
    reduce_extras(whichtablesubtype,iffun,whichd,EOSextra,rho0,du,extras,badlookups);

    int numbad=0;
    for(coli=0;coli<numcols;coli++) if(iffun[coli] && badlookups[coli]) numbad++;
    if(numbad==numcols) *numextrasreturn=0; // indicate that no point in doing more complex computations based upon these results -- can reduce to some simpler estimates

    // set rest to 0 for diagnostics to be uniform
    // These won't be used during calculation
    for(coli=numcols;coli<MAXNUMEXTRAS;coli++){
      extras[coli] = 0.0;
    }


    // DEBUG:
    //    dualfprintf(fail_file,"compute_notallextras1_du_kazfull: rho0=%21.15g du=%21.15g\n",rho0,du);
    //    for(coli=0;coli<numcols;coli++) dualfprintf(fail_file,"notallextras[%d]=%21.15g bad=%d\n",coli,extras[coli],badlookups[coli]);


  }// end if justnum==0


}





// f2c prototype
#include "f2c.h"
#include "tau_neededbyharm.P"
// not linking with libf2c since don't want that dependence and conversion doesn't need it since the original code was simple

// get neutrino terms (and return extras also)
// here function is (rho0,u) -- never used as function of p or chi
int get_extrasprocessed_kazfull(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  FTYPE quant1,quant2;
  // assume all tables have same number of extras or else this doesn't make sense in general
  FTYPE tempk,kbtk,rho0,u,HH,ye,ynu;
  int numextrasreturn,ei;
  FTYPE compute_temp_kazfull(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
  FTYPE unue0,unuebar0,unumu0,
    qtautnueohcm,qtautnuebarohcm,qtautmuohcm,
    qtauanueohcm,qtauanuebarohcm,qtauamuohcm,
    nnue0,nnuebar0,nnumu0,
    ntautnueohcm,ntautnuebarohcm,ntautmuohcm,
    ntauanueohcm,ntauanuebarohcm,ntauamuohcm,
    lambdatot,lambdaintot,
    tauphotonohcm,tauphotonabsohcm,
    nnueth0,nnuebarth0;
  FTYPE qphoton,qm,graddotrhouyl,tthermaltot,tdifftot,rho_nu,p_nu,s_nu,ynulocal,Ynuthermal,Ynuthermal0,enu,enue,enuebar;
  FTYPE qphoton_a[NUMHDIRECTIONS],qm_a[NUMHDIRECTIONS],graddotrhouyl_a[NUMHDIRECTIONS],tthermaltot_a[NUMHDIRECTIONS],tdifftot_a[NUMHDIRECTIONS],rho_nu_a[NUMHDIRECTIONS],p_nu_a[NUMHDIRECTIONS],s_nu_a[NUMHDIRECTIONS],ynulocal_a[NUMHDIRECTIONS],Ynuthermal_a[NUMHDIRECTIONS],Ynuthermal0_a,enu_a[NUMHDIRECTIONS],enue_a[NUMHDIRECTIONS],enuebar_a[NUMHDIRECTIONS];
  FTYPE dquant2;
  FTYPE qarray[NUMINDEPDIMENSMEM];
  int repeatedeos;
  int qi;
  int hi;
  FTYPE frac;
  int notintable;
  int whichd;
  FTYPE testextras[MAXNUMEXTRAS];


 
#if(0)
  // BEGIN DEBUG LOOP
  int loopit,numloops=1000;
  //int loopit,numloops=0;
  FTYPE logynu0;
  for(loopit=0;loopit<=numloops;loopit++){
    //if(doall) EOSextra[YNU0GLOBAL]=1E-15 + (1.0-1E-15)/(1000.0)*loopit;
    if(doall){
      logynu0=-9.0 + (1.0-(-9.0))/(1000.0)*loopit;
      EOSextra[YNU0GLOBAL]=pow(10.0,logynu0);
    }
    else numloops=0;

    // add the below to see *true* functional dependence
    // compute_Hglobal(GLOBALPOINT(EOSextraglobal),GLOBALPOINT(pglobal));


#endif
    // DEBUG:
    //  if((int)EOSextra[IGLOBAL]==0) dualfprintf(fail_file,"doall=%d\n",doall);
    //  else dualfprintf(fail_file,"ijk: %d %d %d\n",(int)EOSextra[IGLOBAL],(int)EOSextra[JGLOBAL],(int)EOSextra[KGLOBAL]);


#if(1)
    // DEBUG:
    // try changing density to a point in the table:
    //    pr[RHO]=3853528593.7105;
    //    pr[RHO]=3199267137.79737;

    // ensure HH is fixed to desired value or comment-out compute_Hglobal() call or assume read-in HH is good if only looking at first iteration
    //    EOSextra[HGLOBAL];
#endif

  
  

    ////////////////
    //
    // checks
    //
    ////////////////
    if(whichdatatype[primarytable]!=4){
      dualfprintf(fail_file,"Neutrino calculation not setup for other datatypes=%d primarytable=%d\n",whichdatatype[primarytable],primarytable);
      myexit(2954654);
    }


    //////////////////////
    //
    // Setup whichd, dquant2, and qarray
    //
    //////////////////////
    whichd=UTOTDIFF; // only one allowed for now
    quant1=pr[RHO];
    quant2=pr[UU];

    // first get dquant2 from quant2
    // avoid getting densities effectively twice, so just get UNUGLOBAL directly
    // Get dquant2 = u - u_\nu
    dquant2 = quant2 - EOSextra[UNUGLOBAL];


    // setup array
    qarray[RHOINDEP]=quant1;
    qarray[TEMPLIKEINDEP]=dquant2;
    // qarray[3+] are always stored in EOSextra
    for(qi=TEMPLIKEINDEP+1;qi<=LASTINDEPDIMENUSED;qi++){
      qarray[qi] = EOSextra[vartypeeosextraarray[qi]];
    }




    ////////////////////////////
    //
    // See if can get old version (repeatedeos)
    //
    // This is to avoid computations related to summing over directions and computing "processed" quantities, not just looking up "extras"
    //
    ////////////////////////////

    if(doallextrasold==doall){
      // check if repeated case (doesn't matter if i,j,k same since result only depends on all qarray values)
      repeatedeos=1;
      for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++){ // only need to repeat used independent variables, not all
        repeatedeos*=(fabs(qarray[qi]-qoldarrayextras[qi])<OLDTOLERANCE);
        //      if((int)EOSextra[IGLOBAL]==0) dualfprintf(fail_file,"repeatedcheck: qi=%d %21.15g %21.15g ig=%d\n",qi,qarray[qi],qoldarrayextras[qi],(int)EOSextra[IGLOBAL]);
      }

    }
    else{
      repeatedeos=0; // can't repeat if doall changed (i.e. if did doall=0, then doall=1 with repeated, then that info isn't there to repeat!)
      doallextrasold=doall; // store old version
    }




    // DEBUG:
    //  dualfprintf(fail_file,"repeatedeos=%d doall=%d doallextrasold=%d ig=%d\n",repeatedeos,doall,doallextrasold,(int)EOSextra[IGLOBAL]);





    ///////////////////////
    //
    // If repatedeos, get old values
    //
    ///////////////////////

    if(repeatedeos){
      if(doall){
        // then repeated case, so just return old result
        for(ei=0;ei<MAXNUMEXTRAS;ei++) extras[ei]=extrasold[ei];
        for(ei=0;ei<MAXPROCESSEDEXTRAS;ei++) processed[ei]=processedold[ei];
      }
      else{
        processed[RHONU]=processedold[RHONU];
        processed[PNU]=processedold[PNU];
        processed[SNU]=processedold[SNU];
      }
    }// end if repeatedeos
    else{


      /////////////////////////////
      //
      // else if not repeated
      //
      /////////////////////////////


      //////////////////////////
      //
      // set extra parameters used for lookup table
      //
      //////////////////////////
      rho0=qarray[RHOINDEP]; // "rho0" IS used
      u=qarray[TEMPLIKEINDEP]; // "u" is NOT used
      ye=qarray[YEINDEP]; // "ye" is NOT used
      if(WHICHEVOLVEYNU==EVOLVEYNURAD){
        ynu=pr[YNU]; // "ynu" is NOT used
      }
      else{
        // ynu is really Ynu0 that is from table lookup EOSextra[]
        ynu=qarray[YNUINDEP];
      }
      HH=qarray[HINDEP]; // "HH" IS used



   

      // not normal temperature(rho0,u), but temp(rho0,d), so direct lookup as tabulated
      // if used temp(rho0,u) then would get densities and this is already being done here, so avoid extra iteration
      // tempk is always stored in dimensionless form of: T[k] k_b/(m_b c^2), so since need kbtk in energy units, then always multiply tempk by m_b c^2 no matter what rho0unittype is
      tempk=compute_temp_whichd_kazfull(whichd, EOSextra, quant1, dquant2);
      kbtk=tempk*mbcsq; // in code energy units




      ///////////////
      //
      // Do lookup of "extra" quantities using to generate "processed" quantities later
      //
      ///////////////

      if(doall){

        // get extras (function of rho0,du as tabulated)
        // this also produces results when out of table
        compute_allextras_du_kazfull(0, EOSextra, quant1, dquant2, &numextrasreturn, extras);

        if(numextrasreturn==0) notintable=1;
        else notintable=0;

#if(0)
        // DEBUG:
        // 1) Use kazeos.loopparms.testynu.dek and edit tau_calc.f to output DEBUG info.
        // 2) Run stellar model generation and ensure first output of last iterated DEBUG output (choose tauiter largest for first set of iterations) is same!
        // 3) Run normal table generation ensure first output of last iterated DEBUG output is same as stellar model.  This confirms stellar model and non-stellar model agree!
        // 4) Copy last iterated result of tau_calc.f's DEBUG output into below structures with units conversions.  
        // 5) Ensure that compute_Hglobal() is commented out so HH calculation is not tested
        // 6) Run HARM with this debug enabled and ensure processed outputs are consistent with stellar model generated DEBUG output
        //
        // for:
        // rhob=3698268044.70451d0
        // tk=8183138022.48367d0
        // ye=0.428493d0
        // HH=0.38730771092310503125D+00008 (only used if stellar model, but processed quantities outputted for non-stellar model use it if you reset hcm (say) inside kazloopfunctions.f to be hcm = that quantity)
        // ynu=0.17028260300336636375E-00005 (only used if non-stellar model)
        // corresponding to:
        // ynu0=0.22703249033971029114E-00001
        //
        //
        // stellar model:
        //DEBUGFORHARM           6
        // extras  3.016208725417535E-009  1.780829277244639E-012  1.528763546663606E-008  1.252471070958565E-012  2.059938688437184E-009  1.314454808328260E-014  2.216369121583122E-009  1.040427344997658E-012  1.119680705812048E-008  7.376261380268049E-013  1.295212495261131E-009  1.099750846149906E-014  2.347702471508822E+026  6.181553876798139E+023  2.968941988288588E+025  5.074545552259924E+031  1.819268265786567E+029  8.338413581438167E+030   447251700.050036        20472315.5009781        3745028921.52609        3110006084.82794       1.117103214733134E+032  4.509919232490958E+028
        // processed   105624263274.740       1.257885287985835E+025  -2621445.40148080       7.332395341431752E-003  1.291919461573330E-003  3.101466000046774E+022  1.033822000015591E+022  3.660157476625981E+028  1.702826030033664E-006 3.759437713997883E-006  7.900220549935097E-006  7.918298450837450E-006  5.769081621748609E-006
        //
        // HELM table at the same HH (hacking kazloopfunctions.f so processed are computed with correct HH) and choosing now Ynu0
        //DEBUGFORHARM          23
        // extras  3.029339808920773E-009  1.602331902357635E-012  1.526768710185669E-008  1.257764285610078E-012  2.070317784546204E-009  1.310244034191246E-014  2.225741157400276E-009  9.353057661507511E-013  1.118029406428476E-008  7.374187614331142E-013  1.301632093683409E-009  1.084654641234010E-014  2.347694110104187E+026  6.181583853271378E+023  2.968941988288588E+025  5.074530210927057E+031  1.819277066697600E+029  8.338413581438167E+030   444961299.277883        20299760.1914007        3704290004.96078        3070417217.24895       1.016747562660904E+032  5.408080673791758E+028
        // processed   106785894260.181       1.132284286057749E+025  -2355920.25380538       7.394723355957793E-003  1.291919461573330E-003  2.793097681097268E+022  9.310325603657561E+021  3.296240345778634E+028  1.530676100627289E-006  3.076387364051629E-006  7.905717444283695E-006  7.925421990296669E-006  5.795089261262683E-006
        // so HELM and stellar model agree even though iterated on different quantity (Ynu vs. Ynu0).  Note that more advanced quasi-thermalization version of stellar model would not agree so is not really consistent to use as input to HARM as initial conditions.


#if(0)
        ei=-1;
        ei++;  testextras[ei]=3.029339808920773E-009/(1.0/Lunit); // 0.000760407388930137
        ei++;  testextras[ei]=1.602331902357635E-012/(1.0/Lunit); // 4.0220810306035e-07
        ei++;  testextras[ei]=1.526768710185669E-008/(1.0/Lunit); // 0.00383240666825727
        ei++;  testextras[ei]=1.257764285610078E-012/(1.0/Lunit); // 3.15716729266852e-07
        ei++;  testextras[ei]=2.070317784546204E-009/(1.0/Lunit); // 0.000519679217289016
        ei++;  testextras[ei]=1.310244034191246E-014/(1.0/Lunit); // 3.28889892763665e-09
        ei++;  testextras[ei]=2.225741157400276E-009/(1.0/Lunit); // 0.000558692695005398
        ei++;  testextras[ei]=9.353057661507511E-013/(1.0/Lunit); // 2.34775053427687e-07
        ei++;  testextras[ei]=1.118029406428476E-008/(1.0/Lunit); // 0.00280641286654554
        ei++;  testextras[ei]=7.374187614331142E-013/(1.0/Lunit); // 1.85102600004856e-07
        ei++;  testextras[ei]=1.301632093683409E-009/(1.0/Lunit); // 0.00032672817317846
        ei++;  testextras[ei]=1.084654641234010E-014/(1.0/Lunit); // 2.72263745784776e-09
        ei++;  testextras[ei]=2.347694110104187E+026/(energyunit/pow(Lunit,3.0)); // 261216.198320417
        ei++;  testextras[ei]=6.181583853271378E+023/(energyunit/pow(Lunit,3.0)); // 687.793962084254
        ei++;  testextras[ei]=2.968941988288588E+025/(energyunit/pow(Lunit,3.0)); // 33033.9346968924
        ei++;  testextras[ei]=5.074530210927057E+031/(1.0/pow(Lunit,3.0)); // 8.0258465358673e+47
        ei++;  testextras[ei]=1.819277066697600E+029/(1.0/pow(Lunit,3.0)); // 2.87735769354505e+45
        ei++;  testextras[ei]=8.338413581438167E+030/(1.0/pow(Lunit,3.0)); // 1.31879848922977e+47
        ei++;  testextras[ei]=444961299.277883/(Lunit); // 1772.65370767648
        ei++;  testextras[ei]=20299760.1914007/(Lunit); // 80.870954904681
        ei++;  testextras[ei]=3704290004.96078/(1.0/Lunit); // 929829490312520
        ei++;  testextras[ei]=3070417217.24895/(1.0/Lunit); // 770718402808099
        ei++;  testextras[ei]=1.016747562660904E+032/(1.0/pow(Lunit,3.0)); // 1.6080818449089e+48
        ei++;  testextras[ei]=5.408080673791758E+028/(1.0/pow(Lunit,3.0)); // 8.5533879467264e+44

        for(ei=0;ei<MAXNUMEXTRAS;ei++) extras[ei]=testextras[ei];
#endif

        // Without temperature correction, HARM produces:
        //repeatedeos=0 doall=1
        //qi=1 qarray[qi] =     3698369293.17202
        //qi=2 qarray[qi] =     30394298.6395949
        //qi=3 qarray[qi] =             0.428493
        //qi=4 qarray[qi] =    0.022703249033971
        //ynu= 1.69552232848888e-06
        //       processed[0]=  1.3478424217287e-15 -> p0*energyunit/Lunit**3/Tunit = 144678132296.455
        //       processed[1]=    0.105487137974405 -> 1.13230462681693e+25
        // processed[2]=    -19.7261783853841 -> p2*energyunit/Lunit**3/Tunit/c**2 = -2355930.23371244
        // processed[3]=     652.032446214747 -> p3*Tunit = 0.00545942420395928
        // processed[4]=     132.573286667649 -> p4*Tunit = 0.00111002729117781
        // processed[5]=     26.3533901908541 -> p5*energyunit/Lunit**3 = 2.36852459113022e+22
        // processed[6]=     8.78446339695136 -> 7.89508197043405e+21
        // processed[7]= 4.57801947359013e+44 -> p7/Lunit**3 = 2.89456047049061e+28
        // processed[8]=  1.3024258110165e-06 (ynulocal = ynu)
        // processed[9]= 2.61741093721727e-06 (ynuthermal)
        // processed[10]= 0.022703249033971   (ynuthermal0)
        // processed[11]= 5.56170975583112e-43 -> 7.90578019862632e-06
        // processed[12]= 5.57557196044043e-43
        // processed[13]=  4.0768867847036e-43
        //
        // So consistent within expected error of temperature.  Implies Ynu[Ynu0] looks good.  Original HELM good, so error must be in lookup table.
        // So error could be in either the degeneracy offset correction or in the fact that I use linear interpolation.
        // Tried 2X more temperature points and did NO better at getting correct Ynu0 even if Newton's method did work and find the correct Ynu with bad Ynu0.
        // Probably not density since new table is restricted to start at 1E7 instead of 1E5 and made no difference.

        // Below is results using "much" lower density by choosing density point on the 100x100 table
        //        qi=1 qarray[qi] =     3199267137.79737
        // qi=2 qarray[qi] =     30394298.6395949
        // qi=3 qarray[qi] =             0.428493
        // qi=4 qarray[qi] =    0.022703249033971
        // ynu= 1.69552232848888e-06
        // TOPLOT     0.022703249033971  1.50561094732914e-06
        // processed[0]=  1.3478424217287e-15
        // processed[1]=    0.105487137974405
        // processed[2]=    -19.7261783853841
        // processed[3]=     652.032446214747
        // processed[4]=     132.573286667649
        // processed[5]=     26.3533901908541
        // processed[6]=     8.78446339695136
        // processed[7]= 3.33354946920567e+44
        // processed[8]= 1.50561094732914e-06
        // processed[9]= 3.02574052771397e-06
        // processed[10]= 5.56170975583112e-43
        // processed[11]= 5.57557196044043e-43
        //        ....
        // processed[13]=  4.0768867847036e-43
        // So Ynu didn't change much with density, implying required change in Ynu0 to match this Ynu wouldn't change things much.
        // Nothing else apart from kbtk, rho, and extras create processed.  So appears look-up of extras is the issue.  Again, could be degeneracy offset or linear interpolation
        //


        // Below is with no extras preset and using original density, etc.


      
        //qi=1 qarray[qi] =     3698369293.17202
        //qi=2 qarray[qi] =     30394298.6395949
        //qi=3 qarray[qi] =             0.428493
        //qi=4 qarray[qi] =    0.022703249033971
        //ynu= 1.69552232848888e-06
        //TOPLOT     0.022703249033971  1.94527528169647e-05
        //extras[0]=  0.00104478340966707 // 3X the correct answer
        //extras[1]= 4.92067922628503e-06 // 10X
        //extras[2]=   0.0036609901318502 // 1X
        //extras[3]= 3.26397637238692e-07 // 1X
        //extras[4]=  0.00048233070306588 // 1X
        //extras[5]= 2.50686585390318e-09 // 1.5X
        //extras[6]= 0.000735490059278506 // 1.5X
        //extras[7]= 3.31139995770982e-06 // 10X
        //extras[8]=  0.00268750225013948 // 1X
        //extras[9]= 2.06638294070713e-07 // 1/10X
        //extras[10]= 0.000303372572818126 // 1X
        //extras[11]= 2.30282715881422e-09 // 1X
        //extras[12]=     279750.222129814 // 1X
        //extras[13]=     361.434906375252 // 1/2X
        //extras[14]=     28636.1055140615 // 1/1.5X
        //extras[15]= 8.31800467820582e+47 // 1X
        //extras[16]= 1.56104644759051e+45 // 1/2X
        //extras[17]= 1.18479626656395e+47 // 1X
        //extras[18]=     1233.32994603728 // 1/1.5X
        //extras[19]=     69.0394429060153 // 1/1.5X
        //extras[20]=      990532959445512 // 1X
        //extras[21]=      828017286074321 // 1X
        //extras[22]= 9.84688169261726e+47 // 1/1.5X
        //extras[23]= 5.82129582930178e+44 // 1/1.5X
        //processed[0]= 1.09679927385663e-15
        //processed[1]=     1.37508287270341
        //processed[2]=      -288.9191413193
        //processed[3]=     763.773349475091
        //processed[4]=     132.573286667649
        //processed[5]=     353.635529001977
        //processed[6]=     117.878509667326
        //processed[7]= 6.14323366595259e+45
        //processed[8]= 1.94527528169647e-05
        //processed[9]= 2.30304648159613e-05
        //processed[10]= 4.99481093689521e-43
        //processed[11]= 4.99543637701315e-43
        //  ....
        //processed[13]= 3.65705255437989e-43

        // Table lookup shows that values are most dependent upon Y_e, so need better table in Y_e and use higher-order interpolation in Y_e

#endif // end debug







        // MAXNUMEXTRAS entries
        // not same order as passed to function: computefinal_fromhcm()
        ei=0;
        qtautnueohcm=extras[ei]; ei++;
        qtauanueohcm=extras[ei]; ei++;
        qtautnuebarohcm=extras[ei]; ei++;
        qtauanuebarohcm=extras[ei]; ei++;
        qtautmuohcm=extras[ei]; ei++;
        qtauamuohcm=extras[ei]; ei++;

        ntautnueohcm=extras[ei]; ei++;
        ntauanueohcm=extras[ei]; ei++;
        ntautnuebarohcm=extras[ei]; ei++;
        ntauanuebarohcm=extras[ei]; ei++;
        ntautmuohcm=extras[ei]; ei++;
        ntauamuohcm=extras[ei]; ei++;

        unue0=extras[ei]; ei++;
        unuebar0=extras[ei]; ei++;
        unumu0=extras[ei]; ei++;

        nnue0=extras[ei]; ei++;
        nnuebar0=extras[ei]; ei++;
        nnumu0=extras[ei]; ei++;

        lambdatot=extras[ei]; ei++; // EXTRA19
        lambdaintot=extras[ei]; ei++; // EXTRA20

        tauphotonohcm=extras[ei]; ei++;
        tauphotonabsohcm=extras[ei]; ei++;

        nnueth0=extras[ei]; ei++;
        nnuebarth0=extras[ei]; ei++;

      }// end if doall
      else{// else if not doing all

        // get extras (only those needed)
        // get extras (function of rho0,du as tabulated)
        // also produces results when out of table
        compute_notallextras1_du_kazfull(0, EOSextra, quant1, dquant2, &numextrasreturn, extras);

        if(numextrasreturn==0) notintable=1;
        else notintable=0;

        // EXTRA1:  qtautnueohcm
        // EXTRA2:  qtauanueohcm
        // EXTRA3:  qtautnuebarohcm
        // EXTRA4:  qtauanuebarohcm
        // EXTRA5:  qtautmuohcm
        // EXTRA6:  qtauamuohcm
        // EXTRA13:  unue0
        // EXTRA14:  unuebar0
        // EXTRA15:  unumu0

        ei=0;
        qtautnueohcm=extras[ei]; ei++;
        qtauanueohcm=extras[ei]; ei++;
        qtautnuebarohcm=extras[ei]; ei++;
        qtauanuebarohcm=extras[ei]; ei++;
        qtautmuohcm=extras[ei]; ei++;
        qtauamuohcm=extras[ei]; ei++;

        ei++;
        ei++;
        ei++;
        ei++;
        ei++;
        ei++;

        unue0=extras[ei]; ei++;
        unuebar0=extras[ei]; ei++;
        unumu0=extras[ei]; ei++;

        ei++;
        ei++;
        ei++;

        ei++; // EXTRA19
        ei++; // EXTRA20

        ei++;
        ei++;

        ei++;
        ei++;
      }// end else if not doing all




  
      //////////////
      //
      // Compute "processed" quantities from "extra" quantities
      //
      // Comput per direction associated with HH
      //
      //////////////


      for(hi=0;hi<NUMHDIRECTIONS;hi++){
        // now process neutrino variables into final form
        HH = EOSextra[vartypeheightarray[hi+1]]; // vartypeheightarray[] index starts at 1

        if(doall){

          if(!notintable){



            // calling f2c generated code, which assumes a certain dimensionless form for inputs
            // rhob/mb should be a number density, and mb*rate/cc should be in rhounit*rate
            // this means if rho0 in c^2 units, then mb should be too
            // Note that Ynuthermal0 is not based upon optical depth, so below just repeats its calculation.
            computefinal_fromhcm__(&Ccode,&mbwithrhounit,&rho0,&kbtk,&HH,
                                   &unue0,&unuebar0,&unumu0,
                                   &qtautnueohcm,&qtautnuebarohcm,&qtautmuohcm,
                                   &qtauanueohcm,&qtauanuebarohcm,&qtauamuohcm,
                                   &nnue0,&nnuebar0,&nnumu0,
                                   &ntautnueohcm,&ntautnuebarohcm,&ntautmuohcm,
                                   &ntauanueohcm,&ntauanuebarohcm,&ntauamuohcm,
                                   &lambdatot,&lambdaintot,
                                   &tauphotonohcm,&tauphotonabsohcm,
                                   &nnueth0,&nnuebarth0,
                                   // outputs below
                                   &qphoton_a[hi],&qm_a[hi],&graddotrhouyl_a[hi],&tthermaltot_a[hi],&tdifftot_a[hi],&rho_nu_a[hi],&p_nu_a[hi],&s_nu_a[hi],&ynulocal_a[hi],&Ynuthermal_a[hi],&Ynuthermal0_a,&enu_a[hi],&enue_a[hi],&enuebar_a[hi]);
          }
          else{
            // then set to approximate optically thin values (i.e. ~0)
            qphoton_a[hi]=qm_a[hi]=graddotrhouyl_a[hi]=rho_nu_a[hi]=p_nu_a[hi]=s_nu_a[hi]=ynulocal_a[hi]=Ynuthermal_a[hi]=Ynuthermal0_a=enu_a[hi]=enue_a[hi]=enuebar_a[hi]=0.0;
            tthermaltot_a[hi]=tdifftot_a[hi]=BIG; // force non-evolution of Ynu
          }
        }// end if doall
        else{// else if not doall
          // DEBUG:
          // dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",qtautnueohcm,qtauanueohcm,qtautnuebarohcm,qtauanuebarohcm,qtautmuohcm,qtauamuohcm,unue0,unuebar0,unumu0);
          // dualfprintf(fail_file,"CCODE: :: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",Ccode,mbwithrhounit,rho0,kbtk,HH,unue0,unuebar0,unumu0,qtautnueohcm,qtautnuebarohcm,qtautmuohcm,qtauanueohcm,qtauanuebarohcm,qtauamuohcm);
          // dualfprintf(fail_file,"mb[cgs]=%21.15g T[K] = %21.15g\n",mbwithrhounit*Munit/Vunit/Vunit,kbtk*energyunit/kb);

          if(!notintable){
            // get rho_nu, p_nu, s_nu from extras
            computefinal_justdensities_fromhcm__(&Ccode,&mbwithrhounit,
                                                 &rho0,&kbtk,&HH,
                                                 &unue0,&unuebar0,&unumu0,
                                                 &qtautnueohcm,&qtautnuebarohcm,&qtautmuohcm,
                                                 &qtauanueohcm,&qtauanuebarohcm,&qtauamuohcm,
                                                 // outputs below
                                                 &rho_nu_a[hi],&p_nu_a[hi],&s_nu_a[hi]); // notice no & since already pointer
          }
          else{
            // not in table
            rho_nu_a[hi]=p_nu_a[hi]=s_nu_a[hi]=0.0;
          }
        }
      }
    


      /////////////////////////////////
      //
      // now get appropriate sum
      // if all h's same then get normal result, otherwise optical depth supression occurs per direction
      // hack for real method
      // GODMARK: a heating method might want access to these different directional quantities
      //
      //////////////////////////////////
      qphoton=qm=graddotrhouyl=tthermaltot=tdifftot=rho_nu=p_nu=s_nu=ynulocal=Ynuthermal=Ynuthermal0=enu=enue=enuebar=0.0;
      frac=1.0/((FTYPE)NUMHDIRECTIONS);
      for(hi=0;hi<NUMHDIRECTIONS;hi++){

        if(doall){
          qphoton += qphoton_a[hi]*frac;
          qm += qm_a[hi]*frac;
          graddotrhouyl += graddotrhouyl_a[hi]*frac;
          tthermaltot += tthermaltot_a[hi]*frac;
          tdifftot += tdifftot_a[hi]*frac;
          ynulocal += ynulocal_a[hi]*frac;
          Ynuthermal += Ynuthermal_a[hi]*frac;
          Ynuthermal0 = Ynuthermal0_a; // not based upon optical depth
          enu += enu_a[hi]*frac;
          enue += enue_a[hi]*frac;
          enuebar += enuebar_a[hi]*frac;
        }
      
        rho_nu += rho_nu_a[hi]*frac;
        p_nu += p_nu_a[hi]*frac;
        s_nu += s_nu_a[hi]*frac;

      }



      /////////////
      //
      // Store final results into processed[]
      //
      /////////////
  
      if(doall){
        // MAXPROCESSEDEXTRAS entries
        ei=-1; // just so easy copy/paste form
        ei++; processed[ei]=qphoton;
        ei++; processed[ei]=qm;
        ei++; processed[ei]=graddotrhouyl;
        ei++; processed[ei]=tthermaltot;
        ei++; processed[ei]=tdifftot;
        ei++; processed[ei]=rho_nu;
        ei++; processed[ei]=p_nu;
        ei++; processed[ei]=s_nu;
        ei++; processed[ei]=ynulocal;
        ei++; processed[ei]=Ynuthermal;
        ei++; processed[ei]=Ynuthermal0;
        // below are energies of *escaping* neutrinos
        ei++; processed[ei]=enu;
        ei++; processed[ei]=enue;
        ei++; processed[ei]=enuebar;
      }
      else{
        processed[RHONU]=rho_nu;
        processed[PNU]=p_nu;
        processed[SNU]=s_nu;
      }

      //  MAXPROCESSEDEXTRAS entries


      ////////////////////////////
      //
      // setup old values
      //
      ////////////////////////////
      for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++) qoldarrayextras[qi]=qarray[qi];
      if(doall){
        for(ei=0;ei<MAXNUMEXTRAS;ei++) extrasold[ei]=extras[ei];
        for(ei=0;ei<MAXPROCESSEDEXTRAS;ei++) processedold[ei]=processed[ei];
      }
      else{
        processedold[RHONU]=processed[RHONU];
        processedold[PNU]=processed[PNU];
        processedold[SNU]=processed[SNU];
      }


    }// end if not repeated lookup





#if(0)
    // DEBUG:
    dualfprintf(fail_file,"repeatedeos=%d doall=%d\n",repeatedeos,doall);

    // qarray[YEINDEP+] are always stored in EOSextra
    for(qi=FIRSTINDEPDIMEN;qi<=LASTINDEPDIMENUSED;qi++){
      dualfprintf(fail_file,"qi=%d qarray[qi] =%21.15g\n",qi,qarray[qi]);
    }

    dualfprintf(fail_file,"primitive ynu=%21.15g\n",pr[YNU]);

    if(doall) dualfprintf(fail_file,"TOPLOT %21.15g %21.15g\n",EOSextra[YNU0GLOBAL],processed[YNULOCAL]);
  
    if(doall){
      for(ei=0;ei<MAXNUMEXTRAS;ei++) dualfprintf(fail_file,"extras[%d]=%21.15g\n",ei,extras[ei]);
      for(ei=0;ei<MAXPROCESSEDEXTRAS;ei++) dualfprintf(fail_file,"processed[%d]=%21.15g\n",ei,processed[ei]);
    }
    else{
      dualfprintf(fail_file,"rhonu=%21.15g pnu=%21.15g snu=%21.15g\n",processed[RHONU],processed[PNU],processed[SNU]);
    }
#if(0)
    //  } // matches loopit loop
#endif
    myexit(0);

#endif


    return(0);


  }








  // assumes this is computed every timestep (or substep) or at least on some timescale that HH changes
  static void compute_upsnu_global(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    FTYPE rho0,u;
    FTYPE rho_nu, p_nu, s_nu;
    int i,j,k;



    if(whichdatatype[primarytable]!=4) return; // doesn't need to be done if !=4

    COMPFULLLOOP{
      compute_and_iterate_upsnu(UTOTDIFF, MAC(EOSextra,i,j,k), MAC(prim,i,j,k));
    }

  }


  // whichd = UTOTDIFF,PTOTDIFF,CHIDIFF
  // get neutrino rho_nu, p_nu, s_nu
  // function is (rho0,u/p/chi)
  // for now assume only needed as function of quant2=u.  Assumes this function called infrequently outside iterative routines
  // Point is that anyways we are doing an approximation, and save memory and gain simplicity by avoiding tabulating other f(rho0,chi) quantities
  static int compute_and_iterate_upsnu(int whichd, FTYPE *EOSextra, FTYPE *pr)
  {
    int get_extrasprocessed_kazfull(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
    int toreturn;
    FTYPE extras[MAXNUMEXTRAS];
    FTYPE processed[MAXPROCESSEDEXTRAS];


    // assume WHICHD==UTOTDIFF only
    toreturn=get_extrasprocessed_kazfull(0, EOSextra, pr, extras, processed);

    // iterate u_nu, p_nu, and s_nu
    iterateupsnu(processed,pr,EOSextra);

    return(toreturn);
  }




  // iterate u_nu, p_nu, and s_nu
  static int iterateupsnu(FTYPE *processed,FTYPE *pr,FTYPE *EOSextra)
  {

    /////////////////
    //
    // Iterate rho_nu, p_nu, and s_nu for table lookup
    //
    // assume DO WANT to update EOSextraglobal[UNUGLOBAL,PNUGLOBAL,SNUGLOBAL][i][j][k] since input is actual value on grid
    EOSextra[UNUGLOBAL] = processed[RHONU];
    EOSextra[PNUGLOBAL] = processed[PNU];
    EOSextra[SNUGLOBAL] = processed[SNU];

    return(0);

  }





  // assumes this is computed every timestep (or substep) or at least on some timescale that HH changes
  static void compute_ynu0_upsnu_global(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    FTYPE rho0,u;
    FTYPE rho_nu, p_nu, s_nu;
    int i,j,k;
    int numberofiterations;
    int iter;


#define NUMT0YNU0PARMSITER 10

    if(nstep==0){
      // then get converged solution beyond first step
      numberofiterations=NUMT0YNU0PARMSITER;
    }
    else{
      numberofiterations=NUMBEROFYNU0ITERATIONS ;
    }


    if(whichdatatype[primarytable]!=4) return; // doesn't need to be done if !=4

    COMPFULLLOOP{
      for(iter=0;iter<numberofiterations;iter++){
        compute_and_iterate_ynu0_upsnu(UTOTDIFF, MAC(EOSextra,i,j,k), MAC(prim,i,j,k));
      }
    }

  }





  // call this if not calling sources so can iterate Ynu0
  static int compute_and_iterate_ynu0_upsnu(int whichd, FTYPE *EOSextra, FTYPE *pr)
  {
    //  int whichd=UTOTDIFF;
    int get_extrasprocessed_kazfull(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
    int toreturn;
    FTYPE extras[MAXNUMEXTRAS];
    FTYPE processed[MAXPROCESSEDEXTRAS];



    // assume WHICHD==UTOTDIFF only
    // 1 below means doall extras and processed
    int doall=1;
    toreturn=get_extrasprocessed_kazfull(doall, EOSextra, pr, extras, processed);

    // iterate u_nu, p_nu, and s_nu
    iterateupsnu(processed,pr,EOSextra);

    // iterate Ynu0[Ynu]
    // only called when old Ynu evolution method
    if(WHICHEVOLVEYNU==EVOLVEYNURAD) iterateynu0(processed,pr,EOSextra);
    else EOSextra[YNUOLDGLOBAL]=processed[YNULOCAL]; // then not really old Ynu, really latest Ynu[Ynu0]


    return(toreturn);

  }


  // update unu,pnu,snu,Ynu0 (or Ynu) over many iterations (used by sources)
  static int get_extrasprocessed_manyiterations(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed)
  {
    int toreturn;
    int iter;


    for(iter=0;iter<NUMBEROFYNU0ITERATIONS;iter++){// iterate like when would otherwise call compute_ynu0_upsnu_global() that isn't called when cooling==COOLEOSGENERAL
    
      // get neutrino and photon quantities
      toreturn=get_extrasprocessed_kazfull(doall, EOSextra, pr, extras, processed); // (rho0,u)
    
      // if here (with cooling==COOLEOSGENERAL), then not doing iteration of upsnu or ynu0 otherwise, so do here
      // iterate u_nu, p_nu, and s_nu
      iterateupsnu(processed,pr,EOSextra);
    
      // iterate Ynu0[Ynu]
      // only called when old Ynu evolution method
      if(WHICHEVOLVEYNU==EVOLVEYNURAD) iterateynu0(processed,pr,EOSextra);
      else EOSextra[YNUOLDGLOBAL]=processed[YNULOCAL]; // then not really old Ynu, really latest Ynu[Ynu0]
    }


    // only care about last call to get_extrasprocessed_kazfull()
    return(toreturn);
  }






  // find Ynu0[Ynu]
  // iterate Ynu0
  // only can be done if doall was used when creating processed[], since otherwise processed[YNULOCAL] has not been recalled or computed
  static int iterateynu0(FTYPE *processed,FTYPE *pr,FTYPE *EOSextra)
  {
    FTYPE damp;
    FTYPE ynu0old=EOSextra[YNU0GLOBAL];
    FTYPE ynu0older=EOSextra[YNU0OLDGLOBAL];
    FTYPE ynuold=processed[YNULOCAL];
    FTYPE ynuolder=EOSextra[YNUOLDGLOBAL];
    FTYPE odRdYnu0,Rvar,errR;



    /////////////////
    //
    // Iterate Ynu0 for table lookup
    //
    
    // Take fake-simple Newton step (assumes monotonically increasing function and that Ynu0~Ynu)
    // Will use below to compute Ynu0[new] = Ynu[new][Ynu0,rho,du,Ye,H] + dYnu, where dYnu = -Ynu[old][Ynu0,rho,du,Ye,H] + Ynu0[old]
    // Recall that ynulocal = Ynu after table lookup (i.e. not Ynu0, but radiative transfer version of Ynu just computed)
    // Recall that ynu=pr[YNU] that is updated (i.e. new) Ynu from last evolution step of Ynu.
    // Recall that EOSextra[YNU0GLOBAL] = Ynu0 used for table lookup
    // Recall that ynulocal is estimate of Ynu using last Ynu0 but current rho,u,Ye.  So it corresponds to OLD value of Ynu
    // EOSextra[YNU0GLOBAL] = EOSextra[YNU0GLOBAL] + (ynu - ynulocal);    
    // EOSextra[YNU0GLOBAL] += (pr[ynu] - processed[YNULOCAL]);


    // TAKE NEWTON STEP:
    // In general, Newton's method is:
    // Resid = Rvar = (-Ynu[new] + Ynu[Ynu0old])
    // damp = starts at 1.0 and can decrease to avoid jumping too far (e.g. out of table)
    // dYnu0 = -damp*Rvar/(dR/dYnu0)
    // But dR/dYnu0 = dYnu[Ynu0old]/dYnu0 \approx (Ynu[Ynu0old] - Ynu[Ynu0older])/(Ynu0old-Ynu0older)
    // So dYnu0 = damp*(Ynu[new] - Ynu[Ynu0old])*(Ynu0old-Ynu0older)/(Ynu[Ynu0old] - Ynu[Ynu0older])
    // So need to store Ynu0older and Ynu[Ynu0older].  Have at first Ynu0old and computed above Ynu[Ynu0old]
    // Then Ynu0new = Ynu0old + dYnu0



  
    //    damp=1.0;
    damp=0.20;
    Rvar = (-pr[YNU] + ynuold);
    errR = fabs(Rvar/(fabs(pr[YNU])+fabs(ynuold)+SMALL));

    if((int)EOSextra[IGLOBAL]==0) dualfprintf(fail_file,"MARK: steppart=%d nstep=%ld Rvar=%21.15g errR=%21.15g\n",steppart,nstep,Rvar,errR);

    if(errR<0.05){
      // don't try to do better if already accurate to less than 5%
      
      if((int)EOSextra[IGLOBAL]==0) dualfprintf(fail_file,"MARK: ynu0old=%21.15g ynu0older=%21.15g : ynuold=%21.15g ynuolder=%21.15g : final=%21.15g\n",ynu0old,ynu0older,ynuold,ynuolder,EOSextra[YNU0GLOBAL]);

    }
    else{// then try to get better Ynu0

      if(ynuold==ynuolder || ynu0old==ynu0older){
        // If just first iteration, so no older information, so just use small slope to change so next iteration can get slope
        // If reached here after first iteration, then don't change Ynu0 much.  Just push so can get slope in case evolves
        // ok that I push it in same sense of Ynu-Ynulocal since just need a push to somewhere else

        // odRdYnu0 = 1/ (dR/dYnu0)
        // can't make guess step too small or else noise in H-integration calculation can trigger bad Ynu[Ynu0]
        odRdYnu0=-0.5*errR;
      }
      else{
        // standard derivative
        // odRdYnu0 = 1/ (dR/dYnu0)
        odRdYnu0=(ynu0old-ynu0older)/(ynuold-ynuolder);
      }

      // iterate
      EOSextra[YNU0GLOBAL] += -damp*Rvar*odRdYnu0;

      if((int)EOSextra[IGLOBAL]==0) dualfprintf(fail_file,"MARK: odRdYnu0=%21.15g ynu0old=%21.15g ynu0older=%21.15g : ynuold=%21.15g ynuolder=%21.15g : finalYnu0=%21.15g\n",odRdYnu0,ynu0old,ynu0older,ynuold,ynuolder,EOSextra[YNU0GLOBAL]);
      
      // need to limit Ynu0 to table (assumes reasonable behavior for above iteration!)
      if(EOSextra[YNU0GLOBAL]>0.999*lineartablelimits[primarytable][YNUEOS][1]) EOSextra[YNU0GLOBAL] = 0.999*lineartablelimits[primarytable][YNUEOS][1];
      if(EOSextra[YNU0GLOBAL]<1.001*lineartablelimits[primarytable][YNUEOS][0]) EOSextra[YNU0GLOBAL] = 1.001*lineartablelimits[primarytable][YNUEOS][0];
    }// end else if error large enough to need correction


    // store old value into older values
    EOSextra[YNU0OLDGLOBAL] = ynu0old;
    EOSextra[YNUOLDGLOBAL] = ynuold;


    return(0);
  }






  // compute sources to equations of motion for HARM
  // duplicate of compute_neutrino() in sources.c
  // Ui and dUother in UNOTHING form
  // Calling this iterates Ynu0 when doing whichdatatype==4 (compute_EOS_parms_kazfull() also iterates Ynu0)
  int compute_sources_EOS_kazfull(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
  {
    int ii,jj,kk,pp;
    int j;
    FTYPE X[NDIM],V[NDIM],r,th,Rcyl;
    FTYPE du;
    FTYPE rho,u,ye,yl,ynu;
    FTYPE cofactor;
    FTYPE dUdtau;
    FTYPE rph,photoncapture;
    FTYPE extras[MAXNUMEXTRAS];
    void compute_allextras_kazfull(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras);
    FTYPE processed[MAXPROCESSEDEXTRAS];
    int ei;
    FTYPE qphoton,qm,graddotrhouyl,tthermaltot,tdifftot,rho_nu,p_nu,s_nu,ynulocal,enu,enue,enuebar;
    FTYPE Ynuthermal,Ynuthermal0,lambdatot,lambdaintot;
    FTYPE graddotrhouynu;
    FTYPE dtau,dtlimit;
    FTYPE gradUU[NDIM];
    FTYPE Unewlocal[NPR];
    FTYPE Uylupdate;
    FTYPE Urhoupdate;
    FTYPE graddotrhouynulimit;
    int numextrasreturn;
    FTYPE Ynuthermaltouse;




    ///////////////////
    //
    // Get coordinate information
    //
    ///////////////////
    ii=geom->i;
    jj=geom->j;
    kk=geom->k;
    pp=geom->p;

    coord_ijk(ii,jj,kk,pp,X) ;
    bl_coord_ijk(ii,jj,kk,pp,V) ;

    r=V[1];
    th=V[2];
    Rcyl = r*sin(th) ;



    ///////////////////
    //
    // Assign primitives
    //
    ///////////////////
    rho = pr[RHO];
    u = pr[UU];

#if(DOYL!=DONOYL)
    //  yl = pr[YL];
    ye = pr[YE]; // now Ye used as primitive and YL as conserved
    // get true yl
    if(WHICHEVOLVEYNU==EVOLVEYNURAD) yl=ye + pr[YNU];
    else yl = ye + EOSextra[YNUOLDGLOBAL];
#else
    yl = 0.0;
#endif

#if(DOYNU!=DONOYNU)
    ynu = pr[YNU];    // if WHICHEVOLVEYNU==EVOLVEYNUNOTRAD, then pr[YNU] really Ynu0, not Ynu, and accounted for below when used.
#else
    ynu = 0.0;
#endif





    // approximately account for photon capture by black hole
    // 50% of masslesss particles from **stationary** isotropic emitters are captured by black hole
    // See Shapiro & Teukolsky (1983) and equation 2.81 and 2.82 in Shibata, Sekiguchi, Takahashi (2007)
    // Note that if MBH=0 then rph=0, so works for NS

    // For now, treat capture process as suppression of cooling/heating rates or any changes in fluid
    // In this way, capture by BH is treated as assuming fluid advected them in since assume if trapping photons then trapping fluid for sure.
    // In reality have to follow those photons/neutrinos to check if this is true, but good assumption.
    // a/MBH could be zero, so resolve so rph is not NaN when no BH so that photoncapture=1
    rph = 2.0*MBH*(1.0 + cos(2.0/3.0*(acos(-a/(fabs(MBH)+SMALL)))));
    photoncapture = (r>rph) ? 1.0 : 0.0;






    if(whichdatatype[primarytable]==4){

      // get extras and processed over potentially many iterations to obtain convergence
      get_extrasprocessed_manyiterations(1,EOSextra, pr, extras, processed);


      // MAXPROCESSEDEXTRAS entries
      qphoton=processed[QPHOTON];
      qm = processed[QNEUTRINO]; 
      graddotrhouyl = processed[GRADDOTRHOUYL]; 
      tthermaltot = processed[TTHERMAL]; 
      tdifftot = processed[TDIFF]; 
      rho_nu = processed[RHONU]; 
      p_nu = processed[PNU]; 
      s_nu = processed[SNU]; 
      ynulocal = processed[YNULOCAL];  // don't actually need this if WHICHEVOLVEYNU==EVOLVEYNURAD
      Ynuthermal = processed[YNUTHERMAL];
      Ynuthermal0 = processed[YNUTHERMAL0];
      // below are energies of *escaping* neutrinos
      enu = processed[ENUAVG]; 
      enue = processed[ENUE]; 
      enuebar = processed[ENUEBAR]; 

      // also store below additional as directly wanted quantities from extras
      lambdatot = extras[whichcolinquantity[EXTRA19]]; // extras[0] is start
      lambdaintot = extras[whichcolinquantity[EXTRA20]];


#if(DOYL!=DONOYL)
      // get corresponding update to yl since get_extrasprocessed_manyiterations() updated EOSextra[YNUOLDGLOBAL] and ynulocal
      ye = pr[YE]; // now Ye used as primitive and YL as conserved
      if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD) yl = ye + ynulocal;
#endif


    }
    else if(whichdatatype[primarytable]==3){

      compute_allextras_kazfull(0, EOSextra, rho, u, &numextrasreturn, extras);

      qphoton = extras[whichcolinquantity[EXTRA1]];
      qm = extras[whichcolinquantity[EXTRA2]];
      graddotrhouyl = extras[whichcolinquantity[EXTRA3]];
      tthermaltot = extras[whichcolinquantity[EXTRA4]];
      tdifftot = extras[whichcolinquantity[EXTRA5]];
      lambdatot = extras[whichcolinquantity[EXTRA6]];
      lambdaintot = extras[whichcolinquantity[EXTRA7]];
      enu = extras[whichcolinquantity[EXTRA8]];
      enue = extras[whichcolinquantity[EXTRA9]];
      enuebar = extras[whichcolinquantity[EXTRA10]];
      Ynuthermal = extras[whichcolinquantity[EXTRA11]];
      //    Ynuthermal0 = extras[whichcolinquantity[EXTRA12]]; // GODMARK: Should add if used.
    }



    ////////////////////////////////
    //
    // pick which Ynu evolution scheme to use
    //
    // note that pr[YNU] changes meaning
    //
    ////////////////////////////////
    if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD){
      // new scheme
      // here pr[YNU] is Ynu0
      Ynuthermaltouse=Ynuthermal0;
    }
    else{
      //  if(WHICHEVOLVEYNU==EVOLVEYNURAD){
      // old scheme
      // here pr[YNU] is Ynu
      Ynuthermaltouse=Ynuthermal;
    }


    // compute some other things
    //  dtau = dt/(q->ucon[TT]);






    /////////////
    //
    // Neutrino thermalization
    // Neutrino thermalization source term (due to neutrinos thermalizing fluid)
    //
    /////////////
    // regardless of the proper Lorentz transformation, this is the right dtau to use to limit ynu to Ynuthermaltouse in the fluid frame
    // well, not really quite true due to advection terms
    // since U[RHO] will get updated, no point in trying to strictly limit using present U[RHO] due to advection terms
    // limit graddotrhouynu so that ynu can't go beyond Ynuthermaltouse (haven't yet accounted for flux, so is not quite right constraint)
    // GODMARK: could add constraints later once flux known -- new function
    //  graddotrhouynu=rho*(Ynuthermaltouse-ynu)/MAX(dtau,tthermaltot);
    //  if(rho*ucon[TT]*ynu+graddotrhouynu*dt>rho*ucon[TT]*Ynuthermaltouse){
    //    graddotrhouynu = rho*(Ynuthermaltouse - ynu)/dtau;
    //  }

    // equation is \nablda_\mu (\rho u^\mu Y_\nu) = \rho dY_\nu/d\tau = \rho (Y_{\nu,thermal} - Y_\nu)/\tau_{thermal} = graddotrhouynu
    //  graddotrhouynu=rho*(Ynuthermaltouse-ynu)/MAX(dtau,tthermaltot);
    Urhoupdate = Ui[RHO] + (dUother[RHO])*dt; // assumes no physics source here for RHO (should add it if there is)

    // then based upon this criterion that doesn't necessarily imply exact dtau, strictly limit
    // Y_nu[new] = Unew[Ynu]/Unew[RHO] = (Ui[Ynu] + dUother[Ynu]*dt + dUsource[Ynu]*dt)/Unew[RHO] and solved for dUsource[Ynu]
    // this guarantees that post-advection an overshoot is truncated and leads to Ynu=Ynuthermaltouse
    graddotrhouynulimit = (Ynuthermaltouse*Urhoupdate - Ui[YNU])/dt - dUother[YNU];
  
    // dtlimit is chosen so the expression for graddotrhouynu is consistent
    // that is graddotrhoynu(original) = graddotrhouynulimit and solve for dtlimit
    dtlimit = rho*(Ynuthermaltouse-ynu)*sign(graddotrhouynulimit)/(fabs(graddotrhouynulimit)+SMALL);

    // GODMARK: should compare dtlimit with estimate of dtau -- for example can dtlimit<0?
    graddotrhouynu = rho*(Ynuthermaltouse-ynu)/MAX(dtlimit,tthermaltot);





    /////////////
    //
    // Lepton losses
    // Lepton source term (due to neutrinos escaping fluid)
    //
    // equation is:
    //
    // \nabla_\mu (\rho u^\mu Y_l) = \rho dY_l/d\tau = graddotrhouyl = m_b (Ndot_{nuebar} - Ndot_{nue})
    // where Ndot>0 means *loss* of that neutrino type.
    // So if Ndot_{nue}>0, then Y_e should decrease.
    //
    // Equation is Unew = Ui + dUother*dt + gradU*dt
    //
    /////////////
    // predicted update
    Uylupdate = Ui[YL] + (dUother[YL] + graddotrhouyl)*dt;

    // limit update so Uylupdate<=0
    // limit Y_l update (Am I sure Y_l can't go negative? GODMARK)
    // here limit is perfect since any change in U[RHO] doesn't matter since U[YL]=0 -> Y_l=0 for all U[RHO]
    if(Uylupdate<0.0){
      graddotrhouyl = -Ui[YL]/dt - dUother[YL]; // so yl=0 after full Riemann + geometry + this source update
    }




    /////////////
    //
    // Lepton+Photon Cooling
    // Note that qm>0 or qphoton>0 means loss of energy, so gradUU = +(qphoton+qm)u_\mu since for energy (j=0) gives overall loss of energy
    //
    /////////////
    // No obvious limit on gradUU since energy can be negative in ergosphere? and momentum and be + or - . GODMARK
    DLOOPA(j){
      gradUU[j]=(qphoton+qm)*(q->ucov[j]);
    }





    /////////////////////////////
    //
    // now apply source terms
    //
    // GODMARK: need a way to track photons that reach infinity vs. black hole
    //*photoncapture
    //
    // Previously had conditional: if(t>0){}else{} around dUcomp's.  Not sure why.
    //
    ////////////////////////////
    DLOOPA(j){
      // provide energy density per second cooled
      // source term has: \detg dU/d\tau u_\mu
      // dUcomp[][j+UU] because j=0 implies UU as required (i.e. dUcomp in NPR form, not NDIM form)
      dUcomp[RADSOURCE][j+UU]=photoncapture*gradUU[j]; // photons+neutrino energy lost in fluid frame.  Some may be captured by black hole, rest reaches infinity
    }


#if(DOYL!=DONOYL)
    j=YL; // accessing NPR-type quantity, so acces dUcomp[][j]
    dUcomp[RADSOURCE][j]=photoncapture*graddotrhouyl;
#endif

#if(DOYNU!=DONOYNU)
    j=YNU; // accessing NPR-type quantity, so acces dUcomp[][j]
    dUcomp[RADSOURCE][j]=photoncapture*graddotrhouynu;
#endif
  


    return(0) ;




  }














  // store parlist into EOSextra
  void store_EOS_parms_kazfull(int numparms, FTYPE *EOSextra, FTYPE *parlist)
  {
    FTYPE TDYNORYEtouse;
    FTYPE YNUtouse,YNUOLDtouse;
    FTYPE Htouse[NUMHDIRECTIONS];
    FTYPE unu,pnu,snu;
    int hi;




    if(numparms!=NUMEOSGLOBALS){
      dualfprintf(fail_file,"numparms=%d NUMEOSGLOBALS=%d\n",numparms,NUMEOSGLOBALS);
      myexit(917692);
    }


    //////////////
    //
    // Constraint EOSextra-related quantities
    //
    //////////////

    // must be right order when this function is used

    // YE assign
    TDYNORYEtouse=parlist[TDYNORYEGLOBAL-FIRSTEOSGLOBAL];

    // YNU0 assign
    YNUtouse=parlist[YNU0GLOBAL-FIRSTEOSGLOBAL];
    YNUOLDtouse=parlist[YNU0OLDGLOBAL-FIRSTEOSGLOBAL];
  
    // H assign
    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      Htouse[hi]=parlist[HGLOBAL-FIRSTEOSGLOBAL+hi]; // 5-3+hi and so for hi=0 gives 2 as required
    }
 
    // YE constrain
    if(TDYNORYEtouse>lineartablelimits[primarytable][YEEOS][1]) TDYNORYEtouse = lineartablelimits[primarytable][YEEOS][1];
    if(TDYNORYEtouse<lineartablelimits[primarytable][YEEOS][0]) TDYNORYEtouse = lineartablelimits[primarytable][YEEOS][0];

    // YNU0 (whichdatatype==4) or YNU (others) constrain
    if(YNUtouse>lineartablelimits[primarytable][YNUEOS][1]) YNUtouse = lineartablelimits[primarytable][YNUEOS][1];
    if(YNUtouse<lineartablelimits[primarytable][YNUEOS][0]) YNUtouse = lineartablelimits[primarytable][YNUEOS][0];
    if(YNUOLDtouse>lineartablelimits[primarytable][YNUEOS][1]) YNUOLDtouse = lineartablelimits[primarytable][YNUEOS][1];
    if(YNUOLDtouse<lineartablelimits[primarytable][YNUEOS][0]) YNUOLDtouse = lineartablelimits[primarytable][YNUEOS][0];

    // H constrain
    if(whichdatatype[primarytable]!=4){
      for(hi=0;hi<NUMHDIRECTIONS;hi++){
        if(Htouse[hi]>lineartablelimits[primarytable][HEOS][1]) Htouse[hi] = lineartablelimits[primarytable][HEOS][1];
        if(Htouse[hi]<lineartablelimits[primarytable][HEOS][0]) Htouse[hi] = lineartablelimits[primarytable][HEOS][0];
      }
    }




    //////////////
    //
    // Set EOSextra
    //
    //////////////


    // YE set  
    EOSextra[TDYNORYEGLOBAL] = TDYNORYEtouse;

    if((int)EOSextra[IGLOBAL]==-100 && WHICHEVOLVEYNU==EVOLVEYNURAD){
      // then not yet iterating on Ynu0, so start iteration

      // YNU0 set
      EOSextra[YNU0GLOBAL] = YNUtouse;
      EOSextra[YNU0OLDGLOBAL] = YNUOLDtouse;
      // YNUold set/assign
      EOSextra[YNUOLDGLOBAL]  = parlist[YNUOLDGLOBAL-FIRSTEOSGLOBAL]; // no particular limit on this quantity
    }

    if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD){
      // YNU0 globals still used for lookup
      EOSextra[YNU0GLOBAL] = YNUtouse;

      EOSextra[YNU0OLDGLOBAL] = YNUOLDtouse; // not used, but assign anyways

      // YNUold set/assign, which is really normal Ynu[Ynu0]
      EOSextra[YNUOLDGLOBAL]  = parlist[YNUOLDGLOBAL-FIRSTEOSGLOBAL]; // no particular limit on this quantity
    }


    // H set
    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      EOSextra[vartypeheightarray[hi+1]]=Htouse[hi];
    }

    // Unu,Pnu,Snu assign/set
    EOSextra[UNUGLOBAL]=parlist[UNUGLOBAL-FIRSTEOSGLOBAL];
    EOSextra[PNUGLOBAL]=parlist[PNUGLOBAL-FIRSTEOSGLOBAL];
    EOSextra[SNUGLOBAL]=parlist[SNUGLOBAL-FIRSTEOSGLOBAL];

    EOSextra[PGASGLOBAL]=parlist[PGASGLOBAL-FIRSTEOSGLOBAL];  // doesn't have to be set always since only used after inversion to check inversion with consistent p(\chi) instead of p(u)

    // IJKGLOBAL assign/set
    EOSextra[IGLOBAL]=parlist[IGLOBAL-FIRSTEOSGLOBAL];
    EOSextra[JGLOBAL]=parlist[JGLOBAL-FIRSTEOSGLOBAL];
    EOSextra[KGLOBAL]=parlist[KGLOBAL-FIRSTEOSGLOBAL];



  }


  // Put EOSextra into parlist
  void get_EOS_parms_kazfull(int*numparms, FTYPE *EOSextra, FTYPE *parlist)
  {
    int hi;

    // expected use of output should be know this ordering
    // resulting array parlist[] should start at index 0 and have "*numparms" elements
    parlist[TDYNORYEGLOBAL-FIRSTEOSGLOBAL]=EOSextra[TDYNORYEGLOBAL];

    parlist[YNU0GLOBAL-FIRSTEOSGLOBAL]=EOSextra[YNU0GLOBAL];

    // YNU0OLDGLOBAL not used if WHICHEVOLVEYNU==EVOLVEYNUNOTRAD, but still assign for consistency
    parlist[YNU0OLDGLOBAL-FIRSTEOSGLOBAL]=EOSextra[YNU0OLDGLOBAL];

    parlist[YNUOLDGLOBAL-FIRSTEOSGLOBAL]=EOSextra[YNUOLDGLOBAL];

    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      parlist[HGLOBAL-FIRSTEOSGLOBAL+hi]=EOSextra[vartypeheightarray[hi+1]];
    }

    parlist[UNUGLOBAL-FIRSTEOSGLOBAL]=EOSextra[UNUGLOBAL];
    parlist[PNUGLOBAL-FIRSTEOSGLOBAL]=EOSextra[PNUGLOBAL];
    parlist[SNUGLOBAL-FIRSTEOSGLOBAL]=EOSextra[SNUGLOBAL];

    parlist[PGASGLOBAL-FIRSTEOSGLOBAL]=EOSextra[PGASGLOBAL]; // probably defined, but not necessarily always.  Just used for p(chi) inversion check

    parlist[IGLOBAL-FIRSTEOSGLOBAL]=EOSextra[IGLOBAL];
    parlist[JGLOBAL-FIRSTEOSGLOBAL]=EOSextra[JGLOBAL];
    parlist[KGLOBAL-FIRSTEOSGLOBAL]=EOSextra[KGLOBAL];

    *numparms=MAXPARLIST; // Only Ye -> ???

  }




  // compute things beyond simple EOS independent variables
  void compute_EOS_parms_kazfull(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    static int firsttime=1; // should be ok to have static here since not inside OpenMP parallel region

    ////////////////////////
    //
    // set grid indices
    //
    ////////////////////////

    compute_IJKglobal(EOSextra);

    /////////////////////////
    //  
    // compute H
    //
    /////////////////////////

    compute_Hglobal(EOSextra,prim);


    /////////////////////////
    //  
    // constrain Tdyn or Ye depending upon whichrnpmethod
    // Also constrain Ynu0
    // GODMARK: Now redundant, so comment out
    // Now Ye and Ynu0 are constrained where computed
    //  
    /////////////////////////

    ////  constrain_TDYNORYE_YNU_global(EOSextra,prim);

    /////////////////////////
    //  
    // 1) Iterate Ynu0 if whichdatatype==4
    // 2) Compute neutriono u_\nu, p_\nu, and s_\nu
    // only do this if no source term.  Otherwise, source term calls processed quantities so use this.
    // This also constrains resulting Ynu0
    //  
    /////////////////////////

    if(cooling!=2 || firsttime) compute_ynu0_upsnu_global(EOSextra,prim);
  

    // indicate no longer firstime here
    firsttime=0;

  }


#define NUMT0FULLPARMSITER 10


  // compute things beyond simple EOS independent variables
  // "full" really means at t=0 for now, so do multiple times.
  void compute_EOS_parms_full_kazfull(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {

    ////////////////////
    //
    // do normal update 10 times so converged solution at t=0
    //
    ////////////////////
    int i;
    for(i=0;i<NUMT0FULLPARMSITER;i++){
      compute_EOS_parms_kazfull(EOSextra,prim);
    }


  }


  static void compute_IJKglobal(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS])
  {
    int i,j,k;
    COMPFULLLOOP{
      MACP0A1(EOSextra,i,j,k,IGLOBAL) = i;
      MACP0A1(EOSextra,i,j,k,JGLOBAL) = j;
      MACP0A1(EOSextra,i,j,k,KGLOBAL) = k;
    }
  }


  static void compute_Hglobal_new(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {

    // 1) Take EOSextra, rho0, u from full overall-CPUs grid to a single smaller grid
    // 2) Ray-trace out from emission zones to determine their optical depth : each CPU gets some fraction of rays
    // 3) Take final per-CPU rays' resulting \tau and feed back to all CPUs.
    // 4) Interpolate back \tau onto hydro grid.


    // get lambdatot
    //get_lambdatot(MAC(EOSextra,i,j,k),rho0,u,&lambdatot);



  }




  // SUPERTODO: Need to integrate up all flavors and absorption/scattering with different summations to get different H's.  Then need to feed these different H's into 2-stream.  Maybe that will help 20X error in Y_e evolution?
  // assumes this is computed every timestep (or substep) or at least on some timescale that H changes
  // EOSextra is not input here -- rather part of output.
  static void compute_Hglobal(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    int i,j,k;
    int ii,jj,kk;
    int pl,pliter;
    FTYPE r,th,phi;
    FTYPE X[NDIM],V[NDIM],dxdxp[NDIM][NDIM];
    FTYPE dr,rdth,rsinthdphi,lambdatot;
    FTYPE rho0,u;
    FTYPE Htest1,Htest2,Htest3,Htest4;
    int lambdatotextra;


    //  return; //assume initial H is good -- use to compare input and output for EOS

    // this calculation makes sense only if using spherical polar coordinates
    if(!ISSPCMCOORD(MCOORD)){
      dualfprintf(fail_file,"Hglobal not setup for anything except spherical polar coordinates: %d\n",ISSPCMCOORD(MCOORD));
      myexit(2673);
    }


    // first store lookups and computations per point
    COMPFULLLOOP{

      //////////////////////////
      // get geometry stuff
      coord_ijk(i,j,k,CENT,X); // doesn't matter if CENT or anyother thing is used since just an estimate
      bl_coord_ijk(i,j,k,CENT,V);
      dxdxprim_ijk(i,j,k,CENT,dxdxp);
      r=V[1];
      th=V[2];
      phi=V[3];
      // assumes original coordinates are r,\theta,\phi
      dr = MAX(fabs(dxdxp[1][1]*dx[1]),fabs(dxdxp[1][2]*dx[2])); // just an estimate -- assume dxdxp[1][1]~1
      rdth = r*MAX(fabs(dxdxp[2][1]*dx[1]),fabs(dxdxp[2][2]*dx[2])); // just an estimate -- assume dxdxp[2][2]~1
#if(0)
      rsinthdphi = r*sin(th)*fabs(dxdxp[3][3]*dx[3]); // just an estimate -- assume dxdxp[3][3]~1
#endif

      rho0=MACP0A1(prim,i,j,k,RHO);
      u=MACP0A1(prim,i,j,k,UU);


      // get lambdatot
      get_lambdatot(MAC(EOSextra,i,j,k),rho0,u,&lambdatot);
      lambdatot=fabs(lambdatot)+SMALL; // ensure positive definite
    


      // get d\tau/dL
      GLOBALMACP0A1(ptemparray,i,j,k,0) = lambdatot;
      GLOBALMACP0A1(ptemparray,i,j,k,1) = dr/lambdatot; // dtau for radial integral
      GLOBALMACP0A1(ptemparray,i,j,k,2) = rdth/lambdatot; // dtau for angular integral
#if(0)
      GLOBALMACP0A1(ptemparray,i,j,k,3) = rsinthddphi/lambdatot; // dtau for radial integral
#endif

    
    }





    // GODMARK : NOT YET FOR MPI
    COMPFULLLOOP{ 


      //////////////////////////////
      //
      // get +- r direction H
      Htest1 = Htest2 = 0.0;
      jj=j;
      kk=k;
      for(ii=i;ii<N1+N1BND;ii++){// outward pointing integral (in MPI, should include other integrals that go to r=infinity)
        Htest1 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,1);
      }
#if(0)
      // convergence of neutrinos to r=0 makes no sense, so only use outgoing rays as estimate
      Htest2 = GLOBALMACP0A1(ptemparray,0,jj,kk,1); // double count ii=0
      for(ii=0;ii<i;ii++){
        Htest2 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,1);
      } 
      MACP0A1(EOSextra,i,j,k,H2GLOBAL) = (Htest1+Htest2*2.0)*GLOBALMACP0A1(ptemparray,i,j,k,0);
#endif
      // This is scale-height used by Kaz's EOS
      MACP0A1(EOSextra,i,j,k,HGLOBAL)  = Htest1*GLOBALMACP0A1(ptemparray,i,j,k,0);
      MACP0A1(EOSextra,i,j,k,H2GLOBAL) = Htest1*GLOBALMACP0A1(ptemparray,i,j,k,0);



      //////////////////////////////
      //
      // get +- \theta direction H (assume spherical polar)
      // GODMARK : NOT YET FOR MPI
      Htest3 = Htest4 = 0.0;
      ii=i;
      kk=k;
      for(jj=j;jj<=N2-1;jj++){
        Htest3 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,2);
      }
      for(jj=j;jj>=0;jj--){
        Htest4 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,2);
      }
      // This is scale-height used by Kaz's EOS
      // add angular part to hack photon trajectory for +-z for disk near BH or NS
      MACP0A1(EOSextra,i,j,k,H3GLOBAL) = (Htest1 + Htest3)*GLOBALMACP0A1(ptemparray,i,j,k,0);

      MACP0A1(EOSextra,i,j,k,H4GLOBAL) = (Htest1 + Htest4)*GLOBALMACP0A1(ptemparray,i,j,k,0);

      //    MACP0A1(EOSextra,i,j,k,H3GLOBAL) = (Htest3 + SHIFT2*MACP0A1(EOSextra,i,SHIFT2,k,HGLOBAL)/(GLOBALMACP0A1(ptemparray,i,SHIFT2,k,0)+SMALL))*GLOBALMACP0A1(ptemparray,i,j,k,0);
      //    MACP0A1(EOSextra,i,j,k,H4GLOBAL) = (Htest4 + SHIFT2*MACP0A1(EOSextra,i,N2-1-SHIFT2,k,HGLOBAL)/(GLOBALMACP0A1(ptemparray,i,N2-1-SHIFT2,k,0)+SMALL))*GLOBALMACP0A1(ptemparray,i,j,k,0);

    }




#if(0)
    //////////////////////////////
    //
    // get +- \phi direction H (assume spherical polar)
    // GODMARK : NOT YET FOR MPI
    COMPFULLLOOP{ 
      Htest1 = Htest2 = 0.0;
      ii=i;
      jj=j;
      for(kk=k;kk<N3;kk++){
        Htest1 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,1);
      }
      for(kk=k;kk>=0;kk--){
        Htest2 +=GLOBALMACP0A1(ptemparray,ii,jj,kk,1);
      }
      // This is scale-height used by Kaz's EOS
      MACP0A1(EOSextra,i,j,k,H5GLOBAL) = Htest1*GLOBALMACP0A1(ptemparray,i,j,k,0);
      MACP0A1(EOSextra,i,j,k,H6GLOBAL) = Htest2*GLOBALMACP0A1(ptemparray,i,j,k,0);
    }
#endif




    ///////////////
    //
    // constrain resulting H
    //
    ///////////////
    constrain_H(EOSextra,prim);



  }



  // get lambdatot(rho0,u)
  static void get_lambdatot(FTYPE *EOSextra, FTYPE rho0, FTYPE u, FTYPE *lambdatotptr)
  {
    int whichd=UTOTDIFF;
    int whichfun;
    FTYPE unu,snu,pnu,du;



    /////////////////////
    //
    // determine what whichfun0 means for this table
    // could put this in lookup with primarytable->whichtable, but expensive to put there and assume all tables same whichdatatype, so this is ok
    //
    /////////////////////
    if(whichdatatype[primarytable]==3){
      whichfun=EXTRA5;
    }
    else if(whichdatatype[primarytable]==4){
      whichfun=EXTRA19;
    }
    else{
      dualfprintf(fail_file,"Shouldn't request whichfun0=%d if primarytable=%d\n",LAMBDATOT,primarytable);
      myexit(46763463);
    }

    if(whichdatatype[primarytable]==4){
      unu = EOSextra[UNUGLOBAL];
      pnu = EOSextra[PNUGLOBAL];
      snu = EOSextra[SNUGLOBAL];
      du  = u - unu;
    }
    else{
      du = u;
    }


    if(getsingle_eos_fromtable(whichfun,whichd,EOSextra,rho0,du,lambdatotptr)){ // uses utotdiff=du
      *lambdatotptr=BIG; // then assume optically thin (worry if outside when rho>>the limit in the table?) GODMARK
    }



  }




  static void constrain_H(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    int hi;
    int i,j,k;

    if(whichdatatype[primarytable]!=4){ // otherwise no need to limit since not using table lookup for H
      for(hi=0;hi<NUMHDIRECTIONS;hi++){
        COMPFULLLOOP{ 
 
 
          // since there exists no information beyond bounds of table, limit H to table values so return value is consistent
          // This is also needed since H is used with other EOS table values and needs to be consistent
          // Note that when using simpletable that these restrictions won't matter
          // offset from strict linear table limits because log-linear conversion leads to error
          if(MACP0A1(EOSextra,i,j,k,HGLOBAL+hi)>0.999*lineartablelimits[primarytable][HEOS][1]) MACP0A1(EOSextra,i,j,k,HGLOBAL+hi) = 0.999*lineartablelimits[primarytable][HEOS][1];
          if(MACP0A1(EOSextra,i,j,k,HGLOBAL+hi)<1.001*lineartablelimits[primarytable][HEOS][0]) MACP0A1(EOSextra,i,j,k,HGLOBAL+hi) = 1.001*lineartablelimits[primarytable][HEOS][0];
        }
      }
    }


  }










  // constrain Y_e,Ynu
  // GODMARK: no longer used since Ye and Ynu0 are constrained where computed
  static void constrain_TDYNORYE_YNU(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
  {
    int i,j,k;
    FTYPE TDYNORYEtouse, YNUtouse;


#if(DOYL==DONOYL || DOYNU==DONOYNU)
    return;
#endif


    COMPFULLLOOP{
    
#if(! (DOYL==DONOYL && DOYNU==DONOYNU) )
      TDYNORYEtouse = MACP0A1(EOSextra,i,j,k,TDYNORYEGLOBAL);

      if(whichdatatype[primarytable]==4){
        // Below really YNU0
        YNUtouse = MACP0A1(EOSextra,i,j,k,YNU0GLOBAL);
      }
      else{
        dualfprintf(fail_file,"whichdatatype==3 method not supported currently.  Not dealing with YNU and constraining YNU correctly yet.\n");
        myexit(9247221);
        YNUtouse = MACP0A1(prim,i,j,k,YNU);
      }
#endif  

      // since there exists no information beyond bounds of table, limit TDYNORYE to table values so return value is consistent
      // This is also needed since TDYNORYE is used with other EOS table values and needs to be consistent
      // Note that when using simpletable that these restrictions won't matter
      // offset from strict linear table limits because log-linear conversion leads to error
      if(TDYNORYEtouse>0.999*lineartablelimits[primarytable][YEEOS][1]) TDYNORYEtouse = 0.999*lineartablelimits[primarytable][YEEOS][1];
      if(TDYNORYEtouse<1.001*lineartablelimits[primarytable][YEEOS][0]) TDYNORYEtouse = 1.001*lineartablelimits[primarytable][YEEOS][0];

      if(YNUtouse>0.999*lineartablelimits[primarytable][YNUEOS][1]) YNUtouse = 0.999*lineartablelimits[primarytable][YNUEOS][1];
      if(YNUtouse<1.001*lineartablelimits[primarytable][YNUEOS][0]) YNUtouse = 1.001*lineartablelimits[primarytable][YNUEOS][0];

      // make final assignment
      MACP0A1(EOSextra,i,j,k,TDYNORYEGLOBAL)=TDYNORYEtouse;
      MACP0A1(EOSextra,i,j,k,YNU0GLOBAL)=YNUtouse;        // Really YNU0 if whichdatatype==4

    }


  }






  // Change the primitives to be limited by constraints on Y_e
  // Do this so if out of bounds, Y_e can still recover eventually if fluid goes out of region where Y_\nu is large
  // assumes this is computed every timestep (or substep) before anything else computed
  // actual computation of Y_e[Y_L,Y_\nu]
  void fix_primitive_eos_scalars_kazfull(FTYPE *EOSextra, FTYPE *pr)
  {
    FTYPE yemin,yemax;
    FTYPE ynu0min,ynu0max;
    FTYPE ylmin,ylmax;
    FTYPE ye;
    FTYPE trueynu,trueyl;




    yemin=1.001*lineartablelimits[primarytable][YEEOS][0];
    yemax=0.999*lineartablelimits[primarytable][YEEOS][1];
    ynu0min=1.001*lineartablelimits[primarytable][YNUEOS][0];
    ynu0max=0.999*lineartablelimits[primarytable][YNUEOS][1];
    // No, yl is not sum of ye and ynu0.  So don't limit yl or ynu
    //  ylmin=yemin+ynumin;
    //  ylmax=yemax+ynumax;


    /////////////
    //
    // update ye from new primitives
    //
    /////////////

    if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD){
      // use latest calculated Ynu[Ynu0] from directly evolve Ynu0
      trueynu=EOSextra[YNUOLDGLOBAL];
    }
    else{
      // then pr[YNU] is really Ynu[Ynu0] already that is Ynu directly evolved
      trueynu=pr[YNU];
    }


    /////////////
    //
    // update Y_e using true Y_l and true Y_\nu
    // Note that Y_e = Y_l - Y_\nu such that Y_e>0 and Y_e<1 as handled in EOS
    // set Y_L
    //
    /////////////
    ye = pr[YE]; // not Y_e is primitive and Y_L is conserved. ///pr[YL] - trueynu;
    trueyl = ye + trueynu;




    ////////////
    //
    // constrain Y_e to table values
    //
    ////////////

    // have trueyl and trueynu share the blame
    if(ye<yemin){
      trueyl  += -0.5*(ye-yemin);
      trueynu += +0.5*(ye-yemin);
      // now effective ye is yemin
      ye=yemin;
      EOSextra[TDYNORYEGLOBAL] = pr[YE] = yemin;
    }
    else if(ye>yemax){
      trueyl  += -0.5*(ye-yemax);
      trueynu += +0.5*(ye-yemax);
      // Now effective ye is yemax
      EOSextra[TDYNORYEGLOBAL] = pr[YE] = yemax;
    }


    ////////////
    //
    // need to limit Ye to table in case changed (really already limited correctly, but let's make sure!)
    //
    ////////////
    if(EOSextra[TDYNORYEGLOBAL]>0.999*lineartablelimits[primarytable][YEEOS][1]) EOSextra[TDYNORYEGLOBAL] = 0.999*lineartablelimits[primarytable][YEEOS][1];
    if(EOSextra[TDYNORYEGLOBAL]<1.001*lineartablelimits[primarytable][YEEOS][0]) EOSextra[TDYNORYEGLOBAL] = 1.001*lineartablelimits[primarytable][YEEOS][0];



    //////////////////
    //
    // update Y_nu
    //
    //////////////////
    if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD){
      // pr[YNU] is Ynu0, so must use generated Ynu[Ynu0] from processed quantities that was stored in EOSextra
      // SUPERTODO: If above constraints changed things, then should obtain new Ynu0[Ynu] so consistent Ynu and Ynu0 used!  Otherwise only Ye has memory of constraint and trueynu value will become lost information relative to Ynu0
      EOSextra[YNU0GLOBAL] = pr[YNU];
    }
    else{
      // pr[YNU] is latest Ynu, so correct update to Y_e from Y_L
      pr[YNU] = trueynu;
      //    EOSextra[YNU0GLOBAL] obtained only after iteration from Ynu0[Ynu]
    }


    // need to limit Ynu0 to table in case changed
    if(EOSextra[YNU0GLOBAL]>0.999*lineartablelimits[primarytable][YNUEOS][1]) EOSextra[YNU0GLOBAL] = 0.999*lineartablelimits[primarytable][YNUEOS][1];
    if(EOSextra[YNU0GLOBAL]<1.001*lineartablelimits[primarytable][YNUEOS][0]) EOSextra[YNU0GLOBAL] = 1.001*lineartablelimits[primarytable][YNUEOS][0];
  

  }



  // get quantity that will be prefactor for mass flux in flux associated with Y_L conservation
  void yl2advect_kazfull(FTYPE *EOSextra, FTYPE pryl, FTYPE prynu, FTYPE *prforadvect)
  {
    FTYPE ye,trueynu,trueyl;



    // pr[YL] is really Y_e now
    // U[YL] is Y_L itself

    // get true ye
    EOSextra[TDYNORYEGLOBAL]=ye=pryl; // pryl=pr[YL]=pr[YE] is really Y_e

    // need to limit Ye to table in case changed due to interpolation or evolution
    if(EOSextra[TDYNORYEGLOBAL]>0.999*lineartablelimits[primarytable][YEEOS][1]) EOSextra[TDYNORYEGLOBAL] = 0.999*lineartablelimits[primarytable][YEEOS][1];
    if(EOSextra[TDYNORYEGLOBAL]<1.001*lineartablelimits[primarytable][YEEOS][0]) EOSextra[TDYNORYEGLOBAL] = 1.001*lineartablelimits[primarytable][YEEOS][0];



    // get true ynu
    if(WHICHEVOLVEYNU==EVOLVEYNURAD){
      trueynu = prynu;
    }
    else{
      // not really old Ynu, actually latest Ynu for this method
      trueynu = EOSextra[YNUOLDGLOBAL];
    }


    // get true Y_L
    trueyl = ye+trueynu;


    // setup Y_L quantity to advect
    *prforadvect=trueyl;
  

  }


  // get primitive associated with conserved quantities
  void advect2yl_kazfull(FTYPE *EOSextra, FTYPE ylforadvect, FTYPE ynuforadvect, FTYPE *ye)
  {
    FTYPE trueynu,trueyl;


    trueyl=ylforadvect;
    if(WHICHEVOLVEYNU==EVOLVEYNURAD){
      trueynu=ynuforadvect;
    }
    else{
      // GODMARK: could update Ynu[Ynu0] here, but not necessary?
      trueynu=EOSextra[YNUOLDGLOBAL];
    }

    // pr[YL] is really Y_e now
    // U[YL] is Y_L itself

    // get true ye
    EOSextra[TDYNORYEGLOBAL]= *ye =trueyl-trueynu;

    // need to limit Ye to table in case changed
    if(EOSextra[TDYNORYEGLOBAL]>0.999*lineartablelimits[primarytable][YEEOS][1]) EOSextra[TDYNORYEGLOBAL] = 0.999*lineartablelimits[primarytable][YEEOS][1];
    if(EOSextra[TDYNORYEGLOBAL]<1.001*lineartablelimits[primarytable][YEEOS][0]) EOSextra[TDYNORYEGLOBAL] = 1.001*lineartablelimits[primarytable][YEEOS][0];



  }



  // get prefactor of massflux for conserved quantity associated with Y_\nu conservation
  void ynu2advect_kazfull(FTYPE *EOSextra, FTYPE pryl, FTYPE prynu, FTYPE *prforadvect)
  {
    FTYPE trueynuevolve;



    // get true ynu that we are evolving
    if(WHICHEVOLVEYNU==EVOLVEYNURAD){
      // Evolving prynu=pr[YNU] is Ynu
      trueynuevolve = prynu;
    }
    else{
      // evolving Ynu0 is prynu=pr[YNU]
      trueynuevolve = EOSextra[YNU0GLOBAL] =  prynu;

      // need to limit Ynu0 to table in case changed
      if(EOSextra[YNU0GLOBAL]>0.999*lineartablelimits[primarytable][YNUEOS][1]) EOSextra[YNU0GLOBAL] = 0.999*lineartablelimits[primarytable][YNUEOS][1];
      if(EOSextra[YNU0GLOBAL]<1.001*lineartablelimits[primarytable][YNUEOS][0]) EOSextra[YNU0GLOBAL] = 1.001*lineartablelimits[primarytable][YNUEOS][0];

    }

    // setup quantity to advect
    *prforadvect=trueynuevolve;


  }




  // get primitive associated with conserved Y_\nu
  void advect2ynu_kazfull(FTYPE *EOSextra, FTYPE ylforadvect, FTYPE ynuforadvect, FTYPE *prynu)
  {
    FTYPE trueynuevolve;


    // in either case, advected quantity is primitive quantity
    if(WHICHEVOLVEYNU==EVOLVEYNURAD){
      // prynu is Ynu
      *prynu = ynuforadvect;
      // GODMARK: could update Ynu0[Ynu]
    }
    else{
      // prynu is Ynu0
      *prynu = ynuforadvect;

      EOSextra[YNU0GLOBAL]=*prynu;

      // need to limit Ynu0 to table in case changed
      if(EOSextra[YNU0GLOBAL]>0.999*lineartablelimits[primarytable][YNUEOS][1]) EOSextra[YNU0GLOBAL] = 0.999*lineartablelimits[primarytable][YNUEOS][1];
      if(EOSextra[YNU0GLOBAL]<1.001*lineartablelimits[primarytable][YNUEOS][0]) EOSextra[YNU0GLOBAL] = 1.001*lineartablelimits[primarytable][YNUEOS][0];
    
    }



  }

