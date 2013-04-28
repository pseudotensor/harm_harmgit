
#include "decs.h"
#include "reconstructeno.h"  

// monoindicator[0]: -1,0,1 for rough, ambiguous, and monotonic
// monoindicator[1]: whether set cell's left value (or central value, for a2c/c2a reconstruction)
// monoindicator[2]: whether set cell's right value
// yout[0][i] is the left interface value for c2e (centered value for a2c & c2a) for grid cell i
// yout[1][i] is the right interface value for c2e for grid cell i, does not make sense for a2c & c2a
void compute_monotonicity_line(
                               int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                               int *minorder, int *maxorder, int *shift,   
                               FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                               FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM] )
{ 

  int i,jj; 

#if( DOMONOINTERP == JMONOINTERP ) 
  extern void compute_jmonotonicity_line(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);
#elif( DOMONOINTERP == SMONOINTERP ) 
#include "interpline.smono.h"
#elif( DOMONOINTERP == NOMONOINTERP )
#else
#error Unknown monointerp type
#endif


  /////////////////
  //
  // initialize monoindicator so can assume 0 if not using mono rather than leaving as "nan" and having to check if 0 or not
  //
  /////////////////

  for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
    monoindicator[MONOINDTYPE][i]=0;
    monoindicator[MONOLEFTSET][i]=0; // or center
    monoindicator[MONORIGHTSET][i]=0;
    yout[0][i]=yout[1][i]=0;
#if(MERGEDC2EA2CMETHOD)
    for(jj=0;jj<MAXSPACEORDER;jj++) youtpolycoef[jj][i]=0;
#endif
  }


  /////////////////
  //
  // compute monotonicity line for JMONO, SMONO, or NOMONO methods
  //
  /////////////////


#if( DOMONOINTERP == JMONOINTERP ) 
  compute_jmonotonicity_line(
                             recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                             minorder, maxorder, shift,   
                             shockindicator, df, 
                             monoindicator, yin,  yout, youtpolycoef );
#elif( DOMONOINTERP == SMONOINTERP ) 
  compute_smonotonicity_line(
                             recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                             minorder, maxorder, shift,   
                             shockindicator, df, 
                             monoindicator, yin,  yout, youtpolycoef );
#elif( DOMONOINTERP == NOMONOINTERP )
  // no monointerp used; so set the monotonicity indicators to zero
  // nothing to do since already initialized indicators
#else
#error Unknown monointerp type
#endif

}




void compute_monotonicity_line_indicatoronly(
                                             int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                             int *minorder, int *maxorder, int *shift,   
                                             FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                             FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM] )
{ 

  int i; 

#if( DOMONOINTERP == JMONOINTERP ) 
#elif( DOMONOINTERP == SMONOINTERP ) 
#include "interpline.smono.h"
#elif( DOMONOINTERP == NOMONOINTERP )
#else
#error Unknown monointerp type
#endif


  /////////////////
  //
  // initialize monoindicator so can assume 0 if not using mono rather than leaving as "nan" and having to check if 0 or not
  //
  /////////////////

  for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
    monoindicator[MONOINDTYPE][i]=0;
    monoindicator[MONOLEFTSET][i]=0; // or center
    monoindicator[MONORIGHTSET][i]=0;
  }

#if( DOMONOINTERP == JMONOINTERP ) 

  dualfprintf(fail_file,"JMONO not setup for split indicator and value\n");
  myexit(1756582);

#elif( DOMONOINTERP == SMONOINTERP ) 
  compute_smonotonicity_line_split(1,0,
                                   recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                                   minorder, maxorder, shift,   
                                   shockindicator, df, 
                                   monoindicator, yin,  yout, youtpolycoef );
#elif( DOMONOINTERP == NOMONOINTERP )
  // no monointerp used; so set the monotonicity indicators to zero
  // do nothing since already set to 0
#else
#error Unknown monointerp type
#endif

}




void compute_monotonicity_line_valueonly(
                                         int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                         int *minorder, int *maxorder, int *shift,   
                                         FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                         FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM] )
{ 

  int i,jj; 

#if( DOMONOINTERP == JMONOINTERP ) 
#elif( DOMONOINTERP == SMONOINTERP ) 
#include "interpline.smono.h"
#elif( DOMONOINTERP == NOMONOINTERP )
#else
#error Unknown monointerp type
#endif

  // initialize
  for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
    yout[0][i]=yout[1][i]=0;
#if(MERGEDC2EA2CMETHOD)
    for(jj=0;jj<MAXSPACEORDER;jj++) youtpolycoef[jj][i]=0;
#endif
  }


#if( DOMONOINTERP == JMONOINTERP ) 

  dualfprintf(fail_file,"JMONO not setup for split indicator and value\n");
  myexit(1756583);

#elif( DOMONOINTERP == SMONOINTERP ) 
  compute_smonotonicity_line_split(0,1,
                                   recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                                   minorder, maxorder, shift,   
                                   shockindicator, df, 
                                   monoindicator, yin,  yout, youtpolycoef );
#elif( DOMONOINTERP == NOMONOINTERP )
  // nothing to do
#else
#error Unknown monointerp type
#endif

}






// multi-pl version
void compute_monotonicity_line_multipl(int stencilvarisnull, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])
{
  void compute_monotonicity_line_indicatoronly(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);
  void compute_monotonicity_line_valueonly(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM],  FTYPE (*monoindicator)[NBIGM] , FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);
  void monotonicity_equalize_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]);
  void monotonicity_minimizeoneboss_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]);
  void monotonicity_minimizeall_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]);
  void monotonicity_minimize_massenergymomentum_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM]);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  int i;
  FTYPE *monostencilvar;


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


#if(0)
  // DEBUG:
  // old unsplit method that still works but isn't used now that want to minimize over weights for all quantities
  NUMPRIMLOOP(pliter,pl) compute_monotonicity_line(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
  return;
#endif


#if(0)
  // DEBUG:
  // assume always call compute_monotonicity_line_indicatoronly() that initializes monoindicator[]
  /////////////////
  //
  // initialize monoindicator so can assume 0 if not using mono rather than leaving as "nan" and having to check if 0 or not
  //
  /////////////////

  NUMPRIMLOOP(pliter,pl) {
    for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
      monoindicator[pl][MONOINDTYPE][i]=0;
      monoindicator[pl][MONOLEFTSET][i]=0; // or center
      monoindicator[pl][MONORIGHTSET][i]=0;
    }
  }
#endif


  /////////////////
  //
  // Choose reference pl in case needed
  //
  /////////////////


  // get starting (reference) pl
  plstart_set(whichquantity,dir,recontype,&plstart);


  /////////////////
  //
  // set mono weights and set values
  //
  /////////////////


  // only operate if a2c or c2a
  if(  MULTIPLTYPE == ENERGY_IS_ALL_WEIGHTS ) {
      
    ////////////////////
    // set pl as plstart
    pl=plstart;

    //////////////////////
    // set stencil variable
    if(stencilvarisnull) monostencilvar=yin[pl][0];
    else monostencilvar=ystencilvar[pl][0];
   
    //////////////////////////////////////
    //
    // get "weight" for SINGLE "pl" for mono
    compute_monotonicity_line_indicatoronly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], ystencilvar[pl][0], yout[pl], youtpolycoef[pl]);

    /////////////////////////////////
    // 
    // now copy over mono weights across "pl"
    monotonicity_equalize_weights(plstart, MULTIPLTYPE,whichquantity,dir,WEIGHT_CALC,recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef);

    /////////////////////////////////////////////////
    //
    // now compute MONO reconstructions given weights for all "pl"
    NUMPRIMLOOP(pliter,pl){
      compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);
    }

  }// end if energy controls
  else if(MULTIPLTYPE == MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS ){

    // get all indicators
    NUMPRIMLOOP(pliter,pl){
      if(stencilvarisnull) monostencilvar=yin[pl][0];
      else monostencilvar=ystencilvar[pl][0];
      compute_monotonicity_line_indicatoronly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], ystencilvar[pl][0], yout[pl], youtpolycoef[pl]);

    }

    // minimize massenergymomentum weights
    monotonicity_minimize_massenergymomentum_weights(plstart, MULTIPLTYPE,whichquantity,dir,WEIGHT_CALC,recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef);

    //////////////////
    //
    // compute value for given weights
    //
    //////////////////
    NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);


  }
  else{


    // any other methods
    NUMPRIMLOOP(pliter,pl){
      if(stencilvarisnull) monostencilvar=yin[pl][0];
      else monostencilvar=ystencilvar[pl][0];
      compute_monotonicity_line_indicatoronly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], ystencilvar[pl][0], yout[pl], youtpolycoef[pl]);

    }

    if(MULTIPLTYPE == ENERGY_CONTROLS_ALL_WEIGHTS){
      monotonicity_minimizeoneboss_weights(plstart, MULTIPLTYPE,whichquantity,dir,WEIGHT_CALC,recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef);
    }
    else if(MULTIPLTYPE == MINIMIZE_ALL_WEIGHTS){
      monotonicity_minimizeall_weights(plstart, MULTIPLTYPE,whichquantity,dir,WEIGHT_CALC,recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef);
    }

    //////////////////
    //
    // compute value for given weights
    //
    //////////////////
    NUMPRIMLOOP(pliter,pl) compute_monotonicity_line_valueonly(recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df[pl], monoindicator[pl], yin[pl][0], yout[pl], youtpolycoef[pl]);

  }// end else if not energy is all


}




// equalize weights
void monotonicity_equalize_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])

{ 
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int pl,pliter;
  int mypl;
  int i;



  if(do_weight_or_recon != WEIGHT_CALC) return; // nothing to do


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  NUMPRIMLOOP(pliter,pl){
    // copy plstart values into all other values
    for(i=bs;i<=be;i++){
      monoindicator[pl][MONOINDTYPE][i]=monoindicator[plstart][MONOINDTYPE][i];
      monoindicator[pl][MONOLEFTSET][i]=monoindicator[plstart][MONOLEFTSET][i];
      monoindicator[pl][MONORIGHTSET][i]=monoindicator[plstart][MONORIGHTSET][i];
    }
  }

 

}







// minimize weights with dominate term used for all + own-pl term as well
void monotonicity_minimizeoneboss_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])

{ 
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int pl,pliter;
  int mypl;
  int i;


  if(do_weight_or_recon != WEIGHT_CALC) return; // nothing to do


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  // copy plstart values into all other values
  for(i=bs;i<=be;i++){
    NUMPRIMLOOP(pliter,pl){
      monoindicator[pl][MONOINDTYPE][i]=MIN(monoindicator[pl][MONOINDTYPE][i],monoindicator[plstart][MONOINDTYPE][i]);
      monoindicator[pl][MONOLEFTSET][i]=MIN(monoindicator[pl][MONOLEFTSET][i],monoindicator[plstart][MONOLEFTSET][i]);
      monoindicator[pl][MONORIGHTSET][i]=MIN(monoindicator[pl][MONORIGHTSET][i],monoindicator[plstart][MONORIGHTSET][i]);
    }
  }

 

}


// minimize weights with dominate term used for all + own-pl term as well
void monotonicity_minimize_massenergymomentum_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])

{ 
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int pl,pliter;
  int mypl;
  int i;
  FTYPE monoall[NUMMONOINDICATORS];
  int odir1,odir2;



  odir1=dir%3+1;
  odir2=(dir+1)%3+1;


  if(do_weight_or_recon != WEIGHT_CALC) return; // nothing to do


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);



  // determine some initial consistent pl
  NUMPRIMLOOP(pliter,pl){
    if(VELTERMSMINIMIZE(pl) || ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl) || PRESSUREMINIMIZE(pl)){// then minimize across all these
      plstart=pl;
      break;
    }
  }


  for(i=bs;i<=be;i++){
    // CHANGINGMARK: This is not quite consistent with reconstructeno.c where we use "unoptimized weights" to determine which quantity dominates stencil choice

    NUMPRIMLOOP(pliter,pl){
      if(VELTERMSMINIMIZE(pl) || ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl) || PRESSUREMINIMIZE(pl)){// then minimize across all these
        monoindicator[pl][MONOINDTYPE][i]=MIN(monoindicator[pl][MONOINDTYPE][i],monoindicator[plstart][MONOINDTYPE][i]);
        monoindicator[pl][MONOLEFTSET][i]=MIN(monoindicator[pl][MONOLEFTSET][i],monoindicator[plstart][MONOLEFTSET][i]);
        monoindicator[pl][MONORIGHTSET][i]=MIN(monoindicator[pl][MONORIGHTSET][i],monoindicator[plstart][MONORIGHTSET][i]);
        // get current global minimized mono weights
        monoall[MONOINDTYPE]=monoindicator[pl][MONOINDTYPE][i];
        monoall[MONOLEFTSET]=monoindicator[pl][MONOLEFTSET][i];
        monoall[MONORIGHTSET]=monoindicator[pl][MONORIGHTSET][i];
      }
    }

    if(emffixedstencil){
      // if not splitting MA and EM, then assume only EMF is treated specially if point method
      NUMPRIMLOOP(pliter,pl){
        if(EMFTERMS(pl)){// copy equal weights
          monoindicator[pl][MONOINDTYPE][i]=1.0;
          monoindicator[pl][MONOLEFTSET][i]=1.0;
          monoindicator[pl][MONORIGHTSET][i]=1.0;
        }
      }
    }

    NUMPRIMLOOP(pliter,pl){
      if(ALLOTHERSMINIMIZE(pl)){// then minimize across all these
        monoindicator[pl][MONOINDTYPE][i]=MIN(monoindicator[pl][MONOINDTYPE][i],monoall[MONOINDTYPE]);
        monoindicator[pl][MONOLEFTSET][i]=MIN(monoindicator[pl][MONOLEFTSET][i],monoall[MONOLEFTSET]);
        monoindicator[pl][MONORIGHTSET][i]=MIN(monoindicator[pl][MONORIGHTSET][i],monoall[MONORIGHTSET]);
      }
    }

    // rest of weights are independent

  }

 

}


// minimize weights over all pl (not recommended since very restrictive)
void monotonicity_minimizeall_weights(int plstart, int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])

{ 
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int pl,pliter;
  int mypl;
  int i;
  FTYPE monoall[NUMMONOINDICATORS];


  if(do_weight_or_recon != WEIGHT_CALC) return; // nothing to do



  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);




  // define single reference and use for all
  for(i=bs;i<=be;i++){
    // use plstart as starting position for simplicity
    monoall[MONOINDTYPE]=monoindicator[plstart][MONOINDTYPE][i];
    monoall[MONOLEFTSET]=monoindicator[plstart][MONOLEFTSET][i];
    monoall[MONORIGHTSET]=monoindicator[plstart][MONORIGHTSET][i];
    NUMPRIMLOOP(pliter,pl){
      monoall[MONOINDTYPE]=MIN(monoall[MONOINDTYPE],monoindicator[pl][MONOINDTYPE][i]);
      monoall[MONOLEFTSET]=MIN(monoall[MONOLEFTSET],monoindicator[pl][MONOLEFTSET][i]);
      monoall[MONORIGHTSET]=MIN(monoall[MONORIGHTSET],monoindicator[pl][MONORIGHTSET][i]);
    }
    NUMPRIMLOOP(pliter,pl){
      monoindicator[pl][MONOINDTYPE][i]=monoall[MONOINDTYPE];
      monoindicator[pl][MONOLEFTSET][i]=monoall[MONOLEFTSET];
      monoindicator[pl][MONORIGHTSET][i]=monoall[MONORIGHTSET];
    }
  }

 

}
