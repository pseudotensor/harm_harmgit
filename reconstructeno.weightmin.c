//#include "decs.h"
//#include "reconstructeno.h"
//#include "reconstructeno_static.h"

// this file is included in reconstructeno.c


////////////////////////////////////////
//
// Functions for minimization routines
//
////////////////////////////////////////

// these were in slope_lim_linetype
static void minimize_weno5_weights( int preforder, int usecurrentlowerorder, FTYPE reductionfactoroffset, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out );
//static  void minimize_weno5_weights( int preforder, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out );
static void minimize_weno5_weights_old( int preforder, int usecurrentlowerorder, FTYPE reductionorder, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out );

static void print_weno5_weights( int preforder, weno_weights_t *weno_weights);

static void print_weno5_weights_allpl( int preforder, int whichquantity, weno_weights_t (*weno_weights)[NBIGM]);


static  void copy_weno5_weights( int preforder, weno_weights_t *weno_weights_source, weno_weights_t *weno_weights_destination );
static  void reset_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set );
static void reset_equal_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set );
static  void compute_lower_order_fraction_weno5_weights( int preforder, weno_weights_t *weno_weights );
static  void rescale_weno5_weights( int preforder, FTYPE rescale_factor, weno_weights_t *weights_array_in, weno_weights_t *weights_array_out );
static void spread_weno5_weights( int preforder, weno_weights_t *weno_weights_in, weno_weights_t *weno_weights_out);
static void interpolate_stiffness(int preforder, FTYPE *indicator, weno_weights_t *weight_stiff, weno_weights_t *weight_notstiff);
static void interpolate_2weightsto1weight(int preforder, int dir, FTYPE (*ystencilvar)[2][NBIGM], weno_weights_t *weightsg1, weno_weights_t *weightsg2, weno_weights_t *weights_out);

static void interpolate_diffweighted_weights(int preforder, int whichquantity, int dir, FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*normalystencilvar)[NBIGM], weno_weights_t (*stencil_weights_array_pl)[NBIGM], struct of_trueijkp *trueijkp);



/////////////////////////////////////
//
// Storage for weights
//
/////////////////////////////////////









// get all "pl" interpolations
void compute_multipl_weno(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, weno_weights_t (*stencil_weights_array_allpl)[NBIGM], int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp)
{
  // for NUMPRIMLOOP:
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  // others:
  int odir1,odir2;

  // OPENMPMARK: Can't leave as static
  FTYPE a_normalystencilvar[NPR][NBIGM];
  FTYPE (*normalystencilvar)[NBIGM];
  weno_weights_t a_weights_array_minimal[NBIGM];  //for storing the minimums of the weights for a one-dimensional array of values; to feed to stencil reduction
  weno_weights_t *weights_array_minimal = &(a_weights_array_minimal[NBIGBND]);

  //weno_weights_t weno_weights_backup[NPR2INTERP][NBIGM]; //for storing the weights for all quantities


  weno_weights_t a_stencil_weights_array_storemin[NBIGM];
  weno_weights_t *stencil_weights_array_storemin = &(a_stencil_weights_array_storemin[NBIGBND]);

  weno_weights_t a_stencil_weights_array_pl[NPR][NBIGM];
  weno_weights_t (*stencil_weights_array_pl)[NBIGM];


#define CODE_FOR_MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD 0

#if(CODE_FOR_MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD)
  // old groups for MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD
  weno_weights_t a_stencil_weights_array_group1[NBIGM];
  weno_weights_t *stencil_weights_array_group1 = &(a_stencil_weights_array_group1[NBIGBND]);
  weno_weights_t a_stencil_weights_array_group2[NBIGM];
  weno_weights_t *stencil_weights_array_group2 = &(a_stencil_weights_array_group2[NBIGBND]);
  weno_weights_t a_stencil_weights_array_group3[NBIGM];
  weno_weights_t *stencil_weights_array_group3 = &(a_stencil_weights_array_group3[NBIGBND]);
  weno_weights_t a_stencil_weights_array_group4[NBIGM];
  weno_weights_t *stencil_weights_array_group4 = &(a_stencil_weights_array_group4[NBIGBND]);
#endif



  /////////////////
  //
  // Setup pointers
  //
  /////////////////

  stencil_weights_array_pl = (weno_weights_t (*)[NBIGM]) (&(a_stencil_weights_array_pl[0][NBIGBND]));

  normalystencilvar = (FTYPE (*)[NBIGM]) (&(a_normalystencilvar[0][NBIGBND]));

  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  ///////////////
  //
  //split version, compute weights in the 1st call and do the reconstruction in 2nd call
  //
  ///////////////

  reset_weno5_weights( preforder, weights_array_minimal );





  if( MULTIPLTYPE == MINIMIZE_ALL_WEIGHTS ) { //a2c on conserved


    NUMPRIMLOOP(pliter,pl){
      pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl],youtpolycoef[pl],trueijkp);
      minimize_weno5_weights_old( preforder, 1, 0.0, stencil_weights_array_allpl[pl], weights_array_minimal, weights_array_minimal );
    }


    NUMPRIMLOOP(pliter,pl) {
      //put the minimal weights back into the weno storage so that weno uses those weights; 1 stands for copy
      copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_allpl[pl] );
      pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
    }
    //SUPERSASMARK : insert a2c/c2a limiting code here


  }
  else if(MULTIPLTYPE == MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS){




    // GODMARK: Note that weight calculation itself has error that is really third order for WENO5 since no third derivatives taken.  This means can't trust weight to some degree of variance.  Hence, start with 1/3 and only go below 1/3 for WENO5 weights if significantly different than 1/3 so don't accumulate lower order fractions that (according to Sasha) generically give only forth order accuracy for generic smooth flows due to the error in the weight calculation.  Should really estimate the weight error for smooth flows (Sasha thinks 1/3 -> 1/30 is range, but that sounds crazy, so I used 1/3 -> 1.3/3 a 30% range to allow the weights to deviate from 1/3 even for smooth flows)
    //


    odir1=dir%3+1;
    odir2=(dir+1)%3+1;





    /////////////////
    //
    // Treat velocity-along-dir terms
    //
    /////////////////

    // Assume weights higher-order to start with
    reset_equal_weno5_weights( preforder, weights_array_minimal);

    NUMPRIMLOOP(pliter,pl){
      if(VELTERMSMINIMIZE(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array_allpl[pl] ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_pl[pl], weights_array_minimal, weights_array_minimal );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, weights_array_minimal,weights_array_minimal);


    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

    // copy to all terms
    NUMPRIMLOOP(pliter,pl){
      if(VELTERMSMINIMIZE(pl)){// then minimize across all these
        copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
      }
    }




    /////////////////
    //
    // treat flux of orthogonal velocity #1
    //
    /////////////////

    // Assume weights higher-order to start with
    reset_equal_weno5_weights( preforder, weights_array_minimal);

    NUMPRIMLOOP(pliter,pl){
      if(ORTHOVEL1TERMSMINIMIZE(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_pl[pl], weights_array_minimal, weights_array_minimal );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, weights_array_minimal,weights_array_minimal);


    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

    // copy to all terms
    NUMPRIMLOOP(pliter,pl){
      if(ORTHOVEL1TERMSMINIMIZE(pl)){// then minimize across all these
        copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
      }
    }



    /////////////////
    //
    // treat flux of orthogonal velocity #2
    //
    /////////////////

    // Assume weights higher-order to start with
    reset_equal_weno5_weights( preforder, weights_array_minimal);

    NUMPRIMLOOP(pliter,pl){
      if(ORTHOVEL2TERMSMINIMIZE(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_pl[pl], weights_array_minimal, weights_array_minimal );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, weights_array_minimal,weights_array_minimal);


    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

    // copy to all terms
    NUMPRIMLOOP(pliter,pl){
      if(ORTHOVEL2TERMSMINIMIZE(pl)){// then minimize across all these
        copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
      }
    }




    /////////////////
    //
    // Treat pressure term or any other terms that are part of main EOMs that determine inversion success
    //
    /////////////////


    // Assume weights higher-order to start with
    reset_equal_weno5_weights( preforder, weights_array_minimal);

    NUMPRIMLOOP(pliter,pl){
      if(PRESSUREMINIMIZE(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_pl[pl], weights_array_minimal, weights_array_minimal );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, weights_array_minimal,weights_array_minimal);

#if(0)
    // linearly interpolate between the independent weight and the minimal weight
    // No need to compute 0,10.0 version of minimal weight -- just switch between totally independent and fully coupled using stiffness parameter
    // stencil_weights_array_storemin vs. stencil_weights_array
    //   interpolate_stiffness(preforder, stiffindicator, stencil_weights_array_storemin, stencil_weights_array);
#endif


    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

    // copy to all terms
    NUMPRIMLOOP(pliter,pl){
      if(PRESSUREMINIMIZE(pl)){// then minimize across all these
        copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
      }
    }




    if(emffixedstencil){
      // if not splitting MA and EM, then assume only EMF is treated specially if point method

      // Assume weights equal
      reset_equal_weno5_weights( preforder, weights_array_minimal);
      
      NUMPRIMLOOP(pliter,pl){
        if(EMFTERMS(pl)){// copy equal weights
          copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
        }
      }
      
    }


    /////////////////
    //
    // NOW treat rest of quantities (avoids EMF terms if emffixedstencil==0 when splitmaem==1)
    //
    /////////////////

    // Assume weights higher-order to start with
    reset_equal_weno5_weights( preforder, weights_array_minimal);

    NUMPRIMLOOP(pliter,pl){
      if(ALLOTHERSMINIMIZE(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_pl[pl], weights_array_minimal, weights_array_minimal );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, weights_array_minimal,weights_array_minimal);

    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, weights_array_minimal );

    // copy to all terms
    NUMPRIMLOOP(pliter,pl){
      if(ALLOTHERSMINIMIZE(pl)){// then minimize across all these
        copy_weno5_weights( preforder, weights_array_minimal, stencil_weights_array_pl[pl]);
      }
    }


    //      print_weno5_weights_allpl( preforder, whichquantity, stencil_weights_array_pl); // CHANGINGMARK



    /////////////////
    //
    // NOW determine final weights
    //
    // Each weight for each quantity is set to be an average of all weights as weighted by all quantities' differences with some extra weight given to the quantity in question
    //
    /////////////////

    // CHANGINGMARK:
    interpolate_diffweighted_weights(preforder, whichquantity, dir, ystencilvar, normalystencilvar, stencil_weights_array_pl,trueijkp);


    //      print_weno5_weights_allpl( preforder, whichquantity, stencil_weights_array_pl); // CHANGINGMARK



    //////////////////////////////
    //
    // NOW reconstruct all
    //
    //////////////////////////////

    //put the minimal weights back into the weno storage so that weno uses those weights; 1 stands for copy

    // now compute all of these quantities using single coupled minimized weights
    NUMPRIMLOOP(pliter,pl) {
      //copy_weno5_weights( preforder, stencil_weights_array_pl[pl], stencil_weights_array );
      pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_pl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
    }




  }



#if(CODE_FOR_MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD)
  else if(MULTIPLTYPE == MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS_OLD){ // older version of this idea


    odir1=dir%3+1;
    odir2=(dir+1)%3+1;


#define WHICHPLTOMINIMIZE1(pl) (VELTERMSMINIMIZE(pl))
#define WHICHPLTOMINIMIZE2(pl) (ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl))
#define WHICHPLTOMINIMIZE3(pl) (PRESSUREMINIMIZE(pl))
#define WHICHPLTOMINIMIZE4(pl) (ALLOTHERSMINIMIZE(pl))


    /////////////////
    //
    // Treat velocity-along-dir terms
    //
    /////////////////

    // Assume weights higher-order to start with

    reset_equal_weno5_weights( preforder, stencil_weights_array_group1);


    ///////////////////////////
    //
    // NOW minimize over mass, energy, and momentum along dir (carrying at first weight from above pressure term)
    //
    // GODMARK: Note that weight calculation itself has error that is really third order for WENO5 since no third derivatives taken.  This means can't trust weight to some degree of variance.  Hence, start with 1/3 and only go below 1/3 for WENO5 weights if significantly different than 1/3 so don't accumulate lower order fractions that (according to Sasha) generically give only forth order accuracy for generic smooth flows due to the error in the weight calculation.  Should really estimate the weight error for smooth flows (Sasha thinks 1/3 -> 1/30 is range, but that sounds crazy, so I used 1/3 -> 1.3/3 a 30% range to allow the weights to deviate from 1/3 even for smooth flows)
    //
    //////////////////////////
    NUMPRIMLOOP(pliter,pl){

      if(WHICHPLTOMINIMIZE1(pl)){// then minimize across all these

        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        //   dualfprintf(fail_file,"pl=%d\n",pl);
        //   print_weno5_weights( preforder, stencil_weights_array ); // CHANGINGMARK

        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_allpl[pl], stencil_weights_array_group1, stencil_weights_array_group1 );
      }
    }

    // spread weno weights
    spread_weno5_weights(preforder, stencil_weights_array_group1,stencil_weights_array_group1);


    // finally set the lower_order_fraction of the weights to match that of the sum of the weno5 weights
    compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array_group1 );





    /////////////////
    //
    // Treat pressure term or any other terms that are part of main EOMs that determine inversion success
    //
    // But override independent weight if in stiff regime between velocity and pressure terms
    //
    /////////////////

    // first assume equal weights
    reset_equal_weno5_weights( preforder, stencil_weights_array_group3);

    // now get weights for this group
    NUMPRIMLOOP(pliter,pl){
      if(WHICHPLTOMINIMIZE3(pl)){
        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

        // minimize weights
        minimize_weno5_weights( preforder, 1, WEIGHTFACTORMINIMIZE, stencil_weights_array_allpl[pl], stencil_weights_array_group3, stencil_weights_array_group3);


        // spread these independent weights before interpolating between weights using stiffness
        //   spread_weno5_weights(preforder, stencil_weights_array,stencil_weights_array); // CHANGINGMARK

#if(1)
        // linearly interpolate between the independent weight and the minimal weight
        // No need to compute 0,10.0 version of minimal weight -- just switch between totally independent and fully coupled using stiffness parameter
        // stencil_weights_array_storemin vs. stencil_weights_array
        interpolate_stiffness(preforder, stiffindicator, stencil_weights_array_storemin, stencil_weights_array);
#endif
      }
    }

    // finally recompute lower order fraction as required
    compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array_group3 );




#if(SPLITPRESSURETERMINFLUX)
    /////////////////
    //
    // NOW determine final common weight for velocity-along-dir and pressure
    //
    /////////////////

    interpolate_2weightsto1weight(preforder, dir, ystencilvar, stencil_weights_array_group1, stencil_weights_array_group3, stencil_weights_array_storemin);
    //      copy_weno5_weights( preforder, stencil_weights_array_group1, stencil_weights_array_storemin);

    //      print_weno5_weights( preforder, stencil_weights_array_storemin); // CHANGINGMARK
#else
    copy_weno5_weights( preforder, stencil_weights_array_group1, stencil_weights_array_storemin);
#endif


    //////////////////////////////
    //
    // NOW reconstruct group1
    //
    //////////////////////////////

    // now compute all of these quantities using single coupled minimized weights
    NUMPRIMLOOP(pliter,pl) {
      //put the minimal weights back into the weno storage so that weno uses those weights; 1 stands for copy
      //      copy_weno5_weights( preforder, stencil_weights_array_group1, stencil_weights_array );
      copy_weno5_weights( preforder, stencil_weights_array_storemin, stencil_weights_array_allpl[pl] );

      if(WHICHPLTOMINIMIZE1(pl)){
        pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
      }
    }


    ///////////////////////////////////
    //
    // NOW reconstruct group3
    //
    ///////////////////////////////////

    // reconstruct group3
    NUMPRIMLOOP(pliter,pl){
      //put the minimal weights back into the weno storage so that weno uses those weights
      //      copy_weno5_weights( preforder, stencil_weights_array_group3, stencil_weights_array );
      copy_weno5_weights( preforder, stencil_weights_array_storemin, stencil_weights_array_allpl[pl] );

      if(WHICHPLTOMINIMIZE3(pl)){
        pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
      }
    }


    /////////////////
    //
    // NOW treat group2
    //
    // treat flux of orthogonal velocity terms starting with above minimized weights but further minimizing
    //
    /////////////////

    NUMPRIMLOOP(pliter,pl) {
      if(WHICHPLTOMINIMIZE2(pl)){
        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
        //otherwise, the common weights have already been precomputed, so modify the current weights to be smallest of the two:  the energy weights and twice the current weights
        minimize_weno5_weights( preforder, 0, 9.0, stencil_weights_array_allpl[pl], stencil_weights_array_storemin, stencil_weights_array_allpl[pl] );

        //spread_weno5_weights(preforder, stencil_weights_array,stencil_weights_array);

        //synchronise the lower_order_fraction with the normalization of the resulting weights
        compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array_allpl[pl] );
        // otherwise do not recompute weights and just use single consistent weight

        pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
      }
    }


    /////////////////
    //
    // NOW treat all other terms
    //
    /////////////////

    NUMPRIMLOOP(pliter,pl) {
      if(WHICHPLTOMINIMIZE4(pl)){
        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
        //otherwise, the common weights have already been precomputed, so modify the current weights to be smallest of the two:  the energy weights and twice the current weights
        minimize_weno5_weights( preforder, 0, 9.0, stencil_weights_array_allpl[pl], stencil_weights_array_storemin, stencil_weights_array_allpl[pl] );

        //   spread_weno5_weights(preforder, stencil_weights_array,stencil_weights_array);

        //synchronise the lower_order_fraction with the normalization of the resulting weights
        compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array_allpl[pl] );
        // otherwise do not recompute weights and just use single consistent weight

        pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
      }
    }


  }
#endif


  else if(  MULTIPLTYPE == ENERGY_CONTROLS_ALL_WEIGHTS ||   MULTIPLTYPE == ENERGY_IS_ALL_WEIGHTS ) {

    // get starting (reference) pl
    plstart_set(whichquantity,dir,recontype,&plstart);



    // had to change this since Sasha assumed certain ordering of pl's
    pl=plstart;
    pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
    //store the weights for this quantity (energy or momentum depending on the interpolation type) as common weights;
    //they will be used for comparison to weights computed for other quantities
    minimize_weno5_weights_old( preforder, 1, 0.0, stencil_weights_array_allpl[pl], weights_array_minimal, weights_array_minimal );
    pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
    //SUPERSASMARK : insert a2c/c2a limiting code here



    NUMPRIMLOOP(pliter,pl) {

      if(pl==plstart) continue;

      // CHANGINGMARK (makes difference on ki-rh42 at least!!)
      // if(pl==B1) continue;
      // if(pl>=B1) continue;


      //      dualfprintf(fail_file,"RECONSTRUCT: pl=%d\n",pl);

      if(  MULTIPLTYPE == ENERGY_CONTROLS_ALL_WEIGHTS ) {
        pass_1d_line_weno_withweights( whichquantity, dir, WEIGHT_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);
        //otherwise, the common weights have already been precomputed, so modify the current weights to be smallest of the two:  the energy weights and twice the current weights
        minimize_weno5_weights_old( preforder, 0, 9.0, stencil_weights_array_allpl[pl], weights_array_minimal, stencil_weights_array_allpl[pl] );
        //synchronise the lower_order_fraction with the normalization of the resulting weights
        compute_lower_order_fraction_weno5_weights( preforder, stencil_weights_array_allpl[pl] );
      }
      // otherwise do not recompute weights and just use single consistent weight

      pass_1d_line_weno_withweights( whichquantity, dir, RECON_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[pl], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl], youtpolycoef[pl],trueijkp);

      // if(pl==B1){
      //   for( i = bs; i <= be; i++ ) {
      //     dualfprintf(fail_file,"ps=%d pe=%d :: i=%d yin=%21.15g %21.15g :: yout=%21.15g %21.15g\n",ps,pe,i,yin[pl][0][i],yin[pl][1][i],yout[pl][0][i],yout[pl][1][i]);
      //     yout[pl][0][i]=yout[pl][1][i]=0.0;
      //   }
      //     }

    }
    //SUPERSASMARK : insert a2c/c2a limiting code here
  }



  return;
}









//NOTE: preforder should be set correctly; if it is not, the behaviour will be wrong.  
//      preforder is used because some of the weights at the boundary are not defined and 
//      for those weights one cannot use the order in the summing procedures.
//      Basically, the order that comes with the weights is copied over but is ignored
//      for the purposes of the summing, etc. in these functions

static void rescale_weno5_weights( int preforder, FTYPE rescale_factor, weno_weights_t *weights_array_in, weno_weights_t *weights_array_out )
{
  int weight_count, i;
  int order = (preforder + 1) / 2;

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weights_array_out[i].weights[weight_count] = rescale_factor * weights_array_in[i].weights[weight_count];  //rescale the weights
    }
    weights_array_out[i].lower_order_fraction = 0.0;
    weights_array_out[i].order = order;
    weights_array_out[i].len = 0;
  }
}



//resets the unoptimized weights to unity; lower_order_fraction to 0.
static void reset_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set )
{
  int weight_count, i;
  int order = (preforder + 1) / 2;

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weno_weights_to_set[i].weights[weight_count] = 1.0;  //set the weights to unity
    }
    weno_weights_to_set[i].lower_order_fraction = 0.0;
    weno_weights_to_set[i].order = order;
    weno_weights_to_set[i].len = 0;
  }
}


//resets the unoptimized weights to unity; lower_order_fraction to 0.
static void reset_equal_weno5_weights( int preforder, weno_weights_t *weno_weights_to_set )
{
  int weight_count, i;
  int order = (preforder + 1) / 2;

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weno_weights_to_set[i].weights[weight_count] = 1.0/order;
    }
    weno_weights_to_set[i].lower_order_fraction = 0.0;
    weno_weights_to_set[i].order = order;
    weno_weights_to_set[i].len = 0;
  }
}

//sets lower_order_fraction
//copies WENO-5 weights, order, len, and lower-order fraction from source to destination (does not modify WENO-3 weights)
//SASMARK: CHANGED BEHAVIOR: previously, would copy structures as a whole but now, since want to preserve the weno-3 weights and the 
//norm of the weights, only copy the relevant information
static void copy_weno5_weights( int preforder, weno_weights_t *weno_weights_source, weno_weights_t *weno_weights_destination )
{
  int weight_count, i;
  int order = (preforder + 1) / 2;
  FTYPE weights_sum;
  int weight_no;

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    //copy the weights over
    //copy the structures as a whole (as single objects)
    //weno_weights_destination[i] = weno_weights_source[i];

    // also "copy" lower order fraction by setting them
    weights_sum = get_sum_of_elements( order, weno_weights_destination[i].weights );
    weno_weights_destination[i].lower_order_fraction = 1.0 - weights_sum;

    // also copy other things
    weno_weights_destination[i].order = weno_weights_source[i].order;
    weno_weights_destination[i].len = weno_weights_source[i].len;

    for( weight_no = 0; weight_no < weno_weights_source[i].len; weight_no++ ) {
      weno_weights_destination[i].weights[weight_no] = weno_weights_source[i].weights[weight_no];
    }

  }

}

static void print_weno5_weights( int preforder, weno_weights_t *weno_weights)
{
  int weight_count, i;
  int order = (preforder + 1) / 2;

  //simply copy the weights over
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d lof=%21.15g ",nstep,steppart,i,weno_weights[i].lower_order_fraction);
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      dualfprintf(fail_file,"w[%d]=%21.15g ",weight_count,weno_weights[i].weights[weight_count]);
    }
    dualfprintf(fail_file,"\n");
  }
}

static void print_weno5_weights_allpl( int preforder, int whichquantity, weno_weights_t (*weno_weights)[NBIGM])
{
  int weight_count, i;
  int order = (preforder + 1) / 2;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  int plother,pllocalother;


  // setup NUMPRIMLOOP
  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  //simply copy the weights over
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d ",nstep,steppart,i);
    dualfprintf(fail_file,"\n");
    NUMPRIMLOOP(pliter,pl){
      dualfprintf(fail_file,"lof[pl=%d]=%21.15g ",pl,weno_weights[pl][i].lower_order_fraction);
      for( weight_count = 0; weight_count < order; weight_count++ ) {
        dualfprintf(fail_file,"w[pl=%d][%d]=%21.15g ",pl,weight_count,weno_weights[pl][i].weights[weight_count]);
      }
      dualfprintf(fail_file,"\n");
    }
    dualfprintf(fail_file,"\n");
  }
}

//computes the sum of the weno5 weights and set the lower order fraction to 1 minus this value
static void compute_lower_order_fraction_weno5_weights( int preforder, weno_weights_t *weno_weights )
{
  FTYPE weights_sum;
  int order = (preforder + 1) / 2;
  int weight_count, i;

  //simply copy the weights over
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    //SASMARK order hard-cored
    weights_sum = get_sum_of_elements( order, weno_weights[i].weights );  //copy the structures as a whole (as single objects)
    weno_weights[i].lower_order_fraction = 1.0 - weights_sum;
  }
}

// use parameter to interpolate weight from one to another weight
// assume indicator: 0=use original weight, 1=use storemin
// assumes notstiff is final output and that no spreading done inside here
static void interpolate_stiffness(int preforder, FTYPE *indicator, weno_weights_t *weight_stiff, weno_weights_t *weight_notstiff)
{
  int order = (preforder + 1) / 2;
  int weight_count, i;

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    // DEBUG:
    //dualfprintf(fail_file,"nstep=%ld steppart=%d :: i=%d stiff=%21.15g\n",nstep,steppart,i,indicator[i]);

    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weight_notstiff[i].weights[weight_count] = (weight_stiff[i].weights[weight_count])*indicator[i] + (weight_notstiff[i].weights[weight_count])*(1.0-indicator[i]);
    }
  }
}



// determines how much to use other quantities to change one quantity in question
#define WEIGHTCOEF (0.5)

// obtain weight as some partially averaged weight
// everything is point-based so never need to worry about duplicated memory
static void interpolate_diffweighted_weights(int preforder, int whichquantity, int dir, FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*normalystencilvar)[NBIGM], weno_weights_t (*stencil_weights_array_pl)[NBIGM], struct of_trueijkp *trueijkp)
{
  int di,dj,dk;
  int i, weight_count;
  int order = (preforder + 1) / 2;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  int stencilpl;
  int plother,pllocalother;
  struct of_geom geom;
  int theotherdir;
  FTYPE myterm,otherterm,dFsum,maxterm;
  int maxpl;
  FTYPE ftemp;
  int minversion;
  int odir1,odir2;
  FTYPE minflat;
  FTYPE ftemp1,ftemp2;


  // setup orthogonal directions to dir-direction
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;

  // setup NUMPRIMLOOP
  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);

  // iterglobal should really be passed
  di=(trueijkp->iter==1);
  dj=(trueijkp->iter==2);
  dk=(trueijkp->iter==3);


  ///////////////////
  //
  // get dF's
  //
  ///////////////////
  NUMPRIMLOOP(pliter,pl){

    // Note FLUXSPLITPMA(dir) pressure term handled normally
    if(VELTERMSMINIMIZE(pl)) stencilpl=UU+dir; // use single quantity for all velocity terms
    else stencilpl=pl;

    // obtain (\delta F)^4 (consistent with p=2)
    for( i = -NBIGBND+1; i < NBIG + NBIGBND-1; i++ ) { // avoid going out of bounds of memory



#if(1) // NEW method


      // finally using fourth power of derivatives, like WENO weight itself
      normalystencilvar[pl][i] = MAX(MAX(fabs(stencil_weights_array_pl[stencilpl][i].smoothness_indicators[0]),fabs(stencil_weights_array_pl[stencilpl][i].smoothness_indicators[1])) , fabs(stencil_weights_array_pl[stencilpl][i].smoothness_indicators[2]));
      normalystencilvar[pl][i] *= normalystencilvar[pl][i]; // square it
      //      ftemp2=normalystencilvar[pl][i];


#else // OLD method


      // centered difference (can miss local jumps)
      normalystencilvar[pl][i] =  0.5*fabs(ystencilvar[stencilpl][0][i+1] - ystencilvar[stencilpl][0][i-1]);
      normalystencilvar[pl][i] *= normalystencilvar[pl][i]; // square it
      normalystencilvar[pl][i] *= normalystencilvar[pl][i]; // square it
      //      ftemp1=normalystencilvar[pl][i];

      // MAX of left and right differences so more robust (NO, less robust!)
      //      normalystencilvar[pl][i] =  MAX(MAX(MAX(fabs(ystencilvar[stencilpl][0][i] - ystencilvar[stencilpl][0][i-1]),fabs(ystencilvar[stencilpl][0][i+1] - ystencilvar[stencilpl][0][i])),fabs(ystencilvar[stencilpl][0][i+2] - ystencilvar[stencilpl][0][i+1])),fabs(ystencilvar[stencilpl][0][i-1] - ystencilvar[stencilpl][0][i-2]));

#endif

      // DEBUG:
      //      dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d pl=%d orig=%21.15g new=%21.15g\n",nstep,steppart,i,pl,ftemp1,ftemp2);



      // avoid 0 but ensure used in way that compared to SMALL later used won't affect equal weights
      normalystencilvar[pl][i] = MAX(normalystencilvar[pl][i],SMALL/NUMEPSILON);

    }

  }


  ///////////////////
  //
  // Process standard advective terms together
  // All pl's for VELTERMSMINIMIZE all have same dF so don't have to treat specially now
  //
  ///////////////////

  NUMPRIMLOOP(pllocal,pl){
    if(VELTERMSMINIMIZE(pl) || ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl) || PRESSUREMINIMIZE(pl)){
      minversion=1;
    }
    else if(ALLOTHERSMINIMIZE(pl)){
      minversion=2;
    }
    else continue; // assume not minimzing that pl (e.g. point field FLUXRECON stag method with splitmaem==0 has EMF terms at constant weights

    if(emffixedstencil==1 && pl>=B1 && pl<=B3){
      dualfprintf(fail_file,"Shouldn't be here in interpolate_diffweighted_weights(): pl=%d\n",pl);
      myexit(893657476);

    }

    
    for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
      for( weight_count = 0; weight_count < order; weight_count++ ) {


        // get other term and sum term
        otherterm = 0.0;
        maxterm = -BIG; //0.0;  -- SASMARKXXX correction
        myterm = 0.0;
        dFsum = 0.0;
        NUMPRIMLOOP(pllocalother,plother){

          if(
             (minversion==1 && (VELTERMSMINIMIZE(plother) || ORTHOVEL1TERMSMINIMIZE(plother) || ORTHOVEL2TERMSMINIMIZE(plother) || PRESSUREMINIMIZE(plother)))
             ||
             (minversion==2 && (ALLOTHERSMINIMIZE(plother) || (VELTERMSMINIMIZE(plother) || ORTHOVEL1TERMSMINIMIZE(plother) || ORTHOVEL2TERMSMINIMIZE(plother) || PRESSUREMINIMIZE(plother)) ) )
             ){

            // get W_i (\delta Fhat)^2
            ftemp = (stencil_weights_array_pl[plother][i].weights[weight_count])*normalystencilvar[plother][i];

            if(normalystencilvar[plother][i]>=maxterm){
              maxterm=normalystencilvar[plother][i];
              maxpl=plother;
            }

            if(plother!=pl){
              otherterm += WEIGHTCOEF*ftemp;
              dFsum += WEIGHTCOEF*normalystencilvar[plother][i]; // over all
            }
            else{
              // then plother==pl
              myterm += ftemp;
              dFsum += normalystencilvar[plother][i]; // over all
            }
          }
        }
        // now interpolate the weight:
        // W_i = [W_i[old] (\delta Fhat_i)^2 + WEIGHTCOEF\Sum_{all others} W_{other} (\delta Fhat_{other})^2]/[\Sum_{all} (\delta F)^2]
        //
        // stencil_weights_array_pl[pl][i].weights[weight_count] = (myterm + otherterm)/(dFsum+SMALL);

        // stencil_weights_array_pl[pl][i].weights[weight_count]
        //   = MIN(10.0*stencil_weights_array_pl[pl][i].weights[weight_count],stencil_weights_array_pl[maxpl][i].weights[weight_count]);
        if(normalystencilvar[UU+dir][i]>0.5*normalystencilvar[maxpl][i]){
          // quasi-energy is boss or energy is all-weights
          maxpl=UU+dir; // override if close enough
        }

        // maxpl=UU+dir; DEBUG

        stencil_weights_array_pl[pl][i].weights[weight_count] = stencil_weights_array_pl[maxpl][i].weights[weight_count]; // most conservative

      }
    }
  }


  ///////////////////////////////
  //
  // get lower order fractions for all pl (assume len and order don't change)
  //
  ///////////////////////////////
  NUMPRIMLOOP(pliter,pl){
    compute_lower_order_fraction_weno5_weights( preforder,  stencil_weights_array_pl[pl]);
  }


}




// interpolate between 2 weights
// purpose is to avoid reducing order of pressure term when v~0 since couple all v terms together always
// done per-point so if out=inputs then still ok
// requires 1 extra boundary cell
static void interpolate_2weightsto1weight(int preforder, int dir, FTYPE (*ystencilvar)[2][NBIGM], weno_weights_t *weightsg1, weno_weights_t *weightsg2, weno_weights_t *weights_out)
{
  int i, weight_count;
  int order = (preforder + 1) / 2;
  FTYPE dg1,dg2;


  // CHANGINGMARK: need 1 more point for a2c to use this (13 total ghost cells)
#if(MAXBND<13)
  // CHANGINGMARK: Should add check for using WENO5BNDPLUSMIN
  dualfprintf(fail_file,"interpolate_2weightsto1weight() cannot be used with MAXBND<13\n");
  myexit(18672160);
#endif


  // loop over spatial positions
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    // get order,len from either g1 or g2
    weights_out[i].order = weightsg1[i].order;
    weights_out[i].len = weightsg1[i].len;

    // compute dg1 and dg2
    // uses +-1 values! -- so need MAXBND==13
    dg1 = fabs(ystencilvar[UU+dir][0][i+1] - ystencilvar[UU+dir][0][i-1]);
    dg2 = fabs(ystencilvar[FLUXSPLITPMA(dir)][0][i+1] - ystencilvar[FLUXSPLITPMA(dir)][0][i-1]);


    //    dg1 = 0.0;
    //    dg2 = 1.0;

    // DEBUG:
    //    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d dg1=%21.15g dg2=%21.15g\n",nstep,steppart,i,dg1,dg2);


    for( weight_count = 0; weight_count < order; weight_count++ ) {
      // interpolate between g1 and g2 weights for final total weight based upon derivative of stencilvar
      weights_out[i].weights[weight_count] = ((dg1 * weightsg1[i].weights[weight_count]) + (dg2 * weightsg2[i].weights[weight_count]))/(dg1+dg2+SMALL);
    }

  }


  // get lower order fraction
  compute_lower_order_fraction_weno5_weights( preforder, weights_out );


}


//compute the minimum of each weight in the weno_weights_current and weno_weights_minimal, and put these minimums into weno_weights_minimal_out
//the weno_weights_current's weights are multipulied by (1-lower_order_fraction) before comparison.
// Only correct for a2c and c2a GODMARK
// usecurrentlowerorder: 0: use reductionfactor to multiply weight  1: use current lower order fraction to multiply weight (as needed to account for fact that weno5 weight normalized so sum_i w_i = 1 but that later RECON_CALC assumes sum of weights from weno3+weno5 =1)
static void minimize_weno5_weights( int preforder, int usecurrentlowerorder, FTYPE reductionfactoroffset, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out )
{
  int i, weight_count;
  int order = (preforder + 1) / 2;
  FTYPE reductionfactor,finalreductionfactor;
  FTYPE lowerfrac;

  reductionfactor = 1.0+reductionfactoroffset;

  
  //keep track of the minimum value of each weight in weno_weights_minimal
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    //    weno_weights_minimal_out[i].lower_order_fraction = -9;  //lower order fraction is not set

    //    finalreductionfactor = ( 1.0 - weno_weights_current[i].lower_order_fraction )*usecurrentlowerorder + (reductionfactor)*(1-usecurrentlowerorder);
    //    finalreductionfactor = ( 1.0 - weno_weights_current[i].lower_order_fraction )*usecurrentlowerorder + (reductionfactor);
    lowerfrac = ( 1.0 - weno_weights_current[i].lower_order_fraction )*usecurrentlowerorder + 1.0*(1-usecurrentlowerorder);
    finalreductionfactor = lowerfrac*(reductionfactor);

    weno_weights_minimal_out[i].order = weno_weights_current[i].order;
    weno_weights_minimal_out[i].len = weno_weights_current[i].len;
    //only keep track of the maximal lower order fraction but 
    //do not copy the order and len since they do not change from quantity to quantity; however, the weights ratios will change, but we don't use them after this point
    //keep track of the *unoptimized* minimal weights and put them into the destination
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weno_weights_minimal_out[i].weights[weight_count] = MIN( 
                                                              weno_weights_current[i].weights[weight_count] * finalreductionfactor, //rescale the weights by the weno5 fraction
                                                              weno_weights_minimal[i].weights[weight_count] );
    }
  }
}

//compute the minimum of each weight in the weno_weights_current and weno_weights_minimal, and put these minimums into weno_weights_minimal_out
//the weno_weights_current's weights are multipulied by (1-lower_order_fraction) before comparison.
// Only correct for a2c and c2a GODMARK
// if 2nd and 4th arguments are same, then de-emphasize current weights and ignore existing lower order fraction
static void minimize_weno5_weights_old( int preforder, int usecurrentlowerorder, FTYPE reductionfactoroffset, weno_weights_t *weno_weights_current, weno_weights_t *weno_weights_minimal, weno_weights_t *weno_weights_minimal_out )
{
  int i, weight_count;
  int order = (preforder + 1) / 2;
  
  //keep track of the minimum value of each weight in weno_weights_minimal
  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {
    weno_weights_minimal_out[i].lower_order_fraction = -9;  //lower order fraction is not set
    weno_weights_minimal_out[i].order = weno_weights_current[i].order;
    weno_weights_minimal_out[i].len = weno_weights_current[i].len;
    //only keep track of the maximal lower order fraction but 
    //do not copy the order and len since they do not change from quantity to quantity; however, the weights ratios will change, but we don't use them after this point
    //keep track of the *unoptimized* minimal weights and put them into the destination
    for( weight_count = 0; weight_count < order; weight_count++ ) {
      weno_weights_minimal_out[i].weights[weight_count] = MIN( 
                                                              weno_weights_current[i].weights[weight_count] * ( 1.0 - weno_weights_current[i].lower_order_fraction ), //rescale the weights by the weno5 fraction
                                                              weno_weights_minimal[i].weights[weight_count] );
    }
  }
}


// spread weight minimization
// CHANGINGMARK: Need to ensure that weight is computed outside where need final weight (i.e. need weight at ps-1 .. pe+1 at least assuming RECON only uses at ps..pe -- what about stencil reduction?)
// Much more robust!
// uses stencil_weights_array_other for temporary space since convolution requires temporary space
static void spread_weno5_weights_old( int preforder, weno_weights_t *weno_weights_in, weno_weights_t *weno_weights_out)
{
  int i, weight_count;
  int order = (preforder + 1) / 2;
  FTYPE left,right;
  int offleft,offright;
  // OPENMPMARK: Can't leave as static
  weno_weights_t a_stencil_weights_array_other[NBIGM];  //for storing the weights for a one-dimensional array of values; used by stencil reduction
  weno_weights_t *stencil_weights_array_other = &(a_stencil_weights_array_other[NBIGBND]); //GODMARK:  what should I use NBIGBND or sth else?


  // first copy over weights since in may equal out
  copy_weno5_weights(preforder, weno_weights_in, stencil_weights_array_other);


  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {

    offleft = (i>-NBIGBND) ? 1 : 0;
    offright = (i<NBIG+NBIGBND-1) ? 1 : 0;

    for( weight_count = 0; weight_count < order; weight_count++ ) {

      left = stencil_weights_array_other[i-offleft].weights[weight_count];
      right = stencil_weights_array_other[i+offright].weights[weight_count];

      weno_weights_out[i].weights[weight_count] = MIN(stencil_weights_array_other[i].weights[weight_count],left);
      weno_weights_out[i].weights[weight_count] = MIN(stencil_weights_array_other[i].weights[weight_count],right);

    }
  }
}




// spread weight minimization
// CHANGINGMARK: Need to ensure that weight is computed outside where need final weight (i.e. need weight at ps-1 .. pe+1 at least assuming RECON only uses at ps..pe -- what about stencil reduction?)
// uses stencil_weights_array_other for temporary space since convolution requires temporary space
// bit more robust but slightly less robust (more oscillatory) than above method -- although this method makes more sense spatially
static void spread_weno5_weights( int preforder, weno_weights_t *weno_weights_in, weno_weights_t *weno_weights_out)
{
  int i, weight_count;
  int order = (preforder + 1) / 2;
  FTYPE left,right,center,mintotal;
  int offleft,offright;
  // OPENMPMARK: Can't leave as static
  weno_weights_t a_stencil_weights_array_other[NBIGM];  //for storing the weights for a one-dimensional array of values; used by stencil reduction
  weno_weights_t *stencil_weights_array_other = &(a_stencil_weights_array_other[NBIGBND]); //GODMARK:  what should I use NBIGBND or sth else?

  // CHANGINGMARK: need 1 more point for a2c to use this (13 total ghost cells)
#if(MAXBND<13 && WENO_EXTRA_A2C_MINIMIZATION==1)
  // CHANGINGMARK: Should add check for using WENO5BNDPLUSMIN
  dualfprintf(fail_file,"spread_weno_weights() cannot be used with MAXBND<13\n");
  myexit(18672159);
#endif

#if(WENO_EXTRA_A2C_MINIMIZATION==0)
  return;
#endif


  // first copy over weights since in may equal out
  copy_weno5_weights(preforder, weno_weights_in, stencil_weights_array_other);


  if(order!=3){
    dualfprintf(fail_file,"spread_weno_weights() not setup for order!=3, order=%d\n",order);
    myexit(2762372);
  }

  for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ) {

    offleft = (i>-NBIGBND) ? 1 : 0;
    offright = (i<NBIG+NBIGBND-1) ? 1 : 0;

    // operating on each spatial positions
    // stencil accessed in space via 2,1,0 NOT 0,1,2
    // name refers to spatial location of centered part of stencil
    left = stencil_weights_array_other[i-offleft].weights[0];
    right = stencil_weights_array_other[i+offright].weights[2];
    center = stencil_weights_array_other[i].weights[1];

    mintotal=MIN(center,MIN(left,right));

    // now assign results
    weno_weights_out[i-offleft].weights[0] = mintotal;
    weno_weights_out[i+offright].weights[2] = mintotal;
    weno_weights_out[i].weights[1] = mintotal;

  }
}
