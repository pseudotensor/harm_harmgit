// NOTES/ASSUMPTIONS OF METHOD:
//
// 1) Assume smoothness indicator is square of sum of derivatives
//    Assume weights are 1 over square of smoothness indicators
//    As discussed in Jiang & Shu (1996), maybe this power p should be p=r=3 for WENO5
// 2) c2e and a2c are independent operations (but aren't really)
// 3) Epsilon matters always in determining reconstruction when close to very flat or switching between very flat and very sharp


///////
// 
//  Changes to the code
//
//  1. Modified the RT problem to evolve half a mushroom
//  2. Modified the makefile to use the -pc64 -unroll -mp for the icc to preserve perfect symmetry
//  3. Corrected format stings in pi2Uavg_specific()
//  4. Use REMOVERESTMASSFROMUU = 2 for test 151 as a try to see if it works better with a2c
//  5. Undo the last change because it does not make any difference
//  6. DOENODEBUG stuff implemented and corrected to give symmetric results. 
//  7. Correct the SMONO asymmetry issue: compute monoindicators over a larger range, [bs+order-1,be-order+1] instead of [ps,pe]
//  8. DOENODEBUG: correct all instances of enodebugarray...[dir]... to enodebugarray...[dir-1]... otherwise memory leaks
//  9. Fixed the number of boundary conditions issue to work correctly with dp/p reduction at the boundaries.  Before, the results were unpredictable
//     for the a2c reduction
// 10. Reduce in interpolation of conserved quantities in one direction if reduced in the other one
// 11. "Energy is the boss" type of reduction for the conserved quantities, "maximum coupling" type of reduction in the fluxes.
// 10. Fixed the analytic connection for the cylindrical geometry
// 11. Flip the sign of the conserved quantities & fluxes in ghost zones near the asymm boundary condition before it gets passed to the a2c/c2a
//     This should avoid a kink in conserved quantities/fluxes due to them being multiplied by detg = |R| that has a kink.  This effectively changes |R| to R.



/////////
//
//  TODO
//
//  1. do a2c at restart to compute the weights for energy is the boss and populate the primtive quantities
//  2. store the a2c weights for the energy and use that for primitives, conserved, and source terms
//  !!!!! IMMEDIATE !!!
//  3. For fluxes, it is the flux of momentum in the direction of interpolation that should control all weights
//  4. WENO_USE_PRIM_REDUCTION:  should distinguish between the averaging of fluxes and source terms


////////
//
//  READ
//
//  1.  Intermediate state due to the Roe solver (seems like they are able to 

/////////
//
//  List of ideas:
//
//  1. SMONO:  require to be monotonic not at [-1,1] but on one of [-1,2] OR [-2,1]. The latter way the requirement is more restrictive and may work better for non-smooth flows
//  2. JMONO vs. SMONO:  JMONO is much less restrictive and can get activated when SMONO is not even thinking about it.  Maybe it is a good idea to re-implement JMONO
//     within SMONO framework with the smooth switch.
//  3. Optimize source terms interpolation (only two of the source terms are non-zero in GR)
//  4. ---(done) pi2Uavg not call when avglim (introduce interporder apart from the number of points required)
//  5. SMONO not switch on when avglim = DONOR
//  6. create NUM_MIN as DBL_MIN or FLOAT_MIN in global.h (see NUM_EPSILON)
//  7. ???(done?) fix diag_fixup_U for FV

//  8. first compute the weights for all quantities, final weight will be some kind of combination of weights

//  9. Can results of test 151 for REMOVERESTMASSFROMUU == 2 be different from those with REMOVERESTMASSFROMUU == 1?

// 10. If a2c leads to a large difference between primitive quantities, maybe should take a linear combination in conserved space?

// Right now ideas
//  1. fractional difference should be modified so as not to reduce when the quantity goes through zero
//  2. +++ great +++ If fractional change is large, then reduce a2c/c2a and for fluxes as well
//  3. Do the limiting of a2c/c2a for primtive variables even if internal energy is negative

//  reduce easier for a2c (lower threshold for reduction)
//  limit the correction for the gamma-factor to 10% due to a2c?

//Reconstruction function updated to contain center to average conversions. Only change to .h file is to add eno_line_c2a() function declaration.
//List of updates
// Replace CVT_C2E -> cvt_type in choose_weno_order (for lower order)
// set w_ratio = BIG to avoid asymmetries and using undefined numbers


// roughness weights (inverse smoothness indicators \omega_j=1/\beta_j^2 normalized)
// optimal weights (d_j's)
// optimized weights (normalized \omega_j d_j)


#include "decs.h"


#include "reconstructeno_global.h"
#include "reconstructeno_defs.h"
#include "reconstructeno_staticfuncs_constvars.h"
#include "reconstructeno_global.funcdeclare.h"


//#if( TESTNUMBER == 153 && N2 == 1 && N3 == 1 )
// static FTYPE a_ucons_bondi[NBIGM][NPR];
// //shift the array index so that it starts at (-NBIGBND)
// static FTYPE (*ucons_bondi)[NPR] = (FTYPE (*)[NPR]) (&(a_ucons_bondi[NPR][NBIGBND]));
//#endif








#include "reconstructeno_set_arrays.c"








int eno_line_a2c( int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                  FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                  FTYPE *Pindicator, FTYPE *yin,  FTYPE *yout, struct of_trueijkp *trueijkp) 
//wrapper function that calls a more general conversion function, eno_line_a2c_c2a
{
  int i;
  int res;
  weno_weights_t one_sided_stencil_weight;
  int weight_no;
 

#if( DO_AVG2CEN == 0 )
  for( i = ps; i <= pf; i++ ) {
    if(monoindicator[CENTYOUT][i]==0) yout[i] = yin[i];
    //yout[i] = yin[i];
  }
  return( 0 );
#else

#if(COUNTCALLS)
  int sum,i;
  sum=0;
  for( i = ps; i <= pf; i++ ) {
    sum+=(monoindicator[CENTYOUT][i]==0);
  }
  dualfprintf(fail_file,"t=%21.15g nstep=%ld pl=%d sum_a2c=%d\n",t,nstep,pl,sum);
#endif

  //#if( TESTNUMBER == 153 && N2 == 1 && N3 == 1 )
  // //for Bondi problem, modify the average values of conserved quantities for X1UP boundary zones to be equal to those at t = 0; this assumes 1d case
  // if( mycpupos[1] == ncpux1 - 1 ) {
  //  for( i = N1; i <= pf; i++ ) {
  //   yin[i] = ucons_bondi[pl][i];
  //  }
  // }
  //#endif

  res = eno_line_reconstruct( whichquantity, do_weight_or_recon, stencil_weights_array, CVT_A2C, whichreduce, 
                              preforder, pl, bs, ps, pf, bf, 
                              minorderit, 
                              maxorderit, 
                              shiftit, 
                              shockindicator,
                              df, dP, monoindicator, Pindicator, yin, yout, NULL, NULL,trueijkp ); 


#if(0) // deprecated with no iterglobal anymore
  /////////////////
  //
  // perform special a2c reconstruction near the real boundaries with special outflow boundary conditions: 
  // avoid using ghost active cells (with artificially evolved conserved quantities) in the a2c reconstruction
  // for a2c reconstruction in active cells for which information in ghostactive cells is required.
  // Instead, use the 5-point stencil that touches the boundary of active cells zone but does not stick out of it.
  // For WENO-5 need to do this only for two points that are closest to the boundary of active cells region.
  //
  // low-r boundary
  if( 0 && iterglobal == 1 && mycpupos[1] == 0 && (BCtype[X1DN] == BONDIMDOTOUTFLOW || BCtype[X1DN] == BONDIINTOUTFLOW) ) {
    one_sided_stencil_weight.order = 5;
    one_sided_stencil_weight.len = 10;

    i = 0;
    //construct the weights that would use the rightmost (0th from right) stencil of five points (order = 5) so that it does not touch ghost cells
    for( weight_no = 0; weight_no < 10; weight_no++ ) {
      one_sided_stencil_weight.weights[weight_no] = (weight_no % 5 == 0);
    }
    eno_cvt( CVT_A2C, bs - i, bf - i, &one_sided_stencil_weight, &yin[i], &yout[i] );

    i = 1;
    //construct the weights that would use the 1st from the right stencil of five points (order = 5)
    for( weight_no = 0; weight_no < 10; weight_no++ ) {
      one_sided_stencil_weight.weights[weight_no] = (weight_no % 5 == 1);
    }
    eno_cvt( CVT_A2C, bs - i, bf - i, &one_sided_stencil_weight, &yin[i], &yout[i] );
  }
  // high-r boundary
  if( 0 && iterglobal == 1 && mycpupos[1] == ncpux1 - 1 && (BCtype[X1UP] == BONDIMDOTOUTFLOW || BCtype[X1UP] == BONDIINTOUTFLOW) ) {
    one_sided_stencil_weight.order = 5;
    one_sided_stencil_weight.len = 10;

    i = N1 - 1;
    //construct the weights that would use the 0th from the left (4th from the right) stencil of five points (order = 5)
    for( weight_no = 0; weight_no < 10; weight_no++ ) {
      one_sided_stencil_weight.weights[weight_no] = (weight_no % 5 == 4);
    }
    eno_cvt( CVT_A2C, bs - i, bf - i, &one_sided_stencil_weight, &yin[i], &yout[i] );

    i = N1 - 2;
    //construct the weights that would use the 1st from the left (3rd from the right) stencil of five points (order = 5)
    for( weight_no = 0; weight_no < 10; weight_no++ ) {
      one_sided_stencil_weight.weights[weight_no] = (weight_no % 5 == 3);
    }
    eno_cvt( CVT_A2C, bs - i, bf - i, &one_sided_stencil_weight, &yin[i], &yout[i] );
  }
  //
  ////////////////////
#endif
  
  return res;

#endif
}









int eno_line_c2a( int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                  FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                  FTYPE *Pindicator, FTYPE *yin,  FTYPE *yout, struct of_trueijkp *trueijkp) 
//wrapper function that calls a more general conversion function, eno_line_reconstruct
{
  int i;
  int res;

#if( DO_CEN2AVG == 0 )
  for( i = ps; i <= pf; i++ ) {
    if(monoindicator[CENTYOUT][i]==0) yout[i] = yin[i];
  }
  return( 0 );
#else

  res = eno_line_reconstruct( whichquantity, do_weight_or_recon, stencil_weights_array, CVT_C2A, whichreduce, preforder, pl, bs, ps, pf, bf, minorderit, maxorderit, shiftit, shockindicator,  df, dP, monoindicator, Pindicator, yin, yout, NULL, NULL,trueijkp ); 

  //#if( TESTNUMBER == 153 && N2 == 1 && N3 == 1 )
  // //for Bondi problem, save the average values of conserved quantities at t = 0; this assumes 1d case
  // if( is_pi2Uavg ) {
  //  for( i = ps; i <= pf; i++ ) {
  //   ucons_bondi[pl][i] = yout[i];
  //  }
  // }
  //#endif

  return res;

#endif
}






int eno_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                  FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM],
                  FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                  FTYPE *Pindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp ) 
{

#if( DO_CENT2EDGE == 0 )
  for( i = ps; i <= pf; i++ ) {
    if(monoindicator[LEFTYOUT][i]==0) yout_left[i] = yin[i];
    if(monoindicator[RIGHTYOUT][i]==0) yout_right[i] = yin[i];
    //yout_left[i] = yin[i];
    //    yout_right[i] = yin[i];
  }
  return( 0 );
#else

#if(COUNTCALLS)
  int sum;
  sum=0;
  for( i = ps; i <= pf; i++ ) {
    sum+=(monoindicator[LEFTYOUT][i]==0)+(monoindicator[RIGHTYOUT][i]==0);
  }
  dualfprintf(fail_file,"t=%21.15g nstep=%ld pl=%d sum_c2e=%d\n",t,nstep,pl,sum/2);
#endif



  return(eno_line_reconstruct( whichquantity, do_weight_or_recon, stencil_weights_array, CVT_C2E, whichreduce, preforder, pl, bs, ps, pf, bf, minorderit, maxorderit, shiftit, shockindicator,  df, dP, monoindicator, Pindicator, yin, yout_left, yout_right, youtpolycoef,trueijkp ));


#endif
}









// define range over which to check VSQ
#if(0)
#define GAMCHECK(i) ( (P[i] > 10 || P[i-1] > 10 || P[i+1] > 10 || P[i+2] > 10 || P[i-2] > 10 || P[i-3] > 10 || P[i+3]>10 || P[i-4]>10 || P[i+4] >10  || P[i-5]>10 || P[i+5]>10  ))
#elif(VSQ!=-100)
#define GAMCHECK(i) ( P[i] > 10 )
#else
#define GAMCHECK(i) (0)
#endif



///Interpolates a set of numbers provided in #yin into one or two sets of numbers #yout_left and, depending on the setting of cvt_type, to #yout_right.
/// \return a set of interpolated values using the interpolation of type #cvt_type
/// \param[in] cvt_type  can be one of  
///       - #CVT_C2E -- reconstruction from centers to edges (yin[] is the values at the centers and yout_left[] and yout_right[] are the values at the edges)
///       - #CVT_A2C -- reconstruction from centered average to centered point quantities.  yin[] and yout_left[] are input and output lines respectively (yout_right[] not used).
///       - #CVT_C2A -- reconstruction from centered point to centered average quantities.  yin[] and yout_left[] are input and output lines respectively (yout_right[] not used).
/// \param[in] whichreduce  can be one of #WENO_REDUCE_TYPE_PPM or #WENO_REDUCE_TYPE_DEFAULT.  Sets the type of stencil reduction to be used.
/// \param[in] bs  first element index where you can read from
/// \param[in] ps  first element index for which a reconstructed value is desired
/// \param[in] pf  last element index for which a reconstructed value is desired
/// \param[in] bf  last element index where you can read from
/// \param[in] minorderit
/// \param[in] maxorderit
/// \param[in] shiftit
/// \param[in] yin  array of quantities to be interpolated
/// \param[out] yout_left  array of interpolated quantities (used for all types of reconstructions)
/// \param[out] yout_right  array of interpolated quantities (used only for those reconstructions that have two outputs, e.g. center to edge reconstruction)
int eno_line_reconstruct(int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int cvt_type, int whichreduce, int preforder, int pl, 
                         int bs, int ps, int pf, int bf, 
                         int *minorderit, int *maxorderit, int *shiftit, 
                         FTYPE (*shockindicator)[NBIGM], 
                         FTYPE (*df)[NBIGM],
                         FTYPE (*dP)[NBIGM],
                         FTYPE (*monoindicator)[NBIGM],
                         FTYPE *Pindicator, 
                         FTYPE *yin,  FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp)
{
  

  int i;
  FTYPE simpleweight;
  int counter;
  int weights_start; //the first point where the weights can be computed
  int weights_finish; //the last point where the weights can be computed
  int maxorder;
  int minorder;
  int order_to_be_used;
  FTYPE yout_delta;
  FTYPE yout_left_temp, yout_right_temp; 
  weno_weights_t *p_stencil_weights_array_to_be_used;
  int max_full_stencil_size;
  int min_full_stencil_size;
  int shift;
  int num_orders;   //the number of reconstructions of different orders whose linear combination is to be taken as the final reconstruction
  int weno_output_count;
  FTYPE *yout[MAX_NO_OF_WENO_OUTPUTS];
  FTYPE output_val;
  FTYPE slope;
  FTYPE interface_value;

  FTYPE fraction_point_value;

  FTYPE monoind_val;

  int derorder;
  int dertype;
  weno_weights_t a_stencil_weights_array_static[NBIGM]; // OPENMPMARK: Can't leave as static
  weno_weights_t *stencil_weights_array_static = &a_stencil_weights_array_static[NBIGBND];

#if(DOENODEBUG || STORE_GAMMA_PRIM_REDUCTION_FRACTION)  //atch enodebug
  int di, dj, dk;
  //  ITERGLOBALDEF  //define di, dj, dk, etc.

  di = (trueijkp->iter==1); dj = (trueijkp->iter==2); dk = (trueijkp->iter==3);
#endif
  
  //////////////////////////////
  //
  // Convert input to internal bounds format and check consistency of boundaries and check symmetry of intput data
  //
  //////////////////////////////

  max_full_stencil_size=preforder;
  min_full_stencil_size=minorderit[ps];
  shift=shiftit[ps];

  assert( bs < - MAXBND, "eno_line_c2f: bs out of bounds" );
  assert( bs > ps || bf < pf, "eno_line_c2f: limits out of order. Should be: bs <= ps <= bf <= pf" );
  assert( max_full_stencil_size < min_full_stencil_size || min_full_stencil_size < 1, "min_full_stencil_size cannot be smaller than min_full_stencil_size or < 1" );
  assert( stencil_weights_array == NULL && do_weight_or_recon == RECON_CALC, "Should provide stencil weights if only doing reconstruction" );

  if( NULL == stencil_weights_array ) {
    //use temporary static space for storing weights
    stencil_weights_array = stencil_weights_array_static;
  }

  //convert the full stencil size to the substencil size
  minorder = ( min_full_stencil_size + 1 ) / 2;
  maxorder = ( max_full_stencil_size + 1 ) / 2;

  assert( maxorder < MIN_CVT_ORDER, "requested order is too small" );

  ////////////////////////
  // Check symmetry of input

#if(SSYMCHECKS)
  //check symmetry of provided solution; if asymmetric, returns 1 and prints out a message
  check_symmetry_of_weno_input( cvt_type, pl, bs, bf, yin );
#endif

  //////////////////////////////
  //
  // Compute range over which weights are computed so all calculations are within boundaries of data
  //
  //////////////////////////////


  weights_start = bs + maxorder - 1;     //the first point where all weno weights can be computed (i.e. default weight calc. will not require the values outside the boundaries)
  weights_finish = bf - (maxorder - 1);  //the last point where all weno weights can be computed

  //adjust the above values if the user requires weights closer to the boundary -- then the weights of those stencils that stick out of the boundary will be zero
  if( weights_start > ps ) { 
    weights_start = ps;
  }

  if( weights_finish < pf ) {
    weights_finish = pf;
  }

  if( do_weight_or_recon == WEIGHT_CALC || do_weight_or_recon == ALL_CALC) { //peform the unoptimized weight calculation
    //////////////////////////////
    //
    // Compute unoptimized weights
    //
    //////////////////////////////

    for( i = weights_start; i <= weights_finish; i++ ) {
      //compute weights for every point and store them in stencil_weights_array[] for future use
      compute_stencil_weights( cvt_type, pl, maxorder, bs - i, bf - i, &monoindicator[MONOYIN][i], &monoindicator[LEFTYOUT][i], &yin[i], &stencil_weights_array[i] );
      
      //compute lower-order (currently, WENO-3 supported) weights and place them into the high-order weights structure (currently, WENO-5 weights)
      compute_lower_order_stencil_weights( cvt_type, pl, maxorder, bs - i, bf - i, &monoindicator[MONOYIN][i], &monoindicator[LEFTYOUT][i], &yin[i], &stencil_weights_array[i] );
      
#if( DO_RESCALE_WEIGHTS )
      rescale_unoptimized_weights( &stencil_weights_array[i] );
#endif
#if( DO_SMOOTH_ADJUSTMENT_UNOPTIMIZED_WEIGHTS )
      mono_adjust_unoptimized_weights( &monoindicator[MONOYIN][i], &stencil_weights_array[i] );
#endif 
      //if( monoindicator[MONOYIN][i] != monoindicator[MONOYIN][N1 - 1 - i] ) {
      // dualfprintf( fail_file, "Asymmetry in monoindicator: monoindicator[MONOYIN][%d] = %21.15g, monoindicator[MONOYIN][%d] = %21.15g\n",
      //   i, monoindicator[MONOYIN][i], N1 - 1 - i, monoindicator[MONOYIN][N1 - 1 - i]);
      //}
    } 

#if(0)
    for( i =0; i < N1; i++ ) {
      if((pl==RHO||pl==U1)&&(steppart==0))
        dualfprintf(fail_file,"%d %d %ld %d %21.15g %21.15g %21.15g %21.15g\n",pl,i,nstep,steppart,stencil_weights_array[i].weights[2],stencil_weights_array[i].weights[1],stencil_weights_array[i].weights[0],yin[i]);
    }
#endif


    for( i = ps-1; i <= pf+1; i++ ) {  //extended to go from ps - 1; pf+1 because additional zones are required by dp/p reduction
      //comutes the weights ratios as well as the lower order fractions for every point based on these ratios
      //note that P[i] contains the value of the lorentz factor (termporary use this array for storing it since no need for pressure)
      if(GAMCHECK(i) &&  cvt_type == CVT_C2E) {
        //use CVT_A2C so that reduction proceeds as the OR reduction
        compute_weights_ratios( CVT_A2C, bs - i, bf - i, &yin[i], &stencil_weights_array[i] ); 
      }
      else {
        compute_weights_ratios( cvt_type, bs - i, bf - i, &yin[i], &stencil_weights_array[i] ); 
      }
    }
    
    for( i = ps; i <= pf; i++ ) {  
      //modify the lower order fraction by applying extra reduction:
      //dp/p reduction, etc. to the weno5 weights and update the lower order values
      apply_additional_reduction_to_weights( cvt_type, whichreduce, maxorder, minorder, i, pl, bs, bf, shockindicator, 
                                             dP, monoindicator[MONOYIN], Pindicator, yin, stencil_weights_array,trueijkp );

      //JCM
      // stencil_weights_array[i].lower_order_fraction;
 
#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )  //SUPERSASMARK
      //reduce the order of reconstruction for interpolation of conserved quantities, source terms and fluxes 
      //if the order of interpolation of gamma has been reduced:
      //SASMARKx: bad idea to use for the interpolation of conserved quantities because gamma is at a previous time step
      if( pl == VSQ && cvt_type == CVT_C2E ) {  

        //if( iterglobal == 2 && stencil_weights_array[i].lower_order_fraction != 0 ) {
        // trifprintf( "gamma reduced, dir = 2, lo_ord_fr = %lg, nstep = %d, steppart = %d, i = %d, j = %d (symj = %d)\n", 
        //  stencil_weights_array[i].lower_order_fraction, nstep, steppart, trueijkp->i, i, N2 - 1 - i );
        //}

        //do this only for gamma, which is the last, NPR (nprend)-th, primitive; STORE_GAMMA_PRIM_REDUCTION_FRACTION == 1 means that gamma actually exists and is being interpolated
        //store the largest value of lower_order_fraction for the current point globally
        //MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) = 
        //   MAX( stencil_weights_array[i].lower_order_fraction,
        //      MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) );  

        //SUPERSASMARK: put in the check for cusp into the lower_order_fraction
        if( monoindicator[MONOYIN][i] == 0 && (monoindicator[LEFTYOUT][i] == 1 || monoindicator[RIGHTYOUT][i] == 1) ) {
          if( fabs(yin[i]-yin[i-1]) + fabs(yin[i+1]-yin[i]) > 0.01 * yin[i] ) {
            MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) = 1.0;
          }
        }
      }
      else if( cvt_type == CVT_A2C ) {  //deaveraging of conserved quantities; use the lower order fractions
        //reduce for a2c if reduced for the gamma
        //stencil_weights_array[i].lower_order_fraction = 
        // MAX( stencil_weights_array[i].lower_order_fraction,
        //     MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) );
      }
      else if( cvt_type == CVT_C2A ) { 
        if( whichquantity == ENOFLUX ) {
          //averaging of fluxes
      
          //reduce for c2a if reduced for interpolation of gamma ...
          //... for the cell to the right of the interface reduced  (no need to convert the indices)
          stencil_weights_array[i].lower_order_fraction = 
            MAX( stencil_weights_array[i].lower_order_fraction,
                 MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) );
          //... for the cell to the left of the interface reduced  (go left in the direction of dirglobal)
          stencil_weights_array[i].lower_order_fraction = 
            MAX( stencil_weights_array[i].lower_order_fraction,
                 MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di - (dirglobal==1),trueijkp->j + i*dj - (dirglobal==2),trueijkp->k + i*dk - (dirglobal==3)) );
        }
        //else { //averaging or primitive quantities or source terms
        // //reduce for c2a if reduced for the gamma
        // stencil_weights_array[i].lower_order_fraction = 
        //  MAX( stencil_weights_array[i].lower_order_fraction,
        //      MACP1A0(weno_prim_lower_order_fraction,dirglobal,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk) );
        //}
      }
      else if( cvt_type != CVT_C2E ) {
        dualfprintf( fail_file, "Unknown reconstruction type\n" );
        myexit(1);
      }
#endif
    }
  } //end if( do_weight_or_recon )

  if( do_weight_or_recon == RECON_CALC || do_weight_or_recon == ALL_CALC ) {
    //////////////////////////////
    //
    // Compute weights
    //
    //////////////////////////////

    for( i = weights_start; i <= weights_finish; i++ ) {
      //compute optimized weights for every point and store them in stencil_weights_array[] for future use
      compute_optimized_stencil_weights( cvt_type, pl, maxorder, bs - i, bf - i, &monoindicator[MONOYIN][i], &yin[i], &stencil_weights_array[i] );
    }

    //////////////////////////////
    //
    // Choose order using weights
    //
    // And then for each point compute the interpolated values
    //
    //////////////////////////////


    //initialize the temporary variables; will be adding pieces from reconstructions of different orders to them
    yout[0] = yout_left;
    yout[1] = yout_right;  //this can be NULL for a2c/c2a recs, so make sure don't write to it unless the reconstruction type requires us to do so

    for( i = ps; i <= pf; i++ ) {
#if( DO_ENO_STENCIL_REDUCTION )
      //for each point: choose the order of interpolation; there can be a linear combination of several orders; the number of orders is returned by the below call
      num_orders = choose_weno_order( cvt_type, whichreduce, maxorder, minorder, i, pl, bs, bf, shockindicator, dP, monoindicator[MONOYIN], monoindicator[LEFTYOUT], Pindicator, yin, stencil_weights_array, &p_stencil_weights_array_to_be_used,trueijkp );
#else
      //stencil reduction is turned off, so only use weights of the default order; they have been computed above
      num_orders = 1;
      p_stencil_weights_array_to_be_used = &stencil_weights_array[i];
#endif
   
   
      for( weno_output_count = 0; weno_output_count < weno_outputs[cvt_type].len; weno_output_count++ ) {
        //init output val with zero if SMONO is not used <-- SASMARK don't need this if yout is well-defined 
        //(not nan's, infitinities, etc. so that 0 * yout[c][i] = 0)
        assert( monoindicator[LEFTYOUT+weno_output_count][i] < 0, "monoindicator[%d][%d] = %lg\n", LEFTYOUT+weno_output_count, i, monoindicator[LEFTYOUT+weno_output_count][i] );
        /////
        //Original version:
        //if( monoindicator[LEFTYOUT+weno_output_count][i] == 0.0 ) {
        /////
        /////
        
        monoind_val = monoindicator[MONOYIN][i];

        //Initialize output values if they have not been assigned according to 5th point stencil in (S/J)SMONO
        if( monoind_val == 0.0 ) {
          yout[weno_output_count][i] = 0.0;
        }
#if(DOENODEBUG)
        else if( dirglobal<=2 && pl <=U2 && monoindicator[LEFTYOUT+weno_output_count][i] > 0.01 ){ //monoindicator is used
          MACP0A4(enodebugarray,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk,dirglobal-1,interporfluxglobal,pl,ENODEBUGPARAM_SMONO)++;
        }
#endif
        ////////
        //The below is the original version that does not assign values in cusps:
        // don't assign value if already assigned before entering here (works for dual or single output)
        //if( monoindicator[LEFTYOUT+weno_output_count][i] >= 0.0 && monoindicator[LEFTYOUT+weno_output_count][i] < 1.0 ){
        ////////

        ////////
        //This is the new version that assigns the values in cusps (uses WENO3 as is set at the top of compute_..._weights())
        ////////
        if( monoind_val >= 0.0 && monoind_val < 1.0 ){
          //the value has not been assigned outside of this function, so assign it
          output_val = 0.0;  //initialize with zero
          //loop over the reconstructions whose linear combination we will be given as the answer
          for( counter = 0; counter < num_orders; counter++ ) {  
            //if( pl < U2 ) dualfprintf( fail_file, "pl = %d, nstep = %ld, steppart = %d, i = %d, order = %d, w0 = %g, w1 = %g, w2 = %g\n",
            // pl, nstep, steppart, i,
            // p_stencil_weights_array_to_be_used[counter].order, 
            // p_stencil_weights_array_to_be_used[counter].weights[0+p_stencil_weights_array_to_be_used[counter].order], 
            // p_stencil_weights_array_to_be_used[counter].weights[1+p_stencil_weights_array_to_be_used[counter].order], 
            // p_stencil_weights_array_to_be_used[counter].weights[2+p_stencil_weights_array_to_be_used[counter].order] );
            eno_cvt( weno_outputs[cvt_type].type[weno_output_count], bs - i, bf - i, &p_stencil_weights_array_to_be_used[counter], &yin[i], &yout_delta );
            output_val += yout_delta;
#if( DOENODEBUG )
            //count the use of WENO3/WENO5 in the SWENO if the fraction of that reconstruction is larger than 1%:
            //
            // WENO3:
            if(  dirglobal<=2 && pl <=U2 && 
                 p_stencil_weights_array_to_be_used[counter].order == 2 && get_sum_of_elements(2, p_stencil_weights_array_to_be_used[counter].weights) > 0.01 ) {
              MACP0A4(enodebugarray,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk,dirglobal-1,interporfluxglobal,pl,ENODEBUGPARAM_WENO3)++;
            }
            // WENO5:
            else if( dirglobal<=2 && pl <=U2 && 
                     p_stencil_weights_array_to_be_used[counter].order == 3 && get_sum_of_elements(3, p_stencil_weights_array_to_be_used[counter].weights) > 0.01 ) {
              MACP0A4(enodebugarray,trueijkp->i + i*di,trueijkp->j + i*dj,trueijkp->k + i*dk,dirglobal-1,interporfluxglobal,pl,ENODEBUGPARAM_WENO5)++;
            }
#endif
          }
          ////////
          //The below is the original version that does not assign values in cusps:
          // don't assign value if already assigned before entering here (works for dual or single output)
          //
          //yout[weno_output_count][i] = ( 1.0 - monoindicator[LEFTYOUT+weno_output_count][i] ) * output_val 
          //                  + monoindicator[LEFTYOUT+weno_output_count][i] * yout[weno_output_count][i];
          ////////

          yout[weno_output_count][i] = ( 1.0 - monoindicator[MONOYIN][i] ) * output_val 
            + monoindicator[MONOYIN][i] * yout[weno_output_count][i];

        }
      } 

    }  //end for

#if( MERGEDC2EA2CMETHOD )
    //Find the values of derivatives due to polynomial reconstruction
    //Do this only for c2e reconstruction and only if youtpolycoef array is non-NULL
    if( cvt_type == CVT_C2E && NULL != youtpolycoef ) {
      //zeroth derivative = value itself
      for( i = ps; i <= pf; i++ ) {
        youtpolycoef[0][i] = yin[i];
      }
      //Fill in 1st - 4th derivative values by calling itself recursively
      for( derorder = 1; derorder <= 4; derorder++ ) {
        dertype = CVT_C2DER1 + derorder - 1;
        eno_line_reconstruct( whichquantity, RECON_CALC, stencil_weights_array, dertype, whichreduce, 
                              preforder, pl, bs, ps, pf, bf, 
                              minorderit, 
                              maxorderit, 
                              shiftit, 
                              shockindicator,
                              df, dP, monoindicator, Pindicator, yin, youtpolycoef[derorder], NULL, NULL,trueijkp ); 
      }
    }
#endif

    ////////////////////////
    //
    // Check symmetry of output

#if(SSYMCHECKS)
    check_symmetry_of_weno_output( cvt_type, monoindicator, pl, ps, pf, maxorder, yout, stencil_weights_array );
#endif

  }  //end doing the reconstruction


  return( 0 );
}








//forces unoptimized weights to be equal if they are different from 
//resets the unoptimized weights to 1/order if the weight is between 1/order - epsilon and 1/order+order*epsilon, lintearly interpolates outside this interval so that 0 -> 0 and 1 -> 1
void rescale_unoptimized_weights( weno_weights_t *weights ) 
{
  int weight_no;
  FTYPE min_weight = 1.;
  FTYPE max_weight = 0.;
  FTYPE frac_equal_weight = 1.;
  FTYPE equal_weight = 1. / weights->order;

  for( weight_no = 0; weight_no < weights->order; weight_no++ ) { 
    //weights->weights[weight_no] = rescale_weight( weights->order, weights->weights[weight_no] );
    min_weight = MIN( weights->weights[weight_no], min_weight );
    max_weight = MAX( weights->weights[weight_no], max_weight );
  }

  frac_equal_weight = MIN( compute_frac_equal_weight( weights->order, min_weight ), frac_equal_weight );
  frac_equal_weight = MIN( compute_frac_equal_weight( weights->order, max_weight ), frac_equal_weight );

  for( weight_no = 0; weight_no < weights->order; weight_no++ ) { 
    weights->weights[weight_no] = equal_weight * frac_equal_weight + (1.0 - frac_equal_weight) * weights->weights[weight_no];
  }
}




FTYPE compute_frac_equal_weight( FTYPE order, FTYPE weight ) 
{
  FTYPE equal_weight = 1./order;
  FTYPE epsilon = 0.5 * equal_weight;
  FTYPE weight_min = equal_weight - epsilon;
  FTYPE weight_start = equal_weight - 0.5 * epsilon;
  FTYPE weight_end = equal_weight + 0.5 * (order - 1) * epsilon;  //so that if all weights except one are equal to weight_start, the left over one will be equal to weight_end
  FTYPE weight_max = equal_weight + (order - 1) * epsilon;
  FTYPE ans;

  if( weight < weight_min ) {
    ans = 0.;
  }
  else if( weight < weight_start ) {
    ans =  (weight - weight_min) / (weight_start - weight_min);
  }
  else if( weight < weight_end ) {
    ans = 1.;
  }
  else if( weight < weight_max ) {
    ans = (weight_max-weight)/(weight_max-weight_end);
  }
  else {
    ans = 0.0;
  }
  
  return( ans );
}





//returns 1/order if the number is between 1/order - epsilon and 1/order+order*epsilon, lintearly interpolates outside this interval so that 0 -> 0 and 1 -> 1
FTYPE rescale_weight( FTYPE order, FTYPE weight ) 
{
  FTYPE epsilon = 0.03;
  FTYPE weight_min = 1.e-3;
  FTYPE weight_start = 1./order - epsilon;
  FTYPE weight_end = 1./order + (order - 1) * epsilon;  //so that if all weights except one are equal to weight_start, the left over one will be equal to weight_end
  FTYPE weight_max = 1.;
  FTYPE ans;

  if( weight < weight_min ) {
    ans = 0.;
  }
  else if( weight < weight_start ) {
    ans =  (weight - weight_min) / ( order*(weight_start - weight_min) );
  }
  else if( weight < weight_end ) {
    ans = 1./order;
  }
  else if( weight < weight_max ) {
    ans = 1./order*( (order-1)*(weight-weight_end)/(weight_max-weight_end)+1.0 );
  }
  else {
    ans = 1.0;
  }
  
  return( ans );
}





void apply_additional_reduction_to_weights( int cvt_type, int whichreduce, int max_order, int min_order, int i0, int pl, int bs, int bf, FTYPE (*shockindicator)[NBIGM], FTYPE (*dP)[NBIGM],
                                            FTYPE *monoindicator, FTYPE *Pindicator, FTYPE *uin, 
                                            weno_weights_t *stencil_weights_array, struct of_trueijkp *trueijkp ) 
{
  FTYPE dPi;
  int counter;
  FTYPE lower_order_fraction;
  FTYPE shock_strength;
  int i;

#if( DOENODEBUG )
  FTYPE lower_order_fraction_old;  //atch debug
#endif 
  int di, dj, dk;
  //  ITERGLOBALDEF  //define values of di, dj, dk, etc.
  di = (trueijkp->iter==1); dj = (trueijkp->iter==2); dk = (trueijkp->iter==3);

  lower_order_fraction = stencil_weights_array[i0].lower_order_fraction;

#if( DO_REDUCE_POST_SHOCK || DO_REDUCE_PRE_SHOCK )
  //Reduce in/against the direction of the pressure growth
  //dPi = dP[0][i0+1] + dP[0][i0];
  //di = ( 2 * (dPi < 0.0) - 1 );
  dPi = dP[0][i0];

#if( DOENODEBUG )
  lower_order_fraction_old = lower_order_fraction;  //save to see if it got changed (for debug info)  //atch debug
#endif 

  shock_strength = MIN( (fabs(dPi) - 0.5 * REDUCE_DP_FRACTION) / (0.5 * REDUCE_DP_FRACTION), 1.);
  lower_order_fraction = MAX( lower_order_fraction, shock_strength * stencil_weights_array[i0+1].lower_order_fraction_fordpp );
  lower_order_fraction = MAX( lower_order_fraction, shock_strength * stencil_weights_array[i0-1].lower_order_fraction_fordpp );

#if( DOENODEBUG )
  //if lower order fraction changed by more than 50% due to dP/P stuff and > 0.01, then count it //atch debug
  if(  dirglobal<=2 && pl <=U2 && lower_order_fraction > MAX(1.5 * lower_order_fraction_old, 0.01)  ) {
    MACP0A4(enodebugarray,trueijkp->i + i0*di,trueijkp->j + i0*dj,trueijkp->k + i0*dk,dirglobal-1,interporfluxglobal,pl,ENODEBUGPARAM_dPP)++;
  }
#endif 

#endif

#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  if( do_weno_lower_order_fraction ) {
    lower_order_fraction = MAX( lower_order_fraction, GLOBALMACP0A1(weno_lower_order_fraction,trueijkp->i + i0*di,trueijkp->j + i0*dj,trueijkp->k + i0*dk,pl) );  //atch
    GLOBALMACP0A1(weno_lower_order_fraction,trueijkp->i + i0*di,trueijkp->j + i0*dj,trueijkp->k + i0*dk,pl) = lower_order_fraction; //store the largest value of lower_order_fraction globally //atch
  }
#endif 

  stencil_weights_array[i0].lower_order_fraction = lower_order_fraction; //update the lower_order_fraction in the weights array so that it can be reused later

}





void mono_adjust_unoptimized_weights( FTYPE *monoindicator, weno_weights_t *stencil_weights_array )
{
  FTYPE simpleweight;
  FTYPE order;
  int i;

  order = stencil_weights_array->order;

  //cases of monindicator == 0 and == 1 are already taken care of inside the function call above
  //SASMARK use of monoindicator[0]
  assert( monoindicator[0] < 0, "monoindicator[0] = %lg\n", monoindicator[0] );
  if( order == 3 && monoindicator[0] > 0.0 && monoindicator[0] < 1.0 ){ 
    //only do this for SWENO-5 and monoindicator strictly between 0 & 1; 
    //if equal to 0, SWENO reconstruction is fully used, so no need to adjust the weights;
    //if equal to 1, SMONO reconstruction is fully used, which means that the unoptimized weights have already been set before (= simpleweight) in compute_stencil_weights
    simpleweight = 1.0 / order;
    //adjust the *only* the unoptimized weights, which determine *only* stencil reduction (there are [order] of them at the beginning):
    //this smooths the jumps in the weights distribution used by stencil reduction when SMONO switches on and off; the reconstruction
    //procedure is so far unaffected because optimazed weights are used there
    for( i = 0; i < order; i++ ) {  
      //in the limit of monoindicator == 1.0 (SMONO fully on), this gives the equal weights (the weights are consistent with 
      //it being on and using the full 5-point polynomial as they should), 
      //in the limit of monoindicator == 0.0 (SMONO off), this leaves the SWENO weights unmodified as it should
      stencil_weights_array->weights[i] = ( 1.0 - monoindicator[0] ) * stencil_weights_array->weights[i] + monoindicator[0] * simpleweight;
    }
  } 
}





int eno_cvt( int cvt_type, int min_index, int max_index, weno_weights_t *stencil_weights_struct, FTYPE *uin, FTYPE *uout )  
//assumes that stencil_weights[] contains num_weights sets of weights, each of length (order_to_be_interpolated)
//the outputs for each of the weights are written to the uout[] array one after another, so uout[] should be large enough to hold num_weights elements
{
  FTYPE *cvtmatrix;
  int shift, symmetric_shift;
  int j;
  int order;
  FTYPE interpolated_value = 0.0;
  FTYPE *stencil_weights;
 
  FTYPE weno_stencil_recs[MAX_CVT_ORDER]; //SASMARK for storing updates due to each stencil
  FTYPE a2cupdate, this_order_weight;
  FTYPE weno_point_updates[MAX_CVT_ORDER];    //atch symmetrize
  FTYPE weno_stencil_updates[MAX_CVT_ORDER];  //atch symmetrize
 


  //////////////////////////////
  //
  // Consistency checks on input values and boundary limits on data
  //
  //////////////////////////////

  assert( NULL == stencil_weights_struct, "eno_cvt: stencil_weights_struct == NULL.\n" );
  assert( NULL == uout, "eno_cvt: FTYPE *uout cannot be a NULL pointer" );
  assert( NULL == uin, "eno_cvt: FTYPE *uin cannot be a NULL pointer" );

  //////////////////////////////
  //
  // Set pointer of weight structure to optimized weights for "order" and cvt_type of reconstruction
  //
  // Get conversion matrix (adding up points within substencil)
  //
  //////////////////////////////


  order = stencil_weights_struct->order;  //assign the order from the weights structure
  stencil_weights = stencil_weights_struct->weights + order * weno_weights_shifts_array[cvt_type];  //points to the start of weight sequence

  cvtmatrix = cvt_matrices[cvt_type][order]; //conversion matrix for that order
  assert( NULL == cvtmatrix, "eno_cvt: conversion matrix not implemented for the order specified\n" );

  //////////////////////////////
  //
  // Add up points within substencil
  //
  // Add up substencils
  //
  //////////////////////////////


  for( shift = 0; shift < order; shift++ ) { 
    //loop through all possible stencils
    //compute the updates due to each point in the stencil
    for( j = 0; j < order; j++ ) {
      //cvt_matrix[shift * order + j] effectively gives cvt_matrix[shift][j]
      weno_point_updates[j] = MATRIX_ROW_COL(cvtmatrix, shift, j, order) * uin[j - shift]; //see (2.11) from Shu Report (1997)
    }

    //initialize with zero to avoid an if-statement later
    weno_stencil_updates[shift] = 0.0;  

    //To get the update due to the whole stencil, need to sum up the above updates; however, this should be performed in a symmetric fashion
    if( stencil_weights[shift] != 0.0 ) {  //SASMARK: if values just outside the boundary are NAN's, then zero time NAN gives a NAN, and the code fails, so use isfinite() or != 0.0 to check for that.
      //reconstructed value due to one stencil
      weno_stencil_recs[shift] = get_sum_of_elements( order, weno_point_updates ); 
      //reconstructed value due to that one stencil multiplied by its weight and ready to be summed up
      weno_stencil_updates[shift] = stencil_weights[shift] * weno_stencil_recs[shift];  
    }
  }

  //sum up weighted contribution from all stencils
  interpolated_value = get_sum_of_elements( order, weno_stencil_updates );

#if( DO_LIMIT_AC_CORRECTION_NEAR_INFLECTIONS )
  //check if the updates due to different stencils have different signs, if yes - do not apply the a2c correction
  if( order == 3 && 
      (cvt_type == CVT_A2C || cvt_type == CVT_C2A ) ) 
    {
      this_order_weight = get_sum_of_elements( order, stencil_weights );
      a2cupdate = interpolated_value - uin[0] * this_order_weight;
      for( shift = 0; shift < order; shift++ ){
        a2cupdate = MINMOD( a2cupdate, (weno_stencil_recs[shift] - uin[0]) * this_order_weight ); //chooses the smallest update possible out of all possible parabolae and the full weno-5
      }
      //assign the interpolated value to the output variable
      uout[0] = uin[0] * this_order_weight + a2cupdate ;
    }
  else {
    uout[0] = interpolated_value;
  }
#else
  uout[0] = interpolated_value;
#endif


  return( 0 );
}


void c2e_simple_eno(int full_order, int is_interpolate_to_right, FTYPE *yin, FTYPE *pout)
//To get the left value, set is_interpolate_to_right = 0
{
  void simple_eno(int full_order, int cvt_type, FTYPE *yin, FTYPE *pout);
  int cvt_type;


  ///////
  // Check input values if in range
  assert( is_interpolate_to_right < 0 || is_interpolate_to_right > 1, "c2e_simple_eno(): is_interpolate_to_right can be 0 or 1\n" );
  ///////

  ///////
  // Set interpolation type: go from the center (i-th point) to left edge (i-th edge)
  cvt_type = CVT_C2L * (1 - is_interpolate_to_right) + CVT_C2R * is_interpolate_to_right;
  ///////

  simple_eno(full_order, cvt_type, yin, pout);
  

}


// performs the conversion on the full line
void compute_c2a_polycoef_simple_eno(int full_order, FTYPE *yin, FTYPE *youtpolycoef )
{

#define TWELFTH 0.083333333333333333333333333333333333L

  ///////
  // Check that type requested exists
  assert( 5 != full_order, "Only order = 5 supported by compute_c2a_polycoef_simple_eno\n" );
  ///////
  
  if( NULL == youtpolycoef ) return;

  //symmetrized expressions for derivatives
  youtpolycoef[0] = yin[0];
  youtpolycoef[1] = TWELFTH * ( (yin[-2] - yin[2]) + 8 * (- yin[-1] + yin[1]) );
  youtpolycoef[2] = TWELFTH * ( (-yin[-2] - yin[2]) + 16*(yin[-1] + yin[1])  - 30*yin[0] );
  youtpolycoef[3] = 0.5 * ( (-yin[-2] + yin[2]) + 2*(yin[-1] - yin[1]) );
  youtpolycoef[4] = ( (yin[-2] + yin[2]) - 4*(yin[-1] + yin[1]) + 6*yin[0] );

  ///////
}



//Reconstructs from the i-th point to the left to the i-th face
//For eno-5 reconstruction, use full_order = 5, etc.
void simple_eno(int full_order, int cvt_type, FTYPE *yin, FTYPE *pout)
{
  FTYPE *cvt_vec;
  FTYPE res;
  int i, is_conv, if_conv;
  FTYPE dres[MAX_CVT_ORDER];

  ///////
  // Check that type requested exists
  assert( cvt_type>NUM_CVT_TYPES, "No such cvt_type for simple_eno\n" );
  ///////


  ///////
  // Get the vector to convolve with the points that gives the reconstructed value
  cvt_vec = get_cvt_vec( cvt_type, full_order );
  ///////

  ///////
  // Convolve the vector with the values

  //initialize the variable to which you add up
  res = 0.0;
 
  // set up the start and end inidices for convolution
  is_conv = - full_order / 2;
  if_conv = is_conv + full_order - 1;

  // convolve the vector of quantities with the reconstruction vector
  for( i = is_conv; i <= if_conv; i++ ) {
    dres[i-is_conv] = cvt_vec[i-is_conv] * yin[i];
#if( !SYMMETRIZE_SIMPLE_WENO )
    res += dres[i-is_conv];
#endif
  }
  ///////

#if( SYMMETRIZE_SIMPLE_WENO )
  res = get_sum_of_elements( full_order, dres );
  //for( i = 0; i < full_order; i++ ) {
  //  res += dres[i];
  //}
#endif

  ///////
  // Return the convolution result -- the reconstructed value
  *pout = res;
  ///////
}


//Reconstructs from the i-th point average centered value to centered point value
//For eno-5 reconstruction, use full_order = 5, etc.
void a2c_simple_eno(int full_order, FTYPE *yin, FTYPE *pout)
{
  simple_eno( full_order, CVT_A2C, yin, pout );
}



//Reconstructs from the i-th centered point value to average centered value
//For eno-5 reconstruction, use full_order = 5, etc.
void c2a_simple_eno(int full_order, FTYPE *yin, FTYPE *pout)
{
  simple_eno( full_order, CVT_C2A, yin, pout );
}





FTYPE *get_cvt_vec( int cvt_type, int full_order )
{
  FTYPE *cvt_matrix;
  FTYPE *cvt_vec;
  int shift;
 
  ///////
  // Pick a matrix for the conversion order that corresponds to the full stencil
  cvt_matrix = cvt_matrices[cvt_type][full_order]; 
  ///////

  ///////
  //Pick the central stencil out of that matrix
 
  //compute the index of the central stencil
  shift = full_order / 2;

  //pick the row of the matrix that corresponds to the reconstruction due to this central stencil
  cvt_vec = &MATRIX_ROW_COL(cvt_matrix, shift, 0, full_order );
  ///////

  return( cvt_vec ); 
}






// returns optimal weights of "order" order
void c2e_simple_weno(int order, int ii, int bs, int be, FTYPE *yin, FTYPE *pleft, FTYPE *pright)
{
  FTYPE simpleweight;
  int i;
  int minindex,maxindex;
  weno_weights_t a_optimal_weights;
  weno_weights_t *optimal_weights;

  optimal_weights=&a_optimal_weights; // so can use -> still


  simpleweight = 1.0/order;

  //monotonicity indicator only implemented for order = 3
  if( order <= 5 ) {  

    //compute the weights structure with the optimal weights to pass to compute_monotonicity_indicator()
    for( i = 0; i < order; i++ ) {

      //roughness weights (inverse smoothness indicators)
      optimal_weights->weights[i] = simpleweight; 

      optimal_weights->order = order;

      // just in case, set these, GODMARK, but shouldn't be needed
      optimal_weights->len = 3 * order;
      //      optimal_weights->do_reduce = 0;
      optimal_weights->w_ratio_min = 0;
      optimal_weights->w_ratio_max = 0;
      optimal_weights->w_ratio_left = 0;
      optimal_weights->w_ratio_right = 0;


      //optimal left
      optimal_weights->weights[i + weno_weights_shifts_array[CVT_C2L] * order] = c2e_optimal_weights_leftface[order-2][i];  

      //optimal right
      optimal_weights->weights[i + weno_weights_shifts_array[CVT_C2R] * order] = c2e_optimal_weights_rightface[order-2][i];
    }
  }
  else{
    dualfprintf(fail_file,"No optimal weights for order=%d\n",order);
    myexit(16);
  }


  minindex=bs-ii;
  maxindex=be-ii;
  
  eno_cvt(CVT_C2L,minindex,maxindex,optimal_weights,yin,pleft);
  eno_cvt(CVT_C2R,minindex,maxindex,optimal_weights,yin,pright);

}






#if( DO_MONOTONICITY_INDICATOR )


#define DEBUGMONO 0

// e.g. WENO5: pr[-2,-1,0,1,2]
void compute_monotonicity_indicator(int order, int minindex, int maxindex, weno_weights_t *optimal_weights, FTYPE *pr, FTYPE *monoindicator)
{
  int i;
  FTYPE dp_a[MAXORDERS],ddp_a[MAXORDERS];
  FTYPE *dp,*ddp;
  FTYPE updp,downdp,upddp,downddp;
  FTYPE monoup,monodown;
  FTYPE pleft, pright;

  dp = (FTYPE (*) ) (&dp_a[order-2]);
  ddp = (FTYPE (*) ) (&ddp_a[order-2]);

  // dp: WENO5: i=-1 0 1 2
  for(i=-order+2;i<=order-1;i++){
    dp[i]=pr[i]-pr[i-1];
#if(DEBUGMONO)
    dualfprintf(fail_file,"dd[%d]=%21.15g\n",i,dp[i]);
#endif
  }

  // dp: WENO5: i=-1 0 1
  for(i=-order+2;i<=order-2;i++){
    ddp[i]=dp[i+1]-dp[i];
#if(DEBUGMONO)
    dualfprintf(fail_file,"ddp[%d]=%21.15g\n",i,ddp[i]);
#endif
  }



  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically increasing
  updp=(dp[-1]>=dp[0])&&(dp[0]>=dp[1])&&(dp[1]>=dp[2]);

  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically decreasing
  downdp=(dp[-1]<=dp[0])&&(dp[0]<=dp[1])&&(dp[1]<=dp[2]);

  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically curved upwards
  upddp=(ddp[-1]>=ddp[0])&&(ddp[0]>=ddp[1]);

  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically curved downwards
  downddp=(ddp[-1]<=ddp[0])&&(ddp[0]<=ddp[1]);

#if(DEBUGMONO)
  dualfprintf(fail_file,"updp=%21.15g downdp=%21.15g upddp=%21.15g downddp=%21.15g\n",updp,downdp,upddp,downddp);
#endif

  if((updp||downdp) && (upddp||downddp)){
    //  if((updp||downdp)){

    //pleft=sasha_dotproduct(order,dl,pr);
    //pright=sasha_dotproduct(order,dr,pr);

    eno_cvt(CVT_C2L,minindex,maxindex,optimal_weights,pr,&pleft);
    eno_cvt(CVT_C2R,minindex,maxindex,optimal_weights,pr,&pright);

    // check monotonicity of final interpolated values
    monoup=((pr[-1]>=pleft)&&(pleft>=pr[0])&&(pr[0]>=pright)&&(pright>=pr[1]));
    monodown=((pr[-1]<=pleft)&&(pleft<=pr[0])&&(pr[0]<=pright)&&(pright<=pr[1]));

    if(monoup||monodown){
      // Tells Sasha to set all weights to same #
      *monoindicator=1;
    }
    else{
      *monoindicator=0;
      // Sasha claims this should never occur, print out if does
      // Sasha wrong, does reach here alot.
      //    dualfprintf(fail_file,"NONMONO @ t=%21.15g\n",t);
    }

  }
  else *monoindicator=0;

#if(DEBUGMONO||1)
  // dualfprintf(fail_file,"t=%21.15g mi=%g\n",t,*monoindicator);
#endif

#if(DEBUGMONO)
  // *monoindicator=0;
#endif

}
#endif 






///A wrapper for computing different types of stencil weights.  Uses an array of pointers to functions that compute weights.  No if statements.
int compute_stencil_weights(  int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
  assert( cvt_type < 0 || cvt_type >= NUM_CVT_TYPES, "compute_stencil_weights(): cvt_type is out of bounds\n" );
  
  return( compute_stencil_weight_func_ptrs[cvt_type](cvt_type, pl, order, min_index, max_index, monoindicator0, monoindicator1, uin, stencil_weights_out) );
}

///A wrapper for computing lower-order stencil weights.  Uses an array of pointers to functions that compute weights.  No if statements.
int compute_lower_order_stencil_weights(  int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
  int lower_order;
  int ret;
  weno_weights_t weights;

  assert( cvt_type < 0 || cvt_type >= NUM_CVT_TYPES, "compute_stencil_weights(): cvt_type is out of bounds\n" );
  
  //////////////////////
  // Preparation: compute lower-order weights and copy them into dedicated lower_order_weights[] array of the weights structure
  //

  //compute lower-order unoptimized and optimized weights
  lower_order = order - 1;
  ret = compute_stencil_weight_func_ptrs[cvt_type](cvt_type, pl, lower_order, min_index, max_index, monoindicator0, monoindicator1, uin, &weights);
  ret += compute_optimized_stencil_weight_func_ptrs[cvt_type](cvt_type, pl, lower_order, min_index, max_index, monoindicator0, uin, &weights);
  extract_weights_into_array( &weights, stencil_weights_out->lower_order_weights );
  //
  //////////////////////
  
  return( ret );
}

void extract_weights_into_array( weno_weights_t *stencil_weights_in, FTYPE *weights_array_out ) 
{
  int weight_no;

  //copy lower-order weights to the dedicated array since the ->weights[] array will be overwritten by next call
  for( weight_no = 0; weight_no < stencil_weights_in->len; weight_no++ ) {
    weights_array_out[weight_no] = stencil_weights_in->weights[weight_no];
  }
}

void create_weights_from_array( int cvt_type, int order, FTYPE *weights_array_in, weno_weights_t *stencil_weights_out ) 
{
  int weight_no;

  stencil_weights_out->order = order;

  //number of weights: for a2c/c2a two sets of weights, for c2e three sets of weights
  if( cvt_type == CVT_A2C || cvt_type == CVT_C2A ) {
    stencil_weights_out->len = 2 * order;
  }
  else {
    stencil_weights_out->len = 3 * order;
  }

  //copy weights to the dedicated array
  for( weight_no = 0; weight_no < stencil_weights_out->len; weight_no++ ) {
    stencil_weights_out->weights[weight_no] = weights_array_in[weight_no];
  }

  // Reset the contents of the structure
  stencil_weights_out->w_ratio_min = 0; 
  stencil_weights_out->w_ratio_max = 0; 
  stencil_weights_out->w_ratio_left = 0; 
  stencil_weights_out->w_ratio_right = 0; 
  stencil_weights_out->w_ratio_combined = 0; 
  stencil_weights_out->lower_order_fraction = 0; 
  stencil_weights_out->lower_order_fraction_fordpp = 0;
}

///A wrapper for computing different types of stencil weights.  Uses an array of pointers to functions that compute weights.  No if statements.
int compute_optimized_stencil_weights(  int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
  assert( cvt_type < 0 || cvt_type >= NUM_CVT_TYPES, "compute_stencil_weights(): cvt_type is out of bounds\n" );

  return( compute_optimized_stencil_weight_func_ptrs[cvt_type](cvt_type, pl, order, min_index, max_index, monoindicator, uin, stencil_weights_out) );
}







/// Computes either the weights for a2c or for c2a conversion depending on whether cvt_type is CVT_A2C or CVT_C2A for the point uin[0].
/// Uses the Shi, Hu, Shu (2002) "A Technique of Treating Negative Weights in WENO Schemes"; see reconstructeno.h for more details on positive and negative weights
/// \param[in]  min_index
/// \param[out] max_index  minimum and maximum values of the indices that can be accessed in #uin array
/// \param[in]  uin  array that contains the quantities to be interpolated.  uin should point to the value around which the reconstruction takes place.
/// \param[out] stencil_weights_out  weights structure for a single point.  The contents are: roughness weights (inverse smoothness indicators) in the first (order) elements of stencil_weights; 
///          the next (order) elements contain the optimized weights.
/// 
//
int compute_ac_ca_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
  int i, min_shift, max_shift, shift; 
  FTYPE um, u0, u1, u2, u3;
  FTYPE second_der, first_der;
  FTYPE c0[3] = { 3.0, 1.0, 1.0 };
  FTYPE c1[3] = { -4.0, .0, -4.0 };
  FTYPE c2[3] = { 1.0, -1.0, 3.0 };
  FTYPE second_der_arr[3];  //for additional stencil reduction for 3rd order

  //assign the local pointers to data arrays
  FTYPE *stencil_weights = &(stencil_weights_out->weights[0]);
  FTYPE *stencil_pos_weights = &(stencil_weights_out->weights[order]);
  FTYPE *stencil_neg_weights = &(stencil_weights_out->weights[2*order]);

  FTYPE smoothness_indicators[MAX_CVT_ORDER];
  FTYPE simpleweight;

  //set the sum of smoothness indicators.  This value
  //is set to 1.0 by default so as to not influence the
  //computation of the sum of the weights if
  //normalization of the smoothness indicators is not performed
  FTYPE sum_smoothness_indicators = 1.0;  
  FTYPE sum_weights;


#if( WENO_AC_REDUCE_NEAR_CUSPS )
  stencil_weights_out->do_reduce = 0;
#endif

  stencil_weights_out->w_ratio_min = 0;
  stencil_weights_out->w_ratio_max = 0;
  stencil_weights_out->w_ratio_left = 0;
  stencil_weights_out->w_ratio_right = 0;

  assert( order < MIN_CVT_ORDER || order > MAX_CVT_ORDER, "compute_ac_stencil_weights: order out of bounds\n" );
  assert( CVT_A2C != cvt_type && CVT_C2A != cvt_type, "compute_ac_stencil_weights: cvt_type can only be CVT_A2C or CVT_C2A\n" );
  assert( max_index - min_index + 1 < order, "compute_ac_stencil_weights: not enough points to interpolate with the required order\n" );
  assert( 0 < min_index || 0 > max_index, "compute_ac_stencil_weights: [min_index, max_index] does not contain initial point" );

  //set the order of the weights that are being computed 
  stencil_weights_out->order = order;
  stencil_weights_out->len = 2 * order;  //store two types of weights: roughness weights and optimized, each of length (order)
  stencil_weights_out->lower_order_fraction = 0.; //initialize the lower order fraction with zero, will be used in compute_weights_ratios() and below in this fcn.
  

  ///////////////////////////////
  //
  // MONOINDICATOR
  //
  assert( *monoindicator0 < 0, "Monodinciator[0] = %lg\n", *monoindicator0 );
  if(*monoindicator0==1){// then equal weights
    simpleweight = 1.0/order;
    //compute the weights structure with the optimal weights to pass to compute_monotonicity_indicator()
    for( i = 0; i < order; i++ ) {
      stencil_weights_out->weights[i] = simpleweight; 
    }
    return(0);
  }
  //else if(monoindicator[0]==-1){// then 0 weight -- does not make sense SASMARK
  //   for( i = 0; i < order; i++ ) {
  //     stencil_weights_out->weights[i] = 0.0; 
  //   }
  //   return(0);
  // }
  //
  // otherwise have to compute weights
  //////////////////////////////////////////
  else if(*monoindicator0==0 && *monoindicator1!=0) { //cusp
    stencil_weights_out->lower_order_fraction = 1.0;
    //*monoindicator1 = 0.;  //reset the monoindicator and take care of mono->MINM reduction inside of WENO. SASMARK CUSPIND
  }



  if( order == 1 ){  //DONOR -- weights are not actually used
    stencil_weights[0] = stencil_weights[1] = 1.;
    return 0;
  }

  if( order == 2 ) {
    //GODMARK: this could be removed if deemed that the if statement is taking longer than the actual computation of the weights
    // for order of 2 and lower it does not matter which stencil you pick because there is no difference between cell-centered and cell-averaged value
    stencil_weights[0] = 0.5; //assign some values, they do not matter anyway due to the above
    stencil_weights[1] = 0.5;
    stencil_weights[2] = 0.5; //assign some values, they do not matter anyway due to the above
    stencil_weights[3] = 0.5;
    return( 0 );
  }


  //not enough zones for interpolation of degree == order
  //assert( max_index - min_index + 1 < order, "compute_ac_stencil_weights: not enough zones for interpolation of degree == order (max_index-min_index+1 < order)" );
  assert( order != 1 && order != 2 && order != 3, "compute_ac_stencil_weights: order can only be 1, 2 or 3" );

  min_shift = order - 1 - ( max_index ); //<= order - 1, but can be <0
  max_shift = - min_index;  //>=0, but can be > order - 1 

  //correct the above values, although, no need to correct because when uncorrected give the right bounds
  //if( min_shift < 0       ) min_shift = 0;
  //if( max_shift > order - 1 ) max_shift = order - 1;

  for( shift = 0; shift < order; shift++ ) {//cycle through all possible stencils
    if( 3 == order ) {
      //u3 =  uin[(3 - shift)];
      u2 =  uin[(2 - shift)];
      u1 =  uin[(1 - shift)];
      u0 =  uin[(0 - shift)];
      //um =  uin[(-1- shift)];

      first_der = (c0[shift] * u0 + c2[shift] * u2) + c1[shift] * u1;
      second_der =  (u0 + u2) - 2.0 * u1;

      second_der_arr[shift] = second_der;

      smoothness_indicators[shift] = 13. / 12. * second_der * second_der + 1. / 4. * first_der * first_der; //this is a smoothness indicator, \beta_r
    }
  }

#if( DO_STORE_SMOOTHNESS_INDICATORS )
  //store a copy of smoothness indicators in the weights structure
  for( shift = 0; shift < order; shift++ ) {
    stencil_weights_out->smoothness_indicators[shift] = smoothness_indicators[shift];
  }
#endif


#if( WENO_AC_REDUCE_NEAR_CUSPS )
  stencil_weights_out->do_reduce = do_reduce_near_cusps( min_index, max_index, uin );
#endif

#if( USE_NORMU )
  desensitise_smoothness_indicators( order, REDUCEEPSILON, uin, smoothness_indicators );
#endif

#if( DO_NORMALIZE_SMOOTHNESS_INDICATORS )
  /////////////////
  ///  
  ///   Normalization of smoothness indicators
  ///
  /////////////////
  for( shift = 0; shift < order; shift++ ) { //reset those smoothness indicators that stick out of the domain to zero
    if( shift < min_shift || shift > max_shift ) {  
      //the stencil would stick outside of the allowed domain -- assign zero to smoothness indicator so that it does not change norm
      smoothness_indicators[shift] = 0.0;
    }
  }

  sum_smoothness_indicators = normalize_array( order, MINNUMREPRESENT, smoothness_indicators ) + MINNUMREPRESENT;
  //now the smoothness indicators are normalized, so it makes perfect sense to use WENO_EPSILON since it is now added to a normalized quantity
#endif 

#if( FORCE_AC_CA_WEIGHT_TO_BE_OPTIMAL )
  for( shift = 0; shift < order; shift++ ) {//cycle through all possible stencils
    stencil_weights[shift] = 1./order;  //make unnormalized roughness weights to be equal to each other
  }
#else
  //computes stencil weights from the values of smoothness indicators, see (2.59) from Shu review, 1997
  compute_weights_from_smoothness_indicators( pl, order, WENO_EPSILON, uin, smoothness_indicators, stencil_weights );
#endif


  //assign zero to those weights that stick out of the domain of assigned values
  for( shift = 0; shift <= order - 1; shift++ ) {
    if( shift < min_shift || shift > max_shift ) {  
      //the stencil would stick outside of the allowed domain -- assign zero weight to it
      stencil_weights[shift] = 0.0;
    }
  }

  //normalize roughness weights
  sum_weights = normalize_array( order, 0, stencil_weights );
  
  //compute the sum of unoptimized weights (correcting for the nomeralization of smoothness indicators, if done above)
  stencil_weights_out->unoptimized_weights_sum = sum_weights / ( sum_smoothness_indicators * sum_smoothness_indicators );

#if( DO_LIMIT_AC_WEIGHTS )
  //obtain the sums of weights in a symmetric way  -- atch symmetrize
  for( shift = 0; shift < order; shift++ ) {
    //if a weight is close to the equilibrium one
    //weight limiting code operates on normalized roughness weights one by one; produces unnormalized weights
    limit_weight( order, stencil_weights[shift], &stencil_weights[shift] );
  }

  //normalize weights
  normalize_array( order, 0, stencil_weights );
#endif  

#if( DESENSITISE_STENCIL_REDUCTION )
  //add a constant to smoothness indicators which effectively makes them less sensitive; note: indicators are not normalized after this
  desensitise_smoothness_indicators( order, DESENISTISE_STENCIL_REDUCTION_EPSILON / sum_smoothness_indicators, uin, smoothness_indicators );
  compute_weights_from_smoothness_indicators( pl, order, WENO_EPSILON, uin, smoothness_indicators, stencil_weights );
  normalize_array( order, 0, stencil_weights );
#endif


  //now stencil_weights[] contains the unoptimized WENO weights, \omega_r (see Shu review (1997), p. 18 - 19)
  //optimized weno weights will be computed outside

  return( 0 );
}


//computes optimized weights assuming that stencil_weights_out array contains the untoptimized weights
int compute_optimized_ac_ca_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
  int i, min_shift, max_shift, shift; 
  FTYPE um, u0, u1, u2, u3, normdpos = 1.0, normdneg = 1.0;
  const FTYPE *dpos, *dneg;

  //assign the local pointers to data arrays
  FTYPE *stencil_weights = &(stencil_weights_out->weights[0]);
  FTYPE *stencil_pos_weights = &(stencil_weights_out->weights[order]);
  FTYPE *stencil_neg_weights = &(stencil_weights_out->weights[2*order]);

  assert( order < MIN_CVT_ORDER || order > MAX_CVT_ORDER, "compute_ac_stencil_weights: order out of bounds\n" );
  assert( CVT_A2C != cvt_type && CVT_C2A != cvt_type, "compute_ac_stencil_weights: cvt_type can only be CVT_A2C or CVT_C2A\n" );
  assert( max_index - min_index + 1 < order, "compute_ac_stencil_weights: not enough points to interpolate with the required order\n" );
  assert( 0 < min_index || 0 > max_index, "compute_ac_stencil_weights: [min_index, max_index] does not contain initial point" );

  ///////////////////////////////
  //
  // MONOINDICATOR
  //
  assert( monoindicator[0] < 0, "Monodindicator[0] = %lg\n", monoindicator[0] );
  if(monoindicator[0]==1){// then equal weights
    //unoptimized weights have already been computed, so nothing to do here
    //no need to adjust optimized weights
    return(0);
  }
  else if(monoindicator[0]==-1){// then 0 weight
    //unoptimized weights have already been computed, so nothing to do here
    //no need to adjust optimized weights
    return(0);
  }
  //
  // otherwise have to compute optimized weights
  //////////////////////////////////////////

  if( order == 1 ){  //DONOR -- weights are not actually used
    return 0;
  }

  if( order == 2 ) {
    //GODMARK: this could be removed if deemed that the if statement is taking longer than the actual computation of the weights
    // for order of 2 and lower it does not matter which stencil you pick because there is no difference between cell-centered and cell-averaged value
    return( 0 );
  }

  //assign the pointers to arrays with optimal weights for the chosen type of conversion (cvt_type) and the order (order)
  dpos = dpos_array[cvt_type][order-2];
  dneg = dneg_array[cvt_type][order-2];

  //SASMARK could be optimized to be normalized in init(), etc.
  //initialization values are already normalized -- see dpos_array & dneg_array
  for( normdpos = 0.0, normdneg = 0.0, shift = 0; shift < order; shift++ ) {
    //normalize positive and negative weights
    normdpos += dpos[shift];
    normdneg += dneg[shift];
  }

  //not enough zones for interpolation of degree == order
  //assert( max_index - min_index + 1 < order, "compute_ac_stencil_weights: not enough zones for interpolation of degree == order (max_index-min_index+1 < order)" );
  assert( order != 1 && order != 2 && order != 3, "compute_ac_stencil_weights: order can only be 1, 2 or 3" );

  //obtain the sums of weights in a symmetric way  -- atch symmetrize
  for( shift = 0; shift < order; shift++ ) {
    stencil_pos_weights[shift]= dpos[shift] * stencil_weights[shift];
    stencil_neg_weights[shift]= dneg[shift] * stencil_weights[shift];
  }

  //normalize weights  SASMARK4 : MINNUMREPRESENT or add something else?
  normalize_array( order, MINNUMREPRESENT, stencil_pos_weights );
  normalize_array( order, MINNUMREPRESENT, stencil_neg_weights );

  for( shift = 0; shift < order; shift++ ) { 
    //In notation of Shi, Hu, Shu (2002), the understandable analog of their (2.4) is:
    //\omega_i = \sigma^{+} \omega^{+} -  \sigma^{-} \omega^{-} 
    //         = \sum [ (\sigma^{+}\gamma_i^{+}/\sum\tilde\omega_i^{+} - \sigma^{-}\gamma_i^{-}/\sum\tilde\omega_i^{-}) / (\beta_i+\epsilon)^2 ].
    //where \sigma^{\pm} = \sum \gamma_j^{\pm}, 
    //\omega_i^{\pm} = \tilde \omega_i^{\pm}/ \sum_j  \tilde \omega_j^{\pm} -- normalized weights for positive (+) and negative (-) groups,
    //\tilde\omega_i^{\pm} = \gamma_i^{\pm}/(\beta_i+\epsilon)^2 -- unnormalized weights
  
    //optimized weights
    stencil_weights[shift + order] = normdpos * stencil_pos_weights[shift] - normdneg * stencil_neg_weights[shift];
  }

#if( DESENSITISE_STENCIL_REDUCTION )
  //add a constant to smoothness indicators which effectively makes them less sensitive; note: indicators are not normalized after this
  desensitise_smoothness_indicators( order, DESENISTISE_STENCIL_REDUCTION_EPSILON / sum_smoothness_indicators, uin, smoothness_indicators );
  compute_weights_from_smoothness_indicators( pl, order, WENO_EPSILON, uin, smoothness_indicators, stencil_weights );
  normalize_array( order, 0, stencil_weights );
#endif


  //now stencil_weights[] contains the WENO weights, \omega_r (see Shu review (1997), p. 18 - 19)

  return( 0 );
}







int compute_cf_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out )
//uin[0] is assumed to be the central value, so the indicies can be both negative and positive
//fills stencil_weights_out array with three sets of weight, of total length (3*order): roughness weights, leftface-optimized weights, rightface-optimized weights
{
  //variable declarations
  int i;
  int index, min_shift, max_shift, shift;
  FTYPE sum_of_nonopt_stencil_weights, norm_weights;
  int order_over_two = order / 2;  //atch symmetrize
  int is_order_odd = order % 2; //atch symmetrize
  int symmetric_shift;  //atch symmetrize
  int counter;

  //assign the local pointers to data arrays
  FTYPE *stencil_weights = &(stencil_weights_out->weights[0]);
  FTYPE *stencil_weights_leftface = &(stencil_weights_out->weights[order]);
  FTYPE *stencil_weights_rightface = &(stencil_weights_out->weights[2*order]);

  FTYPE smoothness_indicators[MAX_CVT_ORDER];

  FTYPE second_der, first_der;

  FTYPE second_der_arr[MAX_STENCIL_LENGTH], first_der_arr[MAX_STENCIL_LENGTH];

  FTYPE *u = &uin[0];  //u used in this function as a temporary pointer to one of the values in the uin[] array 

  FTYPE simpleweight;

  //set the sum of smoothness indicators.  This value
  //is set to 1.0 by default so as to not influence the
  //computation of the sum of the weights if
  //normalization of the smoothness indicators is not performed
  FTYPE sum_smoothness_indicators = 1.0;  
  FTYPE sum_weights;


  //input checks
  assert( order != 1 && order!= 2 && order != 3 && order != 5, "compute_cf_stencil_weights: order can only be 2, 3, 5" );

  assert( 0 < min_index || 0 > max_index || min_index > max_index, "choose_stencil: [min_index, max_index] does not contain initial point" );
    
  //not enough zones for interpolation of degree == order
  assert( max_index - min_index + 1 < order, "compute_cf_stencil_weights: not enough zones for interpolation of degree == order (max_index-min_index+1 < order)" );



  //inits
  stencil_weights_out->w_ratio_min = 0;
  stencil_weights_out->w_ratio_max = 0;
  stencil_weights_out->w_ratio_left = 0;
  stencil_weights_out->w_ratio_right = 0;
  stencil_weights_out->lower_order_fraction = 0.; //initialize the lower order fraction with zero, will be used in compute_weights_ratios()

  //set the order of the weights that are being computed 
  stencil_weights_out->order = order;
  stencil_weights_out->len = 3 * order;  //store three types of weights: roughness and two optimized ones (one for CVT_C2L, another for CVT_C2R), each of length (order)

  min_shift = order - 1 - ( max_index ); //<= order - 1, but can be <0
  max_shift = - min_index;  //>=0, but can be > order - 1 


  ///////////////////////////////
  //
  // MONOINDICATOR
  //
  assert( *monoindicator0 < 0, "Monodinciator[0] = %lg\n", *monoindicator0 );
  if(*monoindicator0==1){// then equal weights
    simpleweight = 1.0/order;
    //compute the weights structure with the optimal weights to pass to compute_monotonicity_indicator()
    for( i = 0; i < order; i++ ) {
      stencil_weights_out->weights[i] = simpleweight; 
    }
    return(0);
  }
  else if(*monoindicator0==0 && *monoindicator1!=0) { //cusp
    stencil_weights_out->lower_order_fraction = 1.0;
    //*monoindicator1 = 0.;  //reset the monoindicator and take care of mono->MINM reduction inside of WENO. SASMARK CUSPIND
  }
  //else if(monoindicator[0]==-1){// then 0 weight  -- does not make sense SASMARK
  //  for( i = 0; i < order; i++ ) {
  //    stencil_weights_out->weights[i] = 0.0; 
  //  }
  //  return(0);
  //}
  //
  // otherwise have to compute weights
  //////////////////////////////////////////


  //SASMARK: note that here we can read from nowhere (no checks!), however the "nowhere" info will not be used below due to checks there

  //setup smoothness indicators for a given order
  if( 1 == order ) {
    stencil_weights[0] = 1.;  //no need to calculate weights for order == 1 (there is only one stencil that contains only one grid cell, this grid cell gets all weight)
    stencil_weights_leftface[0] = 1.;
    stencil_weights_rightface[0] = 1.;
    return 0;  //the above are the weights, so no need for the below; return here.
  }
  else if( 2 == order ) {
    u = &uin[0]; //central value
    first_der_arr[0] = u[1] - u[0];
    first_der_arr[1] = u[0] - u[-1];

    for( shift = 0; shift <= order - 1; shift++ ) {//cycle through all possible stencils
      smoothness_indicators[shift] = first_der_arr[shift] * first_der_arr[shift]; //this is a smoothness indicator, \beta_r
    }

  }
  if( 3 == order ) {
    //note: these are derivatives to an accuracy of a constant factor
    u = &uin[0]; //central value
  
    first_der_arr[0] = - (3 * u[0] + u[2]) + 4 * u[1];  //atch symmetrize
    first_der_arr[1] = u[1] - u[-1];                //atch symmetrize
    first_der_arr[2] = (3 * u[0] + u[-2]) - 4 * u[-1];//atch symmetrize
    second_der_arr[0] = (u[0] + u[2]) - 2 * u[1];     //atch symmetrize
    second_der_arr[1] = (u[-1] + u[1]) - 2 * u[0];    //atch symmetrize
    second_der_arr[2] = (u[0] + u[-2]) - 2 * u[-1];   //atch symmetrize

    for( shift = 0; shift <= order - 1; shift++ ) {//cycle through all possible stencils
      smoothness_indicators[shift] = 13. / 12. * second_der_arr[shift] * second_der_arr[shift]  //SASMARK
        + 1. / 4.  * first_der_arr[shift]* first_der_arr[shift]; //this is a smoothness indicator, \beta_r
    }
  }
  else if( 5 == order ) {
    //note: these weights only work for the case of center to cell face conversion (NOT cell average to cell center, etc.)
    u = &uin[0];

    //dualfprintf( fail_file, "Test2\n" );
 
    smoothness_indicators[0] = (1337954*u[0]*u[0] + 12627689*u[1]*u[1] + 18768339*u[2]*u[2] + 5951369*u[3]*u[3] + 
                                279134*u[4]*u[4] - 21022356*u[2]*u[3] + u[1]*(-30442116*u[2] + 16810942*u[3] - 3568693*u[4]) + 
                                4503117*u[2]*u[4] - 2569471*u[3]*u[4] + u[0]*(-8055511*u[1] + 9424677*u[2] - 5121853*u[3] + 1076779*u[4])) /60480.;
    u--;
    smoothness_indicators[1] = (279134*u[0]*u[0] + 2932409*u[1]*u[1] + 4854159*u[2]*u[2] + 1650569*u[3]*u[3] + 82364*u[4]*u[4] - 
                                5550816*u[2]*u[3] + u[1]*(-7357656*u[2] + 4054702*u[3] - 847303*u[4]) + 1186167*u[2]*u[4] - 
                                725461*u[3]*u[4] + u[0]*(-1714561*u[1] + 2013987*u[2] - 1079563*u[3] + 221869*u[4]))/60480.;
    u--;
    smoothness_indicators[2] = (82364*u[0]*u[0] + 1228889*u[1]*u[1] + 2695779*u[2]*u[2] + 1228889*u[3]*u[3] + 82364*u[4]*u[4] - 
                                3495756*u[2]*u[3] + u[1]*(-3495756*u[2] + 2100862*u[3] - 461113*u[4]) + 799977*u[2]*u[4] - 
                                601771*u[3]*u[4] + u[0]*(-601771*u[1] + 799977*u[2] - 461113*u[3] + 98179*u[4]))/60480.;
    u--;
    smoothness_indicators[3] = (82364*u[0]*u[0] + 1650569*u[1]*u[1] + 4854159*u[2]*u[2] + 2932409*u[3]*u[3] + 279134*u[4]*u[4] - 
                                7357656*u[2]*u[3] + u[1]*(-5550816*u[2] + 4054702*u[3] - 1079563*u[4]) + 2013987*u[2]*u[4] - 
                                1714561*u[3]*u[4] + u[0]*(-725461*u[1] + 1186167*u[2] - 847303*u[3] + 221869*u[4]))/60480.;

    u--;
    smoothness_indicators[4] = (279134*u[0]*u[0] + 5951369*u[1]*u[1] + 18768339*u[2]*u[2] + 12627689*u[3]*u[3] + 
                                1337954*u[4]*u[4] - 30442116*u[2]*u[3] + u[1]*(-21022356*u[2] + 16810942*u[3] - 5121853*u[4]) + 
                                9424677*u[2]*u[4] - 8055511*u[3]*u[4] + u[0]*(-2569471*u[1] + 4503117*u[2] - 3568693*u[3] + 1076779*u[4]))/60480.;

    //dualfprintf( fail_file, "Test3\n" );
    u = &uin[0];

  }

#if( DO_STORE_SMOOTHNESS_INDICATORS )
  //store a copy of smoothness indicators in the weights structure
  for( shift = 0; shift < order; shift++ ) {
    stencil_weights_out->smoothness_indicators[shift] = smoothness_indicators[shift];
  }
#endif


#if( USE_NORMU )
  desensitise_smoothness_indicators( order, REDUCEEPSILON, uin, smoothness_indicators );
#endif


#if( DO_NORMALIZE_SMOOTHNESS_INDICATORS)
  /////////////////
  ///  
  ///   Normalization of smoothness indicators
  ///
  /////////////////
  for( shift = 0; shift < order; shift++ ) { //reset to zero those smoothness indicators that stick out of the domain of assigned values
    if( shift < min_shift || shift > max_shift ) {  
      //the stencil would stick outside of the allowed domain -- assign zero to smoothness indicator so that it does not change norm
      smoothness_indicators[shift] = 0.0;
    }
  }

  sum_smoothness_indicators = normalize_array( order, MINNUMREPRESENT, smoothness_indicators );

  //now the smoothness indicators are normalized, so it makes perfect sense to use WENO_EPSILON since it is now added to a normalized quantity
#endif

#if( FORCE_C2E_WEIGHT_TO_BE_OPTIMAL )
  for( shift = 0; shift < order; shift++ ) {
    stencil_weights[shift] = 1./order;  //make unnormalized roughness weights to be equal to each other
  }
#else
  //computes stencil weights from the values of smoothness indicators, see (2.59) from Shu review, 1997
  compute_weights_from_smoothness_indicators( pl, order, WENO_EPSILON, uin, smoothness_indicators, stencil_weights );
#endif


  //assign zero to those weights that stick out of the domain of assigned values
  for( shift = 0; shift <= order - 1; shift++ ) {//cycle through all possible stencils
    if( shift < min_shift || shift > max_shift ) {  
      //the stencil would stick outside of the allowed domain -- assign zero smoothness ind. to it
      stencil_weights[shift] = 0.0;
    }
  }
 

  //normalize roughness weights
  sum_weights = normalize_array( order, 0, stencil_weights );

  //compute the sum of unoptimized weights (correcting for the normalization of smoothness indicators, if done above)
  stencil_weights_out->unoptimized_weights_sum = sum_weights / ( sum_smoothness_indicators * sum_smoothness_indicators );

#if( DO_LIMIT_C2E_WEIGHTS )
  //compute optimized weights from roughness ones
  for( shift = 0; shift <= order - 1; shift++ ) { 
    //if a weight is close to the equilibrium one
    //weight limiting code operates on normalized roughness weights one by one; produces unnormalized weights
    limit_weight( order, stencil_weights[shift], &stencil_weights[shift] );
  }

  normalize_array( order, 0, stencil_weights );
  //unoptimized weights are now sitting in stencil_weights[0..order-1]; optimized weights are to be computed elsewhere
#endif  

  return( 0 );
}









int compute_optimized_cf_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out )
//uin[0] is assumed to be the central value, so the indicies can be both negative and positive
//fills stencil_weights_out array with three sets of weight, of total length (3*order): roughness weights, leftface-optimized weights, rightface-optimized weights
{
  //variable declarations
  int i;
  int index, min_shift, max_shift, shift;
  FTYPE sum_of_leftface_stencil_weights, sum_of_rightface_stencil_weights, sum_of_nonopt_stencil_weights;
  int order_over_two = order / 2;  //atch symmetrize
  int is_order_odd = order % 2; //atch symmetrize
  int symmetric_shift;  //atch symmetrize
  int counter;

  //assign the local pointers to data arrays
  FTYPE *stencil_weights = &(stencil_weights_out->weights[0]);
  FTYPE *stencil_weights_leftface = &(stencil_weights_out->weights[order]);
  FTYPE *stencil_weights_rightface = &(stencil_weights_out->weights[2*order]);

  FTYPE *u = &uin[0];  //u used in this function as a temporary pointer to one of the values in the uin[] array 

  //input checks
  assert( order != 1 && order!= 2 && order != 3 && order != 5, "compute_cf_stencil_weights: order can only be 2, 3, 5" );

  assert( 0 < min_index || 0 > max_index || min_index > max_index, "choose_stencil: [min_index, max_index] does not contain initial point" );
    
  //not enough zones for interpolation of degree == order
  assert( max_index - min_index + 1 < order, "compute_cf_stencil_weights: not enough zones for interpolation of degree == order (max_index-min_index+1 < order)" );

  //set the order of the weights that are being computed 
  stencil_weights_out->order = order;
  stencil_weights_out->len = 3 * order;  //store three types of weights: roughness and two optimized ones (one for CVT_C2L, another for CVT_C2R), each of length (order)

  ///////////////////////////////
  //
  // MONOINDICATOR
  //
  assert( monoindicator[0] < 0, "Monodindicator[0] = %lg\n", monoindicator[0] );
  if(monoindicator[0]==1){// then equal weights
    return(0);
  }
  else if(monoindicator[0]==-1){// then 0 weight
    return(0);
  }
  //
  // otherwise have to compute weights
  //////////////////////////////////////////


  //compute optimized weights
  for( shift = 0; shift <= order - 1; shift++ ) { 
    stencil_weights_leftface[shift]  = stencil_weights[shift] * c2e_optimal_weights_leftface[order-2][shift]; 
    stencil_weights_rightface[shift] = stencil_weights[shift] * c2e_optimal_weights_rightface[order-2][shift];
  }

  //SASMARK4 : is MINNUMREPRESENT ok here?
  //normalize the optimized weights
  normalize_array( order, MINNUMREPRESENT, stencil_weights_leftface );
  normalize_array( order, MINNUMREPRESENT, stencil_weights_rightface );

  //now the weights are sitting in the stencil_weights_out: (order) of unoptimized weights, (order) of optimized left, and (order) of optimized right weights

  return( 0 );
}





//////
///
///  Make smoothness indicators less sensitive to machine-level noise: add the epsilon on order of machine precision
///  This needs to be done before the smoothness indicators are normalized because that changes the norm of the quantity
///  which is used in the estimate of machine error
void desensitise_smoothness_indicators( int order, FTYPE epsilon, FTYPE *uin, FTYPE *smoothness_indicators )
{

  FTYPE normu; 
  FTYPE error_epsilon;
  int shift, i;
  FTYPE u2[MAX_CVT_ORDER];

  // Get the L_2 norm of the function in the stencil
  for( i = -order+1; i < order; i++ ){
    u2[i+order-1]= uin[i] * uin[i];
  }

  normu = get_sum_of_elements( 2 * order - 1, u2 );

  error_epsilon = normu * epsilon;

  for( shift = 0; shift < order; shift++ ) { 
    smoothness_indicators[shift] += error_epsilon;
  }
}







/// Modifies the weno weights so that the reconstructoin is closer to the one given by the total 5-point stencil in the cases where that would lead
/// to no additional oscillation.  Seems to work well for contacts by leads to a lot of oscillations for the 1D Noh problem.
/// \param[in] order
/// \param[in] u
/// \param[in,out] stencil_weights
/// \param[in,out] stencil_weights_leftface
/// \param[in,out] stencil_weights_rightface
/// \return the weight given to the full 5-point stencil
FTYPE apply_weno5_high_order_central_weight( int order, FTYPE *u, FTYPE *stencil_weights, FTYPE *stencil_weights_leftface, FTYPE *stencil_weights_rightface )
{
  int counter;
  int shift; 

  FTYPE a_sm2[3], a_sm5[3], *sm2 = a_sm2 + 1, *sm5 = a_sm5 + 1;
  FTYPE high_order_central_weight;

  assert( order != 3, "apply_weno5_high_order_central_weight() only implemented for order = 3\n" );

  //reduce to the full 5-point stencil if oscillation due to it is smaller than oscill. due to every 3-point stencil
  sm2[-1] = squareit(-u[-2] + u[0])/4. + (13*squareit(u[-2] - 2*u[-1] + u[0]))/12.;
  sm2[0] = squareit(-u[-1] + u[1])/4. + (13*squareit(u[-1] - 2*u[0] + u[1]))/12.;
  sm2[1] = squareit(-u[0] + u[2])/4. + (13*squareit(u[0] - 2*u[1] + u[2]))/12.;

  //  //note: between sm5[-1] and sm5[1] there is difference in order of the terms; otherwise they are symmetric; see reconstr_para_si3.nb
  sm5[-1] =   //(-1)-centered stencil
    squareit(-u[-2] + u[0])/16. + (21*squareit(u[-2] - 2*u[-1] + u[0]))/40. + 
    (67*squareit(11*u[-2] - 20*u[-1] + 6*u[0] + 4*u[1] - u[2]))/17280. + 
    (229*squareit(-3*u[-2] + 10*u[-1] - 12*u[0] + 6*u[1] - u[2]))/11520. + 
    squareit(-3*u[-2] - 10*u[-1] + 18*u[0] - 6*u[1] + u[2])/192. - 
    (37*squareit(u[-2] - 4*u[-1] + 6*u[0] - 4*u[1] + u[2]))/80640.;


  sm5[0] =   //(0)-centered stencil smoothness indicator using all 5 points of the total stencil
    squareit(-u[-1] + u[1])/16. + (21*squareit(u[-1] - 2*u[0] + u[1]))/40. + 
    squareit(u[-2] - 8*u[-1] + 8*u[1] - u[2])/192. + (67*squareit(-u[-2] + 16*u[-1] - 30*u[0] + 16*u[1] - u[2]))/17280. - 
    (37*squareit(u[-2] - 4*u[-1] + 6*u[0] - 4*u[1] + u[2]))/80640. + (229*squareit(-u[-2] + 2*u[-1] - 2*u[1] + u[2]))/11520.;

  sm5[1] =   //(1)-centered stencil
    squareit(-u[0] + u[2])/16. - (37*squareit(u[-2] - 4*u[-1] + 6*u[0] - 4*u[1] + u[2]))/80640. + 
    (21*squareit(u[0] - 2*u[1] + u[2]))/40. + (229*squareit(u[-2] - 6*u[-1] + 12*u[0] - 10*u[1] + 3*u[2]))/11520. + 
    squareit(-u[-2] + 6*u[-1] - 18*u[0] + 10*u[1] + 3*u[2])/192. + 
    (67*squareit(-u[-2] + 4*u[-1] + 6*u[0] - 20*u[1] + 11*u[2]))/17280.;

  high_order_central_weight = 1.;  //initialize with unity -- maximum allowed value

  for( counter = -1; counter <= 1; counter++ ) {
    sm2[counter] *= sm2[counter];
    sm5[counter] *= sm5[counter];
    high_order_central_weight = 
      MIN(
          100. * (sm2[counter] / (sm2[counter] + sm5[counter] + WENO_EPSILON ) - 0.5), 
          high_order_central_weight
          );
  }

  high_order_central_weight = MAX( high_order_central_weight, 0. );

  //make weno weights more equal so that the 5-point stencil has the above computed weight
  for( shift = 0; shift < order; shift++ ) {
    //w_new_j = w_old_j (1 - w_5) + w_5 optimized_w_j
    //this way when w_5 = 1 the scheme would use weights equal to optimized ones and will be high order
    //in this formula, the w_old_j weights have to be normalized, 0 <= w_5 <= 1 should be true and then w_new_j's are also normalized to machine precision.
    stencil_weights[shift]           *= ( 1. - high_order_central_weight );
    stencil_weights_leftface[shift]  *= ( 1. - high_order_central_weight );
    stencil_weights_rightface[shift] *= ( 1. - high_order_central_weight );

    stencil_weights[shift] += high_order_central_weight / order;  //divide by order because 1/order gives the roughness weight due to the full stencil
    stencil_weights_leftface[shift]  += high_order_central_weight * c2e_optimal_weights_leftface[order-2][shift];
    stencil_weights_rightface[shift] += high_order_central_weight * c2e_optimal_weights_rightface[order-2][shift];
  }

  return( high_order_central_weight );
}


// bad symmetry...and also how used bad symmetry : GODMARK -- needs the sign of epsilon to be flipped
//UPDATE: introduced espilon stuff that seems to be making things much better with test 666 -- oscillations
//go away almost completely
int do_reduce_near_cusps( int min_index, int max_index, FTYPE *u )
//cusp with non-monotonic behaviour aroud -- reduce to lower orders (1st or 2nd)
{
  int do_reduce = 0;
  FTYPE second_der_arr[3];
  FTYPE normu = 0.0;
  FTYPE epsilon;
  int i;

  for( i = -2; i <= 2; i++ ) {
    normu += fabs(u[i]);
  }

  epsilon = normu * WENO_EPSILON;

  if( min_index <= -2 && max_index >= 2 ) {  //only proceed if there is enough data around
    second_der_arr[0] = (u[0] + u[2]) - 2 * u[1];
    second_der_arr[1] = (u[-1] + u[1]) - 2 * u[0];
    second_der_arr[2] = (u[0] + u[-2]) - 2 * u[-1];

    if( second_der_arr[2] * ( u[0] - u[-1] ) > -epsilon &&
        second_der_arr[1] * ( u[0] - u[-1] ) < epsilon)
      {
        if( ( u[2] - u[1] ) * ( u[0] - u[-1] ) < epsilon  ) {
          do_reduce |= 2;  //to 2nd order  //atch correct reduction
        }
        else if(  ( u[1] - u[0] ) * ( u[0] - u[-1] ) < epsilon ) {
          do_reduce |= 1;  //to 1st order  //atch correct reduction
        }
      }

    if( second_der_arr[0] * ( u[0] - u[1] ) > -epsilon &&
        second_der_arr[1] * ( u[0] - u[1] ) < epsilon
        )
      {
        if( ( u[-2] - u[-1] ) * ( u[0] - u[1] ) < epsilon ) {
          do_reduce |= 2;  //to 2nd order  //atch correct reduction
        }
        else if( ( u[-1] - u[0] ) * ( u[0] - u[1] ) < epsilon ) {
          do_reduce |= 1;  //to 1st order  //atch correct reduction
        }
      }
  }

  return( do_reduce );
}







// since same as indicator in monointerp, the do_reduce parameter must be not just set incorrectly above, but ALSO used incorrectly so that symmetry is violated
int do_reduce_near_cusps_jon_version_still_fails_even_though_ok_in_monointerp_file( int min_index, int max_index, FTYPE *u )
//cusp with non-monotonic behaviour aroud -- reduce to lower orders (1st or 2nd)
{
  int do_reduce = 0;
  FTYPE second_der_arr[3];
  FTYPE a_ddf[3],*ddf,a_df[4],*df;
    
  ddf=&a_ddf[1]; // -1 0 1
  df=&a_df[1]; // -1 0 1 2

  if( min_index <= -2 && max_index >= 2 ) {  //only proceed if there is enough data around

    df[-1]=u[-1]-u[-2];
    df[0]=u[0]-u[-1];
    df[1]=u[1]-u[0];
    df[2]=u[2]-u[1];

    ddf[1] = df[2]-df[1];
    ddf[0] = df[1]-df[0];
    ddf[-1] = df[0]-df[-1];


    if(
       (sign(ddf[-1])==-sign(ddf[0]))&&
       (sign(ddf[1])==-sign(ddf[0]))&&
       (sign(ddf[-1])==sign(ddf[1]))&&
       (sign(df[-1])==-sign(df[1]))
       )
      {
        do_reduce |= 2;  //to 2nd order  //atch correct reduction
      }

    //do_reduce|=2;

  }

  return( do_reduce );
}








/// Computes stencil weights from the values of smoothness indicators, see (2.59) from Shu review, 1997.
/// \param[in] order  gives the order of the smoothness indicators (which is the same as their number)
/// \param[in] smoothness_indicators is an array with smoothness indicators, there should be order of them
/// \param[out] stencil_weights  weno weights (roughness and not normalized) computed from the smoothness indicators.  Can share memory space with  #smoothness_indicators
void compute_weights_from_smoothness_indicators( int pl, int order, FTYPE epsilon, FTYPE *uin, FTYPE *smoothness_indicators, FTYPE *stencil_weights )
{
  int shift;
  int i;
  FTYPE *temp = stencil_weights;

  FTYPE timescalefactor;
  FTYPE correctionfactor;

  //timescalefactor = coordparams.timescalefactor;

  ////This rescaling method, that brings the quantity to be interpolated to be of order of unity,
  ////does not work for interpolatin of fluxes (they are rougly conserved quantities times velocity, and
  ////the velocity bring in one more factor of timescalefactor)
  //if( pl >= U1 && pl <= U3 ) {
  // correctionfactor = timescalefactor * timescalefactor;
  //}
  //else if( pl == UU ) {
  // correctionfactor = timescalefactor * timescalefactor * timescalefactor * timescalefactor;
  //}
  //else {
  // correctionfactor = 1.;
  //}

  assert( 0 == USE_NORMU || 0 == DO_NORMALIZE_SMOOTHNESS_INDICATORS, "The epsilon used in WENO reconstruction maybe inconsistently large/small.\n" );

  //Use the usual prescription: do not rescale the epsilon since assume the smoothness indicators have been rescaled
  correctionfactor = 1.;

  for( shift = 0; shift <= order - 1; shift++ ) {  //cycle through all stencils
    //temp[shift] = MAX( smoothness_indicators[shift], epsilon / correctionfactor );
    temp[shift] = smoothness_indicators[shift] + epsilon / correctionfactor;
    temp[shift] *= temp[shift];
    stencil_weights[shift] = 1. / temp[shift];
  }
}








                      

/// \return the sum of elements in #array_to_sum_up evaluated in a symmetric fashion.  If the order of elements in the array were reversed, the result given by this function
/// would be the same.
/// \param[in] number_of_elements  number of elements in the array to sum up
/// \param[in] array_to_sum_up array containing the elements
FTYPE get_sum_of_elements( int number_of_elements, FTYPE *array_to_sum_up )
{
  FTYPE sum;
  int i, symmetric_i;
  int is_number_of_elements_odd;
  int number_of_elements_over_two;

#if( SYMMETRIZE_WENO )
  number_of_elements_over_two = number_of_elements / 2; //SASMARK could be optimized by replacing "/ 2" with ">> 1"
  is_number_of_elements_odd = number_of_elements % 2;  //< equals 1 if the number_of_elements is odd , SASMARK could be optimized by replacing "% 2" with "& 1"
  sum = 0.0;

  for( i = 0; i < number_of_elements_over_two; i++ ) { ///Loop  through array elements
    symmetric_i = number_of_elements - 1 - i;
    sum += array_to_sum_up[i] + array_to_sum_up[symmetric_i]; //pair symmetric elements when summing up
  }

  //add up the leftover un-paired element if it exists
  sum += array_to_sum_up[i] * is_number_of_elements_odd;
#else
  sum = 0.0;
  for( i = 0; i < number_of_elements; i++ ) {
    sum += array_to_sum_up[i];
  }
#endif
  return( sum );
}




FTYPE normalize_array( int number_of_elements, FTYPE epsilon, FTYPE *array_to_normalize )
{
  FTYPE sum;
  int i;
 
  sum = get_sum_of_elements( number_of_elements, array_to_normalize );

  for( i = 0; i < number_of_elements; i++ ) { ///Loop through array elements and divide each one by the sum of elements
    array_to_normalize[i] /= (sum + epsilon); ///adding small epsilon > 0 avoids division by zero if the sum == 0.0
  }

  return( sum );
}






int compute_weights_ratios(  int cvt_type, int bs, int bf, FTYPE *uin, weno_weights_t *stencil_weights_out )
{
#define NUM_POINTS 3  //the number of points at which weights are required by the stencil reduction algorithm

#define LEFT_POINT 0  //indices of points for which the weights are required by the stencil reduction algorithm
#define CENTER_POINT 1
#define RIGHT_POINT 2

  //indices of leftmost and rightmost weights
#define RIGHT_WEIGHT 0
#define LEFT_WEIGHT (order - 1)

  int order;
  int counter;
  int il;
  int ir;
  int i0;

  FTYPE w_ratio_left;
  FTYPE w_ratio_right;
  FTYPE w_ratio_combined;
  FTYPE w_ratio;
  FTYPE w_neighbour;
  FTYPE w_current;
  FTYPE min_rat;
  FTYPE max_rat;
  FTYPE lower_order_fraction;

  weno_weights_t *p_stencil_weights_array[NUM_POINTS];


  i0 = 0;  //the point index for which the weights are computed (relic from the past, now everything is shifted to be around 0)
  order = stencil_weights_out->order;

#if( DO_OR_STENCIL_REDUCTION == 1 )
  w_ratio_left = w_ratio_right = w_ratio = 0.0; //set it to some number so that does not cause reduction near the boundary
#elif(  DO_OR_STENCIL_REDUCTION == -1 )  //combined version of reduction
  w_ratio_left = w_ratio_right = w_ratio = w_ratio_combined = 0.0; //set it to some number so that does not cause reduction near the boundary
#else
  w_ratio_left = w_ratio_right = w_ratio = BIG; //set it to some number so that does not cause reduction near the boundary //atch correct
#endif 

  for( counter = 0; counter < NUM_POINTS; counter++ ) {
    //point the pointers to the arrays with precomputed weights
    //GODMARK: should not get memory violation
    p_stencil_weights_array[counter] = &stencil_weights_out[i0 + (counter - 1) * (order - 1)];  //initially assign the p_stencil_weights_array to point to precomputed weights
  }


  //get indices of two points against whose weights we compare the weights at i = i0
  il = i0 - (order - 1);
  ir = i0 + (order - 1);

  //the leftmost stencil for the i0 point ends at i = il
  if( il >= bs + order - 1 && il <= bf - (order - 1) ) {  //GODMARK: currently, the stencil at the edge will never get reduced; need to clear (order - 1)'s in this line for that to happen
    w_current = p_stencil_weights_array[CENTER_POINT]->weights[LEFT_WEIGHT];

    //this stencil is rightmost one from the point of view of the neighbour to the left
    w_neighbour = p_stencil_weights_array[LEFT_POINT]->weights[RIGHT_WEIGHT];
    w_ratio_left = ( w_current + WENO_EPSILON ) / ( w_neighbour + WENO_EPSILON );
  }
  // else {
  //  dualfprintf( fail_file, "Not enough boundary zones for reduction\n" );
  // }

  //the rightmost stencil for the i0 point ends at i = ir
  if( ir >= bs + order - 1 && ir <= bf - (order - 1) ) {  //GODMARK: currently, the stencil at the edge will never get reduced; need to clear (order - 1)'s in this line for that to happen
    w_current = p_stencil_weights_array[CENTER_POINT]->weights[RIGHT_WEIGHT];

    //this stencil is leftmost one from the point of view of the neighbour to the right
    w_neighbour = p_stencil_weights_array[RIGHT_POINT]->weights[LEFT_WEIGHT];
    w_ratio_right = (w_current + WENO_EPSILON) / ( w_neighbour + WENO_EPSILON );
  }
  // else {
  //  dualfprintf( fail_file, "Not enough boundary zones for reduction\n" );
  // }


  stencil_weights_out->w_ratio_min = MIN( w_ratio_left, w_ratio_right );
  stencil_weights_out->w_ratio_max = MAX( w_ratio_left, w_ratio_right );
  stencil_weights_out->w_ratio_left = w_ratio_left;
  stencil_weights_out->w_ratio_right = w_ratio_right;
 
  //for "(AND+OR)/2" stencil reduction: takes a combination of the weights ratios to the left and to the right; NOTE hard-cored factor of (10.)... <-- SASMARK
  stencil_weights_out->w_ratio_combined = MAX( 
                                              MIN(w_ratio_left, WENO_REDUCE_COMBINED_FACTOR / (p_stencil_weights_array[LEFT_POINT]->weights[RIGHT_WEIGHT] + WENO_EPSILON)), 
                                              MIN( WENO_REDUCE_COMBINED_FACTOR / (p_stencil_weights_array[RIGHT_POINT]->weights[LEFT_WEIGHT] + WENO_EPSILON), w_ratio_right) );


  if( DO_OR_STENCIL_REDUCTION == 1 || cvt_type == CVT_C2A || cvt_type == CVT_A2C ) {
    //"OR" stencil reduction
    w_ratio = stencil_weights_out->w_ratio_max;  
  }
  else if( DO_OR_STENCIL_REDUCTION == -1) { //combined version of reduction, no post/preshock reduction implemented yet (or maybe never)
    w_ratio = stencil_weights_out->w_ratio_combined;  

  }
  else {
    //"AND" (formerly l&r) stencil reduction
    w_ratio = stencil_weights_out->w_ratio_min;
  }

  if( cvt_type == CVT_C2E ) {  //atch modify the transition ratios for a2c/c2a to be smaller
    min_rat = MIN_TRANSITION_RATIO;
    max_rat = MAX_TRANSITION_RATIO;
  }
  else {
    //choose smaller transition ratios so that reduces easier in the a2c/c2a case
    min_rat = MIN_TRANSITION_RATIO_AC;
    max_rat = MAX_TRANSITION_RATIO_AC;
  }

  lower_order_fraction = ( w_ratio - min_rat ) / ( max_rat - min_rat );
  lower_order_fraction = MAX( 0., MIN(lower_order_fraction, 1.0) );

  stencil_weights_out->lower_order_fraction = MAX( lower_order_fraction, stencil_weights_out->lower_order_fraction );
  stencil_weights_out->lower_order_fraction_fordpp = stencil_weights_out->lower_order_fraction; //used for POST/PRESHOK reduction to avoid order-dependent effects
 

  return( 0 );
}







//only weights for cells for i = (bs + order - 1), ..., (bf - order + 1) will be used.
//currently, the stencil is not reduced for the cells sufficiently close to the boundaries: with i = bs, ..., bs + 2 * order - 3, bf - 2 * order + 3, ..., bf.
//GODMARK: this function is an exact copy of choose_ac_ca_weno_order() with the call to compute_ac_stencil_weights being replaced by compute_cf_stencil_weights()
int choose_weno_order( int cvt_type, int whichreduce, int max_order, int min_order, int i0, int pl, int bs, int bf, FTYPE (*shockindicator)[NBIGM], FTYPE (*dP)[NBIGM],
                       FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *Pindicator, FTYPE *uin, 
                       weno_weights_t *stencil_weights_array, 
                       weno_weights_t **pp_stencil_weights_to_be_used, struct of_trueijkp *trueijkp ) 
{
#define NUM_POINTS 3  //the number of points at which weights are required by the stencil reduction algorithm

#define LEFT_POINT 0  //indices of points for which the weights are required by the stencil reduction algorithm
#define CENTER_POINT 1
#define RIGHT_POINT 2

  //indices of leftmost and rightmost weights
#define RIGHT_WEIGHT 0
#define LEFT_WEIGHT (order - 1)

  int order, lower_order;
  int num_weights = 1;  //default value for the number of sets of weights that is returned; currently can be 1 or 2.
  int i, counter, il, ir;
  weno_weights_t stencil_weights_array_static_other[NUM_POINTS]; // OPENMPMARK: Can't leave as static
  weno_weights_t *p_stencil_weights_array_other[NUM_POINTS];  //pointers to actual arrays with weights; used as an intermediary layer
  FTYPE lower_order_fraction;
  FTYPE dPi;
#if( DOENODEBUG )
  FTYPE lower_order_fraction_old;  //atch debug
#endif 
  int di, dj, dk;
  //  ITERGLOBALDEF  //define values of di, dj, dk, etc.
  di = (trueijkp->iter==1); dj = (trueijkp->iter==2); dk = (trueijkp->iter==3);

  //recompute the weights and store them locally in static arrays and
  //re-direct the pointers to point to these local arrays so that they are used
  //GODMARK: this is the only line different from choose_ac_weno_order()

  for( counter = 0; counter < NUM_POINTS; counter++ ) {
    //point the pointers to the arrays with precomputed weights
    //GODMARK: should not get memory violation
    p_stencil_weights_array_other[counter] = &stencil_weights_array[i0 + (counter - 1) * (max_order - 1)];  //initially assign the p_stencil_weights_array_other to point to precomputed weights
  }

#if( WENO_C2E_REDUCE_NEAR_CUSPS || WENO_AC_REDUCE_NEAR_CUSPS)
  if( (stencil_weights_array[i0].do_reduce & 1) != 0) {
    //set the right order
    stencil_weights_array_static_other[CENTER_POINT].do_reduce = 1; //indicates that the reduction to the 1st order was required; probably unnecessary, just set it to be something
    stencil_weights_array_static_other[CENTER_POINT].order = 1;
    stencil_weights_array_static_other[CENTER_POINT].len = 3;
    for( counter = 0; counter < stencil_weights_array_static_other[CENTER_POINT].len; counter++ ) {
      stencil_weights_array_static_other[CENTER_POINT].weights[counter] = 1.;  //set optimized, left and right weights to unity (because reducing to DONOR)
    }

    num_weights = 1;  //reduce all the way to lowest order = 1, therefore no contributions of other orders at all

    //set the pointer to weight arrays: structures that hold weights of two orders follow each other in memory
    p_stencil_weights_array_other[CENTER_POINT] = &stencil_weights_array_static_other[CENTER_POINT];

    *pp_stencil_weights_to_be_used = p_stencil_weights_array_other[CENTER_POINT];  //make sure that in all cases this is the actual pointer to where the weights are

    return( num_weights );
  }
#endif

#if(0)
  if( whichreduce == WENO_REDUCE_TYPE_PPM ) {
    ///SASMARK: NOT KEPT UP TO DATE!!!

    //do PPM reduction instead of default SASHA-type one
    //if WENO_USE_PPM_FLATTENING is enabled, it overrides all other types of flattening/stencil reduction
    lower_order_fraction = shockindicator[EOMSETMHD][i0];

    assert( lower_order_fraction > 1. || lower_order_fraction < 0., "choose_cf_weno_order: shockindicator value is out of bounds [0., 1.]\n" );

    //copy the data over to local static storage for manipulation
    //BE CAREFUL:  make sure this is the first write operation so as not to lose the values of weights for the current order
    //(p_stencil_weights_array_other[CENTER_POINT] may point to the same memory space as stencil_weights_array_static_other[1])
    stencil_weights_array_static_other[CENTER_POINT] = *(p_stencil_weights_array_other[CENTER_POINT]);

    //rescale the weights of current order for output; we have three sets of weights, each of order elements
    for( counter = 0; counter < stencil_weights_array_static_other[CENTER_POINT].len; counter++ ) {
      stencil_weights_array_static_other[CENTER_POINT].weights[counter] *= ( 1. - lower_order_fraction );
    }

    //use a linear mix with 1st order as in PPM flattening
    lower_order = 1; 

    //check for a memory leak
    assert( CENTER_POINT + 1 >= NUM_POINTS,  
            "choose_weno_order: internal error: CENTER_POINT + 1 >= NUM_POINTS, we assumed that there is additional memory space beyond CENTER_POINT.\n" );

    //compute lower-order weights and place them in the next array element after the weights of current order
    compute_stencil_weights( cvt_type, pl, lower_order, bs - i0, bf - i0, &monoindicator0[i0], &monoindicator1[i0], &uin[i0], &stencil_weights_array_static_other[CENTER_POINT+1] );
    compute_optimized_stencil_weights( cvt_type, pl, lower_order, bs - i0, bf - i0, &monoindicator0[i0], &uin[i0], &stencil_weights_array_static_other[CENTER_POINT+1] );

    //rescale lower-order weights so that weights of current order and lower order now add up to unity
    for( counter = 0; counter < stencil_weights_array_static_other[CENTER_POINT+1].len; counter++ ) {
      stencil_weights_array_static_other[CENTER_POINT+1].weights[counter] *= lower_order_fraction;
    }

    //set the pointer to weight arrays: structures that hold weights of two orders follow each other in memory
    p_stencil_weights_array_other[CENTER_POINT] = &stencil_weights_array_static_other[CENTER_POINT];

    *pp_stencil_weights_to_be_used = p_stencil_weights_array_other[CENTER_POINT];  //make sure that in all cases this is the actual pointer to where the weights are

    //the number of sets of weights that is output
    num_weights = 2;  //different from the default value

    return( num_weights );
  }
#endif

  lower_order_fraction = stencil_weights_array[i0].lower_order_fraction;
 
  for( order = max_order; order > min_order; order-- ) {
    
#if( WENO_C2E_REDUCE_NEAR_CUSPS == 1 || WENO_AC_REDUCE_NEAR_CUSPS == 1)
    if( stencil_weights_array[i0].do_reduce == 2 ) {
      w_ratio = max_rat;
    }
#endif

    if( 0.0 < lower_order_fraction ) {   //atch modified to always have two orders
      //the lower order will have the following weight, it is always between 0 and 1:
   
      //copy the data over to local static storage for manipulation
      //BE CAREFUL:  make sure this is the first write operation so as not to lose the values of weights for the current order
      //(p_stencil_weights_array_other[CENTER_POINT] may point to the same memory space as stencil_weights_array_static_other[1])
      stencil_weights_array_static_other[CENTER_POINT] = *(p_stencil_weights_array_other[CENTER_POINT]);

      //rescale the weno5 weights; we have three sets of weights, each of order elements = len elements
      //start from counter = order so that we rescale only the optimized weights; we do not have to rescale the unoptimized weights since
      //they are not used in getting the reconstructed values
      for( counter = order; counter < stencil_weights_array[i0].len; counter++ ) {
        stencil_weights_array_static_other[CENTER_POINT].weights[counter] *= ( 1. - lower_order_fraction );
      }

      //set the right order
      stencil_weights_array_static_other[CENTER_POINT].order = order;

      //the ratio is in the transition region -- use a linear combination of two reconstructions
      //of the current order and the lower order:
      lower_order = order - 1; 

      //check for a memory leak
      assert( CENTER_POINT + 1 >= NUM_POINTS,  
              "choose_weno_order: internal error: CENTER_POINT + 1 >= NUM_POINTS, we assumed that there is additional memory space beyond CENTER_POINT.\n" );

      //NO NEED TO DO COMPUTE WEIGHTS: WEIGHTS HAVE ALREADY BEEN COMPUTED.  INSTEAD, PULL THEM OUT OF EXISTING STRUCTURE (SASMARK)
      //-KILLED-compute lower-order weights and place them in the next array element after the weights of current order
#if(0) // old method of computing weno3 weights in RECON_CALC
      compute_stencil_weights( cvt_type, pl, lower_order, bs - i0, bf - i0, &monoindicator0[i0], &monoindicator1[i0], &uin[i0], &stencil_weights_array_static_other[CENTER_POINT+1] );
      compute_optimized_stencil_weights( cvt_type, pl, lower_order, bs - i0, bf - i0, &monoindicator0[i0], &uin[i0], &stencil_weights_array_static_other[CENTER_POINT+1] );
#else
      // new (current) method of computing weno3 weights in WEIGHT_CALC and pulling them out for use here
      create_weights_from_array( cvt_type, lower_order, stencil_weights_array_static_other[CENTER_POINT].lower_order_weights, &stencil_weights_array_static_other[CENTER_POINT+1] );
#endif

      //rescale lower-order weights so that weights of current order and lower order now add up to unity
      for( counter = 0; counter < stencil_weights_array_static_other[CENTER_POINT+1].len; counter++ ) {
        stencil_weights_array_static_other[CENTER_POINT+1].weights[counter] *= lower_order_fraction;
      }

      //set the pointer to weight arrays: structures that hold weights of two orders follow each other in memory
      p_stencil_weights_array_other[CENTER_POINT] = &stencil_weights_array_static_other[CENTER_POINT];

      //the number of sets of weights that is output
      num_weights = 2;  //different from the default value

    } 

    //if reached here, do not reduce the stencil
    break;
  }  //end for( order = ... )

  *pp_stencil_weights_to_be_used = p_stencil_weights_array_other[CENTER_POINT];  //make sure that in all cases this is the actual pointer to where the weights are

  return( num_weights );
}






FTYPE compute_lower_order_fraction( FTYPE w_ratio_left, FTYPE w_ratio_right ) 
{
  FTYPE w_ratio;
  FTYPE w_ratio_sq;
  const FTYPE w_transition_ratio_sq = 100.;
  const FTYPE w_ratio_min = MIN_TRANSITION_RATIO;
  const FTYPE w_ratio_max = MAX_TRANSITION_RATIO;
  FTYPE lower_order_fraction;
  FTYPE a, b, c, d;
  FTYPE x, x1, x2, x_2, x_3;
  //geometrical mean -- smooth replacement for the w_ratio = MIN( w_ratio_left, w_ratio_right )
  //w_ratio = w_ratio_left * w_ratio_right / ( w_ratio_left + w_ratio_right );
  w_ratio = MIN( w_ratio_left, w_ratio_right );
 
  if( w_ratio <= w_ratio_min ) {
    lower_order_fraction = 0.0;
  }
  else if( w_ratio >= w_ratio_max ) { 
    lower_order_fraction = 1.0;
  }
  else {
    //x = w_ratio; 
    //x_2 = x * x;
    //x_3 = x * x_2;

    //x1 = w_ratio_min;
    //x2 = w_ratio_max;

    //a = (x1 - 3*x2)/(pow(x2,2)*pow(-x1 + x2,3));
    //b = (-2*(pow(x1,2) - 2*x1*x2 - 2*pow(x2,2)))/(pow(x2,2)*pow(-x1 + x2,3));
    //c = -((x1*(pow(x1,2) + x1*x2 - 8*pow(x2,2)))/(pow(x1 - x2,3)*pow(x2,2)));
    //d = (2*(pow(x1,3) - 2*pow(x1,2)*x2))/(pow(x1 - x2,3)*x2);

    //lower_order_fraction = x * ( a * x_3 + b * x_2 + c * x + d );
    lower_order_fraction = ( w_ratio - MIN_TRANSITION_RATIO ) / ( MAX_TRANSITION_RATIO - MIN_TRANSITION_RATIO );
  }

  return( lower_order_fraction );
}


#include "para_and_paraenohybrid.h"


// Pass 1D line to ENO scheme
void pass_1d_line_weno(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp)
{
  void pass_1d_line_weno_withweights(int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);
  weno_weights_t *stencil_weights_array; // shouldn't be used

  stencil_weights_array = NULL;

  pass_1d_line_weno_withweights( whichquantity, dir, do_weight_or_recon, stencil_weights_array, recontype, whichreduce, preforder, 
                                 pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim, ystencilvar, 
                                 yin, yout,youtpolycoef,trueijkp);

}


// Pass 1D line to ENO scheme
void pass_1d_line_weno_withweights(int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp)
{
  FTYPE (*ysend)[NBIGM];

  // These are ENO functions
  //extern int eno_line_c2e(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*shockindicator)[NBIGM], FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, struct of_trueijkp *trueijkp );
  //extern int eno_line_a2c(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*shockindicator)[NBIGM], FTYPE *yin, FTYPE *yout, struct of_trueijkp *trueijkp );
  //extern int eno_line_c2a(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*shockindicator)[NBIGM], FTYPE *yin, FTYPE *yout, struct of_trueijkp *trueijkp );
  //extern int eno_line_a2em(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*shockindicator)[NBIGM], FTYPE *yin, FTYPE *yout, struct of_trueijkp *trueijkp );
  //extern int eno_line_a2ep(int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift, FTYPE (*shockindicator)[NBIGM], FTYPE *yin, FTYPE *yout, struct of_trueijkp *trueijkp );

  // pass yin so that yin[0] is first active point


  /////////////////////////////////
  //
  // Choose which type of quantity to send
  //
  /////////////////////////////////
  if(do_weight_or_recon==WEIGHT_CALC){
    ysend = ystencilvar;
    //ysend = yin;
  }
  else if(do_weight_or_recon==RECON_CALC || do_weight_or_recon==ALL_CALC){
    ysend = yin;
  }




  // DEATHMARK : Sasha, these are your functions
  if(recontype==CVT_C2E){
    // yout filled correctly as 1 quantity interpolated to edges
    //      eno_line_c2e(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,yin[0]+bs,yout,trueijkp);
    // returns left primitive from 0 .. N and right primitive from -1 .. N-1
    // i.e. pleft(0..N) and pright(-1..N-1) as defined below
    // so inclusive on pleft/right need -1..N
    //////////////////////////////////////
    //
    // interpolate primitive using slope
    //
    // |=interface
    // i=zone center of ith zone
    //
    // |         pl(i)|pr(i)    i          |
    // |         Fl(i)|Fr(i)    i          |
    // |         Ul(i)|Ur(i)    i          |
    // |              |pleft(i)   pright(i)|
    // |              |F(i)                |
    //
    //
    //
    //////////////////////////////////////
    //atch: do not shift yin[] -- otherwise need to shift all other indices; spit out two arrays for left and right
    //eno_line_c2e(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,ysend[0], yout[0], yout[1],trueijkp);  

    //GODMARK: pass maxorder and minorder as arrays for each point of interest
    //  eno_line_c2e( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,  shockindicator, df, dP, etai, monoindicator, yprim[UU][0], ysend[0], yout[0], yout[1],trueijkp);  
    if(PARAMODWENO==0){
      eno_line_c2e( whichquantity, dir, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,  shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim[VSQ][0], ysend[0], yout[0], yout[1], youtpolycoef,trueijkp );
    }
    else{
      paraenohybrid_line_c2e( whichquantity, dir, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift,  shockindicator, stiffindicator, Vline, Pline, df, dP, etai, monoindicator, yprim[VSQ][0], ysend[0], yout[0], yout[1], youtpolycoef,trueijkp );
    }
  }
  else if(recontype==CVT_A2C){
    // yout filled correctly as 1 quantity "deconstructed" to 1 quantity
    //  eno_line_a2c(whichreduce,maxorder,minorder,bs,ps,pe,be,shift,ysend[0],yout[0],trueijkp);
    // for ENOFLUXRECONTPE: assume previously bounded fluxes, then returns fluxes from 0 .. N (not just N-1) ( i.e. F(0..N) as defined below)
    // for ENOAVG2CENTTYPE: assume previously bounded average U's from 0-order/2 .. N-1+order/2, then return U's from 0 .. N-1
    //////////////////////////////////////
    //
    // interpolate primitive using slope
    //
    // |=interface
    // i=zone center of ith zone
    //
    // |         pl(i)|pr(i)    i          |
    // |         Fl(i)|Fr(i)    i          |
    // |         Ul(i)|Ur(i)    i          |
    // |              |F(i)                |
    //
    //
    //
    //////////////////////////////////////

#if(0) // DEBUG TEST
    // no i2,j2, etc. present anymore
    //      yiniter=bs;
    //      SUPERGENLOOP(i2,j2,k2,is2,ie2,js2,je2,ks2,ke2,di2,dj2,dk2){
    //
    // yout[0][yiniter]=yin[0][yiniter];
    // yiniter++;
    //      }
#endif
    //  eno_line_a2c( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[UU][0], ysend[0], yout[0],trueijkp ); 
    eno_line_a2c( whichquantity, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[VSQ][0], ysend[0], yout[0],trueijkp ); 
  }
  // if(interporflux==ENOFLUXSPLITTYPE){
  // yout filled as each left/right prepared quantity is interpolated as an average to the left/right edges
  // eno_line_a2em(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,ysend[0],yout[0],trueijkp);
  // eno_line_a2ep(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,ysend[1],yout[1],trueijkp);
  // GODMARK: should return left fluxes from 0 .. N and right fluxes from -1 .. N-1
  // i.e. Fleft(0..N) and Fright(-1..N-1) as defined below
  // so inclusive on Fleft/right need -1..N
  //////////////////////////////////////
  //
  // convert point fluxes to quasi-de-averaged point fluxes at cell interface
  //
  // |=interface
  // i=zone center of ith zone
  //
  // |         pl(i)|pr(i)    i          |
  // |         Fl(i)|Fr(i)    i          |
  // |         Ul(i)|Ur(i)    i          |
  // |              |        F+-(i)      |
  // |              |Fleft(i)   Fright(i)|
  // |              |F(i)                |
  //
  //
  //
  //////////////////////////////////////
  //  }
  else if(recontype==CVT_C2A){
    // yout filled as each left/right prepared quantity is interpolated as an average to the left/right edges
    // eno_line_c2a(1,preforder,bs,ps,pe,be,minorder,maxorder,shift,ysend[0],yout[0],trueijkp);
    // For F1:
    // INPUT  : Inputs fluxes from j=-N2BND .. N2-1+N2BND inclusive for F1(i,j) for a single i
    // OUTPUT : Returns fluxes from j=0 .. N2-1 for F1(i,j) for a single i (i.e. F1(i,j=0..N2-1)
    // where for F1 interporflux==ENOFLUXAVG1TYPE doesn't apply, ENOFLUXAVG2TYPE involves averaging over second direction and ENOFLUXAVG3TYPE involves averaging over third direction
    // the input data is controlled by the loops above
    //////////////////////////////////////
    //
    // convert point flux to surface integrated flux (this call does only 1 of the NDIMEN-1 directions)
    //
    // |=interface
    // i=zone center of ith zone
    // |              |                    |
    // |            F(i)        i          |
    // |              |                    |
    // |              |                    |
    // |              |                    |
    //
    //
    //
    //////////////////////////////////////
    //  eno_line_c2a( whichquantity, do_weight_or_recon, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[UU][0], ysend[0], yout[0],trueijkp );   //atch correct
    eno_line_c2a( whichquantity, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, df, dP, monoindicator, yprim[VSQ][0], ysend[0], yout[0],trueijkp);   //atch correct
  }


}




void multidir_pre_slope_lim_linetype_weno(void)
{
  int i,j,k,pl,pliter;


#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  //tell to SWENO that should look at other directions when performing reduction
  do_weno_lower_order_fraction = 1;
  FULLLOOP PALLLOOP(pl) {  //atch zero out the fractions so that initially use WENO5 everywhere
    GLOBALMACP0A1(weno_lower_order_fraction,i,j,k,pl) = 0.0;
  }
#endif




}





void multidir_post_slope_lim_linetype_weno(void)
{

#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  do_weno_lower_order_fraction = 0;  //do not look at the other directions when performing reduction below this point
#endif

}


#define WHICHETAIPL RHO // temporary test GODMARK SUPERGODMARK
//#define WHICHETAIPL pl // temporary test GODMARK SUPERGODMARK




void pass_1d_line_multipl_weno(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp)
{
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  int myi;
  int is_copy;
  int domultipl;
  int i;
  extern void compute_polycoef_line(
                                    int preforder, int pl, int bs, int ps, int pe, int be, 
                                    FTYPE *yin, FTYPE (*youtpolycoef)[NBIGM]);
  // OPENMPMARK: Can't leave as static
  weno_weights_t a_stencil_weights_array_allpl[MAXNPR][NBIGM];  //for storing the weights for a one-dimensional array of values; used by stencil reduction
  weno_weights_t (*stencil_weights_array_allpl)[NBIGM] = (weno_weights_t (*)[NBIGM]) (&(a_stencil_weights_array_allpl[0][NBIGBND])); //GODMARK:  what should I use NBIGBND or sth else?







  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  // default
  domultipl=0;



  if( (recontype == CVT_A2C || recontype == CVT_C2A) ){
    if(MULTIPLTYPE == NOSPLITA2C) { //do not do any weights minimization for source term integration
      domultipl=0;
      //unsplit version, compute weights and reconstruction in one call
      NUMPRIMLOOP(pliter,pl) pass_1d_line_weno_withweights( whichquantity, dir, ALL_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[WHICHETAIPL], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl],youtpolycoef[pl],trueijkp);
      //SUPERSASMARK : insert a2c/c2a limiting code here
      return;
    }
    else if(MULTIPLTYPE != NOSPLITA2C) {
      domultipl=1;
      // then not done, use below multipl code
    }
    else {
      dualfprintf( fail_file, "Unknown reconstruction split type: MULTIPLTYPE = %d\n", MULTIPLTYPE );
      myexit( 1 );
    }
  }
  else{
    domultipl=0;

    // GODMARK: for c2e, currently no method for using DO_SPLITA2C method.  Could compute conserved quantities and use those for weights for primitive interpolation?

    //unsplit version, compute weights and reconstruction in one call
    NUMPRIMLOOP(pliter,pl) pass_1d_line_weno_withweights( whichquantity, dir, ALL_CALC, stencil_weights_array_allpl[pl], recontype, whichreduce, preforder, pl, bs, ps, pe, be, minorder, maxorder, shift, shockindicator, stiffindicator, Vline, Pline, df[pl], dP, etai[WHICHETAIPL], monoindicator[pl], yprim, ystencilvar[pl], yin[pl], yout[pl],youtpolycoef[pl],trueijkp);

#if(0)  
    //1d test of computing the polynomial coefficients

    // This code was used to verify that WENO and SMONO return correct expressions for 
    // polynomial derivatives. Sasha set up a function that at i = 0 has all derivatives (0th to 4th)
    // equal to 1.  The below code reproduces this, i.e. youtpolycoef[0][0..4][0] = 1 up to machine
    // precision.

    //NOTE: set FORCE_C2E_WEIGHT_TO_BE_OPTIMAL to 1 in reconstructeno.h for the purpose of the test.
    //      This resets the unoptimized weights to be equal.  This test then shows that if the unoptimized
    //      weights are equal, the derivatives of 4th order polynomials are computed exactly up to 4th derivative.
    
    pl = 0; 
    for( i = bs; i < be; i++ ) {
      yin[pl][1][i] = yin[pl][0][i] = 1 + i + pow(i,2)/2. + pow(i,3)/6. + pow(i,4)/24. + 0 * pow(i,5)/120.; 
      ystencilvar[pl][1][i] = ystencilvar[pl][0][i] = yin[pl][0][i];
      monoindicator[pl][0][i] = 0.0;
      monoindicator[pl][1][i] = 0.0;
      monoindicator[pl][2][i] = 0.0;
    }

    // test weno
    pass_1d_line_weno_withweights( whichquantity, 
                                   dir, ALL_CALC, 
                                   stencil_weights_array_allpl[pl], 
                                   recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                                   minorder, maxorder, shift, shockindicator, 
                                   stiffindicator, Vline, Pline, df[pl], 
                                   dP, etai[WHICHETAIPL], 
                                   monoindicator[pl], 
                                   yprim, 
                                   ystencilvar[pl], 
                                   yin[pl], 
                                   yout[pl],
                                   youtpolycoef[pl],trueijkp);
    
    // compare youtpolycoef[0][0..4][0] against 1

    //recompute the same coefficients using simple weno
    compute_polycoef_line( 5, pl, bs, ps, pe, be, 
                           yin[0][pl], youtpolycoef[pl] );

    // compare youtpolycoef[0][0..4][0] against 1
#endif

    return;
  }






  ///////////////
  //
  // Doing multipl version (which right now is used to split weight and value computations)
  //
  ///////////////

  if(domultipl){
    compute_multipl_weno(MULTIPLTYPE, whichquantity, dir, do_weight_or_recon, stencil_weights_array_allpl, recontype, whichreduce, preforder, bs, ps, pe, be, minorder, maxorder, shift,   shockindicator, stiffindicator, Vline,  Pline, df, dP, etai, monoindicator, yprim, ystencilvar, yin, yout, youtpolycoef,trueijkp);
  }


}





/// Limits the weights so that if they are sufficiently close to optimal ones, they will be
/// adjusted to be optimal; if they are close to zero, they will be made zero
/// NOTE: old_weights should be normalized unoptimal weights
static void limit_weight( int order, FTYPE old_weight, FTYPE *new_weight ) 
{
  FTYPE w0, w1, w2, w3, avg_weight, w4;
 
  avg_weight = 1. / order;
  w0 = 0.0;
  w1 = 0.2 * avg_weight;
  w2 = (1.0-WEIGHTSPREADAROUNDEQUAL) * avg_weight;
  w3 = (1.0+WEIGHTSPREADAROUNDEQUAL) * avg_weight;
  w4 = 1.0;

  if( old_weight < w1 ) {
    *new_weight = w0;
  }
  else if( old_weight < w2 ) {
    *new_weight = avg_weight * ( old_weight - w1 ) / ( w2 - w1 );
  }
  else if( old_weight < w3 ) {
    *new_weight = avg_weight;
  }
  else {
    *new_weight = ( old_weight - w3 ) * ( w4 - avg_weight ) / ( w4 - w3 ) + avg_weight;
  }
}

FTYPE squareit( FTYPE base )
{
  return( base * base );
}

FTYPE limit_ac_correction( int order, int pl, int bs, int bf, FTYPE max_frac_difference, FTYPE *yin, FTYPE *yout ) 
{
  FTYPE norm;
  FTYPE ac_correction;
  FTYPE relative_correction;
  FTYPE fraction_point_value;
  FTYPE start_correction;
  FTYPE end_correction;
  int is;
  int ie;
  int i;

  //initialize the difference with zero so that it is smaller than any positive number
  norm = MINNUMREPRESENT;

  //find the norm of the function
  //for( i = is; i <= ie; i++ ) {
  //  diff += fabs(yin[i]);
  //}

  norm += fabs(yin[0]) + fabs(yout[0]);

  ac_correction = yout[0] - yin[0];
  
  //relative correction is defined to be non-negative and normalized by max_frac_difference
  relative_correction = fabs(ac_correction) / ( max_frac_difference * norm );
  
  //continuous triangular-shaped depenendence, maximal (=1) at relative_correction = 0.5; minimal (=0) at relative_correction <= 0 & >= 1.
  fraction_point_value =   MAX( MIN(2. - 2. * relative_correction, 1.0), 0. );   

      
  // apply correction
  yout[0] = yin[0] + ac_correction * fraction_point_value;

  return( fraction_point_value );
}










// included in file instead as separate file so static memory spaces global to this file are global to those files but still private to other files
#include "reconstructeno.weightmin.c"
#include "reconstructeno.debug.c"












