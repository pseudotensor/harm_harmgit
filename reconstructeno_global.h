#ifndef RECONSTRUCTENO_GLOBAL_H
#define RECONSTRUCTENO_GLOBAL_H


#define MAX_NO_OF_WENO_OUTPUTS 2



#define DO_ENO_STENCIL_REDUCTION 1  //atch

#define DO_NORMALIZE_SMOOTHNESS_INDICATORS 1 //atch normbeta
#define DO_RESCALE_WEIGHTS 0

#define USE_NORMU 1 // use WENO epsilon that is proportional to the value of the function

// whether to store smoothness indicators in weight array
#define DO_STORE_SMOOTHNESS_INDICATORS 1

// c2e stuff
#define DO_LIMIT_C2E_CORRECTION 0
#define WENO_C2E_MAX_FRAC_DIFFERENCE (0.5)
#define DO_CENT2EDGE 1 // whether to do c2e; if set to 0, edge values are equal to the center value
#define WENO_C2E_REDUCE_NEAR_CUSPS 0 // let mono handle this

#define COUNTCALLS 0

// a2c/c2a stuff
#define DO_AVG2CEN 1 // whether to do a2c
#define DO_CEN2AVG 1 // whether to do c2a
#define WENO_AC_REDUCE_NEAR_CUSPS 0 // let mono handle this
#define DO_LIMIT_AC_CORRECTION 1
#define WENO_AC_MAX_FRAC_DIFFERENCE (0.5)  //maximum fractional change allowed during the a2c/c2a conversion; should be >= 0 and < 1
#define WENO_DO_FLIP_CONS_SIGN 1
#define WENO_DIR_FLIP_CONS_SIGN_DN ( (ISSPCMCOORDNATIVE(MCOORD)&&WENO_DO_FLIP_CONS_SIGN)?(2):( (MCOORD == CYLMINKMETRIC&&WENO_DO_FLIP_CONS_SIGN)?(1):(0) ) )
#define WENO_DIR_FLIP_CONS_SIGN_UP ( (ISSPCMCOORDNATIVE(MCOORD)&&WENO_DO_FLIP_CONS_SIGN)?(2):( (MCOORD == CYLMINKMETRIC&&WENO_DO_FLIP_CONS_SIGN)?(0):(0) ) )


#define SYMMETRIZE_SIMPLE_WENO 1
#define SYMMETRIZE_WENO 1
#define SSYMCHECKS  0 //atch symmetry checks flag

// JCM GODMARK:below AND causes lack of resolution for high frequency waves that are totally normal
#define DO_OR_STENCIL_REDUCTION (0)  //if equals 1, uses the "OR" version of stencil reduction, if 0 uses the "AND" (formerly l&r version), if (-1) -- the combined (AND+OR)/2 version
#define DO_OR_STENCIL_REDUCTION_REDUCE_AC_EXTRA_POINT (0)  //helps to keep the hard test from RAM paper much more stable by reducing the a2c more, only appies when DO_OR_STENCIL_REDUCTION == 1
#define WENO_REDUCE_COMBINED_FACTOR (0.1)  //only applies for DO_OR_STENCIL_REDUCTION = -1









//whether to smoothly rescale the SWENO weights in accordance with smoothness indicators to allow for smooth change in the behaviour of SWENO stencil reduction
#define DO_SMOOTH_ADJUSTMENT_UNOPTIMIZED_WEIGHTS 1
#define SMONO_DO_SMOOTH_TRANSITION 1
#define SMONO_EPSILON SQRT_WENO_EPSILON


//SMONO settings (if SMONO is enabled; these should not be here but it is nice to have them all in one file)
// The three below settings only control indicator, not value set
#define DO_SMONO_C2A 0
#define DO_SMONO_A2C 0
#define DO_SMONO_C2E 1

#define DO_MONO_1ST_DERIVATIVE 1
#define DO_3RD_DER 1

//#define DO_SMONO_C2E_CUSP_INDICATOR 1
#define DO_SMONO_C2E_CUSP_INDICATOR 0 // JCM GODMARK: causes problems for standard waves resolved with few points -- good for shocks if using AND reduction :: Also, doesn't do as good as expected for Noh type problems (both non-rel and rel) CHANGINGMARK: If 0, then causes problems for Noh but also BH-torus problem
#define DO_SMONO_A2C_CUSP_INDICATOR 0
#define DO_SMONO_C2A_CUSP_INDICATOR 0

////////
//
//  POST-/PRESHOCK REDUCTION
//
// reduce extra point downstream of shock
#define DO_REDUCE_POST_SHOCK 1
// reduce extra point upstream of shock
#define DO_REDUCE_PRE_SHOCK 1
// Fractional pressure jump for which the POST/PRESHOCK reduction fully takes place; it is linearly interpolated to smaller pressure jumps
#define REDUCE_DP_FRACTION (0.5)



////////
//
//  EXPERIMENTAL/DIDN'T WORK
//
#define DO_LIMIT_AC_CORRECTION_NEAR_INFLECTIONS 0
#define DO_MONOTONICITY_INDICATOR 0 // using new mono in interpline.c
#define DESENSITISE_STENCIL_REDUCTION 0  //whether to add epsilon in the calculation of the weights used for stencil reduction
#define WENO_HIGH_ORDER_CENTRAL_WEIGHT 0 // no code with this

///////////////////
// below 2 work quite well for smooth problems and doesn't hurt/help shocks
// only noticed slight noise issue in internal energy for caustic problem when doing full SWENO FV method
#define DO_LIMIT_C2E_WEIGHTS 1  //forces weights to be exactly equal if they are sufficiently close; forces those that are small but nonzero to be exactly zero
#define DO_LIMIT_AC_WEIGHTS 1  //forces weights to be exactly equal if they are sufficiently close; forces those that are small but nonzero to be exactly zero

// 20% variation allowed around equal and then made equal
// accounts for fact that weight calculation itself is inaccurate for (e.g.) 3rd order objects using normal WENO5 weights
// see limit_weight()
#define WEIGHTSPREADAROUNDEQUAL (0.2)

// for minimization, similar purpose to limit_ac_weights
#define WEIGHTFACTORMINIMIZE (0.3)
////////

//   END EXPERIMENTAL CODE
////////////////////////////////////////


#define FORCE_AC_CA_WEIGHT_TO_BE_OPTIMAL 0
#define FORCE_C2E_WEIGHT_TO_BE_OPTIMAL 0


// maximum order of scheme indicators can handle
#define MAXORDERS MAXSPACEORDER


#define MIN_CVT_ORDER 1
#define MAX_CVT_ORDER 7



#define MAX_STENCIL_LENGTH 5  //added by atch (to use as the max. size of an array)
#define MAX_NO_OF_WEIGHTS (3*MAX_STENCIL_LENGTH)  //we maximally store three sets of weights optimized for different points, that's why need to mult. by 3

#define LOWER_ORDER_VALUE (2) //currently chosen as the order of WENO-3, but could be different if use many orders for stencil reduction
#define MAX_NO_OF_LOWER_ORDER_WEIGHTS (3*LOWER_ORDER_VALUE)  //Use for array size to hold weights 
                                                             //instead of MAX_NO_OF_WEIGHTS to save on memory


struct weno_weights_s {
  // weight order in array memory is backwards from spatial order (e.g. WENO5 has 0,1,2 corresponding to right,center,left)
 
  ///////////////////////////////
  //array that holds the weights.  
  //
  //For WENO-5 center to edge there are 9 weights in this array:
  //(unopt_right,unopt_cent,unopt_left; optl_right,optl_cent,optl_left; optr_right,optr_cent,optr_left)
  //
  //For WENO-5 average to center/center to average there are 6 weights in this array:
  //(unopt_right,unopt_cent,unopt_left; opt_right,opt_cent,opt_left)
  ///////////////////////////////
  FTYPE weights[MAX_NO_OF_WEIGHTS];  

#if(1|| DO_STORE_SMOOTHNESS_INDICATORS ) // 1|| is to keep code not specialized when trying to access certain elements of structure
  FTYPE smoothness_indicators[MAX_CVT_ORDER];
#endif


  //total number of weights actually stored in the above array
  int len;                            

  //reconstruction order for which the weights are computed (=3 for WENO-5, =3 for WENO-3)
  //this equals the number of smoothness indicators stored in the smoothness_indicators array
  int order;             

  //the sum of unnormalized unoptimized weights
  FTYPE unoptimized_weights_sum;

  ////////////////////////////////
  //array that holds lower order weights
  //For WENO-3 center to edge there are 6 weights in this array:
  //(unopt_right,unopt,left; optl_right,optl_left; optr_right,optr_left)
  //
  //For WENO-5 average to center/center to average there are 4 trivial (i.e. each is equal to 0.5) 
  //weights in this array:
  //(unopt_right,unopt,left; opt_right,opt_left)
  ///////////////////////////////
  FTYPE lower_order_weights[MAX_NO_OF_LOWER_ORDER_WEIGHTS];

#if( WENO_C2E_REDUCE_NEAR_CUSPS || WENO_AC_REDUCE_NEAR_CUSPS  )
  int do_reduce;
#endif

#if( WENO_HIGH_ORDER_CENTRAL_WEIGHT )
  FTYPE high_order_central_weight;
#endif

  FTYPE w_ratio_min, w_ratio_max, w_ratio_left, w_ratio_right, w_ratio_combined, lower_order_fraction, lower_order_fraction_fordpp;
};





typedef struct weno_weights_s weno_weights_t;





#define LOWER_ORDER_PREFACTOR (1.) //set to zero to switch off stencil reduction
//fewer the number of times when the stencil size is reduced
#define WENO_MIN_NONZERO_WEIGHT (1.) //minimum weight that is considered non-zero

//The interval of weights where a linear combination of reconstructions with different orders are used
//Use different limiting weught ratios for OR and AND versions of stencil reduction.  Because the OR version is more
//sensitive to these ratios, they should be larger for it
#if( DO_OR_STENCIL_REDUCTION )
#define MAX_TRANSITION_RATIO   (15.L)
#define MIN_TRANSITION_RATIO    (9.L)
#define MAX_TRANSITION_RATIO_AC (15.L) //(1.6)  using such small weight thresholds does not work for 2d Noh problem -- lots of negative internal energies in smooth preshock region
#define MIN_TRANSITION_RATIO_AC (9.L)  //(1.3)
#else
#define MAX_TRANSITION_RATIO (1.6L)
#define MIN_TRANSITION_RATIO (1.3L)
#define MAX_TRANSITION_RATIO_AC (15.L)  //using such small weight thresholds does not work for 2d Noh problem -- lots of negative internal energies in smooth preshock region
#define MIN_TRANSITION_RATIO_AC (9.L)
#endif

#define WENO_DELTA_I (1)

#define MONOYIN MONOINDTYPE
#define LEFTYOUT MONOLEFTSET
#define RIGHTYOUT MONORIGHTSET
#define CENTYOUT MONOLEFTSET


//Implemented types of cell-centers to cel avg. reconstruction:
#define USEENO 0
#define USEWENO 1

#define WHICHENO USEWENO  //which type of cell-centers to cel avg. reconstruction to use
#define SQRT_WENO_EPSILON (500.0L*NUMEPSILON)  //used for determining the scale at which the weights change significantly
//#define SQRT_WENO_EPSILON (1.e-13L)  //used for determining the scale at which the weights change significantly
#define WENO_EPSILON (SQRT_WENO_EPSILON*SQRT_WENO_EPSILON)
//#define WENO_EPSILON (1.e-37L)
//#define WENO_EPSILON (1.e-26L)

//#define REDUCEEPSILON (NUMEPSILON * 100.L)
//#define REDUCEEPSILON (1E-10L)
#define REDUCEEPSILON WENO_EPSILON
#define DESENISTISE_STENCIL_REDUCTION_EPSILON (100.0L*NUMSQRTEPSILON)
//#define DESENISTISE_STENCIL_REDUCTION_EPSILON (1E-10L)



#define MATRIX_ROW_COL( matrix_name, index_row, index_col, num_columns ) matrix_name[(index_row) * (num_columns) + (index_col)]






#endif
