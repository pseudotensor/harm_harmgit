#ifndef RECONSTRUCTENO_STATICFUNCS_CONSTVARS_H
#define RECONSTRUCTENO_STATICFUNCS_CONSTVARS_H





static FTYPE normalize_array( int number_of_elements, FTYPE epsilon, FTYPE *array_to_normalize );
static FTYPE get_sum_of_elements( int number_of_elements, FTYPE *array_to_sum_up );
static void compute_weights_from_smoothness_indicators( int pl, int order, FTYPE epsilon, FTYPE *uin, FTYPE *smoothness_indicators, FTYPE *stencil_weights );
static FTYPE get_weno5_high_order_central_weight( FTYPE *u );
static FTYPE apply_weno5_high_order_central_weight( int order, FTYPE *u, FTYPE *stencil_weights, FTYPE *stencil_weights_leftface, FTYPE *stencil_weights_rightface );
static int compute_weights_ratios(  int cvt_type, int bs, int bf, FTYPE *uin, weno_weights_t *stencil_weights_out );
static FTYPE squareit( FTYPE base );

static void compute_monotonicity_indicator(int order, int minindex, int maxindex, weno_weights_t *optimal_weights, FTYPE *pr, int *monoindicator);
static FTYPE *get_cvt_vec( int cvt_type, int full_order );
static void ac_simple_eno_odd(int cvt_type, int order, FTYPE *yin, FTYPE *pout);
static FTYPE rescale_weight( FTYPE order, FTYPE weight ); 
static void rescale_unoptimized_weights( weno_weights_t *weights );
static FTYPE compute_frac_equal_weight( FTYPE order, FTYPE weight );

static void create_weights_from_array( int cvt_type, int order, FTYPE *weights_array_in, weno_weights_t *stencil_weights_out );
static void extract_weights_into_array( weno_weights_t *stencil_weights_in, FTYPE *weights_array_out );

static int check_symmetry_of_weno_input( int cvt_type, int pl, int bs, int bf, FTYPE *yin );
static int check_symmetry_of_weno_output( int cvt_type, FTYPE (*monoindicator)[NBIGM], int pl, int ps, int pf, int maxorder, FTYPE *yout[MAX_NO_OF_WENO_OUTPUTS], weno_weights_t *stencil_weights_array );
static int do_reduce_near_cusps( int min_index, int max_index, FTYPE *uin );

static void limit_weight( int order, FTYPE old_weight, FTYPE *new_weight );
static void mono_adjust_unoptimized_weights( FTYPE *monoindicator, weno_weights_t *stencil_weights_array );

static int eno_cvt( int cvt_type, int min_index, int max_index, weno_weights_t *stencil_weights_struct, FTYPE *uin, FTYPE *uout );
static int compute_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out );
static int compute_lower_order_stencil_weights(  int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out );
static int compute_ac_ca_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out );
static int compute_cf_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator0, FTYPE *monoindicator1, FTYPE *uin, weno_weights_t *stencil_weights_out );

static int compute_optimized_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out );
static int compute_optimized_ac_ca_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out );
static int compute_optimized_cf_stencil_weights( int cvt_type, int pl, int order, int min_index, int max_index, FTYPE *monoindicator, FTYPE *uin, weno_weights_t *stencil_weights_out );

static void apply_additional_reduction_to_weights( int cvt_type, int whichreduce, int max_order, int min_order, int i0, int pl, int bs, int bf, FTYPE (*shockindicator)[NBIGM], FTYPE (*dP)[NBIGM],
                                                   FTYPE *monoindicator, FTYPE *Pindicator, FTYPE *uin, 
                                                   weno_weights_t *stencil_weights_array, struct of_trueijkp *trueijkp );
static int choose_weno_order( int cvt_type, int whichreduce, int max_order, int min_order, int i0, int pl, int bs, int bf, FTYPE (*shockindicator)[NBIGM], FTYPE (*dP)[NBIGM],
                              FTYPE *monoindicator0, FTYPE *monoindicator1, 
                              FTYPE *Pindicator, FTYPE *uin, 
                              weno_weights_t *stencil_weights_array, 
                              weno_weights_t **pp_stencil_weights_to_be_used, struct of_trueijkp *trueijkp  );

static void desensitise_smoothness_indicators( int order, FTYPE epsilon, FTYPE *uin, FTYPE *smoothness_indicators);


static int eno_line_reconstruct( int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int cvt_type, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, 
                                 int *minorderit, int *maxorderit, int *shiftit, 
                                 FTYPE (*shockindicator)[NBIGM], 
                                 FTYPE (*df)[NBIGM],
                                 FTYPE (*dP)[NBIGM],
                                 FTYPE (*monoindicator)[NBIGM],
                                 FTYPE *Pindicator, 
                                 FTYPE *yin,  FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);



static FTYPE compute_lower_order_fraction( FTYPE w_ratio_left, FTYPE w_ratio_right );


// definition
typedef struct weno_output_s  {
  int len;
  int type[MAX_NO_OF_WENO_OUTPUTS];
} weno_outputs_t;


// Here we define what is set in reconstructeno_set_vars.c:
//static int weno_weights_shifts_array[];
//static int (*compute_stencil_weight_func_ptrs[NUM_CVT_TYPES])( int, int, int, int, int, FTYPE *, FTYPE *, FTYPE*, weno_weights_t * );
//static int (*compute_optimized_stencil_weight_func_ptrs[])( int, int, int, int, int, FTYPE *, FTYPE *, weno_weights_t * );
//static weno_outputs_t weno_outputs[];
//static FTYPE dpos_array[][2][3];
//static FTYPE dneg_array[][2][3]
//static FTYPE c2e_optimal_weights_leftface[][5];
//static FTYPE c2e_optimal_weights_rightface[][5];
//static FTYPE a2c_array2[][2];
//static FTYPE c2a_array2[][2];
//static FTYPE a2c_array1[][1];
//static FTYPE c2a_array1[][1];



////////////////////////////////////////
//
// Static [JCM: actually, not written to ever during computation] memory spaces
//
////////////////////////////////////////




static const int weno_weights_shifts_array[] = 
  {
    1,   //CVT_A2C
    1,  //CVT_C2A
    1,  //CVT_C2L
    2,  //CVT_C2R
    1,   //CVT_C2DER1
    1,   //CVT_C2DER2
    1,   //CVT_C2DER3
    1    //CVT_C2DER4
  };



static weno_outputs_t weno_outputs[8] = {
  {1, {CVT_A2C, -1}},   //CVT_A2C  //only centered value
  {1, {CVT_C2A, -1}},  //CVT_C2A  //only centered value
  {2, {CVT_C2L, CVT_C2R}},  //CVT_C2L  //two values: left and right
  {2, {CVT_C2L, CVT_C2R}},   //CVT_C2R  //two values: left and right
  //below for MERGEDA2CMETHOD
  {1, {CVT_C2DER1, -1}},   
  {1, {CVT_C2DER2, -1}},  
  {1, {CVT_C2DER3, -1}},  
  {1, {CVT_C2DER4, -1}}   
};
 

//Using the Shi, Hu, Shu (2002) "A Technique of Treating Negative Weights in WENO Schemes"
//####optimal weights; neeed to rename it to sth more sensible
//positive and negative optimal weights should not be normalized here; because their norm is important in getting the 
//real weights right.  The formula for positive weights is pos_weight = 1/2 ( weight + 3 abs(weight) );  neg_weight = 1/2 ( - weight + 3 abs(weight) )
//Then weight = pos_weight - neg_weight.  Negative weights are taken to be positive; then they always enter with a negative sign as in the preceding eqn.
static const FTYPE dpos_array[8][2][3] = //d[order-2][shift]  
  {
    {  //CVT_A2C
      {1./2., 1./2., 0.0},
      {9./80., 49./20., 9./80.}
    },  //end of CVT_A2C == 0
    {  //CVT_C2A
      {1./2., 1./2., 0.0},
      {17./240., 137./60., 17./240.}
    },  //end of CVT_C2A
    {  //CVT_C2L
      {0, 0, 0},
      {0, 0, 0},
    },  //end of CVT_C2L
    {      //CVT_C2R
      {0, 0, 0},
      {0, 0, 0},
    },  //end of CVT_C2R
    {
      {0.5, 0.5, 0.0},
      {1./3, 4./3, 1./3},
    },
    //below for MERGEDA2CMETHOD
    //CVT_C2DER2
    {
      {0.5, 0.5, 0.0},
      {1./3., 2, 1./6.}
    },
    //CVT_C2DER3
    {
      {0.5, 0.5, 0.0},
      {1, 4, 1}
    },
    //CVT_C2DER4
    {
      {0.5, 0.5, 0.0},
      {2, 2, 2},
    }
  };  

static const FTYPE dneg_array[8][2][3] = //d[order-2][shift]
  {
    { //CVT_A2C
      {1./2., 1./2., 0.0},
      {9./40., 49./40., 9./40.}
    },  //end of CVT_A2C == 0
    { //CVT_C2A
      {1./2., 1./2., 0.0},
      {17./120., 137./120., 17./120.}
    },   //end of CVT_C2A == 1
    {  //CVT_C2L
      {0, 0, 0},
      {0, 0, 0},
    },  //end of CVT_C2L
    {      //CVT_C2R
      {0, 0, 0},
      {0, 0, 0},
    },  //end of CVT_C2R
    //CVT_C2DER1
    {
      {0.5, 0.5, 0.0},
      {1./6, 2./3, 1./6}
    },
    //CVT_C2DER2
    {
      {0.5, 0.5, 0.0},
      {1./6, 1, 1./3}
    },
    //CVT_C2DER3
    {
      {0.5, 0.5, 0.0},
      {2, 2, 2}
    },
    //CVT_C2DER4
    {
      {0.5, 0.5, 0.0},
      {1, 4, 1}
    }
  };

  
static const FTYPE c2e_optimal_weights_leftface[4][5] = //c2e_optimal_weights_leftface[order-2][shift] -- left face optimal weights (a rectangular matrix where unused values are zeroes)
  {{1./4., 3./4, 0.0, 0.0, 0.0},
   { 1./16., 5./8., 5./16., 0.0, 0.0 },
   {0, 0, 0, 0, 0},
   {1./256., 9./64., 63./128., 21./64., 9./256.}
  };  

static const FTYPE c2e_optimal_weights_rightface[4][5] = //c2e_optimal_weights_rightface[order-2][shift] -- right face optimal weights (a rectangular matrix where unused values are zeroes)
  {{3./4., 1./4., 0.0, 0.0, 0.0},
   { 5./16., 5./8., 1./16., 0.0, 0.0 },
   {0, 0, 0, 0, 0},
   {9./256., 21./64., 63./128., 9./64., 1./256.}
  };  

static const FTYPE a2c_array2[2][2] =
  {{ 1., 0. }, { 0., 1. }};

static const FTYPE c2a_array2[2][2] =
  {{ 1., 0. }, { 0., 1. }};


static const FTYPE a2c_array1[1][1] =  //avg to center
  {{ 1.0 }};

static const FTYPE c2a_array1[1][1] = //center to avg
  {{ 1.0 }};

static const FTYPE a2c_array3[3][3] =
  {{0.9583333333333334,0.08333333333333333,
    -0.041666666666666664},
   {-0.041666666666666664,1.0833333333333333,
    -0.041666666666666664},
   {-0.041666666666666664,0.08333333333333333,
    0.9583333333333334}};

static const FTYPE c2a_array3[3][3] =
  {{1.0416666666666667,-0.08333333333333333,
    0.041666666666666664},{0.041666666666666664,
                           0.9166666666666666,0.041666666666666664},
   {0.041666666666666664,-0.08333333333333333,
    1.0416666666666667}};

static const FTYPE a2c_array5[5][5] = {
  {563./640., 57./160., -373./960., 91./480., -71./1920.}, 
  {-71./1920., 511./480., -13./960., -3./160., 3./640.}, 
  {3./640., -29./480., 1067./960., -29./480.,  3./640.}, 
  {3./640., -3./160., -13./960., 511./480., -71./1920.}, 
  {-71./1920., 91./480., -373./960., 57./160., 563./640}};

static const FTYPE c2a_array5[5][5] = {
  {6463./5760., -523./1440., 383./960., -283./1440., 223./5760.}, 
  {223./5760., 1337./1440., 23./960., 17./1440., -17./5760.}, 
  {-17./5760., 77./1440., 863./960., 77./1440., -17./5760.}, 
  {-17./5760., 17./1440., 23./960., 1337./1440., 223./5760.}, 
  {223./5760., -283./1440., 383./960., -523./1440., 6463./5760}};

static const FTYPE a2c_array7[7][7] = 
  {{88069./107520., 36961./53760., -122141./107520., 
    28991./26880., -22327./35840., 
    2173./10752., -3043./107520.}, {-3043./107520., 
                                    10937./10752., 10019./107520., -1303./8960., 
                                    3153./35840., -513./17920., 
                                    143./35840.}, {143./35840., -3023./53760., 
                                                   118379./107520., -1249./26880., -207./35840., 
                                                   15./3584., -5./7168.}, {-5./7168., 
                                                                           159./17920., -7621./107520., 30251./26880., -7621./107520., 
                                                                           159./17920., -5./7168.}, {-5./7168., 
                                                                                                     15./3584., -207./35840., -1249./26880., 
                                                                                                     118379./107520., -3023./53760., 
                                                                                                     143./35840.}, {143./35840., -513./17920., 
                                                                                                                    3153./35840., -1303./8960., 10019./107520., 
                                                                                                                    10937./10752., -3043./107520.}, {-3043./107520., 
                                                                                                                                                     2173./10752., -22327./35840., 
                                                                                                                                                     28991./26880., -122141./107520., 36961./53760., 88069./107520.}};

static const FTYPE c2a_array7[7][7] =
  {{1152511./967680., -7969./10752., 
    134881./107520., -294659./241920., 
    76921./107520., -12629./53760., 32119./967680.}, {32119./967680., 
                                                      154613./161280., -14237./322560., 
                                                      22441./241920., -18157./322560., 
                                                      593./32256., -2489./967680.}, {-2489./967680., 
                                                                                     8257./161280., 291803./322560., 11101./241920., 
                                                                                     883./322560., -367./161280., 
                                                                                     367./967680.}, {367./967680., -281./53760., 6361./107520., 
                                                                                                     215641./241920., 6361./107520., -281./53760., 
                                                                                                     367./967680.}, {367./967680., -367./161280., 883./322560., 
                                                                                                                     11101./241920., 291803./322560., 
                                                                                                                     8257./161280., -2489./967680.}, {-2489./967680., 
                                                                                                                                                      593./32256., -18157./322560., 
                                                                                                                                                      22441./241920., -14237./322560., 154613./161280., 
                                                                                                                                                      32119./967680.}, {32119./967680., -12629./53760., 
                                                                                                                                                                        76921./107520., -294659./241920., 
                                                                                                                                                                        134881./107520., -7969./10752., 1152511./967680.}};

static const FTYPE c2l_array4[4][4] = 
  {{35./16., -35./16., 21./16., -5./16.}, {
      5./16., 15./16., -5./16., 1./16.}, {-1./16.,
                                          9./16., 9./16., -1./16.}, {1./16., -5./16., 15./16., 5./16.}};


static const FTYPE c2r_array4[4][4] = 
  {{5./16., 15./16., -5./16., 1./16.}, {-1./16., 9./16., 
                                        9./16., -1./16.}, {1./16., -5./16., 15./16.,
                                                           5./16.}, {-5./16., 21./16., -35./16., 35./16.}};


static const FTYPE c2l_array5[5][5] = {{315./128., -105./32., 189./64., -45./32., 
                                        35./128.}, {35./128., 35./32., -35./64., 
                                                    7./32., -5./128.}, {-5./128., 15./32., 
                                                                        45./64., -5./32., 3./128.}, {3./128., -5./32., 45./64., 
                                                                                                     15./32., -5./128.}, {-5./128., 7./32., -35./64., 
                                                                                                                          35./32., 35./128.}};

static const FTYPE c2r_array5[5][5] = {{35./128., 35./32., -35./64., 
                                        7./32., -5./128.}, {-5./128., 15./32., 
                                                            45./64., -5./32., 3./128.}, {3./128., -5./32., 45./64., 
                                                                                         15./32., -5./128.}, {-5./128., 7./32., -35./64., 
                                                                                                              35./32., 35./128.}, {35./128., -45./32., 189./64., -105./32., 
                                                                                                                                   315./128.}};

static const FTYPE c2l_array3[3][3] =  //center to left face
  {{15./8., -5./4., 3./8.}, 
   {3./8., 3./4., -1./8.}, 
   {-1./8., 3./4., 3./8.}};

static const FTYPE c2r_array3[3][3] = //center to right face
  {{3./8., 3./4., -1./8.}, 
   {-1./8., 3./4., 3./8.}, 
   {3./8., -5./4., 15./8.}};

static const FTYPE c2l_array2[2][2] =  //center to left face
  {{3./2., -1./2.}, {1./2., 1./2.}};

static const FTYPE c2r_array2[2][2] = //center to right face
  {{1./2., 1./2.}, {-1./2., 3./2.}};

static const FTYPE c2l_array1[1][1] =  //center to left face
  {{ 1.0 }};

static const FTYPE c2r_array1[1][1] = //center to right face
  {{ 1.0 }};

//#if(MERGEDC2EA2CMETHOD)

// FIRST ORDER RECONSTRUCTION (WENO-1, k = 1)
static const FTYPE c2der1_array1[1][1] = 
  {{ 0.0 }};
static const FTYPE c2der2_array1[1][1] = 
  {{ 0.0 }};
static const FTYPE c2der3_array1[1][1] = 
  {{ 0.0 }};
static const FTYPE c2der4_array1[1][1] = 
  {{ 0.0 }};

// THIRD ORDER RECONSTRUCTION (WENO-3, k = 2)
static const FTYPE c2der1_array2[2][2] =  //1st der
  {{-1, 1}, {-1, 1}};
static const FTYPE c2der2_array2[2][2] =  //2nd der -- can in principle be obtained but wrong if weights deviate from optimal weights
  {{0., 0,}, {0., 0.}}; 
static const FTYPE c2der3_array2[2][2] =  //3rd der
  {{0., 0,}, {0., 0.}}; 
static const FTYPE c2der4_array2[2][2] =  //4th der
  {{0., 0,}, {0., 0.}}; 

// FIFTH ORDER RECONSTRUCTION (WENO-5, k = 3)
static const FTYPE c2der1_array3[3][3] =  //1st der
  {{-1.5,2,-0.5},{-0.5,0,0.5},{0.5,-2,1.5}};

static const FTYPE c2der2_array3[3][3] =  //2nd der
  {{-1.5,2,-0.5},{1,-2,1},{0.5,-2,1.5}};

static const FTYPE c2der3_array3[3][3] =  //3rd der
  {{-1.5,2,-0.5},{-0.5,0,0.5},{0.5,-2,1.5}};

static const FTYPE c2der4_array3[3][3] =  //4th der
  {{1,-2,1},{1,-2,1},{1,-2,1}};



///////////////////////////////////
//
// non-const static array of pointers
//
///////////////////////////////////

static int (*compute_stencil_weight_func_ptrs[NUM_CVT_TYPES])( int, int, int, int, int, FTYPE *, FTYPE *, FTYPE*, weno_weights_t * ) = {
  compute_ac_ca_stencil_weights,
  compute_ac_ca_stencil_weights,
  compute_cf_stencil_weights,
  compute_cf_stencil_weights
};


static int (*compute_optimized_stencil_weight_func_ptrs[])( int, int, int, int, int, FTYPE *, FTYPE *, weno_weights_t * ) = {
  compute_optimized_ac_ca_stencil_weights,
  compute_optimized_ac_ca_stencil_weights,
  compute_optimized_cf_stencil_weights,
  compute_optimized_cf_stencil_weights,
  //below for MERGEDA2CMETHOD:
  compute_optimized_ac_ca_stencil_weights,
  compute_optimized_ac_ca_stencil_weights,
  compute_optimized_ac_ca_stencil_weights,
  compute_optimized_ac_ca_stencil_weights
};


static FTYPE *cvt_matrices[8][MAX_CVT_ORDER + 1] = {
  {
    //CVT_A2C
    NULL,
    (FTYPE *)a2c_array1,
    (FTYPE *)a2c_array2, //order 2
    (FTYPE *)a2c_array3, //order 3
    NULL,
    (FTYPE *)a2c_array5,  //order 5
    NULL,
    (FTYPE *)a2c_array7  //order 7
  },
  {
    //CVT_C2A
    NULL,
    (FTYPE *)c2a_array1,
    (FTYPE *)c2a_array2, //order 2
    (FTYPE *)c2a_array3, //order 3
    NULL,
    (FTYPE *)c2a_array5,  //order 5
    NULL,
    (FTYPE *)c2a_array7  //order 7
  },
  {
    //CVT_C2L
    NULL,
    (FTYPE *)c2l_array1,
    (FTYPE *)c2l_array2,
    (FTYPE *)c2l_array3, //order 3
    (FTYPE *)c2l_array4, //order 4
    (FTYPE *)c2l_array5, //order 5
    NULL,
    NULL  //order 7
  },
  {
    //CVT_C2R
    NULL,
    (FTYPE *)c2r_array1,
    (FTYPE *)c2r_array2,
    (FTYPE *)c2r_array3, //order 3
    (FTYPE *)c2r_array4, //order 4
    (FTYPE *)c2r_array5, //order 5
    NULL,
    NULL  //order 7
  },
  //below for merged A2C method
  {
    //CVT_C2DER1
    NULL,
    (FTYPE *)c2der1_array1,
    (FTYPE *)c2der1_array2,
    (FTYPE *)c2der1_array3, //order 3
    NULL, //order 4
    NULL, //order 5
    NULL,
    NULL  //order 7
  },
  {
    //CVT_C2DER2
    NULL,
    (FTYPE *)c2der2_array1,
    (FTYPE *)c2der2_array2,
    (FTYPE *)c2der2_array3, //order 3
    NULL, //order 4
    NULL, //order 5
    NULL,
    NULL  //order 7
  },
  {
    //CVT_C2DER3
    NULL,
    (FTYPE *)c2der3_array1,
    (FTYPE *)c2der3_array2,
    (FTYPE *)c2der3_array3, //order 3
    NULL, //order 4
    NULL, //order 5
    NULL,
    NULL  //order 7
  },
  {
    //CVT_C2DER4
    NULL,
    (FTYPE *)c2der4_array1,
    (FTYPE *)c2der4_array2,
    (FTYPE *)c2der4_array3, //order 3
    NULL, //order 4
    NULL, //order 5
    NULL,
    NULL  //order 7
  },
};





#endif



