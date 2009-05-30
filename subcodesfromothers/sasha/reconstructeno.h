#define MIN_CVT_ORDER 3
#define MAX_CVT_ORDER 5
#define CVT_A2C 0
#define CVT_C2A 1
#define NUM_CVT_TYPES 2

//FTYPE a2c( FTYPE uc[][N2M][NPR], FTYPE ua[][N2M][NPR] );
//FTYPE c2a( FTYPE ut[][N2M][NPR], FTYPE uc[][N2M][NPR] );

int cvtij( int stensil_dir, int cvt_type, int i0, int j0, FTYPE uin[][N2M][NPR], FTYPE *uout );
int choose_stencil( int stencil_dir, int *order, int min_index, int max_index, int i0, int j0, int k0, FTYPE uin[][N2M][NPR] );
FTYPE a2cij( int dir, int i, int j, FTYPE ua[][N2M][NPR], FTYPE *uc );
FTYPE c2aij( int dir, int i, int j, FTYPE ua[][N2M][NPR], FTYPE *uc );
int assert( int val, char *s );
FTYPE newtdiff( int stencil_dir, int order, int i0, int j0, int k0, FTYPE uin[][N2M][NPR] );


FTYPE a2c_array3[3][3] =
   {{0.9583333333333334,0.08333333333333333,
    -0.041666666666666664},
   {-0.041666666666666664,1.0833333333333333,
    -0.041666666666666664},
   {-0.041666666666666664,0.08333333333333333,
    0.9583333333333334}};

FTYPE c2a_array3[3][3] =
    {{1.0416666666666667,-0.08333333333333333,
    0.041666666666666664},{0.041666666666666664,
    0.9166666666666666,0.041666666666666664},
   {0.041666666666666664,-0.08333333333333333,
    1.0416666666666667}};

FTYPE a2c_array5[5][5] = {
  {563./640., 57./160., -373./960., 91./480., -71./1920.}, 
  {-71./1920., 511./480., -13./960., -3./160., 3./640.}, 
  {3./640., -29./480., 1067./960., -29./480.,  3./640.}, 
  {3./640., -3./160., -13./960., 511./480., -71./1920.}, 
  {-71./1920., 91./480., -373./960., 57./160., 563./640}};

FTYPE c2a_array5[5][5] = {
  {6463./5760., -523./1440., 383./960., -283./1440., 223./5760.}, 
  {223./5760., 1337./1440., 23./960., 17./1440., -17./5760.}, 
  {-17./5760., 77./1440., 863./960., 77./1440., -17./5760.}, 
  {-17./5760., 17./1440., 23./960., 1337./1440., 223./5760.}, 
  {223./5760., -283./1440., 383./960., -523./1440., 6463./5760}};

FTYPE *cvt_matrices[2][MAX_CVT_ORDER + 1] = {
  {
    //CVT_A2C
    NULL,
    NULL,
    NULL,
    (FTYPE *)a2c_array3, //order 3
    NULL,
    (FTYPE *)a2c_array5  //order 5
  },
  {
    //CVT_C2A
    NULL,
    NULL,
    NULL,
    (FTYPE *)c2a_array3, //order 3
    NULL,
    (FTYPE *)c2a_array5  //order 5
  }
};

