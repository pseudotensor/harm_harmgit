//#include "global.h"
#include "decs.h"
#include "reconstr.h"


int assert( int val, char *s )
{
    if( 0 != val ) {
	dualfprintf( fail_file, "Assertion failed: %s\n", s );
	myexit( 1 );
    }
    
    return 0;
}

FTYPE newtdiff( int stencil_dir, int order, int i_start, int j_start, int k, FTYPE uin[][N2M][NPR] )
//computes the Newton-divided difference of order 'int order' of k-th quantity 
//returns: order-th degree newton difference with 
//the point (i_start, j_start) being the leftmost in the stencil used to compute the difference; uses recursion
{
    FTYPE res;  //for storing the result
    int di, dj; //defines the unit vector in the direction of the desired Newton 
                //difference given by stencil_dir
    
    assert( stencil_dir != 1 && stencil_dir != 2, "newtdiff: stencil_dir can be 1 or 2" );
    assert( order < 0, "newtdiff: order cannot be negative" );
        
    if( 0 >= order ) {
	//the value itself
	return uin[i_start][j_start][k];
    }

    //order is non-zero, so express the difference by subtracting two lower-order differences (recursion)
    
    //set the unit vector in the dir. of the stencil
    if( 1 == stencil_dir ) {
	//stencil is along the x-direction
	di = 1;
	dj = 0;
    }
    else {
	//stencil is along the y-direction
	di = 0;
	dj = 1;
    }

    //Newton diff. of degree (order) is the difference of the lower-degree Newton differences of degree (order-1) 
    //in the direction stencil_dir
    res = newtdiff( stencil_dir, order - 1, i_start + di, j_start + dj, k, uin ) 
        - newtdiff( stencil_dir, order - 1, i_start,      j_start,      k, uin );
	        	
    return( res );

}

int choose_stencil( int stencil_dir, int *order, int min_index, int max_index, int i0, int j0, int k0, FTYPE uin[][N2M][NPR] )
//This routine assumes that the flow is smooth enough that it is always possible to choose a smooth stencil of length (*order)!!!,
//i.e. that discontinuities are not too close to each other

//(*order) - is the number of points in the stencil; its value maybe reduced in order to satisfy shift < maxshift
//min_index & max_index - max and min values of the index where stencil can be (in the direction of stencil_dir)

//returns:
// the shift of the stencil in the left direction;
//updates: the value of (*order) only in the case if (max_index - min_index + 1 < *order).
{
  //see the Shu review (1997), the code given after (2.46); note that here we operate on v's (derivative of V's), not V's
  //as is done in Shu (1997).

  FTYPE leftdiff; 
  FTYPE rightdiff;

  FTYPE bias_factor = (FTYPE)2.0L;
  FTYPE no_bias_factor = (FTYPE)1.0L;

  FTYPE leftfactor = no_bias_factor;
  FTYPE rightfactor = no_bias_factor;
  FTYPE free_shift_bias_factor = no_bias_factor;

  int can_expand_left = 1;
  int can_expand_right = 1;


#define DIR_INDEX(i,j) ((i) * ei + (j) * ej)    //index of (i,j) point in the direction of stencil_dir

  int is = i0, js = j0; //the left-most point in the stencil for the cell i0 (or j0) for an order-th degree polynomial
  int index_pref; //preferred value of the index of the leftmost point in the stencil (centered stencil)
    
  int m; //iterator
  int shift; //stencil left-shift
    
  int ei, ej; //defines the unit vector in the direction of the desired Newton 
  //difference given by stencil_dir
    
  assert( stencil_dir != 1 && stencil_dir != 2, "choose_stencil: stencil_dir can be 1 or 2" );
  assert( order < 0, "choose_stencil: order cannot be negative" );
        
  //set the unit vector in the dir. of the stencil
  if( 1 == stencil_dir ) {
    //stencil is along the x-direction
    ei = 1;
    ej = 0;
  }
  else {
    //stencil is along the y-direction
    ei = 0;
    ej = 1;
  }
    
  assert( DIR_INDEX(i0, j0) < min_index || DIR_INDEX(i0, j0) > max_index, 
	  "choose_stencil: [min_index, max_index] does not contain initial point" );
    
  //not enough zones for interpolation of degree == order
  if( max_index - min_index + 1 < *order ) {
    (*order) = max_index - min_index + 1; //decrease the interpolation degree
  }

  index_pref = DIR_INDEX(i0,j0) - (*order)/2;

  for( m = 1, shift = 0; m < *order; m++ ) {
    if( DIR_INDEX( is, js ) == min_index ) {
      //can't expand the stencil to the left any more
      can_expand_left = 0;
      continue;      
    }

    if( DIR_INDEX( is, js ) + m - 1 == max_index ) { 
      //can't expand the stencil to the right any more
      can_expand_right = 0;
      is -= 1 * ei;
      js -= 1 * ej;
      shift++;
      continue;
    }

    //calculate the left and right newton differences -- we use their abs. values as indicators of
    //smoothness in the left and right directions respectively
    leftdiff = newtdiff( stencil_dir, m, is - 1 * ei, js - 1 * ej, k0, uin );
    rightdiff = newtdiff( stencil_dir, m, is, js, k0, uin );
	
    //alleviate the ''free adaptation of stencils'' problem, see Shu report (1997), sec. 2.2.2., ~(2.49)
    if( DIR_INDEX(is,js) <= index_pref ) {
      //make expansion to the left less likely since further expansion to the left would lead to an un-preferred stencil
      free_shift_bias_factor = 1 / bias_factor;
    }
    else {
      free_shift_bias_factor = bias_factor;
    }


    if( fabs(leftdiff) < free_shift_bias_factor * fabs(rightdiff) ) {
      //the function is smoother to the left than to the right
      //try to expand the stencil to the left if we can

      is -= 1 * ei;
      js -= 1 * ej;
      shift++;
    }
  }

  ////if everything is successfull, m is equal to (*order), so the following does nothing;
  ////otherwise it limits the interpolation order to the last successful one, m
  //(*order) = m;
    
  return shift;

}

int cvtij( int stencil_dir, int cvt_type, int i0, int j0, FTYPE uin[][N2M][NPR], FTYPE *uout )
{
    int i, j, k, d, shift, rmin, rmax;  // r always >= 0
    int dir_index, max_index, ei, ej;
    int default_order = 3, order;
    FTYPE *cvtmatrix;

    //input checks
    assert( stencil_dir < 1 || stencil_dir > 2, "cvtij: illegal stencil dir" );
    assert( cvt_type != CVT_A2C && cvt_type != CVT_C2A, "cvtij: no such conversion type" );
    assert( i0 < 0 || j0 < 0 || i0 > N1 -  1 || j0 > N2 - 1, "cvtij: (i0, j0) out of bounds" );
    assert( NULL == uout, "cvtij: FTYPE *uout cannot be a NULL pointer" );
    assert( default_order < 0 || default_order > MAX_CVT_ORDER, "cvtij: order of interpolation is out of bounds" );
    
    if( 1 == stencil_dir ) {
        ei = 1;
        ej = 0;
        dir_index = i0;
        max_index = N1 - 1;
    }
    else {
        ei = 0;
        ej = 1;
        dir_index = j0;
        max_index = N2 - 1;
    }


    PLOOP {
      //interpolate each conserved variable separately
      order = default_order; //start out with some initial order
      shift = choose_stencil( stencil_dir, &(order), 0, max_index, i0, j0, k, uin );  //this can change order!!

      if( order < MIN_CVT_ORDER ) {
	//not enough grid cells available -- cannot do anything, so leave values as they are
	uout[k] = uin[i0][j0][k];
	continue;
      }

      do {
	//retrieve the conversion matrix for the chosen transformation type and order
	//and keep lowering the order until the non-NULL matrix found
	cvtmatrix = cvt_matrices[cvt_type][order];
      }
      while( 
	    NULL == cvtmatrix  //if no matrix exists for this order
	 && --order > 0        //decrease the order and
	 && ((shift > order/2)?(shift--):(1))  //adjust the shift to keep stencil close to centered one
	    );

      uout[k] = 0.0;

      //initialize iterators
      i = i0 - shift * ei;
      j = j0 - shift * ej;

      for( d = 0; d < order; d++, i += ei, j += ej )
	{
	  uout[k] += cvtmatrix[shift * order + d] * uin[i][j][k]; //see (2.11) from Shu Report (1997)
	}
    }
    
    return( 0 );
    
}

FTYPE a2cij( int dir, int i, int j, FTYPE ua[][N2M][NPR], FTYPE *uc )
// converts cell averaged to cell centered quantities; output is written to the last argument of the function
{
    cvtij( dir, CVT_A2C, i, j, ua, uc );
    
    return( 0.0 );
}

FTYPE c2aij( int dir, int i, int j, FTYPE uc[][N2M][NPR], FTYPE *ua )
//converts cell-centered quantities to cell averaged; output is written to the last argument of the function
{
    cvtij( dir, CVT_C2A, i, j, uc, ua );
    
    return( 0.0 );
}


FTYPE a2c( FTYPE uc[][N2M][NPR], FTYPE ua[][N2M][NPR] )
// converts cell averaged to cell centered quantities; output is written to the first argument of the function
{
    int i, j;
    
    CZLOOP {
	a2cij( 1, i, j, ua, uc[i][j] );
    }
    return( 0.0 );
}
	
FTYPE c2a( FTYPE ua[][N2M][NPR], FTYPE uc[][N2M][NPR] )
//converts self-centered quantities to cell averaged; output is written to the first argument of the function
{
    int i, j;
    CZLOOP {
	c2aij( 1, i, j, uc, ua[i][j] );
    }
	    
    return(0.);
}

	
	
