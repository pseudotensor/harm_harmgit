
/*! \file math.tools.c
     \brief functions that don't depend upon any global things
*/


#include "decs.h"


#define DO_ASSERTS 0 // similar parameter as in some init.h's

///Description:
///Constructs a second-order polynomial interpolating the set of points {x#, f#}, #=1..3 and
///evaluates it at the point x.
///Interpolation is performed using the standard Lagrange method.
///Arguments:
///x_eval -- abscissa of a point where an interpolated value is to be evaluated
///x1, x2, x3 -- set of abscissas; should all be different; can come in any order
///f1, f2, f3 -- set of function values at the above abscissas, i.e. f# = f( x# ), # = 1..3
FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 ) 
{
  int i, j; //iterator through the abscissa points
#define max_interp_order  6 //should be equal to the number of points input to the function; degree of interpolation is only limited by the number of arguments to the function --
  //so increase that number if require a higher order; also, modify the definitions of x[] and f[] to include the added arguments

  FTYPE f_eval = 0.0; //interpolated value to be returned
  FTYPE x[] = { x1, x2, x3, x4, x5, x6};
  FTYPE f[] = { f1, f2, f3, f4, f5, f6};
  FTYPE x_a, x_b; //temporary variables
  FTYPE delta_f_eval;

  if( order < 1 || order > max_interp_order) {
    trifprintf( "interpn: requested order of interpolation non-positive or too high.\n" );
    myexit( 1 );
  }

  //test if the supplied x-values are valid: there should be no identical values
#if(DO_ASSERTS)
  for( i = 1; i < order; i++ ) {
    for( j = 0; j < i; j++ ) {
      if( x[i] - x[j] == 0.0 ) {
        fprintf( stderr, "interpn: abscissas of values to be interpolated are to be all different\n" );
        return 0.;
      }
    }
  }
#endif

  for( i = 0; i < order; i++ ) {
    //construct a polynomial that would be equal to f[i] at x[i] and equal to 0 at other x-points ( x[(i + j) % order], j = 1..order )
    x_a = x[i];  
    delta_f_eval = f[i];

    for( j = 1; j < order; j++ ) {
      x_b = x[(i + j) % order];
      delta_f_eval *= (x_eval - x_b) / (x_a - x_b);
    }

    //add up the value of this polynomial at x = x_eval to the final result
    f_eval += delta_f_eval;

  }

  return f_eval;
}






/// round to a certain precision in base 10 to avoid round off errors
FTYPE roundprecision(FTYPE value, int precision)
{
  FTYPE fraction,
    temp;
 
 
  value = value*pow(10, precision);
 
  fraction = modf(value, &temp);
  if(fraction>=0.5) temp+=1.0;
  if(fraction<=-0.5) temp-=1.0;
 
  return temp*pow(0.1, precision);
}


/// to use generically (e.g. for parabolic interpolation), call like:
/// interpfun(QUADRATICTYPE,3,1,realposition,array from 0 of positions, array from 0 of values, output of answer);
void interpfun(int interptype, int numpoints, int i, FTYPE pos, FTYPE *xfun, FTYPE *fun, FTYPE *answer)
{
  FTYPE slope,intercept;
  FTYPE slope1,slope2,xminusx0;
  FTYPE f0,f1,f2,x0,x1,x2;
  FTYPE linslope1,linslope2;


  ///////////////////////////////////////
  //
  // First setup points to use.  Restrict if near edge of data
  //
  ///////////////////////////////////////
  if(interptype==LINEARTYPE || interptype==LOGTYPE){
    if(i-1<0){
      f0=fun[i];
      f1=fun[i+1];
      x0=xfun[i];
      x1=xfun[i+1];
    }
    else if(i+1>=numpoints){
      f0=fun[i-1];
      f1=fun[i];
      x0=xfun[i-1];
      x1=xfun[i]; 
    }
    else{
      f0=fun[i-1];
      f1=fun[i];
      x0=xfun[i-1];
      x1=xfun[i];
    }
  }
  else if(interptype==QUADRATICTYPE){

    // quadratically interpolate fun using xfun
    if(i-1<0){
      f0=fun[i];
      f1=fun[i+1];
      f2=fun[i+2];
      x0=xfun[i];
      x1=xfun[i+1];
      x2=xfun[i+2]; 
    }
    else if(i+1>=numpoints){
      f0=fun[i-2];
      f1=fun[i-1];
      f2=fun[i];
      x0=xfun[i-2];
      x1=xfun[i-1];
      x2=xfun[i]; 
    }
    else{
      f0=fun[i-1];
      f1=fun[i];
      f2=fun[i+1];
      x0=xfun[i-1];
      x1=xfun[i];
      x2=xfun[i+1];
    }
  }


  ///////////////////////////////////////
  //
  // Interpolate
  //
  ///////////////////////////////////////


  ///////////////////////////////////////
  // LINEAR
  ///////////////////////////////////////
  if(interptype==LINEARTYPE){
    // linearly interpolate fun using pos
    slope = (f1-f0)/(x1-x0);
    intercept = f0;
    *answer = slope*(pos-x0) + intercept;
  }
  ///////////////////////////////////////
  // Quadratic with limiters
  ///////////////////////////////////////
  else if(interptype==QUADRATICTYPE){

    slope2 = ((f0-f2)/(x0-x2) - (f2-f1)/(x2-x1))/(x0-x1);
    slope1 = (f0-f1)/(x0-x1) + (f0-f2)/(x0-x2) - (f2-f1)/(x2-x1);

    linslope1=(f1-f0)/(x1-x0);
    linslope2=(f2-f1)/(x2-x1);
    
    // MINM truncation:
    if(linslope1*linslope2<0.0){
#if(0)
      if(pos<=x0){
        *answer=f0;
      }
      else if(pos>=x0 && pos<=x1){
        if(fabs(pos-x0)<fabs(pos-x1)) *answer=f0;
        else *answer=f1;
      }
      else if(pos>=x1 && pos<=x2){
        if(fabs(pos-x1)<fabs(pos-x2)) *answer=f1;
        else *answer=f2;
      }
      else *answer=f2;
#elif(1)
      if(fabs(linslope1)<fabs(linslope2)){
        slope = (f1-f0)/(x1-x0);
        intercept = f0;
        *answer = slope*(pos-x0) + intercept;
      }
      else{
        slope = (f2-f1)/(x2-x1);
        intercept = f1;
        *answer = slope*(pos-x1) + intercept;
      }
#elif(0)
      *answer=f1;
#endif
    }
    else{
      xminusx0 = (pos-x0);
      *answer = slope2*pow(xminusx0,2.0) + slope1*xminusx0 + f0;
    }

#if(1)
    // limit to original values' ranges
    if(*answer>f0 && *answer>f1 && *answer>f2){
      if(f0>f1 && f0>f2){
        *answer=f0;
      }
      else if(f1>f0 && f1>f2){
        *answer=f1;
      }
      else if(f2>f0 && f2>f1){
        *answer=f2;
      }
    }
    else if(*answer<f0 && *answer<f1 && *answer<f2){
      if(f0<f1 && f0<f2){
        *answer=f0;
      }
      else if(f1<f0 && f1<f2){
        *answer=f1;
      }
      else if(f2<f0 && f2<f1){
        *answer=f2;
      }
    }
#endif

  }
  ///////////////////////////////////////
  // Linear in Log
  ///////////////////////////////////////
  else if(interptype==LOGTYPE){
    // log interpolate fun using xfun
    slope = log(f1/f0)/log(x1/f0);
    if(fabs(slope)<1E-10) *answer=f1;
    else if(f0<0.0){
      // assume bi-log
      *answer=-exp( slope*log(pos/x0)+log(-f0) );
    }
    else *answer=exp( slope*log(pos/x0)+log(f0) );

    //dualfprintf(fail_file,"ii=%d jj=%d slope=%g myXfun=%g\n",ii,jj,slope,myXfun);
  }



}




