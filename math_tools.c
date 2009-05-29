#include "decs.h"


#define DO_ASSERTS 0 // similar parameter as in some init.h's

FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 ) 
//Arguments:
//x_eval -- abscissa of a point where an interpolated value is to be evaluated
//x1, x2, x3 -- set of abscissas; should all be different; can come in any order
//f1, f2, f3 -- set of function values at the above abscissas, i.e. f# = f( x# ), # = 1..3

//Description:
//Constructs a second-order polynomial interpolating the set of points {x#, f#}, #=1..3 and
//evaluates it at the point x.
//Interpolation is performed using the standard Lagrange method.
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






// round to a certain precision in base 10 to avoid round off errors
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


