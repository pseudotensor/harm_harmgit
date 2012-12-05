/*************************************************************************

                        Mathematica source file

        Copyright 1986 through 1999 by Wolfram Research Inc.


*************************************************************************/

/* C language definitions for use with Mathematica output */


#define Power(x, y)	(powl((long double)(x), (long double)(y)))
#define Sqrt(x)		(sqrtl((long double)(x)))
#define Sqrtl(x)        (sqrtl((long double)(x)))

#define Abs(x)		(fabsl((long double)(x)))

#define Exp(x)		(expl((long double)(x)))
#define Log(x)		(logl((long double)(x)))

#define Sin(x)		(sinl((long double)(x)))
#define Cos(x)		(cosl((long double)(x)))
#define Tan(x)		(tanl((long double)(x)))

#define ArcSin(x)       (asinl((long double)(x)))
#define ArcCos(x)       (acosl((long double)(x)))
#define ArcTan(x)       (atanl((long double)(x)))

#define Sinh(x)          (sinhl((long double)(x)))
#define Cosh(x)          (coshl((long double)(x)))
#define Tanh(x)          (tanhl((long double)(x)))

#define Cot(x)          (1./tanl((long double)(x)))
#define Csc(x)          (1./sinl((long double)(x)))




/** Could add definitions for Random(), SeedRandom(), etc. **/


