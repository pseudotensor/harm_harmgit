/*************************************************************************

                        Mathematica source file

        Copyright 1986 through 1999 by Wolfram Research Inc.


*************************************************************************/

/* C language definitions for use with Mathematica output */


#define Power(x, y)	(pow((x), (y)))
#define Sqrt(x)		(sqrt((x)))

#define Abs(x)		(fabs((x)))

#define Exp(x)		(exp((x)))
#define Log(x)		(log((x)))

#define Sin(x)		(sin((x)))
#define Cos(x)		(cos((x)))
#define Tan(x)		(tan((x)))

#define ArcSin(x)       (asin((x)))
#define ArcCos(x)       (acos((x)))
#define ArcTan(x)       (atan((x)))

#define Sinh(x)          (sinh((x)))
#define Cosh(x)          (cosh((x)))
#define Tanh(x)          (tanh((x)))

#define Cot(x)          (1./tan((x)))
#define Csc(x)          (1./sin((x)))




/** Could add definitions for Random(), SeedRandom(), etc. **/


