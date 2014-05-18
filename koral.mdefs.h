/*! \file koral.mdefs.h
  \brief  Mathematica source file
  C language definitions for use with Mathematica output
  Could add definitions for Random(), SeedRandom(), etc.
*/


#define Power(x, y) (pow((x), (y)))
#define Sqrt(x)  (sqrt((x)))

#define Abs(x)  (fabs((x)))

#define Exp(x)  (exp((x)))
#define Log(x)  (log((x)))

#define Sin(x)  (sin((x)))
#define Cos(x)  (cos((x)))
#define Tan(x)  (tan((x)))

#define ArcSin(x)       (asin((x)))
#define ArcCos(x)       (acos((x)))
#define ArcTan(x)       (atan((x)))

#define Sinh(x)          (sinh((x)))
#define Cosh(x)          (cosh((x)))
#define Tanh(x)          (tanh((x)))

#define Cot(x)          (1./tan((x)))
#define Csc(x)          (1./sin((x)))


#define Conjugate(x) (x) // assume not complex



