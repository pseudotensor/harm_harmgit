//ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
//{
  ldouble  kff=kappaCGS2GU(1.7e-25/massGU2CGS(M_PROTON)/massGU2CGS(M_PROTON)* powl(T,-3.5)*rhoGU2CGS(rho));
  kff*=rho;  
  return kff;

