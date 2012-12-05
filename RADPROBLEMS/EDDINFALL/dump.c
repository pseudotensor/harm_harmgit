/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

v1=Tgas;
v2=Trad;
ldouble tau[3];
calc_tautot(pp,xxx,dx,tau);
v3=tau[0];
calc_tauabs(pp,xxx,dx,tau);
v4=tau[0];

