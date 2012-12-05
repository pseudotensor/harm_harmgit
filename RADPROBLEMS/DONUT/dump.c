/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

 

 analytical_solution(t,ix,iy,iz,uu,pp,vv);
v1=vv[0];
v2=vv[1];
v3=vv[2];
v4=vv[3];
ldouble GG[4][5];
pick_G(ix,iy,iz,GG);
ldouble podpierd=-(GG[0][0]-2.*ELL*GG[0][3]+ELL*ELL*GG[3][3]);
ldouble ulowert=-1./sqrt(podpierd);
if(ulowert<-1. || podpierd<0.)
  {
    v4=0.  ;
  }
 else
   {
     v4=ulowert;
   }
