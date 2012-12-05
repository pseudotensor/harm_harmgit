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

ldouble nx,ny,nz;
  nx=Fx/E;
  ny=Fy/E;
  nz=Fz/E;

 ldouble nlen=sqrtl(nx*nx+ny*ny+nz*nz);
ldouble f;ldouble Rij[4][4];
 
#ifdef EDDINGTON_APR
  f=1./3.;
#else  
  if(nlen>=1.)
    {
      f=1.;
    }
  else //M1
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrtl(4.-3.*(nx*nx+ny*ny+nz*nz)));  
#endif
  
  if(nlen>0) 
    {
      nx/=nlen;
      ny/=nlen;
      nz/=nlen;
    }
  else
    {
      ;
    }

 
  Rij[0][0]=E;
  Rij[0][1]=Rij[1][0]=Fx;
  Rij[0][2]=Rij[2][0]=Fy;
  Rij[0][3]=Rij[3][0]=Fz;

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

v3=Rij[1][1]/E;
v4=Rij[2][2]/E;

v3=f;
