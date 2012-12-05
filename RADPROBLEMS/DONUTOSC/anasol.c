//int
//analytical_solution(ldouble t,int ix,int iy,int iz,ldouble *uu,ldouble *pp,ldouble *vv)
//{

  ldouble gg[4][5],GG[4][5];
  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);
  ldouble podpierd=-(GG[0][0]-2.*ELL*GG[0][3]+ELL*ELL*GG[3][3]);
  ldouble ut=-1./sqrt(podpierd);
  ut/=UTPOT;
  ldouble uint,Vphi,rho;
  ldouble xx=get_x(ix,0);
  if(ut<-1. || podpierd<0. )
    {
      rho=RHOFLOOR;
      uint=UFLOOR;
      Vphi=0.;
    }
  else
    {
		       ldouble D,E,W,eps,uT,uphi,uPhi;
		       /*     uphi=-ELL*ut;
      uT=GG[0][0]*ut+GG[0][3]*uphi;
      eps=1./GAMMA*(-1./ut-1.);
      W=uT/sqrt(-gg[0][0]);
      D=pow(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.))*W;
      E=eps*D*W;
      rho=D/W;
      uint=E/W;
      uPhi=GG[3][3]*uphi+GG[0][3]*ut;
      Vphi=uPhi/uT;*/
		  ldouble h=-1./ut;
		  eps=(h-1.)/GAMMA;
		  rho=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
		  uint=rho*eps;
		  //		  uint=KKK*powl(rho,GAMMA)/(GAMMA-1.);
		  uphi=-ELL*ut;
		  uT=GG[0][0]*ut+GG[0][3]*uphi;
		  uPhi=GG[3][3]*uphi+GG[0][3]*ut;
		  Vphi=uPhi/uT;
		  
    }     

  vv[0]=rho; vv[1]=uint; vv[2]=Vphi; vv[3]=0.;
