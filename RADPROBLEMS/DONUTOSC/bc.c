//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);

/**********************/


  //radius
  if(ix>=NX) //analytical solution at rout only
    {
      ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];
      pick_g(ix,iy,iz,gg);
      pick_T(emuup,ix,iy,iz,eup);
      pick_T(emulo,ix,iy,iz,elo);
      pick_G(ix,iy,iz,GG);
      ldouble podpierd=-(GG[0][0]-2.*ELL*GG[0][3]+ELL*ELL*GG[3][3]);
      ldouble ut=-1./sqrt(podpierd);
      ut/=UTPOT;
      ldouble uint,Vphi,rho,Vr;
      ldouble xx=get_x(ix,0);
      ldouble D,E,W,eps,uT,uphi,uPhi;
      if(ut<-1 || podpierd<0.|| NODONUT)
	{
	  rho=RHO_AMB*pow(xx/2.,-1.5);
	  uint=U_AMB*pow(xx/2.,-5./2.);
	  Vphi=0.;
	  Vr=0.;
	   
	  ldouble r=get_x(ix,0);
	  D=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
	  E=PAR_E/(powl(r*r*sqrt(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
	  ldouble V=sqrtl(2./r)*(1.-2./r);
	  W=1./sqrtl(1.-V*V*gg[1][1]);
	  rho=D/W;
	  uint=E/W;
	  Vr=V;
	 
	  rho=PAR_D/(r*r*sqrtl(2./r));  
	  uPhi=GG[3][3]*uphi+GG[0][3]*ut;

	  //	  Vr=0.;
	  //	  Vphi=0.*uPhi/uT;
	  
	  //corrected rho:
  	  uT=GG[0][0]*ut+GG[0][3]*uphi;

	}
      else
	{
		  ldouble h=-1./ut;
		  ldouble eps=(h-1.)/GAMMA;
		  rho=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
		  uint=rho*eps;
		  //		  uint=KKK*powl(rho,GAMMA)/(GAMMA-1.);
		  uphi=-ELL*ut;
		  uT=GG[0][0]*ut+GG[0][3]*uphi;
		  uPhi=GG[3][3]*uphi+GG[0][3]*ut;
		  Vphi=uPhi/uT;
		  Vr=0.;

		  /*	  uphi=-ELL*ut;
	  uT=GG[0][0]*ut+GG[0][3]*uphi;
	  eps=1./GAMMA*(-1./ut-1.);
	  W=uT/sqrt(-gg[0][0]);
	  D=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.))*W;
	  E=eps*D*W;
	  rho=D/W;
	  uint=E/W;
	  uPhi=GG[3][3]*uphi+GG[0][3]*ut;
	  Vphi=uPhi/uT;
	  Vr=0.;*/
	}     

      pp[0]=rho; pp[1]=uint; pp[4]=Vphi; pp[2]=-Vr; pp[3]=0.;
      pp[5]=calc_Sfromu(rho,uint);

      p2u(pp,uu,gg,eup,elo);

      return 0.;
    }
  else if(ix<0) //outflow near BH
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  
      ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4];

      iix=0;
      iiy=iy;
      iiz=iz;
      pick_g(ix,iy,iz,gg);
      pick_T(emuup,ix,iy,iz,eup);
      pick_T(emulo,ix,iy,iz,elo);
      pick_g(iix,iiy,iiz,ggsrc);
      gdet_src=get_g(g,3,4,iix,iiy,iiz);  
      gdet_bc=get_g(g,3,4,ix,iy,iz);        
      ldouble rsrc=get_x(iix,0);
      ldouble rbc=get_x(ix,0);
      //copying primitives with gdet taken into account
      for(iv=0;iv<NV;iv++)
	{ 
	  if(iv==0 || iv==1 || iv==5)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*gdet_src/gdet_bc;
	  if(iv==2)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.-(rsrc-rbc)/(.5*(rsrc+rbc)));
	  if(iv==3 || iv==4)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.+(rsrc-rbc)/(.5*(rsrc+rbc)));
	  
	  //unchanged primitives
	  //pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}

      /*
      //free-fall acc.to analytical solution from the ix=0 cell
      ldouble Vinfall=sqrtl(2./rbc)*(1.-2./rbc);
      ldouble Dratio=((rsrc*rsrc*sqrtl(2./rsrc*(1.-2./rsrc)))) /
	((rbc*rbc*sqrtl(2./rbc*(1.-2./rbc))));
      ldouble Eratio=(powl(rsrc*rsrc*sqrt(2./rsrc),GAMMA)*powl(1.-2./rsrc,(GAMMA+1.)/4.)) /
	(powl(rbc*rbc*sqrt(2./rbc),GAMMA)*powl(1.-2./rbc,(GAMMA+1.)/4.));
      ldouble rhosrc=get_u(p,0,iix,iiy,iiz);
      ldouble uintsrc=get_u(p,1,iix,iiy,iiz);
      ldouble Vsrc=get_u(p,2,iix,iiy,iiz);
      ldouble Wsrc=1./sqrtl(1.-Vsrc*Vsrc*ggsrc[1][1]);
      ldouble Dsrc=rhosrc*Wsrc;
      ldouble Esrc=uintsrc*Wsrc;
      ldouble V=Vsrc*(sqrtl(2./rbc)*(1.-2./rbc))/(sqrtl(2./rsrc)*(1.-2./rsrc));
      ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
      ldouble D=Dsrc*Dratio;
      ldouble E=Esrc*Eratio;
      ldouble rho=D/W;
      ldouble uint=E/W;      
      pp[0]=rho; pp[1]=uint; pp[2]=V; pp[3]=pp[4]=0.;
	  */
	
      //      pp[5]=calc_Sfromu(pp[0],pp[1]);

      if(ix==-1) //conserved unneccesary for ix=-2 
	p2u(pp,uu,gg,eup,elo);
      return 0;
    }

  //reflections in theta 
  if(iy<0.) //upper
    {      
      iiy=-iy-1;
      iiz=iz;
      iix=ix;
      gdet_src=get_g(g,3,4,iix,iiy,iiz);  
      gdet_bc=get_g(g,3,4,ix,iy,iz);  
      ldouble gg[4][5],eup[4][4],elo[4][4];
      pick_g(ix,iy,iz,gg);
      pick_T(emuup,ix,iy,iz,eup);
      pick_T(emulo,ix,iy,iz,elo);
       for(iv=0;iv<NV;iv++)
	{
	  if(iv==3)
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	  else
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      p2u(pp,uu,gg,eup,elo);
      return 0;
     }
  if(iy>=NY) //lower
    {
      iiy=NY-(iy-NY)-1;
      iiz=iz;
      iix=ix;
      gdet_src=get_g(g,3,4,iix,iiy,iiz);  
      gdet_bc=get_g(g,3,4,ix,iy,iz);  
       ldouble gg[4][5],eup[4][4],elo[4][4];
      pick_g(ix,iy,iz,gg);
      pick_T(emuup,ix,iy,iz,eup);
      pick_T(emulo,ix,iy,iz,elo);
 	  
      for(iv=0;iv<NV;iv++)
	  {
	    if(iv==3)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      p2u(pp,uu,gg,eup,elo); 
      return 0; 
    }
   
  //symmetric in phi:
  iiz=iz;
  iiy=iy;
  iix=ix;
  if(iz<0) iiz=iz+NZ;
  if(iz>NZ-1) iiz=iz-NZ;

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,iix,iiy,iiz);
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  
  return 0;

