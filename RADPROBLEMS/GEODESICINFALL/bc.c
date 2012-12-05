//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);

/**********************/

  //radius
  //  if(ix>=NX || ix<0) //analytical solution on both sides
  if(ix>=NX) //analytical solution at rout only
    {
      //TODO: to lepiej prekalkulowac
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  
 
      ldouble gg[4][5],eup[4][4],elo[4][4];
      pick_g(ix,iy,iz,gg);
      pick_T(emuup,ix,iy,iz,eup);
      pick_T(emulo,ix,iy,iz,elo);

      ldouble r=get_x(ix,0);
      ldouble D=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
      ldouble E=PAR_E/(powl(r*r*sqrt(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
      ldouble V=sqrtl(2./r)*(1.-2./r);
      ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
      ldouble rho=D/W;
      ldouble uint=E/W;

      //corrected rho:
      rho=PAR_D/(r*r*sqrtl(2./r));    

      pp[0]=rho; pp[1]=uint; pp[2]=-V; pp[3]=pp[4]=0.; 
      pp[5]=calc_Sfromu(rho,uint);

      p2u(pp,uu,gg,eup,elo);

      return 0.;
    }
  else if(ix<0) //outflow near BH
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  
      ldouble gg[4][5],eup[4][4],elo[4][4],ggsrc[4][5];

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
      
      /*
      //analytical solution
      ldouble r=get_x(ix,0);
      ldouble D=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
      ldouble E=PAR_E/(powl(r*r*sqrt(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
      ldouble V=sqrtl(2./r)*(1.-2./r);
      ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
      ldouble rho=D/W;
      ldouble uint=E/W;
      pp[0]=rho; pp[1]=uint; pp[2]=-V; pp[3]=pp[4]=0.;
      */
      
      if(ix==-1) //conserved unneccesary for ix=-2 
	p2u(pp,uu,gg,eup,elo);
      return 0;
    }

  //reflections in theta 
  if(iy<0.) //spin axis
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
  if(iy>=NY) //equatorial plane
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
  while(iz<0)    iz+=NZ;
  while(iz>=NZ)    iz-=NZ; 

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,iix,iiy,iiz);
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  
  return 0;
  
