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
  if(ix>=NX) //analytical solution at rout only
    {
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
//mrho=PAR_D/(r*r*sqrtl(2./r));
rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(-sqrtl(2./r))));
uint=E/D*rho;
      

      pp[0]=rho; pp[1]=uint; pp[2]=-V; pp[3]=pp[4]=0.; 
      pp[5]=calc_Sfromu(rho,uint);

      //ldouble Erad=calc_LTE_Efromurho(uint,rho);
ldouble Erad=calc_LTE_EfromT(TAMB);
pp[6]=Erad;
pp[7]=pp[8]=pp[9]=0.;


// prad_zamo2ff(pp,pp,gg,eup);
      p2u(pp,uu,gg,eup,elo);

      return 0.;
    }
  else if(ix<0) //outflow near BH
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  
      ldouble gg[4][5],eup[4][4],elo[4][4],ggsrc[4][5],pp[NV];

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
	  else if(iv==2)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.-(rsrc-rbc)/(.5*(rsrc+rbc)));
	  else if(iv==3 || iv==4)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.+(rsrc-rbc)/(.5*(rsrc+rbc)));
	  else	 
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);

	  pp[iv]=pp[iv];
	}

      //irradiating flux
      ldouble r=rbc;
      ldouble Flux=fluxCGS2GU(LUM*LUMEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)));
      pp[6]=1.01*Flux;
      pp[7]=Flux;
      //    prad_zamo2ff(pp,pp,gg,eup);
      pp[6]=pp[6];
      pp[7]=pp[7];
     
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
	  if(iv==3 || iv==8)
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
	    if(iv==3 || iv==8)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      //       prad_zamo2ff(pp,pp,gg,eup);
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
  
