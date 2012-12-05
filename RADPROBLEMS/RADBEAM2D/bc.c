//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

  gdet_bc=get_g(g,3,4,ix,iy,iz);  
  ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
  pick_g(ix,iy,iz,gg);
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);
  ldouble xx=get_x(ix,0);

//printf("aa\n");

  //radius
  if(iz<0 && xx>BEAML && xx<BEAMR && IFBEAM) //hot boundary
    {
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      iix=ix;
      iiy=iy;
      iiz=0;
      rho=RHOAMB;
      E=calc_LTE_EfromT(TLEFT);
      uint=get_u(p,1,iix,iiy,iiz);
      rho=get_u(p,0,iix,iiy,iiz);
      vx=get_u(p,2,iix,iiy,iiz);
      //      E=calc_LTE_Efromurho(uint,rho);

      Fz=NLEFT*E;
      Fy=Fx=0.;

      //      Fx=Fz/sqrt(2.);
      //      Fz=Fz/sqrt(2.);

      pp[0]=rho;
      pp[1]=uint;
      pp[2]=vx;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;

      prad_zamo2ff(pp,pp,gg,eup);

      p2u(pp,uu,gg,eup,elo);
      return 0.;
    }
  else if(iz<0 )
    {
      /*
      ldouble Fx,Fy,Fz,rho,E,uint,vx,T;
      rho=RHOAMB;
      T=TAMB;
      E=calc_LTE_EfromT(T);
      uint=calc_PEQ_ufromTrho(T,rho);

      Fx=Fy=Fz=0.*E;
      pp[0]=rho;
      pp[1]=uint;
      pp[2]=0;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;
      p2u(pp,uu,gg,eup,elo);
      return 0;
      */
      
      
      iiz=0;
      iiy=iy;
      iix=ix;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      //flux azimuthal only
      pp[7]=0.;
      pp[9]=0.;

      p2u(pp,uu,gg,eup,elo);


      return 0;
      
    }
  
  else if(iz>=NZ ) //copy
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=ix;
      iiy=iy;
      iiz=NZ-1;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      p2u(pp,uu,gg,eup,elo);
      return 0;
    }
  else if(ix<0 && 1) //copy
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=0;
      iiy=iy;
      iiz=iz;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      //       rho=RHOAMB;
       //      E=calc_LTE_EfromT(TAMB);
      //pp[6]=E;

      p2u(pp,uu,gg,eup,elo);
      return 0;
    }
#ifdef FLATBACKGROUND
   else if(ix>NX-1 || ix<0) //copy
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=NX-1;
      iiy=iy;
      iiz=iz;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      //       rho=RHOAMB;
       //      E=calc_LTE_EfromT(TAMB);
      //pp[6]=E;
      p2u(pp,uu,gg,eup,elo);
      return 0;
    }
#endif
   else if(ix>NX-1 || 1 ) //fixed outern boundary
    {
      iix=NX-1;
      iiy=iy;
      iiz=iz;
      //zaczynam jednak od profilu analitycznego:   
      ldouble r=get_x(ix,0);
      ldouble mD=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
      ldouble mE=PAR_E/(powl(r*r*sqrtl(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
      ldouble V=sqrtl(2./r)*(1.-2./r)           ;
      ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
      ldouble rho=PAR_D/(r*r*sqrtl(2./r));
      ldouble T=TAMB;
      ldouble E=calc_LTE_EfromT(T);
      ldouble uint=mE/W;
      ldouble Fx,Fy,Fz;
      E=calc_LTE_Efromurho(uint,rho);
      Fx=Fy=Fz=0.*E;
      pp[0]=rho;
      pp[1]=uint;
      pp[2]=-V;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=pp[8]=pp[9]=0.;
      //copy of radiative primitives
      /*
      pp[6]=get_u(p,6,iix,iiy,iiz);
      pp[7]=get_u(p,7,iix,iiy,iiz);
      pp[8]=get_u(p,8,iix,iiy,iiz);
      pp[9]=get_u(p,9,iix,iiy,iiz);
      */
     prad_zamo2ff(pp,pp,gg,eup);
  p2u(pp,uu,gg,eup,elo);
      return 0;
    }
 


  iix=ix;
  iiz=iz;
  iiy=iy;
  //copy
  while(iix<0)    iix=0.;//iix+=NX;
  while(iix>=NX)    iix=NX-1;//iix-=NX; 
  //periodic
  while(iiy<0)    iiy+=NY;
  while(iiy>=NY)    iiy-=NY; 
 
  for(iv=0;iv<NV;iv++)
    {
      
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  p2u(pp,uu,gg,eup,elo);
return 0;
  
