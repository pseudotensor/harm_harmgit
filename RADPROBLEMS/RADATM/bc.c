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
  if(ix>=NX || ix<0) //analytical solution at rout only
    {
      ldouble Fx,Fy,Fz;
      Fz=Fy=Fx=0.;
ldouble f = (ldouble)KAPPAES*FLUXLEFT*MINX*MINX;

ldouble p0=K_BOLTZ*RHOAMB*TAMB/MU_GAS/M_PROTON;	      
ldouble KKK=p0/powl(RHOAMB,GAMMA);
ldouble C3=GAMMA*KKK/(GAMMA-1.)*powl(RHOAMB,GAMMA-1.)-(1.-f)*(1./MINX+0.*1./MINX/MINX+0.*4./3./MINX/MINX/MINX);

pp[0]=powl(GAMMAM1/GAMMA/KKK*(C3+(1.-f)*(1./xx+0.*1./xx/xx+0.*4./3./xx/xx/xx)),1./GAMMAM1);

ldouble pre=KKK*powl(pp[0],GAMMA);

pp[1]=pre/GAMMAM1;
Fz=Fy=0.;
Fx=FLUXLEFT*(MINX/xx)*(MINX/xx);
ldouble E=Fx/FERATIO;

      pp[2]=0.;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(RHOAMB,pp[1]);
#ifdef RADIATION
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz; 
#endif	    
      p2u(pp,uu,gg,eup,elo);
     
      return 0.;
    }


iix=ix;
  iiz=iz;
  iiy=iy;

  //periodic
  while(iiz<0)    iiz+=NZ;
  while(iiz>=NZ)    iiz-=NZ; 
  //periodic
  while(iiy<0)    iiy+=NY;
  while(iiy>=NY)    iiy-=NY; 


  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,iix,iiy,iiz);
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  
  return 0;
  
