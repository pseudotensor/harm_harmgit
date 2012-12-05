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
      ldouble Fx,Fy,Fz,rho,rho0,Tgas0,E,uint,ur,Tgas,Trad,r,prad,pgas,ut,vx;

      //at outern boundary
      r=MAXX;
      ur=-sqrtl(2./r);
      rho0=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
      Tgas0=TGAS0;
            
      //at given cell
      r=xx;
      ur=-sqrtl(2./r);
      ut=sqrtl((-1.-ur*ur*gg[1][1])/gg[0][0]);
      vx=ur/ut;
      rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
      Tgas=Tgas0*pow(rho/rho0,GAMMA-1.);      
      uint=calc_PEQ_ufromTrho(Tgas,rho);
      pgas=K_BOLTZ*rho*Tgas/MU_GAS/M_PROTON;
      prad=PRADGAS*pgas;
      E=prad*3.;

      Fz=Fy=Fx=0.;
      pp[0]=rho;
      pp[1]=uint;
      pp[2]=vx;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;
#endif
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

      ldouble Fx,Fy,Fz,rho,rho0,Tgas0,E,uint,ur,Tgas,Trad,r,prad,pgas,ut,vx;

      //at outern boundary
      r=rbc;
      ur=-sqrtl(2./r);      
      ut=sqrtl((-1.-ur*ur*gg[1][1])/gg[0][0]);
      vx=ur/ut;

      //copying primitives with gdet taken into account
      for(iv=0;iv<NV;iv++)
	{ 
	  if(iv==2 || iv==7)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.-(rsrc-rbc)/(.5*(rsrc+rbc)));
	  else if(iv==3 || iv==4 || iv==8 || iv==9)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.+(rsrc-rbc)/(.5*(rsrc+rbc)));
	  else 
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);//gdet_src/gdet_bc;

	  //following ~r**-1.5 scaling
	  if(iv==0)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*powl(rsrc/rbc,1.5);
	  if(iv==1)
	    pp[iv]=get_u(p,iv,iix,iiy,iiz)*powl(rsrc/rbc,1.5*GAMMA);
	  if(iv==2)
	    pp[iv]=vx;
	  
	  //unchanged primitives
	  //	  pp[iv]=get_u(p,iv,iix,iiy,iiz);

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
      
      p2u(pp,uu,gg,eup,elo);
      return 0;
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
  
