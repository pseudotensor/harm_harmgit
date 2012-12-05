//KORAL - bc.c
//returns problem specific BC
//**********************
//included in:
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) 
//from problem.c
//**********************

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
/**********************/
/**********************/
if(ix<0 &&  1)
  {

    ldouble rho=1.;
    ldouble uint,E,ux,ut,vx,Fx,Fy,Fz;
    if(NTUBE==1) {uint = 3.e-5 / (GAMMA - 1.); E=1.e-8; Fx=0.e-2*E;ux=0.015;}
    if(NTUBE==2) {uint = 4.e-3 / (GAMMA - 1.);E=2.e-5; Fx=0.e-2*E;ux=0.25;}
    if(NTUBE==3) {uint = 60. / (GAMMA - 1.);E=2.; Fx=0.e-2*E;ux=10.;}
    if(NTUBE==4) {uint = 6.e-3 / (GAMMA - 1.);E=0.18; Fx=0.e-2*E;ux=0.69;}	  
    ut=sqrtl(ux*ux+1.);
    vx=ux/ut;
    Fz=Fy=0.;
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

if(ix>=NX && 1)
  {    
    ldouble rho,uint,E,ux,ut,vx,Fx,Fy,Fz;

    if(NTUBE==1) {rho=2.4;uint = 1.61e-4/ (GAMMA - 1.); E=2.51e-7; Fx=0.e-2*E;ux=6.25e-3;}
    if(NTUBE==2) {rho=3.11;uint = 0.04512 / (GAMMA - 1.);E=3.46e-3; Fx=0.e-2*E;ux=0.0804;}
    if(NTUBE==3) {rho=8.0;uint = 2.34e3 / (GAMMA - 1.);E=1.14e3; Fx=0.e-2*E;ux=1.25;}
    if(NTUBE==4) {rho=3.65;uint =3.59e-2 / (GAMMA - 1.);E=1.30; Fx=0.e-2*E;ux=0.189;}	  
    ut=sqrtl(ux*ux+1.);
    vx=ux/ut;
    Fz=Fy=0.;
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

iix=ix;
iiz=iz;
iiy=iy;
//periodic
while(iiy<0)    iiy+=NY;
while(iiy>=NY)    iiy-=NY; 
//periodic
while(iiz<0)    iiz+=NZ;
while(iiz>=NZ)    iiz-=NZ; 
 
for(iv=0;iv<NV;iv++)
  {
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//converting to conserved
p2u(pp,uu,gg,eup,elo);
  
return 0;
  
