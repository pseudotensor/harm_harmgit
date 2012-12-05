//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

/**********************/
  gdet_bc=get_g(g,3,4,ix,iy,iz);  
  ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
  pick_g(ix,iy,iz,gg);
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);
  ldouble xx=get_x(ix,0);

if(ix<0 ) iix=0;
if(ix>NX-1) iix=NX-1;


  //periodic
  while(iiz<0)    iiz+=NZ;
  while(iiz>=NZ)    iiz-=NZ; 
  //periodic
  while(iiy<0)    iiy+=NY;
  while(iiy>=NY)    iiy-=NY; 


  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }

p2u(pp,uu,gg,eup,elo);

  
  return 0;
  
