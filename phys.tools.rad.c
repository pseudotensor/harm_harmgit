#include "decs.h"

void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G);
void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu);
void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM]);

//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame, returns ultimate deltas ***********************
//******* the fiducial approach *****************************************
//**********************************************************************


//uu0 - original cons. qty
//uu -- current iteration
//f - (returned) errors

//NOTE: uu WILL be changed from inside this fcn
int f_implicit_lab(FTYPE *pp0, FTYPE *uu0,FTYPE *uu,FTYPE realdt, struct of_geom *ptrgeom,  FTYPE *f)
{
  struct of_state q;
  FTYPE pp[NPR];
  int pliter, pl;
  int iv;
  struct of_newtonstats newtonstats;
  int finalstep = 1;  //can choose either 1 or 0 depending on whether want floor-like fixups (1) or not (0).  unclear which one would work best since for Newton method to converge might want to allow negative density on the way to the correct solution, on the other hand want to prevent runaway into rho < 0 region and so want floors.

  // initialize counters
  newtonstats.nstroke=newtonstats.lntries=0;

  PLOOP(pliter,pl)
    pp[pl] = pp0[pl];

  DLOOPA(iv)
    uu[UU+iv] = uu0[UU+iv] - (uu[URAD0+iv]-uu0[URAD0+iv]);
  
  //calculating primitives  
  int corr;
  //u2p(uu,pp,gg,GG,&corr);

  MYFUN(Utoprimgen(finalstep, EVOLVEUTOPRIM, UNOTHING, uu, ptrgeom, pp, &newtonstats),"phys.tools.rad.c:f_implicit_lab()", "Utoprimgen", 1);


  get_state_uconucovonly(pp, ptrgeom, &q);
  get_state_uradconuradcovonly(pp, ptrgeom, &q);

  //radiative covariant four-force
  FTYPE Gd[NDIM];
  calc_Gd(pp, ptrgeom, &q, Gd);
  
  DLOOPA(iv)
    f[iv]=uu[URAD0+iv] - uu0[URAD0+iv] + realdt * Gd[iv];

  return 0;
} 

//SUPERGODMARK only for RK2
FTYPE compute_dt()
{
  if( steppart == 0 ) {
    return( 0.5 * dt );
  }
  else if( steppart == 1 ) {
    return( dt );
  }
  else {
    dualfprintf(fail_file,"Should not get here: compute_dt()\n");
    myexit(12112);
  }

  return(0.0);

}

// compute changes to U (both T and R) using implicit method
void koral_implicit_source_rad(FTYPE *pin, FTYPE *Uin, struct of_geom *ptrgeom, struct of_state *q ,FTYPE (*dUcomp)[NPR])
{
  FTYPE compute_dt();
  int i1,i2,i3,iv,ii,jj,pliter,sc;
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE uu0[NPR],uup[NPR],uu[NPR]; 
  FTYPE f1[NDIM],f2[NDIM],f3[NDIM],x[NDIM];
  FTYPE realdt;
  FTYPE radsource[NPR], deltas[NDIM]; 
  int pl;

  realdt = compute_dt();
 
  //uu0 will hold original vector of conserved
  PLOOP(pliter,iv)
    uu[iv] = uu0[iv] = Uin[iv];
  
  FTYPE EPS = 1.e-6;
  FTYPE CONV = 1.e-6 ;
  
  int iter=0;
  
  do
  {
    iter++;
    
    //vector of conserved at the previous iteration
    PLOOP(pliter,ii)
      uup[ii]=uu[ii];
    
    //values at zero state
    f_implicit_lab(pin, uu0, uu, realdt, ptrgeom, f1);
    
    //calculating approximate Jacobian
    DLOOPA(ii)
      DLOOPA(jj)
      {
	FTYPE del;
	if(uup[jj+URAD0]==0.) 
	  del=EPS*uup[URAD0];
	else 
	  del=EPS*uup[jj+URAD0];

	uu[jj+URAD0]=uup[jj+URAD0]-del;
	
	f_implicit_lab(pin,uu0,uu,realdt,ptrgeom,f2);
	
	J[ii][jj]=(f2[ii] - f1[ii])/(uu[jj+URAD0]-uup[jj+URAD0]);
	
	uu[jj+URAD0]=uup[jj+URAD0];
      }
    
    
    //inversion
    inverse_44matrix(J,iJ);
    
    //updating x
    DLOOPA(ii){
      x[ii]=uup[ii+URAD0];
    }
    
    DLOOPA(ii)
      DLOOPA(jj)
      {
	x[ii]-=iJ[ii][jj]*f1[jj];
      }
    
    DLOOPA(ii)
      uu[ii+URAD0]=x[ii];
    
    //test convergence
    DLOOPA(ii){
      f3[ii]=(uu[ii+URAD0]-uup[ii+URAD0]);
      f3[ii]=fabs(f3[ii]/uup[URAD0]);
    }
    
    if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
      break;
    
    if(iter>50)
    {
      dualfprintf(fail_file,"iter exceeded in solve_implicit_lab()\n");	  
      myexit(21341);
    }
    
  }
  while(1);
  
  DLOOPA(jj){
    deltas[jj]=(uu[URAD0+jj]-uu0[URAD0+jj])/dt;
  }

  PLOOP(pliter,pl){
    radsource[pl] = 0;
  }

  DLOOPA(jj) radsource[UU+jj] = -deltas[jj];
  DLOOPA(jj) radsource[URAD0+jj] = deltas[jj];

  sc = RADSOURCE;

  PLOOP(pliter,pl){
    dUcomp[sc][pl] += radsource[pl];
  }
  
}


// compute changes to U (both T and R) using implicit method
// NOTEMARK: The explicit scheme is only stable if the fluid speed is order the speed of light.  Or, one can force the explicit scheme to use vrad=c.
void koral_explicit_source_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q ,FTYPE (*dUcomp)[NPR])
{
  FTYPE Gd[NDIM], radsource[NPR];
  int pliter, pl, jj, sc;

  calc_Gd(pr, ptrgeom, q, Gd);

  sc = RADSOURCE;
  
  PLOOP(pliter,pl){
    radsource[pl] = 0;
  }

  DLOOPA(jj) radsource[UU+jj] = -Gd[jj];
  DLOOPA(jj) radsource[URAD0+jj] = Gd[jj];

  PLOOP(pliter,pl){
    dUcomp[sc][pl] += radsource[pl];
  }
  
}

void inline koral_source_rad(FTYPE *pin, FTYPE *Uin, struct of_geom *ptrgeom, struct of_state *q ,FTYPE (*dUcomp)[NPR])
{
#if(WHICHRADSOURCEMETHOD==RADSOURCEMETHODEXPLICIT)
  koral_explicit_source_rad( pin, ptrgeom, q, dUcomp);
#elif(WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICIT)
  koral_implicit_source_rad( pin, Uin, ptrgeom, q, dUcomp);
#elif(WHICHRADSOURCEMETHOD==RADSOURCEMETHODNONE)
#endif
}


//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************
//absorption
void calc_kappa(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa)
{

#if(1)
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappa()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  calc_kappa_user(rho,T,xx,yy,zz);
#else
  *kappa = 0.;
#endif  
}

//scattering
void calc_kappaes(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *kappa)
{  
#if(1)
  extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  //user_calc_kappaes()
  FTYPE rho=pr[RHO];
  FTYPE u=pr[UU];
  int ii=ptrgeom->i;
  int jj=ptrgeom->j;
  int kk=ptrgeom->k;
  int loc=ptrgeom->p;
  FTYPE T=compute_temp_simple(ii,jj,kk,loc,rho,u);
  FTYPE V[NDIM],xx,yy,zz;
  bl_coord_ijk(ii,jj,kk,loc,V);
  xx=V[1];
  yy=V[2];
  zz=V[3];
  calc_kappaes_user(rho,T,xx,yy,zz);
#else
  *kappa = 0.;
#endif  
}


void calc_Gd(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *G) 
{
  calc_Gu(pp, ptrgeom, q, G);
  indices_21(G, G, ptrgeom);
}


//**********************************************************************
//****** takes radiative stress tensor and gas primitives **************
//****** and calculates contravariant four-force ***********************
//**********************************************************************
void calc_Gu(FTYPE *pp, struct of_geom *ptrgeom, struct of_state *q ,FTYPE *Gu) 
{
  int i,j,k;
  
  //radiative stress tensor in the lab frame
  FTYPE Rij[NDIM][NDIM];

  //this call returns R^i_j, i.e., the first index is contra-variant and the last index is co-variant
  mhdfull_calc_rad(pp, ptrgeom, q, Rij);
  
  //the four-ve locity of fluid in lab frame
  FTYPE *ucon,*ucov;

  ucon = q->ucon;
  ucov = q->ucov;
  
  
  //gas properties
#if(0)
  FTYPE rho=pp[RHO];
  FTYPE u=pp[UU];
  FTYPE p= pressure_rho0_u_simple(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,rho,u);
  FTYPE T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  FTYPE B = SIGMA_RAD*pow(T,4.)/Pi;
  FTYPE Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
#else
  FTYPE B=1.0; // new HARM units?  KORALTODO SUPERGODMARK
  // Also, doesn't this only allow for thermal black body emission?
#endif
  FTYPE kappa,kappaes;
  calc_kappa(pp,ptrgeom,&kappa);
  calc_kappaes(pp,ptrgeom,&kappaes);
  FTYPE chi=kappa+kappaes;
  
  //contravariant four-force in the lab frame
  
  //R^ab u_a u_b
  FTYPE Ruu=0.;
  DLOOP(i,j)
    Ruu+=Rij[i][j]*ucov[i]*ucon[j];
  
  FTYPE Ru;
  DLOOPA(i)
  {
    Ru=0.;
    DLOOPA(j)
      Ru+=Rij[i][j]*ucon[j];
    Gu[i]=-chi*Ru - (kappaes*Ruu + kappa*4.*Pi*B)*ucon[i];
  }
}

int vchar_all(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxall, FTYPE *vminall,int *ignorecourant)
{
  FTYPE vminmhd,vmaxmhd;
  FTYPE vminrad,vmaxrad;
  
  vchar_each(pr, q, dir, geom, &vmaxmhd, &vminmhd, &vmaxrad, &vminrad, ignorecourant);
  // below correct even if EOMRADTYPE==EOMRADNONE because vchar_each() sets both mhd and rad to be mhd and so below always chooses the mhd values.
  *vminall=MIN(vminmhd,vminrad);
  *vmaxall=MAX(vmaxmhd,vmaxrad);
  

  return(0);
}

int vchar_each(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmaxmhd, FTYPE *vminmhd, FTYPE *vmaxrad, FTYPE *vminrad,int *ignorecourant)
{
  
  vchar(pr, q, dir, geom, vmaxmhd, vminmhd,ignorecourant);
  if(EOMRADTYPE!=EOMRADNONE){
    vchar_rad(pr, q, dir, geom, vmaxrad, vminrad,ignorecourant);
  }
  else{// default as if no other values for wave speeds
    *vmaxrad=*vmaxmhd;
    *vminrad=*vminmhd;
  }

  return(0);
}

int vchar_rad(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin,int *ignorecourant)
{
  extern int simplefast(int dir, struct of_geom *geom,struct of_state *q, FTYPE cms2,FTYPE *vmin, FTYPE *vmax);

  
  // compute kappa
  // Assume computed as d\tau/dorthonormallength as defined by user.
  // Assume \kappa defined in fluid frame (i.e. not radiation frame).
  FTYPE kappa;
  calc_kappa(pr,geom,&kappa);

  //characterisitic wavespeed in the radiation rest frame
  FTYPE vrad2=THIRD;
  FTYPE vrad2limited;

  if(kappa>0.0){// && WHICHRADSOURCEMETHOD==RADSOURCEMETHODIMPLICIT){
    // NOT DOING THIS:
    // compute tautot assuming kappa is optical depth per unit grid dx[1-3].  I.e. calc_kappa() computes grid-based opacity
    // tautot is the total optical depth of the cell in dim dimension
    //  tautot = kappa * dx[dir];

    // DOING THIS:
    // KORALTODO: Approximation to any true path, but approximation is sufficient for approximate wave speeds.
    // \tau_{\rm tot}^2 \approx \kappa^2 [dx^{dir} \sqrt{g_{dirdir}}]^2 
    FTYPE tautotsq,vrad2tau;
    // Note that tautot is frame independent once multiple \kappa by the cell length.  I.e. it's a Lorentz invariant.
    tautotsq = kappa*kappa * dx[dir]*dx[dir]*(geom->gcov[GIND(dir,dir)]);
  
    vrad2tau=(4.0/3.0)*(4.0/3.0)/tautotsq; // KORALTODO: Why 4.0/3.0 ?
    vrad2limited=MIN(vrad2,vrad2tau);
  }
  else{
    vrad2limited=vrad2;
  }


  //need to substitute ucon,ucov with uradcon,uradcov to fool simplefast
  FTYPE ucon[NDIM],ucov[NDIM];
  int ii;
  DLOOPA(ii){
    ucon[ii]=q->ucon[ii];
    ucov[ii]=q->ucov[ii];
    q->ucon[ii]=q->uradcon[ii];
    q->ucov[ii]=q->uradcov[ii];
  }

  //calculating vmin, vmax
  simplefast(dir,geom,q,vrad2limited,vmin,vmax);

#if(0)
  // Cartesian-Minkowski speed-of-light limit of radiation velocity
 FTYPE dxdxp[NDIM][NDIM];
  dxdxprim_ijk(geom->i, geom->j, geom->k, geom->p, dxdxp);
  // characeristic wavespeeds are 3-velocity in lab-frame
  *vmin=-1.0/dxdxp[dir][dir]; // e.g. dxdxp = dr/dx1
  *vmax=+1.0/dxdxp[dir][dir];
#endif

  //restoring gas 4-velocities
  DLOOPA(ii){
    q->ucon[ii]=ucon[ii];
    q->ucov[ii]=ucov[ii];
  }
  
  return(0);
}

// convert primitives to conserved radiation quantities
int p2u_rad(FTYPE *pr, FTYPE *Urad, struct of_geom *ptrgeom, struct of_state *q)
{
  FTYPE Uraddiag;
  
  mhd_calc_rad(pr, TT, ptrgeom, q, Urad);

  // SUPERGODMARK: ensure \detg applied when filling conserved quantity

  return(0);
}

// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_uradconuradcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  void compute_1plusud0(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])

  // urad^\mu
  // ucon_calc() assumes primitive velocities are in U1 through U3, but otherwise calculation is identical for radiation velocity, so just shift effective list of primitives so ucon_calc() operates on U1RAD through U3RAD
  MYFUN(ucon_calc(&pr[URAD1-U1], ptrgeom, q->uradcon,q->othersrad) ,"phys.c:get_state()", "ucon_calc()", 1);
  // urad_\mu
  lower_vec(q->uradcon, ptrgeom, q->uradcov);


  return (0);
}


void mhdfull_calc_rad(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE (*radstressdir)[NDIM])
{
  int jj,kk;

  if(EOMRADTYPE!=EOMRADNONE){
	  DLOOPA(jj) {
  	  mhd_calc_rad( pr, jj, ptrgeom, q, &(radstressdir[jj][0]) );
  	}  
	}
  else DLOOP(jj,kk) radstressdir[jj][kk]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.
}

// compute radiation stres-energy tensor
void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *ptrgeom, struct of_state *q, FTYPE *radstressdir)
{
  int jj;

  if(EOMRADTYPE!=EOMRADNONE){
    // R^{dir}_{jj} radiation stress-energy tensor
    DLOOPA(jj) radstressdir[jj]=THIRD*(4.0*pr[PRAD0]*q->uradcon[dir]*q->uradcov[jj] + pr[PRAD0]*delta(dir,jj));
  }
  else DLOOPA(jj) radstressdir[jj]=0.0; // mhd_calc_rad() called with no condition in phys.tools.c and elsewhere, and just fills normal tempo-spatial components (not RAD0->RAD3), so need to ensure zero.


}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* takes E and F^i from primitives and calculates radiation stress ****
//******* tensor R^ij using M1 closure scheme *****************************
//**********************************************************************
// Use: For initial conditions and dumping
// pp : fluid frame radiation conserved quantities
// Rij : fluid frame radiation stress-energy tensor
int calc_Rij_ff(FTYPE *pp, FTYPE Rij[][NDIM])
{
  FTYPE E=pp[PRAD0];
  FTYPE F[3]={pp[PRAD1],pp[PRAD2],pp[PRAD3]};

  FTYPE nx,ny,nz,nlen,f;

  nx=F[0]/E;
  ny=F[1]/E;
  nz=F[2]/E;

  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
#ifdef EDDINGTON_APR
  f=1./3.;
#else  
  if(nlen>=1.)
    {
      f=1.;
    }
  else //M1
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
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
  Rij[0][1]=Rij[1][0]=F[0];
  Rij[0][2]=Rij[2][0]=F[1];
  Rij[0][3]=Rij[3][0]=F[2];

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}



FTYPE my_min(FTYPE a, FTYPE b)
{
  if(a<b) return a;
  else return b;
}

FTYPE my_sign(FTYPE x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 0.;
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 4by4 matrix
int inverse_44matrix(FTYPE a[][NDIM], FTYPE ia[][NDIM])
{
  FTYPE mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  FTYPE	tmp[12]; FTYPE	src[16]; FTYPE det;
  /* transpose matrix */
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

 
  /* calculate matrix inverse */
  det = 1.0/det; 

  if(isnan(det)){
    dualfprintf(fail_file,"det in inverse 4x4 zero\n");
    myexit(13235);
  }

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      ia[i][j]= dst[i*4+j];

  return 0;
}



/*********************************************************************************/
/****** radiative ortonormal ff primitives (E,F^i) <-> primitives in lab frame  *******/
/*********************************************************************************/
int prad_fforlab(int whichdir, FTYPE *pp1, FTYPE *pp2,struct of_state *q, struct of_geom *ptrgeom)
{
  FTYPE Rij[NDIM][NDIM];
  int pliter,pl;
  int boost22_laborff(int whichdir, FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom);

  //TODO : no distinction basing on whichdir!

  //radiative stress tensor in the fluid frame
  calc_Rij_ff(pp1,Rij);

  //transforming from ortonormal coordinates to code coordinates
  //trans22_on2cc(Rij,Rij,tlo);

  //boosting the tensor to the lab frame
  boost22_laborff(whichdir,Rij,Rij,pp1,q,ptrgeom);

  //R^munu -> R^mu_nu
  indices_2221(Rij,Rij,ptrgeom);

  PLOOP(pliter,pl)
    pp2[pl]=pp1[pl];

  //temporarily store conserved in pp2[]
  pp2[PRAD0]=Rij[0][0];
  pp2[PRAD1]=Rij[0][1];
  pp2[PRAD2]=Rij[0][2];
  pp2[PRAD3]=Rij[0][3];

  //convert to real primitives - conversion does not care about MHD only about radiative conserved
  PFTYPE lpflag=UTOPRIMNOFAIL,lpflagrad=UTOPRIMRADNOFAIL;
  // NOTEMARK: lpflag=UTOPRIMNOFAIL means accept input pp2 primitives for velocity to maybe be used in local reductions to fluid frame.
  u2p_rad(pp2,pp2,ptrgeom, &lpflag, &lpflagrad);

  return 0;
} 


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//calculates general Lorenz matrix for lab <-> ff
//whichdir: [LAB2FF, FF2LAB]
// SUPERGODMARK: What is lab frame velocity?
int
calc_Lorentz_laborff(int whichdir,FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE L[][NDIM])
{
  int ii,jj,kk;
  int verbose=0;

  //"to" frame 4-velocity
  FTYPE ucon[NDIM],ucov[NDIM];
  //"from" frame 4-velocity
  FTYPE wcon[NDIM],wcov[NDIM];

  if(whichdir==LAB2FF)
    {
      //fluid frame
      DLOOPA(ii) 
      {
	ucon[ii]=q->ucon[ii];
	ucov[ii]=q->ucov[ii];
      }
      
      //lab frame
      //TODO: KerrShild???
      wcon[0]=1./sqrt(-ptrgeom->gcov[GIND(0,0)]);
      wcon[1]=wcon[2]=wcon[3]=0.;
      indices_21(wcon,wcov,ptrgeom);
    }
  else if(whichdir==FF2LAB)
    {
      //fluid frame
      DLOOPA(ii) 
      {
	wcon[ii]=q->ucon[ii];
	wcov[ii]=q->ucov[ii];
      }

      //lab frame
      ucon[0]=1./sqrt(-ptrgeom->gcov[GIND(0,0)]);
      ucon[1]=ucon[2]=ucon[3]=0.;
      indices_21(ucon,ucov,ptrgeom);
    }

  //temporary Om matrix
  FTYPE Om[NDIM][NDIM];

  DLOOP(ii,jj)
      Om[ii][jj]=ucon[ii]*wcov[jj]-wcon[ii]*ucov[jj];
  
  //Lorentz factor = -w^mu u_mu
  FTYPE gam=0.;
  DLOOPA(ii)
    gam+=wcon[ii]*ucov[ii];

  FTYPE Omsum;
  //Lorentz matrix components
  DLOOP(ii,jj)
    {
      Omsum=0.;
      DLOOPA(kk)
	Omsum+=Om[ii][kk]*Om[kk][jj];
	
      L[ii][jj]=delta(ii,jj)+1./(1.+gam)*Omsum+Om[ii][jj];
    }

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost lab <-> fluid frame
//whichdir: [LAB2FF, FF2LAB]
int boost22_laborff(int whichdir, FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom)
{ 
  int ii,jj,kk,ll;
  FTYPE Tt[NDIM][NDIM];

  //general Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  calc_Lorentz_laborff(whichdir,pp,q,ptrgeom,L);

  //copying
  DLOOP(ii,jj)
    Tt[ii][jj]=T1[ii][jj];
  
  //boosting
  DLOOP(ii,jj)
    {
      T2[ii][jj]=0.;
      DLOOP(kk,ll)
	T2[ii][jj]+=L[ii][kk]*L[jj][ll]*Tt[kk][ll];
    }

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost from lab to fluid frame
int
boost2_laborff(int whichdir, FTYPE A1[NDIM],FTYPE A2[NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom)
{ 
  int ii,jj,kk,ll;
  FTYPE At[NDIM];

  //general Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  calc_Lorentz_laborff(whichdir,pp,q,ptrgeom,L);

  //copying
  DLOOPA(ii)
    At[ii]=A1[ii];
  
  //boosting
  DLOOPA(ii)
    {
      A2[ii]=0.;
      DLOOPA(kk)
	A2[ii]+=L[ii][kk]*At[kk];
    }

  return 0; 
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//simple Lorentz boosts between ortonormal frames

int boost2_zamo2ff(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(ZAMO2FF, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

int boost2_ff2zamo(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(FF2ZAMO, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost fluid frame -> ZAMO
int boost2_fforzamo(int whichdir, FTYPE A1[NDIM],FTYPE A2[NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{
  // ptrgeom not used for now
  int i,j,k,l;
  FTYPE At[NDIM];
  FTYPE *ulab = q->ucon;
  FTYPE uzamo[NDIM];

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  //proper velocity for ZAMO
  FTYPE vpr[NDIM];
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  if(whichdir==ZAMO2FF){
    vpr[0]*=-1.0;
    vpr[1]*=-1.0;
    vpr[2]*=-1.0;
  }

  //Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  
  //Lorentz factor
  FTYPE vpr2=dot3(vpr,vpr); 
  FTYPE gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2; // NOTEMARK: maybe catestrophic cancellation issue
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }


  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }
  

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }


  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost ZAMO -> fluid frame
int boost22_zamo2ff(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(ZAMO2FF, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

int boost22_ff2zamo(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(FF2ZAMO, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost fluid frame -> ZAMO
int boost22_fforzamo(int whichdir, FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  // ptrgeom not used for now
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];
  FTYPE *ulab = q->ucon;
  FTYPE uzamo[NDIM];

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  //proper velocity for ZAMO
  FTYPE vpr[NDIM];
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  if(whichdir==ZAMO2FF){
    vpr[0]*=-1.0;
    vpr[1]*=-1.0;
    vpr[2]*=-1.0;
  }

  //Lorentz transformation matrix
  FTYPE L[NDIM][NDIM];
  
  //Lorentz factor
  FTYPE vpr2=dot3(vpr,vpr); 
  FTYPE gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2; // NOTEMARK: maybe catestrophic cancellation issue
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }


  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation ZAMO -> lab
int trans22_zamo2lab(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE elo[][NDIM])
{
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=elo[i][k]*elo[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation lab -> ZAMO
int trans22_lab2zamo(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE eup[][NDIM])
{
  int i,j,k,l;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=eup[k][i]*eup[l][j]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation lab -> ZAMO
int trans2_lab2zamo(FTYPE *u1,FTYPE *u2,FTYPE eup[][NDIM])
{
  int i,j,k;
  FTYPE ut[NDIM];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*eup[j][i];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation ZAMO -> lab
int trans2_zamo2lab(FTYPE *u1,FTYPE *u2,FTYPE elo[][NDIM])
{
  int i,j,k;
  FTYPE ut[NDIM];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*elo[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// T^ij -> T^i_j
int indices_2221(FTYPE T1[][NDIM],FTYPE T2[][NDIM], struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE Tt[NDIM][NDIM];

  for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<NDIM;k++)
	    {
	      Tt[i][j]+=T1[i][k]*ptrgeom->gcov[GIND(k,j)];
	    }	  
	}
    }

   for(i=0;i<NDIM;i++)
    {
      for(j=0;j<NDIM;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A^i -> A^_j
int indices_21(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
	{
	  At[i]+=A1[k]*ptrgeom->gcov[GIND(i,k)];
	}	  
    }

   for(i=0;i<NDIM;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A_i -> A^_j
int indices_12(FTYPE A1[NDIM],FTYPE A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  FTYPE At[NDIM];

  for(i=0;i<NDIM;i++)
    {
      At[i]=0.;
      for(k=0;k<NDIM;k++)
	{
	  At[i]+=A1[k]*ptrgeom->gcon[GIND(i,k)];
	}	  
    }

   for(i=0;i<NDIM;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}




#define TOZAMOFRAME 0 // reduce to ZAMO gammarel=1 frame (e.g. in non-GR that would be grid frame or v=0 frame or gammarel=1).
#define TOFLUIDFRAME 1 // reduce to using fluid frame (probably more reasonable in general).

#define M1REDUCE TOFLUIDFRAME // choose

//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//**********************************************************************
//**********************************************************************
//
///////////////
//
// Invert U->direct Primitive for radiation
// OLD (i.e. no longer true): (must come after HD or MHD or whatever sets velocity of fluid, because radiation needs to have updated velocity so that can define fluid frame)
// old code inside utoprimgen.c was:
//    struct of_state qrad;
// this uses new pr to get only ucon and ucov
//get_state_uconucovonly(pr, ptrgeom, &qrad); // OLD
// get new radiation primitives
//
// NEW (currently true): fluid frame no longer needed because go directly from lab-frame conserved quantities to lab-frame primitive quantities.
//
// ERADLIMIT applied here, but (SUPERGODMARK KORALTODO) should in general treat Erf<something as failure in inversion and try to average-out first before reducing to ERf=something and zero relative velocity.  And not sure what that something should be except as relative to other energy densities like how handle UU and BSQ relative to RHO.  Should do that instead and just have very small floor. 
//
// uu: Conserved quantities with URAD0,1,2,3 as radiation conserved quantities
// pp: primitives with PRAD0,1,2,3 as radiation primitive quantities
// ptrgeom: Standard pointer to geometry
// lpflag: see gobal.nondepmnemonics.h .  Tells u2p_rad() if can use/trust fluid velocity.
// lpflagrad: Should be set to indicate success of u2p_rad() inversion
//
// NOTES:
//
// Using *lpflag<=UTOPRIMNOFAIL to check for fluid inversion success rather than a SOFTer condition (e.g. no fail or IFUTOPRIMFAILSOFT==1) because only want to trust fluid as reduction of M1 in case where velocity is accurate with non-negative densities.
//
///////////////
int u2p_rad(FTYPE *uu, FTYPE *pp, struct of_geom *ptrgeom,PFTYPE *lpflag, PFTYPE *lpflagrad)
{
  
  if(WHICHVEL!=VELREL4){
    dualfprintf(fail_file,"u2p_rad() only setup for relative 4-velocity, currently.\n");
    myexit(137432636);
  }


  //////////////////////
  //
  // Prepare inversion from U->p for radiation assuming M1 closure
  //
  //////////////////////

  *lpflagrad=UTOPRIMRADNOFAIL;


  int verbose=0,i;
  FTYPE Rij[NDIM][NDIM];

  //conserved - R^t_mu
  FTYPE Av[NDIM]={uu[URAD0],uu[URAD1],uu[URAD2],uu[URAD3]};
  //indices up - R^tmu
  indices_12(Av,Av,ptrgeom);

  //g_munu R^tmu R^tnu
  FTYPE gRR=ptrgeom->gcov[GIND(0,0)]*Av[0]*Av[0]+ptrgeom->gcov[GIND(0,1)]*Av[0]*Av[1]+ptrgeom->gcov[GIND(0,2)]*Av[0]*Av[2]+ptrgeom->gcov[GIND(0,3)]*Av[0]*Av[3]+
    ptrgeom->gcov[GIND(1,0)]*Av[1]*Av[0]+ptrgeom->gcov[GIND(1,1)]*Av[1]*Av[1]+ptrgeom->gcov[GIND(1,2)]*Av[1]*Av[2]+ptrgeom->gcov[GIND(1,3)]*Av[1]*Av[3]+
    ptrgeom->gcov[GIND(2,0)]*Av[2]*Av[0]+ptrgeom->gcov[GIND(2,1)]*Av[2]*Av[1]+ptrgeom->gcov[GIND(2,2)]*Av[2]*Av[2]+ptrgeom->gcov[GIND(2,3)]*Av[2]*Av[3]+
    ptrgeom->gcov[GIND(3,0)]*Av[3]*Av[0]+ptrgeom->gcov[GIND(3,1)]*Av[3]*Av[1]+ptrgeom->gcov[GIND(3,2)]*Av[3]*Av[2]+ptrgeom->gcov[GIND(3,3)]*Av[3]*Av[3];
 
  //the quadratic equation for u^t of the radiation rest frame (urf[0])
  //supposed to provide two roots for (u^t)^2 of opposite signs
  FTYPE a,b,c,delta,gamma2;
  FTYPE urfcon[NDIM],urfconrel[NDIM],Erf;
  a=16.*gRR;
  b=8.*(gRR*ptrgeom->gcon[GIND(0,0)]+Av[0]*Av[0]);
  c=gRR*ptrgeom->gcon[GIND(0,0)]*ptrgeom->gcon[GIND(0,0)]-Av[0]*Av[0]*ptrgeom->gcon[GIND(0,0)];
  delta=b*b-4.*a*c;
  gamma2=  (-b-sqrt(delta))/2./a; // lab-frame gamma^2
  //if unphysical try the other root
  if(gamma2<0.) gamma2=  (-b+sqrt(delta))/2./a; 

  //cap on u^t
  FTYPE gammamax=GAMMAMAXRAD;

  FTYPE alpha=ptrgeom->alphalapse; //sqrtl(-1./ptrgeom->gcon[GIND(0,0)]);

  FTYPE gammarel2 = gamma2*alpha*alpha;  // /(-ptrgeom->gcon[GIND(TT,TT)]); // relative velocity gammarel^2
 




  //////////////////////
  //
  // Perform regular inversion
  // Or Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////
  int jj,kk;



  //////////////////////
  //
  // First case is if gammarel>gammamax, then set gammarel=gammamax unless Erf<ERADLIMIT (~0) in which case set Erf=ERADLIMIT and gammarel=1.
  // Note, can't set urfcon[0]=gammamax in case gammamax still remains space-like, e.g. inside horizon if gammamax isn't big enough.
  //
  //////////////////////
  if(gammarel2>gammamax*gammamax){

    //    urfcon[0]=gammamax; // ba
    FTYPE gammarel=gammamax;
    gammarel2=gammamax*gammamax;
      
    //proper radiation energy density in the radiation rest frame
    //    Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+ptrgeom->gcon[GIND(0,0)]);
    Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM


    if(Erf<ERADLIMIT){ // Erf too small
      *lpflagrad=UTOPRIMRADFAILCASE1A;
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;

      if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
      else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;

#if(PRODUCTION==0)
      dualfprintf(fail_file,"CASE1A: gammarel>gammamax and Erf<ERADLIMIT: gammarel2=%g gamma2=%g : i=%d j=%d k=%d\n",gammarel2,gamma2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif

    }
    else{ // Erf normal
      
      // check if just near gammamax, in which case don't need to worry about fixing
      //      if(fabs(gammarel2-gammamax*gammamax)>1E-12*gammamax*gammamax){
      //	*lpflagrad=UTOPRIMRADFAILCASE1B;
      //      }
      *lpflagrad=UTOPRIMRADFAILCASE1B;


      if(1){
	// lab-frame radiation relative 4-velocity
	FTYPE Aradrel[NDIM];
	SLOOPA(i) Aradrel[i] = alpha * (Av[i] + 1./3.*Erf*ptrgeom->gcon[GIND(0,i)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);

	// compute \gammarel using this
	FTYPE gammatemp,qsqtemp;
	int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);
	MYFUN(gamma_calc_fromuconrel(Aradrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);

	// now rescale Aradrel[i] so will give desired \gammamax
	SLOOPA(i) Aradrel[i] *= (gammamax/gammatemp);
	
	// copying to urfconrel
	SLOOPA(i) urfconrel[i]=Aradrel[i];

#if(PRODUCTION==0)
	// check that gamma really correctly gammamax
	FTYPE gammatemp2,qsqtemp2;
	MYFUN(gamma_calc_fromuconrel(Aradrel,ptrgeom,&gammatemp2,&qsqtemp2),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);
	dualfprintf(fail_file,"CASE1B: gammarel>gammamax and Erf normal: gammamax=%g gammatemp=%g gammatemp2=%g\n",gammamax,gammatemp,gammatemp2);
#endif


      }
      else if(0){ // Olek way
	
	// lab-frame radiation 4-velocity
	FTYPE Arad[NDIM];
	SLOOPA(i) Arad[i]=(Av[i]-1./3.*Erf*ptrgeom->gcon[GIND(0,i)])/(4./3.*Erf*gammamax);
      
	//is normalized now
	FTYPE Afac;
	a=0.; c=0.; b=0.;
	SLOOPA(i){
	  a+=Arad[i]*Arad[i]*ptrgeom->gcov[GIND(i,i)];
	  b+=2.*Arad[i]*ptrgeom->gcov[GIND(0,i)]*gammamax;
	}

	c=ptrgeom->gcov[GIND(0,0)]*gammamax*gammamax+1.0;
	delta=b*b-4.*a*c;
	Afac= (-b+sqrt(delta))/2./a;

	// lab-frame radiation 4-velocity
	urfcon[0]=gammamax;
	urfcon[1]=Afac*Arad[1];
	urfcon[2]=Afac*Arad[2];
	urfcon[3]=Afac*Arad[3];

	// relative 4-velocity radiation frame
	DLOOPA(i) urfconrel[i]=urfcon[i] - urfcon[0]*ptrgeom->gcon[GIND(0,i)]/ptrgeom->gcon[GIND(0,0)];
	
#if(PRODUCTION==0)
	
	// check that u.u=-1 (seems to be true)
	FTYPE urfcov[NDIM];
	DLOOPA(jj) urfcov[jj]=0.0;
	DLOOP(jj,kk) urfcov[jj] += urfcon[kk]*ptrgeom->gcov[GIND(jj,kk)];
	FTYPE udotu=0.0;
	DLOOPA(jj) udotu += urfcon[jj]*urfcov[jj];
	
	
	dualfprintf(fail_file,"CASE1B(Olek): Erf=%g Afac=%g Arad123=%g %g %g : udotu=%g : i=%d j=%d k=%d\n",Erf,Afac,Arad[1],Arad[2],Arad[3],udotu,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif
      }// end Olek method
	
    }// end else if Erf>0
  }
  //////////////////////
  //
  // Second case is if gammarel<1 or delta<0, then set gammarel=1.  If Erf<ERADLIMIT (~0), then set Erf=ERADLIMIT and gammarel=1.
  // Can't assume this condition is equivalent to large gamma, because if not, then leads to crazy boost of energy.
  //
  //////////////////////
  else if(gammarel2<1. || delta<0.){


    FTYPE gammarel2orig=gammarel2;
    // override
    gammarel2=1.0;
    FTYPE gammarel=1.0;  // use this below


    // get lab-frame 4-velocity u^t
    //    urfcon[0]=gammarel/ptrgeom->alphalapse;

    //radiative energy density in the radiation rest frame
    //    Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+ptrgeom->gcon[GIND(0,0)]);
    Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

    if(Erf<ERADLIMIT){ // JCM
      *lpflagrad=UTOPRIMRADFAILCASE2A;

      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;

      if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
      else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;

#if(PRODUCTION==0)
      dualfprintf(fail_file,"CASE2A: gamma<1 or delta<0 and Erf<ERADLIMIT : gammarel2=%g gamma2=%g : i=%d j=%d k=%d\n",gammarel2,gamma2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif

    }
    else{
      *lpflagrad=UTOPRIMRADFAILCASE2B;
      // relative 4-velocity radiation frame.  Might want to reset Erf to ERADLIMIT to be more strictly avoiding energy runaway problems.
      //      Erf=ERADLIMIT;

      if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
      else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;

#if(PRODUCTION==0)
      dualfprintf(fail_file,"CASE2B: gamma<1 or delta<0 and Erf normal : gammamax=%g gammarel2orig=%g gammarel2=%g gamma2=%g : i=%d j=%d k=%d\n",gammamax,gammarel2orig,gammarel2,gamma2,ptrgeom->i,ptrgeom->j,ptrgeom->k);
#endif

    }

  }
  //////////////////////
  //
  // Third case is if no bad conditions, then try regular calculation.  If Erf<ERADLIMIT, then reset Erf=ERADLIMIT and set gammarel=1
  //
  //////////////////////
  else{
    
    //radiative energy density in the radiation rest frame based upon lab-frame u^t
    //    urfcon[0]=sqrt(gamma2);
    //    Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+ptrgeom->gcon[GIND(0,0)]);
    Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

    if(Erf<ERADLIMIT){ // JCM
      *lpflagrad=UTOPRIMRADFAILCASE3A;
      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
      Erf=ERADLIMIT;

      if(M1REDUCE==TOFLUIDFRAME && *lpflag<=UTOPRIMNOFAIL) SLOOPA(jj) urfconrel[jj]=pp[U1+jj-1];
      else if(M1REDUCE==TOZAMOFRAME) SLOOPA(jj) urfconrel[jj]=0.0;

#if(PRODUCTION==0)
    dualfprintf(fail_file,"CASE3A: normal gamma, but Erf<ERADLIMIT.\n");
#endif

    }
    else{ // no failure yet for *lpflagrad=UTOPRIMRADFAILCASE3B;
      //relative velocity
      FTYPE gammarel=sqrt(gammarel2);
#if(1) // JCM
      SLOOPA(i) urfconrel[i] = alpha * (Av[i] + 1./3.*Erf*ptrgeom->gcon[GIND(0,i)]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
#else // Olek
      SLOOPA(i){
	urfconrel[i]=(3.*Av[i]-Erf*ptrgeom->gcon[GIND(0,i)])/(3.*Av[0]-Erf*ptrgeom->gcon[GIND(0,0)])/alpha+ptrgeom->gcon[GIND(0,i)]/alpha;
	urfconrel[i]*=gammarel;
      }
#endif

#if(PRODUCTION==0)
      //      dualfprintf(fail_file,"CASEnofail: normal gamma and normal Erf (default non-fixed) : Erf=%g : gamma=%g urfconrel123= %g %g %g\n",Erf,gammarel,urfcon[1],urfcon[2],urfcon[3]);
#endif

    }


  }

  
  //new primitives (only uses urfcon[1-3])
  pp[PRAD0]=Erf;
  pp[PRAD1]=urfconrel[1];
  pp[PRAD2]=urfconrel[2];
  pp[PRAD3]=urfconrel[3];

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(FTYPE *uu,FTYPE *pp, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM], FTYPE elo[][NDIM],FTYPE *f)
{
  FTYPE Rij[NDIM][NDIM];
  FTYPE ppp[NPR];

  calc_Rij_ff(pp,Rij);
  boost22_ff2zamo(Rij,Rij,pp,q,ptrgeom,eup);
  trans22_zamo2lab(Rij,Rij,elo);
  indices_2221(Rij,Rij,ptrgeom);

  FTYPE gdet=ptrgeom->gdet;

  // f = error
  f[0]=-Rij[0][0]+uu[URAD0];
  f[1]=-Rij[0][1]+uu[URAD1];
  f[2]=-Rij[0][2]+uu[URAD2];
  f[3]=-Rij[0][3]+uu[URAD3];

  return 0;
} 


// U->P inversion for Eddington approximation using Newton method.
int u2p_rad_num(FTYPE *uu, FTYPE *pp, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM], FTYPE elo[][NDIM])
{
  FTYPE pp0[NPR],pporg[NPR];
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE x[NDIM],f1[NDIM],f2[NDIM],f3[NDIM];
  int i,j,k,iter=0;



  for(i=PRAD0;i<=PRAD3;i++)
    {
      pporg[i]=pp[i];
    }
  
  do
    {
      iter++;
      for(i=PRAD0;i<NPR;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,q,ptrgeom,eup,elo,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+PRAD0]=pp[j+PRAD0]+PRADEPS*pp[PRAD0];
	    
	      f_u2prad_num(uu,pp,q,ptrgeom,eup,elo,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[PRAD0]);

	      pp[j+PRAD0]=pp0[j+PRAD0];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+PRAD0];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}

      for(i=0;i<4;i++)
	{
	  pp[i+PRAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+PRAD0]-pp0[i+PRAD0]);
	  f3[i]=fabs(f3[i]/pp0[PRAD0]);
	}

      if(f3[0]<RADCONV && f3[1]<RADCONV && f3[2]<RADCONV && f3[3]<RADCONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in u2prad_num()\n");
	  
	  for(i=PRAD0;i<NPR;i++)
	    {
	      pp[i]=pporg[i];
	    }
	  
	  return -1;

	  break; // WHY HERE IF UNREACHABLE?
	}
     
    }
  while(1);
  
  if(pp[PRAD0]<EFLOOR) 
    {
      printf("enegative u2prad()\n");
      pp[PRAD0]=EFLOOR;
    }
  

  return 0;
}

 



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives fluid frame -> ZAMO
// Use: Maybe dumping
// pp1 : Full set of primitives with radiation primitives replaced by fluid frame radiation \hat{E} and \hat{F}
// pp2 : Full set of primitives with ZAMO frame for radiation conserved quantities
int prad_ff2zamo(FTYPE *pp1, FTYPE *pp2, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  FTYPE Rij[NDIM][NDIM];
  int i,j;

  // set all (only non-rad needed) primitives
  int pliter,pl;
  PLOOP(pliter,pl) pp2[pl]=pp1[pl];

  calc_Rij_ff(pp1,Rij);
  boost22_ff2zamo(Rij,Rij,pp1,q,ptrgeom,eup);

  // overwrite pp2 with new radiation primitives
  pp2[PRAD0]=Rij[0][0];
  pp2[PRAD1]=Rij[0][1];
  pp2[PRAD2]=Rij[0][2];
  pp2[PRAD3]=Rij[0][3];

  return 0;
} 

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives ZAMO -> fluid frame - numerical solver
// Use: Error (f) for iteration
// ppff : HD primitives (i.e. WHICHVEL velocity) and radiation primitives of \hat{E} and \hat{F} (i.e. fluid frame conserved quantities)
// ppzamo : \hat{E} & \hat{F} -> E & F in ZAMO
int f_prad_zamo2ff(FTYPE *ppff, FTYPE *ppzamo, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM],FTYPE *f)
{
  FTYPE Rij[NDIM][NDIM];
  calc_Rij_ff(ppff,Rij);
  boost22_ff2zamo(Rij,Rij,ppff,q,ptrgeom,eup);

  f[0]=-Rij[0][0]+ppzamo[PRAD0];
  f[1]=-Rij[0][1]+ppzamo[PRAD1];
  f[2]=-Rij[0][2]+ppzamo[PRAD2];
  f[3]=-Rij[0][3]+ppzamo[PRAD3];


  return 0;
} 



// SUPERGODMARK: Is q constant ok?
// numerical iteration to find ff from zamo
// Use: Initial conditions
// ppzamo : E & F in ZAMO -> \hat{E} & \hat{F}
int prad_zamo2ff(FTYPE *ppzamo, FTYPE *ppff, struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{
  FTYPE pp0[NPR],pp[NPR];
  FTYPE J[NDIM][NDIM],iJ[NDIM][NDIM];
  FTYPE x[NDIM],f1[NDIM],f2[NDIM],f3[NDIM];
  int i,j,k,iter=0;


  
  //initial guess
  for(i=0;i<NPR;i++)
    {
      pp[i]=ppzamo[i];
    }


  //debug
  f_prad_zamo2ff(ppzamo,pp,q,ptrgeom,eup,f1);
  for(i=0;i<4;i++)
    {
      x[i]=pp[i+PRAD0];
    }  
 
  do
    {
      iter++;
      for(i=PRAD0;i<NPR;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_prad_zamo2ff(pp,ppzamo,q,ptrgeom,eup,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+PRAD0]=pp[j+PRAD0]+PRADEPS*pp[PRAD0];
	    
	      f_prad_zamo2ff(pp,ppzamo,q,ptrgeom,eup,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[PRAD0]);

	      pp[j+PRAD0]=pp0[j+PRAD0];
	    }
	}
  
      //inversion
      inverse_44matrix(J,iJ);

      //updating unknowns
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+PRAD0];
	}      

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}

      for(i=0;i<4;i++)
	{
	  pp[i+PRAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+PRAD0]-pp0[i+PRAD0]);
	  f3[i]=fabs(f3[i]/pp0[PRAD0]);
	}

      if(f3[0]<PRADCONV && f3[1]<PRADCONV && f3[2]<PRADCONV && f3[3]<PRADCONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in prad_zamo2ff()\n");
	  break;
	}
    }
  while(1);



  //returning prad
  for(i=0;i<NPR;i++)
    {
      ppzamo[i]=pp[i];
    }

  return 0;
}











//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
int calc_tautot(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *dx, FTYPE *tautot)
{
#if(0)
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  extern FTYPE calc_kappaes_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  FTYPE rho=pp[0];
  FTYPE u=pp[1];  
  FTYPE pr=(gamideal-1.)*(u);
  FTYPE T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  FTYPE kappa=calc_kappa_user(rho,T,xx[1],xx[2],xx[3]);
  FTYPE chi=kappa+calc_kappaes_user(rho,T,xx[1],xx[2],xx[3]);
#endif
  //xx[0] holds time
  FTYPE kappa,kappaes,chi;
  calc_kappa(pp,ptrgeom,&kappa);
  calc_kappaes(pp,ptrgeom,&kappaes);
  chi=kappa+kappaes;

  tautot[0]=chi*dx[0];
  tautot[1]=chi*dx[1];
  tautot[2]=chi*dx[2];

  return 0;
}

//**********************************************************************
//******* calculates abs opacity over dx[] ***************************
//**********************************************************************
int
calc_tauabs(FTYPE *pp, struct of_geom *ptrgeom, FTYPE *dx, FTYPE *tauabs)
{
#if(0)
  extern FTYPE calc_kappa_user(FTYPE rho, FTYPE T,FTYPE x,FTYPE y,FTYPE z);
  FTYPE rho=pp[0];
  FTYPE u=pp[1];  
  FTYPE pr=(gamideal-1.)*(u);
  FTYPE T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  FTYPE kappa=calc_kappa_user(rho,T,xx[1],xx[2],xx[3]);
#endif
  FTYPE kappa;
  calc_kappa(pp,ptrgeom,&kappa);

  tauabs[0]=kappa*dx[0];
  tauabs[1]=kappa*dx[1];
  tauabs[2]=kappa*dx[2];

  return 0;
}






//**********************************************************************
//suplementary routines for conversions
// SUPERGODMARK: KORALTODO: These should all be replaced or wrapped with HARM's eos.c functions, but they only appear to be used for initial conditions.  So for now, assume ideal gas.
//**********************************************************************
FTYPE calc_PEQ_ufromTrho(FTYPE T,FTYPE rho)
{
#if(0)
  //  FTYPE p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;
  FTYPE p=rho*T; // /MU_GAS; // KORALTODO  HARMUNITS? SUPERGODMARK
  FTYPE u=p/(gamideal-1.);
#else
  // if use local function function instead of below directly,
  // then assume user doesn't care about position for EOS.
  FTYPE u=u_rho0_T_simple(0, 0, 0, CENT, rho, T);
#endif
  return u;
}

FTYPE calc_PEQ_Tfromurho(FTYPE u,FTYPE rho)
{
#if(0)
  FTYPE p=u*(gamideal-1.);
  //  FTYPE T=p/(K_BOLTZ*rho/MU_GAS/M_PROTON);
  //  FTYPE T=p/(rho/MU_GAS); // KORALTODO  HARMUNITS? SUPERGODMARK
  FTYPE T=p/(rho); // KORALTODO  HARMUNITS? SUPERGODMARK
#else
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
#endif
  return T;
}

FTYPE calc_LTE_EfromT(FTYPE T)
{
  //  return 4.*SIGMA_RAD*T*T*T*T;
  return 4.*T*T*T*T; // KORALTODO: HARMUNITS? SUPERGODMARK
}

FTYPE calc_LTE_TfromE(FTYPE E )
{
  //  return sqrt(sqrt((E/4./SIGMA_RAD)));
  return sqrt(sqrt((E/4.)));
}


FTYPE calc_LTE_Efromurho(FTYPE u,FTYPE rho)
{
#if(0)
  FTYPE p=(gamideal-1.)*(u);
  //  FTYPE T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  //  FTYPE T=p*MU_GAS/rho; // KORALTODO  HARMUNITS? SUPERGODMARK
  FTYPE T=p/rho; // KORALTODO  HARMUNITS? SUPERGODMARK
#else
  FTYPE T=compute_temp_simple(0, 0, 0, CENT, rho, u);
#endif

  return calc_LTE_EfromT(T);
}



