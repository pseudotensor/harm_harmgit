#include "decs.h"

// convert primitives to conserved radiation quantities
int p2u_rad(ldouble *pr, ldouble *Urad, struct of_geom *ptrgeom, struct of_state *q)
{
  FTYPE Uraddiag;
  
  mhd_calc_rad(pr, TT, ptrgeom, q, Urad);

  // SUPERGODMARK: ensure \detg applied when filling conserved quantity

  return(0);
}

// compute radiation stres-energy tensor
void mhd_calc_rad(FTYPE *pr, int dir, struct of_geom *ptrgeom, struct of_state *q, FTYPE *radstressdir)
{
  FTYPE radT[NDIM][NDIM];


  FTYPE (*eup)[4],(*elo)[4];
  eup = GLOBALMETMACP2A0(boostemu,CENT,LAB2ZAMO,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  elo = GLOBALMETMACP2A0(boostemu,CENT,ZAMO2LAB,ptrgeom->i,ptrgeom->j,ptrgeom->k);

  FTYPE Rij[NDIM][NDIM];
  calc_Rij(pr,Rij);

  boost22_ff2zamo(Rij,Rij,pr,q,ptrgeom,eup);
  trans22_zamo2lab(Rij,Rij,elo);
  indices_2221(Rij,radT,ptrgeom);

  int jj;
  DLOOPA(jj) radstressdir[jj] = radT[dir][jj];


}

//**********************************************************************
//******* takes E and F^i from primitives and calculates radiation stress ****
//******* tensor R^ij using M1 closure scheme *****************************
//**********************************************************************
int calc_Rij(ldouble *pp, ldouble Rij[][4])
{
  ldouble E=pp[RAD0];
  ldouble F[3]={pp[RAD1],pp[RAD2],pp[RAD3]};

  ldouble nx,ny,nz,nlen,f;

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



ldouble my_min(ldouble a, ldouble b)
{
  if(a<b) return a;
  else return b;
}

ldouble my_sign(ldouble x)
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
int inverse_44matrix(ldouble a[][4], ldouble ia[][4])
{
  ldouble mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  ldouble	tmp[12]; ldouble	src[16]; ldouble det;
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



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives fluid frame -> ZAMO
int prad_ff2zamo(ldouble *pp1, ldouble *pp2, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4])
{
  ldouble Rij[4][4];
  int i,j;

  // set all (only non-rad needed) primitives
  int pliter,pl;
  PLOOP(pliter,pl) pp2[pl]=pp1[pl];

  calc_Rij(pp1,Rij);
  boost22_ff2zamo(Rij,Rij,pp1,q,ptrgeom,eup);

  // overwrite pp2 with new radiation primitives
  pp2[RAD0]=Rij[0][0];
  pp2[RAD1]=Rij[0][1];
  pp2[RAD2]=Rij[0][2];
  pp2[RAD3]=Rij[0][3];

  return 0;
} 

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//radiative primitives ZAMO -> fluid frame - numerical solver
int f_prad_zamo2ff(ldouble *ppff, ldouble *ppzamo, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4],ldouble *f)
{
  ldouble Rij[4][4];
  calc_Rij(ppff,Rij);
  boost22_ff2zamo(Rij,Rij,ppff,q,ptrgeom,eup);

  f[0]=-Rij[0][0]+ppzamo[RAD0];
  f[1]=-Rij[0][1]+ppzamo[RAD1];
  f[2]=-Rij[0][2]+ppzamo[RAD2];
  f[3]=-Rij[0][3]+ppzamo[RAD3];

  return 0;
} 



// SUPERGODMARK: Is q constant ok?
int prad_zamo2ff(ldouble *ppzamo, ldouble *ppff, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4])
{
  ldouble pp0[NV],pp[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;


  
  //initial guess
  for(i=0;i<NV;i++)
    {
      pp[i]=ppzamo[i];
    }


  //debug
  f_prad_zamo2ff(ppzamo,pp,q,ptrgeom,eup,f1);
  for(i=0;i<4;i++)
    {
      x[i]=pp[i+RAD0];
    }  
 
  do
    {
      iter++;
      for(i=RAD0;i<NV;i++)
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
	      pp[j+RAD0]=pp[j+RAD0]+PRADEPS*pp[RAD0];
	    
	      f_prad_zamo2ff(pp,ppzamo,q,ptrgeom,eup,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[RAD0]);

	      pp[j+RAD0]=pp0[j+RAD0];
	    }
	}
  
      //inversion
      inverse_44matrix(J,iJ);

      //updating unknowns
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+RAD0];
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
	  pp[i+RAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+RAD0]-pp0[i+RAD0]);
	  f3[i]=fabs(f3[i]/pp0[RAD0]);
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
  for(i=0;i<NV;i++)
    {
      ppzamo[i]=pp[i];
    }

  return 0;
}



int boost2_zamo2ff(ldouble A1[],ldouble A2[],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4])
{ 
  boost2_fforzamo(ZAMO2FF, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

int boost2_ff2zamo(ldouble A1[],ldouble A2[],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom,ldouble eup[][4])
{ 
  boost2_fforzamo(FF2ZAMO, A1, A2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost fluid frame -> ZAMO
int boost2_fforzamo(int whichdir, ldouble A1[4],ldouble A2[4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom,ldouble eup[][4])
{
  // ptrgeom not used for now
  int i,j,k,l;
  ldouble At[4];
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
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];

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
int boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4])
{ 
  boost22_fforzamo(ZAMO2FF, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

int boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom,ldouble eup[][4])
{ 
  boost22_fforzamo(FF2ZAMO, T1, T2, pp,q, ptrgeom,  eup);

 return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost fluid frame -> ZAMO
int boost22_fforzamo(int whichdir, ldouble T1[][4],ldouble T2[][4],ldouble *pp,struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4])
{
  // ptrgeom not used for now
  int i,j,k,l;
  ldouble Tt[4][4];
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
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];

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
int trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble elo[][4])
{
  int i,j,k,l;
  ldouble Tt[4][4];

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
int trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble eup[][4])
{
  int i,j,k,l;
  ldouble Tt[4][4];

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
int trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble eup[][4])
{
  int i,j,k;
  ldouble ut[4];

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
int trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble elo[][4])
{
  int i,j,k;
  ldouble ut[4];

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
int indices_2221(ldouble T1[][NDIM],ldouble T2[][NDIM], struct of_geom *ptrgeom)
{
  int i,j,k;
  ldouble Tt[NDIM][NDIM];

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
int indices_21(ldouble A1[NDIM],ldouble A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  ldouble At[NDIM];

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
int indices_12(ldouble A1[NDIM],ldouble A2[NDIM],struct of_geom *ptrgeom)
{
  int i,j,k;
  ldouble At[NDIM];

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




//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of LNRF to transform lab <--> LNRF
int calc_LNRFes(struct of_geom *ptrgeom, ldouble emuup[][4], ldouble emulo[][4])
{
  ldouble e2nu,e2psi,e2mu1,e2mu2,omega;
  ldouble gtt,gtph,gphph,grr,gthth;
  int i,j;
  // recast as [4][4] matrix
  //  ldouble emuup[][4]=(FTYPE (*)[4])(&ptremuup[0]);
  //ldouble emulo[][4]=(FTYPE (*)[4])(&ptremulo[0]);

  // SUPERGODMARK: Only applies for Boyer-Lindquist coordinates

  gtt=ptrgeom->gcov[GIND(0,0)];
  gtph=ptrgeom->gcov[GIND(0,3)];
  gphph=ptrgeom->gcov[GIND(3,3)];
  grr=ptrgeom->gcov[GIND(1,1)];
  gthth=ptrgeom->gcov[GIND(2,2)];

  //Bardeen's 72 coefficients:
  e2nu=-gtt+gtph*gtph/gphph;
  e2psi=gphph;
  e2mu1=grr;
  e2mu2=gthth;
  omega=-gtph/gphph;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	emuup[i][j]=0.;
	emulo[i][j]=0.;
      }

  emuup[0][0]=sqrt(e2nu);
  emuup[1][1]=sqrt(e2mu1);
  emuup[2][2]=sqrt(e2mu2);
  emuup[0][3]=-omega*sqrt(e2psi);
  emuup[3][3]=sqrt(e2psi);

  emulo[3][0]=omega*1./sqrt(e2nu);
  emulo[0][0]=1./sqrt(e2nu);
  emulo[1][1]=1./sqrt(e2mu1);
  emulo[2][2]=1./sqrt(e2mu2);
  emulo[3][3]=1./sqrt(e2psi);

  // need to return Kerr-Schild prime transformations

  return 0;
}





//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//**********************************************************************
//**********************************************************************
int u2p_rad(ldouble *uu, ldouble *pp, struct of_state *q, struct of_geom *ptrgeom)
{
  ldouble Rij[4][4];

  FTYPE (*eup)[4],(*elo)[4];
  eup = GLOBALMETMACP2A0(boostemu,CENT,LAB2ZAMO,ptrgeom->i,ptrgeom->j,ptrgeom->k);
  elo = GLOBALMETMACP2A0(boostemu,CENT,ZAMO2LAB,ptrgeom->i,ptrgeom->j,ptrgeom->k);


  //R^0_{mu}
  ldouble A[4]={uu[RAD0],uu[RAD1],uu[RAD2],uu[RAD3]};
  //indices up
  indices_12(A,A,ptrgeom);

  //covariant formulation
  
  //g_munu R^0mu R^0nu
  ldouble gRR=ptrgeom->gcov[GIND(0,0)]*A[0]*A[0]+ptrgeom->gcov[GIND(0,1)]*A[0]*A[1]+ptrgeom->gcov[GIND(0,2)]*A[0]*A[2]+ptrgeom->gcov[GIND(0,3)]*A[0]*A[3]+
    ptrgeom->gcov[GIND(1,0)]*A[1]*A[0]+ptrgeom->gcov[GIND(1,1)]*A[1]*A[1]+ptrgeom->gcov[GIND(1,2)]*A[1]*A[2]+ptrgeom->gcov[GIND(1,3)]*A[1]*A[3]+
    ptrgeom->gcov[GIND(2,0)]*A[2]*A[0]+ptrgeom->gcov[GIND(2,1)]*A[2]*A[1]+ptrgeom->gcov[GIND(2,2)]*A[2]*A[2]+ptrgeom->gcov[GIND(2,3)]*A[2]*A[3]+
    ptrgeom->gcov[GIND(3,0)]*A[3]*A[0]+ptrgeom->gcov[GIND(3,1)]*A[3]*A[1]+ptrgeom->gcov[GIND(3,2)]*A[3]*A[2]+ptrgeom->gcov[GIND(3,3)]*A[3]*A[3];
 
  //the quadratic equation for u^t of the radiation rest frame (urf[0])
  ldouble a,b,c;
  a=16.*gRR;
  b=8.*(gRR*ptrgeom->gcon[GIND(0,0)]+A[0]*A[0]);
  c=gRR*ptrgeom->gcon[GIND(0,0)]*ptrgeom->gcon[GIND(0,0)]-A[0]*A[0]*ptrgeom->gcon[GIND(0,0)];
  ldouble delta=b*b-4.*a*c;
  ldouble urf[4],Erf;
  urf[0]=sqrt((-b-sqrt(delta))/2./a);
  if(isnan(urf[0])) urf[0]=1.;

  //radiative energy density in the radiation rest frame
  Erf=3.*A[0]/(4.*urf[0]*urf[0]+ptrgeom->gcon[GIND(0,0)]);

  //four-velocity of the rest frame
  urf[1]=3./(4.*Erf*urf[0])*(A[1]-1./3.*Erf*ptrgeom->gcon[GIND(0,1)]);
  urf[2]=3./(4.*Erf*urf[0])*(A[2]-1./3.*Erf*ptrgeom->gcon[GIND(0,2)]);
  urf[3]=3./(4.*Erf*urf[0])*(A[3]-1./3.*Erf*ptrgeom->gcon[GIND(0,3)]);

  //lab frame:
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=4./3.*Erf*urf[i]*urf[j]+1./3.*Erf*ptrgeom->gcon[GIND(i,j)];

  //boosting to ff
  trans22_lab2zamo(Rij,Rij,eup);
  boost22_zamo2ff(Rij,Rij,pp,q,ptrgeom,eup);

  //reading primitives
  pp[RAD0]=Rij[0][0];
  pp[RAD1]=Rij[0][1];
  pp[RAD2]=Rij[0][2];
  pp[RAD3]=Rij[0][3];

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(ldouble *uu,ldouble *pp, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4], ldouble elo[][4],ldouble *f)
{
  ldouble Rij[4][4];
  ldouble ppp[NV];

  calc_Rij(pp,Rij);
  boost22_ff2zamo(Rij,Rij,pp,q,ptrgeom,eup);
  trans22_zamo2lab(Rij,Rij,elo);
  indices_2221(Rij,Rij,ptrgeom);

  ldouble gdet=ptrgeom->gdet;

  f[0]=-Rij[0][0]+uu[RAD0];
  f[1]=-Rij[0][1]+uu[RAD1];
  f[2]=-Rij[0][2]+uu[RAD2];
  f[3]=-Rij[0][3]+uu[RAD3];

  return 0;
} 


int u2p_rad_num(ldouble *uu, ldouble *pp, struct of_state *q, struct of_geom *ptrgeom, ldouble eup[][4], ldouble elo[][4])
{
  ldouble pp0[NV],pporg[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;



  for(i=RAD0;i<NV;i++)
    {
      pporg[i]=pp[i];
    }
  
  do
    {
      iter++;
      for(i=RAD0;i<NV;i++)
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
	      pp[j+RAD0]=pp[j+RAD0]+PRADEPS*pp[RAD0];
	    
	      f_u2prad_num(uu,pp,q,ptrgeom,eup,elo,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(PRADEPS*pp[RAD0]);

	      pp[j+RAD0]=pp0[j+RAD0];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+RAD0];
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
	  pp[i+RAD0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+RAD0]-pp0[i+RAD0]);
	  f3[i]=fabs(f3[i]/pp0[RAD0]);
	}

      if(f3[0]<RADCONV && f3[1]<RADCONV && f3[2]<RADCONV && f3[3]<RADCONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in u2prad_num()\n");
	  
	  for(i=RAD0;i<NV;i++)
	    {
	      pp[i]=pporg[i];
	    }
	  
	  return -1;

	  break; // WHY HERE IF UNREACHABLE?
	}
     
    }
  while(1);
  
  if(pp[RAD0]<EFLOOR) 
    {
      printf("enegative u2prad()\n");
      pp[RAD0]=EFLOOR;
    }
  

  return 0;
}

 
