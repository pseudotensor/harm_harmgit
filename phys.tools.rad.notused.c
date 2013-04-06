
// UNUSED NUMERICAL STUFF
#define RADEPS (1.e-6) // for unused numerical inversion stuff
#define RADCONV (1.e-7) // for unused numerical inversion stuff
#define PRADEPS (1.e-6)  // for unused numerical inversion stuff
#define PRADCONV (1.e-8)  // for unused numerical inversion stuff

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//calculates general Lorenz matrix for lab <-> ff
//whichdir: [LAB2FF, FF2LAB]
// SUPERGODMARK: What is lab frame velocity?
// UNUSED -- KORALTODO :REMOVE?
int calc_Lorentz_laborff(int whichdir,FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE L[][NDIM])
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
int boost2_laborff(int whichdir, FTYPE A1[NDIM],FTYPE A2[NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom)
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
// UNUSED -- KORALTODO :REMOVE?
int boost2_zamo2ff(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(ZAMO2FF, A1, A2, pp,q, ptrgeom,  eup);

  return 0;
}

// UNUSED -- KORALTODO :REMOVE?
int boost2_ff2zamo(FTYPE A1[],FTYPE A2[],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost2_fforzamo(FF2ZAMO, A1, A2, pp,q, ptrgeom,  eup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost fluid frame -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
int boost22_zamo2ff(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom, FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(ZAMO2FF, T1, T2, pp,q, ptrgeom,  eup);

  return 0;
}

// UNUSED -- KORALTODO :REMOVE?
int boost22_ff2zamo(FTYPE T1[][NDIM],FTYPE T2[][NDIM],FTYPE *pp,struct of_state *q, struct of_geom *ptrgeom,FTYPE eup[][NDIM])
{ 
  boost22_fforzamo(FF2ZAMO, T1, T2, pp,q, ptrgeom,  eup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost fluid frame -> ZAMO
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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


//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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
  
  if(pp[PRAD0]<ERADLIMIT) 
    {
      printf("%g=PRAD0<ERADLIMIT=%g u2prad()\n",pp[PRAD0],ERADLIMIT);
      pp[PRAD0]=ERADLIMIT;
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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
// UNUSED -- KORALTODO :REMOVE?
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



