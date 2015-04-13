
#include "decs.h"

// issues:

// 1) Sometimes returns v>c 3-velocity even with P limiter.
//    e.g. aligned pulsar problem with no equatorial AVOIDCS crashes


/*

  Generally:

  \dF^{\alpha\beta} = (1/2) e^{\alpha\beta\mu\nu} F_{\mu\nu}

  F^{\alpha\beta} = (-1/2) e^{\alpha\beta\mu\nu} \dF_{\mu\nu}

  e^{\alpha\beta\gamma\delta} = (-1/\detg) [\alpha\beta\gamma\delta]

  e_{\alpha\beta\gamma\delta} = (\detg) [\alpha\beta\gamma\delta]

  [ijkl] = Signature[{i,j,k,l}]

  With

  F^{\alpha\beta} = A n^\alpha E^\beta + B n^\beta E^\alpha + C B_\gamma n_\delta e^{\alpha\beta\gamma\delta}

  and

  \dF^{\alpha\beta} = D n^\alpha B^\beta + F n^\beta B^\alpha + G E_\gamma n_\delta e^{\alpha\beta\gamma\delta}

  then antisymmetry gives A=-B and D=-F

  where the definition of the dual gives:

  C=-F=D
  A=-G

  which then gives 4 independent equations for 6 unknowns.

  Let A=1=-C=-1

  then:

  F^{\alpha\beta} =  n^\alpha E^\beta - n^\beta E^\alpha - B_\gamma n_\delta e^{\alpha\beta\gamma\delta}

  and

  \dF^{\alpha\beta} = - n^\alpha B^\beta + n^\beta B^\alpha - E_\gamma n_\delta e^{\alpha\beta\gamma\delta}

  As in Ian FFDE paper.

  // Gammie:

  \dF^{\alpha\beta} = b^\alpha u^\beta - b^\beta u^\alpha = -u^\alpha b^\beta + u^\beta b^\alpha

  If \eta->u and then A=1.  Let C=-1 also if setting coefficient of e^\alpha terms.







  //Jon notation:

  // see grmhd-faraday_allforms.nb

  F^{\mu\nu} = -1/\detg {{0,-Eu1,-Eu2,-Eu3},{Eu1,0,-Bd3,Bd2},{Eu2,Bd3,0,-Bd1},{Eu3,-Bd2,Bd1,0}}

  \dF_{\mu\nu} = {{0,Bd1,Bd2,Bd3},{-Bd1,0,Eu3,-Eu2},{-Bd2,-Eu3,0,Eu1},{-Bd3,Eu2,-Eu1,0}}


  F^{ij} = (1/\detg) [ijk] B_k 

  F_{ij}/\detg = B^k [ijk]

  [ijk] F_{ij}/\detg = 2B^k


  \dF^{ij} = [ijk] E_k and \dF^{ij}[ijk] = 2 E_k


  \dF_{ij} = [ijk] E^k


  E_i = F_{it}/\detg

  E^i = F^{ti}\detg // changed sign at some point (now IanPrimitive^i=F^{it} = -E^i/\detg)

  B^i = \dF^{it}  (negative of Ian field Primitive -- IanPrim^i=\dF^{it} = -B^i)

  B_i = \dF_{ti} // changed sign at some point

  This makes E.B=0


  b^\mu = P^\mu_\nu B^\nu/u^t

  b_\mu = P_\mu^\nu B_\nu/u_t






*/


// whether to check if B^2-E^2>0 is satisfied
#define CHECKBSQMINUSESQ 1

// whether to account for dissipation of field
#define ACCOUNTDISS 1

// whether to account for lack of energy conservation in accounting
#define ACCOUNTNONCONS 1

// whether to avoid current sheet dissipation by setting the cross-field drift to vanish.
// currently very model dependent
// need a contact discontinuity detector
#define AVOIDCS 0

// whether to fixup B3 to satisfy u.B=0 even when AVOIDCS=1
#define FIXB3 0 // doesn't work -- causes u^t problems


// 0: komissarov
// 1: jon#1
// 2: jon#2 (paper version BEST)
// 3: jon#3 (trial Lorentz factor CS limit)
#define AVOIDUTINF 2


// whether to remove a T^t_i component and replace it with T^t_t
// due to the number of operations, this apparently reduces accuracy of inversion by a couple orders of magnitude to 10^(-13) instead of 10^(-15).
// actually, can even have order unity relative errors in small velocities if one velocity dominates

// even worse: method completely doesn't work.  Apparently evolution of T^t_t is strongly divergent from evolution of T^t_i
// don't use till more testing in 1D and point-average transformation tests to reduce truncation error.
#define REMOVET 0

// assumes geometry-free U in standard form
// latest version
int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  void TtoEconOalpha(FTYPE *T0, FTYPE *BcovOalpha, FTYPE BsqOalphasq, struct of_geom *geom, FTYPE *EconOalpha);
  int TtoT(FTYPE BsqOalphasq, FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *Tout);
  int i,j,k,pl,pliter;
  FTYPE alphasq;
  FTYPE Fit[NDIM];
  FTYPE EconOalpha[NDIM],EcovOalpha[NDIM],BcovOalpha[NDIM],BconOalpha[NDIM],realBsqOalphasq;
  FTYPE vcon[NDIM],ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE BsqOalphasq,EsqOalphasq,BsqmEsqOalphasq;
  FTYPE gammasq,Gsq;
  int modified;
  int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr);
  void EBcovtovcon(FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE Bsq, struct of_geom *geom, FTYPE*vcon);
  FTYPE prim[NPR];
  int ti,tj,tk;
  FTYPE ucov[NDIM];
  FTYPE Tin[NDIM],Tout[NDIM];
  void U2BBsq(FTYPE *U, struct of_geom *geom, FTYPE *BconOalpha, FTYPE *BcovOalpha, FTYPE *BsqOalphasq);



  ////////////////////////
  //
  // get lapse
  //
  ///////////////////////

  // define \eta_\alpha
  // assume always has 0 value for space components
  alphasq = 1./(- geom->gcon[GIND(0,0)]);  
  
  /////////////
  //
  // get (B^\mu/\alpha)
  //
  /////////////

  U2BBsq(U, geom, BconOalpha, BcovOalpha, &BsqOalphasq);

  ///////////////
  //
  // get T^t_\mu
  //
  /////////////
  DLOOPA(j) Tin[j] = Tout[j] = U[UU+j] ; /* T^t_\mu terms (t up, \mu down) */

#if(REMOVET)
  TtoT(BsqOalphasq,alphasq, BconOalpha, Tin, geom, Tout);
#endif

  /////////////////
  //
  // get (E^\mu/\alpha)
  //
  /////////////////
  TtoEconOalpha(Tout, BcovOalpha, BsqOalphasq, geom, EconOalpha);


  /////////////////////////
  //
  // modify velocity so B^2-E^2>0 always
  //
  ////////////////////////

  // set no modification to energy
  modified=0;

  realBsqOalphasq=BsqOalphasq; // default
#if(CHECKBSQMINUSESQ)


  lower_vec(EconOalpha,geom,EcovOalpha);

  EsqOalphasq = 0.0;
  SLOOPA(j) EsqOalphasq += EconOalpha[j]*EcovOalpha[j];

  BsqmEsqOalphasq = BsqOalphasq-EsqOalphasq;

  
  //check if effective dissipation is occuring

#if(AVOIDUTINF==0)
  if(BsqmEsqOalphasq<0.0){
    //    dualfprintf(fail_file,"You are going to fail! %21.15g\n",BsqmEsqOalphasq);
    // then modify drift velocity
    // see Komissarov (2005)
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq));
    modified=1;
  }
#elif(AVOIDUTINF==1)
  gammasq=BsqOalphasq/BsqmEsqOalphasq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsqOalphasq));
    Gsq=1.0/(GAMMAMAX*GAMMAMAX);
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq)*sqrt(1.0+Gsq));
    // else assume fine
    modified=1;
  }
#elif(AVOIDUTINF==2)
  gammasq=BsqOalphasq/BsqmEsqOalphasq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    Gsq=1.0-1.0/(GAMMAMAX*GAMMAMAX);
    realBsqOalphasq=(sqrt(EsqOalphasq*BsqOalphasq/Gsq));
    modified=1;
  }
#elif(AVOIDUTINF==3)
  // limit Lorentz factor to 1 in current sheets
  // doesn't make sense to limit completely like this since directional dependent
  if((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsqOalphasq));
    Gsq=1.0/(1.0);
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq)*sqrt(1.0+Gsq));
    // else assume fine
    modified=1;
  }
#endif

#endif // end if checking


  /////////////////////
  //
  // Find v^i
  //
  /////////////////////

  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  
  // actually, we removed 1/alpha^2 from Bsq thing, and removed alpha from Bcov and Ecov each (total alpha^2), so left with alpha^2
  // so alpha cancels between Ecov*Bcov/Bsq, so removed versions and versions with it in use same function
  EBcovtovcon(alphasq, EcovOalpha, BcovOalpha, realBsqOalphasq, geom, vcon);



  ///////////////
  //
  // avoid current sheets
  //
  /////////////
  // tried various ways to avoid dissipation in current sheet (CS).
  // 1) use gamma~1 in CS : kinda worked
  // 2) set B2=0 in CS: didn't work
  // 3) set vcon2=0 in CS: works very well!


  // symmetric around equator for centered zones
  ti=startpos[1]+geom->i;
  tj=startpos[2]+geom->j;
  tk=startpos[3]+geom->k;
  if(AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    vcon[TH]=0.0;
    modified=1;
  }

  //////////////////////////
  //
  // convert from v^i -> u^\mu -> primitive velocity
  //
  //////////////////

  //  PLOOP(pliter,pl) prim[pl]=0.0;
  //  prim[U1]=vcon[1];
  //  prim[U2]=vcon[2];
  //  prim[U3]=vcon[3];

  if(vcon2pr(WHICHVEL,vcon,geom,pr)>=1){
    dualfprintf(fail_file, "Utoprim_ffde(vcon2pr): space-like error\n");
    return(1);
  }

  //  ucon_calc_3vel(prim,geom,ucon, others);
  //  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  //  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  //////////////////////////
  //
  // now assign field primitives
  //
  ///////////////////////
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];

  //////////////////
  //
  // rearrange B^\phi so that u.B=0 even after current sheet changes
  //
  /////////////////
  if(FIXB3&&AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    // doesn't work
    PLOOP(pliter,pl) prim[pl]=0.0;
    prim[U1]=vcon[1];
    prim[U2]=vcon[2];
    prim[U3]=vcon[3];
    ucon_calc_3vel(prim,geom,ucon,others);
    lower_vec(ucon,geom,ucov);
    pr[B3]=-ucov[RR]/ucov[PH]*pr[B1];
    // this only works in axisymmetry (else divb!=0)
  }
    


  //////////////////
  //
  // Account for any modifications
  //
  /////////////////

#if(ACCOUNTDISS)
  // account for any energy change
  if(1||modified){ // 1|| because if not conserving energy, then really changed

    // used by diagnostics
    //SLOOPA(j) Econ[j]=-Fit[j];
    Fit[TT]=0.0;
    SLOOPA(j) Fit[j]=-EconOalpha[j];

    diss_account(modified,Fit,U,geom,pr);
  }
#endif
  // done.

  //  dualfprintf(fail_file,"Bsq %21.15g %21.15g\n",BsqOalphasq,realBsqOalphasq);
  //  dualfprintf(fail_file,"vcon: %21.15g %21.15g %21.15g\n",vcon[1],vcon[2],vcon[3]);
  //  dualfprintf(fail_file,"pr: %21.15g %21.15g %21.15g\n",pr[U1],pr[U2],pr[U3]);


  return(0);


}




// assumes geometry-free U in standard form
int Utoprim_ffde_oldnew(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  int i,j,k,pl,pliter;
  int ti,tj,tk;
  void UtoFitMxBsqgeom(FTYPE *U, struct of_geom *geom, FTYPE *Fit, FTYPE*Mx, FTYPE*oBsqgeom);
  FTYPE alphasq,beta[NDIM];
  FTYPE Fit[NDIM],Mx[NDIM],Mt[NDIM],oBsqgeom;
  FTYPE Econ[NDIM],Ecov[NDIM],Bcov[NDIM],Bcon[NDIM],realoBsqgeom;
  FTYPE vcon[NDIM],ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE Bsq,Esq,BsqmEsq;
  FTYPE gammasq,Gsq;
  int modified;
  int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr);
  void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE realoBsqgeom, FTYPE*vcon);
  //  FTYPE prim[NPR];


  // set no modification to energy
  modified=0;

  // lapse
  // define \eta_\alpha
  // assume always has 0 value for space components
  alphasq = 1./(- geom->gcon[GIND(0,0)]);  
  
  SLOOPA(j) beta[j] = geom->gcon[GIND(0,j)]*alphasq ;

  UtoFitMxBsqgeom(U,geom,Fit,Mx,&oBsqgeom);
  // -> B_\mu = -alpha * M^t_\mu
  // Real Bsq = Bsqgeom*alpha^2

  Econ[TT]=0.0;
  // remove alpha since cancels below
  //  SLOOPA(j) Econ[j]=-alpha*Fit[j];
  SLOOPA(j) Econ[j]=-Fit[j];
  lower_vec(Econ,geom,Ecov);

  // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
  // remove alpha since cancels below
  //SLOOPA(j) Bcov[j]=-alpha*Mx[j];
  SLOOPA(j) Bcov[j]=-Mx[j]; // Mx = \dF^t_j

  // remove alpha^2 since cancels below
  //  realoBsqgeom=oBsqgeom/alphasq;
  realoBsqgeom=oBsqgeom;

#if(CHECKBSQMINUSESQ)
  // check B^2-E^2
  Bcon[TT]=0.0;
  // Bcon without alpha
  SLOOPA(j) Bcon[j] = U[B1+j-1] ;  /* Bcon=-M^{ti} , where U[Bi]=B^i=M^{it}*/

  Esq = 0.0;
  SLOOPA(j) Esq += Econ[j]*Ecov[j];

  Bsq = 0.0;
  SLOOPA(j) Bsq += Bcon[j]*Bcov[j];

  BsqmEsq = Bsq-Esq;

  // tried various ways to avoid dissipation in current sheet (CS).
  // 1) use gamma~1 in CS : kinda worked
  // 2) set B2=0 in CS: didn't work
  // 3) set vcon2=0 in CS: works very well!

  
  //check if effective dissipation is occuring

#if(AVOIDUTINF==0)
  if(BsqmEsq<0.0){
    //    dualfprintf(fail_file,"You are going to fail! %21.15g\n",BsqmEsq);
    // then modify drift velocity
    // see Komissarov (2005)
    realoBsqgeom=1.0/(max(Bsq,Esq)*geom->gdet);
    modified=1;
  }
#elif(AVOIDUTINF==1)
  gammasq=Bsq/BsqmEsq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsq));
    Gsq=1.0/(GAMMAMAX*GAMMAMAX);
    realoBsqgeom=1.0/(max(Bsq,Esq)*sqrt(1.0+Gsq)*geom->gdet);
    // else assume fine
    modified=1;
  }
#elif(AVOIDUTINF==2)
  gammasq=Bsq/BsqmEsq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    Gsq=1.0-1.0/(GAMMAMAX*GAMMAMAX);
    realoBsqgeom=1.0/(sqrt(Esq*Bsq/Gsq)*geom->gdet);
    modified=1;
  }
#elif(AVOIDUTINF==3)
  // limit Lorentz factor to 1 in current sheets
  // doesn't make sense to limit completely like this since directional dependent
  ti=startpos[1]+geom->i;
  tj=startpos[2]+geom->j;
  tk=startpos[3]+geom->k;

  if((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsq));
    Gsq=1.0/(1.0);
    realoBsqgeom=1.0/(max(Bsq,Esq)*sqrt(1.0+Gsq)*geom->gdet);
    // else assume fine
    modified=1;
  }
#endif

#endif // end if checking

  /////////////////////
  //
  // Find v^i
  //
  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  
  // actually, we removed 1/alpha^2 from Bsq thing, and removed alpha from Bcov and Ecov each (total alpha^2), so left with alpha^2
  //  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*realoBsqgeom;
  //  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*realoBsqgeom;
  //  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*realoBsqgeom;
  //  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*realoBsqgeom;
  //  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*realoBsqgeom;
  //  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*realoBsqgeom;


  EBcovtovcon_old(beta,alphasq,Ecov,Bcov,realoBsqgeom,vcon);

  // symmetric around equator for centered zones
  ti=startpos[1]+geom->i;
  tj=startpos[2]+geom->j;
  tk=startpos[3]+geom->k;
  if(AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    vcon[2]=0.0;
    modified=1;
  }
  //////////////////////////
  //
  // convert from v^i -> u^\mu
  //  PLOOP(pliter,pl) prim[pl]=0.0;
  //  prim[U1]=vcon[1];
  //  prim[U2]=vcon[2];
  //  prim[U3]=vcon[3];

  if(vcon2pr(WHICHVEL,vcon,geom,pr)>=1){
    dualfprintf(fail_file, "Utoprim_ffde(vcon2pr): space-like error\n");
    return(1);
  }

  //  ucon_calc_3vel(prim,geom,ucon,others);
  //  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  //  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  // now assign field primitives
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];
  


#if(ACCOUNTDISS)
  // account for any energy change
  if(1||modified){ // 1|| because if not conserving energy, then really changed
    diss_account(modified, Fit,U,geom,pr);
  }
#endif
  // done.


  return(0);


}

void EBcovtovcon(FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE Bsq, struct of_geom *geom, FTYPE*vcon)
{
  FTYPE oBsqgeom;
  
  oBsqgeom=1.0/(Bsq*geom->gdet);
  
  vcon[1]=alphasq*(-geom->gcon[GIND(TT,1)] + (Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*oBsqgeom);
  vcon[2]=alphasq*(-geom->gcon[GIND(TT,2)] + (Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*oBsqgeom);
  vcon[3]=alphasq*(-geom->gcon[GIND(TT,3)] + (Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*oBsqgeom);
}

void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE oBsqgeom, FTYPE*vcon)
{
  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*oBsqgeom;
  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*oBsqgeom;
  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*oBsqgeom;
}



// account for dissipation or truncation level lack of conservation
int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr)
{
  int j;
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE (*T)[NDIM]);
  extern int primtoU(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *geom,FTYPE *U, FTYPE *Uabs);
  struct of_state q;
  FTYPE T[NDIM][NDIM];
  FTYPE Ui[NPR],Uf[NPR];


  // do this to get initial T consistent with momentum (excludes energy) equations
  FitUtoT(Fit,U,geom,T);

  // initial energy-momentum and field
  Ui[RHO]=0.0;
  DLOOPA(j) Ui[UU+j]=geom->gdet*T[TT][j];
  SLOOPA(j) Ui[B1+j-1]=geom->gdet*U[B1+j-1];
  
  // get final energy-momentum and field
  MYFUN(get_state(pr, geom, &q),"phys.ffde.c:Utoprim_ffde()", "get_state()", 1);
  primtoU(UDIAG,pr,&q,geom,Uf, NULL);
  Uf[RHO]=0.0;

#if(ACCOUNTNONCONS)
  Uf[UU]-=(geom->gdet*U[UU]-Ui[UU]); // add back in lost energy from nonconservation of energy method
#endif
  
  // account for change
  // diag_fixup_U is forced to account on every time substep since can't step through problems like with negative internal energy.
  if(modified) diag_fixup_U((DOENOFLUX != NOENOFLUX),Ui,Uf,U,geom,0,COUNTUTOPRIMFAILCONV); // real dissipation, absolute limit of Lorentz factor
  else diag_fixup_U((DOENOFLUX != NOENOFLUX),Ui,Uf,U,geom,0,COUNTFLOORACT); // truncation dissipation
   

  return(0);
}  



// input \Omega_F and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  int limit_3vel(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr);


  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  ftemp=(Bccov[TT]+omegaf*Bccov[PH])/Bsq;
  // ftemp2 is set so that no round off error in limit where toroidal field dominates
  // Does this cause problem for frame-dragged fields?
  //ftemp2=omegaf*(1.0 - Bccov[PH]*Bccon[3]/Bsq);
  // below more accurate than above?
  ftemp2=omegaf*(Bccov[RR]*Bccon[RR]+Bccov[TH]*Bccon[TH])/Bsq;

  vcon[1] = -Bccon[1]*ftemp;
  vcon[2] = -Bccon[2]*ftemp;
  //  vcon[3] =  omegaf  - Bccon[3]*ftemp;
  // designed so to avoid catastrophic cancellation like above has problems with
  vcon[3] =  ftemp2 - (Bccon[3]*Bccov[TT]/Bsq);


  // see if 3-velocity is good, otherwise limit to Lorentz factor of <=GAMMAMAX
  // gives final pr as well
  limit_3vel(Bccon, geom, vcon, pr);

  // now can convert to pr  
  //  if(vcon2pr(WHICHVEL,vcon,geom,pr)>=1){
  //    dualfprintf(fail_file, "OBtopr(vcon2pr): space-like error\n");
  //    return(1);
  //  }
  pr[RHO]=pr[UU]=0.0;
  //  pr[B1]=Bccon[B1];
  //  pr[B2]=Bccon[B2];
  //  pr[B3]=Bccon[B3];

  return(0);

}





// check to make sure 3-velocity is good, otherwise will have to limit it
// Bcon and vcon are code versions
// returns code primitives
int limit_3vel(FTYPE *Bcon, struct of_geom *geom, FTYPE *vcon, FTYPE *pr)
{
  int dualf_calc(FTYPE *Bcon, FTYPE *vcon, FTYPE (*dualffull)[NDIM]);
  void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM]);
  FTYPE dualf[NDIM][NDIM], T[NDIM][NDIM];
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  FTYPE U[NPR];

  // get \dF^{\mu\nu}
  dualf_calc(Bcon, vcon, dualf);

  // get T^\mu_\nu
  ffdestresstensor(dualf, geom, T);

  // setup conserved quantities
  U[RHO] = 0;
  U[UU] = T[0][0]; // T^t_t
  U[U1] = T[0][1]; // T^t_x
  U[U2] = T[0][2]; // T^t_y
  U[U3] = T[0][3]; // T^t_z
  U[B1] = Bcon[1];
  U[B2] = Bcon[2];
  U[B3] = Bcon[3];

  // filter through inversion so v^i is limited to have Lorentz factor <=GAMMAMAX
  Utoprim_ffde(U, geom, pr);


  return(0);

}




// input E and B are defined as in GRFFE paper
int EBtopr_2(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE alphasq,beta[NDIM],alpha;
  FTYPE Bsq,realoBsqgeom;
  FTYPE Ecov[NDIM],Bcov[NDIM];
  FTYPE vcon[NDIM];
  void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE realoBsqgeom, FTYPE*vcon);

  

  alphasq = 1./(- geom->gcon[GIND(0,0)]);  
  alpha=sqrt(alphasq);
  SLOOPA(j) beta[j] = geom->gcon[GIND(0,j)]*alphasq ;

  lower_vec(Econ,geom,Ecov);
  lower_vec(Bcon,geom,Bcov);

  Bsq = 0.0;
  SLOOPA(j) Bsq += Bcon[j]*Bcov[j];
  realoBsqgeom=1.0/(Bsq*geom->gdet);

  EBcovtovcon_old(beta,alphasq,Ecov,Bcov,realoBsqgeom,vcon);

  if(vcon2pr(WHICHVEL,vcon,geom,pr)>=1){
    dualfprintf(fail_file, "Utoprim_ffde(vcon2pr): space-like error\n");
    return(1);
  }

  pr[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  pr[B2]=-Bcon[2]/(-alpha);
  pr[B3]=-Bcon[3]/(-alpha);

  return(0);

}


// convert E^\mu and B^\mu to WHICHVEL v^i and B^i (code's primitive quantities)
// input E and B are defined as in GRFFE paper
int EBtopr(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE *pr)
{
  void EBtoT(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE (*T)[NDIM]);
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  FTYPE T[NDIM][NDIM];
  FTYPE U[NPR];
  FTYPE alpha;

  // basically convert T and B -> U

  EBtoT(Econ,Bcon,geom,T);

  alpha = sqrt(1./(- geom->gcon[GIND(0,0)]));  

  U[RHO]=0.0;
  U[UU]=T[TT][TT]; // T^t_\mu
  U[U1]=T[TT][RR];
  U[U2]=T[TT][TH];
  U[U3]=T[TT][PH];
  U[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  U[B2]=-Bcon[2]/(-alpha);
  U[B3]=-Bcon[3]/(-alpha);

  if(Utoprim_ffde(U, geom, pr)>=1){
    dualfprintf(fail_file, "EBtopr(Utoprim_ffde): error1\n");
    return(1);
  }


  return(0);


}

// convert E^\mu and B^\mu to T
// input E and B are defined as in GRFFE paper
void EBtoT(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE (*T)[NDIM])
{
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE (*T)[NDIM]);
  FTYPE Fit[NDIM];
  FTYPE U[NPR];
  FTYPE alpha;
  

  alpha = sqrt(1./(- geom->gcon[GIND(0,0)]));  

  Fit[TT]=0;
  Fit[1]=Econ[1]/(-alpha); // F^{it}
  Fit[2]=Econ[2]/(-alpha);
  Fit[3]=Econ[3]/(-alpha);
  U[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  U[B2]=-Bcon[2]/(-alpha);
  U[B3]=-Bcon[3]/(-alpha);
  
  
  FitUtoT(Fit,U,geom,T);


}

// v.B = KK B^2/\sqrt{B^2}
// component of velocity along field (assumed to be zero in GRFFE formulation)
void computeKK(FTYPE *pr, struct of_geom *geom, FTYPE *KK)
{
  FTYPE ucon[NDIM],ucov[NDIM],Bcon[NDIM],Bcov[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE Bsq;
  int j;

  ucon_calc(pr,geom,ucon,others);
  lower_vec(ucon,geom,ucov);

  SLOOPA(j) Bcon[j]=pr[B1+j-1];
  lower_vec(Bcon,geom,Bcov);

  Bsq=0.0;
  SLOOPA(j) Bsq+=Bcon[j]*Bcov[j];

  *KK=0.0;
  SLOOPA(j) *KK+=(ucov[j]/ucov[TT])*Bcon[j];
  *KK *=sqrt(1.0/Bsq);
  

}


// for E and B in arbitrary frame moving with coodinate 3-velocity veta^i
void EBvetatopr(FTYPE *Econ, FTYPE *Bcon, FTYPE *veta, struct of_geom *geom, FTYPE *pr)
{
  int j,k,l,m;
  int pl,pliter;
  FTYPE T[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM],Fcov[NDIM][NDIM],Fud[NDIM][NDIM];
  FTYPE Mcon[NDIM][NDIM];
  FTYPE prim[NPR];
  FTYPE Bcov[NDIM];
  FTYPE etacon[NDIM],etacov[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE Fsq;
  FTYPE alpha,etacovzamo[NDIM],etaconzamo[NDIM];
  FTYPE beta[NDIM];
  FTYPE U[NPR];

  extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
  void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE (*T)[NDIM]);
  void MtoF(int which, FTYPE (*Max)[NDIM],struct of_geom *geom, FTYPE (*faraday)[NDIM]);
  //\begin{equation}\label{faraday}
  //F^{\alpha\beta} \equiv f\eta^\alpha E^\beta - f\eta^\beta E^\alpha -
  //h B_\gamma \eta_\delta \epsilon^{\alpha\beta\gamma\delta}
  //\end{equation}
  //where $f$ and $h$ are arbitrary independent constants that we set to
  //be $f=h=1$ consistent with the conventions in MTW (see also, e.g.,

  //\begin{equation}\label{tmunu}
  //T^{\mu\nu} = {F^\mu}_\lambda F^{\nu\lambda} - \frac{1}{4}g^{\mu\nu}
  //F^{\lambda\kappa} F_{\lambda\kappa} ,
  //\end{equation}

  PLOOP(pliter,pl) prim[pl]=0.0;
  prim[U1]=veta[1];
  prim[U2]=veta[2];
  prim[U3]=veta[3];

  // arbitrary frame where Bcon and Econ are measured
  ucon_calc_3vel(prim,geom,etacon,others);
  lower_vec(etacon,geom,etacov);

  // B_\mu
  Bcon[TT]=0.0; // force
  lower_vec(Bcon,geom,Bcov);

  Econ[TT]=0.0; // force


  // F^{\mu\nu}
  DLOOP(j,k) Fcon[j][k]=0.0;
  DLOOP(j,k) for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      Fcon[j][k]+=-Bcov[l]*etacov[m]*lc4(1,geom->gdet,j,k,l,m); // \epsilon^{jklm}
    }
  DLOOP(j,k) Fcon[j][k] += etacon[j]*Econ[k]-etacon[k]*Econ[j];

  // F_mu^\nu
  DLOOPA(j) lower_vec(Fcon[j],geom,Fud[j]);

  // F_{\mu\nu}
  lower_A(Fcon,geom,Fcov);



  // T^\mu_\nu
  DLOOP(j,k) T[j][k]=0.0;

  Fsq=0.0;
  DLOOP(j,k) Fsq+=Fcon[j][k]*Fcov[j][k];

  DLOOP(j,k) {
    T[j][k]+=-0.25*delta(j,k)*Fsq;
    for(l=0;l<NDIM;l++) T[j][k]+=(-Fud[l][j])*(Fud[k][l]);
  }


  //  DLOOP(j,k) dualfprintf(fail_file,"Fcov[%d][%d]=%15.7g\n",j,k,Fcov[j][k]);

  MtoF(3,Fcov,geom,Mcon);

  U[RHO]=0.0;
  DLOOPA(j) U[UU+j]=T[TT][j]; // T^t_\mu
  SLOOPA(j) U[B1+j-1]=Mcon[j][TT]; // Jon's B^i = \dF^{it}
  
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"U[%d]=%21.15g\n",pl,U[pl]);

  Utoprim_ffde(U, geom, pr);

}


// needed in case T^t_t isn't consistent with T^i_t
void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE (*T)[NDIM])
{
  FTYPE prffde[NPR];
  FTYPE Mcon[NDIM][NDIM] ; /* contravariant Maxwell */
  void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM]);
  void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE (*Mcon)[NDIM]);


  // get initial U (only from E and B, not from T^t_t)
  prffde[U1]=Fit[1];
  prffde[U2]=Fit[2];
  prffde[U3]=Fit[3];
  prffde[B1]=-U[B1];  // M^{ti} since  U(geomfree)=B^i = \dF^{it} = M^{it}
  prffde[B2]=-U[B2];
  prffde[B3]=-U[B3];
  Max_con(prffde, geom, Mcon);
  /* T^i_t terms (i up, t down) */
  ffdestresstensor(Mcon, geom, T);
}


// 1=x, 2=y, 3=z
#define DEFAULTKILL 2 // choose T^t_y to kill

// convert T^t_t,T^t_y to T^t_x , where y is 2 spatial dimensions and x is one spatial dimension.  Eliminates one spatial dimension for T^t_t.
int TtoT(FTYPE BsqOalphasq, FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *Tout)
//int TtoT(FTYPE Bsq,FTYPE alphasq, FTYPE *Bcon, FTYPE *Uin, struct of_geom *geom, FTYPE *Uout)
{
  void getABC(int which, FTYPE BsqOalphasq,FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *AA, FTYPE*BB, FTYPE*CC);
  int j,k;
  FTYPE subEsq,EsqOBsq;
  FTYPE det,sqrtdet;
  FTYPE AA,BB,CC,qq;
  FTYPE replacea,replaceb;
  int which;

  // only do this if E^2 not close to 0
  subEsq=0.0; SLOOPA(j) subEsq+=Tin[j]*geom->gcon[GIND(TT,j)];
  EsqOBsq=2.0*(-Tin[TT]/alphasq+subEsq)/BsqOalphasq-1.0;
  // see if not enough information to change T
  // also check if violating causality, in which case dissipation is taken care of elsewhere.
  if( ((EsqOBsq)<1E-13)||(EsqOBsq>1.0) ) return(0); 


  // choose which T^t_{which} to remove 
  which=DEFAULTKILL;

  // set default as out=in
  DLOOPA(j) Tout[j]=Tin[j]; // j=0..3

  // get quadratic coefficients for quadratic equation with solution T^t_{which} = (-BB\pm\sqrt{BB^2-4*AA*CC})/(2*AA) .
  getABC(which, BsqOalphasq,alphasq,BconOalpha,Tin,geom,&AA,&BB,&CC);

  det=BB*BB-4.0*AA*CC;
  if(det<0.0){
    dualfprintf(fail_file,"t=%g i=%d j=%d : TtoT det=%g\n",t,geom->i,geom->j,det);
    myexit(1);
  }
  if(AA==0.0){
    dualfprintf(fail_file,"t=%g i=%d j=%d : TtoT AA=%g\n",t,geom->i,geom->j,AA);
    myexit(1);
  }

  sqrtdet=sqrt(det);

  if(BB==0.0){
    replacea=(-BB+sqrtdet)/(2.0*AA);
    replaceb=(-BB-sqrtdet)/(2.0*AA);
  }
  else{
    qq=-0.5*(BB+sign(BB)*sqrtdet);
    if(BB>0.0){
      replacea=CC/qq; // + root
      replaceb=qq/AA; // - root
    }
    else{
      replacea=qq/AA; // + root
      replaceb=CC/qq; // - root
    }
  }

  // figure out which one is closest to original (needed, e.g., inside ergosphere)
  if(fabs(replacea-Tin[which])<fabs(replaceb-Tin[which])){
    Tout[which]=replacea;
  }
  else Tout[which]=replaceb;

  return(0);

}

// apparently this or the above quadratic solution are written to decrease precision of inversion.  Large cancellations?
// even worse if one velocity direction dominates.  Remaining terms have large relative error.

// this gives coefficients to:
// AA (T^t_x)^2 + BB (T^t_x) + CC == 0
// where x=which
void getABC(int which, FTYPE BsqOalphasq,FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *AA, FTYPE*BB, FTYPE*CC)
{
  int j,k;
  FTYPE Q[NDIM][NDIM];
  FTYPE subBB,subCC1,subCC2;



  SSLOOP(j,k) Q[j][k]=(geom->gcon[GIND(j,k)])+(alphasq*geom->gcon[GIND(TT,j)]*geom->gcon[GIND(TT,k)])-BconOalpha[j]*BconOalpha[k]/BsqOalphasq;

  *AA=-1.0/(2.0*BsqOalphasq)*Q[which][which];

  subBB=0; SLOOPA(j) if(j!=which) subBB+=Tin[j]*Q[which][j];
  *BB=alphasq*(geom->gcon[GIND(TT,which)]) - subBB/BsqOalphasq;

  subCC1=0; SLOOPA(j) if(j!=which) subCC1+=geom->gcon[GIND(TT,j)]*Tin[j];
  subCC2=0; SSLOOP(j,k) if((j!=which)&&(k!=which)) subCC2+=Tin[j]*Tin[k]*Q[j][k];
  *CC=-Tin[TT] + alphasq*subCC1 - alphasq*BsqOalphasq*0.5 - subCC2/(2.0*BsqOalphasq);
  
}





// assumes geometry-free U in standard form
int Utoprim_ffde_old(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  int i,j,k;
  FTYPE T[NDIM][NDIM];
  FTYPE Ti0[NDIM] ;   /* stress tensor terms T^i_t */
  FTYPE Fit[NDIM],prffde[NPR];
  FTYPE Mcon[NDIM][NDIM] ; /* contravariant Maxwell */
  FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
  FTYPE ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  FTYPE prim[NPR];
  FTYPE Fcov[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM];
  FTYPE Bsq;
  FTYPE vcon[NDIM];
  void UtoFit(FTYPE *U, struct of_geom *geom, FTYPE *Fit);
  void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE (*Mcon)[NDIM]);
  void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);
  void MtoF(int which, FTYPE (*Max)[NDIM],struct of_geom *geom, FTYPE (*faraday)[NDIM]);
  void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM]);
  void raise_A(FTYPE (*Acov)[NDIM], struct of_geom *geom, FTYPE (*Acon)[NDIM]);
  FTYPE alpha;
  FTYPE betai[NDIM];
  FTYPE etacov[NDIM],etacon[NDIM];
  FTYPE Ecov[NDIM],Bcov[NDIM],Econ[NDIM],Bcon[NDIM];

  // lapse
  // define \eta_\alpha
  // assume always has 0 value for space components
  alpha = 1./sqrt(- geom->gcon[GIND(0,0)]);  

  etacov[TT]=-alpha; // any constant will work.
  SLOOPA(j) etacov[j]=0.0; // must be 0

  // shift
  /* beta */
  //      SLOOPA(j) beta[j] = geom->gcon[GIND(0,j)]*alpha*alpha ;
  // define \eta^\beta
  raise_vec(etacov,geom,etacon);

  SLOOPA(j) betai[j]=-etacon[j]*alpha;


  // get F^{it}
  UtoFit(U,geom,Fit);



  // get M^{\mu\nu}
  // correct order, sign, etc. for Max_con
  prffde[U1]=Fit[1];
  prffde[U2]=Fit[2];
  prffde[U3]=Fit[3];
  prffde[B1]=-U[B1];  // M^{ti} since  U(geomfree)=B^i = \dF^{it} = M^{it}
  prffde[B2]=-U[B2];
  prffde[B3]=-U[B3];

  Max_con(prffde, geom, Mcon);

  // for alternative formulation
  /* T^i_t terms (i up, t down) */
  ffdestresstensor(Mcon, geom, T);
  SLOOPA(j) Ti0[j] = T[j][TT] ; // T^i_t

  /* get covariant Maxwell from contravariant */
  lower_A(Mcon, geom, Mcov) ;

  // get F_{\mu\nu}
  MtoF(0,Mcon,geom,Fcov);

  /* get covariant Maxwell from contravariant */
  raise_A(Fcov, geom, Fcon) ;


  // now can get velocity
  Bsq=0.0;
  DLOOPA(j) Bsq+=-Mcon[j][TT]*Mcov[j][TT]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.

  // v^i = E_j B_k [ijk]/(B^l B_l) in Jon notation

  // S^i = -[ijk] E_j B_k (in Jon notation)
  // S_i = [ijk] E^j B^k = [ijk] \detg F^{tj} \dF^{kt} (Jon notation)
  // T^t_i = - sqrt(-g) * [ijk] * M^tj * F^kt = S_i

  // so v^i = -S^i/(B^l B_l) where S^i = T^i_t (can construct from T^{\mu\nu}) and compare with E/B form


  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  // that is, the below results in \dF_{it} v^i = 0 and F_{it} v^i = 0
  // notice that for v^i that \eta_t constant!=0 doesn't matter (cancels out)
#if(1)
  DLOOPA(j) Ecov[j]=0;
  DLOOP(j,k) Ecov[j] += etacon[k]*Fcov[j][k];

  DLOOPA(j) Bcov[j]=0;
  DLOOP(j,k) Bcov[j] += etacon[k]*Mcov[k][j];

  Bsq=0.0;
  SLOOPA(j) Bsq+=-alpha*prffde[B1+j-1]*Bcov[j];

  vcon[1]=-betai[1]+alpha*alpha*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])/(geom->gdet*Bsq);
  vcon[2]=-betai[2]+alpha*alpha*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])/(geom->gdet*Bsq);
  vcon[3]=-betai[3]+alpha*alpha*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])/(geom->gdet*Bsq);
#elif(0)
  DLOOPA(j) Econ[j]=0;
  DLOOP(j,k) Econ[j] += etacov[k]*Fcon[j][k];
  lower_vec(Econ,geom,Ecov);

  DLOOPA(j) Bcon[j]=0;
  DLOOP(j,k) Bcon[j] += etacov[k]*Mcon[k][j];
  lower_vec(Bcon,geom,Bcov);

  Bsq=0.0;
  DLOOPA(j) Bsq+=Bcon[j]*Bcov[j]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.

  vcon[1]=-betai[1]+alpha*alpha*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])/(geom->gdet*Bsq);
  vcon[2]=-betai[2]+alpha*alpha*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])/(geom->gdet*Bsq);
  vcon[3]=-betai[3]+alpha*alpha*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])/(geom->gdet*Bsq);
#elif(0)
  vcon[1]=(Fcov[TT][2]*Mcov[3][TT]-Fcov[TT][3]*Mcov[2][TT])/(geom->gdet*Bsq);
  vcon[2]=(Fcov[TT][3]*Mcov[1][TT]-Fcov[TT][1]*Mcov[3][TT])/(geom->gdet*Bsq);
  vcon[3]=(Fcov[TT][1]*Mcov[2][TT]-Fcov[TT][2]*Mcov[1][TT])/(geom->gdet*Bsq);
#elif(0) 

  // seems to kinda work

  Bsq=0.0;
  SLOOPA(j) Bsq+=-prffde[B1+j-1]*prffde[B1+j-1];

  vcon[1]=(Fcov[TT][2]*Mcon[3][TT]-Fcov[TT][3]*Mcon[2][TT])/(geom->gdet*Bsq);
  vcon[2]=(Fcov[TT][3]*Mcon[1][TT]-Fcov[TT][1]*Mcon[3][TT])/(geom->gdet*Bsq);
  vcon[3]=(Fcov[TT][1]*Mcon[2][TT]-Fcov[TT][2]*Mcon[1][TT])/(geom->gdet*Bsq);
#elif(0)
  SLOOPA(j) vcon[j]=0;
  SLOOP(j,k) vcon[j]+=((-Fcov[TT][k]/geom->gdet)*geom->gdet*Fcon[j][k])/(Bsq);
#elif(0)
  SLOOPA(j) vcon[j]=-Ti0[j]/Bsq;
#endif

  // convert from v^i -> u^\mu
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  //ucon_calc_3vel(prim,geom,ucon,others);
  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  // now assign field primitives
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];
  

  // done.

  return(0);




}

/*
 * Convert conserved quantities back to primitives
 */

// convert T^t_i & M^{it} -> F^{it}
void UtoFit(FTYPE *U, struct of_geom *geom, FTYPE *Fit)
{
  int k ;
  int j;

  FTYPE Mt[NDIM] ;
  FTYPE Bsq ;
  FTYPE T0[NDIM] ;   /* stress tensor terms T^t_i */
  FTYPE Mx[NDIM] ;   /* mixed Maxwell components M^t_i */


  Mt[TT]=0.0;
  SLOOPA(j) Mt[j] = -U[B1+j-1] ;  /* M^{ti} , where U[Bi]=B^i=M^{it}*/
  // T0's time term not needed
  SLOOPA(j) T0[j] = U[U1+j-1] ; /* T^t_i terms (t up, i down) */

  /* lower second index of Maxwell to get M^t_i */
  lower_vec(Mt,geom,Mx);

  // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
  Bsq=0.0;
  DLOOPA(j) Bsq+=Mx[j]*Mt[j];
      
  /* calculate F^it
   *
   * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
   * cyclic 1,2,3.
   *
   *
   * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
   *
   * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
   *
   * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
   *              = (1/\detg) [ijk] \dF^t_j T^t_k
   *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
   *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
   *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
   */
      
  /* F^{it} */
  Fit[0] = 0;
  Fit[1] = (Mx[2]*T0[3] - Mx[3]*T0[2])/(Bsq*geom->gdet) ;
  Fit[2] = (Mx[3]*T0[1] - Mx[1]*T0[3])/(Bsq*geom->gdet) ;
  Fit[3] = (Mx[1]*T0[2] - Mx[2]*T0[1])/(Bsq*geom->gdet) ;


}      


//  realBsqOalphasq=BsqOalphasq; // default
void U2BBsq(FTYPE *U, struct of_geom *geom, FTYPE *BconOalpha, FTYPE *BcovOalpha, FTYPE *BsqOalphasq)
{
  int j;
  // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
  // remove alpha since cancels below
  //SLOOPA(j) Bcov[j]=-alpha*Mx[j];
  //  SLOOPA(j) BcovOalpha[j]=-Mx[j]; // Mx = \dF^t_j

  // check B^2-E^2
  BconOalpha[TT]=0.0;
  // Bcon without alpha
  SLOOPA(j) BconOalpha[j] = U[B1+j-1] ;  /* Bcon=-M^{ti} , where U[Bi]=B^i=M^{it}*/

  lower_vec(BconOalpha,geom,BcovOalpha);

  // remove alpha^2 since cancels below
  //  realoBsqgeom=oBsqgeom/alphasq;
  *BsqOalphasq=SMALL; // in case not defined
  SLOOPA(j) *BsqOalphasq+=BcovOalpha[j]*BconOalpha[j];



}



void TtoEconOalpha(FTYPE *T0, FTYPE *BcovOalpha, FTYPE BsqOalphasq, struct of_geom *geom, FTYPE *EconOalpha)
{
  FTYPE oBsqOalphasqgeom;

  oBsqOalphasqgeom=1.0/(BsqOalphasq*geom->gdet); // 1/(\detg B^2/\alpha^2)
      
  /* calculate F^it
   *
   * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
   * cyclic 1,2,3.
   *
   *
   * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
   *
   * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
   *
   * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
   *              = (1/\detg) [ijk] \dF^t_j T^t_k
   *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
   *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
   *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
   */
      
  // E^i/\alpha = (\alpha^2/(\detg B^2)) ([ijk] (B_j/\alpha) T^t_k)
  EconOalpha[0] = 0;
  EconOalpha[1] = (BcovOalpha[2]*T0[3] - BcovOalpha[3]*T0[2])*(oBsqOalphasqgeom) ;
  EconOalpha[2] = (BcovOalpha[3]*T0[1] - BcovOalpha[1]*T0[3])*(oBsqOalphasqgeom) ;
  EconOalpha[3] = (BcovOalpha[1]*T0[2] - BcovOalpha[2]*T0[1])*(oBsqOalphasqgeom) ;


}      



void UtoFitMxBsqgeom(FTYPE *U, struct of_geom *geom, FTYPE *Fit, FTYPE*Mx, FTYPE*oBsqgeom)
{
  int k ;
  int j;

  FTYPE Mt[NDIM] ;
  FTYPE T0[NDIM] ;   /* stress tensor terms T^t_i */
  //      FTYPE Mx[NDIM] ;   /* mixed Maxwell components M^t_i */
  // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
  // -> B_\mu = -alpha * M^t_\mu
  // Real Bsq = Bsqgeom*alpha^2
  FTYPE Bsqgeom;

  Mt[TT]=0.0;
  SLOOPA(j) Mt[j] = -U[B1+j-1] ;  /* M^{ti} , where U[Bi]=B^i=M^{it}*/
  // T0's time term not needed
  SLOOPA(j) T0[j] = U[U1+j-1] ; /* T^t_i terms (t up, i down) */

  /* lower second index of Maxwell to get M^t_i */
  lower_vec(Mt,geom,Mx);

  // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
  // RealBsq = Bsq*alpha^2
  Bsqgeom=0.0;
  SLOOPA(j) Bsqgeom+=Mx[j]*Mt[j];
  Bsqgeom*=geom->gdet;
  *oBsqgeom=1.0/Bsqgeom;
      
  /* calculate F^it
   *
   * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
   * cyclic 1,2,3.
   *
   *
   * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
   *
   * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
   *
   * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
   *              = (1/\detg) [ijk] \dF^t_j T^t_k
   *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
   *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
   *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
   */
      
  /* F^{it} */
  Fit[0] = 0;
  Fit[1] = (Mx[2]*T0[3] - Mx[3]*T0[2])*(*oBsqgeom) ;
  Fit[2] = (Mx[3]*T0[1] - Mx[1]*T0[3])*(*oBsqgeom) ;
  Fit[3] = (Mx[1]*T0[2] - Mx[2]*T0[1])*(*oBsqgeom) ;


}      

/* 
 * convert primitives to conserved quantities 
 */

void primtoU_ffde(FTYPE pr[NPR], struct of_geom *geom, FTYPE U[NPR])
{

  int pl,pliter;
  FTYPE F1t,F2t,F3t,Mt1,Mt2,Mt3 ;

  /* load prims */
  F1t = pr[U1] ;   /* F^it -- E components */
  F2t = pr[U2] ;
  F3t = pr[U3] ;
  Mt1 = pr[B1] ;   /* M^ti -- B components */
  Mt2 = pr[B2] ;
  Mt3 = pr[B3] ;

  /* M^ti -- just copy!*/
  U[0] = Mt1 ;
  U[1] = Mt2 ;
  U[2] = Mt3 ;

  /* T^t_i terms( t up i down):
   *
   * T^t_i = - sqrt(-g) * [ijk] * M^tj * F^kt
   *
   */

  U[3] = - geom->gdet * (Mt2*F3t - Mt3*F2t) ;
  U[4] = - geom->gdet * (Mt3*F1t - Mt1*F3t) ;
  U[5] = - geom->gdet * (Mt1*F2t - Mt2*F1t) ;

  /* mult by geometric factor */
  PLOOP(pliter,pl) U[pl] *= geom->gdet ;
}


/*
 * load stress tensor and contravariant Maxwell in structure state
 */

void get_state_ffde(FTYPE pr[NPR], struct of_geom *geom, struct of_state *q)
{
  /* get contravariant Maxwell and stress tensor */
  //      Max_con(pr, geom ,q->Mcon) ;
  //      stresstensor(q->Mcon, geom, q->T) ;
}

// should be able to compare this with mhd_calc
// i.e.
// get_geometry()
// get_state
// mhd_calc(pr,dir,geom,q,mhd,NULL) where T^{dir}_j is returned
//


/*
 * lower both indices of an anti-symmetric  tensor
 *  -- go from contravariant to covariant
 */

/*
 * calculate the contravariant Maxwell
 */
void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE (*Mcon)[NDIM])
{
  int i,j,k ;

  FTYPE Fit[NDIM];
  FTYPE Bcon[NDIM];
  FTYPE Econ[NDIM], Ecov[NDIM] ;          /* E four-vector */
  FTYPE etacon[NDIM],etacov[NDIM];
  FTYPE alpha;

  // lapse
  // define \eta_\alpha
  // assume always has 0 value for space components
  alpha = 1./sqrt(- geom->gcon[GIND(0,0)]);  

  etacov[TT]=-alpha; // any constant will work.
  SLOOPA(j) etacov[j]=0.0; // must be 0

  // shift
  /* beta */
  //      SLOOPA(j) beta[j] = geom->gcon[GIND(0,j)]*alpha*alpha ;
  // define \eta^\beta
  raise_vec(etacov,geom,etacon);

  Fit[0] = 0.0; // Ftt=0
  Fit[1] = prffde[U1] ; /* F^{it} = -E^i/\detg in Jon notation */
  Fit[2] = prffde[U2] ;
  Fit[3] = prffde[U3] ;

  // B^\nu = \eta_\mu \dF^{\mu\nu}
  Bcon[0] = 0.0; // assumes \eta_\alpha={CONST,0,0,0}
  Bcon[1] = etacov[TT]*prffde[B1] ; /* M^{ti} =-B^i in Jon notation */
  Bcon[2] = etacov[TT]*prffde[B2] ;
  Bcon[3] = etacov[TT]*prffde[B3] ;

  /* create E^{\mu} = \eta_\beta F^{\alpha\beta} vector */

  DLOOPA(j) Econ[j]=etacov[TT]*Fit[j];
  lower_vec(Econ, geom, Ecov) ;

  /* diagonal */
  Mcon[0][0] = Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0. ;

  /* out with the old */
  //DLOOP(j,k)  Mcon[j][k] = 0. ;

  /* M^ij terms:
   *
   * M^ij = -beta^i M^tj + beta^j M^ti
   * + alpha * (1/sqrt(-g)) * Ecov[k]
   *
   * use primitive variables for M's
   *
   */

  // again assume \eta_\delta={CONST,0,0,0}
  //
  // \dF^{\alpha\beta} = -eta^\alpha B^\beta + \eta^\beta B^\alpha - E_\gamma \eta_\delta e^{\alpha\beta\gamma\delta}
  //                   = -eta^\alpha B^\beta + \eta^\beta B^\alpha + E_i \eta_t e^{t \alpha\beta i} // 3 sign switches from t
  // \dF^{jk}          = -eta^j B^k + \eta^k B^j - (1/\detg) E_i \eta_t [jki] // 1 sign switch from e^
  // \dF^{jk}          = -eta^j B^k + \eta^k B^j - (\eta_t/\detg) (E_i [jki])
  j=1;k=2;i=3; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->gdet)*Ecov[i];
  j=2;k=3;i=1; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->gdet)*Ecov[i];
  j=3;k=1;i=2; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->gdet)*Ecov[i];


 
  /* copy remaining spacial terms */
   
  Mcon[2][1] = -Mcon[1][2] ;
  Mcon[3][2] = -Mcon[2][3] ;
  Mcon[1][3] = -Mcon[3][1] ;

  /* time terms - easy */
  SLOOPA(j) {
    Mcon[TT][j]=Bcon[j]/etacov[TT];
    Mcon[j][TT]=-Mcon[TT][j];
  }
}


/*
 * calculate the contravariant Maxwell
 */
void Max_con_old(FTYPE prffde[NPR], struct of_geom *geom, FTYPE (*Mcon)[NDIM])
{
  int j ;

  FTYPE F1t, F2t, F3t, Mt1, Mt2, Mt3 ;    /* prims */
  FTYPE Econ[NDIM], Ecov[NDIM] ;          /* E four-vector */
  FTYPE beta[NDIM] ;                      /* shift  3-vector  */
  FTYPE alpha;


  /* lapse */
  alpha = 1./sqrt(- geom->gcon[GIND(0,0)]);  

  /* beta */ 
  // GODMARK: wrongly used below.  Should contract Mt with \eta^\mu = 1/\alpha (1,-\beta^i), so missing 1/\alpha
  // but then left off -alpha on Mt terms in Mcon below, so Mcon below ends up being correct (including signs)
  SLOOPA(j) beta[j] = geom->gcon[GIND(0,j)]*alpha*alpha ;

  F1t = prffde[U1] ; /* F^{it} = -E^i/\detg in Jon notation */
  F2t = prffde[U2] ;
  F3t = prffde[U3] ;
  Mt1 = prffde[B1] ; /* M^{ti} =-B^i in Jon notation */
  Mt2 = prffde[B2] ;
  Mt3 = prffde[B3] ;

  /* create E^{\mu} = \eta_\beta F^{\alpha\beta} vector */
  Econ[0] = 0. ;
  Econ[1] = - alpha * F1t ;
  Econ[2] = - alpha * F2t ;
  Econ[3] = - alpha * F3t ;

  lower_vec(Econ, geom, Ecov) ;

  /* diagonal */
  Mcon[0][0] = Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0. ;

  /* out with the old */
  //DLOOP(j,k)  Mcon[j][k] = 0. ;

  /* M^ij terms:
   *
   * M^ij = -beta^i M^tj + beta^j M^ti
   * + alpha * (1/sqrt(-g)) * Ecov[k]
   *
   * use primitive variables for M's
   *
   */

  // different sign! for E part (GODMARK) (see above comments -- just crazy mixing of terms leads to same results)
  Mcon[1][2] = (-beta[1] * Mt2 + beta[2] * Mt1 + alpha/(geom->gdet) * Ecov[3]) ;
  Mcon[2][3] = (-beta[2] * Mt3 + beta[3] * Mt2 + alpha/(geom->gdet) * Ecov[1]) ;
  Mcon[3][1] = (-beta[3] * Mt1 + beta[1] * Mt3 + alpha/(geom->gdet) * Ecov[2]) ;
  
  /* copy remaining spacial terms */
   
  Mcon[2][1] = -Mcon[1][2] ;
  Mcon[3][2] = -Mcon[2][3] ;
  Mcon[1][3] = -Mcon[3][1] ;

  /* time terms - easy */
  Mcon[0][1] = Mt1 ;
  Mcon[0][2] = Mt2 ;
  Mcon[0][3] = Mt3 ;
  Mcon[1][0] = - Mt1 ;
  Mcon[2][0] = - Mt2 ;
  Mcon[3][0] = - Mt3 ;
}



// see also phys.ffde.debug.c

// GODMARK : only really 2D test
// assumes geometry exists
void testffdeinversion(void)
{
  int i,j,k,pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  FTYPE Bu[NDIM],uu[NDIM];
  int itest,jtest,ktest;


#if(0)
  ktest=0;
  itest=50;
  jtest=32;
  PLOOP(pliter,pl){
    prin[pl]=p[itest][jtest][pl];
    prin[U1]=prin[U3]*0.1;
    prin[U2]=prin[U2]*0.1;
    prin[B3]=prin[B1];
  }

  get_geometry(itest,jtest,ktest,CENT,ptrgeom);
  get_state(prin,ptrgeom,&q);

  primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

  PLOOP(pliter,pl) prin[pl]=prout[pl];

  primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP(pliter,pl) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]);

  k=0;
  i=48;
  j=32;
  PLOOP(pliter,pl){
    prin[pl]=MACP0A1(p,i,j,k,pl);
    prout[pl]=0.0;
  }

  get_geometry(i,j,k,CENT,ptrgeom);
  get_state(prin,ptrgeom,&q);
  primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP(pliter,pl) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]);



  // loop through to test consistency
  PLOOP(pliter,pl){
    prin[pl]=prout[pl];
    prout[pl]=0.0;
  }

  get_geometry(i,j,k,CENT,ptrgeom);
  get_state(prin,ptrgeom,&q);
  primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);

  PLOOP(pliter,pl) dualfprintf(fail_file,"U[%d]=%21.15g\n",pl,U[pl]);
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP(pliter,pl) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]);

#endif
  

  // now do loop over range of relative 4-velocity and field for \rho=u=0

  dualfprintf(fail_file,"realtest\n");




#if(0)
  k=0;
  i=32;
  j=32;
  for(uu[1]=-100.0;uu[1]<=100.0;uu[1]+=10)
    for(uu[2]=-100.0;uu[2]<=100.0;uu[2]+=10)
      for(uu[3]=-100.0;uu[3]<=010.0;uu[3]+=10)

        for(Bu[1]=-10.0;Bu[1]<=10.0;Bu[1]+=1)
          for(Bu[2]=-10.0;Bu[2]<=10.0;Bu[2]+=1)
            for(Bu[3]=-10.0;Bu[3]<=10.0;Bu[3]+=1){

              if((Bu[1]==0)&&(Bu[2]==0)&&(Bu[3]==0)) continue;

              prin[RHO]=prin[UU]=0;

              prin[U1]=uu[1];
              prin[U2]=uu[2];
              prin[U3]=uu[3];

              prin[B1]=Bu[1];
              prin[B2]=Bu[2];
              prin[B3]=Bu[3];

              get_geometry(i,j,k,CENT,ptrgeom);

              if(get_state(prin,ptrgeom,&q)>=1) dualfprintf(fail_file,"getstate failure in realtest\n");
              if(primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL)>=1) dualfprintf(fail_file,"primtoU failure in realtest\n");
       
              Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

              PLOOP(pliter,pl) prin[pl]=prout[pl];

              // only clean solution to test
              get_state(prin,ptrgeom,&q);
              primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);
       
              Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

              // just compare pr in and pr out.
              SLOOPA(j) dualfprintf(fail_file,"realtest uu[%d]=%21.15g\n",j,uu[j]);
              SLOOPA(j) dualfprintf(fail_file,"realtest Bu[%d]=%21.15g\n",j,Bu[j]);
              for(k=U1;k<=B3;k++) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]); 
              fflush(fail_file);

              //       myexit(0);
            }

#endif


#if(1)
  k=0;
  i=0;
  j=32;
  //  uu[1]=uu[2]=uu[3]=-1.9999E3; //yes
  //  uu[1]=uu[2]=uu[3]=-2.4E3; // no
  //  uu[1]=uu[2]=uu[3]=-2.28E3; // top? no?
  //uu[1]=uu[2]=uu[3]=-2.3E3; // works! (with old formulation)
  //  uu[1]=uu[2]=uu[3]=-2.30000099E3; // odd exact machine precision answer
  //uu[1]=uu[2]=uu[3]=-2.300000999E3; // doesn't work.

  // new formulation
  //  uu[1]=uu[2]=uu[3]=-3E4; //yes
  /*
  // I get 5E-5 error for the below with uu2=-1E5
  // 8E-3 error for uu2=-1E6
  // -0.35 error for uu2=-1E7
  uu[1]=0;
  uu[2]=-1E5;
  uu[3]=0;

  Bu[1]=0;
  Bu[2]=0;
  Bu[3]=-1;

  // even uu2=-1E10 works, but uu1=-1E10 doesn't.
  */
  uu[1]=-2270.477;
  uu[2]=-5.407722;
  uu[3]=0;

  Bu[1]=0;
  Bu[2]=-1;
  Bu[3]=0;
  prin[RHO]=prin[UU]=0;

  prin[U1]=uu[1];
  prin[U2]=uu[2];
  prin[U3]=uu[3];

  prin[B1]=Bu[1];
  prin[B2]=Bu[2];
  prin[B3]=Bu[3];

  SLOOPA(j) dualfprintf(fail_file,"realtest uu[%d]=%21.15g\n",j,uu[j]);
  SLOOPA(j) dualfprintf(fail_file,"realtest Bu[%d]=%21.15g\n",j,Bu[j]);

  get_geometry(i,j,k,CENT,ptrgeom);

  if(get_state(prin,ptrgeom,&q)>=1) dualfprintf(fail_file,"getstate failure in realtest\n");
  DLOOPA(pl) dualfprintf(fail_file,"1 uu[%d]=%21.15g\n",k,q.ucon[pl]);
  if(primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL)>=1) dualfprintf(fail_file,"primtoU failure in realtest\n");
       
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion
  for(pl=U1;pl<=B3;pl++) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]); 

  PLOOP(pliter,pl) prin[pl]=prout[pl];

  // only clean solution to test
  get_state(prin,ptrgeom,&q);
  DLOOPA(k) dualfprintf(fail_file,"2 uu[%d]=%21.15g\n",k,q.ucon[pl]);
  primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL);
       
  Utoprim_ffde(U,ptrgeom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  for(pl=U1;pl<=B3;pl++) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]); 
  fflush(fail_file);

  //       myexit(0);
#endif


  myexit(0);




}



