//#include "decs.h"

// u2p_defs already includes decs.h
#include "u2p_defs.h"
#include "utoprim_jon_eos.c"

///////////////////////////////
//
// cold GRMHD limit
//
// 



// questions:
// 1) In cold limit, energy should be derivable from momentum, right?  So can't evolve
// energy separately.  Have to evolve momentums.  Could modify slightly a momentum so consistent with
// how one evolved energy (tried in force-free).
// If no source terms, then both momentum and energy conserved exactly.




// invert from U (conserved quantities in cold GRMHD) to p (primitives in cold GRMHD)
// assumes no geometry in conserved quantity (as if coming from Utoprimgen())
int Utoprim_coldgrmhd(FTYPE *U, struct of_geom *geom, FTYPE *pr, int *positivityproblem)
{
  FTYPE alpha,alphasq;
  FTYPE etacov[NDIM];
  int jj;
  FTYPE Ueta[NPR];
  int construct_Ueta(FTYPE *etacov, FTYPE *U, struct of_geom *geom, FTYPE *Ueta);
  int construct_BsqDsqSSsq(FTYPE *U, struct of_geom *geom, FTYPE *Bsq, FTYPE *D, FTYPE *Dsq, FTYPE *S, FTYPE *Ssq, FTYPE *Msq);
  FTYPE Bsq,D,Dsq,S,Ssq,Msq;
  int get_cold_W_fromMsq(FTYPE Msq, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE *W);
  FTYPE W;
  int returnW;
  FTYPE gammasq;
  int get_cold_gammasq(FTYPE W, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE Msq, FTYPE *gammasq);
  FTYPE vcon[NDIM],vcov[NDIM];
  int get_cold_3vel(FTYPE *Ueta, FTYPE W, FTYPE Bsq, FTYPE S, FTYPE *vcon);
  FTYPE vsq;
  int get_properdensity(FTYPE D, FTYPE gammasq, FTYPE *rho);
  FTYPE rho;
  int get_cold_pressure(FTYPE *EOSextra, FTYPE W, FTYPE Dsq, FTYPE gammasq, FTYPE *pressure);
  FTYPE pressure;
  int get_eta_4vel(FTYPE *vcon, FTYPE gammasq, FTYPE *uetacon);
  FTYPE uetacon[NDIM];
  FTYPE ucon[NDIM];
  FTYPE others[NUMOTHERSTATERESULTS];
  int i,j,k,loc;
  

  i=geom->i;
  j=geom->j;
  k=geom->k;
  loc=geom->p;


  // assume going to be positive
  *positivityproblem=0;

  ////////////////////////
  //
  // get \eta observer
  //
  ///////////////////////

  // define \eta_\alpha
  // assume always has 0 value for space components
  alphasq = 1./(- geom->gcon[GIND(0,0)]);  

  alpha = sqrt(alphasq);
  etacov[TT]=-alpha; // any constant will work.
  SLOOPA(jj) etacov[jj]=0.0; // must be 0 for system to work


  // get Ueta
  construct_Ueta(etacov,U,geom,Ueta);

  // get B^2, D^2, S, S^2, and M^2
  construct_BsqDsqSSsq(Ueta, geom, &Bsq, &D, &Dsq, &S, &Ssq, &Msq);

  // get W
  // in cold case, Msq and EN both contain same and enough information to obtain W, but ultimately to obtain velocity components of m are used.
  //               hence, in general a self-consistent solution is obtained by using Msq and m's.
  //               Indeed, as mentioned in construct_Ueta, T^t_t~EN term is never actually used, so it would be inconsistent to construt Ueta and then get W from EN
  returnW=get_cold_W_fromMsq(Msq, Bsq, Dsq, Ssq, &W); // need to use Msq since that's what is ultimately used to generate the solution
  //  returnW=get_cold_W(Ueta[UU], Bsq, Dsq, Ssq, &W); // Can get cold W only from EN, but evolution probably uses momentum conserved quantities instead of energy
  if(returnW!=0){
    dualfprintf(fail_file,"Utoprim_coldgrmhd: No solution possible since W cannot be found\n");
    dualfprintf(fail_file,"myEN=%21.15g Msq=%21.15g Bsq=%21.15g Dsq=%21.15g Ssq=%21.15g\n",Ueta[UU],Msq,Bsq,Dsq,Ssq);
    *positivityproblem=1;
    return(1);
  }

  get_cold_gammasq(W,Bsq,Dsq,Ssq,Msq,&gammasq);
 

  get_cold_3vel(Ueta,W,Bsq,S,vcon);


  if(gammasq<0.0){
    // deadly
    lower_vec(vcon,geom,vcov);
    vsq =dot(vcon,vcov);

    dualfprintf(fail_file,"Utoprim_coldgrmhd: No solution possible since gammasq=%21.15g < 1.0\n",gammasq);
    dualfprintf(fail_file,"Utoprim_coldgrmhd: No solution possible since vsq=%21.15g < 0.0 or >1.0\n",vsq);
    *positivityproblem=2;
    return(1);
  }
  else if(gammasq<1.0){

    lower_vec(vcon,geom,vcov);
    vsq =dot(vcon,vcov);

    dualfprintf(fail_file,"Utoprim_coldgrmhd: No solution possible since gammasq=%21.15g < 1.0\n",gammasq);
    dualfprintf(fail_file,"Utoprim_coldgrmhd: No solution possible since vsq=%21.15g < 0.0 or >1.0\n",vsq);
    *positivityproblem=3; // not deadly, but screwy
  }


  ///////////////////////////////
  //
  // finally, get primitives
  //
  // proper density
  if(D<0.0) *positivityproblem+=10; // not deadly
  get_properdensity(D,gammasq,&pr[RHO]);

  // notice that pressure can go negative with gammasq>=1.0, but if gammasq<1.0 then pressure is imaginary
  // hence sequence of "failures" is first negative pressure and then gammasq<1.0
  // notice pressure is well-defined (positive) even for gammasq<1.0 but still > 0.0 !
  get_cold_pressure(GLOBALMAC(EOSextraglobal,i,j,k), W, Dsq, gammasq, &pressure);
  // check if perssure >=0
  if(pressure<0.0){
    dualfprintf(fail_file,"Utoprim_coldgrmhd: pressure=%21.15g < 0.0\n",pressure);
    *positivityproblem+=100;
    // but don't fail since negative pressure fixable in general by Runge-Kutta, floors, or other fixups
  }

  // get relative 4-velocity
  get_eta_4vel(vcon,gammasq,uetacon);

  if(WHICHVEL!=VELREL4){
    // then have to convert to 4-velocity and then to correct velocity
    ucon_calc_rel4vel_fromuconrel(uetacon,geom,ucon,others);
    ucon2pr(WHICHVEL,ucon,geom,pr);
  }

  // get code fields
  SLOOPA(jj) pr[B1-1+jj] = U[B1-1+jj]; // no need to convert from Ueta since originally in U


  return(0);
}


// construct energy, momentum, and field strengths as measured by \eta observer
int construct_Ueta(FTYPE *etacov, FTYPE *U, struct of_geom *geom, FTYPE *Ueta)
{
  FTYPE etacon[NDIM];
  FTYPE projetauu[NDIM][NDIM];
  int j,k;


  // get etacon
  raise_vec(etacov,geom,etacon);

  // get projection tensor
  DLOOP(j,k) projetauu[j][k] = geom->gcon[GIND(j,k)] + etacon[j]*etacon[k];


  // D_eta = -\rho u^\mu \eta_\mu = +\gamma \rho
  Ueta[RHO] = -U[RHO]*etacov[TT];

  // h = (\rho+u+p)/\rho
  // E_eta = T^\mu_\nu \eta^\nu \eta_\mu = (\rho+u+p)\gamma^2 - p = D_eta h \gamma - p
  // now start with 0
  Ueta[UU]=0.0;
  DLOOPA(j) Ueta[UU] += etacov[TT]*U[UU+j]*etacon[j];

  // m_eta^\alpha = -T^\mu_\nu \eta_\mu j^{\nu\alpha}
  // m_eta^0 = 0, so ignore time-like term
  // also, T^t_t not actually used here for this projection operator, so U[UU] not really used and could just iterate over space-space (j,k) instead of space-time (j,k)
  // now start with 0
  SLOOPA(j) Ueta[UU+j] = 0.0;
  // j: space, k:time+space
  SDLOOP(j,k) Ueta[UU+j] = -U[UU+k]*etacov[TT]*projetauu[k][j];

  // B_eta^\alpha = \eta_beta \dF^{\beta\alpha} = \eta_t \dF^{t\alpha} = \eta_t (-B^{\alpha t} )
  // j:space
  SLOOPA(j) Ueta[B1-1+j] = -etacov[TT]*U[B1-1+j];

  return(0);
}


// doesn't need Ueta[UU]
int construct_BsqDsqSSsq(FTYPE *Ueta, struct of_geom *geom, FTYPE *Bsq, FTYPE *D, FTYPE *Dsq, FTYPE *S, FTYPE *Ssq, FTYPE *Msq)
{
  FTYPE Betacon[NDIM],Betacov[NDIM]; // contra-field in \eta frame
  FTYPE metacon[NDIM],metacov[NDIM]; // contra-momentum in \eta frame
  int j;
  

  // D_eta^2
  *D = Ueta[RHO];
  *Dsq = (*D) * (*D) ;


  // m_eta^\alpha = D_eta h \tilde{u}^\alpha = D_eta h \gamma \tilde{v}^\alpha
  // \tilde{v}^\alpha = \tilde{u}^\alpha/\gamma
  // gamma = -\eta\cdot u
  metacon[TT] = 0.0;
  SLOOPA(j) metacon[j]=Ueta[U1-1+j];

  // get metacov
  lower_vec(metacon,geom,metacov);
  
  *Msq=dot(metacon,metacov);
  
  // B^\alpha
  Betacon[TT] = 0.0;
  SLOOPA(j) Betacon[j]=Ueta[B1-1+j];

  // get Betacov
  lower_vec(Betacon,geom,Betacov);

  *Bsq=dot(Betacon,Betacov);


  *S = dot(Betacon,metacov); // S = m.B
  *Ssq = (*S) * (*S) ; // S^2
  
  
  return(0);
}


// get cold value of W=\rho h \gamma^2
// h : specific enthalpy  h = (\rho+u+p)/\rho
// see grmhd_positivity.nb
//
// returns 0 if solution exists
// returns 1 if no solution because of discr<0
// returns 2 if no solution because of termundercuberoot<0

// could be problem where cubic solution uses other roots sometimes.  Might be avoided by judicious use of absolute value signs
int get_cold_W_fromEN(FTYPE EN, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE *W)
{
  FTYPE b,s,d;
  FTYPE terma,termb,discr,termundercuberoot;
  FTYPE termc;


  b=Bsq/EN;
  d=Dsq/pow(EN,2.0);
  s=Ssq/pow(EN,3.0);

  terma=27.0*(b*d+s)*0.25;
  termb=-pow( b-1.0 ,3.0);
  discr=terma*(2.0*termb+terma);
  if(discr<0.0){
    dualfprintf(fail_file,"get_cold_W_fromEN discr=%21.15g < 0.0\n",discr);
    return(1);
  }

  termundercuberoot=(termb+terma+sqrt(discr));
  if(termundercuberoot<0.0){
    dualfprintf(fail_file,"get_cold_W_fromEN termundercuberoot=%21.15g < 0.0\n",termundercuberoot);
    return(2);
  }
  
  termc=pow(termundercuberoot,1.0/3.0);

  *W = EN*THIRD*( (b-1.0)*(-1.0+(b-1.0)/termc) + termc );

  return(0);
}



// get cold value of W=\rho h \gamma^2
// h : specific enthalpy  h = (\rho+u+p)/\rho
// see grmhd_positivity.nb
//
// returns 0 if solution exists
// returns 1 if no solution because of terma<0
// returns 2 if no solution because of discr<0
// returns 3 if no solution because of discr2<0

// could be problem with quartic solution using other roots sometimes.  Might be avoided by use of absolute values in some places to mimic all (or just 1 other) roots
// might need to rearrange terms so non-imaginary


int get_cold_W_fromMsq(FTYPE Msq, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE *W)
{
  FTYPE b,s,d;
  FTYPE C1,C2;
  FTYPE terma,termb,termc,termc2,termd,terme,termf,termg;
  FTYPE discr,discr2;
  FTYPE b2,b4,b6,s2;


  b=Bsq/sqrt(Msq);
  d=Dsq/Msq;
  s=Ssq/pow(Msq,1.5);

  b2=b*b;
  b4=b2*b2;
  b6=b2*b4;
  s2=s*s;
  
 
  terma=-432.*(b - 1.*s)*(b*d + s)*(27.*b*(-1. + d)*s + 3.*(1. + (-7. + d)*d)*b2 - 3.*(1. + d)*b4 + b6 - 1.*pow(1. + d,3) + 27.*s2);

  termb=2.*(54.*b*(-1. + d)*s + 3.*(1. + (-16. + d)*d)*b2 - 3.*(1. + d)*b4 + b6 - 1.*pow(1. + d,3) + 54.*s2);

  termc=THIRD*(2. + 2.*d + b2);

  termc2 = 2.0*termc;

  termd=pow(2.0,THIRD)*pow(1. + d - 1.*b2,2);

  // terma>=0.0
  if(terma<0.0){
    dualfprintf(fail_file,"get_cold_W_fromMsq terma=%21.15g < 0.0\n",terma);
    return(1);
  }
  terme = 3.0*pow(termb+sqrt(terma),THIRD);

  termf = 2.*(b*(-1. + d) + 2.*s);
    
  C2 = -1.*b;

  C1 = 1.0/(9.0*pow(2,THIRD));

  discr=termc + termd/terme + C1*terme;
  if(discr<0.0){
    dualfprintf(fail_file,"get_cold_W_fromMsq discr=%21.15g < 0.0\n",discr);
    return(2);
  }
  termg=sqrt(discr);


  discr2=termc2 - (1.*termd)/terme - 1.*C1*terme + termf/termg;
  if(discr2<0.0){
    dualfprintf(fail_file,"get_cold_W_fromMsq discr2=%21.15g < 0.0\n",discr2);
    return(3);
  }

  // if passed 1,2,3 returns, then have solution
  *W=0.5*sqrt(Msq)*(C2 + sqrt(discr2) + termg);


  return(0);
}


// see grmhd_positivity.nb (just take Mignone & Bodo (2006) equation 16 and plug in gamma = W/D (true for cold case)
int get_cold_E_from_W(FTYPE W, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE Msq, FTYPE *EN)
{

  *EN = Bsq - (Bsq*Dsq+Ssq)/(2.0*W*W) + W ;

  return(0);
}


// see grmhd_positivity.nb (just take Mignone & Bodo (2006) equation 17 and plug in gamma = W/D (true for cold case)
int get_cold_Msq_from_W(FTYPE W, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE E, FTYPE *Msq)
{

  *Msq = (1.0-Dsq/(W*W))*(Bsq+W)*(Bsq+W) - Ssq*(Bsq+2.0*W)/(W*W) ;

  return(0);
}


// \gamma^2 = 1- (S^2(2W^2 + B^2) + M^2 W^2)^{1/2} / ( (W+B^2)^2 W^2 )
// Mignone & Bodo (2006) equation 18
// Noble et al. (2006) equation 31 or 28
int get_cold_gammasq(FTYPE W, FTYPE Bsq, FTYPE Dsq, FTYPE Ssq, FTYPE Msq, FTYPE *gammasq)
{
  FTYPE top,bottom,Wsq;

  Wsq=W*W;

  top = Msq*Wsq+Ssq*(Bsq+2.0*W);
  bottom = Wsq*(Bsq+W)*(Bsq+W);

  *gammasq = 1.0-top/bottom;

  return(0);
}


// Mignone & Bodo (2006) equation 23
// Noble et al. (2006) equation 31
int get_cold_3vel(FTYPE *Ueta, FTYPE W, FTYPE Bsq, FTYPE S, FTYPE *vcon)
{
  int j;

  vcon[TT] = 0.0;
  // v^\alpha = 1/(W+B^2) (m_eta^\alpha + S/W B^\alpha)
  // using spatial loop so only have to use Ueta
  SLOOPA(j) vcon[j] = 1.0/(W+Bsq) * ( Ueta[U1+j-1] + S/W * Ueta[B1+j-1] );
  

  return(0);
}


// if D is negative, then density can be negative, which is recoverable.
// but if gammasq<0, then not recoverable
int get_properdensity(FTYPE D, FTYPE gammasq,FTYPE *rho)
{
  FTYPE ftemp;

  *rho = D/sqrt(gammasq);

  return(0);
}


// assumes ideal gas EOS
// Migone & Bodo (2006) equation 19
// Noble et al. (2006) equation 34
int get_cold_pressure(FTYPE *EOSextra, FTYPE W, FTYPE Dsq, FTYPE gammasq, FTYPE *pressure)
{
  //  FTYPE gamr;
  FTYPE vsq;

  //  gamr = gam/(gam-1.0);
  //  *pressure = (W - sqrt(Dsq*gammasq))/(gamr*gammasq);

  vsq=1.0-1.0/gammasq;

  *pressure = pressure_W_vsq(COLDEOS, EOSextra, W,sqrt(fabs(Dsq)),vsq);

  return(0);
}



// same relative 4-velocity as in phys.c
// has properties that:
//
// \tilde{v}^\alpha = \tilde{u}^\alpha / \gamma
// where \gamma = \sqrt{1+u^2}
// So: u^2 = 1/ (1/v^2 - 1)
// So: \gamma = 1/\sqrt{1-v^2} (i.e. just like in special relativity)
//
int get_eta_4vel(FTYPE *vcon, FTYPE gammasq, FTYPE *uetacon)
{
  int j;

  // TT component is 0 and will be set accordingly from vcon[TT]=0
  DLOOPA(j) uetacon[j] = vcon[j] * sqrt(gammasq);

  return(0);
}



// trivial filter
void filter_coldgrmhd(int i, int j, int k, FTYPE *pr)
{
  // kill internal energy
  pr[UU]=0.0;

}


void test_coldgrmhd_inversion(void)
{
  int i,j,k,pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  FTYPE rho,Bu[NDIM],uu[NDIM];
  int itest,jtest,ktest;
  int positivityproblem;

  // pick location
  k=0;
  i=0;
  j=32;

  // choose primitive values
  rho=1.0;
  uu[1]=0.0;
  uu[2]=0.0;
  uu[3]=0.0;

  Bu[1]=0;
  Bu[2]=0;
  Bu[3]=0;

  // assign to prin
  prin[UU] =0.0;
  prin[RHO]=rho;
  prin[U1]=uu[1];
  prin[U2]=uu[2];
  prin[U3]=uu[3];

  prin[B1]=Bu[1];
  prin[B2]=Bu[2];
  prin[B3]=Bu[3];

  dualfprintf(fail_file,"realtest rho=%21.15g\n",rho);
  SLOOPA(j) dualfprintf(fail_file,"realtest uu[%d]=%21.15g\n",j,uu[j]);
  SLOOPA(j) dualfprintf(fail_file,"realtest Bu[%d]=%21.15g\n",j,Bu[j]);

  get_geometry(i,j,k,CENT,ptrgeom);

  if(get_state(prin,ptrgeom,&q)>=1) dualfprintf(fail_file,"getstate failure in realtest\n");
  DLOOPA(pl) dualfprintf(fail_file,"1 uu[%d]=%21.15g\n",k,q.ucon[pl]);
  if(primtoU(UNOTHING,prin,&q,ptrgeom,U, NULL)>=1) dualfprintf(fail_file,"primtoU failure in realtest\n");
       
  Utoprim_coldgrmhd(U,ptrgeom,prout,&positivityproblem); // no need for initial guess since analytic inversion
  for(pl=U1;pl<=B3;pl++) dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",pl,prin[pl],pl,prout[pl],(prin[pl]-prout[pl])/prin[pl]); 



  myexit(0);





}
