
/*! \file rescale_interp.c
     \brief Functions that rescale input primitives/conserves for interpolation purposes, and then functions that unrescale the result back to inputted type of quantities
     // GENERAL PRIMITIVE INTERPOLATION CHANGE OF VARIABLES
     // GODMARK: really should restrict rescale() to npr2interplist since may operate on quantity didn't want to change (to do this, just loop over PINTERPLOOP(pliter,pl) always and put conditional on each pl using if(pl==??) .
*/


#include "decs.h"


/// rescale quantities for interpolation
/// notice that input order is always:
/// which: DORESCALE = rescale UNRESCALE = unrescale
/// dir: which dimension (1,2,3)
/// pr: true primitive
/// p2interp: scaled primitive
/// this also allows any fancy remapping, such as characteristic interpolation -- rest of code is setup to allow any remapping as long as you have an inversion
int rescale(int which, int dir, FTYPE *pr, struct of_geom *ptrgeom,FTYPE *p2interp)
{
  FTYPE scale[NPR2INTERP],r,th,X[NDIM],V[NDIM];
  int ii,jj,kk;
  int pl,pliter;
  struct of_state q;
  extern int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE quasivsq;
  extern int limit_quasivsq(FTYPE quasivsqnew, struct of_geom *geom, FTYPE *pr);
  FTYPE vcon[NDIM],ucon[NDIM];
  int j;
  FTYPE newpr[NPR];
  FTYPE uconrel[NDIM],ucovrel[NDIM];
  FTYPE vconrel[NDIM];
  FTYPE normuconrel,normuconrel_fromui;
  FTYPE gamma,qsq;

#if(VARTOINTERPFIELD==PULSARFIELD)
  extern void getconsts(FTYPE *uconmetp, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconconst);
  extern void undoconsts(FTYPE *uconconst, FTYPE *V, struct of_geom *ptrgeom, FTYPE (*dxdxp)[NDIM],FTYPE *uconmetp);
#endif
  FTYPE Bconin[NDIM],Bconout[NDIM];
  FTYPE Uconin[NDIM],Uconout[NDIM];
  FTYPE dxdxp[NDIM][NDIM];


  int dopl[NPR2INTERP];
  for(pl=0;pl<NPR2INTERP;pl++) dopl[pl]=0; // default don't do
  PINTERPLOOP(pliter,pl) dopl[pl]=1; // only do if in loop range (do this instead of looping over actual rescalings below)

  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;

  coord_ijk(ii,jj,kk,ptrgeom->p,X);
  bl_coord_ijk(ii,jj,kk,ptrgeom->p,V);
  r=V[1]; th=V[2];

#if(VARTOINTERP==PRIMTOINTERP)

  //  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if RESCALEINTERP=%d\n",VARTOINTERP,RESCALEINTERP);
  //  myexit(1);

  // allow multiple types of rescales, and so just identity if not doing this (i.e. PRIMTOINTERP)
  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }




#elif(VARTOINTERP==PRIMTOINTERP_JONRESCALED1)

  if(dir==1){
    // optimized for pole
    if(rescaletype==0){
      scale[RHO]=pow(r,-1.5);
    }
    else if(rescaletype==1){
      scale[RHO]=pow(r,-2.7);
    }
    scale[UU]=scale[RHO]/r;
    scale[U1]=scale[RHO];
    scale[U2]=1.0;
    scale[U3]=1.0/(r*r);
    if(rescaletype==0){
      scale[B1]=scale[U3];
    }
    else if(rescaletype==1){
      scale[B1]=pow(r,-2.4);
    }
    //    if(statpos[2]+jj < 0 || startpos[2]+jj >= totalsize[2]) scale[B1] *= -1. ;
    scale[B2]=scale[B1];
    scale[B3]=scale[B1];

    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else if(dir==2){
    scale[RHO]=1.0;
    scale[UU]=1.0;
    scale[U1]=1.0;
    scale[U2]=1.0;
    scale[U3]=1.0;
    scale[B1]=1.0;
    scale[B2]=1.0;
    scale[B3]=1.0;
    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else{
    dualfprintf(fail_file,"rescale(): no such direction! dir=%d\n",dir);
    myexit(100);
  }


  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl]/scale[pl];
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl]*scale[pl];
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==CONSTOINTERP)
  // this doesn't work at all, even if no bug.
  // doesn't work even if setup better guess, as in interpU code.


  if(which==DORESCALE){ // rescale before interpolation
    MYFUN(get_state(pr, ptrgeom, &q),"interp.c:rescale()", "get_state()", 1);
    MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, p2interp, NULL),"interp.c:rescale()", "primtoU()", 1);
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    struct of_newtonstats newtonstats;   setnewtonstatsdefault(&newtonstats);
    int showmessages=1;
    int allowlocalfailurefixandnoreport=1;
    int eomtype=EOMDEFAULT;
    FTYPE dissmeasure=-1.0; // assume energy inversion ok
    int whichcap=CAPTYPEBASIC;
    int whichmethod=MODEDEFAULT;
    int modprim=0;
    int checkoninversiongas=CHECKONINVERSION;
    int checkoninversionrad=CHECKONINVERSIONRAD;
    MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,OTHERUTOPRIM,UDIAG,p2interp, NULL, ptrgeom, dissmeasure, pr, pr,&newtonstats),"interp.c:rescale()", "Utoprimgen", 1);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERPLGDEN)


  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];
    if(dopl[RHO]) p2interp[RHO]=log(pr[RHO]);
    if(dopl[UU]) p2interp[UU]=log(pr[UU]);
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];
    if(dopl[RHO]) pr[RHO]=exp(p2interp[RHO]);
    if(dopl[UU]) pr[UU]=exp(p2interp[UU]);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }

#elif(VARTOINTERP==PRIMTOINTERP_LGDEN_RHOU)
  // unstable!
  // probably because using LGDEN with RHOU.
  // Between 1E-6 and 1, density would be 1E-3 in log.
  // But v=0 and 1, so rho*v~0 and 1, but then resulting v=10^3 !

  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];
    if(dopl[RHO]) p2interp[RHO]=log(pr[RHO]);
    if(dopl[UU]) p2interp[UU]=log(pr[UU]);

    if(dopl[U1]) p2interp[U1]=pr[RHO]*pr[U1];
    if(dopl[U2]) p2interp[U2]=pr[RHO]*pr[U2];
    if(dopl[U3]) p2interp[U3]=pr[RHO]*pr[U3];
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];
    if(dopl[RHO]) pr[RHO]=exp(p2interp[RHO]);
    if(dopl[UU]) pr[UU]=exp(p2interp[UU]);

    if(dopl[U1]) pr[U1]=p2interp[U1]/pr[RHO];
    if(dopl[U2]) pr[U2]=p2interp[U2]/pr[RHO];
    if(dopl[U3]) pr[U3]=p2interp[U3]/pr[RHO];

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RHOU)
  // kinda works, except for test=7 (peak problem)

  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    if(dopl[U1]) p2interp[U1]=pr[RHO]*pr[U1];
    if(dopl[U2]) p2interp[U2]=pr[RHO]*pr[U2];
    if(dopl[U3]) p2interp[U3]=pr[RHO]*pr[U3];
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    if(dopl[U1]) pr[U1]=p2interp[U1]/pr[RHO];
    if(dopl[U2]) pr[U2]=p2interp[U2]/pr[RHO];
    if(dopl[U3]) pr[U3]=p2interp[U3]/pr[RHO];

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }




#elif(VARTOINTERP==PRIMTOINTERP_VSQ)



#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==DORESCALE){ // before interpolation, get quantities to interpolate
    if(dopl[U1]){
      quasivsq_compute(pr, ptrgeom, &quasivsq);
    }
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];
    
    if(dopl[U1]) p2interp[U1]=pr[RHO]*pr[U1];
    if(dopl[U2]) p2interp[U2]=pr[RHO]*pr[U2];
    if(dopl[U3]) p2interp[U3]=pr[RHO]*pr[U3];
    
    if(dopl[VSQ]){
      // max helps limit oscillatory behavior for non-limiter schemes
      p2interp[VSQ]=max(quasivsq,0.0);
      //p2interp[VSQ]=log(quasivsq); // assumes positive
    }
  }
  else  if(which==UNRESCALE){ // after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    if(dopl[U1]) pr[U1]=p2interp[U1]/pr[RHO];
    if(dopl[U2]) pr[U2]=p2interp[U2]/pr[RHO];
    if(dopl[U3]) pr[U3]=p2interp[U3]/pr[RHO];

    if(dopl[VSQ]){
      // now rescale velocities to agree with quasivsq
      // max helps limit oscillatory behavior for non-limiter schemes
      limit_quasivsq(max(p2interp[VSQ],0.0),ptrgeom,pr);
      //    limit_quasivsq(exp(p2interp[VSQ]),ptrgeom,pr);
    }
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_3VEL_GAMMA)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==DORESCALE){ // before interpolation, get quantities to interpolate

    if(dopl[U1] && VSQ>=0){    
      pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
      // 3-velocity
      SLOOPA(j) vcon[j]=ucon[j]/ucon[TT];
    }

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    for(pl=U1;pl<=U3;pl++) if(dopl[pl]) p2interp[pl]= vcon[pl-U1+1];

    if(dopl[VSQ] && VSQ>=0) p2interp[VSQ]=ucon[TT];
    
  }
  else  if(which==UNRESCALE){ // after interpolation

    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    for(pl=U1;pl<=U3;pl++) if(dopl[pl])  vcon[pl-U1+1] = p2interp[pl]; 

    if(dopl[U1] && VSQ>=0){
      SLOOPA(j) ucon[j] = vcon[j]*p2interp[VSQ];
      ucon2pr(WHICHVEL,ucon,ptrgeom,pr);
    }

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RHOV_GAMMA)

#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==DORESCALE){ // before interpolation, get quantities to interpolate
    if(dopl[U1] && dopl[VSQ] && VSQ>=0){    
      pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
      // 3-velocity
      SLOOPA(j) vcon[j]=ucon[j]/ucon[TT];
    }

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    //comment this out if you do not want to interpolate rho * v
    for(pl=U1;pl<=U3;pl++)     if(dopl[pl]) p2interp[pl]= pr[RHO] * vcon[pl-U1+1];

    if(dopl[VSQ]) p2interp[VSQ]=ucon[TT];

  }
  else  if(which==UNRESCALE){ // after interpolation

    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    //comment this out if you do notwant to interp. rhov
    for(pl=U1;pl<=U3;pl++)  if(dopl[pl]) vcon[pl-U1+1] = p2interp[pl]/pr[RHO]; 

    if(dopl[U1] && dopl[VSQ] && VSQ>=0){
      SLOOPA(j) ucon[j] = vcon[j]*p2interp[VSQ];
      ucon2pr(WHICHVEL,ucon,ptrgeom,pr);
    }

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERP_VELREL4SQ)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif

#if(WHICHVEL==VEL3)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if WHICHVEL=%d\n",VARTOINTERP,WHICHVEL);
  myexit(1);
#endif



  if(which==DORESCALE){ // before interpolation, get quantities to interpolate

    if(dopl[U1]){
      // get relative 4-velocity
      if(WHICHVEL!=VELREL4){
        pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
        ucon2pr(VELREL4,ucon,ptrgeom,newpr);
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=newpr[UU+j];
      }
      else{
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=pr[UU+j];
      }
      
      // get |\tilde{u}|
      lower_vec(uconrel,ptrgeom,ucovrel);
      //normuconrel=sqrt(dot(uconrel,ucovrel));
      normuconrel=dot(uconrel,ucovrel);
    }

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];
    if(dopl[VSQ]) p2interp[VSQ]=normuconrel;

  }
  else  if(which==UNRESCALE){ // after interpolation

    // assign over everything, adjust velocity below
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    if(dopl[U1]){
      // renormalize by relative 4-velocity
      if(WHICHVEL!=VELREL4){
        // get relative 4-velocity from interpolated velocity
        pr2ucon(WHICHVEL,p2interp, ptrgeom ,ucon);
        ucon2pr(VELREL4,ucon,ptrgeom,newpr);
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=newpr[UU+j];
      }
      else{
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=p2interp[UU+j];
      }

      // get |\tilde{u}| from interpolated \tilde{u}^i
      lower_vec(uconrel,ptrgeom,ucovrel);
      //    normuconrel_fromui=sqrt(dot(uconrel,ucovrel));
      normuconrel_fromui=dot(uconrel,ucovrel);

    }

    if(dopl[VSQ] && VSQ>0){
      if(WHICHVEL==VEL3){
        //      pr2ucon(WHICHVEL,p2interp,ptrgeom,ucon);
        //      ucon2pr(VELREL4,ucon,ptrgeom,newpr); // now have relative 4-velocity
        //      uconrel[TT]=0.0;
        //      SLOOPA(j) uconrel[j]=newpr[UU+j];
        // not done, but should be possible
        myexit(1);
      }
      else{ // WHICHVEL==VEL4 or WHICHVEL==VELREL4 : acts same on both
        // directly renormalize primitives
        //      SLOOPA(j) pr[UU+j] = p2interp[UU+j]*p2interp[VSQ]/absuconrel_fromui;
        SLOOPA(j) pr[UU+j] = p2interp[UU+j]*fabs(p2interp[VSQ]/normuconrel_fromui);
        // done!
      }
    }





  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_3VELREL_GAMMAREL)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif


  if(which==DORESCALE){ // before interpolation, get quantities to interpolate

    if(dopl[U1]){
      // get relative 4-velocity
      if(WHICHVEL!=VELREL4){
        pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
        ucon2pr(VELREL4,ucon,ptrgeom,newpr);
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=newpr[UU+j];
      }
      else{
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=pr[UU+j];
      }
      
      // get Lorentz factor w.r.t. relative 4-velocity
      gamma_calc_fromuconrel(uconrel,ptrgeom,&gamma,&qsq);
    }

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    if(dopl[U1]){
      // interpolate relative 3-velocity
      for(pl=U1;pl<=U3;pl++) p2interp[pl]= uconrel[pl-U1+1]/gamma;
    }

    if(dopl[VSQ] && VSQ>=0){
      // interpolate \gamma separately
      p2interp[VSQ]=gamma;
    }

  }
  else  if(which==UNRESCALE){ // after interpolation

    // assign over everything, adjust velocity below
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    if(dopl[U1] && dopl[VSQ] && VSQ>=0){
      // get relative 4-velocity from \gamma and relative 3-velocity
      uconrel[TT]=0;
      SLOOPA(j) uconrel[j]=p2interp[UU+j]*p2interp[VSQ];
      
      // get WHICHVEL velocity
      if(WHICHVEL!=VELREL4){
        pr2ucon(VELREL4,p2interp, ptrgeom ,ucon);
        ucon2pr(WHICHVEL,ucon,ptrgeom,newpr);
        SLOOPA(j) pr[UU+j]=newpr[UU+j];
      }
      else{
        SLOOPA(j) pr[UU+j]=uconrel[j];
      }
      
    }
    

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_3VELREL_GAMMAREL_DXDXP)


#if(DOEXTRAINTERP==0)
  dualfprintf(fail_file,"Shouldn't be trying to do VARTOINTERP=%d if DOEXTRAINTERP=%d\n",VARTOINTERP,DOEXTRAINTERP);
  myexit(1);
#endif

  //compute jacobian
  //dxdxprim(X, V, dxdxp);
  dxdxprim_ijk( ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, dxdxp );
  
  if(which==DORESCALE){ // before interpolation, get quantities to interpolate

    if(dopl[U1]){
      // get relative 4-velocity
      if(WHICHVEL!=VELREL4){
        pr2ucon(WHICHVEL,pr, ptrgeom ,ucon);
        ucon2pr(VELREL4,ucon,ptrgeom,newpr);
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=newpr[UU+j];
      }
      else{
        uconrel[TT]=0.0;
        SLOOPA(j) uconrel[j]=pr[UU+j];
      }
      
      // get Lorentz factor w.r.t. relative 4-velocity
      gamma_calc_fromuconrel(uconrel,ptrgeom,&gamma,&qsq);
    }

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    if(dopl[U1]){
      // interpolate relative 3-velocity
      for(pl=U1;pl<=U3;pl++) p2interp[pl]= uconrel[pl-U1+1]/gamma;
    }
    if(dopl[U2]){
      //Interpolate u^\theta and B^\theta instead of u^2 and B^2:
      p2interp[U2] *= dxdxp[2][2];
    }
    if(dopl[B2]) p2interp[B2] *= dxdxp[2][2];

    if(dopl[VSQ]&&VSQ>=0){
      // interpolate \gamma separately
      p2interp[VSQ]=gamma;
    }

  }
  else  if(which==UNRESCALE){ // after interpolation

    //return back to u^2 and B^2 from u^\theta and B^\theta
    if(dopl[U2]) p2interp[U2] /= dxdxp[2][2];
    if(dopl[B2]) p2interp[B2] /= dxdxp[2][2];
    
    // assign over everything, adjust velocity below
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    if(dopl[U1]){
      // get relative 4-velocity from \gamma and relative 3-velocity
      uconrel[TT]=0;
      SLOOPA(j) uconrel[j]=p2interp[UU+j]*p2interp[VSQ];
      
      
      // get WHICHVEL velocity
      if(WHICHVEL!=VELREL4){
        pr2ucon(VELREL4,p2interp, ptrgeom ,ucon);
        ucon2pr(WHICHVEL,ucon,ptrgeom,newpr);
        SLOOPA(j) pr[UU+j]=newpr[UU+j];
      }
      else{
        SLOOPA(j) pr[UU+j]=uconrel[j];
      }
    }

    

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#elif(VARTOINTERP==PRIMTOINTERP_RAMESH1)


  if(which==DORESCALE){ // rescale before interpolation
    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    // for the fields, do ramesh-like interpolation
    if(dopl[B1]){ pl=B1; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[GIND(1,1)])*pow(r,2.0-nu) ); }
    if(dopl[B2]){ pl=B2; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[GIND(2,2)])*pow(r,2.0-nu) ); }
    if(dopl[B3]){ pl=B3; p2interp[pl] = pr[pl] * (sqrt(ptrgeom->gcov[GIND(3,3)])*pow(r,2.0-nu) ); }
  }
  else if(which==UNRESCALE){ // unrescale after interpolation
    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    // for the fields, do ramesh-like interpolation
    if(dopl[B1]){ pl=B1; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[GIND(1,1)])*pow(r,2.0-nu) ); }
    if(dopl[B2]){ pl=B2; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[GIND(2,2)])*pow(r,2.0-nu) ); }
    if(dopl[B3]){ pl=B3; pr[pl] = p2interp[pl] / (sqrt(ptrgeom->gcov[GIND(3,3)])*pow(r,2.0-nu) ); }

  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERP_GDETFULLVERSION)


  if(which==DORESCALE){ // rescale before interpolation

    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    // if(dopl[U1]) p2interp[U1]=pr[U1]*(ptrgeom->gdet);
    // if(dopl[U2]) p2interp[U2]=pr[U2]*(ptrgeom->gdet);
    // if(dopl[U3]) p2interp[U3]=pr[U3]*(ptrgeom->gdet);

    if(dir==1||dir==3){
      if(dopl[U1]) p2interp[U1]=pr[U1]*(ptrgeom->gdet);
      if(dopl[U2]) p2interp[U2]=pr[U2]*(ptrgeom->gdet);
      if(dopl[U3]) p2interp[U3]=pr[U3]*(ptrgeom->gdet);

      if(URAD1>=0){
        if(dopl[URAD1]) p2interp[URAD1]=pr[URAD1]*(ptrgeom->gdet);
        if(dopl[URAD2]) p2interp[URAD2]=pr[URAD2]*(ptrgeom->gdet);
        if(dopl[URAD3]) p2interp[URAD3]=pr[URAD3]*(ptrgeom->gdet);
      }

      if(dopl[B1]) p2interp[B1]=pr[B1]*(ptrgeom->gdet);
      if(dopl[B2]) p2interp[B2]=pr[B2]*(ptrgeom->gdet);
      if(dopl[B3]) p2interp[B3]=pr[B3]*(ptrgeom->gdet);
    }
    if(0&&dir==2){
      //if(dopl[U1]) p2interp[U1]=pr[U1]*(ptrgeom->gdet);
      //if(dopl[U2]) p2interp[U2]=pr[U2]*(ptrgeom->gdet);
      //if(dopl[U3]) p2interp[U3]=pr[U3]*(ptrgeom->gdet);

      if(dopl[B1]) p2interp[B1]=pr[B1]*(ptrgeom->gdet); //sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=pr[B2]*(ptrgeom->gdet);//*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*pow(V[1],3);
      if(dopl[B3]) p2interp[B3]=pr[B3]*(ptrgeom->gdet);//*pow(V[1],3);
    }
  }
  else if(which==UNRESCALE){ // unrescale after interpolation

    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    //if(dopl[U1]) pr[U1]=p2interp[U1]/(ptrgeom->gdet);
    //if(dopl[U2]) pr[U2]=p2interp[U2]/(ptrgeom->gdet);
    //if(dopl[U3]) pr[U3]=p2interp[U3]/(ptrgeom->gdet);

    if(dir==1||dir==3){
      if(dopl[U1]) pr[U1]=p2interp[U1]/(ptrgeom->gdet);
      if(dopl[U2]) pr[U2]=p2interp[U2]/(ptrgeom->gdet);
      if(dopl[U3]) pr[U3]=p2interp[U3]/(ptrgeom->gdet);

      if(URAD1>=0){
        if(dopl[URAD1]) pr[URAD1]=p2interp[URAD1]/(ptrgeom->gdet);
        if(dopl[URAD2]) pr[URAD2]=p2interp[URAD2]/(ptrgeom->gdet);
        if(dopl[URAD3]) pr[URAD3]=p2interp[URAD3]/(ptrgeom->gdet);
      }

      if(dopl[B1]) pr[B1]=p2interp[B1]/(ptrgeom->gdet);
      if(dopl[B2]) pr[B2]=p2interp[B2]/(ptrgeom->gdet);
      if(dopl[B3]) pr[B3]=p2interp[B3]/(ptrgeom->gdet);
    }
    if(0&&dir==2){
      //if(dopl[U1]) pr[U1]=p2interp[U1]/(ptrgeom->gdet);
      //if(dopl[U2]) pr[U2]=p2interp[U2]/(ptrgeom->gdet);
      //if(dopl[U3]) pr[U3]=p2interp[U3]/(ptrgeom->gdet);

      if(dopl[B1]) pr[B1]=p2interp[B1]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B2]) pr[B2]=p2interp[B2]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(ptrgeom->gdet);///(pow(V[1],3));
    }
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERP==PRIMTOINTERP_GDETFULLVERSION_WALD)


  if(which==DORESCALE){ // rescale before interpolation


    PINTERPLOOP(pliter,pl) p2interp[pl]=pr[pl];

    // if(dopl[U1]) p2interp[U1]=pr[U1]*(ptrgeom->gdet);
    // if(dopl[U2]) p2interp[U2]=pr[U2]*(ptrgeom->gdet);
    // if(dopl[U3]) p2interp[U3]=pr[U3]*(ptrgeom->gdet);

    if(dir==1||dir==3){
      if(dopl[U1]) p2interp[U1]=pr[U1]*(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*(V[1]));
      if(dopl[U2]) p2interp[U2]=pr[U2]*(sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*(V[1]));
      if(dopl[U3]) p2interp[U3]=pr[U3]*(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));

      if(URAD1>=0){
        if(dopl[URAD1]) p2interp[URAD1]=pr[URAD1]*(ptrgeom->gdet);
        if(dopl[URAD2]) p2interp[URAD2]=pr[URAD2]*(ptrgeom->gdet);
        if(dopl[URAD3]) p2interp[URAD3]=pr[URAD3]*(ptrgeom->gdet);
      }

      if(dopl[B1]) p2interp[B1]=pr[B1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)]));
      if(dopl[B2]) p2interp[B2]=pr[B2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]));
      if(dopl[B3]) p2interp[B3]=pr[B3]*(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));
    }
    if(dir==2){
      if(dopl[U1]) p2interp[U1]=pr[U1]*(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*(V[1]));
      if(dopl[U2]) p2interp[U2]=pr[U2]*(sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*(V[1]));
      if(dopl[U3]) p2interp[U3]=pr[U3]*(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));

      if(dopl[B1]) p2interp[B1]=pr[B1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)]));
      if(dopl[B2]) p2interp[B2]=pr[B2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]));
      if(dopl[B3]) p2interp[B3]=pr[B3]*(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));
    }
  }
  else if(which==UNRESCALE){ // unrescale after interpolation

    PINTERPLOOP(pliter,pl) pr[pl]=p2interp[pl];

    //if(dopl[U1]) pr[U1]=p2interp[U1]/(ptrgeom->gdet);
    //if(dopl[U2]) pr[U2]=p2interp[U2]/(ptrgeom->gdet);
    //if(dopl[U3]) pr[U3]=p2interp[U3]/(ptrgeom->gdet);

    if(dir==1||dir==3){
      if(dopl[U1]) pr[U1]=p2interp[U1]/(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*(V[1]));
      if(dopl[U2]) pr[U2]=p2interp[U2]/(sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*(V[1]));
      if(dopl[U3]) pr[U3]=p2interp[U3]/(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));

      if(URAD1>=0){
        if(dopl[URAD1]) pr[URAD1]=p2interp[URAD1]/(ptrgeom->gdet);
        if(dopl[URAD2]) pr[URAD2]=p2interp[URAD2]/(ptrgeom->gdet);
        if(dopl[URAD3]) pr[URAD3]=p2interp[URAD3]/(ptrgeom->gdet);
      }

      if(dopl[B1]) pr[B1]=p2interp[B1]/sqrt(fabs(ptrgeom->gcov[GIND(1,1)]));//(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B2]) pr[B2]=p2interp[B2]/sqrt(fabs(ptrgeom->gcov[GIND(2,2)]));//(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));//(ptrgeom->gdet);///(pow(V[1],3));
    }
    if(dir==2){
      if(dopl[U1]) pr[U1]=p2interp[U1]/(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*(V[1]));
      if(dopl[U2]) pr[U2]=p2interp[U2]/(sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*(V[1]));
      if(dopl[U3]) pr[U3]=p2interp[U3]/(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));

      if(dopl[B1]) pr[B1]=p2interp[B1]/sqrt(fabs(ptrgeom->gcov[GIND(1,1)]));
      if(dopl[B2]) pr[B2]=p2interp[B2]/sqrt(fabs(ptrgeom->gcov[GIND(2,2)]));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(sqrt(fabs(ptrgeom->gcov[GIND(3,3)]))*(V[1]));
    }
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#endif // end over choices for VARTOINTERP








  //////////////////////////////////////////
  //
  // SEPARATE FIELD-ONLY TRANSFORMATIONS.  Assumes VARTOINTERP=PRIMTOINTERP so nothing done otherwise for fields.
  //
  /////////////////////////////////////////


#if(VARTOINTERPFIELD==NOFIELDRESCALE)

  // DO NOTHING, assume VARTOINTERP takes care of any transformation or non-transformation

#elif(VARTOINTERPFIELD==PULSARFIELD)


  if(dir==1){
    if(which==DORESCALE){ // rescale before interpolation

      Bconin[0]=0.0;
      Bconin[1]=pr[B1];
      Bconin[2]=pr[B2];
      Bconin[3]=pr[B3];

      dxdxprim(X, V, dxdxp);

      getconsts(Bconin, V, ptrgeom, dxdxp,Bconout);

      if(dopl[B1]) p2interp[B1]=Bconout[1];
      if(dopl[B2]) p2interp[B2]=Bconout[2];
      if(dopl[B3]) p2interp[B3]=Bconout[3];
    }
    else if(which==UNRESCALE){ // unrescale after interpolation

      Bconin[0]=0.0;
      Bconin[1]=p2interp[B1];
      Bconin[2]=p2interp[B2];
      Bconin[3]=p2interp[B3];

      dxdxprim(X, V, dxdxp);

      undoconsts(Bconin, V, ptrgeom, dxdxp, Bconout);

      if(dopl[B1]) pr[B1]=Bconout[1];
      if(dopl[B2]) pr[B2]=Bconout[2];
      if(dopl[B3]) pr[B3]=Bconout[3];
      
    }
    else{
      dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
      myexit(100);
    }
  }
  // problem when used with dir=2 since near axis end up dividing by 0


#elif(VARTOINTERPFIELD==PULSARFIELD2)


  if(dir==1){
    if(which==DORESCALE){ // rescale before interpolation

      Bconin[0]=0.0;
      Bconin[1]=pr[B1];
      Bconin[2]=pr[B2];
      Bconin[3]=pr[B3];

      dxdxprim(X, V, dxdxp);

      if(dopl[B1]) p2interp[B1]=Bconin[1]*dxdxp[1][1]*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=Bconin[2]*dxdxp[2][2]*pow(V[1],4);
      if(dopl[B3]) p2interp[B3]=Bconin[3]*pow(V[1],3);
    }
    else if(which==UNRESCALE){ // unrescale after interpolation

      dxdxprim(X, V, dxdxp);

      Bconin[0]=0.0;
      Bconin[1]=p2interp[B1]/(dxdxp[1][1]*pow(V[1],3));
      Bconin[2]=p2interp[B2]/(dxdxp[2][2]*pow(V[1],4));
      Bconin[3]=p2interp[B3]/(pow(V[1],3));

      if(dopl[B1]) pr[B1]=Bconin[1];
      if(dopl[B2]) pr[B2]=Bconin[2];
      if(dopl[B3]) pr[B3]=Bconin[3];
      
    }
    else{
      dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
      myexit(100);
    }
  }
  // problem when used with dir=2 since near axis end up dividing by 0


#elif(VARTOINTERPFIELD==PULSARFIELD3)


  if(dir==1){
    if(which==DORESCALE){ // rescale before interpolation

      Bconin[0]=0.0;
      Bconin[1]=pr[B1];
      Bconin[2]=pr[B2];
      Bconin[3]=pr[B3];

      dxdxprim(X, V, dxdxp);

      if(dopl[B1]) p2interp[B1]=Bconin[1]*sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=Bconin[2]*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*pow(V[1],3);
      if(dopl[B3]) p2interp[B3]=Bconin[3]*pow(V[1],3);
    }
    else if(which==UNRESCALE){ // unrescale after interpolation

      dxdxprim(X, V, dxdxp);

      Bconin[0]=0.0;
      Bconin[1]=p2interp[B1]/(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      Bconin[2]=p2interp[B2]/(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      Bconin[3]=p2interp[B3]/(pow(V[1],3));

      if(dopl[B1]) pr[B1]=Bconin[1];
      if(dopl[B2]) pr[B2]=Bconin[2];
      if(dopl[B3]) pr[B3]=Bconin[3];
      
    }
    else{
      dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
      myexit(100);
    }
  }
  // problem when used with dir=2 since near axis end up dividing by 0


#elif(VARTOINTERPFIELD==GDETVERSION)

  //  dxdxprim_ijk( ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, dxdxp );


  if(which==DORESCALE){ // rescale before interpolation

    Bconin[0]=0.0;
    Bconin[1]=pr[B1];
    Bconin[2]=pr[B2];
    Bconin[3]=pr[B3];

    if(dir==1){
      if(dopl[B1]) p2interp[B1]=Bconin[1]*(ptrgeom->gdet); //sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=Bconin[2]*(ptrgeom->gdet);//*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*pow(V[1],3);
      if(dopl[B3]) p2interp[B3]=Bconin[3]*(ptrgeom->gdet);//*pow(V[1],3);
    }
    else if(dir==2){
      if(dopl[B1]) p2interp[B1]=Bconin[1];
      if(dopl[B2]) p2interp[B2]=Bconin[2];
      if(dopl[B3]) p2interp[B3]=Bconin[3];
    }
    else{
      if(dopl[B1]) p2interp[B1]=Bconin[1];
      if(dopl[B2]) p2interp[B2]=Bconin[2];
      if(dopl[B3]) p2interp[B3]=Bconin[3];
    }
  }
  else if(which==UNRESCALE){ // unrescale after interpolation

    Bconin[0]=0.0;
    Bconin[1]=p2interp[B1];
    Bconin[2]=p2interp[B2];
    Bconin[3]=p2interp[B3];
    
    if(dir==1){
      if(dopl[B1]) pr[B1]=p2interp[B1]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B2]) pr[B2]=p2interp[B2]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(ptrgeom->gdet);///(pow(V[1],3));
    }
    else if(dir==2){
      if(dopl[B1]) pr[B1]=Bconin[1];
      if(dopl[B2]) pr[B2]=Bconin[2];
      if(dopl[B3]) pr[B3]=Bconin[3];
    }
    else{
      if(dopl[B1]) pr[B1]=Bconin[1];
      if(dopl[B2]) pr[B2]=Bconin[2];
      if(dopl[B3]) pr[B3]=Bconin[3];
    }
      
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#elif(VARTOINTERPFIELD==GDETFULLVERSION)


  if(which==DORESCALE){ // rescale before interpolation

    if(dir==1||dir==3){
      if(dopl[B1]) p2interp[B1]=pr[B1]*(ptrgeom->gdet); //sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=pr[B2]*(ptrgeom->gdet);//*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*pow(V[1],3);
      if(dopl[B3]) p2interp[B3]=pr[B3]*(ptrgeom->gdet);//*pow(V[1],3);
    }
    if(0&&dir==2){
      if(dopl[B1]) p2interp[B1]=pr[B1]*(ptrgeom->gdet); //sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3);
      if(dopl[B2]) p2interp[B2]=pr[B2]*(ptrgeom->gdet);//*sqrt(fabs(ptrgeom->gcov[GIND(2,2)]))*pow(V[1],3);
      if(dopl[B3]) p2interp[B3]=pr[B3]*(ptrgeom->gdet);//*pow(V[1],3);
    }
  }
  else if(which==UNRESCALE){ // unrescale after interpolation

    if(dir==1||dir==3){
      if(dopl[B1]) pr[B1]=p2interp[B1]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B2]) pr[B2]=p2interp[B2]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(ptrgeom->gdet);///(pow(V[1],3));
    }
    if(0&&dir==2){
      if(dopl[B1]) pr[B1]=p2interp[B1]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B2]) pr[B2]=p2interp[B2]/(ptrgeom->gdet);///(sqrt(fabs(ptrgeom->gcov[GIND(1,1)]))*pow(V[1],3));
      if(dopl[B3]) pr[B3]=p2interp[B3]/(ptrgeom->gdet);///(pow(V[1],3));
    }
      
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }


#endif // end over choices for VARTOINTERPFIELD



  return(0);
}

