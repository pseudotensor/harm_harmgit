
/*! \file metric.c
     \brief Sets metric in some coordinates based upon metric.h type choice for MCOORD
// this file includes metric dependent terms, including for initial
// condition routines for IC coords.

// crucially need to setup analytic form of gcov.  All rest can be done numerically or analytically if wanted.
//  SUPERNOTE: Set DOMIXTHETAPHI to 0 or 1 in definit.h!!
// self-gravity TODO:
// http://www.fftw.org/
// COSMOS++ uses: https://computation.llnl.gov/casc/hypre/software.html
// http://ciera.northwestern.edu/StarCrash/manual/html/node3.html
// http://ccfd-jacob.blogspot.com
// http://farside.ph.utexas.edu/teaching/329/lectures/node60.html
*/

#include "decs.h"




/// obtain gcov in primcoords of whichcoord type metric/coords
/// here ptrgeom is only expected to contain i,j,k,p location
void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc)
{
  void set_gcov_cylminkmetric    (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_spcminkmetric    (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_cartminkmetric   (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_unigravity   (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_htmetric         (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_htmetric_accurate(FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_ksmetric         (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_ks_jp1_metric         (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_ks_bh_tov_metric (FTYPE *X, FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  extern void set_gcov_ks_tov_metric(FTYPE *X, FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  extern void set_gcov_bl_tov_metric(FTYPE *X, FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void set_gcov_blmetric         (FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
  void gcov2gcovprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc, FTYPE *gcovinfuncprim, FTYPE *gcovpertinfuncprim);
  FTYPE gcovselfpert[NDIM];
  extern int set_gcov_selfspcmetric(FTYPE *X, FTYPE *V, FTYPE *gcovselfpert);
  FTYPE V[NDIM],Vmetric[NDIM],Xmetric[NDIM];
  int j,k;
  int presenttime;
  FTYPE phi;
  LOCALMETRICTEMPVARS;






  ////////////////////////////////////
  //
  // determine if asking for metric now or in past
  //
  /////////////////////////////////////
  if(DOEVOLVEMETRIC){

    //    if(ptrgeom->i==7 && nstep==1084){
    //      dualfprintf(fail_file,"X[0]=%21.15g Xmetricold[0]=%21.15g Xmetricnew[0]=%21.15g t=%21.15g p=%d\n",X[0],Xmetricold[0],Xmetricnew[0],t,ptrgeom->p);
    //    }

    if(X[0]>=Xmetricnew[0]-SMALL){ // SMALL trying to account for numerical roundoff error. Assume never want metric at a time difference SMALL before
      // by present we mean last time latest metric was computed
      presenttime=1;
    }
    else{
      // should only ever possibly get here if DOEVOLVEMETRIC==1
      presenttime=0;
#if(DOEVOLVEMETRIC!=1)
      dualfprintf(fail_file,"presenttime=0 and DOEVOLVEMETRIC==0 are incompatible: X[0]=%21.15g t=%21.15g\n",X[0],t);
      myexit(5523);
#endif
      if(ptrgeom->p==NOWHERE){
        dualfprintf(fail_file,"Not capable of obtaining past metric that wasn't stored\n");
        myexit(2356);
        // for example, can't just interpolate metric since then won't necessarily satisfy divg=0
        // although probably not big error for anything requiring this interpolation
        // that is, Connection and most evolution things will use metric at known location
        // so could infact interpolate if this is needed -- wait and see
      }
      if(getprim==0){
        // not designed to feed back anything not in PRIMECOORDS when asking for past metric
        dualfprintf(fail_file,"getprim==0 is not compatible with wanting past metric\n");
        myexit(6246);
      }
    }
  }
  else presenttime=1; // no evolution, so always present




  ///////////////////////////
  //
  // Determine current full X, get V, then get g_{\mu\nu}
  //
  ///////////////////////////
  if(presenttime){


    // if PRIMECOORDS, did store the data before, and choosing a gridded position, then can just grab from memory
    if(whichcoord==PRIMECOORDS && didstoremetricdata==1 && ptrgeom->p!=NOWHERE){

      // get local metric quantities for this loc,i,j,k
      GETLOCALMETRIC(ptrgeom->p,ptrgeom->i,ptrgeom->j,ptrgeom->k);

      DLOOP(j,k){
        gcovinfunc[GIND(j,k)]=localgcov[GIND(j,k)];
      }
      DLOOPA(j){
        gcovpertinfunc[j]=localgcovpert[j];
      }
      
    }
    else{ // then not defined yet or arbitrary position


      // before here, X set by i,j,k or chosen directly
      bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, X, V);
      //bl_coord(X,V);


      //dualfprintf(fail_file,"blcoordcalled: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);


      ////////////////////////////
      //
      // Use V[X] and get Vmetric[V[X]] since all set_gcov's are in terms of Vmetric=r,h,ph (i.e. axisymmetric-type old/original versions)
      //
      ////////////////////////////
      rotate_VtoVmetric(whichcoord,THETAROT,V,Vmetric);
      int jj; DLOOPA(jj) Xmetric[jj]=X[jj]; // no change?  Depends upon how Xmetric used in set_gcov. GODMARK
      

 
      ////////////////////////////
      //
      // Use Vmetric[] and get metric
      //
      ////////////////////////////
      if(whichcoord>=0){
        if(whichcoord==BLCOORDS){
          set_gcov_blmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==KSCOORDS){
          set_gcov_ksmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==KSCARTCOORDS){
          //          set_gcov_kscartmetric(Vmetric, gcovinfunc, gcovpertinfunc); // can just call wrapper around set_gcov_ksmetric()
        }
        else if(whichcoord==KS_JP1_COORDS){
          set_gcov_ks_jp1_metric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==KS_BH_TOV_COORDS){
          set_gcov_ks_bh_tov_metric(Xmetric, Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==KS_TOV_COORDS){
          set_gcov_ks_tov_metric(Xmetric, Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==BL_TOV_COORDS){
          set_gcov_bl_tov_metric(Xmetric, Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==HTMETRIC){
          set_gcov_htmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==HTMETRICACCURATE){
          set_gcov_htmetric_accurate(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==CARTMINKMETRIC || whichcoord==CARTMINKMETRIC2){
          set_gcov_cartminkmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==UNIGRAVITY){
          set_gcov_unigravity(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==CYLMINKMETRIC){
          set_gcov_cylminkmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else if(whichcoord==SPCMINKMETRIC){
          set_gcov_spcminkmetric(Vmetric, gcovinfunc, gcovpertinfunc);
        }
        else{
          dualfprintf(fail_file,"gcov_func(): no such whichcoord=%d\n",whichcoord);
          myexit(12626);
        }

        //      dualfprintf(fail_file,"blcoordcalled2: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);

      }
      else{
        dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcov from prim gcov): getprim=%d whichcoord=%d\n",getprim,whichcoord);
        myexit(26632);
      }

      /////////////////////////
      //
      // Account for self-gravity term.  Add it to g_{\mu\nu}.
      //
      /////////////////////////
#if(DOSELFGRAVVSR)
      if(
         whichcoord==BLCOORDS ||
         whichcoord==KSCOORDS ||
         whichcoord==KS_JP1_COORDS ||
         whichcoord==HTMETRIC ||
         whichcoord==HTMETRICACCURATE ||
         whichcoord==SPCMINKMETRIC
         ){
        // then using spherical polar coordinates and can do self-gravity as designed
        set_gcov_selfspcmetric(Xmetric,Vmetric,gcovselfpert);

        // OLD DEBUG
        //    phi = -0.1/Vmetric[1];
        //gcovinfunc[GIND(TT,TT)] += GINDASSIGNFACTOR(TT,TT)*(-2.0*phi);
        //      gcovpertinfunc[TT] += -2.0*phi;
        //      gcovinfunc[GIND(RR,RR)] += GINDASSIGNFACTOR(RR,RR)*(-2.0*phi);
        //      gcovpertinfunc[RR] += -2.0*phi;

        //    gcovinfunc[GIND(TT,TT)] += GINDASSIGNFACTOR(TT,TT)*gcovselfpert[TT];
        //    gcovpertinfunc[TT] += gcovselfpert[RR];
        //    gcovinfunc[GIND(RR,RR)] += GINDASSIGNFACTOR(RR,RR)*gcovselfpert[TH];
        //    gcovpertinfunc[RR] += gcovselfpert[PH];

        //    if(gcovselfpert[TT]+2.0*phi>0.01){
        // dualfprintf(fail_file,"got TT difference: %21.15g %21.15g\n",gcovselfpert[TT],-2.0*phi);
        //    }

        //    if(gcovselfpert[RR]+2.0*phi>0.01){
        // dualfprintf(fail_file,"got RR difference: %21.15g %21.15g\n",gcovselfpert[RR],-2.0*phi);
        //    }

        //    if(gcovselfpert[TH]>0.01){
        // dualfprintf(fail_file,"got TH difference: %21.15g %21.15g\n",gcovselfpert[TH],0.0);
        //    }

        //    if(gcovselfpert[PH]>0.01){
        // dualfprintf(fail_file,"got PH difference: %21.15g %21.15g\n",gcovselfpert[PH],0.0);
        //    }


        //      DLOOPA(j){
        // dualfprintf(fail_file,"t=%21.15g %ld %d :: X1=%21.15g :: postgcov[%d][%d]=%21.15g :: gcovselfpert[%d]=%21.15g mypert=%2.15g\n",t,steppart,nstep,Xmetric[1],j,j,gcovinfunc[GIND(j,j)],j,gcovselfpert[j],-2.0*phi);
        //      }

        // DLOOPA(j){
        //   dualfprintf(fail_file,"t=%21.15g %ld %d :: X1=%21.15g :: gcovselfpert[%d]=%21.15g\n",t,steppart,nstep,Xmetric[1],j,gcovselfpert[j]);
        // }


        ////////////
        //
        // add self-gravity perturbation to metric
        //
        ////////////
        DLOOPA(j){
          gcovinfunc[GIND(j,j)] += GINDASSIGNFACTOR(j,j)*gcovselfpert[j];
          gcovpertinfunc[j]+=gcovselfpert[j]; // in this way the perturbation part is always resolved
          // dualfprintf(fail_file,"gcovselfpert[%d]=%21.15g\n",j,gcovselfpert[j]);
        }
        if(whichcoord==KSCOORDS || whichcoord==KS_JP1_COORDS){
          // KS-form has these terms
          gcovinfunc[GIND(TT,RR)] += GINDASSIGNFACTOR(TT,RR)*gcovselfpert[TT];
          gcovinfunc[GIND(RR,TT)] = gcovinfunc[GIND(TT,RR)];
        }
      }
      else if(whichcoord==KS_BH_TOV_COORDS || whichcoord==KS_TOV_COORDS || whichcoord==BL_TOV_COORDS){
        // then doing full TOV-like solution that doesn't just come in as a perturbation
      
      }
      else{
        // then not setup for these coordinates
        dualfprintf(fail_file,"DOSELFGRAVVSR=1 will not work with this whichcoord=%d\n",whichcoord);
        myexit(145);
      }
#endif


      //  DLOOP(j,k) { stderrfprintf("1gcov[%d][%d]=%21.15g\n",j,k,gcovinfunc[GIND(j,k)]); fflush(stderr);}



      ///////////
      //
      // If rotating metric, now we perform transformation on the metric differentials
      //
      // X[] inputted is really X, not Xmetric or anything else.
      // i.e., X[] is associated with V, not Vmetric.
      ///////////
      transVmetrictoV(whichcoord,ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, THETAROT, X, V, Xmetric, Vmetric, gcovinfunc,gcovpertinfunc);


      ///////////
      //
      // Convert to prim coords
      //
      ///////////
      if(getprim==1){
        // all the above are analytic, so have to convert to prim coords.
        gcov2gcovprim(ptrgeom, Xmetric, Vmetric, gcovinfunc,gcovpertinfunc, gcovinfunc, gcovpertinfunc);
      }
      //  DLOOP(j,k) { stderrfprintf("2gcov[%d][%d]=%21.15g\n",j,k,gcovinfunc[GIND(j,k)]); fflush(stderr);}

      //    if(ptrgeom->i==7 && nstep==1084){
      //    //      DLOOP(j,k) dualfprintf(fail_file,"present time: gcov[%d][%d]=%21.15g\n",j,k,gcovinfunc[GIND(j,k)]);
      // }
    }// end if needing to recompute metric


  }
  else{

    //////////////
    //
    // then not present time, and assume past time metric stored, so recover
    // also assume never ask for past time at arbitrary location, checked previously
    //
    //////////////

    if(ptrgeom->p!=NOWHERE){

      // get local metric quantities for this loc,i,j,k
      GETLASTLOCALMETRIC(ptrgeom->p,ptrgeom->i,ptrgeom->j,ptrgeom->k);

      DLOOP(j,k){
        gcovinfunc[GIND(j,k)]=localgcov[GIND(j,k)];
      }
      DLOOPA(j){
        gcovpertinfunc[j]=localgcovpert[j];
      }

      //      if(ptrgeom->i==7 && nstep==1084){
      // DLOOP(j,k) dualfprintf(fail_file,"past time: gcov[%d][%d]=%21.15g\n",j,k,gcovinfunc[GIND(j,k)]);
      //      }
    }
    else{
      // then don't have grid at right location and need to interpolate
#if(NEWMETRICSTORAGE)
      interpX_gcov(Xmetric, GLOBALPOINT(compgeomlast), NULL, NULL, gcovinfunc, gcovpertinfunc);
#else
      interpX_gcov(Xmetric, NULL, GLOBALPOINT(gcovlast), GLOBALPOINT(gcovpertlast), gcovinfunc, gcovpertinfunc);
#endif
      // above untested, so fail for now
      dualfprintf(fail_file,"interpX_gcov() unexpectedly called\n");
      myexit(13515);
    } 
  }// end else if doing past time

  //dualfprintf(fail_file,"blcoordcalled3: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);



}


///Note that mathematica's Arctan[x,y] = C's atan2(y,x) (i.e. args are flipped in order)
FTYPE arctanmath(FTYPE x, FTYPE y)
{
  return(atan2(y,x));
}


/////////////////
//// Rotate V[X] (r,thnew,phnew as in bl_coord) to Vmetric (r,th,ph as in set_gcov)
/// That is, when using rotated metric, we assume metric itself still has Vmetric=r,th,ph, while X is mapped to V=rnew,hnew,phnew.
/// So input to set_gcov(X,V) needs to then get V->Vmetric before getting gcov(Vmetric).
///  So real spherical polar coordinates for grid itself is V=rnew,hnew,phnew.
///
/// Use this rather than direct transformation(rnew,hnew,phnew) because that expression becomes too complicated in mathematica.  So just take 3 steps:
/// 1) gcov_func(X,V)
/// 2) rotate_VtoVmetric(V,Vmetric)
/// 3) set_gcov...(Vmetric)
/// 4) transVmetrictoV(gcov) (which internally uses rotate_VtoVmetric() to keep expressions simple (i.e. kept in terms of Vmetric)
///
/// #2 just accounts for Vmetric[V[X]] as far as assignment of metric values so that X[] grid has used correct values of old/metric r,h,ph
/// #4 accounts for differentials in metric so that ds^2 is the same.  This is written in terms of V=rnew,hnew,phnew, so can just feed-in V[X].
///
////////////////
int rotate_VtoVmetric(int whichcoord, FTYPE ROTANGLE, FTYPE *V, FTYPE *Vmetric)
{

  if(ROTANGLE!=0.0 && ALLOWMETRICROT && ISSPCMCOORD(whichcoord)){
    // see metricrot.nb.  
    FTYPE tnew,rnew,hnew,phnew; // true V used in bl_coord() to map to X
 
    tnew=V[0];
    rnew=V[1];
    hnew=V[2];
    phnew=V[3];


    FTYPE told,rold,hold,phold; // what is inside metric functions like set_gcov()
    told=tnew;
    rold=rnew;

    FTYPE b0=ROTANGLE;

    hold=arctanmath (cos(b0)*cos(hnew) - 1.*cos(phnew)*sin(b0)*sin(hnew),pow(pow(cos(hnew)*sin(b0) + cos(b0)*cos(phnew)*sin(hnew),2) + pow(sin(hnew),2)*pow(sin(phnew),2),0.5));

    phold=arctanmath (cos(hnew)*sin(b0) + cos(b0)*cos(phnew)*sin(hnew),sin(hnew)*sin(phnew));


    fix_hp(&hold,&phold);

    Vmetric[0]=told;
    Vmetric[1]=rold;
    Vmetric[2]=hold;
    Vmetric[3]=phold;

  }
  else{
    int jj; DLOOPA(jj) Vmetric[jj]=V[jj];
  }

  
  return(0);

}


/// put \theta,\phi in positive normal region of \theta,\phi.
int fix_hp(FTYPE *h, FTYPE *p)
{
  FTYPE th=*h,ph=*p;

#if(0)
  // keep \theta between 0 and \pi
  if(th<0.0) th+=M_PI;
  if(th>=M_PI) th-=M_PI;

  // keep \phi between 0 and 2\pi
  if(ph<0.0) ph+=2.0*M_PI;
  if(ph>=2.0*M_PI) ph-=2.0*M_PI;
#else

  // keep \phi between 0 and 2\pi.  Can always do that full rotation.
  // assume never more out of phase that 1 full rotation
  if(ph<0.0) ph+=2.0*M_PI;
  if(ph>=2.0*M_PI) ph-=2.0*M_PI;

  // keep \theta between 0 and \pi and \phi between 0 and 2\pi
  // but need to be at same physical SPC location, not arbitrary rotation
  if(th>=0.0 && th<=M_PI){
    // do nothing
  }
  else if(th<0.0){
    th*=-1.0;
    if(ph<=M_PI) ph+=M_PI;
    else if(ph>M_PI) ph-=M_PI;
    else{ dualfprintf(fail_file,"Shouldn't be here1 with th=%g ph=%g\n",th,ph); myexit(1);}
  }
  else if(th>=M_PI){
    th=M_PI-th;
    if(ph<=M_PI) ph+=M_PI;
    else if(ph>M_PI) ph-=M_PI;
    else{ dualfprintf(fail_file,"Shouldn't be here2 with th=%g ph=%g\n",th,ph); myexit(1);}
  }
  else{ dualfprintf(fail_file,"Shouldn't be here3 with th=%g ph=%g\n",th,ph); myexit(1);}
#endif

  *h=th;
  *p=ph;

  return(0);

}

/////////////////////////////
///
/// perform rotation of V (for spherical polar coordinates)
///
/// This takes input of V=(t,r,h,ph) and gives back Vmetric=(t,r,hnew,phnew)
/// where h,ph are in original metric form and hnew,phnew are rotated versions
///
/// So since we want X->V to be mapping V=(t,r,hnew,phnew), this is *not* to be used.
///
/// -----
///
/// This can be used in (e.g.) __init__.py to have python script take data (in Vnew=V) and obtain Vmetric version
///
/// 1) transVtoVmetric(gcovnew) gives gcov[original metric]
/// 2) transVtoVmetric(ucon,bcon,ucov,bcov) or transVmetrictoV(ucon,bcon,ucov,bcov)
/// 3) Rotate actual spatial positions of data, including metrics, so that again axisymmetric so only have to store 1 phi slice!
///
/////////////////////////////
int rotate_VmetrictoV(int whichcoord, FTYPE ROTANGLE, FTYPE *Vmetric, FTYPE *V)
{

  if(ROTANGLE!=0.0 && ALLOWMETRICROT && ISSPCMCOORD(whichcoord)){
    // see metricrot.nb.  Note that mathematica's Arctan[x,y] = C's atan2(y,x) (i.e. args are flipped in order)
    FTYPE t,r,h,ph;
 
    t=Vmetric[0];
    r=Vmetric[1];
    h=Vmetric[2];
    ph=Vmetric[3];

    FTYPE xc,yc,zc,Rold;
    xc=mysin(h)*mycos(ph);
    yc=mysin(h)*mysin(ph);
    zc=mycos(h);
    Rold=sqrt(xc*xc + yc*yc);

    // rotation around y-axis using right-hand rule
    FTYPE xcnew,ycnew,zcnew,Rnew;
    xcnew=xc*mycos(ROTANGLE)-zc*mysin(ROTANGLE);
    ycnew=yc;
    zcnew=zc*mycos(ROTANGLE)+xc*mysin(ROTANGLE);
    Rnew=sqrt(xcnew*xcnew + ycnew*ycnew);

    // Below uses atan2 so gets back correct quadrant
    FTYPE tnew,rnew,hnew,phnew;
    tnew=t;
    rnew=r;
    // these return in range [-\pi,\pi]
    hnew=atan2(Rnew,zcnew);
    phnew=atan2(ycnew,xcnew);


    fix_hp(&hnew,&phnew);


    V[0]=tnew;
    V[1]=rnew;
    V[2]=hnew;
    V[3]=phnew;

  }
  else{
    int jj; DLOOPA(jj) V[jj]=Vmetric[jj];
  }



  return(0);
}



/// interpolate grid-based gcov to arbitrary location
/// used for stored past grid
int interpX_gcov(FTYPE *X, struct of_compgeom PTRDEFMETMACP1A0(compgeom,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3), FTYPE PTRDEFMETMACP1A2(gcovgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM), FTYPE PTRDEFMETMACP1A1(gcovpertgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM), FTYPE *gcov, FTYPE *gcovpert)
{
  int i,j,k;
  int jj,kk;
  extern void icoord(FTYPE *X,int loc, int *i, int *j, int *k);
  int ip,jp,kp;
  FTYPE Xijk[NDIM],Xip[NDIM],Xjp[NDIM],Xkp[NDIM];
  FTYPE gijk,gip,gjp,gkp;
  FTYPE dist[4+1],totaldist;
  int loc;


  // CENT is chosen as reference location for all interpolations
  loc=CENT;

  // use icoord() to get i,j,k from X so can get potential at any location from interpolation of existing potential at certain locations

  // get centered i,j,k
  icoord(X,loc,&i,&j,&k); // truncates to lower i  
  coord_ijk(i,j,k, loc, Xijk);
  // only want effective grid location (i.e. normalized X position) so that differences are 1 between grid cells -- don't care about absolute magnitude of X once normalized
  DLOOPA(jj) Xijk[jj]/=dx[jj];

  /////////////////
  // X1 dir
  //
  // limit lower value of i
  if(i<=-N1BND-SHIFT1){
    i=-N1BND;
  }
  // limit upper value of i
  if(i>=N1+N1BND-SHIFT1){ // -1 is to offset so +1 is in range
    i=N1+N1BND-SHIFT1-SHIFT1; // first -1 is to offset from out of range, and second -1 is so ip is located within array range
  }

  // get X-position of coordinates i and i+1
  ip=i+SHIFT1;
  coord_ijk(ip,j,k, loc, Xip);
  DLOOPA(jj) Xip[jj]/=dx[jj];



  /////////////////
  // X2 dir
  //
  // limit lower value of j
  if(j<=-N2BND-SHIFT2){
    j=-N2BND;
  }
  // limit upper value of j
  if(j>=N2+N2BND-SHIFT2){ // -1 is to offset so +1 is in range
    j=N2+N2BND-SHIFT2-SHIFT2; // first -1 is to offset from out of range, and second -1 is so jp is located within array range
  }

  // get X-position of coordinates j and j+1
  jp=j+SHIFT2;
  coord_ijk(i,jp,k, loc, Xjp);
  DLOOPA(jj) Xjp[jj]/=dx[jj];



  /////////////////
  // X3 dir
  //
  // limit lower value of k
  if(k<=-N3BND-SHIFT3){
    k=-N3BND;
  }
  // limit upper value of k
  if(k>=N3+N3BND-SHIFT3){ // -1 is to offset so +1 is in range
    k=N3+N3BND-SHIFT3-SHIFT3; // first -1 is to offset from out of range, and second -1 is so kp is located within array range
  }

  // get X-position of coordinates k and k+1
  kp=k+SHIFT3;
  coord_ijk(i,j,kp, loc, Xkp);
  DLOOPA(jj) Xkp[jj]/=dx[jj];


  // now use some interpolation of 3 points

  // use this position to interpolate in X (which happens to be the same as interpolating in ijk)
  // distance away from one point toward ijk
  dist[1]=(1-(X[1]-Xijk[1]))*(1-(X[2]-Xijk[2]))*(1-(X[3]-Xijk[3]));
  // away from one point toward ip
  dist[2]=(1-(X[1]-Xip[1]))*(1-(X[2]-Xip[2]))*(1-(X[3]-Xip[3]))*(1.0-dist[1]);
  // away from one point toward jp
  dist[3]=(1-(X[1]-Xjp[1]))*(1-(X[2]-Xjp[2]))*(1-(X[3]-Xjp[3]))*(1.0-dist[2])*(1.0-dist[1]);
  // away from one point toward kp
  dist[4]=(1-(X[1]-Xkp[1]))*(1-(X[2]-Xkp[2]))*(1-(X[3]-Xkp[3]))*(1.0-dist[3])*(1.0-dist[2])*(1.0-dist[1]);
  // normalization of total distance
  totaldist=dist[1]+dist[2]+dist[3]+dist[4];
  


  DLOOP(jj,kk){
#if(NEWMETRICSTORAGE)
    gijk=METMACP1A0(compgeom,loc,i,j,k).gcov[GIND(jj,kk)];
    gip=METMACP1A0(compgeom,loc,ip,j,k).gcov[GIND(jj,kk)];
    gjp=METMACP1A0(compgeom,loc,i,jp,k).gcov[GIND(jj,kk)];
    gkp=METMACP1A0(compgeom,loc,i,j,kp).gcov[GIND(jj,kk)];
#else
    gijk=METMACP1A2(gcovgrid,loc,i,j,k,jj,kk);
    gip=METMACP1A2(gcovgrid,loc,ip,j,k,jj,kk);
    gjp=METMACP1A2(gcovgrid,loc,i,jp,k,jj,kk);
    gkp=METMACP1A2(gcovgrid,loc,i,j,kp,jj,kk);
#endif

    // X's are per unit dX's
    gcov[GIND(jj,kk)] = (gijk*dist[1]+gip*dist[2]+gjp*dist[3]+gkp*dist[4])/totaldist;

  }

  DLOOPA(jj){
#if(NEWMETRICSTORAGE)
    gijk=METMACP1A0(compgeom,loc,i,j,k).gcovpert[jj];
    gip=METMACP1A0(compgeom,loc,ip,j,k).gcovpert[jj];
    gjp=METMACP1A0(compgeom,loc,i,jp,k).gcovpert[jj];
    gkp=METMACP1A0(compgeom,loc,i,j,kp).gcovpert[jj];
#else
    gijk=METMACP1A1(gcovpertgrid,loc,i,j,k,jj);
    gip=METMACP1A1(gcovpertgrid,loc,ip,j,k,jj);
    gjp=METMACP1A1(gcovpertgrid,loc,i,jp,k,jj);
    gkp=METMACP1A1(gcovpertgrid,loc,i,j,kp,jj);
#endif
    // X's are per unit dX's
    gcovpert[jj] = (gijk*dist[1]+gip*dist[2]+gjp*dist[3]+gkp*dist[4])/totaldist;

  }

  return(0);

}



/// obtain prim gcon in primcoords of whichcoord type metric/coords
void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon)
{
  void set_gcon_blmetric(FTYPE *V, FTYPE *gcon);
  void set_gcon_ksmetric(FTYPE *V, FTYPE *gcon);
  void gcon2gconprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE *gcon,FTYPE *gconprim);
  void matrix_inverse_metric(int whichcoord, FTYPE *gcov, FTYPE *gcon);
  int j,k;
  FTYPE V[NDIM];


  if(whichcoord>=0){
    if(whichcoord==BLCOORDS){
      //dualfprintf(fail_file,"mi in BLCOORDS\n");
      if(ANALYTICGCON){
        bl_coord(X, V);
        set_gcon_blmetric(V, gcon);
        if(getprim) gcon2gconprim(ptrgeom, X, V, gcon,gcon);
      }
      // since don't have gcon and want to keep things simple by only having to specify gcov 
      else matrix_inverse_metric(whichcoord,gcov,gcon);

    }
    else if(whichcoord==KSCOORDS){
      //dualfprintf(fail_file,"mi in KSCOORDS\n");
      //DLOOP(j,k) dualfprintf(fail_file,"ks gcov[%d][%d]=%g\n",j,k,gcov[GIND(j,k)]);

      // DEBUG:
      //bl_coord(X, V);
      //dualfprintf(fail_file,"i=%d j=%d k=%d loc=%d :: %21.15g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,V[1]);


      if(ANALYTICGCON){
        bl_coord(X, V);
        set_gcon_ksmetric(V,gcon);
        if(getprim) gcon2gconprim(ptrgeom, X, V,gcon,gcon);
      }
      // since don't have gcon and want to keep things simple by only having to specify gcov 
      else matrix_inverse_metric(whichcoord,gcov,gcon);
    }
    else if(whichcoord==KS_JP1_COORDS || whichcoord==HTMETRIC || whichcoord==HTMETRICACCURATE || whichcoord==CARTMINKMETRIC || whichcoord==CARTMINKMETRIC2 || whichcoord==UNIGRAVITY || whichcoord==CYLMINKMETRIC || whichcoord==SPCMINKMETRIC || whichcoord==KS_BH_TOV_COORDS || whichcoord==KS_TOV_COORDS || whichcoord==BL_TOV_COORDS || whichcoord==KSCARTCOORDS){
      // do not have analytic gcon, so invert numerically
      matrix_inverse_metric(whichcoord,gcov,gcon);
    }
    else{
      dualfprintf(fail_file,"gcon_func(): no such whichcoord=%d\n",whichcoord);
      myexit(2687);
    }
  }
  else{
    dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcon from prim gcon\n");
    myexit(3783);
  }



}









/// connection not simply transformed -- so compute directly from final metric (primcoords)
void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
               FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  void set_conn_general(FTYPE *X, struct of_geom *geom,
                        FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
                         FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
                              FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
                               FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  

  if(whichcoord==KS_JP1_COORDS || whichcoord==BLCOORDS || whichcoord==KS_BH_TOV_COORDS || whichcoord==KS_TOV_COORDS || whichcoord==BL_TOV_COORDS || whichcoord==HTMETRIC || whichcoord==HTMETRICACCURATE || whichcoord==UNIGRAVITY || whichcoord==SPCMINKMETRIC || whichcoord==KSCARTCOORDS){
    set_conn_general(X,geom,conn,conn2);
  }
  else if(whichcoord==KSCOORDS){
    set_conn_ksmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==CYLMINKMETRIC){
    set_conn_cylminkmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==CARTMINKMETRIC || whichcoord==CARTMINKMETRIC2){
    set_conn_cartminkmetric(X,geom,conn,conn2);
  }
  else{
    dualfprintf(fail_file,"conn_func(): no such whichcoord=%d\n",whichcoord);
    myexit(7367);
  }
}



/// obtain eomfunc (f_{(\nu)} factor: see how used in connection below and get_geometry())
/// assumes X set before eomfun_func()
/// should allow for DOEVOLVEMETRIC.  Currently WITHGDET and WITHNOGDET are usable with DOEVOLVEMETRIC
/// do ALL-pl (EOMFUNCMAC(pl))
void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *EOMFUNCNAME)
{

  FTYPE gcovmcoord[SYMMATRIXNDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE r,th;
  int j,k;
  int pl,pliter;
  FTYPE V[NDIM];
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcovpert);
  FTYPE ftemp;


  if(WHICHEOM!=WITHGDET || NEWMETRICSTORAGE==1){
    // then do eomfunc_func() calculation
  }
  else{
    // don't
    // ensures avoid memory leak even if users mistakenly calls this function
    dualfprintf(fail_file,"Don't call eomfunc_func() whichcoord=%d\n",whichcoord);
    return;
  }




  if(WHICHEOM==WITHGDET){
    gcov_func(ptrgeom, getprim, whichcoord,X,gcovmcoord,gcovpertcoord); // actually returns primcoords version of whichcoord
    bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, X, V);
    if(gdet_func_metric(whichcoord,V,gcovmcoord,&ftemp)!=0){
      if(debugfail>=2) dualfprintf(fail_file,"Caught gdet_func_metric() problem in eomfunc_func()\n");
    }
    PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=ftemp;
  }
  else if(WHICHEOM==WITHNOGDET){
    PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=1.0;
  }
  else if(WHICHEOM==WITHSINSQ){ // obviously coordinate dependent
    bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, X, V);

    // rotate consistently with rotation in metric.c
    // that is, anytime take i,j,k -> V for metric stuff, need to first rotate V.
    r=V[1]; th=V[2]; // NOTEMARK :No rotation_V here because th and sin(th)^2 below are for polar grid itself, not the metric.
    ftemp=sin(th)*sin(th);
    PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=ftemp;
  }
  else{
    dualfprintf(fail_file,"eomfunc_func(): no such WHICHEOM=%d\n",WHICHEOM);
    myexit(798436);
  }
  
  if(FORCEGDETPOSITIVE==1){
    PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=fabs(EOMFUNCASSIGN(pl));
  }

}



////////////////////////////
///
/// g_{\mu\nu} (always analytic -- returns some coordinate system)
///
////////////////////////////
/// needs M, J, and Q
/// M: total M+dM mass of star
/// J: total angular momentum
/// Q: mass quadrapole moment
/// external space-time of slowly rotating star, accurate to second order in $\Omega_{\star}$
/// c=G=1
/// NS mass sheds at v=0.8c (eq 28), allowed approx up to v=0.1c
/// WD mass sheds at v=0.01c, allowed approx up to v=0.001c
/// equally restrictive as measured by ratio of rotational energy to gravitational energy
void set_gcov_htmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE P2,Q12, Q22;
  FTYPE z,OO,RO,SS,AA,alphasq,betaphi;
  FTYPE SSpert,OOpert,AApert;
  FTYPE M,J,Q;
  FTYPE r,th;
  FTYPE ftemp;
  // size of r relative to M is controlled by M, not r.
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];


#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

#define HTMETRICTYPE FULLHT

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M
  // NOW: uses MBH as mass of NS (or any compact object within the boundary)
  M=MBH;

  if(HTMETRICTYPE==FULLHT){
    J=a*MBH; // I\Omega = J
    Q=a*a*MBH; // approximately, perhaps, like Kerr geometry, if hard enough EOS
  }
  else if(HTMETRICTYPE==NOSPIN){
    a = 0.0;
    J=a*MBH; // no spin of metric
    Q=a*a*MBH;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=MBH=1E-100; // no mass of star, just flat space (GJ model)
    a=0.0;
    J=a*MBH;
    Q=a*a*MBH;
  }


  // surface radius is not constant radius if Q!=0
  // if Q=0, then surface radius is 3M?  depends on EOS.

  // let's set M=1, such that J is up to M^2 and Q is up to M^3 such that radius is in units of GM/c^2

  // For Sun: M/r<M/R~2E-6, J/r^2<J/R^2~1E-12 (J/M^2=4E-24), Q/r^3<Q/R^3<~1E-10 (Q/R^3=8E-28) and surface roughly spherical

  // Fastest pulsar PSR 1937+21 \tau=1.56ms , 640times/sec
  // Deepto Chakrabarty of MIT
  // Vela: PSR 0833-45 \tau=89.3ms
  // strongest: PSR 0329+54 \tau=0.715s
  // peak: 760times/sec (limit by GR?), 2-3X is theoretical limit but would mass shed.
  // formation spin rate: 30/sec


  // see Hartle & Thorne 1968 or Kim Lee Lee Lee 2005

  // Legendre polynomial
  P2=0.5*(3.0*cos(th)*cos(th)-1.0);

  z=rsharp/M-1.0;

  // Mathematica's LegendreQ (associated Legendre of 2nd kind) is such that:
  //   Q12=Im[Legendre[2,1,z]] and Q22=-Re[Legendre[2,2,z]]
  Q12=sqrt(z*z-1.0)*( (3.0*z*z-2.0)/(z*z-1.0)-3.0/2.0*z*log((z+1)/(z-1)) );

  Q22=(3.0/2.0*(z*z-1.0)*log((z+1)/(z-1))-(3.0*z*z*z-5.0*z)/(z*z-1.0) );



  // metric stuff
  OOpert=-2.0*M/r+2.0*J*J/(r*r*r*r);
  OO=1.0+OOpert;

  RO=1.0+2.0*(J*J/(M*r*r*r)*(1.0+M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  
  SSpert=-2.0*(J*J/(M*r*r*r)*(1.0-5.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  SS=1.0+SSpert;

  AApert=+2.0*(-J*J/(M*r*r*r)*(1.0+2.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*( (2.0*M/(sqrt(r*r*(1.0-2.0*M/r))))*Q12 - Q22  ))*P2;
  AA=1.0+AApert;
  
  alphasq=OO*RO;

  betaphi=-2.0*J/(r*r*r); // note that omega_{zamo}=-betaphi

  gcov[GIND(PH,PH)] = rsharp*rsharp*AA*sin(th)*sin(th) ;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;
  
  gcov[GIND(TT,TT)] = -(alphasq-betaphi*betaphi*gcov[GIND(PH,PH)]) ;
  ftemp = (1.0-alphasq);
  // g_{tt} + 1
  gcovpert[TT] = ftemp + betaphi*betaphi*gcov[GIND(PH,PH)] ;

  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = betaphi*gcov[GIND(PH,PH)] ;
    
  gcov[GIND(RR,TT)] = 0.0 ;
  gcov[GIND(RR,RR)] = SS/OO ;
  // SS/OO = (SSpert+1.0)/OO = SSpert/OO + 1.0/OO , so SS/OO-1.0 = SSpert/OO + (1.0/OO-1.0)
  ftemp=(1.0/OO - 1.0);
  gcovpert[RR] = SSpert/OO + ftemp ; // might help with radial gravity terms

  gcov[GIND(RR,TH)] = 0.0 ;
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = rsharp*rsharp*AA ;
  // r^2 AA - 1  = r^2(AApert+1.0) -1 = r^2 AApert + r^2 -1, so not a big issue with catastrophic cancellation in non-rel limit
  // probably have to define background as being r^2 to have gravity work at large radii (GODMARK)
  gcovpert[TH] = gcov[GIND(TH,TH)] - 1.0;

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;

  
             



}

#undef HTMETRICTYPE

/// from Berti, White, Maniopoulou, and Bruni (2005)
void set_gcov_htmetric_accurate(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE M,J,Q;
  FTYPE u,p,A1,A2,W,F1,F2,H1,H2,L,G1;
  FTYPE r,th;
  FTYPE ftemp;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];



#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

#define HTMETRICTYPE FULLHT

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M
  // NOW: uses MBH as mass of NS (or any compact object within the boundary)
  M=MBH;

  if(HTMETRICTYPE==FULLHT){
    J=a*MBH; // I\Omega = J
    Q=a*a*MBH; // approximately, perhaps, like Kerr geometry, if hard enough EOS
  }
  else if(HTMETRICTYPE==NOSPIN){
    a = 0.0;
    J=a*MBH; // no spin of metric
    Q=a*a*MBH;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=MBH=1E-100; // no mass of star, just flat space (GJ model)
    a=0.0;
    J=a*MBH;
    Q=a*a*MBH;
  }





  L=(80.0*pow(M,6)+8.0*pow(M,4)*r*r+10.0*M*M*M*r*r*r+20.0*M*M*pow(r,4)-45.0*M*pow(r,5)+15.0*pow(r,6));

  u=cos(th);
  p=1.0/(8.0*M*pow(r,4)*(r-2.0*M)); // typo in paper for parenthesis


  A1=(15.0*r*(r-2.0*M)*(1.0-3.0*u*u))/(16.0*M*M)*log(r/(r-2.0*M));
  A2=(15.0*(r*r-2.0*M*M)*(3.0*u*u-1.0))/(16.0*M*M)*log(r/(r-2.0*M));
  
  W=(r-M)*(16.0*pow(M,5)+8.0*pow(M,4)*r-10*M*M*r*r*r-30.0*M*pow(r,4)+15.0*pow(r,5))+u*u*(48.0*pow(M,6)-8.0*pow(M,5)*r-24.0*pow(M,4)*r*r-30.0*M*M*M*r*r*r-60.0*M*M*pow(r,4)+135.0*M*pow(r,5)-45.0*pow(r,6));

  F1=-p*W+A1;
  F2=5.0*r*r*r*p*(3.0*u*u-1.0)*(r-M)*(2.0*M*M+6.0*M*r-3.0*r*r)-A1;

  H1=A2+(1.0/(8.0*M*r*r*r*r))*(1.0-3.0*u*u) * (16.0*pow(M,5)+8.0*pow(M,4)*r-10.0*M*M*r*r*r+15.0*M*pow(r,4)+15.0*pow(r,5)) ;
  H2=-A2 + (1.0/(8.0*M*r))*5.0*(1.0-3.0*u*u)*(2.0*M*M-3.0*M*r-3.0*r*r);

  G1=p*((L-72.0*M*M*M*M*M*r)-3.0*u*u*(L-56.0*M*M*M*M*M*r))-A1;
  




  ftemp  = J*J*F1-Q*F2;
  gcov[GIND(TT,TT)] = -(1.0-2.0*M/r)*(1.0+ftemp);
  // g_{tt} + 1
  //  gcovpert[TT] =  -(J*J*F1-Q*F2) + (2.0*M/r)*(1.0+J*J*F1-Q*F2);
  gcovpert[TT] =  -ftemp +(1.0+ftemp)*(2.0*M/r);



  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = (2.0*J*M*M/r)*sin(th)*sin(th);
    
  gcov[GIND(RR,TT)] = 0.0 ;

  ftemp=J*J*G1+Q*F2;
  gcov[GIND(RR,RR)] = (1.0/(1.0-2.0*M/r))*(1.0+ftemp);
  // g_{rr} - 1
  gcovpert[RR] = (ftemp+2.0*M/r)/(1.0-2.0*M/r); // should help with radial gravity in non-rel regime

  gcov[GIND(RR,TH)] = 0.0 ;
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = rsharp*rsharp*(1.0+J*J*H1-Q*H2) ;
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;
  gcov[GIND(PH,PH)] = gcov[GIND(TH,TH)]*sin(th)*sin(th);
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;

  
             



}
#undef HTMETRICTYPE



/// Cartesian minkowski
/// (t,x,y,z)
void set_gcov_cartminkmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  
  gcov[GIND(TT,TT)] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = 0.0 ;
    
  gcov[GIND(RR,TT)] = 0.0 ;
  gcov[GIND(RR,RR)] = 1.0 ;
  gcovpert[RR] = 0.0;

  gcov[GIND(RR,TH)] = 0.0 ;
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = 1.0 ;
  gcovpert[TH] = 0.0;

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;
  gcov[GIND(PH,PH)] = 1.0 ;
  gcovpert[PH] = 0.0;

}


/// (t,x,y,z)
void set_gcov_unigravity(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{

  FTYPE phi;
  FTYPE x,y,z;
  FTYPE FORCEMAG;

  x=V[1];
  y=V[2];
  z=V[3];

  //this is the hard-cored value for test 28, right? SASMARK
  FORCEMAG=0.1/(coordparams.timescalefactor*coordparams.timescalefactor); // strength of force, scales like v^2

  // uniform force field in +y direction (i.e. dphi/dy = FORCEMAG, so Force = -FORCEMAG)
  phi = FORCEMAG*y;

  gcov[GIND(TT,TT)] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = 0.0 ;
    
  gcov[GIND(RR,TT)] = 0.0 ;
  gcov[GIND(RR,RR)] = 1.0-2.0*phi ;
  gcovpert[RR] = -2.0*phi;

  gcov[GIND(RR,TH)] = 0.0 ;
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = 1.0-2.0*phi ;
  gcovpert[TH] = -2.0*phi;

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;

  gcov[GIND(PH,PH)] = 1.0-2.0*phi ;
  gcovpert[PH] = -2.0*phi;

  //  dualfprintf(fail_file,"tsf=%21.15g :: %21.15g %21.15g :: %21.15g %21.15g %21.15g %21.15g\n",coordparams.timescalefactor,FORCEMAG,phi,gcovpert[TT],gcovpert[RR],gcovpert[TH],gcovpert[PH]);



}



/// (t,R,z,\phi)
void set_gcov_cylminkmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{

  FTYPE r,Mass,RSTAR;
  FTYPE phi;
  FTYPE R,z;



  R = V[1];
  z = V[2];




  Mass=MBH; // NOW: use MBH as any compact object mass (in length units)
  RSTAR=0.1; // could use dxdxp[RR][RR]*dx[RR] to get size
  r=sqrt(R*R+z*z);
  // These SMALL's assume not further squaring or otherwise making smaller to yield numerical 0
  // These SMALL's also assume eventually apply gdet as 0 at axis.
  FTYPE rsq=SMALL+r*r;
  FTYPE Rsq=SMALL+R*R;
  FTYPE zsq=SMALL+z*z;

  if(r<RSTAR){
    phi = (Mass/(2.0*RSTAR))*((r/RSTAR)*(r/RSTAR)-3.0);
  }
  else phi = -Mass/r;

  
  gcov[GIND(TT,TT)] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = 0.0 ;
    
  gcov[GIND(RR,TT)] = 0.0 ;
  gcov[GIND(RR,RR)] = 1.0-2.0*phi*Rsq/(rsq); 
  gcovpert[RR] = -2.0*phi*Rsq/(rsq); 

  gcov[GIND(RR,TH)] = -2.0*phi*R*z/(rsq) ; // order v^4
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = 1.0 -2.0*phi*zsq/(rsq) ;
  gcovpert[TH] =-2.0*phi*zsq/(rsq) ; // order v^4

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;
  gcov[GIND(PH,PH)] = Rsq;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0; // doesn't matter that -1.0 subtracted

  //  dualfprintf(fail_file,"Mass=%g r=%g R=%g z=%g phi=%g gtt=%g\n",Mass,r,R,z,phi,gcov[GIND(TT,TT)]);
  //  int jj,kk;
  //  DLOOP(jj,kk) dualfprintf(fail_file,"jj=%d kk=%d g=%g\n",jj,kk,gcov[GIND(jj,kk)]);

}






/// (t,r,\theta,\phi)
void set_gcov_spcminkmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE r,th;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif




  th=V[2];
  
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }

  // avoid singularity at polar axis
#if(COORDSINGFIX) // just for metric components // hack
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif


  
  gcov[GIND(TT,TT)] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[GIND(TT,RR)] = 0.0 ;
  gcov[GIND(TT,TH)] = 0.0 ;
  gcov[GIND(TT,PH)] = 0.0 ;
    
  gcov[GIND(RR,TT)] = 0.0 ;
  gcov[GIND(RR,RR)] = 1.0; 
  gcovpert[RR] = 0.0;

  gcov[GIND(RR,TH)] = 0.0 ;
  gcov[GIND(RR,PH)] = 0.0 ;
    
  gcov[GIND(TH,TT)] = gcov[GIND(TT,TH)] ;
  gcov[GIND(TH,RR)] = gcov[GIND(RR,TH)] ;
  gcov[GIND(TH,TH)] = rsharp*rsharp ;
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;

  gcov[GIND(TH,PH)] = 0.0 ;
    
  gcov[GIND(PH,TT)] = gcov[GIND(TT,PH)] ;
  gcov[GIND(PH,RR)] = gcov[GIND(RR,PH)] ;
  gcov[GIND(PH,TH)] = gcov[GIND(TH,PH)] ;

  gcov[GIND(PH,PH)] = (rsharp*sin(th))*(rsharp*sin(th));
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;

}


/// (~t,r,\theta,~\phi)
/// mixed KS BH+TOV metric
/// See grmhd-ksksp.nb, grmhd-connectiononly.nb, tov_solution_timedepsol.nb, tov_solution_timeindepsol.nb, ks_from_bl.nb
void set_gcov_ks_bh_tov_metric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);
  FTYPE MSQ;
  FTYPE phi;
  int j,k;
  void set_gcov_ksmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert);
  extern int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  FTYPE gcovtovks[SYMMATRIXNDIM],gcovtovkspert[NDIM];
  FTYPE gcovbhks[SYMMATRIXNDIM],gcovbhkspert[NDIM];
  FTYPE MtotOr,MBHprime,MBHprimeOr;
  SFTYPE MOrself,phiself,vrsqself;
  FTYPE starpot, totalpot;
  FTYPE r2small;
  FTYPE disc;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;



#if(0)
  // debug test
  DLOOP(j,k) gcov[GIND(j,k)] = 0.0;
  DLOOPA(j,j) gcov[GIND(j,j)] = 1.0;
  gcovpert[TT] gcov[GIND(TT,TT)]=-1.0;
  DLOOPA(j)  gcovpert[j]= 0.0;
#endif
  




  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];


  // get TOV ks metric assuming pre-existing phivsr_tot, MvsrOr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_ks_tov_spcmetric(X, V, gcovtovks, gcovtovkspert, &MOrself, &phiself, &vrsqself);

  // get black hole KS metric
  set_gcov_ksmetric(V, gcovbhks, gcovbhkspert);

  // default is KS BH metric
  DLOOP(j,k) gcov[GIND(j,k)] = gcovbhks[GIND(j,k)];
  DLOOPA(j) gcovpert[j] = gcovbhkspert[j];

  // now compute TOV modifications
  //  MBHprime = 2.0*MBH*r*r/(a*a+2.0*r*r+a*a*cos(2.0*th));
  MBHprime = MBH/(1.0 + (a*a + a*a*cos(2.0*th))/(2.0*r2small));
  MBHprimeOr = MBHprime/r;
  MtotOr = MBHprimeOr + MOrself; // assume MOrself has sign built-in
  starpot=vrsqself - 2.0*MOrself;
  totalpot = -2.0*MtotOr;

  // assign mixed metric components
  gcov[GIND(TT,TT)] = gcovbhks[GIND(TT,TT)]*(-gcovtovks[GIND(TT,TT)]);
  gcovpert[TT] = gcov[GIND(TT,TT)] + 1.0; // no obvious way to remove perturbed part

  //  gcov[GIND(RR,RR)] = 1.0 + 2.0*MtotOr-vrsqself;
  gcov[GIND(RR,RR)] = 1.0 + fabs(- totalpot);
  gcovpert[RR] = fabs(-totalpot);

  disc=-gcovtovks[GIND(TT,TT)]/(1.0-fabs(starpot));
  if(disc<0.0){
    disc=0.0; // assume this value never used
    // and can't trust that signature of metric will work out now, so artificially force signature to be negative
    gcov[GIND(TT,TT)]=-fabs(gcov[GIND(TT,TT)]); // should never be used
  }

  gcov[GIND(TT,RR)] = gcov[GIND(RR,TT)] = (-totalpot)*sqrt(disc);

  // DEBUG (overwrite)
  //  DLOOP(j,k) gcov[GIND(j,k)] = gcovtovks[GIND(j,k)];
  //  DLOOPA(j) gcovpert[j] = gcovtovkspert[j];

  //  DLOOP(j,k) dualfprintf(fail_file,"gcovtovks[%d][%d]=%21.15g : %21.15g %21.15g %21.15g\n",j,k,gcovtovks[GIND(j,k)],MOrself,phiself,vrsqself);


}


/// TOV KS metric wrapper
void set_gcov_ks_tov_metric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  extern int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  SFTYPE MOrself,phiself,vrsqself;


  // get TOV ks metric assuming pre-existing phivsr_tot, Mvsr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_ks_tov_spcmetric(X, V, gcov, gcovpert, &MOrself, &phiself, &vrsqself);

}

/// TOV KS metric wrapper
void set_gcov_bl_tov_metric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  extern int set_gcov_bl_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE *gcov, FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  SFTYPE MOrself,phiself,vrsqself;


  // get TOV BL metric assuming pre-existing phivsr_tot, Mvsr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_bl_tov_spcmetric(X, V, gcov, gcovpert, &MOrself, &phiself, &vrsqself);

}





/// (~t,r,\theta,~\phi)
void set_gcov_ksmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3;
  FTYPE rho2,rho2small;
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);
  FTYPE MSQ;
  FTYPE phi;
  int j,k;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;


  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
  //  dualfprintf(fail_file,"got here r=%21.15g drsing=%21.15g\n",r,drsing);
#else
  r2small = r*r + SMALL;
#endif


  //  dualfprintf(fail_file,"r=%21.15g\n",r);

  // theta gives well-defined metric for all th
  th=V[2];


  sth=mysin(th);
  cth=mycos(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = rsharp * rsharp;
  r3 = r2 * r;
  rho2 = r2 + a * a * cth * cth;
  rho2small = r2small + a * a * cth * cth;
  MSQ=(2.*MBH*r-QBH*QBH);

  // really gcov00 diverges at r=0, but force non-divergence to avoid nan's only
#define ks_gcov00 (-1. + MSQ/rho2small)
#define ks_gcov01 (MSQ/rho2small)
#define ks_gcov02 (0)
#define ks_gcov03 (-MSQ*a*s2/rho2small)
#define ks_gcov10 (ks_gcov01)
#define ks_gcov11 (1. + MSQ/rho2small)
#define ks_gcov12 (0)
#define ks_gcov13 (-a*s2*(1. + MSQ/rho2small))
#define ks_gcov20 (0)
#define ks_gcov21 (0)
#define ks_gcov22 (rho2)
#define ks_gcov23 (0)
#define ks_gcov30 (ks_gcov03)
#define ks_gcov31 (ks_gcov13)
#define ks_gcov32 (0)
  // old
  //#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2small)))
  // new
#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + MSQ/rho2small)))
  // Living review format
  //#define ks_gcov33 (s2*(r*r+a*a + a*a*s2*(MSQ/rho2small)))


  gcov[GIND(TT,TT)] = ks_gcov00 ;
  gcovpert[TT] = MSQ/rho2small;


  gcov[GIND(TT,RR)] = ks_gcov01 ;
  gcov[GIND(TT,TH)] = ks_gcov02 ;
  gcov[GIND(TT,PH)] = ks_gcov03 ;
    
  gcov[GIND(RR,TT)] = ks_gcov10 ;

  gcov[GIND(RR,RR)] = ks_gcov11 ;
  gcovpert[RR] = MSQ/rho2small ;

  gcov[GIND(RR,TH)] = ks_gcov12 ;
  gcov[GIND(RR,PH)] = ks_gcov13 ;
    
  gcov[GIND(TH,TT)] = ks_gcov20 ;
  gcov[GIND(TH,RR)] = ks_gcov21 ;

  gcov[GIND(TH,TH)] = ks_gcov22 ;
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;

  gcov[GIND(TH,PH)] = ks_gcov23 ;
    
  gcov[GIND(PH,TT)] = ks_gcov30 ;
  gcov[GIND(PH,RR)] = ks_gcov31 ;
  gcov[GIND(PH,TH)] = ks_gcov32 ;

  gcov[GIND(PH,PH)] = ks_gcov33 ;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;


  //  if(nstep>1051){
  // dualfprintf(fail_file,"MBH=%21.15g a=%21.15g gcov[GIND(TT,TT)]=%21.15g\n",MBH,a,gcov[GIND(TT,TT)]);
    
  //  DLOOP(j,k) dualfprintf(fail_file,"ks gcov[%d][%d]=%g\n",j,k,gcov[GIND(j,k)]);
  //}

}


/// Johannsen & Psaltis 2011 metric
void set_gcov_ks_jp1_metric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);
  FTYPE MSQ;
  FTYPE phi;
  int j,k;
  void set_gcov_ksmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert);
  FTYPE gcovksjp1pert[SYMMATRIXNDIM];
  FTYPE gcovbhks[SYMMATRIXNDIM],gcovbhkspert[NDIM];
  FTYPE r2small;
  FTYPE disc;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;


#if(0)
  // debug test
  DLOOP(j,k) gcov[GIND(j,k)] = 0.0;
  DLOOP(j,k) gcovksjp1pert[GIND(j,k)] = 0.0;
  DLOOPA(j,j) gcov[GIND(j,j)] = 1.0;
  gcovpert[TT] gcov[GIND(TT,TT)]=-1.0;
  DLOOPA(j)  gcovpert[j]= 0.0;
#endif
  

  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];


  // Difference between Johannsen & Psaltis 2011 metric (/data/jon/mathematica/mathematica_merged/timj_kerrschild_form2.nb)
  // See also checkinghowtimjmetricisdifferent.nb for printing of the below
  // only non-zero and non-repeated (due to symmetry) terms:
  DLOOP(j,k) gcovksjp1pert[GIND(j,k)] = 0.0;

  gcovksjp1pert[GIND(0,0)]=-4.*EP3*r*(2.*r*(-2.*MBH + r) + pow(a,2) + cos(2.*th)*pow(a,2))*pow(MBH,3)*pow(pow(a,2) + cos(2.*th)*pow(a,2) + 2.*pow(r,2),-3);

  gcovksjp1pert[GIND(0,1)]=2.*EP3*pow(MBH,4)*pow(r,2)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-3);

  gcovksjp1pert[GIND(1,1)]=-1. - 4.*pow(MBH,2)*pow(r,2)*(EP3*r*pow(MBH,3) + pow(r,4) + 2.*pow(a,2)*pow(r,2)*pow(cos(th),2) + pow(a,4)*pow(cos(th),4))*pow(r*(-2.*MBH + r) + pow(a,2),-1)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-3) - 2.*MBH*r*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-1) + (EP3*r*pow(MBH,3)*(pow(r,2) + pow(a,2)*pow(cos(th),2)) + pow(pow(r,2) + pow(a,2)*pow(cos(th),2),3))*pow((r*(-2.*MBH + r) + pow(a,2))*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2) + EP3*r*pow(a,2)*pow(MBH,3)*pow(sin(th),2),-1) + pow(a,2)*pow(r*(-2.*MBH + r) + pow(a,2),-2)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-3)*pow(sin(th),2)*(-4.*EP3*pow(MBH,5)*pow(r,3) - 4.*pow(MBH,2)*pow(r,2)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2) + (pow(a,2) + pow(r,2))*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),3) + MBH*r*pow(a,2)*(2.*EP3*r*pow(MBH,3) + EP3*pow(MBH,2)*(pow(r,2) + pow(a,2)*pow(cos(th),2)) + 2.*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),2))*pow(sin(th),2));

  gcovksjp1pert[GIND(0,3)]=-2.*a*EP3*pow(MBH,4)*pow(r,2)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-3)*pow(sin(th),2);

  gcovksjp1pert[GIND(1,3)]=0.125*a*EP3*r*pow(MBH,3)*(-4.*r*(2.*MBH + r)*pow(a,2) + 4.*r*(2.*MBH + r)*cos(2.*th)*pow(a,2) - 1.*pow(a,4) + cos(4.*th)*pow(a,4) + 32.*pow(MBH,2)*pow(r,2))*pow(r*(-2.*MBH + r) + pow(a,2),-1)*pow(pow(r,2) + pow(a,2)*pow(cos(th),2),-3)*pow(sin(th),2);

  gcovksjp1pert[GIND(3,3)]=4.*EP3*r*pow(a,2)*(2.*r*(2.*MBH + r) + pow(a,2) + cos(2.*th)*pow(a,2))*pow(MBH,3)*pow(pow(a,2) + cos(2.*th)*pow(a,2) + 2.*pow(r,2),-3)*pow(sin(th),4);



  // get black hole KS metric
  set_gcov_ksmetric(V, gcovbhks, gcovbhkspert);


  // default is KS BH metric
  // and add "perturbation" that is JP1 metric
  // NOTEMARK: when using DLOOP(j,k) over GIND(j,k) for metric, *cannot* use += since that would go over same value twice!
  DLOOP(j,k) gcov[GIND(j,k)] = gcovbhks[GIND(j,k)] + gcovksjp1pert[GIND(j,k)];
  DLOOPA(j) gcovpert[j] = gcovbhkspert[j] + gcovksjp1pert[GIND(j,j)];
  // i.e. JP1 terms are deviation from Kerr.  NOTEMARK: For g_{rr}, this adds and subtracts unity, so not acurate machine accurate.  But that's ok since applications with this metric are always for near the BH.



}





/// (t,r,\theta,\phi)
void set_gcov_blmetric(FTYPE *V, FTYPE *gcov, FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3, DD, mu;
  FTYPE mupert,DDpert;
  FTYPE r,th;
  FTYPE MSQ;
  int j,k;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif



  th=V[2];




  cth = cos(th);
  sth=sin(th);


  MSQ=(2.*MBH*rsharp-QBH*QBH);

  s2 = sth * sth;
  a2 = a * a;
  r2 = rsharp * rsharp;
  r3 = r2 * r;

  DDpert =- MSQ / r2small + a2 / r2small; 
  DD = 1. + DDpert;
  mupert=+ a2 * cth * cth / r2small;
  mu = 1. + mupert;



#define bl_gcov00 (-(1. - MSQ/(r2small*mu)))
#define bl_gcov01 (0)
#define bl_gcov02 (0)
#define bl_gcov03 (-MSQ*a*s2/(r2small*mu))
#define bl_gcov10 (0)
#define bl_gcov11 (mu/DD)
#define bl_gcov12 (0)
#define bl_gcov13 (0)
#define bl_gcov20 (0)
#define bl_gcov21 (0)
#define bl_gcov22 (r2*mu)
#define bl_gcov23 (0)
#define bl_gcov30 (bl_gcov03)
#define bl_gcov31 (0)
#define bl_gcov32 (0)
#define bl_gcov33 (r2*sth*sth*(1. + a2/r2small + MSQ*a2*s2/(r2small*r2small*mu)))

  gcov[GIND(TT,TT)] = bl_gcov00 ;
  gcovpert[TT] = MSQ/(r2small*mu);

  gcov[GIND(TT,RR)] = bl_gcov01 ;
  gcov[GIND(TT,TH)] = bl_gcov02 ;
  gcov[GIND(TT,PH)] = bl_gcov03 ;
    
  gcov[GIND(RR,TT)] = bl_gcov10 ;

  gcov[GIND(RR,RR)] = bl_gcov11 ;
  // mu/DD = (1.0+mupert)/(1.0+DDpert) = 
  gcovpert[RR] = (DDpert - mupert) / (1.0 + DDpert);

  gcov[GIND(RR,TH)] = bl_gcov12 ;
  gcov[GIND(RR,PH)] = bl_gcov13 ;
    
  gcov[GIND(TH,TT)] = bl_gcov20 ;
  gcov[GIND(TH,RR)] = bl_gcov21 ;

  gcov[GIND(TH,TH)] = bl_gcov22 ;
  gcovpert[TH] = gcov[GIND(TH,TH)]-1.0;

  gcov[GIND(TH,PH)] = bl_gcov23 ;
    
  gcov[GIND(PH,TT)] = bl_gcov30 ;
  gcov[GIND(PH,RR)] = bl_gcov31 ;
  gcov[GIND(PH,TH)] = bl_gcov32 ;

  gcov[GIND(PH,PH)] = bl_gcov33 ;
  gcovpert[PH] = gcov[GIND(PH,PH)]-1.0;



  //  DLOOP(j,k) dualfprintf(fail_file,"bl gcov[%d][%d]=%g\n",j,k,gcov[GIND(j,k)]);



}


////////////////////////////
///
/// g^{\mu\nu} (analytic or numerical)
///
////////////////////////////
/// (~t,r,\theta,~\phi)
void set_gcon_ksmetric(FTYPE *V, FTYPE *gcon)
{
  FTYPE sth, cth, s2, a2, r2, r2small,r3;
  FTYPE rho2;
  FTYPE r,th;
  FTYPE MSQ;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif
  // size of r and M controlled by changing M, not rescaling r




  th=V[2];

  cth = cos(th);
  sth=sin(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  rho2 = r2small + a * a * cth * cth; // always divide by this
  MSQ=(2.*MBH*r-QBH*QBH);


#define ks_gcon00 (-(1.+ MSQ/rho2))
#define ks_gcon01 (MSQ/rho2)
#define ks_gcon02 (0)
#define ks_gcon03 (0)
#define ks_gcon10 (ks_gcon01)
  //#define ks_gcon11 ((r*(r-2.)+a*a)/rho2)
#define ks_gcon11 ((r*r+a*a-MSQ)/rho2)
#define ks_gcon12 (0)
#define ks_gcon13 (a/rho2)
#define ks_gcon20 (ks_gcon02)
#define ks_gcon21 (ks_gcon12)
#define ks_gcon22 (1./rho2)
#define ks_gcon23 (0)
#define ks_gcon30 (ks_gcon03)
#define ks_gcon31 (ks_gcon13)
#define ks_gcon32 (ks_gcon23)
#define ks_gcon33 (1./(rho2*s2))


  gcon[GIND(TT,TT)] = ks_gcon00 ;
  gcon[GIND(TT,RR)] = ks_gcon01 ;
  gcon[GIND(TT,TH)] = ks_gcon02 ;
  gcon[GIND(TT,PH)] = ks_gcon03 ;
    
  gcon[GIND(RR,TT)] = ks_gcon10 ;
  gcon[GIND(RR,RR)] = ks_gcon11 ;
  gcon[GIND(RR,TH)] = ks_gcon12 ;
  gcon[GIND(RR,PH)] = ks_gcon13 ;
    
  gcon[GIND(TH,TT)] = ks_gcon20 ;
  gcon[GIND(TH,RR)] = ks_gcon21 ;
  gcon[GIND(TH,TH)] = ks_gcon22 ;
  gcon[GIND(TH,PH)] = ks_gcon23 ;
    
  gcon[GIND(PH,TT)] = ks_gcon30 ;
  gcon[GIND(PH,RR)] = ks_gcon31 ;
  gcon[GIND(PH,TH)] = ks_gcon32 ;
  if(s2!=0.0){
    gcon[GIND(PH,PH)] = ks_gcon33 ;
  }
  else gcon[GIND(PH,PH)]=1.0; // avoid coordinate singularity -- although should never use this value

}

/// (t,r,\theta,\phi)
void set_gcon_blmetric(FTYPE *V, FTYPE *gcon)
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3, DD, mu;
  FTYPE r,th;
  FTYPE MSQ;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  // size of r relative to M is controlled by M, not r.

  th=V[2];

  cth = cos(th);
  sth=sin(th);

  MSQ=(2.*MBH*r-QBH*QBH);

  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - MSQ / r2 + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;


#define bl_gcon00 (-1. - MSQ*(1. + a2/r2small)/(r2small*DD*mu))
#define bl_gcon01 (0)
#define bl_gcon02 (0)
#define bl_gcon03 (-MSQ*a/(r2small*r2small*DD*mu))
#define bl_gcon10 (0)
#define bl_gcon11 (DD/mu)
#define bl_gcon12 (0)
#define bl_gcon13 (0)
#define bl_gcon20 (0)
#define bl_gcon21 (0)
#define bl_gcon22 (1./(r2small*mu))
#define bl_gcon23 (0)
#define bl_gcon30 (bl_gcon03)
#define bl_gcon31 (0)
#define bl_gcon32 (0)
#define bl_gcon33 ((1. - MSQ/(r2small*mu))/(r2small*sth*sth*DD))


  gcon[GIND(TT,TT)] = bl_gcon00 ;
  gcon[GIND(TT,RR)] = bl_gcon01 ;
  gcon[GIND(TT,TH)] = bl_gcon02 ;
  gcon[GIND(TT,PH)] = bl_gcon03 ;
    
  gcon[GIND(RR,TT)] = bl_gcon10 ;
  gcon[GIND(RR,RR)] = bl_gcon11 ;
  gcon[GIND(RR,TH)] = bl_gcon12 ;
  gcon[GIND(RR,PH)] = bl_gcon13 ;
    
  gcon[GIND(TH,TT)] = bl_gcon20 ;
  gcon[GIND(TH,RR)] = bl_gcon21 ;
  gcon[GIND(TH,TH)] = bl_gcon22 ;
  gcon[GIND(TH,PH)] = bl_gcon23 ;
    
  gcon[GIND(PH,TT)] = bl_gcon30 ;
  gcon[GIND(PH,RR)] = bl_gcon31 ;
  gcon[GIND(PH,TH)] = bl_gcon32 ;
  gcon[GIND(PH,PH)] = bl_gcon33 ;

}

////////////////////////////
///
/// CONNECTIONS (numerical or analytic)
///
////////////////////////////
void set_conn_general(FTYPE *X, struct of_geom *geom,
                      FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
                            FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}



void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
                            FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
                            FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcovpert);
  int i,j,k;
  FTYPE gcovmid[SYMMATRIXNDIM];
  FTYPE gcovpertmid[NDIM];
  FTYPE gdetmid;

  FTYPE dxdxp[NDIM][NDIM]; //atch
  FTYPE V[NDIM];           //atch

  if(ANALYTICCONNECTION&&defcoord==UNIFORMCOORDS){ // uniform grid  SUPERSASMARK
    // could directly use gdet in global memory
    // only works for X1=R and X2=z
    gcov_func(geom, 1,CYLMINKMETRIC,X, gcovmid,gcovpertmid);

    bl_coord( X, V );  //actually, dxdxprim() does not use V or X since metric is uniform
    if(gdet_func_metric(MCOORD,V,gcovmid,&gdetmid)!=0){
      if(debugfail>=2) dualfprintf(fail_file,"Caught gdet_func_metric() problem in set_conn_cylminkmetric()\n");
    }



    // see transforms.c and mettometp() and see gcov2gcovprim()  //atch
    dxdxprim(X, V, dxdxp);                                       //atch

    for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
    //conn2[RR]=-1.0/gdetmid;  //wrong as well... shouldn't it be 0 since there is no 2nd connection when WITHGDET?  SUPERSASMARK
 

    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
        for (k = 0; k < NDIM; k++) {
          conn[i][j][k] = 0.;
        }
    conn[PH][RR][PH]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
    conn[PH][PH][RR]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
    conn[RR][PH][PH]= - X[1] * dxdxp[3][3] * dxdxp[3][3]; //-gdetmid;
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }


  
}

/// only works for X1=R and X2=z
void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
                             FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
                            FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);

  int i,j,k;


  if(ANALYTICCONNECTION&&defcoord==UNIFORMCOORDS){// uniform grid
    for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
    
    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
        for (k = 0; k < NDIM; k++) {
          conn[i][j][k] = 0.;
        }
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }

}



void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
                       FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
                            FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2);
  void mks_conn_func(FTYPE *X, struct of_geom *geom,
                     FTYPE (*lconn)[NDIM][NDIM],FTYPE *conn2);

  //  FTYPE DELTA; // for debug below

  // the analytic form can be used with WHICHEOM=WITHNOGDET
  // determine for which we have analytic expressions
  // this currently competes with mks_source_conn in phys.c
  if((!ANALYTICCONNECTION)||(defcoord!=LOGRSINTH)) conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  else mks_conn_func(X,geom,conn,conn2);


  /*
  // some debug stuff
  //  if((geom->i==62)&&(geom->j==32)){ // where c002 is bad
  if(0&&(geom->i==32)&&(geom->j==0)){ // where c023 is bad
  for(DELTA=1E-15;DELTA<1E5;DELTA*=1.01){
  conn_func_numerical1(DELTA,X, geom, conn, conn2);
  dualfprintf(fail_file,"%30.20Lg %30.20Lg\n",DELTA,conn[0][2][3]); fflush(fail_file);
  //      dualfprintf(fail_file,"%21.15g %21.15g\n",DELTA,conn[0][2][3]); fflush(fail_file);
  }
  exit(0);
  }
  else{
  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }
  //mks_conn_func(X,geom,conn,conn2);
  */

}










































//////////////////////////////
//
// below very specific to defcoord==LOGRSINTH and MCOORD=KSCOORDS:
// gives: analytic connection and source
//
/////////////////////////////////



/// jon's MKS connection (and conn2)
/// only applies to defcoord==LOGRSINTH
void mks_conn_func(FTYPE *X, struct of_geom *ptrgeom,
                   FTYPE (*conn)[NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxptrue[NDIM][NDIM];
#ifdef WIN32
  extern FTYPE cot(FTYPE arg);
#endif
  FTYPE dxdxp[NDIM];


  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_conn_func not setup for MBH!=1.0\n");
    myexit(10);
  }


  // get bl coordinates
  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
  r=V[1];
  th=V[2];


  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  sigma=r*r+a*a*cos(th)*cos(th);


  conn[0][0][0]=(-2.*r*sigma + 4.*pow(r,3.))*pow(sigma,-3.);
  conn[0][0][1]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[0][0][2]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
  conn[0][0][3]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[0][1][0]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[0][1][1]=2.*(r + sigma)*pow(dxdxp[1],2.)*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
  conn[0][1][2]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
  conn[0][1][3]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
  conn[0][2][0]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
  conn[0][2][1]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
  conn[0][2][2]=-2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-1.);
  conn[0][2][3]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
  conn[0][3][0]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[0][3][1]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
  conn[0][3][2]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
  conn[0][3][3]=2.*r*pow(sigma,-3.)*pow(sin(th),2.)*
    (-1.*r*pow(sigma,2.) + pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
     pow(sin(th),2.));
  conn[1][0][0]=pow(dxdxp[1],-1.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
  conn[1][0][1]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
  conn[1][0][2]=0.;
  conn[1][0][3]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
  conn[1][1][0]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
  conn[1][1][1]=pow(sigma,-3.)*
    (-1.*dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.)) + 
     pow(sigma,3.) + dxdxp[1]*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
     pow(sin(th),2.));
  conn[1][1][2]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
    ;
  conn[1][1][3]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
                       cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
                       2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[1][2][0]=0.;
  conn[1][2][1]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
    ;
  conn[1][2][2]=-1.*r*pow(dxdxp[1],-1.)*pow(dxdxp[2],2.)*
    pow(sigma,-1.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
  conn[1][2][3]=0.;
  conn[1][3][0]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
  conn[1][3][1]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
                       cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
                       2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[1][3][2]=0.;
  conn[1][3][3]=-1.*pow(dxdxp[1],-1.)*pow(sigma,-3.)*pow(sin(th),2.)*
    (-2.*r + sigma + pow(a,2.)*pow(sin(th),2.))*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[2][0][0]=-1.*r*pow(dxdxp[2],-1.)*pow(a,2.)*pow(sigma,-3.)*
    sin(2.*th);
  conn[2][0][1]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][0][2]=0.;
  conn[2][0][3]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
  conn[2][1][0]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][1][1]=-1.*r*pow(dxdxp[1],2.)*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][1][2]=dxdxp[1]*r*pow(sigma,-1.);
  conn[2][1][3]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
  conn[2][2][0]=0.;
  conn[2][2][1]=dxdxp[1]*r*pow(sigma,-1.);
  conn[2][2][2]=4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*
    pow(M_PI,2.) - 1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)\
    ;
  conn[2][2][3]=0.;
  conn[2][3][0]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
  conn[2][3][1]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
  conn[2][3][2]=0.;
  conn[2][3][3]=-1.*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (pow(sigma,3.) + sigma*(4.*r + sigma)*pow(a,2.)*pow(sin(th),2.) + 
     2.*r*pow(a,4.)*pow(sin(th),4.))*sin(th);
  conn[3][0][0]=a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
  conn[3][0][1]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
    ;
  conn[3][0][2]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
  conn[3][0][3]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[3][1][0]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
    ;
  conn[3][1][1]=a*pow(dxdxp[1],2.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[3][1][2]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
  conn[3][1][3]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[3][2][0]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
  conn[3][2][1]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
  conn[3][2][2]=-1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-1.);
  conn[3][2][3]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
  conn[3][3][0]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[3][3][1]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[3][3][2]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
  conn[3][3][3]=pow(sigma,-3.)*
    (-1.*a*r*pow(sigma,2.)*pow(sin(th),2.) + 
     pow(a,3.)*(-1.*sigma + 2.*pow(r,2.))*pow(sin(th),4.));
  conn2[0]=0.;
  conn2[1]=-1.*pow(sigma,-1.)*(2.*dxdxp[1]*r + pow(r,2.) + 
                               pow(a,2.)*pow(cos(th),2.));
  conn2[2]=-1.*dxdxp[2]*cot(th) + 
    4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
    dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th);
  conn2[3]=0.;


}



/*********************************************************************************************
   Scott's s MKS connection that can be used for any transformation between r,th <-> X1,X2 :
*********************************************************************************************/
#define M (1.)
void mks_conn_func_general(FTYPE *X, struct of_geom *geom, FTYPE (*conn)[NDIM][NDIM] )
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxp[NDIM][NDIM],dxdxp_dxp[NDIM][NDIM][NDIM];

  FTYPE t1,   t10,   t102,   t1024,   t1035,   t1037,   t104,   t11,   t114,   t116,   t119;
  FTYPE  t12,   t121,   t123,   t126,   t129,   t130,   t132,   t14,   t148,   t149,   t15,   t152;
  FTYPE t154,   t156,   t157,   t158,   t159,   t161,   t169,   t17,   t171,   t172,   t175,   t177;
  FTYPE t185,   t2,   t203,   t204,   t208,   t209,   t21,   t210,   t212,   t214,   t22,   t221;
  FTYPE   t222,   t224,   t227,   t23,   t236,   t24,   t240,   t241,   t242,   t243,   t245,   t246;
  FTYPE    t247,   t248,   t25,   t250,   t251,   t258,   t26,   t260,   t261,   t263,   t264,   t271;
  FTYPE   t273,   t275,   t276,   t278,   t28,   t280,   t281,   t283,   t284,   t285,   t286,   t288;
  FTYPE    t289,   t29,   t297,   t299,   t3,   t30,   t300,   t303,   t305,   t306,   t308,   t309;
  FTYPE    t31,   t310,   t313,   t314,   t320,   t325,   t327,   t328,   t329,   t330,   t333,   t336;
  FTYPE    t338,   t34,   t340,   t342,   t344,   t346,   t35,   t358,   t361,   t363,   t366,   t367;
  FTYPE    t368,   t370,   t372,   t375,   t38,   t380,   t381,   t384,   t385,   t387,   t39,   t399;
  FTYPE    t4,   t40,   t400,   t402,   t404,   t405,   t406,   t408,   t409,   t41,   t411,   t412;
  FTYPE    t418,   t42,   t421,   t425,   t428,   t431,   t432,   t433,   t434,   t437,   t440,   t442;
  FTYPE    t448,   t451,   t453,   t454,   t459,   t462,   t467,   t469,   t480,   t481,   t486,   t487;
  FTYPE    t488,   t491,   t492,   t498,   t501,   t504,   t507,   t508,   t510,   t512,   t52,   t521;
  FTYPE    t528,   t53,   t530,   t553,   t556,   t56,   t57,   t588,   t60,   t607,   t627,   t628;
  FTYPE    t63,   t630,   t631,   t632,   t634,   t636,   t637,   t64,   t651,   t652,   t654,   t656;
  FTYPE    t657,   t659,   t661,   t662,   t670,   t673,   t675,   t677,   t686,   t689,   t7,   t712;
  FTYPE    t74,   t748,   t75,   t78,   t793,   t794,   t795,   t799,   t8,   t800,   t801,   t803;
  FTYPE    t806,  t807,   t813,   t816,   t822,   t83,   t831,   t84,   t845,   t86,   t89,   t891; 
  FTYPE    t90,   t91,   t916,   t917,   t920,   t924,   t928,   t940,   t946,   t968,   t97, t970,t991;


  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_conn_func_general not setup for MBH!=1.0\n");
    myexit(11);
  }


  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];

  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim(X,V,dxdxp);
  //  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  //  dx_dxp_calc(X,r,th,dx_dxp);
  //dx_dxp_dxp_calc(X,r,th,dx_dxp_dxp);



  // GODMARK
  // need to set second derivative analytically
  for(i=0;i<NDIM;i++)  for(j=0;j<NDIM;j++)  for(k=0;k<NDIM;k++){
        dxdxp_dxp[i][j][k]=0.0;
      }


  //  t1 = rf(X1,X2);
  t1 = r;
  //  t2 = thf(X1,X2);
  t2 = th;
  t3 = cos(t2);
  t4 = a*t3;
  t7 = (-t1+t4)*(t1+t4);
  t8 = M*M;
  t10 = t1*t1;
  t11 = a*a;
  t12 = t3*t3;
  t14 = t10+t11*t12;
  t15 = t14*t14;
  t17 = 1/t15/t14;
  conn[0][0][0] = -2.0*t7*t8*t1*t17;
  t21 = t10*t10;
  //  t22 = diff(rf(X1,X2),X1);
  t22 = dxdxp[RR][1];
  t23 = t21*t22;
  t24 = t10*t1;
  t25 = M*t24;
  t26 = t25*t22;
  t28 = t24*t3;
  t29 = sin(t2);
  //  t30 = diff(thf(X1,X2),X1);
  t30 = dxdxp[TH][1];
  t31 = t29*t30;
  t34 = M*t1;
  t35 = t22*t12;
  t38 = t12*t12;
  t39 = t38*t22;
  t40 = t29*t1;
  t41 = t12*t3;
  t42 = t41*t30;
  conn[0][0][1] = -M*(-t23-2.0*t26+(2.0*t28*t31+2.0*t34*t35+(t39+2.0*t40*t42)*t11)*t11)*t17;
  //  t52 = diff(rf(X1,X2),X2);
  t52 = dxdxp[RR][2];
  t53 = t21*t52;
  //  t56 = diff(thf(X1,X2),X2);
  t56 = dxdxp[TH][2];
  t57 = t29*t56;
  t60 = t52*t12;
  t63 = t38*t52;
  t64 = t41*t1;
  conn[0][0][2] = -M*(-t53-2.0*t25*t52+(2.0*t28*t57+2.0*t34*t60+(t63+2.0*t64*t57)*t11)*t11)*t17;
  t74 = -1.0+t3;
  t75 = 1.0+t3;
  t78 = t7*t74*t75;
  conn[0][0][3] = -2.0*t78*a*t1*t8*t17;
  conn[0][1][0] = conn[0][0][1];
  t83 = t30*t30;
  t84 = t21*t10;
  t86 = t22*t22;
  t89 = t24*t22;
  t90 = t3*t29;
  t91 = t90*t30;
  t97 = t86*t12;
  t102 = t22*t29;
  t104 = t64*t102;
  conn[0][1][1] = -2.0*(t83*t84-t21*t86-t25*t86+(2.0*t89*t91+2.0*t83*t21*t12+
                                                 t34*t97+(t83*t10*t38+t38*t86+2.0*t104*t30)*t11)*t11)*M*t17;
  t114 = t22*t52;
  t116 = t30*t56;
  t119 = t24*t52;
  t121 = t114*t12;
  t123 = t21*t12;
  t126 = t90*t56;
  t129 = t52*t29;
  t130 = t129*t30;
  t132 = t10*t38;
  conn[0][1][2] = -2.0*(-t25*t114+t116*t84-t23*t52+(t119*t91+t34*t121+2.0*
                                                    t116*t123+t89*t126
                                                    +(t39*t52+t64*t130+t116*t132+t104*t56)*t11)*t11)*M*t17;
  t148 = 2.0*t28*t30;
  t149 = t102*t12;
  t152 = t41*t24;
  t154 = 2.0*t152*t30;
  t156 = 2.0*t64*t30;
  t157 = t39*t29;
  t158 = t30*t1;
  t159 = t38*t3;
  t161 = 2.0*t158*t159;
  t169 = a*t29*t17;
  conn[0][1][3] = -(2.0*t25*t102+t23*t29+(-t148-2.0*t34*t149+t154+(-t156-t157+t161)*t11)*t11)*M*t169;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = conn[0][1][2];
  t171 = t52*t52;
  t172 = t171*M;
  t175 = t56*t56;
  t177 = t1*t12;
  t185 = t52*t41;
  conn[0][2][2] = -2.0*(-t172*t24-t171*t21+t175*t84+(t172*t177+2.0*t119*t126+
                                                     2.0*t175*t21*t12
                                                     +(t171*t38+2.0*t185*t40*t56+t175*t10*t38)*t11)*t11)*M*t17;
  t203 = 2.0*t152*t56;
  t204 = t129*t12;
  t208 = 2.0*t28*t56;
  t209 = t63*t29;
  t210 = t56*t1;
  t212 = 2.0*t210*t159;
  t214 = 2.0*t64*t56;
  conn[0][2][3] = (-2.0*t25*t129-t53*t29+(-t203+2.0*t34*t204+t208+(t209-t212+t214)*t11)*t11)*M*t169;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t221 = t21*t1;
  t222 = t24*t12;
  t224 = t10*t12;
  t227 = t1*t38;
  t236 = (-t221+(-2.0*t222+(-t224+t10)*M+(-t227+(t38-t12)*M)*t11)*t11)*t74*t75;
  conn[0][3][3] = -2.0*t236*t34*t17;
  t240 = t56*M;
  t241 = t240*t24;
  t242 = 2.0*t241;
  t243 = t56*t21;
  t245 = 2.0*t240*t177;
  t246 = t56*t10;
  t247 = t246*t12;
  t248 = t1*t3;
  t250 = 2.0*t248*t129;
  t251 = t56*t12;
  t258 = t30*t52;
  t260 = 1/(-t22*t56+t258);
  t261 = t260*t17;
  conn[1][0][0] = -(-t242+t243+(t245-t247+t250+t246-t251*t11)*t11)*M*t261;
  t263 = M*t22;
  t264 = t56*t38;
  t271 = (-t242+(t245-t247+t250+t246+(-t251+t264)*t11)*t11)*t260*t17;
  conn[1][0][1] = -t263*t271;
  t273 = M*t52;
  conn[1][0][2] = -t273*t271;
  t275 = t24*t29;
  t276 = t240*t275;
  t278 = t119*t3;
  t280 = t243*t29;
  t281 = t40*t12;
  t283 = 2.0*t240*t281;
  t284 = t29*t12;
  t285 = t246*t284;
  t286 = t52*t3;
  t288 = 2.0*t286*t1;
  t289 = t57*t10;
  t297 = a*t29*t260*t17;
  conn[1][0][3] = (-2.0*t276+2.0*t278+t280+(t283-t285+t288+t289-t56*t11*t284)*t11)*M*t297;
  conn[1][1][0] = conn[1][0][1];
  //  t299 = diff(diff(thf(X1,X2),X1),X1);
  t299 = dxdxp_dxp[TH][1][1];
  t300 = t52*t299;
  t303 = t258*t221*t22;
  t305 = t56*t84;
  //  t306 = diff(diff(rf(X1,X2),X1),X1);
  t306 = dxdxp_dxp[RR][1][1] ;
  t308 = t83*t56;
  t309 = t21*t24;
  t310 = t308*t309;
  t313 = 2.0*t308*t84;
  t314 = t56*t24;
  t320 = t30*t22;
  t325 = t258*t89*t12;
  t327 = t308*t221;
  t328 = t52*t83;
  t329 = t90*t21;
  t330 = t328*t329;
  t333 = t306*t12;
  t336 = t221*t12;
  t338 = 2.0*t308*t336;
  t340 = 4.0*t308*t123;
  t342 = t248*t29;
  t344 = 2.0*t52*t86*t342;
  t346 = t56*t86;
  t358 = t306*t38;
  t361 = t1*t22;
  t363 = t258*t361*t38;
  t366 = 2.0*t308*t222;
  t367 = t24*t38;
  t368 = t308*t367;
  t370 = t41*t29*t10;
  t372 = 2.0*t328*t370;
  t375 = 2.0*t308*t132;
  t380 = t159*t29;
  t381 = t380*t56;
  t384 = t308*t227;
  t385 = t38*t12;
  t387 = t328*t380;
  conn[1][1][1] = -(-t300*t84-2.0*t303+t305*t306-t310+(-t243*t86+t313-2.0*t314*t86*M)*M
                    +(-2.0*t320*t3*t280-4.0*t325-t327+t330-3.0*t300*t123+3.0*t243*t333
                      -t338+(t340+t344-t246*t97+t346*t10+2.0*t210*t97*M)*M
                      +(-4.0*t320*t41*t289-3.0*t300*t132+3.0*t246*t358-2.0*t363-t366-t368+t372
                        +(-t346*t12+t375+2.0*t346*t38)*M
                        +(-2.0*t320*t381-t384-t300*t385+t387+t56*t306*t385)*t11)*t11)*t11)*t260*t17;
  //  t399 = diff(diff(thf(X1,X2),X1),X2);
  t399 = dxdxp_dxp[TH][1][2]; 
  t400 = t52*t399;
  //  t402 = diff(diff(rf(X1,X2),X1),X2);
  t402 = dxdxp_dxp[RR][1][2]; 
  t404 = t30*t175;
  t405 = t404*t309;
  t406 = t171*t30;
  t408 = t56*t221;
  t409 = t114*t408;
  t411 = 2.0*t404*t84;
  t412 = t52*t56;
  t418 = t402*t12;
  t421 = t404*t221;
  t425 = t114*t314*t12;
  t428 = 2.0*t404*t336;
  t431 = t22*t175;
  t432 = t431*t329;
  t433 = t10*t22;
  t434 = t433*t12;
  t437 = 4.0*t404*t123;
  t440 = 2.0*t171*t22*t342;
  t442 = t210*M;
  t448 = 2.0*t370*t431;
  t451 = t114*t210*t38;
  t453 = 2.0*t404*t222;
  t454 = t402*t38;
  t459 = t404*t367;
  t462 = 2.0*t404*t132;
  t467 = t404*t227;
  t469 = t431*t380;
  conn[1][1][2] = (t400*t84-t305*t402+t405+t406*t221+t409+(-t411+t412*t23+2.0*t114*t241)*M
                   +(-3.0*t243*t418+t421+3.0*t400*t123+2.0*t425+t428+2.0*t406*t222+
                     t432+(t412*t434-t437-t440-t412*t433-2.0*t121*t442)*M
                     +(t448+t406*t227+t451+t453-3.0*t246*t454+3.0*t400*t132+t459
                       +(t114*t251-t462-2.0*t114*t264)*M+(t467+t400*t385+t469
                                                          -t56*t402*t385)*t11)*t11)*t11)*t260*t17;
  t480 = t286*t21;
  t481 = t57*t221;
  t486 = 2.0*t185*t10;
  t487 = t57*t222;
  t488 = 2.0*t487;
  t491 = t57*t227;
  t492 = t52*t159;
  t498 = (-t491+t492+(-t57*t12+t57*t38)*M)*t11;
  t501 = t480-t481+(2.0*t278-2.0*t276)*M+(t486-t488+(-t285+t289+t288+t283)*M+t498)*t11;
  conn[1][1][3] = t501*t22*t297;
  conn[1][2][0] = conn[1][0][2];
  conn[1][2][1] = conn[1][1][2];
  t504 = t171*t56;
  //  t507 = diff(diff(thf(X1,X2),X2),X2);
  t507 = dxdxp_dxp[TH][2][2];
  t508 = t52*t507;
  //  t510 = diff(diff(rf(X1,X2),X2),X2);
  t510 = dxdxp_dxp[RR][2][2];
  t512 = t175*t56;
  t521 = t512*t221;
  t528 = t52*t175;
  t530 = t510*t12;
  t553 = t510*t38;
  t556 = t512*t24;
  conn[1][2][2] = -(-2.0*t504*t221-t508*t84+t305*t510-t512*t309
                    +(2.0*t512*t84-t504*t21-2.0*t504*t25)*M
                    +(-2.0*t521*t12-3.0*t508*t123-4.0*t504*t222-t521-t528*t329+3.0*t243*t530
                      +(-t504*t224+t504*t10+2.0*t342*t171*t52+4.0*t512*t21*t12+2.0*t12*t171*t442)*M
                      +(-3.0*t508*t132-2.0*t504*t227-2.0*t528*t370+3.0*t246*t553-2.0*t556*t12-t556*t38
                        +(2.0*t512*t10*t38-t504*t12+2.0*t504*t38)*M
                        +(-t528*t380-t512*t1*t38-t508*t385+t56*t510*t385)*t11)*t11)*t11)*t260*t17;
  conn[1][2][3] = t501*t52*t297;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = conn[1][2][3];
  t588 = t84*t29;
  t607 = t29*t38;
  conn[1][3][3] = -(-t56*t309*t29+t286*t84+2.0*t240*t588
                    +(-t481-2.0*t408*t284+t480+2.0*t185*t21+(-4.0*t185*t24+t280+3.0*t243*t284+4.0*t278
                                                             +(-2.0*t314*t29+2.0*t487)*M)*M
                      +(t492*t10+t486-t314*t607-t488+(-2.0*t492*t1+t289+t288-2.0*t285+3.0*t246*t607
                                                      +(-2.0*t491+2.0*t210*t284)*M)*M+t498)*t11)*t11)*t29*t261;
  t627 = t30*t21;
  t628 = t30*M;
  t630 = 2.0*t628*t24;
  t631 = t30*t10;
  t632 = t631*t12;
  t634 = 2.0*t628*t177;
  t636 = 2.0*t361*t90;
  t637 = t30*t12;
  conn[2][0][0] = -(-t627+t630+(t632-t634-t631-t636+t637*t11)*t11)*M*t261;
  t651 = (-t630+(t634+t636+t631-t632+(-t637+t30*t38)*t11)*t11)*t260*t17;
  conn[2][0][1] = t263*t651;
  conn[2][0][2] = t273*t651;
  t652 = t628*t275;
  t654 = t89*t3;
  t656 = t627*t29;
  t657 = t631*t284;
  t659 = 2.0*t628*t281;
  t661 = 2.0*t361*t3;
  t662 = t31*t10;
  conn[2][0][3] = (2.0*t652-2.0*t654-t656+(t657-t659-t661-t662+t30*t11*t284)*t11)*M*t297;
  conn[2][1][0] = conn[2][0][1];
  t670 = t30*t86;
  t673 = t30*t84;
  t675 = t22*t299;
  t677 = t83*t30;
  t686 = t677*t221;
  t689 = t83*t22;
  t712 = t677*t24;
  conn[2][1][1] = -(2.0*t670*t221-t673*t306+t675*t84+t677*t309+(-2.0*t677*t84+t627*t86+2.0*t670*t25)*M
                    +(2.0*t686*t12+t689*t329-3.0*t627*t333+t686+4.0*t670*t222+3.0*t675*t123
                      +(-2.0*t86*t22*t1*t90+t631*t97-t670*t10-4.0*t677*t21*t12-2.0*t637*t86*t34)*M
                      +(2.0*t712*t12+t712*t38+2.0*t670*t227+3.0*t675*t132-3.0*t631*t358+2.0*t689*t370
                        +(t670*t12-2.0*t677*t10*t38-2.0*t670*t38)*M
                        +(t675*t385-t30*t306*t385+t689*t380+t677*t1*t38)*t11)*t11)*t11)*t260*t17;
  t748 = t22*t399;
  conn[2][1][2] = -(t303+t310-t673*t402+t748*t84+t346*t221
                    +(-t313+t258*t23+2.0*t258*t26)*M
                    +(t327+t330+2.0*t325-3.0*t627*t418+2.0*t346*t222+3.0*t748*t123+t338
                      +(-t340-t258*t433-t344+t258*t434-2.0*t60*t30*t361*M)*M
                      +(t346*t227+t372+t363+t366+t368-3.0*t631*t454+3.0*t748*t132
                        +(t258*t35-t375-2.0*t258*t39)*M+(t387+t748*t385+t384
                                                         -t30*t402*t385)*t11)*t11)*t11)*t260*t17;
  t793 = t31*t221;
  t794 = t3*t22;
  t795 = t794*t21;
  t799 = t31*t222;
  t800 = 2.0*t799;
  t801 = t22*t41;
  t803 = 2.0*t801*t10;
  t806 = t31*t227;
  t807 = t22*t159;
  t813 = (-t806+t807+(t31*t38-t31*t12)*M)*t11;
  t816 = -t793+t795+(2.0*t654-2.0*t652)*M+(-t800+t803+(t661-t657+t662+t659)*M+t813)*t11;
  conn[2][1][3] = -t816*t22*t297;
  conn[2][2][0] = conn[2][0][2];
  conn[2][2][1] = conn[2][1][2];
  t822 = t22*t507;
  t831 = t3*t56;
  t845 = t41*t56;
  conn[2][2][2] = -(-t673*t510+2.0*t409+t405+t822*t84+(t406*t21-t411+2.0*t406*t25)*M
                    +(-3.0*t627*t530+t421-t432+2.0*t130*t831*t21+t428+3.0*t822*t123+4.0*t425
                      +(t406*t224-t406*t10-t440-t437-2.0*t406*t177*M)*M
                      +(4.0*t130*t845*t10+2.0*t451+3.0*t822*t132+t453+t459-t448-3.0*t631*t553
                        +(t406*t12-t462-2.0*t406*t38)*M+(2.0*t258*t381+t822*t385-t30*t510*t385
                                                         +t467-t469)*t11)*t11)*t11)*t260*t17;
  conn[2][2][3] = -t816*t52*t297;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = conn[2][2][3];
  t891 = t30*t24;
  conn[2][3][3] = (t794*t84+2.0*t628*t588-t30*t309*t29
                   +(t795-t793-2.0*t30*t221*t284+2.0*t801*t21
                     +(3.0*t627*t284+t656-4.0*t801*t24+4.0*t654+(-2.0*t891*t29+2.0*t799)*M)*M
                     +(t807*t10+t803-t891*t607-t800+(3.0*t631*t607-2.0*t807*t1+t661-2.0*t657+t662
                                                     +(2.0*t158*t284-2.0*t806)*M)*M+t813)*t11)*t11)*t29*t261;
  t917 = a*M;
  conn[3][0][0] = -t7*t917*t17;
  t920 = t102*t10;
  t924 = 1/t29;
  t916 = t924*t17;
  conn[3][0][1] = -t917*(-t920+t148+(t156+t149)*t11)*t916;
  t928 = t129*t10;
  conn[3][0][2] = -t917*(t208-t928+(t214+t204)*t11)*t916;
  conn[3][0][3] = -t78*M*t11*t17;
  conn[3][1][0] = conn[3][0][1];
  t940 = t83*t29;
  t946 = t86*t29;
  t968 = t916;
  conn[3][1][1] = -(t940*t221+2.0*t794*t627+(4.0*t794*t891-t946*t10)*M
                    +(4.0*t801*t631+2.0*t940*t222+(4.0*t361*t42+t946*t12)*M+(t940*t227+2.0*t807*t30)*t11)
                    *t11)*a*t968;
  t970 = t29*t221;
  t991 = t52*t1;
  conn[3][1][2] = -(t116*t970+t286*t627+t794*t243
                    +(2.0*t286*t891-t102*t52*t10+2.0*t794*t314)*M
                    +(2.0*t185*t631+2.0*t801*t246+2.0*t116*t275*t12+(2.0*t361*t845+2.0*t991*t42+t102*t60)*M
                      +(t116*t40*t38+t492*t30+t807*t56)*t11)*t11)*a*t968;
  t1024 = t38*t41;
  conn[3][1][3] = (t3*t30*t84+t970*t22+(3.0*t42*t21+2.0*t275*t35
                                        +(t102*t224+t148-t154-t920)*M
                                        +(t40*t39+3.0*t159*t30*t10+(t156-t161+t149-t157)*M
                                          +t1024*t30*t11)*t11)*t11)*t924*t17;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  t1035 = t175*t29;
  t1037 = t171*t29;
  conn[3][2][2] = -(2.0*t286*t243+t1035*t221+(-t1037*t10+4.0*t286*t314)*M
                    +(4.0*t185*t246+2.0*t1035*t222+(t1037*t12+4.0*t991*t845)*M
                      +(t1035*t227+2.0*t492*t56)*t11)*t11)*a*t968;
  conn[3][2][3] = -(-t970*t52-t831*t84+(-2.0*t275*t60-3.0*t845*t21
                                        +(t203-t208-t129*t224+t928)*M
                                        +(-t40*t63-3.0*t159*t56*t10+(t209+t212-t204-t214)*M
                                          -t1024*t56*t11)*t11)*t11)*t924*t17;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = -t236*a*t17;

  return;

}

#undef M






/// jon's MKS source in mid-simplified form (used to UPDATE the source)
/// only for ideal EOS (gam)
void mks_source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU)
{
  int ii,jj,kk;
  int i=0, j=0, k=0, l=0;
  FTYPE r,th,X[NDIM],V[NDIM],sigma,dxdxptrue[NDIM][NDIM];
#ifdef WIN32
  extern FTYPE cot(FTYPE arg);
#endif
  extern FTYPE csc(FTYPE arg);
  FTYPE b[NDIM],u[NDIM],bsq,en,rho;
  FTYPE dxdxp[NDIM];



  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_source_conn not setup for MBH!=1.0\n");
    myexit(12);
  }


  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;


  bsq = dot(q->bcon, q->bcov);
  u[TT]=q->ucon[TT];
  u[RR]=q->ucon[RR];
  u[TH]=q->ucon[TH];
  u[PH]=q->ucon[PH];

  b[TT]=q->bcon[TT];
  b[RR]=q->bcon[RR];
  b[TH]=q->bcon[TH];
  b[PH]=q->bcon[PH];

  rho=pr[RHO];
  en=pr[UU];

  coord(ptrgeom->i, ptrgeom->j, ptrgeom->k,ptrgeom->p, X);
  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];


  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif
  // set aux vars
  dxdxprim(X,V,dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal


  sigma=r*r+a*a*cos(th)*cos(th);



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU0==1)){
    // see grmhd-fullsource-simplify.nb

    dU[UU]+=pow(sigma,-1.)*(-1.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
                            (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                                    r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
                                       2.*b[PH]*a*pow(sin(th),2.))) + 
                             u[RR]*(bsq + en*gam + rho)*
                             (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                              r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                                 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
                            (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                                    r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                                       2.*b[PH]*a*pow(sin(th),2.))) - 
                             1.*u[TH]*(bsq + en*gam + rho)*
                             (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                              r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                                 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*
                            (-1.*dxdxp[2]*cot(th) + 
                             4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
                             dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


  if((WHICHEOM==WITHNOGDET)&&(NOGDETU1==1)){

    // see grmhd-ksksp-mhd.nb and the other source*simplify*.nb files

    dU[U1]+=0.0625*(16.*dxdxp[1]*r*
                    (0.5*bsq + en*(-1. + gam) - 
                     1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
                     (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
                    pow(sigma,-1.) + 2.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
                    (4.*bsq + 8.*en*(-1. + gam) - 
                     1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                                        8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                                        cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                                        1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                                        8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
                     pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
                     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                      4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.)) - 
                    1.*dxdxp[1]*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
                    (-1.*b[TT]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                                8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                                cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                                1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                                8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
                     u[TT]*(bsq + en*gam + rho)*
                     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                      4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
                    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)) + 
                    8.*a*pow(dxdxp[1],2.)*(-1.*b[RR]*
                                           (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                                            b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                            cos(2.*th)*pow(a,2.)*
                                            (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                                             b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                                            b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
                                           u[RR]*(bsq + en*gam + rho)*
                                           (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                                                     dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                            cos(2.*th)*pow(a,2.)*
                                            (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
                                             u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
                                            u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
                    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
                    8.*dxdxp[1]*a*(-1.*b[TT]*
                                   (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                                    b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                    cos(2.*th)*pow(a,2.)*
                                    (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                                     b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                                    b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
                                   u[TT]*(bsq + en*gam + rho)*
                                   (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                                             dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                    cos(2.*th)*pow(a,2.)*
                                    (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
                                     u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
                                    u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
                    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
                    dxdxp[1]*a*(-1.*b[PH]*
                                (16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
                                 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
                                 4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
                                 b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
                                 4.*b[PH]*a*pow(r,2.)) + 
                                u[PH]*(bsq + en*gam + rho)*
                                (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                                 4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                                 4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                                 u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                                 4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
                    (pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
                     2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
                     cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
                    pow(sin(th),2.) + 8.*dxdxp[1]*pow(sigma,-3.)*
                    (bsq + 2.*en*(-1. + gam) - 
                     1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                               b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                               cos(2.*th)*pow(a,2.)*
                               (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                                b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                               b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*
                     pow(sin(th),2.) + u[PH]*(bsq + en*gam + rho)*
                     (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                      cos(2.*th)*pow(a,2.)*
                      (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
                       u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
                      u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*
                     pow(sin(th),2.))*(r*pow(sigma,2.) + 
                                       pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
                                       ) + 2.*pow(sigma,-3.)*(4.*bsq + 8.*en*(-1. + gam) - 
                                                              1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                                                                                 8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                                                                                 cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                                                                                 1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                                                                                 8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
                                                              pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
                                                              (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                                                               4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                                                               4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                                                               u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                                                               4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*
                    (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
                     dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
                                         pow(a,2.)*pow(cos(th),4.) + 
                                         (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
                    16.*dxdxp[1]*a*pow(sigma,-4.)*
                    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                    (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
                    (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                               2.*b[PH]*a*pow(sin(th),2.))) - 
                     1.*u[PH]*(bsq + en*gam + rho)*
                     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                         u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
                    32.*pow(dxdxp[1],2.)*pow(sigma,-4.)*
                    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
                    (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
                    (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                            r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
                               2.*b[PH]*a*pow(sin(th),2.))) + 
                     u[RR]*(bsq + en*gam + rho)*
                     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                         u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
                    16.*dxdxp[1]*pow(sigma,-3.)*
                    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
                    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                    (0.5*bsq + en*(-1. + gam) + 
                     b[TT]*pow(sigma,-1.)*
                     (b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                      r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                         2.*b[PH]*a*pow(sin(th),2.))) - 
                     1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
                     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                         u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
                    2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*
                    (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                                8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                                cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                                1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                                8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
                     u[TH]*(bsq + en*gam + rho)*
                     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                      4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.)*sin(th) - 
                    8.*dxdxp[1]*dxdxp[2]*a*cos(th)*
                    (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                                b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                cos(2.*th)*pow(a,2.)*
                                (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                                 b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                                b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
                     u[TH]*(bsq + en*gam + rho)*
                     (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                               dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                      cos(2.*th)*pow(a,2.)*
                      (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
                       u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
                      u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
                    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
                    16.*dxdxp[1]*dxdxp[2]*r*
                    (b[TH]*b[TT] - 1.*u[TH]*u[TT]*(bsq + en*gam + rho))*pow(a,2.)*
                    pow(sigma,-2.)*sin(2.*th) + 
                    16.*dxdxp[2]*r*(b[RR]*b[TH] - 
                                    1.*u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[1],2.)*pow(a,2.)*
                    pow(sigma,-2.)*sin(2.*th) - 
                    16.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
                    (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                            r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                               2.*b[PH]*a*pow(sin(th),2.))) - 
                     1.*u[TH]*(bsq + en*gam + rho)*
                     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
                      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
                         u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
                    2.*dxdxp[1]*(-1.*b[TH]*
                                 (16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
                                  4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
                                  4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
                                  b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
                                  4.*b[PH]*a*pow(r,2.)) + 
                                 u[TH]*(bsq + en*gam + rho)*
                                 (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
                                  4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
                                  4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
                                  u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
                                  4.*u[PH]*a*pow(r,2.)))*pow(sigma,-1.)*
                    (-1.*dxdxp[2]*cot(th) + 
                     4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
                     dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)) - 
                    16.*dxdxp[1]*dxdxp[2]*a*
                    (b[PH]*b[TH] - 1.*u[PH]*u[TH]*(bsq + en*gam + rho))*
                    pow(sigma,-2.)*sin(th)*(r*(2. + r)*sigma*cos(th) + 
                                            sigma*pow(a,2.)*pow(cos(th),3.) + r*pow(a,2.)*sin(th)*sin(2.*th)));

  }
  else{

    dU[U1]+=0.5*dxdxp[1]*pow(sigma,-5.)*
      (r*pow(sigma,4.)*(bsq - 2.*en + 2.*en*gam - 
                        2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
                        (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
                        (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.)) + 
       a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
        (r*(2. + 3.*r)*pow(a,2.) + 
         cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
         2.*pow(r,4.)) - 2.*a*
        ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
         0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
         (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
        4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
                      (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
                      (bsq + en*gam + rho)*pow(u[TT],2.)))*pow(sin(th),2.) + 
       dxdxp[1]*a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
                 dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
        pow(dxdxp[1],-1.) + 
        sigma*(r*(2. + 3.*r)*pow(a,2.) + 
               cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
               2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
                              0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
        - 1.*a*pow(dxdxp[1],-1.)*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
        ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
         2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
         (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
         (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       pow(sin(th),2.) - 1.*sigma*
       (4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
         dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*pow(sigma,-1.)*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
        2.*r*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
              (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
              (bsq + en*gam + rho)*pow(u[TT],2.)) - 
        1.*a*(-1.*b[TT]*b[PH] + (bsq + en*gam + rho)*u[TT]*u[PH])*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
       (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
        1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
              1.*dxdxp[1]*b[RR]*b[PH]*
              (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
              dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
              (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
        0.5*(r*(2. + 3.*r)*pow(a,2.) + 
             cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
             2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
                            2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
        )*pow(sin(th),2.)*(r*pow(sigma,2.) + 
                           pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
                           ) + 2.*a*sigma*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
           2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
        - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
        ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
        2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
                pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
                (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.)) + 
       a*sigma*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
                2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
                cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
       pow(sin(th),2.)*(2.*r*(-1.*b[TT]*b[PH] + 
                              (bsq + en*gam + rho)*u[TT]*u[PH]) + 
                        dxdxp[1]*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
                                  0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
                        *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
                        1.*a*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                        (0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*pow(csc(th),2.) - 
                         1.*pow(b[PH],2.) + (bsq + en*gam + rho)*pow(u[PH],2.))*
                        pow(sin(th),2.)) + sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
        2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
        1.*(bsq + en*gam + rho)*
        (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
        4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
        4.*r*(bsq + en*gam + rho)*u[TT]*
        (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))) + 
       2.*sigma*(2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
                       (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
                 dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                 (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
                  ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
                  (bsq + en*gam + rho)*pow(u[RR],2.)) - 
                 1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
                       0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
                 *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*
       (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
        dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
                            pow(a,2.)*pow(cos(th),4.) + 
                            (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
       4.*dxdxp[1]*sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
        2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
        0.5*(bsq + en*gam + rho)*u[RR]*
        (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
         4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))) - 
       1.*dxdxp[2]*a*cos(th)*pow(sigma,2.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
        1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
        (r*(2. + 3.*r)*pow(a,2.) + 
         cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
         2.*pow(r,4.)) + 2.*dxdxp[1]*a*
        (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*sin(th) - 
       2.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,3.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
        dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
        1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
        (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th) + 
       2.*dxdxp[2]*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
       pow(a,2.)*pow(sigma,3.)*sin(2.*th) + 
       2.*dxdxp[1]*dxdxp[2]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(a,2.)*pow(sigma,3.)*
       sin(2.*th) - 2.*dxdxp[2]*r*pow(a,2.)*pow(sigma,2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
        1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
        ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
        2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
        )*sin(2.*th) - 2.*dxdxp[2]*a*
       (b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sigma,3.)*sin(th)*
       (r*(2. + r)*sigma*cos(th) + sigma*pow(a,2.)*pow(cos(th),3.) + 
        r*pow(a,2.)*sin(th)*sin(2.*th)));
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU2==1)){

    dU[U2]+=dxdxp[1]*r*(-1.*b[RR]*b[TH] + 
                        u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[2],2.) - 
      1.*(-1.*b[RR]*b[TH] + u[RR]*u[TH]*(bsq + en*gam + rho))*
      (2.*dxdxp[1]*r + sigma)*pow(dxdxp[2],2.) - 
      0.125*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*
      (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                  8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                  cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                  1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                  8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
       u[TH]*(bsq + en*gam + rho)*
       (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
        4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
        4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
        u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
        4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.) - 
      0.5*a*r*pow(dxdxp[2],2.)*(-1.*b[TH]*
                                (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                                 b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                 cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
                                                       b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
                                 1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
                                 2.*b[PH]*pow(r,4.)) + 
                                u[TH]*(bsq + en*gam + rho)*
                                (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
                                 u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                 cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
                                                       u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
                                 1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
                                 2.*u[PH]*pow(r,4.)))*pow(sigma,-2.)*pow(sin(th),2.) - 
      2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
      (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
              r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                 2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TH]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
        r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
           u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
      2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
      (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
              r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                 2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[PH]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
        r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
           u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*pow(sin(th),3.) - 
      1.*dxdxp[2]*a*r*cos(th)*(-1.*b[TT]*
                               (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                                b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
                                                      b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
                                1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
                                2.*b[PH]*pow(r,4.)) + 
                               u[TT]*(bsq + en*gam + rho)*
                               (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
                                u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                                cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
                                                      u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
                                1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
                                2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*sin(th) - 
      0.125*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*
      (4.*bsq + 8.*en*(-1. + gam) - 
       1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
                          8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
                          cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
                          1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
                          8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
       pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
       (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
        4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
        4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
        u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
        4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*sin(th) - 
      0.5*dxdxp[1]*dxdxp[2]*a*cos(th)*
      (-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                  b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                  cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
                                        b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
                  1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
                  2.*b[PH]*pow(r,4.)) + 
       u[RR]*(bsq + en*gam + rho)*
       (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
        u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
        cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
                              u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
        1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
        2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
      (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
      (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
      (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
       1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) + 
      dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
      (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
              r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
                 2.*b[PH]*a*pow(sin(th),2.))) + 
       u[RR]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
        r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
           u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) - 
      1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
      (0.5*bsq + en*(-1. + gam) + 
       b[TT]*pow(sigma,-1.)*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
                             r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
                                2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
        r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
           u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
      0.5*dxdxp[2]*(bsq + 2.*en*(-1. + gam) - 
                    1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                              b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                              cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
                                                    b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
                              1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
                              2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.) + 
                    u[PH]*(bsq + en*gam + rho)*
                    (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
                     u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                     cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
                                           u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
                     1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
                     2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.))*
      (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)) + 
      (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
      (-1.*dxdxp[2]*cot(th) + 4.*(-1.*M_PI*X[2] + th)*
       pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
       dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th));
  }
  else{


    dU[U2]+=0.5*(-2.*dxdxp[1]*r*
                 (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(dxdxp[2],2.) - 
                 1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-2.)*
                 (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
                  1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
                  (r*(2. + 3.*r)*pow(a,2.) + 
                   cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
                   2.*pow(r,4.)) + 2.*dxdxp[1]*a*
                  (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*pow(sin(th),2.) - 
                 4.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
                 (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
                  1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
                  ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
                  2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
                  ) - 2.*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*pow(sigma,-2.)*
                 (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
                  dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
                  1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
                 4.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
                 (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
                     2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
                  - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
                  ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
                  2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
                          pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
                          (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.))*pow(sin(th),3.)
                 - 2.*dxdxp[2]*a*r*cos(th)*pow(sigma,-4.)*
                 (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
                  (r*(2. + 3.*r)*pow(a,2.) + 
                   cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
                   2.*pow(r,4.)) - 2.*a*
                  ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
                   0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
                   (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
                  4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
                                (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
                                (bsq + en*gam + rho)*pow(u[TT],2.)))*sin(th) - 
                 1.*dxdxp[1]*dxdxp[2]*a*cos(th)*pow(sigma,-4.)*
                 (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                 (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
                           dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
                  pow(dxdxp[1],-1.) + 
                  sigma*(r*(2. + 3.*r)*pow(a,2.) + 
                         cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
                         2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
                                        0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
                  - 1.*a*pow(dxdxp[1],-1.)*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                  ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
                   2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
                   (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
                   (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
                 sin(th) - 2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-2.)*
                 (2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
                        (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
                  dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
                  (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
                   ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
                   (bsq + en*gam + rho)*pow(u[RR],2.)) - 
                  1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
                        0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
                  *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th)\
                 + (bsq - 2.*en + 2.*en*gam - 2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
                    (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
                    (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.))*
                 (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
                  1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) - 
                 1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
                 ((bsq + 2.*en*(-1. + gam))*sigma + 
                  2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
                  1.*(bsq + en*gam + rho)*
                  (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
                  4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
                  4.*r*(bsq + en*gam + rho)*u[TT]*
                  (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))*sin(2.*th) - 
                 2.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
                 (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
                  2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
                  0.5*(bsq + en*gam + rho)*u[RR]*
                  (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
                   4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))))*sin(2.*th) + 
                 dxdxp[2]*pow(sigma,-2.)*
                 (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
                  1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
                        1.*dxdxp[1]*b[RR]*b[PH]*
                        (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
                        dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
                        (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
                  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
                  0.5*(r*(2. + 3.*r)*pow(a,2.) + 
                       cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
                       2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
                                      2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
                  )*pow(sin(th),2.)*(cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)))\
      ;
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU3==1)){

    dU[U3]+=0.5*pow(sigma,-1.)*pow(sin(th),2.)*
      ((-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                   b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                   cos(2.*th)*pow(a,2.)*
                   (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                    b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                   b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
        u[RR]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                  dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
         cos(2.*th)*pow(a,2.)*
         (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
          u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
         u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1. - 2.*dxdxp[1]*r*pow(sigma,-1.)) + 
       (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
                   b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
                   cos(2.*th)*pow(a,2.)*
                   (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
                    b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
                   b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
        u[TH]*(bsq + en*gam + rho)*
        (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
                  dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
         cos(2.*th)*pow(a,2.)*
         (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
          u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
         u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1.*dxdxp[2]*cot(th) + 
        4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
        dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


}



/// GODMARK: remove? Worthless idea? 2D only also.
/// jon's volume diff for uni theta grid and exp r grid (defcoord==LOGRUNITH)
void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol)
{

  /*
    int k, l;
    FTYPE r,th,ph;
    FTYPE r1[2],th1[2],r2[2],th2[2];
    FTYPE X0[NDIM],X1[2][NDIM],X2[2][NDIM];
    //  FTYPE cot(FTYPE arg);


    coord(i, j, CENT, X0);
    coord(i, j, FACE1, X1[0]);
    coord(i+1, j, FACE1, X1[1]);
    coord(i, j, FACE2, X2[0]);
    coord(i, j+1, FACE2, X2[1]);

    // get bl coordinates
    bl_coord(X0,&r,&th);
    bl_coord(X1[0],&r1[0],&th1[0]);
    bl_coord(X1[1],&r1[1],&th1[1]);
    bl_coord(X2[0],&r2[0],&th2[0]);
    bl_coord(X2[1],&r2[1],&th2[1]);

  */




  /* comment out for now until adjust everything


     if(MBH!=1.0){
     dualfprintf(fail_file,"mks_unitheta_idxvold_func not setup for MBH!=1.0\n");
     myexit(13);
     }


     // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
     // follows gcov_func()
     if(POSDEFMETRIC){
     if(th<0.0){ th=-th;}
     if(th>M_PI) { th=M_PI-th; }
     }
     else{
     }
     // avoid singularity at polar axis
     #if(COORDSINGFIX)
     if(fabs(th)<SINGSMALL){
     if(th>=0) th=SINGSMALL;
     if(th<0) th=-SINGSMALL;
     }
     if(fabs(M_PI-th)<SINGSMALL){
     if(th>=M_PI) th=M_PI+SINGSMALL;
     if(th<M_PI) th=M_PI-SINGSMALL;
     }
     #endif
  */

#define IDXR(a,R0,r,th,rl,rh) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/((rh - 1.*rl)*(2.*R0 + rh + rl) + (pow(a,2) + 2.*pow(R0,2))*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl)) + pow(a,2)*cos(2.*th)*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl))))

#define IDXTH(a,R0,r,th,thl,thh) ((-3.*M_PI*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/(cos(thh)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thh)) - 1.*cos(thl)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thl))))


  //#define IDXTH(a,R0,r,th,thl,thh) ((3.*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/((-1.*cos(thh) + cos(thl))*(2.*(pow(a,2) + 3.*pow(r,2)) + pow(a,2)*(cos(2.*thh) + 2.*cos(thh)*cos(thl) + cos(2.*thl)))))

  // fullth integrated -- not used currently
#define FIDXR(a,R0,r,th) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/(4.*R0*(-1.*R0 + r) + pow(-1.*R0 + r,2) + (pow(a,2) + 2.*pow(R0,2) + pow(a,2)*cos(2.*th))*log(-1.*R0 + r)))

#define FIDXTH(a,R0,r,th) ((-3.*pow(a,2)*sin(2.*th) - 6.*pow(r,2)*tan(th))/(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*th)))


  /*

    idxvol[TT]=1.0; // really 1/dt, but changes in time
    // comment out non-volume regularization
    //  idxvol[RR]=1.0/dx[1];
    //idxvol[TH]=1.0/dx[2];
    idxvol[RR]=IDXR(a,R0,r,th,r1[0],r1[1]);
    idxvol[TH]=IDXTH(a,R0,r,th,th2[0],th2[1]);
    idxvol[PH]=1.0/dx[3];

    dualfprintf(fail_file,"%d %d %21.15g %21.15g\n",i,j,idxvol[RR]*dx[1],idxvol[TH]*dx[2]);
  
  */




}


/// geometry only contains i,j,k,p
/// only for spherical polar coordinates with negligible relativistic effects
/// only used if GDETVOLDIFF==1
void gdetvol_func(struct of_geom *ptrgeom, FTYPE *gdettrue, FTYPE *EOMFUNCNAME, FTYPE *gdetvol)
{
  int i,j,k,loc;
  FTYPE Xc[NDIM],Vc[NDIM];
  FTYPE Xim[NDIM],Vim[NDIM],Xip1[NDIM],Vip1[NDIM];
  FTYPE Xjm[NDIM],Vjm[NDIM],Xjp1[NDIM],Vjp1[NDIM];
  FTYPE Xkm[NDIM],Vkm[NDIM],Xkp1[NDIM],Vkp1[NDIM];
  FTYPE dr,dh,dp;
  FTYPE drold,dhold,dpold;
  FTYPE newjacvol,oldjacvol;
  FTYPE dxdxpc[NDIM][NDIM];
  int locarray[4];
  int pl,pliter;


  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;


  // loc==CENT 

  if(loc==CENT){ // e.g. for conserved quantities at CENT
    locarray[0]=loc;
    locarray[1]=FACE1; // left-right
    locarray[2]=FACE2; // up-down
    locarray[3]=FACE3; // in-out
  }
  else if(loc==FACE1){ // e.g. for fluxes F1
    locarray[0]=loc;
    locarray[1]=CENT;
    locarray[2]=CORN3;
    locarray[3]=CORN2;    
  }
  else if(loc==FACE2){ // e.g. for fluxes F2
    locarray[0]=loc;
    locarray[1]=CORN3;
    locarray[2]=CENT;
    locarray[3]=CORN1;    
  }
  else if(loc==FACE3){ // e.g. for fluxes F3
    locarray[0]=loc;
    locarray[1]=CORN2;
    locarray[2]=CORN1;
    locarray[3]=CENT;    
  }
  else{
    // maybe don't need yet :: only for fields -- not used yet.
    locarray[0]=loc;
    locarray[1]=FACE1; // left-right
    locarray[2]=FACE2; // up-down
    locarray[3]=FACE3; // in-out
  }
    

  // rest of code independent of loc    

  coord_ijk(i, j, k, locarray[0], Xc);
  bl_coord_ijk(i, j, k, locarray[0], Vc);
  dxdxprim_ijk(i, j, k, locarray[0], dxdxpc);



 
  ////////////////////////
  //
  // NEW JAC
  coord_ijk(i, j, k, locarray[1], Xim);
  bl_coord_ijk(i, j, k, locarray[1], Vim);

  coord_ijk(i+1, j, k, locarray[1], Xip1); // ok to have i+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i+1, j, k, locarray[1], Vip1);
    
  coord_ijk(i, j, k, locarray[2], Xjm);
  bl_coord_ijk(i, j, k, locarray[2], Vjm);

  coord_ijk(i, j+1, k, locarray[2], Xjp1); // ok to have j+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i, j+1, k, locarray[2], Vjp1);
    
  coord_ijk(i, j, k, locarray[3], Xkm);
  bl_coord_ijk(i, j, k, locarray[3], Vkm);

  coord_ijk(i, j, k+1, locarray[3], Xkp1); // ok to have k+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i, j, k+1, locarray[3], Vkp1);
    
  dr = THIRD*(pow(Vip1[1],3.0)-pow(Vim[1],3.0)); // really \delta(r^3/3) = r^2 dr
  // dh doesn't reduce to Pi, but 2.0 that is the correct answer for any N2 resolution
  dh = -(cos(Vjp1[2])-cos(Vjm[2])); // really \delta(-cos(h)) = sinh dh
  // below is inconsistent with rest of code when N2!=1
  //    dh = (N2==1) ? M_PI : -(cos(Vjp1[2])-cos(Vjm[2])); // really \delta(-cos(h)) = sinh dh
  dp = (Vkp1[3]-Vkm[3]); // just d\phi 
    
  // finite volume of cell
  newjacvol = dr*dh*dp/(dx[1]*dx[2]*dx[3]);


  ////////////////////////
  //
  // OLD JAC
  drold=Vc[1]*Vc[1]*dxdxpc[1][1]*dx[1]; // r^2 dr
  if(totalsize[2]==1) dhold=2.0;
  else dhold=sin(Vc[2])*dxdxpc[2][2]*dx[2]; // sinh dh
  //dhold=sin(Vc[2])*dxdxpc[2][2]*dx[2]; // sinh dh (oldjacvol consistent with gdettrue)
  dpold=2.0*M_PI*dx[3];
  oldjacvol = drold*dhold*dpold/(dx[1]*dx[2]*dx[3]); // only true if diagonal mapping

  // suppose gdettrue already corrected if wanted to be corrected
  //  if(totalsize[2]==1) (*gdettrue) /= (M_PI*0.5); // correct eomfunc and gdettrue
  //  if(WHICHEOM==WITHGDET) PLOOP(pliter,pl) EOMFUNCMAC(pl)=(*gdettrue); // else up to user to make sure right

  //////////////////////////
  //
  // use below, works best.  Only central conserved quantity operated on (source terms too)
  if(loc==CENT){
    PLOOP(pliter,pl) EOMFUNCASSIGN(pl)=newjacvol;
    *gdetvol=(*gdettrue)=newjacvol;
  }
  else *gdetvol=(*gdettrue);


  ///////////////////////////////////
  //
  // SOME OTHER ATTEMPTS for 2 lines "best" above:
  //
  //////////////////////////////////

  // horrible
  //  if(loc==CENT) *gdetvol=newjacvol;
  //  else *gdetvol=(*gdettrue);


  // central region undershoots
  //if(loc==CENT) *gdetvol=newjacvol;
  //  else *gdetvol=(*gdettrue);

  // big overshoot, so probably don't want to use for fluxes
  //  PLOOP(pliter,pl) *gdetvol=(*gdettrue)=EOMFUNCMAC(pl)=newjacvol;

  // with newjac in potential and gdettrue here
  // little bit more of a spike at center
  //  *gdetvol=(*gdettrue);

  // with oldjac in potential and gdettrue here
  // same as above :  little bit more of a spike at center
  //  *gdetvol=(*gdettrue);


  //  *gdetvol = newjacvol;
  //    *gdetvol = dr*dh*dp/5.0;
  //*gdetvol = (*gdettrue); // disable gdetvol but keep memory
  //  PLOOP(pliter,pl) EOMFUNCMAC(pl)=(*gdettrue)=*gdetvol=newjacvol;


  //  *gdetvol=newjacvol;
  //  if(loc==CENT){
  //  PLOOP(pliter,pl) EOMFUNCMAC(pl)=(*gdettrue)=*gdetvol=newjacvol;
  //  PLOOP(pliter,pl) EOMFUNCMAC(pl)=(*gdettrue)=*gdetvol=oldjacvol;
  //  PLOOP(pliter,pl) EOMFUNCMAC(pl)=*gdetvol=(*gdettrue);
  //  }
  //  else{
  //    PLOOP(pliter,pl) EOMFUNCMAC(pl)=*gdetvol=(*gdettrue)=oldjacvol;
  //PLOOP(pliter,pl) EOMFUNCMAC(pl)=*gdetvol=(*gdettrue)=newjacvol;
  //  }


  //  if(loc==CENT){
  //    dualfprintf(fail_file,"i=%d loc=%d :: gdettrue=%15.7g oldjac=%15.7g newjac=%15.7g :: dr=%g drold=%g :: dh=%g dhold=%g :: dp=%g dpold=%g :: %g %g %g\n",i,loc,(*gdettrue),oldjacvol,newjacvol,dr,drold,dh,dhold,dp,dpold,dx[1],dx[2],dx[3]);
  //  }

  //(*gdettrue) = *gdetvol; // replace gdettrue with this version


}

















