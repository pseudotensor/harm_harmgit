
/*! \file interppoint.para.c
     \brief Parabolic/PPM Spatial Interpolation for fluxes based upon providing each point
     /// parabolic interpolation subroutin  
     /// ref. Colella && Woodward's paper
     /// Colella, P., & Woodward, P. R. 1984, J. Comput. Phys., 54, 174-201
     ///
     /// using zone-centered value of 5 continuous zones 
     /// to get left and right value of the middle zone.
     ///  
     /// 
     ///
     // Latest JCM version is para4()/parapl(): 02/25/08

*/


#include "decs.h"

#include "para_and_paraenohybrid.h"


/// initial checks regarding usability of para method
void parainitchecks(void)
{
  int dimen;


  DIMENLOOP(dimen){
    if(PARAGENDQALLOWEXTREMUM && WENOINTERPTYPE(lim[dimen])){
      // if PARAMODWENO==1 implied, but if 0 allow since doesn't matter if doing weno with PARAGENDQALLOWEXTREMUM==1
      // then ok
    }
    else if(PARAMODWENO==1 && PARAGENDQALLOWEXTREMUM==0 && WENOINTERPTYPE(lim[dimen])){
      dualfprintf(fail_file,"WARNING: PARAGENDQALLOWEXTREMUM==0 for hybrid method, will be less accurate for turbulent regions\n");
    }
    else if(PARAGENDQALLOWEXTREMUM==1 && lim[dimen]==PARA){
      dualfprintf(fail_file,"ERROR: PARAGENDQALLOWEXTREMUM==1 but not WENO-based limiter -- will be less stable in stiff regions (e.g. near horizon). Use PARALINE instead.\n");
      myexit(34643463);
    }

    if(lim[dimen]==PARAFLAT){
      dualfprintf(fail_file,"WARNING: PARAFLAT inefficient compared to PARALINE, suggested to use PARALINE instead\n");
    }
  }


}




/// lout/rout is left and right sides of cell
/// note how used in step_ch.c to get appropriate interface value
/// given by Xiaoyue Guan to Scott Noble on Nov 9, 2004, given to me Jan 7, 2005
void para(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}



/// given by Xiaoyue Guan on Jan 9, 2005
void para2(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == PMC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    dq[mm] =Dqc;
#endif
  }
  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  /*
    l=max(min(y[0],y[-1]),l);
    l=min(max(y[0],y[-1]),l);
    r=max(min(y[0],y[1]),r);
    r=min(max(y[0],y[1]),r);
  */

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}




/// 3rd para from Xiaoyue that she bundled with a new TVD-optimal RK3
/// given on 02/17/2005
void para3(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    else dq[mm] = -aDqvanl*sign(Dqc);
    //else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=-min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = -aDqm*sign(Dqc);
    else dq[mm]=-aDqp*sign(Dqc);
#elif(PARA2LIM == NLIM) //w/o slope limiter
    //if(s<=0.) dq[mm] = 0.; // DONOR
    dq[mm] = Dqc;
#endif
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));

  /*
    if (qa <=0. ) {
    l=y[0];
    r=y[0];
    }

    else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
    else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


    *lout=l;   //a_L,j
    *rout=r;
    */

  if (qa <=0. ) {
    *lout=y[0];
    *rout=y[0];
  }  
  else {
    *lout = l;
    *rout = r;
  }

  //2.0 at top/bottom of a steep gradient 
  if (qd*(qd-qe)<0.0) *lout=3.0*y[0]-2.0*r;
  else *lout = l;

  if (qd*(qd+qe)<0.0) *rout=3.0*y[0]-2.0*l;
  else *rout = r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}





/// older method that has no failures problem CHANGINGMARK
void para4_old(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;


#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    //else dq[mm] = aDqvanl*sign(Dqc);
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

#if(0)
    // Jon's version
    dq[mm]=MINMOD(Dqc,MINMOD(Dqm,Dqp));
#else
    // Xioyue's version
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#endif



#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);


#elif(PARA2LIM == MINM) // no steepener, normal MINM

#if(0)
    // Jon's version
    dq[mm] = MINMOD(0.5*Dqm,0.5*Dqp); // no steepening    
#elif(1)
    // Jon's steep version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);
#elif(0)
    // Xioyue's version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = 0.5*aDqm*sign(Dqc);
    else dq[mm]=0.5*aDqp*sign(Dqc);
#endif

#elif(PARA2LIM == NLIM) //w/o slope limiter

    dq[mm] = Dqc;
#endif
  }

#if(JONPARAREDUCE)
  //  if(pl==U1){
  if(pl!=RHO){
    if(
       (fabs(dq[-1]-dq[0])/(fabs(dq[-1])+fabs(dq[0])+SMALL)>0.1)||
       (fabs(dq[1]-dq[0])/(fabs(dq[1])+fabs(dq[0])+SMALL)>0.1)
       ){
      slope_lim_3points(MINM, y[-1], y[0], y[1], dq);
      *lout =y[0] - 0.5* (*dq);
      *rout=y[0] + 0.5* (*dq);
      return;
    }
  }

#endif

  /* CW1.6 */

  // modified as per Matt's paper
  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/8.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/8.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  // modified as per Matt's paper
  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }
  else{

    if (qd*(qd-qe)<0.0) 
      l=3.0*y[0]-2.0*r;


    if (qd*(qd+qe)<0.0) 
      r=3.0*y[0]-2.0*l;
  }


  *lout=l;   //a_L,j
  *rout=r;

  //  *dqleft=dq[-1];
  //  *dqcenter=dq[0];
  //  *dqright=dq[1];

}


#define DO4MONO 0 // whether to add-back-in some trend that one removed -- not working yet in this code

#define OVERRIDEWITHMINM 1

/// Jon's MINM/VANL/MC with parabolic interpolation if 3 dq's would be chosen as centered slope
/// Useful to use cour=0.45 < 0.5 to keep things stable as well
void parajon(int ii, int jj, int kk, int loc, int realisinterp, int dir, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm;
  FTYPE dq0[10];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void jonparasmooth_compute(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *dq1, FTYPE *dq2, FTYPE *lout, FTYPE *rout, int *smooth);
  int smooth=0;
  int a_whichdq[10];
  FTYPE a_y4mono[10];
  FTYPE *y4mono;
  int *whichdq;
  int dqrange=2;
  int usepara;
  FTYPE s0;
  int iii,jjj,kkk;
  FTYPE Dqm4mono,Dqc4mono,Dqp4mono;
  FTYPE s4mono;
  FTYPE aDqm4mono,aDqc4mono,aDqp4mono;
  
  int shifti=(dir==1);
  int shiftj=(dir==2);
  int shiftk=(dir==3);
  

  // shifted dq
  dq=dq0+dqrange;
  whichdq=a_whichdq+dqrange;

#if(DO4MONO)
  y4mono=a_y4mono+dqrange;

  // get true function
  for(mm=-2 ; mm<=2 ; mm++) {
    iii=ii+shifti*mm;
    jjj=jj+shiftj*mm;
    kkk=kk+shiftk*mm;

    y4mono[mm] = y[mm]*GLOBALMACP1A0(pother,RHOCENTA+pl,iii,jjj,kkk);
  }
#endif


  /*CW1.7 */
  s0=1.0;
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);


#if(DO4MONO)
    // 4mono version
    Dqm4mono = 2.0 *(y4mono[mm]-y4mono[mm-1]);
    Dqp4mono = 2.0 *(y4mono[mm+1]-y4mono[mm]);
    Dqc4mono = 0.5 *(y4mono[mm+1]-y4mono[mm-1]);
#else
    // normal version
    Dqm4mono = Dqm;
    Dqp4mono = Dqp;
    Dqc4mono = Dqc;
#endif


    // check what slope is chosen
    if(fabs(Dqc4mono)<=fabs(Dqm4mono) && fabs(Dqc4mono)<=fabs(Dqp4mono)){
      // central picked .. is fine
      whichdq[mm]=0;
    }
    else if(fabs(Dqm4mono)<=fabs(Dqp4mono) && fabs(Dqm4mono)<=fabs(Dqc4mono)){
      // check if really want to steepen
      whichdq[mm]=-1;
    }
    else if(fabs(Dqp4mono)<=fabs(Dqm4mono) && fabs(Dqp4mono)<=fabs(Dqc4mono)){
      whichdq[mm]=+1;
    }
    else{
      whichdq[mm]=+1; // GODMARK: Shouldn't really get here
    }
  



    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    aDqm4mono = fabs(Dqm4mono) ;
    aDqp4mono = fabs(Dqp4mono) ;
    aDqc4mono = fabs(Dqc4mono) ;

    s4mono = Dqm4mono*Dqp4mono;
    if(s4mono<=0) s0=-1;



#if(OVERRIDEWITHMINM || PARA2LIM == MINM) // no steepener, normal MINM
    // in stiff regime, MINM is more stable than other schemes with steepeners
    if (s4mono<=0.) dq[mm] = 0.;
    else if (aDqm4mono<aDqp4mono) dq[mm] = 0.5*Dqm;
    else dq[mm]=0.5*Dqp;
#elif(PARA2LIM == VANL) 
    // NOT DONE for 4mono method
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s4mono <=0.) dq[mm]=0.;
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

    if (s4mono <=0.) dq[mm]=0.;       //CW1.8
    else if (aDqm4mono<aDqp4mono && aDqm4mono<aDqc4mono) dq[mm] = Dqm;
    else if (aDqp4mono<aDqm4mono && aDqp4mono<aDqc4mono) dq[mm] = Dqp;
    else dq[mm]=Dqc;

#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    if (s4mono<=0.) dq[mm] = 0.;
    else if (aDqm4mono<aDqp4mono) dq[mm] = Dqm;
    else dq[mm]=Dqp;


#elif(PARA2LIM == NLIM) //w/o slope limiter
    dq[mm] = Dqc;
#endif
  }




  // no matter what dqrange is, just go from -1..1  
  usepara=0;
  for(mm=-1 ; mm<=1 ; mm++) {
    usepara+=(whichdq[mm]==0);
  }

  //  if(usepara==3 && s0>=0.0){
  //  if(usepara==3 && s0>=0.0 && 0){ // ONLY  MINMOD
  if(usepara==3 && s0>=0.0){
    // then use interior parabola and be done
    *lout=(1./8.)*(6.0*y[0]+3.0*y[-1]-y[1]);
    *rout=(1./8.)*(6.0*y[0]-y[-1]+3.0*y[1]);
  }
  else{
    // then just use 2nd order
    *lout=y[0]-0.5*dq[0];
    *rout=y[0]+0.5*dq[0];
  }


}




/// doesn't use dq
void para4(int realisinterp, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  void para4gen(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout, FTYPE *dq, int *smooth);
  void checkparamonotonicity(int smooth, int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout);
  FTYPE dq0[5];
  FTYPE *dq;
  int dqrange;
  FTYPE a_ddq[7];
  FTYPE *ddq;
  int mm;
  int smooth;


  dqrange=3;
  // shifted dq
  dq=dq0+2;
  ddq=a_ddq+3; // shifted sufficiently


  para4gen(realisinterp,dqrange, pl,y,lout,rout,dq,&smooth);

  for(mm=-dqrange/2+1;mm<=dqrange/2;mm++){
    ddq[mm] = dq[mm] - dq[mm-1];
  }
  checkparamonotonicity(smooth, dqrange, pl, y, ddq, dq, lout, rout, lout, rout);


}




#if(PARAGENDQALLOWEXTREMUM==0)
#define PARAGENMINMOD(a,b) MINMOD(a,b)
#else
#define PARAGENMINMOD(a,b) MINMODB(a,b)
#endif

/// Xiaoyue given on 03/25/05
/// she realized sign error in 1st der's in para3()
/// noted Matt's paper astro-ph/0503420 suggested CW1.6 uses 1/8 rather than 1/6
/// I noticed that Matt uses MC for field variables and PPM+ for hydro variables
/// GODMARK: This step could be done once for entire line of data only once
void para4gen(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout, FTYPE *dq, int *smooth)
{
  int mm ;
  FTYPE a_dq1[10],a_dq2[10];
  FTYPE *dq1,*dq2;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  FTYPE Dqparacenterleft,Dqparacenterright;
  FTYPE Dqsteep;
  FTYPE ddql,ddqr;
  void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont);
  void jonparasmooth_compute(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *dq1, FTYPE *dq2, FTYPE *lout, FTYPE *rout, int *smooth);





  // Procedure:
  // 1) First obtain limited slopes (above)
  // 2) Assume l,r obtained from 3rd order polynomial for points or quartic for averages
  // 3) steepen or flatten solution
  // 4) Check monotonicity: 
  //    Condition 1: If y[0] is a local max or local min compared to l,r, then set l=r=y[0]
  //    Condition 2: If parabola fit through l,y[0],r is non-monotonic, then if non-monotonic region is near l, then modify r until derivative of parabola at l is 0.

  // shifted dq (sufficiently shifted)
  dq1=a_dq1+dqrange;
  dq2=a_dq2+dqrange;
  // setup defaults
  *smooth=0;



#if(JONPARASTEEP)
  /////////////////
  //
  // first check if 3rd order polynomial will have derivative at center of cell that MC says doesn't need steepening
  //
  ////////////////

  mm=0;
  Dqm = 2.0 *(y[mm]-y[mm-1]);   // steepened
  Dqp = 2.0 *(y[mm+1]-y[mm]);   // steepened
  //  Dqc = 0.5 *(y[mm+1]-y[mm-1]); // normal
  // Note that the factors of 3 (and 6) cause asymmetries at roundoff-error
  Dqparacenterleft =  (+0.5*y[0]-y[-1]    +y[-2]/6.0+y[1]/3.0); // used to get a_{j-1/2}
  Dqparacenterright = (-0.5*y[0]-y[-1]/3.0+y[1] -    y[2]/6.0); // used to get a_{j+1/2}
  // see if para will use centered slope that is not steepened enough compared to MC

  // first compare Dqp and Dqm
  Dqsteep=PARAGENMINMOD(Dqm,Dqp);
  //    *dq=PARAGENMINMOD(Dqc,PARAGENMINMOD(Dqm,Dqp));
  // now compare centered with steepened Dqsteep value
  if(fabs(Dqparacenterleft)>fabs(Dqsteep) && fabs(Dqparacenterright)>fabs(Dqsteep)){
    *lout = y[0] - 0.5* Dqsteep;
    *rout = y[0] + 0.5* Dqsteep;
    //dualfprintf(fail_file,"USEDSTEEP\n");
    return;
    // else PARA is ok to use
  }
#endif

  //dualfprintf(fail_file,"USEDPARA\n");



 
  ////////////////////////
  //
  // get slopes
  //
  ////////////////////////

  // Dqm(0) = 2.0*(dq1[0])
  // Dqp(0) = 2.0*(dq1[1]) // so need to go 1 farther to get all needed Dqp's
  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {
    dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face
    dq2[mm] = 0.5 *(y[mm+1]-y[mm-1]); // slope centered at cell center
  }
  mm=dqrange/2+1; // get last dq1 (can't do in loop above since +1 would mean dq2 beyond data range
  dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face
  

  /////////////////
  //
  // Determine monotonized slopes
  //
  /////////////////

  /*CW1.7 */
  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {

    Dqm = 2.0 *dq1[mm];   // steepened
    Dqp = 2.0 *dq1[mm+1]; // steepened
    Dqc = dq2[mm];        // normal


#if(PARA2LIM == VANL) 
    Dqm*=0.5; // desteepen
    Dqp*=0.5; // desteepen
    s = Dqm*Dqp;
    Dqvanl=2.0*s/(Dqm+Dqp);
    //    aDqvanl=fabs(Dqvanl);
    //    aDqm = fabs(Dqm) ;
    //    aDqp = fabs(Dqp) ;
    //    aDqc = fabs(Dqc) ;

    // true VANL:
    if (s <=0.) dq[mm]=0.;
    else dq[mm] = Dqvanl*sign(Dqc);

    // Xioyue was using MC-VANL combo of some kind
    //    else dq[mm] = aDqvanl*sign(Dqc);
    //    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

#if(1)
    // Jon's version
    dq[mm]=PARAGENMINMOD(Dqc,PARAGENMINMOD(Dqm,Dqp));
#else
    // Xioyue's version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#endif



#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);


#elif(PARA2LIM == MINM) // no steepener, normal MINM

#if(1)
    // Jon's version
    dq[mm] = PARAGENMINMOD(0.5*Dqm,0.5*Dqp); // no steepening    
#elif(0)
    // Jon's steep version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);
#elif(0)
    // Xioyue's version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = 0.5*aDqm*sign(Dqc);
    else dq[mm]=0.5*aDqp*sign(Dqc);
#endif

#elif(PARA2LIM == NLIM) //w/o slope limiter

    dq[mm] = Dqc;
#endif
  }// end loop over mm's




#if(JONPARASMOOTH)
  jonparasmooth_compute(realisinterp,dqrange,pl,y,dq1,dq2,lout,rout,smooth);
  if(*smooth) return;
#else
  *smooth=0;
#endif







#if(JONPARAREDUCE)
  //  if(pl==U1){
  if(pl!=RHO){
    if(
       (fabs(dq[-1]-dq[0])/(fabs(dq[-1])+fabs(dq[0])+SMALL)>0.1)||
       (fabs(dq[1]-dq[0])/(fabs(dq[1])+fabs(dq[0])+SMALL)>0.1)
       ){
      slope_lim_3points(MINM, y[-1], y[0], y[1], dq);
      *lout =y[0] - 0.5* (*dq);
      *rout=y[0] + 0.5* (*dq);
      return;
    }
  }

#endif



#if(WHICHLIMITERTOUSEFORLR==0)
  //////////////////////////////////
  //
  // Obtain continuous solution at interface
  //
  // CW1.6 for obtaining a_{j+1/2} using quartic polynomial, but slopes have been replaced with limited slopes
  //
  // GODMARK: This step could be done once for entire line of data
  
  ddql=(dq[0]-dq[-1]);
  ddqr=(dq[1]-dq[0]);
  
  paracont(ddql, &y[0], &l);
  paracont(ddqr, &y[1], &r);





#if(0)
  // doesn't seem to be within CW!
  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);
#endif

  *lout = l;
  *rout = r;



#elif(WHICHLIMITERTOUSEFORLR==1)


  // notice that this results in discontinuous solution at each interface already, unlike para method
  *lout = y[0] - 0.5*dq[0];
  *rout = y[0] + 0.5*dq[0];

  // differs from normal MC,MINM, etc. by computing dq's that will be used later for steepener, fattener, and parabolic check

  // Note that MC is actually 3rd order acccurate (error term O(dx^3)) because linear l,r gives same answer as parabola!


#endif





}




/// check if solution locally without discontinuities as indicated by use of central slope within -1,0,+1 for non-RHO and -2..+2 for RHO (helps moving contact)
/// Required for PARA since PARA otherwise is too non-diffusive when reducing to central slopes and activating 3rd order polynomial at faces that has no jumps and no required diffusion for stability.  PARA otherwise has unstable noise moving supersonically in caustic problem
/// Seems to work in general even without special extra smooth switching!
void jonparasmooth_compute(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *dq1, FTYPE *dq2, FTYPE *lout, FTYPE *rout, int *smooth)
{
  int usepara=0;
  int mm;
  FTYPE Dqm,Dqp,Dqc;
  int a_whichdq[10];
  int *whichdq;


  whichdq=a_whichdq+dqrange;


  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {

    Dqm = 2.0 *dq1[mm];   // steepened
    Dqp = 2.0 *dq1[mm+1]; // steepened
    Dqc = dq2[mm];        // normal

    if(fabs(Dqc)<=fabs(Dqm) && fabs(Dqc)<=fabs(Dqp)){
      // central picked .. is fine
      whichdq[mm]=0;
    }
    else if(fabs(Dqm)<=fabs(Dqp) && fabs(Dqm)<=fabs(Dqc)){
      // check if really want to steepen
      whichdq[mm]=-1;
    }
    else if(fabs(Dqp)<=fabs(Dqm) && fabs(Dqp)<=fabs(Dqc)){
      whichdq[mm]=+1;
    }
    else{
      whichdq[mm]=+1; // GODMARK: Shouldn't really get here
    }
  }


  if(realisinterp==1 && pl==RHO && dqrange>=5){
    // check steepened slopes to ensure really want to steepen
    //  for(mm=-dqrange/2+1 ; mm<=dqrange/2-1 ; mm++) {
    for(mm=-2 ; mm<=2 ; mm++) {
      usepara+=(whichdq[mm]==0);
    }

    if(usepara==5){
      // then use interior parabola and be done
      *lout=(1./8.)*(6.0*y[0]+3.0*y[-1]-y[1]);
      *rout=(1./8.)*(6.0*y[0]-y[-1]+3.0*y[1]);
    
      *smooth=1;
      return;
    }
    // else use para as normal
  }
  else{
    // check steepened slopes to ensure really want to steepen
    //  for(mm=-dqrange/2+1 ; mm<=dqrange/2-1 ; mm++) {
    for(mm=-1 ; mm<=1 ; mm++) {
      usepara+=(whichdq[mm]==0);
    }

    if(usepara==3){
      // then use interior parabola and be done
      *lout=(1./8.)*(6.0*y[0]+3.0*y[-1]-y[1]);
      *rout=(1./8.)*(6.0*y[0]-y[-1]+3.0*y[1]);
    
      *smooth=1;
      return;
    }
    // else use para as normal
  }

  *smooth=0;

}






/// not used right now
void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont)
{
  FTYPE avgpointcoef;


#if(AVGINPUT)
  // /6 is result of passing through points from differential of average values in each cell
  avgpointcoef=1.0/6.0;
  *facecont=0.5*(y[0]+y[-1])-ddq*avgpointcoef;
#else
  avgpointcoef=1.0/8.0;
  // consistent with Matt's paper
  // assume passing quartic polynomial through points y[-2,-1,0,1,2]
  // see PARA_interpolation_checks.nb
  *facecont=0.5*(y[0]+y[-1])-ddq*avgpointcoef;
#endif

}




/// PARA starts with continuous 3rd order polynomial (4th order through averages if used) using monotonized slopes for second derivative term
/// facecont[i] = 0.5*(y[i]+y[i-1]) - ddq[i]   where ddq[0] = dq[0]-dq[-1]
void paracontsmooth(int pl, FTYPE *y, FTYPE *facecont, int *smooth)
{
  FTYPE ddqtest[4];
  int whichmax;
  FTYPE maxddq;
  FTYPE qc,qd,qe;
  int i;

#if(AVGINPUT)
  dualfprintf(fail_file,"NOT SETUP\n");
  myexit(249678346);
#else
  // PARAFLAT has enough points for this
  ddqtest[0]=y[-1]-2*y[-2]+y[-3];
  ddqtest[1]=y[0]-2*y[-1]+y[-2];
  ddqtest[2]=y[1]-2*y[0]+y[-1];
  ddqtest[3]=y[2]-2*y[1]+y[0];
  *smooth=0;

  // get biggest ddq
  maxddq=0.0;
  for(i=0;i<=3;i++){
    if(fabs(ddqtest[i])>maxddq){
      maxddq=ddqtest[i];
      whichmax=i;
    }
  }
  
  // normalize ddq's
  for(i=0;i<=3;i++){
    ddqtest[i]=ddqtest[i]/maxddq; // Then 1.0 will be that ddq that is max
  }

  // value of ddq should be same as largest ddq within some percent error

  *smooth=1;
  for(i=0;i<=3;i++){
    if(ddqtest[i]<-0.1) *smooth=0;
  }

  if(pl==U1){
    if(!*smooth){
      dualfprintf(fail_file,"NOTSMOOTH\n");
      for(i=0;i<=3;i++) dualfprintf(fail_file,"ddqtest[%d]=%g maxddq=%g whichmax=%d\n",i,ddqtest[i],maxddq,whichmax);
    }
  }
    

  //  if( (sign(ddqtest[0])==sign(ddqtest[1]) && sign(ddqtest[1])==sign(ddqtest[2]) && sign(ddqtest[2])==sign(ddqtest[3])){

  // always compute in case used later regardless of smoothness measurement
  *facecont = (1./16.)*(9.*y[0]+9*y[-1]-y[-2]-y[1]);
  //    *smooth=1;
  //}

#if(0)
  if(sign(ddqtest[0])==sign(ddqtest[1]) && sign(ddqtest[1])==sign(ddqtest[2]) && sign(ddqtest[2])==sign(ddqtest[3])){
    if(pl==2){
      dualfprintf(fail_file,"GOOD :: %g %g %g %g :: %d %d %d\n",
                  sign(ddqtest[0]),sign(ddqtest[1]),sign(ddqtest[2]),sign(ddqtest[3])
                  ,sign(ddqtest[0])==sign(ddqtest[1]),sign(ddqtest[1])==sign(ddqtest[2]),sign(ddqtest[2])==sign(ddqtest[3])
                  );
    }
  }
  else{
    if(pl==2){
      dualfprintf(fail_file,"BAD::: %g %g %g %g :: %d %d %d\n",
                  sign(ddqtest[0]),sign(ddqtest[1]),sign(ddqtest[2]),sign(ddqtest[3])
                  ,sign(ddqtest[0])==sign(ddqtest[1]),sign(ddqtest[1])==sign(ddqtest[2]),sign(ddqtest[2])==sign(ddqtest[3])
                  );
    }
  }
#endif




#endif



}




/// check monotonicity of parabola
void checkparamonotonicity(int smooth, int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout)
{
  FTYPE a6COEF;
  FTYPE qa,qb,qd,qe;
  FTYPE r,l;
  int i;
  int numddq;



  

  l = *lin;
  r = *rin;


#if(PARAGENDQALLOWEXTREMUM)

  if(smooth) qb=1;
  else{
    // (dq[0]-dq[-1]) is defined as ddq[0]
    //  ddqr=(dq[1]-dq[0]) is ddq[1]
    qb=0.0;
    numddq=0;
    for(i=-dqrange/2+1;i<=dqrange/2;i++){
      qb+=sign(ddq[i]);
      numddq++;
    }
    // check that ddq's all same sign
    if(fabs(qb)>(FTYPE)numddq-0.1) qb=1.0; // 0.1 is just to avoid machine precision issue
    else qb=-1.0;
  }
#else
  if (smooth){
    *lout=*lin;
    *rout=*rin;
    return;
  }
#endif


  /////////////
  //
  // now perform the limiting procedures that create the discontinuities

  // Condition 1: monotonicity for l,y[0],r. qa>0 if monotonic
  qa=(r-y[0])*(y[0]-l);
  // modify Condition 1 as in Duez et al. (2005)
  // allow nonmonotonic behavior if the second derivative doesn't change sign around the center
  // GODMARK: actually should check 2 derivatives in each direction!  Then need PARAFLAT for extra values


  // CW: \delta a_j
  qd=(r-l);

#if(AVGINPUT)
  // CW: a_{6,j} in CW, which is for averages:
  a6COEF=6.0;
  qe=a6COEF*(y[0]-0.5*(l+r));
#else
  // for points need: see PARA_interpolation_checks.nb
  // parabola has solution y = l + (x-xl)*(qd + a6*(1-(x-xl))) with a6 = 4(y_0-1/2(l+r)) for points
  a6COEF=4.0;
  qe=a6COEF*(y[0]-0.5*(l+r));
#endif





#if(PARAGENDQALLOWEXTREMUM)
  if (qa <=0.0 && qb<=0.0 )
#else
    if (qa <=0.0)
#endif
      { // Condition 1


#if(NONMONOLIM==0)
        l=y[0];
        r=y[0];
#elif(NONMONOLIM==1)
        // makes no sense to reduce all the way to DONOR since to second order can still have monotonic result, so use MONO result in this case and assume flatten result
        // appears to be too speculative as results in  more failures at horizon with PARA2LIM==MC
        l = y[0] - 0.5* dq[0];
        r = y[0] + 0.5* dq[0];
#endif

      }
    else if(1){
      // Condition 2

      // qe can be positive or negative still here even though qa>0
      if     (qd*(qd-qe)<0.0)  l = (-(2.0+a6COEF)*r + 2.0*a6COEF*y[0])/(a6COEF-2.0);
      else if(qd*(qd+qe)<0.0)  r = (-(2.0+a6COEF)*l + 2.0*a6COEF*y[0])/(a6COEF-2.0);
      // else no change needed
      //    else{
      //      dualfprintf(fail_file,"Problem with limiting condition 2\n");
      //    }
    
      // Xiaoyue's verison:
      //    if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
      //    if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;
    }


#if(PARAGENDQALLOWEXTREMUM)
  // recompute monotonicity to confirm didn't screw it up:
  qa=(r-y[0])*(y[0]-l);
  // qb is not changed!
  if (qa <=0. && qb<=0.0 ) { // Condition 1 again
    // with Duez allowance of non-monotonic behavior, this gets triggered if l,r change -- otherwise wouldn't have been triggered
    l=y[0];
    r=y[0];
  }
#endif


  // assign output
  *lout=l;
  *rout=r;



}







/// used when lim=PARAFLAT
void parapl(int i, int j, int k, int loc, int realisinterp, int dir, FTYPE **yrealpl, FTYPE **ypl, FTYPE *loutpl, FTYPE *routpl)
{
  FTYPE dq0[NPR2INTERP][8];
  FTYPE *dq[NPR2INTERP];
  FTYPE *y,*yreal;
  void parasteep(int dir, int pl, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r);
  void getPressure(int whicheom, int i, int j, int k, int loc, FTYPE **yrealpl, FTYPE *P);
  FTYPE a_P[NUMTRUEEOMSETS][10];
  FTYPE *V[NUMTRUEEOMSETS],*P[NUMTRUEEOMSETS];
  FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P);
  int pl,pliter;
  FTYPE Fi[NUMTRUEEOMSETS];
  int dqrange;
  FTYPE a_ddq[7];
  FTYPE *ddq;
  int mm;
  int smooth;
  int whicheom;



  // consistent with PARAFLAT using 7 points
  dqrange = 5; // dq's will exist from -2,-1,0,1,2 and ddq computed from -2,-1,0,1

  // shift P
  for(whicheom=0;whicheom<NUMTRUEEOMSETS;whicheom++){
    P[whicheom]=a_P[whicheom] + 4; // P accessed from -3..3 ( shifted sufficiently)
  }

  // shift dq
  PINTERPLOOP(pliter,pl){
    dq[pl]=dq0[pl]+4; // shifted sufficiently
  }

  ddq=a_ddq+3; // shifted sufficiently


  // assume velocity is istelf
  // KORALTODO: need Ficalc for radiation by itself!
  V[EOMSETMHD] = yrealpl[U1+dir-1];
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
  V[EOMSETRAD] = yrealpl[URAD1+dir-1];
#endif


  // get pressures for all points since needed for reduction or steepening
#if( DOPPMREDUCE || DOPPMCONTACTSTEEP)
  for(whicheom=0;whicheom<NUMTRUEEOMSETS;whicheom++){
    if(whicheom==EOMSETRAD && RADSHOCKFLAT || whicheom==EOMSETMHD) getPressure(whicheom,i, j, k, loc, yrealpl, P[whicheom]);
  }
#endif



  // computed only once for all variables
  for(whicheom=0;whicheom<NUMTRUEEOMSETS;whicheom++){
    if(whicheom==EOMSETRAD && RADSHOCKFLAT || whicheom==EOMSETMHD){
#if( DOPPMREDUCE )
      Fi[whicheom] = Ficalc(dir,V[whicheom],P[whicheom]);
#else
      Fi[whicheom] = 0.0;
#endif
    }
  }


  ///////////////
  //
  // Loop over variables and get interpolated left/right values within cell
  //
  //////////////


  PINTERPLOOP(pliter,pl){

    y=ypl[pl];
    yreal=yrealpl[pl];

    // get continuous solution    
    para4gen(realisinterp,dqrange,pl,y,&loutpl[pl],&routpl[pl],dq[pl],&smooth);

#if(DOPPMCONTACTSTEEP)
    if(RADFULLPL(pl)&&RADSHOCKFLAT) parasteep(dir,pl,V[EOMSETRAD],P[EOMSETRAD],ypl[pl],dq[pl],&loutpl[pl],&routpl[pl]);
    if(!RADFULLPL(pl)) parasteep(dir,pl,V[EOMSETMHD],P[EOMSETMHD],ypl[pl],dq[pl],&loutpl[pl],&routpl[pl]);
#endif


#if( DOPPMREDUCE )
    if(RADFULLPL(pl)&&RADSHOCKFLAT) paraflatten(dir,pl,ypl[pl],Fi[EOMSETRAD],&loutpl[pl],&routpl[pl]);
    if(!RADFULLPL(pl)) paraflatten(dir,pl,ypl[pl],Fi[EOMSETMHD],&loutpl[pl],&routpl[pl]);
#endif


    // finally check monotonicity of the parabola and create discontinuities if non-monotonic
    // FLASH equations 51 -> 53 for points or averages
    for(mm=-dqrange/2+1;mm<=dqrange/2;mm++){
      ddq[mm] = dq[pl][mm] - dq[pl][mm-1];
    }
    checkparamonotonicity(smooth, dqrange, pl, ypl[pl], ddq, dq[pl], &loutpl[pl], &routpl[pl], &loutpl[pl], &routpl[pl]);

#if(NONMONOLIM>0 && DOPPMREDUCE)
    // then flatten again
    if(RADFULLPL(pl)&&RADSHOCKFLAT) paraflatten(dir,pl,ypl[pl],Fi[EOMSETRAD],&loutpl[pl],&routpl[pl]);
    if(!RADFULLPL(pl)) paraflatten(dir,pl,ypl[pl],Fi[EOMSETMHD],&loutpl[pl],&routpl[pl]);
#endif

    
  }



}



/// PPM FLATTENER formula
void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r)
{
  // FLASH Equation 49,50
  *l = Fi * y[0] + ( 1.0 - Fi ) * (*l);
  *r = Fi * y[0] + ( 1.0 - Fi ) * (*r);
}




/// PPM FLATTENER parameter
FTYPE ftilde( int dir, int shift, FTYPE *Vabs, FTYPE *Pabs)
{
  FTYPE Ftilde,Ftilde1,Ftilde2;
  FTYPE Sp;
  FTYPE *V, *P;
  FTYPE P2diff,Pbottom;


  // shift as needed
  P = Pabs + shift;
  V = Vabs + shift;

  // FLASH Equation 43 (pressure jump)
  P2diff=P[2]-P[-2];
  Pbottom=sign(P2diff)/(fabs(P2diff)+SMALL); // singularity avoidance but keeps signature
  Sp = (P[1] - P[-1]) * Pbottom ;




  // FLASH Equation 45 (shock must have sufficient pressure jump)
  Ftilde = max( 0, min( 1.0, 10.0 * (Sp - SP0) ) );

  // FLASH Equation 46 (shock must have pressure jump)
  Ftilde1 = fabs(P[1] - P[-1]) / (min(fabs(P[1]), fabs(P[-1]))+ SMALL );
  Ftilde *= ( (FTYPE)(Ftilde1>=THIRD) );
  //  if(Ftilde1<THIRD) Ftilde=0.0;

  // FLASH Equation 47 (shock must have convergence)
  Ftilde2 = V[1] - V[-1];
  Ftilde *= ( (FTYPE)(Ftilde2<=0.0) );
  //  if(Ftilde2>0.0) Ftilde=0.0;

  //  dualfprintf(fail_file,"Ftilde=%21.15g\n",Ftilde);
  
  return( Ftilde );
}

FTYPE divftilde( int dir, int shift, FTYPE *Vabs, FTYPE *Pabs)
{
  FTYPE Ftilde,Ftilde1,Ftilde2;
  FTYPE Sv;
  FTYPE *V, *P;
  FTYPE V2diff,Vbottom;


  // shift as needed
  P = Pabs + shift;
  V = Vabs + shift;

  // (flow divergence)
  Sv = fabs(V[1] - V[-1]) / ( fabs(V[1]) + fabs(V[-1]) + SMALL );

  // (flow divergence must be sufficient)
  Ftilde = MAX(SV0,MIN(+1.0,Sv));

  Ftilde2 = V[1] - V[-1];
  Ftilde *= ( (FTYPE)(Ftilde2>=0.0) );
  
  return( Ftilde );
}


/// PPM FLATTENERS (final formula)
FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P)
{
  FTYPE ftilde( int dir, int shift, FTYPE *P, FTYPE *V);
  int signdP;
  FTYPE Fi;

  signdP = (P[1] - P[-1] > 0) * 2 - 1;
  // FLASH Equation 48
  Fi = max( ftilde(dir, 0, V,P), ftilde(dir, -signdP, V,P) );

  return(Fi);
}

/// Jon's divergence condition
/// currently Fi check is redundant
FTYPE  Divcalc(int dir, FTYPE Fi, FTYPE *V, FTYPE *P)
{
  FTYPE divftilde( int dir, int shift, FTYPE *Vabs, FTYPE *Pabs);
  FTYPE ftilde( int dir, int shift, FTYPE *P, FTYPE *V);
  int signdV;
  FTYPE Div;

  signdV = (V[1] - V[-1] > 0) * 2 - 1;
  Div = max( divftilde(dir, 0, V,P), divftilde(dir, -signdV, V,P) );

  // now put back actual scale
  Div *= (V[1] - V[-1]);

  return(Div);
}



/// Get pressure
/// Note this is quite inefficient since operating per-point get same pressure for entire line multiple times
/// (SUPERGODMARK: Also no accounting of magnetic field)
void getPressure(int whicheom, int i, int j, int k, int loc, FTYPE **yrealpl, FTYPE *P)
{
  int mm;

#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
  FTYPE tautot[NDIM],tautotmax;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  get_geometry(i,j,k,loc,ptrgeom);
  FTYPE prreal[NPR];
  int pl;
#endif

  // need pressure over full range from -3..3
  for(mm=-interporder[PARAFLAT]/2;mm<=interporder[PARAFLAT]/2;mm++){
    P[mm] = 0.0;

    if(whicheom==EOMSETMHD){
      P[mm] += pressure_rho0_u_simple(i, j, k, loc, yrealpl[RHO][mm],yrealpl[UU][mm]);
    }
#if(RADSHOCKFLAT&&EOMRADTYPE!=EOMRADNONE)
    if(whicheom==EOMSETRAD || whicheom==EOMSETMHD){
      // add radiation pressure to total pressure if optically thick
      PALLREALLOOP(pl) prreal[pl]=yrealpl[pl][mm]; // reorder
      //      calcfull_tautot(prreal, ptrgeom, tautot, &tautotmax);
      calc_tautot(prreal, ptrgeom, NULL, tautot, &tautotmax); // very accurate tautot not necessary, so use Tgas=Trad assumption in opacity

      P[mm] += MIN(tautotmax,1.0)*(4.0/3.0-1.0)*yrealpl[PRAD0][mm]; // KORALNOTE: recall pressure just along diagonal and no velocity in R^\mu_\nu
    }
#endif
  }

}




void parasteep(int dir, int pl, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r)
{
  int odir1,odir2;
  void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  FTYPE etaicalc(int pl, FTYPE *y, FTYPE *V, FTYPE *P);
  FTYPE etai;



#if(DOPPMSTEEPVARTYPE==0)
  if(pl==RHO)
#elif(DOPPMSTEEPVARTYPE==1)
    // define orthogonal directions for field steepening
    odir1=dir%3+1;
  odir2=(dir+1)%3+1;
  if(pl==RHO || pl==B1+odir1-1 || pl==B1+odir2-1)
#endif
    {
      // get contact indicator
      etai=etaicalc(pl,y,V,P);
      // get steepend values
      parasteepgen(pl, etai, V, P, y, dq, l, r);
    }


}



/// steepener
/// doesn't use l,r for any calculations, so if steepens fully then initial l,r values don't matter
/// return etai if needed
void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r)
{
  void pr_contact_compute(int pl, FTYPE *y, FTYPE *dq, FTYPE *prld, FTYPE *prrd);
  FTYPE prld, prrd;
  FTYPE l0,r0;
  FTYPE lmc,rmc;
  FTYPE mceta;


  // compute anti-disspiative left,right values
  pr_contact_compute(pl,y,dq,&prld,&prrd);
  
  // switch to MC for original l,r states if steepening      
  mceta=4.0*max( 0.25 - etai,0.0 );

  lmc = y[0] - 0.5*dq[0];
  rmc = y[0] + 0.5*dq[0];

  l0 = (*l) * mceta  + lmc*(1.0-mceta);
  r0 = (*r) * mceta  + rmc*(1.0-mceta);
  
  // assign steepened density value
  *l = l0 * ( 1.0 - etai ) + prld*etai;
  *r = r0 * ( 1.0 - etai ) + prrd*etai;
  // else make no changes to l,r   

  

}



void pr_contact_compute(int pl, FTYPE *y, FTYPE *dq, FTYPE *prld, FTYPE *prrd)
{

  // equation 33 and 34 in FLASH
  *prld=y[-1]+0.5*dq[-1];
  *prrd=y[+1]-0.5*dq[+1];


}

/// PPM steepener parameter, where etai=1 is steep and etai=0 is normal not steep
/// Acts to modify l and r values so actually nonmonotonic compared to surrounding cells -- This results in diffusion term in HLL or LAXF to actually be an anti-diffusion term causing (e.g. mass) to be sucked back into the cell in order to counter the diffusive flux causing spreading.  In a stationary case the expansive term balances the anti-diffusion term even with jumps at the interface.
/// 
FTYPE etaicalc(int pl, FTYPE *y, FTYPE *V, FTYPE *P)
{
  FTYPE delta2l,delta2r;
  FTYPE etatilde;
  FTYPE etai;
  FTYPE ddcoef;
  int ii;
  FTYPE Pjump,dB,Bmean;
  FTYPE prjumpfactor;
  FTYPE ifinf;
  FTYPE cs2;
  int mm;
  int mmstart,mmend;
  FTYPE max3P,min3P,min3y;

  // y is accessed from y[-2..2]
  // P,dq is accessed via dq[-1,0,1]


#if(AVGINPUT)
  ddcoef=SIXTH;
#else
  ddcoef=1.0/8.0;
#endif

  // equation 35 in FLASH
  ii=-1;
  delta2l=ddcoef*( (y[ii+1]-y[ii]) - (y[ii] - y[ii-1]) );

  // equation 35 in FLASH
  ii=1;
  delta2r=ddcoef*( (y[ii+1]-y[ii]) - (y[ii] - y[ii-1]) );

  // equation 36 in FLASH
  // sign of denominator is important
  // we multiply by conditionall that is 0 or 1
  if(fabs(y[1]-y[-1])>SMALL){
    etatilde=(-(delta2r-delta2l)/(y[1]-y[-1]));
  }
  else etatilde=0.0;

  // below can lead to asymmetries
  //ifinf = ((FTYPE)(fabs(y[1]-y[-1])<0.5*SMALL)); // inf avoidance
  //etatilde=(-(delta2r-delta2l)/(y[1]-y[-1] + ifinf*SMALL));

  

  /////////////////////////////////
  // Pressure jump
  max3P=SMALL;
  min3P=BIG;
  min3y=BIG;
  //  for(mm=-2;mm<=2;mm++){
  // upwinded extension of check
  //  if(V[0]>0.0){ mmstart=-2; mmend=1; }
  //  else if(V[0]<0.0) {  mmstart=-1; mmend=2; }
  //  else {  mmstart=-1; mmend=1; }
  mmstart=-1; mmend=1;
  for(mm=mmstart;mm<=mmend;mm++){
    if(mm==0) continue; // when wanting original method
    max3P=max(max3P,fabs(P[mm]));
    min3P=min(min3P,fabs(P[mm]));
    min3y=min(min3y,fabs(y[mm]));
  }
  min3P+=SMALL;

  //  min3y=min(min(fabs(y[-1]),fabs(y[1])),fabs(y[0]));

  //Pjump=fabs(P[1]-P[-1])/(min(fabs(P[1]),fabs(P[-1]))+SMALL);
  // need more strict and extensive pressure jump check
  Pjump=fabs(max3P - min3P)/min3P;


 

  // consistently use energy density like term
  if(pl>=B1 && pl<=B3){
    dB = y[1]-y[-1];
    Bmean = 0.5*(y[-1] + y[1]); // mean field
    //    prp1sq=0.5*y[1]*y[1]*sign(y[1]); // magnetic pressure
    //    prm1sq=0.5*y[-1]*y[-1]*sign(y[-1]);
    prjumpfactor=fabs(dB)/(fabs(Bmean)+fabs(dB)+SMALL);

    //    dualfprintf(fail_file,"dB=%21.15g Bmean=%21.15g prjumpfactor=%21.15g\n",dB,Bmean,prjumpfactor);
  }
  else{
    prjumpfactor=fabs(y[1]-y[-1])/min3y;

    // check effective c_s given Pressure leaks into steepened density region
    // ensure contributes negligibly to dynamics

    //    cs2 = max3P/min3y;

    // want cs2 dt^2/dx^2 << 1 to avoid steepening leading to velocity comparable to whatever Courant condition is limited by

    // relativistic check
    // don't change energy per baryon by much
    // if pressure is relatively flat
    //    cs2m1check = (cs2 < 5.0 * fabs(P[-1])/fabs(y[-1]));
    //    cs20check = (cs2 < 5.0 * fabs(P[0])/fabs(y[0]));
    //    cs2p1check = (cs2 < 5.0 * fabs(P[1])/fabs(y[1]));
  }




  // equation 39 in FLASH (not a shock)
  //  if( Pjump > 0.1*prjumpfactor) etatilde=0.0;
  etatilde *= (FTYPE)( Pjump <= 0.1*prjumpfactor);



  // equation 37 in FLASH (to avoid triggering on numerical noise)
  //  if(prjumpfactor<0.01) etatilde=0.0;
  etatilde *=  (FTYPE)(prjumpfactor>=0.01);
  //  etatilde *=  (FTYPE)(prjumpfactor>=0.1);

  // equation 38 in FLASH (really a contact)
  //if(delta2l*delta2r>0.0) etatilde=0.0;
  etatilde *=  (FTYPE)(delta2l*delta2r<=0.0);



  // equation 40 in FLASH
  etai=max(0.0,min(20.0*(etatilde-0.05),1.0));

  return(etai);

}















