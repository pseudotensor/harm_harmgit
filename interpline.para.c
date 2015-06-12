
/*! \file interpline.para.c
     \brief Parabolic/PPM Spatial Interpolation for fluxes based upon providing full 1D line information
     // Instead of acting per point, this acts per line to improve memory efficiency.
*/
#include "decs.h"


#include "para_and_paraenohybrid.h"


/// This version has pressure with total pressure, which is more correct than point version
void pass_1d_line_multipl_paraline(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp)
{
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pllocal;
  int numprims;
  int plstart;
  int pl,pliter;
  int mypl;
  int myi;
  int is_copy;
  void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont);
  void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r);
  void checkparamonotonicity(int smooth, int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout);
  void jonparasmooth_compute(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *dq1, FTYPE *dq2, FTYPE *lout, FTYPE *rout, int *smooth);
  FTYPE a_facecont[MAXNPR][NBIGM];
  FTYPE (*facecont)[NBIGM];
  int i;
  int dqrange;
  int odir1,odir2;
  FTYPE left,right;
  FTYPE mymono;
  int smooth;
  int whicheom;


#if(DOPPMREDUCE && SHOCKINDICATOR==0)
#error "Paraline needs SHOCKINDICATOR==1 when DOPPMREDUCE==1"
#endif

#if(DOPPMCONTACTSTEEP && CONTACTINDICATOR==0)
#error "If Steepen in paraline, must turn on CONTACTINDICATOR"
#endif




  // pointer shift
  facecont=(FTYPE (*)[NBIGM]) (&(a_facecont[0][NBIGBND]));


  // define orthogonal directions for field steepening
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;


  /////////////////
  //
  // Define which quantities (pl) to operate on
  //
  /////////////////

  setup_nprlocalist(whichquantity,&nprlocalstart,&nprlocalend,nprlocallist,&numprims);


  /////////////////
  //
  // loop over line first getting continuous solution at face
  //
  ////////////////

  // for( i = ps; i <= pe; i++ ) {
  NUMPRIMLOOP(pliter,pl) for( i = ps; i <= pe+1; i++ ) {
    paracont(df[pl][DF2OFMONO][i], &yin[pl][0][i], &facecont[pl][i]);
  }

  
  //dqrange=preforder-2;
  dqrange=5;

  // default left and right states
  // 1 input and 2 outputs
  NUMPRIMLOOP(pliter,pl) for( i = ps; i <= pe; i++ ) {
    // left from y[i]
    //    yout[pl][0][i]=facecont[pl][i];
    left=facecont[pl][i];
    // right from y[i]
    //    yout[pl][1][i]=facecont[pl][i+1];
    right=facecont[pl][i+1];

    if(RADFULLPL(pl)) whicheom=EOMSETRAD;
    else whicheom=EOMSETMHD;

#if(JONPARASMOOTH)
    int realisinterp=1; // assume not big deal
    jonparasmooth_compute(realisinterp,dqrange,pl,&yin[pl][0][i],&df[pl][DFONESIDED][i],&df[pl][DFCENT][i],&left,&right,&smooth);
#else
    smooth=0;
#endif


#if(DOPPMCONTACTSTEEP)
    if(smooth==0 && (whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT)) parasteepgen(pl,etai[pl][whicheom][i],&Vline[whicheom][i],&Pline[whicheom][i],&yin[pl][0][i],&df[pl][DFMONO][i],&left,&right);
#endif
  
  
#if( DOPPMREDUCE )
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT) paraflatten(dir,pl,&yin[pl][0][i],shockindicator[whicheom][i],&left,&right);
#endif
  

    dqrange=interporder[PARALINE]-2;
    checkparamonotonicity(smooth, dqrange, pl, &yin[pl][0][i], &df[pl][DF2OFMONO][i], &df[pl][DFMONO][i], &left,&right,&left,&right);


    // now see if want to use MONO result and combine with para result
    // doesn't seem to help moving Gresho
#if(PARALINEUSESMONO)

    // in case parafrac==0 and WENO not called and monofrac==0, then need these to be 0 to avoid nan*0=nan    
    if(monoindicator[pl][MONOLEFTSET][i]==0.0) yout[pl][0][i]=0.0;
    if(monoindicator[pl][MONORIGHTSET][i]==0.0) yout[pl][1][i]=0.0;

    mymono=min(max(monoindicator[pl][MONOINDTYPE][i],0.0),1.0);

    //    mymono=1.0;
    
    //if(mymono>0.1) dualfprintf(fail_file,"i=%d mymono=%21.15g\n",i,mymono);

    yout[pl][0][i] = left * (1.0-mymono) + yout[pl][0][i] * mymono;
    yout[pl][1][i] = right * (1.0-mymono) + yout[pl][1][i] * mymono;
#else
    // no mono used, just assign para result
    yout[pl][0][i] = left;
    yout[pl][1][i] = right;
#endif

    // can force MINM for radiation, but don't need to.
    if(0&&whicheom==EOMSETRAD){
      FTYPE extremum=0.0;
      FTYPE theta = 1.0;
      FTYPE Dqm = theta * df[pl][DFONESIDED][i];
      FTYPE Dqp = theta * df[pl][DFONESIDED][i+1];
      FTYPE Dqc = df[pl][DFCENT][i];
      FTYPE mydq = MINMODGEN(extremum,MINMODGEN(extremum,Dqm,Dqc),Dqp);
      
      yout[pl][0][i] = yin[pl][0][i] - 0.5*mydq;
      yout[pl][1][i] = yin[pl][0][i] + 0.5*mydq;
    }

#if(NONMONOLIM>0 && DOPPMREDUCE)
    // then flatten final result
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT) paraflatten(dir,pl,&yin[pl][0][i],shockindicator[whicheom][i],&yout[pl][0][i],&yout[pl][1][i]);
#endif


  }
  
  
}




/// Pass 1D line to PARALINE scheme
void pass_1d_line_paraline(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp)
{


  // assume never need this function
  dualfprintf(fail_file,"pass_1d_line_paraline() not setup\n");
  myexit(2896262);


}



/// get para slopes
void get_limit_slopes_paraline(FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dq)
{

  void get_limit_slopes(int reallim, int extremum, FTYPE *dq1l, FTYPE *dq1r, FTYPE *dq2, FTYPE *dqout);

  get_limit_slopes(PARALINE2LIM, PARALINE2EXTREME, dq1l, dq1r, dq2, dq);

}
