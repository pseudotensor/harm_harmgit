#include "decs.h"

#include "reconstructeno.h"


#include "para_and_paraenohybrid.h"


int paraenohybrid_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array,  int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                            FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM],
                            FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                            FTYPE *Pindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp ) 
{
  // JCM STUFF START
  int paraprocess_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                            FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE *parafrac, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM],
                            FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                            FTYPE *Pindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, struct of_trueijkp *trueijkp  ) ;
  int returnedvalue;
  FTYPE a_parafrac[NBIGM];
  FTYPE *parafrac;
  int allweno,allpara;
  int i;
  FTYPE stiffnonlocal;
  int mm;
  int sum;











  // JCM STUFF FINISH
#if( DO_CENT2EDGE == 0 )
  for( i = ps; i <= pf; i++ ) {
    if(monoindicator[LEFTYOUT][i]==0) yout_left[i] = yin[i];
    if(monoindicator[RIGHTYOUT][i]==0) yout_right[i] = yin[i];
    //yout_left[i] = yin[i];
    //    yout_right[i] = yin[i];
  }
  return( 0 );
#else





#if(COUNTCALLS)
  sum=0;
  for( i = ps; i <= pf; i++ ) {
    sum+=(monoindicator[LEFTYOUT][i]==0)+(monoindicator[RIGHTYOUT][i]==0);
  }
  dualfprintf(fail_file,"t=%21.15g nstep=%ld pl=%d sum_c2e=%d\n",t,nstep,pl,sum/2);
#endif



  // offset pointer
  parafrac=(FTYPE (*)) (&(a_parafrac[NBIGBND]));





  if(FULLHYBRID==1 && PARAMODWENO==1){


    // GODMARK: could also modify bs,ps,pf,bf for each weno/para for range actually doing weno or para
    allweno=allpara=1;
    for( i = bs+1; i <= bf-1; i++ ) { // same range as in  compute_df_line_new() in interpline.c


      // in case parafrac==0 and WENO not called and monofrac==0, then need these to be 0 to avoid nan*0=nan    
      if(monoindicator[LEFTYOUT][i]==0.0) yout_left[i]=0.0;
      if(monoindicator[RIGHTYOUT][i]==0.0) yout_right[i]=0.0;

      // new method that uses para fully up to stiff=0.25 and mix WENO and para up to stiff=0.75 and then use WENO for 0.75 to 1.0
#define PARAEND 0.2
#define MIXEND 0.4

      // use extended maximum stiff indicator to avoid switching to PARA too quickly
      stiffnonlocal=0.0;
      for(mm=-1;mm<=1;mm++){ // allowed since stiffindicator defined over full bs..bf range
        stiffnonlocal=max(stiffnonlocal,stiffindicator[i+mm]);
      }

      if(stiffnonlocal<PARAEND){// assumes stiffindicator already well-truncated (and is)
        parafrac[i]=1.0;
        allweno=0;
      }
      else if(stiffnonlocal<MIXEND){
        parafrac[i] = 1.0 - (stiffnonlocal-PARAEND)/(MIXEND-PARAEND);
        allweno=allpara=0;
      }
      else{
        parafrac[i]=0.0;
        allpara=0;
      }
    }
  }
  else{
    allweno=allpara=0;
  }


  if(FULLHYBRID==1 && !allpara || FULLHYBRID==0 || PARAMODWENO==0){
    returnedvalue=eno_line_c2e( whichquantity, dir, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pf, bf, minorderit, maxorderit, shiftit, shockindicator,  stiffindicator, Vline, Pline, df, dP, etai, monoindicator, Pindicator, yin, yout_left, yout_right, youtpolycoef,trueijkp ); 
  }




  if(PARAMODWENO && MAXBND>=4){ // then assume can store enough dq's for PARA-like steepen/flat/monocheck

    if(FULLHYBRID==1 && !allweno || FULLHYBRID==0){
      paraprocess_line_c2e( whichquantity, dir, do_weight_or_recon, stencil_weights_array, whichreduce, preforder, pl, bs, ps, pf, bf, minorderit, maxorderit, shiftit,  shockindicator, stiffindicator, parafrac, Vline, Pline, df, dP, etai, monoindicator, Pindicator, yin, yout_left, yout_right,trueijkp );  
    }
  }


  return(returnedvalue);
#endif

}







// used to make WENO work better in contact
// works well as setup for isolated stationary and moving contact and blast wave, but problem for TESTNUMBER==4 contact compared to MCSTEEP,PARALINE,PARAFLAT still -- should look at etai's at late time in contact to see why -- maybe shock not flat enough
// GODMARK: If parafrac==0, then must behave *exactly* like weno would have!
int paraprocess_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, 
                          FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE *parafrac, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM],
                          FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], 
                          FTYPE *Pindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, struct of_trueijkp *trueijkp  ) 
{
  extern void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont);
  extern void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  extern void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r);
  extern void checkparamonotonicity(int smooth, int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout);
  void jonparasmooth_compute(int realisinterp, int dqrange, int pl, FTYPE *y, FTYPE *dq1, FTYPE *dq2, FTYPE *lout, FTYPE *rout, int *smooth);
  FTYPE a_face[2][NBIGM];
  FTYPE (*face)[NBIGM];
  int i;
  int dqrange;
  int odir1,odir2;
  FTYPE parafraclocal;
  FTYPE monofrac;
  FTYPE truelowerorderfraction;
  FTYPE myetai,myshock;
  int mm;
  FTYPE mymono;
  FTYPE leftmc,rightmc;
  FTYPE leftweno,rightweno;
  int smooth;


  //////////////////////
  //
  // for now only apply to density (only works if don't introduce switching effects between para and weno, and switching seems to be avoided if use para within -3..3 of contact since WENO is then removed from reaching into contact
  //
  /////////////////////
#if(FULLHYBRID==0)
  if(pl!=RHO) return(0);
#endif


  // stencil_weights_array[i].lower_order_fraction;
  // Still have this local file-global quantity, so let's perform extra PARA-like operations here

  // pointer shift
  face=(FTYPE (*)[NBIGM]) (&(a_face[0][NBIGBND]));

  // define orthogonal directions for field steepening
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;


  //////////////////////////////
  //
  // copy over WENO result since will only use PARA modifications if reduced to WENO3
  //
  /////////////////////////////


#if(0)// choose weno as base for lower order


  for( i = ps; i <= pf; i++ ) {
    face[0][i]=yout_left[i];
    face[1][i]=yout_right[i];
  }


#elif(1) // choose para as base for lower order


  for( i = ps; i <= pf+1; i++ ) {
    paracont(df[DF2OFMONO][i], &yin[i], &face[0][i]);
  }

  
  // default left and right states
  // 1 input and 2 outputs
  for( i = ps-1; i <= pf; i++ ) {
    // left from y[i]
    //      yout_left[i]=face[0][i];
    // right from y[i]
    //      yout_right[i]=face[1][i]=face[0][i+1];

    face[1][i]=face[0][i+1];
  }

#elif(0)// choose weno as base for lower order
  // MC also doesn't work?
  for( i = ps; i <= pf; i++ ) {
    face[0][i]= yin[i] - 0.5*df[DFMONO][i];
    face[1][i]= yin[i] + 0.5*df[DFMONO][i];
  }
#endif



  //dqrange=preforder-2;
  dqrange=5;

  int whicheom;
  if(RADFULLPL(pl)) whicheom=EOMSETRAD;
  else whicheom=EOMSETMHD;



  // default left and right states
  // 1 input and 2 outputs
  for( i = ps; i <= pf; i++ ) {


#if(JONPARASMOOTH)
    int realisinterp=1; // assume not big deal
    jonparasmooth_compute(realisinterp,dqrange,pl,&yin[i],&df[DFONESIDED][i],&df[DFCENT][i],&face[0][i],&face[1][i],&smooth);
#else
    smooth=0;
#endif


#if(DOPPMCONTACTSTEEPMODWENO)

#if(CONTACTINDICATOR==0)
#error If Steepen in paraline, must turn on CONTACTINDICATOR
#endif
    
    if(smooth==0 && (whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT)) parasteepgen(pl,etai[whicheom][i],&Vline[whicheom][i],&Pline[whicheom][i],&yin[i],&df[DFMONO][i],&face[0][i],&face[1][i]);
#endif
  
  
#if( DOPPMREDUCEMODWENO )
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT) paraflatten(dir,pl,&yin[i],shockindicator[whicheom][i],&face[0][i],&face[1][i]);
#endif
  

#if(1)
    // was checking only para result
    checkparamonotonicity(smooth, dqrange, pl, &yin[i], &df[DF2OFMONO][i], &df[DFMONO][i],&face[0][i],&face[1][i],&face[0][i],&face[1][i]);
#endif



    // now see if want to use MONO result and combine with para result
    // doesn't seem to help moving Gresho
#if(1)
    mymono=min(max(monoindicator[MONOINDTYPE][i],0.0),1.0);

    // modify PARA result using MONO if monoindicator>0
    face[0][i] = face[0][i] * (1.0-mymono) +  yout_left[i] * mymono;
    face[1][i] = face[1][i] * (1.0-mymono) +  yout_right[i] * mymono;
#endif


    // now merge WENO and PARA result based upon lower order fraction

    //      if(etai>=0.01) etai=1.0;
    //      if(pl==RHO){
    // dualfprintf(fail_file,"i=%d etai=%21.15g\n",i,etai);
    //      }
    //      parafraclocal = etai;
    //parafraclocal=1.0;
    // GODMARK: Unsure if lower order fraction accounts for MONO

#if(FULLHYBRID==0)
    myetai=0;
    myshock=0;
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT){
      // seem to require -3..3 for blast wave to avoid WENO smearing of contact
      for(mm=-3;mm<=3;mm++){
        myetai=max(myetai,etai[whicheom][i+mm]);
        myshock=max(myshock,shockindicator[whicheom][i+mm]);
      }
    }
    monofrac = min(max(monoindicator[MONOYIN][i],0.0),1.0);

    // monofrac=0.0;
    truelowerorderfraction = (stencil_weights_array[i].lower_order_fraction)*(1.0-monofrac);

    //    truelowerorderfraction = (1.0-monofrac);
    //parafraclocal = (1.0-monofrac)*(1.0-myshock)*myetai;

    //(1.0-monofrac)*(1.0-myshock)*

    //    myetai=1.0;

    //    parafraclocal = min(max(1.5*myetai,0.0),1.0);
    // use para if ANY indication that contact is near
    // 2.0*myetai seems to work fine once set range to be -3..3
    parafraclocal = min(max(2.0*myetai,0.0),1.0);
    //    parafraclocal = min(max(max(2.0*myetai,myshock),0.0),1.0);

    // don't use para if in stiff regime since WENO more robust
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT){
      parafraclocal *= (1.0-stiffindicator[whicheom][i]);
    }

    //    if(pl==RHO) parafraclocal=1.0;

    // don't use para if doing higher order WENO
    //    parafraclocal = min(parafraclocal,truelowerorderfraction);

    // below doesn't work (generates spike at contact) probably because need some monoindicator across all quantities so para or weno used consistently for all quantities in contact
    //    parafraclocal = min(parafraclocal,(1.0-monofrac));

    //dualfprintf(fail_file,"i=%d etai=%21.15g myetai=%21.15g\n",i,etai[i],myetai);


    //    parafraclocal = min(max(truelowerorderfraction,2.0*etai[i]),1.0);
    //    parafraclocal = min(1.0*max(max(truelowerorderfraction,etai),0.0),1.0);
    //    parafraclocal = min(1.0*max(10.0*etai,0.0),1.0);
    //parafraclocal = min(max(1.0*max(truelowerorderfraction,1.1*etai),0.0),1.0);
    //parafraclocal = min(max(1.0*max(truelowerorderfraction,1.1*etai),0.0),1.0);
    //    parafraclocal = 1.0;
    //    parafraclocal = min(max(1.0*max(truelowerorderfraction,etai),0.0),1.0);
    //    if(parafraclocal<0.9) parafraclocal=0.0;
    //    parafraclocal = min(10.0*max(parafraclocal-0.9,0.0),1.0);
    //    parafraclocal = min(max(parafraclocal*(parafraclocal-0.1),0.0),1.0);
    //    if(parafraclocal<0.3) parafraclocal=0;
    //    parafraclocal=0.98;

    //    parafraclocal=0.0;

    //parafraclocal = min(max(etai,0.0),1.0);
    //parafraclocal = 1.0;
    //parafraclocal = 1.0;

    // GODMARK: DEBUG:
    //stencil_weights_array[i].lower_order_fraction=1.0;


#elif(FULLHYBRID==1)

    // if not stiff, then speculate with para.  If stiff, then stay robust with WENO
    //    parafraclocal = (1.0-stiffindicator[i]);

    // new method that uses para or weno fully and only mixes in middle of stiffindicator
    parafraclocal = parafrac[i];

#endif



#if(USEMCFORLOWERORDERWENO)
    leftmc  = yin[i] - 0.5* df[DFMONO][i];
    rightmc = yin[i] + 0.5* df[DFMONO][i];

    // flatten MC in shocks
#if( DOPPMREDUCEMODWENO )
    if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT) paraflatten(dir,pl,&yin[i],shockindicator[whicheom][i],&leftmc,&rightmc);
#endif

    monofrac = min(max(monoindicator[MONOYIN][i],0.0),1.0);
    truelowerorderfraction = min(max((stencil_weights_array[i].lower_order_fraction)*(1.0-monofrac),0.0),1.0);

    leftweno  = yout_left[i]  * (1.0-truelowerorderfraction) + leftmc *truelowerorderfraction;
    rightweno = yout_right[i] * (1.0-truelowerorderfraction) + rightmc*truelowerorderfraction;

#else
    // normal weno
    monofrac = min(max(monoindicator[MONOYIN][i],0.0),1.0);

    leftweno  = yout_left[i];
    rightweno = yout_right[i];
#endif



      
    yout_left[i]  = leftweno  * (1.0-parafraclocal) + face[0][i] * parafraclocal ;
    yout_right[i] = rightweno * (1.0-parafraclocal) + face[1][i] * parafraclocal ;

    //    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d :: etai=%21.15g myetai=%21.15g parafraclocal=%21.15g :: myshock=%21.15g monofrac=%21.15g truelow=%21.15g \n",nstep,steppart,i,etai[i],myetai,parafraclocal,myshock,monofrac,truelowerorderfraction);


#if(0)
    // if using any weno part, then can't modify state so behaves like weno in limit that parafrac==0.  Otherwise, e.g., MPI boundaries not consistently computing fluxes

#if(FULLHYBRID==0 || FULLHYBRID==1)
    //    if(1||parafraclocal>0.01){
    
    if(monofrac<0.1){
      //    if(parafraclocal>0.01){
      // check monotonicity of FINAL result rather than just para result, to ensure monotonic
      dqrange=5;
      checkparamonotonicity(dqrange, pl, &yin[i], &df[DF2OFMONO][i], &df[DFMONO][i],&yout_left[i],&yout_right[i],&yout_left[i],&yout_right[i]);
    }
#endif


#if( DOPPMREDUCEMODWENO )
    // flatten again in case checkparamonotonic is not reducing all the way to DONOR
    // causes major problems
    //     if(whicheom==EOMSETMHD || whicheom==EOMSETRAD&&RADSHOCKFLAT) paraflatten(dir,pl,&yin[i],shockindicator[whicheom][i],&yout_left[i],&yout_right[i]);
#endif


#endif

  }// end loop over i


  return(0);


}


