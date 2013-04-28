
#include "decs.h"
#include "reconstructeno.h"

#define DEBUGMONO 0

#define DEBUGMONONEW 0

// fractional difference has to be monotonic to within this error
#define ERRORNORM (1E3*NUMEPSILON)
//#define ERRORNORM (1E-6)
//#define ERRORNORM (0.0)

// whether to see if inputs are monotonic (perhaps not actually necessary if final interpolated YOUTs are monotonic)
// 0: don't check yin for monotonicity
// 1: use df
// 2: use ddf
// 3: use df&&ddf
#define CHECKYIN 3

#if(CHECKYIN!=0)

#define CHECKYINIF ((!CHECKYIN) ||                                      \
                    ((CHECKYIN==1)&& (dpup||dpdown))||                  \
                    ((CHECKYIN==2)&& (ddpup||ddpdown))||                \
                    ((CHECKYIN==3)&& (dpup||dpdown) && (ddpup||ddpdown)))
#else

#define CHECKYINIF (1)

#endif


// whether to check outputs for monotonicity (crucial)
// 0: don't check
// 1: check function monotonicity
// 2: check only derivatives monotonicity
// 3: check both f and df
// 4: check f, df, and d^f
#define CHECKYOUTC2E 4

// for c2e
// 0: WENO5
// 1: ENO4 // using WENO5 type checks, leads to slightly sharper contacts, but also to a droop at the top due to switching to SWENO5
// 2: new ENO4 doesn't work
#define WHICHENOTYPE 0

// whether to use shock indicator in getting result
#define USESHOCKINDICATOR 0

// whether to use cusp indicator
// 0 : don't
// 1 : old sharp method
// 2 : new detector
#define USECUSPINDICATOR 1

// whether to check outputs for monotonicity (crucial)
// 0: don't check
// 1: check function monotonicity
// 2: check only derivatives monotonicity
// 3: check both f and df
// 4: check f, df, and fractional change
#define CHECKYOUTA2CC2A 4
// fractional change allowed by a2c or c2a
#define FRACCHANGEALLOWED (0.8)



// in the case that WENO5 fails to be monotonic, whether to reduce to WENO3 and see if WENO3 monotonic (if WENO3 isn't monotonic, reduces to ambiguous)
#define MONOFAIL2WENO3 0

// whether to reduce to a LIMITER when YIN's is not monotonic (not activated if CHECKYIN=0)
#define NONMONOTONIC2LIMITER 0

// whether to reduce to a LIMITER when found no good result (normally let Sasha-WENO handle this case)
// Using Sasha-WENO doesn't seem to help at all for test4
#define AMBIGUOUS2LIMITER 0



// monoindicator[0]: -1,0,1 for rough, ambiguous, and monotonic
// monoindicator[1]: whether set cell's left value
// monoindicator[2]: whether set cell's right value
// MERGEDC2EA2CMETHOD  TODO : NOT NECESSARY FOR NOW
void compute_jmonotonicity_line(int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM])
{
  int i;
  int dpup,dpdown,ddpup,ddpdown;
  int monoup,monodown;
  int check_for_cusp(FTYPE *yin, FTYPE *df, FTYPE *ddf);
  int check_for_cusp_indicator(FTYPE *yin, FTYPE *df, FTYPE *ddf, FTYPE *cuspindicator);
  void c2e_simple_limiter(int WHICHLIMITERTOREDUCETO, FTYPE *yin, FTYPE *valueleft, FTYPE *valueright);
  int order;
  void check_mono_input(FTYPE *yin, FTYPE *df, FTYPE *ddf, int *dpdown, int *dpup, int *ddpdown, int *ddpup);
  void check_mono_input_eno4(FTYPE *yin, FTYPE *df, FTYPE *ddf, int *dpdown, int *dpup, int *ddpdown, int *ddpup);
  void compute_check_output(int i, int bs, int be, int recontype, int preforder, FTYPE (*df)[NBIGM], FTYPE *yin, FTYPE (*monoindicator)[NBIGM], FTYPE (*yout)[NBIGM]);
  int pointdi;
  FTYPE pleft,pright,pcenter;
  FTYPE roughnessindicator,cuspindicator;


#if(MERGEDC2EA2CMETHOD)
  dualfprintf(fail_file,"JMONO not setup for MERGEDC2EA2CMETHOD\n");
  myexit(1987526);  
#endif

  if(EOMRADTYPE!=EOMRADNONE){
    dualfprintf(fail_file,"JMONO not setup for koral\n");
    myexit(394343677);      
  }

  if(preforder>5 || preforder<3){
    dualfprintf(fail_file,"JMONO not setup for preforder>5 || preforder<3\n");
    myexit(7176926);
  }

  ///////////////////
  //
  //   set default indicator to be rough and unset
  //
  ////////////////////

  if(recontype!=CVT_C2E){
    for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
      monoindicator[MONOINDTYPE][i]=0;
      monoindicator[MONOLEFTSET][i]=0;
    }
    //if( nstep > 0 || steppart >= 3 ) 
    return;
  }
  else{
    for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
      //      monoindicator[MONOINDTYPE][i]=-1;
      // Sasha-WENO can't handle roughness (-1) yet
      monoindicator[MONOINDTYPE][i]=0;
      monoindicator[MONOLEFTSET][i]=0;
      monoindicator[MONORIGHTSET][i]=0;
    }
    //if( nstep > 0 || steppart >= 3 ) 
    //return;
    // noticed that if always do S-WENO (return above), then normal Sasha-WENO is quite oscillatory for testnumber=9 (shock-entropy interaction)
  }

  // get number of points to jump over each iteration
  //  if(WHICHENOTYPE==0) pointdi=1;
  //  else pointdi=2;
  pointdi=1;


  ///////////////////
  //
  //   get monoindicators and yout if rough (-1) or monotonic (1), but don't set if ambiguous (0)
  //
  ////////////////////
  for(i=ps;i<=pe;i+=pointdi){ // loop over points of interest

    //////////////////////////////
    // rough checks

    // get cuspness indicator
#if(USECUSPINDICATOR==2)
    check_for_cusp_indicator(&yin[i], &df[DFONESIDED][i], &df[DF2OFONESIDED][i], &cuspindicator);
    //    dualfprintf(fail_file,"i=%d cuspindicator=%21.15g\n",i,cuspindicator);
#elif(USECUSPINDICATOR==1)


    
    // check for cusp
    if((check_for_cusp(&yin[i],&df[DFONESIDED][i],&df[DF2OFONESIDED][i]))  ){
      //dualfprintf(fail_file,"shocki[%d]=%21.15g\n",i,shockindicator[EOMSETMHD][i]);
      //if((shockindicator[EOMSETMHD][i]>=0.1)||(check_for_cusp(&df[DFONESIDED][i],&df[DF2OFONESIDED][i],&yin[i]))  ){// does somewhat better than not checking
      //if(0){


      if(recontype==CVT_C2E) c2e_simple_limiter(MINM, &yin[i], &yout[0][i],&yout[1][i]);
      // otherwise do nothing (no a2c/c2a conversion)
      else yout[0][i]=yin[i];

      // set as rough
      //monoindicator[MONOINDTYPE][i]=-1;
      // say unknown roughness
      monoindicator[MONOINDTYPE][i]=0;
      // say already set
      monoindicator[MONOLEFTSET][i]=1;
      monoindicator[MONORIGHTSET][i]=1;      
      //      dualfprintf(fail_file,"ROUGH i=%d\n",i);
    }
    else{ // then not rough
#endif

      //  monoindicator[MONOINDTYPE][i]=0;
      // say already set
      //      monoindicator[MONOLEFTSET][i]=0;
      //      monoindicator[MONORIGHTSET][i]=0;      
      //      continue;

      // check monotonicity of input
#if(CHECKYIN)
      if((WHICHENOTYPE==0)||(WHICHENOTYPE==1)) check_mono_input(&yin[i], &df[DFONESIDED][i], &df[DF2OFONESIDED][i], &dpdown, &dpup, &ddpdown, &ddpup);
      else check_mono_input_eno4(&yin[i], &df[DFONESIDED][i], &df[DF2OFONESIDED][i], &dpdown, &dpup, &ddpdown, &ddpup);
      //      dualfprintf(fail_file,"i=%d :: %d %d %d %d\n",i,dpdown,dpup,ddpdown,ddpup);
#endif

    

      //////////////////////////////
      // finally check that interpolated values are between original values
      if(CHECKYINIF){
        compute_check_output(i, bs, be, recontype, preforder, df, yin, monoindicator, yout);
        // if((pl==UU)&&(recontype==CVT_C2E)) dualfprintf(fail_file,"monoyout=%d\n",monoindicator[MONOINDTYPE][i]);
      
#if(MONOFAIL2WENO3)
        // then try again with lower order
        if(monoindicator[MONOINDTYPE][i]!=1){
          compute_check_output(i, bs, be, recontype, preforder-2, df, yin, monoindicator, yout);
        }
#endif

      }// end if monotonic based upon df's
      else{
        // dualfprintf(fail_file,"C2EINTPUT_NOTMONO: i=%d\n",i);
        monoindicator[MONOINDTYPE][i]=0;
#if(NONMONOTONIC2LIMITER)
        if(recontype==CVT_C2E) c2e_simple_limiter(DONOR, &yin[i], &yout[0][i],&yout[1][i]);
        monoindicator[MONOINDTYPE][i]=0;
        // say already set
        monoindicator[MONOLEFTSET][i]=1;
        monoindicator[MONORIGHTSET][i]=1;      
#endif
      }

#if(USECUSPINDICATOR==1)
    }// end else if not rough
#endif


#if(USESHOCKINDICATOR&&(USECUSPINDICATOR==2))
    roughnessindicator=max(shockindicator[EOMSETMHD][i],cuspindicator);
    //roughnessindicator=shockindicator[EOMSETMHD][i];
#elif(USECUSPINDICATOR==2)
    roughnessindicator=cuspindicator;
#elif(USESHOCKINDICATOR)
    roughnessindicator=shockindicator[EOMSETMHD][i];
#else
    roughnessindicator=0.0;
#endif

    // linearly combine result from DONOR (for roughness) with high order solution
    if((USESHOCKINDICATOR||(USECUSPINDICATOR==2))&&(monoindicator[MONOINDTYPE][i]==1)){
      if(recontype==CVT_C2E){
        c2e_simple_limiter(DONOR, &yin[i], &pleft, &pright);
   
        yout[0][i] = yout[0][i]*(1.0-roughnessindicator)+roughnessindicator*pleft;
        yout[1][i] = yout[1][i]*(1.0-roughnessindicator)+roughnessindicator*pright;
      }
      else{
        pcenter=yin[i]; // DONOR-like
        yout[0][i] = yout[0][i]*(1.0-roughnessindicator)+roughnessindicator*pcenter;
      }
    }
  




#if(AMBIGUOUS2LIMITER)
    if(recontype==CVT_C2E){
      if(monoindicator[MONOINDTYPE][i]==0){
        // fool around with what SWENO might do
        c2e_simple_limiter(DONOR, &yin[i], &yout[0][i],&yout[1][i]);
        monoindicator[MONOINDTYPE][i]=0;
        // say already set
        monoindicator[MONOLEFTSET][i]=1;
        monoindicator[MONORIGHTSET][i]=1;      
      }
    }
    else{
      yout[0][i]=yin[i];
      monoindicator[MONOINDTYPE][i]=0;// unknown roughness
      monoindicator[MONOLEFTSET][i]=1;// set
    }
#endif




  }// end loop over points





  /////////////////////////////////////////
  //
  // check symmetry and stop when broken
#if(0)
  // for test2
  for(i=ps;i<=pe;i++){ // loop over points of interest
    if(monoindicator[MONOINDTYPE][i]!=monoindicator[0][(pe+ps) - i]){
      dualfprintf(fail_file,"problem at i=%d io=%d : %d %d :: recontype=%d\n",i,pe+ps-i, monoindicator[MONOINDTYPE][i],monoindicator[0][pe-i],recontype);
      for(i=ps;i<=pe;i++){
        dualfprintf(fail_file,"mon[%d]=%d :: %21.15g %21.15g :: %21.15g :: %21.15g %21.15g\n",i,monoindicator[MONOINDTYPE][i],yin[i],yout[0][i],df[DF2OFONESIDED][i],df[DFONESIDED][i],df[0][i+1]);
      }
      myexit(1);
    }
  }
#endif




#if(0)
  // GODMARK
  // when S-WENO can handle roughness (-1), need to put -1 in cells that have substencils that cross into nowhere land (i.e. beyond bs or be)
  for(i=bs;i<=bs+preforder/2-1;i++){
    monoindicator[MONOINDTYPE][i]=-1;
  }

  for(i=be;i>=be-preforder/2-1;i--){
    monoindicator[MONOINDTYPE][i]=-1;
  }
#endif



  return;
}











// check monotonicity of 5 point stencil centered on cell center
void check_mono_input(FTYPE *yin, FTYPE *df, FTYPE *ddf, int *dpdown, int *dpup, int *ddpdown, int *ddpup)
{
  FTYPE norm;


  //////////////////////////////
  // mono checks
#if(0)
  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically decreasing
  *dpdown=(df[-1]>=df[0])&&(df[0]>=df[1])&&(df[1]>=df[2]);

  // monotonically increasing
  *dpup=(df[-1]<=df[0])&&(df[0]<=df[1])&&(df[1]<=df[2]);

  // monotonically curved downwards
  *ddpdown=(ddf[-1]>=ddf[0])&&(ddf[0]>=ddf[1]);

  // monotonically curved upwards
  *ddpup=(ddf[-1]<=ddf[0])&&(ddf[0]<=ddf[1]);

#else


#if(DEBUGMONONEW)
  if(pl==RHO){
    dualfprintf(fail_file,"yin[-2]=%21.15g yin[-1]=%21.15g yin[0]=%21.15g yin[1]=%21.15g yin[2]=%21.15g\n",yin[-2],yin[-1],yin[0],yin[1],yin[2]);

    dualfprintf(fail_file,"df[-1]=%21.15g df[0]=%21.15g ddf=%21.15g df[2]=%21.15g\n",df[-1],df[0],df[1],df[2]);
  }
#endif
  // df
  norm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL)*ERRORNORM;  //atch modify: added "-" in the front
      
  *dpdown=          (df[-1]>-norm);
  *dpdown= *dpdown&&(df[0]>-norm);
  *dpdown= *dpdown&&(df[1]>-norm);
  *dpdown= *dpdown&&(df[2]>-norm);
      
  *dpup=        (df[-1]<norm);
  *dpup= *dpup&&(df[0]<norm);
  *dpup= *dpup&&(df[1]<norm);
  *dpup= *dpup&&(df[2]<norm);

  // d^2f
  // using this norm can trigger on machine error
  //  norm=(fabs(df[-1])+fabs(df[0])+fabs(df[1])+fabs(df[2])+SMALL)*ERRORNORM;

  *ddpdown=           (ddf[-1]>-norm);
  *ddpdown= *ddpdown&&(ddf[0]>-norm);
  *ddpdown= *ddpdown&&(ddf[1]>-norm);
      
  *ddpup=         (ddf[-1]<norm);
  *ddpup= *ddpup&&(ddf[ 0]<norm);
  *ddpup= *ddpup&&(ddf[ 1]<norm);
  
#endif

}

// check monotonicity of 5 point stencil centered on cell center
void check_mono_input_eno4(FTYPE *yin, FTYPE *df, FTYPE *ddf, int *dpdown, int *dpup, int *ddpdown, int *ddpup)
{
  FTYPE norm;


  //////////////////////////////
  // mono checks
#if(0)
  // ONLY VALID FOR WENO5 (need to generalize)
  // monotonically decreasing
  *dpdown=(df[-1]>=df[0])&&(df[0]>=df[1]);

  // monotonically increasing
  *dpup=(df[-1]<=df[0])&&(df[0]<=df[1]);

  // monotonically curved downwards
  *ddpdown=(ddf[-1]>=ddf[0]);

  // monotonically curved upwards
  *ddpup=(ddf[-1]<=ddf[0]);

#else

  // df
  norm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;   //atch modify: added "-" in the front
      
  *dpdown=          (df[-1]>-norm);
  *dpdown= *dpdown&&(df[0]>-norm);
  *dpdown= *dpdown&&(df[1]>-norm);
      
  *dpup=        (df[-1]<norm);
  *dpup= *dpup&&(df[0]<norm);
  *dpup= *dpup&&(df[1]<norm);

  // d^2f
  *ddpdown=           (ddf[-1]>-norm);
  *ddpdown= *ddpdown&&(ddf[0]>-norm);
      
  *ddpup=         (ddf[-1]<norm);
  *ddpup= *ddpup&&(ddf[ 0]<norm);
  
#endif

}


//#define WHICHLIMITERTOREDUCETO PARA
//#define WHICHLIMITERTOREDUCETO MINM
//#define WHICHLIMITERTOREDUCETO DONOR

void c2e_simple_limiter(int WHICHLIMITERTOREDUCETO, FTYPE *yin, FTYPE *valueleft, FTYPE *valueright)
{
  FTYPE mydq;
  extern void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void para4(int realisinterp, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout);

  //  *monotype=0; // unknown roughness
  //  *monosetleft=1;
  //  *monosetright=1;

  if(WHICHLIMITERTOREDUCETO<PARA){
    // reduce to MINM
    slope_lim_3points(WHICHLIMITERTOREDUCETO, yin[-1], yin[0], yin[1], &mydq);
    //*left =y[0] - 0.5* mydq;
    //      *right=y[0] + 0.5* mydq;
    *valueleft=yin[0]-0.5*mydq;
    *valueright=yin[0]+0.5*mydq;
  }
  else if(WHICHLIMITERTOREDUCETO==PARA){
    para4(1,RHO,yin,valueleft,valueright);
  }


}




// compute c2e, a2c, or c2a, and check monotonicity of result
void compute_check_output(int i, int bs, int be, int recontype, int preforder, FTYPE (*df)[NBIGM], FTYPE *yin, FTYPE (*monoindicator)[NBIGM], FTYPE (*yout)[NBIGM])
{
  FTYPE pleft, pright, pcenter;
  extern void c2e_simple_weno(int order, int ii, int bs, int be, FTYPE *yin, FTYPE *pleft, FTYPE *pright);
  extern void c2e_simple_eno(int full_order, int is_interpolate_to_left, FTYPE *yin, FTYPE *pout);
  void a2c_simple_eno(int full_order, FTYPE *yin, FTYPE *pout);
  void c2a_simple_eno(int full_order, FTYPE *yin, FTYPE *pout);
  void check_mono_c2e_output(int i, int order, FTYPE *df, FTYPE *ddf, FTYPE pleft, FTYPE pright, FTYPE *yin, FTYPE *monotype, FTYPE *monosetleft, FTYPE *monosetright, FTYPE *valueleft, FTYPE *valueright);
  void check_mono_c2e_output_eno4(int i, int preforder, FTYPE *df, FTYPE *ddf, FTYPE pleft, FTYPE *yin, FTYPE (*monoindicator)[NBIGM], FTYPE (*yout)[NBIGM]);
  void check_mono_a2c_c2a_output(int preforder, FTYPE *df, FTYPE *ddf, FTYPE *pcenter, FTYPE *yin, FTYPE *monotype, FTYPE *monoset, FTYPE *value);
  


  if(recontype==CVT_C2E){
    // get edge values
    if(WHICHENOTYPE==0){
      // old version below
      //    c2e_simple_weno(preforder/2+1,i,bs,be,&yin[i],&pleft,&pright); // no small symmetry problem
      // new version below
      c2e_simple_eno(preforder,0,&yin[i],&pleft); // small asymmetry problems
      c2e_simple_eno(preforder,1,&yin[i],&pright); // small asymmetry problems

      // check if edge values are monotonic
      check_mono_c2e_output(i, preforder, &df[DFONESIDED][i], &df[DF2OFONESIDED][i], pleft, pright, &yin[i], &monoindicator[MONOINDTYPE][i], &monoindicator[MONOLEFTSET][i], &monoindicator[MONORIGHTSET][i],&yout[0][i],&yout[1][i]);
      
    }
    else if(WHICHENOTYPE==1){

      // pleft CENTER pright, and now edges have same value
      c2e_simple_eno(4,0,&yin[i],&pleft); // left face for i
      c2e_simple_eno(4,0,&yin[i+1],&pright); // right face for i

      // check if edge values are monotonic
      check_mono_c2e_output(i, preforder, &df[DFONESIDED][i], &df[DF2OFONESIDED][i], pleft, pright, &yin[i], &monoindicator[MONOINDTYPE][i], &monoindicator[MONOLEFTSET][i], &monoindicator[MONORIGHTSET][i],&yout[0][i],&yout[1][i]);

      if(monoindicator[MONOINDTYPE][i]==1){
        // then set other interface values to be equal
        yout[1][i-1]=pleft; // right face for i-1
        monoindicator[2][i-1]=1;
        monoindicator[0][i-1]=0;

        yout[0][i+1]=pright; // left face for i+1
        monoindicator[1][i+1]=1;
        monoindicator[0][i+1]=0;
      }


    }
    else if(WHICHENOTYPE==2){
      // pleft CENTER pright, and now edges have same value
      c2e_simple_eno(4,0,&yin[i],&pleft); // left face for i
      //      c2e_simple_eno(4,0,&yin[i+1],&pright); // right face for i

      // check if edge values are monotonic
      check_mono_c2e_output_eno4(i, preforder, &df[DFONESIDED][i], &df[DF2OFONESIDED][i], pleft, &yin[i], monoindicator,yout);
    }

  }
  else{
    if(recontype==CVT_A2C){
      a2c_simple_eno(preforder,&yin[i],&pcenter);
    }
    else if(recontype==CVT_C2A){
      c2a_simple_eno(preforder,&yin[i],&pcenter);
    }
    // check if new value is quasi-monotonic (quasi because comparing different types of quantities, averages and points)
    check_mono_a2c_c2a_output(preforder, &df[DFONESIDED][i], &df[DF2OFONESIDED][i], &pcenter, &yin[i], &monoindicator[MONOINDTYPE][i], &monoindicator[MONOLEFTSET][i],&yout[0][i]);
  }
}




#if(CHECKYOUTC2E==0)

// don't question, just accept
#define CHECKYOUTC2EIF (1)

#elif(CHECKYOUTC2E==1)

#define CHECKYOUTC2EIF (mup||mdown)

#elif(CHECKYOUTC2E==2)

#define CHECKYOUTC2EIF (dmup||dmdown)

#elif(CHECKYOUTC2E==3)

#define CHECKYOUTC2EIF ((mup||mdown)&&(dmup||dmdown))

#elif(CHECKYOUTC2E==4)

#define CHECKYOUTC2EIF ((mup||mdown)&&(dmup||dmdown)&&(ddmup||ddmdown))

#endif


// check monotonicity of final interpolated values
void check_mono_c2e_output(int i, int preforder, FTYPE *df, FTYPE *ddf, FTYPE pleft, FTYPE pright, FTYPE *yin, FTYPE *monotype, FTYPE *monosetleft, FTYPE *monosetright, FTYPE *valueleft, FTYPE *valueright)
{
  FTYPE rationorm;
  FTYPE ratios[MAXORDERS];
  int mdown,mup;
  int dmdown,dmup;
  int ddmdown,ddmup;
  FTYPE derivatives[MAXORDERS];
  FTYPE dderf[MAXORDERS];



#if(CHECKYOUTC2E)
    
#if(0)
  mdown=((yin[-1]>=pleft)&&(pleft>=yin[0])&&(yin[0]>=pright)&&(pright>=yin[1]));
  mup=((yin[-1]<=pleft)&&(pleft<=yin[0])&&(yin[0]<=pright)&&(pright<=yin[1]));
#else

#if( (CHECKYOUTC2E==1)||(CHECKYOUTC2E==3)||(CHECKYOUTC2E==4))
  // function monotonicity check
  rationorm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL)*ERRORNORM;
  ratios[0]=(yin[-1]-pleft);
  ratios[1]=(pleft-yin[0]);
  ratios[2]=(yin[0]-pright);
  ratios[3]=(pright-yin[1]);
    
  mdown=       (ratios[0]>-rationorm);
  mdown=mdown&&(ratios[1]>-rationorm);
  mdown=mdown&&(ratios[2]>-rationorm);
  mdown=mdown&&(ratios[3]>-rationorm);
    
  mup=     (ratios[0]<rationorm);
  mup=mup&&(ratios[1]<rationorm);
  mup=mup&&(ratios[2]<rationorm);
  mup=mup&&(ratios[3]<rationorm);

#endif

#if( (CHECKYOUTC2E==2)||(CHECKYOUTC2E==3)||(CHECKYOUTC2E==4))

  // derivative of function monotonicity check
  derivatives[0]=(pleft-yin[-1]);
  derivatives[1]=(yin[0]-pleft);
  derivatives[2]=(pright-yin[0]);
  derivatives[3]=(yin[1]-pright);

  //  dualfprintf(fail_file,"i=%d :: yin[-1]=%21.15g pleft=%21.15g yin[0]=%21.15g pright=%21.15g yin[1]=%21.15g\n",i,yin[-1],pleft,yin[0],pright,yin[1]);

  rationorm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2] + SMALL))*ERRORNORM;
  // some problem with using the below (when values are 0 get derivatives that are always nonzero for der[0] and der[1] for pl=1 test5 nstep=0
  //  rationorm=(fabs(derivatives[0])+fabs(derivatives[1])+fabs(derivatives[2])+fabs(derivatives[3])+SMALL)*ERRORNORM;
  ratios[0]=(derivatives[0]-derivatives[1]);
  ratios[1]=(derivatives[1]-derivatives[2]);
  ratios[2]=(derivatives[2]-derivatives[3]);
  ratios[3]=-(yin[1]-2.0*yin[0]+yin[-1])*0.25;  //atch derivative left to right!!
  
  //  dualfprintf(fail_file,"der[0]=%21.15g der[1]=%21.15g der[2]=%21.15g der[3]=%21.15g\n",derivatives[0],derivatives[1],derivatives[2],derivatives[3]);

  //  dualfprintf(fail_file,"ratios[0]=%21.15g ratios[1]=%21.15g ratios[2]=%21.15g rationorm=%21.15g\n",ratios[0],ratios[1],ratios[2],rationorm);
    
  dmdown=        (ratios[0]>-rationorm);
  dmdown=dmdown&&(ratios[1]>-rationorm);
  dmdown=dmdown&&(ratios[2]>-rationorm);
  dmdown=dmdown&&(ratios[3]>-rationorm);

    
  dmup=      (ratios[0]<rationorm);
  dmup=dmup&&(ratios[1]<rationorm);
  dmup=dmup&&(ratios[2]<rationorm);
  dmup=dmup&&(ratios[3]<rationorm);

#endif // end if checkyout==2

#if(CHECKYOUTC2E==4)

  // derivative of function monotonicity check
  dderf[0]=(derivatives[1]-derivatives[0]);
  dderf[1]=(derivatives[2]-derivatives[1]);
  dderf[2]=(derivatives[3]-derivatives[2]);

  //  dualfprintf(fail_file,"i=%d :: yin[-1]=%21.15g pleft=%21.15g yin[0]=%21.15g pright=%21.15g yin[1]=%21.15g\n",i,yin[-1],pleft,yin[0],pright,yin[1]);

  rationorm= - (fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;
  // some problem with using the below (when values are 0 get dderf that are always nonzero for der[0] and der[1] for pl=1 test5 nstep=0
  //  rationorm=(fabs(dderf[0])+fabs(dderf[1])+fabs(dderf[2])+fabs(dderf[3])+SMALL)*ERRORNORM;
  ratios[0]=(dderf[0]-dderf[1]);
  ratios[1]=(dderf[1]-dderf[2]);
  ratios[2]=- (0.5*0.125*( (yin[1]-3.0*yin[0]+3.0*yin[-1]-yin[-2])+(yin[2]-3.0*yin[1]+3.0*yin[0]-yin[-1]) ) ); //atch 3rd derivative left to right
  
  //  dualfprintf(fail_file,"der[0]=%21.15g der[1]=%21.15g der[2]=%21.15g der[3]=%21.15g\n",dderf[0],dderf[1],dderf[2],dderf[3]);

  //  dualfprintf(fail_file,"ratios[0]=%21.15g ratios[1]=%21.15g ratios[2]=%21.15g rationorm=%21.15g\n",ratios[0],ratios[1],ratios[2],rationorm);
    
  ddmdown=        (ratios[0]>-rationorm);
  ddmdown=ddmdown&&(ratios[1]>-rationorm);
  ddmdown=ddmdown&&(ratios[2]>-rationorm);

    
  ddmup=      (ratios[0]<rationorm);
  ddmup=ddmup&&(ratios[1]<rationorm);
  ddmup=ddmup&&(ratios[2]<rationorm);

#endif // end if checkyout==2
    
#endif // end else if using new machine error non-trigger version

#endif // end if CHECKYOUTC2E>0


  if(CHECKYOUTC2EIF){

    // Tells Sasha to set all weights to same #
    *monotype=1;
    *monosetleft=1;
    *monosetright=1;
    *valueleft=pleft;
    *valueright=pright;
#if(DEBUGMONONEW)
    dualfprintf(fail_file,"MONO: i=%d pl=%d : %21.15g %21.15g pleft=%21.15g %21.15g pright=%21.15g %21.15g %21.15g\n",i,pl,yin[i-2],yin[i-1],pleft,yin[i],pright,yin[i+1],yin[i+2]);
#endif
  }
  else{
    // assume could be set some other way and default is chosen
    *monotype=0;// ambiguous
    //    *monosetleft=*monosetright=0;
    
#if(DEBUGMONONEW)
    //   dualfprintf(fail_file,"C2ENOTMONO: i=%d pl=%d : monoup=%d monodown=%d :: %21.15g %21.15g pleft=%21.15g %21.15g pright=%21.15g %21.15g %21.15g\n",i,pl,monoup,monodown,yin[i-2],yin[i-1],pleft,yin[i],pright,yin[i+1],yin[i+2]);
    dualfprintf(fail_file,"C2ENOTMONO: mup=%d mdown=%d dmup=%d dmdown=%d :: %21.15g %21.15g  %21.15g %21.15g\n",mup,mdown,dmup,dmdown,yin[-1]-pleft,pleft-yin[0],yin[0]-pright,pright-yin[1]);
#endif

  }

}


// check monotonicity of final interpolated values
void check_mono_c2e_output_eno4(int i, int preforder, FTYPE *df, FTYPE *ddf, FTYPE pleft, FTYPE *yin, FTYPE (*monoindicator)[NBIGM], FTYPE (*yout)[NBIGM])
{
  FTYPE rationorm;
  FTYPE ratios[MAXORDERS];
  int mdown,mup;
  int dmdown,dmup;
  int ddmdown,ddmup;
  FTYPE derivatives[MAXORDERS];
  FTYPE dderf[MAXORDERS];



#if(CHECKYOUTC2E)
    
#if(0)
  mdown=((yin[-2]>=yin[-1])&&(yin[-1]>=pleft)&&(pleft>=yin[0])&&(yin[0]>=yin[1]));
  mup=((yin[-2]<=yin[-1])&&(yin[-1]<=pleft)&&(pleft<=yin[0])&&(yin[0]<=yin[1]));
#else

#if( (CHECKYOUTC2E==1)||(CHECKYOUTC2E==3)||(CHECKYOUTC2E==4))
  // function monotonicity check
  rationorm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;
  ratios[0]=(yin[-2]-yin[-1]);
  ratios[1]=(yin[-1]-pleft);
  ratios[2]=(pleft-yin[0]);
  ratios[3]=(yin[0]-yin[1]);
    
  mdown=       (ratios[0]>-rationorm);
  mdown=mdown&&(ratios[1]>-rationorm);
  mdown=mdown&&(ratios[2]>-rationorm);
  mdown=mdown&&(ratios[3]>-rationorm);
    
  mup=     (ratios[0]<rationorm);
  mup=mup&&(ratios[1]<rationorm);
  mup=mup&&(ratios[2]<rationorm);
  mup=mup&&(ratios[3]<rationorm);

#endif

#if( (CHECKYOUTC2E==2)||(CHECKYOUTC2E==3)||(CHECKYOUTC2E==4))

  // derivative of function monotonicity check
  derivatives[0]=(yin[-1]-yin[-2]);
  derivatives[1]=(pleft-yin[-1]);
  derivatives[2]=(yin[0]-pleft);
  derivatives[3]=(yin[1]-yin[0]);

  //  dualfprintf(fail_file,"i=%d :: yin[-2]=%21.15g yin[-1]=%21.15g pleft=%21.15g yin[0]=%21.15g yin[1]=%21.15g\n",i,yin[-2],yin[-1],pleft,yin[0],yin[1]);

  rationorm= - (fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;
  // some problem with using the below (when values are 0 get derivatives that are always nonzero for der[0] and der[1] for pl=1 test5 nstep=0
  //  rationorm=(fabs(derivatives[0])+fabs(derivatives[1])+fabs(derivatives[2])+fabs(derivatives[3])+SMALL)*ERRORNORM;
  ratios[0]=(derivatives[0]-derivatives[1]);
  ratios[1]=(derivatives[1]-derivatives[2]);
  ratios[2]=(derivatives[2]-derivatives[3]);
  
  //  dualfprintf(fail_file,"der[0]=%21.15g der[1]=%21.15g der[2]=%21.15g der[3]=%21.15g\n",derivatives[0],derivatives[1],derivatives[2],derivatives[3]);

  //  dualfprintf(fail_file,"ratios[0]=%21.15g ratios[1]=%21.15g ratios[2]=%21.15g rationorm=%21.15g\n",ratios[0],ratios[1],ratios[2],rationorm);
    
  dmdown=        (ratios[0]>-rationorm);
  dmdown=dmdown&&(ratios[1]>-rationorm);
  dmdown=dmdown&&(ratios[2]>-rationorm);

    
  dmup=      (ratios[0]<rationorm);
  dmup=dmup&&(ratios[1]<rationorm);
  dmup=dmup&&(ratios[2]<rationorm);

#endif // end if checkyout==2

#if(CHECKYOUTC2E==4)

  // derivative of function monotonicity check
  dderf[0]=(derivatives[1]-derivatives[0]);
  dderf[1]=(derivatives[2]-derivatives[1]);
  dderf[2]=(derivatives[3]-derivatives[2]);

  //  dualfprintf(fail_file,"i=%d :: yin[-1]=%21.15g pleft=%21.15g yin[0]=%21.15g pright=%21.15g yin[1]=%21.15g\n",i,yin[-1],pleft,yin[0],pright,yin[1]);

  rationorm= - (fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;
  // some problem with using the below (when values are 0 get dderf that are always nonzero for der[0] and der[1] for pl=1 test5 nstep=0
  //  rationorm=(fabs(dderf[0])+fabs(dderf[1])+fabs(dderf[2])+fabs(dderf[3])+SMALL)*ERRORNORM;
  ratios[0]=(dderf[0]-dderf[1]);
  ratios[1]=(dderf[1]-dderf[2]);
  
  //  dualfprintf(fail_file,"der[0]=%21.15g der[1]=%21.15g der[2]=%21.15g der[3]=%21.15g\n",dderf[0],dderf[1],dderf[2],dderf[3]);

  //  dualfprintf(fail_file,"ratios[0]=%21.15g ratios[1]=%21.15g ratios[2]=%21.15g rationorm=%21.15g\n",ratios[0],ratios[1],ratios[2],rationorm);
    
  ddmdown=        (ratios[0]>-rationorm);
  ddmdown=ddmdown&&(ratios[1]>-rationorm);

    
  ddmup=      (ratios[0]<rationorm);
  ddmup=ddmup&&(ratios[1]<rationorm);

#endif // end if checkyout==2

    
#endif // end else if using new machine error non-trigger version

#endif // end if CHECKYOUTC2E>0


  if(CHECKYOUTC2EIF){

    // Tells Sasha to set all weights to same #
    monoindicator[MONOINDTYPE][i]=1;
    monoindicator[MONOLEFTSET][i]=1;
    monoindicator[2][i-1]=1;
    yout[0][i]=yout[1][i-1]=pleft;
    //    dualfprintf(fail_file,"i=%d pleft=%21.15g\n",i,pleft);
  }
  else{
    monoindicator[MONOINDTYPE][i]=0;
  }

}







#if(CHECKYOUTA2CC2A==0)

// don't question, just accept
#define CHECKYOUTA2CC2AIF (1)

#elif(CHECKYOUTA2CC2A==1)

#define CHECKYOUTA2CC2AIF (mup||mdown)

#elif(CHECKYOUTA2CC2A==2)

#define CHECKYOUTA2CC2AIF (dmup||dmdown)

#elif((CHECKYOUTA2CC2A==3)||(CHECKYOUTA2CC2A==4))

#define CHECKYOUTA2CC2AIF ((mup||mdown)&&(dmup||dmdown)&&(ddmup||ddmdown))

#endif


// check monotonicity of final interpolated values
void check_mono_a2c_c2a_output(int preforder, FTYPE *df, FTYPE *ddf, FTYPE *pcenter, FTYPE *yin, FTYPE *monotype, FTYPE *monoset, FTYPE *value)
{
  FTYPE rationorm;
  FTYPE ratios[MAXORDERS];
  int mdown,mup;
  int dmdown,dmup;
  int ddmdown,ddmup;
  FTYPE fractionalchange;
  FTYPE derivatives[MAXORDERS];
  FTYPE dder[MAXORDERS];

  extern FTYPE limit_ac_correction( int order, int pl, int bs, int bf, FTYPE max_frac_difference, FTYPE *yin, FTYPE *yout );  //atch debug changed return type from int to FTYPE

#if(CHECKYOUTA2CC2A)

    
#if(0)
  mdown=((yin[-1]>=(*pcenter))&&((*pcenter)>=yin[1]));
  mup=((yin[-1]>=(*pcenter))&&((*pcenter)>=yin[1]));
#else

#if( (CHECKYOUTA2CC2A==1)||(CHECKYOUTA2CC2A==3)||(CHECKYOUTA2CC2A==4))
  // function monotonicity check, i.e. check the positivity/negativity of 1st derivative
  rationorm= - (fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL)*ERRORNORM;
  ratios[0]=(yin[-1] - yin[-2]);
  ratios[1]=((*pcenter) - yin[-1] );
  ratios[2]=(yin[1] - (*pcenter));
  ratios[3]=(yin[2] - yin[1]);
    
  // larger than some positive number
  mdown=       (ratios[0]>-rationorm);
  mdown=mdown&&(ratios[1]>-rationorm);
  mdown=mdown&&(ratios[2]>-rationorm);
  mdown=mdown&&(ratios[3]>-rationorm);

  // smaller than some negative number
  mup=     (ratios[0]<rationorm);
  mup=mup&&(ratios[1]<rationorm);
  mup=mup&&(ratios[2]<rationorm);
  mup=mup&&(ratios[3]<rationorm);
#endif

#if( (CHECKYOUTA2CC2A==2)||(CHECKYOUTA2CC2A==3)||(CHECKYOUTA2CC2A==4))

  // derivative of function monotonicity check
  derivatives[0]= ratios[1] - ratios[0];
  derivatives[1]= ratios[2] - ratios[1];
  derivatives[2]= ratios[3] - ratios[2];
  
  // just check sign didn't change
  //dmdown=dmup = (fabs(derivatives[0]-derivatives[1])<ERRORNORM) || (sign(derivatives[0])==sign(derivatives[1]) );

  dmdown=        (derivatives[0]>-rationorm);
  dmdown=dmdown&&(derivatives[1]>-rationorm);
  dmdown=dmdown&&(derivatives[2]>-rationorm);

  dmup=      (derivatives[0]<rationorm);
  dmup=dmup&&(derivatives[1]<rationorm);
  dmup=dmup&&(derivatives[2]<rationorm);

  //3rd der monotonicity check
  dder[0] = derivatives[1] - derivatives[0];
  dder[1] = derivatives[2] - derivatives[1];
 
  ddmdown=         (dder[0]>-rationorm);
  ddmdown=ddmdown&&(dder[1]>-rationorm);

  ddmup=       (dder[0]<rationorm);
  ddmup=ddmup&&(dder[1]<rationorm);


#endif // end if checkyout==2


   
#endif // end else if using new machine error non-trigger version



#endif // end if CHECKYOUT>0

  //  dualfprintf(fail_file,"CHECKYOUTA2CC2A=%d mup=%d mdown=%d dm=%d fc=%21.15g\n",CHECKYOUTA2CC2A,mup,mdown,dmdown,fractionalchange);


  if(CHECKYOUTA2CC2AIF){

#if(CHECKYOUTA2CC2A==4)
    limit_ac_correction(preforder,0,0,0,FRACCHANGEALLOWED,&yin[0],pcenter);
    //  fractionalchange=fabs(yin[0]-(*pcenter))/(fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+SMALL);
#endif

    // Tells Sasha to set all weights to same #
    *monotype=1;
    *monoset=1;
    *value=*pcenter;
    //    *value=yin[0];
    //    *monotype=*monoset=0;

  }
  else{
    // assume default chosen is ok or set some other way, so don't want to reset
    *monotype=0;// ambiguous
    //    *monoset=0;
  }

}


//check for cusp
// too sharp an indicator
int check_for_cusp_old(FTYPE *yin, FTYPE *df, FTYPE *ddf)
{
  FTYPE norm;

  norm=(fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL);
  norm=norm*norm*ERRORNORM;
  

  if( (ddf[-1] * df[0] > norm) && (ddf[0] * df[0] < -norm) ){
    
    if(       df[2] * df[0] < -norm ||  df[1] * df[0] < -norm ) return(1);

  }

  if( (ddf[1] * ( -df[1] ) > norm) && (ddf[0] * ( -df[1] ) < -norm) ){

    if(      ( -df[-1] ) * ( -df[1] ) < -norm  ||   ( -df[0]  ) * ( -df[1] ) < -norm  ) return(1);

  }

  return(0);
}

//check for cusp (cleaned code)
int check_for_cusp(FTYPE *yin, FTYPE *df, FTYPE *ddf)
{
  FTYPE norm;
  FTYPE f1,f2,f3,f4,f5,f6,f7,f8;

  norm=(fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL);
  norm=norm*norm*ERRORNORM;

  f1= ddf[-1] * df[0] -norm;
  f2=-ddf[0]  * df[0] -norm;
  f3= -df[2]  * df[0] -norm;
  f4= -df[1]  * df[0] -norm;

  if(f1>0 && f2>0 && (f3>0 || f4>0) ) return(1);

  f5=-ddf[1]  * df[1] -norm;
  f6= ddf[0]  * df[1] -norm;
  f7= -df[-1] * df[1] -norm;
  f8= -df[0]  * df[1] -norm;

  if(f5>0 && f6>0 && (f7>0 || f8>0) ) return(1);
  

  return(0);
}


#define BADNESSFRAC (1E4)
#define BADNESSBOTTOM (1E-6)

//check for cusp
// indicator from 0 (no cusp) to 1 (super cuspy)
//    check_for_cusp_indicator(&yin[i], &df[DFONESIDED][i], &df[DF2OFONESIDED][i], &cuspindicator);

int check_for_cusp_indicator(FTYPE *yin, FTYPE *df, FTYPE *ddf, FTYPE *cuspindicator)
{
  FTYPE norm;
  FTYPE f1,f2,f3,f4,f5,f6,f7,f8;
  FTYPE badtype1,badtype2,badness;


  norm=(fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL);
  norm=norm*norm*BADNESSBOTTOM;

  f1= ddf[-1] * df[0] -norm;
  f2=-ddf[0]  * df[0] -norm;
  f3= -df[2]  * df[0] -norm;
  f4= -df[1]  * df[0] -norm;

  f5=-ddf[1]  * df[1] -norm;
  f6= ddf[0]  * df[1] -norm;
  f7= -df[-1] * df[1] -norm;
  f8= -df[0]  * df[1] -norm;

  badtype1=(f1>0 && f2>0 && (f3>0 || f4>0) );
  badtype2=(f5>0 && f6>0 && (f7>0 || f8>0) );

  if(badtype1 || badtype2 ){
    badness=(badtype1*fabs(f1+f2+f3+f4)+badtype2*fabs(f5+f6+f7+f8))/(4.0*norm);
    if(badness>100.0){
      *cuspindicator=1.0;
      //      dualfprintf(fail_file,"1 here: %21.15g: %21.15g\n",*cuspindicator,badness);
    }
    else if(badness>1.0){
      *cuspindicator=(badness-1.0)/(100.0-1.0);
      *cuspindicator=1.0;
      //      dualfprintf(fail_file,"2 here: %21.15g: %21.15g\n",*cuspindicator,badness);
    }
    else{
      *cuspindicator=0.0;
      //      dualfprintf(fail_file,"3 here: %21.15g\n",*cuspindicator);
    }
    return(1);
  }
  else{
    //    dualfprintf(fail_file,"4 here: %21.15g\n",*cuspindicator);
    *cuspindicator=0.0;
  }
  //  dualfprintf(fail_file,"END here: %21.15g\n",*cuspindicator);

  //  *cuspindicator=0.0;

  return(0);
}

