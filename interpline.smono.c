
/*! \file interpline.smono.c
     \brief Sasha Mono Spatial Interpolation for fluxes based upon providing full 1D line information
     //Possible bugs:
     //if not reset the whole line of monoindicators to zero, asymmetries develop in the implosion test.  Reason unknown. Probably, ps & pe do not cover a large enough region
     */


#include "decs.h"
#include "reconstructeno.h"  

#include "interpline.smono_static.h"
#include "interpline.smono.h"

#define DEFAULTMONOORDER (5)

/// GODMARK: Sasha created the MONO flags in a confused way.
/// He set things up so that if DO_SMONO_A2C ==0 then still  DO_SMONO_A2C_CUSP_INDICATOR is used and determines monoindicator (left/right) and value
/// I would assume that if DO_SMONO_A2C==0, then nothing is done
void compute_smonotonicity_line(
                                int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                int *minorder, int *maxorder, int *shift,   
                                FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                FTYPE (*monoindicator)[NBIGM], FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM])
{ 
  void compute_smonotonicity_line_split(int setindicator, int setyout,
                                        int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                        int *minorder, int *maxorder, int *shift,   
                                        FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                        FTYPE (*monoindicator)[NBIGM], FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);

  // compute normal monotonicity line
  compute_smonotonicity_line_split(1,1,
                                   recontype, whichreduce, preforder, pl, bs, ps, pe, be, 
                                   minorder, maxorder, shift,   
                                   shockindicator,  df, 
                                   monoindicator, yin, yout, youtpolycoef);

}

void compute_polycoef_line(
                           int preforder, int pl, int bs, int ps, int pe, int be, 
                           FTYPE *yin, FTYPE (*youtpolycoef)[NBIGM])
{
  int i;
  int order=(preforder+1)/2;
  int polyorder;
  FTYPE polycoefarray[MAXSPACEORDER+1];

  for( i = bs+order-1; i <= be-order+1; i++ ) {  
    compute_c2a_polycoef_simple_eno( preforder, &yin[i], polycoefarray );
    //place derivatives into external array
    for( polyorder = 0; polyorder < DEFAULTMONOORDER; polyorder++ ) {
      youtpolycoef[polyorder][i] = polycoefarray[polyorder];
    }
  }  
}


/// monoindicator[0]: -1,0,1 for rough, ambiguous, and monotonic
/// monoindicator[1]: whether set cell's left value (or central value, for a2c/c2a reconstruction)
/// monoindicator[2]: whether set cell's right value
/// yout[0][i] is the left interface value for c2e (centered value for a2c & c2a) for grid cell i
/// yout[1][i] is the right interface value for c2e for grid cell i, does not make sense for a2c & c2a
void compute_smonotonicity_line_split(int setindicator, int setyout,
                                      int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                      int *minorder, int *maxorder, int *shift,   
                                      FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                      FTYPE (*monoindicator)[NBIGM], FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM])
{ 
  int i;
  FTYPE indicator;
  int check_for_cusp_new(FTYPE *yin, FTYPE *df, FTYPE *ddf);
  // int order = (preforder+1)/2;  //the length of the WENO-preforder substencil 
  // JCM: assume want mono result no matter what otherwise order is as long as high enough order so points exist for mono
  int order = (DEFAULTMONOORDER+1)/2;  //the length of the WENO-preforder substencil 
  void set_as_rough_indicator(int recontype, int i, FTYPE (*monoindicator)[NBIGM]);
  void set_as_rough(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], FTYPE (*monoindicator)[NBIGM]);
  void set_as_rough_value(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);
  int setasrough;
#if(MERGEDC2EA2CMETHOD)
  //array for storing polynomial derivatives
  int polyorder;
  FTYPE polycoefarray[MAXSPACEORDER+1];
#endif




  if(setindicator){
    if( preforder < DEFAULTMONOORDER ) {
      //wrong order for SMONO to operate on, therefore initialize the mono indicators in such a way that
      //WENO for that order will be fully used and then return
      for( i = -NBIGBND; i < NBIG + NBIGBND; i++ ){
        monoindicator[MONOINDTYPE][i] = 0.0;
        monoindicator[MONOLEFTSET][i] = 0.0;
        monoindicator[MONORIGHTSET][i] = 0.0;
      }
      return;
    }

#if( DO_ASSERTS )  //this is for checking whether we set enough monoindicators on the grid; the negative values are caught by assert()'s in reconstructeno.c
    //reset the monoindicator line to zeroes --- otherwise asymmetries develop which indicates that ps & pe are not defined to cover large enough region  SASMARK2
    if(recontype!=CVT_C2E){
      for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
        monoindicator[MONOINDTYPE][i]=-1;
        monoindicator[MONOLEFTSET][i]=-1;
      }
    }
    else{
      for(i=-NBIGBND;i<NBIG+NBIGBND;i++){
        monoindicator[MONOINDTYPE][i]=-1;
        monoindicator[MONOLEFTSET][i]=-1;
        monoindicator[MONORIGHTSET][i]=-1;
      }
    }
#endif
  }




  ////////////////////
  //
  //  AVG2CEN interp
  //
  if(recontype==CVT_A2C){
    for( i = bs+order-1; i <= be-order+1; i++ ) {  

#if( DO_SMONO_A2C_CUSP_INDICATOR )
      // check for cusp

      if(setindicator){
        if((check_for_cusp_new(&yin[i],&df[DFONESIDED][i],&df[DF2OFONESIDED][i]))  ){
          set_as_rough_indicator(recontype, i, monoindicator);
        }
        else{
          monoindicator[MONOINDTYPE][i]=0;
          monoindicator[MONOLEFTSET][i]=0;
          monoindicator[MONORIGHTSET][i]=0;
        }
      }

      if(monoindicator[MONOINDTYPE][i]==0 && monoindicator[MONOLEFTSET][i]==1 && monoindicator[MONORIGHTSET][i]==1){
        setasrough=1;
      }
      else setasrough=0;

      // assumes set_as_rough_indicator sets indicator to below settings
      if(setyout && setasrough){
        set_as_rough_value(recontype, i, yin, yout, youtpolycoef);
      }


#else
      setasrough=0;
#endif


      if(setasrough==0){

        if(setindicator){
#if( DO_SMONO_A2C ) 
          indicator = compute_mono_indicator_average_eno5( &yin[i], SQRT_WENO_EPSILON );  //note conversion double -> int; temporary: later this switch will be removed
#else
          indicator = 0;
#endif
#if( SMONO_DO_SMOOTH_TRANSITION )
          monoindicator[MONOINDTYPE][i] = indicator;
#else
          monoindicator[MONOINDTYPE][i] = (int)(indicator + 0.5);
#endif
          monoindicator[MONOLEFTSET][i] = monoindicator[MONOINDTYPE][i];

          // DEBUG
          //   if(i<0 || i>N1){
          //     monoindicator[MONOINDTYPE][i]=sqrt(-1.);
          //     monoindicator[MONOLEFTSET][i]=sqrt(-1.);
          //     monoindicator[MONORIGHTSET][i]=sqrt(-1.);
          //   }
          //   dualfprintf(fail_file,"nstep=%ld steppart=%d mono[%d]=%21.15g yin=%21.15g\n",nstep,steppart,i,monoindicator[MONOINDTYPE][i],yin[i]);
          //   monoindicator[MONOINDTYPE][i]=1;
          //      monoindicator[MONOLEFTSET][i]=1;
          //      monoindicator[MONORIGHTSET][i]=1;

        }

        if(setyout){
          //if( indicator != 0 ){
          if(monoindicator[MONOINDTYPE][i] != 0 ) {
            a2c_simple_eno( DEFAULTMONOORDER, &yin[i], &yout[0][i] ); // MERGEDC2EA2CMETHOD  TODO
          }
        }

      }// end if setasrough
    } // end for loop over i
  }
  //////////////////////
  //
  //    CENT2EDGE interp
  //
  else if(recontype == CVT_C2E) {
    for( i = bs+order-1; i <= be-order+1; i++ ) {  

#if( DO_SMONO_C2E_CUSP_INDICATOR )
      // check for cusp

      if(setindicator){
        if((check_for_cusp_new(&yin[i],&df[DFONESIDED][i],&df[DF2OFONESIDED][i]))  ){
          set_as_rough_indicator(recontype, i, monoindicator);
        }
        else{
          monoindicator[MONOINDTYPE][i]=0;
          monoindicator[MONOLEFTSET][i]=0;
          monoindicator[MONORIGHTSET][i]=0;
        }
      }

      if(monoindicator[MONOINDTYPE][i]==0 && monoindicator[MONOLEFTSET][i]==1 && monoindicator[MONORIGHTSET][i]==1){
        setasrough=1;
      }
      else setasrough=0;


      // assumes set_as_rough_indicator sets indicator to below settings
      if(setyout && setasrough){
        set_as_rough_value(recontype, i, yin, yout, youtpolycoef);
      }

#else
      setasrough=0;
#endif





      if(setasrough==0){
 
        if(setindicator){
#if( DO_SMONO_C2E ) 
          //have a smooth transition between eno5 and weno5 at the level of SQRT_WENO_EPSILON, i.e. at about the same level as that of WENO-5 switching between its stencils.
          indicator = compute_mono_indicator_point_eno5( &yin[i], SMONO_EPSILON );  
#else
          indicator = 0;
#endif
#if( SMONO_DO_SMOOTH_TRANSITION )
          monoindicator[MONOINDTYPE][i] = indicator;
#else
          monoindicator[MONOINDTYPE][i] = (int)(indicator + 0.5);
#endif
          monoindicator[MONOLEFTSET][i] = monoindicator[MONOINDTYPE][i];
          monoindicator[MONORIGHTSET][i] = monoindicator[MONOINDTYPE][i];

        }

        if(setyout){
          if(monoindicator[MONOINDTYPE][i] != 0 ){
            //if(indicator != 0 ){
            //get the left interface value
            c2e_simple_eno( DEFAULTMONOORDER, 0, &yin[i], &yout[0][i] ); // MERGEDC2EA2CMETHOD  TODO
     
            //get the right interface value
            c2e_simple_eno( DEFAULTMONOORDER, 1, &yin[i], &yout[1][i] ); // MERGEDC2EA2CMETHOD  TODO
#if(MERGEDC2EA2CMETHOD)
            compute_c2a_polycoef_simple_eno( DEFAULTMONOORDER, &yin[i], polycoefarray );
            //place derivatives into external array
            for( polyorder = 0; polyorder < DEFAULTMONOORDER; polyorder++ ) {
              youtpolycoef[pl][i];
            }
#endif
          }
        }

      }// end else if smooth
    }
  }
  /////////////////////
  //
  //     CENT2AVG interp
  //
  else if(recontype == CVT_C2A ) {
    for( i = bs+order-1; i <= be-order+1; i++ ) {  
#if( DO_SMONO_C2A_CUSP_INDICATOR )
      // check for cusp


      if(setindicator){
        if((check_for_cusp_new(&yin[i],&df[DFONESIDED][i],&df[DF2OFONESIDED][i]))  ){
          //   dualfprintf(fail_file,"SET ROUGH i=%d\n",i);
          set_as_rough_indicator(recontype, i, monoindicator);
        }
        else{
          monoindicator[MONOINDTYPE][i]=0;
          monoindicator[MONOLEFTSET][i]=0;
          monoindicator[MONORIGHTSET][i]=0;
        }
      }

      if(monoindicator[MONOINDTYPE][i]==0 && monoindicator[MONOLEFTSET][i]==1 && monoindicator[MONORIGHTSET][i]==1){
        setasrough=1;
      }
      else setasrough=0;

      // assumes set_as_rough_indicator sets indicator to below settings
      if(setyout && setasrough){
        set_as_rough_value(recontype, i, yin, yout, youtpolycoef);
      }


#else
      setasrough=0;
#endif




#if(0)
      if(crapdebug||1){
        dualfprintf(fail_file,"i=%d monotype=%21.15g monoleft=%21.15g monoright=%21.15g\n",i,monoindicator[MONOINDTYPE][i],monoindicator[MONOLEFTSET][i],monoindicator[MONORIGHTSET][i]);
      }
#endif


      if(setasrough==0){

        if(setindicator){
#if( DO_SMONO_C2A ) 
          //have a smooth transition between eno5 and weno5 at the level of SQRT_WENO_EPSILON, i.e. at about the same level as that of WENO-5 switching between its stencils.
          indicator = compute_mono_indicator_point_eno5( &yin[i], SMONO_EPSILON );  //note conversion double -> int; temporary: later this switch will be removed
#else
          indicator = 0;
#endif
#if( SMONO_DO_SMOOTH_TRANSITION )
          monoindicator[MONOINDTYPE][i] = indicator;
#else
          monoindicator[MONOINDTYPE][i] = (int)(indicator + 0.5);
#endif
          monoindicator[MONOLEFTSET][i] = monoindicator[MONOINDTYPE][i];
  

        }

        if(setyout){
          if( monoindicator[MONOINDTYPE][i] != 0 ){
            //if( indicator != 0 ){
            c2a_simple_eno( DEFAULTMONOORDER, &yin[i], &yout[0][i] ); // MERGEDC2EA2CMETHOD  TODO
          }
        }

      }// end else if smooth
    }//end for
  }

  ////check symmetry of input
  //for( i = ps; i <= pe; i++ ) {
  // if( yin[i] != yin[N1 - 1 - i] * ((pl==U1)?(-1):(1)) ) {
  //  dualfprintf( fail_file, "Asymmetry in yin: yin[%d] = %21.15g, yin[%d] = %21.15g\n",
  //   i, yin[i], N1 - 1 - i, yin[N1 - 1 - i]);
  // }
  //}

  ////check symmetry of monoindicators
  //for( i = ps; i <= pe; i++ ) {
  // if( monoindicator[MONOYIN][i] != monoindicator[MONOYIN][N1 - 1 - i] ) {
  //  dualfprintf( fail_file, "Asymmetry in monoindicator: monoindicator[MONOYIN][%d] = %21.15g, monoindicator[MONOYIN][%d] = %21.15g\n",
  //   i, monoindicator[MONOYIN][i], N1 - 1 - i, monoindicator[MONOYIN][N1 - 1 - i]);
  // }
  //}
}









void set_as_rough(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], FTYPE (*monoindicator)[NBIGM])
{
  void set_as_rough_indicator(int recontype, int i, FTYPE (*monoindicator)[NBIGM]);
  void set_as_rough_value(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);

  set_as_rough_indicator(recontype, i, monoindicator);
  set_as_rough_value(recontype, i,yin, yout, youtpolycoef);

}

/// if not smooth then reduce to simple limiter
void set_as_rough_indicator(int recontype, int i, FTYPE (*monoindicator)[NBIGM])
{

  // set as rough
  //monoindicator[MONOINDTYPE][i]=-1;
  // say unknown roughness
  monoindicator[MONOINDTYPE][i]=0;
  // say already set
  monoindicator[MONOLEFTSET][i]=1;
  monoindicator[MONORIGHTSET][i]=1;      
  //      dualfprintf(fail_file,"ROUGH i=%d\n",i);

}


/// if not smooth then reduce to simple limiter
void set_as_rough_value(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM])
{
  extern void c2e_simple_limiter(int WHICHLIMITERTOREDUCETO, FTYPE *yin, FTYPE *valueleft, FTYPE *valueright);

  //dualfprintf(fail_file,"shocki[%d]=%21.15g\n",i,shockindicator[EOMSETMHD][i]);
  //if((shockindicator[EOMSETMHD][i]>=0.1)||(check_for_cusp_new(&df[DFONESIDED][i],&df[DF2OFONESIDED][i],&yin[i]))  ){// does somewhat better than not checking
  //if(0){

  // MERGEDC2EA2CMETHOD  TODO
  if(recontype==CVT_C2E) c2e_simple_limiter(MINM, &yin[i], &yout[0][i],&yout[1][i]);
  // otherwise do nothing (no a2c/c2a conversion)
  // MERGEDC2EA2CMETHOD  TODO
  else yout[0][i]=yin[i];

}

///fits a polynomial through the yin[] averages and returns a value between 0 & 1:
/// 0 - the polynomial fit is non-monotonic or about (within epsilon * ||f||) to become non-monotonic (hence do the full WENO-5 routine for this grid cell)
/// 1 - the polynomial fit is safely monotonic, so can set the WENO-5 weights to be equal to simple weights
/// between 0 & 1 - the larger the number is, the closer the WENO weights should be to the optimal ones.  So, do:
///                 new_unoptimized_weight = mono_indicator * simple_weight + (1-monoindicator) * old_unoptimized_weight.
///       NOTE:  only change the unoptimized weights since the weights are req-d to be changed only for the purposes of stencil reduction.
///              The actual summing up of the full-stencil-based reconstruction and the WENO-5 reconstruction is to be performed in eno_cvt.
FTYPE compute_mono_indicator_average_eno5( FTYPE *yin, FTYPE epsilon )
{
#define MAX_DERA2C 3
  FTYPE left_der[MAX_DERA2C+1], right_der[MAX_DERA2C+1], minmod_der[MAX_DERA2C+1];
  FTYPE minabs_der = BIG;
  FTYPE norm;
  FTYPE monoindicator;
  FTYPE cutoff_value;
  int n, maxn;  //order of derivative and max order of derivative whose monotonicity is checked
  FTYPE der2_crit;
  FTYPE acoeff, bcoeff, ccoeff;

  //SASMARK: need to symmetrize the expressions for the derivatives below

  //values of derivatives at i = -1, i.e. the left side of the stencil
  //left_der[1] = ((-9*yin[-2] + 5*yin[2]) + (- 50*yin[-1] - 30*yin[1]) + 84*yin[0])/48.;
  //left_der[2] = ((7*yin[-2] - yin[2]) + (- 12*yin[-1] + 4*yin[1]) + 2*yin[0])/8.;
  //left_der[3] = ((-3*yin[-2] - yin[2]) + (10*yin[-1] + 6*yin[1]) - 12*yin[0])/2.;
  //unsymmetrized expressions at -1
  //left_der[1] = (-9*yin[-2] - 50*yin[-1] + 84*yin[0] - 30*yin[1] + 5*yin[2])/48.;
  //left_der[2] = (7*yin[-2] - 12*yin[-1] + 2*yin[0] + 4*yin[1] - yin[2])/8.;
  //left_der[3] = (-3*yin[-2] + 10*yin[-1] - 12*yin[0] + 6*yin[1] - yin[2])/2.;

  //values of derivatives at i = 1, i.e. the right side of the stencil
  //right_der[1] = ((-5*yin[-2] + 9*yin[2]) + (30*yin[-1] + 50*yin[1]) - 84*yin[0])/48.;
  //right_der[2] = (-yin[-2] + 7*yin[2] + (4*yin[-1] - 12*yin[1]) + 2*yin[0])/8.;
  //right_der[3] = ((yin[-2] + 3*yin[2]) + (- 6*yin[-1] - 10*yin[1]) + 12*yin[0])/2.;
  //unsymmetrized expressions at 1
  //right_der[1] = (-5*yin[-2] + 30*yin[-1] - 84*yin[0] + 50*yin[1] + 9*yin[2])/48.;
  //right_der[2] = (-yin[-2] + 4*yin[-1] + 2*yin[0] - 12*yin[1] + 7*yin[2])/8.;
  //right_der[3] = (yin[-2] - 6*yin[-1] + 12*yin[0] - 10*yin[1] + 3*yin[2])/2.;

  //left at -2
  left_der[1] = (-95*yin[-2] + 174*yin[-1] - 120*yin[0] + 50*yin[1] - 9*yin[2])/48.;
  left_der[2] = (23*yin[-2] - 68*yin[-1] + 74*yin[0] - 36*yin[1] + 7*yin[2])/8.;
  left_der[3] = (-5*yin[-2])/2. + 9*yin[-1] - 12*yin[0] + 7*yin[1] - (3*yin[2])/2.;

  //right at +2
  right_der[1] = (9*yin[-2] - 50*yin[-1] + 120*yin[0] - 174*yin[1] + 95*yin[2])/48.;
  right_der[2] = (7*yin[-2] - 36*yin[-1] + 74*yin[0] - 68*yin[1] + 23*yin[2])/8.;
  right_der[3] = (3*yin[-2])/2. - 7*yin[-1] + 12*yin[0] - 9*yin[1] + (5*yin[2])/2.;

#if( DO_MONO_1ST_DERIVATIVE)
  left_der[1] = MINMOD( left_der[1], yin[-1] - yin[-2] );
  left_der[1] = MINMOD( left_der[1], yin[0] - yin[-1] );
  right_der[1] = MINMOD( right_der[1], yin[2] - yin[1] );
  right_der[1] = MINMOD( right_der[1], yin[1] - yin[0] );
#endif

  //always check the monotonicity of the 3rd derivative for a2c
  maxn = MAX_DERA2C;
  //  maxn = 1;// DEBUG

  //Check the monotonicity of all derivatives of the function. Start with the 3rd order derivative
  //(which is a linear function for 5th order polynomial) and check if it has the same 
  //sign at the ends of the interval. Store the smaller of the two derivatives in absolute value; 
  //store 0.0 for the value if the signs are different.
  //Repeat the same for 2nd and 1st order derivatives and keep track of the minimum of the 
  //absolute value of the derivatives. This value is then used as an indicator of how close
  //the function is to being non-monotonic
  for( n = 1; n <= maxn; n++ ) {
    minmod_der[n] = MINMOD( left_der[n], right_der[n] );
    //keep track of the minimum value of the derivatives on either side of the interval;
    //if the signs of the derivatives differ, set the min value to 0
    minabs_der = MIN( minabs_der, fabs(minmod_der[n]) );  
  }

  //compute the norm of the function in a symmetric way
  norm = fabs(yin[0]) + ( fabs(yin[-1]) + fabs(yin[1]) ) 
    + ( fabs(yin[-2]) + fabs(yin[2]) );

  //cutoff value: if an abs. value of the derivative is smaller than this value, the derivative is considered = 0.
  cutoff_value = norm * epsilon + MINNUMREPRESENT;

  //smooth transition from 0 to 1 depending on the value of minabs_der
  monoindicator = transition_function( minabs_der, 0.0, cutoff_value );

  //if( minabs_der != 0.0 && norm > 1 ) {
  // dualfprintf( fail_file, "monoindicator = %21.15g\n", monoindicator );
  //}

  return( monoindicator );
}

///fits a polynomial through the yin[] points and returns a value between 0 & 1:
/// 0 - the polynomial fit is non-monotonic or about (within epsilon * ||f||) to become non-monotonic (hence do the full WENO-5 routine for this grid cell)
/// 1 - the polynomial fit is safely monotonic, so can set the WENO-5 weights to be equal to simple weights
/// between 0 & 1 - the larger the number is, the closer the WENO weights should be to the optimal ones.  So, do:
///                 new_unoptimized_weight = mono_indicator * simple_weight + (1-monoindicator) * old_unoptimized_weight.
///       NOTE:  only change the unoptimized weights since the weights are req-d to be changed only for the purposes of stencil reduction.
///              The actual summing up of the full-stencil-based reconstruction and the WENO-5 reconstruction is to be performed in eno_cvt.
FTYPE compute_mono_indicator_point_eno5( FTYPE *yin, FTYPE epsilon )
{
#define MAX_DERC2E 3
  FTYPE left_der[MAX_DERC2E+1], right_der[MAX_DERC2E+1], minmod_der[MAX_DERC2E+1];
  FTYPE minabs_der = BIG;
  FTYPE norm;
  FTYPE monoindicator;
  FTYPE cutoff_value;
  int n, maxn;  //order of derivative and max order of derivative whose monotonicity is checked
  FTYPE der2_crit;
  FTYPE acoeff, bcoeff, ccoeff;

  //SASMARK: need to symmetrize the expressions for the derivatives below

  //values of derivatives at i = -2, i.e. the left side of the stencil
  //left_der[1] = (-25*yin[-2])/12. + 4*yin[-1] - 3*yin[0] + (4*yin[1])/3. - yin[2]/4.;
  //left_der[2] = (35*yin[-2] - 104*yin[-1] + 114*yin[0] - 56*yin[1] + 11*yin[2])/12.;
  //left_der[3] = (-5*yin[-2])/2. + 9*yin[-1] - 12*yin[0] + 7*yin[1] - (3*yin[2])/2.;
  //values of derivatives at i = 2, i.e. the right side of the stencil
  //right_der[1] = yin[-2]/4. - (4*yin[-1])/3. + 3*yin[0] - 4*yin[1] + (25*yin[2])/12.;
  //right_der[2] = (11*yin[-2] - 56*yin[-1] + 114*yin[0] - 104*yin[1] + 35*yin[2])/12.;
  //right_der[3] = (3*yin[-2])/2. - 7*yin[-1] + 12*yin[0] - 9*yin[1] + (5*yin[2])/2.;


  //at -1 & 1
  left_der[1] = ( ((-3*yin[-2] + yin[2]) + (- 10*yin[-1] - 6*yin[1])) + 18*yin[0] )/12.;
  left_der[2] = ( ((11*yin[-2] - yin[2]) + (- 20*yin[-1] + 4*yin[1])) + 6*yin[0] )/12.;
  left_der[3] = ( ((-3*yin[-2] - yin[2]) + (10*yin[-1] + 6*yin[1])) - 12*yin[0])/2.;

  right_der[1] = ( ((-yin[-2] + 3*yin[2]) + (6*yin[-1] + 10*yin[1])) - 18*yin[0] )/12.;
  right_der[2] = ( ((-yin[-2] + 11*yin[2]) + (4*yin[-1] - 20*yin[1])) + 6*yin[0] )/12.;
  right_der[3] = ( ((yin[-2] + 3*yin[2]) - (6*yin[-1] + 10*yin[1])) + 12*yin[0])/2.;

  //at -1/2
  //left_der[1] = (yin[-2] - 27*yin[-1] + 27*yin[0] - yin[1])/24.;
  //left_der[2] = (7*yin[-2] + 8*yin[-1] - 42*yin[0] + 32*yin[1] - 5*yin[2])/24.;
  //left_der[3] = -yin[-2] + 3*yin[-1] - 3*yin[0] + yin[1];

  //at 1/2
  //left_der[1] = (yin[-1] - 27*yin[0] + 27*yin[1] - yin[2])/24.;
  //left_der[2] = (-5*yin[-2] + 32*yin[-1] - 42*yin[0] + 8*yin[1] + 7*yin[2])/24.;
  //left_der[3] = -yin[-1] + 3*yin[0] - 3*yin[1] + yin[2];

  // again ders at -1 and 1 for sanity check:
  //left_der[1] = (-3*yin[-2] - 10*yin[-1] + 18*yin[0] - 6*yin[1] + yin[2])/12.;
  //left_der[2] = (11*yin[-2] - 20*yin[-1] + 6*yin[0] + 4*yin[1] - yin[2])/12.;
  //left_der[3] = (-3*yin[-2] + 10*yin[-1] - 12*yin[0] + 6*yin[1] - yin[2])/2.;

  //right_der[1] = (-yin[-2] + 6*yin[-1] - 18*yin[0] + 10*yin[1] + 3*yin[2])/12.;
  //right_der[2] = (-yin[-2] + 4*yin[-1] + 6*yin[0] - 20*yin[1] + 11*yin[2])/12.;
  //right_der[3] = (yin[-2] - 6*yin[-1] + 12*yin[0] - 10*yin[1] + 3*yin[2])/2.;
 

#if( DO_MONO_1ST_DERIVATIVE ) 
  left_der[1] = MINMOD( left_der[1], yin[-1] - yin[-2] );
  right_der[1] = MINMOD( right_der[1], yin[2] - yin[1] );

  //left_der[1] = MINMOD( left_der[1], yin[-0] - yin[-1] );
  //right_der[1] = MINMOD( right_der[1], yin[1] - yin[0] );
#endif

  //whether to check the monotonicity of the 3rd derivative
#if( DO_3RD_DER )
  maxn = MAX_DERC2E;

  //at -1/2
  //left_der[3]  = (-yin[-2] + yin[1]) + (3*yin[-1] - 3*yin[0]);
  //at 1/2
  //right_der[3] = (-yin[-1] + yin[2]) + (3*yin[0] - 3*yin[1]);
  //at -1
  //left_der[3] = (-3*yin[-2] + 10*yin[-1] - 12*yin[0] + 6*yin[1] - yin[2])/2.;
  //at 1
  //right_der[3] = (yin[-2] - 6*yin[-1] + 12*yin[0] - 10*yin[1] + 3*yin[2])/2.;
#else
  maxn = MAX_DERC2E - 1;
#endif

  //Make sure that 2nd der does not change sign
  //in the middle of the interval where it can have a min/max

  //For eno-5 2nd derivative is a parabola:  f'' = a x^2 + b x + c, make sure that if its minimum/maximum is within
  //the interval [-1, 1], the value there is of the same sign as at the ends of the interval:
  //acoeff = ( (yin[-2] + yin[2]) - (4*yin[-1] + 4*yin[1]) + 6*yin[0])/2.;
  //bcoeff = ( (-yin[-2] + yin[2]) + (2*yin[-1] - 2*yin[1]) )/2.;
  //ccoeff = (-yin[-2]/12. - yin[2]/12.) + 4*(yin[-1] + yin[1])/3. - (5*yin[0])/2.;

  //acoeff = (yin[-2] - 4*yin[-1] + 6*yin[0] - 4*yin[1] + yin[2])/2.;
  //bcoeff = (-yin[-2] + 2*yin[-1] - 2*yin[1] + yin[2])/2.;
  //ccoeff = -yin[-2]/12. + (4*yin[-1])/3. - (5*yin[0])/2. + (4*yin[1])/3. - yin[2]/12.;


  //if( fabs(bcoeff) < 2. * fabs(acoeff) ) {  // x_crit = -b/2a, |x_crit| < 1 <=> |b| < 2 |a|
  // //this also makes sure that acoeff != 0, so can divide by it
  // der2_crit = - bcoeff * bcoeff / (4. * acoeff) + ccoeff; //f_crit = - b^2/(4a) + c = -D/(4a)
  // left_der[2] = MINMOD( der2_crit, left_der[2] );  //make sure that near min/max, 2nd der is also of the same sign as at the ends of the interval
  //}


  //Check the monotonicity of all derivatives of the function. Start with the 3rd order derivative
  //(which is a linear function for 5th order polynomial) and check if it has the same 
  //sign at the ends of the interval. Store the smaller of the two derivatives in absolute value; 
  //store 0.0 for the value if the signs are different.
  //Repeat the same for 2nd and 1st order derivatives and keep track of the minimum of the 
  //absolute value of the derivatives. This value is then used as an indicator of how close
  //the function is to being non-monotonic
  for( n = 1; n <= maxn; n++ ) {
    minmod_der[n] = MINMOD( left_der[n], right_der[n] );
    //keep track of the minimum value of the derivatives on either side of the interval;
    //if the signs of the derivatives differ, set the min value to 0
    minabs_der = MIN( minabs_der, fabs(minmod_der[n]) );  
  }

  //compute the norm of the function in a symmetric way
  norm =  fabs(yin[0]) + 
    ( (fabs(yin[-1]) + fabs(yin[1])) + (fabs(yin[-2]) + fabs(yin[2])) );

  //cutoff value: if an abs. value of the derivative is smaller than this value, the derivative is considered = 0.
  cutoff_value = norm * epsilon + MINNUMREPRESENT;

  //smooth transition from 0 to 1 depending on the value of minabs_der
  monoindicator = transition_function( minabs_der, 0.0, cutoff_value );

  //if( minabs_der != 0.0 && norm > 1 ) {
  // dualfprintf( fail_file, "monoindicator = %21.15g\n", monoindicator );
  //}

  return( monoindicator );
}

///returns 0 if x < x1
///        1 if x > x2
///        linearly interpolated value if x1 <= x <= x2 
FTYPE transition_function( FTYPE x, FTYPE x1, FTYPE x2 ) 
{
  FTYPE val;

  assert( x1 == x2, "transition_function(): cannot have transition values x1 equal x2\n" );

  val = MAX( 0., MIN(1., (x - x1)/(x2 - x1)) );

  return( val );
}


/// uses yin[-2..2]
/// uses df[-1..2] uses out to yin[-2..2]
/// uses ddf[-1..1] uses yin[-2..2]
int check_for_cusp_new(FTYPE *yin, FTYPE *df, FTYPE *ddf)
{
  FTYPE norm;
  FTYPE sqrtnorm;
  FTYPE f1a,f5a,f1,f2,f3,f3a,f4,f4a,f5,f6,f7,f7a,f8,f8a;
   
  norm=(fabs(yin[-2])+fabs(yin[-1])+fabs(yin[0])+fabs(yin[1])+fabs(yin[2])+SMALL);
  //norm=norm*norm*SQRT_WENO_EPSILON;  //replaced the ERRORNORM with SQRT_WENO_EPSIOLON
  norm=norm*SQRT_WENO_EPSILON;
   
  f1a = fabs(df[0]) - norm;
  f1= ddf[-1] * sign(df[0]) -norm;
  f2=-ddf[0]  * sign(df[0]) -norm;
  f3= -df[2]  * sign(df[0]) -norm;
  f3a = fabs(-df[2] - df[0]) - 0.25 * (fabs(df[2]) + fabs(df[0]));
  f4= -df[1]  * sign(df[0]) -norm;
  f4a = fabs(-df[1] - df[0]) - 0.25 * (fabs(df[1]) + fabs(df[0]));
       
  //f3 = 0.0;  //SASMARKx

#if(0)
  if(crapdebug||1){
    dualfprintf(fail_file,"f1a=%21.15g f1=%21.15g f2=%21.15g f3=%21.15g f3a=%21.15g f4=%21.15g f4a=%21.15g\n",f1a,f1,f2,f3,f3a,f4,f4a);
   
  } 
#endif          

         
  if(f1a > 0 && f1>0 && f2>0 && ((f3>0&&f3a>0) || (f4>0&&f4a>0)) ) return(1);
        
  f5a = fabs(df[1]) - norm;
  f5=-ddf[1]  * sign(df[1]) -norm;
  f6= ddf[0]  * sign(df[1]) -norm;
  f7= -df[-1] * sign(df[1]) -norm;
  f7a = fabs(-df[-1] - df[1]) - 0.25 * (fabs(df[-1]) + fabs(df[1]));
  f8= -df[0]  * sign(df[1]) -norm;
  f8a = fabs(-df[0] - df[1]) - 0.25 * (fabs(df[0]) + fabs(df[1]));
            
  //f7 = 0.0;  //SASMARKx

#if(0)
  if(crapdebug||1){
    dualfprintf(fail_file,"f5a=%21.15g f5=%21.15g f6=%21.15g f7=%21.15g f7a=%21.15g f8=%21.15g f8a=%21.15g\n",f5a,f5,f6,f7,f7a,f8,f8a);
   
  }           
#endif
              
  if(f5a > 0 && f5>0 && f6>0 && ((f7>0&&f7a>0) || (f8>0&&f8a>0)) ) return(1);
  

                
  return(0);
}
               
