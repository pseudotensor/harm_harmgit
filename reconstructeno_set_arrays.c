// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h

void reconstructeno_set_arrays(void)
{
  int i,j,k,pliter,pl;

  ///////////////////////////
  //
  // COUNTERS and failure checks
  //
  ////////////////////////////

#if( WENOMEMORY && WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  GLOBALPOINT(weno_lower_order_fraction) = (FTYPE PTRMACP0A1(weno_lower_order_fraction,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(weno_lower_order_fraction,N1BND,N2BND,N3BND,0)));  //atch
  
  //init the array with zeroes
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(weno_lower_order_fraction,i,j,k,pl) = 0.0;
#endif
  
#if( WENOMEMORY && STORE_GAMMA_PRIM_REDUCTION_FRACTION )
  GLOBALPOINT(weno_prim_lower_order_fraction) = (FTYPE PTRMACP1A0(weno_prim_lower_order_fraction,FILL,N1M,N2M,N3M)) (&(BASEMACP1A0(weno_prim_lower_order_fraction,0,N1BND,N2BND,N3BND)));  //atch
  
  //init the array with zeroes
  FULLLOOP  DIMENLOOP(dimen) GLOBALMACP1A0(weno_prim_lower_order_fraction,dimen,i,j,k) = 0.0;
#endif





  ////////////////////////////////////////////////
  //
  // DEBUG STUFF USUALLY OFF
  //
  ////////////////////////////////////////////////

#if(DOENODEBUG)

  GLOBALPOINT(enodebugarray) = (CTYPE PTRMACP0A4(enodebugarray,N1M,N2M,N3M,3-1,NUMENOINTERPTYPES,NPR-4,NUMENODEBUGS)) (&(BASEMACP0A4(enodebugarray,N1BND,N2BND,N3BND,0,0,0,0))); //SASMARKK (make dir count from 1 to 2 by changing 0 to -1)  //atch debug


  FULLLOOP DIMENLOOP(dir) INTERPENOTYPELOOP(interpi) PLOOP(pliter,pl) ENODEBUGLOOP(enodebugi){
    if(dir<=2 && pl<=U2){
      GLOBALMACP0A4(enodebugarray,i,j,k,dir-1,interpi,pl,enodebugi)=0;
    }
  }
#endif
}



