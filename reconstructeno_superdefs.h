// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h
// OPENMPNOTE: All these are per global i,j,k, so as with other globals, ok since for each i,j,k separate memory space to be written to when doing loop parallelization

//#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
FTYPE PTRDEFGLOBALMACP0A1(weno_lower_order_fraction,N1M,N2M,N3M,NPR); //atch
//#endif

//#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
FTYPE PTRDEFGLOBALMACP1A0(weno_prim_lower_order_fraction,FILL,N1M,N2M,N3M); //atch
//#endif


//#if(DOENO)
//FTYPE PTRDEFMACP0A1(uenotmp0,N1M,N2M,N3M,NPR); /* for ENO reconstruction of U or dU*/
//FTYPE PTRDEFMACP0A1(uenotmp1,N1M,N2M,N3M,NPR); /* for ENO reconstruction of U or dU*/
//FTYPE PTRDEFMACP0A1(uenotmp2,N1M,N2M,N3M,NPR); /* for ENO reconstruction of U or dU*/
//#endif




#if( WENOMEMORY && WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
FTYPE BASEMACP0A1(weno_lower_order_fraction,N1M,N2M,N3M,NPR); /* space for lower order indicators */  //atch
#endif 

#if( WENOMEMORY && STORE_GAMMA_PRIM_REDUCTION_FRACTION )
FTYPE BASEMACP0A1(weno_prim_lower_order_fraction,NDIM,N1M,N2M,N3M); /* space for lower order fraction of primitive quantities */  //atch
#endif 
