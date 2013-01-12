#include "decs.h"

void set_arrays_multidimen_rad(void)
{
  int i,j,k;
  int valueinit;
  // KORAL


#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif

#if(EOMRADTYPE!=EOMRADNONE)
  GLOBALPOINT(boostemu) = 
    (FTYPE PTRMETMACP2A2(boostemu,BOOSTGRIDPOS,BOOSTDIRS,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM)) 
    (&(BASEMETMACP2A2(boostemu,-CENT,0,N1BND,N2BND,N3BND,0,0)));
  int ll,mm,nn,oo;
  for(ll=CENT;ll<CENT+BOOSTGRIDPOS;ll++) for(mm=0;mm<BOOSTDIRS;mm++) FULLLOOPP1 DLOOP(nn,oo){
      GLOBALMETMACP2A2(boostemu,ll,mm,i,j,k,nn,oo) = valueinit;
    }
#endif


}
