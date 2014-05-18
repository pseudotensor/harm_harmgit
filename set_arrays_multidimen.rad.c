/*! \file set_arrays_multidimen.rad.c
     \brief Sets allocation, pointer shift, and dummy assignment for multi-dimen arrays for RADIATION/Koral code
*/

#include "decs.h"

void set_arrays_multidimen_rad(void)
{
  int i,j,k;
  FTYPE valueinit;
  // KORAL


#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif

#if(EOMRADTYPE!=EOMRADNONE && STORETLAB2ORTHO==1)
  {
    GLOBALPOINT(tlab2ortho) = 
      (FTYPE PTRMETMACP2A2(tlab2ortho,BOOSTGRIDPOS,BOOSTDIRS,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM)) 
      (&(BASEMETMACP2A2(tlab2ortho,-CENT,0,N1BND,N2BND,N3BND,0,0)));
    int ll,mm,nn,oo;
    for(ll=CENT;ll<CENT+BOOSTGRIDPOS;ll++) for(mm=0;mm<BOOSTDIRS;mm++) FULLLOOPP1 DLOOP(nn,oo){
          GLOBALMETMACP2A2(tlab2ortho,ll,mm,i,j,k,nn,oo) = valueinit;
        }
  }
#endif



#if(EOMRADTYPE!=EOMRADNONE)
  {
    GLOBALPOINT(prioritermethod) = 
      (int PTRMACP0A1(prioritermethod,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NUMPRIORITERMETHODINDEX)) 
      (&(BASEMACP0A1(prioritermethod,N1BND,N2BND,N3BND,0)));
    int oo;
    FULLLOOP for(oo=0;oo<NUMPRIORITERMETHODINDEX;oo++){
      GLOBALMACP0A1(prioritermethod,i,j,k,oo) = PRIORITERMETHODNOTSET; //valueinit;
    }
  }
#endif




#if(EOMRADTYPE!=EOMRADNONE)
  {
    GLOBALPOINT(Mradk) = 
      (FTYPE PTRMACP1A1(Mradk,MAXTIMEORDER,N1M,N2M,N3M,NPR)) 
      (&(BASEMACP1A1(Mradk,0,N1BND,N2BND,N3BND,0)));
    int oo,pp;
    FULLLOOP for(oo=0;oo<MAXTIMEORDER;oo++) for(pp=0;pp<NPR;pp++){
        GLOBALMACP1A1(Mradk,oo,i,j,k,pp) = PRIORITERMETHODNOTSET; //valueinit;
      }
  }
#endif


}
