// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h


void set_arrays()
{

  // nothing so far

}



void set_multidimen_arrays()
{
  int i, j, k, pl, pliter, l, m;
  int pl2;
  int ii;
  int floor,pf, tscale,dtstage;
  int dissloop;
  FTYPE valueinit;
  int dir,interpi,enodebugi;
  int dimen;





#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif



  ////////////////////////////////////////////////
  //
  // Basic primitive quantity
  //
  ////////////////////////////////////////////////

  GLOBALPOINT(p) = (FTYPE PTRMACP0A1(p,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(p,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    MACP0A1(p,i,j,k,pl) = valueinit;
  }
  


  GLOBALPOINT(gcon) = (FTYPE PTRMETMACP1A2(gcon,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM)) (&(BASEMETMACP1A2(gcon,0,N1BND,N2BND,N3BND,0,0)));
  
  // rest are always located at CENT
  GLOBALPOINT(conn) = (FTYPE PTRMETMACP0A3(conn,N1M,N2M,N3M,NDIM,NDIM,NDIM)) (&(BASEMETMACP0A3(conn,N1BND,N2BND,N3BND,0,0,0)));

}


