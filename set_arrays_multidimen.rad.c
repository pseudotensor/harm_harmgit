void set_arrays_multidimen_rad(void)
{
  // KORAL


#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif


  GLOBALPOINT(pglobal) = (FTYPE PTRMACP0A1(pglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(pglobal,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(pglobal,i,j,k,pl) = valueinit;
  }


}
