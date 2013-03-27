

#if(SIMULBCCALC==1)
// GODMARK : not setup for 3D and didn't yet work for 2D

// we introduce a stage-dependent loop (non-interfering stages and parts of stages)
// these are completely general and change shape depending upon absolute specified ranges.  All ranges are as specified as if normal non-MPI
#define MIDDLEI(istop) (istop-(N1-1-SAFESIZE)/2)
// safe left i across j
#define STAGECONDITION0(istart,istop,jstart,jstop) ((j>=jstart+SAFESIZE)&&(j<=jstop-SAFESIZE)&&(i>=istart+SAFESIZE)&&(i<=MIDDLEI(istop)))
// safe right i across j
#define STAGECONDITION1(istart,istop,jstart,jstop) ((j>=jstart+SAFESIZE)&&(j<=jstop-SAFESIZE)&&(i>=MIDDLEI(istop)+1)&&(i<=istop-SAFESIZE))
// left unsafe across j
#define STAGECONDITION20(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=jstop)&&(i>=istart)&&(i<=istart+SAFESIZE-1))
// right unsafe across j
#define STAGECONDITION21(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=jstop)&&(i>=istop+1-SAFESIZE)&&(i<=istop))
// upper j unsafe across i
#define STAGECONDITION22(istart,istop,jstart,jstop) ((j>=jstart)&&(j<=SAFESIZE-1)&&(i>=istart+SAFESIZE)&&(i<=istop-SAFESIZE))
// lower j unsafe across i
#define STAGECONDITION23(istart,istop,jstart,jstop) ((j>=N2-SAFESIZE)&&(j<=jstop)&&(i>=istart+SAFESIZE)&&(i<=istop-SAFESIZE))


#define STAGECONDITION2(istart,istop,jstart,jstop) (STAGECONDITION20(istart,istop,jstart,jstop)||STAGECONDITION21(istart,istop,jstart,jstop)||STAGECONDITION22(istart,istop,jstart,jstop)||STAGECONDITION23(istart,istop,jstart,jstop))

#define STAGECONDITION(istart,istop,jstart,jstop) ((stage==STAGEM1)||(stage==STAGE0)&&(STAGECONDITION0(istart,istop,jstart,jstop))||(stage==STAGE1)&&(STAGECONDITION1(istart,istop,jstart,jstop))||(stage==STAGE2)&&(STAGECONDITION2(istart,istop,jstart,jstop)))


#define COMPZLOOP ZLOOP if(STAGECONDITION(0,N1-1,0,N2-1))
// GODMARK: ZSLOOP kept
#define COMPZSLOOP(istart,istop,jstart,jstop) ZSLOOP(istart,istop,jstart,jstop) if(STAGECONDITION(istart,istop,jstart,jstop))
/*
  #define COMPZLOOP ZLOOP
  #define COMPZSLOOP(istart,istop,jstart,jstop) ZSLOOP(istart,istop,jstart,jstop)
*/
// need another set of condition loops for flux calc since need fluxes a bit farther out than values.
// note that flux doesn't need more safe zones! (not even flux_ct)
// flux_ct requires F1 flux to be computed an extra i+1, and j-1,j,j+1
// flux_ct requires F2 flux to be computed an extra j+1, and i-1,i,i+1
// again, these outer fluxes don't use more primitive variables outside the safe zone

// fluxes stored ok
#define COMPFZLOOP(istart,jstart) COMPZSLOOP(istart,N1,jstart,N2)

// emf requires its own storage since don't want to recalculate emf for each stage, only do necessary calculations

// this is already provided by the static array in flux_ct
// i+1 and j+1
#define COMPEMFZLOOP COMPZSLOOP(0,N1,0,N2)
// used inside flux_ct()
// just stores fluxes, which we have F1 and F2 for storage
// i+1
#define COMPF1CTZLOOP COMPZSLOOP(0,N1,0,N2-1,N3-1)
// j+1
#define COMPF2CTZLOOP COMPZSLOOP(0,N1-1,0,N2,N3-1)
// k+1
#define COMPF2CTZLOOP COMPZSLOOP(0,N1-1,0,N2-1,N3)

// also need a loop for dq-like calculations (i.e. for get_bsqflags()) i-1,i+1  and j-1,j+1  this is safe too since stencil is no bigger than loop size
// must keep new dq1 and dq2 in storage since each stage gets new set but calculations require previous stage values
#define COMPDQZLOOP COMPZSLOOP(-1,N1,-1,N2)

// goes over only primitive quantities for rescale()
// affects pr, but undone, so no special storage needed
//#define COMPPREDQZLOOP COMPZSLOOP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND)
#define COMPPREDQZLOOP FULLLOOP


#elif(SIMULBCCALC==2)
#define STAGESETUPM1(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istart;ie=istop;
#define MIDDLEI(istop) (istop-(N1-1-SAFESIZE)/2)
// safe left i across j
#define STAGESETUP0(istart,istop,jstart,jstop,is,ie,js,je) js=jstart+SAFESIZE;je=jstop-SAFESIZE;is=istart+SAFESIZE;ie=MIDDLEI(istop);
// safe right i across j
#define STAGESETUP1(istart,istop,jstart,jstop,is,ie,js,je) js=jstart+SAFESIZE;je=jstop-SAFESIZE;is=MIDDLEI(istop)+1;ie=istop-SAFESIZE;
// left unsafe across j
#define STAGESETUP20(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istart;ie=istart+SAFESIZE-1;
// right unsafe across j
#define STAGESETUP21(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=jstop;is=istop+1-SAFESIZE;ie=istop;
// upper j unsafe across i
#define STAGESETUP22(istart,istop,jstart,jstop,is,ie,js,je) js=jstart;je=SAFESIZE-1;is=istart+SAFESIZE;ie=istop-SAFESIZE;
// lower j unsafe across i
#define STAGESETUP23(istart,istop,jstart,jstop,is,ie,js,je) js=N2-SAFESIZE;je=jstop;is=istart+SAFESIZE;ie=istop-SAFESIZE;

#define STAGECONDITION(istart,istop,jstart,jstop,is,ie,js,je) if(stage==STAGEM1){ STAGESETUPM1(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE0){ STAGESETUP0(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE1){ STAGESETUP1(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE2){ STAGESETUP20(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE3){ STAGESETUP21(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE4){ STAGESETUP22(istart,istop,jstart,jstop,is,ie,js,je) } else if(stage==STAGE5){ STAGESETUP23(istart,istop,jstart,jstop,is,ie,js,je) }

#define STAGELOOP(is,ie,js,je) for(i=is;i<=ie;i++) for(j=js;j<=je;j++)

#define TYPE2 (0)
#define COMPZSLOOP(istart,istop,jstart,jstop,is,ie,js,je) STAGECONDITION(istart,istop,jstart,jstop,is,ie,js,je) STAGELOOP(is,ie,js,je)
#define COMPZLOOP COMPZSLOOP(0,N1-1,0,N2-1,isc,iec,jsc,jec)
#define COMPFZLOOP(istart,jstart) COMPZSLOOP(istart,N1,jstart,N2,isc,iec,jsc,jec)
#define COMPEMFZLOOP CZSLOOP(0,N1,0,N2,isc,iec,jsc,jec)
#define COMPF1CTZLOOP COMPZSLOOP(0,N1,0,N2-1,N3-1,isc,iec,jsc,jec)
#define COMPF2CTZLOOP COMPZSLOOP(0,N1-1,0,N2,N3-1,isc,iec,jsc,jec)
#define COMPF3CTZLOOP COMPZSLOOP(0,N1-1,0,N2-1,N3,isc,iec,jsc,jec)
#define COMPDQZLOOP COMPZSLOOP(-1,N1,-1,N2,isc,iec,jsc,jec)
#define COMPPREDQZLOOP COMPZSLOOP(-N1BND,N1-1+N1BND,-N2BND,N2-1+N2BND,isc,iec,jsc,jec)


/*
  #define TYPE2 (1)
  #define COMPZLOOP STAGELOOP(isc,iec,jsc,jec)
  // ief1 and jef1 are equal to ief2 and jef2, so ok
  #define COMPFZLOOP(istart,jstart) STAGELOOP(istart,ief1,jstart,jef1)
  #define COMPEMFZLOOP STAGELOOP(ise,iee,jse,jee)
  #define COMPF1CTZLOOP STAGELOOP(isf1ct,ief1ct,jsf1ct,jef1ct)
  #define COMPF2CTZLOOP STAGELOOP(isf2ct,ief2ct,jsf2ct,jef2ct)
  #define COMPDQZLOOP STAGELOOP(isdq,iedq,jsdq,jedq)
  #define COMPPREDQZLOOP STAGELOOP(ispdq,iepdq,jspdq,jepdq)
*/










#else

// done in global.h




#endif

