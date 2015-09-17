
/*! \file global.gridsectioning.h
  \brief Definitions and macros for grid sectioning

  // Sasha related enerregion's to grid sectioning, e.g. NUMENERREGIONS
  // See initbase.gridsectioning.c for comments
  // There is no user/problem-dependent code here
 */




#if( DOGRIDSECTIONING )

///corrections to be applied to loops
#define SHIFTX1DN (enerposreg[ACTIVEREGION][X1DN]-0)
#define SHIFTX1UP (enerposreg[ACTIVEREGION][X1UP]-(N1-1))
#define SHIFTX2DN (enerposreg[ACTIVEREGION][X2DN]-0)
#define SHIFTX2UP (enerposreg[ACTIVEREGION][X2UP]-(N2-1))
#define SHIFTX3DN (enerposreg[ACTIVEREGION][X3DN]-0)
#define SHIFTX3UP (enerposreg[ACTIVEREGION][X3UP]-(N3-1))

#else

///no sectioning -- zero corrections
#define SHIFTX1DN (0)
#define SHIFTX1UP (0)
#define SHIFTX2DN (0)
#define SHIFTX2UP (0)
#define SHIFTX3DN (0)
#define SHIFTX3UP (0)

#endif


/// Below WITHINACTIVESECTION() applies to only active computed cells not including their used boundary cells
#if( DOGRIDSECTIONING )
#define WITHINACTIVESECTION(ri,rj,rk) ( ri >= enerposreg[ACTIVEREGION][X1DN] && ri <= enerposreg[ACTIVEREGION][X1UP] \
                                        && rj >= enerposreg[ACTIVEREGION][X2DN] && rj <= enerposreg[ACTIVEREGION][X2UP] \
                                        && rk >= enerposreg[ACTIVEREGION][X3DN] && rk <= enerposreg[ACTIVEREGION][X3UP] )
#else
#define WITHINACTIVESECTION(ri,rj,rk) (ri >=0 && ri<=N1-1 && rj>=0 && rj<=N2-1 && rk>=0 && rk<=N3-1 )  //always within active section since no sectioning (JCM: No, this is only comp cells not BC cells)
#endif


/// Below WITHINACTIVESECTIONEXPAND1() applies to only active computed cells with +-1 boundary cells.  Used, e.g., for setting timestep since often don't want to set timestep in same region flux is computed since just happens that needed p_l and p_r for CT method where fluxes end up not used
#if( DOGRIDSECTIONING )

#define WITHINACTIVESECTIONEXPAND1(ri,rj,rk) ( ri >= enerposreg[ACTIVEREGION][X1DN]-SHIFT1 && ri <= enerposreg[ACTIVEREGION][X1UP]+SHIFT1 \
                                               && rj >= enerposreg[ACTIVEREGION][X2DN]-SHIFT2 && rj <= enerposreg[ACTIVEREGION][X2UP]+SHIFT2 \
                                               && rk >= enerposreg[ACTIVEREGION][X3DN]-SHIFT3 && rk <= enerposreg[ACTIVEREGION][X3UP]+SHIFT3 )

#define LOOPWITHINACTIVESECTIONEXPAND1(ri,rj,rk) GENLOOP(ri,rj,rk,enerposreg[ACTIVEREGION][X1DN]-SHIFT1,enerposreg[ACTIVEREGION][X1UP]+SHIFT1,enerposreg[ACTIVEREGION][X2DN]-SHIFT2,enerposreg[ACTIVEREGION][X2UP]+SHIFT2,enerposreg[ACTIVEREGION][X3DN]-SHIFT3,enerposreg[ACTIVEREGION][X3UP]+SHIFT3)

#define WITHINACTIVESECTIONEXPAND1IS (enerposreg[ACTIVEREGION][X1DN]-SHIFT1)
#define WITHINACTIVESECTIONEXPAND1IE (enerposreg[ACTIVEREGION][X1UP]+SHIFT1)
#define WITHINACTIVESECTIONEXPAND1JS (enerposreg[ACTIVEREGION][X2DN]-SHIFT2)
#define WITHINACTIVESECTIONEXPAND1JE (enerposreg[ACTIVEREGION][X2UP]+SHIFT2)
#define WITHINACTIVESECTIONEXPAND1KS (enerposreg[ACTIVEREGION][X3DN]-SHIFT3)
#define WITHINACTIVESECTIONEXPAND1KE (enerposreg[ACTIVEREGION][X3UP]+SHIFT3)

#else

#define WITHINACTIVESECTIONEXPAND1(ri,rj,rk) (ri >=-SHIFT1 && ri<=N1-1+SHIFT1 && rj>=-SHIFT2 && rj<=N2-1+SHIFT2 && rk>=-SHIFT3 && rk<=N3-1+SHIFT3 )
#define LOOPWITHINACTIVESECTIONEXPAND1(ri,rj,rk) GENLOOP(ri,rj,rk,-SHIFT1,N1-1+SHIFT1,-SHIFT2,N2-1+SHIFT2,-SHIFT3,N3-1+SHIFT3)

#define WITHINACTIVESECTIONEXPAND1IS (-SHIFT1)
#define WITHINACTIVESECTIONEXPAND1IE (N1-1+SHIFT1)
#define WITHINACTIVESECTIONEXPAND1JS (-SHIFT2)
#define WITHINACTIVESECTIONEXPAND1JE (N2-1+SHIFT2)
#define WITHINACTIVESECTIONEXPAND1KS (-SHIFT3)
#define WITHINACTIVESECTIONEXPAND1KE (N3-1+SHIFT3)

#endif


/// Below WITHINACTIVEWITHBNDSECTION() applies to active computed cells including their used boundary cells
#if( DOGRIDSECTIONING )
#define WITHINACTIVEWITHBNDSECTION(ri,rj,rk) ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] \
                                               && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] \
                                               && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] )
#else
#define WITHINACTIVEWITHBNDSECTION(ri,rj,rk) (ri >=-N1BND && ri<=N1-1+N1BND && rj>=-N2BND && rj<=N2-1+N2BND && rk>=-N3BND && rk<=N3-1+N3BND )  //always within active section since no sectioning (JCM: No, this is only comp+bnd cells)
#endif


/// Below WITHINACTIVEWITHBNDSECTION() applies to active computed cells including their used boundary cells
#if( DOGRIDSECTIONING )
#define WITHINACTIVESTAGWITHBNDSECTIONX1(ri,rj,rk) ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= MIN(1+enerposreg[ACTIVEWITHBNDREGION][X1UP],N1-1+N1BND) \
                                                     && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] \
                                                     && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] )
#define WITHINACTIVESTAGWITHBNDSECTIONX2(ri,rj,rk) ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] \
                                                     && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= MIN(1+enerposreg[ACTIVEWITHBNDREGION][X2UP],N2-1+N2BND) \
                                                     && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] )
#define WITHINACTIVESTAGWITHBNDSECTIONX3(ri,rj,rk) ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] \
                                                     && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] \
                                                     && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= MIN(1+enerposreg[ACTIVEWITHBNDREGION][X3UP],N3-1+N3BND) )
#else
/// all same
#define WITHINACTIVESTAGWITHBNDSECTIONX1(ri,rj,rk) (ri >=-N1BND && ri<=N1-1+N1BND && rj>=-N2BND && rj<=N2-1+N2BND && rk>=-N3BND && rk<=N3-1+N3BND )  //always within active section since no sectioning (JCM: No, this is only comp+bnd cells)
#define WITHINACTIVESTAGWITHBNDSECTIONX2(ri,rj,rk) (ri >=-N1BND && ri<=N1-1+N1BND && rj>=-N2BND && rj<=N2-1+N2BND && rk>=-N3BND && rk<=N3-1+N3BND )  //always within active section since no sectioning (JCM: No, this is only comp+bnd cells)
#define WITHINACTIVESTAGWITHBNDSECTIONX3(ri,rj,rk) (ri >=-N1BND && ri<=N1-1+N1BND && rj>=-N2BND && rj<=N2-1+N2BND && rk>=-N3BND && rk<=N3-1+N3BND )  //always within active section since no sectioning (JCM: No, this is only comp+bnd cells)
#endif



/// Below WITHINACTIVEBNDSECTION() applies to only the boundary cells of active computed cells
/// Applies to centered quantities
#if( DOGRIDSECTIONING )
#define WITHINACTIVEBNDSECTION(ri,rj,rk) (                              \
                                          (((ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri < enerposreg[ACTIVEREGION][X1DN]) || (ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && ri > enerposreg[ACTIVEREGION][X1UP])) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP])) \
                                          || (((rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj < enerposreg[ACTIVEREGION][X2DN]) || (rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rj > enerposreg[ACTIVEREGION][X2UP])) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP])) \
                                          || (((rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk < enerposreg[ACTIVEREGION][X3DN]) || (rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] && rk > enerposreg[ACTIVEREGION][X3UP])) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP])) \
                                           )
#else
#define WITHINACTIVEBNDSECTION(ri,rj,rk) (                              \
                                          (((ri >= -N1BND && ri < 0) || (ri <= N1-1+N1BND && ri > N1-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                          || (((rj >= -N2BND && rj < 0) || (rj <= N2-1+N2BND && rj > N2-1)) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                          || (((rk >= -N3BND && rk < 0) || (rk <= N3-1+N3BND && rk > N3-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                           )
#endif

/// Below WITHINACTIVEBNDSECTION() applies to only the boundary cells of active computed cells
/// Applies to staggered quantities (e.g. i=0 is boundary cell for fixed BCs)
#if( DOGRIDSECTIONING )
#define WITHINACTIVESTAGBNDSECTION(ri,rj,rk) (                          \
                                              (((ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEREGION][X1DN]) || (ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && ri > enerposreg[ACTIVEREGION][X1UP])) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP])) \
                                              || (((rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEREGION][X2DN]) || (rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rj > enerposreg[ACTIVEREGION][X2UP])) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP])) \
                                              || (((rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEREGION][X3DN]) || (rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] && rk > enerposreg[ACTIVEREGION][X3UP])) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP])) \
                                               )
#else
#define WITHINACTIVESTAGBNDSECTION(ri,rj,rk) (                          \
                                              (((ri >= -N1BND && ri <= 0) || (ri <= N1-1+N1BND && ri > N1-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              || (((rj >= -N2BND && rj <= 0) || (rj <= N2-1+N2BND && rj > N2-1)) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              || (((rk >= -N3BND && rk <= 0) || (rk <= N3-1+N3BND && rk > N3-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                               )
#endif









/// Below WITHINACTIVEBNDSECTION() applies to only the boundary cells of active computed cells
/// Applies to centered quantities
#if( DOGRIDSECTIONING )
#define WITHINACTIVEBNDSECTIONX1DN(ri,rj,rk) (                          \
                                              ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri < enerposreg[ACTIVEREGION][X1DN] ) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                               )
#define WITHINACTIVEBNDSECTIONX1UP(ri,rj,rk) (                          \
                                              ( ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && ri > enerposreg[ACTIVEREGION][X1UP] ) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                               )

#define WITHINACTIVEBNDSECTIONX2DN(ri,rj,rk) (                          \
                                              ( rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj < enerposreg[ACTIVEREGION][X2DN])  && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                               )
#define WITHINACTIVEBNDSECTIONX2UP(ri,rj,rk) (                          \
                                              ( rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rj > enerposreg[ACTIVEREGION][X2UP]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                               )

#define WITHINACTIVEBNDSECTIONX3DN(ri,rj,rk) (                          \
                                              ( rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk < enerposreg[ACTIVEREGION][X3DN]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP]) \
                                               )
#define WITHINACTIVEBNDSECTIONX3UP(ri,rj,rk) (                          \
                                              ( rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] && rk > enerposreg[ACTIVEREGION][X3UP]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP]) \
  )
#else
#define WITHINACTIVEBNDSECTIONX1DN(ri,rj,rk) (                          \
                                              (((ri >= -N1BND && ri < 0) ) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              )
#define WITHINACTIVEBNDSECTIONX1UP(ri,rj,rk) (                          \
                                              (((ri <= N1-1+N1BND && ri > N1-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              )
#define WITHINACTIVEBNDSECTIONX2DN(ri,rj,rk) (                          \
                                              (((rj >= -N2BND && rj < 0) ) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              )
#define WITHINACTIVEBNDSECTIONX2UP(ri,rj,rk) (                          \
                                              (((rj <= N2-1+N2BND && rj > N2-1)) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                              )
#define WITHINACTIVEBNDSECTIONX3DN(ri,rj,rk) (                          \
                                              (((rk >= -N3BND && rk < 0) ) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                              )
#define WITHINACTIVEBNDSECTIONX3UP(ri,rj,rk) (                          \
                                              (((rk <= N3-1+N3BND && rk > N3-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                              )
#endif



/// Below WITHINACTIVESTAGBNDSECTION() applies to only the boundary cells of active computed cells
/// Applies to staggered quantities (e.g. i=0 is boundary cell for fixed BCs)
#if( DOGRIDSECTIONING )
#define WITHINACTIVESTAGBNDSECTIONX1DN(ri,rj,rk) (                      \
                                                  ( ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEREGION][X1DN] ) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                                   )
#define WITHINACTIVESTAGBNDSECTIONX1UP(ri,rj,rk) (                      \
                                                  ( ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && ri > enerposreg[ACTIVEREGION][X1UP] ) && (rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                                   )

#define WITHINACTIVESTAGBNDSECTIONX2DN(ri,rj,rk) (                      \
                                                  ( rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEREGION][X2DN])  && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                                   )
#define WITHINACTIVESTAGBNDSECTIONX2UP(ri,rj,rk) (                      \
                                                  ( rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP] && rj > enerposreg[ACTIVEREGION][X2UP]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP]) \
                                                   )

#define WITHINACTIVESTAGBNDSECTIONX3DN(ri,rj,rk) (                      \
                                                  ( rk >= enerposreg[ACTIVEWITHBNDREGION][X3DN] && rk <= enerposreg[ACTIVEREGION][X3DN]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP]) \
                                                   )
#define WITHINACTIVESTAGBNDSECTIONX3UP(ri,rj,rk) (                      \
                                                  ( rk <= enerposreg[ACTIVEWITHBNDREGION][X3UP] && rk > enerposreg[ACTIVEREGION][X3UP]) && (ri >= enerposreg[ACTIVEWITHBNDREGION][X1DN] && ri <= enerposreg[ACTIVEWITHBNDREGION][X1UP] && rj >= enerposreg[ACTIVEWITHBNDREGION][X2DN] && rj <= enerposreg[ACTIVEWITHBNDREGION][X2UP]) \
  )
#else
#define WITHINACTIVESTAGBNDSECTIONX1DN(ri,rj,rk) (                      \
                                                  (((ri >= -N1BND && ri <= 0) ) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                                  )
#define WITHINACTIVESTAGBNDSECTIONX1UP(ri,rj,rk) (                      \
                                                  (((ri <= N1-1+N1BND && ri > N1-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                                  )
#define WITHINACTIVESTAGBNDSECTIONX2DN(ri,rj,rk) (                      \
                                                  (((rj >= -N2BND && rj <= 0) ) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                                  )
#define WITHINACTIVESTAGBNDSECTIONX2UP(ri,rj,rk) (                      \
                                                  (((rj <= N2-1+N2BND && rj > N2-1)) && (ri >= -N1BND && ri <= N1-1+N1BND && rk >= -N3BND && rk <= N3-1+N3BND)) \
                                                  )
#define WITHINACTIVESTAGBNDSECTIONX3DN(ri,rj,rk) (                      \
                                                  (((rk >= -N3BND && rk <= 0) ) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                                  )
#define WITHINACTIVESTAGBNDSECTIONX3UP(ri,rj,rk) (                      \
                                                  (((rk <= N3-1+N3BND && rk > N3-1)) && (rj >= -N2BND && rj <= N2-1+N2BND && ri >= -N1BND && ri <= N1-1+N1BND)) \
                                                  )
#endif








extern int setsashawind_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] );
extern int sashawind_set_enerregionupdate(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps);

extern int torus_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] );

extern int jet_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] );
extern int jet_set_myid(void);


/// general functions to be created by users
extern int theproblem_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM]);
extern int theproblem_set_enerregionupdate(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps);
extern int theproblem_set_myid(void);


/// GRIDSECTIONING:
extern int init_gridsectioning(void);
extern int bound_gridsectioning(int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int findandsetactivesection(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
extern int setactivesection(int (*abs)[NDIM], int doprintout);

extern void reset_dothisenerregion(int initialcall);

extern int recompute_fluxpositions(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int setgeneral_enerregion(int (*enerregiondef)[NDIM], int doprintout, int whichregion, int whichbndregion);


extern int setgridsectioning(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );





extern int setflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int setflux_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] );

extern int sethorizonflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int sethorizonflux_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] );

extern int settrueglobalregion(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int settrueglobalregion_set_enerregiondef(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM]);

extern int setjetflux(int initialcall, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );



extern int compute_numcompzones(int (*sectiondef)[NDIM], long long int *localnumcompzones);

