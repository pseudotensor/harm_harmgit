
/*! \file global.bounds.h
  \brief Function declarations (used globally) from bounds.tools.c
 */




extern int bound_x1dn_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_x1up_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

extern int bound_x2dn_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_x2up_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);


extern int bound_x3dn_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_x3up_analytic(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

extern int bound_x1dn_outflow(
                              int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                              int *inboundloop,
                              int *outboundloop,
                              int *innormalloop,
                              int *outnormalloop,
                              int (*inoutlohi)[NUMUPDOWN][NDIM],
                              int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                              int enerregion,
                              int *localenerpos
                              );
extern int bound_x1up_outflow(
                              int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                              int *inboundloop,
                              int *outboundloop,
                              int *innormalloop,
                              int *outnormalloop,
                              int (*inoutlohi)[NUMUPDOWN][NDIM],
                              int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                              int enerregion,
                              int *localenerpos
                              );
extern int bound_x1dn_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x1up_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x2dn_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x2up_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x3dn_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x3up_outflow_simple(
                                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                     int *inboundloop,
                                     int *outboundloop,
                                     int *innormalloop,
                                     int *outnormalloop,
                                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                     int enerregion,
                                     int *localenerpos
                                     );
extern int bound_x1dn_sym(
                          int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                          int *inboundloop,
                          int *outboundloop,
                          int *innormalloop,
                          int *outnormalloop,
                          int (*inoutlohi)[NUMUPDOWN][NDIM],
                          int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                          int enerregion,
                          int *localenerpos
                          );
extern int bound_x2dn_polaraxis_full3d(
                                       int whichcall,
                                       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                       int *inboundloop,
                                       int *outboundloop,
                                       int *innormalloop,
                                       int *outnormalloop,
                                       int (*inoutlohi)[NUMUPDOWN][NDIM],
                                       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                       int enerregion,
                                       int *localenerpos
                                       );
extern int bound_x2dn_polaraxis(
                                int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                int *inboundloop,
                                int *outboundloop,
                                int *innormalloop,
                                int *outnormalloop,
                                int (*inoutlohi)[NUMUPDOWN][NDIM],
                                int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                int enerregion,
                                int *localenerpos
                                );
extern int bound_x2up_polaraxis_full3d(
                                       int whichcall,
                                       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                       int *inboundloop,
                                       int *outboundloop,
                                       int *innormalloop,
                                       int *outnormalloop,
                                       int (*inoutlohi)[NUMUPDOWN][NDIM],
                                       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                       int enerregion,
                                       int *localenerpos
                                       );
extern int bound_x2up_polaraxis(
                                int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                int *inboundloop,
                                int *outboundloop,
                                int *innormalloop,
                                int *outnormalloop,
                                int (*inoutlohi)[NUMUPDOWN][NDIM],
                                int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                int enerregion,
                                int *localenerpos
                                );
extern int bound_x1_periodic(
                             int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                             int *inboundloop,
                             int *outboundloop,
                             int *innormalloop,
                             int *outnormalloop,
                             int (*inoutlohi)[NUMUPDOWN][NDIM],
                             int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                             int enerregion,
                             int *localenerpos
                             );

extern int bound_x2_periodic(
                             int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                             int *inboundloop,
                             int *outboundloop,
                             int *innormalloop,
                             int *outnormalloop,
                             int (*inoutlohi)[NUMUPDOWN][NDIM],
                             int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                             int enerregion,
                             int *localenerpos
                             );

extern int bound_x3_periodic(
                             int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                             int *inboundloop,
                             int *outboundloop,
                             int *innormalloop,
                             int *outnormalloop,
                             int (*inoutlohi)[NUMUPDOWN][NDIM],
                             int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                             int enerregion,
                             int *localenerpos
                             );
extern int bound_x1dn_r0singfixinterior(
                                        int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                                        int *inboundloop,
                                        int *outboundloop,
                                        int *innormalloop,
                                        int *outnormalloop,
                                        int (*inoutlohi)[NUMUPDOWN][NDIM],
                                        int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                                        int enerregion,
                                        int *localenerpos
                                        );
extern int bound_checks1(
                         int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                         int *inboundloop,
                         int *outboundloop,
                         int *innormalloop,
                         int *outnormalloop,
                         int (*inoutlohi)[NUMUPDOWN][NDIM],
                         int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                         int enerregion,
                         int *localenerpos
                         );
extern int extrapfunc(int boundary, int j,int k,
                      int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                      int *inboundloop,
                      int *outboundloop,
                      int *innormalloop,
                      int *outnormalloop,
                      int (*inoutlohi)[NUMUPDOWN][NDIM],
                      int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                      int enerregion,
                      int *localenerpos
                      );
extern int poledeath(int whichx2,
                     int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                     int *inboundloop,
                     int *outboundloop,
                     int *innormalloop,
                     int *outnormalloop,
                     int (*inoutlohi)[NUMUPDOWN][NDIM],
                     int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                     int enerregion,
                     int *localenerpos);
extern int polesmooth(int whichx2,
                      int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                      int *inboundloop,
                      int *outboundloop,
                      int *innormalloop,
                      int *outnormalloop,
                      int (*inoutlohi)[NUMUPDOWN][NDIM],
                      int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
                      int enerregion,
                      int *localenerpos);

extern void user1_adjust_fluxcttoth_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );
extern void user1_adjust_fluxctstag_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern void user1_adjust_fluxcttoth_vpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
extern void user1_adjust_fluxctstag_vpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern void adjust_fluxcttoth_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );
extern void adjust_fluxctstag_emfs(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern void adjust_fluxcttoth_vpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
extern void adjust_fluxctstag_vpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern void check_spc_singularities_user(void);



