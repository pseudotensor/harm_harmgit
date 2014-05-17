
/*! \file bounds.funcdeclare.h
  \brief Boundary conditions declarations

*/


// boundary stuff
extern void set_boundloop(int boundvartype, int *inboundloop, int*outboundloop, int*innormalloop, int*outnormalloop, int (*inoutlohi)[NUMUPDOWN][NDIM], int *riin, int *riout, int *rjin, int *rjout, int *rkin, int *rkout, int *dosetbc);
extern int report_bound_loop(void);
extern void set_numbnd(int boundvartype, int *numbnd, int *numnpr);

// below are for particular purposes
extern int bound_allprim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_evolveprim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_prim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_pstag(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_beforeevolveprim(int boundstage, int finalstep, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);

// below can choose boundvartype
extern int bound_anyallprim(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_anyprim(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_anypstag(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);
extern int bound_uavg(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int doboundmpi);

// only pflag doesn't have boundvartype
extern int bound_pflag(int boundstage, int finalstep, SFTYPE boundtime, PFTYPE (*primbase)[NSTORE2][NSTORE3][NUMPFLAGS], int doboundmpi);


extern int bound_flux(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], int doboundmpi);

extern int bound_vpot(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int doboundmpi, int doboundnonmpi);

// user bounds:
extern int bound_prim_user_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_pstag_user_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int bound_prim_user_after_mpi_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_flux_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int bound_pflag_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS]);
extern int bound_vpot_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int debugspecial3dspc(int which, int whichdir, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);

