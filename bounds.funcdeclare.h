// boundary stuff
extern void set_boundloop(int boundvartype, int *inboundloop, int*outboundloop, int*innormalloop, int*outnormalloop, int (*inoutlohi)[NUMUPDOWN][NDIM], int *riin, int *riout, int *rjin, int *rjout, int *rkin, int *rkout, int *dosetbc);
extern int report_bound_loop(void);
extern void set_numbnd(int boundvartype, int *numbnd, int *numnpr);

// below are for particular purposes
extern int bound_allprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_evolveprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_prim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_pstag(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_beforeevolveprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);

// below can choose boundvartype
extern int bound_anyallprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_anyprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_anypstag(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);
extern int bound_uavg(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);

// only pflag doesn't have boundvartype
extern int bound_pflag(int boundstage, SFTYPE boundtime, PFTYPE (*primbase)[NSTORE2][NSTORE3][NUMPFLAGS], int finalstep, int doboundmpi);


extern int bound_flux(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], int finalstep, int doboundmpi);

// user bounds:
extern int bound_prim_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_pstag_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int bound_prim_user_after_mpi_dir(int boundstage, SFTYPE boundtime, int whichdir, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_flux_user(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int bound_pflag_user(int boundstage, SFTYPE boundtime, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS]);

