
/*! \file initbase.funcdeclare.h
     \brief function declarations for initbase.c functions used globally
*/

// initialize stuff
// specific to init.c's and used in initbase.c and init.c, so leave global
extern int post_init_specific_init(void);
extern int pre_init_specific_init(void);
extern int prepre_init_specific_init(void);
extern int init_consts(void);
extern int init_grid(void);
extern int init_global(void);
extern int init_defcoord(void);

extern int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int copy_prim2panalytic(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR]);


extern int init_vpot(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int get_vpot_fluxctstag_primecoords(SFTYPE time, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *vpot);

extern int init_vpot_justAcov(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int init_vpot_toF(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);



extern int vpot2field(SFTYPE time, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);

extern int update_vpot(int whichmethod, int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*ptrfluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vpot0)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vpotlast)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vpotcum)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
extern int set_emfflux(int whichmethod, int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt,FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpot0)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*vpotlast)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int normalize_field_withnorm(FTYPE norm, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
extern int assign_fieldconservatives_pointvalues(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern void setfdivb(FTYPE *divb, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], int i, int j, int k);

extern int copy_vpot2flux(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int evolve_withvpot(SFTYPE time, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);

extern int init_vpot_user(int *whichcoord, int l, SFTYPE time, int i, int j, int k, int loc, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A);
extern int init_vpot2field_user(SFTYPE time, FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);

extern int field_Bhat_fluxrecon(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR]);


extern int transform_primitive_vB(int whichvel, int whichcoord, int i,int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int assignrough_primitive_pstag(int i,int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern int init_zero_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);

extern int init_postvpot(int i, int j, int k, FTYPE *pr, FTYPE *pstag, FTYPE *ucons);


extern int pi2Uavg(int *fieldfrompotential, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*Upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*Uavg)[NSTORE2][NSTORE3][NPR]);

extern void set_default_nprlists(void);

extern int addremovefieldfromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int interptype, int dir, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*current_in)[NSTORE2][NSTORE3][NPR], FTYPE (*current_out)[NSTORE2][NSTORE3][NPR]);
extern int addremovefromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);

extern int addremovefromanynpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *anynprstart, int *anynprend, int *anynprlist, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);

extern void set_array(void *inbufptr, int num, MPI_Datatype datatype, long double value);
