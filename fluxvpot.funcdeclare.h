
/*! \file fluxvpot.funcdeclare.h
  \brief Function declarations for fluxvpot.c
 */

extern int vpot2field_useflux(int *fieldloc,FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int vpot2field_centeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR]);
extern int vpot2field_staggeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR]);
extern int interpolate_ustag2fieldcent(int stage, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE (*preal)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucent)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR]);
extern int vectorpot_fluxreconorfvavg(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int deaverage_fields_fv(FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
extern int field_integrate_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR]);
extern int vectorpot_useflux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int field_Bhat_fluxrecon(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR]);
extern int ucons2upointppoint(SFTYPE boundtime, FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*ulast)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR]);

extern int deaverage_ustag2pstag(FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*ustag)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);


extern int evolve_vpotgeneral(int whichmethod, int stage,
                              int initialstep, int finalstep,
                              FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                              int *Nvec,
                              FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],
                              FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                              FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime,
                              FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]
                              );


extern void adjust_emfs(SFTYPE time, int whichmethod, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );

extern void adjust_vpot(SFTYPE time, int whichmethod,FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

