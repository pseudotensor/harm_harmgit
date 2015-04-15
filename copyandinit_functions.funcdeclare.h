
/*! \file copyandinit_functions.funcdeclare.h
    \brief copy/init array header

*/

extern void init_3dvpot_fullloopp1(FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern void init_3dvpot(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);


extern void init_3dnpr_fullloop(FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void init_3dnpr_fullloop_flux(FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR+NSPECIAL]);

extern void init_3dnpr_flux(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR+NSPECIAL]);


extern void copy_3dpftype_special(int is, int ie, int js, int je, int ks, int ke,PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3][NUMFAILPFLAGS]);
extern void copy_3dpftype_special_fullloop(PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3][NUMFAILPFLAGS]);


extern void copy_3d_fieldonly_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3d_fieldonly_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3d_nofield_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl_nowait(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl_fullloop_nowait(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);




extern void copy_3d_onepl_fullloop(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3dnpr2interp_2ptrs_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP]);
extern void copy_3dnpr2interp_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP]);

extern void init_3dnpr_3ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR],FTYPE (*dest3)[NSTORE2][NSTORE3][NPR]);

extern void copy_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void init_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR]);
extern void init_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_nofield(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3dnpr_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR]);

extern void copy_tempucum_finalucum(int whichpl, int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR]);
extern void copy_tempucum_finalucum_fieldonly(int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR]);

extern void copy_3dvpot(int is, int ie, int js, int je, int ks, int ke, FTYPE (*source)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
extern void copy_3dvpot_fullloopp1(FTYPE (*source)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
