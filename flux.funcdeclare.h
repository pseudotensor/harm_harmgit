
/*! \file flux.funcdeclare.h
    \brief Function declarations for flux.c

*/


/// stepping
extern int step_ch_full(int truestep,FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);
extern void get_truetime_fluxdt(int numtimeorders, SFTYPE localdt, FTYPE (*CUf)[NUMDTCUFS], FTYPE (*Cunew)[NUMDTCUFS], SFTYPE *fluxdt, SFTYPE *boundtime, SFTYPE *fluxtime, SFTYPE *tstepparti, SFTYPE *tsteppartf);

extern void set_normal_realisinterp(int *realisinterp);

extern int fluxcalc(int stage,
                    int initialstep, int finalstep,
                    FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                    FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
                    FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                    FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                    FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
                    FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
                    FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], 
                    FTYPE *CUf,
                    FTYPE *CUnew,
                    SFTYPE fluxdt,
                    SFTYPE fluxtime,
                    FTYPE *ndt1,
                    FTYPE *ndt2,
                    FTYPE *ndt3
                    );


extern int fluxcalc_fluxctstag(int stage,
                               int initialstep, int finalstep,
                               FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
                               //          FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                               FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                               FTYPE (*wspeed)[COMPDIM][NUMCS][NSTORE1][NSTORE2][NSTORE3],
                               FTYPE (*prc)[NSTORE2][NSTORE3][NPR2INTERP],
                               FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP],
                               FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP],
                               struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
                               struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
                               FTYPE (*geomcornglobal)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                               int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL],
                               FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                               FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime, struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM]);

extern void mergedc2ea2cmethod_compute(int *Nvec,FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
extern int flux_ct(int stage,
                   int initialstep, int finalstep,
                   FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], int *Nvec, FTYPE *CUf, FTYPE *CUnew, SFTYPE fluxdt, SFTYPE fluxtime);

extern void rescale_calc_full(int dir,FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR]);

extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);

extern void set_plpr(int dir, SFTYPE time, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r); // from user boundary routine

extern void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );
extern void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );

extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int loc, int continuous, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);


extern void store_geomcorn(int corner, int odir1, int odir2,FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

