#ifndef RECONSTRUCTENOGFD_H
#define RECONSTRUCTENOGFD_H

extern int eno_line_a2c( int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE (*shockindicator)[NBIGM],    FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *Pindicator, FTYPE *yin,  FTYPE *yout, struct of_trueijkp *trueijkp);
extern int eno_line_c2a( int whichquantity, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE (*shockindicator)[NBIGM],     FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *Pindicator, FTYPE *yin,  FTYPE *yout, struct of_trueijkp *trueijkp);
extern int eno_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit, FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM],    FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE *Pindicator, FTYPE *yin,  FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);
extern int paraenohybrid_line_c2e( int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array, int whichreduce, int preforder, int pl, int bs, int ps, int pf, int bf, int *minorderit, int *maxorderit, int *shiftit,        FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM], FTYPE (*Pline)[NBIGM],       FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM],        FTYPE *Pindicator, FTYPE *yin, FTYPE *yout_left, FTYPE *yout_right, FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp );

extern void compute_c2a_polycoef_simple_eno(int full_order, FTYPE *yin, FTYPE *youtpolycoef );
extern void c2e_simple_eno(int full_order, int is_interpolate_to_right, FTYPE *yin, FTYPE *pout);
extern void a2c_simple_eno(int order, FTYPE *yin, FTYPE *pout);
extern void c2a_simple_eno(int order, FTYPE *yin, FTYPE *pout);
extern void c2e_simple_weno(int order, int ii, int bs, int be, FTYPE *yin, FTYPE *pleft, FTYPE *pright);

extern FTYPE limit_ac_correction( int order, int pl, int bs, int bf, FTYPE max_frac_difference, FTYPE *yin, FTYPE *yout );





extern void pass_1d_line_weno(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);

extern void pass_1d_line_weno_withweights(int whichquantity, int dir, int do_weight_or_recon, weno_weights_t *stencil_weights_array,  int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);

//extern  void pass_1d_line_weno(int whichquantity, int dir, int do_weight_or_recon, int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NBIGM], FTYPE (*monoindicator)[NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[NBIGM], FTYPE (*yin)[NBIGM], FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM], struct of_trueijkp *trueijkp);

extern void compute_multipl_weno(int MULTIPLTYPE, int whichquantity, int dir, int do_weight_or_recon, weno_weights_t (*stencil_weights_array_allpl)[NBIGM], int recontype, int whichreduce, int preforder, int bs, int ps, int pe, int be, int *minorder, int *maxorder, int *shift,   FTYPE (*shockindicator)[NBIGM], FTYPE *stiffindicator, FTYPE (*Vline)[NBIGM],  FTYPE (*Pline)[NBIGM], FTYPE (*df)[NUMDFS][NBIGM], FTYPE (*dP)[NBIGM], FTYPE (*etai)[NUMTRUEEOMSETS][NBIGM], FTYPE (*monoindicator)[NUMMONOINDICATORS][NBIGM], FTYPE (*yprim)[2][NBIGM], FTYPE (*ystencilvar)[2][NBIGM], FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM], struct of_trueijkp *trueijkp);



#endif
