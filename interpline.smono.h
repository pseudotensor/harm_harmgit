
/*! \file interpline.smono.h
     \brief Function declarations for global use of functions in interpline.smono.c

*/

extern void compute_smonotonicity_line(
                                       int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                       int *minorder, int *maxorder, int *shift,   
                                       FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                       FTYPE (*monoindicator)[NBIGM], FTYPE *yin,  FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);


extern void compute_smonotonicity_line_split(int setindicator, int setyout,
                                             int recontype, int whichreduce, int preforder, int pl, int bs, int ps, int pe, int be, 
                                             int *minorder, int *maxorder, int *shift,   
                                             FTYPE (*shockindicator)[NBIGM], FTYPE (*df)[NBIGM], 
                                             FTYPE (*monoindicator)[NBIGM], FTYPE *yin, FTYPE (*yout)[NBIGM], FTYPE (*youtpolycoef)[NBIGM]);
