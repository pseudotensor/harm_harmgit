
FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 );
int bounds_generate( int i, int j, int k, FTYPE *prim );

#if(DONSEMFS)
FTYPE dfluxns( FTYPE r, FTYPE Omega, FTYPE phi, FTYPE th1, FTYPE th2, FTYPE t, FTYPE dt );
FTYPE vpotns_flux( FTYPE r, FTYPE th1, FTYPE th2, FTYPE ph1, FTYPE ph2);
FTYPE get_ns_alpha();
#endif
