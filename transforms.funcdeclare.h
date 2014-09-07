
/*! \file transforms.funcdeclare.h
    \brief function declrations for global use of functions in transforms.c
*/

// transform stuff
extern int bl2met2metp2v(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int bl2met2metp2v_genloc(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk, int loc);
extern int bl2met2metp2v_genloc_fieldonly(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk, int loc);
extern int bl2met2metp2v_gen(int whichvel, int whichcoord, int newwhichvel, int newwhichcoord, FTYPE *pr, int ii, int jj, int kk);

extern int ucov_whichcoord2primecoords(int whichcoord, int ii, int jj, int kk, int loc, FTYPE *ucov);

extern int metp2met2bl(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int metp2met2bl_genloc(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk, int loc);
extern int pr2ucon(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE*ucon);

extern int coordtrans(int whichcoordin, int whichcoordout, int ii, int jj, int kk, int loc, FTYPE*ucon);

extern void transVtoVmetric(int whichcoord, int ii, int jj, int kk, int loc, FTYPE ROTANGLE, FTYPE *X, FTYPE *V, FTYPE *Xmetric, FTYPE *Vmetric, FTYPE*gcov, FTYPE *gcovpert);
extern void transVmetrictoV(int whichcoord, int ii, int jj, int kk, int loc, FTYPE ROTANGLE, FTYPE *X, FTYPE *V, FTYPE *Xmetric, FTYPE *Vmetric, FTYPE*gcov, FTYPE *gcovpert);

extern void transVmetrictoV_trans(FTYPE ROTANGLE, FTYPE *Vmetric, FTYPE (*trans)[NDIM]);
extern void transVmetrictoV_ucov(FTYPE ROTANGLE, FTYPE *Vmetric, FTYPE*ucov);

extern void bltoks_trans(int ii, int jj, int kk, int loc, FTYPE (*bl2ks)[NDIM]);
extern void kstobl_trans(int ii, int jj, int kk, int loc, FTYPE (*ks2bl)[NDIM]);

extern void bltoks(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void kstobl(int ii, int jj, int kk, int loc, FTYPE*ucon);

extern void bltoks_ucov(int ii, int jj, int kk, int loc, FTYPE*ucov);

extern void mettometp(int ii, int jj, int kk, FTYPE*ucon);
extern void metptomet(int ii, int jj, int kk, FTYPE*ucon);
extern void mettometp_genloc(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void metptomet_genloc(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void mettometp_simple(FTYPE (*idxdxp)[NDIM], FTYPE*ucon);
extern void metptomet_simple(FTYPE (*dxdxp)[NDIM], FTYPE*ucon);
extern void metptomet_ucov_simple(FTYPE (*idxdxp)[NDIM], FTYPE*ucon);
extern void mettometp_ucov_simple(FTYPE (*dxdxp)[NDIM], FTYPE*ucon);
extern void metptomet_Tud(int ii, int jj, int kk, FTYPE (*Tud)[NDIM]);
extern void metptomet_simple_Tud(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM], FTYPE (*Tud)[NDIM]);
extern void ucon2pr(int whichvel, FTYPE *ucon, struct of_geom *geom, FTYPE *pr);
extern int vcon2pr(int whichvel, FTYPE *vcon, struct of_geom *geom, FTYPE *pr);
