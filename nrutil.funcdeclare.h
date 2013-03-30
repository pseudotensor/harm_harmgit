/////////////////////////////////////
//
// NR STUFF
//
/////////////////////////////////////
extern FTYPE ranc(int initialize, int seed);

extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
//extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
//       FTYPE tol);


/* NR routines from nrutil.h */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
                         long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
                        long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
                         long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);



extern void free_vector(FTYPE *v, long nl, long nh);
extern FTYPE *vector(long nl, long nh);
extern FTYPE **matrix(long nrl, long nrh, long ncl, long nch);
extern unsigned char **cmatrix(int nrl,int nrh,int ncl,int nch);
extern unsigned char ***c3matrix(int a, int b, int c, int d, int e, int f);
extern unsigned char ****c4matrix(int a, int b, int c, int d, int e, int f, int g, int h);
extern unsigned char *****c5matrix(int a, int b, int c, int d, int e, int f, int g, int h, int i, int j);
extern FTYPE **fmatrix(int nrl,int nrh,int ncl,int nch);
extern FTYPE ***f3matrix(int nzl, int nzh, int nrl, int nrh, int ncl, int nch);
extern FTYPE ****f4matrix(int nql, int nqh, int nzl, int nzh, int nrl, int nrh, int ncl, int nch);
extern FTYPE *****f5matrix(int ncoll, int ncolh, int nql, int nqh, int nzl, int nzh, int nrl, int nrh, int ncl, int nch);

extern void free_matrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);

extern void free_cmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch);
extern void free_c3matrix(unsigned char ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_c4matrix(unsigned char ****m, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_c5matrix(unsigned char *****m, long ncoll, long ncolh, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);


extern void free_fmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);
extern void free_f3matrix(FTYPE ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_f4matrix(FTYPE ****m, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_f5matrix(FTYPE *****m, long ncoll, long ncolh, long nql, long nqh, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);


extern void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
extern void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);
extern void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
extern void rotate(FTYPE **r, FTYPE **qt, int n, int i, FTYPE a, FTYPE b);

extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);

extern void lnsrch(int n, FTYPE parms[], FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], FTYPE *f, FTYPE stpmax, int *check, FTYPE (*func)(FTYPE [], FTYPE []));

extern void lubksb(FTYPE **a, int n, int *indx, FTYPE b[]);
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);

extern void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
extern void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
extern void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);


extern FTYPE nrfmin(FTYPE parms[], FTYPE x[]);



extern void newt(int useanalyticjac
                 ,FTYPE parms[]
                 ,FTYPE x[], int n, int *check
                 ,void (*vecfunc)(int n, FTYPE *parms, FTYPE v[], FTYPE f[])
                 ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
                 );
extern void broydn(int useanalyticjac
                   ,FTYPE parms[]
                   ,FTYPE x[], int n, int *check
                   ,void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[])
                   ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
                   );



extern void bcucof(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE d1, FTYPE d2, FTYPE **c);

extern void bcuint(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE x1l,
                   FTYPE x1u, FTYPE x2l, FTYPE x2u, FTYPE x1, FTYPE x2, FTYPE *ansy,
                   FTYPE *ansy1, FTYPE *ansy2);

extern FTYPE rtbis(FTYPE (*func)(FTYPE,FTYPE*), FTYPE *parms, FTYPE x1, FTYPE x2, FTYPE xacc);
