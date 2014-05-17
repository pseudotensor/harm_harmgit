
/*! \file diag.funcdeclare.h
    \brief diagnostics subroutine declarations

*/


extern int fail(int i, int j, int k, int loc, int fail_type);
extern void setfailresponse(int restartonfail);
extern void setrestart(int*appendold);







extern void init_varstavg(void);
extern void final_varstavg(FTYPE IDT);
extern int set_varstavg(FTYPE tfrac);
extern int average_calc(int doavg);




extern void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
extern void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt);
extern int diag_flux_pureflux(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL],FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL],SFTYPE Dt);
extern int diag_flux_general(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE Dt);


extern int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *prbefore, FTYPE *pr, struct of_newtonstats *newtonstats);

// ENER file stuff
extern int dump_ener(int doener, int dordump, int call_code);

extern int diag(int call_code, FTYPE time, long localnstep, long localrealnstep);

extern void frdotout(void);

extern void appendener(FILE* ener_file,SFTYPE (*pdot_tot)[NPR],SFTYPE*fladd_tot,SFTYPE*sourceadd_tot);

extern void divbmaxavg(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE*ptrdivbmax,FTYPE*ptrdivbavg);
extern void gettotal(int doall, int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void getalltotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void gettotali(int numvars, int* vars[],int*sizes,int*vars_tot[]);
extern int constotal(int enerregion, SFTYPE *vars_tot);
extern int integrate(int numelements, SFTYPE * var,SFTYPE *var_tot,int type, int enerregion);

extern int counttotal(int enerregion, CTYPE *vars_tot, int num);
extern int integratel(int numelements, CTYPE * var,CTYPE *var_tot,int type, int enerregion);
