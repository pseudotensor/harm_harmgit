//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function declarations
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPECIAL USE OF PTRDEF: (and in metric.c):
  //  int interpX_gcov(FTYPE *X, struct of_compgeom (*compgeom)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*gcovgrid)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][SYMMATRIXNDIM], FTYPE (*gcovpertgrid)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][NDIM], FTYPE *gcov, FTYPE *gcovpert);
int interpX_gcov(FTYPE *X, struct of_compgeom PTRDEFMETMACP1A0(compgeom,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3), FTYPE PTRDEFMETMACP1A2(gcovgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM), FTYPE PTRDEFMETMACP1A1(gcovpertgrid,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM), FTYPE *gcov, FTYPE *gcovpert);








extern int main(int argc, char *argv[]);

extern int init(int *argc, char **argv[]);
extern void parainitchecks(void);
extern void myargs(int argc, char *argv[]);



extern int set_dt(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE *dt);




extern void init_3dvpot_fullloopp1(FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern void init_3dvpot(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);


extern void init_3dnpr_fullloop(FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);


extern void copy_3dpftype_special(int is, int ie, int js, int je, int ks, int ke,PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3]);
extern void copy_3dpftype_special_fullloop(PFTYPE (*source)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*destspecial)[NSTORE2][NSTORE3]);


extern void copy_3d_fieldonly_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3d_fieldonly_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_fieldonly(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3d_nofield_nowait(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl_nowait(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl_fullloop_nowait(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);




extern void copy_3d_onepl_fullloop(int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);

extern void copy_3dnpr2interp_2ptrs_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP]);
extern void copy_3dnpr2interp_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR2INTERP]);


extern void copy_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void init_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR]);
extern void init_3dnpr(int is, int ie, int js, int je, int ks, int ke,FTYPE initvalue, FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_nofield(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3d_onepl(int is, int ie, int js, int je, int ks, int ke, int pl, FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3dnpr_fullloop(FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest)[NSTORE2][NSTORE3][NPR]);
extern void copy_3dnpr_2ptrs(int is, int ie, int js, int je, int ks, int ke,FTYPE (*source)[NSTORE2][NSTORE3][NPR],FTYPE (*dest1)[NSTORE2][NSTORE3][NPR],FTYPE (*dest2)[NSTORE2][NSTORE3][NPR]);


// stepping
extern int step_ch_full(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR],FTYPE (*F2)[NSTORE2][NSTORE3][NPR],FTYPE (*F3)[NSTORE2][NSTORE3][NPR],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);
extern void get_truetime_fluxdt(int numtimeorders, SFTYPE localdt, FTYPE (*CUf)[4], FTYPE (*Cunew)[4], SFTYPE *fluxdt, SFTYPE *boundtime, SFTYPE *tstepparti, SFTYPE *tsteppartf);

extern void set_normal_realisinterp(int *realisinterp);

extern int fluxcalc(int stage,
		    FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
		    FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],
		    FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
		    FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
		    FTYPE (*F1)[NSTORE2][NSTORE3][NPR], 
		    FTYPE (*F2)[NSTORE2][NSTORE3][NPR], 
		    FTYPE (*F3)[NSTORE2][NSTORE3][NPR], 
		    FTYPE CUf,
		    FTYPE fluxdt,
		    FTYPE *ndt1,
		    FTYPE *ndt2,
		    FTYPE *ndt3
		    );

extern int fluxcalc_fluxctstag(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],
			       //			       FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
			       FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
			       FTYPE (*wspeed)[NUMCS][NSTORE1][NSTORE2][NSTORE3],
			       FTYPE (*prc)[NSTORE2][NSTORE3][NPR2INTERP],
			       FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP],
			       FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP],
			       struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
			       struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
			       FTYPE (*geomcornglobal)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
			       int *Nvec, FTYPE (*dqvec[NDIM])[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE CUf, struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM]);

extern void mergedc2ea2cmethod_compute(int *Nvec,FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR]);
extern int flux_ct(int stage, FTYPE (*pb)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vconemf)[NSTORE2][NSTORE3][NDIM-1], FTYPE (*dq1)[NSTORE2][NSTORE3][NPR], FTYPE (*dq2)[NSTORE2][NSTORE3][NPR], FTYPE (*dq3)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);

extern void rescale_calc_full(int dir,FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR]);

extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);

extern void set_plpr(int dir, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *p_l, FTYPE *p_r); // from user boundary routine

  extern void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );
  extern void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );

extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);
extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP]);


extern void store_geomcorn(int corner, int odir1, int odir2,FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);


extern void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
extern void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt);
extern int diag_flux_pureflux(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR],FTYPE (*F3)[NSTORE2][NSTORE3][NPR],SFTYPE Dt);
extern int diag_flux_general(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], SFTYPE Dt);

extern int doingmetricsubstep(void);

extern int compute_new_metric_anystep(int whichtime,FTYPE *CUf, FTYPE *Cunew, FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern int compute_new_metric_substep(FTYPE *CUf, FTYPE *Cunew, FTYPE (*pb)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern int compute_new_metric_longsteps(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR],FTYPE (*F2)[NSTORE2][NSTORE3][NPR],FTYPE (*F3)[NSTORE2][NSTORE3][NPR],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);
extern int compute_new_metric_and_prims(int whichtime, SFTYPE MBHpar, SFTYPE apar, SFTYPE QBHpar, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*pl_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pr_ct)[NSTORE1][NSTORE2][NSTORE3][NPR2INTERP],FTYPE (*F1)[NSTORE2][NSTORE3][NPR],FTYPE (*F2)[NSTORE2][NSTORE3][NPR],FTYPE (*F3)[NSTORE2][NSTORE3][NPR],FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);


extern void control_metric_update(void);



extern int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *prbefore, FTYPE *pr, struct of_newtonstats *newtonstats);



extern void copy_tempucum_finalucum(int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR]);
extern void copy_tempucum_finalucum_fieldonly(int *loop, FTYPE (*tempucum)[NSTORE2][NSTORE3][NPR], FTYPE (*ucum)[NSTORE2][NSTORE3][NPR]);

extern void get_inversion_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke);
extern void get_stag_startendindices(int *loop, int dir, int *is,int *ie,int *js,int *je,int *ks,int *ke);
extern void get_flux_startendindices(int *loop, int *is,int *ie,int *js,int *je,int *ks,int *ke);


extern int avg2cen_interp(int *locpl, int *whichpltoavg,  int *ifnotavgthencopy, int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);



extern void set_defaults_performance_checks_prepreinit(void);
extern void set_defaults_performance_checks_preinit(void);
extern void set_file_versionnumbers(void);

extern int timecheck(int whichlocation, SFTYPE comptstart);
extern int gocheck(int whichlocation);
extern int output_steptimedt_info(SFTYPE comptstart);


extern int error_check(int wherefrom);
extern int find_horizon(int fromwhere);

// initialize DUMP stuff
extern int init_dumps(void);
extern void init_dnumcolumns(void);



extern int init_linklists(void);
int setuplinklist(int numcolumns,int which);
extern struct blink * addlink(struct blink * clinkptr);

// ENER file stuff
extern int dump_ener(int doener, int dordump, int call_code);

extern int diag(int call_code, FTYPE time, long localnstep, long localrealnstep);

extern void frdotout(void);

extern void report_systeminfo(FILE * fileout);
extern int IsLittleEndian(void);
extern void *SwapEndian(void* Addr, const int Nb);

extern void makedirs(void);

extern void appendener(FILE* ener_file,SFTYPE (*pdot_tot)[NPR],SFTYPE*fladd_tot,SFTYPE*sourceadd_tot);

extern void divbmaxavg(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE*ptrdivbmax,FTYPE*ptrdivbavg);
extern void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void gettotali(int numvars, int* vars[],int*sizes,int*vars_tot[]);
extern int constotal(int enerregion, SFTYPE *vars_tot);
extern int integrate(int numelements, SFTYPE * var,SFTYPE *var_tot,int type, int enerregion);

extern int counttotal(int enerregion, CTYPE *vars_tot, int num);
extern int integratel(int numelements, CTYPE * var,CTYPE *var_tot,int type, int enerregion);


// DUMP file stuff
extern int isenoughfreespace(unsigned long long need);




// initialize stuff
// specific to init.c's and used in initbase.c and init.c, so leave global
extern int post_init_specific_init(void);
extern int pre_init_specific_init(void);
extern int prepre_init_specific_init(void);
extern int init_consts(void);
extern int init_grid(void);
extern int init_global(void);
extern int init_defcoord(void);

extern int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR],FTYPE (*F1)[NSTORE2][NSTORE3][NPR],FTYPE (*F2)[NSTORE2][NSTORE3][NPR],FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int copy_prim2panalytic(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR]);


extern int init_vpot(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);

extern int update_vpot(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],FTYPE (*ptrfluxvec[NDIM])[NSTORE2][NSTORE3][NPR], FTYPE CUf,FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);

extern int normalize_field_withnorm(FTYPE norm, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);
extern int assign_fieldconservatives_pointvalues(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR]);
extern void setfdivb(FTYPE *divb, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], int i, int j, int k);

extern int copy_vpot2flux(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int evolve_withvpot(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*uconstemp)[NSTORE2][NSTORE3][NPR]);

extern int init_vpot_user(int *whichcoord, int l, int i, int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A);
extern int init_vpot2field_user(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);

extern int field_Bhat_fluxrecon(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR]);


extern int transform_primitive_vB(int whichvel, int whichcoord, int i,int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int transform_primitive_pstag(int whichvel, int whichcoord, int i,int j, int k, FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int init_zero_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR]);

extern int pi2Uavg(int *fieldfrompotential, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*Upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*Uavg)[NSTORE2][NSTORE3][NPR]);

extern void set_default_nprlists(void);

extern int addremovefieldfromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int interptype, int dir, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*current_in)[NSTORE2][NSTORE3][NPR], FTYPE (*current_out)[NSTORE2][NSTORE3][NPR]);
extern int addremovefromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);

extern int addremovefromanynpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *anynprstart, int *anynprend, int *anynprlist, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);

// called in restart.c and initbase.c
extern void set_grid(int whichtime,FTYPE *CUf, FTYPE *Cunew);



extern int assignmetricstorage_new(struct of_compgeom *mygeom, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);
extern int assignmetricstorage_old(int loc, int i, int j, int k, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);
extern int assignmetricstorage_oldlast(int loc, int i, int j, int k, FTYPE **localgcov, FTYPE **localgcon, FTYPE **localgcovpert, FTYPE **localgdet, FTYPE **localgdetvol, FTYPE **localalphalapse, FTYPE **localbetasqoalphasq, FTYPE **beta, FTYPE **localeomfunc);

#if(NEWMETRICSTORAGE)
#define GETLOCALMETRIC(loc,i,j,k)     assignmetricstorage_new(&GLOBALMETMACP1A0(compgeom,loc,i,j,k), &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#define GETLASTLOCALMETRIC(loc,i,j,k) assignmetricstorage_new(&GLOBALMETMACP1A0(compgeomlast,loc,i,j,k), &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#else
// uses globals
#define GETLOCALMETRIC(loc,i,j,k)     assignmetricstorage_old(loc, i, j, k, &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#define GETLASTLOCALMETRIC(loc,i,j,k) assignmetricstorage_oldlast(loc, i, j, k, &localgcov, &localgcon, &localgcovpert, &localgdet, &localgdetvol, &localalphalapse, &localbetasqoalphasq, &localbeta, &localeomfunc)
#endif

//  FTYPE (*localgcov)[NDIM];			\
//  FTYPE (*localgcon)[NDIM];			\

#define LOCALMETRICTEMPVARS \
  FTYPE *localgcov; \
  FTYPE *localgcon;\
  FTYPE *localgcovpert;\
  FTYPE *localgdet,*localgdetvol;\
  FTYPE *localalphalapse;\
  FTYPE *localbetasqoalphasq;\
  FTYPE *localbeta;\
  FTYPE *localeomfunc;



extern int higherorder_set(int whichquantity, int recontype, int*weightsplittype);

extern int get_fluxpldirs(int *Nvec, int dir, int *fluxdir, int* pldir, int *plforflux, FTYPE *signflux);
extern void get_odirs(int dir,int *odir1,int *odir2);
extern int set_location_fluxasemforvpot(int dir, int *numdirs, int *odir1, int *odir2, int *loc);
extern int get_numdirs_fluxasemforvpot(int *numdirs, int *fieldloc);

extern int plstart_set(int whichquantity, int dir, int recontype, int *plstart);


// some physics

extern int sourcephysics(FTYPE *ph, struct of_geom *geom, struct of_state *q, FTYPE *Ugeomfree, FTYPE *dUother, FTYPE (*dUcomp)[NPR]);

extern void postdt(void);
extern int primtoU(int returntype, FTYPE *p, struct of_state *q, struct of_geom *geom,
		   FTYPE *U);

extern int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int ucon_calc_4vel_bothut(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *ucon2, FTYPE *others);

extern int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon, FTYPE *others);
extern int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma, FTYPE *qsq);

#if(RELTYPE==RELEOM)

#if(WHICHVEL==VEL4)
#define ucon_calc ucon_calc_4vel
#define dudp_calc dudp_calc_gen
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#elif(WHICHVEL==VEL3)
#define ucon_calc ucon_calc_3vel
#define dudp_calc dudp_calc_3vel
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#elif(WHICHVEL==VELREL4)
#define ucon_calc ucon_calc_rel4vel
#define dudp_calc dudp_calc_gen
#define compute_1plusud0 compute_1plusud0_rel4vel // uses qsq and gamma from ucon_calc_rel4vel()
#define bsq_calc bsq_calc_rel4vel
#define bsq_calc_fromq bsq_calc_fromq_rel4vel
#elif(RELTYPE==NONRELEOM) // not really right
#define ucon_calc ucon_calc_nonrel
#define dudp_calc dudp_calc_nonrel
#define compute_1plusud0 compute_1plusud0_general
#define bsq_calc bsq_calc_general
#define bsq_calc_fromq bsq_calc_fromq_general
#endif

#endif

extern int ucon_calcother(FTYPE *pr, FTYPE *ucon, FTYPE *others);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,FTYPE *CC, FTYPE *discr);



extern FTYPE ranc(int initialize, int seed);

extern FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 );


// fixup stuff

extern int check_pr(FTYPE *pr, FTYPE *prmodel, FTYPE *ucons, struct of_geom *geom, int modelpos,int finalstep);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
		    FTYPE *ucon);

/* // dudp stuff */

/* extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui); */
/* extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon, */
/* 			FTYPE *b, FTYPE (*dbiduj)[NDIM]); */
/* extern void db2dui_calc(FTYPE (*dbiduj)[NDIM], FTYPE *b, */
/* 			FTYPE *db2dui); */
/* extern void duudud_calc(FTYPE *ucon, FTYPE (*duudud)[NDIM]); */

/* extern void dbsqdui_calc(FTYPE (*dbiduj)[NDIM], FTYPE *b, */
/* 			 FTYPE *dbsqdui); */
/* extern void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi); */
/* extern void duidvj_calc(FTYPE *dgdv,FTYPE (*duidvj)[NDIM]); */
/* extern void dudduu_calc(FTYPE*dutdui, FTYPE (*dudduu)[NDIM]); */
/* extern void dbdiduj_calc(FTYPE (*dbiduj)[NDIM],FTYPE (*dbdiduj)[NDIM]); */
/* extern void ducon_dv3_calc(struct of_state *q,FTYPE (*ducon_dv)[NDIM]); */
extern int sp_stress_calc(FTYPE *pr, FTYPE (*tens_matt)[NDIM],
			  FTYPE (*tens_em)[NDIM], FTYPE *b,
			  FTYPE *ucon);



// log file stuff
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);

// boundary stuff
extern void set_boundloop(int boundvartype, int *inboundloop, int*outboundloop, int*innormalloop, int*outnormalloop, int (*inoutlohi)[NUMUPDOWN][NDIM], int *riin, int *riout, int *rjin, int *rjout, int *rkin, int *rkout, int *dosetbc);
extern int report_bound_loop(void);
extern void set_numbnd(int boundvartype, int *numbnd, int *numnpr);

// below are for particular purposes
extern int bound_allprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_evolveprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_prim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_pstag(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_beforeevolveprim(int boundstage, SFTYPE boundtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

// below can choose boundvartype
extern int bound_anyallprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_anyprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_anypstag(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int bound_uavg(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

// only pflag doesn't have boundvartype
extern int bound_pflag(int boundstage, SFTYPE boundtime, PFTYPE (*primbase)[NSTORE2][NSTORE3][NUMPFLAGS], int finalstep);


extern int bound_flux(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], int finalstep);

// user bounds:
extern int bound_prim_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_pstag_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);
extern int bound_prim_user_after_mpi_dir(int boundstage, SFTYPE boundtime, int whichdir, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_flux_user(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int bound_pflag_user(int boundstage, SFTYPE boundtime, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS]);


extern int inflow_check_4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_3vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_rel4vel(int dir, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep);


// transform stuff
extern int bl2met2metp2v(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int bl2met2metp2v_genloc(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk, int loc);
extern int bl2met2metp2v_gen(int whichvel, int whichcoord, int newwhichvel, int newwhichcoord, FTYPE *pr, int ii, int jj, int kk);

extern int ucov_whichcoord2primecoords(int whichcoord, int ii, int jj, int kk, int loc, FTYPE *ucov);

extern int metp2met2bl(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int pr2ucon(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE*ucon);
extern int coordtrans(int whichcoordin, int whichcoordout, int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void bltoks(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void kstobl(int ii, int jj, int kk, int loc, FTYPE*ucon);
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



// metric stuff
extern int metric_checks(struct of_geom *ptrgeom);
extern void gset_genloc(int getprim, int whichcoord, int i, int j, int k, int loc, struct of_geom *geom);
extern void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *geom);
extern FTYPE gdet_func_metric(int whichcoord, FTYPE *V,FTYPE *gcov);
extern FTYPE gdet_func(int whichcoord, FTYPE *gcov);
extern FTYPE gdet_func_singcheck(int whichcoord, FTYPE *V,FTYPE (*generalmatrixlower)[NDIM]);
extern void metric_sing_check(int whichcoord, FTYPE (*genmatrixlower)[NDIM], int *anglesing, int*centersing, int *truedim);
extern void gdetvol_func(struct of_geom *ptrgeom, FTYPE *gdet, FTYPE *eomfunc, FTYPE *gdetvol);
extern void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *EOMFUNCNAME);

//extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
//extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE *gcov);
//extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE *gcon);
extern void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
		      FTYPE (*lconn)[NDIM][NDIM],FTYPE *conn2);
extern void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol);

extern void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcovinfunc, FTYPE *gcovpertinfunc);
extern void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon);


//extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE *gcov);
//extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon);
extern void matrix_inverse_metric(int whichcoord, FTYPE *gcov, FTYPE *gcon);
extern void matrix_inverse(int whichcoord, FTYPE (*genmatrixlower)[NDIM], FTYPE (*genmatrixupper)[NDIM]);
extern void alphalapse_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE *alphalapse);
extern void betasqoalphasq_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE *betasqoalphasq);
extern void beta_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *gcov, FTYPE *gcon, FTYPE alphalapse, FTYPE *beta);


// coordinate stuff
extern void set_coord_parms(int defcoordlocal);
extern void set_coord_parms_nodeps(int defcoordlocal);
extern void set_coord_parms_deps(int defcoordlocal);
extern void write_coord_parms(int defcoordlocal);
extern void read_coord_parms(int defcoordlocal);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void coord_ijk(int i, int j, int k, int loc, FTYPE *X);
extern void coord_free(int i, int j, int k, int loc, FTYPE *X);

extern void bl_coord(FTYPE *X, FTYPE *V);
extern void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V);
extern void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);

extern void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM]);


extern int setihor(void);
extern FTYPE setRin(int ihor);

extern int is_inside_surface(int dir, int ii, int jj, int kk, int pp);
extern int is_on_surface(int dir, int ii, int jj, int kk, int pp);


// some useful wrappers for user simplicity:
extern FTYPE pressure_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE cs2_compute_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE compute_entropy_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern void get_EOS_parms_simple(int*numparms, int i, int j, int k, int loc, FTYPE *parlist);
extern FTYPE compute_temp_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE u);
extern int get_extrasprocessed_simple(int doall, int i, int j, int k, int loc, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_u_from_entropy_simple(int i, int j, int k, int loc, FTYPE rho0, FTYPE entropy);
extern FTYPE compute_qdot_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdrho0_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);
extern FTYPE dpdu_rho0_u_simple(int i, int j, int k, int loc, FTYPE rho, FTYPE u);

// other wrappers
extern int ufromentropy_calc(struct of_geom *ptrgeom, FTYPE entropy, FTYPE *pr);
extern int entropy_calc(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *entropy);
extern int invertentropyflux_calc(struct of_geom *ptrgeom, FTYPE entropyflux,int dir, struct of_state *q, FTYPE*pr);


// for old inversion methods:
extern FTYPE pressure_rho0_w(FTYPE *EOSextra, FTYPE rho0, FTYPE w);

// eos stuff [should be consistent with defs.general.h ptr's and eos.c assignment to pointers and eos.c,idealgaseos.c,etc. versions of those functions.]
extern int pickeos_eomtype(int whicheos, int whicheom);

extern FTYPE pressure_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_u_from_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE entropy);
extern FTYPE u_rho0_p(FTYPE *EOSextra, FTYPE rho0, FTYPE p);
extern FTYPE dpdu_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE dpdrho0_rho0_u(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE cs2_compute(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_entropy(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdrho(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE compute_dSdu(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern FTYPE pressure_wmrho0(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp(FTYPE *EOSextra, FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_qdot(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern int compute_sources_EOS(FTYPE *EOSextra, FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR]);
extern void compute_allextras(int justnum, FTYPE *EOSextra, FTYPE rho0, FTYPE u,int *numextrasreturn,FTYPE*extras);
extern int get_extrasprocessed(int doall, FTYPE *EOSextra, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_temp(FTYPE *EOSextra, FTYPE rho0, FTYPE u);
extern void compute_EOS_parms(FTYPE (*EOSextra)[NSTORE2][NSTORE3][NUMEOSGLOBALS], FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void store_EOS_parms(int numparms, FTYPE *EOSextra, FTYPE *parlist);
extern void get_EOS_parms(int*numparms, FTYPE *EOSextra, FTYPE *parlist);



// physics stuff
extern int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
extern int set_zamo_ucovuconplus1ud0(struct of_geom *ptrgeom, FTYPE *ucov, FTYPE *ucon, FTYPE *plus1ud0);
extern int set_zamo_ucon(struct of_geom *ptrgeom, FTYPE *ucon);

extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);
extern void bsq_calc_rel4vel_fromq(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *bsq);


extern int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE *gamma, FTYPE *qsq);


extern int dudp_calc_gen(int whichcons, FTYPE *EOSextra, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE **alpha);

extern int dudp_calc_3vel(int whichcons, FTYPE *EOSextra, FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha);


extern int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);

extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
extern void UtoU_evolve2diag(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);



extern int Utoprimgen(int finalstep, int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprimloop(FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*pf)[NSTORE2][NSTORE3][NPR], struct of_newtonstats *newtonstats);
extern int primtoUloop(FTYPE (*pi)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR]);

extern int Utoprim(int entropyeom, FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_1d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_1d_opt(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_2d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_1d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
extern int Utoprim_2d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);
//extern int Utoprim_2d_final_nonrelcompat_inputnorestmass(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);  //wrong function name, corrected by atch, see below
extern int Utoprim_jon_nonrelcompat_inputnorestmass(int eomtype, FTYPE *EOSextra, FTYPE *U, struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE *prim, struct of_newtonstats *newtonstats);
extern int Utoprim_5d2_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr, struct of_newtonstats *newtonstats);

extern int Utoprimdiss(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, PFTYPE *otherfail, struct of_newtonstats *newtonstats);



extern int tetr_func(int inputtype, FTYPE *gcov, FTYPE (*tetr_cov)[NDIM],FTYPE (*tetr_con)[NDIM], FTYPE eigenvalues[]);
extern int tetr_func_frommetric(FTYPE (*dxdxp)[NDIM], FTYPE *gcov, FTYPE (*tetrcov)[NDIM],FTYPE (*tetrcon)[NDIM], FTYPE eigenvalues[]);

//extern void SHOULDNOTREACHHEREEVERBUGYOUHAVE(void);



// presume get_geometry() only necessarily feeds back pointer where geometry is located
// if want hard copy to be created using a memory, should use copy_geometry()
// still required to have pointer point to physical allocated memory in general
#if(NEWMETRICSTORAGE)
// overwrites any prior pointer reference
// better than copying since no point in copying if memory already exists and get_geometry() is always for real computational geometry
// This only works if feed pointer directly into get_geometry.  Can't feed &geom.

#if(MCOORD!=CARTMINKMETRIC)

#define get_geometry(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(compgeom,pp,ii,jj,kk);
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP0A1(gdetgeom,ii,jj,kk,pp);
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,ii,jj,kk);
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,ii,jj,kk);
#define set_igdet(arg)
#define set_igdetsimple(arg)

#else// else if CARTMINKMETRIC

#define get_geometry(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(compgeom,pp,0,0,0);
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP0A1(gdetgeom,0,0,0,pp);
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,0,0,0);
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) ptrgeom=&GLOBALMETMACP1A0(gdetgeomnormal,pp,0,0,0);
#define set_igdet(arg)
#define set_igdetsimple(arg)

#endif // end if CARTMINKMETRIC



#else // else if old metric storage method


#define get_geometry(ii,jj,kk,pp,ptrgeom) get_geometry_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_gdetmix(ii,jj,kk,pp,ptrgeom) get_geometry_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_gdetonly(ii,jj,kk,pp,ptrgeom) get_geometry_gdetonly_old(ii,jj,kk,pp,ptrgeom)
#define get_geometry_geomeonly(ii,jj,kk,pp,ptrgeom) get_geometry_geomeonly_old(ii,jj,kk,pp,ptrgeom)
#define set_igdet(arg) set_igdet_old(arg)
#define set_igdetsimple(arg) set_igdetsimple_old(arg)


#endif





#if(NEWMETRICSTORAGE==0)
extern void get_geometry_old(int i, int j, int k, int loc, struct of_geom *geom);
extern void get_geometry_gdetonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void get_geometry_geomeonly_old(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void set_igdet_old(struct of_geom *geom);
extern void set_igdetsimple_old(struct of_geom *geom);
#else
// only thing required:
extern void set_igdet_old(struct of_geom *geom);
#endif


extern void get_and_copy_geometry(int ii, int jj, int kk, int pp, struct of_geom *ptrgeom);


extern void get_allgeometry(int i, int j, int k, int loc, struct of_allgeom *allgeom, struct of_geom *geom);



extern int get_state(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **q);
extern int get_stateforfluxcalc(int dimen, int isleftright, FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);
extern int get_stateforUdiss(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
extern int get_stateforinterpline(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);
extern int get_stateforglobalwavespeeds(FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);


extern void compute_and_store_fluxstatecent(FTYPE (*pr)[NSTORE2][NSTORE3][NPR]);


extern int primtoflux(int returntype, FTYPE *pa, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *fl);
extern int primtoflux_splitmaem(int returntype, FTYPE *pa, struct of_state *q, int fluxdir, int fundir, struct of_geom *geom, FTYPE *flma, FTYPE *flem);


extern int flux_compute_general(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *ctop);
extern int flux_compute_splitmaem(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *FEM, FTYPE *ctop);
extern int get_global_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *pr,FTYPE *wspeedtemp);
extern int get_global_wavespeeds_full(int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*wspeed)[NUMCS][NSTORE1][NSTORE2][NSTORE3]);
extern int global_vchar(FTYPE (*pointspeed)[NSTORE2][NSTORE3][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE (*wspeed)[2][NSTORE1][NSTORE2][NSTORE3]);


extern int get_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, struct of_state *state_l, struct of_state * state_r, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE *ctop);



extern void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom,
		     struct of_state *q,FTYPE *dU);
extern int source(FTYPE *pa, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUriemann,
		  FTYPE (*Uacomp)[NPR], FTYPE *Ua);

extern FTYPE taper_func(FTYPE R,FTYPE rin) ;
extern FTYPE rhor_calc(int which);
extern FTYPE rmso_calc(int which) ;
extern FTYPE uphi_isco_calc(int which,FTYPE r);

extern int set_atmosphere(int whichcond, int whichvel, struct of_geom *geom, FTYPE *pr);
extern int set_density_floors_default(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int fixup(int stage, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int finalstep);
extern int fixup1zone(FTYPE *pr,FTYPE *ucons, struct of_geom *ptrlgeom, int finalstep);
extern int diag_fixup(FTYPE *pr0, FTYPE *pr, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep,int whocalled);
extern int diag_fixup_U(FTYPE *Ui, FTYPE *Uf, FTYPE *ucons, struct of_geom *ptrgeom, int finalstep,int whocalled);


extern int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);

extern int limit_gamma(FTYPE gammamax, FTYPE*pr, FTYPE *ucons, struct of_geom *geom, int finalstep);

extern int fixup_checksolution(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],int finalstep);
extern int fixup_utoprim(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);
extern int fixup_utoprim_nofixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR], FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], int finalstep);

extern int post_fixup(int stage, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int finalstep);
extern int post_fixup_nofixup(int stageit, SFTYPE boundtime, FTYPE (*pv)[NSTORE2][NSTORE3][NPR],FTYPE (*pbackup)[NSTORE2][NSTORE3][NPR],FTYPE (*ucons)[NSTORE2][NSTORE3][NPR],int finalstep);

extern int pre_fixup(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR]);
extern int get_bsqflags(int stage, FTYPE (*pv)[NSTORE2][NSTORE3][NPR]);


extern int fail(int i, int j, int k, int loc, int fail_type);
extern void setfailresponse(int restartonfail);
extern void setrestart(int*appendold);




extern int vchar(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE (*ecov)[NDIM],
			     FTYPE (*econ)[NDIM]);
extern void transform(FTYPE *vec, FTYPE (*t)[NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE (*t)[NDIM]);

extern void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
extern void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

extern void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
extern void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);

extern int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
		      FTYPE *bcon);

extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
extern void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE (*faraday)[NDIM]);
extern void current_precalc(int which, struct of_geom *geom, struct of_state *q, SFTYPE Dt,FTYPE (*faraday)[3]);
extern void init_varstavg(void);
extern void final_varstavg(FTYPE IDT);
extern int set_varstavg(FTYPE tfrac);
extern void current_calc(FTYPE (*cfaraday)[NSTORE2][NSTORE3][NUMCURRENTSLOTS][3]);
extern int current_doprecalc(int which, FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
extern int average_calc(int doavg);



// interpolation stuff
extern int get_loop(int pointorlinetype, int interporflux, int dir, struct of_loop *loop);
extern int set_interpalltypes_loop_ranges(int pointorlinetype, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);


// line types:
extern void set_interp_loop_gen(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
//extern void set_interp_loop(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
//extern void set_interp_loop_expanded(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);

// point types:
extern int set_interppoint_loop_ranges(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern int set_interppoint_loop_ranges_3Dextended(int interporflux, int *maxis, int *maxie, int *maxjs, int *maxje, int *maxks, int *maxke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_ranges_2D_EMF_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_ranges_geomcorn_formerged(int interporflux, int corner, int odir1, int odir2, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);

extern void set_interppoint_loop(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_expanded(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);


extern int vpot2field_useflux(int *fieldloc,FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int vpot2field_centeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR]);
extern int vpot2field_staggeredfield(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pfield)[NSTORE2][NSTORE3][NPR],FTYPE (*ufield)[NSTORE2][NSTORE3][NPR]);
extern int interpolate_ustag2fieldcent(int stage, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE (*preal)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*ucent)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR]);
extern int vectorpot_fluxreconorfvavg(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int deaverage_fields_fv(FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
extern int field_integrate_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR]);
extern int vectorpot_useflux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
extern int field_Bhat_fluxrecon(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR]);
extern int ucons2upointppoint(SFTYPE boundtime, FTYPE (*pfield)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR],FTYPE (*unew)[NSTORE2][NSTORE3][NPR],FTYPE (*ulast)[NSTORE2][NSTORE3][NPR],FTYPE (*pcent)[NSTORE2][NSTORE3][NPR]);

extern int deaverage_ustag2pstag(FTYPE (*preal)[NSTORE2][NSTORE3][NPR], FTYPE (*ustag)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR]);


// functions for loop stuff
extern void  setup_nprlocalist(int whichprimtype, int *nprlocalstart, int *nprlocalend,int *nprlocallist, int *numprims);

/////////////////////////////////////
//
// NR STUFF
//
/////////////////////////////////////
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
//extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
//		     FTYPE tol);


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

//////////////////////////////
//
// specialty functions
//
//////////////////////////////
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
			FTYPE *Edot);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
			  FTYPE rs, FTYPE urs);
extern void timestep(FTYPE ndtr, FTYPE ndth);
extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);

extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
			  FTYPE r, FTYPE rs, FTYPE urs);
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
			FTYPE *Urs, FTYPE *Edot);
extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);

extern void lower_vec(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void lowerf(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise_vec(FTYPE *v1, struct of_geom *geom, FTYPE *v2);
extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);
extern void set_points(void);
// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;
extern void make_tetr(FTYPE *ucon, FTYPE (*econ)[NDIM]);


extern FTYPE sign_bad(FTYPE a);
extern FTYPE sign_func(FTYPE a);

#ifdef WIN32
// GODMARK: Could refine for a=0
#define sign(a) ((a)>0 ? 1.0 : -1.0)
#else

#if(SUPERLONGDOUBLE)
#define sign(a) (sign_bad(a))
#else
#define sign(a) (copysign(1.0,a))
#endif

#endif



extern FTYPE signavoidzero(FTYPE a);

#ifndef WIN32
extern FTYPE max(FTYPE a, FTYPE b);

extern FTYPE min(FTYPE a, FTYPE b);
#endif

// supplemental trig functions
extern FTYPE mysign(FTYPE x);
extern FTYPE myfabs(FTYPE x);

extern FTYPE mysin(FTYPE th);
extern FTYPE mycos(FTYPE th);

extern FTYPE cot(FTYPE arg);
extern FTYPE csc(FTYPE arg);
extern FTYPE sec(FTYPE arg);
///////////////////////////////////
//
// SUPERLONGDOUBLE declarations
//
///////////////////////////////////

#if(SUPERLONGDOUBLE)
#include "mconf.h"
extern long double ceill ( long double );
extern long double floorl ( long double );
extern long double atan2l ( long double, long double );
extern int signbitl ( long double );
//
extern long double fabsl ( long double );
extern long double sqrtl ( long double );
extern long double cbrtl ( long double );
extern long double expl ( long double );
extern long double logl ( long double );
extern long double tanl ( long double );
extern long double atanl ( long double );
extern long double sinl ( long double );
extern long double asinl ( long double );
extern long double cosl ( long double );
extern long double acosl ( long double );
extern long double powl ( long double, long double );
extern long double tanhl ( long double );
extern long double atanhl ( long double );
extern long double sinhl ( long double );
extern long double asinhl ( long double );
extern long double coshl ( long double );
extern long double acoshl ( long double );
extern long double exp2l ( long double );
extern long double log2l ( long double );
extern long double exp10l ( long double );
extern long double log10l ( long double );
extern long double gammal ( long double );
extern long double lgaml ( long double );
extern long double jnl ( int, long double );
extern long double ynl ( int, long double );
extern long double ndtrl ( long double );
extern long double ndtril ( long double );
extern long double stdtrl ( int, long double );
extern long double stdtril ( int, long double );
extern long double ellpel ( long double );
extern long double ellpkl ( long double );
long double lgammal(long double);
extern int isfinitel ( long double );
#define finite(arg) isfinitel(arg)
//#define isfinite(arg) isfinitel(arg)
#define copysign( a, b ) ( fabsl(a) * sign(b) ) 
extern int merror;
#else



#include <math.h>

#ifdef WIN32
#define finite(arg) _finite(arg)
#define isfinite(arg) _finite(arg)
#endif

#ifndef WIN32
//#if USINGICC==0
#if( !defined(isfinite))
// needed for Sauron
#define isfinite(arg) finite(arg)  //atch -- on mako, in force-free it would complain about multiply-defined __finite() if not include this line
#endif
#endif // end if not defined WIN32

#endif


#ifdef WIN32
#define copysign( a, b ) ( fabs(a) * sign(b) ) 
#endif




#if(!DO_ASSERTS)
#define assert assert_func_empty
#else
#define assert assert_func
#endif


extern int assert_func( int is_bad_val, char *s, ... );
extern int assert_func_empty( int is_bad_val, char *s, ... );




#include "global.funcdeclare.user.h"

