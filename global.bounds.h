extern int bound_x1dn_analytic(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_x1up_analytic(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
extern int bound_x1dn_nssurface(
			 int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
			 int *inboundloop,
			 int *outboundloop,
			 int *innormalloop,
			 int *outnormalloop,
			 int (*inoutlohi)[NUMUPDOWN][NDIM],
			 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
			 int *dosetbc,
			 int enerregion,
			 int *localenerpos
			 );
extern int bound_x1dn_outflow(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       );
extern int bound_x1up_outflow(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
		       );
extern int bound_x1dn_sym(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
		   );
extern int bound_x2dn_polaraxis_full3d(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
				);
extern int bound_x2dn_polaraxis(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
			 );
extern int bound_x2up_polaraxis_full3d(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
				);
extern int bound_x2up_polaraxis(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
			 );
extern int bound_x3_periodic(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
		      );
extern int bound_x1dn_r0singfixinterior(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
				 );
extern int bound_checks1(
		       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
		       int *inboundloop,
		       int *outboundloop,
		       int *innormalloop,
		       int *outnormalloop,
		       int (*inoutlohi)[NUMUPDOWN][NDIM],
		       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		       int enerregion,
		       int *localenerpos
		  );
extern int extrapfunc(int boundary, int j,int k,
	       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
	       int *inboundloop,
	       int *outboundloop,
	       int *innormalloop,
	       int *outnormalloop,
	       int (*inoutlohi)[NUMUPDOWN][NDIM],
	       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
	       int enerregion,
	       int *localenerpos
		      );
extern int poledeath(int whichx2,
	       int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
	       int *inboundloop,
	       int *outboundloop,
	       int *innormalloop,
	       int *outnormalloop,
	       int (*inoutlohi)[NUMUPDOWN][NDIM],
	       int riin, int riout, int rjin, int rjout, int rkin, int rkout, int *dosetbc,
		     int enerregion,
		     int *localenerpos);
extern void adjust_fluxcttoth_emfs(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*emf)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3] );
extern void adjust_fluxctstag_emfs(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR]);
extern void check_spc_singularities_user(void);

extern FTYPE get_omegaf(FTYPE t, FTYPE dt, FTYPE steppart);

