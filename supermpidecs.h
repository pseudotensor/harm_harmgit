// large memory items for mpidecs.h
#if(USEMPI)
extern FTYPE workbca[NUMPACKUNPACK][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
extern PFTYPE workbc_inta[NUMPACKUNPACK][COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
#endif
extern FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];
extern PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM];
