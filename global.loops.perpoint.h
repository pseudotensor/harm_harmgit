
/*! \file global.loops.perpoint.h
    \brief Single-point loop related definitions/macros
    // PER-POINT LOOPS

    
*/



// check for existence in bad form using:
// grep "PLOOP" *.c | grep --invert-match "PLOOP("

// After conversion to PLOOP(pliter,pl) form, regexp to help replace:
// Too unsafe to include () or "\ " in []'s below -- catch afterwards rather than mis-replace
// Note, can't look for " PLOOP(..." since may be tab instead
// 1) PLOOP(\([-_\>a-zA-Z0-9+]+\)) -> PLOOP(pliter,\1)
//    Start with longest name and work to smaller:
//    PFLUXBOUNDLOOP [none apparently]
//    PNOTINTERPLOOP [1]
//    PINTERPLOOP [some]
//    PBOUNDLOOP [several]
//    PDUMPLOOP [several]
//    PINVERTLOOP [none]
//    NUMPRIMLOOP [several]
//    PLOOP
//
// 2) Rename: NUMPRIMLOOPGEN -> NUMPRIMLOOP
//            PLOOPINTERP -> PINTERPLOOP
//            PLOOPNOTINTERP -> PNOTINTERPLOOP
//
// 3) Provide the additional variable:
//
//    regexp:
//             int pl -> int pl,pliter;
//             int dir,pl,sc -> int dir,pl,pliter,sc
//             int i,j,k,pl; -> int i,j,k,pl,pliter;
//             int dir,pl; -> int dir,pl,pliter;
//             int i, j, k, pl, l -> int i, j, k, pl, pliter, l
//             int pl,i,j,k; -> int pl,pliter,i,j,k;
//             int i, j, k, pl; -> int i,j,k,pl,pliter;
//             int i,j,k,l,pl; -> int i,j,k,l,pl,pliter;
//             int i, pl; -> int i,pl,pliter;
//             int i,j,k,pl,l; -> int i,j,k,pl,pliter,l;
//             int l,pl,dir; -> int l,pl,pliter,dir;
//
/// PLOOP controls looping over conserved or primitive -type quantities except during interpolation or other listed below
#if(1) // for now always control interpolated quantities (used to be only for SLPITNPR)

#define PLOOP(pliter,pl) for(pliter=nprstart,pl=nprlist[pliter];pliter<=nprend;pliter++,pl=nprlist[pliter])
//SUPERGODMARK
#define PRIMLOOP(pliter,pl) PLOOP(pliter,pl)
#define PINTERPLOOP(pliter,pl) for(pliter=npr2interpstart,pl=npr2interplist[pliter];pliter<=npr2interpend;pliter++,pl=npr2interplist[pliter])
#define PNOTINTERPLOOP(pliter,pl) for(pliter=npr2notinterpstart,pl=npr2notinterplist[pliter];pliter<=npr2notinterpend;pliter++,pl=npr2notinterplist[pliter])

/// loop over all bounding Primitive variables
#define PBOUNDLOOP(pliter,pl) for(pliter=nprboundstart,pl=nprboundlist[pliter];pliter<=nprboundend;pliter++,pl=nprboundlist[pliter])

/// loop over all bounding flux variables
#define PFLUXBOUNDLOOP(pliter,pl) for(pliter=nprfluxboundstart,pl=nprfluxboundlist[pliter];pliter<=nprfluxboundend;pliter++,pl=nprfluxboundlist[pliter])

/// loop over all dumped Primitive variables
#define PDUMPLOOP(pliter,pd) for(pliter=nprdumpstart,pd=nprdumplist[pliter];pliter<=nprdumpend;pliter++,pd=nprdumplist[pliter])

/// loop over all inversion Primitive variables -- only over 5 quantities (not field)
#define PINVERTLOOP(pliter,pi) for(pliter=nprinvertstart,pl=nprinvertlist[pliter];pliter<=nprinvertend;pliter++,pl=nprinvertlist[pliter])


/// to be used locally:
#define NUMPRIMLOOP(pliter,pl) for(pliter=nprlocalstart,pl=nprlocallist[pliter];pliter<=nprlocalend;pliter++,pl=nprlocallist[pliter])





#else


///  loop over all Primitive variables 
#define PLOOP(pliter,pl) for(pl=0;pl<NPR;pl++) // original
///  loop over all center to edge variables 
#define PINTERPLOOP(pliter,pl) for(pl=0;pl<NPR2INTERP;pl++)
#define PNOTINTERPLOOP(pliter,pl) for(pl=0;pl<0;pl++) // do nothing
///  loop over all bounding Primitive variables 
#define PBOUNDLOOP(pliter,pb) for(pb=0;pb<NPRBOUND;pb++)

///  loop over all dumped Primitive variables 
#define PDUMPLOOP(pliter,pd) for(pd=0;pd<NPRDUMP;pd++)
/// loop over all inversion Primitive variables -- only over 5 quantities (not field)
#define PINVERTLOOP(pliter,pi) for(pi=0;pi<NPRINVERT;pi++)

#define NUMPRIMLOOP(pliter,pl) for(pl=0;pl<numprims;pl++)



#endif




/// always goes over all conserved
#define PALLLOOP(pl) for(pl=0;pl<NPR;pl++)
#define PALLLOOPSPECIAL(pl,special) for(pl=0;pl<NPR+special;pl++)

#define PLOOPSPECIALONLY(pl,special) for(pl=NPR;pl<NPR+special;pl++)


/// always goes over all conserved
#define PALLREALLOOP(pl) for(pl=0;pl<NPRREALSET;pl++)


/// goes over all/any npr lists for copying one list to another
#define PMAXNPRLOOP(pl) for(pl=0;pl<MAXNPR;pl++)

/// always goes over all standard invertable quantities for inversion to operate normally
#define PLOOPALLINVERT(pl) for(pl=0;pl<=B3;pl++)

/// always goes over all interpolatable primitives
#define PLOOPALLINTERP(pl) for(pl=0;pl<NPR2INTERP;pl++)


/// always goes over all primitives
#define PDIAGLOOP(pl) PALLLOOP(pl)





#define PLOOPNOB1(pl) for(pl=0;pl<B1;pl++)
#define PLOOPBONLY(pl) for(pl=B1;pl<=B3;pl++)
#define PLOOPNOB2(pl) for(pl=B3+1;pl<NPR;pl++)
#define PLOOPNOB2SPECIAL(pl,special) for(pl=B3+1;pl<NPR+special;pl++)


#define PLOOPRADONLY(pl) for(pl=PRAD0;pl<=PRAD3;pl++)



///  loop over all Dimensions; second rank loop 
#define DLOOP(j,k) for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
///  loop over all Dimensions; first rank loop 
#define DLOOPA(j) for(j=0;j<NDIM;j++)
///  loop over all Space dimensions; second rank loop 
#define SLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
///  loop over all Space dimensions; first rank loop 
#define SLOOPA(j) for(j=1;j<NDIM;j++)
///  loop over all for j and Space for k; second rank loop 
#define DSLOOP(j,k) for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
///  space-space 
#define SSLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
///  loop over all for k and Space for j; second rank loop 
#define SDLOOP(j,k) for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

/// loop over directions
#define DIRLOOP(dir) for(dir=0;dir<COMPDIM*2;dir++)

#define DIRSIGNLOOP(dirsign) for(dirsign=-1;dirsign<=1;dirsign+=2)


#define DIMENLOOP(dir) for(dir=1;dir<=COMPDIM;dir++)


/// loop over grid positions (used for file writing since not normally done)
#define GRIDLOOP(gridpos) for(gridpos=0;gridpos<NPG;gridpos++)



//////////////////////////////////////
//////////////////////////////////////
///
/// PER-POINT LOOPS (MORE RELATED TO NOT BEING PRIMITIVE OR CONSERVED)
///
//////////////////////////////////////
//////////////////////////////////////


/// loop over fail flag in boundary code
#define FBOUNDLOOP(ff) for(ff=0;ff<NUMPFLAGSBOUND;ff++)

/// loop over jet regions
#define JETLOOP(jetio) for(jetio=0;jetio<NUMJETS;jetio++)

/// loop over ener/flux regions
#define ENERREGIONLOOP(enerregion) for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++) if(dothisenerreg[enerregion])

/// loop over ALL ener/flux regions
#define ENERREGIONALLLOOP(enerregion) for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++)

/// loop over fair/floor types
#define FLOORLOOP(floor) for(floor=0;floor<NUMFAILFLOORFLAGS;floor++)

/// loop over debug time scales
#define TSCALELOOP(tscale) for(tscale=0;tscale<NUMTSCALES;tscale++)


/// loop over ALL sources (including geometry)
#define SCLOOP(sc) for(sc=0;sc<NUMSOURCES;sc++)

/// loop over all sources EXCEPT geometry (assumes GEOMSOURCE==0 or at least nothing before GEOMSOURCE matters)
#define SCPHYSICSLOOP(sc) for(sc=GEOMSOURCE+1;sc<NUMSOURCES;sc++)

/// loop over fluxterms
#define FLLOOP(fl) for(fl=0;fl<NUMFLUXTERMS;fl++)


/// loop over pflag flags
#define PFLAGLOOP(pf) for(pf=0;pf<NUMPFLAGS;pf++)

/// loop over fail pflag flags
#define FAILPFLAGLOOP(pf) for(pf=0;pf<NUMFAILPFLAGS;pf++)

/// for USEMPI&&USEROMIO==1
#define ROMIOCOLLOOP(romiocoliter) for(romiocoliter=0;romiocoliter<romiocloopend;romiocoliter++)

#define BUFFERINIT nextbuf=0
#define COLINIT nextcol=0

/// for mpicombie==0
#define COLLOOP(coliter) for(coliter=0;coliter<numfiles;coliter++)


#define DTSTAGELOOP(dtstage) for(dtstage=0;dtstage<MAXITERDTSTAGES;dtstage++)

#define INTERPENOTYPELOOP(interpi) for(interpi=0;interpi<NUMENOINTERPTYPES;interpi++)

#define ENODEBUGLOOP(enodebugi) for(enodebugi=0;enodebugi<NUMENODEBUGS;enodebugi++)

/// simple loop over 0 and 1 for failfloorcount[]
#define FINALSTEPLOOP(indexfinalstep) for(indexfinalstep=0;indexfinalstep<=1;indexfinalstep++)

#define FAILFLOORLOOP(indexfinalstep,tscale,floor) FINALSTEPLOOP(indexfinalstep) TSCALELOOP(tscale) FLOORLOOP(floor)

