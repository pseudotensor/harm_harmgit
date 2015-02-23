// __WORKINGONIT__

/*! \file set_arrays_multidimen.c
     \brief Sets allocation, pointer shift, and dummy assignment for multi-dimen arrays

// N1M,N2M,N3M here are used correctly w.r.t. global.storage.h
// See also bounds.ns.c [remove simplename=aphifromoutflow or deal with it]
// See also kazfulleos.c [shift consistently there] simplename=EOSextraglobal
// See liaison.c : [shift consistently there]
// Otherwise, all shifts of multi-D arrays should be here [true right now 05/13/09]
*/


#include "decs.h"


/// all arrays that are multi-dimensional that required special reordering considerations when changing ORDERSTORAGE
/// OPENMPOPTMARK: Don't optimize since many different types of arrays
void set_arrays_multidimen()
{
  int i, j, k, pl, pliter, l, m;
  int pl2;
  int ii;
  int jj;
  int qq;
  int indexfinalstep,floor,pf, tscale,dtstage;
  FTYPE valueinit;
  int dir,interpi,enodebugi;
  int dimen;
  int isleftright;
  struct of_state *myfluxstatetempptr;
  int firsttimeinloop;
  extern void set_arrays_multidimen_rad(void);



#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif


  ////////////////////////////////////////////////
  //
  // Basic primitive quantity
  //
  ////////////////////////////////////////////////

  GLOBALPOINT(pglobal) = (FTYPE PTRMACP0A1(pglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(pglobal,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(pglobal,i,j,k,pl) = valueinit;
  }

#if(ANALYTICMEMORY)
  GLOBALPOINT(panalytic) = (FTYPE PTRMACP0A1(panalytic,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(panalytic,N1BND,N2BND,N3BND,0)));
#if(FIELDSTAGMEM)
  GLOBALPOINT(pstaganalytic) = (FTYPE PTRMACP0A1(pstaganalytic,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(pstaganalytic,N1BND,N2BND,N3BND,0)));
#endif

  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(panalytic,i,j,k,pl) = valueinit;
#if(FIELDSTAGMEM)
    GLOBALMACP0A1(pstaganalytic,i,j,k,pl) = valueinit;
#endif
  }
#endif
  

#if(NUMPOTHER>0)
  GLOBALPOINT(pother) = (FTYPE PTRMACP1A0(pother,FILL,N1M,N2M,N3M)) (&(BASEMACP1A0(pother,0,N1BND,N2BND,N3BND)));
  FULLLOOP for(pl=0;pl<NUMPOTHER;pl++) GLOBALMACP1A0(pother,pl,i,j,k) = -1;
#endif
  
  
  // used in fixup.c, higherorder_pointavg.c and kazfulleos.c for compute_Hglobal()  
  GLOBALPOINT(ptemparray) = (FTYPE PTRMACP0A1(ptemparray,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(ptemparray,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(ptemparray,i,j,k,pl) = valueinit;

  // used in advance.c to carefully assign conserved quantity per pl
  GLOBALPOINT(utemparray) = (FTYPE PTRMACP0A1(utemparray,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(utemparray,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(utemparray,i,j,k,pl) = valueinit;


#if(DOEVOLVEMETRIC)
  GLOBALPOINT(ucumformetric) = (FTYPE PTRMACP0A1(ucumformetric,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(ucumformetric,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(ucumformetric,i,j,k,pl) = valueinit;
#endif



#if(STOREFLUXSTATE)
  GLOBALPOINT(fluxstate) = (struct of_state PTRMACP1A1(fluxstate,FILL,N1M,N2M,N3M,NUMLEFTRIGHT)) (&(BASEMACP1A1(fluxstate,-1,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(fluxstatecent) = (struct of_state PTRMACP0A0(fluxstatecent,N1M,N2M,N3M)) (&(BASEMACP0A0(fluxstatecent,N1BND,N2BND,N3BND)));

  for(isleftright=0;isleftright<NUMLEFTRIGHT+1;isleftright++){

    firsttimeinloop=1;
    DIMENLOOP(dimen){
      
      if(isleftright==NUMLEFTRIGHT){ // last entry (related to +1 in above isleftright loop) treated as if for centered fluxstatecent
        if(firsttimeinloop==0) break; // break if dimen iterated since only have fluxstatecent[] that doesn't depend upon dimen
        else firsttimeinloop=0;
      }
      else{
        // then correct to go over all dimen for fluxstate[dimen]
      }

      // loop over i,j,k
      FULLLOOP{

        if(isleftright==NUMLEFTRIGHT){ // last entry treated as if for centered fluxstatecent
          // just way to select centered case so below assignments are in a single code
          myfluxstatetempptr = &GLOBALMACP0A0(fluxstatecent,i,j,k);
        }
        else{
          myfluxstatetempptr = &GLOBALMACP1A1(fluxstate,dimen,i,j,k,isleftright);
        }


#if(MERGEDC2EA2CMETHOD)
        PALLLOOP(pl){
          myfluxstatetempptr->prim[pl]=valueinit;
          myfluxstatetempptr->EOMFUNCMAC(pl)=valueinit;
        }
        DLOOPA(jj){
          myfluxstatetempptr->Blower[jj]=valueinit;
          myfluxstatetempptr->vcon[jj]=valueinit;
          myfluxstatetempptr->gdetBcon[jj]=valueinit;
        }
        myfluxstatetempptr->gdet=valueinit;
        myfluxstatetempptr->overut=valueinit;
#endif
        // always done
        myfluxstatetempptr->pressure=valueinit;
        myfluxstatetempptr->entropy=valueinit;
        myfluxstatetempptr->ifremoverestplus1ud0elseud0=valueinit;
        DLOOPA(jj){
          myfluxstatetempptr->ucon[jj]=valueinit;
          myfluxstatetempptr->ucov[jj]=valueinit;
          myfluxstatetempptr->bcon[jj]=valueinit;
          myfluxstatetempptr->bcov[jj]=valueinit;      
        }

      }
    }
  }

#endif


  // below used in fluct.c and as vector potential storage
  GLOBALPOINT(emf) = (FTYPE PTRMACP1A0(emf,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMACP1A0(emf,0,N1BND,N2BND,N3BND)));// inner shift still same
  for(l=0;l<NDIM;l++) FULLLOOP GLOBALMACP1A0(emf,l,i,j,k) = valueinit;


#if(FIELDTOTHMEM)
  // below only used in fluxct.c
  GLOBALPOINT(vconemf) = (FTYPE PTRMACP0A1(vconemf,N1M,N2M,N3M,NDIM-1)) (&(BASEMACP0A1(vconemf,N1BND,N2BND,N3BND,-U1)));
  FULLLOOP for(l=U1;l<U1+COMPDIM;l++) GLOBALMACP0A1(vconemf,i,j,k,l) = valueinit;
#endif

#if(MODIFYEMFORVPOT==MODIFYVPOT || TRACKVPOT>0 || EVOLVEWITHVPOT>0)
  GLOBALPOINT(vpotarrayglobal) = (FTYPE PTRMACP1A0(vpotarrayglobal,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMACP1A0(vpotarrayglobal,0,N1BND,N2BND,N3BND)));// inner shift still same
  for(l=0;l<NUMVPOT;l++) FULLLOOP{
      GLOBALMACP1A0(vpotarrayglobal,l,i,j,k) = valueinit;
    }

#if(ANALYTICMEMORY)
  GLOBALPOINT(vpotanalytic) = (FTYPE PTRMACP1A0(vpotanalytic,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMACP1A0(vpotanalytic,0,N1BND,N2BND,N3BND)));// inner shift still same

  for(l=0;l<NDIM;l++) FULLLOOP{
      GLOBALMACP1A0(vpotanalytic,l,i,j,k) = valueinit;
    }
#endif

#endif


#if(PERCELLDT)
  GLOBALPOINT(dtijk) = (FTYPE PTRMACP0A1(dtijk,N1M,N2M,N3M,COMPDIM)) (&(BASEMACP0A1(dtijk,N1BND,N2BND,N3BND,-1))); // so access like dtijk[1,2,3]
  FULLLOOP for(l=1;l<=COMPDIM;l++){
    GLOBALMACP0A1(dtijk,i,j,k,l) = -1; //valueinit; // SUPERGODMARK: Can improve
  }
#endif


#if(STOREWAVESPEEDS>0)
  GLOBALPOINT(wspeedtemp) = (FTYPE PTRMACP1A1(wspeedtemp,NUMEOMSETS,N1M,N2M,N3M,NUMCS)) (&(BASEMACP1A1(wspeedtemp,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(wspeed) = (FTYPE PTRMACP3A0(wspeed,FILL,COMPDIM,NUMCS,N1M,N2M,N3M)) (&(BASEMACP3A0(wspeed,0,-1,0,N1BND,N2BND,N3BND))); // shifted so wspeed[1,2,3] accesses the memory
  FULLLOOP for(qq=0;qq<NUMEOMSETS;qq++) for(l=1;l<=COMPDIM;l++) for(m=0;m<NUMCS;m++) GLOBALMACP3A0(wspeed,qq,l,m,i,j,k) = valueinit;
#endif
  
#if(STORESHOCKINDICATOR)
  GLOBALPOINT(shockindicatorarray) = (FTYPE PTRMACP1A0(shockindicatorarray,NUMSHOCKPLS,N1M,N2M,N3M)) (&(BASEMACP1A0(shockindicatorarray,-SHOCKPLDIR1,N1BND,N2BND,N3BND)));
  FULLLOOP for(l=SHOCKPLDIR1;l<SHOCKPLDIR1+NUMSHOCKPLS;l++) GLOBALMACP1A0(shockindicatorarray,l,i,j,k) = valueinit;
#endif
  



  ////////////////////////////////////////////////
  //
  // TIME-STEPPING
  //
  ////////////////////////////////////////////////

  GLOBALPOINT(pk) = (FTYPE PTRMACP1A1(pk,FILL,N1M,N2M,N3M,NPR)) (&(BASEMACP1A1(pk,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(uinitialglobal) = (FTYPE PTRMACP0A1(uinitialglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(uinitialglobal,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(ulastglobal) = (FTYPE PTRMACP0A1(ulastglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(ulastglobal,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(unewglobal) = (FTYPE PTRMACP0A1(unewglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(unewglobal,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(dUgeomarray) = (FTYPE PTRMACP0A1(dUgeomarray,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(dUgeomarray,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    DTSTAGELOOP(dtstage) GLOBALMACP1A1(pk,dtstage,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(uinitialglobal,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(ulastglobal,i,j,k,pl) = valueinit;
    if(FULLOUTPUT==0) GLOBALMACP0A1(unewglobal,i,j,k,pl) = valueinit;
    else  GLOBALMACP0A1(unewglobal,i,j,k,pl) = 0;
    GLOBALMACP0A1(dUgeomarray,i,j,k,pl) = valueinit;
  }      

  
#if(HIGHERORDERMEM||FIELDSTAGMEM) // upoint needed for FV method and STAG for all methods
  GLOBALPOINT(upointglobal) = (FTYPE PTRMACP0A1(upointglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(upointglobal,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(upointglobal,i,j,k,pl) = valueinit;

  GLOBALPOINT(upointglobaluf) = (FTYPE PTRMACP0A1(upointglobaluf,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(upointglobaluf,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(upointglobaluf,i,j,k,pl) = valueinit;

  GLOBALPOINT(oldufstore) = (FTYPE PTRMACP0A1(oldufstore,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(oldufstore,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(oldufstore,i,j,k,pl) = valueinit;
#endif



  ////////////////////////////////////////////////
  //
  // SPATIAL INTERPOLATION
  //
  ////////////////////////////////////////////////

  // Note that PTR does not have N?M+SHIFT? unlike BASE pointer.  This is so flux_ct() can avoid extra loops and still not seg fault when accessing beyind F?'s normal range.
  // Since we want extra space at *bottom* (So can access (e.g.) F1[-N1BND-1] without seg faulting)
  // For FLUXCD, want extra space at top, so in general add 2 extra spaces (one at bottom and one at top)
  // Thus, only do SHIFT+NBND for each since that shifts bottom enough, and top space will be there.
#if(N1>1)
  GLOBALPOINT(F1) = (FTYPE PTRMACP0A1(F1,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F1,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F1,i,j,k,pl) = valueinit;
  }
#endif
#if(N2>1)
  GLOBALPOINT(F2) = (FTYPE PTRMACP0A1(F2,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F2,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F2,i,j,k,pl) = valueinit;
  }
#endif
#if(N3>1)
  GLOBALPOINT(F3) = (FTYPE PTRMACP0A1(F3,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F3,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F3,i,j,k,pl) = valueinit;
  }
#endif

#if(SPLITMAEMMEM)
#if(N1>1)
  GLOBALPOINT(F1EM) = (FTYPE PTRMACP0A1(F1EM,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F1EM,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F1EM,i,j,k,pl) = valueinit;
  }
#endif
#if(N2>1)
  GLOBALPOINT(F2EM) = (FTYPE PTRMACP0A1(F2EM,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F2EM,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F2EM,i,j,k,pl) = valueinit;
  }
#endif
#if(N3>1)
  GLOBALPOINT(F3EM) = (FTYPE PTRMACP0A1(F3EM,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(F3EM,SHIFT1+N1BND,SHIFT2+N2BND,SHIFT3+N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(F3EM,i,j,k,pl) = valueinit;
  }
#endif
#endif


  
#if(SPLITNPR||FIELDSTAGMEM)
  GLOBALPOINT(gp_l) = (FTYPE PTRMACP1A1(gp_l,FILL,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP1A1(gp_l,-1,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(gp_r) = (FTYPE PTRMACP1A1(gp_r,FILL,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP1A1(gp_r,-1,N1BND,N2BND,N3BND,0)));

  FULLLOOP PINTERPLOOP(pliter,pl){
    for(l=1;l<=3;l++){
      GLOBALMACP1A1(gp_l,l,i,j,k,pl) = valueinit;
      GLOBALMACP1A1(gp_r,l,i,j,k,pl) = valueinit;
    }
  }
  
#endif


  GLOBALPOINT(pleft) = (FTYPE PTRMACP0A1(pleft,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(pleft,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(pright) = (FTYPE PTRMACP0A1(pright,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(pright,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(prc) = (FTYPE PTRMACP0A1(prc,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(prc,N1BND,N2BND,N3BND,0)));
  
  FULLLOOP PINTERPLOOP(pliter,pl){
    GLOBALMACP0A1(pleft,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(pright,i,j,k,pl) = valueinit;
    GLOBALMACP0A1(prc,i,j,k,pl) = valueinit;
  }
  



#if(FIELDSTAGMEM)
  //  GLOBALPOINT(wspeedcorn) = (FTYPE PTRMACP2A0(wspeedcorn,FILL,NUMCS,N1M,N2M,N3M)) (&(BASEMACP2A0(wspeedcorn,-1,0,N1BND,N2BND,N3BND))); // shifted so wspeedcorn[1,2,3] accesses the memory
  GLOBALPOINT(pstagglobal) = (FTYPE PTRMACP0A1(pstagglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(pstagglobal,N1BND,N2BND,N3BND,0)));
  
  // -B1 and -U1 are so pbcorninterp[][B1] accesses 0th element of original memory (same for pvcorninterp)
  //  GLOBALPOINT(pbcorninterp) = (FTYPE PTRMACP3A0(pbcorninterp,FILL,COMPDIM,NUMCS,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMACP3A0(pbcorninterp,-1,-B1,0,N1BND,N2BND,N3BND)));
  GLOBALPOINT(pvbcorninterp) = (FTYPE PTRMACP1A3(pvbcorninterp,COMPDIM,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,COMPDIM,NUMCS+1,NUMCS)) (&(BASEMACP1A3(pvbcorninterp,-1,N1BND,N2BND,N3BND,-1,0,0))); // only shift by -1 since holds both U1+dir-1 and B1+dir-1 for dir=1,2,3 that will now just be accessed via "dir" alone

  GLOBALPOINT(geomcornglobal) = (FTYPE PTRMACP1A0(geomcornglobal,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMACP1A0(geomcornglobal,-1,N1BND,N2BND,N3BND))); // geomcornglobal[1,2,3 for CORN1,2,3]

#if(HIGHERORDERMEM)
  // using NPR so can use normal functions
  GLOBALPOINT(Bhatglobal) = (FTYPE PTRMACP0A1(Bhatglobal,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(Bhatglobal,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(Bhatglobal,i,j,k,pl)=valueinit;
  }

#if(ANALYTICMEMORY)
  GLOBALPOINT(Bhatanalytic) = (FTYPE PTRMACP0A1(Bhatanalytic,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(Bhatanalytic,N1BND,N2BND,N3BND,0)));
  FULLLOOP PLOOP(pliter,pl){
    GLOBALMACP0A1(Bhatanalytic,i,j,k,pl)=valueinit;
  }
#endif

#endif


  FULLLOOP PALLLOOP(pl) GLOBALMACP0A1(pstagglobal,i,j,k,pl) = valueinit;

  // corner things that have extra at top
  FULLLOOPP1{
    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=1;pl<=COMPDIM;pl++) for(m=0;m<NUMCS+1;m++) for(l=0;l<NUMCS;l++)  GLOBALMACP1A3(pvbcorninterp,pl2,i,j,k,pl,m,l)=valueinit;
    //    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=B1;pl<=B3;pl++) for(l=0;l<NUMCS;l++)   GLOBALMACP3A0(pbcorninterp,pl2,pl,l,i,j,k)=valueinit;

    for(pl2=1;pl2<=COMPDIM;pl2++) GLOBALMACP1A0(geomcornglobal,pl2,i,j,k)=valueinit;
  }

#endif




  ////////////////////////////////////////////////
  //
  // OLD SPATIAL INTERPOLATION (and new way to store shock indicator for paraflag or paraline)
  //
  ////////////////////////////////////////////////

#if(DODQMEMORY||STORESHOCKINDICATOR)
#if(N1>1)
  GLOBALPOINT(dq1) = (FTYPE PTRMACP0A1(dq1,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(dq1,N1BND,N2BND,N3BND,0)));
  FULLLOOP PINTERPLOOP(pliter,pl) GLOBALMACP0A1(dq1,i,j,k,pl) = valueinit;
#endif
#if(N2>1)
  GLOBALPOINT(dq2) = (FTYPE PTRMACP0A1(dq2,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(dq2,N1BND,N2BND,N3BND,0)));
  FULLLOOP PINTERPLOOP(pliter,pl) GLOBALMACP0A1(dq2,i,j,k,pl) = valueinit;
#endif
#if(N3>1)
  GLOBALPOINT(dq3) = (FTYPE PTRMACP0A1(dq3,N1M,N2M,N3M,NPR2INTERP)) (&(BASEMACP0A1(dq3,N1BND,N2BND,N3BND,0)));
  FULLLOOP PINTERPLOOP(pliter,pl) GLOBALMACP0A1(dq3,i,j,k,pl) = valueinit;
#endif
#endif



  


  ////////////////////////////////////////////////
  //
  // HIGHER ORDER STUFF
  //
  ////////////////////////////////////////////////

#if( HIGHERORDERMEM )
  GLOBALPOINT(fluxvectemp) = (FTYPE PTRMACP0A1(fluxvectemp,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(fluxvectemp,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(Fa) = (FTYPE PTRMACP0A1(Fa,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(Fa,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(Fb) = (FTYPE PTRMACP0A1(Fb,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(Fb,N1BND,N2BND,N3BND,0)));
  trifprintf("Allocating Fa/Fb for UNSPLIT (ALL NOW) flux method\n");
  GLOBALPOINT(stencilvartemp) = (FTYPE PTRMACP0A1(stencilvartemp,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(stencilvartemp,N1BND,N2BND,N3BND,0)));

  // Note that this still references Fa,Fb to save memory
  GLOBALPOINT(a2cin)  = (FTYPE PTRMACP0A1(Fa,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(Fa,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(a2cout) = (FTYPE PTRMACP0A1(Fb,N1M,N2M,N3M,NPR+NSPECIAL)) (&(BASEMACP0A1(Fb,N1BND,N2BND,N3BND,0)));
  trifprintf("Allocating a2cin/a2cout for UNSPLIT (ALL NOW) a2c method\n");

#endif






#if(FLUXDUMP)
  GLOBALPOINT(fluxdump)=(FTYPE PTRMACP0A1(fluxdump,N1M,N2M,N3M,NUMFLUXDUMP)) (&(BASEMACP0A1(fluxdump,N1BND,N2BND,N3BND,0)));

  // normal things
  FULLLOOP for(pl=0;pl<NUMFLUXDUMP;pl++) GLOBALMACP0A1(fluxdump,i,j,k,pl)=0.0; // not used for evolution -- just dumping -- so ok to ignore if leaking(?)

#endif



  ////////////////////////////////////////////////
  //
  // DEBUG STUFF USUALLY ON
  //
  ////////////////////////////////////////////////

  GLOBALPOINT(pflag) = (PFTYPE PTRMACP0A1(pflag,N1M,N2M,N3M,NUMPFLAGS)) (&(BASEMACP0A1(pflag,N1BND,N2BND,N3BND,0)));
  FULLLOOP  PFLAGLOOP(pf) GLOBALMACP0A1(pflag,i,j,k,pf) = NANPFLAG;

  GLOBALPOINT(pflagfailorig) = (PFTYPE PTRMACP0A1(pflagfailorig,N1M,N2M,N3M,NUMFAILPFLAGS)) (&(BASEMACP0A1(pflagfailorig,N1BND,N2BND,N3BND,0)));
  FULLLOOP  FAILPFLAGLOOP(pf) GLOBALMACP0A1(pflagfailorig,i,j,k,pf) = NANPFLAG;

#if(DODEBUG)
  GLOBALPOINT(failfloorcount) = (CTYPE PTRMACP0A3(failfloorcount,N1M,N2M,N3M,2,NUMTSCALES,NUMFAILFLOORFLAGS)) (&(BASEMACP0A3(failfloorcount,N1BND,N2BND,N3BND,0,0,0)));
  FULLLOOP  FAILFLOORLOOP(indexfinalstep,tscale,floor) GLOBALMACP0A3(failfloorcount,i,j,k,indexfinalstep,tscale,floor)=valueinit;
#endif

#if(DOFLOORDIAG)
  GLOBALPOINT(failfloordu) = (FTYPE PTRMACP0A1(failfloordu,N1M,N2M,N3M,NPR)) (&(BASEMACP0A1(failfloordu,N1BND,N2BND,N3BND,0)));
  FULLLOOP{
    PALLLOOP(pl) GLOBALMACP0A1(failfloordu,i,j,k,pl)=valueinit;
  }
#endif
#if(DODISSMEASURE)
  GLOBALPOINT(dissmeasurearray) = (FTYPE PTRMACP0A1(dissmeasurearray,N1M,N2M,N3M,NSPECIAL+1+3*2)) (&(BASEMACP0A1(dissmeasurearray,N1BND,N2BND,N3BND,0)));
  FULLLOOP{
    //    for(pl=0;pl<NSPECIAL+1+3*2;pl++) GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=valueinit;
    for(pl=0;pl<NSPECIAL+1+3*2;pl++) GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0; //  __WORKINGONIT__ For now so can use past measure to help set Fi and Firad
    // for dump to be clean for unused things
    if(N1==1){
      dir=1; pl=NSPECIAL+1+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
      dir=1; pl=NSPECIAL+1+3+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
    }
    if(N2==1){
      dir=2; pl=NSPECIAL+1+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
      dir=2; pl=NSPECIAL+1+3+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
    }
    if(N3==1){
      dir=3; pl=NSPECIAL+1+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
      dir=3; pl=NSPECIAL+1+3+dir-1; GLOBALMACP0A1(dissmeasurearray,i,j,k,pl)=0.0;
    }
  }
#endif

  ////////////////////////////////////////////////
  //
  // other diagnostics
  //
  ////////////////////////////////////////////////


#if(DODISS)
  GLOBALPOINT(dissfunpos) = (FTYPE PTRMACP0A1(dissfunpos,N1M,N2M,N3M,NUMDISSFUNPOS)) (&(BASEMACP0A1(dissfunpos,N1BND,N2BND,N3BND,0)));
#endif


#if(CALCFARADAYANDCURRENTS)
  // this faraday needed for current calculation
  GLOBALPOINT(cfaraday) = (FTYPE PTRMACP0A2(cfaraday,N1M,N2M,N3M,NUMCURRENTSLOTS,3)) (&(BASEMACP0A2(cfaraday,N1BND,N2BND,N3BND,0,0)));
  GLOBALPOINT(fcon) = (FTYPE PTRMACP0A1(fcon,N1M,N2M,N3M,NUMFARADAY)) (&(BASEMACP0A1(fcon,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(jcon) = (FTYPE PTRMACP0A1(jcon,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(jcon,N1BND,N2BND,N3BND,0)));

  FULLLOOP{
    for(pl=0;pl<NUMCURRENTSLOTS;pl++) for(l=0;l<3;l++){
#if(DOGRIDSECTIONING)
        // GODMARK: if grid moves while computing current, then near boundary current will be undefined for time derivative.
        // Could compute J^i using substeps, but harder to do.
        // So for now just assume border-region of current is poorly computed and just avoid nan's
        GLOBALMACP0A2(cfaraday,i,j,k,pl,l)=0.0;
#else
        GLOBALMACP0A2(cfaraday,i,j,k,pl,l)=valueinit;
#endif
      }
    for(pl=0;pl<NUMFARADAY;pl++){
      GLOBALMACP0A1(fcon,i,j,k,pl)=valueinit;
    }
    for(pl=0;pl<NDIM;pl++){
      GLOBALMACP0A1(jcon,i,j,k,pl)=valueinit;
    }
  }

#endif


  ////////////////////////////////////////////////
  //
  // AVG diagnostics
  //
  ////////////////////////////////////////////////
  
  // assume time average stuff gets zeroed in avg routine
#if(DOAVG)
  GLOBALPOINT(normalvarstavg) = (FTYPE PTRMACP0A1(normalvarstavg,N1M,N2M,N3M,NUMNORMDUMP)) (&(BASEMACP0A1(normalvarstavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(anormalvarstavg) = (FTYPE PTRMACP0A1(anormalvarstavg,N1M,N2M,N3M,NUMNORMDUMP)) (&(BASEMACP0A1(anormalvarstavg,N1BND,N2BND,N3BND,0)));
  
  GLOBALPOINT(fcontavg) = (FTYPE PTRMACP0A1(fcontavg,N1M,N2M,N3M,NUMFARADAY)) (&(BASEMACP0A1(fcontavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(fcovtavg) = (FTYPE PTRMACP0A1(fcovtavg,N1M,N2M,N3M,NUMFARADAY)) (&(BASEMACP0A1(fcovtavg,N1BND,N2BND,N3BND,0)));
  
  GLOBALPOINT(afcontavg) = (FTYPE PTRMACP0A1(afcontavg,N1M,N2M,N3M,NUMFARADAY)) (&(BASEMACP0A1(afcontavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(afcovtavg) = (FTYPE PTRMACP0A1(afcovtavg,N1M,N2M,N3M,NUMFARADAY)) (&(BASEMACP0A1(afcovtavg,N1BND,N2BND,N3BND,0)));

  GLOBALPOINT(massfluxtavg) = (FTYPE PTRMACP0A1(massfluxtavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(massfluxtavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(amassfluxtavg) = (FTYPE PTRMACP0A1(amassfluxtavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(amassfluxtavg,N1BND,N2BND,N3BND,0)));

  GLOBALPOINT(othertavg) = (FTYPE PTRMACP0A1(othertavg,N1M,N2M,N3M,NUMOTHER)) (&(BASEMACP0A1(othertavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(aothertavg) = (FTYPE PTRMACP0A1(aothertavg,N1M,N2M,N3M,NUMOTHER)) (&(BASEMACP0A1(aothertavg,N1BND,N2BND,N3BND,0)));

#if(CALCFARADAYANDCURRENTS)
  GLOBALPOINT(jcontavg) = (FTYPE PTRMACP0A1(jcontavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(jcontavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(jcovtavg) = (FTYPE PTRMACP0A1(jcovtavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(jcovtavg,N1BND,N2BND,N3BND,0)));

  GLOBALPOINT(ajcontavg) = (FTYPE PTRMACP0A1(ajcontavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(ajcontavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(ajcovtavg) = (FTYPE PTRMACP0A1(ajcovtavg,N1M,N2M,N3M,NDIM)) (&(BASEMACP0A1(ajcovtavg,N1BND,N2BND,N3BND,0)));
#endif

  GLOBALPOINT(tudtavg) = (FTYPE PTRMACP0A1(tudtavg,N1M,N2M,N3M,NUMSTRESSTERMS)) (&(BASEMACP0A1(tudtavg,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(atudtavg) = (FTYPE PTRMACP0A1(atudtavg,N1M,N2M,N3M,NUMSTRESSTERMS)) (&(BASEMACP0A1(atudtavg,N1BND,N2BND,N3BND,0)));
#endif  





  /* grid functions */
  // GODMARK: for axisymmetric space-times, may want to keep metric functions 2D to save memory
  
  // these have 1 extra value on outer edge.  Shift for real pointer no different
#if(NEWMETRICSTORAGE)
  // new way
  GLOBALPOINT(gdetgeom) = (struct of_gdetgeom PTRMETMACP0A1(gdetgeom,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NPG)) (&(BASEMETMACP0A1(gdetgeom,N1BND,N2BND,N3BND,0)));

  GLOBALPOINT(gdetgeomnormal) = (struct of_gdetgeom PTRMETMACP1A0(gdetgeomnormal,NPG,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(gdetgeomnormal,0,N1BND,N2BND,N3BND)));

  GLOBALPOINT(compgeom) = (struct of_compgeom PTRMETMACP1A0(compgeom,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(compgeom,0,N1BND,N2BND,N3BND)));
#if(DOEVOLVEMETRIC)
  GLOBALPOINT(compgeomlast) = (struct of_compgeom PTRMETMACP1A0(compgeomlast,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(compgeomlast,0,N1BND,N2BND,N3BND)));
#endif
  
#else //else if old way
  GLOBALPOINT(gcon) = (FTYPE PTRMETMACP1A1(gcon,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,SYMMATRIXNDIM)) (&(BASEMETMACP1A1(gcon,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(gcov) = (FTYPE PTRMETMACP1A1(gcov,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,SYMMATRIXNDIM)) (&(BASEMETMACP1A1(gcov,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(gcovpert) = (FTYPE PTRMETMACP1A1(gcovpert,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM)) (&(BASEMETMACP1A1(gcovpert,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(gdet) = (FTYPE PTRMETMACP1A0(gdet,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(gdet,0,N1BND,N2BND,N3BND)));

#if(WHICHEOM!=WITHGDET)
  GLOBALPOINT(eomfunc) = (FTYPE PTRMETMACP1A1(eomfunc,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NPR)) (&(BASEMETMACP1A1(eomfunc,0,N1BND,N2BND,N3BND,0)));
#endif
#if(GDETVOLDIFF)
  GLOBALPOINT(gdetvol) = (FTYPE PTRMETMACP1A0(gdetvol,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(gdetvol,0,N1BND,N2BND,N3BND)));
#endif
  GLOBALPOINT(alphalapse) = (FTYPE PTRMETMACP1A0(alphalapse,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(alphalapse,0,N1BND,N2BND,N3BND)));
  GLOBALPOINT(betasqoalphasq) = (FTYPE PTRMETMACP1A0(betasqoalphasq,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(betasqoalphasq,0,N1BND,N2BND,N3BND)));
  GLOBALPOINT(beta) = (FTYPE PTRMETMACP1A1(beta,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM)) (&(BASEMETMACP1A1(beta,0,N1BND,N2BND,N3BND,0)));

#if(DOEVOLVEMETRIC)
  GLOBALPOINT(gcovlast) = (FTYPE PTRMETMACP1A1(gcovlast,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,SYMMATRIXNDIM)) (&(BASEMETMACP1A1(gcovlast,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(gcovpertlast) = (FTYPE PTRMETMACP1A1(gcovpertlast,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM)) (&(BASEMETMACP1A1(gcovpertlast,0,N1BND,N2BND,N3BND,0)));
  GLOBALPOINT(alphalapselast) = (FTYPE PTRMETMACP1A0(alphalapselast,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3)) (&(BASEMETMACP1A0(alphalapselast,0,N1BND,N2BND,N3BND)));
#endif

#endif // end if old way
  
  
#if(DOSTOREPOSITIONDATA)
  GLOBALPOINT(dxdxpstore) = (FTYPE PTRMETMACP1A2(dxdxpstore,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM)) (&(BASEMETMACP1A2(dxdxpstore,0,N1BND,N2BND,N3BND,0,0)));
  GLOBALPOINT(idxdxpstore) = (FTYPE PTRMETMACP1A2(idxdxpstore,FILL,N1M+SHIFT1,N2M+SHIFT2,N3M+SHIFT3,NDIM,NDIM)) (&(BASEMETMACP1A2(idxdxpstore,0,N1BND,N2BND,N3BND,0,0)));

  GLOBALPOINT(Xstore) = (FTYPE PTRMACP1A1(Xstore,FILL,N1M+SHIFT1*3,N2M+SHIFT2*3,N3M+SHIFT3*3,NDIM)) (&(BASEMACP1A1(Xstore,0,N1BND+SHIFT1,N2BND+SHIFT2,N3BND+SHIFT3,0))); // Xstore is not reduced if CARTMINKMETRIC
  GLOBALPOINT(Vstore) = (FTYPE PTRMACP1A1(Vstore,FILL,N1M+SHIFT1*3,N2M+SHIFT2*3,N3M+SHIFT3*3,NDIM)) (&(BASEMACP1A1(Vstore,0,N1BND+SHIFT1,N2BND+SHIFT2,N3BND+SHIFT3,0))); // Vstore is not reduced if CARTMINKMETRIC
#endif
  
  
  // rest are always located at CENT
  GLOBALPOINT(conn) = (FTYPE PTRMETMACP0A3(conn,N1M,N2M,N3M,NDIM,NDIM,NDIM)) (&(BASEMETMACP0A3(conn,N1BND,N2BND,N3BND,0,0,0)));
  GLOBALPOINT(conn2) = (FTYPE PTRMETMACP0A1(conn2,N1M,N2M,N3M,NDIM)) (&(BASEMETMACP0A1(conn2,N1BND,N2BND,N3BND,0)));
  
#if(VOLUMEDIFF)
  GLOBALPOINT(idxvol) = (FTYPE PTRMETMACP0A1(idxvol,N1M,N2M,N3M,NDIM)) (&(BASEMETMACP0A1(idxvol,N1BND,N2BND,N3BND,0)));
#endif
  
  
  // initialize global pointers for multi-dimensional arrays for radiation
  // KORAL
  set_arrays_multidimen_rad();
  

}
