#include "decs.h"

#define DEBUG 0





// bound all directions
int bound_mpi(int boundstage, int boundtype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR])
{
  int bound_mpi_dir(int boundstage, int whichdir, int boundtype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR]);
  int whichdir;

  whichdir=1; bound_mpi_dir(boundstage, whichdir, boundtype, prim, F1, F2, F3);
  whichdir=2; bound_mpi_dir(boundstage, whichdir, boundtype, prim, F1, F2, F3);
  whichdir=3; bound_mpi_dir(boundstage, whichdir, boundtype, prim, F1, F2, F3);

  return(0);
}





// boundtype specifies whether to bound scalar or to bound vector that is only needed to be bound along that direction
int bound_mpi_dir(int boundstage, int whichdir, int boundtype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR])
{
  FTYPE (*prim2bound[NDIM])[NSTORE2][NSTORE3][NPR];
  int dir;
#if(DEBUG)
  int i,j,k,pl,pliter;
#endif

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  /*
   *
   * 1. do outflow/inflow/etc. boundary conditions (in bounds)
   * 2. pack data into workbc arrays and do send/recv's
   * 3. for each transfer do a wait/unpack
   * 4. do all send waits
   *
   * NOTE: order is important so corner zones are assigned correctly
   *
   * workbc[PACK][][] is data that is being sent
   * workbc[UNPACK][][] is data that is being recvd
   *
   */

  /*
   *  -> higher values of x1
   * 2 -> lower values of x1
   *
   */

  
  /* bounds has already done non-MPI boundary conditions; now do MPI
   * boundary conditions */

  // must go in this order (0,1) then (2,3) then (4,5) or visa versa,
  // a single directions order doesn't matter.  This method of
  // ordering is as opposed to directly transfering the corner zones
  // to the corner CPU.

  // This may be faster since all transfers can proceed at once,
  // although this may be slower since no transfers can occur until
  // packing is completed.  This way packing and transfering occur
  // simultaneously.

  // Although l/r are packed together, since in the end we have to
  // wait for both l/r to complete, so equal time completion is
  // favored //over asynch completion.

  // Also, transfering corner zones with small message sizes increases
  // the importance of latency.

  // for 2D:
  // I choose left-right N?M first, then up/down N?M.  Could
  // just do N? for interior for L/R, but true boundary needs full
  // N?M exchanged since cpu sets boundary using normal bc code
  // which needs to get transfered to neight(i.e. currently if corner
  // has bctype 99/?  then doesn't do corner)

  // GODMARK: Make sure 3D makes sense (no extra things to do)

#if(DEBUG)
  PBOUNDLOOP(pliter,pl){
    FULLLOOP{
      MACP0A1(prim,i,j,k,pl)=-1-pl*100; // should be turned into myid-k*100
    }
    ZLOOP {
      MACP0A1(prim,i,j,k,pl)=myid-pl*100; // differentiates but clear per pr
      //      fprintf(log_file,"%d %d %d %d %21.15g\n",i,j,k,pl,MACP0A1(prim,i,j,k,pl));
    }
  }
#endif


  // this is designed to copy corners by indirectly copying them.  The
  // corners eventually get to corner-related CPUs by passing through
  // another cpu.  This required ordering the left/right and up/down
  // procedures.

  // one could copy the corners directly and get more bandwidth since
  // would transfers 2X as much data, but corners would transfer very
  // slowly alone, and basically has the same number of operations
  // required as does edge transfers.
  //  fprintf(fail_file,"innerboundhere1\n"); fflush(fail_file);





  if(boundtype==BOUNDPRIMTYPE || boundtype==BOUNDPRIMSIMPLETYPE || boundtype==BOUNDPSTAGTYPE || boundtype==BOUNDPSTAGSIMPLETYPE){
    prim2bound[1]=prim;
    prim2bound[2]=prim;
    prim2bound[3]=prim;
  }
  else if(boundtype==BOUNDFLUXTYPE || boundtype==BOUNDFLUXSIMPLETYPE){
    prim2bound[1]=F1;
    prim2bound[2]=F2;
    prim2bound[3]=F3;
  }
  else{
    dualfprintf(fail_file,"No such type of MPI bounding: boundtype=%d\n",boundtype);
    myexit(917616);
  }



  if(whichdir==1){

    ///////////////////
    //
    // x -dir
    //
    /////////////////
    if((boundstage==STAGEM1)||(boundstage==STAGE0)){
      for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound[whichdir],workbc);
      for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
    }
    //fprintf(fail_file,"innerboundhere2\n"); fflush(fail_file);
    if((boundstage==STAGEM1)||(boundstage==STAGE1)){
      for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) recvwait(dir,requests);
      for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound[whichdir]);
      for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendwait(dir,requests);
    }
  }
  else if(whichdir==2){
    ///////////////////
    //
    // y -dir
    //
    /////////////////
    //fprintf(fail_file,"innerboundhere3\n"); fflush(fail_file);
    if((boundstage==STAGEM1)||(boundstage==STAGE2)){
      // now dir=0,1(X1UP,X1DN) is done, so can start 2,3(X2UP,X2DN)
      for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound[whichdir],workbc);
      for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
    }
    //fprintf(fail_file,"innerboundhere4\n"); fflush(fail_file);
    if((boundstage==STAGEM1)||(boundstage==STAGE3)){
      for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) recvwait(dir,requests);
      for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound[whichdir]);
      for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendwait(dir,requests);
    }
    //fprintf(fail_file,"innerboundhere5\n"); fflush(fail_file);
  }
  else if(whichdir==3){
    ///////////////////
    //
    // z -dir
    //
    /////////////////
    //fprintf(fail_file,"innerboundhere3\n"); fflush(fail_file);
    if((boundstage==STAGEM1)||(boundstage==STAGE4)){
      // now dir=0,1,2,3(X1UP,X1DN,X2UP,X2DN) is done, so can start 4,5(X3UP,X3DN)
      for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) pack(dir,boundtype,prim2bound[whichdir],workbc);
      for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendrecv(dir,boundtype,workbc,requests);
    }
    //fprintf(fail_file,"innerboundhere4\n"); fflush(fail_file);
    if((boundstage==STAGEM1)||(boundstage==STAGE5)){
      for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) recvwait(dir,requests);
      for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) unpack(dir,boundtype,workbc,prim2bound[whichdir]);
      for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundtype][dir][DIRIF]) sendwait(dir,requests);
    }
    //fprintf(fail_file,"innerboundhere5\n"); fflush(fail_file);
  }
  else{
    dualfprintf(fail_file,"No such whichdir=%d in boundmpi.c\n",whichdir);
    myexit(1986290386);
  }

#if(DEBUG)
    fprintf(log_file,"\n\nafter\n\n");
    FULLLOOP{
      PBOUNDLOOP(pliter,pl){
        fprintf(log_file,"%d %d %d %d %21.15g\n",i,j,k,pl,MACP0A1(prim,i,j,k,pl));
      }
    }
    myexit(0);
#endif



  // end if mpi
#endif


    return(0);

}	
// end function


// Since pr used to control loop range, PLOOPMPI() must be at outer loop
#define PACKLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk,pr,num) PLOOPMPI(pr,num) SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk)

// packs data for shipment
void pack(int dir, int boundtype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int bci,pr;

  bci=0;
  PACKLOOP(i,j,k
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTART1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTOP1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTART2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTOP2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTART3]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPSTOP3]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPDIR1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPDIR2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRPDIR3]
	   ,pr
	   ,dirgenset[boundtype][dir][DIRNUMPR] // this element not used if SPILTNPR==1
	   ){
    /*
    if(bci>=dirgenset[boundtype][dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirgenset[%d][DIRSIZE]: %d\n",bci,dirgenset[boundtype][dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc[PACK][dir][bci++] = MACP0A1(prim,i,j,k,pr) * primfactor[boundtype][dir][primgridpos[boundtype][dir][pr]][PACK][pr];
  }
}

#if(USEMPI)
void sendrecv(int dir, int boundtype, FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  int truesize;

#if(1)
  // now always controlling range of quantities
  // must be consistent with how PLOOPMPI is setup
  if(boundtype==BOUNDPRIMTYPE || boundtype==BOUNDPRIMSIMPLETYPE || boundtype==BOUNDPSTAGTYPE || boundtype==BOUNDPSTAGSIMPLETYPE) truesize=dirgenset[boundtype][dir][DIRSIZE]/NPRBOUND;
  else if(boundtype==BOUNDFLUXTYPE || boundtype==BOUNDFLUXSIMPLETYPE) truesize=dirgenset[boundtype][dir][DIRSIZE]/NFLUXBOUND;
  else{
    dualfprintf(fail_file,"No such boundtype=%d in sendrecv() in boundmpi.c\n",boundtype);
    myexit(34867364);
  }
  truesize *= (nprboundend-nprboundstart+1);
  
#else
  truesize=dirgenset[boundtype][dir][DIRSIZE];
#endif

  MPI_Irecv(workbc[UNPACK][dir],
	    truesize,
	    MPI_FTYPE,
	    MPIid[dirgenset[boundtype][dir][DIROTHER]],
	    dirgenset[boundtype][dir][DIRTAGR],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc[PACK][dir],
	    truesize,
	    MPI_FTYPE,
	    MPIid[dirgenset[boundtype][dir][DIROTHER]],
	    dirgenset[boundtype][dir][DIRTAGS],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQSEND]);

}
#endif


#if(USEMPI)
void recvwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQRECV], &mpichstatus);

}
#endif


void unpack(int dir, int boundtype, FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  // dir is direction receiving from
  int i,j,k;
  int bci,pr;

  bci=0;
  PACKLOOP(i,j,k
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTART1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTOP1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTART2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTOP2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTART3]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUSTOP3]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUDIR1]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUDIR2]
	   ,dirloopset[boundtype][dir][primgridpos[boundtype][dir][pr]][DIRUDIR3]
	   ,pr
	   ,dirgenset[boundtype][dir][DIRNUMPR] // not used if doing general quantity loop
	   ){
    MACP0A1(prim,i,j,k,pr)=workbc[UNPACK][dir][bci++] * primfactor[boundtype][dir][primgridpos[boundtype][dir][pr]][UNPACK][pr];
  }
}

#if(USEMPI)

void sendwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQSEND], &mpichstatus);
}
#endif
