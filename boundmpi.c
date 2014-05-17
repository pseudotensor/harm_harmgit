#include "decs.h"


/*! \file boundmpi.c
  \brief Boundary conditions using MPI for core-overlapping cells for FTYPE quantities

*/



#define DEBUG 0



static int get_truesize(int dir, int boundvartype);





/// boundvartype specifies whether to bound scalar or to bound vector that is only needed to be bound along that direction
int bound_mpi_dir(int boundstage, int finalstep, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  FTYPE (*prim2bound[NDIM])[NSTORE2][NSTORE3][NPR];
  FTYPE (*flux2bound[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*vpot2bound[NDIM])[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3];
  int dir;
#if(DEBUG)
  int i,j,k,pl,pliter;
#endif

  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
  static int didpostrecvs[COMPDIM*2]={0};


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
      //      logfprintf("%d %d %d %d %21.15g\n",i,j,k,pl,MACP0A1(prim,i,j,k,pl));
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
  //  dualfprintf(fail_file,"innerboundhere1\n"); fflush(fail_file);



  ///////////
  //
  // initialize prim2bound and vpot2bound so interior code knows which to act upon
  //
  ///////////
  {
    int jj;
    for(jj=0;jj<NDIM;jj++){
      prim2bound[jj]=NULL;
      flux2bound[jj]=NULL;
      vpot2bound[jj]=NULL;
    }
  }


  ///////////
  //
  // set prim2bound or flux2bound or vpot2bound
  //
  ///////////
  if(boundvartype==BOUNDPRIMTYPE || boundvartype==BOUNDPRIMSIMPLETYPE || boundvartype==BOUNDPSTAGTYPE || boundvartype==BOUNDPSTAGSIMPLETYPE){
    prim2bound[1]=prim;
    prim2bound[2]=prim;
    prim2bound[3]=prim;
  }
  else if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE){
    flux2bound[1]=F1;
    flux2bound[2]=F2;
    flux2bound[3]=F3;
  }
  else if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE){
    vpot2bound[1]=vpot;
    vpot2bound[2]=vpot;
    vpot2bound[3]=vpot;
  }
  else{
    dualfprintf(fail_file,"No such type of MPI bounding: boundvartype=%d\n",boundvartype);
    myexit(917616);
  }


  int dirstart,dirfinish;
  if(whichdir==-1 || whichdir==1){ dirstart=X1UP; dirfinish=X1DN;}
  if(whichdir==-2 || whichdir==2){ dirstart=X2UP; dirfinish=X2DN;}
  if(whichdir==-3 || whichdir==3){ dirstart=X3UP; dirfinish=X3DN;}


  ////////////
  //
  // pre-post recv's
  //
  // OPTMARK: Could avoid dir-dependent MPI calls in step_ch.c (currently true!) and post all recv's for all dirs at once.
  // OPTMARK: Or setup MPIFLOWCONTROL==2 or SUPERMPI
  //
  ////////////
  if(whichdir==-1 || whichdir==-2 || whichdir==-3){
    if((boundstage==STAGEM1)||(boundstage==STAGE0)){
      for(dir=dirstart;dir<=dirfinish;dir++){
        if(dirgenset[boundvartype][dir][DIRIF]){
          recvonly(dir,boundvartype,workbc,requests);
          didpostrecvs[dir]=1;
        }
      }
    }
  }
  else if(whichdir==1 || whichdir==2 || whichdir==3){

    //////////////////
    //
    // per-whichdir (1,2,3) conditionals for bounding.  Assume whichdir=1 done first, then 2, then 3, so that corners are done correctly as done with normal bounds.
    //
    //////////////////
    
    ///////////////////
    //
    // x or y or z -dir
    // once dir=0,1(X1UP,X1DN) is done, so can start 2,3(X2UP,X2DN)
    //
    /////////////////
    if((boundstage==STAGEM1)||(boundstage==STAGE0&&whichdir==1 || boundstage==STAGE2&&whichdir==2 || boundstage==STAGE4&&whichdir==3)){
      for(dir=dirstart;dir<=dirfinish;dir++){
        if(dirgenset[boundvartype][dir][DIRIF]){
          if(didpostrecvs[dir]==0){
            dualfprintf(fail_file,"Did not post recv and tried to already pack: dir=%d\n",dir);
            myexit(234525155);
          }
          pack(dir,boundvartype,prim2bound[whichdir],flux2bound[whichdir],vpot2bound[whichdir],workbc);
        }
      }
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendonly(dir,boundvartype,workbc,requests);
    }
    if((boundstage==STAGEM1)||(boundstage==STAGE1&&whichdir==1 || boundstage==STAGE3&&whichdir==2 || boundstage==STAGE5&&whichdir==3)){
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]){
          recvwait(dir,requests);
          didpostrecvs[dir]=0; // done with recv's
        }
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack(dir,boundvartype,workbc,prim2bound[whichdir],flux2bound[whichdir],vpot2bound[whichdir]);
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
    }
  }
  else{
    dualfprintf(fail_file,"No such whichdir=%d in boundmpi.c\n",whichdir);
    myexit(1986290386);
  }







#if(DEBUG)
  logfprintf("\n\nafter\n\n");
  FULLLOOP{
    PBOUNDLOOP(pliter,pl){
      logfprintf("%d %d %d %d %21.15g\n",i,j,k,pl,MACP0A1(prim,i,j,k,pl));
    }
  }
  myexit(0);
#endif



  return(0);

} 
// end function


// Since pr used to control loop range, PLOOPMPI() must be at outer loop
#define PACKLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk,pr,num) PLOOPMPI(pr,num) SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk)

#define PACKLOOPORIG(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk,pr,num) PLOOPMPIORIG(pr,num) SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk)


// the last element not used if SPILTNPR==1
#define PACKBLOCK     i,j,k                                             \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART3] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP3] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR3] \
    ,pr                                                                 \
    ,dirgenset[boundvartype][dir][DIRNUMPR]



/// packs data for shipment
void pack(int dir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*flux)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int bci,pr;
  
  bci=0;
  if(prim!=NULL){
    PACKLOOP(i,j,k       \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART3] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP3] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR3] \
             ,pr       \
             ,dirgenset[boundvartype][dir][DIRNUMPR]){
      
      /*
        if(bci>=dirgenset[boundvartype][dir][DIRSIZE]){
        dualfprintf(fail_file,"pack memory leak: bci: %d dirgenset[%d][DIRSIZE]: %d\n",bci,dirgenset[boundvartype][dir][DIRSIZE]);
        myexit(10);
        }
      */

      workbc[PACK][dir][bci++] = MACP0A1(prim,i,j,k,pr) * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][PACK][pr];
      
    }// end of vpot==NULL loop
    
  }// end if prim!=NULL
  else{ // flux or vpot
    PACKLOOPORIG(i,j,k       \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTART3] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPSTOP3] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRPDIR3] \
                 ,pr       \
                 ,dirgenset[boundvartype][dir][DIRNUMPR]){
      if(vpot!=NULL){
        workbc[PACK][dir][bci++] = MACP1A0(vpot,pr,i,j,k) * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][PACK][pr];
      }
      else if(flux!=NULL){
        workbc[PACK][dir][bci++] = MACP1A0(flux,pr,i,j,k) * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][PACK][pr];
      }
      else{
        dualfprintf(fail_file,"No such pack type\n");
        myexit(982359235);
      }
    }// end vpot!=NULL loop
  }// end if vpot!=NULL
  
  
}


/// get true size of NPR-type sized objects
static int get_truesize(int dir, int boundvartype)
{
  int truesize;

#if(1)
  // now always controlling range of quantities
  // must be consistent with how PLOOPMPI is setup
  if(boundvartype==BOUNDPRIMTYPE || boundvartype==BOUNDPRIMSIMPLETYPE || boundvartype==BOUNDPSTAGTYPE || boundvartype==BOUNDPSTAGSIMPLETYPE) truesize=dirgenset[boundvartype][dir][DIRSIZE]/NPRBOUND;
  else if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE) truesize=dirgenset[boundvartype][dir][DIRSIZE]/NFLUXBOUND;
  else if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE) truesize=dirgenset[boundvartype][dir][DIRSIZE]/NDIM;
  else{
    dualfprintf(fail_file,"No such boundvartype=%d in sendrecv() in boundmpi.c\n",boundvartype);
    myexit(34867364);
  }
  truesize *= (nprboundend-nprboundstart+1);
  
#else
  truesize=dirgenset[boundvartype][dir][DIRSIZE];
#endif
  
  return(truesize);
  
}

/// MPI recieve
void recvonly(int dir, int boundvartype, FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  int truesize;

  truesize=get_truesize(dir, boundvartype);

  MPI_Irecv(workbc[UNPACK][dir],
            truesize,
            MPI_FTYPE,
            MPIid[dirgenset[boundvartype][dir][DIROTHER]],
            TAGSTARTBOUNDMPI + dirgenset[boundvartype][dir][DIRTAGR],
            MPI_COMM_GRMHD,
            &requests[dir*2+REQRECV]);
  
}

/// MPI send
void sendonly(int dir, int boundvartype, FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  int truesize;

  truesize=get_truesize(dir, boundvartype);


  if(MPIFLOWCONTROL==1){
    // nominally, good to pre-post recv as here, but other id's Isend might have already been sent a while back causing the unexpected buffer to be filled.  Happens on NICS Kraken.
    // http://climate.ornl.gov/~rmills/pubs/mills-hoffman_CUG2009.pdf
    // The above suggests only pushing the Isend once a 0-byte blocking handshake has been completed from the receiver.  This ensures the relevant Irecv is already posted.

    // See also:
    // http://www.uni-due.de/imperia/md/content/css/ude_3_mpi.pdf
    // http://www.nics.tennessee.edu/computing-resources/kraken/mpi-tips-for-cray-xt5
    // http://www.nersc.gov/assets/Uploads/NERSCXTMPIEnvironment.pdf

    // Could place global barrier here, but then stops not just the pairs of recv/send but all procs.
  
    // So post blocking send of 0-byte size to same place as sending full data.
    int nothingsend=0;
    int nothingrecv=0;
    // maxtag should fit into integer (modern machines is -2147483648 to +2147483647 or maximum numprocs ~ 3.5E8 that is quite large enough for now.
    // otherwise use negative, which would work too.
    int maxtag = numprocs*COMPDIM*2; // see mpi_init.c DIRTAGS,DIRTAGR how set.
    
    MPI_Sendrecv(
                 &nothingsend,0,MPI_INT,
                 MPIid[dirgenset[boundvartype][dir][DIROTHER]],
                 TAGSTARTBOUNDMPI + maxtag + dirgenset[boundvartype][dir][DIRTAGS],

                 &nothingrecv,0,MPI_INT,
                 MPIid[dirgenset[boundvartype][dir][DIROTHER]],
                 TAGSTARTBOUNDMPI + maxtag + dirgenset[boundvartype][dir][DIRTAGR],
   
                 MPI_COMM_GRMHD,MPI_STATUS_IGNORE);
    
  } // end if doing FLOWCONTROL
  
  
  MPI_Isend(workbc[PACK][dir],
            truesize,
            MPI_FTYPE,
            MPIid[dirgenset[boundvartype][dir][DIROTHER]],
            TAGSTARTBOUNDMPI + dirgenset[boundvartype][dir][DIRTAGS],
            MPI_COMM_GRMHD,
            &requests[dir*2+REQSEND]);
}
  

void recvwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQRECV], &mpichstatus);
}


// last element not used if doing general quantity loop
#define UNPACKBLOCK     i,j,k                                           \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART3] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP3] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR1] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR2] \
    ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR3] \
    ,pr                                                                 \
    ,dirgenset[boundvartype][dir][DIRNUMPR]

/// unpack data from MPI transfer
void unpack(int dir, int boundvartype, FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*flux)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  // dir is direction receiving from
  int i,j,k;
  int bci,pr;

  bci=0;
  if(prim!=NULL){
    PACKLOOP(i,j,k \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART3] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP3] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR1] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR2] \
             ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR3] \
             ,pr \
             ,dirgenset[boundvartype][dir][DIRNUMPR]){

      MACP0A1(prim,i,j,k,pr)=workbc[UNPACK][dir][bci++] * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][UNPACK][pr];
    }
  }
  else{
    PACKLOOPORIG(i,j,k \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTART3] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUSTOP3] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR1] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR2] \
                 ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pr]][DIRUDIR3] \
                 ,pr \
                 ,dirgenset[boundvartype][dir][DIRNUMPR]){
      if(vpot!=NULL){
        MACP1A0(vpot,pr,i,j,k)=workbc[UNPACK][dir][bci++] * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][UNPACK][pr];
      }
      else if(flux!=NULL){
        MACP1A0(flux,pr,i,j,k)=workbc[UNPACK][dir][bci++] * primfactor[boundvartype][dir][primgridpos[boundvartype][dir][pr]][UNPACK][pr];
      }
      else{
        dualfprintf(fail_file,"No such unpack type\n");
        myexit(982359236);
      }
    }
  }

}


void sendwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQSEND], &mpichstatus);
}




/// bound all directions
int bound_mpi(int boundstage, int finalstep, int fakedir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int bound_mpi_dir(int boundstage, int finalstep, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3]);
  int whichdir;

  // separate pre-post recv and intermediate user bound and then final other MPI stuff -- implies TAG space is globally separated for these MPI calls and any MPI stuff done by user!
  if(fakedir==-1){
    whichdir=-1; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
    whichdir=-2; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
    whichdir=-3; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
  }
  else if(fakedir==1){
    whichdir=1;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
    whichdir=2;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
    whichdir=3;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
  }
  else if(fakedir==0){
    whichdir=-1; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
    whichdir=1;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
    whichdir=-2; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
    whichdir=2;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
    whichdir=-3; bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot); // pre-post MPI recv's
    whichdir=3;  bound_mpi_dir(boundstage, finalstep, whichdir, boundvartype, prim, F1, F2, F3, vpot);
  }

  return(0);
}


