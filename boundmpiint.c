#include "decs.h"

/*! \file boundmpiint.c
  \brief Boundary conditions using MPI for core-overlapping cells for pflag or other integer quantities

*/


/// see boundmpi.c for more comments
int bound_mpi_int_dir(int boundstage, int finalstep, int whichdir, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS])
{
  int dir;
  static MPI_Request requests[COMPDIM * 2 * 2];
  static int didpostrecvs[COMPDIM*2]={0};

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
          recvonly_int(dir,boundvartype,workbc_int,requests);
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
          pack_int(dir,boundvartype,prim,workbc_int);
        }
      }
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendonly_int(dir,boundvartype,workbc_int,requests);
    }
    if((boundstage==STAGEM1)||(boundstage==STAGE1&&whichdir==1 || boundstage==STAGE3&&whichdir==2 || boundstage==STAGE5&&whichdir==3)){
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]){
          recvwait(dir,requests);
          didpostrecvs[dir]=0; // done with recv's
        }
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
      for(dir=dirstart;dir<=dirfinish;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
    }
  }
  else{
    dualfprintf(fail_file,"No such whichdir=%d in boundmpiint.c\n",whichdir);
    myexit(1986290387);
  }



  // now corner zones will be filled correctly
  // GODMARK: If made fixup_utoprim() and check_solution() not use corner zones could bound all directions at once -- probably not important performance hit


  return(0);

} 
// end function


#define PACKLOOP_INT(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) FBOUNDLOOP(pl) SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk)

/// packs data for shipment
void pack_int(int dir, int boundvartype,PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS],PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int pl,pliter;
  int bci;


  bci=0;
  PACKLOOP_INT(i,j,k
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTART1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTOP1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTART2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTOP2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTART3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPSTOP3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPDIR1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPDIR2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRPDIR3]
               ){
    /*
      if(bci>=dirgenset[boundvartype][dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirgenset[%d][DIRSIZE]: %d\n",bci,dirgenset[boundvartype][dir][DIRSIZE]);
      myexit(10);
      }
    */
    workbc_int[PACK][dir][bci++] = MACP0A1(prim,i,j,k,pl);
  }
}


void recvonly_int(int dir, int boundvartype,PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc_int[UNPACK][dir],
            dirgenset[boundvartype][dir][DIRSIZE],
            MPI_PFTYPE,
            MPIid[dirgenset[boundvartype][dir][DIROTHER]],
            TAGSTARTBOUNDMPIINT + dirgenset[boundvartype][dir][DIRTAGR],
            MPI_COMM_GRMHD,
            &requests[dir*2+REQRECV]);


}


void sendonly_int(int dir, int boundvartype,PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{


  if(MPIFLOWCONTROL==1){
    int nothingsend=0;
    int nothingrecv=0;
    int maxtag = numprocs*COMPDIM*2;
    
    MPI_Sendrecv(
                 &nothingsend,0,MPI_INT,
                 MPIid[dirgenset[boundvartype][dir][DIROTHER]],
                 TAGSTARTBOUNDMPIINT + maxtag + dirgenset[boundvartype][dir][DIRTAGS],

                 &nothingrecv,0,MPI_INT,
                 MPIid[dirgenset[boundvartype][dir][DIROTHER]],
                 TAGSTARTBOUNDMPIINT + maxtag + dirgenset[boundvartype][dir][DIRTAGR],
        
                 MPI_COMM_GRMHD,MPI_STATUS_IGNORE);
  } // end if doing FLOWCONTROL



  MPI_Isend(workbc_int[PACK][dir],
            dirgenset[boundvartype][dir][DIRSIZE],
            MPI_PFTYPE,
            MPIid[dirgenset[boundvartype][dir][DIROTHER]],
            TAGSTARTBOUNDMPIINT + dirgenset[boundvartype][dir][DIRTAGS],
            MPI_COMM_GRMHD,
            &requests[dir*2+REQSEND]);

}


void unpack_int(int dir, int boundvartype,PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS])
{
  // dir is direction receiving from
  int i,j,k;
  int pl,pliter;
  int bci;

  bci=0;
  PACKLOOP_INT(i,j,k
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR3]
               ){
    MACP0A1(prim,i,j,k,pl)=workbc_int[UNPACK][dir][bci++];
  }
}











/// fake MPI bound call so fills same locations with fakevalue
/// no actual MPI calls are made -- just uses the same structures for simplicity
int bound_mpi_int_fakeutoprimmpiinconsisent(int boundstage, int finalstep, int fakedir, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS], int fakevalue)
{
  int dir;


  if(fakedir==+1){ // no need for fakedir==-1
    // only need to "unpack"
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int_fakeutoprimmpiinconsisent(dir,boundvartype,workbc_int,prim,fakevalue);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int_fakeutoprimmpiinconsisent(dir,boundvartype,workbc_int,prim,fakevalue);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int_fakeutoprimmpiinconsisent(dir,boundvartype,workbc_int,prim,fakevalue);
  }
  
  return(0);

} 


/// fake unpack routine that just fills-in MPI boundary cells with fakevalue
void unpack_int_fakeutoprimmpiinconsisent(int dir, int boundvartype,PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS], int fakevalue)
{
  // dir is direction receiving from
  int i,j,k;
  int pl,pliter;
  int bci;

  bci=0;
  PACKLOOP_INT(i,j,k
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTART3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUSTOP3]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR1]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR2]
               ,dirloopset[boundvartype][dir][primgridpos[boundvartype][dir][pl]][DIRUDIR3]
               ){
    MACP0A1(prim,i,j,k,pl)=fakevalue;
  }
}





/// bound all directions
int bound_mpi_int(int boundstage, int finalstep, int fakedir, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS])
{
  int bound_mpi_int_dir(int boundstage, int finalstep, int whichdir, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS]);
  int whichdir;

  // separate pre-post recv and intermediate user bound and then final other MPI stuff -- implies TAG space is globally separated for these MPI calls and any MPI stuff done by user!
  if(fakedir==-1){
    whichdir=-1; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
    whichdir=-2; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
    whichdir=-3; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
  }
  else if(fakedir==1){
    whichdir=1;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
    whichdir=2;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
    whichdir=3;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
  }
  else if(fakedir==0){
    whichdir=-1; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
    whichdir=1;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
    whichdir=-2; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
    whichdir=2;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
    whichdir=-3; bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim); // pre-post MPI recv's
    whichdir=3;  bound_mpi_int_dir(boundstage, finalstep, whichdir, boundvartype, prim);
  }

  return(0);
}


