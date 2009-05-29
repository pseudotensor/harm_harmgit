#include "decs.h"


int bound_mpi_int(int boundstage, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS])
{
  int dir;

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  ///////////////
  // dir=1  
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE1)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }

  ///////////////
  // dir=2
  if((boundstage==STAGE2)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE3)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }

  ///////////////
  // dir=3
  if((boundstage==STAGE4)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE5)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirgenset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }


  // now corner zones will be filled correctly
  // GODMARK: If made fixup_utoprim() and check_solution() not use corner zones could bound all directions at once -- probably not important performance hit


  // end if mpi
#endif

  return(0);

}	
// end function


#define PACKLOOP_INT(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) FBOUNDLOOP(pl) SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk)

// packs data for shipment
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

#if(USEMPI)
void sendrecv_int(int dir, int boundvartype,PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc_int[UNPACK][dir],
	    dirgenset[boundvartype][dir][DIRSIZE],
	    MPI_PFTYPE,
	    MPIid[dirgenset[boundvartype][dir][DIROTHER]],
	    dirgenset[boundvartype][dir][DIRTAGR],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc_int[PACK][dir],
	    dirgenset[boundvartype][dir][DIRSIZE],
	    MPI_PFTYPE,
	    MPIid[dirgenset[boundvartype][dir][DIROTHER]],
	    dirgenset[boundvartype][dir][DIRTAGS],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQSEND]);

}
#endif


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

