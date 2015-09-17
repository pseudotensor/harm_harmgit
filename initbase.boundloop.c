
/*! \file initbase.boundloop.c
     \brief General initialization of code related to loops over boundary conditions.
*/


#include "decs.h"



int report_bound_loop(void)
{
  int inboundloop[NDIM];
  int outboundloop[NDIM];
  int innormalloop[NDIM];
  int outnormalloop[NDIM];
  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  int riin,riout,rjin,rjout,rkin,rkout;
  int dosetbc[COMPDIM*2];
  int ii,jj,dimen,dir;
  int boundvartype;
  int numbnd[NDIM],numnpr;
  


  for(boundvartype=0;boundvartype<NUMBOUNDTYPES;boundvartype++){

    trifprintf("boundvartype=%d\n",boundvartype);

    // get number of boundary cells and number of quantities bounding
    set_numbnd(boundvartype, numbnd, &numnpr);

    trifprintf("numnpr=%d\n",numnpr);
    DIMENLOOP(dimen){
      trifprintf("numbnd[%d]=%d\n",dimen,numbnd[dimen]);
    }

    set_boundloop(boundvartype, inboundloop, outboundloop, innormalloop, outnormalloop, inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);

    // now report
    trifprintf("Note inout/lowhigh iter has POINTDOWN=%d and POINTUP=%d\n",POINTDOWN,POINTUP);
    DIMENLOOP(dimen){
      for(ii=NUMUPDOWN-1;ii>=0;ii--){
        for(jj=NUMUPDOWN-1;jj>=0;jj--){
          trifprintf("inoutlohi[inoutboundary=%d][lowhighrange=%d][dimen=%d]=%d\n",ii,jj,dimen,inoutlohi[ii][jj][dimen]);
        }
      }
      trifprintf("inboundloop[dimen=%d]=%d outboundloop[dimen=%d]=%d\n",dimen,inboundloop[dimen],dimen,outboundloop[dimen]);
      if(dimen==1) trifprintf("riin=%d riout=%d\n",riin,riout);
      else if(dimen==2) trifprintf("rjin=%d rjout=%d\n",rjin,rjout);
      else if(dimen==3) trifprintf("rkin=%d rkout=%d\n",rkin,rkout);
    }
    DIRLOOP(dir){
      trifprintf("dosetbc[dir=%d]=%d\n",dir,dosetbc[dir]);
    }
  } // end over boundvartype


  return(0);

}







// set number of boundary cells and number of quantities for boundary cells for various types of quantities
// sets true number of boundary cells instead of just MAXBND if dimension exists
void set_numbnd(int boundvartype, int *numbnd, int *numnpr)
{


  
  if(boundvartype==BOUNDPRIMTYPE || boundvartype==BOUNDPSTAGTYPE){
    numbnd[1]=N1BND;
    numbnd[2]=N2BND;
    numbnd[3]=N3BND;
    *numnpr = NPRBOUND;
  }
  else if(boundvartype==BOUNDPRIMSIMPLETYPE || boundvartype==BOUNDPSTAGSIMPLETYPE){ // simple means only 1 cell
    // used to bound after fixups that only modify +-1 cell currently
    numbnd[1]=NUMPFLAGBND1;
    numbnd[2]=NUMPFLAGBND2;
    numbnd[3]=NUMPFLAGBND3;      
    *numnpr = NPRBOUND;
  }
  else if(boundvartype==BOUNDINTTYPE){ // always "simple"
    // right now only need 1 boundary cell for averaging or checking by check_solution() and fixup_utoprim()
    numbnd[1]=NUMPFLAGBND1;
    numbnd[2]=NUMPFLAGBND2;
    numbnd[3]=NUMPFLAGBND3;      
    *numnpr = NUMPFLAGSBOUND;
  }
  else if(boundvartype==BOUNDFLUXTYPE){
    numbnd[1]=N1BND;
    numbnd[2]=N2BND;
    numbnd[3]=N3BND;
    *numnpr = NFLUXBOUND;
  }
  else if(boundvartype==BOUNDFLUXSIMPLETYPE){
    // flux only defined with 1 shift
    numbnd[1]=SHIFT1;
    numbnd[2]=SHIFT2;
    numbnd[3]=SHIFT3;
    *numnpr = NFLUXBOUND;
  }
  else if(boundvartype==BOUNDVPOTTYPE){
    numbnd[1]=N1BND;
    numbnd[2]=N2BND;
    numbnd[3]=N3BND;
    *numnpr = NDIM;
  }
  else if(boundvartype==BOUNDVPOTSIMPLETYPE){
    // vpot only defined with one shift
    numbnd[1]=SHIFT1;
    numbnd[2]=SHIFT2;
    numbnd[3]=SHIFT3;
    *numnpr = NDIM;
  }
  else{
    dualfprintf(fail_file,"set_numbnd(): No such boundvartype=%d\n",boundvartype);
    myexit(2436726);
  }


}



// uses global variables local to this file
// determine orthogonal direction loop ranges depending upon presence of real or MPI boundary
// now that bound per direction shouldn't treat real and mpi boundaries differently -- always use full loop
// in order to avoid accessing undefined data, but still fill corner
// zones, the ORDER of boundary LOOPS is as follows:

// X1 in&out: LOOPN2 LOOPN3 LOOPBOUNDIN1 & LOOPBOUNDOUT1
// X2 in&out: LOOPF1 LOOPN3 LOOPBOUNDIN2 & LOOPBOUNDOUT2  // LOOPF1 ok if X1 dimension not there, then LOOPF1->LOOPN1
// X3 in&out: LOOPF1 LOOPF2 LOOPBOUNDIN3 & LOOPBOUNDOUT3  // as above

// above loops replaced with LOOPX1dir LOOPX2dir, LOOPX3dir

// the resulting shifts should include shift due to GRIDSECTION *and* due to normal boundaries
// so specifies really is,ie,js,je,ks,ke for normal loop
// boundvartype is one of NUMBOUNDTYPES in global.nondepmnemonics.h
// This sets up both LOOPBOUND?IN loops *and* LOOPX?dir loops
// Notes:
// 1) Assume set_boundloop() is called each time needed since BCs could change with time
void set_boundloop(int boundvartype, int *inboundloop, int*outboundloop, int*innormalloop, int*outnormalloop, int (*inoutlohi)[NUMUPDOWN][NDIM], int *riin, int *riout, int *rjin, int *rjout, int *rkin, int *rkout, int *dosetbc)
{
  int outflowtype[COMPDIM*2];
  int dir;
  int shifts[NUMUPDOWN][NDIM];
  int ii,jj;
  int dimen;
  int numbnd[NDIM],numnpr;
  int shiftamount[COMPDIM*2];
  const int NxNOT1[NDIM]={0,N1NOT1,N2NOT1,N3NOT1};
  const int Nx[NDIM]={0,N1,N2,N3};
 
  //////////////////////////////
  //
  // Get number of boundary cells and number of quantities bounding per cell
  //
  //////////////////////////////

  set_numbnd(boundvartype, numbnd, &numnpr);

  //////////////////////////////
  //
  // Determine whether shift boundaries with ACTIVEREGION for GRIDSECTIONING
  //
  //////////////////////////////

  shiftamount[X1DN]=SHIFTX1DN;
  shiftamount[X1UP]=SHIFTX1UP;
  shiftamount[X2DN]=SHIFTX2DN;
  shiftamount[X2UP]=SHIFTX2UP;
  shiftamount[X3DN]=SHIFTX3DN;
  shiftamount[X3UP]=SHIFTX3UP;


  if(GLOBALBCMOVEDWITHACTIVESECTION==1){
    DIRLOOP(dir){
      if(
         BCtype[dir]==OUTFLOW
         || BCtype[dir]==HORIZONOUTFLOW
         || BCtype[dir]==FIXEDOUTFLOW
         || BCtype[dir]==OUTFLOWNOINFLOW
         || BCtype[dir]==RAMESHOUTFLOW
         || BCtype[dir]==RESCALEOUTFLOW
         || BCtype[dir]==FREEOUTFLOW
         ){
        outflowtype[dir]=1;
      }
      else outflowtype[dir]=0;
    }

    // USERMARK can choose to force shifts to be always 0 instead of dependending upon BCs if user wants
    // SUPERGODMARK: For example, should add "FIXED" types and shift for those?
    DIRLOOP(dir){
      if(outflowtype[dir]) shifts[POINTFROMDIR(dir)][DIMEN(dir)]=shiftamount[dir];
      else  shifts[POINTFROMDIR(dir)][DIMEN(dir)]=0;
    
      // do set boundary physical conditions
      dosetbc[dir]=1;
    }

  }
  else{
    // then BCs are only true BCs
    // Then never should call COMP type loops in bounds.c
    DIRLOOP(dir){
      // check whether really next to X1DN boundary and need to set physical boundary conditions
      // e.g. SHIFTX1UP*(-1)<=MAXBND [note generically SHIFTX1UP is negative since added in COMP loops] will be true if near that boundary enough for physical boundary conditions to be used in comp loops
      if(shiftamount[dir]*(-DIRSIGN(dir))<=MAXBND) dosetbc[dir]=1;
      else dosetbc[dir]=0;

      // define physical boundary loop range
      if(dosetbc[dir]){
        // standard shift for physical boundary when doing that boundary
        shifts[POINTFROMDIR(dir)][DIMEN(dir)]=0;
      }
      else{
        shifts[POINTFROMDIR(dir)][DIMEN(dir)]=shiftamount[dir];

        // below is wrong:
        // shift so that never access interior of such boundary loops in case not using dosetbc[] to check
        // if(DIRSIGN(dir)==-1) shifts[POINTFROMDIR(dir)][DIMEN(dir)]=totalsize[DIMEN(dir)]+MAXBND;
        // else if(DIRSIGN(dir)==1) shifts[POINTFROMDIR(dir)][DIMEN(dir)]=-MAXBND;
      }

    }
    

  }




  //////////////////////////////
  //
  // For LOOPBOUND?IN/OUT loops:
  //
  //////////////////////////////

  //int shifts[0=inner boundary, 1=outer boundary][0=lower loop index, 1=higher (upper) loop index][dimension];
  if(boundvartype==BOUNDPRIMTYPE || boundvartype==BOUNDPRIMSIMPLETYPE || boundvartype==BOUNDINTTYPE){
    // centered-like quantities

    dimen=1;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][1]=INBOUNDLO1;
    inoutlohi[POINTDOWN][POINTUP][1]=INBOUNDHI1;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][1]=OUTBOUNDLO1;
    inoutlohi[POINTUP][POINTUP][1]=OUTBOUNDHI1;

    dimen=2;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][2]=INBOUNDLO2;
    inoutlohi[POINTDOWN][POINTUP][2]=INBOUNDHI2;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][2]=OUTBOUNDLO2;
    inoutlohi[POINTUP][POINTUP][2]=OUTBOUNDHI2;

    dimen=3;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][3]=INBOUNDLO3;
    inoutlohi[POINTDOWN][POINTUP][3]=INBOUNDHI3;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][3]=OUTBOUNDLO3;
    inoutlohi[POINTUP][POINTUP][3]=OUTBOUNDHI3;
  }
  else if(boundvartype==BOUNDPSTAGTYPE || boundvartype==BOUNDPSTAGSIMPLETYPE){
    // stag-like quantities
    // In general, say for IF3DSPCTHENMPITRANSFERATPOLE, then in transverse direction have to deal with face fields
    // So must include transverse faces on upper edge where face really ends being active cell and where CENT would be first boundary cell

    dimen=1;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][1]=INBOUNDLO1;
    inoutlohi[POINTDOWN][POINTUP][1]=INBOUNDHI1;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][1]=OUTBOUNDLO1;
    //inoutlohi[POINTUP][POINTDOWN][1]=OUTFACEBOUNDLO1;
    inoutlohi[POINTUP][POINTUP][1]=OUTBOUNDHI1;

    dimen=2;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][2]=INBOUNDLO2;
    inoutlohi[POINTDOWN][POINTUP][2]=INBOUNDHI2;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][2]=OUTBOUNDLO2;
    //    inoutlohi[POINTUP][POINTDOWN][2]=OUTFACEBOUNDLO2;
    inoutlohi[POINTUP][POINTUP][2]=OUTBOUNDHI2;

    dimen=3;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][3]=INBOUNDLO3;
    inoutlohi[POINTDOWN][POINTUP][3]=INBOUNDHI3;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][3]=OUTBOUNDLO3;
    //inoutlohi[POINTUP][POINTDOWN][3]=OUTFACEBOUNDLO3;
    inoutlohi[POINTUP][POINTUP][3]=OUTBOUNDHI3;
  }
  else if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE){
    // face-like quantities

    dimen=1;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][1]=INFACEBOUNDLO1;
    inoutlohi[POINTDOWN][POINTUP][1]=INFACEBOUNDHI1;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][1]=OUTFACEBOUNDLO1;
    inoutlohi[POINTUP][POINTUP][1]=OUTFACEBOUNDHI1;

    // dir=2
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][2]=INFACEBOUNDLO2;
    inoutlohi[POINTDOWN][POINTUP][2]=INFACEBOUNDHI2;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][2]=OUTFACEBOUNDLO2;
    inoutlohi[POINTUP][POINTUP][2]=OUTFACEBOUNDHI2;

    // dir=3
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][3]=INFACEBOUNDLO3;
    inoutlohi[POINTDOWN][POINTUP][3]=INFACEBOUNDHI3;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][3]=OUTFACEBOUNDLO3;
    inoutlohi[POINTUP][POINTUP][3]=OUTFACEBOUNDHI3;
  }
  else if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE){
    // face-like quantities for certain directions (GODMARK: Needs updating for true BCs on real boundaries (i.e. non-MPI boundaries), but not used so far)

    dimen=1;
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][1]=INFACEBOUNDLO1;
    inoutlohi[POINTDOWN][POINTUP][1]=INFACEBOUNDHI1;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][1]=OUTFACEBOUNDLO1;
    inoutlohi[POINTUP][POINTUP][1]=OUTFACEBOUNDHI1;

    // dir=2
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][2]=INFACEBOUNDLO2;
    inoutlohi[POINTDOWN][POINTUP][2]=INFACEBOUNDHI2;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][2]=OUTFACEBOUNDLO2;
    inoutlohi[POINTUP][POINTUP][2]=OUTFACEBOUNDHI2;

    // dir=3
    // IN:
    inoutlohi[POINTDOWN][POINTDOWN][3]=INFACEBOUNDLO3;
    inoutlohi[POINTDOWN][POINTUP][3]=INFACEBOUNDHI3;
    // OUT:
    inoutlohi[POINTUP][POINTDOWN][3]=OUTFACEBOUNDLO3;
    inoutlohi[POINTUP][POINTUP][3]=OUTFACEBOUNDHI3;
  }


  ///////////
  //
  // convert original MAXBND number of boundary cells to chosen number of boundary cells
  //
  ///////////
  dimen=1;
  inoutlohi[POINTDOWN][POINTDOWN][1]=INBOUNDLO1+N1BND-numbnd[dimen];
  inoutlohi[POINTUP][POINTUP][1]=OUTBOUNDHI1-N1BND+numbnd[dimen];
  
  dimen=2;
  inoutlohi[POINTDOWN][POINTDOWN][2]=INBOUNDLO2+N2BND-numbnd[dimen];
  inoutlohi[POINTUP][POINTUP][2]=OUTBOUNDHI2-N2BND+numbnd[dimen];
  
  dimen=3;
  inoutlohi[POINTDOWN][POINTDOWN][3]=INBOUNDLO3+N3BND-numbnd[dimen];
  inoutlohi[POINTUP][POINTUP][3]=OUTBOUNDHI3-N3BND+numbnd[dimen];


  
  ////////////
  //
  // now shift the inoutlohi[inner/outer boundary][low/high range][dir] using shifts[][dir]
  //
  ///////////
  // both low/high ranges [jj] get shifted same
  for(ii=0;ii<NUMUPDOWN;ii++){
    for(jj=0;jj<NUMUPDOWN;jj++){
      DIMENLOOP(dimen){
        inoutlohi[ii][jj][dimen] += shifts[ii][dimen];
      }
    }
  }



  // check if really setting BC
  DIRLOOP(dir){
    dimen=DIMEN(dir);
    int point=POINTFROMDIR(dir); // inner and outer boundaries
    if(point==POINTUP){
      if(inoutlohi[point][POINTDOWN][dimen]<-numbnd[dimen] || inoutlohi[point][POINTUP][2]<0){
        dosetbc[dir]=0;
      }
    }
    else if(point==POINTDOWN){
      if(inoutlohi[point][POINTDOWN][dimen]>Nx[dimen]-1 || inoutlohi[point][POINTUP][2]>Nx[dimen]-1+numbnd[dimen]){
        dosetbc[dir]=0;
      }
    }
  }


  ////////////////////
  //
  // Set reference locations user should use to copy boundary cells from
  //
  ////////////////////

  *riin=inoutlohi[POINTDOWN][POINTUP][1]+SHIFT1; // gives 0 as required if no shifts
  *rjin=inoutlohi[POINTDOWN][POINTUP][2]+SHIFT2;
  *rkin=inoutlohi[POINTDOWN][POINTUP][3]+SHIFT3;

  if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE){
    // would be with extra -1 but defined expanded definition of OUTFACEBOUNDLO1
    *riout=inoutlohi[POINTUP][POINTDOWN][1];
    *rjout=inoutlohi[POINTUP][POINTDOWN][2];
    *rkout=inoutlohi[POINTUP][POINTDOWN][3];
  }
  else if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE){
    // would be with extra -1 but defined expanded definition of OUTFACEBOUNDLO1
    *riout=inoutlohi[POINTUP][POINTDOWN][1];
    *rjout=inoutlohi[POINTUP][POINTDOWN][2];
    *rkout=inoutlohi[POINTUP][POINTDOWN][3];
  }
  else{
    *riout=inoutlohi[POINTUP][POINTDOWN][1]-SHIFT1;
    *rjout=inoutlohi[POINTUP][POINTDOWN][2]-SHIFT2;
    *rkout=inoutlohi[POINTUP][POINTDOWN][3]-SHIFT3;
  }


  //////////////////////////////
  //
  // For LOOPX?dir loops:
  //
  //////////////////////////////

  // 1|| below is forced because of how use LOOPX?dir loops in bounds
  // as related to need to set corner zones
  if(1||mycpupos[1]==0) inboundloop[1]=inoutlohi[POINTDOWN][POINTDOWN][1];  //INFULL1+shifts[POINTDOWN][1];
  else inboundloop[1]=0;
  if(1||mycpupos[1]==ncpux1-1) outboundloop[1]=inoutlohi[POINTUP][POINTUP][1];   // OUTFULL1+shifts[POINTUP][1];
  else outboundloop[1]=N1-1;

  if(1||mycpupos[2]==0) inboundloop[2]=inoutlohi[POINTDOWN][POINTDOWN][2]; // INFULL2+shifts[POINTDOWN][2];
  else inboundloop[2]=0;
  if(1||mycpupos[2]==ncpux2-1) outboundloop[2]=inoutlohi[POINTUP][POINTUP][2]; // OUTFULL2+shifts[POINTUP][2];
  else outboundloop[2]=N2-1;

  if(1||mycpupos[3]==0) inboundloop[3]=inoutlohi[POINTDOWN][POINTDOWN][3]; //INFULL3+shifts[POINTDOWN][3];
  else inboundloop[3]=0;
  if(1||mycpupos[3]==ncpux3-1) outboundloop[3]=inoutlohi[POINTUP][POINTUP][3]; // OUTFULL3+shifts[POINTUP][3];
  else outboundloop[3]=N3-1;



  //////////////////////////////
  //
  // For LOOPN? part of LOOPX?dir loops
  //
  //////////////////////////////

  if(boundvartype==BOUNDPRIMTYPE || boundvartype==BOUNDPRIMSIMPLETYPE || boundvartype==BOUNDINTTYPE){
    // centered-like quantities
    innormalloop[1]=0+shifts[POINTDOWN][1];
    outnormalloop[1]=N1-1+shifts[POINTUP][1];

    innormalloop[2]=0+shifts[POINTDOWN][2];
    outnormalloop[2]=N2-1+shifts[POINTUP][2];

    innormalloop[3]=0+shifts[POINTDOWN][3];
    outnormalloop[3]=N3-1+shifts[POINTUP][3];
  }
  else if(boundvartype==BOUNDPSTAGTYPE || boundvartype==BOUNDPSTAGSIMPLETYPE){
    // stag-like quantities
    // SHIFT1/2/3 added so that includes full staggered transverse field
    innormalloop[1]=0+shifts[POINTDOWN][1];
    outnormalloop[1]=N1-1+SHIFT1+shifts[POINTUP][1];

    innormalloop[2]=0+shifts[POINTDOWN][2];
    outnormalloop[2]=N2-1+SHIFT2+shifts[POINTUP][2];

    innormalloop[3]=0+shifts[POINTDOWN][3];
    outnormalloop[3]=N3-1+SHIFT3+shifts[POINTUP][3];

  }
  else if(boundvartype==BOUNDFLUXTYPE || boundvartype==BOUNDFLUXSIMPLETYPE){
    // face-like quantities
    // centered-like quantities
    innormalloop[1]=0+shifts[POINTDOWN][1];
    outnormalloop[1]=N1-1+shifts[POINTUP][1];

    innormalloop[2]=0+shifts[POINTDOWN][2];
    outnormalloop[2]=N2-1+shifts[POINTUP][2];

    innormalloop[3]=0+shifts[POINTDOWN][3];
    outnormalloop[3]=N3-1+shifts[POINTUP][3];

  }
  else if(boundvartype==BOUNDVPOTTYPE || boundvartype==BOUNDVPOTSIMPLETYPE){
    // face-like quantities
    // centered-like quantities
    innormalloop[1]=0+shifts[POINTDOWN][1];
    outnormalloop[1]=N1-1+shifts[POINTUP][1];

    innormalloop[2]=0+shifts[POINTDOWN][2];
    outnormalloop[2]=N2-1+shifts[POINTUP][2];

    innormalloop[3]=0+shifts[POINTDOWN][3];
    outnormalloop[3]=N3-1+shifts[POINTUP][3];

  }

}


int orders_set(void)
{
  int l;

  // maximum number of points *potentially* used for a given interpolation type
  // for a point (i), we assume that the first point potentially used is:
  // si = i - (int)(interporder/2)
  // for WENO4, then if i=2, then si=0 and list of potentially used points is 0,1,2,3 .  This is offset such that the i=2 point is correctly centered at the edge with 2 points on edge side.
  // notice that WENO4 only gives correctly centered left state.  Using WENO4, right state is found from right neighbor.
  // WENO4 useful for avg->point at an interface, but not for point->point at different location, which requires odd WENO.
  // for WENO5, then if i=2, then si=0 and list is 0,1,2,3,4 .  Here, point is assumed to be centered of middle (2) cell.

  // default is 0 to indicate no limiter set
  for(l=-NUMNEGINTERPS;l<NUMPOSINTERPS;l++){
    interporder[l]=0;
  }


  // non-limiters
  interporder[NLIM]=3;
  interporder[NLIMCENT]=3;
  interporder[NLIMUP]=3;
  interporder[NLIMDOWN]=3;

  interporder[DONOR]=1;
  interporder[VANL]=3;
  interporder[MINM]=3;
  interporder[MC]=3;
  interporder[PARA]=5;
  interporder[PARAFLAT]=7;
  interporder[MCSTEEP]=7;
  interporder[CSSLOPE]=5;
  interporder[MP5]=5;
  interporder[EPPM]=5;

  interporder[WENO3]=3;
  interporder[WENO4]=4;
  interporder[WENO5]=5;
  interporder[WENO5BND]=11;   //replace the 5 for WENO5 with 11 = 5 + (2+2 for reduction) + (1+1 for dp/p reduction) to have enough boundary zones for l&r stencil reduction; also see interpline.c:slope_lim_linetype()
  interporder[WENO5BNDPLUSMIN]=13; // added extra minimization of weights JCM
  interporder[WENO5FLAT]=7;
  interporder[WENO6]=6;
  interporder[WENO7]=7;
  interporder[WENO8]=8;
  interporder[WENO9]=9;
  interporder[ENO3]=5;
  interporder[ENO5]=9;

  interporder[PARALINE]=7;


  for(l=-NUMNEGINTERPS;l<NUMPOSINTERPS;l++){
    if(interporder[l]>MAXSPACEORDER){
      dualfprintf(fail_file,"MAXSPACEORDER=%d while interporder[%d]=%d\n",MAXSPACEORDER,l,interporder[l]);
      myexit(13);
    }
  }

  // check to make sure interpolation type is correctly set
  // checked in flux.c from now on for each type of interpolation
  //  if(DOEXTRAINTERP){
  //    if( (VARTOINTERP!=PRIMTOINTERP_VSQ)||(RESCALEINTERP!=1) ){
  //      dualfprintf(fail_file,"Must set VARTOINTERP==PRIMTOINTERP_VSQ and RESCALEINTERP==1 in definit.h/init.h to use DOEXTRAINTERP=%d\n",DOEXTRAINTERP);
  //      myexit(77);
  //    }
  //  }


  /////
  //
  // set whether use dq
  //
  /////

  // default is not to use dq
  for(l=-NUMNEGINTERPS;l<NUMPOSINTERPS;l++){
    usedqarray[l]=0;
  }

  if(LIMADJUST){// if adjusting limiter, cannot use dq since higher order choices cannot use dq
    // nothing to do
  }
  else{ // just change those that use dq

    // 2nd order or lower can/do use dq
    usedqarray[NLIM]=1;
    usedqarray[NLIMCENT]=1;
    usedqarray[NLIMUP]=1;
    usedqarray[NLIMDOWN]=1;
    usedqarray[DONOR]=1;
    usedqarray[VANL]=1;
    usedqarray[MINM]=1;
    usedqarray[MC]=1;
  }


  return(0);
}
