
/*! \file boundsint.c
  \brief User Boundary conditions for PFTYPE for pflag type quantities

  // this is the presently used function
  
  // Regarding GRID SECTIONING, the pflag can be set using LOOPF's not COMPLOOPF since no benefit in moving where define "real" boundary

*/

#include "decs.h"


//  needed as global so all subfunctions have it set
int inboundloop[NDIM];
int outboundloop[NDIM];
int innormalloop[NDIM];
int outnormalloop[NDIM];
int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
int riin,riout,rjin,rjout,rkin,rkout;
int dosetbc[COMPDIM*2];

/// only bounds +-1 cell as required for pflags to exist for check_solution() checking and fixup_utoprim() averaging
int bound_pflag_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, PFTYPE (*prim)[NSTORE2][NSTORE3][NUMPFLAGS])
{
  int i,j,k,pl,pliter;
  int failreturn;
  int ri, rj, rk; // reference i,j,k

  // includes all conditions from global.h except NSSURFACE and FIXED

  // OUTFLOW: direct copy here since value is discrete and has discrete meaning
  // that is, properties are just copied.

  /* if fixed BC: do nothing */
  // GODMARK: no, since doing multi-steps, one should assign some fixed primitive to prim, as done for real primitives
  // general assume flag doesn't care about FIXED, since represents a failure.  One's own FIXED conditions shouldn't fail.

  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);


  /////////////////////
  //
  // PERIODIC
  //
  /////////////////////

  // periodic x1
  if ( (mycpupos[1] == 0)&&(mycpupos[1] == ncpux1 - 1) ) {
    if( (BCtype[X1DN]==PERIODIC)&&(BCtype[X1UP]==PERIODIC) ){
      // just copy from one side to another
 
      LOOPX1dir{

        // copy from upper side to lower boundary zones
        ri=riout;
        rj=j;
        rk=k;
        LOOPBOUND1IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri+1+i,rj,rk,pl);

        // copy from lower side to upper boundary zones
        ri=riin;
        rj=j;
        rk=k;
        LOOPBOUND1OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri+(i-N1),rj,rk,pl);
      }
    }
  }

  // periodic x2
  if ( (mycpupos[2] == 0)&&(mycpupos[2] == ncpux2 - 1) ) {
    if( (BCtype[X2DN]==PERIODIC)&&(BCtype[X2UP]==PERIODIC) ){
      // just copy from one side to another
 
      LOOPX2dir{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=rjout;
        rk=k;
        LOOPBOUND2IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+1+j,rk,pl);

        // copy from lower side to upper boundary zones
        ri=i;
        rj=rjin;
        rk=k;
        LOOPBOUND2OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(j-N2),rk,pl);
      }
    }
  }

  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another
      
      LOOPX3dir{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=j;
        rk=rkout;
        LOOPBOUND3IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+1+k,pl);

        // copy from lower side to upper boundary zones
        ri=i;
        rj=j;
        rk=rkin;
        LOOPBOUND3OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(k-N3),pl);
      }
    }
  }

  /////////////////////
  //
  // OUTFLOW/FIXEDOUTFLOW (for events, just copy)
  //
  /////////////////////


  // outflow inner x1
  if (mycpupos[1] == 0) {
    if( ((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW))||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPX1dir{
        ri=riin;
        rj=j;
        rk=k;
        LOOPBOUND1IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }

  // outflow outer x1
  if (mycpupos[1] == ncpux1 - 1) {
    if( ((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW))||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPX1dir{
        ri=riout;
        rj=j;
        rk=k;
        LOOPBOUND1OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }

  // outflow inner x2
  if (mycpupos[2] == 0) {
    if( ((BCtype[X2DN]==OUTFLOW || BCtype[X2DN]==HORIZONOUTFLOW))||(BCtype[X2DN]==FIXEDOUTFLOW)||(BCtype[X2DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPX2dir{
        ri=i;
        rj=rjin;
        rk=k;
        LOOPBOUND2IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }

  // outflow outer x2
  if (mycpupos[2] == ncpux2 - 1) {
    if( ((BCtype[X2UP]==OUTFLOW || BCtype[X2UP]==HORIZONOUTFLOW))||(BCtype[X2UP]==FIXEDOUTFLOW)||(BCtype[X2UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPX2dir{
        ri=i;
        rj=rjout;
        rk=k;
        LOOPBOUND2OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }

  // outflow inner x3
  if (mycpupos[3] == 0) {
    if( ((BCtype[X3DN]==OUTFLOW || BCtype[X3DN]==HORIZONOUTFLOW))||(BCtype[X3DN]==FIXEDOUTFLOW)||(BCtype[X3DN]==OUTFLOWNOINFLOW) ){
      /* inner r boundary condition: u, just copy */
      LOOPX3dir{
        ri=i;
        rj=j;
        rk=rkin;
        LOOPBOUND3IN FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }

  // outflow outer x3
  if (mycpupos[3] == ncpux3 - 1) {
    if( ((BCtype[X3UP]==OUTFLOW || BCtype[X3UP]==HORIZONOUTFLOW))||(BCtype[X3UP]==FIXEDOUTFLOW)||(BCtype[X3UP]==OUTFLOWNOINFLOW) ){
      /* outer r BC: outflow */
      LOOPX3dir{
        ri=i;
        rj=j;
        rk=rkout;
        LOOPBOUND3OUT FBOUNDLOOP(pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
      }
    }
  }


  /////////////////////
  //
  // POLARAXIS/SYMM/ASYMM (for events (not values) these are the same)
  //
  /////////////////////


  // symmetry on inner x1
  if (mycpupos[1] == 0) {
    if( (BCtype[X1DN]==R0SING)||(BCtype[X1DN]==POLARAXIS)||(BCtype[X1DN]==SYMM)||(BCtype[X1DN]==ASYMM)) {
      LOOPX1dir{
        ri=riin;
        rj=j;
        rk=k;
        // symmetric copy
        LOOPBOUND1IN FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri+(ri-i-1),rj,rk,pl);
      }
    }
  }

  // symmetry on outer x1
  if (mycpupos[1] == ncpux1 - 1) {
    if( (BCtype[X1UP]==R0SING)||(BCtype[X1UP]==POLARAXIS)||(BCtype[X1UP]==SYMM)||(BCtype[X1UP]==ASYMM)) {
      LOOPX1dir{
        ri=riout;
        rj=j;
        rk=k;
        LOOPBOUND1OUT FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri+(ri-i+1),rj,rk,pl);
      }
    }
  }


  // symmetry on inner x2
  if (mycpupos[2] == 0) {
    if( (BCtype[X2DN]==R0SING)||(BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM)) {
      LOOPX2dir{
        ri=i;
        rj=rjin;
        rk=k;
        // symmetric copy
        LOOPBOUND2IN FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j-1),rk,pl);
      }
    }
  }

  // symmetry on outer x2
  if (mycpupos[2] == ncpux2 - 1) {
    if( (BCtype[X2UP]==R0SING)||(BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM)) {
      LOOPX2dir{
        ri=i;
        rj=rjout;
        rk=k;
        LOOPBOUND2OUT FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj+(rj-j+1),rk,pl);
      }
    }
  }

  // symmetry on inner x3
  if (mycpupos[3] == 0) {
    if( (BCtype[X3DN]==R0SING)||(BCtype[X3DN]==POLARAXIS)||(BCtype[X3DN]==SYMM)||(BCtype[X3DN]==ASYMM)) {
      LOOPX3dir{
        ri=i;
        rj=j;
        rk=rkin;
        // symmetric copy
        LOOPBOUND3IN FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(rk-k-1),pl);
      }
    }
  }

  // symmetry on outer x3
  if (mycpupos[3] == ncpux3 - 1) {
    if( (BCtype[X3UP]==R0SING)||(BCtype[X3UP]==POLARAXIS)||(BCtype[X3UP]==SYMM)||(BCtype[X3UP]==ASYMM)) {
      LOOPX3dir{
        ri=i;
        rj=j;
        rk=rkout;
        LOOPBOUND3OUT FBOUNDLOOP(pl)  MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk+(rk-k+1),pl);
      }
    }
  }



  return (0);
}
