

#include "decs.h"

/*! \file boundsflux.c
  \brief User Boundary conditions for fluxes

  Also calls general/frequently-used functions in bounds.tools.c

  // For fluxes, e.g. F1, assume fluxes exist everywhere -- including j/k boundary zones.  Only i-boundary zones need to be bounded.
  // This assumesCOMPZSLOOP(is,ie,js,je,ks,ke) is over boundary zones in flux.c, which in general to be compatible with any flux method (including finite volume) this is how it should be.

  // With fluxes, only need to bound each dir-flux along that direction (as presently used by ENO-type schemes)

  // Assume flux at 0 through N are computed correctly. So only need fluxes in other boundary zones.
  // Self-assigns for 0 or N for simplicity of coding

  // OUTFLOW leaves true edge of boundary unchanged
  // Therefore, if FIXEDOUTFLOW, then extrapolation is always ok.
  // if OUTFLOW, then extrapolation is ok as long as flux is from active zones out of boundary

  // GODMARK: For FLUXCTSTAG method, boundaries of values need (e.g.) F1(B2) = -F2(B1).  In reality I think this means only really need extra bounding of upper N1 or N2 or N3 boundary since lower is fine and real boundary cells are only needed along direction of flux itself.

  */

/// order of outflow extrap
/// 0: none/ copy
/// 1: first order
#define EXTRAP 0 //atch




int inboundloop[NDIM];
int outboundloop[NDIM];
int innormalloop[NDIM];
int outnormalloop[NDIM];
int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
int riin,riout,rjin,rjout,rkin,rkout;
int dosetbc[COMPDIM*2];




int bound_flux_user(int boundstage, int finalstep, SFTYPE boundtime, int boundvartype, FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int i,j,k,pl,pliter;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  int enerregion;
  int *localenerpos;
  int *doflux;


  
  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);
  //  enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
  //  localenerpos=enerposreg[enerregion];
  //  doflux=dofluxreg[enerregion];



  // periodic x1
  if ( (mycpupos[1] == 0)&&(mycpupos[1] == ncpux1 - 1) ) {
    if( (BCtype[X1DN]==PERIODIC)&&(BCtype[X1UP]==PERIODIC) ){
      // just copy from one side to another

      LOOPX1dir{

        // copy from upper side to lower boundary zones
        ri=riout;
        rj=j;
        rk=k;
        LOOPBOUND1IN PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri+1+i,rj,rk,pl); // for i=0 -> i=N1
        if(FLUXB==FLUXCTSTAG){
          if(N2>1) LOOPBOUND1IN PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri+1+i,rj,rk,pl);
          if(N3>1) LOOPBOUND1IN PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri+1+i,rj,rk,pl);
        }

        // copy from lower side to upper boundary zones
        ri=riin;
        rj=j;
        rk=k;
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri+(i-N1),rj,rk,pl); // for i=N1 -> i=0
        if(FLUXB==FLUXCTSTAG){
          if(N2>1) LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri+(i-N1),rj,rk,pl); // for i=N1 -> i=0
          if(N3>1) LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri+(i-N1),rj,rk,pl); // for i=N1 -> i=0
        }
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
        LOOPBOUND2IN PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj+1+j,rk,pl); // for j=0 -> j=N2
        if(FLUXB==FLUXCTSTAG){
          if(N1>1) LOOPBOUND2IN PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj+1+j,rk,pl); // for j=0 -> j=N2
          if(N3>1) LOOPBOUND2IN PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj+1+j,rk,pl); // for j=0 -> j=N2
        }

        // copy from lower side to upper boundary zones
        ri=i;
        rj=rjin;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj+(j-N2),rk,pl); // for j=N2 -> j=0
        if(FLUXB==FLUXCTSTAG){
          if(N1>1) LOOPBOUND2OUT PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj+(j-N2),rk,pl); // for j=N2 -> j=0
          if(N3>1) LOOPBOUND2OUT PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj+(j-N2),rk,pl); // for j=N2 -> j=0
        }
      }
    }
  }


  // periodic x3
  if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
    if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
      // just copy from one side to another

      LOOPF_12{

        // copy from upper side to lower boundary zones
        ri=i;
        rj=j;
        rk=rkout;
        LOOPBOUND3IN PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk+1+k,pl); // for k=0 -> k=N3
        if(FLUXB==FLUXCTSTAG){
          if(N1>1) LOOPBOUND3IN PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk+1+k,pl); // for k=0 -> k=N3
          if(N2>1) LOOPBOUND3IN PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk+1+k,pl); // for k=0 -> k=N3
        }

        // copy from lower side to upper boundary zones
        ri=i;
        rj=j;
        rk=rkin;
        LOOPBOUND3OUT PBOUNDLOOP(pliter,pl) MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk+(k-N3),pl); // for k=N3 -> k=0
        if(FLUXB==FLUXCTSTAG){
          if(N1>1) LOOPBOUND3OUT PBOUNDLOOP(pliter,pl) MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk+(k-N3),pl); // for k=N3 -> k=0
          if(N2>1) LOOPBOUND3OUT PBOUNDLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk+(k-N3),pl); // for k=N3 -> k=0
        }
      }
    }
  }






  ///////////////////////////
  //
  // X1 inner OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////

  // first allow extrapolation
  //
  if (mycpupos[1] == 0) {
    if(((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW))||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner r boundary condition: u, just copy */
      LOOPX1dir{
        ri=riin;
        rj=j;
        rk=k;
        LOOPBOUND1IN PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk,pl)+(MACP0A1(F1,ri+1,rj,rk,pl)-MACP0A1(F1,ri,rj,rk,pl))*(FTYPE)(i-ri);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&((BCtype[X1DN]==OUTFLOW || BCtype[X1DN]==HORIZONOUTFLOW))){ // only for RHO is it obvious what to do without primitives
            if(MACP0A1(F1,i,j,k,pl)>0.0) MACP0A1(F1,i,j,k,pl)=0.0; // GODMARK: hope this is enough to shut-down inflow (what about energy flux?)
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if(BCtype[X1DN]== ASYMM) {  //atch: added this; however unsure if any symmetry exists for fluxes
      LOOPX1dir {
        LOOPBOUND1IN PBOUNDLOOP(pliter,pl)  {
          MACP0A1(F1,i,j,k,pl) = - MACP0A1(F1,-i,j,k,pl);
          if( U1 == pl || B1 == pl ) MACP0A1(F1,i,j,k,pl) *= -1.0;  //resymm them
        }
      }
    }
    else if(BCtype[X1DN]== BCEXTRAP_VEL3 || BCtype[X1DN] == BCEXTRAP) { 
      //extrapolates fluxes with 6th order
      LOOPX1dir{
        LOOPBOUND1IN PBOUNDLOOP(pliter,pl) {
          //     MACP0A1(F1,i,j,k,pl) = interpn( 6, i,  
          //           0, MACP0A1(F1,0,j,k,pl), 
          //           1, MACP0A1(F1,1,j,k,pl), 
          //           2, MACP0A1(F1,2,j,k,pl), 
          //           3, MACP0A1(F1,3,j,k,pl),
          //           4, MACP0A1(F1,4,j,k,pl), 
          //           5, MACP0A1(F1,5,j,k,pl) );
        }
      }
    }
  }// end if mycpupos[1]==0


  ///////////////////////////
  //
  // X1 outer OUTFLOW/FIXEDOUTFLOW
  //
  ///////////////////////////


  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if(((BCtype[X1UP]==OUTFLOW || BCtype[X1UP]==HORIZONOUTFLOW))||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer r BC: outflow */

      LOOPX1dir{
        ri=riout;
        rj=j;
        rk=k;
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F1,i,j,k,pl) = MACP0A1(F1,ri,rj,rk,pl)+(MACP0A1(F1,ri,rj,rk,pl)-MACP0A1(F1,ri-1,rj,rk,pl))*(FTYPE)(i-ri);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&(BCtype[X1UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
            if(MACP0A1(F1,i,j,k,pl)<0.0) MACP0A1(F1,i,j,k,pl)=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if(BCtype[X1UP]== ASYMM) {   //atch: added this; however unsure if any symmetry exists for fluxes
      LOOPX1dir {
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl)  {
          MACP0A1(F1,i,j,k,pl) = - MACP0A1(F1,N1 + N1 - i,j,k,pl);
          if( U1 == pl || B1 == pl ) MACP0A1(F1,i,j,k,pl) *= -1.0;  //resymm them
        }
      }
    }
    else if(BCtype[X1UP]== BCEXTRAP_VEL3 || BCtype[X1UP] == BCEXTRAP) { 
      //extrapolates fluxes with 6th order
      LOOPX1dir{
        LOOPBOUND1OUT PBOUNDLOOP(pliter,pl) {
          //     MACP0A1(F1,i,j,k,pl) = interpn( 6, i, 
          //           N1,   MACP0A1(F1,N1,j,k,pl), 
          //           N1-1, MACP0A1(F1,N1-1,j,k,pl), 
          //           N1-2, MACP0A1(F1,N1-2,j,k,pl), 
          //           N1-3, MACP0A1(F1,N1-3,j,k,pl),
          //           N1-4, MACP0A1(F1,N1-4,j,k,pl), 
          //           N1-5, MACP0A1(F1,N1-5,j,k,pl) );
        }
      }
    }
  }// end if mycpu is correct


  ///////////////////////////
  //
  // X2 inner POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == 0) {
    if(((BCtype[X2DN]==OUTFLOW || BCtype[X2DN]==HORIZONOUTFLOW))||(BCtype[X2DN]==FIXEDOUTFLOW)||(BCtype[X2DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner 2 BC: outflow */

      LOOPX2dir{
        ri=i; // correct for flux
        rj=rjin;
        rk=k;
        LOOPBOUND2IN PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk,pl)+(MACP0A1(F2,ri,rj,rk,pl)-MACP0A1(F2,ri,rj-l,rk,pl))*(FTYPE)(j-rj);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&(BCtype[X2DN]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
            if(MACP0A1(F2,i,j,k,pl)<0.0) MACP0A1(F2,i,j,k,pl)=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
      LOOPX2dir{
        ri=i;
        rj=rjin;
        rk=k;
        LOOPBOUND2IN PBOUNDLOOP(pliter,pl)  MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj+(rj-j),rk,pl); // self-assigns for j=0
      }
    }
    // F2 of U2/B2 is symmetric, so no need to antisymmetrize
    // F2 antisymmetric with respect to all other quantities

    if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPX2dir{
        LOOPBOUND2IN {
          //if(POSDEFMETRIC==0){  //SASMARK need to asymmetrize the fluxes; don't understand why should not asymm. the fluxes if POSDEFMETRIC is 0
          // // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
          // // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
          //}
          //else{
          PLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl)*=-1; // anti-sym all
          MACP0A1(F2,i,j,k,U2)*=-1; // re-sym U2
          MACP0A1(F2,i,j,k,B2)*=-1; // re-sym B2
          //}
        }
      }// end loop 13
    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==0


  ///////////////////////////
  //
  // X2 outer POLARAXIS
  //
  ///////////////////////////


  if (mycpupos[2] == ncpux2-1) {
    if(((BCtype[X2UP]==OUTFLOW || BCtype[X2UP]==HORIZONOUTFLOW))||(BCtype[X2UP]==FIXEDOUTFLOW)||(BCtype[X2UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer 2 BC: outflow */

      LOOPX2dir{
        ri=i; // correct for flux
        rj=rjout;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj,rk,pl)+(MACP0A1(F2,ri,rj,rk,pl)-MACP0A1(F2,ri,rj-l,rk,pl))*(FTYPE)(j-rj);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&(BCtype[X2UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
            if(MACP0A1(F2,i,j,k,pl)>0.0) MACP0A1(F2,i,j,k,pl)=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
    else if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
      LOOPX2dir{
        ri=i;
        rj=rjout;
        rk=k;
        LOOPBOUND2OUT PBOUNDLOOP(pliter,pl)  MACP0A1(F2,i,j,k,pl) = MACP0A1(F2,ri,rj+(rj-j+2),rk,pl); // self-assigns for j=N
      }
    }

    if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

      /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
      LOOPX2dir{
        LOOPBOUND2OUT {
          //if(POSDEFMETRIC==0){  //SASMARK need to asymmetrize the fluxes; don't understand why should not asymm. the fluxes if POSDEFMETRIC is 0
          // // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
          // // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
          //}
          //else{
          PLOOP(pliter,pl) MACP0A1(F2,i,j,k,pl)*=-1; // anti-sym all
          MACP0A1(F2,i,j,k,U2)*=-1; // re-sym U2
          MACP0A1(F2,i,j,k,B2)*=-1; // re-sym B2
          //}
        }
      }// end loop 13
    } // end if POLARXIS or ASYMM
  }// end if mycpupos[2]==ncpux2-1



  //x3 inner
  if ( mycpupos[3] == 0 ) {
    if(((BCtype[X3DN]==OUTFLOW || BCtype[X3DN]==HORIZONOUTFLOW))||(BCtype[X3DN]==FIXEDOUTFLOW)||(BCtype[X3DN]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* inner 3 BC: outflow */

      LOOPX3dir{
        ri=i; // correct for flux
        rj=j;
        rk=rkin;
        LOOPBOUND3IN PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk,pl)+(MACP0A1(F3,ri,rj,rk,pl)-MACP0A1(F3,ri,rj,rk-l,pl))*(FTYPE)(k-rk);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&(BCtype[X3DN]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
            if(MACP0A1(F3,i,j,k,pl)<0.0) MACP0A1(F3,i,j,k,pl)=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
  }

  //x3 outer
  if ( mycpupos[3] == ncpux3 - 1 ) {
    if(((BCtype[X3UP]==OUTFLOW || BCtype[X3UP]==HORIZONOUTFLOW))||(BCtype[X3UP]==FIXEDOUTFLOW)||(BCtype[X3UP]==FIXED)){  //SASMARK: FIXED is not supposed to be here but need to assign sth to fluxes
      /* outer 3 BC: outflow */

      LOOPX3dir{
        ri=i; // correct for flux
        rj=j;
        rk=rkout;
        LOOPBOUND3OUT PBOUNDLOOP(pliter,pl){
#if(EXTRAP==0)
          // zeroth orde extrap
          MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk,pl);
#else
          // linear extrap
          MACP0A1(F3,i,j,k,pl) = MACP0A1(F3,ri,rj,rk,pl)+(MACP0A1(F3,ri,rj,rk,pl)-MACP0A1(F3,ri,rj,rk-l,pl))*(FTYPE)(k-rk);
#endif
          // if OUTFLOW, then disallow extrapolation if flux is from ghost to active zones
          if((pl==RHO)&&(BCtype[X3UP]==FIXEDOUTFLOW)){  //SASMARK:  changed OUTFLOW to FIXEDOUTFLOW here
            if(MACP0A1(F3,i,j,k,pl)>0.0) MACP0A1(F3,i,j,k,pl)=0.0;  //SASMARK:  can be a problem for e.g. Noh problem where there is OUTFLOW BC and the matter actually inflows
          }
        }// end pl and face loop
      }// end 2 3
    }// end if correct bound type
  }

  return (0);
}
