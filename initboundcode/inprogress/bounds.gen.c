
#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: Should try to make this work since clean and general
// code should at least compile, but it's KNOWN not to work correctly.  Definitely not correct in 3D


#define EXTRAP 2
// 0: just copy
// 1: gdet or other extrapolation (no longer supported here, use old bounds.c and boundsint.c)
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][N3M][NPR])
{
  int i, j, k, pl;
  int is,ie,js,je;,ks,ke
  struct of_geom geom,rgeom;
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  int dir;
  int ref;
  
  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    DIRLOOP(dir){
      if((mycpupos[1]==0)&&(dir==X1DN)){
	is=-2;
	ie=-1;
	js=-2;
	je=N2+1;
	ref=1;
      }
      else if((mycpupos[1]==ncpux1-1)&&(dir==X1UP)){
	is=N1;
	ie=N1+1;
	js=-2;
	je=N2+1;
	ref=1;
      }
      else if((mycpupos[2]==0)&&(dir==X2DN)){
	is=-2;
	ie=N1+1;
	js=-2;
	je=-1;
	ref=2;
      }
      else if((mycpupos[2]==ncpux2-1)&&(dir==X2UP)){
	is=-2;
	ie=N1+1;
	js=N2;
	je=N2+1;
	ref=2;
      }
      else continue;

      if(BCtype[dir]==OUTFLOW){
	/* inner r boundary condition: u, just copy */
	ZSLOOP(is,ie,js,je,ks,ke){
	  if(dir==X1DN){
	    ri=0;
	    rj=j;
	  } else if(dir==X1UP){
	    ri=N1-1;
	    rj=j;
	  } else if(dir==X2DN){
	    ri=i;
	    rj=0;
	  } else if(dir==X2UP){
	    ri=i;
	    rj=N2-1;
	  }
#if(EXTRAP==0)
	  PLOOP(pl)	  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
#elif(EXTRAP==2)
	  get_geometry(ri, rj, rk, CENT, &rgeom);
	  rescale(1,ref,prim[ri][rj][rk],&rgeom,prim[ri][rj][rk]);
	  PLOOP(pl)	  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	  get_geometry(i, j, k, CENT, &geom);
	  rescale(-1,ref,prim[i][j][k],&geom,prim[i][j][k]);
	  rescale(-1,ref,prim[ri][rj][rk],&rgeom,prim[ri][rj][rk]);	
#endif
#if(WHICHVEL==VEL4)
	  get_geometry(i, j, k, CENT, &geom);
	  inflow_check_4vel(ref,prim[i][j][k],&geom) ;
#elif(WHICHVEL==VEL3)
	  get_geometry(i, j, k, CENT, &geom);
	  inflow_check_3vel(ref,prim[i][j][k],&geom) ;
	  // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
	  if(jonchecks){
	    //fixup1zone(prim[i][j][k],&geom,0);
	    failreturn=check_pr(prim[i][j][k],prim[i][j][k],&geom,-3);
	    if(failreturn){
	      dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
	      if (fail(FAIL_BCFIX) >= 1) return (1);
	    }
	  }
#elif(WHICHVEL==VELREL4)
	  get_geometry(i,j,k,CENT,&geom) ;
	  inflow_check_rel4vel(ref,prim[i][j][k],&geom) ;
	  if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom)>=1)
	    FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	}
      }
      else if(BCtype[dir]==POLARAXIS){
	ZSLOOP(is,ie,js,je,ks,ke){
	  if(dir==X1DN){
	    if(j==-1){
	      ri=0;
	      rj=j;
	    }
	    else if(j==-2){
	      ri=1;
	      rj=j;
	    }
	  }
	  else if(dir==X1UP){
	    if(j==N2){
	      ri=N1-1;
	      rj=j;
	    }
	    else if(j==N2+1){
	      ri=N1-2;
	      rj=j;
	    }
	  } else if(dir==X2DN){
	    if(i==-1){
	      ri=i;
	      rj=0;
	    }
	    else if(i==-2){
	      ri=i;
	      rj=1;
	    }
	  } else if(dir==X2UP){
	    if(i==N1){
	      ri=i;
	      rj=N2-1;
	    }
	    else if(i==N1+1){
	      ri=i;
	      rj=N2-2;
	    }
	  }
	
	  /* inner polar BC (preserves u^t rho and u) */
	  PLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	  if(POSDEFMETRIC==0){
	    /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	    // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	    // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	    prim[i][j][k][U2] *= -1.;
	    prim[i][j][k][U3] *= 1.;
	    prim[i][j][k][B2] *= -1.;
	    prim[i][j][k][B3] *= 1.;
	  }
	  else{
	    prim[i][j][k][U2] *= -1.;
	    prim[i][j][k][U3] *= 1.;
	    prim[i][j][k][B2] *= -1.;
	    prim[i][j][k][B3] *= 1.;

	  }
	}
      }
    }
  }

  if (USEMPI) bound_mpi(boundstage, prim);


  return (0);
}
