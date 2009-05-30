
#include "decs.h"

/* bound array containing entire set of primitive variables */

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  SFTYPE Dt=-1.0; // no boundary zone should ever be in accounting
  int i, j, k;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j




  if (mycpupos[2] == 0) {
    if((BCtype[X2DN]==OUTFLOW)){
      /* inner r boundary condition: u, just copy */
      for (i = 0; i < N1; i++) {
	ri=i;
	rj=0;
	for(j=-1;j>=-N2BND;j--)	  PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
      }
    }
  }
  
  if (mycpupos[2] == ncpux2-1) {
    if((BCtype[X2UP]==OUTFLOW)){
      /* inner r boundary condition: u, just copy */
      for (i = 0; i < N1; i++) {
	ri=i;
	rj=N2-1;
	for(j=N2;j<=N2-1+N2BND;j++)	  PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
      }
    }
  }
  
  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if((BCtype[X1UP]==OUTFLOW)){
      /* outer r BC: outflow */
      
      for (j = -N2BND; j < N2+N2BND; j++) {
	ri=N1-1;
	rj=j;
	for(i=N1;i<=N1-1+N1BND;i++)	  PBOUNDLOOP prim[i][j][k] = prim[ri][rj][k];
      }
    }
  }

  if (mycpupos[1] == 0) {
    if((BCtype[X1DN]==OUTFLOW)){
      /* inner radial BC (preserves u^t rho and u) */
      for (j = -N2BND; j <=N2-1+N2BND; j++){
	ri=0;
	rj=j;
	for(i=-N1BND;i<=-1;i++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj][k];
      }

    }

    if((BCtype[X1DN]==POLARAXIS)){
      
      /* inner radial BC (preserves u^t rho and u) */
      for (j = -N2BND; j <=N2-1+N2BND; j++){
        ri=0;
        rj=j;
        for(i=-N1BND;i<=-1;i++) PBOUNDLOOP  prim[i][j][k] = prim[ri+(ri-i-1)][rj][k];
      }

      for (j = -N2BND; j <= N2-1+N2BND; j++) {
        for (i = -N1BND; i < 0; i++) {
            prim[i][j][U1] *= -1.;
            prim[i][j][B1] *= -1.;
        }
      }
    }
  }




  if (USEMPI) bound_mpi(boundstage, prim);


  return (0);
}
