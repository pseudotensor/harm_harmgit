
#include "decs.h"

/* bound array containing entire set of primitive variables */

// GODMARK: THIS CODE IS OUTDATED, but SHOULD MAKE IT WORK SINCE NEAT (i.e. CLEAN) AND GENERAL


// should at least compile, but no idea if it works.  Definitely doesn't work in 3D

int bound_pflag(int boundstage, int prim[][N2M][N3M][NUMPFLAGS])
{
  int i, j, k, pl;
  int is,ie,js,je,ks,ke;
  int failreturn;
  int ri, rj, rk; // reference i,j,k
  int dir;
  
  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    DIRLOOP(dir){
      if((mycpupos[1]==0)&&(dir==X1DN)){
	is=-2;
	ie=-1;
	js=-2;
	je=N2+1;
      }
      else if((mycpupos[1]==ncpux1-1)&&(dir==X1UP)){
	is=N1;
	ie=N1+1;
	js=-2;
	je=N2+1;
      }
      else if((mycpupos[2]==0)&&(dir==X2DN)){
	is=-2;
	ie=N1+1;
	js=-2;
	je=-1;
      }
      else if((mycpupos[2]==ncpux2-1)&&(dir==X2UP)){
	is=-2;
	ie=N1+1;
	js=N2;
	je=N2+1;
      }
      else continue;

      if(BCtype[dir]==OUTFLOW){
	/* inner r boundary condition: u, just copy */
	ZSLOOP(is,ie,js,je,ks,ke){
	  if(dir==X1DN){
	    ri=0;
	    rj=j;
	  }
	  else if(dir==X1UP){
	    ri=N1-1;
	    rj=j;
	  } else if(dir==X2DN){
	    ri=i;
	    rj=0;
	  } else if(dir==X2UP){
	    ri=i;
	    rj=N2-1;
	  }
	  FLOOP(pl)	  prim[i][j][k][pl] = prim[ri][rj][rk][pl];
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
	
	  FLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	}
      }
    }
  }

  if (USEMPI) bound_mpi_int(boundstage, prim);


  return (0);
}

