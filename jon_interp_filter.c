#include "decs.h"


// doesn't process t or z directions
// GODMARK3D
void gaussian_filter(int filter,FTYPE sigma, int nt, int nx, int ny, int nz, unsigned char*****oldimage,FTYPE*****olddata)
{
  int i,j;
  int ii,jj;
  FTYPE S,r2,a;
  FTYPE **W,**fmatrix(int a, int b, int c, int d)  ;
  FTYPE **ftemp;
  FTYPE ftemp2;
  unsigned char **ctemp,**cmatrix(int a, int b, int c, int d)  ;
  
  ftemp = fmatrix(0,nx-1,0,ny-1);  

  W = fmatrix(-filter,filter,-filter,filter) ;

  // general weight factors
  S=0;
  r2=2.0*sigma*sigma;
  for(j=-filter;j<=filter;j++)  for(i=-filter;i<=filter;i++){
      a=((FTYPE)i*(FTYPE)i+(FTYPE)j*(FTYPE)j)/r2;
      W[i][j]=exp(-a);
      S+=W[i][j];
    }

  // general image filtering
  for(j=filter;j<ny-filter;j++) for(i=filter;i<nx-filter;i++){
      ftemp2=0.0;
      for(jj=-filter;jj<=filter;jj++)  for(ii=-filter;ii<=filter;ii++){
          if(DATATYPE==0) ftemp2+=W[ii][jj]*(FTYPE)oldimage[0][0][i+ii][j+jj][0]; // GODMARK3D
          else ftemp2+=W[ii][jj]*olddata[0][0][i+ii][j+jj][0]; // GODMARK3D
        }
      ftemp[i][j]=ftemp2/S;
    }
  for(j=filter;j<ny-filter;j++) for(i=filter;i<nx-filter;i++){
      if(DATATYPE==0) oldimage[0][0][i][j][0]=(unsigned char)ftemp[i][j]; // only overwrites filtered parts, assumes rest still there // GODMARK3D
      else olddata[0][0][i][j][0]=ftemp[i][j]; // GODMARK3D
    }

  free_fmatrix(ftemp,0,nx-1,0,ny-1);
  free_fmatrix(W,-filter,filter,-filter,filter);

}

