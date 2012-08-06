#include <stdio.h>
#include <math.h>

int main()
{
  double x,y,z;

  x=0;  y=1;  z=atan2(y,x);  fprintf(stderr,"x=%21.15g y=%21.15g z=%21.15g\n",x,y,z); fflush(stderr);

  x=1;  y=0;  z=atan2(y,x);  fprintf(stderr,"x=%21.15g y=%21.15g z=%21.15g\n",x,y,z); fflush(stderr);

  x=-1;  y=0;  z=atan2(y,x);  fprintf(stderr,"x=%21.15g y=%21.15g z=%21.15g\n",x,y,z); fflush(stderr);

  x=0;  y=-1;  z=atan2(y,x); if(z<0.0) z+=2.0*M_PI; if(z>2.0*M_PI) z-=2.0*M_PI;  fprintf(stderr,"x=%21.15g y=%21.15g z=%21.15g\n",x,y,z); fflush(stderr);

  x=.5;  y=-1;  z=atan2(y,x); if(z<0.0) z+=2.0*M_PI; if(z>2.0*M_PI) z-=2.0*M_PI;  fprintf(stderr,"x=%21.15g y=%21.15g z=%21.15g\n",x,y,z); fflush(stderr);

}
