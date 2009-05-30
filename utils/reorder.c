// gcc -o reorder reorder.c

#include <stdio.h>

#define NUMCOL 13
#define NX 300
#define NY 300

float array[NX][NY][NUMCOL];

int main(void)
{
  FILE *in;
  FILE *out;
  int i,j,col;
  char ch;


  in=fopen("jondisk.dat","r");
  out=fopen("jondiskreordered.dat","w");

  // skip 5 header lines and pipe them to output file
  for(i=1;i<=5;i++){
    //    while(fgetc(in)!='\n');
    while(1){
      ch=fgetc(in);
      if(ch=='\n'){
	fputc(ch,out);
	break;
      }
      else fputc(ch,out);
    }
  }

  fprintf(out,"\n\n");

  // i is associated with fastest quantity, which is R

  for(j=0;j<NY;j++){
    for(i=0;i<NX;i++){
      for(col=0;col<NUMCOL;col++){
	fscanf(in,"%f",&array[i][j][col]);
	//fprintf(stderr,"i=%d j=%d col=%d array=%21.15g\n",i,j,col,array[i][j][col]);
      }
    }
  }

  fclose(in);

  // i make to be slowest quantity, but still attached to R variable

  for(i=-1;i<NX;i++){
    for(j=0;j<NY;j++){
      for(col=0;col<NUMCOL;col++){
	if((i==-1)&&(col==0)) fprintf(out,"%15.7g ",-array[i+1][j][col]);
	else if(i==-1) fprintf(out,"%15.7g ",array[i+1][j][col]);
	else fprintf(out,"%15.7g ",array[i][j][col]);
      }
      fprintf(out,"\n");
    }
  }

  fclose(out);




}
