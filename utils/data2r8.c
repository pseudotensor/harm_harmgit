#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  float outputf;
  unsigned char output;
  float input;
  float min,max;
  int i;
  int snx,sny;

  if(argc<4){
    fprintf(stderr,"not enough arguments\n");
    fprintf(stderr,"e.g. data2r8 25 1950 min max < intputfile > outputfile\n");
    exit(1);
  }

  snx=atoi(argv[1]);
  sny=atoi(argv[2]);
  
  sscanf(argv[3],"%f",&min);
  sscanf(argv[4],"%f",&max);
 
  fprintf(stderr,"snx=%d sny=%d min=%g max=%g\n",snx,sny,min,max); fflush(stderr);

  for(i=0;i<snx*sny;i++){
    fscanf(stdin,"%f",&input);
    //    fprintf(stderr,"input=%f\n",input);
    outputf=(input-min)/(max-min)*255.0;
    output=(unsigned char)(outputf);
    fputc(output,stdout);
    //    fprintf(stderr,"%f -> %c\n",input,output); fflush(stderr);
  }
  return(0);

}
