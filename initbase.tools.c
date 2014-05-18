
/*! \file initbase.tools.c
     \brief Functions that don't depend upon global quantities/functions
*/

#include <stdio.h>
#include <stdlib.h>


/// report system info
void report_systeminfo(FILE * fileout)
{
  int IsLittleEndian(void);
  //////////////////////
  //
  // check endianness of system and MPI compatibility

  if(IsLittleEndian()){
    fprintf(fileout,"System is Little Endian\n");
  }
  else fprintf(fileout,"System is Big Endian\n");

  // inform of data sizes
  fprintf(fileout,"Size of char: %lld\n",(long long int)sizeof(char));
  fprintf(fileout,"Size of short: %lld\n",(long long int)sizeof(short));
  fprintf(fileout,"Size of int: %lld\n",(long long int)sizeof(int));
  fprintf(fileout,"Size of size_t: %lld\n",(long long int)sizeof(size_t));
  fprintf(fileout,"Size of long int: %lld\n",(long long int)sizeof(long int));
  fprintf(fileout,"Size of long long int: %lld\n",(long long int)sizeof(long long int));
  fprintf(fileout,"Size of float: %lld\n",(long long int)sizeof(float));
  fprintf(fileout,"Size of double: %lld\n",(long long int)sizeof(double));
  fprintf(fileout,"Size of long double: %lld\n",(long long int)sizeof(long double));

  if((long long int)sizeof(int)<=4){
    fprintf(fileout,"WARNING: Since size of integer is only 4 bytes, some routines that input and output integer arguments will be limited to 2GB.  For example, ROMIO will fail to work if try to write >2GB file, since buffer sizes must then be <2GB\n");
  }


  ///////////////////
  //
  // Some checks
  //
  ///////////////////

  if(sizeof(char)!=1){
    fprintf(fileout,"sizeof(char) was not 1 byte, void pointer use not going to be correct\n");
    exit(1);
  }

}


// OPENMPMARK: constant static, so ok
static long _TestEndian=1;

int IsLittleEndian(void)
{
  return *(char*)&_TestEndian;
}

///
///  FUNCTION: SwapEndian
///  PURPOSE: Swap the byte order of a structure
///  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
void *SwapEndian(void* Addr, const int Nb)
{
  static char Swapped[16];
  switch (Nb) {
  case 2: Swapped[0]=*((char*)Addr+1);
    Swapped[1]=*((char*)Addr  );
    break;
  case 3: // As far as I know, 3 is used only with RGB images
    Swapped[0]=*((char*)Addr+2);
    Swapped[1]=*((char*)Addr+1);
    Swapped[2]=*((char*)Addr  );
    break;
  case 4: Swapped[0]=*((char*)Addr+3);
    Swapped[1]=*((char*)Addr+2);
    Swapped[2]=*((char*)Addr+1);
    Swapped[3]=*((char*)Addr  );
    break;
  case 8: Swapped[0]=*((char*)Addr+7);
    Swapped[1]=*((char*)Addr+6);
    Swapped[2]=*((char*)Addr+5);
    Swapped[3]=*((char*)Addr+4);
    Swapped[4]=*((char*)Addr+3);
    Swapped[5]=*((char*)Addr+2);
    Swapped[6]=*((char*)Addr+1);
    Swapped[7]=*((char*)Addr  );
    break;
  case 16:Swapped[0]=*((char*)Addr+15);
    Swapped[1]=*((char*)Addr+14);
    Swapped[2]=*((char*)Addr+13);
    Swapped[3]=*((char*)Addr+12);
    Swapped[4]=*((char*)Addr+11);
    Swapped[5]=*((char*)Addr+10);
    Swapped[6]=*((char*)Addr+9);
    Swapped[7]=*((char*)Addr+8);
    Swapped[8]=*((char*)Addr+7);
    Swapped[9]=*((char*)Addr+6);
    Swapped[10]=*((char*)Addr+5);
    Swapped[11]=*((char*)Addr+4);
    Swapped[12]=*((char*)Addr+3);
    Swapped[13]=*((char*)Addr+2);
    Swapped[14]=*((char*)Addr+1);
    Swapped[15]=*((char*)Addr  );
    break;
  }
  return (void*)Swapped;
}

