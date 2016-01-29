#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/////////////////
//
// chooseable:
//
/////////////////
#define HDF (1)
#define V5D (1)





/////////////////
//
// rest NOT chooseable unless you know what you are doing:
//
/////////////////
#ifndef _SWAP_ENDIAN
#define _SWAP_ENDIAN

// Macs and SGIs are Big-Endian; PCs are little endian
// returns TRUE if current machine is little endian
static int IsLittleEndian(void);

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

#define SWAP_SHORT(Var)  Var = *(short*)         SwapEndian((void*)&Var, sizeof(short))
#define SWAP_USHORT(Var) Var = *(unsigned short*)SwapEndian((void*)&Var, sizeof(short))
#define SWAP_LONG(Var)   Var = *(long*)          SwapEndian((void*)&Var, sizeof(long))
#define SWAP_ULONG(Var)  Var = *(unsigned long*) SwapEndian((void*)&Var, sizeof(long))
#define SWAP_RGB(Var)    Var = *(int*)           SwapEndian((void*)&Var, 3)
#define SWAP_FLOAT(Var)  Var = *(float*)         SwapEndian((void*)&Var, sizeof(float))
#define SWAP_DOUBLE(Var) Var = *(double*)        SwapEndian((void*)&Var, sizeof(double))

static void *SwapEndian(void* Addr, const int Nb);

#endif


#define LITTLE_ENDIAN_USER 0
#define BIG_ENDIAN_USER    1


// this program does not convert between data types, only between file formats
// for HDF we have to choose what datatype is read/written, and it's fixed for all output/input types

#if(HDF)
// note, for HDF, dim0,1,2,3,... means the array element array[dim0][dim1][dim2][dim3], which means dimN where N is the dimensionality of the vector, is most quickly iterated.
// so a loop like for(k) for(j) for(i) would be used on array[k][j][i]
// so i is actually dim2, j is dim1, and k is dim0

// SOMETIMES SYSTEM USER OR INSTALLATION TYPE (e.g. ki-rh42,physics-179.umd.edu)
#include "hdf/hdf.h"
#include "hdf/mfhdf.h"

// SOMETIMES USER INSTALLATION TYPE (e.g. ki-rh39)
//#include "hdf.h"
//#include "mfhdf.h"

#define HDFTYPE DFNT_FLOAT32
//#define HDFTYPE DFNT_FLOAT64
// could be DFNT_BYTE DFNT_INT32 or DFNT_FLOAT64 for doubles
#define SDS_NAME      "Jonathan McKinney SDS"
#endif


#if(V5D)
#include <vis5d+/v5d.h>
#include <vis5d+/binio.h>
#endif

// note that for binary, should be concerned with big or little endianness
#define BYTE unsigned char
static void swap(BYTE *x, BYTE size);
#define swaporder16(A)  ( (((A) & 0xff00) >> 8) | (((A) & 0x00ff) << 8) )
#define swaporder32(A)  ( (((A) & 0xff000000) >> 24) | (((A) & 0x00ff0000) >> 8) | (((A) & 0x0000ff00) << 8) | (((A) & 0x000000ff) << 24))

#define DEBUG (0)


// debug compile                       gcc -O0 -g -Wall -o bin2txt bin2txt.c -lmfhdf -ldf -ljpeg -lz  -lv5d

// standard compile on ki-rh42:
// gcc -O2 -Wall -o bin2txt bin2txt.c -I /usr/include/hdf/ -L /usr/lib64/hdf/ -lmfhdf -ldf -ljpeg -lz -lv5d
 
// standard compile on ki-rh39:
// gcc -O2 -Wall -o bin2txt bin2txt.c -I /usr2/local/include/ -L /usr2/local/lib/ -lmfhdf -ldf -ljpeg -lz -lv5d

// JCM 10-4-01

// JCM 10-30-01 added hdf4 capability

#define BUFFERMAP ((k*N2*N1+j*N1+i)*numcolumns+nextbuf++)
#define BUFFERINIT nextbuf=0



#define MAXTERMS 100

int main(
         int argc,
         char *argv[],
         char *envp[]
         )
{
  int machineEndianness(void);
  int argstep=0;
  int NDIM;
  long long int i,j,k; // just so products will stay long long int without cast
  int numrows,numcolumns,numterms,numheader,pipeheader,realnumheader;
  int bytesize,intsize,floatsize,doublesize;
  unsigned char dumb;
  unsigned short shortfordumb;
  int dumi;
  float dumf;
  double dumlf;
  char precision[MAXTERMS];
  int group[MAXTERMS];

  long long int N1,N2,N3,TS; // just so products will stay long long int without cast
  int inputlist,outputlist,input4dblock;

  long long int nextbuf;
  FILE * input;
  FILE * output;
  FILE * inputfilenamelist;
  FILE * outputfilenamelist;
  char INPUT_NAME[400];
  char OUTPUT_NAME[400];
  char INPUTLIST_NAME[400];
  char OUTPUTLIST_NAME[400];
  char V5DHEAD_NAME[400];
  int source,dest,ble;

  short *arrayb;
  int *arrayi;
  float *arrayf;
  double *arrayd;
  void * array;
  float *arrayvisf;
  void * arrayvis;
  float *arrayvisoutput;

#if(HDF)

  // hdf
  int32 sd_id, sds_id, sd_index,istat;
  intn  status;
  int32 rank;
  int32 dim_sizes[4], start[4], edges[4];
#endif
  int it;
#if(V5D)
  // vis5d
  int iv;
  int NumTimes;                      /* number of time steps */
  int NumVars;                       /* number of variables */
  int Nr, Nc, Nl[MAXVARS];           /* size of 3-D grids */
  char VarName[MAXVARS][MAXVARNAME];         /* names of variables */ 
  int TimeStamp[MAXTIMES];           /* real times for each time step */
  int DateStamp[MAXTIMES];           /* real dates for each time step */
  int CompressMode;                  /* number of bytes per grid */
  int Projection;                    /* a projection number */
  float ProjArgs[MAXPROJARGS];       /* the projection parameters */
  int Vertical;                      /* a vertical coord system number */
  float VertArgs[MAXVERTARGS];         /* the vertical coord sys parameters */
  FILE * vis5dheader;
  float pos[3+1][2+1];
  float minmax[2][MAXVARS];
  float a,b;
  long long int ijkvis5d,ijkjon;
#endif
  

  // r8 stuff
  char incommand[300],outcommand[300];
  int lname;
  int inputtype,outputtype;
  char ch;

  // begin

  bytesize=sizeof(unsigned char);
  intsize=sizeof(int);
  floatsize=sizeof(float);
  doublesize=sizeof(double);


  if(argc<11){
    fprintf(stderr,"non-V5D Usage: bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS INPUTNAME OUTPUTNAME <format>\n");
    fprintf(stderr," SOURCE/DEST: 0=r8 (.gz) 1=binary 2=text 3=HDF4 4=HDF5 5=vis5d\n");
    fprintf(stderr," BLE: same endian (0), big to little or little to big (same process) (1), auto-conversion assuming ROMIO binary (-1)\n");
    fprintf(stderr," HEADERLINESTOSKIP: # of header text lines in r8 or text source or binary header (negative values mean pipe text header into new file)\n");
    fprintf(stderr," NDIM/N1/N2/N3/TS: # of dimensions and elements in each dimension (e.g. if 2D, N3's value doesn't matter but needs to be there) and # of timesteps to include TS must be >0\n");
    fprintf(stderr," INPUTNAME/OUTPUTNAME: input and output file names (not used if TS>1)\n");
    fprintf(stderr," <format> type is (for 1-4) b=byte i=integer, f=float, d=double, where you give <type1> #of1 <type2> #of2 ...\n");
    fprintf(stderr,"\ne.g. bin2txt 1 2 0 3 32 32 32 1 data32-3.txt data32-3.bin i 3 f 9\n\n");

    fprintf(stderr,"\n");

    fprintf(stderr,"x=0,1,2 && y=0,1,2 w/ TS>1 then also include:\n");
    fprintf(stderr,"Usage (dest=0,1,2 w/ TS=+|TS|): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS INPUTLIST OUTPUTNAME <format>\n");
    fprintf(stderr,"Usage (source=0,1,2 w/ TS=-|TS|): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS INPUTNAME OUTPUTLIST <format>\n");
    fprintf(stderr,"INTPUTLIST contains: sequence of return delimited input file names for each timeslice\n");
    fprintf(stderr,"OUTPUTLIST contains: sequence of return delimited output file names for each timeslice\n");
    fprintf(stderr,"e.g., to input list and output single 4D file: bin2txtn x y 0 -1 3 256 128 32 +67 dumplist.txt dump.??? d 6\n");
    fprintf(stderr,"e.g., to input 4D file and output many files: bin2txtn x y 0 -1 3 256 128 32 -67 dump4d.??? dumplist.txt d 6\n");

    fprintf(stderr,"\n");

    fprintf(stderr,"Vis5D with TS=1:\n");
    fprintf(stderr,"Usage (dest=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTNAME OUTPUTNAME <format>\n");
    fprintf(stderr,"HEADFILE contains return delimited list of names for all variables and their min and max values (min/max only used if source=0) (names MUST BE <=9 characters!!!)\n");
    fprintf(stderr,"Last line must have:\nx_start x_finish y_start y_finish z_start z_finish\n");
    fprintf(stderr,"e.g. for one variable:\n\ndensity 1E-4 1\n");
    fprintf(stderr,"0 1 0 1 0 1\n");
    fprintf(stderr,"e.g. for fieldline file with 11 columns:\n");
    fprintf(stderr,"density 1E-4 1\n"
            "ug 1E-4 1\n"
            "negudt 1E-4 1\n"
            "mu 1E-4 1\n"
            "uut 1E-4 1\n"
            "vr 1E-4 1\n"
            "vh 1E-4 1\n"
            "vph 1E-4 1\n"
            "Br 1E-4 1\n"
            "Bh 1E-4 1\n"
            "Bp 1E-4 1\n"
            "0 1 0 1 0 1\n"
            );
    fprintf(stderr,"vis5d e.g.: ./bin2txt 0 5 0 3 64 64 64 1 vis5d.image.head imx0-0-0-s1-0000.dat.r8.gz imx0-0-0-s1-0000.dat.v5d b 1\n");

    fprintf(stderr,"\n");

    fprintf(stderr,"V5D w/ TS>1 then also include:\n");
    fprintf(stderr,"Usage (dest=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTLIST OUTPUTNAME <format>\n");
    fprintf(stderr,"Usage (source=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTNAME OUTPUTLIST <format>\n");
    fprintf(stderr,"INTPUTLIST contains: sequence of return delimited input file names for each timeslice\n");
    fprintf(stderr,"OUTPUTLIST contains: sequence of return delimited output file names for each timeslice\n");
    fprintf(stderr,"vis5d: bin2txtn 2 5 0 7 3 32 32 32 101 vis5d.dump.head dumplist.txt dump.v5d d 9\n");
    fprintf(stderr,"vis5d: bin2txtn 2 5 0 0 3 32 32 32 101 vis5d.force.head forcelist.txt forcex.v5d d 8\n");

    exit(1);
  }



  //////////////////
  //
  // BEGIN Read-in arguments
  //
  //////////////////


  // argv[0] is filename of program
  argstep=1;
  source=atoi(argv[argstep++]);
  dest=atoi(argv[argstep++]);

  ble=atoi(argv[argstep++]);
  if(ble==-1){
    if(machineEndianness()==LITTLE_ENDIAN_USER){
      ble=1; // need to convert
    }
    else{
      ble=0; // no need to convert
    }
  }
  // else use user-inputted ble
  if(ble!=0 && ble!=1){
    fprintf(stderr,"ble=%d undefined\n",ble);
    exit(1);
  }

  numheader=atoi(argv[argstep++]);
  
  if(numheader<0){
    pipeheader=1;
    realnumheader=-numheader;
  }
  else{
    pipeheader=0;
    realnumheader=numheader;
  }
  NDIM=atoi(argv[argstep++]);
  N1=atoi(argv[argstep++]);
  N2=atoi(argv[argstep++]);
  N3=atoi(argv[argstep++]);
  TS=atoi(argv[argstep++]);
  if(TS==0){
    TS=1;
    fprintf(stderr,"Assumed by TS=0 you meant TS=%d really\n",TS);
  }
  if(TS<0){
    TS=-TS;
    fprintf(stderr,"Assumed by TS<0 you meant TS=%d but inputting 4D data block\n",TS);
    input4dblock=1;
  }
  else{
    input4dblock=0;
    fprintf(stderr,"outputting 4D data block\n");
  }
  fprintf(stderr,"source: %d dest: %d numheader: %d NDIM: %d N1: %d N2: %d N3: %d TS: %d\n",source,dest,numheader,NDIM,N1,N2,N3,TS);
  // modify N1,N2,N3 for various dimension
  if(NDIM==2){
    N3=1; // force
  }
  else if(NDIM==1){
    N2=N3=1; // force
  }
  else if(NDIM==4){
    fprintf(stderr,"Not yet setup for direct 4D input\n");
    exit(1);
  }

  // set number of rows
  numrows=N1*N2*N3;



  ////////////////////////
  //
  // get input-output file (or list)
  //
  ////////////////////////
  if((TS>1)&&(dest==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"V5DHEAD_NAME: %s INPUTLIST_NAME: %s  OUTPUT_NAME: %s\n",V5DHEAD_NAME,INPUTLIST_NAME,OUTPUT_NAME); fflush(stderr);
  }
  else if((TS>1)&&(source==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    fprintf(stderr,"V5DHEAD_NAME: %s INPUT_NAME: %s OUTPUTLIST: %s \n",V5DHEAD_NAME,INPUT_NAME,OUTPUTLIST_NAME); fflush(stderr);
  }
  else if((dest==5)||(source==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"V5DHEAD_NAME: %s INPUT_NAME: %s OUTPUT_NAME: %s \n",V5DHEAD_NAME,INPUT_NAME,OUTPUT_NAME); fflush(stderr);
  }
  else if((TS>1 && input4dblock==0)&&(dest==0 || dest==1 || dest==2)){
    strcpy(INPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"INPUTLIST_NAME: %s  OUTPUT_NAME: %s\n",INPUTLIST_NAME,OUTPUT_NAME); fflush(stderr);
  }
  else if((TS>1 && input4dblock==1)&&(source==0 || source==1 || source==2)){
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    fprintf(stderr,"INPUT_NAME: %s OUTPUTLIST: %s \n",INPUT_NAME,OUTPUTLIST_NAME); fflush(stderr);
  }
  else{ // then don't need v5d header and don't need multi-file list
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"INPUT_NAME: %s OUTPUT_NAME: %s\n",INPUT_NAME,OUTPUT_NAME); fflush(stderr);
  }



  //////////////////
  //
  // Setup format specifier from command line
  //
  //////////////////
  j=0;
  numcolumns=0;
  while(argstep<argc){  
    sscanf(argv[argstep++],"%c",&precision[j]);
    fprintf(stderr,"precision[%d]=%c\n",j,precision[j]);
    sscanf(argv[argstep++],"%d",&group[j]);
    fprintf(stderr,"group[%d]=%d\n",j,group[j]);
    numcolumns+=group[j];
    j++;
    if(j>MAXTERMS){
      fprintf(stderr,"too many terms!\n");
      exit(1);
    }
  }

  numterms=j;


  //////////////////
  //
  // END Read-in arguments
  //
  //////////////////




  fprintf(stderr,"numrows: %d numcolumns: %d numterms: %d\n",numrows,numcolumns,numterms); fflush(stderr);




  /////////////////////////////
  //
  // vis5d stuff
  //
  /////////////////////////////
#if(V5D)

  // for source==5, don't actually use header, so header can be blank but file name should still be added to command line

  if(dest==5){

    NumTimes=TS;
    NumVars=numcolumns;
    //    //    Nr=N3; // (whine) Avery said use N3 instead of N1
    //    Nr=N2; // (whine) Avery said use N3 instead of N1 // trying
    //    //    Nc=N2;  
    //    Nc=N3; // trying
    Nr=N3;
    Nc=N2;

    fprintf(stderr,"NumTimes=%d  NumVars=%d Nr=%d Nc=%d\n",NumTimes,NumVars,Nr,Nc); fflush(stderr);

    for(i=0;i<NumVars;i++){    
      Nl[i]=N1; // all variables have same # of N3 elements
      // (whine) Avery said use N1 instead of N3
      fprintf(stderr,"Nl[%d]=%d\n",i,Nl[i]);fflush(stderr);
    }

    // read special header file to setup vis5d
    if( (vis5dheader=fopen(V5DHEAD_NAME,"rt"))==NULL){
      fprintf(stderr,"can't open vis5d header file %s\n","vis5d.head");
      exit(1);
    }
    // header has in it:
    // list of names for all variables types in sequence (MUST BE <=9 characters!!!)
    // list of min max for each of the variables
    // x_start x_finish y_start y_finish z_start z_finish
    for(i=0;i<NumVars;i++){
      fscanf(vis5dheader,"%s %f %f",VarName[i],&minmax[0][i],&minmax[1][i]);
      fprintf(stderr,"VarName[%d of %d]=%s min: %g max: %g\n",i,NumVars-1,VarName[i],minmax[0][i],minmax[1][i]);  fflush(stderr);
    }    

    fprintf(stderr,"DONE fscanf1\n"); fflush(stderr);

    fscanf(vis5dheader,"%f",&pos[1][1]);
    fscanf(vis5dheader,"%f",&pos[1][2]);
    fscanf(vis5dheader,"%f",&pos[2][1]);
    fscanf(vis5dheader,"%f",&pos[2][2]);
    fscanf(vis5dheader,"%f",&pos[3][1]);
    fscanf(vis5dheader,"%f",&pos[3][2]);

    fprintf(stderr,"DONE fscanf2: %f %f %f %f %f %f\n",pos[1][1],pos[1][2],pos[2][1],pos[2][2],pos[3][1],pos[3][2]); fflush(stderr);

    while(fgetc(vis5dheader)!='\n'); // skip rest of line

    fprintf(stderr,"DONE fgetc\n"); fflush(stderr);

    for(i=0;i<NumTimes;i++){
      TimeStamp[i]=i%60+((i/60)%60)*100+(i/(60*60))*10000;
      //fprintf(stderr,"ts: %06d\n",TimeStamp[i]); fflush(stderr);
      DateStamp[i]=99036;
    }

  }


  fprintf(stderr,"DONE fscanfs\n"); fflush(stderr);


  if((dest==5)||(source==5)){

    if(source==0 || dest==0) CompressMode=1; // 1,2,4 bytes per grid point
    else CompressMode=4; // 1,2,4 bytes per grid point (assume want more precision if data)
    // in general, CompressMode==2 can be bad for vector tracing when vectors vary alot away from (say) central BH
    //    CompressMode=1; // override

    Projection=0; // 0=linear, rectangular, generic units 1=linear, rectangular,cylindrical-equidistant,2=Lambert Conformal, 3=Stereographic, 4=Rotated
    // For Projection=0, comments below are true for each ProjArgs[?]
    ProjArgs[0]=pos[1][2]; // North bondary of 3D box
    ProjArgs[1]=pos[2][2]; // West boundary of 3D box
    ProjArgs[2]=(pos[1][2]-pos[1][1])/(float)(Nr-1); // increment between rows
    ProjArgs[3]=(pos[2][2]-pos[2][1])/(float)(Nc-1); // increment between columns
    Vertical=0; // 0=equally spaced in generic units 1=equally spaced in km 2=unequally spaced in km 3=unequally spaced in mb
    VertArgs[0]=pos[3][1]; // height of bottom level
    VertArgs[1]=(pos[3][2]-pos[3][1])/(float)(Nl[0]-1); // spacing between levels

    // if use Vertical==3, then do something like:
    //    for (il=0; il< Nl[il]; il++) VertArgs[il] = 1000.0 - 40.0 * il;
    // see convert/test.c in source code
    
    fprintf(stderr,"Args: %g %g : %g %g:  %g %g\n",ProjArgs[0],ProjArgs[1],ProjArgs[2],ProjArgs[3],VertArgs[0],VertArgs[1]); fflush(stderr);

    if(dest==5){
      fclose(vis5dheader); // no longer needed
    }

  }

  fprintf(stderr,"DONE ProJVert\n"); fflush(stderr);


#endif



  /////////////////
  //
  // see if TS>1 and get filename list
  //
  // after done, input and output names will be formed from these lists if inputlist or outputlist are non-zero
  ////////////////

  outputlist=inputlist=0; //default

  if(TS>1){
    // then ignore argument file input/output names and get from list
    // order is ALL input names in 1 file, ALL output names in another file
    if((dest==5 && source!=5) || (input4dblock==0)&&(dest==0 || dest==1 || dest==2 || dest==3 || dest==4)){ // otherwise 1 file
      fprintf(stderr,"Opening %s\n",INPUTLIST_NAME);
      if( (inputfilenamelist=fopen(INPUTLIST_NAME,"rt"))==NULL){
        fprintf(stderr,"can't open input filenamelist header file %s\n",INPUTLIST_NAME);
        exit(1);
      }
      inputlist=1;
      outputlist=0;
    }
    else if((dest!=5 && source==5) || (input4dblock==1)&&(source==0 || source==1 || source==2)){ // otherwise 1 file (can't source 4D from HDF yet)
      fprintf(stderr,"Opening %s\n",OUTPUTLIST_NAME);
      if( (outputfilenamelist=fopen(OUTPUTLIST_NAME,"rt"))==NULL){
        fprintf(stderr,"can't open output filenamelist header file %s\n",OUTPUTLIST_NAME);
        exit(1);
      }
      inputlist=0;
      outputlist=1;
    }
    else{
      fprintf(stderr,"Not setup to process TS=%d>1 with source=%d and dest=%d\n",TS,source,dest);
      exit(1);
    }

  }






  /////////////////////
  //
  // HDF stuff
  //
  /////////////////////
#if(HDF)
  if((source==3)||(dest==3)||(source==4)||(dest==4)){ // use maximal holder array (doubles)
    if(HDFTYPE==DFNT_UINT8){
      arrayi=(int*)malloc(sizeof(int)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayi;
    }
    if(HDFTYPE==DFNT_INT32){
      arrayi=(int32*)malloc(sizeof(int32)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayi;
    }
    if(HDFTYPE==DFNT_FLOAT32){
      arrayf=(float32*)malloc(sizeof(float32)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayf;
    }
    if(HDFTYPE==DFNT_FLOAT64){
      arrayd=(float64*)malloc(sizeof(float64)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayd;
    }
    if(array==NULL){
      fprintf(stderr,"cannot allocate array hdf data\n");
      exit(1);
    }

    if(NDIM==3){
      dim_sizes[3] = numcolumns;
      dim_sizes[2] = N1;
      dim_sizes[1] = N2;
      dim_sizes[0] = N3;
      rank = 1+NDIM;
      
      edges[3] = numcolumns;
      edges[2] = N1;
      edges[1] = N2;
      edges[0] = N3;

      start[0]=start[1]=start[2]=start[3]=0;

    }
    else if(NDIM==2){
      // assum 2d
      dim_sizes[2] = numcolumns;
      dim_sizes[1] = N1;
      dim_sizes[0] = N2;
      rank = 1+NDIM;
      
      edges[2] = numcolumns;
      edges[1] = N1;
      edges[0] = N2;

      start[0]=start[1]=start[2]=0;
    }
    else if(NDIM==1){
      // assum 1d
      dim_sizes[1] = numcolumns;
      dim_sizes[0] = N1;
      rank = 1+NDIM;
      
      edges[1] = numcolumns;
      edges[0] = N1;

      start[0]=start[1]=0;
    }

  }
#endif




  /////////////////////////////
  //
  // vis5d stuff
  //
  /////////////////////////////
#if(V5D)
  if((source==5)||(dest==5)){
    fprintf(stderr,"Allocated memory for source=%d dest=%d\n",source,dest); fflush(stderr);
    arrayvisf=(float*)malloc(sizeof(float)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
    arrayvis=arrayvisf;
    if(arrayvis==NULL){
      fprintf(stderr,"cannot allocate array vis5d data\n");
      exit(1);
    }
    arrayvisoutput=(float*)malloc(sizeof(float)*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
    if(arrayvisoutput==NULL){
      fprintf(stderr,"cannot allocate array vis5d data: arrayvisoutput\n");
      exit(1);
    }
  }
#endif







  ///////////
  //
  //
  // BIG LOOP
  //
  ////////////

  for (it=0;it<TS;it++) { // loop over timeslices


    if(TS>1){ // otherwise already set
      if(inputlist) fscanf(inputfilenamelist,"%s",INPUT_NAME); // otherwise 1 file
      if(outputlist) fscanf(outputfilenamelist,"%s",OUTPUT_NAME); // otherwise 1 file
      fprintf(stderr,"At TS=%d of %d using input file %s and output file %s\n",it,TS,INPUT_NAME,OUTPUT_NAME);
    }

    fprintf(stderr,"Open source=%d\n",source);



    //////////////////////
    //    
    // open source file
    //
    //////////////////////
    if(source==0){
      
      // length of the entire unmodified input file name
      lname=strlen(INPUT_NAME);
      if( (INPUT_NAME[lname-1]=='z')&&(INPUT_NAME[lname-2]=='g')&&(INPUT_NAME[lname-3]=='.') ){
        inputtype=1;
        printf("input flagged as gzip\n");
      }
      else{
        inputtype=0;
      }
      
      if(inputtype==0){
        if( !(input=fopen(INPUT_NAME,"rb"))){
          fprintf(stderr,"trouble opening input file: %s\n",INPUT_NAME);
          exit(1);
        }
      }
      if(inputtype==1){
        sprintf(incommand,"gzip -d < %s",INPUT_NAME);
        if( !(input=popen(incommand,"r"))){
          fprintf(stderr,"trouble opening input file: %s %s\n",INPUT_NAME,incommand);
          exit(1);
        }
      }
      // assume header info provided and # lines of header provided (unlike in r8toras.c and block2tile.c)
      
    }
    else if(source==1){
      if( (input=fopen(INPUT_NAME,"rb"))==NULL){
        fprintf(stderr,"cannot open %s\n",INPUT_NAME);
        exit(1);
      }
    }
    else if(source==2){
      if( (input=fopen(INPUT_NAME,"rt"))==NULL){
        fprintf(stderr,"cannot open %s\n",INPUT_NAME);
        exit(1);
      }
    }
#if(HDF)
    else if(source==3){
      // open HDF file
      sd_id=SDstart(INPUT_NAME,DFACC_READ);
      if (sd_id != FAIL)      printf ("Reading HDF file with READ access\n");
      
      sd_index = 0;
      sds_id = SDselect (sd_id, sd_index);
      
      istat = SDreaddata (sds_id, start, NULL, edges, (VOIDP) array);
    }
    else if(source==4){
      // not yet
      fprintf(stderr,"NOT YET\n"); exit(1);
    }
#endif
#if(V5D)
    else if((source==5)&&(it==it)){ // No longer assume: // assume all input vis5d are multi-timed
      // not yet
      //      fprintf(stderr,"NOT YET\n"); exit(1);

      v5dstruct v;
      int time,var;
      float *data;
      float min, max, sum, sumsum;
      int missing, good;

      // code pulled from v5dstats.c

      fprintf(stderr,"BEGIN: v5dOpenFile\n"); fflush(stderr);
      if (!v5dOpenFile( INPUT_NAME, &v )) {
        printf("Error: couldn't open %s for reading\n", INPUT_NAME );
        exit(0);
      }
      fprintf(stderr,"END: v5dOpenFile\n"); fflush(stderr);

      if(v.NumTimes>1){
        fprintf(stderr,"Not quite setup for multi-timed inputs since need to increase size of array to have time dimension\n");
        fflush(stderr);
        exit(1);
      }

      // assume same size for all times and variables
      int sizeproblem=0;
      for (time=0; time<v.NumTimes; time++) {
        for (var=0; var<v.NumVars; var++) {
          // v5d data size
          int nrncnl;
          nrncnl = v.Nr * v.Nc * v.Nl[var];
          
          if(nrncnl!=N1*N2*N3){
            fprintf(stderr,"Memory created was per-time per-variable of size %d while v5d file had %d\n",N1*N2*N3,nrncnl);
            fflush(stderr);
            sizeproblem++;
          }
          else{
            // debug:
            //      fprintf(stderr,"time=%d var=%d size=%d\n",time,var,nrncnl); fflush(stderr);

            // read data into array
            //      data = (float *) malloc( nrncnl * sizeof(float) );
            data=arrayvisoutput;

            fprintf(stderr,"BEGIN: v5dReadGrid\n"); fflush(stderr);
            if (!v5dReadGrid( &v, time, var, data )) {
              printf("Error while reading grid (time=%d,var=%s)\n", time+1, v.VarName[var] );
              exit(0);
            }
            fprintf(stderr,"END: v5dReadGrid\n"); fflush(stderr);

            //      min = MISSING;
            //      max = -MISSING;
            //      missing = 0;
            //      good = 0;
            //      sum = 0.0;
            //      sumsum = 0.0;
            //
            //      for (i=0;i<nrncnl;i++) {
            //        /*
            //          if (data[i]!=data[i]) {
            //          printf("bad: %g\n", data[i]);
            //          }
            //        */
            //        if ( IS_MISSING(data[i]) ) {
            //          missing++;
            //        }
            //        else {
            //          good++;
            //          if (data[i]<min) {
            //                  min = data[i];
            //          }
            //          if (data[i]>max) {
            //                  max = data[i];
            //          }
            //          sum += data[i];
            //          sumsum += data[i]*data[i];
            //        }
            //      }


            if(dest==0){ // then use min/max conversion for decent legend (at least for fixed scaled data)
              // also assumes linear legend!
              a=minmax[0][var];
              b=minmax[1][var];
            }
            else{ // to cancel change
              a=0.0;
              b=255.0;
            }

            // loop over spatial dimensions //      for (i=0;i<nrncnl;i++)
            for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++){
                  ijkjon=(i+(j+k*N2)*N1)*((long long int)v.NumVars) + (long long int)var;  // so var (columns) is fastest, then i, then j, then k

                  // Nl,Level, Nc,Column, then Nr,Row in speed increasingly.
                  //                  ijkvis5d=k + ((i) + (j) * N1) * N3; // so y,k,N3,Nl fastest, then x,i,N1,Nc then z,j,N2,Nr
                  //                  ijkvis5d=N1-1-i + ((k) + (j) * N3) * N1; // so x,i,N1,Nl fastest, then z,k,N3,Nc then y,j,N2,Nr

                  //                  ijkvis5d=N1-1-i + ((j) + (k) * N2) * N1;


                  //                  ijkvis5d=(N3-1-k) + ((i) + (j) * N1) * N3; // so k fastest, then i then j
                  //ijkvis5d=(i+(j+k*N2)*N1);

                  ijkvis5d=(N3-1-k) + ((j) + (i) * N2) * N3;
                  
                  
                  arrayvisf[ijkjon] = (arrayvisoutput[ijkvis5d]-a)*255.0/(b-a);
                  //              fprintf(stderr,"k=%d j=%d i=%d : var=%d : ijkjon=%d\n",k,j,i,var,ijkjon);

                }

            //      free( data );

            //      if (good==0) {
            //        /* all missing */
            //        printf("%4d  %-8s %-5s  all missing values\n",
            //               time+1, v.VarName[var], v.Units[var] );
            //      }
            //      else {
            //        float mean = sum / good;
            //        float tmp = (sumsum - sum*sum/good) / (good-1);
            //        float sd;
            //        if (tmp<0.0) {
            //          sd = 0.0;
            //        }
            //        else {
            //          sd = sqrt( tmp );
            //        }
            //        printf("%4d  %-8s %-5s %13g%13g%13g%13g  %4d\n",
            //               time+1, v.VarName[var], v.Units[var],
            //               min,  max,  mean, sd,  missing );
            //      }



          }// end else if ok to read i,j,k block of data
        }// end loop over reading variables
      }// end loop over reading times




      if(sizeproblem!=0){
        fprintf(stderr,"Size problem: %d\n",sizeproblem);
        fflush(stderr);
        exit(1);
      }


      fprintf(stderr,"BEGIN: v5dCloseFile\n"); fflush(stderr);
      v5dCloseFile( &v );
      fprintf(stderr,"END: v5dCloseFile\n"); fflush(stderr);
      fprintf(stderr,"Closing vis5d+ file. Done reading vis5d+ file into array\n");  fflush(stderr);

    }
#endif


    //    fprintf(stderr,"Deal with header: pipeheader=%d\n",pipeheader);



    //////////////////////
    //    
    // trial open dest file when processing input list out to 1 file
    //
    //////////////////////
    if(it==0 && (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0)){
      if( !(output=fopen(OUTPUT_NAME,"wt"))){
        fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
        exit(1);
      }
      fclose(output);
    }






    /////////////////////////
    //  
    // deal with header (when doing inputlist, note that this only pipes first input file into any output file)
    //
    if(pipeheader && (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0)){
      if( !(output=fopen(OUTPUT_NAME,"wt"))){
        fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
        exit(1);
      }
    }

    for(i=0;i<realnumheader;i++){
      //printf("headerline: %i\n",i);
      while((ch=fgetc(input))!='\n'){
        if(pipeheader) fputc(ch,output);
        //printf("%c",ch);
      }
    }

    // close header part if opened
    if(pipeheader && (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) ){
      fprintf(output,"\n");
      fclose(output);
    }





    ///////////////////////////
    //
    // open destination file
    //
    ///////////////////////////



    if(dest==0 && ( (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) || outputlist==1) ){
      fprintf(stderr,"Open dest=%d\n",dest);

      lname=strlen(OUTPUT_NAME);
      if( (OUTPUT_NAME[lname-1]=='z')&&(OUTPUT_NAME[lname-2]=='g')&&(OUTPUT_NAME[lname-3]=='.') ){
        outputtype=1;
        printf("output flagged as gzip\n");
      }
      else{
        outputtype=0;
      }
      
      // open output file
      if(outputtype==0){
        if( !(output=fopen(OUTPUT_NAME,"ab"))){
          fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
          exit(1);
        }
      }
      if(outputtype==1){
        sprintf(outcommand,"gzip > %s",OUTPUT_NAME);
        if( !(output=popen(outcommand,"w"))){
          fprintf(stderr,"trouble opening output file: %s %s\n",OUTPUT_NAME,outcommand);
          exit(1);
        }
      }
    }
    else if(dest==1 && ( (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) || outputlist==1) ){
      fprintf(stderr,"Open dest=%d\n",dest);

      if( (output=fopen(OUTPUT_NAME,"at"))==NULL){
        fprintf(stderr,"cannot open %s\n",OUTPUT_NAME);
        exit(1);
      }
    }
    else if(dest==2 && ( (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) || outputlist==1) ){
      if( (output=fopen(OUTPUT_NAME,"at"))==NULL){
        fprintf(stderr,"cannot open %s\n",OUTPUT_NAME);
        exit(1);
      }
    }
#if(HDF)
    else if(dest==3 && ( (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) || outputlist==1) ){
      fprintf(stderr,"Open dest=%d\n",dest);

      // open HDF file
      sd_id = SDstart (OUTPUT_NAME, DFACC_CREATE);
      
      // can change to other output types
      
      sds_id = SDcreate (sd_id, SDS_NAME, HDFTYPE, rank, dim_sizes);
      start[0]=start[1]=start[2]=start[3]=0;
    }
    else if(dest==4 && ( (it==0 && outputlist==0 && inputlist==1 || outputlist==0 && inputlist==0) || outputlist==1) ){
      fprintf(stderr,"Open dest=%d\n",dest);

      // not yet
      fprintf(stderr,"NOT YET\n"); exit(1);
    }
#endif
#if(V5D)
    else if((dest==5)&&(it==0)){ // assume all output vis5d are multi-timed, so only open 1 file for all timeslices
      fprintf(stderr,"Open dest=%d\n",dest);

      fprintf(stderr,"Opening vis5d file: %s %d %d %d %d %d %s\n",OUTPUT_NAME,NumTimes,NumVars,Nr,Nc,Nl[0],VarName[0]); fflush(stderr);


      if(CompressMode<4){
        fprintf(stderr,"WARNING: using CompressMode=%d\n",CompressMode);
        fprintf(stderr,"Compression is per lev, but can make it difficult to measure small values if just happens to be large value somewhere on that lev.  For example, tracing vectors\n");
      }

      /* use the v5dCreate call to create the v5d file and write the header */
      fprintf(stderr,"BEGIN: v5dCreate\n"); fflush(stderr);
      if (!v5dCreate( OUTPUT_NAME, NumTimes, NumVars, Nr, Nc, Nl,
                      (const char (*)[MAXVARNAME]) VarName,
                      TimeStamp, DateStamp, CompressMode,
                      Projection, ProjArgs, Vertical, VertArgs )) {
        printf("Error: couldn't create %s\n", OUTPUT_NAME );
        exit(1);
      }
      fprintf(stderr,"END: v5dCreate\n"); fflush(stderr);
    }
#endif



    
      
    
    ////////////////////////////////////
    //
    // Process file
    //
    ////////////////////////////////////    


    fprintf(stderr,"reading file and putting into other format...\n"); fflush(stderr);
    fprintf(stderr,"newrows=%d numterms=%d group0=%d\n",numrows,numterms,group[0]); fflush(stderr);

    BUFFERINIT;// use to fill/read array if HDF involved
    for(i=0;i<numrows;i++){
#if(DEBUG)
      fprintf(stderr,"rownumber: %d of %d\n",i,numrows-1); fflush(stderr);
#endif
    
      for(j=0;j<numterms;j++){ // over rows and each group
#if(DEBUG)
        fprintf(stderr,"termnumber: %d of %d\n",j,numterms-1); fflush(stderr);
#endif
      
        for(k=0;k<group[j];k++){ // over a group of same kind
#if(DEBUG)
          fprintf(stderr,"groupelements: %d of %d : precision: %c s: %d d: %d\n",k,group[j]-1,precision[j],source,dest); fflush(stderr);
#endif


          //////////////
          // SOURCE
          //////////////
          switch(source){        
          case 0:
            fread(&dumb,bytesize,1,input); precision[j]='b'; //forced
            break;
          case 1:
            if(precision[j]=='b')                   fread(&dumb,bytesize,1,input);
            if(precision[j]=='i')                   fread(&dumi,intsize,1,input);
            if(precision[j]=='f')             fread(&dumf,floatsize,1,input);
            if(precision[j]=='d')             fread(&dumlf,doublesize,1,input);
            break;
          case 2:
            if(precision[j]=='b'){
              fscanf(input,"%hu",&shortfordumb);
              dumb=shortfordumb; // convert short to byte
            }
          
            if(precision[j]=='i')                   fscanf(input,"%d",&dumi);
            if(precision[j]=='f')             fscanf(input,"%f",&dumf);
            if(precision[j]=='d')             fscanf(input,"%lf",&dumlf);
            break;
#if(HDF)
          case 3:
          case 4:
            if(precision[j]=='b') dumb=arrayb[nextbuf++];
            if(precision[j]=='i') dumi=arrayi[nextbuf++];
            if(precision[j]=='f') dumf=arrayf[nextbuf++];
            if(precision[j]=='d') dumlf=arrayd[nextbuf++];
            break;
#endif
#if(V5D)
          case 5:
            // source is always float
            // nextbuf++ just iterates over fastest index, which is columns then i then j then k as normal for HARM
            //      fprintf(stderr,"nextbuf=%d\n",nextbuf); fflush(stderr);
            dumf=arrayvisf[nextbuf++]; precision[j]='f'; // forced
            break;
#endif
          default:
            break;
          }

          //////////////
          // DEST
          //////////////
          switch(dest){
          case 0:
            // always byte size output
            if(precision[j]=='b'){dumb=dumb;  fwrite(&dumb,bytesize,1,output);}
            if(precision[j]=='i'){dumb=dumi;  fwrite(&dumb,bytesize,1,output);}
            if(precision[j]=='f'){dumb=dumf;  fwrite(&dumb,bytesize,1,output);}
            if(precision[j]=='d'){dumb=dumlf; fwrite(&dumb,bytesize,1,output);}
            break;
          case 1:
            if(precision[j]=='b')                 fwrite(&dumb,bytesize,1,output);
            if(precision[j]=='i'){if(ble){ SWAP_LONG(dumi);}     fwrite(&dumi,intsize,1,output);}
            if(precision[j]=='f'){if(ble){ SWAP_FLOAT(dumf);}   fwrite(&dumf,floatsize,1,output);}
            if(precision[j]=='d'){if(ble){ SWAP_DOUBLE(dumlf);} fwrite(&dumlf,doublesize,1,output);}
            break;
          case 2:
            if(precision[j]=='b'){
              if(source==0 || source==1) fprintf(output,"%d ",dumb);
              else fprintf(output,"%c ",dumb);
            }
            if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}     fprintf(output,"%d ",dumi);}
            if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}   fprintf(output,"%17.10g ",dumf);}
            if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);} fprintf(output,"%26.20g ",dumlf);}
            break;
#if(HDF)
          case 3:
          case 4:
            if(precision[j]=='b') arrayb[nextbuf++]=dumb;
            if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}     arrayi[nextbuf++]=dumi;}
            if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}   arrayf[nextbuf++]=dumf;}
            if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);} arrayd[nextbuf++]=dumlf;}
            break;
#endif
#if(V5D)
          case 5:
            // same array, must be float
            if(precision[j]=='b') arrayvisf[nextbuf++]=dumb;
            if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}      arrayvisf[nextbuf++]=dumi;}
            if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}    arrayvisf[nextbuf++]=dumf;}
            if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);}  arrayvisf[nextbuf++]=dumlf;}
            break;
#endif
          default:
            break;
          }
        }//end over group
      }// end over terms
      // skip terms till carraige return (allows parsing out certain extra columns on end)
      if(source==2){ while(fgetc(input)!='\n'); }
      // output carraige return
      if(dest==2) fprintf(output,"\n");

    }// end over rows
  


#if(HDF)
    if(dest==3){
      fprintf(stderr,"Writing HDF file ...\n"); fflush(stderr);
      status = SDwritedata (sds_id, start, NULL, edges, (VOIDP)array); 
      status = SDendaccess (sds_id);
      status = SDend (sd_id);
    }
#endif
    //#if(V5D&&0) // DEBUG GODMARK
#if(V5D)
    if(dest==5){
      fprintf(stderr,"Writing vis5d file ...\n"); fflush(stderr);

      for(iv=0;iv<NumVars;iv++){

        /**
         ** Read your 3-D grid data for timestep it and variable
         ** iv into the array g here.
         ** To help with 3-D array indexing we've defined a macro G.
         ** G(0,0,0) is the north-west-bottom corner, G(Nr-1,Nc-1,Nl-1) is
         ** the south-east-top corner.  If you want a value to be considered
         ** missing, assign it equal to the constant MISSING.  For example:
         ** G(ir,ic,il) = MISSING;
         **/
        //#define G(ROW, COLUMN, LEVEL)   g[ (ROW) + ((COLUMN) + (LEVEL) * Nc) * Nr ]
      
      
        if(source==0){ // then use min/max conversion for decent legend (at least for fixed scaled data)
          // also assumes linear legend!
          a=minmax[0][iv];
          b=minmax[1][iv];
        }
        else{ // to cancel change
          a=0.0;
          b=255.0;
        }
        // convert my format to vis5d format
        for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++){
              ijkjon=(i+(j+k*N2)*N1)*(long long int)NumVars + (long long int)iv;
              //              ijkvis5d=k + ((j) + (i) * N2) * N3;
              //              ijkvis5d=k + ((i) + (j) * N1) * N3; // so y,k,N3,Nl fastest, then x,i,N1,Nc then z,j,N2,Nr
              //              ijkvis5d=N1-1-i + ((k) + (j) * N3) * N1; // so x,i,N1,Nl fastest, then z,k,N3,Nc then y,j,N2,Nr
              //              ijkvis5d=N1-1-i + ((N3-1-k) + (j) * N3) * N1; // so x,i,N1,Nl fastest, then z,k,N3,Nc then y,j,N2,Nr // so x cross y is z instead of -z
              //              ijkvis5d=N1-1-i + ((j) + (k) * N2) * N1;
              

              //              ijkvis5d=(N3-1-k) + ((i) + (j) * N1) * N3; // so k fastest, then i then j
              //ijkvis5d=(i+(j+k*N2)*N1);


              //              ijkvis5d=N1-1-i + ((j) + (k) * N2) * N1;
              
              ijkvis5d=N1-1-i + ((k) + (j) * N3) * N1;

              arrayvisoutput[ijkvis5d]=(arrayvisf[ijkjon]/255.0)*(b-a)+a;

              // DEBUG:
              //        fprintf(stderr,"i=%d j=%d k=%d iv=%d ijkvis5d=%d a=%g b=%g value=%g\n",i,j,k,iv,ijkvis5d,a,b,arrayvisoutput[ijkvis5d]);
            }
      
        /* Write data to v5d file. */
        fprintf(stderr,"BEGIN: v5dWrite: it=%d iv=%d\n",it,iv); fflush(stderr);
        if (!v5dWrite( it+1, iv+1, arrayvisoutput )) {
          printf("Error while writing grid.  Disk full?\n");
          exit(1);
        }
        fprintf(stderr,"END: v5dWrite: it=%d iv=%d\n",it,iv); fflush(stderr);
      }
    }
#endif

  

    /////////////////////
    //
    // close output file
    //
    /////////////////////
    if(it==TS-1 && outputlist==0 && inputlist==1 || outputlist==1){
      if((dest==1)||(dest==2)) fclose(output);
      else if(dest==0){
        if(outputtype==0) fclose(output);
        else if(outputtype==1) pclose(output);
      }
    }


#if(HDF)
    if(source==3){
      /* Terminate access to the array. */
      istat = SDendaccess(sds_id);
    
      /* Terminate access to the SD interface and close the file. */
      istat = SDend(sd_id);
      if (istat != FAIL) printf("... file closed\n\n");
    }
    else if(source==4){
      // not yet
      fprintf(stderr,"NOT YET\n"); exit(1);
    }
#endif


    /////////////////////
    //
    // close input file
    //
    /////////////////////
    if(it==TS-1 && outputlist==1 && inputlist==0 || inputlist==1){
      if((source==1)||(source==2)) fclose(input);
      else if(source==0){
        if(inputtype==0) fclose(input);
        else if(inputtype==1) pclose(input);
      }
    }



  } // END BIG LOOP over all timesteps





  /////////////////
  //
  // vis5d stuff
  //
  /////////////////
#if(V5D)
  // close v5d file which has entire time series in it
  if(dest==5){
    fprintf(stderr,"BEGIN: v5Close\n"); fflush(stderr);
    v5dClose();
    fprintf(stderr,"END: v5Close\n"); fflush(stderr);
    free(arrayvisf);
    free(arrayvisoutput);
  }
  if(source==5){ // close entire time series v5d
    // not yet
    free(arrayvisf);
    free(arrayvisoutput);
  }
#endif


  if(TS>1){
    if(inputlist){ // otherwise 1 file
      fclose(inputfilenamelist);
    }
    if(outputlist){ // otherwise 1 file
      fclose(outputfilenamelist);
    }
  }

  

  fprintf(stderr,"done.\n");
  return(0);
  //  exit(0); // not reachable
}










// apparently doesn't work:
void swap(BYTE *x, BYTE size)
{
  unsigned char c;
  unsigned short s;
  unsigned long l;

  switch (size)
    {
    case 1: /* don't do anything */
      break;
    case 2: /* swap two bytes */
      c = *x;
      *x = *(x+1);
      *(x+1) = c;
      break;
    case 4: /* swap two shorts (2-byte words) */
      s = *(unsigned short *)x;
      *(unsigned short *)x = *((unsigned short *)x + 1);
      *((unsigned short *)x + 1) = s;
      swap ((BYTE *)x, 2);
      swap ((BYTE *)((unsigned short *)x+1), 2);
      break;
    case 8: /* swap two longs (4-bytes words) */
      l = *(unsigned long *)x;
      *(unsigned long *)x = *((unsigned long *)x + 1);
      *((unsigned long *)x + 1) = l;
      swap ((BYTE *)x, 4);
      swap ((BYTE *)((unsigned long *)x+1), 4);
      break;
    }
}




int machineEndianness(void)
{
  int i = 1;
  char *p = (char *) &i;
  if (p[0] == 1){
    // Lowest address contains the least significant byte
    //     fprintf(stderr,"little\n");
    return LITTLE_ENDIAN_USER;
  }
  else{
    //     fprintf(stderr,"big\n");
    return BIG_ENDIAN_USER;
  }

  //   return(0); // not reachable
}



//#include "SwapEndian.h"





static long _TestEndian=1;

int IsLittleEndian(void) {
  return *(char*)&_TestEndian;
}





/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

void *SwapEndian(void* Addr, const int Nb) {
  static char Swapped[16];
  switch (Nb) {
  case 2:       Swapped[0]=*((char*)Addr+1);
    Swapped[1]=*((char*)Addr  );
    break;
  case 3:       // As far as I know, 3 is used only with RGB images
    Swapped[0]=*((char*)Addr+2);
    Swapped[1]=*((char*)Addr+1);
    Swapped[2]=*((char*)Addr  );
    break;
  case 4:       Swapped[0]=*((char*)Addr+3);
    Swapped[1]=*((char*)Addr+2);
    Swapped[2]=*((char*)Addr+1);
    Swapped[3]=*((char*)Addr  );
    break;
  case 8:       Swapped[0]=*((char*)Addr+7);
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
