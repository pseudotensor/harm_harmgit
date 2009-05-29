/*   
	convert "r8" file to sun raster file or view it too
	see USAGE below for usage

	to compile:
	gcc -lm -O3 -Wall -o r8toras r8toras.c
	on alpha:
	ccc -lm -O3 -Wall -o r8toras r8toras.c
        scp r8toras root@metric:/usr/local/bin   

	jcm 07/26/2000
	jcm 08/21/2001
	jcm 09/15/2001
*/

#define MAXFNAME  400		/* max length of filenames */
#define FTYPE float
#include <stdlib.h>  
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rasterfile.h"

int main(int argc, char *argv[])
{
  int itemp;
  int usetmpdir;
  FILE * infile;
  FILE * outfile;
  int inputtype,outputtype;

  int nx, ny, nmap ;
  int nxoutput;
  unsigned char dum ;
  void decint(int *list, int num) ;
  struct rasterfile header ;
  FILE *palfile ;
  char incommand[256];
  char outcommand[256];
  char *outname;
  char newname[256];
  int lname;
  int i,j;
  int view;
  unsigned char firstfew[256];
  int numfirst;
  int numcolors;
  int goodheader;
  int endinput;
  fpos_t beginposition;
  long beginpos,endpos;

  // fuck me, why!  Is outname[256] a constant pointer?
  outname=(char*)malloc(sizeof(char)*256);
  if(outname==NULL){
    printf("death\n");
    exit(1);
  }

#define USAGE \
 fprintf(stderr,"usage: r8toras <outputtype> palfile infile (nx) (ny)\n") ; \
    fprintf(stderr,"or     r8toras <outputtype> palfile infile\n") ; \
    fprintf(stderr,"         if using proper raw header\n") ; \
    fprintf(stderr," 0-> no gzip, no view 1-> gzip, no view 2->no gzip, view 3->gzip,view\n") ;

  if(argc < 4) {
    USAGE
    exit(0) ;
  }
  outputtype=atoi(argv[1]);
  // 0-> no gzip, no view
  // 1-> gzip, no view
  // 2-> no gzip, view
  // 3-> gzip, view
  if((outputtype==2)||(outputtype==3)){
    view=1;
    if(outputtype==2)      outputtype=0;
    if(outputtype==3)      outputtype=1;
  }
  else view=0;

  palfile = fopen(argv[2],"r") ;
  if(palfile == NULL) {
    fprintf(stderr,"trouble opening palette file: %s\n",argv[2]) ;
    exit(0) ;
  }
  
  // length of the entire unmodified input file name
  lname=strlen(argv[3]);
  if( (argv[3][lname-1]=='z')&&(argv[3][lname-2]=='g')&&(argv[3][lname-3]=='.') ){
    inputtype=1;
    endinput=lname-3;
    printf("flagged as gzip\n");
  }
  else{
    inputtype=0;
    endinput=lname;
  }

  // setup output file
  if( (argv[3][endinput-1]=='8')&&(argv[3][endinput-2]=='r')&&(argv[3][endinput-3]=='.') ){
    strcpy(outname,argv[3]);
    outname[endinput-3]='\0'; // truncate ending
    if(outputtype==1) strcat(outname,".ras.gz");
    else strcat(outname,".ras");
    printf("flagged as r8\n");
  }
  else{
    //fprintf(stderr,"input file must have .r8 or .r8.gz extension\n");
    //return(1);
    fprintf(stderr,"Assuming file is a .r8\n");
    strcpy(outname,argv[3]);
    if(outputtype==1) strcat(outname,".ras.gz");
    else strcat(outname,".ras");
  }
  usetmpdir=0;
  // open output file
  if(outputtype==0){
    if( !(outfile=fopen(outname,"wb"))){
      fprintf(stderr,"trouble opening output file: %s, attempting to use tmp dir\n",outname);
      usetmpdir=1;
    }
  }
  if(outputtype==1){
    sprintf(outcommand,"gzip > %s",outname);
    if( !(outfile=popen(outcommand,"w"))){
      fprintf(stderr,"trouble opening output file: %s %s, attempting to use tmp dir\n",outname,outcommand);
      usetmpdir=2;
    }
  }
  if(usetmpdir){
    itemp=strlen(outname);
    for(i=itemp;i>=0;i--){ // start at '\0'
      if(outname[i]=='/') break;
    }
    
    outname+=(i+1); // shift pointer to start after first backslash
    strcat(newname,"/tmp/");
    strcat(newname,outname);
    strcpy(outname,newname);
    if(usetmpdir==1){ // then outname change and reopen
      if( !(outfile=fopen(newname,"wb"))){
	fprintf(stderr,"trouble opening output file: %s\n",newname);
	exit(1);
      }
    }
    if(usetmpdir==2){ // then outname change and reopen
      sprintf(outcommand,"gzip > %s",newname);
      if( !(outfile=popen(outcommand,"w"))){
	fprintf(stderr,"trouble opening output file: %s %s, attempting to use tmp dir\n",newname,outcommand);
	exit(1);
      }
    }
  }


  if(inputtype==0){
    if( !(infile=fopen(argv[3],"rb"))){
      fprintf(stderr,"trouble opening input file: %s\n",argv[3]);
      exit(1);
    }
  }
  if(inputtype==1){
    sprintf(incommand,"gzip -d < %s",argv[3]);
    if( !(infile=popen(incommand,"r"))){
      fprintf(stderr,"trouble opening input file: %s %s\n",argv[3],incommand);
      exit(1);
    }
  }
  // determine what kind of input file we have
  for(i=0;i<5;i++){
    firstfew[i]=fgetc(infile);
  }
  // assume likelyhood of "RAW\n#" appearing in image normally is very near 0
  if( (firstfew[0]=='R')&&(firstfew[1]=='A')&&(firstfew[2]=='W')&&(firstfew[3]=='\n')&&(firstfew[4]=='#') ){
    goodheader=1;
    numfirst=0;
    // read rest of RAW header of input image
    while(fgetc(infile)!='\n'); // skip to next line(skipping comment line)
    fscanf(infile,"%d %d",&nx,&ny); while(fgetc(infile)!='\n'); // skip to next line
    fscanf(infile,"%d",&numcolors); while(fgetc(infile)!='\n'); // skip to next line
  }
  else{
    printf("flagged as totally raw\n");
    goodheader=0; // so assume totally raw
    numfirst=5; // so now need to image these inputs
    if(argc==6){
      nx=atoi(argv[4]); // so size must be on command line
      ny=atoi(argv[5]);
    }
    else{
      if(inputtype==1){
	fprintf(stderr,"Must include nx and ny on command line for this type of input file\n");
	USAGE
	exit(0) ;
      }
      else{
	// if no args given, assume square size
	fprintf(stderr,"no nx ny given, assuming square\n"); fflush(stderr);
	beginpos=ftell(infile);
	fgetpos(infile,&beginposition);
	fseek(infile,0,SEEK_END);
	endpos=ftell(infile);
	fsetpos(infile,&beginposition);
	nx=ny=(int)sqrt((FTYPE)(5)+(FTYPE)(endpos-beginpos+1));      
      }
    }
  }

  // now have all files opened
  fprintf(stderr,"nx: %d ny: %d ... ",nx,ny);

  
  
  nxoutput=nx+nx%2;
  /* set up header */
  if(nxoutput!=nx){
    fprintf(stderr,"nx adjusted to be multiple of 16bits\n"); fflush(stderr);
  }

  header.ras_magic = RAS_MAGIC ;
  header.ras_width = nxoutput ;
  header.ras_height = ny ;
  header.ras_depth = 8 ;
  header.ras_length = nxoutput*ny ;
  header.ras_type = RT_STANDARD ;
  header.ras_maptype = RMT_EQUAL_RGB ;
  header.ras_maplength = 3*256 ;
  
  /* swap & write it */
  decint(&(header.ras_magic),8) ;
  fwrite(&(header.ras_magic), sizeof(int), 8, outfile) ;
  
  /* read & write the palette file */
  nmap = 0 ;
  while(fread(&dum, sizeof(unsigned char), 1, palfile) == 1) {
    fwrite(&dum, sizeof(unsigned char), 1, outfile) ;
    nmap++ ;
  }
  decint(&(header.ras_maplength),1) ;
  if(nmap != header.ras_maplength) {
    fprintf(stderr,"err: header: %d, nmap: %d\n",
	    header.ras_maplength, nmap) ;
  }
  
  /* now read & write the image */
  for(i=0;i<numfirst;i++){
    fwrite(&firstfew[i], sizeof(unsigned char), 1, outfile) ;
  }
  if(nxoutput!=nx){
    i=numfirst; // as just above would do anyways
    j=0;
  }
  while(fread(&dum, sizeof(unsigned char), 1, infile) == 1) {
    fwrite(&dum, sizeof(unsigned char), 1, outfile) ;
    if(nxoutput!=nx){
      i++;
      if(i==nx){
	i=0;
	j++;
	dum=0;
	fwrite(&dum, sizeof(unsigned char), 1, outfile) ; // extra padding to make multiple of 16bits
      }
    }
  }
  
  fprintf(stderr,"done\n") ;
  fclose(palfile);
  if(inputtype==0){
    fclose(infile);
  }
  else pclose(infile);
  if(outputtype==0){
    fclose(outfile);
  }
  else pclose(outfile);

  //now view if asked to
  if(view==1){
    sprintf(outcommand,"display %s\n",outname);
    system(outcommand);
    // assume only want to view
    sprintf(outcommand,"rm %s\n",outname);
    system(outcommand);
  }
  return(0) ;
}

void decint(int *lp, int n)
{
  unsigned int t;
  static unsigned long lmask = 0x00ff0000, rmask = 0x0000ff00;
  
  for(; n--; lp++) {
    t = *lp;
    *lp = (t >> 24) | (t << 24) | ((t << 8) & lmask) |
      ((t >> 8) & rmask);
  }
}

