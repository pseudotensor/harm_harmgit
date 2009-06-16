// compile with:
//   icc -Wall -O3 -axiMKW -unroll -ipo -rcd -tpp7 -march=pentium4 -mcpu=pentium4 -lm -o smcalc smcalc.c
//
// or just:
// gcc -O3 -Wall -o smcalc smcalc.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define sign(a) (copysign(1.0,a))

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#define RHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7
#define SMALL 1E-30
#define THIRD          0.3333333333333333333333333333333333L  /* 1/3 */

#include "global.nondepnmemonics.h"
#define WHICHPARA PARA4
#include "para_and_paraenohybrid.h"

#define NDIM 4
#define NPR 8
#define FTYPE double

#define PLOOP for(k=0;k<NPR;k++)

#define MAP(k,max) (j*nx*max+i*max+(k))
#define MAPGEN(i,j,k,max) ((j)*nx*max+(i)*(max)+(k))
#define MAPIP1(k,max) (j*nx*max+(i+1)*max+k)
#define MAPIM1(k,max) (j*nx*max+(i-1)*max+k)
#define MAPJP1(k,max) ((j+1)*nx*max+i*max+k)
#define MAPJM1(k,max) ((j-1)*nx*max+i*max+k)

#define MAP1 MAP(0,1)
#define MAPDER(k) MAP(k,2)
#define MAPJ(k) MAP(k,4)
#define MAPF(k) MAP(k,6)
#define MAPG(k) MAP(k,numgenvars)
#define MAPFICALC(k) MAP(k,2)

#define MAPP(k) MAP(k,NPR)
#define MAPPGEN(i,j,k) MAPGEN(i,j,k,NPR)
#define MAPPIP1(k) MAPIP1(k,NPR)
#define MAPPIM1(k) MAPIM1(k,NPR)
#define MAPPJP1(k) MAPJP1(k,NPR)
#define MAPPJM1(k) MAPJM1(k,NPR)

#define READHEADER 1 // should always be 1 generally
#define WRITEHEADER 1

#define NUMDUMPCOL1 36 // number of columns including i,j
#define NUMDUMPCOL2 56 // number of columns including i,j
// could be 56 too for jcon/jcov/fcon/fcov stuff if included in dump file

// DUMPTYPE 
#define REALDUMP 0
#define SINGLECOLUMN 1
#define ANYCOLUMN 2   // 2 or greater number of columns
// 0: real dump file
// 1: forfline dump file as in gammie.m (forfldump)
// 2: number of columns given for arbitrary size


//CALCTYPE
#define DOFARADAY 0
#define DOFLINE 1
#define DODER 2 // only use with DUMPTYPE==SINGLECOLUMN
#define DOAVG 3
#define DOFICALC 4


static FTYPE Ficalc(int dimen, FTYPE *V, FTYPE *P);
static FTYPE ftilde( int dimen, int shift, FTYPE *Vabs, FTYPE *Pabs);


int main(
        int argc,
        char *argv[],
        char *envp[]
        )
{
  int k;
  FILE * dumpin;
  FILE *dumpout;
  FILE *colin;
  char filenamein[200];
  char filenameout[200];
  char dumpnamein[200];
  char dumpnameout[200];
  char systemcall[200];
  int i,j;
  int ii,jj;
  int si,sj;
  int iii,jjj;
  FTYPE X[NDIM],r,th,divb,ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],vmin1,vmax1,vmin2,vmax2;
  int totalsize[NDIM];
  long realnstep;
  FTYPE readnstep;
  FTYPE t,startx[NDIM],dx[NDIM],gam,a,R0,Rin,Rout,hslope,dt;
  int defcoord;
  FTYPE *p,*gdet,*aphi,*da,*db,*result1,*faradaydd,*faradayuu,*jcon,*jcov,*faradaydotj,*der;
  FTYPE *Ptot,*ficalc;
  int nx,ny;
  int indexi,indexj;
  FTYPE result2;
  int DUMPTYPE,CALCTYPE,DERTYPE;
  FTYPE *genvar;
  FTYPE ftemp;
  int numgenvars;
  int dumpnum,startdump,enddump;
  int numcol;
  int inow,jnow;
  FTYPE (*Ypl)[20],*V,*P;
  FTYPE a_Ypl[NPR][20],a_V[20],a_P[20];
  int dimen,ismall;


  // shift sufficiently
  Ypl=(FTYPE (*) [20]) (&(a_Ypl[0][10]));
  V=&a_V[10];
  P=&a_P[10];


  if(argc<8){
    fprintf(stderr,"./smcalc DUMPTYPE CALCTYPE DERTYPE nx ny inname outname <startdump> <enddump>\n");
    fprintf(stderr,"DUMPTYPE: 0=realdump 1=single column >=2 any number of columns specified\n");
    fprintf(stderr,"CALCTYPE: 0=faraday 1=fline 2=derivatives 3=average 4=ficalc1,2\n");
    fprintf(stderr,"DERTYPE=0->centered when possible 1->numerical backward if possible 2=forward if possible 10+ of 0 11=+ of 1 12=+ of 2\n");
    exit(1);
  }


  DUMPTYPE=atoi(argv[1]);
  CALCTYPE=atoi(argv[2]);
  DERTYPE=atoi(argv[3]);
  nx=atoi(argv[4]);
  ny=atoi(argv[5]);
  strcpy(filenamein,argv[6]);
  strcpy(filenameout,argv[7]);
  if(CALCTYPE==DOAVG){
    if(argc<10){
      fprintf(stderr,"CALCTYPE==DOAVG\n");
      fprintf(stderr,"./smcalc DUMPTYPE CALCTYPE nx ny inname outname <startdump> <enddump>\n");
      fprintf(stderr,"DUMPTYPE: 0=realdump 1=single column >=2 any number of columns specified\n");
      fprintf(stderr,"CALCTYPE: 0=faraday 1=fline 2=derivatives 3=average\n");
      fprintf(stderr,"DERTYPE=0->centered when possible 1->numerical backward if possible 2=forward if possible 10+ of 0 11=+ of 1 12=+ of 2\n");
      exit(1);
    }
    startdump=atoi(argv[8]);
    enddump=atoi(argv[9]);
  }
  else{
    startdump=0;
    enddump=0;
  }

  fprintf(stderr,"nx=%d ny=%d dumptype=%d calctype=%d dertype=%d\n",nx,ny,DUMPTYPE,CALCTYPE,DERTYPE); fflush(stderr);

  if(CALCTYPE==DOAVG){
    if(DUMPTYPE>=ANYCOLUMN) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*DUMPTYPE);
    else if(DUMPTYPE==SINGLECOLUMN) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*(DUMPTYPE+2));
    else if(DUMPTYPE==REALDUMP) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*NUMDUMPCOL2); // bigger of 2

    if(genvar==NULL){
      fprintf(stderr,"problems allocating genvar memory\n"); fflush(stderr); exit(1);
    }

    if((DUMPTYPE>=ANYCOLUMN)||(DUMPTYPE==SINGLECOLUMN)){
      numgenvars=DUMPTYPE;
    }
    else if(DUMPTYPE==REALDUMP){
      numgenvars=NUMDUMPCOL2;
    }
    fprintf(stderr,"numgenvars=%d\n",numgenvars);

    for(k=0;k<nx*ny*numgenvars;k++)  genvar[k]=0;

    sprintf(dumpnameout,"%s%s",filenameout,"avg");
    if((dumpout=fopen(dumpnameout,"wt"))==NULL){
      fprintf(stderr,"cannot open %s\n",dumpnameout);
      exit(1);
    }
  }
  else{
    if((dumpout=fopen(filenameout,"wt"))==NULL){
      fprintf(stderr,"cannot open %s\n",filenameout);
      exit(1);
    }
  }

  p=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*NPR);
  gdet=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny);
  if(CALCTYPE==DOFLINE){
    aphi=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny);
    da=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny);
    db=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny);
    result1=(FTYPE *)malloc(sizeof(FTYPE)*nx);
  }
  else if(CALCTYPE==DOFARADAY){
    faradaydd=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*6);
    faradayuu=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*6);
  }
  else if(CALCTYPE==DODER){
    der=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*2);
  }
  else if(CALCTYPE==DOFICALC){
    ficalc=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*2);
    Ptot=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*1);
  }

  for(dumpnum=startdump;dumpnum<=enddump;dumpnum++){

    if((CALCTYPE!=DOAVG)||(startdump==enddump)){
      if((dumpin=fopen(filenamein,"rt"))==NULL){
	fprintf(stderr,"problem with file opening %s\n",filenamein);
	exit(1);
      }
      sprintf(systemcall,"cat %s | head -2 | tail -1 | wc -w > numcol.txt",filenamein);
      system(systemcall);
    }
    else{
      sprintf(dumpnamein,"%s%04d",filenamein,dumpnum);
      if((dumpin=fopen(dumpnamein,"rt"))==NULL){
	fprintf(stderr,"problem with file opening %s\n",dumpnamein);
	exit(1);
      }
      sprintf(systemcall,"cat %s | head -2 | tail -1 | wc -w > numcol.txt",dumpnamein);
      system(systemcall);
    }
    colin=fopen("numcol.txt","rt");
    fscanf(colin,"%d",&numcol);
    fclose(colin);
    if(CALCTYPE!=DOFLINE){
      if((DUMPTYPE==SINGLECOLUMN)&&(numcol!=3)){
	fprintf(stderr,"numcol=%d for singlecolumn, should be 3\n",numcol);
	//      exit(1);
      }
    }
    else{
      if((DUMPTYPE==SINGLECOLUMN)&&(numcol!=6)){
	fprintf(stderr,"numcol=%d for singlecolumn/dofline, should be 6\n",numcol);
	//      exit(1);
      }
    }
    if(CALCTYPE==DOAVG){
      if(numcol==NUMDUMPCOL1) numgenvars=NUMDUMPCOL1;
      else if(numcol==NUMDUMPCOL2) numgenvars=NUMDUMPCOL2;
      else{
	fprintf(stderr,"numcol=%d for realdump\n",numcol);
	exit(1);
      }
    }
    
						 
      
#if(READHEADER)
    fprintf(stderr,"start reading header\n"); fflush(stderr);
    fscanf(dumpin, "%lf %d %d %lf %lf %lf %lf %ld %lf %lf %lf %lf %lf %lf %lf %d", &t,&totalsize[1],&totalsize[2],&startx[1],&startx[2],&dx[1],&dx[2],&realnstep,&gam,&a,&R0,&Rin,&Rout,&hslope,&dt,&defcoord);
    
    //realnstep=(long)readnstep;
    while(fgetc(dumpin)!='\n');
#endif
#if(WRITEHEADER)
    if(dumpnum==startdump){
      fprintf(stderr,"start writing header"); fflush(stderr);
      //      fprintf(stderr,"dt=%g defcoord=%d\n",dt,defcoord);
      fprintf(dumpout, "%10.5g %d %d %10.5g %10.5g %10.5g %10.5g %ld %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %d\n", t,totalsize[1], totalsize[2], startx[1], startx[2], dx[1],dx[2],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord); fflush(dumpout);
    }
#endif

    if(CALCTYPE==DOAVG){
      fprintf(stderr,"start reading file: dumpnum=%d dumpname=%s\n",dumpnum,dumpnamein); fflush(stderr);
    }
    else{
      fprintf(stderr,"start reading file: dumpname=%s\n",filenamein); fflush(stderr);
    }

    // LOOP over variables
    for(iii=0;iii<nx*ny;iii++){
      
      // all files assumed to have i and j
      fscanf(dumpin,"%d",&ii); //i
      fscanf(dumpin,"%d",&jj); //j

      // below accounts for if put in boundary zones.  Assumes now iii==0 will be first zone!  This wasn't assumed previously
      if(iii==0){
	si=ii;
	sj=jj;
      }
      i=ii-si;
      j=jj-sj;

      if((i<0)||(j<0)||(i>=nx)||(j>=ny)){
	fprintf(stderr,"data requested larger than said size: i=%d j=%d\n",i,j); fflush(stderr);
	exit(1);
      }
      
      // read in different file types
      if(DUMPTYPE==REALDUMP){
	if(CALCTYPE!=DOAVG){
	  // 2+2+NPR+1+4*4+4+1=34
	  for(k=1;k<=2;k++) fscanf(dumpin,"%lf",&X[k]);
	  fscanf(dumpin,"%lf",&r);
	  fscanf(dumpin,"%lf",&th);
	  PLOOP fscanf(dumpin,"%lf", &p[MAPP(k)]);
	  fscanf(dumpin,"%lf",&divb);
	  
	  for (k = 0; k < NDIM; k++)
	    fscanf(dumpin,"%lf",&(ucon[k]));
	  for (k = 0; k < NDIM; k++)
	    fscanf(dumpin,"%lf",&(ucov[k]));
	  for (k = 0; k < NDIM; k++)
	    fscanf(dumpin,"%lf",&(bcon[k]));
	  for (k = 0; k < NDIM; k++)
	    fscanf(dumpin,"%lf",&(bcov[k]));
	  
	  fscanf(dumpin,"%lf",&vmin1);
	  fscanf(dumpin,"%lf",&vmax1);
	  fscanf(dumpin,"%lf",&vmin2);
	  fscanf(dumpin,"%lf",&vmax2);
	  
	  fscanf(dumpin,"%lf",&gdet[MAP1]);
	}
	else{
	  genvar[MAPG(0)]=i+genvar[MAPG(0)];
	  genvar[MAPG(1)]=j+genvar[MAPG(1)];
	  for(k=2;k<numgenvars;k++){
	    fscanf(dumpin,"%lf", &ftemp);
	    genvar[MAPG(k)]=ftemp+genvar[MAPG(k)];
	  }
	  
	}
      }
      else if(DUMPTYPE==SINGLECOLUMN){
	if(CALCTYPE==DOFLINE){
	  for(k=5;k<=7;k++) fscanf(dumpin,"%lf", &p[MAPP(k)]);
	  //	  for(k=5;k<=7;k++) fprintf(stderr,"%lf ", p[MAPP(k)]);
	  fscanf(dumpin,"%lf",&gdet[MAP1]);
	  //	  fprintf(stderr,"%lf \n",gdet[MAP1]);
	}
	else if(CALCTYPE==DODER){ // dummy variable position
	  for(k=0;k<=0;k++) fscanf(dumpin,"%lf", &p[MAPP(k)]);
	}
	else if(CALCTYPE==DOAVG){ // dummy variable position
	  genvar[MAPG(0)]=i+genvar[MAPG(0)];
	  genvar[MAPG(1)]=j+genvar[MAPG(1)];
	  for(k=2;k<numgenvars;k++){
	    fscanf(dumpin,"%lf", &ftemp);
	    genvar[MAPG(k)]=ftemp+genvar[MAPG(k)];
	  }
	}     
      }
      else if(DUMPTYPE>=ANYCOLUMN){
	if(CALCTYPE==DOAVG){ // dummy variable position
	  //fprintf(stderr,"i=%d j=%d mapg(0)=%d\n",i,j,MAPG(0));
	  genvar[MAPG(0)]=i+genvar[MAPG(0)];
	  genvar[MAPG(1)]=j+genvar[MAPG(1)];
	  for(k=2;k<numgenvars;k++){
	    fscanf(dumpin,"%lf", &ftemp);
	    //fprintf(stderr,"mapg(%d)=%d\n",k,MAPG(k));
	    genvar[MAPG(k)]=ftemp+genvar[MAPG(k)];
	  }
	}
	else if(CALCTYPE==DOFICALC){ // dummy variable position
	  for(k=0;k<NPR;k++){
	    fscanf(dumpin,"%lf", &p[MAPP(k)]);
	    //	    fprintf(stderr,"iii=%d ii=%d jj=%d i=%d j=%d :: p[%d]=%21.15g\n",iii,ii,jj,i,j,k,p[MAPP(k)]);
	  }
	  fscanf(dumpin,"%lf",&Ptot[MAP1]);
	  fscanf(dumpin,"%lf",&gdet[MAP1]);
	}
      }
      

      // perform read-time calculations
      if(CALCTYPE==DOFLINE){
	// field line stuff
	da[MAP1]=-p[MAPP(6)]*gdet[MAP1]*dx[1];
	db[MAP1]=p[MAPP(5)]*gdet[MAP1]*dx[2];
	//fprintf(stdout,"%g %g :: %g %g :: %g %g\n",p[MAPP(6)],p[MAPP(5)],dx[1],dx[2],da[MAP1],db[MAP1]);
      }
      else if(CALCTYPE==DOFARADAY){
	// faraday (down down)
	ftemp=gdet[MAP1];
	faradaydd[MAPF(0)]=ftemp*(bcon[3]*ucon[2]-bcon[2]*ucon[3]);
	faradaydd[MAPF(1)]=ftemp*(bcon[1]*ucon[3]-bcon[3]*ucon[1]);
	faradaydd[MAPF(2)]=ftemp*(bcon[2]*ucon[1]-bcon[1]*ucon[2]);
	faradaydd[MAPF(3)]=ftemp*(bcon[3]*ucon[0]-bcon[0]*ucon[3]);
	faradaydd[MAPF(4)]=ftemp*(bcon[0]*ucon[2]-bcon[2]*ucon[0]);
	faradaydd[MAPF(5)]=ftemp*(bcon[1]*ucon[0]-bcon[0]*ucon[1]);
	
	// faraday (up up)
	ftemp=-1.0/ftemp;
	faradayuu[MAPF(0)]=ftemp*(bcov[3]*ucov[2]-bcov[2]*ucov[3]);
	faradayuu[MAPF(1)]=ftemp*(bcov[1]*ucov[3]-bcov[3]*ucov[1]);
	faradayuu[MAPF(2)]=ftemp*(bcov[2]*ucov[1]-bcov[1]*ucov[2]);
	faradayuu[MAPF(3)]=ftemp*(bcov[3]*ucov[0]-bcov[0]*ucov[3]);
	faradayuu[MAPF(4)]=ftemp*(bcov[0]*ucov[2]-bcov[2]*ucov[0]);
	faradayuu[MAPF(5)]=ftemp*(bcov[1]*ucov[0]-bcov[0]*ucov[1]);
      }

    }// end over variable loop
    
    fclose(dumpin);
    fprintf(stderr,"done reading file\n"); fflush(stderr);
  }// end over dump loop


  fprintf(stderr,"start post-read calculations\n"); fflush(stderr);

  ///////////////////////////////////
  // perform post-read calculations

  if(CALCTYPE==DOAVG){
    fprintf(stderr,"done reading all files\n"); fflush(stderr);

    for(k=0;k<nx*ny*numgenvars;k++){
      genvar[k]=genvar[k]/(enddump-startdump+1);
    }
  }
  else if(CALCTYPE==DOFLINE){
#if(0)
    // first set 0 point
    for(iii=0;iii<nx*ny;iii++)    aphi[iii]=-da[nx*ny-1];
    
    // generate x-direction integral since 1-D
    for(iii=0;iii<nx;iii++){
      result1[iii]=0;    for(jjj=iii;jjj<nx;jjj++)      result1[iii]+=-da[0*nx+jjj];
    }
    // now find i,j'th aphi
    for(iii=0;iii<nx*ny;iii++){
      indexi=(int)(iii%nx);    indexj=(int)(iii/nx);
      
      result2=0;    for(jjj=0;jjj<indexj;jjj++)      result2+=db[jjj*nx+indexi];
      aphi[indexj*nx+indexi]+=result1[indexi]+result2;
    }
#elif(0)
    // first set 0 point
    for(iii=0;iii<nx*ny;iii++)    aphi[iii]=-da[0];
    //for(iii=0;iii<nx*ny;iii++)    aphi[iii]=0;
    
    // generate x-direction integral since 1-D
    for(iii=0;iii<nx;iii++){
      result1[iii]=0;
      // jjj here is index in i.  We are doing integral over dx1
      for(jjj=0;jjj<=iii;jjj++){
	if((jjj==0)||(jjj==iii)) result1[iii]+=da[0*nx+jjj]*0.5;
	else result1[iii]+=da[0*nx+jjj];
      }
    }
    // now find i,j'th aphi
    for(iii=0;iii<nx*ny;iii++){
      indexi=(int)(iii%nx);    indexj=(int)(iii/nx);
      
      // here jjj is over dx2
      result2=0;
      for(jjj=0;jjj<=indexj;jjj++){
	if((jjj==indexj)||(jjj==0)) result2+=db[jjj*nx+indexi]*0.5;
	else result2+=db[jjj*nx+indexi];
      }
      
      aphi[indexj*nx+indexi]+=result1[indexi]+result2;
    }
#elif(0)
    // first set 0 point
    for(iii=0;iii<nx*ny;iii++)    aphi[iii]=-da[nx*(ny-1)];
    //for(iii=0;iii<nx*ny;iii++)    aphi[iii]=0;
    
    // generate x-direction integral since 1-D
    for(iii=0;iii<nx;iii++){
      result1[iii]=0;
      // jjj here is index in i.  We are doing integral over dx1
      for(jjj=0;jjj<=iii;jjj++){
	if((jjj==0)||(jjj==iii)) result1[iii]+=da[nx*(ny-1)+jjj]*0.5;
	else result1[iii]+=da[nx*(ny-1)+jjj];
      }
    }
    // now find i,j'th aphi
    for(iii=0;iii<nx*ny;iii++){
      indexi=(int)(iii%nx);    indexj=(int)(iii/nx);
      
      // here jjj is over dx2
      result2=0;
      for(jjj=ny-1;jjj>=indexj;jjj--){
	if((jjj==indexj)||(jjj==ny-1)) result2+=db[jjj*nx+indexi]*0.5;
	else result2+=db[jjj*nx+indexi];
      }
      
      aphi[indexj*nx+indexi]+=result1[indexi]+result2;
    }

#elif(0) // new Sasha-Jon algorithm
    // first set 0 point
    for(iii=0;iii<nx*ny;iii++)    aphi[iii]=0.0;
    
    // generate x-direction integral since 1-D
    for(iii=0;iii<nx;iii++){
      result1[iii]=-da[0]*0.5;
      // jjj here is index in i.  We are doing integral over dx1
      for(jjj=0;jjj<=iii;jjj++){
	if(jjj==iii) result1[iii]+=da[jjj]*0.5;
	else result1[iii]+=da[jjj];
      }
    }
    // now find i,j'th aphi
    for(iii=0;iii<nx*ny;iii++){
      indexi=(int)(iii%nx);    indexj=(int)(iii/nx);
      
      // here jjj is over dx2
      result2=-db[0*nx+indexi]*0.5;
      for(jjj=0;jjj<=indexj;jjj++){
	if(jjj==indexj) result2+=db[jjj*nx+indexi]*0.5;
	else result2+=db[jjj*nx+indexi];
      }
      
      aphi[indexj*nx+indexi]=result1[indexi]+result2;
    }

#elif(1) // new Sasha-Jon algorithm 2
    // first set 0 point
    for(iii=0;iii<nx*ny;iii++)    aphi[iii]=0.0;
    
    indexi=0;
    // generate x-direction integral since 1-D
    for(jnow=0;jnow<ny;jnow++){
      result1[jnow]=-db[0]*0.5;
      // jjj here is index in i.  We are doing integral over dx1
      for(indexj=0;indexj<=jnow;indexj++){
	if(indexj==jnow) result1[jnow]+=db[nx*indexj+indexi]*0.5;
	else result1[jnow]+=db[nx*indexj+indexi];
      }
    }
    // now find i,j'th aphi
    for(iii=0;iii<nx*ny;iii++){
      indexi=(int)(iii%nx);    indexj=(int)(iii/nx);
      
      // here jjj is over dx2
      result2=-da[indexj*nx+0]*0.5;
      for(inow=0;inow<=indexi;inow++){
	if(inow==indexi) result2+=da[indexj*nx+inow]*0.5;
	else result2+=da[indexj*nx+inow];
      }
      
      aphi[indexj*nx+indexi]=result1[indexj]+result2;
    }

#endif
  }
  else if(CALCTYPE==DODER){
    int signit;
    if(DERTYPE<10) signit=-1;
    else signit=1;

    for(iii=0;iii<nx*ny;iii++){
      i=indexi=(int)(iii%nx);    j=indexj=(int)(iii/nx);
      // x1
      if(nx>1){
	if( ((DERTYPE==2 || DERTYPE==12) || (indexi==0))&&(indexi!=nx-1) ){// forward difference
	  if(signit==-1) der[MAPDER(0)]=(p[MAPPIP1(0)]-p[MAPP(0)])/dx[1];
	  else  der[MAPDER(0)]=(fabs(p[MAPPIP1(0)])+fabs(p[MAPP(0)]))/dx[1];
	}
	else if( ((DERTYPE==1 || DERTYPE==11)||(indexi==nx-1))&&(indexi!=0) ){// backward difference
	  if(signit==-1) der[MAPDER(0)]=(p[MAPP(0)]-p[MAPPIM1(0)])/dx[1];
	  else der[MAPDER(0)]=(fabs(p[MAPP(0)])+fabs(p[MAPPIM1(0)]))/dx[1];
	}
	else if(DERTYPE==0 || DERTYPE==10) {// centered difference
	  if(signit==-1) der[MAPDER(0)]=0.5*(p[MAPPIP1(0)]-p[MAPPIM1(0)])/dx[1];
	  else der[MAPDER(0)]=0.5*(fabs(p[MAPPIP1(0)])+fabs(p[MAPPIM1(0)]))/dx[1];
	}
	else{
	  der[MAPDER(0)]=0.0;
	}
      }
      else{
	der[MAPDER(0)]=0.0;
      }
      // x2
      if(ny>1){
	if( ((DERTYPE==2 || DERTYPE==12) || (indexj==0))&&(indexj!=ny-1) ){// forward difference
	  if(signit==-1) der[MAPDER(1)]=(p[MAPPJP1(0)]-p[MAPP(0)])/dx[2];
	  else der[MAPDER(1)]=(fabs(p[MAPPJP1(0)])+fabs(p[MAPP(0)]))/dx[2];
	}
	else if( ((DERTYPE==1 || DERTYPE==11)||(indexj==ny-1))&&(indexj!=0) ){// backward difference
	  if(signit==-1) der[MAPDER(1)]=(p[MAPP(0)]-p[MAPPJM1(0)])/dx[2];
	  else der[MAPDER(1)]=(fabs(p[MAPP(0)])+fabs(p[MAPPJM1(0)]))/dx[2];
	}
	else if(DERTYPE==0 || DERTYPE==10){// centered difference
	  if(signit==-1) der[MAPDER(1)]=0.5*(p[MAPPJP1(0)]+signit*p[MAPPJM1(0)])/dx[2];
	  else der[MAPDER(1)]=0.5*(fabs(p[MAPPJP1(0)])+fabs(p[MAPPJM1(0)]))/dx[2];
	}
	else{
	  der[MAPDER(0)]=0.0;
	}
      }
      else{
	der[MAPDER(1)]=0.0;
      }
    }
  }
  else if(CALCTYPE==DOFICALC){


    for(dimen=1;dimen<=2;dimen++){
      for(iii=0;iii<nx*ny;iii++){
	i=indexi=(int)(iii%nx);    j=indexj=(int)(iii/nx);
#define NUMBC 3
	if(i>=NUMBC && j>=NUMBC && i<=nx-1-NUMBC && j<=ny-1-NUMBC){
	  // create V,P,Y
	  for(ismall=-NUMBC;ismall<=+NUMBC;ismall++){
	    V[ismall] = p[MAPPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2) , U1+dimen-1)];
	    //	    P[ismall] = p[MAPPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2) , UU)]; // assumes ideal gas and only taking ratios of pressures
	    P[ismall] = Ptot[MAPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2) , 0, 1)]; // true Ptot
	    PLOOP Ypl[k][ismall] = p[MAPPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2) , k)];

	    //	    fprintf(stderr,"i=%d j=%d ismall=%d V=%21.15g P=%21.15g\n",i,j,ismall,V[ismall],P[ismall]);
	  }

	  // now compute Ficalc for all k
	  //	  ficalc[MAPFICALC(dimen-1)]=Ficalc(dimen, V, P, Ypl);
	  ficalc[MAPFICALC(dimen-1)]=Ficalc(dimen, V, P);
	  //	  fprintf(stderr,"dimen=%d ficalc=%21.15g\n",dimen,ficalc[MAPFICALC(dimen-1)]);
	}
	else{
	  // then set to 0
	  ficalc[MAPFICALC(dimen-1)]=0.0;
	}
      }// end over iii
    }// end over dimen
  }
  fprintf(stderr,"done post-read calculations\n"); fflush(stderr);







  //////////////////////////////
  //
  // write calculations to file
  //
  //////////////////////////////

  if(CALCTYPE==DOAVG){
    fprintf(stderr,"start writing file: dumpname=%s\n",dumpnameout); fflush(stderr);
  }
  else{
    fprintf(stderr,"start writing file: dumpname=%s\n",filenameout); fflush(stderr);
  }

  for(iii=0;iii<nx*ny;iii++){
    i=indexi=(int)(iii%nx);    j=indexj=(int)((iii%(nx*ny))/nx);

    //    if(DUMPTYPE==REALDUMP){
    //  fprintf(dumpout,"%d %d ",i,j);
    // }

    if(CALCTYPE==DOFLINE){
      fprintf(dumpout,"%21.15g ",aphi[indexj*nx+indexi]);
    }
    else if(CALCTYPE==DOFARADAY){
      fprintf(dumpout,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "
	      ,faradaydd[indexj*nx*6+indexi*6+0]
	      ,faradaydd[indexj*nx*6+indexi*6+1]
	      ,faradaydd[indexj*nx*6+indexi*6+2]
	      ,faradaydd[indexj*nx*6+indexi*6+3]
	      ,faradaydd[indexj*nx*6+indexi*6+4]
	      ,faradaydd[indexj*nx*6+indexi*6+5]
	      ,faradayuu[indexj*nx*6+indexi*6+0]
	      ,faradayuu[indexj*nx*6+indexi*6+1]
	      ,faradayuu[indexj*nx*6+indexi*6+2]
	      ,faradayuu[indexj*nx*6+indexi*6+3]
	      ,faradayuu[indexj*nx*6+indexi*6+4]
	      ,faradayuu[indexj*nx*6+indexi*6+5]
	      );
    }
    else if(CALCTYPE==DODER){
      fprintf(dumpout,"%21.15g %21.15g "
	      ,der[MAPDER(0)]
	      ,der[MAPDER(1)]);
    }
    else if(CALCTYPE==DOFICALC){
      fprintf(dumpout,"%21.15g %21.15g ",ficalc[MAPFICALC(0)],ficalc[MAPFICALC(1)]);
    }
    else if(CALCTYPE==DOAVG){
      for(k=0;k<numgenvars;k++){
	fprintf(dumpout,"%21.15g ",genvar[MAPG(k)]);
      }
    }
    fprintf(dumpout,"\n");
  }



  //////////////////////////////
  //
  // free memory
  //
  //////////////////////////////
  if(CALCTYPE==DOAVG){
    free(genvar);
  }
  free(p);
  free(gdet);
  if(CALCTYPE==DOFLINE){
    free(aphi);
    free(da);
    free(db);
    free(result1);
  }
  else if(CALCTYPE==DOFARADAY){
    free(faradaydd);
    free(faradayuu);
  }
  else if(CALCTYPE==DODER){
    free(der);
  }
  else if(CALCTYPE==DOFICALC){
    free(ficalc);
  }
  fprintf(stderr,"done with all smcalc: dumpout=%lld\n",(long long int) dumpout); fflush(stderr);
  fclose(dumpout);

}





// PPM FLATTENER parameter
static FTYPE ftilde( int dimen, int shift, FTYPE *Vabs, FTYPE *Pabs)
{
  FTYPE Ftilde,Ftilde1,Ftilde2;
  FTYPE Sp;
  FTYPE *V, *P;
  FTYPE P2diff,Pbottom;


  // shift as needed
  P = Pabs + shift;
  V = Vabs + shift;

  // FLASH Equation 43
  P2diff=P[2]-P[-2];
  Pbottom=sign(P2diff)/(fabs(P2diff)+SMALL); // singularity avoidance but keeps signature
  Sp = (P[1] - P[-1]) * Pbottom ;



  // FLASH Equation 45
  Ftilde = max( 0, min( 1.0, 10.0 * (Sp - SP0) ) );

  //  Ftilde*=Ftilde*Ftilde*Ftilde;

  // FLASH Equation 46
  Ftilde1 = fabs(P[1] - P[-1]) / (min(fabs(P[1]), fabs(P[-1]))+ SMALL );
  Ftilde *= ( (FTYPE)(Ftilde1>=THIRD) );
  //  if(Ftilde1<THIRD) Ftilde=0.0;

  // FLASH Equation 47
  Ftilde2 = V[1] - V[-1];
  Ftilde *= ( (FTYPE)(Ftilde2<=0.0) );
  //  if(Ftilde2>0.0) Ftilde=0.0;
  
  //  if(Sp>0.8){
  //    fprintf(stderr,"Pbottom=%21.15g Sp=%21.15g Ftilde=%21.15g\n",Pbottom,Sp,Ftilde);
  //  }
  
  return( Ftilde );
}


// PPM FLATTENERS (final formula)
static FTYPE  Ficalc(int dimen, FTYPE *V, FTYPE *P)
{
  FTYPE ftilde( int dimen, int shift, FTYPE *P, FTYPE *V);
  int signdP;
  FTYPE Fi;

  signdP = (P[1] - P[-1] > 0) * 2 - 1;
  // FLASH Equation 48
  Fi = max( ftilde(dimen, 0, V,P), ftilde(dimen, -signdP, V,P) );

  return(Fi);
}

